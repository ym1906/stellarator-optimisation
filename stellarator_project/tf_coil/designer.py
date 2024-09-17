# SPDX-FileCopyrightText: 2024-present {{ copyright-holder }} <{{ copyright-holder-email }}>
#
# SPDX-License-Identifier: MIT
"""TF Coil Designer."""

from __future__ import annotations

import functools
import operator
from dataclasses import dataclass
from itertools import starmap
from typing import TYPE_CHECKING

import numpy as np
from bluemira.base.designer import Designer
from bluemira.base.look_and_feel import bluemira_print, bluemira_print_flush
from bluemira.base.parameter_frame import Parameter, ParameterFrame
from bluemira.geometry.parameterisations import GeometryParameterisation
from scipy.optimize import minimize
from simsopt.field import (
    BiotSavart,
    Coil,
    Current,
    apply_symmetries_to_currents,
    apply_symmetries_to_curves,
)
from simsopt.geo import (
    CurveCurveDistance,
    CurveLength,
    create_equally_spaced_curves,
    create_multifilament_grid,
)
from simsopt.objectives import QuadraticPenalty, SquaredFlux

from stellarator_project.tools import (
    NURBSCurveData,
    filament_curves_to_nurbs,
    read_json,
    taylor_test,
)

if TYPE_CHECKING:
    from bluemira.base.parameter_frame.typed import ParameterFrameLike
    from bluemira.base.reactor_config import ConfigParams


@dataclass
class TFCoilDesignerParams(ParameterFrame):
    """Parameters for building a TF coil."""

    R_0: Parameter[float]
    """Major radius for the initial circular coils"""
    R_minor: Parameter[float]
    """Minor radius for the initial circular coils"""
    n_TF_types: Parameter[int]  # noqa: N815
    """Number of unique coil shapes"""
    fourier_modes_cart: Parameter[int]
    """Number of Fourier modes for each Cartesian component of each coil"""
    fourier_modes_rot: Parameter[int]
    """Order of the Fourier expression for the rotation of the filament pack"""
    base_current: Parameter[float]
    min_coil_sep: Parameter[float]
    fil_gap_n: Parameter[float]
    """Gap between filaments in normal direction"""
    fil_gap_b: Parameter[float]
    """Gap between filaments in bi-normal direction"""
    n_fil_n: Parameter[int]
    """Number of filaments in normal direction"""
    n_fil_b: Parameter[int]
    """Number of filaments in bi-normal direction"""


class TFCoilDesigner(Designer[GeometryParameterisation]):
    """TF coil shape designer."""

    params: TFCoilDesignerParams
    param_cls: type[TFCoilDesignerParams] = TFCoilDesignerParams

    def __init__(
        self,
        params: ParameterFrameLike,
        build_config: ConfigParams,
        plasma,
    ):
        super().__init__(params, build_config)
        self.plasma = plasma
        self.coil_data = self.build_config["coil_data"]

    def run(self) -> GeometryParameterisation:  # noqa: PLR0914, PLR0915
        """Run the design of the TF coil."""
        # Create the initial coils
        ncoils = self.params.n_TF_types.value
        r0 = self.params.R_0.value
        r1 = self.params.R_minor.value
        order = self.params.fourier_modes_cart.value

        base_curves = create_equally_spaced_curves(
            ncoils, self.plasma.nfp, stellsym=True, R0=r0, R1=r1, order=order
        )
        base_currents = [Current(self.params.base_current.value) for i in range(ncoils)]
        base_currents[0].fix_all()  # Fix the first current to avoid trivial solution

        # use sum here to concatenate lists
        numfilaments_n = self.params.n_fil_n.value
        numfilaments_b = self.params.n_fil_b.value
        gapsize_n = self.params.fil_gap_n.value
        gapsize_b = self.params.fil_gap_b.value
        rot_order = self.params.fourier_modes_rot.value
        nfil = numfilaments_n * numfilaments_b

        length_pen = self.build_config["optimisation"]["length_penalty"]
        # Threshold and weight for the coil-to-coil distance penalty in the objective
        # function:
        dist_min = self.params.min_coil_sep.value
        dist_pen = self.build_config["optimisation"]["distance_penalty"]

        # Create a dictionary to hold filaments for each base coil
        coil_filament_map = {}

        # Create the filaments and group them by coil name
        for i, base_curve in enumerate(base_curves):
            coil_name = f"coil_{i + 1}"
            filaments = create_multifilament_grid(
                base_curve,
                numfilaments_n,
                numfilaments_b,
                gapsize_n,
                gapsize_b,
                rotation_order=rot_order,
            )
            coil_filament_map[coil_name] = filaments
        base_curves_finite_build = functools.reduce(
            operator.iadd, coil_filament_map.values(), []
        )
        base_currents_finite_build = functools.reduce(
            operator.iadd, [[c] * nfil for c in base_currents], []
        )
        curves_fb = apply_symmetries_to_curves(
            base_curves_finite_build, self.plasma.nfp, stellsym=True
        )
        currents_fb = apply_symmetries_to_currents(
            base_currents_finite_build, self.plasma.nfp, stellsym=True
        )
        coils_fb = list(starmap(Coil, zip(curves_fb, currents_fb, strict=False)))
        bs = BiotSavart(coils_fb)
        bs.set_points(self.plasma.gamma().reshape((-1, 3)))

        # coils = coils_via_symmetries(base_curves, base_currents, self.plasma.nfp, True)
        curves = [c.curve for c in coils_fb]

        # Define the individual terms of the objective function
        # Define the objective function:
        Jf = SquaredFlux(self.plasma, bs)  # noqa: N806
        Jls = [CurveLength(c) for c in base_curves]  # noqa: N806
        Jdist = CurveCurveDistance(curves, dist_min)  # noqa: N806

        # Form the total objective function
        JF = (  # noqa: N806
            Jf
            + length_pen
            * sum(
                QuadraticPenalty(Jls[i], Jls[i].J(), "max")
                for i in range(len(base_curves))
            )
            + dist_pen * Jdist
        )

        def fun(dofs):
            JF.x = dofs
            J = JF.J()  # noqa: N806
            grad = JF.dJ()

            cl_string = ", ".join([f"{J.J():.3f}" for J in Jls])
            mean_absB = np.mean(bs.AbsB())  # noqa: N806
            jf = Jf.J()
            kap_string = ", ".join(f"{np.max(c.kappa()):.1f}" for c in base_curves)
            bluemira_print_flush(
                f"J={J:.3e}, Jflux={jf:.3e},"
                f"sqrt(Jflux)/Mean(|B|)={np.sqrt(jf) / mean_absB:.3e},"
                f" CoilLengths=[{cl_string}], [{kap_string}], "
                f"||∇J||={np.linalg.norm(grad):.3e}"
            )
            return 1e-4 * J, 1e-4 * grad

        taylor_test(fun, JF)
        # Run the optimization
        dofs = JF.x

        bluemira_print("Running the optimization")
        first_opt = self.build_config["optimisation"]["first"]
        res = minimize(
            fun,
            dofs,
            jac=True,
            method=first_opt.get("algorithm", "L-BFGS-B"),
            options=first_opt.get("options", {}),
            tol=first_opt.get("tol", 1e-20),
        )

        # Perform a second optimization with reduced penalty for coil length
        dofs = res.x
        second_opt = self.build_config["optimisation"]["second"]
        res = minimize(
            fun,
            dofs,
            jac=True,
            method=second_opt.get("algorithm", "L-BFGS-B"),
            options=second_opt.get("options", {}),
            tol=second_opt.get("tol", 1e-15),
        )

        # Save the optimized coil shapes and currents
        bs.save(self.build_config["bs_opt"])

        return filament_curves_to_nurbs(
            curves=curves,
            numfil=nfil,
            export_path=self.coil_data,
            plot=True,
        )

    def read(self) -> dict[str, dict[str, NURBSCurveData]]:
        """Read coil data from json."""
        # curve_surface_normals_data = read_json(curve_surface_normals)
        return {
            coil: {
                filament: NURBSCurveData(**f_data) for filament, f_data in data.items()
            }
            for coil, data in read_json(self.coil_data).items()
        }
