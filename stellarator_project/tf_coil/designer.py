# SPDX-FileCopyrightText: 2024-present {{ copyright-holder }} <{{ copyright-holder-email }}>
#
# SPDX-License-Identifier: MIT
"""TF Coil Designer."""

import numpy as np
from bluemira.base.designer import Designer
from bluemira.base.reactor_config import ConfigParams
from bluemira.geometry.optimisation import optimise_geometry
from bluemira.geometry.parameterisations import GeometryParameterisation
from bluemira.geometry.tools import distance_to
from bluemira.geometry.wire import BluemiraWire
from bluemira.utilities.tools import get_class_from_module

from stellarator_project.tools import (
    curves_to_nurbs,
    filament_curves_to_nurbs,
    read_json,
)


class TFCoilDesigner(Designer[GeometryParameterisation]):
    """TF coil shape designer."""

    param_cls = None  # This designer takes no parameters

    def __init__(
        self,
        params: None,
        build_config: ConfigParams,
    ):
        super().__init__(params, build_config)

    def taylor_test(self):
        # Perform a Taylor test
        print("Performing a Taylor test")
        f = fun
        dofs = JF.x
        np.random.seed(1)
        h = np.random.uniform(size=dofs.shape)
        J0, dJ0 = f(dofs)
        dJh = sum(dJ0 * h)
        for eps in [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]:
            J1, _ = f(dofs + eps * h)
            J2, _ = f(dofs - eps * h)
            print("err", (J1 - J2) / (2 * eps) - dJh)

    def run(self) -> GeometryParameterisation:
        """Run the design of the TF coil."""
        self.taylor_test()
        # Run the optimization
        print("Running the optimization")
        res = minimize(
            fun,
            dofs,
            jac=True,
            method="L-BFGS-B",
            options={"maxiter": MAXITER, "maxcor": 400, "gtol": 1e-20, "ftol": 1e-20},
            tol=1e-20,
        )

        # Perform a second optimization with reduced penalty for coil length
        dofs = res.x
        LENGTH_WEIGHT *= 0.1
        res = minimize(
            fun,
            dofs,
            jac=True,
            method="L-BFGS-B",
            options={"maxiter": MAXITER, "maxcor": 300},
            tol=1e-15,
        )

        # Save the optimized coil shapes and currents
        bs.save(OUT_DIR + "biot_savart_opt.json")

        curves_to_nurbs(
            curves=curves,
            export_path=self.build_config.get(
                f"{'finite_' if self.build_config.get('finite') else ''}coil_data"
            ),
            plot=False,
        )
        filament_curves_to_nurbs(
            curves=curves,
            numfil=nfil,
            export_path=self.build_config.get(
                f"{'finite_' if self.build_config.get('finite') else ''}coil_data"
            ),
            plot=False,
        )

        raise NotImplementedError

    def read(self):
        # curve_surface_normals_data = read_json(curve_surface_normals)
        return read_json(self.build_config["coil_data"])
