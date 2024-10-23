# SPDX-FileCopyrightText: 2024-present United Kingdom Atomic Energy Authority
#
# SPDX-License-Identifier: MIT
"""TF Coil Builder."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from bluemira.base.builder import Builder
from bluemira.base.components import Component, PhysicalComponent
from bluemira.base.look_and_feel import bluemira_debug
from bluemira.base.parameter_frame import ParameterFrame
from bluemira.display.palettes import BLUE_PALETTE
from bluemira.geometry.tools import loft_shape, make_bspline

if TYPE_CHECKING:
    from stellarator_project.tools import NURBSCurveData

# The TF coil builder is then passed the centreline from the designer to create
# the Component and the CAD of the TF coil.
# If more TF coils were to be required the build_xyz of `TFCoilBuilder` would
# need to be modified.
#
# Notice that only `TFCoilBuilder` has required parameters in this case.


@dataclass
class TFCoilBuilderParams(ParameterFrame):
    """Parameters for building a TF coil."""


class TFCoilBuilder(Builder):
    """Build a 3D model of a TF Coil from a given centre line."""

    params: TFCoilBuilderParams
    param_cls = TFCoilBuilderParams

    def __init__(self, params: TFCoilBuilderParams, coils_data):
        super().__init__(params, {})
        self.coils_data = coils_data

    def build(self) -> Component:
        """Run the full build for the TF coils."""
        return self.component_tree(
            xz=None,
            xy=None,
            xyz=self.build_xyz(),
        )

    def build_xyz(self) -> list[PhysicalComponent]:
        """Build the xyz Component of the TF coils."""
        coils = []
        for i, (name, filaments) in enumerate(self.coils_data.items()):
            coil = PhysicalComponent(
                f"Winding pack {i}", self.create_coil_from_nurbs(name, filaments)
            )

            coil.display_cad_options.color = BLUE_PALETTE["TF"]
            coils.append(coil)
        return coils

    @staticmethod
    def create_coil_from_nurbs(name: str, filaments: dict[str, NURBSCurveData]):
        """Create coil from nurb filaments."""
        filament_wires = [
            make_bspline(
                poles=curve_dict.poles,
                mults=curve_dict.mults,
                knots=curve_dict.internal_knot_vector,
                degree=curve_dict.degree,
                weights=curve_dict.weights,
                periodic=False,
                check_rational=False,
                label=filament_name,
            )
            for filament_name, curve_dict in filaments.items()
        ]

        # Order wires rotationally. Previous input was in rows to create positive volumes
        # TODO @je-cook: This is not robust and should be improved if possible
        # 7
        size = len(filament_wires) // 2
        filament_wires = filament_wires[:size] + filament_wires[size:][::-1]

        # Create the lofted surface from the filament curves
        coil = loft_shape(
            filament_wires,
            solid=True,
            ruled=True,
            closed=True,
            label=name,
        )
        bluemira_debug(f"L: {coil.length} A: {coil.area} V {coil.volume}")

        return coil
