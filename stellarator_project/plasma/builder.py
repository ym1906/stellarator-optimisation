# SPDX-FileCopyrightText: 2024-present {{ copyright-holder }} <{{ copyright-holder-email }}>
#
# SPDX-License-Identifier: MIT
"""Plasma Builder."""

from __future__ import annotations

from bluemira.base.builder import Builder
from bluemira.base.components import Component, PhysicalComponent
from bluemira.display.palettes import BLUE_PALETTE
from bluemira.geometry.tools import make_bsplinesurface


class PlasmaBuilder(Builder):
    """Build the 3D geometry of a plasma from a given LCFS."""

    param_cls = None

    def __init__(self, surface_data: dict, build_config: dict):
        super().__init__(None, build_config)
        self.s_data = surface_data

    def build(self) -> Component:
        """Run the full build of the Plasma."""
        return self.component_tree(
            xz=None,
            xy=None,
            xyz=[self.build_xyz()],
        )

    def build_xyz(self) -> PhysicalComponent:
        """Build the plasma."""
        # Create a plasma surface from NURBS surface data
        plasma_surface = make_bsplinesurface(
            poles=self.s_data.poles2d,
            mults_u=self.s_data.mults_u,
            mults_v=self.s_data.mults_v,
            knot_vector_u=self.s_data.internal_knot_vector_u,
            knot_vector_v=self.s_data.internal_knot_vector_v,
            degree_u=self.s_data.degree_u,
            degree_v=self.s_data.degree_v,
            weights=self.s_data.weights_reshaped,
            periodic=False,
            check_rational=False,
        )

        component = PhysicalComponent("LCFS", plasma_surface)
        component.display_cad_options.color = BLUE_PALETTE["PL"]
        component.display_cad_options.transparency = 0.5
        return component
