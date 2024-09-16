# SPDX-FileCopyrightText: 2024-present {{ copyright-holder }} <{{ copyright-holder-email }}>
#
# SPDX-License-Identifier: MIT
"""Plasma Builder."""

from __future__ import annotations

from typing import TYPE_CHECKING

from bluemira.base.builder import Builder
from bluemira.base.components import Component, PhysicalComponent
from bluemira.display.palettes import BLUE_PALETTE
from bluemira.geometry.face import BluemiraFace
from bluemira.geometry.tools import (
    revolve_shape,
)

if TYPE_CHECKING:
    from bluemira.geometry.wire import BluemiraWire


class PlasmaBuilder(Builder):
    """Build the 3D geometry of a plasma from a given LCFS."""

    param_cls = None

    def __init__(self, wire: BluemiraWire, build_config: dict):
        super().__init__(None, build_config)
        self.wire = wire

    def build(self) -> Component:
        """Run the full build of the Plasma."""
        xz = self.build_xz()
        return self.component_tree(
            xz=[xz],
            xy=[Component("")],
            xyz=[self.build_xyz(xz.shape)],
        )

    def build_xz(self) -> PhysicalComponent:
        """Build a view of the plasma in the toroidal (xz) plane.

        This generates a ``PhysicalComponent``, whose shape is a face.
        """
        component = PhysicalComponent("LCFS", BluemiraFace(self.wire))
        component.display_cad_options.color = BLUE_PALETTE["PL"]
        component.display_cad_options.transparency = 0.5
        return component

    @staticmethod
    def build_xyz(lcfs: BluemiraFace) -> PhysicalComponent:
        """Build the 3D (xyz) Component of the plasma by revolving
        the given face 360 degrees.
        """
        shape = revolve_shape(lcfs, degree=360)
        component = PhysicalComponent("LCFS", shape)
        component.display_cad_options.color = BLUE_PALETTE["PL"]
        component.display_cad_options.transparency = 0.5
        return component
