# SPDX-FileCopyrightText: 2024-present {{ copyright-holder }} <{{ copyright-holder-email }}>   # noqa: E501
#
# SPDX-License-Identifier: MIT
"""TF Coil Builder."""

from dataclasses import dataclass

from bluemira.base.builder import Builder
from bluemira.base.components import Component, PhysicalComponent
from bluemira.base.parameter_frame import Parameter, ParameterFrame
from bluemira.geometry.face import BluemiraFace
from bluemira.geometry.tools import (
    make_polygon,
    offset_wire,
    sweep_shape,
)
from bluemira.geometry.wire import BluemiraWire

# The TF coil builder is then passed the centreline from the designer to create
# the Component and the CAD of the TF coil.
# If more TF coils were to be required the build_xyz of `TFCoilBuilder` would
# need to be modified.
#
# Notice that only `TFCoilBuilder` has required parameters in this case.


@dataclass
class TFCoilBuilderParams(ParameterFrame):
    """Parameters for building a TF coil."""

    tf_wp_width: Parameter[float]
    tf_wp_depth: Parameter[float]


class TFCoilBuilder(Builder):
    """Build a 3D model of a TF Coil from a given centre line."""

    params: TFCoilBuilderParams
    param_cls = TFCoilBuilderParams

    def __init__(self, params: TFCoilBuilderParams, centreline: BluemiraWire):
        super().__init__(params, {})
        self.centreline = centreline

    def make_tf_wp_xs(self) -> BluemiraWire:
        """Make a wire for the cross-section of the winding pack in xy."""
        width = 0.5 * self.params.tf_wp_width.value
        depth = 0.5 * self.params.tf_wp_depth.value
        return make_polygon(
            {
                "x": [-width, width, width, -width],
                "y": [-depth, -depth, depth, depth],
                "z": 0.0,
            },
            closed=True,
        )

    def build(self) -> Component:
        """Run the full build for the TF coils."""
        return self.component_tree(
            xz=[self.build_xz()],
            xy=[Component("")],
            xyz=[self.build_xyz()],
        )

    def build_xz(self) -> PhysicalComponent:
        """Build the xz Component of the TF coils."""
        inner = offset_wire(
            self.centreline,
            -0.5 * self.params.tf_wp_width.value,
        )
        outer = offset_wire(
            self.centreline,
            0.5 * self.params.tf_wp_width.value,
        )
        return PhysicalComponent("Winding pack", BluemiraFace([outer, inner]))

    def build_xyz(self) -> PhysicalComponent:
        """Build the xyz Component of the TF coils."""
        wp_xs = self.make_tf_wp_xs()
        wp_xs.translate((self.centreline.bounding_box.x_min, 0, 0))
        volume = sweep_shape(wp_xs, self.centreline)
        return PhysicalComponent("Winding pack", volume)
