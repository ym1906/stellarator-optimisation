# SPDX-FileCopyrightText: 2024-present {{ copyright-holder }} <{{ copyright-holder-email }}>
#
# SPDX-License-Identifier: MIT
"""Plasma Designer."""

from dataclasses import dataclass

from bluemira.base.designer import Designer
from bluemira.base.parameter_frame import ParameterFrame
from bluemira.geometry.parameterisations import GeometryParameterisation

from stellarator_project.tools import read_json


@dataclass
class PlasmaDesignerParams(ParameterFrame):
    """Parameters for designing a plasma."""


class PlasmaDesigner(Designer[GeometryParameterisation]):
    """Design a plasma's LCFS using a Johner parametrisation."""

    params: PlasmaDesignerParams
    param_cls = PlasmaDesignerParams

    def run(self) -> GeometryParameterisation:
        """Build the LCFS, returning a closed wire defining its outline."""
        raise NotImplementedError("Not implemented")

    def read(self):
        return read_json(self.build_config["surface_data"])
