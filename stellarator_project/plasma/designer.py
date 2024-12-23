# SPDX-FileCopyrightText: 2024-present United Kingdom Atomic Energy Authority
#
# SPDX-License-Identifier: MIT
"""Plasma Designer."""

from __future__ import annotations

from dataclasses import dataclass

from bluemira.base.designer import Designer
from bluemira.base.parameter_frame import ParameterFrame
from bluemira.geometry.parameterisations import GeometryParameterisation
from simsopt.geo import SurfaceRZFourier

from stellarator_project.tools import NURBSData, read_json, surface_to_nurbs


@dataclass
class PlasmaDesignerParams(ParameterFrame):
    """Parameters for designing a plasma."""


class PlasmaDesigner(Designer[GeometryParameterisation]):
    """Design a plasma's LCFS using a Johner parametrisation."""

    params: PlasmaDesignerParams
    param_cls = PlasmaDesignerParams

    def run(self) -> tuple[NURBSData, SurfaceRZFourier]:
        """Build the LCFS, returning a closed wire defining its outline."""
        s = SurfaceRZFourier.from_wout(
            self.build_config["input_nc3"],
            range=SurfaceRZFourier.RANGE_FIELD_PERIOD,
            nphi=self.build_config.get("nphi", 64),
            ntheta=self.build_config.get("ntheta", 64),
        )
        return surface_to_nurbs(
            simsopt_surface=s,
            export_path=self.build_config.get("surface_data"),
            plot=True,
        ), s

    def read(self) -> tuple[NURBSData, None]:
        """Read surface data from json."""
        return NURBSData(**read_json(self.build_config["surface_data"])), None
