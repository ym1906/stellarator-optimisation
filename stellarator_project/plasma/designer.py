# SPDX-FileCopyrightText: 2024-present {{ copyright-holder }} <{{ copyright-holder-email }}>
#
# SPDX-License-Identifier: MIT
"""Plasma Designer."""

from dataclasses import dataclass

from bluemira.base.designer import Designer
from bluemira.base.parameter_frame import ParameterFrame
from bluemira.geometry.parameterisations import GeometryParameterisation
from simsopt.geo import SurfaceRZFourier

from stellarator_project.tools import read_json, surface_to_nurbs


@dataclass
class PlasmaDesignerParams(ParameterFrame):
    """Parameters for designing a plasma."""


class PlasmaDesigner(Designer[GeometryParameterisation]):
    """Design a plasma's LCFS using a Johner parametrisation."""

    params: PlasmaDesignerParams
    param_cls = PlasmaDesignerParams

    def run(self) -> GeometryParameterisation:
        """Build the LCFS, returning a closed wire defining its outline."""
        s = SurfaceRZFourier.from_wout(
            self.build_config["input_nc3"],
            range="full torus",
            nphi=self.build_config.get("nphi", 64),
            ntheta=self.build_config.get("ntheta", 64),
        )
        surface_to_nurbs(
            simsopt_surface=s,
            export_path=self.build_config.get(
                f"{'finite_' if self.build_config.get('finite') else ''}surface_data"
            ),
            plot=False,
        )

        raise NotImplementedError("Not implemented")

    def read(self):
        return read_json(self.build_config["surface_data"])
