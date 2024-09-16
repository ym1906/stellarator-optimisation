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

from stellarator_project.tools import read_json


class TFCoilDesigner(Designer[GeometryParameterisation]):
    """TF coil shape designer."""

    param_cls = None  # This designer takes no parameters

    def __init__(
        self,
        params: None,
        build_config: ConfigParams,
    ):
        super().__init__(params, build_config)

    def run(self) -> GeometryParameterisation:
        """Run the design of the TF coil."""
        raise NotImplementedError

    def read(self):
        # curve_surface_normals_data = read_json(curve_surface_normals)
        return read_json(self.build_config["coil_data"])
