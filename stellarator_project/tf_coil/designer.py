# SPDX-FileCopyrightText: 2024-present {{ copyright-holder }} <{{ copyright-holder-email }}>
#
# SPDX-License-Identifier: MIT
"""TF Coil Designer."""

import numpy as np
from bluemira.base.designer import Designer
from bluemira.base.reactor_config import ConfigParams
from bluemira.geometry.optimisation import optimise_geometry
from bluemira.geometry.parameterisations import GeometryParameterisation
from bluemira.geometry.tools import (
    distance_to,
)
from bluemira.geometry.wire import BluemiraWire
from bluemira.utilities.tools import get_class_from_module

# And now the TF Coil, in this instance for simplicity we are only making
# one TF coil.
#
# The TF coil designer finds the geometry parameterisation given in
# the `build_config` which should point to a class.
# The parameterisation is then fed into the optimiser that
# minimises the size of the TF coil, whilst keeping at least a meter away
# from the plasma at any point.
# Further information on geometry and geometry optimisations can be found in the
# [geometry tutorial](../geometry/geometry_tutorial.ex.py) and
# [geometry optimisation tutorial](../optimisation/geometry_optimisation.ex.py).
#


class TFCoilDesigner(Designer[GeometryParameterisation]):
    """TF coil shape designer."""

    param_cls = None  # This designer takes no parameters

    def __init__(
        self,
        plasma_lcfs: BluemiraWire,
        params: None,
        build_config: ConfigParams,
    ):
        super().__init__(params, build_config)
        self.lcfs = plasma_lcfs
        self.parameterisation_cls = get_class_from_module(
            self.build_config["param_class"],
            default_module="bluemira.geometry.parameterisations",
        )

    def run(self) -> GeometryParameterisation:
        """Run the design of the TF coil."""
        parameterisation = self.parameterisation_cls(
            var_dict=self.build_config["var_dict"],
        )
        min_dist_to_plasma = 1  # meter
        return self.minimise_tf_coil_size(parameterisation, min_dist_to_plasma)

    def minimise_tf_coil_size(
        self,
        geom: GeometryParameterisation,
        min_dist_to_plasma: float,
    ) -> GeometryParameterisation:
        """Run an optimisation to minimise the size of the TF coil.

        We're minimising the size of the coil whilst always keeping a
        minimum distance to the plasma.
        """
        distance_constraint = {
            "f_constraint": lambda g: self._constrain_distance(
                g,
                min_dist_to_plasma,
            ),
            "tolerance": np.array([1e-6]),
        }
        optimisation_result = optimise_geometry(
            geom=geom,
            f_objective=lambda g: g.create_shape().length,
            opt_conditions={"max_eval": 500, "ftol_rel": 1e-6},
            ineq_constraints=[distance_constraint],
        )
        return optimisation_result.geom

    def _constrain_distance(
        self,
        geom: BluemiraWire,
        min_distance: float,
    ) -> np.ndarray:
        return np.array(
            min_distance - distance_to(geom.create_shape(), self.lcfs)[0],
        )
