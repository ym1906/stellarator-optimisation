# SPDX-FileCopyrightText: 2024-present {{ copyright-holder }} <{{ copyright-holder-email }}>
#
# SPDX-License-Identifier: MIT
"""Plasma Designer."""

from dataclasses import dataclass

from bluemira.base.designer import Designer
from bluemira.base.parameter_frame import Parameter, ParameterFrame
from bluemira.equilibria.shapes import JohnerLCFS
from bluemira.geometry.parameterisations import GeometryParameterisation

# ## ParameterFrames and ComponentMangers
# Firstly we need to define the parameters we're going to use in our reactor
# design for each component.
# See the [ParameterFrame example](../base/params.ex.py) for more information.


# We need to define some `Designers` and `Builders` for our various
# `Components`.
#
# Firstly the plasma.
# The plasma designer will, using its `ParameterFrame`, evaluate a `JohnerLCFS`
# geometry parameterisation, returning a wire representing the plasma's
# last-closed-flux-surface (LCFS).
#
# In this case `PlasmaDesigner` has some required parameters but `PlasmaBuilder`
# does not.


@dataclass
class PlasmaDesignerParams(ParameterFrame):
    """Parameters for designing a plasma."""

    R_0: Parameter[float]
    A: Parameter[float]
    z_0: Parameter[float]
    kappa_u: Parameter[float]
    kappa_l: Parameter[float]
    delta_u: Parameter[float]
    delta_l: Parameter[float]
    phi_neg_u: Parameter[float]
    phi_pos_u: Parameter[float]
    phi_pos_l: Parameter[float]
    phi_neg_l: Parameter[float]


class PlasmaDesigner(Designer[GeometryParameterisation]):
    """Design a plasma's LCFS using a Johner parametrisation."""

    params: PlasmaDesignerParams
    param_cls = PlasmaDesignerParams

    def run(self) -> GeometryParameterisation:
        """Build the LCFS, returning a closed wire defining its outline."""
        return self._build_wire(self.params)

    @staticmethod
    def _build_wire(params: PlasmaDesignerParams) -> GeometryParameterisation:
        return JohnerLCFS(
            var_dict={
                "r_0": {"value": params.R_0.value},
                "z_0": {"value": params.z_0.value},
                "a": {"value": params.R_0.value / params.A.value},
                "kappa_u": {"value": params.kappa_u.value},
                "kappa_l": {"value": params.kappa_l.value},
                "delta_u": {"value": params.delta_u.value},
                "delta_l": {"value": params.delta_l.value},
                "phi_u_pos": {
                    "value": params.phi_pos_u.value,
                    "lower_bound": 0.0,
                },
                "phi_u_neg": {
                    "value": params.phi_neg_u.value,
                    "lower_bound": 0.0,
                },
                "phi_l_pos": {
                    "value": params.phi_pos_l.value,
                    "lower_bound": 0.0,
                },
                "phi_l_neg": {
                    "value": params.phi_neg_l.value,
                    "lower_bound": 0.0,
                    "upper_bound": 90,
                },
            },
        )
