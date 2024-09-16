# SPDX-FileCopyrightText: 2024-present {{ copyright-holder }} <{{ copyright-holder-email }}>
#
# SPDX-License-Identifier: MIT

"""Stellarator Example."""

# %%
from pathlib import Path
from typing import Union

from bluemira.base.parameter_frame import EmptyFrame
from bluemira.base.reactor import Reactor
from bluemira.base.reactor_config import ReactorConfig

from stellarator_project.plasma.builder import PlasmaBuilder
from stellarator_project.plasma.desinger import PlasmaDesigner
from stellarator_project.plasma.manager import Plasma
from stellarator_project.tf_coil.builder import TFCoilBuilder
from stellarator_project.tf_coil.designer import TFCoilDesigner
from stellarator_project.tf_coil.manager import TFCoil

# %% [markdown]
#
# # Simplistic Reactor Design
#
# This example show hows to set up a simple reactor, consisting of a plasma and
# a single TF coil.
# The TF coil will be optimised such that its length is minimised,
# whilst maintaining a minimum distance to the plasma.
#
# To do this we'll run through how to set up the parameters for the build,
# how to define the `Builder`s and `Designer`s
# (including the optimisation problem) for the plasma and TF coil,
# and how to run the build with configurable parameters.
#


# %%
class MyReactor(Reactor):
    """A simple reactor with two components."""

    plasma: Plasma
    tf_coil: TFCoil


def main(build_config: Union[str, Path, dict]) -> MyReactor:  # noqa: FA100
    """Main reactor function."""
    reactor_config = ReactorConfig(build_config, EmptyFrame)

    # %% [markdown]
    #
    # We create our plasma

    # %%
    plasma_designer = PlasmaDesigner(
        reactor_config.params_for("Plasma", "designer"),
        reactor_config.config_for("Plasma", "designer"),
    )
    plasma_parameterisation = plasma_designer.execute()

    plasma_builder = PlasmaBuilder(
        plasma_parameterisation.create_shape(),
        reactor_config.config_for("Plasma"),
    )
    plasma = Plasma(plasma_builder.build())

    # %% [markdown]
    #
    # We create our TF coil

    # %%
    tf_coil_designer = TFCoilDesigner(
        plasma.lcfs(),
        None,
        reactor_config.config_for("TF Coil", "designer"),
    )
    tf_parameterisation = tf_coil_designer.execute()

    tf_coil_builder = TFCoilBuilder(
        reactor_config.params_for("TF Coil", "builder"),
        tf_parameterisation.create_shape(),
    )
    tf_coil = TFCoil(tf_coil_builder.build())

    # %% [markdown]
    #
    # Finally we add the components to the reactor and show the CAD

    # %%
    reactor = MyReactor("Simple Example", n_sectors=1)

    reactor.plasma = plasma
    reactor.tf_coil = tf_coil

    reactor.show_cad(n_sectors=1)
    reactor.show_cad("xz")

    return reactor
