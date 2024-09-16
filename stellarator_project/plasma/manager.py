# SPDX-FileCopyrightText: 2024-present {{ copyright-holder }} <{{ copyright-holder-email }}>
#
# SPDX-License-Identifier: MIT
"""Plasma Manager."""

from bluemira.base.reactor import ComponentManager
from bluemira.geometry.wire import BluemiraWire

# To manage access to properties of the components we need some
# `ComponentManagers`


class Plasma(ComponentManager):
    """Plasma component manager."""

    def lcfs(self) -> BluemiraWire:
        """Get separatrix."""
        return (
            self.component().get_component("xz").get_component("LCFS").shape.boundary[0]
        )
