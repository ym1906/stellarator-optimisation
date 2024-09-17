# SPDX-FileCopyrightText: 2024-present {{ copyright-holder }} <{{ copyright-holder-email }}>
#
# SPDX-License-Identifier: MIT
"""Plasma Manager."""

from bluemira.base.reactor import ComponentManager


class Plasma(ComponentManager):
    """Plasma component manager."""

    def __init__(self, component, surface):
        super().__init__(component)
        self.surface = surface
