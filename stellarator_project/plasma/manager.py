# SPDX-FileCopyrightText: 2024-present United Kingdom Atomic Energy Authority
#
# SPDX-License-Identifier: MIT
"""Plasma Manager."""

from __future__ import annotations

from typing import TYPE_CHECKING

from bluemira.base.reactor import ComponentManager

if TYPE_CHECKING:
    from bluemira.base.component import Component
    from simsopt.geo import SurfaceRZFourier


class Plasma(ComponentManager):
    """Plasma component manager."""

    def __init__(
        self,
        component: Component,
        surface: SurfaceRZFourier | None = None,
        n_field_periods: int | None = None,
    ):
        super().__init__(component)
        self.surface = surface
        self.n_field_periods = n_field_periods or surface.nfp
        if n_field_periods is None:
            raise ValueError("A Surface or number of field periods needs to be provided")
