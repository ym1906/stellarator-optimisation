# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: tags,-all
#     notebook_metadata_filter: -jupytext.text_representation.jupytext_version
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% tags=["remove-cell"]
# SPDX-FileCopyrightText: 2021-present M. Coleman, J. Cook, F. Franza
# SPDX-FileCopyrightText: 2021-present I.A. Maione, S. McIntosh
# SPDX-FileCopyrightText: 2021-present J. Morris, D. Short
#
# SPDX-License-Identifier: LGPL-2.1-or-later

"""
A stellarator build example.
"""

# %% [markdown]
# # Geometry Tutorial
# ## Introduction
#
# Geometry is not plasma physics, but it isn't trivial either. Chances are most of
# your day-to-day interaction with bluemira will revolve around geometry in some form
# or another. Puns intended.
#
# There a few basic concepts you need to familiarise yourself with:
# * Basic objects: [`BluemiraWire`, `BluemiraFace`, `BluemiraShell`, `BluemiraSolid`]
# * Basic properties
# * Matryoshka structure
# * Geometry creation
# * Geometry modification
# * Geometry operations
#
# ## Imports
#
# Let's start out by importing all the basic objects, and some typical tools

# %%
import json
import os
import sys
from typing import Any, List

sys.path.append(os.path.abspath("../bluemira"))


# Some display functionality
from bluemira.codes import _freecadapi as cadapi
from bluemira.display import show_cad

# Basic objects
import Part

# Some useful tools
from bluemira.geometry.tools import (
    make_bsplinesurface,
    save_cad,
    loft_shape,
    make_bspline,
)
from bluemira.geometry.wire import BluemiraWire

# %% [markdown]
# Creating a 3D representation of a stellarator plasma surface generated by Simsopt


def read_json(file_path: str) -> dict[str, Any]:
    """Read JSON data from a file."""
    with open(file_path) as f:
        return json.load(f)


surface_filename = (
    "stellarator_project/data/plasma/finite_plasma_surface_nurbs_data.json"
)
coil_filename = "stellarator_project/data/coils/finite_magnets_nurbs_data.json"
curve_surface_normals = "stellarator_project/data/coils/normals_data.json"
# Read in the json data
surface_data = read_json(surface_filename)
coil_data = read_json(coil_filename)
curve_surface_normals_data = read_json(curve_surface_normals)

# Create a plasma surface from NURBS surface data
plasma_surface = make_bsplinesurface(
    poles=surface_data["poles2d"],
    mults_u=surface_data["mults_u"],
    mults_v=surface_data["mults_v"],
    knot_vector_u=surface_data["internal_knot_vector_u"],
    knot_vector_v=surface_data["internal_knot_vector_v"],
    degree_u=surface_data["degree_u"],
    degree_v=surface_data["degree_v"],
    weights=surface_data["weights_reshaped"],
    periodic=False,
    check_rational=False,
)
coil_curves = []


class PartWrapper:
    def __init__(self, shape):
        self.shape = shape


def create_coil_from_nurbs(json_path: str) -> List[Any]:
    with open(json_path, "r") as f:
        coils_data = json.load(f)

    coil_curves = []
    coil = []
    for coil in coils_data:
        for coil_name, filaments in coil.items():
            filament_curves = []
            filament_wires = []
            for filament_name, curve_dict in filaments.items():
                coil_filament_curve = cadapi.make_bspline(
                    poles=curve_dict["poles"],
                    mults=curve_dict["mults"],
                    knots=curve_dict["internal_knot_vector"],
                    degree=curve_dict["degree"],
                    weights=curve_dict["weights"],
                    periodic=False,
                    check_rational=False,
                )
                coil_wire = make_bspline(
                    poles=curve_dict["poles"],
                    mults=curve_dict["mults"],
                    knots=curve_dict["internal_knot_vector"],
                    degree=curve_dict["degree"],
                    weights=curve_dict["weights"],
                    periodic=False,
                    check_rational=False,
                    label=filament_name,
                )
                filament_curves.append(coil_filament_curve)
                filament_wires.append(coil_wire)

            # Create the lofted surface from the filament curves

            loft = Part.makeLoft(filament_curves, True)
            print("L:", loft.Length, "A:", loft.Area, "V", loft.Volume)
            loft_wire = loft_shape(
                filament_wires, solid=True, ruled=True, closed=True, label=coil_name
            )
            print(loft_wire)
            loft_part = PartWrapper(loft)
            print("lp", loft_part.shape.Volume)
            # filament_curves.append(loft_part)
            coil_curves.append(loft_wire)
            # There is something not right about the volume, I am not using makeloft correctly..
            # The volume is negative...and anyway we need to add filament thickness.
    return coil_curves


coil_curves = create_coil_from_nurbs(coil_filename)
# print(coil_curves)
# Show the CAD of the plasma surface and coils.
show_cad(coil_curves + [plasma_surface])
save_cad(
    coil_curves + [plasma_surface],
    "plasmastellarator.stp",
)
