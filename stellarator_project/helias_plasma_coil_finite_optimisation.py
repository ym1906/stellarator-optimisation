#!/usr/bin/env python

"""Script to import and optimize the HELIAS 5b coil configuration using SIMSOPT."""

import functools
import json
import operator
import os
from itertools import starmap

import h5py
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import minimize
from simsopt.field import (
    BiotSavart,
    Coil,
    Current,
    apply_symmetries_to_currents,
    apply_symmetries_to_curves,
    coils_via_symmetries,
)
from simsopt.geo import (
    CurveCurveDistance,
    CurveLength,
    CurveSurfaceDistance,
    CurveXYZFourier,
    LpCurveCurvature,
    MeanSquaredCurvature,
    SurfaceRZFourier,
    create_equally_spaced_curves,
    create_multifilament_grid,
    curves_to_vtk,
    plot,
)
from simsopt.objectives import QuadraticPenalty, SquaredFlux, Weight
from simsopt.solve import GPMO
from simsopt.util import (
    bluemira_nurbs_utils,
    in_github_actions,
)

# Configuration parameters
nphi = 64  # Number of points in the phi direction
ntheta = 64  # Number of points in the theta direction
input_name_nc3 = "stellarator_project/data/plasma/helias5_vmec.nc3"
coil_filename = "stellarator_project/data/magnets/HELIAS5_coils_fourier.h5"
coil_txt_names = "stellarator_project/data/magnets/HELIAS5_coils_fourier.txt"
ncoils = 6  # Number of unique coil shapes
R0 = 22.0  # Major radius for the initial circular coils
R1 = 12.2  # Minor radius for the initial circular coils
order = 9  # Number of Fourier modes for each Cartesian component of each coil

# Objective function weights
LENGTH_WEIGHT = Weight(0.1)
CC_THRESHOLD = 0.1
CC_WEIGHT = 1000
CS_THRESHOLD = 0.3
CS_WEIGHT = 10
CURVATURE_THRESHOLD = 5.0
CURVATURE_WEIGHT = 1e-6
MSC_THRESHOLD = 5
MSC_WEIGHT = 1e-6
MAXITER = 50 if in_github_actions else 51

OUT_DIR = "./output/"
os.makedirs(OUT_DIR, exist_ok=True)

# Functions to load and process coil data


def load_coils_from_hdf5(filename):
    """Load coil Fourier coefficients from an HDF5 file."""
    coils_data = []
    centers = []
    with h5py.File(filename, "r") as f:
        for coil_name in f:
            group = f[coil_name]
            centre = group["Centre"][:]
            cosine_xyz = group["Cosine_xyz"][:]
            sine_xyz = group["Sine_xyz"][:]
            order = cosine_xyz.shape[0] - 1

            # Create Fourier coefficients array
            fourier_coeffs = np.zeros((order + 1, 6))
            fourier_coeffs[:, 0] = sine_xyz[:, 0]
            fourier_coeffs[:, 1] = cosine_xyz[:, 0]
            fourier_coeffs[:, 2] = sine_xyz[:, 1]
            fourier_coeffs[:, 3] = cosine_xyz[:, 1]
            fourier_coeffs[:, 4] = sine_xyz[:, 2]
            fourier_coeffs[:, 5] = cosine_xyz[:, 2]

            coils_data.append(fourier_coeffs)
            centers.append(centre)

    return np.array(coils_data), np.array(centers), order


def h5_to_fourier_file_format(h5_filename, output_filename):
    """Convert HDF5 file containing Fourier coefficients into a specified format for CurveXYZFourier."""
    with h5py.File(h5_filename, "r") as f:
        coils_data = []
        max_order = 0
        for coil_name in f:
            group = f[coil_name]
            cosine_xyz = group["Cosine_xyz"][:]
            sine_xyz = group["Sine_xyz"][:]
            order = cosine_xyz.shape[0] - 1
            max_order = max(max_order, order)

            fourier_coeffs = np.zeros((order + 1, 6))
            fourier_coeffs[:, 0] = sine_xyz[:, 0]
            fourier_coeffs[:, 1] = cosine_xyz[:, 0]
            fourier_coeffs[:, 2] = sine_xyz[:, 1]
            fourier_coeffs[:, 3] = cosine_xyz[:, 1]
            fourier_coeffs[:, 4] = sine_xyz[:, 2]
            fourier_coeffs[:, 5] = cosine_xyz[:, 2]
            coils_data.append(fourier_coeffs)

        num_coils = len(coils_data)
        full_data = np.zeros((max_order + 1, 6 * num_coils))
        for ic, fourier_coeffs in enumerate(coils_data):
            full_data[: fourier_coeffs.shape[0], 6 * ic : 6 * (ic + 1)] = fourier_coeffs

        np.savetxt(output_filename, full_data, delimiter=",")


def load_curves_from_data(coil_data, centers, order=None, ppp=20):
    """Load Fourier coefficients data into CurveXYZFourier objects."""
    assert coil_data.shape[2] == 6
    num_coils = coil_data.shape[0]
    coils = [CurveXYZFourier(order * ppp, order) for _ in range(num_coils)]
    for ic in range(num_coils):
        dofs = coils[ic].dofs_matrix
        center = centers[ic]
        dofs[0][0] = center[0] + coil_data[ic, 0, 1]
        dofs[1][0] = center[1] + coil_data[ic, 0, 3]
        dofs[2][0] = center[2] + coil_data[ic, 0, 5]

        for io in range(min(order, coil_data.shape[1] - 1)):
            dofs[0][2 * io + 1] = coil_data[ic, io + 1, 0]
            dofs[0][2 * io + 2] = coil_data[ic, io + 1, 1]
            dofs[1][2 * io + 1] = coil_data[ic, io + 1, 2]
            dofs[1][2 * io + 2] = coil_data[ic, io + 1, 3]
            dofs[2][2 * io + 1] = coil_data[ic, io + 1, 4]
            dofs[2][2 * io + 2] = coil_data[ic, io + 1, 5]

        coils[ic].local_x = np.concatenate(dofs)
    return coils


# Initialize the boundary magnetic surface
s = SurfaceRZFourier.from_wout(
    input_name_nc3, range="full torus", nphi=nphi, ntheta=ntheta
)
# Create the initial coils
base_curves = create_equally_spaced_curves(
    ncoils, s.nfp, stellsym=True, R0=R0, R1=R1, order=order
)
base_currents = [Current(1.257e6) for i in range(ncoils)]
base_currents[0].fix_all()  # Fix the first current to avoid trivial solution

# use sum here to concatenate lists
numfilaments_n = 2  # number of filaments in normal direction
numfilaments_b = 3  # number of filaments in bi-normal direction
gapsize_n = 0.02  # gap between filaments in normal direction
gapsize_b = 0.04  # gap between filaments in bi-normal direction
rot_order = 1  # order of the Fourier expression for the rotation of the filament pack, i.e. maximum Fourier mode number
nfil = numfilaments_n * numfilaments_b
LENGTH_PEN = 1e-2

# Threshhold and weight for the coil-to-coil distance penalty in the objective function:
DIST_MIN = 0.1
DIST_PEN = 10
# Number of iterations to perform:
MAXITER = 50 if in_github_actions else 400
# base_curves_finite_build = sum(
#     [
#         create_multifilament_grid(
#             c,
#             numfilaments_n,
#             numfilaments_b,
#             gapsize_n,
#             gapsize_b,
#             rotation_order=rot_order,
#         )
#         for c in base_curves
#     ],
#     [],
# )
# Create a dictionary to hold filaments for each base coil
coil_filament_map = {}

# Create the filaments and group them by coil name
for i, base_curve in enumerate(base_curves):
    coil_name = f"coil_{i + 1}"
    filaments = create_multifilament_grid(
        base_curve,
        numfilaments_n,
        numfilaments_b,
        gapsize_n,
        gapsize_b,
        rotation_order=rot_order,
    )
    coil_filament_map[coil_name] = filaments
base_curves_finite_build = functools.reduce(
    operator.iadd, coil_filament_map.values(), []
)
base_currents_finite_build = functools.reduce(
    operator.iadd, [[c] * nfil for c in base_currents], []
)
curves_fb = apply_symmetries_to_curves(base_curves_finite_build, s.nfp, True)
currents_fb = apply_symmetries_to_currents(base_currents_finite_build, s.nfp, True)
coils_fb = list(starmap(Coil, zip(curves_fb, currents_fb)))
bs = BiotSavart(coils_fb)
bs.set_points(s.gamma().reshape((-1, 3)))


# coils = coils_via_symmetries(base_curves, base_currents, s.nfp, True)
# bs = BiotSavart(coils)
# bs.set_points(s.gamma().reshape((-1, 3)))
curves = [c.curve for c in coils_fb]
# for c in coils_fb:
#     print(c)
#     print(1)
# Define the individual terms of the objective function
# Define the objective function:
Jf = SquaredFlux(s, bs)
Jls = [CurveLength(c) for c in base_curves]
Jdist = CurveCurveDistance(curves, DIST_MIN)

# Form the total objective function
JF = (
    Jf
    + LENGTH_PEN
    * sum(QuadraticPenalty(Jls[i], Jls[i].J(), "max") for i in range(len(base_curves)))
    + DIST_PEN * Jdist
)


def fun(dofs):
    JF.x = dofs
    J = JF.J()
    grad = JF.dJ()
    cl_string = ", ".join([f"{J.J():.3f}" for J in Jls])
    mean_AbsB = np.mean(bs.AbsB())
    jf = Jf.J()
    kap_string = ", ".join(f"{np.max(c.kappa()):.1f}" for c in base_curves)
    print(
        f"J={J:.3e}, Jflux={jf:.3e}, sqrt(Jflux)/Mean(|B|)={np.sqrt(jf) / mean_AbsB:.3e}, CoilLengths=[{cl_string}], [{kap_string}], ||∇J||={np.linalg.norm(grad):.3e}"
    )
    return 1e-4 * J, 1e-4 * grad


# Perform a Taylor test
print("Performing a Taylor test")
f = fun
dofs = JF.x
np.random.seed(1)
h = np.random.uniform(size=dofs.shape)
J0, dJ0 = f(dofs)
dJh = sum(dJ0 * h)
for eps in [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]:
    J1, _ = f(dofs + eps * h)
    J2, _ = f(dofs - eps * h)
    print("err", (J1 - J2) / (2 * eps) - dJh)

# Run the optimization
print("Running the optimization")
res = minimize(
    fun,
    dofs,
    jac=True,
    method="L-BFGS-B",
    options={"maxiter": MAXITER, "maxcor": 400, "gtol": 1e-20, "ftol": 1e-20},
    tol=1e-20,
)


# Perform a second optimization with reduced penalty for coil length
dofs = res.x
LENGTH_WEIGHT *= 0.1
res = minimize(
    fun,
    dofs,
    jac=True,
    method="L-BFGS-B",
    options={"maxiter": MAXITER, "maxcor": 300},
    tol=1e-15,
)

# Save the optimized coil shapes and currents
bs.save(OUT_DIR + "biot_savart_opt.json")
bluemira_nurbs_utils.surface_to_nurbs(
    simsopt_surface=s,
    export_path="stellarator_project/data/plasma/finite_plasma_surface_nurbs_data.json",
    plot=True,
)
bluemira_nurbs_utils.curves_to_nurbs(
    curves=curves,
    export_path="stellarator_project/data/magnets/finite_magnets_nurbs_data.json",
    plot=False,
)

bluemira_nurbs_utils.filament_curves_to_nurbs(
    curves=curves,
    numfil=nfil,
    export_path="stellarator_project/data/magnets/finite_magnets_nurbs_data.json",
    plot=False,
)
