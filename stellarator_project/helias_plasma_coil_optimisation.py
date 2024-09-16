#!/usr/bin/env python

"""Script to import and optimize the HELIAS 5b coil configuration using SIMSOPT."""

import json
import os

import h5py
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import minimize
from simsopt.field import (
    BiotSavart,
    Current,
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
MAXITER = 50 if in_github_actions else 400

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


# Load coils from HDF5 file
coils_data, centers, order = load_coils_from_hdf5(coil_filename)
coils = load_curves_from_data(coil_data=coils_data, centers=centers, order=order)

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

coils = coils_via_symmetries(base_curves, base_currents, s.nfp, True)
bs = BiotSavart(coils)
bs.set_points(s.gamma().reshape((-1, 3)))
curves = [c.curve for c in coils]
curves_to_vtk(curves, OUT_DIR + "curves_init")
pointData = {
    "B_N": np.sum(bs.B().reshape((nphi, ntheta, 3)) * s.unitnormal(), axis=2)[:, :, None]
}
s.to_vtk(OUT_DIR + "surf_init", extra_data=pointData)

# Define the individual terms of the objective function
Jf = SquaredFlux(s, bs)
Jls = [CurveLength(c) for c in base_curves]
Jccdist = CurveCurveDistance(curves, CC_THRESHOLD, num_basecurves=ncoils)
Jcsdist = CurveSurfaceDistance(curves, s, CS_THRESHOLD)
Jcs = [LpCurveCurvature(c, 2, CURVATURE_THRESHOLD) for c in base_curves]
Jmscs = [MeanSquaredCurvature(c) for c in base_curves]

# Form the total objective function
JF = (
    Jf
    + LENGTH_WEIGHT * sum(Jls)
    + CC_WEIGHT * Jccdist
    + CS_WEIGHT * Jcsdist
    + CURVATURE_WEIGHT * sum(Jcs)
    + MSC_WEIGHT * sum(QuadraticPenalty(J, MSC_THRESHOLD, "max") for J in Jmscs)
)


def fun(dofs):
    """Wrapper function for the optimization."""
    JF.x = dofs
    J = JF.J()
    grad = JF.dJ()
    jf = Jf.J()
    BdotN = np.mean(
        np.abs(np.sum(bs.B().reshape((nphi, ntheta, 3)) * s.unitnormal(), axis=2))
    )
    outstr = f"J={J:.1e}, Jf={jf:.1e}, ⟨B·n⟩={BdotN:.1e}"
    cl_string = ", ".join([f"{J.J():.1f}" for J in Jls])
    kap_string = ", ".join(f"{np.max(c.kappa()):.1f}" for c in base_curves)
    msc_string = ", ".join(f"{J.J():.1f}" for J in Jmscs)
    outstr += f", Len=sum([{cl_string}])={sum(J.J() for J in Jls):.1f}, ϰ=[{kap_string}], ∫ϰ²/L=[{msc_string}]"
    outstr += f", C-C-Sep={Jccdist.shortest_distance():.2f}, C-S-Sep={Jcsdist.shortest_distance():.2f}"
    outstr += f", ║∇J║={np.linalg.norm(grad):.1e}"
    return J, grad


# Perform a Taylor test
print("Performing a Taylor test")
dofs = JF.x
np.random.seed(1)
h = np.random.uniform(size=dofs.shape)
J0, dJ0 = fun(dofs)
dJh = sum(dJ0 * h)
for eps in [1e-3, 1e-4, 1e-5, 1e-6, 1e-7]:
    J1, _ = fun(dofs + eps * h)
    J2, _ = fun(dofs - eps * h)
    print("err", (J1 - J2) / (2 * eps) - dJh)

# Run the optimization
print("Running the optimization")
res = minimize(
    fun,
    dofs,
    jac=True,
    method="L-BFGS-B",
    options={"maxiter": MAXITER, "maxcor": 300},
    tol=1e-15,
)
curves_to_vtk(curves, OUT_DIR + "curves_opt_short")
pointData = {
    "B_N": np.sum(bs.B().reshape((nphi, ntheta, 3)) * s.unitnormal(), axis=2)[:, :, None]
}
s.to_vtk(OUT_DIR + "surf_opt_short", extra_data=pointData)

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
curves_to_vtk(curves, OUT_DIR + "curves_opt_long")
pointData = {
    "B_N": np.sum(bs.B().reshape((nphi, ntheta, 3)) * s.unitnormal(), axis=2)[:, :, None]
}
s.to_vtk(OUT_DIR + "surf_opt_long", extra_data=pointData)

# Save the optimized coil shapes and currents
bs.save(OUT_DIR + "biot_savart_opt.json")
plot(coils, engine="plotly", close=True)
for c in coils:
    print(c.volume)
bluemira_nurbs_utils.surface_to_nurbs(
    simsopt_surface=s,
    export_path="/home/graeme/stellarator-project/stellarator_project/data/plasma/plasma_surface_nurbs_data.json",
    plot=False,
)
bluemira_nurbs_utils.curves_to_nurbs(
    curves=curves,
    export_path="/home/graeme/stellarator-project/stellarator_project/data/magnets/magnets_nurbs_data.json",
    plot=False,
)
print(bluemira_nurbs_utils.extract_normals(s, curve_points=curves[1].gamma()))
print(len(bluemira_nurbs_utils.extract_normals(s, curve_points=curves[1].gamma())))
print(curves[1].gamma())


def extract_curve_normals(surface, curves):
    all_data = {}

    for idx, curve in enumerate(curves):
        curve_points = curve.gamma()
        normals = bluemira_nurbs_utils.extract_normals(
            surface, curve_points=curve_points
        )
        # Convert numpy arrays to lists for JSON serialization
        all_data[f"curve_{idx}"] = {
            "curve_points": curve_points.tolist(),
            "normals": normals.tolist(),
        }

    return all_data


# Extracting the normals for all curves
all_normals_data = extract_curve_normals(s, curves)

# Define the output path
output_path = (
    "/home/graeme/stellarator-project/stellarator_project/data/magnets/normals_data.json"
)

# Save the data to a JSON file
with open(output_path, "w") as json_file:
    json.dump(all_normals_data, json_file, indent=4)

print(f"Normals data saved to {output_path}")
