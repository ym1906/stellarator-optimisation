# SPDX-FileCopyrightText: 2024-present {{ copyright-holder }} <{{ copyright-holder-email }}>
#
# SPDX-License-Identifier: MIT

"""Useful tools."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from typing import Any

import h5py
import matplotlib.pyplot as plt
import numpy as np
from bluemira.base.look_and_feel import bluemira_print, bluemira_print_flush
from geomdl import NURBS, helpers
from scipy.spatial import KDTree
from simsopt.geo import CurveXYZFourier


def read_json(file_path: str) -> dict[str, Any]:
    """Read JSON data from a file."""
    with open(file_path) as f:
        return json.load(f)


def taylor_test(fun, JF):  # noqa: N803
    """Perform a Taylor test."""
    bluemira_print("Performing a Taylor test")
    dofs = JF.x
    rng = np.random.default_rng(1)
    h = rng.uniform(size=dofs.shape)
    _J0, dJ0 = fun(dofs)  # noqa: N806
    dJh = sum(dJ0 * h)  # noqa: N806
    for eps in [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]:
        J1, _ = fun(dofs + eps * h)  # noqa: N806
        J2, _ = fun(dofs - eps * h)  # noqa: N806
        bluemira_print(f"err {(J1 - J2) / (2 * eps) - dJh}")


@dataclass
class NURBSData:
    """Nurbs surface data."""

    poles: list[list[float]]
    poles2d: np.ndarray
    num_ctlpts_u: int
    num_ctlpts_v: int
    knot_vector_u: list[float]
    knot_vector_v: list[float]
    internal_knot_vector_u: list[float]
    internal_knot_vector_v: list[float]
    degree_u: int
    degree_v: int
    weights: list[float]
    weights_reshaped: list[list[float]]
    mults_u: list[int]
    mults_v: list[int]

    @classmethod
    def from_surface(cls, surf: NURBS.Surface) -> NURBSData:
        """Extract NURBS data from a geomdl NURBS surface."""
        control_points = surf.ctrlpts
        control_points_size_u = surf.ctrlpts_size_u
        control_points_size_v = surf.ctrlpts_size_v
        control_points_2d = reshape_control_points(
            control_points, control_points_size_u, control_points_size_v
        )
        # poles = [pt[:3] for row in control_points for pt in row]  # Extract x, y, z
        poles = control_points
        poles2d = control_points_2d.tolist()
        degree_u = surf.degree_u
        degree_v = surf.degree_v
        knot_vector_u = surf.knotvector_u
        knot_vector_v = surf.knotvector_v
        internal_knot_vector_u = np.concatenate((
            [0],
            extract_internal_knots(knot_vector_u, degree_u),
            [1],
        )).tolist()
        internal_knot_vector_v = np.concatenate((
            [0],
            extract_internal_knots(knot_vector_v, degree_v),
            [1],
        )).tolist()
        mults_u = get_mults(knot_vector_u, degree_u)
        mults_v = get_mults(knot_vector_v, degree_v)

        # Ensures that the weights are correctly formatted and aligned with
        # the control points grid layout
        weights_reshaped = (
            np.array(surf.weights)
            .reshape(control_points_size_u, control_points_size_v)
            .tolist()
        )

        return cls(
            poles,
            poles2d,
            control_points_size_u,
            control_points_size_v,
            knot_vector_u,
            knot_vector_v,
            internal_knot_vector_u,
            internal_knot_vector_v,
            degree_u,
            degree_v,
            surf.weights,
            weights_reshaped,
            mults_u,
            mults_v,
        )

    def write_to_json(self, file_path):
        """Write NURBSData to a JSON file."""
        with open(file_path, "w") as f:
            json.dump(asdict(self), f)


def get_mults(knot_vector, degree):
    """Get multiplicities from knot vector."""
    mults = np.zeros(len(knot_vector))

    for i, knot in enumerate(knot_vector):
        mults[i] = helpers.find_multiplicity(knot, knot_vector)

    # I am not exactly sure why but it is necessary
    # to remove some of the boundary multiplicities,
    # this is requried by the bsplinesurface freecad (bluemira).
    mults = refine_multiplicities(mults, degree)
    # Convert NumPy arrays to Python lists
    return mults.tolist()


@dataclass
class NURBSCurveData:
    """Nurbs curve data."""

    poles: list[list[float]]
    degree: int
    knot_vector: list[float]
    internal_knot_vector: list[float]
    weights: list[float]
    mults: list[int]

    @classmethod
    def from_curve(cls, curve: NURBS.Curve) -> NURBSCurveData:
        """Extract NURBS curve data from a geomdl NURBS curve."""
        control_points = curve.ctrlpts
        poles = control_points
        weights = curve.weights
        degree = curve.degree
        knot_vector = curve.knotvector
        internal_knot_vector = np.concatenate((
            [0],
            extract_internal_knots(knot_vector, degree),
            [1],
        )).tolist()
        mults = np.zeros(len(knot_vector))
        for i, knot in enumerate(knot_vector):
            mults[i] = helpers.find_multiplicity(knot, knot_vector)
        mults = refine_multiplicities(mults, degree).tolist()

        return cls(poles, degree, knot_vector, internal_knot_vector, weights, mults)

    def write_to_json(self, file_path: str):
        """Write a list of NURBSCurveData to a JSON file."""
        with open(file_path, "w") as f:
            json.dump(asdict(self), f)


def generate_knot_vector(num_points: int, degree: int) -> np.ndarray:
    """Generate a knot vector for NURBS surface."""
    # Number of spans (segments between control points) is one less than
    # the number of control points.
    n = num_points - 1
    # Number of knots required for the knot vector
    num_knots = n + degree + 2
    # Initialise the knot vector with zeroes.
    knot_vector = np.zeros(num_knots)
    # Remove the first and last element.
    internal_knots = np.linspace(0, 1, n - degree + 2)[1:-1]
    # Set the internal knots.
    knot_vector[degree + 1 : n + 1] = internal_knots
    # Sets the last knots equal to one to ensure correct boundary conditions.
    knot_vector[n + 1 :] = 1
    return knot_vector


def make_periodic_knot_vector(num_points: int, degree: int) -> np.ndarray:
    """Generate a periodic knot vector."""
    # Number of spans (segments between control points) is one less than
    # the number of control points.
    n = num_points - 1
    # Number of knots required for the knot vector
    m = n + degree + 1
    # Initialise the knot vector with zeroes.
    knot_vector = np.zeros(m + 1)
    # Fill the knot vector.
    for i in range(m + 1):
        if i <= degree:
            knot_vector[i] = 0
        elif i >= m - degree:
            knot_vector[i] = 1
        else:
            knot_vector[i] = (i - degree) / (m - 2 * degree)
    return knot_vector


def extract_internal_knots(knot_vector: list[float], degree: int) -> np.ndarray:
    """Extract the internal knot vector from a given knot vector
    (excludes boundary knots).
    """
    # The number of repeated boundary knots is equal to the degree + 1
    p = degree
    # Remove the first p+1 elements and the last p+1 elements
    internal_knots = knot_vector[p + 1 : -p - 1]
    return np.array(internal_knots)


def reshape_weights(self) -> list[list[float]]:
    """Reshape weights to a 2D array (list of lists)."""
    num_rows = self.num_ctlpts_u
    num_cols = self.num_ctlpts_v
    return self.weights.reshape(num_rows, num_cols).tolist()


def extract_internal_multiplicities(knot_vector: list[float], degree: int) -> list[int]:
    """Extract the internal multiplicities from a given knot vector."""
    internal_knots = extract_internal_knots(knot_vector, degree)

    # Get the unique internal knots and their counts
    unique_internal_knots, counts = np.unique(internal_knots, return_counts=True)

    # Build the list of multiplicities in the same order as the internal knots
    return [counts[i] for i in range(len(unique_internal_knots))]


def refine_multiplicities(mults: np.ndarray, degree: int) -> np.ndarray:
    """Strip the first and last (degree + 1) multiplicities but retain one instance
    of the boundary values.

    Raises:
        ValueError: If multiplicites is shorter than twice the boundary count
    """
    boundary_count = degree + 1

    if len(mults) < 2 * boundary_count:
        raise ValueError(
            "The length of the multiplicities list is too short for the given degree."
        )

    # Retain one instance of the boundary values and the internal values
    return np.concatenate((
        [mults[0]],
        mults[boundary_count:-boundary_count],
        [mults[-1]],
    ))


def setup_nurbs_surface(
    points: np.ndarray, degree_u: int, degree_v: int, n_phi: int, n_theta: int
) -> NURBS.Surface:
    """Setup a NURBS surface with given parameters."""
    # Create a NURBS surface object
    surface = NURBS.Surface()

    # Set the degrees of the surface in the u and v directions
    surface.degree_u = degree_u
    surface.degree_v = degree_v

    # Check if the last dimension of points is 3 (indicating x, y, z coordinates)
    # If so, add a fourth dimension with a value of 1 (for homogeneous coordinates)
    if points.shape[-1] == 3:  # noqa: PLR2004
        points = np.concatenate(
            (points, np.ones((points.shape[0], points.shape[1], 1))), axis=-1
        )

    # Flatten the points array and convert it to a list
    control_points = points.reshape((-1, 4)).tolist()

    # Set the control points for the NURBS surface
    surface.set_ctrlpts(control_points, n_phi, n_theta)

    # Generate periodic knot vectors for the u and v directions
    knot_vector_u = make_periodic_knot_vector(n_phi, degree_u)
    knot_vector_v = make_periodic_knot_vector(n_theta, degree_v)

    # Set the knot vectors for the NURBS surface
    surface.knotvector_u = knot_vector_u.tolist()
    surface.knotvector_v = knot_vector_v.tolist()

    return surface


def setup_nurbs_curve(points: np.ndarray, degree: int) -> NURBS.Curve:
    """Setup a NURBS curve with given parameters."""
    curve = NURBS.Curve()
    curve.degree = degree
    if points.shape[-1] == 3:  # noqa: PLR2004
        points = np.concatenate((points, np.ones((points.shape[0], 1))), axis=-1)
    control_points = points.tolist()
    curve.set_ctrlpts(control_points)
    knot_vector = make_periodic_knot_vector(points.shape[0], degree)
    curve.knotvector = knot_vector.tolist()
    return curve


def reshape_control_points(control_points, num_ctrlpts_u, num_ctrlpts_v):
    """Reshape the control points of the NURBS surface.
    :param nurbs_data: NURBSData object containing the control points
                       and dimensions.
    :type nurbs_data: NURBSData
    :return: Reshaped 2D array of control points
    :rtype: np.ndarray.
    """
    # Reshape the control points into a 2D array
    return np.reshape(control_points, (num_ctrlpts_u, num_ctrlpts_v, 3))


def extract_points_from_simsopt_surface(simsopt_surface: Any) -> np.ndarray:
    """Extract points from the simsopt surface to model the NURBS surface on.
    The points are concatenated to make the surface fully closed.
    :param simsopt_surface: Simsopt surface object
    :type simsopt_surface: simsopt.surface
    :return: Matrix of points to model
    :rtype: np.ndarray.
    """
    points = simsopt_surface.gamma()
    points = np.concatenate((points, points[:, :1, :]), axis=1)
    return np.concatenate((points, points[:1, :, :]), axis=0)


def extract_points_from_simsopt_curve(simsopt_curve: Any) -> np.ndarray:
    """Extract points from the simsopt surface to model the NURBS surface on.
    The points are concatenated to make the surface fully closed.
    :param simsopt_surface: Simsopt surface object
    :type simsopt_surface: simsopt.surface
    :return: Matrix of points to model
    :rtype: np.ndarray.
    """
    points = simsopt_curve.gamma()
    return np.concatenate((points, points[:1, :]), axis=0)


def extract_normals(surface, curve_points):
    """Extract normal vectors from a surface at the closest points to the given
    curve points.

    Parameters:
    - surface: SurfaceRZFourier object representing the plasma surface.
    - curve_points:
        A numpy array of shape (n, 3) representing points on the curve.

    Returns:
    - normals:
        A numpy array of shape (n, 3) representing the normal vectors
        at the closest points on the surface.
    """
    # Reshape and flatten surface points
    surface_points = surface.gamma().reshape((-1, 3))

    # Build a KDTree for fast nearest neighbor search
    tree = KDTree(surface_points)

    # Find the nearest surface points for each curve point
    _, indices = tree.query(curve_points)

    # Get the normal vectors at the nearest surface points
    return surface.unitnormal().reshape((-1, 3))[indices]


def surface_to_nurbs(
    simsopt_surface: Any,
    export_path: str,
    *,
    return_data: bool = True,
    plot: bool = False,
):
    """Convert a simsopt surface to a NURBS so that it can be used in
    FreeCAD.

    :param simsopt_surface: Simsopt surface object
    :type simsopt_surface: simsopt.geo.surfacerzfourier
    :param export_path: Export path for json file
    :type export_path: str
    :param plot: Plot NURBS with plotly, defaults to False
    :type plot: bool, optional
    """
    points = extract_points_from_simsopt_surface(simsopt_surface)
    # Number of phi and theta points.
    n_phi, n_theta, _ = points.shape
    # Create the NURBS
    nurbs = setup_nurbs_surface(
        points=points, degree_u=3, degree_v=3, n_phi=n_phi, n_theta=n_theta
    )
    nurbs_data = NURBSData.from_surface(nurbs)
    nurbs_data.write_to_json(export_path)
    if plot:
        plot_nurbs(nurbs)
    if return_data:
        return nurbs_data
    return None


def curve_to_nurbs(simsopt_curve: Any, export_path: str, *, plot: bool = False) -> None:
    """Convert a simsopt curve to a NURBS so that it can be used in FreeCAD."""
    points = extract_points_from_simsopt_curve(simsopt_curve)
    nurbs = setup_nurbs_curve(points=points, degree=3)
    nurbs_curve_data = NURBSCurveData.from_curve(nurbs)
    nurbs_curve_data.write_to_json(export_path)
    if plot:
        plot_nurbs(nurbs)


def filament_curves_to_nurbs(
    curves: list[Any],
    numfil: int,
    export_path: str,
    *,
    return_data: bool = True,
    plot: bool = False,
):
    """Convert a list of simsopt curves to NURBS curves and export them as coils
    in a JSON file.

    Raises:
        ValueError: total number of curves not divisible by numfil
    """
    total_curves = len(curves)
    if total_curves % numfil != 0:
        raise ValueError(
            "The number of curves is not a multiple of the number of filaments per coil."
        )

    num_coils = total_curves // numfil

    bluemira_print(
        f"Total number of coils to process: {num_coils}\n"
        f"Total number of filaments per coil: {numfil}\n"
        f"Total number of curves provided: {total_curves}"
    )

    coils_data = {}

    for i in range(num_coils):
        coil = coils_data[f"coil_{i + 1}"] = {}

        bluemira_print_flush(f"Processing coil {i + 1}...")

        for j in range(numfil):
            index = i * numfil + j
            # print(f"  Processing filament {j + 1} for coil {i + 1} (index {index})...")

            points = extract_points_from_simsopt_curve(curves[index])
            nurbs = setup_nurbs_curve(points=points, degree=3)
            nurbs_curve_data = NURBSCurveData.from_curve(nurbs)

            coil[f"filament_{j + 1}"] = nurbs_curve_data

            if plot:
                plot_nurbs(nurbs)

    with open(export_path, "w") as f:
        json.dump(
            {
                c: {f: asdict(dat) for f, dat in fs.items()}
                for c, fs in coils_data.items()
            },
            f,
            indent=4,
        )

    bluemira_print(f"Data successfully written to {export_path}")

    if return_data:
        return coils_data
    return None


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
    """Convert HDF5 file containing Fourier coefficients into a specified
    format for CurveXYZFourier.
    """
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
    """Load Fourier coefficients data into CurveXYZFourier objects.

    Raises:
        ValueError: Coil shape cannot be handled
    """
    if coil_data.shape[2] == 6:  # noqa: PLR2004
        raise ValueError("Coil data shape is not handled")
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


def _mpl_3d_scat(data, title):
    _f, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.scatter(
        xs=data[:, 0],
        ys=data[:, 1],
        zs=data[:, 2],
    )
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title(title)
    plt.show()


def plot_nurbs(surf: NURBS.Surface | NURBS.Curve) -> None:
    """Plot the NURBS surface using Plotly."""
    surf.evaluate()
    evaluated_points = np.array(surf.evalpts)
    _mpl_3d_scat(evaluated_points, "NURBS Surface")
