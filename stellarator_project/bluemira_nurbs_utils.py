from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, List

import numpy as np
import plotly.graph_objects as go
from geomdl import NURBS, helpers
from scipy.spatial import KDTree

# Take a simsopt surface object, and converts it to a NURBS so that
# it is compatible with FreeCAD, which is used in Bluemira.


# A dataclass to hold the NURBS data.
@dataclass
class NURBSData:
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


@dataclass
class NURBSCurveData:
    poles: list[list[float]]
    degree: int
    knot_vector: list[float]
    internal_knot_vector: list[float]
    weights: list[float]
    mults: list[int]


def write_nurbs_data_to_json(nurbs_data: NURBSData, file_path: str) -> None:
    """Write NURBSData to a JSON file."""
    with open(file_path, "w") as f:
        json.dump(nurbs_data.__dict__, f)


def write_nurbs_curve_data_to_json(
    nurbs_curve_data: list[NURBSCurveData], file_path: str
) -> None:
    """Write a list of NURBSCurveData to a JSON file."""
    curves_data_dict = [curve.__dict__ for curve in nurbs_curve_data]
    with open(file_path, "w") as f:
        json.dump(curves_data_dict, f)


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
    if points.shape[-1] == 3:
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


def extract_nurbs_data(surf: NURBS.Surface) -> NURBSData:
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
    weights = surf.weights
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
    # Initialize mult array
    mults_u = np.zeros(len(knot_vector_u))
    # Fill mult array for uknotvector
    for i, knot in enumerate(knot_vector_u):
        mults_u[i] = helpers.find_multiplicity(knot, knot_vector_u)
    # I am not exactly sure why but it is necessary
    # to remove some of the boundary multiplicities,
    # this is requried by the bsplinesurface freecad (bluemira).
    mults_u = refine_multiplicities(mults_u, degree_u)

    # Initialize mult array
    mults_v = np.zeros(len(knot_vector_v))
    # Fill mult array for vknotvector
    for i, knot in enumerate(knot_vector_v):
        mults_v[i] = helpers.find_multiplicity(knot, knot_vector_v)
    mults_v = refine_multiplicities(mults_v, degree_v)
    # Convert NumPy arrays to Python lists
    mults_u = mults_u.tolist()
    mults_v = mults_v.tolist()
    # Ensures that the weights are correctly formatted and aligned with
    # the control points grid layout
    weights_reshaped = np.array(weights)
    weights_reshaped = weights_reshaped.reshape(
        control_points_size_u, control_points_size_v
    ).tolist()

    return NURBSData(
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
        weights,
        weights_reshaped,
        mults_u,
        mults_v,
    )


def extract_nurbs_curve_data(curve: NURBS.Curve) -> NURBSCurveData:
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

    return NURBSCurveData(
        poles,
        degree,
        knot_vector,
        internal_knot_vector,
        weights,
        mults,
    )


def plot_nurbs_surface(surf: NURBS.Surface) -> None:
    """Plot the NURBS surface using Plotly."""
    surf.evaluate()
    evaluated_points = np.array(surf.evalpts)
    fig = go.Figure(
        data=[
            go.Scatter3d(
                x=evaluated_points[:, 0],
                y=evaluated_points[:, 1],
                z=evaluated_points[:, 2],
                mode="markers",
                marker={"size": 2, "color": "blue"},
            )
        ]
    )
    fig.update_layout(
        scene={"xaxis_title": "X", "yaxis_title": "Y", "zaxis_title": "Z"},
        title="NURBS Surface",
    )
    fig.show()


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


def surface_to_nurbs(simsopt_surface: Any, export_path: str, plot: bool = False) -> None:
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
    nurbs_data = extract_nurbs_data(nurbs)
    write_nurbs_data_to_json(nurbs_data, export_path)
    if plot:
        plot_nurbs_surface(nurbs)


def setup_nurbs_curve(points: np.ndarray, degree: int) -> NURBS.Curve:
    """Setup a NURBS curve with given parameters."""
    curve = NURBS.Curve()
    curve.degree = degree
    if points.shape[-1] == 3:
        points = np.concatenate((points, np.ones((points.shape[0], 1))), axis=-1)
    control_points = points.tolist()
    curve.set_ctrlpts(control_points)
    knot_vector = make_periodic_knot_vector(points.shape[0], degree)
    curve.knotvector = knot_vector.tolist()
    return curve


def plot_nurbs_curve(curve: NURBS.Curve) -> None:
    """Plot the NURBS curve using Plotly."""
    curve.evaluate()
    evaluated_points = np.array(curve.evalpts)
    fig = go.Figure(
        data=[
            go.Scatter3d(
                x=evaluated_points[:, 0],
                y=evaluated_points[:, 1],
                z=evaluated_points[:, 2],
                mode="markers",
                marker={"size": 2, "color": "red"},
            )
        ]
    )
    fig.update_layout(
        scene={"xaxis_title": "X", "yaxis_title": "Y", "zaxis_title": "Z"},
        title="NURBS Curve",
    )
    fig.show()


def curve_to_nurbs(simsopt_curve: Any, export_path: str, plot: bool = False) -> None:
    """Convert a simsopt curve to a NURBS so that it can be used in FreeCAD."""
    points = extract_points_from_simsopt_curve(simsopt_curve)
    nurbs = setup_nurbs_curve(points=points, degree=3)
    nurbs_curve_data = extract_nurbs_curve_data(nurbs)
    write_nurbs_curve_data_to_json([nurbs_curve_data], export_path)
    if plot:
        plot_nurbs_curve(nurbs)


def curves_to_nurbs(curves: list[Any], export_path: str, plot: bool = False) -> None:
    """Convert a list of simsopt curves to NURBS curves and export them as a
    single JSON file.
    """
    nurbs_curve_data_list = []
    for curve in curves:
        points = extract_points_from_simsopt_curve(curve)
        nurbs = setup_nurbs_curve(points=points, degree=3)
        nurbs_curve_data = extract_nurbs_curve_data(nurbs)
        nurbs_curve_data_list.append(nurbs_curve_data)
        if plot:
            plot_nurbs_curve(nurbs)
    write_nurbs_curve_data_to_json(nurbs_curve_data_list, export_path)


def filament_curves_to_nurbs(
    curves: list[Any], numfil: int, export_path: str, plot: bool = False
) -> None:
    """Convert a list of simsopt curves to NURBS curves and export them as coils
    in a JSON file.
    """
    print(curves)
    total_curves = len(curves)
    if total_curves % numfil != 0:
        raise ValueError(
            "The number of curves is not a multiple of the number of filaments per coil."
        )

    num_coils = total_curves // numfil

    print(f"Total number of coils to process: {num_coils}")
    print(f"Total number of filaments per coil: {numfil}")
    print(f"Total number of curves provided: {total_curves}")

    coils_data = []

    for i in range(num_coils):
        coil_data = {f"coil_{i + 1}": {}}

        print(f"Processing coil {i + 1}...")

        for j in range(numfil):
            index = i * numfil + j
            # print(f"  Processing filament {j + 1} for coil {i + 1} (index {index})...")

            points = extract_points_from_simsopt_curve(curves[index])
            nurbs = setup_nurbs_curve(points=points, degree=3)
            nurbs_curve_data = extract_nurbs_curve_data(nurbs)

            coil_data[f"coil_{i + 1}"][f"filament_{j + 1}"] = nurbs_curve_data.__dict__

            if plot:
                plot_nurbs_curve(nurbs)

        coils_data.append(coil_data)
        print(f"Finished processing coil {i + 1}")

    with open(export_path, "w") as f:
        json.dump(coils_data, f, indent=4)

    print(f"Data successfully written to {export_path}")
