"""Stellarator TF like shape"""

import numpy as np
from bluemira.geometry.tools import (
    make_polygon,
    make_circle,
    sweep_shape,
    force_wire_to_spline,
)
from bluemira.display import show_cad
from bluemira.geometry.wire import BluemiraWire


def wavey_coils(rad=1, repetitions=(0, 2, 4, 6, 8), size=0.4, discr=1000):
    circ = make_circle(radius=rad, axis=(0, 1, 0))
    coords = circ.discretise(discr)

    tf = {}
    shapes = []
    for i, rep in enumerate(repetitions):
        y = np.cos(np.linspace(0, rep * np.pi, discr)) * size

        coords._array[1] = y
        coords.translate((0, i, 0))
        wire = force_wire_to_spline(make_polygon(coords, closed=True))

        crossection = make_polygon(
            {"x": [-0.05, 0.05, 0.05, -0.05], "y": [-0.05, -0.05, 0.05, 0.05]},
            closed=True,
        )

        crossection.translate(wire.start_point().xyz)
        tf[rep] = sweep_shape(
            crossection, wire, transition=0, solid=False, frenet=False
        )

        shapes.extend([wire, crossection])

    show_cad(list(tf.values()))
    show_cad(shapes)

    return tf


def changing_cross_sections(rad=1):
    path1 = make_circle(radius=rad, start_angle=0, end_angle=90, axis=(0, 1, 0))
    path2 = make_circle(radius=rad, start_angle=90, end_angle=180, axis=(0, 1, 0))
    path3 = make_circle(radius=rad, start_angle=180, end_angle=270, axis=(0, 1, 0))
    path4 = make_circle(radius=rad, start_angle=270, end_angle=360, axis=(0, 1, 0))

    cross_section1 = make_polygon(
        {"x": [-0.05, 0.05, 0.05, -0.05], "y": [-0.05, -0.05, 0.05, 0.05]},
        closed=True,
    )

    cross_section2 = make_circle(radius=0.1)
    cross_section3 = cross_section1.deepcopy()
    cross_section4 = cross_section1.deepcopy()
    cross_section4.rotate(degree=45)

    for cs, pth in zip(
        (cross_section1, cross_section2, cross_section3, cross_section4),
        (path1, path2, path3, path4),
    ):
        cs.translate(pth.start_point().xyz)

    for cs in (cross_section2, cross_section4):
        cs.rotate(degree=90, base=cs.center_of_mass, direction=(0, 1, 0))

    shape = sweep_shape(
        [cross_section1, cross_section2, cross_section3, cross_section4],
        BluemiraWire([path1, path2, path3, path4]),
    )

    show_cad(
        [
            cross_section1,
            cross_section2,
            cross_section3,
            cross_section4,
            path1,
            path2,
            path3,
            path4,
        ]
    )
    show_cad(shape)


def weird_vacuum_vessel(rad=10, discr=1000, size=0.2):
    path1 = make_circle(radius=rad, start_angle=0, end_angle=20, axis=(1, 0, 0))
    path2 = make_circle(radius=rad, start_angle=20, end_angle=40, axis=(1, 0, 0))
    path3 = make_circle(radius=rad, start_angle=40, end_angle=60, axis=(1, 0, 0))

    dx = np.cos(np.linspace(0, 6 * np.pi, discr)) * size
    dz = np.cos(np.linspace(0, 12 * np.pi, discr)) * size

    cross_section1 = make_circle(radius=1)
    cross_section2 = make_circle(radius=1).discretise(discr)
    cross_section3 = make_circle(radius=1)
    cross_section4 = make_circle(radius=1)

    cross_section2._array[0] += dx
    cross_section2._array[1] += dz
    cross_section2 = force_wire_to_spline(make_polygon(cross_section2, closed=True))

    for cs, degree in zip(
        (cross_section2, cross_section3, cross_section4), (20, 40, 60)
    ):
        cs.rotate(degree=degree, direction=(1, 0, 0))

    for cs, pth in zip(
        (cross_section1, cross_section2, cross_section3), (path1, path2, path3)
    ):
        cs.translate(pth.start_point().xyz)

    cross_section4.translate(path3.end_point().xyz)

    show_cad(
        [
            cross_section1,
            cross_section2,
            cross_section3,
            cross_section4,
            path1,
            path2,
            path3,
        ]
    )

    shape = sweep_shape(
        [cross_section1, cross_section2, cross_section3, cross_section4],
        BluemiraWire([path1, path2, path3]),
        solid=False,
        frenet=False,
    )

    show_cad(shape)

    return shape, path1, path2, path3


if __name__ == "__main__":
    tf = wavey_coils(rad=1.5)
    changing_cross_sections()
    vv, *path = weird_vacuum_vessel()

    for cl, ang in zip((tf[0], tf[4], tf[6], tf[8]), (90, 110, 130, 150)):
        cl.rotate(degree=ang, base=cl.center_of_mass, direction=(1, 0, 0))

    tf[0].translate(path[0].start_point().xyz.flat - tf[0].center_of_mass)
    tf[4].translate(path[1].start_point().xyz.flat - tf[4].center_of_mass)
    tf[6].translate(path[2].start_point().xyz.flat - tf[6].center_of_mass)
    tf[8].translate(path[2].end_point().xyz.flat - tf[8].center_of_mass)

    show_cad([*tf.values(), vv])