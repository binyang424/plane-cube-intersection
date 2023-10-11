import numpy as np

if __name__ == '__main__':

    from skspatial.objects import Plane, Points, Line

    # Define the plane.
    plane = Plane([0.031169303653388884, 0.3960435537039571, 0.16455387591698864],
                  [0.07248675268229973, 0.9210315202417607, 0.38268343236508984])

    # Define the line.
    linex1 = Line([1, 0, 0], [1, 0, 0])  # E01
    linex2 = Line([1, 1, 0], [1, 0, 0])  # E25
    linex3 = Line([1, 1, 1], [1, 0, 0])  # E67
    linex4 = Line([1, 0, 1], [1, 0, 0])  # E34

    liney1 = Line([0, 1, 0], [0, 1, 0])  # E02
    liney2 = Line([1, 1, 0], [0, 1, 0])  # E15
    liney3 = Line([1, 1, 1], [0, 1, 0])  # E47
    liney4 = Line([0, 1, 1], [0, 1, 0])  # E36

    linez1 = Line([0, 0, 1], [0, 0, 1])  # E03
    linez2 = Line([1, 0, 1], [0, 0, 1])  # E14
    linez3 = Line([1, 1, 1], [0, 0, 1])  # E57
    linez4 = Line([0, 1, 1], [0, 0, 1])  # E26

    # Find the intersection point.
    pointx1 = plane.intersect_line(linex1)  #
    pointx2 = plane.intersect_line(linex2)  #
    pointx3 = plane.intersect_line(linex3)  #
    pointx4 = plane.intersect_line(linex4)  #

    pointy1 = plane.intersect_line(liney1)  # [0, 0.75, 0]
    pointy2 = plane.intersect_line(liney2)  #
    pointy3 = plane.intersect_line(liney3)  #
    pointy4 = plane.intersect_line(liney4)  #

    pointz1 = plane.intersect_line(linez1)  #
    pointz2 = plane.intersect_line(linez2)  #
    pointz3 = plane.intersect_line(linez3)
    pointz4 = plane.intersect_line(linez4)  #

    pts = [pointx1, pointz2, pointy3, pointy2,
           pointy1, pointx2, pointz3, pointz4,
           pointz1, pointy4, pointx3, pointx4]

    for pt in pts:
        if any(pt < 0) or any(pt > 1):
            pass
        else:
            print("# Intersection point: ", pt)
            # print("-" * 20)

    # Intersection point:  [1.         0.         0.93422714]
    # Intersection point:  [1.         0.38816614 0.        ]
    # Intersection point:  [0.         0.46686784 0.        ]
    # Intersection point:  [0.         0.05137345 1.        ]
    # Intersection point:  [0.65276159 0.         1.        ]

    normal = [0.2198463103929542, 0.2620026302293849, 0.9396926207859084]
    point = [0.2198463103929542, 0.2620026302293849, 0.9396926207859084]
    plane = Plane(point, normal)
    projected = plane.project_point(point)

    p0 = [0, 0, 0]
    distance = plane.distance_point(p0)

    points = [[1., 0., 0.83022222],
              [1., 1., 0.55140484],
              [0., 1., 0.7853604],
              [0., 0.23017853, 1., ],
              [0.27431609, 0., 1.]]

    import pypolylib

    pypolylib.check_principal_axis(points)

    dl = np.linspace(0.1, 1.7, 7, endpoint=True)
    phi = np.linspace(0, 90, 41, endpoint=True)
    theta = np.linspace(90, 0, 41, endpoint=True)

    grid = np.array(np.meshgrid(phi, theta))

    grid = np.array(list(zip(grid[0].flatten(), grid[1].flatten())))
