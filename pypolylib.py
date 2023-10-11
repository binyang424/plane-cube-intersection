import numpy as np


def irregular_polygon_area(points, centroid=False):
    """
     Given the coordinates of the vertices, this routine applies the
     Shoelace formula for the area and centroid of an irregular 2d polygon.
     The vertices may be ordered clockwise or counterclockwise. If they are
     ordered clockwise, the area will be negative but correct in absolute value.

    Parameters
    ----------
    points : array-like
        The points of the polygon in an ordered manner. The points must be in a
        2D plane. The positive direction is counter-clockwise.
    centroid : bool, optional
        If True, the centroid of the polygon is returned as well.

    Returns
    -------
    area : float
        The signed area of the polygon.
    centroid : array-like, optional
        The centroid of the polygon. Only returned if `centroid` is True.

    Examples
    --------
    >>> points = [[0, 0], [1, 0], [1, 1], [0, 1]]
    >>> irregular_polygon_area(points, centroid=True)
    1.0, [0.5, 0.5]
    >>> irregular_polygon_area([[0, 0], [0, 1], [1, 2], [2, 1], [2, 0]], centroid=True, trianglar=False)
    -3.0, [1.0, 0.7777777777777778]

    References
    ----------
    1. [skspatial.measurement — scikit-spatial documentation]
        (https://scikit-spatial.readthedocs.io/en/stable/_modules/skspatial/measurement.html#area_signed)
    2. [An Elegant Proof of the Shoelace Method - YouTube]
        (https://www.youtube.com/watch?v=sNPh8jgngE0&ab_channel=ChenHongming)
    """
    points = np.asarray(points)

    shape = points.shape
    if len(shape) != 2:
        raise ValueError("The points must be 2D.")
    elif shape[0] == 2 and shape[1] != 2:
        points = points.T
    else:
        pass

    n_points = points.shape[0]

    if n_points < 3:
        raise ValueError("The polygon must have at least three points.")

    x = points[:, 0]
    y = points[:, 1]

    indices = np.arange(n_points)
    indices_offset = indices - 1

    # Area of triangles of two consecutive points and the origin.
    area_tri = 0.5 * (x[indices_offset] * y[indices] - x[indices] * y[indices_offset])
    area = np.sum(area_tri)

    if centroid:
        # Centroid of triangles of two consecutive points and the origin.
        centroid_tri_x = (1 / 3) * (x[indices_offset] + x[indices]) * area_tri
        centroid_tri_y = (1 / 3) * (y[indices_offset] + y[indices]) * area_tri

        cent = [np.sum(centroid_tri_x, axis=0) / area, np.sum(centroid_tri_y, axis=0) / area]

        return area, cent

    return area


def polygon_cross_product(a, b):
    """
    Note : Just use np.cross() instead.

    This routine returns the cross product of the vectors (a x b). The cross
    product is defined as the vector perpendicular to both a and b, with a
    direction given by the right hand rule and a magnitude equal to the area
    of the parallelogram that the vectors span.

    Parameters
    ----------
    a : array-like
        The first point.
    b : array-like
        The second point.

    Returns
    -------
    c : array-like
        The cross product of the vectors from a to b and from a to c.
    """
    a = np.asarray(a)
    b = np.asarray(b)

    if a.shape != b.shape:
        raise ValueError("The points must have the same shape.")

    c = np.cross(a, b)

    return c


def get_polygon_normal(v1, v2, v3):
    """
    This routine returns the normal vector of a polygon by three of its vertices.

    Parameters
    ----------
    v1 : array-like
        The first vertex.
    v2 : array-like
        The second vertex.
    v3 : array-like
        The third vertex.

    Returns
    -------
    norm : array-like
        The normal vector of the polygon.
    """

    v1 = np.asarray(v1)
    v2 = np.asarray(v2)
    v3 = np.asarray(v3)

    if v1.shape != v2.shape or v1.shape != v3.shape:
        raise ValueError("The points must have the same shape.")

    v1v2 = v2 - v1
    v1v3 = v3 - v1

    norm = np.cross(v1v2, v1v3)
    norm = norm / np.linalg.norm(norm)

    return norm


def check_principal_axis(vertices, tol=1.0e-12):
    """
    This routine checks if the polygon is perpendicular to a principal axis.

    Parameters
    ----------
    vertices : array-like
        The points of the polygon. the shape of the array must be (n, 3).
        n is the number of points.

    tol : float, optional
        The tolerance for checking if the polygon is perpendicular to a principal axis.

    Returns
    -------
    axis : int
        The principal axis that the polygon is perpendicular to. 0, 1, 2 and -1. 0 for x-axis,
        1 for y-axis, 2 for z-axis and -1 for not perpendicular to any principal axis.
    """

    if not isinstance(vertices, np.ndarray):
        vertices = np.asarray(vertices, dtype=float)

    if vertices.shape[1] != 3:
        raise ValueError("The points must be 3D.")

    if np.all(np.abs(vertices[:, 0] - vertices[0, 0]) < tol):
        return 0
    elif np.all(np.abs(vertices[:, 1] - vertices[0, 1]) < tol):
        return 1
    elif np.all(np.abs(vertices[:, 2] - vertices[0, 2]) < tol):
        return 2
    else:
        return -1


def get_polygon_local_coord(points, sort=True, **kwargs):
    """
    This routine returns the local coordinate system of a 3D planar polygon. The
    unit vectors of the local coordinate system are first defined by the first
    three points of the polygon. The x-axis is defined by the vector from the
    first point to the second point. The z-axis is defined by the cross product
    of the x-axis and the vector from the second point to the third point. The
    y-axis is defined by the cross product of the z-axis and the x-axis. The
    points of the polygon are then projected onto the local coordinate system.

    Parameters
    ----------
    points : array-like
        The points of the polygon in an ordered manner. # TODO : The points can be any order.
    sort : bool, optional
        If True, the points are sorted according to the angle between the x-axis
        and the points in the local coordinate system.

    Returns
    -------
    pts_c : array-like
        The points of the polygon in the local coordinate system.

    Examples
    --------
    >>> points = [[0, 1, 0], [1, 0, 0], [1, 1, 1]]
    >>> get_polygon_local_coord(points)
    array([[-0.70710678,  0.40824829, -0.57735027],
       [ 0.70710678,  0.40824829, -0.57735027],
       [ 0.        ,  1.63299316, -0.57735027]])
    """
    if not isinstance(points, np.ndarray):
        points = np.asarray(points, dtype=np.double)

    # axis = check_principal_axis(points)
    axis = kwargs.get("axis", check_principal_axis(points))

    if axis == 2:
        pts_c = points - points[0]
    elif axis == 0:
        pts_c = points[:, [1, 2, 0]]
        pts_c = pts_c - pts_c[0]
    elif axis == 1:
        pts_c = points[:, [0, 2, 1]]
        pts_c = pts_c - pts_c[0]
    else:
        pts_c = points - points[0]

        i_prime = pts_c[1, :]
        k_prime = np.cross(pts_c[1, :], pts_c[2, :])
        j_prime = np.cross(k_prime, i_prime)

        i_prime = i_prime / np.linalg.norm(i_prime)
        j_prime = j_prime / np.linalg.norm(j_prime)
        k_prime = k_prime / np.linalg.norm(k_prime)

        pts_c[:, 0] = np.dot(points, i_prime)
        pts_c[:, 1] = np.dot(points, j_prime)
        pts_c[:, 2] = np.dot(points, k_prime)

    pts_c = pts_c - np.average(pts_c, axis=0)

    if sort:
        # calculate the angle between the x-axis and the points
        theta = np.arctan2(pts_c[:, 1], pts_c[:, 0])
        # sort the points according to the angle
        idx = np.argsort(theta)
        pts_c = pts_c[idx, :]

        return pts_c[:, : 2], idx

    return pts_c[:, :2]


def pyramid_height(vertices, **kwargs):
    """
    This routine returns the height of a pyramid given the coordinates of its
    vertices (four or more points with the last point being the apex).

    Parameters
    ----------
    vertices : array-like
        The vertices of the pyramid.

    Returns
    -------
    h : float
        The height of the pyramid.

    Examples
    --------
    >>> vertices = [[0, 0, 0], [2, 0, 0], [1, 1.732, 0], [1, 1.732 / 3, 1.732]]
    >>> pyramid_height(vertices)
    1.732
    """
    if not isinstance(vertices, np.ndarray):
        vertices = np.asarray(vertices)

    if vertices.shape[0] < 3:
        raise ValueError("The polygon must have at least three points.")

    if vertices.shape[1] != 3:
        raise ValueError("The point coordinates must be (x, y, z).")

    # axis = check_principal_axis(vertices[: -1, :])
    axis = kwargs.get("axis", check_principal_axis(vertices[: -1, :]))

    if axis == 2:
        h = vertices[-1, 2] - vertices[0, 2]
    elif axis == 0:
        h = vertices[-1, 0] - vertices[0, 0]
    elif axis == 1:
        h = vertices[-1, 1] - vertices[0, 1]
    else:
        v = vertices - vertices[0, :]

        norm = np.cross(v[1, :], v[2, :])
        norm = norm / np.linalg.norm(norm)

        h = np.dot(vertices[-1, :] - vertices[0, :], norm)

    return abs(h)


def pyramid_volume(vertices):
    """
    This routine returns the volume of a pyramid given the coordinates of its
    vertices (four or more points with the last point being the apex).

    Parameters
    ----------
    points : array-like
        The vertices of the pyramid.

    Returns
    -------
    volume : float
        The volume of the pyramid.

    Examples
    --------
    >>> vertices = [[0, 0, 0], [2, 0, 0], [1, 1.732051, 0], [1, 1.732051 / 3, 1.732051]]
    >>> pyramid_volume(vertices)
    1.0000002222003332
    """
    third = 1.0 / 3.0

    local_coord = get_polygon_local_coord(vertices[:-1, :], sort=True)
    area = abs(irregular_polygon_area(local_coord[:, :-1]))

    height = abs(pyramid_height(vertices))

    volume = area * height * third

    return volume


def convex_polyhedron_volume(xyz_vertex, fvertices, nface_vertices, face_mask=None, method="apex", write_stl=False):
    """
    This routine returns the volume of a convex polyhedron given the coordinates
    of its vertices and the vertices of each face. The faces are assumed to be
    convex polygons.

    Parameters
    ----------
    xyz_vertex : array-like
        The coordinates of the vertices of the polyhedron. The columns correspond
        to the x, y and z coordinates.
    fvertices : array-like
        The vertices of each face of the polyhedron. The column index corresponds
        to the face number and the value of each column is the vertex number corresponding
        to the vertex coordinate in `xyz_vertex`.
    face_mask : array-like, optional
        A convenience array for turning faces on or off. By default, all faces are
        considered as part of the polyhedron surface.

    Returns
    -------
    volume : float
        The volume of the polyhedron.
    """
    area = []
    height = []

    # f_orig = fvertices.copy()
    # mask_orig = face_mask.copy()

    if method == "apex":
        apex = xyz_vertex[0, :]
        face_mask[[0, 3, 4]] = 0
    elif method == "center":
        apex = np.average(xyz_vertex, axis=0)

    if face_mask is not None:
        fvertices = fvertices[:, face_mask == 1]
        nface_vertices = nface_vertices[face_mask == 1]

    nfaces = fvertices.shape[1]
    for i in np.arange(nfaces):
        vertices = fvertices[:nface_vertices[i], i]
        coord = xyz_vertex[vertices, :]

        axis = check_principal_axis(coord)
        xyz_local, idx = get_polygon_local_coord(coord, sort=True, axis=axis)

        # update the order of vertices in counter-clockwise direction
        # Note : This is not necessary for computing the volume.
        fvertices[:nface_vertices[i], i] = vertices[idx]

        area.append(irregular_polygon_area(xyz_local))

        height.append(pyramid_height(np.vstack([coord[idx, :], apex]), axis=axis))

    volume = np.sum(np.dot(area, height)) / 3
    if write_stl:
        write_polyhedron_vtk(xyz_vertex, fvertices, nface_vertices, filename="polyhedron.vtk")

    # print("The volume of the polyhedron by pypolylib is {}.".format(volume))

    return volume


def cartisan2spherical(x, y, z, degrees=False):
    """
    Converts Cartesian coordinates to spherical coordinates.

    Parameters
    ----------
    x, y, z : float
        The Cartesian coordinates.
    degrees : bool
        If True, the angles are returned in degrees.

    Returns
    -------
    r : float
        The radius.
    theta : float
        The angle between the projection of the radius vector onto the xy-plane and the x-axis.
    phi : float
        The angle between the radius vector and the z-axis.
    """
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    theta = np.arctan2(y, x)
    phi = np.arccos(z / r)

    if degrees:
        theta = np.degrees(theta)
        phi = np.degrees(phi)

    return r, theta, phi


def write_polyhedron_vtk(xyz_vertex, fvertices, nface_vertices, filename="polyhedron.vtk"):
    """
    This routine writes a polyhedron to a VTK file.

    Parameters
    ----------
    xyz_vertex : array-like
        The coordinates of the vertices of the polyhedron. The columns correspond
        to the x, y and z coordinates.
    fvertices : array-like
        The vertices of each face of the polyhedron. The column index corresponds
        to the face number and the value of each column is the vertex number corresponding
        to the vertex coordinate in `xyz_vertex`.
    face_mask : array-like, optional
        A convenience array for turning faces on or off. By default, all faces are
        considered as part of the polyhedron surface.
    filename : str, optional
        The name of the VTK file.

    Returns
    -------
    None
    """
    header = "# vtk DataFile Version 5.1\nvtk output\nASCII\nDATASET POLYDATA\n"
    pts = "POINTS {} double\n".format(xyz_vertex.shape[0])

    # write to a text file
    with open(filename, "w") as f:
        f.write(header)
        f.write(pts)
        for i in np.arange(xyz_vertex.shape[0]):
            f.write("{} {} {}\n".format(xyz_vertex[i, 0], xyz_vertex[i, 1], xyz_vertex[i, 2]))

        f.write("POLYGONS {} {}\n".format(fvertices.shape[1] + 1, np.sum(nface_vertices)))
        f.write("OFFSETS vtktypeint64\n")
        offset = 0
        for i in nface_vertices:
            f.write("{} ".format(offset))
            offset += i
        f.write("{}\n".format(offset))
        f.write("CONNECTIVITY vtktypeint64\n")

        connectivity = [fvertices[:nface_vertices[i], i] for i in np.arange(fvertices.shape[1])]
        for i in connectivity:
            for j in i:
                f.write(" {}".format(j))

        # save the file
        f.close()

        import pyvista as pv
        mesh = pv.read(filename)
        mesh.triangulate(inplace=True)
        mesh.save(filename[:-4] + ".stl", binary=True)
        # mesh.plot(show_edges=True, color="w", line_width=1.0)
        # delete vtk file
        import os
        os.remove(filename)
        print("The volume of the polyhedron by PyVista is {}.".format(mesh.volume))


def spherical2cartisan(r, theta, phi, degrees=False):
    """
    Converts spherical coordinates to Cartesian coordinates.

    Parameters
    ----------
    r : float
        The radius.
    theta : float
        The angle between the projection of the radius vector onto the xy-plane and the x-axis.
        theta takes value from 0 to 360 degrees or 0 to 2*pi.
    phi : float
        The angle between the radius vector and the z-axis. phi takes value from 0 to 180 degrees or 0 to pi.
    degrees : bool
        If True, the input angles are in degrees.

    Returns
    -------
    x, y, z : float
        The Cartesian coordinates.
    """
    if degrees:
        theta = np.radians(theta)
        phi = np.radians(phi)

    x = r * np.cos(theta) * np.sin(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(phi)

    return x, y, z


def plane_spherical(r, theta, phi, degrees=False):
    """
    Computes the coefficients of the hessian normal form of a plane.

    Parameters
    ----------
    r : float
        The vector radius of the plane. It is the distance from the origin to the plane.
    theta : float
        The angle between the projection of the radius vector onto the xy-plane and the x-axis.
        theta takes value from 0 to 360 degrees or 0 to 2*pi.
    phi : float
        The angle between the radius vector and the z-axis. phi takes value from 0 to 180 degrees or 0 to pi.
    degrees : bool
        If True, the input angles are in degrees.

    Returns
    -------
    a, b, c, d : float
        The coefficients of the hessian normal form of a plane.
    """
    if degrees:
        theta = np.radians(theta)
        phi = np.radians(phi)

    # n = (a, b, c) is the unit normal vector of the plane
    a = np.cos(theta) * np.sin(phi)
    b = np.sin(theta) * np.sin(phi)
    c = np.cos(phi)
    d = r

    return a, b, c, d


def plane_3points(x1, y1, z1, x2, y2, z2, x3, y3, z3):
    """
    Computes the coefficients of the hessian normal form of a plane.
    The plane is defined by three points.
    """
    # TODO : Check if the three points are collinear.

    # vector 1
    v1 = np.array([x2 - x1, y2 - y1, z2 - z1])
    # vector 2
    v2 = np.array([x3 - x1, y3 - y1, z3 - z1])

    # n = (a, b, c) is the unit normal vector of the plane
    n = np.cross(v1, v2)
    a, b, c = n / np.linalg.norm(n)

    d = a * x1 + b * y1 + c * z1

    return a, b, c, d


def plane_plot(a, b, c, d):
    # Define the range of x and y values
    x = np.linspace(-0, 1, 100)
    y = np.linspace(-0, 1, 100)

    # Generate a grid of (x, y) coordinates
    X, Y = np.meshgrid(x, y)
    Z = (d - a * X - b * Y) / c

    return X, Y, Z


def face_number(faces, point_number):
    """
    Given the label of a point, find the face number that the point belongs to.

    Parameters
    ----------
    faces : dict
        The dictionary stores the face label as the key and the label list of the points
        it consists as the value.
    point_number : int
        The label of the point.

    Returns
    -------
    face_number : list of int
        The label of the face that the point belongs to.
    """
    f = [key for key, value in faces.items() if point_number in value]

    return f


def add_vertex_to_face(iface, nfaces, iface_mask, face, nface_vertices, nvert, ivertices):
    """
    adds a face for the plane - cube intersection

    Parameters
    ----------
    iface : int
        The face number specified by cube_paths.pdf.
    nfaces : int
        The number of faces with intersection, a maximum of 7 - six for the cube plus one for the plane.
    iface_mask : array-like
        A convenience array for turning faces on or off. len(face_mask) = 7.
    face : array-like
        The face number of the polyhedra. len(face) = 7.
    nface_vertices : array-like
        How many vertices are with each face. len(nface_vertices) = 7.
    nvert : int
        The number of vertices of the intersected polyhedra in the range of [4, 14],
        whcih includes both cube vertices(8) and intersection points (6 at maximum).
    ivertices : array-like
        The vertices of each face with a shape of (7, 7).

    Returns
    -------
    fvertices : list of int
        The vertices of each face with a shape of (7, 7). The input fvertices is updated.
    """
    if iface_mask[iface] == 0:
        iface_mask[iface] = 1
        nfaces += 1
        face[nfaces] = iface

    nface_vertices[iface] += 1

    ivertices[nface_vertices[iface], iface] = nvert

    return iface_mask, nfaces, face, nface_vertices, ivertices


def check_if_vertex_exists(xyz_intersect, xxyyzz):
    """ √
    checks if a vertex already exists in the list of vertices

    Parameters
    ----------
    xyz_intersect : array-like
        The coordinates of the existing intersection points. shape = [m, 3].
        m is the number of existing intersection points.
    xxyyzz : float
        The x, y and z coordinates of the point to be tested.

    Returns
    -------
    idum : int
        A dummy variable. idum = 1 if the vertex already exists, otherwise idum = 0.
    """
    tiny = 1.0e-12
    # check if xzy_intersect is np.array
    if not isinstance(xyz_intersect, np.ndarray):
        xyz_intersect = np.asarray(xyz_intersect)
    if not isinstance(xxyyzz, np.ndarray):
        xxyyzz = np.array(xxyyzz)

    if np.sum(np.abs(xyz_intersect - xxyyzz), axis=1).min() <= tiny:
        return True
    else:
        return False


def check_edge_intersection(points1, points2, a, b, c, d):
    """ √
    Computes the distance lambda_ along an edge between two vertices where an intersection occurs.

    TODO : Check if there is an intersection at all.

    Parameters
    ----------
    points1 : arraly-like
        The coordinates of the first vertex of the edge.
    points2 : arraly-like
        The coordinates of the second vertex of the edge.
    a, b, c, d : float
        The coefficients of the hessian normal plane ax + by + cz + d = 0.

    Returns
    -------
    lambda_ : float
        The distance along the edge for an intersection. lambda_ will be between 0 and 1
        in the unit cube for a valid intersection.
    xx, yy, zz : float
        The coordinates of the intersection point.
    """
    # TODO : Check if the line is on the plane.

    tiny = 1.0e-12
    eij = np.array(points2) - np.array(points1)
    n = [a, b, c]
    # numerator = d - (a * points1[0] + b * points1[1] + c * points1[2])
    numerator = d - np.dot(n, points1)
    # denominator = a * eij[0] + b * eij[1] + c * eij[2]
    denominator = np.dot(eij, n)

    if denominator != 0.0:
        lambda_ = numerator / denominator
    elif abs(numerator) < tiny and denominator == 0.0:
        lambda_ = 0.0
    elif denominator == 0.0:
        lambda_ = 1000  # for a line parallel to the plane, set lambda_ to a large number as a flag

    if abs(lambda_) < 1.0e-12:
        lambda_ = 0.0
    elif abs(lambda_ - 1.0) < 1.0e-12:
        lambda_ = 1.0

    if 0.0 <= lambda_ <= 1.0:
        xxyyzz = points1 + lambda_ * eij

        return lambda_, xxyyzz
    else:
        return lambda_, None
