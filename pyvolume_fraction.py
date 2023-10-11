import time
import numpy as np
import pypolylib as ppl  # pypolylib.py


def volume_fraction(dl, theta, phi, vxyz, paths, faces_cube):
    """
    Calculates the volume of the polyhedron resulting from the intersection of a plane and a cube.

    Parameters
    ----------
    dl : float
        The length of the radius vector. The radius vector is the distance from the origin to
        the cutting.
    theta : float
           The angle between the radius vector and the x-axis. The angle is measured in degrees.
    phi : float
            The angle between the radius vector and the z-axis. The angle is measured in degrees.
    vxyz : ndarray
        The coordinates of the vertices of the cube.
    paths : dict
        The key is the path number, the value is the path described by the vertex numbers.
    faces_cube : dict
        The key is the vertex number, the value is the faces that attached to the vertex.

    Returns
    -------
    volume : float
        The volume of the polyhedron.
    """
    """ Intersection plane """
    a, b, c, d = ppl.plane_spherical(dl, theta, phi, degrees=True)

    """ Data initialization """
    # vertices
    # 1. the vertices of the cut plane
    xyz_intersect = np.zeros((6, 3), dtype=np.double)
    intersect_mask = np.zeros(6, dtype=int)
    # The number of intersection points; only values of 3, 4, 5, and 6 are possible.

    # 2. the vertices of the intersection polyhedra. 14 vertices at most:
    #    8 vertices of the cube + 6 vertices of the cut plane
    xyz_vertex = np.zeros([14, 3], dtype=np.double)
    ipath = np.full(14, -1, dtype=int)

    # faces
    # number of faces with intersection, a maximum of 7 - six for the cube plus one for the plane
    nfaces = 0
    faces = np.full(7, -1, dtype=int)  # the face number

    nface_vertices = np.zeros(7, dtype=int)  # how many vertices are with each face
    fvertices = np.full([7, 7], -1, dtype=int)  # the vertices of each face

    # a convenience array for turning faces on or off
    face_mask = np.zeros(7, dtype=int)

    """ The intersection test starts with V0 """
    # update vertex information
    xyz_vertex[0, :] = vxyz[0, :]
    nvertex = 0  # number of vertices of the intersection polyhedra (now V0)
    ipath[0] = 0  # the path number of V0

    # update face information
    f = ppl.face_number(faces_cube, 0)
    faces[: len(f)] = f
    face_mask[f] = 1
    nfaces += len(f)
    nface_vertices[f] += 1  # number of vertices of each face
    fvertices[nface_vertices[f] - 1, f] = nvertex

    """ Traverse all the paths """
    for key, path in paths.items():
        for iedge, jedge in zip(path, path[1:]):
            # edge of the path
            points1 = vxyz[iedge, :]
            points2 = vxyz[jedge, :]

            # Intersection test
            lambda_, xxyyzz = ppl.check_edge_intersection(points1, points2, a, b, c, d)
            if xxyyzz is not None:
                # print("lambda_ = {} and iedge = {}".format(lambda_, iedge))

                # check if the intersection point already exists in xyz_vertex
                # if not, add it to xyz_vertex
                if not ppl.check_if_vertex_exists(xyz_vertex, xxyyzz):

                    nvertex += 1

                    # update vertex information
                    xyz_vertex[nvertex, :] = xxyyzz
                    ipath[nvertex] = key

                    # update xyz_intersect
                    xyz_intersect[key, :] = xxyyzz
                    intersect_mask[key] = 1

                    # update face information
                    face_eff = set(ppl.face_number(faces_cube, iedge)).intersection(
                        ppl.face_number(faces_cube, jedge))
                    face_eff = np.array(list(face_eff))

                    mask = face_mask[face_eff] == 0
                    new_faces = face_eff[mask]

                    # print("new_faces = {}".format(new_faces))

                    if len(new_faces) > 0:
                        faces[nfaces: nfaces + len(new_faces)] = new_faces
                        nfaces += len(new_faces)
                        face_mask[new_faces] = 1

                    face_eff = np.append(face_eff, 6)

                    nface_vertices[face_eff] += 1
                    fvertices[nface_vertices[face_eff] - 1, face_eff] = nvertex

                break

            # for the dotted line, only the first case could lead to an intersection. So this
            # if statement is necessary to avoid false intersection detection because the lambda_
            # value is not checked here.
            elif len(path) > 2 and lambda_ > 1:
                if not ppl.check_if_vertex_exists(xyz_vertex, points2):
                    # update vertex information
                    nvertex += 1
                    xyz_vertex[nvertex, :] = points2
                    ipath[nvertex] = key

                    # update face information
                    face_eff = np.array(ppl.face_number(faces_cube, jedge))
                    mask = face_mask[face_eff] == 0
                    new_faces = face_eff[mask]

                    if len(new_faces) > 0:
                        faces[nfaces: nfaces + len(new_faces)] = new_faces
                        nfaces += len(new_faces)
                        face_mask[new_faces] = 1

                        # print("2. new_faces = {}".format(new_faces))

                    nface_vertices[face_eff] += 1
                    fvertices[nface_vertices[face_eff] - 1, face_eff] = nvertex
    # set up the intersection plane

    nintersect = np.sum(intersect_mask)
    if nintersect >= 3:

        faces[-1] = 6  # face 6 is the intersection plane
        face_mask[-1] = 1
        nfaces += 1
        fvertices = fvertices

        """ Data summary """
        xyz_intersect = xyz_intersect[intersect_mask == 1, :]
        nface_vertices[-1] = xyz_intersect.shape[0]
        xyz_vertex = xyz_vertex[:nvertex + 1, :]

        """ Volume of the intersection polyhedra """
        volume = ppl.convex_polyhedron_volume(xyz_vertex, fvertices, nface_vertices, face_mask, method="center")

        return volume, xyz_intersect.shape[0]

    elif nintersect == 0:
        return -1.0, 0


def volume_plane_cube_intersection(ndl, ntheta, nphi, vxyz, paths, faces_cube):
    """  √
    Exercises the routines for calculating the volume of the polyhedron
    resulting from the intersection of a plane and a cube.

    Parameters
    ----------
    ndl : int
        The number of dl, one movie frame per dl.
    ntheta : int
        The number of theta, one movie frame per theta.
    nphi : int
        The number of phi, one movie frame per phi.
    """

    """ loop limits for dl """
    # dl runs from to sqrt(dx^2 + dy^2 + dz^2) = sqrt(3), theta from 0 to 2*pi, phi from 0 to pi.
    dl_hi = 1.7  # high
    dl_lo = 0.1  # low
    dl_steps = np.linspace(dl_lo, dl_hi, ndl)

    """ loop limits for theta """
    theta_hi = 0.  # in degrees
    theta_lo = 90.0  # in degrees
    theta_steps = np.linspace(theta_lo, theta_hi, ntheta)

    """ loop limits for phi """
    phi_hi = 90  # in degrees
    phi_lo = 0  # in degrees
    phi_steps = np.linspace(phi_lo, phi_hi, nphi)

    """ loop over ray length dl - the radius vector """
    for k, dl in enumerate(dl_steps):
        print('working on dl', dl)

        filename = 'volume_' + str(k).zfill(2) + '.dat'

        with open(filename, 'w') as file:
            file.write(str(dl) + '\n')  # write dl to file in the first line

            """ loop over phi """
            # TODO : meshgrid(phi_steps, theta_steps)
            #   coor = np.array(np.meshgrid(phi_steps, theta_steps))

            grid = np.array(np.meshgrid(phi_steps, theta_steps))

            # for phi in phi_steps:
            for theta, phi in zip(grid[0].flatten(), grid[1].flatten()):
                # """ loop over theta """
                # for theta in theta_steps:
                volume, npoints = volume_fraction(dl, theta, phi, vxyz, paths, faces_cube)
                file.write(
                    str(theta) + ' ' + str(phi) + ' ' + str(volume) + \
                    ' ' + str(float(npoints)) + '\n')
            file.write('\n')


if __name__ == "__main__":
    _s = time.time()

    """ Cut Plane: Test Cases """
    # Case 1: three intersection points, parallel to xy-plane.
    case1 = {"dl": 0.4330127, "theta": 45.0, "phi": 54.73561032}

    # Case 2: four intersection points, not parallel to any axis.
    case2 = {"dl": 0.5, "theta": 1, "phi": 89}

    # Case 3: four intersection points, parallel to yz-plane.
    case3 = {"dl": 0.5, "theta": 0, "phi": 90}

    # Case 4: five intersection points
    case4 = {"dl": 1, "theta": 50, "phi": 20}

    # Case 5: six intersection points
    case5 = {"dl": 1, "theta": 50, "phi": 50}

    cases = [case1, case2, case3, case4, case5] * 1000

    """ Unit Cube """
    # coordinates of the unit cube. The vertices are numbered as indicated in cube_paths.pdf
    #                 V0   V1   V2   V3   V4   V5   V6   V7
    vxyz = np.array([[0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0],
                     [0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                     [0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0]], dtype=np.double).T

    # key: ipath, value: path described by the vertex numbers
    paths = {0: [0, 1, 4, 7],  # light gray + x + z + y path (V1 - V2 - V5 - V8) --> p0
             1: [1, 5],  # dotted light gray parallel to y-axis (V2 - V6)        --> p1
             2: [0, 2, 5, 7],  # gray + y + x + z path (V1 - V3 - V6 - V8)       --> p2
             3: [2, 6],  # dotted gray parallel to z-axis (V3 - V7)              --> p3
             4: [0, 3, 6, 7],  # black + z + y + x path (V1 - V4 - V7 - V8)      --> p4
             5: [3, 4]}  # dotted black parallel to x-axis (V4 - V5)             --> p5

    # faces that attached to each vertex
    faces_cube = {0: [0, 2, 5, 1],  # F0
                  1: [1, 4, 7, 5],  # F1
                  2: [4, 3, 6, 7],  # F2
                  3: [0, 3, 6, 2],  # F3
                  4: [0, 1, 4, 3],  # F4
                  5: [2, 5, 7, 6]}  # F5

    # volume = []
    # for case in cases:
    #     """ Test Case """
    #     dl = case["dl"]
    #     theta = case["theta"]  # degree
    #     phi = case["phi"]  # degree
    #
    #     volume.append(volume_fraction(dl, theta, phi, vxyz, paths, faces_cube))

    # print("volume = {}".format(volume[:5]))

    volume_plane_cube_intersection(ndl=32, ntheta=99, nphi=99, vxyz=vxyz, paths=paths, faces_cube=faces_cube)

    _end = time.time()
    print("Time elapsed: {} s".format(_end - _s))
