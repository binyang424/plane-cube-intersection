import numpy as np
import pypolylib
import os


path = "./polyhedra_data/"

filenames = os.listdir(path)
filenames = [filename for filename in filenames if filename.endswith(".dat")]

for filename in filenames:
    print(filename)

    with open(os.path.join(path, filename), "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if i == 1:
            volume_ref = float(line.strip())
            print("volume ref", float(line.strip().split(" ")[0]))

        if i == 2:
            number_vert = int(line.strip())
            print("number vertices", number_vert)

    points = np.zeros((number_vert, 3))

    for i, line in enumerate(lines):
        if 2 < i < 3 + number_vert:
            points[i - 3, :] = [float(x) for x in line.strip().split(" ")]
            # print("points", line.strip().split(" "))

        if i == 3 + number_vert:
            number_faces = int(line.strip())
            print("number faces", number_faces)

    faces = np.zeros((number_faces, 3))

    for i, line in enumerate(lines):
        if i % 2 == 1 and 3 + number_vert < i < 3 + number_vert + 2 * number_faces + 1:
            faces[int((i - 4 - number_vert) / 2), :] = \
                [int(x) for x in line.strip().split(" ")]
            # print("faces", line.strip().split(" "))

    # Python indices start at 0, while the polyhedra data starts at 1.
    # Adapt to Python flavor.
    faces = faces.astype(int) - 1

    """
    The volume should be independent of the common "height" used for each polygon choosing 
    any point on the surface as the height, say the first point, 
    
                            px = x(1)  ; py = y(1) ; pz = z(1)
    
    works fine for convex shapes, but fails for non-convex shapes such as stellated polyhedra. 
    
    Choosing an interior point as the height, say the "centroid", works fine for convex and 
    non-convex polyhedra. If only vertices have weight, then the mean value of the weighted 
    sum of all the points is the "centroid: 
    
                         r_c = sum (r_i * m_i) / nvert. 
    
    For unity mass weightings m, the even simpler:
    
                             r_c = sum (r_i) / nvert. 
    
    We'll use this interior point as the height. This effectly selects the chosen origin point, 
    an interior point, for symmetric polyhedra.
    """
    centroid = np.mean(points, axis=0)
    volume = 0

    for face in faces:
        # Get the three vertices of the face.
        vertices = np.vstack((points[face, :], centroid))
        volume += pypolylib.pyramid_volume(vertices)

    print("volume", volume)

    ratio = abs(volume / volume_ref - 1) * 100  # in percent
    if ratio < 1:
        print("volume correct")
        print("ratio", ratio)
    else:
        print("volume incorrect")
        print("ratio", ratio)
    print("-"*50)

print("Done!")