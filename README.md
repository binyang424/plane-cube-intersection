This Python implementation provides a method to determine plane-cube intersections, originally implemented in Fortran by Dr. [Francis Timmes](https://search.asu.edu/profile/1096683) at Arizona State University. It allows you to calculate the volume of intersection between a plane and a cube using the following specified notations: 

<img src="[drawing.jpg](https://github.com/binyang424/plane-cube-intersection/blob/db0c73025d4e04192d9e3891023a52efcd86455d/111.jpeg)" alt="notations" style="width:400px;"/>

- `pyvolume_fraction.py` is the main function that can be used as an template and a tutorial of this code for volume fraction calculation.
- `pypolylib.py` contains functions for polygon area and polyhedron volume calculations.
- `pypolyhedra_volume.py` script provides examples of how to call the functions in `pypolylib.py` for calculating polygon areas and polyhedron volumes.
- `skspatial_verify.py` employs the `Scikit Spatial` package to calculate the intersection point of a line with a plane. This is done to validate the Python implementation.

This Python implementation is provided as-is and for reference. A copy of the original Fortran version is attached here (`./Fortran`) and a `README` file can be found in this folder (`./Fortran/README.md`), which contains the links point to the website  https://cococubed.com/ maintained by Dr. [Francis Timmes](https://search.asu.edu/profile/1096683).
