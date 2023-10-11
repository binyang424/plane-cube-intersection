./polyhedra_volume.exe polyhedra_data/tetrahedron_01.dat
./polyhedra_volume.exe polyhedra_data/tetrahedron_02.dat
./polyhedra_volume.exe polyhedra_data/tetrahedron_03.dat
./polyhedra_volume.exe polyhedra_data/tetrahedron_04.dat
./polyhedra_volume.exe polyhedra_data/tetrahedron_05.dat
./polyhedra_volume.exe polyhedra_data/tetrahedron_06.dat
./polyhedra_volume.exe polyhedra_data/tetrahedron_07.dat
./polyhedra_volume.exe polyhedra_data/tetrahedron_08.dat
./polyhedra_volume.exe polyhedra_data/tetrahedron_09.dat
./polyhedra_volume.exe polyhedra_data/tetrahedron_10.dat
./polyhedra_volume.exe polyhedra_data/cube_01.dat
./polyhedra_volume.exe polyhedra_data/cube_02.dat
./polyhedra_volume.exe polyhedra_data/cube_03.dat
./polyhedra_volume.exe polyhedra_data/cube_04.dat
./polyhedra_volume.exe polyhedra_data/cube_05.dat
./polyhedra_volume.exe polyhedra_data/octahedron_01.dat
./polyhedra_volume.exe polyhedra_data/octahedron_02.dat
./polyhedra_volume.exe polyhedra_data/octahedron_03.dat
./polyhedra_volume.exe polyhedra_data/octahedron_04.dat
./polyhedra_volume.exe polyhedra_data/octahedron_05.dat
./polyhedra_volume.exe polyhedra_data/dodecahedron_01.dat
./polyhedra_volume.exe polyhedra_data/dodecahedron_02.dat
./polyhedra_volume.exe polyhedra_data/dodecahedron_03.dat
./polyhedra_volume.exe polyhedra_data/dodecahedron_04.dat
./polyhedra_volume.exe polyhedra_data/dodecahedron_05.dat
./polyhedra_volume.exe polyhedra_data/dodecahedron_06.dat
./polyhedra_volume.exe polyhedra_data/icosahedron_01.dat
./polyhedra_volume.exe polyhedra_data/triacontahedron_01.dat
echo ' '
echo ' '
echo 'some mathematica polyhedra:'
./polyhedra_volume.exe polyhedra_data/cube_mathematica_01.dat
./polyhedra_volume.exe polyhedra_data/snub_cube_mathematica_01.dat
echo ' ' 
echo ' ' 
echo 'non-convex test cases:'
./polyhedra_volume.exe polyhedra_data/120hedron_01.dat
./polyhedra_volume.exe polyhedra_data/small_stellated_dodecahedron_mathematica.dat