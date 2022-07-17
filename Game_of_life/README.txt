A modified version of the Game of Life cellular automation rule by J. H. Conway.

The program is parallel and prepared to run using the MPI library. 
It partitions a 2D grid in a toroidal topology and also does
a parallel I/O.

In "tests" the input files and the corresponding outputs can be seen.
Run the program with n_rows*n_cols MPI procceses.