#
#  Copyright (C) 2013, Northwestern University
#  See COPYRIGHT notice in top-level directory.
#
#  $Id: README 4580 2017-11-14 04:39:36Z wkliao $

This software benchmarks the performance of Parallel NetCDF and MPI-IO methods
for the I/O and data partitioning pattern from the NASA's NAS Parallel
Benchmarks (NPB) suite (http://www.nas.nasa.gov/publications/npb.html).

BTIO presents a block-tridiagonal partitioning pattern on a three-dimensional
array across a square number of processes.  Each process is responsible for
multiple Cartesian subsets of the entire data set, whose number increases with
the square root of the number of processes participating in the computation.
Multiple global arrays are consecutively written to a shared file by appending
one after another. The number of global arrays can be adjusted in the input
parameter file.  Each array is stored in the file in a canonical, row-major
order.  The global arrays are later read back, using the same I/O and data
partitioning pattern. The size of global array can also be adjusted from
the input parameter file, named "inputbt.data". For an illustration of data
partitioning pattern, please refer to:
    Wei-keng Liao. "Design and Evaluation of MPI File Domain Partitioning
    Methods under Extent-Based File Locking Protocol", in the IEEE Transactions
    on Parallel and Distributed Systems, 22(2):260-272, February 2011.

Instruction for compiling:
    Edit Makefile and change the following 3 variables
        MPIF90        -- MPI Fortran compiler
        FCFLAGS       -- Compile flag
        PNETCDF_DIR   -- the path of PnetCDF library
                         (1.4.0 and higher is required)

    For example:
        MPIF90      = mpif90
        FCFLAGS     = -O2
        PNETCDF_DIR = ${HOME}/PnetCDF

    Run command "make" to build the executable, named "btio".

Instruction for running:
    The input parameter file named "inputbt.data" is required to run the
    benchmark. Users can set the following parameters in the file.
        w        : IO mode: w for write, r for read
        IO method: for example, 0 for using MPI collective I/O
        number of time steps
        grid_points(1), grid_points(2), grid_points(3)
        directory name for storing the input and output files

    For example, the contents of "inputbt.data" file are:
        w                  # IO mode: w for write, r for read
        3                  # IO method: 0 for MPI collective IO, 1 for MPI independent IO, 2 for PnetCDF blocking I/O, 3 for PnetCDF nonblocking I/O
        40                 # number of writes/reads
        512 512 512        # grid_points(1), grid_points(2), grid_points(3)
        /scratch2/scratchdirs/wkliao/FS_1M_128
    which set 1) w to perform write operations only; 2) the I/O method to
    using Parallel netCDF nonblocking APIs; 3) number of global arrays for
    write and read; 4) the global array size; 5) the input/output directory
    name.

    Note that btio creates a file named btio.nc containing a single 5D array of
    size NUM_DUMPS x Z x Y x X x FIVE_DBL. The data type of variable is double.
    NUM_DUMPS corresponds to the number of writes (or reads for read case). Z,
    Y, and X correspond to grid_points(3), grid_points(2), grid_points(1),
    respectively. FIVE_DBL is the fifth dimension of size 5.

    Command to run the MPI job:
        mpiexec -n 1024 btio
    or
        mpiexec -n 1024 btio inputbt.data
    The only optional command-line argument is a file name. In this example, it
    is "inputbt.data". This argument allows to use a different input file name
    from the default "inputbt.data".

Example output from the standard out:
    -- BT-IO Benchmark (write operation only) --
    Number of MPI processes  :      1024
    Global array size X      :       512
    Global array size Y      :       512
    Global array size Z      :       512
    Number of I/O iterations :        40
    Total I/O amount         :    204800.00 MiB
    Time in sec              :        41.95
    I/O bandwidth                   4882.00 MiB/s
    ------------------------------------------
    Using Parallel netCDF non-blocking I/O method
    output file path         : /scratch2/scratchdirs/wkliao/FS_1M_128
    file striping count      :       128
    file striping size       :   1048576 bytes

Questions/Comments:
    email: wkliao@eecs.northwestern.edu

