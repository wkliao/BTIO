#
#  Copyright (C) 2013, Northwestern University
#  See COPYRIGHT notice in top-level directory.
#
#
# Please change the following variables:
#    MPIF90        -- MPI Fortran compiler
#    FCFLAGS       -- Compile flag
#    PNETCDF_DIR   -- PnetCDF library installation directory
#

SUFFIXES = .o .f90

MPIF90       = mpif90
FCFLAGS      = -O2
PNETCDF_DIR  = $(HOME)/PnetCDF/1.10.0

COMPILE_F90  = $(MPIF90) $(FCFLAGS) $(INC) -c
LINK         = $(MPIF90) $(FCFLAGS)
INC          = -I$(PNETCDF_DIR)/include
LIBS         = -L$(PNETCDF_DIR)/lib -lpnetcdf

SRCS = io_info.f90 header.f90 mpiio_m.f90 make_set.f90 pnetcdf_m.f90 bt.f90

OBJS = $(SRCS:.f90=.o)

TARGET = btio

all: $(TARGET)

%.o:%.f90
	$(COMPILE_F90) $<

$(TARGET): $(OBJS)
	$(LINK) $(OBJS) -o $(TARGET) $(LIBS)

io_info.o:        io_info.f90
header.o:         header.f90
mpiio_m.o:        mpiio_m.f90 header.o
make_set.o:       make_set.f90 header.o
pnetcdf_m.o:      pnetcdf_m.f90 header.o
bt.o:             bt.f90 header.o mpiio_m.o pnetcdf_m.f90

PACKAGE_NAME = btio-pnetcdf-1.1.1

PACKING_LIST = $(SRCS) Makefile README.md COPYRIGHT inputbt.data RELEASE_NOTES

dist:
	/bin/rm -rf $(PACKAGE_NAME) $(PACKAGE_NAME).tar.gz
	mkdir $(PACKAGE_NAME)
	cp $(PACKING_LIST) $(PACKAGE_NAME)
	tar -cf $(PACKAGE_NAME).tar $(PACKAGE_NAME)
	gzip $(PACKAGE_NAME).tar
	/bin/rm -rf $(PACKAGE_NAME)

clean:
	- rm -f *.o core* *.mod $(TARGET)

distclean: clean
	/bin/rm -rf $(PACKAGE_NAME).tar.gz $(PACKAGE_NAME)

