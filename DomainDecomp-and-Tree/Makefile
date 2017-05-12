#-----------------------------------------------------------------
#
# You might be looking for the compile-time Makefile options of the code...
#
# They have moved to a separate file.
#
# To build the code, do the following:
#
#  (1) Copy the file "Template-Config.sh"  to  "Config.sh"
#
#        cp Template-Config.sh Config.sh 
#
#  (2) Edit "Config.sh" as needed for your application
#
#  (3) Run "make"
#
#
#  New compile-time options should be added to the 
#  file "Template-Config.sh" only. Usually, the should be added
#  there in the disabled/default version.
#
#  "Config.sh" should *not* be checked in to the repository
#
#  Note: It is possible to override the default name of the 
#  Config.sh file, if desired, as well as the name of the
#  executable. For example:
#
#   make  CONFIG=MyNewConf.sh  EXEC=Gadget_new
# 
#-----------------------------------------------------------------
#
# You might also be looking for the target system SYSTYPE option
#
# It has also moved to a separate file.
#
# To build the code, do the following:
#
# (A) set the SYSTYPE variable in your .bashrc (or similar file):
#
#        e.g. export SYSTYPE=Magny
# or
#
# (B) set SYSTYPE in Makefile.systype 
#     This file has priority over your shell variable.:
#
#    (1) Copy the file "Template-Makefile.systype"  to  "Makefile.systype"
#
#        cp Template-Makefile.systype Makefile.systype 
#
#    (2) Uncomment your system in  "Makefile.systype".
#
# If you add an ifeq for a new system below, also add that systype to
# Template-Makefile.systype

ifdef SYSTYPE
SYSTYPE := "$(SYSTYPE)"
-include Makefile.systype
else
include Makefile.systype
endif

ifeq ($(wildcard Makefile.systype), Makefile.systype)
INCL = Makefile.systype
else
INCL =
endif


CONFIG   =  Config.sh
PERL     =  /usr/bin/perl

RESULT     := $(shell CONFIG=$(CONFIG) PERL=$(PERL) make -f config-makefile)
CONFIGVARS := $(shell cat gadgetconfig.h)


CC       = mpicc        # sets the C-compiler (default)
CXX       = mpiCC       # sets the C++-compiler (default)

FC 	 = mpif90

OPTIMIZE = -Wall  -g   # optimization and warning flags (default)

MPICHLIB = -lmpich





#----------------------------------------------------------------------------------------------
feq ($(SYSTYPE),"pico")
CC       =   mpicc # sets the gcc C-compiler
#OPTIMIZE =  -O2
OPTIMIZE =  -O0 -g  #-fp-stack-check -check bounds -fpe0
GSL_INCL = -I/gss/gss_work/DRES_murante/local/include
GSL_LIBS = -L/gss/gss_work/DRES_murante/local/lib
FFTW_INCL= -I/gss/gss_work/DRES_murante/local/include
FFTW_LIBS= -L/gss/gss_work/DRES_murante/local/lib
HDF5INCL = -I/gss/gss_work/DRES_murante/local/include
HDF5LIB  = -L/gss/gss_work/DRES_murante/local/lib -lhdf5
IMF_INCL = -I./ImfLib
IMF_LIBS = -L./ImfLib -lIMF
MPICHLIB =
endif

ifeq ($(SYSTYPE),"galileo")
CC       =   mpicc # sets the gcc C-compiler
OPTIMIZE =  -O0 -g
#OPTIMIZE =  -O2
GSL_INCL = -I$(GSL_INCLUDE)
GSL_LIBS = -L$(GSL_LIB)
FFTW_INCL= -I/galileo/home/userexternal/gmurante/local/include
FFTW_LIBS= -L/galileo/home/userexternal/gmurante/local/lib
HDF5INCL = -I$(HDF5_INCLUDE)
HDF5LIB  = -L$(HDF5_LIB) -lhdf5
#HDF5INCL = -I/cineca/prod/libraries/hdf5/1.8.9_par/openmpi--1.4.4-dbg--intel--co-2011.6.233--binary/include
#HDF5LIB  = -L/cineca/prod/libraries/hdf5/1.8.9_par/openmpi--1.4.4-dbg--intel--co-2011.6.233--binary/lib -lhdf5
IMF_INCL = -I./ImfLib
IMF_LIBS = -L./ImfLib -lIMF
MPICHLIB =
#OPT      +=  -DNOCALLSOFSYSTEM
#OPT      +=  -DAUTO_SWAP_ENDIAN_READIC
endif




ifneq (HAVE_HDF5,$(findstring HAVE_HDF5,$(CONFIGVARS)))
ifneq (LT_METAL_COOLING_WAL,$(findstring LT_METAL_COOLING_WAL,$(CONFIGVARS)))
HDF5INCL =
HDF5LIB  =
endif
endif


OPTIONS = $(OPTIMIZE) $(OPT) 

FOPTIONS = $(OPTIMIZE) $(FOPT)

EXEC   = DDT

OBJS  =  domain.o main.o \
	 begrun.o run.o allvars.o io.o\
	 init.o forcetree.o \
	 forcetree_update.o read_ic.o  ngb.o  \
	 system.o allocate.o endrun.o parallel_sort.o\
	 peano.o \
	 mpi_util.o mymalloc.o  \
         compile_time_info.o look_around.o

INCL    += allvars.h proto.h forcetree.h domain.h  \
	 myqsort.h Makefile utilities.h\


CFLAGS = $(OPTIONS) $(GSL_INCL)  


FFTW = $(FFTW_LIBS)  $(FFTW_LIBNAMES) 

LIBS   = -lm $(HDF5LIB) -g $(MPICHLIB) $(GSL_LIBS) -lgsl -lgslcblas $(FFTW)


$(EXEC): $(OBJS) $(FOBJS)  
	$(FC) $(OPTIMIZE) $(OBJS) $(FOBJS) $(LIBS) $(RLIBS) -o $(EXEC)

$(OBJS): $(INCL)  $(CONFIG)  compile_time_info.c


$(FOBJS): $(FINCL)

complile_time_info.c: $(CONFIG)
	$(PERL) prepare-config.perl $(CONFIG)

clean:
	rm -f $(OBJS) $(FOBJS) $(EXEC) gadgetconfig.h


