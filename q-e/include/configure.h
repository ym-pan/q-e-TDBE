!
! Copyright (C) 2006 WanT Group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!
! contains configure infos
!

#define  __HAVE_CONFIG_INFO

#define   __CONF_HOST           "@host@"
#define   __CONF_ARCH           "x86_64"

#define   __CONF_CC             "/shared/apps/.intel/2022/intelmpi/2022.1/mpi/latest/bin/mpiicc"
#define   __CONF_CFLAGS         "-O3"
#define   __CONF_DFLAGS         " -D__FFTW -D__MPI -D__SCALAPACK"
#define   __CONF_CPP            "cpp"
#define   __CONF_CPPFLAGS       "-P -traditional -Uvector"

#define   __CONF_F90            "/shared/apps/intel/2022.1/compiler/latest/linux/bin/ifx"
#define   __CONF_MPIF90         "mpiifort"
#define   __CONF_F90FLAGS       "$(FFLAGS) -nomodule"
#define   __CONF_F77            "@f77@"
#define   __CONF_FFLAGS         "-O2 -assume byterecl -g -traceback"
#define   __CONF_FFLAGS_NOOPT   "-O0 -assume byterecl -g -traceback"
#define   __CONF_PRE_FDFLAGS    "-fpp -allow nofpp_comments "
#define   __CONF_FDFLAGS        "$(DFLAGS) $(MANUAL_DFLAGS)"

#define   __CONF_LD             "mpiifort"
#define   __CONF_LDFLAGS        ""
#define   __CONF_IMOD           "-I"

#define   __CONF_BLAS_LIBS      " -lblas "
#define   __CONF_LAPACK_LIBS    "-L/opt/slurm/current/lib/slurm -llapack  -lblas "
#define   __CONF_FFT_LIBS       ""
#define   __CONF_MASS_LIBS      ""

#define   __CONF_AR             "ar"
#define   __CONF_ARFLAGS        "ruv"
#define   __CONF_RANLIB         "ranlib"


