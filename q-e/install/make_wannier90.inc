#=======================================
# WANNIER90 	
#=======================================

TOPDIR = /home/yp94rucu/codes/q-e

F90=mpiifort
FCOPTS=$(FFLAGS) -nomodule 
LDOPTS=

LIBS = -L/opt/slurm/current/lib/slurm -llapack  -lblas   -lblas  
