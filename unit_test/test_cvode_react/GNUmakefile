PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gnu

USE_MPI    = FALSE
USE_OMP    = FALSE
USE_CUDA   = TRUE

USE_CUDA_CVODE = TRUE
USE_CVODE_CUSOLVER ?= FALSE

USE_REACT = TRUE
USE_REACT_SPARSE_JACOBIAN = TRUE
USE_SPARSE_STOP_ON_OOB = FALSE

EBASE = main

# define the location of the CASTRO top directory
MICROPHYSICS_HOME  := ../..

# This sets the EOS directory in Castro/EOS -- note: gamma_law will not work,
# you'll need to use gamma_law_general
EOS_DIR     := helmholtz

# This sets the network directory in Castro/Networks
NETWORK_DIR := aprox13

CONDUCTIVITY_DIR := stellar

ifeq ($(USE_CUDA), TRUE)
  INTEGRATOR_DIR := VODE90
endif

EXTERN_SEARCH += .

Bpack   := ./Make.package
Blocs   := .

include $(MICROPHYSICS_HOME)/Make.Microphysics


