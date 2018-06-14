ifdef SDC
  F90sources += integrator_sdc.F90
ifndef CUDA
  f90sources += numerical_jacobian_sdc.f90
endif
else
  F90sources += integrator.F90
ifndef CUDA
  F90sources += numerical_jacobian.F90
endif
endif
f90sources += integration_data.f90
F90sources += temperature_integration.F90
F90sources += rpar.F90
