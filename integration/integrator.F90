module integrator_module

  implicit none

  public
  
contains

  subroutine integrator_init()

    use actual_integrator_module, only: actual_integrator_init

    implicit none

    call actual_integrator_init()

  end subroutine integrator_init

  subroutine integrator(state_in, state_out, dt, time)

    !$acc routine seq

    use vode_parameters_module, only: grid_size
    use actual_integrator_module, only: actual_integrator
    use burn_type_module, only: burn_t
    use bl_types, only: dp_t

    implicit none

    type (burn_t),  intent(in   ) :: state_in(grid_size)
    type (burn_t),  intent(inout) :: state_out(grid_size)
    real(dp_t),     intent(in   ) :: dt, time

    call actual_integrator(state_in, state_out, dt, time)

  end subroutine integrator

end module integrator_module
