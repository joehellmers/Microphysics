module integrator_module

  implicit none

  public

contains

  subroutine integrator_init()

#if (INTEGRATOR == 0 || INTEGRATOR == 1)
    use vode_integrator_module, only: vode_integrator_init
    use bs_integrator_module, only: bs_integrator_init
#else
    use actual_integrator_module, only: actual_integrator_init
#endif

    implicit none

#if (INTEGRATOR == 0 || INTEGRATOR == 1)
    call vode_integrator_init()
    call bs_integrator_init()
#else
    call actual_integrator_init()
#endif

  end subroutine integrator_init



  subroutine integrator(state_in, state_out, dt, time)

    !$acc routine seq

#if (INTEGRATOR == 0 || INTEGRATOR == 1)
    use vode_integrator_module, only: vode_integrator
    use bs_integrator_module, only: bs_integrator
#else
    use actual_integrator_module, only: actual_integrator
#endif
    use bl_error_module, only: bl_error
    use bl_constants_module, only: ZERO, ONE
    use burn_type_module, only: burn_t
    use bl_types, only: dp_t
    use integration_data, only: integration_status_t
    use extern_probin_module, only: rtol_spec, rtol_temp, rtol_enuc, &
                                    atol_spec, atol_temp, atol_enuc, &
                                    retry_burn, retry_burn_factor, retry_burn_max_change, &
                                    burning_mode, burning_mode_factor

    implicit none

    type (burn_t),  intent(in   ) :: state_in
    type (burn_t),  intent(inout) :: state_out
    real(dp_t),     intent(in   ) :: dt, time

    real(dp_t) :: limit_factor, t_sound, t_enuc

#if (INTEGRATOR == 0 || INTEGRATOR == 1)
    type (integration_status_t) :: status
    real(dp_t) :: retry_change_factor

    retry_change_factor = ONE

    status % atol_spec = atol_spec
    status % rtol_spec = rtol_spec

    status % atol_temp = atol_temp
    status % rtol_temp = rtol_temp

    status % atol_enuc = atol_enuc
    status % rtol_enuc = rtol_enuc

    do while (retry_change_factor <= retry_burn_max_change)

       status % integration_complete = .true.

#if (INTEGRATOR == 0)
       call vode_integrator(state_in, state_out, dt, time, status)
#elif (INTEGRATOR == 1)
       call bs_integrator(state_in, state_out, dt, time, status)
#endif

       if (status % integration_complete) exit

       ! If we got here, the integration failed, try the next burner.

       status % integration_complete = .true.

#if (INTEGRATOR == 0)
       print *, "Retrying burn with BS integrator"
       call bs_integrator(state_in, state_out, dt, time, status)
#elif (INTEGRATOR == 1)
       print *, "Retrying burn with VODE integrator"
       call vode_integrator(state_in, state_out, dt, time, status)
#endif

       if (status % integration_complete) exit

       ! If we got here, all integrations failed, loosen the tolerances.

       if (.not. retry_burn) then

          call bl_error("ERROR in burner: integration failed")

       else

          print *, "Retrying burn with looser tolerances"

          retry_change_factor = retry_change_factor * retry_burn_factor

          status % atol_spec = status % atol_spec * retry_burn_factor
          status % rtol_spec = status % rtol_spec * retry_burn_factor

          status % atol_temp = status % atol_temp * retry_burn_factor
          status % rtol_temp = status % rtol_temp * retry_burn_factor

          status % atol_enuc = status % atol_enuc * retry_burn_factor
          status % rtol_enuc = status % rtol_enuc * retry_burn_factor

       end if

    end do

    if (.not. status % integration_complete) then

       call bl_error("ERROR in burner: integration failed")

    endif

#else

    call actual_integrator(state_in, state_out, dt, time)

#endif

    ! For burning_mode == 3, limit the change in the state.

    if (burning_mode == 3) then

       t_enuc = min(state_in % e, state_out % e) / max(abs(state_out % e - state_in % e) / dt, 1.d-50)
       t_sound = max(state_out % dx / state_out % cs, 1.d-50)

       limit_factor = min(1.0d0, burning_mode_factor * t_enuc / t_sound)

       state_out % xn = state_in % xn + (state_out % xn - state_in % xn) * limit_factor
       state_out % e = state_in % e + (state_out % e - state_in % e) * limit_factor

    endif

  end subroutine integrator

end module integrator_module
