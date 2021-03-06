module actual_rhs_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module
  use physical_constants, only: N_AVO
  use network
  use table_rates
  use burn_type_module

  implicit none

  type :: rate_eval_t
     real(rt) :: unscreened_rates(4, nrates)
     real(rt) :: screened_rates(nrates)
     real(rt) :: dqweak(nrat_tabular)
     real(rt) :: epart(nrat_tabular)
  end type rate_eval_t
  
contains

  subroutine actual_rhs_init()
    ! STUB FOR MAESTRO'S TEST_REACT. ALL THE INIT IS DONE BY BURNER_INIT
    return
  end subroutine actual_rhs_init


  subroutine update_unevolved_species(state)
    ! STUB FOR INTEGRATOR
    type(burn_t)     :: state

    !$gpu
    
    return
  end subroutine update_unevolved_species


  subroutine zero_rate_eval(rate_eval)

    implicit none

    type(rate_eval_t), intent(inout) :: rate_eval

    !$gpu

    rate_eval % unscreened_rates(i_rate, :) = ZERO
    rate_eval % unscreened_rates(i_drate_dt, :) = ZERO
    rate_eval % unscreened_rates(i_scor, :) = ONE
    rate_eval % unscreened_rates(i_dscor_dt, :) = ZERO
    rate_eval % screened_rates = ZERO
    rate_eval % dqweak = ZERO
    rate_eval % epart = ZERO

  end subroutine zero_rate_eval


  subroutine evaluate_rates(state, rate_eval)
    !$acc routine seq

    use reaclib_rates, only: screen_reaclib, reaclib_evaluate
    use screening_module, only: screen5, plasma_state, fill_plasma_state

    implicit none
    
    type(burn_t)     :: state
    type(rate_eval_t), intent(out) :: rate_eval
    type(plasma_state) :: pstate
    real(rt) :: Y(nspec)
    real(rt) :: raw_rates(4, nrates)
    real(rt) :: reactvec(num_rate_groups+2)
    integer :: i, j
    real(rt) :: rhoy, scor, dscor_dt, dscor_dd

    !$gpu

    Y(:) = state % xn(:) * aion_inv(:)
    rhoy = state % rho * state % y_e

    ! Zero out the rates
    call zero_rate_eval(rate_eval)

    ! Calculate Reaclib rates
    call fill_plasma_state(pstate, state % T, state % rho, Y)
    do i = 1, nrat_reaclib
       call reaclib_evaluate(pstate, state % T, i, reactvec)
       rate_eval % unscreened_rates(:,i) = reactvec(1:4)
    end do

    ! Evaluate screening factors
    if (screen_reaclib) then

      call screen5(pstate, 1, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1) = dscor_dt

    end if


    ! Compute screened rates
    rate_eval % screened_rates = rate_eval % unscreened_rates(i_rate, :) * &
                                 rate_eval % unscreened_rates(i_scor, :)

  end subroutine evaluate_rates


  subroutine actual_rhs(state)
    
    !$acc routine seq

    use extern_probin_module, only: do_constant_volume_burn
    use burn_type_module, only: net_itemp, net_ienuc
    use sneut_module, only: sneut5
    use temperature_integration_module, only: temperature_rhs

    implicit none

    type(burn_t) :: state
    type(rate_eval_t) :: rate_eval
    real(rt) :: Y(nspec), ydot_nuc(nspec)
    real(rt) :: reactvec(num_rate_groups+2)
    integer :: i, j
    real(rt) :: rhoy, ye, enuc
    real(rt) :: sneut, dsneutdt, dsneutdd, snuda, snudz

    !$gpu

    ! Set molar abundances
    Y(:) = state % xn(:) * aion_inv(:)

    call evaluate_rates(state, rate_eval)

    call rhs_nuc(state, ydot_nuc, Y, rate_eval % screened_rates)
    state % ydot(1:nspec) = ydot_nuc

    ! ion binding energy contributions
    call ener_gener_rate(ydot_nuc, enuc)

    ! weak Q-value modification dqweak (density and temperature dependent)
    
    ! weak particle energy generation rates from gamma heating and neutrino loss
    ! (does not include plasma neutrino losses)


    ! Get the neutrino losses
    call sneut5(state % T, state % rho, state % abar, state % zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)

    ! Append the energy equation (this is erg/g/s)
    state % ydot(net_ienuc) = enuc - sneut

    ! Append the temperature equation
    call temperature_rhs(state)

  end subroutine actual_rhs


  subroutine rhs_nuc(state, ydot_nuc, Y, screened_rates)

    !$acc routine seq

    implicit none

    type (burn_t), intent(in) :: state
    real(rt), intent(out) :: ydot_nuc(nspec)
    real(rt), intent(in)  :: Y(nspec)
    real(rt), intent(in)  :: screened_rates(nrates)

    !$gpu



    ydot_nuc(jhe4) = ( &
      -screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4)*state % rho &
       )

    ydot_nuc(jc12) = ( &
      -screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4)*state % rho &
       )

    ydot_nuc(jo16) = ( &
      screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4)*state % rho &
       )


  end subroutine rhs_nuc


  subroutine actual_jac(state)

    !$acc routine seq

    use burn_type_module, only: net_itemp, net_ienuc
    use sneut_module, only: sneut5
    use temperature_integration_module, only: temperature_jac
    use jacobian_sparsity_module, only: get_jac_entry, set_jac_entry, set_jac_zero

    implicit none
    
    type(burn_t) :: state
    type(rate_eval_t) :: rate_eval
    real(rt) :: reactvec(num_rate_groups+2)
    real(rt) :: screened_rates_dt(nrates)
    real(rt) :: Y(nspec), yderivs(nspec)
    real(rt) :: ye, rhoy, b1, scratch
    real(rt) :: sneut, dsneutdt, dsneutdd, snuda, snudz
    integer  :: j, k

    !$gpu

    ! Set molar abundances
    Y(:) = state % xn(:) * aion_inv(:)
    
    call evaluate_rates(state, rate_eval)

    ! Zero out the Jacobian
    call set_jac_zero(state)

    ! Species Jacobian elements with respect to other species
    call jac_nuc(state, Y, rate_eval % screened_rates)

    ! Evaluate the species Jacobian elements with respect to temperature by
    ! calling the RHS using the temperature derivative of the screened rate
    screened_rates_dt = rate_eval % unscreened_rates(i_rate, :) * &
                        rate_eval % unscreened_rates(i_dscor_dt, :) + &
                        rate_eval % unscreened_rates(i_drate_dt, :) * &
                        rate_eval % unscreened_rates(i_scor, :)

    call rhs_nuc(state, yderivs, Y, screened_rates_dt)

    do k = 1, nspec
       call set_jac_entry(state, k, net_itemp, yderivs(k))
    enddo

    ! Energy generation rate Jacobian elements with respect to species
    do j = 1, nspec
       do k = 1, nspec
          call get_jac_entry(state, k, j, yderivs(k))
       enddo
       call ener_gener_rate(yderivs, scratch)
       call set_jac_entry(state, net_ienuc, j, scratch)
    enddo

    ! Account for the thermal neutrino losses
    call sneut5(state % T, state % rho, state % abar, state % zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)

    do j = 1, nspec
       b1 = ((aion(j) - state % abar) * state % abar * snuda + (zion(j) - state % zbar) * state % abar * snudz)
       call get_jac_entry(state, net_ienuc, j, scratch)
       scratch = scratch - b1
       call set_jac_entry(state, net_ienuc, j, scratch)
    enddo

    ! Energy generation rate Jacobian element with respect to temperature
    do k = 1, nspec
       call get_jac_entry(state, k, net_itemp, yderivs(k))
    enddo
    call ener_gener_rate(yderivs, scratch)
    scratch = scratch - dsneutdt    
    call set_jac_entry(state, net_ienuc, net_itemp, scratch)

    ! Temperature Jacobian elements
    call temperature_jac(state)

  end subroutine actual_jac


  subroutine jac_nuc(state, Y, screened_rates)

    !$acc routine seq

    use jacobian_sparsity_module, only: set_jac_entry

    implicit none

    type(burn_t), intent(inout) :: state
    real(rt), intent(in)  :: Y(nspec)
    real(rt), intent(in)  :: screened_rates(nrates)
    real(rt) :: scratch


    !$gpu


    scratch = (&
      -screened_rates(k_he4_c12__o16)*Y(jc12)*state % rho &
       )
    call set_jac_entry(state, jhe4, jhe4, scratch)

    scratch = (&
      -screened_rates(k_he4_c12__o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jc12, scratch)

    scratch = (&
      -screened_rates(k_he4_c12__o16)*Y(jc12)*state % rho &
       )
    call set_jac_entry(state, jc12, jhe4, scratch)

    scratch = (&
      -screened_rates(k_he4_c12__o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jc12, jc12, scratch)

    scratch = (&
      screened_rates(k_he4_c12__o16)*Y(jc12)*state % rho &
       )
    call set_jac_entry(state, jo16, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_c12__o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jo16, jc12, scratch)


  end subroutine jac_nuc

end module actual_rhs_module
