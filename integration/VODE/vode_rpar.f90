! rpar are the real quantities that are passed through the VODE call to the
! RHS and Jacobian routines

module rpar_indices

  implicit none

  integer :: n_rpar_comps = 0

  integer :: irp_dens, irp_cv, irp_cp, irp_dedY, irp_dhdY
  integer :: irp_abar, irp_zbar
  integer :: irp_eta, irp_ye
  integer :: irp_self_heat
  integer :: irp_rates

contains

  function get_next_rpar_index(num) result (next)

    implicit none
    
    ! Return the next starting index for a plotfile quantity,
    ! and increment the counter of plotfile quantities by num.
    
    integer, intent(in) :: num
    integer             :: next

   ! Return the value of the first index of this group of data.

    next = n_rpar_comps + 1

    ! Update our local record of how many variables we are using.
    
    n_rpar_comps = n_rpar_comps + num

  end function get_next_rpar_index


  subroutine init_rpar_indices(nrates, nspec)

    use burn_type_module, only: num_rate_groups

    implicit none
    
    integer, intent(in) :: nrates, nspec

    irp_dens      = get_next_rpar_index(1)
    irp_cp        = get_next_rpar_index(1)
    irp_cv        = get_next_rpar_index(1)
    irp_abar      = get_next_rpar_index(1)
    irp_zbar      = get_next_rpar_index(1)
    irp_eta       = get_next_rpar_index(1)
    irp_ye        = get_next_rpar_index(1)
    irp_self_heat = get_next_rpar_index(1)
    irp_dhdY      = get_next_rpar_index(nspec)
    irp_dedY      = get_next_rpar_index(nspec)
    irp_rates     = get_next_rpar_index(num_rate_groups * nrates)

  end subroutine init_rpar_indices

end module rpar_indices