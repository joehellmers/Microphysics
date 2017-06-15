module react_zones_module
#ifdef CUDA
  use cudafor
#endif
  use bl_types
  use bl_space
  use bl_constants_module, only: ZERO
  use network, only: nspec
  use burn_type_module, only: burn_t, normalize_abundances_burn
  use actual_burner_module, only : actual_burner
  use managed_probin_module, only: cu_tmax
  
  implicit none

  type :: pfidx_t
     integer :: itemp
     integer :: irho
     integer :: ispec
     integer :: ispec_old
     integer :: irodot
     integer :: irho_hnuc
     integer :: endrho
     integer :: endT
     integer :: endX
     integer :: ncomps
  end type pfidx_t
  
contains

#ifdef CUDA  
  attributes(global) &
#endif
  subroutine react_zones_flat(state, pfidx, sdim)
    implicit none
  
    integer, value  :: sdim
    type(pfidx_t)   :: pfidx
    real(kind=dp_t) &
#ifdef CUDA
         , device &
#endif
         :: state(:, :)
    type (burn_t)   :: burn_state_in, burn_state_out
    integer         :: ii, j    

#ifdef CUDA    
    ii = (blockIdx%x - 1) * blockDim % x + threadIdx % x
    if (ii <= sdim) then
#else
    !$OMP PARALLEL DO PRIVATE(ii,j) &
    !$OMP PRIVATE(burn_state_in, burn_state_out) &
    !$OMP SCHEDULE(DYNAMIC,1)
    do ii = 1, sdim
#endif
       burn_state_in % rho = state(ii, pfidx % irho)
       burn_state_in % T = state(ii, pfidx % itemp)
       do j = 1, nspec
          burn_state_in % xn(j) = state(ii, pfidx % ispec_old + j - 1)
       enddo

       call normalize_abundances_burn(burn_state_in)

       ! Integrator doesn't care about the initial internal energy.
       burn_state_in % e = ZERO

       call actual_burner(burn_state_in, burn_state_out, cu_tmax, ZERO)

       do j = 1, nspec
          state(ii, pfidx % ispec + j - 1) = burn_state_out % xn(j)
       enddo

       do j=1, nspec
          ! Explicit loop needed here to keep the GPU happy if running on GPU
          state(ii, pfidx % irodot + j - 1) = &
               (burn_state_out % xn(j) - burn_state_in % xn(j)) / cu_tmax
       enddo

       state(ii, pfidx % irho_hnuc) = state(ii, pfidx % irho) * &
            (burn_state_out % e - burn_state_in % e) / cu_tmax
#ifdef CUDA       
    end if
#else
    enddo
    !$OMP END PARALLEL DO
#endif
  end subroutine react_zones_flat

#ifndef CUDA  
  subroutine react_zones(state, lo, hi, pfidx)
    implicit none
  
    type(pfidx_t)   :: pfidx
    real(kind=dp_t) :: state(:, :, :, :)
    type (burn_t)   :: burn_state_in, burn_state_out
    integer         :: ii, jj, kk, j, lo(3), hi(3)

    !$OMP PARALLEL DO PRIVATE(ii,jj, kk, j) &
    !$OMP PRIVATE(burn_state_in, burn_state_out) &
    !$OMP SCHEDULE(DYNAMIC,1)
    do kk = lo(3), hi(3)
       do jj = lo(2), hi(2)
          do ii = lo(1), hi(1)
             burn_state_in % rho = state(ii, jj, kk, pfidx % irho)
             burn_state_in % T = state(ii, jj, kk, pfidx % itemp)
             do j = 1, nspec
                burn_state_in % xn(j) = state(ii,jj,kk,pfidx % ispec_old+j-1)
             enddo

             call normalize_abundances_burn(burn_state_in)

             ! Integrator doesn't care about the initial internal energy.
             burn_state_in % e = ZERO

             call actual_burner(burn_state_in, burn_state_out, cu_tmax, ZERO)

             do j = 1, nspec
                state(ii, jj, kk, pfidx % ispec + j - 1) = burn_state_out % xn(j)
             enddo

             do j=1, nspec
                state(ii, jj, kk, pfidx % irodot + j - 1) = &
                     (burn_state_out % xn(j) - burn_state_in % xn(j)) / cu_tmax
             enddo

             state(ii, jj, kk, pfidx % irho_hnuc) = state(ii, jj, kk, pfidx % irho) * &
                  (burn_state_out % e - burn_state_in % e) / cu_tmax
          end do
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine react_zones
#endif
  
end module react_zones_module
