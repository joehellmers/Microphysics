! Setup a 3D grid of smoothly varying rho, T, and user-defined X.  Then
! call react_state() on the grid and output the results.

program test_react
#ifdef CUDA
  use cudafor
  use iso_c_binding, only: c_size_t
  use react_zones_module, only: pfidx_t, react_zones_flat  
#else
  use react_zones_module, only: pfidx_t, react_zones
#endif
  
  use BoxLib
  use bl_constants_module
  use bl_types
  use bl_space
  use f2kcli
  use box_util_module
  use ml_layout_module
  use multifab_module
  use variables, only: init_variables, finalize_variables, plot_t
  use probin_module, only: dens_min, dens_max, &
                           temp_min, temp_max, test_set, run_prefix, &
                           small_temp, small_dens, do_acc, &
                           sort_grid
  use runtime_init_module
  use burn_type_module
  use actual_burner_module, only : actual_burner
  use microphysics_module
  use eos_type_module, only : mintemp, mindens
  use network, only: nspec
  use util_module
  use fabio_module
  use build_info_module
  use parallel, only : parallel_wtime
  use linked_list_module

  implicit none

  ! Conventional fluid state multifabs
  type(multifab) , allocatable :: s(:)

  real(kind=dp_t) :: dx(1, MAX_SPACEDIM)

  logical :: pmask(MAX_SPACEDIM)

  type(ml_layout) :: mla
  type(ml_boxarray) :: mba

  integer :: i, j, n
  integer :: ii, jj, kk, mm
  integer :: nrho, nT, nX

  integer :: dm, nlevs

  integer :: n_rhs_min, n_rhs_max, n_rhs_avg

  type(plot_t) :: pf

  integer :: itemp, irho, ispec, ispec_old, irodot, irho_hnuc

  real(kind=dp_t), pointer :: sp(:,:,:,:)

  real(kind=dp_t), allocatable :: state(:,:,:,:)

#ifdef CUDA
  integer, parameter :: cudaNumStreams = 32
  integer :: streamSize, istate, chunkOffset, chunkSize
  integer(kind=cuda_stream_kind) :: streams(cudaNumStreams)
  logical :: pinnedflag
  integer(c_size_t) :: stacksize
  integer :: cuGrid
  integer :: cuThreadBlock = 64
  type(cell_t) :: trunk
  integer :: flatdim
  real(kind=dp_t), allocatable, managed :: flatstate(:,:)
#endif  

  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: domlo(MAX_SPACEDIM), domhi(MAX_SPACEDIM)

  real (kind=dp_t) :: dlogrho, dlogT
  real (kind=dp_t), allocatable :: xn_zone(:, :)

  real (kind=dp_t) :: sum_X

  real (kind=dp_t) :: start_time, end_time

  character (len=256) :: out_name

  type(pfidx_t), &
#ifdef CUDA
       managed, &
#endif
       allocatable :: pfidx

  
  call boxlib_initialize()
  call bl_prof_initialize(on = .true.)

  call runtime_init(.true.)

  ! initialize a grid -- since we are not doing anything hydroy, set the
  ! periodic mask to true
  pmask(:) = .true.
  call read_a_hgproj_grid(mba, test_set)
  call ml_layout_build(mla, mba, pmask)

  nlevs = mla % nlevel
  if (nlevs /= 1) then
     call bl_error("ERROR: only 1 level of refinement currently supported")
  endif

  dm = mla % dim
  if (dm /= 3) then
     call bl_error("ERROR: we require dm = 3")
  endif

  ! we don't care about dx -- we have no physical size
  dx(1,:) = ONE

  ! microphysics
  call microphysics_init(small_temp=mintemp, small_dens=mindens)

  ! we'll store everything in a multifab -- inputs and outputs
  call init_variables(pf)

  allocate(s(nlevs))

  do n = 1,nlevs
    call multifab_build(s(n), mla%la(n), pf % n_plot_comps, 0)
  end do

  nrho = extent(mla%mba%pd(1),1)
  nT = extent(mla%mba%pd(1),2)
  nX = extent(mla%mba%pd(1),3)

  allocate(state(0:nrho-1, 0:nT-1, 0:nX-1, pf % n_plot_comps))

  dlogrho = (log10(dens_max) - log10(dens_min))/(nrho - 1)
  dlogT   = (log10(temp_max) - log10(temp_min))/(nT - 1)

  ! read from the input file to get all the species data
  domlo = lwb(get_pd(get_layout(s(1))))
  domhi = upb(get_pd(get_layout(s(1))))

  allocate(xn_zone(nspec, 0:nX-1))   ! this assumes that lo(3) = 0

  call get_xn(xn_zone, domlo(3), domhi(3))

  ! normalize -- just in case
  do kk = domlo(3), domhi(3)
     sum_X = sum(xn_zone(:, kk))
     xn_zone(:, kk) = xn_zone(:, kk)/sum_X
  enddo

  ! GPU doesn't like derived-types with bound procedures
  allocate(pfidx)
  pfidx % itemp = pf % itemp
  pfidx % irho = pf % irho
  pfidx % ispec = pf % ispec
  pfidx % ispec_old = pf % ispec_old
  pfidx % irodot = pf % irodot
  pfidx % irho_hnuc = pf % irho_hnuc
  pfidx % endrho = nrho - 1
  pfidx % endT = nT - 1
  pfidx % endX = nX - 1
  pfidx % ncomps = pf % n_plot_comps

  n = 1  ! single level assumption

  n_rhs_avg = 0
  n_rhs_max = -100000000
  n_rhs_min = 100000000

  do i = 1, nfabs(s(n))
     sp => dataptr(s(n), i)

     lo = lwb(get_box(s(n), i))
     hi = upb(get_box(s(n), i))

     ! First, construct the input state in a separate loop.

     do kk = lo(3), hi(3)
        do jj = lo(2), hi(2)
           do ii = lo(1), hi(1)

              state(ii, jj, kk, pf % itemp) = 10.0_dp_t**(log10(temp_min) + dble(jj)*dlogT)
              state(ii, jj, kk, pf % irho) = 10.0_dp_t**(log10(dens_min) + dble(ii)*dlogrho)
              state(ii, jj, kk, pf%ispec_old:pf%ispec_old+nspec-1) = max(xn_zone(:, kk), 1.e-10_dp_t)

           enddo
        enddo
     enddo

     write(*,*) 'lo = ', lo
     write(*,*) 'hi = ', hi
     
     ! Set up a timer for the burn.
     start_time = parallel_wtime()
     
#ifdef CUDA
     ! Uncomment to configure Stack Size Limit
     ! stacksize = 64000
     ! istate = cudaDeviceSetLimit(cudaLimitStackSize, stacksize)
     ! write(*,*) 'limiting stack size to ', stacksize, ' with return code ', istate

     flatdim = get_flatdim_2d(lo(1:2), hi(1:2))
     print *, 'Allocating flatstate with size ', flatdim     
     allocate(flatstate(flatdim, pf % n_plot_comps))

     print *, 'Setting up CUDA'

     ! Set up CUDA parameters
     cuGrid = ceiling(real(flatdim)/cuThreadBlock)
     
     write(*,*) 'cuThreadBlock = ', cuThreadBlock
     write(*,*) 'cuGrid = ', cuGrid

     ! Loop over 2-D slices of the 3-D spatial state
     do kk = lo(3), hi(3)
        print *, 'flattening'
        if (sort_grid .eq. 1) then
           ! Sort 2-D slice from state into flatstate
           call sort_and_flatten_2d(state(lo(1):hi(1), lo(2):hi(2), kk, 1:pf % n_plot_comps), &
                                    lo(1:2), hi(1:2), pf % n_plot_comps, &
                                    pf % itemp, trunk, flatstate)
        else
           ! Directly reshape 2-D slice from state into flatstate
           flatstate = reshape(state(lo(1):hi(1),lo(2):hi(2),kk, 1:pf % n_plot_comps), &
                               [flatdim, pf % n_plot_comps])
        endif

        print *, 'Calling react_zones_flat'
        ! React the zones in flatstate using CUDA kernel
        call react_zones_flat<<<cuGrid,cuThreadBlock>>>(flatstate, &
                                                        pfidx, flatdim)

        ! Save reaction results back to state
        if (sort_grid .eq. 1) then
           call restore_2d(state(lo(1):hi(1), lo(2):hi(2), kk, :), &
                           pf % n_plot_comps, trunk, flatstate)
        else
           state(lo(1):hi(1), lo(2):hi(2), kk, :) = reshape(flatstate, &
                                                           [hi(1)-lo(1)+1, &
                                                            hi(2)-lo(2)+1, &
                                                            pf % n_plot_comps])
        end if

        ! Destroy linked list contents in trunk
        call trunk % delete_children
     end do
#else
     call react_zones(state, lo, hi, pfidx)
#endif

     sp(:,:,:,:) = state(:,:,:,:)

     ! End the timer and print the results.     
     end_time = parallel_wtime()
     
     print *, "Net burn time: ", end_time - start_time

  enddo

  ! note: integer division
  ! n_rhs_avg = n_rhs_avg/(nT*nrho*nX)

  ! print *, "RHS stats:"
  ! print *, "  min: ", n_rhs_min
  ! print *, "  avg: ", n_rhs_avg
  ! print *, "  max: ", n_rhs_max

  ! output
  out_name = trim(run_prefix) // "test_react." // trim(integrator_dir)

  call fabio_ml_multifab_write_d(s, mla%mba%rr(:,1), trim(out_name), names=pf%names)

  call write_job_info(out_name, mla%mba)


  ! if you (or a subroutine) built it, destroy it!
  do n = 1,nlevs
    call destroy(s(n))
  end do

  call destroy(mla)

  call finalize_variables(pf)
  
  deallocate(pfidx)
  deallocate(state)
  deallocate(s)
  deallocate(xn_zone)

  call runtime_close()


  call microphysics_finalize()

  ! end boxlib
  call bl_prof_glean("bl_prof_res")
  call bl_prof_finalize()
  call boxlib_finalize()

  call send_success_return_code()
  
end program test_react
