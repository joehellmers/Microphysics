module linked_list_module

  use bl_types, only: dp_t
  use linked_list_data, only: cell_data_t
  
  implicit none

  type :: cell_t
     type(cell_t), pointer :: next=>NULL(), prev=>NULL()
     type(cell_t), pointer :: first_child=>NULL(), last_child=>NULL()
     type(cell_data_t)     :: data
   contains
     procedure :: new_after  => new_after
     procedure :: new_before => new_before
     procedure :: create_only_child => create_only_child
     procedure :: append_new_child  => append_new_child
     procedure :: prepend_new_child => prepend_new_child
     procedure :: delete_cell       => delete_cell
     procedure :: delete_children   => delete_children
  end type cell_t

contains

  subroutine new_after(cell, cell_new)
    class(cell_t), target :: cell
    type(cell_t), pointer :: cell_new
    type(cell_t), pointer :: cell_ptr_next
    cell_ptr_next => cell % next
    nullify(cell % next)
    allocate(cell % next)
    call cell % next % data % init
    cell_new => cell % next
    cell_new % next => cell_ptr_next
    cell_new % prev => cell
    if ( associated(cell_ptr_next) ) then
       cell_new % next % prev => cell_new
    end if
  end subroutine new_after

  subroutine new_before(cell, cell_new)
    class(cell_t), target :: cell
    type(cell_t), pointer :: cell_new
    type(cell_t), pointer :: cell_ptr_prev
    cell_ptr_prev => cell % prev
    nullify(cell % prev)
    allocate(cell % prev)
    call cell % prev % data % init
    cell_new => cell % prev
    cell_new % prev => cell_ptr_prev
    cell_new % next => cell
    if ( associated(cell_ptr_prev) ) then
       cell_new % prev % next => cell_new
    end if
  end subroutine new_before

  ! NOT TESTED:
  ! subroutine insert_after(cell, cell_insert)
  !   class(cell_t), target :: cell
  !   type(cell_t), target  :: cell_insert
  !   cell_insert % next => cell % next
  !   cell_insert % prev => cell
  !   cell % next => cell_insert
  ! end subroutine insert_after

  ! NOT TESTED:
  ! subroutine insert_before(cell, cell_insert)
  !   class(cell_t), target :: cell
  !   type(cell_t), target  :: cell_insert
  !   cell_insert % next => cell
  !   cell_insert % prev => cell % prev
  !   cell % prev => cell_insert
  ! end subroutine insert_before

  subroutine create_only_child(cell, cell_new)
    ! Creates the first child of an otherwise childless cell
    class(cell_t) :: cell
    type(cell_t), pointer :: cell_new
    allocate(cell % first_child)
    call cell % first_child % data % init
    cell_new => cell % first_child
    cell % last_child => cell_new
  end subroutine create_only_child

  subroutine append_new_child(cell, cell_new)
    class(cell_t) :: cell
    type(cell_t), pointer :: cell_new
    if ( associated(cell % last_child) ) then
       call cell % last_child % new_after(cell_new)
       cell % last_child => cell_new
    else
       call cell % create_only_child(cell_new)
    end if
  end subroutine append_new_child

  subroutine prepend_new_child(cell, cell_new)
    class(cell_t) :: cell
    type(cell_t), pointer :: cell_new
    if ( associated(cell % first_child) ) then
       call cell % first_child % new_before(cell_new)
       cell % first_child => cell_new
    else
       call cell % create_only_child(cell_new)
    end if
  end subroutine prepend_new_child

  subroutine delete_cell(cell, cell_next)
    ! Delete this cell's data and its child cells
    class(cell_t) :: cell
    type(cell_t), pointer :: cell_next
!    write(*,*) 'Entered delete_cell'
    call cell % data % terminate
    if ( associated(cell % prev) ) then
       cell % prev % next => cell % next
    end if
    if ( associated(cell % next) ) then
       cell % next % prev => cell % prev
    end if
    if ( associated(cell % first_child) ) then
       call cell % delete_children
    end if
    cell_next => cell % next
!    write(*,*) 'Leaving delete_cell'
  end subroutine delete_cell

  subroutine delete_children(cell)
    ! Delete and deallocate all child cells
    class(cell_t) :: cell
    type(cell_t), pointer    :: cell_ptr
    type(cell_t), pointer    :: cell_next
!    write(*,*) 'Entered delete_children'
!    write(*,*) ''
    cell_ptr => cell % first_child
    do while ( associated(cell_ptr) )
!       write(*,*) 'Deleting cell cell_ptr'
       call cell_ptr % delete_cell(cell_next)
!       write(*,*) 'Deallocating cell_ptr'
       deallocate(cell_ptr)
!       write(*,*) 'Deallocated cell_ptr'
!       write(*,*) ''
       cell_ptr => cell_next
    end do
  end subroutine delete_children

  
  function get_flatdim_2d(lo, hi) result(flatdim)
    integer :: lo(2), hi(2), flatdim
    flatdim = (hi(1)-lo(1)+1) * (hi(2)-lo(2)+1)
  end function get_flatdim_2d

  
  function get_flatdim_3d(lo, hi) result(flatdim)
    integer :: lo(3), hi(3), flatdim
    flatdim = (hi(1)-lo(1)+1) * (hi(2)-lo(2)+1) * (hi(3)-lo(3)+1)
  end function get_flatdim_3d

  
  subroutine restore_2d(state2d, ncomps, trunk, flatstate)
    type(cell_t), intent(inout) :: trunk
    double precision, intent(inout) :: state2d(:,:,:)
    double precision, intent(in) :: flatstate(:,:)
    integer :: ncomps, ii, jj
    type(cell_t), pointer :: aptr
    ii = 1
    aptr => trunk % first_child
    do while ( associated(aptr) )
       do jj = 1, ncomps
          state2d(aptr % data % ni, aptr % data % nj, jj) = flatstate(ii, jj)
       end do
       aptr => aptr % next
       ii = ii + 1
    end do
  end subroutine restore_2d

    
  subroutine sort_and_flatten_2d(state2d, lo, hi, ncomps, &
                                 isort, trunk, flatstate)
    type(cell_t), intent(inout) :: trunk
    double precision, intent(in) :: state2d(:,:,:)
    double precision, intent(inout) :: flatstate(:,:)
    integer :: lo(3), hi(3), ncomps, ii, jj, isort
    type(cell_t), pointer :: aptr, bptr
    real(kind=dp_t) :: tscratch

    do jj = lo(2), hi(2)
       do ii = lo(1), hi(1)
          ! write(*,*) '(i,j): ', ii, jj
          tscratch = state2d(ii,jj,isort)
          if ( associated(trunk % first_child) ) then
             ! There exists a child cell
             aptr => trunk % first_child
             do while ( associated(aptr) )
                if ( tscratch .ge. aptr % data % T ) then
                   ! Add before current cell
                   call aptr % new_before(bptr)
                   if ( .not. associated(bptr % prev) ) then
                      trunk % first_child => bptr
                   end if
                   bptr % data % T = tscratch
                   bptr % data % ni = ii
                   bptr % data % nj = jj
                   bptr % data % nk = 1
                   ! write(*,*) 'Making new cell BEFORE'
                   exit
                else if ( .not. associated(aptr % next) ) then
                   ! Add after current cell
                   call aptr % new_after(bptr)
                   if ( .not. associated(bptr % next) ) then
                      trunk % last_child => bptr
                   end if
                   bptr % data % T = tscratch
                   bptr % data % ni = ii
                   bptr % data % nj = jj
                   bptr % data % nk = 1
                   ! write(*,*) 'Making new cell AFTER'
                   exit
                else
                   aptr => aptr % next
                end if
             end do
          else
             ! Create first child cell
             call trunk % create_only_child(bptr)
             bptr % data % T = tscratch
             bptr % data % ni = ii
             bptr % data % nj = jj
             bptr % data % nk = 1
             ! write(*,*) 'Created first child cell'
          end if
       end do
    end do

    ! write(*,*) 'Done with filling trunk'

    ! Now fill flatstate
    ii = 1
    aptr => trunk % first_child
    do while ( associated(aptr) )
       do jj = 1, ncomps
          flatstate(ii, jj) = state2d(aptr % data % ni, &
               aptr % data % nj, &
               jj)
       end do
       aptr => aptr % next
       ii = ii + 1
    end do
  end subroutine sort_and_flatten_2d

  
  subroutine sort_and_flatten_3d(state3d, lo, hi, ncomps, &
                                  isort, trunk, flatstate)
    type(cell_t), intent(inout) :: trunk
    double precision, intent(in) :: state3d(:,:,:,:)
    double precision, intent(inout) :: flatstate(:,:)
    integer :: lo(3), hi(3), ncomps, ii, jj, kk, isort
    type(cell_t), pointer :: aptr, bptr
    real(kind=dp_t) :: tscratch

    do kk = lo(3), hi(3)
       do jj = lo(2), hi(2)
          do ii = lo(1), hi(1)
             ! write(*,*) '(i,j,k): ', ii, jj, kk
             tscratch = state3d(ii,jj,kk,isort)
             if ( associated(trunk % first_child) ) then
                ! There exists a child cell
                aptr => trunk % first_child
                do while ( associated(aptr) )
                   if ( tscratch .ge. aptr % data % T ) then
                      ! Add before current cell
                      call aptr % new_before(bptr)
                      if ( .not. associated(bptr % prev) ) then
                         trunk % first_child => bptr
                      end if
                      bptr % data % T = tscratch
                      bptr % data % ni = ii
                      bptr % data % nj = jj
                      bptr % data % nk = kk
                      ! write(*,*) 'Making new cell BEFORE'
                      exit
                   else if ( .not. associated(aptr % next) ) then
                      ! Add after current cell
                      call aptr % new_after(bptr)
                      if ( .not. associated(bptr % next) ) then
                         trunk % last_child => bptr
                      end if
                      bptr % data % T = tscratch
                      bptr % data % ni = ii
                      bptr % data % nj = jj
                      bptr % data % nk = kk
                      ! write(*,*) 'Making new cell AFTER'
                      exit
                   else
                      aptr => aptr % next
                   end if
                end do
             else
                ! Create first child cell
                call trunk % create_only_child(bptr)
                bptr % data % T = tscratch
                bptr % data % ni = ii
                bptr % data % nj = jj
                bptr % data % nk = kk
                ! write(*,*) 'Created first child cell'
             end if
          end do
       end do
    end do

    ! write(*,*) 'Done with filling trunk'

    ! Now fill flatstate
    ii = 1
    aptr => trunk % first_child
    do while ( associated(aptr) )
       do jj = 1, ncomps
          flatstate(ii, jj) = state3d(aptr % data % ni, &
               aptr % data % nj, &
               aptr % data % nk, &
               jj)
       end do
       aptr => aptr % next
       ii = ii + 1
    end do
  end subroutine sort_and_flatten_3d

end module linked_list_module
