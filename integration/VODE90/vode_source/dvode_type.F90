module dvode_type_module

  use cudafor
  use bl_types, only: dp_t
  
  implicit none

  type :: dvode_t
     ! Variables previously in common blocks
     real(dp_t), allocatable, dimension(:) :: ACNRM, CCMXJ, CONP, CRATE, DRC
     real(dp_t), allocatable, dimension(:) :: ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1
     real(dp_t), allocatable, dimension(:) :: RC, RL1, TN, UROUND, HU
     
     ! First dimension is size 13
     real(dp_t), allocatable, dimension(:,:) :: EL, TAU

     ! First dimension is size 5
     real(dp_t), allocatable, dimension(:,:) :: TQ

     integer,    allocatable, dimension(:) :: NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
     integer,    allocatable, dimension(:) :: ICF, INIT, IPUP, JCUR, JSTART, JSV
     integer,    allocatable, dimension(:) :: L, LENWM, KFLAG, KUTH
     integer,    allocatable, dimension(:) :: LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL
     integer,    allocatable, dimension(:) :: NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ
     integer,    allocatable, dimension(:) :: NSLP, MXSTEP

     ! Real parameters (first dimension is size n_rpar_comps)
     real(dp_t), allocatable, dimension(:,:) :: RPAR

     ! State flag
     integer,    allocatable, dimension(:)   :: ISTATE

     ! Local time and integration end time
     real(dp_t), allocatable, dimension(:)   :: T, TOUT

     ! Integration vector (first dimension is size VODE_NEQS)
     real(dp_t), allocatable, dimension(:,:) :: Y

     ! Tolerances (size is VODE_NEQS)
     real(dp_t), allocatable, dimension(:) :: RTOL, ATOL

#ifdef CUDA
     attributes(managed) :: ACNRM, CCMXJ, CONP, CRATE, DRC
     attributes(managed) :: ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1
     attributes(managed) :: RC, RL1, TN, UROUND, HU
     attributes(managed) :: EL, TAU
     attributes(managed) :: TQ
     attributes(managed) :: NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
     attributes(managed) :: ICF, INIT, IPUP, JCUR, JSTART, JSV
     attributes(managed) :: L, LENWM, KFLAG, KUTH
     attributes(managed) :: LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL
     attributes(managed) :: NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ
     attributes(managed) :: NSLP, MXSTEP, RPAR, ISTATE, T, TOUT, Y
     attributes(managed) :: RTOL, ATOL
#endif
  end type dvode_t

contains


  subroutine allocate_dvode_state(dvode_state, size)

    use rpar_indices, only: n_rpar_comps

    implicit none

    type (dvode_t), intent(inout) :: dvode_state
    integer,        intent(in   ) :: size

    ! Variables previously in common blocks
    allocate(dvode_state % ACNRM(size), &
             dvode_state % CCMXJ(size), &
             dvode_state % CONP(size), &
             dvode_state % CRATE(size), &
             dvode_state % DRC(size), &
             dvode_state % ETA(size), &
             dvode_state % ETAMAX(size), &
             dvode_state % H(size), &
             dvode_state % HMIN(size), &
             dvode_state % HMXI(size), &
             dvode_state % HNEW(size), &
             dvode_state % HSCAL(size), &
             dvode_state % PRL1(size), &
             dvode_state % RC(size), &
             dvode_state % RL1(size), &
             dvode_state % TN(size), &
             dvode_state % UROUND(size), &
             dvode_state % HU(size), &
             dvode_state % NCFN(size), &
             dvode_state % NETF(size), &
             dvode_state % NFE(size), &
             dvode_state % NJE(size), &
             dvode_state % NLU(size), &
             dvode_state % NNI(size), &
             dvode_state % NQU(size), &
             dvode_state % NST(size), &
             dvode_state % ICF(size), &
             dvode_state % INIT(size), &
             dvode_state % IPUP(size), &
             dvode_state % JCUR(size), &
             dvode_state % JSTART(size), &
             dvode_state % JSV(size), &
             dvode_state % L(size), &
             dvode_state % LENWM(size), &
             dvode_state % KFLAG(size), &
             dvode_state % KUTH(size), &
             dvode_state % LOCJS(size), &
             dvode_state % MAXORD(size), &
             dvode_state % METH(size), &
             dvode_state % MITER(size), &
             dvode_state % MSBJ(size), &
             dvode_state % MXHNIL(size), &
             dvode_state % NEWH(size), &
             dvode_state % NEWQ(size), &
             dvode_state % NHNIL(size), &
             dvode_state % NQ(size), &
             dvode_state % NQNYH(size), &
             dvode_state % NQWAIT(size), &
             dvode_state % NSLJ(size), &
             dvode_state % NSLP(size), &
             dvode_state % MXSTEP(size))
     
    ! First dimension is size 13
    allocate(dvode_state % EL(13, size), dvode_state % TAU(13, size))

    ! First dimension is size 5
    allocate(dvode_state % TQ(5, size))

    ! Real parameters (first dimension is size n_rpar_comps)
    allocate(dvode_state % RPAR(n_rpar_comps, size))

    ! State flag
    allocate(dvode_state % ISTATE(size))

    ! Local time and integration end time
    allocate(dvode_state % T(size), &
             dvode_state % TOUT(size))

    ! Integration vector (first dimension is size VODE_NEQS)
    allocate(dvode_state % Y(VODE_NEQS, size))

    ! Tolerances (size is VODE_NEQS)
    allocate(dvode_state % RTOL(VODE_NEQS), dvode_state % ATOL(VODE_NEQS))

  end subroutine allocate_dvode_state


  subroutine deallocate_dvode_state(dvode_state)

    implicit none

    type (dvode_t), intent(inout) :: dvode_state

    ! Variables previously in common blocks
    deallocate(dvode_state % ACNRM, &
               dvode_state % CCMXJ, &
               dvode_state % CONP, &
               dvode_state % CRATE, &
               dvode_state % DRC, &
               dvode_state % ETA, &
               dvode_state % ETAMAX, &
               dvode_state % H, &
               dvode_state % HMIN, &
               dvode_state % HMXI, &
               dvode_state % HNEW, &
               dvode_state % HSCAL, &
               dvode_state % PRL1, &
               dvode_state % RC, &
               dvode_state % RL1, &
               dvode_state % TN, &
               dvode_state % UROUND, &
               dvode_state % HU, &
               dvode_state % NCFN, &
               dvode_state % NETF, &
               dvode_state % NFE, &
               dvode_state % NJE, &
               dvode_state % NLU, &
               dvode_state % NNI, &
               dvode_state % NQU, &
               dvode_state % NST, &
               dvode_state % ICF, &
               dvode_state % INIT, &
               dvode_state % IPUP, &
               dvode_state % JCUR, &
               dvode_state % JSTART, &
               dvode_state % JSV, &
               dvode_state % L, &
               dvode_state % LENWM, &
               dvode_state % KFLAG, &
               dvode_state % KUTH, &
               dvode_state % LOCJS, &
               dvode_state % MAXORD, &
               dvode_state % METH, &
               dvode_state % MITER, &
               dvode_state % MSBJ, &
               dvode_state % MXHNIL, &
               dvode_state % NEWH, &
               dvode_state % NEWQ, &
               dvode_state % NHNIL, &
               dvode_state % NQ, &
               dvode_state % NQNYH, &
               dvode_state % NQWAIT, &
               dvode_state % NSLJ, &
               dvode_state % NSLP, &
               dvode_state % MXSTEP)
     
    ! First dimension is size 13
    deallocate(dvode_state % EL, dvode_state % TAU)

    ! First dimension is size 5
    deallocate(dvode_state % TQ)

    ! Real parameters (first dimension is size n_rpar_comps)
    deallocate(dvode_state % RPAR)

    ! State flag
    deallocate(dvode_state % ISTATE)

    ! Local time and integration end time
    deallocate(dvode_state % T, &
               dvode_state % TOUT)

    ! Integration vector (first dimension is size VODE_NEQS)
    deallocate(dvode_state % Y)

    ! Tolerances (size is VODE_NEQS)
    deallocate(dvode_state % RTOL, dvode_state % ATOL)

  end subroutine deallocate_dvode_state


#ifndef CUDA  
  subroutine print_state(dvode_state)
    type(dvode_t) :: dvode_state

    write(*,*) 'HU = ', dvode_state % HU
    write(*,*) 'ACNRM = ', dvode_state % ACNRM
    write(*,*) 'CCMXJ = ', dvode_state % CCMXJ
    write(*,*) 'CONP = ', dvode_state % CONP
    write(*,*) 'CRATE = ', dvode_state % CRATE
    write(*,*) 'DRC = ', dvode_state % DRC
    write(*,*) 'EL(1) = ', dvode_state % EL(1)
    write(*,*) 'EL(2) = ', dvode_state % EL(2)
    write(*,*) 'EL(3) = ', dvode_state % EL(3)
    write(*,*) 'EL(4) = ', dvode_state % EL(4)
    write(*,*) 'EL(5) = ', dvode_state % EL(5)
    write(*,*) 'EL(6) = ', dvode_state % EL(6)
    write(*,*) 'EL(7) = ', dvode_state % EL(7)
    write(*,*) 'EL(8) = ', dvode_state % EL(8)
    write(*,*) 'EL(9) = ', dvode_state % EL(9)
    write(*,*) 'EL(10) = ', dvode_state % EL(10)
    write(*,*) 'EL(11) = ', dvode_state % EL(11)
    write(*,*) 'EL(12) = ', dvode_state % EL(12)
    write(*,*) 'EL(13) = ', dvode_state % EL(13)
    write(*,*) 'ETA = ', dvode_state % ETA
    write(*,*) 'ETAMAX = ', dvode_state % ETAMAX
    write(*,*) 'H = ', dvode_state % H
    write(*,*) 'HMIN = ', dvode_state % HMIN
    write(*,*) 'HMXI = ', dvode_state % HMXI
    write(*,*) 'HNEW = ', dvode_state % HNEW
    write(*,*) 'HSCAL = ', dvode_state % HSCAL
    write(*,*) 'PRL1 = ', dvode_state % PRL1
    write(*,*) 'RC = ', dvode_state % RC
    write(*,*) 'RL1 = ', dvode_state % RL1
    write(*,*) 'TAU(1) = ', dvode_state % TAU(1)
    write(*,*) 'TAU(2) = ', dvode_state % TAU(2)
    write(*,*) 'TAU(3) = ', dvode_state % TAU(3)
    write(*,*) 'TAU(4) = ', dvode_state % TAU(4)
    write(*,*) 'TAU(5) = ', dvode_state % TAU(5)
    write(*,*) 'TAU(6) = ', dvode_state % TAU(6)
    write(*,*) 'TAU(7) = ', dvode_state % TAU(7)
    write(*,*) 'TAU(8) = ', dvode_state % TAU(8)
    write(*,*) 'TAU(9) = ', dvode_state % TAU(9)
    write(*,*) 'TAU(10) = ', dvode_state % TAU(10)
    write(*,*) 'TAU(11) = ', dvode_state % TAU(11)
    write(*,*) 'TAU(12) = ', dvode_state % TAU(12)
    write(*,*) 'TAU(13) = ', dvode_state % TAU(13)
    write(*,*) 'TQ(1) = ', dvode_state % TQ(1)
    write(*,*) 'TQ(2) = ', dvode_state % TQ(2)
    write(*,*) 'TQ(3) = ', dvode_state % TQ(3)
    write(*,*) 'TQ(4) = ', dvode_state % TQ(4)
    write(*,*) 'TQ(5) = ', dvode_state % TQ(5)
    write(*,*) 'TN = ', dvode_state % TN
    write(*,*) 'UROUND = ', dvode_state % UROUND
    write(*,*) 'NCFN = ', dvode_state % NCFN
    write(*,*) 'NETF = ', dvode_state % NETF
    write(*,*) 'NFE = ', dvode_state % NFE
    write(*,*) 'NJE = ', dvode_state % NJE
    write(*,*) 'NLU = ', dvode_state % NLU
    write(*,*) 'NNI = ', dvode_state % NNI
    write(*,*) 'NQU = ', dvode_state % NQU
    write(*,*) 'NST = ', dvode_state % NST
    write(*,*) 'ICF = ', dvode_state % ICF
    write(*,*) 'INIT = ', dvode_state % INIT
    write(*,*) 'IPUP = ', dvode_state % IPUP
    write(*,*) 'JCUR = ', dvode_state % JCUR
    write(*,*) 'JSTART = ', dvode_state % JSTART
    write(*,*) 'JSV = ', dvode_state % JSV
    write(*,*) 'KFLAG = ', dvode_state % KFLAG
    write(*,*) 'KUTH = ', dvode_state % KUTH
    write(*,*) 'L = ', dvode_state % L
    write(*,*) 'LENWM = ', dvode_state % LENWM
    write(*,*) 'LOCJS = ', dvode_state % LOCJS
    write(*,*) 'METH = ', dvode_state % METH
    write(*,*) 'MITER = ', dvode_state % MITER
    write(*,*) 'MSBJ = ', dvode_state % MSBJ
    write(*,*) 'MXHNIL = ', dvode_state % MXHNIL
    write(*,*) 'MXSTEP = ', dvode_state % MXSTEP
    write(*,*) 'NEWH = ', dvode_state % NEWH
    write(*,*) 'NEWQ = ', dvode_state % NEWQ
    write(*,*) 'NHNIL = ', dvode_state % NHNIL
    write(*,*) 'NQ = ', dvode_state % NQ
    write(*,*) 'NQNYH = ', dvode_state % NQNYH
    write(*,*) 'NQWAIT = ', dvode_state % NQWAIT
    write(*,*) 'NSLJ = ', dvode_state % NSLJ
    write(*,*) 'NSLP = ', dvode_state % NSLP
  end subroutine print_state
#endif
  
end module dvode_type_module
