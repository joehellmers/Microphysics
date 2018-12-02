module actual_rhs_module

  use network
  use eos_type_module
  use burn_type_module
  use temperature_integration_module, only: temperature_rhs, temperature_jac
  use sneut_module, only: sneut5
  use actual_network, only: nrates
  use rate_type_module

  implicit none

  double precision, parameter :: c54 = 56.0d0/54.0d0

  ! Table interpolation data

  double precision, parameter :: tab_tlo = 6.0d0, tab_thi = 10.0d0
  integer, parameter :: tab_per_decade = 500
  integer, parameter :: nrattab = int(tab_thi - tab_tlo) * tab_per_decade + 1
  integer, parameter :: tab_imax = int(tab_thi - tab_tlo) * tab_per_decade + 1
  double precision, parameter :: tab_tstp = (tab_thi - tab_tlo) / dble(tab_imax - 1)

  double precision, allocatable :: rattab(:,:)
  double precision, allocatable :: drattabdt(:,:)
  double precision, allocatable :: drattabdd(:,:)
  double precision, allocatable :: ttab(:)

contains

  subroutine actual_rhs_init()

    use aprox_rates_module, only: rates_init
    use extern_probin_module, only: use_tables
    use screening_module, only: screening_init
    use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor

    implicit none

    call rates_init()

    call set_up_screening_factors()

    call screening_init()

    if (use_tables) then

       if (parallel_IOProcessor()) then
          print *, ""
          print *, "Initializing pphotcno rate table"
          print *, ""
       endif

       call create_rates_table()

    endif

  end subroutine actual_rhs_init


  subroutine get_rates(state, rr)

    use extern_probin_module, only: use_tables

    implicit none

    type (burn_t), intent(in) :: state
    type (rate_t), intent(out) :: rr

    double precision :: rho, temp
    double precision :: y(nspec)

    double precision :: ratraw(nrates), dratrawdt(nrates), dratrawdd(nrates)
    double precision :: ratdum(nrates), dratdumdt(nrates), dratdumdd(nrates)
    
	!TODO: What should these be?  Do we really even need them?
	!double precision :: dratdumdy1(irsi2ni:irni2si), dratdumdy2(irsi2ni:irni2si)
    double precision :: scfac(nrates),  dscfacdt(nrates),  dscfacdd(nrates)

    ! Get the data from the state

    rho  = state % rho
    temp = state % T
    y    = state % xn * aion_inv

    ! Get the raw reaction rates
    if (use_tables) then
       call pphotcnotab(temp, rho, ratraw, dratrawdt, dratrawdd)
    else
       call iso7rat(temp, rho, ratraw, dratrawdt, dratrawdd)
    endif

    ! Do the screening here because the corrections depend on the composition
    call screen_iso7(temp, rho, y,                 &
                     ratraw, dratrawdt, dratrawdd, &
                     ratdum, dratdumdt, dratdumdd, &
                     dratdumdy1, dratdumdy2,       &
                     scfac,  dscfacdt,  dscfacdd)

    ! Save the rate data, for the Jacobian later if we need it.

    rr % rates(1,:) = ratdum
    rr % rates(2,:) = dratdumdt
    !rr % rates(3,irsi2ni:irni2si) = dratdumdy1
    !rr % rates(4,irsi2ni:irni2si) = dratdumdy2

    rr % T_eval = temp

  end subroutine get_rates



  subroutine pphotcnotab(btemp, bden, ratraw, dratrawdt, dratrawdd)

    implicit none

    double precision :: btemp, bden, ratraw(nrates), dratrawdt(nrates), dratrawdd(nrates)

    integer, parameter :: mp = 4

	integer          i,j,imax,iat,ifirst
	double precision ye,tlo,thi,tstp,bden_sav,btemp_sav,ye_sav
	double precision x,x1,x2,x3,x4,a,b,c,d,e,f,g,h,p,q
    double precision alfa,beta,gama,delt

    double precision :: dtab(nrates)

    ! Set the density dependence array
    dtab(ircag)  = bden
    dtab(iroga)  = 1.0d0
    dtab(ir3a)   = bden*bden
    dtab(irg3a)  = 1.0d0
    dtab(ir1212) = bden
    dtab(ir1216) = bden
    dtab(ir1616) = bden
    dtab(iroag)  = bden
    dtab(irnega) = 1.0d0
    dtab(irneag) = bden
    dtab(irmgga) = 1.0d0
    dtab(irmgag) = bden
    dtab(irsiga) = 1.0d0
    dtab(ircaag) = bden
    dtab(irtiga) = 1.0d0
    dtab(irsi2ni) = 0.0d0
    dtab(irni2si) = 0.0d0

    ! hash locate
    iat = int((log10(btemp) - tab_tlo)/tab_tstp) + 1
    iat = max(1, min(iat - mp / 2 + 1, tab_imax - mp + 1))

    ! setup the lagrange interpolation coefficients for a cubic
    x  = btemp
    x1 = ttab(iat)
    x2 = ttab(iat+1)
    x3 = ttab(iat+2)
    x4 = ttab(iat+3)
    a  = x - x1
    b  = x - x2
    c  = x - x3
    d  = x - x4
    e  = x1 - x2
    f  = x1 - x3
    g  = x1 - x4
    h  = x2 - x3
    p  = x2 - x4
    q  = x3 - x4
    alfa =  b*c*d/(e*f*g)
    beta = -a*c*d/(e*h*p)
    gama =  a*b*d/(f*h*q)
    delt = -a*b*c/(g*p*q)

    ! crank off the raw reaction rates
    do j = 1, nrates

       ratraw(j) = (alfa * rattab(j,iat) &
                    + beta * rattab(j,iat+1) &
                    + gama * rattab(j,iat+2) &
                    + delt * rattab(j,iat+3) ) * dtab(j)

       dratrawdt(j) = (alfa * drattabdt(j,iat) &
                       + beta * drattabdt(j,iat+1) &
                       + gama * drattabdt(j,iat+2) &
                       + delt * drattabdt(j,iat+3) ) * dtab(j)

       !dratrawdd(j) = alfa * drattabdd(j,iat) &
       !             + beta * drattabdd(j,iat+1) &
       !             + gama * drattabdd(j,iat+2) &
       !             + delt * drattabdd(j,iat+3)

    enddo

    ! hand finish the three body reactions
    !dratrawdd(ir3a) = bden * dratrawdd(ir3a)

  end subroutine pphotcnotab



  subroutine create_rates_table()

    implicit none

    ! Allocate memory for the tables
    allocate(rattab(nrates, nrattab))
    allocate(drattabdt(nrates, nrattab))
    allocate(drattabdd(nrates, nrattab))
    allocate(ttab(nrattab))

    call set_iso7rat()

  end subroutine create_rates_table



  subroutine set_iso7rat()

    implicit none

    double precision :: btemp, bden, ratraw(nrates), dratrawdt(nrates), dratrawdd(nrates)
    integer :: i, j

    bden = 1.0d0

    do i = 1, tab_imax

       btemp = tab_tlo + dble(i-1) * tab_tstp
       btemp = 10.0d0**(btemp)

       call iso7rat(btemp, bden, ratraw, dratrawdt, dratrawdd)

       ttab(i) = btemp

       do j = 1, nrates

          rattab(j,i)    = ratraw(j)
          drattabdt(j,i) = dratrawdt(j)
          drattabdd(j,i) = dratrawdd(j)

       enddo

    enddo

  end subroutine set_iso7rat



  subroutine actual_rhs(state)

    implicit none

    ! This routine sets up the system of ODE's for the iso7
    ! nuclear reactions.  This is an alpha chain + heavy ion network
    ! with (a,p)(p,g) links up to silicon, as well as a Si <-> Ni link.
    !
    ! Isotopes: he4,  c12,  o16,  ne20, mg24, si28, ni56

    type (burn_t)    :: state
    type (rate_t)    :: rr

    logical          :: deriva = .false.

    double precision :: sneut, dsneutdt, dsneutdd, snuda, snudz
    double precision :: enuc

    double precision :: rho, temp, abar, zbar
    double precision :: y(nspec)

    ! Get the data from the state

    rho  = state % rho
    temp = state % T
    abar = state % abar
    zbar = state % zbar
    y(:) = state % xn(:) * aion_inv(:)

    call get_rates(state, rr)

    ! Call the RHS to actually get dydt.

    call rhs(y, rr % rates(1,:), rr % rates(1,:), state % ydot(1:nspec), deriva)

    ! Instantaneous energy generation rate -- this needs molar fractions

    call ener_gener_rate(state % ydot(1:nspec), enuc)

    ! Get the neutrino losses

    call sneut5(temp, rho, abar, zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)

    ! Append the energy equation (this is erg/g/s)

    state % ydot(net_ienuc) = enuc - sneut

    ! Append the temperature equation

    call temperature_rhs(state)

  end subroutine actual_rhs



  ! Analytical Jacobian

  subroutine actual_jac(state)

    use amrex_constants_module, only: ZERO

    implicit none

    type (burn_t)    :: state
    type (rate_t)    :: rr

    logical          :: deriva = .true.

    double precision :: b1, sneut, dsneutdt, dsneutdd, snuda, snudz

    integer          :: j

    double precision :: rho, temp, abar, zbar
    double precision :: y(nspec)

    state % jac(:,:) = ZERO

    call get_rates(state, rr)

    ! Get the data from the state

    rho  = state % rho
    temp = state % T
    abar = state % abar
    zbar = state % zbar
    y    = state % xn * aion_inv

    ! Species Jacobian elements with respect to other species

    call dfdy_isotopes_iso7(y, state % jac(1:nspec,1:nspec), rr % rates(1,:), &
                            rr % rates(3,:), rr % rates(4,:))

    ! Energy generation rate Jacobian elements with respect to species

    do j = 1, nspec
       call ener_gener_rate(state % jac(1:nspec,j), state % jac(net_ienuc,j))
    enddo

    ! Account for the thermal neutrino losses

    call sneut5(temp, rho, abar, zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)

    do j = 1, nspec
       b1 = ((aion(j) - abar) * abar * snuda + (zion(j) - zbar) * abar * snudz)
       state % jac(net_ienuc,j) = state % jac(net_ienuc,j) - b1
    enddo

    ! Evaluate the Jacobian elements with respect to temperature by
    ! calling the RHS using d(ratdum) / dT

    call rhs(y, rr % rates(2,:), rr % rates(1,:), state % jac(1:nspec,net_itemp), deriva)

    call ener_gener_rate(state % jac(1:nspec,net_itemp), state % jac(net_ienuc,net_itemp))
    state % jac(net_ienuc,net_itemp) = state % jac(net_ienuc,net_itemp) - dsneutdt

    ! Temperature Jacobian elements

    call temperature_jac(state)

  end subroutine actual_jac



  ! Evaluates the right hand side of the iso7 ODEs

  subroutine rhs(y, rate, ratdum, dydt, deriva)

    use amrex_constants_module, only: ZERO, SIXTH
    use microphysics_math_module, only: esum5, esum6, esum15

    implicit none

    ! deriva is used in forming the analytic Jacobian to get
    ! the derivative wrt A

    logical          :: deriva
    double precision :: y(nspec), rate(nrates), ratdum(nrates), dydt(nspec)

    ! local variables
    integer          :: i
    double precision :: a(15)

    dydt(1:nspec) = ZERO

    ! set up the system of ode's :
    ! 4he reactions
    a(1)  =  3.0d0 * y(ic12) * rate(irg3a)
    a(2)  = -0.5d0 * y(ihe4) * y(ihe4) * y(ihe4) * rate(ir3a)
    a(3)  =  y(io16) * rate(iroga)
    a(4)  = -y(ic12) * y(ihe4) * rate(ircag)
    a(5)  =  0.5d0 * y(ic12) * y(ic12) * rate(ir1212)
    a(6)  =  0.5d0 * y(ic12) * y(io16) * rate(ir1216)
    a(7)  =  0.5d0 * y(io16) * y(io16) * rate(ir1616)
    a(8)  = -y(io16) * y(ihe4) * rate(iroag)
    a(9)  =  y(ine20) * rate(irnega)
    a(10) =  y(img24) * rate(irmgga)
    a(11) = -y(ine20) * y(ihe4) * rate(irneag)
    a(12) =  y(isi28) * rate(irsiga)
    a(13) = -y(img24) * y(ihe4) * rate(irmgag)
    a(14) = -7.0d0 * rate(irsi2ni) * y(ihe4)
    a(15) =  7.0d0 * rate(irni2si) * y(ini56)

    dydt(ihe4) = esum15(a)

    ! 12c reactions
    a(1) =  SIXTH * y(ihe4) * y(ihe4) * y(ihe4) * rate(ir3a)
    a(2) = -y(ic12) * rate(irg3a)
    a(3) =  y(io16) * rate(iroga)
    a(4) = -y(ic12) * y(ihe4) * rate(ircag)
    a(5) = -y(ic12) * y(ic12) * rate(ir1212)
    a(6) = -y(ic12) * y(io16) * rate(ir1216)

    dydt(ic12) = esum6(a)

    ! 16o reactions
    a(1) = -y(io16) * rate(iroga)
    a(2) =  y(ic12) * y(ihe4) * rate(ircag)
    a(3) = -y(ic12) * y(io16) * rate(ir1216)
    a(4) = -y(io16) * y(io16) * rate(ir1616)
    a(5) = -y(io16) * y(ihe4) * rate(iroag)
    a(6) =  y(ine20) * rate(irnega)

    dydt(io16) = esum6(a)

    ! 20ne reactions
    a(1) =  0.5d0 * y(ic12) * y(ic12) * rate(ir1212)
    a(2) =  y(io16) * y(ihe4) * rate(iroag)
    a(3) = -y(ine20) * rate(irnega)
    a(4) =  y(img24) * rate(irmgga)
    a(5) = -y(ine20) * y(ihe4) * rate(irneag)

    dydt(ine20) = esum5(a)

    ! 24mg reactions
    a(1) =  0.5d0 * y(ic12) * y(io16) * rate(ir1216)
    a(2) = -y(img24) * rate(irmgga)
    a(3) =  y(ine20) * y(ihe4) * rate(irneag)
    a(4) =  y(isi28) * rate(irsiga)
    a(5) = -y(img24) * y(ihe4) * rate(irmgag)

    dydt(img24) = esum5(a)

    ! 28si reactions
    a(1) =  0.5d0 * y(ic12) * y(io16) * rate(ir1216)
    a(2) =  0.5d0 * y(io16) * y(io16) * rate(ir1616)
    a(3) = -y(isi28) * rate(irsiga)
    a(4) =  y(img24) * y(ihe4) * rate(irmgag)
    a(5) = -rate(irsi2ni) * y(ihe4)
    a(6) =  rate(irni2si) * y(ini56)

    dydt(isi28) = esum6(a)

    ! ni56 reactions
    a(1) =  rate(irsi2ni) * y(ihe4)
    a(2) = -rate(irni2si) * y(ini56)

    dydt(ini56) = sum(a(1:2))

  end subroutine rhs



  subroutine iso7rat(btemp, bden, ratraw, dratrawdt, dratrawdd)

    ! this routine generates unscreened
    ! nuclear reaction rates for the iso7 network.

    use tfactors_module
    use aprox_rates_module
    use amrex_constants_module, only: ZERO
    use extern_probin_module, only: use_c12ag_deboer17

    double precision :: btemp, bden
    double precision :: ratraw(nrates), dratrawdt(nrates), dratrawdd(nrates)

    integer          :: i
    double precision :: rrate,drratedt,drratedd
    double precision :: ff1,dff1dt,dff1dd,ff2,dff2dt,dff2dd,tot,dtotdt,dtotdd,invtot
    type (tf_t)      :: tf

    do i=1,nrates
       ratraw(i)    = ZERO
       dratrawdt(i) = ZERO
       dratrawdd(i) = ZERO
    enddo

    if (btemp .lt. 1.0d6) return


    ! get the temperature factors
    call get_tfactors(btemp, tf)

    ! Determine which c12(a,g)o16 rate to use
    if (use_c12ag_deboer17) then
    ! deboer + 2017 c12(a,g)o16 rate
       call rate_c12ag_deboer17(tf,bden, &
                    ratraw(ircag),dratrawdt(ircag),dratrawdd(ircag), &
                    ratraw(iroga),dratrawdt(iroga),dratrawdd(iroga))
    else
    ! 1.7 times cf88 c12(a,g)o16 rate
       call rate_c12ag(tf,bden, &
                    ratraw(ircag),dratrawdt(ircag),dratrawdd(ircag), &
                    ratraw(iroga),dratrawdt(iroga),dratrawdd(iroga))
    endif

    ! triple alpha to c12
    call rate_tripalf(tf,bden, &
                      ratraw(ir3a),dratrawdt(ir3a),dratrawdd(ir3a), &
                      ratraw(irg3a),dratrawdt(irg3a),dratrawdd(irg3a))

    ! c12 + c12
    call rate_c12c12(tf,bden, &
                     ratraw(ir1212),dratrawdt(ir1212),dratrawdd(ir1212), &
                     rrate,drratedt,drratedd)

    ! c12 + o16
    call rate_c12o16(tf,bden, &
                     ratraw(ir1216),dratrawdt(ir1216),dratrawdd(ir1216), &
                     rrate,drratedt,drratedd)

    ! 16o + 16o
    call rate_o16o16(tf,bden, &
                     ratraw(ir1616),dratrawdt(ir1616),dratrawdd(ir1616), &
                     rrate,drratedt,drratedd)

    ! o16(a,g)ne20
    call rate_o16ag(tf,bden, &
                    ratraw(iroag),dratrawdt(iroag),dratrawdd(iroag), &
                    ratraw(irnega),dratrawdt(irnega),dratrawdd(irnega))

    ! ne20(a,g)mg24
    call rate_ne20ag(tf,bden, &
                     ratraw(irneag),dratrawdt(irneag),dratrawdd(irneag), &
                     ratraw(irmgga),dratrawdt(irmgga),dratrawdd(irmgga))

    ! mg24(a,g)si28
    call rate_mg24ag(tf,bden, &
                     ratraw(irmgag),dratrawdt(irmgag),dratrawdd(irmgag), &
                     ratraw(irsiga),dratrawdt(irsiga),dratrawdd(irsiga))

    ! ca40(a,g)ti44
    call rate_ca40ag(tf,bden, &
                     ratraw(ircaag),dratrawdt(ircaag),dratrawdd(ircaag), &
                     ratraw(irtiga),dratrawdt(irtiga),dratrawdd(irtiga))

  end subroutine iso7rat



  subroutine screen_iso7(btemp, bden, y, &
                         ratraw, dratrawdt, dratrawdd, &
                         ratdum, dratdumdt, dratdumdd, &
                         !TODO: Do we really need this?
						 !dratdumdy1, dratdumdy2, &
                         scfac, dscfacdt, dscfacdd)

    use amrex_constants_module, only: ZERO, ONE
    use screening_module, only: screen5, plasma_state, fill_plasma_state
    use tfactors_module

    ! this routine computes the screening factors
    ! and applies them to the raw reaction rates,
    ! producing the final reaction rates used by the
    ! right hand sides and jacobian matrix elements

    double precision :: btemp, bden
    double precision :: y(nspec)
    double precision :: ratraw(nrates), dratrawdt(nrates), dratrawdd(nrates)
    double precision :: ratdum(nrates), dratdumdt(nrates), dratdumdd(nrates)
	!TODO: Do we really need this?    
	!double precision :: dratdumdy1(irsi2ni:irni2si), dratdumdy2(irsi2ni:irni2si)
    double precision :: scfac(nrates),  dscfacdt(nrates),  dscfacdd(nrates)

    integer          :: i, jscr
    double precision :: sc1a,sc1adt,sc1add,sc2a,sc2adt,sc2add, &
                        sc3a,sc3adt,sc3add,abar,zbar,ye,z2bar, &
                        t992,t9i92,yeff_ca40,yeff_ca40dt,yeff_ti44,yeff_ti44dt, &
                        denom,denomdt,denomdd,xx,zz

    type (plasma_state) :: pstate
    type (tf_t)         :: tf

    ! initialize
    do i = 1, nrates
       ratdum(i)     = ratraw(i)
       dratdumdt(i)  = dratrawdt(i)
       dratdumdd(i)  = dratrawdd(i)
       scfac(i)      = ONE
       dscfacdt(i)   = ZERO
       dscfacdd(i)   = ZERO
    enddo

	!TODO: Get rid of these?    
	!dratdumdy1(:) = ZERO
    !dratdumdy2(:) = ZERO

    ! get the temperature factors
    call get_tfactors(btemp, tf)

    ! Set up the state data, which is the same for all screening factors.

    call fill_plasma_state(pstate, btemp, bden, y(1:nspec))

    ! first the always fun triple alpha and its inverse
    jscr = 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)

    jscr = jscr + 1
    call screen5(pstate,jscr,sc2a,sc2adt,sc2add)

    sc3a   = sc1a * sc2a
    sc3adt = sc1adt*sc2a + sc1a*sc2adt
    !sc3add = sc1add*sc2a + sc1a*sc2add

    ratdum(ir3a)    = ratraw(ir3a) * sc3a
    dratdumdt(ir3a) = dratrawdt(ir3a)*sc3a + ratraw(ir3a)*sc3adt
    !dratdumdd(ir3a) = dratrawdd(ir3a)*sc3a + ratraw(ir3a)*sc3add

    scfac(ir3a)     = sc3a
    dscfacdt(ir3a)  = sc3adt
    !dscfacdd(ir3a)  = sc3add


    ! c12 to o16
    jscr = jscr + 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)

    ratdum(ircag)     = ratraw(ircag) * sc1a
    dratdumdt(ircag)  = dratrawdt(ircag)*sc1a + ratraw(ircag)*sc1adt
    !dratdumdd(ircag)  = dratrawdd(ircag)*sc1a + ratraw(ircag)*sc1add

    scfac(ircag)      = sc1a
    dscfacdt(ircag)   = sc1adt
    !dscfacdd(ircag)   = sc1add



    ! c12 + c12
    jscr = jscr + 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)

    ratdum(ir1212)    = ratraw(ir1212) * sc1a
    dratdumdt(ir1212) = dratrawdt(ir1212)*sc1a + ratraw(ir1212)*sc1adt
    !dratdumdd(ir1212) = dratrawdd(ir1212)*sc1a + ratraw(ir1212)*sc1add

    scfac(ir1212)     = sc1a
    dscfacdt(ir1212)  = sc1adt
    !dscfacdd(ir1212)  = sc1add



    ! c12 + o16
    jscr = jscr + 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)

    ratdum(ir1216)    = ratraw(ir1216) * sc1a
    dratdumdt(ir1216) = dratrawdt(ir1216)*sc1a + ratraw(ir1216)*sc1adt
    !dratdumdd(ir1216) = dratrawdd(ir1216)*sc1a + ratraw(ir1216)*sc1add

    scfac(ir1216)     = sc1a
    dscfacdt(ir1216)  = sc1adt
    !dscfacdd(ir1216)  = sc1add



    ! o16 + o16
    jscr = jscr + 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)

    ratdum(ir1216)    = ratraw(ir1216) * sc1a
    dratdumdt(ir1216) = dratrawdt(ir1216)*sc1a + ratraw(ir1216)*sc1adt
    !dratdumdd(ir1216) = dratrawdd(ir1216)*sc1a + ratraw(ir1216)*sc1add

    scfac(ir1216)     = sc1a
    dscfacdt(ir1216)  = sc1adt
    !dscfacdd(ir1216)  = sc1add



    ! o16 to ne20
    jscr = jscr + 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)

    ratdum(iroag)    = ratraw(iroag) * sc1a
    dratdumdt(iroag) = dratrawdt(iroag)*sc1a + ratraw(iroag)*sc1adt
    !dratdumdd(iroag) = dratrawdd(iroag)*sc1a + ratraw(iroag)*sc1add

    scfac(iroag)     = sc1a
    dscfacdt(iroag)  = sc1adt
    !dscfacdd(iroag)  = sc1add



    ! o16 to mg24
    jscr = jscr + 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)

    ratdum(irneag)    = ratraw(irneag) * sc1a
    dratdumdt(irneag) = dratrawdt(irneag)*sc1a + ratraw(irneag)*sc1adt
    !dratdumdd(irneag) = dratrawdd(irneag)*sc1a + ratraw(irneag)*sc1add

    scfac(irneag)     = sc1a
    dscfacdt(irneag)  = sc1adt
    !dscfacdd(irneag)  = sc1add


    ! mg24 to si28
    jscr = jscr + 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)

    ratdum(irmgag)    = ratraw(irmgag) * sc1a
    dratdumdt(irmgag) = dratrawdt(irmgag)*sc1a + ratraw(irmgag)*sc1adt
    !dratdumdd(irmgag) = dratrawdd(irmgag)*sc1a + ratraw(irmgag)*sc1add

    scfac(irmgag)     = sc1a
    dscfacdt(irmgag)  = sc1adt
    !dscfacdd(irmgag)  = sc1add



    ! ca40 to ti44
    jscr = jscr + 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)

    ratdum(ircaag)    = ratraw(ircaag) * sc1a
    dratdumdt(ircaag) = dratrawdt(ircaag)*sc1a + ratraw(ircaag)*sc1adt
    !dratdumdd(ircaag) = dratrawdd(ircaag)*sc1a + ratraw(ircaag)*sc1add

    scfac(ircaag)     = sc1a
    dscfacdt(ircaag)  = sc1adt
    !dscfacdd(ircaag)  = sc1add



    ! the publication, timmes, woosley & hoffman apjs, 129, 377
    ! has a typo on page 393, where its says "y(ic12)+y(io16) .gt. 0.004"
    ! it should be less than or equal to, since the idea is this piece
    ! gets activated during silicon buring, after all the c + o from
    ! oxygen burning is gone.


    if (tf%t9 .gt. 2.5 .and. y(ic12)+y(io16) .le. 4.0d-3) then

       t992  = tf%t972 * tf%t9
       t9i92 = 1.0d0/t992

       yeff_ca40   = t9i92 * exp(239.42*tf%t9i - 74.741)
       yeff_ca40dt = -yeff_ca40*(239.42*tf%t9i2 + 4.5d0*tf%t9i)

       yeff_ti44   = t992  * exp(-274.12*tf%t9i + 74.914)
       yeff_ti44dt = yeff_ti44*(274.12*tf%t9i2 + 4.5d0*tf%t9i)

       denom     = (bden * y(ihe4))**3

       ratdum(irsi2ni)     = yeff_ca40*denom*ratdum(ircaag)*y(isi28)
       !TODO: Get rid of these?
	   !dratdumdy1(irsi2ni) = 3.0d0 * ratdum(irsi2ni)/y(ihe4)
       !dratdumdy2(irsi2ni) = yeff_ca40*denom*ratdum(ircaag)
       dratdumdt(irsi2ni)  = (yeff_ca40dt*ratdum(ircaag) &
            + yeff_ca40*dratdumdt(ircaag))*denom*y(isi28)*1.0d-9
       !dratdumdd(irsi2ni)  = 3.0d0*ratdum(irsi2ni)/bden &
       !     + yeff_ca40*denom*dratdumdd(ircaag)*y(isi28)


       if (denom .ne. 0.0) then

          zz     = 1.0d0/denom
          ratdum(irni2si) = min(1.0d10,yeff_ti44*ratdum(irtiga)*zz)

          if (ratdum(irni2si) .eq. 1.0d10) then
             !TODO: Get rid of this?
			 !dratdumdy1(irni2si) = 0.0d0
             dratdumdt(irni2si)  = 0.0d0
             !dratdumdd(irni2si)  = 0.0d0

          else
             !TODO: Get rid of this?
			 !dratdumdy1(irni2si) = -3.0d0 * ratdum(irni2si)/y(ihe4)
             dratdumdt(irni2si)  = (yeff_ti44dt*ratdum(irtiga) &
                  + yeff_ti44*dratdumdt(irtiga))*zz*1.0d-9
             !dratdumdd(irni2si)  = -3.0d0 * ratdum(irni2si)/bden &
             !     + yeff_ti44*dratdumdd(irtiga)*zz
          end if
       endif
    end if

  end subroutine screen_iso7


  !TODO: do we really need dratdumdy1 and dratdumdy2?
  !subroutine dfdy_isotopes_pphotcno(y,dfdy,ratdum,dratdumdy1,dratdumdy2)
  subroutine dfdy_isotopes_pphotcno(y,dfdy,ratdum)

    use network
    use microphysics_math_module, only: esum3, esum4, esum5, esum7, esum8, &
										esum14, esum20

    implicit none

    ! this routine sets up the dense pphotcno jacobian for the isotopes

    double precision :: y(nspec), dfdy(nspec,nspec)
    double precision :: ratdum(nrates)

	!TODO: Do we really need this?	
	!double precision :: dratdumdy1(irsi2ni:irni2si), dratdumdy2(irsi2ni:irni2si)

    double precision :: b(20)


	!!! Hydrogen-1 Elements

    ! d(h1)/d(h1)
    b(1) = -2.0d0*y(ih1)*(ratdum(irpp) + ratdum(irpep))
    b(2) = -y(ih2)*ratdum(irdpg)
    b(3) = -y(ili7)*ratdum(irli7pa) 
    b(4) = -y(ibe7)*ratdum(irbepg)
    b(5) = -y(ihe3)*ratdum(irhep)
    dfdy(ih1,ih1) = esum5(b)

    ! d(h1)/d(h2)
    b(1) =  -y(ih1)*ratdum(irdpg)
    dfdy(ih1,ih2) = b(1)

    ! d(h1)/d(he3)
    b(1) =  2.0d0*y(ihe3)*ratdum(ir33)
    b(2) =  -y(ih1)*ratdum(irhep)
    dfdy(ih1,ihe3) = sum(b(1:2))

    ! d(h1)/d(li7)
    b(1) =  -y(ih1)*ratdum(irli7pa)
	dfdy(ih1,ili7) = b(1)

    ! d(h1)/d(be7)
    b(1) =  -y(ih1)*ratdum(irbepg)
	dfdy(ih1,ibe7) = b(1)

    ! d(h1)/d(b8)
    b(1) =  ratdum(irb8gp)
	dfdy(ih1,ib8) = b(1)


	!!! Hydrogen-2 elements
	! d(h2)/d(h1)    
	b(1) =  -y(ih2)*ratdum(irdpg)
	b(2) =  +y(ih1)*(ratdum(irpp) + ratdum(irpep))
    b(3) =  -y(ih1)*ratdum(irdpg)
	dfdy(ih2,ih1) = esum3(b)

	!!! Helium-3 elements
	! d(he3)/d(h1)    
	b(1) = y(ih2)*ratdum(irdpg)
	b(2) = -y(ihe3)*ratdum(irhep)
    dfdy(ihe3,ih1) = sum(b(1:2))
  
	! d(he3)/d(h2)	
	b(1) = y(ih1)*ratdum(irdpg)
	dfdy(ihe3,ih2) = b(1)
	
	!TODO: What is this about? He-3 going to He-3?
	! d(he3)/d(he3)
	b(1) = -2.0d0*y(ihe3)*ratdum(ir33)
    b(2) = -y(ihe4)*ratdum(irhe3ag)
	b(3) = -y(ih1)*ratdum(irhep)
	dfdy(ihe3,ihe3) = esum3(b)

	! d(he3)/d(he4)
	b(1) = -y(ihe3)*ratdum(irhe3ag)
	dfdy(ihe3,ihe4) = b(1)

	!!! Helium-4 elements
	! d(he4)/d(h1)
    b(1) = 2.0d0*y(ili7)*ratdum(irli7pa)
    b(2) = +y(ihe3)*ratdum(irhep)
	dfdy(ihe4,ih1) = sum(b(1:2))

	! d(he4)/d(he3)
    b(1) = y(ihe3)*ratdum(ir33)
	b(2) = -y(ihe4)*ratdum(irhe3ag)
	b(3) = +y(ih1)*ratdum(irhep)
	dfdy(ihe4,ihe3) = esum3(b)      

	! d(he4)/d(he4)
	b(1) = -y(ihe3)*ratdum(irhe3ag)
	dfdy(ihe4,ihe4) = b(1)

	! d(he4)/d(li7)
    b(1) = 2.0d0*y(ih1)*ratdum(irli7pa)
	dfdy(ihe4,ili7) = b(1)

	! d(he4)/d(b8)
    b(1) = 2.0d0*ratdum(irb8ep)
	dfdy(ihe4,ib8) = b(1)

	!!! Lithium-7 elements
	! d(li7)/d(h1)
    b(1) = -y(ili7)*ratdum(irli7pa)
    dfdy(ili7,ih1) = b(1)

	! d(li7)/d(li7)
	b(1) = -y(ih1)*ratdum(irli7pa)
	dfdy(ili7,ili7) = b(1) 

	! d(li7)/d(be7)
	b(1) = ratdum(irbeec)
	dfdy(ili7,ibe7) = b(1)

	!!! Beryllium-7 elements
	! d(be7)/d(h1)
	b(1) = -y(ibe7)*ratdum(irbepg)
	dfdy(ibe7,ih1) = b(1)

	! d(be7)/d(he3)
    b(1) = y(ihe4)*ratdum(irhe3ag)
	dfdy(ibe7,ihe3) = b(1)

	! d(be7)/d(he4)
    b(1) = y(ihe3)*ratdum(irhe3ag)
	dfdy(ibe7,ihe4) = b(1)

	! d(be7)/d(be7)
	b(1) = -y(ih1)*ratdum(irbepg)
    b(2) = -ratdum(irbeec)
	dfdy(ibe7,ibe7) = sum(b(1:2))

	! d(be7)/d(b8)
    b(1) = ratdum(irb8gp)
	dfdy(ibe7,ib8) = b(1)

	!!! Boron-8 elements
	! d(b8)/d(h1)
	b(1) = y(ibe7)*ratdum(irbepg)
	dfdy(ib8,ih1) = b(1)

	! d(b8)/d(be7)
	b(1) = y(ih1)*ratdum(irbepg)
	dfdy(ib8,ibe7) = b(1)
	
	! d(b8)/d(b8)
	b(1) = -ratdum(irb8ep)
    b(2) = -ratdum(irb8gp)
	dfdy(ib8,ib8) = sum(b(1:2))

!!! HOT CNO CONTRIBUTIONS !!!

	!!! Hot CNO Hydrogen-1

	! d(h1)/d(h1)
	b(1) = dfdy(ih1,ih1)	
	b(2) = -y(ic12)*ratdum(irc12pg)
	b(3) = -y(ic13)*ratdum(irc13pg)
	b(4) = -y(in14)*ratdum(irn14pg)
	b(5) = -y(in15)*ratdum(irn15pa)
	b(6) = -y(in15)*ratdum(irn15pg)
	b(7) = -y(io16)*ratdum(iro16pg)
	b(8) = -y(io17)*ratdum(iro17pa)
	b(9) = -y(io17)*ratdum(iro17pg)
	b(10) = -y(io18)*ratdum(iro18pa)
	b(11) = -y(io18)*ratdum(iro18pg)
	b(12) = -y(if19)*ratdum(irf19pa)
	b(13) = -y(in13)*ratdum(irn13pg)
	b(14) = -y(if17) * ratdum(irf17pa)
	b(15) = -y(if17) * ratdum(irf17pg)
	b(16) = -y(if18) * ratdum(irf18pa)
	b(17) = -3.0d0*y(ine19)*ratdum(irne19pg)
	b(18) = -2.0d0*y(ine20)*ratdum(irne20pg)
	b(19) = -8.0d0 * y(img22)*ratdum(iralam1)*(1.0d0 - ratdum(irdelta1))
	b(20) = -26.0d0 * y(is30)*ratdum(iralam2)*(1.0d0 - ratdum(irdelta2))

	dfdy(ih1,ih1) = esum20(b)

	! d(h1)/d(he4)
    b(1) = dfdy(ih1,ihe4)
	b(2) = +y(ic12)*ratdum(irc12ap)
    b(3) = +y(in14)*ratdum(irn14ap)
    b(4) = +y(in15)*ratdum(irn15ap)
    b(5) = +y(io16)*ratdum(iro16ap)
    b(6) = +y(io14)*ratdum(iro14ap)
    b(7) = +y(io15)*ratdum(iro15ap)
	dfdy(ih1,ihe4) = esum7(b)

	! d(h1)/d(c12)	
	b(1) = -y(ih1)*ratdum(irc12pg)
	b(2) = +y(ihe4)*ratdum(irc12ap)
	dfdy(ih1,ic12) = sum(b(1:2))

	! d(h1)/d(c13)
	b(1) = -y(ih1)*ratdum(irc13pg)
	dfdy(ih1,ic13) = b(1)
	
	! d(h1)/d(n13)    
	b(1) = ratdum(irn13gp)
	b(2) = -y(ih1)*ratdum(irn13pg)
	dfdy(ih1,in13) = sum(b(1:2))

	! d(h1)/d(n14)
	b(1) = ratdum(irn14gp)
	b(2) = -y(ih1)*ratdum(irn14pg)
	b(3) = +y(ihe4)*ratdum(irn14ap)
	dfdy(ih1,in14) = esum3(b)

	! d(h1)/d(n15)    
	b(1) = -y(ih1)*ratdum(irn15pa)
	b(2) = -y(ih1)*ratdum(irn15pg)
	b(3) = +y(ihe4)*ratdum(irn15ap)
	dfdy(ih1,in15) = esum3(b)

	! d(h1)/d(o14)
    b(1) = y(ihe4) * ratdum(iro14ap)
	b(2) = +ratdum(iro14gp)
	dfdy(ih1,io14) = sum(b(1:2))

	! d(h1)/d(o15)
	b(1) = ratdum(iro15gp)
	b(2) = +y(ihe4) * ratdum(iro15ap)
    dfdy(ih1,io15) = sum(b(1:2))

	! d(h1)/d(o16)
    b(1) = ratdum(iro16gp)
	b(2) = -y(ih1)*ratdum(iro16pg)
	b(3) = +y(ihe4)*ratdum(iro16ap)
	dfdy(ih1,io16) = esum3(b)

	! d(h1)/d(o17)
	b(1) = -y(ih1)*ratdum(iro17pa)
	b(2) = -y(ih1)*ratdum(iro17pg)
    dfdy(ih1,io17) = sum(b(1:2))

	! d(h1)/d(o18)
    b(1) = -y(ih1)*ratdum(iro18pa)
	b(2) = -y(ih1)*ratdum(iro18pg)
	dfdy(ih1,io18) = sum(b(1:2))

	! d(h1)/d(f17)
    b(1) = ratdum(irf17gp)
	b(2) = -y(ih1) * ratdum(irf17pa)
	b(3) = -y(ih1) * ratdum(irf17pg)
	dfdy(ih1,if17) = esum3(b)

	! d(h1)/d(f18)
	b(1) = ratdum(irf18gp)
	b(2) = -y(ih1) * ratdum(irf18pa)
	dfdy(ih1,if18) = sum(b(1:2))

	! d(h1)/d(f19)	
	b(1) = ratdum(irf19gp)
	b(2) = -y(ih1)*ratdum(irf19pa)
	dfdy(ih1,if19) = sum(b(1:2))

	! d(h1)/d(ne18)
    b(1) = ratdum(irne18gp)
	dfdy(ih1,ine18) = b(1)

	! d(h1)/d(ne19)  
    b(1) = -3.0d0*y(ih1)*ratdum(irne19pg)
	dfdy(ih1,ine19) = b(1)

	! d(h1)/d(ne20)
	b(1) = -2.0d0*y(ih1)*ratdum(irne20pg)
	dfdy(ih1,ine20) = b(1)

	! d(h1)/d(mg22)
    b(1) = -8.0d0 * y(ih1) * ratdum(iralam1)*(1.0d0 - ratdum(irdelta1))
	dfdy(ih1,img22) = b(1)

	! d(h1)/d(s30)
	b(1) = -26.0d0 * y(ih1) * ratdum(iralam2)*(1.0d0 - ratdum(irdelta2))
    dfdy(ih1,is30) = b(1)

!!! Hot CNO Helium-4

	! d(he4)/d(h1)
    b(1) = dfdy(ihe4,ih1)
	b(2) = +y(in15) * ratdum(irn15pa)
    b(3) = +y(io17) * ratdum(iro17pa)
    b(4) = +y(io18) * ratdum(iro18pa)
    b(5) = +y(if19) * ratdum(irf19pa)
    b(6) = +y(if17) * ratdum(irf17pa)
    b(7) = +y(if18) * ratdum(irf18pa)
	dfdy(ihe4,ih1) = esum7(b)

	! d(he4)/d(he4)
    b(1) =   dfdy(ihe4,ihe4)
    b(2) = -y(ic12) * ratdum(irc12ap)
    b(3) = -y(in14) * ratdum(irn14ap)
    b(4) = -y(in15) * ratdum(irn15ap)
    b(5) = -y(io16) * ratdum(iro16ap)
    b(6) = -y(io14) * ratdum(iro14ap)
    b(7) = -y(io15) * ratdum(iro15ap)
    b(8) = -1.5d0 * y(ihe4) * y(ihe4) * ratdum(ir3a)
    b(9) = -y(io16) * ratdum(iroag)
    b(10) = -2.0d0*y(img22) * ratdum(iralam1) * ratdum(irdelta1)
    b(11) = -6.5d0*y(is30) * ratdum(iralam2) * ratdum(irdelta2)
    b(12) = -y(ic12) * ratdum(ircag)
    b(13) = -y(ine18) * ratdum(irne18ap)
    b(14) = -y(io15) * ratdum(iro15ag)
	dfdy(ihe4,ihe4) = esum14(b)

	! d(he4)/d(c12)      
	b(1) = -y(ihe4)*ratdum(irc12ap)
    b(2) = +3.0d0 * ratdum(irg3a)
    b(3) = -y(ihe4) * ratdum(ircag)
	dfdy(ihe4,ic12) = esum3(b)

	! d(he4)d(in14)
    b(1) = -y(ihe4)*ratdum(irn14ap)
	dfdy(ihe4,in14) = b(1)

	! d(he4)/d(n15)
    b(1) = y(ih1)*ratdum(irn15pa) - y(ihe4)*ratdum(irn15ap)
	dfdy(ihe4,in15) = b(1)
	
	! d(he4)/d(o14)    
	b(1) = -y(ihe4) * ratdum(iro14ap)
	dfdy(ihe4,io14) = b(1)

	! d(he4)/d(o15)
	b(1) = -y(ihe4) * ratdum(iro15ap)
	b(2) = -y(ihe4) * ratdum(iro15ag)
	dfdy(ihe4,io15) = b(2)

	! d(he4)/d(o16)
    b(1) = -y(ihe4)*ratdum(iro16ap)
    b(2) = -y(ihe4) * ratdum(iroag)
    b(3) = +ratdum(iroga)
	dfdy(ihe4,io16) = esum3(b)

	! d(he4)/d(o17)    
  	b(1) = y(ih1)*ratdum(iro17pa)
	dfdy(ihe4,io17) = b(1)

	! d(he4)/d(o18)
    b(1) = y(ih1)*ratdum(iro18pa)
	dfdy(ihe4,io18) = b(1)

	! d(he4)/d(f17)
    b(1) = y(ih1) * ratdum(irf17pa)
	dfdy(ihe4,if17) = b(1)

	! d(he4)/d(f18)
    b(1) = y(ih1) * ratdum(irf18pa)
	dfdy(ihe4,if18) = b(1)

	! d(he4)/d(f19)
    b(1) = y(ih1)*ratdum(irf19pa)
	dfdy(ihe4,if19) = b(1)

	! d(he4)/d(ne18)
    b(1) = -y(ihe4) * ratdum(irne18ap)
	dfdy(ihe4,ine18) = b(1)
	
	! d(he4)/d(mg22)
    b(1) = -2.0d0*y(ihe4)*ratdum(iralam1)*ratdum(irdelta1)
	dfdy(ihe4,img22) = b(1)

	! d(he4)/d(s30)
    b(1) = -6.5d0*y(ihe4)*ratdum(iralam2)*ratdum(irdelta2)
	dfdy(ihe4,is30) = b(1)

!!! Hot CNO Carbon-12 elements
	
	! d(c12)/d(h1)    
  	b(1) = -y(ic12)*ratdum(irc12pg)
	b(2) = +y(in15)*ratdum(irn15pa)
	dfdy(ic12,ih1) = sum(b(1:2))

	! d(c12)/d(he4)
    b(1) = -y(ic12)*ratdum(irc12ap)
    b(2) = +0.5d0 * y(ihe4) * y(ihe4) * ratdum(ir3a)
    b(3) = -y(ic12) * ratdum(ircag)
	dfdy(ic12,ihe4) = esum3(b)

	! d(c12)/d(c12)
    b(1) = -y(ih1)*ratdum(irc12pg)
    b(2) = -y(ihe4)*ratdum(irc12ap)
    b(3) = -ratdum(irg3a)
    b(4) = -y(ihe4) * ratdum(ircag)
	dfdy(ic12,ic12) = esum4(b)
    
	! d(c12)/d(n13)
  	b(1) = ratdum(irn13gp)
	dfdy(ic12,in13) = b(1)

	! d(c12)/d(n15)
    b(1) = y(ih1)*ratdum(irn15pa)
	dfdy(ic12,in15) =  b(1)

	! d(c12)/d(o16)
    b(1) = ratdum(iroga)
	dfdy(ic12,io16) = b(1)

!!! Hot CNO Carbon-13 elements
	
	! d(c13)/d(h1)
    b(1) = -y(ic13)*ratdum(irc13pg)
	dfdy(ic13,ih1) = b(1)

	! d(c13)/d(c13)
	b(1) = -y(ih1)*ratdum(irc13pg)
	dfdy(ic13,ic13) = b(1)

	! d(c13)/d(n13)
    b(1) = ratdum(irn13enu)
	dfdy(ic13,in13) = b(1)

	! d(c13)/d(n14)
	b(1) = ratdum(irn14gp)	    
	dfdy(ic13,in14) = b(1)

!!! Hot CNO Nitrogen-13
	! d(n13)/d(h1)
    b(1) = y(ic12)*ratdum(irc12pg)
    b(2) = -y(in13)*ratdum(irn13pg)
	dfdy(in13,ih1) = sum(b(1:2))

	! d(n13)/d(c12)      
	b(1) = y(ih1)*ratdum(irc12pg)
	dfdy(in13,ic12) = b(1)

	! d(n13)/d(n13)
    b(1) = -ratdum(irn13gp)
    b(2) = -ratdum(irn13enu)
    b(3) = -y(ih1)*ratdum(irn13pg)
	dfdy(in13,in13) = esum3(b)

	! d(n13)/d(o14)
    b(1) = ratdum(iro14gp)
	dfdy(in13,io14) = b(1)

! Hot CNO Nitrogen-14

	! d(n14)/d(h1)
    b(1) = y(ic13)*ratdum(irc13pg)
    b(2) = -y(in14)*ratdum(irn14pg)
    b(3) = +y(io17)*ratdum(iro17pa)
	dfdy(in14,ih1) = esum3(b)

	! d(n14)/d(he4)
	b(1) = -y(in14)*ratdum(irn14ap)
	dfdy(in14,ihe4) = b(1)

	! d(n14)/d(c13)
    b(1) = y(ih1)*ratdum(irc13pg)
	dfdy(in14,ic13) = b(1)
	
	! d(n14)/d(n14)
    b(1) = -ratdum(irn14gp)
    b(2) = -y(ih1)*ratdum(irn14pg)
    b(3) = -y(ihe4)*ratdum(irn14ap)
	dfdy(in14,in14) = esum3(b)
	
	! d(n14)/d(o14)
    b(1) = ratdum(iro14enu)
	dfdy(in14,io14) = b(1)

	! d(n14)/d(o15)
    b(1) = ratdum(iro15gp)
	dfdy(in14,io15) = b(1)

	! d(n14)/d(o17)
    b(1) = y(ih1)*ratdum(iro17pa)
	dfdy(in14,io17) = b(1)

!!! Hot CNO Nitrogen-15 elements

	! d(n15)/d(h1)
    b(1) = -y(in15)*ratdum(irn15pa)
    b(2) = -y(in15)*ratdum(irn15pg)
    b(3) = +y(io18)*ratdum(iro18pa)
	dfdy(in15,ih1) = esum3(b)


	! d(n15)/d(he4)
    b(1) = y(ic12)*ratdum(irc12ap)
    b(2) = -y(in15)*ratdum(irn15ap)
	dfdy(in15,ihe4) = sum(b(1:2))

	! d(n15)/d(c12)
    b(1) = y(ihe4)*ratdum(irc12ap)
	dfdy(in15,ic12) =  b(1)

	! d(n15)/d(n15)
    b(1) = -y(ih1)*ratdum(irn15pa)
    b(2) = -y(ih1)*ratdum(irn15pg)
    b(3) = -y(ihe4)*ratdum(irn15ap)
  	dfdy(in15,in15) = esum3(b)
 
	! d(n15)/d(o15)
    b(1) = ratdum(iro15enu)
	dfdy(in15,io15) = b(1)
	
	! d(n15)/d(o16)
    b(1) = ratdum(iro16gp)
	dfdy(in15,io16) = b(1)  
	
	! d(n15)/d(o18)
	b(1) = y(ih1)*ratdum(iro18pa)
	dfdy(in15,io18) = b(1)  

!!! Hot CNO Oxygen-14 elements

	! d(o14)/d(h1)
    b(1) = y(in13)*ratdum(irn13pg)
    b(2) = +y(if17) * ratdum(irf17pa)
	dfdy(io14,ih1) = sum(b(1:2))

	! d(o14)/d(he4)    
	b(1) = -y(io14) * ratdum(iro14ap)
	dfdy(io14,ihe4) = b(1)

	! d(o14)/d(n13)    
    b(1) = y(ih1)*ratdum(irn13pg)
	dfdy(io14,in13) = b(1)
	
	! d(o14)/d(o14)
	b(1) = -ratdum(iro14gp)
    b(2) = -ratdum(iro14enu)
    b(3) = -y(ihe4) * ratdum(iro14ap)
	dfdy(io14,io14) = esum3(b)
	
	! d(o14)/d(f17)    
    b(1) = y(ih1) * ratdum(irf17pa)
	dfdy(io14,if17) = b(1) 

!!! Hot CNO Oxygen-15 elements
	
	! d(o15)/d(h1)
    b(1) = y(in14)*ratdum(irn14pg)
    b(2) = +y(if18) * ratdum(irf18pa)
	dfdy(io15,ih1) = sum(b(1:2))
  
	! d(o15)/d(he4)
    b(1) = -y(io15) * ratdum(iro15ap)
    b(2) = -y(io15) * ratdum(iro15ag)
	dfdy(io15,ihe4) = sum(b(1:2))

	! d(o15)/d(n14)
    b(1) = y(ih1)*ratdum(irn14pg)
	dfdy(io15,in14) = b(1)
 
	! d(o15)/d(o15)
    b(1) = -ratdum(iro15gp)
    b(2) = -ratdum(iro15enu)
    b(3) = -y(ihe4) * ratdum(iro15ap)
    b(4) = -y(ihe4) * ratdum(iro15ag)
	dfdy(io15,io15) = esum4(b)

	! d(o15)/d(f18)	
	b(1) = y(ih1) * ratdum(irf18pa)
	dfdy(io15,if18) = b(1)  


!!! Hot CNO Oxygen-16 elements
	
	! d(o16)/d(h1)
    b(1) = y(in15)*ratdum(irn15pg)
    b(2) = -y(io16)*ratdum(iro16pg)
    b(3) = +y(if19)*ratdum(irf19pa)
	dfdy(io16,ih1) = esum3(b)

	! d(o16)/d(he4)
    b(1) = -y(io16)*ratdum(iro16ap)
    b(2) = -y(io16) * ratdum(iroag)
    b(3) = +y(ic12) * ratdum(ircag)	
	dfdy(io16,ihe4) = esum3(b)

	! d(o16)/d(c12)
    b(1) = y(ihe4) * ratdum(ircag)
	dfdy(io16,ic12) = b(1)  

	! d(o16)/d(n15)
    b(1) = y(ih1)*ratdum(irn15pg)
  	dfdy(io16,in15) = b(1)
  
	! d(o16)/d(o16)
    b(1) = -ratdum(iro16gp)
    b(2) = -y(ih1)*ratdum(iro16pg)
    b(3) = -y(ihe4)*ratdum(iro16ap)
    b(4) = -y(ihe4) * ratdum(iroag)
    b(5) = -ratdum(iroga)
	dfdy(io16,io16) = esum5(b)

	! d(o16)/d(h1)
    b(1) = ratdum(irf17gp)
	dfdy(io16,if17) = b(1)

	! d(o16)/d(h1)
    b(1) = y(ih1)*ratdum(irf19pa)
	dfdy(io16,if19) = b(1)

!!! Hot CNO Oxygen-17 elements

	! d(io17)/d(h1)      
	b(1) = -y(io17)*ratdum(iro17pa)
    b(2) = -y(io17)*ratdum(iro17pg)
	dfdy(io17,ih1) = sum(b(1:2))

	! d(io17)/d(he4)      
    b(1) = y(in14)*ratdum(irn14ap)
	dfdy(io17,ihe4) = b(1)

	! d(io17)/d(n14)      
    b(1) = y(ihe4)*ratdum(irn14ap)
	dfdy(io17,in14) = b(1)

	! d(io17)/d(o17)      
    b(1) = -y(ih1)*ratdum(iro17pa)
	b(2) = -y(ih1)*ratdum(iro17pg)
	dfdy(io17,io17) = sum(b(1:2))

	! d(io17)/d(f17)      
    b(1) = ratdum(irf17enu)
	dfdy(io17,if17) = b(1)

	! d(io17)/d(f18)      
    b(1) = ratdum(irf18gp)
	dfdy(io17,if18) = b(1)

!!! Hot CNO Oxygen-18 elements

	! d(o18)/d(h1)      
	b(1) = -y(io18)*ratdum(iro18pa)
    b(2) = -y(io18)*ratdum(iro18pg)
	dfdy(io18,ih1) = sum(b(1:2))

	! d(o18)/d(he4)      
	b(1) = y(in15)*ratdum(irn15ap)
	dfdy(io18,ihe4) = b(1)

	! d(o18)/d(n15)      
    b(1) = y(ihe4)*ratdum(irn15ap)
	dfdy(io18,in15) = b(1)

	! d(o18)/d(o18)      
    b(1) = -y(ih1)*ratdum(iro18pa)
    b(2) = -y(ih1)*ratdum(iro18pg)
	dfdy(io18,io18) = sum(b(1:2))

	! d(o18)/d(f18)      
    b(1) = ratdum(irf18enu)
	dfdy(io18,if18) = b(1)

	! d(o18)/d(f19)      
    b(1) = ratdum(irf19gp)
	dfdy(io18,if19) = b(1)

!!! Hot CNO Fluorine-17 elements

	! d(f17)/d(h1)
    b(1) = y(io16) * ratdum(iro16pg)
    b(2) = -y(if17) * ratdum(irf17pa)
    b(3) = -y(if17) * ratdum(irf17pg)
	dfdy(if17,ih1) = esum3(b)

	! d(f17)/d(he4)
    b(1) = y(io14) * ratdum(iro14ap)
	dfdy(if17,ihe4) = b(1)

	! d(f17)/d(o14)
    b(1) = y(ihe4) * ratdum(iro14ap)
	dfdy(if17,io14) = b(1)

	! d(f17)/d(o16)
    b(1) = y(ih1)*ratdum(iro16pg)
	dfdy(if17,io16) = b(1)

	! d(f17)/d(f17)
    b(1) = -ratdum(irf17gp)
    b(2) = -ratdum(irf17enu)
    b(3) = -y(ih1) * ratdum(irf17pa)
    b(4) = -y(ih1) * ratdum(irf17pg)
	dfdy(if17,if17) = esum4(b)

	! d(f17)/d(ne18)
    b(1) = ratdum(irne18gp)
	dfdy(if17,ine18) = b(1)

!!! Hot CNO Fluorine-18 elements

	! d(f18)/d(h1)    
	b(1) = y(io17)*ratdum(iro17pg)
    b(2) = -y(if18) * ratdum(irf18pa)
	dfdy(if18,ih1) = sum(b(1:2))

	! d(f18)/d(he4)    
    b(1) = y(io15) * ratdum(iro15ap)
	dfdy(if18,ihe4) = b(1)

	! d(f18)/d(o15)    
    b(1) = y(ihe4) * ratdum(iro15ap)
	dfdy(if18,io15) = b(1)

	! d(f18)/d(o17)    
    b(1) = y(ih1)*ratdum(iro17pg)
	dfdy(if18,io17) = b(1)

	! d(f18)/d(f18)    
    b(1) = -ratdum(irf18gp)
    b(2) = -ratdum(irf18enu)
    b(3) = -y(ih1) * ratdum(irf18pa)
	dfdy(if18,if18) = esum3(b)

	! d(f18)/d(ne18)    
    b(1) = ratdum(irne18enu)
	dfdy(if18,ine18) = b(1)

!!! Hot CNO Fluorine-19 elements

	! d(f19)/d(h1)      
	b(1) = y(io18)*ratdum(iro18pg)
    b(2) = -y(if19)*ratdum(irf19pa)
	dfdy(if19,ih1) = sum(b(1:2))

	! d(f19)/d(he4)      
    b(1) = y(io16)*ratdum(iro16ap)
	dfdy(if19,ihe4) = b(1)

	! d(f19)/d(o16)      
    b(1) = y(ihe4)*ratdum(iro16ap)
	dfdy(if19,io16) = b(1)

	! d(f19)/d(o18)      
    b(1) = y(ih1)*ratdum(iro18pg)
	dfdy(if19,io18) = b(1)

	! d(f19)/d(f19)      
    b(1) = -ratdum(irf19gp)
    b(2) = -y(ih1)*ratdum(irf19pa)
	dfdy(if19,if19) = sum(b(1:2))

	! d(f19)/d(ne19)      
    b(1) = ratdum(irne19enu)
	dfdy(if19,ine19) = b(1)

!!! Hot CNO Neon-18 elements

	! d(ne18)/d(h1)
    b(1) = y(if17) * ratdum(irf17pg)
	dfdy(ine18,ih1) = b(1)

	! d(ne18)/d(he4)
    b(1) = -y(ine18) * ratdum(irne18ap)
	dfdy(ine18,ihe4) = b(1)

	! d(ne18)/d(f17)
    b(1) = y(ih1) * ratdum(irf17pg)
	dfdy(ine18,if17) = b(1)

	! d(ne18)/d(ne18)
    b(1) = -ratdum(irne18gp)
    b(2) = -ratdum(irne18enu)
    b(3) = -y(ihe4) * ratdum(irne18ap)
	dfdy(ine18,ine18) = esum3(b)

!!! Hot CNO Neon-19 elements

	! d(ne19)/d(h1)	
	b(1) = -y(ine19)*ratdum(irne19pg)
	dfdy(ine19,ih1) = b(1)

	! d(ne19)/d(he4)	
    b(1) = y(io15) * ratdum(iro15ag)
	dfdy(ine19,ihe4) = b(1)

	! d(ne19)/d(o15)	
    b(1) = y(ihe4) * ratdum(iro15ag)
	dfdy(ine19,io15) = b(1)

	! d(ne19)/d(ne19)	
    b(1) = -ratdum(irne19enu)
    b(2) = -y(ih1)*ratdum(irne19pg)
	dfdy(ine19,ine19) = sum(b(1:2))

!!! Hot CNO Neon-20 elements

	! d(ne20)/d(h1)	
	b(1) = -y(ine20)*ratdum(irne20pg)
	dfdy(ine20,ih1) = b(1)

	! d(ne20)/d(he4)	
    b(1) = y(io16) * ratdum(iroag)
	dfdy(ine20,ihe4) = b(1)
	
	! d(ne20)/d(o16)	
	b(1) = y(ihe4) * ratdum(iroag)
	dfdy(ine20,io16) = b(1)

	! d(ne20)/d(ne20)	
    b(1) = -y(ih1) * ratdum(irne20pg)
	dfdy(ine20,ine20) = b(1)

!!! Hot CNO Magnesium-22 elements

	! d(mg22)/d(h1)      
	b(1) = y(ine19) * ratdum(irne19pg)
    b(2) = +y(ine20) * ratdum(irne20pg)
    b(3) = -y(img22) * ratdum(iralam1) * (1.0d0 - ratdum(irdelta1))
	dfdy(img22,ih1) = esum3(b)

	! d(mg22)/d(he4)      
    b(1) = -y(img22)*ratdum(iralam1)*ratdum(irdelta1)
    b(2) = +y(ine18) * ratdum(irne18ap)
	dfdy(img22,ihe4) = sum(b(1:2))

	! d(mg22)/d(ne18)      
    b(1) = y(ihe4) * ratdum(irne18ap)
	dfdy(img22,ine18) = b(1)

	! d(mg22)/d(ne19)      
    b(1) = y(ih1)*ratdum(irne19pg)
	dfdy(img22,ine19) = b(1)

	! d(mg22)/d(ne20)      
    b(1) = y(ih1)*ratdum(irne20pg)
	dfdy(img22,ine20) = b(1)

	! d(mg22)/d(mg22)      
    b(1) = -y(ih1) * ratdum(iralam1) * (1.0d0 - ratdum(irdelta1))
    b(2) = -y(ihe4) * ratdum(iralam1) * ratdum(irdelta1)
	dfdy(img22,img22) = sum(b(1:2))

!!! Hot CNO Sulfur-30 elements

	! d(s30)/d(h1)
    b(1) = y(img22) * ratdum(iralam1) * (1.0d0 - ratdum(irdelta1))
    b(2) = -y(is30)*ratdum(iralam2) * (1.0d0 - ratdum(irdelta2))
	dfdy(is30,ih1) = sum(b(1:2))

	! d(s30)/d(he4)
    b(1) = y(img22) * ratdum(iralam1) * ratdum(irdelta1)
    b(2) = -y(is30) * ratdum(iralam2) * ratdum(irdelta2)
	dfdy(is30,ihe4) = sum(b(1:2))
 
	! d(s30)/d(mg22)
    b(1) = y(ih1)*ratdum(iralam1) * (1.0d0 - ratdum(irdelta1))
    b(2) = + y(ihe4)*ratdum(iralam1)*ratdum(irdelta1)
	dfdy(is30,img22) = sum(b(1:2))
	
	! d(s30)/d(s30)
    b(1) = -y(ih1) * ratdum(iralam2) * (1.0d0 - ratdum(irdelta2))
    b(2) = -y(ihe4) * ratdum(iralam2) * ratdum(irdelta2)
	dfdy(is30,is30) = sum(b(1:2))

!!! Hot CNO Nickel-56 elements

	! d(ni56)/d(h1)     
	b(1) = y(is30)*ratdum(iralam2) * (1.0d0 - ratdum(irdelta2))
	dfdy(ini56,ih1) = b(1)

	! d(ni56)/d(he4)     
    b(1) = y(is30)*ratdum(iralam2)*ratdum(irdelta2)
	dfdy(ini56,ihe4) = b(1)

	! d(ni56)/d(s30)     
    b(1) = y(ih1) * ratdum(iralam2) * (1.0d0 - ratdum(irdelta2))
    b(2) = +y(ihe4) * ratdum(iralam2) * ratdum(irdelta2)
	dfdy(ini56,is30) = sum(b(1:2))

  end subroutine dfdy_isotopes_pphotcno



  ! Computes the instantaneous energy generation rate

  subroutine ener_gener_rate(dydt, enuc)

    use actual_network, only: nspec, mion, enuc_conv2

    implicit none

    double precision :: dydt(nspec), enuc

    ! This is basically e = m c**2

    enuc = sum(dydt(:) * mion(:)) * enuc_conv2

  end subroutine ener_gener_rate



  subroutine set_up_screening_factors()
    ! Compute and store the more expensive screening factors

    use screening_module, only: add_screening_factor
    use network, only: aion, zion

    implicit none

    ! note: we need to set these up in the same order that we evaluate the
    ! rates in actual_rhs.f90 (yes, it's ugly)
    
	!TODO: Need to add the correct screening factors here
	!call add_screening_factor(zion(ihe4),aion(ihe4),zion(ihe4),aion(ihe4))
    !call add_screening_factor(zion(ihe4),aion(ihe4),4.0d0,8.0d0)
    !call add_screening_factor(zion(ic12),aion(ic12),zion(ihe4),aion(ihe4))
    !call add_screening_factor(zion(ic12),aion(ic12),zion(ic12),aion(ic12))
    !call add_screening_factor(zion(ic12),aion(ic12),zion(io16),aion(io16))
    !call add_screening_factor(zion(ic12),aion(ic12),zion(io16),aion(io16))
    !call add_screening_factor(zion(io16),aion(io16),zion(ihe4),aion(ihe4))
    !call add_screening_factor(zion(ine20),aion(ine20),zion(ihe4),aion(ihe4))
    !call add_screening_factor(zion(img24),aion(img24),zion(ihe4),aion(ihe4))
    !call add_screening_factor(20.0d0,40.0d0,zion(ihe4),aion(ihe4))

  end subroutine set_up_screening_factors


  subroutine update_unevolved_species(state)

    implicit none

    type (burn_t)    :: state

  end subroutine update_unevolved_species

end module actual_rhs_module
