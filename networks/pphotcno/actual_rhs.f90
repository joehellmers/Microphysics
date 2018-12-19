module actual_rhs_module

  use network
  use eos_type_module
  use burn_type_module
  use temperature_integration_module, only: temperature_rhs, temperature_jac
  use sneut_module, only: sneut5
  use actual_network, only: nrates
  use rate_type_module

  implicit none

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

	print *,"{JH} actual_rhs_init: Starting..."

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
    
    double precision :: scfac(nrates),  dscfacdt(nrates),  dscfacdd(nrates)

	print *,"{JH} get_rates: Starting..."

    ! Get the data from the state

    rho  = state % rho
    temp = state % T
    y    = state % xn * aion_inv

    ! Get the raw reaction rates
    if (use_tables) then
       call pphotcnotab(temp, rho, ratraw, dratrawdt, dratrawdd)
    else
       call pphotcnorat(temp, rho, ratraw, dratrawdt, dratrawdd, state%y_e)
    endif

    ! Do the screening here because the corrections depend on the composition
    call screen_pphotcno(temp, rho, y,             &
                     ratraw, dratrawdt, dratrawdd, &
                     ratdum, dratdumdt, dratdumdd, &
                     scfac,  dscfacdt,  dscfacdd)

    ! Save the rate data, for the Jacobian later if we need it.

    rr % rates(1,:) = ratdum
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

	print *,"{JH} pphotcnotab: Starting..."

    ! Set the density dependence array
	dtab(irpp)    = bden
	dtab(irdpg)   = bden
	dtab(ir33)    = bden
	dtab(irhe3ag) = bden
	dtab(irbeec)  = bden*ye
	dtab(irbepg)  = bden
	dtab(irb8gp)  = 1.0d0
	dtab(irb8ep)  = 1.0d0
	dtab(irli7pa) = bden
	dtab(irpep)   = ye*bden*bden
	dtab(irhep)   = bden


	dtab(irc12pg)  = bden
	dtab(irn13gp)  = 1.0d0
	dtab(irn13enu) = 1.0d0
	dtab(irc13pg)  = bden
	dtab(irn14gp)  = 1.0d0
	dtab(irn14pg)  = bden
	dtab(iro15gp)  = 1.0d0
	dtab(iro15enu) = 1.0d0
	dtab(irn15pa)  = bden
	dtab(irc12ap)  = bden
	dtab(irn15pg)  = bden
	dtab(iro16gp)  = 1.0d0
	dtab(iro16pg)  = bden
	dtab(irf17gp)  = 1.0d0
	dtab(irf17enu) = 1.0d0
	dtab(iro17pa)  = bden
	dtab(irn14ap)  = bden
	dtab(iro17pg)  = bden
	dtab(irf18gp)  = 1.0d0
	dtab(irf18enu) = 1.0d0
	dtab(iro18pa)  = bden
	dtab(irn15ap)  = bden
	dtab(iro18pg)  = bden
	dtab(irf19gp)  = 1.0d0
	dtab(irf19pa)  = bden
	dtab(iro16ap)  = bden
	dtab(irn13pg)  = bden
	dtab(iro14gp)  = 1.0d0
	dtab(iro14enu) = 1.0d0
	dtab(iro14ap)  = bden
	dtab(irf17pa)  = bden
	dtab(irf17pg)  = bden
	dtab(irne18gp) = 1.0d0
	dtab(irne18enu)= 1.0d0
	dtab(irf18pa)  = bden
	dtab(iro15ap)  = bden
	dtab(ir3a)     = bden*bden
	dtab(irg3a)    = 1.0d0
	dtab(ircag)    = bden
	dtab(iroga)    = 1.0d0
	dtab(iroag)    = bden
	dtab(irne18ap) = bden
	dtab(irne19pg) = bden
	dtab(iro15ag)  = bden
	dtab(irtiap)   = bden
	dtab(irsi26ap) = bden
	dtab(irne19enu)= 1.0d0
	dtab(irne20pg) = bden

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

    call set_pphotcnorat()

  end subroutine create_rates_table



  subroutine set_pphotcnorat()

    implicit none

    double precision :: btemp, bden, ratraw(nrates), dratrawdt(nrates), dratrawdd(nrates)
    integer :: i, j

    bden = 1.0d0

	print *,"{JH} set_pphotcnorat: Starting..."

    do i = 1, tab_imax

       btemp = tab_tlo + dble(i-1) * tab_tstp
       btemp = 10.0d0**(btemp)
	   ! TODO: is 1.0 the correct value of ye to pass in here?
       call pphotcnorat(btemp, bden, ratraw, dratrawdt, dratrawdd,1.0d0)

       ttab(i) = btemp

       do j = 1, nrates

          rattab(j,i)    = ratraw(j)
          drattabdt(j,i) = dratrawdt(j)
          drattabdd(j,i) = dratrawdd(j)

       enddo

    enddo

  end subroutine set_pphotcnorat



  subroutine actual_rhs(state)

    implicit none

    ! This routine sets up the system of ODE's for the pphotcno
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

	print *,"{JH} actual_rhs: Starting..."

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

	print *,"{JH} actual_jac: Starting..."

    call get_rates(state, rr)

    ! Get the data from the state

    rho  = state % rho
    temp = state % T
    abar = state % abar
    zbar = state % zbar
    y    = state % xn * aion_inv

    ! Species Jacobian elements with respect to other species

    call dfdy_isotopes_pphotcno(y, state % jac(1:nspec,1:nspec), rr % rates(1,:))

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

    call rhs(y, rr % rates(1,:), rr % rates(1,:), state % jac(1:nspec,net_itemp), deriva)

    call ener_gener_rate(state % jac(1:nspec,net_itemp), state % jac(net_ienuc,net_itemp))
    state % jac(net_ienuc,net_itemp) = state % jac(net_ienuc,net_itemp) - dsneutdt

    ! Temperature Jacobian elements

    call temperature_jac(state)

  end subroutine actual_jac



  ! Evaluates the right hand side of the pphotcno ODEs

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

	print *,"{JH} rhs: Starting..."

    dydt(1:nspec) = ZERO

!TODO: Do we need to use the a(15) array to sum up the various parts of the reactions?
! set up the system of ode's :
! h1 reactions
      dydt(ih1) =  -y(ih1)*y(ih1)*(rate(irpp) + rate(irpep)) &
                  - y(ih1)*y(ih2)*rate(irdpg) &
                  + y(ihe3)*y(ihe3)*rate(ir33) &
                  - y(ih1)*y(ili7)*rate(irli7pa) &
                  - y(ibe7)*y(ih1)*rate(irbepg) &
                  + y(ib8)*rate(irb8gp) &
                  - y(ihe3)*y(ih1)*rate(irhep)

! h2 reactions
      dydt(ih2) =  -y(ih1)*y(ih2)*rate(irdpg) &
                  + 0.5d0*y(ih1)*y(ih1)*(rate(irpp) + rate(irpep))

! he3 reactions
      dydt(ihe3) =  -y(ihe3)*y(ihe3)*rate(ir33) &
                   + y(ih1)*y(ih2)*rate(irdpg) &
                   - y(ihe3)*y(ihe4)*rate(irhe3ag) &
                   - y(ihe3)*y(ih1)*rate(irhep)

! 4he reactions
      dydt(ihe4) =  0.5d0*y(ihe3)*y(ihe3)*rate(ir33) &
                  - y(ihe3)*y(ihe4)*rate(irhe3ag) &
                  + 2.0d0*y(ili7)*y(ih1)*rate(irli7pa) &
                  + 2.0d0*y(ib8)*rate(irb8ep) &
                  + y(ihe3)*y(ih1)*rate(irhep)

! li7 reactions
      dydt(ili7) =  -y(ili7)*y(ih1)*rate(irli7pa) &
                  + y(ibe7)*rate(irbeec)

! be7 reactions
      dydt(ibe7) =  -y(ibe7)*y(ih1)*rate(irbepg) &
                  + y(ihe3)*y(ihe4)*rate(irhe3ag) &
                  - y(ibe7)*rate(irbeec) &
                  + y(ib8)*rate(irb8gp)

! b8 reactions
      dydt(ib8)  =  y(ibe7)*y(ih1)*rate(irbepg) &
                  - y(ib8)*rate(irb8ep) &
                  - y(ib8)*rate(irb8gp)




! add in the hot cno contributions
! h1 reactions
      dydt(ih1)  = dydt(ih1) &
                   - y(ic12)*y(ih1)*rate(irc12pg) &
                   + y(in13)*rate(irn13gp) &
                   - y(ic13)*y(ih1)*rate(irc13pg) &
                   + y(in14)*rate(irn14gp) &
                   - y(in14)*y(ih1)*rate(irn14pg) &
                   + y(io15)*rate(iro15gp) &
                   - y(in15)*y(ih1)*rate(irn15pa) &
                   + y(ic12)*y(ihe4)*rate(irc12ap) &
                   - y(in15)*y(ih1)*rate(irn15pg) &
                   + y(io16)*rate(iro16gp) &
                   - y(io16)*y(ih1)*rate(iro16pg) &
                   + y(if17)*rate(irf17gp)

      dydt(ih1)  = dydt(ih1) &
                   - y(io17)*y(ih1)*rate(iro17pa) &
                   + y(in14)*y(ihe4)*rate(irn14ap) &
                   - y(io17)*y(ih1)*rate(iro17pg) &
                   + y(if18)*rate(irf18gp) &
                   - y(io18)*y(ih1)*rate(iro18pa) &
                   + y(in15)*y(ihe4)*rate(irn15ap) &
                   - y(io18)*y(ih1)*rate(iro18pg) &
                   + y(if19)*rate(irf19gp) &
                   - y(if19)*y(ih1)*rate(irf19pa) &
                   + y(io16)*y(ihe4)*rate(iro16ap) &
                   - y(in13)*y(ih1)*rate(irn13pg) &
                   + y(io14)*rate(iro14gp)

      dydt(ih1)  = dydt(ih1) &
                   + y(io14) * y(ihe4) * rate(iro14ap) &
                   - y(if17) * y(ih1) * rate(irf17pa) &
                   - y(if17) * y(ih1) * rate(irf17pg) &
                   + y(ine18) * rate(irne18gp) &
                   - y(if18) * y(ih1) * rate(irf18pa) &
                   + y(io15) * y(ihe4) * rate(iro15ap) &
                   - 8.0d0 * y(img22) * y(ih1) * &
                     rate(iralam1)*(1.0d0 - ratdum(irdelta1)) &
                   - 26.0d0 * y(is30) * y(ih1) * &
                     rate(iralam2)*(1.0d0 - ratdum(irdelta2)) &
                   - 3.0d0*y(ine19)*y(ih1)*rate(irne19pg) &
                   - 2.0d0*y(ine20)*y(ih1)*rate(irne20pg)


! he4 reactions
      dydt(ihe4) = dydt(ihe4) &
                   + y(in15)*y(ih1)*rate(irn15pa) &
                   - y(ic12)*y(ihe4)*rate(irc12ap) &
                   + y(io17)*y(ih1)*rate(iro17pa) &
                   - y(in14)*y(ihe4)*rate(irn14ap) &
                   + y(io18)*y(ih1)*rate(iro18pa) &
                   - y(in15)*y(ihe4)*rate(irn15ap) &
                   + y(if19)*y(ih1)*rate(irf19pa) &
                   - y(io16)*y(ihe4)*rate(iro16ap) &
                   - y(io14) * y(ihe4) * rate(iro14ap) &
                   + y(if17) * y(ih1) * rate(irf17pa) &
                   + y(if18) * y(ih1) * rate(irf18pa)

      dydt(ihe4) = dydt(ihe4) &
                   - y(io15) * y(ihe4) * rate(iro15ap) &
                   - 0.5d0 * y(ihe4) * y(ihe4) * y(ihe4) * rate(ir3a) &
                   + 3.0d0 * y(ic12) * rate(irg3a) &
                  -  y(io16) * y(ihe4) * rate(iroag) &
                   - 2.0d0 * y(img22)* y(ihe4) * &
                     rate(iralam1) * ratdum(irdelta1) &
                   - 6.5d0 * y(is30) * y(ihe4) * &
                     rate(iralam2) * ratdum(irdelta2) &
                   + y(io16) * rate(iroga) &
                   - y(ic12) * y(ihe4) * rate(ircag) &
                   - y(ine18) * y(ihe4) * rate(irne18ap) &
                   - y(io15) * y(ihe4) * rate(iro15ag)


! c12 reactions
      dydt(ic12) =  -y(ic12)*y(ih1)*rate(irc12pg) &
                   + y(in13)*rate(irn13gp) &
                   + y(in15)*y(ih1)*rate(irn15pa) &
                   - y(ic12)*y(ihe4)*rate(irc12ap) &
                   + sixth * y(ihe4) * y(ihe4) * y(ihe4) * rate(ir3a) &
                   - y(ic12) * rate(irg3a) &
                   + y(io16) * rate(iroga) &
                   - y(ic12) * y(ihe4) * rate(ircag)

! c13 reactions
      dydt(ic13) =   y(in13)*rate(irn13enu) &
                   - y(ic13)*y(ih1)*rate(irc13pg) &
                   + y(in14)*rate(irn14gp)


! n13 reactions
      dydt(in13) =   y(ic12)*y(ih1)*rate(irc12pg) &
                   - y(in13)*rate(irn13gp) &
                   - y(in13)*rate(irn13enu) &
                   - y(in13)*y(ih1)*rate(irn13pg) &
                   + y(io14)*rate(iro14gp)


! n14 reactions
      dydt(in14) =   y(ic13)*y(ih1)*rate(irc13pg) &
                   - y(in14)*rate(irn14gp) &
                   - y(in14)*y(ih1)*rate(irn14pg) &
                   + y(io15)*rate(iro15gp) &
                   + y(io17)*y(ih1)*rate(iro17pa) &
                   - y(in14)*y(ihe4)*rate(irn14ap) &
                   + y(io14)*rate(iro14enu)


! n15 reactions
      dydt(in15) =   y(io15)*rate(iro15enu) &
                   - y(in15)*y(ih1)*rate(irn15pa) &
                   + y(ic12)*y(ihe4)*rate(irc12ap) &
                   - y(in15)*y(ih1)*rate(irn15pg) &
                   + y(io16)*rate(iro16gp) &
                   + y(io18)*y(ih1)*rate(iro18pa) &
                   - y(in15)*y(ihe4)*rate(irn15ap)


! o14 reactions
      dydt(io14) =   y(in13)*y(ih1)*rate(irn13pg) &
                   - y(io14)*rate(iro14gp) &
                   - y(io14)*rate(iro14enu) &
                   - y(io14) * y(ihe4) * rate(iro14ap) &
                   + y(if17) * y(ih1) * rate(irf17pa)


! o15 reactions
      dydt(io15) =   y(in14)*y(ih1)*rate(irn14pg) &
                   - y(io15)*rate(iro15gp) &
                   - y(io15)*rate(iro15enu) &
                   + y(if18) * y(ih1) * rate(irf18pa) &
                   - y(io15) * y(ihe4) * rate(iro15ap) &
                   - y(io15) * y(ihe4) * rate(iro15ag)


! o16 reactions
      dydt(io16) =   y(in15)*y(ih1)*rate(irn15pg) &
                   - y(io16)*rate(iro16gp) &
                   - y(io16)*y(ih1)*rate(iro16pg) &
                   + y(if17)*rate(irf17gp) &
                   + y(if19)*y(ih1)*rate(irf19pa) &
                   - y(io16)*y(ihe4)*rate(iro16ap) &
                   - y(io16) * y(ihe4) * rate(iroag) &
                   - y(io16) * rate(iroga) &
                   + y(ic12) * y(ihe4) * rate(ircag)


! o17 reactions
      dydt(io17) =   y(if17)*rate(irf17enu) &
                   - y(io17)*y(ih1)*rate(iro17pa) &
                   + y(in14)*y(ihe4)*rate(irn14ap) &
                   - y(io17)*y(ih1)*rate(iro17pg) &
                   + y(if18)*rate(irf18gp)


! o18 reactions
      dydt(io18) =   y(if18)*rate(irf18enu) &
                   - y(io18)*y(ih1)*rate(iro18pa) &
                   + y(in15)*y(ihe4)*rate(irn15ap) &
                   - y(io18)*y(ih1)*rate(iro18pg) &
                   + y(if19)*rate(irf19gp)


! f17 reactions
      dydt(if17) =   y(io16)*y(ih1)*rate(iro16pg) &
                   - y(if17)*rate(irf17gp) &
                   - y(if17)*rate(irf17enu) &
                   + y(io14) * y(ihe4) * rate(iro14ap) &
                   - y(if17) * y(ih1) * rate(irf17pa) &
                   - y(if17) * y(ih1) * rate(irf17pg) &
                   + y(ine18) * rate(irne18gp)


! f18 reactions
      dydt(if18) =   y(io17)*y(ih1)*rate(iro17pg) &
                   - y(if18)*rate(irf18gp) &
                   - y(if18)*rate(irf18enu) &
                   + y(ine18) * rate(irne18enu) &
                   - y(if18) * y(ih1) * rate(irf18pa) &
                   + y(io15) * y(ihe4) * rate(iro15ap)


! f19 reactions
      dydt(if19) =   y(io18)*y(ih1)*rate(iro18pg) &
                   - y(if19)*rate(irf19gp) &
                   - y(if19)*y(ih1)*rate(irf19pa) &
                   + y(io16)*y(ihe4)*rate(iro16ap) &
                   + y(ine19) * rate(irne19enu)


! ne18 reactions
      dydt(ine18) =  y(if17) * y(ih1) * rate(irf17pg) &
                   - y(ine18) * rate(irne18gp) &
                   - y(ine18) * rate(irne18enu) &
                   - y(ine18) * y(ihe4) * rate(irne18ap)


! ne19 reactions
      dydt(ine19) = y(io15) * y(ihe4) * rate(iro15ag) &
                   - y(ine19) * rate(irne19enu) &
                   - y(ine19)*y(ih1)*rate(irne19pg)


! ne20 reactions
      dydt(ine20) = y(io16) * y(ihe4) * rate(iroag) &
                   - y(ine20)*y(ih1)*rate(irne20pg)



! mg22 reactions
! for rp o16(a,g)ne20(p,g)na21(p,g)mg22
!        f17(p,g)ne18(a,g)mg22
!        o15(a,g)ne19(p,g)na20(p,g)mg21(e-nu)na21(p,g)mg22

      dydt(img22) = -y(img22)*y(ih1)* &
                     rate(iralam1)*(1.0d0 - ratdum(irdelta1)) &
                   - y(img22)*y(ihe4)*rate(iralam1)*ratdum(irdelta1) &
                   + y(ine18) * y(ihe4) * rate(irne18ap) &
                   + y(ine19)*y(ih1)*rate(irne19pg) &
                   + y(ine20)*y(ih1)*rate(irne20pg)


! s30 reactions
      dydt(is30) =   y(img22)*y(ih1)* &
                     rate(iralam1)*(1.0d0 - ratdum(irdelta1)) &
                   + y(img22)*y(ihe4)*rate(iralam1)*ratdum(irdelta1) &
                   - y(is30)*y(ih1)* &
                     rate(iralam2)*(1.0d0 - ratdum(irdelta2)) &
                   - y(is30)*y(ihe4)*rate(iralam2)*ratdum(irdelta2)


! ni56 reactions
      dydt(ini56) = y(is30)*y(ih1)* &
                    rate(iralam2)*(1.0d0 - ratdum(irdelta2)) &
                  + y(is30)*y(ihe4)*rate(iralam2)*ratdum(irdelta2)

  end subroutine rhs



  subroutine pphotcnorat(btemp, bden, ratraw, dratrawdt, dratrawdd, ye)

    ! this routine generates unscreened
    ! nuclear reaction rates for the pphotcno network.

    use tfactors_module
    use aprox_rates_module
    use amrex_constants_module, only: ZERO
    use extern_probin_module, only: use_c12ag_deboer17

    double precision :: btemp, bden, ye
    double precision :: ratraw(nrates), dratrawdt(nrates), dratrawdd(nrates)

    integer          :: i
    double precision :: rrate,drratedt,drratedd
	
    double precision :: ff1,dff1dt,dff1dd,ff2,dff2dt,dff2dd,tot,dtotdt,dtotdd,invtot
    type (tf_t)      :: tf

	print *,"{JH} pphotcnorat: Starting...btemp=",btemp

    do i=1,nrates
       ratraw(i)    = ZERO
       dratrawdt(i) = ZERO
       dratrawdd(i) = ZERO
    enddo

    if (btemp .lt. 1.0d6) return

    ! get the temperature factors
    call get_tfactors(btemp, tf)

!!!! START NEW STUFF HERE !!!!

! p(p,e+nu)d
      call rate_pp(tf,bden, &
                   ratraw(irpp),dratrawdt(irpp),dratrawdd(irpp), &
                   rrate,drratedt,drratedd)

! p(e-p,nu)d
      call rate_pep(tf,bden,ye, &
                   ratraw(irpep),dratrawdt(irpep),dratrawdd(irpep), &
                   rrate,drratedt,drratedd)

! d(p,g)3he
      call rate_dpg(tf,bden, &
                   ratraw(irdpg),dratrawdt(irdpg),dratrawdd(irdpg), &
                   rrate,drratedt,drratedd)

! he3(p,e+nu)he4
      call rate_hep(tf,bden, &
                   ratraw(irhep),dratrawdt(irhep),dratrawdd(irhep), &
                   rrate,drratedt,drratedd)

! he3(he3,2p)he4
      call rate_he3he3(tf,bden, &
                   ratraw(ir33),dratrawdt(ir33),dratrawdd(ir33), &
                   rrate,drratedt,drratedd)

! 3he(4he,nu)7be
      call rate_he3he4(tf,bden, &
               ratraw(irhe3ag),dratrawdt(irhe3ag),dratrawdd(irhe3ag), &
               rrate,drratedt,drratedd)

! 7be(e-,nu)7li
      call rate_be7em(tf,bden,ye, &
                   ratraw(irbeec),dratrawdt(irbeec),dratrawdd(irbeec), &
                   rrate,drratedt,drratedd)

! 7be(p,g)8b
      call rate_be7pg(tf,bden, &
                   ratraw(irbepg),dratrawdt(irbepg),dratrawdd(irbepg), &
                   ratraw(irb8gp),dratrawdt(irb8gp),dratrawdd(irb8gp))

! 8b(p=>n)8be=>2 he4  positron decay (half-life = 0.77 sec)
      call rate_b8ep(tf,bden, &
                   ratraw(irb8ep),dratrawdt(irb8ep),dratrawdd(irb8ep), &
                   rrate,drratedt,drratedd)

! 7li(p,g)8be => 2a   and 7li(p,a)a
      call rate_li7pag(tf,bden, &
                 ratraw(irli7pa),dratrawdt(irli7pa),dratrawdd(irli7pa), &
                 rrate,drratedt,drratedd)

! from hotcno
      call rate_c12pg(tf,bden, &
           ratraw(irc12pg),dratrawdt(irc12pg),dratrawdd(irc12pg), &
           ratraw(irn13gp),dratrawdt(irn13gp),dratrawdd(irn13gp))

! 13n(e+nu)13c
      call rate_n13em(tf,bden, &
           ratraw(irn13enu),dratrawdt(irn13enu),dratrawdd(irn13enu), &
           rrate,drratedt,drratedd)

! 13c(p,g)14n
      call rate_c13pg(tf,bden, &
           ratraw(irc13pg),dratrawdt(irc13pg),dratrawdd(irc13pg), &
           ratraw(irn14gp),dratrawdt(irn14gp),dratrawdd(irn14gp))

! 14n(p,g)15o
      call rate_n14pg(tf,bden, &
           ratraw(irn14pg),dratrawdt(irn14pg),dratrawdd(irn14pg), &
           ratraw(iro15gp),dratrawdt(iro15gp),dratrawdd(iro15gp))

! 15o(e+nu)15n
      call rate_o15em(tf,bden, &
           ratraw(iro15enu),dratrawdt(iro15enu),dratrawdd(iro15enu), &
           rrate,drratedt,drratedd)

! 15n(p,a)12c 
      call rate_n15pa(tf,bden, &
           ratraw(irn15pa),dratrawdt(irn15pa),dratrawdd(irn15pa), &
           ratraw(irc12ap),dratrawdt(irc12ap),dratrawdd(irc12ap))

! 15n(p,g)16o
      call rate_n15pg(tf,bden, &
           ratraw(irn15pg),dratrawdt(irn15pg),dratrawdd(irn15pg), &
           ratraw(iro16gp),dratrawdt(iro16gp),dratrawdd(iro16gp))

! 16o(p,g)17f
      call rate_o16pg(tf,bden, &
           ratraw(iro16pg),dratrawdt(iro16pg),dratrawdd(iro16pg), &
           ratraw(irf17gp),dratrawdt(irf17gp),dratrawdd(irf17gp))

! 17f(e+nu)17o
      call rate_f17em(tf,bden, &
           ratraw(irf17enu),dratrawdt(irf17enu),dratrawdd(irf17enu), &
           rrate,drratedt,drratedd)

! 17o(p,a)14n
      call rate_o17pa(tf,bden, &
           ratraw(iro17pa),dratrawdt(iro17pa),dratrawdd(iro17pa), &
           ratraw(irn14ap),dratrawdt(irn14ap),dratrawdd(irn14ap))

! 17o(p,g)18f
      call rate_o17pg(tf,bden, &
           ratraw(iro17pg),dratrawdt(iro17pg),dratrawdd(iro17pg), &
           ratraw(irf18gp),dratrawdt(irf18gp),dratrawdd(irf18gp))

! 18f(e+nu)18o
      call rate_f18em(tf,bden, &
           ratraw(irf18enu),dratrawdt(irf18enu),dratrawdd(irf18enu), &
           rrate,drratedt,drratedd)

! 18o(p,a)15n
      call rate_o18pa(tf,bden, &
           ratraw(iro18pa),dratrawdt(iro18pa),dratrawdd(iro18pa), &
           ratraw(irn15ap),dratrawdt(irn15ap),dratrawdd(irn15ap))

! 18o(p,g)19f
      call rate_o18pg(tf,bden, &
           ratraw(iro18pg),dratrawdt(iro18pg),dratrawdd(iro18pg), &
           ratraw(irf19gp),dratrawdt(irf19gp),dratrawdd(irf19gp))

! 19f(p,a)16o
      call rate_f19pa(tf,bden, &
           ratraw(irf19pa),dratrawdt(irf19pa),dratrawdd(irf19pa), &
           ratraw(iro16ap),dratrawdt(iro16ap),dratrawdd(iro16ap))



! add these for the hot cno cycles
! 13n(p,g)14o
      call rate_n13pg(tf,bden, &
           ratraw(irn13pg),dratrawdt(irn13pg),dratrawdd(irn13pg), &
           ratraw(iro14gp),dratrawdt(iro14gp),dratrawdd(iro14gp))

! 14o(e+nu)14n
      call rate_o14em(tf,bden, &
           ratraw(iro14enu),dratrawdt(iro14enu),dratrawdd(iro14enu), &
           rrate,drratedt,drratedd)

! 14o(a,p)17f cf88 q = 1.191
      call rate_o14ap(tf,bden, &
           ratraw(iro14ap),dratrawdt(iro14ap),dratrawdd(iro14ap), &
           ratraw(irf17pa),dratrawdt(irf17pa),dratrawdd(irf17pa))

! 17f(p,g)18ne
      call rate_f17pg(tf,bden, &
           ratraw(irf17pg),dratrawdt(irf17pg),dratrawdd(irf17pg), &
           ratraw(irne18gp),dratrawdt(irne18gp),dratrawdd(irne18gp))

! 18ne(e+nu)18f
      call rate_ne18em(tf,bden, &
           ratraw(irne18enu),dratrawdt(irne18enu),dratrawdd(irne18enu), &
           rrate,drratedt,drratedd)

! 18f(p,a)15o
      call rate_f18pa(tf,bden, &
           ratraw(irf18pa),dratrawdt(irf18pa),dratrawdd(irf18pa), &
           ratraw(iro15ap),dratrawdt(iro15ap),dratrawdd(iro15ap))

! triple alpha to c12
      call rate_tripalf(tf,bden, &
           ratraw(ir3a),dratrawdt(ir3a),dratrawdd(ir3a), &
           ratraw(irg3a),dratrawdt(irg3a),dratrawdd(irg3a))

! c12(a,g)o16
      call rate_c12ag(tf,bden, &
           ratraw(ircag),dratrawdt(ircag),dratrawdd(ircag), &
           ratraw(iroga),dratrawdt(iroga),dratrawdd(iroga))

! 16o(a,g)20ne
      call rate_o16ag(tf,bden, &
           ratraw(iroag),dratrawdt(iroag),dratrawdd(iroag), &
           rrate,drratedt,drratedd)

! 18ne(a,p)21na
      call rate_ne18ap(tf,bden, &
           ratraw(irne18ap),dratrawdt(irne18ap),dratrawdd(irne18ap), &
           rrate,drratedt,drratedd)

! 19ne(p,g)20na
      call rate_ne19pg(tf,bden, &
           ratraw(irne19pg),dratrawdt(irne19pg),dratrawdd(irne19pg), &
           rrate,drratedt,drratedd)

! 15o(a,g)19ne
      call rate_o15ag(tf,bden, &
           ratraw(iro15ag),dratrawdt(iro15ag),dratrawdd(iro15ag), &
           rrate,drratedt,drratedd)

! 44ti(a,p)47v
      call rate_ti44ap(tf,bden, &
           ratraw(irtiap),dratrawdt(irtiap),dratrawdd(irtiap), &
           rrate,drratedt,drratedd)

!      g47v=1.0+exp(-1.206*t9m1+1.059+9.997d-2*t9)
!      ratraw(irtiap) = ratraw(irtiap)*g47v

! 26si(a,p)29p
      call rate_si26ap(tf,bden, &
           ratraw(irsi26ap),dratrawdt(irsi26ap),dratrawdd(irsi26ap), &
           rrate,drratedt,drratedd)

! 19ne(e+nu)19f
      call rate_ne19em(tf,bden, &
           ratraw(irne19enu),dratrawdt(irne19enu),dratrawdd(irne19enu), &
           rrate,drratedt,drratedd)

! 20ne(p,g)21na
      call rate_ne20pg(tf,bden, &
           ratraw(irne20pg),dratrawdt(irne20pg),dratrawdd(irne20pg), &
           rrate,drratedt,drratedd)

  end subroutine pphotcnorat



  subroutine screen_pphotcno(btemp, bden, y, &
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

	print *,"{JH} screen_pphotcno: Starting..."	
    
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
	!TODO: Do we need to do this first?    

! ppp
	jscr = 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)

    jscr = jscr + 1
    call screen5(pstate,jscr,sc2a,sc2adt,sc2add)

    sc3a   = sc1a * sc2a
    sc3adt = sc1adt*sc2a + sc1a*sc2adt

    ratdum(ir3a)    = ratraw(ir3a) * sc3a
    dratdumdt(ir3a) = dratrawdt(ir3a)*sc3a + ratraw(ir3a)*sc3adt

    scfac(ir3a)     = sc3a
    dscfacdt(ir3a)  = sc3adt

! pp
	jscr = jscr + 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(irpp)    = ratraw(irpp) * sc1a
    dratdumdt(irpp) = dratrawdt(irpp)*sc1a+ratraw(irpp)*sc1adt
    scfac(irpp)     = sc1a
    dscfacdt(irpp)  = sc1adt

	!TODO: What is this for?    
	ratdum(irpep)    = ratraw(irpep) * sc1a
    dratdumdt(irpep) = dratrawdt(irpep)*sc1a + ratraw(irpep)*sc1adt
    scfac(irpep)     = sc1a
    dscfacdt(irpep)  = sc1adt

! d + p
    jscr = jscr + 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(irdpg)   = ratraw(irdpg) * sc1a
    dratdumdt(irdpg)= dratrawdt(irdpg)*sc1a+ratraw(irdpg)*sc1adt
    scfac(irdpg)    = sc1a
    dscfacdt(irdpg) = sc1adt

! he3 + p
	jscr = jscr + 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(irhep)    = ratraw(irhep) * sc1a
	dratdumdt(irhep) = dratrawdt(irhep)*sc1a + ratraw(irhep)*sc1adt
	scfac(irhep)     = sc1a
	dscfacdt(irhep)  = sc1adt

! he3 + he3
	jscr = jscr + 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(ir33)    = ratraw(ir33) * sc1a
	dratdumdt(ir33) = dratrawdt(ir33)*sc1a+ratraw(ir33)*sc1adt
	scfac(ir33)     = sc1a
	dscfacdt(ir33)  = sc1adt

! he3 + he4
	jscr = jscr + 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(irhe3ag)    = ratraw(irhe3ag) * sc1a
	dratdumdt(irhe3ag) = dratrawdt(irhe3ag)*sc1a &
		               + ratraw(irhe3ag)*sc1adt
	scfac(irhe3ag)     = sc1a
	dscfacdt(irhe3ag)  = sc1adt

! be7 + p and inverse
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(irbepg)   = ratraw(irbepg) * sc1a
	dratdumdt(irbepg)= dratrawdt(irbepg)*sc1a+ratraw(irbepg)*sc1adt
	scfac(irbepg)    = sc1a
	dscfacdt(irbepg) = sc1adt

! li7 + p
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(irli7pa)  = ratraw(irli7pa) * sc1a
	dratdumdt(irli7pa)= dratrawdt(irli7pa)*sc1a+ratraw(irli7pa)*sc1adt
	scfac(irli7pa)   = sc1a
	dscfacdt(irli7pa)= sc1adt


! hot cno contributions
! c12 + p
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(irc12pg)    = ratraw(irc12pg) * sc1a
	dratdumdt(irc12pg) =dratrawdt(irc12pg)*sc1a+ratraw(irc12pg)*sc1adt
	scfac(irc12pg)     = sc1a
	dscfacdt(irc12pg)  = sc1adt

! c13 + p
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(irc13pg)    = ratraw(irc13pg) * sc1a
	dratdumdt(irc13pg) =dratrawdt(irc13pg)*sc1a+ratraw(irc13pg)*sc1adt
	scfac(irc13pg)     = sc1a
	dscfacdt(irc13pg)  = sc1adt

! n14 + p
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(irn14pg)    = ratraw(irn14pg) * sc1a
	dratdumdt(irn14pg) =dratrawdt(irn14pg)*sc1a+ratraw(irn14pg)*sc1adt
	scfac(irn14pg)     = sc1a
	dscfacdt(irn14pg)  = sc1adt

! n15 + p
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(irn15pg)    = ratraw(irn15pg) * sc1a
	dratdumdt(irn15pg) =dratrawdt(irn15pg)*sc1a+ratraw(irn15pg)*sc1adt
	scfac(irn15pg)     = sc1a
	dscfacdt(irn15pg)  = sc1adt

	ratdum(irn15pa)    = ratraw(irn15pa) * sc1a
	dratdumdt(irn15pa) =dratrawdt(irn15pa)*sc1a+ratraw(irn15pa)*sc1adt
	scfac(irn15pa)     = sc1a
	dscfacdt(irn15pa)  = sc1adt


! c12 + a
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(irc12ap)    = ratraw(irc12ap) * sc1a
	dratdumdt(irc12ap) =dratrawdt(irc12ap)*sc1a+ratraw(irc12ap)*sc1adt
	scfac(irc12ap)     = sc1a
	dscfacdt(irc12ap)  = sc1adt

	ratdum(ircag)     = ratraw(ircag) * sc1a
	dratdumdt(ircag)  = dratrawdt(ircag)*sc1a + ratraw(ircag)*sc1adt
	scfac(ircag)      = sc1a
	dscfacdt(ircag)   = sc1adt

! o16 + p
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(iro16pg)    = ratraw(iro16pg) * sc1a
	dratdumdt(iro16pg) =dratrawdt(iro16pg)*sc1a+ratraw(iro16pg)*sc1adt
	scfac(iro16pg)     = sc1a
	dscfacdt(iro16pg)  = sc1adt

! o17 + p
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(iro17pg)    = ratraw(iro17pg) * sc1a
	dratdumdt(iro17pg) =dratrawdt(iro17pg)*sc1a+ratraw(iro17pg)*sc1adt
	scfac(iro17pg)     = sc1a
	dscfacdt(iro17pg)  = sc1adt

	ratdum(iro17pa)    = ratraw(iro17pa) * sc1a
	dratdumdt(iro17pa) =dratrawdt(iro17pa)*sc1a+ratraw(iro17pa)*sc1adt
	scfac(iro17pa)     = sc1a
	dscfacdt(iro17pa)  = sc1adt

! n14 + a
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(irn14ap)    = ratraw(irn14ap) * sc1a
	dratdumdt(irn14ap) =dratrawdt(irn14ap)*sc1a+ratraw(irn14ap)*sc1adt
	scfac(irn14ap)     = sc1a
	dscfacdt(irn14ap)  = sc1adt
	dscfacdd(irn14ap)  = sc1add

! o18 + p
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(iro18pg)    = ratraw(iro18pg) * sc1a
	dratdumdt(iro18pg) =dratrawdt(iro18pg)*sc1a+ratraw(iro18pg)*sc1adt
	scfac(iro18pg)     = sc1a
	dscfacdt(iro18pg)  = sc1adt

	ratdum(iro18pa)    = ratraw(iro18pa) * sc1a
	dratdumdt(iro18pa) =dratrawdt(iro18pa)*sc1a+ratraw(iro18pa)*sc1adt
	scfac(iro18pa)     = sc1a
	dscfacdt(iro18pa)  = sc1adt

! n15 + a
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(irn15ap)    = ratraw(irn15ap) * sc1a
	dratdumdt(irn15ap) =dratrawdt(irn15ap)*sc1a+ratraw(irn15ap)*sc1adt
	scfac(irn15ap)     = sc1a
	dscfacdt(irn15ap)  = sc1adt

! f19 + p
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(irf19pa)    = ratraw(irf19pa) * sc1a
	dratdumdt(irf19pa) =dratrawdt(irf19pa)*sc1a+ratraw(irf19pa)*sc1adt
	scfac(irf19pa)     = sc1a
	dscfacdt(irf19pa)  = sc1adt

! o16 + a
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(iro16ap)    = ratraw(iro16ap) * sc1a
	dratdumdt(iro16ap) =dratrawdt(iro16ap)*sc1a+ratraw(iro16ap)*sc1adt
	scfac(iro16ap)     = sc1a
	dscfacdt(iro16ap)  = sc1adt

! n13 + p
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(irn13pg)    = ratraw(irn13pg) * sc1a
	dratdumdt(irn13pg) =dratrawdt(irn13pg)*sc1a+ratraw(irn13pg)*sc1adt
	scfac(irn13pg)     = sc1a
	dscfacdt(irn13pg)  = sc1adt

! f17 + p
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(irf17pg)    = ratraw(irf17pg) * sc1a
	dratdumdt(irf17pg) =dratrawdt(irf17pg)*sc1a+ratraw(irf17pg)*sc1adt
	scfac(irf17pg)     = sc1a
	dscfacdt(irf17pg)  = sc1adt

	ratdum(irf17pa)    = ratraw(irf17pa) * sc1a
	dratdumdt(irf17pa) =dratrawdt(irf17pa)*sc1a+ratraw(irf17pa)*sc1adt
	scfac(irf17pa)     = sc1a
	dscfacdt(irf17pa)  = sc1adt

! o14 + a
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(iro14ap)    = ratraw(iro14ap) * sc1a
	dratdumdt(iro14ap) =dratrawdt(iro14ap)*sc1a+ratraw(iro14ap)*sc1adt
	scfac(iro14ap)     = sc1a
	dscfacdt(iro14ap)  = sc1adt

! f18 + p
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(irf18pa)    = ratraw(irf18pa) * sc1a
	dratdumdt(irf18pa) =dratrawdt(irf18pa)*sc1a+ratraw(irf18pa)*sc1adt
	scfac(irf18pa)     = sc1a
	dscfacdt(irf18pa)  = sc1adt

! o15 + a
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(iro15ap)    = ratraw(iro15ap) * sc1a
	dratdumdt(iro15ap) =dratrawdt(iro15ap)*sc1a+ratraw(iro15ap)*sc1adt
	scfac(iro15ap)     = sc1a
	dscfacdt(iro15ap)  = sc1adt

! o16 to ne20
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(iroag)    = ratraw(iroag) * sc1a
	dratdumdt(iroag) = dratrawdt(iroag)*sc1a + ratraw(iroag)*sc1adt
	scfac(iroag)    = sc1a
	dscfacdt(iroag) = sc1adt

! ne18 to na21
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(irne18ap)   = ratraw(irne18ap) * sc1a
	dratdumdt(irne18ap)= dratrawdt(irne18ap)*sc1a &
		               + ratraw(irne18ap)*sc1adt
	scfac(irne18ap)    = sc1a
	dscfacdt(irne18ap) = sc1adt

! o15 to ne19
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(iro15ag)    = ratraw(iro15ag) * sc1a
	dratdumdt(iro15ag) =dratrawdt(iro15ag)*sc1a+ratraw(iro15ag)*sc1adt
	scfac(iro15ag)     = sc1a
	dscfacdt(iro15ag)  = sc1adt

! ne19 + p
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(irne19pg)   = ratraw(irne19pg) * sc1a
	dratdumdt(irne19pg)= dratrawdt(irne19pg)*sc1a &
		               + ratraw(irne19pg)*sc1adt
	scfac(irne19pg)    = sc1a
	dscfacdt(irne19pg) = sc1adt

! si26 + a
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(irsi26ap)   = ratraw(irsi26ap) * sc1a
	dratdumdt(irsi26ap)= dratrawdt(irsi26ap)*sc1a &
		               + ratraw(irsi26ap)*sc1adt
	scfac(irsi26ap)    = sc1a
	dscfacdt(irsi26ap) = sc1adt
	dscfacdd(irsi26ap) = sc1add

! ti44 + alpha
	jscr = jscr + 1
	call screen5(pstate,jscr,sc1a,sc1adt,sc1add)
	ratdum(irtiap)    = ratraw(irtiap) * sc1a
	dratdumdt(irtiap) = dratrawdt(irtiap)*sc1a + ratraw(irtiap)*sc1adt
	scfac(irtiap)     = sc1a
	dscfacdt(irtiap)  = sc1adt

	print *,"{JH} screen_pphotcno: at end dratdumdt=",dratdumdt

  end subroutine screen_pphotcno


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

	print *,"{JH} dfdy_isotopes_pphotcno: Starting..."
	
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
    
	! TODO: Verify there are screening factors for non-ion related reactions.

	! p(p,e+nu)d
	! p(e-p,nu)d
	! d(p,g)3he
	! he3(p,e+nu)he4
	! he3(he3,2p)he4
	call add_screening_factor(zion(ihe3),aion(ihe3),zion(ihe3),aion(ihe3))
	! 3he(4he,nu)7be
	call add_screening_factor(zion(ihe3),aion(ihe3),zion(ihe4),aion(ihe4))
	! 7be(e-,nu)7li
	! 7be(p,g)8b
	! 8b(p=>n)8be=>2 he4  positron decay (half-life = 0.77 sec)
	! 7li(p,g)8be => 2a   and 7li(p,a)a
	! c12 + p
	! 13n(e+nu)13c
	! 13c(p,g)14n
	! 14n(p,g)15o
	! 15o(e+nu)15n
	! 15n(p,a)12c 
	! 15n(p,g)16o
	! 16o(p,g)17f
	! 17f(e+nu)17o
	! 17o(p,a)14n
	! 17o(p,g)18f
	! 18f(e+nu)18o
	! 18o(p,a)15n
	! 18o(p,g)19f
	! 19f(p,a)16o
	! 13n(p,g)14o
	! 14o(e+nu)14n
	! 14o(a,p)17f cf88 q = 1.191
	call add_screening_factor(zion(io14),aion(io14),zion(ihe4),aion(ihe4))
	! 17f(p,g)18ne
	! 18ne(e+nu)18f
	! 18f(p,a)15o
	! triple alpha to c12
	! c12(a,g)o16
	call add_screening_factor(zion(ic12),aion(ic12),zion(ihe4),aion(ihe4))
	! 16o(a,g)20ne
	call add_screening_factor(zion(io16),aion(io16),zion(ihe4),aion(ihe4))
	! 18ne(a,p)21na
	call add_screening_factor(zion(ine18),aion(ine18),zion(ihe4),aion(ihe4))
	! 19ne(p,g)20na
	call add_screening_factor(zion(ine19),aion(ine19),zion(ihe4),aion(ihe4))
	! 15o(a,g)19ne
	call add_screening_factor(zion(io15),aion(io15),zion(ihe4),aion(ihe4))
	! 44ti(a,p)47v
	call add_screening_factor(22.0d0,44.0d0,zion(ihe4),aion(ihe4))
	! 26si(a,p)29p
	call add_screening_factor(14.0d0,26.0d0,zion(ihe4),aion(ihe4))
	! 19ne(e+nu)19f
	! 20ne(p,g)21na

  end subroutine set_up_screening_factors


  subroutine update_unevolved_species(state)

    implicit none

    type (burn_t)    :: state

  end subroutine update_unevolved_species

end module actual_rhs_module
