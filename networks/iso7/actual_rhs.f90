module actual_rhs_module

  use network
  use network_indices
  use burner_module
  use eos_type_module
  use vode_data, only: net_itemp, net_ienuc
  use rpar_indices

  implicit none

  double precision, parameter :: c54 = 56.0d0/54.0d0

contains

  subroutine actual_rhs(neq,time,y,dydt,rpar)

    use extern_probin_module, only: jacobian

    implicit none

    ! This routine sets up the system of ODE's for the iso7
    ! nuclear reactions.  This is an alpha chain + heavy ion network
    ! with (a,p)(p,g) links up to silicon, as well as a Si <-> Ni link.
    !
    ! Isotopes: he4,  c12,  o16,  ne20, mg24, si28, ni56

    integer          :: neq
    double precision :: time
    double precision :: y(neq), dydt(neq)
    double precision :: rpar(n_rpar_comps)

    integer :: i
    logical :: deriva

    double precision :: ratraw(nrates), dratrawdt(nrates), dratrawdd(nrates)
    double precision :: ratdum(nrates), dratdumdt(nrates), dratdumdd(nrates)
    double precision :: dratdumdy1(nrates), dratdumdy2(nrates)
    double precision :: scfac(nrates),  dscfacdt(nrates),  dscfacdd(nrates)

    double precision :: sneut, dsneutdt, dsneutdd, snuda, snudz
    double precision :: enuc

    double precision :: rho, temp, cv, cp, abar, zbar, dEdY(nspec), dhdY(nspec)

    ! Get the data from rpar and the state

    rho  = rpar(irp_dens)
    temp = y(net_itemp)
    cv   = rpar(irp_cv)
    cp   = rpar(irp_cp)

    dhdY = rpar(irp_dhdY:irp_dhdY+nspec-1)
    dEdY = rpar(irp_dEdY:irp_dEdY+nspec-1)
    abar = rpar(irp_abar)
    zbar = rpar(irp_zbar)

    ! Get the raw reaction rates
    call iso7rat(temp, rho, ratraw, dratrawdt, dratrawdd)

    ! Do the screening here because the corrections depend on the composition
    call screen_iso7(temp, rho, y(1:nspec),        &
                        ratraw, dratrawdt, dratrawdd, &
                        ratdum, dratdumdt, dratdumdd, &
                        dratdumdy1, dratdumdy2,       &
                        scfac,  dscfacdt,  dscfacdd)

    ! Get the right hand side of the ODEs. First, we'll do it
    ! using d(rates)/dT to get the data for the Jacobian and store it.
    ! Then we'll do it using the normal rates.

    if (jacobian == 1) then
       deriva = .true.
       call rhs(y(1:nspec), dratdumdt, ratdum, dydt, deriva)
       rpar(irp_dydt:irp_dydt+nspec-1) = dydt(1:nspec)
       rpar(irp_rates:irp_rates+nrates-1) = ratdum
       rpar(irp_drdy1:irp_drdy1+nrates-1) = dratdumdy1
       rpar(irp_drdy2:irp_drdy2+nrates-1) = dratdumdy2
    endif

    deriva = .false.

    call rhs(y(1:nspec), ratdum, ratdum, dydt, deriva)

    ! Instantaneous energy generation rate
    call ener_gener_rate(dydt, enuc)

    ! Get the neutrino losses
    call sneut5(temp, rho, abar, zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)

    ! Append the energy equation (this is erg/g/s)
    dydt(net_ienuc) = enuc - sneut

    ! Set up the temperature ODE
    call temperature_rhs(neq, y, dydt, rpar)

  end subroutine actual_rhs



  ! Analytical Jacobian

  subroutine actual_jac(neq, t, y, pd, rpar)

    use bl_types
    use bl_constants_module, only: ZERO
    use eos_module
    use extern_probin_module, only: do_constant_volume_burn

    implicit none

    integer         , intent(IN   ) :: neq
    double precision, intent(IN   ) :: y(neq), rpar(n_rpar_comps), t
    double precision, intent(  OUT) :: pd(neq,neq)

    double precision :: ratdum(nrates), dratdumdy1(nrates), dratdumdy2(nrates), ydot(nspec)

    double precision :: b1, sneut, dsneutdt, dsneutdd, snuda, snudz

    integer          :: i, j

    double precision :: rho, temp, cv, cp, abar, zbar, dEdY(nspec), dhdY(nspec)

    pd(:,:) = ZERO

    ! Get the data from rpar and the state

    rho  = rpar(irp_dens)
    temp = y(net_itemp)

    abar = rpar(irp_abar)
    zbar = rpar(irp_zbar)

    ! Note that this RHS has been evaluated using rates = d(ratdum) / dT

    ydot       = rpar(irp_dydt:irp_dydt+nspec-1)
    ratdum     = rpar(irp_rates:irp_rates+nrates-1)
    dratdumdy1 = rpar(irp_drdy1:irp_drdy1+nrates-1)
    dratdumdy2 = rpar(irp_drdy2:irp_drdy2+nrates-1)

    ! Species Jacobian elements with respect to other species

    call dfdy_isotopes_iso7(neq, y, pd, ratdum, dratdumdy1, dratdumdy2)

    ! Energy generation rate Jacobian elements with respect to species

    do j = 1, nspec
       call ener_gener_rate(pd(1:nspec,j), pd(net_ienuc,j))
    enddo

    ! Account for the thermal neutrino losses

    call sneut5(T,rho,abar,zbar,sneut,dsneutdt,dsneutdd,snuda,snudz)

    do j = 1, nspec
       b1 = ((aion(j) - abar) * abar * snuda + (zion(j) - zbar) * abar * snudz)
       pd(net_ienuc,j) = pd(net_ienuc,j) - b1
    enddo

    ! Jacobian elements with respect to temperature

    pd(1:nspec,net_itemp) = ydot

    call ener_gener_rate(pd(1:nspec,net_itemp), pd(net_ienuc,net_itemp))
    pd(net_ienuc,net_itemp) = pd(net_ienuc,net_itemp) - dsneutdt

    ! Temperature Jacobian elements

    call temperature_jac(neq, y, pd, rpar)

  end subroutine actual_jac



  ! Evaluates the right hand side of the iso7 ODEs

  subroutine rhs(y, rate, ratdum, dydt, deriva)

    use bl_constants_module, only: ZERO, SIXTH
    use microphysics_math_module, only: esum

    implicit none

    ! deriva is used in forming the analytic Jacobian to get
    ! the derivative wrt A

    logical          :: deriva
    double precision :: y(nspec),rate(nrates),ratdum(nrates),dydt(nspec)

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

    dydt(ihe4) = esum(a,15)

    ! 12c reactions
    a(1) =  SIXTH * y(ihe4) * y(ihe4) * y(ihe4) * rate(ir3a)
    a(2) = -y(ic12) * rate(irg3a)
    a(3) =  y(io16) * rate(iroga)
    a(4) = -y(ic12) * y(ihe4) * rate(ircag)
    a(5) = -y(ic12) * y(ic12) * rate(ir1212)
    a(6) = -y(ic12) * y(io16) * rate(ir1216)

    dydt(ic12) = esum(a,6)

    ! 16o reactions
    a(1) = -y(io16) * rate(iroga)
    a(2) =  y(ic12) * y(ihe4) * rate(ircag)
    a(3) = -y(ic12) * y(io16) * rate(ir1216)
    a(4) = -y(io16) * y(io16) * rate(ir1616)
    a(5) = -y(io16) * y(ihe4) * rate(iroag)
    a(6) =  y(ine20) * rate(irnega)

    dydt(io16) = esum(a,6)

    ! 20ne reactions
    a(1) =  0.5d0 * y(ic12) * y(ic12) * rate(ir1212)
    a(2) =  y(io16) * y(ihe4) * rate(iroag)
    a(3) = -y(ine20) * rate(irnega)
    a(4) =  y(img24) * rate(irmgga)
    a(5) = -y(ine20) * y(ihe4) * rate(irneag)

    dydt(ine20) = esum(a,5)

    ! 24mg reactions
    a(1) =  0.5d0 * y(ic12) * y(io16) * rate(ir1216)
    a(2) = -y(img24) * rate(irmgga)
    a(3) =  y(ine20) * y(ihe4) * rate(irneag)
    a(4) =  y(isi28) * rate(irsiga)
    a(5) = -y(img24) * y(ihe4) * rate(irmgag)

    dydt(img24) = esum(a,5)

    ! 28si reactions
    a(1) =  0.5d0 * y(ic12) * y(io16) * rate(ir1216)
    a(2) =  0.5d0 * y(io16) * y(io16) * rate(ir1616)
    a(3) = -y(isi28) * rate(irsiga)
    a(4) =  y(img24) * y(ihe4) * rate(irmgag)
    a(5) = -rate(irsi2ni) * y(ihe4)
    a(6) =  rate(irni2si) * y(ini56)

    dydt(isi28) = esum(a,6)

    ! ni56 reactions
    a(1) =  rate(irsi2ni) * y(ihe4)
    a(2) = -rate(irni2si) * y(ini56)

    dydt(ini56) = esum(a,2)

  end subroutine rhs



  subroutine iso7rat(btemp, bden, ratraw, dratrawdt, dratrawdd)

    ! this routine generates unscreened
    ! nuclear reaction rates for the iso7 network.

    use tfactors_module
    use rates_module

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
    tf = get_tfactors(btemp)

    ! c12(a,g)o16
    call rate_c12ag(tf,bden, &
                    ratraw(ircag),dratrawdt(ircag),dratrawdd(ircag), &
                    ratraw(iroga),dratrawdt(iroga),dratrawdd(iroga))

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
                         dratdumdy1, dratdumdy2, &
                         scfac, dscfacdt, dscfacdd)

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
    double precision :: dratdumdy1(nrates), dratdumdy2(nrates)
    double precision :: scfac(nrates),  dscfacdt(nrates),  dscfacdd(nrates)

    integer          :: i, jscr
    double precision :: sc1a,sc1adt,sc1add,sc2a,sc2adt,sc2add, &
                        sc3a,sc3adt,sc3add,abar,zbar,ye,z2bar, &
                        t992,t9i92,yeff_ca40,yeff_ca40dt,yeff_ti44,yeff_ti44dt, &
                        denom,denomdt,denomdd,xx,zz

    type (plasma_state) :: state
    type (tf_t)         :: tf
    
    ! initialize
    do i = 1, nrates
       ratdum(i)     = ratraw(i)
       dratdumdt(i)  = dratrawdt(i)
       dratdumdd(i)  = dratrawdd(i)
       dratdumdy1(i) = ZERO
       dratdumdy2(i) = ZERO
       scfac(i)      = ONE
       dscfacdt(i)   = ZERO
       dscfacdd(i)   = ZERO
    enddo

    ! get the temperature factors
    tf = get_tfactors(btemp)

    ! Set up the state data, which is the same for all screening factors.

    call fill_plasma_state(state, btemp, bden, y(1:nspec))

    ! first the always fun triple alpha and its inverse
    jscr = 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    jscr = jscr + 1
    call screen5(state,jscr,sc2a,sc2adt,sc2add)

    sc3a   = sc1a * sc2a
    sc3adt = sc1adt*sc2a + sc1a*sc2adt
    sc3add = sc1add*sc2a + sc1a*sc2add

    ratdum(ir3a)    = ratraw(ir3a) * sc3a
    dratdumdt(ir3a) = dratrawdt(ir3a)*sc3a + ratraw(ir3a)*sc3adt
    dratdumdd(ir3a) = dratrawdd(ir3a)*sc3a + ratraw(ir3a)*sc3add

    scfac(ir3a)     = sc3a
    dscfacdt(ir3a)  = sc3adt
    dscfacdd(ir3a)  = sc3add


    ! c12 to o16
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ratdum(ircag)     = ratraw(ircag) * sc1a
    dratdumdt(ircag)  = dratrawdt(ircag)*sc1a + ratraw(ircag)*sc1adt
    dratdumdd(ircag)  = dratrawdd(ircag)*sc1a + ratraw(ircag)*sc1add

    scfac(ircag)      = sc1a
    dscfacdt(ircag)   = sc1adt
    dscfacdt(ircag)   = sc1add



    ! c12 + c12
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ratdum(ir1212)    = ratraw(ir1212) * sc1a
    dratdumdt(ir1212) = dratrawdt(ir1212)*sc1a + ratraw(ir1212)*sc1adt
    dratdumdd(ir1212) = dratrawdd(ir1212)*sc1a + ratraw(ir1212)*sc1add

    scfac(ir1212)     = sc1a
    dscfacdt(ir1212)  = sc1adt
    dscfacdd(ir1212)  = sc1add



    ! c12 + o16
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ratdum(ir1216)    = ratraw(ir1216) * sc1a
    dratdumdt(ir1216) = dratrawdt(ir1216)*sc1a + ratraw(ir1216)*sc1adt
    dratdumdd(ir1216) = dratrawdd(ir1216)*sc1a + ratraw(ir1216)*sc1add

    scfac(ir1216)     = sc1a
    dscfacdt(ir1216)  = sc1adt
    dscfacdd(ir1216)  = sc1add



    ! o16 + o16
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ratdum(ir1216)    = ratraw(ir1216) * sc1a
    dratdumdt(ir1216) = dratrawdt(ir1216)*sc1a + ratraw(ir1216)*sc1adt
    dratdumdd(ir1216) = dratrawdd(ir1216)*sc1a + ratraw(ir1216)*sc1add

    scfac(ir1216)     = sc1a
    dscfacdt(ir1216)  = sc1adt
    dscfacdd(ir1216)  = sc1add



    ! o16 to ne20
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ratdum(iroag)    = ratraw(iroag) * sc1a
    dratdumdt(iroag) = dratrawdt(iroag)*sc1a + ratraw(iroag)*sc1adt
    dratdumdd(iroag) = dratrawdd(iroag)*sc1a + ratraw(iroag)*sc1add

    scfac(iroag)     = sc1a
    dscfacdt(iroag)  = sc1adt
    dscfacdd(iroag)  = sc1add



    ! o16 to mg24
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ratdum(irneag)    = ratraw(irneag) * sc1a
    dratdumdt(irneag) = dratrawdt(irneag)*sc1a + ratraw(irneag)*sc1adt
    dratdumdd(irneag) = dratrawdd(irneag)*sc1a + ratraw(irneag)*sc1add

    scfac(irneag)     = sc1a
    dscfacdt(irneag)  = sc1adt
    dscfacdd(irneag)  = sc1add


    ! mg24 to si28
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ratdum(irmgag)    = ratraw(irmgag) * sc1a
    dratdumdt(irmgag) = dratrawdt(irmgag)*sc1a + ratraw(irmgag)*sc1adt
    dratdumdd(irmgag) = dratrawdd(irmgag)*sc1a + ratraw(irmgag)*sc1add

    scfac(irmgag)     = sc1a
    dscfacdt(irmgag)  = sc1adt
    dscfacdd(irmgag)  = sc1add



    ! ca40 to ti44
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ratdum(ircaag)    = ratraw(ircaag) * sc1a
    dratdumdt(ircaag) = dratrawdt(ircaag)*sc1a + ratraw(ircaag)*sc1adt
    dratdumdd(ircaag) = dratrawdd(ircaag)*sc1a + ratraw(ircaag)*sc1add

    scfac(ircaag)     = sc1a
    dscfacdt(ircaag)  = sc1adt
    dscfacdd(ircaag)  = sc1add



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
       dratdumdy1(irsi2ni) = 3.0d0 * ratdum(irsi2ni)/y(ihe4)
       dratdumdy2(irsi2ni) = yeff_ca40*denom*ratdum(ircaag)
       dratdumdt(irsi2ni)  = (yeff_ca40dt*ratdum(ircaag) &
            + yeff_ca40*dratdumdt(ircaag))*denom*y(isi28)*1.0d-9
       dratdumdd(irsi2ni)  = 3.0d0*ratdum(irsi2ni)/bden &
            + yeff_ca40*denom*dratdumdd(ircaag)*y(isi28)


       if (denom .ne. 0.0) then

          zz     = 1.0d0/denom
          ratdum(irni2si) = min(1.0d10,yeff_ti44*ratdum(irtiga)*zz)

          if (ratdum(irni2si) .eq. 1.0d10) then
             dratdumdy1(irni2si) = 0.0d0
             dratdumdt(irni2si)  = 0.0d0
             dratdumdd(irni2si)  = 0.0d0

          else
             dratdumdy1(irni2si) = -3.0d0 * ratdum(irni2si)/y(ihe4)
             dratdumdt(irni2si)  = (yeff_ti44dt*ratdum(irtiga) &
                  + yeff_ti44*dratdumdt(irtiga))*zz*1.0d-9
             dratdumdd(irni2si)  = -3.0d0 * ratdum(irni2si)/bden &
                  + yeff_ti44*dratdumdd(irtiga)*zz
          end if
       endif
    end if

  end subroutine screen_iso7



  subroutine dfdy_isotopes_iso7(neq,y,dfdy,ratdum,dratdumdy1,dratdumdy2)

    use network
    use microphysics_math_module, only: esum

    implicit none

    ! this routine sets up the dense iso7 jacobian for the isotopes

    integer          :: neq
    double precision :: y(neq), dfdy(neq,neq)
    double precision :: ratdum(nrates), dratdumdy1(nrates), dratdumdy2(nrates)

    double precision :: b(8)
    
    ! set up the jacobian
    ! 4he jacobian elements
    ! d(he4)/d(he4)
    b(1) = -1.5d0 * y(ihe4) * y(ihe4) * ratdum(ir3a)
    b(2) = -y(ic12) * ratdum(ircag)
    b(3) = -y(io16) * ratdum(iroag)
    b(4) = -y(ine20) * ratdum(irneag)
    b(5) = -y(img24) * ratdum(irmgag)
    b(6) = -7.0d0 * ratdum(irsi2ni)
    b(7) = -7.0d0 * dratdumdy1(irsi2ni) * y(ihe4)
    b(8) =  7.0d0 * dratdumdy1(irni2si) * y(ini56)

    dfdy(ihe4,ihe4) = esum(b,8)

    ! d(he4)/d(c12)
    b(1) =  3.0d0 * ratdum(irg3a)
    b(2) = -y(ihe4) * ratdum(ircag)
    b(3) =  y(ic12) * ratdum(ir1212)
    b(4) =  0.5d0 * y(io16) * ratdum(ir1216)

    dfdy(ihe4,ic12) = esum(b,4)

    ! d(he4)/d(o16)
    b(1) =  ratdum(iroga)
    b(2) =  0.5d0 * y(ic12) * ratdum(ir1216)
    b(3) =  y(io16) * ratdum(ir1616)
    b(4) = -y(ihe4) * ratdum(iroag)

    dfdy(ihe4,io16) = esum(b,4)

    ! d(he4)/d(ne20)
    b(1) =  ratdum(irnega)
    b(2) = -y(ihe4) * ratdum(irneag)

    dfdy(ihe4,ine20) = esum(b,2)
                      
    ! d(he4)/d(mg24)
    b(1) =  ratdum(irmgga)
    b(2) = -y(ihe4) * ratdum(irmgag)

    dfdy(ihe4,img24) = esum(b,2)

    ! d(he4)/d(si28)
    b(1) =  ratdum(irsiga)
    b(2) = -7.0d0 * dratdumdy2(irsi2ni) * y(ihe4)

    dfdy(ihe4,isi28) = esum(b,2)

    ! d(he4)/d(ni56)
    b(1) =  7.0d0 * ratdum(irni2si)

    dfdy(ihe4,ini56) = esum(b,1)



    ! 12c jacobian elements
    ! d(c12)/d(he4)
    b(1) =  0.5d0 * y(ihe4) * y(ihe4) * ratdum(ir3a)
    b(2) = -y(ic12) * ratdum(ircag)

    dfdy(ic12,ihe4) = esum(b,2)

    ! d(c12)/d(c12)
    b(1) = -ratdum(irg3a)
    b(2) = -y(ihe4) * ratdum(ircag)
    b(3) = -2.0d0 * y(ic12) * ratdum(ir1212)
    b(4) = -y(io16) * ratdum(ir1216)

    dfdy(ic12,ic12) = esum(b,4)

    ! d(c12)/d(o16)
    b(1) =  ratdum(iroga)
    b(2) = -y(ic12) * ratdum(ir1216)

    dfdy(ic12,io16) = esum(b,2)


    ! 16o jacobian elements
    ! d(o16)/d(he4)
    b(1) =  y(ic12) * ratdum(ircag)
    b(2) = -y(io16) * ratdum(iroag)

    dfdy(io16,ihe4) = esum(b,2)

    ! d(o16)/d(c12)
    b(1) =  y(ihe4) * ratdum(ircag)
    b(2) = -y(io16) * ratdum(ir1216)

    dfdy(io16,ic12) = esum(b,2)

    ! d(o16)/d(o16)
    b(1) = -ratdum(iroga)
    b(2) = -y(ic12) * ratdum(ir1216)
    b(3) = -2.0d0 * y(io16) * ratdum(ir1616)
    b(4) = -y(ihe4) * ratdum(iroag)

    dfdy(io16,io16) = esum(b,4)

    ! d(o16)/d(ne20)
    b(1) =  ratdum(irnega)

    dfdy(io16,ine20) = esum(b,1)



    ! 20ne jacobian elements
    ! d(ne20)/d(he4)
    b(1) =  y(io16) * ratdum(iroag) - y(ine20) * ratdum(irneag)

    dfdy(ine20,ihe4) = esum(b,1)

    ! d(ne20)/d(c12)
    b(1) =  y(ic12) * ratdum(ir1212)

    dfdy(ine20,ic12) = esum(b,1)

    ! d(ne20)/d(o16)
    b(1) =  y(ihe4) * ratdum(iroag)

    dfdy(ine20,io16) = esum(b,1)

    ! d(ne20)/d(ne20)
    b(1) = -ratdum(irnega) - y(ihe4) * ratdum(irneag)

    dfdy(ine20,ine20) = esum(b,1)

    ! d(ne20)/d(mg24)
    b(1) =  ratdum(irmgga)

    dfdy(ine20,img24) = esum(b,1)



    ! 24mg jacobian elements
    ! d(mg24)/d(he4)
    b(1) =  y(ine20) * ratdum(irneag)
    b(2) = -y(img24) * ratdum(irmgag)

    dfdy(img24,ihe4) = esum(b,2)

    ! d(mg24)/d(c12)
    b(1) =  0.5d0 * y(io16) * ratdum(ir1216)

    dfdy(img24,ic12) = esum(b,1)

    ! d(mg24)/d(o16)
    b(1) =  0.5d0 * y(ic12) * ratdum(ir1216)

    dfdy(img24,io16) = esum(b,1)

    ! d(mg24)/d(ne20)
    b(1) =  y(ihe4) * ratdum(irneag)

    dfdy(img24,ine20) = esum(b,1)

    ! d(mg24)/d(mg24)
    b(1) = -ratdum(irmgga)
    b(2) = -y(ihe4) * ratdum(irmgag)

    dfdy(img24,img24) = esum(b,2)

    ! d(mg24)/d(si28)
    b(1) =  ratdum(irsiga)

    dfdy(img24,isi28) = esum(b,1)

    ! 28si jacobian elements
    ! d(si28)/d(he4)
    b(1) =  y(img24) * ratdum(irmgag)
    b(2) = -ratdum(irsi2ni)
    b(3) = -dratdumdy1(irsi2ni) * y(ihe4)
    b(4) =  dratdumdy1(irni2si) * y(ini56)

    dfdy(isi28,ihe4) = esum(b,4)

    ! d(si28)/d(c12)
    b(1) =  0.5d0 * y(io16) * ratdum(ir1216)

    dfdy(isi28,ic12) = esum(b,1)

    ! d(si28)/d(o16)
    b(1) =  y(io16) * ratdum(ir1616)
    b(2) =  0.5d0 * y(ic12) * ratdum(ir1216)

    dfdy(isi28,io16) = esum(b,2)

    ! d(si28)/d(mg24)
    b(1) =  y(ihe4) * ratdum(irmgag)

    dfdy(isi28,img24) = esum(b,1)

    ! d(si28)/d(si28)
    b(1) = -ratdum(irsiga)
    b(2) = -dratdumdy2(irsi2ni) * y(ihe4)

    dfdy(isi28,isi28) = esum(b,2)

    ! d(si28)/d(ni56)
    b(1) =  ratdum(irni2si)

    dfdy(isi28,ini56) = esum(b,1)

    ! ni56 jacobian elements
    ! d(ni56)/d(he4)
    b(1) =  ratdum(irsi2ni)
    b(2) =  dratdumdy1(irsi2ni) * y(ihe4)
    b(3) = -dratdumdy1(irni2si) * y(ini56)

    dfdy(ini56,ihe4) = esum(b,3)

    ! d(ni56)/d(si28)
    b(1) = dratdumdy2(irsi2ni) * y(ihe4)

    dfdy(ini56,isi28) = esum(b,1)
    
    ! d(ni56)/d(ni56)
    b(1) = -ratdum(irni2si)

    dfdy(ini56,ini56) = esum(b,1)

  end subroutine dfdy_isotopes_iso7

end module actual_rhs_module