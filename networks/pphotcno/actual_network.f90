module actual_network

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, parameter :: nspec = 26
  integer, parameter :: nspec_evolve = 26
  integer, parameter :: naux  = 0

  integer, parameter :: ih1   = 1
  integer, parameter :: ih2   = 2
  integer, parameter :: ihe3  = 3
  integer, parameter :: ihe4  = 4
  integer, parameter :: ili7  = 5
  integer, parameter :: ibe7  = 6
  integer, parameter :: ib8   = 7
  integer, parameter :: ic12  = 8
  integer, parameter :: ic13  = 9
  integer, parameter :: in13  = 10
  integer, parameter :: in14  = 11
  integer, parameter :: in15  = 12
  integer, parameter :: io14  = 13
  integer, parameter :: io15  = 14
  integer, parameter :: io16  = 15
  integer, parameter :: io17  = 16
  integer, parameter :: io18  = 17
  integer, parameter :: if17  = 18
  integer, parameter :: if18  = 19
  integer, parameter :: if19  = 20
  integer, parameter :: ine18 = 21
  integer, parameter :: ine19 = 22
  integer, parameter :: ine20 = 23
  integer, parameter :: img22 = 24
  integer, parameter :: is30  = 25
  integer, parameter :: ini56 = 26

  double precision, save :: aion(nspec), zion(nspec), nion(nspec)
  double precision, save :: bion(nspec), mion(nspec), wion(nspec)

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  character (len=32), parameter :: network_name = "pphotcno"

  ! Some fundamental physical constants

  double precision, parameter :: avo = 6.0221417930d23
  double precision, parameter :: c_light = 2.99792458d10

  double precision, parameter :: ev2erg  = 1.60217648740d-12
  double precision, parameter :: mev2erg = ev2erg*1.0d6
  double precision, parameter :: mev2gr  = mev2erg/c_light**2

  double precision, parameter :: mn = 1.67492721184d-24
  double precision, parameter :: mp = 1.67262163783d-24
  double precision, parameter :: me = 9.1093821545d-28

  ! Conversion factor for the nuclear energy generation rate.

  double precision, parameter :: enuc_conv2 = -avo*c_light*c_light

  ! Rates data

  integer, parameter :: nrates  = 17
  integer, parameter :: num_rate_groups = 4

  integer, parameter :: ircag   = 1
  integer, parameter :: iroga   = 2
  integer, parameter :: ir3a    = 3
  integer, parameter :: irg3a   = 4
  integer, parameter :: ir1212  = 5
  integer, parameter :: ir1216  = 6
  integer, parameter :: ir1616  = 7
  integer, parameter :: iroag   = 8
  integer, parameter :: irnega  = 9
  integer, parameter :: irneag  = 10
  integer, parameter :: irmgga  = 11
  integer, parameter :: irmgag  = 12
  integer, parameter :: irsiga  = 13
  integer, parameter :: ircaag  = 14
  integer, parameter :: irtiga  = 15
  integer, parameter :: irsi2ni = 16
  integer, parameter :: irni2si = 17

  character (len=20), save :: ratenames(nrates)

contains

  subroutine actual_network_init

    implicit none

    short_spec_names(ih1)  = 'h1'
    short_spec_names(ih2)  = 'h2'
    short_spec_names(ihe3) = 'he3'
    short_spec_names(ihe4) = 'he4'
    short_spec_names(ili7) = 'li7'
    short_spec_names(ibe7) = 'be7'
    short_spec_names(ib8)  = 'b8'
    short_spec_names(ic12) = 'c12'
    short_spec_names(ic13) = 'c13'
    short_spec_names(in13) = 'n13'
    short_spec_names(in14) = 'n14'
    short_spec_names(in15) = 'n15'
    short_spec_names(io14) = 'o14'
    short_spec_names(io15) = 'o15'
    short_spec_names(io16) = 'o16'
    short_spec_names(io17) = 'o17'
    short_spec_names(io18) = 'o18'
    short_spec_names(if17) = 'f17'
    short_spec_names(if18) = 'f18'
    short_spec_names(if19) = 'f19'
    short_spec_names(ine18)= 'ne18'
    short_spec_names(ine19)= 'ne19'
    short_spec_names(ine20)= 'ne20'
    short_spec_names(img22)= 'mg22'
    short_spec_names(is30) = 's30'
    short_spec_names(ini56)= 'ni56'

    spec_names(ih1)  = 'hydrogen-1'
    spec_names(ih2)  = 'hydrogen-2'
    spec_names(ihe3) = 'helium-3'
    spec_names(ihe4) = 'helium-4'
    spec_names(ili7) = 'lithium-7'
    spec_names(ibe7) = 'beryllium-7'
    spec_names(ib8)  = 'boron-8'
    spec_names(ic12) = 'carbon-12'
    spec_names(ic13) = 'carbon-13'
    spec_names(in13) = 'nitrogen-13'
    spec_names(in14) = 'nitrogen-14'
    spec_names(in15) = 'nitrogen-15'
    spec_names(io14) = 'oxygen-14'
    spec_names(io15) = 'oxygen-15'
    spec_names(io16) = 'oxygen-16'
    spec_names(io17) = 'oxygen-17'
    spec_names(io18) = 'oxygen-18'
    spec_names(if17) = 'fluorine-17'
    spec_names(if18) = 'fluorine-18'
    spec_names(if19) = 'fluorine-19'
    spec_names(ine18)= 'neon-18'
    spec_names(ine19)= 'neon-19'
    spec_names(ine20)= 'neon-20'
    spec_names(img22)= 'magnesium-22'
    spec_names(is30) = 'sulfur-30'
    spec_names(ini56)= 'nickel-56'

    ! Set the number of nucleons in the element
    aion(ih1)  = 1.0d0
    aion(ih2)  = 2.0d0
    aion(ihe3) = 3.0d0
    aion(ihe4) = 4.0d0
    aion(ili7) = 7.0d0
    aion(ibe7) = 7.0d0
    aion(ib8)  = 8.0d0
    aion(ic12) = 12.0d0
    aion(ic13) = 13.0d0
    aion(in13) = 13.0d0
    aion(in14) = 14.0d0
    aion(in15) = 15.0d0
    aion(io14) = 14.0d0
    aion(io15) = 15.0d0
    aion(io16) = 16.0d0
    aion(io17) = 17.0d0
    aion(io18) = 18.0d0
    aion(if17) = 17.0d0
    aion(if18) = 18.0d0
    aion(if19) = 19.0d0
    aion(ine18)= 18.0d0
    aion(ine19)= 19.0d0
    aion(ine20)= 20.0d0
    aion(img22)= 22.0d0
    aion(is30) = 30.0d0
    aion(ini56)= 56.0d0

    ! Set the number of protons in the element
    zion(ih1)  = 1.0d0
    zion(ih2)  = 1.0d0
    zion(ihe3) = 2.0d0
    zion(ihe4) = 2.0d0
    zion(ili7) = 3.0d0
    zion(ibe7) = 4.0d0
    zion(ib8)  = 5.0d0
    zion(ic12) = 6.0d0
    zion(ic13) = 6.0d0
    zion(in13) = 7.0d0
    zion(in14) = 7.0d0
    zion(in15) = 7.0d0
    zion(io14) = 8.0d0
    zion(io15) = 8.0d0
    zion(io16) = 8.0d0
    zion(io17) = 8.0d0
    zion(io18) = 8.0d0
    zion(if17) = 9.0d0
    zion(if18) = 9.0d0
    zion(if19) = 9.0d0
    zion(ine18)= 10.0d0
    zion(ine19)= 10.0d0
    zion(ine20)= 10.0d0
    zion(img22)= 12.0d0
    zion(is30) = 16.0d0
    zion(ini56)= 28.0d0

    ! Set the binding energy of the element

    ! set the binding energy of the element
    bion(ih1)   =   0.0d0
    bion(ih2)   =   2.2250d0
    bion(ihe3)  =   7.7204d0
    bion(ili7)  =  39.2440d0
    bion(ibe7)  =  37.6000d0
    bion(ib8)   =  37.7380d0
    bion(ihe4)  =  28.2928d0
    bion(ic12)  =  92.1624d0
    bion(ic13)  =  97.1088d0
    bion(in13)  =  94.1064d0
    bion(in14)  = 104.6598d0
    bion(in15)  = 115.4932d0
    bion(io14)  =  98.7324d0
    bion(io15)  = 111.9558d0
    bion(io16)  = 127.6202d0
    bion(io17)  = 131.7636d0
    bion(io18)  = 139.8080d0
    bion(if17)  = 128.2212d0
    bion(if18)  = 137.3706d0
    bion(if19)  = 147.8020d0
    bion(ine18) = 132.1390d0
    bion(ine19) = 143.7780d0
    bion(ine20) = 160.64788d0
    bion(img22) = 168.57680d0
    bion(is30)  = 243.68660d0
    bion(ini56) = 483.99500d0


    ! Set the number of neutrons
    nion(:) = aion(:) - zion(:)

    ! Set the mass
    mion(:) = nion(:) * mn + zion(:) * (mp + me) - bion(:) * mev2gr

    ! Molar mass
    wion(:) = avo * mion(:)

    ! Common approximation
    wion(:) = aion(:)

    ! set the names of the reaction rates
    ratenames(ircag)   = 'rcag '
    ratenames(iroga)   = 'roga '
    ratenames(ir3a)    = 'r3a  '
    ratenames(irg3a)   = 'rg3a '    ! inverse rate
    ratenames(ir1212)  = 'r1212'
    ratenames(ir1216)  = 'r1216'
    ratenames(ir1616)  = 'r1616'
    ratenames(iroag)   = 'roag '
    ratenames(irnega)  = 'rnega'
    ratenames(irneag)  = 'rneag'
    ratenames(irmgga)  = 'rmgga'
    ratenames(irmgag)  = 'rmgag'
    ratenames(irsiga)  = 'rsiga'
    ratenames(ircaag)  = 'rcaag'
    ratenames(irtiga)  = 'rtiga'

    ratenames(irsi2ni) = 'rsi2ni'
    ratenames(irni2si) = 'rni2si'

  end subroutine actual_network_init



  subroutine actual_network_finalize

    implicit none

    ! Nothing to do here.

  end subroutine actual_network_finalize

end module actual_network
