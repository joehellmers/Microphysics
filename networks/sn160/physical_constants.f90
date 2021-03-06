module physical_constants

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt), parameter :: pi = 3.141592653589793238462643383279502884197d0
  real(rt), parameter :: cm_per_fm = 1.0d-13

  ! From 2014 CODATA recommended values
  real(rt), parameter :: clight  = 299792458d2  ! cm/s

  real(rt), parameter :: N_avo = 6.022140857d23 ! 1/mol

  real(rt), parameter :: gram_per_amu = 1660539.040d-30

  real(rt), parameter :: erg_per_eV  = 1.6021766208d-12
  real(rt), parameter :: erg_per_MeV = erg_per_eV*1.0d6
  real(rt), parameter :: gram_per_MeV = erg_per_MeV/CLIGHT**2

end module physical_constants
