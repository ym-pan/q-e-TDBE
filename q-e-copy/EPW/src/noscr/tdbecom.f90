  !--------------------------------------------------------------------------
  MODULE tdbe_eph_common
  !--------------------------------------------------------------------------
  !!
  !! Global variables for time-dependent Boltzmann simulations
  !! and Semiconductor Bloch equations which will be implemented
  !! in the future.
  !! Implemented by Yiming Pan
  !!
  !--------------------------------------------------------------------------
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !  
  REAL(KIND = DP), ALLOCATABLE :: enk_all(:, :)
  !! Electron energies 
  REAL(KIND = DP), ALLOCATABLE :: wnuq_all(:, :)
  !! Phonon freqencies
  REAL(KIND = DP), ALLOCATABLE :: wkf_all(:) 
  !! weights of k points
  REAL(KIND = DP), ALLOCATABLE :: wqf_all(:) 
  !! weights of q points
  REAL(KIND = DP), ALLOCATABLE :: xkf_all(:,:)
  !! crystal coordinates of k points
  REAL(KIND = DP), ALLOCATABLE :: xqf_all(:,:)
  !! crystal coordinates of q points
  REAL(KIND = DP), ALLOCATABLE :: vnk_all(:,:,:)
  !! Electron velcity
  REAL(KIND = DP), ALLOCATABLE :: vnuq_all(:,:,:)
  !! phonon velcity
  REAL(KIND = DP) :: ef0
  !! Fermi energy
  COMPLEX(KIND = DP), ALLOCATABLE :: uf_all(:,:, :)
  !! Phonon eigenvectors
  !! Fermi energy
  INTEGER:: totq
  !! Total number of q-points within the fsthick window.
  INTEGER, ALLOCATABLE :: selecq(:)
  !! Array of selected q-points for e-ph scattering
  INTEGER, ALLOCATABLE :: bztoibz_dg(:)
  !! Map from k point to its irr-k point
  INTEGER, ALLOCATABLE :: s_bztoibz_dg(:)
  !! The symmetry operation mapping k and irr-k
  INTEGER :: nkirr_dg
  !! Number of irrducible k points on double grid
  REAL(KIND = DP), ALLOCATABLE :: ie_ph(:,:,:)
  !! Collision integral for e-ph interaction
  REAL(KIND = DP), ALLOCATABLE :: iph_e(:,:,:)
  !! Collision integral for e-ph interaction
  REAL(KIND = DP), ALLOCATABLE ::  iph_ph(:,:,:) 
  !! Collision integral for ph-ph interaction
  REAL(KIND = DP), ALLOCATABLE :: ie_e(:,:,:)
  !! Collision integral for el-el interaction
  COMPLEX(KIND = DP), ALLOCATABLE :: rho_e(:,:,:)
  !! Density matrix for electrons rho_e(nbndfst,nbndfst,nktotf)
  ! REAL (KIND = DP), ALLOCATABLE :: Qnu(:), Pnu(:)
  ! !! Cohrent phonon amplitude
  !! For future development of semiconductor Bloch equation
  REAL(KIND = DP), ALLOCATABLE :: felec(:,:)
  !! felec(nbndfst,nkfs) is the electron occupation at time t
  REAL(KIND = DP), ALLOCATABLE :: nphon(:,:) 
  !! nphon(nmodes,totq) is the phnonon occupation at time t
  REAL(KIND = DP), ALLOCATABLE :: nphon_pre(:,:) 
  !! Phonon distribution at the previous ph-ph time step
  INTEGER, ALLOCATABLE :: iq2nqfs(:, :)
  !! map
  INTEGER, ALLOCATABLE :: indx_mapk(:)
  !! map
  INTEGER, ALLOCATABLE :: indx_mapq(:)
  !! map
  INTEGER :: nkfsden
  !! number of k-points
  INTEGER, ALLOCATABLE :: ikfsdx(:)
  !!
  REAL(KIND = DP), PARAMETER :: fs2sec   = 1.0E-15_DP
  !! femtosecond to second
  REAL(KIND = DP), PARAMETER :: ps2sec   = 1.0E-12_DP
  !! picosecond to second
  INTEGER :: nstg
  !! Number of stages for Runge-Kutta mathods
  REAL(KIND = DP) :: b_tab(6)
  !! Butcher Tableaux parameters b_tab
  REAL(KIND = DP) :: c_tab(6)
  !! Butcher Tableaux parameters c_tab
  REAL(KIND = DP) :: a_tab(6,6)
  !! Butcher Tableaux parameters a_tab
  REAL(KIND = DP) :: tstart
  !! The starting time of simulation
  !--------------------------------------------------------------------------
  END MODULE tdbe_eph_common
  !--------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------
  MODULE tdbe_phph_common
  !--------------------------------------------------------------------------
  !!
  !! Global variables for time-dependent Boltzmann simulations
  !! and Semiconductor Bloch equations which will be implemented
  !! in the future.
  !!
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  INTEGER :: nifc
  !! Number of anharmonic interatomic force constants
  INTEGER, ALLOCATABLE :: ind_atm1(:)
  !! indexes of 1st atoms
  INTEGER, ALLOCATABLE :: ind_atm2(:)
  !! indexes of 2nd atoms
  INTEGER, ALLOCATABLE :: ind_atm3(:)
  !! indexes of 3rd atoms
  REAL(KIND = DP), ALLOCATABLE :: Phi(:,:,:,:)
  !! 3rd force constants
  REAL(KIND = DP), ALLOCATABLE :: rvec2(:,:)
  !! vectors of 2nd unit cells in crystal coordinates
  REAL(KIND = DP), ALLOCATABLE :: rvec3(:,:)
  !! vectors of 3rd unit cells in crystal coordinates
  REAL(kind = DP),ALLOCATABLE :: psi2_p(:)
  !! phonon-phonon matrix elements multiplied by the delta function
  REAL(kind = DP),ALLOCATABLE :: psi2_m(:)
  !! phonon-phonon matrix elements multiplied by the delta function
  INTEGER, ALLOCATABLE :: ind_p(:,:)
  !! indexes of the ph-ph matrix elements
  INTEGER, ALLOCATABLE :: ind_m(:,:)
  !! indexes of the ph-ph matrix elements
  REAL(KIND = DP) :: lattv(3,3)
  !! lattice constants
  REAL(KIND = DP) :: rlattv(3,3)
  !! reciprocal lattice constant
  REAL(KIND = DP) :: lfactor
  !! factor that transforms lattice constant to atomic unit
  REAL(KIND = DP) :: Vcell
  !! volume of cell
  REAL(KIND = DP) :: scalebroad
  !!
  REAL(KIND = DP) :: dt_in_ps
  !! time step in picosecond
  INTEGER(KIND = 8) :: nind_p
  !! Number of ph-ph matrix elements
  INTEGER(KIND = 8) :: nind_m
  !! Number of ph-ph matrix elements
  REAL(KIND = DP),allocatable :: ph_lw(:,:)
  !! phonon linewidth 
  INTEGER, parameter :: ntph = 50000
  !!
  REAL(KIND = DP) :: phtemp(ntph), e_latt(ntph)
  !! array of temperatures
  integer :: istart
!   real(kind=DP), parameter :: hbarp=1.05457172647d-22*1e22
!   real(kind=DP), parameter :: massfactor=1.8218779*6.022e-4
!   real(kind=DP), parameter :: Ry2THz=20670.687 !! ry2THz
  !--------------------------------------------------------------------------
  END MODULE tdbe_phph_common
  !--------------------------------------------------------------------------
  !
  !
  !--------------------------------------------------------------------------
  MODULE tdbecom
  !--------------------------------------------------------------------------
  !
  USE tdbe_eph_common
  USE tdbe_phph_common
  !
  !--------------------------------------------------------------------------
  END MODULE tdbecom
  !--------------------------------------------------------------------------
