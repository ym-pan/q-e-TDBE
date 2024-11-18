  !----------------------------------------------------------------------------
  MODULE phph
  !----------------------------------------------------------------------------
  !!
  !! This module contains routines related to phonon-phonon coupling
  !! Codes written by Yiming Pan
  !! Todo : In a serial running, the calling of MPI subroutines need
  !!        to be repalced.
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !--------------------------------------------------------------------------   
    SUBROUTINE psi2_calc()
    !--------------------------------------------------------------------------
    !! This subroutine reads the thrid order force constants produced with 
    !! thirdorder.py and calculate the ph-ph matrix elements on the dense
    !! q grid, according to Eq.(9) of Comp. Phys. Commun. 185, 17471758 (2014)
    !!
    !--------------------------------------------------------------------------
    !
!    USE mpi
    USE kinds,            ONLY : DP
    USE io_global,        ONLY : stdout, ionode,ionode_id
    USE mp_global,        ONLY : inter_pool_comm, npool, my_pool_id
    USE mp_world,         ONLY : mpime
    USE mp,               ONLY : mp_barrier, mp_sum, mp_bcast, mp_gather, &
                                 mp_max, mp_min
    USE ions_base,        ONLY : nat, amass, ityp, tau 
    USE cell_base,        ONLY : at, bg, alat
    USE epwcom,           ONLY : nqf1, nqf2, nqf3, degaussq, phwmin,      &
                                 double_grid, facbroad, tstep, dph
    USE constants_epw,    ONLY : kelvin2Ry, twopi, hbar, ryd2ev, rydcm1,  &
                                 ryd2mev 
    USE constants,        ONLY : pi
    USE elph2,            ONLY : wf,xqf,nqtotf,adapt_smearing
    USE tdbecom,          ONLY : nifc, ind_atm1, ind_atm2, ind_atm3, Phi, &
                                 rvec2, rvec3, uf_all, vnuq_all,          &
                                 ps2sec, dt_in_ps, lattv, rlattv, Vcell,  &
                                 lfactor, scalebroad, nind_p, nind_m,     &
                                 ind_p, ind_m, psi2_p, psi2_m
    USE modes,            ONLY : nmodes
    USE division,         ONLY : fkbounds, fkbounds2
    USE io_var,           ONLY : iunphindx
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_MODE_WRONLY, MPI_MODE_CREATE,MPI_INFO_NULL, &
                                 MPI_OFFSET_KIND, MPI_DOUBLE_PRECISION,          &
                                 MPI_STATUS_IGNORE, MPI_INTEGER, MPI_IN_PLACE,   &
                                 MPI_MODE_RDONLY, MPI_INTEGER8, MPI_SUM,         &
                                 MPI_MODE_DELETE_ON_CLOSE
#endif
    ! 
    !   
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    INTEGER :: ios 
    !! io status
    LOGICAL :: exst
    !! Find if a file exists.
    INTEGER, ALLOCATABLE :: indxp_tmp(:, :)
    !! 
    INTEGER, ALLOCATABLE :: indxm_tmp(:, :)
    !!     
    INTEGER :: ii
    !! direction index in x axis
    INTEGER :: jj
    !! direction index in y axis
    INTEGER :: kk
    !! direction index in z axis
    INTEGER :: iq1
    !! q-index for 1st phonon
    INTEGER :: iq2
    !! q-index for 2nd phonon
    INTEGER :: iq3p 
    !! q-index for 3rd phonon
    INTEGER :: iq3m 
    !! q-index for 3rd phonon
    INTEGER :: nu1
    !! branch index for 1st phonon
    INTEGER :: nu2 
    !! branch index for 2nd phonon
    INTEGER :: nu3 
    !! branch index for 3rd phonon
    INTEGER :: iqx 
    !! direction index in x axis for phonon
    INTEGER :: iqy
    !! direction index in y axis for phonon
    INTEGER :: iqz 
    !! direction index in z axis for phonon
    INTEGER :: ipool 
    !! index for cpu
    INTEGER :: idim
    !! index for dimension
    INTEGER(KIND = 8) :: N_plus(npool) 
    !! Number of ph-ph matrix elements that needs to be
    INTEGER(KIND = 8) :: N_minus(npool) 
    !! Number of ph-ph matrix elements that needs to be
    !! calculated in each pool
    INTEGER(KIND = 8) :: ind1
    !! Counter on phonon modes
    INTEGER(KIND = 8) :: ind2
    !! Counter on phonon modes
    INTEGER(KIND = 8) :: Ntot_p
    !! Total number of ph-ph matrix elements
    INTEGER(KIND = 8) :: Ntot_m
    !! Total number of ph-ph matrix elements
    INTEGER(KIND = MPI_OFFSET_KIND) :: Ntot_p_
    !! Total number of ph-ph matrix elements
    INTEGER(KIND = MPI_OFFSET_KIND) :: Ntot_m_
    !! Total number of ph-ph matrix elements
    INTEGER(KIND = MPI_OFFSET_KIND) :: displ1
    !! displacements in order to balance 
    !! the ph-ph matrix elements between the pools
    INTEGER(KIND = MPI_OFFSET_KIND) :: displ2
    !! displacements in order to balance 
    !! the ph-ph matrix elements between the pools
    INTEGER(KIND = MPI_OFFSET_KIND) :: nindx
    !! number of elements written to file
    INTEGER(KIND = MPI_OFFSET_KIND) :: nindxi 
    !! number of elements read from file
    INTEGER(KIND = MPI_OFFSET_KIND) :: lowerplus
    !! lower bound for parallelization of ph-ph matrix elements
    INTEGER(KIND = MPI_OFFSET_KIND) :: upperplus
    !! upper bound for parallelization of ph-ph matrix elements
    INTEGER(KIND = MPI_OFFSET_KIND) :: lowerminus
    !! lower bound for parallelization of ph-ph matrix elements
    INTEGER(KIND = MPI_OFFSET_KIND) :: upperminus
    !! upper bound for parallelization of ph-ph matrix elements
    INTEGER :: lower_qbnd, upper_qbnd
    !! lower/upper bound for q parallelization
    REAL(KIND = DP) :: wqnu1
    !! frequency of 1st phonon
    REAL(KIND = DP) :: wqnu2 
    !! frequency of 2nd phonon
    REAL(KIND = DP) :: wqnu3p
    !! frequency of 3rd phonon
    REAL(KIND = DP) :: wqnu3m
    !! frequency of 3rd phonon
    REAL(KIND = DP) :: delta
    !! delta function
    REAL(KIND = DP) :: sigmap
    !! smearing parameter
    REAL(KIND = DP) :: sigmam 
    !! smearing parameter
    REAL(KIND = DP) :: q1mq2(3)
    !! q3 = q1 - q2
    REAL(KIND = DP) :: q1pq2(3)
    !! q3 = q1 + q2
    REAL(KIND = DP) :: sigma1max 
    !! smearing parameter
    REAL(KIND = DP) :: sigma1min
    !! smearing parameter
    REAL(KIND = DP) :: sigma2max
    !! smearing parameter
    REAL(KIND = DP) :: sigma2min
    !! smearing parameter
    REAL(KIND = DP) :: psi2
    !! ph-ph matrix element
    CHARACTER(LEN = 256) :: filename
    !! filename for mpi_write
    CHARACTER(LEN = 64) :: idim_ch
    !! Name of the files 
    REAL(KIND = DP) :: topsm1
    !! from ps^-1 to the inverse of Rydberg unit of time 
    !
    topsm1 = (ryd2ev * ps2sec) / hbar
    dt_in_ps = dph * tstep / 1000.d0 
    !! tstep in fs, so devided by 1000 to ps
    !
    CALL start_clock('3-phonon')
    !
    WRITE(stdout, '(5x, a/)') "Phonon-phonon interaction is included, calculate ph-ph coupling matrix elements."
    WRITE(stdout, '(5x, a)') REPEAT('=',67)    
    IF (double_grid) WRITE(stdout, '(/5x, a)') 'Warning : double grid for phonon-phonon &
                                               interaction is not available, use the   & 
                                               grid nqf1, nqf2 and nqf3 to calculate   &
                                               phonon_phonon interaction' 
    !
    IF (mpime == ionode_id) THEN
        INQUIRE(FILE = 'sBTE_config', EXIST = exst)
        IF (exst) THEN
          OPEN(2, file='sBTE_config', status='old',IOSTAT = ios)
          IF (ios /= 0) CALL errore('psi2_calc', 'error opening file sBTE_config',1)
          READ(2, *) lfactor
          DO ii = 1, 3
            READ(2, *) lattv(:,ii)
          ENDDO
          CLOSE(2)
          lattv = lfactor * lattv
        ELSE
          lattv = at * alat  
        ENDIF      
    ENDIF
    !
    !! Now lattice vector is in atomic unit (bohr)
    CALL mp_bcast(lattv, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    ! 
    DO ii = 1, 3
      jj = MOD(ii, 3) + 1
      kk = MOD(jj, 3) + 1
      CALL cross_product(lattv(:,jj), lattv(:,kk), rlattv(:,ii))
    END DO
    !
    Vcell = ABS(DOT_PRODUCT(lattv(:,1), rlattv(:,1)))
    rlattv = rlattv * twopi / Vcell
    scalebroad = facbroad   
    !  
    WRITE(stdout, '(/5x,a)') "Start running psi2_calc"
    CALL read3fc()
    WRITE(stdout, '(/5x,a)') '3rd order force constants have been read from FORCE_CONSTANTS_3RD'
    ! WRITE(stdout, *) 'ityp = ',ityp
    ! WRITE(stdout, *) 'amass = ', amass
    ! WRITE(stdout, *) 'nmodes=',nmodes,"nqtotf=",nqtotf
    CALL fkbounds(nqtotf, lower_qbnd, upper_qbnd)
    ! 
    Ntot_p = 0
    Ntot_m = 0
    N_plus = 0 
    N_minus = 0
!    N_plus_all = 0
!    N_minus_all = 0
    !
    !! The code is modified from ShengBTE/processes.f90
    !! Here we use atomic unit
    IF (adapt_smearing) THEN
      WRITE(stdout,'(/5x,a)') 'Use adaptive smearing in 3-phonon process'
      sigma1min = 1000.0_dp
      sigma2min = 1000.0_dp
      sigma1max = 0.0_dp
      sigma2max = 0.0_dp
    ENDIF
    !
    phwmin = phwmin / ryd2mev ! to Ryd 
    ! 
    WRITE(stdout, '(/5x, a)') REPEAT('=',67) 
    WRITE(stdout, '(/5x, a/)') 'Step 1: Now count the number of allowed 3-phonon scattering processes'
    DO ii = 1, 2
      !! ii = 1 : count and allocate
      !! ii = 2 : save the index while counting
      ! IF (ii == 2 ) THEN
        ind1 = 0
        ind2 = 0
      ! ENDIF !
      DO iq1 = lower_qbnd, upper_qbnd
        DO nu1 = 1, nmodes
          wqnu1 = wf(nu1, iq1)
          IF (wqnu1 < phwmin) CYCLE
          DO iq2 = 1, nqtotf
            DO nu2 = 1, nmodes
              wqnu2 = wf(nu2, iq2)
              IF (wqnu2 < phwmin) CYCLE
              q1pq2(:) = xqf(:, iq1) + xqf(:, iq2)
              iqx = MODULO(INT(q1pq2(1) * nqf1), nqf1) + 1
              iqy = MODULO(INT(q1pq2(2) * nqf2), nqf2) + 1
              iqz = MODULO(INT(q1pq2(3) * nqf3), nqf3) + 1
              iq3p = (iqx - 1) * nqf2 * nqf3 + &
                     (iqy - 1) * nqf3 + iqz
              q1mq2(:) = xqf(:, iq1) - xqf(:, iq2)
              iqx = MODULO(INT(q1mq2(1) * nqf1), nqf1) + 1
              iqy = MODULO(INT(q1mq2(2) * nqf2), nqf2) + 1
              iqz = MODULO(INT(q1mq2(3) * nqf3), nqf3) + 1
              iq3m = (iqx - 1) * nqf2 * nqf3 + &
                     (iqy - 1) * nqf3 + iqz
              DO nu3 = 1, nmodes
                wqnu3p = wf(nu3, iq3p)
                wqnu3m = wf(nu3, iq3m)
                IF (adapt_smearing) THEN
                  sigmap = scalebroad * base_sigma(&
                            vnuq_all(:, nu2, iq2) -&
                            vnuq_all(:, nu3, iq3p))
                  sigmam = scalebroad * base_sigma(&
                            vnuq_all(:, nu2, iq2)-&
                            vnuq_all(:, nu3, iq3m))
                  IF (ii == 1) THEN
                    IF (sigmap < sigma1min) sigma1min = sigmap
                    IF (sigmap > sigma1max) sigma1max = sigmap
                    IF (sigmam < sigma2min) sigma2min = sigmam
                    IF (sigmam > sigma2max) sigma2max = sigmam
                  ENDIF ! ii < 2
                ElSE
                  sigmap = degaussq
                  sigmam = degaussq
                ENDIF ! adapt_smearing
                IF (ABS(wqnu1 + wqnu2 - wqnu3p) <= (2.d0 * sigmap).AND. &
                  wqnu3p > phwmin)  THEN
                  IF (ii == 1) THEN
                    N_plus(my_pool_id + 1)  = N_plus(my_pool_id + 1) + 1
                  ELSE ! ii == 2
                    ind1 = ind1 + 1
                    indxp_tmp(ind1, 1) = iq1
                    indxp_tmp(ind1, 2) = nu1
                    indxp_tmp(ind1, 3) = iq2
                    indxp_tmp(ind1, 4) = nu2
                    indxp_tmp(ind1, 5) = iq3p
                    indxp_tmp(ind1, 6) = nu3
                  ENDIF !! ii == 1 or 2
                ENDIF  !  |wqnu1 + wqnu2 - wqnu3p| < 2 * sigmap                             
                !
                IF (ABS(wqnu1 - wqnu2 - wqnu3m) <= (2.d0 * sigmam).AND. &
                  wqnu3m > phwmin) THEN 
                  IF (ii == 1) THEN
                    N_minus(my_pool_id + 1) = N_minus(my_pool_id + 1) + 1
                  ELSE ! ii == 2
                    ind2 = ind2 + 1
                    indxm_tmp(ind2, 1) = iq1
                    indxm_tmp(ind2, 2) = nu1
                    indxm_tmp(ind2, 3) = iq2
                    indxm_tmp(ind2, 4) = nu2
                    indxm_tmp(ind2, 5) = iq3m
                    indxm_tmp(ind2, 6) = nu3                  
                  ENDIF ! ii
                ENDIF ! |wqnu1 - wqnu2 - wqnu3m| < 2 * sigmam  
              ENDDO ! nu3
            ENDDO ! nu2
          ENDDO ! iq2
        ENDDO ! nu1
      ENDDO ! iq1
      !
      CALL mp_barrier(inter_pool_comm)
      IF (ii == 1) THEN
        !! N_plus and N_minus are kind = 8, can be accepted by mpi_allreduce
        CALL MPI_ALLREDUCE(MPI_IN_PLACE, N_plus,  npool, MPI_INTEGER8, MPI_SUM, inter_pool_comm, ierr)
        CALL MPI_ALLREDUCE(MPI_IN_PLACE, N_minus, npool, MPI_INTEGER8, MPI_SUM, inter_pool_comm, ierr)
        !! N_plus and N_minus are kind = 8, can not be accepted by mp_sum
        ! CALL mp_sum(N_plus, inter_pool_comm)
        ! CALL mp_sum(N_minus, inter_pool_comm)
!        CALL MPI_ALLREDUCE(N_plus,  N_plus_all,  npool, MPI_INTEGER8, MPI_SUM, inter_pool_comm, ierr)
!        CALL MPI_ALLREDUCE(N_minus, N_minus_all, npool, MPI_INTEGER8, MPI_SUM, inter_pool_comm, ierr)
!        N_plus  = N_plus_all
!        N_minus = N_minus_all
!         WRITE(stdout,*) 'N_plus = ', N_plus(:)
!         WRITE(stdout,*) 'N_minus = ', N_minus(:)
        CALL mp_barrier(inter_pool_comm)
        IF (adapt_smearing) THEN
          CALL mp_max(sigma1max, inter_pool_comm)
          CALL mp_min(sigma1min, inter_pool_comm)
          CALL mp_max(sigma2max, inter_pool_comm)
          CALL mp_min(sigma2min, inter_pool_comm)
          WRITE(stdout,'(5x,a,ES20.10,a)') 'sigma1max = ', sigma1max * rydcm1, ' cm^-1'
          WRITE(stdout,'(5x,a,ES20.10,a)') 'sigma2max = ', sigma2max * rydcm1, ' cm^-1'
          WRITE(stdout,'(5x,a,ES20.10,a)') 'sigma1min = ', sigma1min * rydcm1, ' cm^-1'
          WRITE(stdout,'(5x,a,ES20.10,a)') 'sigma2min = ', sigma2min * rydcm1, ' cm^-1'
        ENDIF
        !
        Ntot_p = SUM(N_plus)
        Ntot_m = SUM(N_minus)
        CALL mp_barrier(inter_pool_comm)
        !! Ntot_p_ and Ntot_p_ are MPI_OFFSET_KIND
        Ntot_p_ = INT(Ntot_p, KIND = MPI_OFFSET_KIND)
        Ntot_m_ = INT(Ntot_m, KIND = MPI_OFFSET_KIND)
        WRITE(stdout,'(5x,a,i20)')  "Ntot_p = ", Ntot_p
        WRITE(stdout,'(5x,a,i20)')  "Ntot_m = ", Ntot_m 
        WRITE(stdout,'(/5x,a/)') 'Now find the indexes of three phonon scattering '
        !WRITE(stdout,*) 'Check the memory'
        ALLOCATE(indxp_tmp(N_plus(my_pool_id + 1), 6), STAT = ierr)
        IF (ierr /= 0) CALL errore('psi2_calc', 'Error allocating indxp_tmp', 1)
        ALLOCATE(indxm_tmp(N_minus(my_pool_id + 1),6), STAT = ierr)
        IF (ierr /= 0) CALL errore('psi2_calc', 'Error allocating indxm_tmp', 1)
        WRITE(stdout,'(5x,a)') 'Step 2: indexing according to the 1st q point'
        !
        indxp_tmp = 0
        indxm_tmp = 0
        ! ind1 = 0
        ! ind2 = 0
        !
        CALL mp_barrier(inter_pool_comm)
      ENDIF ! ii == 1
    ENDDO ! ii
    !
    WRITE(stdout, '(/5x,a/)') 'Step 3: distribute 3-phonon indxes'
    !
    CALL fkbounds2(Ntot_p_,  lowerplus, upperplus)
    CALL fkbounds2(Ntot_m_, lowerminus, upperminus)
    !   
    nind_p = upperplus - lowerplus + 1
    nind_m = upperminus - lowerminus + 1
    !
    ! WRITE(stdout, *) 'nind_p =', nind_p
    ! WRITE(stdout, *) 'nind_m =', nind_m
    ALLOCATE(ind_p(nind_p, 6), STAT = ierr)
    IF (ierr /= 0) CALL errore('psi2_calc', 'Error allocating ind_p', 1)
    ALLOCATE(ind_m(nind_m, 6), STAT = ierr)
    IF (ierr /= 0) CALL errore('psi2_calc', 'Error allocating ind_m', 1)
    !
    ind_p = 0
    ind_m = 0
    !
    DO ii = 1, 2 
      !! ii == 1 : redistribute indxp_tmp
      !! ii == 2 : redistribute indxm_tmp
      IF (ii == 1) THEN
        nindx =  INT(N_plus(my_pool_id + 1), KIND = MPI_OFFSET_KIND)
        nindxi = INT(nind_p, KIND = MPI_OFFSET_KIND)
      ELSE
        nindx =  INT(N_minus(my_pool_id + 1), KIND = MPI_OFFSET_KIND)
        nindxi = INT(nind_m, KIND = MPI_OFFSET_KIND)
      ENDIF !! ii == 1 or 2
      CALL mp_barrier(inter_pool_comm)
      !
      displ1 = 0_MPI_OFFSET_KIND
      displ2 = 0_MPI_OFFSET_KIND
      !
      IF (ii == 1) THEN
        IF (my_pool_id > 0) THEN
          DO ipool = 1, my_pool_id
            displ1 = displ1 + INT(N_plus(ipool), KIND = MPI_OFFSET_KIND)
          ENDDO
        ENDIF
        !
        displ1 = displ1 * 4_MPI_OFFSET_KIND
        displ2 = (lowerplus - 1_MPI_OFFSET_KIND) * 4_MPI_OFFSET_KIND
        !
        WRITE(stdout,'(7x,a)') 'Now distribute indx_plus  evenly across the pools'
      ELSE ! ii == 2
        IF (my_pool_id > 0) THEN
          DO ipool = 1, my_pool_id
            displ1 = displ1 + INT(N_minus(ipool), KIND = MPI_OFFSET_KIND)
          ENDDO
        ENDIF
        !
        displ1 = displ1 * 4_MPI_OFFSET_KIND
        displ2 = (lowerminus - 1_MPI_OFFSET_KIND) * 4_MPI_OFFSET_KIND    
        !
        WRITE(stdout,'(7x,a)') 'Now distribute indx_minus evenly across the pools'  
      ENDIF
      CALL mp_barrier(inter_pool_comm)
      !
      IF (ii == 1) THEN    
        DO idim = 1, 6
          WRITE(idim_ch, "(I0)") idim
          filename = 'indph_p'// TRIM(idim_ch)
          ! open file
          CALL MPI_FILE_OPEN(inter_pool_comm, filename, MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, iunphindx, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_OPEN', 1)
          ! Write data to file
          CALL MPI_FILE_WRITE_AT(iunphindx, displ1, indxp_tmp(:,idim), nindx, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_WRITE_AT', 1)
          ! Close file
          CALL MPI_FILE_CLOSE(iunphindx, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_CLOSE', 1)
          ! Wait for all processes to finish writing
          CALL mp_barrier(inter_pool_comm)
          ! Open file for reading
          CALL MPI_FILE_OPEN(inter_pool_comm, filename, &
               MPI_MODE_RDONLY + MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, iunphindx, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_OPEN', 1)
          ! Read data from file
          CALL MPI_FILE_READ_AT(iunphindx, displ2, ind_p(:,idim), nindxi, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_READ_AT', 1)
          ! Close file
          CALL MPI_FILE_CLOSE(iunphindx, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_CLOSE', 1)
          !
          CALL mp_barrier(inter_pool_comm)
        ENDDO  ! idim
      ELSE !! ii == 2
        DO idim = 1, 6
          WRITE(idim_ch, "(I0)") idim
          filename = 'indph_m'// TRIM(idim_ch)
          ! open file
          CALL MPI_FILE_OPEN(inter_pool_comm, filename, &
               MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, iunphindx, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_OPEN', 1)
          ! Write data to file
          CALL MPI_FILE_WRITE_AT(iunphindx, displ1, indxm_tmp(:,idim), nindx, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_WRITE_AT', 1)
          ! Close file
          CALL MPI_FILE_CLOSE(iunphindx, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_CLOSE', 1)
          ! Wait for all processes to finish writing
          CALL mp_barrier(inter_pool_comm)
          ! Open file for reading
          CALL MPI_FILE_OPEN(inter_pool_comm, filename, MPI_MODE_RDONLY + MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, iunphindx, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_OPEN', 1)
          ! Read data from file
          CALL MPI_FILE_READ_AT(iunphindx, displ2, ind_m(:,idim), nindxi, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_READ_AT', 1)
          ! Close file
          CALL MPI_FILE_CLOSE(iunphindx, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_CLOSE', 1)
          !
          CALL mp_barrier(inter_pool_comm)
        ENDDO  ! idim
      ENDIF
    ENDDO !!   
    !
    ! WRITE(stdout, *) 'ind'
    ! DO idim = 1, 6
    ! WRITE(stdout, *) ind_m(10, idim)
    ! WRITE(stdout, *) ind_p(10, idim)
    ! ENDDO
    ! WRITE(stdout, *) 'indxp'
    ! DO idim = 1, 6
    ! WRITE(stdout, *) indxm_tmp(10, idim)
    ! WRITE(stdout, *) indxp_tmp(10, idim)
    ! ENDDO
    DEALLOCATE(indxp_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('psi2_calc', 'Error deallocating indxp_tmp', 1)
    DEALLOCATE(indxm_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('psi2_calc', 'Error deallocating indxm_tmp', 1)
    !
    ALLOCATE(psi2_p(nind_p), STAT = ierr)
    IF (ierr /= 0) CALL errore('psi2_calc', 'Error allocating psi2_p', 1)
    ALLOCATE(psi2_m(nind_m), STAT = ierr)
    IF (ierr /= 0) CALL errore('psi2_calc', 'Error allocating psi2_m', 1)
    psi2_p  = 0.d0
    psi2_m = 0.d0
    WRITE(stdout,'(/5x, a/)') 'Step 4: Calculate ph-ph matrix elements.'
    DO ind1 = 1, nind_p
      iq1  = ind_p(ind1, 1)
      nu1  = ind_p(ind1, 2)
      iq2  = ind_p(ind1, 3)
      nu2  = ind_p(ind1, 4)
      iq3p = ind_p(ind1, 5)
      nu3  = ind_p(ind1, 6)
      CALL psi_interp(nu1, nu2, nu3, iq1, iq2, iq3p, 1.0_dp, psi2)
      IF (adapt_smearing) THEN
        sigmap = scalebroad * base_sigma(          &
                  vnuq_all(:, nu2, iq2) -          &
                  vnuq_all(:, nu3, iq3p))
      ElSE
        sigmap = degaussq
      ENDIF    
      wqnu1  = wf(nu1, iq1)
      wqnu2  = wf(nu2, iq2)
      wqnu3p = wf(nu3, iq3p)
      delta = EXP(- (wqnu1 + wqnu2 - wqnu3p) ** 2  &
           / (sigmap ** 2)) / sigmap / SQRT(pi)
      !
      psi2_p(ind1) = pi / 4.d0 * psi2 * &
                        delta / nqtotf 
      psi2_p(ind1) = psi2_p(ind1) * topsm1
      !! Now the unit of psi2_p is ps^-1
    ENDDO
    DO ind2 = 1, nind_m
      iq1  = ind_m(ind2, 1)
      nu1  = ind_m(ind2, 2)
      iq2  = ind_m(ind2, 3)
      nu2  = ind_m(ind2, 4)
      iq3m = ind_m(ind2, 5)
      nu3  = ind_m(ind2, 6)
      CALL psi_interp(nu1, nu2, nu3, iq1, iq2, iq3m, -1.0_dp, psi2)
      IF (adapt_smearing) THEN
        sigmam=scalebroad * base_sigma(&
               vnuq_all(:, nu2, iq2)-&
               vnuq_all(:, nu3, iq3m))
      ElSE
        sigmam = degaussq
      ENDIF 
      !
      wqnu1  = wf(nu1, iq1)
      wqnu2  = wf(nu2, iq2)
      wqnu3m = wf(nu3, iq3m)
      delta = EXP(-(wqnu1 - wqnu2 - wqnu3m) **2    &
            / (sigmam ** 2)) / sigmam / SQRT(pi)
      psi2_m(ind2) = pi / 4.d0 * psi2 *  &
                          delta / nqtotf
      psi2_m(ind2) = psi2_m(ind2) * topsm1
      !! Now the unit of psi2_m is ps^-1
    ENDDO
    !
    ! WRITE(stdout, *) psi2_m(1:10)
    ! WRITE(stdout, *) psi2_p(1:10)
    WRITE(stdout, '(/5x, a)') REPEAT('=',67)
    WRITE(stdout, '(/5x, a/)') '3-phonon scattering matrix elements have been calculated'
    WRITE(stdout, '(5x, a)') REPEAT('=',67)
    CALL mp_barrier(inter_pool_comm)
    !
    RETURN
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE psi2_calc
    !--------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------- 
    SUBROUTINE iphph_calc()
    !-------------------------------------------------------------------------- 
    !! Calculate phonon-phonon collision integral. More details are given by 
    !! the Eq. (S3) of J. Phys. Chem. Lett. 12, 1734.
    !-------------------------------------------------------------------------- 
    !     
    USE kinds,            ONLY : DP
    USE elph2,            ONLY : nqtotf
    USE mp,               ONLY : mp_barrier, mp_sum, mp_bcast
    USE mp_global,        ONLY : inter_pool_comm, npool, my_pool_id
    USE modes,            ONLY : nmodes
    USE tdbecom,          ONLY : nphon_pre, iph_ph, a_tab, b_tab,c_tab, &
                                 psi2_p, psi2_m, ind_p, ind_m, nstg,    &
                                 dt_in_ps, nind_p, nind_m                       
    !     
    IMPLICIT NONE
    ! 
    INTEGER :: iq1, iq2, iq3
    !! q-indexes of the 3 phonons involved in the ph-ph interaction
    INTEGER :: nu1, nu2, nu3
    !! phonon branch indexes of the 3 phonons involved in the ph-ph interaction
    INTEGER :: istg, pp
    !! Counter on the stages of runge-kutta solver
    INTEGER(KIND = 8) :: indx
    !! Counter on ph-ph coupling matrix element
    REAL(KIND = DP) :: nph1, nph2, nph3
    !! Nonequilibrium phonon distribution of the three phonons involved
    REAL(KIND = DP) :: facp
    !! factor arise from products of phonon distribution
    ! 
    iph_ph = 0.d0
    DO istg = 1, nstg
      DO indx = 1, nind_p
        iq1 = ind_p(indx, 1)
        nu1 = ind_p(indx, 2)
        iq2 = ind_p(indx, 3)
        nu2 = ind_p(indx, 4)
        iq3 = ind_p(indx, 5)
        nu3 = ind_p(indx, 6)
        nph1 = nphon_pre(nu1, iq1)
        nph2 = nphon_pre(nu2, iq2)
        nph3 = nphon_pre(nu3, iq3)
        IF (istg > 1) THEN
          DO pp = 1, istg -1
            nph1 = nph1 + a_tab(istg, pp) * iph_ph(nu1, iq1, pp) * dt_in_ps
            nph2 = nph2 + a_tab(istg, pp) * iph_ph(nu2, iq2, pp) * dt_in_ps
            nph3 = nph3 + a_tab(istg, pp) * iph_ph(nu3, iq3, pp) * dt_in_ps
          ENDDO ! pp
        ENDIF ! istg
!! In Eq. (S3) facp = (1.d0 + nph1) * (1.d0 + nph2) * &
!!      nph3 - nph1 * nph2* (1.d0 + nph3)
!!      It is simplified to the expression below
        facp = - nph1 * (nph2 - nph3) + nph3 * (nph2 + 1.d0)
        iph_ph(nu1, iq1, istg) = iph_ph(nu1, iq1, istg) + psi2_p(indx) * facp
      ENDDO
      !
      DO indx = 1, nind_m
        iq1 = ind_m(indx, 1)
        nu1 = ind_m(indx, 2)
        iq2 = ind_m(indx, 3)
        nu2 = ind_m(indx, 4)
        iq3 = ind_m(indx, 5)
        nu3 = ind_m(indx, 6)
        nph1 = nphon_pre(nu1, iq1)
        nph2 = nphon_pre(nu2, iq2)
        nph3 = nphon_pre(nu3, iq3)
        IF (istg > 1) THEN
          DO pp = 1, istg -1
            nph1 = nph1 + a_tab(istg, pp) * iph_ph(nu1, iq1, pp) * dt_in_ps
            nph2 = nph2 + a_tab(istg, pp) * iph_ph(nu2, iq2, pp) * dt_in_ps
            nph3 = nph3 + a_tab(istg, pp) * iph_ph(nu3, iq3, pp) * dt_in_ps
          ENDDO ! pp
        ENDIF ! istg
!! In Eq. (S3) facp = (1.d0 + nph1) * nph2 * &
!!      nph3 - nph1 * (1.d0 + nph2) * (1.d0 + nph3)
!!      It is simplified to the expression below
        facp = nph2 * nph3 - nph1 * (nph2 + nph3 + 1.d0)
        iph_ph(nu1, iq1, istg) = iph_ph(nu1, iq1, istg) + psi2_m(indx) * facp * 0.5d0
      ENDDO
      CALL mp_sum(iph_ph(:, :, istg), inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
    ENDDO ! istg
    !
    RETURN
    !
    !-------------------------------------------------------------------------- 
    END SUBROUTINE iphph_calc
    !-------------------------------------------------------------------------- 
    !  
    !-------------------------------------------------------------------------- 
    SUBROUTINE inv_tau_rta(temp_min, temp_max, ntemps)
    !-------------------------------------------------------------------------- 
    !! calculate inverse phonon lifetime (linewidth) under RTA, the details of the
    !! expression can be found from many literatures, for example, Eq.(6-8) of 
    !! Comp. Phys. Commun. 185, 17471758 (2014).
    !--------------------------------------------------------------------------     
    ! 
    USE kinds,            ONLY : DP
    USE mp_world,         ONLY : mpime
    USE io_global,        ONLY : stdout, ionode, ionode_id
    USE elph2,            ONLY : nqtotf, wf
    USE mp,               ONLY : mp_barrier, mp_sum, mp_bcast
    USE mp_global,        ONLY : inter_pool_comm, npool, my_pool_id
    USE modes,            ONLY : nmodes
    USE tdbecom,          ONLY : ind_p, ind_m, nind_p, nind_m,      &
                                 psi2_p, psi2_m
    USE constants_epw,    ONLY : kelvin2Ry
    !     
    IMPLICIT NONE
    !
    REAL(KIND=DP),INTENT(in)  :: temp_min, temp_max
    !! lower and upper bound for temperature, in Kelvin
    INTEGER ,INTENT(in) :: ntemps
    !! Number of point in the interval [tmin, tmax]
    REAL(KIND = DP)  :: inv_tau_ph(ntemps, nqtotf,nmodes)
    !! inverse ph-ph lifetime
    REAL(KIND = DP) :: nph2, nph3
    !! Thermal phonon population 
    REAL(KIND = DP) :: temp_in_ry  
    !! Temperature in kelvin (k_B * T)
    INTEGER :: iq1, iq2, iq3
    !! q index for the 3 phonons involved
    INTEGER :: nu1, nu2, nu3
    !! Phonon branch indexes for the 3 phonons involved
    INTEGER(KIND = 8) :: indx
    !! Counter on ph-ph matrix elements
    INTEGER :: itemp
    !! Counter on temperature
    character(len=4) :: temp_in_ch
    !! temperature in string
    character(len=32) ::fname
    !! file name 
    !
    inv_tau_ph = 0.d0
    DO itemp = 1, ntemps
      temp_in_ry = (temp_min + (temp_max - temp_min) / (ntemps - 1)  &
                 * (itemp - 1)) * kelvin2Ry ! 6.33363E-6 is kelvin2Ry
      DO indx = 1, nind_p
        iq1 = ind_p(indx, 1)
        nu1 = ind_p(indx, 2)
        iq2 = ind_p(indx, 3)
        nu2 = ind_p(indx, 4)
        iq3 = ind_p(indx, 5)
        nu3 = ind_p(indx, 6)
        nph2 = 1.0 / (EXP(wf(nu2, iq2) / temp_in_ry) - 1.0)
        ! Bose-Einstein distribution
        nph3 = 1.0 / (EXP(wf(nu3, iq3) / temp_in_ry) - 1.0)
        ! Bose-Einstein distribution
        inv_tau_ph(itemp, iq1, nu1) = inv_tau_ph(itemp, iq1, nu1) +  &
                                      psi2_p(indx) * (nph2 - nph3)
      ENDDO
      !
      DO indx = 1, nind_m
        iq1 = ind_m(indx, 1)
        nu1 = ind_m(indx, 2)
        iq2 = ind_m(indx, 3)
        nu2 = ind_m(indx, 4)
        iq3 = ind_m(indx, 5)
        nu3 = ind_m(indx, 6)
        nph2 = 1.0 / (EXP(wf(nu2, iq2) / temp_in_ry) - 1.0)
        ! Bose-Einstein distribution
        nph3 = 1.0 / (EXP(wf(nu3, iq3) / temp_in_ry) - 1.0)
        ! Bose-Einstein distribution
        inv_tau_ph(itemp, iq1, nu1) = inv_tau_ph(itemp, iq1, nu1) +  &
                  psi2_m(indx) * (nph2 + nph3 + 1.0) / 2.0
      ENDDO
      !
    ENDDO
    CALL mp_sum(inv_tau_ph, inter_pool_comm)
    WRITE(stdout, '(/5x,a/)') 'print the inverse life time of phonon due to phonon-phonon interaction.'
    IF (mpime == ionode_id) THEN
      DO itemp = 1, ntemps
        temp_in_ry = (temp_min + (temp_max - temp_min)/(ntemps - 1) * (itemp - 1))
        WRITE(temp_in_ch, "(I0)") NINT(temp_in_ry)
        fname = "inv_tau_rta_" // TRIM(adjustl(temp_in_ch)) // "K"
        OPEN(unit = 222, FILE = fname)
        DO iq1 = 1, nqtotf
          DO nu1 = 1, nmodes
            WRITE(222, '(2E22.14)') wf(nu1, iq1), inv_tau_ph(itemp, iq1, nu1)
          ENDDO
        ENDDO
        CLOSE(222)
      ENDDO
    ENDIF
    !
    RETURN
    !
    !--------------------------------------------------------------------------------
    END SUBROUTINE inv_tau_rta
    !--------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------- 
    SUBROUTINE phphself(temp_min, temp_max, ntemps)
    !-------------------------------------------------------------------------- 
    !! calculate phonon self energy, the details of the
    !! expression can be found from many literatures, for example, Eq.(2) of 
    !! Phys. Rev. B 68, 220509 (2003).
    !--------------------------------------------------------------------------     
    ! 
    USE kinds,            ONLY : DP
    USE mp_world,         ONLY : mpime
    USE io_global,        ONLY : stdout, ionode, ionode_id
    USE elph2,            ONLY : nqtotf, wf,  xqf
    USE mp,               ONLY : mp_barrier, mp_sum, mp_bcast
    USE mp_global,        ONLY : inter_pool_comm, npool, my_pool_id
    USE modes,            ONLY : nmodes
    USE constants_epw,    ONLY : kelvin2Ry, ryd2mev, twopi 
    USE division,         ONLY : fkbounds
    USE tdbecom,          ONLY : lattv, rlattv, Vcell, lfactor, uf_all
    USE cell_base,        ONLY : at, bg, alat
    USE epwcom,           ONLY : nqf1, nqf2, nqf3, phwmin, phphself_r
    USE elph2,            ONLY : zstar, epsi
    USE ions_base,        ONLY : ityp, tau
    !     
    IMPLICIT NONE
    !
    REAL(KIND=DP),INTENT(in)  :: temp_min, temp_max
    !! lower and upper bound for temperature, in Kelvin
    INTEGER ,INTENT(in) :: ntemps
    !! Number of point in the interval [tmin, tmax]
    REAL(KIND = DP)  :: pi_anh(ntemps, nqtotf, nmodes)
    !! phonon-phonon self energy \Pi_{anh}
    REAL(KIND = DP) :: wqnu1, wqnu2, wqnu3p, wqnu3m
    !! frequency of 3 phonons
    REAL(KIND = DP) :: nph2, nph3m, nph3p
    !! Thermal phonon population 
    REAL(KIND = DP) :: psi2p, psi2m
    !! Square of phonon-phonon coupling matrix elements
    REAL(KIND = DP) :: temp_in_ry  
    !! Temperature in kelvin (k_B * T)
    REAL(KIND = DP) :: eta
    !! small value regularizing the zero denominator
    REAL(KIND = DP) :: q1mq2(3)
    !! q3 = q1 - q2
    REAL(KIND = DP) :: q1pq2(3)
    !! q3 = q1 + q2
    INTEGER :: ierr
    !! Error status
    INTEGER :: ios 
    !! Error status
    INTEGER :: ii, jj, kk
    !! Counter on 3 spacial directions
    INTEGER :: iqx, iqy, iqz
    !! q index on the 3 directions
    INTEGER :: iq1, iq2, iq3p, iq3m
    !! q index for the 3 phonons involved
    INTEGER :: nu1, nu2, nu3
    !! Phonon branch indexes for the 3 phonons involved
    INTEGER :: lower_qbnd
    !! lower q index
    INTEGER :: upper_qbnd
    !! upper q index
    INTEGER :: itemp
    !! Counter on temperature
    LOGICAL :: exst
    !! If the file containing lattice constants exists?
    character(len=4) :: temp_in_ch
    !! temperature in string
    character(len=32) ::fname
    !! file name 
    WRITE(stdout, '(5x, a/)') "Evaluate ph-ph self energy."
    WRITE(stdout, '(5x, a)') REPEAT('=',67)    
    !
    IF (mpime == ionode_id) THEN
      INQUIRE(FILE = 'sBTE_config', EXIST = exst)
      IF (exst) THEN
        OPEN(2, file='sBTE_config', status='old',IOSTAT = ios)
        IF (ios /= 0) CALL errore('psi2_calc', 'error opening file sBTE_config',1)
        READ(2, *) lfactor
        DO ii = 1, 3
          READ(2, *) lattv(:,ii)
        ENDDO
        CLOSE(2)
        lattv = lfactor * lattv
      ELSE
        lattv = at * alat  
      ENDIF      
    ENDIF
    !
    !! Now lattice vector is in atomic unit (bohr)
    CALL mp_bcast(lattv, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    ! 
    DO ii = 1, 3
      jj = MOD(ii, 3) + 1
      kk = MOD(jj, 3) + 1
      CALL cross_product(lattv(:,jj), lattv(:,kk), rlattv(:,ii))
    END DO
    !
    Vcell = ABS(DOT_PRODUCT(lattv(:,1), rlattv(:,1)))
    rlattv = rlattv * twopi / Vcell 
    !  
    CALL read3fc()
    WRITE(stdout, '(/5x,a)') '3rd order force constants have been read from FORCE_CONSTANTS_3RD'
    ! WRITE(stdout, *) 'ityp = ',ityp
    ! WRITE(stdout, *) 'amass = ', amass
    ! WRITE(stdout, *) 'nmodes=',nmodes,"nqtotf=",nqtotf
    !
    pi_anh = 0.d0 
    eta  = 1.0E-6
    !!
    CALL fkbounds(nqtotf, lower_qbnd, upper_qbnd)
    DO iq1 = lower_qbnd, upper_qbnd
      DO nu1 = 1, nmodes 
        wqnu1 = wf(nu1, iq1)
        IF (wqnu1 < phwmin) CYCLE
        DO iq2 = 1, nqtotf
          DO nu2 = 1, nmodes 
            wqnu2 = wf(nu2, iq2)
            IF (wqnu2 < phwmin) CYCLE
            q1pq2(:) = xqf(:, iq1) + xqf(:, iq2)
            iqx = MODULO(INT(q1pq2(1) * nqf1), nqf1) + 1
            iqy = MODULO(INT(q1pq2(2) * nqf2), nqf2) + 1
            iqz = MODULO(INT(q1pq2(3) * nqf3), nqf3) + 1
            iq3p = (iqx - 1) * nqf2 * nqf3 + &
                    (iqy - 1) * nqf3 + iqz
            !
            q1mq2(:) = xqf(:, iq1) - xqf(:, iq2)
            iqx = MODULO(INT(q1mq2(1) * nqf1), nqf1) + 1
            iqy = MODULO(INT(q1mq2(2) * nqf2), nqf2) + 1
            iqz = MODULO(INT(q1mq2(3) * nqf3), nqf3) + 1
            iq3m = (iqx - 1) * nqf2 * nqf3 + &
                    (iqy - 1) * nqf3 + iqz
            !
            DO nu3 = 1, nmodes
              wqnu3p = wf(nu2, iq3p)
              IF (wqnu3p < phwmin) CYCLE
              wqnu3m = wf(nu2, iq3m)
              IF (wqnu3m < phwmin) CYCLE
              CALL psi_interp(nu1, nu2, nu3, iq1, iq2, iq3p,  1.0_dp, psi2p)
              !
              CALL psi_interp(nu1, nu2, nu3, iq1, iq2, iq3m, -1.0_dp, psi2m)
              !
              psi2p = psi2p / 8.d0 / nqtotf 
              psi2m = psi2m / 8.d0 / nqtotf 
              !! The unit of psi2p and psi2m is [Ry]^2
              !
              DO itemp = 1, ntemps
                temp_in_ry = (temp_min + (temp_max - temp_min) / (ntemps - 1)  &
                             * (itemp - 1)) * kelvin2Ry ! 6.33363E-6 is kelvin2Ry
                nph2  = 1.0 / (EXP(wqnu2 / temp_in_ry) - 1.0)
                ! Bose-Einstein distribution
                nph3m = 1.0 / (EXP(wqnu3m / temp_in_ry) - 1.0)
                ! Bose-Einstein distribution
                nph3p = 1.0 / (EXP(wqnu3p / temp_in_ry) - 1.0)
                ! Bose-Einstein distribution
                pi_anh(itemp, iq1, nu1) = pi_anh(itemp, iq1, nu1) - psi2m * (nph2 + nph3m + 1.d0) / 2.d0 * &
                 (( wqnu1 + wqnu2 + wqnu3m) / (( wqnu1 + wqnu2 + wqnu3m)**2 + eta**2) +  &
                  (-wqnu1 + wqnu2 + wqnu3m) / ((-wqnu1 + wqnu2 + wqnu3m)**2 + eta**2))
                !
                pi_anh(itemp, iq1, nu1) = pi_anh(itemp, iq1, nu1) - psi2p * (nph2 - nph3p)        / 2.d0 * &
                 (( wqnu1 - wqnu2 + wqnu3p) / (( wqnu1 - wqnu2 + wqnu3p)**2 + eta**2) +  &
                  (-wqnu1 - wqnu2 + wqnu3p) / ((-wqnu1 - wqnu2 + wqnu3p)**2 + eta**2))
              ENDDO
            ENDDO ! nu3
          ENDDO ! nu2 
        ENDDO ! iq2
      ENDDO ! nu1
    ENDDO ! iq1
    !
    CALL mp_sum(pi_anh, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    !! Start output
    WRITE(stdout, '(/5x,a/)') 'print the self energy of phonon due to phonon-phonon interaction.'
    IF (mpime == ionode_id) THEN
      DO itemp = 1, ntemps
        temp_in_ry = (temp_min + (temp_max - temp_min)/(ntemps - 1) * (itemp - 1))
        WRITE(temp_in_ch, "(I0)") NINT(temp_in_ry)
        fname = "pi_anh_" // TRIM(adjustl(temp_in_ch)) // "K"
        OPEN(unit = 222, FILE = fname)
        DO iq1 = 1, nqtotf
          DO nu1 = 1, nmodes
            WRITE(222, '(2E22.14)') wf(nu1, iq1)* ryd2mev, pi_anh(itemp, iq1, nu1) * ryd2mev
          ENDDO
        ENDDO
        CLOSE(222)
      ENDDO
    ENDIF
    !!! Clear up the arrays allocated
    IF(ALLOCATED(epsi))     DEALLOCATE(epsi)
    IF(ALLOCATED(zstar))    DEALLOCATE(zstar)
    IF(ALLOCATED(ityp))     DEALLOCATE(ityp)
    IF(ALLOCATED(tau))      DEALLOCATE(tau)
    !!!
    DEALLOCATE(wf, STAT = ierr)
    IF (ierr /= 0) CALL errore('phphself', 'Error deallocating wf', 1)
    DEALLOCATE(xqf, STAT = ierr)
    IF (ierr /= 0) CALL errore('phphself', 'Error deallocating xqf', 1)
    DEALLOCATE(uf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('phphself','Error deallocating uf_all',1)
    !
    RETURN
    !
    !--------------------------------------------------------------------------------
    END SUBROUTINE phphself
    !--------------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------------
    SUBROUTINE read3fc()
    !--------------------------------------------------------------------------------
    !! This subroutine is adapted from input.f90 of ShengBTE code 
    !! Read FORCE_CONSTANTS_3RD, generated by thirdorder.py, under ShengBTE package
    !! Regading the format of 3rd force constants, more details can be found in the 
    !! documentations of ShengBTE code :
    !! https://bitbucket.org/sousaw/shengbte/src/master/README.md
    !--------------------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE mp_world,         ONLY : mpime
    USE io_global,        ONLY : stdout, ionode, ionode_id
    USE mp_global,        ONLY : inter_pool_comm, npool, my_pool_id    
    USE mp,               ONLY : mp_barrier, mp_bcast
    USE tdbecom,          ONLY : nifc, Phi, rvec2, rvec3, ind_atm1,   &
                                 ind_atm2, ind_atm3, lattv
    USE io_var,           ONLY : iun3rdfc
    USE constants_epw,    ONLY : bohr2ang, ryd2ev
    !
    implicit none
    !
    INTEGER :: ii, jj, ll, mm, nn, ltem, mtem, ntem, info, P(3)
    !!
    INTEGER :: ierr
    !! Error status
    INTEGER :: ios
    !! I/O status 
    REAL(KIND = DP) :: tmp(3,3)
    !! lattice constants
    !
    IF (mpime == ionode_id) THEN
      OPEN(UNIT = iun3rdfc, FILE = 'FORCE_CONSTANTS_3RD', status="old", IOSTAT = ios)
      IF (ios /= 0) CALL errore('read3fc', 'error opening file FORCE_CONSTANTS_3RD', 1)
      READ(iun3rdfc, *) nifc
      ALLOCATE(ind_atm1(nifc), STAT =ierr)
      IF (ierr /= 0) CALL errore('read3fc', 'Error allocating ind_atm1', 1)
      ALLOCATE(ind_atm2(nifc), STAT =ierr)
      IF (ierr /= 0) CALL errore('read3fc', 'Error allocating ind_atm2', 1)
      ALLOCATE(ind_atm3(nifc), STAT =ierr)
      IF (ierr /= 0) CALL errore('read3fc', 'Error allocating ind_atm3', 1)
      ALLOCATE(Phi(3, 3, 3, nifc), STAT = ierr)
      IF (ierr /= 0) CALL errore('read3fc', 'Error allocating Phi', 1)
      ALLOCATE(rvec2(3, nifc), STAT = ierr)
      IF (ierr /= 0) CALL errore('read3fc', 'Error allocating rvec2', 1)
      ALLOCATE(rvec3(3, nifc), STAT = ierr)
      IF (ierr /= 0) CALL errore('read3fc', 'Error allocating rvec3', 1)
      DO ii = 1, nifc
        READ(iun3rdfc, *) jj
        READ(iun3rdfc, *) rvec2(:,ii)
        READ(iun3rdfc, *) rvec3(:,ii)
        READ(iun3rdfc, *) ind_atm1(ii), ind_atm2(ii), ind_atm3(ii)
        DO ll=1, 3
          DO mm=1, 3
            DO nn=1, 3
              READ(iun3rdfc, *) ltem, mtem, ntem, Phi(ll,mm,nn,ii)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CLOSE(iun3rdfc)
      Phi = Phi * bohr2ang**3 / ryd2ev  ! transform to Rydberg atomic unit
      ! write(stdout,*) 'Read the unit cell coordinates used to calculate FORCE_CONSTANTS_3RD'
      ! open(2,file='sBTE_config',status='old',IOSTAT = ios)
      ! IF (ios /= 0) CALL errore('read3fc', 'error opening file sBTE_config')
      ! READ(2,*) lfactor
      ! DO ll =1,3
      !    READ(2,*) lattv(:,ll)
      ! ENDDO
      ! close(2)
      !lattv=lfactor*lattv
      !! Each vector is rounded to the nearest lattice vector.
      !! rvec2 and rvec3 are transformed to crystal coordinates, which 
      !! is different from the original code ShengBTE/processes.f90 
      tmp = lattv * bohr2ang /10.d0 
      !! to nm, in order to transform unit cell coordinates back to
      !! crystal coordinates
      CALL DGESV(3, nifc, tmp, 3, P, rvec2, 3, info)
      rvec2 = ANINT(rvec2 / 10.d0)     
      tmp = lattv * bohr2ang /10.d0 
      CALL DGESV(3, nifc, tmp, 3, P, rvec3, 3, info)
      rvec3 = ANINT(rvec3 / 10.d0)
    ENDIF ! mpime == ionode_id
    CALL mp_bcast(nifc, ionode_id, inter_pool_comm)
    IF (mpime /= ionode_id) THEN
        ALLOCATE(ind_atm1(nifc), STAT =ierr)
        IF (ierr /= 0) CALL errore('read3fc', 'Error allocating ind_atm1', 1)
        ALLOCATE(ind_atm2(nifc), STAT =ierr)
        IF (ierr /= 0) CALL errore('read3fc', 'Error allocating ind_atm2', 1)
        ALLOCATE(ind_atm3(nifc), STAT =ierr)
        IF (ierr /= 0) CALL errore('read3fc', 'Error allocating ind_atm3', 1)
        ALLOCATE(Phi(3, 3, 3, nifc), STAT = ierr)
        IF (ierr /= 0) CALL errore('read3fc', 'Error allocating Phi', 1)
        ALLOCATE(rvec2(3, nifc), STAT = ierr)
        IF (ierr /= 0) CALL errore('read3fc', 'Error allocating rvec2', 1)
        ALLOCATE(rvec3(3, nifc), STAT = ierr)
        IF (ierr /= 0) CALL errore('read3fc', 'Error allocating rvec3', 1)
    ENDIF
    CALL mp_bcast(ind_atm1, ionode_id, inter_pool_comm)
    CALL mp_bcast(ind_atm2, ionode_id, inter_pool_comm)
    CALL mp_bcast(ind_atm3, ionode_id, inter_pool_comm)
    CALL mp_bcast(Phi,      ionode_id, inter_pool_comm)
    CALL mp_bcast(rvec2,    ionode_id, inter_pool_comm)
    CALL mp_bcast(rvec3,    ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    RETURN
    !
    !--------------------------------------------------------------------------------        
    END SUBROUTINE read3fc
    !--------------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------------    
    SUBROUTINE deallocate_phph()
    !--------------------------------------------------------------------------------  
    !! deallocate the 3-phonon matrix elements and their indxes.
    !-------------------------------------------------------------------------------- 
    !
    USE epwcom,     ONLY : ph_rta, phself_anh
    USE tdbecom,    ONLY : ind_atm1, ind_atm2, ind_atm3, Phi, rvec2,  &
                           rvec3, psi2_p, psi2_m, ind_p, ind_m, ph_lw
    !
    IMPLICIT NONE 
    !
    INTEGER :: ierr
    !
    IF (ph_rta) THEN
      DEALLOCATE(ph_lw, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_phph', 'Error deallocating psi2_p', 1)
    ELSE
      DEALLOCATE(ind_atm1, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_phph', 'Error deallocating ind_atm1', 1)
      DEALLOCATE(ind_atm2, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_phph', 'Error deallocating ind_atm2', 1)
      DEALLOCATE(ind_atm3, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_phph', 'Error deallocating ind_atm3', 1)
      DEALLOCATE(Phi,    STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_phph', 'Error deallocating Phi', 1)
      DEALLOCATE(rvec2,    STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_phph', 'Error deallocating rvec2', 1)
      DEALLOCATE(rvec3,    STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_phph', 'Error deallocating rvec3', 1)
      DEALLOCATE(psi2_p, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_phph', 'Error deallocating psi2_p', 1)
      DEALLOCATE(psi2_m, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_phph', 'Error deallocating psi2_m', 1)
      DEALLOCATE(ind_p, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_phph', 'Error deallocating ind_p', 1)
      DEALLOCATE(ind_m, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_phph', 'Error deallocating ind_m', 1)
    ENDIF
    !
    !--------------------------------------------------------------------------------  
    END SUBROUTINE deallocate_phph
    !-------------------------------------------------------------------------------- 
    !    
    !--------------------------------------------------------------------------------  
    SUBROUTINE psi_interp(nu1, nu2, nu3, iq1, iq2, iq3, pm, psi2)
    !-------------------------------------------------------------------------------- 
    !! Calculate the modulus square of ph-ph matrix element from 3rd force constants
    !! Expressions can be found in Eq.(9) of Comp. Phys. Commun. 185, 17471758 (2014)
    !-------------------------------------------------------------------------------- 
    USE kinds,            ONLY : DP 
    USE modes,            ONLY : nmodes
    USE tdbecom,          ONLY : uf_all, Phi, rvec2, rvec3, ind_atm1,  &
                                 ind_atm2, ind_atm3, nifc, Phi,      &
                                 rvec2, rvec3
    USE elph2,            ONLY : wf, xqf
    USE ions_base,        ONLY : nat, amass, ityp, tau
    USE constants_epw,    ONLY : twopi, ci, eps4
    !
    IMPLICIT NONE 
    !   
    INTEGER, INTENT(IN) :: nu1
    !! branch index of 1st phonon
    INTEGER, INTENT(IN) :: nu2
    !! branch index of 2nd phonon
    INTEGER, INTENT(IN) :: nu3
    !! branch index of 3rd phonon
    INTEGER, INTENT(IN) :: iq1
    !! q-index for 1st phonon
    INTEGER, INTENT(IN) :: iq2
    !! q-index for 2nd phonon
    INTEGER, INTENT(IN) :: iq3 
    !! q-index for 3rd phonon
    REAL(KIND = DP), INTENT(IN) :: pm
    !! plus or minus
    REAL(KIND = DP), INTENT(OUT) :: psi2
    !! results of interpolation
    INTEGER :: ll
    !! index of force constants
    INTEGER :: rr
    !! direction index
    INTEGER :: ss
    !! direction index
    INTEGER :: tt
    !! direction index
    REAL(KIND = DP) :: wqnu1
    !! frequency of 1st phonon
    REAL(KIND = DP) :: wqnu2 
    !! frequency of 2nd phonon
    REAL(KIND = DP) :: wqnu3
    !! frequency of 3rd phonon
    REAL(KIND = DP) :: xxq2(3)
    !! coordinate  of 2nd phonon
    REAL(KIND = DP) :: xxq3(3)
    !! coordinate  of 3rd phonon
    REAL(KIND = DP) :: xxr2(3)
    !! coordinate  of 2nd phonon
    REAL(KIND = DP) :: xxr3(3)
    !! coordinate  of 3rd phonon
    COMPLEX(KIND = DP) :: vectq1(nmodes, nmodes)
    !! eigenmode 1
    COMPLEX(KIND = DP) :: vectq2(nmodes, nmodes)
    !! eigenmode 2
    COMPLEX(KIND = DP) :: vectq3(nmodes, nmodes)
    !! eigenmode 3
    COMPLEX(KIND = DP) :: psi_plus
    !! In case of plus
    COMPLEX(KIND = DP) :: psi_minus
    !! In case of minus
    COMPLEX(KIND = DP) :: prefac, Vp0     
    ! 
    xxq2  = xqf(:, iq2)
    xxq3 =  xqf(:, iq3)
    wqnu1  = wf(nu1, iq1)
    wqnu2  = wf(nu2, iq2)
    wqnu3  = wf(nu3, iq3)
    vectq1 = uf_all(iq1, :, :)
    vectq2 = uf_all(iq2, :, :)
    vectq3 = uf_all(iq3, :, :)
    psi_plus = 0.d0
    psi_minus = 0.d0   
    psi2  = 0.d0         
    IF (ABS(pm - 1.0) < eps4) THEN
      ! 
      DO ll = 1, nifc
        prefac = 1.d0 / SQRT(amass(ityp(ind_atm1(ll)))  *  &
                             amass(ityp(ind_atm2(ll)))  *  &
                             amass(ityp(ind_atm3(ll))))
        xxr2 = rvec2(:,ll)
        xxr3 = rvec3(:,ll)
        prefac = prefac * EXP( twopi * ci * DOT_PRODUCT(xxq2, xxr2)) * &
                          EXP(-twopi * ci * DOT_PRODUCT(xxq3, xxr3))
        Vp0=0.d0
        DO rr = 1, 3
          DO ss = 1, 3
            DO tt = 1, 3
              Vp0 = Vp0 + Phi(tt, ss, rr, ll) * &
                    vectq1(tt + 3 * (ind_atm1(ll) - 1), nu1) * &
                    vectq2(ss + 3 * (ind_atm2(ll) - 1), nu2) * &
              CONJG(vectq3(rr + 3 * (ind_atm3(ll) - 1), nu3))
            END DO
          END DO
        END DO
        psi_plus = psi_plus + prefac * Vp0
      END DO
      psi2= ABS(psi_plus) ** 2 /(wqnu1 * wqnu2 * wqnu3)
      !
    ELSEIF(ABS(pm + 1.0) < eps4) THEN
      ! 
      DO ll = 1, nifc
        prefac = 1.d0 / SQRT(amass(ityp(ind_atm1(ll)))  *  &
                             amass(ityp(ind_atm2(ll)))  *  &
                             amass(ityp(ind_atm3(ll))))
        xxr2 = rvec2(:,ll)
        xxr3 = rvec3(:,ll)
        prefac = prefac * EXP(-twopi * ci * DOT_PRODUCT(xxq2, xxr2)) * &
                          EXP(-twopi * ci * DOT_PRODUCT(xxq3, xxr3))
        Vp0=0.d0
        DO rr = 1, 3
          DO ss = 1, 3
            DO tt = 1, 3
              Vp0 = Vp0 + Phi(tt, ss, rr, ll) * &
                    vectq1(tt + 3 * (ind_atm1(ll) - 1), nu1)  * &
              CONJG(vectq2(ss + 3 * (ind_atm2(ll) - 1), nu2)) * &
              CONJG(vectq3(rr + 3 * (ind_atm3(ll) - 1), nu3))
            END DO
          END DO
        END DO
        psi_minus = psi_minus + prefac * Vp0
      END DO
      psi2= ABS(psi_minus) ** 2 /(wqnu1 * wqnu2 * wqnu3)
      !
    ELSE 
      CALL errore('psi_interp', 'pm must be 1.0 or -1.0', 1)
    ENDIF
    !
    RETURN
    !
    !-------------------------------------------------------------------------------- 
    END SUBROUTINE psi_interp
    !-------------------------------------------------------------------------------- 
    !
    !-------------------------------------------------------------------------------- 
    FUNCTION base_sigma(v)
    !-------------------------------------------------------------------------------- 
    !! Return the base broadening (without prefac) for a mode.
    !-------------------------------------------------------------------------------- 
    USE kinds,            ONLY : DP 
    USE tdbecom,          ONLY : rlattv
    USE constants_epw,    ONLY : ryd2mev 
    USE epwcom,           ONLY : nqf1, nqf2, nqf3
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(IN) :: v(3)
    REAL(KIND = DP):: base_sigma
    REAL(KIND = DP) :: sigma_tmp(3)
    REAL(KIND = DP) :: etamin
    INTEGER :: idim
    !
    etamin = 0.05 / ryd2mev 
    !
    base_sigma=0.
    sigma_tmp(1) = DOT_PRODUCT(rlattv(:,1), v)/ DBLE(nqf1)
    sigma_tmp(2) = DOT_PRODUCT(rlattv(:,2), v)/ DBLE(nqf2)
    sigma_tmp(3) = DOT_PRODUCT(rlattv(:,3), v)/ DBLE(nqf3)
    DO idim = 1,3
      base_sigma = base_sigma + sigma_tmp(idim) ** 2
    END DO
    base_sigma= SQRT(base_sigma / 6.d0)
    !! if base_sigma < 0.05 mev, set 0.05 mev
    IF (base_sigma < etamin ) THEN
      base_sigma = etamin
    ENDIF
    !
    RETURN
    !
    !-------------------------------------------------------------------------------- 
    END function base_sigma
    !-------------------------------------------------------------------------------- 
    !
    !--------------------------------------------------------------------------------     
    SUBROUTINE cross_product(a, b, res)
    !-------------------------------------------------------------------------------- 
    !! cross product of two vectors
    !-------------------------------------------------------------------------------- 
    USE kinds,            ONLY : DP 
    !
    IMPLICIT NONE 
    !
    REAL(KIND = DP), INTENT(IN)  :: a(3), b(3)
    REAL(KIND = DP), INTENT(OUT) :: res(3)
    INTEGER :: i,j,k
    DO i = 1, 3
      j = MOD(i, 3) + 1
      k = MOD(j, 3) + 1
      res(i) = a(j) * b(k) - a(k) * b(j)
    END DO
    !
    RETURN
    !--------------------------------------------------------------------------------     
    END SUBROUTINE cross_product
    !-------------------------------------------------------------------------------- 
    !
    !-------------------------------------------------------------------------------- 
    SUBROUTINE ph_energy()
    !-------------------------------------------------------------------------------- 
    !! Calculate a list of phonon energy for temperatues ranging from temp_ph to temp_el
    !! so that one can use the energy as a reference to find out the effective phonon
    !! temperature
    !-------------------------------------------------------------------------------- 
    !    
    USE kinds,            ONLY : DP 
    USE epwcom,           ONLY : temp_el, temp_ph, eps_acustic
    USE elph2,            ONLY : wf,nqtotf
    use modes,            ONLY : nmodes
    USE division,         ONLY : fkbounds
    USE io_global,        ONLY : stdout
    USE constants_epw,    ONLY : kelvin2Ry
    USE mp,               ONLY : mp_barrier, mp_sum, mp_bcast
    USE mp_global,        ONLY : inter_pool_comm, npool, my_pool_id 
    USE tdbecom,          ONLY : phtemp, e_latt, istart, ntph 
    !
    IMPLICIT NONE 
    !
    INTEGER :: ierr
    !! Error info
    INTEGER :: itemp
    !! Counter of temperature
    INTEGER :: iq 
    !! Counter of q points 
    INTEGER :: nu
    !! Counter of phonon branch
    INTEGER :: lower_qbnd
    !! lower q index
    INTEGER :: upper_qbnd
    !! upper q index
    REAL(KIND = DP) :: dtemp
    !! temperature difference
    REAL(KIND = DP) :: n_eq
    !! Equilibirum phonon population
    !
    WRITE(stdout,'(5x,a)') 'Phonon-phonon interaction is treated with RTA'
    WRITE(stdout,'(5x,a)') 'Phonon energy is evaluated on a list of temperature in the range of Tph,Tel'
    phtemp =0.d0
    e_latt = 0.d0
    dtemp = ABS(temp_el - temp_ph) / DBLE(ntph - 1)
    CALL fkbounds(nqtotf, lower_qbnd, upper_qbnd)
    DO itemp = 1, ntph
      phtemp(itemp) = temp_ph + dtemp*(itemp -1) !! Kelvin
      DO iq = lower_qbnd, upper_qbnd
        DO nu = 1, nmodes
          IF (wf(nu,iq) > eps_acustic) THEN
            n_eq = 1.0 / (EXP(wf(nu, iq) / (phtemp(itemp) * kelvin2Ry)) - 1.0)
            e_latt(itemp) = e_latt(itemp) + n_eq * wf(nu,iq) / nqtotf
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    CALL mp_sum(e_latt, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    istart = 1
    !
    RETURN
    !-------------------------------------------------------------------------------- 
    END SUBROUTINE ph_energy
    !-------------------------------------------------------------------------------- 
    !
    !-------------------------------------------------------------------------------- 
    SUBROUTINE read_ph_lw()
    !-------------------------------------------------------------------------------- 
    !! In case the collision integrals of ph-ph interaction are approximated with
    !! relaxation time approximation, the phonon linewidth can be read from external
    !! files.
    !-------------------------------------------------------------------------------- 
    USE mp_global,        ONLY : inter_pool_comm, npool, my_pool_id
    USE mp_world,         ONLY : mpime
    use modes,            ONLY : nmodes
    USE elph2,            ONLY : wf,nqtotf
    USE io_global,        ONLY : stdout, ionode,ionode_id
    USE mp,               ONLY : mp_barrier, mp_sum, mp_bcast
    USE mp_global,        ONLY : inter_pool_comm, npool, my_pool_id
    USE tdbecom,          ONLY : ph_lw 
    !
    IMPLICIT NONE
    !
    INTEGER :: ios
    !! error info
    INTEGER :: ierr
    !! error info
    INTEGER :: iq
    !! index q
    INTEGER :: nu
    !! index of phonon branch
    WRITE(stdout, '(5x,a)') 'Read phonon linewidth from file ph_lw.dat'
    ALLOCATE (ph_lw(nmodes,nqtotf),STAT =ierr)
    ph_lw = 0.d0
    IF (ierr /= 0) CALL errore('read_ph_lw',' Error allcoating ph_lw',1)
    IF (mpime == ionode_id) THEN
      open(unit = 2, file = 'ph_lw.dat',status="old",IOSTAT = ios)
      DO iq = 1, nqtotf
        DO nu = 1, nmodes
          READ(2,*) ph_lw(nu,iq)
        ENDDO
      ENDDO
      close(2)
      !ph_lw = ph_lw /Ry2THz
    ENDIF
    CALL mp_bcast(ph_lw, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !-------------------------------------------------------------------------------- 
    END SUBROUTINE read_ph_lw
    !-------------------------------------------------------------------------------- 
    !
    !-------------------------------------------------------------------------------- 
    SUBROUTINE iphph_rta_calc()
    !-------------------------------------------------------------------------------- 
    !! This routine calculates ph-ph collision integral 
    !! under relaxation time approximation
    !--------------------------------------------------------------------------------     
    USE kinds,            ONLY : DP 
    USE epwcom,           ONLY : eps_acustic, tstep_write
    USE elph2,            ONLY : wf,nqtotf
    USE modes,            ONLY : nmodes
    USE division,         ONLY : fkbounds
    USE mp,               ONLY : mp_barrier, mp_sum, mp_bcast
    USE mp_global,        ONLY : inter_pool_comm, npool, my_pool_id 
    USE io_global,        ONLY : stdout
    USE tdbecom,          ONLY : nphon_pre, iph_ph, a_tab, b_tab,  &
                                 c_tab, nstg, dt_in_ps, e_latt,    &
                                 phtemp, ph_lw, istart, ntph    
    USE constants_epw,    ONLY : kelvin2Ry
    !
    IMPLICIT NONE
    !
    INTEGER :: lower_qbnd
    !! lower q index
    INTEGER :: upper_qbnd
    !! upper q index
    INTEGER :: iq
    !! q index
    INTEGER :: nu 
    !! phonon branch index
    INTEGER :: istart2
    !! staring index for the guess of temperature 
    INTEGER :: itemp
    !!
    INTEGER :: istg, pp 
    !! Counter on Rugge Kutta solver
    REAL(KIND = DP) :: et
    !! Total energy of the nonequilibrium phonon
    REAL(KIND = DP) :: demin
    !! keep record of min of |et - e_latt(itemp)|
    REAL(KIND = DP) :: detmp 
    !! et - e_latt(itemp)
    REAL(KIND = DP) :: temp_it
    !! Effective lattice temperarure  
    REAL(KIND = DP) :: n_eq
    !! Equilibrium phonon population at temperature
    REAL(KIND = DP) :: nph
    !! Nonqquilibrium phonon population
    !! temp_it
    !
    CALL fkbounds(nqtotf, lower_qbnd, upper_qbnd)
    !
    et = 0.d0
    DO iq = 1, nqtotf
      DO nu = 1, nmodes
        if (wf(nu, iq) > eps_acustic) then
          et = et + wf(nu,iq) * nphon_pre(nu,iq) / nqtotf
        ENDif
      ENDDO ! nu
    ENDDO ! iq

    IF (istart < 6) THEN
      istart2 = 1
    ELSE
      istart2 = istart - 1
    ENDIF
    !
    demin = 1000.d0
    !
    DO itemp = istart2, ntph
      detmp  = et - e_latt(itemp)
      IF (abs(detmp)< demin) THEN
        demin = abs(detmp)
        istart = itemp
      ELSE
        EXIT
      ENDIF
    ENDDO
    !
    temp_it = phtemp(istart)
    WRITE(stdout, '(5x,a,f10.3,a)') 'lattice temperature : ', temp_it, ' K'
    !
    iph_ph = 0.d0
    DO istg = 1, nstg
      DO iq = lower_qbnd, upper_qbnd
        DO nu = 1, nmodes 
          IF (wf(nu, iq) > eps_acustic) then
            n_eq = 1.d0 /(EXP(wf(nu, iq) / (temp_it * kelvin2Ry)) - 1.d0)
            nph = nphon_pre(nu,iq)
            IF (istg > 1) THEN
              DO pp = 1, istg - 1 
                nph = nph + a_tab(istg, pp) * iph_ph(nu, iq, pp) * dt_in_ps
              ENDDO
            ENDIF
            iph_ph(nu, iq, istg) = (n_eq-nph)*ph_lw(nu,iq)
          ENDIF
        ENDDO ! nu
      ENDDO ! iq 
      CALL mp_sum(iph_ph(:, :, istg), inter_pool_comm)    
      CALL mp_barrier(inter_pool_comm)
    ENDDO ! istg
    !
    RETURN
    !
    !-------------------------------------------------------------------------------- 
    END SUBROUTINE iphph_rta_calc
    !-------------------------------------------------------------------------------- 
  !
  !-------------------------------------------------------------------------------- 
  END MODULE phph
  !-------------------------------------------------------------------------------- 
