  !---------------------------------------------------------------------  
  SUBROUTINE tdbe()
  !---------------------------------------------------------------------  
  !!
  !---------------------------------------------------------------------  
  USE epwcom,           ONLY : temp_ph, bte_restart,  bte_phph, mp_mesh_k,        &
                               ph_rta, phself_anh, phphself_r         
  USE tdbe_mod,         ONLY : map2nkfs, double_grid_map,read_g2_tdbe,            &
                               deallocate_tdbe, allocate_tdbe,  &
                               wanninterp, bte_restart_read, &
                               set_initial_dist, write_dist, write_points,        &
                               double_grid_map, read_g2_tdbe, map2nkfs, check_dg, &
                               time_propagation, deallocate_tdbe   
  USE phph,             ONLY : psi2_calc,  inv_tau_rta, ph_energy, read_ph_lw,    &
                               deallocate_phph, phphself
  ! USE phph,             ONLY : ph_energy, read_ph_lw, iphph_rta_calc
  ! USE epwcom,           ONLY : bte_elel
  ! USE elec_elec,        ONLY : elec_energy, iee_rta
  !
  !!
  IMPLICIT NONE
  !
  LOGICAL :: flag
  !
  flag = (phself_anh .AND. (.NOT. phphself_r) .AND. bte_phph)
  !
  CALL start_clock('TDBE')
  CALL wanninterp()
  IF (phself_anh .AND. phphself_r) THEN
    CALL phphself(temp_ph, temp_ph + 600.d0, 6)
    CALL deallocate_phph()
    CALL stop_clock('TDBE')
    RETURN
  ENDIF
  !! Allocate space for time-dependent electron and phonon distribution, and collision inegrals 
  CALL allocate_tdbe() 
  !!
  IF (bte_restart) THEN
    CALL bte_restart_read()  
  ELSE  
    CALL set_initial_dist() 
    CALL write_dist(.TRUE., 0.0d0)
  ENDIF
  !!
  IF (.NOT. flag) THEN
  !!
    CALL write_points()
    CALL double_grid_map()
    CALL read_g2_tdbe()
    CALL map2nkfs() 
    IF (mp_mesh_k) THEN
      CALL check_dg() 
    ENDIF
  ENDIF
  !!
  IF (bte_phph) THEN
    IF (ph_rta) THEN 
      CALL ph_energy()
      CALL read_ph_lw()
    ELSE  
      CALL psi2_calc()
      CALL inv_tau_rta(temp_ph, temp_ph + 600.d0, 6)
    ENDIF
  ENDIF
  ! IF (bte_elel) THEN
  !   ! CALL elec_energy(etf_all,wkf_all,ef0)
  ! ENDIF

  !! Now start the tdBE calulation
  !!
  IF (.NOT. flag) THEN
    CALL time_propagation()
    CALL deallocate_tdbe()
    IF (bte_phph) THEN
      CALL deallocate_phph()
    ENDIF
  ELSE
    CALL deallocate_phph()
  ENDIF
  CALL stop_clock('TDBE')
  !
  !---------------------------------------------------------------------  
  END SUBROUTINE tdbe
  !---------------------------------------------------------------------    

