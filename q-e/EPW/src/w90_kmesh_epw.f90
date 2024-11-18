!-*- mode: F90 -*-!
!------------------------------------------------------------!
! This file is distributed as part of the Wannier90 code and !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the Wannier90      !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the Wannier90 code is www.wannier.org       !
!                                                            !
! The Wannier90 code is hosted on GitHub:                    !
!                                                            !
! https://github.com/wannier-developers/wannier90            !
!------------------------------------------------------------!

module w90_kmesh_epw
  !! Routines adapeted from wannier90 code w90_kmesh.f90
  !! to analyse the regular k-point mesh and determine 
  !! the overlaps neccessary for a finite
  !! difference representation of the spread operator.
  !! These overlaps are defined by a set of vectors (b-vectors) which
  !! connect the Bloch states.
  !! See Eq. B1 in Appendix B of Marzari and
  !!  Vanderbilt  PRB 56 12847 (1997)

  use kinds, only: dp
  use io_global, ONLY : stdout
  use mp,        ONLY : mp_bcast, mp_barrier, mp_sum, mp_max, mp_min
  USE cell_base,        ONLY : bg
  USE mp_global,        ONLY : inter_pool_comm, world_comm,npool, my_pool_id
  implicit none

  private
  ! Definitions of system variables set in this module
  ! nnh     ! the number of b-directions (bka)
  ! nntot   ! total number of neighbours for each k-point
  ! nnlist  ! list of neighbours for each k-point
  ! neigh
  ! nncell  ! gives BZ of each neighbour of each k-point
  ! wbtot
  ! wb      ! weights associated with neighbours of each k-point
  ! bk      ! the b-vectors that go from each k-point to its neighbours
  ! bka     ! the b-directions (not considering inversion) from
  ! 1st k-point to its neighbours

  public :: kmesh_get
!  public :: kmesh_write
  public :: kmesh_dealloc
  public :: kmesh_gradient

  ! public variables 
  integer, public, save              :: num_shells
  integer, public, save              :: nnh           ! the number of b-directions (bka)
  integer, public, save              :: nntot         ! total number of neighbours for each k-point
  integer, public, save, allocatable :: nnlist(:, :)   ! list of neighbours for each k-point
  integer, public, save, allocatable :: neigh(:, :)
  integer, public, save, allocatable :: nncell(:, :, :) ! gives BZ of each neighbour of each k-point
  real(kind=dp), public, save              :: wbtot
  real(kind=dp), public, save, allocatable :: wb(:)         ! weights associated with neighbours of each k-point
  real(kind=dp), public, save, allocatable :: bk(:, :, :)     ! the b-vectors that go from each k-point to its neighbours
  real(kind=dp), public, save, allocatable :: bka(:, :)      ! the b-directions from 1st k-point to its neighbours
  integer, allocatable, public, save :: shell_list(:)

  ! private variables
  integer, parameter :: nsupcell = 3
  integer, parameter :: search_shells = 72
  integer, parameter :: max_shells = 6
  integer, parameter :: num_nnmax = 12
  real(kind=dp), parameter :: kmesh_tol = 0.000001_dp
  !! Size of supercell (of recip cell) in which to search for k-point shells
  integer            :: lmn(3, (2*nsupcell + 1)**3)
  !! Order in which to search the cells (ordered in dist from origin)
  real(kind=dp)      :: bvec_inp(3, num_nnmax, max_shells)
  real(kind=dp), allocatable      :: kpt_cart(:,:)
  integer :: num_kpts
  ! variables for k parallization
  integer :: lower_bnd, upper_bnd
  !! The input b-vectors (only used in the rare case they are read from file)

contains
  !=======================================================
  subroutine kmesh_get(xkf_all,nkpts)
    !=====================================================
    !
    !! Main routine to calculate the b-vectors
    !
    !=====================================================
    use division,  ONLY : fkbounds
    implicit none

    real(kind = dp), intent(in) :: xkf_all(3,nkpts)
    integer, intent(in) :: nkpts 
    ! Variables that are private
    integer :: nlist, nkp, nkp2, l, m, n, ndnn, ndnnx, ndnntot
    integer :: nnsh, nn, nnx, loop, i, j
    integer :: ifound, counter, na, nap, loop_s, loop_b, shell, nbvec, bnum
    integer :: ifpos, ifneg, ierr, multi(search_shells)
    integer :: nnshell(nkpts, search_shells)
    real(kind=dp) :: vkpp(3), vkpp2(3)
    real(kind=dp) :: dist, dnn0, dnn1, bb1, bbn, ddelta
    real(kind=dp), parameter :: eta = 99999999.0_dp    ! eta = very large
    real(kind=dp) :: bweight(max_shells)
    real(kind=dp) :: dnn(search_shells)
    real(kind=dp) :: wb_local(num_nnmax)
    real(kind=dp) :: bk_local(3, num_nnmax, nkpts), kpbvec(3)
    real(kind=dp), allocatable :: bvec_tmp(:, :)
    

    write (stdout, '(/1x,a)') &
      '*---------------------------------- K-MESH ----------------------------------*'
    num_kpts = nkpts
    ALLOCATE(kpt_cart(3,num_kpts), STAT = ierr)
    IF (ierr /= 0) CALL errore('kmesh_get', 'Error allocating kpt_cart', 1)
    kpt_cart = xkf_all
    call fkbounds(num_kpts, lower_bnd, upper_bnd)
    call cryst_to_cart(num_kpts, kpt_cart, bg, 1)
    ! Sort the cell neighbours so we loop in order of distance from the home shell
    call kmesh_supercell_sort

    ! find the distance between k-point 1 and its nearest-neighbour shells
    ! if we have only one k-point, the n-neighbours are its periodic images

    dnn0 = 0.0_dp
    dnn1 = eta
    ndnntot = 0
    do nlist = 1, search_shells
      do nkp = 1, num_kpts
        do loop = 1, (2*nsupcell + 1)**3
          l = lmn(1, loop); m = lmn(2, loop); n = lmn(3, loop)
          !
          vkpp = kpt_cart(:, nkp) + matmul(lmn(:, loop), bg)
          dist = sqrt((kpt_cart(1, 1) - vkpp(1))**2 &
                      + (kpt_cart(2, 1) - vkpp(2))**2 + (kpt_cart(3, 1) - vkpp(3))**2)
          !
          if ((dist .gt. kmesh_tol) .and. (dist .gt. dnn0 + kmesh_tol)) then
            if (dist .lt. dnn1 - kmesh_tol) then
              dnn1 = dist  ! found a closer shell
              counter = 0
            end if
            if (dist .gt. (dnn1 - kmesh_tol) .and. dist .lt. (dnn1 + kmesh_tol)) then
              counter = counter + 1 ! count the multiplicity of the shell
            end if
          end if
        end do ! loop
      enddo ! nkp
      if (dnn1 .lt. eta - kmesh_tol) ndnntot = ndnntot + 1
      dnn(nlist) = dnn1
      multi(nlist) = counter
      dnn0 = dnn1
      dnn1 = eta
    enddo
    write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
    write (stdout, '(1x,a)') '|                    Distance to Nearest-Neighbour Shells                    |'
    write (stdout, '(1x,a)') '|                    ------------------------------------                    |'
    write (stdout, '(1x,a)') '|          Shell             Distance (Bohr^-1)         Multiplicity         |'
    write (stdout, '(1x,a)') '|          -----             ------------------         ------------         |'
    do ndnn = 1, ndnntot
      write (stdout, '(1x,a,11x,i3,17x,f10.6,19x,i4,12x,a)') '|', ndnn, dnn(ndnn), multi(ndnn), '|'
    enddo
    write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
    
    !! we decide the num_shells automatically
    num_shells = 0
    ! if (num_shells == 0) then
    call kmesh_shell_automatic(multi, dnn, bweight)
    ! elseif (num_shells > 0) then
    !   call kmesh_shell_fixed(multi, dnn, bweight)
    ! end if
   
    write (stdout, '(1x,a)', advance='no') '| The following shells are used: '
    do ndnn = 1, num_shells
      if (ndnn .eq. num_shells) then
        write (stdout, '(i3,1x)', advance='no') shell_list(ndnn)
      else
        write (stdout, '(i3,",")', advance='no') shell_list(ndnn)
      endif
    enddo
    do l = 1, 11 - num_shells
      write (stdout, '(4x)', advance='no')
    enddo
    write (stdout, '("|")')

    nntot = 0
    do loop_s = 1, num_shells
      nntot = nntot + multi(shell_list(loop_s))
    end do

    if (nntot > num_nnmax) then
      write (stdout, '(a,i2,a)') ' **WARNING: kmesh has found >', num_nnmax, ' nearest neighbours**'
      write (stdout, '(a)') ' '
      write (stdout, '(a)') ' This is probably caused by an error in your unit cell specification'
      write (stdout, '(a)') ' '
      write (stdout, '(a)') ' If you think this is not the problem; please send your *.win file to the '
      write (stdout, '(a)') ' wannier90 developers'
      write (stdout, '(a)') ' '
      write (stdout, '(a)') ' The problem may be caused by having accidentally degenerate shells of '
      write (stdout, '(a)') ' kpoints. The solution is then to rerun wannier90 specifying the b-vectors '
      write (stdout, '(a)') ' in each shell.  Give devel_flag=kmesh_degen in the *.win file'
      write (stdout, '(a)') ' and create a *.kshell file:'
      write (stdout, '(a)') ' '
      write (stdout, '(a)') ' $>   cat hexagonal.kshell'
      write (stdout, '(a)') ' $>   1 2'
      write (stdout, '(a)') ' $>   5 6 7 8'
      write (stdout, '(a)') ' '
      write (stdout, '(a)') ' Where each line is a new shell (so num_shells in total)'
      write (stdout, '(a)') ' The elements are the bvectors labelled according to the following '
      write (stdout, '(a)') ' list (last column is distance)'
      write (stdout, '(a)') ' '
      allocate (bvec_tmp(3, maxval(multi)), stat=ierr)
      if (ierr /= 0) call errore('kmesh_get','Error allocating bvec_tmp in kmesh_get',1)
      bvec_tmp = 0.0_dp
      counter = 0
      do shell = 1, search_shells
        call kmesh_get_bvectors(multi(shell), 1, dnn(shell), bvec_tmp(:, 1:multi(shell)))
        do loop = 1, multi(shell)
          counter = counter + 1
          write (stdout, '(a,I4,1x,a,2x,3f12.6,2x,a,2x,f12.6,a)') ' | b-vector  ', counter, ': (', &
            bvec_tmp(:, loop), ')', dnn(shell), '  |'
        end do
      end do
      write (stdout, '(a)') ' '
      deallocate (bvec_tmp,stat=ierr)
      if (ierr /= 0) call errore('kmesh_get','Error deallocating bvec_tmp',1)
      call errore('kmesh_get','something wrong, found too many nearest neighbours',1)
    end if

    allocate (nnlist(num_kpts, nntot), stat=ierr)
    if (ierr /= 0) call errore('kmesh_get','Error in allocating nnlist',1)
    allocate (neigh(num_kpts, nntot/2), stat=ierr)
    if (ierr /= 0) call errore('kmesh_get','Error in allocating neigh',1)
    allocate (nncell(3, num_kpts, nntot), stat=ierr)
    if (ierr /= 0) call errore('kmesh_get','Error in allocating nncell',1)
    allocate (wb(nntot), stat=ierr)
    if (ierr /= 0) call errore('kmesh_get','Error in allocating wb',1)
    allocate (bka(3, nntot/2), stat=ierr)
    if (ierr /= 0) call errore('kmesh_get','Error in allocating bka',1)
    allocate (bk(3, nntot, num_kpts), stat=ierr)
    if (ierr /= 0) call errore('kmesh_get','Error in allocating bk',1)

    nnx = 0
    do loop_s = 1, num_shells
      do loop_b = 1, multi(shell_list(loop_s))
        nnx = nnx + 1
        wb_local(nnx) = bweight(loop_s)
      end do
    end do

    ! Now build up the list of nearest-neighbour shells for each k-point.
    ! nnlist(nkp,1...nnx) points to the nnx neighbours (ordered along increa
    ! shells) of the k-point nkp. nncell(i,nkp,nnth) tells us in which BZ is
    ! nnth nearest-neighbour of the k-point nkp. Construct the nnx b-vectors
    ! go from k-point nkp to each neighbour bk(1:3,nkp,1...nnx).
    ! Comment: Now we have bk(3,nntot,num_kps) 09/04/2006
    write (stdout, '(1x,a)') '+----------------------------------------------------------------------------+'
    write (stdout, '(1x,a)') '|                        Shell   # Nearest-Neighbours                        |'
    write (stdout, '(1x,a)') '|                        -----   --------------------                        |'

    ! if (index(devel_flag, 'kmesh_degen') == 0) then
      !
      ! Standard routine
      !
    nnshell = 0
    nncell  = 0
    nnlist  = 0
    neigh   = 0
    bk_local = 0.0_dp
    ! do nkp = 1, num_kpts
    do nkp = lower_bnd, upper_bnd
      nnx = 0
      ok: do ndnnx = 1, num_shells
        ndnn = shell_list(ndnnx)
        do loop = 1, (2*nsupcell + 1)**3
          l = lmn(1, loop); m = lmn(2, loop); n = lmn(3, loop)
          vkpp2 = matmul(lmn(:, loop), bg)
          do nkp2 = 1, num_kpts
            vkpp = vkpp2 + kpt_cart(:, nkp2)
            dist = sqrt((kpt_cart(1, nkp) - vkpp(1))**2 &
                        + (kpt_cart(2, nkp) - vkpp(2))**2 + (kpt_cart(3, nkp) - vkpp(3))**2)
            if ((dist .ge. dnn(ndnn)*(1 - kmesh_tol)) .and. (dist .le. dnn(ndnn)*(1 + kmesh_tol))) then
              nnx = nnx + 1
              nnshell(nkp, ndnn) = nnshell(nkp, ndnn) + 1
              nnlist(nkp, nnx) = nkp2
              nncell(1, nkp, nnx) = l
              nncell(2, nkp, nnx) = m
              nncell(3, nkp, nnx) = n
              bk_local(:, nnx, nkp) = vkpp(:) - kpt_cart(:, nkp)
            endif
            !if we have the right number of neighbours we can exit
            if (nnshell(nkp, ndnn) == multi(ndnn)) cycle ok
          enddo
        enddo
        ! check to see if too few neighbours here
      end do ok
    end do

    call mp_sum(nncell, inter_pool_comm)
    call mp_sum(nnlist, inter_pool_comm)
    call mp_sum(nnshell,inter_pool_comm)
    call mp_sum(bk_local,inter_pool_comm)
    call mp_barrier(inter_pool_comm)

    do ndnnx = 1, num_shells
      ndnn = shell_list(ndnnx)
      write (stdout, '(1x,a,24x,i3,13x,i3,33x,a)') '|', ndnn, nnshell(1, ndnn), '|'
    end do
    write (stdout, '(1x,"+",76("-"),"+")')

    ! do nkp = 1, num_kpts
    do nkp = lower_bnd, upper_bnd
      nnx = 0
      do ndnnx = 1, num_shells
        ndnn = shell_list(ndnnx)
        do nnsh = 1, nnshell(nkp, ndnn)
          bb1 = 0.0_dp
          bbn = 0.0_dp
          nnx = nnx + 1
          do i = 1, 3
            bb1 = bb1 + bk_local(i, nnx, 1)*bk_local(i, nnx, 1)
            bbn = bbn + bk_local(i, nnx, nkp)*bk_local(i, nnx, nkp)
          enddo
          if (abs(sqrt(bb1) - sqrt(bbn)) .gt. kmesh_tol) then
            write (stdout, '(1x,2f10.6)') bb1, bbn
            call errore('kmesh_get','Non-symmetric k-point neighbours!',1)
          endif
        enddo
      enddo
    enddo

    ! now check that the completeness relation is satisfied for every kpoint
    ! We know it is true for kpt=1; but we check the rest to be safe.
    ! Eq. B1 in Appendix  B PRB 56 12847 (1997)
    ! if (.not. skip_B1_tests) then
    ! do nkp = 1, num_kpts
    do nkp = lower_bnd, upper_bnd
      do i = 1, 3
        do j = 1, 3
          ddelta = 0.0_dp
          nnx = 0
          do ndnnx = 1, num_shells
            ndnn = shell_list(ndnnx)
            do nnsh = 1, nnshell(1, ndnn)
              nnx = nnx + 1
              ddelta = ddelta + wb_local(nnx)*bk_local(i, nnx, nkp)*bk_local(j, nnx, nkp)
            enddo
          enddo
          if ((i .eq. j) .and. (abs(ddelta - 1.0_dp) .gt. kmesh_tol)) then
            write (stdout, '(1x,2i3,f12.8)') i, j, ddelta
            call errore('kmesh_get','Eq. (B1) not satisfied in kmesh_get (1)',1)
          endif
          if ((i .ne. j) .and. (abs(ddelta) .gt. kmesh_tol)) then
            write (stdout, '(1x,2i3,f12.8)') i, j, ddelta
            call errore('kmesh_get','Eq. (B1) not satisfied in kmesh_get (2)',1)
          endif
        enddo
      enddo
    enddo
    ! end if

    write (stdout, '(1x,a)') '| Completeness relation is fully satisfied [Eq. (B1), PRB 56, 12847 (1997)]  |'
    write (stdout, '(1x,"+",76("-"),"+")')

    !
    wbtot = 0.0_dp
    nnx = 0
    do ndnnx = 1, num_shells
      ndnn = shell_list(ndnnx)
      do nnsh = 1, nnshell(1, ndnn)
        nnx = nnx + 1
        wbtot = wbtot + wb_local(nnx)
      enddo
    enddo

    nnh = nntot/2
    ! make list of bka vectors from neighbours of first k-point
    ! delete any inverse vectors as you collect them
    na = 0
    do nn = 1, nntot
      ifound = 0
      if (na .ne. 0) then
        do nap = 1, na
          call utility_compar(bka(1, nap), bk_local(1, nn, 1), ifpos, ifneg)
          if (ifneg .eq. 1) ifound = 1
        enddo
      endif
      if (ifound .eq. 0) then
        !         found new vector to add to set
        na = na + 1
        bka(1, na) = bk_local(1, nn, 1)
        bka(2, na) = bk_local(2, nn, 1)
        bka(3, na) = bk_local(3, nn, 1)
      endif
    enddo

    if (na .ne. nnh) call errore('kmesh_get','Did not find right number of bk directions',1)

    write (stdout, '(1x,a)') '|                 b_k Vectors (Bohr^-1) and Weights (Bohr^2)                 |'
    write (stdout, '(1x,a)') '|                 ------------------------------------------                 |'
    write (stdout, '(1x,a)') '|            No.         b_k(x)      b_k(y)      b_k(z)        w_b           |'
    write (stdout, '(1x,a)') '|            ---        --------------------------------     --------        |'
    do i = 1, nntot
      write (stdout, '(1x,"|",11x,i3,5x,3f12.6,3x,f10.6,8x,"|")') &
        i, (bk_local(j, i, 1), j=1, 3), wb_local(i)
    enddo
    write (stdout, '(1x,"+",76("-"),"+")')
    write (stdout, '(1x,a)') '|                           b_k Directions (Bohr^-1)                         |'
    write (stdout, '(1x,a)') '|                           ------------------------                         |'
    write (stdout, '(1x,a)') '|            No.           x           y           z                         |'
    write (stdout, '(1x,a)') '|            ---        --------------------------------                     |'
    do i = 1, nnh
      write (stdout, '(1x,"|",11x,i3,5x,3f12.6,21x,"|")') i, (bka(j, i), j=1, 3)
    enddo
    write (stdout, '(1x,"+",76("-"),"+")')
    write (stdout, *) ' '

    ! find index array
    ! do nkp = 1, num_kpts
    do nkp = lower_bnd, upper_bnd
      do na = 1, nnh
        ! first, zero the index array so we can check it gets filled
        neigh(nkp, na) = 0
        ! now search through list of neighbours of this k-point
        do nn = 1, nntot
          call utility_compar(bka(1, na), bk_local(1, nn, nkp), ifpos, ifneg)
          if (ifpos .eq. 1) neigh(nkp, na) = nn
        enddo
        ! check found
        if (neigh(nkp, na) .eq. 0) then
          write (stdout, *) ' nkp,na=', nkp, na
          call errore('kmesh_get','failed to find neighbours for this kpoint',1)
        endif
      enddo
    enddo
    call mp_sum(neigh, inter_pool_comm)
    call mp_barrier(inter_pool_comm)
    
    !fill in the global arrays from the local ones
    do loop = 1, nntot
      wb(loop) = wb_local(loop)
    end do

    do loop_s = 1, num_kpts
      do loop = 1, nntot
        bk(:, loop, loop_s) = bk_local(:, loop, loop_s)
      end do
    end do
    return

  end subroutine kmesh_get

  !========================================
  subroutine kmesh_dealloc()
    !========================================
    !
    !!  Release memory from the kmesh module
    !   This routine now check to see if arrays
    !   are allocated, as there are some code
    !   paths that will not allocate on all nodes
    !========================================
    implicit none
    integer :: ierr
    ! Deallocate real arrays that are public
    deallocate (bk, stat=ierr)
    if (ierr /= 0) call errore('kmesh_dealloc','Error in deallocating bk',1)
    deallocate (bka, stat=ierr)
    if (ierr /= 0) call errore('kmesh_dealloc','Error in deallocating bka',1)
    deallocate (wb, stat=ierr)
    if (ierr /= 0) call errore('kmesh_dealloc','Error in deallocating wb in',1)
  ! Deallocate integer arrays that are public
    deallocate (neigh, stat=ierr)
    if (ierr /= 0) call errore('kmesh_dealloc','Error in deallocating neigh',1)
    deallocate (nncell, stat=ierr)
    if (ierr /= 0) call errore('kmesh_dealloc','Error in deallocating nncell',1)
    deallocate (nnlist, stat=ierr)
    if (ierr /= 0) call errore('kmesh_dealloc','Error in deallocating nnlist',1)
    return

  end subroutine kmesh_dealloc

  !==================================================================
  subroutine kmesh_supercell_sort
    !==================================================================
    !
    !! We look for kpoint neighbours in a large supercell of reciprocal
    !! unit cells. Done sequentially this is very slow.
    !! Here we order the cells by the distance from the origin.
    !! Doing the search in this order gives a dramatic speed up
    !
    !==================================================================
    implicit none
    integer :: counter, l, m, n, loop
    integer :: lmn_cp(3, (2*nsupcell + 1)**3), indx(1)
    real(kind=dp) :: pos(3)
    real(kind=dp) :: dist((2*nsupcell + 1)**3)
    real(kind=dp) :: dist_cp((2*nsupcell + 1)**3)

    counter = 1
    lmn(:, counter) = 0
    dist(counter) = 0.0_dp
    do l = -nsupcell, nsupcell
      do m = -nsupcell, nsupcell
        do n = -nsupcell, nsupcell
          if (l == 0 .and. m == 0 .and. n == 0) cycle
          counter = counter + 1
          lmn(1, counter) = l; lmn(2, counter) = m; lmn(3, counter) = n
          pos = matmul(lmn(:, counter), bg)
          dist(counter) = sqrt(dot_product(pos, pos))
        end do
      end do
    end do

    do loop = (2*nsupcell + 1)**3, 1, -1
      indx = internal_maxloc(dist)
      dist_cp(loop) = dist(indx(1))
      lmn_cp(:, loop) = lmn(:, indx(1))
      dist(indx(1)) = -1.0_dp
    end do

    lmn = lmn_cp
    dist = dist_cp

  end subroutine kmesh_supercell_sort

  !=============================================================
  subroutine kmesh_get_bvectors(multi, kpt, shell_dist, bvector)
    !=============================================================
    !
    !! Returns the b-vectors for a given shell and kpoint.
    !
    !=============================================================
    implicit none

    integer, intent(in) :: multi   ! the number of kpoints in the shell
    integer, intent(in) :: kpt     ! which kpt is our 'origin'
    real(kind=dp), intent(in) :: shell_dist ! the bvectors
    real(kind=dp), intent(out) :: bvector(3, multi) ! the bvectors
    integer :: loop, nkp2, num_bvec
    real(kind=dp) :: dist, vkpp2(3), vkpp(3)

    bvector = 0.0_dp
    num_bvec = 0
    ok: do loop = 1, (2*nsupcell + 1)**3
      vkpp2 = matmul(lmn(:, loop), bg)
      do nkp2 = 1, num_kpts
        vkpp = vkpp2 + kpt_cart(:, nkp2)
        dist = sqrt((kpt_cart(1, kpt) - vkpp(1))**2 &
                    + (kpt_cart(2, kpt) - vkpp(2))**2 + (kpt_cart(3, kpt) - vkpp(3))**2)
        if ((dist .ge. shell_dist*(1.0_dp - kmesh_tol)) .and. dist .le. shell_dist*(1.0_dp + kmesh_tol)) then
          num_bvec = num_bvec + 1
          bvector(:, num_bvec) = vkpp(:) - kpt_cart(:, kpt)
        endif
        !if we have the right number of neighbours we can exit
        if (num_bvec == multi) cycle ok
      enddo
    enddo ok

    if (num_bvec < multi) call errore('kmesh_get_bvector','Not enough bvectors found',1)
    return

  end subroutine kmesh_get_bvectors

  !==========================================================================
  subroutine kmesh_shell_automatic(multi, dnn, bweight)
    !==========================================================================
    !
    !! Find the correct set of shells to satisfy B1
    !!  The stratagy is:
    !!       1) Take the bvectors from the next shell
    !!       2) Reject them if they are parallel to exisiting b vectors
    !!       3) Test to see if we satisfy B1, if not add another shell and repeat
    !
    !==========================================================================

    use constants_epw, only: eps5, eps6
    implicit none

    integer, intent(in) :: multi(search_shells)   ! the number of kpoints in the shell
    real(kind=dp), intent(in) :: dnn(search_shells) ! the bvectors
    real(kind=dp), intent(out) :: bweight(max_shells)
    real(kind=dp), allocatable     :: bvector(:, :, :) ! the bvectors
    real(kind=dp), dimension(:), allocatable :: singv, tmp1, tmp2, tmp3
    real(kind=dp), dimension(:, :), allocatable :: amat, umat, vmat, smat, tmp0
    integer, parameter :: lwork = max_shells*10
    real(kind=dp) :: work(lwork)
    real(kind=dp), parameter :: target(6) = (/1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
    logical :: b1sat, lpar
    integer :: loop_i, loop_j, loop_bn, loop_b, loop_s, info, cur_shell, ierr
    real(kind=dp) :: delta
    integer :: loop, shell


    allocate (bvector(3, maxval(multi), max_shells), stat=ierr)
    if (ierr /= 0) call errore('kmesh_shell_automatic','Error allocating bvector',1)
    
    bvector = 0.0_dp; bweight = 0.0_dp
    write (stdout, '(1x,a)') '| The b-vectors are chosen automatically                                     |'
    b1sat = .false.
    do shell = 1, search_shells
      cur_shell = num_shells + 1
      ! get the b vectors for the new shell
      call kmesh_get_bvectors(multi(shell), 1, dnn(shell), bvector(:, 1:multi(shell), cur_shell))
      write (stdout, '(1x,a8,1x,I2,a14,1x,I2,49x,a)') '| Shell:', shell, ' Multiplicity:', multi(shell), '|'
      do loop = 1, multi(shell)
        write (stdout, '(1x,a10,I2,1x,a1,4x,3f12.6,5x,a9,9x,a)') '| b-vector ', loop, ':', &
          bvector(:, loop, cur_shell), 'Bohr^-1', '|'
      end do
      ! We check that the new shell is not parrallel to an existing shell (cosine=1)
      lpar = .false.
      if (num_shells > 0) then
        do loop_bn = 1, multi(shell)
          do loop_s = 1, num_shells
            do loop_b = 1, multi(shell_list(loop_s))
              delta = dot_product(bvector(:, loop_bn, cur_shell), bvector(:, loop_b, loop_s))/ &
                      sqrt(dot_product(bvector(:, loop_bn, cur_shell), bvector(:, loop_bn, cur_shell))* &
                           dot_product(bvector(:, loop_b, loop_s), bvector(:, loop_b, loop_s)))
              if (abs(abs(delta) - 1.0_dp) < eps6) lpar = .true.
            end do
          end do
        end do
      end if

      if (lpar) then
        ! if (iprint >= 3 .and. on_root) then
          write (stdout, '(1x,a)') '| This shell is linearly dependent on existing shells: Trying next shell     |'
        ! end if
        cycle
      end if

      num_shells = num_shells + 1
      shell_list(num_shells) = shell
      allocate (tmp0(max_shells, max_shells), stat=ierr)
      if (ierr /= 0) call errore('kmesh_shell_automatic',&
        'Error allocating amat', 1)
      allocate (tmp1(max_shells), stat=ierr)
      if (ierr /= 0) call errore('kmesh_shell_automatic',& 
        'Error allocating amat', 1)
      allocate (tmp2(num_shells), stat=ierr)
      if (ierr /= 0) call errore('kmesh_shell_automatic', &
        'Error allocating amat', 1)
      allocate (tmp3(num_shells), stat=ierr)
      if (ierr /= 0) call errore('kmesh_shell_automatic', &
        'Error allocating amat', 1)
      allocate (amat(max_shells, num_shells), stat=ierr)
      if (ierr /= 0) call errore('kmesh_shell_automatic',&
        'Error allocating amat', 1)
      allocate (umat(max_shells, max_shells), stat=ierr)
      if (ierr /= 0) call errore('kmesh_shell_automatic',&
        'Error allocating umat', 1)
      allocate (vmat(num_shells, num_shells), stat=ierr)
      if (ierr /= 0) call errore('kmesh_shell_automatic',&
        'Error allocating vmat', 1)
      allocate (smat(num_shells, max_shells), stat=ierr)
      if (ierr /= 0) call errore('kmesh_shell_automatic',&
        'Error allocating smat', 1)
      allocate (singv(num_shells), stat=ierr)
      if (ierr /= 0) call errore('kmesh_shell_automatic',&
        'Error allocating singv', 1)
      amat(:, :) = 0.0_dp; umat(:, :) = 0.0_dp; vmat(:, :) = 0.0_dp; smat(:, :) = 0.0_dp; singv(:) = 0.0_dp
      amat = 0.0_dp
      do loop_s = 1, num_shells
        do loop_b = 1, multi(shell_list(loop_s))
          amat(1, loop_s) = amat(1, loop_s) + bvector(1, loop_b, loop_s)*bvector(1, loop_b, loop_s)
          amat(2, loop_s) = amat(2, loop_s) + bvector(2, loop_b, loop_s)*bvector(2, loop_b, loop_s)
          amat(3, loop_s) = amat(3, loop_s) + bvector(3, loop_b, loop_s)*bvector(3, loop_b, loop_s)
          amat(4, loop_s) = amat(4, loop_s) + bvector(1, loop_b, loop_s)*bvector(2, loop_b, loop_s)
          amat(5, loop_s) = amat(5, loop_s) + bvector(2, loop_b, loop_s)*bvector(3, loop_b, loop_s)
          amat(6, loop_s) = amat(6, loop_s) + bvector(3, loop_b, loop_s)*bvector(1, loop_b, loop_s)
        end do
      end do

      info = 0
      call dgesvd('A', 'A', max_shells, num_shells, amat, max_shells, singv, umat, max_shells, vmat, num_shells, work, lwork, info)
      if (info < 0) then
        write (stdout, '(1x,a,1x,I1,1x,a)') 'kmesh_shell_automatic: Argument', abs(info), 'of dgesvd is incorrect'
        call errore('kmesh_shell_automatic',' Problem with Singular Value Decomposition',1)
      else if (info > 0) then
        call errore('kmesh_shell_automatic',' Singular Value Decomposition did not converge',1)
      end if

      if (any(abs(singv) < eps5)) then
        if (num_shells == 1) then
          call errore('kmesh_shell_automatic',' Singular Value Decomposition has found a very small singular value',1)
        else
          write (stdout, '(1x,a)') '| SVD found small singular value, Rejecting this shell and trying the next   |'
          b1sat = .false.
          num_shells = num_shells - 1
          goto 200
        end if
      end if

      smat = 0.0_dp
      do loop_s = 1, num_shells
        smat(loop_s, loop_s) = 1.0_dp/singv(loop_s)
      end do

      ! S. Ponce: The following below is correct but had to be unpacked because of PGI-15
      ! bweight(1:num_shells)=matmul(transpose(vmat),matmul(smat,matmul(transpose(umat),target)))
      tmp0 = transpose(umat)
      tmp1 = matmul(tmp0, target)
      tmp2 = matmul(smat, tmp1)
      tmp3 = matmul(transpose(vmat), tmp2)
      bweight(1:num_shells) = tmp3
      do loop_s = 1, num_shells
        write (stdout, '(1x,a,I2,a,f12.7,5x,a8,36x,a)') '| Shell: ', loop_s, &
          ' w_b ', bweight(loop_s), '(Bohr^2)', '|'
      end do
      !check b1
      b1sat = .true.
      do loop_i = 1, 3
        do loop_j = loop_i, 3
          delta = 0.0_dp
          do loop_s = 1, num_shells
            do loop_b = 1, multi(shell_list(loop_s))
              delta = delta + bweight(loop_s)*bvector(loop_i, loop_b, loop_s)*bvector(loop_j, loop_b, loop_s)
            end do
          end do
          if (loop_i == loop_j) then
            if (abs(delta - 1.0_dp) > kmesh_tol) b1sat = .false.
          end if
          if (loop_i /= loop_j) then
            if (abs(delta) > kmesh_tol) b1sat = .false.
          end if
        end do
      end do

      if (.not. b1sat) then
        if (shell < search_shells) then
          write (stdout, '(1x,a,24x,a1)') '| B1 condition is not satisfied: Adding another shell', '|'
        elseif (shell == search_shells) then
          write (stdout, *) ' '
          write (stdout, '(1x,a,i3,a)') 'Unable to satisfy B1 with any of the first ', search_shells, ' shells'
          write (stdout, '(1x,a)') 'Check that you have specified your unit cell to a high precision'
          write (stdout, '(1x,a)') 'Low precision might cause a loss of symmetry.'
          write (stdout, '(1x,a)') ' '
          write (stdout, '(1x,a)') 'If your cell is very long, or you have an irregular MP grid'
          write (stdout, '(1x,a)') 'Try increasing the parameter search_shells in the win file (default=30)'
          write (stdout, *) ' '
          call errore('kmesh_get_automatic',' ',1)
        end if
      end if

200   continue
      deallocate (tmp0, stat=ierr)
      if (ierr /= 0) call errore('kmesh_shell_automatic','Error deallocating amat',1)
      deallocate (tmp1, stat=ierr)
      if (ierr /= 0) call errore('kmesh_shell_automatic','Error deallocating amat',1)
      deallocate (tmp2, stat=ierr)
      if (ierr /= 0) call errore('kmesh_shell_automatic','Error deallocating amat',1)
      deallocate (tmp3, stat=ierr)
      if (ierr /= 0) call errore('kmesh_shell_automatic','Error deallocating amat',1)
      deallocate (amat, stat=ierr)
      if (ierr /= 0) call errore('kmesh_shell_automatic','Error deallocating amat',1)
      deallocate (umat, stat=ierr)
      if (ierr /= 0) call errore('kmesh_shell_automatic','Error deallocating umat',1)
      deallocate (vmat, stat=ierr)
      if (ierr /= 0) call errore('kmesh_shell_automatic','Error deallocating vmat',1)
      deallocate (smat, stat=ierr)
      if (ierr /= 0) call errore('kmesh_shell_automatic','Error deallocating smat',1)
      deallocate (singv, stat=ierr)
      if (ierr /= 0) call errore('kmesh_shell_automatic','Error deallocating singv',1)
      if (b1sat) exit
    end do

    if (.not. b1sat) then
      write (stdout, *) ' '
      write (stdout, '(1x,a,i3,a)') 'Unable to satisfy B1 with any of the first ', search_shells, ' shells'
      write (stdout, '(1x,a)') 'Your cell might be very long, or you may have an irregular MP grid'
      write (stdout, '(1x,a)') 'Try increasing the parameter search_shells in the win file (default=12)'
      write (stdout, *) ' '
      call errore('kmesh_get_automatic',' ',1)
    end if

    return

  end subroutine kmesh_shell_automatic

  !================================================================
  subroutine kmesh_shell_fixed(multi, dnn, bweight)
    !================================================================
    !
    !!  Find the B1 weights for a set of shells specified by the user
    !
    !================================================================

    !use constants_epw, only: eps7
    ! use w90_io, only: io_error, stdout, io_stopwatch
    implicit none

    integer, intent(in) :: multi(search_shells)   ! the number of kpoints in the shell
    real(kind=dp), parameter :: eps7 = 1.0E-7_DP 
    real(kind=dp), intent(in) :: dnn(search_shells) ! the bvectors
    real(kind=dp), intent(out) :: bweight(max_shells)
    real(kind=dp), allocatable     :: bvector(:, :, :)
    real(kind=dp) :: singv(num_shells)
    real(kind=dp) :: amat(max_shells, num_shells)
    real(kind=dp) :: umat(max_shells, max_shells)
    real(kind=dp) :: vmat(num_shells, num_shells)
    real(kind=dp) :: smat(num_shells, max_shells)
    integer, parameter :: lwork = max_shells*10
    real(kind=dp) :: work(lwork)
    real(kind=dp), parameter :: target(6) = (/1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
    logical :: b1sat
    integer :: ierr, loop_i, loop_j, loop_b, loop_s, info
    real(kind=dp) :: delta
    integer :: loop, shell

    allocate (bvector(3, maxval(multi), num_shells), stat=ierr)
    if (ierr /= 0) call errore('kmesh_shell_fixed','Error allocating bvector',1)
    bvector = 0.0_dp; bweight = 0.0_dp
    amat = 0.0_dp; umat = 0.0_dp; vmat = 0.0_dp; smat = 0.0_dp; singv = 0.0_dp

    write (stdout, '(1x,a)') '| The b-vectors are set in the win file                                      |'

    do shell = 1, num_shells
      ! get the b vectors for this shell
      call kmesh_get_bvectors(multi(shell_list(shell)), 1, dnn(shell_list(shell)), &
                              bvector(:, 1:multi(shell_list(shell)), shell))
    end do


    do shell = 1, num_shells
      write (stdout, '(1x,a8,1x,I2,a14,1x,I2,49x,a)') '| Shell:', shell, ' Multiplicity:', multi(shell_list(shell)), '|'
      do loop = 1, multi(shell_list(shell))
        write (stdout, '(1x,a10,I2,1x,a1,4x,3f12.6,5x,a9,9x,a)') '| b-vector ', loop, ':', &
          bvector(:, loop, shell), '(Bohr^-1)', '|'
      end do
    end do


    do loop_s = 1, num_shells
      do loop_b = 1, multi(shell_list(loop_s))
        amat(1, loop_s) = amat(1, loop_s) + bvector(1, loop_b, loop_s)*bvector(1, loop_b, loop_s)
        amat(2, loop_s) = amat(2, loop_s) + bvector(2, loop_b, loop_s)*bvector(2, loop_b, loop_s)
        amat(3, loop_s) = amat(3, loop_s) + bvector(3, loop_b, loop_s)*bvector(3, loop_b, loop_s)
        amat(4, loop_s) = amat(4, loop_s) + bvector(1, loop_b, loop_s)*bvector(2, loop_b, loop_s)
        amat(5, loop_s) = amat(5, loop_s) + bvector(2, loop_b, loop_s)*bvector(3, loop_b, loop_s)
        amat(6, loop_s) = amat(6, loop_s) + bvector(3, loop_b, loop_s)*bvector(1, loop_b, loop_s)
      end do
    end do

    info = 0
    call dgesvd('A', 'A', max_shells, num_shells, amat, max_shells, singv, umat, max_shells, vmat, num_shells, work, lwork, info)
    if (info < 0) then
      write (stdout, '(1x,a,1x,I1,1x,a)') 'kmesh_shell_fixed: Argument', abs(info), 'of dgesvd is incorrect'
      call errore('kmesh_shell_fixed','Problem with Singular Value Decomposition',1)
    else if (info > 0) then
      call errore('kmesh_shell_fixed','Singular Value Decomposition did not converge',1)
    end if

    if (any(abs(singv) < eps7)) &
      call errore('kmesh_shell_fixed','Singular Value Decomposition has found a very small singular value',1)

    smat = 0.0_dp
    do loop_s = 1, num_shells
      smat(loop_s, loop_s) = 1/singv(loop_s)
    end do

    bweight(1:num_shells) = matmul(transpose(vmat), matmul(smat, matmul(transpose(umat), target)))

    do loop_s = 1, num_shells
!     write(stdout,'(1x,a,I2,a,f12.7,49x,a)') '| Shell: ',loop_s,' w_b ', bweight(loop_s),'|'
      write (stdout, '(1x,a,I2,a,f12.7,5x,a8,36x,a)') '| Shell: ', loop_s, &
        ' w_b ', bweight(loop_s), '(Bohr^2)', '|'
    end do


    !check b1
    b1sat = .true.
!    if (.not. skip_B1_tests) then
      do loop_i = 1, 3
        do loop_j = loop_i, 3
          delta = 0.0_dp
          do loop_s = 1, num_shells
            do loop_b = 1, multi(shell_list(loop_s))
              delta = delta + bweight(loop_s)*bvector(loop_i, loop_b, loop_s)*bvector(loop_j, loop_b, loop_s)
            end do
          end do
          if (loop_i == loop_j) then
            if (abs(delta - 1.0_dp) > kmesh_tol) b1sat = .false.
          end if
          if (loop_i /= loop_j) then
            if (abs(delta) > kmesh_tol) b1sat = .false.
          end if
        end do
      end do
!   end if

    if (.not. b1sat) call errore('kmesh_shell_fixed','B1 condition not satisfied',1)

    return

  end subroutine kmesh_shell_fixed

    !=============================================================
  subroutine kmesh_gradient(f_k, df_k)
    !=============================================================
    !
    !! Returns the k-derivative of a scalar function defined on the
    !! k grid.
    !
    !=============================================================
    implicit none

    integer, intent(in) :: f_k(num_kpts)   
    integer, intent(out) :: df_k(3,num_kpts)
    integer :: nkp, nn, ipol, nkp1

    df_k = 0.0_dp
    do nkp = lower_bnd, upper_bnd
      do nn = 1, nntot
        nkp1 = nnlist(nkp, nn)
        do ipol = 1,3
          df_k(ipol,nkp) = df_k(ipol,nkp)+ wb(nn)*bk(ipol, nn, nkp)*(f_k(nkp1)-f_k(nkp))
        enddo ! ipol
      enddo
    enddo
    call mp_sum(df_k, inter_pool_comm)
    return
  end subroutine kmesh_gradient

  !=================================
  function internal_maxloc(dist)
    !=================================
    !
    !!  A reproducible maxloc function
    !!  so b-vectors come in the same
    !!  order each time
    !=================================

    use constants_epw, only: eps8
    implicit none

    real(kind=dp), intent(in)  :: dist((2*nsupcell + 1)**3)
    !! Distances from the origin of the unit cells in the supercell.
    integer                    :: internal_maxloc
    integer       :: guess(1), loop, counter
    integer       :: list((2*nsupcell + 1)**3)

    list = 0
    counter = 1

    guess = maxloc(dist)
    list(1) = guess(1)
    ! look for any degenerate values
    do loop = 1, (2*nsupcell + 1)**3
      if (loop == guess(1)) cycle
      if (abs(dist(loop) - dist(guess(1))) < eps8) then
        counter = counter + 1
        list(counter) = loop
      endif
    end do
    ! and always return the lowest index
    internal_maxloc = minval(list(1:counter))

  end function internal_maxloc


  !===================================================================
  subroutine utility_compar(a, b, ifpos, ifneg)
    !==================================================================!
    !                                                                  !
    !! Compares two vectors
    !                                                                  !
    !===================================================================
    use constants_epw, only: eps8
    implicit none

    real(kind=dp), intent(in) :: a(3)
    real(kind=dp), intent(in) :: b(3)
    integer, intent(out) :: ifpos, ifneg
    real(kind=dp) :: rrp, rrm

    rrp = (a(1) - b(1))**2 + (a(2) - b(2))**2 + (a(3) - b(3))**2
    rrm = (a(1) + b(1))**2 + (a(2) + b(2))**2 + (a(3) + b(3))**2
    ifpos = 0
    if (abs(rrp) .lt. eps8) ifpos = 1
    ifneg = 0
    if (abs(rrm) .lt. eps8) ifneg = 1

    return

  end subroutine utility_compar

end module w90_kmesh_epw
