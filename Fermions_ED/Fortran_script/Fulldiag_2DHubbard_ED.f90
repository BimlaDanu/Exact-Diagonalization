program hubbard_ed_2d
  implicit none
  !------------------ PARAMETERS ------------------
  integer, parameter ::  numv = 1!number of traget eigenavlues in lapack
  integer, parameter :: Lx = 3, Ly = 2         ! lattice dims (Lx columns, Ly rows)
  integer, parameter :: n = Lx * Ly             ! total sites
  integer, parameter :: N_up = int(n/2)-1, N_dn = int(n/2)-0!N_dn = n-N_up!int(n/2)      ! particle numbers
  real*8, parameter  :: t = 1.d0               ! NN hopping
  real*8, parameter  :: tprime = -0.3d0          ! NNN hopping (diagonal)
  real*8, parameter  :: U = 8.0d0               ! on-site interaction
  real*8, parameter  :: mu = 0.d0              ! chemical potential
  !logical, parameter :: PBC_x = .true., PBC_y = .true.  ! boundary conditions
  logical, parameter :: PBC_x = .false., PBC_y = .false.  ! boundary conditions
  !----------------------------------------------------------------

  integer, parameter :: two_n = 2**(2*n)        ! maximum combined keys (4^n)
  integer :: up_bits, dn_bits, key
  integer :: countu, countd
  integer :: nstates, idx
  integer, allocatable :: state_up(:), state_dn(:)
  integer, allocatable :: bckwindx(:)           ! map combined key -> basis index 
  real*8, allocatable :: H(:,:)
  real*8, allocatable :: eig(:), work(:)
  integer :: i, j, k, site, neigh_count, info, lwork
  integer :: neighbors(8)                       ! up to 4 NN + 4 NNN 
  integer :: nn_list(4), nnn_list(4)
  integer :: nnc, nnnc
  integer :: dest
  integer :: sgn1, sgn2
  integer :: up_new1, up_new2, dn_new1, dn_new2
  integer :: combined_key, jidx
  integer :: sx, sy, nsites
  real*8 :: val
  real*8, allocatable :: psi(:)
  real*8, allocatable :: Szz(:,:), Sxx(:,:), dens(:), dd(:,:), doubl(:)



  !------------------ Build basis mapping ------------------
  print *, "2D Hubbard ED (t, t')"
  print *, "Lx=", Lx, "Ly=", Ly, "n=", n
  print *, "PBC_x=", PBC_x, "PBC_y=", PBC_y
  print *, "N_up=", N_up, "N_dn=", N_dn, "t=", t, "t'=", tprime, "U=", U, "mu=", mu

  ! Count basis states (nstates = C(n,N_up) * C(n,N_dn))
  nstates = 0
  do up_bits = 0, 2**n - 1
    countu = bit_count(up_bits)
    if (countu /= N_up) cycle
    do dn_bits = 0, 2**n - 1
      countd = bit_count(dn_bits)
      if (countd /= N_dn) cycle
      nstates = nstates + 1
    end do
  end do

!  nstates = 0
!  do up_bits = 0, 2**n - 1
!    countu = 0
!    do site = 0, n - 1
!      if (btest(up_bits, site)) countu = countu + 1
!    end do
!    if (countu /= N_up) cycle
!    do dn_bits = 0, 2**n - 1
!      countd = 0
!      do site = 0, n - 1
!        if (btest(dn_bits, site)) countd = countd + 1
!      end do
!      if (countd /= N_dn) cycle
!      nstates = nstates + 1
!    end do
!  end do
!  

  if (nstates == 0) then
    print *, "No basis states for N_up=", N_up, " N_dn=", N_dn
    stop
  end if
  print *, "Basis size (nstates) = ", nstates

  allocate(state_up(nstates))
  allocate(state_dn(nstates))
  allocate(bckwindx(0:two_n-1))
  bckwindx = 0

  idx = 0
  do up_bits = 0, 2**n - 1
    if (bit_count(up_bits) /= N_up) cycle
    do dn_bits = 0, 2**n - 1
      if (bit_count(dn_bits) /= N_dn) cycle
      idx = idx + 1
      state_up(idx) = up_bits
      state_dn(idx) = dn_bits
      key = up_bits + ishft(dn_bits, n)      ! combined key: low n bits = up, high n bits = dn
      bckwindx(key) = idx
    end do
  end do

  ! Allocate Hamiltonian
  allocate(H(nstates, nstates))
  H = 0.d0

  !------------------Hamiltonian matrix ------------------
  do k = 1, nstates
    up_bits = state_up(k)
    dn_bits = state_dn(k)

    ! Diagonal terms: U and chemical potential
    do site = 0, n - 1
      if ( btest(up_bits, site) .and. btest(dn_bits, site) ) H(k,k) = H(k,k) + U
      if ( btest(up_bits, site) ) H(k,k) = H(k,k) - mu
      if ( btest(dn_bits, site) ) H(k,k) = H(k,k) - mu
    end do

    ! For each site, get NN and NNN neighbor lists
    do site = 0, n - 1
      call get_nn_nnn(site, Lx, Ly, PBC_x, PBC_y, nn_list, nnc, nnn_list, nnnc)

      ! --- NN hopping (t) ---
      do i = 1, nnc
        dest = nn_list(i)
        !print*, 'dest,nn_list(i)=',dest,nn_list(i)
        ! spin-up: c^dag_dest c_site  (annihilate at site then create at dest)
        call try_hop_add(up_bits, dn_bits, n, site, dest, t, .true., bckwindx, H, k)
        ! spin-down
        call try_hop_add(dn_bits, up_bits, n, site, dest, t, .false., bckwindx, H, k)
      end do

      ! --- NNN hopping (tprime) ---
      do i = 1, nnnc
        dest = nnn_list(i)
        !print*, 'dest,nnn_list(i)=',dest,nnn_list(i)
        call try_hop_add(up_bits, dn_bits, n, site, dest, tprime, .true., bckwindx, H, k)
        call try_hop_add(dn_bits, up_bits, n, site, dest, tprime, .false., bckwindx, H, k)
      end do
    end do
  end do

  ! Symmetrize H to ensure 
  do i = 1, nstates
    do j = i+1, nstates
      val = 0.5d0 * (H(i,j) + H(j,i))
      H(i,j) = val
      H(j,i) = val
    end do
  end do
  

  !------------------ Diagonalize H  ------------------
  allocate(eig(nstates))
!  lwork = max(1, 3*nstates - 1)
!  allocate(work(lwork))
!  call dsyev('V', 'U', nstates, H, nstates, eig, work, lwork, info)
!  if (info /= 0) then
!    print *, "DSYEV failed with info=", info
!    stop
!  end if
!  print *, "Lowest eigenvalues (first min(10,nstates)):"
!  do i = 1, min(10, nstates)
!    print '(I3,2X,F14.8)', i, eig(i)
!  end do
  allocate(psi(nstates))
  !psi = H(:,1)  ! ground state wavefunction


  !call sym(H, eig, nstates)
  !psi = H(:,1)  ! ground state wavefunction
  !call sym_few(H, eig, nstates)
  call sym_few(H, eig, psi, nstates, numv)
  print*, ('eig(i) =',eig(i), i = 1, numv) 
  print *, "Ground state energy:", eig(1)

  !-------------------------------------------------------
  !correlations
  !-------------------------------------------------------
  nsites = n
  allocate(Szz(nsites,nsites), Sxx(nsites,nsites))
  allocate(dens(nsites), dd(nsites,nsites), doubl(nsites))
  call hubbard_szz_sxx(nsites,nstates,state_up,state_dn,psi,Szz,Sxx)
  call hubbard_density(nsites,nstates,state_up,state_dn,psi,dens,dd,doubl)

  !-------------------------------------------------------
  ! Write output to files
  !-------------------------------------------------------
  call write_correlation_file("2d_Szz.dat",Szz,nsites)
  call write_correlation_file("2d_Sxx.dat",Sxx,nsites)
  call write_correlation_file("2d_density_density.dat",dd,nsites)
  call write_density_file("2d_density_double.dat",dens,doubl,nsites)
  
  

  deallocate(state_up); deallocate(state_dn); deallocate(bckwindx)
  deallocate(eig); deallocate(H)!; deallocate(work)

contains

  !------------------------------
  ! Count bits set to 1 
  !------------------------------
  integer function bit_count(x)
    integer, intent(in) :: x
    integer :: tmp
    tmp = x
    bit_count = 0
    do while (tmp /= 0)
      if (btest(tmp, 0)) bit_count = bit_count + 1
      tmp = ishft(tmp, -1)
    end do
  end function bit_count

  !----------------------------------------------------------------
  ! Build NN and NNN neighbor lists for a 2D lattice 
  ! site = y*Lx + x, x in [0,Lx-1], y in [0,Ly-1]
  ! nn_list returns up to 4 nearest neighbors 
  ! nnn_list returns up to 4 next-nearest (diagonals)
  !----------------------------------------------------------------
  subroutine get_nn_nnn(site, Lx, Ly, PBC_x, PBC_y, nn_list, nnc, nnn_list, nnnc)
    integer, intent(in) :: site, Lx, Ly
    logical, intent(in) :: PBC_x, PBC_y
    integer, intent(out) :: nn_list(4), nnc
    integer, intent(out) :: nnn_list(4), nnnc
    integer :: x, y, xp, yp, s

    x = mod(site, Lx)
    y = site / Lx
    nnc = 0
    nnnc = 0

    ! +x
    if (x + 1 < Lx .or. PBC_x) then
      xp = mod(x + 1, Lx); yp = y
      nnc = nnc + 1
      nn_list(nnc) = yp*Lx + xp
    end if
    ! -x
    if (x - 1 >= 0 .or. PBC_x) then
      xp = mod(x - 1 + Lx, Lx); yp = y
      nnc = nnc + 1
      nn_list(nnc) = yp*Lx + xp
    end if
    ! +y
    if (y + 1 < Ly .or. PBC_y) then
      yp = mod(y + 1, Ly); xp = x
      nnc = nnc + 1
      nn_list(nnc) = yp*Lx + xp
    end if
    ! -y
    if (y - 1 >= 0 .or. PBC_y) then
      yp = mod(y - 1 + Ly, Ly); xp = x
      nnc = nnc + 1
      nn_list(nnc) = yp*Lx + xp
    end if

    ! Diagonal (+1,+1)
    if ( (x+1 < Lx .or. PBC_x) .and. (y+1 < Ly .or. PBC_y) ) then
      xp = mod(x+1, Lx); yp = mod(y+1, Ly)
      nnnc = nnnc + 1
      nnn_list(nnnc) = yp*Lx + xp
    end if
    ! (+1,-1)
    if ( (x+1 < Lx .or. PBC_x) .and. (y-1 >= 0 .or. PBC_y) ) then
      xp = mod(x+1, Lx); yp = mod(y-1+Ly, Ly)
      nnnc = nnnc + 1
      nnn_list(nnnc) = yp*Lx + xp
    end if
    ! (-1,+1)
    if ( (x-1 >= 0 .or. PBC_x) .and. (y+1 < Ly .or. PBC_y) ) then
      xp = mod(x-1+Lx, Lx); yp = mod(y+1, Ly)
      nnnc = nnnc + 1
      nnn_list(nnnc) = yp*Lx + xp
    end if
    ! (-1,-1)
    if ( (x-1 >= 0 .or. PBC_x) .and. (y-1 >= 0 .or. PBC_y) ) then
      xp = mod(x-1+Lx, Lx); yp = mod(y-1+Ly, Ly)
      nnnc = nnnc + 1
      nnn_list(nnnc) = yp*Lx + xp
    end if
  end subroutine get_nn_nnn

  !----------------------------------------------------------------
  ! Try to hop one particle of given spin from 'site' -> 'dest'
  !----------------------------------------------------------------
  subroutine try_hop_add(active_bits, passive_bits, nbits, site, dest, amp, is_up_spin, bckwindx, H, k)
    integer, intent(in) :: active_bits, passive_bits, nbits, site, dest, k
    real*8, intent(in) :: amp
    logical, intent(in) :: is_up_spin
    integer, intent(in) :: bckwindx(0:)
    real*8, intent(inout) :: H(:,:)
    integer :: sgn1, sgn2
    integer :: tmp_bits, new_bits
    integer :: combined_key, jidx

    ! Check occupation at source 
    if (.not. btest(active_bits, site)) return    ! 
    if ( btest(active_bits, dest) ) return        !

    ! sign from annihilating at 'site' on the original bitstring
    sgn1 = fermion_sign(active_bits, site)
    ! state after annihilation
    tmp_bits = ibclr(active_bits, site)
    ! sign for creating at 'dest' on the intermediate state
    sgn2 = fermion_sign(tmp_bits, dest)
    new_bits = ibset(tmp_bits, dest)

    if (is_up_spin) then
      combined_key = new_bits + ishft(passive_bits, nbits)
    else
      combined_key = passive_bits + ishft(new_bits, nbits)
    end if
    jidx = bckwindx(combined_key)
    if (jidx /= 0) then
      H(jidx, k) = H(jidx, k) - amp * dble(sgn1 * sgn2)
    end if
  end subroutine try_hop_add

  !----------------------------------------------------------------
  ! Fermionic sign: +1 if even number of occupied sites with index < pos,
  !                -1 if odd.
  !----------------------------------------------------------------
  integer function fermion_sign(bits, pos)
    integer, intent(in) :: bits, pos
    integer :: p, cnt
    cnt = 0
    if (pos > 0) then
      do p = 0, pos-1
        if (btest(bits, p)) cnt = cnt + 1
      end do
    end if
    if (mod(cnt, 2) == 0) then
      fermion_sign = 1
    else
      fermion_sign = -1
    end if
  end function fermion_sign


  
  
!-----------------------------------------------------------------------
! Find state index in basis
!-----------------------------------------------------------------------
integer function find_state_index(up_bits,dn_bits,state_up,state_dn,nstates)
  implicit none
  integer, intent(in) :: up_bits,dn_bits,nstates
  integer, intent(in) :: state_up(nstates), state_dn(nstates)
  integer :: k
  find_state_index = 0
  do k=1,nstates
     if (state_up(k)==up_bits .and. state_dn(k)==dn_bits) then
        find_state_index = k
        return
     end if
  end do
end function find_state_index  
  
!-----------------------------------------------------------------------
! Spin correlations Szz and Sxx
!-----------------------------------------------------------------------
subroutine hubbard_szz_sxx(nsites,nstates,state_up,state_dn,x,Szz,Sxx)
  implicit none
  integer, intent(in) :: nsites,nstates
  integer, intent(in) :: state_up(nstates), state_dn(nstates)
  real*8, intent(in) :: x(nstates)
  real*8, intent(out) :: Szz(nsites,nsites), Sxx(nsites,nsites)
  integer :: i,j,k,l,tmp_up,tmp_dn,new_up,new_dn,sgn
  integer :: up1,dn1,up2,dn2

  Szz = 0.d0
  Sxx = 0.d0

  ! Szz
  do i=0,nsites-1
     do j=0,nsites-1
        do k=1,nstates
           up1 = ibits(state_up(k),i,1)
           dn1 = ibits(state_dn(k),i,1)
           up2 = ibits(state_up(k),j,1)
           dn2 = ibits(state_dn(k),j,1)
           Szz(i+1,j+1) = Szz(i+1,j+1) + 0.25d0*(up1-dn1)*(up2-dn2)*x(k)**2
        end do
     end do
  end do

  ! Sxx using S^+ S^- + S^- S^+
  do i=0,nsites-1
     do j=0,nsites-1
        do k=1,nstates
           ! S_i^+ S_j^-
           if (btest(state_dn(k),i) .and. .not.btest(state_up(k),i) &
               .and. btest(state_up(k),j) .and. .not.btest(state_dn(k),j)) then
              tmp_up = state_up(k); tmp_dn = state_dn(k)
              new_up = ibset(tmp_up,i); new_up = ibclr(new_up,j)
              new_dn = ibclr(tmp_dn,i); new_dn = ibset(tmp_dn,j)
              sgn = fermion_sign_hubbard(tmp_up,tmp_dn,i,j)
              l = find_state_index(new_up,new_dn,state_up,state_dn,nstates)
              if (l>0) Sxx(i+1,j+1) = Sxx(i+1,j+1) + 0.5d0*sgn*x(k)*x(l)
           end if
           ! S_i^- S_j^+
           if (btest(state_up(k),i) .and. .not.btest(state_dn(k),i) &
               .and. btest(state_dn(k),j) .and. .not.btest(state_up(k),j)) then
              tmp_up = state_up(k); tmp_dn = state_dn(k)
              new_up = ibclr(tmp_up,i); new_up = ibset(new_up,j)
              new_dn = ibset(tmp_dn,i); new_dn = ibclr(new_dn,j)
              sgn = fermion_sign_hubbard(tmp_up,tmp_dn,j,i)
              l = find_state_index(new_up,new_dn,state_up,state_dn,nstates)
              if (l>0) Sxx(i+1,j+1) = Sxx(i+1,j+1) + 0.5d0*sgn*x(k)*x(l)
           end if
        end do
     end do
  end do

end subroutine hubbard_szz_sxx

!-----------------------------------------------------------------------
! Fermion sign for S^+ S^- hopping
!-----------------------------------------------------------------------
integer function fermion_sign_hubbard(up_bits,dn_bits,site_from,site_to)
  implicit none
  integer, intent(in) :: up_bits,dn_bits,site_from,site_to
  integer :: cnt,p
  cnt = 0
  if (site_from < site_to) then
     do p=site_from+1,site_to-1
        if (btest(up_bits,p)) cnt = cnt +1
        if (btest(dn_bits,p)) cnt = cnt +1
     end do
  else if (site_from > site_to) then
     do p=site_to+1,site_from-1
        if (btest(up_bits,p)) cnt = cnt +1
        if (btest(dn_bits,p)) cnt = cnt +1
     end do
  end if
  if (mod(cnt,2)==0) then
     fermion_sign_hubbard = 1
  else
     fermion_sign_hubbard = -1
  end if
end function fermion_sign_hubbard

!-----------------------------------------------------------------------
! Density, density-density, double occupancy
!-----------------------------------------------------------------------
subroutine hubbard_density(nsites,nstates,state_up,state_dn,x,dens,dd_corr,double_occ)
  implicit none
  integer, intent(in) :: nsites,nstates
  integer, intent(in) :: state_up(nstates), state_dn(nstates)
  real*8, intent(in) :: x(nstates)
  real*8, intent(out) :: dens(nsites), dd_corr(nsites,nsites), double_occ(nsites)
  integer :: i,j,k,up,dn

  dens = 0.d0
  dd_corr = 0.d0
  double_occ = 0.d0

  do i=0,nsites-1
     do k=1,nstates
        up = ibits(state_up(k),i,1)
        dn = ibits(state_dn(k),i,1)
        dens(i+1) = dens(i+1) + real(up+dn)*x(k)**2
        double_occ(i+1) = double_occ(i+1) + real(up*dn)*x(k)**2
     end do
  end do

  do i=0,nsites-1
     do j=0,nsites-1
        do k=1,nstates
           dd_corr(i+1,j+1) = dd_corr(i+1,j+1) + real( (ibits(state_up(k),i,1)+ibits(state_dn(k),i,1)) &
                                    *(ibits(state_up(k),j,1)+ibits(state_dn(k),j,1)) )*x(k)**2
        end do
     end do
  end do
end subroutine hubbard_density

!-----------------------------------------------------------------------
! Write spin or density-density correlation files
! Format: i, j, r, value
!-----------------------------------------------------------------------
subroutine write_correlation_file(fname,C,n)
  implicit none
  character(*), intent(in) :: fname
  real*8, intent(in) :: C(n,n)
  integer, intent(in) :: n
  integer :: i,j
  real*8 :: r
  open(unit=11,file=fname,status='replace')
  do i=1,n
     do j=1,n
        r = real(j-i)
       if(abs(r)>0) write(11,'(I4,2X,I4,2X,F6.2,2X,F12.8)') i, j, r, C(i,j)
     end do
  end do
  close(11)
end subroutine write_correlation_file

!-----------------------------------------------------------------------
! Write density and double occupancy file
!-----------------------------------------------------------------------
subroutine write_density_file(fname,dens,double_occ,n)
  implicit none
  character(*), intent(in) :: fname
  real*8, intent(in) :: dens(n), double_occ(n)
  real*8::  sum1, sum2
  integer, intent(in) :: n
  integer :: i
  open(unit=121,file='Total_den_docc_2d.dat',status='replace')
  open(unit=12,file=fname,status='replace')
  sum1 = 0.d0
  sum2 = 0.d0
  do i=1,n
  	  sum1 =  sum1  +  dens(i)
  	  sum2 =  sum2  +  double_occ(i)
     write(12,'(I4,2X,F12.8,2X,F12.8)') i, dens(i), double_occ(i)
  end do
  write(121,'(2X,F12.8,2X,F12.8,2X,F12.8,2X,F12.8)') sum1, sum2, sum1/real(n), sum2/real(n)
  print*, ' sum1, sum2, sum1/real(n), sum2/real(n)=', sum1, sum2, sum1/real(n), sum2/real(n)
  close(12)
  close(121)
end subroutine write_density_file


subroutine sym_few(H, eig, psi0, n, numv)
    implicit none
    integer, intent(in) :: n, numv
    real*8, intent(inout) :: H(n,n)          ! symmetric matrix (will be modified)
    !real*8, intent(out) :: eig               ! ground-state energy
    real*8, intent(out) :: eig(n)
    real*8, intent(out) :: psi0(n)           ! ground-state eigenvector

    ! local variables
    integer :: lda, ldz, il, iu, m, info
    integer, allocatable :: isuppz(:), iwork(:)
    real*8, allocatable :: w(:), z(:,:), work(:)
    real*8 :: abstol
    character(1) :: jobz, range, uplo
    integer :: lwork, liwork

    ! Parameters for DSYEVR
    jobz = 'V'       ! Compute eigenvalues & eigenvectors
    range = 'I'      ! Select eigenvalues by index
    uplo = 'U'       ! Upper triangle of H is stored
    il = 1           ! Smallest eigenvalue
    !iu = 1           ! Only the first (lowest) eigenpair
    iu = min(n, numv)  ! Index of the largest eigenvalue to compute (only 5 eigenpairs)
    lda = n
    ldz = n
    abstol = 1.0d-10

    ! Allocate arrays
    allocate(w(n), z(ldz, iu), isuppz(2*n))

    ! Workspace query
    lwork = -1
    liwork = -1
    allocate(work(1), iwork(1))
    call dsyevr(jobz, range, uplo, n, H, lda, 0.0d0, 0.0d0, il, iu, abstol, &
                m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info)

    if (info /= 0) then
        print *, "Error in DSYEVR workspace query, info =", info
        stop
    endif

    ! Allocate optimal workspace
    lwork = int(work(1))
    liwork = iwork(1)
    deallocate(work, iwork)
    allocate(work(lwork), iwork(liwork))

    ! Actual diagonalization (find lowest eigenpair)
    call dsyevr(jobz, range, uplo, n, H, lda, 0.0d0, 0.0d0, il, iu, abstol, &
                m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info)

    if (info /= 0) then
        print *, "DSYEVR failed, info =", info
        stop
    endif

    ! Return results
    !eig = w(1)           ! ground-state energy
    eig(1:m) = w(1:m)
    psi0 = z(:,1)        ! ground-state vector

    deallocate(work, iwork, w, z, isuppz)
end subroutine sym_few


!For few eigen values
 subroutine sym_few1(H, eig, n)
        implicit none
        integer, intent(in) :: n
        real*8, intent(inout) :: H(n,n)
        real*8, intent(out) :: eig(n)

        integer :: lda, ldz, il, iu, m, info
        integer, allocatable :: isuppz(:), iwork(:)
        double precision, allocatable :: w(:), z(:,:), work(:)
        double precision :: abstol
        character(1) :: jobz, range, uplo
        integer :: lwork, liwork

        ! Set parameters for DSYEVR
        jobz = 'V'    ! Compute eigenvalues and eigenvectors
        range = 'I'   ! 'I' for a subset by index
        uplo = 'U'    ! 'U' for upper triangle of H
        il = 1        ! Index of the smallest eigenvalue to compute
        iu = min(n, 4)  ! Index of the largest eigenvalue to compute (only 5 eigenpairs)
        lda = n
        ldz = n
        abstol = 1.0d-10   ! Tolerance for convergence

        ! Allocate memory for eigenvalues, eigenvectors, and supporting arrays
        allocate(w(n), z(ldz,n), isuppz(2*n))

        ! Workspace query to determine optimal size
        lwork = -1
        liwork = -1
        allocate(work(1))
        allocate(iwork(1))

        ! Query optimal workspace sizes
        call dsyevr(jobz, range, uplo, n, H, lda, 0.0d0, 0.0d0, il, iu, abstol, &
                    m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info)

        ! Check for errors in the query
        if (info /= 0) then
            print *, "Error in DSYEVR workspace query, info =", info
            stop
        endif

        ! Allocate workspace based
        lwork = int(work(1))
        liwork = iwork(1)
        deallocate(work, iwork)
        allocate(work(lwork), iwork(liwork))

        ! Compute eigenvalues and eigenvectors
        call dsyevr(jobz, range, uplo, n, H, lda, 0.0d0, 0.0d0, il, iu, abstol, &
                    m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info)

        ! Check for success
        if (info /= 0) then
            print *, "DSYEVR failed, info =", info
            stop
        endif

        ! Copy the eigenvalues to the output array
        eig(1:m) = w(1:m)

        ! Deallocate allocated memory
        deallocate(work, iwork, w, z, isuppz)
    end subroutine sym_few1

! for full digonalization
	subroutine sym(H,eig,n)
	    implicit none
	    integer n
	    integer lwork,inf
	    real*8::  H(n,n),eig(n),work(max(1,3*n-1))
	    lwork=3*n-1
	    call dsyev('V','U',n,H,n,eig,work,lwork,inf)
	end

end program hubbard_ed_2d
