!-----------------------------------------------------------------------
! 2D Hubbard ED with Lanczos ground-state 
! - Builds basis, implements H * v, Lanczos to obtain GS
! - Returns ground-state energy and ground-state vector psi
!-----------------------------------------------------------------------
program hubbard_ed_2d
  implicit none
  !------------------ PARAMETERS ------------------
  integer, parameter ::  numv = 1      ! number of target eigenvalues (unused here)
  integer, parameter :: Lx = 6, Ly = 2 ! lattice dims (Lx columns, Ly rows)
  integer, parameter :: n = Lx * Ly    ! total sites
  integer, parameter :: N_up = int(n/2)-0, N_dn = int(n/2)-0 ! particle numbers
  real*8, parameter  :: t = 1.d0       ! NN hopping
  real*8, parameter  :: tprime = 0.d0  ! NNN hopping
  real*8, parameter  :: U = 8.0d0      ! on-site interaction
  real*8, parameter  :: mu = 0.d0      ! chemical potential
  logical, parameter :: PBC_x = .false., PBC_y = .false.
  !----------------------------------------------------------------

  integer, parameter :: two_n = 2**(2*n)        ! maximum combined keys (4^n)
  integer :: up_bits, dn_bits, key
  integer :: countu, countd
  integer :: nstates, idx
  integer, allocatable :: state_up(:), state_dn(:)
  integer, allocatable :: bckwindx(:)           ! map combined key -> basis index
  real*8, allocatable :: H(:,:)
  real*8, allocatable :: eig(:), work(:)
   real*8, allocatable ::  gs_energy(:), psi(:)
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
  real*8, allocatable :: Szz(:,:), Sxx(:,:), dens(:), dd(:,:), doubl(:)

  print *, "2D Hubbard ED (t, t')"
  print *, "Lx=", Lx, "Ly=", Ly, "n=", n
  print *, "PBC_x=", PBC_x, "PBC_y=", PBC_y
  print *, "N_up=", N_up, "N_dn=", N_dn, "t=", t, "t'=", tprime, "U=", U, "mu=", mu

  !------------------ Build basis mapping ------------------
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

  allocate(psi(nstates))
  allocate(gs_energy(1))

  ! Run Lanczos ground-state finder
  call lanczos_ground(state_up, state_dn, bckwindx, nstates, 200, 1.d-8, gs_energy(numv), psi,numv)
  print *, 'Ground state energy =', gs_energy(1)
!  do i = 1, min(10, nstates)
!    print '(I6,2X,F12.8)', i, psi(i)
!  end do
!  
!  
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
  call write_correlation_file("2d_Szz_Lanc.dat",Szz,nsites)
  call write_correlation_file("2d_Sxx_Lanc.dat",Sxx,nsites)
  call write_correlation_file("2d_density_density_Lanc.dat",dd,nsites)
  call write_density_file("2d_density_double_Lanc.dat",dens,doubl,nsites)
  
  !clean
  deallocate(state_up); deallocate(state_dn); deallocate(bckwindx)
  deallocate(psi);!deallocate(eig); deallocate(H)

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
  ! Build NN and NNN neighbor lists for a 2D lattice (row-major)
  ! site = y*Lx + x, x in [0,Lx-1], y in [0,Ly-1]
  !----------------------------------------------------------------
  subroutine get_nn_nnn(site, Lx, Ly, PBC_x, PBC_y, nn_list, nnc, nnn_list, nnnc)
    integer, intent(in) :: site, Lx, Ly
    logical, intent(in) :: PBC_x, PBC_y
    integer, intent(out) :: nn_list(4), nnc
    integer, intent(out) :: nnn_list(4), nnnc
    integer :: x, y, xp, yp

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

    ! Diagonals
    if ( (x+1 < Lx .or. PBC_x) .and. (y+1 < Ly .or. PBC_y) ) then
      xp = mod(x+1, Lx); yp = mod(y+1, Ly)
      nnnc = nnnc + 1
      nnn_list(nnnc) = yp*Lx + xp
    end if
    if ( (x+1 < Lx .or. PBC_x) .and. (y-1 >= 0 .or. PBC_y) ) then
      xp = mod(x+1, Lx); yp = mod(y-1+Ly, Ly)
      nnnc = nnnc + 1
      nnn_list(nnnc) = yp*Lx + xp
    end if
    if ( (x-1 >= 0 .or. PBC_x) .and. (y+1 < Ly .or. PBC_y) ) then
      xp = mod(x-1+Lx, Lx); yp = mod(y+1, Ly)
      nnnc = nnnc + 1
      nnn_list(nnnc) = yp*Lx + xp
    end if
    if ( (x-1 >= 0 .or. PBC_x) .and. (y-1 >= 0 .or. PBC_y) ) then
      xp = mod(x-1+Lx, Lx); yp = mod(y-1+Ly, Ly)
      nnnc = nnnc + 1
      nnn_list(nnnc) = yp*Lx + xp
    end if
  end subroutine get_nn_nnn

  !----------------------------------------------------------------
  ! Try to hop one particle and update vec_out accordingly
  ! (used in building H * v)
  !----------------------------------------------------------------
  subroutine try_hop_sparse(active_bits, passive_bits, nbits, site, dest, amp, is_up_spin, vec_in, bckwindx, vec_out, k)
    implicit none
    integer, intent(in) :: active_bits, passive_bits, nbits, site, dest, k
    real*8, intent(in) :: amp
    logical, intent(in) :: is_up_spin
    real*8, intent(in) :: vec_in(:)
    real*8, intent(inout) :: vec_out(:)
    integer, intent(in) :: bckwindx(0:)

    integer :: tmp_bits, new_bits, combined_key, jidx
    integer :: sgn1, sgn2

    if (.not. btest(active_bits, site)) return
    if (btest(active_bits, dest)) return

    sgn1 = fermion_sign(active_bits, site)
    tmp_bits = ibclr(active_bits, site)
    sgn2 = fermion_sign(tmp_bits, dest)
    new_bits = ibset(tmp_bits, dest)

    if (is_up_spin) then
      combined_key = new_bits + ishft(passive_bits, nbits)
    else
      combined_key = passive_bits + ishft(new_bits, nbits)
    end if

    jidx = bckwindx(combined_key)
    if (jidx /= 0) then
      vec_out(jidx) = vec_out(jidx) - amp * dble(sgn1 * sgn2) * vec_in(k)
    end if
  end subroutine try_hop_sparse

  !----------------------------------------------------------------
  ! Hamiltonian action: vec_out = H * vec_in
  !----------------------------------------------------------------
  subroutine hamiltonian_action(vec_in, vec_out, nstates, state_up, state_dn, bckwindx)
    implicit none
    integer, intent(in) :: nstates
    integer, intent(in) :: state_up(nstates), state_dn(nstates)
    integer, intent(in) :: bckwindx(0:)
    real*8, intent(in) :: vec_in(nstates)
    real*8, intent(out) :: vec_out(nstates)

    integer :: k, up_bits, dn_bits, site, dest, i, nnc, nnnc
    integer :: nn_list(4), nnn_list(4)
    real*8 :: val
    real*8, allocatable :: tmp_out(:)

    allocate(tmp_out(nstates))
    tmp_out = 0.d0

    do k = 1, nstates
      up_bits = state_up(k)
      dn_bits = state_dn(k)
      val = 0.d0

      ! Diagonal: U and mu
      do site = 0, n - 1
        if (btest(up_bits, site) .and. btest(dn_bits, site)) val = val + U
        if (btest(up_bits, site)) val = val - mu
        if (btest(dn_bits, site)) val = val - mu
      end do
      tmp_out(k) = tmp_out(k) + val * vec_in(k)

      ! Hopping terms
      do site = 0, n - 1
        call get_nn_nnn(site, Lx, Ly, PBC_x, PBC_y, nn_list, nnc, nnn_list, nnnc)

        do i = 1, nnc
          dest = nn_list(i)
          call try_hop_sparse(up_bits, dn_bits, n, site, dest, t, .true., vec_in, bckwindx, tmp_out, k)
          call try_hop_sparse(dn_bits, up_bits, n, site, dest, t, .false., vec_in, bckwindx, tmp_out, k)
        end do

        do i = 1, nnnc
          dest = nnn_list(i)
          call try_hop_sparse(up_bits, dn_bits, n, site, dest, tprime, .true., vec_in, bckwindx, tmp_out, k)
          call try_hop_sparse(dn_bits, up_bits, n, site, dest, tprime, .false., vec_in, bckwindx, tmp_out, k)
        end do
      end do
    end do

    vec_out = tmp_out
    deallocate(tmp_out)
  end subroutine hamiltonian_action



subroutine lanczos_ground(state_up, state_dn, bckwindx, nstates, max_iter, tol, eigval, psi, numv)
  implicit none
  integer, intent(in) :: state_up(:), state_dn(:), bckwindx(0:)
  integer, intent(in) :: nstates, max_iter, numv
  real*8, intent(in) :: tol
  real*8, intent(out) :: eigval
  real*8, intent(out) :: psi(nstates)

  ! Lanczos internal arrays
  real*8, allocatable :: v0(:), v1(:), v2(:)
  real*8, allocatable :: alpha_list(:), beta_list(:)
  real*8, allocatable :: T(:,:), eigs(:), work(:)
  real*8, allocatable :: V(:,:)         ! store Lanczos basis vectors
  real*8 :: work_query(numv)
  integer :: iter, i, info, lwork
  real*8 :: alpha, beta, rnorm

  ! Allocate working storage
  allocate(v0(nstates), v1(nstates), v2(nstates))
  allocate(alpha_list(max_iter))
  if (max_iter > 1) then
     allocate(beta_list(max_iter-1))
  else
     allocate(beta_list(numv))  ! minimal allocation to avoid zero-size issues
  end if
  allocate(V(nstates, max_iter))

  ! Initialize
  call random_number(v0)
  v0 = v0 / sqrt(max(1.0d-300, dot_product(v0, v0)))  ! normalized initial vector
  v1 = 0.d0
  beta = 0.d0
  iter = 0

  ! Lanczos loop --- build orthonormal basis V(:,1:iter) and tridiagonal coefficients
  do
    iter = iter + 1
    ! store current Lanczos vector
    V(:, iter) = v0

    ! w = H * v0  (uses your provided routine)
    call hamiltonian_action(v0, v2, nstates, state_up, state_dn, bckwindx)

    if (iter > 1) then
      v2 = v2 - beta * v1
    end if

    alpha = dot_product(v0, v2)
    v2 = v2 - alpha * v0

    rnorm = sqrt(max(0.d0, dot_product(v2, v2)))
    alpha_list(iter) = alpha
    if (iter > 1) beta_list(iter-1) = beta

    ! Check convergence / termination
    if (rnorm < tol .or. iter >= max_iter) exit

    !Next step
    v1 = v0
    v0 = v2 / rnorm
    beta = rnorm
  end do

  ! Tridiagonal T of size iter x iter
  allocate(T(iter, iter))
  T = 0.d0
  do i = 1, iter
    T(i,i) = alpha_list(i)
  end do
  do i = 1, iter-1
    T(i,i+1) = beta_list(i)
    T(i+1,i) = beta_list(i)
  end do

  ! Diagonalize T using DSYEV (for small iter it's fine)
  allocate(eigs(iter))
  ! Work space query
  call dsyev('V', 'U', iter, T, iter, eigs, work_query, -1, info)
  if (info /= 0) then
    print *, 'DSYEV workspace query failed, info=', info
    stop
  end if
  lwork = max(1, int(work_query(1)))
  allocate(work(lwork))

  ! Actual diagonalization: T will be overwritten by eigenvectors of T; eigs holds eigenvalues (ascending).
  call dsyev('V', 'U', iter, T, iter, eigs, work, lwork, info)
  if (info /= 0) then
    print *, 'DSYEV failed, info=', info
    stop
  end if

  ! smallest eigenvalue (DSYEV returns eigenvalues in ascending order)
  eigval = eigs(1)

  ! Reconstruct ground-state vector psi = sum_j ( eigenvector_T(j,1) * V(:,j) )
  psi = 0.d0
  do i = 1, iter
    psi = psi + T(i,1) * V(:, i)
  end do

  ! normalize psi
  rnorm = sqrt(max(1.0d-300, dot_product(psi, psi)))
  psi = psi / rnorm

  ! Clean 
  deallocate(v0, v1, v2, alpha_list, beta_list, V, T, eigs, work)

end subroutine lanczos_ground


  !----------------------------------------------------------------
  ! Fermionic sign
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
       !if(abs(r)>0) write(11,'(I4,2X,I4,2X,F6.2,2X,F12.8)') i, j, r, C(i,j)
       if(r>0) write(11,'(I4,2X,I4,2X,F6.2,2X,F12.8)') i, j, r, C(i,j)
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
  open(unit=121,file='Total_den_docc_2d_Lanc.dat',status='replace')
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

end program hubbard_ed_2d
