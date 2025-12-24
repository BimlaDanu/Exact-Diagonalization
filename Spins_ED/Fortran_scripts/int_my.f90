module numbr_spn
integer,parameter::n=6,nb=n-1,nup=n/2,J_ex=1
end module numbr_spn
!-------------------------------!
program main
use numbr_spn
implicit none
integer, allocatable :: state(:)
integer, allocatable :: bckwindx(:)
real*8, allocatable :: H(:,:)
real*8, allocatable :: eig(:)
!--------------------------!
integer*4::i,pos,count,j,fact,nstates,k,bond,di,da,ix,k1,bond1,site1,site2,m
integer, external :: fncup
!print*,fncup(n,nup)
nstates=fncup(n,nup)
!nstates=12870
allocate (state(nstates))
allocate (bckwindx(2**n))
bckwindx = 0
nstates=0
do i =0,2**n-1
count = 0
do pos=0,n-1
if(btest(i,pos).eqv..true.) count = count + 1	
end do
!if (count.eq.nup)then
nstates=nstates+1
state(nstates)=i
bckwindx(i)=nstates
!print*,state(nstates),bckwindx(i)
!end if 
end do
!---------------------------!
allocate(H(nstates,nstates))
H(:,:)=0.d0
do k=1,nstates
do bond =1,nb
site1=2**(bond-1)
bond1=bond+1
!print*,'bond, bond1, nstates= ',bond, bond1,nstates
if(n.gt.2.and.bond1>n)bond1=bond1-n 
!if(bond1>n)bond1=bond1-n 
!print*,'bond, bond1, nstates= ',bond, bond1,nstates
site2=2**(bond1-1)
di=site1+site2
da=iand(state(k),di)
!print*,site1,site2,di,da, nstates
if(da.eq.di.or.da.eq.0)then
H(k,k)=H(k,k)+0.25d0*J_ex
else
H(k,k)=H(k,k)-0.25d0*J_ex
ix=ieor(state(k),di)
j=bckwindx(ix)
H(j,k)=H(j,k)+0.5d0*J_ex
end if
end do
enddo
!do i=1,nstates
!write(*,*)(H(i,j),j=1,nstates)
!end do
!---------------------------!
allocate(eig(nstates))
call sym(H,eig,nstates)
print*,(i,eig(i),i=1,nstates)
do i=1,nstates
!write(*,*)(H(i,j),j=1,nstates)
end do
end
!--------------------------!


function fact(n)
integer*4::fact,j,n
fact=1
do j=2,n
fact=fact*j
end do
end function
!----------------------!

function fncup(n,nup)
integer*4::fncup,nup, n
integer, external :: fact
fncup=nint(dble(fact(n))/dble(fact(nup)*fact(n-nup)))
end function
!-----------------------------------!

subroutine sym(H,eig,n)
implicit none
integer n,lwork,inf
real*8  H(n,n),eig(n),work(max(1,3*n-1))
lwork=3*n-1
call dsyev('V','U',n,H,n,eig,work,lwork,inf)
end subroutine sym




