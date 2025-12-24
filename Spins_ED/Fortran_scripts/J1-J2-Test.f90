module numbr_spn
integer,parameter::n=9,nb=n-1,nup=(n+1)/2
end module numbr_spn

module real_spn
real*8,parameter::J_ex=1.d0,J_exn=0.d0, D=0.0
end module real_spn
!---------------------------------------------------------------------------------!


program main
use numbr_spn
use real_spn
implicit none
integer, allocatable :: state(:)
integer, allocatable :: bckwindx(:)
real*8, allocatable :: H(:,:)
real*8, allocatable :: eig(:)
real*8, allocatable :: x(:)
real*8, allocatable :: x1(:)
real*8, allocatable :: x2(:)
!--------------------------!
real*8 szz(n*n),sxx(n*n),mz(n)
integer*4::i,pos,count,j,fact,nstates,k,bond,di,da,ix,k1,bond1,bond2
integer*4:: site1,site2,site3,m, da1, di1, ix1,j1, indx, ind
integer, external :: fncup


open(1,file='Szz-Corr_L9.dat',status='unknown')
open(2,file='Sxx-Corr_L9.dat',status='unknown')


indx=1


!print*,fncup(n,nup)
!nstates=fncup(n,nup)
!nstates=12870


nstates=2**n

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
!---------------------------------------------------------------------------------!


allocate(H(nstates,nstates))
H(:,:)=0.d0
do k=1,nstates
do bond =1,nb
site1=2**(bond-1)
bond1=bond+1
bond2=bond+2
!print*,'bond, bond1, nstates= ',bond, bond1,bond2,nstates
if(n.gt.2.and.bond1>n)bond1=bond1-n 
if(n.gt.2.and.bond2>n)bond2=bond2-n 
!if(bond1>n)bond1=bond1-n 
!print*,'bond, bond1, nstates= ',bond, bond1,bond2,nstates

!!
site2=2**(bond1-1)
di=site1+site2
da=iand(state(k),di)

di1=site1+site1
da1=iand(state(k),di1)

if(da.eq.di.or.da.eq.0)then
H(k,k)=H(k,k)+(0.25d0*D)!*(-3./4.)
else
H(k,k)=H(k,k)-(0.25d0*D)!*(-3./4.)
endif

if(da.eq.di.or.da.eq.0)then
H(k,k)=H(k,k)+(0.25d0*J_ex)*(1.0/1.0)
else
H(k,k)=H(k,k)-(0.25d0*J_ex)*(1.0/1.0)
ix=ieor(state(k),di)
j=bckwindx(ix)
H(j,k)=H(j,k)+0.5d0*J_ex*1
end if
!!
!!
site3=2**(bond2-1)
di1=site1+site3
da1=iand(state(k),di1)
if(da1.eq.di1.or.da1.eq.0)then
H(k,k)=H(k,k)+0.25d0*J_exn/1.
else
H(k,k)=H(k,k)-0.25d0*J_exn/1.
ix1=ieor(state(k),di1)
j1=bckwindx(ix1)
H(j1,k)=H(j1,k)+0.5d0*J_exn
end if
!!
!!
end do
enddo
!do i=1,nstates
!write(*,*)(H(i,j),j=1,nstates)
!end do
!---------------------------------------------------------------------------------!

allocate(eig(nstates))
allocate (x(nstates))
allocate (x1(nstates))
allocate (x2(nstates))

x(nstates) = 0.d0
x1(nstates) = 0.d0
x2(nstates) = 0.d0


call sym(H,eig,nstates)
print*,(i,eig(i),i=1,6)
!print*,(i,eig(i),i=1,nstates)
do i=1,nstates
!write(*,*)(H(i,j),j=1,nstates)
x(i) = H(i,1)
x1(i) = H(i,1)
x2(i) = H(i,2)
end do


call xcorr(nstates,state,bckwindx,sxx,x,x1,x2,indx)
call  zcorr(nstates,state,szz,x,x1,x2,indx)
!call m_avg(nstates,state,mz,x,x1, x2,indx)
end
!--------------------------------------------------------------------------------!



function fact(n)
integer*4::fact,j,n
fact=1
do j=2,n
fact=fact*j
end do
end function
!--------------------------------------------------------------------------------!

function fncup(n,nup)
integer*4::fncup,nup, n
integer, external :: fact
fncup=nint(dble(fact(n))/dble(fact(nup)*fact(n-nup)))
end function
!--------------------------------------------------------------------------------!

subroutine sym(H,eig,n)
implicit none
integer n,lwork,inf
real*8  H(n,n),eig(n),work(max(1,3*n-1))
lwork=3*n-1
call dsyev('V','U',n,H,n,eig,work,lwork,inf)
end subroutine sym
!--------------------------------------------------------------------------------!


subroutine xcorr(nstates,state,bckwindx,sxx,x,x1,x2,indx)
use numbr_spn
use real_spn
implicit none
real*8 sxx(n*n),x(nstates),x1(nstates),x2(nstates),corr
integer*4::k,i,j,l,di,da,bond,ix,nstates,site1,site2,s1,s2,s
integer*4::state(nstates),a,b,r,bckwindx(2**n),indx

do 10 bond =1,n
i =bond
do j =i+1,n
r =abs(i-j)
s1=2**(i-1)
s2=2**(j-1)

corr=0.d0
s=s1+s2
do 20 k=1,nstates
da=iand(state(k),s)
if(da.eq.0.or.da.eq.s) go to 20
ix=ieor(state(k),s)
l=bckwindx(ix)
!print*,'indx=',indx
if(indx.eq.0)corr=corr+x(k)*x(l)
if(indx.eq.1)corr=corr+x(k)*x(l)
if(indx.eq.2)corr=corr+(x1(k)*x1(l)+x2(k)*x2(l))/2.d0
20 continue
sxx(bond)=corr/4.d0
!print*, 'bond, r, sxx(bond)=',bond, r, sxx(bond)
print*, 'i, j, r, sxx(bond)=',i, j, r, sxx(bond), 4.0*sxx(bond), sxx(bond)/4.0
write(2,26)i,j, r, sxx(bond), 4.0*sxx(bond), sxx(bond)/4.0, 2.0*sxx(bond)
enddo
10	continue
return
26  format(1X,I4, 1X,I4, 1X,I4, 1X,24f12.6)!
end subroutine xcorr
!--------------------------------------------------------------------------------!

subroutine zcorr(nstates,state,szz,x,x1, x2,indx)
use numbr_spn
use real_spn
implicit none
real*8 szz(n*n),x(nstates),x1(nstates),x2(nstates),corr
integer*4:: k,i,j,di,da,bond,ix,nstates,site1,site2,s1,s2,s,factor
integer*4:: state(nstates),r,indx


do 10 bond =1,n
i =bond
do j =i+1,n
r =abs(i-j)
s1=2**(i-1)
s2=2**(j-1)



corr=0.d0
s=s1+s2
do 20 k=1,nstates
!print*,'x(k)=',x(k)
da=iand(state(k),s)
if(da.eq.0.or.da.eq.s)then
factor=1
else
factor=-1
endif
if(indx.eq.0)corr=corr+factor*x(k)**2
if(indx.eq.1)corr=corr+factor*x(k)**2
if(indx.eq.2)corr=corr+factor*(x1(k)**2+x2(k)**2)/2.d0
20 continue

szz(bond)=corr/4.d0
!print*, 'bond, r, szz(bond)=',bond, r, szz(bond)
print*, 'i, j, r, szz(bond)=',i,j, r, szz(bond), 4.0*szz(bond), szz(bond)/4.0
write(1,26)i,j, r, szz(bond), 4.0*szz(bond), szz(bond)*3.d0, szz(bond)*2.d0
enddo
10	continue
return
26  format(1X,I4, 1X,I4, 1X,I4, 1X,24f12.6)!	
end subroutine zcorr




subroutine m_avg1(nstates,state,mz,x,x1, x2,indx)
use numbr_spn
use real_spn
implicit none
real*8 mz(n),x(nstates),x1(nstates),x2(nstates),corr
integer*4:: k,i,j,di,da,bond,ix,nstates,site1,s1,s,factor
integer*4:: state(nstates),r,indx

mz(:) =0.d0
do i  =0,n-1
site1 = 2**i
j =i+1
corr=0.d0
do  k=1,nstates
!if(btest(site1,state(k)).eqv..true.)then
if(btest(i,state(k)).eqv..true.)then
factor=1.d0
else
factor=-1.d0
endif
if(indx.eq.0)corr=corr+factor*x(k)**2/2.d0
!corr=corr+factor*(x1(k)**2+x2(k)**2)/2.d0
end do

mz(j)=corr/2.d0
!print*, 'i, j, mz(j)=',i, j, mz(j)
!print*, 'i, j,  mz(j)=',i,j, mz(j), 4.0*mz(j), mz(j)/4.0
write(1,26)i,j, mz(j), 4.0*mz(j), mz(j)/4.0
end do
return

26  format(1X,I4, 1X,I4, 1X,24f12.6)!	
end subroutine m_avg1






subroutine m_avg(nstates,state,mz,x,x1, x2,indx)
use numbr_spn
use real_spn
implicit none
real*8 mz(n),x(nstates),x1(nstates),x2(nstates),corr
integer*4:: k,i,j,di,da,bond,ix,nstates,site1,s1,s,factor
integer*4:: state(nstates),r,indx

mz(:) =0.d0

do 10 bond =1,n
i =bond
r =i
s1=2**(i-1)

corr=0.d0
s=s1
do 20 k=1,nstates
da=iand(state(k),s)
if(da.eq.s)then
factor=1.d0
else
factor=-1.d0
endif
if(indx.eq.0)corr=corr+factor*x(k)**2
!corr=corr+factor*(x1(k)**2+x2(k)**2)/2.d0
20 continue

mz(bond)=corr/2.d0
!print*, 'bond, r, mz(bond)=',bond, r, mz(bond)
print*, 'i, j, r, mz(bond)=',i,j, r, mz(bond), 4.0*mz(bond), mz(bond)/4.0
write(1,26)i,j, r, mz(bond), 4.0*mz(bond), mz(bond)/4.0
10	continue
return
26  format(1X,I4, 1X,I4, 1X,I4, 1X,24f12.6)!	
end subroutine m_avg


