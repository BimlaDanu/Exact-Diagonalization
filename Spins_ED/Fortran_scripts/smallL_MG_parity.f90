module numbr_spn
integer,parameter::n=4,nup=2,nb=n,nb2=n,szval=n/2-nup,hfbit=2**((n+1)/2)	&
&	,rght=2**((n+1)/2)-1,lft=ieor(2**n-1,rght)
end module numbr_spn
!---------------------------------!
program main
use numbr_spn
implicit none
integer, allocatable :: state(:)
integer, allocatable :: bckwindx(:)
real*8, allocatable :: H(:,:)
real*8, allocatable :: H1(:,:)
real*8, allocatable :: H2(:,:)
real*8, allocatable :: H4(:,:)
real*8, allocatable :: eig(:)
real*8, allocatable :: wk(:,:)
real*8, allocatable :: iwk(:)
real*8, allocatable :: v(:)
real*8, allocatable :: x(:)
real*8, allocatable :: x1(:)
real*8, allocatable :: x2(:)
integer*4::i,j,nstates,np,cont,ja,jb,bjp,ideclr,ne,nvec,k
integer*4::sz,state2(2,0:2**n-1),a,b,jc,npair(2*12),npair1(2*13),npair2(2*13)
real*8:: eps,E(4),hexpec,sxx(16),szz(16),dimcorr(13)
integer, external :: fncup
!-----------------------------------!
eps=1.d-13
nvec=1
ne=4
nstates =fncup(n,nup)
ideclr=nstates
allocate (state(nstates))
allocate (H(nstates,nstates))
allocate (H1(nstates,nstates))
allocate (H2(nstates,nstates))
allocate (H4(nstates,nstates))
allocate (wk(ideclr,8))
allocate (iwk(ideclr))
allocate (v(ideclr))
allocate (x(ideclr))
allocate (x1(ideclr))
allocate (x2(ideclr))
np=n/2+mod(n,2)+(szval+0.001d0)
cont=0
ja=0
jb=0
bjp=0
do 10 i=1,2**n
sz=0
do 20 j=0,n-1
20 sz=sz+mod(i/2**j,2)
if(sz.ne.np)go to 10
cont=cont+1
if (cont.gt.nstates)then
print *,' Incorrect nstates given to sz'
stop
end if
a=iand(i,rght)
b=iand(i,lft)/hfbit
if (b.eq.bjp)then
ja=ja+1
else
bjp=b
ja=1
jb=cont-1
end if
state(cont)=i
state2(1,a)=ja
state2(2,b)=jb
jc=ja+jb
!print*,a,ja,b,jb,jc,state(cont)
10 continue
!if(cont.eq.nstates)return
!stop
call hamiltonian(H1,nstates,state,state2)
call hamiltonian1(H2,nstates,state,state2)
H=H1+H2
!H=-1.d0*H1+1.5d0*H2!0.5*0*H1+H2
!H4=2.d0*H!+N*0.1.5
allocate(eig(nstates))
call sym(H,eig,nstates)

!call spinsquared(H4,nstates,state,state2)
call transform(nstates,H4,H,x)
print*,'x(2)=',x(2)
x(2)=0.5d0*abs(sqrt(1.d0+4.d0*(x(2)+0.75d0*12))-1.d0)
print*,'x(2)=',x(2)
print*,(i,eig(i)/real(n),i=1,5)
!!------------------------------!
do k=1,nstates
x1(k)=H(k,1)
enddo
!data npair /1,2, 3,4, 5,6, 7,8, 9,10, 11,12, 1,3, 3,5, 5,7, 7,9, 9,11	!&
!&	,2,4/!, 4,6, 6,8, 8,10, 10,12, 11,2, 12,1/

#data npair /1,2, 2,3, 3,4, 4,5, 5,6, 6,7, 7,8, 8,9, 9,10, 10,11, 11,12, 12,1/	!&
!data npair /1,3, 2,4, 3,5, 4,6, 5,7, 6,8, 7,9, 8,10, 9,11, 10,12, 11,1, 12,2/

data npair /1,2, 2,3, 3,4, 4,1/	

!npair(1)=1
!npair(2)=2
call xcorr(12,npair,nstates,state,state2,sxx,x1)
call zcorr(12,npair,nstates,state,state2,szz,x1)
!
!data npair1 /1,2, 1,2, 1,2, 1,2, 1,2, 1,2, 1,2, 1,2, 1,2, 1,2, 1,2, 1,2, 1,2/
!data npair2 /3,4, 5,6,  7,8, 9,10, 11,12, 3,5, 5,7, 7,9, 9,11	&
!&	,4,6, 6,8, 8,10, 10,12/!
!call dim2corr(13,npair1,npair2,nstates,state,state2,dimcorr,x1)
stop
end
!---------------------------!

subroutine hamiltonian(H,nstates,state,state2)
use numbr_spn
implicit none
real*8::H(nstates,nstates),diag,wght,bondwt(nb),zrtio(nb)
integer*4::k,j,di,da,bond,ix,nstates,site1,site2,s1,s2,s,cont
integer*4::state(nstates),a,b
integer*4::state2(2,0:2**n-1),ipair(2*nb),j1,new
data bondwt/nb*-1/
data zrtio/nb*1/

!data ipair /1,3, 3,5, 5,7, 7,9, 9,11, 2,4, 4,6, 6,8, 8,10, 10,12/	
!data ipair /1,2, 2,3, 3,4, 4,5, 5,6, 6,7, 7,8, 8,9, 9,10, 10,11, 11,12/!, 12,1/	

data ipair /1,2, 2,3, 3,4, 4,1/


!---------------------------------!
H(:,:)=0
!---------------------------------!
do  30 bond =1,nb
site1=ipair(2*bond-1)-1
site2=ipair(2*bond )-1
s1=2**site1
s2=2**site2
s=s1+s2
wght=bondwt(bond)
diag=0.5*zrtio(bond)
!print*,bondwt(bond),zrtio(bond),site1,site2,wght,diag,-0.5d0*wght,0.5d0*diag
!---------------------------------!
do 30 k=1,nstates
da=iand(state(k),s)
if(da.eq.0.or.da.eq.s)then
H(k,k)=H(k,k)+0.5d0*diag/4.d0
else
H(k,k)=H(k,k)-0.5d0*diag/4.d0
ix=ieor(state(k),s)
a=iand(ix,rght)
b=iand(ix,lft)/hfbit
j=state2(1,a)+state2(2,b)
H(k,j)=-0.5d0*wght
end if
30 continue
return
end subroutine hamiltonian
!---------------------------!


!---------------------------!
subroutine hamiltonian1(H,nstates,state,state2)
use numbr_spn
implicit none
real*8::H(nstates,nstates),diag,wght,bondwt(nb2),zrtio(nb2)
integer*4::k,j,di,da,bond,ix,nstates,site1,site2,s1,s2,s,cont
integer*4::state(nstates),a,b
integer*4::state2(2,0:2**n-1),ipair(2*nb2),j1,new
data bondwt/nb2*-1/
data zrtio/nb2*1/
!data ipair/1,2, 3,4 ,5,6, 7,8, 9,10, 11,12/
!data ipair/1,7, 2,8 ,3,9, 4,10, 5,11, 6,12/
!data ipair /1,3, 2,4, 3,5, 4,6, 5,7, 6,8, 7,9, 8,10, 9,11, 10,12/!, 11,1, 12,2/	

data ipair /1,3, 2,4, 3,1, 4,2/

H(:,:)=0
do  30 bond =1,nb2
site1=ipair(2*bond-1)-1
site2=ipair(2*bond )-1
s1=2**site1
s2=2**site2
s=s1+s2
wght=bondwt(bond)
diag=0.5*zrtio(bond)
!print*,bondwt(bond),zrtio(bond),site1,site2,wght,diag,-0.5d0*wght,0.5d0*diag
!---------------------------------!
do 30 k=1,nstates
da=iand(state(k),s)
if(da.eq.0.or.da.eq.s)then
H(k,k)=H(k,k)+0.5d0*diag/4.d0
else
H(k,k)=H(k,k)-0.5d0*diag/4.d0
ix=ieor(state(k),s)
a=iand(ix,rght)
b=iand(ix,lft)/hfbit
j=state2(1,a)+state2(2,b)
H(k,j)=-0.5d0*wght
end if
30 continue
return
end subroutine hamiltonian1
!---------------------------!


!-------------------------!
function fact(n)
integer*4::fact,j,n
fact=1
do j=2,n
fact=fact*j
end do
end function
!----------------------!

!-----------------------!
function fncup(n,nup)
integer*4::fncup,nup, n
integer, external :: fact
fncup=nint(dble(fact(n))/dble(fact(nup)*fact(n-nup)))
end function
!------------------------!

!------------------------!
subroutine szt(nstates,state,state2)
use numbr_spn
implicit none
integer*4::i,j,np,cont,ja,jb,bjp,a,b,jc,sz
integer*4::nstates,state(nstates),state2(2,0:2**n-1)
cont=0
ja=0
jb=0
bjp=0
do 10 i=1,2**n
sz=0
do 20 j=0,n-1
20 sz=sz+mod(i/2**j,2)
if(sz.ne.nup)go to 10
cont=cont+1
if (cont.gt.nstates)then
print *,' Incorrect nstates given to sz'
stop
end if
a=iand(i,rght)
b=iand(i,lft)/hfbit
if (b.eq.bjp)then
ja=ja+1
else
bjp=b
ja=1
jb=cont-1
end if
state(cont)=i
state2(1,a)=ja
state2(2,b)=jb
jc=ja+jb
print*,a,ja,b,jb,jc,state(cont)
10 continue
!if(cont.eq.nstates)return
!stop
end subroutine szt
!------------------!

subroutine sym(H,eig,n)
implicit none
integer n,lwork,inf
real*8  H(n,n),eig(n),work(max(1,3*n-1))
lwork=3*n-1
call dsyev('V','U',n,H,n,eig,work,lwork,inf)
end subroutine sym
!----------------------!

!-----------------------!
! eigenvalues of a small matrix 
subroutine diag(H,ideclr,nstates,E,v,ne,nvec,eps,wk,iwk)
use numbr_spn
implicit none
real*8:: E(ne),v(ideclr,nvec),H(ideclr,nstates),eps
real*8:: wk(ideclr,8),iwk(ideclr)
integer*4::ne,nvec,nstates,ideclr

if(nvec.lt.0.or.nvec.gt.ne)then
print *,'nvec given to diag out of range'
stop
end if

call hshldr(H,ideclr,nstates,wk(1,1),wk(1,2),wk(1,3),wk(1,4),wk(1,5),wk(1,6))
call bisec(wk(1,1),wk(1,2),nstates,E,ne,eps)
if(nvec.eq.0)return
call vec3(E,H,ideclr,nstates,ne,nvec,wk(1,4),wk(1,5),wk(1,6)	&
&	,wk(1,7),wk(1,8),iwk,wk(1,1),wk(1,2),wk(1,3),v)
end subroutine diag
!---------------------------!

subroutine hshldr(H,ideclr,nstates,alpha,beta,c,w,p,q)
use numbr_spn
implicit none
real*8:: H(ideclr,nstates),alpha(nstates),beta(nstates),c(nstates)
real*8:: w(nstates),p(nstates),q(nstates),s,t
integer*4::i,j,k,nstates,ideclr

do 10 k=1,nstates-2
s=0.d0
do 20 i=k+1,nstates
20	s=s+H(i,k)**2
s=sqrt(s)
if(H(k+1,k).lt.0.d0)s=-s
alpha(k)=H(k,k)
beta(k)=-s
c(k)=0.0d0
if(s**2.lt.1.d-26)goto 10
c(k)=1.d0/(s**2+H(k+1,k)*s)
w(k+1)=H(k+1,k)+s
do 30 i=k+2,nstates
30	w(i)=H(i,k)
H(k+1,k)=w(k+1)
do 40 i=k+1,nstates
t=0.d0
do 50 j=k+1,i
50  t=t+H(i,j)*w(j)
do 55 j=i+1,nstates
55	t=t+H(j,i)*w(j)
p(i)=c(k)*t
40	continue
t=0.d0
do 60 i=k+1,nstates
60	t=t+p(i)*w(i)
s=0.5d0*c(k)*t
do 70 i=k+1,nstates
70	q(i)=p(i)-s*w(i)
do 80 j=k+1,nstates
do 80 i=j,nstates
80	H(i,j)=H(i,j)-w(i)*q(j)-q(i)*w(j)
10   continue
alpha(nstates-1)=H(nstates-1,nstates-1)
alpha(nstates)=H(nstates,nstates)
beta(nstates-1)=H(nstates,nstates-1)
return
end subroutine hshldr
!--------------------------!


!eigenvector of a tridiagonal matrix by inverse iteration ****
subroutine vec3(E,H,ideclr,nstates,ne,nvec,di,bl,bu,bv,cm,lex,alpha,beta,c,v)
use numbr_spn
implicit none
real*8:: E(ne),H(ideclr,nstates)
real*8::di(nstates),bl(nstates),bu(nstates),bv(nstates),cm(nstates),lex(nstates)
real*8:: alpha(nstates),beta(nstates),c(nstates),v(ideclr,nvec),prd,dnorm,s
integer*4::i,j,k,nstates,nvec,km,ideclr,ne,l

if(nvec.gt.ne)then
print *,'nvec given to vec3 out of range'
stop
end if
do 10 k=1,nvec
do 100 j=1,nstates
di(j)=E(k)-alpha(j)
bl(j)=-beta(j)
bu(j)=-beta(j)
100    continue
!LU decomposition
do 110 j=1,nstates-1
if(abs(di(j)).gt.abs(bl(j)))then
!non pivoting
lex(j)=0
if(abs(di(j)).lt.1.d-13)di(j)=1.d-13
cm(j+1)=bl(j)/di(j)
di(j+1)=di(j+1)-cm(j+1)*bu(j)
bv(j)=0.d0
else
!pivoting
lex(j)=1
cm(j+1)=di(j)/bl(j)
di(j)=bl(j)
s=bu(j)
bu(j)=di(j+1)
bv(j)=bu(j+1)
di(j+1)=s-cm(j+1)*bu(j)
bu(j+1)= -cm(j+1)*bv(j)
end if
110    continue
if(abs(di(nstates)).lt.1.d-13)di(nstates)=1.d-13
! initial vector
do 120 j=1,nstates
120    v(j,k)=1.d0/(float(j)*5.d0)

! degeneracy check up
if(k.eq.1)then
km=k
else if(abs(E(k)-E(km)).gt.1.d-13)then
km=k
else
do 130 i=km,k-1
prd=0.d0
do 140 j=1,nstates
140 prd=prd+v(j,i)*v(j,k)
do 150 j=1,nstates
150	v(j,k)=v(j,k)-prd*v(j,i)
130	continue
end if
!inverse iteration
do 160 l=1,k-km+3
if((l.ne.1).or.(k.ne.km))then
! forward substitution
do 170 j=1,nstates-1
if(lex(j).eq.0)then
v(j+1,k)=v(j+1,k)-cm(j+1)*v(j,k)
else
s=v(j,k)
v(j,k)=v(j+1,k)
v(j+1,k)=s-cm(j+1)*v(j,k)
end if
170	continue
end if
! backward substitution
do 180 j=nstates,1,-1
s=v(j,k)
if(j.le.nstates-1)s=s-bu(j)*v(j+1,k)
if(j.le.nstates-2)s=s-bv(j)*v(j+2,k)
v(j,k)=s/di(j)
180	continue
!normalization
dnorm=0.d0
do 190 j=1,nstates
190	dnorm=dnorm+v(j,k)**2
if(dnorm.gt.1.d-13)dnorm=1./sqrt(dnorm)
do 200 j=1,nstates
200	v(j,k)=v(j,k)*dnorm
160	continue
10	continue

!back transformation to the original representation
do 210 k=1,nvec
do 220 i=nstates-2,1,-1
prd=0.d0
do 230 j=i+1,nstates
230	prd=prd+H(j,i)*v(j,k)
s=prd*c(i)
do 240 j=i+1,nstates
240      v(j,k)=v(j,k)-s*H(j,i)
220    continue
210  continue
!orthogonalization for degenerate case
km=1
do 250 k=2,nvec
if(abs(E(k)-E(km)).ge.1.0d-13)then
km=k
else
do 260 i=km,k-1
prd=0.d0
do 270 j=1,nstates
270	prd=prd+v(j,i)*v(j,k)
do 280 j=1,nstates
280	v(j,k)=v(j,k)-prd*v(j,i)
260	continue
dnorm=0.0d0
do 290 j=1,nstates
290	dnorm=dnorm+v(j,k)**2
s=1.d0/sqrt(dnorm)
do 300 j=1,nstates
300	v(j,k)=v(j,k)*s
end if
250  continue
end subroutine vec3
!-----------------------------------!

!-------------------------------!
subroutine bisec(alpha,beta,nstates,E,ne,eps)
use numbr_spn
implicit none
real*8::alpha(nstates),beta(nstates),E(ne),b2(2000),eps,epsabs,range
real*8::a,b,c,g
integer*4::i,j,k,ne,nstates,numneg,ipass
if(nstates.gt.2000)then
print *,' nstates given to bisec exceeds 2000'
stop
end if
if(ne.gt.nstates.or.ne.le.0)then
print *,' ne given to bisec out of range'
stop
end if
!intial bound @@@@@@@@@@@@@@@@@@!

range=abs(alpha(1))+abs(beta(1))
do 10 k=2,nstates-1
10 range=max(range,abs(beta(k-1))+abs(alpha(k))+abs(beta(k)))
range=max(range,abs(beta(nstates-1))+abs(alpha(nstates)))
range=-range

b2(1)=0.d0
do 20 i=2,nstates
20 b2(i)=beta(i-1)**2

epsabs=abs(range)*eps
do 30 i=1,ne
30 E(i)=-range
b=range

!bisection method @@@@@@@@@@@@@@@@@@@@@@@@!
do 100 k=1,ne
a=E(k)
do 110 j=1,110
c=(a+b)/2.d0
if(abs(a-b).lt.epsabs)goto 100
numneg=0
g=1.d0
ipass=0

do 120 i=1,nstates
if(ipass.eq.0)then
g=c-alpha(i)-b2(i)/g
else if(ipass.eq.1)then
ipass=2
else
g=c-alpha(i)
ipass=0
end if
if(ipass.eq.0)then
if(g.le.0.d0)numneg=numneg+1
if(abs(g).le.abs(b2(i)*epsabs*eps))ipass=1
end if
120	continue

numneg=nstates-numneg
if(numneg.lt.k)then
b=c
else
a=c
do 130 i=k,min(numneg,ne)
130	E(i)=c
end if

110	continue
100	continue
end subroutine bisec
!----------------------------------!

!check of the eigenvector and eigenvalue 
subroutine check3(H,nstates,ideclr,x,v,Hexpec)
use numbr_spn
implicit none
real*8::H(ideclr,nstates),hexpec,dnorm,x(nstates),v(nstates),prd
integer*4::i,j,k,nstates,ideclr
!------------------------------!
dnorm=0.d0
do 5 j=1,nstates
5	dnorm=dnorm+x(j)**2
if(dnorm.lt.1.d-30)then
print *,'Null vector given to check3'
return
end if
do 10 i=1,nstates
10   v(i)=0.0d0
do 20 j=1,nstates
do 20 k=1,nstates
v(j)=v(j)+H(j,k)*x(k)
20   continue
prd=0.0d0
do 30 i=1,nstates
30   prd=prd+v(i)*x(i)
Hexpec=prd
print *
print 200
200  format('Information from check3')
print 210,prd
210  format(' <x*H*x> =',1pd16.8)
print 220
220  format(' H*x(j)/x(j) (j=min(nstates/3,13),nstates,max(1,nstates/20))')
print 230,(v(i)/x(i),i=min(nstates/3,13),nstates,max(1,nstates/20))
230  format(4d18.9)
print 240
240  format('-----------')
return
end subroutine check3
!-----------------------------------!

subroutine xcorr(nb1,npair,nstates,state,state2,sxx,x1)
use numbr_spn
implicit none
real*8::sxx(nb1),x1(nstates),x2(nstates),corr
integer*4::k,j,di,da,bond,ix,nstates,site1,site2,s1,s2,s
integer*4::state(nstates),a,b,nb1
integer*4::state2(2,0:2**n-1),npair(2*nb1)
!----------------------------------------!
open(1,file='spncorx.dat',status='unknown')
do 10 bond =1,nb1
!------------------!
site1=npair(2*bond-1)-1
site2=npair(2*bond )-1
s1=2**site1
s2=2**site2
corr=0.d0
s=s1+s2
do 20 k=1,nstates
da=iand(state(k),s)
if(da.eq.0.or.da.eq.s) go to 20
ix=ieor(state(k),s)
a=iand(ix,rght)
b=iand(ix,lft)/hfbit
corr=corr+x1(k)*x1(state2(1,a)+state2(2,b))
20 continue
sxx(bond)=corr/4.d0
print*,site1+1,site2+1,sxx(bond)
write(1,*) site1+1,site2+1,sxx(bond)
10	continue
return
end subroutine xcorr
!------------------------------------------------!

subroutine zcorr(nb1,npair,nstates,state,state2,szz,x1)
use numbr_spn
implicit none
real*8::szz(nb1),x1(nstates),x2(nstates),corr
integer*4::k,j,di,da,bond,ix,nstates,site1,site2,s1,s2,s,factor
integer*4::state2(2,0:2**n-1),npair(2*nb1),state(nstates),a,b,nb1
!---------------------------!
open(1,file='spncorz.dat',status='unknown')
do 10 bond =1,nb1
!------------------!
site1=npair(2*bond-1)-1
site2=npair(2*bond )-1
s1=2**site1
s2=2**site2
corr=0.d0
s=s1+s2
do 20 k=1,nstates
da=iand(state(k),s)
if(da.eq.0.or.da.eq.s)then
factor=1.d0
else
factor=-1.d0
endif
corr=corr+factor*x1(k)**2
20 continue
szz(bond)=corr/4.d0
print*,site1+1,site2+1,szz(bond)
write(1,*) site1+1,site2+1,szz(bond)
10	continue
return
end subroutine zcorr
!----------------------!

subroutine dim2corr(nb1,npair1,npair2,nstates,state,state2,dimcorr,x1)
use numbr_spn
implicit none
real*8::sxy(nb1),x1(nstates),x2(nstates),corr,dimcorr(nb1)
real*8::sxy1(nb1),sxy2(nb1),sxy3(nb1),sxyz(nb1),sxyz1(nb1)
real*8::szz(nb1),szz1(nb1),szz2(nb1),szz3(nb1),szz4(nb1),szz5(nb),sxyz2(nb)
integer*4::k,j,di,da,bond,ix,nstates,site1,site2,s1,s2,s,da1
integer*4::state(nstates),a,b,nb1,is1,is2,is,isite1,isite2,factor
integer*4::state2(2,0:2**n-1),npair1(2*nb1),npair2(2*nb1),ix1,ia,ib
!--------------------!
open(1,file='dim_corln.dat',status='unknown')
do 10 bond =1,nb1
!<sx_i sx_j+sy_i sy_j> nearest neighbor spin-2 correlation function
site1=npair1(2*bond-1)-1
site2=npair1(2*bond )-1
s1=2**site1
s2=2**site2
s=s1+s2
!--------------------------!
corr=0.d0
do 20 k=1,nstates
da=iand(state(k),s)
if(da.eq.0.or.da.eq.s) go to 20
ix=ieor(state(k),s)
a=iand(ix,rght)
b=iand(ix,lft)/hfbit
j=state2(1,a)+state2(2,b)
corr=corr+x1(k)*x1(j)
20 continue
sxy1(bond)=corr/2.d0
10	continue
!---------------------!

do 1 bond =1,nb1
!<sx_k sx_l+sy_k sy_l> nearest neighbor spin-2 correlation function
isite1=npair2(2*bond-1)-1
isite2=npair2(2*bond )-1
is1=2**isite1
is2=2**isite2
is=is1+is2
!------------------------!
corr=0.d0
do 2 k=1,nstates
da1=iand(state(k),is)
if(da1.eq.0.or.da1.eq.is) go to 2
ix1=ieor(state(k),is)
ia=iand(ix1,rght)
ib=iand(ix1,lft)/hfbit
j=state2(1,ia)+state2(2,ib)
corr=corr+x1(k)*x1(j)
2 continue
sxy2(bond)=corr/2.d0
1	continue
!--------------------------!

do 11 bond =1,nb1
!<sz_i sz_j> nearest neighbor spin-2 correlation function
site1=npair1(2*bond-1)-1
site2=npair1(2*bond )-1
s1=2**site1
s2=2**site2
s=s1+s2
!-----------------------!
corr=0.d0
do 12 k=1,nstates
da=iand(state(k),s)
if(da.eq.0.or.da.eq.s)then
factor=1.d0
else
factor=-1.d0
endif
corr=corr+factor*x1(k)**2
12 continue
szz1(bond)=corr/4.d0
11	continue
!-------------------------!


do 3 bond =1,nb1
!<sz_k sz_l> nearest neighbor spin-2 correlation function
isite1=npair2(2*bond-1)-1
isite2=npair2(2*bond )-1
is1=2**isite1
is2=2**isite2
is=is1+is2
!-----------------------!
corr=0.d0
do 4 k=1,nstates
da1=iand(state(k),is)
if(da1.eq.0.or.da1.eq.is)then
factor=1.d0
else
factor=-1.d0
endif
corr=corr+factor*x1(k)**2
4 continue
szz2(bond)=corr/4.d0
3	continue
!---------------------------!

!dimer-dimer corelation
!<sz_i sz_j sz_k sz_l> 
do 21 bond =1,nb1
site1=npair1(2*bond-1)-1
site2=npair1(2*bond )-1
s1=2**site1
s2=2**site2
s=s1+s2
!--------------------------!
isite1=npair2(2*bond-1)-1
isite2=npair2(2*bond )-1
is1=2**isite1
is2=2**isite2
is=is1+is2
corr=0.d0
do 22 k=1,nstates
da=iand(state(k),s)
da1=iand(state(k),is)
if(da==s.and.da1/=is.and.da1/=0)then
factor=-1.d0
elseif(da==0.and.da1/=is.and.da1/=0)then
factor=-1.d0
elseif(da/=s.and.da/=0.and.da1==is)then
factor=-1.d0
elseif(da/=0.and.da/=s.and.da1==0)then
factor=-1.d0
else
factor=1.d0
endif
corr=corr+factor*x1(k)**2
22 continue
szz3(bond)=corr/16.d0
21	continue
!-------------------------!

!----------------------------!
!(<sx_i sx_j+sy_i sy_j)(sx_k sx_l+sy_k sy_l)> 
do 13 bond =1,nb1
site1=npair1(2*bond-1)-1
site2=npair1(2*bond )-1
s1=2**site1
s2=2**site2
s=s1+s2
!--------------------------!
isite1=npair2(2*bond-1)-1
isite2=npair2(2*bond )-1
is1=2**isite1
is2=2**isite2
is=is1+is2
!--------------------------!
corr=0.d0
do 14 k=1,nstates
da=iand(state(k),s)
da1=iand(state(k),is)
if(da==0.or.da==s.or.da1==0.or.da1==is)goto 14
ix1=ieor(state(k),is)
ia=iand(ix1,rght)
ib=iand(ix1,lft)/hfbit
j=state2(1,ia)+state2(2,ib)
!-----------------------!
ix=ieor(state(j),s)
a=iand(ix,rght)
b=iand(ix,lft)/hfbit
j=state2(1,a)+state2(2,b)
corr=corr+x1(k)*x1(j)
14 continue
sxy3(bond)=corr/4.d0
13	continue
!----------------------------!

!----------------------------!
!(<sx_i sx_j+sy_i sy_j)sz_k sz_l> 
do 7 bond =1,nb1
site1=npair1(2*bond-1)-1
site2=npair1(2*bond )-1
s1=2**site1
s2=2**site2
s=s1+s2
!--------------------------!
isite1=npair2(2*bond-1)-1
isite2=npair2(2*bond )-1
is1=2**isite1
is2=2**isite2
is=is1+is2
!--------------------------!
corr=0.d0
do 6 k=1,nstates
da1=iand(state(k),is)
if(da1==0.or.da1==is)goto 6
ix1=ieor(state(k),is)
ia=iand(ix1,rght)
ib=iand(ix1,lft)/hfbit
j=state2(1,ia)+state2(2,ib)
da=iand(state(j),s)
if(da==0.or.da==s)then
factor=1.d0
else
factor=-1.d0
endif
corr=corr+factor*x1(k)*x1(j)
6 continue
szz4(bond)=corr/8.d0
7	continue
!------------------------------!

!----------------------------!
!<(sx_k sx_l+sy_k sy_l)sz_i sz_j> 
do 8 bond =1,nb1
site1=npair1(2*bond-1)-1
site2=npair1(2*bond )-1
s1=2**site1
s2=2**site2
s=s1+s2
!--------------------------!
isite1=npair2(2*bond-1)-1
isite2=npair2(2*bond )-1
is1=2**isite1
is2=2**isite2
is=is1+is2
!--------------------------!
corr=0.d0
do 9 k=1,nstates
da=iand(state(k),s)
if(da==0.or.da==s)goto 9
ix=ieor(state(k),s)
a=iand(ix,rght)
b=iand(ix,lft)/hfbit
j=state2(1,a)+state2(2,b)
da1=iand(state(j),is)
if(da1==0.or.da1==is)then
factor=1.d0
else
factor=-1.d0
endif
corr=corr+factor*x1(k)*x1(j)
9 continue
szz5(bond)=corr/8.d0
!-----------------------!
sxyz1(bond)=sxy1(bond)+szz1(bond)
sxyz2(bond)=sxy2(bond)+szz2(bond)
sxyz(bond)=sxy3(bond)+szz3(bond)+szz4(bond)+szz5(bond)
dimcorr(bond)=sxyz(bond)-sxyz2(bond)*sxyz1(bond)
!print*,'sxyz1(bond),sxyz2(bond),sxyz(bond)=',sxyz1(bond),sxyz2(bond),sxyz(bond),dimcorr(bond)
print*,site1+1,site2+1,isite1+1,isite2+1,dimcorr(bond),sxyz1(bond),sxyz2(bond),sxyz(bond)
write(1,*)site1+1,site2+1,isite1+1,isite2+1,dimcorr(bond)
8	continue
end subroutine dim2corr
!---------------------

subroutine spinsquared(H,nstates,state,state2)
use numbr_spn
implicit none
real*8::H(nstates,nstates)
integer*4::k,j,di,da,bond,ix,nstates,site1,site2,s1,s2,s,ipair1(2*11)
integer*4::state(nstates),a,b,state2(2,0:2**n-1),ipair(2*11),j1,new,k1,i

data ipair1/1,2, 2,3, 3,4, 4,5, 5,6, 6,7, 7,8, 8,9, 9,10, 10,11, 11,12/


H(:,:)=0.d0
do  13 bond =1,11
site1=ipair1(2*bond-1)-1
site2=ipair1(2*bond )-1
s1=2**site1
s2=2**site2
s=s1+s2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do 13 k=1,nstates
da=iand(state(k),s)
H(k,k)=(dfloat(nup)-dfloat(n)/2.d0)**2+dfloat(n)/2.d0
if(da.eq.0.or.da.eq.s)then
!H(k,k)=(dfloat(nup)-dfloat(n)/2.d0)**2+dfloat(n)/2.d0
!H(k,k)=H(k,k)+0.5d0+real(n)*0.75d0/real(11)
else
ix=ieor(state(k),s)
a=iand(ix,rght)
b=iand(ix,lft)/hfbit
j=state2(1,a)+state2(2,b)
!H(k,k)=H(k,k)-0.5d0+real(n)*0.75d0/real(11)
H(k,j)=H(k,k)+1.d0
end if
13 continue
end subroutine spinsquared
!!!!!!!!!!!!!!!!!!!!!!!!

subroutine transform(nstates,H,v1,dia)
implicit none
integer :: i,nstates
real(8) :: H(nstates,nstates),v1(nstates,nstates),dia(nstates)
H=matmul(H,v1)
H=matmul(transpose(v1),H)
do i=1,nstates
dia(i)=H(i,i)
enddo
end subroutine transform
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------------!
!subroutine spinsquared(nu)
!--------------------------!
!use system; implicit none
!integer :: i,j,a,b,sa,sb,nu
!mat(:,:)=0.d0
!do a=1,nst
!sa=state(a)
!mat(a,a)=(dfloat(nu)-dfloat(nn)/2.d0)**2+dfloat(nn)/2.d0
!do i=0,nn-1
!do j=i+1,nn-1
!if (btest(sa,i).neqv.btest(sa,j)) then 
!sb=ieor(sa,2**i+2**j)
!call findstate(sb,b)
!mat(a,b)=1.d0    
!endif
!enddo
!enddo
!enddo
!end subroutine spinsquared




