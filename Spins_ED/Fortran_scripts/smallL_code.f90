module numbr_spn
integer,parameter::n=12,nup=1,nb=12,szval=n/2-nup,hfbit=2**((n+1)/2)	&
&	,rght=2**((n+1)/2)-1,lft=ieor(2**n-1,rght)
end module numbr_spn
!---------------------------------!
program main
use numbr_spn
implicit none
integer, allocatable :: state(:)
integer, allocatable :: bckwindx(:)
real*8, allocatable :: H(:,:)
real*8, allocatable :: eig(:)
real*8, allocatable :: wk(:,:)
real*8, allocatable :: iwk(:)
real*8, allocatable :: v(:)
real*8, allocatable :: x(:)
integer*4::i,j,nstates,np,cont,ja,jb,bjp,ideclr,ne,nvec
integer*4::sz,state2(2,0:2**n-1),a,b,jc,k1,k,npair(2,nb)
real*8:: eps,E(4),hexpec,sxx(1),szz(1),bondwt(nb),zrtio(nb),wght
integer, external :: fncup
!-----------------------------------!
eps=1.d-13
nvec=1
ne=4
nstates =fncup(n,nup)
ideclr=nstates
allocate (state(nstates))
allocate (H(nstates,nstates))
allocate (wk(ideclr,8))
allocate (iwk(ideclr))
allocate (v(ideclr))
allocate (x(ideclr))
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
call hamiltonian(H,nstates,state,state2)
!do k=1,nstates
!write(*,*)(H(k,j),j=1,nstates)
!enddo
!-------------------------------!
!allocate(eig(nstates))
!call sym(H,eig,nstates)
!print*,(i,eig(i),i=1,4)
!do i=1,nstates
!write(*,*)(H(i,j),j=1,1)
!end do
!-----------------------!
call diag(H,nstates,nstates,E,v,ne,nvec,eps,wk,iwk)
print 100,E
100 format(/' [Eigenvalues]  '/2x,4f14.8)
!------------------!
call hamiltonian(H,nstates,state,state2)
!----------------------------!
call check3(H,nstates,nstates,v,wk,Hexpec)
!-------------------------!
npair(1,1)=1
npair(2,1)=2
call xcorr (1,npair,nstates,state,state2,sxx,v)
!------------------------!
call zcorr (1,npair,nstates,state,state2,szz,v)
!--------------------------!
print 130,sxx,szz
130  format('[Nerstnighbrcorrelationfunctions]' 'sxx :',d18.10,',szz :',d18.10)
!-------------------!
!---------------------!
stop
end
!-----------------!

subroutine hamiltonian(H,nstates,state,state2)
use numbr_spn
implicit none
real*8::H(nstates,nstates),diag,wght,bondwt(nb),zrtio(nb)
integer*4::k,j,di,da,bond,ix,nstates,site1,site2,s1,s2,s,cont
integer*4::state(nstates),a,b
integer*4::state2(2,0:2**n-1),ipair(2,nb),j1,new,k1,i
data bondwt/nb*-1/
data zrtio/nb*1/
!----------!
H(:,:)=0
do  30 bond =1,nb
!-----------------------!
!for square lattice with periodic boundary condition
call npair_chain (ipair) !for 1d spin chain 
!call npair (pair)! for 2d square lattice
s1=2**(ipair(1,bond)-1)
s2=2**(ipair(2,bond)-1)
s=s1+s2
!--------------------------!
wght=bondwt(bond)
diag=0.5d0*wght*zrtio(bond)
!--------------------------!
do 30 k=1,nstates
da=iand(state(k),s)
if(da.eq.0.or.da.eq.s)then
H(k,k)=H(k,k)-0.5d0*diag
else
H(k,k)=H(k,k)+0.5d0*diag
ix=ieor(state(k),s)
a=iand(ix,rght)
b=iand(ix,lft)/hfbit
j=state2(1,a)+state2(2,b)
H(k,j)=-0.5d0*wght
end if
30 continue
end subroutine hamiltonian
!---------------------------!

!-----------------------------!
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
!----------------------------------------------------!
subroutine sym(H,eig,n)
implicit none
integer n,lwork,inf
real*8  H(n,n),eig(n),work(max(1,3*n-1))
lwork=3*n-1
call dsyev('V','U',n,H,n,eig,work,lwork,inf)
end subroutine sym

!-------------------------------!



!-----------------------------!
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
integer*4::i,j,k,nstates,ideclr,k2
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
print 200
200  format('Information from check3')
print 210,prd
210  format(' <x*H*x> =',1pd16.8)
print 220
220  format(' H*x(j)/x(j) (j=min(nstates/3,13),nstates,max(1,nstates/20))')
print 230,(v(k2)/x(k2),k2=min(nstates/3,13),nstates,max(1,nstates/20))
230  format(4d18.9)
print 240
240  format('-----------')
return
end subroutine check3
!-----------------------------------!

subroutine xcorr(nb1,npair,nstates,state,state2,sxx,x)
use numbr_spn
implicit none
real*8::H(nstates,nstates),sxx(nb1),x(nstates),corr
integer*4::k,j,di,da,bond,ix,nstates,site1,site2,s1,s2,s
integer*4::state(nstates),a,b,nb1
integer*4::state2(2,0:2**n-1),npair(2,nb1)
!----------------------------------------!
do 10 bond =1,nb1
!------------------!
s1=2**(npair(1,bond)-1)
s2=2**(npair(2,bond)-1)
corr=0.d0
s=s1+s2
!-------------------!
do 20 k=1,nstates
da=iand(state(k),s)
if(da.eq.0.or.da.eq.s) go to 20
ix=ieor(state(k),s)
a=iand(ix,rght)
b=iand(ix,lft)/hfbit
corr=corr+x(k)*x(state2(1,a)+state2(2,b))
20 continue
sxx(bond)=corr/4.d0
!print*,sxx(bond)
10	continue
return
end subroutine xcorr
!------------------------------------------------!
subroutine zcorr(nb1,npair,nstates,state,state2,szz,x)
use numbr_spn
implicit none
real*8::H(nstates,nstates),szz(nb1),x(nstates),corr
integer*4::k,j,di,da,bond,ix,nstates,site1,site2,s1,s2,s,factor
integer*4::state2(2,0:2**n-1),npair(2,nb1),state(nstates),a,b,nb1
!---------------------------!
do 10 bond =1,nb1
!------------------!
s1=2**(npair(1,bond)-1)
s2=2**(npair(2,bond)-1)
corr=0.d0
s=s1+s2
!-----------------------!
do 20 k=1,nstates
da=iand(state(k),s)
if(da.eq.0.or.da.eq.s)then
factor=1.d0
else
factor=-1.d0
endif
corr=corr+factor*x(k)**2
20 continue
szz(bond)=corr/4.d0
!print*,szz(bond)
10	continue
return
end subroutine zcorr
!----------------------!

!--------------------------!
!-----------------------!
subroutine npair (ipair)
use numbr_spn
implicit none
integer :: sa,x1,x2,y1,y2,nn,lx,ly,i,ipair(2,nb)
lx=4
ly=4
nn=lx*ly
do y1=0,ly-1
do x1=0,lx-1
sa=1+x1+y1*lx
x2=mod(x1+1,lx)
y2=y1
ipair(1,sa)=sa
ipair(2,sa)=1+x2+y2*lx
x2=x1
y2=mod(y1+1,ly)
ipair(1,sa+nn)=sa
ipair(2,sa+nn)=y2*lx+x2+1  
enddo
enddo
!do i=1,nb
!print*,ipair(1,i),ipair(2,i)
!enddo
end subroutine npair
!------------------------!


!----------------------------!
subroutine npair_chain (ipair)
use numbr_spn
integer*4::bond,bond1,ipair(2,nb)
do 20 bond =1,nb
!------------------!
ipair(1,bond)=bond
!--------------------!
bond1=bond+1
if(bond1>nb)bond1=bond1-nb
ipair(2,bond)=bond1
!print*,ipair(1,bond),ipair(2,bond)
20	continue
end subroutine npair_chain

