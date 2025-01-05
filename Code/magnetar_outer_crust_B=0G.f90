program atomic_Gibbs
IMPLICIT NONE
!定义变量
INTEGER,PARAMETER :: db=SELECTED_REAL_KIND(12,100)
INTEGER(db) ::i,j,l,ip,rdex=0,rmin=0,count1=0
real(db) ::p,bc=4.41d13,me=0.51099889d0,u,u1,u2,pi=3.14159265358979D0,cbbc=3.40665d-3
real(db)::k,pion=0
real(db),external :: t
type d
real(db) ::z
real(db) ::n
real(db) ::bn
real(db) ::rp
real(db) ::rn
real(db) ::rt
real(db) ::beta
real(db) ::gama
end type d
type struct
real(db) ::z
real(db) ::n
real(db) ::pmax
real(db) ::p1min
real(db) ::nmin
real(db) ::nmax
end type struct
include'sunedf1_0.data'
type(d),dimension(nmass) ::mass1
type(struct),allocatable,save :: str(:)
real(db),ALLOCATABLE ::f1(:),f2(:),pl1(:),pl2(:),pe1(:),pe2(:),xn1(:),ne1(:),xr1(:),xn2(:),ne2(:),xr2(:)
real(db),allocatable ::m(:),a(:),kpl(:),nb(:),ne(:),ue(:),g(:),pl(:),zmin(:),amin(:),nbmin(:),gmin(:)
ALLOCATE(f1(10000),f2(10000),pl1(10000),pl2(10000),pe1(10000),pe2(10000),xn1(10000),ne1(10000),xr1(10000),&
xn2(10000),ne2(10000),xr2(10000))
ALLOCATE(a(10000),m(10000),kpl(10000),nb(10000),ne(10000),ue(10000),g(10000),pl(10000),zmin(10000),&
amin(10000),nbmin(10000),gmin(10000))
ALLOCATE(str(50))
!定义质量数组模块
do i=1,nmass
mass1(i)=pmass(i)
enddo



open(unit=5, file='sunedf1_0.0E+0.txt', status='unknown', action='write', position='append')
iteration5: do rdex=0,9

  rmin=10**rdex


  iteration4: do ip=1,1000,1

p=rmin*1.0d-9+rmin*ip*1.0d-11
!P=5.01E-4
k=me**4/(8.0*pi**2)
!初始化迭代
u=5d0
u1=50d0

!***********************************************************
!开始迭代
Iteration2: do l=1,nmass
 !计算晶格压强系数、质量
  a(l)=mass1(l)%n+mass1(l)%z
  m(l)=938.272088d0*mass1(l)%z+939.56563d0*mass1(l)%n-mass1(l)%bn
  kpl(l)=-1.0/3.0*cbbc*mass1(l)%z**(2.0/3.0)


Iteration: do i=1,100

xn1(i)=sqrt(u**2.0-me**2.0)/me
ne1(i)=(me*xn1(i))**3.0/(3.0*pi**2.0)
xr1(i)=(3*pi**2*ne1(i))**(1.0/3.0)/me
    pe1(i)=k*t(xr1(i))
    pl1(i)=kpl(l)*ne1(i)*sqrt(u**2.0-me**2.0)
  
f1(i)=(pe1(i)+pl1(i))*1.3015949d-7-p


xn2(i)=sqrt(u1**2.0-me**2.0)/me
ne2(i)=(me*xn2(i))**3.0/(3.0*pi**2.0)
xr2(i)=(3*pi**2*ne2(i))**(1.0/3.0)/me
    pe2(i)=k*t(xr2(i))
    pl2(i)=kpl(l)*ne2(i)*sqrt(u1**2.0-me**2.0)
f2(i)=(pe2(i)+pl2(i))*1.3015949d-7-p

u2=(u*f2(i)-u1*f1(i))/(f2(i)-f1(i))

 if(abs(u2-u1)<1.0d-10) then
 !求电子密度、晶格压强、重子密度ne,pl,nb
 ue(l)=u1
  ne(l)=ne2(i)*1.3015949d-7
  nb(l)=a(l)*ne(l)/mass1(l)%z
  pl(l)=pl1(i)*1.3015949d-7
  !write(*,*)  nb(l),i,ue(l)
  exit iteration
 ENDIF
 
 u=u1
 u1=u2
enddo iteration
!*************************************************************
enddo Iteration2

do l=1,nmass
  g(l)=m(l)/a(l)+mass1(l)%z/a(l)*(ue(l)+4.0d0*pl(l)/ne(l))
  !write(*,*) "g=",g(l),mass1(l)%z,a(l)
enddo


!求最小吉布斯自由能
gmin(ip)=g(1)
zmin(ip)=mass1(1)%z
amin(ip)=a(1)
nbmin(ip)=nb(1)


gibssmin:do i=1,nmass
  if(g(i)<gmin(ip)) then
  gmin(ip)=g(i)
  zmin(ip)=mass1(i)%z
  amin(ip)=a(i)
  nbmin(ip)=nb(i)
  endif
enddo gibssmin

!write(*,*) nbmin(ip)
!将得到不同压强下的元素输出
if((zmin(ip)/=zmin(ip-1).or.amin(ip)/=amin(ip-1)).and.(zmin(ip)>10.and.zmin(ip-1)>10)) then
write(*,*) zmin(ip-1),amin(ip-1)-zmin(ip-1),p-rmin*1*1.0d-11,nbmin(ip-1),nbmin(ip)

    count1=count1+1
     str(count1)%z=zmin(ip-1)
     str(count1)%n=amin(ip-1)-zmin(ip-1)
     str(count1)%pmax=p-rmin*1*1.0d-11
     str(count1)%p1min=p
     str(count1)%nmin=nbmin(ip-1)
     str(count1)%nmax=nbmin(ip)
   !  write(*,*)str(count1)


if(count1==1) then
write(5, '(2f6.1, 4es20.2)') str(count1)%z,str(count1)%n,pion,str(count1)%pmax,pion,str(count1)%nmin
else
write(5, '(2f6.1, 4es20.2)/') str(count1)%z,str(count1)%n,str(count1-1)%p1min,str(count1)%pmax,str(count1-1)%nmax,str(count1)%nmin
endif


endif

enddo iteration4

enddo iteration5
close(5)




  DEALLOCATE(f1,f2,pe1,pe2,pl1,pl2,xn1,ne1,xr1,xn2,ne2,xr2,m,a,kpl,nb,ne,ue,g,pl,gmin,amin,zmin,nbmin,str)
  pause
  end program atomic_Gibbs

 !函数
 function t(x) result(t1)
 IMPLICIT NONE
 INTEGER,PARAMETER :: db=SELECTED_REAL_KIND(12,100)
 real(db)::x,t1
 t1=x*(2.0/3.0*x**2.0-1.0)*sqrt(1.0+x**2.0)+log(x+sqrt(1.0+x**2.0))
 end function t

