module global
implicit none
INTEGER,PARAMETER:: db=SELECTED_REAL_KIND(12,100)
real(db),save ::b,b1,bc=4.41d13,me=0.51099889d0,pi=3.14159265358979D0,cbbc=3.40665d-3
real(db),save ::k,pion=0
real(db),allocatable ::m(:),a(:),kpl(:)
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
include'unedf1_a16_8.data'
type(struct),allocatable,save :: str(:)
end module global

!********************************************************************************************************************************
!main:includes OpenMp parallel computing
!********************************************************************************************************************************
program atomic_Gibbs
use global
IMPLICIT NONE
!定义变量
integer ::i,j,l
INTEGER ::i1,ip,rdex,count1=0,vmax_out,gmin_index_local,gmin_index
real(db) ::p,rmin,ne_out,u_out,pl_out,nb_out,gmin_local
character(len=20) :: filename
real(db),allocatable ::nb2(:),ne2(:),u3(:),pl(:),g(:),zmin(:),amin(:),nbmin(:),gmin(:),umin(:)
integer,allocatable ::vmax3(:),vmax4(:)
ALLOCATE(str(50))
ALLOCATE(a(10000),m(10000),pl(10000),kpl(10000),nb2(10000),ne2(10000),u3(10000),g(10000),zmin(10000),&
amin(10000),nbmin(10000),gmin(10000),vmax3(10000),vmax4(10000),umin(10000))

!Calculating the critical magnetic field
b=8.0d16
b1=b/bc
k=b1*me**4/(2.0*pi**2)


!$omp parallel do
  do l=1,nmass
    a(l)=pmass(l)%n+pmass(l)%z
    m(l)=938.272088d0*pmass(l)%z+939.56563d0*pmass(l)%n-pmass(l)%bn
    kpl(l)=-1.0/3.0*cbbc*pmass(l)%z**(2.0/3.0)*(3*pi**2)**(1.0/3.0)*(b1*me**3/(2.0*pi**2))**(4.0/3.0)
  enddo
!$omp end parallel do
!设置循环步长


iteration5: do rdex=0,9

  rmin=10**rdex


  iteration4: do ip=1,1000,1

  p=rmin*1.0d-9+rmin*ip*1.0d-11

!同时计算几个核的自由能
!$OMP PARALLEL DO SCHEDULE(STATIC)DEFAULT(shared) private(l,ne_out,nb_out,u_out,pl_out,vmax_out) 
Iteration3: do l=1,nmass
  call newton_iter(l,p,ne_out,nb_out,u_out,pl_out,vmax_out)
  ne2(l)=ne_out
  nb2(l)=nb_out
  u3(l)=u_out
  pl(l)=pl_out
  vmax3(l)=vmax_out
  g(l)=m(l)/a(l)+pmass(l)%z/a(l)*(u3(l)+4.0d0*pl(l)/ne2(l))
enddo Iteration3
!$omp end parallel do

gmin(ip)=g(1)
gmin_local=g(1)
gmin_index_local = -1  
gmin_index = -1  

 ! OpenMP并行区域  
    !$omp parallel firstprivate(gmin_local, gmin_index_local) shared(g, gmin, gmin_index)  
  
    ! 线程内循环  
    !$omp do schedule(static)  
    do i = 2, nmass  
        if (g(i) < gmin_local) then  
            gmin_local = g(i)  
            gmin_index_local = i  
        end if  
    end do  
    !$omp end do  
  
    ! 同步点来合并结果  
    !$omp critical  
    if (gmin_local < gmin(ip)) then  
        gmin(ip) = gmin_local  
        gmin_index = gmin_index_local  
    end if  
    !$omp end critical  
  
    !$omp end parallel  

      zmin(ip)=pmass(gmin_index)%z
      amin(ip)=a(gmin_index)
      nbmin(ip)=nb2(gmin_index)
      vmax4(ip)=vmax3(gmin_index)
      umin(ip)=u3(gmin_index)


!将得到不同压强下的元素输出
if((zmin(ip)/=zmin(ip-1).or.amin(ip)/=amin(ip-1)).and.(zmin(ip)>10.and.zmin(ip-1)>10)) then
write(*,*) zmin(ip-1),amin(ip-1)-zmin(ip-1),p-rmin*1*1.0d-11,nbmin(ip-1),nbmin(ip),vmax4(ip-1)!,umin(ip-1)

    count1=count1+1
     str(count1)%z=zmin(ip-1)
     str(count1)%n=amin(ip-1)-zmin(ip-1)
     str(count1)%pmax=p-rmin*1*1.0d-11
     str(count1)%p1min=p
     str(count1)%nmin=nbmin(ip-1)
     str(count1)%nmax=nbmin(ip)
   !  write(*,*)str(count1)
    write(filename, '("unedf1_", es7.1, ".txt")') b
open(unit=5, file=filename, status='unknown', action='write', position='append')

if(count1==1) then
write(5, '(2f6.1, 4es20.2)') str(count1)%z,str(count1)%n,pion,str(count1)%pmax,pion,str(count1)%nmin
else
write(5, '(2f6.1, 4es20.2)/') str(count1)%z,str(count1)%n,str(count1-1)%p1min,str(count1)%pmax,str(count1-1)%nmax,str(count1)%nmin
endif
close(5)
endif
enddo iteration4

enddo iteration5


DEALLOCATE(a,m,pl,kpl,nb2,ne2,u3,g,zmin,amin,nbmin,gmin,vmax3,vmax4,umin,str)
pause
end program atomic_Gibbs
 
!***************************************************************************************************************************************
!Newton's method
!****************************************************************************************************************************************
subroutine newton_iter(l,p,ne_p,nb_p,u_p,pl_p,vmax_p)
use global
implicit none
integer,intent(in)::l
real(db),INTENT(in)::p
integer,intent(out)::vmax_p
real(db),INTENT(out)::ne_p,nb_p,u_p,pl_p
integer::vmax,vmax1,i,j,v,v2
real(db)::ne1,ne,vmax2,u,u1,u2,nb
real(db) ::pl1(20),pl2(20),pe1(20),pe2(20),f1(20),f2(20),f3(20),f4(20),pl3(20),pl4(20),pe4(20),pe3(20)
real(db),external :: t


  !初始化迭代
  vmax=200d0
  vmax1=600d0
  pe1(:)=0d0
  pl1(:)=0d0
  pe2(:)=0d0
  pl2(:)=0d0
  !******************************************************
  !开始迭代
  Iteration: do i=1,20
    do v=0,vmax
      if(v==0) then
      pe1(i)=k*t(sqrt(2.0*b1*vmax))
      pl1(i)=sqrt(2*b1*vmax)
      else
      pe1(i)=pe1(i)+k*2.0*(1.0+2.0*v*b1)*t(sqrt((2.0*b1*(vmax-v))/(1.0+2.0*v*b1)))
      pl1(i)=pl1(i)+2*sqrt(2*b1*(vmax-v))
      endif
    enddo
  f1(i)=(pe1(i)+kpl(l)*(pl1(i))**(4.0/3.0))*1.3015949d-7-p   

    do v=0,vmax1
      if(v==0) then
      pe2(i)=k*t(sqrt(2*b1*vmax1))
      pl2(i)=sqrt(2*b1*vmax1)
      else
      pe2(i)=pe2(i)+k*2*(1+2*v*b1)*t(sqrt(2*b1*(vmax1-v)/(1+2*v*b1)))
      pl2(i)=pl2(i)+2*sqrt(2*b1*(vmax1-v))
      endif
    enddo
  f2(i)=(pe2(i)+kpl(l)*(pl2(i))**(4.0/3.0))*1.3015949d-7-p

  vmax2=(vmax*f2(i)-vmax1*f1(i))/(f2(i)-f1(i))

  if(abs(vmax2-vmax1)<0.5) then
  !求电子密度、晶格压强、重子密度ne,pl,nb
  !vmax1=1
    do v=0,vmax1
      if(v==0) then
        ne=b1*me**(3)/(2*pi**2)*sqrt(2.0*b1*vmax1)*1.3015949d-7
      else
        ne=ne+b1*me**(3)/(2*pi**2)*2*sqrt(2.0*b1*(vmax1-v))*1.3015949d-7
      endif
    enddo
    !nb=a*ne/z
    !write(*,*)  vmax2,vmax1
    exit iteration
  ENDIF
  v2=nint(vmax2)
  vmax=vmax1
  vmax1=v2
  enddo iteration
  !*************************************************************
  !初始化迭代
  u=30d0
  u1=40d0
  pe3(:)=0d0
  pl3(:)=0d0
  pe4(:)=0d0
  pl4(:)=0d0
  ne1=0
  !第二次迭代
  !************************************************************

  Iteration1: do j=1,20

  vmax=v2-1

    do v=0,vmax
      if(v==0) then
      pe3(i)=k*t(sqrt(u**2.0-me**2.0)/me)
      pl3(i)=sqrt(u**2.0-me**2.0)/me
      else
      pe3(i)=pe3(i)+k*2.0*(1.0+2.0*v*b1)*t(sqrt((u**2-me**2*(1+2.0*v*b1))/(me**2.0+2.0*v*b1*me**2.0)))
      pl3(i)=pl3(i)+2.0*sqrt(u**2.0-me**2.0*(1.0+2.0*v*b1))/me
      endif
    enddo
  f3(i)=(pe3(i)+kpl(l)*(pl3(i))**(4.0/3.0))*1.3015949d-7-p

    do v=0,vmax
      if(v==0) then
      pe4(i)=k*t(sqrt(u1**2.0-me**2.0)/me)
      pl4(i)=sqrt(u1**2.0-me**2.0)/me
      else
      pe4(i)=pe4(i)+k*2.0*(1.0+2.0*v*b1)*t(sqrt((u1**2.0-me**2.0*(1.0+2.0*v*b1))/(me**2.0+2.0*v*b1*me**2.0)))
      pl4(i)=pl4(i)+2.0*sqrt(u1**2.0-me**2.0*(1.0+2.0*v*b1))/me
      endif
    enddo
  f4(i)=(pe4(i)+kpl(l)*(pl4(i))**(4.0/3.0))*1.3015949d-7-p

  u2=(u*f4(i)-u1*f3(i))/(f4(i)-f3(i))

  if(abs(u2-u1)<1.0d-3) then
  !求电子密度、晶格压强、重子密度ne,pl,nb
    do v=0,vmax
      if(v==0) then
        ne1=b1*me**(3.0)/(2.0*pi**2.0)*sqrt(u1**2.0-me**2.0)/me*1.3015949d-7
      else
        ne1=ne1+b1*me**(3.0)/(2.0*pi**2.0)*2.0*sqrt(u1**2.0-me**2.0*(1.0+2.0*v*b1))/me*1.3015949d-7
      endif
    enddo
ne_p=ne1
nb_p=a(l)*ne1/pmass(l)%z
u_p=u2
pl_p=kpl(l)*(pl4(i))**(4.0/3.0)*1.3015949d-7
vmax_p=vmax
   ! write(*,*)  vmax
    exit iteration1
  ENDIF
  u=u1
  u1=u2
  enddo iteration1
end subroutine newton_iter
!***********************************************************************************************************************************
!***********************************************************************************************************************************

 !函数
 function t(x) result(t1)
 IMPLICIT NONE
 INTEGER,PARAMETER :: db=SELECTED_REAL_KIND(12,100)
 real(db),intent(in)::x
  real(db)::t1
 t1=1.0/2.0*x*sqrt(1+x**2)-1.0/2.0*log(x+sqrt(1+x**2))
 end function t

