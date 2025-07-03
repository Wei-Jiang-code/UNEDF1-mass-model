program drip
IMPLICIT NONE
INTEGER,PARAMETER :: db=SELECTED_REAL_KIND(12,100)
INTEGER(db) ::i,v
real(db)::n,z,a
real(db) ::b,b1,bc=4.41d13,me=0.51099889d0,m,pi=3.14159265358979D0,cbbc=3.40665d-3,c=-1.44423,alfa=1.0/137.036
real(db)::kpl,g,k,pl,bn,pdirp,udrip,pe,pdirp0,p0
real(db),external :: t

b=8d16
z=39d0
n=97d0
bn=915.83d0
m=938.272088d0*z+939.56563d0*n-bn
b1=b/bc
a=n+z
!计算晶格压强系数、p系数
kpl=-1.0/3.0*cbbc*z**(2.0/3.0)*b1*me**4/(2.0*pi**2)
k=b1*me**4/(2.0*pi**2)

udrip=(-m+a*939.56563d0)/z

pdirp=b1*me**2.0*udrip**2.0*(1.0-1.0/3.0*c*alfa*z**(2.0/3.0)*(4.0*b1/pi**2.0)**(1.0/3.0)*(me/udrip)**(2.0/3.0))/(4.0*pi**2.0)&
*1.3015949d-7
pdirp0=(udrip)**4.0*(1.0+4.0*c*alfa*z**(2.0/3.0)/(81d0*pi**2.0)**(1.0/3.0))**(-3.0)/(12*pi**2.0)*1.3015949d-7
p0=me**(4.0)*udrip**(4.0)/(12.0*pi**2.0)*1.3015949d-7
g=m/a+z*udrip/a
write(*,*)"pdrip=",pdirp,udrip,pdirp0!,p0
pause
end program

function t(x) result(t1)
 IMPLICIT NONE
 INTEGER,PARAMETER :: db=SELECTED_REAL_KIND(12,100)
 real(db)::x,t1
 t1=1.0/2.0*x*sqrt(1+x**2)-1.0/2.0*log(x+sqrt(1+x**2))
 end function t
