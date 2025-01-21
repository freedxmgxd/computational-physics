!Curso Fisica Computacional
!Prof: Pedro Alves da Silva Autreto

!Metodo de Euler - pendulo simples

program  pendulo

! Declaracao de variaveis
implicit none
real*8 :: pi, dt, ti,tf,tmax
real*8 :: theta, omega, thetanew, omeganew,theta0
real*8 :: g, L, energia

!Arquivos de saida
open (11, FILE='omega.dat', STATUS='UNKNOWN')
open (12, FILE='theta.dat', STATUS='UNKNOWN')
open (13, FILE='energia.dat', STATUS='UNKNOWN')


!Inicializacao dos parametros
g = 10.0d0
L = 1.0d0
pi = 4*atan(1.0)
dt = 0.1
tmax = 10.0 

! Condicoes iniciais
theta0 = 10.0*pi/180.0
theta = theta0
omega = 0.0d0
ti = 0.0d0

energia = 0.0d0

!Calculo do angulo (theta) e da velocidade angular (omega) via 
! Metodo de Euler
do while (ti <= tmax )
  tf = ti + dt


  thetanew = theta + omega*(tf-ti)
  omeganew = omega  -(g/L)*sin(theta)*(tf-ti)
 
  energia=0.5*L**2*omeganew**2 - g*L*cos(thetanew)
   
  write(11,"(2F16.8)") tf,  omeganew 
  write(12,"(2F16.8)") tf, thetanew
  write(13,"(2F16.8)") tf,  energia
  ti =tf
  theta = thetanew
  omega = omeganew

end do

end program pendulo
