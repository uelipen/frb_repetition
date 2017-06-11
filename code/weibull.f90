  integer, parameter :: n=1024*128,nr=512
  real, parameter :: wk=0.34
  real, dimension(n) :: delta,deltacum
  real, dimension(n*nr) :: rho
  complex, dimension(n*nr) :: crho

  call random_number(delta)

  delta=(-log(1-delta))**(1/wk)

  rho=0
  deltacum(1)=delta(1)
  do i=2,n
     deltacum(i)=deltacum(i-1)+delta(i)
  enddo  
  do i=1,n
     ix=nint(deltacum(i)*nr*n/deltacum(n))
     if (ix<1 .or. ix>n*nr) cycle
     rho(ix)=rho(ix)+1
  end do
  
  write(*,*) deltacum(n)/(n),sum(rho)!,deltacum
  crho=rho/sum(rho)-1
  call four1(crho,n*nr,1)
  crho=crho*conjg(crho)
  crho(1)=0
  call four1(crho,n*nr,-1)
  rho=crho!/(n*nr)
  open(10,file='corr.dat')
  do i=1,n/2
     write(10,*) i-1,rho(i)
  end do
  open(10,file='corrbin.dat')
  i=0
  do
     if (2**(i+1)>n*nr) exit
     write(10,*) 2.**(i+0.5)/nr,sum(rho(2**i:2**(i+1)))/(2**i)
     i=i+1
  end do
end program
