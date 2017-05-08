  integer, parameter :: n=1000,nr=100
  real, parameter :: wk=0.34
  real, dimension(n) :: delta,deltacum
  real, dimension(n*nr/2) :: rho

  call random_number(delta)

  delta=(-log(1-delta))**(1/wk)

  rho=0
  deltacum(1)=delta(1)
  do i=2,n
     deltacum(i)=deltacum(i-1)+delta(i)
     ix=nint(deltacum(i)*nr)
     if (ix<1 .or. ix>n*nr/2) cycle
     rho(ix)=1
  enddo

  
  write(*,*) deltacum(n)/(n),sum(rho)

end program
