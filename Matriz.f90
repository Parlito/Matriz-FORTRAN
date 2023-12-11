program matriz

   implicit none

   double precision :: dist, x(13), y(13), z(13)
   Real*8 :: E(13), gama, H(13,13), Ef
   integer :: i, j, nat 
   character*1:: atom
   
   DOUBLE PRECISION,ALLOCATABLE:: WORK(:)
   INTEGER :: INFO,LWORK

   OPEN(UNIT=1,FILE='triangulene.xyz')

   READ(1,*) nat
   READ(1,*)

   do i = 1, nat
   
    READ(1,*) atom,x(i),y(i),z(i)
    
   end do
   
   CLOSE(UNIT=1)
   
   H=0.0
   
   gama=3.0
   
   do i = 1,13
   
     do j = 1,13  
     
     IF(i.eq.j)THEN
     
     H(i,j)=0.0
     CYCLE
     
     END IF
     
     dist=SQRT((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
     
     IF (dist.LT.1.410D0) THEN
     
      H(i,j) = gama
   
     END IF
     
    end do

   end do
   
   OPEN(UNIT=1,FILE="matriz.dat")
   do i=1,nat
   WRITE(1,*)(H(i,j),j=1,nat)
   end do
   CLOSE(UNIT=1)
    
      !Diagonalizando o Hamiltoniano
      ALLOCATE(WORK(1))
      LWORK=-1
      CALL DSYEV('V','U',nat,H,nat,E,WORK,LWORK,INFO)
      LWORK=WORK(1)
      DEALLOCATE(WORK)
      ALLOCATE(WORK(LWORK))
      CALL DSYEV('V','U',nat,H,nat,E,WORK,LWORK,INFO)

      ! Nível de Fermi
      Ef=(E(nat/2)+E(1+nat/2))/2.0D0
      Print*, E(nat/2),E(1+nat/2),nat/2,1+nat/2
      Print*,'O nível de Fermi é ',Ef
      
      ! Plotando os níveis de Energia
      OPEN(UNIT=1,FILE='energia.dat')
      DO i=1,nat
        WRITE(1,*) 0.0, E(i)
        WRITE(1,*) 1.0, E(i)
        WRITE(1,*) 
      END DO
      CLOSE(UNIT=1)
   
end program matriz
