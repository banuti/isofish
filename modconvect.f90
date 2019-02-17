MODULE modconvect
!USE flowprops

CONTAINS

!***********************************************************************
SUBROUTINE convection !(delt,qinf)

CALL convectinterpolation
!CALL convectcentraldiff


END SUBROUTINE

!=======================================================================

SUBROUTINE convectcentraldiff
!second order central difference 

USE flowprops
USE boundaries

INTEGER ::i,j,k

DIMENSION q_1(imax,jmax,kmax,3) !, qinf(3)


WRITE(*,*) "check: convection (central diff)"




DO i=2,imax-1
  DO j=2,jmax-1
    DO k=2,kmax-1

!x-direction
q_1(i,j,k,1)=q(i,j,k,1)-delt*( q(i,j,k,1)*( q(i+1,j,k,1)-q(i-1,j,k,1) )/(2*delx) +  &
                                q(i,j,k,2)*( q(i,j+1,k,1)-q(i,j-1,k,1) )/(2*dely) + &
                                q(i,j,k,3)*( q(i,j,k+1,1)-q(i,j,k-1,1) )/(2*delz)   )




!y-direction
q_1(i,j,k,2)=q(i,j,k,2)-delt*( q(i,j,k,1)*( q(i+1,j,k,2)-q(i-1,j,k,2) )/(2*delx) +   &
                                q(i,j,k,2)*( q(i,j+1,k,2)-q(i,j-1,k,2) )/(2*dely) +  &
                                q(i,j,k,3)*( q(i,j,k+1,2)-q(i,j,k-1,2) )/(2*delz)    )



!z-direction
q_1(i,j,k,3)=q(i,j,k,3)-delt*( q(i,j,k,1)*( q(i+1,j,k,3)-q(i-1,j,k,3) )/(2*delx) +   &
                                q(i,j,k,2)*( q(i,j+1,k,3)-q(i,j-1,k,3) )/(2*dely) +  &
                                q(i,j,k,3)*( q(i,j,k+1,3)-q(i,j,k-1,3) )/(2*delz)    )


    END DO
  END DO
END DO



!boundary
      CALL boundary_num(q_1) !(imax,jmax,kmax,q_1,qinf)

      q=q_1

!no slip
      CALL impnoslip

RETURN

END SUBROUTINE

!-----------------------------------------------------------------------


SUBROUTINE convectinterpolation


	USE flowprops
	USE boundaries

!
! convection calculation

      REAL nuex, nuey, nuez

      DIMENSION q_1(imax,jmax,kmax,3) !, qinf(3)

      WRITE(*,*) "check: convection (interpolation)"

!------------------
!with convection

      DO 40 i = 2,imax-1
        DO 30 j = 2,jmax-1
          DO 20 k = 2,kmax-1
            IF (F(i,j,k) .GT. 0.) THEN

!sicherstellen, dass Geschwindigkeit in der richtigen angrenzenden Zelle 
!interpoliert wird
            IF (q(i,j,k,1) .GT. 0.) THEN 
              idist = -1
            ELSE
              idist = 1
            END IF
            IF (q(i,j,k,2) .GT. 0.) THEN
              jdist = -1     
            ELSE
              jdist = 1
            END IF
            IF (q(i,j,k,3) .GT. 0.) THEN
              kdist = -1
            ELSE
              kdist = 1
            END IF

            nuex = ABS(q(i,j,k,1))*delt/delx
            nuey = ABS(q(i,j,k,2))*delt/dely
            nuez = ABS(q(i,j,k,3))*delt/delz


!ueberpruefen, ob Interpolationspunkt in benachbarter Zelle liegt    
            IF(nuex .GT. 1. .OR. nuey .GT. 1. .OR. nuez .GT. 1.) THEN
              WRITE(*,*) "interpolationpoint outside neighbourcell"
              WRITE(*,*) nuex, nuey, nuez
            END IF

!velocity after convection / inner domain
!trilinear interpolation
            DO 10 l = 1,3
              q_1(i,j,k,l) = (1-nuex)*(1-nuey)*(1-nuez)*q(i,j,k,l)      &   
     &        +  nuex  *(1-nuey)*(1-nuez)*q(i+idist,j,k,l)				&
     &        +  nuex  *  nuey  *(1-nuez)*q(i+idist,j+jdist,k,l)		&
     &        +(1-nuex)*  nuey  *(1-nuez)*q(i,j+jdist,k,l)				&
     &        +(1-nuex)*(1-nuey)*  nuez  *q(i,j,k+kdist,l)				&
     &        +  nuex  *(1-nuey)*  nuez  *q(i+idist,j,k+kdist,l)		&
     &        +  nuex  *  nuey  *  nuez  *q(i+idist,j+jdist,k+kdist,l)	&
     &        +(1-nuex)*  nuey  *  nuez  *q(i,j+jdist,k+kdist,l)
10          CONTINUE
            END IF
20        CONTINUE
30      CONTINUE
40    CONTINUE       

!boundary
      CALL boundary_num(q_1) !(imax,jmax,kmax,q_1,qinf)

!q = q_1
	q = q_1
!      DO 3000 i = 1,imax
!        DO 2000 j = 1,jmax
!          DO 1000 k = 1,kmax
!            DO 999 l = 1,3
!              q(i,j,k,l) = q_1(i,j,k,l)
!999         CONTINUE
!1000      CONTINUE
!2000    CONTINUE
!3000  CONTINUE
!  
!no slip
      CALL impnoslip
!   
      RETURN
      END SUBROUTINE


END MODULE
