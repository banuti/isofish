MODULE moddiffuse

!USE flowprops

CONTAINS

      SUBROUTINE diffusion !(xmue,delt,qinf)
	USE flowprops
	USE boundaries

!diffusion calculation

!      PARAMETER (imax=189,jmax=71,kmax=101)     

      DIMENSION  q_2(imax,jmax,kmax,3) !, qinf(3)

!      COMMON /level/ F(imax,jmax,kmax)
!      COMMON /grid/ delx,dely,delz
!      COMMON /velo/ q(imax,jmax,kmax,3)
!	COMMON /mode/ nswit

      WRITE(*,*) "check: diffusion"

	  q_2=0.
!------------------
!with diffusion

!inner domain
      DO 40 i = 2,imax-1
        DO 30 j = 2,jmax-1          
          DO 20 k = 2,kmax-1
            IF (F(i,j,k) .GT. 0.) THEN

 !             DO 10 l = 1,3
 !              b_1 = (q(i+1,j,k,l)-2*q(i,j,k,l)		&
 !    &                +q(i-1,j,k,l))/delx**2			
 !               b_2 = (q(i,j+1,k,l)-2*q(i,j,k,l)		&
 !    &                +q(i,j-1,k,l))/dely**2			
 !               b_3 = (q(i,j,k+1,l)-2*q(i,j,k,l)		&
 !    &                 +q(i,j,k-1,l))/delz**2
 !               q_2(i,j,k,l) = q(i,j,k,l)+delt*xmue*(b_1+b_2+b_3)
!10            CONTINUE

q_2(i,j,k,:) = 	q(i,j,k,:)	&
				+delt*xmue*(	&
             (q(i+1,j,k,:)-2*q(i,j,k,:)+q(i-1,j,k,:))/delx**2+	&			
             (q(i,j+1,k,:)-2*q(i,j,k,:)+q(i,j-1,k,:))/dely**2+	&			
             (q(i,j,k+1,:)-2*q(i,j,k,:)+q(i,j,k-1,:))/delz**2		)
                

			
			END IF
20        CONTINUE
30      CONTINUE
40    CONTINUE
      
!boundary
      CALL boundary_num(q_2) !(imax,jmax,kmax,q_2,qinf)

!q = q_2
	              q = q_2     
!	DO 3000 i = 1,imax
!        DO 2000 j = 1,jmax
!          DO 1000 k = 1,kmax
!            DO 999 l = 1,3
!              q(i,j,k,l) = q_2(i,j,k,l)
!999         CONTINUE
!1000      CONTINUE
!2000    CONTINUE
!3000  CONTINUE
!   
!no slip
	IF (nswit.NE. 3.) THEN
		CALL impnoslip
	END IF     
   
      RETURN
END SUBROUTINE

END MODULE