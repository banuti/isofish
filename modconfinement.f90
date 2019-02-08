MODULE modconfinement




CONTAINS


SUBROUTINE vort_conf(n) !(qinf,delt,n)
	USE flowprops
	USE boundaries
!
!vorticity confinement

      REAL neta_x, neta_y, neta_z, mageta, nenner                                                                  

!      PARAMETER (imax=189,jmax=71,kmax=101,
!     &           ipmax=188,jpmax=70,kpmax=100) 

      DIMENSION eta_c(ipmax,jpmax,kpmax),eta(imax,jmax,kmax), s(3)

!      COMMON /grid/ delx,dely,delz
!      COMMON /level/ F(imax,jmax,kmax)
!      COMMON /velo/ q(imax,jmax,kmax,3)
!      COMMON /conf/ eps(imax,jmax,kmax), epsns(imax,jmax,kmax)
!	COMMON /mode/ nswit

      WRITE(*,*) "check: vorticity confinement"

!------------------------
!with confinement

!1. omega and eta in cellcenter
      DO 30 ip = 1,ipmax
        DO 20 jp = 1,jpmax
          DO 10 kp = 1,kpmax
            omega(ip,jp,kp,1) = (q(ip,jp+1,kp,3)-q(ip,jp,kp,3)&
     &           +q(ip,jp+1,kp+1,3)-q(ip,jp,kp+1,3)&
     &           +q(ip+1,jp+1,kp+1,3)-q(ip+1,jp,kp+1,3)&
     &           +q(ip+1,jp+1,kp,3)-q(ip+1,jp,kp,3))/(4*dely)&
     &                         -(q(ip,jp,kp+1,2)-q(ip,jp,kp,2)&
     &           +q(ip,jp+1,kp+1,2)-q(ip,jp+1,kp,2)&
     &           +q(ip+1,jp+1,kp+1,2)-q(ip+1,jp+1,kp,2)&
     &           +q(ip+1,jp,kp+1,2)-q(ip+1,jp,kp,2))/(4*delz)
            omega(ip,jp,kp,2) = (q(ip,jp,kp+1,1)-q(ip,jp,kp,1)&
     &           +q(ip,jp+1,kp+1,1)-q(ip,jp+1,kp,1)&
     &           +q(ip+1,jp+1,kp+1,1)-q(ip+1,jp+1,kp,1)&
     &           +q(ip+1,jp,kp+1,1)-q(ip+1,jp,kp,1))/(4*delz)&
     &                         -(q(ip+1,jp,kp,3)-q(ip,jp,kp,3)&
     &           +q(ip+1,jp,kp+1,3)-q(ip,jp,kp+1,3)&
     &           +q(ip+1,jp+1,kp+1,3)-q(ip,jp+1,kp+1,3)&
     &           +q(ip+1,jp+1,kp,3)-q(ip,jp+1,kp,3))/(4*delx)
            omega(ip,jp,kp,3) = (q(ip+1,jp,kp,2)-q(ip,jp,kp,2)&
     &           +q(ip+1,jp,kp+1,2)-q(ip,jp,kp+1,2)&
     &           +q(ip+1,jp+1,kp+1,2)-q(ip,jp+1,kp+1,2)&
     &           +q(ip+1,jp+1,kp,2)-q(ip,jp+1,kp,2))/(4*delx)&
     &                        -(q(ip,jp+1,kp,1)-q(ip,jp,kp,1)&
     &           +q(ip,jp+1,kp+1,1)-q(ip,jp,kp+1,1)&
     &           +q(ip+1,jp+1,kp+1,1)-q(ip+1,jp,kp+1,1)&
     &           +q(ip+1,jp+1,kp,1)-q(ip+1,jp,kp,1))/(4*dely)
           absomega(ip,jp,kp) = -SQRT(max(0.000001,(omega(ip,jp,kp,1)**2&
     &           +omega(ip,jp,kp,2)**2+omega(ip,jp,kp,3)**2)))
		eta_c(ip,jp,kp) = absomega(ip,jp,kp)
10        CONTINUE
20      CONTINUE
30    CONTINUE

!2. eta = -abs(omega) on gridnode
!inner domain
      DO 60 i = 2,imax-1
        DO 50 j = 2,jmax-1
          DO 40 k = 2,kmax-1
            eta(i,j,k) = (eta_c(i-1,j-1,k-1)+eta_c(i-1,j,k-1)&
     &                   +eta_c(i,j,k-1)+eta_c(i,j-1,k-1)&
     &                   +eta_c(i-1,j-1,k)+eta_c(i-1,j,k)&
     &                   +eta_c(i,j,k)+eta_c(i,j-1,k))/8.
40        CONTINUE
50      CONTINUE
60    CONTINUE

!boundary
      CALL s_extrapolation(eta) !(imax,jmax,kmax,eta)
!
!3. correction 
      DO 100 ip = 2,ipmax-1
        DO 90 jp = 2,jpmax-1
          DO 80 kp = 2,kpmax-1


!SURFACE CONFINEMENT
	IF (F(ip,jp,kp) .LE. fieldRad) THEN

			neta_x = -Fnormp(ip,jp,kp,1)
			neta_y = -Fnormp(ip,jp,kp,2)
			neta_z = -Fnormp(ip,jp,kp,3)


			s(1) = neta_z*omega(ip,jp,kp,2)&
     &				-neta_y*omega(ip,jp,kp,3)
			s(2) = neta_x*omega(ip,jp,kp,3)&
     &				-neta_z*omega(ip,jp,kp,1)
			s(3) = neta_y*omega(ip,jp,kp,1)&
     &				-neta_x*omega(ip,jp,kp,2)

	!correction
            DO l = 1,3
              q(ip+1,jp+1,kp+1,l) = q(ip+1,jp+1,kp+1,l)&
     &              -epsns(ip+1,jp+1,kp+1)*delt*s(l)
              q(ip+1,jp+1,kp,l) = q(ip+1,jp+1,kp,l)&
     &              -epsns(ip+1,jp+1,kp)*delt*s(l)
              q(ip+1,jp,kp,l) = q(ip+1,jp,kp,l)&
     &              -epsns(ip+1,jp,kp)*delt*s(l)
              q(ip+1,jp,kp+1,l) = q(ip+1,jp,kp+1,l)&
     &              -epsns(ip+1,jp,kp+1)*delt*s(l)
              q(ip,jp+1,kp+1,l) = q(ip,jp+1,kp+1,l)&
     &              -epsns(ip,jp+1,kp+1)*delt*s(l)
              q(ip,jp+1,kp,l) = q(ip,jp+1,kp,l)&
     &              -epsns(ip,jp+1,kp)*delt*s(l)
              q(ip,jp,kp,l) = q(ip,jp,kp,l)&
     &              -epsns(ip,jp,kp)*delt*s(l)
              q(ip,jp,kp+1,l) = q(ip,jp,kp+1,l)&
     &              -epsns(ip,jp,kp+1)*delt*s(l)
		  END DO



!FIELD CONFINEMENT
	ELSE

!nabla(eta) in cellcenter
            neta_x = (eta(ip+1,jp,kp)-eta(ip,jp,kp)&
     &           +eta(ip+1,jp,kp+1)-eta(ip,jp,kp+1)&
     &           +eta(ip+1,jp+1,kp+1)-eta(ip,jp+1,kp+1)&
     &           +eta(ip+1,jp+1,kp)-eta(ip,jp+1,kp))/(4*delx)
            neta_y = (eta(ip,jp+1,kp)-eta(ip,jp,kp)&
     &           +eta(ip,jp+1,kp+1)-eta(ip,jp,kp+1)&
     &           +eta(ip+1,jp+1,kp+1)-eta(ip+1,jp,kp+1)&
     &           +eta(ip+1,jp+1,kp)-eta(ip+1,jp,kp))/(4*dely)
            neta_z = (eta(ip,jp,kp+1)-eta(ip,jp,kp)&
     &           +eta(ip,jp+1,kp+1)-eta(ip,jp+1,kp)&
     &           +eta(ip+1,jp+1,kp+1)-eta(ip+1,jp+1,kp)&
     &           +eta(ip+1,jp,kp+1)-eta(ip+1,jp,kp))/(4*delz)
            mageta = SQRT(max(0.00001,(neta_x**2+neta_y**2+neta_z**2)))
!correction value
          IF (mageta .EQ. 0.) THEN

            s(1) = 0.
            s(2) = 0.
            s(3) = 0.
          ELSE          
            s(1) = neta_z/mageta*omega(ip,jp,kp,2)&
     &            -neta_y/mageta*omega(ip,jp,kp,3)
            s(2) = neta_x/mageta*omega(ip,jp,kp,3)&
     &            -neta_z/mageta*omega(ip,jp,kp,1)
            s(3) = neta_y/mageta*omega(ip,jp,kp,1)&
     &            -neta_x/mageta*omega(ip,jp,kp,2)
          END IF  
!weighting values
          nenner = MIN(0.,eta(ip+1,jp+1,kp+1)-eta_c(ip,jp,kp))&
     &             +MIN(0.,eta(ip+1,jp+1,kp)-eta_c(ip,jp,kp)) &
     &             +MIN(0.,eta(ip+1,jp,kp)-eta_c(ip,jp,kp))&
     &             +MIN(0.,eta(ip+1,jp,kp+1)-eta_c(ip,jp,kp))&
     &             +MIN(0.,eta(ip,jp+1,kp+1)-eta_c(ip,jp,kp))&
     &             +MIN(0.,eta(ip,jp+1,kp)-eta_c(ip,jp,kp))&
     &             +MIN(0.,eta(ip,jp,kp)-eta_c(ip,jp,kp))&
     &             +MIN(0.,eta(ip,jp,kp+1)-eta_c(ip,jp,kp))
          IF (nenner .EQ. 0.) THEN
            a_1 = 0.125
            a_2 = 0.125
            a_3 = 0.125
            a_4 = 0.125
            a_5 = 0.125
            a_6 = 0.125
            a_7 = 0.125
            a_8 = 0.125
          ELSE
            a_1 = MIN(0.,eta(ip+1,jp+1,kp+1)-eta_c(ip,jp,kp))/nenner
            a_2 = MIN(0.,eta(ip+1,jp+1,kp)-eta_c(ip,jp,kp))/nenner 
            a_3 = MIN(0.,eta(ip+1,jp,kp)-eta_c(ip,jp,kp))/nenner
            a_4 = MIN(0.,eta(ip+1,jp,kp+1)-eta_c(ip,jp,kp))/nenner
            a_5 = MIN(0.,eta(ip,jp+1,kp+1)-eta_c(ip,jp,kp))/nenner
            a_6 = MIN(0.,eta(ip,jp+1,kp)-eta_c(ip,jp,kp)) /nenner
            a_7 = MIN(0.,eta(ip,jp,kp)-eta_c(ip,jp,kp))/nenner
            a_8 = MIN(0.,eta(ip,jp,kp+1)-eta_c(ip,jp,kp))/nenner 
          END IF
!correction step
            DO l = 1,3
              q(ip+1,jp+1,kp+1,l) = q(ip+1,jp+1,kp+1,l)&
     &              +eps(ip+1,jp+1,kp+1)*delt*a_1*s(l)
              q(ip+1,jp+1,kp,l) = q(ip+1,jp+1,kp,l)&
     &              +eps(ip+1,jp+1,kp)*delt*a_2*s(l)
              q(ip+1,jp,kp,l) = q(ip+1,jp,kp,l)&
     &              +eps(ip+1,jp,kp)*delt*a_3*s(l)
              q(ip+1,jp,kp+1,l) = q(ip+1,jp,kp+1,l)&
     &              +eps(ip+1,jp,kp+1)*delt*a_4*s(l)
              q(ip,jp+1,kp+1,l) = q(ip,jp+1,kp+1,l)&
     &              +eps(ip,jp+1,kp+1)*delt*a_5*s(l)
              q(ip,jp+1,kp,l) = q(ip,jp+1,kp,l)&
     &              +eps(ip,jp+1,kp)*delt*a_6*s(l)
              q(ip,jp,kp,l) = q(ip,jp,kp,l)&
     &              +eps(ip,jp,kp)*delt*a_7*s(l)
              q(ip,jp,kp+1,l) = q(ip,jp,kp+1,l)&
     &              +eps(ip,jp,kp+1)*delt*a_8*s(l)
		  END DO

	END IF !field or surface

80        CONTINUE
90      CONTINUE
100   CONTINUE 

!boundary
      CALL boundary !(qinf)

!no slip for pressure calculation
      IF (nswit .EQ. 0) THEN
        IF (n .GT. 15) THEN
          CALL impnoslip
        ELSE   
          DO 150 i = 1,imax
            DO 140 j = 1,jmax
              DO 130 k = 1,kmax
                IF (F(i,j,k) .LE. 0) THEN
                  DO 125 l = 1,3
                    q(i,j,k,l) = 0.
125               CONTINUE
                END IF
130           CONTINUE
140         CONTINUE
150       CONTINUE
        END IF
      ELSE 
	  IF (nswit.NE. 3.) THEN
          CALL impnoslip      
	  END IF  
      END IF       

      RETURN
END SUBROUTINE

END MODULE