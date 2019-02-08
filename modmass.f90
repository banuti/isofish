MODULE modmass



CONTAINS

SUBROUTINE masscon !(phi)
	USE flowprops

!calculation of phi(n+1) so, that q(n+1) satisfies mass conservation 

!      PARAMETER (imax=189,jmax=71,kmax=101,
!     &           ipmax=188,jpmax=70,kpmax=100)  

DIMENSION marker_2(ipmax,jpmax,kpmax),		&		!phi(ipmax,jpmax,kpmax), 
     &          BDXS(jpmax,kpmax), W(2000),&
     &          BDXF(jpmax,kpmax), BDYS(ipmax,kpmax),&
     &          BDYF(ipmax,kpmax), BDZS(ipmax,jpmax),&
     &          BDZF(ipmax,jpmax)

	REAL::LBDCND,MBDCND,NBDCND 

!      COMMON /level/ F(imax,jmax,kmax)
!      COMMON /grid/ delx,dely,delz
!      COMMON /velo/ q(imax,jmax,kmax,3)

      WRITE(*,*) "check: mass conservation"

!-----------------------------------------------------
!right hand side for Poisson equation

      DO 30 ip = 1,ipmax
        DO 20 jp = 1,jpmax 
          DO 10 kp = 1,kpmax
            rem = 0.
            DO 5 idel = 0,1
              DO 4 jdel = 0,1
                DO 3 kdel = 0,1
                  a = F(ip+idel,jp+jdel,kp+kdel)
                  IF (a .GT. 0.) THEN
                    rem = rem + 1
                  END IF
3               CONTINUE
4             CONTINUE
5           CONTINUE
            IF (rem .GT. 0.) THEN
              marker_2(ip,jp,kp) = 1
              phi(ip,jp,kp) = -(q(ip+1,jp,kp,1)-q(ip,jp,kp,1)&
     &             +q(ip+1,jp,kp+1,1)-q(ip,jp,kp+1,1)&
     &             +q(ip+1,jp+1,kp+1,1)-q(ip,jp+1,kp+1,1)&
     &             +q(ip+1,jp+1,kp,1)-q(ip,jp+1,kp,1))/(4*delx)   &
     &             -(q(ip,jp+1,kp,2)-q(ip,jp,kp,2)&
     &             +q(ip,jp+1,kp+1,2)-q(ip,jp,kp+1,2)&
     &             +q(ip+1,jp+1,kp+1,2)-q(ip+1,jp,kp+1,2)&
     &             +q(ip+1,jp+1,kp,2)-q(ip+1,jp,kp,2))/(4*dely)&
     &             -(q(ip,jp,kp+1,3)-q(ip,jp,kp,3)&
     &             +q(ip,jp+1,kp+1,3)-q(ip,jp+1,kp,3)&
     &             +q(ip+1,jp+1,kp+1,3)-q(ip+1,jp+1,kp,3)&
     &             +q(ip+1,jp,kp+1,3)-q(ip+1,jp,kp,3))/(4*delz)
            ELSE 
              phi(ip,jp,kp) = 0.
            END IF
10        CONTINUE
20      CONTINUE
30    CONTINUE

!values of the normal-derivatives of phi at the boundaries equal zero
      BDXS(1:jpmax,1:kpmax) = 0.0
      BDXF(1:jpmax,1:kpmax) = 0.0
      BDYS(1:ipmax,1:kpmax) = 0.0
      BDYF(1:ipmax,1:kpmax) = 0.0
      BDZS(1:ipmax,1:jpmax) = 0.0
      BDZF(1:ipmax,1:jpmax) = 0.0

!call Poisson solver
!    
      xs = -ipmax/2+0.5
      xf = ipmax/2-0.5
      ys = -jpmax/2+0.5
      yf = jpmax/2-0.5
      zs = -kpmax/2+0.5
      zf = kpmax/2-0.5


	!instream:

!	phi(1:2,:,:)=0.

!boundary conditions
!	LBDCND=3. !3
!	MBDCND=3. !3
!	NBDCND=3. !3
	

CALL HW3CRT(xs,xf,ipmax-1,3,BDXS,BDXF,ys,yf,jpmax-1,3,BDYS,BDYF,zs,zf,kpmax-1,3,BDZS,BDZF,0.,ipmax,jpmax,phi,PERTRB,IERROR,W)

      WRITE(*,*) "finished poisson solver "

!check, if input parameters are valid
      WRITE(*,*) " ierror =" ,ierror
!check, if solution for phi is reasonable
      WRITE(*,*) "pertrb =" ,pertrb

!symmetrisieren
!        DO 3000 jp = 1,jpmax/2
!          DO 2000 ip = 1,ipmax
!            DO 1000 kp = 1,kpmax
!              a = (phi(ip,jp,kp)+phi(ip,jpmax-(jp-1),kp))/2.
!              phi(ip,jp,kp) = a
!              phi(ip,jpmax-(jp-1),kp) = phi(ip,jp,kp)
!1000        CONTINUE
!2000      CONTINUE
!3000    CONTINUE

      WRITE(*,*) "Konvergenz-Kontrolle:"
!     WRITE(*,*) phi(80,35,50)
!	WRITE(*,*) phi(81,35,50)
!      WRITE(*,*) phi(82,35,50)
!      WRITE(*,*) phi(83,35,50)
      WRITE(*,*) phi(85,j0,50)

!	WRITE(*,*) phi(2,2,2)
!	WRITE(*,*) phi(1,1,1)
      RETURN
END SUBROUTINE

!***********************************************************************

SUBROUTINE momentum !(phi,qinf)
	USE flowprops
	USE boundaries

!satisfy momentum conservation

!      PARAMETER (imax=189,jmax=71,kmax=101,
!     &           ipmax=188,jpmax=70,kpmax=100) 

!      DIMENSION phi(ipmax,jpmax,kpmax), qinf(3)

!      COMMON /grid/ delx,dely,delz
 !     COMMON /velo/ q(imax,jmax,kmax,3)

      WRITE(*,*) "check: momentum"
!
!--------------------------------------------------------
!1. inner domain

      DO 30 i = 2,imax-1
        DO 20 j = 2,jmax-1
          DO 10 k = 2,kmax-1
            q(i,j,k,1) = q(i,j,k,1)+(phi(i,j-1,k-1)-phi(i-1,j-1,k-1)&
     &                   +phi(i,j-1,k)-phi(i-1,j-1,k)&
     &                   +phi(i,j,k)-phi(i-1,j,k)&
     &                   +phi(i,j,k-1)-phi(i-1,j,k-1))/(4.*delx)
            q(i,j,k,2) = q(i,j,k,2)+(phi(i-1,j,k-1)-phi(i-1,j-1,k-1)&
     &                   +phi(i-1,j,k)-phi(i-1,j-1,k)&
     &                   +phi(i,j,k)-phi(i,j-1,k)&
     &                   +phi(i,j,k-1)-phi(i,j-1,k-1))/(4.*dely)
            q(i,j,k,3) = q(i,j,k,3)+(phi(i-1,j-1,k)-phi(i-1,j-1,k-1)&
     &                   +phi(i-1,j,k)-phi(i-1,j,k-1)&
     &                   +phi(i,j,k)-phi(i,j,k-1)&
     &                   +phi(i,j-1,k)-phi(i,j-1,k-1))/(4.*delz)
		
10        CONTINUE
20      CONTINUE
30    CONTINUE

!2. boundary
	
      CALL boundary !(qinf)

!no slip
      CALL impnoslip
!  
      RETURN
      
END SUBROUTINE


END MODULE