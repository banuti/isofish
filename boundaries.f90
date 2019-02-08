MODULE boundaries


!TODO
! - assign BC to face
! - choose BC from parafile
! - types:
! -- extrapolation
! -- dirichlet
! -- periodic 


USE flowprops



CONTAINS

!================================================================
SUBROUTINE boundary_num(qq) !(imax,jmax,kmax,q,qinf)

USE flowprops

! INTEGER	::d



! boundary conditions


DIMENSION qq(imax,jmax,kmax,3),q_ik_slice(imax,kmax,3) !qinf(3),

! boundary conditions


! influx
	qq(1:2,:,:,1) = qinf(1)
	qq(1:2,:,:,2) = qinf(2)
	qq(1:2,:,:,3) = qinf(3)





!NOZZLE

!	i = 50
!          DO j = 1,jmax
!            DO k = 1,kmax
!			F = i**2./60**2.
!     &               + j**2./10**2.
!     &               + k**2./10**2. - 1.
!			IF (F.LE.0.0) THEN
!				q(i,j,k,1) = 2.0
!				q(i,j,k,2) = 0.0
!				q(i,j,k,3) = 0.0
!			END IF
!	      END DO
!	    END DO


!OLD

!        i = 1
!        DO 20 j = 1,jmax
!          DO 10 k = 1,kmax
!            q(i,j,k,1) = qinf(1)
!           q(i,j,k,2) = qinf(2)
!            q(i,j,k,3) = qinf(3)
!10        CONTINUE
!20      CONTINUE
!------------------
! BOTTOM
	  kdist_1 = 1
      kdist_2 = 2

      qq(:,:,1,:) = 2*qq(:,:,1+kdist_1,:)-qq(:,:,1+kdist_2,:)



!        k = 1
!        DO 22   i = 2,imax-1
!          DO 11 j = 2,jmax-1
!            q(i,j,k,1) = qinf(1)
!            q(i,j,k,2) = qinf(2)
!            q(i,j,k,3) = qinf(3)
! 11       CONTINUE
! 22     CONTINUE



! TOP
        k = kmax
        kdist_1 = -1
        kdist_2 = -2

        DO i = 2,imax-1
          DO j = 2,jmax-1
            DO l = 1,3
              qq(i,j,k,l) = 2*qq(i,j,k+kdist_1,l)-qq(i,j,k+kdist_2,l)
          END DO
        END DO
      END DO



!------------------
! left,right boundary
        DO 100 j = 1,jmax,jmax-1
          IF (j .EQ. 1) THEN
! right
            jdist_1 = 1
            jdist_2 = 2
          ELSE
! left
            jdist_1 = -1
            jdist_2 = -2
          END IF
!
          DO 90 i = 2,imax-1
            DO 80 k = 1,kmax
              DO 70 l = 1,3
                qq(i,j,k,l) = 2*qq(i,j+jdist_1,k,l)-qq(i,j+jdist_2,k,l)
70            CONTINUE
80          CONTINUE
90        CONTINUE
100     CONTINUE

!	q(:,1:2,:,2)=0.
!	q(:,kmax-1:kmax,:,2)=0.


! LEFT von Neumann
!	q(:,1,:,:)=q(:,1+1,:,:)

	qq(:,1,:,1)=qq(:,1+1,:,1)
	qq(:,1,:,2)=0.
	qq(:,1,:,3)=qq(:,1+1,:,3)


!RIGHT von Neumann
!	q(:,jmax,:,:)=q(:,jmax-1,:,:)

	qq(:,jmax,:,1)=qq(:,jmax-1,:,1)
	qq(:,jmax,:,2)=0.
	qq(:,jmax,:,3)=qq(:,jmax-1,:,3)


!------------------
! outflux

!	i = imax
!          DO j = 1,jmax
!            DO k = 1,kmax
!			F = i**2./60**2.
!     &               + j**2./10**2.
!     &               + k**2./10**2. - 1.
!			IF (F.LE.0.0) THEN
!				q(i,j,k,1) = 0.0
!				q(i,j,k,2) = 0.0
!				q(i,j,k,3) = 0.0
!			ELSE
!				q(i,j,k,1) = q(i-1,j,k,1)
!				q(i,j,k,2) = q(i-1,j,k,2)
!				q(i,j,k,3) = q(i-1,j,k,3)
!			END IF
!	      END DO
!	    END DO

        i = imax
        idist_1 = -1
        idist_2 = -2
        DO 130 j = 1,jmax
          DO 120 k = 1,kmax
            DO 110 l = 1,3
              qq(i,j,k,l) = 2.*qq(i+idist_1,j,k,l)-qq(i+idist_2,j,k,l)
110         CONTINUE
120       CONTINUE
130     CONTINUE



!TOP/BOTTOM BC: W=0.
! qq(:,:,1,3)=0.
! qq(:,:,kmax,3)=0.


!TOP/BOTTOM BC: q=qinf

!DO d=1, 3
!	qq(:,:,1,d)=qinf(d)
!	qq(:,:,kmax,d)=qinf(d)
!END DO

!equalize to q2d

qq(:,:,:,2)=0.


q_ik_slice(:,:,:)=qq(:,j0,:,:)


DO j=1,jmax

qq(:,j,:,:)=q_ik_slice(:,:,:)

END DO




RETURN
END SUBROUTINE


!-------------------------------------------------------------------------------
SUBROUTINE boundary !(qinf)
	USE flowprops


! INTEGER	::d

! boundary conditions

!      PARAMETER(imax=189, jmax=71, kmax=101)

     DIMENSION q_ik_slice(imax,kmax,3)

 !     COMMON /velo/ q(imax,jmax,kmax,3)
!------------------
! influx
	q(1:2,:,:,1) = qinf(1)
	q(1:2,:,:,2) = qinf(2)
	q(1:2,:,:,3) = qinf(3)





!NOZZLE

!	i = 50
!          DO j = 1,jmax
!            DO k = 1,kmax
!			F = i**2./60**2.
!     &               + j**2./10**2.
!     &               + k**2./10**2. - 1.
!			IF (F.LE.0.0) THEN
!				q(i,j,k,1) = 2.0
!				q(i,j,k,2) = 0.0
!				q(i,j,k,3) = 0.0
!			END IF
!	      END DO
!	    END DO


!OLD

!        i = 1
!        DO 20 j = 1,jmax
!          DO 10 k = 1,kmax
!            q(i,j,k,1) = qinf(1)
!           q(i,j,k,2) = qinf(2)
!            q(i,j,k,3) = qinf(3)
!10        CONTINUE
!20      CONTINUE
!------------------
! BOTTOM
	 kdist_1 = 1
       kdist_2 = 2

        q(:,:,1,:) = 2*q(:,:,1+kdist_1,:)-q(:,:,1+kdist_2,:)



!        k = 1
!        DO 22   i = 2,imax-1
!          DO 11 j = 2,jmax-1
!            q(i,j,k,1) = qinf(1)
!            q(i,j,k,2) = qinf(2)
!            q(i,j,k,3) = qinf(3)
! 11       CONTINUE
! 22     CONTINUE



! TOP
        k = kmax
        kdist_1 = -1
        kdist_2 = -2
!
        DO i = 2,imax-1
          DO j = 2,jmax-1
            DO l = 1,3
              q(i,j,k,l) = 2*q(i,j,k+kdist_1,l)-q(i,j,k+kdist_2,l)
          END DO
        END DO
      END DO
!------------------
! left,right boundary
        DO 100 j = 1,jmax,jmax-1
          IF (j .EQ. 1) THEN
! right
            jdist_1 = 1
            jdist_2 = 2
          ELSE
! left
            jdist_1 = -1
            jdist_2 = -2
          END IF
!
          DO 90 i = 2,imax-1
            DO 80 k = 1,kmax
              DO 70 l = 1,3
                q(i,j,k,l) = 2*q(i,j+jdist_1,k,l)-q(i,j+jdist_2,k,l)
70            CONTINUE
80          CONTINUE
90        CONTINUE
100     CONTINUE

!	q(:,1:2,:,2)=0.
!	q(:,kmax-1:kmax,:,2)=0.

! LEFT von Neumann
!	q(:,1,:,:)=q(:,1+1,:,:)

	q(:,1,:,1)=q(:,1+1,:,1)
	q(:,1,:,2)=0.
	q(:,1,:,3)=q(:,1+1,:,3)


!RIGHT von Neumann
!	q(:,jmax,:,:)=q(:,jmax-1,:,:)

	q(:,jmax,:,1)=q(:,jmax-1,:,1)
	q(:,jmax,:,2)=0.
	q(:,jmax,:,3)=q(:,jmax-1,:,3)


!------------------
! outflux

!	i = imax
!          DO j = 1,jmax
!            DO k = 1,kmax
!			F = i**2./60**2.
!     &               + j**2./10**2.
!     &               + k**2./10**2. - 1.
!			IF (F.LE.0.0) THEN
!				q(i,j,k,1) = 0.0
!				q(i,j,k,2) = 0.0
!				q(i,j,k,3) = 0.0
!			ELSE
!				q(i,j,k,1) = q(i-1,j,k,1)
!				q(i,j,k,2) = q(i-1,j,k,2)
!				q(i,j,k,3) = q(i-1,j,k,3)
!			END IF
!	      END DO
!	    END DO

        i = imax
        idist_1 = -1
        idist_2 = -2
        DO 130 j = 1,jmax
          DO 120 k = 1,kmax
            DO 110 l = 1,3
              q(i,j,k,l) = 2.*q(i+idist_1,j,k,l)-q(i+idist_2,j,k,l)
110         CONTINUE
120       CONTINUE
130     CONTINUE




!TOP/BOTTOM BC: W=0.
! q(:,:,1,3)=0.
! q(:,:,kmax,3)=0.


!TOP/BOTTOM BC: q=qinf
!DO d=1,3
!	q(:,:,1,d)=qinf(d)
!	q(:,:,kmax,d)=qinf(d)
!END DO



!equalize to q2d
q(:,:,:,2)=0.

q_ik_slice(:,:,:)=q(:,j0,:,:)


DO j=1,jmax

q(:,j,:,:)=q_ik_slice(:,:,:)

END DO



      RETURN
END SUBROUTINE


!------------------------------------------------------------------------
SUBROUTINE s_extrapolation(s)
	USE flowprops

! boundary treatment with extrapolation for scalar fields on
! velocity grid

      DIMENSION s(imax,jmax,kmax)

!------------------
! inflow,outflow
      DO 30 i = 1,imax,imax-1
        IF (i .EQ. 1) THEN
! inflow
          idist_1 = 1
          idist_2 = 2
        ELSE
! outflow
          idist_1 = -1
          idist_2 = -2
        END IF
        DO 20 j = 2,jmax-1
          DO 10 k = 2,kmax-1
            s(i,j,k) = 2*s(i+idist_1,j,k)-s(i+idist_2,j,k)
10        CONTINUE
20      CONTINUE
30    CONTINUE

!override inflow:

s(1:2,:,:)=0.




!------------------
! left,right boundary
      DO 60 j = 1,jmax,jmax-1
        IF (j .EQ. 1) THEN
! right
          jdist_1 = 1
          jdist_2 = 2
        ELSE
! left
          jdist_1 = -1
          jdist_2 = -2
        END IF
!
        DO 50 i = 1,imax
          DO 40 k = 2,kmax-1
            s(i,j,k) = 2*s(i,j+jdist_1,k)-s(i,j+jdist_2,k)
40        CONTINUE
50      CONTINUE
60    CONTINUE



!------------------
! top,bottom
      DO 90 k = 1,kmax,kmax-1
        IF (k .EQ. 1) THEN
! bottom
          kdist_1 = 1
          kdist_2 = 2
        ELSE
! top
          kdist_1 = -1
          kdist_2 = -2
        END IF
!
        DO 80 i = 1,imax
          DO 70 j = 1,jmax
            s(i,j,k) = 2*s(i,j,k+kdist_1)-s(i,j,k+kdist_2)
70        CONTINUE
80      CONTINUE
90    CONTINUE
!---------------------------
      RETURN
END SUBROUTINE


!--------------------------------------------------------------------------

SUBROUTINE impnoslip
	USE flowprops

! improved enforcement of no slip condition

      REAL nenner

!      PARAMETER (imax =189,jmax =71,kmax =101,
!     &           ipmax=188,jpmax=70,kpmax=100)

      DIMENSION q_i(3)

! values are especially for 6:1 ellipsoid with no angle of attack to grid
      i_min = 1		!(imax+1)/2-71
      i_max = imax	!(imax+1)/2+71
      j_min =	1		!(jmax+1)/2-12
      j_max = jmax	!(jmax+1)/2+12
      k_min = 1		!(kmax+1)/2-12
      k_max = kmax	!(kmax+1)/2+12

      dist_max = 2.

!-----------------------------------------------------------------------
      WRITE(*,*) ' Anfang:impnoslip'

      DO i = i_min,i_max
        DO j = j_min,j_max
          DO k = k_min,k_max



            IF (F(i,j,k) .LE. 0.) THEN

			 q(i,j,k,:) = 0.

            END IF

		END DO
	  END DO
	END DO


		!grid point in boundary layer
	!	c_noslip=2.
	!	IF ((F(i,j,k) .GE. 0.) .AND. (F(i,j,k) .LT. 40.0)) THEN
	!		DO dir=1, 3
	!		 q(i,j,k,dir)=q(i,j,k,dir)*F(i,j,k)/50.0
	!		END DO
	!	END IF

	RETURN
 END SUBROUTINE
!====================================================================


END MODULE
