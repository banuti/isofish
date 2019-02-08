MODULE flowprops

IMPLICIT NONE

!flowfield
	INTEGER	::imax,jmax,kmax,ipmax,jpmax,kpmax,i0,j0,k0
	INTEGER	::nmax,n_step,outoffs,outint,nswit,perturb,geometry
	REAL	::delt,delx,dely,delz,xmue,rho,alpha_deg,alpha,qinfinity,time,D,fieldRad


!	PARAMETER (imax=91,jmax=7,kmax=51,ipmax=90,jpmax=6,kpmax=50)		!D=10
	PARAMETER (imax=181,jmax=7,kmax=101,ipmax=180,jpmax=6,kpmax=100)	!D=20
! PARAMETER (imax=181,jmax=4,kmax=101,ipmax=180,jpmax=3,kpmax=100)	!D=20
	! PARAMETER (imax=361,jmax=7,kmax=201,ipmax=360,jpmax=6,kpmax=200)	!D=40


 !   PARAMETER (imax=189,jmax=71,kmax=101,ipmax=188,jpmax=70,kpmax=100)


!	PARAMETER (imax=201,jmax=81,kmax=161,ipmax=200,jpmax=80,kpmax=160)

  REAL	::F(imax,jmax,kmax), Fnormp(ipmax,jpmax,kpmax,3), Fp(ipmax,jpmax,kpmax)
	REAL	::qinf(3),q(imax,jmax,kmax,3),q_avg(imax,jmax,kmax,3)
	REAL	::phi(ipmax,jpmax,kpmax),cp(ipmax,jpmax,kpmax),pi
	REAL	::eps(imax,jmax,kmax), epsns(imax,jmax,kmax),mark(imax,jmax,kmax)
	REAL	::eps_s, epsns_s
	REAL	::omega(ipmax,jpmax,kpmax,3),absomega(ipmax,jpmax,kpmax)
	REAL	::Cd,Cl,spd(imax,jmax,kmax)

	REAL, DIMENSION(:,:), ALLOCATABLE::cppnts

	INTEGER	::numpnts,numvar

CONTAINS

!===================================================================

SUBROUTINE init

INTEGER	::islice,jslice,kslice,i,j,k,ip,jp,kp,iprime,kprime,jprime
REAL	::Fabs,Fptemp,h

! GRID
h = 1.

delx = h
dely = h
delz = h


jslice=j0
numvar=4

q_avg=0.
Cd=0.
Cl=0.

eps=eps_s
epsns=epsns_s


! switch for continued calculation
! 0 = new calculation
! 1 = continued calculation
! 3 = postprocessing


! Zeitschritt
!      nmax = 4000
!      delt = 0.2
! eps,epsns,mue
!      eps_s = 0.	!0.065
!      epsns_s = 0.	!0.2
!      xmue = 0.1	!0.0025
! Anstroemgeschwindigkeit
      qinfinity = 1.
! angle of attack
!      alpha_deg = 0.
! in Bogenmass
      pi = 3.14159265358979
      alpha = pi/180.*alpha_deg
!
      qinf(1) = COS(alpha)*qinfinity
      qinf(2) = 0.
      qinf(3) = SIN(alpha)*qinfinity


!TRANSFORM LEVELSET FUNCTION TO PRESSURE GRID
WRITE(*,*)'Transform Grid'
	DO ip=1, ipmax
	  DO jp=1, jpmax
	    DO kp=1, kpmax

			Fptemp=0.

			DO i=ip,ip+1
			  DO j=jp, jp+1
				DO k=kp, kp+1
					Fptemp=Fptemp+F(i,j,k)
				END DO
			  END DO
			END DO

			Fp(ip,jp,kp)=Fptemp/8.

		END DO
	  END DO
	END DO



!GET GRADIENT OF LEVELSET IN PRESSURE GRID

	Fnormp=0.

	DO i=2, ipmax-1
	  DO j=2, jpmax-1
	    DO k=2, kpmax-1

		Fnormp(i,j,k,1)=(Fp(i+1,j,k)-Fp(i-1,j,k))/(2*delx)
		Fnormp(i,j,k,2)=(Fp(i,j+1,k)-Fp(i,j-1,k))/(2*dely)
		Fnormp(i,j,k,3)=(Fp(i,j,k+1)-Fp(i,j,k-1))/(2*delz)

		Fabs=SQRT(max(0.00001,(Fnormp(i,j,k,1)**2+Fnormp(i,j,k,2)**2+Fnormp(i,j,k,3)**2)))

		Fnormp(i,j,k,:)=Fnormp(i,j,k,:)/Fabs

		END DO
	  END DO
	END DO

	Fnormp(1,:,:,:)=Fnormp(2,:,:,:)
	Fnormp(ipmax,:,:,:)=Fnormp(ipmax-1,:,:,:)

	Fnormp(:,1,:,:)=Fnormp(:,2,:,:)
	Fnormp(:,jpmax,:,:)=Fnormp(:,jpmax-1,:,:)

	Fnormp(:,:,1,:)=Fnormp(:,:,2,:)
	Fnormp(:,:,kpmax,:)=Fnormp(:,:,kpmax-1,:)





!ESTIMATE NUMBER OF BOUNDARY CELLS
	numpnts=0

	DO k=1,kpmax
	  DO i=1,ipmax
		kprime=k
		iprime=i

		IF (Fp(i,jslice,k).LT. 0.) THEN
		  DO iprime=(i-1), (i+1), 2
			IF (Fp(iprime,jslice,k).GT. 0.) numpnts=numpnts+1
		  END DO
		  DO kprime=(k-1), (k+1), 2
			IF (Fp(i,jslice,kprime).GT. 0.) numpnts=numpnts+1
		  END DO
		END IF
	  END DO
	END DO

	WRITE(*,*)'Number of boundarycells: ',numpnts


	ALLOCATE(cppnts(numpnts, numvar))







END SUBROUTINE

!------------------------------------------------------------------

SUBROUTINE readprops

CHARACTER(LEN=10)	::dummy

WRITE(*,*)'reading properties from file...'

OPEN(10,file='properties.inp')
READ(10,*)dummy,nswit,dummy,eps_s,dummy,epsns_s,dummy,xmue,dummy, &
					delt,dummy,nmax,dummy,outoffs,dummy,outint,dummy,alpha_deg, &
					dummy,perturb,dummy,D,dummy,fieldRad,dummy,geometry
CLOSE(10)

WRITE(*,*)'...done:'
WRITE(*,*)'nswit	   = ',nswit 			!computation mode
WRITE(*,*)'nmax		   = ',nmax			!timesteps
WRITE(*,*)'delt		   = ',delt
WRITE(*,*)'eps_s	   = ',eps_s
WRITE(*,*)'epsns_s	 = ',epsns_s
WRITE(*,*)'xmue		   = ',xmue
WRITE(*,*)'alpha_deg = ',alpha_deg		!angle of attack
WRITE(*,*)'outoffs	 = ',outoffs			!offset start of output
WRITE(*,*)'outint	   = ',outint			!output interval
WRITE(*,*)'perturb	 = ',perturb			!perturbation?
WRITE(*,*)'fieldRad  = ',fieldRad		!radius for field confinement
WRITE(*,*)'geometry  = ',geometry		!levelset fct
WRITE(*,*)


END SUBROUTINE

!------------------------------------------------------------------


SUBROUTINE initial !(qinf,nswit,time)


REAL	::R,deli,delj,delk,rad4,speed
INTEGER	::i,j,k,l,ix,jy,kz

! initial conditions for far field flow around an ellipsoid

!      PARAMETER (imax=189,jmax=71,kmax=101)
!	REAL::F

!      DIMENSION qinf(3)

!      COMMON /velo/ q(imax,jmax,kmax,3)

      WRITE(*,*) "check: initial values"
	  WRITE(*,*)'comp mode nswit: ',nswit



R=D/2.



!------------------------------------------------------------------

!POTENTIAL FLOWFIELD
IF (nswit .EQ. 0) THEN

!	WRITE(*,*)'reading initial potentialflow data file'
!	READ(99) ix,jy,kz,l
!	READ(99) (((q(i,j,k,1),i=1,ix),j=1,jy),k=1,kz),(((q(i,j,k,2),i=1,ix),j=1,jy),k=1,kz),&
!   &            (((q(i,j,k,3),i=1,ix),j=1,jy),k=1,kz)
!	CLOSE(99)

WRITE(*,*)'computing initial potential flowfield'



!CIRCULAR CYLINDER
IF (geometry .EQ. 1) THEN


DO i=1, imax
  DO j=1, jmax
	DO k=1, kmax

	IF (F(i,j,k) .GT. 0.) THEN

		deli=real(i-i0)
		delj=real(j-j0)
		delk=real(k-k0)

		q(i,j,k,1)=qinf(1)*(1+(R**2*(delk**2-deli**2)/(MAX(0.01,(deli**2+delk**2)**2))))
		q(i,j,k,2)=0.
		q(i,j,k,3)=-2*qinf(1)*R**2*deli*delk/(MAX(0.01,(deli**2+delk**2)**2))

	ELSE

		q(i,j,k,:)=0.

	END IF

	END DO
  END DO
END DO



!	OPEN(unit=56,file='potential.dat')
!	WRITE(56,*) 'Title= "potential"'
!	WRITE(56,*) 'Variables=	"i","j","U","W","spd"'
!	WRITE(*,*)'Writing pot File...'
!	WRITE(56,*)'Zone T="n=1",I=',imax,',K=',kmax,',F= Point'

!	DO k=1, kmax
!	    DO i=1, imax

!	speed=SQRT(MAX(0.00001,(q(i,35,k,1)**2+q(i,35,k,2)**2)))
!	WRITE(56,*) i,k,q(i,35,k,1),q(i,35,k,3),speed

!	END DO
!	END DO


END IF


END IF




!CONTINOUS CALCULATION AND POSTPROCESSING
IF ((nswit .EQ. 1) .OR. (nswit .EQ. 3.) )THEN

	WRITE(*,*)'reading previous data file'
	READ(100) time
	CLOSE(100)
	READ(101) ix,jy,kz,l
	READ(101) (((q(i,j,k,1),i=1,ix),j=1,jy),k=1,kz),(((q(i,j,k,2),i=1,ix),j=1,jy),k=1,kz),&
				& (((q(i,j,k,3),i=1,ix),j=1,jy),k=1,kz)
	CLOSE(101)


    time = 0.


END IF




!	q(:,:,:,1) = qinf(1)
!	q(:,:,:,2) = qinf(2)
!	q(:,:,:,3) = qinf(3)





!	  DO i = 1,imax
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
!				q(i,j,k,1) = qinf(1)
!				q(i,j,k,2) = qinf(2)
!				q(i,j,k,3) = qinf(3)
!			END IF
!	      END DO
!	    END DO
!	  END DO




      RETURN
      END SUBROUTINE

!------------------------------------------------------------------

SUBROUTINE perturbation

INTEGER	::i,j,k

DO i=i0, i0+D/2.
  DO k=k0+5, k0+D/1.5

	IF (F(i,j0,k) .GE. 1.0) THEN

		q(i,:,k,1)=1.7
		q(i,:,k,2)=0.
		q(i,:,k,3)=0.

	END IF

  END DO
END DO




END SUBROUTINE

!------------------------------------------------------------------------
SUBROUTINE level_funct


REAL	::deli,delj,delk,R
INTEGER	::i,j,k,l,ix,jy,kz

R=real(D)/2.

WRITE(*,*)'R=',R


!READ LEVELSET FROM FILE
IF (geometry .EQ. 0) THEN

WRITE(*,*)'reading levelset file'

      READ(22) ix,jy,kz,l
      READ(22) (((F(i,j,k),i=1,ix),j=1,jy),k=1,kz)
      CLOSE(22)


	F=F*10.

END IF


!CALCULATION


!X-CYL
IF (geometry .EQ. 1) THEN

WRITE(*,*)'computing x-cyl levelset'

	! i0=3.5*D
	i0 = imax/4. + R
	j0 = jmax/2. - 0.5
	k0 = kmax/2. - 0.5

	WRITE(*,*)'i0=',i0
	WRITE(*,*)'j0=',j0
	WRITE(*,*)'k0=',k0

	 DO i = 1,imax
        DO j = 1,jmax
          DO k = 1,kmax

		  deli=real(i-i0)
		  delj=real(j-j0)
		  delk=real(k-k0)

  F(i,j,k) = (SQRT(deli**2. + delk**2.) - R)

        END DO
      END DO
    END DO

END IF

CONTINUE

RETURN
END SUBROUTINE


!==================================================================
END MODULE
