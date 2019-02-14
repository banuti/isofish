!***********************************************************************
!
!     PROGRAM isofish
!
!***********************************************************************
!
!   Testbed for modeling
!		daniel@banuti.com
!
!
!TODO
! - call boundaries from main loop
! - refactor main loop
! - refactor
! - CFL
! - periodic shear layer
! - passive scalar
! - output pressure
! - 2D minimum width?
! - 1st/2nd order switch
! - modulo screen report
! - write restart file
! - read restart file
! - calculate spectrum
! - run some DNS case?
! - make smoother no-slip boundary
! - square cylinder?
!----------------------------------------------------------------------


USE flowprops
USE modoutput
USE modanalysis
USE boundaries
USE modconvect
USE moddiffuse
USE modmass
USE modconfinement

INTEGER::saved
REAL::qmax


!=================================================================================
WRITE(*,*)
WRITE(*,*)'================================'
WRITE(*,*)'isofish flow simulation   __||_    '
WRITE(*,*)'Daniel Banuti, 2019       \_||_\  '
WRITE(*,*)'daniel@banuti.com          2019  '
WRITE(*,*)'================================'
WRITE(*,*)


!level-set function

!CALL readinput !in flowprops.f90
!CALL initlevelset !in flowprops.f90
!CALL initfield !in flowprops.f90
!CALL initqfield !(qinf,nswit,time)

CALL init


! initial conditions


saved=0

!----------------------------------------------------------------------
CALL outhead

!
!IF (nswit .EQ. 0) THEN
!  IF (perturb .NE. 0) THEN
!
!    CALL perturbation
!
!  END IF
!END IF
!
!



!calculation for every timestep

!*******************LOOP**********************
DO 1000 n = 1,nmax


  time = time+delt
  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) 'T I M E S T E P: ',n,time

!  IF (n .EQ. nmax) THEN
!          deltat = delt
!   END IF





  IF (nswit .NE. 3) THEN

    !convection calculation
    CALL convection

    !add diffusion
    CALL diffusion 

  END IF

!
!add confinement term
    CALL vort_conf(n) !(qinf,delt,n)

!mass conservation
    CALL masscon !(phi)


  IF (nswit.NE. 3.) THEN
    !satisfy momentum conservation
    CALL momentum !(phi,qinf)
  END IF


CALL getspeed

!Velocity-Field
qmax=maxval(q)
WRITE(*,*)'Maximum Velocity-Component: ', qmax
WRITE(*,*)'at: ',maxloc(q)
WRITE(*,*)'Freestream: q(5,j0,k0,1): ',q(5,j0,k0,1)
WRITE(*,*)'Freestream: phi(5,j0,k0): ',phi(5,j0,k0)




!in case of instability
IF ((qmax .GT. 5.) .AND. saved .eq. 0) THEN
  CALL pressure !(deltat) !(phi,cp,deltat,qinfinity)
  CALL output !(time, cp)
  saved=1
END IF



!REPORT
      CALL pressure
CALL drag
CALL lift
CALL outreports(n,qmax)


WRITE(*,*)'Drag: ',Cd
WRITE(*,*)'Lift: ',Cl
WRITE(*,*)'Max vort: ',maxval(absomega)



!from timestep outoffs save every outint'th step
IF ( (n .GE. outoffs) .AND. ( (modulo(n-outoffs,outint) .EQ. 0) &
    .OR. (n-outoffs) .EQ. 0)) THEN

  WRITE(*,*)'Output timestep ',n
  CALL outstep(n,qmax)

END IF



!SUMMATION TO GET AVERAGE VALUES
IF (n .GE. outoffs) THEN
  q_avg=q_avg+q
END IF



1000  CONTINUE
!**********************END LOOP*********************************

!COMPUTE AVG
q_avg=q_avg/real(n-outoffs+1)



CALL pressure !(deltat) !(phi,cp,deltat,qinfinity)

CALL outavg


CALL output !(time,cp)

CALL outend

DEALLOCATE(cppnts)


      END
