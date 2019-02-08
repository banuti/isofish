MODULE modoutput
USE flowprops



CONTAINS
!=====================================================================================
SUBROUTINE outhead


		OPEN(unit=56,file='movie.dat')
		WRITE(56,*) 'Title= "Movie"'
		WRITE(56,*) 'Variables=	"i","k","U","W","spd","vort","F"'
		CALL outinfoheader(56)



	  OPEN(unit=60,file='fielddata.dat')
		WRITE(60,*) 't i k U W spd vort F'


		OPEN(unit=57,file='report.dat')
		WRITE(57,*) 'Title= "Report n,qmax,Cd,Cl"'
		WRITE(57,*) 'Variables="n","qmax","qfree","Cd","Cl","Cpb"'
		CALL outinfoheader(57)


		OPEN(unit=58,file='stream_q.dat')
		WRITE(58,*) 'Title= "Streamwise q"'
		WRITE(58,*) 'Variables="i","U"'	!"time",
		CALL outinfoheader(58)


		OPEN(unit=59,file='cross_q.dat')
		WRITE(59,*) 'Title= "Cross q"'
		WRITE(59,*) 'Variables="k","U 0.0","U 0.6","U 1.0","U 2.0"'	!"time",
		CALL outinfoheader(59)


		OPEN(unit=205, file='xcyl_cp-slice.dat')
		WRITE(205,*) 'Title= "Cp-slice"'
		WRITE(205,*) 'Variables="theta","cp_bnd"'
		CALL outinfoheader(205)




	100 FORMAT(' T="\\n eps = ',F5.3,' \\n mu = ',F5.3,' \\n dt = ',F5.3,' \\n"')




END SUBROUTINE

!-------------------------------------------------------------------------------------

SUBROUTINE outreports(n,qmax)

REAL	::qmax,Cpb,theta,speed
INTEGER	::i,j,k,islice,jslice,kslice,baseindex

islice=i0
jslice=j0
kslice=k0

baseindex=numpnts/2

Cpb=(cppnts(baseindex,2)+cppnts(baseindex+1,2))/2.
theta=(cppnts(baseindex,1)+cppnts(baseindex+1,1))/2.

WRITE(*,*) 'Cpb=',Cpb,' at ',theta

!REPORT-FILE
WRITE(57,*) n, maxval(spd), q(20,j0,50,1),Cd, Cl, Cpb
!10 FORMAT(I5, F12.8, F12.8, F12.8, F12.8)



!MOVIE VELOCITY

IF ((n .GE. outoffs) .AND. ((modulo(n,outint) .EQ. 0) .OR. (n .EQ. 1))) THEN

	WRITE(*,*)'Writing Movie File...'
	WRITE(56,*)'Zone T="n=',n,'",I=',imax-1,',K=',kmax-1,',F= Point'

	DO k=1, kmax-1
	    DO i=1, imax-1

	speed=SQRT(MAX(0.00001,(q(i,jslice,k,1)**2+q(i,jslice,k,2)**2+q(i,jslice,k,3)**2)))
	WRITE(56,101) i,k,q(i,jslice,k,1),q(i,jslice,k,3),speed,absomega(i,jslice,k),F(i,jslice,k)
	101 FORMAT(I4, I4, F7.3, F7.3, F7.3, F8.4, F9.3)

	WRITE(60,111) n,i,k,q(i,jslice,k,1),q(i,jslice,k,3),speed,absomega(i,jslice,k),F(i,jslice,k)
  111 FORMAT(I4, I4, I4, F16.12, F16.12, F16.12, F16.12, F16.12)


		END DO
	END DO

	WRITE(*,*)'Finished'

END IF




END SUBROUTINE

!-------------------------------------------------------------------------------------

SUBROUTINE outstep(n,qmax)

REAL	::qmax!,R
INTEGER	::n,i,j,k,islice,jslice,kslice,jcenter

R=D/2.

islice=i0
jslice=j0
kslice=k0
icenter=i0




!WAKE STREAMWISE
	WRITE(58,*)'Zone T="n=',n,'",I=',int((imax-i0-R)/2+1),' ,F= Point'

	DO i=i0+R, imax, 2

		WRITE(58,*) (i-i0)/D,q(i,jslice,k0,1)

	END DO


!WAKE CROSSWISE
	WRITE(59,*)'Zone T="n=',n,'",K=',int(3.0*D+1),',F= Point'

	DO k=k0-1.5*D, k0+1.5*D

		WRITE(59,*) (k-k0)/D,q(i0,jslice,k,1),q(i0+0.6*D,jslice,k,1),q(i0+D,jslice,k,1),q(i0+2*D,jslice,k,1)

	END DO


!SURFACE PRESSURE SLICE
	WRITE(205,*)'Zone T="n=',n,'",K=',numpnts,' ,F= Point'

	DO i=1,numpnts
		WRITE(205,*) cppnts(i,1), cppnts(i,2), cppnts(i,3)
	END DO

	!WRITE(205,*) ((cppnts(i,j), j=1,numvar), i=1,numpnts)


END SUBROUTINE

!-------------------------------------------------------------------------------------

SUBROUTINE outend

	CLOSE(56)
	CLOSE(57)
	CLOSE(58)
	CLOSE(59)
	CLOSE(60)
	CLOSE(205)

	WRITE(*,*)'Output finished'



END SUBROUTINE

!-------------------------------------------------------------------------------------

SUBROUTINE outavg

INTEGER	::i,j,k,islice,jslice,kslice

icenter=i0
jslice=j0
kslice=k0


!WAKE STREAMWISE
	OPEN(unit=60,file='avg_stream_q.dat')
	WRITE(60,*) 'Title= "Average Streamwise q"'
	WRITE(60,*) 'Variables="i","U"'	!"time",
	CALL outinfoheader(60)


	WRITE(60,*)'Zone T="AVG",I=41 ,F= Point'

	DO i=i0+D/2., imax

		WRITE(60,*) (i-i0)/D, q_avg(i,jslice,kslice,1)

	END DO

	CLOSE(60)


!WAKE CROSSWISE
	OPEN(unit=61,file='avg_cross_q.dat')
	WRITE(61,*) 'Title= "Average Cross q"'
	WRITE(61,*) 'Variables="k","U 0.0","U 0.6","U 1.0","U 2.0"'	!"time",
	CALL outinfoheader(61)


	WRITE(61,*)'Zone T="AVG",K=31 ,F= Point'

	DO k=kslice-int(1.5*D), kslice+int(1.5*D)

		WRITE(61,*) (k-kslice)/D, q_avg(i0+int(0.6*D),jslice,k,1),&
			& q_avg(icenter,jslice,k,1),q_avg(icenter+int(D),jslice,k,1),q_avg(icenter+int(2.*D),jslice,k,1)

	END DO

	CLOSE(61)


!VELOCITY
	WRITE(*,*)'Writing Average Results File...'
	OPEN(unit=62, file='avg_xcyl_vel.dat')
	WRITE(62,*) 'Title= "Average Velocity"'
	WRITE(62,*) 'Variables="i","k","U","W","spd"'
	CALL outinfoheader(62)


	WRITE(62,*)'Zone T="TS",I=',imax,',K=',kmax,',F= Point'

j=j0

	DO k=1, kmax
	    DO i=1, imax

	speed=SQRT(MAX(0.00001,(q_avg(i,j,k,1)**2+q_avg(i,j,k,2)**2+q_avg(i,j,k,3)**2)))
	WRITE(62,101) i,k,q_avg(i,j,k,1),q_avg(i,j,k,3),speed

		END DO
	END DO

	CLOSE(62)
	WRITE(*,*)'Finished'



	101 FORMAT(I4,I4,F7.3,F7.3,F7.3)


END SUBROUTINE

!-------------------------------------------------------------------------------------

SUBROUTINE outinfoheader(unit)

INTEGER	::unit


	WRITE(unit,*) 'TEXT X=90,Y=90,H=2,AN=HEADRIGHT,BX=FILLED,BXF=WHITE,LS=1.5,BXM=50'
	WRITE(unit,100) eps_s,epsns_s,xmue,delt,nmax

	100 FORMAT('T="epsf=',F5.3,'\\n epss=',F5.3,'\\n mu=',F5.3,'\\n dt=',F5.3,'\\n nmax=',I5,'"')



END SUBROUTINE

!-------------------------------------------------------------------------------------

SUBROUTINE output !(time,cp)


    DIMENSION om(3) !, cp(ipmax,jpmax,kpmax)
	REAL::speed !,rabs,theta,thetadeg,deli,delk,dist



      WRITE(*,*) " output: saving....."


!LEVELSET
!	WRITE(*,*)'Writing Levelset File...'
!	OPEN(unit=199, file='levelnorm.dat')
!	WRITE(199,*)'Title= "Levelset"'
!	WRITE(199,*)'Variables="i","j","k","Fnorm1","Fnorm2","Fnorm3","F"'

!	WRITE(199,*)'Zone T="TS",I=',imax,',J=',jmax,'
 !    &			,K=',kmax,',F= Point'

!	DO k=1, kmax
!	  DO j=1, jmax
!	    DO i=1, imax

!	WRITE(199,*) i,j,k,Fnorm(i,j,k,1),Fnorm(i,j,k,2)
 !    &,Fnorm(i,j,k,3),F(i,j,k)

!		END DO
!	  END DO
!	END DO

!	CLOSE(199)
!	WRITE(*,*)'Finished'




!VELOCITY
	WRITE(*,*)'Writing Results File...'
	OPEN(unit=200, file='xcyl_vel.dat')
	WRITE(200,*) 'Title= "Velocity"'
	WRITE(200,*) 'Variables="i","k","U","W","F","spd"'
	CALL outinfoheader(200)

	WRITE(200,*)'Zone T="TS",I=',imax,',J=',jmax,',K=',kmax,',F= Point'


j=j0
	DO k=1, kmax
	    DO i=1, imax

	speed=SQRT(MAX(0.00001,(q(i,j,k,1)**2+q(i,j,k,2)**2+q(i,j,k,3)**2)))

	WRITE(200,101) i,k,q(i,j,k,1),q(i,j,k,3),F(i,j,k),speed
	101 FORMAT(I4,I4,F7.3,F7.3,F7.3,F9.3)

		END DO
	END DO

	CLOSE(200)
	WRITE(*,*)'Finished'




!PRESSURE
	WRITE(*,*)'Writing Pressure File...'
	OPEN(unit=201, file='xcyl_cp.dat')
	WRITE(201,*) 'Title= "Cp"'
	WRITE(201,*) 'Variables="ip", "kp",  "cp", "Fp"'
	CALL outinfoheader(201)



	WRITE(201,*) 'Zone I=',ipmax,',K=',kpmax,',F= Point'


jp=j0
	DO kp=1, kpmax
	    DO ip=1, ipmax
			WRITE(201,*) ip, kp, cp(ip,jp,kp), Fp(ip,jp,kp)
		END DO
	END DO

	CLOSE(201)
	WRITE(*,*)'Finished'


!MAG VORTICITY
	WRITE(*,*)'Writing Vorticity File...'
	OPEN(unit=202, file='xcyl_vor.dat')
	WRITE(202,*) 'Title= "Magnitude of Vorticity"'
	WRITE(202,*) 'Variables="i", "k", "vort"'
	CALL outinfoheader(202)



	WRITE(202,*)'Zone T="TS",I=',ipmax,',K=',kpmax,',F= Point'


jp=j0
	DO kp=1, kpmax
	    DO ip=1, ipmax

		WRITE(202,*) ip, kp, absomega(ip,jp,kp)

		END DO
	END DO

	CLOSE(202)
	WRITE(*,*)'Finished'



!VELOCITY
	IF (nswit.NE. 3.) THEN


	WRITE(*,*)'Write initial-data file'
      WRITE(101) imax,jmax,kmax,3
      WRITE(101) (((q(i,j,k,1),i=1,imax),j=1,jmax),k=1,kmax),	&
     &          (((q(i,j,k,2),i=1,imax),j=1,jmax),k=1,kmax),	&
     &          (((q(i,j,k,3),i=1,imax),j=1,jmax),k=1,kmax)
      CLOSE(101)
!
	WRITE(*,*)'Finished'

!  omega
!      WRITE(41) ipmax,jpmax,kpmax,1
!      WRITE(41) (((omega(ip,jp,kp),ip=1,ipmax),jp=1,jpmax),
!     &             kp=1,kpmax)
!      CLOSE(41)

! cp
!      WRITE(40) ipmax,jpmax,kpmax,1
!      WRITE(40) (((cp(ip,jp,kp),ip=1,ipmax),jp=1,jpmax),
!     &             kp=1,kpmax)
!      CLOSE(40)

!delta t
      WRITE(100) time
      CLOSE(100)
	END IF
!
      RETURN
END SUBROUTINE


END MODULE
