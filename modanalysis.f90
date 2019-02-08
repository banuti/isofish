MODULE modanalysis


CONTAINS
!=========================================================

SUBROUTINE pressure !(deltat) !(phi,cp,deltat,qinfinity)
	USE flowprops

	REAL	::phiinf
	INTEGER	::stagoffs

! calculation of pressure coefficient

!      PARAMETER (imax=189,jmax=71,kmax=101,
!     &           ipmax=188,jpmax=70,kpmax=100)

!      DIMENSION phi(ipmax,jpmax,kpmax), cp(ipmax,jpmax,kpmax)

!	COMMON /grid/ delx,dely,delz

      WRITE(*,*) "check: cp"

	stagoffs=2


!pressure correction

!tuning with gridtrafo
	phiinf=0.5*delt*qinfinity*qinfinity+ &
	&(phi(i0-D/2.,j0,k0)+phi(i0-D/2.-1,j0,k0)+phi(i0-D/2.,j0,k0-1)+phi(i0-D/2.-1,j0,k0-1))/4.

WRITE(*,*)'phiinf=',phiinf
!WRITE(*,*)'phi222=',phi(2,2,2)


!zero
!	phiinf=0.

!dt !!! derived from definition of phi
!	phiinf=delt

	cp(:,:,:) = -2.*(phi(:,:,:)-phiinf)/(qinfinity*qinfinity*delt)

!	WRITE(*,*) 'app stag cp: ',cp(i0-D/2,35,50)



	WRITE(*,*) 'stagnation cp field: ',(cp(i0-D/2.-stagoffs,j0,k0)+cp(i0-D/2.-1-stagoffs,j0,k0)&
	&+cp(i0-D/2.-stagoffs,j0,k0-1)+cp(i0-D/2.-1-stagoffs,j0,k0-1))/4.

	WRITE(*,*) 'stagnation cp surfc: ',(cp(i0-D/2.,j0,k0)+cp(i0-D/2.-1,j0,k0)&
	&+cp(i0-D/2.,j0,k0-1)+cp(i0-D/2.-1,j0,k0-1))/4.

	CALL surfacepressure

      RETURN


END SUBROUTINE

!--------------------------------------------------------

SUBROUTINE getspeed

USE flowprops



	spd(:,:,:)=SQRT(MAX(0.00001,(q(:,:,:,1)**2+q(:,:,:,2)**2+q(:,:,:,3)**2)))

	WRITE(*,*)'Max absolute speed: ',maxval(spd)


END SUBROUTINE

!---------------------------------------------------------

SUBROUTINE surfacepressure

USE flowprops

REAL	::rabs,theta,thetadeg,deli,delk,dist
REAL	::cp_bnd,r,icenter, kcenter
INTEGER	::iprime, kprime
INTEGER	::pntindex,jslice

jslice=j0



!pressure distribution in polar coordinates over
!circular cylinder in j-plane



	icenter=i0-0.5
	kcenter=k0-0.5
	


	WRITE(*,*)'Estimate pressure'


	!WRITE(205,*) 'Zone T="TS",theta=',numpnts,',F= Point'



!	ALLOCATE(cppnts(numpnts, numvar))
	!variables: theta,cp (,r)
	
	pntindex=1

	DO k=1,kpmax
	  DO i=1,ipmax
	
		
		IF (Fp(i,jslice,k).LT. 0.) THEN

			kprime=k
			iprime=i

			
			deli=i-icenter
			delk=k-kcenter

			DO iprime=(i-1), (i+1), 2
				IF (Fp(iprime,jslice,k).GT. 0.) THEN

					dist=delx*abs(Fp(i,jslice,k))/(abs(Fp(i,jslice,k))+abs(Fp(iprime,jslice,k)))

					cp_bnd=cp(i,jslice,k)+(dist/delx)*(cp(iprime,jslice,k)-cp(i,jslice,k))
				
					rabs= SQRT( (deli+(iprime-i)*dist)**2 + delk**2 )
					theta=atan2(delk,-(deli+(iprime-i)*dist))
					thetadeg=(180/pi)*theta

					r= SQRT( deli**2 + delk**2 )
!					theta=atan2(delk,-deli)
!					thetadeg=(180/pi)*theta
					IF (thetadeg.LT. 0.) thetadeg=360+thetadeg
					
					cppnts(pntindex,1)=thetadeg
					cppnts(pntindex,2)=cp_bnd
					!cppnts(pntindex,2)=cp(i,jslice,k)
					cppnts(pntindex,3)=rabs
					cppnts(pntindex,4)=r
					pntindex=pntindex+1

					!WRITE(205,*)thetadeg, cp(i,jslice,k)
				END IF
			END DO


			kprime=k
			iprime=i


			DO kprime=(k-1), (k+1), 2
				
				
				IF (Fp(i,jslice,kprime).GT. 0.) THEN

					dist=delz*abs(Fp(i,jslice,k))/(abs(Fp(i,jslice,k))+abs(Fp(i,jslice,kprime)))
					cp_bnd=cp(i,jslice,k)+(dist/delz)*(cp(i,jslice,kprime)-cp(i,jslice,k))
				
					rabs= SQRT( deli**2 + (delk+(kprime-k)*dist)**2 )
					theta=atan2((delk+(kprime-k)*dist),-deli)
					thetadeg=(180/pi)*theta



					r= SQRT( deli**2 + delk**2 )
!					theta=atan2(delk,-deli)
!					thetadeg=(180/pi)*theta
					IF (thetadeg.LT. 0.) thetadeg=360+thetadeg
					
					cppnts(pntindex,1)=thetadeg
					cppnts(pntindex,2)=cp_bnd
					!cppnts(pntindex,2)=cp(i,jslice,k)
					cppnts(pntindex,3)=rabs
					cppnts(pntindex,4)=r
					pntindex=pntindex+1

					!WRITE(205,*)thetadeg, cp(i,jslice,k)
				END IF
			  END DO



		END IF
	  END DO
	END DO
	
	CALL quicksort(numpnts, numvar, cppnts)
		


END SUBROUTINE

!---------------------------------------------------------

SUBROUTINE drag

USE flowprops

INTEGER	::i
REAL	::cp_loc(numpnts),theta(numpnts),Cdprime


WRITE(*,*)'numpnts: ',numpnts

Cdprime=0.
cp_loc(:)		= cppnts(:,2)
theta(:)	= pi/180.*cppnts(:,1) 



	Cdprime=Cdprime+( cp_loc(1)*COS(theta(1)/2.)*theta(1) )


DO i=1,numpnts-1

	Cdprime=Cdprime+((cp_loc(i)+cp_loc(i+1))/2.*COS((theta(i)+theta(i+1))/2.)*(theta(i+1)-theta(i)))

END DO

	Cdprime=Cdprime+(cp_loc(numpnts)*COS((theta(numpnts)+2*pi)/2.)*(2*pi-theta(numpnts)))

	Cd=0.5*Cdprime 


END SUBROUTINE

!---------------------------------------------------------

SUBROUTINE lift

USE flowprops

!	cppnts(pntindex,1)=thetadeg
!	cppnts(pntindex,2)=cp_bnd

INTEGER	::i
REAL	::cp_loc(numpnts),theta(numpnts),Clprime

Clprime=0.
cp_loc(:)		= cppnts(:,2)
theta(:)	= pi/180.*cppnts(:,1) 



	Clprime=Clprime+( cp_loc(1)*SIN(theta(1)/2.)*theta(1) )


DO i=1,numpnts-1

	Clprime=Clprime+((cp_loc(i)+cp_loc(i+1))/2.*SIN((theta(i)+theta(i+1))/2.)*(theta(i+1)-theta(i)))

END DO

	Clprime=Clprime+(cp_loc(numpnts)*SIN((theta(numpnts)+2*pi)/2.)*(2*pi-theta(numpnts)))

	Cl=-0.5*Clprime !times radius


END SUBROUTINE


!---------------------------------------------------------------------------
SUBROUTINE quicksort(n,var,arr)

	INTEGER n,M,NSTACK,var !n is length,var is number of variables(theta,r,cp.
	REAL arr(n,var)
	PARAMETER (M=7,NSTACK=50)
!Sorts an array arr(1:n) into ascending numerical order using the Quicksort algorithm. n
!is input; arr is replaced on output by its sorted rearrangement.
!Parameters: M is the size of subarrays sorted by straight insertion and NSTACK is the required
!auxiliary storage.
	INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
	REAL a(var),temp(var)
	jstack=0
	l=1
	ir=n
1	if(ir-l.lt.M)then !Insertion sort when subarray small enough.
	do j=l+1,ir
	  a(:)=arr(j,:)
	  do i=j-1,l,-1
	  if(arr(i,1).le. a(1)) goto 2
		arr(i+1,:)=arr(i,:)
	  enddo
	  i=l-1
2	arr(i+1,:)=a(:)
	enddo

	if(jstack.eq.0)return

	ir=istack(jstack) !Pop stack and begin a new round of partitioning.
	l=istack(jstack-1)
	jstack=jstack-2
	else
	k=(l+ir)/2	!Choose median of left, center, and right elements as partitioning
				!element a. Also rearrange so that a(l) . a(l+1) . a(ir).
	temp(:)=arr(k,:)
	arr(k,:)=arr(l+1,:)
	arr(l+1,:)=temp(:)
	if(arr(l,1).gt.arr(ir,1))then
	  temp(:)=arr(l,:)
	  arr(l,:)=arr(ir,:)
	  arr(ir,:)=temp(:)
	endif
	
	if(arr(l+1,1).gt.arr(ir,1))then
	  temp(:)=arr(l+1,:)
	  arr(l+1,:)=arr(ir,:)
	  arr(ir,:)=temp(:)
	endif

	if(arr(l,1).gt.arr(l+1,1))then
	  temp(:)=arr(l,:)
	  arr(l,:)=arr(l+1,:)
	  arr(l+1,:)=temp(:)
	endif

	i=l+1 !Initialize pointers for partitioning.
	j=ir
	a(:)=arr(l+1,:) !Partitioning element.
3	continue !Beginning of innermost loop.

	i=i+1 !Scan up to nd element > a.
	if(arr(i,1).lt.a(1))goto 3
4	continue
	j=j-1 !Scan down to nd element < a.
	
	if(arr(j,1).gt.a(1))goto 4
	if(j.lt.i)goto 5 !Pointers crossed. Exit with partitioning complete.
	
	temp(:)=arr(i,:) !Exchange elements.
	arr(i,:)=arr(j,:)
	arr(j,:)=temp(:)
	goto 3 !End of innermost loop.
5	arr(l+1,:)=arr(j,:) !Insert partitioning element.
	arr(j,:)=a(:)
	jstack=jstack+2
!Push pointers to larger subarray on stack, process smaller subarray immediately.
	if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
	if(ir-i+1.ge.j-l)then
	  istack(jstack)=ir
	  istack(jstack-1)=i
	  ir=j-1
	else
	  istack(jstack)=j-1
	  istack(jstack-1)=l
	  l=i
	endif
	endif
	goto 1

END SUBROUTINE
!=========================================================
END MODULE