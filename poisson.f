C FISHPAK12  FROM PORTLIB                                  03/11/81             
      SUBROUTINE HW3CRT (XS,XF,L,LBDCND,BDXS,BDXF,YS,YF,M,MBDCND,BDYS,          
     1                   BDYF,ZS,ZF,N,NBDCND,BDZS,BDZF,ELMBDA,LDIMF,            
     2                   MDIMF,F,PERTRB,IERROR,W)                               
C                                                                     
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
C     *                                                               *         
C     *                        F I S H P A K                          *         
C     *                                                               *         
C     *                                                               *         
C     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *         
C     *                                                               *         
C     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *         
C     *                                                               *         
C     *                  (VERSION 3.1 , OCTOBER 1980)                 *        
C     *                                                               *         
C     *                             BY                                *         
C     *                                                               *         
C     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *         
C     *                                                               *         
C     *                             OF                                *         
C     *                                                               *         
C     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *         
C     *                                                               *         
C     *                BOULDER, COLORADO  (80307)  U.S.A.             *         
C     *                                                               *         
C     *                   WHICH IS SPONSORED BY                       *         
C     *                                                               *         
C     *              THE NATIONAL SCIENCE FOUNDATION                  *         
C     *                                                               *         
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
C                                                                               
C                                                                               
C    * * * * * * * * *  PURPOSE    * * * * * * * * * * * * * * * * * *          
C                                                                               
C          SUBROUTINE HW3CRT SOLVES THE STANDARD SEVEN-POINT FINITE             
C     DIFFERENCE APPROXIMATION TO THE HELMHOLTZ EQUATION IN CARTESIAN           
C     COORDINATES:                                                              
C                                                                               
C         (D/DX)(DU/DX) + (D/DY)(DU/DY) + (D/DZ)(DU/DZ)                         
C                                                                               
C                    + LAMBDA*U = F(X,Y,Z) .                                    
C                                                                               
C    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *          
C                                                                               
C                                                                               
C    * * * * * * * *    PARAMETER DESCRIPTION     * * * * * * * * * *           
C                                                                               
C                                                                               
C            * * * * * *   ON INPUT    * * * * * *                              
C                                                                               
C     XS,XF                                                                     
C        THE RANGE OF X, I.E. XS .LE. X .LE. XF .                               
C        XS MUST BE LESS THAN XF.                                               
C                                                                               
C     L                                                                         
C        THE NUMBER OF PANELS INTO WHICH THE INTERVAL (XS,XF) IS                
C        SUBDIVIDED.  HENCE, THERE WILL BE L+1 GRID POINTS IN THE               
C        X-DIRECTION GIVEN BY X(I) = XS+(I-1)DX FOR I=1,2,...,L+1,              
C        WHERE DX = (XF-XS)/L IS THE PANEL WIDTH.  L MUST BE AT                 
C        LEAST 5 .                                                              
C                                                                               
C     LBDCND                                                                    
C        INDICATES THE TYPE OF BOUNDARY CONDITIONS AT X = XS AND X = XF.        
C                                                                               
C        = 0  IF THE SOLUTION IS PERIODIC IN X, I.E.                            
C             U(L+I,J,K) = U(I,J,K).                                            
C        = 1  IF THE SOLUTION IS SPECIFIED AT X = XS AND X = XF.                
C        = 2  IF THE SOLUTION IS SPECIFIED AT X = XS AND THE DERIVATIVE         
C             OF THE SOLUTION WITH RESPECT TO X IS SPECIFIED AT X = XF.         
C        = 3  IF THE DERIVATIVE OF THE SOLUTION WITH RESPECT TO X IS            
C             SPECIFIED AT X = XS AND X = XF.                                   
C        = 4  IF THE DERIVATIVE OF THE SOLUTION WITH RESPECT TO X IS            
C             SPECIFIED AT X = XS AND THE SOLUTION IS SPECIFIED AT X=XF.        
C                                                                               
C     BDXS                                                                      
C        A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE VALUES OF THE               
C        DERIVATIVE OF THE SOLUTION WITH RESPECT TO X AT X = XS.                
C        WHEN LBDCND = 3 OR 4,                                                  
C                                                                               
C             BDXS(J,K) = (D/DX)U(XS,Y(J),Z(K)), J=1,2,...,M+1,                 
C                                                K=1,2,...,N+1.                 
C                                                                               
C        WHEN LBDCND HAS ANY OTHER VALUE, BDXS IS A DUMMY VARIABLE.             
C        BDXS MUST BE DIMENSIONED AT LEAST (M+1)*(N+1).                         
C                                                                               
C     BDXF                                                                      
C        A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE VALUES OF THE               
C        DERIVATIVE OF THE SOLUTION WITH RESPECT TO X AT X = XF.                
C        WHEN LBDCND = 2 OR 3,                                                  
C                                                                               
C             BDXF(J,K) = (D/DX)U(XF,Y(J),Z(K)), J=1,2,...,M+1,                 
C                                                K=1,2,...,N+1.                 
C                                                                               
C        WHEN LBDCND HAS ANY OTHER VALUE, BDXF IS A DUMMY VARIABLE.             
C        BDXF MUST BE DIMENSIONED AT LEAST (M+1)*(N+1).                         
C                                                                               
C     YS,YF                                                                     
C        THE RANGE OF Y, I.E. YS .LE. Y .LE. YF.                                
C        YS MUST BE LESS THAN YF.                                               
C                                                                               
C     M                                                                         
C        THE NUMBER OF PANELS INTO WHICH THE INTERVAL (YS,YF) IS                
C        SUBDIVIDED.  HENCE, THERE WILL BE M+1 GRID POINTS IN THE               
C        Y-DIRECTION GIVEN BY Y(J) = YS+(J-1)DY FOR J=1,2,...,M+1,              
C        WHERE DY = (YF-YS)/M IS THE PANEL WIDTH.  M MUST BE AT                 
C        LEAST 5 .                                                              
C                                                                               
C     MBDCND                                                                    
C        INDICATES THE TYPE OF BOUNDARY CONDITIONS AT Y = YS AND Y = YF.        
C                                                                               
C        = 0  IF THE SOLUTION IS PERIODIC IN Y, I.E.                            
C             U(I,M+J,K) = U(I,J,K).                                            
C        = 1  IF THE SOLUTION IS SPECIFIED AT Y = YS AND Y = YF.                
C        = 2  IF THE SOLUTION IS SPECIFIED AT Y = YS AND THE DERIVATIVE         
C             OF THE SOLUTION WITH RESPECT TO Y IS SPECIFIED AT Y = YF.         
C        = 3  IF THE DERIVATIVE OF THE SOLUTION WITH RESPECT TO Y IS            
C             SPECIFIED AT Y = YS AND Y = YF.                                   
C        = 4  IF THE DERIVATIVE OF THE SOLUTION WITH RESPECT TO Y IS            
C             SPECIFIED AT Y = YS AND THE SOLUTION IS SPECIFIED AT Y=YF.        
C                                                                               
C     BDYS                                                                      
C        A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE VALUES OF THE               
C        DERIVATIVE OF THE SOLUTION WITH RESPECT TO Y AT Y = YS.                
C        WHEN MBDCND = 3 OR 4,                                                  
C                                                                               
C             BDYS(I,K) = (D/DY)U(X(I),YS,Z(K)), I=1,2,...,L+1,                 
C                                                K=1,2,...,N+1.                 
C                                                                               
C        WHEN MBDCND HAS ANY OTHER VALUE, BDYS IS A DUMMY VARIABLE.             
C        BDYS MUST BE DIMENSIONED AT LEAST (L+1)*(N+1).                         
C                                                                               
C     BDYF                                                                      
C        A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE VALUES OF THE               
C        DERIVATIVE OF THE SOLUTION WITH RESPECT TO Y AT Y = YF.                
C        WHEN MBDCND = 2 OR 3,                                                  
C                                                                               
C             BDYF(I,K) = (D/DY)U(X(I),YF,Z(K)), I=1,2,...,L+1,                 
C                                                K=1,2,...,N+1.                 
C                                                                               
C        WHEN MBDCND HAS ANY OTHER VALUE, BDYF IS A DUMMY VARIABLE.             
C        BDYF MUST BE DIMENSIONED AT LEAST (L+1)*(N+1).                         
C                                                                               
C     ZS,ZF                                                                     
C        THE RANGE OF Z, I.E. ZS .LE. Z .LE. ZF.                                
C        ZS MUST BE LESS THAN ZF.                                               
C                                                                               
C     N                                                                         
C        THE NUMBER OF PANELS INTO WHICH THE INTERVAL (ZS,ZF) IS                
C        SUBDIVIDED.  HENCE, THERE WILL BE N+1 GRID POINTS IN THE               
C        Z-DIRECTION GIVEN BY Z(K) = ZS+(K-1)DZ FOR K=1,2,...,N+1,              
C        WHERE DZ = (ZF-ZS)/N IS THE PANEL WIDTH.  N MUST BE AT LEAST 5.        
C                                                                               
C     NBDCND                                                                    
C        INDICATES THE TYPE OF BOUNDARY CONDITIONS AT Z = ZS AND Z = ZF.        
C                                                                               
C        = 0  IF THE SOLUTION IS PERIODIC IN Z, I.E.                            
C             U(I,J,N+K) = U(I,J,K).                                            
C        = 1  IF THE SOLUTION IS SPECIFIED AT Z = ZS AND Z = ZF.                
C        = 2  IF THE SOLUTION IS SPECIFIED AT Z = ZS AND THE DERIVATIVE         
C             OF THE SOLUTION WITH RESPECT TO Z IS SPECIFIED AT Z = ZF.         
C        = 3  IF THE DERIVATIVE OF THE SOLUTION WITH RESPECT TO Z IS            
C             SPECIFIED AT Z = ZS AND Z = ZF.                                   
C        = 4  IF THE DERIVATIVE OF THE SOLUTION WITH RESPECT TO Z IS            
C             SPECIFIED AT Z = ZS AND THE SOLUTION IS SPECIFIED AT Z=ZF.        
C                                                                               
C     BDZS                                                                      
C        A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE VALUES OF THE               
C        DERIVATIVE OF THE SOLUTION WITH RESPECT TO Z AT Z = ZS.                
C        WHEN NBDCND = 3 OR 4,                                                  
C                                                                               
C             BDZS(I,J) = (D/DZ)U(X(I),Y(J),ZS), I=1,2,...,L+1,                 
C                                                J=1,2,...,M+1.                 
C                                                                               
C        WHEN NBDCND HAS ANY OTHER VALUE, BDZS IS A DUMMY VARIABLE.             
C        BDZS MUST BE DIMENSIONED AT LEAST (L+1)*(M+1).                         
C                                                                               
C     BDZF                                                                      
C        A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE VALUES OF THE               
C        DERIVATIVE OF THE SOLUTION WITH RESPECT TO Z AT Z = ZF.                
C        WHEN NBDCND = 2 OR 3,                                                  
C                                                                               
C             BDZF(I,J) = (D/DZ)U(X(I),Y(J),ZF), I=1,2,...,L+1,                 
C                                                J=1,2,...,M+1.                 
C                                                                               
C        WHEN NBDCND HAS ANY OTHER VALUE, BDZF IS A DUMMY VARIABLE.             
C        BDZF MUST BE DIMENSIONED AT LEAST (L+1)*(M+1).                         
C                                                                               
C     ELMBDA                                                                    
C        THE CONSTANT LAMBDA IN THE HELMHOLTZ EQUATION. IF                      
C        LAMBDA .GT. 0, A SOLUTION MAY NOT EXIST.  HOWEVER, HW3CRT WILL         
C        ATTEMPT TO FIND A SOLUTION.                                            
C                                                                               
C     F                                                                         
C        A THREE-DIMENSIONAL ARRAY THAT SPECIFIES THE VALUES OF THE             
C        RIGHT SIDE OF THE HELMHOLTZ EQUATION AND BOUNDARY VALUES (IF           
C        ANY).  FOR I=2,3,...,L, J=2,3,...,M, AND K=2,3,...,N                   
C                                                                               
C                   F(I,J,K) = F(X(I),Y(J),Z(K)).                               
C                                                                               
C        ON THE BOUNDARIES F IS DEFINED BY                                      
C                                                                               
C        LBDCND      F(1,J,K)         F(L+1,J,K)                                
C        ------   ---------------   ---------------                             
C                                                                               
C          0      F(XS,Y(J),Z(K))   F(XS,Y(J),Z(K))                             
C          1      U(XS,Y(J),Z(K))   U(XF,Y(J),Z(K))                             
C          2      U(XS,Y(J),Z(K))   F(XF,Y(J),Z(K))   J=1,2,...,M+1             
C          3      F(XS,Y(J),Z(K))   F(XF,Y(J),Z(K))   K=1,2,...,N+1             
C          4      F(XS,Y(J),Z(K))   U(XF,Y(J),Z(K))                             
C                                                                               
C        MBDCND      F(I,1,K)         F(I,M+1,K)                                
C        ------   ---------------   ---------------                             
C                                                                               
C          0      F(X(I),YS,Z(K))   F(X(I),YS,Z(K))                             
C          1      U(X(I),YS,Z(K))   U(X(I),YF,Z(K))                             
C          2      U(X(I),YS,Z(K))   F(X(I),YF,Z(K))   I=1,2,...,L+1             
C          3      F(X(I),YS,Z(K))   F(X(I),YF,Z(K))   K=1,2,...,N+1             
C          4      F(X(I),YS,Z(K))   U(X(I),YF,Z(K))                             
C                                                                               
C        NBDCND      F(I,J,1)         F(I,J,N+1)                                
C        ------   ---------------   ---------------                             
C                                                                               
C          0      F(X(I),Y(J),ZS)   F(X(I),Y(J),ZS)                             
C          1      U(X(I),Y(J),ZS)   U(X(I),Y(J),ZF)                             
C          2      U(X(I),Y(J),ZS)   F(X(I),Y(J),ZF)   I=1,2,...,L+1             
C          3      F(X(I),Y(J),ZS)   F(X(I),Y(J),ZF)   J=1,2,...,M+1             
C          4      F(X(I),Y(J),ZS)   U(X(I),Y(J),ZF)                             
C                                                                               
C        F MUST BE DIMENSIONED AT LEAST (L+1)*(M+1)*(N+1).                      
C                                                                               
C        NOTE:                                                                  
C                                                                               
C        IF THE TABLE CALLS FOR BOTH THE SOLUTION U AND THE RIGHT SIDE F        
C        ON A BOUNDARY, THEN THE SOLUTION MUST BE SPECIFIED.                    
C                                                                               
C     LDIMF                                                                     
C        THE ROW (OR FIRST) DIMENSION OF THE ARRAYS F,BDYS,BDYF,BDZS,           
C        AND BDZF AS IT APPEARS IN THE PROGRAM CALLING HW3CRT. THIS             
C        PARAMETER IS USED TO SPECIFY THE VARIABLE DIMENSION OF THESE           
C        ARRAYS.  LDIMF MUST BE AT LEAST L+1.                                   
C                                                                               
C     MDIMF                                                                     
C        THE COLUMN (OR SECOND) DIMENSION OF THE ARRAY F AND THE ROW (OR        
C        FIRST) DIMENSION OF THE ARRAYS BDXS AND BDXF AS IT APPEARS IN          
C        THE PROGRAM CALLING HW3CRT.  THIS PARAMETER IS USED TO SPECIFY         
C        THE VARIABLE DIMENSION OF THESE ARRAYS.                                
C        MDIMF MUST BE AT LEAST M+1.                                            
C                                                                               
C     W                                                                         
C        A ONE-DIMENSIONAL ARRAY THAT MUST BE PROVIDED BY THE USER FOR          
C        WORK SPACE.  THE LENGTH OF W MUST BE AT LEAST 30 + L + M + 5*N         
C        + MAX(L,M,N) + 7*(INT((L+1)/2) + INT((M+1)/2))                         
C                                                                               
C                                                                               
C            * * * * * *   ON OUTPUT   * * * * * *                              
C                                                                               
C     F                                                                         
C        CONTAINS THE SOLUTION U(I,J,K) OF THE FINITE DIFFERENCE                
C        APPROXIMATION FOR THE GRID POINT (X(I),Y(J),Z(K)) FOR                  
C        I=1,2,...,L+1, J=1,2,...,M+1, AND K=1,2,...,N+1.                       
C                                                                               
C     PERTRB                                                                    
C        IF A COMBINATION OF PERIODIC OR DERIVATIVE BOUNDARY CONDITIONS         
C        IS SPECIFIED FOR A POISSON EQUATION (LAMBDA = 0), A SOLUTION           
C        MAY NOT EXIST.  PERTRB IS A CONSTANT, CALCULATED AND SUBTRACTED        
C        FROM F, WHICH ENSURES THAT A SOLUTION EXISTS.  PWSCRT THEN             
C        COMPUTES THIS SOLUTION, WHICH IS A LEAST SQUARES SOLUTION TO           
C        THE ORIGINAL APPROXIMATION.  THIS SOLUTION IS NOT UNIQUE AND IS        
C        UNNORMALIZED.  THE VALUE OF PERTRB SHOULD BE SMALL COMPARED TO         
C        THE RIGHT SIDE F.  OTHERWISE, A SOLUTION IS OBTAINED TO AN             
C        ESSENTIALLY DIFFERENT PROBLEM.  THIS COMPARISON SHOULD ALWAYS          
C        BE MADE TO INSURE THAT A MEANINGFUL SOLUTION HAS BEEN OBTAINED.        
C                                                                               
C     IERROR                                                                    
C        AN ERROR FLAG THAT INDICATES INVALID INPUT PARAMETERS.  EXCEPT         
C        FOR NUMBERS 0 AND 12, A SOLUTION IS NOT ATTEMPTED.                     
C                                                                               
C        =  0  NO ERROR                                                         
C        =  1  XS .GE. XF                                                       
C        =  2  L .LT. 5                                                         
C        =  3  LBDCND .LT. 0 .OR. LBDCND .GT. 4                                 
C        =  4  YS .GE. YF                                                       
C        =  5  M .LT. 5                                                         
C        =  6  MBDCND .LT. 0 .OR. MBDCND .GT. 4                                 
C        =  7  ZS .GE. ZF                                                       
C        =  8  N .LT. 5                                                         
C        =  9  NBDCND .LT. 0 .OR. NBDCND .GT. 4                                 
C        = 10  LDIMF .LT. L+1                                                   
C        = 11  MDIMF .LT. M+1                                                   
C        = 12  LAMBDA .GT. 0                                                    
C                                                                               
C        SINCE THIS IS THE ONLY MEANS OF INDICATING A POSSIBLY INCORRECT        
C        CALL TO HW3CRT, THE USER SHOULD TEST IERROR AFTER THE CALL.            
C                                                                               
C                                                                               
C    * * * * * * *   PROGRAM SPECIFICATIONS    * * * * * * * * * * * *          
C                                                                               
C     DIMENSION OF   BDXS(MDIMF,N+1),BDXF(MDIMF,N+1),BDYS(LDIMF,N+1),           
C     ARGUMENTS      BDYF(LDIMF,N+1),BDZS(LDIMF,M+1),BDZF(LDIMF,M+1),           
C                    F(LDIMF,MDIMF,N+1),W(SEE ARGUMENT LIST)                    
C                                                                               
C     LATEST         DECEMBER 1, 1978                                           
C     REVISION                                                                  
C                                                                               
C     SUBPROGRAMS    HW3CRT,POIS3D,POS3D1,TRID,RFFTI,RFFTF,RFFTF1,              
C     REQUIRED       RFFTB,RFFTB1,COSTI,COST,SINTI,SINT,COSQI,COSQF,            
C                    COSQF1,COSQB,COSQB1,SINQI,SINQF,SINQB,CFFTI,               
C                    CFFTI1,CFFTB,CFFTB1,PASSB2,PASSB3,PASSB4,PASSB,            
C                    CFFTF,CFFTF1,PASSF1,PASSF2,PASSF3,PASSF4,PASSF,            
C                    PIMACH                                                     
C                                                                               
C     SPECIAL        NONE                                                       
C     CONDITIONS                                                                
C                                                                               
C     COMMON         VALUE                                                      
C     BLOCKS                                                                    
C                                                                               
C     I/O            NONE                                                       
C                                                                               
C     PRECISION      SINGLE                                                     
C                                                                               
C     SPECIALIST     ROLAND SWEET                                               
C                                                                               
C     LANGUAGE       FORTRAN                                                    
C                                                                               
C     HISTORY        WRITTEN BY ROLAND SWEET AT NCAR IN JULY,1977               
C                                                                               
C     ALGORITHM      THIS SUBROUTINE DEFINES THE FINITE DIFFERENCE              
C                    EQUATIONS, INCORPORATES BOUNDARY DATA, AND                 
C                    ADJUSTS THE RIGHT SIDE OF SINGULAR SYSTEMS AND             
C                    THEN CALLS POIS3D TO SOLVE THE SYSTEM.                     
C                                                                               
C     SPACE          7862(DECIMAL) = 17300(OCTAL) LOCATIONS ON THE              
C     REQUIRED       NCAR CONTROL DATA 7600                                     
C                                                                               
C     TIMING AND        THE EXECUTION TIME T ON THE NCAR CONTROL DATA           
C     ACCURACY       7600 FOR SUBROUTINE HW3CRT IS ROUGHLY PROPORTIONAL         
C                    TO L*M*N*(LOG2(L)+LOG2(M)+5), BUT ALSO DEPENDS ON          
C                    INPUT PARAMETERS LBDCND AND MBDCND.  SOME TYPICAL          
C                    VALUES ARE LISTED IN THE TABLE BELOW.                      
C                       THE SOLUTION PROCESS EMPLOYED RESULTS IN A LOSS         
C                    OF NO MORE THAN THREE SIGNIFICANT DIGITS FOR L,M AN        
C                    N AS LARGE AS 32.  MORE DETAILED INFORMATION ABOUT         
C                    ACCURACY CAN BE FOUND IN THE DOCUMENTATION FOR             
C                    SUBROUTINE POIS3D WHICH IS THE ROUTINE THAT ACTUALL        
C                    SOLVES THE FINITE DIFFERENCE EQUATIONS.                    
C                                                                               
C                                                                               
C                       L(=M=N)     LBDCND(=MBDCND=NBDCND)      T(MSECS)        
C                       -------     ----------------------      --------        
C                                                                               
C                         16                  0                    300          
C                         16                  1                    302          
C                         16                  3                    348          
C                         32                  0                   1925          
C                         32                  1                   1929          
C                         32                  3                   2109          
C                                                                               
C     PORTABILITY    AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN.             
C                    THE MACHINE DEPENDENT CONSTANT PI IS DEFINED IN            
C                    FUNCTION PIMACH.                                           
C                                                                               
C     REQUIRED       COS,SIN,ATAN                                               
C     RESIDENT                                                                  
C     ROUTINES                                                                  
C                                                                               
C     REFERENCE      NONE                                                       
C                                                                               
C     REQUIRED         COS,SIN,ATAN                                             
C     RESIDENT                                                                  
C     ROUTINES                                                                  
C                                                                               
C     REFERENCE        NONE                                                     
C                                                                               
C    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *          
C                                                                               
      DIMENSION       BDXS(MDIMF,1)          ,BDXF(MDIMF,1)          ,          
     1                BDYS(LDIMF,1)          ,BDYF(LDIMF,1)          ,          
     2                BDZS(LDIMF,1)          ,BDZF(LDIMF,1)          ,          
     3                F(LDIMF,MDIMF,1)       ,W(1)                              
C                                                                               
C     CHECK FOR INVALID INPUT.                                                  
C                                                                               
      IERROR = 0                                                                
      IF (XF .LE. XS) IERROR = 1                                                
      IF (L .LT. 5) IERROR = 2                                                  
      IF (LBDCND.LT.0 .OR. LBDCND.GT.4) IERROR = 3                              
      IF (YF .LE. YS) IERROR = 4                                                
      IF (M .LT. 5) IERROR = 5                                                  
      IF (MBDCND.LT.0 .OR. MBDCND.GT.4) IERROR = 6                              
      IF (ZF .LE. ZS) IERROR = 7                                                
      IF (N .LT. 5) IERROR = 8                                                  
      IF (NBDCND.LT.0 .OR. NBDCND.GT.4) IERROR = 9                              
      IF (LDIMF .LT. L+1) IERROR = 10                                           
      IF (MDIMF .LT. M+1) IERROR = 11                                           
      IF (IERROR .NE. 0) GO TO 188                                              
      DY = (YF-YS)/FLOAT(M)                                                     
      TWBYDY = 2./DY                                                            
      C2 = 1./(DY**2)                                                           
      MSTART = 1                                                                
      MSTOP = M                                                                 
      MP1 = M+1                                                                 
      MP = MBDCND+1                                                             
      GO TO (104,101,101,102,102),MP                                            
  101 MSTART = 2                                                                
  102 GO TO (104,104,103,103,104),MP                                            
  103 MSTOP = MP1                                                               
  104 MUNK = MSTOP-MSTART+1                                                     
      DZ = (ZF-ZS)/FLOAT(N)                                                     
      TWBYDZ = 2./DZ                                                            
      NP = NBDCND+1                                                             
      C3 = 1./(DZ**2)                                                           
      NP1 = N+1                                                                 
      NSTART = 1                                                                
      NSTOP = N                                                                 
      GO TO (108,105,105,106,106),NP                                            
  105 NSTART = 2                                                                
  106 GO TO (108,108,107,107,108),NP                                            
  107 NSTOP = NP1                                                               
  108 NUNK = NSTOP-NSTART+1                                                     
      LP1 = L+1                                                                 
      DX = (XF-XS)/FLOAT(L)                                                     
      C1 = 1./(DX**2)                                                           
      TWBYDX = 2./DX                                                            
      LP = LBDCND+1                                                             
      LSTART = 1                                                                
      LSTOP = L                                                                 
C                                                                               
C     ENTER BOUNDARY DATA FOR X-BOUNDARIES.                                     
C                                                                               
      GO TO (122,109,109,112,112),LP                                            
  109 LSTART = 2                                                                
      DO 111 J=MSTART,MSTOP                                                     
         DO 110 K=NSTART,NSTOP                                                  
            F(2,J,K) = F(2,J,K)-C1*F(1,J,K)                                     
  110    CONTINUE                                                               
  111 CONTINUE                                                                  
      GO TO 115                                                                 
  112 DO 114 J=MSTART,MSTOP                                                     
         DO 113 K=NSTART,NSTOP                                                  
            F(1,J,K) = F(1,J,K)+TWBYDX*BDXS(J,K)                                
  113    CONTINUE                                                               
  114 CONTINUE                                                                  
  115 GO TO (122,116,119,119,116),LP                                            
  116 DO 118 J=MSTART,MSTOP                                                     
         DO 117 K=NSTART,NSTOP                                                  
            F(L,J,K) = F(L,J,K)-C1*F(LP1,J,K)                                   
  117    CONTINUE                                                               
  118 CONTINUE                                                                  
      GO TO 122                                                                 
  119 LSTOP = LP1                                                               
      DO 121 J=MSTART,MSTOP                                                     
         DO 120 K=NSTART,NSTOP                                                  
            F(LP1,J,K) = F(LP1,J,K)-TWBYDX*BDXF(J,K)                            
  120    CONTINUE                                                               
  121 CONTINUE                                                                  
  122 LUNK = LSTOP-LSTART+1                                                     
C                                                                               
C     ENTER BOUNDARY DATA FOR Y-BOUNDARIES.                                     
C                                                                               
      GO TO (136,123,123,126,126),MP                                            
  123 DO 125 I=LSTART,LSTOP                                                     
         DO 124 K=NSTART,NSTOP                                                  
            F(I,2,K) = F(I,2,K)-C2*F(I,1,K)                                     
  124    CONTINUE                                                               
  125 CONTINUE                                                                  
      GO TO 129                                                                 
  126 DO 128 I=LSTART,LSTOP                                                     
         DO 127 K=NSTART,NSTOP                                                  
            F(I,1,K) = F(I,1,K)+TWBYDY*BDYS(I,K)                                
  127    CONTINUE                                                               
  128 CONTINUE                                                                  
  129 GO TO (136,130,133,133,130),MP                                            
  130 DO 132 I=LSTART,LSTOP                                                     
         DO 131 K=NSTART,NSTOP                                                  
            F(I,M,K) = F(I,M,K)-C2*F(I,MP1,K)                                   
  131    CONTINUE                                                               
  132 CONTINUE                                                                  
      GO TO 136                                                                 
  133 DO 135 I=LSTART,LSTOP                                                     
         DO 134 K=NSTART,NSTOP                                                  
            F(I,MP1,K) = F(I,MP1,K)-TWBYDY*BDYF(I,K)                            
  134    CONTINUE                                                               
  135 CONTINUE                                                                  
  136 CONTINUE                                                                  
C                                                                               
C     ENTER BOUNDARY DATA FOR Z-BOUNDARIES.                                     
C                                                                               
      GO TO (150,137,137,140,140),NP                                            
  137 DO 139 I=LSTART,LSTOP                                                     
         DO 138 J=MSTART,MSTOP                                                  
            F(I,J,2) = F(I,J,2)-C3*F(I,J,1)                                     
  138    CONTINUE                                                               
  139 CONTINUE                                                                  
      GO TO 143                                                                 
  140 DO 142 I=LSTART,LSTOP                                                     
         DO 141 J=MSTART,MSTOP                                                  
            F(I,J,1) = F(I,J,1)+TWBYDZ*BDZS(I,J)                                
  141    CONTINUE                                                               
  142 CONTINUE                                                                  
  143 GO TO (150,144,147,147,144),NP                                            
  144 DO 146 I=LSTART,LSTOP                                                     
         DO 145 J=MSTART,MSTOP                                                  
            F(I,J,N) = F(I,J,N)-C3*F(I,J,NP1)                                   
  145    CONTINUE                                                               
  146 CONTINUE                                                                  
      GO TO 150                                                                 
  147 DO 149 I=LSTART,LSTOP                                                     
         DO 148 J=MSTART,MSTOP                                                  
            F(I,J,NP1) = F(I,J,NP1)-TWBYDZ*BDZF(I,J)                            
  148    CONTINUE                                                               
  149 CONTINUE                                                                  
C                                                                               
C     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.                                     
C                                                                               
  150 CONTINUE                                                                  
      IWB = NUNK+1                                                              
      IWC = IWB+NUNK                                                            
      IWW = IWC+NUNK                                                            
      DO 151 K=1,NUNK                                                           
         I = IWC+K-1                                                            
         W(K) = C3                                                              
         W(I) = C3                                                              
         I = IWB+K-1                                                            
         W(I) = -2.*C3+ELMBDA                                                   
  151 CONTINUE                                                                  
      GO TO (155,155,153,152,152),NP                                            
  152 W(IWC) = 2.*C3                                                            
  153 GO TO (155,155,154,154,155),NP                                            
  154 W(IWB-1) = 2.*C3                                                          
  155 CONTINUE                                                                  
      PERTRB = 0.                                                               
C                                                                               
C     FOR SINGULAR PROBLEMS ADJUST DATA TO INSURE A SOLUTION WILL EXIST.        
C                                                                               
      GO TO (156,172,172,156,172),LP                                            
  156 GO TO (157,172,172,157,172),MP                                            
  157 GO TO (158,172,172,158,172),NP                                            
  158 IF (ELMBDA) 172,160,159                                                   
  159 IERROR = 12                                                               
      GO TO 172                                                                 
  160 CONTINUE                                                                  
      MSTPM1 = MSTOP-1                                                          
      LSTPM1 = LSTOP-1                                                          
      NSTPM1 = NSTOP-1                                                          
      XLP = (2+LP)/3                                                            
      YLP = (2+MP)/3                                                            
      ZLP = (2+NP)/3                                                            
      S1 = 0.                                                                   
      DO 164 K=2,NSTPM1                                                         
         DO 162 J=2,MSTPM1                                                      
            DO 161 I=2,LSTPM1                                                   
               S1 = S1+F(I,J,K)                                                 
  161       CONTINUE                                                            
            S1 = S1+(F(1,J,K)+F(LSTOP,J,K))/XLP                                 
  162    CONTINUE                                                               
         S2 = 0.                                                                
         DO 163 I=2,LSTPM1                                                      
            S2 = S2+F(I,1,K)+F(I,MSTOP,K)                                       
  163    CONTINUE                                                               
         S2 = (S2+(F(1,1,K)+F(1,MSTOP,K)+F(LSTOP,1,K)+F(LSTOP,MSTOP,K))/        
     1                                                          XLP)/YLP        
         S1 = S1+S2                                                             
  164 CONTINUE                                                                  
      S = (F(1,1,1)+F(LSTOP,1,1)+F(1,1,NSTOP)+F(LSTOP,1,NSTOP)+                 
     1    F(1,MSTOP,1)+F(LSTOP,MSTOP,1)+F(1,MSTOP,NSTOP)+                       
     2                                   F(LSTOP,MSTOP,NSTOP))/(XLP*YLP)        
      DO 166 J=2,MSTPM1                                                         
         DO 165 I=2,LSTPM1                                                      
            S = S+F(I,J,1)+F(I,J,NSTOP)                                         
  165    CONTINUE                                                               
  166 CONTINUE                                                                  
      S2 = 0.                                                                   
      DO 167 I=2,LSTPM1                                                         
         S2 = S2+F(I,1,1)+F(I,1,NSTOP)+F(I,MSTOP,1)+F(I,MSTOP,NSTOP)            
  167 CONTINUE                                                                  
      S = S2/YLP+S                                                              
      S2 = 0.                                                                   
      DO 168 J=2,MSTPM1                                                         
         S2 = S2+F(1,J,1)+F(1,J,NSTOP)+F(LSTOP,J,1)+F(LSTOP,J,NSTOP)            
  168 CONTINUE                                                                  
      S = S2/XLP+S                                                              
      PERTRB = (S/ZLP+S1)/((FLOAT(LUNK+1)-XLP)*(FLOAT(MUNK+1)-YLP)*             
     1                                              (FLOAT(NUNK+1)-ZLP))        
      DO 171 I=1,LUNK                                                           
         DO 170 J=1,MUNK                                                        
            DO 169 K=1,NUNK                                                     
               F(I,J,K) = F(I,J,K)-PERTRB                                       
  169       CONTINUE                                                            
  170    CONTINUE                                                               
  171 CONTINUE                                                                  
  172 CONTINUE                                                                  
      NPEROD = 0                                                                
      IF (NBDCND .EQ. 0) GO TO 173                                              
      NPEROD = 1                                                                
      W(1) = 0.                                                                 
      W(IWW-1) = 0.                                                             
  173 CONTINUE                                                                  
      CALL POIS3D (LBDCND,LUNK,C1,MBDCND,MUNK,C2,NPEROD,NUNK,W,W(IWB),          
     1             W(IWC),LDIMF,MDIMF,F(LSTART,MSTART,NSTART),IR,W(IWW))        
C                                                                               
C     FILL IN SIDES FOR PERIODIC BOUNDARY CONDITIONS.                           
C                                                                               
      IF (LP .NE. 1) GO TO 180                                                  
      IF (MP .NE. 1) GO TO 175                                                  
      DO 174 K=NSTART,NSTOP                                                     
         F(1,MP1,K) = F(1,1,K)                                                  
  174 CONTINUE                                                                  
      MSTOP = MP1                                                               
  175 IF (NP .NE. 1) GO TO 177                                                  
      DO 176 J=MSTART,MSTOP                                                     
         F(1,J,NP1) = F(1,J,1)                                                  
  176 CONTINUE                                                                  
      NSTOP = NP1                                                               
  177 DO 179 J=MSTART,MSTOP                                                     
         DO 178 K=NSTART,NSTOP                                                  
            F(LP1,J,K) = F(1,J,K)                                               
  178    CONTINUE                                                               
  179 CONTINUE                                                                  
  180 CONTINUE                                                                  
      IF (MP .NE. 1) GO TO 185                                                  
      IF (NP .NE. 1) GO TO 182                                                  
      DO 181 I=LSTART,LSTOP                                                     
         F(I,1,NP1) = F(I,1,1)                                                  
  181 CONTINUE                                                                  
      NSTOP = NP1                                                               
  182 DO 184 I=LSTART,LSTOP                                                     
         DO 183 K=NSTART,NSTOP                                                  
            F(I,MP1,K) = F(I,1,K)                                               
  183    CONTINUE                                                               
  184 CONTINUE                                                                  
  185 CONTINUE                                                                  
      IF (NP .NE. 1) GO TO 188                                                  
      DO 187 I=LSTART,LSTOP                                                     
         DO 186 J=MSTART,MSTOP                                                  
            F(I,J,NP1) = F(I,J,1)                                               
  186    CONTINUE                                                               
  187 CONTINUE                                                                  
  188 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       

C FISHPAK18  FROM PORTLIB                                  03/11/81             
      SUBROUTINE POIS3D (LPEROD,L,C1,MPEROD,M,C2,NPEROD,N,A,B,C,LDIMF,          
     1                   MDIMF,F,IERROR,W)                                      
C                                                                               
C                                                                               
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
C     *                                                               *         
C     *                        F I S H P A K                          *         
C     *                                                               *         
C     *                                                               *         
C     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *         
C     *                                                               *         
C     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *         
C     *                                                               *         
C     *                  (VERSION 3.1 , OCTOBER 1980)                  *        
C     *                                                               *         
C     *                             BY                                *         
C     *                                                               *         
C     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *         
C     *                                                               *         
C     *                             OF                                *         
C     *                                                               *         
C     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *         
C     *                                                               *         
C     *                BOULDER, COLORADO  (80307)  U.S.A.             *         
C     *                                                               *         
C     *                   WHICH IS SPONSORED BY                       *         
C     *                                                               *         
C     *              THE NATIONAL SCIENCE FOUNDATION                  *         
C     *                                                               *         
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
C                                                                               
C                                                                               
C    * * * * * * * * *  PURPOSE    * * * * * * * * * * * * * * * * * *          
C                                                                               
C     SUBROUTINE POIS3D SOLVES THE LINEAR SYSTEM OF EQUATIONS                   
C                                                                               
C       C1*(X(I-1,J,K)-2.*X(I,J,K)+X(I+1,J,K))                                  
C     + C2*(X(I,J-1,K)-2.*X(I,J,K)+X(I,J+1,K))                                  
C     + A(K)*X(I,J,K-1)+B(K)*X(I,J,K)+C(K)*X(I,J,K+1) = F(I,J,K)                
C                                                                               
C     FOR  I=1,2,...,L , J=1,2,...,M , AND K=1,2,...,N .                        
C                                                                               
C     THE INDICES K-1 AND K+1 ARE EVALUATED MODULO N, I.E.                      
C     X(I,J,0) = X(I,J,N) AND X(I,J,N+1) = X(I,J,1). THE UNKNOWNS               
C     X(0,J,K), X(L+1,J,K), X(I,0,K), AND X(I,M+1,K) ARE ASSUMED TO TAKE        
C     ON CERTAIN PRESCRIBED VALUES DESCRIBED BELOW.                             
C                                                                               
C    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *          
C                                                                               
C                                                                               
C    * * * * * * * *    PARAMETER DESCRIPTION     * * * * * * * * * *           
C                                                                               
C                                                                               
C            * * * * * *   ON INPUT    * * * * * *                              
C                                                                               
C     LPEROD   INDICATES THE VALUES THAT X(0,J,K) AND X(L+1,J,K) ARE            
C              ASSUMED TO HAVE.                                                 
C                                                                               
C              = 0  IF X(0,J,K) = X(L,J,K) AND X(L+1,J,K) = X(1,J,K).           
C              = 1  IF X(0,J,K) = X(L+1,J,K) = 0.                               
C              = 2  IF X(0,J,K) = 0  AND X(L+1,J,K) = X(L-1,J,K).               
C              = 3  IF X(0,J,K) = X(2,J,K) AND X(L+1,J,K) = X(L-1,J,K).         
C              = 4  IF X(0,J,K) = X(2,J,K) AND X(L+1,J,K) = 0.                  
C                                                                               
C     L        THE NUMBER OF UNKNOWNS IN THE I-DIRECTION. L MUST BE AT          
C              LEAST 3.                                                         
C                                                                               
C     C1       THE REAL CONSTANT THAT APPEARS IN THE ABOVE EQUATION.            
C                                                                               
C     MPEROD   INDICATES THE VALUES THAT X(I,0,K) AND X(I,M+1,K) ARE            
C              ASSUMED TO HAVE.                                                 
C                                                                               
C              = 0  IF X(I,0,K) = X(I,M,K) AND X(I,M+1,K) = X(I,1,K).           
C              = 1  IF X(I,0,K) = X(I,M+1,K) = 0.                               
C              = 2  IF X(I,0,K) = 0 AND X(I,M+1,K) = X(I,M-1,K).                
C              = 3  IF X(I,0,K) = X(I,2,K) AND X(I,M+1,K) = X(I,M-1,K).         
C              = 4  IF X(I,0,K) = X(I,2,K) AND X(I,M+1,K) = 0.                  
C                                                                               
C     M        THE NUMBER OF UNKNOWNS IN THE J-DIRECTION. M MUST BE AT          
C              LEAST 3.                                                         
C                                                                               
C     C2       THE REAL CONSTANT WHICH APPEARS IN THE ABOVE EQUATION.           
C                                                                               
C     NPEROD   = 0  IF A(1) AND C(N) ARE NOT ZERO.                              
C              = 1  IF A(1) = C(N) = 0.                                         
C                                                                               
C     N        THE NUMBER OF UNKNOWNS IN THE K-DIRECTION. N MUST BE AT          
C              LEAST 3.                                                         
C                                                                               
C                                                                               
C     A,B,C    ONE-DIMENSIONAL ARRAYS OF LENGTH N THAT SPECIFY THE              
C              COEFFICIENTS IN THE LINEAR EQUATIONS GIVEN ABOVE.                
C                                                                               
C              IF NPEROD = 0 THE ARRAY ELEMENTS MUST NOT DEPEND UPON THE        
C              INDEX K, BUT MUST BE CONSTANT.  SPECIFICALLY,THE                 
C              SUBROUTINE CHECKS THE FOLLOWING CONDITION                        
C                                                                               
C                          A(K) = C(1)                                          
C                          C(K) = C(1)                                          
C                          B(K) = B(1)                                          
C                                                                               
C                  FOR K=1,2,...,N.                                             
C                                                                               
C     LDIMF    THE ROW (OR FIRST) DIMENSION OF THE THREE-DIMENSIONAL            
C              ARRAY F AS IT APPEARS IN THE PROGRAM CALLING POIS3D.             
C              THIS PARAMETER IS USED TO SPECIFY THE VARIABLE DIMENSION         
C              OF F.  LDIMF MUST BE AT LEAST L.                                 
C                                                                               
C     MDIMF    THE COLUMN (OR SECOND) DIMENSION OF THE THREE-DIMENSIONAL        
C              ARRAY F AS IT APPEARS IN THE PROGRAM CALLING POIS3D.             
C              THIS PARAMETER IS USED TO SPECIFY THE VARIABLE DIMENSION         
C              OF F.  MDIMF MUST BE AT LEAST M.                                 
C                                                                               
C     F        A THREE-DIMENSIONAL ARRAY THAT SPECIFIES THE VALUES OF           
C              THE RIGHT SIDE OF THE LINEAR SYSTEM OF EQUATIONS GIVEN           
C              ABOVE.  F MUST BE DIMENSIONED AT LEAST L X M X N.                
C                                                                               
C     W        A ONE-DIMENSIONAL ARRAY THAT MUST BE PROVIDED BY THE             
C              USER FOR WORK SPACE.  THE LENGTH OF W MUST BE AT LEAST           
C              30 + L + M + 2*N + MAX(L,M,N) +                                  
C              7*(INT((L+1)/2) + INT((M+1)/2)).                                 
C                                                                               
C                                                                               
C            * * * * * *   ON OUTPUT   * * * * * *                              
C                                                                               
C     F        CONTAINS THE SOLUTION X.                                         
C                                                                               
C     IERROR   AN ERROR FLAG THAT INDICATES INVALID INPUT PARAMETERS.           
C              EXCEPT FOR NUMBER ZERO, A SOLUTION IS NOT ATTEMPTED.             
C              = 0  NO ERROR                                                    
C              = 1  IF LPEROD .LT. 0 OR .GT. 4                                  
C              = 2  IF L .LT. 3                                                 
C              = 3  IF MPEROD .LT. 0 OR .GT. 4                                  
C              = 4  IF M .LT. 3                                                 
C              = 5  IF NPEROD .LT. 0 OR .GT. 1                                  
C              = 6  IF N .LT. 3                                                 
C              = 7  IF LDIMF .LT. L                                             
C              = 8  IF MDIMF .LT. M                                             
C              = 9  IF A(K) .NE. C(1) OR C(K) .NE. C(1) OR B(I) .NE.B(1)        
C                      FOR SOME K=1,2,...,N.                                    
C              = 10 IF NPEROD = 1 AND A(1) .NE. 0 OR C(N) .NE. 0                
C                                                                               
C              SINCE THIS IS THE ONLY MEANS OF INDICATING A POSSIBLY            
C              INCORRECT CALL TO POIS3D, THE USER SHOULD TEST IERROR            
C              AFTER THE CALL.                                                  
C                                                                               
C                                                                               
C    * * * * * * *   PROGRAM SPECIFICATIONS    * * * * * * * * * * * *          
C                                                                               
C     DIMENSION OF   A(N),B(N),C(N),F(LDIMF,MDIMF,N),                           
C     ARGUMENTS      W(SEE ARGUMENT LIST)                                       
C                                                                               
C     LATEST         DECEMBER 1, 1978                                           
C     REVISION                                                                  
C                                                                               
C     SUBPROGRAMS    POIS3D,POS3D1,TRID,RFFTI,RFFTF,RFFTF1,RFFTB,               
C     REQUIRED       RFFTB1,COSTI,COST,SINTI,SINT,COSQI,COSQF,COSQF1            
C                    COSQB,COSQB1,SINQI,SINQF,SINQB,CFFTI,CFFTI1,               
C                    CFFTB,CFFTB1,PASSB2,PASSB3,PASSB4,PASSB,CFFTF,             
C                    CFFTF1,PASSF1,PASSF2,PASSF3,PASSF4,PASSF,PIMACH,           
C                                                                               
C     SPECIAL        NONE                                                       
C     CONDITIONS                                                                
C                                                                               
C     COMMON         VALUE                                                      
C     BLOCKS                                                                    
C                                                                               
C     I/O            NONE                                                       
C                                                                               
C     PRECISION      SINGLE                                                     
C                                                                               
C     SPECIALIST     ROLAND SWEET                                               
C                                                                               
C     LANGUAGE       FORTRAN                                                    
C                                                                               
C     HISTORY        WRITTEN BY ROLAND SWEET AT NCAR IN JULY,1977               
C                                                                               
C     ALGORITHM      THIS SUBROUTINE SOLVES THREE-DIMENSIONAL BLOCK             
C                    TRIDIAGONAL LINEAR SYSTEMS ARISING FROM FINITE             
C                    DIFFERENCE APPROXIMATIONS TO THREE-DIMENSIONAL             
C                    POISSON EQUATIONS USING THE FOURIER TRANSFORM              
C                    PACKAGE SCLRFFTPAK WRITTEN BY PAUL SWARZTRAUBER.           
C                                                                               
C     SPACE          6561(DECIMAL) = 14641(OCTAL) LOCATIONS ON THE              
C     REQUIRED       NCAR CONTROL DATA 7600                                     
C                                                                               
C     TIMING AND        THE EXECUTION TIME T ON THE NCAR CONTROL DATA           
C     ACCURACY       7600 FOR SUBROUTINE POIS3D IS ROUGHLY PROPORTIONAL         
C                    TO L*M*N*(LOG2(L)+LOG2(M)+5), BUT ALSO DEPENDS ON          
C                    INPUT PARAMETERS LPEROD AND MPEROD.  SOME TYPICAL          
C                    VALUES ARE LISTED IN THE TABLE BELOW WHEN NPEROD=0.        
C                       TO MEASURE THE ACCURACY OF THE ALGORITHM A              
C                    UNIFORM RANDOM NUMBER GENERATOR WAS USED TO CREATE         
C                    A SOLUTION ARRAY X FOR THE SYSTEM GIVEN IN THE             
C                    "PURPOSE" WITH                                             
C                                                                               
C                       A(K) = C(K) = -0.5*B(K) = 1,       K=1,2,...,N          
C                                                                               
C                    AND, WHEN NPEROD = 1                                       
C                                                                               
C                       A(1) = C(N) = 0                                         
C                       A(N) = C(1) = 2.                                        
C                                                                               
C                    THE SOLUTION X WAS SUBSTITUTED INTO THE GIVEN SYS-         
C                    TEM AND, USING DOUBLE PRECISION, A RIGHT SIDE Y WAS        
C                    COMPUTED.  USING THIS ARRAY Y SUBROUTINE POIS WAS          
C                    CALLED TO PRODUCE AN APPROXIMATE SOLUTION Z.  THEN         
C                    THE RELATIVE ERROR, DEFINED AS                             
C                                                                               
C                    E = MAX(ABS(Z(I,J,K)-X(I,J,K)))/MAX(ABS(X(I,J,K)))         
C                                                                               
C                    WHERE THE TWO MAXIMA ARE TAKEN OVER I=1,2,...,L,           
C                    J=1,2,...,M AND K=1,2,...,N, WAS COMPUTED.  THE            
C                    VALUE OF E IS GIVEN IN THE TABLE BELOW FOR SOME            
C                    TYPICAL VALUES OF L,M AND N.                               
C                                                                               
C                                                                               
C                       L(=M=N)   LPEROD    MPEROD    T(MSECS)    E             
C                       ------    ------    ------    --------  ------          
C                                                                               
C                         16        0         0         272     1.E-13          
C                         15        1         1         287     4.E-13          
C                         17        3         3         338     2.E-13          
C                         32        0         0        1755     2.E-13          
C                         31        1         1        1894     2.E-12          
C                         33        3         3        2042     7.E-13          
C                                                                               
C                                                                               
C     PORTABILITY    AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN.             
C                    THE MACHINE DEPENDENT CONSTANT PI IS DEFINED IN            
C                    FUNCTION PIMACH.                                           
C                                                                               
C     REQUIRED       COS,SIN,ATAN                                               
C     RESIDENT                                                                  
C     ROUTINES                                                                  
C                                                                               
C     REFERENCE      NONE                                                       
C                                                                               
C    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *          
C                                                                               
      DIMENSION       A(1)       ,B(1)       ,C(1)       ,                      
     1                F(LDIMF,MDIMF,1)       ,W(1)       ,SAVE(6)               
      LP = LPEROD+1                                                             
      MP = MPEROD+1                                                             
      NP = NPEROD+1                                                             
C                                                                               
C     CHECK FOR INVALID INPUT.                                                  
C                                                                               
      IERROR = 0                                                                
      IF (LP.LT.1 .OR. LP.GT.5) IERROR = 1                                      
      IF (L .LT. 3) IERROR = 2                                                  
      IF (MP.LT.1 .OR. MP.GT.5) IERROR = 3                                      
      IF (M .LT. 3) IERROR = 4                                                  
      IF (NP.LT.1 .OR. NP.GT.2) IERROR = 5                                      
      IF (N .LT. 3) IERROR = 6                                                  
      IF (LDIMF .LT. L) IERROR = 7                                              
      IF (MDIMF .LT. M) IERROR = 8                                              
      IF (NP .NE. 1) GO TO 103                                                  
      DO 101 K=1,N                                                              
         IF (A(K) .NE. C(1)) GO TO 102                                          
         IF (C(K) .NE. C(1)) GO TO 102                                          
         IF (B(K) .NE. B(1)) GO TO 102                                          
  101 CONTINUE                                                                  
      GO TO 104                                                                 
  102 IERROR = 9                                                                
  103 IF (NPEROD.EQ.1 .AND. (A(1).NE.0. .OR. C(N).NE.0.)) IERROR = 10           
  104 IF (IERROR .NE. 0) GO TO 122                                              
      IWYRT = L+1                                                               
      IWT = IWYRT+M                                                             
      IWD = IWT+MAX0(L,M,N)+1                                                   
      IWBB = IWD+N                                                              
      IWX = IWBB+N                                                              
      IWY = IWX+7*((L+1)/2)+15                                                  
      GO TO (105,114),NP                                                        
C                                                                               
C     REORDER UNKNOWNS WHEN NPEROD = 0.                                         
C                                                                               
  105 NH = (N+1)/2                                                              
      NHM1 = NH-1                                                               
      NODD = 1                                                                  
      IF (2*NH .EQ. N) NODD = 2                                                 
      DO 111 I=1,L                                                              
         DO 110 J=1,M                                                           
            DO 106 K=1,NHM1                                                     
               NHPK = NH+K                                                      
               NHMK = NH-K                                                      
               W(K) = F(I,J,NHMK)-F(I,J,NHPK)                                   
               W(NHPK) = F(I,J,NHMK)+F(I,J,NHPK)                                
  106       CONTINUE                                                            
            W(NH) = 2.*F(I,J,NH)                                                
            GO TO (108,107),NODD                                                
  107       W(N) = 2.*F(I,J,N)                                                  
  108       DO 109 K=1,N                                                        
               F(I,J,K) = W(K)                                                  
  109       CONTINUE                                                            
  110    CONTINUE                                                               
  111 CONTINUE                                                                  
      SAVE(1) = C(NHM1)                                                         
      SAVE(2) = A(NH)                                                           
      SAVE(3) = C(NH)                                                           
      SAVE(4) = B(NHM1)                                                         
      SAVE(5) = B(N)                                                            
      SAVE(6) = A(N)                                                            
      C(NHM1) = 0.                                                              
      A(NH) = 0.                                                                
      C(NH) = 2.*C(NH)                                                          
      GO TO (112,113),NODD                                                      
  112 B(NHM1) = B(NHM1)-A(NH-1)                                                 
      B(N) = B(N)+A(N)                                                          
      GO TO 114                                                                 
  113 A(N) = C(NH)                                                              
  114 CONTINUE                                                                  
      CALL POS3D1 (LP,L,MP,M,N,A,B,C,LDIMF,MDIMF,F,W,W(IWYRT),W(IWT),           
     1             W(IWD),W(IWX),W(IWY),C1,C2,W(IWBB))                          
      GO TO (115,122),NP                                                        
  115 DO 121 I=1,L                                                              
         DO 120 J=1,M                                                           
            DO 116 K=1,NHM1                                                     
               NHMK = NH-K                                                      
               NHPK = NH+K                                                      
               W(NHMK) = .5*(F(I,J,NHPK)+F(I,J,K))                              
               W(NHPK) = .5*(F(I,J,NHPK)-F(I,J,K))                              
  116       CONTINUE                                                            
            W(NH) = .5*F(I,J,NH)                                                
            GO TO (118,117),NODD                                                
  117       W(N) = .5*F(I,J,N)                                                  
  118       DO 119 K=1,N                                                        
               F(I,J,K) = W(K)                                                  
  119       CONTINUE                                                            
  120    CONTINUE                                                               
  121 CONTINUE                                                                  
      C(NHM1) = SAVE(1)                                                         
      A(NH) = SAVE(2)                                                           
      C(NH) = SAVE(3)                                                           
      B(NHM1) = SAVE(4)                                                         
      B(N) = SAVE(5)                                                            
      A(N) = SAVE(6)                                                            
  122 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE POS3D1 (LP,L,MP,M,N,A,B,C,LDIMF,MDIMF,F,XRT,YRT,T,D,           
     1                   WX,WY,C1,C2,BB)                                        
      DIMENSION       A(1)       ,B(1)       ,C(1)       ,                      
     1                F(LDIMF,MDIMF,1)       ,XRT(1)     ,YRT(1)     ,          
     2                T(1)       ,D(1)       ,WX(1)      ,WY(1)      ,          
     3                BB(1)                                                     
      PI = PIMACH(DUM)                                                          
      LR = L                                                                    
      MR = M                                                                    
      NR = N                                                                    
C                                                                               
C     GENERATE TRANSFORM ROOTS                                                  
C                                                                               
      LRDEL = ((LP-1)*(LP-3)*(LP-5))/3                                          
      SCALX = LR+LRDEL                                                          
      DX = PI/(2.*SCALX)                                                        
      GO TO (108,103,101,102,101),LP                                            
  101 DI = 0.5                                                                  
      SCALX = 2.*SCALX                                                          
      GO TO 104                                                                 
  102 DI = 1.0                                                                  
      GO TO 104                                                                 
  103 DI = 0.0                                                                  
  104 DO 105 I=1,LR                                                             
         XRT(I) = -4.*C1*(SIN((FLOAT(I)-DI)*DX))**2                             
  105 CONTINUE                                                                  
      SCALX = 2.*SCALX                                                          
      GO TO (112,106,110,107,111),LP                                            
  106 CALL SINTI (LR,WX)                                                        
      GO TO 112                                                                 
  107 CALL COSTI (LR,WX)                                                        
      GO TO 112                                                                 
  108 XRT(1) = 0.                                                               
      XRT(LR) = -4.*C1                                                          
      DO 109 I=3,LR,2                                                           
         XRT(I-1) = -4.*C1*(SIN(FLOAT((I-1))*DX))**2                            
         XRT(I) = XRT(I-1)                                                      
  109 CONTINUE                                                                  
      CALL RFFTI (LR,WX)                                                        
      GO TO 112                                                                 
  110 CALL SINQI (LR,WX)                                                        
      GO TO 112                                                                 
  111 CALL COSQI (LR,WX)                                                        
  112 CONTINUE                                                                  
      MRDEL = ((MP-1)*(MP-3)*(MP-5))/3                                          
      SCALY = MR+MRDEL                                                          
      DY = PI/(2.*SCALY)                                                        
      GO TO (120,115,113,114,113),MP                                            
  113 DJ = 0.5                                                                  
      SCALY = 2.*SCALY                                                          
      GO TO 116                                                                 
  114 DJ = 1.0                                                                  
      GO TO 116                                                                 
  115 DJ = 0.0                                                                  
  116 DO 117 J=1,MR                                                             
         YRT(J) = -4.*C2*(SIN((FLOAT(J)-DJ)*DY))**2                             
  117 CONTINUE                                                                  
      SCALY = 2.*SCALY                                                          
      GO TO (124,118,122,119,123),MP                                            
  118 CALL SINTI (MR,WY)                                                        
      GO TO 124                                                                 
  119 CALL COSTI (MR,WY)                                                        
      GO TO 124                                                                 
  120 YRT(1) = 0.                                                               
      YRT(MR) = -4.*C2                                                          
      DO 121 J=3,MR,2                                                           
         YRT(J-1) = -4.*C2*(SIN(FLOAT((J-1))*DY))**2                            
         YRT(J) = YRT(J-1)                                                      
  121 CONTINUE                                                                  
      CALL RFFTI (MR,WY)                                                        
      GO TO 124                                                                 
  122 CALL SINQI (MR,WY)                                                        
      GO TO 124                                                                 
  123 CALL COSQI (MR,WY)                                                        
  124 CONTINUE                                                                  
      IFWRD = 1                                                                 
      IS = 1                                                                    
  125 CONTINUE                                                                  
C                                                                               
C     TRANSFORM X                                                               
C                                                                               
      DO 141 J=1,MR                                                             
         DO 140 K=1,NR                                                          
            DO 126 I=1,LR                                                       
               T(I) = F(I,J,K)                                                  
  126       CONTINUE                                                            
            GO TO (127,130,131,134,135),LP                                      
  127       GO TO (128,129),IFWRD                                               
  128       CALL RFFTF (LR,T,WX)                                                
            GO TO 138                                                           
  129       CALL RFFTB (LR,T,WX)                                                
            GO TO 138                                                           
  130       CALL SINT (LR,T,WX)                                                 
            GO TO 138                                                           
  131       GO TO (132,133),IFWRD                                               
  132       CALL SINQF (LR,T,WX)                                                
            GO TO 138                                                           
  133       CALL SINQB (LR,T,WX)                                                
            GO TO 138                                                           
  134       CALL COST (LR,T,WX)                                                 
            GO TO 138                                                           
  135       GO TO (136,137),IFWRD                                               
  136       CALL COSQF (LR,T,WX)                                                
            GO TO 138                                                           
  137       CALL COSQB (LR,T,WX)                                                
  138       CONTINUE                                                            
            DO 139 I=1,LR                                                       
               F(I,J,K) = T(I)                                                  
  139       CONTINUE                                                            
  140    CONTINUE                                                               
  141 CONTINUE                                                                  
      GO TO (142,164),IFWRD                                                     
C                                                                               
C     TRANSFORM Y                                                               
C                                                                               
  142 CONTINUE                                                                  
      DO 158 I=1,LR                                                             
         DO 157 K=1,NR                                                          
            DO 143 J=1,MR                                                       
               T(J) = F(I,J,K)                                                  
  143       CONTINUE                                                            
            GO TO (144,147,148,151,152),MP                                      
  144       GO TO (145,146),IFWRD                                               
  145       CALL RFFTF (MR,T,WY)                                                
            GO TO 155                                                           
  146       CALL RFFTB (MR,T,WY)                                                
            GO TO 155                                                           
  147       CALL SINT (MR,T,WY)                                                 
            GO TO 155                                                           
  148       GO TO (149,150),IFWRD                                               
  149       CALL SINQF (MR,T,WY)                                                
            GO TO 155                                                           
  150       CALL SINQB (MR,T,WY)                                                
            GO TO 155                                                           
  151       CALL COST (MR,T,WY)                                                 
            GO TO 155                                                           
  152       GO TO (153,154),IFWRD                                               
  153       CALL COSQF (MR,T,WY)                                                
            GO TO 155                                                           
  154       CALL COSQB (MR,T,WY)                                                
  155       CONTINUE                                                            
            DO 156 J=1,MR                                                       
               F(I,J,K) = T(J)                                                  
  156       CONTINUE                                                            
  157    CONTINUE                                                               
  158 CONTINUE                                                                  
      GO TO (159,125),IFWRD                                                     
  159 CONTINUE                                                                  
C                                                                               
C     SOLVE TRIDIAGONAL SYSTEMS IN Z                                            
C                                                                               
      DO 163 I=1,LR                                                             
         DO 162 J=1,MR                                                          
            DO 160 K=1,NR                                                       
               BB(K) = B(K)+XRT(I)+YRT(J)                                       
               T(K) = F(I,J,K)                                                  
  160       CONTINUE                                                            
            CALL FTRID (NR,A,BB,C,T,D)                                           
            DO 161 K=1,NR                                                       
               F(I,J,K) = T(K)                                                  
  161       CONTINUE                                                            
  162    CONTINUE                                                               
  163 CONTINUE                                                                  
      IFWRD = 2                                                                 
      IS = -1                                                                   
      GO TO 142                                                                 
  164 CONTINUE                                                                  
      DO 167 I=1,LR                                                             
         DO 166 J=1,MR                                                          
            DO 165 K=1,NR                                                       
               F(I,J,K) = F(I,J,K)/(SCALX*SCALY)                                
  165       CONTINUE                                                            
  166    CONTINUE                                                               
  167 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE FTRID (MR,A,B,C,Y,D)                                            
      DIMENSION       A(1)       ,B(1)       ,C(1)       ,Y(1)       ,          
     1                D(1)                                                      
      M = MR                                                                    
      MM1 = M-1                                                                 
      Z = 1./B(1)                                                               
      D(1) = C(1)*Z                                                             
      Y(1) = Y(1)*Z                                                             
      DO 101 I=2,MM1                                                            
         Z = 1./(B(I)-A(I)*D(I-1))                                              
         D(I) = C(I)*Z                                                          
         Y(I) = (Y(I)-A(I)*Y(I-1))*Z                                            
  101 CONTINUE                                                                  
      Z = B(M)-A(M)*D(MM1)                                                      
      IF (Z .NE. 0.) GO TO 102                                                  
      Y(M) = 0.                                                                 
      GO TO 103                                                                 
  102 Y(M) = (Y(M)-A(M)*Y(MM1))/Z                                               
  103 CONTINUE                                                                  
      DO 104 IP=1,MM1                                                           
         I = M-IP                                                               
         Y(I) = Y(I)-D(I)*Y(I+1)                                                
  104 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       

 
C FISHPAK21  FROM PORTLIB                                  03/11/81             
      SUBROUTINE EZFFTF (N,R,AZERO,A,B,WSAVE)                                   
C                                                                               
C                       VERSION 3  JUNE 1979                                    
C                                                                               
      DIMENSION       R(1)       ,A(1)       ,B(1)       ,WSAVE(1)              
      IF (N-2) 101,102,103                                                      
  101 AZERO = R(1)                                                              
      RETURN                                                                    
  102 AZERO = .5*(R(1)+R(2))                                                    
      A(1) = .5*(R(1)-R(2))                                                     
      RETURN                                                                    
  103 DO 104 I=1,N                                                              
         WSAVE(I) = R(I)                                                        
  104 CONTINUE                                                                  
      CALL RFFTF (N,WSAVE,WSAVE(N+1))                                           
      CF = 2./FLOAT(N)                                                          
      CFM = -CF                                                                 
      AZERO = .5*CF*WSAVE(1)                                                    
      NS2 = (N+1)/2                                                             
      NS2M = NS2-1                                                              
      DO 105 I=1,NS2M                                                           
         A(I) = CF*WSAVE(2*I)                                                   
         B(I) = CFM*WSAVE(2*I+1)                                                
  105 CONTINUE                                                                  
      IF (MOD(N,2) .EQ. 0) A(NS2) = .5*CF*WSAVE(N)                              
      RETURN                                                                    
      END                                                                       
      SUBROUTINE EZFFTB (N,R,AZERO,A,B,WSAVE)                                   
      DIMENSION       R(1)       ,A(1)       ,B(1)       ,WSAVE(1)              
      IF (N-2) 101,102,103                                                      
  101 R(1) = AZERO                                                              
      RETURN                                                                    
  102 R(1) = AZERO+A(1)                                                         
      R(2) = AZERO-A(1)                                                         
      RETURN                                                                    
  103 NS2 = (N-1)/2                                                             
      DO 104 I=1,NS2                                                            
         R(2*I) = .5*A(I)                                                       
         R(2*I+1) = -.5*B(I)                                                    
  104 CONTINUE                                                                  
      R(1) = AZERO                                                              
      IF (MOD(N,2) .EQ. 0) R(N) = A(NS2+1)                                      
      CALL RFFTB (N,R,WSAVE(N+1))                                               
      RETURN                                                                    
      END                                                                       
      SUBROUTINE EZFFTI (N,WSAVE)                                               
      DIMENSION       WSAVE(1)                                                  
      IF (N .EQ. 1) RETURN                                                      
      CALL EZFFT1 (N,WSAVE(2*N+1),WSAVE(3*N+1))                                 
      RETURN                                                                    
      END                                                                       
      SUBROUTINE EZFFT1 (N,WA,IFAC)                                             
      DIMENSION       WA(1)      ,IFAC(1)    ,NTRYH(4)                          
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/                         
     1    ,TPI/6.28318530717959/                                                
      NL = N                                                                    
      NF = 0                                                                    
      J = 0                                                                     
  101 J = J+1                                                                   
      IF (J-4) 102,102,103                                                      
  102 NTRY = NTRYH(J)                                                           
      GO TO 104                                                                 
  103 NTRY = NTRY+2                                                             
  104 NQ = NL/NTRY                                                              
      NR = NL-NTRY*NQ                                                           
      IF (NR) 101,105,101                                                       
  105 NF = NF+1                                                                 
      IFAC(NF+2) = NTRY                                                         
      NL = NQ                                                                   
      IF (NTRY .NE. 2) GO TO 107                                                
      IF (NF .EQ. 1) GO TO 107                                                  
      DO 106 I=2,NF                                                             
         IB = NF-I+2                                                            
         IFAC(IB+2) = IFAC(IB+1)                                                
  106 CONTINUE                                                                  
      IFAC(3) = 2                                                               
  107 IF (NL .NE. 1) GO TO 104                                                  
      IFAC(1) = N                                                               
      IFAC(2) = NF                                                              
      ARGH = TPI/FLOAT(N)                                                       
      IS = 0                                                                    
      NFM1 = NF-1                                                               
      L1 = 1                                                                    
      IF (NFM1 .EQ. 0) RETURN                                                   
      DO 111 K1=1,NFM1                                                          
         IP = IFAC(K1+2)                                                        
         L2 = L1*IP                                                             
         IDO = N/L2                                                             
         IPM = IP-1                                                             
         ARG1 = FLOAT(L1)*ARGH                                                  
         CH1 = 1.                                                               
         SH1 = 0.                                                               
         DCH1 = COS(ARG1)                                                       
         DSH1 = SIN(ARG1)                                                       
         DO 110 J=1,IPM                                                         
            CH1H = DCH1*CH1-DSH1*SH1                                            
            SH1 = DCH1*SH1+DSH1*CH1                                             
            CH1 = CH1H                                                          
            I = IS+2                                                            
            WA(I-1) = CH1                                                       
            WA(I) = SH1                                                         
            IF (IDO .LT. 5) GO TO 109                                           
            DO 108 II=5,IDO,2                                                   
               I = I+2                                                          
               WA(I-1) = CH1*WA(I-3)-SH1*WA(I-2)                                
               WA(I) = CH1*WA(I-2)+SH1*WA(I-3)                                  
  108       CONTINUE                                                            
  109       IS = IS+IDO                                                         
  110    CONTINUE                                                               
         L1 = L2                                                                
  111 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE COSTI (N,WSAVE)                                                
      DIMENSION       WSAVE(1)                                                  
      DATA PI /3.14159265358979/                                                
      IF (N .LE. 3) RETURN                                                      
      NM1 = N-1                                                                 
      NP1 = N+1                                                                 
      NS2 = N/2                                                                 
      DT = PI/FLOAT(NM1)                                                        
      FK = 0.                                                                   
      DO 101 K=2,NS2                                                            
         KC = NP1-K                                                             
         FK = FK+1.                                                             
         WSAVE(K) = 2.*SIN(FK*DT)                                               
         WSAVE(KC) = 2.*COS(FK*DT)                                              
  101 CONTINUE                                                                  
      CALL RFFTI (NM1,WSAVE(N+1))                                               
      RETURN                                                                    
      END                                                                       
      SUBROUTINE COST (N,X,WSAVE)                                               
      DIMENSION       X(1)       ,WSAVE(1)                                      
      NM1 = N-1                                                                 
      NP1 = N+1                                                                 
      NS2 = N/2                                                                 
      IF (N-2) 106,101,102                                                      
  101 X1H = X(1)+X(2)                                                           
      X(2) = X(1)-X(2)                                                          
      X(1) = X1H                                                                
      RETURN                                                                    
  102 IF (N .GT. 3) GO TO 103                                                   
      X1P3 = X(1)+X(3)                                                          
      TX2 = X(2)+X(2)                                                           
      X(2) = X(1)-X(3)                                                          
      X(1) = X1P3+TX2                                                           
      X(3) = X1P3-TX2                                                           
      RETURN                                                                    
  103 C1 = X(1)-X(N)                                                            
      X(1) = X(1)+X(N)                                                          
      DO 104 K=2,NS2                                                            
         KC = NP1-K                                                             
         T1 = X(K)+X(KC)                                                        
         T2 = X(K)-X(KC)                                                        
         C1 = C1+WSAVE(KC)*T2                                                   
         T2 = WSAVE(K)*T2                                                       
         X(K) = T1-T2                                                           
         X(KC) = T1+T2                                                          
  104 CONTINUE                                                                  
      MODN = MOD(N,2)                                                           
      IF (MODN .NE. 0) X(NS2+1) = X(NS2+1)+X(NS2+1)                             
      CALL RFFTF (NM1,X,WSAVE(N+1))                                             
      XIM2 = X(2)                                                               
      X(2) = C1                                                                 
      DO 105 I=4,N,2                                                            
         XI = X(I)                                                              
         X(I) = X(I-2)-X(I-1)                                                   
         X(I-1) = XIM2                                                          
         XIM2 = XI                                                              
  105 CONTINUE                                                                  
      IF (MODN .NE. 0) X(N) = XIM2                                              
  106 RETURN                                                                    
      END                                                                       
      SUBROUTINE SINTI (N,WSAVE)                                                
      DIMENSION       WSAVE(1)                                                  
      DATA PI /3.14159265358979/                                                
      IF (N .LE. 1) RETURN                                                      
      NP1 = N+1                                                                 
      NS2 = N/2                                                                 
      DT = PI/FLOAT(NP1)                                                        
      FK = 0.                                                                   
      DO 101 K=1,NS2                                                            
         FK = FK+1.                                                             
         WSAVE(K) = 2.*SIN(FK*DT)                                               
  101 CONTINUE                                                                  
      CALL RFFTI (NP1,WSAVE(NS2+1))                                             
      RETURN                                                                    
      END                                                                       
      SUBROUTINE SINT (N,X,WSAVE)                                               
      DIMENSION       X(1)       ,WSAVE(1)                                      
      DATA SQRT3 /1.73205080756888/                                             
      IF (N-2) 101,102,103                                                      
  101 X(1) = X(1)+X(1)                                                          
      RETURN                                                                    
  102 XH = SQRT3*(X(1)+X(2))                                                    
      X(2) = SQRT3*(X(1)-X(2))                                                  
      X(1) = XH                                                                 
      RETURN                                                                    
  103 NP1 = N+1                                                                 
      NS2 = N/2                                                                 
      X1 = X(1)                                                                 
      X(1) = 0.                                                                 
      DO 104 K=1,NS2                                                            
         KC = NP1-K                                                             
         T1 = X1-X(KC)                                                          
         T2 = WSAVE(K)*(X1+X(KC))                                               
         X1 = X(K+1)                                                            
         X(K+1) = T1+T2                                                         
         X(KC+1) = T2-T1                                                        
  104 CONTINUE                                                                  
      MODN = MOD(N,2)                                                           
      IF (MODN .NE. 0) X(NS2+2) = 4.*X1                                         
      CALL RFFTF (NP1,X,WSAVE(NS2+1))                                           
      X(1) = .5*X(1)                                                            
      DO 105 I=3,N,2                                                            
         XIM1 = X(I-1)                                                          
         X(I-1) = -X(I)                                                         
         X(I) = X(I-2)+XIM1                                                     
  105 CONTINUE                                                                  
      IF (MODN .NE. 0) RETURN                                                   
      X(N) = -X(N+1)                                                            
      RETURN                                                                    
      END                                                                       
      SUBROUTINE COSQI (N,WSAVE)                                                
      DIMENSION       WSAVE(1)                                                  
      DATA PIH /1.57079632679491/                                               
      DT = PIH/FLOAT(N)                                                         
      FK = 0.                                                                   
      DO 101 K=1,N                                                              
         FK = FK+1.                                                             
         WSAVE(K) = COS(FK*DT)                                                  
  101 CONTINUE                                                                  
      CALL RFFTI (N,WSAVE(N+1))                                                 
      RETURN                                                                    
      END                                                                       
      SUBROUTINE COSQF (N,X,WSAVE)                                              
      DIMENSION       X(1)       ,WSAVE(1)                                      
      DATA SQRT2 /1.4142135623731/                                              
      IF (N-2) 102,101,103                                                      
  101 TSQX = SQRT2*X(2)                                                         
      X(2) = X(1)-TSQX                                                          
      X(1) = X(1)+TSQX                                                          
  102 RETURN                                                                    
  103 CALL COSQF1 (N,X,WSAVE,WSAVE(N+1))                                        
      RETURN                                                                    
      END                                                                       
      SUBROUTINE COSQF1 (N,X,W,XH)                                              
      DIMENSION       X(1)       ,W(1)       ,XH(1)                             
      NS2 = (N+1)/2                                                             
      NP2 = N+2                                                                 
      DO 101 K=2,NS2                                                            
         KC = NP2-K                                                             
         XH(K) = X(K)+X(KC)                                                     
         XH(KC) = X(K)-X(KC)                                                    
  101 CONTINUE                                                                  
      MODN = MOD(N,2)                                                           
      IF (MODN .EQ. 0) XH(NS2+1) = X(NS2+1)+X(NS2+1)                            
      DO 102 K=2,NS2                                                            
         KC = NP2-K                                                             
         X(K) = W(K-1)*XH(KC)+W(KC-1)*XH(K)                                     
         X(KC) = W(K-1)*XH(K)-W(KC-1)*XH(KC)                                    
  102 CONTINUE                                                                  
      IF (MODN .EQ. 0) X(NS2+1) = W(NS2)*XH(NS2+1)                              
      CALL RFFTF (N,X,XH)                                                       
      DO 103 I=3,N,2                                                            
         XIM1 = X(I-1)-X(I)                                                     
         X(I) = X(I-1)+X(I)                                                     
         X(I-1) = XIM1                                                          
  103 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE COSQB (N,X,WSAVE)                                              
      DIMENSION       X(1)       ,WSAVE(1)                                      
      DATA TSQRT2 /2.82842712474619/                                            
      IF (N-2) 101,102,103                                                      
  101 X(1) = 4.*X(1)                                                            
      RETURN                                                                    
  102 X1 = 4.*(X(1)+X(2))                                                       
      X(2) = TSQRT2*(X(1)-X(2))                                                 
      X(1) = X1                                                                 
      RETURN                                                                    
  103 CALL COSQB1 (N,X,WSAVE,WSAVE(N+1))                                        
      RETURN                                                                    
      END                                                                       
      SUBROUTINE COSQB1 (N,X,W,XH)                                              
      DIMENSION       X(1)       ,W(1)       ,XH(1)                             
      NS2 = (N+1)/2                                                             
      NP2 = N+2                                                                 
      DO 101 I=3,N,2                                                            
         XIM1 = X(I-1)+X(I)                                                     
         X(I) = X(I)-X(I-1)                                                     
         X(I-1) = XIM1                                                          
  101 CONTINUE                                                                  
      X(1) = X(1)+X(1)                                                          
      MODN = MOD(N,2)                                                           
      IF (MODN .EQ. 0) X(N) = X(N)+X(N)                                         
      CALL RFFTB (N,X,XH)                                                       
      DO 102 K=2,NS2                                                            
         KC = NP2-K                                                             
         XH(K) = W(K-1)*X(KC)+W(KC-1)*X(K)                                      
         XH(KC) = W(K-1)*X(K)-W(KC-1)*X(KC)                                     
  102 CONTINUE                                                                  
      IF (MODN .EQ. 0) X(NS2+1) = W(NS2)*(X(NS2+1)+X(NS2+1))                    
      DO 103 K=2,NS2                                                            
         KC = NP2-K                                                             
         X(K) = XH(K)+XH(KC)                                                    
         X(KC) = XH(K)-XH(KC)                                                   
  103 CONTINUE                                                                  
      X(1) = X(1)+X(1)                                                          
      RETURN                                                                    
      END                                                                       
      SUBROUTINE SINQI (N,WSAVE)                                                
      DIMENSION       WSAVE(1)                                                  
      CALL COSQI (N,WSAVE)                                                      
      RETURN                                                                    
      END                                                                       
      SUBROUTINE SINQF (N,X,WSAVE)                                              
      DIMENSION       X(1)       ,WSAVE(1)                                      
      IF (N .EQ. 1) RETURN                                                      
      NS2 = N/2                                                                 
      DO 101 K=1,NS2                                                            
         KC = N-K                                                               
         XHOLD = X(K)                                                           
         X(K) = X(KC+1)                                                         
         X(KC+1) = XHOLD                                                        
  101 CONTINUE                                                                  
      CALL COSQF (N,X,WSAVE)                                                    
      DO 102 K=2,N,2                                                            
         X(K) = -X(K)                                                           
  102 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE SINQB (N,X,WSAVE)                                              
      DIMENSION       X(1)       ,WSAVE(1)                                      
      IF (N .GT. 1) GO TO 101                                                   
      X(1) = 4.*X(1)                                                            
      RETURN                                                                    
  101 NS2 = N/2                                                                 
      DO 102 K=2,N,2                                                            
         X(K) = -X(K)                                                           
  102 CONTINUE                                                                  
      CALL COSQB (N,X,WSAVE)                                                    
      DO 103 K=1,NS2                                                            
         KC = N-K                                                               
         XHOLD = X(K)                                                           
         X(K) = X(KC+1)                                                         
         X(KC+1) = XHOLD                                                        
  103 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE CFFTI (N,WSAVE)                                                
      DIMENSION       WSAVE(1)                                                  
      IF (N .EQ. 1) RETURN                                                      
      IW1 = N+N+1                                                               
      IW2 = IW1+N+N                                                             
      CALL CFFTI1 (N,WSAVE(IW1),WSAVE(IW2))                                     
      RETURN                                                                    
      END                                                                       
      SUBROUTINE CFFTI1 (N,WA,IFAC)                                             
      DIMENSION       WA(1)      ,IFAC(1)    ,NTRYH(4)                          
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/3,4,2,5/                         
      NL = N                                                                    
      NF = 0                                                                    
      J = 0                                                                     
  101 J = J+1                                                                   
      IF (J-4) 102,102,103                                                      
  102 NTRY = NTRYH(J)                                                           
      GO TO 104                                                                 
  103 NTRY = NTRY+2                                                             
  104 NQ = NL/NTRY                                                              
      NR = NL-NTRY*NQ                                                           
      IF (NR) 101,105,101                                                       
  105 NF = NF+1                                                                 
      IFAC(NF+2) = NTRY                                                         
      NL = NQ                                                                   
      IF (NTRY .NE. 2) GO TO 107                                                
      IF (NF .EQ. 1) GO TO 107                                                  
      DO 106 I=2,NF                                                             
         IB = NF-I+2                                                            
         IFAC(IB+2) = IFAC(IB+1)                                                
  106 CONTINUE                                                                  
      IFAC(3) = 2                                                               
  107 IF (NL .NE. 1) GO TO 104                                                  
      IFAC(1) = N                                                               
      IFAC(2) = NF                                                              
      TPI = 6.28318530717959                                                    
      ARGH = TPI/FLOAT(N)                                                       
      I = 2                                                                     
      L1 = 1                                                                    
      DO 110 K1=1,NF                                                            
         IP = IFAC(K1+2)                                                        
         LD = 0                                                                 
         L2 = L1*IP                                                             
         IDO = N/L2                                                             
         IDOT = IDO+IDO+2                                                       
         IPM = IP-1                                                             
         DO 109 J=1,IPM                                                         
            I1 = I                                                              
            WA(I-1) = 1.                                                        
            WA(I) = 0.                                                          
            LD = LD+L1                                                          
            FI = 0.                                                             
            ARGLD = FLOAT(LD)*ARGH                                              
            DO 108 II=4,IDOT,2                                                  
               I = I+2                                                          
               FI = FI+1.                                                       
               ARG = FI*ARGLD                                                   
               WA(I-1) = COS(ARG)                                               
               WA(I) = SIN(ARG)                                                 
  108       CONTINUE                                                            
            IF (IP .LE. 5) GO TO 109                                            
            WA(I1-1) = WA(I-1)                                                  
            WA(I1) = WA(I)                                                      
  109    CONTINUE                                                               
         L1 = L2                                                                
  110 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE CFFTB (N,C,WSAVE)                                              
      DIMENSION       C(1)       ,WSAVE(1)                                      
      IF (N .EQ. 1) RETURN                                                      
      IW1 = N+N+1                                                               
      IW2 = IW1+N+N                                                             
      CALL CFFTB1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))                             
      RETURN                                                                    
      END                                                                       
      SUBROUTINE CFFTB1 (N,C,CH,WA,IFAC)                                        
      DIMENSION       CH(1)      ,C(1)       ,WA(1)      ,IFAC(1)               
      NF = IFAC(2)                                                              
      NA = 0                                                                    
      L1 = 1                                                                    
      IW = 1                                                                    
      DO 116 K1=1,NF                                                            
         IP = IFAC(K1+2)                                                        
         L2 = IP*L1                                                             
         IDO = N/L2                                                             
         IDOT = IDO+IDO                                                         
         IDL1 = IDOT*L1                                                         
         IF (IP .NE. 4) GO TO 103                                               
         IX2 = IW+IDOT                                                          
         IX3 = IX2+IDOT                                                         
         IF (NA .NE. 0) GO TO 101                                               
         CALL PASSB4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))                      
         GO TO 102                                                              
  101    CALL PASSB4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))                      
  102    NA = 1-NA                                                              
         GO TO 115                                                              
  103    IF (IP .NE. 2) GO TO 106                                               
         IF (NA .NE. 0) GO TO 104                                               
         CALL PASSB2 (IDOT,L1,C,CH,WA(IW))                                      
         GO TO 105                                                              
  104    CALL PASSB2 (IDOT,L1,CH,C,WA(IW))                                      
  105    NA = 1-NA                                                              
         GO TO 115                                                              
  106    IF (IP .NE. 3) GO TO 109                                               
         IX2 = IW+IDOT                                                          
         IF (NA .NE. 0) GO TO 107                                               
         CALL PASSB3 (IDOT,L1,C,CH,WA(IW),WA(IX2))                              
         GO TO 108                                                              
  107    CALL PASSB3 (IDOT,L1,CH,C,WA(IW),WA(IX2))                              
  108    NA = 1-NA                                                              
         GO TO 115                                                              
  109    IF (IP .NE. 5) GO TO 112                                               
         IX2 = IW+IDOT                                                          
         IX3 = IX2+IDOT                                                         
         IX4 = IX3+IDOT                                                         
         IF (NA .NE. 0) GO TO 110                                               
         CALL PASSB5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))              
         GO TO 111                                                              
  110    CALL PASSB5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))              
  111    NA = 1-NA                                                              
         GO TO 115                                                              
  112    IF (NA .NE. 0) GO TO 113                                               
         CALL PASSB (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))                    
         GO TO 114                                                              
  113    CALL PASSB (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))                   
  114    IF (NAC .NE. 0) NA = 1-NA                                              
  115    L1 = L2                                                                
         IW = IW+(IP-1)*IDOT                                                    
  116 CONTINUE                                                                  
      IF (NA .EQ. 0) RETURN                                                     
      N2 = N+N                                                                  
      DO 117 I=1,N2                                                             
         C(I) = CH(I)                                                           
  117 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE PASSB2 (IDO,L1,CC,CH,WA1)                                      
      DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,          
     1                WA1(1)                                                    
      IF (IDO .GT. 2) GO TO 102                                                 
      DO 101 K=1,L1                                                             
         CH(1,K,1) = CC(1,1,K)+CC(1,2,K)                                        
         CH(1,K,2) = CC(1,1,K)-CC(1,2,K)                                        
         CH(2,K,1) = CC(2,1,K)+CC(2,2,K)                                        
         CH(2,K,2) = CC(2,1,K)-CC(2,2,K)                                        
  101 CONTINUE                                                                  
      RETURN                                                                    
  102 DO 104 K=1,L1                                                             
         DO 103 I=2,IDO,2                                                       
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)                               
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)                                       
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)                                     
            TI2 = CC(I,1,K)-CC(I,2,K)                                           
            CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2                                 
            CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2                               
  103    CONTINUE                                                               
  104 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE PASSB3 (IDO,L1,CC,CH,WA1,WA2)                                  
      DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,          
     1                WA1(1)     ,WA2(1)                                        
      DATA TAUR,TAUI /-.5,.866025403784439/                                     
      IF (IDO .NE. 2) GO TO 102                                                 
      DO 101 K=1,L1                                                             
         TR2 = CC(1,2,K)+CC(1,3,K)                                              
         CR2 = CC(1,1,K)+TAUR*TR2                                               
         CH(1,K,1) = CC(1,1,K)+TR2                                              
         TI2 = CC(2,2,K)+CC(2,3,K)                                              
         CI2 = CC(2,1,K)+TAUR*TI2                                               
         CH(2,K,1) = CC(2,1,K)+TI2                                              
         CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))                                       
         CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))                                       
         CH(1,K,2) = CR2-CI3                                                    
         CH(1,K,3) = CR2+CI3                                                    
         CH(2,K,2) = CI2+CR3                                                    
         CH(2,K,3) = CI2-CR3                                                    
  101 CONTINUE                                                                  
      RETURN                                                                    
  102 DO 104 K=1,L1                                                             
         DO 103 I=2,IDO,2                                                       
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)                                       
            CR2 = CC(I-1,1,K)+TAUR*TR2                                          
            CH(I-1,K,1) = CC(I-1,1,K)+TR2                                       
            TI2 = CC(I,2,K)+CC(I,3,K)                                           
            CI2 = CC(I,1,K)+TAUR*TI2                                            
            CH(I,K,1) = CC(I,1,K)+TI2                                           
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))                                
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))                                    
            DR2 = CR2-CI3                                                       
            DR3 = CR2+CI3                                                       
            DI2 = CI2+CR3                                                       
            DI3 = CI2-CR3                                                       
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2                                 
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2                               
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3                                 
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3                               
  103    CONTINUE                                                               
  104 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE PASSB4 (IDO,L1,CC,CH,WA1,WA2,WA3)                              
      DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           ,          
     1                WA1(1)     ,WA2(1)     ,WA3(1)                            
      IF (IDO .NE. 2) GO TO 102                                                 
      DO 101 K=1,L1                                                             
         TI1 = CC(2,1,K)-CC(2,3,K)                                              
         TI2 = CC(2,1,K)+CC(2,3,K)                                              
         TR4 = CC(2,4,K)-CC(2,2,K)                                              
         TI3 = CC(2,2,K)+CC(2,4,K)                                              
         TR1 = CC(1,1,K)-CC(1,3,K)                                              
         TR2 = CC(1,1,K)+CC(1,3,K)                                              
         TI4 = CC(1,2,K)-CC(1,4,K)                                              
         TR3 = CC(1,2,K)+CC(1,4,K)                                              
         CH(1,K,1) = TR2+TR3                                                    
         CH(1,K,3) = TR2-TR3                                                    
         CH(2,K,1) = TI2+TI3                                                    
         CH(2,K,3) = TI2-TI3                                                    
         CH(1,K,2) = TR1+TR4                                                    
         CH(1,K,4) = TR1-TR4                                                    
         CH(2,K,2) = TI1+TI4                                                    
         CH(2,K,4) = TI1-TI4                                                    
  101 CONTINUE                                                                  
      RETURN                                                                    
  102 DO 104 K=1,L1                                                             
         DO 103 I=2,IDO,2                                                       
            TI1 = CC(I,1,K)-CC(I,3,K)                                           
            TI2 = CC(I,1,K)+CC(I,3,K)                                           
            TI3 = CC(I,2,K)+CC(I,4,K)                                           
            TR4 = CC(I,4,K)-CC(I,2,K)                                           
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)                                       
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)                                       
            TI4 = CC(I-1,2,K)-CC(I-1,4,K)                                       
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)                                       
            CH(I-1,K,1) = TR2+TR3                                               
            CR3 = TR2-TR3                                                       
            CH(I,K,1) = TI2+TI3                                                 
            CI3 = TI2-TI3                                                       
            CR2 = TR1+TR4                                                       
            CR4 = TR1-TR4                                                       
            CI2 = TI1+TI4                                                       
            CI4 = TI1-TI4                                                       
            CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2                               
            CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2                                 
            CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3                               
            CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3                                 
            CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4                               
            CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4                                 
  103    CONTINUE                                                               
  104 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE PASSB5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)                          
      DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,          
     1                WA1(1)     ,WA2(1)     ,WA3(1)     ,WA4(1)                
      DATA TR11,TI11,TR12,TI12 /.309016994374947,.951056516295154,              
     1-.809016994374947,.587785252292473/                                       
      IF (IDO .NE. 2) GO TO 102                                                 
      DO 101 K=1,L1                                                             
         TI5 = CC(2,2,K)-CC(2,5,K)                                              
         TI2 = CC(2,2,K)+CC(2,5,K)                                              
         TI4 = CC(2,3,K)-CC(2,4,K)                                              
         TI3 = CC(2,3,K)+CC(2,4,K)                                              
         TR5 = CC(1,2,K)-CC(1,5,K)                                              
         TR2 = CC(1,2,K)+CC(1,5,K)                                              
         TR4 = CC(1,3,K)-CC(1,4,K)                                              
         TR3 = CC(1,3,K)+CC(1,4,K)                                              
         CH(1,K,1) = CC(1,1,K)+TR2+TR3                                          
         CH(2,K,1) = CC(2,1,K)+TI2+TI3                                          
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3                                      
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3                                      
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3                                      
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3                                      
         CR5 = TI11*TR5+TI12*TR4                                                
         CI5 = TI11*TI5+TI12*TI4                                                
         CR4 = TI12*TR5-TI11*TR4                                                
         CI4 = TI12*TI5-TI11*TI4                                                
         CH(1,K,2) = CR2-CI5                                                    
         CH(1,K,5) = CR2+CI5                                                    
         CH(2,K,2) = CI2+CR5                                                    
         CH(2,K,3) = CI3+CR4                                                    
         CH(1,K,3) = CR3-CI4                                                    
         CH(1,K,4) = CR3+CI4                                                    
         CH(2,K,4) = CI3-CR4                                                    
         CH(2,K,5) = CI2-CR5                                                    
  101 CONTINUE                                                                  
      RETURN                                                                    
  102 DO 104 K=1,L1                                                             
         DO 103 I=2,IDO,2                                                       
            TI5 = CC(I,2,K)-CC(I,5,K)                                           
            TI2 = CC(I,2,K)+CC(I,5,K)                                           
            TI4 = CC(I,3,K)-CC(I,4,K)                                           
            TI3 = CC(I,3,K)+CC(I,4,K)                                           
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)                                       
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)                                       
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)                                       
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)                                       
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3                                   
            CH(I,K,1) = CC(I,1,K)+TI2+TI3                                       
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3                                 
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3                                   
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3                                 
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3                                   
            CR5 = TI11*TR5+TI12*TR4                                             
            CI5 = TI11*TI5+TI12*TI4                                             
            CR4 = TI12*TR5-TI11*TR4                                             
            CI4 = TI12*TI5-TI11*TI4                                             
            DR3 = CR3-CI4                                                       
            DR4 = CR3+CI4                                                       
            DI3 = CI3+CR4                                                       
            DI4 = CI3-CR4                                                       
            DR5 = CR2+CI5                                                       
            DR2 = CR2-CI5                                                       
            DI5 = CI2-CR5                                                       
            DI2 = CI2+CR5                                                       
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2                               
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2                                 
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3                               
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3                                 
            CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4                               
            CH(I,K,4) = WA3(I-1)*DI4+WA3(I)*DR4                                 
            CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5                               
            CH(I,K,5) = WA4(I-1)*DI5+WA4(I)*DR5                                 
  103    CONTINUE                                                               
  104 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE PASSB (NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)                  
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,          
     1                C1(IDO,L1,IP)          ,WA(1)      ,C2(IDL1,IP),          
     2                CH2(IDL1,IP)                                              
      IDOT = IDO/2                                                              
      NT = IP*IDL1                                                              
      IPP2 = IP+2                                                               
      IPPH = (IP+1)/2                                                           
      IDP = IP*IDO                                                              
C                                                                               
      IF (IDO .LT. L1) GO TO 106                                                
      DO 103 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         DO 102 K=1,L1                                                          
            DO 101 I=1,IDO                                                      
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)                                 
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)                                
  101       CONTINUE                                                            
  102    CONTINUE                                                               
  103 CONTINUE                                                                  
      DO 105 K=1,L1                                                             
         DO 104 I=1,IDO                                                         
            CH(I,K,1) = CC(I,1,K)                                               
  104    CONTINUE                                                               
  105 CONTINUE                                                                  
      GO TO 112                                                                 
  106 DO 109 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         DO 108 I=1,IDO                                                         
            DO 107 K=1,L1                                                       
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)                                 
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)                                
  107       CONTINUE                                                            
  108    CONTINUE                                                               
  109 CONTINUE                                                                  
      DO 111 I=1,IDO                                                            
         DO 110 K=1,L1                                                          
            CH(I,K,1) = CC(I,1,K)                                               
  110    CONTINUE                                                               
  111 CONTINUE                                                                  
  112 IDL = 2-IDO                                                               
      INC = 0                                                                   
      DO 116 L=2,IPPH                                                           
         LC = IPP2-L                                                            
         IDL = IDL+IDO                                                          
         DO 113 IK=1,IDL1                                                       
            C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)                            
            C2(IK,LC) = WA(IDL)*CH2(IK,IP)                                      
  113    CONTINUE                                                               
         IDLJ = IDL                                                             
         INC = INC+IDO                                                          
         DO 115 J=3,IPPH                                                        
            JC = IPP2-J                                                         
            IDLJ = IDLJ+INC                                                     
            IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP                                  
            WAR = WA(IDLJ-1)                                                    
            WAI = WA(IDLJ)                                                      
            DO 114 IK=1,IDL1                                                    
               C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)                                
               C2(IK,LC) = C2(IK,LC)+WAI*CH2(IK,JC)                             
  114       CONTINUE                                                            
  115    CONTINUE                                                               
  116 CONTINUE                                                                  
      DO 118 J=2,IPPH                                                           
         DO 117 IK=1,IDL1                                                       
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)                                     
  117    CONTINUE                                                               
  118 CONTINUE                                                                  
      DO 120 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         DO 119 IK=2,IDL1,2                                                     
            CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)                                  
            CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)                                 
            CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)                                    
            CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)                                   
  119    CONTINUE                                                               
  120 CONTINUE                                                                  
      NAC = 1                                                                   
      IF (IDO .EQ. 2) RETURN                                                    
      NAC = 0                                                                   
      DO 121 IK=1,IDL1                                                          
         C2(IK,1) = CH2(IK,1)                                                   
  121 CONTINUE                                                                  
      DO 123 J=2,IP                                                             
         DO 122 K=1,L1                                                          
            C1(1,K,J) = CH(1,K,J)                                               
            C1(2,K,J) = CH(2,K,J)                                               
  122    CONTINUE                                                               
  123 CONTINUE                                                                  
      IF (IDOT .GT. L1) GO TO 127                                               
      IDIJ = 0                                                                  
      DO 126 J=2,IP                                                             
         IDIJ = IDIJ+2                                                          
         DO 125 I=4,IDO,2                                                       
            IDIJ = IDIJ+2                                                       
            DO 124 K=1,L1                                                       
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)          
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)            
  124       CONTINUE                                                            
  125    CONTINUE                                                               
  126 CONTINUE                                                                  
      RETURN                                                                    
  127 IDJ = 2-IDO                                                               
      DO 130 J=2,IP                                                             
         IDJ = IDJ+IDO                                                          
         DO 129 K=1,L1                                                          
            IDIJ = IDJ                                                          
            DO 128 I=4,IDO,2                                                    
               IDIJ = IDIJ+2                                                    
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)          
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)            
  128       CONTINUE                                                            
  129    CONTINUE                                                               
  130 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE CFFTF (N,C,WSAVE)                                              
      DIMENSION       C(1)       ,WSAVE(1)                                      
      IF (N .EQ. 1) RETURN                                                      
      IW1 = N+N+1                                                               
      IW2 = IW1+N+N                                                             
      CALL CFFTF1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))                             
      RETURN                                                                    
      END                                                                       
      SUBROUTINE CFFTF1 (N,C,CH,WA,IFAC)                                        
      DIMENSION       CH(1)      ,C(1)       ,WA(1)      ,IFAC(1)               
      NF = IFAC(2)                                                              
      NA = 0                                                                    
      L1 = 1                                                                    
      IW = 1                                                                    
      DO 116 K1=1,NF                                                            
         IP = IFAC(K1+2)                                                        
         L2 = IP*L1                                                             
         IDO = N/L2                                                             
         IDOT = IDO+IDO                                                         
         IDL1 = IDOT*L1                                                         
         IF (IP .NE. 4) GO TO 103                                               
         IX2 = IW+IDOT                                                          
         IX3 = IX2+IDOT                                                         
         IF (NA .NE. 0) GO TO 101                                               
         CALL PASSF4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))                      
         GO TO 102                                                              
  101    CALL PASSF4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))                      
  102    NA = 1-NA                                                              
         GO TO 115                                                              
  103    IF (IP .NE. 2) GO TO 106                                               
         IF (NA .NE. 0) GO TO 104                                               
         CALL PASSF2 (IDOT,L1,C,CH,WA(IW))                                      
         GO TO 105                                                              
  104    CALL PASSF2 (IDOT,L1,CH,C,WA(IW))                                      
  105    NA = 1-NA                                                              
         GO TO 115                                                              
  106    IF (IP .NE. 3) GO TO 109                                               
         IX2 = IW+IDOT                                                          
         IF (NA .NE. 0) GO TO 107                                               
         CALL PASSF3 (IDOT,L1,C,CH,WA(IW),WA(IX2))                              
         GO TO 108                                                              
  107    CALL PASSF3 (IDOT,L1,CH,C,WA(IW),WA(IX2))                              
  108    NA = 1-NA                                                              
         GO TO 115                                                              
  109    IF (IP .NE. 5) GO TO 112                                               
         IX2 = IW+IDOT                                                          
         IX3 = IX2+IDOT                                                         
         IX4 = IX3+IDOT                                                         
         IF (NA .NE. 0) GO TO 110                                               
         CALL PASSF5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))              
         GO TO 111                                                              
  110    CALL PASSF5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))              
  111    NA = 1-NA                                                              
         GO TO 115                                                              
  112    IF (NA .NE. 0) GO TO 113                                               
         CALL PASSF (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))                    
         GO TO 114                                                              
  113    CALL PASSF (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))                   
  114    IF (NAC .NE. 0) NA = 1-NA                                              
  115    L1 = L2                                                                
         IW = IW+(IP-1)*IDOT                                                    
  116 CONTINUE                                                                  
      IF (NA .EQ. 0) RETURN                                                     
      N2 = N+N                                                                  
      DO 117 I=1,N2                                                             
         C(I) = CH(I)                                                           
  117 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE PASSF2 (IDO,L1,CC,CH,WA1)                                      
      DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,          
     1                WA1(1)                                                    
      IF (IDO .GT. 2) GO TO 102                                                 
      DO 101 K=1,L1                                                             
         CH(1,K,1) = CC(1,1,K)+CC(1,2,K)                                        
         CH(1,K,2) = CC(1,1,K)-CC(1,2,K)                                        
         CH(2,K,1) = CC(2,1,K)+CC(2,2,K)                                        
         CH(2,K,2) = CC(2,1,K)-CC(2,2,K)                                        
  101 CONTINUE                                                                  
      RETURN                                                                    
  102 DO 104 K=1,L1                                                             
         DO 103 I=2,IDO,2                                                       
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)                               
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)                                       
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)                                     
            TI2 = CC(I,1,K)-CC(I,2,K)                                           
            CH(I,K,2) = WA1(I-1)*TI2-WA1(I)*TR2                                 
            CH(I-1,K,2) = WA1(I-1)*TR2+WA1(I)*TI2                               
  103    CONTINUE                                                               
  104 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE PASSF3 (IDO,L1,CC,CH,WA1,WA2)                                  
      DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,          
     1                WA1(1)     ,WA2(1)                                        
      DATA TAUR,TAUI /-.5,-.866025403784439/                                    
      IF (IDO .NE. 2) GO TO 102                                                 
      DO 101 K=1,L1                                                             
         TR2 = CC(1,2,K)+CC(1,3,K)                                              
         CR2 = CC(1,1,K)+TAUR*TR2                                               
         CH(1,K,1) = CC(1,1,K)+TR2                                              
         TI2 = CC(2,2,K)+CC(2,3,K)                                              
         CI2 = CC(2,1,K)+TAUR*TI2                                               
         CH(2,K,1) = CC(2,1,K)+TI2                                              
         CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))                                       
         CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))                                       
         CH(1,K,2) = CR2-CI3                                                    
         CH(1,K,3) = CR2+CI3                                                    
         CH(2,K,2) = CI2+CR3                                                    
         CH(2,K,3) = CI2-CR3                                                    
  101 CONTINUE                                                                  
      RETURN                                                                    
  102 DO 104 K=1,L1                                                             
         DO 103 I=2,IDO,2                                                       
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)                                       
            CR2 = CC(I-1,1,K)+TAUR*TR2                                          
            CH(I-1,K,1) = CC(I-1,1,K)+TR2                                       
            TI2 = CC(I,2,K)+CC(I,3,K)                                           
            CI2 = CC(I,1,K)+TAUR*TI2                                            
            CH(I,K,1) = CC(I,1,K)+TI2                                           
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))                                
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))                                    
            DR2 = CR2-CI3                                                       
            DR3 = CR2+CI3                                                       
            DI2 = CI2+CR3                                                       
            DI3 = CI2-CR3                                                       
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2                                 
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2                               
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3                                 
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3                               
  103    CONTINUE                                                               
  104 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE PASSF4 (IDO,L1,CC,CH,WA1,WA2,WA3)                              
      DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           ,          
     1                WA1(1)     ,WA2(1)     ,WA3(1)                            
      IF (IDO .NE. 2) GO TO 102                                                 
      DO 101 K=1,L1                                                             
         TI1 = CC(2,1,K)-CC(2,3,K)                                              
         TI2 = CC(2,1,K)+CC(2,3,K)                                              
         TR4 = CC(2,2,K)-CC(2,4,K)                                              
         TI3 = CC(2,2,K)+CC(2,4,K)                                              
         TR1 = CC(1,1,K)-CC(1,3,K)                                              
         TR2 = CC(1,1,K)+CC(1,3,K)                                              
         TI4 = CC(1,4,K)-CC(1,2,K)                                              
         TR3 = CC(1,2,K)+CC(1,4,K)                                              
         CH(1,K,1) = TR2+TR3                                                    
         CH(1,K,3) = TR2-TR3                                                    
         CH(2,K,1) = TI2+TI3                                                    
         CH(2,K,3) = TI2-TI3                                                    
         CH(1,K,2) = TR1+TR4                                                    
         CH(1,K,4) = TR1-TR4                                                    
         CH(2,K,2) = TI1+TI4                                                    
         CH(2,K,4) = TI1-TI4                                                    
  101 CONTINUE                                                                  
      RETURN                                                                    
  102 DO 104 K=1,L1                                                             
         DO 103 I=2,IDO,2                                                       
            TI1 = CC(I,1,K)-CC(I,3,K)                                           
            TI2 = CC(I,1,K)+CC(I,3,K)                                           
            TI3 = CC(I,2,K)+CC(I,4,K)                                           
            TR4 = CC(I,2,K)-CC(I,4,K)                                           
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)                                       
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)                                       
            TI4 = CC(I-1,4,K)-CC(I-1,2,K)                                       
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)                                       
            CH(I-1,K,1) = TR2+TR3                                               
            CR3 = TR2-TR3                                                       
            CH(I,K,1) = TI2+TI3                                                 
            CI3 = TI2-TI3                                                       
            CR2 = TR1+TR4                                                       
            CR4 = TR1-TR4                                                       
            CI2 = TI1+TI4                                                       
            CI4 = TI1-TI4                                                       
            CH(I-1,K,2) = WA1(I-1)*CR2+WA1(I)*CI2                               
            CH(I,K,2) = WA1(I-1)*CI2-WA1(I)*CR2                                 
            CH(I-1,K,3) = WA2(I-1)*CR3+WA2(I)*CI3                               
            CH(I,K,3) = WA2(I-1)*CI3-WA2(I)*CR3                                 
            CH(I-1,K,4) = WA3(I-1)*CR4+WA3(I)*CI4                               
            CH(I,K,4) = WA3(I-1)*CI4-WA3(I)*CR4                                 
  103    CONTINUE                                                               
  104 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE PASSF5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)                          
      DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,          
     1                WA1(1)     ,WA2(1)     ,WA3(1)     ,WA4(1)                
      DATA TR11,TI11,TR12,TI12 /.309016994374947,-.951056516295154,             
     1-.809016994374947,-.587785252292473/                                      
      IF (IDO .NE. 2) GO TO 102                                                 
      DO 101 K=1,L1                                                             
         TI5 = CC(2,2,K)-CC(2,5,K)                                              
         TI2 = CC(2,2,K)+CC(2,5,K)                                              
         TI4 = CC(2,3,K)-CC(2,4,K)                                              
         TI3 = CC(2,3,K)+CC(2,4,K)                                              
         TR5 = CC(1,2,K)-CC(1,5,K)                                              
         TR2 = CC(1,2,K)+CC(1,5,K)                                              
         TR4 = CC(1,3,K)-CC(1,4,K)                                              
         TR3 = CC(1,3,K)+CC(1,4,K)                                              
         CH(1,K,1) = CC(1,1,K)+TR2+TR3                                          
         CH(2,K,1) = CC(2,1,K)+TI2+TI3                                          
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3                                      
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3                                      
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3                                      
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3                                      
         CR5 = TI11*TR5+TI12*TR4                                                
         CI5 = TI11*TI5+TI12*TI4                                                
         CR4 = TI12*TR5-TI11*TR4                                                
         CI4 = TI12*TI5-TI11*TI4                                                
         CH(1,K,2) = CR2-CI5                                                    
         CH(1,K,5) = CR2+CI5                                                    
         CH(2,K,2) = CI2+CR5                                                    
         CH(2,K,3) = CI3+CR4                                                    
         CH(1,K,3) = CR3-CI4                                                    
         CH(1,K,4) = CR3+CI4                                                    
         CH(2,K,4) = CI3-CR4                                                    
         CH(2,K,5) = CI2-CR5                                                    
  101 CONTINUE                                                                  
      RETURN                                                                    
  102 DO 104 K=1,L1                                                             
         DO 103 I=2,IDO,2                                                       
            TI5 = CC(I,2,K)-CC(I,5,K)                                           
            TI2 = CC(I,2,K)+CC(I,5,K)                                           
            TI4 = CC(I,3,K)-CC(I,4,K)                                           
            TI3 = CC(I,3,K)+CC(I,4,K)                                           
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)                                       
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)                                       
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)                                       
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)                                       
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3                                   
            CH(I,K,1) = CC(I,1,K)+TI2+TI3                                       
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3                                 
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3                                   
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3                                 
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3                                   
            CR5 = TI11*TR5+TI12*TR4                                             
            CI5 = TI11*TI5+TI12*TI4                                             
            CR4 = TI12*TR5-TI11*TR4                                             
            CI4 = TI12*TI5-TI11*TI4                                             
            DR3 = CR3-CI4                                                       
            DR4 = CR3+CI4                                                       
            DI3 = CI3+CR4                                                       
            DI4 = CI3-CR4                                                       
            DR5 = CR2+CI5                                                       
            DR2 = CR2-CI5                                                       
            DI5 = CI2-CR5                                                       
            DI2 = CI2+CR5                                                       
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2                               
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2                                 
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3                               
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3                                 
            CH(I-1,K,4) = WA3(I-1)*DR4+WA3(I)*DI4                               
            CH(I,K,4) = WA3(I-1)*DI4-WA3(I)*DR4                                 
            CH(I-1,K,5) = WA4(I-1)*DR5+WA4(I)*DI5                               
            CH(I,K,5) = WA4(I-1)*DI5-WA4(I)*DR5                                 
  103    CONTINUE                                                               
  104 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE PASSF (NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)                  
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,          
     1                C1(IDO,L1,IP)          ,WA(1)      ,C2(IDL1,IP),          
     2                CH2(IDL1,IP)                                              
      IDOT = IDO/2                                                              
      NT = IP*IDL1                                                              
      IPP2 = IP+2                                                               
      IPPH = (IP+1)/2                                                           
      IDP = IP*IDO                                                              
C                                                                               
      IF (IDO .LT. L1) GO TO 106                                                
      DO 103 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         DO 102 K=1,L1                                                          
            DO 101 I=1,IDO                                                      
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)                                 
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)                                
  101       CONTINUE                                                            
  102    CONTINUE                                                               
  103 CONTINUE                                                                  
      DO 105 K=1,L1                                                             
         DO 104 I=1,IDO                                                         
            CH(I,K,1) = CC(I,1,K)                                               
  104    CONTINUE                                                               
  105 CONTINUE                                                                  
      GO TO 112                                                                 
  106 DO 109 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         DO 108 I=1,IDO                                                         
            DO 107 K=1,L1                                                       
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)                                 
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)                                
  107       CONTINUE                                                            
  108    CONTINUE                                                               
  109 CONTINUE                                                                  
      DO 111 I=1,IDO                                                            
         DO 110 K=1,L1                                                          
            CH(I,K,1) = CC(I,1,K)                                               
  110    CONTINUE                                                               
  111 CONTINUE                                                                  
  112 IDL = 2-IDO                                                               
      INC = 0                                                                   
      DO 116 L=2,IPPH                                                           
         LC = IPP2-L                                                            
         IDL = IDL+IDO                                                          
         DO 113 IK=1,IDL1                                                       
            C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)                            
            C2(IK,LC) = -WA(IDL)*CH2(IK,IP)                                     
  113    CONTINUE                                                               
         IDLJ = IDL                                                             
         INC = INC+IDO                                                          
         DO 115 J=3,IPPH                                                        
            JC = IPP2-J                                                         
            IDLJ = IDLJ+INC                                                     
            IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP                                  
            WAR = WA(IDLJ-1)                                                    
            WAI = WA(IDLJ)                                                      
            DO 114 IK=1,IDL1                                                    
               C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)                                
               C2(IK,LC) = C2(IK,LC)-WAI*CH2(IK,JC)                             
  114       CONTINUE                                                            
  115    CONTINUE                                                               
  116 CONTINUE                                                                  
      DO 118 J=2,IPPH                                                           
         DO 117 IK=1,IDL1                                                       
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)                                     
  117    CONTINUE                                                               
  118 CONTINUE                                                                  
      DO 120 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         DO 119 IK=2,IDL1,2                                                     
            CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)                                  
            CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)                                 
            CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)                                    
            CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)                                   
  119    CONTINUE                                                               
  120 CONTINUE                                                                  
      NAC = 1                                                                   
      IF (IDO .EQ. 2) RETURN                                                    
      NAC = 0                                                                   
      DO 121 IK=1,IDL1                                                          
         C2(IK,1) = CH2(IK,1)                                                   
  121 CONTINUE                                                                  
      DO 123 J=2,IP                                                             
         DO 122 K=1,L1                                                          
            C1(1,K,J) = CH(1,K,J)                                               
            C1(2,K,J) = CH(2,K,J)                                               
  122    CONTINUE                                                               
  123 CONTINUE                                                                  
      IF (IDOT .GT. L1) GO TO 127                                               
      IDIJ = 0                                                                  
      DO 126 J=2,IP                                                             
         IDIJ = IDIJ+2                                                          
         DO 125 I=4,IDO,2                                                       
            IDIJ = IDIJ+2                                                       
            DO 124 K=1,L1                                                       
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)          
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)            
  124       CONTINUE                                                            
  125    CONTINUE                                                               
  126 CONTINUE                                                                  
      RETURN                                                                    
  127 IDJ = 2-IDO                                                               
      DO 130 J=2,IP                                                             
         IDJ = IDJ+IDO                                                          
         DO 129 K=1,L1                                                          
            IDIJ = IDJ                                                          
            DO 128 I=4,IDO,2                                                    
               IDIJ = IDIJ+2                                                    
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)          
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)            
  128       CONTINUE                                                            
  129    CONTINUE                                                               
  130 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE RFFTI (N,WSAVE)                                                
      DIMENSION       WSAVE(1)                                                  
      IF (N .EQ. 1) RETURN                                                      
      CALL RFFTI1 (N,WSAVE(N+1),WSAVE(2*N+1))                                   
      RETURN                                                                    
      END                                                                       
      SUBROUTINE RFFTI1 (N,WA,IFAC)                                             
      DIMENSION       WA(1)      ,IFAC(1)    ,NTRYH(4)                          
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/                         
      NL = N                                                                    
      NF = 0                                                                    
      J = 0                                                                     
  101 J = J+1                                                                   
      IF (J-4) 102,102,103                                                      
  102 NTRY = NTRYH(J)                                                           
      GO TO 104                                                                 
  103 NTRY = NTRY+2                                                             
  104 NQ = NL/NTRY                                                              
      NR = NL-NTRY*NQ                                                           
      IF (NR) 101,105,101                                                       
  105 NF = NF+1                                                                 
      IFAC(NF+2) = NTRY                                                         
      NL = NQ                                                                   
      IF (NTRY .NE. 2) GO TO 107                                                
      IF (NF .EQ. 1) GO TO 107                                                  
      DO 106 I=2,NF                                                             
         IB = NF-I+2                                                            
         IFAC(IB+2) = IFAC(IB+1)                                                
  106 CONTINUE                                                                  
      IFAC(3) = 2                                                               
  107 IF (NL .NE. 1) GO TO 104                                                  
      IFAC(1) = N                                                               
      IFAC(2) = NF                                                              
      TPI = 6.28318530717959                                                    
      ARGH = TPI/FLOAT(N)                                                       
      IS = 0                                                                    
      NFM1 = NF-1                                                               
      L1 = 1                                                                    
      IF (NFM1 .EQ. 0) RETURN                                                   
      DO 110 K1=1,NFM1                                                          
         IP = IFAC(K1+2)                                                        
         LD = 0                                                                 
         L2 = L1*IP                                                             
         IDO = N/L2                                                             
         IPM = IP-1                                                             
         DO 109 J=1,IPM                                                         
            LD = LD+L1                                                          
            I = IS                                                              
            ARGLD = FLOAT(LD)*ARGH                                              
            FI = 0.                                                             
            DO 108 II=3,IDO,2                                                   
               I = I+2                                                          
               FI = FI+1.                                                       
               ARG = FI*ARGLD                                                   
               WA(I-1) = COS(ARG)                                               
               WA(I) = SIN(ARG)                                                 
  108       CONTINUE                                                            
            IS = IS+IDO                                                         
  109    CONTINUE                                                               
         L1 = L2                                                                
  110 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE RFFTB (N,R,WSAVE)                                              
      DIMENSION       R(1)       ,WSAVE(1)                                      
      IF (N .EQ. 1) RETURN                                                      
      CALL RFFTB1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))                           
      RETURN                                                                    
      END                                                                       
      SUBROUTINE RFFTB1 (N,C,CH,WA,IFAC)                                        
      DIMENSION       CH(1)      ,C(1)       ,WA(1)      ,IFAC(1)               
      NF = IFAC(2)                                                              
      NA = 0                                                                    
      L1 = 1                                                                    
      IW = 1                                                                    
      DO 116 K1=1,NF                                                            
         IP = IFAC(K1+2)                                                        
         L2 = IP*L1                                                             
         IDO = N/L2                                                             
         IDL1 = IDO*L1                                                          
         IF (IP .NE. 4) GO TO 103                                               
         IX2 = IW+IDO                                                           
         IX3 = IX2+IDO                                                          
         IF (NA .NE. 0) GO TO 101                                               
         CALL RADB4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))                        
         GO TO 102                                                              
  101    CALL RADB4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))                        
  102    NA = 1-NA                                                              
         GO TO 115                                                              
  103    IF (IP .NE. 2) GO TO 106                                               
         IF (NA .NE. 0) GO TO 104                                               
         CALL RADB2 (IDO,L1,C,CH,WA(IW))                                        
         GO TO 105                                                              
  104    CALL RADB2 (IDO,L1,CH,C,WA(IW))                                        
  105    NA = 1-NA                                                              
         GO TO 115                                                              
  106    IF (IP .NE. 3) GO TO 109                                               
         IX2 = IW+IDO                                                           
         IF (NA .NE. 0) GO TO 107                                               
         CALL RADB3 (IDO,L1,C,CH,WA(IW),WA(IX2))                                
         GO TO 108                                                              
  107    CALL RADB3 (IDO,L1,CH,C,WA(IW),WA(IX2))                                
  108    NA = 1-NA                                                              
         GO TO 115                                                              
  109    IF (IP .NE. 5) GO TO 112                                               
         IX2 = IW+IDO                                                           
         IX3 = IX2+IDO                                                          
         IX4 = IX3+IDO                                                          
         IF (NA .NE. 0) GO TO 110                                               
         CALL RADB5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))                
         GO TO 111                                                              
  110    CALL RADB5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))                
  111    NA = 1-NA                                                              
         GO TO 115                                                              
  112    IF (NA .NE. 0) GO TO 113                                               
         CALL RADBG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))                         
         GO TO 114                                                              
  113    CALL RADBG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))                        
  114    IF (IDO .EQ. 1) NA = 1-NA                                              
  115    L1 = L2                                                                
         IW = IW+(IP-1)*IDO                                                     
  116 CONTINUE                                                                  
      IF (NA .EQ. 0) RETURN                                                     
      DO 117 I=1,N                                                              
         C(I) = CH(I)                                                           
  117 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE RADB2 (IDO,L1,CC,CH,WA1)                                       
      DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,          
     1                WA1(1)                                                    
      DO 101 K=1,L1                                                             
         CH(1,K,1) = CC(1,1,K)+CC(IDO,2,K)                                      
         CH(1,K,2) = CC(1,1,K)-CC(IDO,2,K)                                      
  101 CONTINUE                                                                  
      IF (IDO-2) 107,105,102                                                    
  102 IDP2 = IDO+2                                                              
      DO 104 K=1,L1                                                             
         DO 103 I=3,IDO,2                                                       
            IC = IDP2-I                                                         
            CH(I-1,K,1) = CC(I-1,1,K)+CC(IC-1,2,K)                              
            TR2 = CC(I-1,1,K)-CC(IC-1,2,K)                                      
            CH(I,K,1) = CC(I,1,K)-CC(IC,2,K)                                    
            TI2 = CC(I,1,K)+CC(IC,2,K)                                          
            CH(I-1,K,2) = WA1(I-2)*TR2-WA1(I-1)*TI2                             
            CH(I,K,2) = WA1(I-2)*TI2+WA1(I-1)*TR2                               
  103    CONTINUE                                                               
  104 CONTINUE                                                                  
      IF (MOD(IDO,2) .EQ. 1) RETURN                                             
  105 DO 106 K=1,L1                                                             
         CH(IDO,K,1) = CC(IDO,1,K)+CC(IDO,1,K)                                  
         CH(IDO,K,2) = -(CC(1,2,K)+CC(1,2,K))                                   
  106 CONTINUE                                                                  
  107 RETURN                                                                    
      END                                                                       
      SUBROUTINE RADB3 (IDO,L1,CC,CH,WA1,WA2)                                   
      DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,          
     1                WA1(1)     ,WA2(1)                                        
      DATA TAUR,TAUI /-.5,.866025403784439/                                     
      DO 101 K=1,L1                                                             
         TR2 = CC(IDO,2,K)+CC(IDO,2,K)                                          
         CR2 = CC(1,1,K)+TAUR*TR2                                               
         CH(1,K,1) = CC(1,1,K)+TR2                                              
         CI3 = TAUI*(CC(1,3,K)+CC(1,3,K))                                       
         CH(1,K,2) = CR2-CI3                                                    
         CH(1,K,3) = CR2+CI3                                                    
  101 CONTINUE                                                                  
      IF (IDO .EQ. 1) RETURN                                                    
      IDP2 = IDO+2                                                              
      DO 103 K=1,L1                                                             
         DO 102 I=3,IDO,2                                                       
            IC = IDP2-I                                                         
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)                                      
            CR2 = CC(I-1,1,K)+TAUR*TR2                                          
            CH(I-1,K,1) = CC(I-1,1,K)+TR2                                       
            TI2 = CC(I,3,K)-CC(IC,2,K)                                          
            CI2 = CC(I,1,K)+TAUR*TI2                                            
            CH(I,K,1) = CC(I,1,K)+TI2                                           
            CR3 = TAUI*(CC(I-1,3,K)-CC(IC-1,2,K))                               
            CI3 = TAUI*(CC(I,3,K)+CC(IC,2,K))                                   
            DR2 = CR2-CI3                                                       
            DR3 = CR2+CI3                                                       
            DI2 = CI2+CR3                                                       
            DI3 = CI2-CR3                                                       
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2                             
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2                               
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3                             
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3                               
  102    CONTINUE                                                               
  103 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE RADB4 (IDO,L1,CC,CH,WA1,WA2,WA3)                               
      DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           ,          
     1                WA1(1)     ,WA2(1)     ,WA3(1)                            
      DATA SQRT2 /1.414213562373095/                                            
      DO 101 K=1,L1                                                             
         TR1 = CC(1,1,K)-CC(IDO,4,K)                                            
         TR2 = CC(1,1,K)+CC(IDO,4,K)                                            
         TR3 = CC(IDO,2,K)+CC(IDO,2,K)                                          
         TR4 = CC(1,3,K)+CC(1,3,K)                                              
         CH(1,K,1) = TR2+TR3                                                    
         CH(1,K,2) = TR1-TR4                                                    
         CH(1,K,3) = TR2-TR3                                                    
         CH(1,K,4) = TR1+TR4                                                    
  101 CONTINUE                                                                  
      IF (IDO-2) 107,105,102                                                    
  102 IDP2 = IDO+2                                                              
      DO 104 K=1,L1                                                             
         DO 103 I=3,IDO,2                                                       
            IC = IDP2-I                                                         
            TI1 = CC(I,1,K)+CC(IC,4,K)                                          
            TI2 = CC(I,1,K)-CC(IC,4,K)                                          
            TI3 = CC(I,3,K)-CC(IC,2,K)                                          
            TR4 = CC(I,3,K)+CC(IC,2,K)                                          
            TR1 = CC(I-1,1,K)-CC(IC-1,4,K)                                      
            TR2 = CC(I-1,1,K)+CC(IC-1,4,K)                                      
            TI4 = CC(I-1,3,K)-CC(IC-1,2,K)                                      
            TR3 = CC(I-1,3,K)+CC(IC-1,2,K)                                      
            CH(I-1,K,1) = TR2+TR3                                               
            CR3 = TR2-TR3                                                       
            CH(I,K,1) = TI2+TI3                                                 
            CI3 = TI2-TI3                                                       
            CR2 = TR1-TR4                                                       
            CR4 = TR1+TR4                                                       
            CI2 = TI1+TI4                                                       
            CI4 = TI1-TI4                                                       
            CH(I-1,K,2) = WA1(I-2)*CR2-WA1(I-1)*CI2                             
            CH(I,K,2) = WA1(I-2)*CI2+WA1(I-1)*CR2                               
            CH(I-1,K,3) = WA2(I-2)*CR3-WA2(I-1)*CI3                             
            CH(I,K,3) = WA2(I-2)*CI3+WA2(I-1)*CR3                               
            CH(I-1,K,4) = WA3(I-2)*CR4-WA3(I-1)*CI4                             
            CH(I,K,4) = WA3(I-2)*CI4+WA3(I-1)*CR4                               
  103    CONTINUE                                                               
  104 CONTINUE                                                                  
      IF (MOD(IDO,2) .EQ. 1) RETURN                                             
  105 CONTINUE                                                                  
      DO 106 K=1,L1                                                             
         TI1 = CC(1,2,K)+CC(1,4,K)                                              
         TI2 = CC(1,4,K)-CC(1,2,K)                                              
         TR1 = CC(IDO,1,K)-CC(IDO,3,K)                                          
         TR2 = CC(IDO,1,K)+CC(IDO,3,K)                                          
         CH(IDO,K,1) = TR2+TR2                                                  
         CH(IDO,K,2) = SQRT2*(TR1-TI1)                                          
         CH(IDO,K,3) = TI2+TI2                                                  
         CH(IDO,K,4) = -SQRT2*(TR1+TI1)                                         
  106 CONTINUE                                                                  
  107 RETURN                                                                    
      END                                                                       
      SUBROUTINE RADB5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)                           
      DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,          
     1                WA1(1)     ,WA2(1)     ,WA3(1)     ,WA4(1)                
      DATA TR11,TI11,TR12,TI12 /.309016994374947,.951056516295154,              
     1-.809016994374947,.587785252292473/                                       
      DO 101 K=1,L1                                                             
         TI5 = CC(1,3,K)+CC(1,3,K)                                              
         TI4 = CC(1,5,K)+CC(1,5,K)                                              
         TR2 = CC(IDO,2,K)+CC(IDO,2,K)                                          
         TR3 = CC(IDO,4,K)+CC(IDO,4,K)                                          
         CH(1,K,1) = CC(1,1,K)+TR2+TR3                                          
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3                                      
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3                                      
         CI5 = TI11*TI5+TI12*TI4                                                
         CI4 = TI12*TI5-TI11*TI4                                                
         CH(1,K,2) = CR2-CI5                                                    
         CH(1,K,3) = CR3-CI4                                                    
         CH(1,K,4) = CR3+CI4                                                    
         CH(1,K,5) = CR2+CI5                                                    
  101 CONTINUE                                                                  
      IF (IDO .EQ. 1) RETURN                                                    
      IDP2 = IDO+2                                                              
      DO 103 K=1,L1                                                             
         DO 102 I=3,IDO,2                                                       
            IC = IDP2-I                                                         
            TI5 = CC(I,3,K)+CC(IC,2,K)                                          
            TI2 = CC(I,3,K)-CC(IC,2,K)                                          
            TI4 = CC(I,5,K)+CC(IC,4,K)                                          
            TI3 = CC(I,5,K)-CC(IC,4,K)                                          
            TR5 = CC(I-1,3,K)-CC(IC-1,2,K)                                      
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)                                      
            TR4 = CC(I-1,5,K)-CC(IC-1,4,K)                                      
            TR3 = CC(I-1,5,K)+CC(IC-1,4,K)                                      
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3                                   
            CH(I,K,1) = CC(I,1,K)+TI2+TI3                                       
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3                                 
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3                                   
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3                                 
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3                                   
            CR5 = TI11*TR5+TI12*TR4                                             
            CI5 = TI11*TI5+TI12*TI4                                             
            CR4 = TI12*TR5-TI11*TR4                                             
            CI4 = TI12*TI5-TI11*TI4                                             
            DR3 = CR3-CI4                                                       
            DR4 = CR3+CI4                                                       
            DI3 = CI3+CR4                                                       
            DI4 = CI3-CR4                                                       
            DR5 = CR2+CI5                                                       
            DR2 = CR2-CI5                                                       
            DI5 = CI2-CR5                                                       
            DI2 = CI2+CR5                                                       
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2                             
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2                               
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3                             
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3                               
            CH(I-1,K,4) = WA3(I-2)*DR4-WA3(I-1)*DI4                             
            CH(I,K,4) = WA3(I-2)*DI4+WA3(I-1)*DR4                               
            CH(I-1,K,5) = WA4(I-2)*DR5-WA4(I-1)*DI5                             
            CH(I,K,5) = WA4(I-2)*DI5+WA4(I-1)*DR5                               
  102    CONTINUE                                                               
  103 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE RADBG (IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)                      
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,          
     1                C1(IDO,L1,IP)          ,C2(IDL1,IP),                      
     2                CH2(IDL1,IP)           ,WA(1)                             
      DATA TPI/6.28318530717959/                                                
      ARG = TPI/FLOAT(IP)                                                       
      DCP = COS(ARG)                                                            
      DSP = SIN(ARG)                                                            
      IDP2 = IDO+2                                                              
      NBD = (IDO-1)/2                                                           
      IPP2 = IP+2                                                               
      IPPH = (IP+1)/2                                                           
      IF (IDO .LT. L1) GO TO 103                                                
      DO 102 K=1,L1                                                             
         DO 101 I=1,IDO                                                         
            CH(I,K,1) = CC(I,1,K)                                               
  101    CONTINUE                                                               
  102 CONTINUE                                                                  
      GO TO 106                                                                 
  103 DO 105 I=1,IDO                                                            
         DO 104 K=1,L1                                                          
            CH(I,K,1) = CC(I,1,K)                                               
  104    CONTINUE                                                               
  105 CONTINUE                                                                  
  106 DO 108 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         J2 = J+J                                                               
         DO 107 K=1,L1                                                          
            CH(1,K,J) = CC(IDO,J2-2,K)+CC(IDO,J2-2,K)                           
            CH(1,K,JC) = CC(1,J2-1,K)+CC(1,J2-1,K)                              
  107    CONTINUE                                                               
  108 CONTINUE                                                                  
      IF (IDO .EQ. 1) GO TO 116                                                 
      IF (NBD .LT. L1) GO TO 112                                                
      DO 111 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         DO 110 K=1,L1                                                          
            DO 109 I=3,IDO,2                                                    
               IC = IDP2-I                                                      
               CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)                   
               CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)                  
               CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)                         
               CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)                        
  109       CONTINUE                                                            
  110    CONTINUE                                                               
  111 CONTINUE                                                                  
      GO TO 116                                                                 
  112 DO 115 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         DO 114 I=3,IDO,2                                                       
            IC = IDP2-I                                                         
            DO 113 K=1,L1                                                       
               CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)                   
               CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)                  
               CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)                         
               CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)                        
  113       CONTINUE                                                            
  114    CONTINUE                                                               
  115 CONTINUE                                                                  
  116 AR1 = 1.                                                                  
      AI1 = 0.                                                                  
      DO 120 L=2,IPPH                                                           
         LC = IPP2-L                                                            
         AR1H = DCP*AR1-DSP*AI1                                                 
         AI1 = DCP*AI1+DSP*AR1                                                  
         AR1 = AR1H                                                             
         DO 117 IK=1,IDL1                                                       
            C2(IK,L) = CH2(IK,1)+AR1*CH2(IK,2)                                  
            C2(IK,LC) = AI1*CH2(IK,IP)                                          
  117    CONTINUE                                                               
         DC2 = AR1                                                              
         DS2 = AI1                                                              
         AR2 = AR1                                                              
         AI2 = AI1                                                              
         DO 119 J=3,IPPH                                                        
            JC = IPP2-J                                                         
            AR2H = DC2*AR2-DS2*AI2                                              
            AI2 = DC2*AI2+DS2*AR2                                               
            AR2 = AR2H                                                          
            DO 118 IK=1,IDL1                                                    
               C2(IK,L) = C2(IK,L)+AR2*CH2(IK,J)                                
               C2(IK,LC) = C2(IK,LC)+AI2*CH2(IK,JC)                             
  118       CONTINUE                                                            
  119    CONTINUE                                                               
  120 CONTINUE                                                                  
      DO 122 J=2,IPPH                                                           
         DO 121 IK=1,IDL1                                                       
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)                                     
  121    CONTINUE                                                               
  122 CONTINUE                                                                  
      DO 124 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         DO 123 K=1,L1                                                          
            CH(1,K,J) = C1(1,K,J)-C1(1,K,JC)                                    
            CH(1,K,JC) = C1(1,K,J)+C1(1,K,JC)                                   
  123    CONTINUE                                                               
  124 CONTINUE                                                                  
      IF (IDO .EQ. 1) GO TO 132                                                 
      IF (NBD .LT. L1) GO TO 128                                                
      DO 127 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         DO 126 K=1,L1                                                          
            DO 125 I=3,IDO,2                                                    
               CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)                             
               CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)                            
               CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)                               
               CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)                              
  125       CONTINUE                                                            
  126    CONTINUE                                                               
  127 CONTINUE                                                                  
      GO TO 132                                                                 
  128 DO 131 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         DO 130 I=3,IDO,2                                                       
            DO 129 K=1,L1                                                       
               CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)                             
               CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)                            
               CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)                               
               CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)                              
  129       CONTINUE                                                            
  130    CONTINUE                                                               
  131 CONTINUE                                                                  
  132 CONTINUE                                                                  
      IF (IDO .EQ. 1) RETURN                                                    
      DO 133 IK=1,IDL1                                                          
         C2(IK,1) = CH2(IK,1)                                                   
  133 CONTINUE                                                                  
      DO 135 J=2,IP                                                             
         DO 134 K=1,L1                                                          
            C1(1,K,J) = CH(1,K,J)                                               
  134    CONTINUE                                                               
  135 CONTINUE                                                                  
      IF (NBD .GT. L1) GO TO 139                                                
      IS = -IDO                                                                 
      DO 138 J=2,IP                                                             
         IS = IS+IDO                                                            
         IDIJ = IS                                                              
         DO 137 I=3,IDO,2                                                       
            IDIJ = IDIJ+2                                                       
            DO 136 K=1,L1                                                       
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)          
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)            
  136       CONTINUE                                                            
  137    CONTINUE                                                               
  138 CONTINUE                                                                  
      GO TO 143                                                                 
  139 IS = -IDO                                                                 
      DO 142 J=2,IP                                                             
         IS = IS+IDO                                                            
         DO 141 K=1,L1                                                          
            IDIJ = IS                                                           
            DO 140 I=3,IDO,2                                                    
               IDIJ = IDIJ+2                                                    
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)          
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)            
  140       CONTINUE                                                            
  141    CONTINUE                                                               
  142 CONTINUE                                                                  
  143 RETURN                                                                    
      END                                                                       
      SUBROUTINE RFFTF (N,R,WSAVE)                                              
      DIMENSION       R(1)       ,WSAVE(1)                                      
      IF (N .EQ. 1) RETURN                                                      
      CALL RFFTF1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))                           
      RETURN                                                                    
      END                                                                       
      SUBROUTINE RFFTF1 (N,C,CH,WA,IFAC)                                        
      DIMENSION       CH(1)      ,C(1)       ,WA(1)      ,IFAC(1)               
      NF = IFAC(2)                                                              
      NA = 1                                                                    
      L2 = N                                                                    
      IW = N                                                                    
      DO 111 K1=1,NF                                                            
         KH = NF-K1                                                             
         IP = IFAC(KH+3)                                                        
         L1 = L2/IP                                                             
         IDO = N/L2                                                             
         IDL1 = IDO*L1                                                          
         IW = IW-(IP-1)*IDO                                                     
         NA = 1-NA                                                              
         IF (IP .NE. 4) GO TO 102                                               
         IX2 = IW+IDO                                                           
         IX3 = IX2+IDO                                                          
         IF (NA .NE. 0) GO TO 101                                               
         CALL RADF4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))                        
         GO TO 110                                                              
  101    CALL RADF4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))                        
         GO TO 110                                                              
  102    IF (IP .NE. 2) GO TO 104                                               
         IF (NA .NE. 0) GO TO 103                                               
         CALL RADF2 (IDO,L1,C,CH,WA(IW))                                        
         GO TO 110                                                              
  103    CALL RADF2 (IDO,L1,CH,C,WA(IW))                                        
         GO TO 110                                                              
  104    IF (IP .NE. 3) GO TO 106                                               
         IX2 = IW+IDO                                                           
         IF (NA .NE. 0) GO TO 105                                               
         CALL RADF3 (IDO,L1,C,CH,WA(IW),WA(IX2))                                
         GO TO 110                                                              
  105    CALL RADF3 (IDO,L1,CH,C,WA(IW),WA(IX2))                                
         GO TO 110                                                              
  106    IF (IP .NE. 5) GO TO 108                                               
         IX2 = IW+IDO                                                           
         IX3 = IX2+IDO                                                          
         IX4 = IX3+IDO                                                          
         IF (NA .NE. 0) GO TO 107                                               
         CALL RADF5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))                
         GO TO 110                                                              
  107    CALL RADF5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))                
         GO TO 110                                                              
  108    IF (IDO .EQ. 1) NA = 1-NA                                              
         IF (NA .NE. 0) GO TO 109                                               
         CALL RADFG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))                         
         NA = 1                                                                 
         GO TO 110                                                              
  109    CALL RADFG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))                        
         NA = 0                                                                 
  110    L2 = L1                                                                
  111 CONTINUE                                                                  
      IF (NA .EQ. 1) RETURN                                                     
      DO 112 I=1,N                                                              
         C(I) = CH(I)                                                           
  112 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE RADF2 (IDO,L1,CC,CH,WA1)                                       
      DIMENSION       CH(IDO,2,L1)           ,CC(IDO,L1,2)           ,          
     1                WA1(1)                                                    
      DO 101 K=1,L1                                                             
         CH(1,1,K) = CC(1,K,1)+CC(1,K,2)                                        
         CH(IDO,2,K) = CC(1,K,1)-CC(1,K,2)                                      
  101 CONTINUE                                                                  
      IF (IDO-2) 107,105,102                                                    
  102 IDP2 = IDO+2                                                              
      DO 104 K=1,L1                                                             
         DO 103 I=3,IDO,2                                                       
            IC = IDP2-I                                                         
            TR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)                       
            TI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)                       
            CH(I,1,K) = CC(I,K,1)+TI2                                           
            CH(IC,2,K) = TI2-CC(I,K,1)                                          
            CH(I-1,1,K) = CC(I-1,K,1)+TR2                                       
            CH(IC-1,2,K) = CC(I-1,K,1)-TR2                                      
  103    CONTINUE                                                               
  104 CONTINUE                                                                  
      IF (MOD(IDO,2) .EQ. 1) RETURN                                             
  105 DO 106 K=1,L1                                                             
         CH(1,2,K) = -CC(IDO,K,2)                                               
         CH(IDO,1,K) = CC(IDO,K,1)                                              
  106 CONTINUE                                                                  
  107 RETURN                                                                    
      END                                                                       
      SUBROUTINE RADF3 (IDO,L1,CC,CH,WA1,WA2)                                   
      DIMENSION       CH(IDO,3,L1)           ,CC(IDO,L1,3)           ,          
     1                WA1(1)     ,WA2(1)                                        
      DATA TAUR,TAUI /-.5,.866025403784439/                                     
      DO 101 K=1,L1                                                             
         CR2 = CC(1,K,2)+CC(1,K,3)                                              
         CH(1,1,K) = CC(1,K,1)+CR2                                              
         CH(1,3,K) = TAUI*(CC(1,K,3)-CC(1,K,2))                                 
         CH(IDO,2,K) = CC(1,K,1)+TAUR*CR2                                       
  101 CONTINUE                                                                  
      IF (IDO .EQ. 1) RETURN                                                    
      IDP2 = IDO+2                                                              
      DO 103 K=1,L1                                                             
         DO 102 I=3,IDO,2                                                       
            IC = IDP2-I                                                         
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)                       
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)                       
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)                       
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)                       
            CR2 = DR2+DR3                                                       
            CI2 = DI2+DI3                                                       
            CH(I-1,1,K) = CC(I-1,K,1)+CR2                                       
            CH(I,1,K) = CC(I,K,1)+CI2                                           
            TR2 = CC(I-1,K,1)+TAUR*CR2                                          
            TI2 = CC(I,K,1)+TAUR*CI2                                            
            TR3 = TAUI*(DI2-DI3)                                                
            TI3 = TAUI*(DR3-DR2)                                                
            CH(I-1,3,K) = TR2+TR3                                               
            CH(IC-1,2,K) = TR2-TR3                                              
            CH(I,3,K) = TI2+TI3                                                 
            CH(IC,2,K) = TI3-TI2                                                
  102    CONTINUE                                                               
  103 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE RADF4 (IDO,L1,CC,CH,WA1,WA2,WA3)                               
      DIMENSION       CC(IDO,L1,4)           ,CH(IDO,4,L1)           ,          
     1                WA1(1)     ,WA2(1)     ,WA3(1)                            
      DATA HSQT2 /.7071067811865475/                                            
      DO 101 K=1,L1                                                             
         TR1 = CC(1,K,2)+CC(1,K,4)                                              
         TR2 = CC(1,K,1)+CC(1,K,3)                                              
         CH(1,1,K) = TR1+TR2                                                    
         CH(IDO,4,K) = TR2-TR1                                                  
         CH(IDO,2,K) = CC(1,K,1)-CC(1,K,3)                                      
         CH(1,3,K) = CC(1,K,4)-CC(1,K,2)                                        
  101 CONTINUE                                                                  
      IF (IDO-2) 107,105,102                                                    
  102 IDP2 = IDO+2                                                              
      DO 104 K=1,L1                                                             
         DO 103 I=3,IDO,2                                                       
            IC = IDP2-I                                                         
            CR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)                       
            CI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)                       
            CR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)                       
            CI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)                       
            CR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)                       
            CI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)                       
            TR1 = CR2+CR4                                                       
            TR4 = CR4-CR2                                                       
            TI1 = CI2+CI4                                                       
            TI4 = CI2-CI4                                                       
            TI2 = CC(I,K,1)+CI3                                                 
            TI3 = CC(I,K,1)-CI3                                                 
            TR2 = CC(I-1,K,1)+CR3                                               
            TR3 = CC(I-1,K,1)-CR3                                               
            CH(I-1,1,K) = TR1+TR2                                               
            CH(IC-1,4,K) = TR2-TR1                                              
            CH(I,1,K) = TI1+TI2                                                 
            CH(IC,4,K) = TI1-TI2                                                
            CH(I-1,3,K) = TI4+TR3                                               
            CH(IC-1,2,K) = TR3-TI4                                              
            CH(I,3,K) = TR4+TI3                                                 
            CH(IC,2,K) = TR4-TI3                                                
  103    CONTINUE                                                               
  104 CONTINUE                                                                  
      IF (MOD(IDO,2) .EQ. 1) RETURN                                             
  105 CONTINUE                                                                  
      DO 106 K=1,L1                                                             
         TI1 = -HSQT2*(CC(IDO,K,2)+CC(IDO,K,4))                                 
         TR1 = HSQT2*(CC(IDO,K,2)-CC(IDO,K,4))                                  
         CH(IDO,1,K) = TR1+CC(IDO,K,1)                                          
         CH(IDO,3,K) = CC(IDO,K,1)-TR1                                          
         CH(1,2,K) = TI1-CC(IDO,K,3)                                            
         CH(1,4,K) = TI1+CC(IDO,K,3)                                            
  106 CONTINUE                                                                  
  107 RETURN                                                                    
      END                                                                       
      SUBROUTINE RADF5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)                           
      DIMENSION       CC(IDO,L1,5)           ,CH(IDO,5,L1)           ,          
     1                WA1(1)     ,WA2(1)     ,WA3(1)     ,WA4(1)                
      DATA TR11,TI11,TR12,TI12 /.309016994374947,.951056516295154,              
     1-.809016994374947,.587785252292473/                                       
      DO 101 K=1,L1                                                             
         CR2 = CC(1,K,5)+CC(1,K,2)                                              
         CI5 = CC(1,K,5)-CC(1,K,2)                                              
         CR3 = CC(1,K,4)+CC(1,K,3)                                              
         CI4 = CC(1,K,4)-CC(1,K,3)                                              
         CH(1,1,K) = CC(1,K,1)+CR2+CR3                                          
         CH(IDO,2,K) = CC(1,K,1)+TR11*CR2+TR12*CR3                              
         CH(1,3,K) = TI11*CI5+TI12*CI4                                          
         CH(IDO,4,K) = CC(1,K,1)+TR12*CR2+TR11*CR3                              
         CH(1,5,K) = TI12*CI5-TI11*CI4                                          
  101 CONTINUE                                                                  
      IF (IDO .EQ. 1) RETURN                                                    
      IDP2 = IDO+2                                                              
      DO 103 K=1,L1                                                             
         DO 102 I=3,IDO,2                                                       
            IC = IDP2-I                                                         
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)                       
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)                       
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)                       
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)                       
            DR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)                       
            DI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)                       
            DR5 = WA4(I-2)*CC(I-1,K,5)+WA4(I-1)*CC(I,K,5)                       
            DI5 = WA4(I-2)*CC(I,K,5)-WA4(I-1)*CC(I-1,K,5)                       
            CR2 = DR2+DR5                                                       
            CI5 = DR5-DR2                                                       
            CR5 = DI2-DI5                                                       
            CI2 = DI2+DI5                                                       
            CR3 = DR3+DR4                                                       
            CI4 = DR4-DR3                                                       
            CR4 = DI3-DI4                                                       
            CI3 = DI3+DI4                                                       
            CH(I-1,1,K) = CC(I-1,K,1)+CR2+CR3                                   
            CH(I,1,K) = CC(I,K,1)+CI2+CI3                                       
            TR2 = CC(I-1,K,1)+TR11*CR2+TR12*CR3                                 
            TI2 = CC(I,K,1)+TR11*CI2+TR12*CI3                                   
            TR3 = CC(I-1,K,1)+TR12*CR2+TR11*CR3                                 
            TI3 = CC(I,K,1)+TR12*CI2+TR11*CI3                                   
            TR5 = TI11*CR5+TI12*CR4                                             
            TI5 = TI11*CI5+TI12*CI4                                             
            TR4 = TI12*CR5-TI11*CR4                                             
            TI4 = TI12*CI5-TI11*CI4                                             
            CH(I-1,3,K) = TR2+TR5                                               
            CH(IC-1,2,K) = TR2-TR5                                              
            CH(I,3,K) = TI2+TI5                                                 
            CH(IC,2,K) = TI5-TI2                                                
            CH(I-1,5,K) = TR3+TR4                                               
            CH(IC-1,4,K) = TR3-TR4                                              
            CH(I,5,K) = TI3+TI4                                                 
            CH(IC,4,K) = TI4-TI3                                                
  102    CONTINUE                                                               
  103 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE RADFG (IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)                      
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,          
     1                C1(IDO,L1,IP)          ,C2(IDL1,IP),                      
     2                CH2(IDL1,IP)           ,WA(1)                             
      DATA TPI/6.28318530717959/                                                
      ARG = TPI/FLOAT(IP)                                                       
      DCP = COS(ARG)                                                            
      DSP = SIN(ARG)                                                            
      IPPH = (IP+1)/2                                                           
      IPP2 = IP+2                                                               
      IDP2 = IDO+2                                                              
      NBD = (IDO-1)/2                                                           
      IF (IDO .EQ. 1) GO TO 119                                                 
      DO 101 IK=1,IDL1                                                          
         CH2(IK,1) = C2(IK,1)                                                   
  101 CONTINUE                                                                  
      DO 103 J=2,IP                                                             
         DO 102 K=1,L1                                                          
            CH(1,K,J) = C1(1,K,J)                                               
  102    CONTINUE                                                               
  103 CONTINUE                                                                  
      IF (NBD .GT. L1) GO TO 107                                                
      IS = -IDO                                                                 
      DO 106 J=2,IP                                                             
         IS = IS+IDO                                                            
         IDIJ = IS                                                              
         DO 105 I=3,IDO,2                                                       
            IDIJ = IDIJ+2                                                       
            DO 104 K=1,L1                                                       
               CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)          
               CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)            
  104       CONTINUE                                                            
  105    CONTINUE                                                               
  106 CONTINUE                                                                  
      GO TO 111                                                                 
  107 IS = -IDO                                                                 
      DO 110 J=2,IP                                                             
         IS = IS+IDO                                                            
         DO 109 K=1,L1                                                          
            IDIJ = IS                                                           
            DO 108 I=3,IDO,2                                                    
               IDIJ = IDIJ+2                                                    
               CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)          
               CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)            
  108       CONTINUE                                                            
  109    CONTINUE                                                               
  110 CONTINUE                                                                  
  111 IF (NBD .LT. L1) GO TO 115                                                
      DO 114 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         DO 113 K=1,L1                                                          
            DO 112 I=3,IDO,2                                                    
               C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)                           
               C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)                              
               C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)                                 
               C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)                            
  112       CONTINUE                                                            
  113    CONTINUE                                                               
  114 CONTINUE                                                                  
      GO TO 121                                                                 
  115 DO 118 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         DO 117 I=3,IDO,2                                                       
            DO 116 K=1,L1                                                       
               C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)                           
               C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)                              
               C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)                                 
               C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)                            
  116       CONTINUE                                                            
  117    CONTINUE                                                               
  118 CONTINUE                                                                  
      GO TO 121                                                                 
  119 DO 120 IK=1,IDL1                                                          
         C2(IK,1) = CH2(IK,1)                                                   
  120 CONTINUE                                                                  
  121 DO 123 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         DO 122 K=1,L1                                                          
            C1(1,K,J) = CH(1,K,J)+CH(1,K,JC)                                    
            C1(1,K,JC) = CH(1,K,JC)-CH(1,K,J)                                   
  122    CONTINUE                                                               
  123 CONTINUE                                                                  
C                                                                               
      AR1 = 1.                                                                  
      AI1 = 0.                                                                  
      DO 127 L=2,IPPH                                                           
         LC = IPP2-L                                                            
         AR1H = DCP*AR1-DSP*AI1                                                 
         AI1 = DCP*AI1+DSP*AR1                                                  
         AR1 = AR1H                                                             
         DO 124 IK=1,IDL1                                                       
            CH2(IK,L) = C2(IK,1)+AR1*C2(IK,2)                                   
            CH2(IK,LC) = AI1*C2(IK,IP)                                          
  124    CONTINUE                                                               
         DC2 = AR1                                                              
         DS2 = AI1                                                              
         AR2 = AR1                                                              
         AI2 = AI1                                                              
         DO 126 J=3,IPPH                                                        
            JC = IPP2-J                                                         
            AR2H = DC2*AR2-DS2*AI2                                              
            AI2 = DC2*AI2+DS2*AR2                                               
            AR2 = AR2H                                                          
            DO 125 IK=1,IDL1                                                    
               CH2(IK,L) = CH2(IK,L)+AR2*C2(IK,J)                               
               CH2(IK,LC) = CH2(IK,LC)+AI2*C2(IK,JC)                            
  125       CONTINUE                                                            
  126    CONTINUE                                                               
  127 CONTINUE                                                                  
      DO 129 J=2,IPPH                                                           
         DO 128 IK=1,IDL1                                                       
            CH2(IK,1) = CH2(IK,1)+C2(IK,J)                                      
  128    CONTINUE                                                               
  129 CONTINUE                                                                  
C                                                                               
      IF (IDO .LT. L1) GO TO 132                                                
      DO 131 K=1,L1                                                             
         DO 130 I=1,IDO                                                         
            CC(I,1,K) = CH(I,K,1)                                               
  130    CONTINUE                                                               
  131 CONTINUE                                                                  
      GO TO 135                                                                 
  132 DO 134 I=1,IDO                                                            
         DO 133 K=1,L1                                                          
            CC(I,1,K) = CH(I,K,1)                                               
  133    CONTINUE                                                               
  134 CONTINUE                                                                  
  135 DO 137 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         J2 = J+J                                                               
         DO 136 K=1,L1                                                          
            CC(IDO,J2-2,K) = CH(1,K,J)                                          
            CC(1,J2-1,K) = CH(1,K,JC)                                           
  136    CONTINUE                                                               
  137 CONTINUE                                                                  
      IF (IDO .EQ. 1) RETURN                                                    
      IF (NBD .LT. L1) GO TO 141                                                
      DO 140 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         J2 = J+J                                                               
         DO 139 K=1,L1                                                          
            DO 138 I=3,IDO,2                                                    
               IC = IDP2-I                                                      
               CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)                        
               CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)                       
               CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)                              
               CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)                             
  138       CONTINUE                                                            
  139    CONTINUE                                                               
  140 CONTINUE                                                                  
      RETURN                                                                    
  141 DO 144 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         J2 = J+J                                                               
         DO 143 I=3,IDO,2                                                       
            IC = IDP2-I                                                         
            DO 142 K=1,L1                                                       
               CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)                        
               CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)                       
               CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)                              
               CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)                             
  142       CONTINUE                                                            
  143    CONTINUE                                                               
  144 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       

C FISHPAK24  FROM PORTLIB                                  03/11/81             
      FUNCTION EPMACH (DUM)                                                     
C                                                                               
C     THIS PROGRAM COMPUTES AN APPROXIMATE MACHIINE EPSILON (ACCURACY)          
C                                                                               
      COMMON /VALUE/  V                                                         
      EPS = 1.                                                                  
  101 EPS = EPS/10.                                                             
      CALL STORE (EPS+1.)                                                       
      IF (V-1.) 102,102,101                                                     
  102 EPMACH = 100.*EPS                                                         
      RETURN                                                                    
      END                                                                       
      SUBROUTINE STORE (X)                                                      
      COMMON /VALUE/  V                                                         
      V = X                                                                     
      RETURN                                                                    
      END                                                                       
      FUNCTION PIMACH (DUM)                                                     
C                                                                               
C     THIS SUBPROGRAM SUPPLIES THE VALUE OF THE CONSTANT PI CORRECT TO          
C     MACHINE PRECISION WHERE                                                   
C                                                                               
C     PI=3.1415926535897932384626433832795028841971693993751058209749446        
C                                                                               
      PIMACH = 3.14159265358979                                                 
      RETURN                                                                    
      END                          
