!***BEGIN PROLOGUE  PCHIC                                               
!***DATE WRITTEN   820218   (YYMMDD)                                    
!***REVISION DATE  870707   (YYMMDD)                                    
!***CATEGORY NO.  E1B                                                   
!***KEYWORDS  LIBRARY=SLATEC(PCHIP),                                    
!             TYPE=SINGLE PRECISION(PCHIC-S DPCHIC-D),                  
!             CUBIC HERMITE INTERPOLATION,MONOTONE INTERPOLATION,       
!             PIECEWISE CUBIC INTERPOLATION,                            
!             SHAPE-PRESERVING INTERPOLATION                            
!***AUTHOR  FRITSCH, F. N., (LLNL)                                      
!             MATHEMATICS AND STATISTICS DIVISION                       
!             LAWRENCE LIVERMORE NATIONAL LABORATORY                    
!             P.O. BOX 808  (L-316)                                     
!             LIVERMORE, CA  94550                                      
!             FTS 532-4275, (415) 422-4275                              
!***PURPOSE  Set derivatives needed to determine a piecewise monotone   
!            piecewise cubic Hermite interpolant to given data.         
!            User control is available over boundary conditions and/or  
!            treatment of points where monotonicity switches direction. 
!***DESCRIPTION                                                         
!                                                                       
!         PCHIC:  Piecewise Cubic Hermite Interpolation Coefficients.   
!                                                                       
!     Sets derivatives needed to determine a piecewise monotone piece-  
!     wise cubic interpolant to the data given in X and F satisfying the
!     boundary conditions specified by IC and VC.                       
!                                                                       
!     The treatment of points where monotonicity switches direction is  
!     controlled by argument SWITCH.                                    
!                                                                       
!     To facilitate two-dimensional applications, includes an increment 
!     between successive values of the F- and D-arrays.                 
!                                                                       
!     The resulting piecewise cubic Hermite function may be evaluated   
!     by PCHFE or PCHFD.                                                
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  Calling sequence:                                                    
!                                                                       
!        PARAMETER  (INCFD = ...)                                       
!        INTEGER  IC(2), N, NWK, IERR                                   
!        REAL  VC(2), SWITCH, X(N), F(INCFD,N), D(INCFD,N), WK(NWK)     
!                                                                       
!        CALL  PCHIC (IC, VC, SWITCH, N, X, F, D, INCFD, WK, NWK, IERR) 
!                                                                       
!   Parameters:                                                         
!                                                                       
!     IC -- (input) integer array of length 2 specifying desired        
!           boundary conditions:                                        
!           IC(1) = IBEG, desired condition at beginning of data.       
!           IC(2) = IEND, desired condition at end of data.             
!                                                                       
!           IBEG = 0  for the default boundary condition (the same as   
!                     used by PCHIM).                                   
!           If IBEG.NE.0, then its sign indicates whether the boundary  
!                     derivative is to be adjusted, if necessary, to be 
!                     compatible with monotonicity:                     
!              IBEG.GT.0  if no adjustment is to be performed.          
!              IBEG.LT.0  if the derivative is to be adjusted for       
!                     monotonicity.                                     
!                                                                       
!           Allowable values for the magnitude of IBEG are:             
!           IBEG = 1  if first derivative at X(1) is given in VC(1).    
!           IBEG = 2  if second derivative at X(1) is given in VC(1).   
!           IBEG = 3  to use the 3-point difference formula for D(1).   
!                     (Reverts to the default b.c. if N.LT.3 .)         
!           IBEG = 4  to use the 4-point difference formula for D(1).   
!                     (Reverts to the default b.c. if N.LT.4 .)         
!           IBEG = 5  to set D(1) so that the second derivative is con- 
!              tinuous at X(2). (Reverts to the default b.c. if N.LT.4.)
!              This option is somewhat analogous to the "not a knot"    
!              boundary condition provided by PCHSP.                    
!                                                                       
!          NOTES (IBEG):                                                
!           1. An error return is taken if ABS(IBEG).GT.5 .             
!           2. Only in case  IBEG.LE.0  is it guaranteed that the       
!              interpolant will be monotonic in the first interval.     
!              If the returned value of D(1) lies between zero and      
!              3*SLOPE(1), the interpolant will be monotonic.  This     
!              is **NOT** checked if IBEG.GT.0 .                        
!           3. If IBEG.LT.0 and D(1) had to be changed to achieve mono- 
!              tonicity, a warning error is returned.                   
!                                                                       
!           IEND may take on the same values as IBEG, but applied to    
!           derivative at X(N).  In case IEND = 1 or 2, the value is    
!           given in VC(2).                                             
!                                                                       
!          NOTES (IEND):                                                
!           1. An error return is taken if ABS(IEND).GT.5 .             
!           2. Only in case  IEND.LE.0  is it guaranteed that the       
!              interpolant will be monotonic in the last interval.      
!              If the returned value of D(1+(N-1)*INCFD) lies between   
!              zero and 3*SLOPE(N-1), the interpolant will be monotonic.
!              This is **NOT** checked if IEND.GT.0 .                   
!           3. If IEND.LT.0 and D(1+(N-1)*INCFD) had to be changed to   
!              achieve monotonicity, a warning error is returned.       
!                                                                       
!     VC -- (input) real array of length 2 specifying desired boundary  
!           values, as indicated above.                                 
!           VC(1) need be set only if IC(1) = 1 or 2 .                  
!           VC(2) need be set only if IC(2) = 1 or 2 .                  
!                                                                       
!     SWITCH -- (input) indicates desired treatment of points where     
!           direction of monotonicity switches:                         
!           Set SWITCH to zero if interpolant is required to be mono-   
!           tonic in each interval, regardless of monotonicity of data. 
!             NOTES:                                                    
!              1. This will cause D to be set to zero at all switch     
!                 points, thus forcing extrema there.                   
!              2. The result of using this option with the default boun-
!                 dary conditions will be identical to using PCHIM, but 
!                 will generally cost more compute time.                
!                 This option is provided only to facilitate comparison 
!                 of different switch and/or boundary conditions.       
!           Set SWITCH nonzero to use a formula based on the 3-point    
!              difference formula in the vicinity of switch points.     
!           If SWITCH is positive, the interpolant on each interval     
!              containing an extremum is controlled to not deviate from 
!              the data by more than SWITCH*DFLOC, where DFLOC is the   
!              maximum of the change of F on this interval and its two  
!              immediate neighbors.                                     
!           If SWITCH is negative, no such control is to be imposed.    
!                                                                       
!     N -- (input) number of data points.  (Error return if N.LT.2 .)   
!                                                                       
!     X -- (input) real array of independent variable values.  The      
!           elements of X must be strictly increasing:                  
!                X(I-1) .LT. X(I),  I = 2(1)N.                          
!           (Error return if not.)                                      
!                                                                       
!     F -- (input) real array of dependent variable values to be inter- 
!           polated.  F(1+(I-1)*INCFD) is value corresponding to X(I).  
!                                                                       
!     D -- (output) real array of derivative values at the data points. 
!           These values will determine a monotone cubic Hermite func-  
!           tion on each subinterval on which the data are monotonic,   
!           except possibly adjacent to switches in monotonicity.       
!           The value corresponding to X(I) is stored in                
!                D(1+(I-1)*INCFD),  I=1(1)N.                            
!           No other entries in D are changed.                          
!                                                                       
!     INCFD -- (input) increment between successive values in F and D.  
!           This argument is provided primarily for 2-D applications.   
!           (Error return if  INCFD.LT.1 .)                             
!                                                                       
!     WK -- (scratch) real array of working storage.  The user may wish 
!           to know that the returned values are:                       
!              WK(I)     = H(I)     = X(I+1) - X(I) ;                   
!              WK(N-1+I) = SLOPE(I) = (F(1,I+1) - F(1,I)) / H(I)        
!           for  I = 1(1)N-1.                                           
!                                                                       
!     NWK -- (input) length of work array.                              
!           (Error return if  NWK.LT.2*(N-1) .)                         
!                                                                       
!     IERR -- (output) error flag.                                      
!           Normal return:                                              
!              IERR = 0  (no errors).                                   
!           Warning errors:                                             
!              IERR = 1  if IBEG.LT.0 and D(1) had to be adjusted for   
!                        monotonicity.                                  
!              IERR = 2  if IEND.LT.0 and D(1+(N-1)*INCFD) had to be    
!                        adjusted for monotonicity.                     
!              IERR = 3  if both of the above are true.                 
!           "Recoverable" errors:                                       
!              IERR = -1  if N.LT.2 .                                   
!              IERR = -2  if INCFD.LT.1 .                               
!              IERR = -3  if the X-array is not strictly increasing.    
!              IERR = -4  if ABS(IBEG).GT.5 .                           
!              IERR = -5  if ABS(IEND).GT.5 .                           
!              IERR = -6  if both of the above are true.                
!              IERR = -7  if NWK.LT.2*(N-1) .                           
!             (The D-array has not been changed in any of these cases.) 
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!                                                                       
!***REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE     
!                 CUBIC INTERPOLATION,' SIAM J.NUMER.ANAL. 17, 2 (APRIL 
!                 1980), 238-246.                                       
!               2. F.N.FRITSCH AND J.BUTLAND, 'A METHOD FOR CONSTRUCTING
!                 LOCAL MONOTONE PIECEWISE CUBIC INTERPOLANTS,' SIAM    
!                 J.SCI.STAT.COMPUT.5,2 (JUNE 1984), 300-304.           
!               3. F.N.FRITSCH, 'PIECEWISE CUBIC INTERPOLATION PACKAGE,'
!                 LLNL PREPRINT UCRL-87285 (JULY 1982).                 
!***ROUTINES CALLED  PCHCE,PCHCI,PCHCS,XERROR                           
!***END PROLOGUE  PCHIC                                                 
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  Change record:                                                       
!     82-08-04   Converted to SLATEC library version.                   
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  Programming notes:                                                   
!                                                                       
!     To produce a double precision version, simply:                    
!        a. Change PCHIC to DPCHIC wherever it occurs,                  
!        b. Change PCHCE to DPCHCE wherever it occurs,                  
!        c. Change PCHCI to DPCHCI wherever it occurs,                  
!        d. Change PCHCS to DPCHCS wherever it occurs,                  
      SUBROUTINE PCHIC(IC,VC,SWITCH,N,X,F,D,INCFD,WK,NWK,IERR) 
!        e. Change the real declarations to double precision, and       
!        f. Change the constant  ZERO  to double precision.             
!                                                                       
!  DECLARE ARGUMENTS.                                                   
!                                                                       
      INTEGER  IC(2), N, INCFD, NWK, IERR 
      REAL  VC(2), SWITCH, X(N), F(INCFD,N), D(INCFD,N), WK(NWK) 
!                                                                       
!  DECLARE LOCAL VARIABLES.                                             
!                                                                       
      INTEGER  I, IBEG, IEND, NLESS1 
      REAL  ZERO 
      DATA  ZERO /0./ 
!                                                                       
!  VALIDITY-CHECK ARGUMENTS.                                            
!                                                                       
!***FIRST EXECUTABLE STATEMENT  PCHIC                                   
      IF ( N.LT.2 )  GO TO 5001 
      IF ( INCFD.LT.1 )  GO TO 5002 
      DO 1  I = 2, N 
         IF ( X(I).LE.X(I-1) )  GO TO 5003 
    1 END DO 
!                                                                       
      IBEG = IC(1) 
      IEND = IC(2) 
      IERR = 0 
      IF (IABS(IBEG) .GT. 5)  IERR = IERR - 1 
      IF (IABS(IEND) .GT. 5)  IERR = IERR - 2 
      IF (IERR .LT. 0)  GO TO 5004 
!                                                                       
!  FUNCTION DEFINITION IS OK -- GO ON.                                  
!                                                                       
      NLESS1 = N - 1 
      IF ( NWK .LT. 2*NLESS1 )  GO TO 5007 
!                                                                       
!  SET UP H AND SLOPE ARRAYS.                                           
!                                                                       
      DO 20  I = 1, NLESS1 
         WK(I) = X(I+1) - X(I) 
         WK(NLESS1+I) = (F(1,I+1) - F(1,I)) / WK(I) 
   20 END DO 
!                                                                       
!  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.                        
!                                                                       
      IF (NLESS1 .GT. 1)  GO TO 1000 
      D(1,1) = WK(2) 
      D(1,N) = WK(2) 
      GO TO 3000 
!                                                                       
!  NORMAL CASE  (N .GE. 3) .                                            
!                                                                       
 1000 CONTINUE 
!                                                                       
!  SET INTERIOR DERIVATIVES AND DEFAULT END CONDITIONS.                 
!                                                                       
!     --------------------------------------                            
      CALL PCHCI (N, WK(1), WK(N), D, INCFD) 
!     --------------------------------------                            
!                                                                       
!  SET DERIVATIVES AT POINTS WHERE MONOTONICITY SWITCHES DIRECTION.     
!                                                                       
      IF (SWITCH .EQ. ZERO)  GO TO 3000 
!     ----------------------------------------------------              
      CALL PCHCS (SWITCH, N, WK(1), WK(N), D, INCFD, IERR) 
!     ----------------------------------------------------              
      IF (IERR .NE. 0)  GO TO 5008 
!                                                                       
!  SET END CONDITIONS.                                                  
!                                                                       
 3000 CONTINUE 
      IF ( (IBEG.EQ.0) .AND. (IEND.EQ.0) )  GO TO 5000 
!     -------------------------------------------------------           
      CALL PCHCE (IC, VC, N, X, WK(1), WK(N), D, INCFD, IERR) 
!     -------------------------------------------------------           
      IF (IERR .LT. 0)  GO TO 5009 
!                                                                       
!  NORMAL RETURN.                                                       
!                                                                       
 5000 CONTINUE 
      RETURN 
!                                                                       
!  ERROR RETURNS.                                                       
!                                                                       
 5001 CONTINUE 
!     N.LT.2 RETURN.                                                    
      IERR = -1 
      CALL XERROR ('PCHIC -- NUMBER OF DATA POINTS LESS THAN TWO'       &
     &           , 44, IERR, 1)                                         
      RETURN 
!                                                                       
 5002 CONTINUE 
!     INCFD.LT.1 RETURN.                                                
      IERR = -2 
      CALL XERROR ('PCHIC -- INCREMENT LESS THAN ONE'                   &
     &           , 32, IERR, 1)                                         
      RETURN 
!                                                                       
 5003 CONTINUE 
!     X-ARRAY NOT STRICTLY INCREASING.                                  
      IERR = -3 
      CALL XERROR ('PCHIC -- X-ARRAY NOT STRICTLY INCREASING'           &
     &           , 40, IERR, 1)                                         
      RETURN 
!                                                                       
 5004 CONTINUE 
!     IC OUT OF RANGE RETURN.                                           
      IERR = IERR - 3 
      CALL XERROR ('PCHIC -- IC OUT OF RANGE'                           &
     &           , 24, IERR, 1)                                         
      RETURN 
!                                                                       
 5007 CONTINUE 
!     NWK .LT. 2*(N-1)  RETURN.                                         
      IERR = -7 
      CALL XERROR ('PCHIC -- WORK ARRAY TOO SMALL'                      &
     &           , 29, IERR, 1)                                         
      RETURN 
!                                                                       
 5008 CONTINUE 
!     ERROR RETURN FROM PCHCS.                                          
      IERR = -8 
      CALL XERROR ('PCHIC -- ERROR RETURN FROM PCHCS'                   &
     &           , 32, IERR, 1)                                         
      RETURN 
!                                                                       
 5009 CONTINUE 
!     ERROR RETURN FROM PCHCE.                                          
!   *** THIS CASE SHOULD NEVER OCCUR ***                                
      IERR = -9 
      CALL XERROR ('PCHIC -- ERROR RETURN FROM PCHCE'                   &
     &           , 32, IERR, 1)                                         
      RETURN 
!------------- LAST LINE OF PCHIC FOLLOWS ------------------------------
      END                                           
      SUBROUTINE PCHCE(IC,VC,N,X,H,SLOPE,D,INCFD,IERR) 
!***BEGIN PROLOGUE  PCHCE                                               
!***REFER TO  PCHIC                                                     
!***ROUTINES CALLED  PCHDF,PCHST,XERROR                                 
!***REVISION DATE  870707   (YYMMDD)                                    
!***DESCRIPTION                                                         
!                                                                       
!          PCHCE:  PCHIC End Derivative Setter.                         
!                                                                       
!    Called by PCHIC to set end derivatives as requested by the user.   
!    It must be called after interior derivative values have been set.  
!                      -----                                            
!                                                                       
!    To facilitate two-dimensional applications, includes an increment  
!    between successive values of the D-array.                          
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  Calling sequence:                                                    
!                                                                       
!        PARAMETER  (INCFD = ...)                                       
!        INTEGER  IC(2), N, IERR                                        
!        REAL  VC(2), X(N), H(N), SLOPE(N), D(INCFD,N)                  
!                                                                       
!        CALL  PCHCE (IC, VC, N, X, H, SLOPE, D, INCFD, IERR)           
!                                                                       
!   Parameters:                                                         
!                                                                       
!     IC -- (input) integer array of length 2 specifying desired        
!           boundary conditions:                                        
!           IC(1) = IBEG, desired condition at beginning of data.       
!           IC(2) = IEND, desired condition at end of data.             
!           ( see prologue to PCHIC for details. )                      
!                                                                       
!     VC -- (input) real array of length 2 specifying desired boundary  
!           values.  VC(1) need be set only if IC(1) = 2 or 3 .         
!                    VC(2) need be set only if IC(2) = 2 or 3 .         
!                                                                       
!     N -- (input) number of data points.  (assumes N.GE.2)             
!                                                                       
!     X -- (input) real array of independent variable values.  (the     
!           elements of X are assumed to be strictly increasing.)       
!                                                                       
!     H -- (input) real array of interval lengths.                      
!     SLOPE -- (input) real array of data slopes.                       
!           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!                  H(I) =  X(I+1)-X(I),                                 
!              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.               
!                                                                       
!     D -- (input) real array of derivative values at the data points.  
!           The value corresponding to X(I) must be stored in           
!                D(1+(I-1)*INCFD),  I=1(1)N.                            
!          (output) the value of D at X(1) and/or X(N) is changed, if   
!           necessary, to produce the requested boundary conditions.    
!           no other entries in D are changed.                          
!                                                                       
!     INCFD -- (input) increment between successive values in D.        
!           This argument is provided primarily for 2-D applications.   
!                                                                       
!     IERR -- (output) error flag.                                      
!           Normal return:                                              
!              IERR = 0  (no errors).                                   
!           Warning errors:                                             
!              IERR = 1  if IBEG.LT.0 and D(1) had to be adjusted for   
!                        monotonicity.                                  
!              IERR = 2  if IEND.LT.0 and D(1+(N-1)*INCFD) had to be    
!                        adjusted for monotonicity.                     
!              IERR = 3  if both of the above are true.                 
!                                                                       
!    -------                                                            
!    WARNING:  This routine does no validity-checking of arguments.     
!    -------                                                            
!                                                                       
!  Fortran intrinsics used:  ABS, IABS.                                 
!                                                                       
!***END PROLOGUE  PCHCE                                                 
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,      
!                  Mathematics and Statistics Division,                 
!                  Lawrence Livermore National Laboratory.              
!                                                                       
!  Change record:                                                       
!     82-08-05   Converted to SLATEC library version.                   
!     87-07-07   Minor corrections made to prologue..                   
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  Programming notes:                                                   
!                                                                       
!     1. The function  PCHST(ARG1,ARG2)  is assumed to return zero if   
!        either argument is zero, +1 if they are of the same sign, and  
!        -1 if they are of opposite sign.                               
!     2. To produce a double precision version, simply:                 
!        a. Change PCHCE to DPCHCE wherever it occurs,                  
!        b. Change PCHDF to DPCHCE wherever it occurs,                  
!        c. Change PCHST to DPCHST wherever it occurs,                  
!        d. Change all references to the Fortran intrinsics to their    
!           double presision equivalents,                               
!        e. Change the real declarations to double precision, and       
!        f. Change the constants ZERO, HALF, ... to double precision.   
!     3. One could reduce the number of arguments and amount of local   
!        storage, at the expense of reduced code clarity, by passing in 
!        the array WK (rather than splitting it into H and SLOPE) and   
!        increasing its length enough to incorporate STEMP and XTEMP.   
!     4. The two monotonicity checks only use the sufficient conditions.
!        Thus, it is possible (but unlikely) for a boundary condition to
!        be changed, even though the original interpolant was monotonic.
!        (At least the result is a continuous function of the data.)    
!                                                                       
!  DECLARE ARGUMENTS.                                                   
!                                                                       
      INTEGER  IC(2), N, INCFD, IERR 
      REAL  VC(2), X(N), H(N), SLOPE(N), D(INCFD,N) 
!                                                                       
!  DELCARE LOCAL VARIABLES.                                             
!                                                                       
      INTEGER  IBEG, IEND, IERF, INDEX, J, K 
      REAL  HALF, STEMP(3), THREE, TWO, XTEMP(4), ZERO 
      REAL  PCHDF, PCHST 
!                                                                       
!  INITIALIZE.                                                          
!                                                                       
      DATA  ZERO /0./,  HALF /0.5/,  TWO /2./,  THREE /3./ 
!                                                                       
!***FIRST EXECUTABLE STATEMENT  PCHCE                                   
      IBEG = IC(1) 
      IEND = IC(2) 
      IERR = 0 
!                                                                       
!  SET TO DEFAULT BOUNDARY CONDITIONS IF N IS TOO SMALL.                
!                                                                       
      IF ( IABS(IBEG).GT.N )  IBEG = 0 
      IF ( IABS(IEND).GT.N )  IEND = 0 
!                                                                       
!  TREAT BEGINNING BOUNDARY CONDITION.                                  
!                                                                       
      IF (IBEG .EQ. 0)  GO TO 2000 
      K = IABS(IBEG) 
      IF (K .EQ. 1)  THEN 
!        BOUNDARY VALUE PROVIDED.                                       
         D(1,1) = VC(1) 
      ELSE IF (K .EQ. 2)  THEN 
!        BOUNDARY SECOND DERIVATIVE PROVIDED.                           
         D(1,1) = HALF*( (THREE*SLOPE(1) - D(1,2)) - HALF*VC(1)*H(1) ) 
      ELSE IF (K .LT. 5)  THEN 
!        USE K-POINT DERIVATIVE FORMULA.                                
!        PICK UP FIRST K POINTS, IN REVERSE ORDER.                      
         DO 10  J = 1, K 
            INDEX = K-J+1 
!           INDEX RUNS FROM K DOWN TO 1.                                
            XTEMP(J) = X(INDEX) 
            IF (J .LT. K)  STEMP(J) = SLOPE(INDEX-1) 
   10    CONTINUE 
!                 -----------------------------                         
         D(1,1) = PCHDF (K, XTEMP, STEMP, IERF) 
!                 -----------------------------                         
         IF (IERF .NE. 0)  GO TO 5001 
      ELSE 
!        USE 'NOT A KNOT' CONDITION.                                    
         D(1,1) = ( THREE*(H(1)*SLOPE(2) + H(2)*SLOPE(1))               &
     &             - TWO*(H(1)+H(2))*D(1,2) - H(1)*D(1,3) ) / H(2)      
      ENDIF 
!                                                                       
      IF (IBEG .GT. 0)  GO TO 2000 
!                                                                       
!  CHECK D(1,1) FOR COMPATIBILITY WITH MONOTONICITY.                    
!                                                                       
      IF (SLOPE(1) .EQ. ZERO)  THEN 
         IF (D(1,1) .NE. ZERO)  THEN 
            D(1,1) = ZERO 
            IERR = IERR + 1 
         ENDIF 
      ELSE IF ( PCHST(D(1,1),SLOPE(1)) .LT. ZERO)  THEN 
         D(1,1) = ZERO 
         IERR = IERR + 1 
      ELSE IF ( ABS(D(1,1)) .GT. THREE*ABS(SLOPE(1)) )  THEN 
         D(1,1) = THREE*SLOPE(1) 
         IERR = IERR + 1 
      ENDIF 
!                                                                       
!  TREAT END BOUNDARY CONDITION.                                        
!                                                                       
 2000 CONTINUE 
      IF (IEND .EQ. 0)  GO TO 5000 
      K = IABS(IEND) 
      IF (K .EQ. 1)  THEN 
!        BOUNDARY VALUE PROVIDED.                                       
         D(1,N) = VC(2) 
      ELSE IF (K .EQ. 2)  THEN 
!        BOUNDARY SECOND DERIVATIVE PROVIDED.                           
         D(1,N) = HALF*( (THREE*SLOPE(N-1) - D(1,N-1)) +                &
     &                                           HALF*VC(2)*H(N-1) )    
      ELSE IF (K .LT. 5)  THEN 
!        USE K-POINT DERIVATIVE FORMULA.                                
!        PICK UP LAST K POINTS.                                         
         DO 2010  J = 1, K 
            INDEX = N-K+J 
!           INDEX RUNS FROM N+1-K UP TO N.                              
            XTEMP(J) = X(INDEX) 
            IF (J .LT. K)  STEMP(J) = SLOPE(INDEX) 
 2010    CONTINUE 
!                 -----------------------------                         
         D(1,N) = PCHDF (K, XTEMP, STEMP, IERF) 
!                 -----------------------------                         
         IF (IERF .NE. 0)  GO TO 5001 
      ELSE 
!        USE 'NOT A KNOT' CONDITION.                                    
         D(1,N) = ( THREE*(H(N-1)*SLOPE(N-2) + H(N-2)*SLOPE(N-1))       &
     &             - TWO*(H(N-1)+H(N-2))*D(1,N-1) - H(N-1)*D(1,N-2) )   &
     &                                                         / H(N-2) 
      ENDIF 
!                                                                       
      IF (IEND .GT. 0)  GO TO 5000 
!                                                                       
!  CHECK D(1,N) FOR COMPATIBILITY WITH MONOTONICITY.                    
!                                                                       
      IF (SLOPE(N-1) .EQ. ZERO)  THEN 
         IF (D(1,N) .NE. ZERO)  THEN 
            D(1,N) = ZERO 
            IERR = IERR + 2 
         ENDIF 
      ELSE IF ( PCHST(D(1,N),SLOPE(N-1)) .LT. ZERO)  THEN 
         D(1,N) = ZERO 
         IERR = IERR + 2 
      ELSE IF ( ABS(D(1,N)) .GT. THREE*ABS(SLOPE(N-1)) )  THEN 
         D(1,N) = THREE*SLOPE(N-1) 
         IERR = IERR + 2 
      ENDIF 
!                                                                       
!  NORMAL RETURN.                                                       
!                                                                       
 5000 CONTINUE 
      RETURN 
!                                                                       
!  ERROR RETURN.                                                        
!                                                                       
 5001 CONTINUE 
!     ERROR RETURN FROM PCHDF.                                          
!   *** THIS CASE SHOULD NEVER OCCUR ***                                
      IERR = -1 
      CALL XERROR ('PCHCE -- ERROR RETURN FROM PCHDF'                   &
     &           , 32, IERR, 1)                                         
      RETURN 
!------------- LAST LINE OF PCHCE FOLLOWS ------------------------------
      END                                           
      SUBROUTINE PCHCI(N,H,SLOPE,D,INCFD) 
!***BEGIN PROLOGUE  PCHCI                                               
!***REFER TO  PCHIC                                                     
!***ROUTINES CALLED  PCHST                                              
!***REVISION DATE  870707   (YYMMDD)                                    
!***DESCRIPTION                                                         
!                                                                       
!          PCHCI:  PCHIC Initial Derivative Setter.                     
!                                                                       
!    Called by PCHIC to set derivatives needed to determine a monotone  
!    piecewise cubic Hermite interpolant to the data.                   
!                                                                       
!    Default boundary conditions are provided which are compatible      
!    with monotonicity.  If the data are only piecewise monotonic, the  
!    interpolant will have an extremum at each point where monotonicity 
!    switches direction.                                                
!                                                                       
!    To facilitate two-dimensional applications, includes an increment  
!    between successive values of the D-array.                          
!                                                                       
!    The resulting piecewise cubic Hermite function should be identical 
!    (within roundoff error) to that produced by PCHIM.                 
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  Calling sequence:                                                    
!                                                                       
!        PARAMETER  (INCFD = ...)                                       
!        INTEGER  N                                                     
!        REAL  H(N), SLOPE(N), D(INCFD,N)                               
!                                                                       
!        CALL  PCHCI (N, H, SLOPE, D, INCFD)                            
!                                                                       
!   Parameters:                                                         
!                                                                       
!     N -- (input) number of data points.                               
!           If N=2, simply does linear interpolation.                   
!                                                                       
!     H -- (input) real array of interval lengths.                      
!     SLOPE -- (input) real array of data slopes.                       
!           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!                  H(I) =  X(I+1)-X(I),                                 
!              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.               
!                                                                       
!     D -- (output) real array of derivative values at the data points. 
!           If the data are monotonic, these values will determine a    
!           a monotone cubic Hermite function.                          
!           The value corresponding to X(I) is stored in                
!                D(1+(I-1)*INCFD),  I=1(1)N.                            
!           No other entries in D are changed.                          
!                                                                       
!     INCFD -- (input) increment between successive values in D.        
!           This argument is provided primarily for 2-D applications.   
!                                                                       
!    -------                                                            
!    WARNING:  This routine does no validity-checking of arguments.     
!    -------                                                            
!                                                                       
!  Fortran intrinsics used:  ABS, AMAX1, AMIN1.                         
!                                                                       
!***END PROLOGUE  PCHCI                                                 
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,      
!                  Mathematics and Statistics Division,                 
!                  Lawrence Livermore National Laboratory.              
!                                                                       
!  Change record:                                                       
!     82-06-01   Modified end conditions to be continuous functions of  
!                data when monotonicity switches in next interval.      
!     82-06-02   1. Modified formulas so end conditions are less prone  
!                   to over/underflow problems.                         
!                2. Minor modification to HSUM calculation.             
!     82-08-05   Converted to SLATEC library version.                   
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  Programming notes:                                                   
!                                                                       
!     1. The function  PCHST(ARG1,ARG2)  is assumed to return zero if   
!        either argument is zero, +1 if they are of the same sign, and  
!        -1 if they are of opposite sign.                               
!     2. To produce a double precision version, simply:                 
!        a. Change PCHCI to DPCHCI wherever it occurs,                  
!        b. Change PCHST to DPCHST wherever it occurs,                  
!        c. Change all references to the Fortran intrinsics to their    
!           double presision equivalents,                               
!        d. Change the real declarations to double precision, and       
!        e. Change the constants ZERO and THREE to double precision.    
!                                                                       
!  DECLARE ARGUMENTS.                                                   
!                                                                       
      INTEGER  N, INCFD 
      REAL  H(N), SLOPE(N), D(INCFD,N) 
!                                                                       
!  DECLARE LOCAL VARIBLES.                                              
!                                                                       
      INTEGER  I, NLESS1 
      REAL  DEL1, DEL2, DMAX, DMIN, DRAT1, DRAT2, HSUM, HSUMT3, THREE,  &
     &      W1, W2, ZERO                                                
      REAL  PCHST 
!                                                                       
!  INITIALIZE.                                                          
!                                                                       
      DATA  ZERO /0./,  THREE /3./ 
!***FIRST EXECUTABLE STATEMENT  PCHCI                                   
      NLESS1 = N - 1 
      DEL1 = SLOPE(1) 
!                                                                       
!  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.                        
!                                                                       
      IF (NLESS1 .GT. 1)  GO TO 10 
      D(1,1) = DEL1 
      D(1,N) = DEL1 
      GO TO 5000 
!                                                                       
!  NORMAL CASE  (N .GE. 3).                                             
!                                                                       
   10 CONTINUE 
      DEL2 = SLOPE(2) 
!                                                                       
!  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE        
!     SHAPE-PRESERVING.                                                 
!                                                                       
      HSUM = H(1) + H(2) 
      W1 = (H(1) + HSUM)/HSUM 
      W2 = -H(1)/HSUM 
      D(1,1) = W1*DEL1 + W2*DEL2 
      IF ( PCHST(D(1,1),DEL1) .LE. ZERO)  THEN 
         D(1,1) = ZERO 
      ELSE IF ( PCHST(DEL1,DEL2) .LT. ZERO)  THEN 
!        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.              
         DMAX = THREE*DEL1 
         IF (ABS(D(1,1)) .GT. ABS(DMAX))  D(1,1) = DMAX 
      ENDIF 
!                                                                       
!  LOOP THROUGH INTERIOR POINTS.                                        
!                                                                       
      DO 50  I = 2, NLESS1 
         IF (I .EQ. 2)  GO TO 40 
!                                                                       
         HSUM = H(I-1) + H(I) 
         DEL1 = DEL2 
         DEL2 = SLOPE(I) 
   40    CONTINUE 
!                                                                       
!        SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.                 
!                                                                       
         D(1,I) = ZERO 
         IF ( PCHST(DEL1,DEL2) .LE. ZERO)  GO TO 50 
!                                                                       
!        USE BRODLIE MODIFICATION OF BUTLAND FORMULA.                   
!                                                                       
         HSUMT3 = HSUM+HSUM+HSUM 
         W1 = (HSUM + H(I-1))/HSUMT3 
         W2 = (HSUM + H(I)  )/HSUMT3 
         DMAX = AMAX1( ABS(DEL1), ABS(DEL2) ) 
         DMIN = AMIN1( ABS(DEL1), ABS(DEL2) ) 
         DRAT1 = DEL1/DMAX 
         DRAT2 = DEL2/DMAX 
         D(1,I) = DMIN/(W1*DRAT1 + W2*DRAT2) 
!                                                                       
   50 END DO 
!                                                                       
!  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE        
!     SHAPE-PRESERVING.                                                 
!                                                                       
      W1 = -H(N-1)/HSUM 
      W2 = (H(N-1) + HSUM)/HSUM 
      D(1,N) = W1*DEL1 + W2*DEL2 
      IF ( PCHST(D(1,N),DEL2) .LE. ZERO)  THEN 
         D(1,N) = ZERO 
      ELSE IF ( PCHST(DEL1,DEL2) .LT. ZERO)  THEN 
!        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.              
         DMAX = THREE*DEL2 
         IF (ABS(D(1,N)) .GT. ABS(DMAX))  D(1,N) = DMAX 
      ENDIF 
!                                                                       
!  NORMAL RETURN.                                                       
!                                                                       
 5000 CONTINUE 
      RETURN 
!------------- LAST LINE OF PCHCI FOLLOWS ------------------------------
      END                                           
      SUBROUTINE PCHCS(SWITCH,N,H,SLOPE,D,INCFD,IERR) 
!***BEGIN PROLOGUE  PCHCS                                               
!***REFER TO  PCHIC                                                     
!***ROUTINES CALLED  PCHST,PCHSW                                        
!***REVISION DATE  870707   (YYMMDD)                                    
!***DESCRIPTION                                                         
!                                                                       
!         PCHCS:  PCHIC Monotonicity Switch Derivative Setter.          
!                                                                       
!     Called by  PCHIC  to adjust the values of D in the vicinity of a  
!     switch in direction of monotonicity, to produce a more "visually  
!     pleasing" curve than that given by  PCHIM .                       
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  Calling sequence:                                                    
!                                                                       
!        PARAMETER  (INCFD = ...)                                       
!        INTEGER  N, IERR                                               
!        REAL  SWITCH, H(N), SLOPE(N), D(INCFD,N)                       
!                                                                       
!        CALL  PCHCS (SWITCH, N, H, SLOPE, D, INCFD, IERR)              
!                                                                       
!   Parameters:                                                         
!                                                                       
!     SWITCH -- (input) indicates the amount of control desired over    
!           local excursions from data.                                 
!                                                                       
!     N -- (input) number of data points.  (assumes N.GT.2 .)           
!                                                                       
!     H -- (input) real array of interval lengths.                      
!     SLOPE -- (input) real array of data slopes.                       
!           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!                  H(I) =  X(I+1)-X(I),                                 
!              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.               
!                                                                       
!     D -- (input) real array of derivative values at the data points,  
!           as determined by PCHCI.                                     
!          (output) derivatives in the vicinity of switches in direction
!           of monotonicity may be adjusted to produce a more "visually 
!           pleasing" curve.                                            
!           The value corresponding to X(I) is stored in                
!                D(1+(I-1)*INCFD),  I=1(1)N.                            
!           No other entries in D are changed.                          
!                                                                       
!     INCFD -- (input) increment between successive values in D.        
!           This argument is provided primarily for 2-D applications.   
!                                                                       
!     IERR -- (output) error flag.  should be zero.                     
!           If negative, trouble in PCHSW.  (should never happen.)      
!                                                                       
!    -------                                                            
!    WARNING:  This routine does no validity-checking of arguments.     
!    -------                                                            
!                                                                       
!  Fortran intrinsics used:  ABS, AMAX1, AMIN1.                         
!                                                                       
!***END PROLOGUE  PCHCS                                                 
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,      
!                  Mathematics and Statistics Division,                 
!                  Lawrence Livermore National Laboratory.              
!                                                                       
!  Change record:                                                       
!     82-06-17   Redesigned to (1) fix  problem with lack of continuity 
!                approaching a flat-topped peak (2) be cleaner and      
!                easier to verify.                                      
!                Eliminated subroutines PCHSA and PCHSX in the process. 
!     82-06-22   1. Limited fact to not exceed one, so computed D is a  
!                   convex combination of PCHCI value and PCHSD value.  
!                2. Changed fudge from 1 to 4 (based on experiments).   
!     82-06-23   Moved PCHSD to an inline function (eliminating MSWTYP).
!     82-08-05   Converted to SLATEC library version.                   
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  Programming notes:                                                   
!                                                                       
!     1. The function  PCHST(ARG1,ARG2)  is assumed to return zero if   
!        either argument is zero, +1 if they are of the same sign, and  
!        -1 if they are of opposite sign.                               
!     2. To produce a double precision version, simply:                 
!        a. Change PCHCS to DPCHCS wherever it occurs,                  
!        b. Change PCHSD to DPCHSD wherever it occurs,                  
!        c. Change PCHST to DPCHST wherever it occurs,                  
!        d. Change PCHSW to DPCHSW wherever it occurs,                  
!        e. Change all references to the Fortran intrinsics to their    
!           double precision equivalents,                               
!        f. Change the real declarations to double precision, and       
!        g. Change the constants ZERO and ONE to double precision.      
!                                                                       
!  DECLARE ARGUMENTS.                                                   
!                                                                       
      INTEGER  N, INCFD, IERR 
      REAL  SWITCH, H(N), SLOPE(N), D(INCFD,N) 
!                                                                       
!  DECLARE LOCAL VARIABLES.                                             
!                                                                       
      INTEGER  I, INDX, K, NLESS1 
      REAL  DEL(3), DEXT, DFLOC, DFMX, FACT, FUDGE, ONE, SLMAX,         &
     &      WTAVE(2), ZERO                                              
      REAL  PCHST 
!                                                                       
!  DEFINE INLINE FUNCTION FOR WEIGHTED AVERAGE OF SLOPES.               
!                                                                       
      REAL  PCHSD, S1, S2, H1, H2 
      PCHSD(S1,S2,H1,H2) = (H2/(H1+H2))*S1 + (H1/(H1+H2))*S2 
!                                                                       
!  INITIALIZE.                                                          
!                                                                       
      DATA  ZERO /0./,  ONE /1./ 
      DATA  FUDGE /4./ 
!***FIRST EXECUTABLE STATEMENT  PCHCS                                   
      IERR = 0 
      NLESS1 = N - 1 
!                                                                       
!  LOOP OVER SEGMENTS.                                                  
!                                                                       
      DO 900  I = 2, NLESS1 
         IF ( PCHST(SLOPE(I-1),SLOPE(I)) )  100, 300, 900 
!             --------------------------                                
!                                                                       
  100    CONTINUE 
!                                                                       
!....... SLOPE SWITCHES MONOTONICITY AT I-TH POINT .....................
!                                                                       
!           DO NOT CHANGE D IF 'UP-DOWN-UP'.                            
            IF (I .GT. 2)  THEN 
               IF ( PCHST(SLOPE(I-2),SLOPE(I)) .GT. ZERO)  GO TO 900 
!                   --------------------------                          
            ENDIF 
            IF (I .LT. NLESS1)  THEN 
               IF ( PCHST(SLOPE(I+1),SLOPE(I-1)) .GT. ZERO)  GO TO 900 
!                   ----------------------------                        
            ENDIF 
!                                                                       
!   ....... COMPUTE PROVISIONAL VALUE FOR D(1,I).                       
!                                                                       
            DEXT = PCHSD (SLOPE(I-1), SLOPE(I), H(I-1), H(I)) 
!                                                                       
!   ....... DETERMINE WHICH INTERVAL CONTAINS THE EXTREMUM.             
!                                                                       
            IF ( PCHST(DEXT, SLOPE(I-1)) )  200, 900, 250 
!                -----------------------                                
!                                                                       
  200       CONTINUE 
!              DEXT AND SLOPE(I-1) HAVE OPPOSITE SIGNS --               
!                        EXTREMUM IS IN (X(I-1),X(I)).                  
               K = I-1 
!              SET UP TO COMPUTE NEW VALUES FOR D(1,I-1) AND D(1,I).    
               WTAVE(2) = DEXT 
               IF (K .GT. 1)                                            &
     &            WTAVE(1) = PCHSD (SLOPE(K-1), SLOPE(K), H(K-1), H(K)) 
               GO TO 400 
!                                                                       
  250       CONTINUE 
!              DEXT AND SLOPE(I) HAVE OPPOSITE SIGNS --                 
!                        EXTREMUM IS IN (X(I),X(I+1)).                  
               K = I 
!              SET UP TO COMPUTE NEW VALUES FOR D(1,I) AND D(1,I+1).    
               WTAVE(1) = DEXT 
               IF (K .LT. NLESS1)                                       &
     &            WTAVE(2) = PCHSD (SLOPE(K), SLOPE(K+1), H(K), H(K+1)) 
               GO TO 400 
!                                                                       
  300    CONTINUE 
!                                                                       
!....... AT LEAST ONE OF SLOPE(I-1) AND SLOPE(I) IS ZERO --             
!                     CHECK FOR FLAT-TOPPED PEAK .......................
!                                                                       
            IF (I .EQ. NLESS1)  GO TO 900 
            IF ( PCHST(SLOPE(I-1), SLOPE(I+1)) .GE. ZERO)  GO TO 900 
!                -----------------------------                          
!                                                                       
!           WE HAVE FLAT-TOPPED PEAK ON (X(I),X(I+1)).                  
            K = I 
!           SET UP TO COMPUTE NEW VALUES FOR D(1,I) AND D(1,I+1).       
            WTAVE(1) = PCHSD (SLOPE(K-1), SLOPE(K), H(K-1), H(K)) 
            WTAVE(2) = PCHSD (SLOPE(K), SLOPE(K+1), H(K), H(K+1)) 
!                                                                       
  400    CONTINUE 
!                                                                       
!....... AT THIS POINT WE HAVE DETERMINED THAT THERE WILL BE AN EXTREMUM
!        ON (X(K),X(K+1)), WHERE K=I OR I-1, AND HAVE SET ARRAY WTAVE-- 
!           WTAVE(1) IS A WEIGHTED AVERAGE OF SLOPE(K-1) AND SLOPE(K),  
!                    IF K.GT.1                                          
!           WTAVE(2) IS A WEIGHTED AVERAGE OF SLOPE(K) AND SLOPE(K+1),  
!                    IF K.LT.N-1                                        
!                                                                       
         SLMAX = ABS(SLOPE(K)) 
         IF (K .GT. 1)    SLMAX = AMAX1( SLMAX, ABS(SLOPE(K-1)) ) 
         IF (K.LT.NLESS1) SLMAX = AMAX1( SLMAX, ABS(SLOPE(K+1)) ) 
!                                                                       
         IF (K .GT. 1)  DEL(1) = SLOPE(K-1) / SLMAX 
         DEL(2) = SLOPE(K) / SLMAX 
         IF (K.LT.NLESS1)  DEL(3) = SLOPE(K+1) / SLMAX 
!                                                                       
         IF ((K.GT.1) .AND. (K.LT.NLESS1))  THEN 
!           NORMAL CASE -- EXTREMUM IS NOT IN A BOUNDARY INTERVAL.      
            FACT = FUDGE* ABS(DEL(3)*(DEL(1)-DEL(2))*(WTAVE(2)/SLMAX)) 
            D(1,K) = D(1,K) + AMIN1(FACT,ONE)*(WTAVE(1) - D(1,K)) 
            FACT = FUDGE* ABS(DEL(1)*(DEL(3)-DEL(2))*(WTAVE(1)/SLMAX)) 
            D(1,K+1) = D(1,K+1) + AMIN1(FACT,ONE)*(WTAVE(2) - D(1,K+1)) 
         ELSE 
!           SPECIAL CASE K=1 (WHICH CAN OCCUR ONLY IF I=2) OR           
!                        K=NLESS1 (WHICH CAN OCCUR ONLY IF I=NLESS1).   
            FACT = FUDGE* ABS(DEL(2)) 
            D(1,I) = AMIN1(FACT,ONE) * WTAVE(I-K+1) 
!              NOTE THAT I-K+1 = 1 IF K=I  (=NLESS1),                   
!                        I-K+1 = 2 IF K=I-1(=1).                        
         ENDIF 
!                                                                       
!                                                                       
!....... ADJUST IF NECESSARY TO LIMIT EXCURSIONS FROM DATA.             
!                                                                       
         IF (SWITCH .LE. ZERO)  GO TO 900 
!                                                                       
         DFLOC = H(K)*ABS(SLOPE(K)) 
         IF (K .GT. 1)    DFLOC = AMAX1( DFLOC, H(K-1)*ABS(SLOPE(K-1)) ) 
         IF (K.LT.NLESS1) DFLOC = AMAX1( DFLOC, H(K+1)*ABS(SLOPE(K+1)) ) 
         DFMX = SWITCH*DFLOC 
         INDX = I-K+1 
!        INDX = 1 IF K=I, 2 IF K=I-1.                                   
!        ---------------------------------------------------------------
         CALL PCHSW (DFMX, INDX, D(1,K), D(1,K+1), H(K), SLOPE(K), IERR) 
!        ---------------------------------------------------------------
         IF (IERR .NE. 0)  RETURN 
!                                                                       
!....... END OF SEGMENT LOOP.                                           
!                                                                       
  900 END DO 
!                                                                       
      RETURN 
!------------- LAST LINE OF PCHCS FOLLOWS ------------------------------
      END                                           
      REAL FUNCTION PCHDF(K,X,S,IERR) 
!***BEGIN PROLOGUE  PCHDF                                               
!***REFER TO  PCHCE,PCHSP                                               
!***ROUTINES CALLED  XERROR                                             
!***REVISION DATE  870707   (YYMMDD)                                    
!***DESCRIPTION                                                         
!                                                                       
!          PCHDF:   PCHIP Finite Difference Formula                     
!                                                                       
!     Uses a divided difference formulation to compute a K-point approx-
!     imation to the derivative at X(K) based on the data in X and S.   
!                                                                       
!     Called by  PCHCE  and  PCHSP  to compute 3- and 4-point boundary  
!     derivative approximations.                                        
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!     On input:                                                         
!        K      is the order of the desired derivative approximation.   
!               K must be at least 3 (error return if not).             
!        X      contains the K values of the independent variable.      
!               X need not be ordered, but the values **MUST** be       
!               distinct.  (Not checked here.)                          
!        S      contains the associated slope values:                   
!                  S(I) = (F(I+1)-F(I))/(X(I+1)-X(I)), I=1(1)K-1.       
!               (Note that S need only be of length K-1.)               
!                                                                       
!     On return:                                                        
!        S      will be destroyed.                                      
!        IERR   will be set to -1 if K.LT.2 .                           
!        PCHDF  will be set to the desired derivative approximation if  
!               IERR=0 or to zero if IERR=-1.                           
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  Reference:  Carl de Boor, A Practical Guide to Splines, Springer-    
!              Verlag (New York, 1978), pp. 10-16.                      
!                                                                       
!***END PROLOGUE  PCHDF                                                 
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,      
!                  Mathematics and Statistics Division,                 
!                  Lawrence Livermore National Laboratory.              
!                                                                       
!  Change record:                                                       
!     82-08-05   Converted to SLATEC library version.                   
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  Programming notes:                                                   
!                                                                       
!     To produce a double precision version, simply:                    
!        a. Change PCHDF to DPCHDF wherever it occurs,                  
!        b. Change the real declarations to double precision, and       
!        c. Change the constant ZERO to double precision.               
      INTEGER  K, IERR 
      REAL  X(K), S(K) 
!                                                                       
!  DECLARE LOCAL VARIABLES.                                             
!                                                                       
      INTEGER  I, J 
      REAL  VALUE, ZERO 
      DATA  ZERO /0./ 
!                                                                       
!  CHECK FOR LEGAL VALUE OF K.                                          
!                                                                       
!***FIRST EXECUTABLE STATEMENT  PCHDF                                   
      IF (K .LT. 3)  GO TO 5001 
!                                                                       
!  COMPUTE COEFFICIENTS OF INTERPOLATING POLYNOMIAL.                    
!                                                                       
      DO 10  J = 2, K-1 
         DO 9  I = 1, K-J 
            S(I) = (S(I+1)-S(I))/(X(I+J)-X(I)) 
    9    CONTINUE 
   10 END DO 
!                                                                       
!  EVALUATE DERIVATIVE AT X(K).                                         
!                                                                       
      VALUE = S(1) 
      DO 20  I = 2, K-1 
         VALUE = S(I) + VALUE*(X(K)-X(I)) 
   20 END DO 
!                                                                       
!  NORMAL RETURN.                                                       
!                                                                       
      IERR = 0 
      PCHDF = VALUE 
      RETURN 
!                                                                       
!  ERROR RETURN.                                                        
!                                                                       
 5001 CONTINUE 
!     K.LT.3 RETURN.                                                    
      IERR = -1 
      CALL XERROR ('PCHDF -- K LESS THAN THREE'                         &
     &           , 26, IERR, 1)                                         
      PCHDF = ZERO 
      RETURN 
!------------- LAST LINE OF PCHDF FOLLOWS ------------------------------
      END                                           
      REAL FUNCTION PCHST(ARG1,ARG2) 
!***BEGIN PROLOGUE  PCHST                                               
!***REFER TO  PCHCE,PCHCI,PCHCS,PCHIM                                   
!***ROUTINES CALLED  (NONE)                                             
!***REVISION DATE  870707   (YYMMDD)                                    
!***DESCRIPTION                                                         
!                                                                       
!         PCHST:  PCHIP Sign-Testing Routine.                           
!                                                                       
!                                                                       
!     Returns:                                                          
!        -1. if ARG1 and ARG2 are of opposite sign.                     
!         0. if either argument is zero.                                
!        +1. if ARG1 and ARG2 are of the same sign.                     
!                                                                       
!     The object is to do this without multiplying ARG1*ARG2, to avoid  
!     possible over/underflow problems.                                 
!                                                                       
!  Fortran intrinsics used:  SIGN.                                      
!                                                                       
!***END PROLOGUE  PCHST                                                 
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,      
!                  Mathematics and Statistics Division,                 
!                  Lawrence Livermore National Laboratory.              
!                                                                       
!  Change record:                                                       
!     82-08-05   Converted to SLATEC library version.                   
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  Programming notes:                                                   
!                                                                       
!     To produce a double precision version, simply:                    
!        a. Change PCHST to DPCHST wherever it occurs,                  
!        b. Change all references to the Fortran intrinsics to their    
!           double presision equivalents,                               
!        c. Change the real declarations to double precision, and       
!        d. Change the constants  ZERO  and  ONE  to double precision.  
!                                                                       
!  DECLARE ARGUMENTS.                                                   
!                                                                       
      REAL  ARG1, ARG2 
!                                                                       
!  DECLARE LOCAL VARIABLES.                                             
!                                                                       
      REAL  ONE, ZERO 
      DATA  ZERO /0./,  ONE /1./ 
!                                                                       
!  PERFORM THE TEST.                                                    
!                                                                       
!***FIRST EXECUTABLE STATEMENT  PCHST                                   
      PCHST = SIGN(ONE,ARG1) * SIGN(ONE,ARG2) 
      IF ((ARG1.EQ.ZERO) .OR. (ARG2.EQ.ZERO))  PCHST = ZERO 
!                                                                       
      RETURN 
!------------- LAST LINE OF PCHST FOLLOWS ------------------------------
      END                                           
      SUBROUTINE PCHSW(DFMAX,IEXTRM,D1,D2,H,SLOPE,IERR) 
!***BEGIN PROLOGUE  PCHSW                                               
!***REFER TO  PCHCS                                                     
!***ROUTINES CALLED  R1MACH,XERROR                                      
!***REVISION DATE  870707   (YYMMDD)                                    
!***DESCRIPTION                                                         
!                                                                       
!         PCHSW:  PCHCS Switch Excursion Limiter.                       
!                                                                       
!     Called by  PCHCS  to adjust D1 and D2 if necessary to insure that 
!     the extremum on this interval is not further than DFMAX from the  
!     extreme data value.                                               
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  Calling sequence:                                                    
!                                                                       
!        INTEGER  IEXTRM, IERR                                          
!        REAL  DFMAX, D1, D2, H, SLOPE                                  
!                                                                       
!        CALL  PCHSW (DFMAX, IEXTRM, D1, D2, H, SLOPE, IERR)            
!                                                                       
!   Parameters:                                                         
!                                                                       
!     DFMAX -- (input) maximum allowed difference between F(IEXTRM) and 
!           the cubic determined by derivative values D1,D2.  (assumes  
!           DFMAX.GT.0.)                                                
!                                                                       
!     IEXTRM -- (input) index of the extreme data value.  (assumes      
!           IEXTRM = 1 or 2 .  Any value .NE.1 is treated as 2.)        
!                                                                       
!     D1,D2 -- (input) derivative values at the ends of the interval.   
!           (Assumes D1*D2 .LE. 0.)                                     
!          (output) may be modified if necessary to meet the restriction
!           imposed by DFMAX.                                           
!                                                                       
!     H -- (input) interval length.  (Assumes  H.GT.0.)                 
!                                                                       
!     SLOPE -- (input) data SLOPE on the interval.                      
!                                                                       
!     IERR -- (output) error flag.  should be zero.                     
!           If IERR=-1, assumption on D1 and D2 is not satisfied.       
!           If IERR=-2, quadratic equation locating extremum has        
!                       negative descriminant (should never occur).     
!                                                                       
!    -------                                                            
!    WARNING:  This routine does no validity-checking of arguments.     
!    -------                                                            
!                                                                       
!                                                                       
!  Fortran intrinsics used:  ABS, SIGN, SQRT.                           
!                                                                       
!***END PROLOGUE  PCHSW                                                 
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,      
!                  Mathematics and Statistics Division,                 
!                  Lawrence Livermore National Laboratory.              
!                                                                       
!  Change record:                                                       
!     82-08-05   Converted to SLATEC library version.                   
!     87-07-07   Replaced DATA statement for SMALL with a use of R1MACH.
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  Programming notes:                                                   
!                                                                       
!     To produce a double precision version, simply:                    
!        a. Change PCHSW to DPCHSW wherever it occurs,                  
!        b. Change PCHCS to DPCHCS wherever it occurs,                  
!        c. Change R1MACH to D1MACH wherever it occurs,                 
!        d. Change all references to the Fortran intrinsics to their    
!           double precision equivalents,                               
!        e. Change the real declarations to double precision, and       
!        f. Change constants ZERO, ONE, TWO, ... to double precision.   
!                                                                       
!  DECLARE ARGUMENTS.                                                   
!                                                                       
      INTEGER  IEXTRM, IERR 
      REAL  DFMAX, D1, D2, H, SLOPE 
!                                                                       
!  DECLARE LOCAL VARIABLES.                                             
!                                                                       
      REAL  CP, DMAX, FACT, LAMBDA, NU, ONE, PHI, RADCAL, RHO, SIGMA,   &
     &      SMALL, THAT, THIRD, THREE, TWO, ZERO                        
      DATA  ZERO /0./,  ONE /1./,  TWO /2./,  THREE /3./, FACT /100./ 
!        THIRD SHOULD BE SLIGHTLY LESS THAN 1/3.                        
      DATA  THIRD /0.33333/ 
!                                                                       
!  NOTATION AND GENERAL REMARKS.                                        
!                                                                       
!     RHO IS THE RATIO OF THE DATA SLOPE TO THE DERIVATIVE BEING TESTED.
!     LAMBDA IS THE RATIO OF D2 TO D1.                                  
!     THAT = T-HAT(RHO) IS THE NORMALIZED LOCATION OF THE EXTREMUM.     
!     PHI IS THE NORMALIZED VALUE OF P(X)-F1 AT X = XHAT = X-HAT(RHO),  
!           WHERE  THAT = (XHAT - X1)/H .                               
!        THAT IS, P(XHAT)-F1 = D*H*PHI,  WHERE D=D1 OR D2.              
!     SIMILARLY,  P(XHAT)-F2 = D*H*(PHI-RHO) .                          
!                                                                       
!      SMALL SHOULD BE A FEW ORDERS OF MAGNITUDE GREATER THAN MACHEPS.  
!***FIRST EXECUTABLE STATEMENT  PCHSW                                   
      SMALL = FACT*R1MACH(4) 
!                                                                       
!  DO MAIN CALCULATION.                                                 
!                                                                       
      IF (D1 .EQ. ZERO)  THEN 
!                                                                       
!        SPECIAL CASE -- D1.EQ.ZERO .                                   
!                                                                       
!          IF D2 IS ALSO ZERO, THIS ROUTINE SHOULD NOT HAVE BEEN CALLED.
         IF (D2 .EQ. ZERO)  GO TO 5001 
!                                                                       
         RHO = SLOPE/D2 
!          EXTREMUM IS OUTSIDE INTERVAL WHEN RHO .GE. 1/3 .             
         IF (RHO .GE. THIRD)  GO TO 5000 
         THAT = (TWO*(THREE*RHO-ONE)) / (THREE*(TWO*RHO-ONE)) 
         PHI = THAT**2 * ((THREE*RHO-ONE)/THREE) 
!                                                                       
!          CONVERT TO DISTANCE FROM F2 IF IEXTRM.NE.1 .                 
         IF (IEXTRM .NE. 1)  PHI = PHI - RHO 
!                                                                       
!          TEST FOR EXCEEDING LIMIT, AND ADJUST ACCORDINGLY.            
         DMAX = DFMAX / (H*ABS(PHI)) 
         IF (ABS(D2) .GT. DMAX)  D2 = SIGN (DMAX, D2) 
      ELSE 
!                                                                       
         RHO = SLOPE/D1 
         LAMBDA = -D2/D1 
         IF (D2 .EQ. ZERO)  THEN 
!                                                                       
!           SPECIAL CASE -- D2.EQ.ZERO .                                
!                                                                       
!             EXTREMUM IS OUTSIDE INTERVAL WHEN RHO .GE. 1/3 .          
            IF (RHO .GE. THIRD)  GO TO 5000 
            CP = TWO - THREE*RHO 
            NU = ONE - TWO*RHO 
            THAT = ONE / (THREE*NU) 
         ELSE 
            IF (LAMBDA .LE. ZERO)  GO TO 5001 
!                                                                       
!           NORMAL CASE -- D1 AND D2 BOTH NONZERO, OPPOSITE SIGNS.      
!                                                                       
            NU = ONE - LAMBDA - TWO*RHO 
            SIGMA = ONE - RHO 
            CP = NU + SIGMA 
            IF (ABS(NU) .GT. SMALL)  THEN 
               RADCAL = (NU - (TWO*RHO+ONE))*NU + SIGMA**2 
               IF (RADCAL .LT. ZERO)  GO TO 5002 
               THAT = (CP - SQRT(RADCAL)) / (THREE*NU) 
            ELSE 
               THAT = ONE/(TWO*SIGMA) 
            ENDIF 
         ENDIF 
         PHI = THAT*((NU*THAT - CP)*THAT + ONE) 
!                                                                       
!          CONVERT TO DISTANCE FROM F2 IF IEXTRM.NE.1 .                 
         IF (IEXTRM .NE. 1)  PHI = PHI - RHO 
!                                                                       
!          TEST FOR EXCEEDING LIMIT, AND ADJUST ACCORDINGLY.            
         DMAX = DFMAX / (H*ABS(PHI)) 
         IF (ABS(D1) .GT. DMAX)  THEN 
            D1 = SIGN (DMAX, D1) 
            D2 = -LAMBDA*D1 
         ENDIF 
      ENDIF 
!                                                                       
!  NORMAL RETURN.                                                       
!                                                                       
 5000 CONTINUE 
      IERR = 0 
      RETURN 
!                                                                       
!  ERROR RETURNS.                                                       
!                                                                       
 5001 CONTINUE 
!     D1 AND D2 BOTH ZERO, OR BOTH NONZERO AND SAME SIGN.               
      IERR = -1 
      CALL XERROR ('PCHSW -- D1 AND/OR D2 INVALID'                      &
     &           , 29, IERR, 1)                                         
      RETURN 
!                                                                       
 5002 CONTINUE 
!     NEGATIVE VALUE OF RADICAL (SHOULD NEVER OCCUR).                   
      IERR = -2 
      CALL XERROR ('PCHSW -- NEGATIVE RADICAL'                          &
     &           , 25, IERR, 1)                                         
      RETURN 
!------------- LAST LINE OF PCHSW FOLLOWS ------------------------------
      END                                           
