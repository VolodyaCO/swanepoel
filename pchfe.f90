      SUBROUTINE PCHFE(N,X,F,D,INCFD,SKIP,NE,XE,FE,IERR) 
!***BEGIN PROLOGUE  PCHFE                                               
!***DATE WRITTEN   811020   (YYMMDD)                                    
!***REVISION DATE  870707   (YYMMDD)                                    
!***CATEGORY NO.  E3                                                    
!***KEYWORDS  LIBRARY=SLATEC(PCHIP),                                    
!             TYPE=SINGLE PRECISION(PCHFE-S DPCHFE-D),                  
!             CUBIC HERMITE EVALUATION,HERMITE INTERPOLATION,           
!             PIECEWISE CUBIC EVALUATION                                
!***AUTHOR  FRITSCH, F. N., (LLNL)                                      
!             MATHEMATICS AND STATISTICS DIVISION                       
!             LAWRENCE LIVERMORE NATIONAL LABORATORY                    
!             P.O. BOX 808  (L-316)                                     
!             LIVERMORE, CA  94550                                      
!             FTS 532-4275, (415) 422-4275                              
!***PURPOSE  EVALUATE A PIECEWISE CUBIC HERMITE FUNCTION AT AN ARRAY OF 
!            POINTS.  MAY BE USED BY ITSELF FOR HERMITE INTERPOLATION,  
!            OR AS AN EVALUATOR FOR PCHIM OR PCHIC.                     
!***DESCRIPTION                                                         
!                                                                       
!          PCHFE:  PIECEWISE CUBIC HERMITE FUNCTION EVALUATOR           
!                                                                       
!     EVALUATES THE CUBIC HERMITE FUNCTION DEFINED BY  N, X, F, D  AT   
!     THE POINTS  XE(J), J=1(1)NE.                                      
!                                                                       
!     TO PROVIDE COMPATIBILITY WITH PCHIM AND PCHIC, INCLUDES AN        
!     INCREMENT BETWEEN SUCCESSIVE VALUES OF THE F- AND D-ARRAYS.       
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  CALLING SEQUENCE:                                                    
!                                                                       
!        PARAMETER  (INCFD = ...)                                       
!        INTEGER  N, NE, IERR                                           
!        REAL  X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE)             
!        LOGICAL  SKIP                                                  
!                                                                       
!        CALL  PCHFE (N, X, F, D, INCFD, SKIP, NE, XE, FE, IERR)        
!                                                                       
!   PARAMETERS:                                                         
!                                                                       
!     N -- (INPUT) NUMBER OF DATA POINTS.  (ERROR RETURN IF N.LT.2 .)   
!                                                                       
!     X -- (INPUT) REAL ARRAY OF INDEPENDENT VARIABLE VALUES.  THE      
!           ELEMENTS OF X MUST BE STRICTLY INCREASING:                  
!                X(I-1) .LT. X(I),  I = 2(1)N.                          
!           (ERROR RETURN IF NOT.)                                      
!                                                                       
!     F -- (INPUT) REAL ARRAY OF FUNCTION VALUES.  F(1+(I-1)*INCFD) IS  
!           THE VALUE CORRESPONDING TO X(I).                            
!                                                                       
!     D -- (INPUT) REAL ARRAY OF DERIVATIVE VALUES.  D(1+(I-1)*INCFD) IS
!           THE VALUE CORRESPONDING TO X(I).                            
!                                                                       
!     INCFD -- (INPUT) INCREMENT BETWEEN SUCCESSIVE VALUES IN F AND D.  
!           (ERROR RETURN IF  INCFD.LT.1 .)                             
!                                                                       
!     SKIP -- (INPUT/OUTPUT) LOGICAL VARIABLE WHICH SHOULD BE SET TO    
!           .TRUE. IF THE USER WISHES TO SKIP CHECKS FOR VALIDITY OF    
!           PRECEDING PARAMETERS, OR TO .FALSE. OTHERWISE.              
!           THIS WILL SAVE TIME IN CASE THESE CHECKS HAVE ALREADY       
!           BEEN PERFORMED (SAY, IN PCHIM OR PCHIC).                    
!           SKIP WILL BE SET TO .TRUE. ON NORMAL RETURN.                
!                                                                       
!     NE -- (INPUT) NUMBER OF EVALUATION POINTS.  (ERROR RETURN IF      
!           NE.LT.1 .)                                                  
!                                                                       
!     XE -- (INPUT) REAL ARRAY OF POINTS AT WHICH THE FUNCTION IS TO BE 
!           EVALUATED.                                                  
!                                                                       
!          NOTES:                                                       
!           1. THE EVALUATION WILL BE MOST EFFICIENT IF THE ELEMENTS    
!              OF XE ARE INCREASING RELATIVE TO X;                      
!              THAT IS,   XE(J) .GE. X(I)                               
!              IMPLIES    XE(K) .GE. X(I),  ALL K.GE.J .                
!           2. IF ANY OF THE XE ARE OUTSIDE THE INTERVAL [X(1),X(N)],   
!              VALUES ARE EXTRAPOLATED FROM THE NEAREST EXTREME CUBIC,  
!              AND A WARNING ERROR IS RETURNED.                         
!                                                                       
!     FE -- (OUTPUT) REAL ARRAY OF VALUES OF THE CUBIC HERMITE FUNCTION 
!           DEFINED BY  N, X, F, D  AT THE POINTS  XE.                  
!                                                                       
!     IERR -- (OUTPUT) ERROR FLAG.                                      
!           NORMAL RETURN:                                              
!              IERR = 0  (NO ERRORS).                                   
!           WARNING ERROR:                                              
!              IERR.GT.0  MEANS THAT EXTRAPOLATION WAS PERFORMED AT     
!                 IERR POINTS.                                          
!           "RECOVERABLE" ERRORS:                                       
!              IERR = -1  IF N.LT.2 .                                   
!              IERR = -2  IF INCFD.LT.1 .                               
!              IERR = -3  IF THE X-ARRAY IS NOT STRICTLY INCREASING.    
!              IERR = -4  IF NE.LT.1 .                                  
!             (THE FE-ARRAY HAS NOT BEEN CHANGED IN ANY OF THESE CASES.)
!               NOTE:  THE ABOVE ERRORS ARE CHECKED IN THE ORDER LISTED,
!                   AND FOLLOWING ARGUMENTS HAVE **NOT** BEEN VALIDATED.
!                                                                       
!***REFERENCES  (NONE)                                                  
!***ROUTINES CALLED  CHFEV,XERROR                                       
!***END PROLOGUE  PCHFE                                                 
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  CHANGE RECORD:                                                       
!     82-08-03   MINOR COSMETIC CHANGES FOR RELEASE 1.                  
!     87-07-07   MINOR COSMETIC CHANGES TO PROLOGUE.                    
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  PROGRAMMING NOTES:                                                   
!                                                                       
!     1. TO PRODUCE A DOUBLE PRECISION VERSION, SIMPLY:                 
!        A. CHANGE PCHFE TO DPCHFE, AND CHFEV TO DCHFEV, WHEREVER THEY  
!           OCCUR,                                                      
!        B. CHANGE THE REAL DECLARATION TO DOUBLE PRECISION,            
!                                                                       
!     2. MOST OF THE CODING BETWEEN THE CALL TO CHFEV AND THE END OF    
!        THE IR-LOOP COULD BE ELIMINATED IF IT WERE PERMISSIBLE TO      
!        ASSUME THAT XE IS ORDERED RELATIVE TO X.                       
!                                                                       
!     3. CHFEV DOES NOT ASSUME THAT X1 IS LESS THAN X2.  THUS, IT WOULD 
!        BE POSSIBLE TO WRITE A VERSION OF PCHFE THAT ASSUMES A STRICT- 
!        LY DECREASING X-ARRAY BY SIMPLY RUNNING THE IR-LOOP BACKWARDS  
!        (AND REVERSING THE ORDER OF APPROPRIATE TESTS).                
!                                                                       
!     4. THE PRESENT CODE HAS A MINOR BUG, WHICH I HAVE DECIDED IS NOT  
!        WORTH THE EFFORT THAT WOULD BE REQUIRED TO FIX IT.             
!        IF XE CONTAINS POINTS IN [X(N-1),X(N)], FOLLOWED BY POINTS .LT.
!        X(N-1), FOLLOWED BY POINTS .GT.X(N), THE EXTRAPOLATION POINTS  
!        WILL BE COUNTED (AT LEAST) TWICE IN THE TOTAL RETURNED IN IERR.
!                                                                       
!  DECLARE ARGUMENTS.                                                   
!                                                                       
      INTEGER  N, INCFD, NE, IERR 
      REAL  X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE) 
      LOGICAL  SKIP 
!                                                                       
!  DECLARE LOCAL VARIABLES.                                             
!                                                                       
      INTEGER  I, IERC, IR, J, JFIRST, NEXT(2), NJ 
!                                                                       
!  VALIDITY-CHECK ARGUMENTS.                                            
!                                                                       
!***FIRST EXECUTABLE STATEMENT  PCHFE                                   
      IF (SKIP)  GO TO 5 
!                                                                       
      IF ( N.LT.2 )  GO TO 5001 
      IF ( INCFD.LT.1 )  GO TO 5002 
      DO 1  I = 2, N 
         IF ( X(I).LE.X(I-1) )  GO TO 5003 
    1 END DO 
!                                                                       
!  FUNCTION DEFINITION IS OK, GO ON.                                    
!                                                                       
    5 CONTINUE 
      IF ( NE.LT.1 )  GO TO 5004 
      IERR = 0 
      SKIP = .TRUE. 
!                                                                       
!  LOOP OVER INTERVALS.        (   INTERVAL INDEX IS  IL = IR-1  . )    
!                              ( INTERVAL IS X(IL).LE.X.LT.X(IR) . )    
      JFIRST = 1 
      IR = 2 
   10 CONTINUE 
!                                                                       
!     SKIP OUT OF LOOP IF HAVE PROCESSED ALL EVALUATION POINTS.         
!                                                                       
         IF (JFIRST .GT. NE)  GO TO 5000 
!                                                                       
!     LOCATE ALL POINTS IN INTERVAL.                                    
!                                                                       
         DO 20  J = JFIRST, NE 
            IF (XE(J) .GE. X(IR))  GO TO 30 
   20    CONTINUE 
         J = NE + 1 
         GO TO 40 
!                                                                       
!     HAVE LOCATED FIRST POINT BEYOND INTERVAL.                         
!                                                                       
   30    CONTINUE 
         IF (IR .EQ. N)  J = NE + 1 
!                                                                       
   40    CONTINUE 
         NJ = J - JFIRST 
!                                                                       
!     SKIP EVALUATION IF NO POINTS IN INTERVAL.                         
!                                                                       
         IF (NJ .EQ. 0)  GO TO 50 
!                                                                       
!     EVALUATE CUBIC AT XE(I),  I = JFIRST (1) J-1 .                    
!                                                                       
!       ----------------------------------------------------------------
        CALL CHFEV (X(IR-1),X(IR), F(1,IR-1),F(1,IR), D(1,IR-1),D(1,IR),&
     &              NJ, XE(JFIRST), FE(JFIRST), NEXT, IERC)             
!       ----------------------------------------------------------------
         IF (IERC .LT. 0)  GO TO 5005 
!                                                                       
         IF (NEXT(2) .EQ. 0)  GO TO 42 
!        IF (NEXT(2) .GT. 0)  THEN                                      
!           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(2) TO THE   
!           RIGHT OF X(IR).                                             
!                                                                       
            IF (IR .LT. N)  GO TO 41 
!           IF (IR .EQ. N)  THEN                                        
!              THESE ARE ACTUALLY EXTRAPOLATION POINTS.                 
               IERR = IERR + NEXT(2) 
               GO TO 42 
   41       CONTINUE 
!           ELSE                                                        
!              WE SHOULD NEVER HAVE GOTTEN HERE.                        
               GO TO 5005 
!           ENDIF                                                       
!        ENDIF                                                          
   42    CONTINUE 
!                                                                       
         IF (NEXT(1) .EQ. 0)  GO TO 49 
!        IF (NEXT(1) .GT. 0)  THEN                                      
!           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(1) TO THE   
!           LEFT OF X(IR-1).                                            
!                                                                       
            IF (IR .GT. 2)  GO TO 43 
!           IF (IR .EQ. 2)  THEN                                        
!              THESE ARE ACTUALLY EXTRAPOLATION POINTS.                 
               IERR = IERR + NEXT(1) 
               GO TO 49 
   43       CONTINUE 
!           ELSE                                                        
!              XE IS NOT ORDERED RELATIVE TO X, SO MUST ADJUST          
!              EVALUATION INTERVAL.                                     
!                                                                       
!              FIRST, LOCATE FIRST POINT TO LEFT OF X(IR-1).            
               DO 44  I = JFIRST, J-1 
                  IF (XE(I) .LT. X(IR-1))  GO TO 45 
   44          CONTINUE 
!              NOTE-- CANNOT DROP THROUGH HERE UNLESS THERE IS AN ERROR 
!                     IN CHFEV.                                         
               GO TO 5005 
!                                                                       
   45          CONTINUE 
!              RESET J.  (THIS WILL BE THE NEW JFIRST.)                 
               J = I 
!                                                                       
!              NOW FIND OUT HOW FAR TO BACK UP IN THE X-ARRAY.          
               DO 46  I = 1, IR-1 
                  IF (XE(J) .LT. X(I)) GO TO 47 
   46          CONTINUE 
!              NB-- CAN NEVER DROP THROUGH HERE, SINCE XE(J).LT.X(IR-1).
!                                                                       
   47          CONTINUE 
!              AT THIS POINT, EITHER  XE(J) .LT. X(1)                   
!                 OR      X(I-1) .LE. XE(J) .LT. X(I) .                 
!              RESET IR, RECOGNIZING THAT IT WILL BE INCREMENTED BEFORE 
!              CYCLING.                                                 
               IR = MAX0(1, I-1) 
!           ENDIF                                                       
!        ENDIF                                                          
   49    CONTINUE 
!                                                                       
         JFIRST = J 
!                                                                       
!     END OF IR-LOOP.                                                   
!                                                                       
   50 CONTINUE 
      IR = IR + 1 
      IF (IR .LE. N)  GO TO 10 
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
      CALL XERROR ('PCHFE -- NUMBER OF DATA POINTS LESS THAN TWO'       &
     &           , 44, IERR, 1)                                         
      RETURN 
!                                                                       
 5002 CONTINUE 
!     INCFD.LT.1 RETURN.                                                
      IERR = -2 
      CALL XERROR ('PCHFE -- INCREMENT LESS THAN ONE'                   &
     &           , 32, IERR, 1)                                         
      RETURN 
!                                                                       
 5003 CONTINUE 
!     X-ARRAY NOT STRICTLY INCREASING.                                  
      IERR = -3 
      CALL XERROR ('PCHFE -- X-ARRAY NOT STRICTLY INCREASING'           &
     &           , 40, IERR, 1)                                         
      RETURN 
!                                                                       
 5004 CONTINUE 
!     NE.LT.1 RETURN.                                                   
      IERR = -4 
      CALL XERROR ('PCHFE -- NUMBER OF EVALUATION POINTS LESS THAN ONE' &
     &           , 50, IERR, 1)                                         
      RETURN 
!                                                                       
 5005 CONTINUE 
!     ERROR RETURN FROM CHFEV.                                          
!   *** THIS CASE SHOULD NEVER OCCUR ***                                
      IERR = -5 
      CALL XERROR ('PCHFE -- ERROR RETURN FROM CHFEV -- FATAL'          &
     &           , 41, IERR, 2)                                         
      RETURN 
!------------- LAST LINE OF PCHFE FOLLOWS ------------------------------
      END                                           
      SUBROUTINE CHFEV(X1,X2,F1,F2,D1,D2,NE,XE,FE,NEXT,IERR) 
!***BEGIN PROLOGUE  CHFEV                                               
!***DATE WRITTEN   811019   (YYMMDD)                                    
!***REVISION DATE  870707   (YYMMDD)                                    
!***CATEGORY NO.  E3,H1                                                 
!***KEYWORDS  LIBRARY=SLATEC(PCHIP),                                    
!             TYPE=SINGLE PRECISION(CHFEV-S DCHFEV-D),                  
!             CUBIC HERMITE EVALUATION,CUBIC POLYNOMIAL EVALUATION      
!***AUTHOR  FRITSCH, F. N., (LLNL)                                      
!             MATHEMATICS AND STATISTICS DIVISION                       
!             LAWRENCE LIVERMORE NATIONAL LABORATORY                    
!             P.O. BOX 808  (L-316)                                     
!             LIVERMORE, CA  94550                                      
!             FTS 532-4275, (415) 422-4275                              
!***PURPOSE  EVALUATE A CUBIC POLYNOMIAL GIVEN IN HERMITE FORM AT AN    
!            ARRAY OF POINTS.  WHILE DESIGNED FOR USE BY PCHFE, IT MAY  
!            BE USEFUL DIRECTLY AS AN EVALUATOR FOR A PIECEWISE CUBIC   
!            HERMITE FUNCTION IN APPLICATIONS, SUCH AS GRAPHING, WHERE  
!            THE INTERVAL IS KNOWN IN ADVANCE.                          
!***DESCRIPTION                                                         
!                                                                       
!          CHFEV:  CUBIC HERMITE FUNCTION EVALUATOR                     
!                                                                       
!     EVALUATES THE CUBIC POLYNOMIAL DETERMINED BY FUNCTION VALUES      
!     F1,F2 AND DERIVATIVES D1,D2 ON INTERVAL (X1,X2) AT THE POINTS     
!     XE(J), J=1(1)NE.                                                  
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  CALLING SEQUENCE:                                                    
!                                                                       
!        INTEGER  NE, NEXT(2), IERR                                     
!        REAL  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE)                   
!                                                                       
!        CALL  CHFEV (X1,X2, F1,F2, D1,D2, NE, XE, FE, NEXT, IERR)      
!                                                                       
!   PARAMETERS:                                                         
!                                                                       
!     X1,X2 -- (INPUT) ENDPOINTS OF INTERVAL OF DEFINITION OF CUBIC.    
!           (ERROR RETURN IF  X1.EQ.X2 .)                               
!                                                                       
!     F1,F2 -- (INPUT) VALUES OF FUNCTION AT X1 AND X2, RESPECTIVELY.   
!                                                                       
!     D1,D2 -- (INPUT) VALUES OF DERIVATIVE AT X1 AND X2, RESPECTIVELY. 
!                                                                       
!     NE -- (INPUT) NUMBER OF EVALUATION POINTS.  (ERROR RETURN IF      
!           NE.LT.1 .)                                                  
!                                                                       
!     XE -- (INPUT) REAL ARRAY OF POINTS AT WHICH THE FUNCTION IS TO BE 
!           EVALUATED.  IF ANY OF THE XE ARE OUTSIDE THE INTERVAL       
!           [X1,X2], A WARNING ERROR IS RETURNED IN NEXT.               
!                                                                       
!     FE -- (OUTPUT) REAL ARRAY OF VALUES OF THE CUBIC FUNCTION DEFINED 
!           BY  X1,X2, F1,F2, D1,D2  AT THE POINTS  XE.                 
!                                                                       
!     NEXT -- (OUTPUT) INTEGER ARRAY INDICATING NUMBER OF EXTRAPOLATION 
!           POINTS:                                                     
!            NEXT(1) = NUMBER OF EVALUATION POINTS TO LEFT OF INTERVAL. 
!            NEXT(2) = NUMBER OF EVALUATION POINTS TO RIGHT OF INTERVAL.
!                                                                       
!     IERR -- (OUTPUT) ERROR FLAG.                                      
!           NORMAL RETURN:                                              
!              IERR = 0  (NO ERRORS).                                   
!           "RECOVERABLE" ERRORS:                                       
!              IERR = -1  IF NE.LT.1 .                                  
!              IERR = -2  IF X1.EQ.X2 .                                 
!                (THE FE-ARRAY HAS NOT BEEN CHANGED IN EITHER CASE.)    
!                                                                       
!***REFERENCES  (NONE)                                                  
!***ROUTINES CALLED  XERROR                                             
!***END PROLOGUE  CHFEV                                                 
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  CHANGE RECORD:                                                       
!     82-08-03   MINOR COSMETIC CHANGES FOR RELEASE 1.                  
!                                                                       
! ----------------------------------------------------------------------
!                                                                       
!  PROGRAMMING NOTES:                                                   
!                                                                       
!     TO PRODUCE A DOUBLE PRECISION VERSION, SIMPLY:                    
!        A. CHANGE CHFEV TO DCHFEV WHEREVER IT OCCURS,                  
!        B. CHANGE THE REAL DECLARATION TO DOUBLE PRECISION,            
!        C. CHANGE THE CONSTANT ZERO TO DOUBLE PRECISION, AND           
!        D. CHANGE THE NAMES OF THE FORTRAN FUNCTIONS:  AMAX1, AMIN1.   
!                                                                       
!  DECLARE ARGUMENTS.                                                   
!                                                                       
      INTEGER  NE, NEXT(2), IERR 
      REAL  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE) 
!                                                                       
!  DECLARE LOCAL VARIABLES.                                             
!                                                                       
      INTEGER  I 
      REAL  C2, C3, DEL1, DEL2, DELTA, H, X, XMI, XMA, ZERO 
      DATA  ZERO /0./ 
!                                                                       
!  VALIDITY-CHECK ARGUMENTS.                                            
!                                                                       
!***FIRST EXECUTABLE STATEMENT  CHFEV                                   
      IF (NE .LT. 1)  GO TO 5001 
      H = X2 - X1 
      IF (H .EQ. ZERO)  GO TO 5002 
!                                                                       
!  INITIALIZE.                                                          
!                                                                       
      IERR = 0 
      NEXT(1) = 0 
      NEXT(2) = 0 
      XMI = AMIN1(ZERO, H) 
      XMA = AMAX1(ZERO, H) 
!                                                                       
!  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).                      
!                                                                       
      DELTA = (F2 - F1)/H 
      DEL1 = (D1 - DELTA)/H 
      DEL2 = (D2 - DELTA)/H 
!                                           (DELTA IS NO LONGER NEEDED.)
      C2 = -(DEL1+DEL1 + DEL2) 
      C3 = (DEL1 + DEL2)/H 
!                               (H, DEL1 AND DEL2 ARE NO LONGER NEEDED.)
!                                                                       
!  EVALUATION LOOP.                                                     
!                                                                       
      DO 500  I = 1, NE 
         X = XE(I) - X1 
         FE(I) = F1 + X*(D1 + X*(C2 + X*C3)) 
!          COUNT EXTRAPOLATION POINTS.                                  
         IF ( X.LT.XMI )  NEXT(1) = NEXT(1) + 1 
         IF ( X.GT.XMA )  NEXT(2) = NEXT(2) + 1 
!        (NOTE REDUNDANCY--IF EITHER CONDITION IS TRUE, OTHER IS FALSE.)
  500 END DO 
!                                                                       
!  NORMAL RETURN.                                                       
!                                                                       
      RETURN 
!                                                                       
!  ERROR RETURNS.                                                       
!                                                                       
 5001 CONTINUE 
!     NE.LT.1 RETURN.                                                   
      IERR = -1 
      CALL XERROR ('CHFEV -- NUMBER OF EVALUATION POINTS LESS THAN ONE' &
     &           , 50, IERR, 1)                                         
      RETURN 
!                                                                       
 5002 CONTINUE 
!     X1.EQ.X2 RETURN.                                                  
      IERR = -2 
      CALL XERROR ('CHFEV -- INTERVAL ENDPOINTS EQUAL'                  &
     &           , 33, IERR, 1)                                         
      RETURN 
!------------- LAST LINE OF CHFEV FOLLOWS ------------------------------
      END                                           
