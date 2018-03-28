!***BEGIN PROLOGUE  EFC                                                 
!***DATE WRITTEN   800801   (YYMMDD)                                    
!***REVISION DATE  820801   (YYMMDD)                                    
!***CATEGORY NO.  K1A1A1,K1A2A,L8A3                                     
!***KEYWORDS  B-SPLINES,CONSTRAINED LEAST SQUARES,CURVE FITTING,        
!             LEAST SQUARES                                             
!***AUTHOR  HANSON, R. J., (SNLA)                                       
!***PURPOSE  FITS A PIECE-WISE POLYNOMIAL CURVE TO DISCRETE             
!            DATA.  THE PIECE-WISE POLYNOMIALS ARE REPRESENTED          
!            AS B-SPLINES.  THE FITTING IS DONE IN A WEIGHTED           
!            LEAST SQUARES SENSE.                                       
!***DESCRIPTION                                                         
!                                                                       
!     REVISED 800905-1300                                               
!     REVISED YYMMDD-HHMM                                               
!     DIMENSION XDATA(NDATA),YDATA(NDATA),SDDATA(NDATA),  BKPT(NBKPT),  
!    1  COEFF(NBKPT-NORD),W(*)                                          
!                                                                       
!      THIS SUBPROGRAM FITS A PIECE-WISE POLYNOMIAL CURVE               
!      TO DISCRETE DATA.  THE PIECE-WISE POLYNOMIALS ARE                
!      REPRESENTED AS B-SPLINES.                                        
!      THE FITTING IS DONE IN A WEIGHTED LEAST SQUARES SENSE.           
!                                                                       
!      THE DATA CAN BE PROCESSED IN GROUPS OF MODEST SIZE.              
!      THE SIZE OF THE GROUP IS CHOSEN BY THE USER.  THIS FEATURE       
!      MAY BE NECESSARY FOR PURPOSES OF USING CONSTRAINED CURVE FITTING 
!      WITH SUBPROGRAM FC( ) ON A VERY LARGE DATA SET.                  
!      FOR A DESCRIPTION OF THE B-SPLINES AND USAGE INSTRUCTIONS TO     
!      EVALUATE THEM, SEE                                               
!                                                                       
!      C. W. DE BOOR, PACKAGE FOR CALCULATING WITH B-SPLINES.           
!                     SIAM J. NUMER. ANAL., P. 441, (JUNE, 1977).       
!                                                                       
!      FOR FURTHER DISCUSSION OF (CONSTRAINED) CURVE FITTING USING      
!      B-SPLINES, SEE                                                   
!                                                                       
!      R. J. HANSON, CONSTRAINED LEAST SQUARES CURVE FITTING            
!                   TO DISCRETE DATA USING B-SPLINES, A USER'S          
!                   GUIDE. SANDIA LABS. TECH. REPT. SAND-78-1291,       
!                   DECEMBER, (1978).                                   
!                                                                       
!  INPUT..                                                              
!      NDATA,XDATA(*),                                                  
!      YDATA(*),                                                        
!      SDDATA(*)                                                        
!                         THE NDATA DISCRETE (X,Y) PAIRS AND            
!                         THE Y VALUE STANDARD DEVIATION OR             
!                         UNCERTAINTY, SD, ARE IN THE RESPECTIVE        
!                         ARRAYS XDATA(*), YDATA(*), AND SDDATA(*).     
!                         NO SORTING OF XDATA(*) IS REQUIRED.  ANY      
!                         NON-NEGATIVE VALUE OF NDATA IS ALLOWED.  A    
!                         NEGATIVE VALUE OF NDATA IS AN ERROR.          
!                         A ZERO VALUE FOR ANY ENTRY OF SDDATA(*)       
!                         WILL WEIGHT THAT DATA POINT AS 1.             
!                         OTHERWISE THE WEIGHT OF THAT DATA POINT IS    
!                         THE RECIPROCAL OF THIS ENTRY.                 
!                                                                       
!      NORD,NBKPT,                                                      
!      BKPT(*)                                                          
!                         THE NBKPT KNOTS OF THE B-SPLINE OF ORDER      
!                         NORD ARE IN THE ARRAY BKPT(*).  NORMALLY      
!                         THE PROBLEM DATA INTERVAL WILL BE INCLUDED    
!                         BETWEEN THE LIMITS BKPT(NORD) AND             
!                         BKPT(NBKPT-NORD+1).  THE ADDITIONAL END KNOTS 
!                         BKPT(I),I=1,...,NORD-1 AND I=NBKPT-NORD+2,...,
!                         NBKPT, ARE REQUIRED TO COMPUTE THE FUNCTIONS  
!                         USED TO FIT THE DATA.                         
!                         NO SORTING OF BKPT(*) IS REQUIRED.  INTERNAL  
!                         TO EFC( ) THE EXTREME END KNOTS MAY BE SLIGHT-
!                         LY REDUCED AND INCREASED RESPECTIVELY TO      
!                         ACCOMODATE ANY DATA VALUES THAT ARE EXTERIOR  
!                         TO THE GIVEN KNOT VALUES.  THE CONTENTS OF    
!                         BKPT(*) ARE NOT CHANGED.                      
!                                                                       
!                         THE VALUE OF NORD MUST LIE BETWEEN 1 AND 20.  
!                         THE VALUE OF NBKPT MUST SATISFY NBKPT .GE.    
!                         2*NORD.  OTHER VALUES ARE CONSIDERED ERRORS.  
!                                                                       
!                         (THE ORDER OF THE SPLINE IS ONE MORE          
!                         THAN THE DEGREE OF THE PIECE-WISE POLYNOMIAL  
!                         DEFINED ON EACH INTERVAL.  THIS IS CONSISTENT 
!                         WITH THE B-SPLINE PACKAGE CONVENTION.  FOR    
!                         EXAMPLE, NORD=4 WHEN WE ARE USING PIECE-WISE  
!                         CUBICS.)                                      
!                                                                       
!        MDEIN                                                          
!                         AN INTEGER FLAG, WITH ONE OF TWO POSSIBLE     
!                         VALUES (1 OR 2), THAT DIRECTS THE SUBPROGRAM  
!                         ACTION WITH REGARD TO NEW DATA POINTS PROVIDED
!                         BY THE USER.                                  
!                                                                       
!                         =1  THE FIRST TIME THAT EFC( ) HAS BEEN       
!                         ENTERED.  THERE ARE NDATA POINTS TO PROCESS.  
!                                                                       
!                         =2  THIS IS ANOTHER ENTRY TO EFC( ).  THE SUB-
!                         PROGRAM EFC( ) HAS BEEN ENTERED WITH MDEIN=1  
!                         EXACTLY ONCE BEFORE FOR THIS PROBLEM.  THERE  
!                         ARE NDATA NEW ADDITIONAL POINTS TO MERGE AND  
!                         PROCESS WITH ANY PREVIOUS POINTS.             
!                         (WHEN USING EFC( ) WITH MDEIN=2 IT IS IMPORT- 
!                         ANT THAT THE SET OF KNOTS REMAIN FIXED AT THE 
!                         SAME VALUES FOR ALL ENTRIES TO EFC( ).)       
!       LW                                                              
!                         THE AMOUNT OF WORKING STORAGE ACTUALLY        
!                         ALLOCATED FOR THE WORKING ARRAY W(*).         
!                         THIS QUANTITY IS COMPARED WITH THE            
!                         ACTUAL AMOUNT OF STORAGE NEEDED IN EFC( ).    
!                         INSUFFICIENT STORAGE ALLOCATED FOR W(*) IS    
!                         AN ERROR.  THIS FEATURE WAS INCLUDED IN EFC( )
!                         BECAUSE MISREADING THE STORAGE FORMULA        
!                         FOR W(*) MIGHT VERY WELL LEAD TO SUBTLE       
!                         AND HARD-TO-FIND PROGRAMMING BUGS.            
!                                                                       
!                         THE LENGTH OF THE ARRAY W(*) MUST SATISFY     
!                                                                       
!                         LW .GE. (NBKPT-NORD+3)*(NORD+1)+              
!                                 (NBKPT+1)*(NORD+1)+                   
!                               2*MAX0(NDATA,NBKPT)+NBKPT+NORD**2       
!  OUTPUT..                                                             
!      MDEOUT                                                           
!                         AN OUTPUT FLAG THAT INDICATES THE STATUS      
!                         OF THE CONSTRAINED CURVE FIT.                 
!                                                                       
!                         =-1  A USAGE ERROR OF EFC( ) OCCURRED.  THE   
!                         OFFENDING CONDITION IS NOTED WITH THE SLATEC  
!                         LIBRARY ERROR PROCESSOR, XERROR( ).           
!                         IN CASE THE WORKING ARRAY W(*) IS             
!                         NOT LONG ENOUGH, THE MINIMAL ACCEPTABLE LENGTH
!                         IS PRINTED USING THE ERROR PROCESSING SUBPRO- 
!                         GRAM XERRWV( ).                               
!                                                                       
!                         =1  THE B-SPLINE COEFFICIENTS FOR THE FITTED  
!                         CURVE HAVE BEEN RETURNED IN ARRAY COEFF(*).   
!                                                                       
!                         =2  NOT ENOUGH DATA HAS BEEN PROCESSED TO     
!                         DETERMINE THE B-SPLINE COEFFICIENTS.          
!                         THE USER HAS ONE OF TWO OPTIONS.  CONTINUE    
!                         TO PROCESS MORE DATA UNTIL A UNIQUE SET       
!                         OF COEFFICIENTS IS OBTAINED, OR USE THE       
!                         SUBPROGRAM FC( ) TO OBTAIN A SPECIFIC         
!                         SET OF COEFFICIENTS.  THE USER SHOULD READ    
!                         THE USAGE INSTRUCTIONS FOR FC( ) FOR FURTHER  
!                         DETAILS IF THIS SECOND OPTION IS CHOSEN.      
!      COEFF(*)                                                         
!                         IF THE OUTPUT VALUE OF MDEOUT=1, THIS ARRAY   
!                         CONTAINS THE UNKNOWNS OBTAINED FROM THE LEAST 
!                         SQUARES FITTING PROCESS.  THESE N=NBKPT-NORD  
!                         PARAMETERS ARE THE B-SPLINE COEFFICIENTS.     
!                         FOR MDEOUT=2, NOT ENOUGH DATA WAS PROCESSED TO
!                         UNIQUELY DETERMINE THE B-SPLINE COEFFICIENTS. 
!                         IN THIS CASE, AND ALSO WHEN MDEOUT=-1, ALL    
!                         VALUES OF COEFF(*) ARE SET TO ZERO.           
!                                                                       
!                         IF THE USER IS NOT SATISFIED WITH THE FITTED  
!                         CURVE RETURNED BY EFC( ), THE CONSTRAINED     
!                         LEAST SQUARES CURVE FITTING SUBPROGRAM FC( )  
!                         MAY BE REQUIRED.  THE WORK DONE WITHIN EFC( ) 
!                         TO ACCUMULATE THE DATA CAN BE UTILIZED BY     
!                         THE USER, IF SO DESIRED.  THIS INVOLVES       
!                         SAVING THE FIRST (NBKPT-NORD+3)*(NORD+1)      
!                         ENTRIES OF W(*) AND PROVIDING THIS DATA       
!                         TO FC( ) WITH THE "OLD PROBLEM" DESIGNATION.  
!                         THE USER SHOULD READ THE USAGE INSTRUCTIONS   
!                         FOR SUBPROGRAM FC( ) FOR FURTHER DETAILS.     
!                                                                       
!  WORKING ARRAY..                                                      
!      W(*)                                                             
!                         THIS ARRAY IS TYPED REAL.  ITS                
!                         LENGTH IS SPECIFIED AS AN INPUT               
!                         PARAMETER IN LW AS NOTED ABOVE.               
!                         THE CONTENTS OF W(*) MUST NOT BE MODIFIED     
!                         BY THE USER BETWEEN CALLS TO EFC( ) WITH      
!                         VALUES OF MDEIN=1,2,2,... .  THE FIRST        
!                         (NBKPT-NORD+3)*(NORD+1) ENTRIES OF W(*) ARE   
!                         ACCEPTABLE AS DIRECT INPUT TO FC( ) FOR AN    
!                         "OLD PROBLEM" ONLY WHEN MDEOUT=1 OR 2.        
!                                                                       
!  EVALUATING THE                                                       
!  FITTED CURVE..                                                       
!                         TO EVALUATE DERIVATIVE NUMBER IDER AT XVAL    
!                         USE THE FUNCTION SUBPROGRAM BVALU( ).         
!                                                                       
!                           F = BVALU(BKPT,COEFF,NBKPT-NORD,NORD,IDER,  
!                                      XVAL,INBV,WORKB)                 
!                                                                       
!                         THE OUTPUT OF THIS SUBPROGRAM WILL NOT BE     
!                         DEFINED UNLESS AN OUTPUT VALUE OF MDEOUT=1    
!                         WAS OBTAINED FROM EFC( ), XVAL IS IN THE DATA 
!                         INTERVAL, AND IDER IS NONNEGATIVE AND .LT.    
!                         NORD.                                         
!                         THE FIRST TIME BVALU( ) IS CALLED, INBV=1     
!                         MUST BE SPECIFIED.  THIS VALUE OF INBV IS THE 
!                         OVERWRITTEN BY BVALU( ).  THE ARRAY WORKB(*)  
!                         MUST BE OF LENGTH AT LEAST 3*NORD, AND MUST   
!                         NOT BE THE SAME AS THE W(*) ARRAY USED        
!                         IN THE CALL TO EFC( ).                        
!                                                                       
!                         BVALU( ) EXPECTS THE BREAKPOINT ARRAY BKPT(*  
!                         TO BE SORTED.                                 
!***REFERENCES  HANSON R.J., *CONSTRAINED LEAST SQUARES CURVE FITTING   
!                 TO DISCRETE DATA USING B-SPLINES, A USERS GUIDE*,     
      SUBROUTINE EFC(NDATA,XDATA,YDATA,SDDATA,NORD,NBKPT,BKPT,MDEIN,    &
     &   MDEOUT,COEFF,LW,W)                                             
!                 SAND78-1291, DECEMBER,1978.                           
!***ROUTINES CALLED  EFCMN                                              
!***END PROLOGUE  EFC                                                   
!                                                                       
!+-+-+-+-                                                               
!     THE USER-PROVIDED USAGE INSTRUCTIONS CAN END HERE.                
!      SUBROUTINE           FUNCTION/REMARKS                            
!                                                                       
!      BSPLVN( )          COMPUTE FUNCTION VALUES OF B-SPLINES.  FROM   
!                         THE B-SPLINE PACKAGE OF DE BOOR NOTED ABOVE.  
!                                                                       
!      BNDACC( ),         BANDED LEAST SQUARES MATRIX PROCESSORS.       
!      BNDSOL( )          FROM LAWSON-HANSON, SOLVING LEAST             
!                         SQUARES PROBLEMS.                             
!                                                                       
!      SSORT( )           DATA SORTING SUBROUTINE, FROM THE             
!                         SANDIA MATH. LIBRARY, SAND77-1441.            
!      XERROR( ),         ERROR HANDLING ROUTINES                       
!      XERRWV( )          FOR THE SLATEC MATH. LIBRARY.                 
!                         SEE SAND78-1189, BY R. E. JONES.              
!                                                                       
!      SCOPY( ),SSCAL( )  SUBROUTINES FROM THE BLAS PACKAGE.            
!                         LOCATED ON SANDIA. MATH. LIBRARY. SEE         
!                         SAND77-0898 FOR DESCRIPTION.                  
!                                                                       
!                         WRITTEN BY R. HANSON, SANDIA NATL. LABS.,     
!                         ALB., N. M., AUGUST-SEPTEMBER, 1980.          
      DIMENSION XDATA(1),YDATA(1),SDDATA(1),BKPT(1) 
      DIMENSION COEFF(1),W(1) 
!***FIRST EXECUTABLE STATEMENT  EFC                                     
      MDG=NBKPT+1 
      MDW=NBKPT-NORD+3 
!     LWW=1               USAGE IN EFCMN( ) OF W(*)..                   
!     LWW,...,LG-1        W(*,*)                                        
!                                                                       
!     LG,...,LXTEMP-1     G(*,*)                                        
!                                                                       
!     LXTEMP,...,LPTEMP-1 XTEMP(*)                                      
!                                                                       
!     LPTEMP,...,LBKPT-1  PTEMP(*)                                      
!                                                                       
!     LBKPT,...,LBF       BKPT(*) (LOCAL TO EFCMN( ))                   
!                                                                       
!     LBF,...,LBF+NORD**2 BF(*,*)                                       
!                                                                       
      LWW=1 
      LG=LWW+MDW*(NORD+1) 
      LXTEMP=LG+MDG*(NORD+1) 
      LPTEMP=LXTEMP+MAX0(NDATA,NBKPT) 
      LBKPT=LPTEMP+MAX0(NDATA,NBKPT) 
      LBF=LBKPT+NBKPT 
      CALL EFCMN(NDATA,XDATA,YDATA,SDDATA,                              &
     &         NORD,NBKPT,BKPT,                                         &
     &         MDEIN,MDEOUT,                                            &
     &         COEFF,                                                   &
     &         W(LBF),W(LXTEMP),W(LPTEMP),W(LBKPT),                     &
     &         W(LG),MDG,W(LWW),MDW,LW)                                 
      RETURN 
      END                                           
      SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY) 
!***BEGIN PROLOGUE  SAXPY                                               
!***DATE WRITTEN   791001   (YYMMDD)                                    
!***REVISION DATE  820801   (YYMMDD)                                    
!***CATEGORY NO.  D1A7                                                  
!***KEYWORDS  BLAS,LINEAR ALGEBRA,TRIAD,VECTOR                          
!***AUTHOR  LAWSON, C. L., (JPL)                                        
!           HANSON, R. J., (SNLA)                                       
!           KINCAID, D. R., (U. OF TEXAS)                               
!           KROGH, F. T., (JPL)                                         
!***PURPOSE  S.P. COMPUTATION Y = A*X + Y                               
!***DESCRIPTION                                                         
!                                                                       
!                B L A S  SUBPROGRAM                                    
!    DESCRIPTION OF PARAMETERS                                          
!                                                                       
!     --INPUT--                                                         
!        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)                       
!       SA  SINGLE PRECISION SCALAR MULTIPLIER                          
!       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS                     
!     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX                      
!       SY  SINGLE PRECISION VECTOR WITH N ELEMENTS                     
!     INCY  STORAGE SPACING BETWEEN ELEMENTS OF SY                      
!                                                                       
!     --OUTPUT--                                                        
!       SY  SINGLE PRECISION RESULT (UNCHANGED IF N .LE. 0)             
!                                                                       
!     OVERWRITE SINGLE PRECISION SY WITH SINGLE PRECISION SA*SX +SY.    
!     FOR I = 0 TO N-1, REPLACE  SY(LY+I*INCY) WITH SA*SX(LX+I*INCX) +  
!       SY(LY+I*INCY), WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N 
!       AND LY IS DEFINED IN A SIMILAR WAY USING INCY.                  
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,     
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*, 
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL       
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323 
!***ROUTINES CALLED  (NONE)                                             
!***END PROLOGUE  SAXPY                                                 
!                                                                       
      REAL SX(1),SY(1),SA 
!***FIRST EXECUTABLE STATEMENT  SAXPY                                   
      IF(N.LE.0.OR.SA.EQ.0.E0) RETURN 
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60 
    5 CONTINUE 
!                                                                       
!        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.                   
!                                                                       
      IX = 1 
      IY = 1 
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1 
      DO 10 I = 1,N 
        SY(IY) = SY(IY) + SA*SX(IX) 
        IX = IX + INCX 
        IY = IY + INCY 
   10 END DO 
      RETURN 
!                                                                       
!        CODE FOR BOTH INCREMENTS EQUAL TO 1                            
!                                                                       
!                                                                       
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.   
!                                                                       
   20 M = MOD(N,4) 
      IF( M .EQ. 0 ) GO TO 40 
      DO 30 I = 1,M 
        SY(I) = SY(I) + SA*SX(I) 
   30 END DO 
      IF( N .LT. 4 ) RETURN 
   40 MP1 = M + 1 
      DO 50 I = MP1,N,4 
        SY(I) = SY(I) + SA*SX(I) 
        SY(I + 1) = SY(I + 1) + SA*SX(I + 1) 
        SY(I + 2) = SY(I + 2) + SA*SX(I + 2) 
        SY(I + 3) = SY(I + 3) + SA*SX(I + 3) 
   50 END DO 
      RETURN 
!                                                                       
!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.                  
!                                                                       
   60 CONTINUE 
      NS = N*INCX 
          DO 70 I=1,NS,INCX 
          SY(I) = SA*SX(I) + SY(I) 
   70     CONTINUE 
      RETURN 
      END                                           
      SUBROUTINE SCOPY(N,SX,INCX,SY,INCY) 
!***BEGIN PROLOGUE  SCOPY                                               
!***DATE WRITTEN   791001   (YYMMDD)                                    
!***REVISION DATE  820801   (YYMMDD)                                    
!***CATEGORY NO.  D1A5                                                  
!***KEYWORDS  BLAS,COPY,LINEAR ALGEBRA,VECTOR                           
!***AUTHOR  LAWSON, C. L., (JPL)                                        
!           HANSON, R. J., (SNLA)                                       
!           KINCAID, D. R., (U. OF TEXAS)                               
!           KROGH, F. T., (JPL)                                         
!***PURPOSE  COPY S.P. VECTOR Y = X                                     
!***DESCRIPTION                                                         
!                                                                       
!                B L A S  SUBPROGRAM                                    
!    DESCRIPTION OF PARAMETERS                                          
!                                                                       
!     --INPUT--                                                         
!        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)                       
!       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS                     
!     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX                      
!       SY  SINGLE PRECISION VECTOR WITH N ELEMENTS                     
!     INCY  STORAGE SPACING BETWEEN ELEMENTS OF SY                      
!                                                                       
!     --OUTPUT--                                                        
!       SY  COPY OF VECTOR SX (UNCHANGED IF N .LE. 0)                   
!                                                                       
!     COPY SINGLE PRECISION SX TO SINGLE PRECISION SY.                  
!     FOR I = 0 TO N-1, COPY  SX(LX+I*INCX) TO SY(LY+I*INCY),           
!     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS       
!     DEFINED IN A SIMILAR WAY USING INCY.                              
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,     
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*, 
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL       
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323 
!***ROUTINES CALLED  (NONE)                                             
!***END PROLOGUE  SCOPY                                                 
!                                                                       
      REAL SX(1),SY(1) 
!***FIRST EXECUTABLE STATEMENT  SCOPY                                   
      IF(N.LE.0)RETURN 
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60 
    5 CONTINUE 
!                                                                       
!        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.                    
!                                                                       
      IX = 1 
      IY = 1 
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1 
      DO 10 I = 1,N 
        SY(IY) = SX(IX) 
        IX = IX + INCX 
        IY = IY + INCY 
   10 END DO 
      RETURN 
!                                                                       
!        CODE FOR BOTH INCREMENTS EQUAL TO 1                            
!                                                                       
!                                                                       
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7.   
!                                                                       
   20 M = MOD(N,7) 
      IF( M .EQ. 0 ) GO TO 40 
      DO 30 I = 1,M 
        SY(I) = SX(I) 
   30 END DO 
      IF( N .LT. 7 ) RETURN 
   40 MP1 = M + 1 
      DO 50 I = MP1,N,7 
        SY(I) = SX(I) 
        SY(I + 1) = SX(I + 1) 
        SY(I + 2) = SX(I + 2) 
        SY(I + 3) = SX(I + 3) 
        SY(I + 4) = SX(I + 4) 
        SY(I + 5) = SX(I + 5) 
        SY(I + 6) = SX(I + 6) 
   50 END DO 
      RETURN 
!                                                                       
!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.                  
!                                                                       
   60 CONTINUE 
      NS = N*INCX 
          DO 70 I=1,NS,INCX 
          SY(I) = SX(I) 
   70     CONTINUE 
      RETURN 
      END                                           
      REAL FUNCTION SDOT(N,SX,INCX,SY,INCY) 
!***BEGIN PROLOGUE  SDOT                                                
!***DATE WRITTEN   791001   (YYMMDD)                                    
!***REVISION DATE  820801   (YYMMDD)                                    
!***CATEGORY NO.  D1A4                                                  
!***KEYWORDS  BLAS,INNER PRODUCT,LINEAR ALGEBRA,VECTOR                  
!***AUTHOR  LAWSON, C. L., (JPL)                                        
!           HANSON, R. J., (SNLA)                                       
!           KINCAID, D. R., (U. OF TEXAS)                               
!           KROGH, F. T., (JPL)                                         
!***PURPOSE  S.P. INNER PRODUCT OF S.P. VECTORS                         
!***DESCRIPTION                                                         
!                                                                       
!                B L A S  SUBPROGRAM                                    
!    DESCRIPTION OF PARAMETERS                                          
!                                                                       
!     --INPUT--                                                         
!        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)                       
!       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS                     
!     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX                      
!       SY  SINGLE PRECISION VECTOR WITH N ELEMENTS                     
!     INCY  STORAGE SPACING BETWEEN ELEMENTS OF SY                      
!                                                                       
!     --OUTPUT--                                                        
!     SDOT  SINGLE PRECISION DOT PRODUCT (ZERO IF N .LE. 0)             
!                                                                       
!     RETURNS THE DOT PRODUCT OF SINGLE PRECISION SX AND SY.            
!     SDOT = SUM FOR I = 0 TO N-1 OF  SX(LX+I*INCX) * SY(LY+I*INCY),    
!     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS       
!     DEFINED IN A SIMILAR WAY USING INCY.                              
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,     
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*, 
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL       
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323 
!***ROUTINES CALLED  (NONE)                                             
!***END PROLOGUE  SDOT                                                  
!                                                                       
      REAL SX(1),SY(1) 
!***FIRST EXECUTABLE STATEMENT  SDOT                                    
      SDOT = 0.0E0 
      IF(N.LE.0)RETURN 
      IF(INCX.EQ.INCY) IF(INCX-1)5,20,60 
    5 CONTINUE 
!                                                                       
!        CODE FOR UNEQUAL INCREMENTS OR NONPOSITIVE INCREMENTS.         
!                                                                       
      IX = 1 
      IY = 1 
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1 
      DO 10 I = 1,N 
        SDOT = SDOT + SX(IX)*SY(IY) 
        IX = IX + INCX 
        IY = IY + INCY 
   10 END DO 
      RETURN 
!                                                                       
!        CODE FOR BOTH INCREMENTS EQUAL TO 1                            
!                                                                       
!                                                                       
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.   
!                                                                       
   20 M = MOD(N,5) 
      IF( M .EQ. 0 ) GO TO 40 
      DO 30 I = 1,M 
        SDOT = SDOT + SX(I)*SY(I) 
   30 END DO 
      IF( N .LT. 5 ) RETURN 
   40 MP1 = M + 1 
      DO 50 I = MP1,N,5 
        SDOT = SDOT + SX(I)*SY(I) + SX(I + 1)*SY(I + 1) +               &
     &   SX(I + 2)*SY(I + 2) + SX(I + 3)*SY(I + 3) + SX(I + 4)*SY(I + 4)
   50 END DO 
      RETURN 
!                                                                       
!        CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.                      
!                                                                       
   60 CONTINUE 
      NS=N*INCX 
      DO 70 I=1,NS,INCX 
        SDOT = SDOT + SX(I)*SY(I) 
   70   CONTINUE 
      RETURN 
      END                                           
      SUBROUTINE SSCAL(N,SA,SX,INCX) 
!***BEGIN PROLOGUE  SSCAL                                               
!***DATE WRITTEN   791001   (YYMMDD)                                    
!***REVISION DATE  820801   (YYMMDD)                                    
!***CATEGORY NO.  D1A6                                                  
!***KEYWORDS  BLAS,LINEAR ALGEBRA,SCALE,VECTOR                          
!***AUTHOR  LAWSON, C. L., (JPL)                                        
!           HANSON, R. J., (SNLA)                                       
!           KINCAID, D. R., (U. OF TEXAS)                               
!           KROGH, F. T., (JPL)                                         
!***PURPOSE  S.P. VECTOR SCALE X = A*X                                  
!***DESCRIPTION                                                         
!                                                                       
!                B L A S  SUBPROGRAM                                    
!    DESCRIPTION OF PARAMETERS                                          
!                                                                       
!     --INPUT--                                                         
!        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)                       
!       SA  SINGLE PRECISION SCALE FACTOR                               
!       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS                     
!     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX                      
!                                                                       
!     --OUTPUT--                                                        
!       SX  SINGLE PRECISION RESULT (UNCHANGED IF N .LE. 0)             
!                                                                       
!     REPLACE SINGLE PRECISION SX BY SINGLE PRECISION SA*SX.            
!     FOR I = 0 TO N-1, REPLACE SX(1+I*INCX) WITH  SA * SX(1+I*INCX)    
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,     
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*, 
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL       
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323 
!***ROUTINES CALLED  (NONE)                                             
!***END PROLOGUE  SSCAL                                                 
!                                                                       
      REAL SA,SX(1) 
!***FIRST EXECUTABLE STATEMENT  SSCAL                                   
      IF(N.LE.0)RETURN 
      IF(INCX.EQ.1)GOTO 20 
!                                                                       
!        CODE FOR INCREMENTS NOT EQUAL TO 1.                            
!                                                                       
      NS = N*INCX 
          DO 10 I = 1,NS,INCX 
          SX(I) = SA*SX(I) 
   10     CONTINUE 
      RETURN 
!                                                                       
!        CODE FOR INCREMENTS EQUAL TO 1.                                
!                                                                       
!                                                                       
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.   
!                                                                       
   20 M = MOD(N,5) 
      IF( M .EQ. 0 ) GO TO 40 
      DO 30 I = 1,M 
        SX(I) = SA*SX(I) 
   30 END DO 
      IF( N .LT. 5 ) RETURN 
   40 MP1 = M + 1 
      DO 50 I = MP1,N,5 
        SX(I) = SA*SX(I) 
        SX(I + 1) = SA*SX(I + 1) 
        SX(I + 2) = SA*SX(I + 2) 
        SX(I + 3) = SA*SX(I + 3) 
        SX(I + 4) = SA*SX(I + 4) 
   50 END DO 
      RETURN 
      END                                           
      SUBROUTINE SSWAP(N,SX,INCX,SY,INCY) 
!***BEGIN PROLOGUE  SSWAP                                               
!***DATE WRITTEN   791001   (YYMMDD)                                    
!***REVISION DATE  820801   (YYMMDD)                                    
!***CATEGORY NO.  D1A5                                                  
!***KEYWORDS  BLAS,INTERCHANGE,LINEAR ALGEBRA,VECTOR                    
!***AUTHOR  LAWSON, C. L., (JPL)                                        
!           HANSON, R. J., (SNLA)                                       
!           KINCAID, D. R., (U. OF TEXAS)                               
!           KROGH, F. T., (JPL)                                         
!***PURPOSE  INTERCHANGE S.P VECTORS                                    
!***DESCRIPTION                                                         
!                                                                       
!                B L A S  SUBPROGRAM                                    
!    DESCRIPTION OF PARAMETERS                                          
!                                                                       
!     --INPUT--                                                         
!        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)                       
!       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS                     
!     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX                      
!       SY  SINGLE PRECISION VECTOR WITH N ELEMENTS                     
!     INCY  STORAGE SPACING BETWEEN ELEMENTS OF SY                      
!                                                                       
!     --OUTPUT--                                                        
!       SX  INPUT VECTOR SY (UNCHANGED IF N .LE. 0)                     
!       SY  INPUT VECTOR SX (UNCHANGED IF N .LE. 0)                     
!                                                                       
!     INTERCHANGE SINGLE PRECISION SX AND SINGLE PRECISION SY.          
!     FOR I = 0 TO N-1, INTERCHANGE  SX(LX+I*INCX) AND SY(LY+I*INCY),   
!     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS       
!     DEFINED IN A SIMILAR WAY USING INCY.                              
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,     
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*, 
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL       
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323 
!***ROUTINES CALLED  (NONE)                                             
!***END PROLOGUE  SSWAP                                                 
!                                                                       
      REAL SX(1),SY(1),STEMP1,STEMP2,STEMP3 
!***FIRST EXECUTABLE STATEMENT  SSWAP                                   
      IF(N.LE.0)RETURN 
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60 
    5 CONTINUE 
!                                                                       
!       CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.                     
!                                                                       
      IX = 1 
      IY = 1 
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1 
      DO 10 I = 1,N 
        STEMP1 = SX(IX) 
        SX(IX) = SY(IY) 
        SY(IY) = STEMP1 
        IX = IX + INCX 
        IY = IY + INCY 
   10 END DO 
      RETURN 
!                                                                       
!       CODE FOR BOTH INCREMENTS EQUAL TO 1                             
!                                                                       
!                                                                       
!       CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 3.    
!                                                                       
   20 M = MOD(N,3) 
      IF( M .EQ. 0 ) GO TO 40 
      DO 30 I = 1,M 
        STEMP1 = SX(I) 
        SX(I) = SY(I) 
        SY(I) = STEMP1 
   30 END DO 
      IF( N .LT. 3 ) RETURN 
   40 MP1 = M + 1 
      DO 50 I = MP1,N,3 
        STEMP1 = SX(I) 
        STEMP2 = SX(I+1) 
        STEMP3 = SX(I+2) 
        SX(I) = SY(I) 
        SX(I+1) = SY(I+1) 
        SX(I+2) = SY(I+2) 
        SY(I) = STEMP1 
        SY(I+1) = STEMP2 
        SY(I+2) = STEMP3 
   50 END DO 
      RETURN 
   60 CONTINUE 
!                                                                       
!     CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.                     
!                                                                       
      NS = N*INCX 
        DO 70 I=1,NS,INCX 
        STEMP1 = SX(I) 
        SX(I) = SY(I) 
        SY(I) = STEMP1 
   70   CONTINUE 
      RETURN 
      END                                           
      SUBROUTINE BNDACC(G,MDG,NB,IP,IR,MT,JT) 
!***BEGIN PROLOGUE  BNDACC                                              
!***DATE WRITTEN   790101   (YYMMDD)                                    
!***REVISION DATE  830513   (YYMMDD)                                    
!***CATEGORY NO.  D9                                                    
!***KEYWORDS  CURVE FITTING,LEAST SQUARE                                
!***AUTHOR  LAWSON, C. L., (JPL)                                        
!           HANSON, R. J., (SNLA)                                       
!***PURPOSE  SOLVE THE LEAST SQUARES PROBLEM AX = B FOR BANDED          
!            MATRICES A USING SEQUENTIAL ACCUMULATION OF ROWS OF        
!            THE DATA MATRIX.  EXACTLY ONE RIGHT-HANDED SIDE VECTOR     
!            IS PERMITTED.                                              
!***DESCRIPTION                                                         
!                                                                       
!     THESE SUBROUTINES SOLVE THE LEAST SQUARES PROBLEM AX = B FOR      
!     BANDED MATRICES A USING SEQUENTIAL ACCUMULATION OF ROWS OF THE    
!     DATA MATRIX.  EXACTLY ONE RIGHT-HAND SIDE VECTOR IS PERMITTED.    
!                                                                       
!     THESE SUBROUTINES ARE INTENDED FOR THE TYPE OF LEAST SQUARES      
!     SYSTEMS THAT ARISE IN APPLICATIONS SUCH AS CURVE OR SURFACE       
!     FITTING OF DATA.  THE LEAST SQUARES EQUATIONS ARE ACCUMULATED AND 
!     PROCESSED USING ONLY PART OF THE DATA.  THIS REQUIRES A CERTAIN   
!     USER INTERACTION DURING THE SOLUTION OF AX = B.                   
!                                                                       
!     SPECIFICALLY, SUPPOSE THE DATA MATRIX (A B) IS ROW PARTITIONED    
!     INTO Q SUBMATRICES.  LET (E F) BE THE T-TH ONE OF THESE           
!     SUBMATRICES WHERE E = (0 C 0).  HERE THE DIMENSION OF E IS MT BY N
!     AND THE DIMENSION OF C IS MT BY NB.  THE VALUE OF NB IS THE       
!     BANDWIDTH OF A.  THE DIMENSIONS OF THE LEADING BLOCK OF ZEROS IN E
!     ARE MT BY JT-1.                                                   
!                                                                       
!     THE USER OF THE SUBROUTINE BNDACC PROVIDES MT,JT,C AND F FOR      
!     T=1,...,Q.  NOT ALL OF THIS DATA MUST BE SUPPLIED AT ONCE.        
!                                                                       
!     FOLLOWING THE PROCESSING OF THE VARIOUS BLOCKS (E F), THE MATRIX  
!     (A B) HAS BEEN TRANSFORMED TO THE FORM (R D) WHERE R IS UPPER     
!     TRIANGULAR AND BANDED WITH BANDWIDTH NB.  THE LEAST SQUARES       
!     SYSTEM RX = D IS THEN EASILY SOLVED USING BACK SUBSTITUTION BY    
!     EXECUTING THE STATEMENT CALL BNDSOL(1,...). THE SEQUENCE OF       
!     VALUES FOR JT MUST BE NONDECREASING.  THIS MAY REQUIRE SOME       
!     PRELIMINARY INTERCHANGES OF ROWS AND COLUMNS OF THE MATRIX A.     
!                                                                       
!     THE PRIMARY REASON FOR THESE SUBROUTINES IS THAT THE TOTAL        
!     PROCESSING CAN TAKE PLACE IN A WORKING ARRAY OF DIMENSION MU BY   
!     NB+1.  AN ACCEPTABLE VALUE FOR MU IS                              
!                                                                       
!                       MU = MAX(MT + N + 1),                           
!                                                                       
!     WHERE N IS THE NUMBER OF UNKNOWNS.                                
!                                                                       
!     HERE THE MAXIMUM IS TAKEN OVER ALL VALUES OF MT FOR T=1,...,Q.    
!     NOTICE THAT MT CAN BE TAKEN TO BE A SMALL AS ONE, SHOWING THAT    
!     MU CAN BE AS SMALL AS N+2.  THE SUBPROGRAM BNDACC PROCESSES THE   
!     ROWS MORE EFFICIENTLY IF MU IS LARGE ENOUGH SO THAT EACH NEW      
!     BLOCK (C F) HAS A DISTINCT VALUE OF JT.                           
!                                                                       
!     THE FOUR PRINCIPLE PARTS OF THESE ALGORITHMS ARE OBTAINED BY THE  
!     FOLLOWING CALL STATEMENTS                                         
!                                                                       
!     CALL BNDACC(...)  INTRODUCE NEW BLOCKS OF DATA.                   
!                                                                       
!     CALL BNDSOL(1,...)COMPUTE SOLUTION VECTOR AND LENGTH OF           
!                       RESIDUAL VECTOR.                                
!                                                                       
!     CALL BNDSOL(2,...)GIVEN ANY ROW VECTOR H SOLVE YR = H FOR THE     
!                       ROW VECTOR Y.                                   
!                                                                       
!     CALL BNDSOL(3,...)GIVEN ANY COLUMN VECTOR W SOLVE RZ = W FOR      
!                       THE COLUMN VECTOR Z.                            
!                                                                       
!     THE DOTS IN THE ABOVE CALL STATEMENTS INDICATE ADDITIONAL         
!     ARGUMENTS THAT WILL BE SPECIFIED IN THE FOLLOWING PARAGRAPHS.     
!                                                                       
!     THE USER MUST DIMENSION THE ARRAY APPEARING IN THE CALL LIST..    
!     G(MDG,NB+1)                                                       
!                                                                       
!     DESCRIPTION OF CALLING SEQUENCE FOR BNDACC..                      
!                                                                       
!     THE ENTIRE SET OF PARAMETERS FOR BNDACC ARE                       
!                                                                       
!     INPUT..                                                           
!                                                                       
!     G(*,*)            THE WORKING ARRAY INTO WHICH THE USER WILL      
!                       PLACE THE MT BY NB+1 BLOCK (C F) IN ROWS IR     
!                       THROUGH IR+MT-1, COLUMNS 1 THROUGH NB+1.        
!                       SEE DESCRIPTIONS OF IR AND MT BELOW.            
!                                                                       
!     MDG               THE NUMBER OF ROWS IN THE WORKING ARRAY         
!                       G(*,*).  THE VALUE OF MDG SHOULD BE .GE. MU.    
!                       THE VALUE OF MU IS DEFINED IN THE ABSTRACT      
!                       OF THESE SUBPROGRAMS.                           
!                                                                       
!     NB                THE BANDWIDTH OF THE DATA MATRIX A.             
!                                                                       
!     IP                SET BY THE USER TO THE VALUE 1 BEFORE THE       
!                       FIRST CALL TO BNDACC.  ITS SUBSEQUENT VALUE     
!                       IS CONTROLLED BY BNDACC TO SET UP FOR THE       
!                       NEXT CALL TO BNDACC.                            
!                                                                       
!     IR                INDEX OF THE ROW OF G(*,*) WHERE THE USER IS    
!                       TO PLACE THE NEW BLOCK OF DATA (C F).  SET BY   
!                       THE USER TO THE VALUE 1 BEFORE THE FIRST CALL   
!                       TO BNDACC.  ITS SUBSEQUENT VALUE IS CONTROLLED  
!                       BY BNDACC. A VALUE OF IR .GT. MDG IS CONSIDERED 
!                       AN ERROR.                                       
!                                                                       
!     MT,JT             SET BY THE USER TO INDICATE RESPECTIVELY THE    
!                       NUMBER OF NEW ROWS OF DATA IN THE BLOCK AND     
!                       THE INDEX OF THE FIRST NONZERO COLUMN IN THAT   
!                       SET OF ROWS (E F) = (0 C 0 F) BEING PROCESSED.  
!                                                                       
!     OUTPUT..                                                          
!                                                                       
!     G(*,*)            THE WORKING ARRAY WHICH WILL CONTAIN THE        
!                       PROCESSED ROWS OF THAT PART OF THE DATA         
!                       MATRIX WHICH HAS BEEN PASSED TO BNDACC.         
!                                                                       
!     IP,IR             THE VALUES OF THESE ARGUMENTS ARE ADVANCED BY   
!                       BNDACC TO BE READY FOR STORING AND PROCESSING   
!                       A NEW BLOCK OF DATA IN G(*,*).                  
!                                                                       
!     DESCRIPTION OF CALLING SEQUENCE FOR BNDSOL..                      
!                                                                       
!     THE USER MUST DIMENSION THE ARRAYS APPEARING IN THE CALL LIST..   
!                                                                       
!     G(MDG,NB+1), X(N)                                                 
!                                                                       
!     THE ENTIRE SET OF PARAMETERS FOR BNDSOL ARE                       
!                                                                       
!     INPUT..                                                           
!                                                                       
!     MODE              SET BY THE USER TO ONE OF THE VALUES 1, 2, OR   
!                       3.  THESE VALUES RESPECTIVELY INDICATE THAT     
!                       THE SOLUTION OF AX = B, YR = H OR RZ = W IS     
!                       REQUIRED.                                       
!                                                                       
!     G(*,*),MDG,       THESE ARGUMENTS ALL HAVE THE SAME MEANING AND   
!      NB,IP,IR         CONTENTS AS FOLLOWING THE LAST CALL TO BNDACC.  
!                                                                       
!     X(*)              WITH MODE=2 OR 3 THIS ARRAY CONTAINS,           
!                       RESPECTIVELY, THE RIGHT-SIDE VECTORS H OR W OF  
!                       THE SYSTEMS YR = H OR RZ = W.                   
!                                                                       
!     N                 THE NUMBER OF VARIABLES IN THE SOLUTION         
!                       VECTOR.  IF ANY OF THE N DIAGONAL TERMS ARE     
!                       ZERO THE SUBROUTINE BNDSOL PRINTS AN            
!                       APPROPRIATE MESSAGE.  THIS CONDITION IS         
!                       CONSIDERED AN ERROR.                            
!                                                                       
!     OUTPUT..                                                          
!                                                                       
!     X(*)              THIS ARRAY CONTAINS THE SOLUTION VECTORS X,     
!                       Y OR Z OF THE SYSTEMS AX = B, YR = H OR         
!                       RZ = W DEPENDING ON THE VALUE OF MODE=1,        
!                       2 OR 3.                                         
!                                                                       
!     RNORM             IF MODE=1 RNORM IS THE EUCLIDEAN LENGTH OF THE  
!                       RESIDUAL VECTOR AX-B.  WHEN MODE=2 OR 3 RNORM   
!                       IS SET TO ZERO.                                 
!                                                                       
!     REMARKS..                                                         
!                                                                       
!     TO OBTAIN THE UPPER TRIANGULAR MATRIX AND TRANSFORMED RIGHT-HAND  
!     SIDE VECTOR D SO THAT THE SUPER DIAGONALS OF R FORM THE COLUMNS   
!     OF G(*,*), EXECUTE THE FOLLOWING FORTRAN STATEMENTS.              
!                                                                       
!     NBP1=NB+1                                                         
!                                                                       
!     DO 10 J=1, NBP1                                                   
!                                                                       
!  10 G(IR,J) = 0.E0                                                    
!                                                                       
!     MT=1                                                              
!                                                                       
!     JT=N+1                                                            
!                                                                       
!     CALL BNDACC(G,MDG,NB,IP,IR,MT,JT)                                 
!***REFERENCES  C. L. LAWSON AND R. J. HANSON,                          
!                 SOLVING LEAST SQUARE PROBLEMS,PRENCTICE-HALL, INC     
!                 (1974), CHAPTER 27                                    
!***ROUTINES CALLED  H12,XERROR                                         
!***END PROLOGUE  BNDACC                                                
      DIMENSION G(MDG,1) 
!***FIRST EXECUTABLE STATEMENT  BNDACC                                  
      ZERO=0. 
!                                                                       
!              ALG. STEPS 1-4 ARE PERFORMED EXTERNAL TO THIS SUBROUTINE.
!                                                                       
      NBP1=NB+1 
      IF (MT.LE.0.OR.NB.LE.0) RETURN 
!                                                                       
      IF(.NOT.MDG.LT.IR) GO TO 5 
      NERR=1 
      IOPT=2 
      CALL XERROR( 'BNDACC MDG.LT.IR.. PROBABLE ERROR.',34,NERR,IOPT) 
      RETURN 
    5 CONTINUE 
!                                                                       
!                                             ALG. STEP 5               
      IF (JT.EQ.IP) GO TO 70 
!                                             ALG. STEPS 6-7            
      IF (JT.LE.IR) GO TO 30 
!                                             ALG. STEPS 8-9            
      DO 10 I=1,MT 
        IG1=JT+MT-I 
        IG2=IR+MT-I 
        DO 10 J=1,NBP1 
        G(IG1,J)=G(IG2,J) 
   10 CONTINUE 
!                                             ALG. STEP 10              
      IE=JT-IR 
      DO 20 I=1,IE 
        IG=IR+I-1 
        DO 20 J=1,NBP1 
        G(IG,J)=ZERO 
   20 CONTINUE 
!                                             ALG. STEP 11              
      IR=JT 
!                                             ALG. STEP 12              
   30 MU=MIN0(NB-1,IR-IP-1) 
      IF (MU.EQ.0) GO TO 60 
!                                             ALG. STEP 13              
      DO 50 L=1,MU 
!                                             ALG. STEP 14              
        K=MIN0(L,JT-IP) 
!                                             ALG. STEP 15              
        LP1=L+1 
        IG=IP+L 
        DO 40 I=LP1,NB 
          JG=I-K 
          G(IG,JG)=G(IG,I) 
   40 END DO 
!                                             ALG. STEP 16              
        DO 50 I=1,K 
        JG=NBP1-I 
        G(IG,JG)=ZERO 
   50 CONTINUE 
!                                             ALG. STEP 17              
   60 IP=JT 
!                                             ALG. STEPS 18-19          
   70 MH=IR+MT-IP 
      KH=MIN0(NBP1,MH) 
!                                             ALG. STEP 20              
      DO 80 I=1,KH 
        CALL H12 (1,I,MAX0(I+1,IR-IP+1),MH,G(IP,I),1,RHO,               &
     &            G(IP,I+1),1,MDG,NBP1-I)                               
   80 END DO 
!                                             ALG. STEP 21              
      IR=IP+KH 
!                                             ALG. STEP 22              
      IF (KH.LT.NBP1) GO TO 100 
!                                             ALG. STEP 23              
      DO 90 I=1,NB 
        G(IR-1,I)=ZERO 
   90 END DO 
!                                             ALG. STEP 24              
  100 CONTINUE 
!                                             ALG. STEP 25              
      RETURN 
      END                                           
      SUBROUTINE BNDSOL(MODE,G,MDG,NB,IP,IR,X,N,RNORM) 
!***BEGIN PROLOGUE  BNDSOL                                              
!***DATE WRITTEN   790101   (YYMMDD)                                    
!***REVISION DATE  830513   (YYMMDD)                                    
!***CATEGORY NO.  D9                                                    
!***KEYWORDS  CURVE FITTING,LEAST SQUARE                                
!***AUTHOR  LAWSON, C. L., (JPL)                                        
!           HANSON, R. J., (SNLA)                                       
!***PURPOSE  SOLVE THE LEAST SQUARES PROBLEM AX = B FOR BANDED          
!            MATRICES A USING SEQUENTIAL ACCUMULATION OF ROWS OF        
!            THE DATA MATRIX.  EXACTLY ONE RIGHT-HANDED SIDE VECTOR     
!            IS PERMITTED.                                              
!***DESCRIPTION                                                         
!                                                                       
!     THESE SUBROUTINES SOLVE THE LEAST SQUARES PROBLEM AX = B FOR      
!     BANDED MATRICES A USING SEQUENTIAL ACCUMULATION OF ROWS OF THE    
!     DATA MATRIX.  EXACTLY ONE RIGHT-HAND SIDE VECTOR IS PERMITTED.    
!                                                                       
!     THESE SUBROUTINES ARE INTENDED FOR THE TYPE OF LEAST SQUARES      
!     SYSTEMS THAT ARISE IN APPLICATIONS SUCH AS CURVE OR SURFACE       
!     FITTING OF DATA.  THE LEAST SQUARES EQUATIONS ARE ACCUMULATED AND 
!     PROCESSED USING ONLY PART OF THE DATA.  THIS REQUIRES A CERTAIN   
!     USER INTERACTION DURING THE SOLUTION OF AX = B.                   
!                                                                       
!     SPECIFICALLY, SUPPOSE THE DATA MATRIX (A B) IS ROW PARTITIONED    
!     INTO Q SUBMATRICES.  LET (E F) BE THE T-TH ONE OF THESE           
!     SUBMATRICES WHERE E = (0 C 0).  HERE THE DIMENSION OF E IS MT BY N
!     AND THE DIMENSION OF C IS MT BY NB.  THE VALUE OF NB IS THE       
!     BANDWIDTH OF A.  THE DIMENSIONS OF THE LEADING BLOCK OF ZEROS IN E
!     ARE MT BY JT-1.                                                   
!                                                                       
!     THE USER OF THE SUBROUTINE BNDACC PROVIDES MT,JT,C AND F FOR      
!     T=1,...,Q.  NOT ALL OF THIS DATA MUST BE SUPPLIED AT ONCE.        
!                                                                       
!     FOLLOWING THE PROCESSING OF THE VARIOUS BLOCKS (E F), THE MATRIX  
!     (A B) HAS BEEN TRANSFORMED TO THE FORM (R D) WHERE R IS UPPER     
!     TRIANGULAR AND BANDED WITH BANDWIDTH NB.  THE LEAST SQUARES       
!     SYSTEM RX = D IS THEN EASILY SOLVED USING BACK SUBSTITUTION BY    
!     EXECUTING THE STATEMENT CALL BNDSOL(1,...). THE SEQUENCE OF       
!     VALUES FOR JT MUST BE NONDECREASING.  THIS MAY REQUIRE SOME       
!     PRELIMINARY INTERCHANGES OF ROWS AND COLUMNS OF THE MATRIX A.     
!                                                                       
!     THE PRIMARY REASON FOR THESE SUBROUTINES IS THAT THE TOTAL        
!     PROCESSING CAN TAKE PLACE IN A WORKING ARRAY OF DIMENSION MU BY   
!     NB+1.  AN ACCEPTABLE VALUE FOR MU IS                              
!                                                                       
!                       MU = MAX(MT + N + 1),                           
!                                                                       
!     WHERE N IS THE NUMBER OF UNKNOWNS.                                
!                                                                       
!     HERE THE MAXIMUM IS TAKEN OVER ALL VALUES OF MT FOR T=1,...,Q.    
!     NOTICE THAT MT CAN BE TAKEN TO BE A SMALL AS ONE, SHOWING THAT    
!     MU CAN BE AS SMALL AS N+2.  THE SUBPROGRAM BNDACC PROCESSES THE   
!     ROWS MORE EFFICIENTLY IF MU IS LARGE ENOUGH SO THAT EACH NEW      
!     BLOCK (C F) HAS A DISTINCT VALUE OF JT.                           
!                                                                       
!     THE FOUR PRINCIPLE PARTS OF THESE ALGORITHMS ARE OBTAINED BY THE  
!     FOLLOWING CALL STATEMENTS                                         
!                                                                       
!     CALL BNDACC(...)  INTRODUCE NEW BLOCKS OF DATA.                   
!                                                                       
!     CALL BNDSOL(1,...)COMPUTE SOLUTION VECTOR AND LENGTH OF           
!                       RESIDUAL VECTOR.                                
!                                                                       
!     CALL BNDSOL(2,...)GIVEN ANY ROW VECTOR H SOLVE YR = H FOR THE     
!                       ROW VECTOR Y.                                   
!                                                                       
!     CALL BNDSOL(3,...)GIVEN ANY COLUMN VECTOR W SOLVE RZ = W FOR      
!                       THE COLUMN VECTOR Z.                            
!                                                                       
!     THE DOTS IN THE ABOVE CALL STATEMENTS INDICATE ADDITIONAL         
!     ARGUMENTS THAT WILL BE SPECIFIED IN THE FOLLOWING PARAGRAPHS.     
!                                                                       
!     THE USER MUST DIMENSION THE ARRAY APPEARING IN THE CALL LIST..    
!     G(MDG,NB+1)                                                       
!                                                                       
!     DESCRIPTION OF CALLING SEQUENCE FOR BNDACC..                      
!                                                                       
!     THE ENTIRE SET OF PARAMETERS FOR BNDACC ARE                       
!                                                                       
!     INPUT..                                                           
!                                                                       
!     G(*,*)            THE WORKING ARRAY INTO WHICH THE USER WILL      
!                       PLACE THE MT BY NB+1 BLOCK (C F) IN ROWS IR     
!                       THROUGH IR+MT-1, COLUMNS 1 THROUGH NB+1.        
!                       SEE DESCRIPTIONS OF IR AND MT BELOW.            
!                                                                       
!     MDG               THE NUMBER OF ROWS IN THE WORKING ARRAY         
!                       G(*,*).  THE VALUE OF MDG SHOULD BE .GE. MU.    
!                       THE VALUE OF MU IS DEFINED IN THE ABSTRACT      
!                       OF THESE SUBPROGRAMS.                           
!                                                                       
!     NB                THE BANDWIDTH OF THE DATA MATRIX A.             
!                                                                       
!     IP                SET BY THE USER TO THE VALUE 1 BEFORE THE       
!                       FIRST CALL TO BNDACC.  ITS SUBSEQUENT VALUE     
!                       IS CONTROLLED BY BNDACC TO SET UP FOR THE       
!                       NEXT CALL TO BNDACC.                            
!                                                                       
!     IR                INDEX OF THE ROW OF G(*,*) WHERE THE USER IS    
!                       THE USER TO THE VALUE 1 BEFORE THE FIRST CALL   
!                       TO BNDACC.  ITS SUBSEQUENT VALUE IS CONTROLLED  
!                       BY BNDACC. A VALUE OF IR .GT. MDG IS CONSIDERED 
!                       AN ERROR.                                       
!                                                                       
!     MT,JT             SET BY THE USER TO INDICATE RESPECTIVELY THE    
!                       NUMBER OF NEW ROWS OF DATA IN THE BLOCK AND     
!                       THE INDEX OF THE FIRST NONZERO COLUMN IN THAT   
!                       SET OF ROWS (E F) = (0 C 0 F) BEING PROCESSED.  
!     OUTPUT..                                                          
!                                                                       
!     G(*,*)            THE WORKING ARRAY WHICH WILL CONTAIN THE        
!                       PROCESSED ROWS OF THAT PART OF THE DATA         
!                       MATRIX WHICH HAS BEEN PASSED TO BNDACC.         
!                                                                       
!     IP,IR             THE VALUES OF THESE ARGUMENTS ARE ADVANCED BY   
!                       BNDACC TO BE READY FOR STORING AND PROCESSING   
!                       A NEW BLOCK OF DATA IN G(*,*).                  
!                                                                       
!     DESCRIPTION OF CALLING SEQUENCE FOR BNDSOL..                      
!                                                                       
!     THE USER MUST DIMENSION THE ARRAYS APPEARING IN THE CALL LIST..   
!                                                                       
!     G(MDG,NB+1), X(N)                                                 
!                                                                       
!     THE ENTIRE SET OF PARAMETERS FOR BNDSOL ARE                       
!                                                                       
!     INPUT..                                                           
!                                                                       
!     MODE              SET BY THE USER TO ONE OF THE VALUES 1, 2, OR   
!                       3.  THESE VALUES RESPECTIVELY INDICATE THAT     
!                       THE SOLUTION OF AX = B, YR = H OR RZ = W IS     
!                       REQUIRED.                                       
!                                                                       
!     G(*,*),MDG,       THESE ARGUMENTS ALL HAVE THE SAME MEANING AND   
!      NB,IP,IR         CONTENTS AS FOLLOWING THE LAST CALL TO BNDACC.  
!                                                                       
!     X(*)              WITH MODE=2 OR 3 THIS ARRAY CONTAINS,           
!                       RESPECTIVELY, THE RIGHT-SIDE VECTORS H OR W OF  
!                       THE SYSTEMS YR = H OR RZ = W.                   
!                                                                       
!     N                 THE NUMBER OF VARIABLES IN THE SOLUTION         
!                       VECTOR.  IF ANY OF THE N DIAGONAL TERMS ARE     
!                       ZERO THE SUBROUTINE BNDSOL PRINTS AN            
!                       APPROPRIATE MESSAGE.  THIS CONDITION IS         
!                       CONSIDERED AN ERROR.                            
!                                                                       
!     OUTPUT..                                                          
!                                                                       
!     X(*)              THIS ARRAY CONTAINS THE SOLUTION VECTORS X,     
!                       Y OR Z OF THE SYSTEMS AX = B, YR = H OR         
!                       RZ = W DEPENDING ON THE VALUE OF MODE=1,        
!                       2 OR 3.                                         
!                                                                       
!     RNORM             IF MODE=1 RNORM IS THE EUCLIDEAN LENGTH OF THE  
!                       RESIDUAL VECTOR AX-B.  WHEN MODE=2 OR 3 RNORM   
!                       IS SET TO ZERO.                                 
!                                                                       
!     REMARKS..                                                         
!                                                                       
!     TO OBTAIN THE UPPER TRIANGULAR MATRIX AND TRANSFORMED RIGHT-HAND  
!     SIDE VECTOR D SO THAT THE SUPER DIAGONALS OF R FORM THE COLUMNS   
!     OF G(*,*), EXECUTE THE FOLLOWING FORTRAN STATEMENTS.              
!                                                                       
!     NBP1=NB+1                                                         
!                                                                       
!     DO 10 J=1, NBP1                                                   
!                                                                       
!  10 G(IR,J) = 0.E0                                                    
!                                                                       
!     MT=1                                                              
!                                                                       
!     JT=N+1                                                            
!                                                                       
!     CALL BNDACC(G,MDG,NB,IP,IR,MT,JT)                                 
!***REFERENCES  C. L. LAWSON AND R. J. HANSON,                          
!                 SOLVING LEAST SQUARE PROBLEMS,PRENCTICE-HALL, INC     
!                 (1974), CHAPTER 27                                    
!***ROUTINES CALLED  XERROR                                             
!***END PROLOGUE  BNDSOL                                                
      DIMENSION G(MDG,1),X(N) 
!***FIRST EXECUTABLE STATEMENT  BNDSOL                                  
      ZERO=0. 
!                                                                       
      RNORM=ZERO 
      GO TO (10,90,50), MODE 
!                                   ********************* MODE = 1      
!                                   ALG. STEP 26                        
   10      DO 20 J=1,N 
           X(J)=G(J,NB+1) 
   20 END DO 
      RSQ=ZERO 
      NP1=N+1 
      IRM1=IR-1 
      IF (NP1.GT.IRM1) GO TO 40 
           DO 30 J=NP1,IRM1 
           RSQ=RSQ+G(J,NB+1)**2 
   30 END DO 
      RNORM=SQRT(RSQ) 
   40 CONTINUE 
!                                   ********************* MODE = 3      
!                                   ALG. STEP 27                        
   50      DO 80 II=1,N 
           I=N+1-II 
!                                   ALG. STEP 28                        
           S=ZERO 
           L=MAX0(0,I-IP) 
!                                   ALG. STEP 29                        
           IF (I.EQ.N) GO TO 70 
!                                   ALG. STEP 30                        
           IE=MIN0(N+1-I,NB) 
                DO 60 J=2,IE 
                JG=J+L 
                IX=I-1+J 
                S=S+G(I,JG)*X(IX) 
   60 END DO 
!                                   ALG. STEP 31                        
   70      IF (G(I,L+1)) 80,130,80 
   80      X(I)=(X(I)-S)/G(I,L+1) 
!                                   ALG. STEP 32                        
      RETURN 
!                                   ********************* MODE = 2      
   90      DO 120 J=1,N 
           S=ZERO 
           IF (J.EQ.1) GO TO 110 
           I1=MAX0(1,J-NB+1) 
           I2=J-1 
                DO 100 I=I1,I2 
                L=J-I+1+MAX0(0,I-IP) 
                S=S+X(I)*G(I,L) 
  100 END DO 
  110      L=MAX0(0,J-IP) 
           IF (G(J,L+1)) 120,130,120 
  120      X(J)=(X(J)-S)/G(J,L+1) 
      RETURN 
!                                                                       
  130 CONTINUE 
      NERR=1 
      IOPT=2 
      CALL XERROR (  'BNDSOL A ZERO DIAGONAL TERM IS IN THE N BY N UPPER&
     & TRIANGULAR MATRIX.',69,NERR,IOPT)                                
      RETURN 
      END                                           
      SUBROUTINE BSPLVN(T,JHIGH,INDEX,X,ILEFT,VNIKX) 
!***BEGIN PROLOGUE  BSPLVN                                              
!***REFER TO  FC                                                        
!                                                                       
! CALCULATES THE VALUE OF ALL POSSIBLY NONZERO B-SPLINES AT *X* OF      
!  ORDER MAX(JHIGH,(J+1)(INDEX-1)) ON *T*.                              
!***ROUTINES CALLED  (NONE)                                             
!***END PROLOGUE  BSPLVN                                                
      DIMENSION T(1),VNIKX(1) 
      DIMENSION DELTAM(20),DELTAP(20) 
      DATA J/1/,(DELTAM(I),I=1,20),(DELTAP(I),I=1,20)/40*0./ 
!***FIRST EXECUTABLE STATEMENT  BSPLVN                                  
                                       GO TO (10,20),INDEX 
   10 J = 1 
      VNIKX(1) = 1. 
      IF (J .GE. JHIGH)                GO TO 99 
!                                                                       
   20    IPJ = ILEFT+J 
         DELTAP(J) = T(IPJ) - X 
         IMJP1 = ILEFT-J+1 
         DELTAM(J) = X - T(IMJP1) 
         VMPREV = 0. 
         JP1 = J+1 
         DO 26 L=1,J 
            JP1ML = JP1-L 
            VM = VNIKX(L)/(DELTAP(L) + DELTAM(JP1ML)) 
            VNIKX(L) = VM*DELTAP(L) + VMPREV 
   26       VMPREV = VM*DELTAM(JP1ML) 
         VNIKX(JP1) = VMPREV 
         J = JP1 
         IF (J .LT. JHIGH)             GO TO 20 
!                                                                       
   99                                  RETURN 
      END                                           
      SUBROUTINE EFCMN(NDATA,XDATA,YDATA,SDDATA,NORD,NBKPT,BKPTIN,      &
     &   MDEIN,MDEOUT,COEFF,BF,XTEMP,PTEMP,BKPT,G,MDG,W,MDW,LW)         
!***BEGIN PROLOGUE  EFCMN                                               
!***REFER TO  EFC                                                       
!                                                                       
!     THIS IS A COMPANION SUBPROGRAM TO EFC( ).                         
!     THIS SUBPROGRAM DOES WEIGHTED LEAST SQUARES FITTING               
!     OF DATA BY B-SPLINE CURVES.                                       
!     THE DOCUMENTATION FOR EFC( ) HAS MORE COMPLETE                    
!     USAGE INSTRUCTIONS.                                               
!                                                                       
!     REVISED 800905-1300                                               
!     REVISED YYMMDD-HHMM                                               
!     R. HANSON, SNLA 87185 SEPTEMBER, 1980.                            
!***ROUTINES CALLED  BNDACC,BNDSOL,BSPLVN,SCOPY,SSCAL,SSORT,XERROR,     
!                    XERRWV                                             
!***END PROLOGUE  EFCMN                                                 
      DIMENSION XDATA(NDATA), YDATA(NDATA), SDDATA(NDATA), BKPTIN(NBKPT) 
      DIMENSION COEFF(1), BF(NORD,NORD) 
      DIMENSION XTEMP(1), PTEMP(1), BKPT(NBKPT) 
      DIMENSION G(MDG,1), W(MDW,1) 
!***FIRST EXECUTABLE STATEMENT  EFCMN                                   
      ASSIGN 10 TO NPR001 
      GO TO 40 
   10 ASSIGN 20 TO NPR002 
      GO TO 100 
   20 ASSIGN 30 TO NPR003 
      GO TO 360 
   30 RETURN 
!     PROCEDURE (INITIALIZE-VARIABLES-AND-ANALYZE-INPUT)                
!                                                                       
!     PROCEDURE (INITIALIZE-VARIABLES-AND-ANALYZE-INPUT)                
   40 ZERO = 0.E0 
      ONE = 1.E0 
      L = NBKPT - NORD + 1 
!                                                                       
!     COMPUTE THE NUMBER OF VARIABLES.                                  
      N = NBKPT - NORD 
!                                                                       
!     INITIALLY SET ALL OUTPUT COEFFICIENTS TO ZERO.                    
      CALL SCOPY(N, ZERO, 0, COEFF, 1) 
      NP1 = L 
      IF (.NOT.(.NOT.(1.LE.NORD .AND. NORD.LE.20))) GO TO 50 
      CALL XERROR( 'EFC( ), THE ORDER OF THE B-SPLINE MUST BE 1 THRU 20.&
     &', 52, 3,   1)                                                    
      MDEOUT = -1 
      RETURN 
!                                                                       
   50 IF (.NOT.(NBKPT.LT.2*NORD)) GO TO 60 
      CALL XERROR( 'EFC( ), NUMBER OF KNOTS MUST BE AT LEAST TWICE THE B&
     &-SPLINE ORDER.', 66, 4, 1)                                        
      MDEOUT = -1 
      RETURN 
   60 IF (.NOT.(NDATA.LT.0)) GO TO 70 
      CALL XERROR( 'EFC( ), THE NUMBER OF DATA POINTS MUST BE NONNEGATIV&
     &E.', 54,    5, 1)                                                 
      MDEOUT = -1 
      RETURN 
   70 NB = (NBKPT-NORD+3)*(NORD+1) + (NBKPT+1)*(NORD+1) +               &
     & 2*MAX0(NBKPT,NDATA) + NBKPT + NORD**2                            
      IF (.NOT.(LW.LT.NB)) GO TO 80 
      CALL XERRWV( 'EFC( ). INSUFF. STORAGE FOR W(*). CHECK FORMULA THAT&
     & READS LW.GE. ... . (I1)=NEEDED, (I2)=GIVEN.', 96, 6, 1, 2, NB,   &
     & LW, 0, DUMMY, DUMMY)                                             
      MDEOUT = -1 
      RETURN 
   80 IF (.NOT.(.NOT.(MDEIN.EQ.1 .OR. MDEIN.EQ.2))) GO TO 90 
      CALL XERROR( 'EFC( ), INPUT VALUE OF MDEIN MUST BE 1-2.', 41, 7,  &
     & 1)                                                               
      MDEOUT = -1 
      RETURN 
!                                                                       
!     SORT THE BREAKPOINTS.                                             
   90 CALL SCOPY(NBKPT, BKPTIN, 1, BKPT, 1) 
      CALL SSORT(BKPT, DUMMY, NBKPT, 1) 
      XMIN = BKPT(NORD) 
!                                                                       
!     DEFINE THE INDEX OF RIGHT-MOST INTERVAL-DEFINING KNOT.            
      LAST = L 
      XMAX = BKPT(LAST) 
      NORDM1 = NORD - 1 
      NORDP1 = NORD + 1 
      GO TO NPR001, (10) 
!     PROCEDURE (PROCESS-LEAST-SQUARES-EQUATIONS)                       
!                                                                       
!     SORT DATA AND AN ARRAY OF POINTERS.                               
  100 CALL SCOPY(NDATA, XDATA, 1, XTEMP, 1) 
      I = 1 
      N20019 = NDATA 
      GO TO 120 
  110 I = I + 1 
  120 IF ((N20019-I).LT.0) GO TO 130 
      PTEMP(I) = I 
      GO TO 110 
  130 IF (.NOT.(NDATA.GT.0)) GO TO 140 
      CALL SSORT(XTEMP, PTEMP, NDATA, 2) 
      XMIN = AMIN1(XMIN,XTEMP(1)) 
      XMAX = AMAX1(XMAX,XTEMP(NDATA)) 
!                                                                       
!     FIX BREAKPOINT ARRAY IF NEEDED. THIS SHOULD ONLY INVOLVE VERY     
!     MINOR DIFFERENCES WITH THE INPUT ARRAY OF BREAKPOINTS.            
  140 I = 1 
      N20026 = NORD 
      GO TO 160 
  150 I = I + 1 
  160 IF ((N20026-I).LT.0) GO TO 170 
      BKPT(I) = AMIN1(BKPT(I),XMIN) 
      GO TO 150 
  170 I = LAST 
      N20030 = NBKPT 
      GO TO 190 
  180 I = I + 1 
  190 IF ((N20030-I).LT.0) GO TO 200 
      BKPT(I) = AMAX1(BKPT(I),XMAX) 
      GO TO 180 
!                                                                       
!     INITIALIZE PARAMETERS OF BANDED MATRIX PROCESSOR, BNDACC( ).      
  200 MT = 0 
      IP = 1 
      IR = 1 
      IDATA = 1 
      ILEFT = NORD 
      INTSEQ = 1 
  210 IF (.NOT.(IDATA.LE.NDATA)) GO TO 280 
!                                                                       
!     SORTED INDICES ARE IN PTEMP(*).                                   
      L = PTEMP(IDATA) 
      XVAL = XDATA(L) 
!                                                                       
!     WHEN INTERVAL CHANGES, PROCESS EQUATIONS IN THE LAST BLOCK.       
      IF (.NOT.(XVAL.GE.BKPT(ILEFT+1))) GO TO 250 
      INTRVL = ILEFT - NORDM1 
      CALL BNDACC(G, MDG, NORD, IP, IR, MT, INTRVL) 
      MT = 0 
!                                                                       
!     MOVE POINTER UP TO HAVE BKPT(ILEFT).LE.XVAL, ILEFT.LT.LAST.       
  220 IF (.NOT.(XVAL.GE.BKPT(ILEFT+1) .AND. ILEFT.LT.LAST-1)) GO TO 240 
      IF (.NOT.(MDEIN.EQ.2)) GO TO 230 
!                                                                       
!     DATA IS BEING SEQUENTIALLY ACCUMULATED. TRANSFER                  
!     PREVIOUSLY ACCUMULATED ROWS FROM W(*,*) TO G(*,*)                 
!     AND PROCESS THEM.                                                 
      CALL SCOPY(NORDP1, W(INTSEQ,1), MDW, G(IR,1), MDG) 
      CALL BNDACC(G, MDG, NORD, IP, IR, 1, INTSEQ) 
      INTSEQ = INTSEQ + 1 
  230 ILEFT = ILEFT + 1 
      GO TO 220 
  240 CONTINUE 
!                                                                       
!     OBTAIN B-SPLINE FUNCTION VALUE.                                   
  250 CALL BSPLVN(BKPT, NORD, 1, XVAL, ILEFT, BF) 
!                                                                       
!     MOVE ROW INTO PLACE.                                              
      IROW = IR + MT 
      MT = MT + 1 
      CALL SCOPY(NORD, BF, 1, G(IROW,1), MDG) 
      G(IROW,NORDP1) = YDATA(L) 
!                                                                       
!     SCALE DATA IF UNCERTAINTY IS NONZERO.                             
      IF (.NOT.(SDDATA(L).NE.ZERO)) GO TO 260 
      CALL SSCAL(NORDP1, ONE/SDDATA(L), G(IROW,1), MDG) 
!                                                                       
!     WHEN STAGING WORK AREA IS EXHAUSTED, PROCESS ROWS.                
  260 IF (.NOT.(IROW.EQ.MDG-1)) GO TO 270 
      INTRVL = ILEFT - NORDM1 
      CALL BNDACC(G, MDG, NORD, IP, IR, MT, INTRVL) 
      MT = 0 
  270 IDATA = IDATA + 1 
      GO TO 210 
!                                                                       
!     PROCESS LAST BLOCK OF EQUATIONS.                                  
  280 INTRVL = ILEFT - NORDM1 
      CALL BNDACC(G, MDG, NORD, IP, IR, MT, INTRVL) 
!                                                                       
!     FINISH PROCESSING ANY PREVIOUSLY ACCUMULATED                      
!     ROWS FROM W(*,*) TO G(*,*).                                       
      IF (.NOT.(MDEIN.EQ.2)) GO TO 320 
  290 CALL SCOPY(NORDP1, W(INTSEQ,1), MDW, G(IR,1), MDG) 
      CALL BNDACC(G, MDG, NORD, IP, IR, 1, MIN0(N,INTSEQ)) 
      IF (.NOT.(INTSEQ.EQ.NP1)) GO TO 300 
      GO TO 310 
  300 INTSEQ = INTSEQ + 1 
      GO TO 290 
  310 CONTINUE 
!                                                                       
!     LAST CALL TO ADJUST BLOCK POSITIONING.                            
  320 G(IR,1) = ZERO 
      CALL SCOPY(NORDP1, G(IR,1), 0, G(IR,1), MDG) 
      CALL BNDACC(G, MDG, NORD, IP, IR, 1, NP1) 
!                                                                       
!     TRANSFER ACCUMULATED ROWS FROM G(*,*) TO W(*,*) FOR               
!     POSSIBLE LATER SEQUENTIAL ACCUMULATION.                           
      I = 1 
      N20058 = NP1 
      GO TO 340 
  330 I = I + 1 
  340 IF ((N20058-I).LT.0) GO TO 350 
      CALL SCOPY(NORDP1, G(I,1), MDG, W(I,1), MDW) 
      GO TO 330 
  350 GO TO NPR002, (20) 
!     PROCEDURE (SOLVE-FOR-COEFFICIENTS-WHEN-POSSIBLE)                  
  360 I = 1 
      N20062 = N 
      GO TO 380 
  370 I = I + 1 
  380 IF ((N20062-I).LT.0) GO TO 400 
      IF (.NOT.(G(I,1).EQ.ZERO)) GO TO 390 
      MDEOUT = 2 
      GO TO 410 
  390 CONTINUE 
      GO TO 370 
!                                                                       
!     ALL THE DIAGONAL TERMS IN THE ACCUMULATED TRIANGULAR              
!     MATRIX ARE NONZERO.  THE SOLN. CAN BE COMPUTED BUT                
!     IT MAY BE UNSUITABLE FOR FURTHER USE DUE TO POOR                  
!     CONDITIONING OR THE LACK OF CONSTRAINTS.  NO CHECKING             
!     FOR EITHER OF THESE IS DONE HERE.                                 
  400 MODE = 1 
      CALL BNDSOL(MODE, G, MDG, NORD, IP, IR, COEFF, N, RNORM) 
      MDEOUT = 1 
  410 GO TO NPR003, (30) 
      END                                           
      SUBROUTINE H12(MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV) 
!***BEGIN PROLOGUE  H12                                                 
!***REFER TO  HFTI,LSEI,WNNLS                                           
!                                                                       
!     SUBROUTINE H12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)          
!                                                                       
!     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUN 12 
!     TO APPEAR IN 'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974
!                                                                       
!     MODIFIED AT SANDIA LABS, MAY 1977, TO --                          
!                                                                       
!     1)  REMOVE DOUBLE PRECISION ACCUMULATION, AND                     
!     2)  INCLUDE USAGE OF THE BASIC LINEAR ALGEBRA PACKAGE FOR         
!         VECTORS LONGER THAN A PARTICULAR THRESHOLD.                   
!                                                                       
!     CONSTRUCTION AND/OR APPLICATION OF A SINGLE                       
!     HOUSEHOLDER TRANSFORMATION..     Q = I + U*(U**T)/B               
!                                                                       
!     MODE    = 1 OR 2   TO SELECT ALGORITHM  H1  OR  H2 .              
!     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT.                         
!     L1,M   IF L1 .LE. M   THE TRANSFORMATION WILL BE CONSTRUCTED TO   
!            ZERO ELEMENTS INDEXED FROM L1 THROUGH M.   IF L1 GT. M     
!            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.            
!     U(),IUE,UP    ON ENTRY TO H1 U() CONTAINS THE PIVOT VECTOR.       
!                   IUE IS THE STORAGE INCREMENT BETWEEN ELEMENTS.      
!                                       ON EXIT FROM H1 U() AND UP      
!                   CONTAIN QUANTITIES DEFINING THE VECTOR U OF THE     
!                   HOUSEHOLDER TRANSFORMATION.   ON ENTRY TO H2 U()    
!                   AND UP SHOULD CONTAIN QUANTITIES PREVIOUSLY COMPUTED
!                   BY H1.  THESE WILL NOT BE MODIFIED BY H2.           
!     C()    ON ENTRY TO H1 OR H2 C() CONTAINS A MATRIX WHICH WILL BE   
!            REGARDED AS A SET OF VECTORS TO WHICH THE HOUSEHOLDER      
!            TRANSFORMATION IS TO BE APPLIED.  ON EXIT C() CONTAINS THE 
!            SET OF TRANSFORMED VECTORS.                                
!     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C().      
!     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C().                  
!     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED. IF NCV .LE. 0  
!            NO OPERATIONS WILL BE DONE ON C().                         
!***ROUTINES CALLED  SAXPY,SDOT,SSWAP                                   
!***END PROLOGUE  H12                                                   
      DIMENSION U(IUE,M), C(1) 
!***FIRST EXECUTABLE STATEMENT  H12                                     
      ONE=1. 
!                                                                       
      IF (0.GE.LPIVOT.OR.LPIVOT.GE.L1.OR.L1.GT.M) RETURN 
      CL=ABS(U(1,LPIVOT)) 
      IF (MODE.EQ.2) GO TO 60 
!                            ****** CONSTRUCT THE TRANSFORMATION. ******
          DO 10 J=L1,M 
   10     CL=AMAX1(ABS(U(1,J)),CL) 
      IF (CL) 130,130,20 
   20 CLINV=ONE/CL 
      SM=(U(1,LPIVOT)*CLINV)**2 
          DO 30 J=L1,M 
   30     SM=SM+(U(1,J)*CLINV)**2 
      CL=CL*SQRT(SM) 
      IF (U(1,LPIVOT)) 50,50,40 
   40 CL=-CL 
   50 UP=U(1,LPIVOT)-CL 
      U(1,LPIVOT)=CL 
      GO TO 70 
!            ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ******
!                                                                       
   60 IF (CL) 130,130,70 
   70 IF (NCV.LE.0) RETURN 
      B=UP*U(1,LPIVOT) 
!                       B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN.
!                                                                       
      IF (B) 80,130,130 
   80 B=ONE/B 
      MML1P2=M-L1+2 
      IF (MML1P2.GT.20) GO TO 140 
      I2=1-ICV+ICE*(LPIVOT-1) 
      INCR=ICE*(L1-LPIVOT) 
          DO 120 J=1,NCV 
          I2=I2+ICV 
          I3=I2+INCR 
          I4=I3 
          SM=C(I2)*UP 
              DO 90 I=L1,M 
              SM=SM+C(I3)*U(1,I) 
   90         I3=I3+ICE 
          IF (SM) 100,120,100 
  100     SM=SM*B 
          C(I2)=C(I2)+SM*UP 
              DO 110 I=L1,M 
              C(I4)=C(I4)+SM*U(1,I) 
  110         I4=I4+ICE 
  120     CONTINUE 
  130 RETURN 
  140 CONTINUE 
      L1M1=L1-1 
      KL1=1+(L1M1-1)*ICE 
      KL2=KL1 
      KLP=1+(LPIVOT-1)*ICE 
      UL1M1=U(1,L1M1) 
      U(1,L1M1)=UP 
      IF (LPIVOT.EQ.L1M1) GO TO 150 
      CALL SSWAP(NCV,C(KL1),ICV,C(KLP),ICV) 
  150 CONTINUE 
          DO 160 J=1,NCV 
          SM=SDOT(MML1P2,U(1,L1M1),IUE,C(KL1),ICE) 
          SM=SM*B 
          CALL SAXPY (MML1P2,SM,U(1,L1M1),IUE,C(KL1),ICE) 
          KL1=KL1+ICV 
  160 END DO 
      U(1,L1M1)=UL1M1 
      IF (LPIVOT.EQ.L1M1) RETURN 
      KL1=KL2 
      CALL SSWAP(NCV,C(KL1),ICV,C(KLP),ICV) 
      RETURN 
      END                                           
      SUBROUTINE SSORT(X,Y,N,KFLAG) 
!***BEGIN PROLOGUE  SSORT                                               
!***DATE WRITTEN   761101   (YYMMDD)                                    
!***REVISION DATE  820801   (YYMMDD)                                    
!***CATEGORY NO.  N6A2B1                                                
!***KEYWORDS  QUICKSORT,SINGLETON QUICKSORT,SORT,SORTING                
!***AUTHOR  JONES, R. E., (SNLA)                                        
!           WISNIEWSKI, J. A., (SNLA)                                   
!***PURPOSE  SSORT SORTS ARRAY X AND OPTIONALLY MAKES THE SAME          
!            INTERCHANGES IN ARRAY Y.  THE ARRAY X MAY BE SORTED IN     
!            INCREASING ORDER OR DECREASING ORDER.  A SLIGHTLY MODIFIED 
!            QUICKSORT ALGORITHM IS USED.                               
!***DESCRIPTION                                                         
!                                                                       
!     WRITTEN BY RONDALL E. JONES                                       
!     MODIFIED BY JOHN A. WISNIEWSKI TO USE THE SINGLETON QUICKSORT     
!     ALGORITHM.  DATE 18 NOVEMBER 1976.                                
!                                                                       
!     ABSTRACT                                                          
!         SSORT SORTS ARRAY X AND OPTIONALLY MAKES THE SAME             
!         INTERCHANGES IN ARRAY Y.  THE ARRAY X MAY BE SORTED IN        
!         INCREASING ORDER OR DECREASING ORDER.  A SLIGHTLY MODIFIED    
!         QUICKSORT ALGORITHM IS USED.                                  
!                                                                       
!     REFERENCE                                                         
!         SINGLETON, R. C., ALGORITHM 347, AN EFFICIENT ALGORITHM FOR   
!         SORTING WITH MINIMAL STORAGE, CACM,12(3),1969,185-7.          
!                                                                       
!     DESCRIPTION OF PARAMETERS                                         
!         X - ARRAY OF VALUES TO BE SORTED   (USUALLY ABSCISSAS)        
!         Y - ARRAY TO BE (OPTIONALLY) CARRIED ALONG                    
!         N - NUMBER OF VALUES IN ARRAY X TO BE SORTED                  
!         KFLAG - CONTROL PARAMETER                                     
!             =2  MEANS SORT X IN INCREASING ORDER AND CARRY Y ALONG.   
!             =1  MEANS SORT X IN INCREASING ORDER (IGNORING Y)         
!             =-1 MEANS SORT X IN DECREASING ORDER (IGNORING Y)         
!             =-2 MEANS SORT X IN DECREASING ORDER AND CARRY Y ALONG.   
!***REFERENCES  SINGLETON,R.C., ALGORITHM 347, AN EFFICIENT ALGORITHM   
!                 FOR SORTING WITH MINIMAL STORAGE, CACM,12(3),1969,    
!                 185-7.                                                
!***ROUTINES CALLED  XERROR                                             
!***END PROLOGUE  SSORT                                                 
      DIMENSION X(N),Y(N),IL(21),IU(21) 
!***FIRST EXECUTABLE STATEMENT  SSORT                                   
      NN = N 
      IF (NN.GE.1) GO TO 10 
      CALL XERROR ( 'SSORT- THE NUMBER OF VALUES TO BE SORTED WAS NOT PO&
     &SITIVE.',58,1,1)                                                  
      RETURN 
   10 KK = IABS(KFLAG) 
      IF ((KK.EQ.1).OR.(KK.EQ.2)) GO TO 15 
      CALL XERROR ( 'SSORT- THE SORT CONTROL PARAMETER, K, WAS NOT 2, 1,&
     & -1, OR -2.',62,2,1)                                              
      RETURN 
!                                                                       
! ALTER ARRAY X TO GET DECREASING ORDER IF NEEDED                       
!                                                                       
   15 IF (KFLAG.GE.1) GO TO 30 
      DO 20 I=1,NN 
   20 X(I) = -X(I) 
   30 GO TO (100,200),KK 
!                                                                       
! SORT X ONLY                                                           
!                                                                       
  100 CONTINUE 
      M=1 
      I=1 
      J=NN 
      R=.375 
  110 IF (I .EQ. J) GO TO 155 
  115 IF (R .GT. .5898437) GO TO 120 
      R=R+3.90625E-2 
      GO TO 125 
  120 R=R-.21875 
  125 K=I 
!                                  SELECT A CENTRAL ELEMENT OF THE      
!                                  ARRAY AND SAVE IT IN LOCATION T      
      IJ = I + IFIX (FLOAT (J-I) * R) 
      T=X(IJ) 
!                                  IF FIRST ELEMENT OF ARRAY IS GREATER 
!                                  THAN T, INTERCHANGE WITH T           
      IF (X(I) .LE. T) GO TO 130 
      X(IJ)=X(I) 
      X(I)=T 
      T=X(IJ) 
  130 L=J 
!                                  IF LAST ELEMENT OF ARRAY IS LESS THAN
!                                  T, INTERCHANGE WITH T                
      IF (X(J) .GE. T) GO TO 140 
      X(IJ)=X(J) 
      X(J)=T 
      T=X(IJ) 
!                                  IF FIRST ELEMENT OF ARRAY IS GREATER 
!                                  THAN T, INTERCHANGE WITH T           
      IF (X(I) .LE. T) GO TO 140 
      X(IJ)=X(I) 
      X(I)=T 
      T=X(IJ) 
      GO TO 140 
  135 TT=X(L) 
      X(L)=X(K) 
      X(K)=TT 
!                                  FIND AN ELEMENT IN THE SECOND HALF OF
!                                  THE ARRAY WHICH IS SMALLER THAN T    
  140 L=L-1 
      IF (X(L) .GT. T) GO TO 140 
!                                  FIND AN ELEMENT IN THE FIRST HALF OF 
!                                  THE ARRAY WHICH IS GREATER THAN T    
  145 K=K+1 
      IF (X(K) .LT. T) GO TO 145 
!                                  INTERCHANGE THESE ELEMENTS           
      IF (K .LE. L) GO TO 135 
!                                  SAVE UPPER AND LOWER SUBSCRIPTS OF   
!                                  THE ARRAY YET TO BE SORTED           
      IF (L-I .LE. J-K) GO TO 150 
      IL(M)=I 
      IU(M)=L 
      I=K 
      M=M+1 
      GO TO 160 
  150 IL(M)=K 
      IU(M)=J 
      J=L 
      M=M+1 
      GO TO 160 
!                                  BEGIN AGAIN ON ANOTHER PORTION OF    
!                                  THE UNSORTED ARRAY                   
  155 M=M-1 
      IF (M .EQ. 0) GO TO 300 
      I=IL(M) 
      J=IU(M) 
  160 IF (J-I .GE. 1) GO TO 125 
      IF (I .EQ. 1) GO TO 110 
      I=I-1 
  165 I=I+1 
      IF (I .EQ. J) GO TO 155 
      T=X(I+1) 
      IF (X(I) .LE. T) GO TO 165 
      K=I 
  170 X(K+1)=X(K) 
      K=K-1 
      IF (T .LT. X(K)) GO TO 170 
      X(K+1)=T 
      GO TO 165 
!                                                                       
! SORT X AND CARRY Y ALONG                                              
!                                                                       
  200 CONTINUE 
      M=1 
      I=1 
      J=NN 
      R=.375 
  210 IF (I .EQ. J) GO TO 255 
  215 IF (R .GT. .5898437) GO TO 220 
      R=R+3.90625E-2 
      GO TO 225 
  220 R=R-.21875 
  225 K=I 
!                                  SELECT A CENTRAL ELEMENT OF THE      
!                                  ARRAY AND SAVE IT IN LOCATION T      
      IJ = I + IFIX (FLOAT (J-I) *R) 
      T=X(IJ) 
      TY= Y(IJ) 
!                                  IF FIRST ELEMENT OF ARRAY IS GREATER 
!                                  THAN T, INTERCHANGE WITH T           
      IF (X(I) .LE. T) GO TO 230 
      X(IJ)=X(I) 
      X(I)=T 
      T=X(IJ) 
       Y(IJ)= Y(I) 
       Y(I)=TY 
      TY= Y(IJ) 
  230 L=J 
!                                  IF LAST ELEMENT OF ARRAY IS LESS THAN
!                                  T, INTERCHANGE WITH T                
      IF (X(J) .GE. T) GO TO 240 
      X(IJ)=X(J) 
      X(J)=T 
      T=X(IJ) 
       Y(IJ)= Y(J) 
       Y(J)=TY 
      TY= Y(IJ) 
!                                  IF FIRST ELEMENT OF ARRAY IS GREATER 
!                                  THAN T, INTERCHANGE WITH T           
      IF (X(I) .LE. T) GO TO 240 
      X(IJ)=X(I) 
      X(I)=T 
      T=X(IJ) 
       Y(IJ)= Y(I) 
       Y(I)=TY 
      TY= Y(IJ) 
      GO TO 240 
  235 TT=X(L) 
      X(L)=X(K) 
      X(K)=TT 
      TTY= Y(L) 
       Y(L)= Y(K) 
       Y(K)=TTY 
!                                  FIND AN ELEMENT IN THE SECOND HALF OF
!                                  THE ARRAY WHICH IS SMALLER THAN T    
  240 L=L-1 
      IF (X(L) .GT. T) GO TO 240 
!                                  FIND AN ELEMENT IN THE FIRST HALF OF 
!                                  THE ARRAY WHICH IS GREATER THAN T    
  245 K=K+1 
      IF (X(K) .LT. T) GO TO 245 
!                                  INTERCHANGE THESE ELEMENTS           
      IF (K .LE. L) GO TO 235 
!                                  SAVE UPPER AND LOWER SUBSCRIPTS OF   
!                                  THE ARRAY YET TO BE SORTED           
      IF (L-I .LE. J-K) GO TO 250 
      IL(M)=I 
      IU(M)=L 
      I=K 
      M=M+1 
      GO TO 260 
  250 IL(M)=K 
      IU(M)=J 
      J=L 
      M=M+1 
      GO TO 260 
!                                  BEGIN AGAIN ON ANOTHER PORTION OF    
!                                  THE UNSORTED ARRAY                   
  255 M=M-1 
      IF (M .EQ. 0) GO TO 300 
      I=IL(M) 
      J=IU(M) 
  260 IF (J-I .GE. 1) GO TO 225 
      IF (I .EQ. 1) GO TO 210 
      I=I-1 
  265 I=I+1 
      IF (I .EQ. J) GO TO 255 
      T=X(I+1) 
      TY= Y(I+1) 
      IF (X(I) .LE. T) GO TO 265 
      K=I 
  270 X(K+1)=X(K) 
       Y(K+1)= Y(K) 
      K=K-1 
      IF (T .LT. X(K)) GO TO 270 
      X(K+1)=T 
       Y(K+1)=TY 
      GO TO 265 
!                                                                       
! CLEAN UP                                                              
!                                                                       
  300 IF (KFLAG.GE.1) RETURN 
      DO 310 I=1,NN 
  310 X(I) = -X(I) 
      RETURN 
      END                                           
      FUNCTION BVALU(T,A,N,K,IDERIV,X,INBV,WORK) 
!***BEGIN PROLOGUE  BVALU                                               
!***DATE WRITTEN   800901   (YYMMDD)                                    
!***REVISION DATE  820801   (YYMMDD)                                    
!***CATEGORY NO.  E3,K6                                                 
!***KEYWORDS  B-SPLINE,DATA FITTING,INTERPOLATION,SPLINE                
!***AUTHOR  AMOS, D. E., (SNLA)                                         
!***PURPOSE  EVALUATES THE B-REPRESENTATION OF A B-SPLINE AT X FOR THE  
!            FUNCTION VALUE OR ANY OF ITS DERIVATIVES.                  
!***DESCRIPTION                                                         
!                                                                       
!     WRITTEN BY CARL DE BOOR AND MODIFIED BY D. E. AMOS                
!                                                                       
!     REFERENCE                                                         
!         SIAM J. NUMERICAL ANALYSIS, 14, NO. 3, JUNE, 1977, PP.441-472.
!                                                                       
!     ABSTRACT                                                          
!         BVALU IS THE BVALUE FUNCTION OF THE REFERENCE.                
!                                                                       
!         BVALU EVALUATES THE B-REPRESENTATION (T,A,N,K) OF A B-SPLINE  
!         AT X FOR THE FUNCTION VALUE ON IDERIV = 0 OR ANY OF ITS       
!         DERIVATIVES ON IDERIV = 1,2,...,K-1.  RIGHT LIMITING VALUES   
!         (RIGHT DERIVATIVES) ARE RETURNED EXCEPT AT THE RIGHT END      
!         POINT X=T(N+1) WHERE LEFT LIMITING VALUES ARE COMPUTED.  THE  
!         SPLINE IS DEFINED ON T(K) .LE. X .LE. T(N+1).  BVALU RETURNS  
!         A FATAL ERROR MESSAGE WHEN X IS OUTSIDE OF THIS INTERVAL.     
!                                                                       
!         TO COMPUTE LEFT DERIVATIVES OR LEFT LIMITING VALUES AT A      
!         KNOT T(I), REPLACE N BY I-1 AND SET X=T(I), I=K+1,N+1.        
!                                                                       
!         BVALU CALLS INTRV                                             
!                                                                       
!     DESCRIPTION OF ARGUMENTS                                          
!         INPUT                                                         
!          T       - KNOT VECTOR OF LENGTH N+K                          
!          A       - B-SPLINE COEFFICIENT VECTOR OF LENGTH N            
!          N       - NUMBER OF B-SPLINE COEFFICIENTS                    
!                    N = SUM OF KNOT MULTIPLICITIES-K                   
!          K       - ORDER OF THE B-SPLINE, K .GE. 1                    
!          IDERIV  - ORDER OF THE DERIVATIVE, 0 .LE. IDERIV .LE. K-1    
!                    IDERIV=0 RETURNS THE B-SPLINE VALUE                
!          X       - ARGUMENT, T(K) .LE. X .LE. T(N+1)                  
!          INBV    - AN INITIALIZATION PARAMETER WHICH MUST BE SET      
!                    TO 1 THE FIRST TIME BVALU IS CALLED.               
!                                                                       
!         OUTPUT                                                        
!          INBV    - INBV CONTAINS INFORMATION FOR EFFICIENT PROCESS-   
!                    ING AFTER THE INITIAL CALL AND INBV MUST NOT       
!                    BE CHANGED BY THE USER.  DISTINCT SPLINES REQUIRE  
!                    DISTINCT INBV PARAMETERS.                          
!          WORK    - WORK VECTOR OF LENGTH 3*K.                         
!          BVALU   - VALUE OF THE IDERIV-TH DERIVATIVE AT X             
!                                                                       
!     ERROR CONDITIONS                                                  
!         AN IMPROPER INPUT IS A FATAL ERROR                            
!***REFERENCES  C. DE BOOR, *PACKAGE FOR CALCULATING WITH B-SPLINES*,   
!                 SIAM JOURNAL ON NUMERICAL ANALYSIS, VOLUME 14, NO. 3, 
!                 JUNE 1977, PP. 441-472.                               
!***ROUTINES CALLED  INTRV,XERROR                                       
!***END PROLOGUE  BVALU                                                 
!                                                                       
!                                                                       
      INTEGER I,IDERIV,IDERP1,IHI,IHMKMJ,ILO,IMK,IMKPJ, INBV, IPJ,      &
     & IP1, IP1MJ, J, JJ, J1, J2, K, KMIDER, KMJ, KM1, KPK, MFLAG, N    
      REAL A, FKMJ, T, WORK, X 
!     DIMENSION T(N+K), WORK(3*K)                                       
      DIMENSION T(1), A(N), WORK(1) 
!***FIRST EXECUTABLE STATEMENT  BVALU                                   
      BVALU = 0.0E0 
      IF(K.LT.1) GO TO 102 
      IF(N.LT.K) GO TO 101 
      IF(IDERIV.LT.0 .OR. IDERIV.GE.K) GO TO 110 
      KMIDER = K - IDERIV 
!                                                                       
! *** FIND *I* IN (K,N) SUCH THAT T(I) .LE. X .LT. T(I+1)               
!     (OR, .LE. T(I+1) IF T(I) .LT. T(I+1) = T(N+1)).                   
      KM1 = K - 1 
      CALL INTRV(T, N+1, X, INBV, I, MFLAG) 
      IF (X.LT.T(K)) GO TO 120 
      IF (MFLAG.EQ.0) GO TO 20 
      IF (X.GT.T(I)) GO TO 130 
   10 IF (I.EQ.K) GO TO 140 
      I = I - 1 
      IF (X.EQ.T(I)) GO TO 10 
!                                                                       
! *** DIFFERENCE THE COEFFICIENTS *IDERIV* TIMES                        
!     WORK(I) = AJ(I), WORK(K+I) = DP(I), WORK(K+K+I) = DM(I), I=1.K    
!                                                                       
   20 IMK = I - K 
      DO 30 J=1,K 
        IMKPJ = IMK + J 
        WORK(J) = A(IMKPJ) 
   30 END DO 
      IF (IDERIV.EQ.0) GO TO 60 
      DO 50 J=1,IDERIV 
        KMJ = K - J 
        FKMJ = FLOAT(KMJ) 
        DO 40 JJ=1,KMJ 
          IHI = I + JJ 
          IHMKMJ = IHI - KMJ 
          WORK(JJ) = (WORK(JJ+1)-WORK(JJ))/(T(IHI)-T(IHMKMJ))*FKMJ 
   40   CONTINUE 
   50 END DO 
!                                                                       
! *** COMPUTE VALUE AT *X* IN (T(I),(T(I+1)) OF IDERIV-TH DERIVATIVE,   
!     GIVEN ITS RELEVANT B-SPLINE COEFF. IN AJ(1),...,AJ(K-IDERIV).     
   60 IF (IDERIV.EQ.KM1) GO TO 100 
      IP1 = I + 1 
      KPK = K + K 
      J1 = K + 1 
      J2 = KPK + 1 
      DO 70 J=1,KMIDER 
        IPJ = I + J 
        WORK(J1) = T(IPJ) - X 
        IP1MJ = IP1 - J 
        WORK(J2) = X - T(IP1MJ) 
        J1 = J1 + 1 
        J2 = J2 + 1 
   70 END DO 
      IDERP1 = IDERIV + 1 
      DO 90 J=IDERP1,KM1 
        KMJ = K - J 
        ILO = KMJ 
        DO 80 JJ=1,KMJ 
          WORK(JJ) = (WORK(JJ+1)*WORK(KPK+ILO)+WORK(JJ)                 &
     &              *WORK(K+JJ))/(WORK(KPK+ILO)+WORK(K+JJ))             
          ILO = ILO - 1 
   80   CONTINUE 
   90 END DO 
  100 BVALU = WORK(1) 
      RETURN 
!                                                                       
!                                                                       
  101 CONTINUE 
      CALL XERROR( ' BVALU,  N DOES NOT SATISFY N.GE.K',34,2,1) 
      RETURN 
  102 CONTINUE 
      CALL XERROR( ' BVALU,  K DOES NOT SATISFY K.GE.1',34,2,1) 
      RETURN 
  110 CONTINUE 
      CALL XERROR( ' BVALU,  IDERIV DOES NOT SATISFY 0.LE.IDERIV.LT.K', &
     & 49, 2, 1)                                                        
      RETURN 
  120 CONTINUE 
      CALL XERROR( ' BVALU,  X IS N0T GREATER THAN OR EQUAL TO T(K)',   &
     & 47, 2, 1)                                                        
      RETURN 
  130 CONTINUE 
      CALL XERROR( ' BVALU,  X IS NOT LESS THAN OR EQUAL TO T(N+1)',    &
     & 46, 2, 1)                                                        
      RETURN 
  140 CONTINUE 
      CALL XERROR( ' BVALU,  A LEFT LIMITING VALUE CANN0T BE OBTAINED AT&
     & T(K)',     57, 2, 1)                                             
      RETURN 
      END                                           
      SUBROUTINE INTRV(XT,LXT,X,ILO,ILEFT,MFLAG) 
!***BEGIN PROLOGUE  INTRV                                               
!***DATE WRITTEN   800901   (YYMMDD)                                    
!***REVISION DATE  820801   (YYMMDD)                                    
!***CATEGORY NO.  E3,K6                                                 
!***KEYWORDS  B-SPLINE,DATA FITTING,INTERPOLATION,SPLINE                
!***AUTHOR  AMOS, D. E., (SNLA)                                         
!***PURPOSE  COMPUTES THE LARGEST INTEGER ILEFT IN 1.LE.ILEFT.LE.LXT    
!            SUCH THAT XT(ILEFT).LE.X WHERE XT(*) IS A SUBDIVISION      
!            OF THE X INTERVAL.                                         
!***DESCRIPTION                                                         
!                                                                       
!     WRITTEN BY CARL DE BOOR AND MODIFIED BY D. E. AMOS                
!                                                                       
!     REFERENCE                                                         
!         SIAM J. NUMERICAL ANALYSIS, 14, NO. 3, JUNE 1977, PP. 441-472.
!                                                                       
!     ABSTRACT                                                          
!         INTRV IS THE INTERV ROUTINE OF THE REFERENCE.                 
!                                                                       
!         INTRV COMPUTES THE LARGEST INTEGER ILEFT IN 1 .LE. ILEFT .LE. 
!         LXT SUCH THAT XT(ILEFT) .LE. X WHERE XT(*) IS A SUBDIVISION OF
!         THE X INTERVAL.  PRECISELY,                                   
!                                                                       
!                      X .LT. XT(1)                1         -1         
!         IF  XT(I) .LE. X .LT. XT(I+1)  THEN  ILEFT=I  , MFLAG=0       
!           XT(LXT) .LE. X                         LXT        1,        
!                                                                       
!         THAT IS, WHEN MULTIPLICITIES ARE PRESENT IN THE BREAK POINT   
!         TO THE LEFT OF X, THE LARGEST INDEX IS TAKEN FOR ILEFT.       
!                                                                       
!     DESCRIPTION OF ARGUMENTS                                          
!         INPUT                                                         
!          XT      - XT IS A KNOT OR BREAK POINT VECTOR OF LENGTH LXT   
!          LXT     - LENGTH OF THE XT VECTOR                            
!          X       - ARGUMENT                                           
!          ILO     - AN INITIALIZATION PARAMETER WHICH MUST BE SET      
!                    TO 1 THE FIRST TIME THE SPLINE ARRAY XT IS         
!                    PROCESSED BY INTRV.                                
!                                                                       
!         OUTPUT                                                        
!          ILO     - ILO CONTAINS INFORMATION FOR EFFICIENT PROCESS-    
!                    ING AFTER THE INITIAL CALL, AND ILO MUST NOT BE    
!                    CHANGED BY THE USER.  DISTINCT SPLINES REQUIRE     
!                    DISTINCT ILO PARAMETERS.                           
!          ILEFT   - LARGEST INTEGER SATISFYING XT(ILEFT) .LE. X        
!          MFLAG   - SIGNALS WHEN X LIES OUT OF BOUNDS                  
!                                                                       
!     ERROR CONDITIONS                                                  
!         NONE                                                          
!***REFERENCES  C. DE BOOR, *PACKAGE FOR CALCULATING WITH B-SPLINES*,   
!                 SIAM JOURNAL ON NUMERICAL ANALYSIS, VOLUME 14, NO. 3, 
!                 JUNE 1977, PP. 441-472.                               
!***ROUTINES CALLED  (NONE)                                             
!***END PROLOGUE  INTRV                                                 
!                                                                       
!                                                                       
      INTEGER IHI, ILEFT, ILO, ISTEP, LXT, MFLAG, MIDDLE 
      REAL X, XT 
      DIMENSION XT(LXT) 
!***FIRST EXECUTABLE STATEMENT  INTRV                                   
      IHI = ILO + 1 
      IF (IHI.LT.LXT) GO TO 10 
      IF (X.GE.XT(LXT)) GO TO 110 
      IF (LXT.LE.1) GO TO 90 
      ILO = LXT - 1 
      IHI = LXT 
!                                                                       
   10 IF (X.GE.XT(IHI)) GO TO 40 
      IF (X.GE.XT(ILO)) GO TO 100 
!                                                                       
! *** NOW X .LT. XT(IHI) . FIND LOWER BOUND                             
      ISTEP = 1 
   20 IHI = ILO 
      ILO = IHI - ISTEP 
      IF (ILO.LE.1) GO TO 30 
      IF (X.GE.XT(ILO)) GO TO 70 
      ISTEP = ISTEP*2 
      GO TO 20 
   30 ILO = 1 
      IF (X.LT.XT(1)) GO TO 90 
      GO TO 70 
! *** NOW X .GE. XT(ILO) . FIND UPPER BOUND                             
   40 ISTEP = 1 
   50 ILO = IHI 
      IHI = ILO + ISTEP 
      IF (IHI.GE.LXT) GO TO 60 
      IF (X.LT.XT(IHI)) GO TO 70 
      ISTEP = ISTEP*2 
      GO TO 50 
   60 IF (X.GE.XT(LXT)) GO TO 110 
      IHI = LXT 
!                                                                       
! *** NOW XT(ILO) .LE. X .LT. XT(IHI) . NARROW THE INTERVAL             
   70 MIDDLE = (ILO+IHI)/2 
      IF (MIDDLE.EQ.ILO) GO TO 100 
!     NOTE. IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1         
      IF (X.LT.XT(MIDDLE)) GO TO 80 
      ILO = MIDDLE 
      GO TO 70 
   80 IHI = MIDDLE 
      GO TO 70 
! *** SET OUTPUT AND RETURN                                             
   90 MFLAG = -1 
      ILEFT = 1 
      RETURN 
  100 MFLAG = 0 
      ILEFT = ILO 
      RETURN 
  110 MFLAG = 1 
      ILEFT = LXT 
      RETURN 
      END                                           
      SUBROUTINE XERROR(MESSG,NMESSG,NERR,LEVEL) 
!***BEGIN PROLOGUE  XERROR                                              
!***DATE WRITTEN   790801   (YYMMDD)                                    
!***REVISION DATE  870930   (YYMMDD)                                    
!***CATEGORY NO.  R3C                                                   
!***KEYWORDS  ERROR,XERROR PACKAGE                                      
!***AUTHOR  JONES, R. E., (SNLA)                                        
!***PURPOSE  Processes an error (diagnostic) message.                   
!***DESCRIPTION                                                         
!    From the book "Numerical Methods and Software"                     
!       by  D. Kahaner, C. Moler, S. Nash                               
!           Prentice Hall 1988                                          
!     Abstract                                                          
!        XERROR processes a diagnostic message. It is a stub routine    
!        written for the book above. Actually, XERROR is a sophisticated
!        error handling package with many options, and is described     
!        in the reference below. Our version has the same calling sequen
!        but only prints an error message and either returns (if the    
!        input value of ABS(LEVEL) is less than 2) or stops (if the     
!        input value of ABS(LEVEL) equals 2).                           
!                                                                       
!     Description of Parameters                                         
!      --Input--                                                        
!        MESSG - the Hollerith message to be processed.                 
!        NMESSG- the actual number of characters in MESSG.              
!                (this is ignored in this stub routine)                 
!        NERR  - the error number associated with this message.         
!                NERR must not be zero.                                 
!                (this is ignored in this stub routine)                 
!        LEVEL - error category.                                        
!                =2 means this is an unconditionally fatal error.       
!                =1 means this is a recoverable error.  (I.e., it is    
!                   non-fatal if XSETF has been appropriately called.)  
!                =0 means this is a warning message only.               
!                =-1 means this is a warning message which is to be     
!                   printed at most once, regardless of how many        
!                   times this call is executed.                        
!                 (in this stub routine                                 
!                       LEVEL=2 causes a message to be printed and then 
!                                         stop.                         
!                       LEVEL<2 causes a message to be printed and then 
!                                         return.                       
!                                                                       
!     Examples                                                          
!        CALL XERROR('SMOOTH -- NUM WAS ZERO.',23,1,2)                  
!        CALL XERROR('INTEG  -- LESS THAN FULL ACCURACY ACHIEVED.',     
!                    43,2,1)                                            
!        CALL XERROR('ROOTER -- ACTUAL ZERO OF F FOUND BEFORE INTERVAL F
!    1ULLY COLLAPSED.',65,3,0)                                          
!        CALL XERROR('EXP    -- UNDERFLOWS BEING SET TO ZERO.',39,1,-1) 
!                                                                       
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-    
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,  
!                 1982.                                                 
!***ROUTINES CALLED  XERRWV                                             
!***END PROLOGUE  XERROR                                                
      CHARACTER*(*) MESSG 
!***FIRST EXECUTABLE STATEMENT  XERROR                                  
      CALL XERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.) 
      RETURN 
      END                                           
      SUBROUTINE XERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2) 
!***BEGIN PROLOGUE  XERRWV                                              
!***DATE WRITTEN   800319   (YYMMDD)                                    
!***REVISION DATE  870930   (YYMMDD)                                    
!***CATEGORY NO.  R3C                                                   
!***KEYWORDS  ERROR,XERROR PACKAGE                                      
!***AUTHOR  JONES, R. E., (SNLA)                                        
!***PURPOSE  Processes error message allowing 2 integer and two real    
!            values to be included in the message.                      
!***DESCRIPTION                                                         
!    From the book "Numerical Methods and Software"                     
!       by  D. Kahaner, C. Moler, S. Nash                               
!           Prentice Hall 1988                                          
!     Abstract                                                          
!        XERRWV prints a diagnostic error message.                      
!        In addition, up to two integer values and two real             
!        values may be printed along with the message.                  
!        A stub routine for the book above. The actual XERRWV is describ
!        in the reference below and contains many other options.        
!                                                                       
!     Description of Parameters                                         
!      --Input--                                                        
!        MESSG - the Hollerith message to be processed.                 
!        NMESSG- the actual number of characters in MESSG.              
!                (ignored in this stub)                                 
!        NERR  - the error number associated with this message.         
!                NERR must not be zero.                                 
!                (ignored in this stub)                                 
!        LEVEL - error category.                                        
!                =2 means this is an unconditionally fatal error.       
!                =1 means this is a recoverable error.  (I.e., it is    
!                   non-fatal if XSETF has been appropriately called.)  
!                =0 means this is a warning message only.               
!                =-1 means this is a warning message which is to be     
!                   printed at most once, regardless of how many        
!                   times this call is executed.                        
!                  (in this stub LEVEL=2 causes an error message to be  
!                                          printed followed by a stop,  
!                                LEVEL<2 causes an error message to be  
!                                          printed followed by a return.
!        NI    - number of integer values to be printed. (0 to 2)       
!        I1    - first integer value.                                   
!        I2    - second integer value.                                  
!        NR    - number of real values to be printed. (0 to 2)          
!        R1    - first real value.                                      
!        R2    - second real value.                                     
!                                                                       
!     Examples                                                          
!        CALL XERRWV('SMOOTH -- NUM (=I1) WAS ZERO.',29,1,2,            
!    1   1,NUM,0,0,0.,0.)                                               
!        CALL XERRWV('QUADXY -- REQUESTED ERROR (R1) LESS THAN MINIMUM (
!    1R2).,54,77,1,0,0,0,2,ERRREQ,ERRMIN)                               
!                                                                       
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-    
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,  
!                 1982.                                                 
!***ROUTINES CALLED  (NONE)                                             
!***END PROLOGUE  XERRWV                                                
      CHARACTER*(*) MESSG 
!***FIRST EXECUTABLE STATEMENT  XERRWV                                  
      WRITE(*,*) MESSG 
      IF(NI.EQ.2)THEN 
        WRITE(*,*) I1,I2 
      ELSEIF(NI.EQ.1) THEN 
        WRITE(*,*) I1 
      ENDIF 
      IF(NR.EQ.2) THEN 
        WRITE(*,*) R1,R2 
      ELSEIF(NR.EQ.1) THEN 
        WRITE(*,*) R1 
      ENDIF 
      IF(ABS(LEVEL).LT.2)RETURN 
      STOP 
      END                                           
