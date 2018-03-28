      SUBROUTINE SMOOTH(IERR,TRACE,NORD,N,X,Y,S,SMIN,Z,BADPT,MAXBK,     &
     &                  NBKPTS,BKPT,COEFF,XHALF,LW,W)                   
!                                                                       
! --- This subroutine is a driver for the data-smoothing routine DOFIT. 
!     If the data set to be smoothed is very large, the driver first    
!     smooths a subset of the data and then uses the breakpoints from   
!     this initial fit as a starting place for smoothing the entire     
!     data set.                                                         
!                                                                       
! --- Input Variables                                                   
!     TRACE:         Flag indicating whether to print a trace of the    
!                    routine's progress --                              
!                    TRACE=0 for no trace                               
!                    TRACE=1 for trace of number of data points and     
!                            number of breakpoints in use               
!     NORD:          Order of piecewise polynomials for spline fit      
!                    (NORD=4 for cubic spline)                          
!     N:             Number of data points                              
!     X(N):          Array containing x coordinates of data             
!     Y(N):          Array containing y coordinates of data             
!     S(N):          Array containing y-coordinate error estimates      
!     SMIN:          Minimum allowed error tolerance                    
!     MAXBK:         Maximum allowed number of breakpoints in smoothed  
!                    curve                                              
!     NBKPTS:        Number of initial breakpoints (should be zero if   
!                    data are being smoothed for the first time)        
!     BKPT(NBKPTS):  Array of initial breakpoints (used only if         
!                    NBKPTS>0)                                          
!     LW:            Size of work array W                               
!                                                                       
! --- Output Variables                                                  
!     IERR:          Error indicator --                                 
!                    IERR=0 if no errors were encountered               
!                    IERR=1 if X array is not sorted                    
!                    IERR=2 if an error occurred in subroutine EFC      
!                    IERR=3 if an error occurred in subroutine INTERV   
!                    IERR=4 if the maximum number of breakpoints was    
!                           exceeded                                    
!                    IERR=5 if an error occurred in subroutine INSERT   
!                    IERR=6 if subroutine DOFIT was unable to add new   
!                           breakpoints                                 
!     SMIN:          Minimum allowed error tolerance (changed if input  
!                    was less than or equal to zero)                    
!     Z(N):          Array containing smoothed y coordinates            
!     BADPT(N):      Work array                                         
!     NBKPTS:        Number of breakpoints in smoothed curve            
!     BKPT(NBKPTS):  Array of breakpoints                               
!     COEFF(NBKPTS): Array of b-spline coefficients                     
!     XHALF(NBKPTS): Work array                                         
!     W(LW):         Work array                                         
!                                                                       
! --- Declaration of Calling Sequence Variables                         
      INTEGER IERR,TRACE,NORD,N,MAXBK,NBKPTS,LW 
      REAL X(*),Y(*),S(*),SMIN,Z(*),BADPT(*),BKPT(*),COEFF(*),W(*),     &
     &     XHALF(*)                                                     
!                                                                       
! --- Declaration of Internal Variables                                 
      INTEGER I,NINIT,NTEMP,STEP,INDEX 
      PARAMETER(NINIT=200) 
      REAL XTEMP(NINIT),YTEMP(NINIT),STEMP(NINIT) 
!                                                                       
! --- If data set is large, smooth a subset first.                      
      IF(N.GT.NINIT)THEN 
         STEP=(N-1)/(NINIT-1) 
         IF(MOD(N-1,NINIT-1).NE.0)STEP=STEP+1 
         NTEMP=(N-1)/STEP+2 
         IF(MOD(N-1,STEP).EQ.0)NTEMP=NTEMP-1 
         IF(TRACE.EQ.1)WRITE(*,50)NTEMP 
   50    FORMAT('    Using ',I5,' data points ... ') 
         DO 100 I=1,NTEMP-1 
            INDEX=(I-1)*STEP+1 
            XTEMP(I)=X(INDEX) 
            YTEMP(I)=Y(INDEX) 
            STEMP(I)=S(INDEX) 
  100    CONTINUE 
         XTEMP(NTEMP)=X(N) 
         YTEMP(NTEMP)=Y(N) 
         STEMP(NTEMP)=S(N) 
         CALL DOFIT(IERR,TRACE,NTEMP,XTEMP,YTEMP,STEMP,SMIN,Z,BADPT,    &
     &              NORD,MAXBK,NBKPTS,BKPT,COEFF,XHALF,LW,W)            
         IF((IERR.EQ.1).OR.(IERR.EQ.4))RETURN 
      ENDIF 
!                                                                       
! --- Smooth entire data set.                                           
      IF(TRACE.EQ.1)WRITE(*,200)N 
  200 FORMAT('    Using ',I5,' data points ... ') 
      CALL DOFIT(IERR,TRACE,N,X,Y,S,SMIN,Z,BADPT,NORD,MAXBK,NBKPTS,BKPT,&
     &           COEFF,XHALF,LW,W)                                      
!                                                                       
      RETURN 
      END                                           
!                                                                       
!     ******************************************************************
!                                                                       
      SUBROUTINE DOFIT(IERR,TRACE,N,X,Y,S,SMIN,Z,BADPT,NORD,MAXBK,      &
     &                 NBKPTS,BKPT,COEFF,XHALF,LW,W)                    
!                                                                       
! --- This subroutine smooths a given data set, using the CMLIB b-spline
!     routines EFC and BVALU.  The breakpoints for the b-splines are    
!     determined iteratively, with more breakpoints set in regions where
!     the data vary the most.                                           
!                                                                       
! --- Input Variables                                                   
!     TRACE:         Flag indicating whether to print a trace of the    
!                    routine's progress --                              
!                    TRACE=0 for no trace                               
!                    TRACE=1 for trace of number of data points and     
!                            number of breakpoints in use               
!     N:             Number of data points                              
!     X(N):          Array containing x coordinates of data             
!     Y(N):          Array containing y coordinates of data             
!     S(N):          Array containing y-coordinate error estimates      
!     SMIN:          Minimum allowed error tolerance                    
!     NORD:          Order of piecewise polynomials for spline fit      
!                    (NORD=4 for cubic spline)                          
!     MAXBK:         Maximum allowed number of breakpoints in smoothed  
!                    curve                                              
!     NBKPTS:        Number of initial breakpoints (should be zero if   
!                    data is being smoothed for the first time)         
!     BKPT(NBKPTS):  Array of initial breakpoints (used only if         
!                    NBKPTS>0)                                          
!     LW:            Size of work array W                               
!                                                                       
! --- Output Variables                                                  
!     IERR:          Error indicator --                                 
!                    IERR=0 if no error was encountered                 
!                    IERR=1 if X array is not sorted                    
!                    IERR=2 if an error occurred in subroutine EFC      
!                    IERR=3 if an error occurred in subroutine INTERV   
!                    IERR=4 if the maximum number of breakpoints was    
!                           exceeded                                    
!                    IERR=5 if an error occurred in subroutine INSERT   
!                    IERR=6 if DOFIT was unable to add new breakpoints  
!     SMIN:          Minimum allowed error tolerance (changed if input  
!                    was less than or equal to zero)                    
!     Z(N):          Array containing smoothed y coordinates            
!     BADPT(N):      Work array                                         
!     NBKPTS:        Number of breakpoints in smoothed curve            
!     BKPT(NBKPTS):  Array of breakpoints                               
!     COEFF(NBKPTS): Array of b-spline coefficients                     
!     XHALF(NBKPTS): Work array                                         
!     W(LW):         Work array                                         
!                                                                       
! --- Declaration of Calling Sequence Variables                         
      INTEGER IERR,TRACE,N,NORD,MAXBK,NBKPTS,LW 
      REAL X(*),Y(*),S(*),SMIN,Z(*),BADPT(*),BKPT(*),COEFF(*),W(*),     &
     &     XHALF(*)                                                     
!                                                                       
! --- Declaration of Internal Variables                                 
      LOGICAL ERROR 
      INTEGER I,ITER,MDEOUT,INBV,IDER,BEGIN,END,NPTS,NBAD,OLDBEG,OLDEND,&
     &        NEWBKS,OLDBKS                                             
      REAL BVALU,W2(20) 
!                                                                       
! --- Initialize parameters.                                            
      IERR=0 
      ITER=0 
!                                                                       
! --- Check if x coordinates are in order.                              
      DO 50 I=2,N 
         IF(X(I).LT.X(I-1))THEN 
            IERR=1 
            RETURN 
         ENDIF 
   50 END DO 
!                                                                       
! --- Check if initial breakpoints are supplied.                        
      IF(NBKPTS.GT.0)GOTO 300 
!                                                                       
! --- If not, set initial breakpoints at first and last data points.    
      NBKPTS=2*NORD 
      DO 100 I=1,NORD 
         BKPT(I)=X(1) 
  100 END DO 
      DO 200 I=NORD+1,NBKPTS 
         BKPT(I)=X(N) 
  200 END DO 
!                                                                       
! --- Fit curve through data points.                                    
  300 CONTINUE 
      ITER=ITER+1 
      IF(TRACE.EQ.1)WRITE(*,350)NBKPTS 
  350 FORMAT('       Trying ',I3,' knots ... ') 
      CALL EFC(N,X,Y,S,NORD,NBKPTS,BKPT,1,MDEOUT,COEFF,LW,W) 
      IF(MDEOUT.NE.1)THEN 
         IERR=2 
!        PRINT *,'ERROR in subroutine DOFIT: Error returned from ',     
!    +           'call to EFC.'                                         
!         WRITE(*,380)MDEOUT                                            
! 380    FORMAT('       WARNING: MDEOUT = ',I3,' on return from EFC.')  
         RETURN 
      ENDIF 
!                                                                       
! --- Calculate smoothed data points.                                   
      INBV=1 
      IDER=0 
      DO 400 I=1,N 
         Z(I)=BVALU(BKPT,COEFF,NBKPTS-NORD,NORD,IDER,X(I),INBV,W2) 
  400 END DO 
!                                                                       
! --- Check if smoothed points are within requested accuracy.           
      CALL CHKFIT(N,X,Y,Z,S,SMIN,NBAD,BADPT) 
      IF(NBAD.EQ.0)RETURN 
!                                                                       
! --- If smoothed points are not accurate enough, add breakpoints.      
      OLDBEG=0 
      OLDEND=0 
      OLDBKS=NBKPTS 
      NEWBKS=0 
      DO 500 I=1,NBAD 
         IF(OLDEND.NE.0)THEN 
            IF(BADPT(I).LT.X(OLDEND))GOTO 500 
         ENDIF 
         CALL INTERV(BADPT(I),N,X,NORD,NBKPTS,BKPT,BEGIN,END) 
         IF(BEGIN.EQ.OLDBEG.AND.END.EQ.OLDEND)GOTO 500 
         NPTS=END-BEGIN+1 
         IF(NPTS.LE.0)THEN 
            IERR=3 
!           PRINT *,'ERROR in subroutine DOFIT: Number of points is ',  
!    +              'zero or negative.'                                 
            RETURN 
         ENDIF 
         NEWBKS=NEWBKS+1 
         IF(NEWBKS.GT.MAXBK)THEN 
            IERR=4 
!           PRINT *,'ERROR in subroutine DOFIT: Maximum number of ',    
!    +              'breakpoints has been exceeded --'                  
!           PRINT *,'Increase MAXBK or use a bigger error tolerance.'   
            RETURN 
         ENDIF 
         CALL HALFPT(NPTS,X(BEGIN),Y(BEGIN),XHALF(NEWBKS)) 
         OLDBEG=BEGIN 
         OLDEND=END 
  500 END DO 
      DO 600 I=1,NEWBKS 
         CALL INSERT(ERROR,XHALF(I),MAXBK,NBKPTS,BKPT) 
         IF(ERROR)THEN 
            IERR=5 
            RETURN 
         ENDIF 
  600 END DO 
      CALL DELETE(NORD,NBKPTS,BKPT) 
      IF(NBKPTS.LE.OLDBKS)THEN 
         IERR=6 
!        PRINT *,'ERROR in subroutine DOFIT: After deleting ',          
!    +           'extraneous multiple breakpoints, the total number ',  
!    +           'of breakpoints decreased or stayed the same.'         
         RETURN 
      ENDIF 
      GOTO 300 
!                                                                       
      END                                           
!                                                                       
!     ******************************************************************
!                                                                       
      SUBROUTINE CHKFIT(N,X,Y,Z,S,SMIN,NBAD,BADPT) 
!                                                                       
! --- This subroutine compares each smoothed data point with the        
!     corresponding original data point to check whether the smoothed   
!     curve is within the requested tolerance.                          
!                                                                       
! --- Input Variables                                                   
!     N:           Number of data points                                
!     X(N):        Array containing x coordinates of data points        
!     Y(N):        Array containing y coordinates of data points        
!     Z(N):        Array containing y coordinates of smoothed data      
!     S(N):        Array containing y-coordinate error tolerances       
!     SMIN:        Minimum allowed error tolerance                      
!                                                                       
! --- Output Variables                                                  
!     SMIN:        Minimum allowed error tolerance (changed if input was
!                  less than or equal to zero)                          
!     NBAD:        Number of smoothed points that do not satisfy error  
!                  test                                                 
!     BADPT(NBAD): Array containing x coordinates of smoothed points    
!                  that do not satisfy error test                       
!                                                                       
! --- Declaration of Calling Sequence Variables                         
      INTEGER N,NBAD 
      REAL X(*),Y(*),Z(*),S(*),SMIN,BADPT(*) 
!                                                                       
! --- Declaration of Internal Variables                                 
      INTEGER I 
      REAL TOL,YMIN,YMAX 
!                                                                       
! --- Reset minimum error tolerance, if not greater than zero.          
      IF(SMIN.LE.0.)THEN 
         YMIN=Y(1) 
         YMAX=Y(1) 
         DO 50 I=2,N 
            IF(Y(I).LT.YMIN)YMIN=Y(I) 
            IF(Y(I).GT.YMAX)YMAX=Y(I) 
   50    CONTINUE 
         SMIN=.001*(YMAX-YMIN) 
      ENDIF 
!                                                                       
! --- Compare smoothed points with original data points.                
      NBAD=0 
      DO 100 I=1,N 
         TOL=MAX(ABS(S(I)),SMIN) 
         IF(ABS(Z(I)-Y(I)).GT.TOL)THEN 
            NBAD=NBAD+1 
            BADPT(NBAD)=X(I) 
         ENDIF 
  100 END DO 
!                                                                       
      RETURN 
      END                                           
!                                                                       
!     ******************************************************************
!                                                                       
      SUBROUTINE DELETE(NORD,NBKPTS,BKPT) 
!                                                                       
! --- This subroutine checks for coincident breakpoints and deletes any 
!     over a multiplicity of NORD.  This is done because subroutine EFC 
!     fails with an error "MDEOUT=2" whenever it receives a breakpoint  
!     whose multiplicity is more than NORD.                             
!                                                                       
! --- Input Variables                                                   
!     NORD:         Order of piecewise polynomials for spline fit       
!                   (NORD=4 for cubic spline)                           
!     NBKPTS:       Number of breakpoints (before deletion)             
!     BKPT(NBKPTS): Array containing breakpoints (before deletion)      
!                                                                       
! --- Output Variables                                                  
!     NBKPTS:       Number of breakpoints (after deletion)              
!     BKPT(NBKPTS): Array containing breakpoints (after deletion)       
!                                                                       
! --- Declaration of Calling Sequence Variables                         
      INTEGER NORD,NBKPTS 
      REAL BKPT(*) 
!                                                                       
! --- Declaration of Internal Variables                                 
      INTEGER I,J,INIT 
!                                                                       
! --- Initialize parameters.                                            
      INIT=1 
!                                                                       
! --- Check that number of breakpoints is not degenerate.               
  100 CONTINUE 
      IF(NBKPTS.GT.NORD)THEN 
!                                                                       
! ---    Loop to compare pairs of breakpoints.                          
         DO 300 I=INIT,NBKPTS-NORD 
!                                                                       
! ---       Check if there are more than NORD coincident breakpoints.   
            IF(BKPT(I).EQ.BKPT(I+NORD))THEN 
!                                                                       
! ---          Remove (I+NORD)th breakpoint.                            
               DO 200 J=I+NORD,NBKPTS-1 
                  BKPT(J)=BKPT(J+1) 
  200          CONTINUE 
               NBKPTS=NBKPTS-1 
               INIT=I 
               GOTO 100 
!                                                                       
            ENDIF 
!                                                                       
  300    CONTINUE 
!                                                                       
      ENDIF 
!                                                                       
      RETURN 
      END                                           
!                                                                       
!     ******************************************************************
!                                                                       
      SUBROUTINE HALFPT(N,X,Y,XHALF) 
!                                                                       
! --- This subroutine calculates the (scaled) length of the piecewise   
!     linear interpolant to the given data and returns the x coordinate 
!     of a point approximately half way along the curve.                
!                                                                       
! --- Input Variables                                                   
!     N:     Number of data points                                      
!     X(N):  Array containing x coordinates of data                     
!     Y(N):  Array containing y coordinates of data                     
!                                                                       
! --- Output Variable                                                   
!     XHALF: X coordinate of half-way point                             
!                                                                       
! --- Declaration of Calling Sequence Variables                         
      INTEGER N 
      REAL X(*),Y(*),XHALF 
!                                                                       
! --- Declaration of Internal Variables                                 
      INTEGER I 
      REAL XMIN,XMAX,YMIN,YMAX,XSCALE,YSCALE,LENGTH,XTEMP,YTEMP,DIST 
!                                                                       
! --- Find range of data and set scale factors.                         
      XMIN=X(1) 
      XMAX=X(1) 
      YMIN=Y(1) 
      YMAX=Y(1) 
      DO 100 I=1,N 
         IF(X(I).LT.XMIN)XMIN=X(I) 
         IF(X(I).GT.XMAX)XMAX=X(I) 
         IF(Y(I).LT.YMIN)YMIN=Y(I) 
         IF(Y(I).GT.YMAX)YMAX=Y(I) 
  100 END DO 
      XSCALE=XMAX-XMIN 
      IF(XSCALE.EQ.0.)XSCALE=1. 
      YSCALE=YMAX-YMIN 
      IF(YSCALE.EQ.0.)YSCALE=1. 
!                                                                       
! --- Calculate total length of data curve.                             
      LENGTH=0. 
      DO 200 I=2,N 
         XTEMP=(X(I)-X(I-1))/XSCALE 
         YTEMP=(Y(I)-Y(I-1))/YSCALE 
         LENGTH=LENGTH+SQRT(XTEMP*XTEMP+YTEMP*YTEMP) 
  200 END DO 
!                                                                       
! --- Find x value of half-way point.                                   
      DIST=0. 
      XHALF=X(N) 
      DO 300 I=2,N 
         XTEMP=(X(I)-X(I-1))/XSCALE 
         YTEMP=(Y(I)-Y(I-1))/YSCALE 
         DIST=DIST+SQRT(XTEMP*XTEMP+YTEMP*YTEMP) 
         IF(2*DIST.GE.LENGTH)THEN 
            XHALF=(X(I)+X(I-1))/2. 
            GOTO 400 
         ENDIF 
  300 END DO 
!                                                                       
  400 CONTINUE 
      RETURN 
      END                                           
!                                                                       
!     ******************************************************************
!                                                                       
      SUBROUTINE INSERT(ERROR,NEWBK,MAXBK,NBKPTS,BKPT) 
!                                                                       
! --- This subroutine inserts a new knot in the breakpoint array at a   
!     specified location.                                               
!                                                                       
! --- Input Variables                                                   
!     NEWBK:        X position of new knot to be inserted               
!     MAXBK:        Maximum allowed number of breakpoints               
!     NBKPTS:       Number of breakpoints (before insertion)            
!     BKPT(NBKPTS): Array containing breakpoints (before insertion)     
!                                                                       
! --- Output Variables                                                  
!     ERROR:        Error flag --                                       
!                   ERROR=.TRUE. if an error was encountered            
!                   ERROR=.FALSE. otherwise                             
!     NBKPTS:       Number of breakpoints (after insertion)             
!     BKPT(NBKPTS): Array containing breakpoints (after insertion)      
!                                                                       
! --- Declaration of Calling Sequence Variables                         
      LOGICAL ERROR 
      INTEGER MAXBK,NBKPTS 
      REAL NEWBK,BKPT(*) 
!                                                                       
! --- Declaration of Internal Variables                                 
      INTEGER I,NEWPOS 
!                                                                       
! --- Initialize error flag.                                            
      ERROR=.FALSE. 
!                                                                       
! --- Find position in breakpoint array for new breakpoint.             
      DO 100 I=1,NBKPTS 
         IF(NEWBK.LT.BKPT(I))THEN 
            NEWPOS=I 
            GOTO 200 
         ENDIF 
  100 END DO 
      NEWPOS=NBKPTS+1 
!                                                                       
! --- Insert new breakpoint in breakpoint array.                        
  200 CONTINUE 
      NBKPTS=NBKPTS+1 
      IF(NBKPTS.GT.MAXBK)THEN 
         ERROR=.TRUE. 
!        PRINT *,'ERROR in subroutine INSERT: Maximum number of ',      
!    +           'breakpoints has been exceeded --'                     
!        PRINT *,'Increase MAXBK or use a bigger error tolerance.'      
         RETURN 
      ENDIF 
      IF(NEWPOS.LT.NBKPTS)THEN 
         DO 300 I=NBKPTS,NEWPOS+1,-1 
            BKPT(I)=BKPT(I-1) 
  300    CONTINUE 
      ENDIF 
      BKPT(NEWPOS)=NEWBK 
!                                                                       
      RETURN 
      END                                           
!                                                                       
!     ******************************************************************
!                                                                       
      SUBROUTINE INTERV(XPT,N,X,NORD,NBKPTS,BKPT,BEGIN,END) 
!                                                                       
! --- This subroutine finds the smallest breakpoint interval containing 
!     the given x value and returns the indices of the data points      
!     nearest the endpoints of this interval.                           
!                                                                       
! --- INPUT VARIABLES                                                   
!     XPT:          X value to be located                               
!     N:            Number of data points                               
!     X(N):         Array containing x coordinates of data              
!     NORD:         Order of piecewise polynomials for spline fit       
!                   (NORD=4 for cubic spline)                           
!     NBKPTS:       Number of breakpoints                               
!     BKPT(NBKPTS): Array containing breakpoints                        
!                                                                       
! --- Output Variables                                                  
!     BEGIN:        Index of data point near left breakpoint            
!                   ( X(BEGIN) <= Left Breakpoint <= XPT )              
!     END:          Index of data point near right breakpoint           
!                   ( XPT <= Right Breakpoint <= X(END) )               
!                                                                       
! --- Declaration of Calling Sequence Variables                         
      INTEGER N,NORD,NBKPTS,BEGIN,END 
      REAL XPT,X(*),BKPT(*) 
!                                                                       
! --- Declaration of Internal Variables                                 
      INTEGER I 
      REAL BKPT1,BKPT2 
!                                                                       
! --- Find breakpoint interval containing the given point.              
      DO 100 I=NORD+1,NBKPTS 
         IF(XPT.LE.BKPT(I))THEN 
            BKPT1=BKPT(I-1) 
            BKPT2=BKPT(I) 
            GOTO 200 
         ENDIF 
  100 END DO 
!                                                                       
! --- Find data points nearest the ends of the breakpoint interval.     
  200 CONTINUE 
      DO 300 I=2,N 
         IF(X(I).GT.BKPT1)THEN 
            BEGIN=I-1 
            GOTO 400 
         ENDIF 
  300 END DO 
  400 CONTINUE 
      DO 500 I=BEGIN,N 
         IF(X(I).GT.BKPT2)THEN 
            END=I 
            GOTO 600 
         ENDIF 
  500 END DO 
      END=N 
!                                                                       
  600 CONTINUE 
      RETURN 
      END                                           
