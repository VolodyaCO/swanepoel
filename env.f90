      SUBROUTINE ENV(ERROR,TOP,N,X,Y,DERIV2,NENV,XENV,YENV,NTAN,XTAN,   &
     &               YTAN,D,WK)                                         
!                                                                       
! --- This subroutine calculates the top or bottom envelope of a given  
!     set of data points.                                               
!                                                                       
! --- Written by M. McClain, January 1989.                              
!     Latest revision, October 1992.                                    
!                                                                       
! --- Description of Input Variables                                    
!     TOP:        Flag indicating whether top or bottom envelope is to  
!                 be computed                                           
!     N:          Number of data points                                 
!     X(N):       Array containing x coordinates of data points         
!     Y(N):       Array containing y coordinates of data points         
!     DERIV2(N):  Array containing second derivatives of data points    
!     NTAN:       Number of tangent points where envelope touches data  
!                 curve                                                 
!     XTAN(NTAN): Array containing x coordinates of initial estimate    
!                 of tangent points                                     
!     YTAN(NTAN): Array containing y coordinates of initial estimate    
!                 of tangent points                                     
!                                                                       
! --- Description of Output Variables                                   
!     ERROR:      Error flag                                            
!     NENV:       Number of points in envelope curve                    
!     XENV(NENV): Array containing x coordinates of envelope            
!     YENV(NENV): Array containing y coordinates of envelope            
!     XTAN(NTAN): Array containing x coordinates of tangent points      
!     YTAN(NTAN): Array containing y coordinates of tangent points      
!     D(N):       Work array                                            
!     WK(2*NTAN): Work array                                            
!                                                                       
! --- Declaration of Calling Sequence Variables                         
      LOGICAL TOP,ERROR 
      INTEGER N,NTAN,NENV,IERR,LAST 
      REAL X(*),Y(*),DERIV2(*),XENV(*),YENV(*),XTAN(*),YTAN(*),D(*),    &
     &     WK(*)                                                        
!                                                                       
! --- Declaration of Internal Variables                                 
      LOGICAL SKIP 
      INTEGER I,J,COUNT,INTERV,BEGIN,END,NWK,IC(2) 
      REAL SWITCH,VC(2) 
!                                                                       
! --- If bottom envelope is to be computed, negate y values, second     
!     derivatives, and tangent values.                                  
      IF(.NOT.TOP)THEN 
         DO 100 I=1,N 
            Y(I)=-Y(I) 
            DERIV2(I)=-DERIV2(I) 
  100    CONTINUE 
         DO 150 I=1,NTAN 
            YTAN(I)=-YTAN(I) 
  150    CONTINUE 
      ENDIF 
!                                                                       
! --- Initialize error flag and iteration counter.                      
      ERROR=.FALSE. 
      COUNT=0 
!                                                                       
  200 CONTINUE 
! --- Interpolate a monotone function through the estimated tangent     
!     points.                                                           
      COUNT=COUNT+1 
      IC(1)=0 
      IC(2)=0 
      SWITCH=-1. 
      NWK=2*(NTAN-1) 
      CALL PCHIC(IC,VC,SWITCH,NTAN,XTAN,YTAN,D,1,WK,NWK,IERR) 
      IF(IERR.LT.0)THEN 
         ERROR=.TRUE. 
         PRINT *,'ERROR in subroutine ENV -- IERR = ',IERR,             &
     &           ' on return from PCHIC.'                               
         RETURN 
      ENDIF 
      SKIP=.FALSE. 
      CALL PCHFE(NTAN,XTAN,YTAN,D,1,SKIP,N,X,YENV,IERR) 
      IF(IERR.LT.0)THEN 
         ERROR=.TRUE. 
         PRINT *,'ERROR in subroutine ENV -- IERR = ',IERR,             &
     &           ' on return from PCHFE.'                               
         RETURN 
      ENDIF 
!                                                                       
! --- Check for data points that lie above the envelope.  (Also check   
!     that the data curve is concave at these points.)                  
      INTERV=0 
      LAST=0 
      DO 500 I=2,N-1 
         IF(I.LE.LAST)GOTO 500 
         IF(Y(I).GT.YENV(I).AND.DERIV2(I).LE.0.)THEN 
            DO 300 J=I+1,N 
               IF(Y(J).LE.YENV(J).OR.DERIV2(J).GT.0.)THEN 
                  LAST=J 
                  GOTO 400 
               ENDIF 
  300       CONTINUE 
            LAST=N 
  400       CONTINUE 
            DO 450 J=1,NTAN 
               IF(XTAN(J).GE.X(I-1).AND.XTAN(J).LE.X(LAST))THEN 
                  XTAN(J)=X((I+LAST-1)/2) 
                  YTAN(J)=Y((I+LAST-1)/2) 
               ENDIF 
  450       CONTINUE 
            INTERV=INTERV+1 
         ENDIF 
  500 END DO 
!                                                                       
! --- Stop if no data points are above the envelope, otherwise try again
!     with new tangent point estimates.                                 
      IF(INTERV.GT.0.AND.COUNT.LT.20)GOTO 200 
!     WRITE(*,520)COUNT                                                 
! 520 FORMAT(4X,I2,' iterations.')                                      
!                                                                       
! --- If bottom envelope has been computed, set y values and second     
!     derivatives back to their original values and negate envelope     
!     and tangent values.                                               
      IF(.NOT.TOP)THEN 
         DO 550 I=1,N 
            Y(I)=-Y(I) 
            DERIV2(I)=-DERIV2(I) 
            YENV(I)=-YENV(I) 
  550    CONTINUE 
         DO 580 I=1,NTAN 
            YTAN(I)=-YTAN(I) 
  580    CONTINUE 
      ENDIF 
!                                                                       
! --- Eliminate extrapolation points.                                   
      DO 600 I=1,N 
         IF(X(I).GE.XTAN(1))THEN 
            BEGIN=I 
            GOTO 700 
         ENDIF 
  600 END DO 
  700 CONTINUE 
      DO 800 I=BEGIN,N 
         IF(X(I).GT.XTAN(NTAN))THEN 
            END=I-1 
            GOTO 900 
         ENDIF 
  800 END DO 
      END=NTAN 
  900 CONTINUE 
      NENV=END-BEGIN+1 
      DO 1000 I=1,NENV 
         XENV(I)=X(BEGIN+I-1) 
         YENV(I)=YENV(BEGIN+I-1) 
 1000 END DO 
!                                                                       
      RETURN 
      END                                           
