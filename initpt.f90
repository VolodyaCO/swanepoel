      SUBROUTINE INITPT(ERROR,TOP,N,Y,DERIV2,MAXTAN,NTAN,INDEX) 
!                                                                       
! --- THIS SUBROUTINE CALCULATES INITIAL ESTIMATES FOR POINTS ON THE    
!     DATA CURVE WHICH SHOULD BE TANGENT TO THE TOP OR BOTTOM ENVELOPE  
!     CURVE.  THE METHOD USED IS TO FIND INTERVALS WHERE THE DATA CURVE 
!     IS CONCAVE (FOR THE TOP ENVELOPE) OR CONVEX (FOR THE BOTTOM       
!     ENVELOPE).  THE TANGENT ESTIMATES ARE SET AT THE MIDPOINTS OF     
!     THESE INTERVALS.  CONCAVE (CONVEX) INTERVALS CONTAINING THE FIRST 
!     OR LAST DATA POINT ARE USED ONLY IF THEY ALSO CONTAIN A LOCAL     
!     MAXIMUM (MINIMUM), SINCE OTHERWISE THE DATA FOR THESE INTERVALS   
!     MAY BE INCOMPLETE.                                                
!                                                                       
! --- INPUT VARIABLES                                                   
!     TOP:         FLAG INDICATING WHETHER TANGENTS FOR THE TOP OR THE  
!                  BOTTOM ENVELOPE ARE TO BE COMPUTED                   
!     N:           NUMBER OF DATA POINTS                                
!     Y(N):        ARRAY CONTAINING Y COORDINATES OF DATA POINTS        
!     DERIV2(N):   ARRAY CONTAINING SECOND DERIVATIVES OF DATA POINTS   
!     MAXTAN:      MAXIMUM NUMBER OF TANGENT POINTS ALLOWED             
!                                                                       
! --- OUTPUT VARIABLES                                                  
!     ERROR:       ERROR FLAG                                           
!     NTAN:        NUMBER OF ESTIMATED TANGENT POINTS ON ENVELOPE       
!     INDEX(NTAN): ARRAY CONTAINING INDEXES OF ESTIMATED TANGENT POINTS 
!                  ((X(INDEX(I)),Y(INDEX(I))) IS THE ITH ESTIMATED      
!                  TANGENT POINT.)                                      
!                                                                       
! --- DECLARATION OF CALLING SEQUENCE VARIABLES                         
      LOGICAL ERROR,TOP 
      INTEGER N,MAXTAN,NTAN,INDEX(*) 
      REAL Y(*),DERIV2(*) 
!                                                                       
! --- DECLARATION OF INTERNAL VARIABLES                                 
      INTEGER I,J,K,LAST 
!                                                                       
! --- IF BOTTOM ENVELOPE IS TO BE COMPUTED, NEGATE Y VALUES AND SECOND  
! --- DERIVATIVES                                                       
      IF(.NOT.TOP)THEN 
         DO 20 I=1,N 
            Y(I)=-Y(I) 
            DERIV2(I)=-DERIV2(I) 
   20    CONTINUE 
      ENDIF 
!                                                                       
! --- LOOP TO FIND INITIAL ESTIMATES OF TANGENT POINTS                  
      ERROR=.FALSE. 
      NTAN=0 
      LAST=0 
      DO 200 I=1,N 
         IF(I.LE.LAST)GOTO 200 
         IF(DERIV2(I).LE.0.)THEN 
            DO 100 J=I+1,N 
               IF(DERIV2(J).GT.0..OR.J.EQ.N)THEN 
                  LAST=J 
                  IF(I.EQ.1)THEN 
                     DO 50 K=2,J-1 
                        IF(Y(K).GT.Y(K-1).AND.Y(K).GT.Y(K+1))THEN 
                           NTAN=NTAN+1 
                           IF(NTAN.GT.MAXTAN)THEN 
                              ERROR=.TRUE. 
                              PRINT *,'ERROR in subroutine INITPT: ',   &
     &                                'Too many tangent points --'      
                              PRINT *,'Increase MAXTAN.' 
                              RETURN 
                           ENDIF 
                           INDEX(NTAN)=K 
                        ENDIF 
   50                CONTINUE 
                  ELSEIF(I.GT.1)THEN 
                     IF(J.LT.N)THEN 
                        NTAN=NTAN+1 
                        IF(NTAN.GT.MAXTAN)THEN 
                           ERROR=.TRUE. 
                           PRINT *,'ERROR in subroutine INITPT: ',      &
     &                             'Too many tangent points --'         
                           PRINT *,'Increase MAXTAN.' 
                           RETURN 
                        ENDIF 
                        INDEX(NTAN)=(I+J-1)/2 
                     ELSE 
                        DO 80 K=I,J-1 
                           IF(Y(K).GT.Y(K-1).AND.Y(K).GT.Y(K+1))THEN 
                              NTAN=NTAN+1 
                              IF(NTAN.GT.MAXTAN)THEN 
                                 ERROR=.TRUE. 
                                 PRINT *,'ERROR in subroutine INITPT: ',&
     &                                   'Too many tangent points --'   
                                 PRINT *,'Increase MAXTAN.' 
                                 RETURN 
                              ENDIF 
                              INDEX(NTAN)=K 
                           ENDIF 
   80                   CONTINUE 
                     ENDIF 
                  ENDIF 
                  GOTO 200 
               ENDIF 
  100       CONTINUE 
         ENDIF 
  200 END DO 
!                                                                       
! --- IF BOTTOM ENVELOPE HAS BEEN COMPUTED, SET Y VALUES AND SECOND     
! --- DERIVATIVES BACK TO ORIGINAL VALUES                               
      IF(.NOT.TOP)THEN 
         DO 300 I=1,N 
            Y(I)=-Y(I) 
            DERIV2(I)=-DERIV2(I) 
  300    CONTINUE 
      ENDIF 
!                                                                       
      RETURN 
      END                                           
