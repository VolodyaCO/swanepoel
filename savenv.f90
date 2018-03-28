      SUBROUTINE SAVENV(FILNAM,FULLEN,NAMLEN,TOL1,TOL2,N,X,             &
     &                  NTENV,XTENV,YTENV,NBENV,XBENV,YBENV,            &
     &                  NTTAN,XTTAN,YTTAN,NBTAN,XBTAN,YBTAN)            
!                                                                       
! --- THIS SUBROUTINE SAVES THE RESULTS OF THE ENVELOPE COMPUTATIONS IN 
!     TWO FILES.  A FILE WITH EXTENSION ".ENV" CONTAINS THE TOP AND     
!     BOTTOM ENVELOPE VALUES FOR EACH X VALUE IN THE ORIGINAL DATA FILE 
!     (WITH AN OPTION TO SKIP OVER SOME X VALUES).  A FILE WITH         
!     EXTENSION ".TAN" CONTAINS THE X VALUES FOR EACH TANGENT POINT     
!     BETWEEN THE DATA CURVE AND THE ENVELOPE CURVES, ALONG WITH THE    
!     CORRESPONDING TOP AND BOTTOM ENVELOPE VALUES.                     
!                                                                       
! --- DESCRIPTION OF INPUT VARIABLES                                    
!     FILNAM:       FULL NAME OF USER'S INPUT DATA FILE                 
!     FULLEN:       NUMBER OF CHARACTERS IN FULL FILE NAME              
!     NAMLEN:       NUMBER OF CHARACTERS IN FILE NAME, NOT INCLUDING    
!                   EXTENSION                                           
!     TOL1:         TOLERANCE PARAMETER FOR INITIAL DATA SMOOTHING      
!     TOL2:         TOLERANCE PARAMETER FOR FINAL DATA SMOOTHING        
!     N:            NUMBER OF DATA POINTS IN USER'S INPUT FILE          
!     X(N):         ARRAY CONTAINING X COORDINATES OF USER'S DATA       
!     NTENV:        NUMBER OF CALCULATED POINTS IN TOP ENVELOPE         
!     XTENV(NTENV): ARRAY CONTAINING X COORDINATES OF TOP ENVELOPE      
!                   POINTS                                              
!     YTENV(NTENV): ARRAY CONTAINING Y COORDINATES OF TOP ENVELOPE      
!                   POINTS                                              
!     NBENV:        NUMBER OF CALCULATED POINTS IN BOTTOM ENVELOPE      
!     XBENV(NBENV): ARRAY CONTAINING X COORDINATES OF BOTTOM ENVELOPE   
!                   POINTS                                              
!     YBENV(NBENV): ARRAY CONTAINING Y COORDINATES OF BOTTOM ENVELOPE   
!                   POINTS                                              
!     NTTAN:        NUMBER OF TANGENT POINTS BETWEEN TOP ENVELOPE AND   
!                   DATA CURVE                                          
!     XTTAN(NTTAN): ARRAY CONTAINING X COORDINATES OF TOP ENVELOPE      
!                   TANGENT POINTS                                      
!     YTTAN(NTTAN): ARRAY CONTAINING Y COORDINATES OF TOP ENVELOPE      
!                   TANGENT POINTS                                      
!     NBTAN:        NUMBER OF TANGENT POINTS BETWEEN BOTTOM ENVELOPE AND
!                   DATA CURVE                                          
!     XBTAN(NBTAN): ARRAY CONTAINING X COORDINATES OF BOTTOM ENVELOPE   
!                   TANGENT POINTS                                      
!     YBTAN(NBTAN): ARRAY CONTAINING Y COORDINATES OF BOTTOM ENVELOPE   
!                   TANGENT POINTS                                      
!                                                                       
! --- DECLARATION OF CALLING SEQUENCE VARIABLES                         
      CHARACTER FILNAM*20 
      INTEGER FULLEN,NAMLEN,N,NTENV,NBENV,NTTAN,NBTAN 
      REAL TOL1,TOL2,X(*),XTENV(*),YTENV(*),XBENV(*),YBENV(*),XTTAN(*), &
     &     YTTAN(*),XBTAN(*),YBTAN(*)                                   
!                                                                       
! --- DECLARATION OF INTERNAL VARIABLES                                 
      CHARACTER ANSWER*1,FILE2*20,FILE3*20 
      INTEGER I,STEP,LENGTH,OFFSET,TCOUNT,BCOUNT,TPNTR,BPNTR 
!                                                                       
! --- SAVE ENVELOPE DATA IN ".ENV" FILE                                 
      FILE2=FILNAM(1:NAMLEN)//'.env' 
      LENGTH=NAMLEN+4 
      OPEN(2,FILE=FILE2,STATUS='UNKNOWN') 
      PRINT *,'Do you want to save envelope values for every x ',       &
     &        'value in your '                                          
      PRINT *,'original data file (Y/N)? ' 
      READ(*,2000)ANSWER 
 2000 FORMAT(A) 
      IF(ANSWER.EQ.'Y'.OR.ANSWER.EQ.'y')THEN 
         STEP=1 
      ELSE 
         PRINT *,'Enter an increment for stepping through the x ',      &
     &           'values in your file: '                                
         READ(*,*)STEP 
      ENDIF 
      CALL SEARCH(XTENV(1),N,X,TCOUNT) 
      IF(TCOUNT.EQ.0)THEN 
         PRINT *,'ERROR: X value',XTENV(1),' not found in original ',   &
     &           'data --'                                              
         PRINT *,'       Contact programmer for correction.' 
      ENDIF 
      OFFSET=MOD(TCOUNT,STEP) 
      TCOUNT=MOD(STEP-OFFSET+1,STEP)+1 
      IF(TCOUNT.GT.NTENV)TCOUNT=1 
      CALL SEARCH(XBENV(1),N,X,BCOUNT) 
      IF(BCOUNT.EQ.0)THEN 
         PRINT *,'ERROR: X value',XBENV(1),' not found in original ',   &
     &           'data --'                                              
         PRINT *,'       Contact programmer for correction.' 
      ENDIF 
      OFFSET=MOD(BCOUNT,STEP) 
      BCOUNT=MOD(STEP-OFFSET+1,STEP)+1 
      IF(BCOUNT.GT.NBENV)BCOUNT=1 
      DO 2400 I=1,N,STEP 
         IF(XTENV(TCOUNT).EQ.X(I))THEN 
            IF(XBENV(BCOUNT).EQ.X(I))THEN 
               WRITE(2,2300)X(I),YTENV(TCOUNT),YBENV(BCOUNT) 
 2300          FORMAT(2X,F8.2,12X,F7.4,12X,F7.4) 
               IF(TCOUNT+STEP.LE.NTENV)TCOUNT=TCOUNT+STEP 
               IF(BCOUNT+STEP.LE.NBENV)BCOUNT=BCOUNT+STEP 
            ELSE 
               WRITE(2,2320)X(I),YTENV(TCOUNT) 
 2320          FORMAT(2X,F8.2,12X,F7.4,14X,''' ''') 
               IF(TCOUNT+STEP.LE.NTENV)TCOUNT=TCOUNT+STEP 
            ENDIF 
         ELSE 
            IF(XBENV(BCOUNT).EQ.X(I))THEN 
               WRITE(2,2350)X(I),YBENV(BCOUNT) 
 2350          FORMAT(2X,F8.2,14X,''' ''',14X,F7.4) 
               IF(BCOUNT+STEP.LE.NBENV)BCOUNT=BCOUNT+STEP 
            ELSE 
               WRITE(2,2380)X(I) 
 2380          FORMAT(2X,F8.2,14X,''' ''',16X,''' ''') 
            ENDIF 
         ENDIF 
 2400 END DO 
      WRITE(2,2450) 
 2450 FORMAT('''',5X,'X',12X,'Top Envelope',6X,'Bottom Envelope''') 
      WRITE(2,*) 
      WRITE(2,*)'''Envelope data generated by program ENVELOPE ',       &
     &          'using data from file ',FILNAM(1:FULLEN),''''           
      WRITE(2,2500)TOL1,TOL2 
 2500 FORMAT('''with initial tolerance ',F7.4,' and final ',            &
     &       'tolerance ',F7.4,'''')                                    
      WRITE(2,*)'''(A '' '' entry indicates no envelope point was ',    &
     &          'generated.)'''                                         
      CLOSE(2) 
      PRINT *,'Envelope data saved in file  ',FILE2(1:LENGTH),'.' 
!                                                                       
! --- SAVE TANGENT POINT DATA IN ".TAN" FILE                            
      FILE3=FILNAM(1:NAMLEN)//'.tan' 
      OPEN(3,FILE=FILE3,STATUS='UNKNOWN') 
      TPNTR=1 
      BPNTR=1 
 2800 CONTINUE 
      IF(TPNTR.NE.NTTAN+1)THEN 
         IF((BPNTR.EQ.NBTAN+1).OR.                                      &
     &     ((BPNTR.NE.NBTAN+1).AND.                                     &
     &     (XTTAN(TPNTR).LT.XBTAN(BPNTR))))THEN                         
            CALL SEARCH(XTTAN(TPNTR),NBENV,XBENV,BCOUNT) 
            IF(BCOUNT.EQ.0)THEN 
               WRITE(3,2820)XTTAN(TPNTR),YTTAN(TPNTR) 
 2820          FORMAT(2X,F8.2,12X,F7.4,14X,''' ''') 
            ELSE 
               WRITE(3,2850)XTTAN(TPNTR),YTTAN(TPNTR),                  &
     &            YBENV(BCOUNT)                                         
 2850          FORMAT(2X,F8.2,12X,F7.4,12X,F7.4) 
            ENDIF 
            TPNTR=TPNTR+1 
            GOTO 2800 
         ENDIF 
      ENDIF 
      IF(BPNTR.NE.NBTAN+1)THEN 
         IF((TPNTR.EQ.NTTAN+1).OR.                                      &
     &     ((TPNTR.NE.NTTAN+1).AND.                                     &
     &      (XTTAN(TPNTR).GE.XBTAN(BPNTR))))THEN                        
            CALL SEARCH(XBTAN(BPNTR),NTENV,XTENV,TCOUNT) 
            IF(TCOUNT.EQ.0)THEN 
               WRITE(3,2880)XBTAN(BPNTR),YBTAN(BPNTR) 
 2880          FORMAT(2X,F8.2,14X,''' ''',14X,F7.4) 
            ELSE 
               WRITE(3,2850)XBTAN(BPNTR),YTENV(TCOUNT),                 &
     &            YBTAN(BPNTR)                                          
            ENDIF 
            BPNTR=BPNTR+1 
            GOTO 2800 
         ENDIF 
      ENDIF 
      IF((TPNTR.NE.NTTAN+1).OR.(BPNTR.NE.NBTAN+1))THEN 
         PRINT *,'ERROR: Incorrect merge of tangent points --' 
         PRINT *,'       Contact programmer for correction.' 
      ENDIF 
      WRITE(3,2900) 
 2900 FORMAT('''',5X,'X',17X,'Tmax',15X,'Tmin''') 
      WRITE(3,*) 
      WRITE(3,*)'''Tangent data generated by program ENVELOPE ',        &
     &          'using data from file ',FILNAM(1:FULLEN),''''           
      WRITE(3,2500)TOL1,TOL2 
      CLOSE(3) 
      PRINT *,'Tangent points saved in file ',FILE3(1:LENGTH),'.' 
!                                                                       
      RETURN 
      END                                           
!                                                                       
!     ******************************************************************
!                                                                       
      SUBROUTINE SEARCH(X,N,LIST,POS) 
!                                                                       
! --- THIS SUBROUTINE PERFORMS A BINARY LIST SEARCH.                    
!                                                                       
! --- INPUT VARIABLES                                                   
!     X:       VALUE TO SEARCH FOR                                      
!     N:       NUMBER OF ITEMS IN LIST                                  
!     LIST(N): ARRAY CONTAINING LIST ITEMS, IN INCREASING ORDER         
!                                                                       
! --- OUTPUT VARIABLE                                                   
!     POS:     POSITION OF X IN LIST -- LIST(POS)=X                     
!              (POS=0 IF X WAS NOT FOUND)                               
!                                                                       
! --- DECLARATION OF CALLING SEQUENCE VARIABLES                         
      INTEGER N,POS 
      REAL X,LIST(*) 
!                                                                       
! --- DECLARATION OF INTERNAL VARIABLES                                 
      INTEGER LEFT,RIGHT,MIDPT 
!                                                                       
! --- SEARCH FOR X IN LIST                                              
      POS=0 
      LEFT=1 
      RIGHT=N 
      IF((X.LT.LIST(1)).OR.(X.GT.LIST(N)))RETURN 
  100 CONTINUE 
      IF(LEFT.LT.RIGHT)THEN 
         MIDPT=(LEFT+RIGHT)/2. 
         IF(X.EQ.LIST(MIDPT))THEN 
            POS=MIDPT 
            RETURN 
         ELSEIF(X.LT.LIST(MIDPT))THEN 
            RIGHT=MIDPT 
         ELSE 
            LEFT=MIDPT 
         ENDIF 
         GOTO 100 
      ELSEIF(LEFT.EQ.RIGHT)THEN 
         MIDPT=LEFT 
         IF(X.EQ.LIST(MIDPT))THEN 
            POS=MIDPT 
            RETURN 
         ENDIF 
      ENDIF 
!                                                                       
      RETURN 
      END                                           
