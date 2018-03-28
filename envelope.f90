      PROGRAM ENVELOPE 
!                                                                       
! --- This program calculates the top and bottom envelope curves for a  
!     given set of oscillatory data.  The results are plotted on the    
!     screen (using the VG graphics package), and the user may then save
!     the plot in a file for later printing on a higher-quality output  
!     device and/or save the calculated envelope points and tangent     
!     points in a file.                                                 
!                                                                       
! --- M. McClain, Scientific Computing Division, February 1989.         
!     Latest Revision, September 1993 (Version 2.01).                   
!     Changes since Version 2.0 (October 1992):                         
!     (1) For hardcopy VG output, changed color code from 10 to 5,      
!         because 10 is too dark on Sparc printer.                      
!     (2) Changed text for filename extensions from upper case to       
!         lower case.                                                   
!     Changes since previous version:                                   
!     (1) Removed specific references to transmittance data.            
!     (2) Improved handling of error tolerances for subroutine SMOOTH.  
!     (3) Updated subroutine SMOOTH to October 1992 version.            
!     (4) In subroutine ENV, replaced call to PCHIM with call to PCHIC. 
!                                                                       
! --- Declaration of Variables                                          
      CHARACTER FILNAM*20,ANSWER*1,LETTER*1,TEXT1*7,TEXT2*7,DEV*3,      &
     &          NAME*20                                                 
      LOGICAL FILEOK,ERROR,TOP 
      INTEGER I,N,MAXN,MAXBK,NBKPTS,LW,NORD,INBV,IDER,MAXTAN,NTTAN,     &
     &        NBTAN,NTENV,NBENV,NAMLEN,FULLEN,IERR,TRACE,LOOPIN,LOOPBT  
      REAL XPT,YPT,TOL1,TOL2,BVALU,TEMP,SMIN,TOL,YMIN,YMAX,RANGE 
      PARAMETER(MAXN=1200,MAXBK=200,MAXTAN=50,LW=4*MAXN) 
      INTEGER TINDEX(MAXTAN),BINDEX(MAXTAN) 
      REAL X(MAXN),Y(MAXN),YAVE(MAXN),S(MAXN),Z(MAXN),DERIV2(MAXN),     &
     &     WORK1(MAXN),XTENV(MAXN),YTENV(MAXN),XBENV(MAXN),YBENV(MAXN), &
     &     BKPT(MAXBK),COEFF(MAXBK),WORK2(MAXBK),                       &
     &     XTTAN(MAXTAN),YTTAN(MAXTAN),XBTAN(MAXTAN),YBTAN(MAXTAN),     &
     &     WORK3(2*MAXTAN),                                             &
     &     W(LW),W2(20)                                                 
!                                                                       
! --- Set underflows to zero (for Lahey compiler)                       
!      CALL UNDER0(.TRUE.)                                              
!                                                                       
! --- Print program header and version number.                          
      PRINT *,'Envelope Program Version 2.01' 
      PRINT *,'-----------------------------' 
      PRINT * 
!                                                                       
! --- Open data file and save filename information.                     
   50 CONTINUE 
      PRINT *,'What is the name of your data file? ' 
      READ(*,100)FILNAM 
  100 FORMAT(A) 
      INQUIRE(FILE=FILNAM,EXIST=FILEOK) 
      IF(FILEOK)THEN 
         OPEN(1,FILE=FILNAM,STATUS='OLD') 
      ELSE 
         PRINT *,'File does not exist -- try again.' 
         GOTO 50 
      ENDIF 
      I=0 
  150 CONTINUE 
      I=I+1 
      LETTER=FILNAM(I:I) 
      IF(I.LE.LEN(FILNAM).AND.LETTER.NE.' '.AND.LETTER.NE.'.')GOTO 150 
      NAMLEN=I-1 
      NAME=FILNAM(1:NAMLEN) 
      IF(LETTER.EQ.'.')THEN 
  180    CONTINUE 
         I=I+1 
         LETTER=FILNAM(I:I) 
         IF(I.LE.LEN(FILNAM).AND.LETTER.NE.' ')GOTO 180 
         FULLEN=I-1 
      ELSE 
         FULLEN=NAMLEN 
      ENDIF 
!                                                                       
! --- Read in data.                                                     
      PRINT *,'Reading data ... ' 
      N=0 
  200 CONTINUE 
      READ(1,*,END=300)XPT,YPT 
      N=N+1 
      IF(N.GT.MAXN)THEN 
         PRINT *,'ERROR: Too many data points -- increase MAXN.' 
         STOP 
      ENDIF 
      IF(XPT.LT.X(N-1))THEN 
         PRINT *,'ERROR: X value out of order at point ',N,'.' 
         STOP 
      ENDIF 
      X(N)=XPT 
      Y(N)=YPT 
      GOTO 200 
  300 CONTINUE 
      CLOSE(1) 
      WRITE(*,350)N 
  350 FORMAT(' Found ',I5,' data points.') 
!                                                                       
! --- Perform rough smoothing by averaging every three data points.     
      PRINT * 
      PRINT *,'Averaging data ... ' 
      YAVE(1)=Y(1) 
      DO 400 I=2,N-1 
         YAVE(I)=(Y(I-1)+Y(I)+Y(I+1))/3. 
  400 END DO 
      YAVE(N)=Y(N) 
!                                                                       
! --- Find minimum and maximum data values.                             
      YMIN=YAVE(1) 
      YMAX=YAVE(1) 
      DO 420 I=2,N 
         YMIN=MIN(YMIN,YAVE(I)) 
         YMAX=MAX(YMAX,YAVE(I)) 
  420 END DO 
      RANGE=YMAX-YMIN 
!                                                                       
! --- Make initial pass through smoothing algorithm.                    
      TOL1=.08 
      PRINT * 
      WRITE(*,450)TOL1 
  450 FORMAT(' Enter tolerance factor for initial smoothing ',          &
     &        '(default =',F7.4,'): ')                                  
      READ(*,500)TEMP 
  500 FORMAT(F8.0) 
      IF(TEMP.GT.0.)TOL1=TEMP 
      TOL=TOL1*RANGE 
      DO 600 I=1,N 
         S(I)=TOL 
  600 END DO 
      SMIN=TOL 
      PRINT *,'Smoothing data ... ' 
  610 CONTINUE 
      TRACE=0 
      NORD=4 
      NBKPTS=0 
      CALL SMOOTH(IERR,TRACE,NORD,N,X,YAVE,S,SMIN,Z,WORK1,MAXBK,NBKPTS, &
     &            BKPT,COEFF,WORK2,LW,W)                                
      IF(IERR.NE.0)THEN 
         SMIN=2*SMIN 
         WRITE(*,620)IERR 
  620    FORMAT('    Unable to fit data with current parameters ',      &
     &          '(IERR =',I2,') --')                                    
         PRINT *,'   Doubling tolerance and trying again ... ' 
         GOTO 610 
      ENDIF 
!     WRITE(*,650)NBKPTS                                                
! 650 FORMAT(' Initial smoothing required ',I3,' breakpoints.')         
!                                                                       
! --- Initialize interpolation points for envelope algorithm.           
      PRINT * 
      PRINT *,'Initializing tangent points for envelope ... ' 
      NORD=4 
      INBV=1 
      IDER=2 
      DO 700 I=1,N 
         DERIV2(I)=BVALU(BKPT,COEFF,NBKPTS-NORD,NORD,IDER,X(I),INBV,W2) 
  700 END DO 
      TOP=.TRUE. 
      CALL INITPT(ERROR,TOP,N,Z,DERIV2,MAXTAN,NTTAN,TINDEX) 
      IF(ERROR)STOP 
      WRITE(*,750)NTTAN 
  750 FORMAT(' Found ',I2,' tangent points for top envelope.') 
      IF(NTTAN.LE.1)THEN 
         PRINT *,'ERROR: Not enough tangent points -- cannot proceed.' 
         PRINT *,'       (Try a smaller tolerance factor.)' 
         STOP 
      ENDIF 
      TOP=.FALSE. 
      CALL INITPT(ERROR,TOP,N,Z,DERIV2,MAXTAN,NBTAN,BINDEX) 
      IF(ERROR)STOP 
      WRITE(*,800)NBTAN 
  800 FORMAT(' Found ',I2,' tangent points for bottom envelope.') 
      IF(NBTAN.LE.1)THEN 
         PRINT *,'ERROR: Not enough tangent points -- cannot proceed.' 
         PRINT *,'       (Try a smaller tolerance factor.)' 
         STOP 
      ENDIF 
!                                                                       
! --- Make final pass through smoothing algorithm.                      
      TOL2=.01 
      PRINT * 
      WRITE(*,850)TOL2 
  850 FORMAT(' Enter tolerance factor for final smoothing ',            &
     &        '(default =',F7.4,'): ')                                  
      READ(*,900)TEMP 
  900 FORMAT(F8.0) 
      IF(TEMP.GT.0.)TOL2=TEMP 
      IF(TOL2.LT.TOL1)THEN 
         TOL=TOL2*RANGE 
         DO 1000 I=1,N 
            S(I)=TOL 
 1000    CONTINUE 
         SMIN=TOL 
         PRINT *,'Smoothing data ... ' 
 1010    CONTINUE 
         CALL SMOOTH(IERR,TRACE,NORD,N,X,YAVE,S,SMIN,Z,WORK1,MAXBK,     &
     &               NBKPTS,BKPT,COEFF,WORK2,LW,W)                      
         IF(IERR.NE.0)THEN 
            SMIN=2*SMIN 
            WRITE(*,620)IERR 
            PRINT *,'   Doubling tolerance and trying again ... ' 
            NBKPTS=0 
            GOTO 1010 
         ENDIF 
!        WRITE(*,1020)NBKPTS                                            
!1020    FORMAT(' Final smoothing required ',I3,' breakpoints.')        
      ENDIF 
!                                                                       
! --- Calculate envelope curves.                                        
      PRINT * 
      PRINT *,'Calculating top envelope ... ' 
      NORD=4 
      INBV=1 
      IDER=2 
      DO 1050 I=1,N 
         DERIV2(I)=BVALU(BKPT,COEFF,NBKPTS-NORD,NORD,IDER,X(I),INBV,W2) 
 1050 END DO 
      DO 1100 I=1,NTTAN 
         XTTAN(I)=X(TINDEX(I)) 
         YTTAN(I)=Z(TINDEX(I)) 
 1100 END DO 
      TOP=.TRUE. 
      CALL ENV(ERROR,TOP,N,X,Z,DERIV2,NTENV,XTENV,YTENV,NTTAN,XTTAN,    &
     &         YTTAN,WORK1,work3)                                       
      IF(ERROR)STOP 
      PRINT * 
      PRINT *,'Calculating bottom envelope ... ' 
      DO 1400 I=1,NBTAN 
         XBTAN(I)=X(BINDEX(I)) 
         YBTAN(I)=Z(BINDEX(I)) 
 1400 END DO 
      TOP=.FALSE. 
      CALL ENV(ERROR,TOP,N,X,Z,DERIV2,NBENV,XBENV,YBENV,NBTAN,XBTAN,    &
     &         YBTAN,WORK1,work3)                                       
      IF(ERROR)STOP 
!                                                                       
! --- Plot results on screen.                                           
!      PRINT * 
!      PRINT *,'Press <Enter> to proceed with plotting.' 
!      READ(*,*) 
!      CALL NEWPAG 
! 1800 CONTINUE 
!      CALL PLAC(1.,.1,0.,1.) 
!      CALL SIDTEX(FILNAM,11,'X',11,'Y',11,' ',0) 
!      CALL HOWPLT(23,0,1) 
!      CALL CURV(N,X,Y) 
!      CALL HOWPLT(0,1,15) 
!      CALL CURV(N,X,Z) 
!      CALL HOWPLT(0,1,2) 
!      CALL CURV(NTENV,XTENV,YTENV) 
!      CALL HOWPLT(0,1,2) 
!      CALL CURV(NBENV,XBENV,YBENV) 
!      CALL VG 
!      WRITE(TEXT1,1850)TOL1 
! 1850 FORMAT(F7.4) 
!      CALL HTEX('Initial Tolerance ='//TEXT1,.12,.03,0.,1.,13) 
!      WRITE(TEXT2,1850)TOL2 
!      CALL HTEX('Final Tolerance ='//TEXT2,.55,.03,0.,1.,13) 
!      IF(LOOPIN().EQ.1)GOTO 1800 
!                                                                       
! --- Create plot file for high-quality output.                         
! 1855 CONTINUE 
!      PRINT * 
!      PRINT *,'Do you want to create a high-quality plot file (Y/N)? ' 
!      READ(*,1860)ANSWER 
! 1860 FORMAT(A) 
!      IF(ANSWER.EQ.'Y'.OR.ANSWER.EQ.'y')THEN 
!         PRINT *,'Select an output device from the following --' 
!         PRINT *,'   1. Tektronix' 
!         PRINT *,'   2. HPGL Plotter' 
!         PRINT *,'   3. PostScript' 
!         PRINT *,'   4. QMS Lasergrafix' 
!         PRINT *,'Enter the number of your choice: ' 
!        READ(*,*)TEMP 
!         IF(TEMP.EQ.1)THEN 
!            DEV='tek' 
!         ELSEIF(TEMP.EQ.2)THEN 
!            DEV='hpg' 
!         ELSEIF(TEMP.EQ.3)THEN 
!            DEV='pos' 
!         ELSEIF(TEMP.EQ.4)THEN 
!            DEV='qms' 
!         ELSE 
!            PRINT *,'Invalid selection -- Try again.' 
!            GOTO 1855 
!         ENDIF 
!         CALL NEWPAG 
!         CALL SETDV(DEV) 
!         CALL SETFIL(NAME) 
! 1880    CONTINUE 
!         CALL PLAC(1.,.1,0.,1.) 
!         CALL AXCODE(0,0,5) 
!         CALL SIDTEX(FILNAM,5,'X',5,'Y',5,' ',0) 
!         CALL HOWPLT(23,0,5) 
!         CALL CURV(N,X,Y) 
!         CALL HOWPLT(0,1,5) 
!         CALL CURV(N,X,Z) 
!         CALL HOWPLT(0,1,5) 
!         CALL CURV(NTENV,XTENV,YTENV) 
!         CALL HOWPLT(0,1,5) 
!         CALL CURV(NBENV,XBENV,YBENV) 
!         CALL VG 
!         CALL HTEX('Initial Tolerance ='//TEXT1,.12,.03,0.,1.,5) 
!         CALL HTEX('Final Tolerance ='//TEXT2,.55,.03,0.,1.,5) 
!         IF(LOOPBT().EQ.1)GOTO 1880 
!      ENDIF 
!                                                                       
! --- Save results.                                                     
      PRINT * 
      PRINT *,'Do you want to save the envelope data (Y/N)? ' 
      READ(*,1900)ANSWER 
 1900 FORMAT(A) 
      IF(ANSWER.EQ.'Y'.OR.ANSWER.EQ.'y')THEN 
         CALL SAVENV(FILNAM,FULLEN,NAMLEN,TOL1,TOL2,N,X,                &
     &               NTENV,XTENV,YTENV,NBENV,XBENV,YBENV,               &
     &               NTTAN,XTTAN,YTTAN,NBTAN,XBTAN,YBTAN)               
      ENDIF 
!                                                                       
      STOP 
      END                                           
