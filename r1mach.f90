      REAL FUNCTION R1MACH(I) 
!***BEGIN PROLOGUE  R1MACH                                              
!***DATE WRITTEN   790101   (YYMMDD)                                    
!***REVISION DATE  910131   (YYMMDD)                                    
!***CATEGORY NO.  R1                                                    
!***KEYWORDS  MACHINE CONSTANTS                                         
!***AUTHOR  FOX, P. A., (BELL LABS)                                     
!           HALL, A. D., (BELL LABS)                                    
!           SCHRYER, N. L., (BELL LABS)                                 
!***PURPOSE  Returns single precision machine dependent constants       
!***DESCRIPTION                                                         
!                                                                       
!     This is the CMLIB version of R1MACH, the real machine             
!     constants subroutine originally developed for the PORT library.   
!                                                                       
!     R1MACH can be used to obtain machine-dependent parameters         
!     for the local machine environment.  It is a function              
!     subroutine with one (input) argument, and can be called           
!     as follows, for example                                           
!                                                                       
!          A = R1MACH(I)                                                
!                                                                       
!     where I=1,...,5.  The (output) value of A above is                
!     determined by the (input) value of I.  The results for            
!     various values of I are discussed below.                          
!                                                                       
!  Single-Precision Machine Constants                                   
!  R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.            
!  R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.            
!  R1MACH(3) = B**(-T), the smallest relative spacing.                  
!  R1MACH(4) = B**(1-T), the largest relative spacing.                  
!  R1MACH(5) = LOG10(B)                                                 
!***REFERENCES  FOX, P.A., HALL, A.D., SCHRYER, N.L, *FRAMEWORK FOR     
!                 A PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHE-       
!                 MATICAL SOFTWARE, VOL. 4, NO. 2, JUNE 1978,           
!                 PP. 177-188.                                          
!***ROUTINES CALLED  XERROR                                             
!***END PROLOGUE  R1MACH                                                
!                                                                       
      INTEGER SMALL(2) 
      INTEGER LARGE(2) 
      INTEGER RIGHT(2) 
      INTEGER DIVER(2) 
      INTEGER LOG10(2) 
!                                                                       
      REAL RMACH(5) 
!                                                                       
      EQUIVALENCE (RMACH(1),SMALL(1)) 
      EQUIVALENCE (RMACH(2),LARGE(1)) 
      EQUIVALENCE (RMACH(3),RIGHT(1)) 
      EQUIVALENCE (RMACH(4),DIVER(1)) 
      EQUIVALENCE (RMACH(5),LOG10(1)) 
!                                                                       
!                                                                       
!     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T  
!     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T     
!     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).      
!                                                                       
! === MACHINE = IEEE.MOST-SIG-BYTE-FIRST                                
! === MACHINE = IEEE.LEAST-SIG-BYTE-FIRST                               
! === MACHINE = SUN                                                     
! === MACHINE = 68000                                                   
! === MACHINE = 8087                                                    
! === MACHINE = IBM.PC                                                  
! === MACHINE = ATT.3B                                                  
! === MACHINE = ATT.6300                                                
! === MACHINE = ATT.7300                                                
       DATA SMALL(1) /     8388608 / 
       DATA LARGE(1) /  2139095039 / 
       DATA RIGHT(1) /   864026624 / 
       DATA DIVER(1) /   872415232 / 
!                                                                       
!     MACHINE CONSTANTS FOR AMDAHL MACHINES.                            
!                                                                       
! === MACHINE = AMDAHL                                                  
!      DATA SMALL(1) /    1048576 /                                     
!      DATA LARGE(1) / 2147483647 /                                     
!      DATA RIGHT(1) /  990904320 /                                     
!      DATA DIVER(1) / 1007681536 /                                     
!      DATA LOG10(1) / 1091781651 /                                     
!                                                                       
!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.                  
!                                                                       
! === MACHINE = BURROUGHS.1700                                          
!      DATA RMACH(1) / Z400800000 /                                     
!      DATA RMACH(2) / Z5FFFFFFFF /                                     
!      DATA RMACH(3) / Z4E9800000 /                                     
!      DATA RMACH(4) / Z4EA800000 /                                     
!      DATA RMACH(5) / Z500E730E8 /                                     
!                                                                       
!     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS.       
!                                                                       
! === MACHINE = BURROUGHS.5700                                          
! === MACHINE = BURROUGHS.6700                                          
! === MACHINE = BURROUGHS.7700                                          
!      DATA RMACH(1) / O1771000000000000 /                              
!      DATA RMACH(2) / O0777777777777777 /                              
!      DATA RMACH(3) / O1311000000000000 /                              
!      DATA RMACH(4) / O1301000000000000 /                              
!      DATA RMACH(5) / O1157163034761675 /                              
!                                                                       
!     MACHINE CONSTANTS FOR THE CONVEX C-120 (NATIVE MODE)              
!                                                                       
! === MACHINE = CONVEX.C1                                               
!      DATA RMACH(1) / 2.9387360E-39 /                                  
!      DATA RMACH(2) / 1.7014117E+38 /                                  
!      DATA RMACH(3) / 5.9604645E-08 /                                  
!      DATA RMACH(4) / 1.1920929E-07 /                                  
!      DATA RMACH(5) / 3.0102999E-01 /                                  
!                                                                       
!     MACHINE CONSTANTS FOR THE CONVEX C-120 (NATIVE MODE)              
!     WITH -R8 OPTION                                                   
!                                                                       
! === MACHINE = CONVEX.C1.R8                                            
!      DATA RMACH(1) / 5.562684646268007D-309 /                         
!      DATA RMACH(2) / 8.988465674311577D+307 /                         
!      DATA RMACH(3) / 1.110223024625157D-016 /                         
!      DATA RMACH(4) / 2.220446049250313D-016 /                         
!      DATA RMACH(5) / 3.010299956639812D-001 /                         
!                                                                       
!     MACHINE CONSTANTS FOR THE CONVEX C-120 (IEEE MODE)                
!                                                                       
! === MACHINE = CONVEX.C1.IEEE                                          
!      DATA RMACH(1) / 1.1754945E-38 /                                  
!      DATA RMACH(2) / 3.4028234E+38 /                                  
!      DATA RMACH(3) / 5.9604645E-08 /                                  
!      DATA RMACH(4) / 1.1920929E-07 /                                  
!      DATA RMACH(5) / 3.0102999E-01 /                                  
!                                                                       
!     MACHINE CONSTANTS FOR THE CONVEX C-120 (IEEE MODE)                
!     WITH -R8 OPTION                                                   
!                                                                       
! === MACHINE = CONVEX.C1.IEEE.R8                                       
!      DATA RMACH(1) / 2.225073858507202D-308 /                         
!      DATA RMACH(2) / 1.797693134862315D+308 /                         
!      DATA RMACH(3) / 1.110223024625157D-016 /                         
!      DATA RMACH(4) / 2.220446049250313D-016 /                         
!      DATA RMACH(5) / 3.010299956639812D-001 /                         
!                                                                       
!     MACHINE CONSTANTS FOR THE CYBER 170/180 SERIES USING NOS (FTN5).  
!                                                                       
! === MACHINE = CYBER.170.NOS                                           
! === MACHINE = CYBER.180.NOS                                           
!      DATA RMACH(1) / O"00014000000000000000" /                        
!      DATA RMACH(2) / O"37767777777777777777" /                        
!      DATA RMACH(3) / O"16404000000000000000" /                        
!      DATA RMACH(4) / O"16414000000000000000" /                        
!      DATA RMACH(5) / O"17164642023241175720" /                        
!                                                                       
!     MACHINE CONSTANTS FOR THE CDC 180 SERIES USING NOS/VE             
!                                                                       
! === MACHINE = CYBER.180.NOS/VE                                        
!      DATA RMACH(1) / Z"3001800000000000" /                            
!      DATA RMACH(2) / Z"4FFEFFFFFFFFFFFE" /                            
!      DATA RMACH(3) / Z"3FD2800000000000" /                            
!      DATA RMACH(4) / Z"3FD3800000000000" /                            
!      DATA RMACH(5) / Z"3FFF9A209A84FBCF" /                            
!                                                                       
!     MACHINE CONSTANTS FOR THE CYBER 205                               
!                                                                       
! === MACHINE = CYBER.205                                               
!      DATA RMACH(1) / X'9000400000000000' /                            
!      DATA RMACH(2) / X'6FFF7FFFFFFFFFFF' /                            
!      DATA RMACH(3) / X'FFA3400000000000' /                            
!      DATA RMACH(4) / X'FFA4400000000000' /                            
!      DATA RMACH(5) / X'FFD04D104D427DE8' /                            
!                                                                       
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.                   
!                                                                       
! === MACHINE = CDC.6000                                                
! === MACHINE = CDC.7000                                                
!      DATA RMACH(1) / 00014000000000000000B /                          
!      DATA RMACH(2) / 37767777777777777777B /                          
!      DATA RMACH(3) / 16404000000000000000B /                          
!      DATA RMACH(4) / 16414000000000000000B /                          
!      DATA RMACH(5) / 17164642023241175720B /                          
!                                                                       
!     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.                  
!                                                                       
! === MACHINE = CRAY.46-BIT-INTEGER                                     
! === MACHINE = CRAY.64-BIT-INTEGER                                     
!      DATA RMACH(1) / 200034000000000000000B /                         
!      DATA RMACH(2) / 577767777777777777776B /                         
!      DATA RMACH(3) / 377224000000000000000B /                         
!      DATA RMACH(4) / 377234000000000000000B /                         
!      DATA RMACH(5) / 377774642023241175720B /                         
!                                                                       
!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200              
!                                                                       
!     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -      
!     STATIC RMACH(5)                                                   
!                                                                       
! === MACHINE = DATA_GENERAL.ECLIPSE.S/200                              
!      DATA SMALL/20K,0/,LARGE/77777K,177777K/                          
!      DATA RIGHT/35420K,0/,DIVER/36020K,0/                             
!      DATA LOG10/40423K,42023K/                                        
!                                                                       
!     ELXSI 6400                                                        
!                                                                       
! === MACHINE = ELSXI.6400                                              
!      DATA SMALL(1) / '00800000'X /                                    
!      DATA LARGE(1) / '7F7FFFFF'X /                                    
!      DATA RIGHT(1) / '33800000'X /                                    
!      DATA DIVER(1) / '34000000'X /                                    
!      DATA LOG10(1) / '3E9A209B'X /                                    
!                                                                       
!     MACHINE CONSTANTS FOR THE HARRIS 220                              
!     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7              
!                                                                       
! === MACHINE = HARRIS.220                                              
! === MACHINE = HARRIS.SLASH6                                           
! === MACHINE = HARRIS.SLASH7                                           
!      DATA SMALL(1),SMALL(2) / '20000000, '00000201 /                  
!      DATA LARGE(1),LARGE(2) / '37777777, '00000177 /                  
!      DATA RIGHT(1),RIGHT(2) / '20000000, '00000352 /                  
!      DATA DIVER(1),DIVER(2) / '20000000, '00000353 /                  
!      DATA LOG10(1),LOG10(2) / '23210115, '00000377 /                  
!                                                                       
!     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.              
!     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.              
!                                                                       
! === MACHINE = HONEYWELL.600/6000                                      
! === MACHINE = HONEYWELL.DPS.8/70                                      
!      DATA RMACH(1) / O402400000000 /                                  
!      DATA RMACH(2) / O376777777777 /                                  
!      DATA RMACH(3) / O714400000000 /                                  
!      DATA RMACH(4) / O716400000000 /                                  
!      DATA RMACH(5) / O776464202324 /                                  
!                                                                       
!     MACHINE CONSTANTS FOR THE HP 2100                                 
!     3 WORD DOUBLE PRECISION WITH FTN4                                 
!                                                                       
! === MACHINE = HP.2100.3_WORD_DP                                       
!      DATA SMALL(1), SMALL(2) / 40000B,       1 /                      
!      DATA LARGE(1), LARGE(2) / 77777B, 177776B /                      
!      DATA RIGHT(1), RIGHT(2) / 40000B,    325B /                      
!      DATA DIVER(1), DIVER(2) / 40000B,    327B /                      
!      DATA LOG10(1), LOG10(2) / 46420B,  46777B /                      
!                                                                       
!     MACHINE CONSTANTS FOR THE HP 2100                                 
!     4 WORD DOUBLE PRECISION WITH FTN4                                 
!                                                                       
! === MACHINE = HP.2100.4_WORD_DP                                       
!      DATA SMALL(1), SMALL(2) / 40000B,       1 /                      
!      DATA LARGE91), LARGE(2) / 77777B, 177776B /                      
!      DATA RIGHT(1), RIGHT(2) / 40000B,    325B /                      
!      DATA DIVER(1), DIVER(2) / 40000B,    327B /                      
!      DATA LOG10(1), LOG10(2) / 46420B,  46777B /                      
!                                                                       
!     HP 9000                                                           
!                                                                       
!      R1MACH(1) = 1.17549435E-38                                       
!      R1MACH(2) = 1.70141163E+38                                       
!      R1MACH(3) = 5.960464478E-8                                       
!      R1MACH(4) = 1.119209290E-7                                       
!      R1MACH(5) = 3.01030010E-1                                        
!                                                                       
! === MACHINE = HP.9000                                                 
!      DATA SMALL(1) / 00040000000B /                                   
!      DATA LARGE(1) / 17677777777B /                                   
!      DATA RIGHT(1) / 06340000000B /                                   
!      DATA DIVER(1) / 06400000000B /                                   
!      DATA LOG10(1) / 07646420233B /                                   
!                                                                       
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,                     
!     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86 AND                  
!     THE INTERDATA 3230 AND INTERDATA 7/32.                            
!                                                                       
! === MACHINE = IBM.360                                                 
! === MACHINE = IBM.370                                                 
! === MACHINE = XEROX.SIGMA.5                                           
! === MACHINE = XEROX.SIGMA.7                                           
       DATA LOG10(1) /  1050288283 / 
! === MACHINE = XEROX.SIGMA.9                                           
! === MACHINE = SEL.85                                                  
! === MACHINE = SEL.86                                                  
! === MACHINE = INTERDATA.3230                                          
! === MACHINE = INTERDATA.7/32                                          
!      DATA RMACH(1) / Z00100000 /                                      
!      DATA RMACH(2) / Z7FFFFFFF /                                      
!      DATA RMACH(3) / Z3B100000 /                                      
!      DATA RMACH(4) / Z3C100000 /                                      
!      DATA RMACH(5) / Z41134413 /                                      
!                                                                       
!     MACHINE CONSTANTS FOR THE INTERDATA 8/32                          
!     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.                         
!                                                                       
!     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE                    
!     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.                        
!                                                                       
! === MACHINE = INTERDATA.8/32.UNIX                                     
!      DATA RMACH(1) / Z'00100000' /                                    
!      DATA RMACH(2) / Z'7EFFFFFF' /                                    
!      DATA RMACH(3) / Z'3B100000' /                                    
!      DATA RMACH(4) / Z'3C100000' /                                    
!      DATA RMACH(5) / Z'41134413' /                                    
!                                                                       
!     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR).            
!                                                                       
! === MACHINE = PDP-10.KA                                               
! === MACHINE = PDP-10.KI                                               
!      DATA RMACH(1) / "000400000000 /                                  
!      DATA RMACH(2) / "377777777777 /                                  
!      DATA RMACH(3) / "146400000000 /                                  
!      DATA RMACH(4) / "147400000000 /                                  
!      DATA RMACH(5) / "177464202324 /                                  
!                                                                       
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING                   
!     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).                 
!                                                                       
! === MACHINE = PDP-11.32-BIT                                           
!      DATA SMALL(1) /    8388608 /                                     
!      DATA LARGE(1) / 2147483647 /                                     
!      DATA RIGHT(1) /  880803840 /                                     
!      DATA DIVER(1) /  889192448 /                                     
!      DATA LOG10(1) / 1067065499 /                                     
!                                                                       
!      DATA RMACH(1) / O00040000000 /                                   
!      DATA RMACH(2) / O17777777777 /                                   
!      DATA RMACH(3) / O06440000000 /                                   
!      DATA RMACH(4) / O06500000000 /                                   
!      DATA RMACH(5) / O07746420233 /                                   
!                                                                       
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING                   
!     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL).                
!                                                                       
! === MACHINE = PDP-11.16-BIT                                           
!      DATA SMALL(1),SMALL(2) /   128,     0 /                          
!      DATA LARGE(1),LARGE(2) / 32767,    -1 /                          
!      DATA RIGHT(1),RIGHT(2) / 13440,     0 /                          
!      DATA DIVER(1),DIVER(2) / 13568,     0 /                          
!      DATA LOG10(1),LOG10(2) / 16282,  8347 /                          
!                                                                       
!      DATA SMALL(1),SMALL(2) / O000200, O000000 /                      
!      DATA LARGE(1),LARGE(2) / O077777, O177777 /                      
!      DATA RIGHT(1),RIGHT(2) / O032200, O000000 /                      
!      DATA DIVER(1),DIVER(2) / O032400, O000000 /                      
!      DATA LOG10(1),LOG10(2) / O037632, O020233 /                      
!                                                                       
!     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.                   
!                                                                       
! === MACHINE = SEQUENT.BALANCE.8000                                    
!      DATA SMALL(1) / $00800000 /                                      
!      DATA LARGE(1) / $7F7FFFFF /                                      
!      DATA RIGHT(1) / $33800000 /                                      
!      DATA DIVER(1) / $34000000 /                                      
!      DATA LOG10(1) / $3E9A209B /                                      
!                                                                       
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.                     
!                                                                       
! === MACHINE = UNIVAC.1100                                             
!      DATA RMACH(1) / O000400000000 /                                  
!      DATA RMACH(2) / O377777777777 /                                  
!      DATA RMACH(3) / O146400000000 /                                  
!      DATA RMACH(4) / O147400000000 /                                  
!      DATA RMACH(5) / O177464202324 /                                  
!                                                                       
!     MACHINE CONSTANTS FOR THE VAX 11/780                              
!    (EXPRESSED IN INTEGER AND HEXADECIMAL)                             
!  *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS***              
!                                                                       
! === MACHINE = VAX.11/780                                              
!      DATA SMALL(1) /       128 /                                      
!      DATA LARGE(1) /    -32769 /                                      
!      DATA RIGHT(1) /     13440 /                                      
!      DATA DIVER(1) /     13568 /                                      
!      DATA LOG10(1) / 547045274 /                                      
!                                                                       
!  ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS***      
!                                                                       
!      DATA SMALL(1) / Z00000080 /                                      
!      DATA LARGE(1) / ZFFFF7FFF /                                      
!      DATA RIGHT(1) / Z00003480 /                                      
!      DATA DIVER(1) / Z00003500 /                                      
!      DATA LOG10(1) / Z209B3F9A /                                      
!                                                                       
!                                                                       
!***FIRST EXECUTABLE STATEMENT  R1MACH                                  
      IF (I .LT. 1  .OR.  I .GT. 5)                                     &
     &   CALL XERROR ( 'R1MACH -- I OUT OF BOUNDS',25,1,2)              
!                                                                       
      R1MACH = RMACH(I) 
      RETURN 
!                                                                       
      END                                           
