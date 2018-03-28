To compile the fortran code to find the top and bottom envelopes you will the files: 
envelope.f90, smooth.f90, efc.f90, initpt.f90, env.f90, pchic.f90, pchfe.f90, savenv.f90 
and r1mach.f90 
The aforementioned were modified (using convert.f90) from those that are in the shell archive envelope.shar
(you can find it in https://gams.nist.gov/cgi-bin/serve.cgi/Module/ENVELOPE/ENVELOPE/11673) 
To compile with gfrotran write down in the terminal: 
 gfortran -w envelope.f90 smooth.f90 efc.f90 initpt.f90 env.f90 pchic.f90 pchfe.f90 savenv.f90 r1mach.f90 
This will generate an executable a.out
