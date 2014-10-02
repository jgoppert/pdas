

HIGH ORDER PANEL AERODYNAMICS CODE (PANAIR)         \panair\readme.txt

The files for this program are in the directory \panair on the CD-ROM:
  readme.txt      general description
  panair.exe      executable for Windows
  panair.lnx      executable for Linux
  panair.mac      executable for Macintosh OS X (Intel)
  usermanOCR.pdf  The user guide
  source16.zip    source code, all 518 subroutines in separate files
  arc11398.txt    the original COSMIC program description
  clean502.bat    cleans up all the intermediate files from previous run
  clean502.sh     same as clean502.bat, but for Mac or Linux

Sample cases for this program
  swb.inp         simple wing-body (same as Ch.9 in panair.pdf)
  swb.out         output from swb.inp
  btac.inp        another sample case
  btac.out        output from btac.inp
  .
  .               ditto for ellip,nac6,ppbc,vbc,vepd,vss
  .
  caseslnx.zip    the above cases with Unix end-of-line (Zip archive) 

The reference documents for this program may be accessed
from the web page http://www.pdas.com/panairrefs.html. 

If you need to recompile this program for your machine, use the command
   gfortran  panair.f90 -o panair.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.

This program first asks for the name of the input file. This must
be a file written to conform to the PanAir user guide.  After calculating 
the solution, the program produces a file called panair.out that contains a 
wealth of information concerning the flow problem to be solved. You need 
to read userman.pdf before you try to interpret the output.

One little gotcha associated with this program is that it will not start
if there are left-over files from a previous run -- and it leaves a lot
of debris behind. There is a command file called clean502.bat that will
delete all of these files.  Some of these files could be valuable for 
post-processing, but I don't have any interpreter programs at this time. 
There is an alternate file called clean502.sh that does the same job 
for Mac or Linux folks.


                     DESCRIPTION

The PanAir program (which is widely known by its code name, A502) has
long been recognized as the most accurate and versatile of the long
series of panel codes in use throughout the aeronautical community.
The program was once one of the most demanding tasks for a major
supercomputer center and it is a tribute to the progress of the personal
computer industry that it is feasible to run this program on a PC.

The history of this program goes back to the mid 70's and I am proud to
say that I had a significant part in its development. Numerous people
contributed to the project both within 
the government and the Boeing Company. The concept originated in the 
aerodynamics research group at Boeing and was funded through several studies
leading to a pilot computer program that demonstrated the feasibility of 
the technology. This led to a large contract development effort that was
funded by several agencies in NASA and the Department of Defense. The
effort went on for several years and led to the production version of the 
program. At that time, the program was transmitted to COSMIC and installed
as numerous sites and became rather well known throughout the
aeronautical community.

In the early 1980s, a black cloud appeared over this sunny scene. The
major airplane companies and government labs were changing from Control 
Data computers, on which PanAir was developed, to Cray supercomputers. 
When it was time to convert PanAir to the Cray, it became apparent that
the program had a vast amount of machine-specific code that required a 
significant effort to convert. The level of government funding for panel 
codes was plummeting at the same time. As a result, the program entered 
a period where it was not usable.

The Boeing Company recognized the value of the program and made the 
conversion to the Cray. The new program was made available to the original 
patrons of the project, but contained so much privately funded work that 
it could no longer be considered a publicly funded program. So, the 
distribution was quite limited. This public/private status continued for 
several years and made the agressive distribution of the program impractical.

In 1993, NASA and the Boeing Company entered into a cooperative agreement,
in which NASA gave Boeing access to the Ames Unitary wind tunnels for
a joint research project in exchange for granting rights to NASA to
distribute the current version of the PanAir program. This effectively
restored the program, at what was internally known as the "ht2" level,
to public domain.

This ht2 version of the program is the basis for the program on this
CD-ROM.  I am sure that the program in use at Boeing has undergone
further development since 1993 so those of you who need the very
latest version can pursue this with Boeing and NASA. My intent is to
offer a version that is free from restrictions on distribution rather
than being the very latest in technology.

All of the executable files are made with the gfortran compiler from GNU.
If you wish to recompile the program, you should unzip the file 
panair.zip and use the following command

  gfortran *.f90

If you are using Lahey LF95, you should use
  lf95  -nf95 -nchk  *.f90
to tell your compiler to ignore certain Fortran 95 violations.


