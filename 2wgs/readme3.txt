

CONVERT PANAIR (A502) INPUT FILE TO LAWGS FORMAT          \2wgs\readme3.txt

The files for this program are in directory \2wgs on the CD-ROM and in the 
archive file 2wgs.zip that may be downloaded from the PDAS web site.
  readme3.txt      general description
  a5022wgs.f90     the complete source code
  tm85767.pdf      NASA document describing LaWGS

Sample cases for this program are:
  swb.inp          input data for simple wing-body test case
  swb.wgs          output produced by swb.inp

The description of the input to Panair is part of the user's manual at
   http://www.pdas.com/refs/panairUsermanOCR.pdf.


This program asks for the name of the input file. This must be 
a file written in the input format to PanAir (A502).  After reading 
the input data, the program produces a file called a502.wgs that may 
be used as input to 3view, hlp or wgs2wrl.

To compile this program for your computer, use the command
   gfortran a5022wgs.f90 -o a5022wgs.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.


  a5022wgs.exe     the executable for Windows
  a5022wgs.lnx     the executable for Linux
  a5022wgs.mac     the executable for Macintosh OS X (Intel)

