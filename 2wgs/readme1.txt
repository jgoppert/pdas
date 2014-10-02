

CONVERT WINGBODY INPUT FILE TO LAWGS FORMAT               \2wgs\readme1.txt

The files for this program are in the directory \2wgs on the CD-ROM:
  readme1.txt     general description
  wb2wgs.f90      the complete source code
  tm85767.pdf     NASA document describing LaWGS

Sample cases for this program are:
  tnd4211.inp     input  data for wing-body configuration of NASA TN D-4211
  tnd7505.inp     input  data for wing-body configuration of NASA TN D-7505
  tnd4211.wgs     output data for wing-body configuration of NASA TN D-4211
  tnd7505.wgs     output data for wing-body configuration of NASA TN D-7505

The reference documents for the tnd4211 and tnd7505 configurations may be accessed
from the web page http://www.pdas.com/2wgsrefs.html and the description
of the input to wingbody may be accessed from this page as well.
 
This program asks for the name of the input file. This must be 
a file written in the input format to WingBody. After reading 
the input data, the program produces a file called wb.wgs that may be 
used as input to 3view, hlp or wgs2wrl.

To compile this program for your machine, use the command
   gfortran  wb2wgs.f90 -o wb2wgs.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.

  wb2wgs.exe      the executable for Windows
  wb2wgs.lnx      the executable for Linux
  wb2wgs.mac      the executable for Macintosh OS X (Intel)
