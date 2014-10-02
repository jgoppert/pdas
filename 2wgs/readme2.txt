

CONVERT WAVEDRAG INPUT FILE TO LAWGS FORMAT                   \2wgs\readme2.txt

The files for this program are in the directory \2wgs on the CD-ROM
  readme2.txt     general description.
  wd2wgs.f90      the complete source code
  tm85767.pdf     NASA document describing LaWGS

Sample cases for this program are:
  wdcase1.inp     input data for a WaveDrag test case
  wdcase1.wgs     output data produced by wd2wgs

The description of the input to wavedrag may be accessed
from the web page http://www.pdas.com/2wgsrefs.html and the description
of the Langley Wireframe Geometry Standard may be accessed from this page as well.

This program asks for the name of the input file. This must be a file 
written in the input format to WaveDrag. After reading the input data, 
the program produces a file called wd.wgs that may be used as input to 
3view, hlp or wgs2wrl. A file called wd2wgs.dbg gives a printed description of the 
configuration. It is occasionally useful; it is OK to delete it.

To compile this program for your computer, use the command
   gfortran  wd2wgs.f90 -o wd2wgs.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.

  wd2wgs.exe      the executable for Windows
  wd2wgs.lnx      the executable for Linux
  wd2wgs.mac      the executable for Macintosh OS X (Intel)
