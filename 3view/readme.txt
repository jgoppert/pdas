

THREEVIEW - MAKE PLAN, SIDE, AND REAR VIEWS OF AN OBJECT      \3view\readme.txt

The files for this program are in the directory \3view on the CD-ROM and in the
archive file 3view.zip that may be downloaded from the PDAS web site.
  readme.txt      general description
  3view.f90       the complete source code

Sample cases for this program are:
  tnd7505.wgs     input data for wing-body configuration of NASA TN D-7505
  tnd4211.wgs     input data for wing-body configuration of NASA TN D-4211
  wbf.wgs         simple sample case with vertical fin
  wdcase1.wgs     sample case from wavedrag, case 1

The reference documents for the tnd4211 and tnd7505 configurations may be accessed
from the web page http://www.pdas.com/3viewrefs.html. The NASA report
defining the Langley Wireframe Geometry Standard may be found there as well.
 
This program asks for the name of the input file. This must be 
a file written to conform to the Langley Wire-Frame Geometry Standard.

After reading the input data, the program produces the files
plan.gnu, side.gnu, and rear.gnu which may be used with gnuplot
or with your favorite plotting package.
For example
   gnuplot> plot 'plan.gnu' with lines

You can use the command
   gnuplot> set size ratio -1
to insure that the pictures are not distorted. Otherwize, they will
be stretched both ways to fill your screen.

To compile this program for your machine, use the command
   gfortran  3view.f90 -o 3view.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.

The executables are included on the CD-ROM
  3view.exe       the executable for Windows
  3view.lnx       the executable for Linux
  3view.mac       the executable for Macintosh OS X (Intel)
