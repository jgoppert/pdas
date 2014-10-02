

INPUT PRE-PROCESSOR FOR PANAIR - PANIN                   /panin/readme.txt

The files for this program are in the directory \panin on the CD-ROM and in the
archive file panin.zip that may be downloaded from the PDAS web site.

  readme.txt      this file - general description
  input.txt       instructions for use
  panin.f90       the source code for panin
  case1a.aux      sample input for panin
  case1a.wgs      the associated LaWGS file (input)
  case1a.inp      the resulting input file for PanAir

To compile this program for your computer, use the command
   gfortran  panin.f90 -o panin.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.


Panair Input Format
The input to Panair is described in the user's manual. The input
data is organized in specific columns. Editing an input file is a rather
error-prone process. The PANIN program was written to enable a user to 
select the flow properties and all other program options by editing a short 
free-format file called an auxiliary file. One entry in the auxiliary file 
is the name of a file that contains the geometrical information. The format 
of this file is that of the NASA standard for wireframe geometry as described 
in NASA Report TM 85767. This file is usually referred to as a WGS file, standing for 
Wireframe Geometry Standard. The program reads the various items of control 
information from the auxiliary file and combines this information with the 
panel geometry in the WGS file to produce a combined file that is a properly 
formatted input file for PanAir.
