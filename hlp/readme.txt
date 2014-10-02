

HIDDEN LINE PROGRAM                                           \hlp\readme.txt

The files for this program are in the directory \hlp on the CD-ROM and in the
archive file hlp.zip that may be downloaded from the PDAS web site.

  readme.txt      general description
  hlp.f90         the complete source code
  drh.cmn         you will need it if you want to recompile hlp
  input.txt       instructions for preparing input & interpreting output
  arc12721.txt    the original program description from COSMIC


Sample cases for this program are:
  wbf.wgs         simple sample case with vertical fin
  tnd7505.wgs     wing-body configuration of NASA TN D-7505
  tnd4211.wgs     wing-body configuration of NASA TN D-4211
  f16xl.wgs       F16-XL configurations
  ex1.wgs,ex2.wgs,ex3.wgs,ex4.wgs   examples from the LaWgs document TM 85767

  dpr.xls         Viewing program (from Dan Raymer) (alternate to gnuplot)
  danplot.exe     2008 version of dpr.xls (also from Dan Raymer) (Thanks, Dan)

To use this program, create a directory on your hard disk and copy these 
files to that directory.

This program first asks for the name of the input file. This must
be a file written to conform to the Langley Wire-Frame Geometry Standard.
Next, the program asks for three viewing angles, which correspond to
yaw, roll, and pitch. Enter values in degrees. To complete the viewing
definition, you need to specify the distance of the eye from the object.
Unless you are looking for an extreme perspective effect, give a
number that is large compared to the dimensions of the object. The
image will still be adjusted to fill the screen.

After calculating the scene, the program produces a file called hlp.gnu
that defines the vectors in 2-D that describe the object. Each line in 
hlp.plt defines a point (x,y). Whenever there is a blank line in the file,
this means that the line to the next point is invisible (from the old days 
of pen plotters, this is called a pen-up draw). This file may be displayed 
directly with gnuplot with the command

    gnuplot> plot 'hlp.plt' with lines

You will need to use the command
   gnuplot> set size ratio -1  
to insure that you have a square viewport on your screen to match the
square window of the output from hlp.

There is another file called hlp.ps that is produced by this program.
This is the same picture, but encoded in PostScript format. There is a 
text file on this CD called ps.txt that may give you some tips on using
a PostScript file.

If you do not have gnuplot or a PostScript interpreter installed on your 
system, there is yet another way to visualize the output. Try out the Excel
program called dpr.xls (for Dan Raymer, who wrote the program). This requires
you to have Microsoft Excel on your system, though. Dan wrote a new version
that will work even if you do not have Excel. 

Description -
The program SKETCH (later superseded by SILH), written by David Hedgley 
of the NASA Dryden Flight Research Center is famous as the most widely 
distributed code of the NASA distribution service known as COSMIC. People all
over the world have obtained this program and adapted it to the task of 
displaying three-dimensional objects as wire-frame pictures in perspective 
with hidden lines removed.

Some now say that this program is obsolete because ray-tracing programs can 
create continuous tone color images that are remarkably realistic. But, for 
those of us with black and white laser- or ink-jet printers, this 
procedure still appears to be the best way to display an object like an 
airplane so that it can be visualized. You probably want to learn how to use 
the ray-tracing program so you can make a gee-whiz color figure for a 
presentation, but for the report that will be reproduced on an black-and-white
copy machine, this is the program to use.

I remember jumping at the chance to use this program in 1982, when I was 
looking for a way to check the input for a panel code for aerodynamics known 
as PANAIR. I composed a simple program that read the PANAIR input and used 
SKETCH to view the configuration. Since that time, NASA has settled on a 
single consistent format for describing wireframe objects which is outlined 
in NASA TM 85767, and is referred to as the Langley Wire Frame Geometry 
Standard (LaWGS). This program is an update of the old program PANSKETCH, now
updated to read LaWGS files and use the newer SILH program from David Hedgley.
The output is a plot file, encoded for gnuplot, of the object in question.

