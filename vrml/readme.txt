

MAKE VRML WORLD FROM LaWGS OBJECTS                           \vrml\readme.txt

The files for this program are in the directory \vrml on the CD-ROM:
  readme.txt     this file. General description
  wgs2wrl.f90    the complete source code

Sample cases for this program are:
  tnd7505.wgs     sample wgs file 
  tnd7505.wrl     output .wrl file from this program
  tnd4211.wgs     sample wgs file 
  tnd4211.wrl     output .wrl file from this program
  wdcase1.wgs     sample wgs file 
  wdcase1.wrl     output .wrl file from this program

To use this program, create a directory on your hard disk and copy these 
files to that directory.

This is an experimental program in the initial phases of testing. I hope 
to learn how to define VRML world objects and eventually couple this to 
some of the aerodynamic programs to display streamlines and surface 
pressures, and maybe more.

You are invited to have a good look and perhaps make some comments that will 
aid me in the development. It produces a file in VRML 1.0 format.

The origin of this program came from several user requests that I offer a
version of a program called FAST, developed at NASA Ames, which in turn 
was an outgrowth of a program called PLOT3D. The idea is to display the data
(input and output) from a general class of CFD programs.

Further investigation indicated that there was a large quantity of code 
specific to the Silicon Graphics workstation; I do not want to start 
supporting programs that required expensive workstations! Actually, the
OpenGL part of the program is becoming available on the PC platform, so 
the project might be feasible.

However, there is an alternate approach that seemed more intriguing, namely
capitializing on the VRML (Virtual Reality Modelling Language) work being 
done by many people.  It seems to me all I have to do is write a file in the 
approved VRML format, and all those smart people at Silicon Graphics and 
Platinum and Microsoft will take care of everything concerning viewing 
angles, lighting, navigation in 3D, colors, etc. 

So, what you have on the CD-ROM is the first effort. The program asks for 
the name of the input file. This must be a file in LaWGS format. After 
reading the input data, the program produces a file called wgs.wrl that 
may be used as input to a VRML browser. 

