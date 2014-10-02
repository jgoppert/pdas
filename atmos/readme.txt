

PROPERTIES OF THE STANDARD ATMOSPHERE                       \atmos\readme.txt

The files for this program are in the directory \atmos on the CD-ROM and in the
archive file atmos.zip that may be downloaded from the PDAS web site.
  readme.txt      this file of general information
  coesa.txt       definition of the 1976 US Atmosphere Standard to 86 km.
  at62.for        the original grandfather of all these codes
  press.dpr       the program used to compute pressures at boundaries
  us1.prt         atmosphere table (long form) in U.S. units
  si1.prt         atmosphere table (long form) in SI units
  us2.prt         atmosphere table (short form) in U.S. units
  si2.prt         atmosphere table (short form) in SI units
  tables76.f90    program that produces atmTabs.html
  atmTabs.html     combined tables plus transport properties in HTML format
  atmos76.f90     subroutine to compute atmosphere to 1000 km.
  ussa1976.pas    Steven Pietrobon's program to compute upper atmosphere
  ussa1976.dpr    same as ussa1976.pas, but as a Delphi console application
  ussa1976.exe    executable for ussa1976 (using ussa1976.dpr)
  uasa2000.pas    Pietrobon's program to compute approximate upper atmosphere
  uasa2000.dpr    same as uasa2000.pas, but as a Delphi console application
  uasa2000.exe    executable for uasa2000 (using uasa2000.dpr)
  us76.pdf        a copy of U.S. Standard Atmosphere 1976 (Govt. Printing Office)

  bb.f90          the baseball trajectory program mentioned below
  hotcold.f90     non-standard atmospheres (hot,cold,polar,tropical)





Files associated with the "tables" program
  tables.txt      description of "tables" and notes for each language
  language.txt    discussion of the differences between languages
  rlc1.mak        the project file in Visual Basic
  main.frm        the form description and code for Visual Basic
  tables.bas      the source code for tables in QBASIC
  tables.c        the source code for tables in C
  tables.cpp      the source code for tables in C++
  tables.dpr      the source code for tables in Pascal (as a Delphi console app.)
  tables.f90      the source code for tables in Fortran
  tables.for      the source code for tables in Fortran 77
  tables.py       the source code in Python (by Rich Kwan)
  tables.zip      All of the output files (archived)
  ussa.pro        tables written in IDL (by Martin Shultz)
  tables76.f90    Similar to tables.f90, but writes the output as atmTabs.html to
                     be viewed on a web browser (HTML 5)
  bigtable.f90    similar to tables.f90, but makes tables to 1000 km.
  bigtable.out    output from running bigtable


The program "tables", which is given in several languages, produces four 
atmosphere tables: 
  1. a table is US units from 0 to 280000 ft. by 5000 ft.,
  2. a table in SI units form 0 to 86 km by 2 km
  3. a table in US units from 0 to 65000 ft. by 1000 ft., and 
  4. a table in SI units from 0 to 20 km by 0.5 km.
The goal is to produce exactly the same text in all programs, but there
are small differences.

The Pascal program ussa1976.pas (and its Delphi counterpart ussa1976.dpr)
are contributions by Steven Pietrobon that produce the same tables
as "tables", but also extends from 86 to 1000 km. Another program by
Pietrobon called uasa2000.pas (and Delphi version uasa2000.dpr) compute 
the same quantities but use cubic spline curve fits instead of numerical 
integration and solution of differential equations. uasa is an approximation,
but highly accurate.

The program "bigtable" is similar to uasa2000 in that it uses cubic spline 
curve fits above 86 km to approximate the official results. It also makes a
table, but the atmosphere calculations are enclosed in a module you can
use in your computing projects.

The program VuCalc in this collection uses the Atmosphere procedure
to calculate the flight conditions when you specify both altitude and 
Mach number. 

New for version 10 (January 2005) is the file hotcold.f90.
The purpose of this routine is to compute non-standard atmospheres.
This is still in the final stages of checkout, but I decided to include it
if you want to experiment with it. As always, I appreciate your comments.
Unlike the standard atmosphere, there is no official scientific body that
has put an official stamp of approval on a non-standard atmosphere.
What I have coded is from a Department of Defense document MIL-STD-210A.
The algorithm may be stated simply: For a given altitude, use the standard
atmosphere to compute pressure; then get the temperature from the defined
profile for this non-standard day (either hot,cold,polar, or tropical);
finally, compute density from the pressure and temperature using the
perfect gas law.


DESCRIPTION

Every student of aerodynamics and flight mechanics is introduced to the
atmosphere table, which allows the determination of temperature, pressure, 
and density at any altitude. The purpose of this small program is to compute
and print such a table. The equations used are those adopted 15 October 1976 
by the United States Committee on Extension to the Standard Atmosphere(COESA),
representing 29 U.S. scientific and engineering organizations. The values 
selected in 1976 are slight modifications of those adopted in 1962. The 
equations and parameters used are documented in "U.S. Standard Atmosphere, 
1976," published by the U.S. Government Printing Office, Washington, D.C. 

A summary of the definition of the 1976 atmosphere to 86 km. is contained 
in the file COESA.TXT and the calculation of the pressure at each boundary 
altitude is shown in the program pressure.dpr.

Since nearly every introductory aerodynamics textbook contains such
a table, one may question the value of producing yet another one. However,
by completing this exercise, these routines will be placed in your 
standard toolbox. Then, when you are studying a new vehicle concept or 
flight procedure, you can concentrate on your idea and not on validating 
your atmosphere calculations.

The fundamental procedure is a subroutine called Atmosphere that accepts 
altitude as an input argument and returns non-dimensional values of 
temperature, pressure, and density, which are ratios of the quantity at 
altitude to that at sea-level. The equations are taken directly from the 
official publication. Since the definition of the international standard is 
given in SI units, the altitude is supplied in kilometers. The standard 
atmosphere is defined as a set of layers and the routine determines which 
layer contains the specified altitude. The desired layer is found by 
binary search. If you are not familiar with this search technique, it is 
worth a bit of study.  Binary search converges much faster than sequential 
search and the coding is no more difficult. The desired quantity is then 
computed by interpolation.

The routine Atmosphere implements the first seven layers of the atmosphere,
as defined in the 1976 standard. In order to check the operation of the 
subroutine, a program called "tables" computes and prints a 
formatted page of data showing the atmospheric quantities at various 
altitudes. In addition to the three fundamental non-dimensional numbers, 
it also computes their dimensional values and, to fill out the 
page, supplies the speed of sound and the viscosity. Two versions of this 
program are supplied. One prints the data in metric units from 0 to 86 Km; 
the other prints the data in US Customary units from 0 to 280000 ft. 
Neatly formatting such a table presents a few problems in style. The 
temperature in the atmosphere varies gradually and can easily be printed 
in a fixed format. However, the pressure at sea-level is a million times 
greater than the pressure at 86 km. Printing this quantity in fixed format 
is wasteful, so it is printed in exponential format.

As a practical matter, almost all flight takes place in the first two layers
commonly referred to as the troposphere and the stratosphere. The routine
'Atmosphere,' which does a very elegant binary search through the layers,
will generally use the first or second. An alternate routine called 
SimpleAtmosphere is only correct to 20 km., and makes only one simple test 
for altitude. The table generating program has been modified to use this 
procedure to make alternate tables of the lower atmosphere. It will be left 
as an exercise for the student to make an even simpler routine, perhaps 
called Troposphere, which defines only the first layer. This will be quite 
adequate for work in general aviation, model airplanes, gliders, 
helicopters, etc.

As the short tables print in fewer columns, I have added one more quantity
that I find useful, although I have never seen it in another atmosphere
table. This is the ratio of speed of sound to the kinematic viscosity.
Multiplying this number by the chord length and the Mach number yields
the Reynolds number. For example, if I am cruising at M=0.8 at 34000 ft.
in a jet transport with a chord of 16 ft., then from the table at
34000 ft, this ratio is 2.48 1/ft, and so the Reynolds number is 
2.48*16*0.8=31.7 million (the million is built in). I realize that this
only saves one division, but it seems so much faster, probably
because you don't have to deal with powers of 10. In the metric table
this ratio is given in units of 1/m while it is 1/ft in the US table.

A Helpful hint:  If you work in US customary units and you can remember
that this ratio is about 7 at sea level and about 2-3 for the lower
stratosphere and about 0.5 for SST cruise, then you can estimate Reynolds
number in your head. Astound your colleagues!

APPLICATION - As an example of an interesting program that uses the
atmosphere routine, I have included bb.f90. This program computes the
trajectory of a baseball with a given initial velocity, angle with the
horizontal and altitude of the playing field. With bb, you can answer 
questions as to whether a ball travels further in Denver than it does
at sea level or what is the best initial angle to hit a ball [it is
not 45 degrees].

HTML Output - Also note the program tables76.f90. Here, instead of using 
formatted output and monospaced fonts to achieve neat columns, the output 
is simply put in HTML format as tables. The browser does all the work of 
keeping the numbers in neat columns.

NOTE OF THANKS - I want to thank Steve Pietrobon for his contributions to
this collection. Also, many thanks to Martin Schultz and Rich Kwan for the
codings in IDL and Python.


