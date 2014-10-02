' PROGRAM Tables
' PURPOSE - Make tables of atmosphere properties
' AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
' NOTES - This is a QBASIC program, not Visual Basic
' REVISION HISTORY
'   DATE  VERS PERSON  STATEMENT OF CHANGES
'  4Jul95  1.0   RLC   Translated from versions in other languages
' 12Jul95  1.1   RLC   Added viscosity calculations
'  6Aug95  1.2   RLC   Replaced 1962 tables with 1976 tables
' 25Aug95  1.3   RLC   Replaced some WRITE statements with PRINT
' 19Sep95  1.4   RLC   Added some accuracy to the pressure tables
'
DECLARE SUB LongEnglish ()
DECLARE SUB ShortEnglish ()
DECLARE SUB LongMetric ()
DECLARE SUB ShortMetric ()
DECLARE SUB Atmosphere (alt AS SINGLE, sigma AS SINGLE, delta AS SINGLE, theta AS SINGLE)
DECLARE SUB SimpleAtmosphere (alt AS SINGLE, sigma AS SINGLE, delta AS SINGLE, theta AS SINGLE)
DECLARE FUNCTION LookUp! (n AS INTEGER, x() AS SINGLE, u AS SINGLE)
DECLARE FUNCTION MetricViscosity# (theta AS SINGLE)
 
  CONST REARTH = 6369.1#          ' radius of the Earth (km)
  CONST GMR = 34.163195#          ' gas constant
  
  CONST FT2METERS = .3048#        '   mult. ft. to get meters (exact)
  CONST KELVIN2RANKINE = 1.8#
  CONST PSF2NSM = 47.880258#      ' mult lb/sq.ft to get N/sq.m
  CONST SCF2KCM = 515.379#        ' mult slug/cu.ft to get kg/cu.m


  CONST BETAVISC = 1.458E-06      ' viscosity term, N sec/(sq.m sqrt(deg K)
  CONST SUTH = 110.4              ' Sutherland's constant, deg K
  CONST TZERO = 288.15            ' temperature at sealevel, deg K
  CONST PZERO = 101325!           ' pressure at sealevel, N/sq.m.
  CONST RHOZERO = 1.225           ' density at sealevel, kg/cu.m.
  CONST ASOUNDZERO = 340.294      ' speed of sound at sealevel, m/sec

  DIM altKm AS SINGLE
  DIM sigma, delta, theta AS SINGLE
  DIM temp, pressure, density, speedSound AS SINGLE


  DIM viscosity, kinematicViscosity, vratio   AS SINGLE

  DIM SHARED htab(9) AS SINGLE
  DIM SHARED ttab(9) AS SINGLE
  DIM SHARED ptab(9) AS SINGLE
  DIM SHARED gtab(9) AS SINGLE
  FOR i = 1 TO 8
    READ htab(i)
  NEXT i
  FOR i = 1 TO 8
    READ ttab(i)
  NEXT i
  FOR i = 1 TO 8
    READ ptab(i)
  NEXT i
  FOR i = 1 TO 8
    READ gtab(i)
  NEXT i
  DATA 0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852
  DATA 288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.87
  DATA 1.0, 2.233611E-1, 5.403295E-2, 8.5666784E-3, 1.0945601E-3, 6.6063531E-4, 3.9046834E-5, 3.68501E-6
  DATA -6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0

  PRINT "tables - A QBasic program to compute atmosphere tables"
  PRINT "Ralph Carmichael, Public Domain Aeronautical Software"
  PRINT "Modified by   "   ' Put your name here if make changes
  PRINT "Version 1.4  19Sep95"
  LongEnglish
  ShortEnglish
  LongMetric
  ShortMetric

END   ' ----------------------------------------------- End of Program tables

SUB Atmosphere (alt AS SINGLE, sigma AS SINGLE, delta AS SINGLE, theta AS SINGLE)
  CONST REARTH = 6369.1#   ' radius of the Earth (km)
  CONST GMR = 34.163195#   ' gas constant

  DIM i, j, k AS SINGLE
  DIM h, tgrad, deltah, tbase, tlocal AS SINGLE

  h = alt * REARTH / (alt + REARTH)      ' geometric to geopotential altitude

'  n = LookUp(9, htab(), h)                ' find layer containing this height
'  IF n = 0 THEN n = 1
  
  i = 1                                   ' begin binary search
  j = 8
  DO
    k = (i + j) \ 2                        ' integer division (truncates)
    IF h < htab(k) THEN j = k ELSE i = k
  LOOP WHILE j > i + 1
  
  tgrad = gtab(i)
  tbase = ttab(i)
  deltah = h - htab(i)
  tlocal = tbase + tgrad * deltah
  theta = tlocal / ttab(1)                                ' temperature ratio

  IF tgrad = 0# THEN                                         ' pressure ratio
    delta = ptab(i) * EXP(-GMR * deltah / tbase)
  ELSE
    delta = ptab(i) * (tbase / tlocal) ^ (GMR / tgrad)
  END IF

  sigma = delta / theta                                       ' density ratio

END SUB   ' ------------------------------------ End of Procedure Atmosphere.

SUB LongEnglish
  DIM i AS INTEGER

  OPEN "e1b.prt" FOR OUTPUT AS #1
  PRINT #1, " alt   sigma     delta    theta";
  PRINT #1, "  temp   press     dens       a    visc  k.visc"
  PRINT #1, " Kft                           ";
  PRINT #1, "  degR lb/sq.ft  s/cu.ft     fps s/ft-s sq.ft/s"
 
  FOR i = -1 TO 56   ' 5*i is altitude in thousand ft.
    altKm = 5 * i * FT2METERS
    Atmosphere altKm, sigma, delta, theta
    PRINT #1, USING "####"; 5 * i;
    PRINT #1, USING "##.###^^^^"; sigma;
    PRINT #1, USING "##.###^^^^"; delta;
    PRINT #1, USING "##.####"; theta;
    
    temp = KELVIN2RANKINE * TZERO * theta
    pressure = (PZERO / PSF2NSM) * delta
    density = (RHOZERO / SCF2KCM) * sigma
    speedSound = (ASOUNDZERO / FT2METERS) * SQR(theta)

    PRINT #1, USING "####.#"; temp;
    PRINT #1, USING "##.###^^^^"; pressure;
    PRINT #1, USING "##.###^^^^"; density;
    PRINT #1, USING "#####.#"; speedSound;

    viscosity = (1! / PSF2NSM) * MetricViscosity(theta)
    kinematicViscosity = viscosity / density
    PRINT #1, USING "##.###"; 1000000! * viscosity;
    PRINT #1, USING "##.##^^^^"; kinematicViscosity
  NEXT i
 
  CLOSE #1
  PRINT "File e1b.prt has been added to your directory."
END SUB   ' ------------------------------------------ End of Sub LongEnglish

SUB LongMetric
  OPEN "m1b.prt" FOR OUTPUT AS #1
  PRINT #1, " alt    sigma      delta    theta";
  PRINT #1, "  temp    press     dens     a     visc k.visc"
  PRINT #1, "  Km                             ";
  PRINT #1, "  degK   N/sq.m    kg/cu.m m/sec kg/m-s sq.m/s"
  FOR i = -1 TO 43
    altKm = 2 * i
    Atmosphere altKm, sigma, delta, theta
    PRINT #1, USING "####"; 2 * i;
    PRINT #1, USING "##.####^^^^"; sigma;
    PRINT #1, USING "##.####^^^^"; delta;
    PRINT #1, USING "##.####"; theta;

    temp = TZERO * theta
    pressure = PZERO * delta
    density = RHOZERO * sigma
    speedSound = ASOUNDZERO * SQR(theta)

    PRINT #1, USING "####.#"; temp;
    PRINT #1, USING "##.###^^^^"; pressure;
    PRINT #1, USING "##.###^^^^"; density;
    PRINT #1, USING "####.#"; speedSound;

    viscosity = MetricViscosity(theta)
    kinematicViscosity = viscosity / density
    PRINT #1, USING "###.##"; 1000000! * viscosity;
    PRINT #1, USING "##.##^^^^"; kinematicViscosity
  
  NEXT i
 
  CLOSE #1
  PRINT "File m1b.prt has been added to your directory"
END SUB   ' ------------------------------------------- End of Sub LongMetric

FUNCTION MetricViscosity# (theta AS SINGLE)
  t = TZERO * theta
  MetricViscosity = BETAVISC * SQR(t * t * t) / (t + SUTH)
END FUNCTION   ' ---------------------------- End of Function MetricViscosity

SUB ShortEnglish
  OPEN "e2b.prt" FOR OUTPUT AS #1
  PRINT #1, " alt  sigma  delta  theta";
  PRINT #1, "  temp  press    dens     a     visc   k.visc  ratio "
  PRINT #1, " Kft                     ";
  PRINT #1, "  degR   psf   s/cu.ft   fps s/ft-sec sq.ft/s  1/ft"

  FOR i = -1 TO 65   ' i is altitude in thousand ft.
    altKm = i * FT2METERS
    Atmosphere altKm, sigma, delta, theta
    PRINT #1, USING "####"; i;
    PRINT #1, USING "##.####"; sigma;
    PRINT #1, USING "##.####"; delta;
    PRINT #1, USING "##.####"; theta;
    
    temp = (KELVIN2RANKINE * TZERO) * theta
    pressure = (PZERO / PSF2NSM) * delta
    density = (RHOZERO / SCF2KCM) * sigma
    speedSound = (ASOUNDZERO / FT2METERS) * SQR(theta)

    PRINT #1, USING "####.#"; temp;
    PRINT #1, USING "#####.#"; pressure;
    PRINT #1, USING "##.#######"; density;
    PRINT #1, USING "#####.#"; speedSound;
   
    viscosity = (1! / PSF2NSM) * MetricViscosity(theta)
    kinematicViscosity = viscosity / density
    vratio = speedSound / kinematicViscosity
    PRINT #1, USING "##.###"; 1000000! * viscosity;
    PRINT #1, USING "##.##^^^^"; kinematicViscosity;
    PRINT #1, USING "##.##"; .000001 * vratio
  NEXT i
 
  CLOSE #1
  PRINT "File e2b.prt has been added to your directory."
END SUB   ' ----------------------------------------- End of Sub ShortEnglish

SUB ShortMetric
  OPEN "m2b.prt" FOR OUTPUT AS #1
  PRINT #1, " alt  sigma  delta  theta";
  PRINT #1, "  temp  press  dens   a    visc   k.visc ratio"
  PRINT #1, "  Km                     ";
  PRINT #1, "  degK N/sq.m  kcm  m/sec kg/m-s  sq.m/s  1/m"
 
  FOR i = -1 TO 40
    altKm = .5 * i
    SimpleAtmosphere altKm, sigma, delta, theta
    PRINT #1, USING "##.#"; altKm;
    PRINT #1, USING "##.####"; sigma;
    PRINT #1, USING "##.####"; delta;
    PRINT #1, USING "##.####"; theta;
   
    temp = TZERO * theta
    pressure = PZERO * delta
    density = RHOZERO * sigma
    speedSound = ASOUNDZERO * SQR(theta)

    PRINT #1, USING "####.#"; temp;
    ' PRINT #1, USING "##.###^^^^"; pressure;
    PRINT #1, USING "#######"; pressure;
    PRINT #1, USING "##.###"; density;
    PRINT #1, USING "####.#"; speedSound;
     
    viscosity = MetricViscosity(theta)
    kinematicViscosity = viscosity / density
    vratio = speedSound / kinematicViscosity
    PRINT #1, USING "###.##"; 1000000! * viscosity;
    PRINT #1, USING "##.##^^^^"; kinematicViscosity;
    PRINT #1, USING "###.##"; .000001 * vratio
   
  NEXT i
  CLOSE #1
  PRINT "File m2b.prt has been added to your directory"
END SUB   ' ------------------------------------- End of Sub ShortMetricTable

SUB SimpleAtmosphere (alt AS SINGLE, sigma AS SINGLE, delta AS SINGLE, theta AS SINGLE)
'   sigma   density/sea-level standard density
'   delta   pressure/sea-level standard pressure
'   theta   temperature/sea-level standard temperature
'***********************************************************************
'     L O C A L   C O N S T A N T S                                    *
'***********************************************************************
 CONST REARTH! = 6369.1   '  radius of the Earth (km)
 CONST GMR! = 34.1632   '  gas constants
'***********************************************************************
'     L O C A L   V A R I A B L E S                                    *
'***********************************************************************
 DIM h AS SINGLE   '  geopotential altitude
' -----------------------------------------------------------------------
  h = alt * REARTH / (alt + REARTH) ' geometric to geopotential altitude

  IF (h < 11!) THEN
    theta = 1! + (-6.5 / 288.15) * h
    delta = theta ^ (GMR / 6.5)
  ELSE
    theta = 216.65 / 288.15
    delta = .2233611 * EXP(-GMR * (h - 11!) / 216.65)
  END IF

  sigma = delta / theta
END SUB   ' ------------------------------------- End of Sub SimpleAtmosphere

