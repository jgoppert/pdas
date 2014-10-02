PROGRAM Tables;
{$APPTYPE CONSOLE}
(* ----------------------------------------------------------------------- *)
(* PURPOSE - Make tables of atmospheric properties                         *)
(* AUTHOR  - Ralph L. Carmichael, Public Domain Aeronautical Software      *)
(* REVISION HISTORY                                                        *)
(*   DATE  VERS PERSON  STATEMENT OF CHANGES                               *)
(* 28Feb95  1.0   RLC   Assembled several old codes                        *)
(*  8Jul95  1.1   RLC   Added viscosity calculations                       *)
(*  5Aug95  1.2   RLC   Slight mod to limits for 1976 atmosphere           *)
(* 25Aug95  1.3   RLC   Assembled into one program for distribution        *)
(* 19Sep95  1.4   RLC   Added a little precision to the pressure tables    *)
(* 29Nov96  1.5   RLC   Error mesage if unable to open output files        *)
(* 28Nov00  1.6   RLC   Slight mods to be a Delphi Console App             *)
(* 31Oct01  1.7   RLC   Renamed the tables                                 *)
(* 05Jul12  1.8   RLC   Commented out SysUtils                             *)
(* ----------------------------------------------------------------------- *)

{ older versions of Delphi may need the following statements:}
USES
  SysUtils;   


  CONST
    VERSION = '1.8 (5 July 2012)';
    GREETING = 'Tables - A Pascal/Delphi program to compute atmosphere tables';
    AUTHOR = 'Ralph L. Carmichael, Public Domain Aeronautical Software';
    MODIFIER = ' ';
    FAREWELL = 'Four files added to your directory.';
    FINALMESS = 'Normal termination of tables, version '+VERSION;

  CONST
    TZERO      = 288.15;     { temperature at sealevel, kelvins }
    PZERO      = 101325.0;   { pressure at sealevel, N/sq.m. (Pa) }
    RHOZERO    = 1.2250;     { density at sealevel, kg/cu.m. }
    ASOUNDZERO = 340.294;    { speed of sound at sealevel, m/sec }

    FT2METERS = 0.3048;                   { mult. ft. to get meters (exact) }
    KELVIN2RANKINE = 1.8;                       { mult kelvins to get deg R }
    PSF2NSM = 47.880258;                      { mult lb/sq.ft to get N/sq.m }
    SCF2KCM = 515.379;                    { mult slugs/cu.ft to get kg/cu.m }

(* ======================================================================= *)
PROCEDURE Welcome;
BEGIN
  WriteLn('Executing ', ParamStr(0));
  WriteLn(GREETING);
  WriteLn(AUTHOR);
  IF MODIFIER <> ' ' THEN WriteLn('Modified by ', MODIFIER);
  WriteLn('     version ', VERSION);
END;  (* --------------------------------------- End of Procedure Welcome. *)

(* ======================================================================= *)
FUNCTION Estr(value : EXTENDED;  trail : WORD) : STRING;
  VAR
    i : WORD;
    logval : EXTENDED;
    expt : INTEGER;
    temp : STRING;
    outStr : STRING;
BEGIN
  IF value=0.0 THEN
    IF trail<=0 THEN
      outStr:='0'
    ELSE
      BEGIN
        outStr:='0.';
        FOR i:=1 TO trail DO outStr:=outStr+'0'
      END
  ELSE
    BEGIN
      logval:=Ln(Abs(value))/Ln(10.0);
      IF logval >=0.0 THEN
        expt:=Trunc(logval)
      ELSE
        expt:=Trunc(logval)-1;
      FOR i:=1 TO Abs(expt) DO
                   IF expt>0 THEN value:=value/10.0 ELSE value:=value*10.0;
      Str(value:1:trail, outStr);
      outStr:=outStr+'E';
      Str(expt:1, temp);
      IF expt >=0 THEN temp:='+' + temp;
      outStr:=outStr+temp;
    END;
  Estr:=outStr
END;  (* ------------------------------------------- End of Function Estr. *)


PROCEDURE Atmosphere(alt : SINGLE;
               VAR sigma : SINGLE;
               VAR delta : SINGLE;
               VAR theta : SINGLE);

CONST
  TABLESIZE = 8;
  REARTH = 6369.0;         { radius of the Earth (km)}
  GMR = 34.163195;         { gas constant }

TYPE
  TABLEINDEX = 1..TABLESIZE;
  ATMOSTABLE = ARRAY[TABLEINDEX] OF SINGLE;

CONST
  htab : ATMOSTABLE = (0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852);
  ttab : ATMOSTABLE = (288.15, 216.65, 216.65, 228.65, 270.65,
                       270.65, 214.65, 186.87 );
  ptab : ATMOSTABLE = (1.0, 2.2336110E-1, 5.4032950E-2, 8.5666784E-3,
                     1.0945601E-3, 6.6063531E-4, 3.9046834E-5, 3.685010E-6 );
  gtab : ATMOSTABLE = (-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0);

VAR
   i,j,k : WORD;
   h, tgrad, deltah, tbase, tlocal : SINGLE;
BEGIN
  h:=alt*REARTH/(alt+REARTH);    { Convert geometric to geopotential altitude }

  i:=1;
  j:=TABLESIZE;
  REPEAT                                   { binary search in ordered table }
    k:=(i+j) DIV 2;
    IF h < htab[k] THEN j:=k ELSE i:=k
  UNTIL j <= i+1;

  tgrad:=gtab[i];
  deltah:=h-htab[i];
  tbase:=ttab[i];

  tlocal:=tbase+tgrad*deltah;
  theta:=tlocal/ttab[1];                                { temperature ratio }

  IF tgrad=0.0 THEN                                        { pressure ratio }
    delta:=ptab[i]*Exp(-GMR*deltah/tbase)
  ELSE
    delta:=ptab[i]*Exp(Ln(tbase/tlocal)*GMR/tgrad);

  sigma:=delta/theta;                                       { density ratio }
END;   (* ----------------------------------- End of Procedure Atmosphere. *)

(* ======================================================================= *)
PROCEDURE SimpleAtmosphere(alt : SINGLE;
                     VAR sigma : SINGLE;
                     VAR delta : SINGLE;
                     VAR theta : SINGLE);
  CONST
    REARTH = 6369.0;         { radius of the Earth (km)}
    GMR = 34.163195;         { gas constant }

  VAR
    h : SINGLE;
BEGIN
  h:=alt*REARTH/(alt+REARTH);          { geometric to geopotential altitude }

  IF h<11.0 THEN
    BEGIN                                                     { Troposphere }
      theta:=(288.15-6.5*h)/288.15;
      delta:=Exp((GMR/6.5)*Ln(theta));               { pow(theta, GMR/6.5); }
    END
  ELSE
    BEGIN                                                    { Stratosphere }
      theta:=216.65/288.15;
      delta:=0.2233611*Exp(-GMR*(h-11.0)/216.65);
    END;

  sigma:=delta/theta;
END;   (* ----------------------------- End of Procedure SimpleAtmosphere. *)

(* ======================================================================= *)
FUNCTION MetricViscosity(const theta : SINGLE) : SINGLE;
  CONST
    TZERO      = 288.15;     { temperature at sealevel, kelvins }
    BETAVISC   = 1.458E-6;   { viscosity term, N sec/(sq.m sqrt(kelvins) }
    SUTH       = 110.4;      { Sutherland's constant, kelvins }
  VAR
    t : SINGLE;                             { temperature in kelvins }
BEGIN                                   { returns viscosity in metric units }
  t:=TZERO*theta;
  MetricViscosity:=BETAVISC*Sqrt(t*t*t)/(t+SUTH)     { Sutherland's formula }
END;   (* -------------------------------- End of Function MetricViscosity *)

(* ======================================================================= *)
PROCEDURE LongUSTable;
  VAR
    Itxt : TEXT;
    i : INTEGER;  { altitude in thousand feet }
    altKm : SINGLE;
    sigma,delta,theta : SINGLE;
    temp,pressure,density,asound : SINGLE;
    viscosity,kinematicViscosity : SINGLE;
BEGIN
{$I-}
  Assign(Itxt, 'us1p.prt'); Rewrite(Itxt);
{$I+}
  IF IOresult <> 0 THEN
    BEGIN
      WriteLn('Unable to open e1p.prt for output'); Halt
    END;
  Write(Itxt, ' alt   sigma    delta   theta');
  Write(Itxt, '  temp   press    dens');
  WriteLn(Itxt, '      a    visc  k.visc' );
  Write(Itxt, ' Kft                         ');
  WriteLn(Itxt, '  degR lb/sq.ft s/cu.ft    fps s/ft-s sq.ft/s');

  FOR i:=-1 TO 56 DO
    BEGIN
      altKm:=5*i*FT2METERS;
      Atmosphere(altKm, sigma,delta,theta);
      Write(Itxt, 5*i:4, ' ', FloatToStrF(sigma,ffExponent,4,1), ' ', Estr(delta,3), theta:7:4);
      temp:=(KELVIN2RANKINE*TZERO)*theta;
      Write(Itxt, temp:6:1);
      pressure:=(PZERO/PSF2NSM)*delta;
      density:=(RHOZERO/SCF2KCM)*sigma;
      asound:=(ASOUNDZERO/FT2METERS)*Sqrt(theta);
      Write(Itxt, ' ', Estr(pressure,3), ' ', Estr(density,3), asound:7:1);
      viscosity:=(1/PSF2NSM)*MetricViscosity(theta);
      kinematicViscosity:=viscosity/density;
      WriteLn(Itxt, 1E6*viscosity:6:3, ' ', Estr(kinematicViscosity,2) );
    END;
  Close(Itxt);
END;  (* ----------------------------------- End of Procedure LongUSTable. *)

(* ======================================================================= *)
PROCEDURE ShortUSTable;
  VAR
    Itxt : TEXT;
    i : INTEGER;  { altitude in thousand feet }
    altKm : SINGLE;
    sigma,delta,theta : SINGLE;
    temp,pressure,density,asound : SINGLE;
    viscosity,kinematicViscosity,vratio : SINGLE;
BEGIN
{$I-}
  Assign(Itxt, 'us2p.prt'); Rewrite(Itxt);
{$I+}
  IF IOresult <> 0 THEN
    BEGIN
      WriteLn('Unable to open e1p.prt for output'); Halt
    END;

  Write(Itxt, ' alt  sigma  delta  theta');
  Write(Itxt, '  temp  press    dens');
  WriteLn(Itxt, '     a     visc   k.visc ratio' );
  Write(Itxt, ' Kft                     ');
  WriteLn(Itxt, '  degR   psf   s/cu.ft   fps s/ft-sec sq.ft/s 1/ft');


  FOR i:=-1 TO 65 DO   { alt=i thousand feet }
    BEGIN
      altKm:=i*FT2METERS;
      SimpleAtmosphere(altKm, sigma,delta,theta);
      Write(Itxt, i:4, sigma:7:4, delta:7:4, theta:7:4);

      temp:=(KELVIN2RANKINE*TZERO)*theta;
      pressure:=(PZERO/PSF2NSM)*delta;
      density:=(RHOZERO/SCF2KCM)*sigma;
      asound:=(ASOUNDZERO/FT2METERS)*Sqrt(theta);
      Write(Itxt, temp:6:1, pressure:7:1, density:10:7, asound:7:1);

      viscosity:=(1/PSF2NSM)*MetricViscosity(theta);
      kinematicViscosity:=viscosity/density;
      vratio:=asound/kinematicViscosity;
      WriteLn(Itxt, 1E6*viscosity:6:3, ' ', Estr(kinematicViscosity,2),
                 1E-6*vratio:6:2);

    END;
  Close(Itxt);
END;  (* ---------------------------------- End of Procedure ShortUSTable. *)

(* ======================================================================= *)
PROCEDURE LongSITable;
  VAR
    Itxt : TEXT;
    i : INTEGER;  { altitude in thousand feet }
    altKm : SINGLE;
    sigma,delta,theta : SINGLE;
    temp,pressure,density, asound : SINGLE;
    viscosity,kinematicViscosity : SINGLE;

BEGIN
{$I-}
  Assign(Itxt, 'si1p.prt'); Rewrite(Itxt);
{$I+}
  IF IOresult <> 0 THEN
    BEGIN
      WriteLn('Unable to open m1p.prt for output'); Halt
    END;

  Write(Itxt, ' alt    sigma     delta   theta');
  Write(Itxt, '  temp   press    dens');
  WriteLn(Itxt, '     a    visc  k.visc' );
  Write(Itxt, '  Km                     ');
  WriteLn(Itxt, '          K   N/sq.m   kg/cu.m m/sec kg/m-s sq.m/s');


  FOR i:=-1 TO 43 DO
    BEGIN
      altKm:=2*i;                              { changes integer to floating }
      Atmosphere(altKm, sigma,delta,theta);
      Write(Itxt, 2*i:4, ' ', Estr(sigma,4), ' ', Estr(delta,4), theta:7:4);
      temp:=TZERO*theta;
      pressure:=PZERO*delta;
      density:=RHOZERO*sigma;
      asound:=ASOUNDZERO*Sqrt(theta);
      Write(Itxt, temp:6:1, ' ', Estr(pressure,3), ' ', Estr(density,3));
      Write(Itxt, asound:6:1);
      viscosity:=MetricViscosity(theta);
      kinematicViscosity:=viscosity/density;
      WriteLn(Itxt, 1E6*viscosity:6:2, ' ', Estr(kinematicViscosity,2));

    END;
  Close(Itxt);
END;  (* ----------------------------------- End of Procedure LongSITable. *)

(* ======================================================================= *)
PROCEDURE ShortSITable;
  VAR
    Itxt : TEXT;
    i : INTEGER;  { steps thru altitude }
    altKm : SINGLE;
    sigma,delta,theta : SINGLE;
    temp,pressure,density, asound : SINGLE;
    viscosity,kinematicViscosity,vratio : SINGLE;
BEGIN
{$I-}
  Assign(Itxt, 'si2p.prt'); Rewrite(Itxt);
{$I+}
  IF IOresult <> 0 THEN
    BEGIN
      WriteLn('Unable to open m2p.prt for output'); Halt
    END;

  Write(Itxt, ' alt  sigma  delta  theta');
  Write(Itxt, '  temp  press  dens');
  WriteLn(Itxt, '   a    visc  k.visc ratio' );
  Write(Itxt, '  Km                     ');
  WriteLn(Itxt, '    K  N/sq.m  kcm  m/sec kg/m-s sq.m/s  1/m');

  FOR i:=-1 TO 40 DO
    BEGIN
      altKm:=0.5*i;
      Atmosphere(altKm, sigma,delta,theta);
      Write(Itxt, altKm:4:1, sigma:7:4, delta:7:4, theta:7:4);

      temp:=TZERO*theta;
      pressure:=PZERO*delta;
      density:=RHOZERO*sigma;
      asound:=ASOUNDZERO*Sqrt(theta);
      Write(Itxt, temp:6:1, pressure:7:0,  density:6:3, asound:6:1);
      viscosity:=MetricViscosity(theta);
      kinematicViscosity:=viscosity/density;
      vratio:=asound/kinematicViscosity;
      WriteLn(Itxt, 1.0E6*viscosity:6:2, ' ',
                    Estr(kinematicViscosity,2),
                    1.0E-6*vratio:6:2);
    END;
  Close(Itxt);
END;  (* ---------------------------------- End of Procedure ShortSITable. *)

(* ======================================================================= *)
BEGIN                                          (* Start of Program Tables. *)
  Welcome;
  LongUSTable;
  ShortUSTable;
  LongSITable;
  ShortSITable;
  WriteLn(FAREWELL);
  WriteLn(FINALMESS);

END.
