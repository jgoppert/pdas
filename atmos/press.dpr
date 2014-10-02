PROGRAM Press;                                           { \atmos\press.dpr }
{$APPTYPE CONSOLE}
(* ----------------------------------------------------------------------- *)
(* PURPOSE - Compute the pressure ratio at the altitudes that are the      *)
(*    boundaries between layers of the 1976 Standard Atmosphere            *)
(* AUTHOR  - Ralph L. Carmichael, Public Domain Aeronautical Software      *)
(* NOTES - The type EXTENDED used in these calculations may not be         *)
(*    familiar to everyone. It was included in Turbo Pascal by Borland to  *)
(*    let us use the same 80-bit floating format that is used by the       *)
(*    Intel 8087 processor (and its descendents). The 80-bit numbers will  *)
(*    have even greater accuracy than double precision (IEEE 64 bit).      *)
(* The unit of length used is the kilometer, not meter. The value of the   *)
(*  gas constant R_GAS is given in N km /(kmol K) instead of the value of  *)
(*  8314.32 N m/(kmol K).                                                  *)

(* Running this program should produce the following output:               *)
(* Executing C:press.exe  (will have your full path)                       *)
(* Press - compute pressure at layer boundaries (1976 std)                 *)
(* Ralph L. Carmichael, Public Domain Aeronautical Software                *)
(* Version 1.4 (25Jun12)                                                   *)
(* Hydrostatic constant =      34.163194736310   K/km                      *)
(*                                                                         *)
(*  km      temp           pressure            density                     *)
(*   0   288.15000     1.0000000000000     1.0000000000000                 *)
(*  11   216.65000     0.2233611050922     0.2970759401445                 *)
(*  20   216.65000     0.0540329501078     0.0718651953546                 *)
(*  32   228.65000     0.0085666783593     0.0107959255160                 *)
(*  47   270.65000     0.0010945601338     0.0011653334659                 *)
(*  51   270.65000     0.0006606353133     0.0007033514337                 *)
(*  71   214.65000     0.0000390468337     0.0000524171681                 *)
(*  85   186.94600     0.0000036850095     0.0000056799049                 *)
(* Normal termination of press, version 1.4 (25Jun12)                      *)

(* REVISION HISTORY                                                        *)
(*   DATE  VERS PERSON  STATEMENT OF CHANGES                               *)
(* 10Aug95  1.0   RLC   Original coding                                    *)
(*  3Nov95  1.1   RLC   Added column headers                               *)
(* 21Nov96  1.2   RLC   Writes to a file (easier for Delphi)               *)
(* 28Nov00  1.3   RLC   Made it a Delphi Console Application               *)
(* 25Jun12  1.4   RLC   Commented out the uses SysUtils statement          *)
(* ----------------------------------------------------------------------- *)

(*   if you are running an older version of Delphi or Turbo Pascal,
     you may need the following statement.
uses
  SysUtils;
*)

CONST
    GREETING  = 'Press - compute pressure at layer boundaries (1976 std)';
    AUTHOR    = 'Ralph L. Carmichael, Public Domain Aeronautical Software';
    VERSION   = '1.4 (25Jun12)';
    FINALMESS = 'Normal termination of press, version ' + VERSION;

PROCEDURE ComputePressures(VAR f : TEXT);
  CONST
    GRAVITY : EXTENDED = 9.80665;                       { accel. of gravity }
    MOL_WT  : EXTENDED = 28.9644;                 { molecular weight of air }
    R_GAS   : EXTENDED = 8.31432;             { gas constant N km / (kmol K)}
  VAR
    GMR     : EXTENDED;            { factor used in calc. of pressure, K/km }
    s : EXTENDED;                               { density/sea-level density }
    tgrad : EXTENDED;                              { dT/dh in a given layer }
    t0,t11,t20,t32,t47,t51,t71,t85 : EXTENDED;               { temperatures }
    p0,p11,p20,p32,p47,p51,p71,p85 : EXTENDED; { pressure/sea-level pressure}
                                        { note that the t's are dimensional }
                                          { but the p's are non-dimensional }
BEGIN
  GMR := GRAVITY * MOL_WT / R_GAS;
  WriteLn(f, 'Hydrostatic constant = ', GMR:20:12, '   K/km');
  WriteLn(f);
  WriteLn(f, ' km      temp           pressure            density');

  t0:=288.15;
  p0:=1.0;
  s:=1.0;
  WriteLn(f, '  0', t0:12:5, p0:20:13, s:20:13);

  tgrad:=-6.5;
  t11:=t0+tgrad*11.0;
  p11:=p0*Exp(Ln(t0/t11)*GMR/tgrad);
  s:=p11*t0/t11;
  WriteLn(f, ' 11', t11:12:5, p11:20:13, s:20:13);

  tgrad:=0.0;        { for emphasis. Not used }
  t20:=t11;
  p20:=p11*Exp(-GMR*9.0/t11);                   { special form when tgrad=0 }
  s:=p20*t0/t20;
  WriteLn(f, ' 20', t20:12:5, p20:20:13, s:20:13);

  tgrad:=1.0;
  t32:=t20+tgrad*12.0;
  p32:=p20*Exp(Ln(t20/t32)*GMR/tgrad);
  s:=p32*t0/t32;
  WriteLn(f, ' 32', t32:12:5, p32:20:13, s:20:13);

  tgrad:=2.8;
  t47:=t32+tgrad*15.0;
  p47:=p32*Exp(Ln(t32/t47)*GMR/tgrad);
  s:=p47*t0/t47;
  WriteLn(f, ' 47', t47:12:5, p47:20:13, s:20:13);

  tgrad:=0.0;        { for emphasis. Not used }
  t51:=t47;
  p51:=p47*Exp(-GMR*4.0/t47);                   { special form when tgrad=0 }
  s:=p51*t0/t51;
  WriteLn(f, ' 51', t51:12:5, p51:20:13, s:20:13);

  tgrad:=-2.8;
  t71:=t51+tgrad*20.0;
  p71:=p51*Exp(Ln(t51/t71)*GMR/tgrad);
  s:=p71*t0/t71;
  WriteLn(f, ' 71', t71:12:5, p71:20:13, s:20:13);

  tgrad:=-2.0;
  t85:=t71+tgrad*13.852;
  p85:=p71*Exp(Ln(t71/t85)*GMR/tgrad);
  s:=p85*t0/t85;
  WriteLn(f, ' 85', t85:12:5, p85:20:13, s:20:13);
END;   (* ------------------------------ End of Procedure ComputePressures *)

(* ======================================================================= *)
  VAR
    ff : TEXT;
BEGIN                                           (* Start of Program Press. *)
  Assign(ff, 'press.out');
  Rewrite(ff);
  WriteLn(ff, 'Executing ', ParamStr(0));
  WriteLn(ff, GREETING);
  WriteLn(ff, AUTHOR);
  WriteLn(ff, 'Version ', VERSION);

  ComputePressures(ff);
  WriteLn(ff, FINALMESS);
  Close(ff);
  WriteLn('The file press.out has been added to your directory.')
END.
