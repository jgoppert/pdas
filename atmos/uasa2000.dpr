program uasa2000(input,output);
{Unofficial Australian Standard Atmosphere.
Steven S. Pietrobon.}
{$APPTYPE CONSOLE}
uses
  SysUtils;

const Z0  =       0.; {m}
      date = '3 January 2000';

var P,rho:text;
    Phr,Dh,vs,Z,Zinc,Zmin,Zmax:double;
    ans:char;

const gs =  9.80665; {m/s^2, local gravitational acceleration}
      Re = 6356766.; {m, local radius of planet}

const bmax = 7; {number of geopotential reference levels}
type tarray = array[0..bmax] of double;
const H_:tarray = (0,11000,20000,32000,47000,51000,71000,84852);
        {geopotential height, m'}
      L_:tarray = (-0.0065,0.0,0.0010,0.0028,0.0,-0.0028,-0.0020,0.0);
        {molecular-scale temperature gradiant, K/m'}
      Rs =     8314.32; {Nm/(kmol K), gas constant}
      M0_ =    28.9644; {kg/kmol, mean molecular weight of air}
      g0 =     9.80665; {m/s^2, acceleration of gravity at 45.5425 deg lat.}
var T_,P_:tarray;

const n = 15; {number of regions}
type farray = array[0..n] of double;
     carray = array[0..3,0..n] of double;
const Z_:farray = (86000,100000,115000,130000,150000,175000,200000,250000,
                 300000,400000,500000,600000,700000,800000,900000,1000000);
                 {m, geometric altitude}
      Z7 =      86000.; {m, minimum upper atmosphere geometric height}
      Z12 =   1000000.; {m, maximum upper atmosphere geometric height}
var ln_P,ln_rho:carray;

procedure init_atmosphere;
{Initialse atmosphere parameters}

procedure init_lower;
{Initialise lower atmosphere}
const P0 =      101325.; {Pa, sea level air pressure}
      T0 =       288.15; {K, standard sea-level temperature}
var b:integer;
begin{init lower}
  T_[0] := T0;
  P_[0] := P0;
  for b := 0 to bmax-1 do
    begin{layers}
      T_[b+1] := T_[b] + L_[b]*(H_[b+1]-H_[b]);
      if L_[b] = 0.0
        then P_[b+1] := P_[b]*exp(-(g0*M0_/Rs)*(H_[b+1]-H_[b])/T_[b])
        else P_[b+1] := P_[b]*exp((g0*M0_/Rs)*ln(T_[b]/T_[b+1])/L_[b]);
    end;{layers}
end;{init lower}

procedure init_upper;
{Initialise upper atmosphere}
const P:farray = (3.7338e-1,3.2011e-2,4.0096e-3,1.2505e-3,4.5422e-4,1.7936e-4,
                  8.4736e-5,2.4767e-5,8.7704e-6,1.4518e-6,3.0236e-7,8.2130e-8,
                  3.1908e-8,1.7036e-8,1.0873e-8,7.5138e-9);
                 {Pa, pressure}
      rho:farray =(6.958e-6, 5.604e-7, 4.289e-8, 8.152e-9,2.076e-9,6.339e-10,
                  2.541e-10,6.073e-11,1.916e-11,2.803e-12,5.215e-13,1.137e-13,
                  3.070e-14,1.136e-14,5.759e-15,3.561e-15);
                 {kg/m^3, density}
var I:integer;

procedure csi(x:farray;
              var a:carray);
{Cubic spline interpolation}
var f,h,alpha,l,mu,z:farray;
  i : INTEGER;   { added by RLC }
begin{csi}
  for I := 0 to n-1 do
    h[I] := x[I+1] - x[I];

  for I := 1 to n-1 do
    alpha[I] := 3*(a[0,I+1]*h[I-1]-a[0,I]*(x[I+1]-x[I-1])+a[0,I-1]*h[I])
                 /(h[I-1]*h[I]);
  l[0] := 1;
  mu[0] := 0;
  z[0] := 0;
  for I := 1 to n-1 do
    begin{step 4}
      l[I] := 2*(x[I+1]-x[I-1]) - h[I-1]*mu[I-1];
      mu[I] := h[I]/l[I];
      z[I] := (alpha[I] - h[I-1]*z[I-1])/l[I];
    end;{step 4}

  l[n] := 1;
  z[n] := 0;
  a[2,n] := z[n];
  for I := n-1 downto 0 do
    begin{step 6}
      a[2,I] := z[I] - mu[I]*a[2,I+1];
      a[1,I] := (a[0,I+1]-a[0,I])/h[I] - h[I]*(a[2,I+1]+2*a[2,I])/3;
      a[3,I] := (a[2,I+1]-a[2,I])/(3*h[I]);
    end;{step 6}
end;{csi}

begin{init upper}
  for I := 0 to n do
    begin{ln}
      ln_P[0,I] := ln(P[I]);
      ln_rho[0,I] := ln(rho[I]);
    end;{ln}
  csi(Z_,ln_P);
  csi(Z_,ln_rho);
end;{init upper}

begin{init atmosphere}
  init_lower;
  init_upper;
end;{init atmosphere}

procedure atmosphere(var Zh:double);
{Determine atmosphere pressure, density, and speed of sound.
 Steven S. Pietrobon, 24 Jul 1995. Revised 1 Nov 1999.
 Inputs:  Zh  (m, geometric height above Re)
 Outputs: Phr (Pa, pressure)
          Dh  (kg/m^3, density)
          vs  (m/s, speed of sound)}
const r0 = 6356766.; {m, Earth radius at g0}
      gamma =  1.40; {Ratio of Specific heats for ideal diatomic gas}
var H,Z:double;
    b:integer;

procedure lower_atmosphere;
var Tm:double;
begin{lower atmosphere}
  b := 0;
  while H > H_[b+1] do b := b+1;
  Tm := T_[b] + L_[b]*(H-H_[b]);
  if L_[b] = 0.0
    then Phr := P_[b]*exp(-(g0*M0_/Rs)*(H-H_[b])/T_[b])
    else Phr := P_[b]*exp((g0*M0_/Rs)*ln(T_[b]/Tm)/L_[b]);
  Dh := (Phr/Tm)*(M0_/Rs);
  vs := sqrt((gamma*Rs/M0_)*Tm);
end;{lower atmosphere}

procedure upper_atmosphere;
var dZ:double;
begin{upper atmosphere}
  b := 0;
  while Z > Z_[b+1] do b := b+1;
  dZ := Z - Z_[b];
  Phr := exp(((ln_P[3,b]*dZ+ln_P[2,b])*dZ+ln_P[1,b])*dZ+ln_P[0,b]);
  Dh := exp(((ln_rho[3,b]*dZ+ln_rho[2,b])*dZ+ln_rho[1,b])*dZ+ln_rho[0,b]);
  vs := sqrt((gamma*Rs/M0_)*T_[bmax]);
end;{upper atmosphere}

begin{atmosphere}
  H := (gs/g0)*Re*Zh/(Re+Zh);
  Z := r0*H/(r0-H);
  if Z < Z7
    then lower_atmosphere
    else begin{upper}
           if Z > Z12 then Z := Z12;
           upper_atmosphere;
         end;{upper}
end;{atmosphere}

procedure write_num(var table:text;
                    num:double;
                    width,places:integer);
{Writes reals in scientific form}
var mantissa,expreal:double;
    exponent,wid,dem:longint;
begin{write num}
  if num <> 0
    then begin{non-zero}
           expreal := ln(abs(num))/ln(10);
           if expreal >= 0
             then exponent := trunc(expreal)
             else exponent := trunc(expreal-1);
           mantissa := num*exp(-exponent*ln(10));
           write(table,' ',mantissa:(width+2):(width-1),'e');
           if exponent >= 0
             then write(table,'+')
             else write(table,'-');
           exponent := abs(exponent);
           if places < 1 then places := 1;
           dem := 1;
           for wid := 1 to places-1 do
             dem := 10*dem;
           if dem = 1 then dem := 0;
           wid := places;
           while exponent < dem do
             begin{write zeros}
               write(table,'0');
               dem := dem div 10;
               if dem = 1 then dem := 0;
               wid := wid-1;
             end;{write zeros}
           write(table,exponent:wid);
         end{non-zero}
    else begin{zero}
           write(table,' ',0.0:(width+2):(width-1),'e+');
           for wid := 1 to places do
             write(table,'0');
         end;{zero}
end;{write num}

begin{uasa2000}
  init_atmosphere;
  assign(P,'P.dat');
  rewrite(P);
  assign(rho,'rho.dat');
  rewrite(rho);
  writeln;
  writeln('Unofficial Australian Standard Atmosphere, 2000');
  writeln('Steven S. Pietrobon <steven@sworld.com.au>, ',date);
  writeln('http://www.sworld.com.au/steven/space/atmosphere/');
  writeln;
  repeat
    write('Lower: 0 to 86 km (L), Upper: 86 to 1000 km (U), or Both (B) ');
    write('atmospheres? ');
    readln(ans);
    case ans of
      'l': ans := 'L';
      'u': ans := 'U';
      'b': ans := 'B';
    end;{case}
  until (ans = 'L') or (ans = 'U') or (ans = 'B');
  write('Enter altitude increment in metres: ');
  readln(Zinc);
  writeln('Output: P.DAT (Pressure, Pa), RHO.DAT (Density, kg/m^3)');
  if (ans = 'L') or (ans = 'B') 
    then Zmin := Z0
    else Zmin := Z7;
  if (ans = 'U') or (ans = 'B')
    then Zmax := Z12
    else Zmax := Z7;
  Z := Zmin;
  repeat
    atmosphere(Z);
    write(P,(Z/1000):6:1,' ');
    write_num(P,Phr,6,1);
    writeln(P);
    write(rho,(Z/1000):6:1,' ');
    write_num(rho,Dh,6,1);
    writeln(rho);
    Z := Z+Zinc;
  until Z > Zmax;
  close(P);
  close(rho);
end.{uasa2000}
