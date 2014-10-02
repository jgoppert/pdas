program ussa1976(input,output);
{ Purpose - Determines U.S. Standard Atmosphere 1976 tables. }
{ Steven S. Pietrobon, 30 Dec 1999.}
{ Added to PDAS v6 1 January 2001,  Ralph Carmichael }
{ Variable As changed to Astar to avoid conflict with Delphi keyword }
{$APPTYPE CONSOLE}
uses
  SysUtils;

type gases = (N2,O,O2,Ar,He,H);
     gasarray = array[gases] of double;

const k_ = 1.380622e-23; {Nm/K, Boltzmann constant}
      Na =  6.022169e26; {1/kmol, Avogadro constant}
      Rs =      8314.32; {Nm/(kmol K), gas constant}
      M0 =      28.9644; {kg/kmol, mean molecular weight of air}
      g0 =      9.80665; {m/s^2, acceleration of gravity at 45.5425 deg lat.}
      r0 =     6356766.; {m, Earth radius at g0}
      P0 =      101325.; {Pa, air pressure at g0}
      T0 =       288.15; {K, standard sea-level temperature}
      Td =       273.15; {K, 0 degrees C}
      M_:gasarray = (28.0134,15.9994,31.9988,39.948,4.0026,1.00797);
        {molecular weight, kg/kmol}
      linespace =    10; {number of lines}

var table1,table2,table3,table4{,table}:text;
    rho0,pi,Zinc:double;
    line:longint;
    ans:char;

function gravity(Z:double):double;
{Input:  Z        geometric height, m
Output:  gravity  acceleration of gravity, m/s^2}
begin{gravity}
  gravity := g0*sqr(r0/(r0+Z));
end;{gravity}

procedure write_num(var table:text;
                    num:double;
                    width:integer);
{Writes reals in scientific form}
var mantissa,expreal:double;
    exponent:integer;
begin{write num}
  if num <> 0
    then begin{positive}
           expreal := ln(abs(num))/ln(10);
           if expreal >= 0
             then exponent := trunc(expreal)
             else exponent := trunc(expreal-1);
           mantissa := num*exp(-exponent*ln(10));
           write(table,' ',mantissa:6:4,'e');
           if exponent >= 0
             then write(table,'+')
             else write(table,'-');
           if (abs(exponent) < 10) and (width > 1) then write(table,'0');
           write(table,abs(exponent):1);
         end{positive}
    else write(table,' 0.0000e+00');
end;{write num}

procedure write_real(var table:text;
                     num:double;
                     width,places:integer);
{write(table,num:width:places) doesn't work properly in my version of TMT 
pascal}
var tmp,rem,dem,wid:longint;
begin{write real}
  dem := 1;
  for wid := 1 to places do
    dem := 10*dem;
  tmp := round(dem*abs(num));
  wid := width-places-1;
  if wid < 0 then wid := 0;
  if num >= 0
    then write(table,' ',(tmp div dem):wid)
    else write(table,' ',(-tmp div dem):wid);
  if wid > 0 then
    begin{decimal}
      write(table,'.');
      rem := tmp mod dem;
      tmp := dem div 10;
      if tmp = 1 then tmp := 0;
      wid := places;
      while rem < tmp do
        begin{write zeros}
          write(table,'0');
          tmp := tmp div 10;
          if tmp = 1 then tmp := 0;
          wid := wid-1;
        end;{write zeros}
      write(table,rem:wid);
    end;{decimal}
end;{write real}

procedure write_tables(Z,H,T,Tm,P,rho,g:double);
const sigma = 3.65e-10; {m, effective collision diameter}
var Hp,n,V,L,M:double;
begin{write tables}
  write(table1,round(Z):7);      {geometric height, m}
  write(table1,round(H):7);      {geopotential height, m}
  write_real(table1,T,7,2);      {kinetic temperature, K}
  write_real(table1,T-Td,6,2);   {kinetic temperature, C}
  write_real(table1,Tm,7,2);     {molecular-scale temperature, K}
  write_num(table1,P,1);         {pressure, Pa}
  write_num(table1,P/P0,2);      {normalised pressure}
  write_num(table1,rho,2);       {mass density, kg/m^3}
  write_num(table1,rho/rho0,2);  {normalised mass density}
  writeln(table1);

  write(table2,round(Z):7);      {geometric height, m}
  write(table2,round(H):7);      {geopotential height, m}
  write(table2,'  ',g:6:4);      {acceleration due to gravity, m/s^2}
  Hp := (Rs/M0)*(Tm/g);
  write(table2,'  ',round(Hp):6,' '); {pressure scale height, m}
  n := (Na/Rs)*(P/T);
  write_num(table2,n,2);         {number density, 1/m^3}
  V := sqrt((8*Rs/(pi*M0))*Tm);
  write(table2,' ',V:6:1);       {mean air-particle speed, m/s}
  L := 1.0/((sqrt(2)*pi*sqr(sigma))*n);
  write_num(table2,V/L,2);       {collision frequency, 1/s}
  write_num(table2,L,1);         {mean free path, m}
  M := T*M0/Tm;
  write(table2,'  ',M:6:3);      {mean molecular weight, kg/kmol}
  writeln(table2);

  line := line+1;
  if line = linespace then
    begin{line space}
      writeln(table1);
      writeln(table2);
      line := 0;
    end;{line space}
end;{write tables}

procedure lower_atmosphere;
const bmax = 7; {number of geopotential reference levels}
type temparray = array[0..bmax] of double;
const H_:temparray = (0,11000,20000,32000,47000,51000,71000,84852);
        {geopotential height, m'}
      Lm:temparray = (-0.0065,0.0,0.0010,0.0028,0.0,-0.0028,-0.0020,0.0);
        {molecular-scale temperature gradiant, K/m'}
      Zs = -5000; {m, starting altitude}
var Astar,H,Z,P,Tm,T,M,rho,g:double;
    b,Zo,Zi,Zt:integer;
    P_,T_:temparray;

function mol(Z:double):double;
{mean molecular weight M from 80 to 86 km}
const f:array[0..12] of double = 
            (1.000000,0.999996,0.999989,0.999971,0.999941,0.999909,
             0.999870,0.999829,0.999786,0.999741,0.999694,0.999641,0.999579);
      Zinc = 500.; {m, incremental height}
      Zm = 80000.; {m, initial altitude}
var I:integer;
    Zi:double;
begin{mol}
  if Z < Zm
    then mol := M0
    else begin{change}
           Zi := (Z-Zm)/Zinc;
           I := trunc(Zi);
           mol := M0*((f[I+1]-f[I])*(Zi-I) + f[I]);
         end;{change}
end;{mol}

procedure write_table3(Z,H,T,Tm,rho:double);
const gamma =     1.40; {Ratio of Specific heats for ideal diatomic gas}
      beta  = 1.458e-6; {kg/(s m K^1/2)}
      S     =    110.4; {Sutherland's constant}
var Cs,mu,eta,kt:double;
begin{write table 3}
  write(table3,round(Z):7);      {geometric height, m}
  write(table3,round(H):7);      {geopotential height, m}
  Cs := sqrt((gamma*Rs/M0)*Tm);
  write(table3,' ',Cs:6:2);      {speed of sound, m/s}
  mu := beta*sqrt(T)/(1+S/T);
  write_num(table3,mu,1);        {dynamic viscosity, Ns/m^2}
  eta := mu/rho;
  write_num(table3,eta,1);       {kinematic viscosity, m^2/s}
  kt := 2.64638e-3*sqrt(T)/(1+245.4*exp(-12*ln(10)/T)/T);
  write_num(table3,kt,1);   {coef. of thermal conductivity, W/(m K)}
  writeln(table3);
  if line = linespace-1 then writeln(table3);
end;{write table3}

begin{lower atmosphere}
  assign(table3,'table3.dat');
  rewrite(table3);
  writeln(table3,'              |Speed |         |         |Coef. of');
  writeln(table3,'   Altitude   |  of  | Dynamic |Kinematic|thermal');
  writeln(table3,'              |sound |viscosity|viscosity|conduct.');
  writeln(table3,'-------+------+------+---------+---------+---------');
  writeln(table3,'   Z   |   H  | C_s  |   mu    |   eta   |   k_t');
  writeln(table3,'  (m)  |  (m) |(m/s) |(Ns/m^2) | (m^2/s) |(W/m 1/K)');
  writeln(table3,'-------+------+------+---------+---------+---------');
  T_[0] := T0;
  P_[0] := P0;
  Astar := g0*M0/Rs;
  Zo := round(Zs/Zinc);
  for b := 0 to bmax-1 do
    begin{layers}
      T_[b+1] := T_[b] + Lm[b]*(H_[b+1]-H_[b]);
      if Lm[b] = 0.0
        then P_[b+1] := P_[b]*exp(-Astar*(H_[b+1]-H_[b])/T_[b])
        else P_[b+1] := P_[b]*exp(Astar*ln(T_[b]/T_[b+1])/Lm[b]);
      Zt := trunc(r0*H_[b+1]/((r0-H_[b+1])*Zinc));
      for Zi := Zo to Zt do
        begin{write table}
          Z := Zi*Zinc;
          H := r0*Z/(r0+Z);
          Tm := T_[b] + Lm[b]*(H-H_[b]);
          if Lm[b] = 0.0
            then P := P_[b]*exp(-Astar*(H-H_[b])/T_[b])
            else P := P_[b]*exp(Astar*ln(T_[b]/Tm)/Lm[b]);
          M := mol(Z);
          T := Tm*M/M0;
          rho := (P/Tm)*(M0/Rs);
          g := gravity(Z);
          write_table3(Z,H,T,Tm,rho);
          write_tables(Z,H,T,Tm,P,rho,g);
        end;{write table}
      Zo := Zt+1;
    end;{layers}
  close(table3);
end;{lower atmosphere}

procedure upper_atmosphere;
const Z7  =    86000; {m, minimum geometric height}
      Zm  =   100000; {m, change of N2 layer}
      Zk  =   115000; {m, end of eddy layer}
      Zh  =   150000; {m, start of H layer}
      Z11 =   500000; {m, end of itegration for H layer}
      Z12 =  1000000; {m, maximum geometric height}
      T7  = 186.8673; {K, temperature at Z7}
      T11 = 999.2356; {K, kinetic temperature at Z11}
      n_:gasarray = (1.129794e20,8.6e16,3.030898e19,1.3514e18,7.58173e14,
         8.0e10); {number density for N2..He at 86 km and H at 500 km, 1/m^3}
      alpha:gasarray = (0.0,0.0,0.0,0.0,-0.40,-0.25);
        {thermal-diffusion coefficient, dimensionless}
      phi =      7.2e11; {1/m^2 1/s, vertical flux of H}
      taui = 0.40463343; {dimensionless, tau integral}
      Hi = 9.8776951e10; {1/m^3, H integral}
      Zint =         10; {m, integration region}
      reg =           2; {number of integral regions}
      epsilon =     0.1; {a small number}
var n,fint:gasarray;
    gas:gases;
    Z,T,dT_dZ,g,K:array[1..reg] of double;
    f:array[0..reg] of gasarray;
    I:integer;
    tau:array[0..1] of double;
    tauint:double;

function eddy(Z:double):double;
{Input: Z     geometric height, m
Output: eddy  eddy-diffusion coefficient, m^2/s}
const K7  =  120.0; {m^2/s, eddy diffusion coefficient at 86 km}
      K10 =    0.0; {m^2/s, eddy diffusion coefficient at 115 km}
      Ze  =  95000; {m, geometric height}
      A   = 4.00e8; {m^2}
begin{eddy}
  if Z < Ze
    then eddy := K7
    else if Z < Zk-epsilon
           then eddy := K7*exp(1-A/(A-sqr(Z-Ze)))
           else eddy := K10;
end;{eddy}

procedure temperature(Z:double;
                      var T,dT_dZ:double);
{Input: Z:double      geometric height, m
Output: T:double      kinetic temperature, K
        dT_dZ:double  kinetic temperature gradiant, K/m}
const T9  =     240.0; {K, temperature at 110 km}
      T10 =     360.0; {K, temperature ar 120 km}
      Tinf =   1000.0; {K, defined temperature}
      Z8  =     91000; {m, geometric height}
      Z9  =    110000;
      Z10 =    120000;
      Tc  =  263.1905; {K}
      A   =  -76.3232; {K}
      a_  =   19942.9; {m}
      Lk9 =     0.012; {K/m, kinetic temperature gradient}
      lambda = 1.875e-5; {1/m}
var B,b_,eta:double;
begin{temperature}
  if Z < Z8-epsilon
    then begin
           T := T7;
           dT_dZ := 0.0;
         end
    else if Z < Z9-epsilon
           then begin
                  b_ := (Z-Z8)/a_;
                  B := sqrt(1-sqr(b_));
                  T := Tc + A*B;
                  dT_dZ := -A*b_/(a_*B);
                end
           else if Z < Z10-epsilon
                  then begin
                         T := T9 + Lk9*(Z-Z9);
                         dT_dZ := Lk9;
                       end
                  else begin
                         b_ := (r0+Z10)/(r0+Z);
                         B := (Tinf-T10)*exp(-lambda*(Z-Z10)*b_);
                         T := Tinf - B;
                         dT_dZ := lambda*sqr(b_)*B;
                       end;
end;{temperature}

function flux(Z:double;
              gas:gases):double;
{Input: Z     geometric height, m
        gas   current gas
Output: flux  nu_i/(D_i+K), 1/m}
type fluxarray = array[O..He] of double;
const qO = -3.416248e-12; {1/m^3}
      uO = 97000;         {m}
      wO = 5.008765e-13;  {1/m^3}
      Q:fluxarray = (-5.809644e-13,1.366212e-13,9.434079e-14,-2.457369e-13);
        {coefficient, 1/m^3}
      U:fluxarray = (56903.11,86000,86000,86000);
        {coefficient, m}
      W:fluxarray = (2.706240e-14,8.333333e-14,8.333333e-14,6.666667e-13);
        {coefficient, 1/m^3}
var X:double;
    A,a_:double;
begin{flux}
  a_ := Z-U[gas];
  A := sqr(a_);
  X := Q[gas]*A*exp(-W[gas]*A*a_);
  if (Z < uO) and (gas = O) then
    begin{O}
      a_ := uO-Z;
      A := sqr(a_);
      X := X + qO*A*exp(-wO*A*a_);
    end;{O}
  flux := X;
end;{flux}

procedure weight(Z,T:double;
                 gas:gases;
                 n:gasarray;
                 var M,D:double);
{Input: Z    geometric height, m
        T    kinetic temperature, K
        gas  current gas
        n    number density, 1/m^3
Output: M    mean molecular weight, kg/kmol
        D    molecular-diffusion coefficient, m^2/s}
type diffarray = array[O..H] of double;
const a:diffarray = (6.986e20,4.863e20,4.487e20,1.700e21,3.305e21);
        {coefficient, 1/m 1/s}
      b:diffarray = (0.750,0.750,0.870,0.691,0.500);
        {coefficient, dimensionless}
var maxgas,igas:gases;
    sumn:double;
begin{weight}
  case gas of
    O,O2: maxgas := N2;
    Ar,He: maxgas := O2;
    H: maxgas := He;
  end;{case gas}
  sumn := 0;
  for igas := N2 to maxgas do
    sumn := sumn + n[igas];
  if Z < Zm-epsilon
    then M := M0
    else begin{weight}
           M := 0;
           for igas := N2 to maxgas do
             M := M + n[igas]*M_[igas];
           M := M/sumn;
         end;{weight}
  D := a[gas]*exp(b[gas]*ln(T/Td))/sumn;
end;{weight}

function fcalc(Z,g,K,T,dT_dZ:double;
               gas:gases;
               n:gasarray):double;
var M,D,f:double;
begin{fcalc}
  case gas of
    N2: begin{nitrogen}
          if Z < Zm-epsilon
             then M := M0
             else M := M_[N2];
          fcalc := M*g/(Rs*T);
         end;{nitrogen}
    O..He: begin{Oxygen to Helium}
             if Z < Zk-epsilon
               then begin{diffusion}
                      weight(Z,T,gas,n,M,D);
                      f := D*(M_[gas]+M*K/D+alpha[gas]*Rs*dT_dZ/g)/(D+K);
                    end{diffusion}
               else f := M_[gas]+alpha[gas]*Rs*dT_dZ/g;
             fcalc := g*f/(Rs*T) + flux(Z,gas);
           end;{Oxygen to Helium}
     H: if Z < Zh-epsilon
           then fcalc := 0
           else if Z < Zh+epsilon
                  then begin{init tau}
                         tau[1] := (g/T)*(M_[H]/Rs);
                         tauint := 0;
                         fcalc := 0;
                       end{init tau}
                  else begin{integrate tau}
                         tau[0] := tau[1];
                         tau[1] := (g/T)*(M_[H]/Rs);
                         tauint := tauint + Zint*(tau[0]+tau[1])/(2*reg);
                         if Z < Z11-epsilon
                           then begin{fcalc}
                                  weight(Z,T,H,n,M,D);
                                  fcalc := phi*exp((1+alpha[H])*ln(T/T11) +
                                           tauint - taui)/D;
                                end{fcalc}
                           else fcalc := 0;
                       end;{integrate tau}
   end;{case}
end;{fcalc}

procedure write_data(Z,T,g:double;
                     n:gasarray);
var gas:gases;
    ns,nms,P,rho,H_,Tm:double;
begin{write upper}
  ns := 0;
  nms := 0;
  for gas := N2 to H do
    begin{gas cycle}
      ns := ns + n[gas];
      nms := nms + n[gas]*M_[gas];
    end;{gas cycle}
  H_ := r0*Z/(r0+Z);
  Tm := T*M0*ns/nms;
  P := ns*k_*T;
  rho := nms/Na;
  write(table4,round(Z):7);      {geometric height, m}
  write(table4,round(H_):7);     {geopotential height, m}
  for gas := N2 to H do
    write_num(table4,n[gas],2);
  writeln(table4);
  if line = linespace-1 then writeln(table4);
  write_tables(Z,H_,T,Tm,P,rho,g);

{  write(table,(Z/1000):6:1);
  write_num(table,rho,2);
  writeln(table);}
end;{write lower}

begin{upper atmosphere}
{  assign(table,'table.dat');
  rewrite(table);}
  assign(table4,'table4.dat');
  rewrite(table4);
  write(table4,'   Altitude   |                     Number density (1/m^3)');
  writeln(table4);
  write(table4,'-------+------+----------+----------+----------+----------');
  writeln(table4,'+----------+----------');
  write(table4,' Z (m) | H (m)|    N2    |    O     |    O2    |    Ar    ');
  writeln(table4,'|   He     |    H');
  write(table4,'-------+------+----------+----------+----------+----------');
  writeln(table4,'+----------+----------');
  for I := 1 to reg do
    Z[I] := Z7 - (reg-I)*Zint/reg;
  g[reg] := gravity(Z7);
  K[reg] := eddy(Z7);
  T[reg] := T7;
  dT_dZ[reg] := 0;
  for gas := N2 to He do
    n[gas] := n_[gas];
  n[H] := 0;
  for gas := N2 to H do
    begin{init}
      f[reg][gas] := fcalc(Z[reg],g[reg],K[reg],T[reg],dT_dZ[reg],gas,n);
      fint[gas] := 0;
    end;{init}
  write_data(Z[reg],T[reg],g[reg],n);
  repeat
    for I := 1 to reg do
      begin{calc var}
        Z[I] := Z[I] + Zint;
        g[I] := gravity(Z[I]);
        K[I] := eddy(Z[I]);
        temperature(Z[I],T[I],dT_dZ[I]);
      end;{calc var}
    f[0] := f[reg];
    for gas := N2 to H do
      begin{gas cycle}
        for I := 1 to reg do
          f[I][gas] := fcalc(Z[I],g[I],K[I],T[I],dT_dZ[I],gas,n);
        fint[gas] := fint[gas] 
                     + Zint*(f[0][gas]+4*f[1][gas]+f[2][gas])/(reg*3);
        case gas of
          N2..He: n[gas] := n_[gas]*T7*exp(-fint[gas])/T[reg];
          H: if Z[reg] < Zh-epsilon
               then n[H] := 0
               else n[H] := (n_[H]+Hi-fint[H])*exp((1+alpha[H])*ln(T11/T[reg])
                             + taui - tauint);
        end;{case}
      end;{gas cycle}
    if (round(Z[reg]) mod round(Zinc)) = 0
      then write_data(Z[reg],T[reg],g[reg],n);
  until Z[reg]+epsilon > Z12;
  close(table4);
{  close(table);}
end;{upper atmosphere}

begin{ussa1976}
  assign(table1,'table1.dat');
  rewrite(table1);
  write(table1,'   Altitude   |     Temperature      |     Pressure       ');
  writeln(table1,'|      Density');
  write(table1,'-------+------+-------+------+-------+---------+----------');
  writeln(table1,'+----------+----------');
  write(table1,'  Z (m)| H (m)| T (K) | t (C)|T_M (K)| P (Pa)  |  P/P_0   ');
  writeln(table1,'rho (kg/m^3) rho/rho_0');
  write(table1,'-------+------+-------+------+-------+---------+----------');
  writeln(table1,'+----------+----------');

  assign(table2,'table2.dat');
  rewrite(table2);
  write(table2,'              |Accel. |Pressure|          |Part- ');
  writeln(table2,'|          |  Mean   |');
  write(table2,'   Altitude   |due to | scale  | Number   | icle ');
  writeln(table2,'|Collision |  free   |Molecular');
  write(table2,'              |gravity| height | density  |speed ');
  writeln(table2,'|frequency |  path   | weight');
  write(table2,'-------+------+-------+--------+----------+------');
  writeln(table2,'+----------+---------+---------');
  write(table2,'   Z   |   H  |   g   |  H_P   |    n     |  V   ');
  writeln(table2,'|    nu    |    L    |    M');
  write(table2,'  (m)  |  (m) |(m/s^2)|  (m)   | (1/m^3)  |(m/s) ');
  writeln(table2,'|   (1/s)  |   (m)   |(kg/kmol)');
  write(table2,'-------+------+-------+--------+----------+------');
  writeln(table2,'+----------+---------+---------');

  rho0 := P0*M0/(Rs*T0);
  pi := 4*arctan(1.0);
  line := 0;

  writeln;
  writeln('U.S. Standard Atmosphere, 1976');
  write('Unofficial implementation by Steven S. Pietrobon ');
  writeln('<steven@sworld.com.au>');
  write('30 December 1999      ');
  writeln('http://www.sworld.com.au/steven/space/atmosphere/');
  writeln;
  repeat
    write('Lower: -5 to 86 km (L), Upper: 86 to 1000 km (U), or Both (B) ');
    write('atmospheres? ');
    readln(ans);
    case ans of
      'l': ans := 'L';
      'u': ans := 'U';
      'b': ans := 'B';
    end;{case}
  until (ans = 'L') or (ans = 'U') or (ans = 'B');
  write('Enter altitude increment in metres');
  if (ans = 'L')
    then write(': ')
    else write(' (must be multiple of 10 m): ');
  readln(Zinc);
  writeln('Output:');
  writeln('TABLE1.DAT - Z, H, T, t, T_M, P, P/P_0, rho, rho/rho_0');
  writeln('TABLE2.DAT - Z, H, g, H_P, n, V, nu, L, M');
  writeln('TABLE3.DAT - Z, H, C_s, mu, eta, k_t (lower atmosphere only)');
  write('TABLE4.DAT - Z, H, n of N2, O, O2, Ar, He, H (upper atmosphere ');
  writeln('only)');

  if (ans = 'L') or (ans = 'B') then lower_atmosphere;
  if (ans = 'U') or (ans = 'B') then upper_atmosphere;

  close(table1);
  close(table2);
end.
