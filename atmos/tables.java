// PROGRAM Tables;                                    { atmtos/tables.java
// --------------------------------------------------------------------------
// PURPOSE - Make tables of atmospheric properties
// AUTHOR  - Ralph L. Carmichael, Public Domain Aeronautical Software
// REFERENCE - U.S. Standard Atmosphere, 1976. U.S. Govt Printing Office
// REVISION HISTORY
//   DATE  VERS PERSON  STATEMENT OF CHANGES
// 26Feb95  1.0   RLC   Assembled several old codes
//  7Jul95  1.1   RLC   Added viscosity calculations
//  6Aug95  1.2   RLC   Replaced 1962 tables with 1976
// 25Aug95  1.3   RLC   Added MODIFIER
// 19Sep95  1.4   RLC   Added a little precision to pressure table
// 31Oct01  1.5   RLC   Renamed the output files
// --------------------------------------------------------------------------

import java.io.*;
import java.util.*;

class tables {

// P H Y S I C A L   C O N S T A N T S
  static double FT2METERS      = 0.3048;     // mult. ft. to get meters (exact)
  static double KELVIN2RANKINE = 1.8;              // mult kelvins to get deg R
  static double PSF2NSM        = 47.880258;      // mult lb/sq.ft to get N/sq.m
  static double SCF2KCM        = 515.379;    // mult slugs/cu.ft to get kg/cu.m
  static double TZERO          = 288.15;      // sea level temperature, kelvins
  static double PZERO          = 101325.0;        // sea-level pressure, N/sq.m
  static double RHOZERO        = 1.225;           // sea level density, kg/cu.m
  static double AZERO          = 340.294;    // sea-level speed of sound, m/sec

    

// Compute the ratio of temperature to sea-level temperature in the standard
// atmosphere. Correct to 86 km.  Only approximate thereafter.
// =============================================================================
public static double TemperatureRatio(double  alt) { // geometric altitude, km.

  double REARTH=6369.0;    // radius of the Earth (km)
  double GMR = 34.163195;
  int NTAB = 8;
  int i,j,k;

  double htab[] = {0.0,  11.0, 20.0, 32.0, 47.0,
			     51.0, 71.0, 84.852 };
  double ttab[] = { 288.15, 216.65, 216.65, 228.65, 270.65,
			      270.65, 214.65, 186.946 };
  double ptab[] = { 1.0, 2.2336110E-1, 5.4032950E-2, 8.5666784E-3,
     1.0945601E-3, 6.6063531E-4, 3.9046834E-5, 3.68501E-6 };
  double gtab[] = { -6.5, 0, 1.0, 2.8, 0, -2.8, -2.0, 0 };

  double h=alt*REARTH/(alt+REARTH);     //  geometric to geopotential altitude

  i=0; j=NTAB-1;  // starting values for binary search
  for (i=0, j=NTAB-1; j > i+1; )
    {
      k=(i+j)/2;
      if (h < htab[k]) j=k; else i=k;
    } 

  double tgrad=gtab[i];                      // temp. gradient of local layer
  double tbase=ttab[i];                      // base temp. of local layer
  double deltah=h-htab[i];                   // height above local base
  double tlocal=tbase+tgrad*deltah;          // local temperature
  double theta=tlocal/ttab[0];               // temperature ratio
  
  return theta;
}   // ---------------------------------------- End of function TemperatureRatio

// Compute the ratio of pressure to sea-level pressure in the standard
// atmosphere. Correct to 86 km.  Only approximate thereafter.
// =============================================================================
public static double PressureRatio(double  alt) {     // geometric altitude, km.

  double delta;

  double REARTH=6369.0;    // radius of the Earth (km)
  double GMR = 34.163195;
  int NTAB = 8;
  int i,j,k;

  double htab[] = {0.0,  11.0, 20.0, 32.0, 47.0,
			     51.0, 71.0, 84.852 };
  double ttab[] = { 288.15, 216.65, 216.65, 228.65, 270.65,
			      270.65, 214.65, 186.946 };
  double ptab[] = { 1.0, 2.2336110E-1, 5.4032950E-2, 8.5666784E-3,
     1.0945601E-3, 6.6063531E-4, 3.9046834E-5, 3.68501E-6 };
  double gtab[] = { -6.5, 0, 1.0, 2.8, 0, -2.8, -2.0, 0 };

  double h=alt*REARTH/(alt+REARTH);     //  geometric to geopotential altitude

  i=0; j=NTAB-1;  // starting values for binary search
  for (i=0, j=NTAB-1; j > i+1; )
    {
      k=(i+j)/2;
      if (h < htab[k]) j=k; else i=k;
    } 

  double tgrad=gtab[i];                      // temp. gradient of local layer
  double tbase=ttab[i];                      // base temp. of local layer
  double deltah=h-htab[i];                   // height above local base
  double tlocal=tbase+tgrad*deltah;          // local temperature
  double theta=tlocal/ttab[0];                                  // temperature ratio

  if (0.0 == tgrad)                                         // pressure ratio
    delta=ptab[i]*Math.exp(-GMR*deltah/tbase);
  else
    delta=ptab[i]*Math.pow(tbase/tlocal, GMR/tgrad);

  return delta;
}   // ------------------------------------------- End of function PressureRatio

// Compute the ratio of pressure to sea-level pressure in the standard
// atmosphere. Correct to 86 km.  Only approximate thereafter.
// =============================================================================
public static double DensityRatio(double alt) {
  return PressureRatio(alt)/TemperatureRatio(alt);
}   // -------------------------------------------- End of function DensityRatio

// Compute the ratio of temperature to sea-level temperature in the standard
// atmosphere. Correct to 20 km.  Poor approximation thereafter.
// =============================================================================
public static double SimpleTemperatureRatio(double alt) { // geometric altitude, km.

  double theta;                   // density/sea-level standard density
  double REARTH = 6369.0;    // radius of the Earth (km)

  double h=alt*REARTH/(alt+REARTH);     //  geometric to geopotential altitude

  if (h<11.0)  
    theta=(288.15-6.5*h)/288.15; // Troposphere
  else 
    theta=216.65/288.15;    // Stratosphere

  return theta;
}   // ---------------------------------- End of function SimpleTemperatureRatio

// Compute the ratio of pressure to sea-level pressure in the standard
// atmosphere. Correct to 20 km.  Poor approximation thereafter.
// =============================================================================
public static double SimplePressureRatio(double alt) {   // geometric altitude, km.
	double sigma;                   // density/sea-level standard density
	double delta;                 // pressure/sea-level standard pressure
	double theta;           // temperature/sea-level standard temperature

// Compute the temperature,density, and pressure in the standard atmosphere
// Correct to 20 km.  Only approximate thereafter.

  double REARTH = 6369.0;    // radius of the Earth (km)
  double GMR    = 34.163195;   // gas constant


  double h=alt*REARTH/(alt+REARTH);     //  geometric to geopotential altitude

  if (h<11.0)
    {                                                          // Troposphere
      theta=(288.15-6.5*h)/288.15;
      delta=Math.pow(theta, GMR/6.5);
    }
  else
    {                                                         // Stratosphere
      theta=216.65/288.15;
      delta=0.2233611*Math.exp(-GMR*(h-11.0)/216.65);
    }

  sigma=delta/theta;
  return delta;
}   // ------------------------------------- End of function SimplePressureRatio

// Compute the ratio of density to sea-level density in the standard
// atmosphere. Correct to 20 km.  Poor approximation thereafter.
// =============================================================================
public static double SimpleDensityRatio(double alt) { // geometric altitude, km.
  double sigma;                   // density/sea-level standard density
  sigma = SimplePressureRatio(alt)/SimpleTemperatureRatio(alt);
  return sigma;
}   // -------------------------------------- End of function SimpleDensityRatio

// =============================================================================
public static double MetricViscosity(double theta)
{
  double TZERO = 288.15;               // sea level temperature, kelvins
  double BETAVISC = 1.458E-6;     // constant, N-sec/(sq.m-sqrt(kelvins)
  double SUTH = 110.4;                 // Sutherland's constant, kelvins

  double t=theta*TZERO;                              // temperature in kelvins
  return BETAVISC*Math.sqrt(t*t*t)/(t+SUTH);          // viscosity in kg/(m-sec)
}   // -------------------------------------- End of function MetricViscosity



public static void LongUSTable()
{
  double altKm;
  double sigma,delta,theta;
  double temp,pressure,density,asound;
  double viscosity,kinematicViscosity;
  
  
  try {
    FileWriter fw = new FileWriter("us1j.out");
    PrintWriter pw = new PrintWriter(fw,true);
    
    pw.print(" alt   sigma     delta    theta ");
    pw.println(" temp   press     dens       a    visc  k.visc");
    pw.print(" Kft                            ");
    pw.println(" degR lb/sq.ft  s/cu.ft     fps s/ft-s sq.ft/s\n");

    for (int i=-1; i<=56; i++)  {
      altKm=5*i*FT2METERS;  // 5*i is altitude in thousand feet
      theta=TemperatureRatio(altKm);
      delta=PressureRatio(altKm);
      sigma=DensityRatio(altKm);
      pw.printf(" %4d %10.3e %10.3e %10.3e", 5*i, sigma,delta,theta);

      temp=(KELVIN2RANKINE*TZERO)*theta;
      pressure=(PZERO/PSF2NSM)*delta;
      density=(RHOZERO/SCF2KCM)*sigma;
      asound=(AZERO/FT2METERS)*Math.sqrt(theta);      
      pw.printf(" %5.1f %8.3E %9.3E %6.1f", temp,pressure,density,asound);

      viscosity=(1.0/PSF2NSM)*MetricViscosity(theta);
      kinematicViscosity=viscosity/density;
      pw.printf(" %6.3f %8.2E\n", 1.0E6*viscosity, kinematicViscosity);    
    } // end of for-loop
    pw.close();    
  }  // end of try statement
  catch(IOException ioe) {
    System.out.println("Unable to open file for output. Error # "+ioe); 
  }  
     
}   // ---------------------------------------- End of Procedure LongUSTable.



// ==========================================================================
public static void ShortUSTable()
{
  double altKm;
  double sigma,delta,theta;
  double temp,pressure,density, asound;
  double viscosity,kinematicViscosity,vratio;

    
  try {
    FileWriter fw = new FileWriter("us2j.out");
    PrintWriter pw = new PrintWriter(fw,true);
    
    pw.print(" alt  sigma  delta  theta");
    pw.println("  temp  press    dens     a    visc   k.visc  ratio");
    pw.print(" Kft                     ");
    pw.println("  degR   psf   s/cu.ft   fps  s/ft-s  sq.ft/s  1/ft");

  for (int i=-1; i<=65; i++) {
      altKm=i*FT2METERS;  // i is altitude in Kft
      theta=SimpleTemperatureRatio(altKm);
      delta=SimplePressureRatio(altKm);
      sigma=SimpleDensityRatio(altKm);
      pw.printf("%3d %7.4f %7.4f %7.4f ", i,sigma,delta,theta);
      //Itxt.setf(ios::fixed, ios::floatfield);
      //Itxt.width(4); Itxt << i;
      //Itxt.width(7); Itxt.precision(4); Itxt << sigma;
      //Itxt.width(7); Itxt << delta;
      //Itxt.width(7); Itxt << theta;

      temp=(KELVIN2RANKINE*TZERO)*theta;
      pressure=(PZERO/PSF2NSM)*delta;
      density=(RHOZERO/SCF2KCM)*sigma;
      asound=(AZERO/FT2METERS)*Math.sqrt(theta);
      //Itxt.width(6); Itxt.precision(1); Itxt << temp;
      //Itxt.width(7); Itxt.precision(1); Itxt << pressure;
      //Itxt.width(10); Itxt.precision(7); Itxt << density;
      //Itxt.width(7); Itxt.precision(1); Itxt << asound;

      pw.printf("%6.1f %6.1f %9.7f %6.1f", temp,pressure,density,asound);
      viscosity=(1.0/PSF2NSM)*MetricViscosity(theta);
      kinematicViscosity=viscosity/density;
      vratio=asound/kinematicViscosity;
      pw.printf("%6.3f %8.2E %4.2f\n",
	 1.0E6*viscosity, kinematicViscosity, 1.0E-6*vratio);

     

    }
}
  catch(IOException ioe) {
    System.out.println("Unable to open file for output. Error # "+ioe); 
  }  
       
}   // --------------------------------------- End of Procedure ShortUSTable.

// ==========================================================================
public static void LongSITable()
{
  double altKm;
  double sigma,delta,theta;
  double temp,pressure,density,asound;
  double viscosity,kinematicViscosity;


  try {
    FileWriter fw = new FileWriter("si1j.out");
    PrintWriter pw = new PrintWriter(fw,true);
    

    pw.print(" alt    sigma      delta    theta ");
    pw.println(" temp   press     dens      a   visc   k.visc");
    pw.print("  Km                              ");
    pw.println("   K    N/sq.m   kg/cu.m  m/sec kg/m-s sq.m/s");
  //Itxt.setf(ios::showpoint | ios::uppercase | ios::floatfield);

    for(int i=-1; i<=43; i++)
    {
      altKm=2*i;
      theta=TemperatureRatio(altKm);
      delta=PressureRatio(altKm);
      sigma=DensityRatio(altKm);
      pw.printf("%4d %11.4e %11.4e %7.4f ", 2*i, sigma,delta,theta);
      //Itxt.width(4); Itxt << 2*i;
      //Itxt.setf(ios::scientific);
      //Itxt.width(11); Itxt.precision(4); Itxt << sigma;
      //Itxt.width(11); Itxt << delta;
      //Itxt.setf(ios::fixed);
      //Itxt.width(7); Itxt.precision(4); Itxt << theta;

      temp=TZERO*theta;
      pressure=PZERO*delta;
      density=RHOZERO*sigma;
      asound=AZERO*Math.sqrt(theta);
      //Itxt.width(6); Itxt.precision(1); Itxt << temp;
      //Itxt.setf(ios::scientific);
      //Itxt.width(10); Itxt.precision(3); Itxt << pressure;
      //Itxt.width(10); Itxt.precision(3); Itxt << density;
      //Itxt.setf(ios::fixed, ios::floatfield);
      //Itxt.width(6); Itxt.precision(1); Itxt << asound;

      pw.printf("%6.1f %8.3E %8.3E %5.1f", temp,pressure,density,asound);
      viscosity=MetricViscosity(theta);
      kinematicViscosity=viscosity/density;
      pw.printf("%6.2f %8.2E\n", 1.0E6*viscosity, kinematicViscosity);
 
      
    }
  }
  catch(IOException ioe) {
    System.out.println("Unable to open file for output. Error # "+ioe); 
  }  
    
}   // ---------------------------------------- End of Procedure LongSITable.

// ==========================================================================
public static void ShortSITable()
{
  double altKm;
  double sigma,delta,theta;
  double temp,pressure,density,asound;
  double viscosity,kinematicViscosity,vratio;

  try {
    FileWriter fw = new FileWriter("si2j.out");
    PrintWriter pw = new PrintWriter(fw,true);

    pw.print(" alt  sigma  delta  theta");
    pw.println("  temp  press  dens   a    visc   k.visc ratio");
    pw.print("  Km                     ");
    pw.println("    K  N/sq.m  kcm  m/sec kg/m-s  sq.m/s  1/m");
  //Itxt.setf(ios::showpoint | ios::uppercase);

  for(int i=-1; i<=40; i++)
    {
      altKm=0.5*i;
      theta=SimpleTemperatureRatio(altKm);
      delta=SimplePressureRatio(altKm);
      sigma=SimpleDensityRatio(altKm);
      pw.printf("%4.1f %7.4f %7.4f %7.4f ", altKm, sigma,delta,theta);
      //Itxt.setf(ios::fixed, ios::floatfield);
      //Itxt.width(4); Itxt<< setprecision(1); Itxt << altKm;
      //Itxt.width(7); Itxt.precision(4); Itxt << sigma;
      //Itxt.width(7); Itxt << delta;
      //Itxt.width(7); Itxt << theta;

      temp=TZERO*theta;
      pressure=PZERO*delta;
      density=RHOZERO*sigma;
      asound=AZERO*Math.sqrt(theta);
      //Itxt.width(6); Itxt.precision(1); Itxt << temp;
      //Itxt.width(7); Itxt << long(pressure);              // print as integer
      //Itxt.width(6); Itxt.precision(3); Itxt << density;
      //Itxt.setf(ios::fixed, ios::floatfield);
      //Itxt.width(6); Itxt.precision(1); Itxt << asound;

      pw.printf("%6.1f %6.0f %5.3f %5.1f", temp, pressure,density,asound);
      viscosity=MetricViscosity(theta);
      kinematicViscosity=viscosity/density;
      vratio=asound/kinematicViscosity;
      pw.printf("%6.2f %8.2E %5.2f\n",
	 1.0E6*viscosity, kinematicViscosity, 1.0E-6*vratio);
 
      viscosity=MetricViscosity(theta);
      kinematicViscosity=viscosity/density;
      vratio=asound/kinematicViscosity;
      //Itxt.setf(ios::fixed, ios::floatfield);
      //Itxt.width(6); Itxt.precision(2); Itxt << 1.0E6*viscosity;
      //Itxt.setf(ios::scientific, ios::floatfield);
      //Itxt.width(9); Itxt.precision(2); Itxt << kinematicViscosity;
      //Itxt.setf(ios::fixed, ios::floatfield);
      //Itxt.width(6); Itxt.precision(2); Itxt << 1.0E-6*vratio << endl;
   }
  }
  catch(IOException ioe) {
    System.out.println("Unable to open file for output. Error # "+ioe); 
  }  
}   // --------------------------------------- End of Procedure ShortSITable.



// ==========================================================================
  public static void main(String[] args)
{
  String VERSION = "1.0 (16 December 2012)";
  String GREETING =
	  "Tables - A Java program to compute atmosphere tables";
  String AUTHOR =
	  "Ralph L. Carmichael, Public Domain Aeronautical Software";
  String MODIFIER = "";  // put your name here
  String FAREWELL = "Four files added to your directory.";
  String FINALMESS =
	  "Normal termination of tables, version ";

  System.out.println(GREETING);
  System.out.println(AUTHOR);
  if (MODIFIER.length() > 0) System.out.println("Modified by " + MODIFIER);
  System.out.println("version " + VERSION);
  
  LongUSTable();
  ShortUSTable();
  LongSITable();
  ShortSITable();
  System.out.println(FAREWELL);
  System.out.print(FINALMESS);
  System.out.println(VERSION);
}


} // end of class tables ====================================================

