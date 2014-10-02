// PROGRAM Tables;                                    { \atmtable\Tables.cpp
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
// 31Jul08  1.6   RLC   Adapted to standard library namespace
// --------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
using namespace std;
//#include <string.h>

  const char* const VERSION = "1.6 (31 Jul 2008)";
  const char* const GREETING =
	  "Tables - A C++ program to compute atmosphere tables";
  const char* const AUTHOR =
	  "Ralph L. Carmichael, Public Domain Aeronautical Software";
//  const char* const MODIFIER = "";
  string MODIFIER ("");
  const char* const FAREWELL = "Four files added to your directory.";
  const char* const FINALMESS =
	  "Normal termination of tables, version ";

// P H Y S I C A L   C O N S T A N T S
  const double FT2METERS      = 0.3048;     // mult. ft. to get meters (exact)
  const double KELVIN2RANKINE = 1.8;              // mult kelvins to get deg R
  const double PSF2NSM        = 47.880258;      // mult lb/sq.ft to get N/sq.m
  const double SCF2KCM        = 515.379;    // mult slugs/cu.ft to get kg/cu.m
  const double TZERO          = 288.15;      // sea level temperature, kelvins
  const double PZERO          = 101325.0;        // sea-level pressure, N/sq.m
  const double RHOZERO        = 1.225;           // sea level density, kg/cu.m
  const double AZERO          = 340.294;    // sea-level speed of sound, m/sec

// F U N C T I O N   P R O T O T Y P E S
  void SimpleAtmosphere(const double, double&, double&, double&);
  void Atmosphere(const double, double&, double&, double&);
  void LongUSTable(void);
  void ShortUSTable(void);
  void LongSITable(void);
  void ShortSITable(void);
  double MetricViscosity(const double);

// ==========================================================================
int main(int argc, char* argv[])
{
  std::cout << "Executing " << argv[0] << '\n';
  std::cout << GREETING << '\n' << AUTHOR << '\n';
//  if (strlen(MODIFIER) > 0) std::cout << "Modified by " << MODIFIER << '\n';
  if (MODIFIER.length() > 0) std::cout << "Modified by " << MODIFIER << '\n';

  std::cout << "     version " << VERSION << std::endl;
  LongUSTable();
  ShortUSTable();
  LongSITable();
  ShortSITable();
  std::cout << FAREWELL << '\n' << FINALMESS << VERSION << std::endl;
  return 0;
}   // ----------------------------------------------- End of Program Tables.

// ==========================================================================
void LongUSTable(void)
{
  double altKm;
  double sigma,delta,theta;
  double temp,pressure,density,asound;
  double viscosity,kinematicViscosity;

  std::ofstream Itxt("us1cpp.prt");
  Itxt << "units of alt are Kft; pressure lbs/sq.ft.; density is slugs/cu.ft.\n";
  Itxt << "speed of sound is ft/second; viscosity is slugs per ft-second;\n"; 
  Itxt << "kinematic viscosity is square feet per second.\n";
  Itxt << " alt     sigma      delta   theta ";
  Itxt << " temp    press      dens      a     visc   k.visc" << std::endl;
//  Itxt << " Kft                              ";
//  Itxt << " degR  lb/sq.ft   sl/cu.ft   fps s/ft-s sq.ft/s\n";
  Itxt.setf(std::ios::showpoint | std::ios::uppercase);

  for(int i=-1; i<=56; i++)  // 5*i is altitude in thousand feet }
    {
      altKm=5*i*FT2METERS;
      Atmosphere(altKm, sigma,delta,theta);
      Itxt.width(4); Itxt << 5*i;
      Itxt.setf(std::ios::scientific, std::ios::floatfield);
      Itxt.width(11); Itxt.precision(3); Itxt << sigma;
      Itxt.width(11); Itxt << delta;
      Itxt.setf(std::ios::fixed, std::ios::floatfield);
      Itxt.width(7); Itxt.precision(4); Itxt << theta;

      temp=(KELVIN2RANKINE*TZERO)*theta;
      pressure=(PZERO/PSF2NSM)*delta;
      density=(RHOZERO/SCF2KCM)*sigma;
      asound=(AZERO/FT2METERS)*sqrt(theta);

      Itxt.width(6); Itxt.precision(1); Itxt << temp;
      Itxt.setf(std::ios::scientific, std::ios::floatfield);
      Itxt.width(11); Itxt.precision(3); Itxt << pressure;
      Itxt.width(11); Itxt.precision(3); Itxt << density;
      Itxt.setf(std::ios::fixed, std::ios::floatfield);
      Itxt.width(7); Itxt.precision(1); Itxt << asound;

      viscosity=(1.0/PSF2NSM)*MetricViscosity(theta);
      kinematicViscosity=viscosity/density;
      Itxt.width(6); Itxt.precision(3); Itxt << 1.0E6*viscosity;
      Itxt.setf(std::ios::scientific, std::ios::floatfield);
      Itxt.width(10); Itxt.precision(2); Itxt << kinematicViscosity << std::endl;
    }
}   // ---------------------------------------- End of Procedure LongUSTable.

// ==========================================================================
void ShortUSTable(void)
{
  double altKm;
  double sigma,delta,theta;
  double temp,pressure,density, asound;
  double viscosity,kinematicViscosity,vratio;

  std::ofstream Itxt("us2cpp.prt");
  Itxt << " alt  sigma  delta  theta";
  Itxt << "  temp  press    dens     a    visc   k.visc  ratio\n";
  Itxt << " Kft                     ";
  Itxt << "  degR   psf   s/cu.ft   fps  s/ft-s  sq.ft/s  1/ft\n";
  Itxt.setf(std::ios::showpoint | std::ios::uppercase);

  for(int i=-1;i<=65;i++)
    {
      altKm=i*FT2METERS;
      SimpleAtmosphere(altKm, sigma,delta,theta);
      Itxt.setf(std::ios::fixed, std::ios::floatfield);
      Itxt.width(4); Itxt << i;
      Itxt.width(7); Itxt.precision(4); Itxt << sigma;
      Itxt.width(7); Itxt << delta;
      Itxt.width(7); Itxt << theta;

      temp=(KELVIN2RANKINE*TZERO)*theta;
      pressure=(PZERO/PSF2NSM)*delta;
      density=(RHOZERO/SCF2KCM)*sigma;
      asound=(AZERO/FT2METERS)*sqrt(theta);
      Itxt.width(6); Itxt.precision(1); Itxt << temp;
      Itxt.width(7); Itxt.precision(1); Itxt << pressure;
      Itxt.width(10); Itxt.precision(7); Itxt << density;
      Itxt.width(7); Itxt.precision(1); Itxt << asound;

      viscosity=(1.0/PSF2NSM)*MetricViscosity(theta);
      kinematicViscosity=viscosity/density;
      vratio=asound/kinematicViscosity;
      Itxt.width(6); Itxt.precision(3); Itxt << 1.0E6*viscosity;
      Itxt.setf(std::ios::scientific, std::ios::floatfield);
      Itxt.width(11); Itxt.precision(2); Itxt << kinematicViscosity;
      Itxt.setf(std::ios::fixed, std::ios::floatfield);
      Itxt.width(5); Itxt.precision(2); Itxt << 1.0E-6*vratio << std::endl;

    }
}   // --------------------------------------- End of Procedure ShortUSTable.

// ==========================================================================
void LongSITable()
{
  double altKm;
  double sigma,delta,theta;
  double temp,pressure,density,asound;
  double viscosity,kinematicViscosity;

  std::ofstream Itxt("si1cpp.prt");
  Itxt << "units of alt are Km; pressure N/sq.m.; density is kg/cu.m.\n";
  Itxt << "speed of sound is m/s; viscosity is kg per m-s;\n";
  Itxt << " alt    sigma       delta     theta ";
  Itxt << " temp   press     dens      a   visc   k.visc" << std::endl;
//  Itxt << "kinematic viscosity is square meters per second.\n";
//  Itxt << "  Km                              ";
//  Itxt << "   K    N/sq.m   kg/cu.m  m/sec kg/m-s sq.m/s\n";
  Itxt.setf(std::ios::showpoint | std::ios::uppercase | std::ios::floatfield);

  for(int i=-1; i<=43; i++)
    {
      altKm=2*i;
      Atmosphere(altKm, sigma,delta,theta);
      Itxt.width(4); Itxt << 2*i;
      Itxt.setf(std::ios::scientific, std::ios::floatfield);
      Itxt.width(12); Itxt.precision(4); Itxt << sigma;
      Itxt.width(12); Itxt << delta;
      Itxt.setf(std::ios::fixed, std::ios::floatfield);
      Itxt.width(7); Itxt.precision(4); Itxt << theta;

      temp=TZERO*theta;
      pressure=PZERO*delta;
      density=RHOZERO*sigma;
      asound=AZERO*sqrt(theta);
      Itxt.width(7); Itxt.precision(1); Itxt << temp;
      Itxt.setf(std::ios::scientific, std::ios::floatfield);
      Itxt.width(10); Itxt.precision(3); Itxt << pressure;
      Itxt.width(10); Itxt.precision(3); Itxt << density;
      Itxt.setf(std::ios::fixed, std::ios::floatfield);
      Itxt.width(6); Itxt.precision(1); Itxt << asound;

      viscosity=MetricViscosity(theta);
      kinematicViscosity=viscosity/density;
      Itxt.width(6); Itxt.precision(2); Itxt << 1.0E6*viscosity;
      Itxt.setf(std::ios::scientific, std::ios::floatfield);
      Itxt.width(10); Itxt.precision(2); Itxt << kinematicViscosity << std::endl;
    }
}   // ---------------------------------------- End of Procedure LongSITable.

// ==========================================================================
void ShortSITable()
{
  double altKm;
  double sigma,delta,theta;
  double temp,pressure,density,asound;
  double viscosity,kinematicViscosity,vratio;

  std::ofstream Itxt("si2cpp.prt");
  Itxt << " alt  sigma  delta  theta";
  Itxt << "  temp  press  dens   a    visc   k.visc ratio\n";
  Itxt << "  Km                     ";
  Itxt << "    K  N/sq.m  kcm  m/sec kg/m-s  sq.m/s  1/m\n";
  Itxt.setf(std::ios::showpoint | std::ios::uppercase);

  for(int i=-1; i<=40; i++)
    {
      altKm=0.5*i;
      SimpleAtmosphere(altKm, sigma,delta,theta);
      Itxt.setf(std::ios::fixed, std::ios::floatfield);
      Itxt.width(4); Itxt.precision(1); Itxt << altKm;
      Itxt.width(7); Itxt.precision(4); Itxt << sigma;
      Itxt.width(7); Itxt << delta;
      Itxt.width(7); Itxt << theta;

      temp=TZERO*theta;
      pressure=PZERO*delta;
      density=RHOZERO*sigma;
      asound=AZERO*sqrt(theta);
      Itxt.width(6); Itxt.precision(1); Itxt << temp;
      Itxt.width(7); Itxt << long(pressure);              // print as integer
      Itxt.width(6); Itxt.precision(3); Itxt << density;
      Itxt.setf(std::ios::fixed, std::ios::floatfield);
      Itxt.width(6); Itxt.precision(1); Itxt << asound;

      viscosity=MetricViscosity(theta);
      kinematicViscosity=viscosity/density;
      vratio=asound/kinematicViscosity;
      Itxt.setf(std::ios::fixed, std::ios::floatfield);
      Itxt.width(6); Itxt.precision(2); Itxt << 1.0E6*viscosity;
      Itxt.setf(std::ios::scientific, std::ios::floatfield);
      Itxt.width(10); Itxt.precision(2); Itxt << kinematicViscosity;
      Itxt.setf(std::ios::fixed, std::ios::floatfield);
      Itxt.width(6); Itxt.precision(2); Itxt << 1.0E-6*vratio << std::endl;
   }
}   // --------------------------------------- End of Procedure ShortSITable.

// ==========================================================================
void Atmosphere(const double  alt,                  // geometric altitude, km.
		double& sigma,           // density/sea-level standard density
		double& delta,         // pressure/sea-level standard pressure
		double& theta)   // temperature/sea-level standard temperature
// Compute the temperature,density, and pressure in the standard atmosphere
// Correct to 86 km.  Only approximate thereafter.
{
  const double REARTH=6369.0;    // radius of the Earth (km)
  const double GMR = 34.163195;
  const int NTAB = 8;
  int i,j,k;

  static double htab[NTAB] = {0.0,  11.0, 20.0, 32.0, 47.0,
			     51.0, 71.0, 84.852 };
  static double ttab[NTAB] = { 288.15, 216.65, 216.65, 228.65, 270.65,
			      270.65, 214.65, 186.946 };
  static double ptab[NTAB] = { 1.0, 2.2336110E-1, 5.4032950E-2, 8.5666784E-3,
     1.0945601E-3, 6.6063531E-4, 3.9046834E-5, 3.68501E-6 };
  static double gtab[NTAB] = { -6.5, 0, 1.0, 2.8, 0, -2.8, -2.0, 0 };

  double h=alt*REARTH/(alt+REARTH);     //  geometric to geopotential altitude

  i=0; j=NTAB-1;  // starting values for binary search
  do
    {
      k=(i+j)/2;
      if (h < htab[k]) j=k; else i=k;
    }  while (j > i+1);

  double tgrad=gtab[i];                      // temp. gradient of local layer
  double tbase=ttab[i];                      // base temp. of local layer
  double deltah=h-htab[i];                   // height above local base
  double tlocal=tbase+tgrad*deltah;          // local temperature
  theta=tlocal/ttab[0];                                  // temperature ratio

  if (0.0 == tgrad)                                         // pressure ratio
    delta=ptab[i]*exp(-GMR*deltah/tbase);
  else
    delta=ptab[i]*pow(tbase/tlocal, GMR/tgrad);

  sigma=delta/theta;                                        //  density ratio
}   // ------------------------------------------- End of function Atmosphere

// ==========================================================================
void SimpleAtmosphere(
	const double alt,                           // geometric altitude, km.
	double& sigma,                   // density/sea-level standard density
	double& delta,                 // pressure/sea-level standard pressure
	double& theta)           // temperature/sea-level standard temperature

// Compute the temperature,density, and pressure in the standard atmosphere
// Correct to 20 km.  Only approximate thereafter.
{
  const double REARTH = 6369.0;    // radius of the Earth (km)
  const double GMR    = 34.163195;   // gas constant


  double h=alt*REARTH/(alt+REARTH);     //  geometric to geopotential altitude

  if (h<11.0)
    {                                                          // Troposphere
      theta=(288.15-6.5*h)/288.15;
      delta=pow(theta, GMR/6.5);
    }
  else
    {                                                         // Stratosphere
      theta=216.65/288.15;
      delta=0.2233611*exp(-GMR*(h-11.0)/216.65);
    }

  sigma=delta/theta;
}   // ------------------------------------- End of function SimpleAtmosphere

// ==========================================================================
double MetricViscosity(const double theta)
{
  const double TZERO = 288.15;               // sea level temperature, kelvins
  const double BETAVISC = 1.458E-6;     // constant, N-sec/(sq.m-sqrt(kelvins)
  const double SUTH = 110.4;                 // Sutherland's constant, kelvins

  double t=theta*TZERO;                              // temperature in kelvins
  return BETAVISC*sqrt(t*t*t)/(t+SUTH);            // viscosity in kg/(m-sec)
}   // -------------------------------------- End of function MetricViscosity
