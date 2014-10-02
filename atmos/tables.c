/* PROGRAM Tables;                                      \atmtable\tables.c */
/* ----------------------------------------------------------------------- */
/* PURPOSE - Make tables of atmospheric properties                         */
/* AUTHOR  - Ralph L. Carmichael, Public Domain Aeronautical Software      */
/* REVISION HISTORY                                                        */
/*   DATE  VERS PERSON  STATEMENT OF CHANGES                               */
/*  6Mar95  1.0   RLC   Assembled several old codes                        */
/*  8Jul95  1.1   RLC   Put lookup formula in-line                         */
/*  6Aug95  1.2   RLC   Replaced 1962 tables with 1976 tables              */
/* 19Sep95  1.3   RLC   Added a little precision to pressure tables        */
/* 31Oct01  1.4   RLC   Renamed tables                                     */
/* ----------------------------------------------------------------------- */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define VERSION   "1.4 (31Oct2001)"
#define GREETING  "Tables - A C program to compute atmosphere tables"
#define AUTHOR    "Ralph L. Carmichael, Public Domain Aeronautical Software"
#define MODIFIER  ""
#define FAREWELL  "Four files added to your directory."
#define FINALMESS "Normal termination of tables."

/*    P H Y S I C A L   C O N S T A N T S    */
#define FT2METERS  0.3048               /* mult. ft. to get meters (exact) */
#define KELVIN2RANKINE  1.8                     /* mult deg K to get deg R */
#define PSF2NSM  47.880258                  /* mult lb/sq.ft to get N/sq.m */
#define SCF2KCM  515.379                /* mult slugs/cu.ft to get kg/cu.m */
#define TZERO    288.15                  /* sea-level temperature, kelvins */
#define PZERO    101325.0                    /* sea-level pressure, N/sq.m */
#define RHOZERO  1.225                       /* sea-level density, kg/cu.m */
#define AZERO    340.294                  /* speed of sound at S.L.  m/sec */
#define BETAVISC 1.458E-6                            /* viscosity constant */
#define SUTH     110.4                   /* Sutherland's constant, kelvins */

/*    F U N C T I O N   P R O T O T Y P E S    */
  void SimpleAtmosphere(const float, float*, float*, float*);
  void Atmosphere(const float, float*, float*, float*);
  void LongUSTable(void);
  void ShortUSTable(void);
  void LongSITable(void);
  void ShortSITable(void);
  float MetricViscosity(const float);


/* ======================================================================= */
int main(int argc, char *argv[])
{
  printf("Executing %s\n", argv[0]);
  puts(GREETING);
  puts(AUTHOR);
  if (strlen(MODIFIER) > 0) printf("Modified by %s\n", MODIFIER);
  printf("     version %s\n", VERSION);
  LongUSTable();
  ShortUSTable();
  LongSITable();
  ShortSITable();

  puts(FAREWELL);
  puts(FINALMESS);
  return 0;
}   /* ----------------------------------------------- End of Program Tables. */

/* ======================================================================= */
void LongUSTable(void)
{

  FILE *Itxt;
  int  i;  /* altitude in thousand feet */
  float altKm;
  float sigma,delta,theta;
  float temp,pressure,density,asound;
  float viscosity, kinematicViscosity;

  Itxt=fopen("us1c.prt", "w");
  fprintf(Itxt, " alt   sigma     delta    theta ");
  fprintf(Itxt, " temp   press     dens       a    visc  k.visc\n" );
  fprintf(Itxt, " Kft                            ");
  fprintf(Itxt, " degR lb/sq.ft  s/cu.ft     fps s/ft-s sq.ft/s\n");

  for(i=-1; i<=56; i++)
    {
      altKm=5*i*FT2METERS;
      Atmosphere(altKm, &sigma, &delta, &theta);
      fprintf(Itxt, "%4i %9.3E %9.3E %6.4f ", 5*i,sigma,delta,theta);
      temp=(KELVIN2RANKINE*TZERO)*theta;
      pressure=(PZERO/PSF2NSM)*delta;
      density=(RHOZERO/SCF2KCM)*sigma;
      asound=(AZERO/FT2METERS)*sqrt(theta);
      fprintf(Itxt, "%5.1f %8.3E %9.3E %6.1f", temp,pressure,density,asound);
      viscosity=(1.0/PSF2NSM)*MetricViscosity(theta);
      kinematicViscosity=viscosity/density;
      fprintf(Itxt, "%6.3f %8.2E\n", 1.0E6*viscosity, kinematicViscosity);
    }
  fclose(Itxt);
}   /* ------------------------------------- End of Procedure LongUSTable. */

/* ======================================================================= */
void ShortUSTable(void)
{
  FILE *Itxt;
  float altKm;
  float sigma,delta,theta;
  float temp,pressure,density, asound;
  float viscosity, kinematicViscosity, vratio;
  int i;

  Itxt=fopen("us2c.prt", "w");
  fprintf(Itxt, " alt  sigma  delta  theta");
  fprintf(Itxt, "  temp  press    dens     a     visc   k.visc ratio\n" );
  fprintf(Itxt, " Kft                     ");
  fprintf(Itxt, "  degR   psf   s/cu.ft   fps s/ft-sec sq.ft/s   1/ft\n");

  for(i=-1;i<=65;i++)
    {
      altKm=i*FT2METERS;
      SimpleAtmosphere(altKm, &sigma, &delta, &theta);
      fprintf(Itxt, "%4i %6.4f %6.4f %6.4f", i, sigma, delta, theta);
      temp=KELVIN2RANKINE*TZERO*theta;
      pressure=PZERO*delta/47.88;
      density=RHOZERO*sigma/515.379;
      asound=(AZERO/FT2METERS)*sqrt(theta);
      fprintf(Itxt, "%6.1f %6.1f %9.7f %6.1f", temp,pressure,density,asound);
      viscosity=(1.0/PSF2NSM)*MetricViscosity(theta);
      kinematicViscosity=viscosity/density;
      vratio=asound/kinematicViscosity;
      fprintf(Itxt, "%6.3f %8.2E %4.2f\n",
	 1.0E6*viscosity, kinematicViscosity, 1.0E-6*vratio);

    }
  fclose(Itxt);

}   /* ------------------------------------ End of Procedure ShortUSTable. */

/* ======================================================================= */
void LongSITable()
{

  FILE *Itxt;
  int i;
  float altKm;
  float sigma,delta,theta;
  float temp,pressure,density, asound;
  float viscosity, kinematicViscosity;

  Itxt=fopen("si1c.prt", "w");
  fprintf(Itxt, " alt    sigma      delta    theta");
  fprintf(Itxt, "  temp   press      dens     a   visc   k.visc\n" );
  fprintf(Itxt, "  Km                             ");
  fprintf(Itxt, "    K   N/sq.m    kg/cu.m  m/sec kg/m-s sq.m/s\n");

  for(i=-1; i<=43; i++)
    {
      altKm=2*i;                            /* changes integer to floating */
      Atmosphere(altKm, &sigma, &delta, &theta);
      fprintf(Itxt, "%4i %8.4E %8.4E %6.4f", 2*i, sigma,delta,theta);
      temp=TZERO*theta;
      pressure=PZERO*delta;
      density=RHOZERO*sigma;
      asound=AZERO*sqrt(theta);
      fprintf(Itxt, "%6.1f %8.3E %8.3E %5.1f", temp,pressure,density,asound);
      viscosity=MetricViscosity(theta);
      kinematicViscosity=viscosity/density;
      fprintf(Itxt, "%6.2f %8.2E\n", 1.0E6*viscosity, kinematicViscosity);
    }
  fclose(Itxt);
}   /* ------------------------------------- End of Procedure LongSITable. */

/* ======================================================================= */
void ShortSITable()
{
  FILE *Itxt;
  int i;
  float altKm;
  float sigma,delta,theta;
  float temp,pressure,density, asound;
  float viscosity, kinematicViscosity, vratio;

  Itxt=fopen("si2c.prt", "w");
  fprintf(Itxt, " alt  sigma  delta  theta");
  fprintf(Itxt, "  temp  press  dens   a    visc   k.visc ratio\n" );
  fprintf(Itxt, "  Km                     ");
  fprintf(Itxt, "    K  N/sq.m  kcm  m/sec kg/m-s  sq.m/s  1/m\n");

  for(i=-1; i<=40; i++)
    {
      altKm=0.5*i;
      SimpleAtmosphere(altKm, &sigma, &delta, &theta);
      fprintf(Itxt, "%4.1f %6.4f %6.4f %6.4f", altKm, sigma, delta, theta);
      temp=TZERO*theta;
      pressure=PZERO*delta;
      density=RHOZERO*sigma;
      asound=AZERO*sqrt(theta);
      fprintf(Itxt, "%6.1f %6.0f %5.3f %5.1f", temp, pressure,density,asound);
      viscosity=MetricViscosity(theta);
      kinematicViscosity=viscosity/density;
      vratio=asound/kinematicViscosity;
      fprintf(Itxt, "%6.2f %8.2E %5.2f\n",
	 1.0E6*viscosity, kinematicViscosity, 1.0E-6*vratio);
    }
  fclose(Itxt);

}   /* ------------------------------------ End of Procedure ShortSITable. */

/* ======================================================================= */
void Atmosphere(const float  alt,              /* geometric altitude, km.  */
		float *sigma,        /* density/sea-level standard density */
		float *delta,      /* pressure/sea-level standard pressure */
		float *theta)    /* temperature/sea-level std. temperature */

/* Compute the temperature,density,and pressure in the standard atmosphere */
/* Correct to 86 km.  Only approximate thereafter.    */
{
  float REARTH=6369.0;                         /* radius of the Earth (km) */
  float GMR = 34.163195;                                   /* gas constant */

  int i,j,k;
  float h, tgrad, deltah, tbase, tlocal;

#define NTAB 8
  static float htab[NTAB] = {0.0, 11.0, 20.0, 32.0, 47.0,
			    51.0, 71.0, 84.852 };
  static float ttab[NTAB] = { 288.15, 216.65, 216.65, 228.65, 270.65,
			      270.65, 214.65, 186.946 };
  static float ptab[NTAB] = { 1.0, 2.2336110E-1, 5.4032950E-2, 8.5666784E-3,
		      1.0945601E-3, 6.6063531E-4, 3.9046834E-5, 3.68501E-6 };
  static float gtab[NTAB] = { -6.5, 0, 1.0, 2.8, 0, -2.8, -2.0, 0 };

  h=alt*REARTH/(alt+REARTH);        /* geometric to geopotential altitude  */

  i=0; j=NTAB-1;                      /* starting values for binary search */
  do
    {
      k=(i+j)/2;
      if (h < htab[k]) j=k; else i=k;
    }  while (j > i+1);

  tgrad=gtab[i];
  tbase=ttab[i];
  deltah=h-htab[i];
  tlocal=tbase+tgrad*deltah;
  *theta=tlocal/ttab[0];                            /*  temperature ratio  */

  if (tgrad == 0.0)                            /*  compute pressure ratio  */
    *delta=ptab[i]*exp(-GMR*deltah/tbase);
  else
    *delta=ptab[i]*pow(tbase/tlocal, GMR/tgrad);

  *sigma= *delta / *theta;                              /*  density ratio  */
}   /*  --------------------------------------- End of function Atmosphere */

/* ======================================================================= */
void SimpleAtmosphere(
	const float  alt,                      /* geometric altitude, km.  */
	float *sigma,                /* density/sea-level standard density */
	float *delta,              /* pressure/sea-level standard pressure */
	float *theta)        /* temperature/sea-level standard temperature */

/* Compute the temperature,density,and pressure in the standard atmosphere */
/* Correct to 20 km.  Only approximate thereafter.    */
{
  float REARTH = 6369.0;                       /* radius of the Earth (km) */
  float GMR    = 34.163195;                                /* gas constant */

  float h;

  h=alt*REARTH/(alt+REARTH);        /* geometric to geopotential altitude  */

  if (h<11.0)
    {                                                       /* Troposphere */
      *theta=(288.15-6.5*h)/288.15;
      *delta=pow(*theta, GMR/6.5);
    }
  else
    {                                                      /* Stratosphere */
      *theta=216.65/288.15;
      *delta=0.2233611*exp(-GMR*(h-11.0)/216.65);
    }

  *sigma= *delta / *theta;
}   /* ---------------------------------- End of function SimpleAtmosphere */

/* ======================================================================= */
float MetricViscosity(const float theta)
{
  float t=theta*TZERO;                           /* temperature in kelvins */
  return BETAVISC*sqrt(t*t*t)/(t+SUTH);         /* viscosity in kg/(m-sec) */
}   /* ----------------------------------- End of function MetricViscosity */
