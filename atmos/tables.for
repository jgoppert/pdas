*+
      PROGRAM Tables                              ! \atmtable\tables.for
*   --------------------------------------------------------------------
*     PURPOSE - Make tables of atmospheric properties
*
*     AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
*
*     REVISION HISTORY
*   DATE  VERS PERSON  STATEMENT OF CHANGES
* 28Feb95  1.0   RLC   Assembled several old codes
*  7Jul95  1.1   RLC   Renamed routines, coded viscosity in US units
*  5Aug95  1.2   RLC   Replaced 1962 tables with 1976 tables
* 25Aug95  1.3   RLC   Added MODIFIER
* 19Sep95  1.4   RLC   Added a little accuracy to the pressure tables
* 29Nov96  1.5   RLC   Added message if unable to open output files
* 31Oct01  1.6   RLC   Renamed table names
*
*     NOTES-
*
      IMPLICIT NONE
*
      EXTERNAL  LongUSTable, ShortUSTable
      EXTERNAL  LongSITable, ShortSITable
************************************************************************
*     C O N S T A N T S                                                *
************************************************************************
      CHARACTER GREETING*60, AUTHOR*60, VERSION*30, FAREWELL*60
      CHARACTER MODIFIER*63
      PARAMETER (GREETING=
     &   ' tables - A Fortran program to make atmosphere tables')
      PARAMETER (AUTHOR=
     &   ' Ralph L. Carmichael, Public Domain Aeronautical Software')
      PARAMETER (MODIFIER=' ')
      PARAMETER (VERSION=' 1.6 (31Oct01)' )
      PARAMETER (FAREWELL=
     & 'Four files added to your directory.')
*-----------------------------------------------------------------------
      WRITE (*,*) GREETING, AUTHOR, VERSION
*
      CALL LongUSTable
      CALL ShortUSTable
      CALL LongSITable
      CALL ShortSITable
*
      WRITE (*,*) FAREWELL
      WRITE(*,*) 'tables has terminated successfully.'
      STOP
*
      END ! --------------------------------- End of main program tables
*+
      SUBROUTINE LongUSTable
* ----------------------------------------------------------------------
*     PURPOSE - Make an atmosphere table that varies altitude from
*       -5000 feet to 280,000 feet by 5000 foot intervals
      IMPLICIT NONE
*
      REAL MetricViscosity
      EXTERNAL Atmosphere, MetricViscosity
************************************************************************
*     L O C A L   C O N S T A N T S                                    *
************************************************************************
      CHARACTER HEAD8*78, HEAD9*78
      PARAMETER (HEAD8=' alt   sigma    delta   theta  temp'//
     & '   press    dens      a    visc  k.visc')
      PARAMETER (HEAD9=' Kft                           degR'//
     & ' lb/sq.ft s/cu.ft    fps s/ft-s sq.ft/s')

      REAL FT2METERS                   ! mult. ft. to get meters (exact)
      REAL KELVIN2RANKINE                    ! mult kelvins to get deg R
      REAL PSF2NSM                         ! mult lb/sq.ft to get N/sq.m
      REAL SCF2KCM                      ! mult slug/cu.ft to get kg/cu.m
      REAL TZERO                      ! temperature at sealevel, kelvins
      REAL PZERO                          ! pressure at sealevel, N/sq.m
      REAL RHOZERO                        ! density at sealevel, kg/cu.m
      REAL ASOUNDZERO                ! speed of sound at sealevel, m/sec
      PARAMETER(FT2METERS = 0.3048)
      PARAMETER(KELVIN2RANKINE = 1.8)
      PARAMETER(PSF2NSM=47.880258)
      PARAMETER(SCF2KCM=515.379)
      PARAMETER(TZERO=288.15)
      PARAMETER(PZERO=101325.0)
      PARAMETER(RHOZERO=1.2250)
      PARAMETER(ASOUNDZERO=340.294)
      INTEGER ITXT
      PARAMETER (ITXT=1)
************************************************************************
*     L O C A L   V A R I A B L E S                                    *
************************************************************************
      INTEGER i                              ! altitude in thousand feet
      REAL altKm                                       !  altitude in Km
      REAL sigma,delta,theta       ! density,pressure,temperature ratios
      REAL temp,pressure,density
      REAL asound                             ! speed of sound in ft/sec
      REAL viscosity                       ! viscosity in slugs/(ft*sec)
      REAL kinematicViscosity                                ! sq.ft/sec
!!!      REAL vratio        ! ratio of asound to kinematic viscosity (1/ft)
      INTEGER errCode
*-----------------------------------------------------------------------
      OPEN(UNIT=ITXT, FILE='us1f77.prt',                                &
     & STATUS='UNKNOWN', IOSTAT=errCode)
      IF (errCode .NE. 0) THEN
        WRITE(*,*) 'Unable to open output file us1f77.prt'
        STOP
      END IF
      WRITE(ITXT,*) HEAD8
      WRITE(ITXT,*) HEAD9

      DO i=-1,56   ! vary alt (ft) from -5000 to 280,000 by 5000
         altKm=REAL(5*i)*FT2METERS
         CALL Atmosphere(altKm, sigma,delta,theta)
         temp=KELVIN2RANKINE*TZERO*theta
         pressure=(PZERO/PSF2NSM)*delta
         density=(RHOZERO/SCF2KCM)*sigma
         asound=(ASOUNDZERO/FT2METERS)*Sqrt(theta)
         viscosity=(1.0/PSF2NSM)*MetricViscosity(theta)
         kinematicViscosity=viscosity/density
         WRITE(ITXT,5000) 5*i, sigma, delta, theta, temp, pressure,
     &      density, asound, 1E6*viscosity, kinematicViscosity
      END DO
      Close(1)
 5000 FORMAT(I4,1P2E9.3E1,0PF7.4,F6.1,1P2E9.3E1,0PF7.1,0PF6.3,1PE8.2E1)
      END   ! ----------------------- End of Subroutine LongEnglishTable
*+
      SUBROUTINE ShortUSTable
* ----------------------------------------------------------------------
*     PURPOSE - Make an atmosphere table that varies altitude from
*       -1000 feet to 65000 feet by 1000 foot intervals
      IMPLICIT NONE
*
      REAL MetricViscosity
      EXTERNAL SimpleAtmosphere, MetricViscosity
************************************************************************
*     L O C A L   C O N S T A N T S                                    *
************************************************************************
      CHARACTER HEAD8*78, HEAD9*78
      PARAMETER (HEAD8=' alt  sigma  delta  theta  temp  press'//
     & '    dens     a    visc  k.visc  ratio')
      PARAMETER (HEAD9=' Kft                       degR   psf  '//
     & ' s/cu.ft   fps s/ft-sec sq.ft/s 1/ft')
      REAL FT2METERS                    ! mult ft. to get meters (exact)
      REAL KELVIN2RANKINE                    ! mult kelvins to get deg R
      REAL PSF2NSM                         ! mult lb/sq.ft to get N/sq.m
      REAL SCF2KCM                      ! mult slug/cu.ft to get kg/cu.m
      REAL TZERO, PZERO, RHOZERO, ASOUNDZERO
      PARAMETER(FT2METERS = 0.3048)
      PARAMETER(KELVIN2RANKINE = 1.8)
      PARAMETER(PSF2NSM=47.880258)
      PARAMETER(SCF2KCM=515.379)
      PARAMETER(TZERO=288.15)      !    temperature at sealevel, kelvins
      PARAMETER(PZERO=101325.0)      !     pressure at sealevel, N/sq.m.
      PARAMETER(RHOZERO=1.2250)       !    density at sealevel, kg/cu.m.
      PARAMETER(ASOUNDZERO=340.294)  ! speed of sound at sealevel, m/sec
      INTEGER ITXT
      PARAMETER (ITXT=1)
************************************************************************
*     L O C A L   V A R I A B L E S                                    *
************************************************************************
      INTEGER i !     altitude in thousand feet
      REAL altKm  !
      REAL sigma,delta,theta
      REAL temp,pressure,density,asound
      REAL viscosity,kinematicViscosity,vratio
      INTEGER errCode
*-----------------------------------------------------------------------
      OPEN(UNIT=ITXT, FILE='us2f77.prt',                                &
     &  STATUS='UNKNOWN', IOSTAT=errCode)
      IF (errCode .NE. 0) THEN
        WRITE(*,*) 'Unable to open output file us2f77.prt'
        STOP
      END IF
      WRITE(ITXT,*) HEAD8
      WRITE(ITXT,*) HEAD9
      DO i=-1,65   ! vary alt (ft) from -1000 to 65000 by 1000
         altKm=REAL(i)*FT2METERS
         CALL SimpleAtmosphere(altKm, sigma,delta,theta)
         temp=(KELVIN2RANKINE*TZERO)*theta
         pressure=(PZERO/PSF2NSM)*delta
         density=(RHOZERO/SCF2KCM)*sigma
         asound=(ASOUNDZERO/FT2METERS)*Sqrt(theta)
         viscosity=(1.0/PSF2NSM)*MetricViscosity(theta)
         kinematicViscosity=viscosity/density
         vratio=1.0E-6*asound/kinematicViscosity
         WRITE(ITXT,5000) i, sigma, delta, theta, temp, pressure,
     &      density, asound, 1E6*viscosity, kinematicViscosity, vratio
      END DO
      Close(Itxt)
 5000 FORMAT(I4,3F7.4,F6.1,F7.1,F10.7,F7.1,F6.3,1PE8.2E1,0PF6.2)
      END ! ------------------------ End of subroutine ShortEnglishTable
*+
      SUBROUTINE LongSITable
*   --------------------------------------------------------------------
*     PURPOSE - Make an atmosphere table that varies altitude from
*       -2000 meters to 86000 meters by 2000 meter intervals
      IMPLICIT NONE
*
      REAL MetricViscosity
      EXTERNAL Atmosphere, MetricViscosity
************************************************************************
*     L O C A L   C O N S T A N T S                                    *
************************************************************************
      CHARACTER HEAD8*78, HEAD9*78
      PARAMETER (HEAD8=' alt    sigma     delta   theta  temp'//
     & '   press    dens     a    visc  k.visc')
      PARAMETER (HEAD9='  Km                               K '//
     & '  N/sq.m   kg/cu.m m/sec kg/m-s sq.m/s')
      REAL TZERO, PZERO, RHOZERO, ASOUNDZERO
      PARAMETER(TZERO=288.15)      !    temperature at sealevel, kelvins
      PARAMETER(PZERO=101325.0)      !     pressure at sealevel, N/sq.m.
      PARAMETER(RHOZERO=1.2250)     !      density at sealevel, kg/cu.m.
      PARAMETER(ASOUNDZERO=340.294)  ! speed of sound at sealevel, m/sec
      INTEGER ITXT
      PARAMETER (ITXT=1)
************************************************************************
*     L O C A L   V A R I A B L E S                                    *
************************************************************************
      INTEGER i
      REAL altKm                                ! altitude in kilometers
      REAL sigma,delta,theta
      REAL temp,pressure,density, asound
      REAL viscosity,kinematicViscosity
      INTEGER errCode
*-----------------------------------------------------------------------
      OPEN(UNIT=ITXT, FILE='si1f77.prt',                                &
     &  STATUS='UNKNOWN', IOSTAT=errCode)
      IF (errCode .NE. 0) THEN
        WRITE(*,*) 'Unable to open output file si1f77.prt'
        STOP
      END IF
      WRITE(ITXT,*) HEAD8
      WRITE(ITXT,*) HEAD9


      DO i=-2,86,2
         altKm=REAL(i)
         CALL Atmosphere(altKm, sigma,delta,theta)
         temp=TZERO*theta
         pressure=PZERO*delta
         density=RHOZERO*sigma
         asound=ASOUNDZERO*Sqrt(theta)
         viscosity=MetricViscosity(theta)
         kinematicViscosity=viscosity/density
         WRITE(ITXT,5000) i, sigma, delta, theta, temp, pressure,
     &      density, asound,
     &      1.0E6*viscosity, kinematicViscosity
      END DO
      Close(Itxt)
 5000 FORMAT(I4,1P2E10.4E1,0PF7.4,F6.1,1P2E9.3E1,0PF6.1,F6.2,1P2E8.2E1)
      END ! -------------------------- End of subroutine LongMetricTable
*+
      SUBROUTINE ShortSITable
*   --------------------------------------------------------------------
*     PURPOSE - Make an atmosphere table that varies altitude from
*       -1000 meters to 20000 meters by 500 meter intervals

      IMPLICIT NONE
*
      REAL MetricViscosity
      EXTERNAL SimpleAtmosphere, MetricViscosity
************************************************************************
*     L O C A L   C O N S T A N T S                                    *
************************************************************************
      CHARACTER HEAD8*78, HEAD9*78
      PARAMETER (HEAD8=' alt  sigma  delta  theta  temp  press'//
     & '  dens   a    visc  k.visc ratio')
      PARAMETER (HEAD9='  Km                         K  N/sq.m'//
     & '  kcm  m/sec kg/m-s sq.m/s  1/m')

      REAL TZERO, PZERO, RHOZERO, ASOUNDZERO
      PARAMETER(TZERO=288.15)         ! temperature at sealevel, kelvins
      PARAMETER(PZERO=101325.0)          ! pressure at sealevel, N/sq.m.
      PARAMETER(RHOZERO=1.2250)          ! density at sealevel, kg/cu.m.
      PARAMETER(ASOUNDZERO=340.294)  ! speed of sound at sealevel, m/sec
      INTEGER ITXT                          ! unit number of output file
      PARAMETER (ITXT=1)
************************************************************************
*     L O C A L   V A R I A B L E S                                    *
************************************************************************
      INTEGER i                                   !  steps thru altitude
      REAL altKm
      REAL sigma,delta,theta
      REAL temp,pressure,density, asound
      REAL viscosity,kinematicViscosity,vratio
      INTEGER errCode
*-----------------------------------------------------------------------
      OPEN(UNIT=ITXT, FILE='si2f77.prt',                                &
     &   STATUS='UNKNOWN', IOSTAT=errCode)
      IF (errCode .NE. 0) THEN
        WRITE(*,*) 'Unable to open output file si2f77.prt'
        STOP
      END IF
      WRITE(ITXT,*) HEAD8
      WRITE(ITXT,*) HEAD9

      DO i=-1,40
         altKm=0.5*REAL(i)
         CALL SimpleAtmosphere(altKm, sigma,delta,theta)
         temp=TZERO*theta
         pressure=PZERO*delta
         density=RHOZERO*sigma
         asound=ASOUNDZERO*Sqrt(theta)
         viscosity=MetricViscosity(theta)
         kinematicViscosity=viscosity/density
         vratio=asound/kinematicViscosity
         WRITE(ITXT,5000) altKm, sigma, delta, theta, temp,
     &      INT(pressure),  density, asound,
     &      1.0E6*viscosity, kinematicViscosity, 1.0E-6*vratio
      END DO
      Close(Itxt)
 5000 FORMAT(F4.1,3F7.4,F6.1,I7,F6.3,F6.1,F6.2,1PE8.2E1,0PF6.2)
      END ! ------------------------- End of subroutine ShortMetricTable
*+
      SUBROUTINE Atmosphere(alt, sigma, delta, theta)
*   --------------------------------------------------------------------
*     PURPOSE -
*
      IMPLICIT NONE
************************************************************************
*     A R G U M E N T S                                                *
************************************************************************
      REAL alt    ! geometric altitude, km.                           IN
      REAL sigma  ! density/sea-level standard density               OUT
      REAL delta  ! pressure/sea-level standard pressure             OUT
      REAL theta  ! temperature/sea-level standard temperature       OUT
************************************************************************
*     L O C A L   C O N S T A N T S                                    *
************************************************************************
      REAL REARTH, GMR
      PARAMETER(REARTH=6369.0)                ! radius of the Earth (km)
      PARAMETER(GMR = 34.163195)                          ! gas constant
      INTEGER NTAB            ! number of entries in the defining tables
      PARAMETER(NTAB=8)
************************************************************************
*     L O C A L   V A R I A B L E S                                    *
************************************************************************
      INTEGER i,j,k                                           ! counters
      REAL h                                ! geopotential altitude (km)
      REAL tgrad, deltah, tbase, tlocal
************************************************************************
*     L O C A L   A R R A Y S                                          *
************************************************************************
      REAL htab(NTAB),ttab(NTAB),ptab(NTAB),gtab(NTAB)
      DATA htab/0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852/
      DATA ttab/288.15, 216.65, 216.65, 228.65, 270.65,
     &                  270.65, 214.65, 186.946/
      DATA ptab/1.0, 2.233611E-1, 5.403295E-2, 8.5666784E-3, 
     & 1.0945601E-3, 6.6063531E-4, 3.9046834E-5, 3.68501E-6/
      DATA gtab/-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0/
*-----------------------------------------------------------------------
      h=alt*REARTH/(alt+REARTH)   !  geometric to geopotential altitude

      i=1
      j=NTAB
   10 k=(i+j)/2                                       ! integer division
      IF (h .LT. htab(k)) THEN
         j=k
      ELSE
         i=k
      ENDIF
      IF (j .GT. i+1) GOTO 10

      tgrad=gtab(i)                            ! i will be in 1...NTAB-1
      deltah=h-htab(i)
      tbase=ttab(i)

      tlocal=tbase+tgrad*deltah
      theta=tlocal/ttab(1)                           ! temperature ratio

      IF (tgrad .EQ. 0.0) THEN                          ! pressure ratio
         delta=ptab(i)*EXP(-GMR*deltah/tbase)
      ELSE
         delta=ptab(i)*(tbase/tlocal)**(GMR/tgrad)
      ENDIF

      sigma=delta/theta                                  ! density ratio
      END ! ------------------------------- End of subroutine Atmosphere
*+
      SUBROUTINE SimpleAtmosphere(alt,sigma,delta,theta)
*   --------------------------------------------------------------------
*     PURPOSE - Compute the characteristics of the lower atmosphere.
*
*     NOTES-Correct to 20 km. Only approximate above there
*
      IMPLICIT NONE
************************************************************************
*     A R G U M E N T S                                                *
************************************************************************
      REAL alt    ! geometric altitude, km.                           IN
      REAL sigma  ! density/sea-level standard density               OUT
      REAL delta  ! pressure/sea-level standard pressure             OUT
      REAL theta  ! temperature/sea-level standard temperature       OUT
************************************************************************
*     L O C A L   C O N S T A N T S                                    *
************************************************************************
      REAL REARTH                             ! radius of the Earth (km)
      REAL GMR                                            ! gas constant
      PARAMETER (REARTH=6369.0, GMR=34.163195)
************************************************************************
*     L O C A L   V A R I A B L E S                                    *
************************************************************************
      REAL h   ! geopotential altitude
*-----------------------------------------------------------------------
      h=alt*REARTH/(alt+REARTH)     ! geometric to geopotential altitude
*
      IF (h .LT. 11.0) THEN
         theta=1.0+(-6.5/288.15)*h                         ! Troposphere
         delta=theta**(GMR/6.5)
      ELSE
         theta=216.65/288.15                              ! Stratosphere
         delta=0.2233611*EXP(-GMR*(h-11.0)/216.65)
      ENDIF
*
      sigma=delta/theta
      END ! ------------------------- End of subroutine SimpleAtmosphere
*+
      REAL FUNCTION MetricViscosity(theta)
*   --------------------------------------------------------------------
*     PURPOSE - Compute viscosity using Sutherland's formula.
*        Returns viscosity in kg/(meter-sec)
*
      IMPLICIT NONE
************************************************************************
*     A R G U M E N T S                                                *
************************************************************************
      REAL theta             ! temperature/sea-level temperature      IN
************************************************************************
*     L O C A L   C O N S T A N T S                                    *
************************************************************************
      REAL BETAVISC         ! viscosity term,  N sec/(sq.m sqrt(kelvins)
      REAL SUTH                         ! Sutherland's constant, kelvins
      REAL TZERO                      ! temperature at sealevel, kelvins
      PARAMETER(BETAVISC= 1.458E-6, SUTH=110.4, TZERO= 288.15)
************************************************************************
*     L O C A L   V A R I A B L E S                                    *
************************************************************************
      REAL temp                                 ! temperature in kelvins
*-----------------------------------------------------------------------
      temp=TZERO*theta
      MetricViscosity=BETAVISC*Sqrt(temp*temp*temp)/(temp+SUTH)
      END ! ---------------------------- End of Function MetricViscosity
