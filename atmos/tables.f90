MODULE PhysicalConstants
  IMPLICIT NONE

  REAL,PARAMETER:: REARTH = 6369.0                 ! radius of the Earth (km)
  REAL,PARAMETER:: GMR = 34.163195                             ! gas constant

  REAL,PARAMETER:: FT2METERS = 0.3048               ! mult. ft. to get meters (exact)
  REAL,PARAMETER:: KELVIN2RANKINE = 1.8             ! mult deg K to get deg R
  REAL,PARAMETER:: PSF2NSM = 47.880258          ! mult lb/sq.ft to get N/sq.m
  REAL,PARAMETER:: SCF2KCM = 515.379         ! mult slug/cu.ft to get kg/cu.m
  REAL,PARAMETER:: TZERO = 288.15            ! temperature at sealevel, deg K
  REAL,PARAMETER:: PZERO = 101325.0            ! pressure at sealevel, N/sq.m
  REAL,PARAMETER:: RHOZERO = 1.2250            ! density at sealevel, kg/cu.m
  REAL,PARAMETER:: ASOUNDZERO = 340.294   ! speed of sound at sealevel, m/sec

  REAL,PARAMETER:: BETAVISC = 1.458E-6 ! viscosity term, N sec/(sq.m sqrt(deg K)
  REAL,PARAMETER:: SUTH = 110.4              ! Sutherland's constant, deg K

END MODULE PhysicalConstants   ! --------------------------------------------

!+
PROGRAM Tables                                         ! \atmtable\tables.f90
! ---------------------------------------------------------------------------
! PURPOSE - Make tables of atmospheric properties
!
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
!
!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 28Feb95  1.0   RLC   Assembled several old codes
!  7Jul95  1.1   RLC   Renamed routines, coded viscosity in US units
! 26Jul95  1.2   RLC   Converted f77 program to f90 style.
!                      Made Module PhysicalConstants
!  6Aug95  1.3   RLC   Replaced 1962 tables with 1976 tables
! 25Aug95  1.4   RLC   Added MODIFIER
! 19Sep95  1.5   RLC   Added a little precision to pressure tables
!  3Oct95  1.6   RLC   Numerous style changes to run under Elf90
! 29Nov96  1.7   RLC   Error message if unable to open output file
! 29Oct01  1.8   RLC   Renamed output files
!
! NOTES-

  IMPLICIT NONE
!============================================================================
!     C O N S T A N T S                                                     |
!============================================================================
  CHARACTER(LEN=*),PARAMETER::GREETING =  &
        ' tables - A Fortran 95 program to make atmosphere tables'
  CHARACTER(LEN=*),PARAMETER::AUTHOR=  &
    ' Ralph L. Carmichael, Public Domain Aeronautical Software'
  CHARACTER(LEN=*),PARAMETER::MODIFIER=' '   ! put your name here if you mod it
  CHARACTER(LEN=*),PARAMETER::VERSION=' 1.8 (29Oct01)'
  CHARACTER(LEN=*),PARAMETER::FAREWELL= &
    'Four files have been added to your directory.'
!-----------------------------------------------------------------------
  WRITE(*,*) GREETING
  WRITE(*,*) AUTHOR
  IF (MODIFIER /= ' ') WRITE(*,*) 'Modified by '//MODIFIER
  WRITE(*,*) VERSION

  CALL LongUSTable()
  CALL ShortUSTable()
  CALL LongSItable()
  CALL ShortSItable()

  WRITE(*,*) FAREWELL
  WRITE(*,*) 'tables has terminated successfully.'
  STOP

CONTAINS

!+
SUBROUTINE LongUSTable()
! ---------------------------------------------------------------------------
  USE PhysicalConstants

  CHARACTER(LEN=*),PARAMETER::HEAD1=' alt   sigma    delta   theta  temp'// &
      '   press    dens      a    visc  k.visc'
  CHARACTER(LEN=*),PARAMETER::HEAD2=' Kft                           degR'// &
      ' lb/sq.ft s/cu.ft    fps s/ft-s sq.ft/s'

  INTEGER,PARAMETER::ITXT=1
  CHARACTER(LEN=*),PARAMETER:: FMT = &
    '(I4,2ES9.3E1,F7.4,F6.1,2ES9.3E1,F7.1,F6.3,ES8.2E1)'

  INTEGER:: errCode
  INTEGER:: i                              ! altitude in thousand feet
  REAL:: altKm                                       !  altitude in Km
  REAL:: sigma,delta,theta       ! density,pressure,temperature ratios
  REAL:: temp,pressure,density
  REAL:: asound                             ! speed of sound in ft/sec
  REAL:: viscosity                       ! viscosity in slugs/(ft*sec)
  REAL:: kinematicViscosity                                ! sq.ft/sec
!----------------------------------------------------------------------------
  OPEN(UNIT=ITXT, FILE='us1f.prt', &
     IOSTAT=errCode, STATUS='REPLACE', ACTION='WRITE')
  IF (errCode /= 0) THEN
    WRITE(*,*) "Unable to open output file us1f.prt"
    STOP
  END IF
  WRITE(ITXT,*) HEAD1
  WRITE(ITXT,*) HEAD2

  DO i=-1,56   ! vary alt (ft) from -5000 to 280,000 by 5000
    altKm=REAL(5*i)*FT2METERS
    CALL Atmosphere(altKm, sigma,delta,theta)
    temp=KELVIN2RANKINE*TZERO*theta
    pressure=(PZERO/PSF2NSM)*delta
    density=(RHOZERO/SCF2KCM)*sigma
    asound=(ASOUNDZERO/FT2METERS)*Sqrt(theta)
    viscosity=(1.0/PSF2NSM)*MetricViscosity(theta)
    kinematicViscosity=viscosity/density
    WRITE(ITXT, FMT) 5*i, sigma, delta, theta, temp, pressure, &
      density, asound, 1E6*viscosity, kinematicViscosity
  END DO
  Close(UNIT=ITXT)
  RETURN
 
END Subroutine LongUStable   ! ----------------------------------------------

!+
SUBROUTINE ShortUStable()
!   --------------------------------------------------------------------
  USE PhysicalConstants

! L O C A L   C O N S T A N T S
  CHARACTER(LEN=*),PARAMETER::HEAD1=' alt  sigma  delta  theta  temp'//  &
      '  press    dens     a    visc  k.visc  ratio'
  CHARACTER(LEN=*),PARAMETER::HEAD2=' Kft                       degR   psf  '//  &
      ' s/cu.ft   fps s/ft-sec sq.ft/s 1/ft'
  CHARACTER(LEN=*),PARAMETER:: FMT = &
      '(I4,3F7.4,F6.1,F7.1,F10.7,F7.1,F6.3,ES8.2E1,F6.2)'
  INTEGER,PARAMETER::ITXT=1

! L O C A L   V A R I A B L E S
  INTEGER:: errCode
  INTEGER:: i !     altitude in thousand feet
  REAL:: altKm  !
  REAL:: sigma,delta,theta
  REAL:: temp,pressure,density,asound
  REAL:: viscosity,kinematicViscosity,vratio
!----------------------------------------------------------------------------
  OPEN(UNIT=ITXT, FILE='us2f.prt', &
    IOSTAT=errCode, STATUS='REPLACE', ACTION='WRITE')
  IF (errCode /= 0) THEN
    WRITE(*,*) "Unable to open output file us2f.prt"
    STOP
  END IF
  WRITE(ITXT,*) HEAD1
  WRITE(ITXT,*) HEAD2
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
    WRITE(ITXT,FMT) i, sigma, delta, theta, temp, pressure,    &
      density, asound, 1E6*viscosity, kinematicViscosity, vratio
  END DO
  Close(Itxt)
  RETURN
END Subroutine ShortUSTable   ! ---------------------------------------------

!+
SUBROUTINE LongSItable()
! ---------------------------------------------------------------------------
  USE PhysicalConstants

! L O C A L   C O N S T A N T S
  CHARACTER(LEN=*),PARAMETER::HEAD1=' alt    sigma     delta   theta  temp'// &
    '   press    dens     a    visc  k.visc'
  CHARACTER(LEN=*),PARAMETER::HEAD2='  Km                             degK'// &
    '  N/sq.m   kg/cu.m m/sec kg/m-s sq.m/s'
  CHARACTER(LEN=*),PARAMETER:: FMT = &
    '(I4,2ES10.4E1,F7.4,F6.1,2ES9.3E1,F6.1,F6.2,2ES8.2E1)'
  INTEGER,PARAMETER::ITXT=1

! L O C A L   V A R I A B L E S
  INTEGER:: errCode
  INTEGER:: i
  REAL:: altKm                                ! altitude in kilometers
  REAL:: sigma,delta,theta
  REAL:: temp,pressure,density, asound
  REAL:: viscosity,kinematicViscosity
!----------------------------------------------------------------------------
  OPEN(UNIT=ITXT, FILE='si1f.prt', &
    IOSTAT=errCode, STATUS='REPLACE', ACTION='WRITE')
  IF (errCode /= 0) THEN
    WRITE(*,*) "Unable to open output file si1f.prt"
    STOP
  END IF
  WRITE(ITXT,*) HEAD1
  WRITE(ITXT,*) HEAD2

  DO i=-2,86,2
    altKm=REAL(i)
    CALL Atmosphere(altKm, sigma,delta,theta)
    temp=TZERO*theta
    pressure=PZERO*delta
    density=RHOZERO*sigma
    asound=ASOUNDZERO*Sqrt(theta)
    viscosity=MetricViscosity(theta)
    kinematicViscosity=viscosity/density
    WRITE(ITXT,FMT) i, sigma, delta, theta, temp, pressure,  &
      density, asound, 1.0E6*viscosity, kinematicViscosity
  END DO
  Close(Itxt)
  RETURN
END Subroutine LongSItable   ! ----------------------------------------------

!+
SUBROUTINE ShortSItable()
! ---------------------------------------------------------------------------
  USE PhysicalConstants

! L O C A L   C O N S T A N T S
  CHARACTER(LEN=*),PARAMETER::HEAD1=' alt  sigma  delta  theta  temp'//   &
    '  press  dens   a    visc  k.visc ratio'
  CHARACTER(LEN=*),PARAMETER::HEAD2='  Km                       degK N/sq.m'// &
    '  kcm  m/sec kg/m-s sq.m/s  1/m'
  CHARACTER(LEN=*),PARAMETER:: FMT = &
    '(F4.1,3F7.4,F6.1,I7,F6.3,F6.1,F6.2,ES8.2E1,F6.2)'
  INTEGER,PARAMETER::ITXT=1                   ! unit number of output file

! L O C A L   V A R I A B L E S
  INTEGER:: errCode
  INTEGER:: i                                   !  steps thru altitude
  REAL:: altKm
  REAL:: sigma,delta,theta
  REAL:: temp,pressure,density, asound
  REAL:: viscosity,kinematicViscosity,vratio
!----------------------------------------------------------------------------
  OPEN(UNIT=ITXT, FILE='si2f.prt', &
    IOSTAT=errCode, STATUS='REPLACE', ACTION='WRITE')
  IF (errCode /= 0) THEN
    WRITE(*,*) "Unable to open output file si2f.prt"
    STOP
  END IF
  WRITE(ITXT,*) HEAD1
  WRITE(ITXT,*) HEAD2

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
    WRITE(ITXT,FMT) altKm, sigma, delta, theta, temp, INT(pressure),  &
      density, asound, 1.0E6*viscosity, kinematicViscosity, 1.0E-6*vratio
  END DO
  Close(Itxt)
  RETURN
 
END Subroutine ShortSITable   ! ---------------------------------------------

!+
SUBROUTINE Atmosphere(alt, sigma, delta, theta)
!   -------------------------------------------------------------------------
! PURPOSE - Compute the properties of the 1976 standard atmosphere to 86 km.

  IMPLICIT NONE
!============================================================================
!     A R G U M E N T S                                                     |
!============================================================================
  REAL,INTENT(IN)::  alt    ! geometric altitude, km.                        
  REAL,INTENT(OUT):: sigma  ! density/sea-level standard density              
  REAL,INTENT(OUT):: delta  ! pressure/sea-level standard pressure           
  REAL,INTENT(OUT):: theta  ! temperature/sea-level standard temperature
!============================================================================
!     L O C A L   C O N S T A N T S                                         |
!============================================================================
  REAL,PARAMETER:: REARTH = 6369.0                 ! radius of the Earth (km)
  REAL,PARAMETER:: GMR = 34.163195                             ! gas constant
  INTEGER,PARAMETER:: NTAB=8       ! number of entries in the defining tables
!============================================================================
!     L O C A L   V A R I A B L E S                                         |
!============================================================================
  INTEGER:: i,j,k                                                  ! counters
  REAL:: h                                       ! geopotential altitude (km)
  REAL:: tgrad, tbase      ! temperature gradient and base temp of this layer
  REAL:: tlocal                                           ! local temperature
  REAL:: deltah                             ! height above base of this layer
!============================================================================
!     L O C A L   A R R A Y S   ( 1 9 7 6   S T D.  A T M O S P H E R E )   |
!============================================================================
  REAL,DIMENSION(NTAB),PARAMETER:: htab= &
                          (/0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852/)
  REAL,DIMENSION(NTAB),PARAMETER:: ttab= &
          (/288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946/)
  REAL,DIMENSION(NTAB),PARAMETER:: ptab= &
               (/1.0, 2.233611E-1, 5.403295E-2, 8.5666784E-3, 1.0945601E-3, &
                                     6.6063531E-4, 3.9046834E-5, 3.68501E-6/)
  REAL,DIMENSION(NTAB),PARAMETER:: gtab= &
                                (/-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0/)
!----------------------------------------------------------------------------
  h=alt*REARTH/(alt+REARTH)      ! convert geometric to geopotential altitude

  i=1 
  j=NTAB                                  ! setting up for binary search
  DO
    k=(i+j)/2                                              ! integer division
    IF (h < htab(k)) THEN
      j=k
    ELSE
      i=k
    END IF   
    IF (j <= i+1) EXIT
  END DO

  tgrad=gtab(i)                                     ! i will be in 1...NTAB-1
  tbase=ttab(i)
  deltah=h-htab(i)
  tlocal=tbase+tgrad*deltah
  theta=tlocal/ttab(1)                                    ! temperature ratio

  IF (tgrad == 0.0) THEN                                     ! pressure ratio
    delta=ptab(i)*EXP(-GMR*deltah/tbase)
  ELSE
    delta=ptab(i)*(tbase/tlocal)**(GMR/tgrad)
  END IF

  sigma=delta/theta                                           ! density ratio
  RETURN
END Subroutine Atmosphere   ! -----------------------------------------------

!+
SUBROUTINE SimpleAtmosphere(alt,sigma,delta,theta)
!   -------------------------------------------------------------------------
! PURPOSE - Compute the characteristics of the lower atmosphere.

! NOTES-Correct to 20 km. Only approximate above there

  IMPLICIT NONE
!============================================================================
!     A R G U M E N T S                                                     |
!============================================================================
  REAL,INTENT(IN)::  alt    ! geometric altitude, km.
  REAL,INTENT(OUT):: sigma  ! density/sea-level standard density             
  REAL,INTENT(OUT):: delta  ! pressure/sea-level standard pressure            
  REAL,INTENT(OUT):: theta  ! temperature/sea-level standard temperature   
!============================================================================
!     L O C A L   C O N S T A N T S                                         |
!============================================================================
  REAL,PARAMETER:: REARTH = 6369.0                ! radius of the Earth (km)
  REAL,PARAMETER:: GMR = 34.163195                            ! gas constant
!============================================================================
!     L O C A L   V A R I A B L E S                                         |
!============================================================================
  REAL:: h   ! geopotential altitude
!----------------------------------------------------------------------------
  h=alt*REARTH/(alt+REARTH)      ! convert geometric to geopotential altitude

  IF (h < 11.0) THEN
    theta=1.0+(-6.5/288.15)*h                                   ! Troposphere
    delta=theta**(GMR/6.5)
  ELSE
    theta=216.65/288.15                                        ! Stratosphere
    delta=0.2233611*EXP(-GMR*(h-11.0)/216.65)
  END IF

  sigma=delta/theta
  RETURN
END Subroutine SimpleAtmosphere   ! -----------------------------------------

!+
FUNCTION MetricViscosity(theta) RESULT(visc)
!   -------------------------------------------------------------------------
! PURPOSE - Compute viscosity using Sutherland's formula.
!        Returns viscosity in kg/(meter-sec)
  USE PhysicalConstants

  IMPLICIT NONE
!============================================================================
!     A R G U M E N T S                                                     |
!============================================================================
  REAL,INTENT(IN) :: theta                ! temperature/sea-level temperature  
!============================================================================
!     L O C A L   V A R I A B L E S                                         |
!============================================================================
  REAL:: visc
  REAL:: temp                              ! temperature in deg Kelvin
!----------------------------------------------------------------------------
  temp=TZERO*theta
  visc=BETAVISC*Sqrt(temp*temp*temp)/(temp+SUTH)
  RETURN
END Function MetricViscosity   ! --------------------------------------------

END Program Tables   ! ------------------------------------------------------
