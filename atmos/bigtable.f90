!+
!PROGRAM BigTable                                      ! \atmos\V16\bigtable.f90
! ------------------------------------------------------------------------------
! PURPOSE - Make tables of atmospheric properties from 0 to 1000 km.
!
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software


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
!  1Jan01  1.8   RLC   Converted to BigTable, going to 1000 km.
!-------------------------------------------------------------------------------



!+
MODULE Atmosphere1976
! ------------------------------------------------------------------------------
! PURPOSE - Compute properties of the U.S. Standard Atmosphere 1976
! AUTHORS - Steven S. Pietrobon, Small World
!           Ralph L. Carmichael, Public Domain Aeronautical Software
!
!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 28Feb95  0.1   RLC   Assembled several old codes
!  1Aug00  0.2   RLC   Copied from old Tables76
! 23Aug00  0.3   RLC   Added NitrogenNumber using QUANC8 (later removed)
! 24Aug00  0.4   RLC   Added KineticTemperatureGradient
! 30Aug00  0.5   RLC   Corrected numerous errors
! 30Dec00  0.6   RLC   Adapted UpperAtmosphere from Pietrobon's Unofficial
!                        Australian Standard Atmosphere
! 29Dec10  0.7   RLC   Restructured code; made tables non-dimensional
!-------------------------------------------------------------------------------
IMPLICIT NONE
  CHARACTER(LEN=*),PUBLIC,PARAMETER:: ATM76_VERSION = "0.75 (6 January 2011)"
  REAL,PRIVATE,PARAMETER:: PI = 3.14159265
  REAL,PARAMETER,PRIVATE:: POLAR_RADIUS = 6356.7523
  REAL,PARAMETER,PRIVATE:: EQUATORIAL_RADIUS = 6378.1370
  REAL,PARAMETER,PRIVATE:: RADIUS_AT_45DEG_LAT = 6367.4895
  REAL,PARAMETER,PRIVATE:: AUTHALIC_RADIUS = 6371.0012
  REAL,PARAMETER,PRIVATE:: VOLUMETRIC_RADIUS = 6371.0008
  REAL,PARAMETER:: REARTH = POLAR_RADIUS              ! radius of the Earth (km)
  REAL,PARAMETER:: GZERO = 9.80665                 !  accel. of gravity, m/sec^2

  REAL,PARAMETER:: TZERO = 288.15                ! temperature at sealevel, K
  REAL,PARAMETER:: PZERO = 101325.0            ! pressure at sealevel, N/sq.m
  REAL,PARAMETER:: RHOZERO = 1.2250            ! density at sealevel, kg/cu.m
!  REAL,PARAMETER:: RSTAR = 8314.32       ! perfect gas constant, N-m/(kmol-K)
  REAL,PARAMETER:: ASOUNDZERO = 340.294   ! speed of sound at sealevel, m/sec



!  REAL,PARAMETER:: MZERO      = 28.9644 ! molecular weight of air at sealevel
!  REAL,PARAMETER:: AVOGADRO =  6.022169E26        ! 1/kmol, Avogadro constant
!  REAL,PARAMETER:: BOLTZMANN = 1.380622E-23        ! Nm/K, Boltzmann constant


  PUBLIC:: Atmosphere
  PRIVATE:: EvaluateCubic
  PRIVATE:: KineticTemperature
  PRIVATE:: LowerAtmosphere
  PRIVATE:: UpperAtmosphere
  PUBLIC:: ThermalConductivity
  PUBLIC:: Viscosity

CONTAINS

!+
PURE FUNCTION EvaluateCubic(a,fa,fpa, b,fb,fpb, u) RESULT(fu)
! ---------------------------------------------------------------------------
! PURPOSE - Evaluate a cubic polynomial defined by the function and the
!   1st derivative at two points
  REAL,INTENT(IN):: u   ! point where function is to be evaluated
  REAL,INTENT(IN):: a,fa,fpa   ! a, f(a), f'(a)  at first point
  REAL,INTENT(IN):: b,fb,fpb   ! b, f(b), f'(b)  at second point
  REAL:: fu                    ! computed value of f(u)

  REAL:: d,t,p
!----------------------------------------------------------------------------
  d=(fb-fa)/(b-a)
  t=(u-a)/(b-a)
  p=1.0-t

  fu = p*fa + t*fb - p*t*(b-a)*(p*(d-fpa)-t*(d-fpb))
  RETURN
END Function EvaluateCubic   ! ----------------------------------------------

!+
FUNCTION KineticTemperature(z) RESULT(t)
! ------------------------------------------------------------------------------
! PURPOSE - Compute kinetic temperature above 86 km.

  REAL,INTENT(IN)::  z     ! geometric altitude, km.                        
  REAL:: t     ! kinetic temperature, K

! PARAMETERS FOR THE DEFINITION OF KINETIC TEMPERATURE FROM 86 km to 1000 km
  REAL,PARAMETER:: Z7 =  86.0,  T7=186.8673
  REAL,PARAMETER:: z8 =  91.0,  T8=T7
  REAL,PARAMETER:: Z9 = 110.0,  T9=240.0
  REAL,PARAMETER:: Z10= 120.0, T10=360.0
!  REAL,PARAMETER:: Z11= 500.0, T11=999.2356   ! not used
  REAL,PARAMETER:: Z12=1000.0, T12=1000.0

  REAL,PARAMETER:: C1 = -76.3232  ! uppercase A in document
  REAL,PARAMETER:: C2 = 19.9429   ! lowercase a in document
  REAL,PARAMETER:: C3 = 12.0
  REAL,PARAMETER:: C4 = 0.01875   ! lambda in document
  REAL,PARAMETER:: TC = 263.1905

  REAL:: xx,yy
!-------------------------------------------------------------------------------
  IF (z <= Z8) THEN
    t=T7
  ELSE IF (z < Z9) THEN  
    xx=(z-Z8)/C2                        ! from Appendix B, p.223
    yy=SQRT(1.0-xx*xx)
    t=TC+C1*yy
  ELSE IF (z <= Z10) THEN
    t=T9+C3*(z-Z9)
  ELSE
    xx=(REARTH+Z10)/(REARTH+z)
    yy=(T12-T10)*EXP(-C4*(z-Z10)*xx)
    t=T12-yy
  END IF

  RETURN
END Function KineticTemperature   ! -----------------------------------------


!+
SUBROUTINE UpperAtmosphere(alt, sigma, delta, theta)
!   -------------------------------------------------------------------------
! PURPOSE - Compute the properties of the 1976 standard atmosphere from
!   86 km. to 1000 km.
! NOTE - The results for pressure and density are approximate and are
!   computed from a curve fit of selected points from the more precise
!   program ussa1976.

IMPLICIT NONE
  REAL,INTENT(IN)::  alt    ! geometric altitude, km.                        
  REAL,INTENT(OUT):: sigma  ! density/sea-level standard density              
  REAL,INTENT(OUT):: delta  ! pressure/sea-level standard pressure           
  REAL,INTENT(OUT):: theta  ! temperature/sea-level standard temperature

  INTEGER,PARAMETER:: NTABLE = 25

! altitude table (km)
  REAL,PARAMETER,DIMENSION(NTABLE):: Z = (/      &
     86.0,  93.0, 100.0, 107.0, 114.0, &
    121.0, 128.0, 135.0, 142.0, 150.0, &
    160.0, 170.0, 180.0, 190.0, 200.0, &
    220.0, 260.0, 300.0, 400.0, 500.0, &
    600.0, 700.0, 800.0, 900.0,1000.0 /)

! pressure table  (DELTA = P/P0)
  REAL,PARAMETER,DIMENSION(SIZE(Z)):: DELTA_TABLE = (/                  &
    3.6850E-6, 1.0660E-6, 3.1593E-7, 1.0611E-7, 4.3892E-8,  &
    2.3095E-8, 1.3997E-8, 9.2345E-9, 6.4440E-9, 4.4828E-9,  &
    2.9997E-9, 2.0933E-9, 1.5072E-9, 1.1118E-9, 8.3628E-10, &
    4.9494E-10, 1.9634E-10, 8.6557E-11, 1.4328E-11, 2.9840E-12, &
    8.1056E-13, 3.1491E-13, 1.6813E-13, 1.0731E-13, 7.4155E-14 /)
    

! density table  (SIGMA = RHO/RHO0)
  REAL,PARAMETER,DIMENSION(SIZE(Z)):: SIGMA_TABLE = (/                              &
    5.680E-6, 1.632E-6, 4.575E-7, 1.341E-7, 4.061E-8, &
    1.614E-8, 7.932E-9, 4.461E-9, 2.741E-9, 1.694E-9, &
    1.007E-9, 6.380E-10, 4.240E-10, 2.923E-10, 2.074E-10, &
    1.116E-10, 3.871E-11, 1.564E-11, 2.288E-12, 4.257E-13, &
    9.279E-14, 2.506E-14, 9.272E-15, 4.701E-15, 2.907E-15 /)
    
  REAL,PARAMETER,DIMENSION(SIZE(Z)):: LOGDELTA = LOG(DELTA_TABLE)
  REAL,PARAMETER,DIMENSION(SIZE(Z)):: LOGSIGMA = LOG(SIGMA_TABLE)

  REAL,PARAMETER,DIMENSION(SIZE(Z)):: DLOGDELTA = (/   & 
    -0.174061, -0.177924, -0.167029, -0.142755, -0.107859,   &
    -0.079322, -0.064664, -0.054879, -0.048260, -0.042767,   &
    -0.037854, -0.034270, -0.031543, -0.029384, -0.027632,   &
    -0.024980, -0.021559, -0.019557, -0.016735, -0.014530,   &
    -0.011314, -0.007677, -0.005169, -0.003944, -0.003612 /)
  REAL,PARAMETER,DIMENSION(SIZE(Z)):: DLOGSIGMA = (/   & 
    -0.172421, -0.182258, -0.178090, -0.176372, -0.154322,   &
    -0.113750, -0.090582, -0.075033, -0.064679, -0.056067,   &
    -0.048461, -0.043042, -0.038869, -0.035648, -0.033063,   &
    -0.029164, -0.024220, -0.021336, -0.017686, -0.016035,   &
    -0.014327, -0.011631, -0.008248, -0.005580, -0.004227 /)



  INTEGER:: i,j,k                                                  ! counters

  REAL:: p,rho
!----------------------------------------------------------------------------

  IF (alt >= Z(SIZE(Z))) THEN          ! trap altitudes greater than 1000 km.
    delta=7.42E-14                     ! ~value at 1000 km
    sigma=2.907E-15                    ! ~value at 1000 km
    theta=1000.0/TZERO                 ! value at 1000 km
    RETURN
  END IF

  i=1 
  j=SIZE(Z)                                    ! setting up for binary search
  DO
    k=(i+j)/2                                              ! integer division
    IF (alt < Z(k)) THEN
      j=k
    ELSE
      i=k
    END IF   
    IF (j <= i+1) EXIT
  END DO

  delta=EXP(EvaluateCubic(Z(i), LOGDELTA(i),   DLOGDELTA(i),          &
                        Z(i+1), LOGDELTA(i+1), DLOGDELTA(i+1), alt))

  sigma=EXP(EvaluateCubic(Z(i), LOGSIGMA(i),   DLOGSIGMA(i),          &
                        Z(i+1), LOGSIGMA(i+1), DLOGSIGMA(i+1), alt))

  theta=KineticTemperature(alt)/TZERO
  RETURN
END Subroutine UpperAtmosphere   ! ------------------------------------------

!+
SUBROUTINE LowerAtmosphere(alt, sigma, delta, theta)
!   -------------------------------------------------------------------------
! PURPOSE - Compute the properties of the 1976 standard atmosphere to 86 km.

  IMPLICIT NONE

  REAL,INTENT(IN)::  alt    ! geometric altitude, km.
  REAL,INTENT(OUT):: sigma  ! density/sea-level standard density              
  REAL,INTENT(OUT):: delta  ! pressure/sea-level standard pressure           
  REAL,INTENT(OUT):: theta  ! temperature/sea-level standard temperature

! use the value of REARTH from module specifications.
! If you extract this code for other use, uncomment next line
!  REAL,PARAMETER:: REARTH = 6371.0                ! radius of the Earth (km)
  REAL,PARAMETER:: GMR = 34.163195                     ! hydrostatic constant

  INTEGER:: i,j,k                                                  ! counters
  REAL:: h                                       ! geopotential altitude (km)
  REAL:: tgrad, tbase      ! temperature gradient and base temp of this layer
  REAL:: tlocal                                           ! local temperature
  REAL:: deltah                             ! height above base of this layer

  INTEGER,PARAMETER:: NTAB=8       ! number of entries in the defining tables
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
END Subroutine LowerAtmosphere   ! ------------------------------------------

!+  
PURE FUNCTION ThermalConductivity(t) RESULT(k)   
!   -------------------------------------------------------------------------
! PURPOSE - Compute the coefficient of thermal conductivity at a
!    given temperature

  REAL,INTENT(IN) :: t                ! temperature, K

  REAL:: k  ! coefficient of thermal conductivity, watts per meter per kelvin
  REAL,PARAMETER:: C1=2.64638E-3, C2=245.4
! NOTE - This empirical equation is given in the reference document.
!----------------------------------------------------------------------------
  k=C1*SQRT(t*t*t)/(t+C2*10.0**(-12.0/t))
  RETURN
END Function ThermalConductivity   ! ----------------------------------------



!+
PURE FUNCTION Viscosity(theta) RESULT(visc)
!   ----------------------------------------------------------------------------
! PURPOSE - Compute viscosity using Sutherland's formula.
!        Returns viscosity in kg/(meter-sec)
  REAL,INTENT(IN) :: theta                   ! temperature/sea-level temperature  
  REAL:: visc
  REAL:: temp                                              ! temperature, kelvin
  REAL,PARAMETER:: BETAVISC = 1.458E-6       ! viscosity term, N s/(sq.m sqrt(K)
  REAL,PARAMETER:: SUTH = 110.4                       ! Sutherland's constant, K

!-------------------------------------------------------------------------------
  temp=TZERO*theta
  visc=BETAVISC*Sqrt(temp*temp*temp)/(temp+SUTH)
  RETURN
END Function Viscosity   ! -----------------------------------------------------

!+
SUBROUTINE Atmosphere(alt,sigma,delta,theta)
!   ----------------------------------------------------------------------------
! PURPOSE - Compute the characteristics of the U.S. Standard Atmosphere 1976

  IMPLICIT NONE
  REAL,INTENT(IN)::  alt    ! geometric altitude, km.
  REAL,INTENT(OUT):: sigma  ! density/sea-level standard density             
  REAL,INTENT(OUT):: delta  ! pressure/sea-level standard pressure            
  REAL,INTENT(OUT):: theta  ! temperature/sea-level standard temperature   
!------------------------------------------------------------------------------
  IF (alt > 86.0) THEN
    CALL UpperAtmosphere(alt,sigma,delta,theta)
  ELSE
    CALL LowerAtmosphere(alt,sigma,delta,theta)
  END IF
  RETURN
END Subroutine Atmosphere   ! --------------------------------------------------

END Module Atmosphere1976   ! ==================================================



!+
PROGRAM BigTable                                        ! \atmos\bigtable.f90
! ---------------------------------------------------------------------------
! PURPOSE - Make tables of atmospheric properties from 0 to 1000 km.
!
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
USE Atmosphere1976
IMPLICIT NONE

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
!  1Jan01  1.8   RLC   Converted to BigTable, going to 1000 km.
! 31Dec10  1.9   RLC   Big revision of module atmosphere
! 24Apr13  1.96  RLC   Aligned headers a little better
!
! NOTES-

  CHARACTER(LEN=*),PARAMETER::GREETING =  &
        'bigTable - A Fortran 90 program to make atmosphere tables'
  CHARACTER(LEN=*),PARAMETER::AUTHOR=  &
    ' Ralph L. Carmichael, Public Domain Aeronautical Software'
  CHARACTER(LEN=*),PARAMETER::MODIFIER=' '   ! put your name here if you mod it
  CHARACTER(LEN=*),PARAMETER::VERSION=' 1.96 (24 April 2013)'
  CHARACTER(LEN=*),PARAMETER::FAREWELL= &
    'File bigtable.out has been added to your directory.'
    
!  REAL,PARAMETER:: FT2METERS = 0.3048       ! mult. ft. to get meters (exact)
!  REAL,PARAMETER:: KELVIN2RANKINE = 1.8             ! mult K to get deg R
!  REAL,PARAMETER:: PSF2NSM = 47.880258          ! mult lb/sq.ft to get N/sq.m
!  REAL,PARAMETER:: SCF2KCM = 515.379         ! mult slug/cu.ft to get kg/cu.m
    
!-------------------------------------------------------------------------------
  WRITE(*,*) GREETING
  WRITE(*,*) AUTHOR
  IF (MODIFIER /= ' ') WRITE(*,*) 'Modified by '//MODIFIER
  WRITE(*,*) VERSION

  CALL LongTable()

  WRITE(*,*) FAREWELL
  WRITE(*,*) 'bigTable has terminated successfully.'
  STOP

CONTAINS


!+
SUBROUTINE LongTable()
! ------------------------------------------------------------------------------

  CHARACTER(LEN=*),PARAMETER::HEAD1=' alt    sigma     delta    theta   temp'// &
    '    press      dens    a    visc  k.visc'
  CHARACTER(LEN=*),PARAMETER::HEAD2=' Km                                  K '// &
    '   N/sq.m   kg/cu.m   m/s  kg/m-s sq.m/s'
  CHARACTER(LEN=*),PARAMETER:: FMT = &
    '(I4,2ES11.4,F7.4,F7.1,2ES10.3,F6.1,F6.2,2ES9.2)'
  INTEGER,PARAMETER::ITXT=1
  INTEGER:: errCode
  INTEGER:: i
  REAL:: altKm                                ! altitude in kilometers
  REAL:: sigma,delta,theta
  REAL:: temp,pressure,density, asound
  REAL:: visc,kinematicViscosity
!-------------------------------------------------------------------------------
  OPEN(UNIT=ITXT, FILE='bigtable.out', STATUS='REPLACE', &
    IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')
  IF (errCode /= 0) THEN
    WRITE(*,*) "Unable to open output file bigtable.out"
    STOP
  END IF
  WRITE(ITXT,*) HEAD1
  WRITE(ITXT,*) HEAD2

  DO i=0,1000,5
    altKm=REAL(i)
    CALL Atmosphere(altKm, sigma,delta,theta)
    temp=TZERO*theta
    pressure=PZERO*delta
    density=RHOZERO*sigma
    asound=ASOUNDZERO*Sqrt(theta)
    visc=Viscosity(theta)
    kinematicViscosity=visc/density
    WRITE(ITXT,FMT) i, sigma, delta, theta, temp, pressure,  &
      density, asound, 1.0E6*visc, kinematicViscosity
  END DO
  Close(Itxt)
  RETURN
END Subroutine LongTable   ! ---------------------------------------------------

END Program BigTable   ! =======================================================
