!+
! PROGRAM TablesHTML
! ------------------------------------------------------------------------------
! PURPOSE - Make tables of atmospheric properties. Instead of making
!   uniformly spaced tables using monospaced font printing, this program
!   makes an HTML file and uses the <tables> construct to keep the
!   columns lined up properly. This program makes the same four tables as the
!   programs tables.f90, tables.cpp, etc. All four tables are in one HTML5 file.
!
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
!
!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 28Feb95  1.0   RLC   Assembled several old codes
!  7Jul95  1.1   RLC   Renamed routines, coded viscosity in US units
! 26Jul95  1.2   RLC   Converted f77 program to f90 style.
!  6Aug95  1.3   RLC   Replaced 1962 tables with 1976 tables
! 25Aug95  1.4   RLC   Added MODIFIER
! 19Sep95  1.5   RLC   Added a little precision to pressure tables
!  3Oct95  1.6   RLC   Numerous style changes to run under Elf90
! 29Nov96  1.7   RLC   Error message if unable to open output file
! 22Apr99  1.8   RLC   Made HTML output file
! 22Dec00  1.9   RLC   Corrected two small errors and removed 86+
! 29Oct01  2.0   RLC   Fixed it up for release 7 of PDAS.
! 09Dec10  2.1   RLC   Updated HTML to HTML 5
! 16Jul12  2.2   RLC   Restructured into modules
!

!+
MODULE AtmosphereProcedures
! ------------------------------------------------------------------------------
! PURPOSE - Collect the procedures for calculating pressure, density, 
! temperature, and viscosity in the standard atmosphere.

IMPLICIT NONE

! The standard assumes a spherical earth. From the actual dimensions of the
! ellipsoidal approximation of the earth, one may choose various values for the
! single value of REARTH. The value of POLAR_RADIUS is chosen because that is
! what was used for the tables in the big yellow book.
  REAL,PARAMETER,PRIVATE:: POLAR_RADIUS = 6356.766  ! also have seen 6356.7523
  REAL,PARAMETER,PRIVATE:: EQUATORIAL_RADIUS = 6378.1370
  REAL,PARAMETER,PRIVATE:: RADIUS_AT_45DEG_LAT = 6367.4895
  REAL,PARAMETER,PRIVATE:: AUTHALIC_RADIUS = 6371.0012    ! same surface area
  REAL,PARAMETER,PRIVATE:: VOLUMETRIC_RADIUS = 6371.0008  ! same volume
! the next line is where you decide what value to use for the radius of the earth
  REAL,PARAMETER:: REARTH = POLAR_RADIUS              ! radius of the Earth (km)

! The constant GMR (gravity times molecular weight divided by gas constant) has
!  units of kelvins per meter. Since I have coded the atmosphere routines to
!  use altitude in kilometers rather than meters, the factor of 1000 is needed
!  so that GMR will have units of kelvins per kilometer.

  REAL,PARAMETER:: GZERO = 9.80665 !  accel. of gravity, m/sec^2
  REAL,PARAMETER:: MZERO = 28.9644 !  molecular weight of air at sealevel
  REAL,PARAMETER:: RSTAR = 8314.32 !  perfect gas constant
  REAL,PARAMETER:: GMR = 1000*GZERO*MZERO/RSTAR
!  REAL(WP),PARAMETER:: GMR = 34.163195                   ! hydrostatic constant

CONTAINS


SUBROUTINE Atmosphere(alt, sigma, delta, theta)
!   -------------------------------------------------------------------------
! PURPOSE - Compute the properties of the 1976 standard atmosphere to 86 km.

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
!  REAL,PARAMETER:: GMR = 34.163195              ! hydrostatic constant, K/Km
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
                                              ! REARTH is global
  i=1 
  j=NTAB                                       ! setting up for binary search
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
    delta=ptab(i)*EXP(-GMR*deltah/tbase)                      ! GMR is global
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

  REAL,INTENT(IN)::  alt    ! geometric altitude, km.
  REAL,INTENT(OUT):: sigma  ! density/sea-level standard density             
  REAL,INTENT(OUT):: delta  ! pressure/sea-level standard pressure            
  REAL,INTENT(OUT):: theta  ! temperature/sea-level standard temperature   

  REAL:: h   ! geopotential altitude
!----------------------------------------------------------------------------
  h=alt*REARTH/(alt+REARTH)      ! convert geometric to geopotential altitude
                                                ! REARTH is a module variable
  IF (h < 11.0) THEN
    theta=1.0+(-6.5/288.15)*h                                   ! Troposphere
    delta=theta**(GMR/6.5)                         ! GMR is a module variable
  ELSE
    theta=216.65/288.15                                        ! Stratosphere
    delta=0.2233611*EXP(-GMR*(h-11.0)/216.65)
  END IF

  sigma=delta/theta
  RETURN
END Subroutine SimpleAtmosphere   ! -----------------------------------------

!+
FUNCTION DynamicViscosity(theta) RESULT(visc)
! ------------------------------------------------------------------------------
! PURPOSE - Compute viscosity using Sutherland's formula.
!        Returns viscosity in kg/(meter-sec)

  REAL,INTENT(IN) :: theta                ! temperature/sea-level temperature  
  REAL:: visc
  REAL:: temp                              ! temperature in deg Kelvin
  REAL,PARAMETER:: BETAVISC = 1.458E-6 ! viscosity term, N s/(sq.m sqrt(kelvins)
  REAL,PARAMETER:: SUTH = 110.4              ! Sutherland's constant, kelvins
  REAL,PARAMETER:: TZERO = 288.15            ! temperature at sealevel, deg K
!-------------------------------------------------------------------------------
  temp=TZERO*theta
  visc=BETAVISC*SQRT(temp*temp*temp)/(temp+SUTH)
  RETURN
END Function DynamicViscosity   ! ----------------------------------------------

END Module AtmosphereProcedures   ! ============================================

!+
MODULE TableMakingProcedures
! ------------------------------------------------------------------------------
! PURPOSE - Collect the procedures used in making the tables. Everything is to
!  be standard HTML 5.

USE AtmosphereProcedures, ONLY: Atmosphere, SimpleAtmosphere, DynamicViscosity, REARTH
IMPLICIT NONE

  CHARACTER(LEN=*),PARAMETER:: TD1="<td>", TD2="</td>"
  CHARACTER(LEN=*),PARAMETER:: TDC1='<td class="center">'

  REAL,PARAMETER:: PI = 3.14159265

  REAL,PARAMETER:: FT2METERS = 0.3048       ! mult. ft. to get meters (exact)
  REAL,PARAMETER:: KELVIN2RANKINE = 1.8           ! mult kelvins to get deg R
  REAL,PARAMETER:: PSF2NSM = 47.880258          ! mult lb/sq.ft to get N/sq.m
  REAL,PARAMETER:: SCF2KCM = 515.379         ! mult slug/cu.ft to get kg/cu.m
  REAL,PARAMETER:: TZERO = 288.15            ! temperature at sealevel, deg K
  REAL,PARAMETER:: PZERO = 101325.0            ! pressure at sealevel, N/sq.m
  REAL,PARAMETER:: RHOZERO = 1.2250            ! density at sealevel, kg/cu.m
  REAL,PARAMETER:: ASOUNDZERO = 340.294   ! speed of sound at sealevel, m/sec

  INTEGER,PARAMETER:: HTML=1

  INTEGER:: errCode

CONTAINS

!+
SUBROUTINE CreateHTMLfile()
! ------------------------------------------------------------------------------
! PURPOSE - Open the file 
  OPEN(UNIT=HTML, FILE='atmTabs.html', STATUS='REPLACE', &
    IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')
  IF (errCode/=0) THEN
    WRITE(*,*) " Unable to open atmTabs.html for output. Abort"
    STOP
  END IF
 
  WRITE(HTML,*) '<!DOCTYPE html>'
  WRITE(HTML,*) '<html lang="en">'
  WRITE(HTML,*) '<head>'
  WRITE(HTML,*) '<meta charset="utf-8" />'
  WRITE(HTML,*) '<title>Tables of U.S. Standard Atmosphere 1976</title>'

  WRITE(HTML,*) '<style type="text/css">'
  WRITE(HTML,*) 'body {margin: 5px; color: black; background: white;', &
    ' font-family: Arial,helvetica,sans-serif;}'
  WRITE(HTML,*) 'h1 {font-size: 200%;}'
  WRITE(HTML,*) 'h2 {font-size: 170%;}'

  WRITE(HTML,*) 'table {border-collapse: collapse; ', & 
    'black; cellspacing: 0; cellpadding: 2px;}'
  WRITE(HTML,*) 'td {border: 1px solid black; text-align: right; ', &
    'padding: 2px; font-size: 80%;}'
  WRITE(HTML,*) 'td.center {text-align:center}'
  WRITE(HTML,*) 'td.right {text-align:right}'
  WRITE(HTML,*) 'tr.gray {background-color: #f0f0f0;}'   ! displays vary. 
                                                  ! You may need to tweak this.
  
  WRITE(HTML,*) 'footer {font-size: 90%; border-top: 3px solid red; }'
  WRITE(HTML,*) '</style>'

  WRITE(HTML,*) '</head>'
  WRITE(HTML,*) '<body>'
  WRITE(HTML,*) '<h1>Tables of the U.S. Standard Atmosphere, 1976</h1>'
  
  WRITE(HTML,*) '<nav>'
  WRITE(HTML,*) '<a id="top"></a>'
  WRITE(HTML,*) '<br /><a href="#Table1">Table 1, US units to 280000 ft.</a>'
  WRITE(HTML,*) '<br /><a href="#Table2">Table 2, US units to 65000 ft.</a>'
  WRITE(HTML,*) '<br /><a href="#Table3">Table 3, SI units to 86 km.</a>'
  WRITE(HTML,*) '<br /><a href="#Table4">Table 4, SI units to 20 km.</a>'
  WRITE(HTML,*) '</nav>'

  RETURN
END Subroutine CreateHTMLfile   ! ==============================================

!+
SUBROUTINE CloseHTMLfile()
! ------------------------------------------------------------------------------
! PURPOSE - Write the appropriate HTML commands to complete the HTML statements
!  and then close the file.
!-------------------------------------------------------------------------------
  WRITE(HTML,*) "<footer>"
  WRITE(HTML,*) '<br />Created by program tablesHTML from '// &
    '<a href="http://www.pdas.com/index.html">'// &
    'Public Domain Aeronautical Software</a>.'
  WRITE(HTML,*) "</footer>"
  WRITE(HTML,*) " </body></html>"  
  CLOSE(UNIT=HTML)
  RETURN
END Subroutine CloseHTMLfile   ! -----------------------------------------------
 
!+
SUBROUTINE Table1()
! ---------------------------------------------------------------------------
! PURPOSE - Make an atmosphere table from 0 to 280000 ft at 5000 foot intervals.
!   Use procedure Atmosphere for calculations; print answers in US units

  INTEGER:: i                              ! altitude in thousand feet
  REAL:: altKm                                       !  altitude in Km
  REAL:: sigma,delta,theta       ! density,pressure,temperature ratios
  REAL:: temp,pressure,density
  REAL:: asound                             ! speed of sound in ft/sec
  REAL:: viscosity                       ! viscosity in slugs/(ft*sec)
  REAL:: kinematicViscosity                                ! sq.ft/sec
! ---------------------------------------------------------------------------
  WRITE(HTML,*) '<a id="Table1"> </a>'

  WRITE(HTML,*) "<h2>Atmosphere to 280 Kft by 5000 ft (US units)</h2>"
  WRITE(HTML,*) "<table>"
  WRITE(HTML,*) "<tr>"
  WRITE(HTML,*) TDC1//"alt"//TD2
  WRITE(HTML,*) TDC1//"sigma"//TD2
  WRITE(HTML,*) TDC1//"delta"//TD2
  WRITE(HTML,*) TDC1//"theta"//TD2
  WRITE(HTML,*) TDC1//"temp"//TD2
  WRITE(HTML,*) TDC1//"pressure"//TD2
  WRITE(HTML,*) TDC1//"density"//TD2
  WRITE(HTML,*) TDC1//"V sound"//TD2
  WRITE(HTML,*) TDC1//"dynamic<br /> viscosity"//TD2
  WRITE(HTML,*) TDC1//"kinematic<br /> viscosity"//TD2
  WRITE(HTML,*) "</tr>"

  WRITE(HTML,*) "<tr>"
  WRITE(HTML,*) TDC1//"Kft"//TD2
  WRITE(HTML,*) "<td> </td>"
  WRITE(HTML,*) "<td> </td>"
  WRITE(HTML,*) "<td> </td>"
  WRITE(HTML,*) TDC1//"deg R"//TD2
  WRITE(HTML,*) TDC1//"psf"//TD2
  WRITE(HTML,*) TDC1//"slug/ft<sup>3</sup>"//TD2
  WRITE(HTML,*) TDC1//"ft/sec"//TD2
  WRITE(HTML,*) TDC1//"slugs/ft-sec"//TD2
  WRITE(HTML,*) TDC1//"ft<sup>2</sup>/sec"
  WRITE(HTML,*) "</tr>"

  DO i=-1,56   ! vary alt (ft) from -5000 to 280,000 by 5000
    IF (MODULO(i,2)==0) THEN
      WRITE(HTML,*) '<tr>'
    ELSE
      WRITE(HTML,*) '<tr class="gray">'
    END IF

    altKm=5*i*FT2METERS
    WRITE(HTML,'(A,I4,A)') TD1, 5*i, TD2

    CALL Atmosphere(altKm, sigma,delta,theta)
    WRITE(HTML,'(A,ES12.4E1,A)') TD1, sigma, TD2
    WRITE(HTML,'(A,ES12.4E1,A)') TD1, delta, TD2
    WRITE(HTML,'(A,F7.4,A)') TD1, theta, TD2

    temp=KELVIN2RANKINE*TZERO*theta
    WRITE(HTML,'(A,F7.2,A)') TD1, temp, TD2

    pressure=(PZERO/PSF2NSM)*delta
    WRITE(HTML,'(A,ES12.4E1,A)') TD1, pressure, TD2

    density=(RHOZERO/SCF2KCM)*sigma
    WRITE(HTML,'(A,ES12.4E1,A)') TD1, density, TD2

    asound=(ASOUNDZERO/FT2METERS)*Sqrt(theta)
    WRITE(HTML,'(A,F9.2,A)') TD1, asound, TD2

    viscosity=(1.0/PSF2NSM)*DynamicViscosity(theta)
    WRITE(HTML,'(A,ES12.4E1,A)') TD1, viscosity, TD2

    kinematicViscosity=viscosity/density
    WRITE(HTML,'(A,ES12.4E1,A)') TD1, kinematicViscosity, TD2


    WRITE(HTML,*) "</tr>"
  END DO

  WRITE(HTML,*) "</table>"
  WRITE(HTML,*) '<br />Go to the <a href="#Top">top</a> of the page'
  
RETURN
 
END Subroutine Table1   ! ---------------------------------------------------


!+
SUBROUTINE Table2()
! ----------------------------------------------------------------------------
! PURPOSE - Make an atmosphere table from 0 to 65000 ft at 1000 foot intervals. 
!   Use Simple Atmosphere for calculations; print answers in US units.

  INTEGER:: i !     altitude in thousand feet
  REAL:: altKm  !
  REAL:: sigma,delta,theta
  REAL:: temp,pressure,density,asound
  REAL:: viscosity,kinematicViscosity,vratio
!----------------------------------------------------------------------------
  WRITE(HTML,*) '<a id="Table2"> </a>'

  WRITE(HTML,*) "<h2>Atmosphere to 65000 ft by 1000 ft (US units)</h2>"
  WRITE(HTML,*) "<table>"
  WRITE(HTML,*) "<tr>"
  WRITE(HTML,*) TDC1//"alt"//TD2
  WRITE(HTML,*) TDC1//"sigma"//TD2
  WRITE(HTML,*) TDC1//"delta"//TD2
  WRITE(HTML,*) TDC1//"theta"//TD2
  WRITE(HTML,*) TDC1//"temp"//TD2
  WRITE(HTML,*) TDC1//"pressure"//TD2
  WRITE(HTML,*) TDC1//"density"//TD2
  WRITE(HTML,*) TDC1//"V sound"//TD2
  WRITE(HTML,*) TDC1//"dynamic<br />viscosity"//TD2
  WRITE(HTML,*) TDC1//"kinematic<br />viscosity"//TD2
  WRITE(HTML,*) TDC1//"ratio"//TD2
  WRITE(HTML,*) "</tr>"

  WRITE(HTML,*) "<tr>"
  WRITE(HTML,*) TDC1//"kft"//TD2
  WRITE(HTML,*) "<td> </td>"
  WRITE(HTML,*) "<td> </td>"
  WRITE(HTML,*) "<td> </td>"
  WRITE(HTML,*) TDC1//"deg R"//TD2
  WRITE(HTML,*) TDC1//"psf"//TD2
  WRITE(HTML,*) TDC1//"slugs/ft<sup>3</sup>"//TD2
  WRITE(HTML,*) TDC1//"ft/sec"//TD2
  WRITE(HTML,*) TDC1//"slugs/ft-sec"//TD2
  WRITE(HTML,*) TDC1//"ft<sup>2</sup>/sec"//TD2
  WRITE(HTML,*) TDC1//"ft<sup>-1</sup>"//TD2
  WRITE(HTML,*) "</tr>"

  DO i=-1,65   ! vary alt (ft) from -1000 to 65000 by 1000
    IF (MODULO(i,2)==0) THEN
      WRITE(HTML,*) '<tr>'
    ELSE
      WRITE(HTML,*) '<tr class="gray">'
    END IF

    altKm=i*FT2METERS
    WRITE(HTML,'(A,I4,A)') TD1, i, TD2

    CALL SimpleAtmosphere(altKm, sigma,delta,theta)
    WRITE(HTML,'(A,F7.4,A)') TD1, sigma, TD2
    WRITE(HTML,'(A,F7.4,A)') TD1, delta, TD2
    WRITE(HTML,'(A,F7.4,A)') TD1, theta, TD2

    temp=(KELVIN2RANKINE*TZERO)*theta
    WRITE(HTML,'(A,F7.2,A)') TD1, temp, TD2

    pressure=(PZERO/PSF2NSM)*delta
    WRITE(HTML,'(A,ES12.4E1,A)') TD1, pressure, TD2

    density=(RHOZERO/SCF2KCM)*sigma
    WRITE(HTML,'(A,ES12.4E1,A)') TD1, density, TD2

    asound=(ASOUNDZERO/FT2METERS)*Sqrt(theta)
    WRITE(HTML,'(A,F9.2,A)') TD1, asound, TD2

    viscosity=(1.0/PSF2NSM)*DynamicViscosity(theta)
    WRITE(HTML,'(A,ES12.4E1,A)') TD1, viscosity, TD2

    kinematicViscosity=viscosity/density
    WRITE(HTML,'(A,ES12.4E1,A)') TD1, kinematicViscosity, TD2

    vratio=1.0E-6*asound/kinematicViscosity
    WRITE(HTML,'(A,F9.3,A)') TD1, vratio, TD2
    WRITE(HTML,*) "</tr>"
  END DO

  WRITE(HTML,*) "</table>"
  WRITE(HTML,*) '<br>go to the <a href="#Top">top</a> of the page'

  RETURN
END Subroutine Table2   ! ---------------------------------------------------

!+
SUBROUTINE Table3()
! ---------------------------------------------------------------------------
! PURPOSE - Make an atmosphere table from 0 to 86 km at 2 km intervals.
!   Use procedure Atmosphere for calculations; print answers in SI units

  INTEGER:: i
  REAL:: altKm                                ! altitude in kilometers
  REAL:: sigma,delta,theta
  REAL:: temp,pressure,density, asound
  REAL:: viscosity,kinematicViscosity
!----------------------------------------------------------------------------
  WRITE(HTML,*) '<a id="Table3"> </a>'

  WRITE(HTML,*) "<h2>Atmosphere to 86 Km by 2 Km (SI units)</h2>"
  WRITE(HTML,*) "<table>"
  WRITE(HTML,*) "<tr>"
  WRITE(HTML,*) TDC1//"alt"//TD2
  WRITE(HTML,*) TDC1//"sigma"//TD2
  WRITE(HTML,*) TDC1//"delta"//TD2
  WRITE(HTML,*) TDC1//"theta"//TD2
  WRITE(HTML,*) TDC1//"temp"//TD2
  WRITE(HTML,*) TDC1//"pressure"//TD2
  WRITE(HTML,*) TDC1//"density"//TD2
  WRITE(HTML,*) TDC1//"V sound"//TD2
  WRITE(HTML,*) TDC1//"dynamic<br /> viscosity"//TD2
  WRITE(HTML,*) TDC1//"kinematic<br /> viscosity"//TD2
  WRITE(HTML,*) "</tr>"

  WRITE(HTML,*) "<tr>"
  WRITE(HTML,*) TDC1//"km"//TD2
  WRITE(HTML,*) "<td> </td>"
  WRITE(HTML,*) "<td> </td>"
  WRITE(HTML,*) "<td> </td>"
  WRITE(HTML,*) TDC1//"K"//TD2
  WRITE(HTML,*) TDC1//"N/m<sup>2</sup>"//TD2
  WRITE(HTML,*) TDC1//"kg/m<sup>3</sup>"//TD2
  WRITE(HTML,*) TDC1//"m/s"//TD2
  WRITE(HTML,*) TDC1//"kg-m<sup>-1</sup>-s<sup>-1</sup>"//TD2
  WRITE(HTML,*) TDC1//"m<sup>2</sup> s<sup>-1</sup>"//TD2
  WRITE(HTML,*) "</tr>"

  DO i=-2,86,2
    IF (MODULO(i,4)==0) THEN
      WRITE(HTML,*) '<tr>'
    ELSE
      WRITE(HTML,*) '<tr class="gray">'
    END IF

    altKm=REAL(i)
    WRITE(HTML,'(A,I4,A)') TD1, i, TD2

    CALL Atmosphere(altKm, sigma,delta,theta)
    WRITE(HTML,'(A,ES12.4E1,A)') TD1, sigma, TD2
    WRITE(HTML,'(A,ES12.4E1,A)') TD1, delta, TD2
    WRITE(HTML,'(A,F9.4,A)') TD1, theta, TD2

    temp=TZERO*theta
    WRITE(HTML,'(A,F7.2,A)') TD1, temp, TD2

    pressure=PZERO*delta
    WRITE(HTML,'(A,ES12.4E1,A)') TD1, pressure, TD2

    density=RHOZERO*sigma
    WRITE(HTML,'(A,ES12.4E1,A)') TD1, density, TD2

    asound=ASOUNDZERO*Sqrt(theta)
    WRITE(HTML,'(A,F9.2,A)') TD1, asound, TD2

    viscosity=DynamicViscosity(theta)
    WRITE(HTML,'(A,ES12.4E1,A)') TD1, viscosity, TD2

    kinematicViscosity=viscosity/density
    WRITE(HTML,'(A,ES12.4E1,A)') TD1, kinematicViscosity, TD2


    WRITE(HTML,*) "</tr>"

  END DO

  WRITE(HTML,*) "</table>"
  WRITE(HTML,*) '<br />Go to the <a href="#Top">top</a> of the page'

  RETURN
END Subroutine Table3   ! ---------------------------------------------------

!+
SUBROUTINE Table4()
! ---------------------------------------------------------------------------
! PURPOSE - Make an atmosphere table from 0 to 20 km at 0.5 km intervals.
!   Use Simple Atmosphere for calculations; print answers in SI units

  INTEGER:: i                                   !  steps thru altitude
  REAL:: altKm
  REAL:: sigma,delta,theta
  REAL:: temp,pressure,density, asound
  REAL:: viscosity,kinematicViscosity,vratio
!----------------------------------------------------------------------------
  WRITE(HTML,*) '<a id="Table4"></a>'

  WRITE(HTML,*) "<h2>Atmosphere to 20 Km by 0.5 Km (SI units)</h2>"
  WRITE(HTML,*) "<table>"
  WRITE(HTML,*) "<tr>"     
  WRITE(HTML,*) TDC1//'alt'//TD2
  WRITE(HTML,*) TDC1//'sigma'//TD2
  WRITE(HTML,*) TDC1//'delta'//TD2
  WRITE(HTML,*) TDC1//'theta'//TD2
  WRITE(HTML,*) TDC1//'temp'//TD2
  WRITE(HTML,*) TDC1//'pressure'//TD2
  WRITE(HTML,*) TDC1//'density'//TD2
  WRITE(HTML,*) TDC1//'v sound'//TD2
  WRITE(HTML,*) TDC1//'dynamic<br />viscosity'//TD2
  WRITE(HTML,*) TDC1//'kinematic<br /> viscosity'//TD2
  WRITE(HTML,*) TDC1//'a/nu'//TD2
  WRITE(HTML,*) "</tr>"

  WRITE(HTML,*) "<tr>"
  WRITE(HTML,*) TDC1//'km'//TD2
  WRITE(HTML,*) "<td> </td>"
  WRITE(HTML,*) "<td> </td>"
  WRITE(HTML,*) "<td> </td>"
  WRITE(HTML,*) TDC1//'K'//TD2
  WRITE(HTML,*) TDC1//'N/m<sup>2</sup>'//TD2
  WRITE(HTML,*) TDC1//'kg/m<sup>3</sup>'//TD2
  WRITE(HTML,*) TDC1//'m/s'//TD2
  WRITE(HTML,*) TDC1//'N s m<sup>-2</sup>'//TD2
  WRITE(HTML,*) TDC1//'m<sup>2</sup> s<sup>-1</sup>'//TD2
  WRITE(HTML,*) TDC1//"m<sup>-1</sup>"//TD2
  WRITE(HTML,*) "</tr>"

  DO i=-1,40
    IF (MODULO(i,2)==0) THEN
      WRITE(HTML,*) '<tr>'
    ELSE
      WRITE(HTML,*) '<tr class="gray">'
    END IF
    altKm=0.5*i
    WRITE(HTML,'(A,F0.1,A)') TD1, altKm, TD2

    CALL SimpleAtmosphere(altKm, sigma,delta,theta)
    WRITE(HTML,'(A,F0.4,A)') TD1, sigma, TD2
    WRITE(HTML,'(A,F0.4,A)') TD1, delta, TD2
    WRITE(HTML,'(A,F0.4,A)') TD1, theta, TD2

    temp=TZERO*theta
    WRITE(HTML,'(A,F0.2,A)') TD1, temp, TD2

    pressure=PZERO*delta
    WRITE(HTML,'(A,ES12.4E1,A)') TD1, pressure, TD2

    density=RHOZERO*sigma
    WRITE(HTML,'(A,ES12.4E1,A)') TD1, density, TD2

    asound=ASOUNDZERO*Sqrt(theta)
    WRITE(HTML,'(A,F0.2,A)') TD1, asound, TD2

    viscosity=DynamicViscosity(theta)
    WRITE(HTML,'(A,ES12.4E1,A)') TD1, viscosity, TD2

    kinematicViscosity=viscosity/density
    WRITE(HTML,'(A,ES12.4E1,A)') TD1, kinematicViscosity, TD2

    vratio=1E-6*asound/kinematicViscosity
    WRITE(HTML,'(A,F0.2,A)') TD1, vratio, TD2

    WRITE(HTML,*) "</tr>"
  END DO
  WRITE(HTML,*) "</table>"
  WRITE(HTML,*) '<br>Go to the <a href="#Top">top</a> of the page'

  RETURN 
END Subroutine Table4   ! ---------------------------------------------------

END Module TableMakingProcedures   ! ========================================


!+
PROGRAM TablesHTML
! ---------------------------------------------------------------------------
! PURPOSE - 

USE TableMakingProcedures,ONLY: CreateHTMLfile,CloseHTMLfile, &
  Table1,Table2,Table3,Table4
USE ISO_FORTRAN_ENV
IMPLICIT NONE

  CHARACTER(LEN=*),PARAMETER::GREETING =  &
        ' tables76 - A Fortran program to make atmosphere tables with HTML'
  CHARACTER(LEN=*),PARAMETER::AUTHOR=  &
    ' Ralph L. Carmichael, Public Domain Aeronautical Software'
  CHARACTER(LEN=*),PARAMETER::MODIFIER=' '   ! put your name here if you mod it
  CHARACTER(LEN=*),PARAMETER::VERSION=' 2.21 (17 July 2012)'
  CHARACTER(LEN=*),PARAMETER::FAREWELL= &
    'The file atmTabs.html has been added to your directory.'

!----------------------------------------------------------------------------
  WRITE(*,*) GREETING
  WRITE(*,*) AUTHOR
  IF (MODIFIER /= ' ') WRITE(*,*) 'Modified by '//MODIFIER
  WRITE(*,*) VERSION
!  WRITE(*,'(2A/2A)') ' This file was compiled by ', compiler_version(), &
!                     ' using the options ', compiler_options()

  CALL CreateHTMLfile()

  CALL Table1()
  CALL Table2()
  CALL Table3()
  CALL Table4()
  
  CALL CloseHTMLfile()


  WRITE (*,*) FAREWELL
  STOP
END Program TablesHTML   ! --------------------------------------------------
