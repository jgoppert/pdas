!+
! PROGRAM BaseBall                                              ! \atmos\bb.f90
! ---------------------------------------------------------------------------
! PURPOSE - Calculate the trajectory of a baseball with air drag
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software

! REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
!  7Aug00  1.0   RLC   Initial release
! 22Dec00  1.1   RLC   Removed all graphics for maximum portability
! 04Mar13  1.2   RLC   Small fixes, comments

!+
MODULE BaseBallAuxiliary
! ---------------------------------------------------------------------------
! PURPOSE - Collect all of the functions and subroutines used by the main
!  program. This guarantees explicit interfaces. Subroutine Trajectory is 
!  called from the main program and all others called from Trajectory.
!  All are declared Private except Trajectory, which is Public.
IMPLICIT NONE

  TYPE:: TrajectoryPoint
    REAL:: time
    REAL,DIMENSION(2):: position, velocity, acceleration
  END TYPE TrajectoryPoint

! MODULE CONSTANTS:
  REAL,PRIVATE,PARAMETER:: PI=3.14159265
  REAL,PRIVATE,PARAMETER:: FT2METERS = 0.3048   ! convert feet to meters
  REAL,PRIVATE,PARAMETER:: G = 9.8066   ! acceleration of gravity
  REAL,PRIVATE,PARAMETER:: LBS2KG = 0.45359237  ! weight in pounds to mass in kg.
  REAL,PRIVATE,PARAMETER:: RHOZERO = 1.2250    ! density of air at sealevel, kg/cu.m

! BASEBALL CONSTANTS - by the rules, the circumference of a baseball must
! be no less than 9 inches and no more than 9.25 inches. The weight of a
! baseball must be no less than 5.00 ounces and no more than 5.25 ounces.
! I assume my ideal ball lies exactly in the middle.
! In SI units:
  REAL,PARAMETER:: DIAM=(9.125/(12*PI))*FT2METERS ! diameter of a baseball (m)
  REAL,PARAMETER:: MASS=(5.125/16)*LBS2KG     ! mass of a baseball (kg)
  REAL,PARAMETER:: SREF = 0.25*PI*DIAM*DIAM   ! frontal area (sq.m)


  PRIVATE:: Accel
  PRIVATE:: BaseballKutta 
  PRIVATE:: CorrectFinalPosition
  PRIVATE:: Viscosity
  PRIVATE:: SimpleAtmosphere
  PUBLIC :: Trajectory

CONTAINS

!+
FUNCTION Accel(time, position, velocity) RESULT(acceleration)      ! not PURE
! ---------------------------------------------------------------------------
! PURPOSE - Compute the acceleration (vector) for a spherical projectile
!    moving through a viscous medium. Assume Mach number is small enough
!    that wave drag may be neglected. Ignore added mass term.
! NOTE - position has units of meters, but first argument to SimpleAtmosphere
!   is in kilometers. Be sure to remember to multiply by 0.001

  REAL,INTENT(IN):: time
  REAL,INTENT(IN),DIMENSION(:):: position, velocity

  REAL,DIMENSION(2):: acceleration

!  CHARACTER(LEN=*),PARAMETER:: FMT = &
!    "(2F7.2, F9.2, 2F7.2, F7.0, ES9.2E1, F7.4,6F7.2)"
  REAL,PARAMETER,DIMENSION(2):: VERTICAL = (/0.0, 1.0/)

  REAL:: cd
  REAL:: density
  REAL,DIMENSION(2):: drag
  REAL:: dragMagnitude
  REAL:: q      ! dynamic pressure
  REAL:: reynolds
  REAL:: sigma,delta,theta
  REAL,DIMENSION(2):: unitVelocity
  REAL:: vmag
  REAL:: vsq
!----------------------------------------------------------------------------
  vsq=SUM(velocity**2)
  vmag=SQRT(vsq)
  unitVelocity=velocity/vmag
  CALL SimpleAtmosphere(0.001*position(2), sigma,delta,theta )
                                        ! first arg is altitude in kilometers
  density=sigma*RHOZERO
  q=0.5*density*vsq
  reynolds=density*vmag*DIAM/Viscosity(theta)
  cd=CdSphere(reynolds)
  dragMagnitude=cd*q*SREF
  drag=-dragMagnitude*unitVelocity
  acceleration=drag/MASS - G*VERTICAL

!  WRITE(DBG, FMT) time, position, &
!     velocity, q, reynolds, cd, dragMagnitude, drag, acceleration/G
  RETURN
END Function Accel   ! ------------------------------------------------------

!+
FUNCTION CdSphere(r) RESULT(d)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the drag coefficient of a sphere as a function of
!   Reynolds number. Assumes Mach number is small.
!   Taken from Chow, "Computational Aerodynamics"
  REAL,INTENT(IN):: r   ! Reynolds number
  REAL:: d              ! drag coefficient based on cross-section area
!----------------------------------------------------------------------------
  IF (r <= 0) THEN
    d=0
  ELSE IF (r <= 1.0) THEN
    d=24.0/r
  ELSE IF (r <= 400.0) THEN
    d=24*r**(-0.646)
  ELSE IF (r <= 3E5) THEN
    d=0.5
  ELSE IF (r <= 2E6) THEN
    d=3.66E-4*r**0.4275
  ELSE
    d=0.18
  END IF

  RETURN
END Function CdSphere   ! ---------------------------------------------------

!+
SUBROUTINE CorrectFinalPosition(initialAltitude, a1, a2)
! ---------------------------------------------------------------------------
! PURPOSE - Correct the final point of a trajectory so that its  y-coordinate
!  is exactly equal to initial altitude. Assume that the altitudes of a1
!  and a2 straddle initialAltitude.

  REAL,INTENT(IN):: initialAltitude
  TYPE(TrajectoryPoint),INTENT(IN):: a1
  TYPE(TrajectoryPoint),INTENT(IN OUT):: a2

  REAL:: fraction
!----------------------------------------------------------------------------
  fraction=(initialAltitude-a1%position(2))/(a2%position(2)-a1%position(2))
  a2%time    =a1%time     + fraction*(a2%time    -a1%time)
  a2%position=a1%position + fraction*(a2%position-a1%position)
  a2%velocity=a1%velocity + fraction*(a2%velocity-a1%velocity)
  a2%acceleration=a1%acceleration + &
                            fraction*(a2%acceleration-a1%acceleration)

  RETURN
END Subroutine CorrectFinalPosition   ! -------------------------------------

!+
FUNCTION Viscosity(theta) RESULT(visc)
!   -------------------------------------------------------------------------
! PURPOSE - Compute viscosity using Sutherland's formula.
!        Returns viscosity in kg/(meter-sec)

  REAL,INTENT(IN) :: theta                ! temperature/sea-level temperature  
  REAL:: visc

  REAL,PARAMETER:: BETAVISC = 1.458E-6 ! viscosity term, N sec/(sq.m sqrt(deg K)
  REAL,PARAMETER:: SUTH = 110.4              ! Sutherland's constant, deg K
  REAL,PARAMETER:: TZERO = 288.15            ! temperature at sealevel, deg K

  REAL:: temp                              ! temperature in deg Kelvin
!----------------------------------------------------------------------------
  temp=TZERO*theta
  visc=BETAVISC*Sqrt(temp*temp*temp)/(temp+SUTH)
  RETURN
END Function Viscosity   ! --------------------------------------------------

!+
SUBROUTINE SimpleAtmosphere(alt,sigma,delta,theta)
!   -------------------------------------------------------------------------
! PURPOSE - Compute the characteristics of the lower atmosphere.

! NOTES-Correct to 20 km. Only approximate above there

  REAL,INTENT(IN)::  alt    ! geometric altitude, km.
  REAL,INTENT(OUT):: sigma  ! density/sea-level standard density             
  REAL,INTENT(OUT):: delta  ! pressure/sea-level standard pressure            
  REAL,INTENT(OUT):: theta  ! temperature/sea-level standard temperature   

  REAL,PARAMETER:: REARTH = 6369.0                ! radius of the Earth (km)
  REAL,PARAMETER:: GMR = 34.163195                            ! gas constant

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
SUBROUTINE BaseBallKutta(p1, h, p2)
! ---------------------------------------------------------------------------
! PURPOSE - Advance one time-like step in a trajectory. This is a system
!  of four first order ordinary differential equations. Use fourth-order
!  Runge-Kutta equation to advance one time step.

  TYPE(TrajectoryPoint),INTENT(IN):: p1   ! current position
  REAL,INTENT(IN):: h                     ! step to be taken in time
  TYPE(TrajectoryPoint),INTENT(OUT):: p2  ! next position

  REAL,PARAMETER:: HALF=0.5, SIXTH=1.0/6.0

  REAL:: t
  REAL,DIMENSION(2):: x,v,a
  REAL,DIMENSION(2):: dx1,dx2,dx3,dx4
  REAL,DIMENSION(2):: dv1,dv2,dv3,dv4
!----------------------------------------------------------------------------
  t=p1%time
  x=p1%position
  v=p1%velocity

  a=Accel(t, x, v)
  dx1 = h*v
  dv1 = h*a 

  a=Accel(t+HALF*h, x+HALF*dx1, v+HALF*dv1)
  dx2 = h*(v + HALF*dv1)
  dv2 = h*a
  
  a=Accel(t+HALF*h, x+HALF*dx2, v+HALF*dv2)
  dx3 = h*(v + HALF*dv2)
  dv3 = h*a 

  a=Accel(t+h, x+dx3, v+dv3)
  dx4 = h*(v + dv3)   ! oops, was dx3.   Thanks, Craig Bushman
  dv4 = h*a

  p2%time = t + h
  p2%position = p1%position + SIXTH*(dx1+dx2+dx2+dx3+dx3+dx4)
  p2%velocity = p1%velocity + SIXTH*(dv1+dv2+dv2+dv3+dv3+dv4)

  p2%acceleration=Accel(p2%time, p2%position, p2%velocity)

  RETURN
END SUBROUTINE BaseBallKutta   ! --------------------------------------------

!+
SUBROUTINE Trajectory(initialAltitude, initialVelocity, initialTheta, &
  dt, normalized, npts, history)
! ---------------------------------------------------------------------------
! PURPOSE - Compute a trajectory, performing numerical solution of a set of
!   ordinary differential equations with a fixed time step. Halt the
!   calculation when the altitude is less than the initial altitude and
!   correct the final point to have the same altitude as the initial altitude.

  REAL,INTENT(IN):: initialAltitude, initialVelocity, initialTheta
  REAL,INTENT(IN):: dt
  LOGICAL,INTENT(IN):: normalized
  INTEGER,INTENT(OUT):: npts
  TYPE(TrajectoryPoint),INTENT(OUT),DIMENSION(:):: history

  INTEGER:: k
  REAL:: t

  REAL,DIMENSION(2):: position,velocity,acceleration

!----------------------------------------------------------------------------
  t=0.0
  position(1)=0.0
  position(2)=initialAltitude
  velocity(1)=initialVelocity*COS(initialTheta*PI/180)
  velocity(2)=initialVelocity*SIN(initialTheta*PI/180)
  acceleration=Accel(t, position, velocity)  
  history(1)%time=t
  history(1)%position=position
  history(1)%velocity=velocity
  history(1)%acceleration=acceleration

  DO k=2,SIZE(history)
    CALL BaseballKutta(history(k-1), dt, history(k))
    IF (history(k)%position(2) < initialAltitude) EXIT
  END DO


  CALL CorrectFinalPosition(initialAltitude, history(k-1), history(k) )

  IF (normalized) THEN
    history(1:k)%position(2)=history(1:k)%position(2)-history(1)%position(2)
  END IF

  npts=k
  RETURN
END Subroutine Trajectory   ! -----------------------------------------------

!+
SUBROUTINE VacuumTrajectory(initialVelocity, initialTheta, dt, npts, history)
! ---------------------------------------------------------------------------
! PURPOSE - Same as subroutine Trajectory, but for the case where there is
!   no air resistance.

  REAL,INTENT(IN):: initialVelocity, initialTheta
  REAL,INTENT(IN):: dt
  INTEGER,INTENT(OUT):: npts
  TYPE(TrajectoryPoint),INTENT(OUT),DIMENSION(:):: history

  INTEGER:: k
  REAL:: t
  REAL:: vZeroX, vZeroZ

!  CHARACTER(LEN=30):: identifier
!----------------------------------------------------------------------------
  vZeroX=initialVelocity*COS(initialTheta*PI/180)
  vZeroZ=initialVelocity*SIN(initialTheta*PI/180)

  t=0.0
  DO k=1,SIZE(history)
    history(k)%time=t
    history(k)%position(1)=vZeroX*t
    history(k)%position(2)=t*(vZeroZ - 0.5*G*t)
    history(k)%velocity(1)=vZeroX
    history(k)%velocity(2)=vZeroZ-G*t
    history(k)%acceleration(1)=0.0  ! not used
    history(k)%acceleration(2)=-G   ! not used
    IF (history(k)%position(2) < 0.0 ) EXIT
    t=t+DT
  END DO

  CALL CorrectFinalPosition(0.0, history(k-1), history(k) )
  npts=k
  RETURN
END Subroutine VacuumTrajectory   ! -----------------------------------------

END MODULE BaseBallAuxiliary   ! ============================================


!+
PROGRAM BaseBall                                              ! \atmos\bb.f90
! ---------------------------------------------------------------------------
! PURPOSE - Calculate the trajectory of a baseball with air drag
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software


USE BaseBallAuxiliary
IMPLICIT NONE

  REAL,PARAMETER:: DT=0.1   ! time step
  TYPE(TrajectoryPoint),DIMENSION(1000):: history
  REAL:: initialAltitude  ! meters
  REAL:: initialAngle     ! degrees from horizontal
  REAL:: initialVelocity  ! meters per second
  INTEGER:: k,n
  LOGICAL,PARAMETER:: NORMALIZED=.TRUE.

!----------------------------------------------------------------------------
  WRITE(*,*) "Baseball trajectory from Public Domain Aeronautical Software."
  WRITE(*,*) "Check the web site at www.pdas.com"

! Instructions for use of this program: put your values for the three
!  initial conditions on the next three lines. Recompile. Run bb.
! Hints: Denver=1609  Mexico City=2420  La Paz=3650  Everest=8850
  initialAltitude = 1609.0  ! meters
  initialAngle = 40.0        ! deg
  initialVelocity = 35.0     ! m/s

  WRITE(*,*) &
   "     t          x          y         vx          vy         ax        ay"
  CALL Trajectory(initialAltitude, initialVelocity, initialAngle, &
      DT, NORMALIZED, n, history)
  DO k=1,n
    WRITE(*, '(F9.2,6F11.2)' ) history(k)
  END DO
  WRITE(*,*) "End of BaseBall"
  STOP
End Program BaseBall   ! ====================================================

