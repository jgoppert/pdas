!+
! PROGRAM ThreeView                                          ! \3view\3view.f90 
! ---------------------------------------------------------------------------
! PURPOSE - Plot three orthogonal views of an object defined in LaWGS
!   (Langley Wireframe Geometry Standard) format.
! AUTHOR  - Ralph L. Carmichael, Public Domain Aeronautical Software
! REVISION HISTORY                                                        
!   DATE  VERS PERSON  STATEMENT OF CHANGES                               
! 19Sep95  0.1   RLC   Converted old program ( C language)              
! 25Sep95  0.2   RLC   Got most bugs out                                
! 30Nov95  1.0   RLC   Made Fig have STATUS='DELETE'                    
! 18Dec02  1.5   RLC   Complete rewrite of the imaging code.
! 28Dec02  1.6   RLC   Fixed several errors in the imaging logic
! 02Dec08  1.7   RLC   Read file name from command line
! 14Dec08  1.8   RLC   Moved procedures to module
! 27Jan13  1.9   RLC   Omitted icomment from CopyData; also generalTitle



!+
MODULE NetworkProcedures
! ---------------------------------------------------------------------------
! PURPOSE - Define the derived type Network and collect the procedures that
!  operate on data with the type Network. 

! AUTHORS  - Ralph L. Carmichael, Public Domain Aeronautical Software
!            Charlotte Craidon, NASA Langley, coordinator of the Langley
!               wireframe geometry standard as documented in NASA TM 85767. 

! REVISION HISTORY                                                  
!   DATE  VERS PERSON  STATEMENT OF CHANGES   
!   1977         CC    Published NASA TM 85767
!   1989   0.1   RLC   Various routines to read networks in LaWgs format
! 23Nov02  0.5   RLC   Rewrote completely (Fortran 95)


IMPLICIT NONE

  TYPE:: Network
    CHARACTER(LEN=80):: name
    INTEGER:: id
    INTEGER:: rows,cols
    INTEGER:: symLocal, symGlobal
    INTEGER:: localImage  ! index of original network; =0 if an original
    INTEGER:: globalImage ! index of original network; =0 if an original
    REAL,POINTER,DIMENSION(:,:):: x,y,z
    REAL,DIMENSION(3):: rotate,translate,scale
  END TYPE Network
  

  INTEGER,PARAMETER,PRIVATE:: DBG = 3
  CHARACTER(LEN=*),PARAMETER:: NETWORK_VERSION = '0.64 (18Dec02)'
  
  PRIVATE:: BuildRotationMatrix
  PUBLIC::  CreateGlobalImageNetworks
  PUBLIC::  CreateLocalImageNetworks  
  PUBLIC::  DeallocateNetworks
  PUBLIC::  DumpNetworkGridPoints
  PUBLIC::  ReadOneNetwork
  PUBLIC::  ScanNetworks
  PRIVATE:: StripFortranQuotes
  PUBLIC:: TransformNetworks
! ---------------------------------------------------------------------------

CONTAINS

!+
FUNCTION BuildRotationMatrix(r) RESULT(rot)
! ---------------------------------------------------------------------------
! PURPOSE - Form the rotation matrix used for geometrical transformations
!  Taken from NASA TM 85767 defining LaWGS

  REAL,INTENT(IN),DIMENSION(:):: r   ! rotation angles, degrees
  REAL,DIMENSION(3,3):: rot

  REAL,PARAMETER:: PI=3.14159265, RAD=PI/180

  REAL:: cPhi,cTheta,cPsi
  REAL:: sPhi,sTheta,sPsi
  INTRINSIC:: COS,SIN
!----------------------------------------------------------------------------
  cPhi=COS(RAD*r(1))
  sPhi=SIN(RAD*r(1))
  cTheta=COS(RAD*r(2))
  sTheta=SIN(RAD*r(2))
  cPsi=COS(RAD*r(3))
  sPsi=SIN(RAD*r(3))

  rot(1,1)= cTheta*cPsi
  rot(2,1)= cTheta*sPsi
  rot(3,1)=-sTheta

  rot(1,2)=-sPsi*cPhi + sTheta*cPsi*sPhi
  rot(2,2)= cPsi*cPhi + sTheta*sPsi*sPhi
  rot(3,2)= cTheta*sPhi

  rot(1,3)= sPsi*sPhi + sTheta*cPsi*cPhi
  rot(2,3)=-cPsi*sPhi + sTheta*sPsi*cPhi
  rot(3,3)= cTheta*cPhi

  RETURN
END FUNCTION BuildRotationMatrix   ! ----------------------------------------

!+
SUBROUTINE CreateGlobalImageNetworks(n,a)
! ---------------------------------------------------------------------------
! PURPOSE - If any networks in the original file have the symGlobal flag = 1,
!  make an image network of it. Note that n and a will be modified and that a
!  must be dimensioned large enough to hold the new networks. Fatal if not.

  INTEGER,INTENT(IN OUT):: n
  TYPE(Network),INTENT(IN OUT),DIMENSION(:):: a 
 
  INTEGER:: i,j,k
  INTEGER:: nOriginal
  INTRINSIC:: SIZE, TRIM
! ---------------------------------------------------------------------------
  nOriginal=n
  DO k=1,nOriginal
    IF (a(k)%symGlobal==0) CYCLE
    IF (n+1 > SIZE(a)) THEN
      WRITE(*,*) 'Too many networks (2)'
      STOP
    END IF
    n=n+1
    a(n)%name=TRIM('global image of '//a(k)%name)
    a(n)%id=10000+a(k)%id
    a(n)%rows=a(k)%rows
    a(n)%cols=a(k)%cols
    a(n)%symLocal=a(k)%symLocal   ! but it won't be used
    a(n)%symGlobal=1   ! but it won't be used
    a(n)%localImage=0 
    a(n)%globalImage=k  ! keep track of where it came from
    i=a(n)%rows
    j=a(n)%cols
    ALLOCATE(a(n)%x(i,j), a(n)%y(i,j), a(n)%z(i,j) )
    a(n)%x=a(k)%x
    a(n)%y=-a(k)%y   ! this is where we do the imaging
    a(n)%z=a(k)%z
  END DO
  RETURN
END Subroutine CreateGlobalImageNetworks   ! ---------------------------------

!+
SUBROUTINE CreateLocalImageNetworks(n,a)
! ---------------------------------------------------------------------------
! PURPOSE - If any networks in the original file have the symLocal flag = 1,
!  make an image network of it. Note that n and a will be modified and that a
!  must be dimensioned large enough to hold the new networks. Fatal if not.

  INTEGER,INTENT(IN OUT):: n
  TYPE(Network),INTENT(IN OUT),DIMENSION(:):: a 
 
  INTEGER:: i,j,k
  INTEGER:: nOriginal
  INTRINSIC:: SIZE, TRIM
! ---------------------------------------------------------------------------
  nOriginal=n
  DO k=1,nOriginal
    IF (a(k)%symLocal==0) CYCLE
    IF (n+1 > SIZE(a)) THEN
      WRITE(*,*) 'Too many networks (1)'
      STOP
    END IF
    n=n+1
    a(n)%name=TRIM('local image of '//a(k)%name)
    a(n)%id=1000+a(k)%id
    a(n)%rows=a(k)%rows
    a(n)%cols=a(k)%cols
    a(n)%symLocal=1   ! but it won't be used
    a(n)%symGlobal=a(k)%symGlobal   ! we still may need to image it again!
    a(n)%localImage=k   ! keep track of where it came from
    a(n)%globalImage=0   ! will be set by CreateGlobalImageNetworks
    i=a(n)%rows
    j=a(n)%cols
    ALLOCATE(a(n)%x(i,j), a(n)%y(i,j), a(n)%z(i,j) )
    a(n)%x=a(k)%x
    a(n)%y=-a(k)%y   ! this is where we do the imaging
    a(n)%z=a(k)%z
    a(n)%rotate=a(k)%rotate
    a(n)%translate=a(k)%translate
    a(n)%scale=a(k)%scale
  END DO
  RETURN
END Subroutine CreateLocalImageNetworks   ! ---------------------------------

!+
SUBROUTINE DeAllocateNetworks(a)
! ---------------------------------------------------------------------------
! PURPOSE - Deallocate the memory allocated for gridpoints in each network
  TYPE(Network),INTENT(IN OUT),DIMENSION(:):: a
  
  INTEGER:: k
  INTRINSIC:: SIZE
!----------------------------------------------------------------------------
  DO k=1,SIZE(a)
    DEALLOCATE(a(k)%x, a(k)%y, a(k)%z)
  END DO
  RETURN
END Subroutine DeAllocateNetworks   ! ---------------------------------------

!+
SUBROUTINE DumpNetworkGridPoints(efu, a)
! ---------------------------------------------------------------------------
! PURPOSE - Print the defined networks
  INTEGER,INTENT(IN):: efu
  TYPE(Network),INTENT(IN),DIMENSION(:):: a
  
  CHARACTER(LEN=*),PARAMETER:: FMT = '(2I4,3F12.6)'
  INTEGER:: i,j,n
  INTRINSIC:: SIZE
!----------------------------------------------------------------------------
  DO n=1,SIZE(a)
    WRITE(efu,*) 'NETWORK',n
    DO j=1,a(n)%cols
      DO i=1,a(n)%rows
        WRITE(efu,FMT) i,j, a(n)%x(i,j), a(n)%y(i,j), a(n)%z(i,j)
      END DO
    END DO    
    WRITE(efu,*)   ! blank line
  END DO
  RETURN
END Subroutine DumpNetworkGridPoints   ! ------------------------------------

!+
SUBROUTINE GetXYZlimits(a, xmin,xmax, ymin,ymax, zmin,zmax)
! ---------------------------------------------------------------------------
! PURPOSE - Find the smallest box that will contain all the gridpoints in an
!  array of networks
  TYPE(Network),INTENT(IN),DIMENSION(:):: a
  REAL,INTENT(OUT):: xmin,xmax, ymin,ymax, zmin,zmax
  
  INTEGER:: k
  INTRINSIC:: HUGE, MAX,MIN, MAXVAL,MINVAL
!----------------------------------------------------------------------------
  xmin=HUGE(xmin)
  ymin=xmin
  zmin=xmin
  xmax=-xmin
  ymax=xmax
  zmax=xmax
  
  DO k=1,SIZE(a)
    xmin=MIN(xmin, MINVAL(a(k)%x) )
    xmax=MAX(xmax, MAXVAL(a(k)%x) )
    ymin=MIN(ymin, MINVAL(a(k)%y) )
    ymax=MAX(ymax, MAXVAL(a(k)%y) )
    zmin=MIN(zmin, MINVAL(a(k)%z) )
    zmax=MAX(zmax, MAXVAL(a(k)%z) )
  END DO
  
  RETURN
END Subroutine GetXYZlimits   ! ---------------------------------------------  

!+
SUBROUTINE ReadOneNetwork(efu, a)
! ---------------------------------------------------------------------------
! PURPOSE - Read one network in LaWgs format. Allocate the memory for 
!  variables x,y,z and load them with the grid point coordinates.
  INTEGER,INTENT(IN):: efu
  TYPE(NETWORK),INTENT(OUT):: a
  
  CHARACTER(LEN=132):: buffer
  INTEGER:: errCode
  INTEGER:: i,j,k,n
  INTEGER:: nobj
  CHARACTER(LEN=132):: title
  INTEGER:: rows,cols
  INTEGER::symLocal,symGlobal
  REAL,DIMENSION(3):: rotate,translate,scale
  
  REAL,ALLOCATABLE,DIMENSION(:):: inpX,inpY,inpZ 
  INTRINSIC:: TRIM
!----------------------------------------------------------------------------
  READ(efu, '(A)') title
  CALL StripFortranQuotes(title)
  WRITE(*,*) "Reading network "//Trim(title)
  
  rotate=0.0
  translate=0.0
  scale=1.0
  READ(efu,'(A)') buffer
  READ(buffer,*,IOSTAT=errCode) &
    nobj, cols, rows, symLocal, rotate,translate,scale, symGlobal
  IF (errCode /= 0) THEN
    WRITE(DBG,*) 'Some parameters on the col,row record missing'  
    WRITE(DBG,*) 'buffer:'//TRIM(buffer)
  END IF  
  n=rows*cols
  ALLOCATE(inpX(n),inpY(n),inpZ(n))
  Read(efu,*) (inpX(k), inpY(k), inpZ(k), k=1,n)

  ALLOCATE(a%x(rows,cols), a%y(rows,cols), a%z(rows,cols) )

  a%name = TRIM(title)
  a%id = nobj
  a%cols = cols
  a%rows = rows
  a%symLocal=symLocal
  a%symGlobal=symGlobal
  a%rotate=rotate
  a%translate=translate
  a%scale=scale
  a%localImage=0
  a%globalImage=0

  k=0
  DO j=1,cols
    DO i=1,rows
      k=k+1
      a%x(i,j)=inpX(k)
      a%y(i,j)=inpY(k)
      a%z(i,j)=inpZ(k)
    END DO
  END DO    
      
  DEALLOCATE(inpX, inpY, inpZ)
  
  RETURN
END Subroutine ReadOneNetwork   ! -------------------------------------------

!+
SUBROUTINE ScanNetworks(efu, nets,ngrid,npanel,maxRows,maxCols,maxGrid)
! ---------------------------------------------------------------------------
! PURPOSE - Read the Wgs file that is already open as unit efu to determine
!  the number of elements defined by the file. Read networks until 
!  end-of-file or a I/O error is encountered or a net with zero rows or cols.
!  You will probably get warning messages that errCode, nobj, x,y,z are
!  never used. True, but we have to read past them in the file. 
  INTEGER,INTENT(IN):: efu
  INTEGER,INTENT(OUT),OPTIONAL:: nets,ngrid,npanel
  INTEGER,INTENT(OUT),OPTIONAL:: maxRows,maxCols,maxGrid    ! maybe optional???

  INTEGER:: k 
  INTEGER:: nobj 
  CHARACTER(LEN=132):: title, buffer
  INTEGER:: rows,cols
  INTEGER:: errCode
  REAL:: x,y,z
  INTRINSIC:: MAX
!----------------------------------------------------------------------------
  IF (Present(nets))    nets=0
  IF (Present(ngrid))   ngrid=0
  IF (Present(npanel))  npanel=0
  IF (Present(maxRows)) maxRows=0
  IF (Present(maxCols)) maxCols=0
  IF (Present(maxGrid)) maxGrid=0

  REWIND(efu)
  READ(efu, '(A)', IOSTAT=errCode) title     ! general title
  IF (errCode /= 0) WRITE(*,*) "Unable to read general title"

  DO
    READ(efu, '(A)', IOSTAT=errCode) title   ! network title
    IF (errCode /= 0) EXIT
    CALL StripFortranQuotes(title)
    WRITE(*,*) "Scanning network "//Trim(title)
    READ(efu, '(A)', IOSTAT=errCode) buffer
    IF (errCode /= 0) EXIT
    READ(buffer, *, IOSTAT=errCode) nobj, cols, rows
    IF (errCode /= 0) EXIT
    IF (rows <= 0) EXIT 
    IF (cols <= 0) EXIT
    READ(efu,*,IOSTAT=errCode) (x,y,z, k=1,rows*cols)
    nobj=nobj+1   ! this and the following 3 statements do absolutely nothing
    x=x+1         ! but they do keep my compiler from constantly nagging me
    y=y+1         ! about the fact that nobj,x,y,z are input but never used
    z=z+1

    IF (Present(nets))    nets=nets+1
    IF (Present(ngrid))   ngrid=ngrid+rows*cols
    IF (Present(npanel))  npanel=npanel+(rows-1)*(cols-1) 
    IF (Present(maxRows)) maxRows=MAX(rows,maxRows)
    IF (Present(maxCols)) maxCols=MAX(cols,maxCols)
    IF (Present(maxGrid)) maxGrid=MAX(rows*cols,maxGrid)
  END DO

  REWIND(efu)
  RETURN
END Subroutine ScanNetworks   ! ---------------------------------------------

!+
SUBROUTINE StripFortranQuotes(t)
! ---------------------------------------------------------------------------
! PURPOSE - If the first non-blank character and the last non-blank
!   character in a string are both APOSTROPHE, then change them to
!   blanks and adjust the resulting string to the left.

  CHARACTER(LEN=*),INTENT(IN OUT):: t

  CHARACTER(LEN=1),PARAMETER:: APOSTROPHE = "'"
  CHARACTER(LEN=LEN(t)):: s   ! I prefer this, but it fails on Absoft F90
!  CHARACTER(LEN=256):: s   ! use this if above is not OK
  INTEGER:: k
  INTRINSIC:: ADJUSTL, LEN_TRIM, TRIM
!----------------------------------------------------------------------------
  s=TRIM(ADJUSTL(t))   ! get rid of leading and trailing blanks
  k=LEN_TRIM(s)
  IF (k==0) RETURN   ! just forget it
  IF( (s(k:k) == APOSTROPHE) .AND. (s(1:1)==APOSTROPHE) ) THEN
    s(k:k)=" "
    s(1:1)=" "
  END IF
  t = TRIM(ADJUSTL(s))
  RETURN
END Subroutine StripFortranQuotes   ! ---------------------------------------

!+
SUBROUTINE TransformNetworks(a) 
!   -------------------------------------------------------------------------
! PURPOSE - Transform the gridpoints in a n array of networks according to
!  the transformation equations. From the LaWgs document, NASA TM-85767, 
!  the order of conversion is rotation, translation, scale.

  TYPE(Network),INTENT(IN OUT),DIMENSION(:):: a
!  REAL,INTENT(IN),DIMENSION(:):: r   ! rotation angles (deg)
!  REAL,INTENT(IN),DIMENSION(:):: t   ! translation distances
!  REAL,INTENT(IN),DIMENSION(:):: s   ! scale factors

  REAL,DIMENSION(3):: before,after
  INTEGER:: i,j,k
  REAL,DIMENSION(3,3):: rot
  INTRINSIC:: MAXVAL, MATMUL
!----------------------------------------------------------------------------
  DO k=1,SIZE(a)
    IF (MAXVAL(ABS(a(k)%rotate))==0.0    .AND.   &
        MAXVAL(ABS(a(k)%translate))==0.0 .AND.   &
        MAXVAL(ABS(a(k)%scale-1.0))==0.0 ) CYCLE
        
    rot=BuildRotationMatrix(a(k)%rotate)

    DO j=1,a(k)%cols
      DO i=1,a(k)%rows
        before(1)=a(k)%x(i,j)
        before(2)=a(k)%y(i,j)
        before(3)=a(k)%z(i,j)
        after=MATMUL(rot,before)
        after=after + a(k)%translate
        after=after * a(k)%scale
        a(k)%x(i,j)=after(1)
        a(k)%y(i,j)=after(2)
        a(k)%z(i,j)=after(3)
      END DO  
    END DO
    
  END DO  
  RETURN
END Subroutine TransformNetworks   ! ----------------------------------------

END Module NetworkProcedures   ! ============================================

!+
MODULE ThreeViewProcedures
!-------------------------------------------------------------------------------
! PURPOSE -
USE NetworkProcedures
IMPLICIT NONE
  CHARACTER(LEN=15):: dateTimeStr
  CHARACTER(LEN=*),PARAMETER:: VERSION = "1.81 (16 Dec 2008)"
  CHARACTER(LEN=132):: wgsFileName
  INTEGER,PARAMETER:: GNU=1, FIG=2

  PUBLIC:: BigWindow
  PUBLIC:: CopyData
  PUBLIC:: CreateGnuFiles
  PUBLIC:: CreateMetafile
  PUBLIC:: GetDateTimeStr
  PRIVATE:: PlotOneNetwork
  PUBLIC:: Welcome

CONTAINS

!+
SUBROUTINE CopyData(efu1,efu2)
! ------------------------------------------------------------------------------
! PURPOSE - Copy a502 input data from unit efu1 to unit efu2.
!   Change CR or LF characters to blanks.
!   Omit A502 style comment records.
IMPLICIT NONE
  INTEGER,INTENT(IN):: efu1   ! external file unit number of the input file
  INTEGER,INTENT(IN):: efu2   ! external file unit number of the output file

  CHARACTER(LEN=132):: buffer
  INTEGER:: errCode
  INTEGER:: idata
  INTEGER:: k
  INTRINSIC:: IACHAR, TRIM
!-------------------------------------------------------------------------------
  idata=0

  DO
    READ(efu1, '(A)', IOSTAT=errCode) buffer
    IF (errCode < 0) EXIT
    
    idata=idata+1
    DO k=1,LEN_TRIM(buffer)
      IF (IACHAR(buffer(k:k))==10) buffer(k:k)=' '
      IF (IACHAR(buffer(k:k))==13) buffer(k:k)=' '
    END DO
    WRITE(efu2, '(A)' ) TRIM(buffer)
  END DO

  WRITE(*, '(I5,A,I9,A)' ) idata, ' total input records '
  RETURN
END Subroutine CopyData   !-------------------------------------------------

!!+
!SUBROUTINE AnnotateFrame(efu,view)
!! ---------------------------------------------------------------------------
!! PURPOSE - 
!  INTEGER,INTENT(IN):: efu
!  INTEGER,INTENT(IN):: view
!  CHARACTER(LEN=9),PARAMETER,DIMENSION(3):: viewTitle = &
!    (/ 'Plan view','Side view','Rear view' /)
!!----------------------------------------------------------------------------
!  WRITE(efu,*) "FONT 6 20"
!  WRITE(efu,*) "SCMOVE 0.5 0.71"
!  WRITE(efu,*) "CTEXT "//viewTitle(view)//" of "//Trim(wgsFileName)
!  WRITE(efu,*) "FONT 5 8"
!  WRITE(efu,*) "SCMOVE 0.99 0.01"
!  WRITE(efu,*) "RTEXT "//dateTimeStr
!  WRITE(efu,*) "ENDFRAME  "//viewTitle(view)//" ----------------------------"
!  RETURN
!END Subroutine AnnotateFrame   ! --------------------------------------------

!+                                                                      
SUBROUTINE BigWindow(vminx, vmaxx, vminy, vmaxy,                  &
  xmin, xmax, ymin, ymax, xmin2, xmax2, ymin2, ymax2) 
! ---------------------------------------------------------------------------
! PURPOSE - Define a window that contains the data in the rectangle
!   [xmin,xmax]X[ymin,ymax] and has the same aspect ratio as the
!   viewport [VMINX,VMAXX]x[VMINY,VMAXY]
!                                                                       
! NOTES-Returns a zero window size if XMIN=XMAX and YMIN=YMAX       
                                                                      
  REAL,INTENT(IN)::  vminx ! left boundary of viewport
  REAL,INTENT(IN)::  vmaxx ! right boundary of viewport
  REAL,INTENT(IN)::  vminy ! bottom boundary of viewport
  REAL,INTENT(IN)::  vmaxy ! top boundary of viewport
  REAL,INTENT(IN)::   xmin ! left boundary of data to be drawn
  REAL,INTENT(IN)::   xmax ! right boundary of data to be drawn
  REAL,INTENT(IN)::   ymin ! bottom boundary of data to be drawn
  REAL,INTENT(IN)::   ymax ! top boundary of data to be drawn
  REAL,INTENT(OUT):: xmin2 ! left boundary of resulting window
  REAL,INTENT(OUT):: xmax2 ! right boundary of resulting window
  REAL,INTENT(OUT):: ymin2 ! bottom boundary of resulting window
  REAL,INTENT(OUT):: ymax2 ! top boundary of resulting window

  REAL:: dx ! window size in x-direction if y is dominant
  REAL:: dy ! window size in y-direction if x is dominant
!----------------------------------------------------------------------------
  IF ((xmax-xmin)*(vmaxy-vminy) > (ymax-ymin)*(vmaxx-vminx)) THEN 
    xmin2=xmin 
    xmax2=xmax 
    dy=(xmax-xmin)*(VMAXY-VMINY)/(VMAXX-VMINX) 
    ymin2=0.5*(ymax+ymin-dy) 
    ymax2=0.5*(ymax+ymin+dy) 
  ELSE 
    dx=(ymax-ymin)*(VMAXX-VMINX)/(VMAXY-VMINY) 
    xmin2=0.5*(xmax+xmin-dx) 
    xmax2=0.5*(xmax+xmin+dx) 
    ymin2=ymin 
    ymax2=ymax 
  END IF
  RETURN
END Subroutine BigWindow  ! -------------------------------------------------

!+
SUBROUTINE CreateGnuFiles()
! ---------------------------------------------------------------------------
! PURPOSE - Read the metafile  describing the three views and create the
!  three gnuplot files showing the plan, side, and rear views.

  CHARACTER(LEN=80):: buffer
  INTEGER:: errCode

  INTEGER:: k,n
  REAL:: x,y
!----------------------------------------------------------------------------
  OPEN(UNIT=FIG, FILE='3view.fig', STATUS='OLD', &
    IOSTAT=errCode, ACTION='READ', POSITION='REWIND')

  OPEN(UNIT=GNU, FILE='plan.gnu', STATUS='REPLACE', &
    IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')
  DO
    READ(FIG,'(A)',IOSTAT=errCode) buffer
    IF (errCode < 0) EXIT
    IF (buffer(1:4)=="DATA") THEN
      READ(buffer(6:10), '(I5)') n
      DO k=1,n
        READ(FIG,*) x,y
        WRITE(GNU,'(2ES15.4)') x,y
      END DO
      WRITE(GNU,'(A)')
    END IF  
    IF (buffer(1:8)=="ENDFRAME") EXIT
  END DO
  CLOSE(UNIT=GNU)
  
  OPEN(UNIT=GNU, FILE='side.gnu', STATUS='REPLACE', &
    IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')
  DO
    READ(FIG,'(A)',IOSTAT=errCode) buffer
    IF (errCode < 0) EXIT
    IF (buffer(1:4)=="DATA") THEN
      READ(buffer(6:10), '(I5)') n
      DO k=1,n
        READ(FIG,*) x,y
        WRITE(GNU,'(2ES15.4)') x,y
      END DO
      WRITE(GNU,'(A)')
    END IF  
    IF (buffer(1:8)=="ENDFRAME") EXIT
  END DO
  CLOSE(UNIT=GNU)
  
  OPEN(UNIT=GNU, FILE='rear.gnu', STATUS='REPLACE', &
    IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')
  DO
    READ(FIG,'(A)',IOSTAT=errCode) buffer
    IF (errCode < 0) EXIT
    IF (buffer(1:4)=="DATA") THEN
      READ(buffer(6:10), '(I5)') n
      DO k=1,n
        READ(FIG,*) x,y
        WRITE(GNU,'(2ES15.4)') x,y
      END DO
      WRITE(GNU,'(A)')
    END IF  
    IF (buffer(1:8)=="ENDFRAME") EXIT
  END DO
  CLOSE(UNIT=GNU)
  CLOSE(UNIT=FIG)

  RETURN
END Subroutine CreateGnuFiles   ! ----------------------------------------

!+
SUBROUTINE CreateMetaFile(a)
! ---------------------------------------------------------------------------
  TYPE(Network),INTENT(IN),DIMENSION(:):: a

  CHARACTER(LEN=*),PARAMETER:: FMTA = "(A)"
  CHARACTER(LEN=*),PARAMETER:: FMTWIN = "('WINDOW ', 4ES15.4)"
  REAL,PARAMETER:: VMINX=0.01, VMAXX=0.99, VMINY=0.01, VMAXY=0.70

  INTEGER:: errCode
  INTEGER:: k
  REAL:: xmin,xmax, ymin,ymax, zmin,zmax
  REAL:: xmin2,xmax2, ymin2,ymax2

!----------------------------------------------------------------------------
  OPEN(UNIT=FIG, FILE='3view.fig', STATUS='REPLACE', &
    IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')
  IF (errCode == 0) THEN
    WRITE(FIG,*) "--- Created by 3view, version "//VERSION//" on "//dateTimeStr
    WRITE(FIG,*) "--- Input file: "//Trim(wgsFileName)
    WRITE(FIG,*) "LANDSCAPE"
    WRITE(FIG,*) "VIEWPORT  0.01 0.99  0.01 0.70"
  ELSE
    WRITE(*,*) "Unable to open graphics metafile"
    STOP
  END IF

  CALL GetXYZlimits(a, xmin,xmax, ymin,ymax, zmin,zmax)

!... Plan View (X-Y)
  CALL BigWindow(VMINX,VMAXX,VMINY,VMAXY, &
    xmin,xmax,ymin,ymax, xmin2,xmax2,ymin2,ymax2)
  WRITE(FIG,FMTWIN) xmin2,xmax2,ymin2,ymax2

  DO k=1,SIZE(a)
    CALL PlotOneNetwork(FIG,a(k),1)
  END DO
  WRITE(FIG,'(A)') "FONT 6 20"
  WRITE(FIG,'(A)') "SCMOVE 0.5 0.71"
  WRITE(FIG,'(A)') "CTEXT Plan View of "//Trim(wgsFileName)
  WRITE(FIG,'(A)') "ENDFRAME  (PLAN VIEW) ------------------------------------"

!... Side View (X-Y)
  CALL BigWindow(VMINX,VMAXX,VMINY,VMAXY, &
    xmin,xmax,zmin,zmax, xmin2,xmax2,ymin2,ymax2)
  WRITE(FIG,FMTWIN) xmin2,xmax2,ymin2,ymax2

  DO k=1,SIZE(a)
    CALL PlotOneNetwork(FIG,a(k),2)
  END DO
  WRITE(FIG,'(A)') "FONT 6 20"
  WRITE(FIG,'(A)') "SCMOVE 0.5 0.71"
  WRITE(FIG,'(A)') "CTEXT Side View of "//Trim(wgsFileName)
  WRITE(FIG,'(A)') "ENDFRAME  (SIDE VIEW) ------------------------------------"


!... Rear View (X-Y)
  CALL BigWindow(VMINX,VMAXX,VMINY,VMAXY, &
    ymin,ymax,zmin,zmax, xmin2,xmax2,ymin2,ymax2)
  WRITE(FIG,FMTWIN) xmin2,xmax2,ymin2,ymax2

  DO k=1,SIZE(a)
    CALL PlotOneNetwork(FIG, a(k), 3)
  END DO
  WRITE(FIG,FMTA) "FONT 6 20"
  WRITE(FIG,FMTA) "SCMOVE 0.5 0.71"
  WRITE(FIG,FMTA) "CTEXT Rear View of "//Trim(wgsFileName)
  WRITE(FIG,FMTA) "FONT 5 8"
  WRITE(FIG,FMTA) "SCMOVE 0.99 0.01"
  WRITE(FIG,FMTA) "RTEXT "//dateTimeStr
  WRITE(FIG,FMTA) "ENDFRAME  (REAR VIEW) ------------------------------------"

  CLOSE(UNIT=FIG)   ! FIG has ACTION='WRITE'; reopen it again with 'READ'
                    ! We could have opened it 'READWRITE' and REWIND at this
                    ! point, but I prefer this method for added safety. 
  RETURN
END Subroutine CreateMetaFile   ! -------------------------------------------


!+
FUNCTION GetDateTimeStr() RESULT(s)
! ---------------------------------------------------------------------------
! PURPOSE - Return a string with the current date and time
  CHARACTER(LEN=*),PARAMETER:: MONTH="JanFebMarAprMayJunJulAugSepOctNovDec"
  CHARACTER(LEN=*),PARAMETER:: FMT = "(I2.2,A1,I2.2,I3,A3,I4)"
  CHARACTER(LEN=15):: s
  INTEGER,DIMENSION(8):: v
!----------------------------------------------------------------------------
  CALL DATE_AND_TIME(VALUES=v)

  WRITE(s,FMT) v(5), ':', v(6), v(3), MONTH(3*v(2)-2:3*v(2)), v(1)
  RETURN
END FUNCTION GetDateTimeStr   ! ---------------------------------------------

!+
SUBROUTINE PlotOneNetwork(efu,a,view)
! ---------------------------------------------------------------------------
! PURPOSE - Plot the rows and columns of the grid of points

  INTEGER,INTENT(IN):: efu
  TYPE(Network),INTENT(IN):: a

  INTEGER,INTENT(IN):: view
  
!  REAL,INTENT(IN),DIMENSION(:,:):: grid
  INTEGER:: rows,cols
  CHARACTER(LEN=*),PARAMETER:: FMTPTS= "(2F15.6)"
  CHARACTER(LEN=*),PARAMETER:: FMTROW= &
    "('DATA ',I5,' for row ',I3)"
  CHARACTER(LEN=*),PARAMETER:: FMTCOL= &
    "('DATA ',I5,' for col ',I3)"

  INTEGER:: i,j
!----------------------------------------------------------------------------
  rows=a%rows
  cols=a%cols

!... Plot the columns
  DO j=1,cols
    IF (j==1 .OR. j==cols) THEN
      WRITE(efu,*) "LINEWIDTH 0.5"
    ELSE
      WRITE(efu,*) "LINEWIDTH 0.2"
    END IF

    WRITE(efu,FMTCOL) rows, j
    SELECT CASE(view)
      CASE(1)
        WRITE(efu,FMTPTS) (a%x(i,j),a%y(i,j),i=1,rows)
      CASE(2)
        WRITE(efu,FMTPTS) (a%x(i,j),a%z(i,j),i=1,rows)
      CASE(3)
        WRITE(efu,FMTPTS) (a%y(i,j),a%z(i,j),i=1,rows)
    END SELECT
  END DO

!... Plot the rows
  DO i=1,rows
    IF (i==1 .OR. i==rows) THEN
      WRITE(efu,*) "LINEWIDTH 0.5"
    ELSE
      WRITE(efu,*) "LINEWIDTH 0.2"
    END IF

    WRITE(efu,FMTROW) cols, i
    SELECT CASE(view)
      CASE(1)
        WRITE(efu,FMTPTS) (a%x(i,j),a%y(i,j),j=1,cols)
      CASE(2)
        WRITE(efu,FMTPTS) (a%x(i,j),a%z(i,j),j=1,cols)
      CASE(3)
        WRITE(efu,FMTPTS) (a%y(i,j),a%z(i,j),j=1,cols)
    END SELECT
  END DO

  RETURN
END Subroutine PlotOneNetwork   ! -------------------------------------------

!+
SUBROUTINE Welcome(efu)
! ---------------------------------------------------------------------------
! PURPOSE -
IMPLICIT NONE
  INTEGER,INTENT(IN):: efu  
  CHARACTER(LEN=*),PARAMETER:: GREETING = "New ThreeView"
  CHARACTER(LEN=*),PARAMETER:: AUTHOR = &
    "Ralph L. Carmichael, Public Domain Aeronautical Software"
  CHARACTER(LEN=*),PARAMETER:: MODIFIER = "none"

  CHARACTER(LEN=132):: fileName
  INTEGER:: errCode
  INTRINSIC:: GET_COMMAND_ARGUMENT,LEN_TRIM
!----------------------------------------------------------------------------


  WRITE(*,*) GREETING
  WRITE(*,*) "Version "//VERSION                   ! VERSION is module variable
  WRITE(*,*) AUTHOR
  WRITE(*,*) "Modified by "//Trim(MODIFIER)

  errCode=99
  CALL GET_COMMAND_ARGUMENT(1,fileName)
  IF (LEN_TRIM(fileName) > 0) THEN
    OPEN(UNIT=efu,FILE=fileName,STATUS='OLD',IOSTAT=errCode,ACTION='READ') 
  END IF

  IF (errCode /= 0) THEN
    DO
      WRITE(*,*) 'Enter the name of the LaWgs file: '
      READ(*,'(A)') fileName   ! fileName is global
      IF (Len_Trim(fileName)==0) STOP
      OPEN(UNIT=efu, FILE=fileName, STATUS='OLD', &
        IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
      IF (errCode==0) EXIT  
      OPEN(UNIT=efu, FILE=Trim(fileName)//'.wgs', STATUS='OLD', &
        IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
      IF (errCode==0) EXIT  
      WRITE(*,*) 'Unable to open this file. Try again.'
    END DO
  END IF
  INQUIRE(UNIT=efu, NAME=wgsFileName)                  ! fileName is global
  WRITE(*,*) "Reading from "//Trim(wgsFileName)

  RETURN
END Subroutine Welcome   ! --------------------------------------------------


END Module ThreeViewProcedures   ! =============================================

!+
PROGRAM ThreeView                                         ! \3view\v14\3view.f90
! ------------------------------------------------------------------------------
! PURPOSE - Plot three orthogonal views of an object defined in LaWGS
!   (Langley Wireframe Geometry Standard) format.
! AUTHOR  - Ralph L. Carmichael, Public Domain Aeronautical Software
! REVISION HISTORY                                                        
!   DATE  VERS PERSON  STATEMENT OF CHANGES                               
! 19Sep95  0.1   RLC   Converted old program ( C language)              
! 25Sep95  0.2   RLC   Got most bugs out                                
! 30Nov95  1.0   RLC   Made Fig have STATUS='DELETE'                    
! 18Dec02  1.5   RLC   Complete rewrite of the imaging code.
! 28Dec02  1.6   RLC   Fixed several errors in the imaging logic
! 14Dec08  1.7   RLC   Moved procedures to module. Removed PostScript logic
USE NetworkProcedures
USE ThreeViewProcedures
IMPLICIT NONE

 
  INTEGER,PARAMETER:: WGS=1, DBG=3, TMP=4      ! units

  INTEGER:: errCode
  CHARACTER(LEN=80):: generalTitle
  INTEGER:: knet
  TYPE(Network),ALLOCATABLE,DIMENSION(:):: networks
  INTEGER:: nnets
 
!-------------------------------------------------------------------------------
  dateTimeStr=GetDateTimeStr()
  OPEN(UNIT=DBG, FILE='3view.dbg', STATUS='REPLACE', &
    IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')
  IF (errCode == 0) THEN
    WRITE(DBG,*) "--- Created by 3view, version "//VERSION
    WRITE(DBG,*) "--- Date: "//dateTimeStr
  ELSE
    WRITE(*,*) "Unable to open log file. Fatal error"
    STOP
  END IF

  CALL Welcome(WGS)

  WRITE(DBG,*) "--- Reading from file: "//Trim(wgsFileName)
  OPEN(UNIT=TMP, FILE='3view.tmp',STATUS='REPLACE',ACTION='WRITE')
  CALL CopyData(WGS,TMP) ! makes a local copy of the WGS file
  CLOSE(UNIT=WGS)
  CLOSE(UNIT=TMP)   ! never write to TMP again
  OPEN(UNIT=TMP,FILE='3view.tmp',STATUS='OLD',ACTION='READ')
  CALL CopyData(TMP,DBG)
  REWIND(UNIT=TMP)
  WRITE(DBG,'(A//)') 'END OF THE LaWgs INPUT FILE'
  
  CALL ScanNetworks(TMP, NETS=nnets)
  ALLOCATE(networks(4*nnets))   ! worst case: each net has a local image and
                                ! each of these has a global image.
  
  READ(TMP,'(A)') generalTitle
  WRITE(*,*) 'The general title of the WGS file is'
  WRITE(*,*) Trim(generalTitle)
  WRITE(DBG,*) "--- General title:"//Trim(generalTitle)
  DO knet=1,nnets
    CALL ReadOneNetwork(TMP,networks(knet))
  END DO  
  WRITE(*,*) nnets, ' networks have been read.'
  CLOSE(UNIT=TMP)
  
  CALL CreateLocalImageNetworks(nnets,networks)
  WRITE(*,*) 'After local imaging, there are', nnets, ' networks.'
  WRITE(DBG,*) 'LIST OF NETWORKS AFTER LOCAL IMAGING'
  CALL DumpNetworkGridPoints(DBG, networks(1:nnets))
 
  CALL TransformNetworks(networks(1:nnets))
    
  CALL CreateGlobalImageNetworks(nnets,networks)
  WRITE(*,*) 'After global imaging, there are', nnets, ' networks.'
  WRITE(DBG,*) 'LIST OF NETWORKS AFTER TRANSFORMATION AND GLOBAL IMAGING'
  CALL DumpNetworkGridPoints(DBG, networks(1:nnets))
 
  CALL CreateMetaFile(networks(1:nnets))
  CALL CreateGnuFiles()
!  CALL CreatePsFile()   ! may get added later

  STOP
END Program ThreeView   ! ======================================================

