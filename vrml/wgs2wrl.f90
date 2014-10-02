!+
MODULE WgsUtilities                                          ! \wgs\wgsstuff.f90
! ---------------------------------------------------------------------------
! PURPOSE - Utility routines for programs using WGS format files          
!    The geometry must be defined in LaWGS                                
!                           (Langley Wire Frame Geometry Standard) format 
! AUTHOR  - Ralph L. Carmichael, Public Domain Aeronautical Software
! REVISION HISTORY                                                        
!   DATE  VERS PERSON  STATEMENT OF CHANGES                               
!  6Sep91  0.1   RLC   Original coding (based on old MODNETS)             
! 26Sep91  0.2   RLC   Revised structure of objects                       
! 26Jun92  0.3   RLC   Objects restructured again                         
! 23Jul92  0.4   RLC   Added the ReadNetworks & ScanNetworks routines     
! 17Aug92  0.5   RLC   Added CONFIGURATION, MAXNETS, LoadNetworkId        
! 17Aug92  0.52  RLC   Added precision to WriteToWGSfile                  
! 16Mar93  0.6   RLC   Dump all old stuff.                                
!  6Jul96  0.7   RLC   Fortran90. Types Network and NetworkNode
!  8Jul96  0.75  RLC   ReadNetworks and ReadNetworksIntoList
! 13Jul96  0.8   RLC   Added PrintSummary
! 22Jun98  0.81  RLC   Revised field widths in PrintSummary
! 26Jun98  0.82  RLC   Added CreateVRMLobject
! 27Jun98  0.83  RLC   Added CountNetworks
! 03Jan05  0.84  RLC   Changed length of variable s in StripFortranQuotes
!----------------------------------------------------------------------------

IMPLICIT NONE


  REAL,PARAMETER,PRIVATE :: PI=3.14159265
  REAL,PARAMETER,PRIVATE:: TWOPI=PI+PI, HALFPI=0.5*PI, RAD=180/PI
  CHARACTER(LEN=*),PARAMETER:: WGSSTUFF_MODULE_VERSION = "0.84 (3 January 2005)"
  CHARACTER(LEN=*),PARAMETER,PRIVATE:: APOSTROPHE = "'"

  TYPE:: Network
    CHARACTER(LEN=132):: name
    INTEGER:: id
    INTEGER:: rows,cols
    REAL,POINTER,DIMENSION(:):: x,y,z
  END TYPE Network

  TYPE:: NetworkNode
    TYPE(Network):: data
    TYPE(NetworkNode),POINTER:: next
  END TYPE NetworkNode

  PUBLIC:: CountNetworksInList
  PUBLIC:: EdgeEndPoints
  PUBLIC:: GetDateTimeStr
  PUBLIC:: IndexOfEdgePoint
  PUBLIC:: NetworkLimits
  PUBLIC:: ScanNetworks
  PUBLIC:: StripFortranQuotes
  PUBLIC:: AppendNetsFromFile
  PUBLIC:: OpenInput
  PUBLIC:: ReadNetworks
  PUBLIC:: ReadNetworksIntoList
  PUBLIC:: ReadNonWakeNetworksIntoList
  PUBLIC:: PrintSummary
  PUBLIC:: PrintData
  PUBLIC:: PrintNetworks502
  PUBLIC:: PrintSize
  PUBLIC:: PrintLimits
  PUBLIC:: GetGlobalLimits
  PUBLIC:: GetNonWakeGlobalLimits
  PRIVATE:: UpCase




CONTAINS


!+
FUNCTION CountNetworksInList(a) RESULT(k)
! ---------------------------------------------------------------------------
! PURPOSE - Return an integer counting the number of networks in a list.
  TYPE(NetworkNode),POINTER:: a
  INTEGER:: k

  TYPE(networkNode),POINTER:: trav
!----------------------------------------------------------------------------
  trav => a   ! copy to avoid destroying original

  k=0
  DO
    IF ( .NOT. Associated(trav) ) EXIT
    k=k+1
    trav => trav%next
  END DO
  RETURN
END Function CountNetworksInList   ! ----------------------------------------

!+
SUBROUTINE EdgeEndPoints(edge, rows,cols, npt1,npt2)
! ---------------------------------------------------------------------------
  INTEGER,INTENT(IN):: edge
  INTEGER,INTENT(IN):: rows,cols
  INTEGER,INTENT(OUT):: npt1,npt2
!----------------------------------------------------------------------------
  SELECT CASE(edge)
    CASE (1)
      npt1=1
      npt2=1+(cols-1)*rows  
    CASE (2)
      npt1=1+(cols-1)*rows
      npt2=rows*cols         
    CASE (3)
      npt1=rows*cols
      npt2=rows              
    CASE (4)
      npt1=rows
      npt2=1                 
  END SELECT 
  RETURN
END Subroutine EdgeEndPoints   ! --------------------------------------------


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
FUNCTION IndexOfEdgePoint(edge,iz, rows, cols, norm) RESULT (i)
! ---------------------------------------------------------------------------
! PURPOSE -Calculate gridpoint index of a mesh point, given its
!        index along a specified side.
! NOTES-Adapted from pilot PANAIR routine MSHIND.                     *)
!
!
!                --->
!           _________________
!           |    edge 1     |
!         ^ |e             e|
!         | |d             d| |
!         | |g             g| |
!         | |e             e| |
!           |               | v
!           |4             2|
!           |    edge 3     |
!           |_______________|
!                <-----
!
!

  INTEGER,INTENT(IN):: edge    ! The old INDEDG function }
  INTEGER,INTENT(IN):: iz 
  INTEGER,INTENT(IN):: rows, cols 
  LOGICAL,INTENT(IN):: norm 
  INTEGER:: i
  INTEGER:: row,col  !  don't confuse with rows, cols 
!----------------------------------------------------------------------------
  SELECT CASE (edge)
    CASE (1)
      row=1
      IF (norm) THEN
        col=iz
      ELSE
        col=cols-iz+1
      END IF
    CASE (2)
      col=cols
      IF (norm) THEN
        row=iz
      ELSE
        row=rows-iz+1
      END IF
    CASE (3)
      row=rows
      IF (norm) THEN
        col=cols-iz+1
      ELSE
        col=iz
      END IF
    CASE (4)
      col=1
      IF (norm) THEN
        row=rows-iz+1
      ELSE
        row=iz
      END IF
  END SELECT 
  i=row+rows*(col-1)
  RETURN
END Function IndexOfEdgePoint   ! -------------------------------------------

!+
SUBROUTINE NetworkLimits(net, xmin,xmax, ymin,ymax, zmin,zmax)
! ---------------------------------------------------------------------------
! PURPOSE - Determine the minimum 3-D box that will contain all the points
!   of one network.

  TYPE(Network),INTENT(IN):: net
  REAL,INTENT(OUT):: xmin,xmax, ymin,ymax, zmin,zmax
!----------------------------------------------------------------------------
  xmin=HUGE(xmin)
  xmax=-xmin
  ymin=xmin
  ymax=xmax
  zmin=xmin
  zmax=xmax

  xmin=MINVAL(net%x)
  xmax=MAXVAL(net%x)
  ymin=MINVAL(net%y)
  ymax=MAXVAL(net%y)
  zmin=MINVAL(net%z)
  zmax=MAXVAL(net%z)

  RETURN
END Subroutine NetworkLimits   ! --------------------------------------------

!+
SUBROUTINE ScanNetworks(efu, nets,ngrid,npanel,maxRows,maxCols,maxGrid)
! ---------------------------------------------------------------------------
  INTEGER,INTENT(IN):: efu
  INTEGER,INTENT(OUT):: nets,ngrid,npanel
  INTEGER,INTENT(OUT):: maxRows,maxCols,maxGrid    ! maybe optional???

  INTEGER:: i 
  INTEGER:: nobj 
  CHARACTER(LEN=132):: title, buffer
  INTEGER:: rows,cols
  INTEGER:: errCode
  REAL:: x,y,z
!----------------------------------------------------------------------------
  nets=0
  ngrid=0
  npanel=0
  maxRows=0
  maxCols=0
  maxGrid=0

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
    nobj=nobj+1   ! USELESS
    READ(efu,*,IOSTAT=errCode) (x,y,z, i=1,rows*cols)
      x=x+1       ! USELESS
      y=y+1       ! USELESS
      z=z+1       ! USELESS
    IF (errCode /= 0) EXIT

    IF (rows <= 0) EXIT 
    IF (cols <= 0) EXIT

    nets=nets+1
    ngrid=ngrid+rows*cols
    npanel=npanel+(rows-1)*(cols-1) 
    IF (rows > maxRows) maxRows=rows
    IF (cols > maxCols) maxcols=cols
    IF (rows*cols > maxGrid) maxGrid=rows*cols

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

!  CHARACTER(LEN=LEN(t)):: s
  CHARACTER(LEN=256):: s
  INTEGER:: k
!----------------------------------------------------------------------------
  s=Trim(AdjustL(t))   ! get rid of leading and trailing blanks
  k=Len_Trim(s)
  IF (k==0) RETURN   ! just forget it
  IF( (s(k:k) == APOSTROPHE) .AND. (s(1:1)==APOSTROPHE) ) THEN
    s(k:k)=" "
    s(1:1)=" "
  END IF
  t = Trim(AdjustL(s))
  RETURN
END Subroutine StripFortranQuotes   ! ---------------------------------------

!+
SUBROUTINE AppendNetsFromFile(efu, efuLocal, s)
! ---------------------------------------------------------------------------
  INTEGER,INTENT(IN):: efu, efuLocal
  CHARACTER(LEN=*),INTENT(IN):: s
  
  INTEGER:: errCode
  CHARACTER(LEN=255):: s1
!----------------------------------------------------------------------------
  OPEN(UNIT=efuLocal, FILE=s, &
     IOSTAT=errCode, STATUS='OLD', ACTION='WRITE', POSITION='REWIND')
  IF (errCode /= 0) THEN
    WRITE(*,*) "Unable to open file "//Trim(s)
    RETURN
  END IF

  READ(efuLocal, '(A)' ) s1   ! discard the title line
  DO
    READ(efuLocal, '(A)', IOSTAT=errCode ) s1
    IF (errCode < 0) EXIT
    WRITE(efu, '(A)' ) Trim(s1)
  END DO

  Close(UNIT=efuLocal)
  WRITE(*,*) "Nets copied from "//s
  RETURN
END Subroutine AppendNetsFromFile   ! =======================================

!+
SUBROUTINE OpenInput(efu)
! ---------------------------------------------------------------------------
  INTEGER,INTENT(IN):: efu

  INTEGER:: code
  CHARACTER(LEN=80):: fileName
!----------------------------------------------------------------------------
  DO
    WRITE(*,*) 'Enter the name of the input file: '
    READ(*, '(A)') fileName
    IF (fileName == ' ') STOP
    OPEN (UNIT=efu, FILE=fileName, STATUS='OLD', IOSTAT=code, ACTION='READ')
    IF (code == 0) EXIT
    OPEN (UNIT=efu, FILE=Trim(fileName)//'.wgs', &
                                   STATUS='OLD', IOSTAT=code, ACTION='READ')
    IF (code == 0) EXIT
    WRITE(*,*) ' Unable to open this file.  Try again.'
  END DO

  RETURN
END SUBROUTINE OpenInput   ! ================================================

!+
SUBROUTINE ReadNetworks(efu, np, maxGrid)
! ---------------------------------------------------------------------------
  INTEGER,INTENT(IN):: efu
  TYPE(NETWORK),INTENT(OUT),DIMENSION(:):: np   !   VAR np : PCONFIGURATION;
  INTEGER,INTENT(IN):: maxGrid
  
  INTEGER:: errCode
  INTEGER:: i,k
  INTEGER:: nobj
  CHARACTER(LEN=132):: title
  INTEGER:: rows,cols
  REAL,ALLOCATABLE,DIMENSION(:):: inpX,inpY,inpZ 
!----------------------------------------------------------------------------
  REWIND (efu)
  READ(efu, '(A)' ) title    ! read the general title 

  ALLOCATE(inpX(maxGrid), inpY(maxGrid), inpZ(maxGrid) )

  DO k=1,SIZE(np)
    READ(efu, '(A)', IOSTAT=errCode) title
    IF (errCode /= 0) RETURN
    CALL StripFortranQuotes(title)
    np(k)%name = title
    WRITE(*,*) "Reading network "//Trim(title)
    READ(efu,*) nobj, cols, rows
    np(k)%id = nobj
    np(k)%cols = cols
    np(k)%rows = rows
    Read(efu,*) (inpX(i), inpY(i), inpZ(i), i=1,rows*cols)
    !!! now load np
    ALLOCATE(np(k)%x(rows*cols))
    ALLOCATE(np(k)%y(rows*cols))
    ALLOCATE(np(k)%z(rows*cols))
    np(k)%x=inpX(1:rows*cols)
    np(k)%y=inpY(1:rows*cols)
    np(k)%z=inpZ(1:rows*cols)
  END DO

  DEALLOCATE(inpX, inpY, inpZ)
  RETURN
END Subroutine ReadNetworks   ! =============================================


!+
SUBROUTINE ReadNetworksIntoList(efu, root)
! ---------------------------------------------------------------------------
  INTEGER,INTENT(IN):: efu   ! external file unit to read from. Assumed open.
  TYPE(NetworkNode),POINTER:: root  ! pointer
  
  INTEGER:: errCode
  INTEGER:: i,n
  INTEGER:: nobj
  CHARACTER(LEN=132):: title
  INTEGER:: rows,cols
  REAL,ALLOCATABLE,DIMENSION(:):: inpX,inpY,inpZ 
  TYPE(NetworkNode),POINTER:: previousNetwork, nextNetwork
!----------------------------------------------------------------------------
  NULLIFY(root)   ! just in case there are no valid networks

  REWIND (efu)
  READ(efu, '(A)' ) title    ! read the general title 
!
!..... Read the first network and start the list
  READ(efu, '(A)', IOSTAT=errCode) title
  IF (errCode /= 0) RETURN
  CALL StripFortranQuotes(title)
  WRITE(*,*) "Reading network "//Trim(title)
  READ(efu,*,IOSTAT=errCode) nobj, cols, rows
  IF (errCode /= 0) THEN
    WRITE(*,*) "Error reading nobj,cols,rows record"
    RETURN
  END IF
!  a%id = nobj
!  a%cols = cols
!  a%rows = rows
  n=rows*cols
  ALLOCATE(inpX(n), inpY(n), inpZ(n), STAT=errCode)
  IF (errCode /= 0) THEN
    WRITE(*,*) "Unable to allocate. Error code=", errCode
    RETURN
  END IF
  Read(efu,*,IOSTAT=errCode) (inpX(i), inpY(i), inpZ(i), i=1,rows*cols)
  IF (errCode /= 0) THEN
    WRITE(*,*) "Error reading data records, i=", i
    DEALLOCATE(inpX, inpY, inpZ)
    RETURN
  END IF
!
  ALLOCATE(root)
  root%data%name = title
  root%data%id = nobj
  root%data%cols = cols
  root%data%rows = rows
  ALLOCATE(root%data%x(n), root%data%y(n), root%data%z(n) )
  root%data%x = inpX
  root%data%y = inpY
  root%data%z = inpZ
  NULLIFY(root%next)
  previousNetwork => root
  DEALLOCATE(inpX, inpY, inpZ)
!
!..... Read the remainder of the networks and add them to the list
  DO 
    READ(efu, '(A)', IOSTAT=errCode) title
    IF (errCode /= 0) RETURN
    CALL StripFortranQuotes(title)
    WRITE(*,*) "Reading network "//Trim(title)
    READ(efu,*,IOSTAT=errCode) nobj, cols, rows
    IF (errCode /= 0) THEN
      WRITE(*,*) "Error reading nobj,cols,rows record"
      RETURN
    END IF
    n=cols*rows
    ALLOCATE(inpX(n), inpY(n), inpZ(n) )
    Read(efu,*,IOSTAT=errCode) (inpX(i), inpY(i), inpZ(i), i=1,n)
    IF (errCode /= 0) THEN
      WRITE(*,*) "Error reading data, i=", i
      RETURN
    END IF

    ALLOCATE(nextNetwork)
    nextNetwork%data%name = title
    nextNetwork%data%id = nobj
    nextNetwork%data%cols = cols
    nextNetwork%data%rows = rows
    ALLOCATE(nextNetwork%data%x(n) )
    ALLOCATE(nextNetwork%data%y(n) )
    ALLOCATE(nextNetwork%data%z(n) )
    nextNetwork%data%x=inpX
    nextNetwork%data%y=inpY
    nextNetwork%data%z=inpZ
    DEALLOCATE(inpX, inpY, inpZ)
    NULLIFY(nextNetwork%next)
    previousNetwork%next => nextNetwork
    previousNetwork => nextNetwork
  END DO

  RETURN
END Subroutine ReadNetworksIntoList   ! -------------------------------------


!+
SUBROUTINE ReadNonWakeNetworksIntoList(efu, root)
! ---------------------------------------------------------------------------
  INTEGER,INTENT(IN):: efu   ! external file unit to read from. Assumed open.
  TYPE(NetworkNode),POINTER:: root  ! pointer
  
  INTEGER:: errCode
  INTEGER:: i,n
  INTEGER:: nobj
  CHARACTER(LEN=132):: title
  INTEGER:: rows,cols
  REAL,ALLOCATABLE,DIMENSION(:):: inpX,inpY,inpZ 
  TYPE(NetworkNode),POINTER:: previousNetwork, nextNetwork
!----------------------------------------------------------------------------
  NULLIFY(root)   ! just in case there are no valid networks

  REWIND (efu)
  READ(efu, '(A)' ) title    ! read the general title 
!
!..... Read the first network and start the list
  READ(efu, '(A)', IOSTAT=errCode) title
  IF (errCode /= 0) RETURN
  CALL StripFortranQuotes(title)
  WRITE(*,*) "Reading network "//Trim(title)
  READ(efu,*,IOSTAT=errCode) nobj, cols, rows
  IF (errCode /= 0) RETURN
!  a%id = nobj
!  a%cols = cols
!  a%rows = rows
  n=rows*cols
  ALLOCATE(inpX(n), inpY(n), inpZ(n) )
  Read(efu,*,IOSTAT=errCode) (inpX(i), inpY(i), inpZ(i), i=1,rows*cols)
  IF (errCode /= 0) THEN
    DEALLOCATE(inpX, inpY, inpZ)
    RETURN
  END IF
!
  ALLOCATE(root)
  root%data%name = title
  root%data%id = nobj
  root%data%cols = cols
  root%data%rows = rows
  ALLOCATE(root%data%x(n), root%data%y(n), root%data%z(n) )
  root%data%x = inpX
  root%data%y = inpY
  root%data%z = inpZ
  NULLIFY(root%next)
  previousNetwork => root
  DEALLOCATE(inpX, inpY, inpZ)
!
!..... Read the remainder of the networks and add them to the list
  DO 
    READ(efu, '(A)', IOSTAT=errCode) title
    IF (errCode /= 0) RETURN
    CALL StripFortranQuotes(title)
    WRITE(*,*) "Reading network "//Trim(title)
    READ(efu,*,IOSTAT=errCode) nobj, cols, rows
    IF (errCode /= 0) RETURN
    n=cols*rows
    ALLOCATE(inpX(n), inpY(n), inpZ(n) )
    Read(efu,*,IOSTAT=errCode) (inpX(i), inpY(i), inpZ(i), i=1,rows*cols)
    IF (errCode /= 0) RETURN

    ALLOCATE(nextNetwork)
    nextNetwork%data%name = title
    nextNetwork%data%id = nobj
    nextNetwork%data%cols = cols
    nextNetwork%data%rows = rows
    ALLOCATE(nextNetwork%data%x(n) )
    ALLOCATE(nextNetwork%data%y(n) )
    ALLOCATE(nextNetwork%data%z(n) )
    nextNetwork%data%x=inpX
    nextNetwork%data%y=inpY
    nextNetwork%data%z=inpZ
    DEALLOCATE(inpX, inpY, inpZ)
    NULLIFY(nextNetwork%next)
    previousNetwork%next => nextNetwork
    previousNetwork => nextNetwork
  END DO

  RETURN
END Subroutine ReadNonWakeNetworksIntoList   ! ------------------------------

!+
SUBROUTINE PrintSummary(a)
! ---------------------------------------------------------------------------
  TYPE(NetworkNode),POINTER:: a

  TYPE(networkNode),POINTER:: trav
  INTEGER:: k
  REAL:: xmin,xmax, ymin,ymax, zmin,zmax
  REAL:: xminGlobal, xmaxGlobal
  REAL:: yminGlobal, ymaxGlobal
  REAL:: zminGlobal, zmaxGlobal
!----------------------------------------------------------------------------

  xminGlobal=HUGE(xminGlobal)
  xmaxGlobal= -xminGlobal
  yminGlobal= xminGlobal
  ymaxGlobal= xmaxGlobal
  zminGlobal= xminGlobal
  zmaxGlobal= xmaxGlobal

  trav => a   ! copy to avoid destroying original

  k=0
  WRITE(*,*) " net  rows cols     ID"
  DO
    IF ( .NOT. Associated(trav) ) EXIT
    k=k+1
    WRITE(*,'(3I5,3X,A)') k, &
       trav%data%rows, trav%data%cols, Trim(trav%data%name)
    trav => trav%next
  END DO
!!!  PauseForUser();

  trav => a   ! start back at the root
  k=0
  WRITE(*,*) "net rows cols     "// &
             "xmin     xmax     ymin     ymax     zmin     zmax"
  DO
    IF ( .NOT.Associated(trav) ) EXIT
    k=k+1
    xmin = MINVAL(trav%data%x)
    xmax = MAXVAL(trav%data%x)
    ymin = MINVAL(trav%data%y)
    ymax = MAXVAL(trav%data%y)
    zmin = MINVAL(trav%data%z)
    zmax = MAXVAL(trav%data%z)
    xminGlobal = MIN(xminGlobal, xmin)
    xmaxGlobal = MAX(xmaxGlobal, xmax)
    yminGlobal = MIN(yminGlobal, ymin)
    ymaxGlobal = MAX(ymaxGlobal, ymax)
    zminGlobal = MIN(zminGlobal, zmin)
    zmaxGlobal = MAX(zmaxGlobal, zmax)
    WRITE(*, '(I3,2I5,6F10.3)' ) k, trav%data%rows, trav%data%cols, &
      xmin,xmax, ymin,ymax, zmin,zmax
    trav => trav%next
  END DO

!  WRITE(*, '(A,I4,A,I4,A,I4,A)' ) "There are", nets, " networks with", &
!     npanel, " panels and", ngrid, " grid points in this file."
  WRITE(*,*) "The data limits are"
  WRITE(*,*) "               min         max"
  WRITE(*, '(A,2F13.3)' ) "   x:", xminGlobal, xmaxGlobal
  WRITE(*, '(A,2F13.3)' ) "   y:", yminGlobal, ymaxGlobal
  WRITE(*, '(A,2F13.3)' ) "   z:", zminGlobal, zmaxGlobal

  RETURN
END Subroutine PrintSummary   ! ---------------------------------------------


!+
SUBROUTINE PrintData(efu,a)
! ---------------------------------------------------------------------------
! PURPOSE - Print the coordinates of each network to a file

  INTEGER,INTENT(IN):: efu   ! external file unit for output
  TYPE(NetworkNode),POINTER:: a

  TYPE(networkNode),POINTER:: trav
  INTEGER:: i,k,n
!----------------------------------------------------------------------------
  trav => a   ! copy to avoid destroying original

  k=0
  DO
    IF ( .NOT. Associated(trav) ) EXIT
    k=k+1
    n=trav%data%rows*trav%data%cols
    WRITE(efu,*) "DATA for "//Trim(trav%data%name)
    WRITE(efu,'(3ES15.5)') (trav%data%x(i), &
                            trav%data%y(i),trav%data%z(i),i=1,n)
    trav => trav%next
  END DO
  RETURN
END Subroutine PrintData   ! ------------------------------------------------

!+
SUBROUTINE PrintNetworks502(efu,a,prec,kt,nlopt1,nropt1,nlopt2,nropt2,mnswitch)
! ---------------------------------------------------------------------------
! PURPOSE - Print the coordinates of each network to a file in the
!   format for input to A502

  INTEGER,INTENT(IN):: efu   ! external file unit for output
  TYPE(NetworkNode),POINTER:: a
  INTEGER,INTENT(IN):: prec
  INTEGER,INTENT(IN),DIMENSION(:):: kt
  INTEGER,INTENT(IN),DIMENSION(:):: nlopt1,nropt1,nlopt2,nropt2
  INTEGER,INTENT(IN),DIMENSION(:):: mnswitch

  INTEGER:: c,r
  CHARACTER(LEN=8):: fmt
  CHARACTER(LEN=*),PARAMETER:: FMT1="(2F10.1,20X,2F10.1)"
  CHARACTER(LEN=*),PARAMETER:: FMT2="( F10.1,30X,2F10.1)"

  TYPE(networkNode),POINTER:: trav
  INTEGER:: j,j1,j2,k,n
  REAL:: nts,ntd
  REAL,ALLOCATABLE,DIMENSION(:):: xyz
!----------------------------------------------------------------------------
  fmt="(6F10.4)"
  WRITE(fmt(7:7), '(I1)') prec
  trav => a   ! copy to avoid destroying original

  k=0
  DO
    IF ( .NOT. Associated(trav) ) EXIT
    k=k+1
    c=trav%data%cols
    r=trav%data%rows
    n=c*r
    WRITE(efu,'(A,I5,3X,A)' ) "$POINTS    for network #", &
      k, Trim(trav%data%name)
    WRITE(efu,*) "1.0"    ! nets per group
    SELECT CASE(kt(k))
      CASE(-18)
        WRITE(efu,FMT1) 18.0, 1.0, REAL(mnswitch(k)), 0.0
      CASE(30)
        IF(nlopt1(k)==0) THEN
          nts=0.0
        ELSE
          nts=1.0
        END IF
        IF(nlopt2(k)==0) THEN
          ntd=0.0
        ELSE
          ntd=12.0
        END IF
        WRITE(efu,'(6F10.1)') REAL(kt(k)), nts,ntd, 0.0, REAL(mnswitch(k)), 0.0
        WRITE(efu,'(6F10.1)') REAL(nlopt1(k)),REAL(nropt1(k)), &
                              REAL(nlopt2(k)),REAL(nropt2(k))
      CASE DEFAULT
        WRITE(efu,FMT2) REAL(kt(k)), REAL(mnswitch(k)), 0.0
    END SELECT
    WRITE(efu,'(2F10.1,40X,A)') REAL(r),REAL(c),trav%data%name(1:10)
    ALLOCATE (xyz(3*n))
    xyz(1:3*n-2:3)=trav%data%x
    xyz(2:3*n-1:3)=trav%data%y
    xyz(3:3*n  :3)=trav%data%z
    DO j=1,c
      j2=3*r*j
      j1=j2-3*r+1
      WRITE(efu,fmt) xyz(j1:j2)
    END DO
    DEALLOCATE(xyz)
!    WRITE(efu,'(3ES15.5)') (trav%data%x(i), &
!                            trav%data%y(i),trav%data%z(i),i=1,n)
    WRITE(*,*) "Printed "//Trim(trav%data%name)//" in a502 format"
    trav => trav%next
  END DO
  RETURN
END Subroutine PrintNetworks502   ! ------------------------------------------------


!+
SUBROUTINE PrintSize(efu,a)
! ---------------------------------------------------------------------------
! PURPOSE - Print the ID,rows,cols and name of each network to a file

  INTEGER,INTENT(IN):: efu   ! external file unit for output
  TYPE(NetworkNode),POINTER:: a

  CHARACTER(LEN=*),PARAMETER:: HEADER= " net ID rows cols  name"
  TYPE(networkNode),POINTER:: trav
  INTEGER:: k
!----------------------------------------------------------------------------
  trav => a   ! copy to avoid destroying original

  k=0
  WRITE(efu,*) HEADER
  DO
    IF ( .NOT. Associated(trav) ) EXIT
    k=k+1
    WRITE(efu,'(4I5,3X,A)') k, trav%data%id, &
       trav%data%rows, trav%data%cols, Trim(trav%data%name)
    trav => trav%next
  END DO
  RETURN
END Subroutine PrintSize   ! ------------------------------------------------


!+
SUBROUTINE PrintLimits(efu,a)
! ---------------------------------------------------------------------------
! PURPOSE - Print the min and max coordinate (in x,y,z) of each network

  INTEGER,INTENT(IN):: efu   ! external file unit for output
  TYPE(NetworkNode),POINTER:: a  

  CHARACTER(LEN=*),PARAMETER:: HEADER = &
    "net rows cols     xmin     xmax     ymin     ymax     zmin     zmax"

  TYPE(networkNode),POINTER:: trav
  INTEGER:: k
  REAL:: xmin,xmax, ymin,ymax, zmin,zmax
!----------------------------------------------------------------------------
  trav => a   ! copy to avoid destroying original
  k=0
  WRITE(efu,*) HEADER
  DO
    IF ( .NOT.Associated(trav) ) EXIT
    k=k+1
    xmin = MINVAL(trav%data%x)
    xmax = MAXVAL(trav%data%x)
    ymin = MINVAL(trav%data%y)
    ymax = MAXVAL(trav%data%y)
    zmin = MINVAL(trav%data%z)
    zmax = MAXVAL(trav%data%z)
    WRITE(efu, '(I3,2I5,6F10.3)' ) k, trav%data%rows, trav%data%cols, &
      xmin,xmax, ymin,ymax, zmin,zmax
    trav => trav%next
  END DO

  RETURN
END Subroutine PrintLimits   ! ---------------------------------------------

!+
SUBROUTINE GetGlobalLimits(a,xmin,xmax,ymin,ymax,zmin,zmax)
! ---------------------------------------------------------------------------
! PURPOSE - Determine the maximum and minimum dimensions of a list of nets.

  TYPE(NetworkNode),POINTER:: a                        ! the root of the list
  REAL,INTENT(OUT):: xmin,xmax, ymin,ymax, zmin,zmax

  TYPE(networkNode),POINTER:: trav
!----------------------------------------------------------------------------
  xmin=HUGE(xmin)
  xmax= -xmin
  ymin= xmin
  ymax= xmax
  zmin= xmin
  zmax= xmax

  trav => a   ! copy to avoid destroying original

  DO
    IF ( .NOT.Associated(trav) ) EXIT
    xmin = MIN(xmin, MINVAL(trav%data%x))
    xmax = MAX(xmax, MAXVAL(trav%data%x))
    ymin = MIN(ymin, MINVAL(trav%data%y))
    ymax = MAX(ymax, MAXVAL(trav%data%y))
    zmin = MIN(zmin, MINVAL(trav%data%z))
    zmax = MAX(zmax, MAXVAL(trav%data%z))
    trav => trav%next
  END DO

  RETURN
END Subroutine GetGlobalLimits   ! ------------------------------------------


!+
SUBROUTINE GetNonWakeGlobalLimits(a,xmin,xmax,ymin,ymax,zmin,zmax)
! ---------------------------------------------------------------------------
! PURPOSE - Determine the maximum and minimum dimensions of a list of nets,
!   except that any networks with the text string 'WAKE' in the name
!   are ignored. Case-insensitive.

  TYPE(NetworkNode),POINTER:: a                        ! the root of the list
  REAL,INTENT(OUT):: xmin,xmax, ymin,ymax, zmin,zmax

  INTEGER:: k
  TYPE(networkNode),POINTER:: trav
  CHARACTER(LEN=132):: tempName
!----------------------------------------------------------------------------
  xmin=HUGE(xmin)
  xmax= -xmin
  ymin= xmin
  ymax= xmax
  zmin= xmin
  zmax= xmax

  trav => a   ! copy to avoid destroying original

  DO
    IF ( .NOT.Associated(trav) ) EXIT
    tempName=trav%data%name
    CALL UpCase(tempName)
    k=INDEX(tempName, "WAKE")
    IF (k == 0) THEN
      xmin = MIN(xmin, MINVAL(trav%data%x))
      xmax = MAX(xmax, MAXVAL(trav%data%x))
      ymin = MIN(ymin, MINVAL(trav%data%y))
      ymax = MAX(ymax, MAXVAL(trav%data%y))
      zmin = MIN(zmin, MINVAL(trav%data%z))
      zmax = MAX(zmax, MAXVAL(trav%data%z))
    END IF
    trav => trav%next
  END DO

  RETURN
END Subroutine GetNonWakeGlobalLimits   ! -----------------------------------

!+
SUBROUTINE UpCase(a)
! ---------------------------------------------------------------------------
! PURPOSE - Convert lower case elements of a character variable to upper
  CHARACTER(LEN=*),INTENT(IN OUT):: a
  INTEGER:: j,k
!----------------------------------------------------------------------------
  DO k=1,LEN(a)
    j=IACHAR(a(k:k))
    IF (j>96 .AND. j<123) a(k:k)=ACHAR(j-32)    ! assumes ASCII character set
  END DO
  RETURN
END Subroutine UpCase   ! ---------------------------------------------------

END Module WgsUtilities   ! =================================================



!+
PROGRAM Wgs2Wrl                                  !   \wgs\wgs2wrl\wgs2wrl.f90
! ---------------------------------------------------------------------------
! PURPOSE -  Create a VRML world of the objects (networks) in a WGS file.  
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software

! REVISION HISTORY
!   DATE   VERS PERSON  STATEMENT OF CHANGES
! 26Jun98  0.1    RLC   Original coding (from wgssumm2 and wgs3view)
!  2Jan00  0.5    RLC   Made compatible with ELF

! NOTES-

USE WgsUtilities

IMPLICIT NONE

  INTEGER,PARAMETER:: DAT=1, WRL=2                         ! I/O unit numbers
  CHARACTER(LEN=*),PARAMETER:: VERSION="0.3 (3Jan2005)"

  CHARACTER (LEN=15) :: dateTimeStr        ! current date and time
  INTEGER:: errCode
  TYPE(NetworkNode),POINTER:: root
!----------------------------------------------------------------------------
  CALL Welcome()
  CALL ReadNetworksIntoList(DAT, root)

  OPEN(UNIT=WRL, FILE='wgs.wrl', IOSTAT=errCode, &
    STATUS='REPLACE', ACTION='WRITE', POSITION='REWIND')
  IF (errCode==0) THEN
    WRITE(WRL,'(A)') "#VRML V1.0 ascii"
  ELSE
    WRITE(*,*) "Unable to open wgs.wrl for output"
    STOP
  END IF

  CALL CreateWorld(WRL,root)
  CALL Goodbye()

  WRITE(*,*) "wgs2wrl has terminated successfully."
  STOP

CONTAINS


!+
SUBROUTINE CreateVRMLobject(efu, a)
! ---------------------------------------------------------------------------
! PURPOSE - Write the description of a network on unit efu using the
!   Coordinate3 and IndexedFaceSet constructs of VRML.
  INTEGER,INTENT(IN):: efu
  TYPE(Network),INTENT(IN):: a

  CHARACTER(LEN=*),PARAMETER:: COMMA=','
  INTEGER,PARAMETER:: PREC=4   ! number of digits of precision

  INTEGER:: i,j,k
  INTEGER,DIMENSION(5):: m
  CHARACTER(LEN=132):: s
  REAL,DIMENSION(3):: xyz
!----------------------------------------------------------------------------
  WRITE(efu,*) "Coordinate3  {"
  WRITE(efu,*) "      point  ["

  k=0
  DO j=1,a%cols
    WRITE(efu,*) "              # column ", j
    DO i=1,a%rows
      k=k+1
      xyz(1)=a%x(k)
      xyz(2)=a%y(k)
      xyz(3)=a%z(k)
      s=RealsToString(xyz,PREC)
      IF (j < a%rows  .OR.  i < a%cols) s=Trim(s)//COMMA
      WRITE(efu,*) Trim(s)
    END DO
  END DO
  WRITE(efu,*) "      ]"
  WRITE(efu,*) "}"

  WRITE(efu,*) "IndexedFaceSet  {"
  WRITE(efu,*) "    coordIndex  ["
  m(5)=-1
  
  DO j=2,a%cols
    WRITE(efu,*) "      # panel column ", j-1
    DO i=2,a%rows
      m(1)=i-2+(j-2)*a%rows
      m(2)=m(1)+1
      m(3)=m(2)+a%rows
      m(4)=m(3)-1
      s=IntegersToStringComma(m)
      s=Trim(s)//COMMA
      WRITE(efu,*) Trim(s)
    END DO

  END DO
  WRITE(efu,*) "      ]"
  WRITE(efu,*) "}"

  RETURN
END Subroutine CreateVRMLobject   ! -----------------------------------------


!+
SUBROUTINE CreateWorld(efu, a)
! ---------------------------------------------------------------------------
  INTEGER,INTENT(IN):: efu
  TYPE(NetworkNode),POINTER:: a   ! root of the list of networks

  TYPE(networkNode),POINTER:: trav

!----------------------------------------------------------------------------

  trav => a   ! copy to avoid destroying original

!  k=0
  DO
    IF ( .NOT. Associated(trav) ) EXIT
    CALL CreateVRMLobject(efu, trav%data)
    trav => trav%next
  END DO
  
  RETURN
END SUBROUTINE CreateWorld  ! -----------------------------------------------

!+
SUBROUTINE GoodBye()
! ---------------------------------------------------------------------------
  CHARACTER (LEN=*), PARAMETER :: FAREWELL= &
     "File wgs.wrl added to your directory."
!----------------------------------------------------------------------------
  WRITE (*,*) FAREWELL
  RETURN
END SUBROUTINE GoodBye  ! ---------------------------------------------------

!+
FUNCTION IntegersToString(a) RESULT(s)
! ---------------------------------------------------------------------------
! PURPOSE - Convert an array of integers to a string with one space
!   between successive integers

  INTEGER,INTENT(IN),DIMENSION(:):: a
  CHARACTER(LEN=255):: s   ! no test is made for possibility of overflowing

  INTEGER:: k
  CHARACTER(LEN=15):: t
!----------------------------------------------------------------------------
  s=" "
  DO k=1,SIZE(a)
    WRITE(t,'(I15)') a(k)
    s=Trim(s)//" "//Trim(AdjustL(t))
  END DO
  RETURN
END Function IntegersToString   ! -------------------------------------------

!+
FUNCTION IntegersToStringComma(a) RESULT(s)
! ---------------------------------------------------------------------------
! PURPOSE - Similar to IntegersToString except that the fields are
!   separated by commas instead of blanks. No comma at end.
  INTEGER,INTENT(IN),DIMENSION(:):: a
  CHARACTER(LEN=255):: s

  INTEGER:: k
  CHARACTER(LEN=15):: t
!----------------------------------------------------------------------------
  WRITE(t,'(I15)') a(1)
  s=AdjustL(t)
  DO k=2,SIZE(a)
    WRITE(t,'(I15)') a(k)
    s=Trim(s)//","//Trim(AdjustL(t))
  END DO
  RETURN
END Function IntegersToStringComma   ! --------------------------------------


!+
FUNCTION RealsToString(a,prec) RESULT(s)
! ---------------------------------------------------------------------------
! PURPOSE - Like IntegersToString, but with reals

  REAL,INTENT(IN),DIMENSION(:):: a
  INTEGER,INTENT(IN):: prec
  CHARACTER(LEN=255):: s

  CHARACTER(LEN=10):: FMT
  INTEGER:: k
  CHARACTER(LEN=15):: t
!----------------------------------------------------------------------------
  WRITE(t,'(I2)') prec
  FMT="(F15."//Trim(AdjustL(t))//")"

  s=" "
  DO k=1,SIZE(a)
    WRITE(t,FMT) a(k)
    s=Trim(s)//" "//Trim(AdjustL(t))
  END DO
  RETURN
END Function RealsToString   ! ----------------------------------------------


!+
SUBROUTINE Welcome()
! ---------------------------------------------------------------------------
  CHARACTER (LEN=*), PARAMETER :: GREETING = &
    "wgssumm - show a summary of a LaWGS file"
  CHARACTER (LEN=*), PARAMETER :: AUTHOR   = &
    "Ralph L. Carmichael, Public Domain Aeronautical Software"
  CHARACTER (LEN=*), PARAMETER :: MODIFIER = " "  ! put your name here
!----------------------------------------------------------------------------
  dateTimeStr=GetDateTimeStr()
!
  WRITE(*,*) GREETING
  WRITE(*,*) "It is now "// dateTimeStr
  WRITE(*,*) AUTHOR
  IF (MODIFIER /= " ") WRITE(*,*) "Modified by "//MODIFIER
  WRITE(*,*) "Version "//VERSION

  CALL OpenInput(Dat)
  RETURN
END SUBROUTINE Welcome  ! ---------------------------------------------------



END Program Wgs2Wrl   ! =====================================================
