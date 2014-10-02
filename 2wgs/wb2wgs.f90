!+                                                                      
! PROGRAM WingBodyToLaWGS 
! ------------------------------------------------------------------------------
! PURPOSE - Make a LaWGS file corresponding to one defined by a
!   namelist input to WingBody

! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software

!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
!  2Aug95  0.1   RLC   Original coding
! 25Sep95  0.2   RLC   Further debugging
! 18Oct95  0.3   RLC   Removed old equivalence stuff in input.
! 13Dec95  1.0   RLC   Production version
! 22Dec96  2.0   RLC   Protect against failure to open output file
! 30Dec01  2.1   RLC   F90 format
! 22Jul08  2.2   RLC   Moved procedures into module. File name on command line
! 02Dec08  2.3   RLC   GET_COMMAND_ARGUMENT to read command line

!-------------------------------------------------------------------------------

!+
MODULE WingBody2LaWgsProcedures
! ------------------------------------------------------------------------------
  CHARACTER(LEN=*),PARAMETER:: VERSION='2.3 (2 December 2008)' 

  PRIVATE:: BodyCamber
  PRIVATE:: BodyRadius
  PUBLIC:: BodyToWgs
  PUBLIC:: CopyData
  PRIVATE:: Ellipsoid
  PRIVATE:: FillArray
  PUBLIC:: GetDateTimeStr
  PRIVATE:: Parabolic
  PRIVATE:: SearsHaack
  PRIVATE:: TAINT
  PRIVATE:: VonKarmanOgive
  PUBLIC:: WingsToWgs


CONTAINS


!+
SUBROUTINE BodyToWgs(efu1,efu2)
! ---------------------------------------------------------------------------
! PURPOSE - Read the BODY namelist and create the body

IMPLICIT NONE 
  INTEGER,INTENT(IN):: efu1,efu2

  REAL,PARAMETER:: PI=3.14159265, HALFPI=0.5*PI

  INTEGER:: errCode
  INTEGER:: i,j 
  INTEGER:: nb=2,nxbody=21,nrows=8,bcode=0
  REAL:: lnose=0.5,lbody=1.0,ltail=0.0
  REAL:: xnose=0.0,rnose=0.0,radius=0.0,rbase=0.0
  REAL:: znose=0.0,zbase=0.0
  INTEGER:: ntab 
  REAL:: xstart=0.0,sref,cbar,refmom 
  REAL:: mach,rnl 
  INTEGER:: oc,opt, bpcode 
  REAL,DIMENSION(20):: bpfs
  REAL:: per

  REAL,DIMENSION(51):: x,r,rprime,z
  REAL,DIMENSION(101):: xbody,rbody,zbody 
  REAL,DIMENSION(12):: sinth,costh 

  NAMELIST/BODY/nb,nxbody,nrows,bcode, ntab,x,r,rprime,z,           &
     & lnose, lbody, ltail, rnose,radius, rbase, znose,zbase, xnose,    &
     & xstart,sref,cbar,refmom,mach,oc,opt, rnl,bpcode,bpfs,per

!-----------------------------------------------------------------------
  READ(UNIT=efu1, NML=BODY, IOSTAT=errCode)
  IF (errCode < 0) RETURN
                                                                        
  IF (nb .LT. 1) THEN 
    WRITE(*,*) 'nb must be 1 or greater. Setting it to 2' 
    nb=2 
  END IF 
  IF (nb .GT. 5) THEN 
    WRITE(*,*) 'nb must be 5 or less. Setting it to 5' 
    nb=5 
  END IF 
  IF (nxbody .LT. 2) THEN 
    WRITE(*,*) 'nxbody must be 2 or greater. Setting it to 21' 
    nxbody=21 
  END IF 
  IF (nxbody .GT. 101) THEN 
    WRITE(*,*) 'nxbody must be 101 or less. Setting it to 101' 
    nxbody=101 
  END IF 

!!!      NB2=NB+NB
!!!      NBODY=NROWS*NB2
  IF (bcode .LT. 0) THEN 
    xnose=x(1) 
    lbody=x(ntab)-x(1) 
  END IF 

  DO i=0,nb+nb 
    sinth(i+1)=SIN(HALFPI*REAL(i)/REAL(nb)) 
    costh(i+1)=COS(HALFPI*REAL(i)/REAL(nb)) 
  END DO 

  IF (bcode .GE. 0) THEN 
    CALL FillArray(xnose,xnose+lbody,xbody(1:nxbody) )

    DO i=1,nxbody                  ! load the rbody, rpbody arrays
      rbody(i)=BodyRadius(bcode,xbody(i),lnose,lbody,ltail,rnose,rbase,radius)
      zbody(i)=BodyCamber(xbody(i), lnose,lbody,ltail, znose,zbase) 
    END DO 
  ELSE 
    CALL FillArray(x(1),x(ntab),xbody(1:nxbody) )
    DO i=1,nxbody 
      CALL TAINT(x,r, xbody(i),rbody(i), ntab,2) 
      CALL TAINT(x,z, xbody(i),zbody(i), ntab,2) 
    END DO 
  END IF 

  IF (xstart .GE. xbody(nxbody)) nrows=0 
                                                    ! load the bpfs arra
  IF (bpcode .EQ. 0) THEN 
    CALL FillArray(xstart,xbody(nxbody),bpfs(1:nrows) ) ! currently not plotted
  END IF 

  WRITE(efu2,*) '''BODY''' 
  WRITE(efu2,'(3I4,A)' ) 1, nb+nb+1, nxbody, ' 0   0 0 0   0 0 0   1 1 1   0'               
  DO j=1,nb+nb+1 
    DO i=1,nxbody 
      WRITE(efu2, '(3F15.5)' ) xbody(i), rbody(i)*sinth(j),zbody(i)+rbody(i)*costh(j)                                  
    END DO 
  END DO 

  RETURN
END Subroutine BodyToWgs   ! ------------------------------------------------

!+
SUBROUTINE WingsToWgs(efu1,efu2)
! ---------------------------------------------------------------------------
! PURPOSE - Read the WING namelist and create the wing
! NOTE - for now, this only makes the mean surface with no thickness,
!  mo camber, no twist.

  IMPLICIT NONE
  INTEGER,INTENT(IN):: efu1,efu2
  INTEGER:: errCode

  REAL:: rootle=0.0,rootte=1.0,rooty=0.0,rootz=0.0
  REAL:: tiple=0.0,tipte=1.0,tipy=1.0,tipz=0.0 
  INTEGER:: rows=6,cols=8, type=4 
  REAL,DIMENSION(51):: f,g,p,shear 
  REAL:: tcroot,tctip 
  INTEGER:: sect 
  INTEGER:: oc 
  REAL:: mach,sref,refmom,cbar 
  INTEGER:: opt 
  REAL:: rnl 
  REAL:: per 

  INTEGER:: i,j 
  REAL:: theta 

  NAMELIST/WING/ROOTLE,ROOTTE,ROOTY,ROOTZ,                          &
     & TIPLE,TIPTE,TIPY,TIPZ,                                           &
     & ROWS,COLS, TYPE, F,G,P,SHEAR,                                    &
     & TCROOT,TCTIP,  SECT,                                             &
     & OC,MACH,SREF,REFMOM,CBAR, OPT, RNL, PER                          

!----------------------------------------------------------------------------
  DO 
    READ(UNIT=efu1, NML=WING, IOSTAT=errCode)
    IF (errCode < 0) EXIT

    IF (type==3 .OR. type==4) THEN 
      CALL FillArray(rootle,rootte,f(1:rows+1)) 
      CALL FillArray(tiple,tipte,g(1:rows+1)) 
    END IF 

    IF (type==2 .OR. type==4) THEN 
      CALL FillArray(rooty,tipy,p(1:cols+1)) 
      CALL FillArray(rootz,tipz,shear(1:cols+1)) 
    END IF 

    WRITE(efu2, '(A)' ) '''WING''' 
    WRITE(efu2,'(3I4,A)' ) 11,                                       &
     &    cols+1, rows+1, ' 0   0 0 0   0 0 0   1 1 1   0'              
    DO j=1,cols+1 
      theta=(p(j)-p(1))/(p(cols+1)-p(1)) 
      DO i=1,rows+1 
        WRITE(efu2, '(3F15.5)' ) theta*g(i) + (1-theta)*f(i), p(j), shear(j)               
      END DO 
    END DO 

  END DO 

  RETURN
END Subroutine WingsToWgs   ! -----------------------------------------------

!+
FUNCTION BodyCamber(x,lnose,lbody,ltail,znose,zbase) RESULT(z)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the z-coordinate of the center of the fuselage
!  cross-section as a function of longitudinal station x

IMPLICIT NONE 
  REAL,INTENT(IN):: x,lnose,lbody,ltail,znose,zbase 
  REAL:: z

  REAL:: xx
!----------------------------------------------------------------------------
  IF (x < lnose) THEN 
    xx=(lnose-x)/lnose                                     ! nose region
    z=znose*xx*xx
  ELSE IF (x .GT. lbody-ltail) THEN
    xx=(x-lbody+ltail)/ltail                          ! afterbody region
    z=zbase*xx*xx
  ELSE
    z=0.0
  END IF

  RETURN
END Function BodyCamber   ! -------------------------------------------------


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
  INTEGER:: idata,icomment
  INTEGER:: k
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
    WRITE(efu2, '(A)' ) Trim(buffer)
  END DO

  WRITE(*, '(I5,A,I9,A)' ) idata, ' total input records '
  RETURN
END Subroutine CopyData   !-------------------------------------------------


!+
SUBROUTINE TAINT(XTAB,FTAB,X,FX,N,K) 
! ---------------------------------------------------------------------------
! PURPOSE - Table interpolation
  IMPLICIT NONE 

  REAL,INTENT(IN),DIMENSION(:):: xtab,ftab
  REAL,INTENT(IN):: x 
  REAL,INTENT(OUT):: fx 
  INTEGER,INTENT(IN):: n 
  INTEGER,INTENT(IN):: k 


  REAL,DIMENSION(10):: c
      INTEGER:: i,j 
      INTEGER:: kp1 
      INTEGER:: l,m 
      REAL,DIMENSION(10):: t
!----------------------------------------------------------------------------
  IF (n<=1) THEN
    fx=ftab(1)
    RETURN
  END IF

  DO i=1,n
    IF (x.GT.xtab(i)) CYCLE  ! sequential search, tsk,tsk!         
    j=i 
    GO TO 18 
  END DO
  j=n


   18 J=J-(K+1)/2 
      IF(J.LE.0)J=1 
   20 M=J+K 
      IF(M.LE.N)GO TO 21 
      J=J-1 
      GO TO 20 
   21 KP1=K+1 

      DO  L=1,K+1        ! Lagrange polynomial of order K
      C(L)=X-XTAB(J) 
      T(L)=FTAB(J) 
      J=J+1
      END DO


      DO J=1,K 
      I=J+1 
   25 T(I)=(C(J)*T(I)-C(I)*T(J))/(C(J)-C(I)) 
      I=I+1 
      IF(I.LE.KP1)GO TO 25 
      END DO

      FX=T(K+1) 
  RETURN 
END Subroutine Taint   ! ----------------------------------------------------

!+
SUBROUTINE FillArray(start,end, array, spacingCode)
! ---------------------------------------------------------------------------
! PURPOSE - fill an array from start to end. The intermediate points are
!    computed according to various spacing rules.

  REAL,INTENT(IN):: start,end
  REAL,INTENT(OUT),DIMENSION(:):: array
  INTEGER,INTENT(IN),OPTIONAL:: spacingCode
                            ! =2 full cosine
                            ! =3 half cosine
                            ! =4 half sine
                            ! anything else (or nothing)= uniform spacing
  INTEGER:: k,n
  REAL,PARAMETER:: PI=3.14159265, HALFPI=0.5*PI
  REAL,ALLOCATABLE,DIMENSION(:):: temp
!----------------------------------------------------------------------------
  n=SIZE(array)
  IF (n <= 0) RETURN
  array(n)=end
  array(1)=start
  IF (n <= 2) RETURN
  ALLOCATE(temp(n-2))
  temp= (/ (REAL(k), k=1,n-2) /) / REAL(n-1)

  IF (Present(spacingCode)) THEN
    SELECT CASE(spacingCode)
      CASE (2)
        temp=0.5*(1.0-COS(PI*temp))   ! full cosine, dense near both ends
      CASE (3)
        temp=1.0-COS(HALFPI*temp)         ! half cosine, dense near start
      CASE (4)
        temp=SIN(HALFPI*temp)                 ! half sine, dense near end
    END SELECT
  END IF

  array(2:n-1)=start + (end-start)*temp

  DEALLOCATE(temp)
  RETURN
END Subroutine FillArray   ! ------------------------------------------------


!+
FUNCTION Parabolic(x) RESULT(z)
! ---------------------------------------------------------------------------
  REAL,INTENT(IN):: x
  REAL:: z
!----------------------------------------------------------------------------
  z=x*(2.0-x)
  RETURN
END Function Parabolic   ! --------------------------------------------------

FUNCTION SearsHaack(x) RESULT(z)
! ---------------------------------------------------------------------------
  REAL,INTENT(IN):: x
  REAL:: z
!----------------------------------------------------------------------------
  z=(x*(2.0-x))**0.75
  RETURN
END Function SearsHaack   ! -------------------------------------------------

!+
FUNCTION VonKarmanOgive(x) RESULT(z)
! ---------------------------------------------------------------------------
  REAL,INTENT(IN):: x
  REAL:: z

  REAL,PARAMETER:: PI=3.14159265
  REAL:: area 
!----------------------------------------------------------------------------
  area=2.0*(ASIN(SQRT(x)) - (1.0-x-x)*SQRT(x*(1.0-x))) 
  z=SQRT(area/PI)
  RETURN
END Function VonKarmanOgive   ! ---------------------------------------------

!+
FUNCTION Ellipsoid(x) RESULT(z)
! ---------------------------------------------------------------------------
  REAL,INTENT(IN):: x
  REAL:: z
!----------------------------------------------------------------------------
  z=SQRT(x*(2.0-x)) 
  RETURN
END Function Ellipsoid   ! --------------------------------------------------

!+                                                                      
FUNCTION BodyRadius(bcode,x,lnose,lbody,ltail,rnose,rtail,radius) RESULT(r)
! ---------------------------------------------------------------------------
! PURPOSE - Calculate radius at a given point in a
!   nose-cylinder-afterbody fuselage.

IMPLICIT NONE 

  INTEGER,INTENT(IN):: bcode     ! =0   parabolic body
                                 ! =1   Sears-Haack-Adams body
                                 ! =2   Von Karman ogive
                                 ! =3   ellipsoidal body
                                 ! =4   conical body

  REAL,INTENT(IN):: x 
  REAL,INTENT(IN):: lnose,lbody,ltail 
  REAL,INTENT(IN):: rnose,rtail,radius 

  REAL:: r,xx 
!-------------------------------------------------------------------------------
  IF (x .LT. lnose) THEN 
    xx=x/lnose                                 ! non-dimensional nose coordinate
    SELECT CASE (bcode)
      CASE(0) 
        r=rnose+(radius-rnose)*Parabolic(xx)                    ! parabolic nose
      CASE(1) 
        r=rnose+(radius-rnose)*SearsHaack(xx)                  ! SearsHaack nose
      CASE(2)                                                       
        r=rnose+(radius-rnose)*VonKarmanOgive(xx)               ! VonKarman nose
      CASE(3) 
        r=rnose+(radius-rnose)*Ellipsoid(xx)                  ! ellipsoidal nose
      CASE DEFAULT 
        r=rnose+(radius-rnose)*xx                                 ! conical nose
    END SELECT 
  ELSE IF (x .GT. lbody-ltail) THEN
    xx=(lbody-x)/ltail                    ! non-dimensional afterbody coordinate
    SELECT CASE (bcode)
      CASE(0) 
        r=rtail+(radius-rtail)*Parabolic(xx)                    ! parabolic tail
      CASE(1) 
        r=rtail+(radius-rtail)*SearsHaack(xx)                  ! SearsHaack tail
      CASE(2) 
        r=rtail+(radius-rtail)*VonKarmanOgive(xx)               ! VonKarman tail
      CASE(3) 
        r=rtail+(radius-rtail)*Ellipsoid(xx)                  ! ellipsoidal tail
      CASE DEFAULT 
        r=rtail+(radius-rtail)*xx                                 ! conical tail
    END SELECT 
  ELSE
    r=radius                                                ! cylindrical region
  END IF 
  RETURN 
END Function BodyRadius   ! -------------------------------------------------

!+
FUNCTION GetDateTimeStr() RESULT(s)
! ---------------------------------------------------------------------------
! PURPOSE - Return a string with the current date and time
  CHARACTER(LEN=*),PARAMETER:: MONTH='JanFebMarAprMayJunJulAugSepOctNovDec'
  CHARACTER(LEN=*),PARAMETER:: FMT = '(I2.2,A1,I2.2,I3,A3,I4)'
  CHARACTER(LEN=15):: s
  INTEGER,DIMENSION(8):: v
!----------------------------------------------------------------------------
  CALL DATE_AND_TIME(VALUES=v)

  WRITE(s,FMT) v(5), ':', v(6), v(3), MONTH(3*v(2)-2:3*v(2)), v(1)
  RETURN
END FUNCTION GetDateTimeStr   ! ---------------------------------------------



END Module WingBody2LaWgsProcedures   ! ==============================================


!+
PROGRAM WingBodyToLaWgs
! ------------------------------------------------------------------------------

USE WingBody2LaWgsProcedures
IMPLICIT NONE

  INTEGER,PARAMETER:: IN=1,WGS=2,TMP=3 


  
!----------------------------------------------------------------------------
  CALL Welcome()

  WRITE (*,*) 'Looking for BODY...' 
  CALL BodyToWgs(TMP,WGS) 
  REWIND(UNIT=TMP)

  WRITE(*,*) 'Looking for WINGs ...' 
  CALL WingsToWgs(TMP,WGS) 
  CLOSE(IN) 
  CLOSE(WGS) 
  WRITE(*,*) 'File wb.wgs has been added to your directory.'
  WRITE(*,*) 'wb2wgs has terminated successfully.'
  STOP

CONTAINS



!+
SUBROUTINE Welcome()
! ------------------------------------------------------------------------------
! PURPOSE - Greet user, get file name, open input and output files

  CHARACTER(LEN=*),PARAMETER:: GREETING = &
    'wb2wgs - convert wingbody input to LaWGS'
  CHARACTER(LEN=*),PARAMETER:: AUTHOR = &
    'Ralph L. Carmichael, Public Domain Aeronautical Software'
  CHARACTER(LEN=*),PARAMETER:: MODIFIER=' '  ! put your name here 

  CHARACTER(LEN=15):: dateTimeStr
  INTEGER:: errCode
  CHARACTER(LEN=132):: fileName
  INTRINSIC:: GET_COMMAND_ARGUMENT,LEN_TRIM
!-------------------------------------------------------------------------------
  dateTimeStr=GetDateTimeStr()
  WRITE(*,*) GREETING
  WRITE(*,*) "Version " // VERSION 
  WRITE(*,*) AUTHOR 
  IF (MODIFIER /= ' ') WRITE(*,*) 'Modified by '//MODIFIER 

  errCode=99
  CALL GET_COMMAND_ARGUMENT(1,fileName)
  IF (Len_Trim(fileName) > 0) THEN
    OPEN(UNIT=IN,FILE=fileName,STATUS='OLD',IOSTAT=errCode,ACTION='READ') 
  END IF

  IF (errCode /= 0) THEN
    DO
      WRITE(*,*) 'Enter the name of the WingBody input file: '
      READ(*, '(A)' ) fileName 
      IF (fileName == ' ') STOP 
      OPEN(UNIT=IN,FILE=fileName,STATUS='OLD',IOSTAT=errCode,ACTION='READ') 
      IF (errCode == 0) EXIT
      OPEN(UNIT=IN,FILE=Trim(fileName)//'.inp', &
        STATUS='OLD',IOSTAT=errCode,ACTION='READ')
      IF (errCode==0) EXIT
      OPEN(UNIT=IN,FILE=Trim(fileName)//'.nml', &
        STATUS='OLD',IOSTAT=errCode,ACTION='READ')
      IF (errCode==0) EXIT
      WRITE (*,*) 'Unable to open this file. Try again' 
    END DO
  END IF

  INQUIRE(UNIT=IN,NAME=fileName)
  WRITE(*,*) "Reading from "//Trim(fileName)
  OPEN(UNIT=TMP,FILE='wb2wgs.tmp',STATUS='REPLACE',IOSTAT=errCode,ACTION='READWRITE') 
  CALL CopyData(IN,TMP)
  CLOSE(UNIT=IN)
  REWIND(UNIT=TMP)

  OPEN(UNIT=WGS,FILE='wb.wgs',STATUS='REPLACE',ACTION='WRITE') 
  WRITE(WGS,*) 'Created by wb2wgs from '//Trim(fileName)//"    "//dateTimeStr
  RETURN
END Subroutine Welcome   ! --------------------------------------------------

END Program WingBodyToLaWgs   ! ==============================================
