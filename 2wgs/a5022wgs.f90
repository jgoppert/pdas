!+
!  PROGRAM A5022WGS
! ------------------------------------------------------------------------------
!     PURPOSE - Convert a A502 input file to a WGS file
!
!     AUTHOR - Ralph L. Carmichael, NASA/AMES
!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 11Feb87  0.1   RLC   Original coding of VHELM2NET
!  6Aug91  0.2   RLC   Rework for A502 and WGS
! 13Oct96  0.3   RLC   Fortran 90
! 23Dec01  0.4   RLC   Look for netname in cols 71-80
! 22Jul08  0.5   RLC   Moved procedures into module; file name from command line
! 02Dec08  0.6   RLC   GET_COMMAND_ARGUMENT to read command line

!+
MODULE A502ToWgsProcedures
! ------------------------------------------------------------------------------
IMPLICIT NONE

  CHARACTER(LEN=*),PARAMETER:: VERSION= "0.6 (2 Dec 2008)"

  REAL,PARAMETER:: PI = 3.14159265   ! used in two subroutines

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

!+
SUBROUTINE UpCase(a)
! ---------------------------------------------------------------------------
! PURPOSE - Make all elements of a character variable upper case

  CHARACTER(LEN=*), INTENT(IN OUT):: a
  INTEGER:: i,j
!----------------------------------------------------------------------------
  DO i=1,LEN(a)
    j=IACHAR(a(i:i))
    IF (j>96 .AND. j<123) a(i:i)=ACHAR(j-32)  ! assumes the ASCII char. set
  END DO
  RETURN
END Subroutine UpCase  ! ----------------------------------------------------

!+
SUBROUTINE PointsToNets(I502, WGS)
! ---------------------------------------------------------------------------
! PURPOSE - Read the points and write them as a WGS network

  INTEGER,INTENT(IN):: I502 ! unit number of the 502 input file
  INTEGER,INTENT(IN):: WGS   ! unit number of the PANAIR output file

  CHARACTER(LEN=*),PARAMETER:: FMT1 = '(6F10.0)'
  INTEGER,PARAMETER:: NDIM=300
  INTEGER:: i,j,k   ! counters
  INTEGER:: nrows, ncols  ! number of rows and columns
  INTEGER:: kn
  REAL:: akn, akt, amn, ann, amnswch
  CHARACTER(LEN=80):: dummy
  CHARACTER(LEN=10):: netName
  CHARACTER(LEN=*),PARAMETER:: QUOTE = "'"
  REAL,DIMENSION(NDIM):: x,y,z
!-----------------------------------------------------------------------
  READ(I502,FMT1) AKN
  KN=INT(AKN)
  READ(I502,'(A)') dummy
  READ(dummy(1:10),FMT1) akt
  READ(dummy(41:50),FMT1) amnswch  ! does not use this right now
!  READ(I502, '(6F10.0)' ) AKT, DUM, DUM, DUM, AMNSWCH

  DO K=1,KN
    READ(I502,'(A)') dummy
    READ(dummy(1:20),FMT1) amn,ann
    netName=dummy(71:80)
!    READ (I502, '(2F10.0,40X,A10)' ) AMN, ANN, NETNAME
    WRITE(WGS,'(A)') QUOTE//netName//QUOTE
!    WRITE(WGS, '(A1,A10,A1)' ) '''', NETNAME, ''''
    nrows=INT(amn)
    ncols=INT(ann)
    WRITE (WGS, '(3I3,A)' ) KN, NCOLS, NROWS, &
                               '   0    0 0 0   0 0 0   1 1 1   1'
    DO J=1,NCOLS
      READ (I502,FMT1) (x(I), y(I), z(I), I=1,NROWS)
      WRITE (WGS, '(6F13.4)' ) (x(I), y(I), z(I), I=1,NROWS)
    END DO
  END DO
  RETURN
END Subroutine PointsToNets   !----------------------------------------------

!+
SUBROUTINE CircToNets(I502, WGS)
! ---------------------------------------------------------------------------
! PURPOSE - Read the 502 points and write them as a PANAIR network

  INTEGER,INTENT(IN):: I502 ! unit number of the 502 input file
  INTEGER,INTENT(IN):: WGS   ! unit number of the PANAIR output file

  INTEGER,PARAMETER:: NDIM= 300
  INTEGER:: I, J, K   ! counters
  INTEGER:: NROWS, NCOLS  ! number of rows and columns
  INTEGER:: KN
  REAL:: AKN, AKT, AMN, ANN, AMNSWCH, ADNSMSH, ANOPT, DUM
  REAL:: S, C
  CHARACTER(LEN=60):: dummy
  CHARACTER(LEN=10):: netName
  REAL,DIMENSION(NDIM):: x,r,deltaz,theta
!----------------------------------------------------------------------------
  READ (I502, '(F10.0)' ) AKN
  KN=INT(AKN)
  READ (I502, '(6F10.0,A10)' ) AKT, DUM, DUM, DUM, AMNSWCH,netName

  DO K=1,KN
    READ (I502, '(6F10.0)' ) ANOPT
    READ (I502, '(6F10.0)' ) AMN
    NROWS=INT(AMN)
    READ (I502, '(6F10.0)' ) (X(I), R(I), I=1,NROWS)
    IF (ANOPT .NE. 0.0) THEN
      READ (I502, '(6F10.0)' ) (DELTAZ(I), I=1,NROWS)
    ELSE
      DO I=1,NROWS
        DELTAZ(I)=0.0
      END DO
    END IF
    READ (I502, '(6F10.0)' ) ANN
    NCOLS=INT(ANN)
    READ (I502, '(6F10.0)' ) (THETA(I), I=1,NCOLS)

    WRITE(WGS, '(A1,A10,A1)' ) '''', NETNAME, ''''
    WRITE(WGS, '(3I3,A)' ) KN, NCOLS, NROWS, &
                               '   0    0 0 0   0 0 0   1 1 1   1'
    DO J=1,NCOLS
      S=SIN(THETA(J)*PI/180.0)
      C=COS(THETA(J)*PI/180.0)
      WRITE (WGS, '(6F13.4)' )   &
                      (X(I), R(I)*C, DELTAZ(I)+R(I)*S, I=1,NROWS)
    END DO
  END DO
  RETURN
END Subroutine CircToNets   !------------------------------------------------

!+
SUBROUTINE QuadToNets(I502, WGS)
! ---------------------------------------------------------------------------
! PURPOSE - Read the 502 points and write them as a PANAIR network
  INTEGER,INTENT(IN):: I502 ! unit number of the 502 input file
  INTEGER,INTENT(IN):: WGS   ! unit number of the PANAIR output file
  
      INTEGER,PARAMETER:: NDIM=300
  INTEGER:: I, J, K   ! counters
  INTEGER:: NROWS, NCOLS  ! number of rows and columns
  INTEGER:: KN
  REAL:: AKN, AKT, AMN, ANN, AMNSWCH, ADNSMSH, ANOPT, DUM
  REAL:: XLE,YLE,ZLE, XTE,YTE,ZTE
  CHARACTER(LEN=60):: dummy
  CHARACTER(LEN=10):: netName
  REAL,DIMENSION(4):: scx,scy,scz
  REAL,DIMENSION(NDIM):: rowpc,colpc
!----------------------------------------------------------------------------
      READ (I502, '(F10.0)' ) AKN
      KN=INT(AKN)
      READ (I502, '(6F10.0,A10)' ) AKT, DUM, DUM, DUM, AMNSWCH, netName

      DO K=1,KN
         READ (I502, '(6F10.0)' ) (SCX(I), SCY(I), SCZ(I), I=1,4)
         READ (I502, '(6F10.0)' ) AMN
         NROWS=INT(AMN)
         READ (I502, '(6F10.0)' ) (ROWPC(I), I=1,NROWS)
         READ (I502, '(6F10.0)' ) ANN
         NCOLS=INT(ANN)
         READ (I502, '(6F10.0)' ) (COLPC(I), I=1,NCOLS)

         WRITE(WGS, '(A1,A10,A1)' ) '''', NETNAME, ''''
         WRITE (WGS, '(3I3,A)' ) KN, NCOLS, NROWS, &
                                    '   0    0 0 0   0 0 0   1 1 1   1'
         DO J=1,NCOLS
            XLE=COLPC(J)*SCX(1)+(1.0-COLPC(J))*SCX(4)
            YLE=COLPC(J)*SCY(1)+(1.0-COLPC(J))*SCY(4)
            ZLE=COLPC(J)*SCZ(1)+(1.0-COLPC(J))*SCZ(4)
            XTE=COLPC(J)*SCX(2)+(1.0-COLPC(J))*SCX(3)
            YTE=COLPC(J)*SCY(2)+(1.0-COLPC(J))*SCY(3)
            ZTE=COLPC(J)*SCZ(2)+(1.0-COLPC(J))*SCZ(3)
            WRITE (WGS, '(6F13.4)' )                        &
                (ROWPC(I)*XLE+(1.0-ROWPC(I))*XTE,          &
                 ROWPC(I)*YLE+(1.0-ROWPC(I))*YTE,          &
                 ROWPC(I)*ZLE+(1.0-ROWPC(I))*ZTE, I=1,NROWS)
         END DO
      END DO
  RETURN
END Subroutine QuadToNets   !------------------------------------------------

!+
SUBROUTINE EllToNets(I502, WGS)
! ---------------------------------------------------------------------------
! PURPOSE - Read the 502 points and write them as a PANAIR network

  INTEGER,INTENT(IN):: I502 ! unit number of the 502 input file
  INTEGER,INTENT(IN):: WGS   ! unit number of the PANAIR output file

      INTEGER,PARAMETER:: NDIM=300
      INTEGER:: I, J, K   ! counters
      INTEGER:: NROWS, NCOLS  ! number of rows and columns
      INTEGER:: KN
      REAL:: AKN, AKT, AMN, ANN, AMNSWCH, ADNSMSH, ANOPT, DUM
      REAL:: S, C
      CHARACTER(LEN=60):: dummy
      CHARACTER(LEN=10):: netName
      REAL,DIMENSION(NDIM):: x,ry,rz,deltaz,theta
!----------------------------------------------------------------------------
      READ (I502, '(F10.0)' ) AKN
      KN=INT(AKN)
      READ (I502, '(6F10.0,A10)' ) AKT, DUM, DUM, DUM, AMNSWCH,netName

      DO K=1,KN
         READ (I502, '(6F10.0)' ) ANOPT
         READ (I502, '(6F10.0)' ) AMN
         NROWS=INT(AMN)
         READ (I502, '(6F10.0)' ) (X(I), RY(I), RZ(I), I=1,NROWS)
         IF (ANOPT .NE. 0.0) THEN
            READ (I502, '(6F10.0)' ) (DELTAZ(I), I=1,NROWS)
         ELSE
            DO I=1,NROWS
               DELTAZ(I)=0.0
            END DO
         END IF
         READ (I502, '(6F10.0)' ) ANN
         NCOLS=INT(ANN)
         READ (I502, '(6F10.0)' ) (THETA(I), I=1,NCOLS)

         WRITE(WGS, '(A1,A10,A1)' ) '''', NETNAME, ''''
         WRITE (WGS, '(3I3,A)' ) KN, NCOLS, NROWS, &
                                    '   0    0 0 0   0 0 0   1 1 1   1'
         DO J=1,NCOLS
            S=SIN(THETA(J)*PI/180.0)
            C=COS(THETA(J)*PI/180.0)
            WRITE (WGS, '(3F12.6)' ) &
                          (X(I), RY(I)*C, DELTAZ(I)+RZ(I)*S, I=1,NROWS)
         END DO
      END DO
  RETURN
END Subroutine EllToNets   !-------------------------------------------------



END Module A502ToWgsProcedures

!+
PROGRAM A5022wgs
! ------------------------------------------------------------------------------
USE A502ToWgsProcedures
IMPLICIT NONE
  INTEGER,PARAMETER:: I502=1 ! unit number for the 502 file
  INTEGER,PARAMETER:: WGS=2  ! unit number for the network file
  INTEGER,PARAMETER:: IALT=3 ! unit number for the scratch file

  CHARACTER(LEN=80):: dummy
  INTEGER:: errCode
!-----------------------------------------------------------------------
  CALL Welcome()   ! opens all files, etc.

  DO
    READ(IALT, '(A)', IOSTAT=errCode) dummy
    IF (errCode < 0) EXIT
    CALL UpCase(dummy(1:4))
    IF (dummy(1:4) .EQ. '$POI') CALL PointsToNets(IALT, WGS)
    IF (dummy(1:4) .EQ. '$CIR') CALL CircToNets(IALT, WGS)
    IF (dummy(1:4) .EQ. '$QUA') CALL QuadToNets(IALT, WGS)
    IF (dummy(1:4) .EQ. '$ELL') CALL EllToNets(IALT, WGS)
  END DO

  CLOSE(UNIT=IALT)
  CLOSE(UNIT=WGS)
  WRITE(*,*) 'File a502.wgs has been created. Normal termination.'
  STOP

CONTAINS

!+
SUBROUTINE Welcome()
! ---------------------------------------------------------------------------
! PURPOSE - Greet user, get input file name, copy all input to scratch file,
!   omitting comments and open the output file.


  CHARACTER(LEN=*),PARAMETER:: GREETING =  &
        "a5022wgs - A program to convert an a502 input file to wgs."
  CHARACTER(LEN=*),PARAMETER:: AUTHOR= &
        " Ralph L. Carmichael, Public Domain Aeronautical Software"
  CHARACTER(LEN=*),PARAMETER:: MODIFIER=' '  ! put your name here 

  CHARACTER(LEN=15):: dateTimeString
  INTEGER:: errCode
  CHARACTER(LEN=132):: fileName
  INTRINSIC:: GET_COMMAND_ARGUMENT, LEN_TRIM
!----------------------------------------------------------------------------
  dateTimeString = GetDateTimeStr()
                                               
  WRITE(*,*) GREETING
  WRITE(*,*) "Version " // VERSION
  WRITE(*,*) AUTHOR
  IF (MODIFIER /= ' ') WRITE(*,*) 'Modified by '//MODIFIER 

  errCode=99
  CALL GET_COMMAND_ARGUMENT(1,fileName)
  IF (LEN_TRIM(fileName) > 0) THEN
    OPEN(UNIT=I502,FILE=fileName,STATUS='OLD',IOSTAT=errCode,ACTION='READ') 
  END IF

  IF (errCode /= 0) THEN
    DO
      WRITE(*,*) "Enter the name of the input file: "
      READ(*,'(A)') fileName
      IF (fileName == ' ') STOP
      OPEN(UNIT=I502, FILE=fileName, IOSTAT=errCode, &
        STATUS='OLD', ACTION='READ', POSITION='REWIND')
      IF (errCode == 0) EXIT
      WRITE(*,*) "Unable to open this file. Try again."
    END DO
  END IF
  INQUIRE (UNIT=I502, NAME=fileName)
  WRITE(*,*) 'Reading from '//fileName

!..... Open the temporary file and copy the input data, removing all
!.     comment records...
  OPEN(UNIT=IALT, IOSTAT=errCode, STATUS='SCRATCH', ACTION='READWRITE')
  IF (errCode /= 0) THEN
    WRITE(*,*) 'Unable to open temporary file'
    STOP
  END IF
  CALL CopyData(I502, IALT)
  CLOSE(UNIT=I502)                 ! never use the original input again
  REWIND (IALT)

!..... Open the wgs file for output and write the general title for LaWgs....
  OPEN(UNIT=WGS, FILE='a502.wgs', IOSTAT=errCode,  &
    STATUS='REPLACE', ACTION='WRITE', POSITION='REWIND')
  IF (errCode /= 0) THEN
   WRITE(*,*) "Unable to open a502.wgs for output"
   STOP
  END IF
  WRITE(WGS,*) &
     "Conversion of file "//Trim(fileName)//" to WGS     "//dateTimeString

  WRITE(*,*) 'Converting 502 records ...'
  RETURN 
END Subroutine Welcome   !---------------------------------------------------

END Program A5022WGS !=======================================================
