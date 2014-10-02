!+
PROGRAM PanAirInput
!   --------------------------------------------------------------------
! PURPOSE - Prepare an input file for PanAir (A502). The goal is to be
! able to prepare a readable file with the boundary conditions, flight
! conditions and the name of the geometry file (in LaWgs format).
! This program this data into the efficient, but obscure A502 input
! format.
!
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
!
! REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
!  7Jun89  0.1   RLC   Copied from Pascal version
! 22Jul89  0.3   RLC   Changes for UNICOS A502
! 10Aug89  0.5   RLC   Added creation of shell script
! 11Aug89  0.6   RLC   New keyword "PATH"
! 25Sep89  0.7   RLC   Fixed error in OUTCODES list
! 15Feb90  0.71  RLC   Fixed error in type 30 output list
! 16Feb90  0.72  RLC   Solved several problems in streamline section
!  1Jun90  0.75  RLC   Revised script to do rcp from ra-iris
!  9Jul90  0.76  RLC   Save/restore from columbia, not AFS
! 12Jul90  0.8   RLC   XYZ points from file; kt=-18 for base wakes
! 13Jul90  0.81  RLC   Separated subroutine SHELL; diff versions vax/iris
!  4Sep90  0.84  RLC   tape13, ignore, average
! 28Sep90  0.85  RLC   bet1 values for kt=4,9, or 14
! 15Oct90  0.86  RLC   STREAM points from file
!  3Jan91  0.87  RLC   CUTS for sectional properties from a file
! 10Jan91  0.88  RLC   Get add'l files when using CUTS
! 18Jan91  0.89  RLC   Put network name on the $POINTS record
! 15Feb91  0.90  RLC   Added MEMORY,MEM keywords and filesize estimate
! 12Mar91  0.91  RLC   Added OFFNETS, DIVERT
! 10Jul91  0.92  RLC   Changed memory size; get exec from columbia
! 31Jul91  0.93  RLC   New names for some A502 files (version H)
!  1Aug91               (another minor fix to names)
! 26Aug91  0.94  RLC   Allow for two planes of symmetry
!  9Sep91  0.95  RLC   Added WAKE; read BETAs on kt=30 networks (if needed)
! 24Sep91  0.96  RLC   wake nets add to count of nets
!  9Oct91  0.97  RLC   added COPYDATA; allows a502 style comments
! 17Oct91  0.971 RLC   added NOF
! 24Oct91  0.98  RLC   fixed problems with kt=30 source or doublet only
! 11Aug92  0.995 RLC   put PEA data after the network data
!  1Sep92  1.0   RLC   fixed the last bug
!  9Nov92  1.1   RLC   ignored and average nets into sectional properties
! 14Jan93  1.2   RLC   added BWAKE
!  8Feb93  1.3   RLC   omit writing number of nets. A502 no longer needs it
! 26Sep93  1.4   RLC   increased linelength in auxiliary file
! 18Jan94  1.41  RLC   minor mod to shellsub,eaglesub
! 27Jan94  1.42  RLC   returning iflggp instead of ifil63
!  4Feb94  1.5   RLC   change way of making scratch space on eagle
! 15Feb94  1.6   RLC   allow multiple CUTS commands
!  1Aug98  0.99  RLC   total conversion to f90. Defined types
! 20Aug98  0.991 RLC   added PrintSummary after reading nets
! 17Dec98  0.992 RLC   last touches using lf95
!  4Jan00  1.0   RLC   adapted to ELF 

IMPLICIT NONE


  TYPE:: PanAirNetwork
    CHARACTER(LEN=132):: name
    INTEGER:: id
    INTEGER:: rows,cols
    INTEGER:: kt
    INTEGER:: nlopt1,nropt1, nlopt2,nropt2
    INTEGER:: mnSwitch
    INTEGER:: ignore
    INTEGER:: average
    REAL,POINTER,DIMENSION(:):: x,y,z
  END TYPE PanAirNetwork


!***********************************************************************
!     G L O B A L   C O N S T A N T S
!***********************************************************************
  INTEGER, PARAMETER:: IWGS=2                  ! unit number of the WGS file
  INTEGER, PARAMETER:: IAUX=1! unit number of the copy of the auxiliary file
  INTEGER, PARAMETER:: I502=3              ! unit number of the a502.in file
  INTEGER, PARAMETER:: IAUXTEMP=8  ! unit number of the actual auxiliary file
  INTEGER, PARAMETER:: DBG=9   ! debug output file

  INTEGER, PARAMETER:: MAXNETS = 150  ! maximum number of networks
  INTEGER, PARAMETER:: MAXSTRM = 100 ! maximum number of streamlines
  INTEGER, PARAMETER:: MAXOFB  = 1000 !
  INTEGER, PARAMETER:: MAXCUTS = 50 !
  INTEGER, PARAMETER:: MAXABUTS = 50
  INTEGER, PARAMETER:: MAXEDGES = 10

  CHARACTER(LEN=*),PARAMETER :: VERSION='1.0 (4Jan2000)'
  CHARACTER(LEN=*),PARAMETER :: &
           FAREWELL=" Normal termination of panin, version "// VERSION

!============================================================================
!  G L O B A L    V A R I A B L E S                                                       !
!============================================================================
  CHARACTER(LEN=80):: a502File
  REAL:: alpc = 0.0   ! compressibility angle of attack
  REAL:: ABSERR=0.01
  REAL:: betc = 0.0   ! compressibility angle of sideslip

  INTEGER:: CASES=1
  REAL:: CBAR=1.0   ! reference chord, used for pitching moment
  INTEGER:: DATACHK = 0
  CHARACTER(LEN=15):: dateTimeStr                     ! current date and time
  REAL:: direction=1.0
  INTEGER:: errCode

  REAL:: HMIN=0.1
  REAL:: HMAX=0.5

  INTEGER:: IDEFAULT=0
  INTEGER:: IVCORR = 0

  INTEGER:: k

  REAL:: mach = 0.25
  REAL:: MAXORDR=6.0
  INTEGER:: MISYMM = 1   !
  INTEGER:: MJSYMM = 0   !

  CHARACTER(LEN=20) :: NAME = 'a502' ! job name; used for forming file names

  INTEGER:: nets=0          ! number of networks in the file     !
  INTEGER:: NETSIGN=0
  
  INTEGER:: NETSAVG=0
  INTEGER:: NETSOFF=0
  INTEGER:: NWAKES=0
  INTEGER:: nbwakes=0


  INTEGER:: NGRIDS=0
  INTEGER:: NCUTS=0

  INTEGER:: NUMPTS=0
  INTEGER:: NXYZ=0   ! off-body and streamline variables ................

  REAL:: pointsPerStreamline=100.0
  INTEGER:: PREC = 4
  INTEGER:: RESTART = 0
  REAL:: SPAN=1.0   ! reference span, used for yawing and rolling
  REAL:: SREF=1.0   ! reference area
  
  CHARACTER(LEN=80):: title1 = "INPUT FILE FOR PanAir (A502-ht2)"
  CHARACTER(LEN=80):: title2 = "Created by PANIN"

  REAL:: TOL = 0.0


  CHARACTER(LEN=80):: wgsFile 
  LOGICAL:: writeNOF=.FALSE. ! if TRUE, puts $NOF flag in input file



  REAL:: XREF=0.0   ! x-coor of the moment reference point
  REAL:: YREF=0.0   ! y-coor of the moment reference point
  REAL:: ZREF=0.0   ! z-coor of the moment reference point
  



!==============================================================================
!     A R R A Y S                                                             !
!==============================================================================
  INTEGER,DIMENSION(17):: outCodes =(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
  REAL,DIMENSION(4):: ALPHA = (/ 0.0, 0.0, 0.0, 0.0 /)
  REAL,DIMENSION(4):: BETA  = (/ 0.0, 0.0, 0.0, 0.0 /)
  TYPE(PanAirNetwork), DIMENSION(MAXNETS):: np
!  INTEGER, DIMENSION(MAXNETS) :: KT
!  INTEGER, DIMENSION(MAXNETS) :: nlopt1,nropt1, nlopt2,nropt2
!  INTEGER, DIMENSION(MAXNETS) :: MNSWITCH
  INTEGER, DIMENSION(MAXNETS) :: NETIGN, NETAVG, OFFNETS
  INTEGER, DIMENSION(MAXNETS) :: WAKE ! number of net to which wake is attached
  INTEGER, DIMENSION(MAXNETS) :: WAKEEDGE ! edge to which wake is attached
  REAL, DIMENSION(MAXNETS) :: WAKEXTE      ! x coor of downstream edge of wake
!!! analysis says bwake and bwakeedge never used
  INTEGER, DIMENSION(MAXNETS) :: bwake ! number of net to which wake is attached
  INTEGER, DIMENSION(MAXNETS) :: bWakeEdge ! edge to which wake is attached
  REAL, DIMENSION(MAXNETS) :: bWakeXte      ! x coor of downstream edge of wake
  REAL, DIMENSION(MAXNETS) :: BETA1   ! never used???
  REAL, DIMENSION(MAXNETS) :: BETA2  ! never used??

  REAL, DIMENSION(MAXOFB) :: XOFF,YOFF,ZOFF, XORIGIN,YORIGIN,ZORIGIN
  REAL, DIMENSION(MAXOFB) :: DXGRID,DYGRID,DZGRID
  INTEGER, DIMENSION(MAXOFB) :: NXGRID,NYGRID,NZGRID

  REAL, DIMENSION(MAXSTRM) :: XSTR,YSTR,ZSTR, DXSTR,DYSTR,DZSTR

  REAL, DIMENSION(MAXCUTS) :: XCUT,YCUT,ZCUT
  REAL, DIMENSION(MAXCUTS) :: NXCUT,NYCUT,NZCUT

  REAL,DIMENSION(4,MAXEDGES,MAXABUTS):: edgeDef   ! never used???


!!!   REAL,DIMENSION(MAXGRID):: x,y,z 
!
!!!   CHARACTER (LEN=80), DIMENSION(MAXNETS) :: NETNAMES
!----------------------------------------------------------------------------
  CALL Welcome()   ! opens the auxiliary file as unit IAUX

  CALL FindWgsFileName(IAUX, wgsFile)   ! Find the name of the WGS file 
  IF(Len_Trim(wgsFile)==0) THEN
    WRITE(*,*) "No WGS file specified in the input file."
    STOP
  ELSE
    WRITE(*,*) " Geometry data to be read from "//wgsFile
  END IF

!... Read the gridpoints from the WGS file ...........................
  OPEN(UNIT=IWGS, FILE=wgsFile, STATUS='OLD', &
    IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
  IF(errCode /= 0) THEN
    WRITE(*,*) "Cannot open the WGS file"
    STOP
  END IF
  WRITE(*,*) "Reading WGS file..."
  CALL ReadNetworks(IWGS,np,nets)
  CLOSE(UNIT=IWGS)
!  CALL PrintSize(DBG,netRoot)
!  CALL PrintLimits(DBG,netRoot)
!  CALL PrintData(DBG,netRoot)

!  nets=CountNetworksInList(netRoot)
  WRITE(DBG,*) "There are ", nets, " networks."
  DO k=1,nets
    np(k)%mnswitch=0
  END DO
  np(1:nets)%mnswitch=0
  netign(1:nets)=0
  netavg(1:nets)=0

!... Read the flight conditions reference dimensions, output codes,
!     boundary conditions, etc. from the auxiliary file .................
  WRITE(*,*) "Reading input file instructions..."
  CALL ProcessCommands()
  CLOSE(UNIT=IAUX)

!..... Write the A502 input file (all but the networks) ................
  a502file=Trim(AdjustL(name))//".in"
  OPEN(UNIT=I502, FILE=a502File, STATUS='REPLACE', &
    IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')
  IF (errCode /= 0) THEN
    WRITE(*,*) "Unable to open the new PanAir input file"
    STOP
  END IF
  WRITE(*,*) "Writing PanAir input file..."

  WRITE(I502, '(A)') '$TITLE    Created by panin, version '//   &
                                VERSION//' on '//dateTimeStr
  CALL WRITE502()   ! write the flight conditions, etc.

!..... Write the network gridpoints to the A502 file ...................
! (I502, netRoot, prec)
!  CALL PrintNetworks502(I502,netRoot,prec,kt,nlopt1,nropt1,nlopt2,nropt2,mnswitch)

  CALL WriteNets(I502, PREC, np(1:nets))   

!..... Write any wakes defined by WAKE records .........................
   CALL WriteWakes(I502, NWAKES, WAKE, WAKEEDGE, WAKEXTE)
   WRITE(I502, '(A)' ) '$END OF INPUT DATA'
   CLOSE(UNIT=I502)

!..... Say goodbye and ride off into the sunset ........................
  WRITE(*,*) " Files "//Trim(A502FILE)//" added to your directory."
  WRITE(*,*) "Also, file panin.dbg"
  WRITE(*, '(A)' ) FAREWELL
  WRITE(*,*) "Normal termination of panin"
  STOP

CONTAINS


!+
SUBROUTINE COPYDATA (IN, OUT)
!   --------------------------------------------------------------------
! PURPOSE - Copy data records from file attached as unit IN to the file
!   attached as unit OUT. All characters including and following a
!   comment character are changed to blanks. The comment characters
!   are ! = and /
!   Blank lines are not copied to OUT.
!
!
! NOTES- Special treatment for lines beginning with TITLE1 or TITLE2
!
  INTEGER,INTENT(IN):: IN   ! unit number of the  input file
  INTEGER,INTENT(IN):: OUT  ! unit number of the output file
                            ! these are assumed to be opened

  CHARACTER(LEN=1),PARAMETER:: BLANK = ' '
  CHARACTER(LEN=132):: buffer
  INTEGER:: errCode
  CHARACTER(LEN=4):: test
  INTEGER:: ngood=0  ! the number of data records
  INTEGER:: nrec=0   ! the total number of records (data & comments)
  INTEGER:: I
  INTEGER:: k
  CHARACTER(LEN=4),PARAMETER,DIMENSION(1):: WORDS = (/ "TITL" /)
!-----------------------------------------------------------------------
  DO
    READ(IN, '(A)', IOSTAT=errCode) buffer
    IF (errCode < 0) EXIT
    IF (errCode > 0) THEN
      WRITE(*,*) "Fatal error reading input file"
      STOP            ! that's really fatal
    END IF
    nrec=nrec+1
    IF (Len_Trim(buffer)==0) CYCLE

    TEST=AdjustL(buffer)            ! first 4 chars of buffer token
    k=TestKeyWord(test, WORDS)
    IF (k==0) THEN
      I=INDEX(BUFFER, '/')
      IF (i>0) BUFFER(I:)=BLANK

      I=INDEX(BUFFER, '=')
      IF (i>0) BUFFER(I:)=BLANK

      I=INDEX(BUFFER, '!')
      IF (i>0) BUFFER(I:)=BLANK

    END IF
    IF (Len_Trim(buffer) > 0) THEN
      WRITE(OUT,'(A)') Trim(buffer)
      ngood=ngood+1
    END IF
  END DO

  WRITE(*,*) NREC,  ' records copied from auxiliary file.'
  WRITE(*,*) ngood, ' records in the internal data file.'
  RETURN
END Subroutine CopyData   ! ----------------------------------------------------


!+
SUBROUTINE FindWgsFileName(IAUX, WGSFILE)
! ---------------------------------------------------------------------------
! PURPOSE - Find the name of the WGS file in the auxiliary file
!
  INTEGER,INTENT(IN):: IAUX
  CHARACTER(LEN=*),INTENT(OUT):: WGSFILE

  CHARACTER(LEN=3),PARAMETER,DIMENSION(1):: DICT = (/ "WGS" /)

  CHARACTER(LEN=80):: buffer
  INTEGER:: k
  INTEGER:: nfound
  INTEGER,DIMENSION(2):: startToken, endToken
!-----------------------------------------------------------------------
  DO
    READ(IAUX, '(A)', IOSTAT=errCode) buffer
    IF (errCode < 0) EXIT
    IF (errCode > 0) THEN
      WRITE(*,*) "Fatal error 1 looking for name of WGS file"
      STOP
    END IF

    CALL ParseLine(buffer, nfound, startToken,endToken)
    k=TestKeyWord(buffer(startToken(1):endToken(1)), DICT)
    IF (k==1) THEN
      IF (nfound < 2) THEN
        WRITE(*,*) "Fatal error 2 looking for name of WGS file"
        STOP
      END IF
      wgsFile=buffer(startToken(2):endToken(2))
      RETURN
    END IF
  END DO

  wgsFile=" "
  REWIND(UNIT=IAUX)
  RETURN
END Subroutine FindWgsFileName   ! ------------------------------------------


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
SUBROUTINE ParseLine(string, nfound, startToken,endToken)
! ---------------------------------------------------------------------------
! PURPOSE - Parse a line to delimit the tokens
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
!          from Sergio Gelato
!
!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 19Apr82  0.1   RLC   Original coding
!  6May83  1.0   RLC   Replaced SCAN_STRING with SCAN
!  4Dec84  1.1   RLC   Added test for INC>NFOUND
! 28Jul86  1.2   RLC   IMPLICIT NONE
!  2Dec91  1.3   RLC   END DO
!  4Feb95  1.4   RLC   Renamed Scan to ScanString
!  1Jun95  1.5   RLC   Added LenString (N2 parameter in ScanString)
!                         SCAN is an intrinsic in Fortran 90
! 10Jul96  1.6   RLC   Recoded, inspired by a code fragment posted to
!                         comp.lang.fortran by Sergio Gelato
!                         gelato@astrosun.tn.cornell.edu
!                      Also, omitted nmax from args. Use length of start,end

  CHARACTER(LEN=*),INTENT(IN):: string   ! character variable to be analyzed
  INTEGER,INTENT(OUT):: nfound           ! number of tokens found in string
  INTEGER,INTENT(OUT),DIMENSION(:):: startToken ! starting position 
                                                ! (in string) of each token
  INTEGER,INTENT(OUT),DIMENSION(:):: endToken ! ending position (in string)
                                              ! of each token
                               ! token i is string(startToken(i):endToken(i))

  INTEGER:: k,nmax
  LOGICAL:: inword
!----------------------------------------------------------------------------
  inword=.FALSE.
  nfound=0
  nmax=MIN(SIZE(startToken), SIZE(endToken))
  IF (nmax <= 0) RETURN

  DO k=1,LEN_TRIM(string)
    IF (string(k:k) == " ") THEN
      IF (inword) THEN
        endToken(nfound)=k-1
        IF (nfound == nmax) RETURN
        inword=.FALSE.
      END IF
    ELSE
      IF (.NOT.inword) THEN
        IF (nfound+1 > nmax) RETURN
        nfound=nfound+1
        startToken(nfound)=k
        inword=.TRUE.
      END IF
    END IF
  END DO
  IF (inword) endToken(nFound)=LEN_TRIM(string)

  RETURN
END Subroutine ParseLine   ! ------------------------------------------------


!+
SUBROUTINE ProcessCommands()
! ---------------------------------------------------------------------------
! PURPOSE - Read the auxiliary file and update the values of the global
!   control variables. In some cases, an external file must be opened
!   and read.
!
  INTEGER,PARAMETER:: ITMP=7         ! unit number of the file of XYZ points
  INTEGER,PARAMETER:: NMAX=40           ! maximum number of tokens on an input line


  CHARACTER(LEN=9),PARAMETER,DIMENSION(63):: WORDS = (/                 &
     'MACH     ', 'SREF     ', 'XREF     ', 'YREF     ', 'ZREF     ',   &
     'CBAR     ', 'SPAN     ', 'CASES    ', 'ALPC     ', 'BETC     ',   &
     'ALPHA    ', 'BETA     ', 'SYMM     ', 'ISINGS   ', 'IGEOMP   ',   &
     'ISINGP   ', 'ICONTP   ', 'IBCONP   ', 'IEDGEP   ', 'IPRAIC   ',   &
     'NEXDGN   ', 'IOUTPR   ', 'IFMCPR   ', 'IVCORR   ', 'HMIN     ',   &
     'HMAX     ', 'EAT      ', 'IGEOIN   ', 'IGEOUT   ', 'ABUTPR   ',   &
     'TITLE1   ', 'TITLE2   ', 'XYZ      ', 'STREAM   ', 'BOUN     ',   &
     'FORCE    ', 'GRID     ', 'IGNORE   ', 'AVERAGE  ', 'MNSWITCH ',   &
     'MNSWCH   ', 'PRECISION', 'SAVE     ', 'RESTART  ', 'NAME     ',   &
     'TIME     ', 'CPF      ', 'CHECK    ', 'SUBMITTER', 'WGS      ',   &
     'PATH     ', 'TAPE13   ', 'CUTS     ', 'MEMORY   ', 'MEM      ',   &
     'DIVERT   ', 'OFFNETS  ', 'WAKE     ', 'MAXSTREAM', 'UPSTREAM ',   &
     'NOF      ', 'BWAKE    ', 'PEA      ' /)



  CHARACTER(LEN=255):: arg1,buffer
  INTEGER:: NFOUND
  INTEGER:: I,J,K

  INTEGER,DIMENSION(NMAX)::  startToken,endToken
  REAL:: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
  INTEGER:: itmp1,itmp2,itmp3
!-----------------------------------------------------------------------

  DO
    READ(IAUX, '(A)', IOSTAT=errCode) buffer
    IF (errCode < 0) RETURN
    CALL ParseLine(buffer, nfound, startToken, endToken)
    IF (nfound==0) CYCLE
    k=TestKeyWord(buffer(startToken(1):endToken(1)), WORDS)
    IF (k==0) THEN
      WRITE(*,*) "***** Comment: "//Trim(buffer)
      WRITE(DBG,*) "***** Comment: "//Trim(buffer)
      CYCLE
    ELSE
      WRITE(*, '(A,I2,1X,A)' ) ' Command ', k, Trim(buffer)
      WRITE(DBG, '(A,I2,1X,A)' ) ' Command ', k, Trim(buffer)
      IF (nfound > 1) arg1=buffer(startToken(2):endToken(2))
    END IF

!..... Input record is a valid command. Process it......................
  SELECT CASE(k)
    CASE(1)
      READ(arg1, '(F80.0)' ) MACH               ! MACH
      WRITE(DBG,*) "Mach=", mach
    CASE(2)
      READ(arg1, '(F80.0)' ) SREF               ! SREF
    CASE(3)
      READ(arg1, '(F80.0)' ) XREF               ! XREF
    CASE(4)
      READ(arg1, '(F80.0)' ) YREF               ! YREF
    CASE(5)
      READ(arg1, '(F80.0)' ) ZREF               ! ZREF
    CASE(6)
      READ(arg1, '(F80.0)' ) CBAR               ! CBAR
    CASE(7)
      READ(arg1, '(F80.0)' ) SPAN               ! SPAN
    CASE(8)
      
    CASE(9)
      READ(arg1, '(F80.0)' ) ALPC               ! ALPC
    CASE(10)
      READ(arg1, '(F80.0)' ) BETC               ! BETC
    CASE(11)
      CASES=MIN(4, NFOUND-1)                                          ! ALPHA
      DO I=1,CASES
         READ(buffer(startToken(I+1):endToken(I+1)), '(F80.0)' ) ALPHA(I)
      END DO
    CASE(12)
      CASES=MIN(4, NFOUND-1)                                           ! BETA
      DO I=1,CASES
         READ(buffer(startToken(I+1):endToken(I+1)), '(F80.0)' ) BETA(I)
      END DO
    CASE(13)
       READ(arg1, '(I80)' ) I                    ! SYMM
       IF (I == 0) THEN
          MISYMM=0
          MJSYMM=0
       ELSE IF (I == 2) THEN
          MISYMM=1
          MJSYMM=1
       ELSE
          MISYMM=1
          MJSYMM=0
       END IF
    CASE(14)
      READ(arg1, '(I80)' ) OUTCODES(1)        ! ISINGS
    CASE(15)
      READ(arg1, '(I80)' ) OUTCODES(2)        ! IGEOMP
    CASE(16)
      READ(arg1, '(I80)' ) OUTCODES(3)        ! ISINGP
    CASE(17)
      READ(arg1, '(I80)' ) OUTCODES(4)        ! ICONTP
    CASE(18)
      READ(arg1, '(I80)' ) OUTCODES(5)        ! IBCONP
    CASE(19)
      READ(arg1, '(I80)' ) OUTCODES(6)        ! IEDGEP
    CASE(20)
      READ(arg1, '(I80)' ) OUTCODES(7)        ! IPRAIC
    CASE(21)
      READ(arg1, '(I80)' ) OUTCODES(8)        ! NEXDGN
    CASE(22)
      READ(arg1, '(I80)' ) OUTCODES(9)        ! IOUTPR
    CASE(23)
      READ(arg1, '(I80)' ) OUTCODES(10)       ! IFMCPR
    CASE(24)
      READ(arg1, '(I80)' ) OUTCODES(12)       ! IVCORR
    CASE(25)
      READ(arg1, '(F80.0)' ) HMIN               ! HMIN
    CASE(26)
      READ(arg1, '(F80.0)' ) HMAX               ! HMAX
    CASE(27)
      READ(arg1, '(F80.0)' ) TOL           ! TOLERANCE
    CASE(28)
      READ(arg1, '(I80)' ) OUTCODES(13)       ! IGEOIN
    CASE(29)
      READ(arg1, '(I80)' ) OUTCODES(14)       ! IGEOUT
    CASE(30)
      READ(arg1, '(I80)' ) OUTCODES(15)       ! ABUTPR
    CASE(31)
      IF (NFOUND .GT. 1) THEN                                   ! TITLE1
        TITLE1=buffer(startToken(2):endToken(NFOUND))
      ELSE
        TITLE1=' '
      END IF
    CASE(32)
      IF (NFOUND .GT. 1) THEN                                   ! TITLE2
        TITLE2=buffer(startToken(2):endToken(NFOUND))
      ELSE
        TITLE2=' '
      END IF

    CASE(33)                                                       ! XYZ
      OPEN(UNIT=ITMP, FILE=arg1, STATUS='OLD', &
        IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
      IF (errCode /= 0) THEN
        WRITE(*,*) "Unable to find XYZ file "//Trim(arg1)
        CYCLE
      END IF
      DO
        READ(ITMP, *, IOSTAT=errCode) tmp1,tmp2,tmp3
        IF (errCode /= 0) EXIT
        IF (nxyz >= MAXOFB) EXIT
        NXYZ=NXYZ+1
        xoff(nxyz)=tmp1
        yoff(nxyz)=tmp2
        zoff(nxyz)=tmp3
      END DO
      CLOSE(UNIT=ITMP)

    CASE(34)                                                  ! STREAM   
      OPEN(UNIT=ITMP, FILE=arg1, STATUS='OLD', &
        IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
      IF (errCode /= 0) THEN
        WRITE(*,*) "Unable to find streamline file "//Trim(arg1)
        CYCLE
      END IF
      DO
        READ(ITMP,*, IOSTAT=errCode) tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
        IF (errCode /= 0) EXIT
        IF (numpts >= MAXSTRM) EXIT
        NUMPTS=NUMPTS+1
          XSTR(NUMPTS)=tmp1
          YSTR(NUMPTS)=tmp2
          ZSTR(NUMPTS)=tmp3
          DXSTR(NUMPTS)=tmp4
          DYSTR(NUMPTS)=tmp5
          DZSTR(NUMPTS)=tmp6
      END DO
      CLOSE(UNIT=ITMP)

    CASE(35)                                                      ! BOUN
      I=0
      J=2
      DO i=1,nets
        IF (j > nfound) THEN
          READ(IAUX, '(A)') buffer
          CALL ParseLine(buffer, nfound, startToken, endToken)
          j=1
        END IF
        READ(buffer(startToken(j):endToken(j)), '(I80)' ) np(i)%kt
        j=j+1
        IF (np(i)%kt .EQ. 30) THEN
          J=J+1
          READ(buffer(startToken(J):endToken(J)), '(I80)' ) np(i)%nlopt1
          J=J+1
          READ(buffer(startToken(J):endToken(J)), '(I80)' ) np(i)%nropt1
          J=J+1
          READ(buffer(startToken(J):endToken(J)), '(I80)' ) np(i)%nlopt2
          J=J+1
          READ(buffer(startToken(J):endToken(J)), '(I80)' ) np(i)%nropt2
!          IF (NROPT1(I)==1 .OR. NROPT1(I)==7 .OR. NROPT1(I)==7) THEN
!            J=J+1
!            READ(buffer(startToken(J):endToken(J)), '(F80.0)' ) BETA1(I)
!          END IF
!          IF(NROPT2(I)==1 .OR. NROPT2(I)==7 .OR. NROPT2(I)==7) THEN
!            J=J+1
!            READ(buffer(startToken(J):endToken(J)), '(F80.0)' ) BETA2(I)
!          END IF
          WRITE(DBG,*) "Net ", i, "   kt=30 and ", &
            np(i)%nlopt1, np(i)%nropt1, np(i)%nlopt2, np(i)%nropt2
        ELSE
          WRITE(DBG,*) "Net ", i, "   kt=", np(i)%kt
        END IF
      END DO

    CASE(36)                                                     ! FORCE

    CASE(37)
      IF (ngrids > MAXOFB) CYCLE
      NGRIDS=NGRIDS+1                                             ! GRID
      READ(buffer(startToken(2):endToken(2)),  '(F80.0)' ) XORIGIN(NGRIDS)
      READ(buffer(startToken(3):endToken(3)),  '(F80.0)' ) YORIGIN(NGRIDS)
      READ(buffer(startToken(4):endToken(4)),  '(F80.0)' ) ZORIGIN(NGRIDS)
      READ(buffer(startToken(5):endToken(5)),  '(F80.0)' )  DXGRID(NGRIDS)
      READ(buffer(startToken(6):endToken(6)),  '(F80.0)' )  DYGRID(NGRIDS)
      READ(buffer(startToken(7):endToken(7)),  '(F80.0)' )  DZGRID(NGRIDS)
      READ(buffer(startToken(8):endToken(8)),  '(I80)'   )  NXGRID(NGRIDS)
      READ(buffer(startToken(9):endToken(9)),  '(I80)'   )  NYGRID(NGRIDS)
      READ(buffer(startToken(10):endToken(10)),'(I80)'   )  NZGRID(NGRIDS)

    CASE(38)                                                     ! IGNORE
      DO i=2,nfound
        IF (netsIgn >= MAXNETS) EXIT
        NETSIGN=NETSIGN+1
        READ(buffer(startToken(I):endToken(I)), '(I80)' ) NETIGN(NETSIGN)
      END DO    
    CASE(39)                                                    ! AVERAGE
      DO i=2,nfound
        IF (netsAvg >= MAXNETS) EXIT
        NETSAVG=NETSAVG+1
        READ(buffer(startToken(i):endToken(i)), '(I80)' ) NETAVG(NETSAVG)
      END DO  
    
    CASE(40:41)           ! MNSWITCH
      DO i=2,nfound
        READ(buffer(startToken(i):endToken(i)), '(I80)' ) J
        IF (j>0 .AND. j<= MAXNETS) np(j)%MNSWITCH=1
      END DO    
    CASE(42)
      READ(arg1, '(I80)' ) PREC            ! PRECISION
    CASE(43)
!!!      SAVE=.TRUE.                                                 ! SAVE
    CASE(44)
      READ(arg1, '(I80)' ) RESTART           ! RESTART
    CASE(45)
!      IF (NFOUND.GT.1) THEN                                       ! NAME
!         NAME=arg1
   !!!      NNAME=endToken(2)-startToken(2)+1
!      END IF
    CASE(46)
    CASE(47)
    CASE(48)
      READ(arg1, '(I80)' ) DATACHK           ! DATACHK
    CASE(49)
    CASE(50)      
    CASE(51)
    CASE(52)
!!!      GET13=.TRUE.                                              ! TAPE13
!
    CASE(53)                                                      ! CUTS
      OPEN(UNIT=ITMP, FILE=arg1, STATUS='OLD', &
        IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
      IF(errCode /= 0) THEN
        WRITE(*,*) "Unable to find CUTS file "//Trim(arg1)
        CYCLE
      END IF
      DO
        READ(ITMP,*,IOSTAT=errCode)  tmp1,tmp2,tmp3, tmp4,tmp5,tmp6
        IF (errCode < 0) EXIT
        IF (ncuts >= MAXCUTS) EXIT
        NCUTS=NCUTS+1
        XCUT(NCUTS)=tmp1
        YCUT(NCUTS)=tmp2
        ZCUT(NCUTS)=tmp3
        NXCUT(NCUTS)=tmp4
        NYCUT(NCUTS)=tmp5
        NZCUT(NCUTS)=tmp6
      END DO
      CLOSE(UNIT=ITMP)
    CASE(54)                                                    ! MEMORY
    CASE(55)                                                       ! MEM
    CASE(56)
    CASE(57)                                               ! OFFNETS
      DO i=2,nfound                                     
        READ(buffer(startToken(I):endToken(I)), '(I80)',IOSTAT=errCode) itmp1
        IF (errCode /= 0) EXIT
        IF (netsOff >= MAXNETS) EXIT
        NETSOFF=NETSOFF+1
        offNets(netsoff)=itmp1
      END DO

    CASE(58)
      IF (NFOUND .LT. 4) THEN                                  ! WAKE
         WRITE(*,*) 'Insufficient parameters on WAKE record'
      ELSE
        READ(arg1, '(I80)', IOSTAT=errCode) itmp1
        IF (errCode /= 0) CYCLE
        READ(buffer(startToken(3):endToken(3)), '(I80)', IOSTAT=errCode) itmp2
        IF (errCode /= 0) CYCLE
        READ(buffer(startToken(4):endToken(4)), '(F80.0)', IOSTAT=errCode) tmp1
        IF (errCode /= 0) CYCLE
        IF (nwakes >= MAXNETS) CYCLE
        nwakes=nwakes+1
        wake(nwakes)=itmp1
        wakeEdge(nwakes)=itmp2
        wakeXte(nwakes)=tmp1
      END IF

    CASE(59)
      READ(arg1, '(F80.0)' ) pointsPerStreamline ! MAXSTREAM
    CASE(60)
      direction=-1.0                                              ! UPSTREAM
    CASE(61)
      writeNOF=.TRUE.                                              ! NOF

    CASE(62)
      IF (NFOUND .LT. 4) THEN                                        ! BWAKE
         WRITE(*,*) 'Insufficient parameters on BWAKE record'
      ELSE
        NBWAKES=NBWAKES+1
        READ(arg1, '(I80)', IOSTAT=errCode) itmp1
        IF (errCode /= 0) CYCLE
        READ(buffer(startToken(3):endToken(3)), '(I80)', IOSTAT=errCode) itmp2
        IF (errCode /= 0) CYCLE
        READ(buffer(startToken(4):endToken(4)), '(F80.0)', IOSTAT=errCode) tmp1
        IF (errCode /= 0) CYCLE
        IF (nbWakes >= MAXNETS) CYCLE
        nbwakes=nbwakes+1
        bWake(nbwakes)=itmp1
        bWakeEdge(nbwakes)=itmp2
        bWakeXte(nbwakes)=tmp1
      END IF

    CASE(63)


  END SELECT

  END DO

  RETURN
END Subroutine ProcessCommands   ! --------------------------------------------------

!+
SUBROUTINE ReadNetworks(efu, np, nets)
! ---------------------------------------------------------------------------
  INTEGER,INTENT(IN):: efu
  TYPE(PanAirNetwork),INTENT(OUT),DIMENSION(:):: np   
  INTEGER,INTENT(OUT):: nets

  INTEGER,PARAMETER:: MAXGRID = 10000    ! a kludge. Get rid of this.
  INTEGER:: errCode
  INTEGER:: i,k
  INTEGER:: nobj
  CHARACTER(LEN=132):: title
  INTEGER:: n
  INTEGER:: localRows, localCols
  REAL,ALLOCATABLE,DIMENSION(:):: inpX,inpY,inpZ 
!----------------------------------------------------------------------------
  REWIND (efu)
  READ(efu, '(A)' ) title    ! read the general title 
  WRITE(DBG,*) "Reading wgs file general title:"//Trim(title)

  ALLOCATE(inpX(maxGrid), inpY(maxGrid), inpZ(maxGrid) )

  nets=0
  DO k=1,SIZE(np)
    READ(efu, '(A)', IOSTAT=errCode) title
    IF (errCode /= 0) RETURN
    CALL StripFortranQuotes(title)
    np(k)%name = title
    WRITE(*,*) "Reading network "//Trim(title)
    READ(efu,*) nobj, localCols, localRows
    np(k)%id = nobj
    np(k)%cols = localCols
    np(k)%rows = localRows
    n=localCols*localRows
    Read(efu,*) (inpX(i), inpY(i), inpZ(i), i=1,n)
    !!! now load np
    ALLOCATE(np(k)%x(n))
    ALLOCATE(np(k)%y(n))
    ALLOCATE(np(k)%z(n))
    np(k)%x=inpX(1:n)
    np(k)%y=inpY(1:n)
    np(k)%z=inpZ(1:n)
    WRITE(DBG,*) "Completed reading network:"//Trim(title)
    nets=nets+1
  END DO

  DEALLOCATE(inpX, inpY, inpZ)
  RETURN
END Subroutine ReadNetworks   ! =============================================

!+
SUBROUTINE StripFortranQuotes(t)
! ---------------------------------------------------------------------------
! PURPOSE - If the first non-blank character and the last non-blank
!   character in a string are both APOSTROPHE, then change them to
!   blanks and adjust the resulting string to the left.

  CHARACTER(LEN=*),INTENT(IN OUT):: t

  CHARACTER(LEN=1),PARAMETER:: APOSTROPHE = "'" ! 
!  CHARACTER(LEN=LEN(t)):: s   ! elegant, works on Lahey, but not Absoft
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
FUNCTION TestKeyWord(test,dict) RESULT(k)
! ---------------------------------------------------------------------------
  CHARACTER(LEN=*),INTENT(IN):: test
  CHARACTER(LEN=*),INTENT(IN),DIMENSION(:):: dict
  INTEGER:: k

!  CHARACTER(LEN=WORDLENGTH):: testWord
  CHARACTER(LEN=LEN(dict(1))):: testWord
  INTEGER:: j
!----------------------------------------------------------------------------
  k=0
  testWord=test         ! makes all the tests between strings of equal length
  CALL UpCase(testWord) ! assumes dict entries are all upper case

  DO j=1,SIZE(dict)
    IF (testWord == dict(j)) THEN
      k=j
      RETURN
    END IF
  END DO
  RETURN
END Function TestKeyWord   ! --------------------------------------------------

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




!+
SUBROUTINE WRITE502()
! ---------------------------------------------------------------------------
! PURPOSE - Write all data to the A502 input file except the grid points.

!!!   INTEGER, INTENT(IN) :: I502  ! unit number of the A502 input file

  REAL,PARAMETER:: ZERO=0.0, ONE=1.0
   INTEGER :: i
   CHARACTER (LEN=8) :: FMT = '(6F10. )'
   CHARACTER(LEN=*),PARAMETER:: FMT1="(7F10.4)", FMT2="(A/7F10.4)"
!-----------------------------------------------------------------------
   WRITE(FMT(7:7), '(I1)' ) PREC
   WRITE (I502, '(A)' ) TITLE1
   WRITE (I502, '(A)' ) TITLE2
   WRITE (I502,FMT2) "$DATACHECK", REAL(DATACHK)
   IF (RESTART > 0) WRITE(I502,'(A)') "$RESTART"
   IF (RESTART == 1) THEN
      WRITE (I502, '(3F10.1)' ) 0.0, 1.0, 1.0
   ELSE IF (RESTART == 2) THEN
      WRITE (I502, '(3F10.1)' ) 1.0, 1.0, 1.0
   ELSE IF (RESTART == 3) THEN
      WRITE (I502, '(3F10.1)' ) 1.0, 0.0, 1.0
   END IF
!
!..... Flow properties .................................................
   WRITE(I502,FMT2) '$SYM', REAL(MISYMM), REAL(MJSYMM)
   WRITE(I502,FMT2) '$MACH NUMBER', MACH
   WRITE(I502,FMT2) '$CASES', REAL(CASES)
   WRITE(I502,FMT2) '$ANGLES OF ATTACK', ALPC
   WRITE(I502,FMT1) alpha(1:cases) 
   WRITE(I502,FMT2) '$YAW ANGLES', BETC
   WRITE(I502,FMT1) beta(1:cases) 

!..... Reference dimensions ............................................
  WRITE(I502,'(A)') "$REFERENCE DATA FOR 3-D FORCES AND MOMENTS"
  WRITE(I502,FMT1) XREF,YREF,ZREF, ONE
  WRITE(I502,FMT) SREF, SPAN, CBAR, SPAN

!..... Printout control.................................................
   WRITE(I502,'(A)') "$PRINTOUT CONTROL"
   WRITE(I502, '(6F10.1)' ) (REAL(OUTCODES(I)),I=1,12)
   WRITE(I502,'(A)') "$EAT"
   WRITE(I502, '(F10.6,5F10.1)' ) TOL, (REAL(OUTCODES(I)),I=13,17)
   WRITE(I502,FMT2) "$VELOCITY CORRECTION", REAL(IVCORR)

!... Sectional properties ............................................
  IF (NCUTS .GT. 0) THEN
    WRITE(I502,'(A)') "$SECTIONAL PROPERTIES"
    WRITE(I502,'(A)') "1.0"
    WRITE(I502,'(A)') "*CUT"
    WRITE(I502, '(7F10.1)' ) 0.0,0.0,1.0,1.0,0.0,0.0,1.0
    WRITE(I502, '(F10.1)' ) REAL(NCUTS)
    WRITE(I502,fmt) (XCUT(I), YCUT(I), ZCUT(I),   &
         NXCUT(I), NYCUT(I), NZCUT(I), I=1,NCUTS)
  END IF

!... NO-filaments for Trefftz plane drag analysis ....................
   IF (writeNOF) WRITE(I502, '(A)' ) '$NOF'

!... Flow field ......................................................
   IF (NXYZ.GT.0 .OR. NUMPTS.GT.0 .OR. NGRIDS.GT.0 .OR. NETSOFF.GT.0)   &
     WRITE(I502, '(A/F10.1)' )  '$FLOW-FIELDS', 1.0
!
!..... Network selection ...............................................
   IF (NETSOFF .GT. 0) THEN
      WRITE(I502, '(A)' )  '$NWL NETWORK SELECTION FOR FLOW-FIELD PROP'
      WRITE(I502, '(F10.1)' ) REAL(NETSOFF)
      WRITE(I502,fmt) (REAL(OFFNETS(I)), I=1,NETSOFF)
   END IF
!
!..... Offbody points...................................................
   IF (NXYZ .GT. 0) THEN
      WRITE(I502, '(A)' ) '$XYZ'
      WRITE(I502, '(F10.1)' ) REAL(NXYZ)
      WRITE(I502, FMT) (XOFF(I), YOFF(I), ZOFF(I), I=1,NXYZ)
   END IF
!
   IF (NGRIDS .GT. 0) THEN
      WRITE(I502, '(A/F10.1)' ) '$GRID', REAL(NGRIDS)
      DO I=1,NGRIDS
         WRITE(I502, FMT)                             &
                    XORIGIN(I),       YORIGIN(I),      ZORIGIN(I),   &
                1.0+XORIGIN(I),       YORIGIN(I),      ZORIGIN(I),   &
                    XORIGIN(I),   1.0+YORIGIN(I),      ZORIGIN(I),   &
                    XORIGIN(I),       YORIGIN(I),  1.0+ZORIGIN(I),   &
                     DXGRID(I),        DYGRID(I),       DZGRID(I),   &
              REAL(NXGRID(I)), REAL(NYGRID(I)), REAL(NZGRID(I))
      END DO
   END IF
!
!..... Streamlines......................................................
   IF (NUMPTS .GT. 0) THEN
      WRITE(I502, '(A)' ) '$STREAMLINES', '1.0'
      WRITE(I502, '(F10.1,F10.5,F10.2,2F10.1,F10.4)' ) REAL(NUMPTS),   &
                      HMIN, HMAX, pointsPerStreamline, MAXORDR, ABSERR
      WRITE(I502, '(3F10.1)' ) 10000.0,  0.0,  0.0

      DO I=1,NUMPTS
         WRITE(I502, FMT)  XSTR(I), YSTR(I), ZSTR(I),   &
                          DXSTR(I),DYSTR(I),DZSTR(I), direction
      END DO
   END IF
!
!..... Network selection for force & moment calculations................
   IF (NETSIGN.GT.0 .OR. NETSAVG.GT.0) WRITE(I502, '(A)' )   &
      '$FORCE AND MOMENT SELECTION'
   IF (NETSIGN .GT. 0) THEN
      WRITE(I502, '(A)' ) '*NETWORK SELECTION'
      WRITE(I502, '(2F10.1)' ) REAL(-NETSIGN), REAL(IDEFAULT)
      WRITE(I502, '(F10.1)' ) (NETIGN(I),I=1,NETSIGN)
   END IF
!
   IF (NETSAVG .GT. 0) THEN
      WRITE(I502, '(A)' ) '*NETWORK SELECTION'
      WRITE(I502, '(2F10.1)' ) REAL(NETSAVG), REAL(IDEFAULT)
      WRITE(I502, '(F10.1)' ) (NETAVG(I),I=1,NETSAVG)
      WRITE(I502, '(A)' ) '*PRESSURE SURFACE'
      WRITE(I502, '(2F10.1)' ) REAL(NETSAVG), REAL(IDEFAULT)
      WRITE(I502, '(F10.1)' ) (NETAVG(I),I=1,NETSAVG)
   END IF
   RETURN
END Subroutine Write502   ! -------------------------------------------------

!+
SUBROUTINE WriteWakes(I502, NWAKES, WAKE,WAKEEDGE,WAKEXTE)
! ---------------------------------------------------------------------------
!  PURPOSE - Write the wake definitions to the 502 input file.

   INTEGER,INTENT(IN):: I502
   INTEGER,INTENT(IN):: NWAKES
   INTEGER,INTENT(IN),DIMENSION(:):: WAKE
   INTEGER,INTENT(IN),DIMENSION(:):: WAKEEDGE
   REAL,INTENT(IN),DIMENSION(:):: WAKEXTE

  CHARACTER(LEN=*),PARAMETER:: FMT = "(2F10.1,F10.3)"
  INTEGER:: k
!----------------------------------------------------------------------------
  IF (NWAKES <= 0) RETURN

  WRITE(I502,'(A)') "$TRAILING  WAKES FROM WINGS"
  WRITE(I502,FMT) REAL(NWAKES)
  WRITE(I502,'(A)') "18.0"

  DO K=1,NWAKES
    WRITE(I502,FMT) REAL(WAKE(K)), REAL(WAKEEDGE(K)), WAKEXTE(K)
  END DO
  RETURN
END Subroutine WriteWakes   ! -------------------------------------------------

!+
SUBROUTINE Welcome()
! ---------------------------------------------------------------------------
! PURPOSE -

  CHARACTER(LEN=*),PARAMETER :: GREETING="Prepare input for PanAir"
  CHARACTER(LEN=*),PARAMETER:: AUTHOR= &
    "Ralph L. Carmichael, Public Domain Aeronautical Software"

  CHARACTER(LEN=80):: auxFile
  INTEGER:: errCode
!----------------------------------------------------------------------------
!!    CALL TIMEID(TIMESTR)    DateTimeStr
   dateTimeStr=" "

  WRITE(*,*) GREETING
  WRITE(*,*) ' Version '//VERSION
  WRITE(*,*) AUTHOR

  OPEN(UNIT=DBG, FILE="panin.dbg", STATUS='REPLACE', &
    IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')
  IF (errCode==0) THEN
    WRITE(DBG,*) "Debug output from panin on ", dateTimeStr
  ELSE
    WRITE(*,*) "Unable to open debugging file"
    STOP
  END IF
!
!
!..... Get the name of the auxiliary file and copy contents to IAUX ..........
  DO
    WRITE(*,*) "Enter the name of the auxiliary file: "
    READ(*, '(A)') AUXFILE
    IF (Len_trim(auxFile)==0) STOP
    OPEN(UNIT=IauxTemp, FILE=auxFile, STATUS='OLD',   &
      IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
    IF(errCode==0) EXIT
    OPEN(UNIT=IauxTemp, FILE=Trim(auxFile)//'.aux', STATUS='OLD',   &
      IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
    IF(errCode==0) EXIT
    WRITE(*,*) "Unable to open this file. Try again."
  END DO

  OPEN(UNIT=IAUX, FILE='panin.tmp', STATUS='REPLACE', &
    IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')
  IF (errCode==0) THEN
    WRITE(DBG,*) "Reading auxiliary data from "//Trim(auxFile)
  ELSE
    WRITE(*,*) "Unable to open temp file for AUX"
    STOP
  END IF

!... Copy the contents of the auxiliary file to a new file, omitting
!    the comment records and blank records. Use this new file for
!    the actual calculations. 
  CALL COPYDATA(IAUXTEMP, IAUX)
  CLOSE(UNIT=IAUX)             ! re-open later for reading
  CLOSE(UNIT=IAUXTEMP)

!... Now re-open the copy of the auxiliary file and read the instructions ...
  OPEN(UNIT=IAUX, FILE='panin.tmp', STATUS='OLD', &
    IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
  IF (errCode /=0) THEN
    WRITE(*,*) "Unable to reopen temp file for AUX to read it"
    STOP
  END IF
  

  RETURN

END Subroutine Welcome   ! --------------------------------------------------

!+
SUBROUTINE WriteNets(I502, prec, np)
! ---------------------------------------------------------------------------
! PURPOSE - Write the gridpoints and boundary condition information
!   for each network to the results file.

  INTEGER,INTENT(IN):: I502
  INTEGER,INTENT(IN):: prec
  TYPE(PanAirNetwork),INTENT(IN),DIMENSION(:):: np

  REAL:: amn
  INTEGER:: OFFSET,I,J,K
  INTEGER:: N1,N2,N3
  INTEGER:: i1,i2
  REAL:: NTS,NTD
  CHARACTER(LEN=8):: FMT= '(6F10. )'
!-----------------------------------------------------------------------
  WRITE(FMT(7:7), '(I1)' ) PREC
  DO k=1,SIZE(np)
    amn=REAL(np(k)%mnSwitch)
    WRITE(I502, '(A,I3,3X,A)' ) '$POINTS for network #', k, np(k)%name(1:50)
    WRITE(I502, '(A)' ) '       1.0'  ! NETS PER GROUP =1
    IF (np(k)%kt == -18) THEN
      WRITE(I502, '(2F10.1,20X,2F10.1)' ) 18.0,1.0, amn, 0.0
    ELSE IF (np(k)%kt /= 30) THEN
      WRITE(I502, '(F10.1,30X,2F10.1)' ) REAL(np(k)%kt), amn, 0.0
    ELSE
      IF (np(k)%nlopt1==0) THEN
        NTS=0.0
      ELSE
        NTS=1.0
      END IF
      IF (np(k)%nlopt2==0) THEN
        NTD=0.0
      ELSE
        NTD=12.0
      END IF
      WRITE(I502, '(6F10.1)' ) REAL(np(k)%kt),nts,ntd,0.0,amn,0.0
      WRITE(I502, '(6F10.1)' ) REAL(np(k)%nlopt1), REAL(np(k)%nropt1), &
                               REAL(np(k)%nlopt2), REAL(np(k)%nropt2)
    END IF

    WRITE(I502, '(2F10.1,50X,A10)' ) &
      REAL(np(k)%rows), REAL(np(k)%cols), np(k)%name(1:10)

    DO J=1,np(k)%cols
      i1=1+(j-1)*np(k)%rows
      i2=j*np(k)%rows
      WRITE(I502, FMT)  (np(k)%x(i),np(k)%y(i),np(k)%z(i),i=i1,i2)
    END DO

    IF (np(k)%kt==4 .OR. np(k)%kt==9 .OR. np(k)%kt==14) THEN
      DO J=1,(np(k)%rows-1)*(np(k)%cols-1)
        WRITE(I502,'(A)' ) "       0.0       0.0       0.0      0.0"
      END DO
    END IF

!    IF (np(k)%kt==30 .AND.  &
!     &     (np(k)%nropt1==1 .OR. np(k)%nropt11==7 .OR. NROPT1(K).EQ.8))
!     &     WRITE(I502,'(F10.5)') (BETA1(K),I=1,(ROWS(K)-1)*(COLS(K)-1))

!         IF (KT(K).EQ.30 .AND.
!     &     (NROPT2(K).EQ.1 .OR. NROPT2(K).EQ.7 .OR. NROPT2(K).EQ.8))
!     &     WRITE(I502,'(F10.5)') (BETA2(K),I=1,(ROWS(K)+1)*(COLS(K)+1))

  END DO
  RETURN
END Subroutine WriteNets   ! ------------------------------------------------

!+
!SUBROUTINE WriteWakes(I502, NWAKES, WAKE,WAKEEDGE,WAKEXTE)
! ---------------------------------------------------------------------------
! PURPOSE - Write the wake definitions to the 502 input file.

!  INTEGER,INTENT(IN):: I502
!  INTEGER,INTENT(IN):: NWAKES
!  INTEGER,INTENT(IN),DIMENSION(:):: WAKE
!  INTEGER,INTENT(IN),DIMENSION(:):: WAKEEDGE
!  REAL,INTENT(IN),DIMENSION(:):: WAKEXTE

!  INTEGER:: K
!----------------------------------------------------------------------------
!  IF (NWAKES <= 0) RETURN
!
!  WRITE(I502, '(A)' ) "$TRAILING  WAKES FROM WINGS"
!  WRITE(I502, '(F10.1)' ) REAL(NWAKES)
!  WRITE(I502, '(A)' ) "18.0"
!
!  DO K=1,NWAKES
!    WRITE(I502, '(2F10.1, F10.3)' ) REAL(WAKE(K)),
!     &      REAL(WAKEEDGE(K)), WAKEXTE(K)
!  END DO
!  RETURN
!END Subroutine WriteWakes   ! -----------------------------------------------

!+
SUBROUTINE WriteBWakes(I502, NBWAKES, BWAKE,BWAKEEDGE,BWAKEXTE)
! ---------------------------------------------------------------------------
! PURPOSE - Write the wake definitions to the 502 input file.

  INTEGER,INTENT(IN):: I502
  INTEGER,INTENT(IN):: NBWAKES
  INTEGER,INTENT(IN),DIMENSION(:):: BWAKE
  INTEGER,INTENT(IN),DIMENSION(:):: BWAKEEDGE
  REAL,INTENT(IN),DIMENSION(:):: BWAKEXTE

  INTEGER:: K
!----------------------------------------------------------------------------
  IF (NBWAKES <= 0) RETURN

  WRITE(I502, '(A)') "$TRAILING  WAKES FROM BODIES, NACELLES OR UPSTREAM WAKES"
  WRITE(I502, '(F10.1)' ) REAL(NBWAKES)
  WRITE(I502, '(A)' ) "18.0      1.0"

  DO K=1,NBWAKES
    WRITE(I502, '(2F10.1, F10.3)' ) REAL(BWAKE(K)), &
      REAL(BWAKEEDGE(K)), BWAKEXTE(K)
  END DO
  RETURN
END Subroutine WriteBWakes   ! ----------------------------------------------

!+
SUBROUTINE WritePEA(I502, NABUTS, EDGES, EDGEDEF)
! ---------------------------------------------------------------------------
! PURPOSE - Write the partial edge abutment data to the 502 file

  INTEGER,INTENT(IN):: I502
  INTEGER,INTENT(IN):: NABUTS
!!!  INTEGER,INTENT(IN):: MAXEDGES
  INTEGER,INTENT(IN),DIMENSION(:):: EDGES
  INTEGER,INTENT(IN),DIMENSION(:,:,:):: EDGEDEF

  INTEGER:: I,J,K
!----------------------------------------------------------------------------
  IF (NABUTS > 0) THEN
    WRITE(I502, '(A)' ) "$PEA - PARTIAL OR FULL EDGE ABUTMENT"
    WRITE(I502, '(3F10.0)' ) REAL(NABUTS), 0.0, 0.0
    DO K=1,NABUTS
      WRITE(I502, '(F10.0)' ) REAL(EDGES(K))
      DO J=1,EDGES(K)
        WRITE(I502, '(4F10.0)' ) (REAL(EDGEDEF(I,J,K)),I=1,4)
      END DO
    END DO
  END IF

  RETURN
END Subroutine WritePEA   ! -------------------------------------------------

END Program PanAirInput   ! =================================================
