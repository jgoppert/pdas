!+
! PROGRAM HypersonicArbitraryBodyToLangleyWireframeGeometry
! ----------------------------------------------------------------------
! PURPOSE - Convert geometry data from the Hypersonic Arbitrary Body
!  Program into Langley Wireframe Geometry format. The goal is to read
!  old data files that were used for the SHABP and end up with Wgs files
!  that can be used in the new hypersonic program or in PanAir.

! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
! ORIGINAL AUTHORS (SHABP) -
!   Arvel E. Gentry, Douglas N. Smyth, Wayne R. Oliver, Douglas Aircraft
! REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
!   1973   xxx   AEG   Original release of SHABP, Mark 4
! 05Dec02  0.1   RLC   Original coding, based on main prog & GEOM in SHA

MODULE HypersonicArbitraryBodyToLangleyWireframeGeometryProcedures
! ----------------------------------------------------------------------
! PURPOSE - 

  CHARACTER(LEN=*),PARAMETER:: VERSION = "(0.2 22 July 2008)"

CONTAINS      

!+
SUBROUTINE GEOM()
! ----------------------------------------------------------------------
! PURPOSE - Control the execution of the geometry generation routines

!!!      include 'comp.cmn'
!
! *** This file defines the common block 'comp' which is used to save co
!     names, panel names (new longer version), and which components are
!     each save set.
!
!     maxcom = maximum number of components
!     maxppc = maximum points per component
!     maxpan = maximum number of panels
!     mxpapc = maximum number of panels/component
!     maxsav = maximum number of save sets that can be later summed

      INTEGER,PARAMETER:: MAXCOM=20,MXPAPC=10,MAXPAN=25,MAXSAV=50

!     comnam - array of component names (character array)
!     idpan  - array of panel ids for each component.  The first index i
!              the component, the second index is the panel numbers
!     npanls - array of the numbers of panels in each component
!     ncomps - is the number of components in file
!     pannam - array of panel names
!     savnam - name of component for each save set including deflection
!              angle of applicable

      CHARACTER(LEN=16),DIMENSION(MAXCOM):: comnam
      CHARACTER(LEN=16),DIMENSION(MAXPAN):: pannam
      CHARACTER(LEN=32),DIMENSION(MAXSAV):: savnam
      INTEGER,DIMENSION(MAXCOM,MXPAPC):: idpan


!      common /comp/ comnam(maxcom),idpan(maxcom,mxpapc),npanls(maxcom),
!     &              ncomps,pannam(maxpan),savnam(maxsav)

!     character comnam*16,pannam*16,savnam*32


!      COMMON /EXEC/CASE,TITLE,PAGE,ERROR
      !!! also needs /GDATA/ ??? RLC
      COMMON /GFLAG/IOUT,ISTAT3,IORIN,COMPIN
!     COMMON /TAPE/TAPEIN,TAPEOT,TAPEA,TAPEB,TAPEC,TAPED,TAPEE,TAPEF,   &
!     &             TAPEG,TAPEH,TAPEI,TAPEJ,TAPEK
!      COMMON /TAGS/ITAG4,ITAG9,ITAG10
!      common /units/ iunit
!
      DIMENSION  XA(250),XB(250),YA(250),YB(250),ZA(250),ZB(250),       &      
     & XI(4),ETA(4),XIN(4),YIN(4),ZIN(4),TITLE(15),XPA(4),YPA(4),ZPA(4)
      DIMENSION CONFIG(15),IGEOM(10),E(25),EO(25),EP(25)
      DIMENSION ELEM(25),PANEL(100),IORN(100),SYMFCT(100),IFA(100)
!
      REAL NX,NY,NZ
      LOGICAL  RFLAG, AFLAG, BFLAG
!      INTEGER  TAPEIN,TAPEOT,TAPEA,TAPEB,TAPEC,TAPED,TAPEE,TAPEF,       &  
!     &         TAPEG,TAPEH,TAPEI,TAPEJ,TAPEK
      INTEGER  STAT, STATT, PAGE, TYPE, ERROR, PRINTS, SYMFCT, SEQ

!----------------------------------------------
      ISTAT3 = 0
!
!  READ GEOMETRY CONTROL CARD
      READ(TAPEIN,10)IOUT,IREW,PRINTS,IQUAD,I3MAX,NEW,NPMAX,CONFIG,   &
     &                 iunit
   10  FORMAT (I2,5I1,I2,15A4,i1)
!  INITIALIZE ELEMENT DATA STORAGE UNIT IF REQUIRED
!      DO 13 I=1,25
!      E (I) = 0.0
!      EP(I) = 0.0
!   13  EO(I) = 0.0
!      IF (NEW .EQ. 1) GO TO 15
!      IG4 = 1
!     WRITE (4:IG4) CONFIG
!!!!      CALL WRITMS (4,CONFIG,15,IG4)
      IF (NPMAX .EQ. 0) NPMAX = 50
!  NREM IS THE NUMBER OF RECORDS PER ELEMENT
!      NREM = 3
!      NEXT = (NPMAX+1)*5
!      NP = 0
!      E(1) = NEXT
!      E(2) = NP
!      E(3) = NPMAX
!      E(4) = NREM
!      GO TO 19
!   15 IG4 = 1
!     READ (4:IG4) CONFIG
!!!!      CALL READMS (4,CONFIG,15,IG4)
 !     IG4 = 2
!     READ (4:IG4) E
!!!!      CALL READMS (4,E,25,IG4)
!      NEXT = E(1)
!      NP   = E(2)
!      NPMAX= E(3)
!      NREM = E(4)
!   19 CONTINUE
!
!  CHECK IF QUADRILATERALS ONLY ARE TO BE CALCULATED
!      IF (IQUAD .EQ. 1) GO TO 110
!      IF (IREW.EQ.0 .AND. IOUT.NE.5 .AND. IOUT.NE.1) REWIND IOUT
!  READ Panel IDENTIFICATION CARD
   20 ISTAT3 = ISTAT3 + 1
      READ (TAPEIN,30) PANEL(ISTAT3),LAST,IORN(ISTAT3),IGEOM,            &
     &    SYMFCT(ISTAT3),IFA(ISTAT3),NADJ1,NADJ2,NADJ3,NADJ4,ivis_typ,   &
     &    pannam(istat3)
   30  FORMAT (A4,2I1,10I1,2I1,4I2,i1,2x,a)
      IORIEN = IORN(ISTAT3)
      COMPIN = PANEL(ISTAT3)
!
!
! ****************************
!  READ OR GENERATE ELEMENTS FOR EACH VEHICLE Panel
      DO 90 I=1,10
       IF (IGEOM(I) .EQ. 0) GO TO 100
      IG = IGEOM(I)
      GO TO (40,50,60,70), IG
   40 CALL IELE
      GO TO 80
!
   50 CALL ELLIP
      GO TO 80
   60 CALL CUBIC
      GO TO 80
                                        
   70 WRITE(*,*) 'disabled. use wd2wgs'    !!! CALL AIRCFT
!
   80 IF (ERROR .NE. 0) THEN
        WRITE(*,*) 'Fatal error occured'
        STOP
      END IF  
   90 END DO
!
  100 IF (LAST .EQ. 0) GO TO 20
  440 FORMAT (/,' ***** SURFACE DATA ROUTINE HAS ATTEMPTED TO READ A'   &
     &' NON SURFACE CARD - CHECK YOUR CARDS ***** ' ,a4,i4)
      ERROR = 1
      GO TO 490
!
  450 FORMAT (/,7X, I4, 1P4E14.5,0PF10.6,1P2E14.5,I6,2(/12X,4E14.5,     &
     & 0PF10.6,1P2E14.5)  )
!
  460 FORMAT (/,3X, 2I4,1P4E14.5,0PF10.6,1P2E14.5,I6,2(/12X,4E14.5,     &
     & 0PF10.6,1P2E14.5)  )
!
  470 FORMAT (/,' Panel: ',A,/,8x,'Number of Sections = ',i2,4x,        &
     & 'Number of Elements = ',I5,/,8x,'Total Surface Area = ', F12.3,  &
     & 4X,'Total Volume  = ',F12.3,/,3(10X,'**********'))
!
  480 FORMAT(/,'  INPUT SURFACE ELEMENT DATA',/,6X,1HN3X,1HM7X,1HX,     &
     & 3(13X,1HX),11X,2HNX,9X,5HXCENT,9X,4HAREA8X,1HL ,/1H ,5X, 4(13X,  &
     & 1HY),11X,2HNY9X,5HYCENT,7X,7HDELTA V,/1H ,5X,4(13X,1HZ),11X,2HNZ,&
     & 9X,5HZCENT ,7X,6HVOLUME,/ )
!
  490 RETURN
! 490 CONTINUE
      END Subroutine Geom


      SUBROUTINE IELE()
! ---------------------------------------------------------------------------
! PURPOSE - Convert one batch of type 3 input records to a network

IMPLICIT NONE
!
 !     COMMON /TAPE/TAPEIN,TAPEOT,TAPEA,TAPEB,TAPEC,TAPED,TAPEE,TAPEF,   &
 !    &             TAPEG,TAPEH,TAPEI,TAPEJ,TAPEK
 !     COMMON /EXEC/CASE,TITLE,PAGE,ERROR
!      COMMON /GFLAG/IOUT,ISTAT3,IORIN,COMPIN
      INTEGER  STAT,STATT,TYPE,SEQ
!      INTEGER  TAPEIN,TAPEOT,TAPEA,TAPEB,TAPEC,TAPED,TAPEE,TAPEF,       &
!     &         TAPEG,TAPEH,TAPEI,TAPEJ,TAPEK

  INTEGER:: cols,rows
  INTEGER:: i3max,in,irew,i3
  INTEGER:: is3
  INTEGER:: k,n
  INTEGER,DIMENSION(1000):: statusFlag
  REAL,DIMENSION(1000):: x,y,z
  REAL,DIMENSION(3):: xyz1,xyz2
!----------------------------------------------------------------------------
    
!
! the only item used on the next record is I3max
      READ(TAPEIN,'(2I2,2I1)') I3MAX,IN,IREW,I3   !  READ ELEMENT CONTROL CARD
!      IF (IN .EQ. 0) IN = TAPEIN
!      IF (IN.NE.5 .AND. IREW.EQ.1) REWIND IN
!
!      WRITE(TAPEOT,20) IN,IOUT,I3MAX
!   20 FORMAT (/,' TYPE 3 CARDS ARE BEING TRANSFERED FROM UNIT ',I2,     &
!     &  ' TO UNIT ',I2,5X,'I3MAX = ',I4)

      IF (I3MAX .EQ. 0) I3MAX = 1

  k=0
  is3=0
  DO
    READ(TAPEIN,40) XYZ1,STAT,XYZ2,STATT   !!!,CASE,SECT,TYPE,SEQ
   40 FORMAT (3F10.4,I1,3F10.4,I1,2X,I2,1A4,I2,4X,I4)
      IF (STAT .LE. 0) STAT = 0
      IF (STATT .LE. 0) STATT = 0
      k=k+1
      x(k)=-xyz1(1)
      y(k)=xyz1(2)
      z(k)=xyz1(3)
      statusFlag(k)=stat
      IF (stat==3) is3=is3+1
      IF (is3==i3max) EXIT
      k=k+1
      x(k)=-xyz2(1)
      y(k)=xyz2(2)
      z(k)=xyz2(3)
      statusFlag(k)=statt
      IF (statt==3) is3=is3+1
      IF (IS3==I3MAX) EXIT
  END DO
  n=k
  WRITE(*,*) n, 'grid points'
  
!...Now find the first index with statusFlag = 1. This should be one greater than the number
! of rows and it should divide n evenly. Points rows+1, 2*rows+1, 3*rows+1...
! should all have statusFlag = 1.
  DO k=1,n
    IF (statusFlag(k)==1) EXIT
  END DO
  rows=k-1
  cols=n/rows
  networks=networks+1
  
  WRITE(WGS,'(A,I4,A)') "'Network ", networks, " '"
  WRITE(WGS,'(3I4,A)') networks, cols, rows, ' 0   0 0 0   0 0 0   1 1 1   1'
  WRITE(WGS, '(3F15.6)') (x(k),y(k),z(k),k=1,n)  
  
  RETURN
END Subroutine Iele

FUNCTION Radd(BBI,AAI,THP) RESULT(r)
! -----------------------
! PURPOSE - Compute the radius of an ellipse
IMPLICIT NONE
  REAL,INTENT(IN):: bbi,aai
  REAL,INTENT(IN):: thp   ! angle, radians
  REAL:: r
!----------------------------------------------------------------------------
  r= SQRT(BBI*BBI*COS(THP)*COS(THP) + AAI*AAI*SIN(THP)*SIN(THP) )
  RETURN
END Function Radd   ! -------------------------------------------------------  

!+
SUBROUTINE FillArray(start,end, array, spacingCode)
! ---------------------------------------------------------------------------
! PURPOSE - fill an array from start to end. The intermediate points are
!    computed according to various spacing rules.

  REAL,INTENT(IN):: start,end
  REAL,INTENT(OUT),DIMENSION(:):: array
  INTEGER,INTENT(IN),OPTIONAL:: spacingCode  ! =2 full cosine
                                             ! =3 half cosine
                                             ! =4 half sine
                                             ! anything else = uniform spacing
                                             ! if omitted, uniform spacing

  INTEGER:: i,n
  REAL,PARAMETER:: PI=3.14159265, HALFPI=PI/2
  REAL,ALLOCATABLE,DIMENSION(:):: temp
!----------------------------------------------------------------------------
  n=SIZE(array)
  IF (n <= 0) RETURN

  array(n)=end
  array(1)=start
  IF (n <= 2) RETURN

  ALLOCATE(temp(n-2))
  temp(:)= (/ (REAL(i),i=1,n-2) /)   !  1. 2. 3. ...
  temp=temp*(1.0/REAL(n-1))

  IF (Present(spacingCode)) THEN
    SELECT CASE(spacingCode)
      CASE (2)
        temp=0.5*(1.0-COS(PI*temp))       ! full cosine, dense near both ends
      CASE (3)
        temp=1.0-COS(HALFPI*temp)             ! half cosine, dense near start
      CASE (4)
        temp=SIN(HALFPI*temp)                     ! half sine, dense near end
    END SELECT
  END IF

  array(2:n-1)=start + (end-start)*temp

  DEALLOCATE(temp)

  RETURN
END Subroutine FillArray   ! ------------------------------------------------


      SUBROUTINE ELLIP()
!     THIS SUBROUTINE PREPARES THE REQUIRED SURFACE ELEMENTS FOR
!     CIRCULAR OR ELLIPTICAL ARC SECTIONS.  EACH CROSS-SECTION
!     IS CONSIDERED SEPARATELY. DUMMY POINTS ARE COMPUTED SO THAT EACH
!     SECTION IS FORCED TO HAVE AN EVEN NUMBER OF POINTS AND SO THAT
!     POINTS IN A ROW ARE CORRECTLY MATCHED WITH POINTS IN AN ADJACENT
!     ROW WHEN THESE ROWS CONTAIN AN UNEQUAL NUMBER OF POINTS.
!
!  THE PARAMETER DISCON WHICH IS SPECIFIED BY THE PROGRAMMER IS VALUED
!     DEPENDING ON HOW THE POINTS ARE TO BE MATCHED
!          DISCON= 1  ALL THETAO AND THETAL ARE THE SAME. DELTHE MUST
!                     DIVIDE THE ANGULAR INCREMENT THETAL - THETAO
!                     EVENLY.
!                = 2  ALL THETAL ARE EQUAL BUT THETAO VARIES
!                = 3  ALL THETAO ARE EQUAL BUT THETAL VARIES
 !     COMMON /EXEC/CASE,TITLE,PAGE,ERROR
 !     COMMON /TAPE/TAPEIN,TAPEOT,TAPEA,TAPEB,TAPEC,TAPED,TAPEE,TAPEF,   &
 !    &             TAPEG,TAPEH,TAPEI,TAPEJ,TAPEK
      DIMENSION AX(100),THETOX(100),THETLX(100),DELTHX(100),  &
     &   NN(100),DELZX(100),DELYX(100),AA(100),BB(100)
  CHARACTER(LEN=4):: sect   
  INTEGER:: L  
  INTEGER:: last   
  REAL,PARAMETER:: PI = 3.14159265
  INTEGER:: rows
      INTEGER  STAT,STATT,STATD,STATC
!      INTEGER  TAPEIN,TAPEOT,TAPEA,TAPEB,TAPEC,TAPED,TAPEE,TAPEF,       &
!     &         TAPEG,TAPEH,TAPEI,TAPEJ,TAPEK
      INTEGER  STATA, STATB,PAGE,SEQ,TYPE,DISCON
  REAL,ALLOCATABLE,DIMENSION(:):: theta,y,z      

!----------------------------------------------------------------------------
  WRITE(*,*)'ELLIPTICAL GEOMETRY DATA IS BEING GENERATED'
!     SET COUNTERS
      TYPE = 3
      NREC = 0
!  READ IN TITLE CARD
   20 READ (TAPEIN,30) title(1:48),DISCON,IPRINT,CASE,SECT,ITYPE
   30 FORMAT (A,11X,2I1,3X,I2,1A4,I2)
      IF (ITYPE .NE. 4) THEN
        WRITE(*,*) 'Record type not 4 where required.'
        STOP
      END IF  
      LINE = 100
      SEQ  = 1
!  READ IN ALL DATA CARDS FOR THE SECTION
!      REWIND TAPEB
  DO i=1,SIZE(ax)
   40 READ (TAPEIN,50) ax(i),THETO,THETL,NN(I),aa(i),bb(i),DELZ,DELY,LAST,ITYPE
   50 FORMAT (F10.0,2F6.0,I3,2F10.0,2F7.0,I1,10X,I2)
      IF (ITYPE .NE. 5) THEN
        WRITE(*,*) 'Record type not 5 where required.'
        STOP
      END IF  
      DELTH = (THETL-THETO)/REAL(NN(I))
      THETOX(I) = (PI/180)*THETO
      THETLX(I) = (PI/180)*THETL
      DELTHX(I) = (PI/180)*DELTH
!      AA(I) = A
!      BB(I) = B
      DELYX(I) = DELY
      DELZX(I) = DELZ
      IF (LAST /= 0) THEN
        N = I
        EXIT
      END IF
  END DO
  
!... This program assumes that all NN's are the same. This is rows.
! But, each cross-section can have its own theta start and end  
  rows=nn(1)
  ALLOCATE(theta(rows),y(rows),z(rows))
  
!  110 GO TO (120,180,270), DISCON
!
!  120 M = M + 1
!
  WRITE(WGS,'(A,I4,A)') "'Network ", networks, " '"
  WRITE(WGS,'(3I4,A)') networks, n,rows, " 0   0 0 0   0 0 0   1 1 1   1"
  DO I=1,N   ! number of columns
    CALL FillArray(thetox(i),thetlx(i),theta)
    XA = AX(I)

    DO J=1,rows
      RAD = RADD(BB(I),AA(I),ABS(theta(j)-PI/2) )
      IF (RAD .NE. 0.0)  RAD = AA(I)*BB(I) / RAD
      Y(j) = RAD * SIN(THETA(j))
      Z(j) =-RAD * COS(THETA(j))
    END DO
    y=y+delyx(i)
    z=z+delzx(i)  
    WRITE(WGS,'(3F15.6)') (-ax(i),y(j),z(j),j=1,rows)
  END DO  
  DEALLOCATE(z,y,theta)
  RETURN
END Subroutine Ellip

      SUBROUTINE CUBIC()
! THIS SUBPROGRAM CALCULATES THE QUADRILATERAL DATA FOR A SURFACE GIVEN
! BY THE COONS MIT SURFACE FIT TECHNIQUE.
!
      DIMENSION XA(20),XB(20),YA(20),YB(20),ZA(20),ZB(20),XB1(4,20),YB1(&
     &4,20),ZB1(4,20),NPTS(4),D(4,9),TITLE(15),SECT(1)
      REAL L21,L31,L32,L1,M1,N1,L2,M2,N2,LN,MN,NN
      INTEGER STAT,STATT,TYPE,SEQ,CASE,PAGE,ERROR

      WRITE(*,*) 'PARAMETRIC CUBIC GEOMETRY DATA IS BEING GENERATED'

      TYPE=3
      SEQ=1
      LINE=100
      NREC=0
!
! *******READ IN BOUNDARY CURVE DATA
!     SET UP STARTING CONSTANTS
!
   20 CONTINUE
      N=-1
      L=0
   30 II=1
      ITRUE=0
      IFALSE=1
      READ (TAPEIN,40) (TITLE(K),K=1,12),NOU,NOW,LAST,ISOVR,IPRINT,CASE,&
     &     SECT,ITYPE
   40 FORMAT (12A4,1X,I3,1X,I3,3X,3I1,2X,I2,1A4,I2)
      IF (ITYPE .NE. 6) THEN
        WRITE(*,*) 'Type not 6 where required'
        STOP
      END IF  
!
!     READ IN BOUNDARY CURVE DATA
!
   50 CONTINUE
      READ (TAPEIN,60) X,Y,Z,ISTAT,XX,YY,ZZ,ISTATT,ITYPE
   60 FORMAT(3F10.4,I1,3F10.4,I1,8X,I2)
      IF (ITYPE .NE. 7) THEN
        WRITE(*,*) 'Type not 7 where required (1).'
        STOP
      END IF  
      IRFLAG=IFALSE
      GO TO 150
   70 IF(IRFLAG)80,90,80
   80 IRFLAG=ITRUE
      X=XX
      Y=YY
      Z=ZZ
      ISTAT=ISTATT
      GO TO 100
   90 IRFLAG=IFALSE
      READ (TAPEIN,60) X,Y,Z,ISTAT,XX,YY,ZZ,ISTATT,ITYPE
      IF (ITYPE .NE. 7) THEN
        WRITE(*,*) 'Type not 7 where required (2).'
        STOP
      END IF  
  100 IF(ISTAT)110,250,110
  110 IF(ISTAT-3)120,250,120
  120 IF(ISTAT-2)130,270,130
  130 IF(IAFLAG-1)140,270,140
  140 MC=M
  150 M=1
      IF(ISTAT-2)160,230,160
  160 IF(IBFLAG-1)170,200,170
  170 DO 180 J=1,MC
      XA(J)=XB(J)
      YA(J)=YB(J)
  180 ZA(J)=ZB(J)
  190 XB(1)=X
      YB(1)=Y
      ZB(1)=Z
      GO TO 70
  200 IF(IAFLAG)210,220,210
  210 IBFLAG=0
      GO TO 170
  220 IAFLAG=1
      GO TO 190
  230 IAFLAG=0
      IBFLAG=1
      N=N+1
  240 XA(M)=X
      YA(M)=Y
      ZA(M)=Z
      GO TO 70
  250 M=M+1
      IF(IAFLAG)260,240,260
  260 XB(M)=X
      YB(M)=Y
      ZB(M)=Z
      IF(ISTAT-3)70,270,70
  270 ML=MC
      MC=M
  280 N=N+1
!
!  SET UP BOUNDARY CURVE COORDINATE ARRAYS
  290 CONTINUE
      IF(II-1)300,300,320
  300 DO 310 I=1,ML
      XB1(II,I)=XA(I)
      YB1(II,I)=YA(I)
  310 ZB1(II,I)=ZA(I)
      NPTS(II)=ML
  320 II=II+1
      DO 330 I=1,MC
      XB1(II,I)=XB(I)
      YB1(II,I)=YB(I)
  330 ZB1(II,I)=ZB(I)
      NPTS(II)=MC
      IF(II-4)150,340,340
  340 CONTINUE
      IF(ISTAT-3)30,350,350
  350 CONTINUE
!
!
! *******CALCULATE BOUNDARY CURVE CONSTANTS
!  CALCULATE ARC LENGTH S ON BOUNDARY
      NB=1
  360 S=0.0
      K=NPTS(NB)-2
      DO 370 I=2,K
  370 S=S+SQRT((XB1(NB,I+1)-XB1(NB,I))**2+(YB1(NB,I+1)-YB1(NB,I))**2+(ZB&
     &1(NB,I+1)-ZB1(NB,I))**2)
!  CALCULATE TANGENT VECTORS AT THE START OF THE BOUNDARY
      IFLAG1=0
      J1=1
      J2=2
      J3=3
  380 X2X1=XB1(NB,J2)-XB1(NB,J1)
      X3X1=XB1(NB,J3)-XB1(NB,J1)
      X3X2=XB1(NB,J3)-XB1(NB,J2)
      Y2Y1=YB1(NB,J2)-YB1(NB,J1)
      Y3Y1=YB1(NB,J3)-YB1(NB,J1)
      Y3Y2=YB1(NB,J3)-YB1(NB,J2)
      Z2Z1=ZB1(NB,J2)-ZB1(NB,J1)
      Z3Z1=ZB1(NB,J3)-ZB1(NB,J1)
      Z3Z2=ZB1(NB,J3)-ZB1(NB,J2)
      L21=SQRT(X2X1*X2X1+Y2Y1*Y2Y1+Z2Z1*Z2Z1)
      L31=SQRT(X3X1*X3X1+Y3Y1*Y3Y1+Z3Z1*Z3Z1)
      L32=SQRT(X3X2*X3X2+Y3Y2*Y3Y2+Z3Z2*Z3Z2)
      L1=X3X1/L31
      M1=Y3Y1/L31
      N1=Z3Z1/L31
      L2=X2X1/L21
      M2=Y2Y1/L21
      N2=Z2Z1/L21
      LN=-(N1*(L1*N2-L2*N1)+M1*(L1*M2-L2*M1))
      MN=-(N1*(M1*N2-M2*N1)+L1*(L1*M2-L2*M1))
      NN=M1*(M1*N2-M2*N1)+L1*(L1*N2-L2*N1)
      COSEP1=(X2X1*X3X1+Y2Y1*Y3Y1+Z2Z1*Z3Z1)/(L21*L31)
      IF(COSEP1-0.999999)400,400,390
  390 EPS1=0.0
      GO TO 430
  400 IF(COSEP1+0.999999)410,420,420
  410 EPS1=0.0
      GO TO 430
  420 EPS1=ACOS(COSEP1)
  430 COSEP2=(X3X2*X3X1+Y3Y2*Y3Y1+Z3Z2*Z3Z1)/(L32*L31)
      IF(COSEP2-0.999999)450,450,440
  440 EPS2=0.0
      GO TO 480
  450 IF(COSEP2+0.999999)460,470,470
  460 EPS2=0.0
      GO TO 480
  470 EPS2=ACOS(COSEP2)
  480 DELTA=EPS1+EPS2
      TX=L1*COS(DELTA)+LN*SIN(DELTA)
      TY=M1*COS(DELTA)+MN*SIN(DELTA)
      TZ=N1*COS(DELTA)+MN*SIN(DELTA)
!
!  CALCULATE END POINT DERIVATIVES
      IF(IFLAG1)490,490,500
  490 X1VOO=TX*S
      Y1VOO=TY*S
      Z1VOO=TZ*S
      J1=NPTS(NB)-2
      J2=NPTS(NB)-1
      J3=NPTS(NB)
      IFLAG1=1
      GO TO 380
  500 X1VO1=TX*S
      Y1VO1=TY*S
      Z1VO1=TZ*S
!
!
! *********CALCULATE CONSTANTS FOR BOUNDARY CURVE
      D(NB,1)=2.0*(XB1(NB,2)-XB1(NB,J2))+X1VOO+X1VO1
      D(NB,2)=3.0*(XB1(NB,J2)-XB1(NB,2))-2.0*X1VOO-X1VO1
      D(NB,3)=X1VOO
      D(NB,4)=2.0*(YB1(NB,2)-YB1(NB,J2))+Y1VOO+Y1VO1
      D(NB,5)=3.0*(YB1(NB,J2)-YB1(NB,2))-2.0*Y1VOO-Y1VO1
      D(NB,6)=Y1VOO
      D(NB,7)=2.0*(ZB1(NB,2)-ZB1(NB,J2))+Z1VOO+Z1VO1
      D(NB,8)=3.0*(ZB1(NB,J2)-ZB1(NB,2))-2.0*Z1VOO-Z1VO1
      D(NB,9)=Z1VOO
      NB=NB+1
      IF(NB-4)360,360,510
!
! *********CALCULATE PATCH DATA
  510 NOW=NOW/2*2+1
      DELU=1.0/FLOAT(NOU)
      DELW=1.0/FLOAT(NOW)
      NOU=NOU+1
      NOW=NOW+1
      STATT=0
      U=0.0
!
  WRITE(WGS,'(A,I4,A)' ) "'network ", networks, " '"
  WRITE(WGS,'(3I4,A)' ) networks, nou,now, " 0   0 0 0   0 0 0   1 1 1   0"

      DO 560 I=1,NOU
      STAT=1
      INU=0
      W=0.0
!
      DO 550 K=1,NOW
!
!
      W3=W**3
      W2=W**2
      U3=U**3
      U2=U**2
!  CALCULATE BLENDING FUNCTIONS
      F1U=3.0*U2-2.0*U3
      FOU=1.0-3.0*U2+2.0*U3
      F1W=3.0*W2-2.0*W3
      FOW=1.0-3.0*W2+2.0*W3
!
!  CALCULATE POINTS ON BOUNDARY CURVES
      XOW=D(1,1)*W3+D(1,2)*W2+D(1,3)*W+XB1(1,2)
      YOW=D(1,4)*W3+D(1,5)*W2+D(1,6)*W+YB1(1,2)
      ZOW=D(1,7)*W3+D(1,8)*W2+D(1,9)*W+ZB1(1,2)
      X1W=D(2,1)*W3+D(2,2)*W2+D(2,3)*W+XB1(2,2)
      Y1W=D(2,4)*W3+D(2,5)*W2+D(2,6)*W+YB1(2,2)
      Z1W=D(2,7)*W3+D(2,8)*W2+D(2,9)*W+ZB1(2,2)
      XUO=D(3,1)*U3+D(3,2)*U2+D(3,3)*U+XB1(3,2)
      YUO=D(3,4)*U3+D(3,5)*U2+D(3,6)*U+YB1(3,2)
      ZUO=D(3,7)*U3+D(3,8)*U2+D(3,9)*U+ZB1(3,2)
      XU1=D(4,1)*U3+D(4,2)*U2+D(4,3)*U+XB1(4,2)
      YU1=D(4,4)*U3+D(4,5)*U2+D(4,6)*U+YB1(4,2)
      ZU1=D(4,7)*U3+D(4,8)*U2+D(4,9)*U+ZB1(4,2)
      NPT1=NPTS(1)-1
      NPT2=NPTS(2)-1
!
!  CALCULATE POSITION OF A POINT ON THE SURFACE
      XS = XOW*FOU+X1W*F1U+XUO*FOW+XU1*F1W-XB1(1,2)*FOU*FOW-XB1(1,NPT1)*&
     &  FOU*F1W-XB1(2,2)*F1U*FOW-XB1(2,NPT2)*F1U*F1W
      YS = YOW*FOU+Y1W*F1U+YUO*FOW+YU1*F1W-YB1(1,2)*FOU*FOW-YB1(1,NPT1)*&
     &  FOU*F1W-YB1(2,2)*F1U*FOW-YB1(2,NPT2)*F1U*F1W
      ZS = ZOW*FOU+Z1W*F1U+ZUO*FOW+ZU1*F1W-ZB1(1,2)*FOU*FOW-ZB1(1,NPT1)*&
     &  FOU*F1W-ZB1(2,2)*F1U*FOW-ZB1(2,NPT2)*F1U*F1W
     
     WRITE(WGS,'(3F15.6)') -xs,ys,zs
     
      IF(INU-1)520,530,530
  520 XXS = XS
      YYS = YS
      ZZS = ZS
      INU=1
      GO TO 540
  530 IF (I.EQ.NOU .AND. K.EQ.NOW .AND. LAST.EQ.1) STATT = 3
      IF (I.EQ.NOU .AND. K.EQ.NOW .AND. LAST.EQ.3) STATT = 3
      IF (I.EQ.NOU .AND. K.EQ.NOW .AND. LAST.EQ.4) STATT = 3
      IF (I.EQ.1 .AND. K.EQ.2 .AND. ISOVR.EQ.0) STAT = 2
!      CALL PUNCH (XXS,YYS,ZZS,STAT,XS,YS,ZS,STATT,SECT,TYPE,LINE,SEQ,   &
!     &    LAST,IPRINT,NREC)
      INU=0
      STAT=0
  540 W=W+DELW
  550 END DO
      U=U+DELU
  560 END DO
  
  IF (LAST==2) GO TO 20
  
!  570 GO TO 590
!  580 ERROR = 1
!  590 RETURN
! 590 CONTINUE

  RETURN
END Subroutine Cubic

      SUBROUTINE PUNCH (X1,Y1,Z1,NSTAT1,X2,Y2,Z2,NSTAT3,SECT,TYPE,      &
     &     LINE,SEQ,LAST,IPRINT,NREC)
!
!  THIS SUBROUTINE PREPARES VEHICLE GEOMETRY DATA IN THE PROPER FORM
!   FOR USE BY SDATA ROUTINE
!
      COMMON /EXEC/CASE,TITLE,PAGE,ERROR
      COMMON /GFLAG/IOUT,ISTAT3,IORIN,COMPIN
      COMMON /TAPE/TAPEIN,TAPEOT,TAPEA,TAPEB,TAPEC,TAPED,TAPEE,TAPEF,   &
     &             TAPEG,TAPEH,TAPEI,TAPEJ,TAPEK
      DIMENSION  TITLE(15),SECT(1)
!
      INTEGER  PAGE,SEQ,ERROR,CASE
      INTEGER  TAPEIN,TAPEOT,TAPEA,TAPEB,TAPEC,TAPED,TAPEE,TAPEF,       &
     &         TAPEG,TAPEH,TAPEI,TAPEJ,TAPEK
      INTEGER TYPE
!
      NSTAT2 = NSTAT3
!
!  CHECK IF THIS IS THE LAST POINT OF THE ENTIRE VEHICLE
      IF (NSTAT3.EQ.3 .AND. LAST.EQ.1) NSTAT2 = 0
      IF (IPRINT .EQ. 0) GO TO 40
!
      IF  (LINE .LT. 50 )  GO TO 20
!
!  WRITE PAGE HEADER FOR STANDARD OUTPUT TAPE
      WRITE (TAPEOT,10) CASE,(TITLE(L),L=1,12),PAGE
   10 FORMAT (1H1,5X,36HANALYTICALLY GENERATED ELEMENT DATA ,/          &
     &   1H0,6H  CASE,3X,I2,17X,12A4,17X,5HPAGE ,I4,/                   &
     &  1H0,5X,1HX,9X,1HY9X,1HZ4X,1HS5X,1HX,9X,1HY8X,1HZ5X,1HS10H CASE S&
     &ECT, 6X,3HSEQ )
!
      PAGE = PAGE + 1
      LINE = 5
!
!  WRITE GEOMETRY CARDS ON STANDARD OUTPUT TAPE
   20 WRITE (TAPEOT,30) X1,Y1,Z1,NSTAT1,X2,Y2,Z2,NSTAT2,CASE,SECT,      &
     &      TYPE,SEQ
   30 FORMAT (1H0,3F10.4,I1,3F10.4,I1,2X,I2,A4,1X,I1,4HAERO,I4 )
!
      LINE = LINE + 2
!
!  WRITE GEOMETRY DATA ON GEOMETRY TAPE
   40 WRITE (IOUT,50)   X1,Y1,Z1,NSTAT1,X2,Y2,Z2,NSTAT2,CASE,SECT,      &
     &      TYPE,SEQ
   50 FORMAT(3F10.4,I1,3F10.4,I1,2X,I2,A4,1X,I1,4HAERO,I4 )
!
      NREC = NREC + 1
      SEQ = SEQ + 1
!
      RETURN
      END Subroutine Punch

END Module HypersonicArbitraryBodyToLangleyWireframeGeometryProcedures

!+
PROGRAM HypersonicArbitraryBodyToLangleyWireframeGeometry
! ----------------------------------------------------------------------
! PURPOSE - 


USE HypersonicArbitraryBodyToLangleyWireframeGeometryProcedures
IMPLICIT NONE

      INTEGER:: errCode
      CHARACTER(LEN=80):: fileName
      INTEGER,PARAMETER:: TAPEIN=1, WGS=2, DBG=3, TAPEOT=4
      INTEGER:: networks
      CHARACTER(LEN=60):: title

!      COMMON /EXEC/CASE,TITLE,PAGE,ERROR


      DIMENSION  IPG(20)
      INTEGER::      PAGE,ERROR
!----------------------------------------------------------------------------
  networks=0
  
      DO
        WRITE(*,*) 'Enter input file name:'
        READ(*,'(A)') fileName
        IF (Len_Trim(fileName)==0) STOP
        OPEN(UNIT=TAPEIN,FILE=fileName,STATUS='OLD',                    &
     &   IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
        IF (errCode==0) EXIT
        OPEN(UNIT=TAPEIN,FILE=Trim(fileName)//'.mk4',STATUS='OLD',      &
     &   IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
        IF (errCode==0) EXIT
        OPEN(UNIT=TAPEIN,FILE=Trim(fileName)//'.mk5',STATUS='OLD',      &
     &   IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
        IF (errCode==0) EXIT
        WRITE(*,*) 'Unable to find this file. Try again.'
      END DO
      INQUIRE(UNIT=TAPEIN, NAME=fileName)
      WRITE(*,*) 'Reading from ',Trim(fileName)
      open(UNIT=WGS,FILE='shabp.wgs',STATUS='REPLACE',ACTION='WRITE')
      OPEN(UNIT=DBG,FILE='hab2wgs.dbg',STATUS='REPLACE',ACTION='WRITE')
      open(unit=tapeot,file='SHABP.OUT',status='REPLACE',ACTION='WRITE')

!
!**   INITIALIZE PARAMETERS
!
      ERROR  =  0
                                                 
  READ(TAPEIN,'(2I1,1X,A)') IEROR,INMONT,title(1:60) ! Executive Flag Card
  WRITE(WGS,*) "'"//title//"'"
 !      IF (INMONT .EQ. 0) CALL MONITR(iecho)
                                   
      READ (TAPEIN,'(20I1)') IPG  ! System Control Card

! *** Determine the number of program options entered
      I = 0
      NPROG = 0
   35 I = I+1
      IF(IPG(I) .GT. 0 .AND. IPG(I) .LT. 4)THEN
        NPROG = NPROG + 1
        GO TO 35
      else if(ipg(i).eq.9) then
        continue
      ELSE IF(IPG(I) .GT. 3)THEN
        WRITE(TAPEOT,155)IPG(I)
  155   FORMAT(' *** ERROR:',/,                                         &
     &  ' UNKNOWN PROGRAM OPTION ',I1,'.  PROGRAM TERMINATING')
        STOP
      ENDIF
      IF(NPROG .EQ. 0)THEN
        WRITE (TAPEOT,230)
  230   FORMAT (' *** ERROR:',/,' **** FIRST PHASE OPTION IS ZERO ****')
        STOP
      ENDIF
!
! *** Write program options selected by user to output file
!
      WRITE (TAPEOT,30)
   30 FORMAT (////,7x,'*** SUPERSONIC-HYPERSONIC ARBITRARY-BODY',       &
     & ' AERODYNAMIC SYSTEM ***',///,                                   &
     & 10X,'PROGRAM OPTIONS ARE IN THE FOLLOWING ORDER....   ',/)
!
      DO 80 I=1,NPROG
        IF (IPG(I) .EQ. 1) WRITE (TAPEOT,40) I
   40     FORMAT (/,15X,I2,' GEOMETRY PROGRAM     (OPTION 1)',/)
        IF (IPG(I) .EQ. 2) WRITE (TAPEOT,50) I
   50     FORMAT (/,15X,I2,' AERODYNAMIC PROGRAM  (OPTION 2)',/)
        IF (IPG(I) .EQ. 3) WRITE (TAPEOT,70) I
   70     FORMAT (/,15X,I2,' AUXILIARY PROGRAMS   (OPTION 3)',/)
   80 END DO
!
! *** Call program options selected by user
!
      DO 500 I=1,NPROG
        IPROG = IPG(I)
        IF(IPROG .EQ. 1)THEN
          CALL GEOM()
        ENDIF
        IF (ERROR .NE. 0) THEN
          WRITE(TAPEOT,170) I, ERROR
  170     FORMAT (/,' *** ERROR:',/,                                    &
     &    ' ***** AN ERROR HAS OCCURRED IN OPTION ',I2,/,               & 
     &    ' AND WAS DETECTED AFTER RETURN TO THE MAIN EXECUTIVE',       &
     &    ' ROUTINE., ERROR = ',I2)
          WRITE (TAPEOT,210)
  210     FORMAT (//,' ******** FATAL ERROR *** PROGRAM STOPPED')
          STOP
        ENDIF
  500 END DO

  WRITE(*,*) '*** NORMAL TERMINATION ***'
  STOP

END PROGRAM HypersonicArbitraryBodyToLangleyWireframeGeometry

