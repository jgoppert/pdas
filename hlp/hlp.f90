MODULE HiddenLineProcedures
! ------------------------------------------------------------------------------

CONTAINS

!+
      SUBROUTINE Plot(x,y,ipen) 
!   --------------------------------------------------------------------
!     PURPOSE - Write appropriate information to the temporary file 

      IMPLICIT NONE 
!***********************************************************************
!     A R G U M E N T S                                                *
!***********************************************************************
      REAL,INTENT(IN):: x         ! x-coordinate
      REAL,INTENT(IN):: y         ! y-coordinate
      INTEGER,INTENT(IN):: ipen   ! =2 for draw; =3 for move
!***********************************************************************
!     C O N S T A N T S                                                *
!***********************************************************************
      INTEGER PLT 
      PARAMETER (PLT=2)  ! note that this MUST agree with the main       
!-----------------------------------------------------------------------
      IF (ipen .EQ. 3) THEN 
        WRITE(PLT,'(/2ES15.4)')  x,y   ! move to this point
      ELSE 
        WRITE(PLT, '(2ES15.4)' ) x,y   ! draw to this point
      ENDIF 
      RETURN
      END Subroutine Plot   ! --------------------------------------------------

!+
! The following programs are part of COSMIC program ARC-12721           
!    SILHOUETTE - Hidden Line Computer Code with                        
!        Generalized Silhouette Solution                                
!          by David R. Hedgley, Jr.                                     
!           Ames Research Center                                        
!            Dryden Flight Research Facility                            
!             Edwards, California                                       
!                                                                       
!                                                                       
!                                                                       
! MODIFICATIONS BY PUBLIC DOMAIN AERONAUTICAL SOFTWARE:                 
!                                                                       
!   1. The common block /DRH/ has been replaced with the statement      
!      INCLUDE 'drh.cmn'                                                
!      at all occurences in the program. This allows the size           
!      to be modified with one change, rather than trying to find       
!      them all.                                                        
!   2. /GO/ is currently set for 1000 elements.                         
!   3. In subroutine LIN, the dimensions of array ICCT have been        
!      changed from (1:200) to (0:200). When testing the program        
!      I frequently encountered messages that ICCT(IV1) was             
!      being accessed with IV1=0. There may still be a problem,         
!      but the error messages are gone.                                 
!   4. The routines VSRTR and VSRT1 that are included with the          
!      COSMIC distribution have been replaced with equivalent           
!      routines that have no copyright restrictions.                    
                                                                        
                                                                        
                                                                        
      SUBROUTINE STATUM(OJ,TMJ,XXX,TGM,RV,RVI,TGI,ZM,NNO,II,            &
     & H,IM,JXT,ZJ,NC,ZMI,CCC,LZ)                                        
      DIMENSION X(50),Y(50) 
      DIMENSION CCC(*),XXX(*),ZMI(*),TGM(*),RV(*),RVI(*),TGI(*)
      DIMENSION NNO(*),H(*)                                            
      dimension zm(*)   !!! ???           
      COMMON/GO3/L0,L1,L00,L01,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13 
      include 'drh.cmn' 
!!!      COMMON/DRH/ISILH(1)                                            
      DATA GGK,F,EI/.005,.045,.2/ 
      OJ1=OJ 
      TMJ1=TMJ 
      IM=0 
      J=1 
      T=2.*F/J 
      Y(1)=TMJ+F 
      X(1)=OJ-F 
      DO 97 I=1,J 
      X(I+1)=X(I)+T 
      Y(I+1)=TMJ+F 
   97 END DO 
      DO 98 I=1,J 
      X(J+1+I)=OJ+F 
      Y(J+1+I)=Y(J+1+I-1)-T 
   98 END DO 
      N=2*J+1 
      DO 99 I=1,J 
      Y(N+I)=Y(N) 
      X(N+I)=X(N+I-1)-T 
   99 END DO 
      N=3*J+1 
      DO 100 I=1,J 
      X(N+I)=X(N) 
      Y(N+I)=Y(N+I-1)-T 
  100 END DO 
      KS=4*J 
    9 CONTINUE 
      DO 80 L=1,KS 
      OJ=X(L) 
      TMJ=Y(L) 
   12 CONTINUE 
      D=EI*OJ-TMJ 
      DO 65 JO=1,II 
      JG=NNO(L4+JO) 
      IF(ISILH(JXT).NE.ISILH(JG))GO TO 60 
      IF((TMJ.GE.RV(L7+JG)).OR.(TMJ.LE.RVI(L8+JG)))GO TO 60 
      IF((OJ.GE.TGI(L6+JG)).OR.(OJ.LE.TGM(L5+JG)))GO TO 60 
      JS=L13+(JG-1)*LZ 
      JT=L12+(JG-1)*5 
      JN=0 
      IF(JG.EQ.JXT)JN=H(6) 
      NS=XXX(5+JT)-JN 
      IB=NS*5 
      DO 20 J=1,IB,5 
      JJ=J+JS 
      IF(CCC(JJ).NE.0)GO TO 15 
      S=TMJ-CCC(JJ+3) 
      S1=TMJ-CCC(JJ+4) 
      DY=(-CCC(JJ+2)/CCC(JJ+1))-OJ 
      GO TO 16 
   15 CONTINUE 
      S=OJ-CCC(JJ+3) 
      S1=OJ-CCC(JJ+4) 
      DY=(-CCC(JJ+2)-CCC(JJ+1)*OJ)-TMJ 
   16 CONTINUE 
      IF((ABS(DY).GE.GGK).OR.(S*S1.GT.0.))GO TO 20 
      GO TO 80 
   20 END DO 

      I=0 
      DO 40 J=1,IB,5 
      JJ=J+JS 
      R=EI*CCC(JJ)+CCC(JJ+1) 
      IF(R.EQ.0.)GO TO 40 
      T=(-CCC(JJ+2)+CCC(JJ)*D)/R 
      IF(T.LT.OJ)GO TO 40 
   34 CONTINUE 
      IF(CCC(JJ).NE.0.)GO TO 30 
      T=EI*T-D 
   30 CONTINUE 
      IF((T-CCC(JJ+3))*(T-CCC(JJ+4)).GT.0.)GO TO 40 
      I=I+1 
   40 END DO 

      IF(I-(I/2)*2.EQ.0)GO TO 60 
      IM=1 
      GO TO 80 
   60 CONTINUE 
   65 END DO 

      IM=0 
      GO TO 90 
  800 CONTINUE 
   80 END DO 
      IM=1 
   90 CONTINUE 
      TMJ=TMJ1 
      OJ=OJ1 
      RETURN 
      END Subroutine Statum   ! --------------------------------------------------

      SUBROUTINE SKETCH(X,Y,Z,NP,NC) 



!     THIS SUBROUTINE SETS UP PEN MOTION INDICATORS.                    


      DIMENSION X(*),Y(*),Z(*) 
      DIMENSION X1(160),Y1(160),Z1(160) 
      COMMON/SCALAR/SCX,YAW,ROL,PIT,LZ,VP,JJJ,ICORE 

      x1(:)=0.0   ! these 3 statements added by RLC, 24 Dec 2004
      y1(:)=0.0
      z1(:)=0.0

      L=NP 
      LI=NP 
      IF(L.LE.2)GO TO 50 
      LX=1 
      NPX=NP 
    1 NPX=NPX-1 
      I=LX 
      DO 8 M=I,NPX 
      RX=0 
      A=X(M+1)-X(M) 
      B=Y(M+1)-Y(M) 
      C=Z(M+1)-Z(M) 
      IF(A.NE.0.)GO TO 8 
      IF(B.NE.0.)GO TO 8 
      IF(C.NE.0.)GO TO 8 
      IX=M 
      IX1=NPX 
      DO 4 MX=IX,IX1 
      X(MX)=X(MX+1) 
      Y(MX)=Y(MX+1) 
      Z(MX)=Z(MX+1) 
    4 END DO 
      RX=1 
      LX=M 
      IF(LX.EQ.NPX)GO TO 10 
      GO TO 1 
    8 END DO 
   10 CONTINUE 
      IF(RX.EQ.1.)NPX=NPX-1 
      NP=NPX+1 
      LI=NP 
      IF(NP.LE.2)GO TO 50 
      IX=0 
      M1=0 
      M=1 
      IS=NP-1 
   20 CONTINUE 
      M=M+IX 
      M1=M1+IX+1 
      IF(M-1.EQ.LI)GO TO 70 


!     SEARCH FOR MATCHING COORDINATES.                                  


      DO 40 J=M,IS 
      T=X(J+1)-X(M) 
      U=Z(J+1)-Z(M) 
      W=Y(J+1)-Y(M) 
      IF(T.NE.0.)GO TO 40 
      IF(W.NE.0.)GO TO 40 
      IF(U.NE.0.)GO TO 40 
      NP=NP+1 

!     MATCH FOUND.....STORE COORDINATES AND SET SWITCH TO LIFT PEN      
!     AND/OR END SET.                                                   


      IX=J+2-M 
      IX1=J-IS+1 
      I2=M1-1 
      I3=M-1 
      DO 30 IK=1,IX 
      X1(I2+IK)=X(I3+IK) 
      Y1(I2+IK)=Y(I3+IK) 
      Z1(I2+IK)=Z(I3+IK) 
   30 END DO 
      Z1(M1+IX)=-ISIGN(1,IX1)*9999. 
      GO TO 20 
   40 END DO 
   50 CONTINUE 
      DO 60 J=1,LI 
      X1(J)=X(J) 
      Y1(J)=Y(J) 
      Z1(J)=Z(J) 
   60 END DO 
      NP=NP+1 
      Z1(NP)=-9999. 
   70 CONTINUE 
      CALL LIN(X1,Y1,Z1,NP,NC)   ! the only call to Lin
      NP=L 


!     RESET VALUE FOR MAXIMUM NUMBER OF EDGES IF ARGUMENT IS            
!     COMPLETED.                                                        


      IF(VP.GT.0.)LZ=LZ/5 
      RETURN 
      END Subroutine Sketch   ! ----------------------------------

!+
      SUBROUTINE COEF(X,Y,Z,XXX,JXX,NC,NS,CCC,LZ) 
! -----------------------------------------------------------------------
! PURPOSE - This subroutine determines equation of lines and planes.

      DIMENSION CCC(*),XXX(*),X(*),Y(*),Z(*)    ! were 1s. RLC
      DIMENSION ZZ(30),ZZZ(30) 
      DIMENSION IB(5) 
      DOUBLE PRECISION T,T1,E,F,A1,B1,C1,A2,B2,C2,COE(8) 
      COMMON/GO3/L0,L1,L00,L01,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13 
!!!      COMMON/DRH/ISILH(1)                                            
      INCLUDE 'drh.cmn' 
      DATA EPSS/.00001/ 
      J1=JXX-1 
      LE=0 
      JA=L13+J1*LZ 
      JF=L12+J1*5 
      I=0 
      J=1 
   10 CONTINUE 


!     SEARCH FOR MATCHING COORDINATES.                                  


      I=I+1 
      JN=I+1 
      IF(X(JN)-X(I).NE.0.)GO TO 20 
      IF(Y(JN)-Y(I).NE.0.)GO TO 20 
      IF(Z(JN)-Z(I).NE.0.)GO TO 20 


!     MATCH FOUND.....PROCEED IF LIST IS NOT EXHAUSTED.                 


      I=I+2 
   20 CONTINUE 
      IF(I.GT.NS)GO TO 70 


!     DETERMINE EQUATION OF LINE-SEGMENTS.                              


      T=X(I+1)-X(I) 
      T1=Y(I+1)-Y(I) 
      IF((T.EQ.0.).AND.(T1.EQ.0.))GO TO 10 
      IF(T.NE.0.)GO TO 30 
   29 CONTINUE 
      CCC(J+JA)=0 
      CCC(J+1+JA)=1 
      CCC(J+2+JA)=-X(I) 
      GO TO 40 
   30 CONTINUE 
      CCC(J+JA)=1 
      E=T1/T 
      F=E*X(I)-Y(I) 
      IF(DABS(E).GE.10000.)GO TO 29 
      IF(DABS(E).LT.100.)GO TO 220 
      YO=E*X(I+1)-F 
      G=ABS(Y(I+1)-YO) 
      IF(G.GT.EPSS)GO TO 29 
  220 CONTINUE 
      CCC(J+1+JA)=-E 
      CCC(J+2+JA)=F 
   40 CONTINUE 
      IF(CCC(J+JA).NE.0.)GO TO 50 
      CCC(J+3+JA)=Y(I) 
      CCC(J+4+JA)=Y(I+1) 
      GO TO 60 
   50 CONTINUE 
      CCC(J+3+JA)=X(I) 
      CCC(J+4+JA)=X(I+1) 
   60 CONTINUE 
      J=J+5 
      LE=LE+1 
      ZZ(LE)=Z(I) 
      ZZZ(LE)=Z(I+1) 
      IF(LE.GT.3)GO TO 10 
      IB(LE)=I 
      GO TO 10 
   70 CONTINUE 


!     DETERMINE EQUATION OF PLANE.                                      

      J=(J-1)/5 
      XXX(JF+5)=J 
      IF(NS.LE.3)GO TO 120 
      A1=X(3)-X(1) 
      B1=Y(3)-Y(1) 
      C1=Z(3)-Z(1) 
      A2=X(2)-X(1) 
      B2=Y(2)-Y(1) 
      C2=Z(2)-Z(1) 
      COE(1)=B1*C2-B2*C1 
      COE(2)=C1*A2-C2*A1 
      COE(3)=A1*B2-A2*B1 
      COE(4)=COE(1)*X(1)+COE(2)*Y(1)+COE(3)*Z(1) 
      COE(4)=-COE(4) 
      DO 110 J=1,4 
  110 XXX(JF+J)=COE(J) 
      IF(COE(3).NE.0.)GO TO 140 
      J=1 
      DO 25 K=1,LE 
      CCC(JA+J)=ZZ(K) 
      CCC(JA+J+1)=ZZZ(K) 
      J=J+5 
   25 END DO 
      IF(COE(1).NE.0.)I=1 
      IF(COE(2).NE.0.)I=2 
      P=COE(I) 
      DO 26 K=1,4 
   26 XXX(JF+K)=XXX(JF+K)/P 
      GO TO 140 
  120 CONTINUE 
      XXX(JF+5)=1 
      DO 130 IX=1,2 
  130 XXX(JF+IX)=Z(IX) 
      XXX(JF+3)=0 
  140 CONTINUE 
      RETURN 
      END Subroutine Coef   ! --------------------------------------------------

!+
      SUBROUTINE LIN(X,Y,Z,NP,NC) 
! ------------------------------------------------------------------------------
! PURPOSE - This subroutine is the executive. 

      DIMENSION X(*),Y(*),Z(*)
      INCLUDE 'drh.cmn' 
! the following dimension statements do not reflect the actual usage.
! they are all equivalenced to WORK and are a rather dangerous form of
! programming where all of the variables use the same storage. The 
! programmer must be sure that a quantity is not needed later if it is
! commanded to be overwritten. All were originally written as being
! of dimension 1, so there was never any possibility of checking if the
! array usage exceeded that which was expected.
      INTEGER XIND(16*MAXPANELS)
      DIMENSION XXX(30*MAXPANELS),CCC(26*MAXPANELS) 
      DIMENSION TGM(MAXPANELS),IN(14*MAXPANELS),ZM(5*MAXPANELS)
      DIMENSION ZMI(6*MAXPANELS)
      DIMENSION TGMT(13*MAXPANELS),TGI(2*MAXPANELS),NNO(7*MAXPANELS)
      DIMENSION RV(3*MAXPANELS),RVI(4*MAXPANELS),NOCT(12*MAXPANELS)                        
      DIMENSION ZMIN(MAXPANELS),YMIN(10*MAXPANELS)
      DIMENSION IN1(11*MAXPANELS),IN2(MAXPANELS)
      DIMENSION COORD(69*MAXPANELS) 
      DIMENSION SNDT(7*MAXPANELS) 
      DIMENSION IND(15*MAXPANELS) 
      DIMENSION NEH(69*MAXPANELS),KEEP(20*MAXPANELS)         ! changed by RLC
 
      DIMENSION RRX(20) 
      DIMENSION RCT(200) 
      DIMENSION NGX(15),IADR(200) 
      DIMENSION X21(500),Y21(500),Z21(500),IIA(500) 
      DIMENSION XI(1000),YI(1000),ZI(1000),DI(1000) 
      DIMENSION U(6),V(6),W(6),X1(10),Y1(10),Z1(10),H(15) 
      DIMENSION ICOUNT(150) 
      DIMENSION XE(150),YE(150),XU(150),YU(150) 
      DIMENSION I2(2),I3(2) 
      DIMENSION REX(20) 
      DIMENSION IBEG(200),IEND(200),ICT(200),ICCT(0:200) 
      DIMENSION XB(70),YB(70),ZB(70) 
      DIMENSION YSUM(3),SUM1(3,3) 
!      REAL work(69*25000) 
!      COMMON/GO/WORK 
      COMMON/GO3/L0,L1,L00,L01,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13 
      COMMON/SCALAR/SCX,YAW,ROLL,PIT,LZ, VP,JJJ,ICORE 
      COMMON/HEDG/JAT,ME,JT,D4,D2,D1,D3,NS 
!!!      COMMON/DRH/ISILH(1) ! added by RLC                             
!!!      INCLUDE 'go.cmn'       ! added by RLC                          
      COMMON/ISIL/ISIL 
      COMMON/INDX/JT1,JO 
      EQUIVALENCE(WORK(1),XXX(1),CCC(1),IN(1),TGMT(1)) 
      EQUIVALENCE(WORK(1),ZM(1),ZMI(1),NNO(1))
      EQUIVALENCE(WORK(1),TGM(1),TGI(1),RV(1),RVI(1),NOCT(1))                               
      EQUIVALENCE(WORK(1),IN1(1),IN2(1),YMIN(1),ZMIN(1)) 
      EQUIVALENCE(WORK(1),COORD(1)) 
      EQUIVALENCE(WORK(1),SNDT(1)) 
      EQUIVALENCE(WORK(1),IND(1),XIND(1)) 
      EQUIVALENCE(WORK(1),KEEP(1),NEH(1)) 

      IF(VP.LT.0.)GO TO 20 
      HXX=.005 
      AVA=0 
      ISIL=0 
      HX1=.001 
      LC=10**6 
      IXXX=0 
      IF(SCX .LT. 0.)IXXX=1 
      SCX=ABS(SCX) 

!     INITIALIZE VARIABLES.                                             

      LZ=LZ*5 
      AXMIN=10.**6 
      AXMAX=-10.**6 
      SW1=0 
      VTX=VP+10 
      SW=0 
      IDAV=0 
      I1=0   ! added by RLC  20Dec2004

!     CALCULATE MAXIMUM ALLOWABLE ELEMENTS.                             

      IABC=ICORE/(25+LZ+4*JJJ) 
      ISAVE=NC 
      NC=IABC 
      L5=0 
      L6=NC 
      L7=2*NC 
      L8=3*NC 
      L2=4*NC 
      L3=5*NC 
      L4=6*NC 
      L00=7*NC 
      L01=8*NC 
      L1=9*NC 
      L0=10*NC 
      L9=11*NC 
      L10=12*NC 
      L11=13*NC 
      L15=14*NC 
      L16=15*NC 
      L17=16*NC 
      L18=19*NC 
      L12=20*NC 
      L13=25*NC 
      L14=L13+LZ*NC 

      DO 10 J=1,NC 
      RVI(L8+J)=10**6 
      TGM(L5+J)=10**6 
      RV(L7+J)=-RVI(L8+J) 
      TGI(L6+J)=-TGM(L5+J) 
      NOCT(L9+J)=0 
      ZM(L2+J)=RV(L7+J) 
      ZMI(L3+J)=RVI(L8+J) 
      XIND(L16+J)=0 
      IND(L15+J)=J 
      KEEP(L18+J)=0 
   10 END DO 

      NC=ISAVE 
      IK=0 
      IKT=0 
      PI=3.1416/180. 
      IJBB=0 
      VP=-VP 
!                                                                       
!     STORE EULERIAN ANGLES.                                            
!                                                                       
      XX=YAW*PI 
      YY=ROLL*PI 
      ZZ=PIT*PI 
      COSY=COS(YY) 
      SINY=SIN(YY) 
      COSZ=COS(ZZ) 
      SINZ=SIN(ZZ) 
      COSX=COS(XX) 
      SINX=SIN(XX) 
   20 CONTINUE 
      NT=NP-1 
      IKK=IK+1 
      IK=IK+1 


!     SET ERROR CODES, IF NECESSARY.                                    


      IF(IKK.LE.IABC)GO TO 30 
      SW=1 
   30 CONTINUE 

      IF(NC.EQ.0)GO TO 40 
      IDAV=1 
      NC=-SW1 
      IF(SW.EQ.0.)GO TO 50 
      ICORE=(25+LZ+4*JJJ)*IKK 
      NC=-(SW+SW1) 
   40 CONTINUE 
   50 CONTINUE 

      DO 60 J=1,NP 
      X21(J)=X(J) 
      Y21(J)=Y(J) 
      Z21(J)=Z(J) 
      IIA(J)=0 
   60 END DO 


!     STORE COORDINATES AND SET PEN POSITION WHENEVER ABS(Z)=9999.      


      DO 70 J=1,NT 
      IF(Z21(J).NE.9999.)GO TO 70 
      IIA(J)=1 
      IXU=J-2 
      IBB=J-ISIGN(1,IXU) 
      X21(J)=X21(IBB) 
      Y21(J)=Y21(IBB) 
      Z21(J)=Z21(IBB) 
   70 END DO 

      IIA(NP)=1 
      Z21(NP)=Z21(NT) 
      Y21(NP)=Y21(NT) 
      X21(NP)=X21(NT) 
      JXX=IKK 
      I=1 
      VL=ABS(VP) 


!     LOOP THAT DOES THE THREE DIMENSIONAL TRANSFORMATION ON THE        
!     COORDINATES.                                                      

      JV=L14+(IKK-1)*4*JJJ 
      JT=1 
      JX2=JXX+L2 
      JX3=JXX+L3 
      JX4=JXX+L4 
      JX5=JXX+L5 
      JX6=JXX+L6 
      JX7=JXX+L7 
      JX8=JXX+L8 
      JV1=JV+1 
      JV2=JV+2 
      JV3=JV+3 

      DO 90 J=1,NP 
      XJ=X21(J) 
      YJ=Y21(J) 
      ZJ=Z21(J) 
      X21(J)=ZJ*(COSY*SINX)+XJ*(COSY*COSX)-YJ*SINY 
      TW=YJ*COSY*COSZ 
      TZ=XJ*(SINZ*SINX+SINY*COSZ*COSX) 
      TY=ZJ*(-SINZ*COSX+SINY*COSZ*SINX) 
      Y21(J)=TZ+TW+TY 
      PT=YJ*COSY*SINZ 
      PK=ZJ*(COSZ*COSX+SINY*SINZ*SINX) 
      PS=XJ*(-COSZ*SINX+SINY*SINZ*COSX) 
      Z21(J)=PK+PS+PT 
      RV(JX7)=AMAX1(RV(JX7),Y21(J)) 
      RVI(JX8)=AMIN1(RVI(JX8),Y21(J)) 
      TGI(JX6)=AMAX1(TGI(JX6),X21(J)) 
      TGM(JX5)=AMIN1(TGM(JX5),X21(J)) 
      ZM(JX2)=AMAX1(ZM(JX2),Z21(J)) 
      ZMI(JX3)=AMIN1(ZMI(JX3),Z21(J)) 
      COORD(JV+JT)=X21(J) 
      COORD(JV1+JT)=Y21(J) 
      COORD(JV2+JT)=Z21(J) 
      COORD(JV3+JT)=IIA(J) 
      JT=JT+4 
   90 END DO 

      NOCT(L9+IKK)=NOCT(L9+IKK)+NP 
      NS=NP 
      AVA=AVA+(TGI(JX6)-TGM(JX5))*(RV(JX7)-RVI(JX8)) 
      IF(IXXX.EQ.1)GO TO 95 


!     CALL SUBROUTINE WHICH CALCULATES BOTH THE EQUATIONS OF THE LINE   
!     SEGMENTS AND POLYGONS.                                            

      CALL COEF(X21,Y21,Z21,XXX,JXX,NC,NS,CCC,LZ) 

!     CHECKS TO SEE IF ALL ELEMENTS(SETS) HAVE BEEN PASSED.             

   95 CONTINUE 
      IF(IDAV.EQ.1)GO TO 100 
      GO TO 400 
  100 CONTINUE 
      AVA=AVA/IKK 

      DO 1301 J=1,200 
      ICCT(J)=0 
      ICT(J)=0 
      RCT(J)=J-1 
      IBEG(J)=1 
      IEND(J)=0 
 1301 END DO 

      AMAXX=-999999. 
      AMAXY=-999999. 
      AMINX=999999. 
      AMINY=999999. 
      DO 1400 J=1,IKK 
      AMAXX=AMAX1(AMAXX,TGI(L6+J)) 
      AMAXY=AMAX1(AMAXY,RV(L7+J)) 
      AMINX=AMIN1(AMINX,TGM(L5+J)) 
      AMINY=AMIN1(AMINY,RVI(L8+J)) 
 1400 END DO 

      IAUG=50+((IKK/10000)*2) 
      TMAX=(AMAXX-AMINX)*(AMAXY-AMINY) 
      IBL=TMAX/AVA 
      IBL=IBL/4 

!     DETERMINES THE NUMBER OF GRID POINTS IN THE GRID.                 


      IF(IXXX.EQ.1)GO TO 19990 
      EN=IKK 
      K=(ALOG(EN)/ALOG(2.))+.01 
      K=K+IAUG 
      K=MIN0(K,IBL) 
      IF(K.LE.1)K=1 
      IF(K.GE.199)K=199 
      K=MAX0(K,4) 
      T=K 
      R=(T**.5) 
      KS=R+.5 
      S=T/KS 
      MS=S+.5 
      N=KS*MS 
      MND=N+1 
      XMD=MND 
      T=3./(MND-1) 
      IGY=T*IKK 
      K=KS 
      K1=MS 
      CRX=(AMAXX-AMINX)/K 
      CRY=(AMAXY-AMINY)/K1 


!     DETERMINES THE RELEVANT ELEMENTS VIA THE GRID BLOCKS.             

      DO 93 J=1,IKK 
      IA=0 
      XMAT=TGI(L6+J) 
      XMIT=TGM(L5+J) 
      YMAT=RV(L7+J) 
      YMIT=RVI(L8+J) 
      M=0 
      DO 91 I=1,K1 
      A1=YMAT-(I*CRY+AMINY) 
      A=A1+CRY 
      B1=YMIT-(I*CRY+AMINY) 
      B=B1+CRY 
      A2=A*A1 
      B2=B*B1 
      DO 92 L=1,K 
      M=M+1 
      S1=XMAT-(L*CRX+AMINX) 
      S=S1+CRX 
      R1=XMIT-(L*CRX+AMINX) 
      R=R1+CRX 
      IF((S.LT.0.).OR.(R1.GT.0.))GO TO 92 
      IF((A.LT.0.).OR.(B1.GT.0.))GO TO 92 
      IF((S*S1.GT.0.).OR.(R*R1.GT.0.))GO TO 94 
      IF((A2.GT.0.).OR.(B2.GT.0.))GO TO 94 
      GO TO 93 
   94 CONTINUE 
      ICCT(M)=ICCT(M)+1 
   92 END DO 
   91 END DO 
   93 END DO 

      IADR(1)=0 
      MM=K*K1 
      I8=3*IKK 
      DO 8977 I=2,MM 
      IADR(I)=ICCT(I-1)+IADR(I-1) 
      IF(IADR(I)+ICCT(I).LE.I8) GO TO 8977 

      DO 6666 K9=I,MM 
      ICCT(K9)=-1 
 6666 END DO 
      GO TO 5555 
 8977 END DO 

 5555 CONTINUE 
      DO 191 J=1,MM 
      IF(ICCT(J).GE.0)ICCT(J)=0 
  191 END DO  

      DO 3 J=1,IKK 
      IA=0 
      XMAT=TGI(L6+J) 
      XMIT=TGM(L5+J) 
      YMAT=RV(L7+J) 
      YMIT=RVI(L8+J) 
      J16=L16+J 
      M=0 
      DO 1 I=1,K1 
      A1=YMAT-(I*CRY+AMINY) 
      A=A1+CRY 
      B1=YMIT-(I*CRY+AMINY) 
      B=B1+CRY 
      A2=A*A1 
      B2=B*B1 
      DO 2 L=1,K 
      M=M+1 
      S1=XMAT-(L*CRX+AMINX) 
      S=S1+CRX 
      R1=XMIT-(L*CRX+AMINX) 
      R=R1+CRX 
      IF((S.LT.0.).OR.(R1.GT.0.))GO TO 2 
      IF((A.LT.0.).OR.(B1.GT.0.))GO TO 2 
      IF((S*S1.GT.0.).OR.(R*R1.GT.0.))GO TO 4 
      IF((A2.GT.0.).OR.(B2.GT.0.))GO TO 4 
      XIND(J16)=M 
      GO TO 3 
    4 CONTINUE 
      IA=IA+1 
      IF(IA.LE.4)GO TO 8000 
      XIND(J16)=0 
      GO TO 8001 
 8000 CONTINUE 
      XIND(J16)=XIND(J16)+M*(MND**(IA-1)) 
 8001 CONTINUE 
      IF(ICCT(M).LT.0)GO TO 2 
      ICCT(M)=ICCT(M)+1 
      JK=IADR(M)+ICCT(M)+L17 
      NEH(JK)=J 
    2 END DO 
    1 END DO 
    3 END DO 

!!!      CALL VSRT1(XIND(L16+1),IK,IND(L15+1)) 
      CALL VSRT1(XIND(L16+1:L16+IK),IK,IND(L15+1:L15+IK)) 
      SW=0 
      L=1 
      DO 5 I=1,IKK 
   11 CONTINUE 
      IF(XIND(L16+I).NE.RCT(L))GO TO 6 
      SW=SW+1 
      IF(SW.EQ.1.)LT=I 
      ICT(L)=ICT(L)+1 
      GO TO 5 
    6 CONTINUE 
      IF(SW.NE.0.)GO TO 8 
      L=L+1 
      GO TO 11 
    8 CONTINUE 
      IBEG(L)=LT 
      IEND(L)=LT+ICT(L)-1 
      SW=0 
      IF(XIND(L16+I).GE.MND)GO TO 2110 
      L=L+1 
      GO TO 11 
    5 END DO 

      IBEG(L)=LT 
      IEND(L)=LT+ICT(L)-1 
 2110 CONTINUE 
      DO 2111 J=1,IKK 
      SNDT(L4+J)=IND(L15+J) 
 2111 END DO 
!!!      CALL VSRTR(SNDT(L4+1),IK,XIND(L16+1)) 
      CALL VSRTR(SNDT(L4+1:L4+IK),IK,XIND(L16+1:L16+IK)) 

      EN=IKK 
      IGX=(ALOG(EN)/ALOG(2.))+1. 
      DO 105 J=1,IGX 
      RRX(J)=2**(IGX-J) 
  105 END DO 

19990 CONTINUE 
      IJ=0 
      X1(3)=(AMAXX+AMINX)/2 
      Y1(3)=(AMAXY+AMINY)/2 
      X1(4)=SCX 
      Y1(4)=SCX 
      IF(IXXX.EQ.1)GO TO 18880 
      DO 115 J=1,IKK 
      IN(L11+J)=J 
      IN1(L0+J)=J 
      TGMT(L10+J)=TGM(L5+J) 
      YMIN(L1+J)=RVI(L8+J) 
  115 END DO 

!     CALL SUBROUTINE WHICH WILL SORT ON X,Y AND Z.                     

!!!      CALL VSRTR(TGMT(L10+1),IK,IN(L11+1)) 
!!!      CALL VSRTR(YMIN(L1+1),IK,IN1(L0+1)) 
      CALL VSRTR(TGMT(L10+1:L10+IK),IK,IN(L11+1:L11+IK)) 
      CALL VSRTR(YMIN(L1+1:L1+IK),IK,IN1(L0+1:L0+IK)) 
      H(8)=0 
18880 CONTINUE 
      L10W=L10-1 
      L1W=L1-1 
      L01W=L01-1 
      TZQ=.05*IKK
 
      DO 395 J=1,IKK   ! start of a really big loop
      N9=0 
      IF(ISILH(J).NE.0)N9=2 
      IF((IJBB.EQ.0).AND.(J.GE.TZQ))N9=0 
      N7=J-1 
      JXT=J 
      KS=IKK 
      JJ=L14+N7*4*JJJ 
      JH=1 
      II=0 
      IXR=NOCT(L9+J) 
      NIT=0 
      JT=L12+5*N7 
      JT1=JT 
      JO=L13+LZ*N7 
      IF(IXXX.EQ.1)GO TO 200 
      NS=XXX(5+JT) 
      DO 1910 IC=1,NS 
 1910 IIA(400+IC)=0 
      D1=XXX(JT+1) 
      D2=XXX(JT+2) 
      D3=XXX(JT+3) 
      D4=XXX(JT+4) 
      IF(D3.EQ.0.)GO TO 108 
      Q1=D1/D3 
      Q2=D2/D3 
      Q4=D4/D3 
  108 CONTINUE 
      I9=0 
      NG=NS*5 
      I=0 
      JD=1 
      JW=JJ+1 
      JW1=JJ+2 
      DO 121 I=1,NS 
      XE(I)=COORD(JJ+JD) 
      YE(I)=COORD(JW+JD) 
      DI(I)=COORD(JW1+JD) 
      JD=JD+4 
  121 END DO 
  122 CONTINUE 
!     THOSE OTHER ELEMENTS WHICH COULD POSSIBLY HIDE SOME PORTION       
!     OF THE GIVEN ELEMENT.                                             



      K=2**IGX 
      K1=K 
      K2=K 

!     DO LOGARITHMIC SEARCH TO DETERMINE RELEVANT ELEMENTS.             

      TGII=TGI(L6+J) 
      RVV=RV(L7+J) 
      ZMII=ZMI(L3+J) 
      TGMM=TGM(L5+J) 
      RVII=RVI(L8+J) 
      S=-1 
      DO 131 I=1,IGX 
      K=K+SIGN(RRX(I),S) 
      IF(K.GT.IKK)K=IKK 
      S=TGII-TGMT(L10+K) 
      IF(S*(TGII-TGMT(L10W+K)).LT.0.)GO TO 132 
  131 END DO 
      K=IKK 
  132 CONTINUE 
      S=-1 
      DO 133 I=1,IGX 
      K1=K1+SIGN(RRX(I),S) 
      IF(K1.GT.IKK)K1=IKK 
      S=RVV-YMIN(L1+K1) 
      IF(S*(RVV-YMIN(L1W+K1)).LT.0.)GO TO 134 
  133 END DO 
      K1=IKK 
  134 CONTINUE 

!     RETRIEVE THE RELEVANT ELEMENTS DETERMINED FROM SCHEME 1.          

 8181 CONTINUE 
      IR=XIND(L16+J) 
      IF(IR.EQ.0)GO TO 1270 
      VX=IR 
      T=ALOG(VX) 
      IF(IR.LE.LC)GO TO 1800 
      E=LC 
      LG=IR/LC 
      MU=MOD(IR,LC) 
      UX=LG+(MU/E) 
      T=ALOG(UX)+ALOG(E) 
 1800 CONTINUE 
      IXT=0 
      IEXP=(T/ALOG(XMD))+1 
      DO 8004 L=1,IEXP 
      JCZ=MND**(IEXP-L) 
      IV=IR/JCZ 
      IR=IR-IV*JCZ 
      IV=IV+1 
      IV1=IV-1 
      IF(ICCT(IV1).EQ.0)GO TO 4000 
      IF(ICCT(IV1).GT.0)GO TO 4001 
      GO TO 1270 
 4001 CONTINUE 
      KE=ICCT(IV1) 
      IL=0 
      JTT=IADR(IV1)+L17 
      JJG=L4+IXT 
      DO 4003 I=1,KE 
      KV=NEH(I+JTT) 
      IF(KEEP(L18+KV).EQ.J)GO TO 4003 
      IL=IL+1 
      NNO(JJG+IL)=KV 
      KEEP(L18+KV)=J 
 4003 END DO 
      IXT=IXT+IL 
 4000 CONTINUE 
      IX=IBEG(IV) 
      IX1=IEND(IV) 
      LJ=L4+IXT-IX+1 
      DO 1170 I=IX,IX1 
 1170 NNO(LJ+I)=IND(L15+I) 
      IXT=IXT+IX1-IX+1 
 8004 END DO 
      KS=IXT 
 1270 CONTINUE 
      IM=MIN0(K,K1) 

!     PICK MINIMUM COUNT FROM BOTH SCHEMES.                             

      IF(KS.LT.IM)GO TO 129 
      IF(IM.EQ.K)GO TO 1001 
      IF(IM.EQ.K1)GO TO 1002 
                                    ! I1 has never been defined - RLC 
      KS=I1                         ! this statement may not be correct  - RLC
      IJ1=L00+IKK+1 
      DO 1003 I=1,KS 
 1003 NNO(L4+I)=IN2(IJ1-I) 
      GO TO 129 
 1001 CONTINUE 
      KS=K 
      DO 1004 I=1,KS 
 1004 NNO(L4+I)=IN(L11+I) 
      GO TO 129 
 1002 CONTINUE 
      KS=K1 
      DO 1006 I=1,KS 
 1006 NNO(L4+I)=IN1(L0+I) 
  129 CONTINUE 
      DO 170 I=1,KS    ! start of big loop
      JB=NNO(L4+I) 
      IF(ISILH(J).EQ.0)GO TO 265 
      IF(ISILH(J).NE.ISILH(JB))GO TO 265 
      IF((RVV.LT.RVI(L8+JB)).OR.(RVII.GT.RV(L7+JB)))GO TO 170 
      IF((TGMM.GT.TGI(L6+JB)).OR.(TGII.LT.TGM(L5+JB)))GO TO 170 
      JS=L12+(JB-1)*5 
      IF(XXX(3+JS).EQ.0)GO TO 170 
      IF(J.EQ.JB)GO TO 166 
      GO TO 2399 
  265 CONTINUE 
      IF((TGMM.GE.TGI(L6+JB)).OR.(TGII.LE.TGM(L5+JB)))GO TO 170 
      IF((RVV.LE.RVI(L8+JB)).OR.(RVII.GE.RV(L7+JB)))GO TO 170 
      IF(ZMII.GE.ZM(L2+JB))GO TO 170 
      IF(J.EQ.JB)GO TO 170 
      JS=L12+(JB-1)*5 
      IF(XXX(3+JS).EQ.0.)GO TO 170 
 2399 CONTINUE 
      IF(D3.EQ.0.)GO TO 166 
      NV=XXX(5+JS) 

!     TEST TO SEE IF ALL VERTICES LIE EITHER BEHIND OR IN FRONT OF      
!     THE GIVEN POLYGON.                                                

      IT=0 
      JD=1 
      JI=L14+(JB-1)*4*JJJ 
      JI1=JI+1 
      JI2=JI+2 
      DO 145 M=1,NV 
      XU(M)=COORD(JI+JD) 
      YU(M)=COORD(JI1+JD) 
      XI(M)=COORD(JI2+JD) 
      JD=JD+4 
      ZS1=-(Q4+Q2*YU(M)+Q1*XU(M)) 
      IF (ABS(XI(M)-ZS1).LT.HXX) CYCLE   !!! GO TO 145 
      IT=IT+1 
      ICOUNT(IT)=0 
      IF(XI(M).GT.ZS1)ICOUNT(IT)=1 
  145 END DO 


!     TESTS FOR SEMI-RELEVANT PLANES.  THAT IS,NEGATIVE INDEXES         
!     INDICATE ELEMENT IS TO BE USED FOR VISIBILITY TEST, BUT NOT FOR   
!     INTERSECTION LINE DETERMINATION.                                  


!     IF ALL EDGES HAVE ALREADY BEEN DRAWN,THEN DISREGARD THE FOLLOWING 

      IF((ISILH(J).NE.ISILH(JB)).AND.(ISILH(J).NE.0))GO TO 1212 
      IF(J.LT.JB)GO TO 1650 
      IF(I9.EQ.NS)GO TO 1650 
      MG=0 
      NV1=NV-1 
      NS1=NS-1 
      DO 1630 L=1,NS 
      DO 1580 IX=1,NV 
      IF (ABS(XI(IX)-DI(L)).GT.HXX) CYCLE   !!! GO TO 1580 
      IF ((ABS(YU(IX)-YE(L)).GT.HXX).OR.(ABS(XU(IX)-XE(L)).GT. HXX)) CYCLE   !!! GO TO 1580                                                   
      IF(MG.EQ.2)GO TO 1640 
      MG=MG+1 
      I2(MG)=L 
      I3(MG)=IX 
      GO TO 1630 
 1580 END DO 
 1630 END DO 
 1640 CONTINUE 
      IF(MG.NE.2)GO TO 1650 
      IR=IABS(I3(2)-I3(1)) 
      IR1=IABS(I2(2)-I2(1)) 
      IF((IR.NE.1).AND.(IR.NE.NV1))GO TO 1650 
      IF((IR1.NE.1).AND.(IR1.NE.NS1))GO TO 1650 
      IX=MAX0(I2(1),I2(2)) 
      IF(IR1.EQ.1)IX=IX-1 
      IIA(400+IX)=1 
      I9=I9+1 
 1650 CONTINUE 
 1212 CONTINUE 
      L=0 
      DO M=1,IT 
        L=L+ICOUNT(M) 
      END DO

      IF((ISILH(J).NE.0).AND.(ISILH(JB).EQ.ISILH(J)))GO TO 165 
      IF(L.EQ.0)GO TO 170 
      IF(II.GT.N9)GO TO 165 


!     INTERROGATE THE RELATIONSHIP OF THE CANDIDATE POLYGON TO THE      
!     GIVEN POLYGON BY DETERMINING IF THE PROJECTION OF ONE POLYGON     
!     CAN BE SEPARATED BY AN EDGE FROM THE OTHER]S PROJECTION           


      JK=L13+(JB-1)*LZ 
      E3=XXX(3+JS) 
      E1=XXX(1+JS) 
      E4=XXX(4+JS) 
      E2=XXX(2+JS) 
      SD=0 
      I3(1)=JK 
      I3(2)=JO 
      I2(1)=NV*5 
      I2(2)=NS*5 
      DO 164 KU=1,2 
      IS=I3(KU) 
      IB=I2(KU) 
      DO 163 LL=1,IB,5 
  151 CONTINUE 
      IF(SD.EQ.1.)GO TO 152 
      A=D2*E3-E2*D3 
      B=D1*E3-E1*D3 
      C=D4*E3-E4*D3 
      GO TO 153 
  152 CONTINUE 
      A=CCC(LL+IS) 
      B=CCC(LL+IS+1) 
      C=CCC(LL+IS+2) 
  153 CONTINUE 
      IF((A.EQ.0.).AND.(B.EQ.0.))GO TO 162 
      IF(A.NE.0.)GO TO 154 
      A=0 
      C=C/B 
      B=1 
      GO TO 155 
  154 CONTINUE 
      B=B/A 
      C=C/A 
      A=1 
  155 CONTINUE 
      M=0 
      R1=0 
      DO 158 IX=1,NV 
      M=M+1 
      YG=YU(M) 
      IF(A.NE.0.)GO TO 156 
      DY=-C/B 
      YG=XU(M) 
      GO TO 157 
  156 CONTINUE 
      DY=-C-B*XU(M) 
  157 IF(ABS(DY-YG).LT.HXX)GO TO 158 
      R=YG-DY 
      IF(R*R1.LT.0.)GO TO 162 
      R1=R 
  158 END DO 
      M=0 
      R2=0 
      DO 161 IX=1,NS 
      M=M+1 
      YG=YE(M) 
      IF(A.NE.0.)GO TO 159 
      DY=-C/B 
      YG=XE(M) 
      GO TO 160 
  159 CONTINUE 
      DY=-C-B*XE(M) 
  160 IF(ABS(DY-YG).LT.HXX)GO TO 161 
      R=YG-DY 
      IF(R*R2.LT.0.)GO TO 162 
      R2=R 
  161 END DO 
      IF(R1*R2.LT.0.)GO TO 170 
  162 CONTINUE 
      IF(SD.NE.0.)GO TO 163 
      SD=1 
      GO TO 151 
  163 END DO 
  164 END DO 
  165 CONTINUE 
      IF((L.EQ.IT).OR.(L.EQ.0))JB=-JB 
  166 CONTINUE 
      II=II+1 
      NNO(L4+II)=JB 
  170 END DO                           ! end of big loop

      IF(II.LE.2)IJBB=1 
      JAT=-4 
      IF(IXR.LE.3)GO TO 200 
      IF(II.EQ.0)GO TO 190 


!     CALL SUBROUTINE WHICH SOLVES FOR THE LINES OF INTERSECTION,IF ANY,
!     OF THE JTH ELEMENT WITH OTHER ELEMENTS.                           


      CALL SOLVE(IXR,J,XXX,CCC,II,NNO,NIT,X21,Y21,Z21,IIA,NC,ZM,ZMI,LZ) 
  190 CONTINUE 
  200 CONTINUE 
      IJ1=JJ+1 
      IJ2=JJ+2 
      IJ3=JJ+3 
      DO 210 JM=1,IXR 
      X21(JM)=COORD(JH+JJ) 
      Y21(JM)=COORD(JH+IJ1) 
      Z21(JM)=COORD(JH+IJ2) 
      IIA(JM)=COORD(JH+IJ3) 
      JH=JH+4 
  210 END DO 
      IXR=IXR+3*NIT 
      IF(II.EQ.0)GO TO 220 
      IF(IXXX.NE.1)GO TO 240 
  220 CONTINUE 
      DO 230 JM=1,IXR 
      IF(IXXX.EQ.1)GO TO 1993 
      IF(JM.GE.IXR-1)GO TO 1993 
      IF(IIA(400+JM).EQ.1)IIA(JM+1)=1 
 1993 CONTINUE 
      X1(2)=X21(JM) 
      Y1(2)=Y21(JM) 
      IM=IABS(IIA(JM)) 
      CALL PLT(X1,Y1,IJ,IM) 
      IF((JM.NE.IXR).AND.(IIA(JM+1) .EQ.0))CALL PLT(X1,Y1,IJ,0) 
  230 END DO 
      GO TO 390 
  240 CONTINUE 
      JX=1 
  250 CONTINUE 


!     PLOTS IF IIA(JX+1) IS EQUAL TO 1.                                 

      IF((IIA(JX).EQ.0).AND.(IIA(JX+1).EQ.0))GO TO 260 
      IF(IIA(JX+1).NE.-1)GO TO 251 
      IIA(JX+1)=0 
      IM=1 
      JAT=JAT+5 
      JX=JX-1 
      GO TO 252 
  251 CONTINUE 
      IM=IIA(JX+1) 
  252 CONTINUE 
      X1(2)=X21(JX+1) 
      Y1(2)=Y21(JX+1) 
      CALL PLT(X1,Y1,IJ,IM) 
      JX=JX+2 
      IF(JX.GE.IXR)GO TO 390 
      GO TO 250 
  260 CONTINUE 
      IJ=0 
      JAT=JAT+5 
      ME=0 
      IF(JX.GT.NS)GO TO 255 
      IF(IIA(400+JX).EQ.1)GO TO 381 
  255 CONTINUE 

!     CALL SUBROUTINE WHICH DETERMINES THE POINTS OF INTERSECTIONS      
!     OF THE LINES OF THE JTH SET WITH THE RELEVANT LINES AND PLANES    
!     OF OTHER ELEMENTS.                                                

      H(6)=NIT 
      CALL CHECK(XXX,CCC,NNO,J,II,NC,XI,YI,NGX,ZM,ZMI,RV,RVI,TGM,TGI,   &
     &ZI,LZ)                                                            
      NG=NGX(1)+2 
      IMB=JX+1 
      XI(1)=X21(JX) 
      YI(1)=Y21(JX) 
      ZI(1)=Z21(JX) 
      XI(NG)=X21(IMB) 
      YI(NG)=Y21(IMB) 
      ZI(NG)=Z21(IMB) 
  880 CONTINUE 
 1101 CONTINUE 


!     THE FOLLOWING CODE SORTS THE INTERSECTION POINTS IN ASCENDING     
!     ORDER OF OCCURENCE AND THEN SHRINKS THE LIST IF REDUNDANCY EXIST. 

      IF(NG.LE.3)GO TO 340 
      NI=NG-2 
      NII=NI 
      DO 270 M=1,NG 
      DI(M)=(XI(M)-XI(1))**2+(YI(M)-YI(1))**2 
  270 END DO 
      DO 290 M=2,NI 
      DO 280 MX=2,NII 
      IF (DI(MX).LE.DI(MX+1)) CYCLE   !!! GO TO 280 
      IW=MX+1 
      HOLD=DI(MX) 
      HOLD1=XI(MX) 
      HOLD2=YI(MX) 
      HOLD3=ZI(MX) 
      XI(MX)=XI(IW) 
      YI(MX)=YI(IW) 
      ZI(MX)=ZI(IW) 
      DI(MX)=DI(IW) 
      DI(IW)=HOLD 
      XI(IW)=HOLD1 
      YI(IW)=HOLD2 
      ZI(IW)=HOLD3 
  280 END DO 
      NII=NII-1 
  290 END DO 
      LX=1 
      NPX=NG 
  300 NPX=NPX-1 
      I=LX 
      DO 320 M=I,NPX 
      RX=0 
      T=SQRT((XI(M)-XI(M+1))**2+(YI(M)-YI(M+1))**2) 
      IF(T.GT.HX1)GO TO 320 
      IX=M 
      IX1=NPX 
      DO 310 MX=IX,IX1 
      IW=MX+1 
      XI(MX)=XI(IW) 
      YI(MX)=YI(IW) 
      ZI(MX)=ZI(IW) 
  310 END DO 
      RX=1 
      LX=M 
      IF(LX.EQ.NPX)GO TO 330 
      GO TO 300 
  320 END DO 
  330 CONTINUE 
      IF(RX.EQ.1.)NPX=NPX-1 
      NG=NPX+1 
  340 CONTINUE 


!     THIS CODE DETERMINES THE STATUS(VISIBILITY) OF EVERY OTHER POINT  
!     AS SUGGESTED BY THE THEOREM IN THE TECHNICAL REPORT.              

      DO 350 L=1,NG,2 
      OJ=XI(L) 
      TMJ=YI(L) 
      ZJ=ZI(L) 
      CALL STATUS(OJ,TMJ,XXX,TGM,RV,RVI,TGI,ZM,NNO,II,H,IM,JXT,         &
     &ZJ,NC,ZMI,CCC,LZ)                                                 
      DI(L)=IM 
  350 END DO 
      JII=NG-1 
      DO 370 L=1,NG,2 
      IF(L.EQ.NG)GO TO 370 
      IF(L.EQ.JII)GO TO 360 
      IF(DI(L)+DI(L+2).NE.2.)GO TO 360 
      DI(L+1)=DI(L) 
      GO TO 370 
  360 CONTINUE 
      MN=L+1 
      OJ=XI(MN) 
      TMJ=YI(MN) 
      ZJ=ZI(MN) 
      CALL STATUS(OJ,TMJ,XXX,TGM,RV,RVI,TGI,ZM,NNO,II,H,IM,             &
     &JXT,ZJ,NC,ZMI,CCC,LZ)                                             
      DI(MN)=IM 
  370 END DO 


!     THE FOLLOWING CODE ACTUALLY PLOTS THE POINTS ON A GIVEN LINE      
!     GOVERNED BY THE VALUE(IM) RETURNED BY STATUS SUBROUTINE.          
!     1 MEANS INVISIBLE,...0 MEANS VISIBLE.                             


  379 CONTINUE 
      DO 380 L=1,NG 
      X1(2)=XI(L) 
      Y1(2)=YI(L) 
      IM=DI(L) 
      CALL PLT(X1,Y1,IJ,IM) 
      IF(L.EQ.NG)GO TO 380 
      IF(DI(L)+DI(L+1).GT.0.)GO TO 380 
      H(8)=1 
      MN=L+1 
      OJ=(XI(L)+XI(MN))/2 
      TMJ=(YI(L)+YI(MN))/2 
      ZJ=(ZI(L)+ZI(MN))/2 
      CALL STATUS(OJ,TMJ,XXX,TGM,RV,RVI,TGI,ZM,NNO,II,H,IM,JXT,ZJ,NC,   &
     &ZMI,CCC,LZ)                                                       
      H(8)=0 
      IF(IM.EQ.0)GO TO 380 
      X1(2)=OJ 
      Y1(2)=TMJ 
      CALL PLT(X1,Y1,IJ,IM) 
  380 END DO 
  381 CONTINUE 
      JX=JX+1 
      GO TO 250 
  390 CONTINUE 

!     DECREMENTS THE COUNT OF THE NUMBER OF LINES IN THE JTH SET        
!     SINCE THE LINES OF INTERSECTIONS WERE ADDED TO THIS ELEMENT       
!     BY THE SUBROUTINE SOLVE.

      XXX(5+JT)=XXX(5+JT)-NIT 
  395 END DO    ! end of a really big loop



  400 CONTINUE 
      RETURN 
      END Subroutine Lin   ! --------------------------------

!+
      SUBROUTINE STATUS(OJ,TMJ,XXX,TGM,RV,RVI,                          &
     & TGI,ZM,NNO,II,H,IM,JXT,ZJ,NC,ZMI,CCC,LZ)                          
! ------------------------------------------------------------------------------
! PURPOSE - This subroutine determines the visibility of an arbitrary point by
!  drawing a line from the point in question to infinity and counting the
!  number of times it crosses the boundaries of a relevant element.

      DIMENSION CCC(*),XXX(*) 
      DIMENSION ZMI(*),TGM(*),RV(*),RVI(*)
      DIMENSION TGI(*),ZM(*),NNO(*),H(*)                                          
      COMMON/GO3/L0,L1,L00,L01,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13 
!!!      COMMON/DRH/ISILH(1)                                            
      INCLUDE 'drh.cmn' 
      DATA GGK,GGJ,EI/.005,.015,.2/ 
      IM=0 
      D=EI*OJ-TMJ 
      DO 60 JO=1,II 
      JG=NNO(L4+JO) 

!     PRELIMINARY CHECK TO SEE IF THE POINT IS OUTSIDE THE BOUNDARY     
!     BOXES IN THE X,Y,Z DIMENSIONS.                                    

      IF((TMJ.GE.RV(L7+JG)).OR.(TMJ.LE.RVI(L8+JG)))GO TO 60 
      IF((OJ.GE.TGI(L6+JG)).OR.(OJ.LE.TGM(L5+JG)))GO TO 60 
      IF(ZJ.GE.ZM(L2+JG))GO TO 60 
      JT=L12+(JG-1)*5 
      ZS=-(XXX(4+JT)+XXX(2+JT)*TMJ+XXX(1+JT)*OJ)/XXX(JT+3) 
      IF(ABS(ZJ-ZS).LT.GGJ)GO TO 60 
      IF((ZJ.GT.ZS).OR.(JXT.EQ.JG))GO TO 60 
      IB=XXX(5+JT)*5 
      JS=L13+(JG-1)*LZ 
      DO 20 J=1,IB,5 
      JJ=JS+J 
      IF(CCC(JJ).EQ.0)GO TO 15 
      S=OJ-CCC(JJ+3) 
      S1=OJ-CCC(JJ+4) 
      DY=(-CCC(JJ+2)-CCC(JJ+1)*OJ)-TMJ 
      GO TO 16 
   15 CONTINUE 
      DY=(-CCC(JJ+2)/CCC(JJ+1))-OJ 
      S=TMJ-CCC(JJ+3) 
      S1=TMJ-CCC(JJ+4) 
   16 CONTINUE 
      IF((ABS(DY).LT.GGK).AND.(S*S1.LE.0.))GO TO 61 
   20 END DO 


!     THE FOLLOWING CODE COUNTS THE INTERSECTIONS OF BOUNDARIES         
!     OF A GIVEN ELEMENT WITH THE INFINITE LINE AND CHECKS,IF INSIDE    
!     OF THE BOUNDARY, WHETHER OR NOT THE POINT IS BEHIND OR IN FRONT   
!     OF THE ELEMENT.                                                   


      I=0 
      DO 40 J=1,IB,5 
      JJ=JS+J 
      R=EI*CCC(JJ)+CCC(JJ+1) 
      IF(R.EQ.0.)GO TO 40 
      T=(-CCC(JJ+2)+CCC(JJ)*D)/R 
      IF(T.LT.OJ)GO TO 40 
      IF(CCC(JJ).NE.0.)GO TO 30 
      T=EI*T-D 
   30 CONTINUE 
      IF((T-CCC(JJ+3))*(T-CCC(JJ+4)).GT.0.)GO TO 40 
   35 CONTINUE 
      I=I+1 
   40 END DO 
      IF(I-(I/2)*2.EQ.0)GO TO 60 
   50 CONTINUE 
      IM=1 
      GO TO 70 
   61 CONTINUE 
      IF(H(8).NE.1)GO TO 60 
      IF((ISILH(JXT).EQ.ISILH(JG)).AND.(ISILH(JXT).NE.0))GO TO 60 
      IF(ZJ.LT.ZMI(L3+JG))GO TO 50 
   60 END DO 
      IF((ISILH(JXT).EQ.0).OR.(H(8).NE.1.))GO TO 70 
      CALL STATUM(OJ,TMJ,XXX,TGM,RV,RVI,TGI,ZM,NNO,II,H,IM,JXT,         &
     &ZJ,NC,ZMI,CCC,LZ)                                                 
   70 CONTINUE 
      RETURN 
      END Subroutine Status   ! -----------------------------------------------

!+
      SUBROUTINE PLT(X1,Y1,IJ,IM) 
! ------------------------------------------------------------------------------
! PURPOSE - Plots points governed by the value of im.

!        NOTE THAT CALL PLOT(X,Y,2) MEANS MOVE PEN FROM THE CURRENT     
!     POSITION TO THE POINT,(X,Y),WITH THE PEN DOWN.                    


!     CALL PLOT(X,Y,3) MEANS MOVE THE PEN FROM THE CURRENT POSITION     
!     TO THE POINT,(X,Y), WITH THE PEN UP.                              

      DIMENSION X1(4),Y1(4)  ! changed from 1 to 4 by RLC  17Aug95    
!-------------------------------------------------------------------------------
      IF(IM.EQ.1)GO TO 20 
      IF(IJ.EQ.0)GO TO 10 

!     HERE IS WHERE THE POINTS ARE DRAWN. YOU MUST ALSO CHECK TO BE     
!     SURE THAT THE POINT AT '10 CONTINUE'IS INCLUDED WHENEVER          
!     POINTS ARE DRAWN;IT WILL BE THE FIRST IN THE SEQUENCE.            

      CALL PLOT((X1(2)-X1(3))/X1(4),(Y1(2)-Y1(3))/Y1(4),2) 
      GO TO 30 
   10 CONTINUE 
      CALL PLOT((X1(2)-X1(3))/X1(4),(Y1(2)-Y1(3))/Y1(4),3) 
      IJ=1 
      GO TO 30 
   20 CONTINUE 
      IJ=0 
   30 CONTINUE 
      RETURN 
      END Subroutine Plt   ! ---------------------------------------------------

!+
      SUBROUTINE CHECK(XXX,CCC,NNO,J,II,NC,XI,YI,                       &
     &                   NGX,ZM,ZMI,RV,RVI,TGM,TGI,ZI,LZ) 
! ------------------------------------------------------------------------------
! PURPOSE - This subroutine solves for the points of intersection on the      
!     lines of the jth element with other lines and planes(relevant)    


      DIMENSION CCC(:),XXX(:) 
      DIMENSION RV(:),RVI(:),TGM(:),TGI(:),ZM(:),ZMI(:)
      DIMENSION NNO(:),NGX(15),XI(:),YI(:),ZI(:)                                  
      DOUBLE PRECISION T8 
      COMMON/DAVE/XCC(800) 
      COMMON/HEDG/JS,M,JT,VX,VX1,VX2,VX3,NN 
      COMMON/GO3/L0,L1,L00,L01,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13 
!!!      COMMON/DRH/ISILH(1)                                            
      INCLUDE 'drh.cmn' 
      DATA EEX,EXP/.015,.005/ 
!-------------------------------------------------------------------------------
      NGX(1)=0 
      XCC3=XCC(JS+3) 
      XCC4=XCC(JS+4) 
      XCC1=XCC(JS+1) 
      XCC2=XCC(JS+2) 
      IF(NN.EQ.1)GO TO 5 
      IF(VX3.NE.0.)GO TO 5 
      A=XXX(JT+2) 
      B=XXX(JT+1) 
      C=XXX(JT+4) 
      Z1=XCC(JS) 
      Z2=XCC1 
      IF(A.EQ.0.)GO TO 1 
      Y1=-XCC3*B-C 
      Y2=-XCC4*B-C 
      X1=XCC3 
      X2=XCC4 
      GO TO 15 
    1 CONTINUE 
      Y1=XCC3 
      Y2=XCC4 
      X1=-C 
      X2=X1 
      GO TO 15 
    5 CONTINUE 
      A=XCC(JS) 
      B=XCC1 
      C=XCC2 
      IF(A.EQ.0.)GO TO 20 
      Y1=-XCC3*B-C 
      Y2=-XCC4*B-C 
      X1=XCC3 
      X2=XCC4 
      GO TO 30 
   20 CONTINUE 
      Y1=XCC3 
      Y2=XCC4 
      X1=-XCC2 
      X2=X1 
   30 CONTINUE 
      IF(NN.NE.1)GO TO 40 
      Z1=XXX(1+JT) 
      Z2=XXX(2+JT) 
      GO TO 50 
   40 CONTINUE 
      Z1=-(VX+VX1*Y1+VX2*X1)/VX3 
      Z2=-(VX+VX1*Y2+VX2*X2)/VX3 
   50 CONTINUE 
   15 CONTINUE 
      AL=X2-X1 
      BL=Y2-Y1 
      CL=Z2-Z1 
      EG=AMIN1(Z1,Z2) 
      EGX=AMAX1(X1,X2) 
      EGX1=AMIN1(X1,X2) 
      EGY=AMAX1(Y1,Y2) 
      EGY1=AMIN1(Y1,Y2) 
      EGZ=AMAX1(Z1,Z2) 
      IF(AL.EQ.0.)GO TO 51 
      BLA=BL/AL 
      BLAX=(BLA)*X1 
      CLA=CL/AL 
      CLAX=(CLA)*X1 
      GO TO 52 
   51 CONTINUE 
      IF(BL.EQ.0.)GO TO 52 
      ALB=AL/BL 
      ALBY=(ALB)*Y1 
      CLB=CL/BL 
      CLBY=(CLB)*Y1 
   52 CONTINUE 


!     THIS CODE DETERMINES THE POINTS OF INTERSECTIONS ON THE LINES OF  
!     JTH ELEMENT RESULTING FROM THE INTERSECTION OF THE PLANES WITH    
!     THESE LINES.                                                      


      DO 170 JR=1,II 
      KM=L4+JR 
      LG=NNO(KM) 
      NNO(KM)=IABS(LG) 
      LE=NNO(KM) 
      IF(J.EQ.LE)GO TO 170 
      IF(EGX.LT.TGM(LE+L5))GO TO 170 
      IF(EGX1.GT.TGI(LE+L6))GO TO 170 
      IF(EGY.LT.RVI(LE+L8))GO TO 170 
      IF(EGY1.GT.RV(LE+L7))GO TO 170 
      IF((ISILH(J).EQ.ISILH(LE)).AND.(ISILH(J).NE.0))GO TO 1172 
      IF(EG.GT.ZM(L2+LE))GO TO 170 
 1172 CONTINUE 
      JE=L13+LZ*(LE-1) 
      JU=L12+5*(LE-1) 
      AC=XXX(1+JU) 
      BC=XXX(2+JU) 
      CC=XXX(3+JU) 
      D=XXX(4+JU) 
      G=1./CC 
      NK=XXX(5+JU)*5 
      IF((LG.LT.0).OR.(EGZ.LE.ZMI(L3+LE)))GO TO 80 
      IF((AL.EQ.0).AND.(BL.EQ.0))GO TO 80 
      IF(AL.EQ.0)GO TO 60 
      VU=AC+BC*BLA+CC*CLA 
      IF(VU.EQ.0.)GO TO 80 
      XP=BC*BLAX+CC*CLAX-D 
      XP=XP-BC*Y1-CC*Z1 
      XP=XP/VU 
      T1=(XP-X1)/AL 
      YP=T1*BL+Y1 
      GO TO 70 
   60 CONTINUE 
      VU=BC+AC*ALB+CC*CLB 
      IF(VU.EQ.0.)GO TO 80 
      YP=AC*ALBY+CC*CLBY-D 
      YP=YP-CC*Z1-AC*X1 
      YP=YP/VU 
      T1=(YP-Y1)/BL 
      XP=T1*AL+X1 
   70 CONTINUE 
      IF((XP-TGM(LE+L5))*(XP-TGI(LE+L6)).GT.0.)GO TO 80 
      IF((YP-RV(LE+L7))*(YP-RVI(LE+L8)).GT.0.)GO TO 80 
      T=XP 
      IF(A.EQ.0.)T=YP 
      IF((T-XCC3)*(T-XCC4).GE.0.)GO TO 80 
      ZP=T1*CL+Z1 
      S=ZP-ZM(L2+LE) 
      S1=ZP-ZMI(L3+LE) 
      IF((ABS(S).LT.EEX).OR.(ABS(S1).LT.EEX))GO TO 56 
      IF(S*S1.GT.0.)GO TO 80 
   56 CONTINUE 

!     STORES INTERSECTIONS.                                             

      M=M+1 
      NQ=M+1 
      XI(NQ)=XP 
      YI(NQ)=YP 
      ZI(NQ)=ZP 
   80 CONTINUE 

!     THIS CODE DETERMINES INTERSECTION POINTS OF LINES WITH LINES.     

      JE1=JE+1 
      JE2=JE+2 
      JE3=JE+3 
      JE4=JE+4 
      DO 160 JV=1,NK,5 
      B1=CCC(JV+JE1) 
      A1=CCC(JV+JE) 
      T8=A1*B-B1*A 
      IF(T8.EQ.0.)GO TO 160 
      C1=CCC(JV+JE2) 
      XO=(C1*A-C*A1)/T8 
      IF(A.NE.0.)GO TO 90 
      YO=-C1-B1*XO 
      GO TO 100 
   90 CONTINUE 
      YO=-C-B*XO 
  100 CONTINUE 
      T=XO 
      IF(A.EQ.0.)T=YO 
      IF((T-XCC3)*(T-XCC4).GE.0.)GO TO 160 
      T=XO 
      IF(A1.EQ.0.)T=YO 
      S1=T-CCC(JV+JE4) 
      S=T-CCC(JV+JE3) 
      IF((ABS(S).LE.EXP).OR.(ABS(S1).LE.EXP))GO TO 110 
      IF(S*S1.GT.0.)GO TO 160 
  110 CONTINUE 
      IF(VX3.NE.0.)GO TO 130 
      TSZ=CL 
      TSX=AL 
      VT=XO-X1 
      IF(TSX.NE.0.)GO TO 120 
      VT=YO-Y1 
      TSX=BL 
  120 CONTINUE 
      ZX1=(TSZ/TSX)*VT+Z1 
      GO TO 140 
  130 CONTINUE 
      ZX1=-(VX+VX1*YO+VX2*XO)/VX3 
  140 CONTINUE 
      IF(ISILH(J).EQ.0)GO TO 148 
      IF(ISILH(J).EQ.ISILH(LE))GO TO 150 
  148 CONTINUE 
      ZX=-(AC*XO+BC*YO+D)*G 
      IF(ABS(ZX-ZX1).LT.EXP)GO TO 150 
      IF(ZX1.GT.ZX)GO TO 160 
  150 CONTINUE 
      M=M+1 
      NQ=M+1 

!     STORES INTERSECTIONS.

      XI(NQ)=XO 
      YI(NQ)=YO 
      ZI(NQ)=ZX1 
  160 END DO 
  170 END DO 
      NGX(1)=M 
  190 RETURN 
      END Subroutine Check   ! -------------------------------------------------

      SUBROUTINE SOLVE(IXR,J,XXX,CCC,II,NNO,NIT,                        &
     & X21,Y21,Z21,IIA,NC,ZM,ZMI,LZ)                                     

!     THIS SUBROUTINE SOLVES FOR THE LINES OF INTERSECTION RESULTING    
!     FROM THE INTERSECTIONS OF THE JTH ELEMENT WITH THE OTHER          
!     RELEVANT ELEMENTS.



      DIMENSION XXX(:),CCC(:),NNO(:),ZM(:)                                   
      DIMENSION ZMI(:),X21(:),Y21(:),Z21(:),IIA(:),IV(2)                    
      DIMENSION XA(50),YA(50),ZA(50) 
      COMMON/DAVE/XCC(800) 
      COMMON/GO3/L0,L1,L00,L01,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13 
!!!      COMMON/DRH/ISILH(1)
      INCLUDE 'drh.cmn' 
      COMMON/INDX/JT,JB 
      DATA EXX,ERS/.001,.015/ 
      C3=XXX(3+JT) 
      IF(C3.EQ.0.)GO TO 80 
      A3=XXX(1+JT) 
      B3=XXX(2+JT) 
      D3=XXX(4+JT) 
      ZQ=ZM(L2+J) 
      DO 70 L=1,II 
      K=NNO(L4+L) 

!     CHECKS TO SEE IF THIS RELEVANT ELEMENT IS TO BE CONSIDERED FOR    
!     INTERSECTION

      IF((K.LT.0).OR.(K.LT.J))GO TO 70 
      IF(ZQ.LT.ZMI(L3+K))GO TO 70 
      IF((ISILH(J).EQ.ISILH(K)).AND.(ISILH(J).NE.0))GO TO 70 
      JX=L12+(K-1)*5 
      MT=0 
      A4=XXX(1+JX) 
      B4=XXX(2+JX) 
      C4=XXX(3+JX) 
      D4=XXX(4+JX) 

!     DETERMINES THE EQUATION OF LINE OF INTERSECTION.                  

      B=A3*C4-A4*C3 
      A=B3*C4-B4*C3 
      C=D3*C4-D4*C3 
      IF((A.EQ.0.).AND.(B.EQ.0.))GO TO 70 
      IF(A.NE.0.)GO TO 10 
      A=0 
      C=C/B 
      B=1 
      GO TO 20 
   10 CONTINUE 
      B=B/A 
      C=C/A 
      A=1 
   20 CONTINUE 
      IV(1)=J 
      IV(2)=K 
      S3=2*A+2*B+C 
      DO 60 M=1,2 
      I=IV(M) 
      JJ=L13+(I-1)*LZ 
      IG=5+L12+(I-1)*5 
      NK=XXX(IG)*5 
      JJ1=JJ+1 
      JJ2=JJ+2 
      JJ3=JJ+3 
      JJ4=JJ+4 
      DO 50 JV=1,NK,5 
      A1=CCC(JV+JJ) 
      B1=CCC(JV+JJ1) 
      C1=CCC(JV+JJ2) 

!     CHECK TO BE SURE LINE OF INTERSECTION IS NOT BOUNDARY LINE        
!     OF THE JTH SET.                                                   

      S2=A1*2.+B1*2.+C1 
      IF(ABS(S2-S3).LT.EXX)GO TO 70 


!     DETERMINES THE POINTS OF INTERSECTIONS OF THE LINE OF INTERSECTION
!     WITH OTHER LINES OF RELEVANT ELEMENTS.                            


      T8=A1*B-B1*A 
      IF(ABS(T8).LE.ERS)GO TO 50 
      XO=(C1*A-C*A1)/T8 
      IF(A.NE.0.)GO TO 30 
      YO=-C1-B1*XO 
      GO TO 40 
   30 CONTINUE 
      YO=-C-B*XO 
   40 CONTINUE 
      T=XO 
      IF(A1.EQ.0.)T=YO 
      IF((T-CCC(JV+JJ4))*(T-CCC(JV+JJ3)).GT.0.)GO TO 50 
      MT=MT+1 

!     STORE THE PTS OF INTERSECTIONS. 

      XA(MT)=XO 
      YA(MT)=YO 
      ZA(MT)=-(D3+A3*XO+B3*YO)/C3 
      ZT=-(D4+A4*XO+B4*YO)/C4 
      IF(ABS(ZT-ZA(MT)).GT.EXX)GO TO 70 
   50 END DO 
   60 END DO 
      CALL STAT2(MT,NIT,IXR,X21,Y21,Z21,IIA,IV,A,B,C,J,XA,YA,ZA,CCC,XXX,NC,LZ)
   70 END DO
   80 CONTINUE 
      NR=5*XXX(5+JT) 
      DO 90 IS=1,NR 
      XCC(IS)=CCC(IS+JB) 
   90 END DO 
      XXX(5+JT)=XXX(5+JT)+NIT 
      RETURN 
      END Subroutine Solve

!+
      SUBROUTINE STAT2(MT,NIT,IXR,X21,Y21,Z21,IIA,IV,A,B,C,IK,XA,YA,ZA,CCC,XXX,NC,LZ)
! ------------------------------------------------------------------------------
! PURPOSE - This subroutine takes the pts of intersection determined by
!     subroutine Solve and picks the coordinates with the max and
!     min x coordinates provided they lie on the interior/boundary
!     of both elements.

      DIMENSION XXX(:),CCC(:),X21(:)
      DIMENSION Y21(:),Z21(:),IIA(:),IV(:),XA(:),YA(:),ZA(:)                      
      COMMON/DAVE/XCC(800) 
      COMMON/GO3/L0,L1,L00,L01,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13 
      DATA EXX/.015/ 

      NX=0 
      IF(NIT.GE.120)GO TO 160 
      IF(MT.EQ.0)GO TO 160 
      DO 50 JX=1,MT 
      EI=0 
   10 EI=EI+.1 
      IF(EI.GE..5)GO TO 160 
      D=EI*XA(JX)-YA(JX) 
      DO 40 JO=1,2 
      M=IV(JO)-1 
      JC=L13+(M)*LZ 
      JXC=L12+(M)*5 
      NK=XXX(5+JXC) 
      I=0 
      IB=NK*5 


!     DETERMINE IF THE PROJECTION OF THE POINT OF INTERSECTION          
!     BELONGS TO THE INTERIOR OF BOTH PLANES.                           


      DO 30 J=1,IB,5 
      VE=XA(JX) 
      J6=J+JC 
      IF(CCC(J6).EQ.0.)VE=YA(JX) 
      T=CCC(J6)*YA(JX)+CCC(J6+1)*XA(JX)+CCC(J6+2) 
      IF((ABS(T).LT.EXX).AND.((VE-CCC(J6+3))*(VE-CCC(J6+4)).LE.0.))     &
     &GO TO 40                                                          
      R=EI*CCC(J6)+CCC(J6+1) 
      IF(R.EQ.0.)GO TO 30 
      T=(-CCC(J6+2)+CCC(J6)*D)/R 
      IF(T.LT.XA(JX))GO TO 30 
      IF(CCC(J6).NE.0.)GO TO 20 
      T=EI*T-D 
   20 CONTINUE 
      IF((T.EQ.CCC(J6+3)).OR.(T.EQ.CCC(J6+4)))GO TO 10 
      IF((T-CCC(J6+3))*(T-CCC(J6+4)).GT.0.)GO TO 30 
      I=I+1 
   30 END DO 
      IF(I-(I/2)*2.EQ.0)GO TO 50 
   40 END DO 
      NX=NX+1 
      XA(NX)=XA(JX) 
      YA(NX)=YA(JX) 
      ZA(NX)=ZA(JX) 
   50 END DO 
      IF(NX.EQ.0)GO TO 160 


!     THIS CODE FINDS THE MAX/MIN X-COORDINATES(Y-COORDINATES) AND
!     STORES THEM. FUTHERMORE BOTH THE EQUATION OF LINE AND POINTS(2)
!     ARE TREATED LIKE ADDITIONAL EDGES. IN THIS WAY, THE ALGORITHM NEED
!     NOT BE DISTURBED. ESSENTIALLY,THEN,THIS TRICK IS TRANSPARENT TO
!     THE REST OF THE PROGRAM.


      AMAXX=-(10**6) 
      AMINX=-AMAXX 
      AMAXY=AMAXX 
      AMINY=AMINX 
      IS=5+(IK-1)*5+L12 
      IS=XXX(IS) 
      DO 110 JI=1,NX 
      IF(A.EQ.0.)GO TO 80 
      IF(XA(JI).GE.AMINX)GO TO 60 
      AMINX=XA(JI) 
      YI=YA(JI) 
      ZI=ZA(JI) 
   60 IF(XA(JI).LE.AMAXX)GO TO 70 
      AMAXX=XA(JI) 
      YII=YA(JI) 
      ZII=ZA(JI) 
   70 CONTINUE 
      GO TO 110 
   80 CONTINUE 
      IF(YA(JI).GE.AMINY)GO TO 90 
      AMINY=YA(JI) 
      XI=XA(JI) 
      ZI=ZA(JI) 
   90 CONTINUE 
      IF(YA(JI).LE.AMAXY)GO TO 100 
      XII=XA(JI) 
      AMAXY=YA(JI) 
      ZII=ZA(JI) 
  100 CONTINUE 
  110 END DO 
      NIT=NIT+1 
      K=5*(NIT-1+IS)+1 
      XCC(K)=A 
      XCC(K+1)=B 
      XCC(K+2)=C 
      IF(A.EQ.0.)GO TO 120 
      XCC(K+3)=AMINX 
      XCC(K+4)=AMAXX 
      AMIN=AMINX 
      AMAX=AMAXX 
      YE=YII 
      ZE=ZII 
      GO TO 130 
  120 CONTINUE 
      XCC(K+3)=AMINY 
      XCC(K+4)=AMAXY 
      AMIN=XI 
      AMAX=XII 
      YI=AMINY 
      YE=AMAXY 
      ZE=ZII 
  130 CONTINUE 
      IG=IXR+NIT*3 
      I8=IG-2 
      X21(I8)=AMIN 
      Y21(I8)=YI 
      Z21(I8)=ZI 
      DO 140 JK=1,2 
      IE=IG-JK+1 
      X21(IE)=AMAX 
      Y21(IE)=YE 
      Z21(IE)=ZE 
  140 END DO 
      DO 150 JK=1,2 
      IIA(IG-JK)=0 
  150 END DO 
      IIA(IG)=1 
      TX=(AMAX-AMIN)**2 
      TY=(YE-YI)**2 
      DX=(TX+TY)**.5 
      IF(DX.LT..1)NIT=NIT-1 
  160 CONTINUE 
      IF(EI.GE..5)NIT=0 
      RETURN 
      END Subroutine Stat2
                                                                        
!   The following routines VSRTR and VSRT1 are supplied as substitutes
!   for the routines of the same name in the NASA distribution from 
!   COSMIC labelled ARC-12721.                                                 
!   The NASA distribution used routines that were copyrighted by IMSL,In
!   and were distributed with the kind permission of IMSL. The IMSL     
!   routines are possibly more robust and efficient than these equivalen
!   ones and if you have access to the IMSL libraries, you may want to  
!   put them back.                                                      

!   The standard practice of Public Domain Aeronautical Software is to  
!   distribute only code that is truly free of legal restrictions.      

!+                                                                      
      SUBROUTINE VSRTR(a, n, ix) 
!   --------------------------------------------------------------------
!     PURPOSE - Sort and index a real array                             

!     AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software 
!       with great assistance from the authors of "Numerical Recipes,   
!       The Art of Scientific Computing".                               

!     NOTES - This is a replacement for the routine of the same name    
!       in the IMSL library. This routine is not copyrighted. The       
!       routine uses the heapsort algorithm. The coding follows the     
!       examples of the routines SORT2 and INDEXX in "Numerical Recipes"

!       On input, a contains the vector of length n to be sorted.       
!       On input, ix contains an identifying integer for each element   
!       of a. Usually, one sets ix=1,2,3,4,....                         
!       For real generality, the values of ix may be anything that is   
!       useful and meaningful to the programmer. Whenever two values of 
!       the array a are swapped, the corresponding entries in the array 
!       ix are also swapped.                                            
!       After execution, the a-array will be in increasing order and the
!       array ix will hold a history of the transformations             

!       The hidden-line program (subroutine LIN) is an example of the   
!       use of an array other than 1,2,3,4,... for ix.                  

      IMPLICIT NONE 
!***********************************************************************
!     A R G U M E N T S                                                *
!***********************************************************************
      INTEGER,INTENT(IN):: n 
      REAL,INTENT(IN OUT),DIMENSION(:):: a
      INTEGER,INTENT(IN OUT),DIMENSION(:):: ix
!***********************************************************************
!     L O C A L   V A R I A B L E S                                    *
!***********************************************************************
      INTEGER:: i,j,k, ir, ixsave 
      REAL:: asave 
!-----------------------------------------------------------------------
      IF (n .LE. 1) RETURN 

      k=(n/2)+1 
      ir=n 

   20 CONTINUE 
      IF (k .GT. 1) THEN 
        k=k-1 
        asave=a(k) 
        ixsave=ix(k) 
      ELSE 
        asave=a(ir) 
        ixsave=ix(ir) 
        a(ir)=a(1) 
        ix(ir)=ix(1) 
        ir=ir-1 
        IF (ir .EQ. 1) THEN 
          a(1)=asave 
          ix(1) = ixsave 
          RETURN   ! this is the only way out of here          
        END IF 
      END IF 

      i = k 
      j = k+k 
   30 CONTINUE 
      IF ( j .LE. ir) THEN 
        IF (j .LT. ir) THEN 
          IF ( a(j) .LT. a(j+1) ) j=j+1 
        END IF 

        IF ( asave .LT. a(j) ) THEN 
          a(i)=a(j) 
          ix(i) = ix(j) 
          i = j 
          j = j+j 
        ELSE 
          j = ir+1 
        END IF 
        GO TO 30 
      END IF 

      a(i)=asave 
      ix(i) = ixsave 
      GOTO 20           
      END Subroutine Vsrtr   ! ----------------------- End of Subroutine SortAndIndexReal

!+                                                                      
      SUBROUTINE VSRT1(a, n, ix) 
!   --------------------------------------------------------------------
!     PURPOSE - Sort and index an integer array                         

!     AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software 
!       with great assistance from the authors of "Numerical Recipes,   
!       The Art of Scientific Computing".                               

!     NOTES - This is a replacement for the routine of the same name    
!       in the IMSL library. This routine is not copyrighted. The       
!       routine uses the heapsort algorithm. The coding follows the     
!       examples of the routines SORT2 and INDEXX in "Numerical Recipes"

!       On input, a contains the vector of length n to be sorted.       
!       On input, ix contains an identifying integer for each element   
!       of a. Usually, one sets ix=1,2,3,4,....                         
!       For real generality, the values of ix may be anything that is   
!       useful and meaningful to the programmer. Whenever two values of 
!       the array a are swapped, the corresponding entries in the array 
!       ix are also swapped.                                            
!       After execution, the a-array will be in increasing order and the
!       array ix will hold a history of the transformations             

      IMPLICIT NONE 
!***********************************************************************
!     A R G U M E N T S                                                *
!***********************************************************************
      INTEGER,INTENT(IN):: n 
                    
      INTEGER,INTENT(IN OUT),DIMENSION(:):: a  ! just like VSRTR but a is integer                 
      INTEGER,INTENT(IN OUT),DIMENSION(:):: ix
!***********************************************************************
!     L O C A L   V A R I A B L E S                                    *
!***********************************************************************
      INTEGER i,j,k, ir, ixsave 
      INTEGER asave 
!-----------------------------------------------------------------------
      IF (n .LE. 1) RETURN 

      k=(n/2)+1 
      ir=n 

   20 CONTINUE 
      IF (k .GT. 1) THEN 
        k=k-1 
        asave=a(k) 
        ixsave=ix(k) 
      ELSE 
        asave=a(ir) 
        ixsave=ix(ir) 
        a(ir)=a(1) 
        ix(ir)=ix(1) 
        ir=ir-1 
        IF (ir .EQ. 1) THEN 
          a(1)=asave 
          ix(1) = ixsave 
          RETURN    ! this is the only way out of here          
        END IF 
      END IF 

      i = k 
      j = k+k 
   30 CONTINUE 
      IF ( j .LE. ir) THEN 
        IF (j .LT. ir) THEN 
          IF ( a(j) .LT. a(j+1) ) j=j+1 
        ENDIF 

        IF ( asave .LT. a(j) ) THEN 
          a(i)=a(j) 
          ix(i) = ix(j) 
          i = j 
          j = j+j 
        ELSE 
          j = ir+1 
        ENDIF 
        GOTO 30 
      ENDIF 

      a(i)=asave 
      ix(i) = ixsave 
      GOTO 20 
      END Subroutine Vsrt1  ! -------------------- End of Subroutine SortAndIndexInteger                                           


!+
      SUBROUTINE BigWindow(vminx,vmaxx,vminy,vmaxy,                     &  
     &    xmin,xmax, ymin,ymax, xmin2,xmax2,ymin2,ymax2)
! ----------------------------------------------------------------------
! PURPOSE - Define a window that contains the data in the rectangle
!    [xmin,xmax]x[ymin,ymax] and has the same aspect ratio as the
!    viewport [vminx,vmaxx]x[vminy,vmaxy]
      IMPLICIT NONE

! the four quantities in each line below are left,right,bottom,top
      REAL  vminx,vmaxx, vminy,vmaxy  ! viewport
      REAL  xmin,xmax, ymin,ymax      ! data window
      REAL  xmin2,xmax2,ymin2,ymax2   ! resulting window

      REAL dx              ! window size in x-direction if y is dominant
      REAL dy              ! window size in y-direction if x is dominant
!-----------------------------------------------------------------------
      IF ((xmax-xmin)*(vmaxy-vminy) .GT. (ymax-ymin)*(vmaxx-vminx)) THEN
        xmin2=xmin
        xmax2=xmax
        dy=(xmax-xmin)*(vmaxy-vminy)/(vmaxx-vminx)
        ymin2=0.5*(ymax+ymin-dy)
        ymax2=0.5*(ymax+ymin+dy)
      ELSE
        dx=(ymax-ymin)*(vmaxx-vminx)/(vmaxy-vminy)
        xmin2=0.5*(xmax+xmin-dx)
        xmax2=0.5*(xmax+xmin+dx)
        ymin2=ymin
        ymax2=ymax
      END IF

      RETURN
      END Subroutine BigWindow   ! =================================

      SUBROUTINE ScanGnu(efu, xmin,xmax, ymin,ymax)
      IMPLICIT NONE
      INTEGER efu

      INTEGER code
      REAL x,y
      REAL xmin,xmax, ymin,ymax
      CHARACTER dummy*80

      xmin=1E37
      xmax=-xmin
      ymin=xmin
      ymax=xmax

      REWIND(efu)
      DO
        READ(efu,'(A)', IOSTAT=code) dummy
        IF (code .LT. 0) EXIT
        IF (dummy(1:1) .EQ. '#') CYCLE
        IF (dummy .EQ. ' ') CYCLE

        READ(dummy,*,IOSTAT=code) x,y
        IF (code .NE. 0) CYCLE

        xmax=MAX(xmax,x)
        xmin=MIN(xmin,x)
        ymax=MAX(ymax,y)
        ymin=MIN(ymin,y)
      END DO

      RETURN
      END Subroutine ScanGnu   ! ====================================

      SUBROUTINE CreateFig(gnu,fig)
      IMPLICIT NONE
      INTEGER gnu,fig

      CHARACTER dummy*80
      CHARACTER*80 data(1000)
      INTEGER i,k ,code


      Write(fig,*) 'LINEWIDTH 0.5'

      k=0
      DO
        Read(gnu,'(A)',IOSTAT=code) dummy
        IF (code .LT. 0) EXIT
        IF (dummy(1:1) .EQ. '#') CYCLE
        IF (dummy .EQ. ' ') THEN
          IF (k .GT. 1) THEN
            Write(fig, '(A,I5)')  'DATA ', k
            WRITE(fig,'(A)') (Trim(data(i)),i=1,k)
          END IF
          k=0
        ELSE
          k=k+1
          data(k) = dummy
        END IF
      END DO

      Write(fig,'(A)') 'ENDFRAME ------------------------------------'

      RETURN
      END  Subroutine CreateFig   ! =================================



      SUBROUTINE CreatePs(gnu,ps)
      IMPLICIT NONE
      INTEGER efu
      INTEGER,PARAMETER:: FIG=7

      REAL:: xmin,xmax, ymin,ymax
!      REAL:: xmin2,xmax2, ymin2,ymax2
      REAL:: xlow,xhigh,ylow,yhigh
      REAL,PARAMETER:: VMINX=0.01, VMAXX=0.99, VMINY=0.03, VMAXY=0.73

      INTEGER gnu,ps
      REAL MARGINPS,SCALEPS, MARGINPCL,SCALEPCL
      PARAMETER (MARGINPS=36.0, SCALEPS=720.0)
      PARAMETER (MARGINPCL=0.0, SCALEPCL=10160.0)

      CHARACTER buffer*132
      INTEGER code
      INTEGER i,n
      REAL x,y, sx,sy
      REAL sclx,scly
      INTEGER ix,iy

      REWIND(UNIT=gnu)
      CALL ScanGnu(gnu, xmin,xmax, ymin,ymax)
      REWIND(UNIT=gnu)
      CALL BigWindow(vminx,vmaxx,vminy,vmaxy,                           &  
     &    xmin,xmax, ymin,ymax, xlow,xhigh, ylow,yhigh)
      sclx=(vmaxx - vminx) / (xhigh - xlow)
      scly=(vmaxy - vminy) / (yhigh - ylow)

      OPEN(UNIT=FIG,FILE='hlp.fig',STATUS='REPLACE',ACTION='WRITE')
      CALL CreateFig(gnu,FIG)
      CLOSE(UNIT=GNU)
      CLOSE(UNIT=FIG)
      OPEN(UNIT=FIG,FILE='hlp.fig',STATUS='OLD',ACTION='READ')

      WRITE(PS,*) '%!PS-Adobe'
      WRITE(PS,*) '612 0 translate 90 rotate'
      WRITE(PS,*) '0.5 setlinewidth'

      DO
        Read(fig, '(A)', IOSTAT=code) buffer
        IF (code .LT. 0) EXIT

        IF (buffer(1:4) .EQ. 'DATA') THEN
          READ(buffer(6:10), '(I5)' ) n
          Read(fig,*) x, y
          sx=VMINX + (x-xlow)*sclx
          sy=VMINY + (y-ylow)*scly
          sx=MARGINPS +SCALEPS*sx
          sy=MARGINPS +SCALEPS*sy
          WRITE(Ps,'(2F9.1,A)' )  sx, sy, ' moveto'
          DO i=2,n
            Read(fig,*) x, y
            sx=VMINX + (x-xlow)*sclx
            sy=VMINY + (y-ylow)*scly
            sx=MARGINPS +SCALEPS*sx
            sy=MARGINPS +SCALEPS*sy
            Write(Ps, '(2F9.1,A)' ) sx, sy,  ' lineto'
          END DO
        ELSE IF (buffer(1:8) .EQ. 'ENDFRAME') THEN
          Write(Ps,*) 'stroke showpage'
        END IF
      END DO

      CLOSE(UNIT=gnu)
      Close(UNIT=PS)

      RETURN
      END Subroutine CreatePs   ! ----------------------------------------------



END Module HiddenLineProcedures   ! ============================================

!+                                                                      
      PROGRAM HiddenLineProgram 
                                                       !    \hlp\hlp.f90
!   --------------------------------------------------------------------
!     PURPOSE - Plot a perspective view with hidden lines removed of an 
!        object defined in WGS format                                   
!                                                                       
!     AUTHOR - David Hedgley, NASA/AMES (Dryden)                        
!              Ralph L. Carmichael, Public Domain Aeronautical Software 
!                                                                       
!     REVISION HISTORY                                                  
!   DATE  VERS PERSON  STATEMENT OF CHANGES                             
!   1980   xxx   DRH   Original SKETCH and HIDLIN programs              
!   1982   xxx   RLC   Old PANSKETCH program. Input in PANAIR format    
!   1988   xxx   DRH   Upgrade to Silhouette                            
! 13Dec88  0.1   RLC   Original adaptation from PANSKETCH. WGS format   
! 30Oct92  0.2   RLC   Extensive rewrite; incorporated Silhouette       
! 14Dec93  0.3   RLC   Raised MAXGRID to 100000 (later removed)         
! 16Dec94  0.4   RLC   Adapted to MS-DOS environment                    
! 20Dec94  0.5   RLC   Use external file instead of virtual memory      
! 17Aug95  0.6   RLC   Changed output to be usable by gnuplot           
! 23Aug95  0.7   RLC   Replaced IMSL routines VSRTR and VSRT1 with PD on
!  5Dec95  1.0   RLC   Took the E1 out of Format in Plot                
! 31Dec96  2.0   RLC   Protect against inability to open output file
!  5Jan99  2.1   RLC   Raised MAXPANELS=25000. Called it bighlp.
! 20Dec04  2.2   RLC   Set a value for I1 in LIN. Otherwise it is undefined.
!                      Open PLT with 'REPLACE' instead of 'UNKNOWN'
!                      Increased MAXPANELS to 50000, MAXROWCOL=1000
! 09Feb09  2.3   RLC   Changed name of subroutine Stat to Stat2 to avoid
!                         conflict with new Fortran intrinsic
!                                                                       
!     NOTES- The common blocks /GO/ and /DRH/ must have the same length 
!       consistently throughout the program and the length is related   
!       to the value of MAXPANELS. /DRH/ must have a length of MAXPANELS
!       reals and /GO/ must have a length of 69*MAXPANELS reals.        
!       Since /DRH/ appears several times, it is INCLUDEd from an       
!       external file called drh.cmn. The block /GO/ only appears       
!       twice (in the main program and in LIN) and must be changed      
!       whenever MAXPANELS is modified.                                 
!                                                                       
USE HiddenLineProcedures
IMPLICIT NONE 
                                                                       
!     BUG LIST -                                                        
!       * WGS transformation parameters are not processed      
!       * if you compile with bounds checking turned ON, you always 
!         abort at run time due to an array being addressed outside 
!         its bounds. It is in LIN on the arrays that are
!         equivalenced to WORK
!***********************************************************************
!     C O N S T A N T S                                                *
!***********************************************************************
      INTEGER,PARAMETER:: WGS=1, GNU=2, PS=4   ! units
      CHARACTER AUTHOR*63, MODIFIER*63 
      CHARACTER VERSION*30, FAREWELL*60 
  CHARACTER(LEN=*),PARAMETER:: GREETING='hlp - draw configuration views'
      PARAMETER (AUTHOR=                                                &
     & 'Ralph L. Carmichael, Public Domain Aeronautical Software')      
      PARAMETER (MODIFIER=' ')  ! put your name here if you modify 
      PARAMETER (VERSION= '2.3 (9 February 2009)' ) 
      PARAMETER (FAREWELL= 'File hlp.gnu added to your directory.') 
                       
      INTEGER MAXROWCOL  ! max gridpoints in ONE network                        
!      INTEGER MAXPANELS  ! max number of panels in configuration(&image)
!      PARAMETER (MAXROWCOL=1000, MAXPANELS=50000)
      PARAMETER (MAXROWCOL=1000)
!***********************************************************************
!     V A R I A B L E S                                                *
!***********************************************************************
      INTEGER errCode 
      CHARACTER*80 fileName 
      INTEGER i,j,k 
      INTEGER irow 
      INTEGER jcol 
      INTEGER nc 
      INTEGER nets 
      CHARACTER*80 netTitle 
      INTEGER np 
      INTEGER npanels 
      CHARACTER*1 yorn 
      LOGICAL yImages 
!***********************************************************************
!     A R R A Y S                                                      *
!***********************************************************************
      REAL x(5),y(5),z(5) 
      REAL xnet(MAXROWCOL),ynet(MAXROWCOL),znet(MAXROWCOL) 
!***********************************************************************
!     C O M M O N   B L O C K S                                        *
!***********************************************************************
!!!      INTEGER isilh(MAXPANELS)                                       
!!!      COMMON /DRH/ isilh                                             
      INCLUDE 'drh.cmn' 

!      REAL work(69*MAXPANELS) 
!      COMMON /GO/ work 

      INTEGER ICORE 
      INTEGER MNE 
      INTEGER MNP 
      REAL dv,scf 
      REAL roll,pitch,yaw 
      COMMON /SCALAR/ scf, yaw, roll, pitch, MNE, dv, MNP, icore 
!-----------------------------------------------------------------------
      MNE=4 
      MNP=6 
                                          
      ICORE=(25 + 5*MNE + 4*MNP)*MAXPANELS ! =69*MAXPANELS for this case 
      DO i=1,MAXPANELS 
        isilh(i)=0 
      END DO 
!                                                                       
      WRITE(*,*) GREETING 
      WRITE(*,*) AUTHOR 
      WRITE(*,*) 'Modified by '//MODIFIER 
      WRITE(*,*) 'Version '//VERSION 
      WRITE(*,*) 

      DO                                                                       
        WRITE(*,*) 'Enter the name of the wgs input file: ' 
        READ(*, '(A)' ) FILENAME 
        IF (FILENAME .EQ. ' ') STOP 
        OPEN (UNIT=WGS, FILE=FILENAME, STATUS='OLD', IOSTAT=errCode) 
        IF (errCode==0) EXIT
        OPEN (UNIT=WGS, FILE=Trim(FILENAME)//'.wgs', &
          STATUS='OLD', IOSTAT=errCode) 
        IF (errCode==0) EXIT
        WRITE (*, '(A)' ) ' Unable to open this file. Try again...' 
      END DO
 
!..... Determine which options are to be used...........................
      WRITE (*, '(A)' ) ' If only the right half of the vehicle is',    &
     & ' defined on the WGS file, this program may be directed to',     &
     & ' create the image panels in order to see the complete vehicle.',&
     & ' Do you want to make y-images of each network (Y/N) ? '         
      READ (*, '(A)' ) yorn 
      yImages=yorn.EQ.'Y' .OR. yorn.EQ.'y' 
      IF (yImages) THEN 
         WRITE (*,*) ' y-images will be drawn.' 
      ELSE 
         WRITE (*,*) ' No image networks will be drawn.' 
      END IF 

    9 CONTINUE                ! come back here if more views are made

!..... Determine the viewing position...................................
      WRITE(*,*) ' Enter the first rotation angle (degrees): ' 
      READ(*,*) ROLL 
      WRITE(*,*) ' Enter the second rotation angle (degrees): ' 
      READ(*,*) PITCH 
      WRITE(*,*) ' Enter the third rotation angle (degrees): ' 
      READ(*,*) YAW 
      WRITE(*,*) ' Enter distance from observer to origin: ' 
      READ(*,*) DV 
                       
      SCF=1.0     ! i have never used anything but 1


!..... Open the plot file. This example may be used for gnuplot.        
      OPEN(UNIT=GNU, FILE='hlp.gnu', &
        STATUS='REPLACE', IOSTAT=errCode, ACTION='WRITE') 
      IF (errCode .NE. 0) THEN
        WRITE(*,*) 'Unable to open plot file' 
        STOP
      END IF

!!!      WRITE(PLT,*) '#Created by hlp'                                 
!!!      WRITE(PLT,*) '#FILE '//fileName                                
!!!      WRITE(PLT,'(A,F12.1)') '#ROLL ' , roll                         
!!!      WRITE(PLT,'(A,F12.1)') '#PITCH ', pitch                        
!!!      WRITE(PLT,'(A,F12.1)') '#YAW '  , yaw                          
!!!      WRITE(PLT,'(A,F15.1)') '#DV '   , dv                           
!!!..... gnuplot comments have # in col 1                               
!                                                                       
!                                                                       
!.....Read the networks one by one, calling Sketch for each one (and    
!.       and its image if yImages is .TRUE.). Then call Sketch again to   
!.       compute the picture.                                           
      REWIND WGS 
      WRITE(*,*) ' Reading the network data...' 
      nets=0 
      npanels=0 

      nc=0              ! this is a flag for Sketch
      np=5 
                                            
      READ (WGS, '(A)') netTitle     ! read the general title

   10 CONTINUE    ! keep coming back to read more nets until EOF
   
      READ(WGS, '(A)', END=900) netTitle 
      write(*,*) ' reading '//netTitle(1:69) 

      READ(WGS, *, END=900) j, jcol, irow   ! ignore transformation
      IF (jcol*irow .EQ. 0) THEN 
         WRITE(*,*) ' Zero panels in this network' 
         GOTO 10 
      END IF 
      IF (jcol*irow .GT. MAXROWCOL) THEN
        WRITE(*,*) 'Too many gridpoints in this net.'
        STOP
      END IF
      READ (WGS, *, END=900) (xnet(i),ynet(i),znet(i),i=1,irow*jcol) 

      DO j=2,jcol 
        DO i=2,irow 
          k = i-1 + (j-2)*irow 
          x(1)=xnet(k) 
          y(1)=ynet(k) 
          z(1)=znet(k) 
          k = i-1 + (j-1)*irow 
          x(2)=xnet(k) 
          y(2)=ynet(k) 
          z(2)=znet(k) 
          k = i + (j-1)*irow 
          x(3)=xnet(k) 
          y(3)=ynet(k) 
          z(3)=znet(k) 
          k = i + (j-2)*irow 
          x(4)=xnet(k) 
          y(4)=ynet(k) 
          z(4)=znet(k) 
          x(5)=x(1) 
          y(5)=y(1) 
          z(5)=z(1) 

          npanels=npanels+1 
          IF (npanels .GT. MAXPANELS) STOP 'Too many panels. ABORT!!' 

          CALL Sketch(X, Y, Z, NP, NC)  ! send one panel to sketch
          IF (yImages) THEN 
            DO k=1,5 
              y(k)=-y(k) 
            END DO 
            npanels=npanels+1 
            IF (npanels .GT. MAXPANELS) STOP 'Too many panels' 
            CALL Sketch(x, y, z, np, nc) 
          END IF 
        END DO 
      END DO 
      nets=nets+1 
      GOTO 10    ! go back and look for another network

!..... After all networks have been sent to SKETCH, call it again.......
!      This time with nc=1 

  900 CONTINUE 
      WRITE(*,*) nets,    ' networks with a total of' 
      WRITE(*,*) npanels, ' panels read from the input file.' 
      WRITE(*,*) ' Computing the scene...' 
      nc=1 
      CALL Sketch(x, y, z, np, nc) 


!..... See if there are any more views to make..........................
!!!!!      WRITE (*, '(A)') '$Any more views? '
!!!!!      READ  (*,'(A)') yorn
!!!!!      IF (yorn.EQ.'Y' .OR. yorn.EQ.'y') GOTO 9

      CLOSE(UNIT=WGS) 
      CLOSE(UNIT=GNU) 

! call PostScript writer here...
      OPEN(UNIT=GNU, FILE='hlp.gnu',STATUS='OLD',ACTION='READ')
      OPEN(UNIT=PS, FILE='hlp.ps',STATUS='REPLACE',ACTION='WRITE')
      WRITE(*,*) 'Scene computed. Generating PostScript file.'
      CALL CreatePs(GNU,PS)

      WRITE (*,*) FAREWELL 
      WRITE(*,*) 'hlp has terminated successfully.' 
      STOP
      END Program HiddenLineProgram   ! ----------------------------------------

