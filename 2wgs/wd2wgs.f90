!+
! PROGRAM WaveDragToWgs 
!   --------------------------------------------------------------------
! PURPOSE - Convert a configuration defined in WaveDrag format to
!   Langley Wireframe Geometry Standard (LaWGS) format.

! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software

!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 21Jun95  0.1   RLC   Original coding
! 19Sep95  0.2   RLC   Added MODIFIER. errCode
! 13Dec95  1.0   RLC   Production version
! 22Dec96  2.0   RLC   Protect against failure to open output file
! 31Dec01  2.1   RLC   F90
! 22Jul08  2.2   RLC   Moved procedures into module; file name from command line
! 02Dec08  2.3   RLC   GET_COMMAND_ARGUMENT for reading command line
!-------------------------------------------------------------------------------

!+
MODULE Wd2WgsCommon
! ------------------------------------------------------------------------------
! PURPOSE - Replace the common blocks in original code

IMPLICIT NONE
  CHARACTER(LEN=*),PARAMETER:: VERSION = '2.3 (2 December 2008)' 

  INTEGER:: j0,j1,j2,j3,j4,j5,j6
  INTEGER:: nwaf,nwafor,nfus,nradx(4),nforx(4),np,npodor 
  INTEGER:: nf,nfinor,ncan,ncanor

  REAL,DIMENSION(30):: xaf
  REAL,DIMENSION(20,4):: waforg,w
  REAL,DIMENSION(20,3,30):: waford
  REAL,DIMENSION(20,30):: tzord
  REAL,DIMENSION(20,2):: ordmax
!      COMMON/WING/xaf,waforg,w,waford,tzord,ordmax

  REAL,DIMENSION(30,4):: xfus,zfus,fusard,fusrad
  REAL,DIMENSION(30,30,8):: sfus
!      COMMON/FUSE/xfus,zfus,fusard,fusrad,sfus 


  REAL,DIMENSION(9,3):: podorg
  REAL,DIMENSION(9,30):: xpod,podord
!      COMMON /POD/ podorg,xpod,podord

  REAL,DIMENSION(6,2,4):: finorg
  REAL,DIMENSION(6,10):: xfin
  REAL,DIMENSION(6,2,10):: finord, finx2,finx3 
  REAL,DIMENSION(6):: finmx1,finmx2, finth1,finth2 
!      COMMON /FIN/ finorg,xfin,finord,                            &
!     & finx2,finx3,finmx1,finmx2,finth1,finth2

  REAL,DIMENSION(2,2,4):: canorg
  REAL,DIMENSION(2,10):: xcan
  REAL,DIMENSION(2,2,10):: canord,canor1,canorx
  REAL,DIMENSION(2,2,2):: canmax
!      COMMON /CAN/ canorg,xcan,canord,canor1,canorx,canmax

  REAL,PARAMETER:: PI=3.14159265
END Module Wd2WgsCommon   ! ====================================================

!+
MODULE WaveDragToWgsProcedures
! PURPOSE -

USE Wd2WgsCommon

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
SUBROUTINE ReadConfigurationData(lun) 
! --------------------------------------------------------------------
! PURPOSE - Read the configuration data.                            
!                                                                       
! NOTES -Data in /PARAMS/ has already been read.                    
!        The info read here goes into /WING/, /FUSE/, /PODS/,           
!        /FINS/, and /CAN/                                              
!                                                                       
IMPLICIT NONE

  INTEGER,INTENT(IN):: lun       ! logical unit number of the input file

  INTEGER:: i,j,k
  INTEGER:: ln 
  INTEGER:: n,nfu,nn 

!      REAL podorg(9,3),xpod(9,30),podord(9,30)
!      COMMON /POD/ podorg,xpod,podord 
!                                                                       
!-----------------------------------------------------------------------
  WRITE(*,*) 'Begin subroutine ReadConfigurationData' 

   20 FORMAT (10F7.0)            ! read everything here with this format

  IF (j0 == 1) READ(lun,20)  ! read, ignore reference area        

  IF (j1 .NE. 0) THEN                                   ! START WING
    READ(lun,20) (xaf(i),i=1,ABS(nwafor)) 
    DO i=1,nwaf 
      READ(lun,20) (waforg(i,j),j=1,4) 
    END DO 

    IF (j1 .LT. 0) THEN                    ! read or zero the camber
      DO i=1,nwaf 
        DO k=1,ABS(nwafor) 
          tzord(i,k)=0.0 
        END DO 
      END DO 
    ELSE 
      DO nn=1,nwaf 
        READ(lun,20) (tzord(nn,i),i=1,ABS(nwafor)) 
      END DO 
    END IF 

                    
    DO nn=1,nwaf    ! ordinates (z) 
      READ(lun,20) (waford(nn,1,i),i=1,ABS(nwafor)) 
      IF (nwafor.GE.0) THEN 
        DO k=1,ABS(nwafor) 
          waford(nn,2,k)=waford(nn,1,k) 
        END DO 
      ELSE 
        READ(lun,20) (waford(nn,2,i),i=1,ABS(nwafor)) 
      END IF 
    END DO                                                    ! END WING

  END IF 

  IF (j2.NE.0) THEN                                 ! START FUSELAGE
    DO nfu=1,nfus 
      n=nforx(nfu)    ! number of defining stations                  
      READ(lun,20) (xfus(i,nfu),i=1,n) 
      IF (j6.EQ.0 .AND. j2.EQ.-1) THEN       ! cambered fuselage    
        READ(lun,20) (zfus(i,nfu),i=1,nforx(nfu)) 
      ELSE 
!        DO i=1,nforx(nfu) 
!          zfus(i,nfu)=0.0 
!        END DO
        zfus(:,nfu)=0.0
      END IF 

      IF (j2.EQ.1) THEN 
        DO ln=1,nforx(nfu) 
          READ(lun,20) (sfus(i,ln,nfu+nfu-1),i=1,nradx(nfu)) 
          READ(lun,20) (sfus(i,ln,nfu+nfu),  i=1,nradx(nfu)) 
        END DO 
      ELSE 
        READ(lun,20) fusard(1:n,nfu)
        fusrad(1:n,nfu)=SQRT(fusard(1:n,nfu)/PI)
!        DO i=1,nforx(nfu) 
!          fusrad(i,nfu)=SQRT(fusard(i,nfu)/PI) 
!        END DO 
      END IF 
    END DO                   
  END IF                                              ! END FUSELAGE

                             
  IF (j3.NE.0) THEN                                     ! START PODS
    DO nn=1,np   ! all input is dimensional                         
      READ(lun,20) (podorg(nn,i),i=1,3) 
      READ(lun,20) (xpod(nn,i),i=1,npodor) 
      READ(lun,20) (podord(nn,i),i=1,npodor) 
    END DO 
  END IF                                                  ! END PODS

  IF (j4.NE.0) THEN                                     ! START FINS
    DO nn=1,ABS(nf) 
      READ(lun,20) ((finorg(nn,i,j),j=1,4),i=1,2) 
      READ(lun,20) (xfin(nn,i),i=1,nfinor) 
      READ(lun,20) (finord(nn,1,j),j=1,nfinor) 
      IF (nf.LT.0) THEN 
        READ(lun,20) (finord(nn,2,j), j=1,nfinor) 
      ELSE
        finord(nn,2,1:nfinor)=finord(nn,1,1:nfinor)
!        DO j=1,nfinor 
!          finord(nn,2,j)=finord(nn,1,j) 
!        END DO 
      END IF 
    END DO 
  END IF                                                  ! END FINS

  IF (j5.NE.0) THEN                                  ! START CANARDS
    n=ABS(ncanor) 
    DO nn=1,ncan 
      READ(lun,20) ((canorg(nn,i,j),j=1,4),i=1,2) 
      READ(lun,20) (xcan(nn,i),i=1,ABS(ncanor)) 
      READ(lun,20) (canord(nn,1,j),j=1,ABS(ncanor)) 
      canord(nn,2,1:n)=canord(nn,1,1:n)
!      DO j=1,ABS(ncanor)
!        canord(nn,2,j)=canord(nn,1,j) 
!      END DO 
      IF (ncanor >= 0) THEN 
        DO j=1,ncanor 
          canor1(nn,1,j)=canord(nn,1,j) 
          canor1(nn,2,j)=canord(nn,2,j) 
        END DO 
      ELSE 
        READ(lun,20) (canor1(nn,1,j),j=1,ABS(ncanor)) 
        DO j=1,ABS(ncanor) 
          canor1(nn,2,j)=canor1(nn,1,j) 
        END DO 
      END IF 
    END DO 
  END IF                                            ! END CANARDS

  WRITE (*,*) 'Configuration Data has been read'
  RETURN
END Subroutine ReadConfigurationData   ! ------------------------------------

!+                                                                      
SUBROUTINE TransformConfigurationData()
!   --------------------------------------------------------------------
! PURPOSE - Convert non-dimensional input data to actual units      
!          Transform wing coordinates from pct-chord to actual units    
!          of length, referred to common origin of problem. compute     
!          maximum ordinate of each airfoil.                            
!                                                                       
!     NOTES-                                                            
!                                                                       
      IMPLICIT NONE 
!***********************************************************************
!     L O C A L   V A R I A B L E S                                    *
!***********************************************************************
      REAL:: e,e2,e3,ee 
      INTEGER:: i,j,k 
      INTEGER:: lq 
      INTEGER:: nn 
                                                                       
!      REAL xaf(30),waforg(20,4),w(20,4) 
!      REAL waford(20,3,30),tzord(20,30),ordmax(20,2) 
!      COMMON/WING/xaf,waforg,w,waford,tzord,ordmax 
!                                                                       
!      REAL finorg(6,2,4),xfin(6,10) 
!      REAL finord(6,2,10),finx2(6,2,10),finx3(6,2,10) 
!      REAL finmx1(6),finmx2(6),finth1(6),finth2(6) 
!      COMMON /FIN/ finorg,xfin,finord,                                  &
!     & finx2,finx3,finmx1,finmx2,finth1,finth2                          
!                                                                       
!      REAL canorg(2,2,4),xcan(2,10),canord(2,2,10) 
!      REAL canor1(2,2,10),canorx(2,2,10),canmax(2,2,2) 
!      COMMON /CAN/ canorg,xcan,canord,canor1,canorx,canmax 
!-----------------------------------------------------------------------
  WRITE(*,*) 'Begin subroutine TransformConfigurationData' 

  IF (j1 .NE. 0) THEN                                   ! START WING
    DO I=1,nwaf 
      e=0.01*waforg(i,4) 
      DO j=1,ABS(nwafor) 
        waford(i,1,j)=  e*waford(i,1,j) + waforg(i,3) + tzord(i,j) 
        waford(i,2,j)= -e*waford(i,2,j) + waforg(i,3) + tzord(i,j) 
        waford(i,3,j)=waforg(i,1)+e*xaf(j) 
      END DO 
    END DO 
  END IF                                                  ! END WING

  IF (j4.NE.0) THEN                                     ! START FINS
    DO lq=1,ABS(nf) 
      DO j=2,1,-1                              
        e=.01*finorg(lq,j,4)   ! chord                              
        e2=finorg(lq,j,2)      ! y-coor                             
        DO k=1,nfinor 
          ee=finord(lq,j,k)*e           ! delta-y                    
          finord(lq,j,k)=e2+ee         ! outboard y                 
          finx2(lq,j,k)=e2-ee          ! inboard  y                 
          finx3(lq,j,k)=finorg(lq,j,1) + e*xfin(lq,k)   ! x-coor    
        END DO 
      END DO 
    END DO 
  END IF                                                  ! END FINS

  IF (j5.NE.0) THEN                                  ! START CANARDS
    DO nn=1,ncan 
      DO i=2,1,-1    ! why?                                         
        e=.01*canorg(nn,i,4) 
        e3=canorg(nn,i,3) 
        DO j=1,ABS(ncanor) 
          canord(nn,i,j)=e*canord(nn,i,j)+e3      ! upper z         
          canor1(nn,i,j)=-e*canor1(nn,i,j)+e3     ! lower z         
          canorx(nn,i,j)=canorg(nn,i,1)+e*xcan(nn,j)  ! x-coor       
        END DO
      END DO 
    END DO 
  END IF                                               ! END CANARDS

  WRITE (*,*) 'Configuration Data has been transformed'
  RETURN
END Subroutine TransformConfigurationData   ! -------------------------------

!+                                                                      
      SUBROUTINE ComputeConfigurationBounds()
!   --------------------------------------------------------------------
!     PURPOSE - Compute limits                                          
!                                                                       
!     NOTES-     ***********  NOT NEEDED ANYMORE ***********************
!                                                                       
      IMPLICIT NONE 
!***********************************************************************
!     L O C A L   V A R I A B L E S                                    *
!***********************************************************************
      INTEGER:: i,j,k 
      INTEGER:: lq 
      INTEGER:: nn 

!      REAL xaf(30),waforg(20,4),w(20,4)
!      REAL waford(20,3,30),tzord(20,30),ordmax(20,2) 
!      COMMON/WING/xaf,waforg,w,waford,tzord,ordmax 
!                                                                       
!      REAL finorg(6,2,4),xfin(6,10) 
!      REAL finord(6,2,10),finx2(6,2,10),finx3(6,2,10) 
!      REAL finmx1(6),finmx2(6),finth1(6),finth2(6) 
!      COMMON /FIN/ finorg,xfin,finord,                                  &
!     & finx2,finx3,finmx1,finmx2,finth1,finth2                          
!                                                                       
!      REAL canorg(2,2,4),xcan(2,10),canord(2,2,10) 
!      REAL canor1(2,2,10),canorx(2,2,10),canmax(2,2,2) 
!      COMMON /CAN/ canorg,xcan,canord,canor1,canorx,canmax 
!-----------------------------------------------------------------------
      WRITE(*,*) 'Begin subroutine ComputeConfigurationBounds' 
!                                                                       
                                                            ! START WING
      IF (j1 .NE. 0) THEN 
        DO I=1,NWAF 
          ordmax(i,1)=waford(i,1,1) 
          ordmax(i,2)=waford(i,2,1) 
          DO j=2,nwafor 
            ordmax(i,1)=MAX(ordmax(i,1),waford(i,1,j)) 
            ordmax(i,2)=MIN(ordmax(i,2),waford(i,2,j)) 
          END DO 
        END DO 
                                                              ! END WING
      END IF 
!                                                                       
                                                            ! START FINS
      IF (j4.NE.0) THEN 
         DO lq=1,ABS(nf) 
           finmx1(lq)=0. 
           finmx2(lq)=0. 
           DO k=1,ABS(nfinor) 
             finmx1(lq)=MAX(finmx1(lq),finord(lq,1,k)) 
             finmx2(lq)=MAX(finmx2(lq),finord(lq,2,k)) 
           END DO 
           finth1(lq)=2.0*(finmx1(lq)-finorg(lq,1,2)) 
           finth2(lq)=2.0*(finmx2(lq)-finorg(lq,2,2)) 
        END DO 
                                                              ! END FINS
      END IF 
!                                                                       
                                                         ! START CANARDS
      IF (j5.NE.0) THEN 
         DO nn=1,ncan 
           DO i=1,2 
             DO j=2,ABS(ncanor) 
               k=j-1 
               IF (canord(nn,i,k).GE.canord(nn,i,j)) GO TO 470 
             END DO 
  470        canmax(nn,i,1)=canord(nn,i,k) 
             DO j=2,ABS(ncanor) 
               k=j-1 
               IF (canor1(nn,i,k).LE.canor1(nn,i,j)) GO TO 480 
             END DO 
  480        canmax(nn,i,2)=canor1(nn,i,k) 
           END DO 
         END DO 
                                                           ! END CANARDS
      END IF 
!                                                                       
      WRITE (*,*) 'Configuration Bounds have been computed'
  RETURN
END Subroutine ComputeConfigurationBounds   ! -------------------------------

!+
      SUBROUTINE DumpConfigurationData(lun) 
!   --------------------------------------------------------------------
!     PURPOSE - Dump output to debug file                               
      IMPLICIT NONE

  INTEGER,INTENT(IN):: lun   ! logical unit number of the input file

  INTEGER:: i,j,k

!      REAL xaf(30),waforg(20,4),w(20,4) 
!      REAL waford(20,3,30),tzord(20,30),ordmax(20,2) 
!      COMMON /WING/ xaf,waforg,w,waford,tzord,ordmax 
!                                                                       
!      REAL xfus(30,4),zfus(30,4),fusard(30,4),fusrad(30,4),sfus(30,30,8) 
!      COMMON /FUSE/ xfus,zfus,fusard,fusrad,sfus 
!                                                                       
!      REAL podorg(9,3),xpod(9,30),podord(9,30) 
!      COMMON /POD/ podorg,xpod,podord 
!                                                                       
!      REAL finorg(6,2,4),xfin(6,10) 
!      REAL finord(6,2,10),finx2(6,2,10),finx3(6,2,10) 
!      REAL finmx1(6),finmx2(6),finth1(6),finth2(6) 
!      COMMON /FIN/ finorg,xfin,finord,                                  &
!     & finx2,finx3,finmx1,finmx2,finth1,finth2                          
!                                                                       
!      REAL canorg(2,2,4),xcan(2,10),canord(2,2,10) 
!      REAL canor1(2,2,10),canorx(2,2,10),canmax(2,2,2) 
!      COMMON /CAN/ canorg,xcan,canord,canor1,canorx,canmax 
!-----------------------------------------------------------------------
  WRITE(*,*) 'Begin subroutine DumpConfigurationData' 
  WRITE(lun, '(7I5)') j0,j1,j2,j3,j4,j5,j6 
  WRITE(lun,*) nwaf,nwafor 
  WRITE(lun,*) nfus 
  WRITE(lun,*) nradx 
  WRITE(lun,*) nforx 
  WRITE(lun,*) np,npodor 
  WRITE(lun,*) nf,nfinor 
  WRITE(lun,*) ncan,ncanor 

  IF (j1 .NE. 0) THEN                                         ! WING
    DO j=1,ABS(nwaf) 
      WRITE(lun, '(A,I2,A,I2)' ) 'Airfoil #', j, '  of ', ABS(nwaf) 
      WRITE(lun, '(8X,A,8X,A,5X,A,2X,A)' )                         &
                    'x','y','z-upper','z-lower'                        
      WRITE(lun, '(I3,5F9.4)' )                                     &
           (i, waford(j,3,i), waforg(j,2), waford(j,1,i),              &
               waford(j,2,i), tzord(j,i), i=1,ABS(nwafor) )            
    END DO 

    WRITE(lun,*) 'ORDMAX' 
    WRITE(lun, '(2F15.6)' ) (ordmax(i,1),ordmax(i,2), i=1,nwaf) 
  END IF 

  IF (j2 .NE. 0) THEN                                     ! FUSELAGE
    DO k=1,nfus 
      WRITE(lun, '(A,I2,A,I2)' ) 'Fuselage Segment #',k, '  of',nfus 
      IF (j2.EQ.1) THEN 
        WRITE(lun, '(13X,A,9X,A,9X,A)' ) 'x','y','z' 
        DO j=1,nforx(k) 
          WRITE(lun, '(2I4,3F10.5)' ) (j,i, xfus(j,k),              &
                        sfus(i,j,k+k-1), sfus(i,j,k+k), i=1,nradx(k))
        END DO 
      ELSE 
        WRITE(lun, '(9X,A,9X,A,9X,A,9X,A)' ) 'x','z','S','r' 
        WRITE(lun, '(I4,4F10.5)' ) (j, xfus(j,k), zfus(j,k),        &
                              fusard(j,k), fusrad(j,k), j=1,nforx(k)) 
      END IF 
    END DO 
  END IF 

  IF (j3.NE.0) THEN                                           ! PODS
    DO k=1,np 
      WRITE(lun, '(A,I2,A,I2)' ) ' Pod #',k, '   of ', np 
      WRITE(lun, '(9X,A,9X,A,9X,A,9X,A)' ) 'x','y','z','r' 
      WRITE(lun, '(I4,4F10.5)' ) (j, xpod(k,j)+podorg(k,1),        &
                podorg(k,2), podorg(k,3), podord(k,j), j=1,npodor)
    END DO 
  END IF 

  IF (j4.NE.0) THEN                                           ! FINS
    DO k=1,ABS(nf) 
      WRITE(lun, '(A,I2,A,I2)' ) 'Fin #', k, '   of ', ABS(nf) 
      WRITE(lun, '(15X,A,35X,A)' ) 'root','tip' 
      WRITE(lun, '(8X,A,6X,A,5X,A,10X,A,6X,A,5X,A)' )              &
       'x','y-inbd  y-outbd','z',                                   &
       'x','y-inbd  y-outbd','z'                                    
      WRITE(lun, '(I3,4F9.4,F11.4,3F9.4)' ) (j,                    &
        finx3(k,1,j), finx2(k,1,j), finord(k,1,j), finorg(k,1,3),   &
        finx3(k,2,j), finx2(k,2,j), finord(k,2,j), finorg(k,2,3),   &
        j=1,ABS(nfinor))                                            
    END DO 
  END IF 

  IF (j5.NE.0) THEN                                        ! CANARDS
    DO k=1,ABS(ncan) 
      WRITE(lun, '(A,I2,A,I2)' ) 'Canard #', k, '   of ', ABS(ncan) 
      WRITE(lun, '(15X,A,35X,A)' ) 'root','tip' 
      WRITE(lun, '(8X,A,8X,A,5X,A,7X,A,8X,A,5X,A)' )                &
        'x','y','z-upper  z-lower',                                  &
        'x','y','z-upper  z-lower'                                   
      WRITE(lun, '(I3,4F9.4,F11.4,3F9.4)' ) (j,                     &
        canorx(k,1,j), canorg(k,1,2), canord(k,1,j), canor1(k,1,j),  &
        canorx(k,2,j), canorg(k,2,2), canord(k,2,j), canor1(k,2,j),  &
        j=1,ABS(ncanor))                                             
    END DO 
  END IF 

  WRITE (*,*) 'Configuration Data has been dumped'
  RETURN
END Subroutine DumpConfigurationData   ! ------------------------------------

!+                                                                      
SUBROUTINE CountPanels(lun) 
!   --------------------------------------------------------------------
! PURPOSE - Count the number of panels (including images)           
      IMPLICIT NONE 

  INTEGER,INTENT(IN):: lun   ! logical unit number of the input file

  INTEGER,PARAMETER:: NSECP=6     ! MUST agree with PodsToSketch ???

  INTEGER:: j,k
  INTEGER:: nwing,nfuselage,npods,nfins,ncanards 
  INTEGER:: ntotal 
!                                                                       
!      REAL podorg(9,3),xpod(9,30),podord(9,30) 
!      COMMON /POD/ podorg,xpod,podord 
!                                                                       
!      REAL finorg(6,2,4),xfin(6,10) 
!      REAL finord(6,2,10),finx2(6,2,10),finx3(6,2,10) 
!      REAL finmx1(6),finmx2(6),finth1(6),finth2(6) 
!      COMMON /FIN/ finorg,xfin,finord,                                  &
!     & finx2,finx3,finmx1,finmx2,finth1,finth2                          
!                                                                       
!-----------------------------------------------------------------------
  WRITE(*,*) 'Begin subroutine CountPanels' 

  nwing= 2*(nwaf-1)*(ABS(nwafor)-1)                         ! WING  

  nfuselage=0
  IF (j2.NE.0) THEN                                       ! FUSELAGE
    DO k=1,nfus 
      nfuselage=nfuselage + (nradx(k)-1)*(nforx(k)-1) 
    END DO 
  END IF 

  npods=0
  IF (j3.NE.0) THEN                                           ! PODS
    DO j=1,np 
      npods=npods+1 
      IF (podorg(j,2) .NE. 0.0) npods=npods+1                ! image
    END DO 
    npods=npods*NSECP*(npodor-1) 
  END IF 

  nfins=0                                                 
  IF (j4.NE.0) THEN                                           ! FINS
    DO j=1,ABS(nf) 
      IF (finorg(j,1,2).NE.0.0 .OR. finorg(j,2,2).NE.0.0) THEN 
        nfins=nfins+2                         ! inboard and outboard
      ELSE 
        nfins=nfins+1                               ! centerline fin
      END IF 
    END DO 
    nfins=nfins*(ABS(nfinor)-1) 
  END IF 
                                                
  ncanards=2*ABS(ncan)*(ABS(ncanor)-1)                     ! CANARDS

  ntotal=nwing+nfuselage+npods+nfins+ncanards 
  WRITE(*,*) nwing,     ' wing panels' 
  WRITE(*,*) nfuselage, ' fuselage panels' 
  WRITE(*,*) npods,     ' pod panels' 
  WRITE(*,*) nfins,     ' fin panels' 
  WRITE(*,*) ncanards,  ' canard panels' 
  WRITE(*,*) ntotal,    ' total panels' 
  WRITE(*,*) 'Panels have been counted' 
  WRITE(lun,*) nwing,     ' wing panels' 
  WRITE(lun,*) nfuselage, ' fuselage panels' 
  WRITE(lun,*) npods,     ' pod panels' 
  WRITE(lun,*) nfins,     ' fin panels' 
  WRITE(lun,*) ncanards,  ' canard panels' 
  WRITE(lun,*) ntotal,    ' total panels' 
  WRITE(lun,*) 'Panels have been counted'
  RETURN
END Subroutine CountPanels   ! ----------------------------------------------

!+                                                                      
      SUBROUTINE TestForValidity()
!   --------------------------------------------------------------------
!     PURPOSE - Check the data for validity. All errors are fatal.      
      IMPLICIT NONE 

      INTEGER:: j,nsum
!-----------------------------------------------------------------------
  IF (j0.LT.0  .OR. j0.GT.1) THEN
    WRITE(*,*) 'Invalid j0'
    STOP
  END IF
  IF (j1.LT.-1 .OR. j1.GT.1) THEN
    WRITE(*,*) 'Invalid j1'
    STOP
  END IF
  IF (j2.LT.-1 .OR. j2.GT.1) THEN
    WRITE(*,*) 'Invalid j2'
    STOP
  END IF
  IF (j3.LT.0  .OR. j3.GT.1) THEN
    WRITE(*,*) 'Invalid j3'
    STOP
  END IF
  IF (j4.LT.0  .OR. j4.GT.1) THEN
    WRITE(*,*) 'Invalid j4'
    STOP
  END IF
  IF (j5.LT.0  .OR. j5.GT.1) THEN
    WRITE(*,*) 'Invalid j5'
    STOP
  END IF
  IF (j6.LT.-1 .OR. j6.GT.1) THEN
    WRITE(*,*) 'Invalid j6'
    STOP
  END IF
  IF (j1 .EQ. 0) THEN 
    IF (nwaf.NE.0) THEN
      WRITE(*,*) 'j1=0 and nwaf !=0 is inconsistent'
      STOP
    END IF

    IF (nwaf.LT.2 .OR. nwaf.GT.20) THEN
      WRITE(*,*) 'Invalid NWAF'
      STOP
    END IF

    IF (nwafor.LT.3 .OR. nwafor.GT.30) THEN
      WRITE(*,*) 'Invalid NWAFOR'
      STOP
    END IF
  END IF

  IF (j2 .EQ. 0) THEN 
    IF (nfus.NE.0) THEN
      WRITE(*,*) 'j2=0 and nfus !=0 is inconsistent'
      STOP
    END IF
  ELSE 
    IF (nfus.LT.1 .OR. nfus.GT.4) THEN
      WRITE(*,*) 'Invalid NFUS'
      STOP
    END IF
    nsum=1 
    DO j=1,nfus 
       IF (nradx(j).LT.3 .OR. nradx(j).GT.50) THEN
         WRITE(*,*) 'Invalid NRADX'
         STOP
       END IF
       nsum=nsum+nforx(j) 
    END DO 
    IF (nsum>101) THEN
      WRITE(*,*) 'Too many fuselage points.'
      STOP
    END IF
    IF(j2.EQ.1 .AND. j6.EQ.1) THEN
      WRITE(*,*) 'j2 and j6 cannot both be 1'
      STOP
    END IF
  END IF 

  IF (j3 .EQ. 0) THEN 
    IF (np.NE.0) THEN
      WRITE(*,*) 'j3=0 and np !=0 is inconsistent'
      STOP
    END IF
    IF (npodor.NE.0) THEN
      WRITE(*,*) 'j3=0 and npodor !=0 is inconsistent'
      STOP
    END IF
  ELSE 
    IF (np.LT.1 .OR. np.GT.9) THEN
      WRITE(*,*) 'Invalid NP'
      STOP
    END IF
    IF (npodor.LT.4 .OR. npodor.GT.30) THEN
      WRITE(*,*) 'Invalid NPODOR'
      STOP
    END IF
  END IF 

  IF (j4 .EQ. 0) THEN
    IF (nf.NE.0) THEN
      WRITE(*,*) 'j4=0 and nf !=0 is inconsistent'
      STOP
    END IF
    IF (nfinor.NE.0) THEN
      WRITE(*,*) 'j4=0 and nfinor !=0 is inconsistent'
      STOP
    END IF
  ELSE 
    IF (ABS(nf).LT.1 .OR. ABS(nf).GT.6) THEN
      WRITE(*,*) 'Invalid NF'
      STOP
    END IF
    IF (nfinor.LT.3 .OR. nfinor.GT.10) THEN
      WRITE(*,*) 'Invalid NFINOR'
      STOP
    END IF
  END IF 

  IF (j5 .EQ. 0) THEN 
    IF (ncan.NE.0) THEN
      WRITE(*,*) 'j5=0 and ncan !=0 is inconsistent'
      STOP
    END IF
    IF (ncanor.NE.0) THEN
      WRITE(*,*) 'j5=0 and ncanor !=0 is inconsistent'
      STOP
    END IF
  ELSE 
    IF (ABS(ncan).LT.1 .OR. ABS(ncan).GT.2) THEN
      WRITE(*,*) 'Invalid NCAN'
      STOP
    END IF
    IF (ABS(ncanor).LT.3.OR.ABS(ncanor).GT.10) THEN
      WRITE(*,*) 'Invalid NCANOR'
      STOP
    END IF
  END IF 

  WRITE(*,*) 'Valid data.'
  RETURN
END Subroutine TestForValidity   ! ------------------------------------------

!+
SUBROUTINE WingToWgs(lun) 
!   --------------------------------------------------------------------
! PURPOSE - Write the wing data in WGS format                       
!                                                                       
!     NOTES- Each wing surface has (nwaf-1)*(nwafor-1) panels...........
!        (There are 2 surfaces: upper,lower)                            
!        So a total of 2*(nwaf-1)*(nwafor-1) panels                     
!                                                                       
IMPLICIT NONE 
  INTEGER,INTENT(IN):: lun 

  INTEGER:: i,j
!                                                                       
!      REAL xaf(30),waforg(20,4),w(20,4) 
!      REAL waford(20,3,30),tzord(20,30),ordmax(20,2) 
!      COMMON/WING/xaf,waforg,w,waford,tzord,ordmax 
!-----------------------------------------------------------------------

      IF (j1 .EQ. 0) RETURN 
!..... Upper surface ... (nwaf-1)*(nwafor-1) panels.....................
      WRITE(lun,*) '''WING UPPER''' 
      WRITE(lun,*) 1, nwaf, nwafor, '  0   0 0 0   0 0 0   1 1 1   1' 
      DO j=1,nwaf 
        DO i=1,nwafor 
          WRITE(lun,*) waford(j,3,i), waforg(j,2), waford(j,1,i) 
        END DO 
      END DO 
!..... Lower surface ...................................................
      WRITE(lun,*) '''WING LOWER''' 
      WRITE(lun,*) 2, nwaf, nwafor, '  0   0 0 0   0 0 0   1 1 1   1' 
      DO j=1,nwaf 
        DO i=1,nwafor 
          WRITE(lun,*) waford(j,3,i), waforg(j,2), waford(j,2,i) 
        END DO 
      END DO 
      WRITE(*,*) 'Wing data sent to Wgs.'
  RETURN
END Subroutine WingToWgs   ! ------------------------------------------------

!+                                                                      
SUBROUTINE FuselageToWgs(lun) 
!   --------------------------------------------------------------------
! PURPOSE - Write the fuselage data in WGS format                   
      IMPLICIT NONE 

      INTEGER,INTENT(IN):: lun

      INTEGER:: i,j,k
      REAL:: theta,c,s,r 

!      REAL xfus(30,4),zfus(30,4),fusard(30,4),fusrad(30,4),sfus(30,30,8) 
!      COMMON/FUSE/xfus,zfus,fusard,fusrad,sfus 
!-----------------------------------------------------------------------
      IF (j2 .EQ. 0) RETURN 
      DO k=1,nfus 
        WRITE(lun,*) '''FUSELAGE SEGMENT #', k, '''' 
        WRITE(lun,*)                                                    &
     &          20+k, nradx(k), nforx(k), ' 0   0 0 0   0 0 0   1 1 1 0'
        IF (j2 .EQ. 1) THEN 
          DO j=1,nradx(k) 
            DO i=1,nforx(k) 
              WRITE(lun,*) xfus(i,k), sfus(j,i,k+k-1), sfus(j,i,k+k) 
            END DO 
          END DO 
        ELSE 
          DO j=1,nradx(k) 
            theta=REAL(j-1)*PI/REAL(nradx(k)-1) 
            s=SIN(theta) 
            c=COS(theta) 
            DO i=1,nforx(k) 
              r=SQRT(fusard(i,k)/PI) 
              WRITE(lun,*) xfus(i,k), r*s, zfus(i,k)+r*c 
            END DO 
          END DO 
        END IF 
      END DO 
      WRITE(*,*) 'Fuselage data sent to Wgs.'
  RETURN
END Subroutine FuselageToWgs   ! --------------------------------------------

!+                                                                      
SUBROUTINE PodsToWgs(lun) 
!   --------------------------------------------------------------------
! PURPOSE - Write the pod data in WGS format                        
      IMPLICIT NONE 
      INTEGER,INTENT(IN):: lun 

!      INTEGER NSECP
      INTEGER,PARAMETER:: NSECP=6
!      REAL PI 
!      PARAMETER (PI=3.14159265) 
  REAL,PARAMETER:: DELTHETA=PI/NSECP

      INTEGER:: i,j,k,n
      REAL:: theta 
      REAL:: c,s 

!-----------------------------------------------------------------------
      IF (j3 .EQ. 0) RETURN 
      DO k=1,np 
        IF (podorg(k,2) .EQ. 0) THEN 
          n=1+NSECP 
        ELSE 
          n=1+NSECP+NSECP 
        END IF 
!                                                                       
        WRITE(lun,*) '''POD #', k, '''' 
        WRITE(lun,*) 30+k, n, npodor, ' 0   0 0 0   0 0 0   1 1 1 0' 
        DO j=1,n 
          theta=REAL(j)*DELTHETA 
          s=SIN(theta) 
          c=COS(theta) 
!                                                                       
          DO i=1,npodor 
            WRITE(lun,*) podorg(k,1) + xpod(k,i),                       &
     &                   podorg(k,2) + podord(k,i)*s,                   &
     &                   podorg(k,3) + podord(k,i)*c                    
          END DO 
        END DO 
      END DO 
      WRITE(*,*) 'Pod data sent to Wgs.'
  RETURN
END Subroutine PodsToWgs   ! ------------------------------------------------

!+                                                                      
SUBROUTINE FinsToWgs(lun) 
!   --------------------------------------------------------------------
! PURPOSE - Write the fin data in WGS format                        
      IMPLICIT NONE 

      INTEGER,INTENT(IN):: lun

      INTEGER:: i,k

!      REAL finorg(6,2,4),xfin(6,10) 
!      REAL finord(6,2,10),finx2(6,2,10),finx3(6,2,10) 
!      REAL finmx1(6),finmx2(6),finth1(6),finth2(6) 
!      COMMON /FIN/ finorg,xfin,finord,                                  &
!     & finx2,finx3,finmx1,finmx2,finth1,finth2                          
!-----------------------------------------------------------------------
      IF (j4 .EQ. 0) RETURN 
      DO k=1,ABS(nf) 
        WRITE(lun,*) '''FIN #', k, '   OUTBOARD''' 
        WRITE(lun,*) 39+k+k, 2, nfinor, ' 0   0 0 0   0 0 0   1 1 1 0' 
        DO i=1,nfinor 
          WRITE(lun,*) finx3(k,1,i), finord(k,1,i), finorg(k,1,3) 
        END DO 
        DO i=1,nfinor 
          WRITE(lun,*) finx3(k,2,i), finord(k,2,i), finorg(k,2,3) 
        END DO 

        IF (finorg(k,1,2).NE.0.0 .OR. finorg(k,2,2).NE.0.0) THEN 
          WRITE(lun,*) '''FIN #', k, '   INBOARD''' 
          WRITE(lun,*) 40+k+k, 2, nfinor, ' 0   0 0 0   0 0 0   1 1 1 0' 
          DO i=1,nfinor 
            WRITE(lun,*) finx3(k,1,i), finx2(k,1,i), finorg(k,1,3) 
          END DO 
          DO i=1,nfinor 
            WRITE(lun,*) finx3(k,2,i), finx2(k,2,i), finorg(k,2,3) 
          END DO 
        END IF 
      END DO 
      WRITE(*,*) 'Fin data sent to Wgs.'
  RETURN
END Subroutine FinsToWgs   ! ------------------------------------------------

!+                                                                      
SUBROUTINE CanardsToWgs(lun) 
!   --------------------------------------------------------------------
! PURPOSE - Write the canard data in WGS format                     
      IMPLICIT NONE 
      INTEGER,INTENT(IN):: lun 

      INTEGER:: i,k

!      REAL canorg(2,2,4),xcan(2,10),canord(2,2,10) 
!      REAL canor1(2,2,10),canorx(2,2,10),canmax(2,2,2) 
!      COMMON /CAN/ canorg,xcan,canord,canor1,canorx,canmax 
!-----------------------------------------------------------------------
      IF (j5 .EQ. 0) RETURN 
      DO k=1,ncan 
        WRITE(lun,*) '''CANARD #', k, '  UPPER''' 
        WRITE(lun,*) 59+k+k, 2, ncanor, ' 0   0 0 0   0 0 0   1 1 1 0' 
        DO i=1,ncanor 
          WRITE(lun,*) canorx(k,1,i), canorg(k,1,2), canord(k,1,i) ! root 
        END DO

        DO i=1,ncanor 
          WRITE(lun,*) canorx(k,2,i), canorg(k,2,2), canord(k,2,i) ! tip 
        END DO 

        WRITE(lun,*) '''CANARD #', k, '  LOWER''' 
        WRITE(lun,*) 60+k+k, 2, ncanor, ' 0   0 0 0   0 0 0   1 1 1 0' 
        DO i=1,ncanor 
          WRITE(lun,*) canorx(k,1,i), canorg(k,1,2), canor1(k,1,i) !root 
        END DO 
        DO i=1,ncanor
          WRITE(lun,*) canorx(k,2,i), canorg(k,2,2), canor1(k,2,i) ! tip 
        END DO 
      END DO 
                                                                        
      WRITE(*,*) 'Canard data sent to Wgs.'
  RETURN
END Subroutine CanardsToWgs   ! ------------------------------------------------

END Module WaveDragToWgsProcedures   ! =========================================

!+
PROGRAM WaveDragToWgs 
! ------------------------------------------------------------------------------
USE WaveDragToWgsProcedures
IMPLICIT NONE 
                     
  INTEGER,PARAMETER:: IN=1,WGS=2,DBG=3,TMP=4    ! units                                        

  CHARACTER(LEN=*),PARAMETER:: FAREWELL=                                &
     &   'File wd.wgs has been added to your directory.'
                                                                       
  CHARACTER(LEN=80):: abc
  INTEGER:: i 
!-------------------------------------------------------------------------------
  CALL Welcome()    
  OPEN(UNIT=TMP,FILE='wd2wgs.tmp',STATUS='REPLACE',ACTION='READWRITE')
  CALL CopyData(IN,TMP)
  CLOSE(UNIT=IN)
  REWIND(UNIT=TMP)

  WRITE(*,*) 'Computing configuration: ' 
  READ(TMP,'(A)') abc       ! first card
  WRITE(*,*) abc 
  WRITE(DBG,*) abc 

  READ(TMP,'(24I3)') j0,j1,j2,j3,j4,j5,j6,nwaf,nwafor,nfus,          &
     & (nradx(i),nforx(i),i=1,4),np,npodor,nf,nfinor,ncan,ncanor        

  CALL TestForValidity()
  CALL ReadConfigurationData(TMP)   ! input the configuration data   
  CALL TransformConfigurationData()
  CALL ComputeConfigurationBounds()
  CALL DumpConfigurationData(DBG) 
  CALL CountPanels(DBG) 

  CALL WingToWgs(WGS)
  CALL FuselageToWgs(WGS) 
  CALL PodsToWgs(WGS) 
  CALL FinsToWgs(WGS) 
  CALL CanardsToWgs(WGS) 

  WRITE(*,*) FAREWELL
  CLOSE(UNIT=WGS) 
  CLOSE(UNIT=DBG) 
  CLOSE(UNIT=TMP) 
  WRITE(*,*) 'wd2wgs has terminated successfully.'
  STOP

CONTAINS
!+
SUBROUTINE Welcome()
! ---------------------------------------------------------------------------
! PURPOSE - Greet user, get file name and open all files

  CHARACTER(LEN=*),PARAMETER:: GREETING = &
    'wd2wgs - convert wd configuration to WGS'
  CHARACTER(LEN=*),PARAMETER:: AUTHOR = &  
    ' Ralph L. Carmichael, Public Domain Aeronautical Software'
  CHARACTER(LEN=*),PARAMETER:: MODIFIER=' '   ! your name here

  INTEGER:: errCode
  CHARACTER(LEN=132):: fileName
  INTRINSIC:: GET_COMMAND_ARGUMENT, LEN_TRIM
!----------------------------------------------------------------------------
  WRITE(*,*) GREETING
  WRITE(*,*) AUTHOR 
  IF (MODIFIER .NE. ' ') WRITE(*,*) 'Modified by '//MODIFIER 
  WRITE(*,*) 'Version '//VERSION

  errCode=99
  CALL GET_COMMAND_ARGUMENT(1,fileName)
  IF (LEN_TRIM(fileName) > 0) THEN
    OPEN(UNIT=IN,FILE=fileName,STATUS='OLD',IOSTAT=errCode,ACTION='READ') 
  END IF

  IF (errCode /= 0) THEN
    DO
      WRITE(*,*) 'Enter the name of the WaveDrag input file: ' 
      READ  (*, '(A)' ) fileName 
      IF (fileName .EQ. ' ') STOP 
      OPEN(UNIT=IN,FILE=fileName,STATUS='OLD',IOSTAT=errCode,ACTION='READ')
      IF (errCode == 0) EXIT
      OPEN(UNIT=IN,FILE=Trim(fileName)//'.inp', &
        STATUS='OLD',IOSTAT=errCode,ACTION='READ')
      IF (errCode == 0) EXIT
      WRITE (*,*) ' Unable to open this file. Try again' 
    END DO
  END IF
  INQUIRE(UNIT=IN, NAME=fileName)
                                                                      
  OPEN(UNIT=WGS, FILE='wd.wgs', STATUS='REPLACE', ACTION='WRITE') 
  WRITE(WGS,*) 'Created by wd2wgs from: '// Trim(fileName)
                                                                        
  OPEN(UNIT=DBG,FILE='wd2wgs.dbg',STATUS='REPLACE',ACTION='WRITE') 
  WRITE(DBG,*) 'Created by wd2wgs from: '//Trim(fileName)

  RETURN
END Subroutine Welcome   ! --------------------------------------------------


END Program WaveDragToWGS   ! ===============================================
