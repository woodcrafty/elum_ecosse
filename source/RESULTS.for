C------------------------------------------------------------
C Rothamsted Nitrogen Turnover Model
C Results Subroutines
C     N.J.Bradbury
C
C Version - 93-7
C     J.U.Smith
C
C Link in Submodels NMMAIN,NMSET,NMCROP,NMORG,NMNIT,NMWAT,NMAPP,NMRES
C
C---------------------------------------------------------------------
C

      SUBROUTINE SETRES(BALANCE,CO2,SX,SXORGN)
C
C Subroutine to set results
C
C Common Blocks
C
      INCLUDE 'G.FOR'
      INCLUDE 'TOTVAR2.FOR'
C
C Dimension Statement
C
      DIMENSION BALANCE(20)
	REAL CO2(2)
C
C Initialize BALANCE array
C
      DO 50 I=1,14
       BALANCE(I)=0
 50   CONTINUE
C
C Set BALANCE(1) = C in HUM Pool
C Set BALANCE(2) = C in BIO Pool
C Set BALANCE(3) = C in RO Pool
C
      BALANCE(1)=HCARB0(1)+HCARB0(2)
      BALANCE(2)=BCARB0(1)+BCARB0(2)
      BALANCE(3)=DPMCARB0(1)+DPMCARB0(2)+
     &           RPMCARB0(1)+RPMCARB0(2)
c balance(11), balance(12) are initial dpm & rpm amount in pool
	BALANCE(11)=DPMCARB0(1)+DPMCARB0(2)
	BALANCE(12)=RPMCARB0(1)+RPMCARB0(2)
c      BALANCE(3)=RCARB0(1)+RCARB0(2)
C
C Initialize variables
C
      CALL SETVAR1
      SX=SXORGN
      TRNIN=0
	CO2(1)=0
	CO2(2)=0
C
C Leave SETRES
C
      RETURN
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE SETVAR1
C
C Subroutine to initalize variables
C
C Common Blocks
C
      INCLUDE 'TOTVAR2.FOR'
C
C Initialize Variables
C
       TCINP=0
       TRNIN=0
       TRNIN15=0
       IR_TS=0
       TC_TS=0
       ITROOT=0
C
C Leave SETVAR1
C
       RETURN
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE STARTRES(IRYEAR,IK,RSTART,DPMNIT0,RPMNIT0,
     &                    BSTART,BNIT0,HNIT0,
     &                    RSTART15,DPMNLAB0,RPMNLAB0,
     &                    BSTART15,BNLAB0,HNLAB0,
     &                    TOTAST,AMMN,TOTAST15,AMMN15,
     &                    TOTNST,SOILN,TOTNST15,SOIL15)
     &                    
C
C Subroutine to save results at start of timestep
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      INTEGER MAXLAYER
	PARAMETER (MAXLAYER=60)
	INTEGER IL
C
C Variables passed to/from calling subroutine
C
      INTEGER IRYEAR,IK
	REAL AMMN(MAXLAYER),AMMN15(MAXLAYER)
	REAL SOILN(MAXLAYER),SOIL15(MAXLAYER)
	REAL TOTAST,TOTNST,TOTAST15,TOTNST15
	REAL RSTART,BSTART,DPMNIT0(MAXLAYER),RPMNIT0(MAXLAYER)
	REAL BNIT0(MAXLAYER),HNIT0(MAXLAYER)
	REAL RSTART15,BSTART15,DPMNLAB0(MAXLAYER),RPMNLAB0(MAXLAYER)
	REAL BNLAB0(MAXLAYER),HNLAB0(MAXLAYER)
C
C Set last week's end results to this week's start result
C
	RSTART=DPMNIT0(1)+RPMNIT0(1)+DPMNIT0(2)+RPMNIT0(2)
	BSTART=BNIT0(1)+BNIT0(2)+HNIT0(1)+HNIT0(2)
      RSTART15=DPMNLAB0(1)+RPMNLAB0(1)+DPMNLAB0(2)+RPMNLAB0(2)
      BSTART15=BNLAB0(1)+BNLAB0(2)+HNLAB0(1)+HNLAB0(2)
	TOTAST=0
	TOTAST15=0
	TOTNST=0
	TOTNST15=0
      DO 20 IL=1,3
	  TOTAST=TOTAST+AMMN(IL)
        TOTAST15=TOTAST15+AMMN15(IL)
20    CONTINUE
      DO 21 IL=1,12
        TOTNST=TOTNST+SOILN(IL)
        TOTNST15=TOTNST15+SOIL15(IL)
21    CONTINUE
      END
C
C------------------------------------------------------------------
C
      SUBROUTINE TS_RES1(IYEAR,IRYEAR,NXYEARS,SECONDS, 
     &                   IK,N_TS,SUM_TS,FIXEND,NSOW,ISOWN,MEND,
     &                   SX,SXORGN15,SORGN,SORGN15,RNIN15,
     &                   CACT,CACT15,CACTOT,CATOT15,
     &                   CLOSSX,CLOSSX15,VOLAT,VOLAT15,
     &                   ATM,ATM15,TDNIT,T15,
     &                   THISFERT,THISFERT15,WLEACH,WLEACH15,
     &                   DN15,                   
     &                   FYMFERT,FYMFERT15,
     &                   IANTHES,CONC,CONC15,SLEACH,
     &                   SOILN,AMMN,TOTAST,TOTNST,TOTAST15,TOTNST15,
     &                   DPMCARB0,RPMCARB0,BCARB0,HCARB0,
     &                   DPMNIT0,RPMNIT0,HNIT0,BNIT0,DENIT,
     &                   BSTART,RSTART,BSTART15,RSTART15,RNIN,
     &                   ICROP,SEEDIN,GNN2O,GNNO,GPNN2O,G15NN2O,G15NNO,
     &                   G15PNN2O,GDN2,GDN2O,G15DN2,G15DN2O)
C
C Subroutine to record weeks results
C
C
C Parameters
C
      PARAMETER (MAXCROP=36,MAXSOIL=50,MAXLAYER=60)
C
C Define variables
C
      INTEGER MAXLAYER1
	DATA MAXLAYER1 /60/
	INTEGER MAXDEPTH
	DATA MAXDEPTH /300/ 
      REAL FTEMP,FTEMP15,SLEACH(MAXLAYER)
      INTEGER SUM_TS,FIXEND,ISOWN,ICROP
      REAL SECONDS
	REAL SOILN(MAXLAYER),AMMN(MAXLAYER)
	REAL HCARB0(MAXLAYER),DPMCARB0(MAXLAYER),RPMCARB0(MAXLAYER)
      REAL DPMNIT0(MAXLAYER),RPMNIT0(MAXLAYER),HNIT0(MAXLAYER)
	REAL BCARB0(MAXLAYER),BNIT0(MAXLAYER),DENIT
	REAL TOTAST,TOTNST,TOTAST15,TOTNST15
	REAL RNIN,SEEDIN
	REAL GNN2O,GNNO,GPNN2O,G15NN2O,G15NNO,G15PNN2O					
	REAL GDN2,GDN2O,G15DN2,G15DN2O									
C
C Common block
C
      INCLUDE 'CROPPAR.FOR'
      INCLUDE 'A.FOR'
      INCLUDE 'FIXER.FOR'
      INCLUDE 'LASTHARV.FOR'
      INCLUDE 'LEVEL.FOR'
      INCLUDE 'N15.FOR'
      INCLUDE 'POOL.FOR'
      INCLUDE 'RES.FOR'
      INCLUDE 'TOTVAR.FOR'
      INCLUDE 'TOTVAR15.FOR'
      INCLUDE 'TOTVAR2.FOR'
	SAVE
C
C Initialize variables at sowing
C
      IF(IK.LT.NSOW)THEN
        CACL=CACTOT
      ELSE
        CACL=0
      ENDIF
C
      START=RSTART+BSTART+TOTNST+TOTAST
C
C Sum inputs during week
C
      WATM=ATM+FIXN
      SUBTOTIN=THISFERT+WATM+RNIN+FYMFERT
C
C Sum N at start
C
      TOTIN=START+SUBTOTIN
C
C Calculate total nitrate and ammonia at end of week
C
      TOTNIT=0
      TOTAMM=0
C
C Write out balance sheet 0-150cm
C Nitrate...
C
c      DO 400 IL=1,(150*MAXLAYER1/MAXDEPTH)
      DO 400 IL=1,(300*MAXLAYER1/MAXDEPTH)
        TOTNIT=TOTNIT+SOILN(IL)
400   CONTINUE
C
C Ammonium...
C
c      DO 500 IL=1,(150*MAXLAYER1/MAXDEPTH)
      DO 500 IL=1,(300*MAXLAYER1/MAXDEPTH)
        TOTAMM=TOTAMM+AMMN(IL)
500   CONTINUE
C
C Plant debris and biomas + humus
C
      ROOTN=0
	ROOTDPM=0
	ROOTRPM=0
	BIOHUM=0
      ROOTN15=0
c      DO 550 IL=1,(150*MAXLAYER1/MAXDEPTH)
      DO 550 IL=1,(300*MAXLAYER1/MAXDEPTH)
        ROOTN=ROOTN+DPMNIT0(IL)+RPMNIT0(IL)
	  ROOTDPM=ROOTDPM+DPMNIT0(IL)
	  ROOTRPM=ROOTRPM+RPMNIT0(IL)
        BIOHUM=BIOHUM+BNIT0(IL)+HNIT0(IL)
550   CONTINUE
      ENDRES=ROOTN+BIOHUM+TOTNIT+TOTAMM
C
C Sum outputs during week
C
      SUBTOTOUT=DENIT+CACT+WLEACH+VOLAT-CLOSSX
C
C Sum N at end of week
C
      TOTOUT=ENDRES+SUBTOTOUT
C
C Write out Weekly N Balance sheet
C
      BALIN=BALIN+THISFERT+FYMFERT+RNIN+WATM+SEEDIN
	BALOUT=BIOHUM+ROOTDPM+ROOTRPM+TOTNIT+TOTAMM
	BALOUT=BALOUT+DENIT+VOLAT-CLOSSX+CACT+WLEACH
	DIFF=BALIN-BALOUT
      WRITE(58,5801) N_TS,IK,IYEAR, BIOHUM,ROOTN,
     1    ROOTDPM,ROOTRPM,TOTNIT,TOTAMM,THISFERT+FYMFERT,RNIN,
     2    WATM,SEEDIN,TDNIT,DENIT,GNN2O,GNNO,GPNN2O,G15NN2O,G15NNO,
     3    G15PNN2O,GDN2,GDN2O,G15DN2,G15DN2O,VOLAT-CLOSSX,CACT,WLEACH,
     4    CLOSSX,BALIN,BALOUT,DIFF
	BALIN=BIOHUM+ROOTDPM+ROOTRPM+TOTNIT+TOTAMM

 5801 FORMAT (2X,3I5,1X,F7.1,5F8.2,16F8.2,4F8.2,3F12.1)
C
      CALL GETCAC(CACTOT,CACT,CATOT15,CACT15)
      RETURN
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE CBALANCE(IYEAR,CO2,CH4,CH4TOAIR,IK,TFYMC,N_TS,TCINP,
     &                   HCARB0,BCARB0,DPMCARB0,RPMCARB0,IOM,PIANN,TEMP)
C
C Subroutine to write out the carbon balance
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER
	PARAMETER (MAXLAYER=60)
      INTEGER IL
	REAL TOTDPM,TOTRPM,TOTBIO,TOTHUM,TOTCO2,TOTCH4,TOTIOM
	REAL ANNCO2, PIANN
C
C Variables passed to/from calling subroutine
C
	REAL CH4(MAXLAYER)			! IN:CH4 emitted (kgC/ha/layer)
	REAL CH4TOAIR				! IN:CH4 emitted to atmosphere (kgC/ha/layer)
      REAL BCARB0(MAXLAYER),HCARB0(MAXLAYER),IOM(MAXLAYER)
	REAL DPMCARB0(MAXLAYER),RPMCARB0(MAXLAYER),CO2(MAXLAYER)
	INTEGER IYEAR,IK,N_TS
	REAL TCINP,TFYMC,TEMP
C
C Work out the totals
C
	TOTDPM=0
	TOTRPM=0
	TOTBIO=0
	TOTHUM=0
	TOTCO2=0
	TOTCH4=0
      TOTIOM=0
	DO 100 IL=1,MAXLAYER
	  TOTDPM=TOTDPM+DPMCARB0(IL)
	  TOTRPM=TOTRPM+RPMCARB0(IL)
	  TOTBIO=TOTBIO+BCARB0(IL)
	  TOTHUM=TOTHUM+HCARB0(IL)
	  TOTCO2=TOTCO2+CO2(IL)
	  TOTCH4=TOTCH4+CH4(IL)
        TOTIOM=TOTIOM+IOM(IL)
100   CONTINUE
      IF(IK.EQ.1)ANNCO2=0
	ANNCO2=ANNCO2+TOTCO2
C
C Write out the C results
C
      WRITE(57,5701)N_TS,IK,IYEAR,
     &             (DPMCARB0(IL),IL=1,10),
     &             (RPMCARB0(IL),IL=1,10),
     &             (BCARB0(IL),IL=1,10),
     &             (HCARB0(IL),IL=1,10),
     &             (CO2(IL),IL=1,10),
     &             (CH4(IL),IL=1,10),
     &              TOTDPM,TOTRPM,TOTBIO,TOTHUM,TOTCO2,TOTCH4,CH4TOAIR,
     &              TOTIOM,TCINP+TFYMC,ANNCO2,
     &              TOTDPM+TOTRPM+TOTBIO+TOTHUM+TOTIOM,PIANN,TEMP
 5701 FORMAT(1X,3I5,74(F12.2,1X))
	TFYMC = 0.0
      RETURN
      END
