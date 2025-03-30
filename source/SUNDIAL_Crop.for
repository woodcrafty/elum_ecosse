C------------------------------------------------------------
C
C Rothamsted Carbon and Nitrogen Turnover Model
C Soil Crop Routines
C
C Adapted from SUNDIAL (MAGEC)
C by Jo Smith & Kevin Coleman
C 02/03/05
C
C Modified for spatial applications
C by Jo Smith 
C Started 01/08/06
C
C-------------------------------------------------------------
C
C EXTERNAL SUBROUTINES
C 1. INIT_SUNDIAL_CROP
C 2. RUN1_SUNDIAL_CROP
C 3. RUN2_SUNDIAL_CROP
C
C-------------------------------------------------------------
C
      SUBROUTINE INIT_SUNDIAL_CROP(IYEAR, LCROP, ICROP, 
     &                             INSTRAW, IEND, ISTHARV,   
     &                             ISOWN, MEND, NSOW, JSTOP, L_TSS,
     &                             YLD, PREYLD, EXYLD,  
     &                             SORGC,SORGN, SXORGC, SXORGN,  
     &                             HZ1, TORGC, TORGN,TXORGC, TXORGN, 
     &                             OCROPN, RORGN,DDAYS,
     &                             ORGC, ORGN, XORGC, XORGN,
     &                             NXYEARS,NDATE,FIXEND,LHARV,IHARV,
     &                             NORGM,IORGM,NFERT,IFERT,IANTHES,
     &							 SECONDS,
C introduced for SWAT water routines
     &				             PLANTBIOMASS,PLANTHEIGHT,LAI,C1,T1,
     &							 canMax,CONVER_F,RRG,IROCKS)
C
C Subroutine to initialise SUNDIAL crop characteristics
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      INTEGER MAXLAYER
	PARAMETER (MAXLAYER=60)
      INTEGER MAXCROP,MAXORGM,MAXFERT
	PARAMETER (MAXCROP=36,MAXORGM=52,MAXFERT=5)
	INTEGER MCROP,IROCKS(0:MAXCROP),ISTART,ISAVE
	REAL ROOTS,CREQN,TOTC,CINP,ANINP,TCNEED,BREQN
      REAL CAO(5,0:MAXCROP),UR(3,0:MAXCROP),UT(3,0:MAXCROP) 
	REAL INC(2,0:MAXCROP),CRATE(0:MAXCROP),ANRATE(0:MAXCROP)
      REAL CFACT(0:MAXCROP),SEED(0:MAXCROP),C1(0:MAXCROP),T1(0:MAXCROP)
      REAL RRG(0:MAXCROP),RRX(0:MAXCROP)
      REAL CSC(4,0:MAXCROP),STRAW(4,0:MAXCROP)
C
C Variables passed to/from calling subroutine
C
      INTEGER IYEAR, LCROP, ICROP,INSTRAW,IEND,N_REMAIN
      INTEGER ISTHARV,N_STEPS,ISOWN,MEND,NSOW,JSTOP
	INTEGER L_TSS(0:MAXCROP)
	INTEGER IANTHES,NXYEARS,NDATE,FIXEND,LHARV,IHARV
	INTEGER NORGM,IORGM(MAXORGM),NFERT,IFERT(MAXFERT)
	REAL SECONDS
      REAL CONVER_F
	REAL YLD,PREYLD,EXYLD,SORGC,SORGN
      REAL SXORGC,SXORGN,HZ1(MAXLAYER),TORGC,TORGN,TXORGC
      REAL TXORGN,OCROPN,RORGN
	REAL DDAYS
	REAL ORGC, ORGN, XORGC, XORGN
C
C Variables required for SWAT water routines
C
	REAL PLANTBIOMASS		! Plant biomass at harvest [ kg ha-1]
	REAL PLANTHEIGHT		! Plant height at harvest [cm]
	REAL LAI				! Plant LAI at harvest
	REAL canMAX				! maximum amount of water that can be trapped in the canopy when the canopy is fully developed (mm H2O)
C
C Set crop parameters from file
C
      CALL PAR_CROP(L_TSS, IROCKS, CAO, UR, UT, INC, CRATE, 
     &                    ANRATE, CFACT, SEED, C1, T1, RRG, RRX,  
     &                    CSC, STRAW)
C
C Initialise variables
C
      CALL SETVARS_CROP()
C
C ...specifiy timing of management operations relative to previous harvest (or anthesis in previous year)
C
      CALL SETTIME_SUNDIAL_CROP(IYEAR,NXYEARS,NDATE,MEND,FIXEND,
     &                  ISTHARV,LHARV,LCROP,IANTHES,
     &                  SECONDS,CONVER_F,N_STEPS,N_REMAIN,L_TSS,
     &                  NSOW,ISOWN,IHARV,
     &                  NORGM,IORGM,NFERT,IFERT)
C
C Call subroutine MOREFIX to set parameters dependant on inputs
C from sites and to call SETC to set C inputs from previous crop
C
      CALL MOREFIX(IYEAR,MCROP,LCROP,ICROP,INSTRAW,
     &             YLD,PREYLD,EXYLD,SORGC,SORGN, 
     &             SXORGC,SXORGN,HZ1,TORGC,TORGN,TXORGC, 
     &             TXORGN,CAO,CSC,STRAW,
     &             ORGC,ORGN,XORGC,XORGN,
C introduced for SWAT water routines
     &				             PLANTBIOMASS,PLANTHEIGHT,LAI,C1,T1,
     &							 canMax,CONVER_F)
C
C In previous year to first growing season
C
      IF(IYEAR.EQ.0)THEN
C
C Set the crop N requirement
C
        CALL SETCROPN(IYEAR,INSTRAW,MCROP,LCROP,ICROP,
     &                    SXORGN,SORGN,YLD,PREYLD,CREQN, 
     &                    OCROPN,ROOTS,EXYLD,TOTC,
     &                    CFACT,STRAW,INC,SEED,UR,UT,XORGN,ORGN)
        XORGC=TXORGC-SXORGC
C
C RORGN is yearly input from roots and stubble N from harvest to
C harvest.
C In initial year, assume that previous stubble input, SXORGN, (i.e.
C from two crops before) is same as current stubble input.
C
        RORGN=XORGN+SXORGN
C
C Set rate of return of C (CRATE) and N (ANRATE) from dead roots etc.
C
        MCROP=LCROP
C
C Set the end of the growing season
C
        CALL SETEND(IYEAR,IEND,ISTHARV,N_STEPS,ISOWN,MEND,NSOW)
C
C Set degree-days to 2000 at anthesis of initial year
C
        DDAYS=2000
      ENDIF
C
C Set ready for senescence
C
      JSTOP=0
C
C Save crop characteristics
C
      ISAVE=1
      CALL SAVE_CROP(MCROP,CAO,UR,UT,INC,CRATE,ANRATE,CFACT,SEED,
     &                     C1,T1,RRG,RRX,CSC,STRAW,IROCKS,ROOTS,CREQN,
     &                     TOTC,ISTART,CINP,ANINP,TCNEED,BREQN,
     &                     CONVER_F,N_STEPS,N_REMAIN,
     &                     ISAVE)
C
C Leave INIT_CROP
C
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE RUN1_SUNDIAL_CROP(IYEAR,IK,MEND,IS_TS,JSTOP,
     &                             IEND,IANTHES,I_TS,
     &                             NSOW,ISTHARV,INSTRAW,ICROP,LCROP,
     &                             ISOWN,OCROPN,YLD,PREYLD,EXYLD,C_TS,
     &                             SXORGN,SORGN,TORGC,SORGC,RORGN,ORGC,
     &                             SXORGN15,SXORGC,XORGC,XORGN,ORGN,
     &                             TCINP,TRNIN,RNIN,RNIN15,
     &                             WR,CTOT,
     &                             CACTOT,CATOT15,ICOVER,
C required for SWAT water
     &				             PLANTBIOMASS,PLANTHEIGHT,LAI,
     &							 CANMAX,CULTIVATE,NCULT,N_TS)      
C
C Subroutine to Calculate the total C returns for next crop and set up crop at sowing
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      INTEGER MAXCROP
	PARAMETER (MAXCROP = 36)
	INTEGER MCROP,IROCKS(0:MAXCROP),ISTART,ISAVE
	REAL CINP,ANINP,ROOTS,CREQN,TOTC,TCNEED,BREQN 
      REAL CAO(5,0:MAXCROP),UR(3,0:MAXCROP),UT(3,0:MAXCROP) 
	REAL INC(2,0:MAXCROP),CRATE(0:MAXCROP),ANRATE(0:MAXCROP)
      REAL CFACT(0:MAXCROP),SEED(0:MAXCROP),C1(0:MAXCROP),T1(0:MAXCROP)
      REAL RRG(0:MAXCROP),RRX(0:MAXCROP)
      REAL CSC(4,0:MAXCROP),STRAW(4,0:MAXCROP)
C
C Variables passed to/from calling subroutine
C
      INTEGER CULTIVATE			! IN(CALL): Code to cultivate this week (1=yes, 0=no)
	INTEGER NCULT				! IN(CALL): Number of cultivations
	INTEGER IYEAR,IK,NSOW,ISTHARV,IEND,MEND
	INTEGER IS_TS,INSTRAW,ICROP
	INTEGER LCROP,N_STEPS,ISOWN,JSTOP,N_REMAIN
      INTEGER IANTHES,I_TS,N_TS
	REAL OCROPN,SXORGN,SORGN,CACTOT,CATOT15,TORGC,SORGC,RORGN,RNIN15
      REAL WR,YLD,PREYLD,CTOT,EXYLD,ORGC
      REAL RNIN,CONVER_F,TRNIN,SXORGN15,C_TS,SXORGC,TCINP
	REAL XORGN,ORGN,XORGC
	INTEGER ICOVER
C
C Variables required for SWAT water routines
C
	REAL PLANTBIOMASS		! Plant biomass at harvest [ kg ha-1]
	REAL PLANTHEIGHT		! Plant height at harvest [cm]
	REAL LAI				! Plant LAI at harvest
	REAL canMAX				! maximum amount of water that can be trapped in the canopy when the canopy is fully developed (mm H2O)
	REAL aVP					! IN: water vapor pressure of air at height z (kPa)
	REAL z						! elevation above sea level [m]
C								  in the canopy when the canopy is fully developed (mm H2O)
	REAL CURRENTPLBIOM			! Current plant biomass [ kg ha-1]
	REAL CURRENTPLHEI			! Current plant height [cm]
	REAL CURRENTLAI				! Current LAI
	REAL canStor				! amount of free water held in the canopy on a given day (mm H2O)
	INTEGER ROOTLAYER
C
C
C Retrieve crop characteristics
C
      ISAVE=0
      CALL SAVE_CROP(MCROP,CAO,UR,UT,INC,CRATE,ANRATE,CFACT,SEED,
     &                     C1,T1,RRG,RRX,CSC,STRAW,IROCKS,ROOTS,CREQN,
     &                     TOTC,ISTART,CINP,ANINP,TCNEED,BREQN,
     &                     CONVER_F,N_STEPS,N_REMAIN,
     &                     ISAVE)
C
C Count number of weeks since sowing
C
      CALL COUNTW(IYEAR,IK,ISTART,N_STEPS,ISOWN,IANTHES,I_TS,NSOW)
C
C Add any crop debris
C 
      CALL RESIDU(I_TS,ISTART,IEND,IYEAR,CINP,ANINP,IK,
     &                  CACTOT,CREQN,TRNIN,RNIN,RNIN15,CATOT15,
     &                  ICROP,SXORGN,SXORGN15,C_TS,SXORGC,TCINP,
     &                  XORGC,XORGN,ORGC,ORGN,CRATE,ANRATE,MCROP,
     &                  CONVER_F,CULTIVATE,NCULT,N_TS)
C
C In week of sowing set up parameters for the new sown crop
C
 
      IF(IYEAR.GT.0.AND.IK.EQ.NSOW)THEN
        CALL SOWING(IYEAR,INSTRAW,MCROP,LCROP,ICROP,IEND,ISTHARV,
     &              N_STEPS, ISOWN, MEND, NSOW, IS_TS, JSTOP,
     &              OCROPN, SXORGN, SORGN, YLD, CREQN, PREYLD,
     &              TORGC, CTOT, EXYLD, TOTC, SORGC, RORGN, ORGC,
     &              CACTOT, CATOT15, TCNEED, BREQN, UR,
     &              UT, SEED, INC, CFACT, XORGN, ORGN, ROOTS,STRAW)
C
C Extra crop parameters required for SWAT water (not interferring with main programme)
C
	 CALL SETEXTRACROPPARAM(YLD,PLANTBIOMASS,PLANTHEIGHT,LAI,canMax)
C
       END IF
C
C Cycle root N from total crop N back to RO pool
C  N and C is cycled back from dead roots to RO pool, starting
C in the week after crop uptake starts (i.e. IK.GT.NSOW)
C 
      IF(ICROP.NE.0)THEN
        IF(IYEAR.GT.0.AND.IK.GT.NSOW.AND.IK.LE.MEND)THEN
          CALL RCYCLE(MEND, IK, N_REMAIN, ICROP, 
     &                  WR, CACTOT, CATOT15, RNIN, RNIN15, CONVER_F)
        ENDIF
      ENDIF
C
C Work out if there is crop cover or not
C
	IF(IYEAR.GT.0.AND.IK.GT.NSOW.AND.IK.LE.MEND)THEN
		ICOVER=1
      ELSE IF (IYEAR==0) THEN
		ICOVER=1
	ELSE
		ICOVER=0
      ENDIF
C
C Save crop characteristics
C
      ISAVE=1
      CALL SAVE_CROP(MCROP,CAO,UR,UT,INC,CRATE,ANRATE,CFACT,SEED,
     &                     C1,T1,RRG,RRX,CSC,STRAW,IROCKS,ROOTS,CREQN,
     &                     TOTC,ISTART,CINP,ANINP,TCNEED,BREQN,
     &                     CONVER_F,N_STEPS,N_REMAIN,
     &                     ISAVE)
C
C Leave RUN1_SUNDIAL_CROP
C
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE RUN2_SUNDIAL_CROP(IYEAR,IK,NSOW,ICROP,JSTOP,MEND,IS_TS,
     &                             NSOIL,IL_TSS,AIRTEMP,DDAYS,CACT,
     &                             CACT15,ATM,ATM15,CACTOT,TC,
     &                             L_TSS,TACTOT,SRNIN,TRNIN,
     &                             RNIN,OCROPN,CLOSSX,CLOSSX15,RORGN,
     &                             VOLAT,VOLAT15,CATOT15,CUPTN,SOILN, 
     &                             TIN,SOIL15,AMMN15,AMMN,TAM,CRIT,CTOT,
     &                             SEEDIN,XORGN,
C Variables required for SWAT water
     &				             PLANTBIOMASS,PLANTHEIGHT,LAI,
     &							 CURRENTPLBIOM,CURRENTPLHEI,CURRENTLAI,
     &							 ROOTLAYER,PLANTUP)	
C
C Subroutine to take up N and return C and N during season
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER
	PARAMETER (MAXLAYER=60)
      INTEGER MAXCROP,MAXSOIL
	PARAMETER (MAXCROP=36)
	PARAMETER (MAXSOIL=50)
C
C Variables local to subroutine
C

C
C Variables passed to/from this subroutine
C
	REAL AIRTEMP
	REAL AMMN15(MAXLAYER)
	REAL AMMN(MAXLAYER)
	REAL ANINP
	REAL ANRATE(0:MAXCROP)
	REAL ATM
	REAL ATM15
	REAL BREQN 
	REAL C1(0:MAXCROP)
      REAL CAO(5,0:MAXCROP)
	REAL CACT
	REAL CACT15
	REAL CACTOT
	REAL CATOT15
      REAL CFACT(0:MAXCROP)
	REAL CINP
      REAL CLOSSX
	REAL CLOSSX15
      REAL CONVER_F
	REAL CRATE(0:MAXCROP)
	REAL CREQN
	REAL CRIT(MAXSOIL,MAXLAYER)
      REAL CSC(4,0:MAXCROP)
	REAL CTOT
	REAL CUPTN
	REAL CURRENTLAI			! Current LAI
	REAL CURRENTPLBIOM		! Current plant biomass [ kg ha-1]
	REAL CURRENTPLHEI		! Current plant height [cm]
	REAL DDAYS
	INTEGER ICROP
	INTEGER IK
	INTEGER IL_TSS
	REAL INC(2,0:MAXCROP)
	INTEGER IROCKS(0:MAXCROP)
	INTEGER IS_TS
	INTEGER ISTART
	INTEGER ISAVE
	INTEGER IYEAR
	INTEGER JSTOP
	INTEGER L_TSS(0:MAXCROP)
	REAL LAI				! Plant LAI at harvest
	INTEGER MEND
	INTEGER MCROP
	INTEGER N_REMAIN
	INTEGER N_STEPS
	INTEGER NSOIL
	INTEGER NSOW
	REAL OCROPN
	REAL PLANTBIOMASS		! Plant biomass at harvest [ kg ha-1]
	REAL PLANTHEIGHT		! Plant height at harvest [cm]
	REAL PLANTUP(MAXLAYER)	! IN(CROPGROW)/OUT:Plant uptake from each soil layer
	INTEGER ROOTLAYER
	REAL ROOTS
      REAL RRG(0:MAXCROP)
	REAL RRX(0:MAXCROP)
	REAL SEED(0:MAXCROP)
	REAL SEEDIN
	REAL SOIL15(MAXLAYER)
      REAL SOILN(MAXLAYER)
	REAL SRNIN
	REAL STRAW(4,0:MAXCROP)
	REAL T1(0:MAXCROP)
	REAL TACTOT
	REAL TAM(MAXLAYER)
	REAL TC
	REAL TCNEED
	REAL TIN(MAXLAYER)
	REAL TOTC
	REAL TRNIN
	REAL RNIN
	REAL RORGN
	REAL UR(3,0:MAXCROP)
	REAL UT(3,0:MAXCROP) 
	REAL VOLAT
	REAL VOLAT15
	REAL XORGN
C
C Retrieve crop characteristics
C
      ISAVE=0
      CALL SAVE_CROP(MCROP,CAO,UR,UT,INC,CRATE,ANRATE,CFACT,SEED,
     &                     C1,T1,RRG,RRX,CSC,STRAW,IROCKS,ROOTS,CREQN,
     &                     TOTC,ISTART,CINP,ANINP,TCNEED,BREQN,
     &                     CONVER_F,N_STEPS,N_REMAIN,
     &                     ISAVE)

C
C Caclulate the number of degree days since sowing
C
	CALL DAYD(IYEAR,IK,NSOW,DDAYS,AIRTEMP,CONVER_F) 
C
C If freezing or not in growing season set crop uptake parameters to zero
C
      IF(IYEAR.EQ.0.OR.ICROP.EQ.0.OR.AIRTEMP.LE.0.0.OR.
     &  (IYEAR.GT.0.AND.IK.LT.NSOW).OR.JSTOP.EQ.1)THEN
        CALL ZEROUP(CUPTN,CACT,CACT15)
C
C Otherwise calculate the N/C turnover associated with growing the crop
C
      ELSE
        CALL CROPGROW(MCROP,IYEAR,IK,NSOW,MEND,IS_TS,NSOIL, 
     &                CACT15,SOILN,TIN,SOIL15,AMMN15,AMMN,TAM,
     &                DDAYS,CUPTN,ROOTS,CONVER_F,ATM,ATM15, 
     &                T1,C1,RRG,IROCKS,RRX,CRIT,CREQN,CACT,CTOT,TCNEED,
     &				PLANTBIOMASS,PLANTHEIGHT,LAI,
     &				CURRENTPLBIOM,CURRENTPLHEI,CURRENTLAI,
     &				ROOTLAYER,PLANTUP)
      END IF
C
C In sowing week add seed N
C
      IF(IYEAR.GT.0.AND.IK.EQ.NSOW)THEN
	  CACT=CACT+SEED(MCROP)
	  SEEDIN=SEED(MCROP)
	ELSE
	  SEEDIN=0
	ENDIF
C
C Sum N in crop
C
      TC=TC+CACT
C
C l_tss weeks before harvest....
C

      IF(ICROP.NE.0)THEN
        IL_TSS = NINT(L_TSS(MCROP) * (1 / CONVER_F))
        IF(JSTOP.EQ.0.AND.MEND-IK.EQ.IL_TSS-1)THEN
C
C UPdate this week's crop uptake total
C
          TACTOT=CACTOT+CACT
          SRNIN=TRNIN+RNIN
C
C start senescence
C
          CALL SENES(IYEAR, MCROP, JSTOP,OCROPN,SRNIN,RORGN,TACTOT,
     &                 CONVER_F,CLOSSX,L_TSS, XORGN)
        ENDIF
C
C Required for SWAT water
C
	  IF(JSTOP.EQ.1.AND.MEND-IK.LE.IL_TSS-1)THEN
C
C Calculate current status of extra crop parameters (biomass,height,LAI)
C
		CALL CURRENTPLANTPARAM(PLANTBIOMASS,PLANTHEIGHT,LAI, 
     &							C1,T1,DDAYS,MCROP,
     &							CURRENTPLBIOM,CURRENTPLHEI,CURRENTLAI)
	 END IF
C
C Set current crop parameters to ZERO
C
	  IF (JSTOP.EQ.1.AND.IK.EQ.MEND) THEN
			CURRENTPLBIOM=0.
			CURRENTPLHEI =0.
			CURRENTLAI   =0.
			ROOTLAYER	 =0.
	  END IF
C
C End of requirement for SWAT water
C
C If stopping take off losses due to senescence
C
        IF(JSTOP.EQ.1)THEN
          CLOSSX15=CLOSSX*CATOT15/CACTOT
          CACTOT=CACTOT-CLOSSX
          CATOT15=CATOT15-CLOSSX15
          VOLAT=VOLAT+CLOSSX
          VOLAT15=VOLAT15+CLOSSX15
        ENDIF
      ENDIF
C
C Save crop characteristics
C
      ISAVE=1
      CALL SAVE_CROP(MCROP,CAO,UR,UT,INC,CRATE,ANRATE,CFACT,SEED,
     &                     C1,T1,RRG,RRX,CSC,STRAW,IROCKS,ROOTS,CREQN,
     &                     TOTC,ISTART,CINP,ANINP,TCNEED,BREQN,
     &                     CONVER_F,N_STEPS,N_REMAIN,
     &                     ISAVE)
      CALL GETCAC(CACTOT,CACT,CATOT15,CACT15)
C
C Leave RUN2_SUNDIAL_CROP
C
      END
C
C-----------------------------------------------------------------------
C
C INTERNAL SUBROUTINES
C
C------------------------------------------------------------------
C
      SUBROUTINE ADDCARB(CINP,RNIN,CACTOT,RNIN15,CATOT15,ANINP,ICROP,
     &                   C_TS)
C
C Subroutine to add plant rsidue additions from the growing crop
C
      IMPLICIT NONE
C
C Variable local to this subroutine
C
C
C Variables passed to/from calling subroutine
C
      INTEGER ICROP				! IN:Current crop number
	REAL ANINP					! IN:N added to soil as non-cartable deb. 
								!        in this timestep (kgN/ha/timestep)
	REAL CINP					! IN:C added to soil as non-cartable deb.
								!        in this timestep (kgN/ha/timestep)
	REAL RNIN					! OUT:Deb.N input in this timestep (kgN/ha)
	REAL RNIN15					! OUT:Deb.N15 input in timestep (kgN15/ha)
	REAL CACTOT					! IN:N taken up by crop (kgN/ha)
	REAL CATOT15				! IN:N15 taken up by crop (kgN15/ha)
	REAL C_TS					! OUT:Deb.C input in this timestep (kgC/ha)
C
C  Set C,N and N15 additions
C
	C_TS=CINP
      RNIN=ANINP
      IF(CACTOT.GT.0)THEN
        RNIN15=RNIN*CATOT15/CACTOT
C      ELSEIF(ICROP.GT.0)THEN
C	  WRITE(15,*)'Error ! Unexpected zero crop offtake in ADDCARB'
      ENDIF
C
C Leave ADDCARB
C
      RETURN
      END
C
C--------------------------------------------------------------
C
      SUBROUTINE ADDUP(BREQN, CACT, CREQN, TCNEED)
C
C Subroutine to the crop N requirement for this week
C
C Set the crop N requirement for this week
C
      IMPLICIT NONE

C REAL passed from subroutine
      REAL BREQN, CACT, CREQN, TCNEED

      BREQN=BREQN-CACT
      CREQN=CREQN+TCNEED
      IF(CREQN.LE.0)CREQN=0.0001
      TCNEED=0
C
C Leave ADDUP
C
      RETURN
      END
C
C----------------------------------------------------------
C
      SUBROUTINE ATMFIX(CACT, CFIXNEED, ATM, CACT15, ATM15, FIXN)
C
C Subroutine to fix nitrogen from the atmosphere
C
      IMPLICIT NONE
C
C REAL passed from subroutine
      REAL CACT, CFIXNEED, ATM, CACT15, ATM15, FIXN 
C
C Fix required nitrogen from the atmosphere
C
      CACT=CACT+CFIXNEED
      IF(ATM.GT.0)CACT15=CACT15+(CFIXNEED*ATM15/ATM)
C
C Record how much N is fixed for weekly flow chart and balance sheets
C
      FIXN=CFIXNEED
C
C Leave atmfix
C
      END
C
C-----------------------------------------------------------
C
      SUBROUTINE CROPGROW(MCROP,IYEAR,IK,NSOW,MEND,IS_TS,NSOIL, 
     &                    CACT15,SOILN,TIN,SOIL15,AMMN15,AMMN,TAM,
     &                    DDAYS,CUPTN,ROOTS,CONVER_F,ATM,ATM15,
     &                    T1,C1,RRG,IROCKS,RRX,CRIT,CREQN,CACT,CTOT,
     &                    TCNEED,
     &				PLANTBIOMASS,PLANTHEIGHT,LAI,
     &				CURRENTPLBIOM,CURRENTPLHEI,CURRENTLAI,
     &				ROOTLAYER,PLANTUP)
C
C Subroutine to calculate the C/N turnover associated with the growing crop
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER
	PARAMETER (MAXLAYER=60)
      INTEGER MAXCROP, MAXSOIL
	PARAMETER (MAXCROP=36)
	PARAMETER (MAXSOIL=50)
	REAL CNEED
	REAL CFIXNEED
      REAL FIXFRAC(0:MAXCROP)
      INTEGER M
C
C Variables passed to/from calling subroutine
C
	REAL AMMN15(MAXLAYER)
	REAL AMMN(MAXLAYER)
	REAL ATM
	REAL ATM15  
	REAL BREQN
	REAL C1(0:MAXCROP)
	REAL CACT
	REAL CACT15
	REAL CACTOT
	REAL CONVER_F
	REAL CREQN
	REAL CRIT(MAXSOIL,MAXLAYER)
	REAL CUPTN
	REAL CURRENTLAI			! Current LAI
	REAL CURRENTPLBIOM		! Current plant biomass [ kg ha-1]
	REAL CURRENTPLHEI		! Current plant height [cm]
	REAL DDAYS
	REAL CTOT
	REAL FIXN
      INTEGER IROCKS(0:MAXCROP)
	INTEGER IS_TS
	INTEGER IYEAR
	INTEGER IK
	REAL LAI				! Plant LAI at harvest
	INTEGER MEND
	INTEGER MCROP
	INTEGER NSOIL
	INTEGER NSOW
	REAL PLANTBIOMASS		! Plant biomass at harvest [ kg ha-1]
	REAL PLANTHEIGHT		! Plant height at harvest [cm]
	REAL PLANTUP(MAXLAYER)	! IN(EXTRACT):OUT:Plant uptake from each soil layer
	INTEGER ROOTLAYER		! Layer until which roots have grown
	REAL ROOTS
	REAL RRG(0:MAXCROP)
      REAL RRX(0:MAXCROP)
	REAL SOIL15(MAXLAYER)
      REAL SOILN(MAXLAYER)
      REAL T1(0:MAXCROP)
	REAL TAM(MAXLAYER)
	REAL TCNEED
	REAL TIN(MAXLAYER)

C
C Set FIXFRAC (REPLACE WITH FIXFRAC READ FROM FILE AND BLOCK DATA
C
      FIXFRAC(MCROP)=0
c      IF(MCROP.EQ.4.OR.MCROP.EQ.8.OR.MCROP.EQ.13.OR.MCROP.EQ.14
c     &   .OR.MCROP.EQ.20.OR.MCROP.EQ.23)THEN
c        FIXFRAC(MCROP)=0.8
c      ENDIF
C
C In first year or in growing season set the crop N requirement for this week
C
      IF(IYEAR.EQ.0.OR.(IYEAR.GT.0.AND.IK.GE.NSOW.AND.
     &   IK.LE.MEND))CALL ADDUP(BREQN,CACT,CREQN,TCNEED)
C
C In first year or first week of the growing season
C set the inital rate of crop N uptake
C
      IF(IYEAR.EQ.0.AND.IK.EQ.1)CALL INITUP(CTOT,CREQN,BREQN,CACTOT)
C
C Calculate crop N uptake
C
      CALL UPTAKE(DDAYS,CTOT,T1,C1,CREQN,MCROP,CUPTN)
      CNEED=CUPTN
      IS_TS=IS_TS+1
C
C Required for SWAT water
C Calculate current status of extra crop parameters (biomass,height,LAI)
C
	CALL CURRENTPLANTPARAM(PLANTBIOMASS,PLANTHEIGHT,LAI, 
     &							C1,T1,DDAYS,MCROP,
     &							CURRENTPLBIOM,CURRENTPLHEI,CURRENTLAI)
C
C Calc. root depth at 5cm/week since sowing
C
      ROOTS=IS_TS*RRG(MCROP)*CONVER_F
      CALL ROOTL(M,MCROP,IROCKS,ROOTS,RRX)
	ROOTLAYER=M				!required for SWAT water
C
C Proportion of N15 in SOILN and AMMN remains same for each layer
C after extraction by roots and crop. Minimum value of SOILN and
C AMMN prevents N15 going negative.
C
      CUPTN=(1-FIXFRAC(MCROP))*CUPTN
      CALL EXTRACT(M,CACT,NSOIL,CUPTN,SOILN,CRIT,AMMN,PLANTUP)
      CALL EXTR15(CACT15,SOILN,TIN,SOIL15,AMMN15,AMMN,TAM)
      CFIXNEED=CNEED-CUPTN
C
C Fix nitrogen from the atmosphere for legumes
C
      CALL ATMFIX(CACT,CFIXNEED,ATM,CACT15,ATM15,FIXN)
C
C Sum the total N still required by the crop
C
      TCNEED=TCNEED+(CNEED-CACT)
C
C Leave CROPGROW
C
      RETURN
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE CROPN(CREQN,YLD, UT, MCROP)
C
C Calculate crop N requirement
C
C
C Parameters
C
      IMPLICIT NONE

      INTEGER MAXCROP, MCROP

	PARAMETER (MAXCROP=36)

C REAL passed from subroutine
      REAL CREQN, YLD, UT(3,0:MAXCROP)
C
C Set tops N uptake
C
       CREQN=UT(1,MCROP)*(EXP(UT(2,MCROP)*YLD)-UT(3,MCROP))
C
C Leave CROPN
C
       RETURN
       END
C
C-----------------------------------------------------------------
C
      SUBROUTINE COUNTW(IYEAR,IK,ISTART,N_STEPS,ISOWN,IANTHES,I_TS,NSOW)
C
C Subroutine to count the number of weeks since sowing, I_TS
C
      IMPLICIT NONE
C
C Variables passed to/from calling subroutine
C
      INTEGER IYEAR,IK,ISTART,N_STEPS,ISOWN,IANTHES,I_TS,NSOW
C
C Count time steps
C
      IF(IYEAR.EQ.0.AND.IK.EQ.1)THEN
       ISTART=(N_STEPS-ISOWN)+IANTHES
       I_TS=ISTART
      ELSE IF(IYEAR.GT.0.AND.IK.EQ.NSOW)THEN
       ISTART=1
       I_TS=ISTART
      ELSE
       I_TS = I_TS + 1
      END IF
C
C Leave COUNTW
C
      RETURN
      END
C
C------------------------------------------------------------
C
	SUBROUTINE CURRENTPLANTPARAM(PLANTBIOMASS,PLANTHEIGHT,LAI, 
     &							C1,T1,DDAYS,MCROP,
     &							CURRENTPLBIOM,CURRENTPLHEI,CURRENTLAI)
C
C Subroutine to calcluate current plant height, biomass and LAI according the 
C cumulative N uptake
C
	IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXCROP				! Max.no.of crop types allows
      PARAMETER (MAXCROP=36)
C
C Variables passed to calling subroutine
C
	REAL PLANTBIOMASS		! Plant biomass at harvest [ kg ha-1]
	REAL PLANTHEIGHT		! Plant height at harvest [cm]
	REAL LAI				! Plant LAI at harvest
	REAL C1(0:MAXCROP)
	REAL T1(0:MAXCROP)
	REAL DDAYS
	INTEGER MCROP

C
C Variables passed from calling subroutine
C
	REAL CURRENTPLBIOM		! Current plant biomass [ kg ha-1]
	REAL CURRENTPLHEI		! Current plant height [cm]
	REAL CURRENTLAI			! Current LAI
C
	IF(T1(MCROP).GT.0)THEN
      CURRENTPLBIOM=((100**(-1/T1(MCROP))+
     &			  EXP(-C1(MCROP)*DDAYS))**(-T1(MCROP)))/
     &              100*PLANTBIOMASS
	CURRENTPLHEI=((100**(-1/T1(MCROP))+
     &			 EXP(-C1(MCROP)*DDAYS))**(-T1(MCROP)))/100*PLANTHEIGHT
	CURRENTLAI=((100**(-1/T1(MCROP))+
     &			EXP(-C1(MCROP)*DDAYS))**(-T1(MCROP)))/100*LAI
      ELSE
        WRITE(*,*)'Error! Unexpected zero in crop parameters
     &									 (biomass/height/LAI)'
        WRITE(15,*)'Error! Unexpected zero in crop parameters
     &									 (biomass/height/LAI)'
      ENDIF
C
C Leave CURRENTPLANTPARAM
C
      RETURN
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE DAYD(IYEAR,IK,NSOW,DDAYS,AIRTEMP,CONVER_F)
C
C Subroutine to calculate the number of degree days
C
      IMPLICIT NONE
C
C Variables passed to/from calling subroutine
C
	REAL AIRTEMP				! IN:Air temperature (deg.C)
	REAL CONVER_F				! IN:Conversion between this timestep & weeks
      REAL DDAYS					! OUT:Degree days since sowing (deg.C)
	INTEGER IK					! IN:No.timesteps from prev.harv. to current
	INTEGER IYEAR				! IN:Current growing season number
	INTEGER NSOW				! IN:Timesteps from prev.harv. to sowing date
C
C Calc. the number of degree days since sowing
C
      IF(IYEAR.GT.0.AND.IK.EQ.NSOW)DDAYS=0
      IF(AIRTEMP.GT.0.0)DDAYS=DDAYS+AIRTEMP*7*CONVER_F
C
C Leave DAYD
C
      RETURN
      END
C
C------------------------------------------------------------------
C
      SUBROUTINE EXTRACT(MC,ACTUPT,NSOIL,CUPTN,SOILN,CRIT,AMMN,PLANTUP)
C
C Subroutine to extract mineral nitrogen required by the crop
C from the soil profile
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER MAXSOIL				! Max.no.of soil types
	PARAMETER (MAXSOIL=50)

	REAL ACRIT
	REAL AUPTN
      INTEGER IN
	INTEGER IUPT
	INTEGER N
      REAL SCRIT
	REAL UPTN
	REAL UPTN1
C
C Variables passed to/from calling subroutine
C
	REAL ACTUPT					! IN:Crop N uptake (kgN/ha/timestep)
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL CRIT(MAXSOIL,MAXLAYER)	! IN:Critical N content (kgN/ha/layer)
	REAL CUPTN					! IN:Crop N uptake (kgN/ha/timestep)
      INTEGER MC					! IN:Bottom layer of roots
	INTEGER NSOIL				! IN:Soil code number
	REAL PLANTUP(MAXLAYER)		! OUT:Plant uptake from each soil layer
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
C
C Crop takes up NO3 and NH4 in proportion, depletes layer of
C both before taking from lower layer
C
      IN=0
      ACTUPT=0
      SCRIT=0
      IUPT=0
      UPTN=CUPTN
      DO 44 N=1,MC
        IN=0
        IF(IUPT.EQ.0)THEN
C
C SCRIT = available NO3. ACRIT = available NH4
C
          SCRIT=SOILN(N)-(CRIT(NSOIL,N)*(MAXDEPTH/(MAXLAYER1*50.)))
          ACRIT=AMMN(N)-(CRIT(NSOIL,N)*(MAXDEPTH/(MAXLAYER1*50.)))
C
C if uptake exceeds available N, take all NH4 , set IN = 1 (take all
C NO3: see below)
C
          IF(UPTN.GE.(SCRIT+ACRIT))THEN
            ACTUPT=ACTUPT+SCRIT+ACRIT
	      PLANTUP(N)=SCRIT+ACRIT
            AMMN(N)=CRIT(NSOIL,N)*(MAXDEPTH/(MAXLAYER1*50.))
	      SOILN(N)=CRIT(NSOIL,N)*(MAXDEPTH/(MAXLAYER1*50.))
            UPTN=UPTN-(SCRIT+ACRIT)
            IN=1
C
C if uptake less than available N, take NH4 and NO3 in proprtion to NH4+NO3
C
          ELSE
            AUPTN=UPTN*(ACRIT/(ACRIT+SCRIT))
            UPTN1=UPTN*(SCRIT/(SCRIT+ACRIT))
            ACTUPT=ACTUPT+AUPTN+UPTN1
            AMMN(N)=AMMN(N)-AUPTN
	      SOILN(N)=SOILN(N)-UPTN1
	      PLANTUP(N)=AUPTN+UPTN1
C
C stop loop if N requirement has been extracted
C
            IUPT=1
          ENDIF
        END IF
   44 CONTINUE
C
C Leave EXTRACT
C
      RETURN
      END
C
C-------------------------------------------------------------------------
C
      SUBROUTINE EXTR15(CACT15,SOILN,TIN,SOIL15,AMMN15,AMMN,TAM)
C
C Subroutine to extract Min-N15 required by the crop from the soil profile
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER
	INTEGER MAXLAYER1
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/
      INTEGER IL
C
C Variables passed to/from calling subroutine
C
      REAL CACT15,SOILN(MAXLAYER),TIN(MAXLAYER),SOIL15(MAXLAYER)
	REAL AMMN15(MAXLAYER),AMMN(MAXLAYER),TAM(MAXLAYER)
C
C Extract 15N in proportion to uptake of total NO3 and NH4
C
      CACT15=0
      DO 100 IL=1,MAXLAYER1
        CACT15=CACT15+(SOIL15(IL)-SOILN(IL)*TIN(IL))
        SOIL15(IL)=SOILN(IL)*TIN(IL)
        CACT15=CACT15+(AMMN15(IL)-AMMN(IL)*TAM(IL))
        AMMN15(IL)=AMMN(IL)*TAM(IL)
  100 CONTINUE
C
C Leave EXTR15
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE GETCAC(CACTOT,CACT,CATOT15,CACT15)
C
C Subroutine to get crop uptake
C
      CACTOT=CACTOT+CACT
      CATOT15=CATOT15+CACT15
      END
C
C----------------------------------------------------
C
C
      SUBROUTINE INITUP(CTOT, CREQN, BREQN, CACTOT)
C
C Subroutine to set the initial rate of crop N uptake
C
      IMPLICIT NONE

C REAL passed from subroutine
      REAL CTOT, CREQN, BREQN, CACTOT
C
C CTOT=uptake to the current week
C Set to the maximum crop requirement
C
      CTOT=0
      BREQN=CREQN-CTOT
      CACTOT=CTOT
C
C Leave INITUP
C
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MOREFIX(IYEAR,MCROP,LCROP,ICROP,INSTRAW,
     &                   YLD,PREYLD,EXYLD,SORGC,SORGN, 
     &                   SXORGC,SXORGN,HZ1,TORGC,TORGN,TXORGC, 
     &                   TXORGN,CAO,CSC,STRAW,
     &                   ORGC,ORGN,XORGC,XORGN,
C introduced for SWAT water routines
     &				             PLANTBIOMASS,PLANTHEIGHT,LAI,C1,T1,
     &							 canMax,CONVER_F)
C
C Subroutine to set parameters dependant on inputs
C from sites
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      INTEGER MAXCROP,MAXSOIL
	PARAMETER (MAXCROP=36)
      PARAMETER (MAXSOIL=50)
	INTEGER MAXLAYER
	PARAMETER (MAXLAYER=60)
C
C Variables passed from calling subroutine
C
      INTEGER IYEAR, MCROP, LCROP, ICROP, INSTRAW
      REAL YLD, PREYLD, EXYLD, SORGC, SORGN, SXORGC, SXORGN
	REAL HZ1(MAXLAYER), TORGC, TORGN, TXORGC, TXORGN
	REAL CAO(5,0:MAXCROP), CSC(4,0:MAXCROP), STRAW(4,0:MAXCROP)
	REAL ORGC, ORGN, XORGC, XORGN
C
C Variables required for SWAT water routines
C
	REAL PLANTBIOMASS		! Plant biomass at harvest [ kg ha-1]
	REAL PLANTHEIGHT		! Plant height at harvest [cm]
	REAL LAI				! Plant LAI at harvest
	REAL canMAX				! maximum amount of water that can be trapped in the canopy when the canopy is fully developed (mm H2O)
	REAL aVP					! IN: water vapor pressure of air at height z (kPa)
	REAL z						! elevation above sea level [m]
C								  in the canopy when the canopy is fully developed (mm H2O)
	REAL CURRENTPLBIOM			! Current plant biomass [ kg ha-1]
	REAL CURRENTPLHEI			! Current plant height [cm]
	REAL CURRENTLAI				! Current LAI
	REAL canStor				! amount of free water held in the canopy on a given day (mm H2O)
	REAL C1(0:MAXCROP)
	REAL T1(0:MAXCROP)
	INTEGER ROOTLAYER
	REAL SOILFC(MAXLAYER)		! IN: soil water content at FC (mm/layer)
	REAL CONVER_F				! IN:Conversion between timestep & weeks

C
C In first year calc Total and Stubble Carbon return from initial crop
C
      IF(IYEAR.LE.1)THEN
       MCROP=LCROP
       YLD=PREYLD
       CALL SETC(INSTRAW,MCROP,TXORGC,SXORGC,YLD,CAO,CSC,STRAW)
C
C Extra crop parameters required for SWAT water (not interferring with main programme)
C
	 CALL SETEXTRACROPPARAM(YLD,PLANTBIOMASS,PLANTHEIGHT,LAI,canMax)
C
C
C In subsequent years...
C
      ELSE
C
C At Anthesis, growing crop's C and N relabelled as they will
C be input as stubble at harvest.
C
       XORGN=ORGN
       XORGC=ORGC
       SXORGN=SORGN
       SXORGC=SORGC
      END IF
C
C Calc total C returns for next crop
C
      MCROP=ICROP
      YLD=EXYLD
      CALL SETC(INSTRAW,MCROP,TORGC,SORGC,YLD,CAO,CSC,STRAW)
C
C Calc N inputs for previous and current year
C
      IF(IYEAR.EQ.0)THEN
       TXORGN=HZ1(1)*TXORGC
      ELSE
       TXORGN=TORGN
      END IF
      TORGN=HZ1(1)*TORGC
C
C Leave MOREFIX
C
      RETURN
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE PAR_CROP(L_TSS, IROCKS, CAO, UR, UT, INC, CRATE, 
     &                    ANRATE, CFACT, SEED, C1, T1, RRG, RRX,  
     &                    CSC, STRAW)
C
C Subroutine to reset parameters from files PARAM.OUT and PARAM.DAT
C Parameters
C
      IMPLICIT NONE
      INTEGER MAXCROP
	PARAMETER (MAXCROP=36)

C INTEGER passed from subroutine    
	INTEGER L_TSS(0:MAXCROP),IROCKS(0:MAXCROP)

C REAL passed from subroutine
      REAL CAO(5,0:MAXCROP),UR(3,0:MAXCROP),UT(3,0:MAXCROP) 
	REAL INC(2,0:MAXCROP),CRATE(0:MAXCROP),ANRATE(0:MAXCROP)
      REAL CFACT(0:MAXCROP),SEED(0:MAXCROP),C1(0:MAXCROP),T1(0:MAXCROP)
      REAL RRG(0:MAXCROP),RRX(0:MAXCROP)
      REAL CSC(4,0:MAXCROP),STRAW(4,0:MAXCROP)

C INTEGER local to subroutine
      INTEGER NUMCROP, I, J

C CHARACTER local to subroutine (KWC: I think cname is only in crop_par)
      CHARACTER*40 TEMP, CNAME(0:MAXCROP)
C
C Set default array descriptors
C
      I=1
      J=1
C
C Set default number of crops to number of crops previously parameterised
C
      NUMCROP=4
C
C Try to open CROP_SUN.DAT
C
      OPEN(43,FILE='CROP_SUN.DAT',STATUS='OLD',ERR=111)
C
C Read in parameters from CROP_SUN.DAT
C
101   CONTINUE
      READ(43,9,ERR=111,END=111)TEMP,I
      READ(43,10,ERR=112,END=111)(CAO(J,I),J=1,5),
     &                          (UR(J,I),J=1,3),(UT(J,I),J=1,3),
     &                          (INC(J,I),J=1,2),CRATE(I),ANRATE(I),
     &                          CFACT(I),SEED(I),C1(I),T1(I),
     &                          L_TSS(I),RRG(I),RRX(I),IROCKS(I),
     &                          (CSC(J,I),J=1,4),
     &                          (STRAW(J,I),J=1,4)
C
C Set crop name and number of crops
C
      CNAME(I)=TEMP
      NUMCROP=I
      GOTO 101
C
C Format statements for parlis...
C
C Line 1: Crop Name; Line 2: Crop Number
9     FORMAT(A40/I3)
C Line 3: Cao (5 values)
10    FORMAT(F10.0,2X,F10.1,2X,F10.2,2X,F10.2,2X,F10.3,2X/
C Line 4: Ur (3 values) and Ut (3 values)
     &       F10.0,2X,F10.2,2X,F10.3,2X,F10.0,2X,F10.3,2X,F10.3,2X/
C Line 5: Inc (2 values), Crate, Anrate
     &       2(F10.2,2X),2(F10.3,2X)/
C Line 6: Cfact, Seed, C1, T1
     &       F10.2,2X,F10.1,2X,F10.3,2X,F10.1,2X/
C Line 7: L_TSs, RRG,RRX,IROCKS
     &       I10,2X,2(F10.1,2X),I10/
C Line 8: Csc (4 values)
     &       F10.0,2X,F10.2,2X,F10.3,2X,F10.3,2X/
C Line 9: Straw (4 values)
     &       F10.2,2X,F10.0,2X,F10.2,2X,F10.3,2X)
C
C Record error in the format of CROP_SUN.DAT
C
112   CONTINUE
      WRITE(*,*)'Warning! Error in general crop parameters!'
      WRITE(*,*)'Check format of crop parameter file, CROP_SUN.DAT'
      WRITE(15,*)'Warning! Error in general crop parameters!'
      WRITE(15,*)'Check format of crop parameter file, CROP_SUN.DAT'
C
C Try to open PARLIS.DAT
C
111   CONTINUE
      OPEN(44,FILE='PARLIS.DAT',STATUS='UNKNOWN',ERR=222)
C
C Read in parameters from PARLIS.DAT
C
      READ(44,9,ERR=222,END=222)TEMP,I
      READ(44,10,ERR=223,END=223)(CAO(J,I),J=1,5),
     &                          (UR(J,I),J=1,3),(UT(J,I),J=1,3),
     &                          (INC(J,I),J=1,2),CRATE(I),ANRATE(I),
     &                          CFACT(I),SEED(I),C1(I),T1(I),
     &                          L_TSS(I),RRG(I),RRX(I),IROCKS(I),
     &                          (CSC(J,I),J=1,4),
     &                          (STRAW(J,I),J=1,4)
C
C Set crop name
C
      CNAME(I)=TEMP
      GOTO 222
C
C Record error in the format of PARLIS.DAT
C
223   CONTINUE
      WRITE(*,*)'Warning! Error in current crop parameters!'
      WRITE(*,*)'Check format of crop parameter file, PARLIS.DAT'
      WRITE(15,*)'Warning! Error in current crop parameters!'
      WRITE(15,*)'Check format of crop parameter file, PARLIS.DAT'
C
C Write out parameters used to CROP_SUN.DAT
C
222   CONTINUE
      REWIND 43
      DO 100 I=1,NUMCROP
      WRITE(43,9)CNAME(I),I
      WRITE(43,10)(CAO(J,I),J=1,5),
     &            (UR(J,I),J=1,3),(UT(J,I),J=1,3),
     &            (INC(J,I),J=1,2),CRATE(I),ANRATE(I),
     &            CFACT(I),SEED(I),C1(I),T1(I),
     &            L_TSS(I),RRG(I),RRX(I),IROCKS(I),
     &            (CSC(J,I),J=1,4),
     &            (STRAW(J,I),J=1,4)
100   CONTINUE
C
C Close channels 43 and 44 
C
      CLOSE(43)
      CLOSE(44)
C
C Leave RESPAR
C
      END
C
C------------------------------------------------------------
C
      SUBROUTINE RCYCLE(MEND, IK, N_REMAIN, ICROP, 
     &                  WR, CACTOT, CATOT15, RNIN, RNIN15, CONVER_F)
C
C Subroutine to recycle C and N from crop debris back to RO pool
C

      IMPLICIT NONE

C INTEGER passed from subroutine
      INTEGER MEND, IK, N_REMAIN, ICROP

C REAL passed from subroutine
      REAL WR, CACTOT, CATOT15, RNIN, RNIN15, CONVER_F

C INTEGER local to subroutine
      INTEGER IROOT, NROOT

C REAL local to subroutine (WR15 not initialized before use)
      REAL W, W15, WR15, RNIN1, RN15
C
C Save variable values
C
      SAVE
c
C if RNIN (root N addition) exceeds CACTOT (crop n), allow roots to take
C 20% of crop N. Cumulate excess requirement (WR) until CACTOT large enough
C to subtract from. Subtract excess in up to 10 equal weekly amounts (W,W15).
C If occurs near harvest (usually only if CACTOT has been limited for some
C reason) check number of weeks remaining.
C
C if excess root N cumulated
C
         IF(WR.GT.0)THEN
c
C if week's root reqirement + cumulated excess < crop N
C
          IF(RNIN+WR.LT.CACTOT)THEN
           IROOT=1
C
C Check number of weeks remaining (NROOT)
C
           IF((MEND-IK)+1.LT.N_REMAIN)THEN
            NROOT=((MEND-IK)+1)
           ELSE
            NROOT=N_REMAIN
           END IF
C
C
C Calc weekly addition to root requirement
C
           W=WR/NROOT
           W15=WR15/NROOT
           WR=0
           WR15=0
          END IF
         END IF
C
C
c if adding excess
C
         IF(IROOT.GE.1)THEN
          IROOT=IROOT+1
          RNIN=RNIN+W
          RNIN15=RNIN15+W15
          IF(IROOT.EQ.NROOT+1)THEN
           IROOT=0
           W=0
           W15=0
          END IF
         END IF
C
         RNIN=RNIN+WR
         WR=0
C
         RNIN15=RNIN15+WR15
         WR15=0
C
C
C If week's root requirement greater than crop N, re-set root req. to 20%
C of crop total. Calculate excess (WR,WR15)
C

          IF(RNIN.GT.CACTOT)THEN
           RNIN1=RNIN
           RN15=RNIN15
           RNIN=0.2*CACTOT*CONVER_F
           IF(CACTOT.NE.0)THEN
             RNIN15=RNIN*CATOT15/CACTOT
C           ELSEIF(ICROP.GT.0)THEN
C             WRITE(15,*)'Error! Unexpected zero crop offtake in RCYCLE'
           ENDIF
           WR=RNIN1-RNIN
           WR15=RN15-RNIN15
          END IF
C
C Subtract root req. from crop total
C
         CACTOT=CACTOT-RNIN
         CATOT15=CATOT15-RNIN15
C
C leave RCYCLE
C
      RETURN
      END
C
C-----------------------------------------------------------------
C
      SUBROUTINE RESIDU(I_TS,ISTART,IEND,IYEAR,CINP,ANINP,IK,
     &                  CACTOT,CREQN,TRNIN,RNIN,RNIN15,CATOT15,
     &                  ICROP,SXORGN,SXORGN15,C_TS,SXORGC,TCINP,
     &                  XORGC,XORGN,ORGC,ORGN,CRATE,ANRATE,MCROP,
     &                  CONVER_F,CULTIVATE,NCULT,N_TS)
C
C Subroutine to add crop debris
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXCROP				! Max.no.of crop types allows
      PARAMETER (MAXCROP=36)
C
C Variables passed to/from calling subroutine
C
	INTEGER IK					! IN:No.timesteps from prev.harv. to current
	INTEGER IEND				! IN:Timesteps from sowing date to harvest
	INTEGER ISTART				! IN:Sowing month
      INTEGER I_TS	,N_TS			! Number of weeks since sowing, N_TS is day of the year
	INTEGER IYEAR				! IN:Current growing season number
      INTEGER ICROP				! IN:Current crop number
	INTEGER MCROP				! IN:Crop code number

	REAL ANINP					! IN/OUT:N added to soil as non-cartable deb. 
								!        in this timestep (kgN/ha/timestep)
	REAL ANRATE(0:MAXCROP)		! IN:Rate factor for C in non-cartable deb
	REAL CACTOT					! IN/OUT:N taken up by crop (kgN/ha)
	REAL CATOT15				! IN/OUT:N15 taken up by crop (kgN15/ha)
	REAL CINP					! IN/OUT:C added to soil as non-cartable deb.
								!        in this timestep (kgN/ha/timestep)
	REAL CONVER_F				! IN:Conversion between timestep & weeks
      INTEGER CULTIVATE			! IN(CALL): Code to cultivate this week (1=yes, 0=no)
	REAL CRATE(0:MAXCROP)		! IN:Rate factor for C in non-cartable deb
	REAL CREQN					! IN/OUT:Amount of N required by plant each 
	                            !        year (kgN/ha/year)
	REAL C_TS					! IN/OUT:Deb.C input in this timestep (kgC/ha)
	INTEGER NCULT				! IN(CALL): Number of cultivations
	REAL ORGC					! IN:Total org.C input (kgC/ha) 
	REAL ORGN					! IN:Total org.N input (kgN/ha) 
	REAL RNIN					! IN/OUT:Deb.N input in this timestep (kgN/ha)
	REAL RNIN15					! IN/OUT:Deb.N15 input in timestep (kgN15/ha)
	REAL SXORGC					! IN:Total org.N15 input (kgN/ha) 
	REAL SXORGN					! IN:Total org.N15 input (kgN/ha) 
	REAL SXORGN15				! IN/OUT:Total org.N15 input (kgN/ha) 
	REAL TCINP					! IN/OUT:Total org.C input (kgC/ha) 
	REAL TRNIN					! IN/OUT:Total litter N input (kgN/ha)
	REAL XORGC					! IN:Total org.C input (kgC/ha) 
	REAL XORGN					! IN:Total org.N input (kgN/ha) 
C
C If crop is growing set and add the carbon from crop debris
C
      IF(I_TS.GE.ISTART.AND.I_TS.LE.IEND)THEN
        CALL SETCARB(IYEAR,XORGC,XORGN,ORGC,ORGN,I_TS,ISTART,
     &                CRATE,MCROP,IEND,CONVER_F,ANRATE,CINP,ANINP,N_TS)
C
        IF(IYEAR.EQ.0.AND.IK.EQ.2)THEN
	    CACTOT=CREQN-TRNIN
        END IF
C
        CALL ADDCARB(CINP,RNIN,CACTOT,RNIN15,CATOT15,ANINP,ICROP,
     &                   C_TS)
C
C If in next cultivation month/week/day of growing season add the stubble N
C
      ELSEIF((NCULT.EQ.0.AND.IYEAR.GT.0.AND.IK.EQ.1).OR.
     &       (CULTIVATE.EQ.1))THEN
        RNIN=SXORGN
        IF(CACTOT.GT.0)THEN
c Stubble+chaff
           SXORGN15=SXORGN*0.27*CATOT15/CACTOT
C        ELSEIF(ICROP.GT.0)THEN
           SXORGN15=0
        ENDIF
        RNIN15=SXORGN15
    
        C_TS=SXORGC
C
C If crop is not growing add no residues
C
      ELSE
        C_TS=0
        RNIN=0
        RNIN15=0
      END IF
C
C Sum C additions
C
      TCINP=TCINP+C_TS
C
C Leave RESIDU
C
      RETURN
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE ROOTL(M,MCROP,IROCKS,ROOTS,RRX)
C
C Subroutine to calculate the layer that the roots have reached
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER,MAXDEPTH,MAXLAYER1
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
      INTEGER MAXCROP
	PARAMETER (MAXCROP=36)
C
C Variables passed to/from calling subroutine
C
      INTEGER M,MCROP,IROCKS(0:MAXCROP)
      REAL ROOTS,RRX(0:MAXCROP)
C
C Set layer using root exploration delay distance
C
	M = NINT((ROOTS)*((REAL(MAXLAYER1)/REAL(MAXDEPTH))))
	
C
C Do not allow roots to grow below maximum rooting depth
C
      IF(M.GT.IROCKS(MCROP)*(REAL(MAXLAYER1)/REAL(MAXDEPTH)))
     &  M=IROCKS(MCROP)*(REAL(MAXLAYER1)/REAL(MAXDEPTH))
C
C Leave ROOTL
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE ROOTN(ROOTS,YLD, UR, MCROP)
C
C Calculate root N requirement for various crops
C
C
C
      IMPLICIT NONE

      INTEGER MAXCROP

	PARAMETER (MAXCROP=36)

C INTEGER passed from subroutine
      INTEGER MCROP

C REAL passed from subroutine
      REAL ROOTS, YLD, UR(3,0:MAXCROP)
C
C Set root N uptake
C
      ROOTS=UR(1,MCROP)*(UR(3,MCROP)-EXP(UR(2,MCROP)*YLD))
C
C For zero yield set root N uptake to zero
C
      IF(NINT(YLD).EQ.0)ROOTS=0
C
C Leave ROOTN
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE SAVE_CROP(MCROP,CAO,UR,UT,INC,CRATE,ANRATE,CFACT,SEED,
     &                     C1,T1,RRG,RRX,CSC,STRAW,IROCKS,ROOTS,CREQN,
     &                     TOTC,ISTART,CINP,ANINP,TCNEED,BREQN,
     &                     CONVER_F,N_STEPS,N_REMAIN,
     &                     ISAVE)
C
C Subroutine to save and retrieve crop parameters 
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      INTEGER MAXCROP,MAXCROP1
      PARAMETER (MAXCROP=36)
	DATA MAXCROP1 /36/
      INTEGER IROCKS1(0:MAXCROP),MCROP1,ISTART1
      REAL CAO1(5,0:MAXCROP),UR1(3,0:MAXCROP),UT1(3,0:MAXCROP) 
	REAL INC1(2,0:MAXCROP),CRATE1(0:MAXCROP),ANRATE1(0:MAXCROP)
      REAL CFACT1(0:MAXCROP),SEED1(0:MAXCROP),C11(0:MAXCROP)
	REAL T11(0:MAXCROP),RRG1(0:MAXCROP),RRX1(0:MAXCROP)
      REAL CSC1(4,0:MAXCROP),STRAW1(4,0:MAXCROP)
	REAL ROOTS1,CREQN1,TOTC1,CINP1,ANINP1,TCNEED1,BREQN1
	REAL CONVER_F1
	INTEGER N_STEPS1,N_REMAIN1
	INTEGER I,J
C
C Variables passed to/from calling subroutine
C
      INTEGER IROCKS(0:MAXCROP),MCROP,ISAVE,ISTART
      REAL CAO(5,0:MAXCROP), UR(3,0:MAXCROP),UT(3,0:MAXCROP) 
	REAL INC(2,0:MAXCROP), CRATE(0:MAXCROP), ANRATE(0:MAXCROP)
      REAL CFACT(0:MAXCROP),SEED(0:MAXCROP),C1(0:MAXCROP),T1(0:MAXCROP)
      REAL RRG(0:MAXCROP),RRX(0:MAXCROP)
      REAL CSC(4,0:MAXCROP), STRAW(4,0:MAXCROP)
	REAL ROOTS,CREQN,TOTC,CINP,ANINP,TCNEED,BREQN
	REAL CONVER_F
	INTEGER N_STEPS,N_REMAIN
C
C Save variable
C
	SAVE
C
C Save parameters
C
	IF(ISAVE.EQ.1)THEN
	  N_REMAIN1=N_REMAIN
	  N_STEPS1=N_STEPS
	  CONVER_F1=CONVER_F
	  MCROP1=MCROP
	  ROOTS1=ROOTS
	  CREQN1=CREQN
	  TOTC1=TOTC
	  CINP1=CINP
	  ANINP1=ANINP
	  ISTART1=ISTART
	  TCNEED1=TCNEED
	  BREQN1=BREQN
	  DO 100 I=0,MAXCROP1
	    DO 200 J=1,5
            CAO1(J,I)=CAO(J,I)
200       CONTINUE
          DO 300 J=1,3
		  UR1(J,I)=UR(J,I)
		  UT1(J,I)=UT(J,I)
300       CONTINUE
          DO 400 J=1,2
		  INC1(J,I)=INC(J,I)
400       CONTINUE
          CRATE1(I)=CRATE(I)
		ANRATE1(I)=ANRATE(I)
		CFACT1(I)=CFACT(I)
		SEED1(I)=SEED(I)
		C11(I)=C1(I)
		T11(I)=T1(I)
		RRG1(I)=RRG(I)
		RRX1(I)=RRX(I)
	    IROCKS1(I)=IROCKS(I)
	    DO 500 J=1,4
		  CSC1(J,I)=CSC(J,I)
		  STRAW1(J,I)=STRAW(J,I)
500       CONTINUE
100     CONTINUE	  
C
C Retrieve parameterd
C
      ELSEIF(ISAVE.EQ.0)THEN
	  N_REMAIN=N_REMAIN1
	  N_STEPS=N_STEPS1
	  CONVER_F=CONVER_F1
	  MCROP=MCROP1
	  ROOTS=ROOTS1
	  CREQN=CREQN1
        TOTC=TOTC1
	  ISTART=ISTART1
	  CINP=CINP1
	  ANINP=ANINP1
        TCNEED=TCNEED1
	  BREQN=BREQN1
	  DO 600 I=0,MAXCROP1
	    DO 700 J=1,5
            CAO(J,I)=CAO1(J,I)
700       CONTINUE
          DO 800 J=1,3
		  UR(J,I)=UR1(J,I)
		  UT(J,I)=UT1(J,I)
800       CONTINUE
          DO 900 J=1,2
		  INC(J,I)=INC1(J,I)
900       CONTINUE
          CRATE(I)=CRATE1(I)
		ANRATE(I)=ANRATE1(I)
		CFACT(I)=CFACT1(I)
		SEED(I)=SEED1(I)
		C1(I)=C11(I)
		T1(I)=T11(I)
		RRG(I)=RRG1(I)
		RRX(I)=RRX1(I)
	    IROCKS(I)=IROCKS1(I)
	    DO 1000 J=1,4
		  CSC(J,I)=CSC1(J,I)
		  STRAW(J,I)=STRAW1(J,I)
1000      CONTINUE
600     CONTINUE	  
	ENDIF
C
C Leave SAVE_CROP
C
	END
C
C------------------------------------------------------------
C
      SUBROUTINE SENES(IYEAR, MCROP, JSTOP,OCROPN1,SRNIN,RORGN,TACTOT,
     &                 CONVER_F,CLOSSX,L_TSS, XORGN)
C
C Subroutine to set loss by senescence
C
      IMPLICIT NONE

      INTEGER MAXCROP

	PARAMETER (MAXCROP=36)

C INTEGER passed from subroutine
      INTEGER IYEAR, MCROP, JSTOP,L_TSS(0:MAXCROP)

C REAL passed from subroutine
      REAL OCROPN1,SRNIN,RORGN,TACTOT
	REAL CONVER_F,CLOSSX

	REAL XORGN

C REAL local to subroutine
      REAL XN
C
C Calc difference between annual root N requirement (RORGN) and
C root N already recovered(TRNIN).
C
C Note that : i) RORGN=SXORGN+ORGN i.e. N from previous crop's stubble
C                + current crop's total root N requirement.
C            ii) TRNIN=SXORGN+current root N uptake
C
C therefore RORGN-TRNIN=(SXORGN+ORGN)-(SXORGN+root N uptake) so that
C XN = net amount of N still to be taken up by roots.
C
C NJB
      CLOSSX=0
	IF(IYEAR.EQ.0)THEN
	XN=XORGN-SRNIN
	ELSE
	XN=RORGN-SRNIN
	END IF
C
C if total crop N is greater than final crop N requirement+amount
C still required by roots : lose excess amount by senescence
C
      IF(TACTOT.GT.OCROPN1+XN)THEN
       CLOSSX=(TACTOT-(OCROPN1+XN))/(NINT(L_TSS(MCROP)*(1/CONVER_F)))
      ELSE
       CLOSSX=0
      END IF
	
	JSTOP=1
C
C Leave SENES
C
      RETURN
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE SETC(INSTRAW, MCROP, TC, SC, YLD, CAO, CSC, STRAW)
C
C Subroutine to Calculate the total C returns for next crop
C
      IMPLICIT NONE

      INTEGER MAXCROP

	PARAMETER (MAXCROP=36)

C INTEGER passed from subroutine
      INTEGER INSTRAW, MCROP

C REAL passed from subroutine
      REAL TC, SC, YLD
	REAL CAO(5,0:MAXCROP),CSC(4,0:MAXCROP),STRAW(4,0:MAXCROP)

C REAL local to subroutine
      REAL RC, STRAWC
C
C Set total C from crop (=TC) and stubble & chaff (=SC)
C
       TC=CAO(2,MCROP)*(CAO(5,MCROP)-EXP(CAO(3,MCROP)*YLD))
       TC=(TC+CAO(4,MCROP))*CAO(1,MCROP)
       SC=CSC(4,MCROP)-(CSC(2,MCROP)*EXP(CSC(3,MCROP)*YLD))
       SC=SC*CSC(1,MCROP)
C
C Straw incorporation
C Set Straw yield: STRAW(1,MCROP)=HARVEST INDEX
C                  STRAW(2,MCROP)=FRACTION DM
C                  STRAW(3,MCROP)=FRACTION C
C                  STRAW(4,MCROP)=FRACTION N
C
        IF(INSTRAW.EQ.1)THEN
         RC=TC-SC
         IF(STRAW(1,MCROP).GT.0)THEN
          STRAWC=STRAW(2,MCROP)*YLD
          STRAWC=STRAWC*((1/STRAW(1,MCROP))-1)*STRAW(3,MCROP)
         ELSE
          STRAWC=0
         ENDIF
         TC=RC+STRAWC
         SC=STRAWC
        ENDIF

C
C Leave SETC
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE SETCARB(IYEAR,XORGC,XORGN,ORGC,ORGN,I_TS,ISTART,
     &                   CRATE,MCROP,IEND,CONVER_F,ANRATE,CINP,ANINP
     &					,n_ts)
C
C Subroutine to set the plant residues additions from the growing crop
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXCROP				! Max.no.of crop types allows
      PARAMETER (MAXCROP=36)
      CHARACTER*40 TEMP			! Temporary character string
	REAL CADD					! Total non-cartable C debris added 
								! through the growing season (kgC/ha/season)
	REAL CADD1					! C added up to end of prev.timestep (kgC/ha)
	REAL CADD2					! C added up to end of this timestep (kgC/ha)
	REAL ANADD					! Total non-cartable N debris added 
								! through the growing season (kgN/ha/season)
	REAL ANADD1					! N added up to end of prev.timestep (kgN/ha)
	REAL ANADD2					! N added up to end of this timestep (kgN/ha)
	INTEGER I					! Local counter
	INTEGER J					! Local counter
C
C Variables passed to/from calling subroutine
C
	INTEGER IEND				! IN:Timesteps from sowing date to harvest
	INTEGER ISTART				! IN:Sowing month
      INTEGER I_TS,N_TS			! Number of weeks since sowing
	INTEGER IYEAR				! IN:Current growing season number
      INTEGER MCROP				! IN:Crop code number

	REAL ANINP					! OUT:N added to soil as non-cartable deb. 
								!        in this timestep (kgN/ha/timestep)
	REAL ANRATE(0:MAXCROP)		! IN:Rate factor for C in non-cartable deb
	REAL CINP					! OUT:C added to soil as non-cartable deb.
								!        in this timestep (kgN/ha/timestep)
	REAL CONVER_F				! IN:Conversion between timestep & weeks
	REAL CRATE(0:MAXCROP)		! IN:Rate factor for C in non-cartable deb
	REAL ORGC					! IN:Total org.C input (kgC/ha) 
	REAL ORGN					! IN:Total org.N input (kgN/ha) 
	REAL XORGC					! IN:Total org.C input (kgC/ha) 
	REAL XORGN					! IN:Total org.N input (kgN/ha) 
C
C Save Variable values
C
      SAVE
C
C Set total plant debris
C
      IF(IYEAR.EQ.0)THEN
        CADD=XORGC
        ANADD=XORGN
      ELSE
        CADD=ORGC
        ANADD=ORGN
      END IF
C
C Calculate additions from up to Previous month/week/day
C
      IF(I_TS.EQ.ISTART+1)THEN
        IF(IYEAR.EQ.0)THEN
          CADD2=CADD*EXP(-CRATE(MCROP)*(IEND-(I_TS-1))*CONVER_F)
          ANADD2=ANADD*EXP(-ANRATE(MCROP)*(IEND-(I_TS-1))*CONVER_F)
        ELSE
          CADD2=0
          ANADD2=0
        END IF
      END IF
C
C Calculate additions up to current month/week/day
C
      CADD1=CADD*EXP(-CRATE(MCROP)*(IEND-I_TS)*CONVER_F)
	IF(I_TS.EQ.1)then
	  CINP= 0.0
	ELSE
        CINP=CADD1-CADD2
	ENDIF
      CADD2=CADD1
      
	ANADD1=ANADD*EXP(-ANRATE(MCROP)*(IEND-I_TS)*CONVER_F)
	IF(I_TS.EQ.1)then
	  ANINP=0.0
	ELSE
        ANINP=ANADD1-ANADD2
	ENDIF
      ANADD2=ANADD1
C
C Leave SETCARB
C
        RETURN
        END
C
C---------------------------------------------------------------
C
      SUBROUTINE SETCROPN(IYEAR, INSTRAW, MCROP, LCROP, ICROP,
     &                    SXORGN, SORGN, YLD, PREYLD, CREQN, 
     &                    OCROPN, ROOTS, EXYLD, TOTC,
     &                    CFACT, STRAW, INC, SEED, UR, UT, XORGN, ORGN)
C
C Subroutine to set the crop N requirement
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXCROP				! Max.no.of crop types allows
	PARAMETER (MAXCROP=36)
	REAL CLOSS
      INTEGER IFIT
      REAL STRAWN
C
C Variables passed to/from calling subroutine
C
	REAL CFACT(0:MAXCROP)		! IN:Fraction of N lost by senescence 
	REAL CREQN					! IN:Amount of N required by plant each yr
	                            !        (kgN/ha/year)
	REAL EXYLD					! IN:Yield of current crop (t/ha)
	INTEGER ICROP				! IN:Current crop code
      REAL INC(2,0:MAXCROP)
	INTEGER INSTRAW				! IN:Straw incorporated? 0=No 1=Yes
	INTEGER IYEAR				! IN:Current growing season number
	INTEGER LCROP				! IN:Previous crop code
	INTEGER MCROP				! IN:Crop type code number
	REAL OCROPN					! IN:N uptake of crop (kgN/ha) 
								! (0->calculate within model)
	REAL ORGN					! IN:Total org.N input (kgN/ha) 
	REAL PREYLD					! IN:Yield of previous crop (t/ha)
	REAL ROOTS					! OUT:N in below ground plant (kgN/ha)
	REAL SEED(0:MAXCROP)
	REAL SORGN					! IN:Total org.N15 input (kgN/ha) 
	REAL STRAW(4,0:MAXCROP)  
	REAL SXORGN					! IN:Total org.N15 input (kgN/ha) 
	REAL TOTC
	REAL UR(3,0:MAXCROP)
	REAL UT(3,0:MAXCROP)
	REAL XORGN					! IN:Total org.N input (kgN/ha) 
	REAL YLD					! IN:Yield (t/ha) 
C
C In first year or if crop N requirement has not been set
C code IFIT for calculation of crop N uptake
C
      IF(OCROPN.LT.0.01.OR.IYEAR.EQ.0)THEN
       IFIT=0
      ELSE
       IFIT=1
      END IF
C
C Set crop index number for current crop (ICROP) or initial crop
C (LCROP) if running from anthesis to harvest of first year.
C
      IF(IYEAR.EQ.0)THEN
       MCROP=LCROP
       YLD=PREYLD
      ELSE
       MCROP=ICROP
       YLD=EXYLD
      END IF
C
C Set crop and root N
C
      CALL ROOTN(ROOTS,YLD, UR, MCROP)
C
C If crop N uptake has been set... 
C
      IF(IFIT.EQ.1)THEN
C
C ...proportion root N requirement according to N in tops
C
        CALL CROPN(CREQN,YLD, UT, MCROP)
	  ROOTS=ROOTS*OCROPN/CREQN
C
C ...enter measured crop N
C
        CREQN=OCROPN
        IF(CREQN.LE.0)CREQN=0.0001
C
C If crop N has not been set estimate crop N in relation to yield
C
      ELSE
         CALL CROPN(CREQN,YLD, UT, MCROP)
         OCROPN=CREQN
         IF(CREQN.LE.0)CREQN=0.0001
      END IF
C
C Allow fraction of crop N which may be lost by senescence after anthesis
C (e.g.5%)
C
      CLOSS=CFACT(MCROP)*CREQN
C
C In first year always incorporate chaff
C
      IF(IYEAR.EQ.0)THEN
        STRAWN=0
        IF(INSTRAW.EQ.1)THEN
          IF(STRAW(1,MCROP).GT.0)THEN
            STRAWN=STRAW(2,MCROP)*PREYLD
            STRAWN=STRAWN*((1/STRAW(1,MCROP))-1)*STRAW(4,MCROP)
          ENDIF
        ENDIF
      SXORGN=INC(1,MCROP)*CREQN+STRAWN
      XORGN=INC(2,MCROP)*ROOTS
      CREQN=CREQN+CLOSS-SEED(MCROP)+XORGN
      IF(CREQN.LE.0)CREQN=0.0001
C
C In subsequent years...
C
      ELSE
        STRAWN=0
        IF(INSTRAW.EQ.1)THEN
          IF(STRAW(1,MCROP).GT.0)THEN
            STRAWN=STRAW(2,MCROP)*EXYLD
            STRAWN=STRAWN*((1/STRAW(1,MCROP))-1)*STRAW(4,MCROP)
          ENDIF
        ENDIF
       SORGN=INC(1,MCROP)*CREQN+STRAWN
C
C Set ORGN = N requirement of the roots
C
       ORGN=INC(2,MCROP)*ROOTS
C
C Set CREQN = N requirement f the tops
C
       CREQN=CREQN+CLOSS-SEED(MCROP)+ORGN
       IF(CREQN.LE.0)CREQN=0.0001
      END IF
C
      TOTC=CREQN
C
C Leave SETCROPN
C
      RETURN
      END
C
C------------------------------------------------------------------
C
      SUBROUTINE SETEND(IYEAR,IEND,ISTHARV,N_STEPS,ISOWN,MEND,NSOW)
C
C Subroutine to set the end of the growing season
C
      IMPLICIT NONE

C INTEGER passed to subrontine
      INTEGER IYEAR, IEND, ISTHARV, N_STEPS, ISOWN, MEND, NSOW
C
C In first year...
C
       IF(IYEAR.EQ.0)THEN
        IEND=(ISTHARV+N_STEPS)-ISOWN
C
C In subsequent years...
       ELSE
        IEND=(MEND-NSOW)+1
       END IF
C
C Leave SETEND
C
       RETURN
       END
C
C----------------------------------------------------------------------------
C
	SUBROUTINE SETEXTRACROPPARAM(YLD,PLANTBIOMASS,PLANTHEIGHT,LAI,
     &							canMax)
C
C This subroutine estimates biomass [kg ha-1], LAI and plant height [cm] from yield
C according to fitted functions for maize by Pia Gottschalk
C Fitted equations:	a) plant height = 33.655 * ln(YLD)+154.08 ; R2=0.88
C					b) LAI			= 0.9049 * ln(YLD) + 1.4738 ; R2=0.89
C					c) biomass		= 4431.1 * EXP(0.1948*YLD) ; 0.99
C
	IMPLICIT NONE
C
C Variables passed to calling subroutine
C
	REAL YLD
C
C Variables passed from calling subroutine
C
	REAL PLANTBIOMASS		! Plant biomass at harvest [ kg ha-1]
	REAL PLANTHEIGHT		! Plant height at harvest [cm]
	REAL LAI				! Plant LAI at harvest
	REAL canMAX				! maximum amount of water that can be trapped in the canopy when the canopy is fully developed (mm H2O)
C
	PLANTBIOMASS = 4431.1 * EXP(0.1948*YLD)
	PLANTHEIGHT = 33.655 * LOG(YLD)+154.08
	LAI = 0.9049 * LOG(YLD) + 1.4738
C
C Setting canMax
	canMax=2.8
C Leave SETEXTRACROPPARAM
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE SETTIME_SUNDIAL_CROP(IYEAR,NXYEARS,NDATE,MEND,FIXEND,
     &                    ISTHARV,LHARV,LCROP,IANTHES,
     &                    SECONDS,CONVER_F,N_STEPS,N_REMAIN,L_TSS,
     &                    NSOW,ISOWN,IHARV,
     &                    NORGM,IORGM,NFERT,IFERT)
C
C Subroutine to set the timing of field operations with respect to previous harvest
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXFERT,MAXORGM,MAXCROP
	PARAMETER(MAXFERT=5,MAXORGM=52,MAXCROP=36)
	INTEGER NF
C
C Variables passed to/from calling subroutine
C
      INTEGER IYEAR,NXYEARS,NDATE,MEND,FIXEND,ISTHARV,LHARV,LCROP
	INTEGER IANTHES,L_TSS(0:MAXCROP),NSOW,ISOWN,IHARV
	INTEGER NORGM,IORGM(MAXORGM),NFERT,IFERT(MAXFERT),N_STEPS,N_REMAIN
	REAL CONVER_F
	REAL SECONDS
C
C Set time factors from SECONDS
C
      N_STEPS=(365.25*24*60*60)/SECONDS
	CONVER_F=SECONDS/(7*24*60*60)
	N_REMAIN=10/CONVER_F
C
C Set timing of events wrt previous harvest
C ...in previous year, set sowing to anthesis 
C    and length of season to time from anthesis to harvest

	IF(IYEAR.EQ.0)THEN
        IANTHES=ISTHARV - NINT(L_TSS(LCROP) * (1 / CONVER_F))
        ISOWN=IANTHES
        NDATE=2
        MEND=ISTHARV-IANTHES
        MEND=MEND+1
C
C ...after last growing season set length of season to time from harvest 
C    to the fixed end of simulation
C
      ELSEIF(IYEAR.EQ.NXYEARS+1)THEN
        MEND=FIXEND-LHARV-1
C
C ...in other growing seasons, set season length to time from previous harvest 
C    to harvest and calculate timing of operations
C
      ELSE
        NDATE=2
        MEND=IHARV-LHARV
C
C Calc. no. of weeks between sowing and harvest (NSOW)
C
        NSOW=ISOWN-LHARV	
C
C Calc no. of weeks between harvest and fertilizer appln.
C
        DO 21 NF=1,NFERT
          IFERT(NF)=IFERT(NF)-LHARV
21      CONTINUE
C
C Calc no. of weeks between harvest and organic manure appln.
C
        DO 121 NF=1,NORGM
          IORGM(NF)=IORGM(NF)-LHARV
121     CONTINUE
C
       END IF
C
C Leave SETTIME_CROP
C
	END
C
C-----------------------------------------------------------
C
      SUBROUTINE SETVAR3(IS_TS,CACTOT,CATOT15, CTOT, TCNEED, JSTOP)
C
C Subroutine to initalize variables
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER IS_TS				! OUT:Crop counter
	INTEGER JSTOP

C
C Variables passed to/from calling subroutine
C
	REAL CTOT					! IN:Total litter C input (kgC/ha)
	REAL TCNEED, CACTOT, CATOT15
C
C Initialize variables
C
      CTOT=0
      TCNEED=0
      IS_TS=0
      JSTOP=0
      CACTOT=0
      CATOT15=0
C
C Leave SETVAR3
C
      RETURN
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE SETVARS_CROP()
C
C Subroutine to initialise variables for crop calculations
C

C
C Initialise variables
C
      END
C
C------------------------------------------------------------
C
      SUBROUTINE SOWING(IYEAR,INSTRAW,MCROP,LCROP,ICROP,IEND,ISTHARV,
     &                  N_STEPS,ISOWN,MEND,NSOW,IS_TS,JSTOP,
     &                  OCROPN,SXORGN,SORGN,YLD,CREQN,PREYLD,
     &                  TORGC,CTOT,EXYLD,TOTC,SORGC,RORGN,ORGC, 
     &                  CACTOT,CATOT15,TCNEED,BREQN,UR,
     &                  UT,SEED,INC,CFACT,XORGN,ORGN,ROOTS,STRAW)


      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXCROP				! Max.no.of crop types allows
	PARAMETER (MAXCROP=36)
C
C Variables passed to/from calling subroutine
C
	REAL BREQN
	REAL CACTOT					! OUT:N taken up by crop (kgN/ha)
	REAL CATOT15				! OUT:N15 taken up by crop (kgN15/ha)
	REAL CFACT(0:MAXCROP)		! IN:Fraction of N lost by senescence 
	REAL CREQN					! IN:Amount of N required by plant each yr
	                            !        (kgN/ha/year)
	REAL CTOT					! OUT:Total litter C input (kgC/ha)
	REAL EXYLD					! IN:Yield of current crop (t/ha)
	INTEGER ICROP				! IN:Current crop code
	INTEGER IEND				! IN:Timesteps from sowing date to harvest
      REAL INC(2,0:MAXCROP)
	INTEGER INSTRAW				! IN:Straw incorporated? 0=No 1=Yes
	INTEGER IS_TS				! IN:Crop counter
	INTEGER ISOWN				! IN:Timesteps from 01/01/01 to sowing date 
	INTEGER	ISTHARV				! IN:Timesteps from 01/01/01 to first harvest
	INTEGER IYEAR				! IN:Current growing season number
	INTEGER JSTOP				! IN:Code to stop simulation 0=No 1=Yes
	INTEGER LCROP				! IN:Previous crop type
	INTEGER MCROP				! IN:Crop type code number
	INTEGER MEND				! IN:Number of timesteps from prev.harv.-harv.
	INTEGER N_STEPS				! IN:No.timesteps in a year
	INTEGER NSOW				! IN:Timesteps from prev.harv. to sowing date
	REAL OCROPN					! IN:N uptake of crop (kgN/ha) 
								! (0->calculate within model)
	REAL ORGC					! IN:Total org.C input (kgC/ha) 
	REAL ORGN					! IN:Total org.N input (kgN/ha) 
	REAL PREYLD					! IN:Yield of previous crop (t/ha)
	REAL ROOTS
	REAL RORGN					! IN:Total org.N input (kgN/ha) DONT PASS?
	REAL SEED(0:MAXCROP)
	REAL SORGC					! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SORGN					! IN:Total org.N15 input (kgN/ha) 
	REAL STRAW(4,0:MAXCROP)  
	REAL SXORGN					! IN:Total org.N15 input (kgN/ha) 
	REAL TCNEED					! OUT:Total N requirement (kgN/ha)
      REAL TORGC					! IN:Total org.C input (kgC/ha) 
	REAL TOTC
	REAL UR(3,0:MAXCROP)
	REAL UT(3,0:MAXCROP)
	REAL YLD					! IN:Yield (t/ha) 
	REAL XORGN					! IN:Total org.N input (kgN/ha)
C
C Set N uptake for the sown crop
C
      CALL SETCROPN(IYEAR, INSTRAW, MCROP, LCROP, ICROP,
     &                    SXORGN, SORGN, YLD, PREYLD, CREQN, 
     &                    OCROPN, ROOTS, EXYLD, TOTC,
     &                    CFACT, STRAW, INC, SEED, UR, UT, XORGN, ORGN)
C
C Set Organic C and N returned as debris from the previous crop
C
      ORGC=TORGC-SORGC
      RORGN=SXORGN+ORGN
C
C Reset crop code number, MCROP
C
      MCROP=ICROP
C
C Set the end of the growing season for the sown crop
C
      CALL SETEND(IYEAR, IEND, ISTHARV, N_STEPS, ISOWN, MEND, NSOW)
C
C Initialize crop variables
C
      CALL SETVAR3(IS_TS,CACTOT,CATOT15,CTOT,TCNEED,JSTOP)
      BREQN=TOTC
C
C Leave SOWING
C
      RETURN
      END
C
C----------------------------------------------------------
C
C
      SUBROUTINE UPTAKE(DDAYS, CTOT, T1, C1, CREQN, MCROP, CUPTN)
C
C Subroutine to calculate the crop N uptake
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      INTEGER MAXCROP
	PARAMETER (MAXCROP=36)
	REAL TOT1
C
C Variables passed to/from calling subroutine
C
      INTEGER MCROP
      REAL DDAYS, CTOT, T1(0:MAXCROP), C1(0:MAXCROP), CREQN, CUPTN
C
C Set TOT1 to uptake up to last week
C
      TOT1=CTOT
C
C Calculate CTOT, uptake up to this week
C
      IF(T1(MCROP).GT.0)THEN
        CTOT=(CREQN**(-1/T1(MCROP))+EXP(-C1(MCROP)*DDAYS))**(-T1(MCROP))
      ELSE
        WRITE(*,*)'Error! Unexpected zero in crop parameters'
        WRITE(15,*)'Error! Unexpected zero in crop parameters'
      ENDIF
C
C Caculate the amount of n that was taken up this week
C
      CUPTN=CTOT-TOT1
	IF(CUPTN.LT.0)CUPTN=0
C
C Leave UPTAKE
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE ZEROUP(CUPTN,CACT,CACT15)
C
C Subroutine to zero crop uptake parameters
C
      IMPLICIT NONE
C
C Variables passed to/from calling subroutine
C
	REAL CACT				! Actual N uptake (kgN/ha)
	REAL CACT15				! Actual N15 uptake (kgN15/ha)
      REAL CUPTN				! Crop N requirement (kgN/ha)
C
      CUPTN=0
      CACT=0
      CACT15=0
C
C Leave ZEROUP
C
      RETURN
      END