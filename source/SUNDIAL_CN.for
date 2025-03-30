C-----------------------------------------------------------
C This is new file
C 
C Rothamsted Carbon and Nitrogen Turnover Model
C Soil C and N Routines
C
C Adapted from SUNDIAL (MAGEC)
C by Jo Smith (Univ.Aberdeen) & Kevin Coleman (Rothamsted Research)
C 02/03/05
C
C Modified for highly organic soils
C by Jo Smith, Bente Foereid & Matt Aitkenhead (Univ.Aberdeen)
C Started 01/04/05
C
C Modified for tropical soils
C by Jo Smith & Pia Gottschalk (Univ.Aberdeen)
C Started 01/04/05
C
C Modified for spatial applications
C by Jo Smith (Univ.Aberdeen)
C Started 01/08/06
C 
C-------------------------------------------------------------
C EXTERNAL ROUTINES FOR SITE SPECIFIC RUNS
C 1. ADD_SUNDIAL_SOILCN
C 2. DUMP_SOILCN
C 3. GET_CTON_FALLOON
C 4. GET_CTON_FOEREID
C 5. GET_CTON_BRADBURY
C 6. GET_SUNDIAL_SOILPARS
C 7. INIT_SUNDIAL_SOILCN
C 8. INIT_SUNDIAL_SOILCN_FIXED
C 9. INIT_SUNDIAL_SOILCN_NOPARS
C 10. MICROBIAL_SUNDIAL_SOILCN
C 11. PHYSICAL_SUNDIAL_SOILCN
C 12. RESTART_SUNDIAL_SOILCN
C 13.RUN_ROTHC_AND_N_TO_EQUILIBRIUM
C
C-------------------------------------------------------------
C
      SUBROUTINE ADD_PLANT_DIST(C_TS,RNIN,RNIN15,PI_CEQ_MON,THISMON,
     &                          PI_C,PI_N,PI_N15) 
C
C     to add C, N and N15 input as distributed in PI_CEQ_MON
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/

      INTEGER IL					! Layer counter
	REAL PI_CEQ					! Total equilibrium plant C input this month (kgC/ha)
C
C Variables passed to/from calling subroutine
C
	REAL C_TS					! IN: Litter C input in this timestep (kgC/ha)
      REAL PI_C(MAXLAYER)			! OUT: Plant input C to soil (kgC/ha/step)
	REAL PI_CEQ_MON(12,MAXLAYER)! IN: Equilibrium plant C input each month
								!         in each layer (kgC/ha/month/layer)
      REAL PI_N(MAXLAYER)			! OUT: Plant input N to soil (kgC/ha/step)
      REAL PI_N15(MAXLAYER)		! OUT: Plant input N15 to soil(kgC/ha/step)
								!(1=arable,2=grassland,3=forestry,4=semi-natural)
	REAL RNIN					! IN:Litter N input in this timestep (kgN/ha)
	REAL RNIN15					! IN:Litter N15 input in timestep (kgN15/ha)
	INTEGER THISMON				! IN: This month
C
C Cumulate plant input distribution in this month
C
      PI_CEQ=0
      DO 100 IL=1,MAXLAYER1
        PI_CEQ=PI_CEQ+PI_CEQ_MON(THISMON,IL)
100   CONTINUE
C
C Distribute the calculated plant input this month
C
      DO 200 IL=1,MAXLAYER1
        IF(PI_CEQ.GT.0)THEN
	    PI_C(IL)=C_TS*PI_CEQ_MON(THISMON,IL)/PI_CEQ
	    PI_N(IL)=RNIN*PI_CEQ_MON(THISMON,IL)/PI_CEQ
          PI_N15(IL)=RNIN15*PI_CEQ_MON(THISMON,IL)/PI_CEQ
	  ELSE
	    PI_C(IL)=0
	    PI_N(IL)=0
          PI_N15(IL)=0
	  ENDIF
200    CONTINUE
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE ADD_SUNDIAL_SOILCN(
C INPUTS: time, environment, organic manure, atmospheric, 
     &                              SECONDS,
     &                              NSOIL,RAIN,
     &                              ORGMANF,JORGMNF,IOLABNF,      
     &							  FERTNF,TFERTNF,ILABNF,IVOLNF,
     &                              ATM,								  
C OUTPUTS: organic manure, atmospheric
     &                              FYMFERT,FYMFERT15,TFYMC,			  
     &                              VOLAT,VOLAT15,					  
C INPUTS/OUTPUTS: Soil C&N
     &                              SOILN,SOIL15,AMMN,AMMN15,			  
     &                              DPMCARB0,DPMNIT0,DPMNLAB0,		  
     &                              RPMCARB0,RPMNIT0,RPMNLAB0,		  
     &                              HCARB0,HNIT0,HNLAB0,TIN,TAM)		  
C
C Subroutine to add stuff to SUNDIAL soil C and N
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXORGM				! Max.no.of org.manure applications allowed
	INTEGER MAXSOIL				! Max.no.of soil types
	INTEGER MAXFERT				! Max.no.of fertiliser applications allowed 
	PARAMETER(MAXSOIL=50,MAXORGM=52,MAXFERT=5)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXFERTSTEPS		! Maximum number of timesteps 
								!	over which fertiliser is added
      PARAMETER (MAXFERTSTEPS=30*60*24*7*3)	! Timesteps = 3 wks (max.30 secs)
      INTEGER ISAVE				! Code to save or retrieve variables
      INTEGER NUMSOIL				! Number of soils defined
	REAL BYRAIN					! IN/OUT:Rain lost by bypass flow mm/timestep	
	REAL DFACT(MAXLAYER)		! IN/OUT: Denitrification factor
	REAL ANIT					! Ammonium immobilised (kgN/ha)
	REAL ANIT15					! Ammonium N15 nitrified (kgN/ha)
	REAL FANIT					! Fertiliser N nitrified (kgN/ha)
	REAL FANIT15				! Fertiliser N15 nitrified (kgN/ha)
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: % clay content in this layer
      INTEGER LUARRAY(MAXSOIL)	! IN/OUT:Land use before equilibrium 
C
C Variables passed to/from calling subroutine
C ...Timing factors
C
	REAL SECONDS				! IN:Number of seconds in one timestep
      REAL CONVER_F				! IN/OUT:Conversion to correct timestep
C
C ...Weather factors
C
	REAL RAIN					! IN: Rainfall (mm/timestep)
C
C ...Manure factors
C
      REAL FPROPX(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
	REAL FYMCX(MAXORGM)			! IN/OUT:Proportion of C in FYM
	REAL FYMNX(MAXORGM)			! IN/OUT:Proportion of Organic N in FYM
	REAL FYMAX(MAXORGM)			! IN/OUT:Proportion of ammonium-N in FYM
	REAL FYMWX(MAXORGM)			! IN/OUT:Amount of water added in FYM
	REAL FYMLOSS(MAXORGM)		! IN/OUT:Prop.FYM lost by volatn/timestep
	REAL FYMPOS(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
      REAL FYMSTART(MAXORGM)		! IN/OUT:Amount of rainfall in 1 week, 
								!		 below which volatilisation will occur
	INTEGER JORGMNF				! IN:Type of manure application
	INTEGER IOLABNF				! IN:Labelling on manure? 0=No 1=Yes
	REAL ORGMANF				! IN:Amount of manure applied (t/ha)
C
C ...Fertiliser factors
C
	REAL FERTNF					! INOUT:Amount of fertiliser applied (kgN/ha)
	REAL TFERTNF(3)				! INOUT:Prop.NO3,NH4,urea in fertiliser
	INTEGER ILABNF				! INOUT:Labelling on fertiliser? 0=No 1=Yes
	INTEGER IVOLNF				! INOUT:Ammonium salt (0=sulphate, 1=other)
	REAL FERTADD(MAXFERTSTEPS,4)	! INOUT:Fertiliser N added in this time 
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
	REAL FERTADD15(MAXFERTSTEPS,4)	! INOUT:Fertiliser N15 added in this time
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
C
C ...Soil factors
C
	CHARACTER*40 SNAME(MAXSOIL)	! IN/OUT:Soil name
	INTEGER NSOIL				! IN:Soil code number
	INTEGER NFERT				! IN:No.fertiliser applications to this crop
	REAL ATM					! IN:Atmospheric N input (kgN/ha/timestep)
	REAL VOLAT					! IN:N lost by volatilisation (kgN/ha)
	REAL VOLAT15				! IN:N15 lost by volatilisation (kgN15/ha)
	REAL FYMFERT				! IN:N input by fertiliser & FYM (kgN/ha)
	REAL FYMFERT15				! IN:N15 input by fertiliser & FYM (kgN15/ha)
	REAL TFYMC					! Total C added in FYM (kgC/ha)
  	REAL WLEACH					! IN:N leached (kgN/ha)
	REAL WLEACH15				! IN:N15 leached (kgN15/ha)
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN:N in soil humus (kgN/ha/layer)
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN:N15 in soil biomass (kgN15/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
 	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
	REAL TAM(MAXLAYER)			! IN:15N/N in the ammonium in the layer
	REAL TIN(MAXLAYER)			! IN:15N/N in the nitrate in the layer
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: pH above which decomp.max.
	REAL HY(MAXSOIL,MAXLAYER)	! IN/OUT:Stable N:C ratio of BIO & HUM pools
	REAL BPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL HPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL BPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from biomass decompositn
	REAL HPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from humus decompositn
	REAL ALPHA(MAXLAYER)			! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL BETA(MAXLAYER)				! IN/OUT: Prop.HUM produced on BIO decompn
	REAL GAMMA(MAXLAYER)			! IN/OUT: Prop.BIO produced on HUM decompn 
	REAL DELTA(MAXLAYER)			! IN/OUT: Prop.HUM produced on HUM decompn 
	REAL HZ1(MAXLAYER)				! IN/OUT: N:C ratio for steady state
	REAL BIOP(MAXSOIL,MAXLAYER)		! IN/OUT: BIO/TOC 
	REAL DPMRATE(MAXLAYER)			! IN/OUT: Rate constant for DPM decompn
	REAL RPMRATE(MAXLAYER)			! IN/OUT: Rate constant for RPM decompn
	REAL BRATE(MAXLAYER)			! IN/OUT: Rate constant for HUM decompn
	REAL HRATE(MAXLAYER)			! IN/OUT: Rate constant for BIO decompn
	REAL CRIT(MAXSOIL,MAXLAYER)		! IN/OUT: Minimum nitrate level (Sozanska)
	REAL BIORATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for BIO decompn/yr
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for HUM decompn/yr
	REAL IOMARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium IOM (kgC/ha/layer)
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium TOC (kgC/ha/layer)
  	REAL ICFACTOR(MAXLAYER)			! IN:Adjustment in rate needed to	
									!    achieve measured NPP and TOC	
C
C Retrieve soil CN characteristics 
C
      ISAVE=0
      CALL SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,CONVER_F,
     &             SECONDS,NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)
	CALL SAVE_FERTADD(SECONDS,FERTADD,FERTADD15,ISAVE)
C
C Add fertiliser
C
      CALL ADDFERT(FERTNF,TFERTNF,ILABNF,IVOLNF,FERTADD,FERTADD15,
     &             SECONDS)
C
C  Add any FYM
C
      CALL ADDFYM(FPROPX,JORGMNF,FYMCX,FYMNX,FYMAX,FYMWX,
     &                  RAIN,FYMSTART,FYMLOSS,ORGMANF,HY,
     &                  FYMPOS,HCARB0,HNIT0,DPMCARB0,RPMCARB0,
     &                  DPMNIT0,RPMNIT0,HNLAB0,DPMNLAB0,RPMNLAB0,
     &                  AMMN,VOLAT,IOLABNF,AMMN15,VOLAT15,NSOIL,
     &                  FYMFERT,FYMFERT15,CONVER_F)
C
C Add any atmospheric dry deposition
C
      CALL ADDATM(ATM,SOILN)
C
C
C calculate the proportion of N15 added
C
      CALL PROPN15(SOILN,SOIL15,TIN,AMMN,AMMN15,TAM)
C
C Save soil CN characteristics 
C
      ISAVE=1
      CALL SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,CONVER_F,
     &             SECONDS,NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)
	CALL SAVE_FERTADD(SECONDS,FERTADD,FERTADD15,ISAVE)
      END
C
C--------------------------------------------------------

	 SUBROUTINE DUMP_SOILCN(IDUMP)

C Subroutine to write initial C & N parameters to file (Initial_CN_pars.DAT)
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/
	INTEGER MAXSOIL				! Max.no.of soil types
	INTEGER MAXORGM				! Max.no.of org.manure applications allowed
      PARAMETER (MAXSOIL=50,MAXORGM=52)
      INTEGER IL					! Layer counter
      INTEGER ISAVE				! Code to save or retrieve variables
	INTEGER I					! Local counter
	INTEGER J					! Local counter
C
C Variables passed to/from calling subroutine
C
      INTEGER IDUMP				! IN:Dump data to file or retrieve 
								! from file. (1=dump;0=retrieve)
C
C ...Timing factors
C
	REAL SECONDS				! IN:Number of seconds in one timestep
	REAL CONVER_F
C
C ...Weather factors
C
	REAL BYRAIN					! IN/OUT:Rain lost by bypass flow mm/timestep	
C
C ...Soil factors
C
      INTEGER NUMSOIL					! Number of soils defined
	REAL HY(MAXSOIL,MAXLAYER)		! IN/OUT:Stable N:C ratio of BIO&HUM pools
	REAL BPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL HPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL BPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from biomass decompositn
	REAL HPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from humus decompositn
	REAL ALPHA(MAXLAYER)			! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL BETA(MAXLAYER)				! IN/OUT: Prop.HUM produced on BIO decompn
	REAL GAMMA(MAXLAYER)			! IN/OUT: Prop.BIO produced on HUM decompn 
	REAL DELTA(MAXLAYER)			! IN/OUT: Prop.HUM produced on HUM decompn 
	REAL HZ1(MAXLAYER)				! IN/OUT: N:C ratio for steady state
	REAL BIOP(MAXSOIL,MAXLAYER)		! IN/OUT: BIO/TOC 
	REAL DPMRATE(MAXLAYER)			! IN/OUT: Rate constant for DPM decompn
	REAL RPMRATE(MAXLAYER)			! IN/OUT: Rate constant for RPM decompn
	REAL BRATE(MAXLAYER)			! IN/OUT: Rate constant for HUM decompn
	REAL HRATE(MAXLAYER)			! IN/OUT: Rate constant for BIO decompn
	REAL CRIT(MAXSOIL,MAXLAYER)		! IN/OUT: Minimum nitrate level (Sozanska)
	REAL BIORATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for BIO decompn/yr
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for HUM decompn/yr
	REAL IOMARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium IOM (kgC/ha/layer)
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium TOC (kgC/ha/layer)
	CHARACTER*40 SNAME(MAXSOIL)		! IN/OUT:Soil name
	REAL DFACT(MAXLAYER)			! IN/OUT: Denitrification factor
	REAL ANIT						! Ammonium immobilised (kgN/ha)
	REAL ANIT15						! Ammonium N15 nitrified (kgN/ha)
	REAL FANIT						! Fertiliser N nitrified (kgN/ha)
	REAL FANIT15					! Fertiliser N15 nitrified (kgN/ha)
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)! IN/PUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)! IN/PUT: pH above which decomp.max.
      REAL ICFACTOR(MAXLAYER)
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: % clay content in this layer
      INTEGER LUARRAY(MAXSOIL)	! IN/OUT:Land use before equilibrium 
C
C ...Manure factors
C
      REAL FPROPX(MAXORGM)			! IN/OUT:Prop.of FYM added to top 25cm
	REAL FYMCX(MAXORGM)				! IN/OUT:Prop.of C in FYM
	REAL FYMNX(MAXORGM)				! IN/OUT:Prop.of Organic N in FYM
	REAL FYMAX(MAXORGM)				! IN/OUT:Proportion of ammonium-N in FYM
	REAL FYMWX(MAXORGM)				! IN/OUT:Amount of water added in FYM
	REAL FYMLOSS(MAXORGM)			! IN/OUT:Prop.FYM lost by volatn/timestep
	REAL FYMPOS(MAXORGM)			! IN/OUT:Prop.of FYM added to top 25cm
	REAL FYMSTART(MAXORGM)			! IN/OUT:Amount of rainfall in 1 week, 
									!		 below which volatilisation will occur
C
C Open dump file
C
      OPEN(91,FILE='Initial_CN_Pars.DAT',STATUS='UNKNOWN',ERR=111)
C
C Dump soil CN characteristics to file
C
      IF(IDUMP.EQ.1)THEN
C
C ...Retrieve soil CN characteristics from soil CN routines
C
        ISAVE=0
        CALL SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,CONVER_F,
     &             SECONDS,NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)

C
C ...Save parameters to file
C
	  WRITE(91,10)CONVER_F,ANIT,ANIT15,FANIT,FANIT15,BYRAIN,
     &             NUMSOIL
10      FORMAT(6(F15.4,2X),I4)
	  DO 100 I=1,MAXSOIL
          WRITE(91,20)SNAME(I),LUARRAY(I)
20        FORMAT(A40,2X,I4)
	    DO 200 J=1,MAXLAYER1
            WRITE(91,30)HY(I,J),BPART(I,J),HPART(I,J),BPROP(I,J),
     &                 HPROP(I,J),BIOP(I,J),CRIT(I,J),BIORATE(I,J),
     &                 HUMRATE(I,J),IOMARRAY(I,J),TOCARRAY(I,J),
     &                 PHARRAY(I,J),PHP1ARRAY(I,J),PHP2ARRAY(I,J),
     &                 CLARRAY(I,J)
30          FORMAT(15(F15.4,2X))
200       CONTINUE
100     CONTINUE	  
        DO 300 I=1,7
          WRITE(91,40)FPROPX(I),FYMCX(I),FYMNX(I),FYMAX(I),
     &                FYMWX(I),FYMLOSS(I),FYMPOS(I),FYMSTART(I)
40        FORMAT(8(F15.4,2X))
300     CONTINUE
        DO 400 I=1,MAXLAYER1
	    WRITE(91,50)HZ1(I),DFACT(I),BRATE(I),HRATE(I),
     &               DPMRATE(I),RPMRATE(I),ALPHA(I),BETA(I),GAMMA(I),
     &               DELTA(I),ICFACTOR(I)
50        FORMAT(11(F15.4,2X))
400     CONTINUE
C
C Retrieve soil CN characteristics from dump file
C
      ELSEIF(IDUMP.EQ.0)THEN
C
C ...Retrieve parameters to dump file
C
	  READ(91,10,ERR=111)CONVER_F,ANIT,ANIT15,FANIT,FANIT15,BYRAIN,
     &                    NUMSOIL
	  DO 500 I=1,MAXSOIL
          READ(91,20,ERR=111)SNAME(I),LUARRAY(I)
	    DO 600 J=1,MAXLAYER1
            READ(91,30,ERR=111)HY(I,J),BPART(I,J),HPART(I,J),BPROP(I,J),
     &                 HPROP(I,J),BIOP(I,J),CRIT(I,J),BIORATE(I,J),
     &                 HUMRATE(I,J),IOMARRAY(I,J),TOCARRAY(I,J),
     &                 PHARRAY(I,J),PHP1ARRAY(I,J),PHP2ARRAY(I,J),
     &                 CLARRAY(I,J)
600       CONTINUE
500     CONTINUE	  
        DO 700 I=1,7
          READ(91,40,ERR=111)FPROPX(I),FYMCX(I),FYMNX(I),FYMAX(I),
     &               FYMWX(I),FYMLOSS(I),FYMPOS(I),FYMSTART(I)
700     CONTINUE
        DO 800 I=1,MAXLAYER1
	    READ(91,50,ERR=111)HZ1(I),DFACT(I),BRATE(I),HRATE(I),
     &               DPMRATE(I),RPMRATE(I),ALPHA(I),BETA(I),GAMMA(I),
     &               DELTA(I),ICFACTOR(I)
800     CONTINUE
C
C ...Save soil CN characteristics for use in soil CN routines
C
        ISAVE=1
        CALL SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,CONVER_F,
     &             SECONDS,NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)
	ENDIF
      GOTO 101
C
C Errors in open or read
C
111   CONTINUE
      PRINT*,'Error in opening or reading the dump file, ',
     &       'Initial_CN_Pars.DAT'
C
C Close dump file
C
101   CONTINUE
      CLOSE(91)
	END

C
C------------------------------------------------------------
C
      SUBROUTINE GET_CTON_FALLOON(DPMCTON,RPMCTON,LUCODE)
C
C Subroutine to get C:N ratio of DPM and RPM
C
      IMPLICIT NONE
C
C Variables local to this routine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLU				! Maximum number of land use types
	PARAMETER (MAXLU=6)	

      INTEGER IL					! Layer counter
      INTEGER ILU					! Land use counter
	REAL DCTON(MAXLU)			! C:N ratio of DPM
C                                    Ara  Gra  For  Nat  Mis  SRC
	DATA (DCTON(ILU),ILU=1,MAXLU) /25.0,25.0,55.0,40.0,25.0,55.0/ 
	REAL RCTON(MAXLU)			! C:N ratio of RPM
C                                    Ara  Gra  For  Nat  Mis  SRC
	DATA (RCTON(ILU),ILU=1,MAXLU) /25.0,25.0,55.0,40.0,25.0,55.0/
C
C Variables passed to/from this subroutine
C
      REAL DPMCTON(MAXLAYER)		! OUT:DPM C:N ratio
	INTEGER LUCODE				! IN:Land use code 1=Ara;2=Gra;3=For;4=Natural
      REAL RPMCTON(MAXLAYER)		! OUT:RPM C:N ratio
C
C Get C:N ratio for this LU from array
C
      DO 100 IL=1,MAXLAYER
        DPMCTON(IL)=DCTON(LUCODE)
	  RPMCTON(IL)=RCTON(LUCODE)
100   CONTINUE
C
C Leave GET_CTON_FALLOON
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE GET_CTON_BRADBURY(DPMCTON,RPMCTON,HZ1)
C
C Subroutine to get C:N ratio of DPM and RPM
C
      IMPLICIT NONE
C
C Variables local to this routine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLU				! Maximum number of land use types
	PARAMETER (MAXLU=6)	

      INTEGER IL					! Layer counter
      INTEGER ILU					! Land use counter
C
C Variables passed to/from this subroutine
C
      REAL DPMCTON(MAXLAYER)		! OUT:DPM C:N ratio
	REAL HZ1(MAXLAYER)		   	! IN: N:C ratio for steady state
      REAL RPMCTON(MAXLAYER)		! OUT:RPM C:N ratio
C
C Get C:N ratio for this LU from array
C
      DO 100 IL=1,MAXLAYER
        DPMCTON(IL)=1.0/HZ1(IL)
	  RPMCTON(IL)=100
100   CONTINUE
C
C Leave GET_CTON_BRADBURY
C
      RETURN
      END
C
C-------------------------------------------------------------
C
	SUBROUTINE GET_DECOMP_EFF(CLAY_IL,DECOMP_EFF)
C
C Subroutine to get the decomposition efficiency
C Decomposition efficiency = (BIO + HUM) / (Total decomposition)
C
      IMPLICIT NONE
C
C Variables passes to/from this subroutine
C
      REAL CLAY_IL		! IN: Clay content of the layer (%)
      REAL DECOMP_EFF		! OUT: Decomposition efficiency
C
C Calculate decomp.eff.from equation from Coleman & Jenkinson, 1996.
C
      DECOMP_EFF=1/(1+1.67*(1.85+1.6*EXP(-0.0786*CLAY_IL)))   
C
C Leave GET_DECOMP_EFF
C
      END
C
C-------------------------------------------------------------
C
	SUBROUTINE GETEQUIVALENTS(SOILN_IL,AMMN_IL,MOBDOC_IL,SOILW_IL,
     &                          CNO3_IL,CNH4_IL,CMOBDOC_IL)
C
C Subroutine to convert kg/ha/mm to equivalents/m3
C	 
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      REAL MWNO3					! Molecular weight of nitrate g/mol
	DATA MWNO3 /62.01/
	REAL MWNH4					! Molecular weight of ammonium g/mol
	DATA MWNH4 /18.042/

C
C Variables passed to/from this subroutine
C
	REAL AMMN_IL				! IN:Soil ammonium-N (kgN/ha/layer)
	REAL CMOBDOC_IL				! OUT: Concentration of mobile DOC in this layer (mg/l)
	REAL CNH4_IL				! OUT: Concentration of ammonium in this layer (eq/m3)
      REAL CNO3_IL				! OUT: Concentration of nitrate in this layer (eq/m3)
	REAL MOBDOC_IL				! IN:Mobile DOC in each layer (kgC/ha)
	REAL SOILN_IL				! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOILW_IL				! IN:Available water (mm/layer)
      
C
C Convert to equivalents
C
      IF(SOILW_IL.GT.0)THEN
	  CNO3_IL=SOILN_IL*100/(SOILW_IL*MWNO3)
	  CNH4_IL=AMMN_IL*100/(SOILW_IL*MWNH4)
	  CMOBDOC_IL=MOBDOC_IL*10000/SOILW_IL
	ELSEIF(SOILW_IL.LE.0)THEN
	  CNO3_IL=0
	  CNH4_IL=0
	  CMOBDOC_IL=0
	ENDIF
	END
C
C--------------------------------------------------------
C
      SUBROUTINE GET_PLANT_DIST(PI_ANN,PI_CEQ_MON,LUCODE) 
C
C Subroutine to get plant input distribution
C
      IMPLICIT NONE
C
C Variables local to this routine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER MAXDEPTH			! Maximum depth of the soil profile
	DATA MAXDEPTH /300/
	INTEGER MAXLU				! Maximum number of land use types
	PARAMETER (MAXLU=10)	
	INTEGER MAXPILAY			! Maximum number of PI layers
	PARAMETER (MAXPILAY=3)
	REAL PLADD(12,MAXLU,MAXPILAY)	! Plant input (t/ha/month)
      INTEGER IL					! Local layer counter
      INTEGER IMON				! Local month counter
      INTEGER M					! Local month counter
	REAL DEPTH					! Depth
	INTEGER ILAY				! Counter for PI layer number
	REAL PISPL					! Split of PI layer into SUNDIAL LAYERS
	REAL PL_TOT					! Total standard plant input 
C
C Variables passed to/from this routine
C ...Plant factors
C
      REAL PI_ANN						! IN:Annual plant input C (kgC/ha/year)
	REAL PI_CEQ_MON(12,MAXLAYER)	! OUT:Equilibrium plant C input each 
									!         month in each layer (kgC/ha/month/layer)
	INTEGER LUCODE					! IN:Code for land use type


C Initialise variables
C (1) Arable crops
C (1.1) 0-30cm
C
	DATA (PLADD(M,1,1),M=1,12)/0.00,0.00,0.30,0.30,0.30,
     &                           0.60,1.87,0.00,0.00,0.00,0.00,0.00/
C
C (1.2) 30-100cm
C
	DATA (PLADD(M,1,2),M=1,12)/0.00,0.00,0.30,0.30,0.30,
     &                           0.30,1.17,0.00,0.00,0.00,0.00,0.00/
C
C (1.3) >100cm
C
	DATA (PLADD(M,1,3),M=1,12)/0.00,0.00,0.30,0.30,0.30,
     &                           0.30,1.17,0.00,0.00,0.00,0.00,0.00/
C
C (2) Grassland
C (2.1) 0-30cm
C
      DATA (PLADD(M,2,1),M=1,12)/0.25,0.25,0.25,0.25,0.25,0.25,
     &                         0.25,0.89,0.25,0.25,0.25,0.25/
C
C (2.2) 30-100cm
C
      DATA (PLADD(M,2,2),M=1,12)/0.2083,0.2083,0.2083,0.2083,0.2083,
     &                           0.2083,0.2083,0.2083,0.2083,0.2083,
     &                           0.2083,0.2083/
C
C (2.3) >100cm
C
      DATA (PLADD(M,2,3),M=1,12)/0.2083,0.2083,0.2083,0.2083,0.2083,
     &                           0.2083,0.2083,0.2083,0.2083,0.2083,
     &                           0.2083,0.2083/
C
C (3) Forestry
C (3.1) 0-30cm
C
      DATA (PLADD(M,3,1),M=1,12)/0.2375,0.2375,0.2375,0.2375,0.2375,
     &                           0.2375,0.2375,0.2375,0.2375,0.2375,
     &                           0.2375,0.2375/
C
C (3.2) 30-100cm
C
      DATA (PLADD(M,3,2),M=1,12)/0.2375,0.2375,0.2375,0.2375,0.2375,
     &                           0.2375,0.2375,0.2375,0.2375,0.2375,
     &                           0.2375,0.2375/
C
C (3.3) >100cm
C
      DATA (PLADD(M,3,3),M=1,12)/0.2375,0.2375,0.2375,0.2375,0.2375,
     &                           0.2375,0.2375,0.2375,0.2375,0.2375,
     &                           0.2375,0.2375/
C
C (4) Semi-natural
C (4.1) 0-30cm
C
      DATA (PLADD(M,4,1),M=1,12)/0.475,0.475,0.475,0.475,0.475,0.475,
     &                           0.475,0.475,0.475,0.475,0.475,0.475/
C
C (4.2) 30-100cm
C
      DATA (PLADD(M,4,2),M=1,12)/0.25,0.25,0.25,0.25,0.25,0.25,
     &                           0.25,0.25,0.25,0.25,0.25,0.25/
C>>> Temp change - Add for accumulating peats JUS 17/04/08 >>>>
C      DATA (PLADD(M,4,2),M=1,12)/0,0,0,0,0,0,0,0,0,0,0,0/
c<<< Temp change JUS 17/04/08 >>>>
C
C (4.3) >100cm
C
      DATA (PLADD(M,4,3),M=1,12)/0.25,0.25,0.25,0.25,0.25,0.25,
     &                           0.25,0.25,0.25,0.25,0.25,0.25/
C>>> Temp change - Add for accumulating peats JUS 17/04/08 >>>>
c      DATA (PLADD(M,4,3),M=1,12)/0,0,0,0,0,0,0,0,0,0,0,0/
c<<< Temp change JUS 17/04/08 >>>>
C
C (5) Miscanthus
C (5.1) 0-30cm
C
      DATA (PLADD(M,5,1),M=1,12)/0.25,0.25,0.25,0.25,0.25,0.25,
     &                         0.25,0.89,0.25,0.25,0.25,0.25/
C
C (5.2) 30-100cm
C
      DATA (PLADD(M,5,2),M=1,12)/0.2083,0.2083,0.2083,0.2083,0.2083,
     &                           0.2083,0.2083,0.2083,0.2083,0.2083,
     &                           0.2083,0.2083/
C
C (5.3) >100cm
C
      DATA (PLADD(M,5,3),M=1,12)/0.2083,0.2083,0.2083,0.2083,0.2083,
     &                           0.2083,0.2083,0.2083,0.2083,0.2083,
     &                           0.2083,0.2083/
C
C (6) SRC
C (6.1) 0-30cm
C
      DATA (PLADD(M,6,1),M=1,12)/0.2375,0.2375,0.2375,0.2375,0.2375,
     &                           0.2375,0.2375,0.2375,0.2375,0.2375,
     &                           0.2375,0.2375/
C
C (6.2) 30-100cm
C
      DATA (PLADD(M,6,2),M=1,12)/0.2375,0.2375,0.2375,0.2375,0.2375,
     &                           0.2375,0.2375,0.2375,0.2375,0.2375,
     &                           0.2375,0.2375/
C
C (6.3) >100cm
C
      DATA (PLADD(M,6,3),M=1,12)/0.2375,0.2375,0.2375,0.2375,0.2375,
     &                           0.2375,0.2375,0.2375,0.2375,0.2375,
     &                           0.2375,0.2375/
C (7) Sugar beet
C (7.1) 0-30cm
C
	DATA (PLADD(M,7,1),M=1,12)/0.00,0.00,0.30,0.30,0.30,
     &                           0.60,1.87,0.00,0.00,0.00,0.00,0.00/
C
C (7.2) 30-100cm
C
	DATA (PLADD(M,7,2),M=1,12)/0.00,0.00,0.30,0.30,0.30,
     &                           0.30,1.17,0.00,0.00,0.00,0.00,0.00/
C
C (7.3) >100cm
C
	DATA (PLADD(M,7,3),M=1,12)/0.00,0.00,0.30,0.30,0.30,
     &                           0.30,1.17,0.00,0.00,0.00,0.00,0.00/
C (8) Oilseed rape
C (8.1) 0-30cm
C
	DATA (PLADD(M,8,1),M=1,12)/0.00,0.00,0.30,0.30,0.30,
     &                           0.60,1.87,0.00,0.00,0.00,0.00,0.00/
C
C (8.2) 30-100cm
C
	DATA (PLADD(M,8,2),M=1,12)/0.00,0.00,0.30,0.30,0.30,
     &                           0.30,1.17,0.00,0.00,0.00,0.00,0.00/
C
C (8.3) >100cm
C
	DATA (PLADD(M,8,3),M=1,12)/0.00,0.00,0.30,0.30,0.30,
     &                           0.30,1.17,0.00,0.00,0.00,0.00,0.00/
C (9) SRF
C (9.1) 0-30cm
C
      DATA (PLADD(M,9,1),M=1,12)/0.2375,0.2375,0.2375,0.2375,0.2375,
     &                           0.2375,0.2375,0.2375,0.2375,0.2375,
     &                           0.2375,0.2375/
C
C (9.2) 30-100cm
C
      DATA (PLADD(M,9,2),M=1,12)/0.2375,0.2375,0.2375,0.2375,0.2375,
     &                           0.2375,0.2375,0.2375,0.2375,0.2375,
     &                           0.2375,0.2375/
C
C (9.3) >100cm
C
      DATA (PLADD(M,9,3),M=1,12)/0.2375,0.2375,0.2375,0.2375,0.2375,
     &                           0.2375,0.2375,0.2375,0.2375,0.2375,
     &                           0.2375,0.2375/        
C (10) Wheat
C (10.1) 0-30cm
C
	DATA (PLADD(M,10,1),M=1,12)/0.00,0.00,0.30,0.30,0.30,
     &                           0.60,1.87,0.00,0.00,0.00,0.00,0.00/
C
C (10.2) 30-100cm
C
	DATA (PLADD(M,10,2),M=1,12)/0.00,0.00,0.30,0.30,0.30,
     &                           0.30,1.17,0.00,0.00,0.00,0.00,0.00/
C
C (10.3) >100cm
C
	DATA (PLADD(M,10,3),M=1,12)/0.00,0.00,0.30,0.30,0.30,
     &                           0.30,1.17,0.00,0.00,0.00,0.00,0.00/
C
C For each layer down the profile
C
      PL_TOT=0
      DO 100 ILAY=1,MAXPILAY
	  DO 200 IMON=1,12
	    PL_TOT=PL_TOT+PLADD(IMON,LUCODE,ILAY)
200     CONTINUE
100   CONTINUE
      IF(PI_ANN.LT.1)THEN				! THIS ENABLES SPINUP TO IDENTIFY WHEN ICOVER EQUALS 1 OR ZERO
		PI_ANN=1						! IT IS THEN OVER WRITEN LATER WHEN PIANN IS GREATER THAN 1
      ENDIF

C
C For each layer down the profile
C
      DO 300 IL=1,MAXLAYER1
C
C ...work out depth
C
	  DEPTH=IL*MAXDEPTH/(MAXLAYER1*1.)    
C
C ...set layer for SOM calculation and split of SOM layer into SUNDIAL layers (SOMSPL)
C
        IF(DEPTH.LE.30)THEN
	    ILAY=1
	    PISPL=(MAXDEPTH/MAXLAYER1)/30.
	  ELSEIF(DEPTH.GT.30.AND.DEPTH.LE.100)THEN
	    ILAY=2
	    PISPL=(MAXDEPTH/MAXLAYER1)/(100.-30.)
	  ELSEIF(DEPTH.GT.100)THEN
	    ILAY=3
	    PISPL=(MAXDEPTH/MAXLAYER1)/(MAXDEPTH-100.)
	  ENDIF
C
C ...distribute plant input down the soil profile and across the months
C
        DO 400 IMON=1,12
          PI_CEQ_MON(IMON,IL)=PISPL*PLADD(IMON,LUCODE,ILAY)/PL_TOT
          PI_CEQ_MON(IMON,IL)=PI_ANN*PI_CEQ_MON(IMON,IL)
400     CONTINUE
300   CONTINUE
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE GET_SUNDIAL_SOILPARS(NSOIL,TOC,IOM,CRIT,ILU,
     &                                CLAY,BULKDENS,SPARMODEL)
C
C Subroutine to get soil parameters
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXSOIL				! Max.no.of soil types
	INTEGER MAXORGM				! Max.no.of org.manure applications allowed
	PARAMETER (MAXSOIL=50)
	PARAMETER (MAXORGM=52)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER IL					! Layer counter
	INTEGER SPARFILE			! Soil parameters read in from file
	INTEGER SPARCALC			! Soil parameters calculated from TOC etc
	DATA SPARFILE,SPARCALC /1,2/
C
C Variables passed to/from this subroutine
C
      INTEGER ISAVE				! OUT:Code to save or retrieve variables
C
C ...Timing factors
C
      REAL CONVER_F				! IN/OUT:Conversion to correct timestep
	REAL SECONDS				! IN:Number of seconds in one timestep
C
C ...Model descriptors
C
	INTEGER SPARMODEL			! IN:Soil parameter model 
  	REAL ICFACTOR(MAXLAYER)			! IN:Adjustment in rate needed to	
									!    achieve measured NPP and TOC	
C
C ...Soil factors
C
	INTEGER NSOIL				! IN:Soil code number
      INTEGER NUMSOIL				! IN/OUT:Number of soils defined
	REAL HZ1(MAXLAYER)				! IN/OUT: N:C ratio for steady state
	REAL TOC(MAXLAYER)			! IN:Meas.total organic C (kgC/ha/layer)
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
      REAL CLAY(MAXLAYER)			! IN:Clay content of the layer (%)
	REAL BULKDENS(MAXLAYER)		! IN:Bulk density of the layer (g/cm3)
	REAL HY(MAXSOIL,MAXLAYER)	! IN/OUT:Stable N:C ratio of BIO & HUM pools
	REAL BPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL HPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL BPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from biomass decompositn
	REAL HPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from humus decompositn
	REAL BIOP(MAXSOIL,MAXLAYER)		! IN/OUT: BIO/TOC 
	REAL CRIT(MAXSOIL,MAXLAYER)		! IN/OUT: Minimum nitrate level (Sozanska)
	REAL BIORATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for BIO decompn/yr
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for HUM decompn/yr
	REAL IOMARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium IOM (kgC/ha/layer)
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium TOC (kgC/ha/layer)
	CHARACTER*40 SNAME(MAXSOIL)	! IN/OUT:Soil name
	REAL DFACT(MAXLAYER)		! IN/OUT: Denitrification factor
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: pH above which decomp.max.
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: % clay content in this layer
      INTEGER LUARRAY(MAXSOIL)	! IN/OUT:Land use before equilibrium 
	INTEGER ILU						! IN:Counter for land use 
	REAL ALPHA(MAXLAYER)			! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL BETA(MAXLAYER)				! IN/OUT: Prop.HUM produced on BIO decompn
	REAL GAMMA(MAXLAYER)			! IN/OUT: Prop.BIO produced on HUM decompn 
	REAL DELTA(MAXLAYER)			! IN/OUT: Prop.HUM produced on HUM decompn 
	REAL DPMRATE(MAXLAYER)			! IN/OUT: Rate constant for DPM decompn
	REAL RPMRATE(MAXLAYER)			! IN/OUT: Rate constant for RPM decompn
	REAL BRATE(MAXLAYER)			! IN/OUT: Rate constant for HUM decompn
	REAL HRATE(MAXLAYER)			! IN/OUT: Rate constant for BIO decompn
	REAL ANIT					! IN/OUT:NH4-N nitrified(kgN/ha)
	REAL ANIT15					! IN/OUT:NH4-N15 nitrified (kgN/ha)
	REAL FANIT					! IN/OUT:Fertiliser NH4-N nitrified(kgN/ha)
	REAL FANIT15				! IN/OUT:Fertiliser NH4-N15 nitrified(kgN/ha)
C
C ...Fertiliser factors
C
	REAL BYRAIN					! IN/OUT:Rain lost by bypass flow mm/timestep	
C
C ...Manure factors
C
      REAL FPROPX(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
	REAL FYMCX(MAXORGM)			! IN/OUT:Proportion of C in FYM
	REAL FYMNX(MAXORGM)			! IN/OUT:Proportion of Organic N in FYM
	REAL FYMAX(MAXORGM)			! IN/OUT:Proportion of ammonium-N in FYM
	REAL FYMWX(MAXORGM)			! IN/OUT:Amount of water added in FYM
	REAL FYMLOSS(MAXORGM)		! IN/OUT:Prop.FYM lost by volatn/timestep
	REAL FYMPOS(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
      REAL FYMSTART(MAXORGM)		! IN/OUT:Amount of rainfall in 1 week, 
								!		 below which volatilisation will occur
C
C Retrieve soil CN characteristics 
C
      ISAVE=0
      CALL SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,CONVER_F,
     &             SECONDS,NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)
C
C Set Parameters
C
      IF(SPARMODEL.EQ.SPARFILE)THEN
        CALL PAR_CN(NUMSOIL,HY,BPART,BPROP,HPROP,BIOP,CRIT,
     &           BIORATE,HUMRATE,BULKDENS,
     &           TOCARRAY,IOMARRAY,LUARRAY,CLARRAY,SNAME,DFACT,
     &           FYMSTART,FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,FYMLOSS,FYMPOS,
     &           PHARRAY,PHP1ARRAY,PHP2ARRAY)

	ELSEIF(SPARMODEL.EQ.SPARCALC)THEN
	  CALL MAKE_PAR_CN(NSOIL,TOC,IOM,ILU,CLAY,BULKDENS,
     &                   HY,BPART,BPROP,HPART,HPROP,
     &				   BIOP,CRIT,BIORATE,HUMRATE,TOCARRAY,IOMARRAY,
     &				   SNAME,CLARRAY,LUARRAY,
     &				   PHARRAY,PHP1ARRAY,PHP2ARRAY,DFACT,
     &				   FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,FYMLOSS,FYMPOS,
     &				   FYMSTART)
	ENDIF
      ILU=LUARRAY(NSOIL)
C
C Save soil CN characteristics 
C
      ISAVE=1
      CALL SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,CONVER_F,
     &             SECONDS,NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE INIT_SUNDIAL_SOILCN(
C INPUTS: Time, environment	
     &                               ICMODEL,INMODEL,DOCMODEL,EQMODEL,
     &                               SECONDS,NSOIL,
C OUTPUTS: Soil C&N     
     &                               TOC,IOM,HZ1,CRIT,
C INPUTS/OUTPUTS: Soil C&N
     &                               SOILN,SOIL15,AMMN,AMMN15,
     &                               DPMCARB0,DPMNIT0,DPMNLAB0,
     &                               RPMCARB0,RPMNIT0,RPMNLAB0,
     &                               BCARB0,BNIT0,BNLAB0,
     &                               HCARB0,HNIT0,HNLAB0,
     &                               MOBDOC,MOBDON,SPARMODEL,
     &                               ILU,CLAY,BULKDENS,PI_CEQ_MON,
     &                               LTA_AWC,LTA_TEMP,ITFUNC,IMFUNC,
     &                               WMAX,WSAT,
     &                               ICFACTOR,DPMCTON,RPMCTON,SOILPH)
C
C Subroutine to initialise soil C and N
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLU				! Maximum number of land use types
	PARAMETER (MAXLU=6)	
	INTEGER MAXSOIL				! Max.no.of soil types
	INTEGER MAXORGM				! Max.no.of org.manure applications allowed
	PARAMETER(MAXSOIL=50,MAXORGM=52)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER IL					! Layer counter
	INTEGER MAXFERTSTEPS		! Maximum number of timesteps 
								!	over which fertiliser is added
      PARAMETER (MAXFERTSTEPS=30*60*24*7*3)	! Timesteps = 3 wks (max.30 secs)
C
C ...Model descriptors
C
	INTEGER DOC_OFF				! DOC model off
	INTEGER DOC_ON				! DOC model on
	DATA DOC_ON,DOC_OFF /1,2/
	INTEGER EQNPP				! Model initialised using plant inputs 
								!    from measured NPP
	INTEGER EQTOC				! Model initialised using measured TOC
	INTEGER EQNPPTOC			! Model initialised using both plant inputs
								!    from measured NPP and measured TOC
	INTEGER EQHILLIER,EQJONES	! Model initialised using TOC and
								!    the Hillier solver to get C pools and plant inputs
	DATA EQNPP,EQTOC,EQNPPTOC,EQHILLIER,EQJONES /1,2,3,4,5/
      INTEGER SPARFILE			! Soil parameters read in from file
      INTEGER SPARCALC			! Soil parameters calculated from TOC etc	
      DATA SPARFILE,SPARCALC /1,2/										
C
C Variables passed to/from this subroutine
C
      INTEGER ISAVE				! OUT:Code to save or retrieve variables
C
C ...Weather factors
C
	REAL LTA_AWC(12,MAXLAYER)	! IN:Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)	! IN:Average air temp this month (deg.C)
C
C ...Model descriptors
C
	INTEGER ICMODEL				! OUT:Type of C initialisation model 
	INTEGER ICFIXED				! Initialisation of C is fixed
	INTEGER ICROTHCEQ			! Initialisation of C by RothC 
								!  equilibrium run
	DATA ICFIXED,ICROTHCEQ /1,2/
	INTEGER INMODEL				! OUT:Type of N initialisation model
	INTEGER INSTABLECN			! Initisalisation of N by assuming 
								!  steady state (after Bradbury)
	INTEGER INPASSCN			! Initialisation of N by passing the
								!  C:N ratio of DPM and RPM
	DATA INSTABLECN,INPASSCN /1,2/
	INTEGER DOCMODEL			! IN:DOC model (on or off)
	INTEGER EQMODEL				! IN:Type of equilibrium run (NPP,TOC or both)
	INTEGER ITFUNC				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IMFUNC				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
      INTEGER SPARMODEL			! IN:Soil parameter model
C
C ...Timing factors
C
      REAL CONVER_F				! IN/OUT:Conversion to correct timestep
	REAL SECONDS				! IN:Number of seconds in one timestep
C
C ...Soil factors
C
	REAL ALPHA(MAXLAYER)			! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
	REAL ANIT					! IN/OUT:NH4-N nitrified(kgN/ha)
	REAL ANIT15					! IN/OUT:NH4-N15 nitrified (kgN/ha)
	REAL BCARB0(MAXLAYER)		! IN:C in soil biomass (kgC/ha/layer)
	REAL BETA(MAXLAYER)				! IN/OUT: Prop.HUM produced on BIO decompn
	REAL BIOP(MAXSOIL,MAXLAYER)		! IN/OUT: BIO/TOC 
	REAL BIORATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for BIO decompn/yr
      REAL BNIT0(MAXLAYER)		! IN:N in soil biomass (kgN/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN:N15 in soil humus (kgN15/ha/layer)
	REAL BPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL BPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from biomass decompositn
	REAL BRATE(MAXLAYER)			! IN/OUT: Rate constant for HUM decompn
	REAL BULKDENS(MAXLAYER)		! IN:Bulk density of the layer (g/cm3)
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: % clay content in this layer
      REAL CLAY(MAXLAYER)			! IN:Clay content of the layer (%)
	REAL CRIT(MAXSOIL,MAXLAYER)		! IN/OUT: Minimum nitrate level (Sozanska)
	REAL DELTA(MAXLAYER)			! IN/OUT: Prop.HUM produced on HUM decompn 
	REAL DFACT(MAXLAYER)			! IN/OUT: Denitrification factor
      REAL DPM_RPM_RATIO			! OUT:Ratio of DPM:RPM. If set to 
								!     zero, this will be worked out using a call to SET_DPMRPMRATIO
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
      REAL DPMCTON(MAXLAYER)		! IN/OUT:Ratio of C:N in DPM
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL DPMRATE(MAXLAYER)			! IN/OUT: Rate constant for DPM decompn
	REAL FANIT					! IN/OUT:Fertiliser NH4-N nitrified(kgN/ha)
	REAL FANIT15				! IN/OUT:Fertiliser NH4-N15 nitrified(kgN/ha)
	REAL GAMMA(MAXLAYER)			! IN/OUT: Prop.BIO produced on HUM decompn 
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN:N in soil humus (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN:N15 in soil biomass (kgN15/ha/layer)
	REAL HPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL HPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from humus decompositn
	REAL HRATE(MAXLAYER)			! IN/OUT: Rate constant for BIO decompn
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for HUM decompn/yr
	REAL HY(MAXSOIL,MAXLAYER)	! IN/OUT:Stable N:C ratio of BIO & HUM pools
	REAL HZ1(MAXLAYER)				! IN/OUT: N:C ratio for steady state
  	REAL ICFACTOR(MAXLAYER)			! IN:Adjustment in rate needed to	
									!    achieve measured NPP and TOC	
	INTEGER ILU						! IN:Counter for land use 
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
	REAL IOMARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium IOM (kgC/ha/layer)
      INTEGER LUARRAY(MAXSOIL)	! IN/OUT:Land use before equilibrium 
	INTEGER NSOIL				! IN:Soil code number
      INTEGER NUMSOIL				! IN/OUT:Number of soils defined
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)! IN/OUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)! IN/OUT: pH above which decomp.max.
	REAL PI_CEQ_MON(12,MAXLAYER)	! IN/OUT: Equilibrium plant C input each 
									!         month to this layer (kgC/ha/month/layer)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
      REAL RPMCTON(MAXLAYER)		! IN/OUT:Ratio of C:N in RPM
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL RPMRATE(MAXLAYER)			! IN/OUT: Rate constant for RPM decompn
	CHARACTER*40 SNAME(MAXSOIL)	! IN/OUT:Soil name
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOILPH(MAXLAYER)		! IN/OUT:pH of soil in this layer
	REAL TOC(MAXLAYER)			! IN:Total organic C (kgC/ha/layer)
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium TOC (kgC/ha/layer)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
C
C ...Manure factors
C
      REAL FPROPX(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
	REAL FYMCX(MAXORGM)			! IN/OUT:Proportion of C in FYM
	REAL FYMNX(MAXORGM)			! IN/OUT:Proportion of Organic N in FYM
	REAL FYMAX(MAXORGM)			! IN/OUT:Proportion of ammonium-N in FYM
	REAL FYMWX(MAXORGM)			! IN/OUT:Amount of water added in FYM
	REAL FYMLOSS(MAXORGM)		! IN/OUT:Prop.FYM lost by volatn/timestep
	REAL FYMPOS(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
	REAL FYMSTART(MAXORGM)		! IN/OUT:Amount of rainfall in 1 week, 
								!		 below which volatilisation will occur
C
C ...Fertiliser factors
C
	REAL BYRAIN					! IN/OUT:Rain lost by bypass flow mm/timestep	
	REAL FERTADD(MAXFERTSTEPS,4)	! INOUT:Fertiliser N added in this time 
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
	REAL FERTADD15(MAXFERTSTEPS,4)	! INOUT:Fertiliser N15 added in this time
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
C
C ...Dissolved organic matter factors
C
	REAL MOBDOC(MAXLAYER)		! IN:Mobile DOC in each layer (kgC/ha)
	REAL MOBDON(MAXLAYER)		! IN:Mobile DON in each layer (kgN/ha)
C
C Get soil parameters
C
  	CALL GET_SUNDIAL_SOILPARS(NSOIL,TOC,IOM,CRIT,ILU,
     &                          CLAY,BULKDENS,SPARMODEL)

C
C Retrieve soil CN characteristics 
C
      ISAVE=0
      CALL SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,CONVER_F,
     &             SECONDS,NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)
	CALL SAVE_FERTADD(SECONDS,FERTADD,FERTADD15,ISAVE)
C
C Set time factors according to SECONDS
C
      CALL SETTIME_CN(SECONDS,CONVER_F)
C
C Set starting value of organic matter
C
      CALL STARTORG(NSOIL,TOC,TOCARRAY,IOM,IOMARRAY)
C
C Set BIO and HUM decompostion parameters
C
      CALL OMPROP(NSOIL,HY,HZ1,ALPHA,BETA,GAMMA,DELTA,
     &                  BPART,HPART,BPROP,HPROP,PHARRAY)
C
C Set NO3 and NH4 residual
C
      CALL SOILFIX(NSOIL,CRIT,AMMN,AMMN15,SOILN,SOIL15)
C
C Set the residual N and C in the RO,BIO and HUM pools
C
      IF(INMODEL.EQ.INSTABLECN)THEN
        CALL GET_CTON_BRADBURY(DPMCTON,RPMCTON,HZ1)
	ENDIF
	IF(ICMODEL.EQ.ICROTHCEQ)THEN
        CALL SETORGM_EQRUN(NSOIL,DPMRATE,RPMRATE,TOC,IOM,CONVER_F,
     &                   BRATE,BIORATE,HRATE,HUMRATE,HZ1,HY,
     &                   BCARB0,HCARB0,DPMCARB0,RPMCARB0,
     &                   BNIT0,HNIT0,DPMNIT0,RPMNIT0,
     &                   HNLAB0,BNLAB0,DPMNLAB0,RPMNLAB0,
     &                   LUARRAY,CLARRAY,
     &                   PHARRAY,PHP1ARRAY,PHP2ARRAY,
     &                   EQMODEL,PI_CEQ_MON,ICFACTOR,LTA_AWC,LTA_TEMP,
     &                   ITFUNC,IMFUNC,WMAX,WSAT,
     &                   DPM_RPM_RATIO,DPMCTON,RPMCTON,ALPHA,BETA,CLAY)	
      ELSEIF(ICMODEL.EQ.ICFIXED)THEN
        CALL SETORGM_FIXED(NSOIL,DPMRATE,RPMRATE,TOC,IOM,CONVER_F,
     &                   BRATE,BIORATE,HRATE,HUMRATE,HZ1,HY,
     &                   BCARB0,HCARB0,DPMCARB0,RPMCARB0,
     &                   BNIT0,HNIT0,DPMNIT0,RPMNIT0,
     &                   HNLAB0,BNLAB0,DPMNLAB0,RPMNLAB0,ICFACTOR,
     &                   DPM_RPM_RATIO,DPMCTON,RPMCTON)	
	ENDIF 
C
C Set residual DOC content
C
      IF(DOCMODEL.EQ.DOC_ON)CALL SETDOC(NSOIL,MOBDOC,MOBDON)
C
C Set initial fertiliser additions
C
      CALL SETFERT(FERTADD,FERTADD15)
C
C Pass soil pH out of routine
C
      DO 100 IL=1,MAXLAYER
	  SOILPH(IL)=PHARRAY(NSOIL,IL)
100   CONTINUE
C
C Save soil CN characteristics 
C
      ISAVE=1
      CALL SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,CONVER_F,
     &             SECONDS,NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)
	CALL SAVE_FERTADD(SECONDS,FERTADD,FERTADD15,ISAVE)
C
C Leave INIT_CN
C
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE INIT_SUNDIAL_SOILCN_FIXED(
C INPUTS: Time, environment	
     &                               INMODEL,DOCMODEL,SECONDS,
     &                               NSOIL,
C OUTPUTS: Soil C&N     
     &                               TOC,IOM,HZ1,CRIT,
C INPUTS/OUTPUTS: Soil C&N
     &                               SOILN,SOIL15,AMMN,AMMN15,
     &                               DPMCARB0,DPMNIT0,DPMNLAB0,
     &                               RPMCARB0,RPMNIT0,RPMNLAB0,
     &                               BCARB0,BNIT0,BNLAB0,
     &                               HCARB0,HNIT0,HNLAB0,
     &                               MOBDOC,MOBDON,SPARMODEL,
     &                               ILU,CLAY,BULKDENS,ICFACTOR,
     &                               DPM_RPM_RATIO,DPMCTON,RPMCTON,
     &                               SOILPH)	 
C
C Subroutine to initialise soil C and N
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLU				! Maximum number of land use types
	PARAMETER (MAXLU=6)	
	INTEGER MAXSOIL				! Max.no.of soil types
	INTEGER MAXORGM				! Max.no.of org.manure applications allowed
	PARAMETER(MAXSOIL=50,MAXORGM=52)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER IL					! Layer counter
	INTEGER MAXFERTSTEPS		! Maximum number of timesteps 
								!	over which fertiliser is added
      PARAMETER (MAXFERTSTEPS=30*60*24*7*3)	! Timesteps = 3 wks (max.30 secs)
C
C ...Model descriptors
C
	INTEGER INSTABLECN			! Initisalisation of N by assuming 
								!  steady state (after Bradbury)
	INTEGER INPASSCN			! Initialisation of N by passing the
								!  C:N ratio of DPM and RPM
	DATA INSTABLECN,INPASSCN /1,2/
	INTEGER DOC_OFF				! DOC model off
	INTEGER DOC_ON				! DOC model on
	DATA DOC_ON,DOC_OFF /1,2/
	INTEGER EQNPP				! Model initialised using plant inputs 
								!    from measured NPP
	INTEGER EQTOC				! Model initialised using measured TOC
	INTEGER EQNPPTOC			! Model initialised using both plant inputs
								!    from measured NPP and measured TOC
	INTEGER EQHILLIER,EQJONES	! Model initialised using TOC and
								!    the Hillier solver to get C pools and plant inputs
	DATA EQNPP,EQTOC,EQNPPTOC,EQHILLIER,EQJONES /1,2,3,4,5/
      INTEGER SPARFILE			! Soil parameters read in from file
      INTEGER SPARCALC			! Soil parameters calculated from TOC etc	
      DATA SPARFILE,SPARCALC /1,2/										
C
C Variables passed to/from this subroutine
C
      INTEGER ISAVE				! OUT:Code to save or retrieve variables
C
C ...Model descriptors
C
	INTEGER DOCMODEL			! IN:DOC model (on or off)
	INTEGER INMODEL				! IN:Type of N initialisation model
      INTEGER SPARMODEL			! IN:Soil parameter model
C
C ...Timing factors
C
      REAL CONVER_F				! IN/OUT:Conversion to correct timestep
	REAL SECONDS				! IN:Number of seconds in one timestep
C
C ...Soil factors
C
      REAL DPM_RPM_RATIO			! IN/OUT:Ratio of DPM:RPM. If set to 
								!     zero, this will be worked out using a call to SET_DPMRPMRATIO
      REAL DPMCTON(MAXLAYER)		! IN/OUT:Ratio of C:N in DPM
      REAL RPMCTON(MAXLAYER)		! IN/OUT:Ratio of C:N in RPM
	REAL TOC(MAXLAYER)			! IN:Total organic C (kgC/ha/layer)
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
      REAL CLAY(MAXLAYER)			! IN:Clay content of the layer (%)
	REAL BULKDENS(MAXLAYER)		! IN:Bulk density of the layer (g/cm3)
	INTEGER NSOIL				! IN:Soil code number
      INTEGER NUMSOIL				! IN/OUT:Number of soils defined
	CHARACTER*40 SNAME(MAXSOIL)	! IN/OUT:Soil name
	REAL ANIT					! IN/OUT:NH4-N nitrified(kgN/ha)
	REAL ANIT15					! IN/OUT:NH4-N15 nitrified (kgN/ha)
	REAL FANIT					! IN/OUT:Fertiliser NH4-N nitrified(kgN/ha)
	REAL FANIT15				! IN/OUT:Fertiliser NH4-N15 nitrified(kgN/ha)
	REAL BCARB0(MAXLAYER)		! IN:C in soil biomass (kgC/ha/layer)
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
      REAL BNIT0(MAXLAYER)		! IN:N in soil biomass (kgN/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN:N in soil humus (kgN/ha/layer)
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN:N15 in soil biomass (kgN15/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN:N15 in soil humus (kgN15/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILPH(MAXLAYER)		! IN/OUT:pH of soil in this layer
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)! IN/OUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)! IN/OUT: pH above which decomp.max.
      INTEGER LUARRAY(MAXSOIL)	! IN/OUT:Land use before equilibrium 
	REAL HY(MAXSOIL,MAXLAYER)	! IN/OUT:Stable N:C ratio of BIO & HUM pools
	REAL BPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL HPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL BPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from biomass decompositn
	REAL HPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from humus decompositn
	REAL BIOP(MAXSOIL,MAXLAYER)		! IN/OUT: BIO/TOC 
	REAL CRIT(MAXSOIL,MAXLAYER)		! IN/OUT: Minimum nitrate level (Sozanska)
	REAL BIORATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for BIO decompn/yr
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for HUM decompn/yr
	REAL ALPHA(MAXLAYER)			! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL BETA(MAXLAYER)				! IN/OUT: Prop.HUM produced on BIO decompn
	REAL GAMMA(MAXLAYER)			! IN/OUT: Prop.BIO produced on HUM decompn 
	REAL DELTA(MAXLAYER)			! IN/OUT: Prop.HUM produced on HUM decompn 
	REAL HZ1(MAXLAYER)				! IN/OUT: N:C ratio for steady state
	REAL DFACT(MAXLAYER)			! IN/OUT: Denitrification factor
	REAL DPMRATE(MAXLAYER)			! IN/OUT: Rate constant for DPM decompn
	REAL RPMRATE(MAXLAYER)			! IN/OUT: Rate constant for RPM decompn
	REAL BRATE(MAXLAYER)			! IN/OUT: Rate constant for HUM decompn
	REAL HRATE(MAXLAYER)			! IN/OUT: Rate constant for BIO decompn
	REAL IOMARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium IOM (kgC/ha/layer)
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium TOC (kgC/ha/layer)
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: % clay content in this layer
	INTEGER ILU						! IN:Counter for land use 
  	REAL ICFACTOR(MAXLAYER)			! IN:Adjustment in rate needed to	
									!    achieve measured NPP and TOC	
C
C ...Manure factors
C
      REAL FPROPX(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
	REAL FYMCX(MAXORGM)			! IN/OUT:Proportion of C in FYM
	REAL FYMNX(MAXORGM)			! IN/OUT:Proportion of Organic N in FYM
	REAL FYMAX(MAXORGM)			! IN/OUT:Proportion of ammonium-N in FYM
	REAL FYMWX(MAXORGM)			! IN/OUT:Amount of water added in FYM
	REAL FYMLOSS(MAXORGM)		! IN/OUT:Prop.FYM lost by volatn/timestep
	REAL FYMPOS(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
	REAL FYMSTART(MAXORGM)		! IN/OUT:Amount of rainfall in 1 week, 
								!		 below which volatilisation will occur
C
C ...Fertiliser factors
C
	REAL BYRAIN					! IN/OUT:Rain lost by bypass flow mm/timestep	
	REAL FERTADD(MAXFERTSTEPS,4)	! INOUT:Fertiliser N added in this time 
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
	REAL FERTADD15(MAXFERTSTEPS,4)	! INOUT:Fertiliser N15 added in this time
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
C
C ...Dissolved organic matter factors
C
	REAL MOBDOC(MAXLAYER)		! IN:Mobile DOC in each layer (kgC/ha)
	REAL MOBDON(MAXLAYER)		! IN:Mobile DON in each layer (kgN/ha)
C
C Get soil parameters
C
  	CALL GET_SUNDIAL_SOILPARS(NSOIL,TOC,IOM,CRIT,ILU,
     &                          CLAY,BULKDENS,SPARMODEL)
C
C Retrieve soil CN characteristics 
C
      ISAVE=0
      CALL SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,CONVER_F,
     &             SECONDS,NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)
	CALL SAVE_FERTADD(SECONDS,FERTADD,FERTADD15,ISAVE)
C
C Set time factors according to SECONDS
C
      CALL SETTIME_CN(SECONDS,CONVER_F)
C
C Set starting value of organic matter
C
      CALL STARTORG(NSOIL,TOC,TOCARRAY,IOM,IOMARRAY)
C
C Set BIO and HUM decompostion parameters
C
      CALL OMPROP(NSOIL,HY,HZ1,ALPHA,BETA,GAMMA,DELTA,
     &                  BPART,HPART,BPROP,HPROP,PHARRAY)
C
C Set NO3 and NH4 residual
C
      CALL SOILFIX(NSOIL,CRIT,AMMN,AMMN15,SOILN,SOIL15)
C
C Set the residual N and C in the RO,BIO and HUM pools
C
      IF(INMODEL.EQ.INSTABLECN)THEN
        CALL GET_CTON_BRADBURY(DPMCTON,RPMCTON,HZ1)
	ENDIF
      CALL SETORGM_FIXED(NSOIL,DPMRATE,RPMRATE,TOC,IOM,CONVER_F,
     &                   BRATE,BIORATE,HRATE,HUMRATE,HZ1,HY,
     &                   BCARB0,HCARB0,DPMCARB0,RPMCARB0,
     &                   BNIT0,HNIT0,DPMNIT0,RPMNIT0,
     &                   HNLAB0,BNLAB0,DPMNLAB0,RPMNLAB0,ICFACTOR,
     &                   DPM_RPM_RATIO,DPMCTON,RPMCTON)	
C
C Set residual DOC content
C
      IF(DOCMODEL.EQ.DOC_ON)CALL SETDOC(NSOIL,MOBDOC,MOBDON)
C
C Set initial fertiliser additions
C
      CALL SETFERT(FERTADD,FERTADD15)
C
C Pass soil pH out of routine
C
      DO 100 IL=1,MAXLAYER
	  SOILPH(IL)=PHARRAY(NSOIL,IL)
100   CONTINUE
C
C Save soil CN characteristics 
C
      ISAVE=1
      CALL SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,CONVER_F,
     &             SECONDS,NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)
	CALL SAVE_FERTADD(SECONDS,FERTADD,FERTADD15,ISAVE)
C
C Leave INIT_CN
C
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE INIT_SUNDIAL_SOILCN_NOPARS(
C INPUTS: Time, environment	
     &                               ICMODEL,INMODEL,DOCMODEL,
     &                               EQMODEL,SECONDS,
     &                               NSOIL,
C OUTPUTS: Soil C&N     
     &                               TOC,IOM,HZ1,CRIT,
C INPUTS/OUTPUTS: Soil C&N
     &                               SOILN,SOIL15,AMMN,AMMN15,
     &                               DPMCARB0,DPMNIT0,DPMNLAB0,
     &                               RPMCARB0,RPMNIT0,RPMNLAB0,
     &                               BCARB0,BNIT0,BNLAB0,
     &                               HCARB0,HNIT0,HNLAB0,
     &                               MOBDOC,MOBDON,
     &                               ILU,CLAY,BULKDENS,PI_CEQ_MON,
     &                               LTA_AWC,LTA_TEMP,ITFUNC,IMFUNC,
     &                               WMAX,WSAT,ICFACTOR,
     &                               DPM_RPM_RATIO,DPMCTON,RPMCTON,
     &                               SOILPH)
C
C Subroutine to initialise soil C and N
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLU				! Maximum number of land use types
	PARAMETER (MAXLU=6)	
	INTEGER MAXSOIL				! Max.no.of soil types
	INTEGER MAXORGM				! Max.no.of org.manure applications allowed
	PARAMETER(MAXSOIL=50,MAXORGM=52)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER IL					! Layer counter
	INTEGER MAXFERTSTEPS		! Maximum number of timesteps 
								!	over which fertiliser is added
      PARAMETER (MAXFERTSTEPS=30*60*24*7*3)	! Timesteps = 3 wks (max.30 secs)
C
C ...Model descriptors
C
	INTEGER DOC_OFF				! DOC model off
	INTEGER DOC_ON				! DOC model on
	DATA DOC_ON,DOC_OFF /1,2/
	INTEGER EQNPP				! Model initialised using plant inputs 
								!    from measured NPP
	INTEGER EQTOC				! Model initialised using measured TOC
	INTEGER EQNPPTOC			! Model initialised using both plant inputs
								!    from measured NPP and measured TOC
	INTEGER EQHILLIER,EQJONES			! Model initialised using TOC and
								!    the Hillier solver to get C pools and plant inputs
	DATA EQNPP,EQTOC,EQNPPTOC,EQHILLIER,EQJONES /1,2,3,4,5/
C
C Variables passed to/from this subroutine
C
      INTEGER ISAVE				! OUT:Code to save or retrieve variables
C
C ...Weather factors
C
	REAL LTA_AWC(12,MAXLAYER)	! IN:Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)	! IN:Average air temp this month (deg.C)
C
C ...Model descriptors
C
	INTEGER ICMODEL				! OUT:Type of C initialisation model 
	INTEGER ICFIXED				! Initialisation of C is fixed
	INTEGER ICROTHCEQ			! Initialisation of C by RothC 
								!  equilibrium run
	DATA ICFIXED,ICROTHCEQ /1,2/
	INTEGER INMODEL				! OUT:Type of N initialisation model
	INTEGER INSTABLECN			! Initisalisation of N by assuming 
								!  steady state (after Bradbury)
	INTEGER INPASSCN			! Initialisation of N by passing the
								!  C:N ratio of DPM and RPM
	DATA INSTABLECN,INPASSCN /1,2/
	INTEGER DOCMODEL			! IN:DOC model (on or off)
	INTEGER EQMODEL				! IN:Type of equilibrium run (NPP,TOC or both)
	INTEGER ITFUNC				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IMFUNC				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
C
C ...Timing factors
C
      REAL CONVER_F				! IN/OUT:Conversion to correct timestep
	REAL SECONDS				! IN:Number of seconds in one timestep
C
C ...Soil factors
C
	REAL ALPHA(MAXLAYER)			! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
	REAL ANIT					! IN/OUT:NH4-N nitrified(kgN/ha)
	REAL ANIT15					! IN/OUT:NH4-N15 nitrified (kgN/ha)
	REAL BCARB0(MAXLAYER)		! IN:C in soil biomass (kgC/ha/layer)
	REAL BETA(MAXLAYER)				! IN/OUT: Prop.HUM produced on BIO decompn
	REAL BIOP(MAXSOIL,MAXLAYER)		! IN/OUT: BIO/TOC 
	REAL BIORATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for BIO decompn/yr
      REAL BNIT0(MAXLAYER)		! IN:N in soil biomass (kgN/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN:N15 in soil humus (kgN15/ha/layer)
	REAL BPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL BPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from biomass decompositn
	REAL BRATE(MAXLAYER)			! IN/OUT: Rate constant for HUM decompn
	REAL BULKDENS(MAXLAYER)		! IN:Bulk density of the layer (g/cm3)
      REAL CLAY(MAXLAYER)			! IN:Clay content of the layer (%)
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: % clay content in this layer
	REAL CRIT(MAXSOIL,MAXLAYER)		! IN/OUT: Minimum nitrate level (Sozanska)
	REAL DELTA(MAXLAYER)			! IN/OUT: Prop.HUM produced on HUM decompn 
	REAL DFACT(MAXLAYER)			! IN/OUT: Denitrification factor
      REAL DPM_RPM_RATIO			! OUT:Ratio of DPM:RPM. If set to 
								!     zero, this will be worked out using a call to SET_DPMRPMRATIO
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
      REAL DPMCTON(MAXLAYER)		! IN/OUT:Ratio of C:N in DPM
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL DPMRATE(MAXLAYER)			! IN/OUT: Rate constant for DPM decompn
	REAL FANIT					! IN/OUT:Fertiliser NH4-N nitrified(kgN/ha)
	REAL FANIT15				! IN/OUT:Fertiliser NH4-N15 nitrified(kgN/ha)
	REAL GAMMA(MAXLAYER)			! IN/OUT: Prop.BIO produced on HUM decompn 
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN:N in soil humus (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN:N15 in soil biomass (kgN15/ha/layer)
	REAL HPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL HPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from humus decompositn
	REAL HRATE(MAXLAYER)			! IN/OUT: Rate constant for BIO decompn
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for HUM decompn/yr
	REAL HY(MAXSOIL,MAXLAYER)	! IN/OUT:Stable N:C ratio of BIO & HUM pools
	REAL HZ1(MAXLAYER)				! IN/OUT: N:C ratio for steady state
  	REAL ICFACTOR(MAXLAYER)			! IN:Adjustment in rate needed to	
									!    achieve measured NPP and TOC	
	INTEGER ILU						! IN:Counter for land use 
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
	REAL IOMARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium IOM (kgC/ha/layer)
      INTEGER LUARRAY(MAXSOIL)	! IN/OUT:Land use before equilibrium 
	INTEGER NSOIL				! IN:Soil code number
      INTEGER NUMSOIL				! IN/OUT:Number of soils defined
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)	! IN/PUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)	! IN/PUT: pH above which decomp.max.
	REAL PI_CEQ_MON(12,MAXLAYER)	! IN/OUT: Equilibrium plant C input each 
									!         month to this layer (kgC/ha/month/layer)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
      REAL RPMCTON(MAXLAYER)		! IN/OUT:Ratio of C:N in RPM
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL RPMRATE(MAXLAYER)			! IN/OUT: Rate constant for RPM decompn
	CHARACTER*40 SNAME(MAXSOIL)	! IN/OUT:Soil name
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILPH(MAXLAYER)		! IN/OUT:pH of soil in this layer
	REAL TOC(MAXLAYER)			! IN:Total organic C (kgC/ha/layer)
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium TOC (kgC/ha/layer)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
C
C ...Manure factors
C
      REAL FPROPX(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
	REAL FYMCX(MAXORGM)			! IN/OUT:Proportion of C in FYM
	REAL FYMNX(MAXORGM)			! IN/OUT:Proportion of Organic N in FYM
	REAL FYMAX(MAXORGM)			! IN/OUT:Proportion of ammonium-N in FYM
	REAL FYMWX(MAXORGM)			! IN/OUT:Amount of water added in FYM
	REAL FYMLOSS(MAXORGM)		! IN/OUT:Prop.FYM lost by volatn/timestep
	REAL FYMPOS(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
	REAL FYMSTART(MAXORGM)		! IN/OUT:Amount of rainfall in 1 week, 
								!		 below which volatilisation will occur
C
C ...Fertiliser factors
C
	REAL BYRAIN					! IN/OUT:Rain lost by bypass flow mm/timestep	
	REAL FERTADD(MAXFERTSTEPS,4)	! INOUT:Fertiliser N added in this time 
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
	REAL FERTADD15(MAXFERTSTEPS,4)	! INOUT:Fertiliser N15 added in this time
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
C
C ...Dissolved organic matter factors
C
	REAL MOBDOC(MAXLAYER)		! IN:Mobile DOC in each layer (kgC/ha)
	REAL MOBDON(MAXLAYER)		! IN:Mobile DON in each layer (kgN/ha)
C
C Retrieve soil CN characteristics 
C
      ISAVE=0
      CALL SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,CONVER_F,
     &             SECONDS,NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)
	CALL SAVE_FERTADD(SECONDS,FERTADD,FERTADD15,ISAVE)
C
C Set time factors according to SECONDS
C
      CALL SETTIME_CN(SECONDS,CONVER_F)
C
C Set starting value of organic matter
C
      CALL STARTORG(NSOIL,TOC,TOCARRAY,IOM,IOMARRAY)
C
C Set BIO and HUM decompostion parameters
C
      CALL OMPROP(NSOIL,HY,HZ1,ALPHA,BETA,GAMMA,DELTA,
     &                  BPART,HPART,BPROP,HPROP,PHARRAY)
C
C Set NO3 and NH4 residual
C
      CALL SOILFIX(NSOIL,CRIT,AMMN,AMMN15,SOILN,SOIL15)
C
C Set the residual N and C in the RO,BIO and HUM pools
C
      IF(INMODEL.EQ.INSTABLECN)THEN
        CALL GET_CTON_BRADBURY(DPMCTON,RPMCTON,HZ1)
	ENDIF
	IF(ICMODEL.EQ.ICROTHCEQ)THEN
        CALL SETORGM_EQRUN(NSOIL,DPMRATE,RPMRATE,TOC,IOM,CONVER_F,
     &                   BRATE,BIORATE,HRATE,HUMRATE,HZ1,HY,
     &                   BCARB0,HCARB0,DPMCARB0,RPMCARB0,
     &                   BNIT0,HNIT0,DPMNIT0,RPMNIT0,
     &                   HNLAB0,BNLAB0,DPMNLAB0,RPMNLAB0,
     &                   LUARRAY,CLARRAY,
     &                   PHARRAY,PHP1ARRAY,PHP2ARRAY,
     &                   EQMODEL,PI_CEQ_MON,ICFACTOR,LTA_AWC,LTA_TEMP,
     &                   ITFUNC,IMFUNC,WMAX,WSAT,
     &                   DPM_RPM_RATIO,DPMCTON,RPMCTON,ALPHA,BETA,CLAY)	
      ELSEIF(ICMODEL.EQ.ICFIXED)THEN
        CALL SETORGM_FIXED(NSOIL,DPMRATE,RPMRATE,TOC,IOM,CONVER_F,
     &                   BRATE,BIORATE,HRATE,HUMRATE,HZ1,HY,
     &                   BCARB0,HCARB0,DPMCARB0,RPMCARB0,
     &                   BNIT0,HNIT0,DPMNIT0,RPMNIT0,
     &                   HNLAB0,BNLAB0,DPMNLAB0,RPMNLAB0,ICFACTOR,
     &                   DPM_RPM_RATIO,DPMCTON,RPMCTON)	
	ENDIF 
C
C Set residual DOC content
C
      IF(DOCMODEL.EQ.DOC_ON)CALL SETDOC(NSOIL,MOBDOC,MOBDON)
C
C Set initial fertiliser additions
C
      CALL SETFERT(FERTADD,FERTADD15)
C
C Pass soil pH out of routine
C
      DO 100 IL=1,MAXLAYER
	  SOILPH(IL)=PHARRAY(NSOIL,IL)
100   CONTINUE
C
C Save soil CN characteristics 
C
      ISAVE=1
      CALL SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,CONVER_F,
     &             SECONDS,NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)
	CALL SAVE_FERTADD(SECONDS,FERTADD,FERTADD15,ISAVE)
C
C Leave INIT_CN
C
	END
C
C-------------------------------------------------------------
C
	SUBROUTINE MIAMI(NPP,AVERAIN,AVETEMP)
C
C Subroutine to calculate NPP
C Equations for NPP based on annual rainfall and temp 
C Leith, H., 1972. Modelling the primary productivity of the world. 
C Nature and Resources, UNESCO, VIII, 2:5-10.
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      REAL ANNRAIN				! Annual rainfall (mm)
	REAL ANNTEMP				! Annual temp (deg.C)
	INTEGER IMON				! Month counter
	REAL NPPTEMP				! Temp limited NPP (g DM / m2 / year)
	REAL NPPRAIN				! Rain limited NPP (g DM / m2 / year)
	REAL CTODM					! Assumed C:DM content of plant material
	DATA CTODM /0.4/
C
C Variables passed to/from calling subroutine
C ...Plant factors
C
	REAL NPP					! IN/OUT:Plant C input calculated using 
								!	  ROTHC equilibrium run for each LU 
								!     (kgC/ha/yr)
C
C ...Met factors
C
      REAL AVERAIN(12)			! Long term average rainfall(mm/month)
      REAL AVETEMP(12)			! Long term average monthly average 
								!	air temp (deg.C)
C
C Get annual rainfall and air temp
C
      ANNRAIN=0
	ANNTEMP=0
      DO 100 IMON=1,12
        ANNRAIN=ANNRAIN+AVERAIN(IMON)
	  ANNTEMP=ANNTEMP+AVETEMP(IMON)
100   CONTINUE
      ANNTEMP=ANNTEMP/12
C
C Calculate NPP
C
	NPPRAIN=3000*(1-EXP(-0.000664*ANNRAIN))
	NPPTEMP=3000*1/(1+EXP(1.315-(0.119*ANNTEMP)))
      IF(NPPRAIN.LT.NPPTEMP)THEN
	  NPP=CTODM*10*NPPRAIN
	ELSE
	  NPP=CTODM*10*NPPTEMP
	ENDIF
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE MICROBIAL_SUNDIAL_SOILCN(
C INPUTS: time, environment, organic manure, fertilisers, crop debris, soil water
     &                                 DMODEL,DOCMODEL,CNMODEL,                                 
     &                                 ITFUNC,IMFUNC,CH4MODEL,INMODEL,
     &                                 SECONDS,ICOVER,
     &                                 NSOIL,SOILTEMP,RAIN,
     &								 PI_C,PI_N,PI_N15,DRRAT,                                                                     
     &                                 WMAX,SOILW,WSAT,
     &                                 WILTPOINT,SATWATCONT,
C OUTPUTS: fertilisers, atmospheric, mineralisation, leaching
     &                                 THISFERT,THISFERT15,
     &                                 CO2,VOLAT,VOLAT15,VOLATN,
     &                                 CH4,CH4_PROD,CH4_SOIL_OX,
     &							     CH4_ATMOS_OX, CH4_FLUX,
     &                                 DENITRIFN,NITRIFN,
     &                                 DENIT,DN15,GNN2O,GNNO,GPNN2O,	
     &                                 G15NN2O,G15NNO,G15PNN2O,		
     &                                 GDN2,GDN2O,G15DN2,G15DN2O,		
     &                                 DNIT,DNIT15,TDNIT,T15,
     &                                 WLEACH,WLEACH15,
C INPUTS/OUTPUTS: Soil C&N
     &                                 SOILN,SOIL15,AMMN,AMMN15,
     &                                 DPMCARB0,DPMNIT0,DPMNLAB0,
     &                                 RPMCARB0,RPMNIT0,RPMNLAB0,
     &                                 BCARB0,BNIT0,BNLAB0,BULKDENS,
     &                                 TOC,CH4TOAIR,HCARB0,HNIT0,HNLAB0,
     &                                 TIN,TAM,MOBDOC,MOBDON,CO2FROMDOC,
     &                                 ICFACTOR,DPMCTON,RPMCTON,
C Required for Monte Carlo - 
     &					DPMCARB,RPMCARB,BCARB,HCARB,
     &					WRATEDM,TRATEM,SOILPH,MEASLAY)
C
C Subroutine to run SUNDIAL soil C and N
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXFERTSTEPS		! Maximum number of timesteps 
								!	over which fertiliser is added
      PARAMETER (MAXFERTSTEPS=30*60*24*7*3)	! Timesteps = 3 wks (max.30 secs)
 	INTEGER MAXORGM				! Max.no.of organic manure applications allowed
	PARAMETER (MAXORGM=52)
	INTEGER MAXSOIL				! Max.no.of soil types
	PARAMETER (MAXSOIL=50)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/

	INTEGER CH4_OFF				! CH4 model off
	INTEGER CH4_RICHARDS	    ! Richards CH4 model on		
	INTEGER CH4_AITKENHEAD      ! Aitkenhead CH4 model on
	DATA CH4_OFF,CH4_RICHARDS,CH4_AITKENHEAD /0,1,2/
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: % clay content in this layer
	INTEGER CNFOEREID			! C:N ratio obtained by method of Foereid
	INTEGER CNMAGEC				! C:N ratio obtained by method of MAGEC
	DATA CNMAGEC,CNFOEREID /1,2/
	REAL CO2DPM(MAXLAYER)       ! CO2 emitted from the DPM pool [kgC/ha/layer]
	REAL CO2RPM(MAXLAYER)       ! CO2 emitted from the RPM pool [kgC/ha/layer]
	REAL CO2BIO(MAXLAYER)       ! CO2 emitted from the BIO pool [kgC/ha/layer]
	REAL CO2HUM(MAXLAYER)       ! CO2 emitted from the HUM pool [kgC/ha/layer]
	INTEGER DBRADBURY			! Bradbury model for denitrification
	INTEGER DNEMIS				! Nemis model for denitrification
	DATA DBRADBURY,DNEMIS /1,2/
	INTEGER DOC_OFF				! DOC model off
	INTEGER DOC_ON				! DOC model on
	DATA DOC_ON,DOC_OFF /1,2/
      INTEGER IL					! Local layer counter
	INTEGER INSTABLECN			! Initisalisation of N by assuming 
								!  steady state (after Bradbury)
	INTEGER INPASSCN			! Initialisation of N by passing the
								!  C:N ratio of DPM and RPM
	DATA INSTABLECN,INPASSCN /1,2/
      INTEGER ISAVE				! OUT:Code to save or retrieve variables

      INTEGER LUARRAY(MAXSOIL)	! IN/OUT:Land use before equilibrium 
      INTEGER NUMSOIL				! IN/OUT:Number of soils defined
C
C Variables passed to/from this subroutine
C
	REAL ALPHA(MAXLAYER)			! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL AMMN(MAXLAYER)				! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)			! IN:Soil ammonium-N15 (kgN15/ha/layer)
	REAL ANIT						! IN/OUT:NH4-N nitrified(kgN/ha)
	REAL ANIT15						! IN/OUT:NH4-N15 nitrified (kgN/ha)
	REAL BCARB(MAXLAYER)			! IN:C in BIO at end of timestep (MC)
  	REAL BCARB0(MAXLAYER)			! IN:C in soil biomass (kgC/ha/layer)
	REAL BETA(MAXLAYER)				! IN/OUT: Prop.HUM produced on BIO decompn
	REAL BIOP(MAXSOIL,MAXLAYER)		! IN/OUT: BIO/TOC 
	REAL BIORATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for BIO decompn/yr
      REAL BNIT0(MAXLAYER)			! IN:N in soil biomass (kgN/ha/layer)
	REAL BNLAB0(MAXLAYER)			! IN:N15 in soil humus (kgN15/ha/layer)
	REAL BPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL BPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from biomass decompositn
	REAL BRATE(MAXLAYER)			! IN/OUT: Rate constant for HUM decompn
	REAL BULKDENS(MAXLAYER)			! OUT: Bulk density of soil (g cm-3)
	REAL BYRAIN						! IN/OUT:Rain lost by bypass flow mm/timestep	
	REAL CH4(MAXLAYER)				! IN:CH4 emitted (kgC/ha/layer)
	INTEGER CH4MODEL				! IN:Methane model (off, Richards or Aitkenhead model on) 
	REAL CH4TOAIR					! OUT:CH4 release to atmosphere (kgC/ha)
	REAL CH4_PROD(MAXLAYER)  		! OUT: CH4 produced in each soil layer [kgC/ha/layer] (Richards CH4 model)
	REAL CH4_SOIL_OX				! OUT: CH4 oxidised from methane produced in the soil [kgC/ha] (Richards CH4 model)
      REAL CH4_ATMOS_OX				! OUT: CH4 oxidised from atmospheric CH4 [kgC/ha] (Richards CH4 model)
      REAL CH4_FLUX					! OUT: Net CH4 flux from soil column [kgC/ha] (Richards CH4 model)
	INTEGER CNMODEL					! OUT:Type of Model for C:N ratio of DPM and 
									!     RPM (MAGEC or Foereid)
	REAL CO2(MAXLAYER)				! IN:CO2 emitted (kgC/ha/layer)
	REAL CO2FROMDOC(MAXLAYER)		! IN:CO2 production rate (kgC/ha/layer/day)
      REAL CONVER_F					! IN/OUT:Conversion to correct timestep
	REAL CRIT(MAXSOIL,MAXLAYER)		! IN/OUT: Minimum nitrate level (Sozanska)
	REAL DELTA(MAXLAYER)			! IN/OUT: Prop.HUM produced on HUM decompn 
	REAL DENIT						! IN:N lost by denitrification (kgN/ha)
	REAL DENITRIFN(MAXLAYER)		! IN/OUT: N lost by denitrification (kgN/ha/layer
	REAL DFACT(MAXLAYER)			! IN/OUT: Denitrification factor
	INTEGER DMODEL					! OUT:Type denitrification model chosen 
									!     (Bradbury or NEMIS)
	REAL DN15						! IN:N15 lost by denitrification (kgN15/ha)
	REAL DNIT(MAXLAYER)				! IN:Net Mineralised N (kgN/ha/lay/timestep)
	REAL DNIT15(MAXLAYER)			! IN:Net Mineralised N15 (kgN15/ha/lay/time)
	INTEGER DOCMODEL				! IN:DOC model (on or off)
     	REAL DPMCARB(MAXLAYER)			! IN:C in DPM at end of timestep (MC) 
	REAL DPMCARB0(MAXLAYER)			! IN:C in decomposable PM (kgC/ha/layer)
	REAL DPMNIT0(MAXLAYER)			! IN:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)			! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL DPMRATE(MAXLAYER)			! IN/OUT: Rate constant for DPM decompn
	REAL DRRAT						! IN:DPM:RPM ratio
      REAL DPMCTON(MAXLAYER)			! IN:Ratio of C:N in DPM
	REAL FANIT						! IN/OUT:Fertiliser NH4-N nitrified(kgN/ha)
	REAL FANIT15					! IN/OUT:Fertiliser NH4-N15 nitrified(kgN/ha)
	REAL FERTADD(MAXFERTSTEPS,4)	! INOUT:Fertiliser N added in this time 
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
	REAL FERTADD15(MAXFERTSTEPS,4)	! INOUT:Fertiliser N15 added in this time
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
      REAL FPROPX(MAXORGM)			! OUT:Proportion of FYM added to top 25cm
	REAL FYMCX(MAXORGM)				! OUT:Proportion of C in FYM
	REAL FYMNX(MAXORGM)				! OUT:Proportion of Organic N in FYM
	REAL FYMAX(MAXORGM)				! OUT:Proportion of ammonium-N in FYM
	REAL FYMWX(MAXORGM)				! OUT:Amount of water added in FYM
	REAL FYMLOSS(MAXORGM)			! OUT:Prop.FYM lost by volatn/timestep
	REAL FYMPOS(MAXORGM)			! OUT:Proportion of FYM added to top 25cm
	REAL FYMSTART(MAXORGM)			! OUT:Amount of rainfall in 1 week, 								! 	 below which volatilisation will occur
	REAL G15DN2						! IN:15N2 lost by denitrification (kgN15/ha)
	REAL G15DN2O					! IN:15N2O lost by denitrification (kgN15/ha)					
	REAL G15NN2O					! IN:15N2O lost by nitrification (kgN15/ha)
	REAL G15NNO						! IN:15NO lost by nitrification (kgN15/ha)			
	REAL G15PNN2O					! IN:15N2O lost by part.nitrif.(kgN15/ha)			
	REAL GAMMA(MAXLAYER)			! IN/OUT: Prop.BIO produced on HUM decompn 
	REAL GDN2						! IN:N2 lost by denitrification (kgN/ha)
	REAL GDN2O						! IN:N2O lost by denitrification (kgN/ha)
	REAL GNN2O						! IN:N2O lost by nitrification (kgN/ha)
	REAL GNNO						! IN:NO lost by nitrification (kgN/ha)
	REAL GPNN2O						! IN:N2O lost by part.nitrification (kgN/ha)
	REAL HCARB(MAXLAYER)			! IN:C in HUM at end of timestep (MC)
	REAL HCARB0(MAXLAYER)			! IN:C in soil humus (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)			! IN:N in soil humus (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)			! IN:N15 in soil biomass (kgN15/ha/layer)
	REAL HPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL HPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from humus decompositn
	REAL HRATE(MAXLAYER)			! IN/OUT: Rate constant for BIO decompn
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for HUM decompn/yr
	REAL HY(MAXSOIL,MAXLAYER)		! IN/OUT:Stable N:C ratio of BIO & HUM pools
	REAL HZ1(MAXLAYER)				! IN/OUT: N:C ratio for steady state
	REAL ICFACTOR(MAXLAYER)			! IN:Adjustment in rate needed to 
									!    achieve measured NPP and TOC
      INTEGER ICOVER					! IN/OUT:Crop cover 1=Covered 0=Bare
	INTEGER ITFUNC					! IN:Choice of temperature rate modifier 
									! (ROTHC or HADLEY)
	INTEGER IMFUNC					! IN:Choice of moisture rate modifier 
									! (ROTHC or HADLEY)
	INTEGER INMODEL					! OUT:Type of N initialisation model
	REAL IOMARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium IOM (kgC/ha/layer)
	INTEGER MEASLAY					! Layer that soil is measured to
	REAL MOBDOC(MAXLAYER)			! IN:Total mobile DOC pool (kgC/ha/layer)
	REAL MOBDON(MAXLAYER)			! IN:Total mobile DON pool (kgN/ha/layer)
	REAL NITRIFN(MAXLAYER)			! IN/OUT:Nitrification from this layer (kg/ha/layer/timestep)
	INTEGER NSOIL					! IN:Soil code number
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)	! IN/PUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)	! IN/PUT: pH above which decomp.max.
      REAL PI_C(MAXLAYER)				! IN/OUT: Plant input C to soil (kgC/ha/step)
      REAL PI_N(MAXLAYER)				! IN/OUT: Plant input N to soil (kgC/ha/step)
      REAL PI_N15(MAXLAYER)			! IN/OUT: Plant input N15 to soil (kgC/ha/step)
	REAL RAIN					! IN: Rainfall (mm/timestep)
	REAL RPMCARB(MAXLAYER)			! IN:C in RPM at end of timestep (MC)
	REAL RPMCARB0(MAXLAYER)			! IN:C in resistant PM (kgC/ha/layer)
      REAL RPMCTON(MAXLAYER)			! IN:Ratio of C:N in RPM
	REAL RPMNIT0(MAXLAYER)			! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)			! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL RPMRATE(MAXLAYER)			! IN/OUT: Rate constant for RPM decompn
	REAL SATWATCONT(MAXLAYER)	    ! IN:Total water content at saturation (mm/layer)
	REAL SECONDS					! IN:Number of seconds in one timestep
	CHARACTER*40 SNAME(MAXSOIL)		! IN/OUT:Soil name
	REAL SOIL15(MAXLAYER)			! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILN(MAXLAYER)			! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOILPH(MAXLAYER)			! IN/OUT:pH of soil in this layer
	REAL SOILTEMP(MAXLAYER)			! IN: Soil temperature (deg.C/timestep)
	REAL SOILW(MAXLAYER)			! IN:Available water (mm/layer)
	REAL T15						! IN:Net Mineralised N (kgN/ha/layer)
	REAL TAM(MAXLAYER)				! IN:15N/N in the ammonium in the layer
	REAL TDNIT						! IN:Net Mineralised N (kgN/ha/layer)
	REAL THISFERT					! IN:N input by fertiliser (kgN/ha)
	REAL THISFERT15					! IN:N15 input by fertiliser (kgN15/ha)
	REAL TIN(MAXLAYER)				! IN:15N/N in the nitrate in the layer
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium TOC (kgC/ha/layer)
	REAL TOC(MAXLAYER)						! TOC (kgC/ha/layer)
	REAL TRATEM						! IN:Temperature rate modifier (MC)
	REAL VOLAT						! IN:N lost by volatilisation (kgN/ha)
	REAL VOLAT15					! IN:N15 lost by volatilisation (kgN15/ha)
	REAL VOLATN(MAXLAYER)		! IN/OUT:Volatilisation from this layer (kg/ha/layer/timestep)
	REAL WILTPOINT(MAXLAYER)	! IN:Water content at wilting point (mm/layer)
	REAL WLEACH					! IN:N leached (kgN/ha)
	REAL WLEACH15				! IN:N15 leached (kgN15/ha)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WRATEDM				! IN:Moisture rate modifier (MC)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)

C
C Denitrification, DOC and CN model types
C
	DBRADBURY=1
	DNEMIS=2
	DOC_ON=1
	DOC_OFF=2
	CNMAGEC=1
	CNFOEREID=2
C
C Retrieve soil CN characteristics 
C
      ISAVE=0
      CALL SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,CONVER_F,
     &             SECONDS,NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)
	CALL SAVE_FERTADD(SECONDS,FERTADD,FERTADD15,ISAVE)
C
C Get CTON ratios
C
      IF(INMODEL.EQ.INSTABLECN)THEN
        CALL GET_CTON_BRADBURY(DPMCTON,RPMCTON,HZ1)
	ENDIF
C
C Pass soil pH into routine
C
      DO 100 IL=1,MAXLAYER
	  PHARRAY(NSOIL,IL)=SOILPH(IL)
100   CONTINUE
C
C Calculate mineralisation
C TDNIT = Mineralisation for this timestep
C
      IF(CNMODEL.EQ.CNFOEREID)THEN
        CALL MINER1_FOEREID(NSOIL,ICOVER,ITFUNC,IMFUNC,
     &                    PI_C,PI_N,PI_N15,DRRAT,
     &                    SOILW,WMAX,WSAT,SOILTEMP,
     &                    ALPHA,BETA,GAMMA,DELTA,HY,
     &                    BRATE,HRATE,DPMRATE,RPMRATE,
     &                    AMMN,AMMN15,SOILN,SOIL15,CRIT,
     &                    CO2,
     &                    DNIT,DNIT15,TDNIT,T15,
     &                    BCARB0,BNIT0,BNLAB0,
     &                    HCARB0,HNIT0,HNLAB0,
     &                    DPMCARB0,DPMNIT0,DPMNLAB0,
     &                    RPMCARB0,RPMNIT0,RPMNLAB0,
     &					PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &                    DPMCTON,RPMCTON,
C Required for Monte Carlo - 
     &					DPMCARB,RPMCARB,BCARB,HCARB,
     &					WRATEDM,TRATEM)
      ELSE
c        CALL MINER1_BRADBURY(NSOIL,ICOVER,ITFUNC,IMFUNC,
c     &                    PI_C,PI_N,PI_N15,DRRAT,
c     &                    SOILW,WMAX,WSAT,SOILTEMP,
c     &                    ALPHA,BETA,GAMMA,DELTA,HY,
c     &                    BRATE,HRATE,DPMRATE,RPMRATE,
c     &                    AMMN,AMMN15,SOILN,SOIL15,CRIT,
c     &                    CO2, CO2DPM, CO2RPM, CO2BIO, CO2HUM,
c     &                    DNIT,DNIT15,TDNIT,T15,
c     &                    BCARB0,BNIT0,BNLAB0,
c     &                    HCARB0,HNIT0,HNLAB0,
c     &                    DPMCARB0,DPMNIT0,DPMNLAB0,
c     &                    RPMCARB0,RPMNIT0,RPMNLAB0,
c     &					PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
c     &                    DPMCTON,RPMCTON,
C Required for Monte Carlo - 
c     &					DPMCARB,RPMCARB,BCARB,HCARB,
c     &					WRATEDM,TRATEM,MEASLAY)
        CALL ADD_PI_TO_DPMRPM(MAXLAYER, MEASLAY, DRRAT, PI_C, PI_N,
     &                 ICOVER,DPMCARB0, DPMNIT0, RPMCARB0, RPMNIT0)
        
        CALL MINER1_RICHARDS(MAXLAYER, MAXDEPTH, MEASLAY,
     &                 NSOIL, MAXSOIL, ICOVER,
     &                 DPMRATE, RPMRATE, BRATE, HRATE,
     &                 ALPHA, BETA, DELTA, GAMMA, HY,
     &                 ITFUNC, IMFUNC, ICFACTOR,
     &                 SOILW, WMAX, WSAT, SOILTEMP, CRIT,
     &                 WILTPOINT, SATWATCONT,
     &                 PHARRAY, PHP1ARRAY, PHP2ARRAY,
     &                 SOILN, AMMN,
     &                 DPMCARB0 , DPMNIT0, RPMCARB0, RPMNIT0,
     &                 BCARB0, BNIT0, HCARB0, HNIT0,
     &                 CO2, CO2DPM, CO2RPM, CO2BIO, CO2HUM,
     &                 DNIT, TDNIT, WRATEDM, TRATEM,PI_C)
	ENDIF
C
C Calculate methane emissions
C
	IF (CH4MODEL.EQ.CH4_RICHARDS) THEN 
	    CALL METHANE_RICHARDS(MAXLAYER, MAXDEPTH, MEASLAY, MAXSOIL,
     &						  CH4_PROD, CH4_SOIL_OX, CH4_ATMOS_OX,
     &						  CH4_FLUX, SOILTEMP, SOILW, WILTPOINT,
     &                          WMAX, WSAT, SATWATCONT, BULKDENS,
     &                          CO2, CO2DPM, CO2RPM,
     &                          CO2HUM, CO2BIO, DPMCARB0, RPMCARB0,
     &                          BCARB0, HCARB0, DPMNIT0, RPMNIT0, BNIT0,
     &                          HNIT0, AMMN, SOILN, NSOIL, HY, CRIT,
     &                          ALPHA, BETA, DELTA, GAMMA, PHARRAY,
     &                          SECONDS)
      ELSEIF (CH4MODEL.EQ.CH4_AITKENHEAD) THEN 
C          CALL EMIT_METHANE(SOILTEMP,NSOIL,
C     &                    BCARB0,HCARB0,DPMCARB0,RPMCARB0,
C     &                    BNIT0,HNIT0,DPMNIT0,RPMNIT0,
C     &                    BNLAB0,HNLAB0,DPMNLAB0,RPMNLAB0,
C     &                    DPMRATE,RPMRATE,BRATE,HRATE,
C     &                    HY,SOILW,WMAX,WSAT,
C     &                    ALPHA,BETA,GAMMA,DELTA,
C     &                    SOILN,SOIL15,AMMN,AMMN15,CRIT,
C     &                    PHARRAY,CO2,CH4,CLARRAY,BULKDENS,TOC,
C     &                    SECONDS,CH4TOAIR)
	ENDIF
C
C Calculate fertilizer additions after bypass flow
C
      CALL BYPASSFLOW(
C ...Weather factors			
     &                      RAIN,BYRAIN,
C ...Timing factors
     &					  SECONDS,CONVER_F,
C ...Fertiliser factors
     &					  FERTADD,FERTADD15,THISFERT,THISFERT15,
C ...Soil factors
     &					  SOILN,SOIL15,AMMN,AMMN15,
     & 					  TAM,TIN,WLEACH,WLEACH15)
C
C Calculate the proportion of N15
C
      CALL PROPN15(SOILN,SOIL15,TIN,AMMN,AMMN15,TAM)
C
C Calculate nitrification
C
      IF(DMODEL.EQ.DNEMIS)THEN
        CALL NITRIF_NEMIS(SOILTEMP,RAIN,NSOIL,ITFUNC,IMFUNC,
     &                  SOILW,WMAX,WSAT,
     &                  SOILN,AMMN,SOIL15,AMMN15,CRIT,TAM,
     &                  FERTADD,FERTADD15,THISFERT,THISFERT15,
     &                  FANIT,FANIT15,VOLAT,VOLAT15,CONVER_F,
     &                  PHARRAY,PHP1ARRAY,PHP2ARRAY,
     &                  GNN2O,GNNO,GPNN2O,G15NN2O,G15NNO,G15PNN2O,
     &                  NITRIFN,VOLATN) 
	ELSE
        CALL NITRIF_BRADBURY(SOILTEMP,RAIN,NSOIL,ITFUNC,IMFUNC,
     &                  SOILW,WMAX,WSAT,
     &                  SOILN,AMMN,SOIL15,AMMN15,CRIT,TAM,
     &                  FERTADD,FERTADD15,THISFERT,THISFERT15,
     &                  FANIT,FANIT15,VOLAT,VOLAT15,CONVER_F,
     &                  PHARRAY,PHP1ARRAY,PHP2ARRAY,NITRIFN,VOLATN)
	ENDIF
C
C Reset the proportion N15:N
C
      CALL PROPN15(SOILN,SOIL15,TIN,AMMN,AMMN15,TAM)
C
C Calculate denitrification
C
      IF(DMODEL.EQ.DNEMIS)THEN
        CALL DENITRIF3_NEMIS(DFACT,DENIT,DN15,RAIN,CH4,
     &                     WMAX,CO2,SOILW,SOILN,CRIT,NSOIL,TIN,SOIL15,
     &                     GDN2,GDN2O,G15DN2,G15DN2O,DENITRIFN,BULKDENS)
	ELSE
        CALL DENITRIF3_BRADBURY(DFACT,DENIT,DN15,RAIN,CH4,CO2,
     &                          WMAX,SOILW,SOILN,CRIT,NSOIL,TIN,SOIL15,
     &                          DENITRIFN)
	ENDIF
C
C Calculate exchange of C and N into dissolved organic C and N
C
      IF(DOCMODEL.EQ.DOC_ON)CALL MAKE_DOC(PHP1ARRAY,PHP2ARRAY,PHARRAY,
     &                    NSOIL,ICOVER,ITFUNC,IMFUNC,
     &                    WMAX,WSAT,SOILW,SOILTEMP,
     &                    DPMCARB0,RPMCARB0,BCARB0,HCARB0,
     &                    DPMNIT0,RPMNIT0,BNIT0,HNIT0,
     &                    DPMNLAB0,RPMNLAB0,BNLAB0,HNLAB0,
     &                    AMMN,AMMN15,
     &                    MOBDOC,MOBDON,CO2FROMDOC,HY)
C
C Save soil CN characteristics 
C
      ISAVE=1
      CALL SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,CONVER_F,
     &             SECONDS,NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)
	CALL SAVE_FERTADD(SECONDS,FERTADD,FERTADD15,ISAVE)
C
C Leave RUN1_SUNDIAL_SOILCN
C
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE PARTITION_PLANTINPUT(C_TS,RNIN,RNIN15,PI_C,PI_N,PI_N15)
C
C Subroutine to partition plant input down the soil profile
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)	
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/
	DATA MAXLAYER1 /60/
      INTEGER IDEPTH				! Depth in current layer
	REAL SPL					! Split of input into the layer
	INTEGER IL					! Layer counter 
C
C Variables passed to/from this subroutine
C
	REAL C_TS					! IN:Litter C input in this timestep (kgC/ha)
	REAL RNIN					! IN:Litter N input in this timestep (kgN/ha)
	REAL RNIN15					! IN:Litter N15 input in timestep (kgN15/ha)
      REAL PI_C(MAXLAYER)			! OUT: Plant input C to soil (kgC/ha/step)
      REAL PI_N(MAXLAYER)			! OUT: Plant input N to soil (kgC/ha/step)
      REAL PI_N15(MAXLAYER)		! OUT: Plant input N15 to soil (kgC/ha/step)
C
C 80%:20% partitioning in top 50cm as used in earlier versions of SUNDIAL (Bradbury et al, 1993)
C
      DO 100 IL=1,MAXLAYER1
	  IDEPTH=IL*MAXDEPTH/(MAXLAYER1*1.)
	  IF(IDEPTH.LE.25)THEN
	    SPL=0.8*(MAXLAYER1*1.)/(MAXDEPTH)
	  ELSEIF(IDEPTH.GT.25.AND.IDEPTH.LE.50)THEN 
	    SPL=0.2*(MAXLAYER1*1.)/MAXDEPTH
	  ELSE
	    SPL=0
	  ENDIF
        PI_C(IL)=C_TS*SPL
        PI_N(IL)=RNIN*SPL
        PI_N15(IL)=RNIN15*SPL
100   CONTINUE
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE PHYSICAL_SUNDIAL_SOILCN(
C INPUTS: time, environment	
     &                                   DOCMODEL,SECONDS,
     &                                   NSOIL,SOILTEMP,
     &                                   WMAX,SOILW,DRAINW,REFIL,
C OUTPUTS: Soil C&N
     &                                   SLEACH,SLEA15,
     &                                   WLEACH,WLEACH15,CONC,CONC15,
C INPUTS/OUTPUTS: Soil C&N
     &                                   SOILN,SOIL15,AMMN,AMMN15,
     &                                   TIN,TAM,MOBDOC,MOBDON,
     &                                   LEACHDOC,LEACHDON)
C
C Subroutine to run SUNDIAL soil C and N
C
      IMPLICIT NONE
C
C Variable local to this subroutine
C
	INTEGER MAXSOIL				! Max.no.of soil types
	INTEGER MAXORGM				! Max.no.of organic manure applications allowed
  	INTEGER MAXFERT				! Max.no.of fertiliser applications allowed 
	PARAMETER(MAXSOIL=50,MAXORGM=52,MAXFERT=5)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/
	INTEGER IL					! Local layer counter
      INTEGER NUMSOIL				! Number of soils defined
      INTEGER ISAVE				! Code to save or retrieve variables
	REAL ANIT					! Ammonium immobilised (kgN/ha)
	REAL ANIT15					! Ammonium N15 nitrified (kgN/ha)
	REAL FANIT					! Fertiliser N nitrified (kgN/ha)
	REAL FANIT15				! Fertiliser N15 nitrified (kgN/ha)
C
C Variables passed to/from calling subroutine
C ...Timing factors
C
	REAL SECONDS				! IN:Number of seconds in one timestep
      REAL CONVER_F				! IN/OUT:Conversion to correct timestep
C
C ...Manure factors
C
      REAL FPROPX(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
	REAL FYMCX(MAXORGM)			! IN/OUT:Proportion of C in FYM
	REAL FYMNX(MAXORGM)			! IN/OUT:Proportion of Organic N in FYM
	REAL FYMAX(MAXORGM)			! IN/OUT:Proportion of ammonium-N in FYM
	REAL FYMWX(MAXORGM)			! IN/OUT:Amount of water added in FYM
	REAL FYMLOSS(MAXORGM)		! IN/OUT:Prop.FYM lost by volatn/timestep
	REAL FYMPOS(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
	REAL FYMSTART(MAXORGM)		! IN/OUT:Amount of rainfall in 1 week, 
								!		 below which volatilisation will occur
C
C ...Weather factors
C
	REAL SOILTEMP(MAXLAYER)		! IN: Soil temperature (deg.C/timestep)
C
C ...Soil factors
C
	REAL ALPHA(MAXLAYER)			! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
	REAL BETA(MAXLAYER)				! IN/OUT: Prop.HUM produced on BIO decompn
	REAL BIOP(MAXSOIL,MAXLAYER)		! IN/OUT: BIO/TOC 
	REAL BIORATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for BIO decompn/yr
	REAL BPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL BPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from biomass decompositn
	REAL BRATE(MAXLAYER)			! IN/OUT: Rate constant for HUM decompn
	REAL BYRAIN					! IN/OUT:Rain lost by bypass flow mm/timestep	
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: % clay content in this layer
      REAL CONC					! Concentration of N leached from the profile
								!		kgN/dm3
	REAL CONC15					! Concentration of N15 leached from profile
								!		kgN15/dm3
	REAL CRIT(MAXSOIL,MAXLAYER)		! IN/OUT: Minimum nitrate level (Sozanska)
	REAL DELTA(MAXLAYER)			! IN/OUT: Prop.HUM produced on HUM decompn 
	REAL DFACT(MAXLAYER)			! IN/OUT: Denitrification factor
	INTEGER DOCMODEL			! IN:DOC model (on or off)
	INTEGER DOC_OFF				! DOC model off
	INTEGER DOC_ON				! DOC model on
	DATA DOC_ON,DOC_OFF /1,2/

	REAL DPMRATE(MAXLAYER)			! IN/OUT: Rate constant for DPM decompn
	REAL DRAINW(MAXLAYER)		! IN:Water drained from this layer (mm/layer)
	REAL GAMMA(MAXLAYER)			! IN/OUT: Prop.BIO produced on HUM decompn 
	REAL HPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL HPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from humus decompositn
	REAL HRATE(MAXLAYER)			! IN/OUT: Rate constant for BIO decompn
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for HUM decompn/yr
	REAL HY(MAXSOIL,MAXLAYER)	! IN:Stable N:C ratio of BIO & HUM pools
	REAL HZ1(MAXLAYER)				! IN/OUT: N:C ratio for steady state
	REAL ICFACTOR(MAXLAYER)			! IN:Adjustment in rate needed to 
									!    achieve measured NPP and TOC
	REAL IOMARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium IOM (kgC/ha/layer)
	REAL LEACHDOC(MAXLAYER)		! IN/OUT:DOC leached from the layer (kgC/ha/layer/timestep)
	REAL LEACHDON(MAXLAYER)		! IN/OUT:DON leached from the layer (kgN/ha/layer/timestep)
      INTEGER LUARRAY(MAXSOIL)	! IN/OUT:Land use before equilibrium 
	REAL MOBDOC(MAXLAYER)		! IN:Mobile DOC in each layer (kgC/ha)
	REAL MOBDON(MAXLAYER)		! IN:Mobile DON in each layer (kgN/ha)
 	INTEGER NSOIL				! IN:Soil code number
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)	! IN/PUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)	! IN/PUT: pH above which decomp.max.
	REAL REFIL(MAXLAYER)		! Water deficit in the layer mm/layer
	REAL RPMRATE(MAXLAYER)			! IN/OUT: Rate constant for RPM decompn
	REAL SLEACH(MAXLAYER)		! INOUT: Nitrate N leached from this layer 
								!		 to the next (kgN/ha/layer/timestep)
	REAL SLEA15(MAXLAYER)		! INOUT: Nitrate N15 leached from this layer 
								!		 to the next (kgN15/ha/layer/timestep)
	CHARACTER*40 SNAME(MAXSOIL)	! IN/OUT:Soil name
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
 	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL TAM(MAXLAYER)			! IN:15N/N in the ammonium in the layer
	REAL TIN(MAXLAYER)			! IN:15N/N in the nitrate in the layer
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium TOC (kgC/ha/layer)
	REAL WLEACH					! IN:N leached from profile (kgN/ha)
	REAL WLEACH15				! IN:N15 leached from profile (kgN15/ha)
								!		 to the next (kgN15/ha/layer/timestep)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
C
C DOC model types
C
	DOC_ON=1
	DOC_OFF=2
C
C Retrieve soil CN characteristics 
C
      ISAVE=0
      CALL SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,CONVER_F,
     &             SECONDS,NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)
C
C Calculate leaching of nitrate
C
	WLEACH=0
      CALL LEACH_NITRATE(SOILN,SOIL15,SLEACH,SLEA15,DRAINW,WMAX,
     &                         SOILW,CRIT,NSOIL,REFIL)
C
C Calculate leaching of DOC
C
      IF(DOCMODEL.EQ.DOC_ON)
     &       CALL LEACH_DOC(MOBDOC,MOBDON,LEACHDOC,LEACHDON,
     &                      DRAINW,WMAX,SOILW,REFIL,
     &                      PHARRAY,NSOIL)
C
C Calculate leaching of ammonium
C
      CALL LEACH_NH4(AMMN,AMMN15,SLEACH,SLEA15,DRAINW,WMAX,
     &                         SOILW,CRIT,NSOIL,REFIL)
C
C Proportion N15:N
C
      CALL PROPN15(SOILN,SOIL15,TIN,AMMN,AMMN15,TAM)
C
C Sum leached N and N15
C
      DRAINW(MAXLAYER1)=DRAINW(MAXLAYER1)+BYRAIN
      WLEACH=WLEACH+SLEACH(MAXLAYER1)
      WLEACH15=WLEACH15+SLEA15(MAXLAYER1)
      IF(DRAINW(MAXLAYER1).NE.0)THEN
       CONC=WLEACH*100/DRAINW(MAXLAYER1)
       CONC15=WLEACH15*100/DRAINW(MAXLAYER1)
      ELSE
       CONC=0
       CONC15=0
      ENDIF
C
C Save soil CN characteristics 
C
      ISAVE=1
      CALL SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,CONVER_F,
     &             SECONDS,NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)
C
C Leave PHYSICAL_SUNDIAL_SOILCN
C
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE RESTART_SUNDIAL_SOILCN(
C INPUTS: Time, environment	
     &                               EQMODEL,DOCMODEL,SECONDS,
     &                               SPARMODEL,NSOIL,
C OUTPUTS: Soil C&N     
     &                               TOC,IOM,HZ1,CRIT,
C INPUTS/OUTPUTS: Soil C&N
     &                               MOBDOC,MOBDON,
     &                               ILU,CLAY,BULKDENS,
     &							   ICFACTOR,SOILPH)
C
C Subroutine to initialise soil C and N
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXSOIL				! Max.no.of soil types
	INTEGER MAXORGM				! Max.no.of organic manure applications allowed
	INTEGER MAXFERT				! Max.no.of fertiliser applications allowed 
	PARAMETER(MAXSOIL=50,MAXORGM=52,MAXFERT=5)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER IL					! Layer counter
	INTEGER MAXFERTSTEPS		! Maximum number of timesteps 
								!	over which fertiliser is added
      PARAMETER (MAXFERTSTEPS=30*60*24*7*3)	! Timesteps = 3 wks (max.30 secs)
C
C ...Model descriptors
C
	INTEGER DOC_OFF				! DOC model off
	INTEGER DOC_ON				! DOC model on
	DATA DOC_ON,DOC_OFF /1,2/
	INTEGER SPARFILE			! Soil parameters read in from file
	INTEGER SPARCALC			! Soil parameters calculated from TOC etc
	DATA SPARFILE,SPARCALC /1,2/
	INTEGER EQMODEL				! IN:Type of equilibrium run (NPP,TOC or both)
	INTEGER EQNPP				! Model initialised using plant inputs 
								!    from measured NPP
	INTEGER EQTOC				! Model initialised using measured TOC
	INTEGER EQNPPTOC			! Model initialised using both plant inputs
								!    from measured NPP and measured TOC
	INTEGER EQHILLIER,EQJONES			! Model initialised using TOC and
								!    the Hillier solver to get C pools and plant inputs
	DATA EQNPP,EQTOC,EQNPPTOC,EQHILLIER,EQJONES /1,2,3,4,5/
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: % clay content in this layer
      INTEGER LUARRAY(MAXSOIL)	! IN/OUT:Land use before equilibrium 
C
C Variables passed to/from this subroutine
C
      INTEGER ISAVE				! OUT:Code to save or retrieve variables
C
C ...Model descriptors
C
	INTEGER DOCMODEL			! IN:DOC model (on or off)
	INTEGER SPARMODEL			! IN:Soil parameter model (from file or calc)
C
C ...Timing factors
C
      REAL CONVER_F				! IN/OUT:Conversion to correct timestep
	REAL SECONDS				! IN:Number of seconds in one timestep
C
C ...Soil factors
C
	REAL TOC(MAXLAYER)			! IN:Total organic C (kgC/ha/layer)
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
      REAL CLAY(MAXLAYER)			! IN:Clay content of the layer (%)
	REAL BULKDENS(MAXLAYER)		! IN:Bulk density of the layer (g/cm3)
	INTEGER NSOIL				! IN:Soil code number
      INTEGER NUMSOIL				! IN/OUT:Number of soils defined
	CHARACTER*40 SNAME(MAXSOIL)	! IN/OUT:Soil name
	REAL ANIT					! IN/OUT:NH4-N nitrified(kgN/ha)
	REAL ANIT15					! IN/OUT:NH4-N15 nitrified (kgN/ha)
	REAL FANIT					! IN/OUT:Fertiliser NH4-N nitrified(kgN/ha)
	REAL FANIT15				! IN/OUT:Fertiliser NH4-N15 nitrified(kgN/ha)
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)	! IN/PUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)	! IN/PUT: pH above which decomp.max.
	REAL HY(MAXSOIL,MAXLAYER)	! IN/OUT:Stable N:C ratio of BIO & HUM pools
	REAL BPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL HPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL BPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from biomass decompositn
	REAL HPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from humus decompositn
	REAL BIOP(MAXSOIL,MAXLAYER)		! IN/OUT: BIO/TOC 
	REAL CRIT(MAXSOIL,MAXLAYER)		! IN/OUT: Minimum nitrate level (Sozanska)
	REAL BIORATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for BIO decompn/yr
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for HUM decompn/yr
	REAL ALPHA(MAXLAYER)			! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL BETA(MAXLAYER)				! IN/OUT: Prop.HUM produced on BIO decompn
	REAL GAMMA(MAXLAYER)			! IN/OUT: Prop.BIO produced on HUM decompn 
	REAL DELTA(MAXLAYER)			! IN/OUT: Prop.HUM produced on HUM decompn 
	REAL HZ1(MAXLAYER)				! IN/OUT: N:C ratio for steady state
	REAL DFACT(MAXLAYER)			! IN/OUT: Denitrification factor
	REAL DPMRATE(MAXLAYER)			! IN/OUT: Rate constant for DPM decompn
	REAL RPMRATE(MAXLAYER)			! IN/OUT: Rate constant for RPM decompn
	REAL BRATE(MAXLAYER)			! IN/OUT: Rate constant for HUM decompn
	REAL HRATE(MAXLAYER)			! IN/OUT: Rate constant for BIO decompn
	REAL IOMARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium IOM (kgC/ha/layer)
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium TOC (kgC/ha/layer)
	INTEGER ILU						! IN:Counter for land use 
	REAL ICFACTOR(MAXLAYER)			! IN:Adjustment in rate needed to 
									!    achieve measured NPP and TOC
	REAL SOILPH(MAXLAYER)		! IN/OUT:pH of soil in this layer
C
C ...Manure factors
C
      REAL FPROPX(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
	REAL FYMCX(MAXORGM)			! IN/OUT:Proportion of C in FYM
	REAL FYMNX(MAXORGM)			! IN/OUT:Proportion of Organic N in FYM
	REAL FYMAX(MAXORGM)			! IN/OUT:Proportion of ammonium-N in FYM
	REAL FYMWX(MAXORGM)			! IN/OUT:Amount of water added in FYM
	REAL FYMLOSS(MAXORGM)		! IN/OUT:Prop.FYM lost by volatn/timestep
	REAL FYMPOS(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
	REAL FYMSTART(MAXORGM)		! IN/OUT:Amount of rainfall in 1 week, 
								!		 below which volatilisation will occur
C
C ...Fertiliser factors
C
	REAL BYRAIN					! IN/OUT:Rain lost by bypass flow mm/timestep	
	REAL FERTADD(MAXFERTSTEPS,4)	! INOUT:Fertiliser N added in this time 
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
	REAL FERTADD15(MAXFERTSTEPS,4)	! INOUT:Fertiliser N15 added in this time
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
C
C ...Dissolved organic matter factors
C
	REAL MOBDOC(MAXLAYER)		! IN:Mobile DOC in each layer (kgC/ha)
	REAL MOBDON(MAXLAYER)		! IN:Mobile DON in each layer (kgN/ha)
C
C Get soil parameters
C
  	CALL GET_SUNDIAL_SOILPARS(NSOIL,TOC,IOM,CRIT,ILU,
     &                          CLAY,BULKDENS,SPARMODEL)
C
C Retrieve soil CN characteristics 
C
      ISAVE=0
      CALL SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,CONVER_F,
     &             SECONDS,NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)
	CALL SAVE_FERTADD(SECONDS,FERTADD,FERTADD15,ISAVE)
C
C Set time factors according to SECONDS
C
      CALL SETTIME_CN(SECONDS,CONVER_F)
C
C Set starting value of organic matter
C
      CALL STARTORG(NSOIL,TOC,TOCARRAY,IOM,IOMARRAY)
C
C Set BIO and HUM decompostion parameters
C
      CALL OMPROP(NSOIL,HY,HZ1,ALPHA,BETA,GAMMA,DELTA,
     &                  BPART,HPART,BPROP,HPROP,PHARRAY)
C
C Set Rate constants for decomposition of RO (=R), BIO (=B) and HUM (=H)
C
	CALL GET_DECOMP_RATECONSTS(CONVER_F,NSOIL,BIORATE,HUMRATE,
     &                           DPMRATE,RPMRATE,BRATE,HRATE)	  
C
C Set residual DOC content
C
      IF(DOCMODEL.EQ.DOC_ON)CALL SETDOC(NSOIL,MOBDOC,MOBDON)
C
C Set initial fertiliser additions
C
      CALL SETFERT(FERTADD,FERTADD15)									
C
C Pass soil pH out of routine
C
      DO 100 IL=1,MAXLAYER
	  SOILPH(IL)=PHARRAY(NSOIL,IL)
100   CONTINUE
C
C Save soil CN characteristics 
C
      ISAVE=1
      CALL SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,CONVER_F,
     &             SECONDS,NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)
	CALL SAVE_FERTADD(SECONDS,FERTADD,FERTADD15,ISAVE)
C
C Leave RESTART_SOIL_SUNDIALCN
C
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE RUN_ROTHC_AND_N_TO_EQUILIBRIUM(
C INPUTS: Time, environment	
     &                               ICMODEL,INMODEL,EQMODEL,SECONDS,
     &                               NSOIL,
C OUTPUTS: Soil C&N     
     &                               TOC,IOM,HZ1,CRIT,
C INPUTS/OUTPUTS: Soil C&N
     &                               DPMCARB0,DPMNIT0,DPMNLAB0,
     &                               RPMCARB0,RPMNIT0,RPMNLAB0,
     &                               BCARB0,BNIT0,BNLAB0,
     &                               HCARB0,HNIT0,HNLAB0,
     &                               PI_CEQ_MON,
     &                               LTA_AWC,LTA_TEMP,ITFUNC,IMFUNC,
     &                               WMAX,WSAT,ICFACTOR,
     &                               DPM_RPM_RATIO,DPMCTON,RPMCTON,
     &                               SOILPH)
C
C Subroutine to equilibrium run for sundial soil C and N
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXORGM				! Max.no.of org.manure applications allowed
	PARAMETER(MAXORGM=52)
	INTEGER MAXSOIL				! Max.no.of soil types
	PARAMETER(MAXSOIL=50)

	INTEGER ICFIXED				! Initialisation of C is fixed
	INTEGER ICROTHCEQ			! Initialisation of C by RothC 
								!  equilibrium run
	DATA ICFIXED,ICROTHCEQ /1,2/
	INTEGER INSTABLECN			! Initisalisation of N by assuming 
								!  steady state (after Bradbury)
      INTEGER IL					! Local layer counter
	INTEGER INPASSCN			! Initialisation of N by passing the
								!  C:N ratio of DPM and RPM
	DATA INSTABLECN,INPASSCN /1,2/
C
C Variables passed to/from this subroutine
C
	REAL ALPHA(MAXLAYER)			! IN/OUT:Prop.BIO produced on BIO decompn 
	REAL ANIT					! IN/OUT:NH4-N nitrified(kgN/ha)
	REAL ANIT15					! IN/OUT:NH4-N15 nitrified (kgN/ha)
	REAL BCARB0(MAXLAYER)		! IN:C in soil biomass (kgC/ha/layer)
	REAL BETA(MAXLAYER)				! IN/OUT:Prop.HUM produced on BIO decompn
	REAL BIOP(MAXSOIL,MAXLAYER)		! IN/OUT:BIO/TOC 
	REAL BIORATE(MAXSOIL,MAXLAYER)	! IN/OUT:Rate constant for BIO decompn/yr
      REAL BNIT0(MAXLAYER)		! IN:N in soil biomass (kgN/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN:N15 in soil humus (kgN15/ha/layer)
	REAL BPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL BPROP(MAXSOIL,MAXLAYER)	! IN/OUT:BIO/HUM from biomass decompositn
	REAL BRATE(MAXLAYER)			! IN/OUT:Rate constant for HUM decompn
	REAL BYRAIN						! IN/OUT:Rain lost by bypass flow mm/timestep	
      REAL CLAY(MAXLAYER)				! OUT:Clay content of the layer (%)
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:% clay content in this layer
      REAL CONVER_F					! IN/OUT:Conversion to correct timestep
	REAL CRIT(MAXSOIL,MAXLAYER)		! IN/OUT:Minimum nitrate level (Sozanska)
	REAL DELTA(MAXLAYER)			! IN/OUT:Prop.HUM produced on HUM decompn 
	REAL DFACT(MAXLAYER)			! IN/OUT:Denitrification factor
      REAL DPM_RPM_RATIO				! IN/OUT:Ratio of DPM:RPM. If set 
									!        to zero, this will be worked out using a call to SET_DPMRPMRATIO
	REAL DPMCARB0(MAXLAYER)			! IN:C in decomposable PM (kgC/ha/layer)
      REAL DPMCTON(MAXLAYER)			! IN:Ratio of C:N in DPM
	REAL DPMNIT0(MAXLAYER)			! IN:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)			! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL DPMRATE(MAXLAYER)			! IN/OUT:Rate constant for DPM decompn
	INTEGER EQMODEL					! IN:Type of equilibrium run (NPP,TOC or both)
	REAL FANIT						! IN/OUT:Fertiliser NH4-N nitrified(kgN/ha)
	REAL FANIT15					! IN/OUT:Fertiliser NH4-N15 nitrified(kgN/ha)
      REAL FPROPX(MAXORGM)			! IN/OUT:Proportion of FYM added to top 25cm
	REAL FYMAX(MAXORGM)				! IN/OUT:Proportion of ammonium-N in FYM
	REAL FYMCX(MAXORGM)				! IN/OUT:Proportion of C in FYM
	REAL FYMLOSS(MAXORGM)			! IN/OUT:Prop.FYM lost by volatn/timestep
	REAL FYMNX(MAXORGM)				! IN/OUT:Proportion of Organic N in FYM
	REAL FYMPOS(MAXORGM)			! IN/OUT:Proportion of FYM added to top 25cm
	REAL FYMSTART(MAXORGM)			! IN/OUT:Amount of rainfall in 1 week, 
									!		 below which volatilisation will occur
	REAL FYMWX(MAXORGM)				! IN/OUT:Amount of water added in FYM
	REAL GAMMA(MAXLAYER)			! IN/OUT: Prop.BIO produced on HUM decompn 
	REAL HCARB0(MAXLAYER)			! IN:C in soil humus (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)			! IN:N in soil humus (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)			! IN:N15 in soil biomass (kgN15/ha/layer)
	REAL HPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL HPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from humus decompositn
	REAL HRATE(MAXLAYER)			! IN/OUT: Rate constant for BIO decompn
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for HUM decompn/yr
	REAL HY(MAXSOIL,MAXLAYER)		! IN/OUT:Stable N:C ratio of BIO & HUM pools
	REAL HZ1(MAXLAYER)				! IN/OUT: N:C ratio for steady state
  	REAL ICFACTOR(MAXLAYER)			! IN:Adjustment in rate needed to	
									!    achieve measured NPP and TOC	
      INTEGER IDUMP				    ! OUT:Dump data to file or retrieve 
							    	! from file. (1=dump;0=retrieve)
	INTEGER IMFUNC					! IN:Choice of moisture rate modifier 
									! (ROTHC or HADLEY)
	INTEGER ICMODEL				! IN:Type of C initialisation model 
	INTEGER INMODEL				! IN:Type of N initialisation model
	REAL IOM(MAXLAYER)				! IN:Inert organic C (kgC/ha/layer)
	REAL IOMARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium IOM (kgC/ha/layer)
	INTEGER ITFUNC					! IN:Choice of temperature rate modifier 
									! (ROTHC or HADLEY)
      INTEGER ISAVE					! OUT:Code to save or retrieve variables
	REAL LTA_AWC(12,MAXLAYER)		! IN:Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)		! IN:Average air temp this month (deg.C)
      INTEGER LUARRAY(MAXSOIL)		! IN/OUT:Land use before equilibrium 
	INTEGER NSOIL					! IN:Soil code number
      INTEGER NUMSOIL					! IN/OUT:Number of soils defined
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)! IN/OUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)! IN/OUT: pH above which decomp.max.
	REAL PI_CEQ_MON(12,MAXLAYER)	! IN/OUT: Equilibrium plant C input each 
									!         month to this layer (kgC/ha/month/layer)
	REAL RPMCARB0(MAXLAYER)			! IN:C in resistant PM (kgC/ha/layer)
      REAL RPMCTON(MAXLAYER)					! IN/OUT:Ratio of C:N in RPM
	REAL RPMNIT0(MAXLAYER)			! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)			! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL RPMRATE(MAXLAYER)			! IN/OUT: Rate constant for RPM decompn
	REAL SECONDS					! IN:Number of seconds in one timestep
	CHARACTER*40 SNAME(MAXSOIL)		! IN/OUT:Soil name
	REAL SOILPH(MAXLAYER)		! IN/OUT:pH of soil in this layer
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL TOC(MAXLAYER)				! IN:Total organic C (kgC/ha/layer)
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium TOC (kgC/ha/layer)
	REAL WMAX(MAXLAYER)				! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)				! IN:Available water at saturation (mm/layer)
C
C Retrieve soil CN characteristics 
C
      ISAVE=0
      CALL SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,CONVER_F,
     &             SECONDS,NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)
C
C Pass soil pH out of routine and set clay content
C
	DO 100 IL=1,MAXLAYER
	  CLAY(IL)=CLARRAY(NSOIL,IL)
	  SOILPH(IL)=PHARRAY(NSOIL,IL)
100   CONTINUE
C
C Set the residual N and C in the RO,BIO and HUM pools
C
      IF(INMODEL.EQ.INSTABLECN)THEN
        CALL GET_CTON_BRADBURY(DPMCTON,RPMCTON,HZ1)
	ENDIF
	IF(ICMODEL.EQ.ICROTHCEQ)THEN
        CALL SETORGM_EQRUN(NSOIL,DPMRATE,RPMRATE,TOC,IOM,CONVER_F,
     &                   BRATE,BIORATE,HRATE,HUMRATE,HZ1,HY,
     &                   BCARB0,HCARB0,DPMCARB0,RPMCARB0,
     &                   BNIT0,HNIT0,DPMNIT0,RPMNIT0,
     &                   HNLAB0,BNLAB0,DPMNLAB0,RPMNLAB0,
     &                   LUARRAY,CLARRAY,
     &                   PHARRAY,PHP1ARRAY,PHP2ARRAY,
     &                   EQMODEL,PI_CEQ_MON,ICFACTOR,LTA_AWC,LTA_TEMP,
     &                   ITFUNC,IMFUNC,WMAX,WSAT,
     &                   DPM_RPM_RATIO,DPMCTON,RPMCTON,ALPHA,BETA,CLAY)	
      ELSEIF(ICMODEL.EQ.ICFIXED)THEN
        CALL SETORGM_FIXED(NSOIL,DPMRATE,RPMRATE,TOC,IOM,CONVER_F,
     &                   BRATE,BIORATE,HRATE,HUMRATE,HZ1,HY,
     &                   BCARB0,HCARB0,DPMCARB0,RPMCARB0,
     &                   BNIT0,HNIT0,DPMNIT0,RPMNIT0,
     &                   HNLAB0,BNLAB0,DPMNLAB0,RPMNLAB0,ICFACTOR,
     &                   DPM_RPM_RATIO,DPMCTON,RPMCTON)	
	ENDIF 
C
C Save soil CN characteristics 
C
      ISAVE=1
      CALL SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,CONVER_F,
     &             SECONDS,NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)
C
C Dump to file so have equilibrium position when model is restarted
C
	IDUMP = 1
	CALL DUMP_SOILCN(IDUMP)
C
C Leave equilibrium run for sundial soil N and N
C
	END
C
C-------------------------------------------------------------
C
C INTERNAL SUBROUTINES
C
C------------------------------------------------------------
C
      SUBROUTINE ADDATM(ATM,SOILN)
C
C Subroutine to add atmospheric dry deposition
C
      IMPLICIT NONE
C
C Variable local to this subroutine
C
	INTEGER MAXLAYER
	PARAMETER (MAXLAYER=60)
C
C Variables passed to/from calling subroutine
C
      REAL SOILN(MAXLAYER),ATM
C
C Add NO3 at a rate ATM
C
      SOILN(1)=SOILN(1)+ATM
C
C Leave ADDATM
C
      RETURN
      END
C
C---------------------------------------------------------
C
      SUBROUTINE ADDFERT(FERTNF,TFERTNF,ILABNF,IVOLNF,FERTADD,FERTADD15,
     &                   SECONDS)
C
C Subroutine to add fertilizer
C
      IMPLICIT NONE
C
C Variable local to this subroutine
C
	INTEGER MAXFERT				! Max.no.of fertiliser applications allowed 
	PARAMETER(MAXFERT=5)
	INTEGER MAXFERTSTEPS		! Maximum number of timesteps 
								!	over which fertiliser is added
      PARAMETER (MAXFERTSTEPS=30*60*24*7*3)	! Timesteps = 3 wks (max.30 secs)
	INTEGER IT					! Local counter for fertiliser type
	INTEGER IS					! Local counter for steps of fert.addition
	INTEGER INSTEPS				! Timesteps over which fert.is added (= 3 wks)
	REAL FERTVOL(4)				! Fertiliser N added  
								!   of type 1=nitrate,2=amm.sulphate,
								!           3=other amm.salt,4=Urea
C
C Variables passed to/from calling subroutine
C
	REAL FERTADD(MAXFERTSTEPS,4)	! INOUT:Fertiliser N added in this time 
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
	REAL FERTADD15(MAXFERTSTEPS,4)	! INOUT:Fertiliser N15 added in this time
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
	REAL FERTNF					! IN:Amount of fertiliser applied (kgN/ha)
	REAL TFERTNF(3)				! IN:Prop.NO3,NH4,urea in fertiliser
	INTEGER ILABNF				! IN:Labelling on fertiliser? 0=No 1=Yes
	INTEGER IVOLNF				! IN:Ammonium salt (0=sulphate, 1=other)
	REAL SECONDS				! IN:Number of seconds in one timestep
C
C Add fertiliser to pool of fertiliser not yet incorporated in the soil
C
      INSTEPS=INT(3*(7*24*60*60)/SECONDS)
	IF(INSTEPS.LT.1)INSTEPS=1
	FERTVOL(1)=FERTNF*TFERTNF(1)/100
	IF(IVOLNF.EQ.0)FERTVOL(2)=FERTNF*TFERTNF(2)/100
	IF(IVOLNF.EQ.1)FERTVOL(3)=FERTNF*TFERTNF(2)/100
	FERTVOL(4)=FERTNF*TFERTNF(3)/100
      DO 100 IS=1,INSTEPS   
        FERTADD(IS,1)=FERTADD(IS,1)+FERTVOL(1)/INSTEPS
        IF(ILABNF.EQ.1)FERTADD15(IS,1)=FERTADD15(IS,1)
     &					+FERTVOL(1)/INSTEPS
100   CONTINUE
      DO 200 IT=2,4
	  FERTADD(1,IT)=FERTVOL(IT)
200   CONTINUE
C
C Leave ADDFERT
C
      RETURN
      END
C
C--------------------------------------------------------------
C
      SUBROUTINE ADDFYM(FPROPX,JORGMNF,FYMCX,FYMNX,FYMAX,FYMWX,
     &                  RAIN,FYMSTART,FYMLOSS,ORGMANF,HY,
     &                  FYMPOS,HCARB0,HNIT0,DPMCARB0,RPMCARB0,
     &                  DPMNIT0,RPMNIT0,HNLAB0,DPMNLAB0,RPMNLAB0,
     &                  AMMN,VOLAT,IOLABNF,AMMN15,VOLAT15,NSOIL,
     &                  FYMFERT,FYMFERT15,CONVER_F)
C
C Subroutine to add FYM
C
      IMPLICIT NONE
C
C Variable local to this subroutine
C
	INTEGER MAXSOIL				! Max.no.of soil types
	INTEGER MAXORGM				! Max.no.of org.manure applications allowed
	PARAMETER(MAXSOIL=50,MAXORGM=52)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	INTEGER MAXDEPTH			! Maximum depth of the soil profile
	INTEGER IDEPTH				! Depth of this layer
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
	REAL FPROP					! Proportion of FYM added to top 25cm
	REAL FYMC					! Proportion of C in FYM
	REAL FYMN					! Proportion of Organic N in FYM
	REAL FYMA					! Proportion of ammonium-N in FYM
	REAL FLOSS					! Prop.FYM lost by volatn/timestep
	REAL FYMW					! Amount of water added in FYM
	REAL TFYMN					! Total N added in FYM (kgN/ha)
	REAL TFYMC					! Total C added in FYM (kgC/ha)
	REAL TFYMA					! Total ammonium added in FYM (kgN/ha)
	REAL HFYMC					! Humus C added in FYM (kgC/ha)
	REAL HFYMN					! Humus N added in FYM (kgN/ha)
	REAL RFYMC					! Plant debris C added in FYM (kgC/ha)
	REAL RFYMN					! Plant debris N added in FYM (kgN/ha)
	REAL PLAYER					! Size of layer compared to 25cm
	INTEGER I					! Local counter variable
	REAL SPL					! Split of manure into this layer
	REAL SPLIT					! Split of manure into this layer
C
C Variables passed to/from calling subroutine
C ... Weather factors
C
	REAL RAIN					! IN: Rainfall (mm/timestep)
C
C ...Manure factors
C
	INTEGER JORGMNF				! IN:Type of manure application
	INTEGER IOLABNF				! IN:Labelling on manure? 0=No 1=Yes
	REAL ORGMANF				! IN:Amount of manure applied (t/ha)
	INTEGER NSOIL				! IN:Soil code number
      REAL FPROPX(MAXORGM)		! IN/OUT:Proportion of FYM that is decomposed organic matter
	REAL FYMCX(MAXORGM)			! IN/OUT:Proportion of C in FYM
	REAL FYMNX(MAXORGM)			! IN/OUT:Proportion of Organic N in FYM
	REAL FYMAX(MAXORGM)			! IN/OUT:Proportion of ammonium-N in FYM
	REAL FYMWX(MAXORGM)			! IN/OUT:Amount of water added in FYM
	REAL FYMLOSS(MAXORGM)		! IN/OUT:Prop.FYM lost by volatn/timestep
	REAL FYMPOS(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
	REAL FYMSTART(MAXORGM)		! IN/OUT:Amount of rainfall in 1 week, 
								!		 below which volatilisation will occur
C
C ...Soil factors
C
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN:N in soil humus (kgN/ha/layer)
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN:N15 in soil biomass (kgN15/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
	REAL VOLAT					! IN:N lost by volatilisation (kgN/ha)
	REAL VOLAT15				! IN:N15 lost by volatilisation (kgN15/ha)
	REAL FYMFERT				! IN:N input by fertiliser & FYM (kgN/ha)
	REAL FYMFERT15				! IN:N15 input by fertiliser & FYM (kgN15/ha)
      REAL CONVER_F				! IN/OUT:Conversion to correct timestep
	REAL HY(MAXSOIL,MAXLAYER)	! IN/OUT:Stable N:C ratio of BIO & HUM pools
C
C Set parameters of current manure application
C
      IF(JORGMNF.LE.0)RETURN
      FPROP=FPROPX(JORGMNF)
      FYMC=FYMCX(JORGMNF)
      FYMN=FYMNX(JORGMNF)
      FYMA=FYMAX(JORGMNF)
      FYMW=FYMWX(JORGMNF)
C
C Adjust amount lost by volatilisation according to rainfall (excluding
C water added with FYM)
C
      IF(RAIN.LT.(FYMSTART(JORGMNF)*CONVER_F))THEN
       FLOSS=FYMLOSS(JORGMNF)
      ELSE
       FLOSS=0.0
      ENDIF
C
C Add water in FYM to incoming rainfall
C
      RAIN=RAIN+FYMW*ORGMANF
C
C Calculate manure applcation rates...
C TFYMN : Total organic N in FYM
C TFYMC : Total organic C in FYM
C HFYMN : N entering HUM
C HFYMC : C entering HUM
C RFYMN : N entering RO
C RFYMC : C entering RO
C
      TFYMN=FYMN*ORGMANF
      TFYMC=FYMC*ORGMANF
	TFYMA=FYMA*ORGMANF
      HFYMN=FPROP*TFYMN
      HFYMC=HFYMN/HY(NSOIL,1)
      RFYMN=TFYMN-HFYMN
      RFYMC=TFYMC-HFYMC
C
C Divide organic part of application between top 50cm
C Layer 1 = 0-25cm; Layer 2 = 25-50cm=0%
C
      DO 20 I=1,MAXLAYER1
	  IDEPTH=I*MAXDEPTH/(MAXLAYER1*1.)
        SPL=FYMPOS(JORGMNF)
	  PLAYER=MAXDEPTH/(MAXLAYER1*25.)
        IF(IDEPTH.LE.25)THEN
         SPLIT=SPL*PLAYER
        ELSEIF(IDEPTH.GT.25.AND.IDEPTH.LE.50)THEN
         SPLIT=(1.0-SPL)*PLAYER
        END IF
C
C Put C and N into HUM and RO pools
C
        IF(IDEPTH.LE.50)THEN
          HCARB0(I)=HCARB0(I)+SPLIT*HFYMC
          HNIT0(I)=HNIT0(I)+SPLIT*HFYMN
          DPMCARB0(I)=DPMCARB0(I)+SPLIT*RFYMC*0.5
          RPMCARB0(I)=RPMCARB0(I)+SPLIT*RFYMC*0.5
          DPMNIT0(I)=DPMNIT0(I)+SPLIT*RFYMN*0.5
          RPMNIT0(I)=RPMNIT0(I)+SPLIT*RFYMN*0.5
C
C Add any label to HUM and RO pools
C
          IF(IOLABNF.EQ.1)THEN
            HNLAB0(I)=HNLAB0(I)+SPLIT*HFYMN
            DPMNLAB0(I)=DPMNLAB0(I)+SPLIT*RFYMN*0.5
            RPMNLAB0(I)=RPMNLAB0(I)+SPLIT*RFYMN*0.5
          ENDIF
	  ENDIF
20    CONTINUE
C
C Put N into ammounium pool
C
      AMMN(1)=AMMN(1)+((TFYMA-(FLOSS*TFYMA)))
      VOLAT=VOLAT+(FLOSS*TFYMA)
      FYMFERT=FYMFERT+(TFYMN+TFYMA)
C
C Add any label to ammounium pool
C
      IF(IOLABNF.EQ.1)THEN
        AMMN15(1)=AMMN15(1)+TFYMA-(FLOSS*TFYMA)
        VOLAT15=VOLAT15+(FLOSS*TFYMA)
        FYMFERT15=FYMFERT15+(TFYMN+TFYMA)
      ENDIF
C
C Leave ADDFYM
C
      RETURN
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE BYPASSFLOW(
C ...Weather factors			
     &                      RAIN,BYRAIN,
C ...Timing factors
     &					  SECONDS,CONVER_F,
C ...Fertiliser factors
     &					  FERTADD,FERTADD15,THISFERT,THISFERT15,
C ...Soil factors
     &					  SOILN,SOIL15,AMMN,AMMN15,
     & 					  TAM,TIN,WLEACH,WLEACH15)

C
C Subroutine to calculate fertilizer additions into soil and leaching loss via bypass flow.
C
      IMPLICIT NONE
C
C Variable local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXFERTSTEPS		! Maximum number of timesteps 
								!	over which fertiliser is added
      PARAMETER (MAXFERTSTEPS=30*60*24*7*3)	! Timesteps = 3 wks (max.30 secs)

      REAL RBY					! Fraction of added fertiliser
								!   lost by bypass flow
	REAL RCRIT					! Critical level of rainfall above which 
								!   bypass flow will occur
	REAL BYFRAC					! Fraction of rainfall above critical level 
								!   that will be lost by bypass flow
      PARAMETER (RBY=0.015, RCRIT=15, BYFRAC=0.5)
	INTEGER INSTEPS				! Number of timesteps in 3 weeks 
	REAL BYPASS					! N lost by bypass flow (kgN/ha/timestep)
	REAL BYPASS15				! N15 lost by bypass flow (kgN15/ha/timestep)
	INTEGER IS					! Local counter for steps of fert.addition
C
C Variables passed to/from calling subroutine
C
C ...Weather factors
C
	REAL RAIN					! IN: Rainfall (mm/timestep)
	REAL BYRAIN					! IN/OUT:Rain lost by bypass flow mm/timestep	
C
C ...Timing factors
C
	REAL SECONDS				! IN:Number of seconds in one timestep
      REAL CONVER_F				! IN/OUT:Conversion to correct timestep
C
C ...Fertiliser factors
C
	REAL FERTADD(MAXFERTSTEPS,4)	! INOUT:Fertiliser N added in this time 
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
	REAL FERTADD15(MAXFERTSTEPS,4)	! INOUT:Fertiliser N15 added in this time
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
	REAL THISFERT				! OUT:N input by fertiliser (kgN/ha)
	REAL THISFERT15				! OUTN15 input by fertiliser (kgN15/ha)
C
C ...Soil factors
C
 	REAL SOILN(MAXLAYER)		! INOUT:Soil nitrate-N (kgN/ha/layer)
	REAL SOIL15(MAXLAYER)		! INOUT:Soil nitrate-N15 (kgN15/ha/layer)
	REAL AMMN(MAXLAYER)			! INOUT:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! INOUT:Soil ammonium-N15 (kgN15/ha/layer)
	REAL TAM(MAXLAYER)			! INOUT:15N/N in the ammonium in the layer
	REAL TIN(MAXLAYER)			! INOUT:15N/N in the nitrate in the layer
  	REAL WLEACH					! INOUT:N leached (kgN/ha)
	REAL WLEACH15				! INOUT:N15 leached (kgN15/ha)
C
C Saved variables
C
      SAVE 
C
C Initialize variables
C
      BYRAIN=0.0
      THISFERT=0.0
      THISFERT15=0.0
      INSTEPS=NINT(3*(7*24*60*60)/SECONDS)
C
C If rainfall has exceeded the critical level...
C
      IF(RAIN.GT.(RCRIT*CONVER_F))THEN
C
C ...calculate rainfall lost by bypass flow (= fraction (BYFRAC) of rainfall above critical minimum)
C
        IF(BYRAIN.LE.0.0001)BYRAIN=BYFRAC*(RAIN-(RCRIT*CONVER_F))
C 
C ...loose fertiliser by bypassflow, and wash the rest into the soil
C
        DO 100 IS=1,INSTEPS
C
C ...calculate the N lost by bypass flow (= fraction RBY of added fertilizer * bypass flow water)
C
          BYPASS=MIN(RBY*BYRAIN,1.0)*FERTADD(IS,1)
          BYPASS15=MIN(RBY*BYRAIN,1.0)*FERTADD15(IS,1)
C
C ...take off bypassed nitrate fertilizer and add the rest
C
          FERTADD(IS,1)=FERTADD(IS,1)-BYPASS
          SOILN(1)=SOILN(1)+FERTADD(IS,1)
          FERTADD15(1,1)=FERTADD15(IS,1)-BYPASS15
          SOIL15(1)=SOIL15(1)+FERTADD15(IS,1)
C
C ...save results
C
          THISFERT=THISFERT+FERTADD(IS,1)
          THISFERT15=THISFERT15+FERTADD15(IS,1)
          WLEACH=WLEACH+BYPASS
          WLEACH15=WLEACH+BYPASS15
	    FERTADD(IS,1)=0
	    FERTADD15(IS,1)=0
100     CONTINUE
C
C Otherwise, if rainfall has not exceeded critical level...
C
      ELSE
C
C ...add this timesteps fertiliser nitrate
C
        SOILN(1)=SOILN(1)+FERTADD(1,1)
        SOIL15(1)=SOIL15(1)+FERTADD15(1,1)
C
C ...save results
C
        THISFERT=THISFERT+FERTADD(1,1)
        THISFERT15=THISFERT15+FERTADD15(1,1)
        WLEACH=WLEACH+BYPASS
        WLEACH15=WLEACH+BYPASS15
C
C ...save the rest of the nitrate for adding next timestep
C
        DO 200 IS=1,INSTEPS-1
          FERTADD(IS,1)=FERTADD(IS+1,1)
          FERTADD15(IS,1)=FERTADD15(IS+1,1)
200     CONTINUE
        FERTADD(INSTEPS,1)=0
        FERTADD15(INSTEPS,1)=0
      ENDIF
C
C Reset the proportion of N:N15
C
      CALL PROPN15(SOILN,SOIL15,TIN,AMMN,AMMN15,TAM)
C >>> Commented out by Mark Richards, 15/02/2013 because BYRAIN is more or less the same water as 
C that lost due to flowprop in the SUNDIAL water model which results in the water being lost twice.
C      RAIN=RAIN-BYRAIN
C End of commenting out <<<
C	WRITE(33,*)BYRAIN
C
C Leave BYPASSFLOW
C
      END
C
CC--------------------------------------------------------
C
      SUBROUTINE CULTIV(LU1,DPMCARB0,DPMNIT0,DPMNLAB0,
     &                  RPMCARB0,RPMNIT0,RPMNLAB0,
     &                  BCARB0,BNIT0,BNLAB0,
     &                  HCARB0,HNIT0,HNLAB0,CULTDEPTH,VIGOUR,INVERT)

C
C Subroutine to cultivate the soil
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/
	INTEGER MAXLAYER			! No.of layers in the soil profile
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/

      REAL CTONDPM				! C:N ratio of DPM
      REAL CTONRPM				! C:N ratio of RPM
      REAL CTONBIO				! C:N ratio of BIO
      REAL CTONHUM				! C:N ratio of HUM
      INTEGER ICULT				! Bottom layer of cultivation				
      REAL DEPTH					! Depth of current layer (cm)
	REAL HUMTOPM				! Amount of humus passed to plant 
								! material on cultivation (kg C / ha)
      INTEGER IL					! Layer counter
	REAL PASS					! Proportion of humus passed to biomass, dpm and rpm
	DATA PASS /0.5/				! Soil parameter maximum proportion observed by West & Post (2002) = 0.5. 
								! More mechanistic description needed to obtain a more generic model.
	REAL PROPB					! Proportion bio/(bio&hum*dpm)
	REAL PROPD					! Proportion dpm/(bio&hum*dpm)
	REAL PROPR					! Proportion rpm/(bio&hum*dpm)
	REAL SOCNT					! Total SOC in this layer before tillage (kg C / ha)
	REAL N15TON					! N15 to N ratio of given pool
	REAL TOTDPMC				! Total DPM	C in tillage zone (kgC/ha)
	REAL TOTRPMC				! Total RPM C in tillage zone (kgC/ha)
	REAL TOTBIOC				! Total Bio C in tillage zone (kgC/ha)
	REAL TOTHUMC				! Total Hum C in tillage zone (kgC/ha)
	REAL TOTDPMN				! Total DPM	N in tillage zone (kgC/ha)
	REAL TOTRPMN				! Total RPM N in tillage zone (kgC/ha)
	REAL TOTBION				! Total Bio N in tillage zone (kgC/ha)
	REAL TOTHUMN				! Total Hum N in tillage zone (kgC/ha)
	REAL TOTDPMN15				! Total DPM	N15 in tillage zone (kgC/ha)
	REAL TOTRPMN15				! Total RPM N15 in tillage zone (kgC/ha)
	REAL TOTBION15				! Total Bio N15 in tillage zone (kgC/ha)
	REAL TOTHUMN15				! Total Hum N15 in tillage zone (kgC/ha)
C
C Variables passed to/from calling subroutine
C
	INTEGER INVERT				! IN:Invert / mix soil on cultivation? 0=No 1=Yes
      INTEGER ISAVE				! Code to save or retrieve variables
	REAL CULTDEPTH				! IN:Cultivation depth (cm)
	REAL DRRAT					! IN:DPM:RPM ratio
	INTEGER LU1					! IN:Land use before land use change
	REAL BCARB0(MAXLAYER)		! IN:C in soil biomass (kgC/ha/layer)
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
      REAL BNIT0(MAXLAYER)		! IN:N in soil biomass (kgN/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN:N in soil humus (kgN/ha/layer)
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN:N15 in soil biomass (kgN15/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN:N15 in soil humus (kgN15/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL VIGOUR					! IN: Vigour of cultivation (0-1) Determines the proportion of humus released to biomass,DPM and RPM pools
C
C Work out the bottom layer of the cultivation
C
      ICULT=((CULTDEPTH*MAXLAYER1)/MAXDEPTH)      
C
C Work out total C, N and N15 in the different pools to cultivation depth
C
	TOTDPMC=0
	TOTRPMC=0
	TOTBIOC=0
	TOTHUMC=0
	TOTDPMN=0
	TOTRPMN=0
	TOTBION=0
	TOTHUMN=0
	TOTDPMN15=0
	TOTRPMN15=0
	TOTBION15=0
	TOTHUMN15=0
C
C Mix up soil to depth of cultivation
C
      IF(INVERT.EQ.1)THEN
        DO 100 IL=1,ICULT
	    TOTDPMC=TOTDPMC+DPMCARB0(IL)
	    TOTRPMC=TOTRPMC+RPMCARB0(IL)
	    TOTBIOC=TOTBIOC+BCARB0(IL)
	    TOTHUMC=TOTHUMC+HCARB0(IL)
	    TOTDPMN=TOTDPMN+DPMNIT0(IL)
	    TOTRPMN=TOTRPMN+RPMNIT0(IL)
	    TOTBION=TOTBION+BNIT0(IL)
	    TOTHUMN=TOTHUMN+HNIT0(IL)
	    TOTDPMN15=TOTDPMN15+DPMNLAB0(IL)
	    TOTRPMN15=TOTRPMN15+RPMNLAB0(IL)
	    TOTBION15=TOTBION15+BNLAB0(IL)
	    TOTHUMN15=TOTHUMN15+HNLAB0(IL)
100     CONTINUE
        DO 200 IL=1,ICULT
	    DPMCARB0(IL)=TOTDPMC/ICULT
	    RPMCARB0(IL)=TOTRPMC/ICULT
	    BCARB0(IL)=TOTBIOC/ICULT
	    HCARB0(IL)=TOTHUMC/ICULT
	    DPMNIT0(IL)=TOTDPMN/ICULT
	    RPMNIT0(IL)=TOTRPMN/ICULT
	    BNIT0(IL)=TOTBION/ICULT
	    HNIT0(IL)=TOTHUMN/ICULT
	    DPMNLAB0(IL)=TOTDPMN15/ICULT
	    RPMNLAB0(IL)=TOTRPMN15/ICULT
	    BNLAB0(IL)=TOTBION15/ICULT
	    HNLAB0(IL)=TOTHUMN15/ICULT
200     CONTINUE
      ENDIF
C
C For each depth down to bottom of cultivation
C
      DO 300 IL=1,ICULT
C
C ...work out the c:N ratios of the pools
C
        CTONDPM=DPMCARB0(IL)/DPMNIT0(IL)
        CTONRPM=RPMCARB0(IL)/RPMNIT0(IL)
        CTONBIO=BCARB0(IL)/BNIT0(IL)
        CTONHUM=HCARB0(IL)/HNIT0(IL)        
C
C ...work out total carbon and depth of this layer
C
	  SOCNT=DPMCARB0(IL)+RPMCARB0(IL)+BCARB0(IL)+HCARB0(IL)
	  DEPTH=(IL*MAXDEPTH/MAXLAYER1)
C
C ...calculate amount of humus passed back to plant material: 
C    according to West & Post, 2002. Soil Sci Am J, 66, 1930-1946
C    SOC(till) - SOC(no till) in 0-7cm depth 
C      = 10 x (((SOC(no till)/10)-255.12)/1.2) 
C
        IF(DEPTH.LE.7)THEN
	    HUMTOPM=10*((SOCNT/10)-(((SOCNT/10)-255.12)/1.2))
C
C    SOC(till) - SOC(no till) in 0-7cm depth 
C      = 10 x (((SOC(no till)/10)-255.12)/1.2) 
C
	  ELSEIF(DEPTH.GT.7)THEN
	    HUMTOPM=10*((SOCNT/10)-(((SOCNT/10)-181.36)/0.93))
	  ENDIF
C
C ...Set the proportion of humus going to different pools
C
	  PROPB=0
	  CALL SET_DPMRPMRATIO(LU1,DRRAT)
	  PROPD=DRRAT/(1+DRRAT)
	  PROPR=1/(1+DRRAT)

C
C ... Calculate the proportion of humus passed to plant material
C
        PASS=HUMTOPM/HCARB0(IL)
	  IF(PASS.GT.1)PASS=1
C
C ... Pass C and N to appropriate pools
C
        IF(VIGOUR.GT.0.AND.VIGOUR.LT.1)PASS=VIGOUR

	  BCARB0(IL)=BCARB0(IL)+(PROPB*PASS*HCARB0(IL))
	  N15TON=BNLAB0(IL)/BNIT0(IL)
	  BNIT0(IL)=BNIT0(IL)+(PROPB*PASS*HCARB0(IL)/CTONBIO)
	  BNLAB0(IL)=BNLAB0(IL)+(N15TON*PROPB*PASS*HCARB0(IL)/CTONBIO)
	  HNIT0(IL)=HNIT0(IL)-(PROPB*PASS*HCARB0(IL)/CTONBIO)
	  HNLAB0(IL)=HNLAB0(IL)-(N15TON*PROPB*PASS*HCARB0(IL)/CTONBIO)

	  DPMCARB0(IL)=DPMCARB0(IL)+(PROPD*PASS*HCARB0(IL))
	  N15TON=DPMNLAB0(IL)/DPMNIT0(IL)
	  DPMNIT0(IL)=DPMNIT0(IL)+(PROPD*PASS*HCARB0(IL)/CTONHUM)
	  DPMNLAB0(IL)=DPMNLAB0(IL)+(N15TON*PROPD*PASS*HCARB0(IL)/CTONHUM)
	  HNIT0(IL)=HNIT0(IL)-(PROPD*PASS*HCARB0(IL)/CTONHUM)
	  HNLAB0(IL)=HNLAB0(IL)-(N15TON*PROPD*PASS*HCARB0(IL)/CTONHUM)

	  RPMCARB0(IL)=RPMCARB0(IL)+(PROPR*PASS*HCARB0(IL))
	  N15TON=RPMNLAB0(IL)/RPMNIT0(IL)
	  RPMNIT0(IL)=RPMNIT0(IL)+(PROPR*PASS*HCARB0(IL)/CTONHUM)
	  RPMNLAB0(IL)=RPMNLAB0(IL)+(N15TON*PROPR*PASS*HCARB0(IL)/CTONHUM)
	  HNIT0(IL)=HNIT0(IL)-(PROPR*PASS*RPMCARB0(IL)/CTONRPM)
	  HNLAB0(IL)=HNLAB0(IL)-(N15TON*PROPR*PASS*RPMCARB0(IL)/CTONRPM)

	  HCARB0(IL)=HCARB0(IL)-(PASS*HCARB0(IL))
300   CONTINUE
C
C Leave cultiv
C
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE DENITRIF3_NEMIS(DFACT,DENIT,DN15,RAIN,CH4,
     &                     WMAX,CO2,SOILW,SOILN,CRIT,NSOIL,TIN,SOIL15,
     &                     GDN2,GDN2O,G15DN2,G15DN2O,DENITRIFN,BULKDENS)
C
C Subroutine to calculate denitrification
C
      IMPLICIT NONE
C
C Variable local to this subroutine
C
	INTEGER MAXSOIL				! Max.no.of soil types
      PARAMETER(MAXSOIL=50)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	INTEGER MAXDEPTH			! Maximum depth of the soil profile
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
      REAL PN2FC					! Proportion of N2 produced at field capacity
	REAL PN2ZERO				! Proportion of N2 produced at wilting point
	DATA PN2FC /0.55/
	DATA PN2ZERO /10/
      REAL PARTH2O
	REAL FTEMP
	REAL PARTNO3
	REAL FWATEFF
	REAL FNITEFF
	REAL SCRIT
 	REAL DEFI(MAXLAYER)			! Calculate the amount of water before layer reaches field capacity
	REAL XRAIN
	REAL DPROP
	REAL TRAT
	REAL DL15
	REAL DLAYER
	INTEGER IL
	INTEGER NVERSION
	INTEGER NBRADBURY
	INTEGER NEMIS
	REAL CD1					! Amount of nitrate at which denitrification is 50% of its full potential (kg N / ha)
	DATA CD1 /16.5/
	REAL CD2					! Fitted constant describing water modifier for denitrification
	DATA CD2 /0.62/
	REAL CD3					! Fitted constant describing water modifier for denitrification
	DATA CD3 /0.38/
	REAL CD4					! Fitted constant describing water modifier for denitrification
	DATA CD4 /1.74/
C
C Variables passed to/from calling subroutine
C Soil factors
C
	REAL DFACT(MAXLAYER)		! IN/OUT: Denitrification factor
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL CO2(MAXLAYER)			! IN:CO2 emitted (kgC/ha/layer)
	REAL CH4(MAXLAYER)			! IN:CH4 emitted (kgC/ha/layer)
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL TIN(MAXLAYER)			! IN:15N/N in the nitrate in the layer
	REAL CRIT(MAXSOIL,MAXLAYER)	! IN/OUT: Minimum nitrate level (Sozanska)
	REAL DENIT					! IN:N lost by denitrification (kgN/ha)
	REAL DENITRIFN(MAXLAYER)	! OUT: N lost by denitrification (kgN/ha/layer
	REAL DN15					! IN:N15 lost by denitrification (kgN15/ha)
	REAL RAIN					! IN: Rainfall (mm/timestep)
	INTEGER NSOIL				! IN:Soil code number
	REAL GDN2					! IN:N2 lost by denitrification (kgN/ha)
	REAL GDN2O					! IN:N2O lost by denitrification (kgN/ha)
	REAL G15DN2					! IN:15N2 lost by denitrification (kgN15/ha)
	REAL G15DN2O				! IN:15N2O lost by denitrification (kgN15/ha)	
	REAL FNIT
	REAL FCO2
	REAL FNO3				
	REAL RATION2N2O
      REAL BULKDENS(MAXLAYER)
C
C Initialize Variables
C DENIT=Total N denitrified from top 50cm this week
C DN15=Total N15 denitrified from top 50cm this week
C XRAIN=rainfall this week
C
      DENIT=0
      DN15=0
      XRAIN=RAIN
      GDN2=0
      GDN2O=0
      G15DN2=0
      G15DN2O=0
C
C For each layer down the profile...
C
c      DO 20 IL=1,MAXLAYER1

      DO 20 IL=1,5
C
C Calculate the amount of water before layer reaches field capacity = DEFI()
C
      
        IF(SOILW(IL).LT.WMAX(IL))THEN
          DEFI(IL)=WMAX(IL)-SOILW(IL)
        ELSE
          DEFI(IL)=0.0
        END IF
C
C NEMIS DENITRIFICATION MODEL
C
        IF(SOILW(IL).LT.WMAX(IL))THEN
	    FWATEFF=(((SOILW(IL)/WMAX(IL)) - CD2)/CD3)
	  ELSE
	    FWATEFF=((1 - CD2)/CD3)
        ENDIF
  	  IF(FWATEFF.LT.0)FWATEFF=0
	  FWATEFF=FWATEFF**CD4
c	  FNITEFF=SOILN(IL)/(CD1 + SOILN(IL))
        FNITEFF=1/(CD1 + SOILN(IL))	
	  IF(FNITEFF.LT.0)FNITEFF=0						
	  DLAYER=(CH4(IL)+CO2(IL))*DFACT(IL)*SOILN(IL)*FWATEFF*FNITEFF		
C
C Pass excess water down to the next layer
C
        XRAIN=XRAIN-DEFI(IL)
        IF(XRAIN.LT.0)XRAIN=0
	  SCRIT=SOILN(IL)-CRIT(NSOIL,IL)*MAXDEPTH/(50.*MAXLAYER1)
        IF(DLAYER.GT.SCRIT)THEN
          DLAYER=SCRIT
        END IF
C
C Calculate the N15 denitrified
C
        TRAT=SOIL15(IL)/SOILN(IL)
        DL15=DLAYER*TRAT
        DENIT=DENIT+DLAYER
	  DENITRIFN(IL)=DLAYER
        DN15=DN15+DL15
        SOILN(IL)=SOILN(IL)-DLAYER
        SOIL15(IL)=SOIL15(IL)-DL15
        IF(SOILN(IL).GT.0)THEN
          TIN(IL)=SOIL15(IL)/SOILN(IL)
        ELSE
          TIN(IL)=0
        END IF
C Partition denitrification losses into N2 and N20 ##needs calibrating##
C
        IF(DLAYER.GT.0)THEN
	      IF(SOILW(IL).LT.WMAX(IL))THEN
 	        PARTH2O=PN2FC*SOILW(IL)/WMAX(IL)
	      ELSE
 	        PARTH2O=PN2FC
	      ENDIF	
		  PARTNO3=(1-SOILN(IL)/(200+SOILN(IL)))
	  IF (SOILN(IL).LT.(DLAYER*PN2ZERO))THEN
          GDN2=GDN2+0
		  GDN2O=GDN2O+DLAYER 
      ELSE
           GDN2=GDN2+(PARTH2O*PARTNO3*DLAYER)
	      GDN2O=GDN2O+(DLAYER-(PARTH2O*PARTNO3*DLAYER))
	  ENDIF 
		  FTEMP=PARTNO3*PARTH2O*DL15
	      G15DN2=G15DN2+FTEMP
	      G15DN2O=G15DN2O+DL15-FTEMP
	  ENDIF		
	  
	  					
C     Partioning of N2 and N2O, equations taken from Parton et al 1996 needs diffusion however
c      FNIT=1.4/(13**(17/13**(2.2*(SOILW(IL)/WMAX(IL)))))
c
c	FNO3=(1-(0.5+
c     &(1*ATAN(3.14*0.01*(10*SOILN(IL)/1.1-190))/3.14)))*25
c
c	FCO2=13+30.78*ATAN(3.14*0.07*((CH4(IL)+CO2(IL))-13))/3.14
c      IF (FNO3.GT.FCO2)THEN
c      RATION2N2O=FNIT*FCO2
c	ELSE
c     RATION2N2O=FNIT*FNO3
c	ENDIF
c
c	GDN2O=GDN2O+DLAYER/(1+RATION2N2O)
c	GDN2=GDN2+DLAYER/(1+1/RATION2N2O)
C
C Go back to repeat the calculation for the next layer down
C
  20  CONTINUE
C
C Leave DENITRIF3
C
      RETURN
      END
C
C--------------------------------------------------------------
C
      SUBROUTINE DIFFUSION(GAS,GASATM,LAYWIDTH,CLARRAY,BULKDENS,TOC,
     &                     SOILW,WSAT,AIRDIFF,WATDIFF,SECONDS,GASOUT)
C
C Subroutine to calculate gas diffusion between soil layers
C
c *** Note to Matt (delete after reading)- implicit none needed to check all variables are declared***
      IMPLICIT NONE
C
C Declare variables and arrays 
C (**Note to Matt (delete after reading), separate into variable passed in, out and local variables as done in other subroutines***)
C
      INTEGER IL						! Local layer counter
	INTEGER MAXLAY					! No.of layers in the soil profile
	INTEGER MAXLAYER
	PARAMETER (MAXLAYER=60)
	DATA MAXLAY /60/
	REAL LAYWIDTH					! Layer width (m)
	REAL AIRDIFF					! Coefficient of diffusion in air
	REAL WATDIFF					! Coefficient of diffusion in water
	REAL GAS(MAXLAYER)				! Gas density (kg ha-1 layer-1)
	REAL GASNEW(MAXLAYER)			! Used in calculation of gas release to atmos
	REAL GASATM						! Gas concentration in atmosphere
									! (kg ha-1 in 5cm layer)
	REAL CLARRAY(MAXLAYER)			! Clay content of layer (%)
	REAL BULKDENS(MAXLAYER)			! Bulk density of soil (g cm-3)
	REAL TOC(MAXLAYER)				! Total carbon (kg ha-1 layer-1)
	REAL SOILW(MAXLAYER)			! Available water (mm/layer)
	REAL WSAT(MAXLAYER)				! Available water at saturation (mm/layer)
	REAL GASGRAD(MAXLAYER)			! Gradient between layers
	REAL TORTUOUS(MAXLAYER)			! Soil tortuosity variable
	REAL PORE(MAXLAYER)				! Porosity
	REAL POREAIR(MAXLAYER)			! Air-filled porosity
	REAL POREWAT(MAXLAYER)			! Water-filled porosity
	REAL GASDA(MAXLAYER)			! Rate of diffusion in air
	REAL GASDW(MAXLAYER)			! Rate of diffusion in water
	REAL GASMOVE(MAXLAYER)			! Overall diffusion rate
	REAL SECONDS					! Number of seconds in a timestep
	REAL GASOUT						! Diffusion to or from atmosphere
	REAL GASOUTSTEP                 ! Timestep diffusion to or from atmosphere
	REAL TSTEP                      ! Time step for iterating diffusion
	REAL TSTEPS                     ! Number of time steps for iterating diffusion
      REAL DIFFSTEP                   ! Loop to allow diffusion up and down
C
C Determine gas gradients between layers & gradient to atmosphere
C (Marshall & Holmes, Soil Physics, any edition)

C
C Iteration to prevent crashing with large time step
C
	TSTEPS=SECONDS/3600.0
	GASOUT=0.0
C
C
	DO 500 TSTEP=1,TSTEPS
C
	DO 600 DIFFSTEP=1,2
C
	DO 800 IL=1,MAXLAY
        GASNEW(IL)=GAS(IL)
800   CONTINUE
C
	DO 300 IL=1,MAXLAY
	  if(DIFFSTEP.EQ.1)then
          if(IL.EQ.1)then
            GASGRAD(IL)=(GAS(IL)-GASATM)/LAYWIDTH
  	    else
            GASGRAD(IL)=(GAS(IL)-GAS(IL-1))/LAYWIDTH
	    endif
	  else
          if(IL.EQ.MAXLAY)then
            GASGRAD(IL)=0.
  	    else
            GASGRAD(IL)=(GAS(IL)-GAS(IL+1))/LAYWIDTH
	    endif
	  endif
300   CONTINUE
C
	DO 400 IL=1,MAXLAY
C Determine soil topology value
        TORTUOUS(IL)=1.5+(0.015*CLARRAY(IL))
C Calculate porosity
	  PORE(IL)=0.4+(0.4*TOC(IL)/55000.)+(CLARRAY(1)/1000.)
	  if(PORE(IL).GT.1.)then
	    PORE(IL)=1.0
        endif
C Calculate air-filled porosity
        if(WSAT(IL).GT.0.)then
	    POREWAT(IL)=PORE(IL)*(SOILW(IL)/WSAT(IL))
	  else
          POREWAT(IL)=0.0
	  endif
C Calculate water-filled porosity
        POREAIR(IL)=PORE(IL)-POREWAT(IL)
C Calculate diffusion rate in air
	  GASDA(IL)=AIRDIFF*(POREAIR(IL)**(TORTUOUS(IL)+1.))
C Calculate diffusion rate in water
        GASDW(IL)=WATDIFF*(POREWAT(IL)**(TORTUOUS(IL)+1.))
C Eventual proportional move equals proportional move raised to the power of 3600
        GASMOVE(IL)=1.-((GASDA(IL)+GASDW(IL))/LAYWIDTH)
        GASMOVE(IL)=GASMOVE(IL)**3600.
C Calculate size of move to be made
	  GASMOVE(IL)=GASMOVE(IL)*GASGRAD(IL)*LAYWIDTH*0.5
        if(DIFFSTEP.EQ.1)then
	    GASNEW(IL)=GASNEW(IL)-GASMOVE(IL)
	    if(IL.EQ.1)then
	      GASOUTSTEP=GASMOVE(IL)
	    else
	      GASNEW(IL-1)=GASNEW(IL-1)+GASMOVE(IL)
	    endif
	  else
	    GASNEW(IL)=GASNEW(IL)-GASMOVE(IL)
	    if(IL.LT.MAXLAY)then
	      GASNEW(IL+1)=GASNEW(IL+1)+GASMOVE(IL)
	    endif
	  endif
C
400   CONTINUE
C
	DO 700 IL=1,MAXLAY
        GAS(IL)=GASNEW(IL)
	  if(GAS(IL).LT.0)then
          GAS(IL)=0.0
	  endif
700   CONTINUE
C
      GASOUT=GASOUT+GASOUTSTEP
C
600   CONTINUE
C
500   CONTINUE

	END
C
C--------------------------------------------------------------
C

      SUBROUTINE METHANE_AITKENHEAD(SOILTEMP,NSOIL,
     &                        BCARB0,HCARB0,DPMCARB0,RPMCARB0,
     &                        BNIT0,HNIT0,DPMNIT0,RPMNIT0,
     &                        BNLAB0,HNLAB0,DPMNLAB0,RPMNLAB0,
     &                        DPMRATE,RPMRATE,BRATE,HRATE,
     &                        HY,SOILW,WMAX,WSAT,
     &                        ALPHA,BETA,GAMMA,DELTA,
     &                        SOILN,SOIL15,AMMN,AMMN15,CRIT,
     &                        PHARRAY,CO2,CH4,CLARRAY,BULKDENS,TOC,
     &                        SECONDS,CH4TOAIR)		
C   
C Subroutine to calculate the amount of methane emissions
C
      IMPLICIT NONE
C
C Declare local  variables
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	INTEGER MAXDEPTH			! Maximum depth of the soil profile
	INTEGER MAXSOIL				! Max.no.of soil types
	PARAMETER(MAXSOIL=50)
	INTEGER IDEPTH				! Depth of this layer
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
	REAL CH4INAIR				! Concentration of methane in air at sea level
								! (kg/ha/5cm layer)
	DATA CH4INAIR /0.0011/
	REAL ANDECDPM				! Anaerobic decomposition of DPM (kg C / ha)
	REAL ANDECRPM				! Anaerobic decomposition of RPM (kg C / ha)
	REAL ANDECBIO				! Anaerobic decomposition of BIO (kg C / ha)
	REAL ANDECHUM				! Anaerobic decomposition of HUM (kg C / ha)
	REAL BIOCPROD				! Biomass C produced by anaer.decomp. 
								! (kg C / ha)
	REAL HUMCPROD				! Humus C produced by anaer.decomp.
								! (kg C / ha)
	REAL BIONPROD				! Biomass C produced by anaer.decomp. 
								! (kg C / ha)
	REAL HUMNPROD				! Humus C produced by anaer.decomp.
								! (kg C / ha)
	REAL CH4DPM					! CH4 C produced by anaer.decomp. from DPM
								! (kg C / ha)
	REAL CH4RPM					! CH4 C produced by anaer.decomp. from RPM
								! (kg C / ha)
	REAL CH4BIO					! CH4 C produced by anaer.decomp. from BIO
								! (kg C / ha)
	REAL CH4HUM					! CH4 C produced by anaer.decomp. from HUM
								! (kg C / ha)
      REAL VFACT					! Soil dependent factor accounting for 
								!  diffusion and oxidation rates
      REAL TFACT					! transport factor 
								! (non-transporters t=0; transporting non-sedges t=0.25; 
								!  sedges t=1)
	REAL CH4OXID				! CH4 oxidised
      REAL NTOCDPM				! N:C ratio of DPM
      REAL NTOCRPM				! N:C ratio of RPM
      REAL NTOCBIO				! N:C ratio of BIO
      REAL NTOCHUM				! N:C ratio of HUM
	REAL BIONNEED				! N needed to maintain stable C:N ratio of BIO
	REAL HUMNNEED				! N needed to maintain stable C:N ratio of HUM
      REAL AVAILAMMN				! Ammonium N available kgN/ha/layer
      REAL AVAILSOILN				! Nitrate N available kgN/ha/layer
	REAL BIONIMM				! N immobilised by biomass pool kgN/ha/layer
	REAL HUMNIMM				! N immobilised by humus pool kgN/ha/layer
	INTEGER IL					! Layer counter
	REAL BIOCPROD_DPM			! Biomass C produced by anaer.decomp.of DPM 
	REAL BIONPROD_DPM			! Biomass N produced by anaer.decomp.of DPM 
	REAL BIONNEED_DPM			! N needed by BIO for anaer.decomp.of DPM
	REAL HUMCPROD_DPM			! Humus C produced by anaer.decomp.of DPM 
	REAL HUMNPROD_DPM			! Humus N produced by anaer.decomp.of DPM 
	REAL HUMNNEED_DPM			! N needed by HUM for anaer.decomp.of DPM
	REAL BIOCPROD_RPM			! Biomass C produced by anaer.decomp.of RPM 
	REAL BIONPROD_RPM			! Biomass N produced by anaer.decomp.of RPM 
	REAL BIONNEED_RPM			! N needed by BIO for anaer.decomp.of RPM
	REAL HUMCPROD_RPM			! Humus C produced by anaer.decomp.of RPM 
	REAL HUMNPROD_RPM			! Humus N produced by anaer.decomp.of RPM 
	REAL HUMNNEED_RPM			! N needed by HUM for anaer.decomp.of RPM
	REAL FALLOON_FACTOR			! Ratio of anaerobic to aerobic decomposition 
	REAL SECONDS                ! The number of seconds in a timestep

C>>> Temp change JUS 28/05/08 >>>>
      DATA FALLOON_FACTOR /0.3/	! Falloon_factor - adjusted to get
c     DATA FALLOON_FACTOR /0.165/	! Falloon_factor - adjusted to get
C<<< Temp change JUS 28/05/08 <<<<
	REAL METHADIF               ! Diffusion of methane in air
	REAL METHWDIF               ! Diffusion of methane in water
	DATA METHADIF /0.00002/
	DATA METHWDIF /0.0000000025/
	REAL LAYWIDTH               ! Width of soil layer
	DATA LAYWIDTH /0.05/
	REAL OXMAX					! Maximum proportional methane oxidation rate
	DATA OXMAX /0.05/
C
C Declare variables passed to/from this subroutine
C ...Weather descriptors
C
	REAL SOILTEMP(MAXLAYER)		! IN: Soil temperature (deg.C/timestep)
C
C ...Soil descriptors
C
	INTEGER NSOIL				! IN:Soil code number
  	REAL BCARB0(MAXLAYER)		! IN:C in soil biomass (kgC/ha/layer)
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
      REAL BNIT0(MAXLAYER)		! IN:N in soil biomass (kgN/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN:N in soil humus (kgN/ha/layer)
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN:N15 in soil biomass (kgN15/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN:N15 in soil humus (kgN15/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL DPMRATE(MAXLAYER)		! IN/OUT: Rate constant for DPM decompn
	REAL RPMRATE(MAXLAYER)		! IN/OUT: Rate constant for RPM decompn
	REAL BRATE(MAXLAYER)		! IN/OUT: Rate constant for HUM decompn
	REAL HRATE(MAXLAYER)		! IN/OUT: Rate constant for BIO decompn
	REAL HY(MAXSOIL,MAXLAYER)	! IN/OUT:Stable N:C ratio of BIO & HUM pools
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
	REAL ALPHA(MAXLAYER)		! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL BETA(MAXLAYER)			! IN/OUT: Prop.HUM produced on BIO decompn
	REAL GAMMA(MAXLAYER)		! IN/OUT: Prop.BIO produced on HUM decompn 
	REAL DELTA(MAXLAYER)		! IN/OUT: Prop.HUM produced on HUM decompn 
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
	REAL CRIT(MAXSOIL,MAXLAYER)		! IN:Minimum nitrate level (kgN/ha/50cm)
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN:pH of soil in this layer
      REAL PHRATE					! IN:pH rate modifier (proportion)
	REAL TRATE					! IN:temperature rate modifying factor (prop.)
	REAL WRATE					! IN:Moisture rate modifyer for (prop)
	REAL CO2(MAXLAYER)			! IN:CO2 emitted (kgC/ha/layer)
	REAL CH4(MAXLAYER)			! IN:CH4 produced (kgC/ha/layer)
	REAL CH4GEN(MAXLAYER)	    ! IN:CH4 generated in soil (kgC/ha/layer)
	REAL CH4GRAD(MAXLAYER)      ! IN:CH4 gradient to neighbouring layers
	                            ! (kgC/ha/layer/metre)
      REAL TORTUOUS(MAXLAYER)     ! IN:Topology variable controlling diffusion
	REAL PORE(MAXLAYER)         ! IN:Total porosity
	REAL POREAIR(MAXLAYER)      ! IN:Air-filled porosity
	REAL POREWAT(MAXLAYER)      ! IN:Water-filled porosity
	REAL CH4DA(MAXLAYER)        ! IN:Diffusion rate in air
	REAL CH4DW(MAXLAYER)        ! IN:Diffusion rate in water
	REAL CH4MOVE(MAXLAYER)      ! IN:Diffusion rate between layers
      REAL CLARRAY(MAXLAYER)		! IN:Clay content of the layer (%)
	REAL BULKDENS(MAXLAYER)		! IN:Bulk density of the layer (g/cm3)
	REAL TOC(MAXLAYER)			! IN:Total organic C (kgC/ha/layer)
      REAL CH4TOAIR				! OUT:CH4 released to atmosphere
C
C Calculate anaerobic decomposition
C Dan = C(pool) x exp(k(pool) x m(water) x m(temp) x m(pH)
C 1. Rate constants assumed to be the same as for aerobic decomposition
C    k(dpm)=10yr-1; k(rpm)=0.3yr-1; k(bio)=0.66yr-1; k(hum)=0.02yr-1 reduced by falloon_factor
C 2. m(water)=0 at awc/wsat=0 (awc = available water in 5cm layer)
C                             (wsat = available water at saturation in 5cm layer)
C             increases to 1 at awc/wsat=0.2, decreasing again to 0 at awc/wsat=0.6
C 3. m(temp) = T/25      (T = temp of 5cm layer)
C 4. m(pH) = 1-((pH-4)/3)
C                        (pH = pH of the 5cm layer)
C 5. m(crop) assumed to be 1
C
      DO 100 IL=1,MAXLAYER1
	  CH4(IL)=0
        CALL MODFACTS_ANDEC(SOILW(IL),WMAX(IL),WSAT(IL),WRATE,
     &                SOILTEMP(IL),TRATE,
     &			    PHARRAY(NSOIL,IL),PHRATE)
	  ANDECDPM=(1-EXP(-FALLOON_FACTOR*DPMRATE(IL)*WRATE*TRATE*PHRATE))
	  ANDECDPM=ANDECDPM*DPMCARB0(IL)
	  ANDECRPM=(1-EXP(-FALLOON_FACTOR*RPMRATE(IL)*WRATE*TRATE*PHRATE))
	  ANDECRPM=ANDECRPM*RPMCARB0(IL)
	  ANDECBIO=(1-EXP(-FALLOON_FACTOR*BRATE(IL)*WRATE*TRATE*PHRATE))
	  ANDECBIO=ANDECBIO*BCARB0(IL)
	  ANDECHUM=(1-EXP(-FALLOON_FACTOR*HRATE(IL)*WRATE*TRATE*PHRATE))
	  ANDECHUM=ANDECHUM*HCARB0(IL)
C
C Proportion anaerobic decomposition into biomass, humus and methane
C 1. BIOprod = alpha x Dan (alpha assumed to be the same as for aerobic decomposition)
C 2. HUMprod = beta x Dan (beta assumed to be the same as for aerobic decomposition)
C 3. CH4prod = (1-alpha-beta) x Dan
C
        BIOCPROD=ALPHA(IL)*(ANDECDPM+ANDECRPM+ANDECBIO)
	  BIOCPROD=BIOCPROD+(GAMMA(IL)*ANDECHUM)
	  HUMCPROD=BETA(IL)*(ANDECDPM+ANDECRPM+ANDECBIO)
	  HUMCPROD=HUMCPROD+(DELTA(IL)*ANDECHUM)
C
C Work out how much N is needed to maintain C:N of BIO and HUM pools
C
        IF(DPMCARB0(IL).LE.0)THEN
	    NTOCDPM=0
	  ELSE 
	    NTOCDPM=DPMNIT0(IL)/DPMCARB0(IL)
	  ENDIF
        IF(RPMCARB0(IL).LE.0)THEN
	    NTOCRPM=0
	  ELSE 
          NTOCRPM=RPMNIT0(IL)/RPMCARB0(IL)
	  ENDIF
        IF(BCARB0(IL).LE.0)THEN
	    NTOCBIO=0
	  ELSE 
          NTOCBIO=BNIT0(IL)/BCARB0(IL)
	  ENDIF
        IF(HCARB0(IL).LE.0)THEN
	    NTOCHUM=0
	  ELSE 
          NTOCHUM=HNIT0(IL)/HCARB0(IL)
	  ENDIF
        BIONPROD=ALPHA(IL)*ANDECDPM*NTOCDPM
        BIONPROD=BIONPROD+(ALPHA(IL)*ANDECRPM*NTOCRPM)
        BIONPROD=BIONPROD+(ALPHA(IL)*ANDECBIO*NTOCBIO)
        BIONPROD=BIONPROD+(GAMMA(IL)*ANDECHUM*NTOCHUM)
        HUMNPROD=BETA(IL)*ANDECDPM*NTOCDPM
        HUMNPROD=HUMNPROD+(BETA(IL)*ANDECRPM*NTOCRPM)
        HUMNPROD=HUMNPROD+(BETA(IL)*ANDECBIO*NTOCBIO)
        HUMNPROD=HUMNPROD+(DELTA(IL)*ANDECHUM*NTOCHUM)
	  BIONNEED=HY(NSOIL,IL)*(BCARB0(IL)+BIOCPROD)
	  BIONNEED=BIONNEED-BNIT0(IL)-BIONPROD
	  HUMNNEED=HY(NSOIL,IL)*(HCARB0(IL)+HUMCPROD)
	  HUMNNEED=HUMNNEED-HNIT0(IL)-HUMNPROD
C
C If pools need to lose N, add to ammonium pool
C

        AVAILAMMN=AMMN(IL)-CRIT(NSOIL,IL)*MAXDEPTH/(50.*MAXLAYER1)
        AVAILSOILN=SOILN(IL)-CRIT(NSOIL,IL)*MAXDEPTH/(50.*MAXLAYER1)
        IF(BIONNEED+HUMNNEED.LT.0)THEN
	    AMMN(IL)=AMMN(IL)-BIONNEED-HUMNNEED
	    BNIT0(IL)=BNIT0(IL)+BIONNEED
	    HNIT0(IL)=HNIT0(IL)+HUMNNEED
C
C ...or if pools need to gain N, take it from the ammonium pool, 
C
	  ELSEIF(BIONNEED+HUMNNEED.LT.AVAILAMMN)THEN
	    AMMN(IL)=AMMN(IL)-BIONNEED-HUMNNEED
	    BNIT0(IL)=BNIT0(IL)+BIONNEED
	    HNIT0(IL)=HNIT0(IL)+HUMNNEED
C
C ...then when that runs out, take it from the nitrate pool as well,
C
	  ELSEIF(BIONNEED+HUMNNEED.LT.AVAILAMMN+AVAILSOILN)THEN
	    AMMN(IL)=CRIT(NSOIL,IL)*MAXDEPTH/(50.*MAXLAYER1)
	    IF(BIONNEED+HUMNNEED.GT.0)THEN
	      BIONIMM=AVAILAMMN*BIONNEED/(BIONNEED+HUMNNEED)
	      HUMNIMM=AVAILAMMN*HUMNNEED/(BIONNEED+HUMNNEED)
	    ELSE
	      BIONIMM=0
	      HUMNIMM=0
          ENDIF
	    BIONNEED=BIONNEED-BIONIMM
	    HUMNNEED=HUMNNEED-HUMNIMM
	    SOILN(IL)=SOILN(IL)-BIONNEED-HUMNNEED
	    BIONIMM=BIONIMM+BIONNEED
	    HUMNIMM=HUMNIMM+HUMNNEED
	    BNIT0(IL)=BNIT0(IL)+BIONIMM
	    HNIT0(IL)=HNIT0(IL)+HUMNIMM
C
C ...or if there is not enough N, only decompose the fraction that is possible
C
	  ELSEIF(BIONNEED+HUMNNEED.GT.AVAILAMMN+AVAILSOILN)THEN
	    AMMN(IL)=CRIT(NSOIL,IL)*MAXDEPTH/(50.*MAXLAYER1)
	    SOILN(IL)=CRIT(NSOIL,IL)*MAXDEPTH/(50.*MAXLAYER1)
	    IF(BIONNEED+HUMNNEED.GT.0)THEN
	      BIONIMM=(AVAILAMMN+AVAILSOILN)*BIONNEED/(BIONNEED+HUMNNEED)
	      HUMNIMM=(AVAILAMMN+AVAILSOILN)*HUMNNEED/(BIONNEED+HUMNNEED)
	    ELSE
	      BIONIMM=0
	      HUMNIMM=0
          ENDIF
	    BIOCPROD_DPM=ALPHA(IL)*ANDECDPM
	    BIONPROD_DPM=NTOCDPM*BIOCPROD_DPM
	    BIONNEED_DPM=(HY(NSOIL,IL)*BIOCPROD_DPM)-BIONPROD_DPM
	    HUMCPROD_DPM=BETA(IL)*ANDECDPM
	    HUMNPROD_DPM=NTOCDPM*HUMCPROD_DPM
	    HUMNNEED_DPM=(HY(NSOIL,IL)*HUMCPROD_DPM)-HUMNPROD_DPM
	    BIOCPROD_RPM=ALPHA(IL)*ANDECRPM
	    BIONPROD_RPM=NTOCRPM*BIOCPROD_RPM
	    BIONNEED_RPM=(HY(NSOIL,IL)*BIOCPROD_RPM)-BIONPROD_RPM
	    HUMCPROD_RPM=BETA(IL)*ANDECRPM
	    HUMNPROD_RPM=NTOCRPM*HUMCPROD_RPM
	    HUMNNEED_RPM=(HY(NSOIL,IL)*HUMCPROD_RPM)-HUMNPROD_RPM
	    IF(BIONNEED_DPM.GT.0.OR.HUMNNEED_DPM.GT.0) 
     &	  ANDECDPM=ANDECDPM*(BIONIMM+HUMNIMM)/(BIONNEED+HUMNNEED)
	    IF(BIONNEED_RPM.GT.0.OR.HUMNNEED_RPM.GT.0) 
     &	  ANDECRPM=ANDECRPM*(BIONIMM+HUMNIMM)/(BIONNEED+HUMNNEED)
	    BIONNEED=BIONNEED-BIONIMM
	    HUMNNEED=HUMNNEED-HUMNIMM
	    BNIT0(IL)=BNIT0(IL)+BIONIMM
	    HNIT0(IL)=HNIT0(IL)+HUMNIMM
        ENDIF
C
C Adjust pools to account for C lost by anaerobic decomposition
C
        DPMCARB0(IL)=DPMCARB0(IL)-ANDECDPM
	  RPMCARB0(IL)=RPMCARB0(IL)-ANDECRPM
        BCARB0(IL)=BCARB0(IL)-ANDECBIO+BIONPROD
	  HCARB0(IL)=HCARB0(IL)-ANDECHUM+HUMNPROD
	  CH4DPM=(1-ALPHA(IL)-BETA(IL))*ANDECDPM
	  CH4RPM=(1-ALPHA(IL)-BETA(IL))*ANDECRPM
	  CH4BIO=(1-ALPHA(IL)-BETA(IL))*ANDECBIO
	  CH4HUM=(1-GAMMA(IL)-DELTA(IL))*ANDECHUM
	  CH4GEN(IL)=CH4DPM+CH4RPM+CH4BIO+CH4HUM
C
C Calculate the amount of CH4 that is oxidised to CO2
C 1. CH4oxid = CH4prod with temperature, moisture and pH variation
C
C Assumption that oxidation rate maximum equals some constant value 'OXMAX'
C (Teh, Silver & Conrad, 2005)
C (Tamai, Takenaka, Ishizuka, 2007)
C
C Temperate rate goes from 0 at 0 Centigrade to maximum of 1 at 25C
C (Blodau, Roulet, Heitmann, Stewart, Beer, Lafleur & Moore, 2007)
C 
C Moisture rate goes from 0 at 0 water to maximum at 0.2, then back to 0 at sat
C (Von Arnold, Weslien, Nilsson, Svensson & Klemedtsson, 2005)
C (Bradford, Wookey, Ineson & Lappin-Scott, 2001)
C
C PH rate goes from 0 at PH=2 to 1 at PH=7
C (Regina, Pihlatie, Esala, Alakukku, 2007)
C (Kravchenko & Yu, 2006)
C
        CH4OXID=OXMAX*CH4(IL)
	  if(SOILTEMP(IL).LT.0)then
	    CH4OXID=0
	  endif
	  if(SOILTEMP(IL).LT.25)then
	    CH4OXID=CH4OXID*(SOILTEMP(IL)/25.)
	  endif

	  if(SOILW(IL).GT.0)then
          if((SOILW(IL)/WSAT(IL)).GT.0.2)then
            if((SOILW(IL)/WSAT(IL)).GT.0.6)then
  	        CH4OXID=CH4OXID*2.5*(1.-(SOILW(IL)/WSAT(IL)))
	      else
  	        CH4OXID=CH4OXID*1.
	      endif
  	    else
            if((SOILW(IL)/WSAT(IL)).LT.0)then
              CH4OXID=0
	      else
  	        CH4OXID=CH4OXID*5.*(SOILW(IL)/WSAT(IL))
	      endif
          endif
	  else
          CH4OXID=0
	  endif

        if(PHARRAY(NSOIL,IL).LT.2)then
          CH4OXID=0
	  else
          if(PHARRAY(NSOIL,IL).GT.7)then
            CH4OXID=CH4OXID*1.
	    else
            CH4OXID=CH4OXID*(PHARRAY(NSOIL,IL)-2.0)/5.0
	    endif
	  endif
C
C Adjust pools to account for methane oxidation
C
	  CO2(IL)=CO2(IL)+CH4OXID
	  CH4(IL)=CH4(IL)-CH4OXID
100   CONTINUE
C
C Calculate addition to methane stored in soil
C
      DO 200 IL=1,MAXLAYER1
        CH4(IL)=CH4(IL)+CH4GEN(IL)
200   CONTINUE
C
C If timestep is less than weekly, calculate methane movement between soil layers
C
      IF(SECONDS.LT.7*24*60*60)THEN
        CALL DIFFUSION(CH4,CH4INAIR,LAYWIDTH,CLARRAY,BULKDENS,TOC,
     &               SOILW,WSAT,METHADIF,METHWDIF,SECONDS,CH4TOAIR)
	ENDIF
	END
C
C--------------------------------------------------------------
C
      SUBROUTINE GET_DECOMP_RATECONSTS(CONVER_F,NSOIL,BIORATE,HUMRATE,
     &                           DPMRATE,RPMRATE,BRATE,HRATE)
C
C Set Rate constants for decomposition of RO (=R), BIO (=B) and HUM (=H)
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXSOIL				! Max.no.of soil types
	PARAMETER(MAXSOIL=50)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER IL					! Layer counter
C
C Variables passed to/from this subroutine
C
	REAL BIORATE(MAXSOIL,MAXLAYER)	! IN: Rate constant for BIO decompn/yr
	REAL BRATE(MAXLAYER)		! OUT: Rate constant for HUM decompn
      REAL CONVER_F				! IN:Conversion to correct timestep
	REAL DPMRATE(MAXLAYER)		! OUT: Rate constant for DPM decompn
	REAL HRATE(MAXLAYER)		! OUT: Rate constant for BIO decompn
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN: Rate constant for HUM decompn/yr
	INTEGER NSOIL				! IN:Soil code number
	REAL RPMRATE(MAXLAYER)		! OUT: Rate constant for RPM decompn

C
C Set the residual N and C in the RO,BIO and HUM pools
C
      DO 100 IL=1,MAXLAYER1
        DPMRATE(IL)=10.0/52.0*CONVER_F
        RPMRATE(IL)=0.3/52*CONVER_F
        BRATE(IL)=BIORATE(NSOIL,IL)*CONVER_F
        HRATE(IL)=HUMRATE(NSOIL,IL)*CONVER_F
100   CONTINUE
C
C Leave GET_DECOMP_RATECONSTS
C
      END
C
C-------------------------------------------------------------
C
	SUBROUTINE GET_IOM_FROM_FALLOON_EQN(TOC_IL,IOM_IL)
C
C Subroutine to get IOM from the Falloon equation
C
      IMPLICIT NONE
C
C Variables passed to/from this routine
C
	REAL IOM_IL				! OUT:Inert organic C (kgC/ha/layer)
	REAL TOC_IL				! IN:Meas.total organic C (kgC/ha/layer)
C
C Calculate IOM from the Falloon equation
C      
	IOM_IL=1000*(0.049*((TOC_IL/1000)**1.139))
C
C Leave GET_IOM_FROM_FALLOON_EQN
C
      END
C
C--------------------------------------------------------------
C
      SUBROUTINE GETSTABLECN(PH,STABCN)
C
C Subroutine to get the stable C:N ratio from pH
C
      IMPLICIT NONE
C
C Variables internal to this routine
C
	REAL BACCN,FUNCN,BOTPH,BOTBAC,TOPPH,TOPBAC,FRACBAC,FRACFUN
C
C Variables passed to/from this routine
C
	REAL PH,STABCN
C
C Initialise variables (linear relationship from pH 4-5.5, bacteria = 20% to 50%)
C
      BACCN=5.5
	FUNCN=11.5
	BOTPH=4
	BOTBAC=0.2
	TOPPH=5.5
	TOPBAC=0.5
C
C Caculate the fraction of bacteria according to pH
C
	IF(PH.LE.BOTPH)THEN
	  FRACBAC=BOTBAC
	ELSEIF(PH.GT.BOTPH.AND.PH.LE.TOPPH)THEN
	  FRACBAC=BOTBAC+((TOPBAC-BOTBAC)*(PH-BOTPH)/(TOPPH-BOTPH))
	ELSEIF(PH.GT.TOPPH)THEN
	  FRACBAC=TOPBAC
	ENDIF
C
C Calculate the fraction of fungi by difference
C
	FRACFUN=1-FRACBAC
C
C Calculate the stable C:N ratio from the proportion of bacteria and fungi and the stable C:N ratio of each pool
C
	STABCN=(FRACBAC*BACCN)+(FRACFUN*FUNCN)
C
C Leave GETSTABLECN
C
	RETURN
	END
C
C------------------------------------------------------------
C
      SUBROUTINE LEACH_DOC(MOBDOC,MOBDON,LEACHDOC,LEACHDON,
     &                     DRAINW,WMAX,SOILW,REFIL,
     &                     PHARRAY,NSOIL)			
C
C Subroutine to calculate denitrification
C
      IMPLICIT NONE
C
C Variable local to this subroutine
C MAXLAYER = Maximum number of layers in the soil profile
C MAXLAYER1 = Maximum number of layers in the soil profile
C MAXDEPTH = Maximum depth of the soil profile
C IL = layer counter
C MOBPART1 = partitioning of mobile DOC into dissolved and sorbed
C DOCSORB() = Amount of DOC pool that is sorbed (kg C / ha / layer)
C DOCSOLV() = Amount of DOC pool that is in solution (kg C / ha / layer)
C DOCIN() = DOC input to the layer (kg C / ha / layer)
C DOCOUT() = DOC output from the layer (kg C / ha / layer)
C DOCSTOP = DOC that is leached into this layer and stops here (kg C / ha / layer)
C DONSORB() = Amount of DON pool that is sorbed (kg N / ha / layer)
C DONSOLV() = Amount of DON pool that is in solution (kg N / ha / layer)
C DONIN() = DON input to the layer (kg N / ha / layer)
C DONOUT() = DON output from the layer (kg N / ha / layer)
C DONSTOP = DON that is leached into this layer and stops here (kg N / ha / layer)
C
	INTEGER MAXLAYER,MAXDEPTH,MAXLAYER1,MAXSOIL
	PARAMETER (MAXLAYER=60)
	PARAMETER (MAXSOIL=50)
	DATA MAXLAYER1 /60/ 
	DATA MAXDEPTH /300/
	INTEGER IL
	REAL MOBPART1,PH
	REAL DOCIN(MAXLAYER)
	REAL DOCOUT(MAXLAYER)
      REAL DOCSOLV(MAXLAYER)
	REAL DOCSORB(MAXLAYER)
	REAL DOCSTOP
	REAL DONIN(MAXLAYER)
	REAL DONOUT(MAXLAYER)
      REAL DONSOLV(MAXLAYER)
	REAL DONSORB(MAXLAYER)
	REAL DONSTOP
	
C
C Variables passed to/from calling subroutine
C LEACHDOC() = DOC leached from each layer
C LEACHDON() = DON leached from each layer
C WMAX() = Maximum water content of this layer (mm / layer)
C SOILW() = Available water content of this layer (mm / layer)
C MOBDOC() = Total amount of mobile DOC pool (kg C / ha / layer)
C MOBDON() = Total amount of mobile DON pool (kg N / ha / layer)
C DRAINW() = Amount of water draining from this layer to the next (mm / layer)
C REFIL() = Amount of space left in layer for water (mm / layer)			
C PHARRAY()
C NSOIL()
C 
	REAL LEACHDOC(MAXLAYER),LEACHDON(MAXLAYER)
      REAL WMAX(MAXLAYER),SOILW(MAXLAYER),DRAINW(MAXLAYER)
	REAL REFIL(MAXLAYER),MOBDOC(MAXLAYER),MOBDON(MAXLAYER)
	REAL PHARRAY(MAXSOIL)
	INTEGER NSOIL
C
C Calculation for every layer down the profile...
C
      DO 100 IL=1,MAXLAYER1
C
C
        PH=PHARRAY(NSOIL)
        MOBPART1=(PH-3)/4
	  IF (MOBPART1.LT.0)THEN 
	     MOBPART1=0.0
	  ENDIF
	  IF (MOBPART1.GT.1)THEN
	     MOBPART1=1
	  ENDIF
        DOCSOLV(IL)=MOBDOC(IL)*MOBPART1
        DOCSORB(IL)=MOBDOC(IL)*(1-MOBPART1)
        DONSOLV(IL)=MOBDON(IL)*MOBPART1
        DONSORB(IL)=MOBDON(IL)*(1-MOBPART1)
C
C If not in the top layer pass leached DOC down from the layer above
C
        IF(IL.GT.1)THEN
	    DOCIN(IL)=DOCOUT(IL-1)
	    DONIN(IL)=DONOUT(IL-1)
C
C If in top layer do not pass from the layer above
C
        ELSE
          DOCIN(IL)=0.0
          DONIN(IL)=0.0
        END IF
C
C Calculate the amount of DOC leaching from the proportion of the total
C water moving to the next layer
C Assumptions:	1. Only non-sorbed DOC is leached
C				2. DOC is leached by piston flow according to concentration
C 1. Pick up DOC from this layer
C a. If all water in layer is drained out, leach all soluble DOC
C
        IF(DRAINW(IL).GT.SOILW(IL))THEN
	    DOCOUT(IL)=DOCSOLV(IL)
	    DOCSOLV(IL)=0
	    DONOUT(IL)=DONSOLV(IL)
	    DONSOLV(IL)=0
C
C b. If only a proportion of water in layer is drained out, 
C     leach the same proportion of soluble DOC.
C
        ELSEIF(DRAINW(IL).LE.SOILW(IL))THEN
	    IF(SOILW(IL).GT.0)THEN
		  DOCOUT(IL)=DOCSOLV(IL)*DRAINW(IL)/SOILW(IL)
	      DOCSOLV(IL)=DOCSOLV(IL)*(1-(DRAINW(IL)/SOILW(IL)))
		  DONOUT(IL)=DONSOLV(IL)*DRAINW(IL)/SOILW(IL)
	      DONSOLV(IL)=DONSOLV(IL)*(1-(DRAINW(IL)/SOILW(IL)))
          ELSEIF(SOILW(IL).LE.0)THEN
            DOCOUT(IL)=0
            DONOUT(IL)=0
	    ENDIF
        END IF
C
C 2. Add DOC leached from layer above, and pass remainder down to next layer
C
        IF((DRAINW(IL)+REFIL(IL)).GT.0)THEN
          DOCSTOP=DOCIN(IL)*REFIL(IL)/(DRAINW(IL)+REFIL(IL))
          DOCSOLV(IL)=DOCSOLV(IL)+DOCSTOP
	    DOCOUT(IL)=DOCOUT(IL)+DOCIN(IL)-DOCSTOP
          DONSTOP=DONIN(IL)*REFIL(IL)/(DRAINW(IL)+REFIL(IL))
          DONSOLV(IL)=DONSOLV(IL)+DONSTOP
	    DONOUT(IL)=DONOUT(IL)+DONIN(IL)-DONSTOP
	  ELSEIF((DRAINW(IL)+REFIL(IL)).LE.0)THEN
          DOCSTOP=0
	    DOCOUT(IL)=DOCOUT(IL)+DOCIN(IL)
          DONSTOP=0
	    DONOUT(IL)=DONOUT(IL)+DONIN(IL)
        ENDIF
C
C 3. Add up mobile and leached DOC
C
        MOBDOC(IL)=DOCSORB(IL)+DOCSOLV(IL)
	  MOBDON(IL)=DONSORB(IL)+DONSOLV(IL)
	  LEACHDOC(IL)=DOCOUT(IL)
	  LEACHDON(IL)=DONOUT(IL)
100   CONTINUE
      END
C
C------------------------------------------------------------
C
      SUBROUTINE LEACH_NITRATE(SOILN,SOIL15,SLEACH,SLEA15,DRAINW,WMAX,
     &                         SOILW,CRIT,NSOIL,REFIL)
C
C Subroutine to calculate the nitrate leaching
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXSOIL				! Max.no.of soil types
      PARAMETER(MAXSOIL=50)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	INTEGER MAXDEPTH			! Maximum depth of the soil profile
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
      INTEGER IL					! Layer counter
	INTEGER J					! Local counter
	REAL SOL(MAXLAYER)			! Temporary nitrate N in layer kgN/ha
	REAL SOL15(MAXLAYER)		! Temporary nitrate N15 in layer kgN/ha
	REAL T(MAXLAYER)			! N15/N ratio of nitrate pool
	REAL CRITL					! Crit.minimum nitrate in the 
								!  layer (kgN/ha/layer)
	REAL S
	REAL S15
C
C Variables passed from calling subroutine
C
	INTEGER NSOIL				! IN:Soil code number
 	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL CRIT(MAXSOIL,MAXLAYER)	! IN/OUT: Minimum nitrate level (Sozanska)
	REAL DRAINW(MAXLAYER)		! IN:Water drained from this layer (mm/layer)

	REAL SLEACH(MAXLAYER)		! INOUT: Nitrate N leached from this layer 
								!		 to the next (kgN/ha/layer/timestep)
	REAL SLEA15(MAXLAYER)		! INOUT: Nitrate N15 leached from this layer 
								!		 to the next (kgN15/ha/layer/timestep)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL REFIL(MAXLAYER)		! Water deficit in the layer mm/layer
C
C For each layer set NO3 = SOL()
C
      DO 25 IL=1,MAXLAYER1
        SOL(IL)=SOILN(IL)
        SOL15(IL)=SOIL15(IL)
        IF(SOL(IL).GT.0)THEN
          T(IL)=SOL15(IL)/SOL(IL)
        ELSE
          T(IL)=0
        END IF
        SLEACH(IL)=0
        SLEA15(IL)=0
25    CONTINUE
C
C set the minimum amount of nitrate for each layer, CRITL
C
      DO 30 IL=1,MAXLAYER1
        CRITL=CRIT(NSOIL,IL)*MAXDEPTH/(50.*MAXLAYER1)
C
C If not in the top layer pass leached NO3 down from the layer above
C
        IF(IL.GT.1)THEN
          J=IL-1
          S=SLEACH(J)
          S15=SLEACH(J)*T(J)
C
C If in top layer do not pass from the layer above
C
        ELSE
          S=0.0
          S15=0.0
        END IF
C
C If no rainfall do not allow any leaching
C
        IF(DRAINW(IL).LT.0.0001)THEN
          SOL(IL)=SOL(IL)+S
          SOL15(IL)=SOL15(IL)+S15
          S=0
          S15=0
          SLEACH(IL)=0
          SLEA15(IL)=0
          IF(SOL(IL).GT.0)THEN
            T(IL)=SOL15(IL)/SOL(IL)
          ELSE
            T(IL)=0
          END IF
          GOTO 31
        END IF
C
C If Soil NO3 is not zero calculate the proportion N15
C
        IF(SOL(IL)+S.GT.0)THEN
          T(IL)=(SOL15(IL)+S15)/(SOL(IL)+S)
        ELSE
          T(IL)=0
        END IF
C
C Calculate the amount of NO3 leaching from the proportion of the total
C water moving to the next layer
C
        IF(DRAINW(IL)+REFIL(IL).GT.0)THEN
          SOL(IL)=SOL(IL)+S*REFIL(IL)/(DRAINW(IL)+REFIL(IL))
          S=S*(1-REFIL(IL)/(DRAINW(IL)+REFIL(IL)))
	  ENDIF
C
C If the NO3 leached does not exceed the critical minimum
C set the amount left to leach, SLEACH to 0
C
        IF(SOL(IL).LT.CRITL)THEN
          SLEACH(IL)=0
C
C Otherwise, if the amount exceeds the critical level,
C leach nitrate down to the next layer in proportion to the 
C fraction of soil water drained...
C
        ELSE
          IF(DRAINW(IL).GT.SOILW(IL))THEN
            SLEACH(IL)=SOL(IL)
          ELSE
            SLEACH(IL)=SOL(IL)*DRAINW(IL)/SOILW(IL)
          END IF
C
C ...but allow it to drain only to the critical level.
C
          IF(SOL(IL)-(SLEACH(IL)-S).LT.CRITL)THEN
           SLEACH(IL)=(SOL(IL)+S)-CRITL
          END IF
        END IF
C
C Take leached nitrate out of soil layer
C
        SOL(IL)=SOL(IL)-(SLEACH(IL)-S)
        SOL15(IL)=SOL(IL)*T(IL)
        SLEA15(IL)=SLEACH(IL)*T(IL)
   31   CONTINUE
   30 CONTINUE
C
C Reset nitrate in soil layers according to leaching 
C
      DO 400 IL=1,MAXLAYER1
	  SOILN(IL)=SOL(IL)
400   CONTINUE
C
C Leave LEACH_NITRATE
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE LEACH_NH4(AMMN,AMMN15,SLEACH,SLEA15,DRAINW,WMAX,
     &                         SOILW,CRIT,NSOIL,REFIL)
C
C Subroutine to calculate the nitrate leaching
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXSOIL				! Max.no.of soil types
      PARAMETER(MAXSOIL=50)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	INTEGER MAXDEPTH			! Maximum depth of the soil profile
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
      INTEGER IL					! Layer counter
	INTEGER J					! Local counter
	REAL SOL(MAXLAYER)			! Temporary nitrate N in layer kgN/ha
	REAL SOL15(MAXLAYER)		! Temporary nitrate N15 in layer kgN/ha
	REAL T(MAXLAYER)			! N15/N ratio of nitrate pool
	REAL CRITL					! Crit.minimum nitrate in the 
								!  layer (kgN/ha/layer)
	REAL S
	REAL S15
C
C Variables passed from calling subroutine
C
	INTEGER NSOIL				! IN:Soil code number
 	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
	REAL CRIT(MAXSOIL,MAXLAYER)	! IN/OUT: Minimum nitrate level (Sozanska)
	REAL DRAINW(MAXLAYER)		! IN:Water drained from this layer (mm/layer)

	REAL SLEACH(MAXLAYER)		! INOUT: Nitrate N leached from this layer 
								!		 to the next (kgN/ha/layer/timestep)
	REAL SLEA15(MAXLAYER)		! INOUT: Nitrate N15 leached from this layer 
								!		 to the next (kgN15/ha/layer/timestep)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL REFIL(MAXLAYER)		! Water deficit in the layer mm/layer
C
C For each layer set AMMN = SOL()
C
      DO 25 IL=1,MAXLAYER1
        SOL(IL)=AMMN(IL)
        SOL15(IL)=AMMN15(IL)
        IF(SOL(IL).GT.0)THEN
          T(IL)=SOL15(IL)/SOL(IL)
        ELSE
          T(IL)=0
        END IF
        SLEACH(IL)=0
        SLEA15(IL)=0
25    CONTINUE
C
C set the minimum amount of ammonium for each layer, CRITL (temp.set to 10 kg/ha/5cm)
C
      DO 30 IL=1,MAXLAYER1
        CRITL=10 !!temp
C
C If not in the top layer pass leached NH4 down from the layer above
C
        IF(IL.GT.1)THEN
          J=IL-1
          S=SLEACH(J)
          S15=SLEACH(J)*T(J)
C
C If in top layer do not pass from the layer above
C
        ELSE
          S=0.0
          S15=0.0
        END IF
C
C If no rainfall do not allow any leaching
C
        IF(DRAINW(IL).LT.0.0001)THEN
          SOL(IL)=SOL(IL)+S
          SOL15(IL)=SOL15(IL)+S15
          S=0
          S15=0
          SLEACH(IL)=0
          SLEA15(IL)=0
          IF(SOL(IL).GT.0)THEN
            T(IL)=SOL15(IL)/SOL(IL)
          ELSE
            T(IL)=0
          END IF
          GOTO 31
        END IF
C
C If Soil NH4 is not zero calculate the proportion N15
C
        IF(SOL(IL)+S.GT.0)THEN
          T(IL)=(SOL15(IL)+S15)/(SOL(IL)+S)
        ELSE
          T(IL)=0
        END IF
C
C Calculate the amount of NH4 leaching from the proportion of the total
C water moving to the next layer
C
        IF(DRAINW(IL)+REFIL(IL).GT.0)THEN
          SOL(IL)=SOL(IL)+S*REFIL(IL)/(DRAINW(IL)+REFIL(IL))
          S=S*(1-REFIL(IL)/(DRAINW(IL)+REFIL(IL)))
	  ENDIF
C
C If the NH4 leached does not exceed the critical minimum
C set the amount left to leach, SLEACH to 0
C
        IF(SOL(IL).LT.CRITL)THEN
          SLEACH(IL)=0
C
C Otherwise, if the amount exceeds the critical level,
C leach nitrate down to the next layer in proportion to the 
C fraction of soil water drained...
C
        ELSE
          IF(DRAINW(IL).GT.SOILW(IL))THEN
            SLEACH(IL)=SOL(IL)
          ELSE
            SLEACH(IL)=SOL(IL)*DRAINW(IL)/SOILW(IL)
          END IF
C
C ...but allow it to drain only to the critical level.
C
          IF(SOL(IL)-(SLEACH(IL)-S).LT.CRITL)THEN
           SLEACH(IL)=(SOL(IL)+S)-CRITL
          END IF
        END IF
C
C Take leached ammonium out of soil layer
C
        SOL(IL)=SOL(IL)-(SLEACH(IL)-S)
        SOL15(IL)=SOL(IL)*T(IL)
        SLEA15(IL)=SLEACH(IL)*T(IL)
   31   CONTINUE
   30 CONTINUE
C
C Reset ammonium in soil layers according to leaching 
C
      DO 400 IL=1,MAXLAYER1
	  AMMN(IL)=SOL(IL)
400   CONTINUE
C
C Leave LEACH_NH4
C
      RETURN
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE MAKE_DOC(PHP1ARRAY,PHP2ARRAY,PHARRAY,
     &                    NSOIL,ICOVER,ITFUNC,IMFUNC,
     &                    WMAX,WSAT,SOILW,SOILTEMP,
     &                    DPMCARB0,RPMCARB0,BCARB0,HCARB0,
     &                    DPMNIT0,RPMNIT0,BNIT0,HNIT0,
     &                    DPMNLAB0,RPMNLAB0,BNLAB0,HNLAB0,
     &                    AMMN,AMMN15,
     &                    MOBDOC,MOBDON,CO2FROMDOC,HY)
C
C Subroutine to calculate denitrification
C
      IMPLICIT NONE
C
C Variable local to this subroutine
C MAXLAYER = Maximum number of layers in the soil profile
C MAXLAYER1 = Maximum number of layers in the soil profile
C MAXDEPTH = Maximum depth of the soil profile
C IL = layer counter
C MAXSOIL = Maximum number of defined soil types
C DPMRATE = 2e-4 = Rate of change of decomposable plant material into mobile DOC (kg C / ha / day)
C RPMRATE = 4e-5 = Rate of change of resistant plant material into mobile DOC (kg C / ha / day)
C BRATE = 1e-4 = Rate of change of biomass into mobile DOC (kg C / ha / day)
C HRATE = 4e-6 = Rate of change of humus into mobile DOC (kg C / ha / day)
C DOCBIO = 1e-2 = Rate of change of mobile DOC into biomass (kg C / ha / day)
C DOCCO2 = 1e-5 = production of CO2 from mobile DOC (kg C / ha / day)
C RATEMOD = product of each of the rate modifiers previously calculated
C DOCPROD = Change in amount of mobile DOC due to pool decomposition (kg C / ha / layer / day)
C DONPROD = Change in amount of mobile DON due to pool decomposition (kg N / ha / layer / day)
C DOCUP = Change in amount of mobile DOC due to mobile DOC uptake (kg C / ha / layer / day)
C DONUP = Change in amount of mobile DON due to mobile DON uptake (kg N / ha / layer / day)
C TRATE = Rate modifier according to soil temperature
C WRATE = Rate modifier according to soil water
C PHRATE = Rate modifier according to soil pH
C CRRATE = Rate modifier according to crop cover
C PH = pH of this soil
C PHP1 = pH below which decomposition stops
C PHP2 = pH above which decomposition rate is not modified by pH
C T_GERM = timesteps needed for germination
C
	INTEGER MAXLAYER,MAXDEPTH,MAXLAYER1
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
	INTEGER IL
      INTEGER MAXSOIL
      PARAMETER(MAXSOIL=50)
      REAL*8 DPMRATE,RPMRATE,BRATE,HRATE,DOCBIO,DOCCO2
	DATA DPMRATE, RPMRATE, BRATE,  HRATE,   DOCBIO,DOCCO2 
     &    /0.0002,0.00004,0.0001,0.000004,0.005, 0.00001/
	REAL RATEMOD,MOBPART1,MOBPART2,MOBPART3
      REAL DOCPROD,DOCUP,DONUP,DONPROD
      REAL CRRATE,TRATE,WRATED,WRATER,PHRATE,PH,PHP1,PHP2
C
C Variables passed to/from calling subroutine
C WMAX() = Maximum water content of this layer (mm / layer)
C SOILW() = Available water content of this layer (mm / layer)
C AIRTEMP = Air temperature (deg.C)
C PHARRAY() = Array holding pH of defined soils
C PHP1ARRAY() = Array holding pH below which decomposition stops for defined soils
C PHP2ARRAY() = Array holding pH above which decomposition rate is not modified by pH for defined soils
C NSOIL = Soil type number code
C HY = inverse of C:N ratio
C DPMCARB0() =  Amount of decomposable plant material C in layer (kg C / ha / layer)
C RPMCARB0() =  Amount of resistant plant material C in layer (kg C / ha / layer)
C BCARB0() = Amount of biomass C in layer (kg C / ha / layer)
C HCARB0() = Amount of humus C in layer (kg C / ha / layer)
C MOBDOC() = Total amount of mobile DOC pool (kg C / ha / layer)
C MOBDON() = Total amount of mobile DON pool (kg N / ha / layer)
C CO2FROMDOC() = CO2 production rate of layer (kg C / ha / layer / day)
C
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER),PHP2ARRAY(MAXSOIL,MAXLAYER)
	REAL PHARRAY(MAXSOIL,MAXLAYER)
	INTEGER NSOIL,ICOVER
      REAL WMAX(MAXLAYER),SOILW(MAXLAYER),SOILTEMP(MAXLAYER)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
	REAL HY(MAXSOIL,MAXLAYER)
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
	REAL DPMCARB0(MAXLAYER),RPMCARB0(MAXLAYER)
	REAL BCARB0(MAXLAYER),HCARB0(MAXLAYER)
 	REAL DPMNIT0(MAXLAYER),RPMNIT0(MAXLAYER)
	REAL BNIT0(MAXLAYER),HNIT0(MAXLAYER)
 	REAL DPMNLAB0(MAXLAYER),RPMNLAB0(MAXLAYER)
	REAL BNLAB0(MAXLAYER),HNLAB0(MAXLAYER)
      REAL MOBDOC(MAXLAYER),CO2FROMDOC(MAXLAYER),MOBDON(MAXLAYER)
C
C ...Model descriptors
C
	INTEGER ITFUNC				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IMFUNC				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
C
C Calculation for every layer down the profile...
C
      DO 100 IL=1,MAXLAYER1
C
C Calculate rate modifiers (Source Bradbury et al, 1993)
C
        CALL MODFACTS_MINER(SOILW(IL),WMAX(IL),WSAT(IL),
     &                WRATED,WRATER,
     &                SOILTEMP(IL),TRATE,								
     &			    PHARRAY(NSOIL,IL),PHP1ARRAY(NSOIL,IL),
     &                PHP2ARRAY(NSOIL,IL),PHRATE,
     &                ICOVER,CRRATE,ITFUNC,IMFUNC)									
C
C 5. Overall rate modifier
C      
        RATEMOD=((WRATED+WRATER)/2)*TRATE*PHRATE*CRRATE
C
C Breakdown of RPM into mobile DOC 
C (after Michalzic et al, 2003. Biogeochem. 66,241-264)
C Assumptions:	1. RPM breakdown produces DOC only
C                 2. A constant fraction breaks down each timestep, modified 
C                     by environmental factors 
C				3. The rate modifiers used for the first order process of microbial 
C                     decomposition can be used unchanged to describe rate modification.
C
        DOCPROD=(RPMCARB0(IL)*RPMRATE*RATEMOD)
        DONPROD=(RPMNIT0(IL)*RPMRATE*RATEMOD)
        RPMCARB0(IL)=RPMCARB0(IL)-DOCPROD
        RPMNLAB0(IL)=RPMNLAB0(IL)-(RPMNLAB0(IL)/RPMNIT0(IL))*DONPROD
        RPMNIT0(IL)=RPMNIT0(IL)-DONPROD
        MOBDOC(IL)=MOBDOC(IL)+DOCPROD
        MOBDON(IL)=MOBDON(IL)+DONPROD
C
C Breakdown of DPM into mobile DOC 
C (after Michalzic et al, 2003. Biogeochem. 66,241-264)
C Assumptions:	1. DPM breakdown produces DOC only
C                 2. A constant fraction breaks down each timestep, modified 
C                     by environmental factors
C				3. The rate modifiers used for the first order process of microbial 
C                     decomposition can be used unchanged to describe rate modification.
C
        DOCPROD=(DPMCARB0(IL)*DPMRATE*RATEMOD)
        DONPROD=(DPMNIT0(IL)*DPMRATE*RATEMOD)
        DPMCARB0(IL)=DPMCARB0(IL)-DOCPROD
        DPMNLAB0(IL)=DPMNLAB0(IL)-(DPMNLAB0(IL)/DPMNIT0(IL))*DONPROD
        DPMNIT0(IL)=DPMNIT0(IL)-DONPROD
        MOBDOC(IL)=MOBDOC(IL)+DOCPROD
        MOBDON(IL)=MOBDON(IL)+DONPROD
C
C Breakdown of biomass into mobile DOC 
C (after Michalzic et al, 2003. Biogeochem. 66,241-264)
C Assumptions:	1. Biomass breakdown produces DOC only
C                 2. A constant fraction breaks down each timestep, modified 
C                     by environmental factors
C				3. The rate modifiers used for the first order process of microbial 
C                     decomposition can be used unchanged to describe rate modification.
C
        DOCPROD=(BCARB0(IL)*BRATE*RATEMOD)
        DONPROD=(BNIT0(IL)*BRATE*RATEMOD)
        BCARB0(IL)=BCARB0(IL)-DOCPROD
        BNLAB0(IL)=BNLAB0(IL)-(BNLAB0(IL)/BNIT0(IL))*DONPROD
        BNIT0(IL)=BNIT0(IL)-DONPROD
        MOBDOC(IL)=MOBDOC(IL)+DOCPROD
        MOBDON(IL)=MOBDON(IL)+DONPROD
C
C Breakdown of humus into mobile DOC 
C (after Michalzic et al, 2003. Biogeochem. 66,241-264)
C Assumptions:	1. Humus breakdown produces DOC only
C                 2. A constant fraction breaks down each timestep, modified 
C                     by environmental factors
C				3. The rate modifiers used for the first order process of microbial 
C                     decomposition can be used unchanged to describe rate modification.
C
        DOCPROD=(HCARB0(IL)*HRATE*RATEMOD)
        DONPROD=(HNIT0(IL)*HRATE*RATEMOD)
        HCARB0(IL)=HCARB0(IL)-DOCPROD
        HNLAB0(IL)=HNLAB0(IL)-(HNLAB0(IL)/HNIT0(IL))*DONPROD
        HNIT0(IL)=HNIT0(IL)-DONPROD
        MOBDOC(IL)=MOBDOC(IL)+DOCPROD
        MOBDON(IL)=MOBDON(IL)+DONPROD
C
C Uptake by biomass of mobile DOC
C Assumptions 1. A constant fraction breaks down each timestep, modified 
C                     by environmental factors
C             2. Only biomass takes up DOC (note this is a route from humus to biomass)
C
        DOCUP=(MOBDOC(IL)*DOCBIO*RATEMOD)
        DONUP=(MOBDON(IL)*DOCBIO*RATEMOD)
        MOBDOC(IL)=MOBDOC(IL)-DOCUP
        MOBDON(IL)=MOBDON(IL)-DONUP
        BCARB0(IL)=BCARB0(IL)+DOCUP
        BNLAB0(IL)=BNLAB0(IL)+(BNLAB0(IL)/BNIT0(IL))*DONUP
        BNIT0(IL)=BNIT0(IL)+DONUP
C
C Production of CO2 from mobile DOC
C Assumptions: 1. A constant fraction is taken up by biomass at each timestep 
C
        CO2FROMDOC(IL)=MOBDOC(IL)*DOCCO2
        MOBDOC(IL)=MOBDOC(IL)*(1-DOCCO2)
	  MOBDON(IL)=MOBDON(IL)*(1-DOCCO2)
	  AMMN15(IL)=AMMN15(IL)+(AMMN15(IL)/AMMN(IL))*(MOBDON(IL)*DOCCO2)
	  AMMN(IL)=AMMN(IL)+(MOBDON(IL)*DOCCO2)
100   CONTINUE
      END
C
C-------------------------------------------------------------
C
	SUBROUTINE MAKE_PAR_CN(NSOIL,TOC,IOM,ILU,CLAY,BULKDENS,
     &                   HY,BPART,BPROP,HPART,HPROP,
     &				   BIOP,CRIT,BIORATE,HUMRATE,TOCARRAY,IOMARRAY,
     &				   SNAME,CLARRAY,LUARRAY,
     &				   PHARRAY,PHP1ARRAY,PHP2ARRAY,DFACT,
     &				   FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,FYMLOSS,FYMPOS,
     &				   FYMSTART)
C
C Subroutine to reset parameters from files PARAM.OUT and PARAM.DAT
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXSOIL				! Max.no.of soil types
	INTEGER MAXORGM				! Max.no.of organic manure applications allowed
	PARAMETER(MAXSOIL=50,MAXORGM=52)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	INTEGER MAXDEPTH			! Maximum depth of the soil profile
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
	INTEGER IL					! Layer counter
C
C Variables passed to/from calling subroutine
C ...Soil factors
C
	INTEGER NSOIL				! IN:Soil code number
	REAL TOC(MAXLAYER)			! IN:Total organic C (kgC/ha/layer)
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
      REAL CLAY(MAXLAYER)			! IN:Clay content of the layer (%)
	REAL BULKDENS(MAXLAYER)		! IN:Bulk density of the layer (g/cm3)
	INTEGER ILU					! IN:Counter for land use 
								!    (1=arable;2=grass;3=forestry;4=semi-nat)
	REAL HY(MAXSOIL,MAXLAYER)	! OUT:Stable N:C ratio of BIO & HUM pools
	REAL BPART(MAXSOIL,MAXLAYER)! OUT:Decomposition efficiency 
								! = (BIO + HUM) / (Total decomposition)
	REAL BPROP(MAXSOIL,MAXLAYER)! OUT: BIO/HUM from biomass decompositn
	REAL HPART(MAXSOIL,MAXLAYER)! OUT:Decomposition efficiency 
								! = (BIO + HUM) / (Total decomposition)
	REAL HPROP(MAXSOIL,MAXLAYER)! OUT: BIO/HUM from humus decompositn
	REAL BIOP(MAXSOIL,MAXLAYER)	! OUT: BIO/TOC 
	REAL CRIT(MAXSOIL,MAXLAYER)	! OUT: Minimum nitrate level / 50cm 
								! (Sozanska eqn)
	REAL BIORATE(MAXSOIL,MAXLAYER)	! OUT: Rate constant for BIO decompn/yr
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! OUT: Rate constant for HUM decompn/yr
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! OUT: Equilibrium TOC (kgC/ha/layer)
	REAL IOMARRAY(MAXSOIL,MAXLAYER)	! OUT: Equilibrium IOM (kgC/ha/layer)
	CHARACTER*40 SNAME(MAXSOIL)	! OUT:Soil name
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! OUT: % clay content in this layer
      INTEGER LUARRAY(MAXSOIL)	! OUT:Land use before equilibrium 
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)! OUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)! OUT: pH above which decomp.max.
	REAL DFACT(MAXLAYER)		! OUT: Denitrification factor
C
C ...Manure factors
C
      REAL FPROPX(MAXORGM)		! OUT:Amount of N in fresh manure that is decomposed organic matter (kg N / t manure)
	REAL FYMCX(MAXORGM)			! OUT:Amount of C in fresh manure (kg C / t manure)
	REAL FYMNX(MAXORGM)			! OUT:Amount of N in fresh manure (kg N / t manure)
	REAL FYMAX(MAXORGM)			! OUT:Amount of NH4+ in fresh manure (kg N / t manure)
	REAL FYMWX(MAXORGM)			! OUT:Amount of water in fresh manure (mm water / t manure)
	REAL FYMLOSS(MAXORGM)		! OUT:Proportion of organic waste lost by volatn/timestep
	REAL FYMPOS(MAXORGM)		! OUT:Proportion of FYM added to top 25cm
	REAL FYMSTART(MAXORGM)		! OUT:Amount of rainfall in 1 week, 
								! 	 below which volatilisation will occur
C
C Soil Parameters
C
      SNAME(NSOIL)='Calculated soil parameters'
C
C Layered parameters
C
      DO 100 IL=1,MAXLAYER1
	  TOCARRAY(NSOIL,IL)=TOC(IL)
C IOM from Falloon equation
        IF(IOM(IL).LT.0)THEN
          CALL GET_IOM_FROM_FALLOON_EQN(TOCARRAY(NSOIL,IL),
     &                                  IOMARRAY(NSOIL,IL))
	  ELSE
	    IOMARRAY(NSOIL,IL)=IOM(IL)
	  ENDIF
C Clay content
        CLARRAY(NSOIL,IL)=CLAY(IL)
C HY(MAXSOIL,MAXLAYER) = Stable N:C ratio of biomass & humus
        HY(NSOIL,IL)=0.118
C BPART(MAXSOIL,MAXLAYER) = decomposition efficiency = (BIO + HUM) / (Total decomposition)
C HPART(MAXSOIL,MAXLAYER) = decomposition efficiency = (BIO + HUM) / (Total decomposition)
c If BPART or HPART are changed also change in BPCALC
        CALL GET_DECOMP_EFF(CLAY(IL),BPART(NSOIL,IL))
c       BPART(NSOIL,IL)=1/(1+1.67*(1.85+1.6*EXP(-0.0786*CLAY(IL))))      
        HPART(NSOIL,IL)=BPART(NSOIL,IL)
C BPROP(MAXSOIL,MAXLAYER) = BIO/HUM from biomass decomposition
c        BPROP(NSOIL,IL)=1.1
        BPROP(NSOIL,IL)=0.85
C HPROP(MAXSOIL,MAXLAYER) = BIO/HUM from humus decomposition
c        HPROP(NSOIL,IL)=1.1
        HPROP(NSOIL,IL)=0.85
C BIOP(MAXSOIL,MAXLAYER) = BIO/TOC 
        BIOP(NSOIL,IL)=0.028
C CRIT(MAXSOIL,MAXLAYER) = Minimum nitrate level / 50cm (Sozanska eqn)
c        CRIT(NSOIL,IL)=5+(CLAY(IL)*10./65.)
        CRIT(NSOIL,IL)=0
C BIORATE(MAXSOIL,MAXLAYER) = Rate constant for biomass decomposition (/year)
        BIORATE(NSOIL,IL)=0.6604/52.
C HUMRATE(MAXSOIL,MAXLAYER) = Rate constant for humus decomposition (/year)
        HUMRATE(NSOIL,IL)=0.0208/52.
C pH = currently set to neutral - not used
        PHARRAY(NSOIL,IL)=7
	  PHP1ARRAY(NSOIL,IL)=1.0
	  PHP2ARRAY(NSOIL,IL)=4.5
C Denitrification factor
        DFACT(IL)=0.005
100   CONTINUE
C
C Non-layered parameters
C Land use change arrays        
	LUARRAY(NSOIL)=ILU
C
C Organic Manure Parameters
C Manure types
C (1) Cattle FYM
C (2) Pig FYM
C (3) Layer Manure
C (4) Broiler/Turkey Manure
C (5) Sewage sludge cake (undigested)
C (6) Sewage sludge cake (digested)
C (7) Dairy slurry
C (8) Beef slurry
C (9) Pig Slurry
C (10) Strainer box seperated slurry
C (11) Weeping wall seperated slurry
C (12) Mechanically seperated slurry
C (13) Sewage sludge liquid (undigested)
C (14) Sewage sludge liquid (digested)
C
C FPROPX - Amount of N in fresh manure that is decomposed organic matter (kg N / t manure)
C
      FPROPX(1)=100
      FPROPX(2)=100
      FPROPX(3)=250
      FPROPX(4)=350
      FPROPX(5)=350
      FPROPX(6)=350
      FPROPX(7)=200
      FPROPX(8)=200
      FPROPX(9)=320
      FPROPX(10)=200
      FPROPX(11)=250
      FPROPX(12)=200
      FPROPX(13)=150
      FPROPX(14)=150

C     FPROPX(1)=0.2
C     FPROPX(2)=0.2
C     FPROPX(3)=0.2
C     FPROPX(4)=0.2	
C     FPROPX(5)=0.2
C     FPROPX(6)=0.2
C     FPROPX(7)=0.2

C
C FYMCX - C in fresh manure (kg C / t manure)
C
	FYMCX(1)=52.5
	FYMCX(2)=62.5
	FYMCX(3)=96.0
	FYMCX(4)=180.0
	FYMCX(5)=62.5
	FYMCX(6)=62.5
	FYMCX(7)=12.0
	FYMCX(8)=12.0
	FYMCX(9)=15.0
	FYMCX(10)=4.5
	FYMCX(11)=9.0
	FYMCX(12)=12.0
	FYMCX(13)=12.5
	FYMCX(14)=12.0

C	FYMCX(1)=80.0
C	FYMCX(2)=80.0
C	FYMCX(3)=80.0
C	FYMCX(4)=12.0
C	FYMCX(5)=80.0
C	FYMCX(6)=12.0
C	FYMCX(7)=24.0

C
C FYMNX - N in fresh manure (kg N / t manure)
C
      FYMNX(1)=6.0
      FYMNX(2)=7.0
      FYMNX(3)=15.0
      FYMNX(4)=28.8
      FYMNX(5)=7.5
      FYMNX(6)=7.5
      FYMNX(7)=3.0
      FYMNX(8)=2.28
      FYMNX(9)=4.98
      FYMNX(10)=1.5
      FYMNX(11)=3.0
      FYMNX(12)=4.0
      FYMNX(13)=1.8
      FYMNX(14)=2.0

C     FYMNX(1)=5.5
C     FYMNX(2)=5.5
C     FYMNX(3)=5.5
C     FYMNX(4)=5.5
C     FYMNX(5)=5.5
C     FYMNX(6)=5.5
C     FYMNX(7)=8.2

C
C FYMAX - NH4+ in fresh manure (kg N / t manure)
C
	FYMAX(1)=1.2
	FYMAX(2)=1.4
	FYMAX(3)=7.5
	FYMAX(4)=14.4
	FYMAX(5)=1.5
	FYMAX(6)=1.125
	FYMAX(7)=1.71
	FYMAX(8)=1.14
	FYMAX(9)=3.486
	FYMAX(10)=0.75
	FYMAX(11)=1.5
	FYMAX(12)=2
	FYMAX(13)=0.81
	FYMAX(14)=0.9

C     FYMAX(1)=1.5
C     FYMAX(2)=1.5
C     FYMAX(3)=1.5
C     FYMAX(4)=1.5
C     FYMAX(5)=1.5
C     FYMAX(6)=1.5
C     FYMAX(7)=4.4

C
C FYMWX - Water in fresh manure (mm water / t manure)
C
	FYMWX(1)=0.075
	FYMWX(2)=0.075
	FYMWX(3)=0.07
	FYMWX(4)=0.04
	FYMWX(5)=0.075
	FYMWX(6)=0.075
	FYMWX(7)=0.094
	FYMWX(8)=0.094
	FYMWX(9)=0.094
	FYMWX(10)=0.0985
	FYMWX(11)=0.097
	FYMWX(12)=0.096
	FYMWX(13)=0.095
	FYMWX(14)=0.096

C     FYMWX(1)=0.077
C     FYMWX(2)=0.077
C     FYMWX(3)=0.077
C     FYMWX(4)=0.097
C     FYMWX(5)=0.077
C     FYMWX(6)=0.097
C     FYMWX(7)=0.094

C
C FYMLOSS Proportion of organic waste lost by volatn/timestep
C

      FYMLOSS(1)=0.05
      FYMLOSS(2)=0.05
      FYMLOSS(3)=0.15
      FYMLOSS(4)=0.15
      FYMLOSS(5)=0.05
      FYMLOSS(6)=0.0375
      FYMLOSS(7)=0.228
      FYMLOSS(8)=0.2
      FYMLOSS(9)=0.294
      FYMLOSS(10)=0.1
      FYMLOSS(11)=0.1
      FYMLOSS(12)=0.1
      FYMLOSS(13)=0.18
      FYMLOSS(14)=0.18

C     FYMLOSS(1)=0.15
C     FYMLOSS(2)=0.15
C     FYMLOSS(3)=0.15
C     FYMLOSS(4)=0.15
C     FYMLOSS(5)=0.15
C     FYMLOSS(6)=0.15
C     FYMLOSS(7)=0.15

C
C FYMPOS Proportion of FYM added to top 25cm 
C 
	FYMPOS(1)=1.0
	FYMPOS(2)=1.0
	FYMPOS(3)=1.0
	FYMPOS(4)=1.0
	FYMPOS(5)=1.0
	FYMPOS(6)=1.0
	FYMPOS(7)=0.8
	FYMPOS(8)=0.8
	FYMPOS(9)=0.8
	FYMPOS(10)=0.8
	FYMPOS(11)=0.8
	FYMPOS(12)=0.8
	FYMPOS(13)=0.8
	FYMPOS(14)=0.8

C     FYMPOS(1)=1.0
C     FYMPOS(2)=1.0
C     FYMPOS(3)=1.0
C     FYMPOS(4)=0.8
C     FYMPOS(5)=1.0
C     FYMPOS(6)=0.8
C     FYMPOS(7)=0.8

C
C FYMSTART Amount of rainfall in 1 week, below which volatilisation will occur
C
      FYMSTART(1)=5
      FYMSTART(2)=5
      FYMSTART(3)=5
      FYMSTART(4)=5
      FYMSTART(5)=5
      FYMSTART(6)=5
      FYMSTART(7)=5
      FYMSTART(8)=5
      FYMSTART(9)=5
      FYMSTART(10)=5
      FYMSTART(11)=5
      FYMSTART(12)=5
      FYMSTART(13)=5
      FYMSTART(14)=5
C
C Leave MAKEPAR_CN
C
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE MINER1_FOEREID(NSOIL,ICOVER,ITFUNC,IMFUNC,
     &                    PI_C,PI_N,PI_N15,DRRAT,								
     &                    SOILW,WMAX,WSAT,SOILTEMP,
     &                    ALPHA,BETA,GAMMA,DELTA,HY,
     &                    BRATE,HRATE,DPMRATE,RPMRATE,
     &                    AMMN,AMMN15,SOILN,SOIL15,CRIT,
     &                    CO2,
     &                    DNIT,DNIT15,TDNIT,T15,
     &                    BCARB0,BNIT0,BNLAB0,
     &                    HCARB0,HNIT0,HNLAB0,
     &                    DPMCARB0,DPMNIT0,DPMNLAB0,
     &                    RPMCARB0,RPMNIT0,RPMNLAB0,
     &                    PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &                    DPMCTON,RPMCTON,
C Required for Monte Carlo - 
     &					DPMCARB,RPMCARB,BCARB,HCARB,
     &					WRATEDM,TRATEM)
C
C Subroutine to calculate mineralisation and immobilisation
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXSOIL				! Max.no.of soil types
      PARAMETER(MAXSOIL=50)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	INTEGER MAXDEPTH			! Maximum depth of the soil profile
	INTEGER IDEPTH				! Depth of this layer
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/

	REAL ALOW					! Limited parameter for biomass decomposition
	REAL AN						! Amount ammonium avail.for immob.(kgN/ha)
 	REAL AN15					! Amount ammonium-15 avail.for immob.(kgN/ha)
	REAL ANIT					! Ammonium immobilised (kgN/ha)
	REAL BCARB(MAXLAYER)		! C in soil biomass at end (kgC/ha/layer)
      REAL BNIT(MAXLAYER)			! N in soil biomass at end (kgN/ha/layer)
	REAL BNLAB(MAXLAYER)		! N15 in soil humus at end (kgN15/ha/layer)
	REAL CARB(MAXLAYER)			! Carbon input to layer (kgC/ha)
	REAL CLOW					! Limited parameter for humus decomposition
	REAL CRITL					! Critical minimum N in the layer (kgN/ha)
	REAL DPMCARB(MAXLAYER)		! C in decomposable PM at end (kgC/ha/layer)
	REAL DPMNIT(MAXLAYER)		! N in decomposable PM at end (kgN/ha/layer)
	REAL DPMNLAB(MAXLAYER)		! N15 in decomposable PM at end (kgN15/ha/layer)
      REAL DPMNRAT				! Ratio of plant matter N in DPM
	REAL*8 EB					! Proportion of biomass decomposed 
	REAL*8 EDPM					! Proportion of DPM decomposed 
	REAL*8 EH					! Proportion of humus decomposed 
	REAL*8 ERPM					! Proportion of RPM decomposed 
	REAL HCARB(MAXLAYER)		! C in soil humus at end (kgC/ha/layer)
	REAL HNIT(MAXLAYER)			! N in soil humus at end (kgN/ha/layer)
      REAL HNLAB(MAXLAYER)		! N15 in soil biomass at end (kgN15/ha/layer)
	INTEGER LIMITED				! Marker to indicate immobilisation is limited
	INTEGER M					! Layer counter
	REAL P(MAXLAYER)			! N15/N of humus
	REAL PH						! Current pH of the soil
	REAL PHP1					! pH below which rate is zero
	REAL PHP2					! pH above which rate is optimum
	REAL PROPIMM				! Proportion of mineral N available for immob.												
	REAL Q(MAXLAYER)			! N15/N of PM
	REAL RADD(MAXLAYER)			! Nitrogen input to layer (kgN/ha)
	REAL RADD15(MAXLAYER)		! N15 input to layer (kgN15/ha)
	REAL RPMCARB(MAXLAYER)		! C in resistant PM at end (kgC/ha/layer)
	REAL RPMNIT(MAXLAYER)		! N in resistant PM at end (kgN/ha/layer)
	REAL RPMNLAB(MAXLAYER)		! N15 in resistant PM at end (kgN15/ha/layer)
      REAL RPMNRAT				! Ratio of plant matter N in RPM
	REAL S(MAXLAYER)			! N15/N of biomass
	REAL SN						! Amount nitrate avail.for immob.(kgN/ha)
	REAL SPL(MAXLAYER)			! Split of input into the layer
	REAL SN15					! Amount nitrate-15 avail.for immob.(kgN/ha)
	REAL SNIT					! Nitrate immobilised (kgN/ha)
	REAL T(MAXLAYER)			! N15/N of mineral N
	REAL TA(MAXLAYER)			! N15/N of ammonium
	REAL TN(MAXLAYER)			! N15/N of nitrate
	REAL U(MAXLAYER)			! Variable used in calculation of 15N loss
	REAL XN						! Variable used in calc.of changes in pools
	REAL Z(MAXLAYER)			! N:C of PM
C
C Variables passed to/from calling subroutine
C ...Weather factors
C
	REAL SOILTEMP(MAXLAYER)		! IN: Soil temperature (deg.C/timestep)
C
C ...Model descriptors
C
	INTEGER ITFUNC				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IMFUNC				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
C
C ... Soil factors
C
	REAL ALPHA(MAXLAYER)			! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
  	REAL BCARB0(MAXLAYER)		! IN:C in soil biomass (kgC/ha/layer)
	REAL BETA(MAXLAYER)				! IN/OUT: Prop.HUM produced on BIO decompn
      REAL BNIT0(MAXLAYER)		! IN:N in soil biomass (kgN/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN:N15 in soil humus (kgN15/ha/layer)
	REAL BRATE(MAXLAYER)			! IN/OUT: Rate constant for HUM decompn
	REAL CO2(MAXLAYER)			! IN:CO2 emitted (kgC/ha/layer)
	REAL CRIT(MAXSOIL,MAXLAYER)		! IN/OUT: Minimum nitrate level (Sozanska)
      REAL CRRATE					! IN:crop modifying factor (proportion)
	REAL DELTA(MAXLAYER)			! IN/OUT: Prop.HUM produced on HUM decompn 
	REAL DNIT(MAXLAYER)			! IN:Net Mineralised N (kgN/ha/lay/timestep)
	REAL DNIT15(MAXLAYER)		! IN:Net Mineralised N15 (kgN15/ha/lay/time)
	REAL DRRAT					! IN:DPM:RPM ratio
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
      REAL DPMCTON(MAXLAYER)		! IN:DPM C:N ratio
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL DPMRATE(MAXLAYER)			! IN/OUT: Rate constant for DPM decompn
	REAL GAMMA(MAXLAYER)			! IN/OUT: Prop.BIO produced on HUM decompn 
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN:N in soil humus (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN:N15 in soil biomass (kgN15/ha/layer)
	REAL HRATE(MAXLAYER)			! IN/OUT: Rate constant for BIO decompn
	REAL HY(MAXSOIL,MAXLAYER)	! IN/OUT:Stable N:C ratio of BIO & HUM pools
  	REAL ICFACTOR(MAXLAYER)		! IN:Adjustment in rate needed to	
								!    achieve measured NPP and TOC	
      INTEGER ICOVER				! IN:Code for crop cover: 1=covered; 2=bare
	INTEGER NSOIL				! IN:Soil code number
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)	! IN/PUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)	! IN/PUT: pH above which decomp.max.
      REAL PHRATE					! IN:pH rate modifier (proportion)
      REAL PI_C(MAXLAYER)			! IN/OUT: Plant input C to soil (kgC/ha/step)
      REAL PI_N(MAXLAYER)			! IN/OUT: Plant input N to soil (kgC/ha/step)
      REAL PI_N15(MAXLAYER)		! IN/OUT: Plant input N15 to soil (kgC/ha/step)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
      REAL RPMCTON(MAXLAYER)		! IN:RPM C:N ratio
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL RPMRATE(MAXLAYER)			! IN/OUT: Rate constant for RPM decompn
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL T15					! IN:Net Mineralised N (kgN/ha/layer)
	REAL TDNIT					! IN:Net Mineralised N (kgN/ha/layer)
	REAL TRATE					! IN:temperature rate modifying factor (prop.)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WRATED					! IN:Moisture rate modifyer for DPM & BIO(prop)
	REAL WRATER					! IN:Moisture rate modifyer for RPM & HUM(prop)
	REAL TRATEM					! Temperature rate modifier summed for layers
	REAL WRATEDM				! Moisture rate modifier summed for layers
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
C
C Save Variable Values
C
      SAVE
C
C  Calculate the amount added from plant debris this week of
C  C(=CARB(M)), N(=RADD(M)) and N15(=RADD15(M)) to layer M
C
      DO 24 M=1,MAXLAYER1
	  IF(PI_C(M).LT.0)PI_C(M)=0
	  IF(PI_N(M).LT.0)PI_N(M)=0
	  IF(PI_N15(M).LT.0)PI_N15(M)=0
        CARB(M)=PI_C(M)
       IF(PI_C(M).GT.0)THEN
	 ICOVER=1
	 ELSE
	 ICOVER=0
	 ENDIF
        RADD(M)=PI_N(M)
        RADD15(M)=PI_N15(M)
   24 CONTINUE
C
C For each layer...
C
C
	WRATEDM=0
	TRATEM=0
C
      DO 25 M=1,MAXLAYER1
        CALL MODFACTS_MINER(SOILW(M),WMAX(M),WSAT(M),
     &                WRATED,WRATER,
     &                SOILTEMP(M),TRATE,
     &			    PHARRAY(NSOIL,M),PHP1ARRAY(NSOIL,M),
     &                PHP2ARRAY(NSOIL,M),PHRATE,
     &                ICOVER,CRRATE,ITFUNC,IMFUNC)
C
C
C
  		IF (M.LT.11) THEN
		WRATEDM=WRATEDM+WRATED
		TRATEM=TRATEM+TRATE
		ENDIF
C
C
C Calculate the exponential factors for decomposition of
C the RO(=ER), BIO(=EB) and HUM(=EH) pools
C
        EDPM=EXP(-DPMRATE(M)*WRATED*TRATE*CRRATE*PHRATE*ICFACTOR(M))
	  ERPM=EXP(-RPMRATE(M)*WRATER*TRATE*CRRATE*PHRATE*ICFACTOR(M))
        EB=EXP(-BRATE(M)*WRATED*TRATE*CRRATE*PHRATE*ICFACTOR(M))
        EH=EXP(-HRATE(M)*WRATER*TRATE*CRRATE*PHRATE*ICFACTOR(M))
C
C Add to the RO pool in the layer
C
c------------------------------------------------------------
c Previous approach N -> DPM:RPM in 1:1 ratio
c        DPMCARB0(M)=DPMCARB0(M)+CARB(M)*(DRRAT/(1+DRRAT))
c        RPMCARB0(M)=RPMCARB0(M)+CARB(M)*(1/(1+DRRAT))
c        DPMNIT0(M)=DPMNIT0(M)+RADD(M)*(DRRAT/(1+DRRAT))
c        RPMNIT0(M)=RPMNIT0(M)+RADD(M)*(1/(1+DRRAT))
c------------------------------------------------------------
        DPMCARB0(M)=DPMCARB0(M)+CARB(M)*(DRRAT/(1+DRRAT))
        RPMCARB0(M)=RPMCARB0(M)+CARB(M)*(1/(1+DRRAT))
	  IF(CARB(M).GT.0)THEN
	   DPMNRAT=CARB(M)*((DRRAT/(1+DRRAT))/DPMCTON(M))
	   RPMNRAT=CARB(M)*((1/(1+DRRAT))/RPMCTON(M))
         DPMNIT0(M)=DPMNIT0(M)+(RADD(M)*(DPMNRAT/(DPMNRAT+RPMNRAT)))
         RPMNIT0(M)=RPMNIT0(M)+(RADD(M)*(RPMNRAT/(DPMNRAT+RPMNRAT)))
         DPMNLAB0(M)=DPMNLAB0(M)+(RADD15(M)*(DPMNRAT/(DPMNRAT+RPMNRAT)))
         RPMNLAB0(M)=RPMNLAB0(M)+(RADD15(M)*(RPMNRAT/(DPMNRAT+RPMNRAT)))
	  ENDIF
	  IF((DPMCARB0(M)+RPMCARB0(M)).GT.0)THEN
          Z(M)=(DPMNIT0(M)+RPMNIT0(M))/(DPMCARB0(M)+RPMCARB0(M))
	  ELSE
	    Z(M)=0
	  ENDIF
C
C Calc N15/N for NH4(=TA()), NO3(=TN()), MinN(=T()), RO(=Q()), BIO(=S()) and HUM(=P())
C
        TA(M)=AMMN15(M)/AMMN(M)
        TN(M)=SOIL15(M)/SOILN(M)
        T(M)=(AMMN15(M)+SOIL15(M))/(AMMN(M)+SOILN(M))
	  Q(M)=0
	  S(M)=0
        P(M)=0
	  IF((DPMNIT0(M)+RPMNIT0(M)).GT.0)
     &    Q(M)=(DPMNLAB0(M)+RPMNLAB0(M))/(DPMNIT0(M)+RPMNIT0(M))
	  IF(BNIT0(M).GT.0)S(M)=BNLAB0(M)/BNIT0(M)
	  IF(HNIT0(M).GT.0)P(M)=HNLAB0(M)/HNIT0(M)
C
C Calculation of Pool Changes...
C
        LIMITED=0
  100   CONTINUE
C
C Adjust C(=RCARB()), N(=RNIT()) and N15(=RNLAB()) in RO pool
C due to biological activity
C
        DPMCARB(M)=DPMCARB0(M)*EDPM
        RPMCARB(M)=RPMCARB0(M)*ERPM
        DPMNIT(M)=DPMNIT0(M)*EDPM
        RPMNIT(M)=RPMNIT0(M)*ERPM
        DPMNLAB(M)=DPMNLAB0(M)*EDPM
        RPMNLAB(M)=RPMNLAB0(M)*ERPM
C
C Adjust C(=BCARB()) in BIO pool due to biological activity
C
        BCARB(M)=(BCARB0(M)*EB)+
     &         (ALPHA(M)*BCARB0(M)*(1-EB))+
     &         (GAMMA(M)*HCARB0(M)*(1-EH))+
     &         (ALPHA(M)*DPMCARB0(M)*(1-EDPM))+
     &         (ALPHA(M)*RPMCARB0(M)*(1-ERPM))
C
C Adjust C(=HCARB()) in HUM pool due to biological activity
C
        HCARB(M)=(HCARB0(M)*EH)+
     &         (DELTA(M)*HCARB0(M)*(1-EH))+
     &         (BETA(M)*BCARB0(M)*(1-EB))+
     &         (BETA(M)*DPMCARB0(M)*(1-EDPM))+
     &         (BETA(M)*RPMCARB0(M)*(1-ERPM))
C
C Calculate the CO2(=CO2()) released due to biological activity
C
        CO2(M)=((1-(ALPHA(M)+BETA(M)))*DPMCARB0(M)*(1-EDPM))+
     &       ((1-(ALPHA(M)+BETA(M)))*RPMCARB0(M)*(1-ERPM))+
     &       ((1-(ALPHA(M)+BETA(M)))*BCARB0(M)*(1-EB))+
     &       ((1-(GAMMA(M)+DELTA(M)))*HCARB0(M)*(1-EH))
C
C Calculate the change in nitrogen in BIO(=BNIT()) and HUM(=HNIT())
C pools from N:C ratio ( = HY)
C
        BNIT(M)=BCARB(M)*HY(NSOIL,M)
        HNIT(M)=HCARB(M)*HY(NSOIL,M)
C
C Calculate the change in N15 in BIO(=BNLAB()) and HUM(=HNLAB())
C pools from N:C ratio ( = HY) amd N15:N ratio
C
        IF(Z(M).GT.(ALPHA(M)+BETA(M))*HY(NSOIL,M))THEN
          BNLAB(M)=(S(M)*HY(NSOIL,M)*BCARB0(M)*EB)+
     &           (S(M)*HY(NSOIL,M)*ALPHA(M)*BCARB0(M)*(1-EB))+
     &           (P(M)*HY(NSOIL,M)*GAMMA(M)*HCARB0(M)*(1-EH))+
     &           (Q(M)*HY(NSOIL,M)*ALPHA(M)*DPMCARB0(M)*(1-EDPM))+
     &           (Q(M)*HY(NSOIL,M)*ALPHA(M)*RPMCARB0(M)*(1-ERPM))
          HNLAB(M)=(P(M)*HY(NSOIL,M)*HCARB0(M)*EH)+
     &           (P(M)*HY(NSOIL,M)*DELTA(M)*HCARB0(M)*(1-EH))+
     &           (S(M)*HY(NSOIL,M)*BETA(M)*BCARB0(M)*(1-EB))+
     &           (Q(M)*HY(NSOIL,M)*BETA(M)*DPMCARB0(M)*(1-EDPM))+
     &           (Q(M)*HY(NSOIL,M)*BETA(M)*RPMCARB0(M)*(1-ERPM))
        ELSE
          U(M)=(SOIL15(M)+AMMN15(M)+
     &      (S(M)*HY(NSOIL,M)*(1-(ALPHA(M)+BETA(M)))*BCARB0(M)*(1-EB))+
     &      (P(M)*HY(NSOIL,M)*(1-(GAMMA(M)+DELTA(M)))*HCARB0(M)*(1-EH))-
     &      (0.5*T(M)*(HY(NSOIL,M)*(ALPHA(M)+BETA(M))-Z(M))*DPMCARB0(M)*
     &        (1-EDPM))-
     &      (0.5*T(M)*(HY(NSOIL,M)*(ALPHA(M)+BETA(M))-Z(M))*RPMCARB0(M)*
     &        (1-ERPM)))/
     &      (SOILN(M)+AMMN(M)+
     &      (0.5*(Z(M)-HY(NSOIL,M)*(ALPHA(M)+BETA(M)))*DPMCARB0(M)*
     &        (1-EDPM))+
     &      (0.5*(Z(M)-HY(NSOIL,M)*(ALPHA(M)+BETA(M)))*RPMCARB0(M)*
     &        (1-ERPM))+
     &      (HY(NSOIL,M)*(1-(ALPHA(M)+BETA(M)))*BCARB0(M)*(1-EB))+
     &      (HY(NSOIL,M)*(1-(GAMMA(M)+DELTA(M)))*HCARB0(M)*(1-EH)))
          BNLAB(M)=(S(M)*HY(NSOIL,M)*BCARB0(M)*EB)+
     &      (S(M)*HY(NSOIL,M)*ALPHA(M)*BCARB0(M)*(1-EB))+
     &      (P(M)*HY(NSOIL,M)*GAMMA(M)*HCARB0(M)*(1-EH))+
     &      ((ALPHA(M)/(ALPHA(M)+BETA(M)))*Q(M)*Z(M)*DPMCARB0(M)*
     &        (1-EDPM))+
     &      ((ALPHA(M)/(ALPHA(M)+BETA(M)))*Q(M)*Z(M)*RPMCARB0(M)*
     &        (1-ERPM))+
     &      ((ALPHA(M)/(ALPHA(M)+BETA(M)))*0.5*(T(M)+U(M))*
     &        (HY(NSOIL,M)*(ALPHA(M)+BETA(M))-Z(M))*
     &      DPMCARB0(M)*(1-EDPM))+
     &      ((ALPHA(M)/(ALPHA(M)+BETA(M)))*0.5*(T(M)+U(M))*
     &        (HY(NSOIL,M)*(ALPHA(M)+BETA(M))-Z(M))*
     &      RPMCARB0(M)*(1-ERPM))
          HNLAB(M)=(P(M)*HY(NSOIL,M)*HCARB0(M)*EH)+
     &      (P(M)*HY(NSOIL,M)*DELTA(M)*HCARB0(M)*(1-EH))+
     &      (S(M)*HY(NSOIL,M)*BETA(M)*BCARB0(M)*(1-EB))+
     &      ((BETA(M)/(ALPHA(M)+BETA(M)))*Q(M)*Z(M)*DPMCARB0(M)*
     &        (1-EDPM))+
     &      ((BETA(M)/(ALPHA(M)+BETA(M)))*Q(M)*Z(M)*RPMCARB0(M)*
     &        (1-ERPM))+
     &      ((BETA(M)/(ALPHA(M)+BETA(M)))*0.5*(T(M)+U(M))*
     &        (HY(NSOIL,M)*(ALPHA(M)+BETA(M))-Z(M))*
     &      DPMCARB0(M)*(1-EDPM))+
     &      ((ALPHA(M)/(ALPHA(M)+BETA(M)))*0.5*(T(M)+U(M))*
     &        (HY(NSOIL,M)*(ALPHA(M)+BETA(M))-Z(M))*
     &      RPMCARB0(M)*(1-ERPM))
        ENDIF
C
C Calculate the net mineralization = DNIT() and DNIT15()
C
        DNIT(M)=(DPMNIT0(M)-DPMNIT(M))+
     &        (RPMNIT0(M)-RPMNIT(M))+
     &        (BNIT0(M)-BNIT(M))+
     &        (HNIT0(M)-HNIT(M))
        DNIT15(M)=(DPMNLAB0(M)-DPMNLAB(M))+
     &          (RPMNLAB0(M)-RPMNLAB(M))+
     &          (BNLAB0(M)-BNLAB(M))+
     &          (HNLAB0(M)-HNLAB(M))
C
C Calculate amount available for immobilization
C of nitrate-N(=SN) and ammonium-N(=AN)
C
        CRITL=CRIT(NSOIL,M)*MAXDEPTH/(50.*MAXLAYER1)
        SN=SOILN(M)-CRITL
        IF(SN.LT.0)SN=0
        AN=AMMN(M)
        IF(AN.LT.0)AN=0
        SN15=SOIL15(M)-CRITL*TN(M)
        IF(SN15.LT.0)SN15=0
        AN15=AMMN15(M)*TA(M)
        IF(AN15.LT.0)AN15=0
C
C If amount of N immobilized exceeds the N available as nitrate and ammonium
C stop the decomposition of RO (i.e. ER=1) and recalculate the pool changes
C
        IF(LIMITED.EQ.1)THEN
	    EDPM=1
	    ERPM=1
	    EB=1
	    EH=1
	  ENDIF
        IF(-DNIT(M).GT.AN+SN)THEN
C
C Calculate increase in fraction lost as CO2, and recalculate
C pool changes
C
	PROPIMM=(-DPMNIT0(M)+DPMNIT0(M)*EDPM-RPMNIT0(M)+
     & RPMNIT0(M)*ERPM-BNIT0(M)+HY(NSOIL,M)*BCARB0(M)*EB+
     & HY(NSOIL,M)*ALPHA(M)*BCARB0(M)-HY(NSOIL,M)*ALPHA(M)*BCARB0(M)*EB+
     & HY(NSOIL,M)*GAMMA(M)*HCARB0(M)-
     & HY(NSOIL,M)*GAMMA(M)*HCARB0(M)*EH-AN-HNIT0(M)+
     & HY(NSOIL,M)*HCARB0(M)*EH+
     & HY(NSOIL,M)*DELTA(M)*HCARB0(M)-HY(NSOIL,M)*DELTA(M)*HCARB0(M)*EH+
     & HY(NSOIL,M)*BETA(M)*BCARB0(M)-
     & HY(NSOIL,M)*BETA(M)*BCARB0(M)*EB-SN)/(HY(NSOIL,M)*
     & (-BETA(M)*DPMCARB0(M)+
     & BETA(M)*DPMCARB0(M)*EDPM-BETA(M)*RPMCARB0(M)+
     & BETA(M)*RPMCARB0(M)*ERPM-
     & ALPHA(M)*DPMCARB0(M)+ALPHA(M)*DPMCARB0(M)*EDPM-
     & ALPHA(M)*RPMCARB0(M)+ALPHA(M)*RPMCARB0(M)*ERPM))
C
C Set microbial efficiency to max in the beginning, Bente, 16.06.06
C then calculate new value
C
	    ALOW=ALPHA(M)
	    CLOW=BETA(M)
          ALOW=ALOW*PROPIMM
	    CLOW=CLOW*PROPIMM
C							
C Changed to repeat code rather then going back to avoid going through 
C the loop many times. Bente, 08.08.06
C
C Adjust C(=RCARB()), N(=RNIT()) and N15(=RNLAB()) in RO pool
C due to biological activity
C
	    DPMCARB(M)=DPMCARB0(M)*EDPM
          RPMCARB(M)=RPMCARB0(M)*ERPM
          DPMNIT(M)=DPMNIT0(M)*EDPM
          RPMNIT(M)=RPMNIT0(M)*ERPM
          DPMNLAB(M)=DPMNLAB0(M)*EDPM
          RPMNLAB(M)=RPMNLAB0(M)*ERPM
C
C Adjust C(=BCARB()) in BIO pool due to biological activity
C
          BCARB(M)=(BCARB0(M)*EB)+
     &         (ALPHA(M)*BCARB0(M)*(1-EB))+
     &         (GAMMA(M)*HCARB0(M)*(1-EH))+
     &         (ALOW*DPMCARB0(M)*(1-EDPM))+
     &         (ALOW*RPMCARB0(M)*(1-ERPM))
C
C Adjust C(=HCARB()) in HUM pool due to biological activity
C
          HCARB(M)=(HCARB0(M)*EH)+
     &         (DELTA(M)*HCARB0(M)*(1-EH))+
     &         (BETA(M)*BCARB0(M)*(1-EB))+
     &         (CLOW*DPMCARB0(M)*(1-EDPM))+
     &         (CLOW*RPMCARB0(M)*(1-ERPM))
C
C Calculate the CO2(=CO2()) released due to biological activity
C
          CO2(M)=((1-(ALOW+CLOW))*DPMCARB0(M)*(1-EDPM))+
     &       ((1-(ALOW+CLOW))*RPMCARB0(M)*(1-ERPM))+
     &       ((1-(ALPHA(M)+BETA(M)))*BCARB0(M)*(1-EB))+
     &       ((1-(GAMMA(M)+DELTA(M)))*HCARB0(M)*(1-EH))
C
C Calculate the change in nitrogen in BIO(=BNIT()) and HUM(=HNIT())
C pools from N:C ratio ( = HY)
C
          BNIT(M)=BCARB(M)*HY(NSOIL,M)
          HNIT(M)=HCARB(M)*HY(NSOIL,M)
C
C Calculate the change in N15 in BIO(=BNLAB()) and HUM(=HNLAB())
C pools from N:C ratio ( = HY) amd N15:N ratio
C This part has not been thoroughly checked, modified model needs 
C more scrutiny before it is used with 15N (Bente,04.08.06)
C
          IF(Z(M).GT.(ALOW+CLOW)*HY(NSOIL,M))THEN
            BNLAB(M)=(S(M)*HY(NSOIL,M)*BCARB0(M)*EB)+
     &           (S(M)*HY(NSOIL,M)*ALPHA(M)*BCARB0(M)*(1-EB))+
     &           (P(M)*HY(NSOIL,M)*GAMMA(M)*HCARB0(M)*(1-EH))+
     &           (Q(M)*HY(NSOIL,M)*ALOW*DPMCARB0(M)*(1-EDPM))+
     &           (Q(M)*HY(NSOIL,M)*ALOW*RPMCARB0(M)*(1-ERPM))
            HNLAB(M)=(P(M)*HY(NSOIL,M)*HCARB0(M)*EH)+
     &           (P(M)*HY(NSOIL,M)*DELTA(M)*HCARB0(M)*(1-EH))+
     &           (S(M)*HY(NSOIL,M)*BETA(M)*BCARB0(M)*(1-EB))+
     &           (Q(M)*HY(NSOIL,M)*CLOW*DPMCARB0(M)*(1-EDPM))+
     &           (Q(M)*HY(NSOIL,M)*CLOW*RPMCARB0(M)*(1-ERPM))
          ELSE
            U(M)=(SOIL15(M)+AMMN15(M)+
     &      (S(M)*HY(NSOIL,M)*(1-(ALPHA(M)+BETA(M)))*BCARB0(M)*(1-EB))+
     &      (P(M)*HY(NSOIL,M)*(1-(GAMMA(M)+DELTA(M)))*HCARB0(M)*(1-EH))-
     &      (0.5*T(M)*(HY(NSOIL,M)*(ALOW+CLOW)-Z(M))*DPMCARB0(M)*
     &        (1-EDPM))-
     &      (0.5*T(M)*(HY(NSOIL,M)*(ALOW+CLOW)-Z(M))*RPMCARB0(M)*
     &        (1-ERPM)))/
     &      (SOILN(M)+AMMN(M)+
     &      (0.5*(Z(M)-HY(NSOIL,M)*(ALOW+CLOW))*DPMCARB0(M)*(1-EDPM))+
     &      (0.5*(Z(M)-HY(NSOIL,M)*(ALOW+CLOW))*RPMCARB0(M)*(1-ERPM))+
     &      (HY(NSOIL,M)*(1-(ALPHA(M)+BETA(M)))*BCARB0(M)*(1-EB))+
     &      (HY(NSOIL,M)*(1-(GAMMA(M)+DELTA(M)))*HCARB0(M)*(1-EH)))
           BNLAB(M)=(S(M)*HY(NSOIL,M)*BCARB0(M)*EB)+
     &      (S(M)*HY(NSOIL,M)*ALPHA(M)*BCARB0(M)*(1-EB))+
     &      (P(M)*HY(NSOIL,M)*GAMMA(M)*HCARB0(M)*(1-EH))+
     &      ((ALOW/(ALOW+CLOW))*Q(M)*Z(M)*DPMCARB0(M)*(1-EDPM))+
     &      ((ALOW/(ALOW+CLOW))*Q(M)*Z(M)*RPMCARB0(M)*(1-ERPM))+
     &      ((ALOW/(ALOW+CLOW))*0.5*(T(M)+U(M))*
     &	  (HY(NSOIL,M)*(ALOW+CLOW)-Z(M))*
     &      DPMCARB0(M)*(1-EDPM))+
     &      ((ALOW/(ALOW+CLOW))*0.5*(T(M)+U(M))*
     &	  (HY(NSOIL,M)*(ALOW+CLOW)-Z(M))*
     &      RPMCARB0(M)*(1-ERPM))
           HNLAB(M)=(P(M)*HY(NSOIL,M)*HCARB0(M)*EH)+
     &      (P(M)*HY(NSOIL,M)*DELTA(M)*HCARB0(M)*(1-EH))+
     &      (S(M)*HY(NSOIL,M)*BETA(M)*BCARB0(M)*(1-EB))+
     &      ((CLOW/(ALOW+CLOW))*Q(M)*Z(M)*DPMCARB0(M)*(1-EDPM))+
     &      ((CLOW/(ALOW+CLOW))*Q(M)*Z(M)*RPMCARB0(M)*(1-ERPM))+
     &      ((CLOW/(ALOW+CLOW))*0.5*(T(M)+U(M))*
     &	  (HY(NSOIL,M)*(ALOW+CLOW)-Z(M))*DPMCARB0(M)*(1-EDPM))+
     &      ((CLOW/(ALOW+CLOW))*0.5*(T(M)+U(M))
     &	  *(HY(NSOIL,M)*(ALOW+CLOW)-Z(M))*
     &      RPMCARB0(M)*(1-ERPM))
        ENDIF
C
C Calculate the net mineralization = DNIT() and DNIT15()
C
        DNIT(M)=(DPMNIT0(M)-DPMNIT(M))+
     &        (RPMNIT0(M)-RPMNIT(M))+
     &        (BNIT0(M)-BNIT(M))+
     &        (HNIT0(M)-HNIT(M))
        DNIT15(M)=(DPMNLAB0(M)-DPMNLAB(M))+
     &          (RPMNLAB0(M)-RPMNLAB(M))+
     &          (BNLAB0(M)-BNLAB(M))+
     &          (HNLAB0(M)-HNLAB(M))
C
C Calculate amount available for immobilization
C of nitrate-N(=SN) and ammonium-N(=AN)
C
        CRITL=CRIT(NSOIL,M)*MAXDEPTH/(50.*MAXLAYER1)
        SN=SOILN(M)-CRITL
        IF(SN.LT.0)SN=0
        AN=AMMN(M)-CRITL
        IF(AN.LT.0)AN=0
        SN15=SOIL15(M)-CRITL*TN(M)
        IF(SN15.LT.0)SN15=0
        AN15=AMMN15(M)-CRITL*TA(M)
        IF(AN15.LT.0)AN15=0
C
C If amount of N immobilized exceeds the N available as ammonium take
C the excess N from nitrate pool
C
        ELSEIF(-DNIT(M).GT.AN)THEN
          SOILN(M)=SOILN(M)+AN+DNIT(M)
          AMMN(M)=CRITL
          SNIT=-AN-DNIT(M)
          ANIT=AN
C
C If amount of N immobilized does not exceed N available as ammonium take
C all the required N from ammonium pool
C
        ELSE
          AMMN(M)=AMMN(M)+DNIT(M)
          ANIT=-DNIT(M)
        END IF
C
C Calculate changes in N15 in nitrate and ammonium pools
C
        IF(-DNIT15(M).GT.(SN15+AN15))THEN								
          XN=SN15+AN15+DNIT15(M)
          DNIT15(M)=SN15+AN15
          BNLAB(M)=BNLAB(M)+XN*(ALPHA(M)/(ALPHA(M)+BETA(M)))
          HNLAB(M)=HNLAB(M)+XN*(ALPHA(M)/(ALPHA(M)+BETA(M)))
          AMMN15(M)=AMMN15(M)-AN15
          SOIL15(M)=SOIL15(M)-SN15
	  ELSEIF(-DNIT15(M).GT.AN15)THEN
          SOIL15(M)=SOIL15(M)+AN15+DNIT15(M)
          AMMN15(M)=CRITL*TA(M)
        ELSE
          AMMN15(M)=AMMN15(M)+DNIT15(M)
        END IF
C
C Calculate the ratio of mineral N15:N
C
        T(M)=(AMMN15(M)+SOIL15(M))/(AMMN(M)+SOILN(M))
C
C Reset RO, HUM and BIO pools
C
        DPMCARB0(M)=DPMCARB(M)
        RPMCARB0(M)=RPMCARB(M)
        BCARB0(M)=BCARB(M)
        HCARB0(M)=HCARB(M)
        DPMNIT0(M)=DPMNIT(M)
        RPMNIT0(M)=RPMNIT(M)
        BNIT0(M)=BNIT(M)
        HNIT0(M)=HNIT(M)
        DPMNLAB0(M)=DPMNLAB(M)
        RPMNLAB0(M)=RPMNLAB(M)
        BNLAB0(M)=BNLAB(M)
        HNLAB0(M)=HNLAB(M)
C
C Go back and calculate pool changes for next layer
C
   25 CONTINUE
C
C Total N mineralization in 0-50cm layer = TDNIT
C Total N15 mineralization in 0-50cm layer = T15
C
      TDNIT=0
	T15=0
	DO 26 M=1,MAXLAYER
	  TDNIT=TDNIT+DNIT(M)
	  T15=T15+DNIT15(M)
26    CONTINUE
C
C Leave MINER1
C
      RETURN
      END
C
C-------------------------------------------------------------

C-------------------------------------------------------------
C
      SUBROUTINE MINER1_BRADBURY(NSOIL,ICOVER,ITFUNC,IMFUNC,
     &                    PI_C,PI_N,PI_N15,DRRAT,								
     &                    SOILW,WMAX,WSAT,SOILTEMP,
     &                    ALPHA,BETA,GAMMA,DELTA,HY,					
     &                    BRATE,HRATE,DPMRATE,RPMRATE,				
     &                    AMMN,AMMN15,SOILN,SOIL15,CRIT,
     &                    CO2, CO2DPM, CO2RPM, CO2BIO, CO2HUM,
     &                    DNIT,DNIT15,TDNIT,T15,
     &                    BCARB0,BNIT0,BNLAB,
     &                    HCARB0,HNIT0,HNLAB,
     &                    DPMCARB0,DPMNIT0,DPMNLAB0,
     &                    RPMCARB0,RPMNIT0,RPMNLAB0,
     &                    PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &                    DPMCTON,RPMCTON,
     &					DPMCARB,RPMCARB,BCARB,HCARB,
     &					WRATEDM,TRATEM,MEASLAY)
C
C Subroutine to calculate mineralisation and immobilisation
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXSOIL				! Max.no.of soil types
      PARAMETER(MAXSOIL=50)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	INTEGER MAXDEPTH			! Maximum depth of the soil profile
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/

	REAL ALOW					! Limited parameter for biomass decomposition
	REAL AN						! Amount ammonium avail.for immob.(kgN/ha)
	REAL AN15					! Amount ammonium-15 avail.for immob.(kgN/ha)
	REAL ANIT					! Ammonium immobilised (kgN/ha)
	REAL BCARB(MAXLAYER)		! C in soil biomass at end (kgC/ha/layer)
      REAL BNIT(MAXLAYER)			! N in soil biomass at end (kgN/ha/layer)
	REAL BNLAB(MAXLAYER)		! N15 in soil humus at end (kgN15/ha/layer)
	REAL CARB(MAXLAYER)			! Carbon input to layer (kgC/ha)
	REAL CLOW					! Limited parameter for humus decomposition
	REAL CRITL					! Critical minimum N in the layer (kgN/ha)
	REAL DPMCARB(MAXLAYER)		! C in decomposable PM at end (kgC/ha/layer)
	REAL DPMNIT(MAXLAYER)		! N in decomposable PM at end (kgN/ha/layer)
	REAL DPMNLAB(MAXLAYER)		! N15 in decomposable PM at end (kgN15/ha/layer)
      REAL DPMNRAT				! Ratio of plant matter N in DPM
	REAL*8 EB					! Proportion of biomass decomposed 
	REAL*8 EDPM					! Proportion of DPM decomposed 
	REAL*8 EH					! Proportion of humus decomposed 
	REAL*8 ERPM					! Proportion of RPM decomposed 
	REAL HCARB(MAXLAYER)		! C in soil humus at end (kgC/ha/layer)
	REAL HNIT(MAXLAYER)			! N in soil humus at end (kgN/ha/layer)
      REAL HNLAB(MAXLAYER)		! N15 in soil biomass at end (kgN15/ha/layer)
	INTEGER IDEPTH				! Depth of this layer
	INTEGER ISAVE				! Save = 1, retrieve = 0
	INTEGER LIMITED				! Marker to indicate immobilisation is limited
	INTEGER M					! Layer counter
	REAL MODALPHA				! Proportion of decomposed BIO that goes to BIO
	REAL MODBETA				! Proportion of decomposed BIO that goes to HUM
	REAL MODDELTA				! Proportion of decomposed HUM that goes to HUM
	REAL MODGAMMA				! Proportion of decomposed HUM that goes to BIO
	REAL MODBPART				! Efficiency of decomposition (fraction retained in soil)
	REAL MODBPROP				! BIO/HUM
	REAL P(MAXLAYER)			! N15/N of humus
	REAL PH						! Current pH of the soil
	REAL PHP1					! pH below which rate is zero
	REAL PHP2					! pH above which rate is optimum
	REAL PROPIMM				! Proportion of mineral N available for immob.												
	REAL Q(MAXLAYER)			! N15/N of PM
	REAL RADD(MAXLAYER)			! Nitrogen input to layer (kgN/ha)
	REAL RADD15(MAXLAYER)		! N15 input to layer (kgN15/ha)
	REAL RPMCARB(MAXLAYER)		! C in resistant PM at end (kgC/ha/layer)
	REAL RPMNIT(MAXLAYER)		! N in resistant PM at end (kgN/ha/layer)
	REAL RPMNLAB(MAXLAYER)		! N15 in resistant PM at end (kgN15/ha/layer)
      REAL RPMNRAT				! Ratio of plant matter N in RPM
	REAL S(MAXLAYER)			! N15/N of biomass
	REAL SN						! Amount nitrate avail.for immob.(kgN/ha)
	REAL SN15					! Amount nitrate-15 avail.for immob.(kgN/ha)
	REAL SNIT					! Nitrate immobilised (kgN/ha)
	REAL SPL(MAXLAYER)			! Split of input into the layer
	REAL T(MAXLAYER)			! N15/N of mineral N
	REAL TA(MAXLAYER)			! N15/N of ammonium
	REAL TN(MAXLAYER)			! N15/N of nitrate
	REAL U(MAXLAYER)			! Variable used in calculation of 15N loss
	REAL XN						! Variable used in calc.of changes in pools
	REAL Z(MAXLAYER)			! N:C of PM
C
C Variables passed to/from calling subroutine
C
	INTEGER MEASLAY					! Layer that soil is measured to     
C
C ...Weather factors
C
	REAL SOILTEMP(MAXLAYER)		! IN: Soil temperature (deg.C/timestep)
C
C ...Model descriptors
C
	INTEGER ITFUNC				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IMFUNC				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
C
C ... Soil factors
C
	REAL ALPHA(MAXLAYER)		! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
  	REAL BCARB0(MAXLAYER)		! IN:C in soil biomass (kgC/ha/layer)
	REAL BETA(MAXLAYER)			! IN/OUT: Prop.HUM produced on BIO decompn
      REAL BNIT0(MAXLAYER)		! IN:N in soil biomass (kgN/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN:N15 in soil humus (kgN15/ha/layer)
	REAL BRATE(MAXLAYER)		! IN/OUT: Rate constant for HUM decompn
	REAL CO2(MAXLAYER)			! IN:CO2 emitted (kgC/ha/layer)
	REAL CO2DPM(MAXLAYER)       ! OUT:CO2 emitted from the DPM pool [kgC/ha/layer]
	REAL CO2RPM(MAXLAYER)       ! OUT:CO2 emitted from the RPM pool [kgC/ha/layer]
	REAL CO2BIO(MAXLAYER)       ! OUT:CO2 emitted from the BIO pool [kgC/ha/layer]
	REAL CO2HUM(MAXLAYER)       ! OUT:CO2 emitted from the HUM pool [kgC/ha/layer]
	REAL CRIT(MAXSOIL,MAXLAYER)	! IN/OUT: Minimum nitrate level (Sozanska)
      REAL CRRATE					! IN:crop modifying factor (proportion)
	REAL DELTA(MAXLAYER)		! IN/OUT: Prop.HUM produced on HUM decompn 
	REAL DNIT(MAXLAYER)			! IN:Net Mineralised N (kgN/ha/lay/timestep)
	REAL DNIT15(MAXLAYER)		! IN:Net Mineralised N15 (kgN15/ha/lay/time)
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
      REAL DPMCTON(MAXLAYER)		! IN:DPM C:N ratio
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL DPMRATE(MAXLAYER)		! IN/OUT: Rate constant for DPM decompn
	REAL DRRAT					! IN:DPM:RPM ratio
	REAL GAMMA(MAXLAYER)		! IN/OUT: Prop.BIO produced on HUM decompn 
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN:N in soil humus (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN:N15 in soil biomass (kgN15/ha/layer)
	REAL HRATE(MAXLAYER)		! IN/OUT: Rate constant for BIO decompn
	REAL HY(MAXSOIL,MAXLAYER)	! IN/OUT:Stable N:C ratio of BIO & HUM pools
  	REAL ICFACTOR(MAXLAYER)		! IN:Adjustment in rate needed to	
								!    achieve measured NPP and TOC	
      INTEGER ICOVER				! IN:Code for crop cover: 1=covered; 2=bare
	INTEGER NSOIL				! IN:Soil code number
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)	! IN/PUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)	! IN/PUT: pH above which decomp.max.
      REAL PHRATE					! IN:pH rate modifier (proportion)
      REAL PI_C(MAXLAYER)			! IN/OUT: Plant input C to soil (kgC/ha/step)
      REAL PI_N(MAXLAYER)			! IN/OUT: Plant input N to soil (kgC/ha/step)
      REAL PI_N15(MAXLAYER)		! IN/OUT: Plant input N15 to soil(kgC/ha/step)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
      REAL RPMCTON(MAXLAYER)		! IN:RPM C:N ratio
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL RPMRATE(MAXLAYER)		! IN/OUT: Rate constant for RPM decompn
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL T15					! IN:Net Mineralised N (kgN/ha/layer)
	REAL TDNIT					! IN:Net Mineralised N (kgN/ha/layer)
	REAL TRATE					! IN:temperature rate modifying factor (prop.)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WRATED					! IN:Moisture rate modifyer for DPM&BIO(prop)
	REAL WRATER					! IN:Moisture rate modifyer for RPM&HUM(prop)
	REAL WRATEDM
	REAL TRATEM
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
C
C Save Variable Values
C
      SAVE
C
C  Calculate the amount added from plant debris this week of
C  C(=CARB(M)), N(=RADD(M)) and N15(=RADD15(M)) to layer M
C
      DO 24 M=1,MEASLAY
	  IF(PI_C(M).LT.0)PI_C(M)=0
	  IF(PI_N(M).LT.0)PI_N(M)=0
	  IF(PI_N15(M).LT.0)PI_N15(M)=0
        CARB(M)=PI_C(M)
	 IF(PI_C(M).GT.0)THEN
	 ICOVER=1
	 ELSE
	 ICOVER=0
	 ENDIF
        RADD(M)=PI_N(M)
        RADD15(M)=PI_N15(M)
   24 CONTINUE
C
C For each layer...
C
	WRATEDM=0
	TRATEM=0

      DO 25 M=1,MEASLAY
        CALL MODFACTS_MINER(SOILW(M),WMAX(M),WSAT(M),
     &                WRATED,WRATER,
     &                SOILTEMP(M),TRATE,
     &			    PHARRAY(NSOIL,M),PHP1ARRAY(NSOIL,M),
     &                PHP2ARRAY(NSOIL,M),PHRATE,
     &                ICOVER,CRRATE,ITFUNC,IMFUNC)
C
		IF (M.LT.11) THEN
		WRATEDM=WRATEDM+WRATED
		TRATEM=TRATEM+TRATE
		ENDIF
C
C
C Calculate the exponential factors for decomposition of
C the RO(=ER), BIO(=EB) and HUM(=EH) pools
C
        EDPM=EXP(-DPMRATE(M)*WRATED*TRATE*CRRATE*PHRATE*ICFACTOR(M))
	  ERPM=EXP(-RPMRATE(M)*WRATER*TRATE*CRRATE*PHRATE*ICFACTOR(M))
        EB=EXP(-BRATE(M)*WRATED*TRATE*CRRATE*PHRATE*ICFACTOR(M))
        EH=EXP(-HRATE(M)*WRATER*TRATE*CRRATE*PHRATE*ICFACTOR(M))
C
C Add to the RO pool in the layer
C
        DPMNRAT=((DRRAT/(1+DRRAT))/DPMCTON(M))
	  RPMNRAT=((1/(1+DRRAT))/RPMCTON(M))

        IF(CARB(M).GT.0)THEN
!	    DPMNRAT=CARB(M)*((DRRAT/(1+DRRAT))/DPMCTON(M))
!	    RPMNRAT=CARB(M)*((1/(1+DRRAT))/RPMCTON(M))
           DPMCARB0(M)=DPMCARB0(M)+CARB(M)*(DRRAT/(1+DRRAT))
           RPMCARB0(M)=RPMCARB0(M)+CARB(M)*(1/(1+DRRAT))
           DPMNIT0(M)=DPMNIT0(M)+(RADD(M)*(DPMNRAT/(DPMNRAT+RPMNRAT)))
           RPMNIT0(M)=RPMNIT0(M)+(RADD(M)*(RPMNRAT/(DPMNRAT+RPMNRAT)))
           DPMNLAB0(M)=DPMNLAB0(M)+
     &                 (RADD15(M)*(DPMNRAT/(DPMNRAT+RPMNRAT)))
           RPMNLAB0(M)=RPMNLAB0(M)+
     &                 (RADD15(M)*(RPMNRAT/(DPMNRAT+RPMNRAT)))
        ENDIF
	  IF((DPMCARB0(M)+RPMCARB0(M)).GT.0)THEN
          Z(M)=(DPMNIT0(M)+RPMNIT0(M))/(DPMCARB0(M)+RPMCARB0(M))
	  ELSE
	    Z(M)=0
	  ENDIF
C
C Calc N15/N for NH4(=TA()), NO3(=TN()), MinN(=T()), RO(=Q()), BIO(=S()) and HUM(=P())
C
        TA(M)=AMMN15(M)/AMMN(M)
        TN(M)=SOIL15(M)/SOILN(M)
        T(M)=(AMMN15(M)+SOIL15(M))/(AMMN(M)+SOILN(M))
	  Q(M)=0
	  S(M)=0
        P(M)=0
	  IF((DPMNIT0(M)+RPMNIT0(M)).GT.0)
     &    Q(M)=(DPMNLAB0(M)+RPMNLAB0(M))/(DPMNIT0(M)+RPMNIT0(M))
	  IF(BNIT0(M).GT.0)S(M)=BNLAB0(M)/BNIT0(M)
	  IF(HNIT0(M).GT.0)P(M)=HNLAB0(M)/HNIT0(M)
C
C Calculation of Pool Changes...
C
        LIMITED=0
	  MODALPHA=ALPHA(M)
	  MODBETA=BETA(M)
	  MODDELTA=DELTA(M)
	  MODGAMMA=GAMMA(M)
  100   CONTINUE
C
C Adjust C(=RCARB()), N(=RNIT()) and N15(=RNLAB()) in RO pool
C due to biological activity
C
        DPMCARB(M)=DPMCARB0(M)*EDPM
        RPMCARB(M)=RPMCARB0(M)*ERPM
        DPMNIT(M)=DPMNIT0(M)*EDPM
        RPMNIT(M)=RPMNIT0(M)*ERPM
        DPMNLAB(M)=DPMNLAB0(M)*EDPM
        RPMNLAB(M)=RPMNLAB0(M)*ERPM
C
C Adjust C(=BCARB()) in BIO pool due to biological activity
C
        BCARB(M)=(BCARB0(M)*EB)+
     &         (MODALPHA*BCARB0(M)*(1-EB))+
     &         (MODGAMMA*HCARB0(M)*(1-EH))+
     &         (MODALPHA*DPMCARB0(M)*(1-EDPM))+
     &         (MODALPHA*RPMCARB0(M)*(1-ERPM))
C
C Adjust C(=HCARB()) in HUM pool due to biological activity
C
        HCARB(M)=(HCARB0(M)*EH)+
     &         (MODDELTA*HCARB0(M)*(1-EH))+
     &         (MODBETA*BCARB0(M)*(1-EB))+
     &         (MODBETA*DPMCARB0(M)*(1-EDPM))+
     &         (MODBETA*RPMCARB0(M)*(1-ERPM))
C
C Calculate the CO2(=CO2()) released due to biological activity
C
c >>> Change by Mark Richards required for methane_richards()
c        CO2(M)=((1-(MODALPHA+MODBETA))*DPMCARB0(M)*(1-EDPM))+
c     &       ((1-(MODALPHA+MODBETA))*RPMCARB0(M)*(1-ERPM))+
c     &       ((1-(MODALPHA+MODBETA))*BCARB0(M)*(1-EB))+
c     &       ((1-(MODGAMMA+MODDELTA))*HCARB0(M)*(1-EH))
       CO2DPM(M) = (1-(MODALPHA+MODBETA))*DPMCARB0(M)*(1-EDPM)
       CO2RPM(M) = (1-(MODALPHA+MODBETA))*RPMCARB0(M)*(1-ERPM)
       CO2BIO(M) = (1-(MODALPHA+MODBETA))*BCARB0(M)*(1-EB)
       CO2HUM(M) = (1-(MODGAMMA+MODDELTA))*HCARB0(M)*(1-EH)
       CO2(M) = CO2DPM(M) + CO2RPM(M) + CO2BIO(M) + CO2HUM(M)
c <<< EOF change by Mark Richards
C
C Calculate the change in nitrogen in BIO(=BNIT()) and HUM(=HNIT())
C pools from N:C ratio ( = HY)
C
        BNIT(M)=BCARB(M)*HY(NSOIL,M)
        HNIT(M)=HCARB(M)*HY(NSOIL,M)
C
C Calculate the change in N15 in BIO(=BNLAB()) and HUM(=HNLAB())
C pools from N:C ratio ( = HY) amd N15:N ratio
C
        IF(Z(M).GT.(MODALPHA+MODBETA)*HY(NSOIL,M))THEN
          BNLAB(M)=(S(M)*HY(NSOIL,M)*BCARB0(M)*EB)+
     &           (S(M)*HY(NSOIL,M)*MODALPHA*BCARB0(M)*(1-EB))+
     &           (P(M)*HY(NSOIL,M)*MODGAMMA*HCARB0(M)*(1-EH))+
     &           (Q(M)*HY(NSOIL,M)*MODALPHA*DPMCARB0(M)*(1-EDPM))+
     &           (Q(M)*HY(NSOIL,M)*MODALPHA*RPMCARB0(M)*(1-ERPM))
          HNLAB(M)=(P(M)*HY(NSOIL,M)*HCARB0(M)*EH)+
     &           (P(M)*HY(NSOIL,M)*MODDELTA*HCARB0(M)*(1-EH))+
     &           (S(M)*HY(NSOIL,M)*MODBETA*BCARB0(M)*(1-EB))+
     &           (Q(M)*HY(NSOIL,M)*MODBETA*DPMCARB0(M)*(1-EDPM))+
     &           (Q(M)*HY(NSOIL,M)*MODBETA*RPMCARB0(M)*(1-ERPM))
        ELSE
          U(M)=(SOIL15(M)+AMMN15(M)+
     &     (S(M)*HY(NSOIL,M)*(1-(MODALPHA+MODBETA))*BCARB0(M)
     &      *(1-EB))+
     &     (P(M)*HY(NSOIL,M)*(1-(MODGAMMA+MODDELTA))*HCARB0(M)*
     &     (1-EH))-(0.5*T(M)*(HY(NSOIL,M)*(MODALPHA+MODBETA)-Z(M))
     &      *DPMCARB0(M)*(1-EDPM))-
     &     (0.5*T(M)*(HY(NSOIL,M)*(MODALPHA+MODBETA)-Z(M))*
     &      RPMCARB0(M)*(1-ERPM)))/
     &     (SOILN(M)+AMMN(M)+
     &     (0.5*(Z(M)-HY(NSOIL,M)*(MODALPHA+MODBETA))*DPMCARB0(M)*
     &     (1-EDPM))+
     &     (0.5*(Z(M)-HY(NSOIL,M)*(MODALPHA+MODBETA))*RPMCARB0(M)*
     &     (1-ERPM))+
     &     (HY(NSOIL,M)*(1-(MODALPHA+MODBETA))*BCARB0(M)*(1-EB))+
     &     (HY(NSOIL,M)*(1-(MODGAMMA+MODDELTA))*HCARB0(M)*(1-EH)))
          
		BNLAB(M)=(S(M)*HY(NSOIL,M)*BCARB0(M)*EB)+
     &      (S(M)*HY(NSOIL,M)*MODALPHA*BCARB0(M)*(1-EB))+
     &      (P(M)*HY(NSOIL,M)*MODGAMMA*HCARB0(M)*(1-EH))+
     &      ((MODALPHA/(MODALPHA+MODBETA))*Q(M)*Z(M)*
     &       DPMCARB0(M)*(1-EDPM))+
     &      ((MODALPHA/(MODALPHA+MODBETA))*Q(M)*Z(M)*
     &       RPMCARB0(M)*(1-ERPM))+
     &      ((MODALPHA/(MODALPHA+MODBETA))*0.5*(T(M)+U(M))*
     &        (HY(NSOIL,M)*(MODALPHA+MODBETA)-Z(M))*
     &      DPMCARB0(M)*(1-EDPM))+
     &      ((MODALPHA/(MODALPHA+MODBETA))*0.5*(T(M)+U(M))*
     &        (HY(NSOIL,M)*(MODALPHA+MODBETA)-Z(M))*
     &      RPMCARB0(M)*(1-ERPM))
          HNLAB(M)=(P(M)*HY(NSOIL,M)*HCARB0(M)*EH)+
     &      (P(M)*HY(NSOIL,M)*MODDELTA*HCARB0(M)*(1-EH))+
     &      (S(M)*HY(NSOIL,M)*MODBETA*BCARB0(M)*(1-EB))+
     &      ((MODBETA/(MODALPHA+MODBETA))*Q(M)*Z(M)*
     &       DPMCARB0(M)*(1-EDPM))+
     &      ((MODBETA/(MODALPHA+MODBETA))*Q(M)*Z(M)*
     &       RPMCARB0(M)*(1-ERPM))+
     &      ((MODBETA/(MODALPHA+MODBETA))*0.5*(T(M)+U(M))*
     &        (HY(NSOIL,M)*(MODALPHA+MODBETA)-Z(M))*
     &      DPMCARB0(M)*(1-EDPM))+
     &      ((MODALPHA/(MODALPHA+MODBETA))*0.5*(T(M)+U(M))*
     &        (HY(NSOIL,M)*(MODALPHA+MODBETA)-Z(M))*
     &      RPMCARB0(M)*(1-ERPM))
        ENDIF
C
C Calculate the net mineralization = DNIT() and DNIT15()
C
        DNIT(M)=(DPMNIT0(M)-DPMNIT(M))+
     &        (RPMNIT0(M)-RPMNIT(M))+
     &        (BNIT0(M)-BNIT(M))+
     &        (HNIT0(M)-HNIT(M))
        DNIT15(M)=(DPMNLAB0(M)-DPMNLAB(M))+
     &          (RPMNLAB0(M)-RPMNLAB(M))+
     &          (BNLAB0(M)-BNLAB(M))+
     &          (HNLAB0(M)-HNLAB(M))
C
C Calculate amount available for immobilization
C of nitrate-N(=SN) and ammonium-N(=AN)
C
        CRITL=CRIT(NSOIL,M)*MAXDEPTH/(50.*MAXLAYER1)
        SN=SOILN(M)-CRITL
        IF(SN.LT.0)SN=0
        AN=AMMN(M)-CRITL
        IF(AN.LT.0)AN=0
        SN15=SOIL15(M)-CRITL*TN(M)
        IF(SN15.LT.0)SN15=0
        AN15=AMMN15(M)-CRITL*TA(M)
        IF(AN15.LT.0)AN15=0
C
C If amount of N immobilized exceeds the N available as nitrate and ammonium
C stop the decomposition of RO (i.e. ER=1) and recalculate the pool changes
C
        IF(-DNIT(M).GT.(AN+SN))THEN
	    IF(LIMITED.EQ.0)THEN
c            IF(AN+SN.GT.0)THEN
c 	        PROPIMM=(AN+SN)/(-1*DNIT(M))									
c	        EDPM = EDPM*(1-PROPIMM)
c	        ERPM = ERPM*(1-PROPIMM)
c	        EB = EB*(1-PROPIMM)
c	        EH = EB*(1-PROPIMM)
c            ELSE
c              EDPM=1
c 	        ERPM=1
c 	        EB=1
c 	        EH=1
c	      ENDIF
********230909         
            MODBPROP=ALPHA(M)/BETA(M)
C            MODBPART=(AN+SN+
C     &                (DPMNIT0(M)*(1-EDPM))+
C     &                (RPMNIT0(M)*(1-ERPM))+
C     &                (HY(NSOIL,M)*(BCARB0(M)*(1-EB)))+
C     &                (HY(NSOIL,M)*(HCARB0(M)*(1-EH))))/
C     &               (HY(NSOIL,M)*(DPMCARB0(M)*(1-EDPM))+
C     &                HY(NSOIL,M)*(RPMCARB0(M)*(1-ERPM))+
C     &                HY(NSOIL,M)*(BCARB0(M)*(1-EB))+
C     &                HY(NSOIL,M)*(HCARB0(M)*(1-EH)))

	      MODBPART=(MODALPHA+MODBETA)*0.9 
		  IF(MODBPART.LE.0.0001.OR.AN+SN.LT.0.0001)THEN
              EDPM=1
 	        ERPM=1
 	        EB=1
 	        EH=1
	        LIMITED=1
            ENDIF
	      MODBETA=MODBPART/(1+MODBPROP)
	      MODALPHA=MODBPART-MODBETA
	      MODDELTA=MODBPART/(1+MODBPROP)
	      MODGAMMA=MODBPART-MODBETA
		  DPMNIT(M)=DPMNIT0(M)
		  RPMNIT(M)=RPMNIT0(M)
		  BNIT(M)=BNIT0(M)
		  HNIT(M)=HNIT0(M)
		  DPMNLAB(M)=DPMNLAB0(M)
		  RPMNLAB(M)=RPMNLAB0(M)
		  BNLAB(M)=BNLAB0(M)
		  HNLAB(M)=HNLAB0(M)
		  DPMCARB(M)=DPMCARB0(M)
		  RPMCARB(M)=RPMCARB0(M)
		  BCARB(M)=BCARB0(M)
		  HCARB(M)=HCARB0(M)
	      DNIT(M)=0
	      DNIT15(M)=0
            GOTO 100														
	    ENDIF
C	      
C
C If amount of N immobilized exceeds the N available as ammonium take
C the excess N from nitrate pool	erpm	0.9985460

C
        ELSEIF(-DNIT(M).GT.AN)THEN
          SOILN(M)=SOILN(M)+AN+DNIT(M)
          AMMN(M)=CRITL
          SNIT=-AN-DNIT(M)
          ANIT=AN
C
C If amount of N immobilized does not exceed N available as ammonium take
C all the required N from ammonium pool
C
        ELSE
          AMMN(M)=AMMN(M)+DNIT(M)
          ANIT=-DNIT(M)
        END IF
	  LIMITED=0
C
C Calculate changes in N15 in nitrate and ammonium pools
C
        IF(-DNIT15(M).GT.(SN15+AN15))THEN								
          XN=SN15+AN15+DNIT15(M)
          DNIT15(M)=SN15+AN15
          BNLAB(M)=BNLAB(M)+XN*(ALPHA(M)/(ALPHA(M)+BETA(M)))
          HNLAB(M)=HNLAB(M)+XN*(ALPHA(M)/(ALPHA(M)+BETA(M)))
          AMMN15(M)=AMMN15(M)-AN15
          SOIL15(M)=SOIL15(M)-SN15
	  ELSEIF(-DNIT15(M).GT.AN15)THEN
          SOIL15(M)=SOIL15(M)+AN15+DNIT15(M)
          AMMN15(M)=CRITL*TA(M)
        ELSE
          AMMN15(M)=AMMN15(M)+DNIT15(M)
        END IF
C
C Calculate the ratio of mineral N15:N
C
        T(M)=(AMMN15(M)+SOIL15(M))/(AMMN(M)+SOILN(M))
C
C Reset RO, HUM and BIO pools
C
        DPMCARB0(M)=DPMCARB(M)
        RPMCARB0(M)=RPMCARB(M)
        BCARB0(M)=BCARB(M)
        HCARB0(M)=HCARB(M)
        DPMNIT0(M)=DPMNIT(M)
        RPMNIT0(M)=RPMNIT(M)
        BNIT0(M)=BNIT(M)
        HNIT0(M)=HNIT(M)
        DPMNLAB0(M)=DPMNLAB(M)
        RPMNLAB0(M)=RPMNLAB(M)
        BNLAB0(M)=BNLAB(M)
        HNLAB0(M)=HNLAB(M)
C
C Go back and calculate pool changes for next layer
C
25    CONTINUE
C
C Total N mineralization in 0-50cm layer = TDNIT
C Total N15 mineralization in 0-50cm layer = T15
C
      TDNIT=0
	T15=0
	DO 26 M=1,MAXLAYER
	  TDNIT=TDNIT+DNIT(M)
	  T15=T15+DNIT15(M)
26    CONTINUE
C
C Leave MINER1
C
       RETURN
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE MODFACTS_ANDEC(LAYERSOILW,LAYERWMAX,LAYERWSAT,WRATE,
     &                          THISTEMP,TRATE,
     &			              PH,PHRATE)
C
C Subroutine to calculate modifying factors for rate of anaerobic decomposition
C
	IMPLICIT NONE
C
C Variables local to this subroutine
C ...Profile descriptors
C
      INTEGER MAXLAYER1	! Maximum number of layers in the profile
	INTEGER MAXDEPTH	! Maximum depth of the profile
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
C
C ...Parameters 
C
      REAL C1				! Constant for calculation of moisture modifier
      REAL C3				! Constant for calculation of pH modifier
      REAL C4				! Constant for calculation of pH modifier
C
C Variables passed to /from this subroutine
C ... Weather variables
C
      REAL THISTEMP		! IN:Soil temperature in this layer(deg.C)
C
C ...Soil descriptors
C
	REAL LAYERSOILW		! IN:Available soil water (mm/layer)
	REAL LAYERWMAX		! IN:Maximum available soil water (mm/layer)
	REAL LAYERWSAT		! IN:Available water at saturation (mm/layer)
	REAL PH				! IN:Current pH of the soil
      REAL PHRATE			! IN:pH rate modifier (proportion)
	REAL TRATE			! IN:temperature rate modifying factor (prop.)
	REAL WRATE			! IN:Moisture rate modifyer (prop.)
C
C m(water) = linear increase from 0.2 at no water to 1 at 0.5 of maximum
C
	if(LAYERSOILW.GT.0)then
        if((LAYERSOILW/LAYERWSAT).GT.0.5)then
          WRATE=1.
  	  else
          if((LAYERSOILW/LAYERWSAT).LT.0)then
            WRATE=0.2
	    else
  	      WRATE=0.2+(1.6*(LAYERSOILW/LAYERWSAT))
	    endif
        endif
	else
        WRATE=0.2
	endif
C
C m(temp) = T/25 (T = temp of 5cm layer), zero below 0C and 1 above 25C
C
	if(THISTEMP.LT.0)then
	  TRATE=0.
	endif
	if(THISTEMP.LT.25)then
	  TRATE=THISTEMP/25.
	else
	  TRATE=1.
	endif
C
C m(pH) = 0.5 below 2.0, 1 above 7.0 and increases linearly from 2.0 to 7.0
C
	if(PH.LT.2)then
        PHRATE=0.5
	else
	  if(PH.GT.7)then
          PHRATE=1.
	  else
          PHRATE=(PH+3.0)/10.0
	  endif
	endif
C
C m(crop) assumed to be 1
C
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE MODFACTS_MINER(LAYERSOILW,LAYERWMAX,LAYERWSAT,
     &                    WRATED,WRATER,
     &					THISTEMP,TRATE,
     &					PH,PHP1,PHP2,PHRATE,
     &					ICOVER,CRRATE,ITFUNC,IMFUNC)
C
C Subroutine to calculate the rate modifying factors
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C ...Profile descriptors
C
      INTEGER MAXLAYER1	! Maximum number of layers in the profile
	INTEGER MAXDEPTH	! Maximum depth of the profile
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
C
C ...ROTHC variables
C
      REAL PHRATE_MIN		! Minimum pH rate modifier
      REAL WMOD			! Water deficit (mm/layer)
	REAL ONEBAR			! Water deficit at -1bar (20mm / 25cm)
	REAL WMIN			! Minimum water level as fraction of the maximum
	REAL WMINR			! Minimum water level as fraction of the maximum
	REAL WMIND			! Minimum water level as fraction of the maximum
C	DATA WMIN,WMINR,WMIND /0.2,0.2,0.6/		! Bente's values
	DATA WMIN,WMINR,WMIND /0.2,0.2,0.2/		! Values used in RothC
C
C ...Hadley centre variables
C
	REAL STH_SOIL		! Soil moisture as a fraction of saturation 
						! (m3/m3) (value is between (~0.3 and 1.0)
	REAL STH_WILT		! Wilting soil moisture as a fraction of saturation.
	REAL STH_OPT		! Fractional soil moisture at which respiration is 
						! maximum
	REAL Q10			! Q10 value = 2.0
	DATA Q10 /2.0/
	REAL TSOIL			! Soil temperature (in degrees kelvin) 
	REAL WILTP			! Wilting point in mm
C
C Variables passed to /from this subroutine
C ... Weather variables
C
      REAL THISTEMP		! IN:Soil temperature in this layer(deg.C)
C
C ...Soil descriptors
C
	REAL LAYERSOILW		! IN:Available soil water (mm/layer)
	REAL LAYERWMAX		! IN:Maximum available soil water (mm/layer)
	REAL LAYERWSAT		! IN:Available water at saturation (mm/layer)
	REAL PH				! IN:Current pH of the soil
	REAL PHP1			! IN:pH below which rate is zero
	REAL PHP2			! IN:pH above which rate is optimum
      INTEGER ICOVER		! IN:Code for crop cover: 1=covered; 2=bare

      REAL CRRATE			! IN:crop modifying factor (proportion)
      REAL PHRATE			! IN:pH rate modifier (proportion)
	REAL TRATE			! IN:temperature rate modifying factor (prop.)
	REAL WRATED			! IN:Moisture rate modifyer for DPM and BIO (prop.)
	REAL WRATER			! IN:Moisture rate modifyer for RPM and HUM (prop.)
C
C ...Model descriptors
C
	INTEGER ITFUNC				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER ITFUNC_ROTHC		! ROTHC temperature rate modifier
	INTEGER ITFUNC_HADLEY		! HADLEY temperature rate modifier
      DATA ITFUNC_ROTHC,ITFUNC_HADLEY /0,1/
	INTEGER IMFUNC				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IMFUNC_ROTHC		! ROTHC moisture rate modifier
	INTEGER IMFUNC_HADLEY		! HADLEY moisture rate modifier
      DATA IMFUNC_ROTHC,IMFUNC_HADLEY /0,1/
C---------------------------------
C CROP MODIFICATION FACTOR, CRRATE
C---------------------------------
      IF(ICOVER.EQ.1)THEN
	  CRRATE = 0.6
	ELSE
	  CRRATE= 1.0
	ENDIF
C----------------------------------------------
C MOISTURE MODIFICATION FACTOR, WRATER & WRATED
C----------------------------------------------
C ROTH-C moisture rate modification factors
C
      IF(IMFUNC.EQ.IMFUNC_ROTHC)THEN
C
C Calcualate the water modifying factor, WRATE,
C Changed include reduction of water modifier above field capacity
C Based on Devevre & Horwarth (2000) and the assumption that resistant
C and decomposable material react differently. 
C
        WMOD=LAYERWMAX-LAYERSOILW
        ONEBAR=20.0*MAXDEPTH/(MAXLAYER1*25.)
        IF(WMOD.LE.0)THEN					
          WRATER=((WMINR-1)*LAYERSOILW)+LAYERWSAT-(WMINR*LAYERWMAX)
	    WRATER=WRATER/(LAYERWSAT-LAYERWMAX)
	    WRATED=((WMIND-1)*LAYERSOILW)+LAYERWSAT-(WMIND*LAYERWMAX)
	    WRATED=WRATED/(LAYERWSAT-LAYERWMAX)
        ELSEIF(WMOD.LE.ONEBAR.AND.WMOD.GT.0)THEN
          WRATER=1.0
          WRATED=1.0
        ELSEIF(WMOD.GT.ONEBAR)THEN
          WRATER=1.0-((1-WMIN)*(WMOD-ONEBAR))/(LAYERWMAX-ONEBAR)
          WRATED=1.0-((1-WMIN)*(WMOD-ONEBAR))/(LAYERWMAX-ONEBAR)
        ENDIF
   
        IF(WRATER.LT.0.2)then
        WRATER=0.2
        WRATED=0.2
        endif
C
C Hadley Centre moisture rate modification factors
C
	ELSEIF(IMFUNC.EQ.IMFUNC_HADLEY)THEN
C Note: STH_WILT set to 0.3 for all soils. 
C Check: Does this introduces an error? 
C        If it does, must pass wilting point from wrapper
        STH_WILT=0.3
	  WILTP=STH_WILT*LAYERWSAT/(1+STH_WILT)
	  STH_SOIL=(LAYERSOILW+WILTP)/(LAYERWSAT+WILTP)
        STH_OPT=0.5*(1+STH_WILT)                                  
        IF(STH_SOIL.LE.STH_WILT)THEN                             
          WRATER=0.2                                                     
          WRATED=0.2                                                     
        ELSEIF(STH_SOIL.GT.STH_WILT.AND.STH_SOIL.LE.STH_OPT)THEN
          WRATER=0.2+0.8*((STH_SOIL-STH_WILT)/(STH_OPT-STH_WILT))                      
          WRATED=0.2+0.8*((STH_SOIL-STH_WILT)/(STH_OPT-STH_WILT))                      
        ELSEIF(STH_SOIL.GT.STH_OPT)THEN                          
          WRATER=1-0.8*(STH_SOIL-STH_OPT)                      
          WRATED=1-0.8*(STH_SOIL-STH_OPT)                      
        ENDIF
	ENDIF
C---------------------------------------
C TEMPERATURE MODIFICATION FACTOR, TRATE
C---------------------------------------
C ROTH-C temperature rate modification factors
C
      IF(ITFUNC.EQ.ITFUNC_ROTHC)THEN
C
C Calculate the temperature modifying factor, TRATE
C
        IF(THISTEMP.LT.-10)THEN
	    TRATE=0
	  ELSE
	    TRATE=47.91/(1+EXP(106.06/(THISTEMP+18.27)))
	  ENDIF
C
C Hadley Centre temperature rate modification factors
C
	ELSEIF(ITFUNC.EQ.ITFUNC_HADLEY)THEN
        TSOIL=THISTEMP+273.15
        TRATE=Q10**(0.1*(TSOIL-298.15))                      
	ENDIF
C-------------------------------
C PH MODIFICATION FACTOR, PHRATE
C-------------------------------
      PHRATE_MIN=0.2
      IF(PH.LE.PHP1)THEN
	  PHRATE=PHRATE_MIN
	ELSEIF(PH.GT.PHP1.AND.PH.LT.PHP2)THEN
	  PHRATE=PHRATE_MIN+(1-PHRATE_MIN)*(PH-PHP1)/(PHP2-PHP1)
	ELSEIF(PH.GE.PHP2)THEN
	  PHRATE=1
	ENDIF
C
C Leave MODFACTS
C
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE MODFACTS_NITRIF(LAYERSOILW,LAYERWMAX,LAYERWSAT,WRATE,
     &					THISTEMP,TRATE,
     &					PH,PHP1,PHP2,PHRATE,
     &					ICOVER,CRRATE,ITFUNC,IMFUNC)
C
C Subroutine to calculate the rate modifying factors
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C ...Layer descriptors
C
      INTEGER MAXLAYER1	! Maximum number of layers in the profile
	INTEGER MAXDEPTH	! Maximum depth of the profile
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
C
C ...ROTH-C variables
C
      REAL WMOD			! Water deficit (mm/layer)
	REAL ONEBAR			! Water deficit at -1bar (20mm / 25cm)
	REAL WMIN			! Minimum water level as fraction of the maximum
	DATA WMIN /0.6/
	REAL NWATER			! Fitted constant describing water modifier for nitrification
	DATA NWATER /0.8/
	REAL NTEMP1			! Fitted constant describing temperature modifier for nitrification
	DATA NTEMP1 /47.9/
	REAL NTEMP2			! Fitted constant describing temperature modifier for nitrification
	DATA NTEMP2 /106.06/
	REAL NTEMP3			! Fitted constant describing temperature modifier for nitrification
	DATA NTEMP3 /18.27/
C
C ...Hadley centre variables
C
c	REAL V_SAT			! Volumetric soil moisture content at saturation 
						! (m3 H2O/m3 soil)(typical value 0.458150)
	REAL STH_SOIL		! Soil moisture as a fraction of saturation 
						! (m3/m3) (value is between (~0.3 and 1.0)
	REAL STH_WILT		! Wilting soil moisture as a fraction of saturation.
	REAL STH_OPT		! Fractional soil moisture at which respiration is 
						! maximum
	REAL WILTP			! Wilting point in mm
	REAL Q10			! Q10 value = 2.0
	DATA Q10 /2.0/
	REAL TSOIL			! Soil temperature (in degrees kelvin) 
C
C Variables passed to /from this subroutine
C ...Weather descriptors
C
      REAL THISTEMP		! IN:Soil temperature in this layer(deg.C)
C
C ...Soil descriptors
C
	REAL LAYERSOILW		! IN:Available soil water (mm/layer)
	REAL LAYERWMAX		! IN:Maximum available soil water (mm/layer)
	REAL LAYERWSAT		! IN:Available water at saturation (mm/layer)
	REAL PH				! IN:Current pH of the soil
	REAL PHP1			! IN:pH below which rate is zero
	REAL PHP2			! IN:pH above which rate is optimum
      INTEGER ICOVER		! IN:Code for crop cover: 1=covered; 2=bare

      REAL CRRATE			! OUT crop modifying factor (proportion)
      REAL PHRATE			! OUT pH rate modifier (proportion)
	REAL TRATE			! OUT temperature rate modifying factor (prop.)
	REAL WRATE			! OUT Moisture rate modifying factor (proportion)
C
C ...Model descriptors
C
	INTEGER ITFUNC				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER ITFUNC_ROTHC		! ROTHC temperature rate modifier
	INTEGER ITFUNC_HADLEY		! HADLEY temperature rate modifier
      DATA ITFUNC_ROTHC,ITFUNC_HADLEY /0,1/
	INTEGER IMFUNC				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IMFUNC_ROTHC		! ROTHC moisture rate modifier
	INTEGER IMFUNC_HADLEY		! HADLEY moisture rate modifier
      DATA IMFUNC_ROTHC,IMFUNC_HADLEY /0,1/
C---------------------------------
C CROP MODIFICATION FACTOR, CRRATE
C---------------------------------
      IF(ICOVER.EQ.1)THEN
	  CRRATE = 0.6
	ELSE
	  CRRATE= 1.0
	ENDIF
C----------------------------------------------
C MOISTURE MODIFICATION FACTOR, WRATER & WRATED
C----------------------------------------------
C ROTH-C moisture rate modification factors
C Calcualate the water modifying factor, WRATE,
C The modification of anaerobic nitrification is set
C so that at saturation, no nitrification takes place
C as it is an aerobic process
C
      IF(IMFUNC.EQ.IMFUNC_ROTHC)THEN
        WMOD=LAYERWMAX-LAYERSOILW
        ONEBAR=20.0*MAXDEPTH/(MAXLAYER1*25.)
  	  IF(LAYERSOILW.GT.LAYERWMAX)THEN
	    WRATE=(LAYERWSAT-LAYERSOILW)/(LAYERWSAT-LAYERWMAX)
        ELSEIF(WMOD.LE.ONEBAR)THEN
          WRATE=1.0
        ELSE
          WRATE=1.0-((1-WMIN)*(WMOD-ONEBAR))/(LAYERWMAX-ONEBAR)
        ENDIF
C
C Hadley Centre moisture rate modification factors
C
	ELSEIF(IMFUNC.EQ.IMFUNC_HADLEY)THEN
C Note: STH_WILT set to 0.3 for all soils. 
C Check: Does this introduces an error? 
C        If it does, must pass wilting point from wrapper
        STH_WILT=0.3
	  WILTP=STH_WILT*LAYERWSAT/(1+STH_WILT)
	  STH_SOIL=(LAYERSOILW+WILTP)/(LAYERWSAT+WILTP)
        STH_OPT=0.5*(1+STH_WILT)                                  
        IF(STH_SOIL.LE.STH_WILT)THEN                             
          WRATE=0.2                                                     
        ELSEIF(STH_SOIL.GT.STH_WILT.AND.STH_SOIL.LE.STH_OPT)THEN
          WRATE=0.2+NWATER*((STH_SOIL-STH_WILT)/(STH_OPT-STH_WILT))                      
        ELSEIF(STH_SOIL.GT.STH_OPT)THEN                          
          WRATE=1-NWATER*(STH_SOIL-STH_OPT)                      
        ENDIF
	ENDIF
C---------------------------------------
C TEMPERATURE MODIFICATION FACTOR, TRATE
C---------------------------------------
C ROTH-C temperature rate modification factors
C
      IF(ITFUNC.EQ.ITFUNC_ROTHC)THEN
        IF(THISTEMP.LT.-10)THEN
  	    TRATE=0
	  ELSE
	    TRATE=NTEMP1/(1+EXP(NTEMP2/(THISTEMP+NTEMP3)))
	  ENDIF
C
C Hadley Centre moisture rate modification factors
C
	ELSEIF(ITFUNC.EQ.ITFUNC_HADLEY)THEN
        TSOIL=THISTEMP+273.15
        TRATE=Q10**(0.1*(TSOIL-298.15))                      
	ENDIF
C-------------------------------
C PH MODIFICATION FACTOR, PHRATE
C-------------------------------
	PHRATE=0.333*PH-1.333						! Ste-Marie & Pare, 1999
	PHRATE=0.56+(ATAN(3.14*0.45*(-5.+PH)))/3.14	! Parton et al, 1996
	IF(PHRATE.LT.0)THEN 
	  PHRATE=0.
	ElSEIF(PHRATE.GT.1.5)THEN
	  PHRATE=1.5
	ENDIF
C
C Leave MODFACTS
C
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE NITRIF_NEMIS(SOILTEMP,RAIN,NSOIL,ITFUNC,IMFUNC,
     &                  SOILW,WMAX,WSAT,
     &                  SOILN,AMMN,SOIL15,AMMN15,CRIT,TAM,
     &                  FERTADD,FERTADD15,THISFERT,THISFERT15,
     &                  FANIT,FANIT15,VOLAT,VOLAT15,CONVER_F,
     &                  PHARRAY,PHP1ARRAY,PHP2ARRAY,
     &                  GNN2O,GNNO,GPNN2O,G15NN2O,G15NNO,G15PNN2O,
     &                  NITRIFN,VOLATN)	

C
      IMPLICIT NONE




C
C Variables local to this subroutine
C
	INTEGER MAXFERTSTEPS		! Maximum number of timesteps 
								!	over which fertiliser is added
      PARAMETER (MAXFERTSTEPS=30*60*24*7*3)	! Timesteps = 3 wks (max.30 secs)
	INTEGER MAXSOIL				! Max.no.of soil types
      PARAMETER (MAXSOIL=50)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	INTEGER MAXDEPTH			! Maximum depth of the soil profile
	INTEGER IDEPTH				! Depth of this layer
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
	REAL FCGAS					! Prop.N2O loss by partial nitrificn at FC
	REAL PE1					! Prop.gas.loss by nitrificn
	REAL PE2					! Prop.N2O loss by part.nitrificn at this 
								! water content
	REAL PE3					! Prop.gas.loss by nitrificn that is NO
	DATA PE1 /0.02/
	DATA FCGAS /0.02/
      DATA PE3 /0.4/
	REAL CRITL					! Critical minimum N in the layer (kgN/ha)
	REAL AMM1					! Ammonium N remaining after nitrificn(kgN/ha)
	REAL AMM2					! Fert.N remaining after nitrificn(kgN/ha)
	REAL ANIT					! Ammonium N nitrified (kgN/ha)
	REAL ANIT15					! Ammonium N15 nitrified (kgN/ha)
	REAL FANIT					! Fertiliser N nitrified (kgN/ha)
	REAL FANIT15				! Fertiliser N15 nitrified (kgN/ha)
	REAL FLOSS					! Fraction lost by volatn.if rain < 5mm / week
	REAL ULOSS					! Urea lost by vol.if rain < 5mm/week (kg/ha)
	REAL UFERT					! Urea in fertiliser (kgN/ha)
	REAL ALOSS					! Ammn.lost by vol.if rain < 5mm/week (kg/ha)
	REAL AMFERT					! Ammn.in fertiliser (kgN/ha)
	REAL AFERT					! Ammn.and urea in fertiliser (kgN/ha)
	REAL ULOSS15				! Urea lost by vol.if rain < 5mm/week (kg/ha)
	REAL UFERT15				! Urea in fertiliser (kgN/ha)
	REAL ALOSS15				! Ammn.lost by vol.if rain < 5mm/week (kg/ha)
	REAL AMFERT15				! Ammn.in fertiliser (kgN/ha)
	REAL AFERT15				! Ammn.and urea in fertiliser (kgN/ha)
	INTEGER IL					! Layer counter
	REAL TEMPN1					! Sum of N added to nitrate pool
	REAL TEMPN15				! Sum of N15 added to nitrate pool
	REAL KNITRIF				! Rate constant for nitrification (/week)
	DATA KNITRIF /0.6/			! Rate constant for nitrification

C Variables passed to/from calling subroutine
C ....Timing factors
C
      REAL CONVER_F				! IN/OUT:Conversion to correct timestep
C
C ...Model descriptors
C
	INTEGER ITFUNC				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IMFUNC				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
C
C ...Weather factors
C
	REAL SOILTEMP(MAXLAYER)		! IN: Soil temperature (deg.C/timestep)
	REAL RAIN					! IN: Rainfall (mm/timestep)
C
C ...Soil factors
C
	INTEGER NSOIL				! IN:Soil code number
	REAL VOLAT					! IN:N lost by volatilisation (kgN/ha)
	REAL VOLAT15				! IN:N15 lost by volatilisation (kgN15/ha)
	REAL VOLATN(MAXLAYER)		! OUT:Volatilisation from this layer (kg/ha/layer/timestep)
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
	REAL CRIT(MAXSOIL,MAXLAYER)	! IN/OUT: Minimum nitrate level (Sozanska)
	REAL TAM(MAXLAYER)			! IN:15N/N in the ammonium in the layer
	REAL GNN2O					! IN:N2O lost by nitrification (kgN/ha)
	REAL GNNO					! IN:NO lost by nitrification (kgN/ha)
	REAL GPNN2O					! IN:N2O lost by part.nitrification (kgN/ha)
	REAL G15NN2O				! IN:15N2O lost by nitrification (kgN15/ha)
	REAL G15NNO					! IN:15NO lost by nitrification (kgN15/ha)			
	REAL G15PNN2O				! IN:15N2O lost by part.nitrif.(kgN15/ha)			
	REAL NITRIFN(MAXLAYER)		! OUT:Nitrification from this layer (kg/ha/layer/timestep)
	REAL WRATE					! IN:Moisture rate modifyer (prop.)
      REAL CRRATE					! IN:crop modifying factor (proportion)
      REAL PHRATE					! IN:pH rate modifier (proportion)
	REAL TRATE					! IN:temperature rate modifying factor (prop.)
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)! IN:pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)! IN:pH above which decomp.max.
C
C ... Fertiliser factors
C
	REAL THISFERT				! IN:N input by fertiliser (kgN/ha)
	REAL THISFERT15				! IN:N15 input by fertiliser (kgN15/ha)
	REAL FERTADD(MAXFERTSTEPS,4)	! INOUT:Fertiliser N added in this time 
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
	REAL FERTADD15(MAXFERTSTEPS,4)	! INOUT:Fertiliser N15 added in this time
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
C
C Initialize Variables
C
	GNN2O=0															
	GNNO=0
	GPNN2O=0
	G15NN2O=0
	G15NNO=0
	G15PNN2O=0
      FANIT=0
      FANIT15=0
C
C For each layer in the profile to 25cm calculate nitrification
C
      DO 100 IL=1,5
        CALL MODFACTS_NITRIF(SOILW(IL),WMAX(IL),WSAT(IL),WRATE,
     &                SOILTEMP(IL),TRATE,								
     &			    PHARRAY(NSOIL,IL),PHP1ARRAY(NSOIL,IL),
     &                PHP2ARRAY(NSOIL,IL),PHRATE,
     &                0,CRRATE,ITFUNC,IMFUNC)										  														
C
C Calculate ANIT = Week's nitrification
C
C       CRITL=CRIT(NSOIL,IL)*MAXDEPTH/(50.*MAXLAYER1)
C       AMM1=AMMN(IL)*EXP(-KNITRIF*TRATE*WRATE)
C       ANIT=(AMMN(IL)-AMM1)*PHRATE*CONVER_F
C       IF(ANIT.GT.(AMMN(IL)-CRITL))THEN
C       ANIT=AMMN(IL)-CRITL
C       END IF
C	  NITRIFN(IL)=ANIT
C	  VOLATN(IL)=0
C
        CRITL=CRIT(NSOIL,IL)*MAXDEPTH/(50.*MAXLAYER1)
c       Changed the conversion factor to change from weekly rate to daily rate ##possibly fudge
        ANIT=AMMN(IL)*(1-EXP(-0.6*TRATE*WRATE*PHRATE*CONVER_F))
c       added a saturation function to simulate limiting nitrifiers##did the same for fertiliser addition below
     &  *(1-AMMN(IL)/(50+AFERT+AMMN(IL)))
c        ANIT=(AMMN(IL)-AMM1)
c        IF(ANIT.GT.(AMMN(IL)-CRITL))THEN
c          ANIT=AMMN(IL)-CRITL
c        END IF
        NITRIFN(IL)=ANIT
	  VOLATN(IL)=0
C
C Take nitrified N from the NH4-N pool
C
        ANIT15=ANIT*TAM(IL)
        AMMN(IL)=AMMN(IL)-ANIT
        AMMN15(IL)=AMMN15(IL)-ANIT15
C
C ....Calculate the gaseous losses associated with nitrification
C
	  PE2=FCGAS*SOILW(IL)/WMAX(IL)
	  GPNN2O=GPNN2O+(PE2*ANIT)
	  GNNO=GNNO+(PE3*PE1*ANIT)
	  GNN2O=GNN2O+(PE1*(1-PE3)*ANIT)
        G15PNN2O=G15PNN2O+(PE2*ANIT15)
        G15NNO=G15NNO+(PE3*PE1*ANIT15)
        G15NN2O=G15NN2O+(PE1*(1-PE3)*ANIT15)
C
C Add nitrified N to the nitrate pool
C TEMPN1=sum of N added to nitrate pool
C TEMPN15=sum of N15 added to nitrate pool
C
        TEMPN1=0
        TEMPN15=0
        SOILN(IL)=SOILN(IL)+ANIT
        SOIL15(IL)=SOIL15(IL)+ANIT15
        TEMPN1=TEMPN1+SOILN(IL)
        TEMPN15=TEMPN15+SOIL15(IL)
 100  CONTINUE
C
C For this weeks fertilizer addition...
C Set the proportion of the added fertilizer that is NH4 and the fraction
C that is lost by volatalization (FLOSS) if rain < 5mm in week.
C
         IF(RAIN.LT.(5*CONVER_F))THEN
          FLOSS=0.15
         ELSE
          FLOSS=0.0
         END IF
C
         ULOSS=FLOSS*FERTADD(1,4)
	   VOLATN(1)=VOLATN(1)+ULOSS
         UFERT=FERTADD(1,4)
         ULOSS15=FLOSS*FERTADD15(1,4)
         UFERT15=FERTADD15(1,4)
         THISFERT=THISFERT+FERTADD(1,4)
         THISFERT15=THISFERT15+FERTADD15(1,4)
C
         UFERT=UFERT-ULOSS
         UFERT15=UFERT15-ULOSS15
C
C If NH4 fertilizer is ammonium sulphate, (IVOL=0), also lose by
C volatilization
C
         ALOSS=FLOSS*FERTADD(1,2)
         ALOSS15=FLOSS*FERTADD15(1,2)
	   VOLATN(1)=VOLATN(1)+ALOSS
C
C Calc. amount of N as NH4
C
         AMFERT=FERTADD(1,2)+FERTADD(1,3)
         AMFERT15=FERTADD15(1,2)+FERTADD15(1,3)
         THISFERT=THISFERT+FERTADD(1,2)+FERTADD(1,3)
         THISFERT15=THISFERT15+FERTADD15(1,2)+FERTADD15(1,3)
         AMFERT=AMFERT-ALOSS
         AMFERT15=AMFERT15-ALOSS
C
C Calculate the amount of added fertilizer lost by volatalization
C
         AFERT=AMFERT+UFERT
         AFERT15=AMFERT15+UFERT15
         VOLAT=VOLAT+ULOSS+ALOSS
         VOLAT15=VOLAT15+ULOSS15+ALOSS15
C
C Calculate the amount of NH4 fertilizer nitrified
C
        CALL MODFACTS_NITRIF(SOILW(1),WMAX(1),WSAT(1),WRATE,
     &                SOILTEMP(1),TRATE,
     &			    PHARRAY(NSOIL,1),PHP1ARRAY(NSOIL,1),
     &                PHP2ARRAY(NSOIL,1),PHRATE,
     &                0,CRRATE,ITFUNC,IMFUNC)
	   FANIT=AFERT*(1-EXP(-0.6*WRATE*TRATE*PHRATE*CONVER_F))
     &   *(1-AFERT/(50+AFERT+AMMN(1)))
c         AMM2=AFERT*EXP(-0.6*WRATE*TRATE)
c         FANIT=(AFERT-AMM2)*PHRATE*CONVER_F
c	   FANIT=AFERT-AMM2
	   NITRIFN(1)=NITRIFN(1)+FANIT
C
C ....Calculate the gaseous losses associated with FERTILISER nitrification
C
	   GPNN2O=GPNN2O+(PE2*FANIT)
	   GNNO=GNNO+(PE3*PE1*FANIT)
	   GNN2O=GNN2O+(PE1*(1-PE3)*FANIT)
         AFERT=AFERT-FANIT
         AMMN(1)=AMMN(1)+AFERT
         SOILN(1)=SOILN(1)+FANIT
         TEMPN1=TEMPN1+FANIT
         IF(FERTADD15(1,2)+FERTADD15(1,3)+FERTADD15(1,4).GT.0)THEN
          FANIT15=FANIT
          G15PNN2O=G15PNN2O+(PE2*FANIT15)
          G15NNO=G15NNO+(PE3*PE1*FANIT15)
	    G15NN2O=G15NN2O+(PE1*(1-PE3)*FANIT15)
          AMMN15(1)=AMMN15(1)+AFERT
          SOIL15(1)=SOIL15(1)+FANIT15
          TEMPN15=TEMPN15+FANIT15
         END IF
C
C Go back to next fertilizer addition
C
  400 CONTINUE
C
C Leave NITRIF
C
      RETURN
      END
C
C----------------------------------------------------------
C
      SUBROUTINE OMPROP(NSOIL,HY,HZ1,ALPHA,BETA,GAMMA,DELTA,
     &                  BPART,HPART,BPROP,HPROP,PHARRAY)
C
C Subroutine to set BIO and HUM decompostion parameters
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      INTEGER MAXSOIL,MAXLAYER
      PARAMETER (MAXSOIL=50,MAXLAYER=60)
	INTEGER IL
	REAL PH(MAXLAYER),STABCN
C
C Variables passed to/from calling subroutine
C
      INTEGER NSOIL
      REAL HY(MAXSOIL,MAXLAYER),HZ1(MAXLAYER)
	REAL ALPHA(MAXLAYER),BETA(MAXLAYER)
	REAL GAMMA(MAXLAYER),DELTA(MAXLAYER)
      REAL BPART(MAXSOIL,MAXLAYER),HPART(MAXSOIL,MAXLAYER)
	REAL BPROP(MAXSOIL,MAXLAYER),HPROP(MAXSOIL,MAXLAYER)
	REAL PHARRAY(MAXSOIL,MAXLAYER)
C
C For each layer down the profile set the parameters
C
      DO 100 IL=1,MAXLAYER
C
C Set N:C Ratio of N in HUM and ROOT pools
C
	  PH(IL)=PHARRAY(NSOIL,IL)
	  CALL GETSTABLECN(PH(IL),STABCN)
	  HY(NSOIL,IL)=1/STABCN
C
C BPART and HPART = proportion of Carbon going to BIO + HUM (=1-EFFICIENCY)
C
        HPART(NSOIL,IL)=BPART(NSOIL,IL)
C
C ALPHA = Proportion of decomposing biomass that is transformed to biomass
C BETA = Proportion of decomposing biomass that is transformed to humus
C 
        BETA(IL)=BPART(NSOIL,IL)/(1+BPROP(NSOIL,IL))
        ALPHA(IL)=BPART(NSOIL,IL)-BETA(IL)
C
C GAMMA = Proportion of decomposing humus that is transformed to biomass
C DELTA = Proportion of decomposing humus that is transformed to humus
C
        DELTA(IL)=HPART(NSOIL,IL)/(1+HPROP(NSOIL,IL))
        GAMMA(IL)=HPART(NSOIL,IL)-DELTA(IL)

C
C Calc. N:C of input plant material for steady state conditions = HZ1
C
        HZ1(IL)=HY(NSOIL,IL)*(ALPHA(IL)+BETA(IL))
100   CONTINUE
C
C Leave OMPROP
C
      RETURN
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE PAR_CN(NUMSOIL,HY,BPART,BPROP,HPROP,BIOP,CRIT,
     &           BIORATE,HUMRATE,BULKDENS,
     &           TOCARRAY,IOMARRAY,LUARRAY,CLARRAY,SNAME,DFACT,
     &           FYMSTART,FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,FYMLOSS,FYMPOS,
     &           PHARRAY,PHP1ARRAY,PHP2ARRAY)
C
C Subroutine to reset parameters from files PARAM.OUT and PARAM.DAT
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	INTEGER MAXDEPTH			! Maximum depth of the soil profile
	INTEGER IDEPTH				! Depth of this layer
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
	INTEGER MAXSOIL				! Max.no.of soil types
	INTEGER MAXORGM				! Max.no.of org.manure applications allowed
      PARAMETER (MAXSOIL=50,MAXORGM=52)
      CHARACTER*40 TEMP
	INTEGER I					! Local counter 
	INTEGER J					! Local counter 
	INTEGER IL					! Local layer counter
	REAL FNULL					! Unwanted real number
C
C Variables passed to/from calling subroutine
C ... Soil factors
C
      INTEGER NUMSOIL				! Number of soils defined
	REAL HY(MAXSOIL,MAXLAYER)	! IN/OUT:Stable N:C ratio of BIO & HUM pools
	REAL CRIT(MAXSOIL,MAXLAYER)	! IN/OUT: Minimum nitrate level (Sozanska)
	REAL BPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL HPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL BPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from biomass decompositn
	REAL HPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from humus decompositn
	REAL BIOP(MAXSOIL,MAXLAYER)		! IN/OUT: BIO/TOC 
	REAL BIORATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for BIO decompn/yr
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for HUM decompn/yr
	REAL IOMARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium IOM (kgC/ha/layer)
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium TOC (kgC/ha/layer)
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: % clay content in this layer
      INTEGER LUARRAY(MAXSOIL)	! OUT:Land use before equilibrium 
	REAL DEPTH
	CHARACTER*40 SNAME(MAXSOIL)	! OUT:Soil name
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)! OUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)! OUT: pH above which decomp.max.
	REAL DFACT(MAXLAYER)		! OUT: Denitrification factor
	REAL BULKDENS(MAXLAYER)     ! OUT: Bulk density of each layer [g/cm3]
C
C ... Manure factors
C
      REAL FPROPX(MAXORGM)		! OUT:Proportion of FYM added to top 25cm
	REAL FYMCX(MAXORGM)			! OUT:Proportion of C in FYM
	REAL FYMNX(MAXORGM)			! OUT:Proportion of Organic N in FYM
	REAL FYMAX(MAXORGM)			! OUT:Proportion of ammonium-N in FYM
	REAL FYMWX(MAXORGM)			! OUT:Amount of water added in FYM
	REAL FYMLOSS(MAXORGM)		! OUT:Prop.FYM lost by volatn/timestep
	REAL FYMPOS(MAXORGM)		! OUT:Proportion of FYM added to top 25cm
	REAL FYMSTART(MAXORGM)		! OUT:Amount of rainfall in 1 week, 
								! 	 below which volatilisation will occur
C
C Soil Parameters
C
C Soil name
      SNAME(1)='Sand'
	SNAME(2)='Loam'
	SNAME(3)='Clay'
	SNAME(4)='User'
C
C Get parameters for each layer
C
      DO 100 IL=1,MAXLAYER1
C Rate of biomass decomposition
	  BIORATE(1,IL)=0.0127*50
	  BIORATE(2,IL)=0.0127*50
	  BIORATE(3,IL)=0.0127*50
C Rate of humus decomposition
        HUMRATE(1,IL)=0.0004*50
        HUMRATE(2,IL)=0.0004*50
        HUMRATE(3,IL)=0.0004*50
C N:C ratio of BIO and HUM
        HY(1,IL)=0.118
	  HY(2,IL)=0.118
	  HY(3,IL)=0.118
	  HY(4,IL)=0.083
C Fraction of BIO + HUM formed (alpha+beta)
        BPART(1,IL)=0.35
	  BPART(2,IL)=0.4
	  BPART(3,IL)=0.42
	  BPART(4,IL)=0.27
C Propn. BIO:HUM (alpha/beta) from BIO
        BPROP(1,IL)=1.1
        BPROP(2,IL)=1.1
        BPROP(3,IL)=1.1
        BPROP(4,IL)=1.1
C Propn. BIO:HUM (alpha/beta) from HUM
        HPROP(1,IL)=1.1
        HPROP(2,IL)=1.1
        HPROP(3,IL)=1.1
        HPROP(4,IL)=1.1
C 8 Fraction of BIO in total organic C and N 
        BIOP(1,IL)=0.028
        BIOP(2,IL)=0.028
        BIOP(3,IL)=0.028
        BIOP(4,IL)=0.028
C  Minimum Level of Nitrate in soil
        CRIT(1,IL)=5.0
	  CRIT(2,IL)=10.0
	  CRIT(3,IL)=15.0
	  CRIT(4,IL)=10.0
C Denitrification factor
        DFACT(IL)=0.005
100   CONTINUE
C
C Organic Manure Parameters
C Manure types
C (1) Cattle FYM
C (2) Pig FYM
C (3) Layer Manure
C (4) Broiler/Turkey Manure
C (5) Sewage sludge cake (undigested)
C (6) Sewage sludge cake (digested)
C (7) Dairy slurry
C (8) Beef slurry
C (9) Pig Slurry
C (10) Strainer box seperated slurry
C (11) Weeping wall seperated slurry
C (12) Mechanically seperated slurry
C (13) Sewage sludge liquid (undigested)
C (14) Sewage sludge liquid (digested)
C
C FPROPX - Amount of N in fresh manure that is decomposed organic matter (kg N / t manure)
C
      FPROPX(1)=100
      FPROPX(2)=100
      FPROPX(3)=250
      FPROPX(4)=350
      FPROPX(5)=350
      FPROPX(6)=350
      FPROPX(7)=200
      FPROPX(8)=200
      FPROPX(9)=320
      FPROPX(10)=200
      FPROPX(11)=250
      FPROPX(12)=200
      FPROPX(13)=150
      FPROPX(14)=150

C     FPROPX(1)=0.2
C     FPROPX(2)=0.2
C     FPROPX(3)=0.2
C     FPROPX(4)=0.2	
C     FPROPX(5)=0.2
C     FPROPX(6)=0.2
C     FPROPX(7)=0.2

C
C FYMCX - C in fresh manure (kg C / t manure)
C
	FYMCX(1)=52.5
	FYMCX(2)=62.5
	FYMCX(3)=96.0
	FYMCX(4)=180.0
	FYMCX(5)=62.5
	FYMCX(6)=62.5
	FYMCX(7)=12.0
	FYMCX(8)=12.0
	FYMCX(9)=15.0
	FYMCX(10)=4.5
	FYMCX(11)=9.0
	FYMCX(12)=12.0
	FYMCX(13)=12.5
	FYMCX(14)=12.0

C	FYMCX(1)=80.0
C	FYMCX(2)=80.0
C	FYMCX(3)=80.0
C	FYMCX(4)=12.0
C	FYMCX(5)=80.0
C	FYMCX(6)=12.0
C	FYMCX(7)=24.0

C
C FYMNX - N in fresh manure (kg N / t manure)
C
      FYMNX(1)=6.0
      FYMNX(2)=7.0
      FYMNX(3)=15.0
      FYMNX(4)=28.8
      FYMNX(5)=7.5
      FYMNX(6)=7.5
      FYMNX(7)=3.0
      FYMNX(8)=2.28
      FYMNX(9)=4.98
      FYMNX(10)=1.5
      FYMNX(11)=3.0
      FYMNX(12)=4.0
      FYMNX(13)=1.8
      FYMNX(14)=2.0

C     FYMNX(1)=5.5
C     FYMNX(2)=5.5
C     FYMNX(3)=5.5
C     FYMNX(4)=5.5
C     FYMNX(5)=5.5
C     FYMNX(6)=5.5
C     FYMNX(7)=8.2

C
C FYMAX - NH4+ in fresh manure (kg N / t manure)
C
	FYMAX(1)=1.2
	FYMAX(2)=1.4
	FYMAX(3)=7.5
	FYMAX(4)=14.4
	FYMAX(5)=1.5
	FYMAX(6)=1.125
	FYMAX(7)=1.71
	FYMAX(8)=1.14
	FYMAX(9)=3.486
	FYMAX(10)=0.75
	FYMAX(11)=1.5
	FYMAX(12)=2
	FYMAX(13)=0.81
	FYMAX(14)=0.9

C     FYMAX(1)=1.5
C     FYMAX(2)=1.5
C     FYMAX(3)=1.5
C     FYMAX(4)=1.5
C     FYMAX(5)=1.5
C     FYMAX(6)=1.5
C     FYMAX(7)=4.4

C
C FYMWX - Water in fresh manure (mm water / t manure)
C
	FYMWX(1)=0.075
	FYMWX(2)=0.075
	FYMWX(3)=0.07
	FYMWX(4)=0.04
	FYMWX(5)=0.075
	FYMWX(6)=0.075
	FYMWX(7)=0.094
	FYMWX(8)=0.094
	FYMWX(9)=0.094
	FYMWX(10)=0.0985
	FYMWX(11)=0.097
	FYMWX(12)=0.096
	FYMWX(13)=0.095
	FYMWX(14)=0.096

C     FYMWX(1)=0.077
C     FYMWX(2)=0.077
C     FYMWX(3)=0.077
C     FYMWX(4)=0.097
C     FYMWX(5)=0.077
C     FYMWX(6)=0.097
C     FYMWX(7)=0.094

C
C FYMLOSS Proportion of organic waste lost by volatn/timestep
C

      FYMLOSS(1)=0.05
      FYMLOSS(2)=0.05
      FYMLOSS(3)=0.15
      FYMLOSS(4)=0.15
      FYMLOSS(5)=0.05
      FYMLOSS(6)=0.0375
      FYMLOSS(7)=0.228
      FYMLOSS(8)=0.2
      FYMLOSS(9)=0.294
      FYMLOSS(10)=0.1
      FYMLOSS(11)=0.1
      FYMLOSS(12)=0.1
      FYMLOSS(13)=0.18
      FYMLOSS(14)=0.18

C     FYMLOSS(1)=0.15
C     FYMLOSS(2)=0.15
C     FYMLOSS(3)=0.15
C     FYMLOSS(4)=0.15
C     FYMLOSS(5)=0.15
C     FYMLOSS(6)=0.15
C     FYMLOSS(7)=0.15

C
C FYMPOS Proportion of FYM added to top 25cm 
C 
	FYMPOS(1)=1.0
	FYMPOS(2)=1.0
	FYMPOS(3)=1.0
	FYMPOS(4)=1.0
	FYMPOS(5)=1.0
	FYMPOS(6)=1.0
	FYMPOS(7)=0.8
	FYMPOS(8)=0.8
	FYMPOS(9)=0.8
	FYMPOS(10)=0.8
	FYMPOS(11)=0.8
	FYMPOS(12)=0.8
	FYMPOS(13)=0.8
	FYMPOS(14)=0.8

C     FYMPOS(1)=1.0
C     FYMPOS(2)=1.0
C     FYMPOS(3)=1.0
C     FYMPOS(4)=0.8
C     FYMPOS(5)=1.0
C     FYMPOS(6)=0.8
C     FYMPOS(7)=0.8

C
C FYMSTART Amount of rainfall in 1 week, below which volatilisation will occur
C
      FYMSTART(1)=5
      FYMSTART(2)=5
      FYMSTART(3)=5
      FYMSTART(4)=5
      FYMSTART(5)=5
      FYMSTART(6)=5
      FYMSTART(7)=5
      FYMSTART(8)=5
      FYMSTART(9)=5
      FYMSTART(10)=5
      FYMSTART(11)=5
      FYMSTART(12)=5
      FYMSTART(13)=5
      FYMSTART(14)=5
C
C Set default array descriptors
C
      I=1
      J=1
C
C Set default number of soils to number of soils previously parameterised
C
      NUMSOIL=3
C
C Read in parameters from SOIL_PAR.DAT
C
	REWIND(46)
303   CONTINUE
      READ(46,29,ERR=334,END=333)TEMP,I
      READ(46,*,ERR=334,END=333)FNULL,HY(I,1),FNULL,FNULL,
     &           BPART(I,1),BPROP(I,1),HPROP(I,1),BIOP(I,1),CRIT(I,1),				
     &           BIORATE(I,1),HUMRATE(I,1),TOCARRAY(I,1),IOMARRAY(I,1),			
     &           LUARRAY(I),CLARRAY(I,1),DEPTH,						
     &           PHARRAY(I,1),PHP1ARRAY(I,1),PHP2ARRAY(I,1)
	BIORATE(I,1)=BIORATE(I,1)/52
	HUMRATE(I,1)=HUMRATE(I,1)/52
C
C Fill up parameters in the soil profile
C
	TOCARRAY(I,1)=(MAXDEPTH/MAXLAYER1)*(TOCARRAY(I,1)/DEPTH)
	IOMARRAY(I,1)=(MAXDEPTH/MAXLAYER1)*(IOMARRAY(I,1)/DEPTH)
      DO 400 IL=2,MAXLAYER1
	  HY(I,IL)=HY(I,1)
	  BPART(I,IL)=BPART(I,1)
	  BPROP(I,IL)=BPROP(I,1)
	  HPROP(I,IL)=HPROP(I,1)
	  BIOP(I,IL)=BIOP(I,1)
	  CRIT(I,IL)=CRIT(I,1)
	  BIORATE(I,IL)=BIORATE(I,1)
	  HUMRATE(I,IL)=HUMRATE(I,1)
	  CLARRAY(I,IL)=CLARRAY(I,1)
	  PHARRAY(I,IL)=PHARRAY(I,1)
	  PHP1ARRAY(I,IL)=PHP1ARRAY(I,1)
  	  PHP2ARRAY(I,IL)=PHP2ARRAY(I,1)					
        IF(IL.LE.DEPTH*MAXLAYER/MAXDEPTH)THEN
	    TOCARRAY(I,IL)=TOCARRAY(I,1)
	    IOMARRAY(I,IL)=IOMARRAY(I,1)
	  ELSE
	    TOCARRAY(I,IL)=0
	    IOMARRAY(I,IL)=0
	  ENDIF
400   CONTINUE
C
C Set soil name and number of soils
C
      SNAME(I)=TEMP
      NUMSOIL=I
      GOTO 303
29    FORMAT(A40/I3)
      PRINT*,'SOIL: CN parameters from SOIL_PARS.DAT'
C
C Record error in the format of SOIL_PAR.DAT
C
334   CONTINUE
C
C If can open file SOIL.DAT, read in soil parameters from this
C
      REWIND(46)
      READ(46,10,ERR=222)TEMP,I
10    FORMAT(A40/I3)
	READ(46,20,ERR=222)FNULL,HY(I,1),BPART(I,1),BPROP(I,1),
     &                   HPROP(I,1),BIOP(I,1),CRIT(I,1),BIORATE(I,1),
     &                   HUMRATE(I,1),TOCARRAY(I,1),IOMARRAY(I,1),
     &                   LUARRAY(I)
      READ(46,30,ERR=222)CLARRAY(I,1),DEPTH,
     &                   PHARRAY(I,1),PHP1ARRAY(I,1),PHP2ARRAY(I,1),
     &                   BULKDENS(1)
20    FORMAT(11(F9.0/),I9)
30    FORMAT(5(F9.0/),F9.0)									
	BIORATE(I,1)=BIORATE(I,1)/52
	HUMRATE(I,1)=HUMRATE(I,1)/52
C
C Fill up parameters in the soil profile
C
	TOCARRAY(I,1)=(MAXDEPTH/MAXLAYER1)*(TOCARRAY(I,1)/DEPTH)
	IOMARRAY(I,1)=(MAXDEPTH/MAXLAYER1)*(IOMARRAY(I,1)/DEPTH)
      DO 401 IL=2,MAXLAYER1
	  HY(I,IL)=HY(I,1)
	  BPART(I,IL)=BPART(I,1)
	  BPROP(I,IL)=BPROP(I,1)
	  HPROP(I,IL)=HPROP(I,1)
	  BIOP(I,IL)=BIOP(I,1)
	  CRIT(I,IL)=CRIT(I,1)
	  BIORATE(I,IL)=BIORATE(I,1)
	  HUMRATE(I,IL)=HUMRATE(I,1)
	  CLARRAY(I,IL)=CLARRAY(I,1)
	  PHARRAY(I,IL)=PHARRAY(I,1)
	  PHP1ARRAY(I,IL)=PHP1ARRAY(I,1)
  	  PHP2ARRAY(I,IL)=PHP2ARRAY(I,1)	
        BULKDENS(IL)=BULKDENS(1)
        IF(IL.LE.DEPTH*MAXLAYER/MAXDEPTH)THEN
	    TOCARRAY(I,IL)=TOCARRAY(I,1)
	    IOMARRAY(I,IL)=IOMARRAY(I,1)
	  ELSE
	    TOCARRAY(I,IL)=0
	    IOMARRAY(I,IL)=0
	  ENDIF
401	CONTINUE
C
C Set soil name and number of soils
C
      SNAME(I)=TEMP
      NUMSOIL=I
C
      PRINT*,'SOIL: CN parameters from SOIL.DAT'
	RETURN
222   CONTINUE
      WRITE(*,*)'Warning! Error in general soil parameters!'
      WRITE(*,*)'Check format of soil parameter file, SOILLIS.OUT'
      WRITE(15,*)'Warning! Error in general soil parameters!'
      WRITE(15,*)'Check format of soil parameter file, SOILLIS.OUT'
333   CONTINUE
C
C Leave RESPAR
C
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE PROPN15(SOILN,SOIL15,TIN,AMMN,AMMN15,TAM)
C
C Subroutine to calculate the proportion of N15 added
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER,MAXLAYER1
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/
      INTEGER IL
C
C Variables passed to/from calling subroutine
C
      REAL SOILN(MAXLAYER),SOIL15(MAXLAYER),TIN(MAXLAYER)
	REAL AMMN(MAXLAYER),AMMN15(MAXLAYER),TAM(MAXLAYER)
C
C TIN()=proportion of N15:N in NO3 pool of the layer
C
      DO 4151 IL=1,MAXLAYER1
        IF(SOILN(IL).GT.0)THEN
	    TIN(IL)=SOIL15(IL)/SOILN(IL)
        ELSE
          TIN(IL)=0
        END IF
 4151 CONTINUE
C
C TAM()=proportion of N15:N in NH4 pool of the layer
C
        DO 4152 IL=1,MAXLAYER1
         IF(AMMN(IL).GT.0)THEN
          TAM(IL)=AMMN15(IL)/AMMN(IL)
         ELSE
          TAM(IL)=0
         END IF
 4152   CONTINUE
C
C Leave PROPN15
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,
     &             CONVER_F,SECONDS,
     &             NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)
C
C Subroutine to save internal SUNDIAL_CN parameters for use in other routines
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER INSTEPS				! Timesteps over which fert.is added(=3 weeks)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/
	INTEGER MAXSOIL				! Max.no.of soil types
	INTEGER MAXORGM				! Max.no.of org.manure applications allowed
	INTEGER MAXFERT				! Max.no.of fertiliser applications allowed 
	INTEGER MAXFERT1			! Max.no.of fertiliser applications allowed 
      PARAMETER (MAXSOIL=50,MAXORGM=52,MAXFERT=5)
	INTEGER I,J
	REAL ICFACTOR1(MAXLAYER)
C
C ...Timing factors
C
	REAL CONVER_F1
C
C ...Weather factors
C
	REAL BYRAIN1					! IN/OUT:Rain lost by bypass mm/timestep	
C
C ...Soil factors
C
      INTEGER NUMSOIL1				! Number of soils defined
	REAL HY1(MAXSOIL,MAXLAYER)		! IN/OUT:Stable N:C ratio of BIO&HUM pools
	REAL BPART1(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL HPART1(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL BPROP1(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from biomass decompositn
	REAL HPROP1(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from humus decompositn
	REAL ALPHA1(MAXLAYER)			! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL BETA1(MAXLAYER)			! IN/OUT: Prop.HUM produced on BIO decompn
	REAL GAMMA1(MAXLAYER)			! IN/OUT: Prop.BIO produced on HUM decompn 
	REAL DELTA1(MAXLAYER)			! IN/OUT: Prop.HUM produced on HUM decompn 
	REAL HZ11(MAXLAYER)				! IN/OUT: N:C ratio for steady state
	REAL BIOP1(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/TOC 
	REAL DPMRATE1(MAXLAYER)			! IN/OUT: Rate constant for DPM decompn
	REAL RPMRATE1(MAXLAYER)			! IN/OUT: Rate constant for RPM decompn
	REAL BRATE1(MAXLAYER)			! IN/OUT: Rate constant for HUM decompn
	REAL HRATE1(MAXLAYER)			! IN/OUT: Rate constant for BIO decompn
	REAL CRIT1(MAXSOIL,MAXLAYER)	! IN/OUT: Minimum nitrate level (Sozanska)
	REAL BIORATE1(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for BIO decompn/yr
	REAL HUMRATE1(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for HUM decompn/yr
	REAL IOMARRAY1(MAXSOIL,MAXLAYER)! IN/OUT: Equilibrium IOM (kgC/ha/layer)
	REAL TOCARRAY1(MAXSOIL,MAXLAYER)! IN/OUT: Equilibrium TOC (kgC/ha/layer)
	CHARACTER*40 SNAME1(MAXSOIL)	! IN/OUT:Soil name
	REAL DFACT1(MAXLAYER)			! IN/OUT: Denitrification factor
	REAL ANIT1						! Ammonium immobilised (kgN/ha)
	REAL ANIT151					! Ammonium N15 nitrified (kgN/ha)
	REAL FANIT1						! Fertiliser N nitrified (kgN/ha)
	REAL FANIT151					! Fertiliser N15 nitrified (kgN/ha)
	REAL PHARRAY1(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY1(MAXSOIL,MAXLAYER)! IN/PUT: pH below which decomp.zero
	REAL PHP2ARRAY1(MAXSOIL,MAXLAYER)! IN/PUT: pH above which decomp.max.
	REAL CLARRAY1(MAXSOIL,MAXLAYER)	! IN/OUT: % clay content in this layer
      INTEGER LUARRAY1(MAXSOIL)	! IN/OUT:Land use before equilibrium 
C
C ...Manure factors
C
      REAL FPROPX1(MAXORGM)			! IN/OUT:Prop.of FYM added to top 25cm
	REAL FYMCX1(MAXORGM)			! IN/OUT:Proportion of C in FYM
	REAL FYMNX1(MAXORGM)			! IN/OUT:Proportion of Organic N in FYM
	REAL FYMAX1(MAXORGM)			! IN/OUT:Proportion of ammonium-N in FYM
	REAL FYMWX1(MAXORGM)			! IN/OUT:Amount of water added in FYM
	REAL FYMLOSS1(MAXORGM)			! IN/OUT:Prop.FYM lost by volatn/timestep
	REAL FYMPOS1(MAXORGM)			! IN/OUT:Prop.of FYM added to top 25cm
	REAL FYMSTART1(MAXORGM)			! IN/OUT:Amount of rainfall in 1 week, 
									!		 below which volatilisation will occur
C
C Variables passed to/from calling subroutine
C ...Timing factors
C
	REAL SECONDS				! IN:Number of seconds in one timestep
	REAL CONVER_F
C
C ...Weather factors
C
	REAL BYRAIN					! IN/OUT:Rain lost by bypass flow mm/timestep	
C
C ...Soil factors
C
      INTEGER NUMSOIL					! Number of soils defined
      INTEGER ISAVE					! Code to save or retrieve variables
	REAL HY(MAXSOIL,MAXLAYER)		! IN/OUT:Stable N:C ratio of BIO&HUM pools
	REAL BPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL HPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL BPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from biomass decompositn
	REAL HPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from humus decompositn
	REAL ALPHA(MAXLAYER)			! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL BETA(MAXLAYER)				! IN/OUT: Prop.HUM produced on BIO decompn
	REAL GAMMA(MAXLAYER)			! IN/OUT: Prop.BIO produced on HUM decompn 
	REAL DELTA(MAXLAYER)			! IN/OUT: Prop.HUM produced on HUM decompn 
	REAL HZ1(MAXLAYER)				! IN/OUT: N:C ratio for steady state
	REAL BIOP(MAXSOIL,MAXLAYER)		! IN/OUT: BIO/TOC 
	REAL DPMRATE(MAXLAYER)			! IN/OUT: Rate constant for DPM decompn
	REAL RPMRATE(MAXLAYER)			! IN/OUT: Rate constant for RPM decompn
	REAL BRATE(MAXLAYER)			! IN/OUT: Rate constant for HUM decompn
	REAL HRATE(MAXLAYER)			! IN/OUT: Rate constant for BIO decompn
	REAL CRIT(MAXSOIL,MAXLAYER)		! IN/OUT: Minimum nitrate level (Sozanska)
	REAL BIORATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for BIO decompn/yr
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for HUM decompn/yr
	REAL IOMARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium IOM (kgC/ha/layer)
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium TOC (kgC/ha/layer)
	CHARACTER*40 SNAME(MAXSOIL)		! IN/OUT:Soil name
	REAL DFACT(MAXLAYER)			! IN/OUT: Denitrification factor
	REAL ANIT						! Ammonium immobilised (kgN/ha)
	REAL ANIT15						! Ammonium N15 nitrified (kgN/ha)
	REAL FANIT						! Fertiliser N nitrified (kgN/ha)
	REAL FANIT15					! Fertiliser N15 nitrified (kgN/ha)
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)! IN/PUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)! IN/PUT: pH above which decomp.max.
      REAL ICFACTOR(MAXLAYER)
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: % clay content in this layer
      INTEGER LUARRAY(MAXSOIL)	! IN/OUT:Land use before equilibrium 
C
C ...Manure factors
C
      REAL FPROPX(MAXORGM)			! IN/OUT:Prop.of FYM added to top 25cm
	REAL FYMCX(MAXORGM)				! IN/OUT:Prop.of C in FYM
	REAL FYMNX(MAXORGM)				! IN/OUT:Prop.of Organic N in FYM
	REAL FYMAX(MAXORGM)				! IN/OUT:Proportion of ammonium-N in FYM
	REAL FYMWX(MAXORGM)				! IN/OUT:Amount of water added in FYM
	REAL FYMLOSS(MAXORGM)			! IN/OUT:Prop.FYM lost by volatn/timestep
	REAL FYMPOS(MAXORGM)			! IN/OUT:Prop.of FYM added to top 25cm
	REAL FYMSTART(MAXORGM)			! IN/OUT:Amount of rainfall in 1 week, 
									!		 below which volatilisation will occur
C
C Save parameters
C
      SAVE
      INSTEPS=INT(3*(7*24*60*60)/SECONDS)
	IF(ISAVE.EQ.1)THEN
	  CONVER_F1=CONVER_F
	  ANIT1=ANIT
	  ANIT151=ANIT15
	  FANIT1=FANIT
	  FANIT151=FANIT15
	  BYRAIN1=BYRAIN
	  NUMSOIL1=NUMSOIL
	  DO 100 I=1,MAXSOIL
          SNAME1(I)=SNAME(I)
		LUARRAY1(I)=LUARRAY(I)
	    DO 200 J=1,MAXLAYER1
            HY1(I,J)=HY(I,J)
	      BPART1(I,J)=BPART(I,J)
	      HPART1(I,J)=HPART(I,J)
	      BPROP1(I,J)=BPROP(I,J)
	      HPROP1(I,J)=HPROP(I,J)
	      BIOP1(I,J)=BIOP(I,J)
	      CRIT1(I,J)=CRIT(I,J)
	      BIORATE1(I,J)=BIORATE(I,J)
	      HUMRATE1(I,J)=HUMRATE(I,J)
	      IOMARRAY1(I,J)=IOMARRAY(I,J)
	      TOCARRAY1(I,J)=TOCARRAY(I,J)
	      PHARRAY1(I,J)=PHARRAY(I,J)
	      PHP1ARRAY1(I,J)=PHP1ARRAY(I,J)
	      PHP2ARRAY1(I,J)=PHP2ARRAY(I,J)
		  CLARRAY1(I,J)=CLARRAY(I,J)
200       CONTINUE
100     CONTINUE	  
        DO 300 I=1,7
          FPROPX1(I)=FPROPX(I)
		FYMCX1(I)=FYMCX(I)
		FYMNX1(I)=FYMNX(I)
		FYMAX1(I)=FYMAX(I)
		FYMWX1(I)=FYMWX(I)
		FYMLOSS1(I)=FYMLOSS(I)
		FYMPOS1(I)=FYMPOS(I)
          FYMSTART1(I)=FYMSTART(I)
300     CONTINUE
        DO 500 I=1,MAXLAYER1
	    HZ11(I)=HZ1(I)
	    DFACT1(I)=DFACT(I)
          BRATE1(I)=BRATE(I)
	    HRATE1(I)=HRATE(I)
	    DPMRATE1(I)=DPMRATE(I)
	    RPMRATE1(I)=RPMRATE(I)
	    ALPHA1(I)=ALPHA(I)
	    BETA1(I)=BETA(I)
	    GAMMA1(I)=GAMMA(I)
	    DELTA1(I)=DELTA(I)
          ICFACTOR1(I)=ICFACTOR(I)
500     CONTINUE

C
C Retrieve parameters
C
	ELSEIF(ISAVE.EQ.0)THEN
	  CONVER_F=CONVER_F1
	  ANIT=ANIT1
	  ANIT15=ANIT151
	  FANIT=FANIT1
	  FANIT15=FANIT151
	  BYRAIN=BYRAIN1
	  NUMSOIL=NUMSOIL1
	  DO 700 I=1,MAXSOIL
          SNAME(I)=SNAME1(I)
		LUARRAY(I)=LUARRAY1(I)
	    DO 800 J=1,MAXLAYER1
            HY(I,J)=HY1(I,J)
	      BPART(I,J)=BPART1(I,J)
	      HPART(I,J)=HPART1(I,J)
	      BPROP(I,J)=BPROP1(I,J)
	      HPROP(I,J)=HPROP1(I,J)
	      BIOP(I,J)=BIOP1(I,J)
	      CRIT(I,J)=CRIT1(I,J)
	      BIORATE(I,J)=BIORATE1(I,J)
	      HUMRATE(I,J)=HUMRATE1(I,J)
	      IOMARRAY(I,J)=IOMARRAY1(I,J)
	      TOCARRAY(I,J)=TOCARRAY1(I,J)
	      PHARRAY(I,J)=PHARRAY1(I,J)
	      PHP1ARRAY(I,J)=PHP1ARRAY1(I,J)
	      PHP2ARRAY(I,J)=PHP2ARRAY1(I,J)
		  CLARRAY(I,J)=CLARRAY1(I,J)
800       CONTINUE
700     CONTINUE	  
        DO 900 I=1,7
          FPROPX(I)=FPROPX1(I)
		FYMCX(I)=FYMCX1(I)
		FYMNX(I)=FYMNX1(I)
		FYMAX(I)=FYMAX1(I)
		FYMWX(I)=FYMWX1(I)
		FYMLOSS(I)=FYMLOSS1(I)
		FYMPOS(I)=FYMPOS1(I)
	    FYMSTART(I)=FYMSTART1(I)
900     CONTINUE
        DO 1100 I=1,MAXLAYER1
          HZ1(I)=HZ11(I)
	    DFACT(I)=DFACT1(I)
          BRATE(I)=BRATE1(I)
	    HRATE(I)=HRATE1(I)
	    DPMRATE(I)=DPMRATE1(I)
	    RPMRATE(I)=RPMRATE1(I)
	    ALPHA(I)=ALPHA1(I)
	    BETA(I)=BETA1(I)
	    GAMMA(I)=GAMMA1(I)
	    DELTA(I)=DELTA1(I)
!          ICFACTOR(I)=ICFACTOR1(I) !Commented out so ICFACTOR can be change in dynamic spin up
1100    CONTINUE
	ENDIF
C
C Leave SAVE_CN
C
      END

C
C------------------------------------------------------------
C
      SUBROUTINE SAVE_FERTADD(SECONDS,FERTADD,FERTADD15,ISAVE)
C
C Subroutine to save fertiliser additions
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXFERTSTEPS			! Maximum number of timesteps 
									!	over which fertiliser is added
      PARAMETER (MAXFERTSTEPS=30*60*24*7*3)	! Timesteps = 3 wks (max.30 secs)
	INTEGER INSTEPS					! Timesteps over which fert.is added(=3 weeks)
	INTEGER I						! Local counter
	INTEGER J						! Local counter
	REAL FERTADD1(MAXFERTSTEPS,4)	! INOUT:Fertiliser N added in this time 
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
	REAL FERTADD151(MAXFERTSTEPS,4)	! INOUT:Fertiliser N15 added in this time
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
C
C Variables passed to/from calling subroutine
C
	REAL FERTADD(MAXFERTSTEPS,4)	! INOUT:Fertiliser N added in this time 
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
	REAL FERTADD15(MAXFERTSTEPS,4)	! INOUT:Fertiliser N15 added in this time
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
      INTEGER ISAVE				! Code to save or retrieve variables
	REAL SECONDS				! IN:Number of seconds in one timestep
C
C Save parameters
C
      SAVE
      INSTEPS=INT(3*(7*24*60*60)/SECONDS)
	IF(INSTEPS.EQ.0)INSTEPS=1											
	IF(ISAVE.EQ.1)THEN
	  DO 400 I=1,INSTEPS
	    DO 450 J=1,4
	      FERTADD1(I,J)=FERTADD(I,J)
	      FERTADD151(I,J)=FERTADD15(I,J)
450     CONTINUE
400     CONTINUE

C
C Retrieve parameters
C
	ELSEIF(ISAVE.EQ.0)THEN
	  DO 1000 I=1,INSTEPS
	    DO 1050 J=1,4
	     FERTADD(I,J)=FERTADD1(I,J)
	     FERTADD15(I,J)=FERTADD151(I,J)
1050      CONTINUE
1000    CONTINUE
	ENDIF
C
C Leave SAVE_FERTADD
C
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE SETDOC(NSOIL,MOBDOC,MOBDON)
C
C Subroutine to calculate denitrification
C
      IMPLICIT NONE
C
C Variables local to this routine
C
	INTEGER MAXLAYER,MAXLAYER1
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/
      INTEGER IL
C
C Variables passed to/from this routine
C
      REAL MOBDOC(MAXLAYER)
      REAL MOBDON(MAXLAYER)
	INTEGER NSOIL
C
C Initialise mobile DOC down profile 
C Assumed to be 0 - correct when more information available
C
      DO 100 IL=1,MAXLAYER1
	  MOBDOC(IL)=0
	  MOBDON(IL)=0
100   CONTINUE
      END
C
C----------------------------------------------------------
C
      SUBROUTINE SETFACTORS(NSOIL,LUCODE,INPUTFACTORS,ATM)
C
C Subroutine to set BIO and HUM decompostion parameters
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXORGM				! Max.no.of org.manure applications allowed
	PARAMETER(MAXORGM=52)
	INTEGER MAXSOIL				! Max.no.of soil types
      PARAMETER (MAXSOIL=50)
	INTEGER IL
	REAL PH(MAXLAYER),STABCN
	INTEGER INPUTFACTORS			! Input factors for decomposition process. If = 1, asks user for value
	REAL INPUTATMN					! Input atmospheric N input: Def = 0 kg N / ha
	REAL INPUTDECEFF				! Input decomposition efficiency: Def = 46%
	REAL INPUTSTABLECN				! Input stable C:N: Def=8
C	
C Variables passed to/from calling subroutine
C
	REAL ALPHA(MAXLAYER)			! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL ANIT						! Ammonium immobilised (kgN/ha)
	REAL ANIT15						! Ammonium N15 nitrified (kgN/ha)
      REAL ATM						! IN/OUT:Atmospheric N input (kgN/ha/timestep)
	REAL BETA(MAXLAYER)				! IN/OUT: Prop.HUM produced on BIO decompn
	REAL BIOP(MAXSOIL,MAXLAYER)		! IN/OUT: BIO/TOC 
	REAL BIORATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for BIO decompn/yr
	REAL BPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL BPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from biomass decompositn
	REAL BRATE(MAXLAYER)			! IN/OUT: Rate constant for HUM decompn
	REAL BYRAIN						! IN/OUT:Rain lost by bypass flow mm/timestep	
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: % clay content in this layer
      REAL CONVER_F					! IN/OUT:Conversion to correct timestep
	REAL CRIT(MAXSOIL,MAXLAYER)		! IN/OUT: Minimum nitrate level (Sozanska)
	REAL DELTA(MAXLAYER)			! IN/OUT: Prop.HUM produced on HUM decompn 
	REAL DFACT(MAXLAYER)			! IN/OUT: Denitrification factor
	REAL DPMRATE(MAXLAYER)			! IN/OUT: Rate constant for DPM decompn
	REAL FANIT						! Fertiliser N nitrified (kgN/ha)
	REAL FANIT15					! Fertiliser N15 nitrified (kgN/ha)
      REAL FPROPX(MAXORGM)			! IN/OUT:Proportion of FYM added to top 25cm
	REAL FYMAX(MAXORGM)				! IN/OUT:Proportion of ammonium-N in FYM
	REAL FYMCX(MAXORGM)				! IN/OUT:Proportion of C in FYM
	REAL FYMLOSS(MAXORGM)			! IN/OUT:Prop.FYM lost by volatn/timestep
	REAL FYMNX(MAXORGM)				! IN/OUT:Proportion of Organic N in FYM
	REAL FYMPOS(MAXORGM)			! IN/OUT:Proportion of FYM added to top 25cm
      REAL FYMSTART(MAXORGM)			! IN/OUT:Amount of rainfall in 1 week, 
									!		 below which volatilisation will occur
	REAL FYMWX(MAXORGM)				! IN/OUT:Amount of water added in FYM
	REAL GAMMA(MAXLAYER)			! IN/OUT: Prop.BIO produced on HUM decompn 
	REAL HPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL HPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from humus decompositn
	REAL HRATE(MAXLAYER)			! IN/OUT: Rate constant for BIO decompn
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for HUM decompn/yr
	REAL HY(MAXSOIL,MAXLAYER)		! IN/OUT:Stable N:C ratio of BIO & HUM pools
	REAL HZ1(MAXLAYER)				! IN/OUT: N:C ratio for steady state
  	REAL ICFACTOR(MAXLAYER)			! IN:Adjustment in rate needed to	
									!    achieve measured NPP and TOC	
	REAL IOMARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium IOM (kgC/ha/layer)
      INTEGER ISAVE					! Code to save or retrieve variables
      INTEGER LUARRAY(MAXSOIL)		! IN/OUT:Land use before equilibrium 
	INTEGER LUCODE					! IN: Land use code
	INTEGER NSOIL					! IN:Soil code number
      INTEGER NUMSOIL					! Number of soils defined
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: pH above which decomp.max.
	REAL RPMRATE(MAXLAYER)			! IN/OUT: Rate constant for RPM decompn
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium TOC (kgC/ha/layer)
	REAL SECONDS					! IN:Number of seconds in one timestep
 	CHARACTER*40 SNAME(MAXSOIL)		! IN/OUT:Soil name
C
C Retrieve soil CN characteristics 
C
      ISAVE=0
      CALL SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,CONVER_F,
     &             SECONDS,NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)
C
C Set factors according to input values
C ......get set atmospheric N input
C
      IF(INPUTFACTORS.EQ.1)THEN
	  CALL GETSTABLECN(PHARRAY(NSOIL,1),STABCN)
	  WRITE(*,10)LUCODE,STABCN
10      FORMAT('LAND USE ',I3,
     &         ' Enter stable C:N ratio of BIO & HUM'/
     &         ' (Default =  ',F5.1,' Range = 4-30.',
     &         ' Enter -1 to use default): ')
	  READ(*,*)INPUTSTABLECN
	  WRITE(*,20)BPART(NSOIL,1)*100
20      FORMAT(/' Enter efficiency of decomp. (%)'/
     &          ' (Default = ',F5.1,' Range = 15-100.',
     &          ' Enter -1 to use default): ')
	  READ(*,*)INPUTDECEFF
	  INPUTDECEFF=INPUTDECEFF/100.
	  WRITE(*,30)
30      FORMAT(/' Enter atm N input (kg/ha/year)',
     &          ' Range = 0-50. Enter -1 to use default): ')
	  READ(*,*)INPUTATMN
	  IF(INPUTATMN.GE.0)ATM=INPUTATMN
	  WRITE(*,40)
40      FORMAT(/' ====================================================',
     &           '=='/)
	ENDIF
C
C For each layer down the profile set the parameters
C
      DO 100 IL=1,MAXLAYER
C
C Set N:C Ratio of N in HUM and ROOT pools
C
	  PH(IL)=PHARRAY(NSOIL,IL)
	  CALL GETSTABLECN(PH(IL),STABCN)
	  IF(INPUTFACTORS.EQ.1.AND.INPUTSTABLECN.GT.0)STABCN=INPUTSTABLECN
	  HY(NSOIL,IL)=1/STABCN
C
C BPART and HPART = proportion of Carbon going to BIO + HUM (=1-EFFICIENCY)
C
        IF(INPUTFACTORS.EQ.1.AND.INPUTDECEFF.GT.0)
     &     BPART(NSOIL,IL)=INPUTDECEFF
        HPART(NSOIL,IL)=BPART(NSOIL,IL)
C
C ALPHA = Proportion of decomposing biomass that is transformed to biomass
C BETA = Proportion of decomposing biomass that is transformed to humus
C 
        BETA(IL)=BPART(NSOIL,IL)/(1+BPROP(NSOIL,IL))
        ALPHA(IL)=BPART(NSOIL,IL)-BETA(IL)
		
C
C GAMMA = Proportion of decomposing humus that is transformed to biomass
C DELTA = Proportion of decomposing humus that is transformed to humus
C
        DELTA(IL)=HPART(NSOIL,IL)/(1+HPROP(NSOIL,IL))
        GAMMA(IL)=HPART(NSOIL,IL)-DELTA(IL)
C
C Calc. N:C of input plant material for steady state conditions = HZ1
C
        HZ1(IL)=HY(NSOIL,IL)*(ALPHA(IL)+BETA(IL))
100   CONTINUE
C
C Save soil CN characteristics 
C
      ISAVE=1
      CALL SAVE_CN(HY,BPART,HPART,BPROP,HPROP,
     &             ALPHA,BETA,GAMMA,DELTA,
     &             HZ1,BIOP,DPMRATE,RPMRATE,BRATE,HRATE,CONVER_F,
     &             SECONDS,NUMSOIL,CRIT,BIORATE,HUMRATE,
     &             IOMARRAY,TOCARRAY,SNAME,
     &             FPROPX,FYMCX,FYMNX,FYMAX,FYMWX,
     &             FYMLOSS,FYMPOS,FYMSTART,BYRAIN,
     &             ANIT,ANIT15,FANIT,FANIT15,DFACT,ISAVE,
     &             PHARRAY,PHP1ARRAY,PHP2ARRAY,ICFACTOR,
     &             CLARRAY,LUARRAY)
C
C Leave SETFACTORS
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE SETORGM_FIXED(NSOIL,DPMRATE,RPMRATE,TOC,IOM,CONVER_F,
     &                   BRATE,BIORATE,HRATE,HUMRATE,HZ1,HY,
     &                   BCARB0,HCARB0,DPMCARB0,RPMCARB0,
     &                   BNIT0,HNIT0,DPMNIT0,RPMNIT0,
     &                   HNLAB0,BNLAB0,DPMNLAB0,RPMNLAB0,ICFACTOR,
     &                   DPM_RPM_RATIO,DPMCTON,RPMCTON)	
C
C Subroutine to set the residual N and C in the RO,BIO and HUM pools
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	INTEGER MAXDEPTH			! Maximum depth of the soil profile
	INTEGER MAXSOMLAY			! Maximum number of SOM layers
	PARAMETER (MAXLAYER=60)								
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
	INTEGER MAXSOIL				! Max.no.of soil types
      PARAMETER (MAXSOIL=50)

	REAL ACT_OM_FIXED			! Fixed amount of active organic matter
								! (used to adjust ratios)
	REAL BIO_OM_FIXED			! Fixed amount of biomass (used to get ratios)
	REAL DPM_OM_FIXED			! Fixed amount of DPM (used to get ratios)
      REAL DPM_PLUS_RPM			! Fixed C in DPM and RPM (used to get ratios)
	REAL HUM_OM_FIXED			! Fixed amount of humus (used to get ratios)
	INTEGER IL					! Local counter for layers
	INTEGER IMON				! Local counter for months
	REAL RPM_OM_FIXED			! Fixed amount of RPM (used to get ratios)

	DATA ACT_OM_FIXED,DPM_OM_FIXED,RPM_OM_FIXED,BIO_OM_FIXED,
     &     HUM_OM_FIXED /90.71,0.81,13.21,1.95,74.74/
C
C Variables passed to/from calling subroutine
C
  	REAL BCARB0(MAXLAYER)		! IN:C in soil biomass (kgC/ha/layer)
	REAL BIORATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for BIO decompn/yr
      REAL BNIT0(MAXLAYER)		! IN:N in soil biomass (kgN/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN:N15 in soil humus (kgN15/ha/layer)
	REAL BORGM					! IN:C in soil biomass (kgC/ha/layer)
	REAL BRATE(MAXLAYER)		! IN/OUT: Rate constant for HUM decompn
      REAL CONVER_F				! IN/OUT:Conversion to correct timestep
	REAL DORGM					! IN:C in soil DPM (kgC/ha/layer)
      REAL DPM_RPM_RATIO			! OUT:Ratio of DPM:RPM. If set to 
								!     zero, this will be worked out using a call to SET_DPMRPMRATIO
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
      REAL DPMCTON(MAXLAYER)		! IN:Ratio of C:N in DPM
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL DPMRATE(MAXLAYER)		! IN/OUT: Rate constant for DPM decompn
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN:N in soil humus (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN:N15 in soil biomass (kgN15/ha/layer)
	REAL HORGM					! IN:C in soil humus (kgC/ha/layer)
	REAL HRATE(MAXLAYER)		! IN/OUT: Rate constant for BIO decompn
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for HUM decompn/yr
	REAL HY(MAXSOIL,MAXLAYER)	! IN/OUT:Stable N:C ratio of BIO & HUM pools
	REAL HZ1(MAXLAYER)			! IN/OUT: N:C ratio for steady state
	REAL ICFACTOR(MAXLAYER)		! IN:Adjustment in rate needed to 
								!    achieve measured NPP and TOC
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
	INTEGER NSOIL				! IN:Soil code number
	REAL RORGM					! IN:C in soil RPM (kgC/ha/layer)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
      REAL RPMCTON(MAXLAYER)		! IN:Ratio of C:N in RPM
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL RPMRATE(MAXLAYER)		! IN/OUT: Rate constant for RPM decompn
	REAL TOC(MAXLAYER)			! IN:Total organic C (kgC/ha/layer)
C
C Set Rate constants for decomposition of RO (=R), BIO (=B) and HUM (=H)
C
	CALL GET_DECOMP_RATECONSTS(CONVER_F,NSOIL,BIORATE,HUMRATE,
     &                           DPMRATE,RPMRATE,BRATE,HRATE)	  
C
C Set the amounts in initial organic matter
C of biomass(=BORGM) and humus (=HORGM)
C and partition down the soil profile
C
      DO 100 IL=1,MAXLAYER1
C
C Set the amounts in initial organic matter
C ...of biomass(=BORGM) and humus (=HORGM)
C
        BORGM = (TOC(IL)-IOM(IL))*BIO_OM_FIXED/ACT_OM_FIXED
	  HORGM = (TOC(IL)-IOM(IL))*HUM_OM_FIXED/ACT_OM_FIXED
        BCARB0(IL)=BORGM
        HCARB0(IL)=HORGM
C
C ...and of DPM and RPM
C
        IF(DPM_RPM_RATIO.LT.0.0001)THEN
	    DPMCARB0(IL)=(TOC(IL)-IOM(IL))*DPM_OM_FIXED/ACT_OM_FIXED
	    RPMCARB0(IL)=(TOC(IL)-IOM(IL))*RPM_OM_FIXED/ACT_OM_FIXED
	  ELSEIF(DPM_RPM_RATIO.GE.0.0001)THEN
	    DPM_PLUS_RPM=(DPM_OM_FIXED+RPM_OM_FIXED)/ACT_OM_FIXED
          RPMCARB0(IL)=(TOC(IL)-IOM(IL))*DPM_PLUS_RPM/(1+DPM_RPM_RATIO)
	    DPMCARB0(IL)=(TOC(IL)-IOM(IL))*(DPM_PLUS_RPM-RPMCARB0(IL))
	  ENDIF
C
C Set the initial N15 in organic pools to 0
C
        HNLAB0(IL)=0
        BNLAB0(IL)=0
        DPMNLAB0(IL)=0
        RPMNLAB0(IL)=0
C
C Set the amount of nitrogen in BIO and HUM pools from N:C ratio of BIO&HUM
C
        BNIT0(IL)=BCARB0(IL)*HY(NSOIL,IL)
        HNIT0(IL)=HCARB0(IL)*HY(NSOIL,IL)
C
C ... and of DPM and RPM from passed in ratio
C
        DPMNIT0(IL)=DPMCARB0(IL)/DPMCTON(IL)
        RPMNIT0(IL)=RPMCARB0(IL)/RPMCTON(IL)
C
C Set ICFACTOR to 1
C
        ICFACTOR(IL)=1
100   CONTINUE
C
C Leave SETORGM
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE SETORGM_EQRUN(NSOIL,DPMRATE,RPMRATE,TOC,IOM,CONVER_F,
     &                   BRATE,BIORATE,HRATE,HUMRATE,HZ1,HY,
     &                   BCARB0,HCARB0,DPMCARB0,RPMCARB0,
     &                   BNIT0,HNIT0,DPMNIT0,RPMNIT0,
     &                   HNLAB0,BNLAB0,DPMNLAB0,RPMNLAB0,
     &                   LUARRAY,CLARRAY,
     &                   PHARRAY,PHP1ARRAY,PHP2ARRAY,
     &                   EQMODEL,PI_CEQ_MON,ICFACTOR,LTA_AWC,LTA_TEMP,
     &                   ITFUNC,IMFUNC,WMAX,WSAT,
     &                   DPM_RPM_RATIO,DPMCTON,RPMCTON,ALPHA,BETA,CLAY)	
C
C Subroutine to set the residual N and C in the RO,BIO and HUM pools
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	INTEGER MAXDEPTH			! Maximum depth of the soil profile
	INTEGER MAXSOMLAY			! Maximum number of SOM layers
	PARAMETER (MAXLAYER=60)								
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
	INTEGER MAXSOIL				! Max.no.of soil types
      PARAMETER (MAXSOIL=50)
	INTEGER MAXLU				! Max no. of land use types
      PARAMETER (MAXLU=6)
	REAL AWC_IL(12)				! OUT:Long term ave.soil moist.def.(mm)
	REAL DEPTH					! Depth of this layer
	INTEGER IL					! Local counter for layers
	INTEGER IMON				! Local counter for months
	REAL PI_IL(12)              ! IN/OUT: Equilibrium plant C input each month
	REAL TEMP_IL(12)			! OUT:Long term ave.soil moist.def.(mm)   
	REAL THISPI(12)			    ! Equilibrium plant C input each month 
								! to this layer (t C / ha / month)
C
C Variables passed to/from calling subroutine
C
	REAL ALPHA(MAXLAYER)			! IN/OUT: Prop.BIO produced on BIO decompn 
  	REAL BCARB0(MAXLAYER)		! IN:C in soil biomass (kgC/ha/layer)
	REAL BETA(MAXLAYER)				! IN/OUT: Prop.HUM produced on BIO decompn
	REAL BIORATE(MAXSOIL,MAXLAYER)	! IN: Rate constant for BIO decompn/yr
      REAL BNIT0(MAXLAYER)		! IN:N in soil biomass (kgN/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN:N15 in soil humus (kgN15/ha/layer)
	REAL BORGM					! IN:C in soil biomass (kgC/ha/layer)
	REAL BRATE(MAXLAYER)		! IN: Rate constant for HUM decompn
      REAL CLAY(MAXLAYER)			! OUT:Clay content of the layer (%)
      REAL CLAY_IL				! IN: Clay content of the layer (%)
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN: % clay content in this layer
      REAL CONVER_F				! IN:Conversion to correct timestep
	REAL DORGM					! IN:C in soil DPM (kgC/ha/layer)
      REAL DPM_RPM_RATIO			! IN/OUT:Ratio of DPM:RPM. If set to 
								!        zero, this will be worked out using a call to SET_DPMRPMRATIO
	REAL DPMCARB0(MAXLAYER)		! IN/OUT:C in decomposable PM (kgC/ha/layer)
      REAL DPMCTON(MAXLAYER)		! IN:Ratio of C:N in DPM
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL DPMRATE(MAXLAYER)		! IN: Rate constant for DPM decompn
	INTEGER EQMODEL				! IN:Type of equilibrium run (NPP,TOC or both)
	INTEGER EQNPP				! Model initialised using plant inputs 
								!    from measured NPP
	INTEGER EQTOC				! Model initialised using measured TOC
	INTEGER EQNPPTOC			! Model initialised using both plant inputs
								!    from measured NPP and measured TOC
	INTEGER EQHILLIER,EQJONES			! Model initialised using TOC and
								!    the Hillier solver to get C pools and plant inputs
	DATA EQNPP,EQTOC,EQNPPTOC,EQHILLIER,EQJONES /1,2,3,4,5/
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN:N in soil humus (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN:N15 in soil biomass (kgN15/ha/layer)
	REAL HORGM					! IN:C in soil humus (kgC/ha/layer)
	REAL HRATE(MAXLAYER)		! IN: Rate constant for BIO decompn
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN: Rate constant for HUM decompn/yr
	REAL HY(MAXSOIL,MAXLAYER)	! IN:Stable N:C ratio of BIO & HUM pools
	REAL HZ1(MAXLAYER)			! IN: N:C ratio for steady state
	REAL ICFACTOR(MAXLAYER)		! IN:Adjustment in rate needed to 
								!    achieve measured NPP and TOC
	INTEGER IMFUNC				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)     
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
	INTEGER ITFUNC				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	REAL LTA_AWC(12,MAXLAYER)	! IN:Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)	! IN:Average air temp this month (deg.C)

      INTEGER LUARRAY(MAXSOIL)	! IN:Land use before equilibrium 
	INTEGER NSOIL				! IN:Soil code number
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)	! IN/PUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)	! IN/PUT: pH above which decomp.max.
	REAL PI_CEQ_MON(12,MAXLAYER)! IN/OUT: Equilibrium plant C input each month
	REAL RORGM					! IN:C in soil RPM (kgC/ha/layer)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
      REAL RPMCTON(MAXLAYER)				! IN/OUT:Ratio of C:N in RPM
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL RPMRATE(MAXLAYER)		! IN: Rate constant for RPM decompn
	REAL TOC(MAXLAYER)			! IN:Total organic C (kgC/ha/layer)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)

C
C Set Rate constants for decomposition of RO (=R), BIO (=B) and HUM (=H)
C
	CALL GET_DECOMP_RATECONSTS(CONVER_F,NSOIL,BIORATE,HUMRATE,
     &                           DPMRATE,RPMRATE,BRATE,HRATE)	  
C
C
C Set the amounts in initial organic matter
C of biomass(=BORGM) and humus (=HORGM)
C and partition down the soil profile
C
      DO 100 IL=1,MAXLAYER1
C
C Set clay content
C
	  CLAY(IL)=CLARRAY(NSOIL,IL)

C Get depth of this layer
C
        DEPTH=IL*MAXDEPTH/(MAXLAYER1*1.)
C
C Set the amounts in initial organic matter
C of biomass(=BORGM) and humus (=HORGM)
C
        
        DO 200 IMON=1,12
	    THISPI(IMON)=PI_CEQ_MON(IMON,IL)/1000
200     CONTINUE
        IF(EQMODEL.EQ.EQNPP.OR.
     &     EQMODEL.EQ.EQTOC.OR.
     &     EQMODEL.EQ.EQNPPTOC)THEN
          CALL GETEQOM_FROM_ROTHC(TOC(IL),IOM(IL),
     &             BORGM,HORGM,RORGM,DORGM,
     &             CLARRAY(NSOIL,IL),DEPTH,
     &             PHARRAY(NSOIL,IL),PHP1ARRAY(NSOIL,IL),
     &             PHP2ARRAY(NSOIL,IL),
     &             EQMODEL,THISPI,ICFACTOR(IL),LTA_AWC,LTA_TEMP,
     &		     ITFUNC,IMFUNC,WMAX,WSAT,DPM_RPM_RATIO)       
	  ELSEIF(EQMODEL.EQ.EQHILLIER)THEN
          DO 250 IMON=1,12
		  AWC_IL(IMON)=LTA_AWC(IMON,IL)
	      TEMP_IL(IMON)=LTA_TEMP(IMON,IL)
		  PI_IL(IMON)=PI_CEQ_MON(IMON,IL)
250       CONTINUE
	    CALL HILLIER_SOLVER(
C Input 
     &                   ALPHA(IL),AWC_IL,BETA(IL),CLAY(IL),
     &                   DPM_RPM_RATIO,ITFUNC,IMFUNC,
     &                   DPMRATE(IL),RPMRATE(IL),BRATE(IL),HRATE(IL),
     &                   TOC(IL),PHARRAY(NSOIL,IL),
     &                   PHP1ARRAY(NSOIL,IL),PHP2ARRAY(NSOIL,IL),
     &                   TEMP_IL,WMAX(IL),WSAT(IL),
C Input / Output
     &                   PI_IL,ICFACTOR(IL),
C Output
     &                   BORGM,DORGM,HORGM,RORGM,IL)
          DO 260 IMON=1,12
	      THISPI(IMON)=PI_IL(IMON)/1000
260       CONTINUE

         
          ELSEIF(EQMODEL.EQ.EQJONES)THEN
      
            DO 251 IMON=1,12
		  AWC_IL(IMON)=LTA_AWC(IMON,IL)
	      TEMP_IL(IMON)=LTA_TEMP(IMON,IL)
		  PI_IL(IMON)=PI_CEQ_MON(IMON,IL)
251         CONTINUE
            

	DPMRATE=52*DPMRATE/CONVER_F !rates per year
	RPMRATE=52*RPMRATE/CONVER_F
	BRATE=52*BRATE/CONVER_F
	HRATE=52*HRATE/CONVER_F

 	CALL ED_SOLVER_1(
C Input 
     &                   ALPHA(IL),AWC_IL,BETA(IL),CLAY(IL),
     &                   DPM_RPM_RATIO,ITFUNC,IMFUNC,
     &                   DPMRATE(IL),RPMRATE(IL),BRATE(IL),HRATE(IL),
     &                   TOC(IL),PHARRAY(NSOIL,IL),
     &                   PHP1ARRAY(NSOIL,IL),PHP2ARRAY(NSOIL,IL),
     &                   TEMP_IL,WMAX(IL),WSAT(IL),
C Input / Output
     &                   PI_IL,ICFACTOR(IL),
C Output
     &                   BORGM,DORGM,HORGM,RORGM)

	DPMRATE=DPMRATE*CONVER_F/52 !convert back
	RPMRATE=RPMRATE*CONVER_F/52
	BRATE=BRATE*CONVER_F/52
	HRATE=HRATE*CONVER_F/52
! ICFACTOR(IL)=12*ICFACTOR(IL)*CONVER_F/52

          DO 261 IMON=1,12
	      THISPI(IMON)=PI_IL(IMON)/1000
261       CONTINUE
	  ENDIF
C
C Partition the biomass, humus, DPM and RPM in the SOM layer 
C between the equal soil layers
C
        BCARB0(IL)=BORGM
        HCARB0(IL)=HORGM
        DPMCARB0(IL)=DORGM											
        RPMCARB0(IL)=RORGM											
C
C Set the initial N15 in organic pools to 0
C
        HNLAB0(IL)=0
        BNLAB0(IL)=0
        DPMNLAB0(IL)=0
        RPMNLAB0(IL)=0
C
C Set the amount of nitrogen in BIO and HUM pools from N:C ratio of BIO&HUM
C
        BNIT0(IL)=BCARB0(IL)*HY(NSOIL,IL)
        HNIT0(IL)=HCARB0(IL)*HY(NSOIL,IL)

C
C C/N ratio of DPM and RPM changed to be more in equilibrium with those during simulation
C
        DPMNIT0(IL)=DPMCARB0(IL)/DPMCTON(IL)
        RPMNIT0(IL)=RPMCARB0(IL)/RPMCTON(IL)
100   CONTINUE
C
C Leave SETORGM
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE SETORGM_ROTHC(NSOIL,DPMRATE,RPMRATE,TOC,IOM,CONVER_F,
     &                   BRATE,BIORATE,HRATE,HUMRATE,HZ1,HY,
     &                   BCARB0,HCARB0,DPMCARB0,RPMCARB0,
     &                   BNIT0,HNIT0,DPMNIT0,RPMNIT0,
     &                   HNLAB0,BNLAB0,DPMNLAB0,RPMNLAB0,
     &                   LUARRAY,CLARRAY,
     &                   PHARRAY,PHP1ARRAY,PHP2ARRAY,
     &                   EQMODEL,PI_CEQ_MON,ICFACTOR,LTA_AWC,LTA_TEMP,
     &                   ITFUNC,IMFUNC,WMAX,WSAT,DPM_RPM_RATIO)
C
C Subroutine to set the residual N and C in the RO,BIO and HUM pools
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	INTEGER MAXDEPTH			! Maximum depth of the soil profile
	INTEGER MAXSOMLAY			! Maximum number of SOM layers
	PARAMETER (MAXLAYER=60)								
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
	INTEGER MAXSOIL				! Max.no.of soil types
      PARAMETER (MAXSOIL=50)
	INTEGER MAXLU				! Max no. of land use types
      PARAMETER (MAXLU=6)
	REAL DEPTH					! Depth of this layer
	INTEGER IL					! Local counter for layers
	INTEGER IMON				! Local counter for months
	REAL THISPI(12)			    ! Equilibrium plant C input each month 
								! to this layer (t C / ha / month)
C
C Variables passed to/from calling subroutine
C ...Model descriptors
C
	INTEGER ITFUNC				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IMFUNC				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)     
C
C ...Weather factors
C
	REAL LTA_AWC(12,MAXLAYER)	! IN:Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)	! IN:Average air temp this month (deg.C)
C
C ...Soil descriptors
C
      REAL DPM_RPM_RATIO			! IN/OUT:Ratio of DPM:RPM. If set to 
								!        zero, this will be worked out using a call to SET_DPMRPMRATIO
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
	REAL BORGM					! IN:C in soil biomass (kgC/ha/layer)
	REAL HORGM					! IN:C in soil humus (kgC/ha/layer)
	REAL RORGM					! IN:C in soil RPM (kgC/ha/layer)
	REAL DORGM					! IN:C in soil DPM (kgC/ha/layer)
	INTEGER NSOIL				! IN:Soil code number
	REAL DPMRATE(MAXLAYER)		! IN: Rate constant for DPM decompn
	REAL RPMRATE(MAXLAYER)		! IN: Rate constant for RPM decompn
	REAL BRATE(MAXLAYER)		! IN: Rate constant for HUM decompn
	REAL HRATE(MAXLAYER)		! IN: Rate constant for BIO decompn
      REAL CONVER_F				! IN:Conversion to correct timestep
	REAL HZ1(MAXLAYER)			! IN: N:C ratio for steady state
	REAL HY(MAXSOIL,MAXLAYER)	! IN:Stable N:C ratio of BIO & HUM pools
  	REAL BCARB0(MAXLAYER)		! IN:C in soil biomass (kgC/ha/layer)
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
      REAL BNIT0(MAXLAYER)		! IN:N in soil biomass (kgN/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN:N in soil humus (kgN/ha/layer)
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN:N15 in soil biomass (kgN15/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN:N15 in soil humus (kgN15/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL TOC(MAXLAYER)			! IN:Total organic C (kgC/ha/layer)
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
	REAL BIORATE(MAXSOIL,MAXLAYER)	! IN: Rate constant for BIO decompn/yr
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN: Rate constant for HUM decompn/yr
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN: % clay content in this layer
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)	! IN/PUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)	! IN/PUT: pH above which decomp.max.
      INTEGER LUARRAY(MAXSOIL)	! IN:Land use before equilibrium 
	INTEGER EQMODEL				! IN:Type of equilibrium run (NPP, TOC or both)
	INTEGER EQNPP				! Model initialised using plant inputs 
								!    from measured NPP
	INTEGER EQTOC				! Model initialised using measured TOC
	INTEGER EQNPPTOC			! Model initialised using both plant inputs
								!    from measured NPP and measured TOC
	INTEGER EQHILLIER,EQJONES			! Model initialised using TOC and
								!    the Hillier solver to get C pools and plant inputs
	DATA EQNPP,EQTOC,EQNPPTOC,EQHILLIER,EQJONES /1,2,3,4,5/
	REAL PI_CEQ_MON(12,MAXLAYER)! IN/OUT: Equilibrium plant C input each month
	REAL ICFACTOR(MAXLAYER)		! IN:Adjustment in rate needed to 
								!    achieve measured NPP and TOC
C
C Set Rate constants for decomposition of RO (=R), BIO (=B) and HUM (=H)
C
	CALL GET_DECOMP_RATECONSTS(CONVER_F,NSOIL,BIORATE,HUMRATE,
     &                           DPMRATE,RPMRATE,BRATE,HRATE)	  
C
C Set the amounts in initial organic matter
C of biomass(=BORGM) and humus (=HORGM)
C and partition down the soil profile
C
      DO 100 IL=1,MAXLAYER1
C
C Get depth of this layer
C
        DEPTH=IL*MAXDEPTH/(MAXLAYER1*1.)
C
C Set the amounts in initial organic matter
C of biomass(=BORGM) and humus (=HORGM)
C
        DO 200 IMON=1,12
	    THISPI(IMON)=PI_CEQ_MON(IMON,IL)/1000
200     CONTINUE
        CALL GETEQOM_FROM_ROTHC(TOC(IL),IOM(IL),
     &             BORGM,HORGM,RORGM,DORGM,
     &             CLARRAY(NSOIL,IL),DEPTH,
     &             PHARRAY(NSOIL,IL),PHP1ARRAY(NSOIL,IL),
     &             PHP2ARRAY(NSOIL,IL),
     &             EQMODEL,THISPI,ICFACTOR(IL),LTA_AWC,LTA_TEMP,
     &		     ITFUNC,IMFUNC,WMAX,WSAT,DPM_RPM_RATIO)
C
C Partition the biomass, humus, DPM and RPM in the SOM layer 
C between the equal soil layers
C
        BCARB0(IL)=BORGM
        HCARB0(IL)=HORGM
        DPMCARB0(IL)=DORGM											
        RPMCARB0(IL)=RORGM											
        DPMNIT0(IL)=DPMCARB0(IL)*HZ1(IL)
        RPMNIT0(IL)=RPMCARB0(IL)*0.01
C
C Set the initial N15 in organic pools to 0
C
        HNLAB0(IL)=0
        BNLAB0(IL)=0
        DPMNLAB0(IL)=0
        RPMNLAB0(IL)=0
C
C Set the amount of nitrogen in BIO and HUM pools from N:C ratio of BIO&HUM
C
        BNIT0(IL)=BCARB0(IL)*HY(NSOIL,IL)
        HNIT0(IL)=HCARB0(IL)*HY(NSOIL,IL)
100   CONTINUE
C
C Leave SETORGM
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE SETSIM(MCSTART,SOIL,DRRAT,
     &                 FPLANT,FFYM,FDEC1,FPDEC1,DECOMP,
     &                 TPAR,TPPAR,TFPAR)
C
C Sets up simulation for simple phase
C
      IMPLICIT NONE
C
C Variables local to this routine
C
      INTEGER I
C
C Variables passed to/from this routine
C
 	REAL SOIL(0:6)
      INTEGER MCSTART
	REAL FPLANT(5),FFYM(5),FDEC1(5),FPDEC1(5),DECOMP(5)
	REAL DRRAT,TPAR,TPPAR,TFPAR
C
C Set fixed parameters for current land-use
C
      DO 10 I = 1,5
        FPLANT(I) = 0.0
 10   CONTINUE
      FFYM(1) = 0.49
      FFYM(2) = 0.49
      FFYM(3) = 0.0
      FFYM(4) = 0.0
      FFYM(5) = 0.02
      FDEC1(3) = 0.46
      FDEC1(5) = 0.54
      FPDEC1(4) = 0.46
      FPDEC1(5) = 0.54 
      DECOMP(1) = 10.0
      DECOMP(2) = 0.3
      DECOMP(3) = 0.66
      DECOMP(4) = 0.66
      DECOMP(5) = 0.02
      FPLANT(1)=DRRAT/(1+DRRAT)
      FPLANT(2)=1/(1+DRRAT)
      TPAR=FDEC1(3)+FDEC1(5)
      TPPAR=FPDEC1(4)+FPDEC1(5)
      TFPAR=FPLANT(1)+FPLANT(2)
      END
C
C------------------------------------------------------------
C
      SUBROUTINE SETTIME_CN(SECONDS,CONVER_F)
C
C Subroutine to set time factors from SECONDS
C
      IMPLICIT NONE
C
C Variables passed to/from calling subroutine
C
	REAL SECONDS
      REAL CONVER_F
C
C Set factors according to timestep
C
	CONVER_F=SECONDS/(7*24*60*60)
	END
C
C------------------------------------------------------------
C
      SUBROUTINE SETFERT(FERTADD,FERTADD15)
C
C Set initial fertiliser additions
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER IS					! Local timesteps counter
	INTEGER IT					! Local fertiliser type counter
	INTEGER MAXFERTSTEPS		! Maximum number of timesteps 
								!	over which fertiliser is added
      PARAMETER (MAXFERTSTEPS=30*60*24*7*3)	! Timesteps = 3 wks (max.30 secs)
C
C Variables passed to/from calling subroutine
C
	REAL FERTADD(MAXFERTSTEPS,4)	! INOUT:Fertiliser N added in this time 
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
	REAL FERTADD15(MAXFERTSTEPS,4)	! INOUT:Fertiliser N15 added in this time
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
C
C Initialise fertiliser addition to zero
C
      DO 100 IS=1,MAXFERTSTEPS
	  DO 200 IT=1,4
          FERTADD(IS,IT)=0
          FERTADD15(IS,IT)=0
200     CONTINUE
100   CONTINUE
      END
C
C------------------------------------------------------------
C
      SUBROUTINE SOILFIX(NSOIL,CRIT,AMMN,AMMN15,SOILN,SOIL15)
C
C Subroutine to set NO3 and NH4 residual
C            : 0.5kg NO3/ 5cm layer/ha (SOILN(1...10))
C            : 5.0kg NO3/50cm layer/ha (SOILN(11...12))
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER,MAXLAYER1,MAXDEPTH
	PARAMETER (MAXLAYER=60)
	DATA MAXDEPTH /300/
	DATA MAXLAYER1 /60/
      INTEGER MAXSOIL
      PARAMETER (MAXSOIL=50)
	INTEGER IL
	REAL SPL
C
C Variables passed to/from calling subroutine
C
      INTEGER NSOIL
      REAL CRIT(MAXSOIL,MAXLAYER)
	REAL AMMN(MAXLAYER),AMMN15(MAXLAYER)
	REAL SOILN(MAXLAYER),SOIL15(MAXLAYER)
C
C ..Ammonium-N
C
	SPL=MAXDEPTH/(MAXLAYER1*50.)
      DO 100 IL=1,MAXLAYER1
        AMMN(IL)=CRIT(NSOIL,IL)*SPL
        AMMN15(IL)=0
C
C ..Nitrate-N
C
        SOILN(IL)=CRIT(NSOIL,IL)*SPL
        SOIL15(IL)=0
  100 CONTINUE
C
C Leave SOILFIX
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE STARTORG(NSOIL,TOC,TOCARRAY,IOM,IOMARRAY)
C
C Subroutine to set starting value of organic matter
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER,MAXDEPTH,IDEPTH,MAXLAYER1
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
      INTEGER MAXSOIL
      PARAMETER (MAXSOIL=50)
	INTEGER IL
C
C Variables passed to/from calling subroutine
C
      INTEGER NSOIL
	REAL TOCARRAY(MAXSOIL,MAXLAYER),IOMARRAY(MAXSOIL,MAXLAYER)
	REAL TOC(MAXLAYER),IOM(MAXLAYER)
C
C Set TOC and IOM
C
      DO 100 IL=1,MAXLAYER1
        TOC(IL)=TOCARRAY(NSOIL,IL)
	  IOM(IL)=IOMARRAY(NSOIL,IL)
	  IF(IOM(IL).LT.0.0)THEN
	    CALL GET_IOM_FROM_FALLOON_EQN(TOC(IL),IOM(IL))
	  ENDIF
100   CONTINUE
C
C Leave STARTORG
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE SUM_PLANT_INPUT(TOTPIC_LU,PI_CEQ_MON,LU1)
C
C Subroutine to set starting value of organic matter
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER IL					! Local layer counter
	INTEGER IMON				! Local month counter
C
C Variables passed to/from calling subroutine
C
	INTEGER LU1					! IN:First land use
	REAL PI_CEQ_MON(12,MAXLAYER)! IN: Equilibrium plant C input each month
								!         in each layer (kgC/ha/month/layer)
	REAL TOTPIC_LU				! IN/OUT:Plant C input calculated using 
								!	  ROTHC equilibrium run for this LU (kgC/ha/yr)
C
C Sum plant input over the whole year and the whole profile
C
      TOTPIC_LU=0
      DO 100 IMON=1,12
	  DO 200 IL=1,MAXLAYER
	    TOTPIC_LU=TOTPIC_LU+PI_CEQ_MON(IMON,IL)
200     CONTINUE
100   CONTINUE
C
C Leave SUM_PLANT_INPUT
C
      END

!----------------------------------------------------------------------------------------------

      subroutine add_pi_to_dpmrpm(maxlayer, measlay, dpm2rpm,
     &                            pi_c, pi_n, icover, dpm_c, dpm_n,
     &                            rpm_c, rpm_n)
      
      ! Adds plant inputs to the DPM & RPM pools      
      
      implicit none

      ! Scalar arguments with intent(in)
      real, intent(in) :: dpm2rpm     ! DPM:RPM ratio of plant inputs
      integer, intent(in) :: maxlayer ! Number of layers in the soil profile
      integer, intent(in) :: measlay  ! Layer to which soil is measured

      ! Array arguments with intent(in)
      real, intent(in) :: pi_c(maxlayer)        ! Plant input C to soil [kgC/ha/timestep]
      real, intent(in) :: pi_n(maxlayer)        ! Plant input N to soil [kgN/ha/timestep]

      ! Array arguments with intent(inout)
      real, intent(inout) :: dpm_c(maxlayer)    ! C in Decomposable Plant Matter pool before decomposition [kgC/ha/layer]
      real, intent(inout) :: dpm_n(maxlayer)    ! N in Decomposable Plant Matter pool before decomposition [kgN/ha/layer]
      real, intent(inout) :: rpm_c(maxlayer)    ! C in Resistant Plant Matter pool before decomposition [kgC/ha/layer]
      real, intent(inout) :: rpm_n(maxlayer)    ! N in Resistant Plant Matter pool before decomposition [kgN/ha/layer]

      ! Local scalar variables
      real :: dpm_frac    ! Proportion of the plant inputs that go to DPM [proportion] 
      integer :: i, icover        ! Soil layer loop counter
      real :: rpm_frac    ! Proportion of the plant inputs that go to RPM [proportion] 

      ! Convert DPM:RPM ratio to fractions
      dpm_frac = dpm2rpm / ( 1 + dpm2rpm)
      rpm_frac = 1 / (1 + dpm2rpm)
      icover=0
      do i=1, measlay
          ! Add the plant C & N to DPM & RPM pools
          if (pi_c(i) > 0) then
              icover=1
              dpm_c(i) = dpm_c(i) + pi_c(i) * dpm_frac
              rpm_c(i) = rpm_c(i) + pi_c(i) * rpm_frac
              dpm_n(i) = dpm_n(i) + pi_n(i) * dpm_frac
              rpm_n(i) = rpm_n(i) + pi_n(i) * rpm_frac
          endif          
      enddo  ! layers
      
      end subroutine add_pi_to_dpmrpm
