     
C---------------------------------------------------------------------------------
C---------------------------------------------------------------------------------
C 
     
     
      SUBROUTINE NITRIF_BRADBURY(SOILTEMP,RAIN,NSOIL,ITFUNC,IMFUNC,
     &                  SOILW,WMAX,WSAT,
     &                  SOILN,AMMN,SOIL15,AMMN15,CRIT,TAM,
     &                  FERTADD,FERTADD15,THISFERT,THISFERT15,
     &                  FANIT,FANIT15,VOLAT,VOLAT15,CONVER_F,
     &                  PHARRAY,PHP1ARRAY,PHP2ARRAY,NITRIFN,VOLATN)
C
C Subroutine to calculate nitrification
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
C
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
	REAL NITRIFN(MAXLAYER)		! OUT:Nitrification from this layer (kg/ha/layer/timestep)
	INTEGER NSOIL				! IN:Soil code number
	REAL VOLAT					! IN:N lost by volatilisation (kgN/ha)
	REAL VOLAT15				! IN:N15 lost by volatilisation (kgN15/ha)
	REAL VOLATN(MAXLAYER)		! OUT:Volatikisation from this layer (kg/ha/layer/timestep)
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
	REAL CRIT(MAXSOIL,MAXLAYER)	! IN/OUT: Minimum nitrate level (Sozanska)
	REAL TAM(MAXLAYER)			! IN:15N/N in the ammonium in the layer
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
      FANIT=0
      FANIT15=0
C
C For each layer in the profile to 25cm calculate nitrification
C
      DO 100 IL=1,MAXLAYER1
        CALL MODFACTS_NITRIF(SOILW(IL),WMAX(IL),WSAT(IL),WRATE,
     &                SOILTEMP(IL),TRATE,
     &			    PHARRAY(NSOIL,IL),PHP1ARRAY(NSOIL,IL),
     &                PHP2ARRAY(NSOIL,IL),PHRATE,
     &                0,CRRATE,ITFUNC,IMFUNC)
C
C Calculate ANIT = Week's nitrification
C
        CRITL=CRIT(NSOIL,IL)*MAXDEPTH/(50.*MAXLAYER1)
        AMM1=AMMN(IL)*EXP(-0.6*TRATE*WRATE*PHRATE*CONVER_F)
        ANIT=(AMMN(IL)-AMM1)
        IF(ANIT.GT.(AMMN(IL)-CRITL))THEN
          ANIT=AMMN(IL)-CRITL
        END IF
        NITRIFN(IL)=ANIT
	  VOLATN(IL)=0
C
C Take nitrified N from the NH4-N pool
C
        ANIT15=ANIT*TAM(IL)
        AMMN(IL)=AMMN(IL)-ANIT
        AMMN15(IL)=AMMN15(IL)-ANIT15
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
C Set the proportion of the added fertilizer that is urea and the fraction
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
C Calculate the total amount of added fertilizer lost by volatalization
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
         AMM2=AFERT*EXP(-0.6*WRATE*TRATE*PHRATE*CONVER_F)
         FANIT=(AFERT-AMM2)
	   NITRIFN(1)=NITRIFN(1)+FANIT
         AFERT=AFERT-FANIT
         AMMN(1)=AMMN(1)+AFERT
         SOILN(1)=SOILN(1)+FANIT
         TEMPN1=TEMPN1+FANIT
         IF(FERTADD15(1,2)+FERTADD15(1,3)+FERTADD15(1,4).GT.0)THEN
          FANIT15=FANIT
          AMMN15(1)=AMMN15(1)+AFERT
          SOIL15(1)=SOIL15(1)+FANIT15
          TEMPN15=TEMPN15+FANIT15
         END IF
C
C Leave NITRIF
C
      RETURN
      END
C

C---------------------------------------------------------------------------------
C---------------------------------------------------------------------------------
C
      SUBROUTINE MINER1_BRADBURY_OLD(NSOIL,ICOVER,ITFUNC,IMFUNC,
     &                    PI_C,PI_N,PI_N15,DRRAT,								
     &                    SOILW,WMAX,WSAT,SOILTEMP,
     &                    ALPHA,BETA,GAMMA,DELTA,HY,					
     &                    BRATE,HRATE,DPMRATE,RPMRATE,				
     &                    AMMN,AMMN15,SOILN,SOIL15,CRIT,
     &                    CO2,
     &                    DNIT,DNIT15,TDNIT,T15,
     &                    BCARB0,BNIT0,BNLAB,
     &                    HCARB0,HNIT0,HNLAB,
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
	REAL*8 EB						! Proportion of biomass decomposed 
	REAL*8 EDPM					! Proportion of DPM decomposed 
	REAL*8 EH						! Proportion of humus decomposed 
	REAL*8 ERPM					! Proportion of RPM decomposed 
	REAL HCARB(MAXLAYER)		! C in soil humus at end (kgC/ha/layer)
	REAL HNIT(MAXLAYER)			! N in soil humus at end (kgN/ha/layer)
      REAL HNLAB(MAXLAYER)		! N15 in soil biomass at end (kgN15/ha/layer)
	INTEGER IDEPTH				! Depth of this layer
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
	REAL CRIT(MAXSOIL,MAXLAYER)		! IN/OUT: Minimum nitrate level (Sozanska)
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
      DO 24 M=1,MAXLAYER1
	  IF(PI_C(M).LT.0)PI_C(M)=0
	  IF(PI_N(M).LT.0)PI_N(M)=0
	  IF(PI_N15(M).LT.0)PI_N15(M)=0
        CARB(M)=PI_C(M)
        RADD(M)=PI_N(M)
        RADD15(M)=PI_N15(M)
   24 CONTINUE
C
C For each layer...
C
	WRATEDM=0
	TRATEM=0

      DO 25 M=1,MAXLAYER1
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
              EDPM=1
 	        ERPM=1
 	        EB=1
 	        EH=1
c	      ENDIF
		  LIMITED=1			
            GOTO 100														
	    ELSEIF(LIMITED.EQ.1)THEN
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
C---------------------------------------------------------------------------------
C---------------------------------------------------------------------------------

      SUBROUTINE BYPASSFLOW_old(ILAB,SECONDS,IK,FERTADD,BYRAIN,THISFERT,
     &                      THISFERT15,NFERT,IFERT,FERT,TFERT,
     &                      SOILN,WLEACH,SOIL15,WLEACH15,RAIN,
     &                      TIN,TAM,AMMN,AMMN15,CONVER_F)
C
C Subroutine to calculate fertilizer additions into soil and leaching loss via bypass flow.
C
      IMPLICIT NONE
C
C Variable local to this subroutine
C
      INTEGER MAXFERT
      PARAMETER (MAXFERT=5)
	INTEGER MAXLAYER
	PARAMETER (MAXLAYER=60)
      REAL RBY,RCRIT,BYFRAC,TS_ADD(MAXFERT),TS_ADD15(MAXFERT)
      PARAMETER (RBY=0.015, RCRIT=15, BYFRAC=0.5)
	INTEGER NF,N_TIMESTEPS
	REAL BYPASS(MAXFERT),FERTMIX(MAXFERT),FMIX
C
C Variables passed to/from calling subroutine
C
	INTEGER MAXFERTSTEPS		! Maximum number of timesteps 
								!	over which fertiliser is added
      PARAMETER (MAXFERTSTEPS=30*60*24*7*3)	! Timesteps = 3 wks (max.30 secs)
	REAL FERTADD(MAXFERTSTEPS,4)	! INOUT:Fertiliser N added in this time 
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
	REAL FERTADD15(MAXFERTSTEPS,4)	! INOUT:Fertiliser N15 added in this time
									!   of type 1=nitrate,2=amm.sulphate,
									!           3=other amm.salt,4=Urea
      REAL BYRAIN,THISFERT,THISFERT15,RAIN,CONVER_F
	REAL SOILN(MAXLAYER),AMMN(MAXLAYER),TIN(MAXLAYER),TAM(MAXLAYER)
	REAL WLEACH,WLEACH15,SOIL15(MAXLAYER),AMMN15(MAXLAYER)
	REAL TFERT(MAXFERT,3),FERT(MAXFERT)
	REAL SECONDS
	INTEGER ILAB(MAXFERT)
	INTEGER IK
	INTEGER IMIX(MAXFERT),IBY(MAXFERT),NFERT,IFERT(MAXFERT)
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
C
C For each fertilizer application...
C
      DO 82 NF=1,NFERT
        TS_ADD(NF)=0.0
        TS_ADD15(NF)=0.0
C
C in week of fertilizer application and 2 weeks following, add fertilizer to the top of the soil (not yet incorporated)
C
        N_TIMESTEPS=INT(3*(7*24*60*60)/SECONDS)-1
        IF(IK.GE.IFERT(NF).AND.IK.LE.IFERT(NF)+N_TIMESTEPS)THEN
          IF(IK.EQ.IFERT(NF))THEN
            FERTADD(NF,1)=FERT(NF)*(TFERT(NF,1)/100)
            IMIX(NF)=0
            IBY(NF)=0
          ENDIF
C
C If bypass flow has not yet occurred...
C
          IF(IBY(NF).EQ.0)THEN
C
C and if rainfall has exceeded the critical level...
C

            IF(RAIN.GT.(RCRIT*CONVER_F))THEN
              IBY(NF)=1
C
C Calculate rainfall lost by bypass flow (= fraction (BYFRAC) of rainfall above critical minimum)
C
              IF(BYRAIN.LE.0.0001)BYRAIN=BYFRAC*(RAIN-(RCRIT*CONVER_F))
C
C Calculate the N lost by bypass flow (= fraction RBY of added fertilizer * bypass flow water)
C
              BYPASS(NF)=MIN(RBY*BYRAIN,1.0)*FERTADD(NF,1)
C
C Take off bypassed fertilizer and add the rest
C
              FERTMIX(NF)=FERTADD(NF,1)-BYPASS(NF)
              FERTADD(NF,1)=0.0
              SOILN(1)=SOILN(1)+FERTMIX(NF)
C
C Save results
C
              TS_ADD(NF)=FERTMIX(NF)+BYPASS(NF)
              WLEACH=WLEACH+BYPASS(NF)
C
C Calculate labelled N lost by bypassflow
C
              IF(ILAB(NF).EQ.1)THEN
                SOIL15(1)=SOIL15(1)+FERTMIX(NF)
C
C Save results
C
                TS_ADD15(NF)=FERTMIX(NF)+BYPASS(NF)
                WLEACH15=WLEACH15+BYPASS(NF)
              ENDIF
C
C Otherwise, if rainfall has not exceeded critical level...
C
            ELSE
C
C ...set fertilizer lost by bypassflow to 0
C
              IMIX(NF)=IMIX(NF)+1
              FMIX=0
C
C Add only 1/3 fertilizer to labelled pool
C
              IF(IMIX(NF).LT.N_TIMESTEPS+2)  
     &                   FMIX=1/FLOAT(N_TIMESTEPS+2-IMIX(NF))
              FERTMIX(NF)=FMIX*FERTADD(NF,1)
              FERTADD(NF,1)=FERTADD(NF,1)-FERTMIX(NF)
              SOILN(1)=SOILN(1)+FERTMIX(NF)
C
C Save results
C
              TS_ADD(NF)=FERTMIX(NF)
C
C If fertilizer labelled, add 1/3 fertilizer to labelled pool
C
              IF(ILAB(NF).EQ.1)THEN
                SOIL15(1)=SOIL15(1)+FERTMIX(NF)
C
C Save results
C
                TS_ADD15(NF)=FERTMIX(NF)
              ENDIF
            ENDIF
          ENDIF
C
C Reset the proportion of N:N15
C
        CALL PROPN15(SOILN,SOIL15,TIN,AMMN,AMMN15,TAM)
C
C Save results
C
        THISFERT=THISFERT+TS_ADD(NF)
        IF(ILAB(NF).EQ.1)THISFERT15=THISFERT15+TS_ADD15(NF)
      ENDIF
   82 CONTINUE
      RAIN=RAIN-BYRAIN
C
C Leave BYPASSFLOW
C
      END
C
C---------------------------------------------------------------------------------
C---------------------------------------------------------------------------------

      SUBROUTINE DENITRIF3_BRADBURY(DFACT,DENIT,DN15,RAIN,CH4,CO2,
     &                     WMAX,SOILW,SOILN,CRIT,NSOIL,TIN,SOIL15,
     &                     DENITRIFN)
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
 	REAL DEFI(MAXLAYER)			! Calculate the amount of water before layer reaches field capacity
	REAL XRAIN					! Rainfall this timestep (mm)
	REAL DPROP					! Proportion (DPROP) of max.denitrification 
								! depending on how far below field capacity
	REAL TRAT					! N15:N in nitrate pool
	REAL DL15					! N15 lost by denitrification from this layer (kg N / ha)
	REAL DLAYER					! N lost by denitrification from this layer (kg N / ha)
	INTEGER IL					! Layer counter
	REAL SCRIT					! Minimum amount of nitrate in the layer (kg N / ha)
C
C Variables passed to/from calling subroutine
C
	REAL CO2(MAXLAYER)			! IN:CO2 emitted (kgC/ha/layer)
	REAL CH4(MAXLAYER)			! IN:CH4 emitted (kgC/ha/layer)
	REAL DFACT(MAXLAYER)		! IN: Denitrification factor
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL TIN(MAXLAYER)			! IN:15N/N in the nitrate in the layer
	REAL CRIT(MAXSOIL,MAXLAYER)		! IN/OUT: Minimum nitrate level (Sozanska)
	REAL DENIT					! IN:N lost by denitrification (kgN/ha)
	REAL DENITRIFN(MAXLAYER)	! OUT: N lost by denitrification (kgN/ha/layer
	REAL DN15					! IN:N15 lost by denitrification (kgN15/ha)
	REAL RAIN					! IN: Rainfall (mm/timestep)
	INTEGER NSOIL				! IN:Soil code number
C
C Initialize Variables
C DENIT=Total N denitrified from top 50cm this week
C DN15=Total N15 denitrified from top 50cm this week
C XRAIN=rainfall this week
C
      DENIT=0
      DN15=0
      XRAIN=RAIN
C
C For each layer down the profile...
C
      DO 20 IL=1,MAXLAYER1
C
C Calculate the amount of water before layer reaches field capacity = DEFI()
C
        IF(SOILW(IL).LT.WMAX(IL))THEN
          DEFI(IL)=WMAX(IL)-SOILW(IL)
        ELSE
          DEFI(IL)=0.0
        END IF
C
C Calculate the denitrification in this layer = DLAYER...
C ...if rainfall exceeds capacity of the layer allow full denitrifcation
C as indicated by the CO2 (ie biological activity) in the layer
C
        IF(XRAIN.GT.DEFI(IL))THEN
          DLAYER=(CO2(IL)+CH4(IL))*DFACT(IL)*SOILN(IL)
C
C ...if rainfall is less than capacity of the layer allow only a
C proportion (DPROP) of denitrification depending on how far below
C the capacity of the layer
C
        ELSE
	    IF(SOILW(IL).GE.WMAX(IL))THEN
	      DPROP=1
	    ELSEIF(SOILW(IL).LT.WMAX(IL))THEN
            IF(WMAX(IL).GT.0)THEN
		    DPROP=(SOILW(IL)+XRAIN)/WMAX(IL)
	      ELSE
	        DPROP=0
	      ENDIF
	    ENDIF
          DLAYER=(CO2(IL)+CH4(IL))*DFACT(IL)*SOILN(IL)*DPROP
        END IF
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

C---------------------------------------------------------------------------------
C---------------------------------------------------------------------------------
C
      SUBROUTINE GET_CTON_FOEREID(DPMCTON,RPMCTON,LUCODE)
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
	PARAMETER (MAXLU=9)	

      INTEGER IL					! Layer counter
      INTEGER ILU					! Land use counter
	REAL DCTON(MAXLU)			! C:N ratio of DPM
C                                    Ara  Gra  For  Nat  Mis  SRC
	DATA (DCTON(ILU),ILU=1,MAXLU) /14.3,14.3,49.0,30.0,14.3,49.0,
     &                               14.3,14.3,49.0/ 
	REAL RCTON(MAXLU)			! C:N ratio of RPM
C                                    Ara  Gra  For  Nat  Mis  SRC
	DATA (RCTON(ILU),ILU=1,MAXLU) /41.9,41.9,49.0,30.0,14.3,49.0,
     &                               41.9,41.9,41.9/ 
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
C Leave GET_CTON_FOEREID
C
      RETURN
      END