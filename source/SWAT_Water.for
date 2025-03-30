C------------------------------------------------------------
C
C Rothamsted Carbon and Nitrogen Turnover Model
C Soil Water Routines
C
C Adapted from SWAT
C (Reference: Neitsch, S.L., Arnold, J.G., Kiniry, J.R., Williams, J.R., 2002. 
C Soil Water Assessment Tool - User Manual - Version 2000. TWRI Report. 
C Texas Water Resources Institute, College Station, Texas, Temple, Texas, p. 506.)
C
C by Pia Gottschalk, Jessica Ballerby & Jo Smith
C University of Aberdeen
C 31/08/07
C
C crop PARAMETERS hardcoded: RL = minimum effective stomatal resistance of a single leaf (s m-1)
C-------------------------------------------------------------
C
C EXTERNAL SUBROUTINES
C 1. INTERCEPTION_SWAT_WATER
C 2. EVAP_SWAT_WATER
C
C Modified external subroutines for SWAT
C
C 3. GETWEATHER_AT_FC_SWAT
C 4. GETWEATHER_SWAT
C 5. INIT_SUNDIAL_WATER_SWAT
C
C-------------------------------------------------------------
C
	SUBROUTINE EVAP_SWAT_WATER(
C INPUTS
     &						T,Z,AVP,LATD,JD,CURRENTPLHEI,
     &						CURRENTLAI,CANSTOR,SOILW,FIELDCAP,
     &						WMAX,WILTPOINT,ROOTLAYER,
     &						CURRENTPLBIOM,
C OUTPUTS
     &						EVAPW,
C INPUT/OUTPUTS
     &						WATCONT)
C
C Subroutine to calculate actual evapotranspiration (SWAT Chapter 7.3 actual evapotranspiration)
C
C For well-watered plants under neutral atmospheric stability and assuming
C logarithmic wind profiles, the Penman-Monteith equation may be written (Jensen
C et al., 1990):
C
C L*Et = (DELTAVP*(RN-G)+PSC*K1*(0.622*L*pair/P)*(SVP - AVP)/RA)/(DELTAVP+PSC*(1+RC/RA))
C
C where L is the latent heat of vaporization (MJ kg-1), 
C Et is the maximum transpiration rate (mm d-1), 
C DELTAVP is the slope of the saturation vapor pressure-temperature curve, de/dT (kPa DegreeC-1)
C RN is the net radiation (MJ m-2 d-1),
C G is the heat flux density to the ground (MJ m-2 d-1)
C PSC is the psychrometric constant (kPa DegreeC-1)
C K1 is a dimension coefficient needed to ensure the two terms in the numerator have the same units (for uz in m s-1, K1 = 8.64 x 104), 
C pair is the air density (kg m-3),
C and P is the atmospheric pressure (kPa)
C SVP is the saturation vapor pressure of air at height Z (kPa),
C AVP is the water vapor pressure of air at height Z (kPa),
C RA is the diffusion resistance of the air layer (aerodynamic resistance) (s m-1).
C RC is the plant canopy resistance (s m-1),
C
	IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	PARAMETER (MAXDEPTH=300)
	REAL L						! latent heat of vaporization (MJ kg-1)
	REAL PSC					! sychrometric constant (kPa DegreeC-1)
	REAL SVP					! saturation vapor pressure of air at height Z (kPa)
	REAL ALBEDO					! albedo
	REAL WS						! wind speed at 2 m height [m s-1]
	REAL DELTAVP				! slope of the saturation vapor pressure-
								! temperature curve, de/dT (kPa DegreeC-1)
	REAL G						! heat flux density to the ground (MJ m-2 d-1)
	REAL ETPOT					! Penman-Monteith potential evapotranspiration [mm day-1]
	REAL ET						! Penman-Monteith potential evapotranspiration for actual crop [mm day-1]
	REAL ET_DEMAND				! evaporative demand after intercepted water is evaporated [mm day-1]
	REAL ACTPLANTTRANSP			! actual plant transpiration (mm H2O)
	REAL ES						! total soil evaporation (mm H2O)
	REAL ACTPLANTTRANSP50		! actual plant transpiration from 50 cm soil profile (mm H2O)
	REAL EVAPW					! total water evapotranspirated from 50 cm soil profile (mm H2O)
C
C Variables passed to/from calling subroutine
C
	REAL LATD					! IN: latitude in degree
	REAL CANSTOR				! IN: amount of free water held in the canopy on a given day (mm H2O)
	REAL CURRENTLAI				! IN: leaf area index for a given day
	INTEGER ROOTLAYER			! IN: Layer until which roots have grown
	REAL T						! IN: mean daily air temperature (at 2 m height) [�C]
	REAL AVP					! IN: water vapor pressure of air at height Z (kPa)
	INTEGER JD					! IN: Julian day
	REAL Z						! IN: elevation above sea level [m]
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL WILTPOINT(MAXLAYER)	! IN: Total water conent of the soil
								!	  at WP (mm / layer)
	REAL WMAX(MAXLAYER)			! IN: Available water at field cap. (mm/layer)
	REAL FIELDCAP(MAXLAYER)		! IN: soil water content at FC (mm/layer)
	REAL CURRENTPLHEI			! IN: canopy hight in cm
	REAL CURRENTPLBIOM			! IN: Current plant biomass [ kg ha-1]
      REAL WATCONT(MAXLAYER)		! IN/OUT: Total water content of the soil
								!		layer (mm / layer)
C
C Get WATCONT from SOILW and WILTPOINT
C
	CALL GET_WATCONT_FROM_SOILW(SOILW,WILTPOINT,WATCONT)
C
C Calculation of climate variables valid for PET (alfalfa) & ET (actual crop)
C
C ... L latent heat of vaporization (MJ kg-1)
C
	CALL LATHEAT(T,L)
C
C ... psychrometric constant [kPa �C-1]
C
	CALL PSYCHO_D(Z,PSC)
C
C ... SVP - AVP saturation vapour pressure deficit [kPa]
C
	CALL SATVAPOURPRESSURE_D(AVP,T,SVP)
C
C ... slope vapour pressure curve [kPa �C-1]
C
	CALL SLOPEVAPOURPRESSURE_D(T,SVP,DELTAVP)
C
C ... soil heat flux density [MJ m-2 day-1]
C
	G=0
C
C ... wind speed at 2 m height [m s-1]
C
	WS=2
C
C Calculation of potential evapotranspiration: reference crop: alfalfa  at 40 cm height (ETPOT)
C
	CALL PET_ALFALFA(LATD,JD,T,AVP,SVP,L,WS,DELTAVP,G,PSC,ETPOT)
C
C Calculation of potential evapotranspiration: actual crop (ET)
C
	IF (CURRENTPLHEI.GT.0.0) THEN
		CALL PET(LATD,JD,T,AVP,SVP,L,WS,DELTAVP,G,PSC,ALBEDO,
     &			CURRENTPLHEI,CURRENTLAI,CURRENTPLBIOM,ETPOT,ET)
	ELSE
		ET=0
	END IF
C
C For running bare sites, set ET=0 throughout
C	ET=0
C
C Calculation of evaporative demand after evaporation of intercepted water (ET_DEMAND)
C
	CALL EVAP_INTERCEPT(ETPOT,CANSTOR,ET_DEMAND)
C
C For running bare sites, set ET_DEMAND=ETPOT throughout
C	ET_DEMAND=ETPOT
C
C Calculation of actual plant transpiration (actET)
C 
	IF (CURRENTPLHEI.GT.0.0) THEN
		CALL PLANT_TRANSPIRATION(WMAX,SOILW,CURRENTLAI,ET,ROOTLAYER,
     & ACTPLANTTRANSP,ACTPLANTTRANSP50,WATCONT)
	ELSE
		ACTPLANTTRANSP=0.0
	END IF
C
C For running bare sites, set ACTPLANTTRANSP=0 & CURRENTPLBIOM=0throughout
C	ACTPLANTTRANSP=0
C	CURRENTPLBIOM=0
C
C Calculation of soil evaporation
C
	CALL SOIL_EVAP(WATCONT,FIELDCAP,WILTPOINT,CURRENTPLBIOM,
     &			SOILW,ET_DEMAND,ACTPLANTTRANSP,ES)
C
	EVAPW=CANSTOR+ACTPLANTTRANSP50+ES
	WRITE(32,10)JD,ETPOT-ET_DEMAND,ETPOT,ET,ACTPLANTTRANSP,ES
10	FORMAT(I6,5F6.2)
C
C Leave EVAP_SWAT_WATER
C
	RETURN
	END
C
C------------------------------------------------------------
C
	SUBROUTINE INTERCEPTION_SWAT_WATER(
C INPUTS
     &						LAI,CURRENTLAI,CANMAX,
C INPUTS/OUTPUTS
     &						RAIN,
C OUTPUTS
     &						CANSTOR)
C
C Calculation of how much water is intercepted by the plants (Chapter 7.1, pg. 118)
C
C Variables local to this subroutine
C
	REAL CANDAYMAX		!maximum amount of water that can be trapped in the 
						!canopy on a given day (mm H2O)
C
C Variables passed to/from calling subroutine
C
	REAL LAI			! IN: maximum leaf area index for the plant
	REAL CURRENTLAI		! IN: leaf area index for a given day
	REAL CANMAX			! IN: maximum amount of water that can be trapped in 
						! the canopy when the canopy is fully developed (mm H2O),
	INTEGER I
	REAL RAIN			! IN/OUT: precipitation on a given day before and 
						! after canopy interception is removed (mm H2O)
	REAL CANSTOR		! OUT: amount of free water held in the canopy on a 
						! given day (mm H2O)
C
C Calculation of CANDAYMAX
C
	CANDAYMAX = CANMAX * CURRENTLAI/LAI					! equation 7.1.1 (SWAT 
														! Theoretical Doc.)
C
C Calculation of how much water is held in the canopy and how much precipitation 
C reaches the ground ( -> precipitation is reduced by what is held in the canopy storage)
C
	if (RAIN .LE. CANDAYMAX-CANSTOR) then 
		CANSTOR = CANSTOR + RAIN						! equation 7.1.2 (SWAT 
														! Theoretical Doc.)
		RAIN = 0.
	else
		RAIN = RAIN - (CANDAYMAX - CANSTOR)	
		CANSTOR = CANDAYMAX								! equation 7.1.3 (SWAT 
														! Theoretical Doc.)
	end if
C
C Leave INTERCEPTION_SWAT_WATER
C
	RETURN
	END
C
C-------------------------------------------------------------
C INTERNAL SUBROUTINES
C-------------------------------------------------------------
C
	SUBROUTINE AERODYRES(WS,CURRENTPLHEI,RA)

C Subroutine to calculate the aerodynamic resistance to sensible heat and vapor transfer
C SWAT Theoretical Doc., Chapter 7.2.1.2, p. 121
C
C RA = [ln[(ZW-d)/ZOM]/ln[(ZP-d)/ZOV]]/K**2*WS
C ZW is the height of the wind speed measurement (cm)
C ZP is the height of the humidity (psychrometer) and temperature measurements (cm)
C D is the zero plane displacement of the wind profile (cm)
C ZOM is the roughness length for momentum transfer (cm)
C ZOV is the roughness length for vapor transfer (cm)
C K is the von K�rm�n constant
C and WS is the wind speed at height ZW (m s-1)
C CURRENTPLHEI is canopy hight (cm)
C
	IMPLICIT NONE
C
C Variables local to this subroutine
C
	REAL ZW						! height of the wind speed measurement (cm)
	REAL ZP						! height of the humidity (psychrometer) and temperature measurements (cm)
	REAL D						! zero plane displacement of the wind profile (cm)
	REAL ZOM					! roughness length for momentum transfer (cm)
	REAL ZOV					! roughness length for vapor transfer (cm)
	REAL K						! K�rm�n constant
C
C Variables passed to/from calling subroutine
C
	REAL WS						! IN: wind speedd (m s-1)
	REAL CURRENTPLHEI			! IN: canopy hight in cm
	REAL RA						! OUT: diffusion resistance of the air layer (aerodynamic resistance) (s m-1)
C
C The von K�rm�n constant is considered to be a universal constant
C in turbulent flow. Its value has been calculated to be near 0.4 with a range
C of 0.36 to 0.43 (Jensen et al., 1990). A value of 0.41 is used by SWAT for
C the von K�rm�n constant.
C
	K=0.41
C
C Calculation of ZOM:
C
	IF (CURRENTPLHEI .LE. 200) THEN
		ZOM=0.123*CURRENTPLHEI				! equation 7.2.4 (SWAT 
											! Theoretical Doc.)
	ELSE
		ZOM=0.058*CURRENTPLHEI**1.19		! equation 7.2.5 (SWAT 
											! Theoretical Doc.)
	END IF
C
C Calculation of ZOV
C
	ZOV=0.1*ZOM								! equation 7.2.6 (SWAT 
											! Theoretical Doc.)
C
C Calculation of D
C
	D= 2/3*CURRENTPLHEI						! equation 7.2.7 (SWAT 
											! Theoretical Doc.)
C
C Assumptions of ZW and ZP:
C
	ZW = 170
	ZP = 170
C
C Calculation of RA
C
	RA = (log((ZW-D)/ZOM)/log((ZP-D)/ZOV))/K**2*WS	! equation 7.2.3 (SWAT 
													! Theoretical Doc.)
C
C Leave AERODYRES
C
	RETURN
	END 
C
C-------------------------------------------------------------
C
	SUBROUTINE CANOPYRES(LAI,RC)
C
C Subroutine to calculate canopy resistance (SWAT Theoretical Doc., Chapter 7.2.1.3, p. 122)
C
C RC= RL/(0.5* LAI)
C RL is is the minimum effective stomatal resistance of a single leaf (s m-1)
C RL = 250 for maize, according to Mihailovic & Eitzinger (2007), Ecological modelling
C
	IMPLICIT NONE
C
C Variables local to this subroutine
C
	REAL RL					! minimum effective stomatal resistance of a single leaf (s m-1)
C
C Variables passed to/from calling subroutine
C
	REAL LAI				! IN: maximum leaf area index for the plant
	REAL RC					! plant canopy resistance (s m-1)
C
	RL = 250
	RC = RL/(0.5*LAI)		! equation 7.2.8 (SWAT 
							! Theoretical Doc.)
C
C Leave CANOPYRES
C
	RETURN
	END
C
C-------------------------------------------------------------
C
	SUBROUTINE EVAP_INTERCEPT(ETPOT,CANSTOR,ET_DEMAND)
C
C Subroutine to calculate evaporative demand after intercepted water is evaporated
C (SWAT Theoretical Doc., Chapter 7.3.1, p. 128)
C
	IMPLICIT NONE
C
C Variables passed to/from calling subroutine
	REAL ETPOT			! IN: Penman-Monteith potential evapotranspiration [mm day-1]
	REAL CANSTOR		! IN: amount of free water held in the canopy on a 
						! given day (mm H2O)
	REAL ET_DEMAND		! OUT: evaporative demand after intercepted water is evaporated [mm day-1]
C
C Calculation of evapotranspirative demand after all the water is evaporated from the canopy
C
	IF (ETPOT .LE. CANSTOR) THEN
		CANSTOR = CANSTOR - ETPOT
		ET_DEMAND=0
	ELSE
		ET_DEMAND=ETPOT-CANSTOR
		CANSTOR=0
	END IF
C
C Leave EVAP_INTERCEPT
C
	RETURN
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE DEF_SWAT(SUM_TS,SOILW,WMAX,IDATEFC,IANTHES,
     &               K_TS,RAIN,EVAP,AIRTEMP,RDD,TMMN,TMMX,VP,WN,
     &               MODTYPE,WSAT,FLOWPROP,IYEAR,ISTHARV, ISOWNJ,
     &				N_STEPS,LAI,CURRENTLAI,CANMAX,CANSTOR,CONVER_F,
     &                PLANTBIOMASS,PLANTHEIGHT,LCROP,C1,T1,
     &				CURRENTPLBIOM,CURRENTPLHEI,DDAYS,
     &				Z,AVP,LAT,FIELDCAP,WATCONT,SATWATCONT,
     &				WILTPOINT,MCROP,IROCKS,RRG)
C
C Subroutine to calculate soil moisture deficit in mid-July.
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXCROP				! Max.no.of crop types allows
      PARAMETER (MAXCROP=36)
	INTEGER MAXGROW			! Max.no.of growing seasons allowed
	PARAMETER (MAXGROW=300)
	INTEGER MAXSOIL			! Max.no.of soil types
      PARAMETER(MAXSOIL=50)
	INTEGER MAXLAYER		! No.of layers in the soil profile
	INTEGER MAXLAYER1		! No.of layers in the soil profile
	INTEGER MAXDEPTH		! Maximum depth of the soil profile
	INTEGER IDEPTH			! Depth of this layer
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
	REAL REFIL(MAXLAYER)	! Water deficit in the layer mm/layer
	REAL DRAIN(MAXLAYER)	! Water drained from this layer (mm/layer)
      REAL AVAILW				! Water available in this layer (mm/layer)
	INTEGER IL				! Layer counter
	INTEGER NN				! Date counter
	INTEGER IBARE			! Soil covered / frozen (1=no,0=yes)
	INTEGER IS_TS			! days since sowing (for root length calculation
	INTEGER M				! layer until which roots have grown
	REAL ROOTS				! root length (cm)
C
C Variables passed to/from calling subroutine
C
      INTEGER SUM_TS
	INTEGER IDATEFC			! IN:Date of field capacity 
	INTEGER IANTHES			! IN:Date of anthesis
	INTEGER K_TS			! IN:Date counter
	INTEGER MODTYPE			! Type of crop model (1=MAGEC, 2=SUNDIAL)
	INTEGER IYEAR			! IN:Current growing season number
	INTEGER ISTHARV			! IN:Timesteps from 01/01/01 to first harvest
	INTEGER ISOWNJ(0:MAXGROW) ! IN:Timesteps from 01/01/01 to sowing date 
	INTEGER N_STEPS			! IN:No.timesteps in a year
	REAL SOILW(MAXLAYER)	! IN:Available water (mm/layer)
      REAL FIELDCAP(MAXLAYER)		! IN/OUT: Total water in the layer at 
	                            !		field capacity (mm / layer)
	REAL SATWATCONT(MAXLAYER)	! IN: Total water content at 
								!		saturation (mm / layer)
      REAL WATCONT(MAXLAYER)		! IN: Total water content of the soil
								!		layer (mm / layer)
	REAL WILTPOINT(MAXLAYER)	! IN: Total water conent of the soil
								!		at WP layer (mm / layer)
	REAL WMAX(MAXLAYER)		! IN:Available water at field cap. (mm/layer)
	REAL RAIN				! IN: Rainfall (mm/timestep)
	REAL EVAP				! IN: Potential evap. (mm/timestep)
	REAL AIRTEMP			! IN: Air temperature (deg.C/timestep)
	REAL RDD
	REAL TMMN
	REAL TMMX
	REAL VP
	REAL WN
	REAL FLOWPROP			! IN:Proportion of flow needed to achieve
							!	     water table at depth WTABLE
	REAL WSAT(MAXLAYER)		! IN:Available water at saturation (mm/layer)
	REAL Z					! elevation above sea level [m]
	REAL AVP				! IN: water vapor pressure of air at height Z (kPa)
	REAL LAT				! IN: Latitude
	REAL LAI				! maximum leaf area index for the plant
	INTEGER ROOTLAYER		! Layer until which roots have grown
	REAL CURRENTLAI			! leaf area index for a given day
	REAL CANMAX				! canmx is the maximum amount of water that can be trapped in the canopy when the canopy is fully developed (mm H2O),
	REAL CANSTOR			! amount of free water held in the canopy on a given day (mm H2O)
	REAL DDAYS				! OUT:Degree days since sowing (deg.C)
	REAL CONVER_F
	REAL PLANTBIOMASS
	REAL PLANTHEIGHT
	REAL C1(0:MAXCROP)
	REAL T1(0:MAXCROP)
	INTEGER LCROP
	REAL CURRENTPLBIOM
	REAL CURRENTPLHEI
	REAL EVAPW
	INTEGER MCROP				! IN: Crop type code number
	INTEGER IROCKS(0:MAXCROP)	! IN: Maximum rooting depth / impermeable layer
	REAL RRG(0:MAXCROP)			! IN: Rate of root growth (cm/week)
C
C Set the variable for summing weeks
C
      SUM_TS=0
C
C Set soil water to max
C
      DO 1 IL=1,MAXLAYER1
       SOILW(IL)=WMAX(IL)
    1 CONTINUE
C
C Initialise canopy storage
	CANSTOR = 0.
C
C Calculate water deficit from soil at field capacity up to
C anthesis of previous crop
C (simulation proper started at anthesis of previous crop)
C
      DO 20 NN=IDATEFC,IANTHES-1
C
C Read in weather data
C
        IF(MODTYPE.EQ.1)THEN
          READ(3,*)K_TS,RAIN,EVAP,AIRTEMP, RDD, TMMN, TMMX, VP, WN
	  ELSE
          READ(3,*)K_TS,RAIN,EVAP,AIRTEMP
		READ(3000,*) K_TS,AVP
        ENDIF
        SUM_TS=SUM_TS+1
        IF(EVAP.LT.0)EVAP=0
C
C If freezing soil acts as if it is bare
C Soil is set to covered during "previous" growing season which is estimated to last
C from sowing date of the actual first season simulated minus total timesteps per year to
C previous harvest and soil is set to bare from the start of the simulation to the estimated
C sowing date
C
        IF(AIRTEMP.LT.0.0) THEN
		IBARE=1
	  ELSE IF ((K_TS.LE.ISTHARV) .AND.
     &	(K_TS.GT.1))THEN						! Grignon, France
!     &   (K_TS.GT.(ISOWNJ(2)-N_STEPS)))THEN				! Ochinga_M (bi-annual cropping)
!     &   (K_TS.GT.(ISOWNJ(1)-N_STEPS)))THEN				! Monte Carlo,MZ12 (annual cropping),Paulinenaue P5b
	    IBARE=0
        ELSE
          IBARE=1
        END IF
C
C
C Calculation of current crop parameters for actual evapotranspiration
C
C ... 1. calculation of degree days
C ...		IDATEFC (date where soil is at FC corresponds to DOY=JD=ik
C		DDAYS (degree days)
C ...     IS_TS (days since sowing)
C ... 1a. calculation of degree days from sowing of previous crop to IDATEFC
C					assuming that the sowing data of the previous crop is the same
C					as the first proper crop (second drop if bi-annual system) and 
C					also assuming that AIRTEMP of first day where weather is read 
C					in is valid for the growing period since sowing of previous crop
	   IF(IDATEFC==K_TS) THEN
!			DDAYS = (IDATEFC-(ISOWNJ(2)-365))*	! Ochinga_M
!     &                   AIRTEMP*7*CONVER_F		! Ochinga_M
!			IS_TS =	ABS(IDATEFC-(ISOWNJ(2)-365))		! Ochinga_M
!			DDAYS = (ABS(ISOWNJ(1)-N_STEPS))*AIRTEMP*7*CONVER_F		! MZ12
!			IS_TS = (ABS(ISOWNJ(1)-N_STEPS))						! MZ12
!			DDAYS = 61*AIRTEMP*7*CONVER_F		! Mafu
!			IS_TS = 61							! Mafu
!			DDAYS = (IDATEFC-(ISOWNJ(1)-365))* ! Monte Carlo
!     &                   AIRTEMP*7*CONVER_F		! Monte Carlo
			DDAYS =	1188						! Grignon
			IS_TS =	121							! Grignon
!			DDAYS = 0							! Paulinenaue P5b, here 0 because the previous growing seaons starts after date soil reaches FC
!			IS_TS = 0							! Paulinenaue P5b

	   ELSE
C
C ... 1b. calculation of degree days for each day from IDATEFC to ANTHES-1
C
		 IF(AIRTEMP.GT.0.0)THEN							!normally!
!		 IF(AIRTEMP.GT.0.0.AND.K_TS.GE.(ISOWNJ(1)-365))THEN			!Paulinenaue P5b
				DDAYS=DDAYS+AIRTEMP*7*CONVER_F
				IS_TS=IS_TS+(7*CONVER_F)
		 END IF
		END IF
C
C ... 2. calculation of current crop height, biomass and LAI
C
	  CALL CURRENTPLANTPARAM(PLANTBIOMASS,PLANTHEIGHT,LAI, 
     &							C1,T1,DDAYS,LCROP,
     &							CURRENTPLBIOM,CURRENTPLHEI,CURRENTLAI)
C
C Calculation of intercepted water by plants
C
	  CALL INTERCEPTION_SWAT_WATER(LAI,CURRENTLAI,CANMAX,RAIN,
     &					CANSTOR)
C
C Evapourate and drain
C
C	  CALL DRAIN_SUNDIAL_WATER(RAIN, WMAX, SOILW,DRAIN, REFIL,
C     &                            WSAT, FLOWPROP)
        CALL DRAIN_SUNDIAL_WATER2(RAIN, WMAX, SOILW, WSAT, FLOWPROP,
     &                            DRAIN, REFIL)
C
C Calculation of actual evapotranspiration (SWAT) (until anthesis) 
C (alternative to EVAP_SUNDIAL_WATER)
C
      ROOTS=IS_TS*RRG(MCROP)*CONVER_F
      CALL ROOTL_SWAT(M,MCROP,IROCKS,ROOTS)
      ROOTLAYER=M
C
	  CALL EVAP_SWAT_WATER(AIRTEMP,Z,AVP,LAT,K_TS,CURRENTPLHEI,
     &						CURRENTLAI,CANSTOR,SOILW,FIELDCAP,WMAX,
     &                        WILTPOINT,ROOTLAYER,
     &						CURRENTPLBIOM,EVAPW,WATCONT)
        CALL GET_AVAILWAT_FROM_TOTWAT(FIELDCAP,WATCONT,SATWATCONT,
     &                                WILTPOINT,WMAX,SOILW,WSAT)		 
C
C Go back and do next week
C
   20 CONTINUE
C
C Leave DEF
C
      RETURN
      END
C
C-------------------------------------------------------------
C
	SUBROUTINE GETWEATHER_AT_FC_SWAT(IFILE,IDATEFC,N_STEPS,IRYEAR,K_TS,
     &                            RAIN,EVAP,AIRTEMP,RDD,TMMN,TMMX,VP,WN,
     &                            FIXFILE,FIXFILEaVP,MODTYPE,SOILTEMP,
     &							AVP)
C
C Subroutine to move weather data file to the week when the soil is at field capacity
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      INTEGER MAXLAYER
	PARAMETER(MAXLAYER=60)
	INTEGER MAXLAYER1
	DATA MAXLAYER1 /60/
	INTEGER IL
      INTEGER MAXWEATH
	PARAMETER (MAXWEATH=300)
      INTEGER NN,IK,MEND
C
C Variables passed from calling subroutine
C
      INTEGER IFILE,IDATEFC,N_STEPS,IRYEAR,K_TS,MODTYPE
	REAL RAIN,EVAP,AIRTEMP,SOILTEMP(MAXLAYER),RDD,TMMN,TMMX,VP,WN
	REAL AVP		! IN: water vapor pressure of air at height Z (kPa)
	CHARACTER*100 FIXFILE(MAXWEATH)
	CHARACTER*100 FIXFILEaVP(MAXWEATH)	! IN: Weather files (actual vapour pressure) (kPa)
C
C Open weather file
C
      IFILE=IFILE+1
      OPEN(3,FILE=FIXFILE(IFILE),STATUS='OLD')
	OPEN(3000,FILE=FIXFILEaVP(IFILE),STATUS='OLD',ACTION='READ')
C
C Read lines to date of field capacity
C
      DO 15 NN=1,IDATEFC-1
C
C If file is at end, open next weather file
C
        IF(NN.EQ.N_STEPS)THEN
          IFILE=IFILE+1
          IK=1
          MEND=2
          CALL NEXTYEAR_SWAT(IRYEAR,IFILE,IK,MEND,FIXFILE,FIXFILEaVP)
        END IF
C
C read weather data
C
	  IF(MODTYPE.EQ.1)THEN
          READ(3,*)IK,RAIN,EVAP,AIRTEMP, RDD, TMMN, TMMX, VP, WN
	  ELSE
	    READ(3,*)IK,RAIN,EVAP,AIRTEMP
		READ(3000,*) IK, AVP
	  ENDIF
        IF(EVAP.LT.0.0)EVAP=0.0
15    CONTINUE
C
C Set soil temperature to air temperature
C
      DO 100 IL=1,MAXLAYER1
	  SOILTEMP(IL)=AIRTEMP
100   CONTINUE
C
C leave GETWEATHER_AT_FC
C
	RETURN
      END
C
C-------------------------------------------------------------
C
	SUBROUTINE GET_WATCONT_FROM_SOILW(SOILW,WILTPOINT,WATCONT)
	
	IMPLICIT NONE
C
C Variable declaration
C
C internal parameters
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	PARAMETER (MAXDEPTH=300)
	INTEGER IL					! Layer counter
C
C ... input
C
	REAL WILTPOINT(MAXLAYER)	! IN: Total water conent of the soil
								!		layer at WP (mm / layer)
	REAL SOILW(MAXLAYER)
C
C ... output
C
	REAL WATCONT(MAXLAYER)		! IN: Total water content of the soil
								!		layer (mm / layer)	
C
C Translate available soil water into total soil water
C
      DO 100 IL=1,MAXLAYER1
		IF(WATCONT(IL).GT.WILTPOINT(IL))THEN
			WATCONT(IL)=SOILW(IL)+WILTPOINT(IL)
		ELSE
			WATCONT(IL)=WATCONT(IL)+SOILW(IL)
		END IF
100   CONTINUE
C
C Leave GET_WATCONT_FROM_SOILW
C
	RETURN
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE GETWEATHER_SWAT(N_TS,N_STEPS,IFILE,IRYEAR,IK,MEND,
     &						SUM_TS,
     &                      RAIN,EVAPW,AIRTEMP,RDD,TMMN,TMMX,VP,WN,
     &                      FIXFILE,FIXFILEaVP,MODTYPE,SOILTEMP,AVP,
C required for Monte carlo - need to pass on original value ofprecipiation
     &						PRECIP)
C
C Subroutine to deliver met data
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      INTEGER MAXLAYER
	PARAMETER(MAXLAYER=60)
	INTEGER MAXLAYER1
	DATA MAXLAYER1 /60/
	INTEGER IL
      INTEGER MAXWEATH
	PARAMETER(MAXWEATH=300)
C
C Variables passed to/from calling subroutine
C
	INTEGER SUM_TS,N_TS,N_STEPS,IFILE,IRYEAR,IK,MEND,MODTYPE
      REAL RAIN,EVAPW,AIRTEMP,SOILTEMP(MAXLAYER),RDD,TMMN,TMMX,VP,WN
	REAL AVP					! IN: water vapor pressure of air at height Z (kPa)
	REAL PRECIP !original value of rain needs to be passed on from file for MCMODE
	CHARACTER*100 FIXFILE(MAXWEATH)
	CHARACTER*100 FIXFILEaVP(MAXWEATH)	! IN: Weather files (actual vapour pressure) (kPa)	
C
C In last week of year open next years weather data file
C
      IF(N_TS.EQ.N_STEPS)THEN
        IFILE=IFILE+1
        CALL NEXTYEAR_SWAT(IRYEAR,IFILE,IK,MEND,FIXFILE,FIXFILEaVP)
      END IF
C
C Read in weather data for this week
C
      SUM_TS = SUM_TS + 1
	IF(MODTYPE.EQ.1)THEN
        READ(3,*)N_TS,RAIN,EVAPW,AIRTEMP, RDD, TMMN, TMMX, VP, WN
	ELSE
	  READ(3,*)N_TS,RAIN,EVAPW,AIRTEMP
	  READ(3000,*) N_TS,AVP
	ENDIF
      IF(EVAPW.LT.0.0)EVAPW=0.0
C
C Set soil temperature to air temperature
C
      DO 100 IL=1,MAXLAYER1
	  SOILTEMP(IL)=AIRTEMP
100   CONTINUE
C
C required for Monte carlo - need to pass on original value of precipiation
C
	PRECIP=RAIN
C
C Leave Getweather
C
	RETURN
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE INIT_SUNDIAL_WATER_SWAT(WATCONT,WILTPOINT,ROOT,IROCK,
     &                              SUM_TS,IDATEFC,IANTHES,
     &                              IRYEAR,IFILE,FIELDCAP,K_TS,RAIN,
     &                              EVAP,AIRTEMP,RDD,TMMN,TMMX,VP,WN,
     &                              IAWC,NSOIL,MODTYPE,SPARMODEL,
     &                              TOC,CLAY,SAND,SILT,
     &                              BULKDENS,SATWATCONT,WTABLE,
     &                              FLOWPROP,AVERAIN,AVEPET,AVETEMP,
     &                              IYEAR,ISTHARV,ISOWNJ,N_STEPS,
C required for SWAT water
     &							  LAI,CURRENTLAI,CANMAX,CANSTOR,
     &                              CONVER_F,PLANTBIOMASS,PLANTHEIGHT,
     &							  LCROP,C1,T1,
     &						  CURRENTPLBIOM,CURRENTPLHEI,DDAYS,
     &							Z,AVP,LAT,MCROP,RRG,
     &							  IROCKS,SOILFC,WMAX,SOILW,SECONDS)
C
C Subroutine to initialise soil water 
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXCROP				! Max.no.of crop types allows
      PARAMETER (MAXCROP=36)
	INTEGER MAXGROW			! Max.no.of growing seasons allowed
	PARAMETER (MAXGROW=300)
	INTEGER MAXSOIL			! Max.no.of soil types
	PARAMETER (MAXSOIL=50)
	INTEGER MAXLAYER		! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER SPARFILE		! IN:Soil parameters read in from file
	INTEGER SPARCALC		! IN:Soil parameters calculated from TOC etc
	DATA SPARFILE,SPARCALC /1,2/
C
C Variables passed to/from calling subroutine
C ...Model descriptors
C
      INTEGER MODTYPE			! IN:Crop model type: 0=SUNDIAL, 1=MAGEC
C
C ...Timing factors
C
	INTEGER IDATEFC			! IN:Date of field capacity (1=01/01; 2=01/06)
	INTEGER IANTHES			! IN:No.timesteps from 01/01 to anthesis 
	INTEGER IRYEAR			! IN:Current weather year number
	INTEGER K_TS			! IN:No.timesteps since start of year(UNUSED?)
	INTEGER SUM_TS			! IN:Total number of timesteps passed
	INTEGER N_STEPS			! IN:No.timesteps in a year
	INTEGER IYEAR			! IN:Current growing season number
	INTEGER ISTHARV			! IN:Timesteps from 01/01/01 to first harvest
	INTEGER ISOWNJ(0:MAXGROW) ! IN:Timesteps from 01/01/01 to sowing date 
	REAL CONVER_F				! IN:Conversion between timestep & weeks
	REAL SECONDS			! Number of seconds in one timestep
C
C ...Soil factors
C
      INTEGER ISAVE			! OUT:Code to save or retrieve variables
	INTEGER IAWC			! IN:Water movement code number
	INTEGER IROCK			! IN:No.layers: 1=50cm, 2=100cm , 3=150cm
      INTEGER NUMSOIL			! IN/OUT: Number of soils defined
	INTEGER NSOIL			! IN:Soil code number
	REAL SOILW(MAXLAYER)	! IN/OUT:Available water (mm/layer)
	REAL WILTPOINT(MAXLAYER)	! IN: Total water conent of the soil
								!		layer at WP (mm / layer)
	REAL WMAX(MAXLAYER)		! IN:Available water at field cap. (mm/layer)
	REAL MAXWAT(MAXSOIL,MAXLAYER)	! IN/OUT: Maximum water content (mm/layer)
	REAL WATSAT(MAXSOIL,MAXLAYER)	! IN:Avail.water at saturation (mm/layer)
	REAL SATWATCONT(MAXLAYER)	! IN: Total water content at 
								!		saturation (mm / layer)
	REAL WSAT(MAXLAYER)		! IN:Available water at saturation (mm/layer)
	REAL FIELDCAP(MAXLAYER)	! OUT: soil water content at FC (mm/layer)
	REAL SOILFC(MAXLAYER)	! OUT: soil water content at FC (mm/layer)
	REAL WATCONT(MAXLAYER)	! IN: water content (mm/layer)
	REAL WTABLE				! IN:Water table depth in cm
	REAL WILTPOINT5(MAXSOIL,MAXLAYER)	! OUT: soil water content at wiltint point (mm/layer)
	CHARACTER*40 SNAME(MAXSOIL)		! IN/OUT: Soil name
	INTEGER SPARMODEL		! IN:Soil parameter model (from file or calc)
	REAL TOC(MAXLAYER)		! IN:Total organic C (kgC/ha/layer)
      REAL CLAY(MAXLAYER)		! IN:Clay content of the layer (%)
	REAL SAND(MAXLAYER)		! IN:Sand content of the layer (%)
	REAL SILT(MAXLAYER)		! IN:Silt content of the layer (%)
	REAL BULKDENS(MAXLAYER)	! IN:Bulk density of the layer (g/cm3)
	REAL FLOWPROP			! IN/OUT:Proportion of flow needed to achieve
							!	     water table at depth WTABLE
C
C ... Crop factors
C
	REAL ROOT				! IN:Rooting depth according to restriction (cm)
	REAL LAI			! maximum leaf area index for the plant
	REAL CURRENTLAI		! leaf area index for a given day
	REAL CANMAX			! canmx is the maximum amount of water that can be trapped in the canopy when the canopy is fully developed (mm H2O),
	REAL CANSTOR		! amount of free water held in the canopy on a given day (mm H2O)
	REAL PLANTBIOMASS		! Plant biomass at harvest [ kg ha-1]
	REAL PLANTHEIGHT		! Plant height at harvest [cm]
	INTEGER LCROP
	REAL C1(0:MAXCROP)
	REAL T1(0:MAXCROP)
	REAL CURRENTPLBIOM
	REAL CURRENTPLHEI
	REAL DDAYS				! OUT:Degree days since sowing (deg.C)
	INTEGER MCROP				! IN: Crop type code number
	INTEGER IROCKS(0:MAXCROP)	! IN: Maximum rooting depth / impermeable layer
	REAL RRG(0:MAXCROP)			! Rate of root growth (cm/week)
C
C ...Weather factors
C
	INTEGER IFILE			! IN:Current weather file
	REAL RAIN				! IN: Rainfall (mm/timestep)
	REAL EVAP				! IN: Potential evap. (mm/timestep)
	REAL AIRTEMP			! IN: Air temperature (deg.C/timestep)
	REAL Z					! elevation above sea level [m]
	REAL AVP				! IN: water vapor pressure of air at height Z (kPa)
	REAL LAT				! IN: Latitude
	REAL RDD				! IN: Weather data used by MAGEC
	REAL TMMN				! IN: Weather data used by MAGEC
	REAL TMMX				! IN: Weather data used by MAGEC
	REAL VP					! IN: Weather data used by MAGEC
	REAL WN					! IN: Weather data used by MAGEC
	REAL AVERAIN(12)		! IN/OUT: Long term average rainfall (mm)
	REAL AVETEMP(12)		! IN/OUT: Long term average temperature (deg.C)
	REAL AVEPET(12)			! IN/OUT: Long term average PET (mm)
      REAL PI_CEQ_MON(12,MAXLAYER)
C
C Read soil water parameters from input files
C
	CALL PAR_WAT_SWAT(NUMSOIL,MAXWAT,SNAME,WATSAT,WILTPOINT5)
C
C Save parameters for future use
C      
	CALL SAVE_WAT(NUMSOIL,MAXWAT,SNAME,ISAVE,WATSAT)
C
C Set the available water at field capacity and saturation
C
      CALL SWATER_SWAT(IAWC,MAXWAT,WSAT,WATSAT,WMAX,WILTPOINT,
     &				FIELDCAP,WATCONT,WILTPOINT5,SATWATCONT)
C
C Calculate the proportion restriction occuring given the minimum water table depth
C
      CALL RESTRICT_DRAINAGE(FLOWPROP,WTABLE,WSAT,WMAX,PI_CEQ_MON,
     &                       AVERAIN,AVEPET,AVETEMP,SOILW,SECONDS)
C
C Set the rooting depth
C
      CALL WATFIX(IROCK,ROOT,NSOIL)
C
C Calculate soil moisture deficit in mid-July.
C
      CALL DEF_SWAT(SUM_TS,SOILW,WMAX,IDATEFC,IANTHES,
     &               K_TS,RAIN,EVAP,AIRTEMP,RDD,TMMN,TMMX,VP,WN,
     &               MODTYPE,WSAT,FLOWPROP,IYEAR,ISTHARV, ISOWNJ,
     &				N_STEPS,LAI,CURRENTLAI,CANMAX,CANSTOR,CONVER_F,
     &                PLANTBIOMASS,PLANTHEIGHT,LCROP,C1,T1,
     &				CURRENTPLBIOM,CURRENTPLHEI,DDAYS,
     &				Z,AVP,LAT,FIELDCAP,WATCONT,SATWATCONT,
     &				WILTPOINT,MCROP,IROCKS,RRG)
C
C Leave INIT_WAT
C
	END
C
C-------------------------------------------------------------
C
	SUBROUTINE LATHEAT(T,L)
C
C Subroutine to calculate latent heat of vaporization
C
	IMPLICIT NONE
C
C Variables passed to/from calling subroutine
C
	REAL T					! IN: mean daily air temperature (at 2 m height) [�C]
	REAL L					! OUT: latent heat of vaporization (MJ kg-1)
C
	L = 2.501 - 2.362e-3 * T							! equation 3.3.6 (SWAT 
														! Theoretical Doc.)
C Leave LATHEAT
C
	RETURN
	END
C
C-------------------------------------------------------------
C
	SUBROUTINE NETRADIATION_D(JD,LATD,T,AVP,ALBEDO,RN)
C
C Subroutine to calculate solar net radiation, which is net short wave radiation minus
C net long wave radiation
C
	IMPLICIT NONE
C
C Variables local to this subroutine
C
	REAL LATT_RAD			! Latitude in Radiants
	REAL SD					! solar declination in radiants
	REAL DD					! relative distance sun to earth
	REAL CH
	REAL H
	REAL DAYL
	REAL YS
	REAL YC
	REAL MAXRAD				! maximum possible daily short wave radiation
	REAL RALB				! net short wave radiation
	REAL ALBEDO
	REAL RBO				! net emittance 2.2.20
	REAL RTO				! cloud cover adjustment factor
	REAL ROUT				! net long-wave radiation
C
C Variables passed to/from calling subroutine
C
	REAL LATD		! IN: latitude in degree
	REAL T			! IN: mean daily air temperature (at 2 m height) [�C]
	REAL AVP		! IN: water vapor pressure of air at height Z (kPa)
	INTEGER JD		! IN: Julian day
	REAL RN			! OUT: net radiation per day (net short + net long)
C
C Location, conversion of geographical degree into radiants
C
	LATT_RAD=LATD*(3.14/180)
C
C solar declination in radiants eq 2.1.2 in SWAT Theoretical Doc.
C
	SD=0.
	SD=ASIN(0.4*SIN((REAL(JD) - 82.)/58.09))  
C
C relative distance sun to earth  2.1.1 in SWAT Theoretical Doc.
C
	DD=0.
	DD=1.0+0.033*COS(REAL(JD) / 58.09)
C
C daylength calculation
C
	CH=0.
	H=0.
	CH=-SIN(LATT_RAD) * TAN(SD) / COS(LATT_RAD)
	IF (CH > 1.) THEN
		H = 0.000000000000001
	ELSE IF(CH >= -1.) THEN
		H=ACOS(CH)
	ELSE
		H=3.1416
	END IF

	DAYL= 7.6394 * H
C
C max possible daily radiation
C
	YS = 0.
	YC = 0.
	YS = SIN((LATT_RAD)) * SIN(SD)
	YC = COS((LATT_RAD)) * COS(SD)
	MAXRAD = 30. * DD * (H*YS+YC * SIN(H)) 
C
C calculate net short-wave radiation for PET all SWAT Theoretical Doc.
C
	RALB = 0.   ! net short wave radiation
C			
			!  if (avtmp(month).le.0.0) ALBEDO=0.5 !! snow
			!  if (avtmp(month).le.5.0.and.avtmp(month).gt.0.1) ALBEDO=0.3 !! bare soil
	RALB = MAXRAD * (1.0 - ALBEDO) 
C
C calculate net long-wave radiation
C net emissivity  equation 2.2.20 in SWAT manual
C AVP actual vapure pressure kPa	
	RBO = 0.
	RBO = -(0.39 - 0.158 * SQRT(AVP))
C
C cloud cover factor equation 2.2.19 from SWAT
C
	RTO = 0.
	RTO = 1.2 * (MAXRAD / MAXRAD) + (-0.2)
C
C net long-wave radiation equation 2.2.21
	ROUT = 0.
	ROUT = RBO * RTO * 4.9e-9 * ((T + 273.3)**4)
C
C calculate net radiation
C
	RN = 0.
	RN = RALB + ROUT  ! net solar radiation for PET
C
C Leave NETRADIATION_D
C
	RETURN	
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE NEXTYEAR_SWAT(IRYEAR,NFILE,IK,MEND,FIXFILE,FIXFILEaVP)
C
C Subroutine to open next years weather data file
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      INTEGER MAXWEATH
      PARAMETER(MAXWEATH=300)
C
C Variables passed from calling subroutine
C
      INTEGER IRYEAR,NFILE,IK,MEND
      CHARACTER*100 FIXFILE(MAXWEATH)
	CHARACTER*100 FIXFILEaVP(MAXWEATH)	! IN: Weather files (actual vapour pressure) (kPa)
C
C Close previously opened weather data file on channel 3
C
      CLOSE(3)
	CLOSE(3000)
C
C Open new met. data file on channel 3
C
	 IF(IK.LT.MEND)THEN
	 OPEN(3,FILE=FIXFILE(NFILE),STATUS='OLD',ERR=111)
	 OPEN(3000,FILE=FIXFILEaVP(NFILE),STATUS='OLD',ERR=111)
	 GOTO 101
111      CONTINUE
	 WRITE(15,10)FIXFILE(NFILE)
10       FORMAT('Error in weather data file!'/
     &          'Check format of file ',A12/
     &          'SIMULATION NOT COMPLETED!')
	 STOP
101      CONTINUE
C
C IRYEAR is no. of calendar years starting from the first Jan.
C
	IRYEAR=IRYEAR+1
	ENDIF
C
C Leave NEXTYEAR
C
	RETURN
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE PAR_WAT_SWAT(NUMSOIL,MAXWAT,SNAME,WATSAT,WILTPOINT5)
C
C Subroutine to reset parameters from files PARAM.OUT and PARAM.DAT
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXSOIL				! Max.no.of soil types
	PARAMETER(MAXSOIL=50)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXDEPTH			! Maximum depth of the soil profile
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
	CHARACTER*40 TEMP			! Local character string
	INTEGER I,J,N,IL,K			! Local counter variables
	REAL FNULL
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! Equilibrium TOC (kgC/ha/layer)
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! % clay content in this layer
	REAL WILTPOINT25(MAXSOIL,MAXLAYER)	! soil water content at wiltint point (mm/0-25cm)
C
C Variables passed to/from calling subroutine
C
	REAL WATSAT(MAXSOIL,MAXLAYER) ! OUT:Available water at saturation (mm/layer)
      INTEGER NUMSOIL				! OUT:Number of soils defined
	CHARACTER*40 SNAME(MAXSOIL)	! OUT:Soil name
	REAL MAXWAT(MAXSOIL,MAXLAYER)	! IN/OUT:Maximum water content (mm/layer)
	REAL WILTPOINT5(MAXSOIL,MAXLAYER)	! OUT: soil water content at wiltint point (mm/layer)
C
C Soil Parameters
C Soil name
C
      SNAME(1)='Sand'
	SNAME(2)='Loam'
	SNAME(3)='Clay'
	SNAME(4)='User'
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
C Try to open SOIL_PAR.DAT
C
C
C Read in parameters from SOIL_PAR.DAT
	REWIND(46)
303   CONTINUE
3009  READ(46,29,ERR=334,END=333)TEMP,I
      READ(46,*,ERR=334,END=333)MAXWAT(I,1),FNULL,FNULL,FNULL,
     &           FNULL,FNULL,FNULL,FNULL,FNULL,				
     &           FNULL,FNULL,TOCARRAY(I,1),FNULL,			
     &           FNULL,CLARRAY(I,1),FNULL,
     &           FNULL,FNULL,FNULL,
     &		   WATSAT(I,1),WILTPOINT25(I,1)
	MAXWAT(I,1)=MAXWAT(I,1)*(MAXDEPTH/(MAXLAYER1*1.))/25.
	WATSAT(I,1)=WATSAT(I,1)*(MAXDEPTH/(MAXLAYER1*1.))/25.
	WILTPOINT5(I,1)=WILTPOINT25(I,1)*(MAXDEPTH/(MAXLAYER1*1.))/25.
	DO 200 IL=2,MAXLAYER1
	  MAXWAT(I,IL)=MAXWAT(I,1)
	  WATSAT(I,IL)=WATSAT(I,1)
	  WILTPOINT5(I,IL)=WILTPOINT5(I,1)
200   CONTINUE
C
C Set soil name and number of soils
C
      SNAME(I)=TEMP
      NUMSOIL=I
      GOTO 303
29    FORMAT(A40/I3)
      PRINT*,'SOIL: Water parameters from SOIL_PARS.DAT'
C
C Record error in the format of SOIL_PAR.DAT
C
334   CONTINUE
C
C If can open file SOIL.DAT, read in soil parameters from this
C
	REWIND(46)
      READ(46,29,ERR=222)TEMP,I
	READ(46,*,ERR=222)MAXWAT(I,1)
	DO K=1,16
		READ(46,*,ERR=222)
	END DO
	READ(46,*,ERR=222)WATSAT(I,1)
	READ(46,*,ERR=222)WILTPOINT25(I,1)
	MAXWAT(I,1)=MAXWAT(I,1)*(MAXDEPTH/(MAXLAYER1*1.))/25.
	WATSAT(I,1)=WATSAT(I,1)*(MAXDEPTH/(MAXLAYER1*1.))/25.
	WILTPOINT5(I,1)=WILTPOINT25(I,1)*(MAXDEPTH/(MAXLAYER1*1.))/25.
20    FORMAT(F9.1/,16/,F9.1/,F9.1)
	DO 100 IL=2,MAXLAYER1
	  MAXWAT(I,IL)=MAXWAT(I,1)
	  WATSAT(I,IL)=WATSAT(I,1)
	  WILTPOINT5(I,IL)=WILTPOINT5(I,1)
100   CONTINUE
      SNAME(I)=TEMP
      NUMSOIL=I
      REWIND(46)
	RETURN
      PRINT*,'SOIL: Water parameters from SOIL.DAT'
333   CONTINUE
222   CONTINUE
      WRITE(15,*)'Warning! Error in general soil parameters!'
      WRITE(15,*)'Check format of soil parameter file'
C
C Leave PAR_WAT_SWAT
C
      END
C
C-------------------------------------------------------------
C
	SUBROUTINE PET(LATD,JD,T,AVP,SVP,L,WS,DELTAVP,G,PSC,ALBEDO,
     &				CURRENTPLHEI,LAI,CURRENTPLBIOM,ETPOT,ET)
C
C Subroutine to calculate potential evapotranspiration of actual crop (ET)
C
	IMPLICIT NONE
C
C Variables local to this subroutine
C
	REAL L					! latent heat of vaporization (MJ kg-1)
	REAL RN					! net radiation at the crop surface [MJ m-2 day-1]
	REAL G					! soil heat flux density [MJ m-2 day-1]
	REAL WS					! wind speed at 2 m height [m s-1]
	REAL SVP				! saturation vapour pressure [kPa] calculated from mean temperature
	REAL DELTAVP			! slope vapour pressure curve [kPa �C-1]
	REAL PSC				! psychrometric constant [kPa �C-1]
	REAL RA					! diffusion resistance of the air layer (aerodynamic resistance) (s m-1)
	REAL RC					! plant canopy resistance (s m-1)
	REAL COMBTERM			! combind Term K1*(0.622*L*pair/P) = 1710 - 6.85 * T
	REAL ALBEDO				! Albedo 
C 
C Variables passed to/from calling subroutine
C
	REAL LATD				! IN: latitude in degree (user input)
	INTEGER JD				! IN: Julian day
	REAL T					! IN: mean daily air temperature (at 2 m height?) from CRU [�C]
	REAL AVP				! IN: actual vapour pressure from CRU [kPa]
	REAL Z					! IN: elevation above sea level [m] (user input)
	REAL CURRENTPLHEI		! IN: canopy hight in cm
	REAL LAI				! IN: LAI
	REAL CURRENTPLBIOM		! IN: Current plant biomass [ kg ha-1]
	REAL ETPOT				! IN: Penman-Monteith potential evapotranspiration [mm day-1]
	REAL ET					! OUT: Penman-Monteith potential evapotranspiration for actual crop [mm day-1]
C
C ... calculation of ALBEDO (if CURRENTPLBIOM was zero, albedo=0.3 = bare soil)
C
	ALBEDO=0.23 * (1-(EXP(-5.*0.0001*CURRENTPLBIOM)))+		! equation 2.2.14 - 2.2.16
     & 0.3*EXP(-5.*0.0001*CURRENTPLBIOM)						! (SWAT Theoretical Doc.)
C
C calculation of net radiation at the crop surface [MJ m-2 day-1]
C
	CALL NETRADIATION_D(JD,LATD,T,AVP,ALBEDO,RN)
C
C ... diffusion resistance of the air layer (aerodynamic resistance) (s m-1)
C
	CALL AERODYRES(WS,CURRENTPLHEI,RA)
C
C ... plant canopy resistance (s m-1)
C
	CALL CANOPYRES(LAI,RC)
C
C ... K1*(0.622*L*pair/P) = 1710 - 6.85 * T (equation pg. 121 (SWAT Theoretical Doc.) 
C
	COMBTERM = 1710 - 6.85 * T				
C
C calculation of potential evapotranspiration of actual crop (ET)
C
	ET = ((DELTAVP*(RN-G)+PSC*COMBTERM*(SVP - AVP)/RA)/
     &		(DELTAVP+PSC*(1+RC/RA)))/L
C	WRITE(4001,*) DELTAVP,RN,G,PSC,COMBTERM,SVP,AVP,RA,RC
	IF (ET.LT.0) ET = 0.
      ET = MIN(ET,ETPOT)
C
C Leave PET
C
	RETURN
	END
C
C-------------------------------------------------------------
C
	SUBROUTINE PET_ALFALFA(
C INPUTS
     &				LATD,JD,T,AVP,SVP,L,WS,DELTAVP,G,PSC,
C OUTPUTS
     &				ETPOT)
C 
C Subroutine to calculate potential evapotranspiration of a well watered plant (ETPOT)
C SWAT Chapter 7.2.1 (Penman-Monteith method)
C
	IMPLICIT NONE
C
C Variables local to this subroutine
C
	REAL RN					! net radiation at the crop surface [MJ m-2 day-1]
	REAL RA					! diffusion resistance of the air layer (aerodynamic resistance) (s m-1)
	REAL RC					! plant canopy resistance (s m-1)
	REAL COMBTERM			! combind Term K1*(0.622*L*pair/P) = 1710 - 6.85 * T
	REAL ALBEDO				! Albedo
C
C Variables passed to/from calling subroutine
C
	REAL LATD				! latitude in degree (user input)
	INTEGER JD				! Julian day
	REAL T					! mean daily air temperature (at 2 m height) [�C]
	REAL Tmin				! minimum Temperature [�C] if available
	REAL Tmax				! maximum Temperature [�C] if available
	REAL AVP				! actual vapour pressure from CRU [kPa]
	REAL L					! latent heat of vaporization (MJ kg-1)
	REAL WS					! wind speed at 2 m height [m s-1]
	REAL DELTAVP			! slope vapour pressure curve [kPa �C-1]
	REAL G					! soil heat flux density [MJ m-2 day-1]
	REAL PSC				! psychrometric constant [kPa �C-1]
	REAL SVP				! saturation vapour pressure [kPa] calculated from mean temperature
	REAL ETPOT				! Penman-Monteith potential evapotranspiration [mm day-1]
C
	ALBEDO = 0.23
C
C ... net radiation at the crop surface [MJ m-2 day-1]
C
	CALL NETRADIATION_D(JD,LATD,T,AVP,ALBEDO,RN)
C
C ... diffusion resistance of the air layer (aerodynamic resistance) (s m-1)
C
	RA = 114/WS								! equation 7.2.20 (SWAT 
											! Theoretical Doc.)
C ... plant canopy resistance (s m-1)
C
	RC = 49. / (1.4 - 0.4 * 380 / 330.)		! equation 7.2.22 (SWAT 
											! Theoretical Doc.)
C ... K1*(0.622*L*pair/P) = 1710 - 6.85 * T
C
	COMBTERM = 1710 - 6.85 * T				! equation 7.2.19 (SWAT 
											! Theoretical Doc.)
C
	ETPOT = ((DELTAVP*(RN-G)+PSC*COMBTERM*(SVP - AVP)/RA)/
     &		(DELTAVP+PSC*(1+RC/RA)))/L
C
	ETPOT = Max(0., ETPOT)
C
C Leave PET_ALFALFA
C
	RETURN
	END 
C
C-------------------------------------------------------------
C
	SUBROUTINE PLANT_TRANSPIRATION(WMAX,SOILW,CURRENTLAI,ET,ROOTLAYER,
     &			ACTPLANTTRANSP,ACTPLANTTRANSP50,WATCONT)
C 
C Subroutine to calculate actual plant water uptake (SWAT Theoretical Doc., Chapter 7.3 actual evapotranspiration)
C
C DETAILS:
C SOILWUPTAKE(pot)=(ET/(1-EXP(-BETA))*(1-(EXP(-BETA*(Z/ZROOT))))
C
	IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER					! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1					! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER MAXDEPTH					! Maximum depth of the soil profile (cm)
	PARAMETER (MAXDEPTH=300)
	INTEGER I							! internal counter
	REAL SOILWUPTAKE_D(0:MAXDEPTH)		! potential water uptake from the soil surface to a specified depth (mm H2O)
	REAL BETA							! water-use distribution parameter 
	REAL Z								! depth from the soil surface (mm)
	REAL ZROOT							! depth of root development in the soil (mm). = ROOTLAYER*LAYERDEPTH
	REAL LAYERDEPTH						! depth of soil layer (mm)
	REAL EPCO							! plant uptake compensation factor (see explaination below)
	REAL PLANTWUPTAKE(0:MAXLAYER)		! actual water uptake per layer (mm H2O)
	REAL TOTPLANTWUPTAKE				! total actual water uptake (mm H2O)
C
C Variables passed to/from calling subroutine
C
	REAL CURRENTLAI						! IN: leaf area index for a given day
	REAL ET								! IN: Penman-Monteith potential evapotranspiration for actual crop [mm day-1]
	INTEGER ROOTLAYER					! IN: layer until which roots have grown
	REAL WMAX(MAXLAYER)					! IN: Available water at field cap. (mm/layer)
	REAL WILTPOINT(MAXLAYER)			! IN: Total water conent of the soil
										!	  at WP (mm / layer)
	REAL SOILW(MAXLAYER)				! IN: Available water (mm/layer)
	REAL WATCONT(MAXLAYER)				! IN/OUT: Total water content of the soil
										!		layer (mm / layer)
	REAL ACTPLANTTRANSP					! IN/OUT: actual plant transpiration (mm H2O)
	REAL ACTPLANTTRANSP50				! IN/OUT: actual plant transpiration from 50 cm soil profile (mm H2O)
C
C Calculate layer depth (LAYERDEPTH) [mm]
C
	LAYERDEPTH = (MAXDEPTH/MAXLAYER)*10
C
C The water-use distribution parameter, BETA, is set to 10 in SWAT. 
C With this value, 50% of the water uptake will occur in the upper 6% of
C the root zone. Figure 18-3 graphically displays the uptake of water at different
C depths in the root zone.
C
	BETA =10.
C
C	The parameter EPCO: If upper layers in the soil profile do not contain enough water to meet the
C	potential water uptake calculated with equation 18.2.2, users may allow lower
C	layers to compensate. EPCO is the plant uptake compensation factor. The plant uptake compensation
C	factor can range from 0.01 to 1.00 and is set by the user. As EPCO approaches 1.0,
C	the model allows more of the water uptake demand to be met by lower layers in
C	the soil. As EPCO approaches 0.0, the model allows less variation from the depth
C	distribution described by previous equation for potential uptake to take place.
C
	EPCO =0.
C
C 2. Calculate POTENTIAL plant water uptake per layer (lower boundery depth - upper boundery depth)
C
	SOILWUPTAKE_D(0)=0.0
	PLANTWUPTAKE(0)=0.0
C loop over layers:
C
	do I=1,ROOTLAYER	! loop over profile down to rooting depth
C
C	 1. Calculation of pot. water uptake to specific depth (e.g. 5 cm, 10 cm...):
C
	SOILWUPTAKE_D(I)=(ET/(1-EXP(-BETA))*(1-(EXP(-BETA*((I*LAYERDEPTH)/		! equation 18.2.1 (SWAT 
     &	(ROOTLAYER*LAYERDEPTH))))))											! Theoretical Doc.)	
C
C	 2. Calculate plant water uptake according to availabe water in the soil (ACTUAL plant water uptake):
C
	! Details: Subtract pot. uptake down to previous from down to current + 
	!				excess demand (if compensation allowed)
	!				excess = difference between pot. uptake and total water taken up up to previous layer
C
	PLANTWUPTAKE(I)=SUM(SOILWUPTAKE_D(1:I))-SUM(SOILWUPTAKE_D(1:I-1))+
     &				(SUM(SOILWUPTAKE_D(1:I-1))-
     &				SUM(PLANTWUPTAKE(1:I-1)))*EPCO							! equation 18.2.3 (SWAT 
C																			! Theoretical Doc.)
C	Adjust uptake if soil water is less than 25% of plant available water:
C	As the water content of the soil decreases, the water in the soil is held
C	more and more tightly by the soil particles and it becomes increasingly difficult
C	for the plant to extract water from the soil. To reflect the decrease in the
C	efficiency of the plant in extracting water from dryer soils, the potential water
C	uptake is modified using the following equations:
C 
	IF(WATCONT(I).LT.(0.25*WMAX(I)))THEN
		PLANTWUPTAKE(I)=PLANTWUPTAKE(I)*EXP(5*(((WATCONT(I))/
     &	 (0.25*WMAX(I)))-1))												! equation 18.2.4 (SWAT 
																			! Theoretical Doc.)										
	ELSE
		PLANTWUPTAKE(I)=PLANTWUPTAKE(I)										! equation 18.2.5 (SWAT 
																			! Theoretical Doc.)	
	END IF
																	
																			
C	3.update soil water status
C
	IF(PLANTWUPTAKE(I).GE.(WATCONT(I)-WILTPOINT(I)))THEN
		PLANTWUPTAKE(I)=WATCONT(I)-WILTPOINT(I)
	END IF
	
	WATCONT(I)=MAX(WATCONT(I)-PLANTWUPTAKE(I),WILTPOINT(I))
	END DO
C
C 4. Calculate plant actual transpiration
C
C 4a) calculate total plant water uptake over all the layers
C
	TOTPLANTWUPTAKE=SUM(PLANTWUPTAKE(1:ROOTLAYER))
	ACTPLANTTRANSP50=SUM(PLANTWUPTAKE(1:10))
C
C 4b) calculate actual plant transpiration
C
	ACTPLANTTRANSP=TOTPLANTWUPTAKE																		
C
C Leave PLANT_TRANSPIRATION
C
	RETURN
	END
C
C-------------------------------------------------------------
C
	SUBROUTINE PSYCHO_D(Z,PSC)

C Subroutine to calculate the psychrometric constant [kPa �C-1]
C
	IMPLICIT NONE
C
C Variables local to this subroutine
C
	REAL P			! atmospheric pressure [kPa]
	REAL L			! latent heat of vaporization, 2.45 [MJ kg-1]
	REAL cp			! specific heat at constant pressure, 1.013 10-3 [MJ kg-1 �C-1]
	REAL e			! ratio molecular weight of water vapour/dry air = 0.622.
C
C Variables passed to/from calling subroutine
C
	REAL Z			! IN: elevation above sea level [m]
	REAL PSC		! OUT: psychrometric constant [kPa �C-1]
C
	cp = 1.013e-3
	e=0.622
	L=2.45
C	
	P=101.3*((293-0.0065*Z)/293)**5.26	! FAO, Chapter 3, 
										! http://www.fao.org/docrep/X0490E/x0490e07.htm#atmospheric%20pressure%20(p)
C
	PSC= (cp*P)/(e*L)					! equation 3.3.7 (SWAT 
										! Theoretical Doc.)
C
C Leave PSYCHO_D
C
	RETURN
	END
C
C-------------------------------------------------------------
C
	SUBROUTINE READWATER(WMAX,FIELDCAP)

C Subroutine to read in max. available water and wilting point

C
C Variables local to this subroutine
C
	INTEGER MAXLAYER		! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1		! No.of layers in the soil profile
	INTEGER MAXDEPTH		! Depth of profile (cm)
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
	INTEGER I				! Local counter variable
	REAL SOILWP(MAXLAYER)	! soil water content at WP [mm H2O / 5 cm layers]
	REAL SOIL2(MAXLAYER)	! plant available water (2.0-4.2 pF) [mm H2O / 5 cm layers]
	REAL SOIL23(MAXLAYER)	! plant available water (2.3-4.2 pF) [mm H2O / 5 cm layers]
	REAL SOIL25(MAXLAYER)	! plant available water (2.5-4.2 pF) [mm H2O / 5 cm layers]
C
C Variables passed to/from calling subroutine
C
	REAL WMAX(MAXLAYER)		! IN:Available water at field cap. (mm/layer)
	REAL FIELDCAP(MAXLAYER)	! OUT: soil water content at FC (mm/layer)
C
	OPEN(45,FILE='SOILWATER.prn', ACTION="READ",STATUS="OLD",ERR=555)
	DO I=1,MAXLAYER
		READ(45,*)	SOIL2(I), SOIL23(I), SOIL25(I), SOILWP(I)
		WMAX(I) = SOIL2(I)
		FIELDCAP(I) = SOIL2(I)+SOILWP(I)
	END DO
C
C Leave READWATER
C
      RETURN
555   END
C
C-------------------------------------------------------------------
C
      SUBROUTINE ROOTL_SWAT(M,MCROP,IROCKS,ROOTS)
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
      REAL ROOTS
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
C-------------------------------------------------------------
C
	SUBROUTINE SATVAPOURPRESSURE_D(AVP,T,SVP)
C
C Subroutine to calculate saturation vapour pressure [kPa]
C
	IMPLICIT NONE
C
C Variables passed to/from calling subroutine
C
	REAL T				! IN: mean daily air temperature (at 2 m height) [�C]
	REAL AVP			! IN: water vapor pressure of air at height Z (kPa)
	REAL SVP			! OUT: saturation vapor pressure of air at height Z (kPa)

C As saturation vapour pressure is related to air temperature, it can be calculated 
C from the air temperature. The relationship is expressed by:
C
	SVP = 0.6108 * EXP((17.27 * T)/(T+237.3))	! FAO, Chapter 3,
												! http://www.fao.org/docrep/X0490E/x0490e07.htm#air%20humidity
C
C Due to the non-linearity of the above equation, the mean saturation vapour pressure 
C for a day, week, decade or month should be computed as the mean between the saturation 
C vapour pressure at the mean daily maximum and minimum air temperatures for that period
C
C sVPminmaxT = ((0.6108 * EXP((17.27 * Tmin)/(Tmin+237.3)))+(0.6108 * EXP((17.27 * Tmax)/(Tmax+237.3))))/2
C
C Leave SATVAPOURPRESSURE_D
C
	RETURN
	END 
C
C-------------------------------------------------------------
C
	SUBROUTINE SLOPEVAPOURPRESSURE_D(T,SVP,DELTAVP)
C
C Subroutine to calculate the slope of the relationship between saturation vapour pressure and temperature
C
	IMPLICIT NONE
C
C Variables passed to/from calling subroutine
C
	REAL T					! IN: mean daily air temperature (at 2 m height) [�C]
	REAL SVP				! IN: saturation vapor pressure of air at height Z (kPa)
	REAl DELTAVP			! OUT: slope vapour pressure curve [kPa �C-1]
C
	DELTAVP = (4098 * SVP) / ((T+237.3)**2)
C
C ... if Tmin and Tmax are available!!!
C
C DELTAVP = (4098 * sVPminmaxT) / ((T+237.3)**2)
C
C Leave SLOPEVAPOURPRESSURE_D
C
	RETURN
	END
C
C-------------------------------------------------------------
C
	SUBROUTINE SOIL_EVAP(WATCONT,FIELDCAP,WILTPOINT,CURRENTPLBIOM,
     &			SOILW,ET_DEMAND,ACTPLANTTRANSP,ES)
C 
C Subroutine to calculate soil evapotranspiration (SWAT Theoretical Doc., Chapter 7.3.3)
C
	IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	PARAMETER (MAXDEPTH=300)
	INTEGER I					! internal counter
	REAL ESMAX					! maximum sublimation/soil evaporation on a given day (mm H2O)
	REAL ESMAX_ADJ				! maximum sublimation/soil evaporation adjusted for plant water use on a given day (mm H2O)
	REAL LAYERDEPTH				! depth of soil layer (mm)			
	REAL SOILCOVER				! soil cover index
	REAL ESMAX_D(0:MAXLAYER)	! evaporative demand at depth Z (mm H2O)
	REAL ESMAX_L(MAXLAYER)		! evaporative demand per layer (mm H2O)
	REAL ESCO					! soil evaporation compensation coefficient
	REAL ES_ADJ_L(MAXLAYER)		! evaporative demand for layer ly adjusted for water content (mm H2O)
	REAL ES_L(MAXLAYER)			! amount of water removed from layer by evaporation per layer (mm H2O)
	REAL ES						! total soil evaporation (mm H2O)
	REAL ESMAXProfile			! SUM(ESMAX_L)
	REAL ES_ADJ_PROFILE			! SUM(ES_ADJ_L)
C
C Variables passed to/from calling subroutine
C
	REAL SOILW(MAXLAYER)		! IN: Available water (mm/layer)
	REAL ET_DEMAND				! IN: evaporative demand after intercepted water is evaporated [mm day-1]
	REAL ACTPLANTTRANSP			! IN: actual plant transpiration (mm H2O)
	REAL CURRENTPLBIOM			! IN: Current plant biomass [ kg ha-1]
	REAL FIELDCAP(MAXLAYER)		! IN: soil water content at FC [mm]
	REAL WILTPOINT(MAXLAYER)	! IN: Total water conent of the soil
								!		layer at WP (mm / layer)
	REAL WATCONT(MAXLAYER)		! IN/OUT: Total water content at 
								!		saturation (mm / layer)
C
C Calculate layer depth (LAYERDEPTH) [mm]
C
	LAYERDEPTH = (MAXDEPTH/MAXLAYER)*10
C
C Calculation of FIELDCAP for testing (Lena data set), WP=12.8v%
C 1. Available water = 57.49 (25cm)
C 2. WP = 12.8%/100*250 mm= 32
C 3. WP+availabe water = 89.49
C
C 1. The amount of sublimation and soil evaporation will be impacted by the
C degree of shading. The maximum amount of sublimation/soil evaporation on a
C given day is calculated as:
C where SOILCOVER soil the cover index. The soil cover index is calculated:
C
	SOILCOVER = EXP(-5.*0.0001*CURRENTPLBIOM)
C
	ESMAX = ET_DEMAND* SOILCOVER
C
	ESMAX_ADJ=MIN(ESMAX,(ESMAX*ET_DEMAND)/(ESMAX+ACTPLANTTRANSP+1.e-10))
	ESMAX_ADJ=MAX(ESMAX_ADJ,0.)
C
C
C 2. Partitioning soil evaporative demand between layers. First, the maximum amount 
C of water that is allowed to be evaporated is determined :
C The coefficients in this equation were selected so
C that 50% of the evaporative demand is extracted from the top 10 mm of
C soil and 95% of the evaporative demand is extracted from the top 100 mm
C of soil.
C SWAT does not allow a different layer to compensate for the
C inability of another layer to meet its evaporative demand. The evaporative
C demand not met by a soil layer results in a reduction in actual
C evapotranspiration.
C A coefficient has been incorporated into equation 7.3.16 to allow
C the user to modify the depth distribution used to meet the soil evaporative
C demand.
C
C Setting ESCO (As the value for ESCO is reduced, the model is able to extract more
C of the evaporative demand from lower levels.)
C
	ESCO=1
C
	ESMAX_D(0)=0.0
	DO I=1,MAXLAYER
	ESMAX_D(I) = ESMAX_ADJ*((I*50)/((I*50)+EXP(2.374-0.00713*(I*50))))
	ESMAX_L(I) = ESMAX_D(I)-(ESMAX_D(I-1)*ESCO)					! equation 7.3.17 (SWAT 
																! Theoretical Doc.)
	END DO
	ESMAXProfile=SUM(ESMAX_L)
C
C When the water content of a soil layer is below field capacity, the
C evaporative demand for the layer is reduced according to the following
C equations:

	DO I=1,MAXLAYER
	IF (WATCONT(I).LT. FIELDCAP(I)) THEN
		ES_ADJ_L(I)=ESMAX_L(I) * EXP(((2.5*(WATCONT(I)-FIELDCAP(I)))
     &												/FIELDCAP(I)))		! equation 7.3.18 (SWAT 
																		! Theoretical Doc.)
C Error in SWATuserTheory, pg. 133, above equation changed according to equation of SWAT code in etact.f
	ELSE
		ES_ADJ_L(I)=ESMAX_L(I)											! equation 7.3.19 (SWAT 
																		! Theoretical Doc.)
	END IF
	END DO
	ES_ADJ_PROFILE=sum(ES_ADJ_L)
C
C In addition to limiting the amount of water removed by
C evaporation in dry conditions, SWAT defines a maximum value of water
C that can be removed at any time. This maximum value is 80% of the plant
C available water on a given day
C
	DO I=1,MAXLAYER
C	ES_L(I) = MIN(ES_ADJ_L(I),0.8*(WATCONT(I)-WILTPOINT(I)))
	ES_L(I) = ES_ADJ_L(I)
	WATCONT(I)=WATCONT(I)-ES_L(I)
	END DO
	ES = sum(ES_L)
C
C Leave SOIL_EVAP
C
	RETURN
	END
C
C---------------------------------------------------------
C
      SUBROUTINE SWATER_SWAT(IAWC,MAXWAT,WSAT,WATSAT,WMAX,WILTPOINT,
     &				FIELDCAP,WATCONT,WILTPOINT5,SATWATCONT)
C
C Subroutine to set the water holding capacity
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXSOIL			! Max.no.of soil types
	PARAMETER (MAXSOIL=50)
	INTEGER MAXLAYER		! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1		! No.of layers in the soil profile
	INTEGER MAXDEPTH		! Depth of profile (cm)
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
	INTEGER IL,IL1,IL2		! Local counter variable
	REAL CHWAT				! Change in water content with depth
C
C Variables passed to/from calling subroutine
C
	INTEGER*4 IAWC			! IN:Water movement code number
	REAL WMAX(MAXLAYER)		! IN:Available water at field cap. (mm/layer)
	REAL WILTPOINT5(MAXSOIL,MAXLAYER) ! IN: Total water conent of the soil
								!		layer at WP (mm / layer)
	REAL WILTPOINT(MAXLAYER)	! IN: Total water conent of the soil
								!		layer at WP (mm / layer)
	REAL FIELDCAP(MAXLAYER)	! OUT: soil water content at FC (mm/layer)
	REAL WATCONT(MAXLAYER)
	REAL MAXWAT(MAXSOIL,MAXLAYER)	! IN/OUT: Maximum water content (mm/layer)
	REAL WATSAT(MAXSOIL,MAXLAYER)	! IN:Avail.water at saturation (mm/layer)
	REAL WSAT(MAXLAYER)		! IN:Available water at saturation (mm/layer)
	REAL SATWATCONT(MAXLAYER)	! IN: Total water content at 
								!		saturation (mm / layer)
C
C
C Set values of available water : Soil layers 1...10 (0-5,5-10...45-50cm),
C                                 set at max. water holding capacity/10.
C                                 Soil layers 11 & 12 (50-100,100-150cm),
C
C Set values of available water 
C
      DO 101 IL=1,MAXLAYER1
        WMAX(IL)=(MAXWAT(IAWC,IL))
        WSAT(IL)=(WATSAT(IAWC,IL))
	  FIELDCAP(IL)=MAXWAT(IAWC,IL)+WILTPOINT5(IAWC,IL)
	  WATCONT(IL)=MAXWAT(IAWC,IL)+WILTPOINT5(IAWC,IL)
	  WILTPOINT(IL)= WILTPOINT5(IAWC,IL)
	  SATWATCONT(IL)=WSAT(IL)+WILTPOINT(IL)
101   CONTINUE
C
C	CALL READWATER(WMAX,FIELDCAP)
C
C Leave SWATER
C
      RETURN
      END
C
C-------------------------------------------------------------
C