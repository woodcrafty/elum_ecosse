C------------------------------------------------------------
C
C Site routines for Rothamsted Carbon and Nitrogen Turnover Model
C
C Jo Smith 
C Started 16/01/07
C
C
C
C*************************************************************
C MAIN RUN ROUTINES
C*************************************************************
C
C EXTERNAL SUBROUTINES
C ECOSSE_SITE_RUN()
C-------------------------------------------------------------
C
      SUBROUTINE ECOSSE_SITE_RUN(ISWAIT)
C
C Subroutine to run ECOSSE for SITE
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
	INTEGER MAXCROP				! Max.no.of crop types allows
      PARAMETER (MAXCROP=36)
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/
	INTEGER MAXFERT				! Max.no.of fertiliser applications allowed 
	PARAMETER (MAXFERT=5)
	INTEGER MAXGROW				! Max.no.of growing seasons allowed
	PARAMETER (MAXGROW=300)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXORGM				! Max.no.of organic manure applications allowed
	PARAMETER (MAXORGM=52)
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXSOIL				! Max.no.of soil types
	PARAMETER (MAXSOIL=50)
	INTEGER MAXSOMLAY			! Max.no.of soil layers entered			
	PARAMETER (MAXSOMLAY=10)
	INTEGER MAXWEATH			! Max.no.of years allowed
	PARAMETER (MAXWEATH=300)

	INTEGER IL					! Local counter variable
	INTEGER NF					! No.fertiliser application
	INTEGER I1TEMP
C
C Variables passed to/from other subroutines
C (IN = originates in other subroutine; 
C  OUT = originates in this subroutine or subroutine called here; 
C  Name in brackets shows first definition of term;
C  CALL indicates passed to/from calling subroutine)
C
	REAL AE						! IN(RUN2_MAGEC_CROP): Weather data used by MAGEC
	REAL AIRTEMP				! IN(GETWEATHER_AT_FC): Air temperature (deg.C/timestep)	
	REAL AMMN(MAXLAYER)			! IN(INIT_SUNDIAL_SOILCN_NOPARS):Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN(INIT_SUNDIAL_SOILCN_NOPARS):Soil ammonium-N15 (kgN15/ha/layer)
      REAL ATM					! IN(SETFILE):Atmospheric N input (kgN/ha/timestep)
	REAL ATM15					! IN(RUN2_SUNDIAL_CROP):Atmospheric N15 input (kgN15/ha/time)
      REAL AVEPET(12)				! IN(GET_SITE_SOIL):Long term average PET (mm/month)
      REAL AVERAIN(12)			! IN(GET_SITE_SOIL):Long term average rainfall(mm/month)
      REAL AVETEMP(12)			! IN(GET_SITE_SOIL):Long term average monthly average 
	REAL AVP					! IN(GETWEATHER_AT_FC_SWAT): water vapor pressure of air at height z (kPa)
      REAL BALANCE(20)			! IN(SETRES):C Results
	REAL BBRADS					! IN(RUN2_MAGEC_CROP): Weather data used by MAGEC
	REAL BCARB(MAXLAYER)		! IN(MICROBIAL_SUNDIAL_SOILCN):C in BIO at start of timestep
	REAL BCARB0(MAXLAYER)		! IN(INIT_SUNDIAL_SOILCN_NOPARS):C in soil biomass (kgC/ha/layer)
      REAL BNIT0(MAXLAYER)		! IN(INIT_SUNDIAL_SOILCN_NOPARS):N in soil biomass (kgN/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN(INIT_SUNDIAL_SOILCN_NOPARS):N15 in soil humus (kgN15/ha/layer)
	REAL BSTART					! IN(STARTRES):Initial N in biomass pool (kgN/ha)
	REAL BSTART15				! IN(STARTRES):Initial N in biomass pool (kgN15/ha)
	REAL BULKDENS(MAXLAYER)		! IN(GET_SITE_SOIL):Bulk density of the layer (g/cm3)
	REAL C_TS					! IN(RUN1_SUNDIAL_CROP):Litter C input in this timestep (kgC/ha)
	REAL C1(0:MAXCROP)			! IN(INIT_SUNDIAL_CROP):Crop parameter for N uptake
      REAL CACCUM					! IN(GET_SITE_SOIL): Measured annual C accumulation 
								!     in the soil kgC/ha/yr
	REAL CACT					! IN(RUN2_SUNDIAL_CROP):Crop N uptake (kgN/ha/timestep)
	REAL CACT15					! IN(RUN2_SUNDIAL_CROP):Crop N15 uptake (kgN15/ha/timestep)
	REAL CACTOT					! IN(RUN1_SUNDIAL_CROP):N taken up by crop (kgN/ha)
	REAL CANMAX					! IN(INIT_SUNDIAL_CROP):Maximum amount of water that can be trapped in the canopy when the canopy is fully developed (mm H2O)
	REAL CANSTOR				! IN(INIT_SUNDIAL_WATER_SWAT):amount of free water held in the canopy on a given day (mm H2O)
	REAL CATOT15				! IN(RUN1_SUNDIAL_CROP):N15 taken up by crop (kgN15/ha)
	REAL CFACT(0:MAXCROP)		! IN(GETYRSET):Crop parameter used to estimate yield
	REAL CH4(MAXLAYER)			! IN(MICROBIAL_SUNDIAL_SOILCN):CH4 emitted (kgC/ha/layer) in the Aitkenhead methane model

	INTEGER CH4MODEL			! IN(READ_MODEL):Methane model (off, Richards or Aitkenhead model) 
	INTEGER CH4_OFF				! CH4 model off
	INTEGER CH4_RICHARDS    	! Richards CH4 model on
	INTEGER CH4_AITKENHEAD   	! Aitkenhead CH4 model on		
      DATA CH4_OFF,CH4_RICHARDS,CH4_AITKENHEAD /0,1,2/
	REAL CH4TOAIR				! IN(MICROBIAL_SUNDIAL_SOILCN):CH4 released to atmosphere (kgC/ha) from the Aitkenhead methane model
	REAL CH4_PROD(MAXLAYER)  	! CH4 produced in each soil layer [kgC/ha/layer] (Richards CH4 model)
	REAL CH4_SOIL_OX			! CH4 oxidised from methane produced in the soil [kgC/ha] (Richards CH4 model)
      REAL CH4_ATMOS_OX			! CH4 oxidised from atmospheric CH4 [kgC/ha] (Richards CH4 model)
      REAL CH4_FLUX				! Net CH4 flux from soil column [kgC/ha] (Richards CH4 model)
      REAL CLAY(MAXLAYER)			! IN(GET_SITE_SOIL):Clay content of the layer (%)
	INTEGER CNMODEL				! IN(READ_MODEL):Type of Model for C:N ratio of DPM and 
	                            !     RPM (MAGEC or Foereid)
	 INTEGER CNFOEREID			!   C:N ratio obtained by method of Foereid
	 INTEGER CNMAGEC			!   C:N ratio obtained by method of MAGEC
	 DATA CNMAGEC,CNFOEREID /1,2/
	REAL CLOSSX					! IN(SETWEEKVARS):N lost by senescence (kgN/ha) 
	REAL CLOSSX15				! IN(SETWEEKVARS):N15 lost by senescence (kgN15/ha) 	
	REAL CO2(MAXLAYER)			! IN(MICROBIAL_SUNDIAL_SOILCN):CO2 emitted (kgC/ha/layer)
	REAL CO2FROMDOC(MAXLAYER)	! IN(MICROBIAL_SUNDIAL_SOILCN):CO2 output from DOC in each layer(kgC/ha)
	REAL CONC					! IN(PHYSICAL_SUNDIAL_SOILCN):Conc.N leached (kgN/ha/mm/layer)
	REAL CONC15					! IN(PHYSICAL_SUNDIAL_SOILCN):Conc.N15 leached (kgN/ha/mm/layer)
	REAL CONVER_F				! IN(INIT_SUNDIAL_CROP):Conversion between timestep & weeks
	REAL CRIT(MAXSOIL,MAXLAYER)	! IN(GET_SUNDIAL_SOILPARS):Critical N content (kgN/ha/layer)
	REAL CRITMIN				! IN(RUN2_MAGEC_CROP): Data used by MAGEC
	REAL CTOT					! IN(RUN1_SUNDIAL_CROP):Total litter C input (kgC/ha)
	REAL CULTDEPTH				! OUT(CULTIV):Cultivation depth (cm)
      INTEGER CULTIVATE			! OUT(RUN1_SUNDIAL_CROP): Code to cultivate this week (1=yes, 0=no)
	REAL CUPTN					! IN(RUN2_SUNDIAL_CROP):Crop N uptake (kgN/ha/timestep)
	REAL CURRENTLAI				! IN(INIT_SUNDIAL_WATER_SWAT):Current LAI
	REAL CURRENTPLBIOM			! IN(INIT_SUNDIAL_WATER_SWAT):Current plant biomass [ kg ha-1]
	REAL CURRENTPLHEI			! IN(INIT_SUNDIAL_WATER_SWAT):Current plant height [cm]
	REAL DAVTMP					! IN(RUN2_MAGEC_CROP): Weather data used by MAGEC
	REAL DDAYS					! IN(RUN2_MAGEC_CROP): Weather data used by MAGEC
	REAL DENITRIFN(MAXLAYER)	! IN(MICROBIAL_SUNDIAL_SOILCN): N lost by denitrification (kgN/ha/layer
	REAL DENIT					! IN(MICROBIAL_SUNDIAL_SOILCN):N lost by denitrification (kgN/ha)
	INTEGER DMODEL				! IN(READ_MODEL):Type denitrification model chosen 
	                            ! (Bradbury or NEMIS)
	 INTEGER DBRADBURY			!   Bradbury model for denitrification
	 INTEGER DNEMIS				!   Nemis model for denitrification
	 DATA DBRADBURY,DNEMIS /1,2/
	REAL DN15					! IN(MICROBIAL_SUNDIAL_SOILCN):N15 lost by denitrification (kgN15/ha)
	REAL DNIT(MAXLAYER)			! IN(MICROBIAL_SUNDIAL_SOILCN):Net Mineralised N (kgN/ha/lay/timestep)
	REAL DNIT15(MAXLAYER)		! IN(MICROBIAL_SUNDIAL_SOILCN):Net Mineralised N15 (kgN15/ha/lay/time)
 	INTEGER DOCMODEL			! IN(READ_MODEL):DOC model (on or off)
	 INTEGER DOC_OFF			!   DOC model off
	 INTEGER DOC_ON				!   DOC model on
	 DATA DOC_ON,DOC_OFF /1,2/
	REAL DOMSOILIMPDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN(GET_SITE_SOIL):Depth of impermeable layer if present (cm)
	INTEGER DOMSOILISIMP(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN(GET_SITE_SOIL):Flag for impermeable layer (0=No; 1=Yes)
      REAL DPM_RPM_RATIO			! IN(SET_DPMRPMRATIO):Ratio of DPM:RPM. If set to 
								!        zero, this will be worked out using a call to SET_DPMRPMRATIO
	REAL DPMCARB(MAXLAYER)		! IN(MICROBIAL_SUNDIAL_SOILCN):C in DPM at start of timestep
	REAL DPMCARB0(MAXLAYER)		! IN(INIT_SUNDIAL_SOILCN_NOPARS):C in decomposable PM (kgC/ha/layer)
      REAL DPMCTON(MAXLAYER)		! IN(GET_CTON_FOEREID):DPM C:N ratio 
	REAL DPMNIT0(MAXLAYER)		! IN(INIT_SUNDIAL_SOILCN_NOPARS):N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN(INIT_SUNDIAL_SOILCN_NOPARS):N15 in decomposable PM (kgN15/ha/layer)
	REAL DRAINW(MAXLAYER)		! IN(DRAIN_SUNDIAL_WATER):Water drained from this layer (mm/layer)
	REAL DRRAT					! IN(SET_DPMRPMRATIO):DPM:RPM ratio
	REAL DTR					! IN(RUN2_MAGEC_CROP): Weather data used by MAGEC
	REAL DTRJM2					! IN(RUN2_MAGEC_CROP): Weather data used by MAGEC
	REAL DVP					! IN(RUN2_MAGEC_CROP): Weather data used by MAGEC
	INTEGER EC_EQRUN			! IN(READ_MODEL):Initialisation using a full ECOSSE equilibrium run (on or off)
	 INTEGER EC_EQRUN_OFF		!   Initialisation using a full ECOSSE equilibrium run is off
	 INTEGER EC_EQRUN_ON		!   Initialisation using a full ECOSSE equilibrium run is on
	 DATA EC_EQRUN_OFF,EC_EQRUN_ON /0,1/
	INTEGER ENTER_WTD			! OUT: Code to use read in water table depth 1 = read in, 0=dont read in
	INTEGER EQMODEL				! IN(READ_MODEL):Type of equilibrium run 
	                            !   (NPP, TOC or both)
	 INTEGER EQNPP				!   Model initialised using plant inputs 
								!      from measured NPP
	 INTEGER EQTOC				!   Model initialised using measured TOC
	 INTEGER EQNPPTOC			!   Model initialised using both plant inputs
								!      from measured NPP and measured TOC
	 INTEGER EQHILLIER,EQJONES			!   Model initialised using TOC and
								!      the Hillier solver to get C pools and plant inputs
	 DATA EQNPP,EQTOC,EQNPPTOC,EQHILLIER,EQJONES /1,2,3,4,5/
	REAL EVAP					! IN(GETWEATHER_AT_FC): Potential evap. (mm/timestep)
	REAL EVAPO					! IN(RUN2_MAGEC_CROP): Data used by MAGEC
	REAL EVAPW					! IN(GETWEATHER): Potential evap. (mm/timestep)
	REAL EXYLD					! IN(SETFILE):Yield of current crop (t/ha)
      REAL EXYLDJ(0:MAXGROW)		! IN(SETFILE):Yield of current crop (t/ha)
								! (0->calculate within model)
	REAL FERT(MAXFERT)			! IN(GETYRSET):Amount of fertiliser applied (kgN/ha)
	REAL FERTJ(0:MAXGROW,MAXFERT)	! IN(SETFILE):Amount of fertiliser applied (kgN/ha)
	REAL FERTNF					! IN(GETFERT):Amount of fertiliser applied (kgN/ha)
      REAL FIELDCAP(MAXLAYER)		! IN(INIT_SUNDIAL_WATER): Total water in the layer at 
	                            !		field capacity (mm/layer)
	INTEGER FIXEND				! IN(SETFILE):Fixed end? 0=No 1=Yes
	CHARACTER*100 FIXFILE(MAXWEATH)	! IN(SETFILE): Weather files
	CHARACTER*100 FIXFILEAVP(MAXWEATH)	! IN(SETFILE): Weather files (actual vapour pressure) (kPa)
	REAL FIXN					! IN(SETWEEKVARS):N fixed from atmosphere (kgN/ha)
	REAL FLOWPROP				! IN(INIT_SUNDIAL_WATER):Proportion of flow needed to achieve
							    !	     water table at depth WTABLE
	REAL FYMFERT				! IN(GETFERT):N input by fertiliser & FYM (kgN/ha)
	REAL FYMFERT15				! IN(GETFERT):N15 input by fertiliser & FYM (kgN15/ha)
	REAL GNN2O					! IN(MICROBIAL_SUNDIAL_SOILCN):N2O lost by nitrification (kgN/ha)
	REAL GNNO					! IN(MICROBIAL_SUNDIAL_SOILCN):NO lost by nitrification (kgN/ha)
	REAL GPNN2O					! IN(MICROBIAL_SUNDIAL_SOILCN):N2O lost by part.nitrification (kgN/ha)
	REAL G15NN2O				! IN(MICROBIAL_SUNDIAL_SOILCN):15N2O lost by nitrification (kgN15/ha)
	REAL G15NNO					! IN(MICROBIAL_SUNDIAL_SOILCN):15NO lost by nitrification (kgN15/ha)			
	REAL G15PNN2O				! IN(MICROBIAL_SUNDIAL_SOILCN):15N2O lost by part.nitrif.(kgN15/ha)			
	REAL GDN2					! IN(MICROBIAL_SUNDIAL_SOILCN):N2 lost by denitrification (kgN/ha)
	REAL GDN2O					! IN(MICROBIAL_SUNDIAL_SOILCN):N2O lost by denitrification (kgN/ha)
	REAL G15DN2					! IN(MICROBIAL_SUNDIAL_SOILCN):15N2 lost by denitrification (kgN15/ha)
	REAL G15DN2O				! IN(MICROBIAL_SUNDIAL_SOILCN):15N2O lost by denitrification (kgN15/ha)					
	REAL HCARB(MAXLAYER)		! IN(MICROBIAL_SUNDIAL_SOILCN):C in HUM at start of timestep
	REAL HCARB0(MAXLAYER)		! IN(INIT_SUNDIAL_SOILCN_NOPARS):C in soil humus (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN(INIT_SUNDIAL_SOILCN_NOPARS):N in soil humus (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN(INIT_SUNDIAL_SOILCN_NOPARS):N15 in soil biomass (kgN15/ha/layer)
	REAL HZ1(MAXLAYER)			! IN(INIT_SUNDIAL_SOILCN_NOPARS):N:C of input PM for steady state 
	INTEGER I_TS				! IN(RUN1_SUNDIAL_CROP):No.timesteps since sowing
	INTEGER IANTHES				! IN(INIT_SUNDIAL_CROP):No.timesteps from 01/01 to anthesis 
	INTEGER IAWC				! IN(SETFILE):Water movement code number
	INTEGER IAWCJ				! IN(SETFILE):Water movement code number
  	REAL ICFACTOR(MAXLAYER)		! IN(INIT_SUNDIAL_SOILCN_NOPARS):Adjustment in rate needed to	
								!    achieve measured NPP and TOC	
	INTEGER ICMODEL				! IN(READ_MODEL):Type of C initialisation model 
	 INTEGER ICFIXED			!   Initialisation of C is fixed
	 INTEGER ICROTHCEQ			!   Initialisation of C by RothC 
								!    equilibrium run
	 DATA ICFIXED,ICROTHCEQ /1,2/
      INTEGER ICOVER				! IN(RUN1_SUNDIAL_CROP):Crop cover 1=Covered 0=Bare
	INTEGER ICROP				! IN(GETYRSET):Current crop code
	INTEGER ICROPJ(0:MAXGROW)	! IN(SETFILE):Crop codes
	INTEGER ICULT(MAXGROW)		! IN(SETFILE): Timesteps from 01/01 when cultivation occured
	INTEGER IDAG					! IN(RUN2_MAGEC_CROP): Weather data used by MAGEC
	INTEGER IDATEFC				! IN(SETFILE):Date of field capacity (1=01/01; 2=01/06)
	INTEGER IDRAINJ				! IN(SETFILE):Water movement code number 
	INTEGER IEND				! IN(INIT_SUNDIAL_CROP):Timesteps from sowing date to harvest
	INTEGER IFERT(MAXFERT)		! IN(GETYRSET):No.timesteps to split fert.application
      INTEGER IFERTJ(0:MAXGROW,MAXFERT)	! IN(SETFILE):No.timesteps to fert.application
	INTEGER IFILE				! IN(GETYRSET):Current weather file
	INTEGER IHARV				! IN(GETYRSET):Timesteps from 01/01/01 to harvest date 
	INTEGER IHARVJ(0:MAXGROW)	! IN(SETFILE):Timesteps from 01/01/01 to harvest date 
	INTEGER IK					! OUT(GETWEATHER):No.timesteps from prev.harv. to current
	INTEGER IL_TSS				! IN(RUN2_SUNDIAL_CROP):No.timesteps before harvest when senesces
	INTEGER ILAB(MAXFERT)		! IN(GETYRSET):Labelling on fertiliser? 0=No 1=Yes
	INTEGER ILABJ(0:MAXGROW,MAXFERT)	! IN(SETFILE):Lab.on fertiliser? 0=No 1=Yes
	INTEGER ILABNF				! IN(GETFERT):Labelling on fertiliser? 0=No 1=Yes
	INTEGER ILAST_TS			! IN(SETFILE):Last simulated timestep since 01/01/01
	INTEGER ILU					! IN(GET_SITE_SOIL):Counter for land use 
	INTEGER IMFUNC				! IN(READ_MODEL):Choice of moisture rate modifier 
								!   (ROTHC or HADLEY)
	 INTEGER IMFUNC_ROTHC		!   ROTHC moisture rate modifier
	 INTEGER IMFUNC_HADLEY		!   HADLEY moisture rate modifier
       DATA IMFUNC_ROTHC,IMFUNC_HADLEY /0,1/
	INTEGER INMODEL				! IN(READ_MODEL):Type of N initialisation model
	 INTEGER INSTABLECN			!   Initisalisation of N by assuming 
								!    steady state (after Bradbury)
	 INTEGER INPASSCN			!   Initialisation of N by passing the
								!    C:N ratio of DPM and RPM
	INTEGER INSTRAW				! IN(GETYRSET):Straw incorporated? 0=No 1=Yes
	INTEGER INVERT				! OUT(CULTIV):Invert / mix soil on cultivation? 0=No 1=Yes
	INTEGER IOLAB(MAXORGM)		! IN(GETYRSET):Labelling on manure? 0=No 1=Yes
      INTEGER IOLABJ(0:MAXGROW,MAXORGM)	! IN(SETFILE):Labelling on manure? 0=No 1=Yes
	INTEGER IOLABNF				! IN(GETFYM):Labelling on manure? 0=No 1=Yes
	REAL IOM(MAXLAYER)			! IN(GET_SITE_SOIL):Inert organic C (kgC/ha/layer)
	INTEGER IORGMJ(0:MAXGROW,MAXORGM)	! IN(SETFILE):No.timesteps to manure applic.
	INTEGER IORGM(MAXORGM)		! IN(GETYRSET):No.timesteps to split manure application
	INTEGER IROCK				! IN(GETYRSET):No.layers: 1=50cm, 2=100cm , 3=150cm
	INTEGER IROCKJ				! IN(SETFILE):No.layers: 1=50cm, 2=100cm , 3=150cm
	INTEGER IROCKS(0:MAXCROP)	! IN(INIT_SUNDIAL_CROP)Maximum rooting depth (cm)
	INTEGER IRYEAR				! IN(GETYRSET):Current weather year number
	INTEGER IS_TS				! IN(GETYRSET):Crop counter
 	INTEGER ISERIES				! IN(GET_SITE_SOIL):Counter for soil series
	INTEGER ISOWN				! IN(GETYRSET):Timesteps from 01/01/01 to sowing date 
	INTEGER ISOWNJ(0:MAXGROW)	! IN(SETFILE):Timesteps from 01/01/01 to sowing date 
	INTEGER ISPINUP				! IN(READ_MODEL):Is N limitation spin-up used?
	 INTEGER ISPINUP_OFF		!   N limitation spin-up is not used
	 INTEGER ISPINUP_ON			!   N limitation spin-up is used
	 DATA ISPINUP_OFF,ISPINUP_ON /0,1/
	 DATA INSTABLECN,INPASSCN /1,2/
      INTEGER ISTART_TS			! IN(SETFILE):First simulated timestep since 01/01/01 
	INTEGER	ISTHARV				! IN(SETFILE):Timesteps from 01/01/01 to first harvest
	INTEGER ISTYR				! IN(SETFILE):First year in simulation (eg. 2001)
	INTEGER ISWAIT				! IN(CALL):Code to wait for key press (1) or not (0)
	INTEGER ITFUNC				! IN(READ_MODEL):Choice of temperature rate modifier 
								!   (ROTHC or HADLEY)
	 INTEGER ITFUNC_ROTHC		!   ROTHC temperature rate modifier
	 INTEGER ITFUNC_HADLEY		!   HADLEY temperature rate modifier
       DATA ITFUNC_ROTHC,ITFUNC_HADLEY /0,1/
	INTEGER IVOL(MAXFERT)		! IN(GETYRSET):N volatilised from fertiliser (kgN/ha)
      INTEGER IVOLJ(0:MAXGROW,MAXFERT)	! IN(SETFILE):Does fert.contain ammonium salts 
										! other than ammonium sulphate? (0=No; 1=Yes)
	INTEGER IVOLNF				! IN(GETFERT):N volatilised from fertiliser (kgN/ha)
	INTEGER IYEAR				! OUT(GETYRSET):Current growing season number
	INTEGER JCULT(MAXGROW) ! IN(SETFILE):Type of cultivation 
						   ! 0=Zero tillage (with or without mulching) may also be included
						   ! 1=Minimum tillage (or eco-tillage/conservation): Non inversion but depth of cultivation up to 5-10cm, normally this tillage makes a furrow for sowing of seeds and fertilizer application. This practice is applicable for both spring and cover/winter crops.
						   ! 2=Reduced tillage (or non-inversion): Cover both type i.e. non-inversion but depth of cultivation up to 15-20 cm and less number of tillage practices compared to conventional one, where applicable.
						   ! 3=Conventional  tillage (inversion): Depth of cultivation up to 20-30 cm, no limitation of tillage intensity and plough type.
	REAL JFERTJ(0:MAXGROW,MAXFERT,3)	! IN(SETFILE):Prop.NO3,NH4,urea in fertiliser
      INTEGER MODTYPE				! IN(SETFILE):Crop model type: 0=SUNDIAL, 1=MAGEC
	 INTEGER IS_SUNDIAL			!   Crop model type: 0=SUNDIAL
	 INTEGER IS_MAGEC			!   Crop model type: 1=MAGEC
	 DATA IS_SUNDIAL,IS_MAGEC /0,1/
	INTEGER JORGM(MAXORGM)		! IN(GETYRSET):Type of manure application
	INTEGER JORGMJ(0:MAXGROW,MAXORGM)	! IN(SETFILE):Type of manure application
	INTEGER JORGMNF				! IN(GETFYM):Type of manure application
	INTEGER JSTOP				! IN(INIT_SUNDIAL_CROP):Code to stop simulation 0=No 1=Yes
	INTEGER K_TS				! IN(GETWEATHER_AT_FC):No.timesteps since start of year(UNUSED?)
	INTEGER L_TSS(0:MAXCROP)	! IN(INIT_SUNDIAL_CROP):No.timesteps before harvest
	REAL LAI					! IN(INIT_SUNDIAL_CROP):Plant LAI at harvest
	REAL LAT					! IN(SETFILE): Latitude
	INTEGER LCROP				! IN(SETFILE):Previous crop type
	REAL LEACHDOC(MAXLAYER)		! IN/OUT:DOC leached from the layer (kgC/ha/layer/timestep)
	REAL LEACHDON(MAXLAYER)		! IN/OUT:DON leached from the layer (kgN/ha/layer/timestep)
	INTEGER LHARV				! IN(GETYRSET):No.timesteps from 01/01 to prev. harvest
	REAL LTA_AWC(12,MAXLAYER)	! IN(GETLTA):Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)	! IN(GETLTA):Average air temp this month (deg.C)
	INTEGER LUCODE				! OUT(GET_CTON_FOEREID):LU code for equilibrium run 
	INTEGER MCROP				! IN(RUN2_MAGEC_CROP):Crop type code number
	INTEGER MEASLAY				! OUT:Layer that soil is measured to
	INTEGER MEND				! IN(INIT_SUNDIAL_CROP):Number of timesteps from prev.harv.-harv.
	REAL MOBDOC(MAXLAYER)		! IN(INIT_SUNDIAL_SOILCN_NOPARS):Mobile DOC in each layer (kgC/ha)
	REAL MOBDON(MAXLAYER)		! IN(INIT_SUNDIAL_SOILCN_NOPARS):Mobile DON in each layer (kgN/ha)
	INTEGER N_STEPS				! IN(SETFILE):No.timesteps in a year
	INTEGER N_TS				! OUT(GETWEATHER):No.timesteps since start of year(UNUSED?)
	REAL NAVTMP					! IN(RUN2_MAGEC_CROP): Weather data used by MAGEC
	INTEGER NCULT				! IN(SETFILE): Number of cultivations
	INTEGER NDATE				! IN(INIT_SUNDIAL_CROP):Not used?
	INTEGER NFERT				! IN(GETYRSET):No.fertiliser applications to this crop
      INTEGER NFERTJ(0:MAXGROW)	! IN(SETFILE):No.fertiliser applications to crops
	INTEGER NGROW				! IN(SETFILE):No.growing seasons simulated
	REAL NITRIFN(MAXLAYER)		! IN(MICROBIAL_SUNDIAL_SOILCN):Nitrification from this layer (kg/ha/layer/timestep)
	INTEGER NORGM				! IN(GETYRSET):No.manure applications to this crop
	INTEGER NORGMJ(0:MAXGROW)	! IN(SETFILE):No.manure applications to crops
	REAL NRADS					! IN(RUN2_MAGEC_CROP): Weather data used by MAGEC		
      INTEGER NRESJ(0:MAXGROW)	! IN(SETFILE):Straw incorporated? 0=No 1=Yes
	INTEGER NSOIL				! IN(GETYRSET):Soil code number
	INTEGER NSOILJ				! IN(SETFILE):Soil code number
	INTEGER NSOMLAY(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN(GET_SITE_SOIL):Number of SOM layers
	INTEGER NSOW				! IN(INIT_SUNDIAL_CROP):Timesteps from prev.harv. to sowing date
	INTEGER NXYEARS				! IN(SETFILE):No.growing seasons simulated
	INTEGER NYEARS				! IN(SETFILE):No.weather years included in the sim.
	REAL OCROPN					! IN(GETYRSET):N uptake of crop (kgN/ha) 
	REAL OCROPNJ(0:MAXGROW)		! IN(SETFILE):N uptake of crop (kgN/ha) 
	REAL ORGC					! IN(INIT_SUNDIAL_CROP):Total org.C input (kgC/ha) 
	REAL ORGMA(MAXORGM)			! IN(GETYRSET):Amount of manure applied (t/ha)
      REAL ORGMJ(0:MAXGROW,MAXORGM)	! IN(SETFILE):Amount of manure applied (t/ha)
	REAL ORGMANF				! IN(GETFYM):Amount of manure applied (t/ha)
	REAL ORGN					! IN(INIT_SUNDIAL_CROP):Total org.N input (kgN/ha) 
	REAL PE						! IN(RUN2_MAGEC_CROP): Weather data used by MAGEC
	REAL PENMD					! IN(RUN2_MAGEC_CROP): Weather data used by MAGEC
	REAL PENMRS					! IN(RUN2_MAGEC_CROP): Weather data used by MAGEC
	INTEGER PH_MODEL			! IN(READ_MODEL):How is pH calculated?
	 INTEGER PH_PARAM			! pH set to neutral or passed from parameter file
       INTEGER PH_STATIC			! pH is read in from input
								!   file and does not change
	 INTEGER PH_FROM_VSD		! pH is calculated using VSD (a very 
								!   simple dynamic version of the MAGIC model
								!   by Ed Rowe & Chris Evans, CEH, Bangor)
	 DATA PH_PARAM,PH_STATIC,PH_FROM_VSD /0,1,2/			 
	REAL PLANTBIOMASS			! IN(INIT_SUNDIAL_CROP):Plant biomass at harvest [ kg ha-1]
	REAL PLANTHEIGHT			! IN(INIT_SUNDIAL_CROP):Plant height at harvest [cm]
	REAL PLANTUP(MAXLAYER)		! IN(RUN2_GIS_CROP)/OUT(TEST1_RES):Plant uptake from each soil layer
	REAL PNREQ(MAXLU)			! IN(GET_SITE_SOIL):Amount of N required by plant each yr
	                            !        (kgN/ha/year)
	REAL PI_ANN					! IN(GET_SITE_SOIL):Annual plant input (kgC/ha/year)
      REAL PI_C(MAXLAYER)			! IN(PARTITION_PLANTINPUT): Plant input C to soil (kgC/ha/step)
	REAL PI_CEQ_MON(12,MAXLAYER)! IN(GET_PLANT_DIST): Equilibrium plant C input each month
								!         in each layer (kgC/ha/month/layer)
      REAL PI_N(MAXLAYER)			! IN(PARTITION_PLANTINPUT): Plant input N to soil (kgC/ha/step)
      REAL PI_N15(MAXLAYER)		! IN(PARTITION_PLANTINPUT): Plant input N15 to soil(kgC/ha/step)
								!(1=arable,2=grassland,3=forestry,4=semi-natural)
	REAL PRECIP					! IN(GETWEATHER):original value of rain needs to be passed on from file for MCMODE
	REAL PREYLD					! IN(SETFILE):Yield of previous crop (t/ha)
	REAL RAIN					! IN(GETWEATHER_AT_FC): Rainfall (mm/timestep)
	REAL RDD					! IN(GETWEATHER_AT_FC): Weather data used by MAGEC
	REAL REFIL(MAXLAYER)		! IN(DRAIN_SUNDIAL_WATER):Water needed to fill to FC (mm/lay)
	REAL RLWNS					! IN(RUN2_MAGEC_CROP): Weather data used by MAGEC
	REAL RNIN					! IN(RUN1_SUNDIAL_CROP):Litter N input in this timestep (kgN/ha)
	REAL RNIN15					! IN(RUN1_SUNDIAL_CROP):Litter N15 input in timestep (kgN15/ha)
	REAL ROOT					! IN(INIT_SUNDIAL_WATER):Rooting depth according to restriction (cm)
	INTEGER ROOTLAYER			! IN(RUN2_SUNDIAL_CROP):Rooting layer 
	REAL RORGN					! IN(INIT_SUNDIAL_CROP):Total org.N input (kgN/ha) DONT PASS?
	REAL RPMCARB(MAXLAYER)		! IN(MICROBIAL_SUNDIAL_SOILCN):C in RPM at start of timestep
	REAL RPMCARB0(MAXLAYER)		! IN(INIT_SUNDIAL_SOILCN_NOPARS):C in resistant PM (kgC/ha/layer)
      REAL RPMCTON(MAXLAYER)		! IN(GET_CTON_FOEREID):RPM C:N ratio 
	REAL RPMNIT0(MAXLAYER)		! IN(INIT_SUNDIAL_SOILCN_NOPARS):N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN(INIT_SUNDIAL_SOILCN_NOPARS):N15 in resistant PM (kgN15/ha/layer)
	REAL RRG(0:MAXCROP)			! IN(INIT_SUNDIAL_CROP):Rate of root growth (cm/week)
	REAL RSTART					! IN(STARTRES):Initial N in debris pool (kgN/ha)
	REAL RSTART15				! IN(STARTRES):Initial N in debris pool (kgN15/ha)
	REAL SAND(MAXLAYER)			! IN(GET_SITE_SOIL):Sand content of the layer (%)
	REAL SATWATCONT(MAXLAYER)	! IN(INIT_SUNDIAL_WATER):Total water content at saturation (mm/layer)
	REAL SECONDS				! IN(SETFILE):Number of seconds in one timestep
	REAL SEED(0:MAXCROP)				! IN(RUN2_MAGEC_CROP): Data used by MAGEC
	REAL SEEDIN					! IN(RUN2_MAGEC_CROP): Data used by MAGEC
	REAL SEEDN_S				! IN(RUN2_MAGEC_CROP): Data used by MAGEC
	REAL SILT(MAXLAYER)			! IN(GET_SITE_SOIL):Silt content of the layer (%)
	CHARACTER SITE*80			! IN(OPENCHAN): Name of file containing soil characteristics at site 
	REAL SLEA15(MAXLAYER)		! IN(PHYSICAL_SUNDIAL_SOILCN):Nitrate-N15 leached (kgN15/ha)
	REAL SLEACH(MAXLAYER)		! IN(PHYSICAL_SUNDIAL_SOILCN):Nitrate-N leached (kgN/ha)
	REAL SLOPES					! IN(RUN2_MAGEC_CROP): Weather data used by MAGEC
	REAL SOIL15(MAXLAYER)		! IN(INIT_SUNDIAL_SOILCN_NOPARS):Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILFC(MAXLAYER)		! IN(INIT_SUNDIAL_WATER): soil water content at FC (mm/layer)
	REAL SOILN(MAXLAYER)		! IN(INIT_SUNDIAL_SOILCN_NOPARS):Soil nitrate-N (kgN/ha/layer)
	REAL SOILPH(MAXLAYER)		! IN(GET_SITE_SOIL):pH of soil in this layer
	REAL SOILTEMP(MAXLAYER)		! IN(GETWEATHER_AT_FC): Soil temperature (deg.C/timestep)
	REAL SOILW(MAXLAYER)		! IN(GET_AVAILWAT_FROM_TOTWAT):Available water (mm/layer)
	REAL SOMDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN(GET_SITE_SOIL):depth of soil organic matter layers
	REAL SORGC					! IN(INIT_SUNDIAL_CROP):Total org.N15 input (kgN/ha) DONT PASS?
	REAL SORGN					! IN(INIT_SUNDIAL_CROP):Total org.N15 input (kgN/ha) DONT PASS?
	REAL SORGN15				! IN(INIT_SUNDIAL_CROP):Total org.N15 input (kgN/ha) DONT PASS?
	INTEGER SPARMODEL			! IN(READ_MODEL):Soil parameter model (from file or calc)
	 INTEGER SPARFILE			!   Soil parameters read in from file
	 INTEGER SPARCALC			!   Soil parameters calculated from TOC etc
	 DATA SPARFILE,SPARCALC /1,2/
	REAL SRNIN					! IN(RUN2_SUNDIAL_CROP):Total litter N input (kgN/ha)
	INTEGER SUM_TS				! IN(GETWEATHER)/OUT(SETWTD):Total number of timesteps passed
	REAL SVPS					! IN(RUN2_MAGEC_CROP): Weather data used by MAGEC
	REAL SX						! IN(SETRES):Total org.N15 input (kgN/ha) DONT PASS?
	REAL SXORGC					! IN(INIT_SUNDIAL_CROP):Total org.N15 input (kgN/ha) DONT PASS?
	REAL SXORGN					! IN(INIT_SUNDIAL_CROP):Total org.N15 input (kgN/ha) DONT PASS?
	REAL SXORGN15				! IN(INIT_SUNDIAL_CROP):Total org.N15 input (kgN/ha) DONT PASS?
	REAL T1(0:MAXCROP)			! IN(INIT_SUNDIAL_CROP):Crop parameter for N uptake
	REAL T15					! IN(MICROBIAL_SUNDIAL_SOILCN):Net Mineralised N (kgN/ha/layer)
	REAL TACTOT					! IN(RUN2_SUNDIAL_CROP):Total crop N uptake (kgN/ha)
	REAL TAVS					! IN(RUN2_MAGEC_CROP): Weather data used by MAGEC
	REAL TAM(MAXLAYER)			! IN(ADD_SUNDIAL_SOILCN):15N/N in the ammonium in the layer
	REAL TC						! IN(RUN2_SUNDIAL_CROP):Total org.N15 input (kgN/ha) DONT PASS?
	REAL TCINP					! IN(RUN1_SUNDIAL_CROP):Total org.C input (kgC/ha) DONT PASS?
	REAL TDNIT					! IN(MICROBIAL_SUNDIAL_SOILCN):Net Mineralised N (kgN/ha/layer)
	REAL TFERT(MAXFERT,3)		! IN(GETYRSET):Prop.NO3,NH4,urea in fertiliser
	REAL TFERTNF(3)				! IN(GETYRSET):Prop.NO3,NH4,urea in fertiliser
	REAL TFYMC					! IN(ADD_SUNDIAL_SOILCN):Total FYM C input (kgC/ha) DONT PASS?
	REAL THISFERT				! IN(SETWEEKVARS):N input by fertiliser (kgN/ha)
	REAL THISFERT15				! IN(SETWEEKVARS):N15 input by fertiliser (kgN15/ha)
	REAL THISSOMDEPTH			! OUT(TEST1_RES):Depth of maximum SOM layer
	REAL TIN(MAXLAYER)			! IN:15N/N in the nitrate in the layer
	REAL TMAX					! IN(RUN2_MAGEC_CROP): Weather data used by MAGEC
	REAL TMIN					! IN(RUN2_MAGEC_CROP): Weather data used by MAGEC
	REAL TMMN					! IN(GETWEATHER_AT_FC): Weather data used by MAGEC
	REAL TMMX					! IN(GETWEATHER_AT_FC): Weather data used by MAGEC
	REAL TOC(MAXLAYER)			! IN(GET_SITE_SOIL):Total organic C (kgC/ha/layer)
      REAL TORGC					! IN(INIT_SUNDIAL_CROP):Total org.C input (kgC/ha) DONT PASS?
	REAL TORGN					! IN(INIT_SUNDIAL_CROP):Total org.N input (kgN/ha) DONT PASS?
	REAL TOTAST					! IN(STARTRES):Total ammonium in soil (kgN/ha)
	REAL TOTNST					! IN(STARTRES):Total nitrate in soil (kgN/ha)
	REAL TOTAST15				! IN(STARTRES):Total ammonium in soil (kgN15/ha)
	REAL TOTNST15				! IN(STARTRES):Total nitrate in soil (kgN15/ha)
	REAL TOTPIC(MAXLU)			! IN(GET_SITE_SOIL):Plant C input calculated using 
								!	  ROTHC equilibrium run for each LU (kgC/ha/yr)
	REAL TRATEM					! IN(MICROBIAL_SUNDIAL_SOILCN):Temperature modifier
	REAL TRNIN					! IN(RUN1_SUNDIAL_CROP):Total litter N input (kgN/ha)
	REAL TXORGC					! IN(INIT_SUNDIAL_CROP):Total org.C input (kgC/ha) DONT PASS?
	REAL TXORGN					! IN(INIT_SUNDIAL_CROP):Total org.N input (kgN/ha) DONT PASS?
	REAL UT(3,0:MAXCROP)		! IN(GETYRSET):Crop parameter used to estimate yield
	REAL VIGOUR					! OUT: Vigour of cultivation (0-1) Determines the proportion of humus released to biomass,DPM and RPM pools
	REAL VOLAT					! IN(SETWEEKVARS):N lost by volatilisation (kgN/ha)
	REAL VOLAT15				! IN(SETWEEKVARS):N15 lost by volatilisation (kgN15/ha)
	REAL VOLATN(MAXLAYER)		! IN(MICROBIAL_SUNDIAL_SOILCN):Volatilisation from this layer (kg/ha/layer/timestep)
	REAL VP						! IN(GETWEATHER_AT_FC): Weather data used by MAGEC
      REAL WATCONT(MAXLAYER)		! IN(MICROBIAL_SUNDIAL_SOILCN): Total water content of the soil
								!		layer (mm / layer)
	INTEGER WATERMODEL			! IN(READ_MODEL):Water model
	 INTEGER SUNDIAL_WATER		!   SUNDIAL water model
	 INTEGER SWAT_WATER			!	SWAT water model
	 DATA SUNDIAL_WATER,SWAT_WATER /0,1/
	REAL WDF					! IN(RUN2_MAGEC_CROP): Weather data used by MAGEC
	REAL WILTPOINT(MAXLAYER)	! IN(MICROBIAL_SUNDIAL_SOILCN): Total water conent of the soil layer (mm/layer)
	REAL WLEACH					! IN(SETWEEKVARS):N leached (kgN/ha)
	REAL WLEACH15				! IN(SETWEEKVARS):N15 leached (kgN15/ha)
	REAL WMAX(MAXLAYER)			! IN(SETWEEKVARS):Available water at field cap. (mm/layer)
	REAL WN,PIANN,TEMP			! IN(GETWEATHER_AT_FC): Weather data used by MAGEC
	REAL WR						! IN(RUN1_SUNDIAL_CROP):Root requirement (kgN/ha)
	REAL WRATEDM				! IN(MICROBIAL_SUNDIAL_SOILCN):Moisture modifier
	REAL WSAT(MAXLAYER)			! IN(SETWEEKVARS):Available water at saturation (mm/layer)
	REAL WTABLE					! IN(SETFILE):Water table depth in cm
	REAL XORGC					! IN(INIT_SUNDIAL_CROP):Total org.C input (kgC/ha) 
	REAL XORGN					! IN(INIT_SUNDIAL_CROP):Total org.N input (kgN/ha) 
	REAL YLD					! IN(INIT_SUNDIAL_CROP):Yield (t/ha) DONT PASS?
	REAL Z						! IN(GETWEATHER_AT_FC_SWAT):elevation above sea level [m]
C								  in the canopy when the canopy is fully developed (mm H2O)
	LOGICAL FULL_OUTPUT         ! Flags whether to write detailed output files
	LOGICAL SUMMARY_OUTPUT      ! Flags whether to write the summary output file
C
C Data Statements
C
      DATA IYEAR/0/
C
C************************************************
C Set model type FOR SITE SPECIFIC RUNS here!!!!
C 1. Denitrification model, DMODEL: DNEMIS = NEMIS model; 
C    
	DMODEL=DNEMIS
C
C 2. Initialisation of SOM pools, ICMODEL and INMODEL
C
	INMODEL=INSTABLECN
      ICMODEL=ICROTHCEQ
      ISPINUP=ISPINUP_OFF
C
C 3. DOC model, DOCMODEL: DOC_OFF	= DOC model off; DOC_ON = DOC model on
C
	DOCMODEL=DOC_ON
C
C 4. CN model, CNMODEL: CNFOEREID = C:N ratio obtained by method of Foereid; C:N ratio obtained by method of MAGEC
C
	CNMODEL=CNMAGEC
C
C 5. Soil parameter model, SPARMODEL: SPARFILE = Soil parameters read in from file; SPARCALC = Soil parameters calculated from TOC etc
C
	SPARMODEL=SPARFILE
C
C 6. Type of equilibrium run: EQNPP = Model initialised using plant inputs from measured NPP;
C							EQTOC = Model initialised using measured TOC;
C							EQNPPTOC = Model initialised using both plant inputs.
C
	EQMODEL=EQHILLIER
	EQMODEL=EQJONES
C
C 7. Calculation of moisture rate modifiers, IMFUNC_ROTHC=0, IMFUNC_HADLEY=1
C    Calculation of temp.rate modifiers, ITFUNC_ROTHC=0, ITFUNC_HADLEY=1
C
      IMFUNC=IMFUNC_ROTHC
      ITFUNC=ITFUNC_ROTHC
C
C 8. CH4 model, CH4MODEL: CH4_OFF	= CH4 model off; CH4_RICHARDS = Richards CH4 model on' CH4_AITKENHEAD = Aitkenhead CH4 model on
C
      CH4MODEL=CH4_RICHARDS  
C
C 9. Water model, WATERMODEL: SUNDIAL_WATER = original Sundial water routines
C							SWAT_WATER = water rountines from SWAT implemented
C
C This is still in development and currently only works for Tropical maize!
C
	WATERMODEL=SUNDIAL_WATER
C
C 10. Use equilibrium run of ECOSSE to initialise or not
C EC_EQRUN_ON = yes, EC_EQRUN_OFF = no
C
	EC_EQRUN=EC_EQRUN_OFF
C
C 11. Read in water table depth? 0=no, 1=yes
C
      ENTER_WTD=0
C
C Land use type 
C
	LUCODE=1			! Note: LUCODE set to arable because this 
						! version of the model is arable crops only
	ILU=1
C
C
C Output options
C
	FULL_OUTPUT=.TRUE.
	SUMMARY_OUTPUT=.TRUE.
C
C************************************************
C Set time factors
C
      IYEAR=0
	C_TS=0.0
      N_TS=0
C 
C 1. Open Channels for file input/output
C
	CALL OPENCHAN(SITE,ENTER_WTD)
	CALL TEST1_OPENCHAN(WATERMODEL,FULL_OUTPUT,SUMMARY_OUTPUT)
C
C 2. Read in mode of model use
C
      CALL READ_MODEL(DMODEL,ICMODEL,INMODEL,CNMODEL,DOCMODEL,SPARMODEL,
     &  		        EQMODEL,IMFUNC,ITFUNC,CH4MODEL,EC_EQRUN)
      CALL TELL_MODEL(DMODEL,ICMODEL,INMODEL,CNMODEL,DOCMODEL,SPARMODEL,
     &  		        EQMODEL,IMFUNC,ITFUNC,CH4MODEL,ISPINUP,ISWAIT)
C
C 3. Get filenames and simulation inputs
C
      CALL SETFILE(ATM,IDATEFC,NSOILJ,IDRAINJ,IROCKJ,LCROP,
     &                   PREYLD,SECONDS,MODTYPE,N_STEPS,
     &                   IAWCJ,NYEARS,
     &                   ISTHARV,ISTYR,ISTART_TS,ILAST_TS,FIXEND,LAT,
     &                   FIXFILE,NGROW,NXYEARS,ICROPJ,ISOWNJ,OCROPNJ,
     &                   IHARVJ,EXYLDJ,NRESJ,NFERTJ,NORGMJ,FERTJ,IFERTJ,
     &                   JFERTJ,IVOLJ,ILABJ,ORGMJ,IORGMJ,JORGMJ,
     &                   IOLABJ,WTABLE,
     &                   NCULT,ICULT,JCULT,
C Variable required for SWAT water
     &					FIXFILEAVP,WATERMODEL,SWAT_WATER,Z)
C
C
C 4. Read in TOC, IOM, CLAY, SILT, SAND, BULKDENS, NPP
C
      MEASLAY=MAXLAYER1
      IF(SPARMODEL.EQ.SPARCALC)THEN
	  CALL GET_SITE_SOIL(AVEPET,AVERAIN,AVETEMP,BULKDENS,CACCUM,
     &                     CLAY,DOMSOILISIMP,DOMSOILIMPDEPTH,EC_EQRUN,
     &                     EQMODEL,ICMODEL,ILU,IOM,ISERIES,NSOMLAY,
     &                     PI_ANN,PNREQ,SAND,SILT,SITE,SOILPH,SOMDEPTH,
     &                     TOC,TOTPIC)
        MEASLAY=SOMDEPTH(ISERIES,ILU,NSOMLAY(ISERIES,ILU))
    	  MEASLAY=MEASLAY*MAXLAYER1/MAXDEPTH
	  PH_MODEL=PH_STATIC
	ENDIF
C
C For each growing season... 
C        
      DO 100 IYEAR=0,NXYEARS
C
C Get setup information for this year
C
	  CALL GETYRSET(IYEAR,NXYEARS,NSOIL,NSOILJ,IAWC,IAWCJ,
     &                IROCK,IROCKJ,INSTRAW,NRESJ,
     &                ICROP,ICROPJ,ISOWN,ISOWNJ,OCROPN,OCROPNJ,
     &                IHARV,IHARVJ,EXYLD,EXYLDJ,NFERT,NFERTJ,
     &                NORGM,NORGMJ,FERT,FERTJ,IFERT,IFERTJ,
     &                ILAB,ILABJ,IVOL,IVOLJ,TFERT,JFERTJ,
     &                ORGMA,ORGMJ,IORGM,IORGMJ,IOLAB,IOLABJ,
     &                JORGM,JORGMJ,LCROP,CFACT,UT,
     &                LHARV,SECONDS,IS_TS,IRYEAR,IFILE,ISTHARV)
C
C Initialise current crop
C 
        CALL INIT_SUNDIAL_CROP(IYEAR, LCROP, ICROP, 
     &                             INSTRAW, IEND, ISTHARV,   
     &                             ISOWN, MEND, NSOW, JSTOP, L_TSS,
     &                             YLD, PREYLD, EXYLD,  
     &                             SORGC,SORGN, SXORGC, SXORGN,  
     &                             HZ1, TORGC, TORGN,TXORGC, TXORGN, 
     &                             OCROPN, RORGN,DDAYS,
     &                             ORGC, ORGN, XORGC, XORGN,
     &                             NXYEARS,NDATE,FIXEND,LHARV,IHARV,
     &                             NORGM,IORGM,NFERT,IFERT,IANTHES,
     &                             SECONDS,
C introduced for SWAT water routines
     &				             PLANTBIOMASS,PLANTHEIGHT,LAI,C1,T1,
     &							 CANMAX,CONVER_F,RRG,IROCKS)
C
C Initialise soil water and C&N in year 0
C
        IF(IYEAR.EQ.0)THEN
C
C Required for SWAT water
C
		IF(WATERMODEL.EQ.SWAT_WATER)THEN
		 CALL GETWEATHER_AT_FC_SWAT(IFILE,IDATEFC,N_STEPS,IRYEAR,K_TS,
     &                          RAIN,EVAP,AIRTEMP,RDD,TMMN,TMMX,VP,WN,
     &                          FIXFILE,FIXFILEAVP,MODTYPE,SOILTEMP,AVP)
		 CALL AVEMET(AVERAIN,AVEPET,AVETEMP) ! 
           CALL INIT_SUNDIAL_WATER_SWAT(WATCONT,WILTPOINT,ROOT,IROCK,
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
     &						      CURRENTPLBIOM,CURRENTPLHEI,DDAYS,
     &							  z,AVP,LAT,LCROP,RRG,
     &							  IROCKS,SOILFC,WMAX,SOILW,SECONDS)

           CALL GET_AVAILWAT_FROM_TOTWAT(FIELDCAP,WATCONT,SATWATCONT,
     &                                   WILTPOINT,WMAX,SOILW,WSAT)		 
           CALL GETLTA_SWAT(LTA_AWC,LTA_TEMP,WMAX,WSAT,FLOWPROP,
     &                  AVERAIN,AVEPET,AVETEMP)

          ELSE
		 CALL GETWEATHER_AT_FC(IFILE,IDATEFC,N_STEPS,IRYEAR,K_TS,
     &                          RAIN,EVAP,AIRTEMP,RDD,TMMN,TMMX,VP,WN,
     &                          FIXFILE,MODTYPE,SOILTEMP)
		 CALL AVEMET(AVERAIN,AVEPET,AVETEMP)
          CALL MIAMI_DYCE(LUCODE,sum(AVETEMP)/12,sum(AVERAIN),PI_ANN)
           CALL INIT_SUNDIAL_WATER(SOILW,ROOT,IROCK,SUM_TS,IDATEFC,
     &                              IANTHES,IRYEAR,IFILE,WMAX,K_TS,RAIN,
     &                              EVAP,AIRTEMP,RDD,TMMN,TMMX,VP,WN,
     &                              IAWC,NSOIL,MODTYPE,SPARMODEL,
     &                              TOC,CLAY,SAND,SILT,BULKDENS,WSAT,
     &                              WTABLE,
     &                              FLOWPROP,AVERAIN,AVEPET,AVETEMP,
     &                              IYEAR,ISTHARV,ISOWNJ,N_STEPS,
     &                              WILTPOINT, SATWATCONT,SECONDS)
		 CALL GETLTA(LTA_AWC,LTA_TEMP,WMAX,WSAT,FLOWPROP,
     &                  AVERAIN,AVEPET,AVETEMP)
          ENDIF
C
C
C
  	    CALL GET_SUNDIAL_SOILPARS(NSOIL,TOC,IOM,CRIT,ILU,
     &                            CLAY,BULKDENS,SPARMODEL)
	    CALL GET_PLANT_DIST(PI_ANN,PI_CEQ_MON,ILU)
          IF(INMODEL.EQ.INPASSCN)
     &      CALL GET_CTON_FOEREID(DPMCTON,RPMCTON,ILU)
c     &      CALL GET_CTON_FALLOON(DPMCTON,RPMCTON,ILU)
	    IF(ICMODEL.EQ.ICROTHCEQ)
     &      CALL SET_DPMRPMRATIO(ILU,DPM_RPM_RATIO)
	    CALL INIT_SUNDIAL_SOILCN_NOPARS(ICMODEL,INMODEL,DOCMODEL,
     &                               EQMODEL,SECONDS,
     &                               NSOIL,TOC,IOM,HZ1,CRIT,
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
          CALL STARTRES(IRYEAR,IK,RSTART,DPMNIT0,RPMNIT0,
     &                    BSTART,BNIT0,HNIT0,
     &                    RSTART15,DPMNLAB0,RPMNLAB0,
     &                    BSTART15,BNLAB0,HNLAB0,
     &                    TOTAST,AMMN,TOTAST15,AMMN15,
     &                    TOTNST,SOILN,TOTNST15,SOIL15)
C
C ... equilibrium run using ECOSSE
C
          IF(EC_EQRUN.EQ.EC_EQRUN_ON)THEN
	      CALL ECOSSE_EQRUN(AMMN,AMMN15,ATM,AVEPET,AVERAIN,AVETEMP,
     &              BCARB0,BNIT0,BNLAB0,CACCUM,CH4MODEL,CNMODEL,CTOT,
     &			  CLOSSX,CLOSSX15,CRIT,DMODEL,DOCMODEL,DPMCARB0,
     &			  DPMNIT0,DPMNLAB0,DRAINW,FERT,FIXN,FLOWPROP,
     &              FYMFERT,FYMFERT15,HCARB0,HNIT0,HNLAB0,
     &              IANTHES,IEND,IFERT,IHARV,ILAB,IMFUNC,ICMODEL,
     &              INMODEL,IOLAB,IOM,IORGM,ISERIES,ISOWN,ITFUNC,IVOL,
     &              JORGM,JSTOP,LHARV,ILU,MEND,NFERT,NORGM,
     &              NSOIL,NSOW,ORGMA,PNREQ,REFIL,RPMCARB0,
     &              RPMNIT0,RPMNLAB0,SOIL15,SOILN,SOILW,TAM,
     &              TFERT,THISFERT,THISFERT15,TIN,TOC,TOTPIC,VOLAT,
     &              VOLAT15,WLEACH,WLEACH15,WSAT,WMAX,SOILPH,NXYEARS,
     &              EQMODEL,CLAY,PI_CEQ_MON,DOMSOILISIMP,
     &              DOMSOILIMPDEPTH,SOMDEPTH,NSOMLAY,LTA_AWC,LTA_TEMP,
     &              PH_MODEL)
C
C Reinitialise current crop
C 
            CALL INIT_SUNDIAL_CROP(IYEAR, LCROP, ICROP, 
     &                             INSTRAW, IEND, ISTHARV,   
     &                             ISOWN, MEND, NSOW, JSTOP, L_TSS,
     &                             YLD, PREYLD, EXYLD,  
     &                             SORGC,SORGN, SXORGC, SXORGN,  
     &                             HZ1, TORGC, TORGN,TXORGC, TXORGN, 
     &                             OCROPN, RORGN,DDAYS,
     &                             ORGC, ORGN, XORGC, XORGN,
     &                             NXYEARS,NDATE,FIXEND,LHARV,IHARV,
     &                             NORGM,IORGM,NFERT,IFERT,IANTHES,
     &                             SECONDS,
     &				             PLANTBIOMASS,PLANTHEIGHT,LAI,C1,T1,
     &							 CANMAX,CONVER_F,RRG,IROCKS)
	    ENDIF
	  ENDIF
C
C For each timestep till end of growing season...
C
        DO 200 IK=1,MEND
C
C Initialize THISFERT=current weeks fertilizer addition
C            WLEACH=current weeks water leaching
C      
          CALL SETWEEKVARS(JSTOP,CLOSSX,CLOSSX15,FYMFERT,FYMFERT15,
     &                     FIXN,THISFERT,THISFERT15,WLEACH,WLEACH15,
     &                     VOLAT,VOLAT15)
C
C
C Get weather data
C
C	Required for SWAT water
C
		IF(WATERMODEL.EQ.SWAT_WATER)THEN
            CALL GETWEATHER_SWAT(N_TS,N_STEPS,IFILE,IRYEAR,IK,MEND,
     &                    SUM_TS,RAIN,EVAPW,AIRTEMP,RDD,TMMN,TMMX,
     &                    VP,WN,FIXFILE,FIXFILEAVP,MODTYPE,SOILTEMP,
C required for Monte carlo - need to pass on original value ofprecipiation
     &					AVP,PRECIP)
		ELSE
            CALL GETWEATHER(N_TS,N_STEPS,IFILE,IRYEAR,IK,MEND,SUM_TS,
     &                    RAIN,EVAPW,AIRTEMP,RDD,TMMN,TMMX,VP,WN,
     &                    FIXFILE,MODTYPE,SOILTEMP,
C required for Monte carlo - need to pass on original value ofprecipiation
     &						precip)
		ENDIF
C
C Read in water table depth
C
          CALL SETWTD(SUM_TS,ENTER_WTD,SOILW,WSAT)
C
C If current week is first in the growing season set results
C
          CALL SETRES(BALANCE,CO2,SX,SXORGN)
          IS_SUNDIAL=0
	    IS_MAGEC=1
C Work out what cultivations are hapenning this week
C	0=Zero tillage (with or without mulching) may also be included
C	1=Minimum tillage (or eco-tillage/conservation): Non inversion but depth of cultivation up to 5-10cm, normally this tillage makes a furrow for sowing of seeds and fertilizer application. This practice is applicable for both spring and cover/winter crops.
C	2=Reduced tillage (or non-inversion): Cover both type i.e. non-inversion but depth of cultivation up to 15-20 cm and less number of tillage practices compared to conventional one, where applicable.
C	3=Conventional  tillage (inversion): Depth of cultivation up to 20-30 cm, no limitation of tillage intensity and plough type.
C
          CULTIVATE=0
          DO 300 IL=1,NCULT
	      IF(IK.EQ.ICULT(IL)-LHARV)THEN
              IF(JCULT(IL).EQ.0)THEN
			  CULTDEPTH=0
	          INVERT=0
			  VIGOUR=0.0
			ELSEIF(JCULT(IL).EQ.1)THEN
			  CULTDEPTH=10
	          INVERT=0
			  VIGOUR=0.0
	        ELSEIF(JCULT(IL).EQ.2)THEN
	          CULTDEPTH=15
	          INVERT=0
	          VIGOUR=0.0
	        ELSEIF(JCULT(IL).EQ.3)THEN
	          CULTDEPTH=20
	          INVERT=1
	          VIGOUR=0.0
              ENDIF
   		    PRINT*,'Please enter the vigour of cultivation ',IL,
     &               ' on timestep',ICULT(IL)
	        READ(*,*)VIGOUR
		    CULTIVATE=1	        
	      ENDIF
300       CONTINUE
C
C SUNDIAL Crop Model...
C
          IF(MODTYPE.EQ.IS_SUNDIAL)THEN
C
C ...Calculate crop C and N returns and N offtake 
C
	      CALL RUN1_SUNDIAL_CROP(IYEAR,IK,MEND,IS_TS,JSTOP,
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
      
	    ELSEIF(MODTYPE.EQ.IS_MAGEC)THEN
C
C ...Calculate crop C and N returns 
C
	      IF((IYEAR.EQ.0).OR.((IYEAR.EQ.1).AND.(IK.EQ.1)))then
	        CALL RUN1_MAGEC_CROP(I_TS,IEND,IYEAR,IK,SECONDS,
     &                  CACTOT,TRNIN,RNIN,RNIN15,CATOT15,
     &                  ICROP,SXORGN,SXORGN15,C_TS,SXORGC,TCINP,
     &                  XORGC,XORGN,ORGC,ORGN,
     &                  ISOWN,IANTHES,NSOW,ICOVER,CULTIVATE,NCULT)
            ENDIF
   		ENDIF
C
C ...Add stuff to soil
C
          CALL GETFYM(IK,NORGM,IORGM,ORGMA,JORGM,IOLAB,
     &                ORGMANF,JORGMNF,IOLABNF) 
          CALL GETFERT(IK,NFERT,IFERT,FERT,TFERT,ILAB,IVOL,
     &                  FERTNF,TFERTNF,ILABNF,IVOLNF) 
          CALL ADD_SUNDIAL_SOILCN(
C INPUTS: time, environment, organic manure, fertiliser, atmospheric, 
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
C Cultivate soil
C
          IF(CULTIVATE.EQ.1)THEN 
      	  CALL CULTIV(LUCODE,
     &                  DPMCARB0,DPMNIT0,DPMNLAB0,
     &                  RPMCARB0,RPMNIT0,RPMNLAB0,
     &                  BCARB0,BNIT0,BNLAB0,
     &                  HCARB0,HNIT0,HNLAB0,
     &                  CULTDEPTH,VIGOUR,INVERT)
	    ENDIF
C
C Run SUNDIAL Microbial Soil C and N routines
C 
          CALL PARTITION_PLANTINPUT(C_TS,RNIN,RNIN15,PI_C,PI_N,PI_N15)
          IF(INMODEL.EQ.INPASSCN)
     &      CALL GET_CTON_FOEREID(DPMCTON,RPMCTON,LUCODE)
c     &      CALL GET_CTON_FALLOON(DPMCTON,RPMCTON,LUCODE)
	    CALL SET_DPMRPMRATIO(LUCODE,DRRAT)
          CALL MICROBIAL_SUNDIAL_SOILCN(DMODEL,DOCMODEL,CNMODEL,
     &                                 ITFUNC,IMFUNC,CH4MODEL,INMODEL,
     &                                 SECONDS,ICOVER,
     &                                 NSOIL,SOILTEMP,RAIN,
     &								 PI_C,PI_N,PI_N15,DRRAT,                                                                     
     &                                 WMAX,SOILW,WSAT,
     &                                 WILTPOINT,SATWATCONT,
     &                                 THISFERT,THISFERT15,
     &                                 CO2,VOLAT,VOLAT15,VOLATN,
     &                                 CH4,CH4_PROD,CH4_SOIL_OX,
     &	    					     CH4_ATMOS_OX,CH4_FLUX,
     &                                 DENITRIFN,NITRIFN,
     &                                 DENIT,DN15,GNN2O,GNNO,GPNN2O,	
     &                                 G15NN2O,G15NNO,G15PNN2O,		
     &                                 GDN2,GDN2O,G15DN2,G15DN2O,		
     &                                 DNIT,DNIT15,TDNIT,T15,
     &                                 WLEACH,WLEACH15,
     &                                 SOILN,SOIL15,AMMN,AMMN15,
     &                                 DPMCARB0,DPMNIT0,DPMNLAB0,
     &                                 RPMCARB0,RPMNIT0,RPMNLAB0,
     &                                 BCARB0,BNIT0,BNLAB0,BULKDENS,
     &                                 TOC,CH4TOAIR,HCARB0,HNIT0,HNLAB0,
     &                                 TIN,TAM,MOBDOC,MOBDON,CO2FROMDOC,
     &                                 ICFACTOR,DPMCTON,RPMCTON,
     &					             DPMCARB,RPMCARB,BCARB,HCARB,
     &					             WRATEDM,TRATEM,SOILPH,MEASLAY)
C
C Run SUNDIAL within season crop routines
C
          IF(MODTYPE.EQ.IS_SUNDIAL)THEN
            CALL RUN2_SUNDIAL_CROP(IYEAR,IK,NSOW,ICROP,JSTOP,MEND,
     &                             IS_TS,NSOIL,IL_TSS,AIRTEMP,DDAYS,
     &                             CACT,CACT15,ATM,ATM15,CACTOT,TC,
     &                             L_TSS,TACTOT,SRNIN,TRNIN,
     &                             RNIN,OCROPN,CLOSSX,CLOSSX15,RORGN,
     &                             VOLAT,VOLAT15,CATOT15,CUPTN,SOILN, 
     &                             TIN,SOIL15,AMMN15,AMMN,TAM,CRIT,CTOT,
     &                             SEEDIN,XORGN,
C Variables required for SWAT water
     &				             PLANTBIOMASS,PLANTHEIGHT,LAI,
     &							 CURRENTPLBIOM,CURRENTPLHEI,CURRENTLAI,
     &							 ROOTLAYER,PLANTUP)
          ELSEIF(MODTYPE.EQ.IS_MAGEC)THEN
            CALL RUN2_MAGEC_CROP(CRIT,IYEAR,ICROP,IK,NSOW,MEND,C_TS,
     &                          RNIN,DTR,RDD,DVP,VP,TMAX,
     &                          TMMX,TMIN,TMMN,DAVTMP,NAVTMP,
     &                          DTRJM2,TAVS,BBRADS,SVPS,SLOPES,
     &                          RLWNS,NRADS,PENMRS,WDF,WN,PENMD,
     &                          PE,AE,SOILW,EVAPW,DDAYS,IDAG,
     &                          WMAX,N_TS,RAIN,
     &                          LAT,SOILN,AMMN,CRITMIN,IROCKJ,INSTRAW,
     &                          EVAPO,CACT,SEEDN_S,SEED,MCROP,NSOIL)
	    ENDIF
C
C Water routines
C
C Check if SWAT watermodel should be used
C
	    IF(WATERMODEL.EQ.SWAT_WATER)THEN
C
C Intercepted rainfall
C
		  CALL INTERCEPTION_SWAT_WATER(LAI,CURRENTLAI,CANMAX,RAIN,
     &								CANSTOR)
	    ENDIF		
C
C Leaching routines
C
C          CALL DRAIN_SUNDIAL_WATER(RAIN,WMAX,SOILW,DRAINW,REFIL,
C     &                             WSAT,FLOWPROP)
	    CALL DRAIN_SUNDIAL_WATER2(RAIN, WMAX, SOILW, WSAT, FLOWPROP,
     &                              DRAINW, REFIL)
C
C Evapotranspiration 
C
	    IF(WATERMODEL.EQ.SWAT_WATER)THEN
C
C Calculation of actual evapotranspiration (SWAT)
C
		  CALL EVAP_SWAT_WATER(AIRTEMP,Z,AVP,LAT,N_TS,CURRENTPLHEI,
     &						CURRENTLAI,CANSTOR,SOILW,FIELDCAP,
     &                        WMAX,WILTPOINT,ROOTLAYER,
     &                        CURRENTPLBIOM,EVAPW,WATCONT)
            CALL GET_AVAILWAT_FROM_TOTWAT(FIELDCAP,WATCONT,SATWATCONT,
     &                                   WILTPOINT,WMAX,SOILW,WSAT)		 
	    ELSE
C
C Evapotranspiration from input file
C
	      CALL EVAP_SUNDIAL_WATER(EVAPW,AIRTEMP,ICOVER,
     &                            WMAX,SOILW)
	    ENDIF
C
C End of water routines
C

C
C Run SUNDIAL physical routines
C
          CALL PHYSICAL_SUNDIAL_SOILCN(DOCMODEL,SECONDS,
     &                                 NSOIL,SOILTEMP,
     &                                 WMAX,SOILW,
     &                                 DRAINW,REFIL,
     &                                 SLEACH,SLEA15,
     &                                 WLEACH,WLEACH15,CONC,CONC15,
     &                                 SOILN,SOIL15,AMMN,AMMN15,
     &                                 TIN,TAM,MOBDOC,MOBDON,
     &                                 LEACHDOC,LEACHDON)
C
C Record the results at the end of the week
C
	    THISSOMDEPTH=MEASLAY*MAXDEPTH/MAXLAYER1
          IF(IYEAR.GT.0)
     &       CALL TEST1_RES(IYEAR,IRYEAR,NXYEARS,LHARV,
     &                   SECONDS,IK,N_TS,SUM_TS,FIXEND,NSOW,ISOWN,MEND,
     &                   SX,SXORGN15,SORGN,SORGN15,RNIN15,
     &                   CACT,CACT15,CACTOT,CATOT15,
     &                   CLOSSX,CLOSSX15,VOLAT,VOLAT15,
     &                   ATM,ATM15,DNIT,TDNIT,T15,
     &                   THISFERT,THISFERT15,WLEACH,WLEACH15,
     &                   DN15,                   
     &                   FYMFERT,FYMFERT15,
     &                   IANTHES,CONC,CONC15,SLEACH,
     &                   SOILN,AMMN,TOTAST,TOTNST,TOTAST15,TOTNST15,
     &                   DPMCARB0,RPMCARB0,BCARB0,HCARB0,
     &                   DPMNIT0,RPMNIT0,HNIT0,BNIT0,DENIT,
     &                   BSTART,RSTART,BSTART15,RSTART15,RNIN,
     &                   ICROP,SEEDIN,GNN2O,GNNO,GPNN2O,G15NN2O,G15NNO,
     &                   G15PNN2O,GDN2,GDN2O,G15DN2,G15DN2O,
     &                   SOILW,MOBDOC,MOBDON,CO2FROMDOC,CO2,CH4,TFYMC,
     &                   CH4MODEL,CH4_PROD,CH4_SOIL_OX,
     &                   CH4_ATMOS_OX,CH4_FLUX,
     &                   TCINP,IOM,CH4TOAIR,WMAX,WSAT,C_TS,
     &                   DENITRIFN,NITRIFN,VOLATN,
     &                   LEACHDOC,LEACHDON,
c     &				   SOMDEPTH(ISERIES,ILU,NSOMLAY(ISERIES,ILU)),
     &				   THISSOMDEPTH,
     &                   PLANTUP,FULL_OUTPUT,SUMMARY_OUTPUT)
          CALL PUTRES(IYEAR,IK,SOILN,AMMN,SOILW,MOBDOC,CO2FROMDOC,
     &				WILTPOINT,WATCONT,WATERMODEL,
     &				CURRENTPLBIOM,CURRENTPLHEI,CURRENTLAI)
c          CALL TS_RES1(IYEAR,IRYEAR,NXYEARS,SECONDS, 
c     &                 IK,N_TS,SUM_TS,FIXEND,NSOW,ISOWN,MEND,
c     &                 SX,SXORGN15,SORGN,SORGN15,RNIN15,
c     &                 CACT,CACT15,CACTOT,CATOT15,
c     &                 CLOSSX,CLOSSX15,VOLAT,VOLAT15,
c     &                 ATM,ATM15,TDNIT,T15,
c     &                 THISFERT,THISFERT15,
c     &                 WLEACH,WLEACH,
c     &                 DN15,                   
c     &                 FYMFERT,FYMFERT15,
c     &                 IANTHES,CONC,CONC15,SLEACH,
c     &                 SOILN,AMMN,TOTAST,TOTNST,TOTAST15,TOTNST15,
c     &                 DPMCARB0,RPMCARB0,BCARB0,HCARB0,
c     &                 DPMNIT0,RPMNIT0,HNIT0,BNIT0,DENIT,
c     &                 BSTART,RSTART,BSTART15,RSTART15,RNIN,
c     &                 ICROP,SEEDIN,GNN2O,GNNO,GPNN2O,G15NN2O,G15NNO,
c     &                 G15PNN2O,GDN2,GDN2O,G15DN2,G15DN2O)
C
C Record the carbon balance
C
          IF (FULL_OUTPUT) THEN
            CALL CBALANCE(IYEAR,CO2,CH4,CH4TOAIR,IK,TFYMC,N_TS,TCINP,
     &                 HCARB0,BCARB0,DPMCARB0,RPMCARB0,IOM,PIANN,C_TS)
	    ENDIF
C
C If at the fixed end jump out of simulation
C
          IF(SUM_TS.EQ.FIXEND) GOTO 7200
C
C Go back and do calculation for the next week in the growing season
C
200     CONTINUE
C
C Go back and set up the next crop
C
100   CONTINUE
C
C Fixed End
C
7200  CONTINUE
C 
C Close Channels for file input/output
C
	CALL CLOSECHAN()
C
C Close the fileS
C
      CLOSE(3)
C
C Loop back for multiple runs
C
      CLOSE(4)
C
C Record successful completion of simulation in noerror.msg
C
      CLOSE(15)
      OPEN(15,FILE='NOERROR.MSG',STATUS='UNKNOWN')
      WRITE(15,*)'SIMULATION COMPLETED!'
      CLOSE(15)
	CLOSE(21)

      WRITE(*,*)'****************************************************'
      WRITE(*,*)'SIMULATION SUCCESSFULLY COMPLETED!!!'
	WRITE(*,*)
	WRITE(*,*)'...check N results in N_Balance.out'		
	WRITE(*,*)
      WRITE(*,*)'****************************************************'
c      IF(ISWAIT)THEN
c	  WRITE(*,*)'                       ....Press any key to continue'
c	  READ(*,*)
c	ENDIF
      END
C*******************************************************************
C INTERNAL ROUTINES
C
C-------------------------------------------------------------------
C
	SUBROUTINE ADDYEAR(IFERT,IHARV,IORGM,ISOWN,IYEAR,LHARV,SECONDS)
C
C to add year cro management operations
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXFERT				! Max.no.of fertiliser applications allowed 
	PARAMETER (MAXFERT=5)
	INTEGER MAXORGM				! Max.no.of organic manure applications allowed
	PARAMETER (MAXORGM=52)

      INTEGER NF					! Fertiliser / organic manure counter
	INTEGER STEPSIN1YR			! Time steps in 1 year
C
C Variables passed to/from this subroutine
C
	INTEGER IFERT(MAXFERT)		! IN(CALL)/OUT(CALL):No.timesteps to split fert.application
	INTEGER IHARV				! IN(CALL)/OUT(CALL):Timesteps from 01/01/01 to harvest date 
	INTEGER IORGM(MAXORGM)		! IN(CALL)/OUT(CALL):No.timesteps to split manure application
	INTEGER ISOWN				! IN(CALL)/OUT(CALL):Timesteps from 01/01/01 to sowing date 
	INTEGER IYEAR				! IN(CALL)/OUT(CALL):Current growing season number
	INTEGER LHARV				! IN(CALL)/OUT(CALL):No.timesteps from 01/01 to prev. harvest
	REAL SECONDS				! IN(CALL):Number of seconds in one timestep
C
C Work out number of timesteps in each year
C
      STEPSIN1YR=(365.25*24*60*60)/SECONDS
C
C Add the number of timesteps in 1 year to each variable
C
      DO 100 NF=1,MAXFERT
	  IFERT(NF)=IFERT(NF)+STEPSIN1YR
100   CONTINUE
	IHARV=IHARV+STEPSIN1YR
      DO 200 NF=1,MAXORGM
	  IORGM(NF)=IORGM(NF)+STEPSIN1YR
200   CONTINUE
	ISOWN=ISOWN+STEPSIN1YR
	IYEAR=IYEAR+1
	LHARV=LHARV+STEPSIN1YR
C
C leave ADDYEAR
C
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE AVEMET(AVERAIN,AVEPET,AVETEMP)
C
C To average weather data 1900-1930 for equilibrium run <1900
C
      IMPLICIT NONE
C
C Variables local to this routine
C
	INTEGER MONTH				! Month as read in from a file
	INTEGER M					! Month counter
C
C Variables passed to/from this routine
C
      REAL AVERAIN(12),AVEPET(12),AVETEMP(12)
C
C Read in weather data from averaged weather data file
C
      OPEN(50,FILE='AVEMET.DAT',STATUS='UNKNOWN',ERR=111)
	GOTO 101
111   CONTINUE
      PRINT*,'Error in file AVEMET.DAT'
	STOP
101   CONTINUE
	DO 100 M=1,12
	  READ(50,*)MONTH,AVERAIN(M),AVEPET(M),AVETEMP(M)
100   CONTINUE
      CLOSE(50)
	END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CLOSECHAN()
C
C Subroutine to close channels opened for site specific calculations
C
      IMPLICIT NONE
      CLOSE(4)
      CLOSE(10)
      CLOSE(15)
	CLOSE(20)
	CLOSE(21)
	CLOSE(23)
	CLOSE(46)
      CLOSE(52)
      CLOSE(53)
      CLOSE(54)
	CLOSE(55)
	CLOSE(56)
	CLOSE(57)
	CLOSE(58)
	CLOSE(64)
	CLOSE(97)
	CLOSE(99)
	CLOSE(98)
	CLOSE(99)
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE GET_AVAILWAT_FROM_TOTWAT(FIELDCAP,WATCONT,SATWATCONT,
     &                                    WILTPOINT,WMAX,SOILW,WSAT)
C
C Subroutine to translate total soil water into available soil water 
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER IL					! Layer counter
C
C Variables passed to/from this subroutine
C
      REAL FIELDCAP(MAXLAYER)		! IN: Total water in the layer at 
	                            !		field capacity (mm / layer)
	REAL SATWATCONT(MAXLAYER)	! IN: Total water content at 
								!		saturation (mm / layer)
	REAL SOILW(MAXLAYER)		! OUT:Available water (mm/layer)
      REAL WATCONT(MAXLAYER)		! IN: Total water content of the soil
								!		layer (mm / layer)
	REAL WILTPOINT(MAXLAYER)	! IN: Total water conent of the soil
								!		layer (mm / layer)
	REAL WMAX(MAXLAYER)			! OUT:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! OUT:Available water at saturation (mm/layer)
C
C Translate total soil water into available soil water
C
      DO 100 IL=1,MAXLAYER
        WMAX(IL)=FIELDCAP(IL)-WILTPOINT(IL)
	  SOILW(IL)=WATCONT(IL)-WILTPOINT(IL)
	  WSAT(IL)=SATWATCONT(IL)-WILTPOINT(IL)
	  IF(SOILW(IL).GT.WSAT(IL))SOILW(IL)=WSAT(IL)
	  IF(SOILW(IL).LT.0)SOILW(IL)=0
	  IF(WMAX(IL).GT.WSAT(IL))WMAX(IL)=WSAT(IL)
	  IF(WMAX(IL).LT.0)WMAX(IL)=0
	  IF(WSAT(IL).LT.SOILW(IL))WSAT(IL)=SOILW(IL)
	  IF(WSAT(IL).LT.0)WSAT(IL)=0
100   CONTINUE
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GET_SITE_SOIL(AVEPET,AVERAIN,AVETEMP,BULKDENS,CACCUM,
     &                     CLAY,DOMSOILISIMP,DOMSOILIMPDEPTH,EC_EQRUN,
     &                     EQMODEL,ICMODEL,ILU,IOM,ISERIES,NSOMLAY,
     &                     PI_ANN,PNREQ,SAND,SILT,SITE,SOILPH,SOMDEPTH,
     &                     TOC,TOTPIC)
C
C To get soil characteristics at site
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/
	INTEGER MAXGROW				! Max.no.of growing seasons allowed
	PARAMETER (MAXGROW=300)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXSOMLAY			! Max.no.of soil layers entered			
	PARAMETER (MAXSOMLAY=10)

	INTEGER IL					! Local counter variable
C
C Variables passed to/from this subroutine
C
      REAL AVEPET(12)				! IN(SETFILE_SITE):Long term average PET (mm/month)
      REAL AVERAIN(12)			! IN(SETFILE_SITE):Long term average rainfall(mm/month)
      REAL AVETEMP(12)			! IN(SETFILE_SITE):Long term average monthly average 
	REAL BULKDENS(MAXLAYER)		! OUT(CALL):Bulk density of the layer (g/cm3)
      REAL CACCUM					! OUT(CALL): Measured annual C accumulation 
								!     in the soil kgC/ha/yr
      REAL CLAY(MAXLAYER)			! OUT(CALL):Clay content of the layer (%)
	REAL DOMSOILC(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN(SETFILE_SITE):Soil C in 5 
															! major soil series under different LU 
															! in SOM layers (kgC/ha)
	REAL DOMSOILCLAY(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN(SETFILE_SITE):% Soil clay in 5 
															! major soil series under different LU 
	REAL DOMSOILIMPDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN(SETFILE_SITE): Depth of impermeable layer if present (cm)
	INTEGER DOMSOILISIMP(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN(SETFILE_SITE): Flag for impermeable layer (0=No; 1=Yes)
	REAL DOMSOILSILT(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN(SETFILE_SITE):% Soil silt in 5 
															! major soil series under different LU 
	REAL DOMSOILSAND(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN(SETFILE_SITE):% Soil sand in 5 
															! major soil series under different LU 
	REAL DOMSOILBD(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN(SETFILE_SITE):Soil BD in 5 
															! major soil series under different LU 
															! in SOM layers (g/cm3)
	REAL DOMSOILPH(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN(SETFILE_SITE):Soil PH in 5 
												! major soil series under different LU 
												! in SOM layers (g/cm3)
      INTEGER DRAINCLASS			! IN(SETFILE_SITE):Drainage class (not currently used)
	INTEGER EC_EQRUN			! IN(CALL):Initialisation using a full ECOSSE equilibrium run (on or off)
	 INTEGER EC_EQRUN_OFF		!   Initialisation using a full ECOSSE equilibrium run is off
	 INTEGER EC_EQRUN_ON		!   Initialisation using a full ECOSSE equilibrium run is on
	 DATA EC_EQRUN_OFF,EC_EQRUN_ON /0,1/
	
	INTEGER EQMODEL				! IN(CALL):Type of equilibrium run 
	                            !   (NPP, TOC or both)
	 INTEGER EQNPP				!   Model initialised using plant inputs 
								!      from measured NPP
	 INTEGER EQTOC				!   Model initialised using measured TOC
	 INTEGER EQNPPTOC			!   Model initialised using both plant inputs
								!      from measured NPP and measured TOC
	 INTEGER EQHILLIER			!   Model initialised using TOC and
								!      the Hillier solver to get C pools and plant inputs
	 DATA EQNPP,EQTOC,EQNPPTOC,EQHILLIER /1,2,3,4/
      REAL GRIDLAT				! IN(SETFILE_SITE):Average latitude of this 20km2 grid cell - NOT USED
	INTEGER ICMODEL				! IN(READ_MODEL):Type of C initialisation model 
	 INTEGER ICFIXED			!   Initialisation of C is fixed
	 INTEGER ICROTHCEQ			!   Initialisation of C by RothC 
								!    equilibrium run
	INTEGER ILU					! OUT(CALL):Counter for land use 
	REAL IOM(MAXLAYER)			! OUT(CALL):Inert organic C (kgC/ha/layer)
	INTEGER ISDYNPI				! IN(SETFILE_SITE):Code for dynamic PI (adjusted by external factors) - NOT USED
	  INTEGER ISDYNPI_OFF			! Adjustment of PI by external factors is off
	  INTEGER ISDYNPI_ON			! Adjustment of PI by external factors is on
	DATA ISDYNPI_OFF, ISDYNPI_ON /0,1/
 	INTEGER ISERIES				! OUT(CALL):Counter for soil series
      INTEGER ISERROR				! IN(GETSOIL):Code for error in soil 1=error 0=noerror
      INTEGER LUSEQ(MAXGROW)		! IN(SETFILE_SITE):Land use sequence
	INTEGER NSOILNULL			! IN(SETFILE_SITE):Soil code number - NOT USED
	INTEGER NSOMLAY(CELLSTOGRID*MAXSERIES/2,MAXLU)	! OUT(CALL):Number of SOM layers
	INTEGER NXYEARS				! IN(SETFILE_SITE):No.growing seasons simulated
	REAL PI_ANN					! OUT(CALL):Annual plant input (kgC/ha/year)
	                            !        (kgN/ha/year)
	REAL PNREQ(MAXLU)			! IN(SETFILE_SITE):Amount of N required by plant each yr
	                            !        (kgN/ha/year)
	REAL SAND(MAXLAYER)			! OUT(CALL):Sand content of the layer (%)
	REAL SILT(MAXLAYER)			! OUT(CALL):Silt content of the layer (%)
	CHARACTER SITE*80			! IN(OPENCHAN): Name of file containing soil characteristics at site 
	REAL SOILPH(MAXLAYER)		! OUT(CALL):pH of soil in this layer
	REAL SOMDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! OUT(CALL):depth of soil organic matter layers
	REAL TOC(MAXLAYER)			! OUT(CALL):Total organic C (kgC/ha/layer)
	REAL TOTPIC(MAXLU)			! OUT(CALL):Plant C input calculated using 
								!	  ROTHC equilibrium run for each LU (kgC/ha/yr)
C
C Get simulation inputs
C
      CALL SETFILE_SITE(ICMODEL,ISDYNPI,EQMODEL,DRAINCLASS,
     &                 SOMDEPTH,LUSEQ,NXYEARS,
     &                 DOMSOILC,DOMSOILBD,DOMSOILPH,
     &                 DOMSOILCLAY,DOMSOILSILT,DOMSOILSAND,
     &                 DOMSOILISIMP,DOMSOILIMPDEPTH,
     &				 TOTPIC,AVERAIN,AVEPET,AVETEMP,
     &                 CACCUM,GRIDLAT,NSOMLAY,SITE)
      ILU=LUSEQ(1)
C
C ......get soil characteristics 
C
      ISERIES=1
	PNREQ(ILU)=0 ! Set properly
	PI_ANN=TOTPIC(ILU)
	CALL GETSOIL(ISERIES,ILU,SOMDEPTH,NSOMLAY,							
     &             DOMSOILC,DOMSOILCLAY,DOMSOILSILT,				
     &             DOMSOILSAND,DOMSOILBD,DOMSOILPH,
     &             DOMSOILISIMP,DOMSOILIMPDEPTH,
     &             TOC,IOM,CLAY,SILT,SAND,BULKDENS,SOILPH,ISERROR)
	IF(ISERROR .eq. 1)THEN
	  PRINT*,'ERROR! Series ',ISERIES,' not defined for LU ',ILU
	  PRINT*,'Press any key to continue...'
        READ(*,*)
	  STOP
	ENDIF
C
C Leave GET_SITE_SOIL
C
      return
C
C Read soil characteristics in from console
C
	PRINT*,'----------------------------------------------'
	PRINT*,'Enter TOC to 25cm (eg. 195000 kgC/ha)'
      TOC(1)=195000
	READ(*,*)TOC(1)
	TOC(1)=TOC(1)*(MAXDEPTH/MAXLAYER1)/25
	PRINT*,'----------------------------------------------'
	PRINT*,'Enter IOM to 25cm (eg. 27000 kgC/ha)'
	IOM(1)=27061
	READ(*,*)IOM(1)
	IOM(1)=IOM(1)*(MAXDEPTH/MAXLAYER1)/25
	PRINT*,'----------------------------------------------'
303   CONTINUE
	PRINT*,'Enter the steady state land use'
	PRINT*,'       1 = Arable'
	PRINT*,'       2 = Grassland'
	PRINT*,'       3 = Forestry'
	PRINT*,'       4 = Semi-Natural'
	ILU=1
      READ(*,*)ILU
	IF(ILU.LT.1.OR.ILU.GT.4)GOTO 303
	PRINT*,'----------------------------------------------'
	PRINT*,'Enter percent clay in top 25cm (eg 7%)'
	CLAY(1)=70
      READ(*,*)CLAY(1)
	PRINT*,'----------------------------------------------'
	PRINT*,'Enter percent silt in top 25cm (eg 23%)'
	SILT(1)=10
 	READ(*,*)SILT(1)
	PRINT*,'----------------------------------------------'
	PRINT*,'Enter bulk density in top 25cm (eg. 0.34g/cm3)'
	BULKDENS(1)=1.1
	READ(*,*)BULKDENS(1)
	PRINT*,'----------------------------------------------'

      DO 300 IL=1,MAXLAYER1
	  IF((IL*MAXDEPTH/MAXLAYER1).LE.25)THEN
	    TOC(IL)=TOC(1)
	    IOM(IL)=IOM(1)
        ELSE
	    TOC(IL)=0
	    IOM(IL)=0
	  ENDIF
        CLAY(IL)=CLAY(1)
	  SILT(IL)=SILT(1)
	  BULKDENS(IL)=BULKDENS(1)
300   CONTINUE
      IF(EQMODEL.EQ.EQNPP.OR.EQMODEL.EQ.EQNPPTOC.OR.
     &   EC_EQRUN.EQ.EC_EQRUN_ON)THEN
	  PRINT*,'Enter equilibrium plant C input / year (kgC/ha/year)',
     &         ' e.g. 5000 kgC/ha/month'
	  READ(*,*)PI_ANN
	  PRINT*,'----------------------------------------------'
	ELSE
	  PI_ANN=1000
	ENDIF
	IF(EC_EQRUN.EQ.EC_EQRUN_ON)THEN
        PRINT*,'Enter rate of C accumulation in top 25cm (kgC/ha/yr)',
     &         ' +ve: accumulation; -ve: degradation; 0: equilibrium',
     &         ' (eg. 154 kg C / ha / yr)'
	  READ(*,*)CACCUM
	ENDIF
      ISERIES=1
	PNREQ(ILU)=0 ! Set properly
	TOTPIC(ILU)=PI_ANN
	DOMSOILISIMP(ISERIES,ILU)=0
	DOMSOILIMPDEPTH(ISERIES,ILU)=MAXDEPTH
	SOMDEPTH(ISERIES,ILU,1)=25
	NSOMLAY(ISERIES,ILU)=1
C
C Leave GET_SITE_SOIL
C
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GETDEF(ICROP,EXYLD)
C
C Subroutine to return defaults
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      INTEGER MAXCROP
      PARAMETER(MAXCROP=36)
	REAL DEXYLD(MAXCROP)
	INTEGER I
C
C Variables passed to/from calling subroutine
C 
      INTEGER ICROP
      REAL EXYLD
C
C Default values
C Expected Yield: Average yields according to Nix,J., 1995.
C Crop Name                         WW   WB   WOSR WBEA SW   SB   SOSR
      DATA (DEXYLD(I),I=1,MAXCROP)/7.25,6.00,3.00,3.60,5.25,5.00,2.10,
C Crop Name                         SBEA POT  SBEE ROOT LIN  FPEA FLEG
     &                             3.40,40.00,42.5,0.0,1.75,3.70, 0.0,
C Crop Name                         FMAI GMAI SOAT TRIT CERE VPEA BRASS
     &                              0.0, 0.0,4.75,5.25, 0.0,4.50, 0.0,
C Crop Name                         MUST LUPI LENT RICE SCAN SETA FRAP
     &                              0.0, 0.0, 4.5, 0.0, 0.0, 0.0, 0.0,
C Crop Name                         RYEG RYEC UNKN UNKN UNKN UNKN UNKN
     &                             5.25,4.25, 0.0, 0.0, 0.0, 0.0, 0.0,
C Crop Name                         UNKN
     &                              0.0/
      EXYLD=DEXYLD(ICROP)
      END
C
C--------------------------------------------------------
C
      SUBROUTINE GETFERT(IK,NFERT,IFERT,FERT,TFERT,ILAB,IVOL,
     &                  FERTNF,TFERTNF,ILABNF,IVOLNF) 
C
C Subroutine to return any manure application in this timestep
C	
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXFERT				! Max.no.of fertiliser applications allowed 
      PARAMETER (MAXFERT=5)
      INTEGER NF					! Fertiliser counter
	INTEGER IT					! Local counter for fertiliser type
C
C Variables passed to/from this subroutine
C
	INTEGER IK					! IN:No.timesteps from prev.harv. to current
C
C ... Fertiliser factors
C
	INTEGER NFERT				! IN:No.fertiliser applications to this crop
	INTEGER IFERT(MAXFERT)		! IN:No.timesteps to split fert.application
	REAL FERT(MAXFERT)			! IN:Amount of fertiliser applied (kgN/ha)
	REAL TFERT(MAXFERT,3)		! IN:Prop.NO3,NH4,urea in fertiliser
	INTEGER ILAB(MAXFERT)		! IN:Labelling on fertiliser? 0=No 1=Yes
	INTEGER IVOL(MAXFERT)		! IN:N volatilised from fertiliser (kgN/ha)
	REAL FERTNF					! OUT:Amount of fertiliser applied (kgN/ha)
	REAL TFERTNF(3)				! OUT:Prop.NO3,NH4,urea in fertiliser
	INTEGER ILABNF				! OUT:Labelling on fertiliser? 0=No 1=Yes
	INTEGER IVOLNF				! OUT:N volatilised from fertiliser (kgN/ha)
C
      
      DO 100 NF=1,NFERT												
        IF(IK.EQ.IFERT(NF))THEN
          FERTNF=FERT(NF)
		ILABNF=ILAB(NF)
	    IVOLNF=IVOL(NF)
	    DO 200 IT=1,3
            TFERTNF(IT)=TFERT(NF,IT)
200       CONTINUE
          RETURN
        ELSE
	    FERTNF=0
          ILABNF=0
	    DO 300 IT=1,3
            TFERTNF(IT)=0
300       CONTINUE
	  ENDIF											
100   CONTINUE	
      if(NFERT.LE.0)then
       FERTNF=0
          ILABNF=0 
          endif
      END
C
C--------------------------------------------------------
C
      SUBROUTINE GETFYM(IK,NORGM,IORGM,ORGMA,JORGM,IOLAB,
     &                  ORGMANF,JORGMNF,IOLABNF) 
C
C Subroutine to return any manure application in this timestep
C	
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXORGM				! Max.no.of organic manure applications allowed
      PARAMETER (MAXORGM=52)
      INTEGER NF					! Fertiliser counter
C
C Variables passed to/from this subroutine
C
	INTEGER IK					! IN:No.timesteps from prev.harv. to current
	INTEGER NORGM				! IN:No.manure applications to this crop
	INTEGER IORGM(MAXORGM)		! IN:No.timesteps to split manure application
	INTEGER JORGM(MAXORGM)		! IN:Type of manure application
	INTEGER IOLAB(MAXORGM)		! IN:Labelling on manure? 0=No 1=Yes
	REAL ORGMA(MAXORGM)			! IN:Amount of manure applied (t/ha)
	INTEGER JORGMNF				! IN:Type of manure application
	INTEGER IOLABNF				! IN:Labelling on manure? 0=No 1=Yes
	REAL ORGMANF				! IN:Amount of manure applied (t/ha)
C
	norgm=maxorgm
      DO 100 NF=1,NORGM												
        IF(IK.EQ.IORGM(NF))THEN
          ORGMANF=ORGMA(NF)
          JORGMNF=JORGM(NF)
		IOLABNF=IOLAB(NF)      
        ELSE
	    ORGMANF=0
          JORGMNF=0
		IOLABNF=0
	  ENDIF											
100   CONTINUE															
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GETLTA_LIM(LTA_AWC,LTA_TEMP,WMAX,WSAT,FLOWPROP,
     &                  AVERAIN,AVEPET,AVETEMP,ICROP)
C
C Subroutine to return the long term average SMD
C
      IMPLICIT NONE
C 
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	INTEGER MAXDEPTH			! Maximum depth of the soil profile
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
	INTEGER IL					! Local layer counter 
	INTEGER M					! Local month counter
	REAL AVAILWAT(MAXLAYER)		! Available water in the layer (mm/5cm)
	REAL RAINTHISMON			! Rainfall this month (mm)
	REAL REF(MAXLAYER)			! Refill space in the layer (mm/5cm)
	REAL DRAIN(MAXLAYER)		! Drainage from the layer (mm/5cm)
      INTEGER ICOVER,J				! Crop cover 1=Covered 0=Bare
C
C Variables passed to / from this subroutine
C
      REAL AVERAIN(12)			! IN:Average rainfall (mm/month)
	REAL AVEPET(12)				! IN:Average pot.evapotranspirn (mm/month)
	REAL AVETEMP(12)			! IN:Average air temp this month (deg.C)
	REAL LTA_AWC(12,MAXLAYER)	! OUT:Long term ave.soil moist.def.(mm)
	REAL LTA_TEMP(12,MAXLAYER)	! OUT:Long term ave.soil moist.def.(mm)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL FLOWPROP				! IN:Proportion of flow needed to achieve
								!	     water table at depth WTABLE
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
	REAL ICROP(12)
C
C Set long term average soil temp, and initialise available water to FC in Jan
C
      DO 10 IL=1,MAXLAYER1
	  AVAILWAT(IL)=WMAX(IL)
        DO 20 M=1,12
          LTA_TEMP(M,IL)=AVETEMP(M)
20      CONTINUE
10    CONTINUE
C
C Work out available water from rain and pet data
C
      DO j=1,10
	DO 30 M=1,12
      if(icrop(m).gt.0)then 
	  ICOVER=1
         else
       icover=0
         endif    ! Assume soil covered - results in soil being slightly drier 
	  RAINTHISMON=AVERAIN(M)
        CALL RUN_SUNDIAL_WATER_SUBMONTHLY(RAINTHISMON,WMAX,
     &                                        AVAILWAT,WSAT,FLOWPROP,
     &                                        DRAIN,REF,AVEPET(m),
     &                                        AVETEMP(m),ICOVER)
	  DO 40 IL=1,MAXLAYER1
          LTA_AWC(M,IL)=AVAILWAT(IL)
40      CONTINUE
30    CONTINUE
	ENDDO
	END
	 !-----------------    -----------------------
      SUBROUTINE GETLTA(LTA_AWC,LTA_TEMP,WMAX,WSAT,FLOWPROP,
     &                  AVERAIN,AVEPET,AVETEMP)
C
C Subroutine to return the long term average SMD
C
      IMPLICIT NONE
C 
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	INTEGER MAXDEPTH			! Maximum depth of the soil profile
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
	INTEGER IL					! Local layer counter 
	INTEGER M					! Local month counter
	REAL AVAILWAT(MAXLAYER)		! Available water in the layer (mm/5cm)
	REAL RAINTHISMON			! Rainfall this month (mm)
	REAL REF(MAXLAYER)			! Refill space in the layer (mm/5cm)
	REAL DRAIN(MAXLAYER)		! Drainage from the layer (mm/5cm)
      INTEGER ICOVER				! Crop cover 1=Covered 0=Bare
C
C Variables passed to / from this subroutine
C
      REAL AVERAIN(12)			! IN:Average rainfall (mm/month)
	REAL AVEPET(12)				! IN:Average pot.evapotranspirn (mm/month)
	REAL AVETEMP(12)			! IN:Average air temp this month (deg.C)
	REAL LTA_AWC(12,MAXLAYER)	! OUT:Long term ave.soil moist.def.(mm)
	REAL LTA_TEMP(12,MAXLAYER)	! OUT:Long term ave.soil moist.def.(mm)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL FLOWPROP				! IN:Proportion of flow needed to achieve
								!	     water table at depth WTABLE
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
C
C Set long term average soil temp, and initialise available water to FC in Jan
C
      DO 10 IL=1,MAXLAYER1
	  AVAILWAT(IL)=WMAX(IL)
        DO 20 M=1,12
          LTA_TEMP(M,IL)=AVETEMP(M)
20      CONTINUE
10    CONTINUE
C
C Work out available water from rain and pet data
C
      DO 30 M=1,12
	  ICOVER=1	! Assume soil covered - results in soil being slightly drier 
	  RAINTHISMON=AVERAIN(M)
C        CALL DRAIN_SUNDIAL_WATER(RAINTHISMON,WMAX,AVAILWAT,DRAIN,REF,
C     &                               WSAT,FLOWPROP)
        CALL DRAIN_SUNDIAL_WATER2(RAINTHISMON, WMAX, AVAILWAT, WSAT,
     &                            FLOWPROP, DRAIN, REF)
   	  CALL EVAP_SUNDIAL_WATER(AVEPET(M),AVETEMP(M),ICOVER,WMAX,
     &                          AVAILWAT)
         CALL RUN_SUNDIAL_WATER_SUBMONTHLY(RAINTHISMON,WMAX,
     &                                        AVAILWAT,WSAT,FLOWPROP,
     &                                        DRAIN,REF,AVEPET(m),
     &                                        AVETEMP(m),ICOVER)



	  DO 40 IL=1,MAXLAYER1
          LTA_AWC(M,IL)=AVAILWAT(IL)
40      CONTINUE
30    CONTINUE
	END
C
C-------------------------------------------------------------------
C
      SUBROUTINE GETLTA_SWAT(LTA_AWC,LTA_TEMP,WMAX,WSAT,FLOWPROP,
     &                  AVERAIN,AVEPET,AVETEMP)
C
C Subroutine to return the long term average SMD
C
      IMPLICIT NONE
C 
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	INTEGER MAXDEPTH			! Maximum depth of the soil profile
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
	INTEGER IL					! Local layer counter 
	INTEGER M					! Local month counter
	REAL AVAILWAT(MAXLAYER)		! Available water in the layer (mm/5cm)
	REAL REF(MAXLAYER)			! Refill space in the layer (mm/5cm)
	REAL DRAIN(MAXLAYER)		! Drainage from the layer (mm/5cm)
      INTEGER ICOVER				! Crop cover 1=Covered 0=Bare
C
C Variables passed to / from this subroutine
C
      REAL AVERAIN(12)			! IN:Average rainfall (mm/month)
	REAL AVEPET(12)				! IN:Average pot.evapotranspirn (mm/month)
	REAL AVETEMP(12)			! IN:Average air temp this month (deg.C)
	REAL LTA_AWC(12,MAXLAYER)	! OUT:Long term ave.soil moist.def.(mm)
	REAL LTA_TEMP(12,MAXLAYER)	! OUT:Long term ave.soil moist.def.(mm)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL FLOWPROP				! IN:Proportion of flow needed to achieve
								!	     water table at depth WTABLE
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
	REAL Z					! elevation above sea level [m]
	REAL AVP				! IN: water vapor pressure of air at height z (kPa)
	REAL LAT				! IN: Latitude
	INTEGER ROOTLAYER		! Layer until which roots have grown
	REAL CURRENTLAI			! leaf area index for a given day
	REAL CANSTOR			! amount of free water held in the canopy on a given day (mm H2O)
	REAL CURRENTPLBIOM
	REAL CURRENTPLHEI
	REAL SOILFC(MAXLAYER)	! OUT: soil water content at FC (mm/layer)
	REAL SOILW(MAXLAYER)	! IN:Available water (mm/layer)

C
C Set long term average soil temp, and initialise available water to FC in Jan
C
      DO 10 IL=1,MAXLAYER1
	  AVAILWAT(IL)=WMAX(IL)
        DO 20 M=1,12
          LTA_TEMP(M,IL)=AVETEMP(M)
20      CONTINUE
10    CONTINUE
C
C Work out available water from rain and pet data
C
      DO 30 M=1,12
	  ICOVER=1	! Assume soil covered - results in soil being slightly drier 
C        CALL DRAIN_SUNDIAL_WATER(AVERAIN(M),WMAX,AVAILWAT,DRAIN,REF,
C     &                               WSAT,FLOWPROP)
	  CALL DRAIN_SUNDIAL_WATER2(AVERAIN(M),WMAX,AVAILWAT,WSAT,FLOWPROP,
     &                            DRAIN,REF)

  	  CALL EVAP_SUNDIAL_WATER(AVEPET(M),AVETEMP(M),ICOVER,WMAX,
     &                          AVAILWAT)
	  DO 40 IL=1,MAXLAYER1
          LTA_AWC(M,IL)=AVAILWAT(IL)
40      CONTINUE
30    CONTINUE
	END
C
C--------------------------------------------------------
C
      SUBROUTINE GETWEATHER(N_TS,N_STEPS,IFILE,IRYEAR,IK,MEND,SUM_TS,
     &                      RAIN,EVAPW,AIRTEMP,RDD,TMMN,TMMX,VP,WN,
     &                      FIXFILE,MODTYPE,SOILTEMP,
C required for Monte carlo - need to pass on original value ofprecipiation
     &						precip)
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
	REAL precip !original value of rain needs to be passed on from file for MCMODE
      REAL RAIN,EVAPW,AIRTEMP,SOILTEMP(MAXLAYER),RDD,TMMN,TMMX,VP,WN
	CHARACTER*100 FIXFILE(MAXWEATH)
C
C In last week of year open next years weather data file
C
      IF(N_TS.EQ.N_STEPS)THEN
        IFILE=IFILE+1
        CALL NEXTYEAR(IRYEAR,IFILE,IK,MEND,FIXFILE)
      END IF
C
C Read in weather data for this week
C
      SUM_TS = SUM_TS + 1
	IF(MODTYPE.EQ.1)THEN
        READ(3,*)N_TS,RAIN,EVAPW,AIRTEMP, RDD, TMMN, TMMX, VP, WN
	ELSE
	  READ(3,*)N_TS,RAIN,EVAPW,AIRTEMP
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
	precip=RAIN
C
C
C Leave Getweather
C
	END
C
C-----------------------------------------------------------------------
C
	SUBROUTINE GETWEATHER_AT_FC(IFILE,IDATEFC,N_STEPS,IRYEAR,K_TS,
     &                            RAIN,EVAP,AIRTEMP,RDD,TMMN,TMMX,VP,WN,
     &                            FIXFILE,MODTYPE,SOILTEMP)
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
	CHARACTER*100 FIXFILE(MAXWEATH)
C
C Variables required for SWAT water
C
	INTEGER WATERMODEL
	INTEGER SWAT_WATER
	REAL AVP	
C
C Open weather file
C
      IFILE=IFILE+1
      OPEN(3,FILE=FIXFILE(IFILE),STATUS='OLD',ERR=111)
	GOTO 101
111   CONTINUE
      WRITE(*,10)FIXFILE(IFILE)
10	FORMAT('Cannot open weather file ',A100)
	PRINT*,'... press any key to continue'
	READ(*,*)
	STOP
101   CONTINUE
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
          CALL NEXTYEAR(IRYEAR,IFILE,IK,MEND,FIXFILE)
        END IF
C
C read weather data
C
	  IF(MODTYPE.EQ.1)THEN
          READ(3,*)IK,RAIN,EVAP,AIRTEMP, RDD, TMMN, TMMX, VP, WN
	  ELSE
	    READ(3,*)IK,RAIN,EVAP,AIRTEMP
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
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE GETYRSET(IYEAR,NXYEARS,NSOIL,NSOILJ,IAWC,IAWCJ,
     &                IROCK,IROCKJ,INSTRAW,NRESJ,
     &                ICROP,ICROPJ,ISOWN,ISOWNJ,OCROPN,OCROPNJ,
     &                IHARV,IHARVJ,EXYLD,EXYLDJ,NFERT,NFERTJ,
     &                NORGM,NORGMJ,FERT,FERTJ,IFERT,IFERTJ,
     &                ILAB,ILABJ,IVOL,IVOLJ,TFERT,JFERTJ,
     &                ORGMA,ORGMJ,IORGM,IORGMJ,IOLAB,IOLABJ,
     &                JORGM,JORGMJ,LCROP,CFACT,UT,
     &                LHARV,SECONDS,IS_TS,IRYEAR,IFILE,ISTHARV)
C
C Subroutine to set simulation parameters for this year
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      INTEGER MAXCROP,MAXFERT,MAXORGM,MAXWEATH,MAXGROW
      PARAMETER (MAXCROP=36,MAXFERT=5,MAXORGM=52,MAXWEATH=300)
      PARAMETER (MAXGROW=300)
      REAL RPLANT
	INTEGER I,J
C
C Variables passed from calling subroutine
C
	INTEGER NSOIL,NSOILJ,IAWC,IAWCJ,IROCK,IROCKJ,IYEAR
      INTEGER NXYEARS,INSTRAW,ICROP,ISOWN,IHARV,NFERT,NORGM
	INTEGER LHARV,IS_TS,IRYEAR,IFILE,ISTHARV
      INTEGER NRESJ(0:MAXGROW)
      INTEGER ICROPJ(0:MAXGROW),ISOWNJ(0:MAXGROW),IHARVJ(0:MAXGROW)
      INTEGER NFERTJ(0:MAXGROW),NORGMJ(0:MAXGROW),IFERT(MAXFERT)
      INTEGER ILAB(MAXFERT),IFERTJ(0:MAXGROW,MAXFERT)
      INTEGER ILABJ(0:MAXGROW,MAXFERT),IVOL(MAXFERT)
      INTEGER IVOLJ(0:MAXGROW,MAXFERT),IORGM(MAXORGM)
      INTEGER IORGMJ(0:MAXGROW,MAXORGM),IOLAB(MAXORGM)
      INTEGER JORGMJ(0:MAXGROW,MAXORGM),IOLABJ(0:MAXGROW,MAXORGM)
      INTEGER JORGM(MAXORGM),LCROP
	REAL SECONDS
      REAL OCROPN,EXYLD,OCROPNJ(0:MAXGROW),EXYLDJ(0:MAXGROW)
      REAL FERT(MAXFERT),FERTJ(0:MAXGROW,MAXFERT),TFERT(MAXFERT,3)
      REAL CFACT(0:MAXCROP),UT(3,0:MAXCROP),ORGMA(MAXORGM)
      REAL ORGMJ(0:MAXGROW,MAXORGM),JFERTJ(0:MAXGROW,MAXFERT,3)
C
C In first year, set crop counter for mid-July (anthesis of previous crop)
C 
      IF(IYEAR.EQ.0)THEN
        IS_TS=INT((237*24*60*60)/SECONDS)
        IRYEAR=0
        IFILE=0
        INSTRAW=0
        ICROP=0
        ISOWN=0
        OCROPN=0
        IHARV=0
        EXYLD=0
        NFERT=0
        NORGM=0
	ENDIF
C
C Save current harvest date
C
      IF(IYEAR.GT.1)LHARV=IHARV
C
C Reset Soil Parameters
C
      NSOIL=NSOILJ
      IAWC=IAWCJ
      IROCK=IROCKJ
C
C After last crop
C
      IF(IYEAR.EQ.NXYEARS+1)THEN
       INSTRAW=0
       ICROP=0
       ISOWN=0
       OCROPN=0
       IHARV=0
       EXYLD=0
       NFERT=0
       NORGM=0
C
C Current Crop
C Set next crop to current crop
C
      ELSEIF(IYEAR.GT.0)THEN
	INSTRAW=NRESJ(IYEAR)
	ICROP=ICROPJ(IYEAR)
	ISOWN=ISOWNJ(IYEAR)
	OCROPN=OCROPNJ(IYEAR)
	IHARV=IHARVJ(IYEAR)
	EXYLD=EXYLDJ(IYEAR)
	NFERT=NFERTJ(IYEAR)
	NORGM=NORGMJ(IYEAR)
	DO 100 I=1,NFERT
	 FERT(I)=FERTJ(IYEAR,I)
	 IFERT(I)=IFERTJ(IYEAR,I)
         ILAB(I)=ILABJ(IYEAR,I)
         IVOL(I)=IVOLJ(IYEAR,I)
	 DO 200 J=1,3
	   TFERT(I,J)=JFERTJ(IYEAR,I,J)
200      CONTINUE
100     CONTINUE
C
C If crop IS present, and EXYLD is 0, set to default
C
       IF(ICROP.GT.0.AND.EXYLD.LE.0)THEN
	 IF(OCROPN.GT.0)THEN
	   RPLANT=OCROPN/(1+CFACT(ICROP))
	   EXYLD=LOG((RPLANT/UT(1,ICROP))+1)/UT(2,ICROP)
	 ELSEIF(OCROPN.LE.0)THEN
	   CALL GETDEF(ICROP,EXYLD)
	 ENDIF
       ENDIF
C
C Organic Manure
C
      DO 300 I=1,NORGM
       ORGMA(I)=ORGMJ(IYEAR,I)
       IORGM(I)=IORGMJ(IYEAR,I)
       IOLAB(I)=IOLABJ(IYEAR,I)
       JORGM(I)=JORGMJ(IYEAR,I)
300   CONTINUE
C
C .. Previous Crop
C
      ELSE
       ICROP=LCROP
      ENDIF
C
C Set the very first harvest
C
      IF(IYEAR.EQ.1)LHARV=ISTHARV
C
C Leave GETYRSET
C
      RETURN
      END
C
C-----------------------------------------------------------
C
      SUBROUTINE NEXTYEAR(IRYEAR,NFILE,IK,MEND,FIXFILE)
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
C
C Close previously opened weather data file on channel 3
C
      CLOSE(3)
C
C Open new met. data file on channel 3
C
	 IF(IK.LT.MEND)THEN
	 OPEN(3,FILE=FIXFILE(NFILE),STATUS='OLD',ERR=111)
	 GOTO 101
111      CONTINUE
	 WRITE(15,10)FIXFILE(NFILE)
	 WRITE(*,10)FIXFILE(NFILE)
10       FORMAT('Error in weathered data file!'/
     &          'Check format of file ',A12/
     &          'SIMULATION NOT COMPLETED!')
	 PRINT*,' ....Press any key to continue'
	 READ(*,*)
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
C--------------------------------------------------------
C
	SUBROUTINE OPENCHAN(SITE,ENTER_WTD)
C
C Subroutine to open channels
C
      IMPLICIT NONE
C
C Local variable
C
	CHARACTER SETUP*80
	CHARACTER SOIL*80
C
C Variables passed to/from this subroutine
C (IN = originates in other subroutine; 
C  OUT = originates in this subroutine or subroutine called here; 
C  Name in brackets shows first definition of term;
C  CALL indicates passed to/from calling subroutine)
C
	INTEGER ENTER_WTD			! IN(CALL): Code to use read in water table depth 1 = read in, 0=dont read in
	CHARACTER SITE*80			! OUT(CALL): Name of file containing soil characteristics at site 
C
C Open channels 
C
      OPEN(4,FILE='FNAMES.DAT',STATUS='OLD',ERR=111)
C
C Delete message in ERROR.MSG file if one is present from a previous run
C
	OPEN(15,FILE='ERROR.MSG',STATUS='UNKNOWN')
	CLOSE(15,STATUS='DELETE')
C
      OPEN(15,FILE='ERROR.MSG',STATUS='UNKNOWN',ERR=444)
      GOTO 1221
111   CONTINUE
      WRITE(*,110)'FNAMES.DAT'
      GOTO 1221
444   CONTINUE
      WRITE(*,110)'ERROR.MSG'
	READ(*,*)
      STOP
1001  CONTINUE
      WRITE(*,110)'SOILN.OUT'
      GOTO 1221
1002  CONTINUE
      WRITE(*,110)'SOILW.OUT'
1221  CONTINUE
C
C Open Channels
C
c      OPEN(52, FILE='crop_1.out',STATUS='UNKNOWN')
c      OPEN(53, FILE='crop_2.out',STATUS='UNKNOWN')
c      OPEN(54, FILE='crop_3.out',STATUS='UNKNOWN')
c      OPEN(55, FILE='crop_4.out',STATUS='UNKNOWN')
      OPEN(56, FILE='INPUTS.OUT',STATUS='UNKNOWN')
c      OPEN(57, FILE='BALANCE_C.OUT',STATUS='UNKNOWN')
c      OPEN(58, FILE='BALANCE_N.OUT',STATUS='UNKNOWN')
c	OPEN(64, FILE='EVAPO.OUT',  STATUS='UNKNOWN')
c	OPEN(21, FILE='SOILN.OUT',  STATUS='UNKNOWN')
c	OPEN(20, FILE='SOILW.OUT',  STATUS='UNKNOWN')
c	OPEN(23, FILE='DOC.OUT',  STATUS='UNKNOWN')
c	OPEN(24, FILE='SOILWC.OUT',  STATUS='UNKNOWN',ACTION='WRITE')
C
C Get the name of the Simulation Input File from FNAMES.DAT
C
      READ(4,*,ERR=444,END=333)SETUP,SOIL,SITE
C
C Open the Setup File
C
	OPEN(46,FILE=SOIL,STATUS='OLD',ERR=445)
	IF(ENTER_WTD.EQ.1)OPEN(47,FILE='WTD.DAT',STATUS='OLD',ERR=446)
      OPEN(10,FILE=SETUP,STATUS='OLD',ERR=447)
C
C Record SETUP file in results files
C 
	WRITE(*,*)'Setup file = ',SETUP
	WRITE(*,*)'Soil parameter file = ',SOIL
	WRITE(*,*)'File of site characteristics = ',SITE
C	WRITE(57,*)'Setup file = ',SETUP,' Soil parameter file = ',SOIL
	IF(ENTER_WTD.EQ.1)THEN
	  WRITE(*,*)'Water table depth read from file WTD.DAT'
C	  WRITE(57,*)'Water table depth read from file WTD.DAT'
	ENDIF

C	WRITE(57,*)'Soil parameter file = ',SOIL
C	WRITE(57,*)'File of site characteristics = ',SITE
c	WRITE(58,*)'Setup file = ',SETUP,' Soil parameter file = ',SOIL
c	WRITE(58,*)'Soil parameter file = ',SOIL
c	WRITE(58,*)'File of site characteristics = ',SITE
c	WRITE(20,*)'Setup file = ',SETUP,' Soil parameter file = ',SOIL
c	WRITE(20,*)'Soil parameter file = ',SOIL
c	WRITE(20,*)'File of site characteristics = ',SITE
c	WRITE(21,*)'Setup file = ',SETUP,' Soil parameter file = ',SOIL
c	WRITE(21,*)'Soil parameter file = ',SOIL
c	WRITE(21,*)'File of site characteristics = ',SITE
c	WRITE(23,*)'Setup file = ',SETUP,' Soil parameter file = ',SOIL
c	WRITE(23,*)'Soil parameter file = ',SOIL
c	WRITE(23,*)'File of site characteristics = ',SITE
c	WRITE(24,*)'Setup file = ',SETUP,' Soil parameter file = ',SOIL
c	WRITE(24,*)'Soil parameter file = ',SOIL
c	WRITE(24,*)'File of site characteristics = ',SITE
C
C Titles in results files
C
c      WRITE(21,210)
210   FORMAT("Season Step      NO3       NO3       NO3       NO3",
     &       "       NO3       NO3       NO3       NO3",
     &       "       NO3       NO3       NO3       NO3",
     &       "       NO3       NO3       NO3       NO3",
     &       "       NO3       NO3       NO3       NO3",
     &       "       NO3       NO3       NO3       NO3",
     &       "       NO3       NO3       NO3       NO3",
     &       "       NO3       NO3",
     &       "       NH4       NH4       NH4       NH4",
     &       "       NH4       NH4       NH4       NH4",
     &       "       NH4       NH4       NH4       NH4",
     &       "       NH4       NH4       NH4       NH4",
     &       "       NH4       NH4       NH4       NH4",
     &       "       NH4       NH4       NH4       NH4",
     &       "       NH4       NH4       NH4       NH4",
     &       "       NH4       NH4"/
     &       "                0-5cm     5-10cm   10-15cm   15-20cm",
     &       "   20-25cm   25-30cm   30-35cm   35-40cm",
     &       "   40-45cm   45-50cm   50-55cm   55-60cm",
     &       "   60-65cm   65-70cm   70-75cm   75-80cm",
     &       "   80-85cm   85-90cm   90-95cm   95-100cm",
     &       " 100-105cm 105-110cm 110-115cm  115-120cm",
     &       " 120-125cm 125-130cm 130-135cm  135-140cm",
     &       " 140-145cm 145-150cm",
     &       " 0-5cm     5-10cm   10-15cm   15-20cm",
     &       "   20-25cm   25-30cm   30-35cm   35-40cm",
     &       "   40-45cm   45-50cm   50-55cm   55-60cm",
     &       "   60-65cm   65-70cm   70-75cm   75-80cm",
     &       "   80-85cm   85-90cm   90-95cm   95-100cm",
     &       " 100-105cm 105-110cm 110-115cm  115-120cm",
     &       " 120-125cm 125-130cm 130-135cm  135-140cm",
     &       " 140-145cm 145-150cm")
c      WRITE(23,230)
230   FORMAT("Season Step      DOC       DOC       DOC       DOC",
     &       "       DOC       DOC       DOC       DOC",
     &       "       DOC       DOC       DOC       DOC",
     &       "       DOC       DOC       DOC       DOC",
     &       "       DOC       DOC       DOC       DOC",
     &       "       DOC       DOC       DOC       DOC",
     &       "       DOC       DOC       DOC       DOC",
     &       "       DOC       DOC",
     &       "       CO2       CO2       CO2       CO2",
     &       "       CO2       CO2       CO2       CO2",
     &       "       CO2       CO2       CO2       CO2",
     &       "       CO2       CO2       CO2       CO2",
     &       "       CO2       CO2       CO2       CO2",
     &       "       CO2       CO2       CO2       CO2",
     &       "       CO2       CO2       CO2       CO2",
     &       "       CO2       CO2 from DOC"/
     &       "                0-5cm     5-10cm   10-15cm   15-20cm",
     &       "   20-25cm   25-30cm   30-35cm   35-40cm",
     &       "   40-45cm   45-50cm   50-55cm   55-60cm",
     &       "   60-65cm   65-70cm   70-75cm   75-80cm",
     &       "   80-85cm   85-90cm   90-95cm   95-100cm",
     &       " 100-105cm 105-110cm 110-115cm  115-120cm",
     &       " 120-125cm 125-130cm 130-135cm  135-140cm",
     &       " 140-145cm 145-150cm",
     &       " 0-5cm     5-10cm   10-15cm   15-20cm",
     &       "   20-25cm   25-30cm   30-35cm   35-40cm",
     &       "   40-45cm   45-50cm   50-55cm   55-60cm",
     &       "   60-65cm   65-70cm   70-75cm   75-80cm",
     &       "   80-85cm   85-90cm   90-95cm   95-100cm",
     &       " 100-105cm 105-110cm 110-115cm  115-120cm",
     &       " 120-125cm 125-130cm 130-135cm  135-140cm",
     &       " 140-145cm 145-150cm")
c      write(52,*)"  DOY   FNSH   FWSH    PT    AE  PCANA     SLA    NDEM
c     &    NSUP    NUPT"

c	write(53,*)"  DOY    WLV   WLVD    WST    WSO    WRT   WRTD   WST
c     &R  WLVDS    WSH"

c	write(54, *)"  DOY    DS NTOT  NSO  NRT  NSH  NLV  NST 
c     & LAI LAIN   SLN  TSOIL   TSUM"

c	write(55,*)"  DOY   PNC   HNC   RNC   LNC   ONC   NNC  LITC  LITN  
c     &DERI  RTSA  SHSA  ADIF"

C	WRITE(57,571) 
571   FORMAT("    DOY  H-H YEAR      ",
     &       "DPM        DPM        DPM        DPM        DPM        ",
     &       "DPM        DPM        DPM        DPM        DPM        ",
     &       "RPM        RPM        RPM        RPM        RPM        ",
     &       "RPM        RPM        RPM        RPM        RPM        ",
     &       "BIO        BIO        BIO        BIO        BIO        ",
     &       "BIO        BIO        BIO        BIO        BIO        ",
     &       "HUM        HUM        HUM        HUM        HUM        ",
     &       "HUM        HUM        HUM        HUM        HUM        ",
     &       "CO2        CO2        CO2        CO2        CO2        ",
     &       "CO2        CO2        CO2        CO2        CO2        ",
     &       "CH4        CH4        CH4        CH4        CH4        ",
     &       "CH4        CH4        CH4        CH4        CH4        ",
     &       "Tot.DPM    Tot.RPM    Tot.BIO    Tot.HUM    Tot.CO2    ",
     &       "Tot.CH4    CH4.ATM    Tot.IOM  Tot.Input  CO2_for_year",
     &       "Total_C    PLant_Input_C"      )
C	WRITE(57,572),
572	FORMAT("                     0-5cm     5-10cm    10-15cm    ",
     &       "15-20cm    20-25cm    25-30cm    30-35cm    35-40cm    ",
     &       "40-45cm    45-50cm    ",
     &       "0-5cm     5-10cm    10-15cm    15-20cm    20-25cm    ",
     &       "25-30cm    30-35cm    35-40cm    40-45cm    45-50cm    ",
     &     "    0-5cm     5-10cm    10-15cm    15-20cm    20-25cm    ",
     &       "25-30cm    30-35cm    35-40cm    40-45cm    45-50cm    ",
     &     "    0-5cm     5-10cm    10-15cm    15-20cm    20-25cm    ",
     &       "25-30cm    30-35cm    35-40cm    40-45cm    45-50cm    ",
     &     "    0-5cm     5-10cm    10-15cm    15-20cm    20-25cm    ",
     &       "25-30cm    30-35cm    35-40cm    40-45cm    45-50cm    ",
     &     "    0-5cm     5-10cm    10-15cm    15-20cm    20-25cm    ",
     &       "25-30cm    30-35cm    35-40cm    40-45cm    45-50cm    ")
c	WRITE(58,22)"Timestep H-H Year BioHum    PM      DPM     RPM  ",
c     &            "   NO3     NH4     Fert    Stub    ATM     Seed   M",
c     &            "iner   Denit  NitN2O   NitNO PNitN2O Nit15N2O Nit15",
c     &            "NO PNit15N2O DenN2 DenN2O Den15N2 Den15N2O Vol   Cr",
c     &            "op Leach Senes  Tot.In Tot.Out  Diff               "
22    FORMAT(5A51)

c      WRITE(64,*)"  DOY    EVAPO      AE      AT"
C
C If no error, skip error messages
C
      GOTO 1222
C
C Error Messages for file opening
C
110   FORMAT('Warning! Unable to open file ',A12/
     &       'Check computer memory or file format')
222   CONTINUE
      WRITE(15,110)'CYEAR.OUT'
      GOTO 1222
333   CONTINUE
      WRITE(15,110)'N_TS.OUT'
      GOTO 1222
445   CONTINUE
      WRITE(15,110)'SOIL_FILE'
	READ(*,*)
	STOP
      GOTO 1222
446   CONTINUE
      WRITE(15,110)'WTD.DAT'
	WRITE(*,*)'WTD.DAT'
	ENTER_WTD=0
      GOTO 1222
447   CONTINUE
      WRITE(15,110)'SETUP_FILE'
      GOTO 1222
353   CONTINUE
      WRITE(15,110)'NBAL100.OUT'
      GOTO 1222
555   CONTINUE
      WRITE(15,110)'NYEAR.OUT'
      GOTO 1222
666   CONTINUE
      WRITE(15,110)'NFLOW.PLT'
      GOTO 1222
777   CONTINUE
      WRITE(15,110)'N_TS.PLT'
      GOTO 1222
888   CONTINUE
      WRITE(15,110)'NFERT.PLT'
      GOTO 1222
999   CONTINUE
      WRITE(15,110)'NHARV.PLT'
      GOTO 1222
1111  CONTINUE
      WRITE(15,110)'N15_TS.PLT'
1222  CONTINUE
      END

C
C-------------------------------------------------------------------
C
      SUBROUTINE PUTRES(IYEAR,IK,SOILN,AMMN,SOILW,MOBDOC,CO2FROMDOC,
     &					WILTPOINT,WATCONT,WATERMODEL,
     &					CURRENTPLBIOM,CURRENTPLHEI,CURRENTLAI)
	IMPLICIT NONE
	INTEGER MAXLAYER,MAXLAYER1
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/
	INTEGER IYEAR,IK,I,WATERMODEL
	INTEGER SUNDIAL_WATER
	INTEGER SWAT_WATER
	DATA SUNDIAL_WATER,SWAT_WATER /0,1/
	REAL SOILN(MAXLAYER),SOILW(MAXLAYER),AMMN(MAXLAYER)
	REAL MOBDOC(MAXLAYER),CO2FROMDOC(MAXLAYER),WILTPOINT(MAXLAYER)
	REAL WATCONT(MAXLAYER),WATCONT_OUT(MAXLAYER)
	REAL CURRENTPLBIOM			! Current plant biomass [ kg ha-1]
	REAL CURRENTPLHEI			! Current plant height [cm]
	REAL CURRENTLAI				! Current LAI
	WRITE(21,10)IYEAR,IK,(SOILN(I),I=1,MAXLAYER1),
     &                     (AMMN(I),I=1,MAXLAYER1)
c	WRITE(20,10)IYEAR,IK,(SOILW(I),I=1,MAXLAYER1)
c	WRITE(23,10)IYEAR,IK,(MOBDOC(I),I=1,MAXLAYER1),
c     &                     (CO2FROMDOC(I),I=1,MAXLAYER1)
	IF(WATERMODEL.EQ.SWAT_WATER)THEN
	DO I=1,MAXLAYER1,1
		IF(WATCONT(I).GE.WILTPOINT(I))THEN
			WATCONT_OUT(I)=WATCONT(I)-WILTPOINT(I)
		ELSE
C			WATCONT_OUT(I)=SOILW(I)-(WILTPOINT(I)-WATCONT(I))
			WATCONT_OUT(I)=WILTPOINT(I)-WATCONT(I)
		END IF
	END DO
	WRITE(30,10)IYEAR,IK,(WATCONT_OUT(I),I=1,MAXLAYER1)
	WRITE(31,10)IYEAR,IK,CURRENTPLBIOM,CURRENTPLHEI,CURRENTLAI

	ELSE
	END IF
10    FORMAT(I5,I5,120F10.2)
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE READ_MODEL(DMODEL,ICMODEL,INMODEL,CNMODEL,DOCMODEL,
     &                SPARMODEL,EQMODEL,IMFUNC,ITFUNC,CH4MODEL,EC_EQRUN)
C
C Subroutine to get filenames and simulation inputs
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
C
C Variables passed to/from calling subroutine
C
	INTEGER CH4MODEL			! OUT:Methane model (off, Richards or Aitkenhead model) 
	INTEGER CNMODEL				! OUT:Type of Model for C:N ratio of DPM and 
	                            !     RPM (MAGEC or Foereid)
	INTEGER DMODEL				! OUT:Type denitrification model chosen 
	                            ! (Bradbury or NEMIS)
	INTEGER DOCMODEL			! OUT:DOC model (on or off)
	INTEGER EC_EQRUN			! OUT:Initialisation using a full ECOSSE equilibrium run (on or off)
	INTEGER EQMODEL				! OUT:Type of equilibrium run 
	                            ! (NPP, TOC or both)
	INTEGER ICMODEL				! OUT:Type of C initialisation model 
	INTEGER INMODEL				! OUT:Type of N initialisation model
	INTEGER IMFUNC				! OUT:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
	INTEGER ITFUNC				! OUT:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER SPARMODEL			! OUT:Soil parameter model (from file or calc)
C
C Get mode of simulation if specified
C 1. Denitrification model, DMODEL: DNEMIS = NEMIS model; DBRADBURY = Bradbury model
C
      READ(4,*,ERR=111,END=111)DMODEL
C
C 2a. Initialisation of SOM pools, ICMODEL: ICFIXED = fixed initialisation; 
C                                           ICROTHCEQ = RothC equilibrium run.
C
	READ(4,*,ERR=111,END=111)ICMODEL
C
C 2b. Initialisation of SOM pools, INMODEL: INSTABLECN = Bradbury's assumption of stable C:N ratio
C                                           INPASSCN = passed C:N ratio of DPM and RPM
C
	READ(4,*,ERR=111,END=111)INMODEL
C
C 3. DOC model, DOCMODEL: DOC_OFF	= DOC model off; DOC_ON = DOC model on
C
	READ(4,*,ERR=111,END=111)DOCMODEL
C
C 4. CN model, CNMODEL: CNFOEREID = C:N ratio obtained by method of Foereid; C:N ratio obtained by method of MAGEC
C
	READ(4,*,ERR=111,END=111)CNMODEL
C
C 5. Soil parameter model, SPARMODEL: SPARFILE = Soil parameters read in from file; SPARCALC = Soil parameters calculated from TOC etc
C
	READ(4,*,ERR=111,END=111)SPARMODEL
C
C 6. Type of equilibrium run: EQNPP = Model initialised using plant inputs from measured NPP;
C							EQTOC = Model initialised using measured TOC;
C							EQNPPTOC = Model initialised using both plant inputs.
C
	READ(4,*,ERR=111,END=111)EQMODEL
C
C 7. Calculation of moisture rate modifiers, IMFUNC_ROTHC=0, IMFUNC_HADLEY=1
C    Calculation of temp.rate modifiers, ITFUNC_ROTHC=0, ITFUNC_HADLEY=1
C
      READ(4,*,ERR=111,END=111)IMFUNC
      READ(4,*,ERR=111,END=111)ITFUNC
C
C 8. CH4 model, CH4MODEL: CH4_OFF	= CH4 model off; CH4_RICHARDS = Richards CH4 model on; CH4_AITKENHEAD = Aitkenhead CH4 model on
C
      READ(4,*,ERR=111,END=111)CH4MODEL
C
C 9. ECOSSE equilibrium run,EC_EQRUN
C
      READ(4,*,ERR=111,END=111)EC_EQRUN
C
C If model type not specified in file, jump to end
C
111   CONTINUE
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SAVE_SITE_SOIL(ISAVE,LU1,
     &                        SOILN,SOIL15,AMMN,AMMN15,
     &                        SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
     &                        DPMCARB0,DPMNIT0,DPMNLAB0,
     &                        RPMCARB0,RPMNIT0,RPMNLAB0,
     &                        BCARB0,BNIT0,BNLAB0,
     &                        HCARB0,HNIT0,HNLAB0,
     &                        IOM,PI_CEQ_MON)
C
C Subroutine to save/retrieve soil characteristics
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)

	REAL AMMN1(MAXLAYER,MAXLU)		! Soil ammonium-N (kgN/ha/layer)
	REAL AMMN151(MAXLAYER,MAXLU)	! Soil ammonium-N15 (kgN15/ha/layer)
	REAL BCARB01(MAXLAYER,MAXLU)	! C in soil biomass (kgC/ha/layer)
      REAL BNIT01(MAXLAYER,MAXLU)		! N in soil biomass (kgN/ha/layer)
	REAL BNLAB01(MAXLAYER,MAXLU)	! N15 in soil humus (kgN15/ha/layer)
	REAL DPMCARB01(MAXLAYER,MAXLU)	! C in decomposable PM (kgC/ha/layer)
	REAL DPMNIT01(MAXLAYER,MAXLU)	! N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB01(MAXLAYER,MAXLU)	! N15 in decomposable PM (kgN15/ha/layer)
	REAL FLOWPROP1(MAXLU)			! Proportion of flow needed to achieve
									!	water table at depth WTABLE
	REAL HCARB01(MAXLAYER,MAXLU)	! C in soil humus (kgC/ha/layer)
	REAL HNIT01(MAXLAYER,MAXLU)		! N in soil humus (kgN/ha/layer)
      REAL HNLAB01(MAXLAYER,MAXLU)	! N15 in soil biomass (kgN15/ha/layer)
	INTEGER IL						! Local layer counter variable
	INTEGER IMON					! Local month counter
	REAL IOM1(MAXLAYER,MAXLU)		! Inert organic C (kgC/ha/layer)
	REAL LTA_AWC1(12,MAXLAYER,MAXLU)! IN:Long term available water (mm)
	REAL LTA_TEMP1(12,MAXLAYER,MAXLU) ! IN:Average air temp this month (deg.C)
	REAL PI_CEQ_MON1(12,MAXLAYER,MAXLU)	! Equilibrium plant C input each month
									! in each layer (kgC/ha/month/layer)
	REAL RPMCARB01(MAXLAYER,MAXLU)	! C in resistant PM (kgC/ha/layer)
	REAL RPMNIT01(MAXLAYER,MAXLU)	! N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB01(MAXLAYER,MAXLU)	! N15 in resistant PM (kgN15/ha/layer)
	REAL SOIL151(MAXLAYER,MAXLU)	! Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILN1(MAXLAYER,MAXLU)		! Soil nitrate-N (kgN/ha/layer)
	REAL SOILW1(MAXLAYER,MAXLU)		! Available water (mm/layer)
	REAL WSAT1(MAXLAYER,MAXLU)			! Available water at saturation (mm/layer)
	REAL WMAX1(MAXLAYER,MAXLU)			! Available water at field capacity (mm/layer)
C
C Variables passed to / from this routine
C 
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
	REAL BCARB0(MAXLAYER)		! IN:C in soil biomass (kgC/ha/layer)
      REAL BNIT0(MAXLAYER)		! IN:N in soil biomass (kgN/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN:N15 in soil humus (kgN15/ha/layer)
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL FLOWPROP				! IN:Proportion of flow needed to achieve
							    !	     water table at depth WTABLE
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN:N in soil humus (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN:N15 in soil biomass (kgN15/ha/layer)
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
      INTEGER ISAVE				! IN:Code to save or retrieve variables
	REAL LTA_AWC(12,MAXLAYER)	! IN:Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)	! IN:Average air temp this month (deg.C)
      INTEGER LU1					! IN: Land use code
	REAL PI_CEQ_MON(12,MAXLAYER)! IN/OUT: Equilibrium plant C input each month
								!         in each layer (kgC/ha/month/layer)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
	REAL WMAX(MAXLAYER)			! IN:Available water at field capacity(mm/layer)
C
C Save soil 
C
      IF(ISAVE.EQ.1)THEN
	  FLOWPROP1(LU1)=FLOWPROP
	  DO 100 IL=1,MAXLAYER1
	    AMMN1(IL,LU1)=AMMN(IL)
		AMMN151(IL,LU1)=AMMN15(IL)
		BCARB01(IL,LU1)=BCARB0(IL)
		BNIT01(IL,LU1)=BNIT0(IL)
		BNLAB01(IL,LU1)=BNLAB0(IL)
		DPMCARB01(IL,LU1)=DPMCARB0(IL)
		DPMNIT01(IL,LU1)=DPMNIT0(IL)
		DPMNLAB01(IL,LU1)=DPMNLAB0(IL)
		HCARB01(IL,LU1)=HCARB0(IL)
		HNIT01(IL,LU1)=HNIT0(IL)
		HNLAB01(IL,LU1)=HNLAB0(IL)
		RPMCARB01(IL,LU1)=RPMCARB0(IL)
		RPMNIT01(IL,LU1)=RPMNIT0(IL)
		RPMNLAB01(IL,LU1)=RPMNLAB0(IL)
	    IOM1(IL,LU1)=IOM(IL)
		SOIL151(IL,LU1)=SOIL15(IL)
		SOILN1(IL,LU1)=SOILN(IL)
		SOILW1(IL,LU1)=SOILW(IL)
	    WSAT1(IL,LU1)=WSAT(IL)
	    WMAX1(IL,LU1)=WMAX(IL)
	    DO 200 IMON=1,12
	 	  LTA_AWC1(IMON,IL,LU1)=LTA_AWC(IMON,IL)
	      LTA_TEMP1(IMON,IL,LU1)=LTA_TEMP(IMON,IL)
		  PI_CEQ_MON1(IMON,IL,LU1)=PI_CEQ_MON(IMON,IL)
200       CONTINUE
100     CONTINUE
C
C Retrieve soil
C
	ELSEIF(ISAVE.EQ.0)THEN
	  FLOWPROP=FLOWPROP1(LU1)
	  DO 300 IL=1,MAXLAYER1
	    AMMN(IL)=AMMN1(IL,LU1)
		AMMN15(IL)=AMMN151(IL,LU1)
		BCARB0(IL)=BCARB01(IL,LU1)
		BNIT0(IL)=BNIT01(IL,LU1)
		BNLAB0(IL)=BNLAB01(IL,LU1)
		DPMCARB0(IL)=DPMCARB01(IL,LU1)
		DPMNIT0(IL)=DPMNIT01(IL,LU1)
		DPMNLAB0(IL)=DPMNLAB01(IL,LU1)
		HCARB0(IL)=HCARB01(IL,LU1)
		HNIT0(IL)=HNIT01(IL,LU1)
		HNLAB0(IL)=HNLAB01(IL,LU1)
		RPMCARB0(IL)=RPMCARB01(IL,LU1)
		RPMNIT0(IL)=RPMNIT01(IL,LU1)
		RPMNLAB0(IL)=RPMNLAB01(IL,LU1)
	    IOM(IL)=IOM1(IL,LU1)
		SOIL15(IL)=SOIL151(IL,LU1)
		SOILN(IL)=SOILN1(IL,LU1)
		SOILW(IL)=SOILW1(IL,LU1)
	    WSAT(IL)=WSAT1(IL,LU1)
	    WMAX(IL)=WMAX1(IL,LU1)
	    DO 400 IMON=1,12
	 	  LTA_AWC(IMON,IL)=LTA_AWC1(IMON,IL,LU1)
	      LTA_TEMP(IMON,IL)=LTA_TEMP1(IMON,IL,LU1)
		  PI_CEQ_MON(IMON,IL)=PI_CEQ_MON1(IMON,IL,LU1)
400       CONTINUE
300     CONTINUE
	ENDIF
C
C Leave SAVE_SITE_SOIL
C
      END

C
C-------------------------------------------------------------
C
      SUBROUTINE SETFILE(ATM,IDATEFC,NSOILJ,IDRAINJ,IROCKJ,LCROP,
     &                   PREYLD,SECONDS,MODTYPE,N_STEPS,
     &                   IAWCJ,NYEARS,
     &                   ISTHARV,ISTYR,ISTART_TS,ILAST_TS,FIXEND,LAT,
     &                   FIXFILE,NGROW,NXYEARS,ICROPJ,ISOWNJ,OCROPNJ,
     &                   IHARVJ,EXYLDJ,NRESJ,NFERTJ,NORGMJ,FERTJ,IFERTJ,
     &                   JFERTJ,IVOLJ,ILABJ,ORGMJ,IORGMJ,JORGMJ,
     &                   IOLABJ,WTABLE,
     &                   NCULT,ICULT,JCULT,
C Variable required for SWAT water
     &					FIXFILEAVP,WATERMODEL,SWAT_WATER,Z)
C
C Subroutine to get filenames and simulation inputs
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      INTEGER MAXCROP,MAXFERT,MAXORGM,MAXWEATH,MAXGROW
      PARAMETER (MAXCROP=36,MAXFERT=5,MAXORGM=52,MAXWEATH=300)
      PARAMETER (MAXGROW=300)
	INTEGER I,J,IGRASSJ
	INTEGER TIMESTEP
C
C Variables passed from calling subroutine
C
	CHARACTER*100 FIXFILE(MAXWEATH)
	INTEGER IDATEFC,NSOILJ,IDRAINJ,IROCKJ,LCROP
      INTEGER MODTYPE,N_STEPS,IAWCJ,NYEARS,ISTHARV,ISTYR
      INTEGER ISTART_TS,ILAST_TS,FIXEND,NGROW,NXYEARS
      INTEGER ICROPJ(0:MAXGROW),ISOWNJ(0:MAXGROW),IHARVJ(0:MAXGROW)
      INTEGER NRESJ(0:MAXGROW),NFERTJ(0:MAXGROW),NORGMJ(0:MAXGROW)
      INTEGER IFERTJ(0:MAXGROW,MAXFERT)
      INTEGER IVOLJ(0:MAXGROW,MAXFERT),ILABJ(0:MAXGROW,MAXFERT)
      INTEGER IORGMJ(0:MAXGROW,MAXORGM),JORGMJ(0:MAXGROW,MAXORGM)
      INTEGER IOLABJ(0:MAXGROW,MAXORGM)
	REAL SECONDS
	REAL ATM,PREYLD,OCROPNJ(0:MAXGROW)
      REAL EXYLDJ(0:MAXGROW),FERTJ(0:MAXGROW,MAXFERT)
      REAL ORGMJ(0:MAXGROW,MAXORGM),LAT,JFERTJ(0:MAXGROW,MAXFERT,3)
	REAL WTABLE					! IN:Water table depth in cm
	CHARACTER*100 FIXFILEAVP(MAXWEATH) !required for SWAT water
	INTEGER WATERMODEL
	INTEGER SWAT_WATER
	REAL Z !altitude of the site, required for SWAT
	INTEGER NCULT		   ! Number of cultivations
	INTEGER ICULT(MAXGROW) ! Timesteps from 01/01 when cultivation occured
	INTEGER JCULT(MAXGROW) ! Type of cultivation 
						   ! 0=Zero tillage (with or without mulching) may also be included
						   ! 1=Minimum tillage (or eco-tillage/conservation): Non inversion but depth of cultivation up to 5-10cm, normally this tillage makes a furrow for sowing of seeds and fertilizer application. This practice is applicable for both spring and cover/winter crops.
						   ! 2=Reduced tillage (or non-inversion): Cover both type i.e. non-inversion but depth of cultivation up to 15-20 cm and less number of tillage practices compared to conventional one, where applicable.
						   ! 3=Conventional  tillage (inversion): Depth of cultivation up to 20-30 cm, no limitation of tillage intensity and plough type.
C
C Read in Simulation Data...
C ...Default values
C
      ATM=40
      IDATEFC=1
C
C Values read from SETUP file
C
	READ(10,*,ERR=111)NSOILJ,IDRAINJ,IROCKJ,LCROP,PREYLD,IGRASSJ,
     &                  ATM,IDATEFC,TIMESTEP,MODTYPE
	if(modtype.eq.0)then
	  write(*,*)' You are using Sundials old crop model'
	else
	  write(*,*)' You are using Sundial with Yins crop model'
	endif
      IF(TIMESTEP.EQ.0)THEN
	  SECONDS=30*60
	ELSEIF(TIMESTEP.EQ.1)THEN
	  SECONDS=24*60*60
	ELSEIF(TIMESTEP.EQ.2)THEN
	  SECONDS=7*24*60*60
	ELSEIF(TIMESTEP.EQ.3)THEN
	  SECONDS=(365.25/12)*24*60*60
	ENDIF
      N_STEPS=365.25*24*60*60/SECONDS
      IAWCJ=NSOILJ
      ATM=ATM/N_STEPS
C
C ...Weather Data
C
      READ(10,*,ERR=111)NYEARS,ISTHARV,ISTYR,ISTART_TS,ILAST_TS,
     &                  FIXEND,LAT,WTABLE
 	IF(WATERMODEL.EQ.SWAT_WATER)THEN
	  READ(10,*,ERR=111) Z
	ENDIF

      IF(FIXEND.EQ.0)FIXEND=ILAST_TS
      DO 100 I=1,NYEARS,1
        READ(10,*,ERR=111)FIXFILE(I)
	IF(WATERMODEL.EQ.SWAT_WATER)THEN
	  READ(10,*,ERR=111)FIXFILEAVP(I)
	ENDIF
100   CONTINUE
C
C ...Crop Data
C
      READ(10,*,ERR=111)NGROW
      NXYEARS=NGROW
      DO 200 I=1,NGROW,1
      READ(10,*,ERR=111)ICROPJ(I),ISOWNJ(I),OCROPNJ(I),IHARVJ(I),
     &                 EXYLDJ(I),NRESJ(I),NFERTJ(I),NORGMJ(I)
      DO 300 J=1,NFERTJ(I)
      READ(10,*,ERR=111)FERTJ(I,J),IFERTJ(I,J),
     &                 JFERTJ(I,J,1),JFERTJ(I,J,2),JFERTJ(I,J,3),
     &                 IVOLJ(I,J),ILABJ(I,J)
300   CONTINUE
C
C ...Organic Manure Data
C
      DO 400 J=1,NORGMJ(I)
       READ(10,*,ERR=111)ORGMJ(I,J),IORGMJ(I,J),JORGMJ(I,J),IOLABJ(I,J)
400   CONTINUE
200   CONTINUE
C
C Close the Parameter Simulation File
C
      CLOSE(10)
C
C Bypass reading in from MANAGEMENT.DAT if no error found
C
      GOTO 222
C
C Go straight here if error in SETUP file and read from MANAGEMENT.DAT 
C style file instead
C
111   CONTINUE
      REWIND(10)
C
C Value read in from MANAGEMENT.DAT
C
C (1) Soil type, drainage class, depth to impermeable layer, previous crop, previous yield, 
C     atmospheric N input (kg N / ha / yr), date of field capacity (1=01/01; 2=01/06), timestep (0 = 30 minute; 1 = daily; 2 = weekly; 3 = monthly), Crop model type (0=SUNDIAL; 1=MAGEC)			
C     number of years included in the simulation, Timesteps from 01/01/01 to harvest of previous crop, First year of simulation, End of simulation (in timesteps), 
C     fixed end of simulation (0=No; 1=Yes), Latitude, Water Table of soil, if > 300 cm than there is no effect
C
      READ(10,30,ERR=444)NSOILJ,IDRAINJ,IROCKJ,LCROP,PREYLD,
     &                   ATM,IDATEFC,TIMESTEP,MODTYPE,
     &                   NYEARS,ISTHARV,ISTYR,ILAST_TS,
     &                   FIXEND,LAT,WTABLE
30    FORMAT(4(I10/),2(F10.0/),8(I10/),(F10.0/),F10.0)
      IF(TIMESTEP.EQ.0)THEN
	  SECONDS=30*60
	ELSEIF(TIMESTEP.EQ.1)THEN
	  SECONDS=24*60*60
	ELSEIF(TIMESTEP.EQ.2)THEN
	  SECONDS=7*24*60*60
	ELSEIF(TIMESTEP.EQ.3)THEN
	  SECONDS=(365.25/12)*24*60*60
	ENDIF
      N_STEPS=365.25*24*60*60/SECONDS
C
C (2) blank line  for SWAT
C
	IF(WATERMODEL.EQ.SWAT_WATER)THEN
	  READ(10,*,ERR=444)
	ENDIF
C
C (3) Read in weather files
C
      DO 450 I=1,NYEARS,1
        READ(10,*,ERR=444)FIXFILE(I)
		IF(WATERMODEL.EQ.SWAT_WATER)THEN
	      READ(10,*,ERR=111)FIXFILEAVP(I)
	    ENDIF
450   CONTINUE
C
C (4) Read in cropping details
C (4a) Number of seasons
C
      READ(10,35,ERR=444)NGROW
35    FORMAT(I10)
      NXYEARS=NGROW
C
C (4b) For each season...
C
      DO 500 I=1,NGROW,1
C
C      ...Crop type, sowing date, N offtake (kg N / ha), harvest date
C      expected yield (t/ha), straw incorporation (0=no,1=yes), number of fertiliser application, number of manure applications
C
      READ(10,40,ERR=444)ICROPJ(I),ISOWNJ(I),OCROPNJ(I),IHARVJ(I),
     &                 EXYLDJ(I),NRESJ(I),NFERTJ(I),NORGMJ(I)
40    FORMAT(2(I10/),F10.0/I10/F10.0/2(I10/),I10)
C
C (4c) For each fertiliser application...
C
      DO 600 J=1,NFERTJ(I)
C
C      ... amount of fertiliser(kg N / ha), date applied,
C      prop. nitrate, prop. ammonium, prop. urea,
C      does fert.contain ammonium salts other than ammonium sulphate? (0=No; 1=Yes), Is fertiliser labelled (0=No; 1=Yes)
C
        READ(10,50,ERR=444)FERTJ(I,J),IFERTJ(I,J),
     &                 JFERTJ(I,J,1),JFERTJ(I,J,2),JFERTJ(I,J,3),
     &                 IVOLJ(I,J),ILABJ(I,J)
50      FORMAT(F10.0/I10/3(F10.0/),I10/I10)
600   CONTINUE
C
C (4d) For each manure application...
C
      DO 700 J=1,NORGMJ(I)
C
C      ... amount of manure (t manure / ha), date applied, type manure, Is manure labelled (0=No; 1=Yes)
C
       READ(10,60,ERR=444)ORGMJ(I,J),IORGMJ(I,J),JORGMJ(I,J),IOLABJ(I,J)
60     FORMAT(F10.0/2(I10/),I10)
700   CONTINUE      
500   CONTINUE
C
C (4e) Number of cultivations
C
      READ(10,35,END=555,ERR=666)NCULT
      GOTO 777
666   PRINT*,'Error in entered number of cultivations'
555   CONTINUE
      NCULT=0
777   CONTINUE
C
C (4f) For each cultivation...
C
      DO 800 I=1,NCULT
	  READ(10,70,ERR=888)ICULT(I),JCULT(I)
800   CONTINUE     
      GOTO 999
888   CONTINUE
      PRINT*,'Error in format of cultivation ',I
	READ(*,*)
	STOP
999   CONTINUE
70    FORMAT(I10/I10)
C
C Finish up
C
      IAWCJ=NSOILJ
      ATM=ATM/N_STEPS
	CLOSE(10)
	RETURN
444   CONTINUE
      WRITE(15,10)
	WRITE(*,10)
10    FORMAT('Error in simulation setup!'/
     &       'Check the format of file'/
     &       'SIMULATION NOT COMPLETED!'//
     &       '....press any to continue')
	READ(*,*)
      STOP
222   CONTINUE
333   CONTINUE
C
C Leave SETFILE
C
      RETURN
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE SETFILE_SITE(
     &                 ICMODEL,ISDYNPI,EQMODEL,DRAINCLASS,
     &                 SOMDEPTH,LUSEQ,NXYEARS,
     &                 DOMSOILC,DOMSOILBD,DOMSOILPH,
     &                 DOMSOILCLAY,DOMSOILSILT,DOMSOILSAND,
     &                 DOMSOILISIMP,DOMSOILIMPDEPTH,
     &				 TOTPIC,AVERAIN,AVEPET,AVETEMP,
     &                 CACCUM,GRIDLAT,NSOMLAY,SITE)
C
C Subroutine to get filenames and simulation inputs
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/
	INTEGER MAXGROW				! Max.no.of growing seasons allowed
	PARAMETER (MAXGROW=300)
	INTEGER MAXLU				! Max.no.of land use types
      PARAMETER (MAXLU=6)
	INTEGER MAXLU1				! Max.no.of land use types
	DATA MAXLU1 /6/
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXSOMLAY			! Maximum number of SOM layers
	PARAMETER (MAXSOMLAY=10)			

	INTEGER ILU					! Counter for land use 
	INTEGER IMON				! Counter for months
      real INMODE				! Input mode of model run 
								!		1=use plant inputs; 
								!		2=use toc; 
								!		3=use plant inputs and toc; 
								!		4=use C accumulation;
								!		5=use fixed C pools;
								!		6=use toc and the Hillier solver 
								!		  (an analytical solution of RothC =brium run)
	CHARACTER*50 INFILE			! Input file name
	INTEGER ISOMLAY				! Current SOM layer
	INTEGER IYEAR				! Current growing season number
	REAL TEMP					! Temporary real number

C
C Variables passed to/from calling subroutine
C ...Model factors
C
	INTEGER EQMODEL				! OUT:Type of equilibrium run 
	                            ! (NPP, TOC or both)
	INTEGER EQNPP				! Model initialised using plant inputs 
								!    from measured NPP
	INTEGER EQTOC				! Model initialised using measured TOC
	INTEGER EQNPPTOC			! Model initialised using both plant inputs
								!    from measured NPP and measured TOC
	INTEGER EQHILLIER			! Model initialised using TOC and
								!    the Hillier solver to get C pools and plant inputs
	DATA EQNPP,EQTOC,EQNPPTOC,EQHILLIER /1,2,3,4/
	INTEGER ISCHECKED(MAXLU)
	INTEGER ICMODEL				! OUT:Type of C initialisation model 
	INTEGER ICFIXED				! Initialisation of C is fixed
	INTEGER ICROTHCEQ			! Initialisation of C by RothC 
								!  equilibrium run
	DATA ICFIXED,ICROTHCEQ /1,2/
	INTEGER ISDYNPI				! OUT(CALL):Code for dynamic PI (adjusted by external factors)
	 INTEGER ISDYNPI_OFF		!   Adjustment of PI by external factors is off
	 INTEGER ISDYNPI_ON			!   Adjustment of PI by external factors is on
	 DATA ISDYNPI_OFF, ISDYNPI_ON /0,1/
	CHARACTER SITE*80			! IN(CALL): Name of file containing soil characteristics at site 
C
C ...Soil factors
C
      REAL CACCUM					! OUT: Measured annual C accumulation 
								!     in the soil kgC/ha/yr
	REAL DOMSOILC(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! OUT:Soil C in 5 
												! major soil series under different LU 
												! in SOM layers (kgC/ha)
	REAL DOMSOILCLAY(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! OUT:% Soil clay in 5 
												! major soil series under different LU 
	REAL DOMSOILIMPDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN: Depth of impermeable layer if present (cm)
	INTEGER DOMSOILISIMP(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN: Flag for impermeable layer (0=No; 1=Yes)
	REAL DOMSOILSILT(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! OUT:% Soil silt in 5 
												! major soil series under different LU 
	REAL DOMSOILSAND(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! OUT:% Soil sand in 5 
												! major soil series under different LU 
	REAL DOMSOILBD(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! OUT:Soil BD in 5 
												! major soil series under different LU 
												! in SOM layers (g/cm3)
	REAL DOMSOILPH(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! OUT:Soil PH in 5 
												! major soil series under different LU 
												! in SOM layers (g/cm3)
      INTEGER DRAINCLASS			! OUT:Drainage class
      INTEGER ISERROR				! IN:Code for error in soil 1=error 0=noerror
      INTEGER ISMISS				! IN:Code for missing values in soil 1=missing 0=nomissing
      INTEGER LUSEQ(MAXGROW)		! OUT:Land use sequence
	INTEGER NSOMLAY(CELLSTOGRID*MAXSERIES/2,MAXLU)	! OUT:Number of SOM layers
	INTEGER NXYEARS				! OUT:No.growing seasons simulated
      REAL SOMDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! OUT:Depth of each SOM layer
C
C ...Plant factors
C
	REAL TOTPIC(MAXLU)			! OUT:Plant C input calculated using 
								!	  ROTHC equilibrium run for each LU (kgC/ha/yr)
C
C ...Met factors
C
      REAL AVERAIN(12)			! OUT:Long term average rainfall(mm/month)
      REAL AVEPET(12)				! OUT:Long term average PET (mm/month)
      REAL AVETEMP(12)			! OUT:Long term average monthly average 
								!	air temp (deg.C)
      REAL GRIDLAT				! IN/OUT:Average latitude of this 20km2 grid cell
C
C Soil type
C
      ISDYNPI=ISDYNPI_OFF
C
C Open data file
C
111   CONTINUE
      OPEN(42,FILE=SITE,STATUS='OLD',ERR=112)
      GOTO 113
C
C Get input file name
C
112   CONTINUE
      PRINT*,'Enter the name of the file containing input data'
      PRINT*,'======================================================'
	PRINT*,'  Mode of equilibrium run: calculate C pools...'
	PRINT*,'          1=using known PI' 
	PRINT*,'          2=adjusting PI to give measured TOC'
	PRINT*,'          3=adjusting decomp.rate to give meas.PI and TOC' 
      PRINT*,'          4=using C accumulation'
      PRINT*,'          5=using fixed proportions of C pools'
      PRINT*,'          6=using TOC and the Hillier solver ',
     & '(an analytical soln of RothC =brium run) to get PI & C pools'
	PRINT*,'          7=adjusting PI to give measured TOC at ',
     & 'equilibrium, and further adjusting PI according to weather'
	PRINT*,'  n = Number of layers (max. 10)'
      PRINT*,'  Depth of top of SOM layer 1 (cm)'
      PRINT*,'  Depth of top of SOM layer 2 (cm)'
      PRINT*,'  ...Depth of top of SOM layer n (cm)'
	PRINT*,'  For this soil under arable:'
c	PRINT*,'   Impermeable layer? (0=No, 1=Yes)'
c     PRINT*,'   Depth of impermeable layer (if present) (cm)'
	PRINT*,'   In SOM layer 1'
	PRINT*,'    C content (kg C / ha)'
	PRINT*,'    Bulk density (g/cm3)'
	PRINT*,'    Soil pH'
	PRINT*,'    % clay by weight'
	PRINT*,'    % silt by weight'
	PRINT*,'    % sand weight'
	PRINT*,'   In SOM layer 2...'
	PRINT*,'   In SOM layer n...'
	PRINT*,'  For this soil under grassland:'
c	PRINT*,'   Impermeable layer? (0=No, 1=Yes)'
c     PRINT*,'   Depth of impermeable layer (if present) (cm)'
	PRINT*,'   In SOM layer 1'
	PRINT*,'    C content (kg C / ha)'
	PRINT*,'    Bulk density (g/cm3)'
	PRINT*,'    Soil pH'
	PRINT*,'    % clay by weight'
	PRINT*,'    % silt by weight'
	PRINT*,'    % sand weight'
	PRINT*,'   In SOM layer 2...'
	PRINT*,'   In SOM layer n...'
	PRINT*,'  For this soil under forestry:'
c	PRINT*,'   Impermeable layer? (0=No, 1=Yes)'
c     PRINT*,'   Depth of impermeable layer (if present) (cm)'
	PRINT*,'   In SOM layer 1'
	PRINT*,'    C content (kg C / ha)'
	PRINT*,'    Bulk density (g/cm3)'
	PRINT*,'    Soil pH'
	PRINT*,'    % clay by weight'
	PRINT*,'    % silt by weight'
	PRINT*,'    % sand weight'
	PRINT*,'   In SOM layer 2...'
	PRINT*,'   In SOM layer n...'
	PRINT*,'  For this soil under natural/seminatural:'
c	PRINT*,'   Impermeable layer? (0=No, 1=Yes)'
c     PRINT*,'   Depth of impermeable layer (if present) (cm)'
	PRINT*,'   In SOM layer 1'
	PRINT*,'    C content (kg C / ha)'
	PRINT*,'    Bulk density (g/cm3)'
	PRINT*,'    Soil pH'
	PRINT*,'    % clay by weight'
	PRINT*,'    % silt by weight'
	PRINT*,'    % sand weight'
	PRINT*,'   In SOM layer 2...'
	PRINT*,'   In SOM layer n...'
      PRINT*,'  Arable plant C input (kg C /ha/yr)'
      PRINT*,'  Grassland plant C input (kg C /ha/yr)'
      PRINT*,'  Forestry plant C input (kg C /ha/yr)'
      PRINT*,'  Natural plant C input (kg C /ha/yr)'
      PRINT*,'  Long term average rainfall (mm/month)'
	PRINT*,'    January'
	PRINT*,'    February'
	PRINT*,'    March'
	PRINT*,'    April'
	PRINT*,'    May'
	PRINT*,'    June'
	PRINT*,'    July'
	PRINT*,'    August'
	PRINT*,'    September'
	PRINT*,'    October'
	PRINT*,'    November'
	PRINT*,'    December'
      PRINT*,'  Long term average temperature (deg.C/month)'
	PRINT*,'    January'
	PRINT*,'    February'
	PRINT*,'    March'
	PRINT*,'    April'
	PRINT*,'    May'
	PRINT*,'    June'
	PRINT*,'    July'
	PRINT*,'    August'
	PRINT*,'    September'
	PRINT*,'    October'
	PRINT*,'    November'
	PRINT*,'    December'
	PRINT*,'  Latitude'
      PRINT*,'  Water table depth at start (cm)'
	PRINT*,'  Drainage class'
      PRINT*,'  C accum.before change (kgC/ha/yr)'
      PRINT*,'  CH4 emission before change (kgC/ha/yr)'
      PRINT*,'  CO2 emission before change (kgC/ha/yr)'
      PRINT*,'  DOC loss before change (kgC/ha/yr)'
      PRINT*,'  Number of growing seasons, n'
      PRINT*,'  For each growing season: Land use' 
	PRINT*,'                            1=arable'
	PRINT*,'                            2=grassland'
	PRINT*,'                            3=forestry'
	PRINT*,'                            4=natural/seminatural'
      PRINT*,'  Name of annual weather file'
	PRINT*,'   (if used). If not use long term average.'
	PRINT*,'   Format: Month,Rain(mm/mon),PET(mm/mon),Temp(deg.C/mon)'
	PRINT*,'   Note: need quotes round file name'
	PRINT*,'======================================================'
	PRINT*,''
      PRINT*,'Enter the name of the file containing input data here:'
	PRINT*,'....'
      READ(*,10)SITE
	PRINT*,'======================================================'
10    FORMAT(A50)
      GOTO 111
113   CONTINUE
C
C Read in data
C
C INMODE: 1=use plant inputs; 2=use toc; 3=use plant inputs and toc; 4=use C accumulation
      READ(42,*)INMODE
	IF(INMODE.EQ.1)THEN			
	  EQMODEL=EQNPP
	  ICMODEL=ICROTHCEQ
	ELSEIF(INMODE.EQ.2)THEN
	  EQMODEL=EQTOC
	  ICMODEL=ICROTHCEQ
	ELSEIF(INMODE.EQ.3)THEN
	  EQMODEL=EQNPPTOC
	  ICMODEL=ICROTHCEQ
	ELSEIF(INMODE.EQ.4)THEN
	  EQMODEL=EQTOC
	  ICMODEL=ICROTHCEQ
	ELSEIF(INMODE.EQ.5)THEN
	  EQMODEL=EQTOC
	  ICMODEL=ICFIXED
	ELSEIF(INMODE.EQ.6)THEN
	  EQMODEL=EQHILLIER
	  ICMODEL=ICROTHCEQ
	ELSEIF(INMODE.EQ.7)THEN
	  EQMODEL=EQTOC
	  ICMODEL=ICROTHCEQ
	  ISDYNPI=ISDYNPI_ON
	ENDIF
	READ(42,*)NSOMLAY(1,1)
	DO 100 ISOMLAY=1,NSOMLAY(1,1)
        READ(42,*)SOMDEPTH(1,1,ISOMLAY)
100   CONTINUE
      DO 200 ILU=1,MAXLU1
	  NSOMLAY(1,ILU)=NSOMLAY(1,1)
	  DOMSOILISIMP(1,ILU)=0
	  DOMSOILIMPDEPTH(1,ILU)=300
C
C Read in soil layers
C
C	  READ(42,*)DOMSOILISIMP(1,ILU)
C	  READ(42,*)DOMSOILIMPDEPTH(1,ILU)
	  DO 300 ISOMLAY=1,NSOMLAY(1,ILU)
	    SOMDEPTH(1,ILU,ISOMLAY)=SOMDEPTH(1,1,ISOMLAY)
          READ(42,*)DOMSOILC(1,ILU,ISOMLAY)
          READ(42,*)DOMSOILBD(1,ILU,ISOMLAY)
          READ(42,*)DOMSOILPH(1,ILU,ISOMLAY)
          READ(42,*)DOMSOILCLAY(1,ILU,ISOMLAY)
          READ(42,*)DOMSOILSILT(1,ILU,ISOMLAY)
          READ(42,*)DOMSOILSAND(1,ILU,ISOMLAY)
		IF(DOMSOILCLAY(1,ILU,ISOMLAY).LT.0.0001)
     &          DOMSOILCLAY(1,ILU,ISOMLAY)=0.0001
		IF(DOMSOILSILT(1,ILU,ISOMLAY).LT.0.0001)
     &          DOMSOILSILT(1,ILU,ISOMLAY)=0.0001
		IF(DOMSOILSAND(1,ILU,ISOMLAY).LT.0.0001)
     &          DOMSOILSAND(1,ILU,ISOMLAY)=0.0001
300     CONTINUE
200   CONTINUE
      DO 400 ILU=1,MAXLU1
        READ(42,*)TOTPIC(ILU)
400   CONTINUE
      DO 500 IMON=1,12
	  READ(42,*)AVERAIN(IMON)
500   CONTINUE
      DO 600 IMON=1,12
	  READ(42,*)AVETEMP(IMON)
600   CONTINUE
	READ(42,*)GRIDLAT
	READ(42,*)
	READ(42,*)DRAINCLASS
	READ(42,*)CACCUM
	READ(42,*)
	READ(42,*)
	READ(42,*)
	READ(42,*)NXYEARS
C
C Read in LU data and check profile information for each LU included
C
      DO 650 ILU=1,MAXLU1
	  ISCHECKED(ILU)=0
650   CONTINUE
	DO 700 IYEAR=1,NXYEARS
	  READ(42,*)LUSEQ(IYEAR)
        IF(ISCHECKED(LUSEQ(IYEAR)).EQ.0)THEN
          CALL CHECK_PROFILE(LUSEQ(IYEAR),NSOMLAY,
     &                       SOMDEPTH,DOMSOILC,DOMSOILCLAY,
     &                       DOMSOILIMPDEPTH,DOMSOILISIMP,DOMSOILSILT,
     &                       DOMSOILSAND,DOMSOILBD,DOMSOILPH)
	    ISCHECKED(LUSEQ(IYEAR))=1
	  ENDIF
700   CONTINUE
	DO 800 IYEAR=1,NXYEARS
        READ(42,*,END=801)
800   CONTINUE
801   CONTINUE
C
C If NPP entered as zero, calculate NPP values using the MIAMI model
C
      DO 900 ILU=1,MAXLU1
	  IF(TOTPIC(ILU).LT.0.0001.AND.TOTPIC(ILU).GT.-0.0001)THEN
	    CALL MIAMI(TOTPIC(ILU),AVERAIN,AVETEMP)
	  ENDIF
900   CONTINUE
C
C Get PET from Thornthwaite equation
C
	CALL GET_PET(GRIDLAT,AVETEMP,AVEPET)
C
C Close data file
C
      CLOSE(42)
C
C Leave SETFILE_SITE
C
      RETURN
      END
C
C--------------------------------------------------------
C
	SUBROUTINE SETWEEKVARS(JSTOP,CLOSSX,CLOSSX15,FYMFERT,FYMFERT15,
     &                       FIXN,THISFERT,THISFERT15,WLEACH,WLEACH15,
     &                       VOLAT,VOLAT15)
C
C Subroutine to reset weekly variables
C      
      IMPLICIT NONE
C
C Variables passed to/from calling subroutine
C
	INTEGER JSTOP
	REAL CLOSSX,CLOSSX15,FYMFERT,FYMFERT15,FIXN,THISFERT,THISFERT15
      REAL WLEACH,WLEACH15,VOLAT,VOLAT15
C
C Rezero weekly variables
C
      IF(JSTOP.EQ.0)THEN
       CLOSSX=0
       CLOSSX15=0
      END IF
	FYMFERT=0
      FYMFERT15=0
      FIXN=0
      THISFERT=0
      THISFERT15=0
      WLEACH=0
      WLEACH15=0
      VOLAT=0
      VOLAT15=0
C
C Leave SETWEEKVARS
C
      END
C
C--------------------------------------------------------
C
      SUBROUTINE SETWTD(SUM_TS,ENTER_WTD,SOILW,WSAT)
C
C Subroutine to get filenames and simulation inputs
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER IL					! Local counter variable
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	REAL TIMESTEP				! Timesteps since beginning of simulation
	REAL WTD					! Water table depth (cm)
	INTEGER WTL					! Layer number of water table depth

C
C Variables passed to/from this subroutine
C
	INTEGER ENTER_WTD			! IN(CALL): Code to use read in water table depth 1 = read in, 0=dont read in
	REAL SOILW(MAXLAYER)		! IN(CALL):Available water (mm/layer)
	INTEGER SUM_TS				! IN(CALL):Total number of timesteps passed
	REAL WSAT(MAXLAYER)			! IN(CALL):Available water at saturation (mm/layer)
C
C Save last timestep and water table depth read in
C
      SAVE
C
C If not coded to read in water table depth, leave subroutine
C
      IF(ENTER_WTD.EQ.0)RETURN
C
C If the last timestep read in is the current timestep...
C
      IF(SUM_TS.EQ.TIMESTEP.OR.SUM_TS.EQ.1)THEN
C
C ...set all layers below WTD to saturated water content
C
        WTL=(WTD*MAXLAYER/MAXDEPTH)+1
	  IF(WTL.GT.MAXLAYER)WTL=MAXLAYER
	  IF(WTL.LT.1)WTL=1
	  DO 100 IL=WTL,MAXLAYER
	    SOILW(IL)=WSAT(IL)
100     CONTINUE
C
C ...read in next water table depth
C
        READ(47,*,ERR=111,END=112)TIMESTEP,WTD
	  GOTO 101
C
C ......error in file reached
C
111     CONTINUE
        WRITE(*,*)'Error in file WTD.DAT'
C
C ......end of WTD.DAT file reached
C
112     CONTINUE
        ENTER_WTD=0
	  CLOSE(47)
C
C .....continue reading
C
101   CONTINUE
	ENDIF
C
C Leave SETWTD
C
      END