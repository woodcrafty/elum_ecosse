C------------------------------------------------------------
C
C Monte Carlo analysis subroutines
C
C developed by Pia Gottschalk for sensitivity analysis within the 
C AfricaNuances project (2005 - 2007)
C incorporated into the main code by Jessica Bellarby
C
C*************************************************************
C MAIN RUN ROUTINES
C*************************************************************
C-------------------------------------------------------------
C
C EXTERNAL SUBROUTINES
C 1. OPENCHAN_MC
C 2. MC_input_prep
C 3. SETFILE_MC
C 4. SETCROPPARAM_MC
C 5. MOREFIX_MC
C 6. INIT_SUNDIAL_CROP_MC
C 7. SETCROPN_MC
C 8. SETC_MC
C 9. TS_RES1_MC
C 10.MC_RESULTS
C
C------------------------------------------------------------
C
C-------------------------------------------------------------
C
      SUBROUTINE ECOSSE_MONTE_CARLO()
C
C Subroutine to run ECOSSE for SITE in MONTE CARLO MODE
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXCROP				! Max.no.of crop types allows
      PARAMETER (MAXCROP=36)
	INTEGER MAXFERT				! Max.no.of fertiliser applications allowed 
	PARAMETER (MAXFERT=5)
	INTEGER MAXORGM				! Max.no.of organic manure applications allowed
	PARAMETER (MAXORGM=52)
	INTEGER MAXSOIL				! Max.no.of soil types
	PARAMETER (MAXSOIL=50)
	INTEGER MAXWEATH			! Max.no.of years allowed
	PARAMETER (MAXWEATH=300)
	INTEGER MAXGROW				! Max.no.of growing seasons allowed
	PARAMETER (MAXGROW=300)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/

	INTEGER IL					! Local counter variable
	INTEGER ILU					! Counter for land use 
	INTEGER MEASLAY				! Layer that soil is measured to
C
C Variables passed to/from other subroutines
C ...Model descriptors
C
      INTEGER MODTYPE				! OUT:Crop model type: 0=SUNDIAL, 1=MAGEC
	INTEGER IS_SUNDIAL			! Crop model type: 0=SUNDIAL
	INTEGER IS_MAGEC			! Crop model type: 1=MAGEC
	DATA IS_SUNDIAL,IS_MAGEC /0,1/
	INTEGER CNMODEL				! OUT:Type of Model for C:N ratio of DPM and 
	                            !     RPM (MAGEC or Foereid)
	INTEGER CNFOEREID			! C:N ratio obtained by method of Foereid
	INTEGER CNMAGEC				! C:N ratio obtained by method of MAGEC
	DATA CNMAGEC,CNFOEREID /1,2/
	INTEGER DMODEL				! OUT:Type denitrification model chosen 
	                            ! (Bradbury or NEMIS)
	INTEGER DBRADBURY			! Bradbury model for denitrification
	INTEGER DNEMIS				! Nemis model for denitrification
	DATA DBRADBURY,DNEMIS /1,2/
	INTEGER DOCMODEL			! OUT:DOC model (on or off)
	INTEGER DOC_OFF				! DOC model off
	INTEGER DOC_ON				! DOC model on
	DATA DOC_ON,DOC_OFF /1,2/
	INTEGER EC_EQRUN			! IN(READ_MODEL):Initialisation using a full ECOSSE equilibrium run (on or off)
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
	INTEGER SPARMODEL			! OUT:Soil parameter model (from file or calc)
	INTEGER SPARFILE			! Soil parameters read in from file
	INTEGER SPARCALC			! Soil parameters calculated from TOC etc
	DATA SPARFILE,SPARCALC /1,2/
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
	INTEGER ITFUNC				! OUT:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER ITFUNC_ROTHC		! ROTHC temperature rate modifier
	INTEGER ITFUNC_HADLEY		! HADLEY temperature rate modifier
      DATA ITFUNC_ROTHC,ITFUNC_HADLEY /0,1/
	INTEGER IMFUNC				! OUT:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IMFUNC_ROTHC		! ROTHC moisture rate modifier
	INTEGER IMFUNC_HADLEY		! HADLEY moisture rate modifier
      DATA IMFUNC_ROTHC,IMFUNC_HADLEY /0,1/
	INTEGER CH4MODEL			! OUT:Methane model (off, Richards or Aitkenhead) 
	INTEGER CH4_OFF				! CH4 model off
	INTEGER CH4_RICHARDS		! Richards CH4 model on		
	INTEGER CH4_AITKENHEAD      ! Aitkenhead CH4 model on
      DATA CH4_OFF,CH4_RICHARDS,CH4_AITKENHEAD /0,1,2/
C
C ...Timing factors
C
	INTEGER FIXEND				! IN:Fixed end? 0=No 1=Yes
	INTEGER I_TS				! IN:No.timesteps since sowing
	INTEGER IANTHES				! IN:No.timesteps from 01/01 to anthesis 
	INTEGER IDATEFC				! IN:Date of field capacity (1=01/01; 2=01/06)
	INTEGER IEND				! IN:Timesteps from sowing date to harvest
	INTEGER IFILE				! IN:Current weather file
	INTEGER IK					! OUT:No.timesteps from prev.harv. to current
	INTEGER IL_TSS				! IN:No.timesteps before harvest when senesces
	INTEGER ILAST_TS			! IN:Last simulated timestep since 01/01/01
	INTEGER IRYEAR				! IN:Current weather year number
	INTEGER IS_TS				! IN:Crop counter
      INTEGER ISTART_TS			! IN:First simulated timestep since 01/01/01 
	INTEGER ISTYR				! IN:First year in simulation (eg. 2001)
	INTEGER IYEAR				! IN:Current growing season number
	INTEGER JSTOP				! IN:Code to stop simulation 0=No 1=Yes
	INTEGER K_TS				! IN:No.timesteps since start of year(UNUSED?)
	INTEGER L_TSS(0:MAXCROP)	! IN:No.timesteps before harvest
	INTEGER LHARV				! IN:No.timesteps from 01/01 to prev. harvest
	INTEGER MEND				! IN:Number of timesteps from prev.harv.-harv.
	INTEGER N_TS				! IN:No.timesteps since start of year(UNUSED?)
	INTEGER NDATE				! IN:Not used?
	INTEGER NF					! IN:No.fertiliser application
	INTEGER NYEARS				! IN:No.weather years included in the sim.
	INTEGER N_STEPS				! IN:No.timesteps in a year
	REAL SECONDS				! IN:Number of seconds in one timestep
	INTEGER SUM_TS				! IN:Total number of timesteps passed
C
C ...Soil factors
C
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
      REAL ATM					! IN:Atmospheric N input (kgN/ha/timestep)
	REAL ATM15					! IN:Atmospheric N15 input (kgN15/ha/time)
      REAL BALANCE(20)			! IN:C Results
	REAL BCARB0(MAXLAYER)		! IN:C in soil biomass (kgC/ha/layer)
      REAL BNIT0(MAXLAYER)		! IN:N in soil biomass (kgN/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN:N15 in soil humus (kgN15/ha/layer)
	REAL BSTART					! IN:Initial N in biomass pool (kgN/ha)
	REAL BSTART15				! IN:Initial N in biomass pool (kgN15/ha)
	REAL BULKDENS(MAXLAYER)		! IN:Bulk density of the layer (g/cm3)
	REAL C_TS					! IN:Litter C input in this timestep (kgC/ha)
	REAL CH4(MAXLAYER)			! IN:CH4 emitted (kgC/ha/layer) (Aitkenhead CH4 model)
	REAL CH4TOAIR				! OUT:CH4 released to atmosphere (kgC/ha) (Aitkenhead CH4 model)
	REAL CH4_PROD(MAXLAYER)  	! OUT: CH4 produced in each soil layer [kgC/ha/layer] (Richards CH4 model)
	REAL CH4_SOIL_OX			! OUT: CH4 oxidised from methane produced in the soil [kgC/ha] (Richards CH4 model)
      REAL CH4_ATMOS_OX			! OUT: CH4 oxidised from atmospheric CH4 [kgC/ha] (Richards CH4 model)
      REAL CH4_FLUX				! Net CH4 flux from soil column [kgC/ha] (Richards CH4 model)
      REAL CLAY(MAXLAYER)			! IN:Clay content of the layer (%)
	REAL CO2(MAXLAYER)			! IN:CO2 emitted (kgC/ha/layer)
	REAL CONC					! IN:Conc.N leached (kgN/ha/mm/layer)
	REAL CONC15					! IN:Conc.N15 leached (kgN/ha/mm/layer)
	REAL CRIT(MAXSOIL,MAXLAYER)	! IN:Critical N content (kgN/ha/layer)
	REAL CTOT					! IN:Total litter C input (kgC/ha)
	REAL DENITRIFN(MAXLAYER)	! IN/OUT: N lost by denitrification (kgN/ha/layer
	REAL DENIT					! IN:N lost by denitrification (kgN/ha)
	REAL DN15					! IN:N15 lost by denitrification (kgN15/ha)
	REAL DNIT(MAXLAYER)			! IN:Net Mineralised N (kgN/ha/lay/timestep)
	REAL DNIT15(MAXLAYER)		! IN:Net Mineralised N15 (kgN15/ha/lay/time)
      REAL DPM_RPM_RATIO			! IN/OUT:Ratio of DPM:RPM. If set to 
								!        zero, this will be worked out using a call to SET_DPMRPMRATIO
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
      REAL DPMCTON(MAXLAYER)		! OUT:DPM C:N ratio 
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL DRAINW(MAXLAYER)		! IN:Water drained from this layer (mm/layer)
      REAL FIELDCAP(MAXLAYER)		! IN/OUT: Total water in the layer at 
	                            !		field capacity (mm / layer)
	REAL FLOWPROP				! IN:Proportion of flow needed to achieve
							    !	     water table at depth WTABLE
	REAL GNN2O					! IN:N2O lost by nitrification (kgN/ha)
	REAL GNNO					! IN:NO lost by nitrification (kgN/ha)
	REAL GPNN2O					! IN:N2O lost by part.nitrification (kgN/ha)
	REAL G15NN2O				! IN:15N2O lost by nitrification (kgN15/ha)
	REAL G15NNO					! IN:15NO lost by nitrification (kgN15/ha)			
	REAL G15PNN2O				! IN:15N2O lost by part.nitrif.(kgN15/ha)			
	REAL GDN2					! IN:N2 lost by denitrification (kgN/ha)
	REAL GDN2O					! IN:N2O lost by denitrification (kgN/ha)
	REAL G15DN2					! IN:15N2 lost by denitrification (kgN15/ha)
	REAL G15DN2O				! IN:15N2O lost by denitrification (kgN15/ha)					
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN:N in soil humus (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN:N15 in soil biomass (kgN15/ha/layer)
	REAL HZ1(MAXLAYER)			! IN:N:C of input PM for steady state 
	INTEGER IAWC				! IN:Water movement code number
	INTEGER IAWCJ				! IN:Water movement code number
  	REAL ICFACTOR(MAXLAYER)		! IN:Adjustment in rate needed to	
								!    achieve measured NPP and TOC	
	INTEGER IDRAINJ				! IN:Water movement code number 
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
	INTEGER IROCK				! IN:No.layers: 1=50cm, 2=100cm , 3=150cm
	INTEGER IROCKJ				! IN:No.layers: 1=50cm, 2=100cm , 3=150cm
	REAL LEACHDOC(MAXLAYER)		! IN/OUT:DOC leached from the layer (kgC/ha/layer/timestep)
	REAL LEACHDON(MAXLAYER)		! IN/OUT:DON leached from the layer (kgN/ha/layer/timestep)
	INTEGER LUCODE				! IN/OUT:LU code for equilibrium run 
	REAL DRRAT					! IN/OUT:DPM:RPM ratio
	REAL NITRIFN(MAXLAYER)		! IN/OUT:Nitrification from this layer (kg/ha/layer/timestep)
	INTEGER NSOIL				! IN:Soil code number
	INTEGER NSOILJ				! IN:Soil code number
								!(CHECK USAGE OF IDRAINJ AND IAWCJ)
	REAL ORGC					! IN:Total org.C input (kgC/ha) 
	REAL ORGN					! IN:Total org.N input (kgN/ha) 
	REAL PI_ANN					! IN:Annual plant input (kgC/ha/year)
      REAL PI_C(MAXLAYER)			! IN/OUT: Plant input C to soil (kgC/ha/step)
	REAL PI_CEQ_MON(12,MAXLAYER)! IN/OUT: Equilibrium plant C input each month
								!         in each layer (kgC/ha/month/layer)
      REAL PI_N(MAXLAYER)			! IN/OUT: Plant input N to soil (kgC/ha/step)
      REAL PI_N15(MAXLAYER)		! IN/OUT: Plant input N15 to soil(kgC/ha/step)
								!(1=arable,2=grassland,3=forestry,4=semi-natural)
	REAL PLANTUP(MAXLAYER)		! IN(RUN2_GIS_CROP)/OUT(TEST1_RES):Plant uptake from each soil layer
	REAL REFIL(MAXLAYER)		! IN:Water needed to fill to FC (mm/lay)
	REAL RNIN					! IN:Litter N input in this timestep (kgN/ha)
	REAL RNIN15					! IN:Litter N15 input in timestep (kgN15/ha)
	REAL RORGN					! IN:Total org.N input (kgN/ha) DONT PASS?
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
      REAL RPMCTON(MAXLAYER)		! OUT:RPM C:N ratio 
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL RSTART					! IN:Initial N in debris pool (kgN/ha)
	REAL RSTART15				! IN:Initial N in debris pool (kgN15/ha)
	REAL SAND(MAXLAYER)			! IN:Sand content of the layer (%)
	REAL SATWATCONT(MAXLAYER)	! IN:Total water content at saturation (mm/layer)
	REAL SC						! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SILT(MAXLAYER)			! IN:Silt content of the layer (%)
	REAL SLEA15(MAXLAYER)		! IN:Nitrate-N15 leached (kgN15/ha)
	REAL SLEACH(MAXLAYER)		! IN:Nitrate-N leached (kgN/ha)
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOILPH(MAXLAYER)		! IN/OUT:pH of soil in this layer
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL SORGC					! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SORGN					! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SORGN15				! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SRNIN					! IN:Total litter N input (kgN/ha)
	REAL SX						! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SXORGC					! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SXORGN					! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SXORGN15				! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL T15					! IN:Net Mineralised N (kgN/ha/layer)
	REAL TAM(MAXLAYER)			! IN:15N/N in the ammonium in the layer
	REAL TC						! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL TCINP					! IN:Total org.C input (kgC/ha) DONT PASS?
	REAL TDNIT					! IN:Net Mineralised N (kgN/ha/layer)
	REAL TFYMC					! IN:Total FYM C input (kgC/ha) DONT PASS?
	REAL TIN(MAXLAYER)			! IN:15N/N in the nitrate in the layer
	REAL TOC(MAXLAYER)			! IN:Total organic C (kgC/ha/layer)
      REAL TORGC					! IN:Total org.C input (kgC/ha) DONT PASS?
	REAL TORGN					! IN:Total org.N input (kgN/ha) DONT PASS?
	REAL TOTAST					! IN:Total ammonium in soil (kgN/ha)
	REAL TOTNST					! IN:Total nitrate in soil (kgN/ha)
	REAL TOTAST15				! IN:Total ammonium in soil (kgN15/ha)
	REAL TOTNST15				! IN:Total nitrate in soil (kgN15/ha)
	REAL TRNIN					! IN:Total litter N input (kgN/ha)
	REAL TXORGC					! IN:Total org.C input (kgC/ha) DONT PASS?
	REAL TXORGN					! IN:Total org.N input (kgN/ha) DONT PASS?
	REAL VOLAT					! IN:N lost by volatilisation (kgN/ha)
	REAL VOLAT15				! IN:N15 lost by volatilisation (kgN15/ha)
	REAL VOLATN(MAXLAYER)		! IN/OUT:Volatilisation from this layer (kg/ha/layer/timestep)
      REAL WATCONT(MAXLAYER)		! IN: Total water content of the soil layer (mm/layer)
	REAL WILTPOINT(MAXLAYER)	! IN:Total water conent of the soil layer (mm/layer)
	REAL WLEACH					! IN:N leached (kgN/ha)
	REAL WLEACH15				! IN:N15 leached (kgN15/ha)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
	REAL WTABLE					! IN:Water table depth in cm
	REAL XORGC					! IN:Total org.C input (kgC/ha) 
	REAL XORGN					! IN:Total org.N input (kgN/ha) 
C
C ... Crop factors
C
	REAL CACT					! IN:Crop N uptake (kgN/ha/timestep)
	REAL CACT15					! IN:Crop N15 uptake (kgN15/ha/timestep)
	REAL CACTOT					! IN:N taken up by crop (kgN/ha)
	REAL CATOT15				! IN:N15 taken up by crop (kgN15/ha)
	REAL CFACT(0:MAXCROP)		! IN:Crop parameter used to estimate yield
	REAL CLOSSX					! IN:N lost by senescence (kgN/ha) 
	REAL CLOSSX15				! IN:N15 lost by senescence (kgN15/ha) 	
	REAL CUPTN					! IN:Crop N uptake (kgN/ha/timestep)
	REAL EXYLD					! IN:Yield of current crop (t/ha)
      REAL EXYLDJ(0:MAXGROW)		! IN:Yield of current crop (t/ha)
								! (0->calculate within model)
	REAL FIXN					! IN:N fixed from atmosphere (kgN/ha)
      INTEGER ICOVER				! OUT:Crop cover 1=Covered 0=Bare
	INTEGER MCROP				! IN:Crop type code number
	INTEGER NGROW				! IN:No.growing seasons simulated
	INTEGER NXYEARS				! IN:No.growing seasons simulated
	INTEGER ICROP				! IN:Current crop code
	INTEGER ICROPJ(0:MAXGROW)	! IN:Crop codes
	INTEGER ISOWN				! IN:Timesteps from 01/01/01 to sowing date 
	INTEGER ISOWNJ(0:MAXGROW)	! IN:Timesteps from 01/01/01 to sowing date 
	INTEGER	ISTHARV				! IN:Timesteps from 01/01/01 to first harvest
	INTEGER IHARV				! IN:Timesteps from 01/01/01 to harvest date 
	INTEGER IHARVJ(0:MAXGROW)	! IN:Timesteps from 01/01/01 to harvest date 
	INTEGER INSTRAW				! IN:Straw incorporated? 0=No 1=Yes
	INTEGER LCROP				! IN:Previous crop type
	INTEGER NFERT				! IN:No.fertiliser applications to this crop
      INTEGER NFERTJ(0:MAXGROW)	! IN:No.fertiliser applications to crops
	INTEGER NORGM				! IN:No.manure applications to this crop
	INTEGER NORGMJ(0:MAXGROW)	! IN:No.manure applications to crops
      INTEGER NRESJ(0:MAXGROW)	! IN:Straw incorporated? 0=No 1=Yes
	INTEGER NSOW				! IN:Timesteps from prev.harv. to sowing date
	REAL PREYLD					! IN:Yield of previous crop (t/ha)
	REAL OCROPN					! IN:N uptake of crop (kgN/ha) 
	REAL OCROPNJ(0:MAXGROW)		! IN:N uptake of crop (kgN/ha) 
	REAL ROOT					! IN:Rooting depth according to restriction (cm)
	REAL TACTOT					! IN:Total crop N uptake (kgN/ha)
	REAL UT(3,0:MAXCROP)		! IN:Crop parameter used to estimate yield
	REAL YLD					! IN:Yield (t/ha) DONT PASS?
	REAL WR						! IN:Root requirement (kgN/ha)
C
C ... Fertiliser factors
C
	REAL FERT(MAXFERT)			! IN:Amount of fertiliser applied (kgN/ha)
	REAL FERTJ(0:MAXGROW,MAXFERT)	! IN:Amount of fertiliser applied (kgN/ha)
	REAL FERTNF					! INOUT:Amount of fertiliser applied (kgN/ha)
	REAL FYMFERT				! IN:N input by fertiliser & FYM (kgN/ha)
	REAL FYMFERT15				! IN:N15 input by fertiliser & FYM (kgN15/ha)
	INTEGER IFERT(MAXFERT)		! IN:No.timesteps to split fert.application
      INTEGER IFERTJ(0:MAXGROW,MAXFERT)	! IN:No.timesteps to fert.application
	INTEGER ILAB(MAXFERT)		! IN:Labelling on fertiliser? 0=No 1=Yes
	INTEGER ILABJ(0:MAXGROW,MAXFERT)	! IN:Lab.on fertiliser? 0=No 1=Yes
	INTEGER ILABNF				! INOUT:Labelling on fertiliser? 0=No 1=Yes
	INTEGER IVOL(MAXFERT)		! IN:N volatilised from fertiliser (kgN/ha)
      INTEGER IVOLJ(0:MAXGROW,MAXFERT)	! IN:Does fert.contain ammonium salts 
										! other than ammonium sulphate? (0=No; 1=Yes)
	INTEGER IVOLNF				! INOUT:N volatilised from fertiliser (kgN/ha)
	REAL JFERTJ(0:MAXGROW,MAXFERT,3)	! IN:Prop.NO3,NH4,urea in fertiliser
	REAL TFERT(MAXFERT,3)		! IN:Prop.NO3,NH4,urea in fertiliser
	REAL TFERTNF(3)				! INOUT:Prop.NO3,NH4,urea in fertiliser
	REAL THISFERT				! IN:N input by fertiliser (kgN/ha)
	REAL THISFERT15				! IN:N15 input by fertiliser (kgN15/ha)
C
C ... Manure factors
C     
	INTEGER IORGM(MAXORGM)		! IN:No.timesteps to split manure application
	INTEGER JORGM(MAXORGM)		! IN:Type of manure application
	INTEGER IOLAB(MAXORGM)		! IN:Labelling on manure? 0=No 1=Yes
	INTEGER IORGMJ(0:MAXGROW,MAXORGM)	! IN:No.timesteps to manure applic.
	INTEGER JORGMJ(0:MAXGROW,MAXORGM)	! IN:Type of manure application
      INTEGER IOLABJ(0:MAXGROW,MAXORGM)	! IN:Labelling on manure? 0=No 1=Yes
	REAL ORGMA(MAXORGM)			! IN:Amount of manure applied (t/ha)
      REAL ORGMJ(0:MAXGROW,MAXORGM)	! IN:Amount of manure applied (t/ha)
	INTEGER IOLABNF				! IN:Labelling on manure? 0=No 1=Yes
	INTEGER JORGMNF				! IN:Type of manure application
	REAL ORGMANF				! IN:Amount of manure applied (t/ha)
C
C ...Weather factors
C
      REAL AVEPET(12)				! Long term average PET (mm/month)
      REAL AVERAIN(12)			! Long term average rainfall(mm/month)
      REAL AVETEMP(12)			! Long term average monthly average 
	REAL LTA_AWC(12,MAXLAYER)	! IN:Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)	! IN:Average air temp this month (deg.C)
	CHARACTER*100 FIXFILE(MAXWEATH)	! IN: Weather files
	REAL LAT					! IN: Latitude
	REAL RAIN					! IN: Rainfall (mm/timestep)
	REAL EVAP					! IN: Potential evap. (mm/timestep)
	REAL EVAPW					! IN: Potential evap. (mm/timestep)
								!	CONSIDER COMBINING EVAP AND EVAPW
	REAL AIRTEMP				! IN: Air temperature (deg.C/timestep)	
	REAL SOILTEMP(MAXLAYER)		! IN: Soil temperature (deg.C/timestep)
	REAL RDD					! IN: Weather data used by MAGEC
	REAL TMMN					! IN: Weather data used by MAGEC
	REAL TMMX					! IN: Weather data used by MAGEC
	REAL VP						! IN: Weather data used by MAGEC
	REAL WN						! IN: Weather data used by MAGEC
	REAL DTR					! IN: Weather data used by MAGEC
	REAL DVP					! IN: Weather data used by MAGEC
	REAL TMAX					! IN: Weather data used by MAGEC
	REAL TMIN					! IN: Weather data used by MAGEC
	REAL DAVTMP					! IN: Weather data used by MAGEC
	REAL NAVTMP					! IN: Weather data used by MAGEC
	REAL DTRJM2					! IN: Weather data used by MAGEC
	REAL TAVS					! IN: Weather data used by MAGEC
	REAL BBRADS					! IN: Weather data used by MAGEC
	REAL SVPS					! IN: Weather data used by MAGEC
	REAL SLOPES					! IN: Weather data used by MAGEC
	REAL RLWNS					! IN: Weather data used by MAGEC
	REAL NRADS					! IN: Weather data used by MAGEC		
	REAL PENMRS					! IN: Weather data used by MAGEC
	REAL WDF					! IN: Weather data used by MAGEC
	REAL PENMD					! IN: Weather data used by MAGEC
	REAL PE						! IN: Weather data used by MAGEC
	REAL AE						! IN: Weather data used by MAGEC
	REAL DDAYS					! IN: Weather data used by MAGEC
	INTEGER IDAG					! IN: Weather data used by MAGEC
	REAL CRITMIN				! IN: Data used by MAGEC
	REAL EVAPO					! IN: Data used by MAGEC
	REAL SEEDN_S				! IN: Data used by MAGEC
	REAL SEED(0:MAXCROP)					! IN: Data used by MAGEC
	REAL SEEDIN					! IN: Data used by MAGEC
C
C ...Dissolved organic matter factors
C
	REAL CO2FROMDOC(MAXLAYER)	! IN:CO2 output from DOC in each layer(kgC/ha)
	REAL MOBDOC(MAXLAYER)		! IN:Mobile DOC in each layer (kgC/ha)
	REAL MOBDON(MAXLAYER)		! IN:Mobile DON in each layer (kgN/ha)
C
C
C
C ... variables required for MONTE CARLO runs 
C
C ... variables for input (Monte carlo)
C
	INTEGER ios				!internal counter
	INTEGER DOY(365)		!day of the year
	INTEGER I				!internal counter
C
C Weather (Monte carlo)
C
	REAL samples(19)
	REAL precip		!daily precipitation
	REAL temp(365)			!daily temperature
C
C Crops (Monte carlo)
C
	REAL carbonTotRInitYear			!Total carbon return in the initial year
	REAL carbonTotRGrowSeason		!Total carbon return in the growing season
	REAL carbonRetStubInitYear		!Carbon return as stubble in the initial year
	REAL carbonRetStubGrowSeason	!Carbon return as stubble in the growing season
	REAL carbonRetStrawInitYear		!Carbon return as straw in the initial year
	REAL carbonRetStrawGrowSeason	!Carbon return as straw in the grwoing season
	REAL nitrogenReqAboveInitYear	!N requirement aboveground in the initial year
	REAL nitrogenReqBelowInitYear	!N requirement belowground in the inital year
	REAL nitrogenReqAboveGrowS		!N requirement above ground in the growing season
	REAL nitrogenReqBelowGrowS		!N requirement below ground in the growing season
	REAL nitrogenRetStrawInitYear	!N return as straw in the initial year
	REAL nitrogenRetStrawGrowSeason	!N return as straw in the growing season
	REAL DPM_RPM_SOILN_ratio		!Ration of soilN DPM/RPM
C
C	SETup file (Monte carlo)
C
	CHARACTER(LEN=100) text(20)		!?
	INTEGER cropNo					!crop number
	INTEGER plantingdate			!planting date
	INTEGER NRequire				!N requirement
	INTEGER harvestDate				!harvest date
	INTEGER strawFlag				!straw???
	INTEGER noFert					!No fertilizer
	INTEGER noManure				!No manure
	INTEGER dateFert				!date of fertilizer application
	INTEGER ureaFlag				!urea??
	INTEGER n15Flag					!?
	REAL amountFert					!Amount of fertilizer used
	REAL no3Fert					!Fertilizer as nitrate
	REAL nh4Fert					!Fertilizer as ammonium
	REAL ureaFert					!Fertilizer as urea
	REAL expYld						!Expected yield
	REAL amountManure				!Amount of manure
	INTEGER dateManure				!Date of manure application
	INTEGER typeManure				!Type of manure
	INTEGER lableManure				!Labelled manure
	CHARACTER CLIMATEFILE			!Name of climatefile
	CHARACTER SETUP					!setup file
	CHARACTER CROP_MC*50				!crop parameter file for Monte Carlo
	CHARACTER TEMPFIL				!tempfileC			
C
C...additional variables for output (Monte carlo)
C
	REAL SUM_FYMNH4			!Sum of farm manure additions
	REAL SUM_FERT			!Sum of fertilizer additions
	REAL CREQN				!crop requirement for N
	REAL TCNEED				!
	REAL ACTUPT				!Actual plant uptake
	REAL SUM_ACTUPT			!Sum of actual plant uptake
	REAL SumDPMDecomp		!Sum of carbon loss in the DPM pool
	REAL SumRPMDecomp		!Sum of carbon loss in the RPM pool
	REAL SumBDecomp			!Sum of carbon loss in the Humus pool
	REAL SumHDecomp			!Sum of carbon loss in the Biomass pool
	REAL DPMCARB(MAXLAYER)
	REAL RPMCARB(MAXLAYER)
	REAL BCARB(MAXLAYER)
	REAL HCARB(MAXLAYER)
	REAL DPMCARB_BEGINN		!Carbon in the DPM pool at NSOW
	REAL RPMCARB_BEGINN		!Carbon in the RPM pool at NSOW
	REAL HCARB_BEGINN		!Carbon in the Humus pool at NSOW
      REAL BCARB_BEGINN		!Carbon in the biomass pool at NSOW
	REAL TOTCARB_BEGINN		!Total carbon at NSOW
	REAL SOILN_BEGINN		!Soil nitrate N at NSOW
	REAL AMMN_BEGINN		!Soil ammonium N at NSOW
	REAL TD_NIT				!
	REAL SUM_TDNIT			!
	REAL TOTNIT_BEGINN		!Total nitrogen at NSOW
	REAL DPMNIT_BEGINN		!Nitrogen in the DPM pool at NSOW
	REAL RPMNIT_BEGINN		!Nitrogen in the RPM pool at NSOW
	REAL HNIT_BEGINN		!Nitrogen in the Humus pool at NSOW
	REAL BNIT_BEGINN		!Nitrogen in the Biomass pool at NSOW
	REAL HCARB_CHANGE		!Relative change of carbon in the Humus pool during the season
	REAL BCARB_CHANGE		!Relative change of carbon in the biomass pool during the season
	REAL DPMCARB_CHANGE		!Relative change of carbon in the DPM pool during the season
	REAL RPMCARB_CHANGE		!Relative change of carbon in the RPM pool during the season
	REAL TOTCARB_CHANGE		!Relative change of total carbon during the season
	REAL SUM_DENIT			!Sum of denitrification
	REAL SUM_LEACH			!Sum of leaching
	INTEGER TIMESTEP
	REAL nitrogenReqBelow
      REAL nitrogenReqAbove,nitrogenRetStraw
	CHARACTER CLIMATE_PREP*50
	CHARACTER SETUP_PREP*50
C
	INTEGER delayPlanting
C
C Variables required for SWAT water routines
C
	REAL PLANTBIOMASS		! Plant biomass at harvest [ kg ha-1]
	REAL PLANTHEIGHT		! Plant height at harvest [cm]
	REAL LAI				! Plant LAI at harvest
	REAL canMAX				! maximum amount of water that can be trapped in the canopy when the canopy is fully developed (mm H2O)
	CHARACTER*100 FIXFILEaVP(MAXWEATH)	! IN: Weather files (actual vapour pressure) (kPa)
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
	INTEGER WATERMODEL
	INTEGER SUNDIAL_WATER
	INTEGER SWAT_WATER
	DATA SUNDIAL_WATER,SWAT_WATER /0,1/
C
C Setting the depth of profile for the MC_analysis
C
	INTEGER MC_DEPTH
	REAL TRATEM					! Temperature rate modifier
	REAL WRATEDM					! Moisture rate modifier

C
C Data Statements
C
      DATA IYEAR/0/

C RUN1_SUNDIAL_CROP
C Variables local to this subroutine
C
	INTEGER IROCKS(0:MAXCROP),ISTART,ISAVE
	REAL CINP,ANINP,ROOTS,TOTC,BREQN 
      REAL CAO(5,0:MAXCROP),UR(3,0:MAXCROP)
	REAL INC(2,0:MAXCROP),CRATE(0:MAXCROP),ANRATE(0:MAXCROP)
      REAL RRG(0:MAXCROP),RRX(0:MAXCROP)
      REAL CSC(4,0:MAXCROP),STRAW(4,0:MAXCROP)
C
C Variables passed to/from calling subroutine
C
	INTEGER N_REMAIN
C
C************************************************
C Set model type FOR SITE SPECIFIC RUNS here!!!!
C 1. Denitrification model, DMODEL: DNEMIS = NEMIS model; DBRADBURY = Bradbury model
C    
	DMODEL=DBRADBURY
C
C 2. Initialisation of SOM pools, ICMODEL and INMODEL
C
	INMODEL=INPASSCN
      ICMODEL=ICROTHCEQ
C
C 3. DOC model, DOCMODEL: DOC_OFF	= DOC model off; DOC_ON = DOC model on
C
	DOCMODEL=DOC_OFF
C
C 4. CN model, CNMODEL: CNFOEREID = C:N ratio obtained by method of Foereid; C:N ratio obtained by method of MAGEC
C
	CNMODEL=CNFOEREID
C
C 5. Soil parameter model, SPARMODEL: SPARFILE = Soil parameters read in from file; SPARCALC = Soil parameters calculated from TOC etc
C
	SPARMODEL=SPARFILE
C
C 6. Type of equilibrium run: EQNPP = Model initialised using plant inputs from measured NPP;
C							EQTOC = Model initialised using measured TOC;
C							EQNPPTOC = Model initialised using both plant inputs.
C
	EQMODEL=EQTOC
C
C 7. Calculation of moisture rate modifiers, IMFUNC_ROTHC=0, IMFUNC_HADLEY=1
C    Calculation of temp.rate modifiers, ITFUNC_ROTHC=0, ITFUNC_HADLEY=1
C
      IMFUNC=IMFUNC_ROTHC
      ITFUNC=ITFUNC_ROTHC
C
C 8. CH4 model, CH4MODEL: CH4_OFF = CH4 model off; CH4_RICHARDS = Richards CH4 model on, CH4_AITKENHEAD = Aitkenhead CH4 model on
C
      CH4MODEL=CH4_OFF
C
C 9. Water model, WATERMODEL: SUNDIAL_WATER = original Sundial water routines
C							SWAT_WATER = water rountines from SWAT implemented
C
C This is still in development and currently only works for Tropical maize!
C
	WATERMODEL=SWAT_WATER
C
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
	CALL OPENCHAN_MC(CLIMATE_PREP,SETUP_PREP)
C
C 2. Read in mode of model use
C
      CALL READ_MODEL(DMODEL,ICMODEL,INMODEL,CNMODEL,DOCMODEL,SPARMODEL,
     &  		        EQMODEL,IMFUNC,ITFUNC,CH4MODEL,EC_EQRUN)
	EC_EQRUN=0
C
C
C - prepare input files for Monte Carlo
C
	CALL MC_input_prep (delayPlanting,WATERMODEL,SWAT_WATER)
C
C 3. Get filenames and simulation inputs
C
      CALL SETFILE_MC(ATM,IDATEFC,NSOILJ,IDRAINJ,IROCKJ,LCROP,
     &                   PREYLD,SECONDS,MODTYPE,N_STEPS,
     &                   IAWCJ,NYEARS,
     &                   ISTHARV,ISTYR,ISTART_TS,ILAST_TS,FIXEND,LAT,
     &                   FIXFILE,NGROW,NXYEARS,ICROPJ,ISOWNJ,OCROPNJ,
     &                   IHARVJ,EXYLDJ,NRESJ,NFERTJ,NORGMJ,FERTJ,IFERTJ,
     &                   JFERTJ,IVOLJ,ILABJ,ORGMJ,IORGMJ,JORGMJ,
     &                   IOLABJ,WTABLE,
C Variable required for SWAT water
     &					FIXFILEavp,WATERMODEL,SWAT_WATER,z)
C
C
C 4. Read in TOC, IOM, CLAY, SILT, SAND, BULKDENS, NPP
C
      IF(SPARMODEL.EQ.SPARCALC)THEN
	  PRINT*,'----------------------------------------------'
	  PRINT*,'SPARCALC CANNOT BE DONE FOR MONTE CARLO!!!!'
	  PRINT*,'CHANGE SETTING AND START AGAIN!'
	  PRINT*,'----------------------------------------------'
      ENDIF

      IF(EQMODEL.EQ.EQNPP.OR.EQMODEL.EQ.EQNPPTOC)THEN
	  PRINT*,'EQNPP AND EQNPPTOC CANNOT BE DONE IN MONTE CARLO!'
	  PRINT*,'CHANGE SETTING AND START AGAIN!'
	  PRINT*,'----------------------------------------------'
	ELSE
	  PI_ANN=1000
	ENDIF
C
C For each growing season... 
C        
      DO 2 IYEAR=0,NXYEARS

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
        CALL INIT_SUNDIAL_CROP_MC(IYEAR, LCROP, ICROP, 
     &                             INSTRAW, IEND, ISTHARV,   
     &                             ISOWN, MEND, NSOW, JSTOP, L_TSS,
     &                             YLD, PREYLD, EXYLD,  
     &                             SORGC,SORGN, SXORGC, SXORGN,  
     &                             HZ1, TORGC, TORGN,TXORGC, TXORGN, 
     &                             OCROPN, RORGN,DDAYS,
     &                             ORGC, ORGN, XORGC, XORGN,
     &                             NXYEARS,NDATE,FIXEND,LHARV,IHARV,
     &                             NORGM,IORGM,NFERT,IFERT,IANTHES,
     &                             SECONDS,CREQN,
     &			nitrogenReqAboveGrowS,nitrogenReqBelowGrowS,
     &			nitrogenRetStrawGrowSeason,DPM_RPM_SOILN_ratio,
C introduced for SWAT water routines
     &				             PLANTBIOMASS,PLANTHEIGHT,LAI,C1,T1,
     &							 canMax,CONVER_F,RRG,IROCKS)
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
     &                          FIXFILE,FIXFILEaVP,MODTYPE,SOILTEMP,aVP)
		 CALL AVEMET(AVERAIN,AVEPET,AVETEMP)
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
           CALL INIT_SUNDIAL_WATER(SOILW,ROOT,IROCK,SUM_TS,IDATEFC,
     &                              IANTHES,IRYEAR,IFILE,WMAX,K_TS,RAIN,
     &                              EVAP,AIRTEMP,RDD,TMMN,TMMX,VP,WN,
     &                              IAWC,NSOIL,MODTYPE,SPARMODEL,
     &                              TOC,CLAY,SAND,SILT,BULKDENS,WSAT,
     &                              WTABLE,
     &                              FLOWPROP,AVERAIN,AVEPET,AVETEMP,
     &                              IYEAR,ISTHARV,ISOWNJ,N_STEPS,
     &                              WILTPOINT,SATWATCONT,SECONDS)
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
	  ENDIF
C
C For each week till end of growing season...
C
        DO 5 IK=1,MEND
C
C Initialize THISFERT=current weeks fertilizer addition
C            WLEACH=current weeks water leaching
C      
          CALL SETWEEKVARS(JSTOP,CLOSSX,CLOSSX15,FYMFERT,FYMFERT15,
     &                     FIXN,THISFERT,THISFERT15,WLEACH,WLEACH15,
     &                     VOLAT,VOLAT15)
C
C Get weather data
C
C	Required for SWAT water
C
		IF(WATERMODEL.EQ.SWAT_WATER)THEN
          CALL GETWEATHER_SWAT(N_TS,N_STEPS,IFILE,IRYEAR,IK,MEND,SUM_TS,
     &                    RAIN,EVAPW,AIRTEMP,RDD,TMMN,TMMX,VP,WN,
     &                    FIXFILE,FIXFILEaVP,MODTYPE,SOILTEMP,aVP,
C required for Monte carlo - need to pass on original value ofprecipiation
     &						precip)
		ELSE
           CALL GETWEATHER(N_TS,N_STEPS,IFILE,IRYEAR,IK,MEND,SUM_TS,
     &                    RAIN,EVAPW,AIRTEMP,RDD,TMMN,TMMX,VP,WN,
     &                    FIXFILE,MODTYPE,SOILTEMP,
C required for Monte carlo - need to pass on original value ofprecipiation
     &						precip)
		ENDIF
C
C If current week is first in the growing season set results
C
          CALL SETRES(BALANCE,CO2,SX,SXORGN)
          IS_SUNDIAL=0
	    IS_MAGEC=1
C
C SUNDIAL Crop Model...
C
          IF(MODTYPE.EQ.IS_SUNDIAL)THEN
C
C ...Calculate crop C and N returns and N offtake 
C
	      CALL RUN1_SUNDIAL_CROP_MC(IYEAR,IK,MEND,IS_TS,JSTOP,
     &                             IEND,IANTHES,I_TS,
     &                             NSOW,ISTHARV,INSTRAW,ICROP,LCROP,
     &                             ISOWN,OCROPN,YLD,PREYLD,EXYLD,C_TS,
     &                             SXORGN,SORGN,TORGC,SORGC,RORGN,ORGC,
     &                             SXORGN15,SXORGC,XORGC,XORGN,ORGN,
     &                             TCINP,TRNIN,RNIN,RNIN15,
     &                             WR,CTOT,
     &                             CACTOT,CATOT15,ICOVER,
C required for Monte carlo
     &				    nitrogenReqBelowGrowS,
     &					nitrogenReqAboveGrowS,
     &					nitrogenRetStrawGrowSeason,CREQN,
C required for SWAT water
     &				             PLANTBIOMASS,PLANTHEIGHT,LAI,
     &							 canMax)      
     
      
	    ELSEIF(MODTYPE.EQ.IS_MAGEC)THEN
C
C ...Calculate crop C and N returns 
C
	      IF((IYEAR.EQ.0).OR.((IYEAR.EQ.1).AND.(IK.EQ.1)))then
	        CALL RUN1_MAGEC_CROP(I_TS,IEND,IYEAR,IK,SECONDS,
     &                  CACTOT,TRNIN,RNIN,RNIN15,CATOT15,
     &                  ICROP,SXORGN,SXORGN15,C_TS,SXORGC,TCINP,
     &                  XORGC,XORGN,ORGC,ORGN,
     &                  ISOWN,IANTHES,NSOW,ICOVER,0,0)
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
C Run SUNDIAL Soil C and N routines
C
          CALL PARTITION_PLANTINPUT(C_TS,RNIN,RNIN15,PI_C,PI_N,PI_N15)
	    LUCODE=1		! Note: LUCODE set to arable because this 
						! version of the model is arable crops only
          IF(INMODEL.EQ.INPASSCN)
     &      CALL GET_CTON_FOEREID(DPMCTON,RPMCTON,LUCODE)
c     &      CALL GET_CTON_FALLOON(DPMCTON,RPMCTON,LUCODE)
	    CALL SET_DPMRPMRATIO(LUCODE,DRRAT)
	    MEASLAY=MAXLAYER
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
     &							       CH4_ATMOS_OX,CH4_FLUX,
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
C Required for Monte Carlo - 
     &					DPMCARB,RPMCARB,BCARB,HCARB,
     &					WRATEDM,TRATEM,SOILPH,MEASLAY)			
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
		  CALL INTERCEPTION_SWAT_WATER(LAI,CURRENTLAI,canMax,RAIN,
     &								canStor)
	    ENDIF		
C
C Leaching routines
C
          CALL DRAIN_SUNDIAL_WATER(RAIN,WMAX,SOILW,DRAINW,REFIL,
     &                             WSAT,FLOWPROP)
C
C Evapotranspiration 
C
	    IF(WATERMODEL.EQ.SWAT_WATER)THEN
C
C Calculation of actual evapotranspiration (SWAT)
C
		  CALL EVAP_SWAT_WATER(AIRTEMP,z,aVP,LAT,N_TS,CURRENTPLHEI,
     &						CURRENTLAI,canStor,SOILW,FIELDCAP,
     &                        WMAX,WILTPOINT,ROOTLAYER,
     &						CURRENTPLBIOM,EVAPW,WATCONT)
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

          CALL PHYSICAL_SUNDIAL_SOILCN(DOCMODEL,SECONDS,
     &                                 NSOIL,SOILTEMP,
     &                                 WMAX,SOILW,DRAINW,REFIL,
     &                                 SLEACH,SLEA15,
     &                                 WLEACH,WLEACH15,CONC,CONC15,
     &                                 SOILN,SOIL15,AMMN,AMMN15,
     &                                 TIN,TAM,MOBDOC,MOBDON,
     &                                 LEACHDOC,LEACHDON)
C
C Record only summarised results for Monte Carlo
C
      	CALL TS_RES1_MC(IYEAR,IRYEAR,NXYEARS,SECONDS, 
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
      
		CALL MC_RESULTS (IYEAR,IK,NSOW,MEND,NSOIL,
     &					MAXLAYER,precip,EVAPW,AIRTEMP,
     &                   SXORGC,TORGC,SORGC,SEEDIN,
     &					SORGN,SXORGN,ORGN,CREQN,
     &                   DPMCARB0,RPMCARB0,BCARB0,HCARB0,
     &					DPMCARB,RPMCARB,BCARB,HCARB,
     &                   DPMNIT0,RPMNIT0,HNIT0,BNIT0,
     &                   SOILN,AMMN,CACT,
     &                    VOLAT,DENIT,ATM,TDNIT,THISFERT,
     &                   FYMFERT,SLEACH,delayPlanting,SOILW,
     &					WRATEDM,TRATEM)
C
C
C If at the fixed end jump out of simulation
C
          IF(SUM_TS.EQ.FIXEND) GOTO 7200
C
C Go back and do calculation for the next week in the growing season
C
   5    CONTINUE
C
C Go back and set up the next crop
C
    2 CONTINUE
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
      END
C
C #######################################################################
C
C         SUBROUTINES
C
C #######################################################################


	SUBROUTINE OPENCHAN_MC (CLIMATE_PREP,SETUP_PREP)
C	
C Subroutine to open channels
C
      IMPLICIT NONE
	integer ios
	CHARACTER SOIL_MC*50
	CHARACTER SETUP_PREP*50
	CHARACTER CLIMATE_PREP*50
C
C Open report files
C
	OPEN(15,FILE='ERROR.MSG',STATUS='UNKNOWN',ERR=444)
C
C Open and get the name of the Simulation Input File from MCNAMES.DAT
C
	OPEN(4,FILE='MC_NAMES.DAT',STATUS='OLD',ERR=444)
      READ(4,*,ERR=444,END=333)SETUP_PREP,CLIMATE_PREP,SOIL_MC
C
C Open files for Monte Carlo input preparation
C
	OPEN(21,FILE=SETUP_PREP,STATUS='OLD',ERR=444)
	OPEN(16,FILE=CLIMATE_PREP,STATUS='OLD',ERR=444)
	OPEN(46,FILE=SOIL_MC,STATUS='OLD',ERR=444)
C
C If there is an error stop programme
C	
      GOTO 333
444   CONTINUE
      WRITE(*,*)'ERROR.MSG'
      STOP
	CONTINUE
C
C Leave OPENCHANNEL_MC    
C
333	CONTINUE	
	END
C

C
C----------------------------------------------------------------------------------
C
	SUBROUTINE MC_input_prep (delayPlanting,WATERMODEL,SWAT_WATER)
C
C Subroutine to prepare input files - all variables are written to a file so no need to call any
C
	IMPLICIT NONE
C
C Variable declarations
C
	INTEGER ios				!internal counter
	INTEGER DOY				!day of the year
	INTEGER I				!internal counter
C	
C Weather
C
	REAL precip				!daily precipitation
	REAL evap				!daily evaporation
	REAL temp				!daily temperature
C
C Crops
C
	REAL samples(19)				!samples created for Monte Carlo by Simlab
	REAL carbonTotRInitYear			!Total carbon return in the initial year
	REAL carbonTotRGrowSeason		!Total carbon return in the growing season
	REAL carbonRetStubInitYear		!Carbon return as stubble in the initial year
	REAL carbonRetStubGrowSeason	!Carbon return as stubble in the growing season
	REAL carbonRetStrawInitYear		!Carbon return as straw in the initial year
	REAL carbonRetStrawGrowSeason	!Carbon return as straw in the grwoing season
	REAL nitrogenReqAboveInitYear	!N requirement aboveground in the initial year
	REAL nitrogenReqBelowInitYear	!N requirement belowground in the inital year
	REAL nitrogenReqAboveGrowS		!N requirement above ground in the growing season
	REAL nitrogenReqBelowGrowS		!N requirement below ground in the growing season
	REAL nitrogenRetStrawInitYear	!N return as straw in the initial year
	REAL nitrogenRetStrawGrowSeason	!N return as straw in the growing season
	REAL DPM_RPM_SOILN_ratio		!Ration of soilN DPM/RPM
C
C	SETup file
C
	CHARACTER(LEN=100) text(20)		!?
	INTEGER NSOIL					!Soil number
	INTEGER delayPlanting
	INTEGER cropNo					!crop number
	INTEGER plantingdate			!planting date
	INTEGER NRequire				!N requirement
	INTEGER harvestDate				!harvest date
	INTEGER strawFlag				!straw???
	INTEGER noFert					!No fertilizer
	INTEGER noManure				!No manure
	INTEGER dateFert				!date of fertilizer application
	INTEGER ureaFlag				!urea??
	INTEGER n15Flag					!?
	REAL amountFert					!Amount of fertilizer used
	REAL no3Fert					!Fertilizer as nitrate
	REAL nh4Fert					!Fertilizer as ammonium
	REAL ureaFert					!Fertilizer as urea
	REAL expYld						!Expected yield
	REAL amountManure				!Amount of manure
	INTEGER dateManure				!Date of manure application
	INTEGER typeManure				!Type of manure
	INTEGER lableManure				!Labelled manure
	INTEGER WATERMODEL
	INTEGER SWAT_WATER
C
C Open factor.tmp, which contains input parameters from MC_samples2 for the respective run
C
	OPEN (51, FILE="factor.tmp" , STATUS= 'OLD', ACTION = "READ", 
     &			IOSTAT=ios,ERR=1500)
	if (ios .NE. 0) then
		if (ios .LT. 0) then
			write(*,*) "End-of-file or end-of-record in factor.tmp"
		else
			write(*,*) "Error condition occured in factor.tmp"
		end if
	end if
C
C Open climate file to be manipulated for Monte Carlo
C
	OPEN (17, FILE="CLIMATEx.prn",STATUS='UNKNOWN',ERR=1600)
	OPEN (19, FILE="CROP_MC.txt" , STATUS= 'UNKNOWN',ERR=1800)
	OPEN (10, FILE="SETUP.SET" , STATUS='REPLACE',ERR=2000)
C
C READ sample file from factor.tmp
C
	REWIND(51)
	READ(51,*) samples
C
C Make an integer from sample(7) = soil number
C		and sample(17) = delay in planting date
C
	samples(7) = INT(samples(7))
	delayPlanting = INT(samples(17))  
C
C ---------------------------------------------------------------------------------
C Weather file preparation: read CLIMATE_PREP file, modify and write to CLIMATEX file
C
	DO i=1,365,1
	READ(16,*) DOY, precip, evap, temp
	precip=precip+samples(9)*precip
	temp=temp+samples(8)
	WRITE(17,"(I8, F14.1,F10.1, F13.2)") DOY, precip, evap, 
     &				temp
	END DO
	CLOSE(17)
C
C END: Weather file preparation
C -------------------------------------------------------------------------------------
C 
C Crop parameter file preparation from MC_samples2
C
	carbonTotRInitYear=samples(1)
	carbonTotRGrowSeason=samples(10)
	carbonRetStubInitYear=samples(1)*samples(2)
	carbonRetStubGrowSeason=samples(10)*samples(11)
!	carbonRetStrawInitYear=carbonTotRInitYear*samples(6)
	carbonRetStrawInitYear=samples(6)
! carbonRetStrawGrowSeason=carbonTotRGrowSeason*samples(15)
	carbonRetStrawGrowSeason=samples(15)
	nitrogenReqAboveInitYear=samples(3)
	nitrogenReqBelowInitYear=samples(4)
	nitrogenReqAboveGrowS=samples(12)
	nitrogenReqBelowGrowS=samples(13)
!	nitrogenRetStrawInitYear=(nitrogenReqAboveInitYear+nitrogenReqBelowInitYear)*samples(5)
	nitrogenRetStrawInitYear=samples(5)
! nitrogenRetStrawGrowSeason=(nitrogenReqAboveGrowS+nitrogenReqBelowGrowS)*samples(14)
	nitrogenRetStrawGrowSeason=samples(14)
	DPM_RPM_SOILN_ratio=samples(16)
C
C Write crop parameters
C
	WRITE(19,*)
	WRITE(19,*) carbonTotRInitYear
	WRITE(19,*) carbonTotRGrowSeason
	WRITE(19,*) carbonRetStubInitYear
	WRITE(19,*) carbonRetStubGrowSeason
	WRITE(19,*) carbonRetStrawInitYear
	WRITE(19,*) carbonRetStrawGrowSeason
	WRITE(19,*) nitrogenReqAboveInitYear
	WRITE(19,*) nitrogenReqBelowInitYear
	WRITE(19,*) nitrogenReqAboveGrowS
	WRITE(19,*) nitrogenReqBelowGrowS
	WRITE(19,*) nitrogenRetStrawInitYear
	WRITE(19,*) nitrogenRetStrawGrowSeason
	WRITE(19,*) DPM_RPM_SOILN_ratio
C
C END: Crop parameter file preparation
C ---------------------------------------------------------------------------------
C
C SETUP file preparation
C
C ... first line
C
	READ(21,"(I1,A50)") NSOIL, text(1)
	NSOIL = samples(7)
	WRITE(10,"(I2,A50)") NSOIL, text(1)
C
C ... second line
C
	READ(21,"(A50)") text(1)
	write(10,"(A50)") text(1)
C
C ... third line if SWAT water option is on
C	
	IF(WATERMODEL.EQ.SWAT_WATER)THEN
	READ(21,"(A50)") text(1)
	write(10,"(A50)") text(1)
	ENDIF

C
C ... climate files
C
	READ(21,"(A50)") text(1)
	write(10,*) "' CLIMATEx.prn ' 0"
	IF(WATERMODEL.EQ.SWAT_WATER)THEN
	  READ(21,*)text(1)
	  WRITE(10,*)text(1)
	ENDIF
C
	READ(21,"(A50)") text(1)
	write(10,*) "' CLIMATEx.prn ' 0"
	 IF(WATERMODEL.EQ.SWAT_WATER)THEN
	  READ(21,*)text(1)
	  WRITE(10,*)text(1)
	 ENDIF
C
	READ(21,"(A50)") text(1)
	write(10,"(A50)") text(1)
C
C Modify planting date to mimic delay in planting
C
	READ(21,*) cropNo, plantingdate, NRequire, harvestDate, expYld, 
     &				strawFlag, noFert, noManure
	plantingdate=plantingdate+samples(17)
	WRITE(10,"(I3,I5,I2,I5,F4.1,I2,2I3)") cropNo, plantingdate, 
     &	NRequire, harvestDate, expYld, strawFlag, noFert, noManure
C
C Modify date of fertilizer application to match delay in planting
C Modify amount of fertilizer
C
	DO i=1,noFert,1
	read(21,*) amountFert, dateFert, no3Fert, nh4Fert, ureaFert,
     &				ureaFlag, n15Flag
	dateFert=dateFert+samples(17)
	amountFert=amountFert+amountFert*samples(18)
	write(10,"(F5.1,I5,3F6.1,2I2)") amountFert, dateFert, no3Fert, 
     &						nh4Fert, ureaFert,ureaFlag, n15Flag
	END DO
C
C Modify date of manure application to match delay in planting
C Modify amount of manure
C
	DO i=1,noManure,1
	read(21,*) amountManure, dateManure, typeManure, lableManure
	dateManure=dateManure+samples(17)
	amountManure=amountManure+amountManure*samples(19)
	write(10,"(F5.1,I5,I3,I2)") amountManure, dateManure, typeManure, 
     &								lableManure
	END DO
C
	CLOSE(10)
	GOTO 2200
C
C END: SETUP file preparation
C--------------------------------------------------------------------------------
C
C 'Program terminated normally'

!***  ERROR-Branches
 1500 CONTINUE

	STOP "ERROR - something is wrong with file factor.tmp !!"

 1600 CONTINUE

	STOP "ERROR - something is wrong with file CLIMATE_PREP !!"

 1800 CONTINUE

	STOP "ERROR - something is wrong with file CROP_MC.txt !!"

 2000 CONTINUE

	STOP "ERROR - something is wrong with file SETUP_PREP !!"
C 
C Leave Monte Carlo preparation
C
2200	END 
C
C------------------------------------------------------------------------------------
C 
      SUBROUTINE SETFILE_MC(ATM,IDATEFC,NSOILJ,IDRAINJ,IROCKJ,LCROP,
     &                   PREYLD,SECONDS,MODTYPE,N_STEPS,
     &                   IAWCJ,NYEARS,
     &                   ISTHARV,ISTYR,ISTART_TS,ILAST_TS,FIXEND,LAT,
     &                   FIXFILE,NGROW,NXYEARS,ICROPJ,ISOWNJ,OCROPNJ,
     &                   IHARVJ,EXYLDJ,NRESJ,NFERTJ,NORGMJ,FERTJ,IFERTJ,
     &                   JFERTJ,IVOLJ,ILABJ,ORGMJ,IORGMJ,JORGMJ,
     &                   IOLABJ,WTABLE,
C Variable required for SWAT water
     &					FIXFILEavp,WATERMODEL,SWAT_WATER,z)


C
C Subroutine to get filenames and simulation inputs
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      INTEGER MAXCROP,MAXFERT,MAXORGM,MAXWEATH,MAXGROW,IOS
      PARAMETER (MAXCROP=36,MAXFERT=5,MAXORGM=52,MAXWEATH=300)
      PARAMETER (MAXGROW=300)
      CHARACTER TEMPFIL*20
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
	REAL WTABLE
	CHARACTER*100 FIXFILEaVP(MAXWEATH) !required for SWAT water
	INTEGER WATERMODEL
	INTEGER SWAT_WATER
	REAL z ! altitude required for SWAT
C
C Read in Simulation Data...
C ...Soil Data
C
C Default values
C
      ATM=40
      IDATEFC=1
C
      OPEN (10, FILE="SETUP.SET" , STATUS='OLD',ERR=111)
C
C Values read from file
C
      READ(10,*,ERR=111)NSOILJ,IDRAINJ,IROCKJ,LCROP,PREYLD,IGRASSJ,
     &                  ATM,IDATEFC,TIMESTEP,MODTYPE
C
C ... set timesteps
C	
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
	  READ(10,*,ERR=111) z
	ENDIF

      IF(FIXEND.EQ.0)FIXEND=ILAST_TS
      DO 100 I=1,NYEARS,1
        READ(10,*,ERR=111)FIXFILE(I)
	 	IF(WATERMODEL.EQ.SWAT_WATER)THEN
	  READ(10,*,ERR=111)FIXFILEaVP(I)
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
C Write planting and harvest date to file
C
	OPEN(1122,FILE='DATES.txt',STATUS='REPLACE')
	WRITE(1122,*) ISOWNJ(1)-365,IHARVJ(1)-365
	CLOSE(1122)
C
C Close the Parameter Simulation File
C
      
C
C Bypass quit if no error found
C
      GOTO 222
C
C Go straight here if error in input file
C
111   CONTINUE
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
C--------------------------------------------------------
C
	SUBROUTINE SETCROPPARAM_MC(carbonTotRInitYear,carbonTotRGrowSeason,
     &				carbonRetStubInitYear,carbonRetStubGrowSeason,
     &				carbonRetStrawInitYear,carbonRetStrawGrowSeason,
     &				nitrogenReqAboveInitYear,nitrogenReqBelowInitYear,
     &				nitrogenReqAboveGrowS,nitrogenReqBelowGrowS,
     &				nitrogenRetStrawInitYear,
     &				nitrogenRetStrawGrowSeason,DPM_RPM_SOILN_ratio)

C Subroutine to read in general crop parameters from File GeneralCropParameters.txt

	REAL carbonTotRInitYear,carbonTotRGrowSeason
      REAL carbonRetStubInitYear,carbonRetStubGrowSeason
      REAL carbonRetStrawInitYear,carbonRetStrawGrowSeason
      REAL nitrogenReqAboveInitYear,nitrogenReqBelowInitYear
      REAL nitrogenReqAboveGrowS,nitrogenReqBelowGrowS
	REAL nitrogenRetStrawInitYear
      REAL nitrogenRetStrawGrowSeason
	REAL DPM_RPM_SOILN_ratio

	OPEN(19, FILE='CROP_MC.txt', ACTION='READ',STATUS='OLD')

	REWIND(19)
	READ(19,*)
	READ(19,*) carbonTotRInitYear
	READ(19,*) carbonTotRGrowSeason
	READ(19,*) carbonRetStubInitYear
	READ(19,*) carbonRetStubGrowSeason
	READ(19,*) carbonRetStrawInitYear
	READ(19,*) carbonRetStrawGrowSeason
	READ(19,*) nitrogenReqAboveInitYear
	READ(19,*) nitrogenReqBelowInitYear
	READ(19,*) nitrogenReqAboveGrowS
	READ(19,*) nitrogenReqBelowGrowS
	READ(19,*) nitrogenRetStrawInitYear
	READ(19,*) nitrogenRetStrawGrowSeason
	READ(19,*) DPM_RPM_SOILN_ratio

	REWIND(19)

C
C Leave SETCROPPARAM
C
      RETURN
      END
C	
C-------------------------------------------------------------------------------
C	
	SUBROUTINE MOREFIX_MC(IYEAR,MCROP,LCROP,ICROP,INSTRAW,
     &                   YLD,PREYLD,EXYLD,SORGC,SORGN, 
     &                   SXORGC,SXORGN,HZ1,TORGC,TORGN,TXORGC, 
     &                   TXORGN,CAO,CSC,STRAW,
     &                   ORGC,ORGN,XORGC,XORGN,
     &				carbonTotRInitYear,carbonTotRGrowSeason,
     &				carbonRetStubInitYear,carbonRetStubGrowSeason,
     &				carbonRetStrawInitYear,carbonRetStrawGrowSeason,
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
      INTEGER MAXCROP,MAXSOIL,MAXLAYER
	PARAMETER (MAXCROP=36)
      PARAMETER (MAXSOIL=50)
       PARAMETER (MAXLAYER=60)
C
C Variables passed from calling subroutine
C
      INTEGER IYEAR, MCROP, LCROP, ICROP, INSTRAW
      REAL YLD, PREYLD, EXYLD, SORGC, SORGN, SXORGC(MAXLAYER), SXORGN
	REAL HZ1(MAXLAYER)
      Real TORGC(MAXLAYER)
      Real TORGN(MAXLAYER)
      Real TXORGC(MAXLAYER)
      Real TXORGN(MAXLAYER)
	REAL CAO(5,0:MAXCROP), CSC(4,0:MAXCROP), STRAW(4,0:MAXCROP)
	REAL ORGC, ORGN, XORGC, XORGN
C 
C Variables specific for Monte carlo
C
	REAL carbonTotRInitYear,carbonTotRGrowSeason
      REAL carbonRetStubInitYear,carbonRetStubGrowSeason
      REAL carbonRetStrawInitYear,carbonRetStrawGrowSeason
	REAL STRAWC
C
C Variables required for SWAT water routines
C
	REAL PLANTBIOMASS		! Plant biomass at harvest [ kg ha-1]
	REAL PLANTHEIGHT		! Plant height at harvest [cm]
	REAL LAI				! Plant LAI at harvest
	REAL canMAX				! maximum amount of water that can be trapped in the canopy when the canopy is fully developed (mm H2O)
	REAL CONVER_F				! IN:Conversion between timestep & weeks
	REAL C1(0:MAXCROP)
	REAL T1(0:MAXCROP)
C 
C In first year calc Total and Stubble Carbon return from initial crop
C
      IF(IYEAR.LE.1)THEN
       MCROP=LCROP
       YLD=PREYLD
	 TXORGC=carbonTotRInitYear
	 SXORGC=carbonRetStubInitYear
	 STRAWC=carbonRetStrawInitYear
       CALL SETC_MC(TXORGC,SXORGC,STRAWC)
C
C Extra crop parameters required for SWAT water (not interferring with main programme)
C
	 CALL SETEXTRACROPPARAM(YLD,PLANTBIOMASS,PLANTHEIGHT,LAI,canMax)
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
	 STRAWC=carbonRetStrawGrowSeason
      END IF
C
C Calc total C returns for next crop
C
      MCROP=ICROP
      YLD=EXYLD
	TORGC=carbonTotRGrowSeason
	SORGC=carbonRetStubGrowSeason

      CALL SETC_MC(TXORGC,SXORGC,STRAWC)
C
C Calc N inputs for previous and current year
C
C Seems to be not used anymore! Pia, 20/09/06
      IF(IYEAR.EQ.0)THEN
       TXORGN=HZ1*TXORGC
      ELSE
       TXORGN=TORGN
      END IF
      TORGN=HZ1*TORGC
C
C Leave MOREFIX
C
      RETURN
      END
C
C--------------------------------------------------------------------------------------
C
	SUBROUTINE INIT_SUNDIAL_CROP_MC(IYEAR, LCROP, ICROP, 
     &                             INSTRAW, IEND, ISTHARV,   
     &                             ISOWN, MEND, NSOW, JSTOP, L_TSS,
     &                             YLD, PREYLD, EXYLD,  
     &                             SORGC,SORGN, SXORGC, SXORGN,  
     &                             HZ1, TORGC, TORGN,TXORGC, TXORGN, 
     &                             OCROPN, RORGN,DDAYS,
     &                             ORGC, ORGN, XORGC, XORGN,
     &                             NXYEARS,NDATE,FIXEND,LHARV,IHARV,
     &                             NORGM,IORGM,NFERT,IFERT,IANTHES,
     &                             SECONDS,CREQN,
     &					nitrogenReqAboveGrowS,nitrogenReqBelowGrowS,
     &				nitrogenRetStrawGrowSeason,DPM_RPM_SOILN_ratio,
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
      INTEGER MAXCROP,MAXORGM,MAXFERT,IK
	PARAMETER (MAXCROP=36,MAXORGM=52,MAXFERT=5)
	INTEGER MCROP,IROCKS(0:MAXCROP),ISTART,ISAVE
	REAL ROOTS,CREQN,TOTC,CINP,ANINP,TCNEED,BREQN
      REAL CAO(5,0:MAXCROP),UR(3,0:MAXCROP),UT(3,0:MAXCROP) 
	REAL INC(2,0:MAXCROP),CRATE(0:MAXCROP),ANRATE(0:MAXCROP)
      REAL CFACT(0:MAXCROP),SEED(0:MAXCROP),C1(0:MAXCROP),T1(0:MAXCROP)
      REAL RRG(0:MAXCROP),RRX(0:MAXCROP)
      REAL CSC(4,0:MAXCROP),STRAW(4,0:MAXCROP)

C ******************************************************************************
C Input for MC (response curve analysis), Pia, 12/09/06
C Declaration of general crop parameters
	REAL carbonTotRInitYear,carbonTotRGrowSeason
      REAL carbonRetStubInitYear,carbonRetStubGrowSeason
      REAL carbonRetStrawInitYear,carbonRetStrawGrowSeason
      REAL nitrogenReqAboveInitYear,nitrogenReqBelowInitYear
	REAL nitrogenReqAboveGrowS,nitrogenReqBelowGrowS
      REAL nitrogenRetStrawInitYear
      REAL nitrogenRetStrawGrowSeason,DPM_RPM_SOILN_ratio
C ******************************************************************************
C
C Variables passed to/from calling subroutine
C
      INTEGER IYEAR, LCROP, ICROP,INSTRAW,IEND,N_REMAIN
      INTEGER ISTHARV,N_STEPS,ISOWN,MEND,NSOW,JSTOP
	INTEGER L_TSS(0:MAXCROP)
	INTEGER IANTHES,TIMESTEP,NXYEARS,NDATE,FIXEND,LHARV,IHARV
	INTEGER NORGM,IORGM(MAXORGM),NFERT,IFERT(MAXFERT),MAXLAYER
      parameter (MAXLAYER=60)
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
C READ in parameters for general crop description
C 
	CALL SETCROPPARAM_MC(carbonTotRInitYear,carbonTotRGrowSeason,
     &				carbonRetStubInitYear,carbonRetStubGrowSeason,
     &				carbonRetStrawInitYear,carbonRetStrawGrowSeason,
     &				nitrogenReqAboveInitYear,nitrogenReqBelowInitYear,
     &				nitrogenReqAboveGrowS,nitrogenReqBelowGrowS,
     &				nitrogenRetStrawInitYear,
     &				nitrogenRetStrawGrowSeason,DPM_RPM_SOILN_ratio)
C
C Call subroutine MOREFIX to set parameters dependant on inputs
C from sites and to call SETC to set C inputs from previous crop
C
c      CALL MOREFIX_MC(IYEAR,MCROP,LCROP,ICROP,INSTRAW,
c     &             YLD,PREYLD,EXYLD,SORGC,SORGN, 
c     &             SXORGC,SXORGN,HZ1,TORGC,TORGN,TXORGC, 
c     &             TXORGN,CAO,CSC,STRAW,
c     &             ORGC,ORGN,XORGC,XORGN,
c     &				carbonTotRInitYear,carbonTotRGrowSeason,
c     &				carbonRetStubInitYear,carbonRetStubGrowSeason,
c     &				carbonRetStrawInitYear,carbonRetStrawGrowSeason,
C introduced for SWAT water routines
c     &				             PLANTBIOMASS,PLANTHEIGHT,LAI,C1,T1,
c     &							 canMax,CONVER_F)				
C
C In previous year to first growing season
C
      IF(IYEAR.EQ.0)THEN
C
C Set the crop N requirement
C
        CALL SETCROPN_MC(IYEAR,INSTRAW,MCROP,LCROP,ICROP,
     &                    SXORGN,SORGN,YLD,PREYLD,CREQN, 
     &                    OCROPN,ROOTS,EXYLD,TOTC,
     &                    CFACT,STRAW,INC,SEED,UR,UT,XORGN,ORGN,
     &					nitrogenReqBelowInitYear,
     &					nitrogenReqAboveInitYear,
     &					nitrogenRetStrawInitYear)
C 
C Set the crop C return without stubble & chaff = XORGC
C
        XORGC=TXORGC-SXORGC
C	write(*,*) "XORGC=TXORGC-SXORGC", XORGC,TXORGC,SXORGC
C
C RORGN is yearly input from roots and stubble N from harvest to
C harvest.
C In initial year, assume that previous stubble input, SXORGN, (i.e.
C from two crops before) is same as current stubble input.
C
        RORGN=XORGN+SXORGN
C	write(*,*) "RORGN=XORGN+SXORGN", RORGN,XORGN,SXORGN
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
C Leave INIT_CROP
C
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE SETCROPN_MC(IYEAR, INSTRAW, MCROP, LCROP, ICROP,
     &                    SXORGN, SORGN, YLD, PREYLD, CREQN, 
     &                    OCROPN, ROOTS, EXYLD, TOTC,
     &                    CFACT, STRAW, INC, SEED, UR, UT, XORGN, ORGN,
     &					nitrogenReqBelow,
     &					nitrogenReqAbove,nitrogenRetStraw)
C
C Subroutine to set the crop N requirement
C
      IMPLICIT NONE

      INTEGER MAXCROP
	PARAMETER (MAXCROP=36)
C
C INTEGER passed from subroutine
C 
      INTEGER IYEAR, INSTRAW, MCROP, LCROP, ICROP
C
C REAL passed from subroutine
C
      REAL SXORGN, SORGN, YLD, PREYLD, CREQN, OCROPN, ROOTS, EXYLD, TOTC
	REAL CFACT(0:MAXCROP), STRAW(4,0:MAXCROP), INC(2,0:MAXCROP) 
	REAL SEED(0:MAXCROP), UR(3,0:MAXCROP), UT(3,0:MAXCROP)
	REAL XORGN, ORGN
C
C Input for MC analysis
C
	REAL nitrogenReqBelow,nitrogenReqAbove,nitrogenRetStraw
C
C INTEGER local to subroutine
C 
      INTEGER IFIT
C
C REAL local to subroutine
C 
      REAL STRAWN, CLOSS
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
C  Nitrogen requirements below and above ground are given as inputs for Monte Carlo
C
     	  ROOTS=nitrogenReqBelow
	  CREQN=nitrogenReqAbove
c
C Allow fraction of crop N which may be lost by senescence after anthesis
C (e.g.5%)
C
      CLOSS=CFACT(MCROP)*CREQN
C
C Incorporation of stubble and chaff is an input for Monte Carlo
C In first year always incorporate chaff
C
       IF(IYEAR.EQ.0)THEN
        STRAWN=nitrogenRetStraw
        SXORGN=INC(1,MCROP)*CREQN+STRAWN
        XORGN=INC(2,MCROP)*ROOTS
        CREQN=CREQN+CLOSS-SEED(MCROP)+XORGN
         IF(CREQN.LE.0)CREQN=0.0001
C
C In subsequent years...
C
      ELSE
       STRAWN=nitrogenRetStraw
       SORGN=INC(1,MCROP)*CREQN+STRAWN
C
C Set ORGN = N requirement of the roots
C
       ORGN=INC(2,MCROP)*ROOTS
C 
C Set CREQN = N requirement of the tops
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
      SUBROUTINE SETC_MC(TC,SC,STRAWC)
C
C Subroutine to Calculate the total C returns for next crop
C
      IMPLICIT NONE
      integer MAXLAYER
       PARAMETER (MAXLAYER=60)
C REAL passed from subroutine
C 
      REAL TC(MAXLAYER) 
      REAL SC(MAXLAYER)
C
C REAL local to subroutine
C
      REAL RC(MAXLAYER), STRAWC
C
C  
C
        IF(STRAWC.GT. 0.)THEN
         RC=TC-SC
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
      SUBROUTINE TS_RES1_MC(IYEAR,IRYEAR,NXYEARS,SECONDS, 
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
	DATA MAXLAYER1 /30/
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
      DO 400 IL=1,MAXLAYER1
        TOTNIT=TOTNIT+SOILN(IL)
400   CONTINUE
C
C Ammonium...
C
      DO 500 IL=1,MAXLAYER1
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
      DO 550 IL=1,MAXLAYER1
        ROOTN=ROOTN+DPMNIT0(IL)+RPMNIT0(IL)
	  ROOTDPM=ROOTDPM+DPMNIT0(IL)
	  ROOTRPM=ROOTRPM+RPMNIT0(IL)
        BIOHUM=BIOHUM+BNIT0(IL)+HNIT0(IL)
550   CONTINUE
      END=ROOTN+BIOHUM+TOTNIT+TOTAMM
C
C Sum outputs during week
C
      SUBTOTOUT=DENIT+CACT+WLEACH+VOLAT-CLOSSX
C
C Sum N at end of week
C
      TOTOUT=END+SUBTOTOUT
C
      CALL GETCAC(CACTOT,CACT,CATOT15,CACT15)
      RETURN
      

	END
C-------------------------------------------------------------------------------------
	 SUBROUTINE MC_RESULTS(IYEAR,IK,NSOW,MEND,NSOIL,
     &					MAXLAYER,precip,EVAPW,AIRTEMP,
     &                   SXORGC,TORGC,SORGC,SEEDIN,
     &					SORGN,SXORGN,ORGN,CREQN,
     &                   DPMCARB0,RPMCARB0,BCARB0,HCARB0,
     &					DPMCARB,RPMCARB,BCARB,HCARB,
     &                   DPMNIT0,RPMNIT0,HNIT0,BNIT0,
     &                   SOILN,AMMN,CACT,
     &                    VOLAT,DENIT,ATM,TDNIT,THISFERT,
     &                   FYMFERT,SLEACH,delayPlanting,SOILW,
     &					WRATEDM,TRATEM)
C
	IMPLICIT NONE
C
C Define variables called from other subroutines
C
	INTEGER MAXLAYER
      REAL SLEACH(MAXLAYER)
	REAL SOILN(MAXLAYER)
	REAL AMMN(MAXLAYER)
	REAL HCARB0(MAXLAYER)
	REAL DPMCARB0(MAXLAYER)
	REAL RPMCARB0(MAXLAYER)
      REAL DPMNIT0(MAXLAYER)
	REAL RPMNIT0(MAXLAYER)
	REAL HNIT0(MAXLAYER)
	REAL BCARB0(MAXLAYER)
	REAL BNIT0(MAXLAYER)
	REAL DENIT
	REAL SEEDIN
	REAL CACT
	REAL TDNIT
	REAL CLARRAY(30)

	INTEGER IYEAR
	INTEGER IK
	INTEGER NSOW
	INTEGER MEND
	INTEGER IHARVJ
	INTEGER ISOWNJ
	REAL SXORGC
	REAL SXORGN
	REAL TORGC
	REAL SORGC
	REAL ORGN
	REAL CREQN
	REAL ATM
	REAL VOLAT
	REAL FYMFERT
	REAL THISFERT

	REAL CreturnInitYear
	REAL NreturnInitYear
	REAL CreturnHarvestInitYear
	REAL NreturnHarvestInitYear
	REAL NRequirement
	REAL TOT_CN
	REAL DPM_CN
	REAL RPM_CN
	REAL BIO_CN
	REAL HUM_CN
	REAL CreturnDuringGrowS
	REAL NreturnDuringGrowS
	REAL DPMCARB(MAXLAYER)
	REAL RPMCARB(MAXLAYER)
	REAL BCARB(MAXLAYER)
	REAL HCARB(MAXLAYER)
	REAL TOCARRAY(30)
	REAL MAXWAT(30)
	REAL TON(30)
	INTEGER delayPlanting
	REAL SORGN
	INTEGER harvestdate
	INTEGER sowingdate
	REAL FNULL
	REAL EVAPW
	REAL precip
	REAL AIRTEMP
C	
C Variables local to this subroutine
C	
C Climate
	REAL avgTemp
	REAL waterBalance
	REAL sumprecip
	REAL sumtemp
	REAL sumevap
C
C
C
	INTEGER ios
	REAL SUM_ATM
	REAL SUM_VOLAT
	REAL SUM_ACTUPT
	REAL SUM_NinCROP
	REAL SumDPMDecomp
	REAL SumRPMDecomp
	REAL SumBDecomp
	REAL SumHDecomp
	REAL DPMCARB_BEGINN
	REAL RPMCARB_BEGINN
	REAL HCARB_BEGINN
	REAL BCARB_BEGINN
	REAL TOTCARB_BEGINN
	REAL SOILN_BEGINN
	REAL AMMN_BEGINN
	REAL SUM_TDNIT
	REAL TOTNIT_BEGINN
	REAL DPMNIT_BEGINN
	REAL RPMNIT_BEGINN
	REAL HNIT_BEGINN
	REAL BNIT_BEGINN
	REAL HCARB_CHANGE
	REAL BCARB_CHANGE
	REAL DPMCARB_CHANGE
	REAL RPMCARB_CHANGE
	REAL TOTCARB_CHANGE
	REAL SUM_DENIT
	REAL SUM_LEACH
	REAL SUM_FERTN
	REAL SUM_FYMN
	REAL SUM_SOILW
	REAL avgSOILW
	REAL avgWRATEDM
	REAL avgTRATEM
	REAL SUM_WRATEDM
	REAL SUM_TRATEM

	INTEGER D1,D2,DE1,DE2,DELTA1,DELTA2
	INTEGER D1B,D2B,DE1B,DE2B,DELTA1B,DELTA2B
	REAL WATDEF



	INTEGER NSOIL
	INTEGER run
	REAL sampels
	INTEGER I
	INTEGER MC_DEPTH
	REAL WRATEDM
	REAL TRATEM
	REAL SOILW(MAXLAYER)
	MC_DEPTH=10
C								
C
C
C-----------------------------------------------------------------------------
C 
	OPEN (20, FILE="results" , STATUS= 'REPLACE', ACTION = "WRITE")
	OPEN (21, FILE="control_calc.txt", STATUS='OLD', ACTION = "WRITE", 
     &	POSITION = 'APPEND',IOSTAT=ios,ERR=2100)
	if (ios .NE. 0) then
		if (ios .LT. 0) then
			write(*,*) "End-of-file or end-of-record in factor.tmp"
		else
			write(*,*) "Error condition occured in factor.tmp"
		end if
	end if

C
C END: File connections and error messages
C----------------------------------------------------------------------------------
C
	IF (IYEAR == 0 .AND. IK == 1) THEN

	REWIND(46)
c	DO WHILE (.NOT.EOF(46))
c	READ(46,*) 
c	READ(46,*) I
c      READ(46,*)MAXWAT(I),FNULL,FNULL,FNULL,
c     &           FNULL,FNULL,FNULL,FNULL,FNULL,				
c     &           FNULL,FNULL,TOCARRAY(I),FNULL,			
c     &           FNULL,CLARRAY(I),FNULL,						
c     &           FNULL,FNULL,FNULL,FNULL,TON(I)
c	ENDDO
C
C VARIABLES AT HARVEST IN INITIAL YEAR:
C
C Carbon returned to the field  (stubble, chaff, straw)
C Nitrogen returned to the field  (stubble, chaff, straw)
C
	CreturnHarvestInitYear=SXORGC
	NreturnHarvestInitYear=SXORGN
	END IF

C 
C------------------------------------------------------------------------------
C VARIABLES AT START OF CURRENT GROWING SEASON (IN 50 CM PROFILE):
C
C Output for carbon pool changes
C save carbon pools at the beginning of growing season
C
	IF (IYEAR .GT. 0 .AND. IK == NSOW) THEN
		DPMCARB_BEGINN=SUM(DPMCARB0(1:MC_DEPTH))
		RPMCARB_BEGINN=SUM(RPMCARB0(1:MC_DEPTH))
		BCARB_BEGINN=SUM(BCARB0(1:MC_DEPTH))
		HCARB_BEGINN=SUM(HCARB0(1:MC_DEPTH))
		TOTCARB_BEGINN=DPMCARB_BEGINN+RPMCARB_BEGINN+BCARB_BEGINN
     &					+HCARB_BEGINN
C
C save soil mineral N status at the beginning of the growing season
C
		SOILN_BEGINN=SUM(SOILN(1:MC_DEPTH))
		AMMN_BEGINN=SUM(AMMN(1:MC_DEPTH))
		DPMNIT_BEGINN=SUM(DPMNIT0(1:MC_DEPTH))
		RPMNIT_BEGINN=SUM(RPMNIT0(1:MC_DEPTH))
		BNIT_BEGINN=SUM(BNIT0(1:MC_DEPTH))
		HNIT_BEGINN=SUM(HNIT0(1:MC_DEPTH))
		TOTNIT_BEGINN=DPMNIT_BEGINN+RPMNIT_BEGINN+BNIT_BEGINN
     &					+HNIT_BEGINN
C
C write to file
C
	TOT_CN=TOTCARB_BEGINN/TOTNIT_BEGINN
	DPM_CN=DPMCARB_BEGINN/DPMNIT_BEGINN
	RPM_CN=RPMCARB_BEGINN/RPMNIT_BEGINN
	BIO_CN=BCARB_BEGINN/BNIT_BEGINN
	HUM_CN=HCARB_BEGINN/HNIT_BEGINN

	END IF
C -----------------------------------------------------------------------------
C VARIABLES DURING GROWING SEASON WITH EXCEPTIONAL CALCULATIONS:
C
C Carbon returned to the field during growing season
C
	IF (IYEAR .GT. 0 .AND. IK==1) THEN
	CReturnDuringGrowS=TORGC-SORGC !ORGC
	END IF
C 
C Nitrogen returned to the field during growing season
C Crop requirement N in initial year
C
	IF (IYEAR .GT. 0 .AND. IK==NSOW) THEN
	NReturnDuringGrowS=ORGN
	NRequirement=CREQN
	END IF
C 
C --------------------------------------------------------------------------------
C VARIABLES DURING GROWING SEASON:
C
C Calculation of 
C	atmospheric N additions - SUM_ATM
C	volatilization - SUM_VOLAT
C	ammonium N from FYM - SUM_FYMNH4
C	N fertilizer additions - SUM_FERT
C	Net mineralization - TDNIT
C	Denitrification - SUM_DENIT
C	Leaching - SUM_LEACH
C	Actual plant uptake - SUM_ACTUPT
C	Soil water - SUM_SOILW AS PERCENT OF MAXWAT
C	water rate modifier - SUM_WRATEDM
C	temperature rate modifier - SUM_TRATEM
C
	IF (IYEAR .GT. 0 .AND. IK .GE. NSOW) THEN
	 SUM_ATM=SUM_ATM+ATM
	 SUM_VOLAT=SUM_VOLAT+VOLAT
	 SUM_FYMN=SUM_FYMN+FYMFERT
	 SUM_FERTN=SUM_FERTN+THISFERT
	 SUM_TDNIT=SUM_TDNIT+TDNIT
	 SUM_DENIT=SUM_DENIT+DENIT
	 SUM_LEACH=SUM_LEACH+SLEACH(MC_DEPTH)
	 SUM_NinCROP=SUM_NinCROP+CACT 
       SUM_ACTUPT=SUM_ACTUPT+CACT-SEEDIN
	SUM_SOILW=SUM_SOILW+(SUM(SOILW(1:MC_DEPTH)))/(2*MAXWAT(NSOIL))*100
	 SUM_WRATEDM=SUM_WRATEDM+WRATEDM/10	! AVERAGE OF 10 LAYERS
	 SUM_TRATEM=SUM_TRATEM+TRATEM/10     !AVERAGE OF 10 LAYERS

	END IF
C
C calculation of absolute losses
C
	IF (IYEAR .GT. 0 .AND. IK .GE. NSOW) THEN
	 SumDPMDecomp=SumDPMDecomp+(SUM(DPMCARB0(1:MC_DEPTH))-
     &			SUM(DPMCARB(1:MC_DEPTH)))
	 SumRPMDecomp=SumRPMDecomp+(SUM(RPMCARB0(1:MC_DEPTH))-
     &			SUM(RPMCARB(1:MC_DEPTH)))
	 SumBDecomp = SumBDecomp + (SUM(BCARB0(1:MC_DEPTH))-
     &			SUM(BCARB(1:MC_DEPTH)))
	 SumHDecomp = SumHDecomp + (SUM(HCARB0(1:MC_DEPTH))-
     &			SUM(HCARB(1:MC_DEPTH)))
	END IF
C
	IF (IYEAR .GT. 0 .AND. IK==MEND) THEN
C
C Calcualte relative change in carbon pools from beginning to en of growing season
C	
		DPMCARB_CHANGE=SumDPMDecomp/DPMCARB_BEGINN
		RPMCARB_CHANGE=SumRPMDecomp/RPMCARB_BEGINN
		BCARB_CHANGE=SumBDecomp/BCARB_BEGINN
		HCARB_CHANGE=SumHDecomp/HCARB_BEGINN		
	END IF
	
C	 Climate averages over growing season
C Calculate distance between fertilizer application and rainfall
C
	IF (IYEAR.GT.0.AND.IK==NSOW)THEN
		D1=1
	ENDIF

	IF (IYEAR.GT.0.AND.IK.GT.NSOW)THEN
		DE1=DE1+D1
		IF (precip.GT.0)THEN
			DELTA1=DE1
			D1=0
		ENDIF
	ENDIF
	
	IF (IYEAR.GT.0.AND.IK==NSOW+42)THEN
		D2=1
	ENDIF

	IF (IYEAR.GT.0.AND.IK.GT.NSOW+42)THEN
		DE2=DE2+D2
		IF (precip.GT.0)THEN
			DELTA2=DE2
			D2=0
		ENDIF
	ENDIF
C
C Calculate distance between fertilizer application and rainfall
C
	IF (IYEAR.GT.0.AND.IK==NSOW)THEN
		D1B=1
	ENDIF

	IF (IYEAR.GT.0.AND.IK.GT.NSOW)THEN
		DE1B=DE1B+D1B
		IF (precip.GT.0)THEN
			DELTA1B=DE1B
			D1B=0
		ENDIF
	ENDIF
	
	IF (IYEAR.GT.0.AND.IK==NSOW+42)THEN
		D2B=1
	ENDIF

	IF (IYEAR.GT.0.AND.IK.GT.NSOW+42)THEN
		DE2B=DE2B+D2B
	WATDEF=MAXWAT(NSOIL)*2-(SUM(SOILW(1:MC_DEPTH)))
		IF (precip.GT.WATDEF)THEN
			DELTA2B=DE2B
			D2B=0
		ENDIF
	ENDIF
C

C
C Sum up weather data
C
	IF (IYEAR .GT. 0 .AND. IK .GE. NSOW) THEN
	  sumtemp=sumtemp+AIRTEMP
	  sumprecip=sumprecip+precip
	  sumevap=sumevap+EVAPW
	ENDIF

	IF (IYEAR .GT. 0 .AND. IK==MEND) THEN
	  avgTemp=sumtemp/(MEND-NSOW)
		avgSOILW=SUM_SOILW/(MEND-NSOW)
	avgWRATEDM=SUM_WRATEDM/(MEND-NSOW)
	avgTRATEM=SUM_TRATEM/(MEND-NSOW)
	  waterBalance=sumprecip-sumevap
	ENDIF
C
C
C Calculate average temperature and water balance during growing season
C avgTemp=SUM(temp(sowingDate:harvestDate))/(harvestDate-sowingDate)
C waterBalance=SUM(precip(sowingDate:harvestDate))-SUM(evap(sowingDate:harvestDate))
C	
C
C



	CLOSE(25)
	OPEN (25, FILE="Run.txt", STATUS='OLD', ACTION='READ')
	read(25,*) run



	IF (IYEAR.EQ.1.AND.IK.EQ.MEND)THEN
		write (*,*) "run", run
	WRITE(20,3010) run, NSOIL,MEND-NSOW,      
     &	DELTA1,DELTA2,DELTA1B,DELTA2B, 
C 2F10.1
     &	claRRAY(NSOIL),MAXWAT(NSOIL)*100/250,

C 23F10.3
     & CreturnHarvestInitYear,NreturnHarvestInitYear, 
     & CReturnDuringGrowS,NReturnDuringGrowS,NRequirement,TOT_CN,DPM_CN, 
     & RPM_CN, BIO_CN,HUM_CN,SUM_ATM,SUM_TDNIT,DPMCARB_Change,
     & RPMCARB_Change,BCARB_Change,HCARB_Change, 
     & SOILN_BEGINN,AMMN_BEGINN,SUM_VOLAT,SUM_DENIT, 
     & SUM_NinCROP,SUM_ACTUPT, SUM_LEACH,
C 9F10.1
     & SUM_FYMN,SUM_FERTN,TOCARRAY(NSOIL), 
     & TON(NSOIL),avgTemp,waterBalance,avgSOILW,
     & sumevap,sumprecip,
C 10F10.3
     & TOTCARB_BEGINN,DPMCARB_BEGINN,
     & RPMCARB_BEGINN,BCARB_BEGINN,HCARB_BEGINN,
     & TOTNIT_BEGINN,DPMNIT_BEGINN,RPMNIT_BEGINN,BNIT_BEGINN,
     & HNIT_BEGINN,avgWRATEDM,avgTRATEM,
C 1I5
     & delayPlanting
C
3010	FORMAT(7I10,2F10.1,23F10.3,9F10.1,12F10.3,1I5)
c
c	WRITE(20,*)
c	WRITE(20,*) run, "run"
c	WRITE(20,*) NSOIL, "soil number"
c	WRITE(20,*) claRRAY(NSOIL), "clay"
c	WRITE(20,*) CreturnHarvestInitYear, "CretunrHarvestinityear"
c	WRITE(20,*) NreturnHarvestInitYear, "Nreturnharvestinityear"
c	WRITE(20,*) CReturnDuringGrowS, "CReturnDuringGrowS" 
c	WRITE(20,*) NReturnDuringGrowS, "NReturnDuringGrowS"
c	WRITE(20,*) NRequirement,"NRequirement"
c	WRITE(20,*) TOT_CN, "TOT_CN"
c	WRITE(20,*) DPM_CN, "DPM_CN"
c	WRITE(20,*) RPM_CN, "RPM_CN"
c	WRITE(20,*) BIO_CN, "BIO_CN"
c	WRITE(20,*) HUM_CN, "HUM_CN"
c	WRITE(20,*) SUM_ATM, "SUM_ATM"
c	WRITE(20,*) SUM_TDNIT,"SUM_TDNIT"
c	WRITE(20,*) DPMCARB_Change, "DPMCARB_Change"
c	WRITE(20,*) RPMCARB_Change, "RPMCARB_Change"
c	WRITE(20,*) BCARB_Change, "BCARB_Change"
c	WRITE(20,*) HCARB_Change, "HCARB_Change"
c	WRITE(20,*) SOILN_BEGINN, "SOILN_BEGINN"
c	WRITE(20,*) AMMN_BEGINN, "AMMN_BEGINN"
c	WRITE(20,*) SUM_VOLAT, "SUM_VOLAT"
c	WRITE(20,*) SUM_DENIT, "SUM_DENIT"
c	WRITE(20,*) TXORGC, "TXORGC"
c	WRITE(20,*) SUM_ACTUPT, "SUM_ACTUPT"
c	WRITE(20,*) SUM_LEACH, "SUM_LEACH"
c	WRITE(20,*) SUM_FYMN, "SUM_FYMN"
c	WRITE(20,*) SUM_FERTN,"SUM_FERTN"
c	WRITE(20,*) TOCARRAY(NSOIL), "TOCARRAY(soilNo), "
c	WRITE(20,*) TON(NSOIL), "TON(soilNo),"
c	WRITE(20,*) avgTemp, "avgTemp"
c	WRITE(20,*) waterBalance, "waterBalance"
c	WRITE(20,*) TOTCARB_BEGINN, "TOTCARB_BEGINN"
c	WRITE(20,*) DPMCARB_BEGINN,"DPMCARB_BEGINN"
c	WRITE(20,*) RPMCARB_BEGINN, "RPMCARB_BEGINN"
c	WRITE(20,*) BCARB_BEGINN, "BCARB_BEGINN"
c	WRITE(20,*) HCARB_BEGINN, "HCARB_BEGINN"
c	WRITE(20,*) TOTNIT_BEGINN, "TOTNIT_BEGINN"
c	WRITE(20,*) DPMNIT_BEGINN, "DPMNIT_BEGINN"
c	WRITE(20,*) RPMNIT_BEGINN, "RPMNIT_BEGINN"
c	WRITE(20,*) BNIT_BEGINN, "BNIT_BEGINN"
c	WRITE(20,*) HNIT_BEGINN, "HNIT_BEGINN"
c	WRITE(20,*) delayPlanting, "delayPlanting"
	ENDIF

	GOTO 3012
	CONTINUE
	
1500  CONTINUE

	STOP "ERROR - something is wrong with file CLIMATEx.prn !!"

1600  CONTINUE

	STOP "ERROR - something is wrong with file MC_output.txt !!"

1900  CONTINUE

	STOP "ERROR - something is wrong with file SoilNumber.txt !!"

2100  CONTINUE

	STOP "ERROR - something is wrong with file control_calc.txt !!"

2200  CONTINUE

	STOP "ERROR - something is wrong with file DATES.txt !!"
C Leave MC_RESULTS
3012	CONTINUE

	END
C

C------------------------------------------------------------
C
      SUBROUTINE RUN1_SUNDIAL_CROP_MC(IYEAR,IK,MEND,IS_TS,JSTOP,
     &                             IEND,IANTHES,I_TS,
     &                             NSOW,ISTHARV,INSTRAW,ICROP,LCROP,
     &                             ISOWN,OCROPN,YLD,PREYLD,EXYLD,C_TS,
     &                             SXORGN,SORGN,TORGC,SORGC,RORGN,ORGC,
     &                             SXORGN15,SXORGC,XORGC,XORGN,ORGN,
     &                             TCINP,TRNIN,RNIN,RNIN15,
     &                             WR,CTOT,
     &                             CACTOT,CATOT15,ICOVER,
C required for Monte carlo
     &				nitrogenReqBelowGrowS,
     &					nitrogenReqAboveGrowS,
     &					nitrogenRetStrawGrowSeason,CREQN,
C required for SWAT water
     &				             PLANTBIOMASS,PLANTHEIGHT,LAI,
     &							 canMax)      
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
	INTEGER IYEAR,IK,NSOW,ISTHARV,IEND,MEND
	INTEGER IS_TS,INSTRAW,ICROP
	INTEGER LCROP,N_STEPS,ISOWN,JSTOP,N_REMAIN
      INTEGER IANTHES,I_TS
	REAL OCROPN,SXORGN,SORGN,CACTOT,CATOT15,TORGC,SORGC,RORGN,RNIN15
      REAL WR,YLD,PREYLD,CTOT,EXYLD,ORGC
      REAL RNIN,CONVER_F,TRNIN,SXORGN15,C_TS,SXORGC,TCINP
	REAL XORGN,ORGN,XORGC
	INTEGER ICOVER

C Crops FOR MONTE CARLO
C
	REAL nitrogenReqAboveGrowS		!N requirement above ground in the growing season
	REAL nitrogenReqBelowGrowS		!N requirement below ground in the growing season
	REAL nitrogenRetStrawInitYear	!N return as straw in the initial year
	REAL nitrogenretstrawgrowseason !N return as straw in the growing season
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
     &                  CONVER_F,0,0)
C
C
C In week of sowing set up parameters for the new sown crop
C
 
      IF(IYEAR.GT.0.AND.IK.EQ.NSOW)THEN
        CALL SOWING_MC(IYEAR,INSTRAW,MCROP,LCROP,ICROP,IEND,ISTHARV,
     &              N_STEPS, ISOWN, MEND, NSOW, IS_TS, JSTOP,
     &              OCROPN, SXORGN, SORGN, YLD, CREQN, PREYLD,
     &              TORGC, CTOT, EXYLD, TOTC, SORGC, RORGN, ORGC,
     &              CACTOT, CATOT15, TCNEED, BREQN, UR,
     &              UT, SEED, INC, CFACT, XORGN, ORGN, ROOTS,STRAW,
C required for Monte Carlo routines
     &				nitrogenReqBelowGrowS,
     &					nitrogenReqAboveGrowS,
     &					nitrogenRetStrawGrowSeason)
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
C Leave RUN1_SUNDIAL_CROP_MC
C
	END
C
C -----------------------------------------------------------------------
C
      SUBROUTINE SOWING_MC(IYEAR,INSTRAW,MCROP,LCROP,ICROP,IEND,ISTHARV,
     &                  N_STEPS,ISOWN,MEND,NSOW,IS_TS,JSTOP,
     &                  OCROPN,SXORGN,SORGN,YLD,CREQN,PREYLD,
     &                  TORGC,CTOT,EXYLD,TOTC,SORGC,RORGN,ORGC, 
     &                  CACTOT,CATOT15,TCNEED,BREQN,UR,
     &                  UT,SEED,INC,CFACT,XORGN,ORGN,ROOTS,STRAW,
C Variables required for Monte carlo
     &					nitrogenReqBelowGrowS,
     &					nitrogenReqAboveGrowS,
     &					nitrogenRetStrawGrowSeason)


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
C Variables specific for MONTE_CARLO
C
	REAL nitrogenReqBelowGrowS
      REAL nitrogenReqAboveGrowS
      REAL nitrogenRetStrawGrowSeason
C
C If Monte carlo is used set crop N according to external inputs
C
	CALL SETCROPN_MC(IYEAR, INSTRAW, MCROP, LCROP, ICROP,
     &                    SXORGN, SORGN, YLD, PREYLD, CREQN, 
     &                    OCROPN, ROOTS, EXYLD, TOTC,
     &                    CFACT, STRAW, INC, SEED, UR, UT, XORGN, ORGN,
     &					nitrogenReqBelowGrowS,
     &					nitrogenReqAboveGrowS,
     &					nitrogenRetStrawGrowSeason)
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
C Leave SOWING_MC
C
      RETURN
      END
C
