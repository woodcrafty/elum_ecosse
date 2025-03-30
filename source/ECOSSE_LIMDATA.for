C------------------------------------------------------------
C
C Site routines for ECOSSE Carbon and Nitrogen Turnover Model
C using on limited data to run the model
C
C Jo Smith 
C Started 28/06/07
C---------------------------------------------------
C
C
C*************************************************************
C MAIN RUN ROUTINES
C*************************************************************
C
C EXTERNAL SUBROUTINES
C ECOSSE_LIM_RUN()
C-------------------------------------------------------------
C
      SUBROUTINE ECOSSE_LIM_RUN(ISWAIT)
C
C Subroutine to run ECOSSE for GIS
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
	INTEGER MAXLUSOIL   		! Max.no.of land use types for soils MLR
	PARAMETER (MAXLUSOIL=6)
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=10)         ! MLR
	INTEGER MAXLU1, N			! Max.no.of land use types, N counter in DO loop
	DATA MAXLU1 /6/
	INTEGER MAXORGM				! Max.no.of organic manure applications allowed
	PARAMETER (MAXORGM=52)
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXSOIL				! Max.no.of soil types
	PARAMETER (MAXSOIL=50)
	INTEGER MAXSOMLAY			! Max.no.of soil layers entered			
	PARAMETER (MAXSOMLAY=10)

	INTEGER CH4_OFF				! CH4 model off
	INTEGER CH4_RICHARDS		! Richards CH4 model on		
	INTEGER CH4_AITKENHEAD      ! Aitkenhead CH4 model on
	INTEGER CNFOEREID			! C:N ratio obtained by method of Foereid
	INTEGER CNMAGEC				! C:N ratio obtained by method of MAGEC
      logical cultivate           ! Flags whether a cultivation event should take place
	INTEGER DBRADBURY			! Bradbury model for denitrification
	INTEGER DNEMIS				! Nemis model for denitrification
	INTEGER DOC_OFF				! DOC model off
	INTEGER DOC_ON				! DOC model on
	INTEGER EC_EQRUN_OFF		! Initialisation using a full ECOSSE equilibrium run is off
	INTEGER EC_EQRUN_ON			! Initialisation using a full ECOSSE equilibrium run is on
	INTEGER EQNPP				! Model initialised using plant inputs 
								!    from measured NPP
	INTEGER EQTOC				! Model initialised using measured TOC
	INTEGER EQNPPTOC			! Model initialised using both plant inputs
								!    from measured NPP and measured TOC
	INTEGER EQHILLIER, EQJONES	! Model initialised using TOC and
								!    the Hillier solver to get C pools and plant inputs
	! MLR inserted frac as quick fix to convert yield data into plant inputs
      real:: frac(10)	! soil input as fraction of npp for each land cover type
      real:: rootshoot(10)	! soil input as fraction of npp for each land cover type
      INTEGER ICFIXED				! Initialisation of C is fixed
	INTEGER ICROTHCEQ			! Initialisation of C by RothC 
								!  equilibrium run
	INTEGER IL					! Local counter variable
	INTEGER IMFUNC_ROTHC		! ROTHC moisture rate modifier
	INTEGER IMFUNC_HADLEY		! HADLEY moisture rate modifier
	INTEGER IMODDEC_ON			! Decomposition modification is on
	INTEGER IMODDEC_OFF			! Decomposition modification is off
	INTEGER IMODPI_ON			! Plant input modification is on
	INTEGER IMODPI_OFF			! Plant input modification is off
	INTEGER INPASSCN			! Initialisation of N by passing the
								!  C:N ratio of DPM and RPM
	INTEGER INPUTFACTORS		! Set factors for decomposition process. If = 1, asks user for value
	INTEGER IS_SUNDIAL			! Crop model type: 0=SUNDIAL
	INTEGER IS_MAGEC			! Crop model type: 1=MAGEC
	INTEGER ISDYNPI_OFF			! Adjustment of PI by external factors is off
	INTEGER ISDYNPI_ON			! Adjustment of PI by external factors is on
	INTEGER ISPINUP_OFF			! N limitation spin-up is not used
	INTEGER ISPINUP_ON			! N limitation spin-up is used
	INTEGER ITFUNC_ROTHC		! ROTHC temperature rate modifier
	INTEGER ITFUNC_HADLEY		! HADLEY temperature rate modifier
      integer :: lastlu           ! Land use code at the previsous timestep
      integer :: lu_change_yr     ! The last year when the land use changed
	INTEGER PH_PARAM			! pH set to neutral or passed from parameter file
      INTEGER PH_STATIC			! pH is read in from input
								!   file and does not change
	INTEGER PH_FROM_VSD			! pH is calculated using VSD (a very 
								!   simple dynamic version of the MAGIC model
								!   by Ed Rowe & Chris Evans, CEH, Bangor)
	INTEGER SPARFILE			! Soil parameters read in from file
	INTEGER SPARCALC			! Soil parameters calculated from TOC etc
	INTEGER SUNDIAL_WATER		!   SUNDIAL water model
	INTEGER SWAT_WATER			!	SWAT water model
	LOGICAL FULL_OUTPUT         ! Flags whether to write detailed output files
	LOGICAL SUMMARY_OUTPUT      ! Flags whether to write summary output file
C
C Data values
C
      DATA CH4_OFF,CH4_RICHARDS,CH4_AITKENHEAD /0,1,2/
	DATA CNMAGEC,CNFOEREID /1,2/
	DATA DBRADBURY,DNEMIS /1,2/
	DATA DOC_ON,DOC_OFF /1,2/
	DATA EC_EQRUN_OFF,EC_EQRUN_ON /0,1/
	DATA EQNPP,EQTOC,EQNPPTOC,EQHILLIER,EQJONES /1,2,3,4,5/
	DATA ICFIXED,ICROTHCEQ /1,2/
      DATA IMFUNC_ROTHC,IMFUNC_HADLEY /0,1/
	DATA IMODDEC_OFF,IMODDEC_ON /0,1/
	DATA IMODPI_OFF,IMODPI_ON /0,1/
	DATA INSTABLECN,INPASSCN /1,2/
	DATA INPUTFACTORS /0/
	DATA IS_SUNDIAL,IS_MAGEC /0,1/
	DATA ISDYNPI_OFF, ISDYNPI_ON /0,1/
	DATA ISPINUP_OFF,ISPINUP_ON /0,1/
      DATA ITFUNC_ROTHC,ITFUNC_HADLEY /0,1/
	DATA PH_PARAM,PH_STATIC,PH_FROM_VSD /0,1,2/			 
	DATA SPARFILE,SPARCALC /1,2/
	DATA SUNDIAL_WATER,SWAT_WATER /0,1/
C
C Variables passed to / from this subroutine
C
	REAL AIRTEMP				! IN: Air temperature (deg.C/timestep)	
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
      REAL ANCH4					! IN/OUT: Measured annual CH4 emissions
								!     kgC/ha/yr
	REAL ANCO2					! IN/OUT: Measured annual CO2 emissions
								!     kgC/ha/yr
	REAL ANDOC					! IN/OUT: Measured annual DOC loss
								!     kgC/ha/yr
      REAL ATM					! IN:Atmospheric N input (kgN/ha/timestep)
	REAL ATM15					! IN:Atmospheric N15 input (kgN15/ha/time)
      REAL AVEPET(12)				! IN:Long term average PET (mm/month)
      REAL AVERAIN(12)			! IN:Long term average rainfall(mm/month)
      REAL AVETEMP(12)			! IN:Long term average monthly average 
      REAL BALANCE(20)			! IN:C Results
	REAL BCARB(MAXLAYER)
	REAL BCARB0(MAXLAYER)		! IN:C in soil biomass (kgC/ha/layer)
      REAL BNIT0(MAXLAYER)		! IN:N in soil biomass (kgN/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN:N15 in soil humus (kgN15/ha/layer)
	REAL BSTART					! IN:Initial N in biomass pool (kgN/ha)
	REAL BSTART15				! IN:Initial N in biomass pool (kgN15/ha)
	REAL BULKDENS(MAXLAYER)		! IN:Bulk density of the layer (g/cm3)
	REAL C_TS					! IN:Litter C input in this timestep (kgC/ha)
      REAL CACCUM					! IN/OUT: Measured annual C accumulation 
								!     in the soil kgC/ha/yr
	REAL CACT					! IN:Crop N uptake (kgN/ha/timestep)
	REAL CACT15					! IN:Crop N15 uptake (kgN15/ha/timestep)
	REAL CACTOT					! IN:N taken up by crop (kgN/ha)
	REAL CAL(MAXLAYER)			! IN/OUT: Concentration of Al3+ in layer (units?)
	REAL CATOT15				! IN:N15 taken up by crop (kgN15/ha)
	REAL CDOCIN					! IN/OUT: Concentration of mobile DOC in this layer (mg/l)
	REAL CFACT(0:MAXCROP)		! IN:Crop parameter used to estimate yield
	REAL CH4(MAXLAYER)			! IN:CH4 emitted (kgC/ha/layer)
	INTEGER CH4MODEL			! OUT:Methane model (off, Richards or Aitkenhead model) 
	REAL CH4TOAIR				! OUT:CH4 released to atmosphere (Aitkenhead CH4 model) (kgC/ha)
	REAL CH4_PROD(MAXLAYER)  	! OUT: CH4 produced in each soil layer [kgC/ha/layer] (Richards CH4 model)
	REAL CH4_SOIL_OX			! OUT: CH4 oxidised from methane produced in the soil [kgC/ha] (Richards CH4 model)
      REAL CH4_ATMOS_OX			! OUT: CH4 oxidised from atmospheric CH4 [kgC/ha] (Richards CH4 model)
      REAL CH4_FLUX				! OUT: Net CH4 flux from soil column [kgC/ha] (Richards CH4 model)
      REAL CLAY(MAXLAYER)			! IN:Clay content of the layer (%)
	REAL CLOSSX					! IN:N lost by senescence (kgN/ha) 
	REAL CLOSSX15				! IN:N15 lost by senescence (kgN15/ha) 	
	REAL CNH4IN					! IN/OUT: Concentration of ammonium in this layer (eq/m3)
	INTEGER CNMODEL				! OUT:Type of Model for C:N ratio of DPM and 
	                            !     RPM (MAGEC or Foereid)
      REAL CNO3IN					! IN/OUT: Concentration of nitrate in this layer (eq/m3)
	REAL CO2(MAXLAYER)			! IN:CO2 emitted (kgC/ha/layer)
	REAL CO2FROMDOC(MAXLAYER)	! IN:CO2 output from DOC in each layer(kgC/ha)
	REAL CONC					! IN:Conc.N leached (kgN/ha/mm/layer)
	REAL CONC15					! IN:Conc.N15 leached (kgN/ha/mm/layer)
	REAL CRIT(MAXSOIL,MAXLAYER)	! IN:Critical N content (kgN/ha/layer)
	REAL CTOT					! IN:Total litter C input (kgC/ha)
	REAL CULTDEPTH				! OUT:Cultivation depth (cm)
	REAL CUPTN					! IN:Crop N uptake (kgN/ha/timestep)
	REAL DENITRIFN(MAXLAYER)	! IN/OUT: N lost by denitrification (kgN/ha/layer
	REAL DENIT					! IN:N lost by denitrification (kgN/ha)
	INTEGER DMODEL				! OUT:Type denitrification model chosen 
	                            ! (Bradbury or NEMIS)
	REAL DN15					! IN:N15 lost by denitrification (kgN15/ha)
	REAL DNIT(MAXLAYER)			! IN:Net Mineralised N (kgN/ha/lay/timestep)
	REAL DNIT15(MAXLAYER)		! IN:Net Mineralised N15 (kgN15/ha/lay/time)
	INTEGER DOCMODEL			! OUT:DOC model (on or off)
	REAL DOMSOILC(CELLSTOGRID*MAXSERIES/2,MAXLUSOIL,MAXSOMLAY)	! IN(SETFILE_LIM):Soil C in 5 
															! major soil series under different LU 
															! in SOM layers (kgC/ha)
	REAL DOMSOILCLAY(CELLSTOGRID*MAXSERIES/2,MAXLUSOIL,MAXSOMLAY)	! IN(SETFILE_LIM):% Soil clay in 5 
															! major soil series under different LU 
	REAL DOMSOILIMPDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLUSOIL)	! IN(SETFILE_LIM): Depth of impermeable layer if present (cm)
	INTEGER DOMSOILISIMP(CELLSTOGRID*MAXSERIES/2,MAXLUSOIL)	! IN(SETFILE_LIM): Flag for impermeable layer (0=No; 1=Yes)
	REAL DOMSOILSILT(CELLSTOGRID*MAXSERIES/2,MAXLUSOIL,MAXSOMLAY)	! IN(SETFILE_LIM):% Soil silt in 5 
															! major soil series under different LU 
	REAL DOMSOILSAND(CELLSTOGRID*MAXSERIES/2,MAXLUSOIL,MAXSOMLAY)	! IN(SETFILE_LIM):% Soil sand in 5 
															! major soil series under different LU 
	REAL DOMSOILBD(CELLSTOGRID*MAXSERIES/2,MAXLUSOIL,MAXSOMLAY)	! IN(SETFILE_LIM):Soil BD in 5 
															! major soil series under different LU 
															! in SOM layers (g/cm3)
	REAL DOMSOILPH(CELLSTOGRID*MAXSERIES/2,MAXLUSOIL,MAXSOMLAY)	! IN(SETFILE_LIM):Soil PH in 5 
												! major soil series under different LU 
												! in SOM layers (g/cm3)
      REAL DPM_RPM_RATIO			! IN/OUT:Ratio of DPM:RPM. If set to 
								!        zero, this will be worked out using a call to SET_DPMRPMRATIO
     	REAL DPMCARB(MAXLAYER)
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
      REAL DPMCTON(MAXLAYER)		! IN/OUT:DPM C:N ratio
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
      INTEGER DRAINCLASS			! IN(SETFILE_LIM):Drainage class
	REAL DRAINW(MAXLAYER)		! IN:Water drained from this layer (mm/layer)
	REAL DRRAT					! IN/OUT:DPM:RPM ratio
	INTEGER EC_EQRUN			! IN(SETFILE_LIM):Initialisation using a full ECOSSE equilibrium run (on or off)
	INTEGER EQMODEL				! OUT:Type of equilibrium run 
	                            ! (NPP, TOC or both)
	REAL EVAP					! IN: Potential evap. (mm/timestep)
	REAL EVAPW					! IN: Potential evap. (mm/timestep)
	REAL EXYLD					! IN:Yield of current crop (t/ha)
	REAL FERT(MAXFERT)			! IN:Amount of fertiliser applied (kgN/ha)
	REAL FERTNF					! INOUT:Amount of fertiliser applied (kgN/ha)
	REAL FIELDCAP(MAXLAYER)	    ! IN:Soil water content at field capacity (mm/layer)
	REAL FIXN					! IN:N fixed from atmosphere (kgN/ha)
	INTEGER FIXEND				! IN:Fixed end? 0=No 1=Yes
	REAL FLOWPROP				! IN:Proportion of flow needed to achieve
							    !	     water table at depth WTABLE
	REAL FYMFERT				! IN:N input by fertiliser & FYM (kgN/ha)
	REAL FYMFERT15				! IN:N15 input by fertiliser & FYM (kgN15/ha)
	REAL G15DN2					! IN:15N2 lost by denitrification (kgN15/ha)
	REAL G15DN2O				! IN:15N2O lost by denitrification (kgN15/ha)					
	REAL G15NN2O				! IN:15N2O lost by nitrification (kgN15/ha)
	REAL G15NNO					! IN:15NO lost by nitrification (kgN15/ha)			
	REAL G15PNN2O				! IN:15N2O lost by part.nitrif.(kgN15/ha)			
	REAL GDN2					! IN:N2 lost by denitrification (kgN/ha)
	REAL GDN2O					! IN:N2O lost by denitrification (kgN/ha)
	INTEGER GISOK				! IN:Check GIS data 0=WRONG, 1=OK
 	REAL GNN2O					! IN:N2O lost by nitrification (kgN/ha)
	REAL GNNO					! IN:NO lost by nitrification (kgN/ha)
	REAL GPNN2O					! IN:N2O lost by part.nitrification (kgN/ha)
      REAL GRIDLAT				! IN(SETFILE_LIM):Average latitude of this 20km2 grid cell
	REAL HCARB(MAXLAYER)
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN:N in soil humus (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN:N15 in soil biomass (kgN15/ha/layer)
	REAL HZ1(MAXLAYER)			! IN:N:C of input PM for steady state 
	INTEGER I_TS				! IN:No.timesteps since sowing
	INTEGER IANTHES				! IN:No.timesteps from 01/01 to anthesis 
	INTEGER IAWC				! IN(SETFILE_LIM):Water movement code number
 	REAL ICFACTOR(MAXLAYER)		! IN:Adjustment in rate needed to	
								!    achieve measured NPP and TOC	
	INTEGER ICMODEL				! OUT:Type of C initialisation model 
	INTEGER ICROP				! IN:Current crop code
      INTEGER ICOVER				! OUT:Crop cover 1=Covered 0=Bare
	INTEGER IDATEFC				! IN:Date of field capacity (1=01/01; 2=01/06)
      INTEGER IDEC				! OUT:Decade counter 
	INTEGER IEND				! IN:Timesteps from sowing date to harvest
	INTEGER IFERT(MAXFERT)		! IN:No.timesteps to split fert.application
	INTEGER IFILE				! IN:Current weather file
	INTEGER IHARV				! IN:Timesteps from 01/01/01 to harvest date 
	INTEGER IK					! OUT:No.timesteps from prev.harv. to current
	INTEGER IL_TSS				! IN:No.timesteps before harvest when senesces
	INTEGER ILAB(MAXFERT)		! IN:Labelling on fertiliser? 0=No 1=Yes
	INTEGER ILABNF				! INOUT:Labelling on fertiliser? 0=No 1=Yes
	INTEGER ILAST_TS			! IN:Last simulated timestep since 01/01/01
	INTEGER ILU					! IN:Counter for land use 
	INTEGER IMFUNC				! OUT:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
      INTEGER IMODDEC				! Modify decomposition according to 
	                            !   full run using long term average weather data
      INTEGER IMODPI				! Modify plant inputs according to 
	                            !   full run using long term average weather data
	REAL INMODE					! SPINUP mode number
	INTEGER INMODEL				! OUT:Type of N initialisation model
	INTEGER INSTABLECN			! Initisalisation of N by assuming 
								!  steady state (after Bradbury)
	INTEGER INSTRAW				! IN:Straw incorporated? 0=No 1=Yes
	INTEGER INVERT				! OUT(CULTIV):Invert / mix soil on cultivation? 0=No 1=Yes
	INTEGER IOLAB(MAXORGM)		! IN:Labelling on manure? 0=No 1=Yes
	INTEGER IOLABNF				! IN:Labelling on manure? 0=No 1=Yes
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
	INTEGER IORGM(MAXORGM)		! IN:No.timesteps to split manure application
	INTEGER IROCK				! IN:No.layers: 1=50cm, 2=100cm , 3=150cm
	INTEGER IRYEAR				! IN:Current weather year number
	INTEGER IS_TS				! IN:Crop counter
      INTEGER ISAVE				! OUT:Code to save or retrieve variables
								!		0=save; 1=retrieve
	INTEGER ISDYNPI				! Code for dynamic PI (adjusted by external factors)
	INTEGER ISERIES				! IN:Counter for soil series
      INTEGER ISERROR				! IN:Code for error in soil 1=error 0=noerror
	INTEGER ISOPEN(MAXGROW)		! Is the metfile open (1=yes; 0=no)
	INTEGER ISOWN				! IN(GETYRSET):Timesteps from 01/01/01 to sowing date 
	INTEGER ISPINUP				! Is N limitation spin-up used?
      INTEGER ISTART_TS			! IN:First simulated timestep since 01/01/01 
	INTEGER	ISTHARV				! IN:Timesteps from 01/01/01 to first harvest
	INTEGER ISTYR				! IN:First year in simulation (eg. 2001)
	INTEGER ISWAIT				! IN/OUT:Code to wait for key press (1) or not (0)
	INTEGER ITFUNC				! OUT:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IVOL(MAXFERT)		! IN:N volatilised from fertiliser (kgN/ha)
	INTEGER IVOLNF				! INOUT:N volatilised from fertiliser (kgN/ha)
	INTEGER IYEAR				! IN:Current growing season number
	INTEGER JORGM(MAXORGM)		! IN:Type of manure application
	INTEGER JORGMNF				! IN:Type of manure application
	INTEGER JSTOP				! IN:Code to stop simulation 0=No 1=Yes
	INTEGER K_TS				! IN:No.timesteps since start of year(UNUSED?)
	INTEGER L_TSS(0:MAXCROP)	! IN:No.timesteps before harvest
	REAL LAT					! IN: Latitude
	INTEGER LCROP				! IN:Previous crop type
	REAL LEACHDOC(MAXLAYER)		! IN/OUT:DOC leached from the layer (kgC/ha/layer/timestep)
	REAL LEACHDON(MAXLAYER)		! IN/OUT:DON leached from the layer (kgN/ha/layer/timestep)
	INTEGER LHARV				! IN:No.timesteps from 01/01 to prev. harvest
	REAL LTA_AWC(12,MAXLAYER)	! IN:Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)	! IN:Average air temp this month (deg.C)
      INTEGER LUSEQ(MAXGROW)		! IN(SETFILE_LIM):Land use sequence
	INTEGER LU1					! OUT:First land use
	INTEGER LU2					! OUT:Land use changed to
	INTEGER MCROP				! IN:Crop type code number
	INTEGER MEASLAY				! Layer that soil is measured to
 	INTEGER MEND				! IN:Number of timesteps from prev.harv.-harv.
	CHARACTER*40 METFILE(MAXGROW)	! IN(SETFILE_LIM):Met.file  
	Real Miscticker				! Ticker to track miscanthus growth
	REAL MOBDOC(MAXLAYER)		! IN:Mobile DOC in each layer (kgC/ha)
	REAL MOBDON(MAXLAYER)		! IN:Mobile DON in each layer (kgN/ha)
      INTEGER MODTYPE				! OUT:Crop model type: 0=SUNDIAL, 1=MAGEC
	INTEGER N_STEPS				! IN:No.timesteps in a year
	INTEGER N_TS				! IN:No.timesteps since start of year(UNUSED?)
	INTEGER NDATE				! IN:Not used?
	INTEGER NF					! IN:No.fertiliser application
	INTEGER NFERT				! IN:No.fertiliser applications to this crop
	INTEGER NGROW				! IN:No.growing seasons simulated
	REAL NITRIFN(MAXLAYER)		! IN/OUT:Nitrification from this layer (kg/ha/layer/timestep)
	INTEGER NORGM				! IN:No.manure applications to this crop
	INTEGER NSOIL				! IN(SETFILE_LIM):Soil code number
	INTEGER NSOMLAY(CELLSTOGRID*MAXSERIES/2,MAXLUSOIL)	! IN:Number of SOM layers
	INTEGER NSOW				! IN:Timesteps from prev.harv. to sowing date
	INTEGER NXYEARS				! IN(SETFILE_LIM):No.growing seasons simulated
	INTEGER NYEARS				! IN:No.weather years included in the sim.
	REAL ORGC					! IN:Total org.C input (kgC/ha) 
	REAL ORGMA(MAXORGM)			! IN:Amount of manure applied (t/ha)
	REAL ORGMANF				! IN:Amount of manure applied (t/ha)
      REAL ORGMJ(0:MAXGROW,MAXORGM)	! IN:Amount of manure applied (t/ha)
	REAL ORGN					! IN:Total org.N input (kgN/ha) 
	INTEGER PH_MODEL			! How is pH calculated?
	REAL PIANN					! OUT: Total annual plant C input 
	                            !      (kg C / ha / yr) 
								!      - used as temprary value to allow adjustment of PI according 
								!        to external factors
      REAL PI_C(MAXLAYER)			! IN/OUT: Plant input C to soil (kgC/ha/step)
	REAL PI_CEQ_MON(12,MAXLAYER)! IN/OUT: Equilibrium plant C input each month
	Real Equ_plant(12,MAXLAYER) !         in each layer (kgC/ha/month/layer)
      REAL PI_N(MAXLAYER)			! IN/OUT: Plant input N to soil (kgC/ha/step)
      REAL PI_N15(MAXLAYER)		! IN/OUT: Plant input N15 to soil(kgC/ha/step)
								!(1=arable,2=grassland,3=forestry,4=semi-natural)
	REAL PI_SPEC(MAXGROW)		! IN(SETFILE_LIM): Total annual plant C input (kg C / ha / yr) 
	REAL PLANTUP(MAXLAYER)		! IN(RUN2_GIS_CROP)/OUT(TEST1_RES):Plant uptake from each soil layer
	REAL PNREQ(MAXLU)			! IN/OUT:Amount of N required by plant each yr
	                            !        (kgN/ha/year)
	REAL PREYLD					! IN:Yield of previous crop (t/ha)
	REAL OCROPN					! IN:N uptake of crop (kgN/ha) 
	REAL RAIN					! IN: Rainfall (mm/timestep)
	REAL REFIL(MAXLAYER)		! IN:Water needed to fill to FC (mm/lay)
	REAL RNIN					! IN:Litter N input in this timestep (kgN/ha)
	REAL RNIN15					! IN:Litter N15 input in timestep (kgN15/ha)
	REAL ROOT					! IN:Rooting depth according to restriction (cm)
	REAL RORGN					! IN:Total org.N input (kgN/ha) DONT PASS?
	REAL RPMCARB(MAXLAYER)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
      REAL RPMCTON(MAXLAYER)		! IN/OUT:RPM C:N ratio
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL RSTART					! IN:Initial N in debris pool (kgN/ha)
	REAL RSTART15				! IN:Initial N in debris pool (kgN15/ha)
	INTEGER RUNMODE				! Run using normal run mode for a 
								!  specified number of years
								! or using C accumulation
								!  & C losses to estimate plant input
								!  and spinnig up from 0 C to current 
								!  measured state
	REAL SAND(MAXLAYER)			! IN:Sand content of the layer (%)
	REAL SATWATCONT(MAXLAYER)	! IN:Total water content at saturation (mm/layer)
	REAL SC						! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SECONDS				! IN:Number of seconds in one timestep
	REAL SEEDIN					! IN: Data used by MAGEC
	REAL SILT(MAXLAYER)			! IN:Silt content of the layer (%)
	REAL SLEA15(MAXLAYER)		! IN:Nitrate-N15 leached (kgN15/ha)
	REAL SLEACH(MAXLAYER)		! IN:Nitrate-N leached (kgN/ha)
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOILPH(MAXLAYER)		! IN/OUT:pH of soil in this layer
	REAL SOILTEMP(MAXLAYER)		! IN: Soil temperature (deg.C/timestep)
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL SOMDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLUSOIL,MAXSOMLAY)	! IN:depth of soil organic matter layers
	REAL SORGC					! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SORGN					! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SORGN15				! IN:Total org.N15 input (kgN/ha) DONT PASS?
	INTEGER SPARMODEL			! OUT:Soil parameter model (from file or calc)
	REAL SRNIN					! IN:Total litter N input (kgN/ha)
      logical(4) START				! OUT: Starting pH run? 0=No, 1=Yes
	INTEGER SUM_TS				! IN:Total number of timesteps passed
	REAL SX						! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SXORGC					! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SXORGN					! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SXORGN15				! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL T15					! IN:Net Mineralised N (kgN/ha/layer)
	REAL TACTOT					! IN:Total crop N uptake (kgN/ha)
	REAL TAM(MAXLAYER)			! IN:15N/N in the ammonium in the layer
	REAL TC						! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL TCINP					! IN:Total org.C input (kgC/ha) DONT PASS?
	REAL TDNIT					! IN:Net Mineralised N (kgN/ha/layer)
	REAL TFERT(MAXFERT,3)		! IN:Prop.NO3,NH4,urea in fertiliser
	REAL TFERTNF(3)				! INOUT:Prop.NO3,NH4,urea in fertiliser
	REAL TFYMC					! IN:Total FYM C input (kgC/ha) DONT PASS?
      REAL THICK					! OUT: Thickness of the layer (m)
	REAL THISFERT				! IN:N input by fertiliser (kgN/ha)
	REAL THISFERT15				! IN:N15 input by fertiliser (kgN15/ha)
	INTEGER THISLU				! OUT:LU in the current growing season
	REAL TIN(MAXLAYER)			! IN:15N/N in the nitrate in the layer
	REAL TOC(MAXLAYER)			! IN:Measured total organic C (kgC/ha/layer)
      REAL TORGC					! IN:Total org.C input (kgC/ha) DONT PASS?
	REAL TORGN					! IN:Total org.N input (kgN/ha) DONT PASS?
	REAL TOTAST					! IN:Total ammonium in soil (kgN/ha)
	REAL TOTNST					! IN:Total nitrate in soil (kgN/ha)
	REAL TOTAST15				! IN:Total ammonium in soil (kgN15/ha)
	REAL TOTNST15				! IN:Total nitrate in soil (kgN15/ha)
	REAL TOTPIC(MAXLU)			! IN/OUT:Plant C input calculated using 
								!	  ROTHC equilibrium run for each LU (kgC/ha/yr)
	REAL TRATEM
	REAL TRNIN					! IN:Total litter N input (kgN/ha)
	REAL TXORGC					! IN:Total org.C input (kgC/ha) DONT PASS?
	REAL TXORGN					! IN:Total org.N input (kgN/ha) DONT PASS?
	REAL VIGOUR					! OUT: Vigour of cultivation (0-1) Determines the proportion of humus released to biomass,DPM and RPM pools
	REAL VOLAT					! IN:N lost by volatilisation (kgN/ha)
	REAL VOLAT15				! IN:N15 lost by volatilisation (kgN15/ha)
	REAL VOLATN(MAXLAYER)		! IN/OUT:Volatilisation from this layer (kg/ha/layer/timestep)
	REAL VSD_STATE_IL(10)		! IN: VSD State layer IL
	REAL WATERIN				! OUT: Water moving into the layer (m/timestep)
	INTEGER WATERMODEL			! IN(READ_MODEL):Water model
	REAL WATEROUT				! OUT: Water moving out of the layer (m/timestep)
	REAL WILTPOINT(MAXLAYER)	! IN:Water content at wilting point (mm/layer)
	REAL WLEACH					! IN:N leached (kgN/ha)
	REAL WLEACH15				! IN:N15 leached (kgN15/ha)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
	REAL WR						! IN:Root requirement (kgN/ha)
	REAL WRATEDM
	REAL WTABLE					! IN:Water table depth in cm
	REAL XORGC					! IN:Total org.C input (kgC/ha) 
	REAL XORGN					! IN:Total org.N input (kgN/ha) 
	REAL YLD,MAXNITER			! IN:Yield (t/ha) DONT PASS?
	Real ANNRain,anntemp
	integer IMON, NITER, P		!Counters in DO loops
	REAL ACTTOC(MAXLAYER)		!Actual measured TOC
	REAL ACTTOCDASH(MAXLAYER)	!Est. TOC
	REAL CROP(12)            ! Crop cover during each month, 1=Covered, 0=Bare
      REAL DIFF			!Conversion of est PI in all LU types with spinup 6.1
	REAL DDAYS

      ! Land use types: 1 ara, 2 gra, 3 for, 4 nat, 5 mis, 6 src, 7 sug, 8 osr, 9 srf, 10 whe)
      frac = [0.53, 0.71, 0.88, 0.71, 0.36, 0.35, 0.25, 0.75, 0.4, 0.53]  ! MLR quick fix for elum
      rootshoot = [0.0, 0.0, 0.0, 0.0,
     &             0.35, 0.25, 0.15, 0.15, 0.25, 0.15]
      rootshoot = 1 + (rootshoot / (1 - rootshoot))
C
C Set flag for opening metfiles
C
      DO 25 IL=1,MAXGROW
	  ISOPEN(IL)=0
25    CONTINUE
C
C Set restrictions on calculation
C Get parameters from TOC,IOM,CLAY,SILT,SAND,BULKDENS,	
C no restriction due to drainage,
C initialisation using ROTHC model & Bente Foereid's C:N ratios
C and get pH / water restriction on decomposition from NPP and TOC
C Use only 1 decade calculation
C
      DMODEL=2
	MODTYPE=IS_SUNDIAL
	SPARMODEL=SPARCALC
	WTABLE=300	
      ISPINUP=ISPINUP_OFF
	IMODPI=IMODPI_OFF
	IMODDEC=IMODDEC_OFF
      IMFUNC=IMFUNC_ROTHC
      ITFUNC=ITFUNC_ROTHC
	EC_EQRUN=EC_EQRUN_OFF
	INMODEL=INPASSCN
C      PH_MODEL=PH_FROM_VSD
      PH_MODEL=PH_STATIC
      CH4MODEL=CH4_RICHARDS !RICHARDS
 	DOCMODEL=DOC_OFF
      FULL_OUTPUT=.FALSE.
      SUMMARY_OUTPUT=.TRUE.
      CALL TELL_MODEL(DMODEL,ICMODEL,INMODEL,CNMODEL,DOCMODEL,SPARMODEL,
     &				EQMODEL,IMFUNC,ITFUNC,CH4MODEL,ISPINUP,ISWAIT)
C
C Open Channels for file input/output
C
	CALL TEST1_OPENCHAN(WATERMODEL,FULL_OUTPUT,SUMMARY_OUTPUT)
C
C Loop back to next cell
C
101   CONTINUE
C
C Get simulation inputs
C
      CALL SETFILE_LIM(ICMODEL,ISDYNPI,EQMODEL,EC_EQRUN,DRAINCLASS,
     &                 NSOIL,IAWC,SOMDEPTH,LUSEQ,NXYEARS,
     &                 DOMSOILC,DOMSOILBD,DOMSOILPH,INMODE,
     &                 DOMSOILCLAY,DOMSOILSILT,DOMSOILSAND,
     &                 DOMSOILISIMP,DOMSOILIMPDEPTH,
     &				 TOTPIC,AVERAIN,AVEPET,AVETEMP,METFILE,WTABLE,
     &                 CACCUM,ANCH4,ANCO2,ANDOC,GRIDLAT,NSOMLAY,PI_SPEC)
      LU1=LUSEQ(1)
C
C ...Set time factors (timestep=monthly)
C
	C_TS=0.0
      N_TS=0
	SECONDS=(365.25/12)*24*60*60
C
C For each land use 
C
      DO 50 P=1,MAXLU+1 !MLR 1+1	!EOJ change so spin up is for landuse in year 1
C
C ... check if this land use is found
C
       LU1=P
        DO IYEAR=1,NXYEARS
	 IF(LUSEQ(IYEAR).EQ.LU1)GOTO 70
	ENDDO
	IF(LU1.EQ.MAXLU+1)THEN 
	LU1=LUSEQ(1) !MAKES SPIN UP DO FIRST ONE LAST
	GOTO 70
	ENDIF
      GOTO 50
C
C... if this land use is found then...
C
70    CONTINUE
  
      IF(INMODE.GE.9.AND.INMODE.LT.10)THEN
          LU1=LUSEQ(1)
      endif
C
C ......get soil characteristics 
C
        ISERIES=1
	  CALL GETSOIL(ISERIES,LU1,SOMDEPTH,NSOMLAY,							
     &             DOMSOILC,DOMSOILCLAY,DOMSOILSILT,				
     &             DOMSOILSAND,DOMSOILBD,DOMSOILPH,
     &             DOMSOILISIMP,DOMSOILIMPDEPTH,
     &             TOC,IOM,CLAY,SILT,SAND,BULKDENS,SOILPH,ISERROR)
	  IF(ISERROR.EQ.1)THEN
	    PRINT*,'ERROR! Series ',ISERIES,' not defined for LU ',LU1
	    IF(ISWAIT.EQ.1)THEN
	      PRINT*,'Press any key to continue...'
            READ(*,*)
	    ENDIF
	    STOP
	  ENDIF
        CALL GET_PLANT_DIST(TOTPIC(LU1),PI_CEQ_MON,LU1)
C
C ......initialise water
C              
	  CALL INIT_GIS_WATER(BULKDENS,CLAY,SAND,FLOWPROP,IAWC,IROCK,
     &                        NSOIL,SILT,TOC,WMAX,WSAT,WTABLE,ROOT,
     &                        AVEPET,AVERAIN,AVETEMP,SOILW,PI_CEQ_MON,
     &                        WILTPOINT,FIELDCAP,SATWATCONT,SECONDS)
C
C ......get plant distribution for this land use
C

C ......get weather data
	DO IMON=1,12
	IF(PI_CEQ_MON(IMON,LU1).GT.0)THEN
	CROP(IMON)=1
	ELSE
	CROP(IMON)=0
	ENDIF
	ENDDO
C
	  CALL GETLTA_LIM(LTA_AWC,LTA_TEMP,WMAX,WSAT,FLOWPROP,
     &                AVERAIN,AVEPET,AVETEMP,CROP)
C
C ......get soil C and N parameters
C
	  CALL GET_SUNDIAL_SOILPARS(NSOIL,TOC,IOM,CRIT,LU1,
     &                          CLAY,BULKDENS,SPARMODEL)
C
C ......get plant distribution for this land use
C
        CALL GET_PLANT_DIST(TOTPIC(LU1),PI_CEQ_MON,LU1)
C
C ......initialise the soil for this land use
C
        IF(INMODEL.EQ.INPASSCN)
     &    CALL GET_CTON_FOEREID(DPMCTON,RPMCTON,LU1)
	  IF(ICMODEL.EQ.ICROTHCEQ)
     &    CALL SET_DPMRPMRATIO(LU1,DPM_RPM_RATIO)
        
        equ_plant=pi_ceq_mon
 !       IF(EQMODEL.EQ.EQJONES)then
                 iyear=1
              CALL INIT_GIS_CROP(IYEAR,LU1,SECONDS,
     &                     TOTPIC(LU1),PNREQ(LU1),
     &                     N_STEPS,ISOWN,NSOW,IHARV,LHARV,IANTHES,
     &                     MEND,IEND,NFERT,IFERT,ILAB,FERT,TFERT,
     &                     JSTOP,CTOT)
              
       do  IK=1,12
              iyear=1
                  CALL RUN1_GIS_CROP(IYEAR,IK,LU1,IEND,MEND,NSOW,
     &                          N_STEPS,IS_TS,TOTPIC(LU1),
     &                          PNREQ(LU1),CACTOT,CATOT15,
     &                          ICOVER,SXORGN,ORGN,
     &                          TRNIN,RNIN,RNIN15,C_TS,TCINP,WR)   
          iyear=0
          if(sum(PI_CEQ_MON(ik,:)).GT.0)then
          Equ_plant(ik,:)= C_TS*PI_CEQ_MON(ik,:)/sum(PI_CEQ_MON(ik,:))      
           else
               Equ_plant(ik,:)=0
           endif
      Enddo
!199          continue    
 !            endif
             PI_CEQ_MON=Equ_plant
             
	  CALL INIT_GIS_SOILCN_NOPARS(ICMODEL,INMODEL,DOCMODEL,
     &                                  EQMODEL,SECONDS,
     &                                  NSOIL,TOC,IOM,HZ1,CRIT,
     &                                  SOILN,SOIL15,AMMN,AMMN15,
     &                                  DPMCARB0,DPMNIT0,DPMNLAB0,
     &                                  RPMCARB0,RPMNIT0,RPMNLAB0,
     &                                  BCARB0,BNIT0,BNLAB0,
     &                                  HCARB0,HNIT0,HNLAB0,
     &                                  MOBDOC,MOBDON,SATWATCONT,
     &                                  CLAY,BULKDENS,Equ_plant,
     &                                  LTA_AWC,LTA_TEMP,ITFUNC,IMFUNC,
     &                                  WMAX,WSAT,ICFACTOR,wiltpoint,
     &                                  DPM_RPM_RATIO,DPMCTON,RPMCTON,
     &                                  PH_MODEL,SOILPH,
     &                                  CACCUM,ANCH4,ANCO2,ANDOC)
                 
             
	  DO 100 IL=1,MAXLAYER1
	    IF(BCARB0(IL).LT.0.OR.HCARB0(IL).LT.0.OR.
     &     RPMCARB0(IL).LT.0.OR.DPMCARB0(IL).LT.0.OR.
     &     BNIT0(IL).LT.0.OR.HNIT0(IL).LT.0.OR.
     &     RPMNIT0(IL).LT.0.OR.DPMNIT0(IL).LT.0.OR.
     &     SOILN(IL).LT.0.OR.AMMN(IL).LT.0)THEN
	      PRINT*,'ERROR! Equilibrium not found!'
	      PRINT*,'    .....Press any key to continue'
	      READ(*,*)
	      STOP
	    ENDIF
100     CONTINUE
C
C ......save annual plant input
C
        CALL SUM_PLANT_INPUT(TOTPIC(LU1),PI_CEQ_MON,LU1)
C
C .....initialise pH calculation
C
        IF(PH_MODEL.EQ.PH_FROM_VSD)THEN
	    START=1
          DO 200 IL=1,MAXLAYER1
	      CALL GETEQUIVALENTS(SOILN(IL),AMMN(IL),MOBDOC(IL),SOILW(IL),
     &                        CNO3IN,CNH4IN,CDOCIN)
            IF(IL.EQ.1)THEN
	        WATERIN=AVERAIN(1)/1000.
	        WATEROUT=DRAINW(IL)/1000.
	      ELSEIF(IL.GT.1)THEN
	        WATERIN=DRAINW(IL-1)/1000.
	        WATEROUT=DRAINW(IL)/1000.
	      ENDIF
	      THICK=MAXDEPTH/(MAXLAYER*100.)
            CALL RUN_VSD(START,CNO3IN,CNH4IN,CDOCIN,
     &                 WATERIN,WATEROUT,THICK,
     &                 SOILPH(IL),CAL(IL),ISERIES,VSD_STATE_IL)
	      ISAVE=1
	      CALL SAVE_VSD_STATE(ISAVE,VSD_STATE_IL,IL)
200       CONTINUE
        ENDIF
C
C ......get rate modifier associated with SUNDIAL routines,
C
        ISPINUP=ISPINUP_OFF
        !ATM = 1.5 ! MLR temp fix
        IF(ISPINUP.EQ.ISPINUP_ON)THEN
          CALL GET_NLIM(AMMN,AMMN15,ATM,AVEPET,AVERAIN,AVETEMP,
     &            BCARB0,BNIT0,BNLAB0,CH4MODEL,CNMODEL,CTOT,
     &            CLOSSX,CLOSSX15,CRIT,DMODEL,DOCMODEL,DPMCARB0,
     &		    DPMNIT0,DPMNLAB0,DRAINW,FERT,FIXN,FLOWPROP,
     &            FYMFERT,FYMFERT15,HCARB0,HNIT0,HNLAB0,IANTHES,
     &            IEND,IFERT,IHARV,ILAB,IMFUNC,ICMODEL,INMODEL,IOLAB,
     &            IOM,IORGM,ISERIES,ISOWN,ITFUNC,IVOL,JORGM,JSTOP,
     &            LHARV,LU1,MEND,NFERT,NORGM,N_STEPS,NSOIL,NSOW,ORGMA,
     &            PNREQ,REFIL,RPMCARB0,RPMNIT0,RPMNLAB0,SECONDS,
     &            SOIL15,SOILN,SOILW,TAM,TFERT,THISFERT,
     &            THISFERT15,TIN,TOC,TOTPIC,VOLAT,VOLAT15,WLEACH,
     &            WLEACH15,WSAT,WMAX,SOILPH,PI_CEQ_MON)		  
	  ENDIF
C
C ......get modification to plant inputs due to full run
C
        IMODPI=IMODPI_OFF
        IF(IMODPI.EQ.IMODPI_ON)THEN
	    CALL GET_MODPI_LAYERED(AMMN,AMMN15,ATM,AVEPET,AVERAIN,AVETEMP,
     &              BCARB0,BNIT0,BNLAB0,CH4MODEL,CNMODEL,CTOT,
     &			  CLOSSX,CLOSSX15,CRIT,DMODEL,DOCMODEL,DPMCARB0,
     &			  DPMNIT0,DPMNLAB0,DRAINW,FERT,FIXN,FLOWPROP,
     &              FYMFERT,FYMFERT15,HCARB0,HNIT0,HNLAB0,IANTHES,
     &              IEND,IFERT,IHARV,ILAB,IMFUNC,ICMODEL,INMODEL,IOLAB,
     &              IOM,IORGM,ISERIES,ISOWN,ITFUNC,IVOL,JORGM,JSTOP,
     &              LHARV,LU1,MEND,NFERT,NORGM,N_STEPS,NSOIL,NSOW,ORGMA,
     &              PNREQ,REFIL,RPMCARB0,RPMNIT0,RPMNLAB0,SECONDS,
     &              SOIL15,SOILN,SOILW,TAM,TFERT,THISFERT,
     &              THISFERT15,TIN,TOC,TOTPIC,VOLAT,VOLAT15,WLEACH,
     &              WLEACH15,WSAT,WMAX,SOILPH,NXYEARS,EQMODEL,CLAY,
     &              PI_CEQ_MON,DOMSOILISIMP,DOMSOILIMPDEPTH,
     &              SOMDEPTH,NSOMLAY,LTA_AWC,LTA_TEMP,PH_MODEL)
	  ENDIF
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
     &              JORGM,JSTOP,LHARV,LU1,MEND,NFERT,NORGM,
     &              NSOIL,NSOW,ORGMA,PNREQ,REFIL,RPMCARB0,
     &              RPMNIT0,RPMNLAB0,SOIL15,SOILN,SOILW,TAM,
     &              TFERT,THISFERT,THISFERT15,TIN,TOC,TOTPIC,VOLAT,
     &              VOLAT15,WLEACH,WLEACH15,WSAT,WMAX,SOILPH,NXYEARS,
     &              EQMODEL,CLAY,PI_CEQ_MON,DOMSOILISIMP,
     &              DOMSOILIMPDEPTH,SOMDEPTH,NSOMLAY,LTA_AWC,LTA_TEMP,
     &              PH_MODEL)
	  ENDIF
C
C ......get modification to decomposition due to full run
C
        IF(IMODDEC.EQ.IMODDEC_ON)THEN
	   CALL GET_MODDEC_LAYERED(AMMN,AMMN15,ATM,AVEPET,AVERAIN,AVETEMP,
     &              BCARB0,BNIT0,BNLAB0,CH4MODEL,CNMODEL,CTOT,
     &			  CLOSSX,CLOSSX15,CRIT,DMODEL,DOCMODEL,DPMCARB0,
     &			  DPMNIT0,DPMNLAB0,DRAINW,FERT,FIXN,FLOWPROP,
     &              FYMFERT,FYMFERT15,HCARB0,HNIT0,HNLAB0,IANTHES,
     &              IEND,IFERT,IHARV,ILAB,IMFUNC,ICMODEL,INMODEL,IOLAB,
     &              IOM,IORGM,ISERIES,ISOWN,ITFUNC,IVOL,JORGM,JSTOP,
     &              LHARV,LU1,MEND,NFERT,NORGM,N_STEPS,NSOIL,NSOW,ORGMA,
     &              PNREQ,REFIL,RPMCARB0,RPMNIT0,RPMNLAB0,SECONDS,
     &              SOIL15,SOILN,SOILW,TAM,TFERT,THISFERT,
     &              THISFERT15,TIN,TOC,TOTPIC,VOLAT,VOLAT15,WLEACH,
     &              WLEACH15,WSAT,WMAX,SOILPH,NXYEARS,EQMODEL,CLAY,
     &              PI_CEQ_MON,DOMSOILISIMP,DOMSOILIMPDEPTH,
     &              SOMDEPTH,NSOMLAY,LTA_AWC,LTA_TEMP,PH_MODEL)		  
	  ENDIF
C
C ... Save state of soil for this land use
C
        ISAVE=1 ! Save soil 
        CALL SAVE_GIS_SOIL(ISAVE,LU1,
     &                     SOILN,SOIL15,AMMN,AMMN15,
     &                     SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
     &                     DPMCARB0,DPMNIT0,DPMNLAB0,
     &                     RPMCARB0,RPMNIT0,RPMNLAB0,
     &                     BCARB0,BNIT0,BNLAB0,
     &                     HCARB0,HNIT0,HNLAB0,
     &                     IOM,PI_CEQ_MON)
      
50    CONTINUE
C
C ... Retrieve state of soil for first land use
C
      LU1=LUSEQ(1)
	THISLU=LUSEQ(1)
	ISERIES=1
	CALL GETSOIL(ISERIES,LU1,SOMDEPTH,NSOMLAY,							
     &             DOMSOILC,DOMSOILCLAY,DOMSOILSILT,				
     &             DOMSOILSAND,DOMSOILBD,DOMSOILPH,
     &             DOMSOILISIMP,DOMSOILIMPDEPTH,
     &             TOC,IOM,CLAY,SILT,SAND,BULKDENS,SOILPH,ISERROR)
      ISAVE=0 ! Retrieve soil 
      CALL SAVE_GIS_SOIL(ISAVE,LU1,
     &                   SOILN,SOIL15,AMMN,AMMN15,
     &                   SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
     &                   DPMCARB0,DPMNIT0,DPMNLAB0,
     &                   RPMCARB0,RPMNIT0,RPMNLAB0,
     &                   BCARB0,BNIT0,BNLAB0,
     &                   HCARB0,HNIT0,HNLAB0,
     &                   IOM,PI_CEQ_MON)
    
	IYEAR=1
      lastlu = thislu
      lu_change_yr = 0

      DO WHILE (IYEAR.LE.NXYEARS)
      THISLU=LUSEQ(IYEAR)
	  
	  CALL GET_PLANT_DIST(TOTPIC(THISLU),PI_CEQ_MON,THISLU) !added by eoj, 6/3/12, adds plant carbon each year for specific land use

C
C Initialise soil and crop in first year
C	
        IF(IYEAR.EQ.1)THEN
C
C ...initialise the soil for initial land use,
C

			CALL SET_DPMRPMRATIO(LU1,DPM_RPM_RATIO)
              
            
	      CALL GETSOIL(ISERIES,LU1,SOMDEPTH,NSOMLAY,							
     &                 DOMSOILC,DOMSOILCLAY,DOMSOILSILT,				
     &                 DOMSOILSAND,DOMSOILBD,DOMSOILPH,
     &                 DOMSOILISIMP,DOMSOILIMPDEPTH,
     &                 TOC,IOM,CLAY,SILT,SAND,BULKDENS,SOILPH,ISERROR)   
  
          CALL INIT_GIS_SOILCN_NOOPT(DOCMODEL,
     &                               EQMODEL,SECONDS,
     &                               NSOIL,TOC,IOM,HZ1,CRIT,
     &                               SOILN,SOIL15,AMMN,AMMN15,
     &                               DPMCARB0,DPMNIT0,DPMNLAB0,
     &                               RPMCARB0,RPMNIT0,RPMNLAB0,
     &                               BCARB0,BNIT0,BNLAB0,
     &                               HCARB0,HNIT0,HNLAB0,
     &                               MOBDOC,MOBDON,THISLU,
     &                               CLAY,BULKDENS,PI_CEQ_MON,
     &                               LTA_AWC,LTA_TEMP,ITFUNC,IMFUNC,
     &                               WMAX,WSAT,ICFACTOR,PH_MODEL,SOILPH)
          IF(ISPINUP.EQ.ISPINUP_ON)THEN
            CALL SET_NLIM(THISLU,ICFACTOR,
     &                  SOILN,SOIL15,AMMN,AMMN15,SOILW,
     &                  DPMCARB0,DPMNIT0,DPMNLAB0,
     &                  RPMCARB0,RPMNIT0,RPMNLAB0,
     &                  BCARB0,BNIT0,BNLAB0,
     &                  HCARB0,HNIT0,HNLAB0,IOM,PI_CEQ_MON)
	    ENDIF
C
C Set up results arrays
C
          CALL STARTRES(IRYEAR,IK,RSTART,DPMNIT0,RPMNIT0,
     &                  BSTART,BNIT0,HNIT0,
     &                  RSTART15,DPMNLAB0,RPMNLAB0,
     &                  BSTART15,BNLAB0,HNLAB0,
     &                  TOTAST,AMMN,TOTAST15,AMMN15,
     &                  TOTNST,SOILN,TOTNST15,SOIL15)
          CALL SETFACTORS(NSOIL,THISLU,INPUTFACTORS,ATM)
        ENDIF
C
C Initialise current crop
      IF(IYEAR.GT.1.)THEN
      IF(((LUSEQ(IYEAR).NE.LUSEQ(IYEAR-1)).or.miscticker.eq.20).AND.
     &                (LUSEQ(IYEAR).EQ.5.or.LUSEQ(IYEAR).EQ.6
     &                .or.LUSEQ(IYEAR).EQ.9))THEN
	    MISCTICKER=1
      endif
      endif
C 
      IF(PI_SPEC(IYEAR).GT.0)Then         !if not plant input specified use estimated land use specific miami
		TOTPIC(THISLU)=PI_SPEC(IYEAR)
          ! MLR: quick fix for ELUM runs to change above-ground yield into plant input:
          ! Convert above-ground biomass to total biomass by multiplying by rootshoot
          ! Convert DM to C by multiplying by 0.5
          ! Convert total biomass to plant inputs to soil by multiplying by frac
          TOTPIC(THISLU) = TOTPIC(THISLU) * rootshoot(thislu) * 0.5 *
     &                     frac(thislu)
          ! MLR end
		PIANN=TOTPIC(THISLU)  
          
	ELSE
	    ANNRAIN=0
	    ANNTEMP=0
          rewind(3)
	    DO 660 IMON=1,12 
              CALL GETWEATHER_LIM(ISOPEN,IMON,LHARV,AVERAIN,AVEPET,
     &                        AVETEMP,METFILE,GRIDLAT,
     &                        RAIN,EVAPW,SOILTEMP,AIRTEMP,IYEAR)
              ANNRAIN=ANNRAIN+RAIN
	        ANNTEMP=ANNTEMP+AIRTEMP
660       CONTINUE
	    ANNTEMP=ANNTEMP/12 
	    IF(INMODE.LT.6.OR.INMODE.GE.7)THEN
	        CALL MIAMI_DYCE(THISLU,ANNTEMP,ANNRAIN,PIANN)
	        TOTPIC(THISLU)=PIANN
	    ELSE
	        PIANN=TOTPIC(THISLU)
	    ENDIF
          rewind(3)
      endif

      call adj_pi_by_lu_age(iyear, ik, thislu, lastlu, lu_change_yr,
     &                      piann)
      
      CALL INIT_GIS_CROP_LIM(IYEAR,THISLU,SECONDS,
     &                     TOTPIC(THISLU),PNREQ(THISLU),
     &                     N_STEPS,ISOWN,NSOW,IHARV,LHARV,IANTHES,
     &                     MEND,IEND,NFERT,IFERT,ILAB,FERT,TFERT,
     &                     JSTOP,CTOT,miscticker)
	 
      IF(ISPINUP.EQ.ISPINUP_ON)THEN
	    CALL SET_NLIM_RATE(THISLU,ICFACTOR)
	ENDIF
C
C Cultivate soil for gra->ara, for->ara, for->gra, nat->ara, for->gra, nat->for
C
      Miscticker=(miscticker+1) 
          if (iyear > 1) then
              cultivate = .false.
              if (luseq(iyear) /= luseq(iyear -1)) then
                  ! Assume that no "additional" cultivation required to
                  ! convert arable, oilseed rape and beet to other land uses
                  
                  ! Grass to...
                  if (luseq(iyear - 1) == 2 .or. luseq(iyear - 1) == 4) then
c                      if (luseq(iyear) == 1) .or.  ! ...arable
c     &                   (luseq(iyear) == 5) .or.  ! ...miscanthus
c     &                   (luseq(iyear) == 7) .or.  ! ...sugar beet
c     &                   (luseq(iyear) == 8) .or.  ! ...oilseed rape
c                              then
                      cultivate = .true.
                      cultdepth = 30
                      vigour = 0.25
                      invert = 1                     

                  ! Forest to...
                  elseif (luseq(iyear - 1) == 3) then
                      cultivate = .true.
                      cultdepth = 50
                      vigour = 0.25
                      invert = 0                                      

                  ! Natural to...
c                  elseif (luseq(iyear - 1) == 4) then
c                      cultivate = .true.
c                      cultdepth = 30
c                      vigour = 0.25
c                      invert = 1                                      
                      
                  ! Miscanthus to...
                  elseif (luseq(iyear - 1) == 5) then
                      cultivate = .true.
                      cultdepth = 30
                      vigour = 0.25
                      invert = 1                                      

                  ! Short Rotation Coppice to
                  elseif (luseq(iyear - 1) == 6) then
                      cultivate = .true.
                      cultdepth = 30
                      vigour = 0.25
                      invert = 1                                      

                  ! Short Rotation Forest to...
                  elseif (luseq(iyear - 1) == 9) then
                      cultivate = .true.
                      cultdepth = 40
                      vigour = 0.25
                      invert = 1                                   
                  endif
              endif
              if (cultivate) then
                  call cultiv(luseq(iyear - 1),
     &                   dpmcarb0, dpmnit0, dpmnlab0,
     &                   rpmcarb0, rpmnit0, rpmnlab0,
     &                   bcarb0, bnit0, bnlab0,
     &                   hcarb0, hnit0, hnlab0,
     &                   cultdepth, vigour, invert)              
              endif
          endif
!          IF(
!!     &        (LUSEQ(IYEAR-1).EQ.1.AND.LUSEQ(IYEAR).EQ.1).OR.
!     &        (LUSEQ(IYEAR-1).EQ.7.AND.LUSEQ(IYEAR).EQ.1).OR.
!     &        (LUSEQ(IYEAR-1).EQ.8.AND.LUSEQ(IYEAR).EQ.1).OR.     
!     &	    (LUSEQ(IYEAR-1).EQ.7.AND.LUSEQ(IYEAR).EQ.7).OR.
!     &        (LUSEQ(IYEAR-1).EQ.8.AND.LUSEQ(IYEAR).EQ.8))THEN  !SRF->SRF
!	         CULTDEPTH=30
!	         VIGOUR=0.2
!	         INVERT=1
!                          	   CALL CULTIV(LUSEQ(IYEAR-1),
!     &                   DPMCARB0,DPMNIT0,DPMNLAB0,
!     &                   RPMCARB0,RPMNIT0,RPMNLAB0,
!     &                   BCARB0,BNIT0,BNLAB0,
!     &                   HCARB0,HNIT0,HNLAB0, 
!     &                   CULTDEPTH,VIGOUR,INVERT)
!      ENDif
!      endif
!	
!       IF(IYEAR.GT.1.)THEN
!	     IF((LUSEQ(IYEAR).NE.LUSEQ(IYEAR-1)).OR.(MISCTICKER.EQ.20.
!     &	   AND.(LUSEQ(IYEAR).EQ.5.or.LUSEQ(IYEAR).EQ.6)))THEN
!	         CULTDEPTH=30 
!	         VIGOUR=0.5
!               INVERT=1
!              
!c                         	   CALL CULTIV(LUSEQ(IYEAR-1),
!c     &                   DPMCARB0,DPMNIT0,DPMNLAB0,
!c     &                   RPMCARB0,RPMNIT0,RPMNLAB0,
!c     &                   BCARB0,BNIT0,BNLAB0,
!c     &                   HCARB0,HNIT0,HNLAB0, 
!c     &                   CULTDEPTH,VIGOUR,INVERT)
!                           !  Endif
!	     ELSEIF(
!     &	     (LUSEQ(IYEAR-1).EQ.2.AND.LUSEQ(IYEAR).EQ.1).OR.
!     &         (LUSEQ(IYEAR-1).EQ.3.AND.LUSEQ(IYEAR).EQ.1).OR.
!     &         (LUSEQ(IYEAR-1).EQ.4.AND.LUSEQ(IYEAR).EQ.1).OR.
!     &         (LUSEQ(IYEAR-1).EQ.5.AND.LUSEQ(IYEAR).EQ.1).OR.
!     &         (LUSEQ(IYEAR-1).EQ.6.AND.LUSEQ(IYEAR).EQ.1).OR.
!     &         (LUSEQ(IYEAR-1).EQ.4.AND.LUSEQ(IYEAR).EQ.2).OR.  !For->SRC/F
!     &         (LUSEQ(IYEAR-1).EQ.4.AND.LUSEQ(IYEAR).EQ.3).OR.  !For->Misc
!     &         (LUSEQ(IYEAR-1).EQ.4.AND.LUSEQ(IYEAR).EQ.4).OR.  !For->Misc
!     &         (LUSEQ(IYEAR-1).EQ.4.AND.LUSEQ(IYEAR).EQ.5).OR.  !For->Misc
!     &         (LUSEQ(IYEAR-1).EQ.4.AND.LUSEQ(IYEAR).EQ.6))THEN
!	         CULTDEPTH=50
!	         VIGOUR=0.5
!	         INVERT=1
!           ELSEIF(
!     &              (LUSEQ(IYEAR-1).EQ.3.AND.LUSEQ(IYEAR).EQ.2).OR.
!     &			  (LUSEQ(IYEAR-1).EQ.3.AND.LUSEQ(IYEAR).EQ.6).OR.  !Gra->SRC/F
!     &			  (LUSEQ(IYEAR-1).EQ.3.AND.LUSEQ(IYEAR).EQ.5).OR.  !Gra->Misc
!     &	          (LUSEQ(IYEAR-1).EQ.5.AND.LUSEQ(IYEAR).EQ.3).OR.  !Misc->Gra   
!     &	          (LUSEQ(IYEAR-1).EQ.5.AND.LUSEQ(IYEAR).EQ.5).OR.  !Misc->Misc
!     &              (LUSEQ(IYEAR-1).EQ.3.AND.LUSEQ(IYEAR).EQ.4))THEN
!	         CULTDEPTH=30
!	         VIGOUR=0.5
!	         INVERT=1
!	     ELSEIF(
!     &              (LUSEQ(IYEAR-1).EQ.4.AND.LUSEQ(IYEAR).EQ.3).OR.    
!     &	          (LUSEQ(IYEAR-1).EQ.5.AND.LUSEQ(IYEAR).EQ.6))THEN  !SRF->SRF
!	         CULTDEPTH=10
!	         VIGOUR=0.1
!	         INVERT=1
!           ELSEIF(
!     &              (LUSEQ(IYEAR-1).EQ.1.AND.LUSEQ(IYEAR).EQ.2).OR.     
!     &	          (LUSEQ(IYEAR-1).EQ.1.AND.LUSEQ(IYEAR).EQ.3).OR.     
!     &	          (LUSEQ(IYEAR-1).EQ.1.AND.LUSEQ(IYEAR).EQ.4).OR.     
!     &	          (LUSEQ(IYEAR-1).EQ.1.AND.LUSEQ(IYEAR).EQ.5).OR.     
!     &	          (LUSEQ(IYEAR-1).EQ.1.AND.LUSEQ(IYEAR).EQ.6))THEN  !SRF->SRF
!	         CULTDEPTH=5
!	         VIGOUR=0.01
!	         INVERT=0
!           ELSEIF(
!     &              (LUSEQ(IYEAR-1).EQ.2.AND.LUSEQ(IYEAR).EQ.1).OR.     
!     &	          (LUSEQ(IYEAR-1).EQ.2.AND.LUSEQ(IYEAR).EQ.3).OR.     
!     &	          (LUSEQ(IYEAR-1).EQ.2.AND.LUSEQ(IYEAR).EQ.4).OR.     
!     &	          (LUSEQ(IYEAR-1).EQ.2.AND.LUSEQ(IYEAR).EQ.5).OR.     
!     &	          (LUSEQ(IYEAR-1).EQ.2.AND.LUSEQ(IYEAR).EQ.6))THEN  !SRF->SRF
!	         CULTDEPTH=30
!	         VIGOUR=0.3
!	         INVERT=1
!               
!
!           	   CALL CULTIV(LUSEQ(IYEAR-1),
!     &                   DPMCARB0,DPMNIT0,DPMNLAB0,
!     &                   RPMCARB0,RPMNIT0,RPMNLAB0,
!     &                   BCARB0,BNIT0,BNLAB0,
!     &                   HCARB0,HNIT0,HNLAB0, 
!     &                   CULTDEPTH,VIGOUR,INVERT)
!	    ENDIF	         
!      ENDIF
C If annual plant inputs are specified, and mode is not 7 (adjust PI due to weather conditions), set PI to specified value
C
!        IF(PI_SPEC(IYEAR).GT.0)Then         !if not plant input specified use estimated land use specific miami
!		TOTPIC(THISLU)=PI_SPEC(IYEAR)
!		PIANN=TOTPIC(THISLU)
!			  ELSE
!     
!	ANNRAIN=0
!	ANNTEMP=0
!	
!      rewind(3)
!	DO 660 IMON=1,12
!          
!         CALL GETWEATHER_LIM(ISOPEN,IMON,LHARV,AVERAIN,AVEPET,AVETEMP, 
!     &                        METFILE,GRIDLAT,
!     &                        RAIN,EVAPW,SOILTEMP,AIRTEMP,IYEAR)
!        ANNRAIN=ANNRAIN+RAIN
!	  ANNTEMP=ANNTEMP+AIRTEMP
!660   CONTINUE
!	ANNTEMP=ANNTEMP/12 
!	IF(INMODE.LT.6.OR.INMODE.GE.7)THEN
!	  CALL MIAMI_DYCE(THISLU,ANNTEMP,ANNRAIN,PIANN)
!	  TOTPIC(THISLU)=PIANN
!	ELSE
!	PIANN=TOTPIC(THISLU)
!	ENDIF
!	
!
!      rewind(3)
!	endif
C
C For each week till end of growing season...
C           
        DO 700 IK=1,12  !MEND  Runs for whole year not just till end of harvest
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
C
C Set soil temperature to air temperature
C
          CALL GETWEATHER_LIM(ISOPEN,IK,LHARV,AVERAIN,AVEPET,AVETEMP, 
     &                        METFILE,GRIDLAT,
     &                        RAIN,EVAPW,SOILTEMP,AIRTEMP,IYEAR)
C
C Adjust plant inputs for dynamic weather conditions
C
       !   PIANN=TOTPIC(THISLU)
 !         IF(ISDYNPI.EQ.ISDYNPI_ON)
  !   &      CALL ADJUST_PI(AVERAIN,AVETEMP,IK,LHARV,PIANN,RAIN,SOILTEMP)
C
C If current week is first in the growing season set results
C
       CALL SETRES(BALANCE,CO2,SX,SXORGN)
       
C
C ...Calculate crop C and N returns and N offtake 
C
	    CALL RUN1_GIS_CROP(IYEAR,IK,THISLU,IEND,MEND,NSOW,
     &                          N_STEPS,IS_TS,PIANN,
     &                          PNREQ(THISLU),CACTOT,CATOT15,
     &                          ICOVER,SXORGN,ORGN,
     &                          TRNIN,RNIN,RNIN15,C_TS,TCINP,WR)  
          
          IF(C_TS.GT.0)Then
              ICOVER=1
          else
              icover=0
          endif
          
           CALL RUN_SUNDIAL_WATER_SUBMONTHLY(RAIN,WMAX,
     &                                      SOILW,WSAT,FLOWPROP,
     &                                      DRAINW,REFIL,EVAPW,
     &                                      SOILTEMP(1),ICOVER)

C
C ...Add stuff to soil
C
          CALL GETFYM(IK,NORGM,IORGM,ORGMA,JORGM,IOLAB,
     &                ORGMANF,JORGMNF,IOLABNF) 
          CALL GETFERT(IK,NFERT,IFERT,FERT,TFERT,ILAB,IVOL,
     &                FERTNF,TFERTNF,ILABNF,IVOLNF) 
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
c          CALL PARTITION_PLANTINPUT(C_TS,RNIN,RNIN15,
c     &                                 PI_C,PI_N,PI_N15)
          CALL ADD_PLANT_DIST(C_TS,RNIN,RNIN15,PI_CEQ_MON,IK,
     &                        PI_C,PI_N,PI_N15) 
          IF(INMODEL.EQ.INPASSCN)
     &      CALL GET_CTON_FOEREID(DPMCTON,RPMCTON,THISLU)
	    CALL SET_DPMRPMRATIO(THISLU,DRRAT)
c	    IF(DOMSOILISIMP(ISERIES,THISLU))
c     &      CALL LIMIT_DRAIN(DOMSOILIMPDEPTH(ISERIES,THISLU),SOILW)
          MEASLAY=SOMDEPTH(ISERIES,LU1,NSOMLAY(ISERIES,LU1))
    	    MEASLAY=MEASLAY*MAXLAYER1/MAXDEPTH
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
     &					             DPMCARB,RPMCARB,BCARB,HCARB,
     &					             WRATEDM,TRATEM,
     &                                 SOILPH,MEASLAY)
C
C Run within season crop routines
C
          CALL RUN2_GIS_CROP(SOILTEMP(1),AMMN,AMMN15,ATM,ATM15,
     &                         CACT,CACT15,CACTOT,CATOT15,CLOSSX,
     &                         CLOSSX15,PNREQ(THISLU),CRIT,CTOT,
     &                         CUPTN,IK,IS_TS,IYEAR,JSTOP,MEND,NSOIL,
     &                         NSOW,OCROPN,ORGN,RNIN,SOIL15,SOILN,SRNIN,
     &                         SXORGN,TACTOT,TAM,TC,TIN,TRNIN,VOLAT,
     &                         VOLAT15,XORGN,SEEDIN,PLANTUP,DDAYS)	
C
C Leaching routines
C

          CALL PHYSICAL_SUNDIAL_SOILCN(DOCMODEL,SECONDS,
     &                                   NSOIL,SOILTEMP,
     &                                   WMAX,SOILW,DRAINW,REFIL,
     &                                   SLEACH,SLEA15,
     &                                   WLEACH,WLEACH15,CONC,CONC15,
     &                                   SOILN,SOIL15,AMMN,AMMN15,
     &                                   TIN,TAM,MOBDOC,MOBDON,
     &                                   LEACHDOC,LEACHDON)
C
C pH calculation
C
          IF(PH_MODEL.EQ.PH_FROM_VSD)THEN
	      START=0
            DO 401 IL=1,MAXLAYER1
	        ISAVE=0
	        CALL SAVE_VSD_STATE(ISAVE,VSD_STATE_IL,IL)
	        CALL GETEQUIVALENTS(SOILN(IL),AMMN(IL),MOBDOC(IL),
     &                            SOILW(IL),CNO3IN,CNH4IN,CDOCIN)
              IF(IL.EQ.1)THEN
		      WATERIN=AVERAIN(IK)/1000.
	          WATEROUT=DRAINW(IL)/1000.
	        ELSEIF(IL.GT.1)THEN
	          WATERIN=DRAINW(IL-1)/1000.
	          WATEROUT=DRAINW(IL)/1000.
	        ENDIF
	        THICK=MAXDEPTH/(MAXLAYER*100.)
              CALL RUN_VSD(START,CNO3IN,CNH4IN,CDOCIN,
     &                     WATERIN,WATEROUT,THICK,
     &                     SOILPH(IL),CAL(IL),ISERIES,VSD_STATE_IL)
	        ISAVE=1
	        CALL SAVE_VSD_STATE(ISAVE,VSD_STATE_IL,IL)
401         CONTINUE
          ENDIF
C
          IF (FULL_OUTPUT) THEN
            CALL CBALANCE(IYEAR,CO2,CH4,CH4TOAIR,IK,TFYMC,N_TS,TCINP,
     &                 HCARB0,BCARB0,DPMCARB0,RPMCARB0,IOM,PIANN,C_TS)
	    ENDIF

C Output results in normal mode only
C
          IF(IYEAR.GT.0)
     &      CALL TEST1_RES(IYEAR,IRYEAR,NXYEARS,LHARV,
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
     &                   SOMDEPTH(ISERIES,LU1,NSOMLAY(ISERIES,LU1)),
     &                   PLANTUP,FULL_OUTPUT,SUMMARY_OUTPUT)
700     CONTINUE
C
C Test the whether simulated C matches measured C in no steady state run
C
        IYEAR=IYEAR+1
C
C Go back and set up the next crop
C
        END DO
C 
C Close Channels for file input/output
C
	CALL CLOSECHAN2()
C
C Normal termination message
C
      WRITE(*,*)'****************************************************'
      WRITE(*,*)'SIMULATION SUCCESSFULLY COMPLETED!!!'
      WRITE(*,*)'****************************************************'
      IF(ISWAIT.EQ.1)THEN
	  WRITE(*,*)'                       ....Press any key to continue'
	  READ(*,*)
	ENDIF
	END

C
C-------------------------------------------------------------
C
      SUBROUTINE CHECK_PROFILE(
C Inputs: Land use and number of SOM layers
     &                         ILU,NSOMLAY,
C     ... Profile information     
     &                         SOMDEPTH,DOMSOILC,DOMSOILCLAY,
     &                         DOMSOILIMPDEPTH,DOMSOILISIMP,DOMSOILSILT,
     &                         DOMSOILSAND,DOMSOILBD,DOMSOILPH)

C
C Subroutine to check profile information
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXSOMLAY			! Max.no.of soil layers entered			
	PARAMETER (MAXSOMLAY=10)

	INTEGER ISOMLAY				! Current SOM layer
	REAL TEMP					! Temporary real number
C
C Variables passed to/from this subroutine
C
	INTEGER ILU					! IN: Land use code
      INTEGER ISERROR				! IN:Code for error in soil 1=error 0=noerror
      INTEGER ISMISS				! IN:Code for missing values in soil 1=missing 0=nomissing
	REAL DOMSOILC(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! Soil C in 5 
												! major soil series under different LU 
												! in SOM layers (kgC/ha)
	REAL DOMSOILCLAY(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! % Soil clay in 5 
												! major soil series under different LU 
	REAL DOMSOILIMPDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN: Depth of impermeable layer if present (cm)
	INTEGER DOMSOILISIMP(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN: Flag for impermeable layer (0=No; 1=Yes)
	REAL DOMSOILSILT(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! % Soil silt in 5 
												! major soil series under different LU 
	REAL DOMSOILSAND(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! % Soil sand in 5 
												! major soil series under different LU 
	REAL DOMSOILBD(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! Soil BD in 5 
												! major soil series under different LU 
												! in SOM layers (g/cm3)
	REAL DOMSOILPH(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! OUT:Soil PH in 5 
												! major soil series under different LU 
												! in SOM layers (g/cm3)
	INTEGER NSOMLAY(CELLSTOGRID*MAXSERIES/2,MAXLU)	! OUT:Number of SOM layers
	REAL SOMDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN:depth of soil organic matter layers
C
C Check profile information
C
      DO 100 ISOMLAY=1,NSOMLAY(1,ILU)
        if (ilu < 7) then
        IF(ISOMLAY.GT.1)THEN
          CALL CHECK_DOMSOIL(ISERROR,ISMISS,
     &                       SOMDEPTH(1,ILU,ISOMLAY),
     &                       SOMDEPTH(1,ILU,ISOMLAY-1),
     &                       DOMSOILC(1,ILU,ISOMLAY),
     &                       DOMSOILCLAY(1,ILU,ISOMLAY),
     &                       DOMSOILSILT(1,ILU,ISOMLAY),
     &                       DOMSOILSAND(1,ILU,ISOMLAY),
     &                       DOMSOILBD(1,ILU,ISOMLAY),
     &                       DOMSOILPH(1,ILU,ISOMLAY),
     &                       DOMSOILISIMP(1,ILU),
     &                       DOMSOILIMPDEPTH(1,ILU))
	    IF(ISMISS.EQ.1)THEN
	      IF(SOMDEPTH(1,ILU,ISOMLAY).LT.
     &	     SOMDEPTH(1,ILU,ISOMLAY-1))
     &        SOMDEPTH(1,ILU,ISOMLAY)=MAXDEPTH
	      IF(ISOMLAY.EQ.1)TEMP=0
	      IF(ISOMLAY.GT.1)TEMP=SOMDEPTH(1,ILU,ISOMLAY-1)
            CALL GET_DOMSOIL_LAYER(DOMSOILC(1,ILU,ISOMLAY),
     &                             DOMSOILC(1,ILU,ISOMLAY-1),
     &                             DOMSOILCLAY(1,ILU,ISOMLAY),
     &                             DOMSOILCLAY(1,ILU,ISOMLAY-1),
     &                             DOMSOILSILT(1,ILU,ISOMLAY),
     &                             DOMSOILSILT(1,ILU,ISOMLAY-1),
     &                             DOMSOILSAND(1,ILU,ISOMLAY),
     &                             DOMSOILSAND(1,ILU,ISOMLAY-1),
     &                             DOMSOILBD(1,ILU,ISOMLAY),
     &                             DOMSOILBD(1,ILU,ISOMLAY-1),
     &                             DOMSOILPH(1,ILU,ISOMLAY),
     &                             DOMSOILPH(1,ILU,ISOMLAY-1),
     &                             SOMDEPTH(1,ILU,ISOMLAY),
     &                             TEMP)
	      ISMISS=0
	      ISERROR=0
	    ENDIF
        ELSEIF(ISOMLAY.LE.1)THEN
	    TEMP=0
          CALL CHECK_DOMSOIL(ISERROR,ISMISS,
     &                       SOMDEPTH(1,ILU,ISOMLAY),
     &                       TEMP,
     &                       DOMSOILC(1,ILU,ISOMLAY),
     &                       DOMSOILCLAY(1,ILU,ISOMLAY),
     &                       DOMSOILSILT(1,ILU,ISOMLAY),
     &                       DOMSOILSAND(1,ILU,ISOMLAY),
     &                       DOMSOILBD(1,ILU,ISOMLAY),
     &                       DOMSOILPH(1,ILU,ISOMLAY),
     &                       DOMSOILISIMP(1,ILU),
     &                       DOMSOILIMPDEPTH(1,ILU))
        ENDIF
        endif
	  IF(ISERROR.EQ.1)THEN
	    PRINT*,'Error in profile data. LU ',ILU,' layer ',ISOMLAY
	    PRINT*,'Press any key to continue'
	    READ(*,*)
	    STOP
	  ENDIF
	  IF(ISOMLAY.EQ.NSOMLAY(1,ILU).AND.
     &     SOMDEPTH(1,ILU,ISOMLAY).LT.MAXDEPTH)THEN
          SOMDEPTH(1,ILU,ISOMLAY+1)=MAXDEPTH
          CALL GET_DOMSOIL_LAYER(DOMSOILC(1,ILU,ISOMLAY+1),
     &                           DOMSOILC(1,ILU,ISOMLAY),
     &                           DOMSOILCLAY(1,ILU,ISOMLAY+1),
     &                           DOMSOILCLAY(1,ILU,ISOMLAY),
     &                           DOMSOILSILT(1,ILU,ISOMLAY+1),
     &                           DOMSOILSILT(1,ILU,ISOMLAY),
     &                           DOMSOILSAND(1,ILU,ISOMLAY+1),
     &                           DOMSOILSAND(1,ILU,ISOMLAY),
     &                           DOMSOILBD(1,ILU,ISOMLAY+1),
     &                           DOMSOILBD(1,ILU,ISOMLAY),
     &                           DOMSOILPH(1,ILU,ISOMLAY+1),
     &                           DOMSOILPH(1,ILU,ISOMLAY),
     &                           SOMDEPTH(1,ILU,ISOMLAY),
     &                           TEMP)
	  ENDIF
100   CONTINUE
C
C Leave CHECK_PROFILE
C
      END
C
C-----------------------------------------------------------
C
      SUBROUTINE GETWEATHER_LIM(ISOPEN,IK,LHARV,AVERAIN,AVEPET,AVETEMP, 
     &                           METFILE,GRIDLAT,
     &                           RAIN,EVAPW,SOILTEMP,AIRTEMP,IYEAR)
C
C Subroutine to get this months weather for GIS run
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXGROW				! Max.no.of growing seasons allowed
	PARAMETER (MAXGROW=300)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/

	INTEGER IL,IYEAR			! Local counter for months, year counter
	INTEGER IMON				! Month counter
	INTEGER ISOPEN(MAXGROW)		! Is the metfile open (1=yes; 0=no)
	REAL*8 RMON					! Working for current month
	INTEGER THISMON				! Current month
	INTEGER THISYEAR			! Current year
C
C Variables passed from calling subroutine
C
	REAL AIRTEMP				! IN: Air temperature (deg.C/timestep)	
      REAL AVEPET(12)				! Long term average PET (mm/month)
      REAL AVERAIN(12)			! Long term average rainfall(mm/month)
      REAL AVETEMP(12)			! Long term average monthly average 
	REAL EVAPW					! IN: Potential evap. (mm/timestep)
      REAL GRIDLAT				! IN:Average latitude of this 20km2 grid cell
	INTEGER IK					! IN:No.timesteps from prev.harv. to current
	INTEGER LHARV				! IN:No.timesteps from 01/01 to prev. harvest
	CHARACTER*40 METFILE(MAXGROW)	! IN:Met.file  
	REAL RAIN					! IN: Rainfall (mm/timestep)
	REAL SOILTEMP(MAXLAYER)		! IN: Soil temperature (deg.C/timestep)
C
C Get the current month
C
	THISMON=AINT((IK+LHARV)/12.)
	RMON=(REAL(IK+LHARV)/12.)
	THISMON=NINT(12.*(RMON-THISMON))
	IF(THISMON.EQ.0)THISMON=12
	THISMON=IK  !easy way of doing month?
C
C Get current year
C
	THISYEAR=1+AINT((IK+LHARV-1)/12.)
	THISYEAR=IYEAR !Easier to use year counter?
	IF(THISYEAR.GT.MAXGROW)THEN
	  PRINT*,'Error in number of growing year included.'
	  PRINT*,'Maximum allowed = ',MAXGROW
	  READ(*,*)
	  STOP
	ENDIF
C
C Get weather for current month
C
      RAIN=AVERAIN(THISMON)
	EVAPW=AVEPET(THISMON)
	AIRTEMP=AVETEMP(THISMON)
C
C Open met.file
C
      IF(ISOPEN(THISYEAR).NE.1)THEN
	  IF(THISYEAR.GT.1)THEN 
	    CLOSE(3)
	    ISOPEN(THISYEAR)=0
	  ENDIF
	  OPEN(3,FILE=METFILE(THISYEAR),STATUS='OLD',ERR=111)
	  GOTO 101
111     CONTINUE
        IF(THISMON.EQ.1)
     &    WRITE(*,10)THISYEAR,METFILE(THISYEAR)
10        FORMAT('Note LTA data used in year ',I4,
     &           ' Cannot open met file ',A20)
	  GOTO 102
101     CONTINUE
        ISOPEN(THISYEAR)=1
C
C Move to line before current timestep
C
        DO 100 IMON=1,THISMON-1,1
	    READ(3,*,ERR=112)
	    GOTO 103
112       CONTINUE
          PRINT*,'Note! Data missing in ',METFILE(THISYEAR),' LTA used'
	    GOTO 104
103       CONTINUE
100     CONTINUE
	ENDIF
C
C Read in this month
C
      READ(3,*,ERR=113)IMON,RAIN,EVAPW,AIRTEMP
      GOTO 105
113   CONTINUE
      PRINT*,'Note! Data missing in ',METFILE(THISYEAR),' LTA used ',
     &       'for month ',THISMON
	GOTO 104
105   CONTINUE
104   CONTINUE
C
C Error check read in data
C
      IF(RAIN.LT.0.OR.RAIN.GT.5000)THEN
	  PRINT*,'Error in rain: ',RAIN,'mm in month ',THISMON,
     &         ' LTA used'
	  RAIN=AVERAIN(THISMON)
	ENDIF
	IF(AIRTEMP.LT.-100.OR.AIRTEMP.GT.50)THEN
	  PRINT*,'Error in temp: ',AIRTEMP,'deg.C in month ',THISMON,
     &         ' LTA used'
	  AIRTEMP=AVETEMP(THISMON)
	ENDIF
      IF(EVAPW.LT.-998.AND.EVAPW.GT.-1000)
c     &  CALL GET_THIS_PET(GRIDLAT,AIRTEMP,EVAPW,THISMON)
     &  EVAPW=AVEPET(THISMON)
	IF(EVAPW.LT.0.OR.EVAPW.GT.500)THEN
	  PRINT*,'Error in evap: ',EVAPW,'mm in month ',THISMON,
     &         ' LTA used'
	  EVAPW=AVEPET(THISMON)
	ENDIF
C
C Set soil temperature from air temperature
C Assume soil temperature is same thoughout profile
C
102   CONTINUE
      DO 200 IL=1,MAXLAYER1
	  SOILTEMP(IL)=AIRTEMP
200   CONTINUE
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE SAVE_VSD_STATE(ISAVE,VSD_STATE_IL,IL)
C
C Subroutine to save the state of VSD in layer IL
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	REAL VSD_STATE(MAXLAYER,10)	! VSD state in layers
	INTEGER IS					! Item in state array
C
C Variables passed to/from this subroutine
C
	INTEGER IL					! IN: Layer to be saved/retrieved
	INTEGER ISAVE				! IN: Save=1; Retrieve=0
      REAL VSD_STATE_IL(10)		! IN/OUT: State of VSD in layer IL
C
C Save state within this routine
C
      SAVE VSD_STATE
C
C Save state of layer IL
C
      IF(ISAVE.EQ.1)THEN
	  DO 100 IS=1,10
	    VSD_STATE(IL,IS)=VSD_STATE_IL(IS)
100     CONTINUE
C
C Retrieve state of layer IL
C
      ELSEIF(ISAVE.EQ.0)THEN
	  DO 200 IS=1,10
	    VSD_STATE_IL(IS)=VSD_STATE(IL,IS)
200     CONTINUE
      ENDIF
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE SETFILE_LIM(
C OUTPUTS: entered mode of run, drainage class                        
     &                 ICMODEL,ISDYNPI,EQMODEL,EC_EQRUN,DRAINCLASS,
C          soil type, available water type, soil organic matter layers, 
C          land uses, number of growing seasons
     &                 NSOIL,IAWC,SOMDEPTH,LUSEQ,NXYEARS,
C          soil C, %clay, %silt, %sand, BD, pH, impermeable layer / depth
     &                 DOMSOILC,DOMSOILBD,DOMSOILPH,INMODE,
     &                 DOMSOILCLAY,DOMSOILSILT,DOMSOILSAND,
     &                 DOMSOILISIMP,DOMSOILIMPDEPTH,
C		 plant input, long term average rain, PET, temp, metfiles
     &				 TOTPIC,AVERAIN,AVEPET,AVETEMP,METFILE,WTABLE,
C          annual C accumulation, CH4 emission, CO2 emission, DOC emission, PI
     &                 CACCUM,ANCH4,ANCO2,ANDOC,GRIDLAT,NSOMLAY,PI_SPEC)
C
C Subroutine to get filenames and simulation inputs
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	REAL AVE_ANNTEMP            ! Long term average annual temperature [degC]
	REAL AVE_ANNRAIN            ! Long term average annual precipitation [mm/year]
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/
	INTEGER MAXGROW				! Max.no.of growing seasons allowed
	PARAMETER (MAXGROW=300)
	INTEGER MAXLUSOIL   		! Max.no.of land use types
      PARAMETER (MAXLUSOIL=6)
	INTEGER MAXLU				! Max.no.of land use types
      PARAMETER (MAXLU=10)
	INTEGER MAXLU1				! Max.no.of land use types
	DATA MAXLU1 /6/
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXSOMLAY			! Maximum number of SOM layers
	PARAMETER (MAXSOMLAY=10)			

	INTEGER ILU					! Counter for land use 
	INTEGER IMON				! Counter for months
      REAL INMODE				! Input mode of model run 
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
	INTEGER EC_EQRUN			! OUT:Initialisation using a full ECOSSE equilibrium run (on or off)
	INTEGER EQMODEL				! OUT:Type of equilibrium run 
	                            ! (NPP, TOC or both)
	INTEGER EQNPP				! Model initialised using plant inputs 
								!    from measured NPP
	INTEGER EQTOC				! Model initialised using measured TOC
	INTEGER EQNPPTOC			! Model initialised using both plant inputs
								!    from measured NPP and measured TOC
	INTEGER EQHILLIER,EQJONES	! Model initialised using TOC and
								!    the Hillier solver to get C pools and plant inputs
	DATA EQNPP,EQTOC,EQNPPTOC,EQHILLIER,EQJONES /1,2,3,4,5/
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
C
C ...Soil factors
C
      REAL ANCH4					! OUT: Measured annual CH4 emissions
								!     kgC/ha/yr
	REAL ANCO2					! OUT: Measured annual CO2 emissions
								!     kgC/ha/yr
	REAL ANDOC					! OUT: Measured annual DOC loss
								!     kgC/ha/yr
      REAL CACCUM					! OUT: Measured annual C accumulation 
								!     in the soil kgC/ha/yr
	REAL DOMSOILC(CELLSTOGRID*MAXSERIES/2,MAXLUSOIL,MAXSOMLAY)	! OUT:Soil C in 5 
												! major soil series under different LU 
												! in SOM layers (kgC/ha)
	REAL DOMSOILCLAY(CELLSTOGRID*MAXSERIES/2,MAXLUSOIL,MAXSOMLAY)	! OUT:% Soil clay in 5 
												! major soil series under different LU 
	REAL DOMSOILIMPDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLUSOIL)	! IN: Depth of impermeable layer if present (cm)
	INTEGER DOMSOILISIMP(CELLSTOGRID*MAXSERIES/2,MAXLUSOIL)	! IN: Flag for impermeable layer (0=No; 1=Yes)
	REAL DOMSOILSILT(CELLSTOGRID*MAXSERIES/2,MAXLUSOIL,MAXSOMLAY)	! OUT:% Soil silt in 5 
												! major soil series under different LU 
	REAL DOMSOILSAND(CELLSTOGRID*MAXSERIES/2,MAXLUSOIL,MAXSOMLAY)	! OUT:% Soil sand in 5 
												! major soil series under different LU 
	REAL DOMSOILBD(CELLSTOGRID*MAXSERIES/2,MAXLUSOIL,MAXSOMLAY)	! OUT:Soil BD in 5 
												! major soil series under different LU 
												! in SOM layers (g/cm3)
	REAL DOMSOILPH(CELLSTOGRID*MAXSERIES/2,MAXLUSOIL,MAXSOMLAY)	! OUT:Soil PH in 5 
												! major soil series under different LU 
												! in SOM layers (g/cm3)
      INTEGER DRAINCLASS			! OUT:Drainage class
	INTEGER IAWC				! OUT:Available water code - set to 1
      INTEGER ISERROR				! IN:Code for error in soil 1=error 0=noerror
      INTEGER ISMISS				! IN:Code for missing values in soil 1=missing 0=nomissing
      INTEGER LUSEQ(MAXGROW)		! OUT:Land use sequence
      INTEGER NSOIL				! OUT:Soil code - set to 1
	INTEGER NSOMLAY(CELLSTOGRID*MAXSERIES/2,MAXLUSOIL)	! OUT:Number of SOM layers
	INTEGER NXYEARS				! OUT:No.growing seasons simulated
	REAL PI_SPEC(MAXGROW)		! OUT: Total annual plant C input (kg C / ha / yr) 
      REAL SOMDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLUSOIL,MAXSOMLAY)	! OUT:Depth of each SOM layer
	REAL WTABLE,ANNTEMP,ANNRAIN					! OUT:Water table depth in cm
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
      CHARACTER*40 METFILE(MAXGROW)	! OUT:Met.file  
C
C Soil type
C
      ISDYNPI=ISDYNPI_OFF
      NSOIL=1
	IAWC=1
C
C Get input file name
C
111   CONTINUE
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
	PRINT*,'          8=adjusting PI to give measured TOC at ',
     & 'equilibrium, and further adjusting PI according to weather. ',
     & 'No ECOSSE equilibration' 
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
      READ(*,10)INFILE
	PRINT*,'======================================================'
10    FORMAT(A50)
C
C Open data file
C
      OPEN(42,FILE=INFILE,STATUS='OLD',ERR=111)
C
C Read in data
C
C INMODE: 1=use plant inputs; 2=use toc; 3=use plant inputs and toc; 4=use C accumulation
      READ(42,*)INMODE
	EC_EQRUN=0     
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
	ELSEIF(INMODE.GE.6.AND.INMODE.LT.7)THEN
	  EQMODEL=EQHILLIER
	  ICMODEL=ICROTHCEQ	
	ELSEIF(INMODE.EQ.7)THEN
	  EQMODEL=EQTOC
	  ICMODEL=ICROTHCEQ
	  ISDYNPI=ISDYNPI_ON
	ELSEIF(INMODE.EQ.8)THEN
	  EQMODEL=EQTOC
	  ICMODEL=ICROTHCEQ
	  ISDYNPI=ISDYNPI_ON
	ELSEIF(INMODE.GE.9.OR.INMODE.LT.10)THEN
	  EQMODEL=EQJONES
	  ICMODEL=ICROTHCEQ
	  EC_EQRUN=0
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
	READ(42,*)WTABLE
	READ(42,*)DRAINCLASS
	READ(42,*)CACCUM
	READ(42,*)ANCH4
	READ(42,*)ANCO2
	READ(42,*)ANDOC
	READ(42,*)NXYEARS
C
C Estimate each land use's long term average plant inputs using MIAMI_DYCE
C if they have not been specified in the input file
C
      AVE_ANNRAIN = 0
      AVE_ANNTEMP = 0
      DO IMON=1,12
          AVE_ANNRAIN = AVE_ANNRAIN + AVERAIN(IMON)
          AVE_ANNTEMP = AVE_ANNTEMP + AVETEMP(IMON)
      ENDDO
      AVE_ANNTEMP = AVE_ANNTEMP / 12
      DO ILU=1,MAXLU
          IF(TOTPIC(ILU).EQ.0)THEN    
              CALL MIAMI_DYCE(ILU,AVE_ANNTEMP,AVE_ANNRAIN,TOTPIC(ILU))
              PRINT *,'Plant inputs brought to you by Miami Dyce'
			  PRINT*, ILU, TOTPIC(ILU)
          ENDIF
      ENDDO
C
C Read in LU data and check profile information for each LU included
C
      DO 650 ILU=1,MAXLU1
	  ISCHECKED(ILU)=0
650   CONTINUE
	DO 700 IYEAR=1,NXYEARS
	  PI_SPEC(IYEAR)=0
	  READ(42,*)LUSEQ(IYEAR),PI_SPEC(IYEAR)
        IF(ISCHECKED(LUSEQ(IYEAR)).EQ.0)THEN
          if (LUSEQ(IYEAR) < 7) then
              CALL CHECK_PROFILE(LUSEQ(IYEAR),NSOMLAY,
     &                       SOMDEPTH,DOMSOILC,DOMSOILCLAY,
     &                       DOMSOILIMPDEPTH,DOMSOILISIMP,DOMSOILSILT,
     &                       DOMSOILSAND,DOMSOILBD,DOMSOILPH)
          endif          
	    ISCHECKED(LUSEQ(IYEAR))=1
	  ENDIF
700   CONTINUE
	DO 800 IYEAR=1,NXYEARS
        READ(42,*,END=801)METFILE(IYEAR)
800   CONTINUE
801   CONTINUE
C
C If NPP entered as zero, calculate NPP values using the MIAMI model
C
      DO 900 ILU=1,MAXLU
	 IF(TOTPIC(ILU).LT.0.0001.AND.TOTPIC(ILU).GT.-0.0001)THEN
	 ANNRAIN=0
	 ANNTEMP=0
	 DO 1100 IMON=1,12
        ANNRAIN=ANNRAIN+AVERAIN(IMON)
	  ANNTEMP=ANNTEMP+AVETEMP(IMON)
1100   CONTINUE
      ANNTEMP=ANNTEMP/12
	IF(INMODE.LT.6.OR.INMODE.GT.7)THEN
	 CALL MIAMI_DYCE(ILU,ANNTEMP,ANNRAIN,TOTPIC(ILU))
!      CALL MIAMI(TOTPIC(ILU),AVERAIN,AVETEMP)
	  ENDIF
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
C Leave SETFILE_LIM
C
      RETURN
      END