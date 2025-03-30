C------------------------------------------------------------
C
C Site routines for ECOSSE Carbon and Nitrogen Turnover Model
C using on limited data to run the model
C
C Jo Smith 
C Started 29/01/08
C
C
C
C*************************************************************
C MAIN RUN ROUTINES
C*************************************************************
C
C EXTERNAL SUBROUTINES
C TEST_ECOSSE_FOR_VSD()
C INIT_ECOSSE_FOR_VSD()
C RUN_ECOSSE_FOR_VSD()
C
C-------------------------------------------------------------
C
      SUBROUTINE TEST_ECOSSE_FOR_VSD(ISWAIT)
C
C Subroutine to initialise ECOSSE for VSD
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
	INTEGER MAXFERT				! Max.no.of fertiliser applications allowed 
	PARAMETER (MAXFERT=5)
	INTEGER MAXGROW				! Max.no.of growing seasons allowed
	PARAMETER (MAXGROW=300)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXLU1				! Max.no.of land use types
	DATA MAXLU1 /6/
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXSOMLAY			! Max.no.of soil layers entered			
	PARAMETER (MAXSOMLAY=10)

	INTEGER IL					! Local counter variable
C
C Variables passed to / from this subroutine
C ...Model factors
C
	INTEGER ICMODEL				! IN/OUT:Type of C initialisation model 
	INTEGER EC_EQRUN			! IN(SETFILE_LIM):Initialisation using a full ECOSSE equilibrium run (on or off)
	INTEGER EQMODEL				! IN/OUT:Type of equilibrium run 
	                            ! (NPP, TOC or both)     
	INTEGER ISDYNPI				! Code for dynamic PI (adjusted by external factors)
	INTEGER ISDYNPI_OFF			! Adjustment of PI by external factors is off
	INTEGER ISDYNPI_ON			! Adjustment of PI by external factors is on
	DATA ISDYNPI_OFF, ISDYNPI_ON /0,1/
      INTEGER ISERROR				! IN:Code for error in soil 1=error 0=noerror
	INTEGER ISWAIT				! IN/OUT:Code to wait for key press (1) or not (0)
C
C ...Timing factors
C
	INTEGER IK					! OUT:No.timesteps from prev.harv. to current
	INTEGER IYEAR				! OUT:Current growing season number
	INTEGER MEND				! IN/OUT:Number of timesteps from prev.harv.-harv.
	INTEGER SUM_TS				! IN:Total number of timesteps passed
C
C ...Soil factors
C
      REAL ANCH4					! OUT: Measured annual CH4 emissions
								!     kgC/ha/yr
	REAL ANCO2					! OUT: Measured annual CO2 emissions
								!     kgC/ha/yr
	REAL ANDOC					! OUT: Measured annual DOC loss
								!     kgC/ha/yr
	REAL BULKDENS(MAXLAYER)		! IN:Bulk density of the layer (g/cm3)
      REAL CACCUM					! OUT: Measured annual C accumulation 
								!     in the soil kgC/ha/yr
      REAL CLAY(MAXLAYER)			! IN:Clay content of the layer (%)
	REAL CO2(MAXLAYER)			! IN:CO2 emitted (kgC/ha/layer)
	REAL CMOBDOC(MAXLAYER)		! OUT: Concentration of mobile DOC in this layer (eq/m3)
	REAL CNH4(MAXLAYER)			! OUT: Concentration of ammonium in this layer (eq/m3)
      REAL CNO3(MAXLAYER)			! OUT: Concentration of nitrate in this layer (eq/m3)
	REAL DOMSOILC(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN:Soil C in 5 
															! major soil series under different LU 
															! in SOM layers (kgC/ha)
	REAL DOMSOILCLAY(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN:% Soil clay in 5 
															! major soil series under different LU 
	REAL DOMSOILIMPDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN: Depth of impermeable layer if present (cm)
	INTEGER DOMSOILISIMP(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN: Flag for impermeable layer (0=No; 1=Yes)
	REAL DOMSOILSILT(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN:% Soil silt in 5 
															! major soil series under different LU 
	REAL DOMSOILSAND(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN:% Soil sand in 5 
															! major soil series under different LU 
	REAL DOMSOILBD(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN:Soil BD in 5 
															! major soil series under different LU 
															! in SOM layers (g/cm3)
	REAL DOMSOILPH(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! OUT:Soil PH in 5 
												! major soil series under different LU 
												! in SOM layers (g/cm3)
      INTEGER DRAINCLASS			! IN:Drainage class
	REAL DRAINW(MAXLAYER)		! IN:Water drained from this layer (mm/layer)
	INTEGER IAWC				! IN:Water movement code number
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
	INTEGER NSOIL				! IN:Soil code number
	INTEGER NSOMLAY(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN:Number of SOM layers
	REAL SAND(MAXLAYER)			! IN:Sand content of the layer (%)
	REAL SILT(MAXLAYER)			! IN:Silt content of the layer (%)
	REAL SOILPH(MAXLAYER)		! IN/OUT:pH of soil in this layer
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL TOC(MAXLAYER)			! IN:Measured total organic C (kgC/ha/layer)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
	REAL WTABLE	,INMODE				! IN:Water table depth in cm
C
C ...Plant factors
C
      INTEGER LUSEQ(MAXGROW)		! Land use sequence
	INTEGER NXYEARS			! IN:No.growing seasons simulated
	REAL PI_SPEC(MAXGROW)		! IN(SETFILE_LIM): Total annual plant C input (kg C / ha / yr) 
	INTEGER THISLU				! OUT:LU in the current growing season
	REAL TOTPIC(MAXLU)			! IN/OUT:Plant C input calculated using 
								!	  ROTHC equilibrium run for each LU (kgC/ha/yr)
C
C ...Weather factors
C
      REAL AVEPET(12)				! IN:Long term average PET (mm/month)
      REAL AVERAIN(12)			! IN:Long term average rainfall(mm/month)
      REAL AVETEMP(12)			! IN:Long term average monthly average 
      REAL GRIDLAT				! IN/OUT:Average latitude of this 20km2 grid cell
C
C ... GIS specific variables
C
	INTEGER ILU					! IN:Counter for land use 
	INTEGER ISERIES				! IN:Counter for soil series
	CHARACTER*40 METFILE(MAXGROW)	! IN:Met.file  

	REAL SOMDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN:depth of soil organic matter layers
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
C
C Get soil characteristics 
C
      ISERIES=1
	CALL GETSOIL(ISERIES,LUSEQ(1),SOMDEPTH,NSOMLAY,							
     &             DOMSOILC,DOMSOILCLAY,DOMSOILSILT,				
     &             DOMSOILSAND,DOMSOILBD,DOMSOILPH,
     &             DOMSOILISIMP,DOMSOILIMPDEPTH,
     &             TOC,IOM,CLAY,SILT,SAND,BULKDENS,SOILPH,ISERROR)
	IF(ISERROR.EQ.1)THEN
	  PRINT*,'ERROR IN SOIL PARAMETERS'
        GOTO 606
	ENDIF
C
C Initialise ECOSSE for VSD
C
C***********************************************************************
C INIT_ECOSSE_FOR_VSD = Subroutine to initialise ECOSSE for VSD
C Argument list:	ICMODEL		=	Type of C initialisation model 
C 								(1=Initialisation of C is fixed;
C	                             2=Initialisation of C by RothC)
C				EQMODEL		=	Type of equilibrium run 
C								(1=Model initialised using plant inputs (measured NPP)
C								 2=Model initialised using measured TOC
C								 3=Model initialised using both plant inputs and measured TOC)
C                 ISERIES		=	Soil series number
C				LUSEQ(1)	=	Land use in year 1
C				TOC()		=	Measured total organic C (kgC/ha/layer)
C				CLAY()		=	Clay content of the layer (%)
C				SILT()		=	Silt content of the layer (%)
C				SAND()		=	Sand content of the layer (%)
C				BULKDENS()	=	Bulk density of the layer (g/cm3)
C				AVEPET()	=	Long term average PET (mm/month)
C				AVERAIN()	=	Long term average rainfall(mm/month)
C				AVETEMP()	=	Long term average temperature (deg.C/month)
C				TOTPIC()	=	Plant C input (kgC/ha/yr)
C				WMAX()		=	Available water at field cap. (mm/layer)
C				WSAT()		=	Available water at saturation (mm/layer)
C				SOILW()		=	Available water (mm/layer)
C                 SOILPH()	=   pH of soil in the layer
C
      CALL INIT_ECOSSE_FOR_VSD(ICMODEL,EQMODEL,ISERIES,LUSEQ(1),
     &                         TOC,CLAY,SILT,SAND,BULKDENS,
     &                         AVERAIN,AVEPET,AVETEMP,TOTPIC,
     &                         WMAX,WSAT,SOILW,SOILPH,WTABLE)
C***********************************************************************
C
C For each growing season... 
C        
      DO 600 IYEAR=1,NXYEARS
C
C If annual plant inputs are specified, and mode is not 7 (adjust PI due to weather conditions), set PI to specified value
C
        IF(PI_SPEC(IYEAR).GT.0)TOTPIC(LUSEQ(IYEAR))=PI_SPEC(IYEAR)
C
C For each month till end of growing season...
C           
        IK=1
700     CONTINUE
C
C Run ECOSSE for this month
C
C***********************************************************************
C RUN_ECOSSE_FOR_VSD = Subroutine to run ECOSSE for this month
C Argument list:	ICMODEL		=	Type of C initialisation model 
C 								(1=Initialisation of C is fixed;
C	                             2=Initialisation of C by RothC)
C				EQMODEL		=	Type of equilibrium run 
C								(1=Model initialised using plant inputs (measured NPP)
C								 2=Model initialised using measured TOC
C								 3=Model initialised using both plant inputs and measured TOC)
C				ISERIES		=	Soil series number
C				IYEAR		=	Current growing season number
C				IK			=	No.timesteps from prev.harv.to current month
C				MEND		=	Number of timesteps from prev.harv.to next harv.
C				LUSEQ()		=	Land use sequence
C				TOC()		=	Measured total organic C (kgC/ha/layer)
C				CLAY()		=	Clay content of the layer (%)
C				BULKDENS()	=	Bulk density of the layer (g/cm3)
C				AVERAIN()	=	Long term average rainfall(mm/month)
C				AVEPET()	=	Long term average PET (mm/month)
C				AVETEMP()	=	Long term average monthly average 
C				TOTPIC()	=	Plant C input (kgC/ha/yr)
C				WMAX()		=	Available water at field cap. (mm/layer)
C				WSAT()		=	Available water at saturation (mm/layer)
C				SOILW()		=	Available water (mm/layer)
C				CNO3()		=	Concentration of nitrate in this layer (eq/m3)
C				CNH4()		=	Concentration of ammonium in this layer (eq/m3)
C				CMOBDOC()	=	Concentration of mobile DOC in this layer (eq/m3)
C				CO2()		=	CO2 emitted (kgC/ha/layer)
C	            SOILPH()    =   pH of soil in this layer
C                 WTABLE		=   Water table depth (cm)
C                 ISIMP		=	Impermeable layer? 0=No; 1=Yes
C                 IMPDEPTH	=	Depth of impermeable layer (cm)
C
          CALL RUN_ECOSSE_FOR_VSD(ICMODEL,EQMODEL,
     &                            ISERIES,IYEAR,IK,MEND,LUSEQ,
     &							TOC,CLAY,BULKDENS,
     &                            AVERAIN,AVEPET,AVETEMP,TOTPIC,
     &                            WMAX,WSAT,SOILW,DRAINW,
     &                            CNO3,CNH4,CMOBDOC,CO2,SOILPH,WTABLE,
     &                            DOMSOILISIMP(ISERIES,LUSEQ(1)),
     &                            DOMSOILIMPDEPTH(ISERIES,LUSEQ(1)),
     &                            ISDYNPI)
C***********************************************************************

	    IK=IK+1
        IF(IK.LE.MEND)GOTO 700
C
C Go back and set up the next crop
C
600    CONTINUE
606   CONTINUE
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
      WRITE(*,*)'                       ....Press any key to continue'
	READ(*,*)
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE INIT_ECOSSE_FOR_VSD(ICMODEL,EQMODEL,ISERIES,LU1,
     &             TOC,CLAY,SILT,SAND,BULKDENS,
     &             AVERAIN,AVEPET,AVETEMP,TOTPIC,
     &             WMAX,WSAT,SOILW,SOILPH,WTABLE)
C
C Subroutine to initialise ECOSSE for VSD
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
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
	INTEGER MAXLU1				! Max.no.of land use types
	DATA MAXLU1 /6/
	INTEGER MAXORGM				! Max.no.of organic manure applications allowed
	PARAMETER (MAXORGM=52)
	INTEGER MAXSOIL				! Max.no.of soil types
	PARAMETER (MAXSOIL=50)

	INTEGER IL					! Local counter variable
	INTEGER ISTEST				! Code for test results (1) or not (0)
C
C Variables passed to / from this subroutine
C
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
	INTEGER ICMODEL				! OUT:Type of C initialisation model 
	INTEGER ICFIXED				! Initialisation of C is fixed
	INTEGER ICROTHCEQ			! Initialisation of C by RothC 
								!  equilibrium run
	DATA ICFIXED,ICROTHCEQ /1,2/
	INTEGER ISPINUP				! Is N limitation spin-up used?
	INTEGER ISPINUP_OFF			! N limitation spin-up is not used
	INTEGER ISPINUP_ON			! N limitation spin-up is used
	DATA ISPINUP_OFF,ISPINUP_ON /0,1/
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
	INTEGER CH4MODEL			! OUT:Methane model (on, Richards or Aitkenhead) 
	INTEGER CH4_OFF             ! CH4 model off
	INTEGER CH4_RICHARDS    	! Richards CH4 model on
	INTEGER CH4_AITKENHEAD   	! Aitkenhead CH4 model on		
      DATA CH4_OFF,CH4_RICHARDS,CH4_AITKENHEAD /0,1,2/
	INTEGER PH_MODEL			! How is pH calculated?
	INTEGER PH_PARAM			! pH set to neutral or passed from parameter file
      INTEGER PH_STATIC			! pH is read in from input
								!   file and does not change
	INTEGER PH_FROM_VSD			! pH is calculated using VSD (a very 
								!   simple dynamic version of the MAGIC model
								!   by Ed Rowe & Chris Evans, CEH, Bangor)
	DATA PH_PARAM,PH_STATIC,PH_FROM_VSD /0,1,2/			 
C
C ...Timing factors
C
	INTEGER FIXEND				! IN:Fixed end? 0=No 1=Yes
	INTEGER IANTHES				! IN:No.timesteps from 01/01 to anthesis 
	INTEGER IEND				! IN:Timesteps from sowing date to harvest
	INTEGER IK					! OUT:No.timesteps from prev.harv. to current
	INTEGER IRYEAR				! IN:Current weather year number
	INTEGER IYEAR				! IN:Current growing season number
	INTEGER JSTOP				! IN:Code to stop simulation 0=No 1=Yes
	INTEGER LHARV				! IN:No.timesteps from 01/01 to prev. harvest
	INTEGER MEND				! IN:Number of timesteps from prev.harv.-harv.
	INTEGER N_STEPS				! IN:No.timesteps in a year
	REAL SECONDS				! IN:Number of seconds in one timestep
	INTEGER SUM_TS				! IN:Total number of timesteps passed
C
C ...Soil factors
C
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
	REAL BCARB0(MAXLAYER)		! IN:C in soil biomass (kgC/ha/layer)
      REAL BNIT0(MAXLAYER)		! IN:N in soil biomass (kgN/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN:N15 in soil humus (kgN15/ha/layer)
	REAL BSTART					! IN:Initial N in biomass pool (kgN/ha)
	REAL BSTART15				! IN:Initial N in biomass pool (kgN15/ha)
	REAL BULKDENS(MAXLAYER)		! IN:Bulk density of the layer (g/cm3)
      REAL CACCUM					! IN/OUT: Measured annual C accumulation 
								!     in the soil kgC/ha/yr
	REAL C_TS					! IN:Litter C input in this timestep (kgC/ha)
      REAL CLAY(MAXLAYER)			! IN:Clay content of the layer (%)
	REAL CO2(MAXLAYER)			! IN:CO2 emitted (kgC/ha/layer)
	REAL CRIT(MAXSOIL,MAXLAYER)	! IN:Critical N content (kgN/ha/layer)
	REAL CTOT					! IN:Total litter C input (kgC/ha)
      REAL DPM_RPM_RATIO			! IN/OUT:Ratio of DPM:RPM. If set to 
								!        zero, this will be worked out using a call to SET_DPMRPMRATIO
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
      REAL DPMCTON(MAXLAYER)		! IN/OUT:DPM C:N ratio
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL DRAINW(MAXLAYER)		! IN:Water drained from this layer (mm/layer)
	REAL FIELDCAP(MAXLAYER)	    ! IN:Soil water content at field capacity (mm/layer)
	REAL FLOWPROP				! IN:Proportion of flow needed to achieve
							    !	     water table at depth WTABLE
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN:N in soil humus (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN:N15 in soil biomass (kgN15/ha/layer)
	REAL HZ1(MAXLAYER)			! IN:N:C of input PM for steady state 
	INTEGER IAWC				! IN:Water movement code number
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
	INTEGER IROCK				! IN:No.layers: 1=50cm, 2=100cm , 3=150cm
	INTEGER ISERIES				! IN:Counter for soil series
	INTEGER NSOIL				! IN:Soil code number
      REAL PI_C(MAXLAYER)			! IN/OUT: Plant input C to soil (kgC/ha/step)
	REAL PI_CEQ_MON(12,MAXLAYER)! IN/OUT: Equilibrium plant C input each month
								!         in each layer (kgC/ha/month/layer)
      REAL PI_N(MAXLAYER)			! IN/OUT: Plant input N to soil (kgC/ha/step)
      REAL PI_N15(MAXLAYER)		! IN/OUT: Plant input N15 to soil(kgC/ha/step)
								!(1=arable,2=grassland,3=forestry,4=semi-natural)
	REAL REFIL(MAXLAYER)		! IN:Water needed to fill to FC (mm/lay)
	REAL RNIN					! IN:Litter N input in this timestep (kgN/ha)
	REAL RNIN15					! IN:Litter N15 input in timestep (kgN15/ha)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
      REAL RPMCTON(MAXLAYER)		! IN/OUT:RPM C:N ratio
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL RSTART					! IN:Initial N in debris pool (kgN/ha)
	REAL RSTART15				! IN:Initial N in debris pool (kgN15/ha)
	REAL SAND(MAXLAYER)			! IN:Sand content of the layer (%)
	REAL SATWATCONT(MAXLAYER)	! IN:Total water content at saturation (mm/layer)
	REAL SILT(MAXLAYER)			! IN:Silt content of the layer (%)
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOILPH(MAXLAYER)		! IN/OUT:pH of soil in this layer
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL SORGN15				! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL T15					! IN:Net Mineralised N (kgN/ha/layer)
	REAL TAM(MAXLAYER)			! IN:15N/N in the ammonium in the layer
	REAL TFYMC					! IN:Total FYM C input (kgC/ha) DONT PASS?
	REAL TIN(MAXLAYER)			! IN:15N/N in the nitrate in the layer
	REAL TOC(MAXLAYER)			! IN:Measured total organic C (kgC/ha/layer)
	REAL TOTAST					! IN:Total ammonium in soil (kgN/ha)
	REAL TOTNST					! IN:Total nitrate in soil (kgN/ha)
	REAL TOTAST15				! IN:Total ammonium in soil (kgN15/ha)
	REAL TOTNST15				! IN:Total nitrate in soil (kgN15/ha)
	REAL VOLAT					! IN:N lost by volatilisation (kgN/ha)
	REAL VOLAT15				! IN:N15 lost by volatilisation (kgN15/ha)
	REAL WILTPOINT(MAXLAYER)	! IN:Water conent at wilting point (mm/layer)
	REAL WLEACH					! IN:N leached (kgN/ha)
	REAL WLEACH15				! IN:N15 leached (kgN15/ha)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
	REAL WTABLE					! IN:Water table depth in cm
C
C ... Crop factors
C
	REAL CACT					! IN:Crop N uptake (kgN/ha/timestep)
	REAL CACT15					! IN:Crop N15 uptake (kgN15/ha/timestep)
	REAL CACTOT					! IN:N taken up by crop (kgN/ha)
	REAL CATOT15				! IN:N15 taken up by crop (kgN15/ha)
	REAL CLOSSX					! IN:N lost by senescence (kgN/ha) 
	REAL CLOSSX15				! IN:N15 lost by senescence (kgN15/ha) 	
	REAL FIXN					! IN:N fixed from atmosphere (kgN/ha)
	INTEGER ISOWN				! IN:Timesteps from 01/01/01 to sowing date 
	INTEGER IHARV				! IN:Timesteps from 01/01/01 to harvest date 
      INTEGER LUSEQ(MAXGROW)		! Land use sequence
	INTEGER LU1					! OUT:First land use
	INTEGER NFERT				! IN:No.fertiliser applications to this crop
	INTEGER NORGM				! IN:No.manure applications to this crop
	INTEGER NSOW				! IN:Timesteps from prev.harv. to sowing date
	INTEGER NXYEARS				! IN:No.growing seasons simulated
	REAL ROOT					! IN:Rooting depth according to restriction (cm)
	INTEGER THISLU				! OUT:LU in the current growing season
	REAL TOTPIC(MAXLU)			! IN/OUT:Plant C input calculated using 
								!	  ROTHC equilibrium run for each LU (kgC/ha/yr)
	REAL WR						! IN:Root requirement (kgN/ha)
C
C ... Fertiliser factors
C
	REAL FERT(MAXFERT)			! IN:Amount of fertiliser applied (kgN/ha)
	REAL FERTNF					! INOUT:Amount of fertiliser applied (kgN/ha)
	REAL FYMFERT				! IN:N input by fertiliser & FYM (kgN/ha)
	REAL FYMFERT15				! IN:N15 input by fertiliser & FYM (kgN15/ha)
	INTEGER IFERT(MAXFERT)		! IN:No.timesteps to split fert.application
	INTEGER ILAB(MAXFERT)		! IN:Labelling on fertiliser? 0=No 1=Yes
	INTEGER ILABNF				! INOUT:Labelling on fertiliser? 0=No 1=Yes
	INTEGER IVOL(MAXFERT)		! IN:N volatilised from fertiliser (kgN/ha)
	INTEGER IVOLNF				! INOUT:N volatilised from fertiliser (kgN/ha)
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
	REAL ORGMA(MAXORGM)			! IN:Amount of manure applied (t/ha)
      REAL ORGMJ(0:MAXGROW,MAXORGM)	! IN:Amount of manure applied (t/ha)
	INTEGER IOLABNF				! IN:Labelling on manure? 0=No 1=Yes
	INTEGER JORGMNF				! IN:Type of manure application
	REAL ORGMANF				! IN:Amount of manure applied (t/ha)
C
C ...Weather factors
C
	REAL LTA_AWC(12,MAXLAYER)	! IN:Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)	! IN:Average air temp this month (deg.C)
	REAL LAT					! IN: Latitude
	REAL RAIN					! IN: Rainfall (mm/timestep)
	REAL EVAP					! IN: Potential evap. (mm/timestep)
	REAL EVAPW					! IN: Potential evap. (mm/timestep)
								!	CONSIDER COMBINING EVAP AND EVAPW
      INTEGER IDEC				! OUT:Decade counter 
	REAL SOILTEMP(MAXLAYER)		! IN: Soil temperature (deg.C/timestep)
	REAL SEEDIN					! IN: Data used by MAGEC
C
C ...Dissolved organic matter factors
C
	REAL CO2FROMDOC(MAXLAYER)	! IN:CO2 output from DOC in each layer(kgC/ha)
	REAL MOBDOC(MAXLAYER)		! IN:Mobile DOC in each layer (kgC/ha)
	REAL MOBDON(MAXLAYER)		! IN:Mobile DON in each layer (kgN/ha)
C
C ... GIS specific variables
C
      REAL AVEPET(12)				! IN:Long term average PET (mm/month)
      REAL AVERAIN(12)			! IN:Long term average rainfall(mm/month)
      REAL AVETEMP(12)			! IN:Long term average monthly average 
  	REAL ICFACTOR(MAXLAYER)		! IN:Adjustment in rate needed to	
								!    achieve measured NPP and TOC	
	INTEGER ILU					! IN:Counter for land use 
      INTEGER ISAVE				! OUT:Code to save or retrieve variables
								!		0=save; 1=retrieve
	REAL PNREQ(MAXLU)			! IN/OUT:Amount of N required by plant each yr
	                            !        (kgN/ha/year)
C
C Required for Monte Carlo - 
C
     	REAL DPMCARB(MAXLAYER)
	REAL RPMCARB(MAXLAYER)
	REAL BCARB(MAXLAYER)
	REAL HCARB(MAXLAYER)
	REAL WRATEDM
	REAL TRATEM
C
C Test results
C Set up results arrays
C
      ISTEST=1
      IF(ISTEST.EQ.1)THEN
        CALL STARTRES(IRYEAR,IK,RSTART,DPMNIT0,RPMNIT0,
     &                  BSTART,BNIT0,HNIT0,
     &                  RSTART15,DPMNLAB0,RPMNLAB0,
     &                  BSTART15,BNLAB0,HNLAB0,
     &                  TOTAST,AMMN,TOTAST15,AMMN15,
     &                  TOTNST,SOILN,TOTNST15,SOIL15)
	ENDIF

C
C Set restrictions on calculation
C Get parameters from TOC,IOM,CLAY,SILT,SAND,BULKDENS,	
C no restriction due to drainage,
C initialisation using ROTHC model & Bente Foereid's C:N ratios
C and get pH / water restriction on decomposition from NPP and TOC
C Use only 1 decade calculation
C
      MODTYPE=IS_SUNDIAL
	SPARMODEL=SPARCALC
C      ISPINUP=ISPINUP_ON
      ISPINUP=ISPINUP_OFF
	INMODEL=INPASSCN
	ITFUNC=ITFUNC_ROTHC
	IMFUNC=IMFUNC_ROTHC
	PH_MODEL=PH_FROM_VSD
	NORGM=0
C
C Set IOM from TOC
C
      DO 100 IL=1,MAXLAYER1
        CALL GET_IOM_FROM_FALLOON_EQN(TOC(IL),IOM(IL))
100   CONTINUE
C
C ...Set time factors (timestep=monthly)
C
	C_TS=0.0
	SECONDS=(365.25/12)*24*60*60
C
C ...Initialise water
C              
      NSOIL=1
	IAWC=1
	IROCK=3
      
      CALL GET_PLANT_DIST(TOTPIC(LU1),PI_CEQ_MON,LU1)    !enables icover to be worked out properly   
      
	CALL INIT_GIS_WATER(BULKDENS,CLAY,SAND,FLOWPROP,IAWC,IROCK,
     &                        NSOIL,SILT,TOC,WMAX,WSAT,WTABLE,ROOT,
     &                        AVEPET,AVERAIN,AVETEMP,SOILW,PI_CEQ_MON,
     &                        WILTPOINT,FIELDCAP,SATWATCONT,SECONDS)
C
C ...Get weather data
C
	CALL GETLTA(LTA_AWC,LTA_TEMP,WMAX,WSAT,FLOWPROP,
     &                AVERAIN,AVEPET,AVETEMP)
C
C ...Get soil C and N parameters
C
	CALL GET_SUNDIAL_SOILPARS(NSOIL,TOC,IOM,CRIT,LU1,
     &                          CLAY,BULKDENS,SPARMODEL)
C
C ...Get plant distribution for this land use
C
     	CALL GET_PLANT_DIST(TOTPIC(LU1),PI_CEQ_MON,LU1)
C
C ...Initialise the soil for this land use
C
      IF(INMODEL.EQ.INPASSCN)
     &  CALL GET_CTON_FOEREID(DPMCTON,RPMCTON,LU1)
	IF(ICMODEL.EQ.ICROTHCEQ)
     &  CALL SET_DPMRPMRATIO(LU1,DPM_RPM_RATIO)
	CALL INIT_GIS_SOILCN_NOPARS(ICMODEL,INMODEL,DOCMODEL,
     &                                EQMODEL,SECONDS,
     &                                NSOIL,TOC,IOM,HZ1,CRIT,
     &                                SOILN,SOIL15,AMMN,AMMN15,
     &                                DPMCARB0,DPMNIT0,DPMNLAB0,
     &                                RPMCARB0,RPMNIT0,RPMNLAB0,
     &                                BCARB0,BNIT0,BNLAB0,
     &                                HCARB0,HNIT0,HNLAB0,
     &                                 MOBDOC,MOBDON,SATWATCONT,
     &                                CLAY,BULKDENS,PI_CEQ_MON,
     &                                LTA_AWC,LTA_TEMP,ITFUNC,IMFUNC,
     &                                WMAX,WSAT,ICFACTOR,wiltpoint,
     &                                DPM_RPM_RATIO,DPMCTON,RPMCTON,
     &                                PH_MODEL,SOILPH,
     &                                CACCUM,ANCH4,ANCO2,ANDOC)
	DO 200 IL=1,MAXLAYER1
	  IF(BCARB0(IL).LT.0.OR.HCARB0(IL).LT.0.OR.
     &     RPMCARB0(IL).LT.0.OR.DPMCARB0(IL).LT.0.OR.
     &     BNIT0(IL).LT.0.OR.HNIT0(IL).LT.0.OR.
     &     RPMNIT0(IL).LT.0.OR.DPMNIT0(IL).LT.0.OR.
     &     SOILN(IL).LT.0.OR.AMMN(IL).LT.0)THEN
	     PRINT*,'ERROR! Negative state variables found!'
	     PRINT*,'...Press any key to continue'
	     READ(*,*)
	     STOP
	  ENDIF
200   CONTINUE
C
C ...Get rate modifier associated with SUNDIAL routines,
C
      IF(ISPINUP.EQ.ISPINUP_ON)THEN
        CALL GET_NLIM(AMMN,AMMN15,ATM,AVEPET,AVERAIN,AVETEMP,
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
     &              WLEACH15,WSAT,WMAX,SOILPH,PI_CEQ_MON)		  
	ENDIF
C
C Save results for VSD run
C
      ISAVE=1
      CALL SAVE_FOR_VSD(ISAVE,
     &                  SOILN,SOIL15,AMMN,AMMN15,
     &                  DPMCARB0,DPMNIT0,DPMNLAB0,
     &                  RPMCARB0,RPMNIT0,RPMNLAB0,
     &                  BCARB0,BNIT0,BNLAB0,
     &                  HCARB0,HNIT0,HNLAB0,
     &                  NFERT,IFERT,ILAB,FERT,TFERT,IVOL,
     &                  C_TS,FLOWPROP,HZ1,RNIN,RNIN15,
     &                  PI_C,PI_CEQ_MON,PI_N,PI_N15,PNREQ,
     &	              CACT,CACT15,CACTOT,CATOT15,NSOW,ICFACTOR)
C
C Leave INIT_ECOSSE_FOR_VSD
C
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE RUN_ECOSSE_FOR_VSD(ICMODEL,EQMODEL,
     &                            ISERIES,IYEAR,IK,MEND,LUSEQ,
     &							TOC,CLAY,BULKDENS,
     &                            AVERAIN,AVEPET,AVETEMP,TOTPIC,
     &                            WMAX,WSAT,SOILW,DRAINW,
     &                            CNO3,CNH4,CMOBDOC,CO2,SOILPH,WTABLE,
     &                            ISIMP,IMPDEPTH,ISDYNPI)
C
C Subroutine to run ECOSSE for GIS
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXFERT				! Max.no.of fertiliser applications allowed 
	PARAMETER (MAXFERT=5)
	INTEGER MAXGROW				! Max.no.of growing seasons allowed
	PARAMETER (MAXGROW=300)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXLU1				! Max.no.of land use types
	DATA MAXLU1 /6/
	INTEGER MAXORGM				! Max.no.of organic manure applications allowed
	PARAMETER (MAXORGM=52)
	INTEGER MAXSOIL				! Max.no.of soil types
	PARAMETER (MAXSOIL=50)

	INTEGER IL					! Local counter variable
	INTEGER ISTEST				! Code for test results (1) or not (0)
C
C Variables passed to / from this subroutine
C
C ...Model descriptors
C
	INTEGER MEASLAY				! Layer that soil is measured to
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
	INTEGER ICMODEL				! OUT:Type of C initialisation model 
	INTEGER ICFIXED				! Initialisation of C is fixed
	INTEGER ICROTHCEQ			! Initialisation of C by RothC 
								!  equilibrium run
	DATA ICFIXED,ICROTHCEQ /1,2/
	INTEGER ISDYNPI				! Code for dynamic PI (adjusted by external factors)
	INTEGER ISDYNPI_OFF			! Adjustment of PI by external factors is off
	INTEGER ISDYNPI_ON			! Adjustment of PI by external factors is on
	DATA ISDYNPI_OFF, ISDYNPI_ON /0,1/
	INTEGER ISIMP				! Is impermeable layer? 0=No, 1=Yes
	REAL IMPDEPTH				! Depth of impermeable layer (cm)
	INTEGER ISPINUP				! Is N limitation spin-up used?
	INTEGER ISPINUP_OFF			! N limitation spin-up is not used
	INTEGER ISPINUP_ON			! N limitation spin-up is used
	DATA ISPINUP_OFF,ISPINUP_ON /0,1/
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
	INTEGER INVERT				! OUT(CULTIV):Invert / mix soil on cultivation? 0=No 1=Yes
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
	INTEGER CH4_RICHARDS    	! Richards CH4 model on
	INTEGER CH4_AITKENHEAD   	! Aitkenhead CH4 model on		
      DATA CH4_OFF,CH4_RICHARDS,CH4_AITKENHEAD /0,1,2/
	INTEGER PH_MODEL			! How is pH calculated?
	INTEGER PH_PARAM			! pH set to neutral or passed from parameter file
      INTEGER PH_STATIC			! pH is read in from input
								!   file and does not change
	INTEGER PH_FROM_VSD			! pH is calculated using VSD (a very 
								!   simple dynamic version of the MAGIC model
								!   by Ed Rowe & Chris Evans, CEH, Bangor)
	DATA PH_PARAM,PH_STATIC,PH_FROM_VSD /0,1,2/			 
C
C ...Timing factors
C
	INTEGER FIXEND				! IN:Fixed end? 0=No 1=Yes
	INTEGER IANTHES				! IN:No.timesteps from 01/01 to anthesis 
	INTEGER IEND				! IN:Timesteps from sowing date to harvest
	INTEGER IK					! OUT:No.timesteps from prev.harv. to current
	INTEGER IRYEAR				! IN:Current weather year number
	INTEGER IS_TS				! IN:Crop counter
	INTEGER IYEAR				! IN:Current growing season number
	INTEGER JSTOP				! IN:Code to stop simulation 0=No 1=Yes
	INTEGER LHARV				! IN:No.timesteps from 01/01 to prev. harvest
	INTEGER MEND				! IN:Number of timesteps from prev.harv.-harv.
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
	REAL CH4_PROD(MAXLAYER)  	! OUT:CH4 produced in each soil layer [kgC/ha/layer] (Richards CH4 model)
	REAL CH4_SOIL_OX			! OUT:CH4 oxidised from methane produced in the soil [kgC/ha] (Richards CH4 model)
      REAL CH4_ATMOS_OX			! OUT:CH4 oxidised from atmospheric CH4 [kgC/ha] (Richards CH4 model)
      REAL CH4_FLUX				! OUT:Net CH4 flux from soil column [kgC/ha] (Richards CH4 model)
      REAL CLAY(MAXLAYER)			! IN:Clay content of the layer (%)
	REAL CO2(MAXLAYER)			! IN:CO2 emitted (kgC/ha/layer)
	REAL CONC					! IN:Conc.N leached (kgN/ha/mm/layer)
	REAL CONC15					! IN:Conc.N15 leached (kgN/ha/mm/layer)
	REAL CRIT(MAXSOIL,MAXLAYER)	! IN:Critical N content (kgN/ha/layer)
	REAL CTOT					! IN:Total litter C input (kgC/ha)
	REAL CULTDEPTH				! OUT:Cultivation depth (cm)
	REAL DENIT					! IN:N lost by denitrification (kgN/ha)
	REAL DENITRIFN(MAXLAYER)	! IN/OUT: N lost by denitrification (kgN/ha/layer
	REAL DN15					! IN:N15 lost by denitrification (kgN15/ha)
	REAL DNIT(MAXLAYER)			! IN:Net Mineralised N (kgN/ha/lay/timestep)
	REAL DNIT15(MAXLAYER)		! IN:Net Mineralised N15 (kgN15/ha/lay/time)
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
      REAL DPMCTON(MAXLAYER)		! IN/OUT:DPM C:N ratio
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL DRAINW(MAXLAYER)		! IN:Water drained from this layer (mm/layer)
	REAL DRRAT					! IN/OUT:DPM:RPM ratio
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
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
	REAL LEACHDOC(MAXLAYER)		! IN/OUT:DOC leached from the layer (kgC/ha/layer/timestep)
	REAL LEACHDON(MAXLAYER)		! IN/OUT:DON leached from the layer (kgN/ha/layer/timestep)
	REAL NITRIFN(MAXLAYER)		! OUT:Nitrification from this layer (kg/ha/layer/timestep)
	INTEGER NSOIL				! IN:Soil code number
	REAL PIANN					! OUT: Total annual plant C input 
	                            !      (kg C / ha / yr) 
								!      - used as temprary value to allow adjustment of PI according 
								!        to external factors
      REAL PI_C(MAXLAYER)			! IN/OUT: Plant input C to soil (kgC/ha/step)
	REAL PI_CEQ_MON(12,MAXLAYER)! IN/OUT: Equilibrium plant C input each month
								!         in each layer (kgC/ha/month/layer)
      REAL PI_N(MAXLAYER)			! IN/OUT: Plant input N to soil (kgC/ha/step)
      REAL PI_N15(MAXLAYER)		! IN/OUT: Plant input N15 to soil(kgC/ha/step)
								!(1=arable,2=grassland,3=forestry,4=semi-natural)
	REAL REFIL(MAXLAYER)		! IN:Water needed to fill to FC (mm/lay)
	REAL RNIN					! IN:Litter N input in this timestep (kgN/ha)
	REAL RNIN15					! IN:Litter N15 input in timestep (kgN15/ha)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
      REAL RPMCTON(MAXLAYER)		! IN/OUT:RPM C:N ratio
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL RSTART					! IN:Initial N in debris pool (kgN/ha)
	REAL RSTART15				! IN:Initial N in debris pool (kgN15/ha)
	REAL SAND(MAXLAYER)			! IN:Sand content of the layer (%)
	REAL SATWATCONT(MAXLAYER)	! IN:Total water content at saturation (mm/layer)
	REAL SILT(MAXLAYER)			! IN:Silt content of the layer (%)
	REAL SLEA15(MAXLAYER)		! IN:Nitrate-N15 leached (kgN15/ha)
	REAL SLEACH(MAXLAYER)		! IN:Nitrate-N leached (kgN/ha)
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOILPH(MAXLAYER)		! IN/OUT:pH of soil in this layer
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL SRNIN					! IN:Total litter N input (kgN/ha)
	REAL SX						! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SXORGN					! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SXORGN15				! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL T15					! IN:Net Mineralised N (kgN/ha/layer)
	REAL TAM(MAXLAYER)			! IN:15N/N in the ammonium in the layer
	REAL TC						! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL TCINP					! IN:Total org.C input (kgC/ha) DONT PASS?
	REAL TDNIT					! IN:Net Mineralised N (kgN/ha/layer)
	REAL TFYMC					! IN:Total FYM C input (kgC/ha) DONT PASS?
	REAL TIN(MAXLAYER)			! IN:15N/N in the nitrate in the layer
	REAL TOC(MAXLAYER)			! IN:Measured total organic C (kgC/ha/layer)
	REAL TOTAST					! IN:Total ammonium in soil (kgN/ha)
	REAL TOTNST					! IN:Total nitrate in soil (kgN/ha)
	REAL TOTAST15				! IN:Total ammonium in soil (kgN15/ha)
	REAL TOTNST15				! IN:Total nitrate in soil (kgN15/ha)
	REAL TRNIN					! IN:Total litter N input (kgN/ha)
	REAL VIGOUR					! OUT: Vigour of cultivation (0-1) Determines the proportion of humus released to biomass,DPM and RPM pools
	REAL VOLAT					! IN:N lost by volatilisation (kgN/ha)
	REAL VOLAT15				! IN:N15 lost by volatilisation (kgN15/ha)
	REAL VOLATN(MAXLAYER)		! IN/OUT:Volatilisation from this layer (kg/ha/layer/timestep)
	REAL WILTPOINT(MAXLAYER)	! IN:Water content at wilting point (mm/layer)
	REAL WLEACH					! IN:N leached (kgN/ha)
	REAL WLEACH15				! IN:N15 leached (kgN15/ha)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
	REAL WTABLE					! IN:Water table depth in cm
C
C ... Crop factors
C
	REAL CACT					! IN:Crop N uptake (kgN/ha/timestep)
	REAL CACT15					! IN:Crop N15 uptake (kgN15/ha/timestep)
	REAL CACTOT					! IN:N taken up by crop (kgN/ha)
	REAL CATOT15				! IN:N15 taken up by crop (kgN15/ha)
	REAL CLOSSX					! IN:N lost by senescence (kgN/ha) 
	REAL CLOSSX15				! IN:N15 lost by senescence (kgN15/ha) 	
	REAL CUPTN					! IN:Crop N uptake (kgN/ha/timestep)
	REAL FIXN					! IN:N fixed from atmosphere (kgN/ha)
      INTEGER ICOVER				! OUT:Crop cover 1=Covered 0=Bare
	INTEGER ICROP				! IN:Current crop code
	INTEGER ISOWN				! IN:Timesteps from 01/01/01 to sowing date 
	INTEGER IHARV				! IN:Timesteps from 01/01/01 to harvest date 
      INTEGER LUSEQ(MAXGROW)		! Land use sequence
	INTEGER LU1					! OUT:First land use
	INTEGER NFERT				! IN:No.fertiliser applications to this crop
	INTEGER NORGM				! IN:No.manure applications to this crop
	INTEGER NSOW				! IN:Timesteps from prev.harv. to sowing date
	INTEGER NXYEARS				! IN:No.growing seasons simulated
	REAL OCROPN					! IN:N uptake of crop (kgN/ha) 
	REAL ORGN					! IN:Total org.N input (kgN/ha) 
	REAL PLANTUP(MAXLAYER)		! IN(RUN2_GIS_CROP)/OUT(TEST1_RES):Plant uptake from each soil layer
	REAL ROOT					! IN:Rooting depth according to restriction (cm)
	REAL SORGN					! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SORGN15				! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL TACTOT					! IN:Total crop N uptake (kgN/ha)
	INTEGER THISLU				! OUT:LU in the current growing season
	REAL TOTPIC(MAXLU)			! IN/OUT:Plant C input calculated using 
								!	  ROTHC equilibrium run for each LU (kgC/ha/yr)
	REAL WR						! IN:Root requirement (kgN/ha)
	REAL XORGN					! IN:Total org.N input (kgN/ha) 
C
C ... Fertiliser factors
C
	REAL FERT(MAXFERT)			! IN:Amount of fertiliser applied (kgN/ha)
	REAL FERTNF					! INOUT:Amount of fertiliser applied (kgN/ha)
	REAL FYMFERT				! IN:N input by fertiliser & FYM (kgN/ha)
	REAL FYMFERT15				! IN:N15 input by fertiliser & FYM (kgN15/ha)
	INTEGER IFERT(MAXFERT)		! IN:No.timesteps to split fert.application
	INTEGER ILAB(MAXFERT)		! IN:Labelling on fertiliser? 0=No 1=Yes
	INTEGER ILABNF				! INOUT:Labelling on fertiliser? 0=No 1=Yes
	INTEGER IVOL(MAXFERT)		! IN:N volatilised from fertiliser (kgN/ha)
	INTEGER IVOLNF				! INOUT:N volatilised from fertiliser (kgN/ha)
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
	REAL ORGMA(MAXORGM)			! IN:Amount of manure applied (t/ha)
      REAL ORGMJ(0:MAXGROW,MAXORGM)	! IN:Amount of manure applied (t/ha)
	INTEGER IOLABNF				! IN:Labelling on manure? 0=No 1=Yes
	INTEGER JORGMNF				! IN:Type of manure application
	REAL ORGMANF				! IN:Amount of manure applied (t/ha)
C
C ...Weather factors
C
	REAL AIRTEMPC				! IN: Air temperature (deg.C/timestep)	
	REAL LTA_AWC(12,MAXLAYER)	! IN:Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)	! IN:Average air temp this month (deg.C)
	REAL LAT					! IN: Latitude
	REAL RAIN					! IN: Rainfall (mm/timestep)
	REAL EVAP					! IN: Potential evap. (mm/timestep)
	REAL EVAPW					! IN: Potential evap. (mm/timestep)
      INTEGER IDEC				! OUT:Decade counter 
	REAL SOILTEMP(MAXLAYER)		! IN: Soil temperature (deg.C/timestep)
	REAL SEEDIN					! IN: Data used by MAGEC
C
C ...Dissolved organic matter factors
C
	REAL CO2FROMDOC(MAXLAYER)	! IN:CO2 output from DOC in each layer(kgC/ha)
	REAL MOBDOC(MAXLAYER)		! IN:Mobile DOC in each layer (kgC/ha)
	REAL MOBDON(MAXLAYER)		! IN:Mobile DON in each layer (kgN/ha)
C
C ... GIS specific variables
C
      REAL AVEPET(12)				! IN:Long term average PET (mm/month)
      REAL AVERAIN(12)			! IN:Long term average rainfall(mm/month)
      REAL AVETEMP(12)			! IN:Long term average monthly average 
  	REAL ICFACTOR(MAXLAYER)		! IN:Adjustment in rate needed to	
								!    achieve measured NPP and TOC	
	INTEGER ILU					! IN:Counter for land use 
      INTEGER ISAVE				! OUT:Code to save or retrieve variables
								!		0=save; 1=retrieve
	INTEGER ISERIES				! IN:Counter for soil series
	REAL PNREQ(MAXLU)			! IN/OUT:Amount of N required by plant each yr
	                            !        (kgN/ha/year)
C
C Results for VSD
C
	REAL CMOBDOC(MAXLAYER)		! IN: Concentration of mobile DOC in this layer (eq/m3)
	REAL CNH4(MAXLAYER)			! IN: Concentration of ammonium in this layer (eq/m3)
      REAL CNO3(MAXLAYER)			! IN: Concentration of nitrate in this layer (eq/m3)
C
C Required for Monte Carlo - 
C
     	REAL DPMCARB(MAXLAYER)
	REAL RPMCARB(MAXLAYER)
	REAL BCARB(MAXLAYER)
	REAL HCARB(MAXLAYER)
	REAL WRATEDM
	REAL TRATEM
	REAL DDAYS
CC
C Set IOM from TOC
C
      DO 100 IL=1,MAXLAYER
        CALL GET_IOM_FROM_FALLOON_EQN(TOC(IL),IOM(IL))
100   CONTINUE

C Retrieve state variables etc, that have been saved internally 
C
	SECONDS=(365.25/12)*24*60*60
      ISAVE=0
      CALL SAVE_FOR_VSD(ISAVE,
     &                  SOILN,SOIL15,AMMN,AMMN15,
     &                  DPMCARB0,DPMNIT0,DPMNLAB0,
     &                  RPMCARB0,RPMNIT0,RPMNLAB0,
     &                  BCARB0,BNIT0,BNLAB0,
     &                  HCARB0,HNIT0,HNLAB0,
     &                  NFERT,IFERT,ILAB,FERT,TFERT,IVOL,
     &                  C_TS,FLOWPROP,HZ1,RNIN,RNIN15,
     &                  PI_C,PI_CEQ_MON,PI_N,PI_N15,PNREQ,
     &	              CACT,CACT15,CACTOT,CATOT15,NSOW,ICFACTOR)
C
C Set restrictions on calculation
C Get parameters from TOC,IOM,CLAY,SILT,SAND,BULKDENS,	
C no restriction due to drainage,
C initialisation using ROTHC model & Bente Foereid's C:N ratios
C and get pH / water restriction on decomposition from NPP and TOC
C Use only 1 decade calculation
C
      MODTYPE=IS_SUNDIAL
	SPARMODEL=SPARCALC
C      ISPINUP=ISPINUP_ON
      ISPINUP=ISPINUP_OFF
	INMODEL=INPASSCN
	ITFUNC=ITFUNC_ROTHC
	IMFUNC=IMFUNC_ROTHC
	PH_MODEL=PH_FROM_VSD
	NORGM=0
	NSOIL=1
	IAWC=1
C
C Initialize THISFERT=current weeks fertilizer addition
C            WLEACH=current weeks water leaching
C      
      CALL SETWEEKVARS(JSTOP,CLOSSX,CLOSSX15,FYMFERT,FYMFERT15,
     &                     FIXN,THISFERT,THISFERT15,WLEACH,WLEACH15,
     &                     VOLAT,VOLAT15)
C
C Initialise soil and crop in first year
C	
      THISLU=LUSEQ(IYEAR)
      IF(IYEAR.EQ.1.AND.IK.EQ.1)THEN
C
C ...initialise the soil for initial land use,
C
        CALL INIT_GIS_SOILCN_NOOPT(DOCMODEL,
     &                             EQMODEL,SECONDS,
     &                             NSOIL,TOC,IOM,HZ1,CRIT,
     &                             SOILN,SOIL15,AMMN,AMMN15,
     &                             DPMCARB0,DPMNIT0,DPMNLAB0,
     &                             RPMCARB0,RPMNIT0,RPMNLAB0,
     &                             BCARB0,BNIT0,BNLAB0,
     &                             HCARB0,HNIT0,HNLAB0,
     &                             MOBDOC,MOBDON,THISLU,
     &                             CLAY,BULKDENS,PI_CEQ_MON,
     &                             LTA_AWC,LTA_TEMP,ITFUNC,IMFUNC,
     &                             WMAX,WSAT,ICFACTOR,PH_MODEL,SOILPH)
        IF(ISPINUP.EQ.ISPINUP_ON)THEN
          CALL SET_NLIM(THISLU,ICFACTOR,
     &                  SOILN,SOIL15,AMMN,AMMN15,SOILW,
     &                  DPMCARB0,DPMNIT0,DPMNLAB0,
     &                  RPMCARB0,RPMNIT0,RPMNLAB0,
     &                  BCARB0,BNIT0,BNLAB0,
     &                  HCARB0,HNIT0,HNLAB0,IOM,PI_CEQ_MON)
	  ENDIF
      ENDIF
C
C In first year of season....
C
      IF(IK.EQ.1)THEN
C
C Initialise current crop
C 
        CALL INIT_GIS_CROP(IYEAR,THISLU,SECONDS,
     &                     TOTPIC(THISLU),PNREQ(THISLU),
     &                     N_STEPS,ISOWN,NSOW,IHARV,LHARV,IANTHES,
     &                     MEND,IEND,NFERT,IFERT,ILAB,FERT,TFERT,
     &                     JSTOP,CTOT)
        IF(ISPINUP.EQ.ISPINUP_ON)THEN
	    CALL SET_NLIM_RATE(THISLU,ICFACTOR)
	  ENDIF
C
C Cultivate soil for gra->ara, for->ara, for->gra, nat->ara, for->gra, nat->for
C
	  IF(IYEAR.GT.1.)THEN
	    IF((LUSEQ(IYEAR).NE.LUSEQ(IYEAR-1)))THEN
	      CULTDEPTH=0
	      VIGOUR=0
              INVERT=0
	      IF(
     &	     (LUSEQ(IYEAR-1).EQ.2.AND.LUSEQ(IYEAR).EQ.1).OR.
     &         (LUSEQ(IYEAR-1).EQ.3.AND.LUSEQ(IYEAR).EQ.1).OR.
     &         (LUSEQ(IYEAR-1).EQ.4.AND.LUSEQ(IYEAR).EQ.1))THEN
	         CULTDEPTH=50
	         VIGOUR=0.5
	         INVERT=1
             ELSEIF(
     &              (LUSEQ(IYEAR-1).EQ.3.AND.LUSEQ(IYEAR).EQ.2).OR.
     &              (LUSEQ(IYEAR-1).EQ.3.AND.LUSEQ(IYEAR).EQ.4))THEN
	         CULTDEPTH=30
	         VIGOUR=0.5
	         INVERT=1
	       ELSEIF(
     &              (LUSEQ(IYEAR-1).EQ.4.AND.LUSEQ(IYEAR).EQ.3))THEN
	         CULTDEPTH=0
	         VIGOUR=0
	         INVERT=0
	       ENDIF
      	   CALL CULTIV(LUSEQ(IYEAR-1),
     &                   DPMCARB0,DPMNIT0,DPMNLAB0,
     &                   RPMCARB0,RPMNIT0,RPMNLAB0,
     &                   BCARB0,BNIT0,BNLAB0,
     &                   HCARB0,HNIT0,HNLAB0,
     &                   CULTDEPTH,VIGOUR,INVERT)
	    ENDIF	         
	  ENDIF
	ENDIF
C
C ...Calculate crop C and N returns and N offtake 
C
C
C Set soil temperature to air temperature
C
      CALL GETWEATHER_GIS(IK,LHARV,AVERAIN,AVEPET,AVETEMP, 
     &                            RAIN,EVAPW,SOILTEMP,AIRTEMPC)
C
C Adjust plant inputs for dynamic weather conditions
C
          PIANN=TOTPIC(THISLU)
          IF(ISDYNPI.EQ.ISDYNPI_ON)
     &      CALL ADJUST_PI(AVERAIN,AVETEMP,IK,LHARV,PIANN,RAIN,SOILTEMP)
C
C ...Calculate crop C and N returns and N offtake 
C
	CALL RUN1_GIS_CROP(IYEAR,IK,THISLU,IEND,MEND,NSOW,
     &                          N_STEPS,IS_TS,PIANN,
     &                          PNREQ(THISLU),CACTOT,CATOT15,
     &                          ICOVER,SXORGN,ORGN,
     &                          TRNIN,RNIN,RNIN15,C_TS,TCINP,WR)      
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
      CALL PARTITION_PLANTINPUT(C_TS,RNIN,RNIN15,
     &                                 PI_C,PI_N,PI_N15)
      IF(INMODEL.EQ.INPASSCN)
     &  CALL GET_CTON_FOEREID(DPMCTON,RPMCTON,THISLU)
	CALL SET_DPMRPMRATIO(THISLU,DRRAT)
	IF(ISIMP.EQ.1)CALL LIMIT_DRAIN(IMPDEPTH,SOILW)
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
     &					             DPMCARB,RPMCARB,BCARB,HCARB,
     &					             WRATEDM,TRATEM,SOILPH,MEASLAY)
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
      CALL DRAIN_SUNDIAL_WATER(RAIN,WMAX,SOILW,DRAINW,REFIL,
     &                              WSAT,FLOWPROP)
	CALL EVAP_SUNDIAL_WATER(EVAPW,SOILTEMP(1),ICOVER,
     &                            WMAX,SOILW)
      CALL PHYSICAL_SUNDIAL_SOILCN(DOCMODEL,SECONDS,
     &                                   NSOIL,SOILTEMP,
     &                                   WMAX,SOILW,DRAINW,REFIL,
     &                                   SLEACH,SLEA15,
     &                                   WLEACH,WLEACH15,CONC,CONC15,
     &                                   SOILN,SOIL15,AMMN,AMMN15,
     &                                   TIN,TAM,MOBDOC,MOBDON,
     &                                   LEACHDOC,LEACHDON)
C
C Calculate results as a concentration 
C
      DO 200 IL=1,MAXLAYER
	  CALL GETEQUIVALENTS(SOILN(IL),AMMN(IL),MOBDOC(IL),SOILW(IL),
     &                      CNO3(IL),CNH4(IL),CMOBDOC(IL))
200   CONTINUE
C
C Save state variables etc internally 
C
      ISAVE=1
      CALL SAVE_FOR_VSD(ISAVE,
     &                  SOILN,SOIL15,AMMN,AMMN15,
     &                  DPMCARB0,DPMNIT0,DPMNLAB0,
     &                  RPMCARB0,RPMNIT0,RPMNLAB0,
     &                  BCARB0,BNIT0,BNLAB0,
     &                  HCARB0,HNIT0,HNLAB0,
     &                  NFERT,IFERT,ILAB,FERT,TFERT,IVOL,
     &                  C_TS,FLOWPROP,HZ1,RNIN,RNIN15,
     &                  PI_C,PI_CEQ_MON,PI_N,PI_N15,PNREQ,
     &	              CACT,CACT15,CACTOT,CATOT15,NSOW,ICFACTOR)
C
C Leave RUN_ECOSSE_FOR_VSD
C
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE SAVE_FOR_VSD(ISAVE,
     &                  SOILN,SOIL15,AMMN,AMMN15,
     &                  DPMCARB0,DPMNIT0,DPMNLAB0,
     &                  RPMCARB0,RPMNIT0,RPMNLAB0,
     &                  BCARB0,BNIT0,BNLAB0,
     &                  HCARB0,HNIT0,HNLAB0,
     &                  NFERT,IFERT,ILAB,FERT,TFERT,IVOL,
     &                  C_TS,FLOWPROP,HZ1,RNIN,RNIN15,
     &                  PI_C,PI_CEQ_MON,PI_N,PI_N15,PNREQ,
     &	              CACT,CACT15,CACTOT,CATOT15,NSOW,ICFACTOR)
C
C Subroutine to save / retrieve results for VSD run
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXFERT				! Max.no.of fertiliser applications allowed 
	PARAMETER (MAXFERT=5)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)

	REAL AMMN1(MAXLAYER)		! Soil ammonium-N (kgN/ha/layer)
	REAL AMMN151(MAXLAYER)		! Soil ammonium-N15 (kgN15/ha/layer)
	REAL BCARB01(MAXLAYER)		! C in soil biomass (kgC/ha/layer)
      REAL BNIT01(MAXLAYER)		! N in soil biomass (kgN/ha/layer)
	REAL BNLAB01(MAXLAYER)		! N15 in soil humus (kgN15/ha/layer)
	REAL CACT1					! Crop N uptake (kgN/ha/timestep)
	REAL CACT151				! Crop N15 uptake (kgN15/ha/timestep)
	REAL CACTOT1				! N taken up by crop (kgN/ha)
	REAL CATOT151				! N15 taken up by crop (kgN15/ha)
	REAL C_TS1					! Litter C input in this timestep (kgC/ha)
	REAL DPMCARB01(MAXLAYER)	! C in decomposable PM (kgC/ha/layer)
	REAL DPMNIT01(MAXLAYER)		! N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB01(MAXLAYER)	! N15 in decomposable PM (kgN15/ha/layer)	
	REAL FERT1(MAXFERT)			! Amount of fertiliser applied (kgN/ha)
	REAL FLOWPROP1				! Proportion of flow needed to achieve
							    !	     water table at depth WTABLE
	REAL HCARB01(MAXLAYER)		! C in soil humus (kgC/ha/layer)
	REAL HNIT01(MAXLAYER)		! N in soil humus (kgN/ha/layer)
      REAL HNLAB01(MAXLAYER)		! N15 in soil biomass (kgN15/ha/layer)
	REAL HZ11(MAXLAYER)			! N:C of input PM for steady state 
  	REAL ICFACTOR1(MAXLAYER)	! Adjustment in rate needed to	
								!    achieve measured NPP and TOC	
	INTEGER IFC					! Fertiliser counter
	INTEGER IFERT1(MAXFERT)		! No.timesteps to split fert.application
	INTEGER ILAB1(MAXFERT)		! Labelling on fertiliser? 0=No 1=Yes
	INTEGER IFP					! Fertiliser parts counter
	INTEGER IVOL1(MAXFERT)		! N volatilised from fertiliser (kgN/ha)
	INTEGER IL					! Local layer counter
	INTEGER ILU					! Local land use counter
	INTEGER IMON				! Local month counter
	REAL LTA_AWC1(12,MAXLAYER)	! Long term available water (mm)
	REAL LTA_TEMP1(12,MAXLAYER)	! Average air temp this month (deg.C)
	INTEGER NFERT1				! No.fertiliser applications to this crop
	INTEGER NSOW1				! Timesteps from prev.harv. to sowing date
      REAL PI_C1(MAXLAYER)		! Plant input C to soil (kgC/ha/step)
	REAL PI_CEQ_MON1(12,MAXLAYER)! Equilibrium plant C input each month
								!  in each layer (kgC/ha/month/layer)
      REAL PI_N1(MAXLAYER)		! Plant input N to soil (kgC/ha/step)
      REAL PI_N151(MAXLAYER)		! Plant input N15 to soil(kgC/ha/step)
								!(1=arable,2=grassland,3=forestry,4=semi-natural)
	REAL PNREQ1(MAXLU)			! IN/OUT:Amount of N required by plant each yr
	                            !        (kgN/ha/year)
	REAL RNIN1					! Litter N input in this timestep (kgN/ha)
	REAL RNIN151				! Litter N15 input in timestep (kgN15/ha)
	REAL RPMCARB01(MAXLAYER)	! C in resistant PM (kgC/ha/layer)
	REAL RPMNIT01(MAXLAYER)		! N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB01(MAXLAYER)	! N15 in resistant PM (kgN15/ha/layer)
	REAL SOIL151(MAXLAYER)		! Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILN1(MAXLAYER)		! Soil nitrate-N (kgN/ha/layer)
	REAL TFERT1(MAXFERT,3)		! Prop.NO3,NH4,urea in fertiliser
C
C Variables passed to/from this subroutine
C
	REAL AMMN(MAXLAYER)			! IN/OUT:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN/OUT:Soil ammonium-N15 (kgN15/ha/layer)
	REAL BCARB0(MAXLAYER)		! IN/OUT:C in soil biomass (kgC/ha/layer)
      REAL BNIT0(MAXLAYER)		! IN/OUT:N in soil biomass (kgN/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN/OUT:N15 in soil humus (kgN15/ha/layer)
	REAL CACT					! IN:Crop N uptake (kgN/ha/timestep)
	REAL CACT15					! IN:Crop N15 uptake (kgN15/ha/timestep)
	REAL CACTOT					! IN:N taken up by crop (kgN/ha)
	REAL CATOT15				! IN:N15 taken up by crop (kgN15/ha)
	REAL C_TS					! IN/OUT:Litter C input in this timestep (kgC/ha)
	REAL DPMCARB0(MAXLAYER)		! IN/OUT:C in decomposable PM (kgC/ha/layer)
	REAL DPMNIT0(MAXLAYER)		! IN/OUT:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN/OUT:N15 in decomposable PM (kgN15/ha/layer)	
	REAL FERT(MAXFERT)			! IN:Amount of fertiliser applied (kgN/ha)
	REAL FLOWPROP				! IN:Proportion of flow needed to achieve
							    !	     water table at depth WTABLE
	REAL HCARB0(MAXLAYER)		! IN/OUT:C in soil humus (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN/OUT:N in soil humus (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN/OUT:N15 in soil biomass (kgN15/ha/layer)
	REAL HZ1(MAXLAYER)			! IN:N:C of input PM for steady state 
  	REAL ICFACTOR(MAXLAYER)		! IN:Adjustment in rate needed to	
								!    achieve measured NPP and TOC	
	INTEGER IFERT(MAXFERT)		! IN:No.timesteps to split fert.application
	INTEGER ILAB(MAXFERT)		! IN:Labelling on fertiliser? 0=No 1=Yes
	INTEGER IVOL(MAXFERT)		! IN:N volatilised from fertiliser (kgN/ha)
      INTEGER ISAVE				! IN:Code to save or retrieve variables
								!		0=save; 1=retrieve
	REAL LTA_AWC(12,MAXLAYER)	! IN:Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)	! IN:Average air temp this month (deg.C)
	INTEGER NFERT				! IN:No.fertiliser applications to this crop
	INTEGER NSOW				! IN:Timesteps from prev.harv. to sowing date
      REAL PI_C(MAXLAYER)			! IN/OUT: Plant input C to soil (kgC/ha/step)
	REAL PI_CEQ_MON(12,MAXLAYER)! IN/OUT: Equilibrium plant C input each month
								!         in each layer (kgC/ha/month/layer)
      REAL PI_N(MAXLAYER)			! IN/OUT: Plant input N to soil (kgC/ha/step)
      REAL PI_N15(MAXLAYER)		! IN/OUT: Plant input N15 to soil(kgC/ha/step)
								!(1=arable,2=grassland,3=forestry,4=semi-natural)
	REAL PNREQ(MAXLU)			! IN/OUT:Amount of N required by plant each yr
	                            !        (kgN/ha/year)
	REAL RNIN					! IN:Litter N input in this timestep (kgN/ha)
	REAL RNIN15					! IN:Litter N15 input in timestep (kgN15/ha)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL TFERT(MAXFERT,3)		! IN:Prop.NO3,NH4,urea in fertiliser
C
C Save parameters
C
      SAVE
	IF(ISAVE.EQ.1)THEN
	  DO 100 IL=1,MAXLAYER
	    AMMN1(IL)=AMMN(IL)
		AMMN151(IL)=AMMN15(IL)
		BCARB01(IL)=BCARB0(IL)
		BNIT01(IL)=BNIT0(IL)		
		BNLAB01(IL)=BNLAB0(IL)
		DPMCARB01(IL)=DPMCARB0(IL)
		DPMNIT01(IL)=DPMNIT0(IL)
		DPMNLAB01(IL)=DPMNLAB0(IL)
		HCARB01(IL)=HCARB0(IL)
		HNIT01(IL)=HNIT0(IL)
		HNLAB01(IL)=HNLAB0(IL)
		HZ11(IL)=HZ1(IL)
		PI_C1(IL)=PI_C(IL)
	    DO 200 IMON=1,12 
	      LTA_AWC1(IMON,IL)=LTA_AWC(IMON,IL)
	      LTA_TEMP1(IMON,IL)=LTA_TEMP(IMON,IL)
		  PI_CEQ_MON1(IMON,IL)=PI_CEQ_MON(IMON,IL)
200       CONTINUE
		PI_N1(IL)=PI_N(IL)
		PI_N151(IL)=PI_N15(IL)							
		RPMCARB01(IL)=RPMCARB0(IL)
		RPMNIT01(IL)=RPMNIT0(IL)
		RPMNLAB01(IL)=RPMNLAB0(IL)
		SOIL151(IL)=SOIL15(IL)
		SOILN1(IL)=SOILN(IL)
	    ICFACTOR1(IL)=ICFACTOR(IL)
100     CONTINUE
	  DO 300 IFC=1,MAXFERT
	    FERT1(IFC)=FERT(IFC)
	    IFERT1(IFC)=IFERT(IFC)
	    ILAB1(IFC)=ILAB(IFC)
	    IVOL1(IFC)=IVOL(IFC)
	    DO 400 IFP=1,3
	      TFERT1(IFC,IFP)=TFERT(IFC,IFP)
400       CONTINUE
300     CONTINUE
        DO 500 ILU=1,MAXLU
	    PNREQ1(ILU)=PNREQ(ILU)
500     CONTINUE
        C_TS1=C_TS
	  FLOWPROP1=FLOWPROP
	  RNIN1=RNIN
	  RNIN151=RNIN15
	  CACT1=CACT
	  CACT151=CACT15
	  CACTOT1=CACTOT
	  CATOT151=CATOT15
	  NFERT1=NFERT
	  NSOW1=NSOW
C
C Retrieve parameters
C   
      ELSEIF(ISAVE.EQ.0)THEN
	  DO 600 IL=1,MAXLAYER
	    AMMN(IL)=AMMN1(IL)
		AMMN15(IL)=AMMN151(IL)
		BCARB0(IL)=BCARB01(IL)
		BNIT0(IL)=BNIT01(IL)		
		BNLAB0(IL)=BNLAB01(IL)
		DPMCARB0(IL)=DPMCARB01(IL)
		DPMNIT0(IL)=DPMNIT01(IL)
		DPMNLAB0(IL)=DPMNLAB01(IL)
		HCARB0(IL)=HCARB01(IL)
		HNIT0(IL)=HNIT01(IL)
		HNLAB0(IL)=HNLAB01(IL)
		HZ1(IL)=HZ11(IL)
		PI_C(IL)=PI_C1(IL)
	    DO 700 IMON=1,12 
	      LTA_AWC(IMON,IL)=LTA_AWC1(IMON,IL)
	      LTA_TEMP(IMON,IL)=LTA_TEMP1(IMON,IL)
		  PI_CEQ_MON(IMON,IL)=PI_CEQ_MON1(IMON,IL)
700       CONTINUE
		PI_N(IL)=PI_N1(IL)
		PI_N15(IL)=PI_N151(IL)							
		RPMCARB0(IL)=RPMCARB01(IL)
		RPMNIT0(IL)=RPMNIT01(IL)
		RPMNLAB0(IL)=RPMNLAB01(IL)
		SOIL15(IL)=SOIL151(IL)
		SOILN(IL)=SOILN1(IL)
	    ICFACTOR(IL)=ICFACTOR1(IL)
600     CONTINUE
	  DO 800 IFC=1,MAXFERT
	    FERT(IFC)=FERT1(IFC)
	    IFERT(IFC)=IFERT1(IFC)
	    ILAB(IFC)=ILAB1(IFC)
	    IVOL(IFC)=IVOL1(IFC)
	    DO 900 IFP=1,3
	      TFERT(IFC,IFP)=TFERT1(IFC,IFP)
900       CONTINUE
800     CONTINUE
        C_TS=C_TS1
	  FLOWPROP=FLOWPROP1
	  RNIN=RNIN1
	  RNIN15=RNIN151
	  CACT=CACT1
	  CACT15=CACT151
	  CACTOT=CACTOT1
	  CATOT15=CATOT151
	  NFERT=NFERT1
	  NSOW=NSOW1
        DO 1000 ILU=1,MAXLU
	    PNREQ(ILU)=PNREQ1(ILU)
1000    CONTINUE
	ENDIF
C
C Leave SAVE_FOR_VSD
C
      END


