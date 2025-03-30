C------------------------------------------------------------
C
C GIS routines for ECOSSE Carbon and Nitrogen Turnover Model
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
C ECOSSE_EQRUN()
C ECOSSE_GIS_RUN()
C-------------------------------------------------------------
C
      SUBROUTINE ECOSSE_EQRUN(AMMN,AMMN15,ATM,AVEPET,AVERAIN,AVETEMP,
     &                BCARB0,BNIT0,BNLAB0,CACCUM,CH4MODEL,CNMODEL,CTOT,
     &			    CLOSSX,CLOSSX15,CRIT,DMODEL,DOCMODEL,DPMCARB0,
     &			    DPMNIT0,DPMNLAB0,DRAINW,FERT,FIXN,FLOWPROP,
     &                FYMFERT,FYMFERT15,HCARB0,HNIT0,HNLAB0,
     &                IANTHES,IEND,IFERT,IHARV,ILAB,IMFUNC,ICMODEL,
     &                INMODEL,IOLAB,IOM,IORGM,ISERIES,ISOWN,ITFUNC,IVOL,
     &                JORGM,JSTOP,LHARV,LU1,MEND,
     &                NFERT,NORGM,NSOIL,NSOW,ORGMA,
     &                PNREQ,REFIL,RPMCARB0,RPMNIT0,RPMNLAB0,
     &                SOIL15,SOILN,SOILW,TAM,TFERT,THISFERT,
     &                THISFERT15,TIN,TOC,TOTPIC,VOLAT,VOLAT15,WLEACH,
     &                WLEACH15,WSAT,WMAX,SOILPH,NXYEARS,EQMODEL,CLAY,
     &                PI_CEQ_MON,DOMSOILISIMP,DOMSOILIMPDEPTH,
     &                SOMDEPTH,NSOMLAY,LTA_AWC,LTA_TEMP,PH_MODEL)		  
C
C Modify plant inputs to achieve steady state in full ECOSSE run
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/
	INTEGER MAXFERT				! Max.no.of fertiliser applications allowed 
	PARAMETER (MAXFERT=5)
	INTEGER MAXGROW				! Max.no.of growing seasons allowed
	PARAMETER (MAXGROW=300)
	INTEGER MAXITER1			! Maximum allowed number of iterations
      DATA MAXITER1 /1000/
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER MAXLUSOIL   		! Max.no.of land use types
	PARAMETER (MAXLUSOIL=6)
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=10)
	INTEGER MAXORGM				! Max.no.of organic manure applications allowed
	PARAMETER (MAXORGM=52)
	INTEGER MAXSEQ				! Max.no of sequences included in GIS
	PARAMETER (MAXSEQ=MAXLU*MAXLU)						
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXSOIL				! Max.no.of soil types
	PARAMETER (MAXSOIL=50)
	INTEGER MAXSOMLAY			! Max.no.of soil layers entered			
	PARAMETER (MAXSOMLAY=10)

      REAL*8 DC1(MAXLAYER)		! Change in C over simulation with MSS = 1 (kg C / ha)
      REAL*8 DC2(MAXLAYER)		! Change in C over simulation with calculated MSS (kg C / ha)
	REAL*8 DCTOT				! Change in C over simulation over whole measured profile (kg C / ha)
	REAL DELTAC1				! Change in C during last time step (kgC/ha/yr)
	REAL DELTAC2				! Change in C during this time step (kgC/ha/yr)
	REAL ENDC					! C at end of year(kg C / ha)
      REAL ENDB					! C in BIO at end of year(kg C / ha)
	REAL ENDD					! C in DPM at end of year(kg C / ha)
	REAL ENDH					! C in HUM at end of year(kg C / ha)
	REAL ENDR					! C in RPM at end of year(kg C / ha)
      REAL FNULL					
      REAL FNULLARRAY(MAXLAYER)					
	INTEGER ICFIXED				! Initialisation of C is fixed
	INTEGER ICROTHCEQ			! Initialisation of C by RothC 
								!  equilibrium run
	DATA ICFIXED,ICROTHCEQ /1,2/
	INTEGER IK					! No.timesteps from prev.harv. to current
      INTEGER IL					! Local layer counter
	INTEGER INSTABLECN			! Initisalisation of N by assuming 
								!  steady state (after Bradbury)
	INTEGER INPASSCN			! Initialisation of N by passing the
								!  C:N ratio of DPM and RPM
	DATA INSTABLECN,INPASSCN /1,2/
	INTEGER INULL				
	INTEGER ISATEQ				! Code to indicate whether equilibrium 
								! has been found 0=No 1=Yes
	INTEGER ISPEAT				! 1=peat, 0=non-peat
	INTEGER ISTEST				! Output balance sheet files 
								! 0 = Do not output test files
	                            ! 1 = Output balance sheet files (BALANCE_N.OUT, BALANCE_C.OUT & SOILW.OUT)
	                            ! 2 = C Change  from one land use to aother
	                            ! 3 = Put solution results
	INTEGER IYEAR				! Current growing season number - set to 1
	REAL LINEA					! Parameter a in fitted line (y = a+bx) for optimisation 
	REAL LINEB					! Parameter b in fitted line (y = a+bx) for optimisation 
	REAL MEASC					! Measured C at end of optimisation iteration (kg C / ha)
	INTEGER MEASLAY				! Layer that soil is measured to
	INTEGER METMON				! Month selected for met data
	INTEGER METYEAR				! Year selected for met data
	INTEGER NITER				! No.of iterations in plant input optimisation
	INTEGER N_STEPS				! No.timesteps in a year
	REAL PC_TS(12)				! Proportion of total C litter input added each timestep (kgC/ha)
	REAL PDPM					! Proportion of DPM
	REAL PRPM					! Proportion of RPM
      REAL PI_CEQ_IN(12,MAXLAYER)	! Total annual plant input before modification (kg C / ha)
      REAL PIC_LAST(MAXLAYER)		! Total annual plant input before modification in last iteration (kg C / ha)
	REAL PINC					! Proportion increase in soil C for degrading soils
	REAL PITOT					! Total annual plant input (kg C /ha)
	INTEGER REPYR				! Reporting year during optimisation
	REAL SECONDS				! Number of seconds in one timestep
	REAL SIMC(MAXLAYER)			! Simulated C at end of optimisation iteration (kg C / ha)
	REAL SIMC_LAST				! Simulated C at end of last optimisation iteration (kg C / ha)	
	REAL STARTC					! C at start of year (kg C / ha)
	REAL STARTC1				! C at start of year 1 (kg C / ha)
	REAL STARTD					! C in DPM at start of year 1 (kg C / ha)
	REAL STARTR					! C in RPM at start of year 1 (kg C / ha)
      REAL STARTPIC(MAXLAYER)		! Total annual plant input before modification (kg C / ha)
	REAL SUMC_TS				! Total C litter input added each year (kgC/ha/year)
      INTEGER TESTLU				! LU to be output
      INTEGER THISLU				! This LU
      INTEGER VECTOR_OPT			! Optimise by vectors (1) or using proportions (0)
C
C Variables passed to/from other subroutines
C
	REAL AIRTEMP				! IN: Air temperature (deg.C/timestep)	
	REAL ALPHA(MAXLAYER)			! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
	REAL ANIT					! Ammonium immobilised (kgN/ha)
	REAL ANIT15					! Ammonium N15 nitrified (kgN/ha)
      REAL ATM					! IN:Atmospheric N input (kgN/ha/timestep)
	REAL ATM15					! IN:Atmospheric N15 input (kgN15/ha/time)
      REAL AVEPET(12)				! IN:Long term average PET (mm/month)
      REAL AVERAIN(12)			! IN:Long term average rainfall(mm/month)
      REAL AVETEMP(12)			! IN:Long term average monthly average 
      REAL BALANCE(20)			! IN:C Results
	REAL BCARB(MAXLAYER)
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
	REAL BULKDENS(MAXLAYER)				! Soil bulk density (g cm-3)
	REAL BYRAIN					! IN/OUT:Rain lost by bypass flow mm/timestep	
	REAL C_TS					! IN:Litter C input in this timestep (kgC/ha)
      REAL CACCUM					! IN: Measured annual C accumulation 
								!     in the soil kgC/ha/yr
	REAL CACT					! IN:Crop N uptake (kgN/ha/timestep)
	REAL CACT15					! IN:Crop N15 uptake (kgN15/ha/timestep)
	REAL CACTOT					! IN:N taken up by crop (kgN/ha)
	REAL CATOT15				! IN:N15 taken up by crop (kgN15/ha)
	REAL CH4(MAXLAYER)			! IN:CH4 emitted (kgC/ha/layer)
	INTEGER CH4MODEL			! OUT:Methane model (off, Richards or Aitkenhead model) 
	REAL CH4TOAIR				! OUT:CH4 released to atmosphere (kgC/ha)
	REAL CH4_PROD(MAXLAYER)  	! OUT: CH4 produced in each soil layer [kgC/ha/layer] (Richards CH4 model)
	REAL CH4_SOIL_OX			! OUT: CH4 oxidised from methane produced in the soil [kgC/ha] (Richards CH4 model)
      REAL CH4_ATMOS_OX			! OUT: CH4 oxidised from atmospheric CH4 [kgC/ha] (Richards CH4 model)
      REAL CH4_FLUX				! OUT: Net CH4 flux from soil column [kgC/ha] (Richards CH4 model)
      REAL CLAY(MAXLAYER)			! IN:Clay content of the layer (%)
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: % clay content in this layer
	REAL CLOSSX					! IN:N lost by senescence (kgN/ha) 
	REAL CLOSSX15				! IN:N15 lost by senescence (kgN15/ha) 	
	INTEGER CNMODEL				! IN:Type of Model for C:N ratio of DPM and 
	                            !     RPM (MAGEC or Foereid)
	REAL CO2(MAXLAYER)			! IN:CO2 emitted (kgC/ha/layer)
	REAL CO2FROMDOC(MAXLAYER)	! IN:CO2 output from DOC in each layer(kgC/ha)
	REAL CONC					! IN:Conc.N leached (kgN/ha/mm/layer)
	REAL CONC15					! IN:Conc.N15 leached (kgN/ha/mm/layer)
      REAL CONVER_F					! IN/OUT:Conversion to correct timestep
	REAL CRIT(MAXSOIL,MAXLAYER)	! IN:Critical N content (kgN/ha/layer)
	REAL CTOT					! IN:Total litter C input (kgC/ha)
	REAL CUPTN					! IN:Crop N uptake (kgN/ha/timestep)
	REAL DELTA(MAXLAYER)			! IN/OUT: Prop.HUM produced on HUM decompn 
	REAL DENIT					! IN:N lost by denitrification (kgN/ha)
	REAL DENITRIFN(MAXLAYER)	! IN/OUT: N lost by denitrification (kgN/ha/layer
	REAL DFACT(MAXLAYER)		! IN/OUT: Denitrification factor
	INTEGER DMODEL				! IN:Type denitrification model chosen 
	                            ! (Bradbury or NEMIS)
	REAL DN15					! IN:N15 lost by denitrification (kgN15/ha)
	REAL DNIT(MAXLAYER)			! IN:Net Mineralised N (kgN/ha/lay/timestep)
	REAL DNIT15(MAXLAYER)		! IN:Net Mineralised N15 (kgN15/ha/lay/time)
	INTEGER DOCMODEL			! IN:DOC model (on or off)
	REAL DOMSOILIMPDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN: Depth of impermeable layer if present (cm)
	INTEGER DOMSOILISIMP(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN: Flag for impermeable layer (0=No; 1=Yes)
     	REAL DPMCARB(MAXLAYER)
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
      REAL DPMCTON(MAXLAYER)		! OUT:DPM C:N ratio
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL DPMRATE(MAXLAYER)			! IN/OUT: Rate constant for DPM decompn
	REAL DRAINW(MAXLAYER)		! IN:Water drained from this layer (mm/layer)
	REAL DRRAT						! IN:DPM:RPM ratio
	INTEGER EQMODEL				! IN:Type of equilibrium run 
	                            ! (NPP, TOC or both)
	REAL EVAPW					! IN: Potential evap. (mm/timestep)
	REAL FANIT					! Fertiliser N nitrified (kgN/ha)
	REAL FANIT15				! Fertiliser N15 nitrified (kgN/ha)
	REAL FERT(MAXFERT)			! IN:Amount of fertiliser applied (kgN/ha)
	REAL FERTNF					! INOUT:Amount of fertiliser applied (kgN/ha)
	REAL FIXN					! IN:N fixed from atmosphere (kgN/ha)
	REAL FLOWPROP				! IN:Proportion of flow needed to achieve
							    !	     water table at depth WTABLE
      REAL FPROPX(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
	REAL FYMAX(MAXORGM)			! IN/OUT:Proportion of ammonium-N in FYM
	REAL FYMCX(MAXORGM)			! IN/OUT:Proportion of C in FYM
	REAL FYMFERT				! IN:N input by fertiliser & FYM (kgN/ha)
	REAL FYMFERT15				! IN:N15 input by fertiliser & FYM (kgN15/ha)
	REAL FYMLOSS(MAXORGM)		! IN/OUT:Prop.FYM lost by volatn/timestep
	REAL FYMNX(MAXORGM)			! IN/OUT:Proportion of Organic N in FYM
	REAL FYMPOS(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
      REAL FYMSTART(MAXORGM)		! IN/OUT:Amount of rainfall in 1 week, 
								!		 below which volatilisation will occur
	REAL FYMWX(MAXORGM)			! IN/OUT:Amount of water added in FYM

	REAL GAMMA(MAXLAYER)			! IN/OUT: Prop.BIO produced on HUM decompn 
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
	REAL HCARB(MAXLAYER)
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
	INTEGER IANTHES				! IN:No.timesteps from 01/01 to anthesis 
  	REAL ICFACTOR(MAXLAYER)			! IN:Adjustment in rate needed to	
									!    achieve measured NPP and TOC	
	INTEGER ICMODEL				! OUT:Type of C initialisation model 
  	REAL ICNULL(MAXLAYER)			! IN:Adjustment in rate needed to	
									!    achieve measured NPP and TOC - not used	
      INTEGER ICOVER				! IN:Crop cover 1=Covered 0=Bare
	INTEGER IEND				! IN:Timesteps from sowing date to harvest
	INTEGER IFERT(MAXFERT)		! IN:No.timesteps to split fert.application
	INTEGER IHARV				! IN:Timesteps from 01/01/01 to harvest date 
	INTEGER ILAB(MAXFERT)		! IN:Labelling on fertiliser? 0=No 1=Yes
	INTEGER ILABNF				! INOUT:Labelling on fertiliser? 0=No 1=Yes
	INTEGER IMFUNC				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
	INTEGER INMODEL				! OUT:Type of N initialisation model
	INTEGER IOLAB(MAXORGM)		! IN:Labelling on manure? 0=No 1=Yes
	INTEGER IOLABNF				! IN:Labelling on manure? 0=No 1=Yes
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
	REAL IOMARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium IOM (kgC/ha/layer)
	INTEGER IORGM(MAXORGM)		! IN:No.timesteps to split manure application
	INTEGER IS_TS				! IN:Crop counter
      INTEGER ISAVE				! OUT:Code to save or retrieve variables
								!		0=save; 1=retrieve
	INTEGER ISERIES				! IN:Counter for soil series
	INTEGER ISOWN				! IN:Timesteps from 01/01/01 to sowing date 
	INTEGER ITFUNC				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IVOL(MAXFERT)		! IN:N volatilised from fertiliser (kgN/ha)
	INTEGER IVOLNF				! INOUT:N volatilised from fertiliser (kgN/ha)
	INTEGER JORGMNF				! IN:Type of manure application
	INTEGER JORGM(MAXORGM)		! IN:Type of manure application
	INTEGER JSTOP				! IN:Code to stop simulation 0=No 1=Yes
	INTEGER LHARV				! IN:No.timesteps from 01/01 to prev. harvest
	REAL LEACHDOC(MAXLAYER)		! IN/OUT:DOC leached from the layer (kgC/ha/layer/timestep)
	REAL LEACHDON(MAXLAYER)		! IN/OUT:DON leached from the layer (kgN/ha/layer/timestep)
      INTEGER LU1					! IN:Land use code
      INTEGER LUARRAY(MAXSOIL)	! IN/OUT:Land use before equilibrium 
	REAL LTA_AWC(12,MAXLAYER)	! IN:Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)	! IN:Average air temp this month (deg.C)
	INTEGER MEND				! IN:Number of timesteps from prev.harv.-harv.
	REAL MOBDOC(MAXLAYER)		! IN:Mobile DOC in each layer (kgC/ha)
	REAL MOBDON(MAXLAYER)		! IN:Mobile DON in each layer (kgN/ha)
	INTEGER NFERT				! IN:No.fertiliser applications to this crop
	REAL NITRIFN(MAXLAYER)		! OUT:Nitrification from this layer (kg/ha/layer/timestep)
	INTEGER NORGM				! IN:No.manure applications to this crop
	INTEGER NSOIL				! IN:Soil code number
	INTEGER NSOMLAY(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN:Number of SOM layers
	INTEGER NSOW				! IN:Timesteps from prev.harv. to sowing date
      INTEGER NUMSOIL					! Number of soils defined
	INTEGER NXYEARS				! IN:No.growing seasons simulated
	REAL OCROPN					! IN:N uptake of crop (kgN/ha) 
	REAL ORGMA(MAXORGM)			! IN:Amount of manure applied (t/ha)
	REAL ORGMANF				! IN:Amount of manure applied (t/ha)
	REAL ORGN					! IN:Total org.N input (kgN/ha) 
	INTEGER PH_MODEL			! IN:How is pH calculated?
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: pH above which decomp.max.

      REAL PI_C(MAXLAYER)			! IN: Plant input C to soil (kgC/ha/step)
	REAL PI_CEQ_MON(12,MAXLAYER)! IN/OUT: Equilibrium plant C input each month
								!         in each layer (kgC/ha/month/layer)
      REAL PI_N(MAXLAYER)			! IN: Plant input N to soil (kgC/ha/step)
      REAL PI_N15(MAXLAYER)		! IN: Plant input N15 to soil (kgC/ha/step)
	REAL PLANTUP(MAXLAYER)		! IN(RUN2_GIS_CROP)/OUT(TEST1_RES):Plant uptake from each soil layer
	REAL PNREQ(MAXLU)			! IN:Amount of N required by plant each yr
	                            !        (kgN/ha/year)
	REAL RAIN					! IN: Rainfall (mm/timestep)
	REAL REFIL(MAXLAYER)		! IN:Water needed to fill to FC (mm/lay)
	REAL RNIN					! IN:Litter N input in this timestep (kgN/ha)
	REAL RNIN15					! IN:Litter N15 input in timestep (kgN15/ha)
	REAL RPMCARB(MAXLAYER)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
      REAL RPMCTON(MAXLAYER)		! OUT:RPM C:N ratio
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL RPMRATE(MAXLAYER)			! IN/OUT: Rate constant for RPM decompn
	REAL SATWATCONT(MAXLAYER)	    ! IN:Total water content at saturation (mm/layer)
	REAL SEEDIN					! IN: Data used by MAGEC
	REAL SLEA15(MAXLAYER)		! IN:Nitrate-N15 leached (kgN15/ha)
	REAL SLEACH(MAXLAYER)		! IN:Nitrate-N leached (kgN/ha)
	CHARACTER*40 SNAME(MAXSOIL)		! IN/OUT:Soil name
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOILPH(MAXLAYER)		! IN/OUT:pH of soil in this layer
	REAL SOILTEMP(MAXLAYER)		! IN: Soil temperature (deg.C/timestep)
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL SOMDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN:depth of soil organic matter layers
	INTEGER SPINCYCLE,temp			! OUT:Number of years in the spin-up
	REAL SRNIN					! IN:Total litter N input (kgN/ha)
	REAL SX						! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SXORGN					! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL T15					! IN:Net Mineralised N (kgN/ha/layer)
	REAL TACTOT					! IN:Total crop N uptake (kgN/ha)
	REAL TAM(MAXLAYER)			! IN:15N/N in the ammonium in the layer
	REAL TC						! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL TCINP					! IN:Total org.C input (kgC/ha) DONT PASS?
	REAL TDNIT					! IN:Net Mineralised N (kgN/ha/layer)
	REAL TFERT(MAXFERT,3)		! IN:Prop.NO3,NH4,urea in fertiliser
	REAL TFERTNF(3)				! INOUT:Prop.NO3,NH4,urea in fertiliser
	REAL TFYMC					! IN:Total FYM C input (kgC/ha) DONT PASS?
	REAL THISFERT				! IN:N input by fertiliser (kgN/ha)
	REAL THISFERT15				! IN:N15 input by fertiliser (kgN15/ha)
	REAL TIN(MAXLAYER)			! IN:15N/N in the nitrate in the layer
	REAL TOC(MAXLAYER)			! IN:Measured total organic C (kgC/ha/layer)
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium TOC (kgC/ha/layer)
	REAL TOTPIC(MAXLU)			! IN:Plant C input calculated using 
	REAL TRATEM
	REAL TRNIN					! IN:Total litter N input (kgN/ha)
	REAL VOLAT					! IN:N lost by volatilisation (kgN/ha)
	REAL VOLAT15				! IN:N15 lost by volatilisation (kgN15/ha)
	REAL VOLATN(MAXLAYER)		! IN/OUT:Volatilisation from this layer (kg/ha/layer/timestep)
	REAL WILTPOINT(MAXLAYER)	! IN:Water content at wilting point (mm/layer)
	REAL WLEACH					! IN:N leached (kgN/ha)
	REAL WLEACH15				! IN:N15 leached (kgN15/ha)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WR						! IN:Root requirement (kgN/ha)
	REAL WRATEDM
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
	REAL XORGN					! IN:Total org.N input (kgN/ha) 
	REAL DDAYS					! Degree days
C
C Set as peat (ie all DPM and RPM) (1=yes, 0=no)
C
202   CONTINUE
      IF(CACCUM.LT.0)THEN
	  PRINT*,'How do you want to start the optimisation?'
	  PRINT*,'    1. Assume most C is undecomposed plant material'
	  PRINT*,'    2. Assume initial pool ratios calculated by RothC'
	  READ(*,*)ISPEAT
	  IF(ISPEAT.LT.1.AND.ISPEAT.GT.2)THEN
	    PRINT*,'Error in response'
	    GOTO 202
	  ENDIF
      ENDIF
C
C ......retrieve soil parameters
C	
      IF(ISPEAT.EQ.1)THEN	
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
	ENDIF
C
C Set timestep to monthly
C
      SECONDS=(365.25/12)*60*60*24
	N_STEPS=12
	CONVER_F=SECONDS/(7*24*60*60)
C
C Calculate plant input at start of equilibrium run
C
	MEASLAY=SOMDEPTH(ISERIES,LU1,NSOMLAY(ISERIES,LU1))
	MEASLAY=MEASLAY*MAXLAYER1/MAXDEPTH
	DO 200 IL=1,MAXLAYER1
	  STARTPIC(IL)=0
	  DO 300 IK=1,12
	    IF(IL.GT.MEASLAY)THEN
	      PI_CEQ_MON(IK,IL)=0
	      DPMCARB0(IL)=0
	      DPMNIT0(IL)=0
	      DPMNLAB0(IL)=0
	      RPMCARB0(IL)=0
	      RPMNIT0(IL)=0
	      RPMNLAB0(IL)=0
	      BCARB0(IL)=0
	      BNIT0(IL)=0
	      BNLAB0(IL)=0
	      HCARB0(IL)=0
	      HNIT0(IL)=0
	      HNLAB0(IL)=0
	      AMMN(IL)=0
	      AMMN15(IL)=0
	      SOILN(IL)=0
	      SOIL15(IL)=0
	    ENDIF
	    STARTPIC(IL)=STARTPIC(IL)+PI_CEQ_MON(IK,IL)
	    PI_CEQ_IN(IK,IL)=PI_CEQ_MON(IK,IL)
300     CONTINUE
200   CONTINUE
C
C Calculate measured C
C
      MEASC=0
	DO 400 IL=1,MEASLAY
	  MEASC=MEASC+TOC(IL)-IOM(IL)
400   CONTINUE
C
C If agrading, start off with empty soil organic matter pools, and allow it to agrade 
C If degrading, start with RothC pools plus annual degradation
C
      NITER=0
      DO 100 IL=1,MAXLAYER
	  IF(CACCUM.GE.0.0001)THEN
	    DPMCARB0(IL)=0
	    RPMCARB0(IL)=0
	    BCARB0(IL)=0
	    HCARB0(IL)=0
	    DPMNIT0(IL)=0
	    RPMNIT0(IL)=0
	    BNIT0(IL)=0
	    HNIT0(IL)=0
	    DPMNLAB0(IL)=0
	    RPMNLAB0(IL)=0
	    BNLAB0(IL)=0
	    HNLAB0(IL)=0
	  ELSEIF(CACCUM.LT.-0.0001)THEN
	    PINC=(MEASC-(CACCUM))/MEASC
          DPMCARB0(IL)=PINC*DPMCARB0(IL)
          DPMNIT0(IL)=PINC*DPMNIT0(IL)
		DPMNLAB0(IL)=PINC*DPMNLAB0(IL)
          RPMCARB0(IL)=PINC*RPMCARB0(IL)
		RPMNIT0(IL)=PINC*RPMNIT0(IL)
		RPMNLAB0(IL)=PINC*RPMNLAB0(IL)
          BCARB0(IL)=PINC*BCARB0(IL)
		BNIT0(IL)=PINC*BNIT0(IL)
		BNLAB0(IL)=PINC*BNLAB0(IL)
          HCARB0(IL)=PINC*HCARB0(IL)
		HNIT0(IL)=PINC*HNIT0(IL)
		HNLAB0(IL)=PINC*HNLAB0(IL)
	  ENDIF
100   CONTINUE
      PRINT*,'ECOSSE initialisation'
	PRINT*,'=============================='
C
C Save current state of the soil 
C
      ISAVE=1 ! Save soil 
      CALL SAVE_GIS_SOIL(ISAVE,LU1,
     &                   SOILN,SOIL15,AMMN,AMMN15,
     &                   SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
     &                   DPMCARB0,DPMNIT0,DPMNLAB0,
     &                   RPMCARB0,RPMNIT0,RPMNLAB0,
     &                   BCARB0,BNIT0,BNLAB0,
     &                   HCARB0,HNIT0,HNLAB0,
     &                   IOM,PI_CEQ_MON)
C
C ...Calculate crop C and N returns and N offtake 
C
      IYEAR=1
	MEND=IHARV
      THISLU=LU1
C******************************
C Redistribute plant inputs according to SUNDIAL calculations
C ...Initialise current crop
C 
      CALL INIT_GIS_CROP(IYEAR,THISLU,SECONDS,
     &                     TOTPIC(THISLU),PNREQ(THISLU),
     &                     N_STEPS,ISOWN,NSOW,IHARV,LHARV,IANTHES,
     &                     MEND,IEND,NFERT,IFERT,ILAB,FERT,TFERT,
     &                     JSTOP,CTOT)
C
C ...Calculate crop C and N returns and N offtake 
C
      DO 500 IK=1,MEND
	  CALL RUN1_GIS_CROP(IYEAR,IK,THISLU,IEND,MEND,NSOW,
     &                          N_STEPS,IS_TS,TOTPIC(THISLU),
     &                          PNREQ(THISLU),CACTOT,CATOT15,
     &                          ICOVER,SXORGN,ORGN,
     &                          TRNIN,RNIN,RNIN15,C_TS,TCINP,WR)      
C
C ...Save carbon added this timestep
C
        PC_TS(IK)=C_TS
C
C ...Sum carbon added over year
C
        IF(IK.EQ.1)SUMC_TS=0
	  SUMC_TS=SUMC_TS+C_TS
500   CONTINUE
C
C ...Work out proportion of total C added each timestep
C
      DO 600 IK=1,MEND
	  PC_TS(IK)=PC_TS(IK)/SUMC_TS
C
C ...Correct distribution over year to match SUNDIAL equations
C
        DO 700 IL=1,MAXLAYER1
	    PI_CEQ_IN(IK,IL)=STARTPIC(IL)*PC_TS(IK)
	    PI_CEQ_MON(IK,IL)=STARTPIC(IL)*PC_TS(IK)
700     CONTINUE
600   CONTINUE
      
	PITOT=0
      DO 750 IL=1,MEASLAY
        PITOT=PITOT+STARTPIC(IL)
750   CONTINUE
      PRINT*,'Start: Meas. C = ',MEASC,' dC=',CACCUM
C******************************
C Full ECOSSE simulation using
C
101   CONTINUE
      NITER=NITER+1
	WRITE(*,20)NITER,PITOT
20       FORMAT(/'Iteration ',I3,' PI ',F10.0/'   YR       DPM        '
     &          'RPM        BIO        HUM        TOT(kgC/ha) dC    ')
C
C Adjust plant input according to modifier
C
      TOTPIC(LU1)=0
      DO 800 IL=1,MAXLAYER1
	  DO 900 IK=1,12
	    TOTPIC(LU1)=TOTPIC(LU1)+PI_CEQ_MON(IK,IL)
900     CONTINUE
800   CONTINUE
C
C Set maximum number of cycles in spin-up
C
      SPINCYCLE=20000
C
C For each growing season... 
C
      IYEAR=1
      DO WHILE (IYEAR.LE.SPINCYCLE)
C
C Calculate C at start of year
C 
        STARTC=0
	  IF(ISPEAT.EQ.1)STARTD=0
        IF(ISPEAT.EQ.1)STARTR=0
	  DO 1000 IL=1,MEASLAY
	    STARTC=STARTC+DPMCARB0(IL)+RPMCARB0(IL)+BCARB0(IL)+HCARB0(IL)
	    IF(ISPEAT.EQ.1)STARTD=STARTD+DPMCARB0(IL)
	    IF(ISPEAT.EQ.1)STARTR=STARTR+RPMCARB0(IL)
1000    CONTINUE
        IF(IYEAR.EQ.1)STARTC1=STARTC
C
C Initialise soil and crop in first year
C	
        IF(IYEAR.EQ.1)THEN
C
C ...initialise the soil for initial land use,
C
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
	  ENDIF
C
C Initialise current crop
C 
        CALL INIT_GIS_CROP(IYEAR,THISLU,SECONDS,
     &                     TOTPIC(THISLU),PNREQ(THISLU),
     &                     N_STEPS,ISOWN,NSOW,IHARV,LHARV,IANTHES,
     &                     MEND,IEND,NFERT,IFERT,ILAB,FERT,TFERT,
     &                     JSTOP,CTOT)
C
C For each week till end of growing season...
C           
        DO 1100 IK=1,MEND
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
	    RAIN=AVERAIN(IK)
	    AIRTEMP=AVETEMP(IK)
	    EVAPW=AVEPET(IK)
          DO 1200 IL=1,MAXLAYER1
	      SOILTEMP(IL)=AIRTEMP
1200      CONTINUE
C
C ...Calculate crop C and N returns and N offtake 
C
	    CALL RUN1_GIS_CROP(IYEAR,IK,THISLU,IEND,MEND,NSOW,
     &                          N_STEPS,IS_TS,TOTPIC(THISLU),
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
          CALL ADD_PLANT_DIST(C_TS,RNIN,RNIN15,PI_CEQ_MON,IK,
     &                        PI_C,PI_N,PI_N15)
          IF(INMODEL.EQ.INPASSCN)
     &      CALL GET_CTON_FOEREID(DPMCTON,RPMCTON,THISLU)
	    CALL SET_DPMRPMRATIO(THISLU,DRRAT)
c	    IF(DOMSOILISIMP(ISERIES,THISLU))
c     &      CALL LIMIT_DRAIN(DOMSOILIMPDEPTH(ISERIES,THISLU),SOILW)
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
          CALL RUN_SUNDIAL_WATER_SUBMONTHLY(RAIN,WMAX,
     &                                      SOILW,WSAT,FLOWPROP,
     &                                      DRAINW,REFIL,EVAPW,
     &                                      SOILTEMP(1),ICOVER)
          CALL PHYSICAL_SUNDIAL_SOILCN(DOCMODEL,SECONDS,
     &                                   NSOIL,SOILTEMP,
     &                                   WMAX,SOILW,DRAINW,REFIL,
     &                                   SLEACH,SLEA15,
     &                                   WLEACH,WLEACH15,CONC,CONC15,
     &                                   SOILN,SOIL15,AMMN,AMMN15,
     &                                   TIN,TAM,MOBDOC,MOBDON,
     &                                   LEACHDOC,LEACHDON)
1100    CONTINUE
C
C Calculate change in soil C. If soil C is not changing, stop optimisation
C 
        ENDC=0
	  ENDB=0
	  ENDD=0
	  ENDH=0
	  ENDR=0
        FNULL=0
	  DO 1300 IL=1,MEASLAY
	    ENDC=ENDC+DPMCARB0(IL)+RPMCARB0(IL)+BCARB0(IL)+HCARB0(IL)
	    ENDD=ENDD+DPMCARB0(IL)
	    ENDR=ENDR+RPMCARB0(IL)
	    ENDB=ENDB+BCARB0(IL)
	    ENDH=ENDH+HCARB0(IL)
	    FNULL=FNULL+AMMN(IL)
1300     CONTINUE
         REPYR=500
	   IF(IYEAR.EQ.1)THEN
          WRITE(*,10)IYEAR,ENDD,ENDR,ENDB,ENDH,ENDC
	   ENDIF
         IF((IYEAR.NE.SPINCYCLE.AND.
     &       REPYR*AINT(REAL(IYEAR/REPYR)).EQ.IYEAR))THEN
	     WRITE(*,10)IYEAR,ENDD,ENDR,ENDB,ENDH,ENDC,DELTAC2
10       FORMAT(I5,6(1X,F10.0))
	   ENDIF
	   DELTAC1=DELTAC2
	   DELTAC2=ENDC-STARTC
C
C Jump out of spin-up phase if current soil conditions have been reached or passed
C
************ Temporary code to skip optimisaton
C          DELTAC2=CACCUM
************ End of temporary code to skip optimisaton
C         IF(IYEAR.GT.1)THEN        
         IF(IYEAR.GE.1)THEN        
C If increasing curve (agrading site) and change in C is less than or equal to C accumulation,...
           IF(CACCUM.GE.0.0001)THEN
		   IF(DELTAC2.LT.CACCUM+0.0001)THEN
               WRITE(*,10)IYEAR,ENDD,ENDR,ENDB,ENDH,ENDC,DELTAC2
		     GOTO 1303
	       ENDIF
C ...or if steady state and change in C is approximately 0,...
           ELSEIF(CACCUM.GT.-0.0001.AND.CACCUM.LT.0.0001)THEN
		   IF(DELTAC2.GT.-0.0001.AND.DELTAC2.LT.0.0001)THEN
               WRITE(*,10)IYEAR,ENDD,ENDR,ENDB,ENDH,ENDC,DELTAC2
  		     GOTO 1303
	       ENDIF
C ...or if decreasing curve (degrading site) and change in C is greater than (less negative than) C accumulation,...
           ELSEIF(CACCUM.LE.-0.0001)THEN
             IF(DELTAC2.GT.CACCUM+0.0001)THEN
               WRITE(*,10)IYEAR,ENDD,ENDR,ENDB,ENDH,ENDC,DELTAC2
		     GOTO 1303
	       ENDIF
	     ENDIF
C ... if reached maximum umber of iterations and still not there, stop calculation
           IF(IYEAR.EQ.SPINCYCLE)THEN
	       PRINT*,'Steady state not achieved. Press any key to continue...'
             READ(*,*)	
	       STOP
	     ENDIF
	   ENDIF	
C
C For peat, jump out after 1 year
C
         IF(ISPEAT.EQ.1)THEN
	     IF(NITER.EQ.1)THEN
	       GOTO 1303
	     ELSE
	       GOTO 1700
	     ENDIF
	   ENDIF
C
C Go back and set up the next crop
C
        IYEAR=IYEAR+1
      END DO
C
C Jump out of dynamic simulation loop
C
1303  CONTINUE
C End of ECOSSE simulation
C
C If is peat, calculate RPM and DPM using assumption that DPMC + RPMC = TOTC
C Therefore, pDPM = ((dTot/Tot)-rRPM)/(rDPM-rRPM) 
C				where	dTot is the measured change in total C (kg C / ha / yr)
C						Tot is the measred total C (kg C / ha)
C						rRPM is the rate of change of RPM C / unit C   (/yr)
C						rDPM is the rate of change of DPM C / unit C   (/yr)
C
      IF(ISPEAT.EQ.1)THEN
	  PDPM=(CACCUM/(MEASC*(1-ALPHA(1)-BETA(1))))
	  PDPM=PDPM-((ENDR-STARTR)/STARTR)
	  PDPM=PDPM/(((ENDD-STARTD)/STARTD)-((ENDR-STARTR)/STARTR))
	  PRPM=1-PDPM
	  DO 1350 IL=1,MEASLAY
	    DPMCARB0(IL)=(DPMCARB0(IL)+RPMCARB0(IL)+BCARB0(IL)+HCARB0(IL))
	    DPMNIT0(IL)=(DPMNIT0(IL)+RPMNIT0(IL)+BNIT0(IL)+HNIT0(IL))
	    DPMNLAB0(IL)=(DPMNLAB0(IL)+RPMNLAB0(IL)+BNLAB0(IL)+HNLAB0(IL))
		DPMCARB0(IL)=DPMCARB0(IL)*PDPM
		DPMNIT0(IL)=DPMNIT0(IL)*PDPM
		DPMNLAB0(IL)=DPMNLAB0(IL)*PDPM
		RPMCARB0(IL)=(1-PDPM)*DPMCARB0(IL)/PDPM
		RPMNIT0(IL)=(1-PDPM)*DPMNIT0(IL)/PDPM
		RPMNLAB0(IL)=(1-PDPM)*DPMNLAB0(IL)/PDPM
	    BCARB0(IL)=0
	    HCARB0(IL)=0
	    BNIT0(IL)=0
	    HNIT0(IL)=0
	    BNLAB0(IL)=0
	    HNLAB0(IL)=0
1350    CONTINUE
        ISAVE=1 ! Save soil 
        CALL SAVE_GIS_SOIL(ISAVE,LU1,
     &                   SOILN,SOIL15,AMMN,AMMN15,
     &                   SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
     &                   DPMCARB0,DPMNIT0,DPMNLAB0,
     &                   RPMCARB0,RPMNIT0,RPMNLAB0,
     &                   BCARB0,BNIT0,BNLAB0,
     &                   HCARB0,HNIT0,HNLAB0,
     &                   IOM,PI_CEQ_MON)      
	  GOTO 101
	ENDIF
C
C Revise plant inputs
C
c      VECTOR_OPT=0
      IF(ENDC-MEASC.GT.MEASC/100000.OR.ENDC-MEASC.LT.-MEASC/100000)THEN
C	  IF(ENDC-MEASC.LT.MEASC/10000.OR.ENDC-MEASC.GT.-MEASC/10000)THEN
c    &     VECTOR_OPT=1
        DO 1400 IL=1,MEASLAY
	    SIMC_LAST=SIMC(IL)
	    SIMC(IL)=DPMCARB0(IL)+RPMCARB0(IL)+BCARB0(IL)+HCARB0(IL)
c	    IF(VECTOR_OPT.EQ.0.OR.NITER.EQ.1)THEN  
	      PIC_LAST(IL)=STARTPIC(IL)
	      STARTPIC(IL)=STARTPIC(IL)*(TOC(IL)-IOM(IL))/SIMC(IL)
c	      STARTPIC(IL)=STARTPIC(IL)*DELTAC2/CACCUM
c	    ELSEIF(VECTOR_OPT.EQ.1.AND.NITER.GT.1)THEN
c	      LINEB=SIMC(IL)-SIMC_LAST/STARTPIC(IL)-PIC_LAST(IL)
c	      LINEA=SIMC(IL)-(LINEB*STARTPIC(IL))
c	      PIC_LAST(IL)=STARTPIC(IL)
c	      STARTPIC(IL)=(TOC(IL)-IOM(IL)-LINEA)/LINEB
c	      IF(STARTPIC(IL).LT.0)STARTPIC(IL)=0.0001/LINEB
c	    ENDIF
1400    CONTINUE
	ENDIF
C
C For a soil not at equilibrium steady state, retrieve initial soil characteristics and redo simulation with revised plant input 
C
************ Temporary code to skip optimisaton
c      ENDC=MEASC
************ End of temporary code to skip optimisaton
      IF(ENDC-MEASC.GT.MEASC/100000.OR.ENDC-MEASC.LT.-MEASC/100000)THEN
C      IF(ENDC-MEASC.GT.MEASC/10000.OR.ENDC-MEASC.LT.-MEASC/10000)THEN
        IF(CACCUM.LE.-0.0001.OR.CACCUM.GE.0.0001)THEN
          ISAVE=0 ! Retrieve soil 
          CALL SAVE_GIS_SOIL(ISAVE,LU1,
     &                   SOILN,SOIL15,AMMN,AMMN15,
     &                   SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
     &                   DPMCARB0,DPMNIT0,DPMNLAB0,
     &                   RPMCARB0,RPMNIT0,RPMNLAB0,
     &                   BCARB0,BNIT0,BNLAB0,
     &                   HCARB0,HNIT0,HNLAB0,
     &                   IOM,PI_CEQ_MON)      
	  ENDIF
c        IF(CACCUM.LE.-0.0001)THEN
c	    IF(DELTAC2.GT.0.8*CACCUM)THEN
c		  PINC=CACCUM/DELTAC2
C		  PINC=MEASC/ENDC
c	      PRINT*,' Starting soil C changed by proportion ',PINC
c	      DO 1450 IL=1,MEASLAY
c              DPMCARB0(IL)=PINC*DPMCARB0(IL)
c		    DPMNIT0(IL)=PINC*DPMNIT0(IL)
c		    DPMNLAB0(IL)=PINC*DPMNLAB0(IL)
c              RPMCARB0(IL)=PINC*RPMCARB0(IL)
c		    RPMNIT0(IL)=PINC*RPMNIT0(IL)
c		    RPMNLAB0(IL)=PINC*RPMNLAB0(IL)
c              BCARB0(IL)=PINC*BCARB0(IL)
c		    BNIT0(IL)=PINC*BNIT0(IL)
c		    BNLAB0(IL)=PINC*BNLAB0(IL)
c              HCARB0(IL)=PINC*HCARB0(IL)
c		    HNIT0(IL)=PINC*HNIT0(IL)
c		    HNLAB0(IL)=PINC*HNLAB0(IL)
c              ISAVE=1 ! Save soil 
c              CALL SAVE_GIS_SOIL(ISAVE,LU1,
c     &                   SOILN,SOIL15,AMMN,AMMN15,
c     &                   SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
c     &                   DPMCARB0,DPMNIT0,DPMNLAB0,
c     &                   RPMCARB0,RPMNIT0,RPMNLAB0,
c     &                   BCARB0,BNIT0,BNLAB0,
c     &                   HCARB0,HNIT0,HNLAB0,
c     &                   IOM,PI_CEQ_MON)      
c1450        CONTINUE
c		ENDIF
c        ENDIF
C
C ...Work out proportion of total C added each timestep
C
        PITOT=0
        DO 1500 IK=1,MEND
C
C ...Correct distribution over year to match SUNDIAL equations
C
          DO 1600 IL=1,MAXLAYER1
	      PI_CEQ_MON(IK,IL)=STARTPIC(IL)*PC_TS(IK)
	      IF(IL.GE.MEASLAY)PITOT=PITOT+STARTPIC(IL)
1600      CONTINUE
1500    CONTINUE
	  GOTO 101
C
C For degrading soil, if rate of degradation not within 80% of measured, start at a high soil C content
C
      ELSEIF(CACCUM.LE.-0.0001.AND.
     &         (DELTAC2.GT.0.8*CACCUM))THEN
C
C Retrieve initial soil characteristics and redo simulation with 10% more SOM
C
        ISAVE=0 ! Retrieve soil 
        CALL SAVE_GIS_SOIL(ISAVE,LU1,
     &                   SOILN,SOIL15,AMMN,AMMN15,
     &                   SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
     &                   DPMCARB0,DPMNIT0,DPMNLAB0,
     &                   RPMCARB0,RPMNIT0,RPMNLAB0,
     &                   BCARB0,BNIT0,BNLAB0,
     &                   HCARB0,HNIT0,HNLAB0,
     &                   IOM,PI_CEQ_MON)      
	  DO 1800 IL=1,MEASLAY
	    DPMCARB0(IL)=DPMCARB0(IL)*1.1
	    RPMCARB0(IL)=RPMCARB0(IL)*1.1
	    BCARB0(IL)=BCARB0(IL)*1.1
	    HCARB0(IL)=HCARB0(IL)*1.1
	    DPMNIT0(IL)=DPMNIT0(IL)*1.1
	    RPMNIT0(IL)=RPMNIT0(IL)*1.1
	    BNIT0(IL)=BNIT0(IL)*1.1
	    HNIT0(IL)=HNIT0(IL)*1.1
	    DPMNLAB0(IL)=DPMNLAB0(IL)*1.1
	    RPMNLAB0(IL)=RPMNLAB0(IL)*1.1
	    BNLAB0(IL)=BNLAB0(IL)*1.1
	    HNLAB0(IL)=HNLAB0(IL)*1.1
1800    CONTINUE
C
C Save current state of the soil 
C
        ISAVE=1 ! Save soil 
        CALL SAVE_GIS_SOIL(ISAVE,LU1,
     &                   SOILN,SOIL15,AMMN,AMMN15,
     &                   SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
     &                   DPMCARB0,DPMNIT0,DPMNLAB0,
     &                   RPMCARB0,RPMNIT0,RPMNLAB0,
     &                   BCARB0,BNIT0,BNLAB0,
     &                   HCARB0,HNIT0,HNLAB0,
     &                   IOM,PI_CEQ_MON)
        NITER=0
	  PRINT*,'Increasing starting soil C by 10%'
        GOTO 101   
	ENDIF 
C
C	When best fit with these pools found, save result
C
1700  CONTINUE
      ISAVE=1 ! Save soil 
      CALL SAVE_GIS_SOIL(ISAVE,LU1,
     &                   SOILN,SOIL15,AMMN,AMMN15,
     &                   SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
     &                   DPMCARB0,DPMNIT0,DPMNLAB0,
     &                   RPMCARB0,RPMNIT0,RPMNLAB0,
     &                   BCARB0,BNIT0,BNLAB0,
     &                   HCARB0,HNIT0,HNLAB0,
     &                   IOM,PI_CEQ_MON) 
      PRINT*,'ECOSSE initialisation complete'
	PRINT*,'=============================='
	PRINT*,'Dynamic simulation phase...'
C
C Leave ECOSSE_EQRUN
C
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE ECOSSE_GIS_RUN(ISWAIT)
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
      INTEGER MAXDEC				! Max.no.of decades included in calculation
	PARAMETER(MAXDEC=9)		
	INTEGER MAXDEC1				! Max.no.of decades included in calculation
	DATA MAXDEC1 /9/
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
	INTEGER MAXLU1				! Max.no.of land use types
	DATA MAXLU1 /6/
      INTEGER MAXMET				! No.of met years
	PARAMETER (MAXMET=90)
	INTEGER MAXORGM				! Max.no.of organic manure applications allowed
	PARAMETER (MAXORGM=52)
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXSEQ				! Max.no of sequences included in GIS
	PARAMETER (MAXSEQ=MAXLU*MAXLU)						
	INTEGER MAXSEQ1				! Max.no of sequences included in GIS
	DATA MAXSEQ1 /MAXSEQ/
	INTEGER MAXSOIL				! Max.no.of soil types
	PARAMETER (MAXSOIL=50)
	INTEGER MAXSOMLAY			! Max.no.of soil layers entered			
	PARAMETER (MAXSOMLAY=10)
	INTEGER MAXWC				! Max.no of wetness classes
      PARAMETER (MAXWC=6)
	INTEGER MAXWEATH			! Max.no.of years allowed
	PARAMETER (MAXWEATH=300)

	INTEGER J					! Local counter variable
	INTEGER IL					! Local counter variable
	INTEGER TESTLU1				! Initial land use to be tested
	INTEGER TESTLU2				! Land use change to be tested
	INTEGER TOTSEQ				! Total number of sequences
      CHARACTER*20 GISERROR		! Error Message for GIS
	INTEGER WATERMODEL			! IN(READ_MODEL):Water model
C
C Variables passed to / from this subroutine
C
	REAL AE						! IN: Weather data used by MAGEC
	REAL AIRTEMP				! IN: Air temperature (deg.C/timestep)	
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
      REAL ANCH4					! OUT: Measured annual CH4 emissions
								!     kgC/ha/yr
	REAL ANCO2					! OUT: Measured annual CO2 emissions
								!     kgC/ha/yr
	REAL ANDOC					! OUT: Measured annual DOC loss
								!     kgC/ha/yr
      REAL ATM					! IN:Atmospheric N input (kgN/ha/timestep)
	REAL ATM15					! IN:Atmospheric N15 input (kgN15/ha/time)
      REAL AVEPET(12)				! Long term average PET (mm/month)
      REAL AVERAIN(12)			! Long term average rainfall(mm/month)
      REAL AVETEMP(12)			! Long term average monthly average 
      REAL BALANCE(20)			! IN:C Results
	REAL BBRADS					! IN: Weather data used by MAGEC
	REAL BCARB(MAXLAYER)
	REAL BCARB0(MAXLAYER)		! IN:C in soil biomass (kgC/ha/layer)
      REAL BNIT0(MAXLAYER)		! IN:N in soil biomass (kgN/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN:N15 in soil humus (kgN15/ha/layer)
	REAL BSTART					! IN:Initial N in biomass pool (kgN/ha)
	REAL BSTART15				! IN:Initial N in biomass pool (kgN15/ha)
	REAL BULKDENS(MAXLAYER)		! IN:Bulk density of the layer (g/cm3)
	REAL C_TS					! IN:Litter C input in this timestep (kgC/ha)
      REAL CACCUM					! OUT: Measured annual C accumulation 
								!     in the soil kgC/ha/yr
	REAL CACT					! IN:Crop N uptake (kgN/ha/timestep)
	REAL CACT15					! IN:Crop N15 uptake (kgN15/ha/timestep)
	REAL CACTOT					! IN:N taken up by crop (kgN/ha)
	REAL CAL(MAXLAYER)			! IN/OUT: Concentration of Al3+ in layer (units?)
	REAL CATOT15				! IN:N15 taken up by crop (kgN15/ha)
	REAL CCHANGE10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	!IN:Change 
								! in C content over the decade (kgC/ha/decade)
	REAL CDOCIN					! IN/OUT: Concentration of mobile DOC in this layer (mg/l)
	INTEGER CELLCOUNT			! Count of cells simulated
	REAL CFACT(0:MAXCROP)		! IN:Crop parameter used to estimate yield
	REAL CH4(MAXLAYER)			! IN:CH4 emitted (kgC/ha/layer) (Aitkenhead model)
	INTEGER CH4MODEL			! OUT:Methane model (off, Richards or Aitkenhead model) 
	INTEGER CH4_OFF				! CH4 model off
	INTEGER CH4_RICHARDS    	! Richards CH4 model on
	INTEGER CH4_AITKENHEAD   	! Aitkenhead CH4 model on		
      DATA CH4_OFF,CH4_RICHARDS,CH4_AITKENHEAD /0,1,2/
      REAL CH4TOAIR				! OUT:CH4 released to atmosphere (kgC/ha) (Aitkenhead CH4 model)
	REAL CH4_PROD(MAXLAYER)  	! OUT: CH4 produced in each soil layer [kgC/ha/layer] (Richards CH4 model)
	REAL CH4_SOIL_OX			! OUT: CH4 oxidised from methane produced in the soil [kgC/ha] (Richards CH4 model)
      REAL CH4_ATMOS_OX			! OUT: CH4 oxidised from atmospheric CH4 [kgC/ha] (Richards CH4 model)
      REAL CH4_FLUX				! OUT: Net CH4 flux from soil column [kgC/ha] (Richards CH4 model)
	REAL CH4C10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! IN:CH4-C 
								! emitted over the decade (kgC/ha/decade)
      REAL CLAY(MAXLAYER)			! IN:Clay content of the layer (%)
	REAL CLOSSX					! IN:N lost by senescence (kgN/ha) 
	REAL CLOSSX15				! IN:N15 lost by senescence (kgN15/ha) 	
	REAL CNH4IN					! IN/OUT: Concentration of ammonium in this layer (eq/m3)
	INTEGER CNMODEL				! OUT:Type of Model for C:N ratio of DPM and 
	                            !     RPM (MAGEC or Foereid)
	INTEGER CNFOEREID			! C:N ratio obtained by method of Foereid
	INTEGER CNMAGEC				! C:N ratio obtained by method of MAGEC
	DATA CNMAGEC,CNFOEREID /1,2/
      REAL CNO3IN					! IN/OUT: Concentration of nitrate in this layer (eq/m3)
	REAL CO2(MAXLAYER)			! IN:CO2 emitted (kgC/ha/layer)
	REAL CO2FROMDOC(MAXLAYER)	! IN:CO2 output from DOC in each layer(kgC/ha)
	REAL CO2C10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! IN:CO2-C  
								! emitted over the decade (kgC/ha/decade)
	REAL CONC					! IN:Conc.N leached (kgN/ha/mm/layer)
	REAL CONC15					! IN:Conc.N15 leached (kgN/ha/mm/layer)
	REAL CRIT(MAXSOIL,MAXLAYER)	! IN:Critical N content (kgN/ha/layer)
	REAL CRITMIN				! IN: Data used by MAGEC
	REAL CTOT					! IN:Total litter C input (kgC/ha)
	REAL CULTDEPTH				! OUT:Cultivation depth (cm)
	REAL CUPTN					! IN:Crop N uptake (kgN/ha/timestep)
	REAL DAVTMP					! IN: Weather data used by MAGEC
	INTEGER DBRADBURY			! Bradbury model for denitrification
	REAL DDAYS					! IN: Weather data used by MAGEC
	INTEGER DECEND				! OUT:Code to save change at end of decade
	INTEGER DECSEQ				! One decade simulation only
      INTEGER DECSTART			! OUT:Code to save results at start of decade
	REAL DENITRIFN(MAXLAYER)	! IN/OUT: N lost by denitrification (kgN/ha/layer
	REAL DENIT					! IN:N lost by denitrification (kgN/ha)
	INTEGER DMODEL				! OUT:Type denitrification model chosen 
	                            ! (Bradbury or NEMIS)
	REAL DN15					! IN:N15 lost by denitrification (kgN15/ha)
	INTEGER DNEMIS				! Nemis model for denitrification
	DATA DBRADBURY,DNEMIS /1,2/
	REAL DNIT(MAXLAYER)			! IN:Net Mineralised N (kgN/ha/lay/timestep)
	REAL DNIT15(MAXLAYER)		! IN:Net Mineralised N15 (kgN15/ha/lay/time)
 	INTEGER DOC_OFF				! DOC model off
	INTEGER DOC_ON				! DOC model on
	INTEGER DOCMODEL			! OUT:DOC model (on or off)
	DATA DOC_ON,DOC_OFF /1,2/
      INTEGER DOMSOIL(CELLSTOGRID,MAXSERIES)	! 5 major soil series in the 1km 		
								! square cell x no.cells in 20km2 grid
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
      REAL DPM_RPM_RATIO			! IN/OUT:Ratio of DPM:RPM. If set to 
								!        zero, this will be worked out using a call to SET_DPMRPMRATIO
	  REAL TOC_LU(MAXLU,MAXLAYER)
      REAL IOM_LU(MAXLU,MAXLAYER) 
     	REAL DPMCARB(MAXLAYER)
	REAL DPMCARB0(MAXLAYER)		! OUT:C in decomposable PM (kgC/ha/layer)
      REAL DPMCTON(MAXLAYER)		! IN:Ratio of C:N in DPM
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL DRAINW(MAXLAYER)		! IN:Water drained from this layer (mm/layer)
	REAL DRRAT					! IN/OUT:DPM:RPM ratio
	REAL DTRJM2					! IN: Weather data used by MAGEC
	REAL DTR					! IN: Weather data used by MAGEC
	REAL DVP					! IN: Weather data used by MAGEC
	REAL EAST(CELLSTOGRID)		! National grid easting x 1000
	INTEGER EC_EQRUN			! OUT:Initialisation using a full ECOSSE equilibrium run (on or off)
	INTEGER EC_EQRUN_OFF		! Initialisation using a full ECOSSE equilibrium run is off
	INTEGER EC_EQRUN_ON			! Initialisation using a full ECOSSE equilibrium run is on
	DATA EC_EQRUN_OFF,EC_EQRUN_ON /0,1/
	INTEGER EQMODEL				! OUT:Type of equilibrium run 
	                            ! (NPP, TOC or both)
	INTEGER EQNPP				! Model initialised using plant inputs 
								!    from measured NPP
	INTEGER EQTOC				! Model initialised using measured TOC
	INTEGER EQNPPTOC			! Model initialised using both plant inputs
								!    from measured NPP and measured TOC
	INTEGER EQHILLIER,EQJONES			! Model initialised using TOC and
								!    the Hillier solver to get C pools and plant inputs
	DATA EQNPP,EQTOC,EQNPPTOC,EQHILLIER,EQJONES /1,2,3,4,5/
	REAL EVAP					! IN: Potential evap. (mm/timestep)
	REAL EVAPO					! IN: Data used by MAGEC
	REAL EVAPW					! IN: Potential evap. (mm/timestep)
								!	CONSIDER COMBINING EVAP AND EVAPW
	REAL EXYLD					! IN:Yield of current crop (t/ha)
      REAL EXYLDJ(0:MAXGROW)		! IN:Yield of current crop (t/ha)
								! (0->calculate within model)
	REAL FAREA(CELLSTOGRID)		! IN:Fraction of 20km2 in this 1km2 cell
	REAL FERT(MAXFERT)			! IN:Amount of fertiliser applied (kgN/ha)
	REAL FERTJ(0:MAXGROW,MAXFERT)	! IN:Amount of fertiliser applied (kgN/ha)
	REAL FERTNF					! INOUT:Amount of fertiliser applied (kgN/ha)
	REAL FIELDCAP(MAXLAYER)	    ! IN:Soil water content at field capacity (mm/layer)
	INTEGER FIXEND				! IN:Fixed end? 0=No 1=Yes
	CHARACTER*100 FIXFILE(MAXWEATH)	! IN: Weather files
	REAL FIXN					! IN:N fixed from atmosphere (kgN/ha)
	REAL FLOWPROP				! IN:Proportion of flow needed to achieve
							    !	     water table at depth WTABLE
	REAL FRACLU50(CELLSTOGRID,MAXLU+2) ! Fraction of 1km cell in LU1 in dec.
	INTEGER FULLSEQ				! Full sequences used (all options included)
      REAL FUTNPP(MAXMET,12,MAXLU)! IN(SETFILE_GIS):Future monthly rainfall(mm/month)
      REAL FUTRAIN(MAXMET,12)		! Future monthly rainfall(mm/month)
      REAL FUTPET(MAXMET,12)		! Future monthly PET (mm/month)
      REAL FUTTEMP(MAXMET,12)		! Future monthly average 
								!	air temp (deg.C)
	REAL FYMFERT				! IN:N input by fertiliser & FYM (kgN/ha)
	REAL FYMFERT15				! IN:N15 input by fertiliser & FYM (kgN15/ha)
	REAL G15DN2					! IN:15N2 lost by denitrification (kgN15/ha)
	REAL G15DN2O				! IN:15N2O lost by denitrification (kgN15/ha)					
	REAL G15NN2O				! IN:15N2O lost by nitrification (kgN15/ha)
	REAL G15NNO					! IN:15NO lost by nitrification (kgN15/ha)			
	REAL G15PNN2O				! IN:15N2O lost by part.nitrif.(kgN15/ha)			
	REAL GDN2					! IN:N2 lost by denitrification (kgN/ha)
	REAL GDN2O					! IN:N2O lost by denitrification (kgN/ha)
      INTEGER GIS_INDATA			! OUT(OPENCHAN_GIS): Type of data is used in this run
      INTEGER GIS_JULES			! JULES data
	INTEGER GIS_NONSSKIB		! Scottish non-SSKIB
	INTEGER GIS_SSKIB			! Scottish SSKIB
	DATA GIS_NONSSKIB,GIS_SSKIB,GIS_JULES /1,2,3/
	INTEGER GISOK				! IN:Check GIS data 0=WRONG, 1=OK
	REAL GNN2O					! IN:N2O lost by nitrification (kgN/ha)
	REAL GNNO					! IN:NO lost by nitrification (kgN/ha)
	REAL GPNN2O					! IN:N2O lost by part.nitrification (kgN/ha)
      REAL GRIDLAT				! OUT:Average latitude of this 20km2 grid cell
	REAL GRIDNPP				! Average NPP in the 1km cell (kgC/ha)
	REAL HCARB(MAXLAYER)
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN:N in soil humus (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN:N15 in soil biomass (kgN15/ha/layer)
	REAL HZ1(MAXLAYER)			! IN:N:C of input PM for steady state 
	INTEGER I_TS				! IN:No.timesteps since sowing
	INTEGER IANTHES				! IN:No.timesteps from 01/01 to anthesis 
	INTEGER IAWC				! IN:Water movement code number
	INTEGER IAWCJ				! IN:Water movement code number
  	REAL ICFACTOR(MAXLAYER)		! IN:Adjustment in rate needed to	
								!    achieve measured NPP and TOC	
	INTEGER ICFIXED				! Initialisation of C is fixed
	INTEGER ICMODEL				! OUT:Type of C initialisation model 
      INTEGER ICOVER				! OUT:Crop cover 1=Covered 0=Bare
	INTEGER ICROP				! IN:Current crop code
	INTEGER ICROPJ(0:MAXGROW)	! IN:Crop codes
	INTEGER ICROTHCEQ			! Initialisation of C by RothC 
								!  equilibrium run
	DATA ICFIXED,ICROTHCEQ /1,2/
	INTEGER IDAG					! IN: Weather data used by MAGEC
	INTEGER IDATEFC				! IN:Date of field capacity (1=01/01; 2=01/06)
      INTEGER IDEC				! OUT:Decade counter 
	INTEGER IDRAINJ				! IN:Water movement code number 
	INTEGER IEND				! IN:Timesteps from sowing date to harvest
	INTEGER IFERT(MAXFERT)		! IN:No.timesteps to split fert.application
      INTEGER IFERTJ(0:MAXGROW,MAXFERT)	! IN:No.timesteps to fert.application
	INTEGER IFILE				! IN:Current weather file
	INTEGER IHARV				! IN:Timesteps from 01/01/01 to harvest date 
	INTEGER IHARVJ(0:MAXGROW)	! IN:Timesteps from 01/01/01 to harvest date 
	INTEGER IK					! OUT:No.timesteps from prev.harv. to current
	INTEGER IL_TSS				! IN:No.timesteps before harvest when senesces
	INTEGER ILAB(MAXFERT)		! IN:Labelling on fertiliser? 0=No 1=Yes
	INTEGER ILABJ(0:MAXGROW,MAXFERT)	! IN:Lab.on fertiliser? 0=No 1=Yes
	INTEGER ILABNF				! INOUT:Labelling on fertiliser? 0=No 1=Yes
	INTEGER ILAST_TS			! IN:Last simulated timestep since 01/01/01
	INTEGER ILU					! Counter for land use 
	INTEGER IMFUNC				! OUT:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IMFUNC_HADLEY		! HADLEY moisture rate modifier
	INTEGER IMFUNC_ROTHC		! ROTHC moisture rate modifier
      DATA IMFUNC_ROTHC,IMFUNC_HADLEY /0,1/
      INTEGER IMODDEC				! Modify decomposition according to 
	                            !   full run using long term average weather data
	INTEGER IMODDEC_ON			! Decomposition modification is on
	INTEGER IMODDEC_OFF			! Decomposition modification is off
	DATA IMODDEC_OFF,IMODDEC_ON /0,1/
      INTEGER IMODPI				! Modify plant inputs according to 
	                            !   full run using long term average weather data
	INTEGER IMODPI_ON			! Plant input modification is on
	INTEGER IMODPI_OFF			! Plant input modification is off
	DATA IMODPI_OFF,IMODPI_ON /0,1/
	INTEGER INMODEL				! OUT:Type of N initialisation model
	INTEGER INPASSCN			! Initialisation of N by passing the
								!  C:N ratio of DPM and RPM
	INTEGER INSTABLECN			! Initisalisation of N by assuming 
								!  steady state (after Bradbury)
	DATA INSTABLECN,INPASSCN /1,2/
	INTEGER INSTRAW				! IN:Straw incorporated? 0=No 1=Yes
	INTEGER INVERT				! OUT(CULTIV):Invert / mix soil on cultivation? 0=No 1=Yes
	INTEGER IOLAB(MAXORGM)		! IN:Labelling on manure? 0=No 1=Yes
      INTEGER IOLABJ(0:MAXGROW,MAXORGM)	! IN:Labelling on manure? 0=No 1=Yes
	INTEGER IOLABNF				! IN:Labelling on manure? 0=No 1=Yes
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
	INTEGER IORGM(MAXORGM)		! IN:No.timesteps to split manure application
	INTEGER IORGMJ(0:MAXGROW,MAXORGM)	! IN:No.timesteps to manure applic.
	INTEGER IROCK				! IN:No.layers: 1=50cm, 2=100cm , 3=150cm
	INTEGER IROCKJ				! IN:No.layers: 1=50cm, 2=100cm , 3=150cm
	INTEGER IRYEAR				! IN:Current weather year number
	INTEGER IS_MAGEC			! Crop model type: 1=MAGEC
	INTEGER IS_SUNDIAL			! Crop model type: 0=SUNDIAL
	DATA IS_SUNDIAL,IS_MAGEC /0,1/
	INTEGER IS_TS				! IN:Crop counter
      INTEGER ISAVE				! OUT:Code to save or retrieve variables
								!		0=save; 1=retrieve
	INTEGER ISDYNPI				! Code for dynamic PI (adjusted by external factors)
	INTEGER ISDYNPI_OFF			! Adjustment of PI by external factors is off
	INTEGER ISDYNPI_ON			! Adjustment of PI by external factors is on
	DATA ISDYNPI_OFF, ISDYNPI_ON /0,1/
	INTEGER ISERIES				! Counter for soil series
      INTEGER ISERROR				! IN:Code for error in soil 1=error 0=noerror
	INTEGER ISOWN				! IN:Timesteps from 01/01/01 to sowing date 
	INTEGER ISOWNJ(0:MAXGROW)	! IN:Timesteps from 01/01/01 to sowing date 
	INTEGER ISPINUP				! OUT:Is N limitation spin-up used?
	INTEGER ISPINUP_OFF			! N limitation spin-up is not used
	INTEGER ISPINUP_ON			! N limitation spin-up is used
	DATA ISPINUP_OFF,ISPINUP_ON /0,1/
	INTEGER	ISTHARV				! IN:Timesteps from 01/01/01 to first harvest
	INTEGER ISWAIT				! IN/OUT:Code to wait for key press (1) or not (0)
      INTEGER ISTART_TS			! IN:First simulated timestep since 01/01/01 
	INTEGER ISTEST				! Output balance sheet files 
								! 0 = Do not output test files
	                            ! 1 = Output balance sheet files (BALANCE_N.OUT, BALANCE_C.OUT & SOILW.OUT)
								! 2 = Output carbon losses in t/ha 
	INTEGER ISTYR				! IN:First year in simulation (eg. 2001)
	INTEGER ITFUNC				! OUT:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER ITFUNC_ROTHC		! ROTHC temperature rate modifier
	INTEGER ITFUNC_HADLEY		! HADLEY temperature rate modifier
      DATA ITFUNC_ROTHC,ITFUNC_HADLEY /0,1/
	INTEGER IVOL(MAXFERT)		! IN:N volatilised from fertiliser (kgN/ha)
      INTEGER IVOLJ(0:MAXGROW,MAXFERT)	! IN:Does fert.contain ammonium salts 
										! other than ammonium sulphate? (0=No; 1=Yes)
	INTEGER IVOLNF				! INOUT:N volatilised from fertiliser (kgN/ha)
	INTEGER IWC					! OUT: Counter for wetness classes
	INTEGER IYEAR,temp,iiyear	! IN:Current growing season number, temp ,iiyear
								! for changing start date of the run in
	REAL JFERTJ(0:MAXGROW,MAXFERT,3)	! IN:Prop.NO3,NH4,urea in fertiliser
	INTEGER JORGM(MAXORGM)		! IN:Type of manure application
	INTEGER JORGMJ(0:MAXGROW,MAXORGM)	! IN:Type of manure application
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
	INTEGER LU1					! OUT:First land use
	INTEGER LU2					! OUT:Land use changed to
	REAL LU1TOLU2(CELLSTOGRID,MAXDEC,MAXLU+2,MAXLU+2)	! Fracn of LU1 changed 
														! to LU2 in dec
	INTEGER LUCODE				! IN/OUT:LU code for equilibrium run 
      INTEGER LUSEQ(MAXGROW)		! Land use sequence
	CHARACTER*10 KM20GRIDID		! IN:20km sqare identifier
	INTEGER MCROP				! IN:Crop type code number
      INTEGER MEASLAY				! OUT(TEST1_RES): Measured layer numbers
	INTEGER MEND				! IN:Number of timesteps from prev.harv.-harv.
	REAL MOBDOC(MAXLAYER)		! IN:Mobile DOC in each layer (kgC/ha)
	REAL MOBDON(MAXLAYER)		! IN:Mobile DON in each layer (kgN/ha)
      INTEGER MODTYPE				! OUT:Crop model type: 0=SUNDIAL, 1=MAGEC
	INTEGER N_STEPS				! IN:No.timesteps in a year
	INTEGER N_TS				! IN:No.timesteps since start of year(UNUSED?)
	REAL N2ON10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! IN:N2O-N emitted 
								! over the decade (kgN/ha/decade)
	REAL NAVTMP					! IN: Weather data used by MAGEC
	INTEGER NDATE				! IN:Not used?
	INTEGER NF					! IN:No.fertiliser application
	INTEGER NFERT				! IN:No.fertiliser applications to this crop
      INTEGER NFERTJ(0:MAXGROW)	! IN:No.fertiliser applications to crops
	INTEGER NGROW				! IN:No.growing seasons simulated
	REAL NITRIFN(MAXLAYER)		! IN/OUT:Nitrification from this layer (kg/ha/layer/timestep)
	INTEGER NORGM				! IN:No.manure applications to this crop
	INTEGER NORGMJ(0:MAXGROW)	! IN:No.manure applications to crops
	REAL NORTH(CELLSTOGRID)		! National grid northing x 1000
	REAL NRADS					! IN: Weather data used by MAGEC		
      INTEGER NRESJ(0:MAXGROW)	! IN:Straw incorporated? 0=No 1=Yes
	INTEGER NSEQ				! Counter for number of land use sequences 
	INTEGER NSOIL				! IN:Soil code number
	INTEGER NSOILJ				! IN:Soil code number
								!(CHECK USAGE OF IDRAINJ AND IAWCJ)
	INTEGER NSOMLAY(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN:Number of SOM layers
	INTEGER NSOW				! IN:Timesteps from prev.harv. to sowing date
	INTEGER NUMCELLSIN20KM2		! IN:Number of cells in 20km2 cell
	INTEGER NXYEARS				! IN:No.growing seasons simulated
	INTEGER NYEARS				! IN:No.weather years included in the sim.
	REAL ORGC					! IN:Total org.C input (kgC/ha) 
	REAL ORGMA(MAXORGM)			! IN:Amount of manure applied (t/ha)
      REAL ORGMJ(0:MAXGROW,MAXORGM)	! IN:Amount of manure applied (t/ha)
	REAL ORGMANF				! IN:Amount of manure applied (t/ha)
	REAL ORGN					! IN:Total org.N input (kgN/ha) 
	REAL PE						! IN: Weather data used by MAGEC
	REAL PENMD					! IN: Weather data used by MAGEC
	REAL PENMRS					! IN: Weather data used by MAGEC
	REAL PERSOIL(CELLSTOGRID,MAXSERIES)		! % of each major soil type in the
								!	1km square cell
	INTEGER PH_MODEL			! How is pH calculated?
	INTEGER PH_PARAM			! pH set to neutral or passed from parameter file
      INTEGER PH_STATIC			! pH is read in from input
								!   file and does not change
	INTEGER PH_FROM_VSD			! pH is calculated using VSD (a very 
								!   simple dynamic version of the MAGIC model
								!   by Ed Rowe & Chris Evans, CEH, Bangor)
	DATA PH_PARAM,PH_STATIC,PH_FROM_VSD /0,1,2/			 
      REAL PI_C(MAXLAYER)			! IN/OUT: Plant input C to soil (kgC/ha/step)
	REAL PI_CEQ_MON(12,MAXLAYER)! IN/OUT: Equilibrium plant C input each month
	REAL Equ_plant(12,MAXLAYER)	!  in each layer (kgC/ha/month/layer)
	INTEGER PI_FROMTOC			! Calculate PI from TOC during initialisation
	INTEGER PI_INPUT			! Input PI and use for initialisation
	DATA PI_INPUT,PI_FROMTOC /1,2/
      REAL PI_N(MAXLAYER)			! IN/OUT: Plant input N to soil (kgC/ha/step)
      REAL PI_N15(MAXLAYER)		! IN/OUT: Plant input N15 to soil(kgC/ha/step)
								!(1=arable,2=grassland,3=forestry,4=semi-natural)
	INTEGER PI_SOURCE			! Source of PI 1=input 2=calc from TOC
	REAL PIANN					! OUT: Total annual plant C input 
	                            !      (kg C / ha / yr) 
								!      - used as temprary value to allow adjustment of PI according 
								!        to external factors
	REAL PLANTUP(MAXLAYER)		! IN(RUN2_GIS_CROP)/OUT(TEST1_RES):Plant uptake from each soil layer
	REAL PNREQ(MAXLU)			! IN/OUT:Amount of N required by plant each yr
	                            !        (kgN/ha/year)
	REAL PREYLD					! IN:Yield of previous crop (t/ha)
	REAL OCROPN					! IN:N uptake of crop (kgN/ha) 
	REAL OCROPNJ(0:MAXGROW)		! IN:N uptake of crop (kgN/ha) 
	REAL RAIN					! IN: Rainfall (mm/timestep)
	REAL RDD					! IN: Weather data used by MAGEC
	REAL REFIL(MAXLAYER)		! IN:Water needed to fill to FC (mm/lay)
	REAL REPORTDEPTH			! IN/OUT:Depth of reporting (cm)
	REAL RLWNS					! IN: Weather data used by MAGEC
	REAL RNIN					! IN:Litter N input in this timestep (kgN/ha)
	REAL RNIN15					! IN:Litter N15 input in timestep (kgN15/ha)
	REAL ROOT					! IN:Rooting depth according to restriction (cm)
	REAL RORGN					! IN:Total org.N input (kgN/ha) DONT PASS?
	REAL RPMCARB(MAXLAYER)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
      REAL RPMCTON(MAXLAYER)		! OUT:Ratio of C:N in RPM
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL RSTART					! IN:Initial N in debris pool (kgN/ha)
	REAL RSTART15				! IN:Initial N in debris pool (kgN15/ha)
	INTEGER RUNFUTMET			! IN(OPENCHAN_GIS):Integer code for 
								!					future met.run (0=No, 1=Yes)
	REAL SAND(MAXLAYER)			! IN:Sand content of the layer (%)
	REAL SATWATCONT(MAXLAYER)	! IN:Total water content at saturation (mm/layer)
	REAL SC						! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SECONDS				! IN:Number of seconds in one timestep
	REAL SEED(MAXCROP)					! IN: Data used by MAGEC
	REAL SEEDIN					! IN: Data used by MAGEC
	REAL SEEDN_S				! IN: Data used by MAGEC
	INTEGER SEQTYPE				! Type of sequence to be used
	REAL SILT(MAXLAYER)			! IN:Silt content of the layer (%)
	REAL SLEA15(MAXLAYER)		! IN:Nitrate-N15 leached (kgN15/ha)
	REAL SLEACH(MAXLAYER)		! IN:Nitrate-N leached (kgN/ha)
	REAL SLOPES					! IN: Weather data used by MAGEC
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
      INTEGER SOILID(CELLSTOGRID*MAXSERIES/2)	! IN:Soil series integer codes
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOILPH(MAXLAYER)		! IN/OUT:pH of soil in this layer
	REAL SOILTEMP(MAXLAYER)		! IN: Soil temperature (deg.C/timestep)
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL SOMDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN:depth of soil organic matter layers
	REAL SORGC					! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SORGN					! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SORGN15				! IN:Total org.N15 input (kgN/ha) DONT PASS?
	INTEGER SPARMODEL			! OUT:Soil parameter model (from file or calc)
	INTEGER SPARFILE			! Soil parameters read in from file
	INTEGER SPARCALC			! Soil parameters calculated from TOC etc
	DATA SPARFILE,SPARCALC /1,2/
	REAL SQUID(CELLSTOGRID)		! 1km square identifier
	REAL SRNIN					! IN:Total litter N input (kgN/ha)
      logical(4) START				! OUT: Starting pH run? 0=No, 1=Yes
	INTEGER STYPES				! IN:Number of soil types in 20km2 grid
	INTEGER SUM_TS				! IN:Total number of timesteps passed
	REAL SVPS					! IN: Weather data used by MAGEC
	REAL SX						! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SXORGC					! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SXORGN					! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SXORGN15				! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL T15					! IN:Net Mineralised N (kgN/ha/layer)
	REAL TACTOT					! IN:Total crop N uptake (kgN/ha)
	REAL TAM(MAXLAYER)			! IN:15N/N in the ammonium in the layer
	REAL TAVS					! IN: Weather data used by MAGEC
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
	REAL TMAX					! IN: Weather data used by MAGEC
	REAL TMIN					! IN: Weather data used by MAGEC
	REAL TMMN					! IN: Weather data used by MAGEC
	REAL TMMX					! IN: Weather data used by MAGEC
	REAL TOC(MAXLAYER)			! IN:Measured total organic C (kgC/ha/layer)
      REAL TORGC					! IN:Total org.C input (kgC/ha) DONT PASS?
	REAL TORGN					! IN:Total org.N input (kgN/ha) DONT PASS?
	REAL TOTAST					! IN:Total ammonium in soil (kgN/ha)
	REAL TOTNST					! IN:Total nitrate in soil (kgN/ha)
	REAL TOTAST15				! IN:Total ammonium in soil (kgN15/ha)
	REAL TOTNST15				! IN:Total nitrate in soil (kgN15/ha)
	REAL TOTPIC(MAXLU)			! IN/OUT:Plant C input calculated using 
								!	  ROTHC equilibrium run for each LU (kgC/ha/yr)
      INTEGER TRAPERR(CELLSTOGRID*MAXSERIES/2,MAXSEQ) 
								! IN:Error trapping
								!		0=No error
								!		1=Error in NPP data
								!		2=Error in weather data
								!		3=Error in LU data
								!		4=Error in LU change data
								!		5=Soil error
								!		6=Equilibrium not found
	REAL TRATEM
	REAL TRNIN					! IN:Total litter N input (kgN/ha)
	REAL TXORGC					! IN:Total org.C input (kgC/ha) DONT PASS?
	REAL TXORGN					! IN:Total org.N input (kgN/ha) DONT PASS?
	REAL UT(3,0:MAXCROP)		! IN:Crop parameter used to estimate yield
	REAL VIGOUR					! OUT: Vigour of cultivation (0-1) Determines the proportion of humus released to biomass,DPM and RPM pools
	REAL VOLAT					! IN:N lost by volatilisation (kgN/ha)
	REAL VOLAT15				! IN:N15 lost by volatilisation (kgN15/ha)
	REAL VOLATN(MAXLAYER)		! IN/OUT:Volatilisation from this layer (kg/ha/layer/timestep)
	REAL VP						! IN: Weather data used by MAGEC
	REAL VSD_STATE_IL(10)		! IN: VSD State layer IL
	REAL WATERIN				! OUT: Water moving into the layer (m/timestep)
	REAL WATEROUT				! OUT: Water moving out of the layer (m/timestep)
      INTEGER WC_SOIL(CELLSTOGRID*MAXSERIES/2,MAXWC)	! OUT:Wetness classes for this soil series
	REAL WDF					! IN: Weather data used by MAGEC
      INTEGER WETCLASS(CELLSTOGRID,MAXSERIES) !OUT:Wetness class for 
								! each major soil type in the 1km2 square cell 
								! x 20km2 cell
	REAL WILTPOINT(MAXLAYER)	! IN:Water content at wilting point (mm/layer)
	REAL WLEACH					! IN:N leached (kgN/ha)
	REAL WLEACH15				! IN:N15 leached (kgN15/ha)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WN						! IN: Weather data used by MAGEC
	REAL WR						! IN:Root requirement (kgN/ha)
	REAL WRATEDM
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
	REAL WTABLE					! IN:Water table depth in cm
      INTEGER WTYPES(CELLSTOGRID*MAXSERIES/2)	! OUT:Number of wetness classes for each soil type
	REAL XORGC					! IN:Total org.C input (kgC/ha) 
	REAL XORGN					! IN:Total org.N input (kgN/ha) 
	REAL YLD					! IN:Yield (t/ha) DONT PASS?
	LOGICAL FULL_OUTPUT  	    ! Flags whether to write detailed results files
	LOGICAL SUMMARY_OUTPUT  	! Flags whether to write a summary results file
	integer imon				!counter for months for crop
	real crop(12)				!array for icover calcs cover =1 no cover =0
	real gley(MAXLU)

	DATA FULLSEQ,DECSEQ /1,2/
	DATA DECSTART,DECEND /1,2/
C
C Set restrictions on calculation
C Get parameters from TOC,IOM,CLAY,SILT,SAND,BULKDENS,	
C no restriction due to drainage,
C initialisation using ROTHC model & Bente Foereid's C:N ratios
C and get pH / water restriction on decomposition from NPP and TOC
C Use only 1 decade calculation
C
      ISTEST=3
	TESTLU1=1
	TESTLU2=1
	PI_SOURCE=PI_INPUT 
      CELLCOUNT=0
	ICMODEL=ICROTHCEQ
	ISPINUP=ISPINUP_OFF
      IMFUNC=IMFUNC_ROTHC
      ITFUNC=ITFUNC_ROTHC
	INMODEL=INPASSCN
      SEQTYPE=DECSEQ
	IF(PI_SOURCE.EQ.PI_FROMTOC)THEN
	  EQMODEL=EQTOC		! Select to optimise PI against meas.TOC when NPP is not reliable
	  IMODPI=IMODPI_ON
	  IMODDEC=IMODDEC_OFF
	ELSEIF(PI_SOURCE.EQ.PI_INPUT)THEN
	  EQMODEL=EQJONES	! Select to use PI to initialise pools and ICFACTOR when NPP is reliable
	  IMODPI=IMODPI_OFF
	  IMODDEC=IMODDEC_OFF
	ENDIF
	ISDYNPI=ISDYNPI_OFF
	CNMODEL=1 !CNFOEREID  Uses MINER1_RICHARDS
	DMODEL=DBRADBURY
	NUMCELLSIN20KM2=CELLSTOGRID
	ATM=35./52.
      MODTYPE=IS_SUNDIAL
	SPARMODEL=SPARCALC
	WTABLE=300	
      ISPINUP=ISPINUP_OFF
	INMODEL=INPASSCN
      PH_MODEL=PH_STATIC
      CH4MODEL=CH4_OFF
 	DOCMODEL=DOC_OFF 

99    CONTINUE
      PRINT*,'What type of data is used in this run?'
      PRINT*,'	1 = Scottish non-SSKIB'
      PRINT*,'	2 = Scottish SSKIB'
      PRINT*,'	3 = JULES data'
	READ(*,*)GIS_INDATA
	IF(GIS_INDATA.LT.GIS_NONSSKIB.OR.GIS_INDATA.GT.GIS_JULES)GOTO 99

      CALL TELL_MODEL(DMODEL,ICMODEL,INMODEL,CNMODEL,DOCMODEL,SPARMODEL,
     &				EQMODEL,IMFUNC,ITFUNC,CH4MODEL,ISPINUP,ISWAIT)
C
C Open Channels for file input/output
C
	CALL OPENCHAN_GIS(GIS_INDATA,REPORTDEPTH,RUNFUTMET)
      IF(ISTEST.EQ.1)THEN
	  CALL TEST1_OPENCHAN(WATERMODEL,FULL_OUTPUT, SUMMARY_OUTPUT)
	  PRINT*,'Test results output in main simulation'
	ENDIF
	IF(ISTEST.EQ.2)CALL TEST2_OPENCHAN()
	IF(ISTEST.EQ.3)CALL TEST3_OPENCHAN()
C
C Loop back to next cell
C
      PRINT*,'---------------------------------------------------------'
	PRINT*,' '
      IF(RUNFUTMET.EQ.0)THEN
	  PRINT*,'SIMULATION OF THE IMPACTS OF LAND USE CHANGE'
      ELSEIF(RUNFUTMET.EQ.1)THEN
	  PRINT*,'SIMULATION OF THE IMPACTS OF CLIMATE CHANGE'
	ENDIF
101   CONTINUE
C
C Get simulation inputs
C
c      CALL SETFILE_GIS(ISWAIT,GIS_INDATA,NSOIL,IAWC,SOMDEPTH,NSOMLAY,
c     &                 KM20GRIDID,NUMCELLSIN20KM2,
c     &                 SQUID,EAST,NORTH,DOMSOIL,PERSOIL,FAREA,
c     &                 DOMSOILC,DOMSOILCLAY,DOMSOILSILT,				
c     &                 DOMSOILSAND,DOMSOILBD,DOMSOILPH,
c     &                 DOMSOILISIMP,DOMSOILIMPDEPTH,
c     &                 STYPES,SOILID,WTYPES,WC_SOIL,WETCLASS,
c     &				 GRIDNPP,FRACLU50,LU1TOLU2,
c     &                 AVERAIN,AVEPET,AVETEMP,
c     &                 RUNFUTMET,FUTNPP,FUTRAIN,FUTPET,FUTTEMP,
c     &                 PI_SOURCE,GLEY)
	CELLCOUNT=CELLCOUNT+1
C
C Set total number of LU sequences
C
	TOTSEQ=(MAXLU1*MAXLU1)
C
C Record cell on screen
C
	WRITE(*,10)CELLCOUNT,KM20GRIDID,STYPES
10    FORMAT(
     &    /' ---------------------------------------------------------'
     &    /' GRID SQUARE ',I6,' ID ',A10/
     &     ' Number of valid soil types = ',I4/
     &     ' ---------------------------------------------------------')
C
C Initially show no error in this grid cell
C 
      CALL SETNOERROR(TRAPERR)
C
C If cell data not valid (NPP and weather), miss out completely
C
      GISOK=0
      CALL CHECK_GIS_GRID(GISOK,TRAPERR,STYPES,TOTSEQ,
     &                    GRIDNPP,AVERAIN,AVEPET,AVETEMP)
	IF(GISOK.EQ.0)THEN
	  GOTO 101
      ENDIF
C
C Run simulation for all combinations of soil, wetness class and LU type in the cell
C
      DO 300 ISERIES=1,STYPES
        WRITE(*,20)CELLCOUNT,ISERIES,SOILID(ISERIES)
20      FORMAT(I4,'.'I4,'. SERIES ',I6)
	  DO 350 IWC=1,WTYPES(ISERIES)
C
C Get the plant inputs for different LU types on this soil series
C
C
 !          IF(PI_SOURCE.EQ.PI_INPUT)THEN
 !            CALL GET_GIS_PI_FROM_ROTHC(GRIDNPP,FRACLU50,ISERIES,
 !    &                    SOMDEPTH,DOMSOILC,DOMSOILCLAY,DOMSOILSILT,				
 !    &                    DOMSOILSAND,DOMSOILBD,DOMSOILPH,
 !   &                    IAWC,IROCK,NSOIL,
 !    &                    WTABLE,ROOT,AVEPET,AVETEMP,AVERAIN,
 !    &                    ITFUNC,IMFUNC,SPARMODEL,TOTPIC,
 !    &                    NUMCELLSIN20KM2,NSOMLAY)
  !         ELSEIF(PI_SOURCE.EQ.PI_FROMTOC)THEN
!	       DO 375 ILU=1,MAXLU1
!	         CALL MIAMI(TOTPIC(ILU),AVERAIN,AVETEMP)
!375          CONTINUE		! ILU=1,MAXLU1
!           ENDIF Use EqJones instead
C
C Use soil wetness class
C
          WRITE(*,25)WC_SOIL(ISERIES,IWC)
25        FORMAT('           Wetness class ',I2)
          CALL USE_WETNESS_CLASS(WC_SOIL(ISERIES,IWC))
          
C
C Initialise soil for each land use type
C
          DO 400 LU1=1,MAXLU1
C
C ...Set time factors (timestep=monthly)
C
	      C_TS=0.0
            N_TS=0
	      SECONDS=(365.25/12)*24*60*60
C
C ...Get soil characteristics 
C
	      CALL GETSOIL(ISERIES,LU1,SOMDEPTH,NSOMLAY,							
     &                 DOMSOILC,DOMSOILCLAY,DOMSOILSILT,				
     &                 DOMSOILSAND,DOMSOILBD,DOMSOILPH,
     &                 DOMSOILISIMP,DOMSOILIMPDEPTH,
     &                 TOC,IOM,CLAY,SILT,SAND,BULKDENS,SOILPH,ISERROR)
						DO 1979 IL=1,MAXLAYER
                        TOC_LU(LU1,IL)=TOC(IL)
                        IOM_LU(LU1,IL)=IOM(IL)
1979                      continue


	      IF(ISERROR.eq.1)THEN
	        WRITE(*,30)LU1
30            FORMAT('            Exclude LU ',I1,'. Soil not defined')
	        GOTO 400
	      ENDIF
	      
C
C
C ...Get plant distribution for this land use
C
      CALL GET_PLANT_DIST(TOTPIC(LU1),PI_CEQ_MON,LU1)
      IF(GLEY(LU1).GE.0)THEN
      WTABLE=GLEY(LU1)*100
      ENDIF

C
C ...Initialise water
C              
		  CALL INIT_GIS_WATER(BULKDENS,CLAY,SAND,FLOWPROP,IAWC,IROCK,
     &                        NSOIL,SILT,TOC,WMAX,WSAT,WTABLE,ROOT,
     &                        AVEPET,AVERAIN,AVETEMP,SOILW,PI_CEQ_MON,
     &                        WILTPOINT,FIELDCAP,SATWATCONT,SECONDS)
C
C ...Get weather data
C
	DO IMON=1,12
	IF(PI_CEQ_MON(IMON,LU1).GT.0)THEN
	CROP(IMON)=1
	ELSE
	CROP(IMON)=0
	ENDIF
	ENDDO
	  CALL GETLTA_LIM(LTA_AWC,LTA_TEMP,WMAX,WSAT,FLOWPROP,  !Gets long term average water, temp etc. uses crop to denote cover,  
     &                AVERAIN,AVEPET,AVETEMP,CROP)            !thus matches water in spinup and the runs. 
C
C ...Get soil C and N parameters
C
	      CALL GET_SUNDIAL_SOILPARS(NSOIL,TOC,IOM,CRIT,LU1,
     &                              CLAY,BULKDENS,SPARMODEL)

      IF(PI_SOURCE.EQ.PI_INPUT)THEN
   	        EQMODEL=EQNPPTOC
			EQMODEL=EQJONES
      Do 660 THISLU=1,6
          CALL MIAMI_DYCE(THISLU,sum(AVETEMP)/12,
     &    sum(AVERAIN),PIANN)
        TOTPIC(THISLU)=PIANN
 660   continue
			
            ELSEIF(PI_SOURCE.EQ.PI_FROMTOC)THEN
	        EQMODEL=EQTOC
            EQMODEL=EQHILLIER
	      ENDIF
		  
C
C ...Get plant distribution for this land use
C
      CALL GET_PLANT_DIST(TOTPIC(LU1),PI_CEQ_MON,LU1)

C
C ...Initialise the soil for this land use
C
            IF(INMODEL.EQ.INPASSCN)
     &        CALL GET_CTON_FOEREID(DPMCTON,RPMCTON,LU1)
	      IF(ICMODEL.EQ.ICROTHCEQ.AND.EQMODEL.NE.EQJONES)
     &        CALL SET_DPMRPMRATIO(LU1,DPM_RPM_RATIO)
            
             IF(EQMODEL.NE.EQJONES)then
                 iyear=1
              CALL INIT_GIS_CROP(IYEAR,LU1,SECONDS,
     &                     TOTPIC(LU1),PNREQ(LU1),
     &                     N_STEPS,ISOWN,NSOW,IHARV,LHARV,IANTHES,
     &                     MEND,IEND,NFERT,IFERT,ILAB,FERT,TFERT,
     &                     JSTOP,CTOT)
              
       do 1999 IK=1,12
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
1999          continue    
        
		CALL INIT_GIS_SOILCN_NOPARS(ICMODEL,INMODEL,DOCMODEL,
     &                                EQMODEL,SECONDS,
     &                                NSOIL,TOC,IOM,HZ1,CRIT,
     &                                SOILN,SOIL15,AMMN,AMMN15,
     &                                DPMCARB0,DPMNIT0,DPMNLAB0,
     &                                RPMCARB0,RPMNIT0,RPMNLAB0,
     &                                BCARB0,BNIT0,BNLAB0,
     &                                HCARB0,HNIT0,HNLAB0,
     &                                MOBDOC,MOBDON,SATWATCONT,
     &                                CLAY,BULKDENS,Equ_plant,
     &                                LTA_AWC,LTA_TEMP,ITFUNC,IMFUNC,
     &                                WMAX,WSAT,ICFACTOR,wiltpoint,
     &                                DPM_RPM_RATIO,DPMCTON,RPMCTON,
     &                                PH_MODEL,SOILPH,
     &                                CACCUM,ANCH4,ANCO2,ANDOC)
             endif
             PI_CEQ_MON=Equ_plant
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
	          CALL GETEQUIVALENTS(SOILN(IL),AMMN(IL),MOBDOC(IL),
     &                       SOILW(IL),CNO3IN,CNH4IN,CDOCIN)
                IF(IL.EQ.1)THEN
	            WATERIN=AVERAIN(1)/1000.
	            WATEROUT=DRAINW(IL)/1000.
	          ELSEIF(IL.GT.1)THEN
	            WATERIN=DRAINW(IL-1)/1000.
	            WATEROUT=DRAINW(IL)/1000.
	          ENDIF
	          THICK=MAXDEPTH/(MAXLAYER*100.)
                CALL RUN_VSD(START,CNO3IN,CNH4IN,CDOCIN,
     &               WATERIN,WATEROUT,THICK,
     &               SOILPH(IL),CAL(IL),ISERIES,VSD_STATE_IL)
	          ISAVE=1
	          CALL SAVE_VSD_STATE(ISAVE,VSD_STATE_IL,IL)
200           CONTINUE		! IL=1,MAXLAYER1
            ENDIF
C
C ...Get rate modifier associated with SUNDIAL routines,
C
            ISPINUP=ISPINUP_OFF
            IF(ISPINUP.EQ.ISPINUP_ON.AND.ICMODEL.EQ.ICROTHCEQ)THEN
              CALL GET_NLIM(AMMN,AMMN15,ATM,AVEPET,AVERAIN,AVETEMP,
     &                BCARB0,BNIT0,BNLAB0,CH4MODEL,CNMODEL,CTOT,
     &                CLOSSX,CLOSSX15,CRIT,DMODEL,DOCMODEL,DPMCARB0,
     &			    DPMNIT0,DPMNLAB0,DRAINW,FERT,FIXN,FLOWPROP,
     &                FYMFERT,FYMFERT15,HCARB0,HNIT0,HNLAB0,IANTHES,
     &                IEND,IFERT,IHARV,ILAB,IMFUNC,ICMODEL,INMODEL,
     &                IOLAB,IOM,IORGM,ISERIES,ISOWN,ITFUNC,IVOL,
     &                JORGM,JSTOP,LHARV,LU1,MEND,NFERT,NORGM,N_STEPS,
     &                NSOIL,NSOW,ORGMA,PNREQ,REFIL,
     &                RPMCARB0,RPMNIT0,RPMNLAB0,
     &                SECONDS,SOIL15,SOILN,SOILW,TAM,TFERT,THISFERT,
     &                THISFERT15,TIN,TOC,TOTPIC,VOLAT,VOLAT15,WLEACH,
     &                WLEACH15,WSAT,WMAX,SOILPH,PI_CEQ_MON)		  
            ENDIF
C
C ......get modification to plant inputs due to full run
C
            IMODPI=IMODPI_OFF
            IF(IMODPI.EQ.IMODPI_ON)THEN
              NXYEARS=(MAXDEC1*10)
	        CALL GET_MODPI_LAYERED(AMMN,AMMN15,ATM,AVEPET,AVERAIN,
     &              AVETEMP,BCARB0,BNIT0,BNLAB0,CH4MODEL,CNMODEL,CTOT,
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
            EC_EQRUN=EC_EQRUN_OFF


            IF(EC_EQRUN.EQ.EC_EQRUN_ON)THEN
c	    CALL ECOSSE_EQRUN(AMMN,AMMN15,ATM,AVEPET,AVERAIN,AVETEMP,
c     &              BCARB0,BNIT0,BNLAB0,CACCUM,CH4MODEL,CNMODEL,CTOT,
c     &			  CLOSSX,CLOSSX15,CRIT,DMODEL,DOCMODEL,DPMCARB0,
c     &			  DPMNIT0,DPMNLAB0,DRAINW,FERT,FIXN,FLOWPROP,
c     &              FYMFERT,FYMFERT15,HCARB0,HNIT0,HNLAB0,
c     &              IANTHES,IEND,IFERT,IHARV,ILAB,IMFUNC,ICMODEL,
c     &              INMODEL,IOLAB,IOM,IORGM,ISERIES,ISOWN,ITFUNC,IVOL,
c     &              JORGM,JSTOP,LHARV,LU1,MEND,NFERT,NORGM,
c     &              NSOIL,NSOW,ORGMA,PNREQ,REFIL,RPMCARB0,
c     &              RPMNIT0,RPMNLAB0,SOIL15,SOILN,SOILW,TAM,
c     &              TFERT,THISFERT,THISFERT15,TIN,TOC,TOTPIC,VOLAT,
c     &              VOLAT15,WLEACH,WLEACH15,WSAT,WMAX,SOILPH,NXYEARS,
c     &              EQMODEL,CLAY,PI_CEQ_MON,DOMSOILISIMP,
c     &              DOMSOILIMPDEPTH,SOMDEPTH,NSOMLAY,LTA_AWC,LTA_TEMP,
c     &              PH_MODEL)
	      ENDIF
C
C ......get modification to decomposition due to full run
C
            IF(IMODDEC.EQ.IMODDEC_ON)THEN
	        CALL GET_MODDEC_LAYERED(AMMN,AMMN15,ATM,AVEPET,AVERAIN,
     &              AVETEMP,BCARB0,BNIT0,BNLAB0,CH4MODEL,CNMODEL,CTOT,
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
     &                         SOILN,SOIL15,AMMN,AMMN15,
     &                        SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
     &                        DPMCARB0,DPMNIT0,DPMNLAB0,
     &                        RPMCARB0,RPMNIT0,RPMNLAB0,
     &                        BCARB0,BNIT0,BNLAB0,
     &                        HCARB0,HNIT0,HNLAB0,
     &                        IOM,PI_CEQ_MON)
      
400       CONTINUE	! LU1=1,MAXLU1
C
C For each possible sequence...
C
          DO 500 NSEQ=1,TOTSEQ
            GPNN2O=0
            GDN2O=0
            GISERROR=' '
C
C Get land use sequence
C
            CALL GET_LUSEQ(LUSEQ,NGROW,NSEQ,SEQTYPE)
C
C For simulation using long term average met data...
C
            IF(RUNFUTMET.EQ.0)THEN
C
C ...If this sequence not occuring on this soil, miss out this sequence
C
	        CALL CHECK_GIS_LU(GISOK,LUSEQ,LU1TOLU2,
     &                      PERSOIL,DOMSOIL,SOILID,
     &                      NUMCELLSIN20KM2,TRAPERR,ISERIES,NSEQ)
C
C For simulation using future met data...
C
            ELSEIF(RUNFUTMET.EQ.1)THEN
C
C ...only run for no land use change (initially)
C
              IF(LUSEQ(1).EQ.LUSEQ(6))THEN
	          GISOK=1
 	        ELSE
	          GISOK=0
	        ENDIF
	      ENDIF
C
C Check soil data for this sequence
C
c		  IF(GISOK)CALL GETSOIL(ISERIES,LUSEQ(1),SOMDEPTH,NSOMLAY,							
c     &                 DOMSOILC,DOMSOILCLAY,DOMSOILSILT,				
c     &                 DOMSOILSAND,DOMSOILBD,DOMSOILPH,
c     &                 DOMSOILISIMP,DOMSOILIMPDEPTH,
c     &                 TOC,IOM,CLAY,SILT,SAND,BULKDENS,SOILPH,ISERROR)
c	      IF(ISERROR)GISOK=0
c		  IF(GISOK)CALL GETSOIL(ISERIES,LUSEQ(6),SOMDEPTH,NSOMLAY,							
c     &                 DOMSOILC,DOMSOILCLAY,DOMSOILSILT,				
c     &                 DOMSOILSAND,DOMSOILBD,DOMSOILPH,
c     &                 DOMSOILISIMP,DOMSOILIMPDEPTH,
c     &                 TOC,IOM,CLAY,SILT,SAND,BULKDENS,SOILPH,ISERROR)
c	      IF(ISERROR)GISOK=0
	      IF(GISOK.EQ.0)THEN
              CALL SETERROR(FYMFERT,FYMFERT15,TFYMC,VOLAT,VOLAT15,					  
     &                    SOILN,SOIL15,AMMN,AMMN15,			  
     &                    DPMCARB0,DPMNIT0,DPMNLAB0,		  
     &                    RPMCARB0,RPMNIT0,RPMNLAB0,		  
     &                    HCARB0,HNIT0,HNLAB0,TIN,TAM,
     &                    CCHANGE10,CO2C10,CH4C10,N2ON10,
     &                    ISERIES,NSEQ,WC_SOIL(ISERIES,IWC))
	        GOTO 202
            ENDIF	     
C
C Set time factors (timestep=monthly)
C
	      C_TS=0.0
            N_TS=0
	      SECONDS=(365.25/12)*24*60*60
C
C For each growing season... 
C        
            NXYEARS=(MAXDEC1*10)
            DO 600 iiyear=1,NXYEARS
	        IYEAR=iiyear
	        LU1=LUSEQ(1)
	        LU2=LUSEQ(6)
	        IF(IYEAR.GT.1)THEN !NGROW  Changes to make LU from year 1
	          THISLU=LU2
	        ELSE
	          THISLU=LUSEQ(IYEAR)
	        ENDIF


	CALL GET_PLANT_DIST(TOTPIC(THISLU),PI_CEQ_MON,THISLU)
C
C Initialise soil and crop in first year
C	
              IF(IYEAR.EQ.1)THEN
			  IF(GLEY(LU1).GE.0)THEN
                WTABLE=GLEY(LU1)*100
              ELSE
                  WTABLE=300
         
              ENDIF
			   CALL INIT_GIS_WATER(BULKDENS,CLAY,SAND,FLOWPROP,IAWC,IROCK,
     &                        NSOIL,SILT,TOC,WMAX,WSAT,WTABLE,ROOT,
     &                        AVEPET,AVERAIN,AVETEMP,SOILW,PI_CEQ_MON,
     &                        WILTPOINT,FIELDCAP,SATWATCONT,SECONDS)
			IF(EQMODEL.EQ.EQJONES)THEN       
            CALL SET_DPMRPMRATIO(LU1,DPM_RPM_RATIO)
            
	      CALL GETSOIL(ISERIES,LU1,SOMDEPTH,NSOMLAY,							
     &                 DOMSOILC,DOMSOILCLAY,DOMSOILSILT,				
     &                 DOMSOILSAND,DOMSOILBD,DOMSOILPH,
     &                 DOMSOILISIMP,DOMSOILIMPDEPTH,
     &                 TOC,IOM,CLAY,SILT,SAND,BULKDENS,SOILPH,ISERROR)
			DO 899 il=1,maxlayer
            TOC(il)=TOC_LU(LU1,il)
            IOM(il)=IOM_LU(LU1,il)
899           continue
        	   CALL INIT_GIS_WATER(BULKDENS,CLAY,SAND,FLOWPROP,IAWC,IROCK,
     &                        NSOIL,SILT,TOC,WMAX,WSAT,WTABLE,ROOT,
     &                        AVEPET,AVERAIN,AVETEMP,SOILW,PI_CEQ_MON,
     &                        WILTPOINT,FIELDCAP,SATWATCONT,SECONDS)
              
            CALL INIT_GIS_CROP(IYEAR,THISLU,SECONDS,
     &            TOTPIC(THISLU),PNREQ(THISLU),
     &            N_STEPS,ISOWN,NSOW,IHARV,LHARV,IANTHES,
     &            MEND,IEND,NFERT,IFERT,ILAB,FERT,TFERT,
     &            JSTOP,CTOT)
              do 999 IK=1,12
             CALL RUN1_GIS_CROP(IYEAR,IK,THISLU,IEND,MEND,NSOW,
     &                          N_STEPS,IS_TS,TOTPIC(THISLU),
     &                          PNREQ(THISLU),CACTOT,CATOT15,
     &                          ICOVER,SXORGN,ORGN,
     &                          TRNIN,RNIN,RNIN15,C_TS,TCINP,WR)   
           if(sum(PI_CEQ_MON(ik,:)).GT.0)then
          Equ_plant(ik,:)= C_TS*PI_CEQ_MON(ik,:)/sum(PI_CEQ_MON(ik,:))      
           else
               Equ_plant(ik,:)=0
               endif
999           continue
      CALL GETLTA_LIM(LTA_AWC,LTA_TEMP,WMAX,WSAT,FLOWPROP,  !Gets long term average water, temp etc. uses crop to denote cover,  
     &          AVERAIN,AVEPET,AVETEMP,CROP)
            CALL INIT_GIS_SOILCN_NOPARS(ICMODEL,INMODEL,DOCMODEL,
     &                                EQMODEL,SECONDS,
     &                                NSOIL,TOC,IOM,HZ1,CRIT,
     &                                SOILN,SOIL15,AMMN,AMMN15,
     &                                DPMCARB0,DPMNIT0,DPMNLAB0,
     &                                RPMCARB0,RPMNIT0,RPMNLAB0,
     &                                BCARB0,BNIT0,BNLAB0,
     &                                HCARB0,HNIT0,HNLAB0,
     &                                MOBDOC,MOBDON,SATWATCONT,
     &                                CLAY,BULKDENS,Equ_plant,
     &                                LTA_AWC,LTA_TEMP,ITFUNC,IMFUNC,
     &                                WMAX,WSAT,ICFACTOR,wiltpoint,
     &                                DPM_RPM_RATIO,DPMCTON,RPMCTON,
     &                                PH_MODEL,SOILPH,
     &                                CACCUM,ANCH4,ANCO2,ANDOC)
            ENDIF
C
C ...Get soil characteristics 
C
	          CALL GETSOIL(ISERIES,LU1,SOMDEPTH,NSOMLAY,							
     &                 DOMSOILC,DOMSOILCLAY,DOMSOILSILT,				
     &                 DOMSOILSAND,DOMSOILBD,DOMSOILPH,
     &                 DOMSOILISIMP,DOMSOILIMPDEPTH,
     &                 TOC,IOM,CLAY,SILT,SAND,BULKDENS,SOILPH,ISERROR)
c	          IF(ISERROR)THEN
c                  CALL SETERROR(FYMFERT,FYMFERT15,TFYMC,VOLAT,VOLAT15,					  
c     &                      SOILN,SOIL15,AMMN,AMMN15,			  
c     &                      DPMCARB0,DPMNIT0,DPMNLAB0,		  
c     &                      RPMCARB0,RPMNIT0,RPMNLAB0,		  
c     &                      HCARB0,HNIT0,HNLAB0,TIN,TAM,
c     &                      CCHANGE10,CO2C10,CH4C10,N2ON10,
c     &                      ISERIES,NSEQ,WC_SOIL(ISERIES,IWC))
c	            GOTO 202
c	          ENDIF
                ISAVE=0 ! Retrieve soil 
				if(EQMODEL.EQ.EQJONES)
     & ISAVE=1    

                CALL SAVE_GIS_SOIL(ISAVE,LU1,
     &                        SOILN,SOIL15,AMMN,AMMN15,
     &                        SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
     &                        DPMCARB0,DPMNIT0,DPMNLAB0,
     &                        RPMCARB0,RPMNIT0,RPMNLAB0,
     &                        BCARB0,BNIT0,BNLAB0,
     &                        HCARB0,HNIT0,HNLAB0,
     &                        IOM,PI_CEQ_MON)
C
C Record sequence
C
	          IF(SEQTYPE.EQ.FULLSEQ)THEN
                  WRITE(*,40)NSEQ,LUSEQ(1),LUSEQ(6),LUSEQ(16),
     &                 LUSEQ(26),GISERROR
                ELSE
                  WRITE(*,50)NSEQ,LUSEQ(1),LUSEQ(6),GISERROR
	          ENDIF
40              FORMAT('           Seq.',I3,': ',4(I2,1X),A20)
50              FORMAT('           Seq.',I3,': ',2(I2,1X),A20)
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
     &                             MOBDOC,MOBDON,LU1,
     &                             CLAY,BULKDENS,PI_CEQ_MON,
     &                             LTA_AWC,LTA_TEMP,ITFUNC,IMFUNC,
     &                             WMAX,WSAT,ICFACTOR,PH_MODEL,SOILPH)
                IF(ISPINUP.EQ.ISPINUP_ON)THEN
		        CALL SET_NLIM(LU1,ICFACTOR,
     &                           SOILN,SOIL15,AMMN,AMMN15,SOILW,
     &                           DPMCARB0,DPMNIT0,DPMNLAB0,
     &                           RPMCARB0,RPMNIT0,RPMNLAB0,
     &                           BCARB0,BNIT0,BNLAB0,
     &                           HCARB0,HNIT0,HNLAB0,IOM,PI_CEQ_MON)
	          ENDIF
C
C Check the initialised soil data - if it is not correct miss out this sequence
C
                CALL CHECK_GIS_SOIL(GISOK,SOILN,SOIL15,AMMN,AMMN15,
     &			  SOILW,DPMCARB0,DPMNIT0,DPMNLAB0,RPMCARB0,RPMNIT0,
     &			  RPMNLAB0,BCARB0,BNIT0,BNLAB0,HCARB0,HNIT0,HNLAB0,
     &              TRAPERR,ISERIES,NSEQ)
	          IF(GISOK.EQ.0)THEN
	            PRINT*,'SOIL NOT PROPERLY DEFINED'
                  CALL SETERROR(FYMFERT,FYMFERT15,TFYMC,VOLAT,VOLAT15,					  
     &                      SOILN,SOIL15,AMMN,AMMN15,			  
     &                      DPMCARB0,DPMNIT0,DPMNLAB0,		  
     &                      RPMCARB0,RPMNIT0,RPMNLAB0,		  
     &                      HCARB0,HNIT0,HNLAB0,TIN,TAM,
     &                      CCHANGE10,CO2C10,CH4C10,N2ON10,
     &                      ISERIES,NSEQ,WC_SOIL(ISERIES,IWC))
	            GOTO 202
                ENDIF		  
C
C Set up results arrays
C
                CALL STARTRES(IRYEAR,IK,RSTART,DPMNIT0,RPMNIT0,
     &                    BSTART,BNIT0,HNIT0,
     &                    RSTART15,DPMNLAB0,RPMNLAB0,
     &                    BSTART15,BNLAB0,HNLAB0,
     &                    TOTAST,AMMN,TOTAST15,AMMN15,
     &                    TOTNST,SOILN,TOTNST15,SOIL15)
	          CALL INIT_GIS_RES(CO2C10,CH4C10,N2ON10,ISERIES,NSEQ,
     $                            WC_SOIL(ISERIES,IWC))
              ENDIF
C
C Initialise current crop
C 
              CALL INIT_GIS_CROP(IYEAR,THISLU,SECONDS,
     &                        TOTPIC(THISLU),PNREQ(THISLU),
     &                        N_STEPS,ISOWN,NSOW,IHARV,LHARV,IANTHES,
     &                        MEND,IEND,NFERT,IFERT,ILAB,FERT,TFERT,
     &                        JSTOP,CTOT)
              IF(ISPINUP.EQ.ISPINUP_ON)THEN
	          CALL SET_NLIM_RATE(THISLU,ICFACTOR)
	        ENDIF
C
C Cultivate soil for gra->ara, for->ara, for->gra, nat->ara, for->gra, nat->for
C
              temp=iyear
	        IF(IYEAR.eq.2)THEN
                  
                  iyear=6
	          IF((LUSEQ(IYEAR).NE.LUSEQ(IYEAR-1)))THEN
	            CULTDEPTH=0
	            VIGOUR=0
	            IF(
     &	           (LUSEQ(IYEAR-1).EQ.2.AND.LUSEQ(IYEAR).EQ.1).OR.
     &               (LUSEQ(IYEAR-1).EQ.3.AND.LUSEQ(IYEAR).EQ.1).OR.
     &               (LUSEQ(IYEAR-1).EQ.4.AND.LUSEQ(IYEAR).EQ.1))THEN
	             CULTDEPTH=50
	             VIGOUR=0.5
	             INVERT=1
                  ELSEIF(
     &               (LUSEQ(IYEAR-1).EQ.3.AND.LUSEQ(IYEAR).EQ.2).OR.
     &               (LUSEQ(IYEAR-1).EQ.3.AND.LUSEQ(IYEAR).EQ.4))THEN
	             CULTDEPTH=30
	             VIGOUR=0.5
	             INVERT=1
               ELSEIF(
     &               (LUSEQ(IYEAR-1).EQ.2.AND.LUSEQ(IYEAR).EQ.3))Then
                        IF(WTABLE.LE.50)THEN 
           PRINT*,'Cultivation to forestry',
     & ' lowering the water table from', wtable,'cm to 50cm'  
                       WTABLE=50
            CALL INIT_GIS_WATER(BULKDENS,CLAY,SAND,FLOWPROP,IAWC,IROCK,
     &                        NSOIL,SILT,TOC,WMAX,WSAT,WTABLE,ROOT,
     &                        AVEPET,AVERAIN,AVETEMP,SOILW,PI_CEQ_MON,
     &                        WILTPOINT,FIELDCAP,SATWATCONT,SECONDS)
                          
                   ENDIF
                  
	            ELSEIF(
     &               (LUSEQ(IYEAR-1).EQ.4.AND.LUSEQ(IYEAR).EQ.3))THEN
	             CULTDEPTH=30
	             VIGOUR=0.1
	             INVERT=1
			IF(WTABLE.LE.50)THEN 
      PRINT*,'Cultivation to forestry',
     & ' lowering the water table from', wtable,'cm to 50cm'  
                       WTABLE=50
            CALL INIT_GIS_WATER(BULKDENS,CLAY,SAND,FLOWPROP,IAWC,IROCK,
     &                        NSOIL,SILT,TOC,WMAX,WSAT,WTABLE,ROOT,
     &                        AVEPET,AVERAIN,AVETEMP,SOILW,PI_CEQ_MON,
     &                        WILTPOINT,FIELDCAP,SATWATCONT,SECONDS)
                          
                   ENDIF
     
	            ENDIF
      	        CALL CULTIV(LUSEQ(IYEAR-1),
     &                        DPMCARB0,DPMNIT0,DPMNLAB0,
     &                        RPMCARB0,RPMNIT0,RPMNLAB0,
     &                        BCARB0,BNIT0,BNLAB0,
     &                        HCARB0,HNIT0,HNLAB0,
     &                        CULTDEPTH,VIGOUR,INVERT)
			  ENDIF	         
	        ENDIF
C
C For each timestep till end of growing season...
          iyear=temp

      PIANN=TOTPIC(THISLU)
          
C           
       DO 700 IK=1,12 !Mend, runs till end of the year (12 months)
C
C Initialize THISFERT=current weeks fertilizer addition
C            WLEACH=current weeks water leaching
C      
                CALL SETWEEKVARS(JSTOP,CLOSSX,CLOSSX15,FYMFERT,
     &                           FYMFERT15,FIXN,THISFERT,THISFERT15,
     &                           WLEACH,WLEACH15,VOLAT,VOLAT15)
C
C Get weather data and adjust PI
C
			  
                IF(RUNFUTMET.EQ.0)THEN
                  CALL GETWEATHER_GIS(IK,LHARV,AVERAIN,AVEPET,AVETEMP, 
     &                              RAIN,EVAPW,SOILTEMP,AIRTEMP)
	          ELSEIF(RUNFUTMET.EQ.1)THEN
	            IF(GIS_INDATA.EQ.GIS_JULES)THEN
	              CALL GETFUTWEATHER_JULES(IK,THISLU,IYEAR,LHARV,
     &                              FUTNPP,FUTRAIN,FUTPET,FUTTEMP,
     &                              GRIDNPP,RAIN,EVAPW,SOILTEMP,AIRTEMP)
                    CALL ADJUST_PI_BYNPP(FUTNPP,IK,IYEAR,THISLU,PIANN)
				ELSE
	              CALL GETFUTWEATHER(IK,IYEAR,LHARV,
     &                               FUTRAIN,FUTPET,FUTTEMP,
     &                               RAIN,EVAPW,SOILTEMP,AIRTEMP)
				  IF(ISDYNPI.EQ.ISDYNPI_ON)
     &			    CALL ADJUST_PI(AVERAIN,AVETEMP,IK,LHARV,PIANN,
     &							   RAIN,SOILTEMP)
				ENDIF
	          ENDIF
C
C If current week is first in the growing season set results
C
                CALL SETRES(BALANCE,CO2,SX,SXORGN)
               

C
C ...Calculate crop C and N returns and N offtake 
C
	          CALL RUN1_GIS_CROP(IYEAR,IK,THISLU,IEND,MEND,NSOW,
     &                            N_STEPS,IS_TS,PIANN,
     &                            PNREQ(THISLU),CACTOT,CATOT15,
     &                            ICOVER,SXORGN,ORGN,
     &                            TRNIN,RNIN,RNIN15,C_TS,TCINP,WR)   
                
       IF(C_TS.GT.0)Then
              ICOVER=1
          else
              icover=0
       endif
                
             CALL RUN_SUNDIAL_WATER_SUBMONTHLY(RAIN,WMAX,
     &                                         SOILW,WSAT,FLOWPROP,
     &                                         DRAINW,REFIL,EVAPW,
     &                                         SOILTEMP(1),ICOVER)
C
C ...Add stuff to soil
C
                CALL GETFYM(IK,NORGM,IORGM,ORGMA,JORGM,IOLAB,
     &                  ORGMANF,JORGMNF,IOLABNF) 
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
                CALL ADD_PLANT_DIST(C_TS,RNIN,RNIN15,PI_CEQ_MON,IK,
     &                        PI_C,PI_N,PI_N15) 
                IF(INMODEL.EQ.INPASSCN)
     &            CALL GET_CTON_FOEREID(DPMCTON,RPMCTON,THISLU)
	          CALL SET_DPMRPMRATIO(THISLU,DRRAT)
c	          IF(DOMSOILISIMP(ISERIES,THISLU))
c     &           CALL LIMIT_DRAIN(DOMSOILIMPDEPTH(ISERIES,THISLU),SOILW)
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
                  DO 800 IL=1,MAXLAYER1
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
800               CONTINUE	! IL=1,MAXLAYER1
                ENDIF
C
C Output results 
C
       IF(IK.EQ.12)THEN
	            IDEC=INT(IYEAR/10)+1
	            IF(IYEAR.EQ.(10*(IDEC-1))+1)THEN
	              CALL GET_GIS_CHANGE(REPORTDEPTH,ISERIES,NSEQ,
     &                            WC_SOIL(ISERIES,IWC),
     &                            CCHANGE10,CO2C10,CH4C10,LU1,LU2,IDEC,
     &                            DPMCARB0,RPMCARB0,BCARB0,HCARB0)
			    ENDIF
     	          ENDIF
	          IF(LU1.EQ.TESTLU1.AND.LU2.EQ.TESTLU2)THEN ! Change equator 
						! to LU1 and LU2 to output different LU change results	
      	        IF(ISTEST.EQ.1.AND.IYEAR.GT.0)
     &              CALL TEST1_RES(IYEAR,IRYEAR,NXYEARS,
     &                   LHARV,SECONDS,IK,N_TS,SUM_TS,FIXEND,NSOW,ISOWN,
     &                   MEND,SX,SXORGN15,SORGN,SORGN15,RNIN15,
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
	          ENDIF
C
C For now assume N2O emissions are 10% of total denitrification
C
c                IF(DMODEL.EQ.DBRADBURY)THEN
c	            GNN2O=0
c                  GPNN2O=0
c                  GDN2O=0.005*DENIT
c	          ENDIF
                IDEC=INT(IYEAR/10)+1
     	          CALL GET_GIS_EMISSIONS(REPORTDEPTH,CH4,GNN2O,GPNN2O,
     &                                 GDN2O,CH4C10,N2ON10,
     &                                 ISERIES,NSEQ,IDEC,
     &                                 WC_SOIL(ISERIES,IWC))
			  IF(ISTEST.EQ.3)THEN
	            CALL TEST3_RES(KM20GRIDID,IK,IYEAR,LU1,
     &                           SLEACH,LEACHDOC,DRAINW)
	          ENDIF
700           CONTINUE		! IK=1,MEND
C
C Go back and set up the next crop
C
600         CONTINUE			! IYEAR=1,NXYEARS
	      IF(ISTEST.EQ.2)CALL TEST2_RES(KM20GRIDID,SOILID,NSEQ,
     &                                  WC_SOIL(ISERIES,IWC),
     &                                  SOMDEPTH,ISERIES,CCHANGE10,
     &                                  CO2C10,CH4C10,N2ON10,
     &                                  LU1TOLU2,
     &                                  DOMSOILC,DOMSOILCLAY,DOMSOILBD,
     &                                  NUMCELLSIN20KM2,
     &	   			                  LU1,LU2,TRAPERR)
202         CONTINUE
500       CONTINUE			! NSEQ=1,TOTSEQ
350     CONTINUE				! IWC=1,WTYPES(ISERIES)
300   CONTINUE				! ISERIES=1,STYPES
C
C Output results for this cell
C
      IF(RUNFUTMET.EQ.0)THEN
c        CALL PUT_GIS_RES(SQUID,EAST,NORTH,NUMCELLSIN20KM2,KM20GRIDID,
c     &                       DOMSOIL,PERSOIL,SOILID,STYPES,WETCLASS,
c     &                       TOTSEQ,LUSEQ,SEQTYPE,
c     &                       LU1TOLU2,CCHANGE10,CH4C10,CO2C10,N2ON10,
c     &                       TRAPERR,DOMSOILC,DOMSOILBD,DOMSOILCLAY,
c     &                       SOMDEPTH,gley,PIANN,lu1)	 	 
	ELSEIF(RUNFUTMET.EQ.1)THEN
c	  CALL PUT_GIS_FUTRES(SQUID,EAST,NORTH,NUMCELLSIN20KM2,KM20GRIDID,
c     &                       DOMSOIL,PERSOIL,SOILID,STYPES,WETCLASS,
c     &                       TOTSEQ,LUSEQ,SEQTYPE,FRACLU50,
c     &                       CCHANGE10,CH4C10,CO2C10,N2ON10,
c     &                       TRAPERR,DOMSOILC,DOMSOILBD,SOMDEPTH,gley)
	ENDIF

C
C Spatial simulation of cells
C
      GOTO 101
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
c      IF(ISWAIT)THEN
c	  WRITE(*,*)'                       ....Press any key to continue'
c	  READ(*,*)
c	ENDIF
	END

C*************************************************************
C CROP ROUTINES
C*************************************************************
C
C EXTERNAL SUBROUTINES
C 1. GET_LUSEQ
C 2. GET_NLIM
C 3. INIT_GIS_CROP
C 4. RUN1_GIS_CROP
C 5. RUN2_GIS_CROP
C 6. SET_NLIM
C
C-------------------------------------------------------------
C
      SUBROUTINE GET_LUSEQ(LUSEQ,NGROW,NSEQ,SEQTYPE)
C
C Get land use sequence
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLU1				! Max.no.of land use types
	DATA MAXLU1 /6/						
	INTEGER MAXGROW				! Max.no.of growing seasons allowed
	PARAMETER (MAXGROW=300)
	INTEGER ILU1				! Land use counter
	INTEGER ILU2				! Land use counter
	INTEGER IGROW				! Counter for growing season 
	INTEGER TOTAL				! Number used in calculating sequence
	INTEGER NDIV				! Number used in calculating sequence
      INTEGER LU40S				! Land use in the 40s
      INTEGER LU50S				! Land use in the 50s
	INTEGER LU60S				! Land use in the 60s
	INTEGER LU70S				! Land use in the 70s
	INTEGER LU80S				! Land use in the 80s
	INTEGER LU90S				! Lans use in the 90s
C
C Variables passed to/from other subroutines
C
      INTEGER LUSEQ(MAXGROW)			! Land use 
									! sequence over each growing season
	INTEGER NSEQ					! Counter for number of land use sequences 
	INTEGER NGROW					! Number of growing seasons
	INTEGER SEQTYPE				! Type of sequence to be used
	INTEGER FULLSEQ				! Full sequences used (all options included)
	INTEGER DECSEQ				! One decade simulation only
	DATA FULLSEQ,DECSEQ /1,2/
C
C FULL SEQUENCE CALCULATION
C Work out land use change in the 50s, 60s and 70s
C
      IF(SEQTYPE.EQ.FULLSEQ)THEN 
        NGROW=30
        TOTAL=NSEQ-1
        NDIV=(MAXLU1)*(MAXLU1)*(MAXLU1)
        LU40S=TOTAL/NDIV
	  TOTAL=TOTAL-(LU40S*NDIV)
	  NDIV=(MAXLU1)*(MAXLU1)
	  LU50S=TOTAL/NDIV
	  TOTAL=TOTAL-(LU50S*NDIV)
	  NDIV=(MAXLU1)
	  LU60S=TOTAL/NDIV
	  TOTAL=TOTAL-(LU60S*NDIV)
        NDIV=1
	  LU70S=TOTAL/NDIV
	  TOTAL=TOTAL-(LU70S*NDIV)
C
C CALCULATION FOR 1 DECADE ONLY
C
      ELSEIF(SEQTYPE.EQ.DECSEQ)THEN
        NGROW=11
        TOTAL=NSEQ-1
	  NDIV=(MAXLU1)
	  LU40S=TOTAL/NDIV
	  TOTAL=TOTAL-(LU40S*NDIV)
        NDIV=1
	  LU50S=TOTAL/NDIV
	  TOTAL=TOTAL-(LU50S*NDIV)
	ENDIF
C
C For each growing year
C      
      DO 100 IGROW=1,NGROW
	  IF(IGROW.LE.5)THEN
          LUSEQ(IGROW)=LU40S+1
	  ELSEIF(IGROW.GT.5.AND.IGROW.LE.15)THEN
	    LUSEQ(IGROW)=LU50S+1
	  ELSEIF(IGROW.GT.15.AND.IGROW.LE.25)THEN
	    LUSEQ(IGROW)=LU60S+1
	  ELSEIF(IGROW.GT.25)THEN
	    LUSEQ(IGROW)=LU70S+1
	  ENDIF
100   CONTINUE
      END


C
C-------------------------------------------------------------
C
      SUBROUTINE GET_MODPI_LAYERED(AMMN,AMMN15,ATM,AVEPET,AVERAIN,
     &                AVETEMP,BCARB0,BNIT0,BNLAB0,CH4MODEL,CNMODEL,CTOT,
     &			    CLOSSX,CLOSSX15,CRIT,DMODEL,DOCMODEL,DPMCARB0,
     &			    DPMNIT0,DPMNLAB0,DRAINW,FERT,FIXN,FLOWPROP,
     &                FYMFERT,FYMFERT15,HCARB0,HNIT0,HNLAB0,IANTHES,
     &                IEND,IFERT,IHARV,ILAB,IMFUNC,ICMODEL,INMODEL,
     &                IOLAB,IOM,IORGM,ISERIES,ISOWN,ITFUNC,IVOL,
     &                JORGM,JSTOP,LHARV,LU1,MEND,
     &                NFERT,NORGM,N_STEPS,NSOIL,NSOW,ORGMA,
     &                PNREQ,REFIL,RPMCARB0,RPMNIT0,RPMNLAB0,SECONDS,
     &                SOIL15,SOILN,SOILW,TAM,TFERT,THISFERT,
     &                THISFERT15,TIN,TOC,TOTPIC,VOLAT,VOLAT15,WLEACH,
     &                WLEACH15,WSAT,WMAX,SOILPH,NXYEARS,EQMODEL,CLAY,
     &                PI_CEQ_MON,DOMSOILISIMP,DOMSOILIMPDEPTH,
     &                SOMDEPTH,NSOMLAY,LTA_AWC,LTA_TEMP,PH_MODEL)		  
C
C Modify plant inputs to achieve steady state in full ECOSSE run
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/
	INTEGER MAXFERT				! Max.no.of fertiliser applications allowed 
	PARAMETER (MAXFERT=5)
	INTEGER MAXGROW				! Max.no.of growing seasons allowed
	PARAMETER (MAXGROW=300)
	INTEGER MAXITER1			! Maximum allowed number of iterations
      DATA MAXITER1 /1000/
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXORGM				! Max.no.of organic manure applications allowed
	PARAMETER (MAXORGM=52)
	INTEGER MAXSEQ				! Max.no of sequences included in GIS
	PARAMETER (MAXSEQ=MAXLU*MAXLU)						
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXSOIL				! Max.no.of soil types
	PARAMETER (MAXSOIL=50)
	INTEGER MAXSOMLAY			! Max.no.of soil layers entered			
	PARAMETER (MAXSOMLAY=10)

      REAL*8 DC1(MAXLAYER)		! Change in C over simulation with MSS = 1 (kg C / ha)
      REAL*8 DC2(MAXLAYER)		! Change in C over simulation with calculated MSS (kg C / ha)
	REAL*8 DCTOT				! Change in C over simulation over whole measured profile (kg C / ha)
	REAL BTOT, HTOT, RTOT, DTOT ! 230909
	REAL ENDC(MAXLAYER)			! C at end of full simulation (kg C / ha)
	REAL ENDC0(MAXLAYER)		! C at end of simulation with no PI (kg C / ha)
	REAL ENDC90(MAXLAYER)		! C at end of simulation with 90% PI (kg C / ha)
	REAL ENDC100(MAXLAYER)		! C at end of simulation with 100% PI (kg C / ha)
      REAL FNULL					
      REAL FNULLARRAY(MAXLAYER)					
	INTEGER IYEAR				! Current growing season number - set to 1
	INTEGER IK					! No.timesteps from prev.harv. to current
      INTEGER IL					! Local layer counter
	INTEGER INULL				
	INTEGER ICFIXED				! Initialisation of C is fixed
	INTEGER ICROTHCEQ			! Initialisation of C by RothC 
								!  equilibrium run
	DATA ICFIXED,ICROTHCEQ /1,2/
	INTEGER INSTABLECN			! Initisalisation of N by assuming 
								!  steady state (after Bradbury)
	INTEGER INPASSCN			! Initialisation of N by passing the
								!  C:N ratio of DPM and RPM
	DATA INSTABLECN,INPASSCN /1,2/
	INTEGER ISATEQ				! Code to indicate whether equilibrium 
								! has been found 0=No 1=Yes
	INTEGER ISTEST				! Output balance sheet files 
								! 0 = Do not output test files
	                            ! 1 = Output balance sheet files (BALANCE_N.OUT, BALANCE_C.OUT & SOILW.OUT)
	INTEGER MEASLAY				! Layer that soil is measured to
	REAL*8 MSS(MAXLAYER)			! Modifier needed to acieve steady state
	REAL*8 MSS2					! Plant input modifier in second cycle
	INTEGER NITER				! No.of iterations in plant input optimisation
	INTEGER NITER2				! No.of iterations in soil pool optimisation
	INTEGER OPTCYCLE			! Number of years in the optimisation cycle
	REAL PC_TS(12)				! Proportion of total C litter input added each timestep (kgC/ha)
      REAL PI_CEQ_IN(12,MAXLAYER)	! Total annual plant input before modification
	REAL STARTC(MAXLAYER)		! C at start of full simulation (kg C / ha)
      REAL STARTPIC(MAXLAYER)		! Total annual plant input before modification
	REAL SUMC_TS				! Total C litter input added each year (kgC/ha/year)
      INTEGER TESTLU				! LU to be output
      INTEGER THISLU				! This LU

C
C Variables passed to/from other subroutines
C
	REAL AIRTEMP				! IN: Air temperature (deg.C/timestep)	
	REAL ALPHA(MAXLAYER)			! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL ANIT					! Ammonium immobilised (kgN/ha)
	REAL ANIT15					! Ammonium N15 nitrified (kgN/ha)
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
      REAL ATM					! IN:Atmospheric N input (kgN/ha/timestep)
	REAL ATM15					! IN:Atmospheric N15 input (kgN15/ha/time)
      REAL AVEPET(12)				! IN:Long term average PET (mm/month)
      REAL AVERAIN(12)			! IN:Long term average rainfall(mm/month)
      REAL AVETEMP(12)			! IN:Long term average monthly average 
	REAL BCARB(MAXLAYER)
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
	REAL BULKDENS(MAXLAYER)				! Soil bulk density (g cm-3)
	REAL BYRAIN					! IN/OUT:Rain lost by bypass flow mm/timestep	
	REAL C_TS					! IN:Litter C input in this timestep (kgC/ha)
	REAL CACT					! IN:Crop N uptake (kgN/ha/timestep)
	REAL CACT15					! IN:Crop N15 uptake (kgN15/ha/timestep)
	REAL CACTOT					! IN:N taken up by crop (kgN/ha)
	REAL CATOT15				! IN:N15 taken up by crop (kgN15/ha)
	REAL CH4(MAXLAYER)			! IN:CH4 emitted (kgC/ha/layer) (Aitkenhead CH4 model)
	INTEGER CH4MODEL			! OUT:Methane model (off, Richards or Aitkenhead model) 
	REAL CH4TOAIR				! OUT:CH4 released to atmosphere (kgC/ha) (Aitkenhead CH4 model)
	REAL CH4_PROD(MAXLAYER)  	! OUT: CH4 produced in each soil layer [kgC/ha/layer] (Richards CH4 model)
	REAL CH4_SOIL_OX			! OUT: CH4 oxidised from methane produced in the soil [kgC/ha] (Richards CH4 model)
      REAL CH4_ATMOS_OX			! OUT: CH4 oxidised from atmospheric CH4 [kgC/ha] (Richards CH4 model)
      REAL CH4_FLUX				! OUT: Net CH4 flux from soil column [kgC/ha] (Richards CH4 model)
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: % clay content in this layer
      REAL CLAY(MAXLAYER)			! IN:Clay content of the layer (%)
	REAL CLOSSX					! IN:N lost by senescence (kgN/ha) 
	REAL CLOSSX15				! IN:N15 lost by senescence (kgN15/ha) 	
	INTEGER CNMODEL				! IN:Type of Model for C:N ratio of DPM and 
	                            !     RPM (MAGEC or Foereid)
	REAL CO2(MAXLAYER)			! IN:CO2 emitted (kgC/ha/layer)
	REAL CO2FROMDOC(MAXLAYER)	! IN:CO2 output from DOC in each layer(kgC/ha)
	REAL CONC					! IN:Conc.N leached (kgN/ha/mm/layer)
	REAL CONC15					! IN:Conc.N15 leached (kgN/ha/mm/layer)
      REAL CONVER_F					! IN/OUT:Conversion to correct timestep
	REAL CRIT(MAXSOIL,MAXLAYER)	! IN:Critical N content (kgN/ha/layer)
	REAL CTOT					! IN:Total litter C input (kgC/ha)
	REAL CUPTN					! IN:Crop N uptake (kgN/ha/timestep)
 	REAL DELTA(MAXLAYER)			! IN/OUT: Prop.HUM produced on HUM decompn 
	REAL DENITRIFN(MAXLAYER)	! IN/OUT: N lost by denitrification (kgN/ha/layer
	REAL DENIT					! IN:N lost by denitrification (kgN/ha)
	REAL DFACT(MAXLAYER)		! IN/OUT: Denitrification factor
	INTEGER DMODEL				! IN:Type denitrification model chosen 
	                            ! (Bradbury or NEMIS)
	INTEGER DOCMODEL			! IN:DOC model (on or off)
	REAL DOMSOILIMPDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN: Depth of impermeable layer if present (cm)
	INTEGER DOMSOILISIMP(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN: Flag for impermeable layer (0=No; 1=Yes)
	REAL DN15					! IN:N15 lost by denitrification (kgN15/ha)
	REAL DNIT(MAXLAYER)			! IN:Net Mineralised N (kgN/ha/lay/timestep)
	REAL DNIT15(MAXLAYER)		! IN:Net Mineralised N15 (kgN15/ha/lay/time)
     	REAL DPMCARB(MAXLAYER)
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
      REAL DPMCTON(MAXLAYER)		! OUT:DPM C:N ratio
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL DPMRATE(MAXLAYER)			! IN/OUT: Rate constant for DPM decompn
	REAL DRAINW(MAXLAYER)		! IN:Water drained from this layer (mm/layer)
	REAL DRRAT						! IN:DPM:RPM ratio
	REAL ENDTOTC				! IN:Total organic C at end (kgC/ha)
	INTEGER EQMODEL				! IN:Type of equilibrium run 
	                            ! (NPP, TOC or both)
	REAL EVAPW					! IN: Potential evap. (mm/timestep)
	REAL FANIT					! Fertiliser N nitrified (kgN/ha)
	REAL FANIT15				! Fertiliser N15 nitrified (kgN/ha)
	REAL FERT(MAXFERT)			! IN:Amount of fertiliser applied (kgN/ha)
	REAL FERTNF					! INOUT:Amount of fertiliser applied (kgN/ha)
	REAL FIXN					! IN:N fixed from atmosphere (kgN/ha)
	REAL FLOWPROP				! IN:Proportion of flow needed to achieve
							    !	     water table at depth WTABLE
      REAL FPROPX(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
	REAL FYMAX(MAXORGM)			! IN/OUT:Proportion of ammonium-N in FYM
	REAL FYMCX(MAXORGM)			! IN/OUT:Proportion of C in FYM
	REAL FYMFERT				! IN:N input by fertiliser & FYM (kgN/ha)
	REAL FYMFERT15				! IN:N15 input by fertiliser & FYM (kgN15/ha)
	REAL FYMLOSS(MAXORGM)		! IN/OUT:Prop.FYM lost by volatn/timestep
	REAL FYMNX(MAXORGM)			! IN/OUT:Proportion of Organic N in FYM
	REAL FYMPOS(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
      REAL FYMSTART(MAXORGM)		! IN/OUT:Amount of rainfall in 1 week, 
								!		 below which volatilisation will occur
	REAL FYMWX(MAXORGM)			! IN/OUT:Amount of water added in FYM
	REAL GAMMA(MAXLAYER)			! IN/OUT: Prop.BIO produced on HUM decompn 
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
	REAL HCARB(MAXLAYER)
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
	INTEGER IANTHES				! IN:No.timesteps from 01/01 to anthesis 
  	REAL ICFACTOR(MAXLAYER)			! IN:Adjustment in rate needed to	
									!    achieve measured NPP and TOC	
	INTEGER ICMODEL				! OUT:Type of C initialisation model 
  	REAL ICNULL(MAXLAYER)			! IN:Adjustment in rate needed to	
									!    achieve measured NPP and TOC - not used	
      INTEGER ICOVER				! IN:Crop cover 1=Covered 0=Bare
	INTEGER IEND				! IN:Timesteps from sowing date to harvest
	INTEGER IFERT(MAXFERT)		! IN:No.timesteps to split fert.application
	INTEGER IHARV				! IN:Timesteps from 01/01/01 to harvest date 
	INTEGER ILAB(MAXFERT)		! IN:Labelling on fertiliser? 0=No 1=Yes
	INTEGER ILABNF				! INOUT:Labelling on fertiliser? 0=No 1=Yes
	INTEGER IMFUNC				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
	INTEGER INMODEL				! OUT:Type of N initialisation model
	INTEGER IOLAB(MAXORGM)		! IN:Labelling on manure? 0=No 1=Yes
	INTEGER IOLABNF				! IN:Labelling on manure? 0=No 1=Yes
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
	REAL IOMARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium IOM (kgC/ha/layer)
	INTEGER IORGM(MAXORGM)		! IN:No.timesteps to split manure application
	INTEGER IS_TS				! IN:Crop counter
      INTEGER ISAVE				! OUT:Code to save or retrieve variables
								!		0=save; 1=retrieve
	INTEGER ISERIES				! IN:Counter for soil series
	INTEGER ISOWN				! IN:Timesteps from 01/01/01 to sowing date 
	INTEGER ITFUNC				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IVOL(MAXFERT)		! IN:N volatilised from fertiliser (kgN/ha)
	INTEGER IVOLNF				! INOUT:N volatilised from fertiliser (kgN/ha)
	INTEGER JORGMNF				! IN:Type of manure application
	INTEGER JORGM(MAXORGM)		! IN:Type of manure application
	INTEGER JSTOP				! IN:Code to stop simulation 0=No 1=Yes
	INTEGER LHARV				! IN:No.timesteps from 01/01 to prev. harvest
	REAL LEACHDOC(MAXLAYER)		! IN/OUT:DOC leached from the layer (kgC/ha/layer/timestep)
	REAL LEACHDON(MAXLAYER)		! IN/OUT:DON leached from the layer (kgN/ha/layer/timestep)
      INTEGER LU1					! IN:Land use code
      INTEGER LUARRAY(MAXSOIL)	! IN/OUT:Land use before equilibrium 
	REAL LTA_AWC(12,MAXLAYER)	! IN:Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)	! IN:Average air temp this month (deg.C)
	INTEGER MEND				! IN:Number of timesteps from prev.harv.-harv.
	REAL MOBDOC(MAXLAYER)		! IN:Mobile DOC in each layer (kgC/ha)
	REAL MOBDON(MAXLAYER)		! IN:Mobile DON in each layer (kgN/ha)
	INTEGER N_STEPS				! IN:No.timesteps in a year
	INTEGER NFERT				! IN:No.fertiliser applications to this crop
	REAL NITRIFN(MAXLAYER)		! IN/OUT:Nitrification from this layer (kg/ha/layer/timestep)
	INTEGER NORGM				! IN:No.manure applications to this crop
	INTEGER NSOIL				! IN:Soil code number
	INTEGER NSOMLAY(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN:Number of SOM layers
	INTEGER NSOW				! IN:Timesteps from prev.harv. to sowing date
      INTEGER NUMSOIL					! Number of soils defined
	INTEGER NXYEARS				! IN:No.growing seasons simulated
	REAL OCROPN					! IN:N uptake of crop (kgN/ha) 
	REAL ORGMA(MAXORGM)			! IN:Amount of manure applied (t/ha)
	REAL ORGMANF				! IN:Amount of manure applied (t/ha)
	REAL ORGN					! IN:Total org.N input (kgN/ha) 
	INTEGER PH_MODEL			! IN:How is pH calculated?
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: pH above which decomp.max.
      REAL PI_C(MAXLAYER)			! IN: Plant input C to soil (kgC/ha/step)
	REAL PI_CEQ_MON(12,MAXLAYER)! IN/OUT: Equilibrium plant C input each month
								!         in each layer (kgC/ha/month/layer)
      REAL PI_N(MAXLAYER)			! IN: Plant input N to soil (kgC/ha/step)
      REAL PI_N15(MAXLAYER)		! IN: Plant input N15 to soil (kgC/ha/step)
	REAL PLANTUP(MAXLAYER)		! IN(RUN2_GIS_CROP)/OUT(TEST1_RES):Plant uptake from each soil layer
	REAL PNREQ(MAXLU)			! IN:Amount of N required by plant each yr
	                            !        (kgN/ha/year)
	REAL RAIN					! IN: Rainfall (mm/timestep)
	REAL REFIL(MAXLAYER)		! IN:Water needed to fill to FC (mm/lay)
	REAL RNIN					! IN:Litter N input in this timestep (kgN/ha)
	REAL RNIN15					! IN:Litter N15 input in timestep (kgN15/ha)
	REAL RPMCARB(MAXLAYER)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
      REAL RPMCTON(MAXLAYER)		! OUT:RPM C:N ratio
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL RPMRATE(MAXLAYER)			! IN/OUT: Rate constant for RPM decompn
	REAL SATWATCONT(MAXLAYER)	    ! IN:Total water content at saturation (mm/layer)
	REAL SECONDS				! IN:Number of seconds in one timestep
	REAL SEEDIN					! IN: Data used by MAGEC
	REAL SLEA15(MAXLAYER)		! IN:Nitrate-N15 leached (kgN15/ha)
	REAL SLEACH(MAXLAYER)		! IN:Nitrate-N leached (kgN/ha)
	CHARACTER*40 SNAME(MAXSOIL)		! IN/OUT:Soil name
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOILPH(MAXLAYER)		! IN/OUT:pH of soil in this layer
	REAL SOILTEMP(MAXLAYER)		! IN: Soil temperature (deg.C/timestep)
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL SOMDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN:depth of soil organic matter layers
	REAL SRNIN					! IN:Total litter N input (kgN/ha)
	REAL STARTTOTC				! IN:Total organic C at start (kgC/ha)
	REAL SX						! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SXORGN					! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL T15					! IN:Net Mineralised N (kgN/ha/layer)
	REAL TACTOT					! IN:Total crop N uptake (kgN/ha)
	REAL TAM(MAXLAYER)			! IN:15N/N in the ammonium in the layer
	REAL TC						! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL TCINP					! IN:Total org.C input (kgC/ha) DONT PASS?
	REAL TDNIT					! IN:Net Mineralised N (kgN/ha/layer)
	REAL TFERT(MAXFERT,3)		! IN:Prop.NO3,NH4,urea in fertiliser
	REAL TFERTNF(3)				! INOUT:Prop.NO3,NH4,urea in fertiliser
	REAL TFYMC					! IN:Total FYM C input (kgC/ha) DONT PASS?
	REAL THISFERT				! IN:N input by fertiliser (kgN/ha)
	REAL THISFERT15				! IN:N15 input by fertiliser (kgN15/ha)
	REAL TIN(MAXLAYER)			! IN:15N/N in the nitrate in the layer
	REAL TOC(MAXLAYER)			! IN:Measured total organic C (kgC/ha/layer)
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium TOC (kgC/ha/layer)

	REAL TOTPIC(MAXLU)			! IN:Plant C input calculated using 
	REAL TRATEM
	REAL TRNIN					! IN:Total litter N input (kgN/ha)
	REAL VOLAT					! IN:N lost by volatilisation (kgN/ha)
	REAL VOLAT15				! IN:N15 lost by volatilisation (kgN15/ha)
	REAL VOLATN(MAXLAYER)		! IN/OUT:Volatilisation from this layer (kg/ha/layer/timestep)
	REAL WILTPOINT(MAXLAYER)	! IN:Water content at wilting point (mm/layer)
	REAL WLEACH					! IN:N leached (kgN/ha)
	REAL WLEACH15				! IN:N15 leached (kgN15/ha)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WR						! IN:Root requirement (kgN/ha)
	REAL WRATEDM
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
	REAL XORGN					! IN:Total org.N input (kgN/ha) 
	REAL DDAYS
C
C Save current state of the soil 
C
      ISAVE=1 ! Save soil 
      CALL SAVE_GIS_SOIL(ISAVE,LU1,
     &                   SOILN,SOIL15,AMMN,AMMN15,
     &                   SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
     &                   DPMCARB0,DPMNIT0,DPMNLAB0,
     &                   RPMCARB0,RPMNIT0,RPMNLAB0,
     &                   BCARB0,BNIT0,BNLAB0,
     &                   HCARB0,HNIT0,HNLAB0,
     &                   IOM,PI_CEQ_MON)
C
C Calculate C at start of simulation
C
      NITER=0
      NITER2=0
	DO 100 IL=1,MAXLAYER1
	  STARTPIC(IL)=0
        MSS(IL)=1
	  DO 200 IK=1,12
	    STARTPIC(IL)=STARTPIC(IL)+PI_CEQ_MON(IK,IL)
	    PI_CEQ_IN(IK,IL)=PI_CEQ_MON(IK,IL)
200     CONTINUE
        STARTC(IL)=BCARB0(IL)+HCARB0(IL)+DPMCARB0(IL)+RPMCARB0(IL)
100   CONTINUE
      print*,'Start PIC= ',startpic(1)
C
C ...Calculate crop C and N returns and N offtake 
C
      IYEAR=1
	MEND=IHARV
      THISLU=LU1
C******************************
C Redistribute plant inputs according to SUNDIAL calculations
C ...Initialise current crop
C 
      CALL INIT_GIS_CROP(IYEAR,THISLU,SECONDS,
     &                     TOTPIC(THISLU),PNREQ(THISLU),
     &                     N_STEPS,ISOWN,NSOW,IHARV,LHARV,IANTHES,
     &                     MEND,IEND,NFERT,IFERT,ILAB,FERT,TFERT,
     &                     JSTOP,CTOT)
C
C ...Calculate crop C and N returns and N offtake 
C
      DO 250 IK=1,MEND
	  CALL RUN1_GIS_CROP(IYEAR,IK,THISLU,IEND,MEND,NSOW,
     &                          N_STEPS,IS_TS,TOTPIC(THISLU),
     &                          PNREQ(THISLU),CACTOT,CATOT15,
     &                          ICOVER,SXORGN,ORGN,
     &                          TRNIN,RNIN,RNIN15,C_TS,TCINP,WR)      
C
C ...Save carbon added this timestep
C
        PC_TS(IK)=C_TS
C
C ...Sum carbon added over year
C
        IF(IK.EQ.1)SUMC_TS=0
	  SUMC_TS=SUMC_TS+C_TS
250   CONTINUE
C
C ...Work out proportion of total C added each timestep
C
      DO 350 IK=1,MEND
	  PC_TS(IK)=PC_TS(IK)/SUMC_TS
C
C ...Correct distribution over year to match SUNDIAL equations
C
        DO 450 IL=1,MAXLAYER1
	    PI_CEQ_IN(IK,IL)=STARTPIC(IL)*PC_TS(IK)
	    PI_CEQ_MON(IK,IL)=STARTPIC(IL)*PC_TS(IK)
450     CONTINUE
350   CONTINUE
C******************************
C Full ECOSSE simulation using
C 1. 1 x PIeq
C 2. 0.9 x PIeq
C
101   CONTINUE
      NITER=NITER+1
C
C Adjust plant input according to modifier
C
      TOTPIC(LU1)=0
      DO 300 IL=1,MAXLAYER1
	  DO 400 IK=1,12
	    TOTPIC(LU1)=TOTPIC(LU1)+PI_CEQ_MON(IK,IL)
400     CONTINUE
300   CONTINUE
C
C For each growing season... 
C
      OPTCYCLE=NXYEARS 
c      OPTCYCLE=1000
      IYEAR=1
      DO WHILE (IYEAR.LE.OPTCYCLE)
C
C Initialise soil and crop in first year
C	
        IF(IYEAR.EQ.1)THEN
C
C ...initialise the soil for initial land use,
C
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
	    ENDIF
C
C Initialise current crop
C 
        CALL INIT_GIS_CROP(IYEAR,THISLU,SECONDS,
     &                     TOTPIC(THISLU),PNREQ(THISLU),
     &                     N_STEPS,ISOWN,NSOW,IHARV,LHARV,IANTHES,
     &                     MEND,IEND,NFERT,IFERT,ILAB,FERT,TFERT,
     &                     JSTOP,CTOT)
C
C For each week till end of growing season...
C           
        DO 500 IK=1,MEND
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
	    RAIN=AVERAIN(IK)
	    AIRTEMP=AVETEMP(IK)
	    EVAPW=AVEPET(IK)
          DO 600 IL=1,MAXLAYER1
	      SOILTEMP(IL)=AIRTEMP
600       CONTINUE
C
C ...Calculate crop C and N returns and N offtake 
C
	    CALL RUN1_GIS_CROP(IYEAR,IK,THISLU,IEND,MEND,NSOW,
     &                          N_STEPS,IS_TS,TOTPIC(THISLU),
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
          CALL ADD_PLANT_DIST(C_TS,RNIN,RNIN15,PI_CEQ_MON,IK,
     &                        PI_C,PI_N,PI_N15)
          IF(INMODEL.EQ.INPASSCN)
     &      CALL GET_CTON_FOEREID(DPMCTON,RPMCTON,THISLU)
	    CALL SET_DPMRPMRATIO(THISLU,DRRAT)
	    IF(DOMSOILISIMP(ISERIES,THISLU).eq.1)
     &      CALL LIMIT_DRAIN(DOMSOILIMPDEPTH(ISERIES,THISLU),SOILW)
********Temp change to remove N limitation************
c          do 9 il=1,maxlayer ! 230909
c 	    soiln(1)=soiln(il)+1000 ! 230909
c9         continue !230909
******************************************************
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
          CALL RUN_SUNDIAL_WATER_SUBMONTHLY(RAIN,WMAX,
     &                                      SOILW,WSAT,FLOWPROP,
     &                                      DRAINW,REFIL,EVAPW,
     &                                      SOILTEMP(1),ICOVER)
          CALL PHYSICAL_SUNDIAL_SOILCN(DOCMODEL,SECONDS,
     &                                   NSOIL,SOILTEMP,
     &                                   WMAX,SOILW,DRAINW,REFIL,
     &                                   SLEACH,SLEA15,
     &                                   WLEACH,WLEACH15,CONC,CONC15,
     &                                   SOILN,SOIL15,AMMN,AMMN15,
     &                                   TIN,TAM,MOBDOC,MOBDON,
     &                                   LEACHDOC,LEACHDON)
C
C Save results
C
        IF(IYEAR.EQ.OPTCYCLE.AND.IK.EQ.1)THEN
C          IF(NITER.GT.1)THEN ! 230909 
          IF(NITER.EQ.2)THEN	! 230909
C
C ...Calculate change from previous iteration
C
            DO 700 IL=1,MEASLAY
 	        DC1(IL)=ENDC(IL)-STARTC(IL)
700         CONTINUE
          ENDIF
          DO 850 IL=1,MEASLAY
            ENDC(IL)=BCARB0(IL)+HCARB0(IL)+DPMCARB0(IL)+RPMCARB0(IL)
850       CONTINUE
          IF(NITER.GT.1)THEN
            DO 750 IL=1,MEASLAY
              DC2(IL)=ENDC(IL)-STARTC(IL)
750         CONTINUE
	    ENDIF
	  ENDIF
C
C Write out results
C
        BTOT=0 ! 230909
	  HTOT=0 ! 230909
	  RTOT=0 ! 230909
	  DTOT=0 ! 230909
        DO 501 IL=1,MEASLAY ! 230909
          BTOT=BTOT+BCARB0(IL) ! 230909
		HTOT=HTOT+HCARB0(IL) ! 230909
		DTOT=DTOT+DPMCARB0(IL) ! 230909
		RTOT=RTOT+RPMCARB0(IL) ! 230909 	    
501     CONTINUE ! 230909
10      FORMAT(4F8.0) ! 230909

500     CONTINUE
        IYEAR=IYEAR+1
C
C Go back and set up the next crop
C
        END DO
C End of ECOSSE simulation
C******************************
C
C In first iteration (uses provided steady state plant inputs without adjustment)
C
      print*,'mss = ',mss(1)
      DCTOT=0 ! 230909
	DO 1600 IL=1,MEASLAY ! 230909
	  DCTOT=DCTOT+DC2(IL) ! 230909
1600  CONTINUE ! 230909
c	IF(NITER.LE.2.OR.DCTOT.GT.0.1.OR.DCTOT.LT.-0.1)THEN ! 230909
	IF(NITER.LE.2)THEN ! 230909

        ISAVE=0 ! Retrieve soil 
        CALL SAVE_GIS_SOIL(ISAVE,LU1,
     &                   SOILN,SOIL15,AMMN,AMMN15,
     &                   SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
     &                   DPMCARB0,DPMNIT0,DPMNLAB0,
     &                   RPMCARB0,RPMNIT0,RPMNLAB0,
     &                   BCARB0,BNIT0,BNLAB0,
     &                   HCARB0,HNIT0,HNLAB0,
     &                   IOM,PI_CEQ_MON)      

	ENDIF ! 230909
C
C After run using MSS=1,... 
C
      IF(NITER.EQ.1)THEN
C
C ... and reset layer multiplier to 0.9
C
        MSS2=0.9
	  DO 900 IL=1,MEASLAY
          MSS(IL)=MSS2
900     CONTINUE
        DO 1000 IL=MEASLAY+1,MAXLAYER1
	    MSS(IL)=MSS2
1000    CONTINUE
C
C After run using MSS=0.9, calculate change in C using modifier MSS=0.9 (DC2)
C and calculate adjustment to plant inputs according to PI = mss x PIeq
C where MSS = 1-(DC1 x (1-0.9))/(DC1-DC2)
C       DC1 = change in carbon using equilibrium run plant inputs
C   and DC2 = change in carbon using 0.9 x equilibrium run plant inputs
C
      ELSEIF(NITER.GT.1)THEN
	  DO 1100 IL=1,MEASLAY
	    MSS(IL)=(1-(DC1(IL)*(1-MSS(IL)))/(DC1(IL)-DC2(IL))) ! 230909
1100    CONTINUE
	ENDIF
C
C Work out total plant input according to previous proportions
C
      TOTPIC(THISLU)=0
      DO 1200 IL=1,MEASLAY
	  DO 1300 IK=1,12
	    IF(STARTC(IL).GT.0.AND.MSS(IL).GT.0)THEN
            PI_CEQ_MON(IK,IL)=MSS(IL)*PI_CEQ_IN(IK,IL)
	    ELSE
	      PI_CEQ_MON(IK,IL)=0
          ENDIF
	    TOTPIC(THISLU)=TOTPIC(THISLU)+PI_CEQ_MON(IK,IL)
1300    CONTINUE
1200  CONTINUE
      DO 1400 IL=MEASLAY+1,MAXLAYER1
	  MSS(IL)=1
	  DO 1500 IK=1,12
          PI_CEQ_MON(IK,IL)=MSS(IL)*PI_CEQ_IN(IK,IL)
	    TOTPIC(THISLU)=TOTPIC(THISLU)+PI_CEQ_MON(IK,IL)
1500    CONTINUE
1400  CONTINUE
C
C Redo calculation with recalculated plant input
C
C      DCTOT=0
C	DO 1600 IL=1,MEASLAY
C	  DCTOT=DCTOT+DC2(IL)
C1600  CONTINUE
      IF(NITER.EQ.1)PRINT*,'CALCULATION OF STEADY STATE' ! 230909
      IF(NITER.GT.1)PRINT*,NITER,'. dC = ',DCTOT! 230909
	IF(NITER.LE.2)THEN ! 230909
	  IF(NITER.LE.2)THEN 
	    GOTO 101 
	  ELSE ! 230909
	    IF(DCTOT.GT.0.1.OR.DCTOT.LT.-0.1)THEN ! 230909
	      GOTO 101 ! 230909
	    ENDIF ! 230909
	  ENDIF ! 230909
	ENDIF
C
C	When best fit with these pools found, save result and redo fit with new pools
C
      ISAVE=1 ! Save soil ! 230909
      CALL SAVE_GIS_SOIL(ISAVE,LU1,
     &                   SOILN,SOIL15,AMMN,AMMN15,
     &                   SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
     &                   DPMCARB0,DPMNIT0,DPMNLAB0,
     &                   RPMCARB0,RPMNIT0,RPMNLAB0,
     &                   BCARB0,BNIT0,BNLAB0,
     &                   HCARB0,HNIT0,HNLAB0,
     &                   IOM,PI_CEQ_MON) ! 230909
	IF(NITER2.LT.2)THEN
        NITER2=NITER2+1
	  NITER=0
	  DO 1700 IL=1,MEASLAY
	    STARTPIC(IL)=0
          MSS(IL)=1
	    DO 1800 IK=1,12
	      STARTPIC(IL)=STARTPIC(IL)+PI_CEQ_MON(IK,IL)
	      PI_CEQ_IN(IK,IL)=PI_CEQ_MON(IK,IL)
1800      CONTINUE
         STARTC(IL)=BCARB0(IL)+HCARB0(IL)+DPMCARB0(IL)+RPMCARB0(IL)
1700    CONTINUE
      print*,'End PIC= ',startpic(1)
        GOTO 101
	ENDIF
C
C Leave GET_MODPI_LAYERED2
C
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE GET_MODDEC_LAYERED(AMMN,AMMN15,ATM,AVEPET,AVERAIN,
     &                AVETEMP,BCARB0,BNIT0,BNLAB0,CH4MODEL,CNMODEL,CTOT,
     &			    CLOSSX,CLOSSX15,CRIT,DMODEL,DOCMODEL,DPMCARB0,
     &			    DPMNIT0,DPMNLAB0,DRAINW,FERT,FIXN,FLOWPROP,
     &                FYMFERT,FYMFERT15,HCARB0,HNIT0,HNLAB0,IANTHES,
     &                IEND,IFERT,IHARV,ILAB,IMFUNC,ICMODEL,INMODEL,
     &                IOLAB,IOM,IORGM,ISERIES,ISOWN,ITFUNC,IVOL,
     &                JORGM,JSTOP,LHARV,LU1,MEND,
     &                NFERT,NORGM,N_STEPS,NSOIL,NSOW,ORGMA,
     &                PNREQ,REFIL,RPMCARB0,RPMNIT0,RPMNLAB0,SECONDS,
     &                SOIL15,SOILN,SOILW,TAM,TFERT,THISFERT,
     &                THISFERT15,TIN,TOC,TOTPIC,VOLAT,VOLAT15,WLEACH,
     &                WLEACH15,WSAT,WMAX,SOILPH,NXYEARS,EQMODEL,CLAY,
     &                PI_CEQ_MON,DOMSOILISIMP,DOMSOILIMPDEPTH,
     &                SOMDEPTH,NSOMLAY,LTA_AWC,LTA_TEMP,PH_MODEL)		  
C
C Modify plant inputs to achieve steady state in full ECOSSE run
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/
	INTEGER MAXFERT				! Max.no.of fertiliser applications allowed 
	PARAMETER (MAXFERT=5)
	INTEGER MAXGROW				! Max.no.of growing seasons allowed
	PARAMETER (MAXGROW=300)
	INTEGER MAXITER1			! Maximum allowed number of iterations
      DATA MAXITER1 /1000/
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXORGM				! Max.no.of organic manure applications allowed
	PARAMETER (MAXORGM=52)
	INTEGER MAXSEQ				! Max.no of sequences included in GIS
	PARAMETER (MAXSEQ=MAXLU*MAXLU)						
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXSOIL				! Max.no.of soil types
	PARAMETER (MAXSOIL=50)
	INTEGER MAXSOMLAY			! Max.no.of soil layers entered			
	PARAMETER (MAXSOMLAY=10)

      REAL*8 DC1(MAXLAYER)		! Change in C over simulation with MSS = 1 (kg C / ha)
      REAL*8 DC2(MAXLAYER)		! Change in C over simulation with MSS = 0.9 (kg C / ha)
	REAL ENDC(MAXLAYER)			! C at end of full simulation (kg C / ha)
	REAL ENDC0(MAXLAYER)		! C at end of simulation with no PI (kg C / ha)
	REAL ENDC90(MAXLAYER)		! C at end of simulation with 90% PI (kg C / ha)
	REAL ENDC100(MAXLAYER)		! C at end of simulation with 100% PI (kg C / ha)
      REAL FNULL					
      REAL FNULLARRAY(MAXLAYER)					
	INTEGER IYEAR				! Current growing season number - set to 1
	INTEGER IK					! No.timesteps from prev.harv. to current
      INTEGER IL					! Local layer counter
	INTEGER INULL				
	INTEGER ISATEQ				! Code to indicate whether equilibrium 
								! has been found 0=No 1=Yes
	INTEGER ISTEST				! Output balance sheet files 
								! 0 = Do not output test files
	                            ! 1 = Output balance sheet files (BALANCE_N.OUT, BALANCE_C.OUT & SOILW.OUT)
	INTEGER MEASLAY				! Layer that soil is measured to
	REAL*8 MSS(MAXLAYER)			! Modifier needed to acieve steady state
	REAL*8 MSS2					! Plant input modifier in second cycle
	INTEGER NITER				! No.of iterations
	INTEGER OPTCYCLE			! Number of years in the optimisation cycle
	REAL PC_TS(12)				! Proportion of total C litter input added each timestep (kgC/ha)
      REAL PI_CEQ_IN(12,MAXLAYER)	! Total annual plant input before modification
	REAL STARTC(MAXLAYER)		! C at start of full simulation (kg C / ha)
      REAL STARTPIC(MAXLAYER)		! Total annual plant input before modification
	REAL SUMC_TS				! Total C litter input added each year (kgC/ha/year)
      INTEGER TESTLU				! LU to be output
      INTEGER THISLU				! This LU

	REAL BIN(MAXLAYER)			! C in soil biomass at start of optimisation (kgC/ha/layer)
	REAL DIN(MAXLAYER)			! C in soil dpm at start of optimisation (kgC/ha/layer)
	REAL HIN(MAXLAYER)			! C in soil humus at start of optimisation (kgC/ha/layer)
	REAL RIN(MAXLAYER)			! C in soil rpm at start of optimisation (kgC/ha/layer)
  	REAL ICFACTIN(MAXLAYER)		! Starting value for adjustment in rate needed to	
								!    achieve measured NPP and TOC	
	INTEGER ICFIXED				! Initialisation of C is fixed
	INTEGER ICROTHCEQ			! Initialisation of C by RothC 
								!  equilibrium run
	DATA ICFIXED,ICROTHCEQ /1,2/
	INTEGER INSTABLECN			! Initisalisation of N by assuming 
								!  steady state (after Bradbury)
	INTEGER INPASSCN			! Initialisation of N by passing the
								!  C:N ratio of DPM and RPM
	DATA INSTABLECN,INPASSCN /1,2/
C
C Variables passed to/from other subroutines
C
	REAL AIRTEMP				! IN: Air temperature (deg.C/timestep)	
	REAL ATM15					! IN:Atmospheric N15 input (kgN15/ha/time)
      REAL AVEPET(12)				! IN:Long term average PET (mm/month)
      REAL AVERAIN(12)			! IN:Long term average rainfall(mm/month)
      REAL AVETEMP(12)			! IN:Long term average monthly average 
      REAL BALANCE(20)			! IN:C Results
	REAL BULKDENS(MAXLAYER)				! Soil bulk density (g cm-3)
	REAL C_TS					! IN:Litter C input in this timestep (kgC/ha)
	REAL CACT					! IN:Crop N uptake (kgN/ha/timestep)
	REAL CACT15					! IN:Crop N15 uptake (kgN15/ha/timestep)
	REAL CACTOT					! IN:N taken up by crop (kgN/ha)
	REAL CATOT15				! IN:N15 taken up by crop (kgN15/ha)
	REAL CH4(MAXLAYER)			! IN:CH4 emitted (kgC/ha/layer) (Aitkenhead CH4 model)
	REAL CH4TOAIR				! OUT:CH4 released to atmosphere (kgC/ha) (Aitkenhead CH4 model)
	REAL CH4_PROD(MAXLAYER)  	! OUT: CH4 produced in each soil layer [kgC/ha/layer] (Richards CH4 model)
	REAL CH4_SOIL_OX			! OUT: CH4 oxidised from methane produced in the soil [kgC/ha] (Richards CH4 model)
      REAL CH4_ATMOS_OX			! OUT: CH4 oxidised from atmospheric CH4 [kgC/ha] (Richards CH4 model)
      REAL CH4_FLUX				! Net CH4 flux from soil column [kgC/ha] (Richards CH4 model)
      REAL CLAY(MAXLAYER)			! IN:Clay content of the layer (%)
	REAL CO2(MAXLAYER)			! IN:CO2 emitted (kgC/ha/layer)
	REAL CO2FROMDOC(MAXLAYER)	! IN:CO2 output from DOC in each layer(kgC/ha)
	REAL CONC					! IN:Conc.N leached (kgN/ha/mm/layer)
	REAL CONC15					! IN:Conc.N15 leached (kgN/ha/mm/layer)
	REAL DENITRIFN(MAXLAYER)	! IN/OUT: N lost by denitrification (kgN/ha/layer
	REAL DENIT					! IN:N lost by denitrification (kgN/ha)
	REAL DOMSOILIMPDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN: Depth of impermeable layer if present (cm)
	INTEGER DOMSOILISIMP(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN: Flag for impermeable layer (0=No; 1=Yes)
	REAL DN15					! IN:N15 lost by denitrification (kgN15/ha)
	REAL DNIT(MAXLAYER)			! IN:Net Mineralised N (kgN/ha/lay/timestep)
	REAL DNIT15(MAXLAYER)		! IN:Net Mineralised N15 (kgN15/ha/lay/time)
	REAL ENDTOTC				! IN:Total organic C at end (kgC/ha)
	INTEGER EQMODEL				! IN:Type of equilibrium run 
	                            ! (NPP, TOC or both)
	REAL EVAPW					! IN: Potential evap. (mm/timestep)
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
      INTEGER ICOVER				! IN:Crop cover 1=Covered 0=Bare
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
	INTEGER IS_TS				! IN:Crop counter
      INTEGER ISAVE				! OUT:Code to save or retrieve variables
								!		0=save; 1=retrieve
	INTEGER ISERIES				! IN:Counter for soil series
      INTEGER LUARRAY(MAXSOIL)	! IN/OUT:Land use before equilibrium 
	REAL MOBDOC(MAXLAYER)		! IN:Mobile DOC in each layer (kgC/ha)
	REAL MOBDON(MAXLAYER)		! IN:Mobile DON in each layer (kgN/ha)
	INTEGER NSOMLAY(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN:Number of SOM layers
	REAL OCROPN					! IN:N uptake of crop (kgN/ha) 
	REAL ORGMA(MAXORGM)			! IN:Amount of manure applied (t/ha)
	REAL ORGMANF				! IN:Amount of manure applied (t/ha)
	REAL ORGN					! IN:Total org.N input (kgN/ha) 
	INTEGER PH_MODEL			! IN:How is pH calculated?
      REAL PI_C(MAXLAYER)			! IN: Plant input C to soil (kgC/ha/step)
	REAL PI_CEQ_MON(12,MAXLAYER)! IN/OUT: Equilibrium plant C input each month
								!         in each layer (kgC/ha/month/layer)
      REAL PI_N(MAXLAYER)			! IN: Plant input N to soil (kgC/ha/step)
      REAL PI_N15(MAXLAYER)		! IN: Plant input N15 to soil (kgC/ha/step)
	REAL PLANTUP(MAXLAYER)		! IN(RUN2_GIS_CROP)/OUT(TEST1_RES):Plant uptake from each soil layer
	REAL RAIN					! IN: Rainfall (mm/timestep)
	REAL RNIN					! IN:Litter N input in this timestep (kgN/ha)
	REAL RNIN15					! IN:Litter N15 input in timestep (kgN15/ha)
	REAL SEEDIN					! IN: Data used by MAGEC
	REAL SLEA15(MAXLAYER)		! IN:Nitrate-N15 leached (kgN15/ha)
	REAL SLEACH(MAXLAYER)		! IN:Nitrate-N leached (kgN/ha)
	REAL SOILTEMP(MAXLAYER)		! IN: Soil temperature (deg.C/timestep)
	REAL SRNIN					! IN:Total litter N input (kgN/ha)
	REAL STARTTOTC				! IN:Total organic C at start (kgC/ha)
	REAL SX						! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SXORGN					! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL T15					! IN:Net Mineralised N (kgN/ha/layer)
	REAL TACTOT					! IN:Total crop N uptake (kgN/ha)
	REAL TC						! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL TCINP					! IN:Total org.C input (kgC/ha) DONT PASS?
	REAL TDNIT					! IN:Net Mineralised N (kgN/ha/layer)
	REAL TFYMC					! IN:Total FYM C input (kgC/ha) DONT PASS?
	REAL TOC(MAXLAYER)			! IN:Measured total organic C (kgC/ha/layer)
	REAL TRNIN					! IN:Total litter N input (kgN/ha)
	REAL WLEACH					! IN:N leached (kgN/ha)
	REAL WLEACH15				! IN:N15 leached (kgN15/ha)
	REAL WR						! IN:Root requirement (kgN/ha)
	REAL XORGN					! IN:Total org.N input (kgN/ha) 

	REAL HY(MAXSOIL,MAXLAYER)	! IN/OUT:Stable N:C ratio of BIO & HUM pools
	REAL BPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL HPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL BPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from biomass decompositn
	REAL HPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from humus decompositn
C
C Soil factors
C
	REAL ALPHA(MAXLAYER)			! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL BETA(MAXLAYER)				! IN/OUT: Prop.HUM produced on BIO decompn
	REAL GAMMA(MAXLAYER)			! IN/OUT: Prop.BIO produced on HUM decompn 
	REAL DELTA(MAXLAYER)			! IN/OUT: Prop.HUM produced on HUM decompn 
	REAL HZ1(MAXLAYER)				! IN/OUT: N:C ratio for steady state
	REAL BIOP(MAXSOIL,MAXLAYER)		! IN/OUT: BIO/TOC 
	REAL DPMRATE(MAXLAYER)			! IN/OUT: Rate constant for DPM decompn
	REAL DRRAT						! IN:DPM:RPM ratio
	REAL RPMRATE(MAXLAYER)			! IN/OUT: Rate constant for RPM decompn
	REAL BRATE(MAXLAYER)			! IN/OUT: Rate constant for HUM decompn
	REAL HRATE(MAXLAYER)			! IN/OUT: Rate constant for BIO decompn
      REAL CONVER_F					! IN/OUT:Conversion to correct timestep
	REAL NITRIFN(MAXLAYER)			! IN/OUT:Nitrification from this layer (kg/ha/layer/timestep)
      INTEGER NUMSOIL					! Number of soils defined
	REAL BIORATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for BIO decompn/yr
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for HUM decompn/yr
	REAL IOMARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium IOM (kgC/ha/layer)
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium TOC (kgC/ha/layer)
  	REAL ICFACTOR(MAXLAYER)			! IN:Adjustment in rate needed to	
									!    achieve measured NPP and TOC	
  	REAL ICNULL(MAXLAYER)			! IN:Adjustment in rate needed to	
									!    achieve measured NPP and TOC - not used	
	CHARACTER*40 SNAME(MAXSOIL)		! IN/OUT:Soil name
      REAL FPROPX(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
	REAL FYMCX(MAXORGM)			! IN/OUT:Proportion of C in FYM
	REAL FYMNX(MAXORGM)			! IN/OUT:Proportion of Organic N in FYM
	REAL FYMAX(MAXORGM)			! IN/OUT:Proportion of ammonium-N in FYM
	REAL FYMWX(MAXORGM)			! IN/OUT:Amount of water added in FYM
	REAL FYMLOSS(MAXORGM)		! IN/OUT:Prop.FYM lost by volatn/timestep
	REAL FYMPOS(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
      REAL FYMSTART(MAXORGM)		! IN/OUT:Amount of rainfall in 1 week, 
								!		 below which volatilisation will occur
	REAL BYRAIN					! IN/OUT:Rain lost by bypass flow mm/timestep	
	REAL ANIT					! Ammonium immobilised (kgN/ha)
	REAL ANIT15					! Ammonium N15 nitrified (kgN/ha)
	REAL FANIT					! Fertiliser N nitrified (kgN/ha)
	REAL FANIT15				! Fertiliser N15 nitrified (kgN/ha)
	REAL DFACT(MAXLAYER)		! IN/OUT: Denitrification factor
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: pH above which decomp.max.
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: % clay content in this layer
	REAL VOLATN(MAXLAYER)		! IN/OUT:Volatilisation from this layer (kg/ha/layer/timestep)
C
C Variables passed from calling routines
C
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
      REAL ATM					! IN:Atmospheric N input (kgN/ha/timestep)
	REAL BCARB0(MAXLAYER)		! IN:C in soil biomass (kgC/ha/layer)
      REAL BNIT0(MAXLAYER)		! IN:N in soil biomass (kgN/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN:N15 in soil humus (kgN15/ha/layer)
	INTEGER CH4MODEL			! OUT:Methane model (off, Richards or Aitkenhead model) 
	REAL CLOSSX					! IN:N lost by senescence (kgN/ha) 
	REAL CLOSSX15				! IN:N15 lost by senescence (kgN15/ha) 	
	INTEGER CNMODEL				! IN:Type of Model for C:N ratio of DPM and 
	                            !     RPM (MAGEC or Foereid)
	REAL CRIT(MAXSOIL,MAXLAYER)	! IN:Critical N content (kgN/ha/layer)
	REAL CTOT					! IN:Total litter C input (kgC/ha)
	REAL CUPTN					! IN:Crop N uptake (kgN/ha/timestep)
	INTEGER DMODEL				! IN:Type denitrification model chosen 
	                            ! (Bradbury or NEMIS)
	INTEGER DOCMODEL			! IN:DOC model (on or off)
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
      REAL DPMCTON(MAXLAYER)		! OUT:DPM C:N ratio
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL DRAINW(MAXLAYER)		! IN:Water drained from this layer (mm/layer)
	REAL FERT(MAXFERT)			! IN:Amount of fertiliser applied (kgN/ha)
	REAL FERTNF					! INOUT:Amount of fertiliser applied (kgN/ha)
	REAL FIXN					! IN:N fixed from atmosphere (kgN/ha)
	REAL FLOWPROP				! IN:Proportion of flow needed to achieve
							    !	     water table at depth WTABLE
	REAL FYMFERT				! IN:N input by fertiliser & FYM (kgN/ha)
	REAL FYMFERT15				! IN:N15 input by fertiliser & FYM (kgN15/ha)
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN:N in soil humus (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN:N15 in soil biomass (kgN15/ha/layer)
	INTEGER IANTHES				! IN:No.timesteps from 01/01 to anthesis 
	INTEGER ICMODEL				! OUT:Type of C initialisation model 
	INTEGER IEND				! IN:Timesteps from sowing date to harvest
	INTEGER IFERT(MAXFERT)		! IN:No.timesteps to split fert.application
	INTEGER IHARV				! IN:Timesteps from 01/01/01 to harvest date 
	INTEGER ILAB(MAXFERT)		! IN:Labelling on fertiliser? 0=No 1=Yes
	INTEGER ILABNF				! INOUT:Labelling on fertiliser? 0=No 1=Yes
	INTEGER IMFUNC				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
	INTEGER INMODEL				! OUT:Type of N initialisation model
	INTEGER IOLAB(MAXORGM)		! IN:Labelling on manure? 0=No 1=Yes
	INTEGER IOLABNF				! IN:Labelling on manure? 0=No 1=Yes
	INTEGER IORGM(MAXORGM)		! IN:No.timesteps to split manure application
	INTEGER ISOWN				! IN:Timesteps from 01/01/01 to sowing date 
	INTEGER ITFUNC				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IVOL(MAXFERT)		! IN:N volatilised from fertiliser (kgN/ha)
	INTEGER IVOLNF				! INOUT:N volatilised from fertiliser (kgN/ha)
	INTEGER JORGMNF				! IN:Type of manure application
	INTEGER JORGM(MAXORGM)		! IN:Type of manure application
	INTEGER JSTOP				! IN:Code to stop simulation 0=No 1=Yes
	REAL LEACHDOC(MAXLAYER)		! IN/OUT:DOC leached from the layer (kgC/ha/layer/timestep)
	REAL LEACHDON(MAXLAYER)		! IN/OUT:DON leached from the layer (kgN/ha/layer/timestep)
	INTEGER LHARV				! IN:No.timesteps from 01/01 to prev. harvest
      INTEGER LU1					! IN:Land use code
	REAL LTA_AWC(12,MAXLAYER)	! IN:Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)	! IN:Average air temp this month (deg.C)
	INTEGER MEND				! IN:Number of timesteps from prev.harv.-harv.
	INTEGER N_STEPS				! IN:No.timesteps in a year
	INTEGER NFERT				! IN:No.fertiliser applications to this crop
	INTEGER NORGM				! IN:No.manure applications to this crop
	INTEGER NSOIL				! IN:Soil code number
	INTEGER NSOW				! IN:Timesteps from prev.harv. to sowing date
	INTEGER NXYEARS				! IN:No.growing seasons simulated
	REAL PNREQ(MAXLU)			! IN:Amount of N required by plant each yr
	                            !        (kgN/ha/year)
	REAL REFIL(MAXLAYER)		! IN:Water needed to fill to FC (mm/lay)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
      REAL RPMCTON(MAXLAYER)		! OUT:RPM C:N ratio
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL SATWATCONT(MAXLAYER)	    ! IN:Total water content at saturation (mm/layer)
	REAL SECONDS				! IN:Number of seconds in one timestep
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOILPH(MAXLAYER)		! IN/OUT:pH of soil in this layer
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL SOMDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN:depth of soil organic matter layers
	REAL TAM(MAXLAYER)			! IN:15N/N in the ammonium in the layer
	REAL TFERT(MAXFERT,3)		! IN:Prop.NO3,NH4,urea in fertiliser
	REAL TFERTNF(3)				! INOUT:Prop.NO3,NH4,urea in fertiliser
	REAL THISFERT				! IN:N input by fertiliser (kgN/ha)
	REAL THISFERT15				! IN:N15 input by fertiliser (kgN15/ha)
	REAL TIN(MAXLAYER)			! IN:15N/N in the nitrate in the layer
	REAL TOTPIC(MAXLU)			! IN:Plant C input calculated using 
	REAL VOLAT					! IN:N lost by volatilisation (kgN/ha)
	REAL VOLAT15				! IN:N15 lost by volatilisation (kgN15/ha)
	REAL WILTPOINT(MAXLAYER)	! IN:Water content at wilting point (mm/layer)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
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
C
C Save current state of the soil 
C
      ISAVE=1 ! Save soil 
      CALL SAVE_GIS_SOIL(ISAVE,LU1,
     &                   SOILN,SOIL15,AMMN,AMMN15,
     &                   SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
     &                   DPMCARB0,DPMNIT0,DPMNLAB0,
     &                   RPMCARB0,RPMNIT0,RPMNLAB0,
     &                   BCARB0,BNIT0,BNLAB0,
     &                   HCARB0,HNIT0,HNLAB0,
     &                   IOM,PI_CEQ_MON)
C
C Calculate C at start of simulation
C
      NITER=0
	MEASLAY=SOMDEPTH(ISERIES,LU1,NSOMLAY(ISERIES,LU1))
	MEASLAY=MEASLAY*MAXLAYER1/MAXDEPTH
	DO 100 IL=1,MEASLAY
	  STARTPIC(IL)=0
        MSS(IL)=1
	  DO 200 IK=1,12
	    STARTPIC(IL)=STARTPIC(IL)+PI_CEQ_MON(IK,IL)
	    PI_CEQ_IN(IK,IL)=PI_CEQ_MON(IK,IL)
200     CONTINUE
        STARTC(IL)=BCARB0(IL)+HCARB0(IL)+DPMCARB0(IL)+RPMCARB0(IL)
c >> temp change >>
	  ICFACTIN(IL)=ICFACTOR(IL)
c       BIN(IL)=BCARB0(IL)
c	  HIN(IL)=HCARB0(IL)
c	  DIN(IL)=DPMCARB0(IL)
c	  RIN(IL)=RPMCARB0(IL)
c << temp change <<
100   CONTINUE
C
C ...Calculate crop C and N returns and N offtake 
C
      IYEAR=1
	MEND=IHARV
      THISLU=LU1
C******************************
C Redistribute plant inputs according to SUNDIAL calculations
C ...Initialise current crop
C 
      CALL INIT_GIS_CROP(IYEAR,THISLU,SECONDS,
     &                     TOTPIC(THISLU),PNREQ(THISLU),
     &                     N_STEPS,ISOWN,NSOW,IHARV,LHARV,IANTHES,
     &                     MEND,IEND,NFERT,IFERT,ILAB,FERT,TFERT,
     &                     JSTOP,CTOT)
C
C ...Calculate crop C and N returns and N offtake 
C
      DO 250 IK=1,MEND
	  CALL RUN1_GIS_CROP(IYEAR,IK,THISLU,IEND,MEND,NSOW,
     &                          N_STEPS,IS_TS,TOTPIC(THISLU),
     &                          PNREQ(THISLU),CACTOT,CATOT15,
     &                          ICOVER,SXORGN,ORGN,
     &                          TRNIN,RNIN,RNIN15,C_TS,TCINP,WR)      
C
C ...Save carbon added this timestep
C
        PC_TS(IK)=C_TS
C
C ...Sum carbon added over year
C
        IF(IK.EQ.1)SUMC_TS=0
	  SUMC_TS=SUMC_TS+C_TS
250   CONTINUE
C
C ...Work out proportion of total C added each timestep
C
      DO 350 IK=1,MEND
	  PC_TS(IK)=PC_TS(IK)/SUMC_TS
C
C ...Correct distribution over year to match SUNDIAL equations
C
        DO 450 IL=1,MAXLAYER1
	    PI_CEQ_IN(IK,IL)=STARTPIC(IL)*PC_TS(IK)
4	    PI_CEQ_MON(IK,IL)=STARTPIC(IL)*PC_TS(IK)
450     CONTINUE
350   CONTINUE
C******************************
C Full ECOSSE simulation using
C 1. 1 x PIeq
C 2. 0.9 x PIeq
C
101   CONTINUE
      NITER=NITER+1
C
C Adjust plant input according to modifier
C
      TOTPIC(LU1)=0
      DO 300 IL=1,MAXLAYER1
	  DO 400 IK=1,12
	    TOTPIC(LU1)=TOTPIC(LU1)+PI_CEQ_MON(IK,IL)
400     CONTINUE
300   CONTINUE
C
C For each growing season... 
C
      OPTCYCLE=NXYEARS
      IYEAR=1
      DO WHILE (IYEAR.LE.OPTCYCLE)
C
C Initialise soil and crop in first year
C	
        IF(IYEAR.EQ.1)THEN
C
C ...initialise the soil for initial land use,
C
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
	    ENDIF
C
C Initialise current crop
C 
        CALL INIT_GIS_CROP(IYEAR,THISLU,SECONDS,
     &                     TOTPIC(THISLU),PNREQ(THISLU),
     &                     N_STEPS,ISOWN,NSOW,IHARV,LHARV,IANTHES,
     &                     MEND,IEND,NFERT,IFERT,ILAB,FERT,TFERT,
     &                     JSTOP,CTOT)
C
C For each week till end of growing season...
C           
        DO 500 IK=1,MEND
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
	    RAIN=AVERAIN(IK)
	    AIRTEMP=AVETEMP(IK)
	    EVAPW=AVEPET(IK)
          DO 600 IL=1,MAXLAYER1
	      SOILTEMP(IL)=AIRTEMP
600       CONTINUE
C
C ...Calculate crop C and N returns and N offtake 
C
	    CALL RUN1_GIS_CROP(IYEAR,IK,THISLU,IEND,MEND,NSOW,
     &                          N_STEPS,IS_TS,TOTPIC(THISLU),
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
          CALL ADD_PLANT_DIST(C_TS,RNIN,RNIN15,PI_CEQ_MON,IK,
     &                        PI_C,PI_N,PI_N15)
          IF(INMODEL.EQ.INPASSCN)
     &      CALL GET_CTON_FOEREID(DPMCTON,RPMCTON,THISLU)
	    CALL SET_DPMRPMRATIO(THISLU,DRRAT)
c	    IF(DOMSOILISIMP(ISERIES,THISLU))
c     &      CALL LIMIT_DRAIN(DOMSOILIMPDEPTH(ISERIES,THISLU),SOILW)
********Temp change to remove N limitation************
          do 9 il=1,maxlayer ! 230909
 	    soiln(1)=soiln(il)+1000 ! 230909
9         continue !230909
******************************************************
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
          CALL RUN_SUNDIAL_WATER_SUBMONTHLY(RAIN,WMAX,
     &                                      SOILW,WSAT,FLOWPROP,
     &                                      DRAINW,REFIL,EVAPW,
     &                                      SOILTEMP(1),ICOVER)
          CALL PHYSICAL_SUNDIAL_SOILCN(DOCMODEL,SECONDS,
     &                                   NSOIL,SOILTEMP,
     &                                   WMAX,SOILW,DRAINW,REFIL,
     &                                   SLEACH,SLEA15,
     &                                   WLEACH,WLEACH15,CONC,CONC15,
     &                                   SOILN,SOIL15,AMMN,AMMN15,
     &                                   TIN,TAM,MOBDOC,MOBDON,
     &                                   LEACHDOC,LEACHDON)
C
C In first month of last year in optimisation cycle...
C
        IF(IYEAR.EQ.OPTCYCLE.AND.IK.EQ.1)THEN
          IF(NITER.GT.1)THEN
C
C ...Calculate change when modification factor = 1
C
            DO 700 IL=1,MEASLAY
 	        DC1(IL)=ENDC(IL)-STARTC(IL)
700         CONTINUE
          ENDIF
C
C ... Save end carbon content in each layer
C
          DO 850 IL=1,MEASLAY
            ENDC(IL)=BCARB0(IL)+HCARB0(IL)+DPMCARB0(IL)+RPMCARB0(IL)
850       CONTINUE
C
C ... and calculate change when modification factor = 0.9
C
          IF(NITER.GT.1)THEN
            DO 750 IL=1,MEASLAY
              DC2(IL)=ENDC(IL)-STARTC(IL)
750         CONTINUE
	    ENDIF
	  ENDIF
500     CONTINUE
C
C Go back and set up the next crop
C
        IYEAR=IYEAR+1
        END DO
C End of ECOSSE simulation
C******************************
C
C In first iteration (uses provided steady state plant inputs without adjustment)
C
C
C Retrieve soil characteristics from start of simulation
C
      ISAVE=0 ! Retrieve soil 
      CALL SAVE_GIS_SOIL(ISAVE,LU1,
     &                   SOILN,SOIL15,AMMN,AMMN15,
     &                   SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
     &                   DPMCARB0,DPMNIT0,DPMNLAB0,
     &                   RPMCARB0,RPMNIT0,RPMNLAB0,
     &                   BCARB0,BNIT0,BNLAB0,
     &                   HCARB0,HNIT0,HNLAB0,
     &                   IOM,PI_CEQ_MON)
      
C
C After run using MSS=1,... 
C
      MSS2=0.9
      IF(NITER.EQ.1)THEN
C
C ... and reset layer multiplier to 0.9
C
	  DO 900 IL=1,MEASLAY
          MSS(IL)=MSS2
900     CONTINUE
        DO 1000 IL=MEASLAY,MAXLAYER1
	    MSS(IL)=MSS2
1000    CONTINUE
C
C After run using MSS=0.9, calculate change in C using modifier MSS=0.9 (DC2)
C and calculate adjustment to plant inputs according to PI = mss x PIeq
C where MSS = 1-(DC1 x (1-0.9))/(DC1-DC2)
C       DC1 = change in carbon using equilibrium run plant inputs
C   and DC2 = change in carbon using 0.9 x equilibrium run plant inputs
C
      ELSEIF(NITER.GT.1)THEN
	  DO 1100 IL=1,MEASLAY
	    MSS(IL)=(1-(DC1(IL)*(1-MSS(IL)))/(DC1(IL)-DC2(IL)))
1100    CONTINUE
	ENDIF
C
C Work out total plant input according to previous proportions
C
      TOTPIC(THISLU)=0
      DO 1200 IL=1,MEASLAY
	  DO 1300 IK=1,12
	    IF(STARTC(IL).GT.0.AND.MSS(IL).GT.0)THEN
c*** >> temp change          
		  ICFACTOR(IL)=MSS(IL)*ICFACTIN(IL)
c	      BCARB0(IL)=MSS(IL)*BIN(IL)
c	      HCARB0(IL)=HIN(IL)-(BCARB0(IL)-BIN(IL))
c	      IF(HCARB0(IL).LT.0)THEN
c	        BCARB0(IL)=BCARB0(IL)+HCARB0(IL)
c	        HCARB0(IL)=0
c	      ENDIF
c	      DPMCARB0(IL)=MSS(IL)*DIN(IL)
c	      RPMCARB0(IL)=RIN(IL)-(DPMCARB0(IL)-DIN(IL))
c	      IF(RPMCARB0(IL).LT.0)THEN
c	        DPMCARB0(IL)=DPMCARB0(IL)+RPMCARB0(IL)
c	        RPMCARB0(IL)=0
c	      ENDIF
	    ELSE
	      ICFACTOR(IL)=ICFACTIN(IL)
          ENDIF
c*** << temp change
	    TOTPIC(THISLU)=TOTPIC(THISLU)+PI_CEQ_MON(IK,IL)
1300    CONTINUE
1200  CONTINUE
      DO 1400 IL=MEASLAY,MAXLAYER1
	  MSS(IL)=1
	  DO 1500 IK=1,12
          ICFACTOR(IL)=MSS(IL)*ICFACTIN(IL)
	    TOTPIC(THISLU)=TOTPIC(THISLU)+PI_CEQ_MON(IK,IL)
1500    CONTINUE
1400  CONTINUE
C
C Redo calculation with recalculated plant input
C
      IF(NITER.LT.2)THEN
	  GOTO 101
	ENDIF
C
C Save current state of the soil 
C
      ISAVE=1 ! Save soil 
      CALL SAVE_GIS_SOIL(ISAVE,LU1,
     &                   SOILN,SOIL15,AMMN,AMMN15,
     &                   SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
     &                   DPMCARB0,DPMNIT0,DPMNLAB0,
     &                   RPMCARB0,RPMNIT0,RPMNLAB0,
     &                   BCARB0,BNIT0,BNLAB0,
     &                   HCARB0,HNIT0,HNLAB0,
     &                   IOM,PI_CEQ_MON)
C
C Leave GET_MODDEC_LAYERED
C
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE GET_MODPI_UNLAYERED(AMMN,AMMN15,ATM,AVEPET,AVERAIN,
     &                AVETEMP,BCARB0,BNIT0,BNLAB0,CH4MODEL,CNMODEL,CTOT,
     &			    CLOSSX,CLOSSX15,CRIT,DMODEL,DOCMODEL,DPMCARB0,
     &			    DPMNIT0,DPMNLAB0,DRAINW,FERT,FIXN,FLOWPROP,
     &                FYMFERT,FYMFERT15,HCARB0,HNIT0,HNLAB0,IANTHES,
     &                IEND,IFERT,IHARV,ILAB,IMFUNC,ICMODEL,INMODEL,
     &                IOLAB,IOM,IORGM,ISERIES,ISOWN,ITFUNC,IVOL,
     &                JORGM,JSTOP,LHARV,LU1,MEND,
     &                NFERT,NORGM,N_STEPS,NSOIL,NSOW,ORGMA,
     &                PNREQ,REFIL,RPMCARB0,RPMNIT0,RPMNLAB0,SECONDS,
     &                SOIL15,SOILN,SOILW,TAM,TFERT,THISFERT,
     &                THISFERT15,TIN,TOC,TOTPIC,VOLAT,VOLAT15,WLEACH,
     &                WLEACH15,WSAT,WMAX,SOILPH,NXYEARS,EQMODEL,CLAY,
     &                PI_CEQ_MON,DOMSOILISIMP,DOMSOILIMPDEPTH,
     &                SOMDEPTH,NSOMLAY,LTA_AWC,LTA_TEMP,PH_MODEL)		  
C
C Modify plant inputs to achieve steady state in full ECOSSE run
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/
	INTEGER MAXFERT				! Max.no.of fertiliser applications allowed 
	PARAMETER (MAXFERT=5)
	INTEGER MAXGROW				! Max.no.of growing seasons allowed
	PARAMETER (MAXGROW=300)
	INTEGER MAXITER1			! Maximum allowed number of iterations
      DATA MAXITER1 /1000/
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXORGM				! Max.no.of organic manure applications allowed
	PARAMETER (MAXORGM=52)
	INTEGER MAXSEQ				! Max.no of sequences included in GIS
	PARAMETER (MAXSEQ=MAXLU*MAXLU)						
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXSOIL				! Max.no.of soil types
	PARAMETER (MAXSOIL=50)
	INTEGER MAXSOMLAY			! Max.no.of soil layers entered			
	PARAMETER (MAXSOMLAY=10)

      REAL DC1					! Change in C over simulation with MSS = 1 (kg C / ha)
      REAL DC2					! Change in C over simulation with MSS = 0.9 (kg C / ha)
	REAL ENDC					! C at end of full simulation (kg C / ha)
      REAL FNULL					
      REAL FNULLARRAY(MAXLAYER)					
	INTEGER IYEAR				! Current growing season number - set to 1
	INTEGER IK					! No.timesteps from prev.harv. to current
      INTEGER IL					! Local layer counter
	INTEGER INULL				
	INTEGER ISATEQ				! Code to indicate whether equilibrium 
								! has been found 0=No 1=Yes
	INTEGER ISTEST				! Output balance sheet files 
								! 0 = Do not output test files
	                            ! 1 = Output balance sheet files (BALANCE_N.OUT, BALANCE_C.OUT & SOILW.OUT)
	INTEGER MEASLAY				! Layer that soil is measured to
	REAL MSS					! Modifier needed to acieve steady state
	INTEGER NITER				! No.of iterations
	REAL STARTC					! C at start of full simulation (kg C / ha)
      REAL STARTPIC				! Total annual plant input before modification
      INTEGER TESTLU				! LU to be output
      INTEGER THISLU				! This LU

	INTEGER ICFIXED				! Initialisation of C is fixed
	INTEGER ICROTHCEQ			! Initialisation of C by RothC 
								!  equilibrium run
	DATA ICFIXED,ICROTHCEQ /1,2/
	INTEGER INSTABLECN			! Initisalisation of N by assuming 
								!  steady state (after Bradbury)
	INTEGER INPASSCN			! Initialisation of N by passing the
								!  C:N ratio of DPM and RPM
	DATA INSTABLECN,INPASSCN /1,2/
C
C Variables passed to/from other subroutines
C
	REAL AIRTEMP				! IN: Air temperature (deg.C/timestep)	
	REAL ATM15					! IN:Atmospheric N15 input (kgN15/ha/time)
      REAL AVEPET(12)				! IN:Long term average PET (mm/month)
      REAL AVERAIN(12)			! IN:Long term average rainfall(mm/month)
      REAL AVETEMP(12)			! IN:Long term average monthly average 
      REAL BALANCE(20)			! IN:C Results
	REAL BULKDENS(MAXLAYER)				! Soil bulk density (g cm-3)
	REAL C_TS					! IN:Litter C input in this timestep (kgC/ha)
	REAL CACT					! IN:Crop N uptake (kgN/ha/timestep)
	REAL CACT15					! IN:Crop N15 uptake (kgN15/ha/timestep)
	REAL CACTOT					! IN:N taken up by crop (kgN/ha)
	REAL CATOT15				! IN:N15 taken up by crop (kgN15/ha)
	REAL CH4(MAXLAYER)			! IN:CH4 emitted (kgC/ha/layer)  (Aitkenhead CH4 model)
	REAL CH4TOAIR				! OUT: CH4 released to atmosphere (kgC/ha) (Aitkenheaf CH4 model)
	REAL CH4_PROD(MAXLAYER)  	! OUT: CH4 produced in each soil layer [kgC/ha/layer] (Richards CH4 model)
	REAL CH4_SOIL_OX			! OUT: CH4 oxidised from methane produced in the soil [kgC/ha] (Richards CH4 model)
      REAL CH4_ATMOS_OX			! OUT: CH4 oxidised from atmospheric CH4 [kgC/ha] (Richards CH4 model)
      REAL CH4_FLUX				! OUT: Net CH4 flux from soil column [kgC/ha] (Richards CH4 model)
      REAL CLAY(MAXLAYER)			! IN:Clay content of the layer (%)
	REAL CO2(MAXLAYER)			! IN:CO2 emitted (kgC/ha/layer)
	REAL CO2FROMDOC(MAXLAYER)	! IN:CO2 output from DOC in each layer(kgC/ha)
	REAL CONC					! IN:Conc.N leached (kgN/ha/mm/layer)
	REAL CONC15					! IN:Conc.N15 leached (kgN/ha/mm/layer)
	REAL DENITRIFN(MAXLAYER)	! IN/OUT: N lost by denitrification (kgN/ha/layer
	REAL DENIT					! IN:N lost by denitrification (kgN/ha)
	REAL DOMSOILIMPDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN: Depth of impermeable layer if present (cm)
	INTEGER DOMSOILISIMP(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN: Flag for impermeable layer (0=No; 1=Yes)
	REAL DN15					! IN:N15 lost by denitrification (kgN15/ha)
	REAL DNIT(MAXLAYER)			! IN:Net Mineralised N (kgN/ha/lay/timestep)
	REAL DNIT15(MAXLAYER)		! IN:Net Mineralised N15 (kgN15/ha/lay/time)
	REAL ENDTOTC				! IN:Total organic C at end (kgC/ha)
	INTEGER EQMODEL				! IN:Type of equilibrium run 
	                            ! (NPP, TOC or both)
	REAL EVAPW					! IN: Potential evap. (mm/timestep)
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
      INTEGER ICOVER				! IN:Crop cover 1=Covered 0=Bare
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
	INTEGER IS_TS				! IN:Crop counter
      INTEGER ISAVE				! OUT:Code to save or retrieve variables
								!		0=save; 1=retrieve
	INTEGER ISERIES				! IN:Counter for soil series
      INTEGER LUARRAY(MAXSOIL)	! IN/OUT:Land use before equilibrium 
	REAL MOBDOC(MAXLAYER)		! IN:Mobile DOC in each layer (kgC/ha)
	REAL MOBDON(MAXLAYER)		! IN:Mobile DON in each layer (kgN/ha)
	INTEGER NSOMLAY(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN:Number of SOM layers
	REAL OCROPN					! IN:N uptake of crop (kgN/ha) 
	REAL ORGMA(MAXORGM)			! IN:Amount of manure applied (t/ha)
	REAL ORGMANF				! IN:Amount of manure applied (t/ha)
	REAL ORGN					! IN:Total org.N input (kgN/ha) 
	INTEGER PH_MODEL			! IN:How is pH calculated?
      REAL PI_C(MAXLAYER)			! IN: Plant input C to soil (kgC/ha/step)
	REAL PI_CEQ_MON(12,MAXLAYER)! IN/OUT: Equilibrium plant C input each month
								!         in each layer (kgC/ha/month/layer)
      REAL PI_N(MAXLAYER)			! IN: Plant input N to soil (kgC/ha/step)
      REAL PI_N15(MAXLAYER)		! IN: Plant input N15 to soil (kgC/ha/step)
	REAL RAIN					! IN: Rainfall (mm/timestep)
	REAL RNIN					! IN:Litter N input in this timestep (kgN/ha)
	REAL RNIN15					! IN:Litter N15 input in timestep (kgN15/ha)
	REAL SEEDIN					! IN: Data used by MAGEC
	REAL SLEA15(MAXLAYER)		! IN:Nitrate-N15 leached (kgN15/ha)
	REAL SLEACH(MAXLAYER)		! IN:Nitrate-N leached (kgN/ha)
	REAL SOILTEMP(MAXLAYER)		! IN: Soil temperature (deg.C/timestep)
	REAL SRNIN					! IN:Total litter N input (kgN/ha)
	REAL STARTTOTC				! IN:Total organic C at start (kgC/ha)
	REAL SX						! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SXORGN					! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL T15					! IN:Net Mineralised N (kgN/ha/layer)
	REAL TACTOT					! IN:Total crop N uptake (kgN/ha)
	REAL TC						! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL TCINP					! IN:Total org.C input (kgC/ha) DONT PASS?
	REAL TDNIT					! IN:Net Mineralised N (kgN/ha/layer)
	REAL TFYMC					! IN:Total FYM C input (kgC/ha) DONT PASS?
	REAL TOC(MAXLAYER)			! IN:Measured total organic C (kgC/ha/layer)
	REAL TRNIN					! IN:Total litter N input (kgN/ha)
	REAL WLEACH					! IN:N leached (kgN/ha)
	REAL WLEACH15				! IN:N15 leached (kgN15/ha)
	REAL WR						! IN:Root requirement (kgN/ha)
	REAL XORGN					! IN:Total org.N input (kgN/ha) 

	REAL HY(MAXSOIL,MAXLAYER)	! IN/OUT:Stable N:C ratio of BIO & HUM pools
	REAL BPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL HPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL BPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from biomass decompositn
	REAL HPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from humus decompositn
C
C Soil factors
C
	REAL ALPHA(MAXLAYER)			! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL BETA(MAXLAYER)				! IN/OUT: Prop.HUM produced on BIO decompn
	REAL GAMMA(MAXLAYER)			! IN/OUT: Prop.BIO produced on HUM decompn 
	REAL DELTA(MAXLAYER)			! IN/OUT: Prop.HUM produced on HUM decompn 
	REAL HZ1(MAXLAYER)				! IN/OUT: N:C ratio for steady state
	REAL BIOP(MAXSOIL,MAXLAYER)		! IN/OUT: BIO/TOC 
	REAL DPMRATE(MAXLAYER)			! IN/OUT: Rate constant for DPM decompn
	REAL DRRAT						! IN:DPM:RPM ratio
	REAL RPMRATE(MAXLAYER)			! IN/OUT: Rate constant for RPM decompn
	REAL BRATE(MAXLAYER)			! IN/OUT: Rate constant for HUM decompn
	REAL HRATE(MAXLAYER)			! IN/OUT: Rate constant for BIO decompn
      REAL CONVER_F					! IN/OUT:Conversion to correct timestep
      INTEGER NUMSOIL					! Number of soils defined
	REAL BIORATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for BIO decompn/yr
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for HUM decompn/yr
	REAL IOMARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium IOM (kgC/ha/layer)
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium TOC (kgC/ha/layer)
  	REAL ICFACTOR(MAXLAYER)			! IN:Adjustment in rate needed to	
									!    achieve measured NPP and TOC	
  	REAL ICNULL(MAXLAYER)			! IN:Adjustment in rate needed to	
									!    achieve measured NPP and TOC - not used	
	CHARACTER*40 SNAME(MAXSOIL)		! IN/OUT:Soil name
      REAL FPROPX(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
	REAL FYMCX(MAXORGM)			! IN/OUT:Proportion of C in FYM
	REAL FYMNX(MAXORGM)			! IN/OUT:Proportion of Organic N in FYM
	REAL FYMAX(MAXORGM)			! IN/OUT:Proportion of ammonium-N in FYM
	REAL FYMWX(MAXORGM)			! IN/OUT:Amount of water added in FYM
	REAL FYMLOSS(MAXORGM)		! IN/OUT:Prop.FYM lost by volatn/timestep
	REAL FYMPOS(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
      REAL FYMSTART(MAXORGM)		! IN/OUT:Amount of rainfall in 1 week, 
								!		 below which volatilisation will occur
	REAL BYRAIN					! IN/OUT:Rain lost by bypass flow mm/timestep	
	REAL ANIT					! Ammonium immobilised (kgN/ha)
	REAL ANIT15					! Ammonium N15 nitrified (kgN/ha)
	REAL FANIT					! Fertiliser N nitrified (kgN/ha)
	REAL FANIT15				! Fertiliser N15 nitrified (kgN/ha)
	REAL DFACT(MAXLAYER)		! IN/OUT: Denitrification factor
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: pH above which decomp.max.
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: % clay content in this layer
	REAL NITRIFN(MAXLAYER)		! IN/OUT:Nitrification from this layer (kg/ha/layer/timestep)
	REAL VOLATN(MAXLAYER)		! IN/OUT:Volatilisation from this layer (kg/ha/layer/timestep)
C
C Variables passed from calling routines
C
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
      REAL ATM					! IN:Atmospheric N input (kgN/ha/timestep)
	REAL BCARB0(MAXLAYER)		! IN:C in soil biomass (kgC/ha/layer)
      REAL BNIT0(MAXLAYER)		! IN:N in soil biomass (kgN/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN:N15 in soil humus (kgN15/ha/layer)
	INTEGER CH4MODEL			! OUT:Methane model (off, Richards or Aitkenhead model) 
	REAL CLOSSX					! IN:N lost by senescence (kgN/ha) 
	REAL CLOSSX15				! IN:N15 lost by senescence (kgN15/ha) 	
	INTEGER CNMODEL				! IN:Type of Model for C:N ratio of DPM and 
	                            !     RPM (MAGEC or Foereid)
	REAL CRIT(MAXSOIL,MAXLAYER)	! IN:Critical N content (kgN/ha/layer)
	REAL CTOT					! IN:Total litter C input (kgC/ha)
	REAL CUPTN					! IN:Crop N uptake (kgN/ha/timestep)
	INTEGER DMODEL				! IN:Type denitrification model chosen 
	                            ! (Bradbury or NEMIS)
	INTEGER DOCMODEL			! IN:DOC model (on or off)
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
      REAL DPMCTON(MAXLAYER)		! OUT:DPM C:N ratio
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL DRAINW(MAXLAYER)		! IN:Water drained from this layer (mm/layer)
	REAL FERT(MAXFERT)			! IN:Amount of fertiliser applied (kgN/ha)
	REAL FERTNF					! INOUT:Amount of fertiliser applied (kgN/ha)
	REAL FIXN					! IN:N fixed from atmosphere (kgN/ha)
	REAL FLOWPROP				! IN:Proportion of flow needed to achieve
							    !	     water table at depth WTABLE
	REAL FYMFERT				! IN:N input by fertiliser & FYM (kgN/ha)
	REAL FYMFERT15				! IN:N15 input by fertiliser & FYM (kgN15/ha)
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN:N in soil humus (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN:N15 in soil biomass (kgN15/ha/layer)
	INTEGER IANTHES				! IN:No.timesteps from 01/01 to anthesis 
	INTEGER ICMODEL				! OUT:Type of C initialisation model 
	INTEGER IEND				! IN:Timesteps from sowing date to harvest
	INTEGER IFERT(MAXFERT)		! IN:No.timesteps to split fert.application
	INTEGER IHARV				! IN:Timesteps from 01/01/01 to harvest date 
	INTEGER ILAB(MAXFERT)		! IN:Labelling on fertiliser? 0=No 1=Yes
	INTEGER ILABNF				! INOUT:Labelling on fertiliser? 0=No 1=Yes
	INTEGER IMFUNC				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
	INTEGER INMODEL				! OUT:Type of N initialisation model
	INTEGER IOLAB(MAXORGM)		! IN:Labelling on manure? 0=No 1=Yes
	INTEGER IOLABNF				! IN:Labelling on manure? 0=No 1=Yes
	INTEGER IORGM(MAXORGM)		! IN:No.timesteps to split manure application
	INTEGER ISOWN				! IN:Timesteps from 01/01/01 to sowing date 
	INTEGER ITFUNC				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IVOL(MAXFERT)		! IN:N volatilised from fertiliser (kgN/ha)
	INTEGER IVOLNF				! INOUT:N volatilised from fertiliser (kgN/ha)
	INTEGER JORGMNF				! IN:Type of manure application
	INTEGER JORGM(MAXORGM)		! IN:Type of manure application
	INTEGER JSTOP				! IN:Code to stop simulation 0=No 1=Yes
	REAL LEACHDOC(MAXLAYER)		! IN/OUT:DOC leached from the layer (kgC/ha/layer/timestep)
	REAL LEACHDON(MAXLAYER)		! IN/OUT:DON leached from the layer (kgN/ha/layer/timestep)
	INTEGER LHARV				! IN:No.timesteps from 01/01 to prev. harvest
      INTEGER LU1					! IN:Land use code
	REAL LTA_AWC(12,MAXLAYER)	! IN:Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)	! IN:Average air temp this month (deg.C)
	INTEGER MEND				! IN:Number of timesteps from prev.harv.-harv.
	INTEGER N_STEPS				! IN:No.timesteps in a year
	INTEGER NFERT				! IN:No.fertiliser applications to this crop
	INTEGER NORGM				! IN:No.manure applications to this crop
	INTEGER NSOIL				! IN:Soil code number
	INTEGER NSOW				! IN:Timesteps from prev.harv. to sowing date
	INTEGER NXYEARS				! IN:No.growing seasons simulated
	REAL PNREQ(MAXLU)			! IN:Amount of N required by plant each yr
	                            !        (kgN/ha/year)
	REAL REFIL(MAXLAYER)		! IN:Water needed to fill to FC (mm/lay)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
      REAL RPMCTON(MAXLAYER)		! OUT:RPM C:N ratio
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL SATWATCONT(MAXLAYER)	    ! IN:Total water content at saturation (mm/layer)
	REAL SECONDS				! IN:Number of seconds in one timestep
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOILPH(MAXLAYER)		! IN/OUT:pH of soil in this layer
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL SOMDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN:depth of soil organic matter layers
	REAL TAM(MAXLAYER)			! IN:15N/N in the ammonium in the layer
	REAL TFERT(MAXFERT,3)		! IN:Prop.NO3,NH4,urea in fertiliser
	REAL TFERTNF(3)				! INOUT:Prop.NO3,NH4,urea in fertiliser
	REAL THISFERT				! IN:N input by fertiliser (kgN/ha)
	REAL THISFERT15				! IN:N15 input by fertiliser (kgN15/ha)
	REAL TIN(MAXLAYER)			! IN:15N/N in the nitrate in the layer
	REAL TOTPIC(MAXLU)			! IN:Plant C input calculated using 
	REAL VOLAT					! IN:N lost by volatilisation (kgN/ha)
	REAL VOLAT15				! IN:N15 lost by volatilisation (kgN15/ha)
	REAL WILTPOINT(MAXLAYER)	! IN:Water content at wilting point (mm/layer)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
C
C Required for Monte Carlo - 
C
     	REAL DPMCARB(MAXLAYER)
	REAL RPMCARB(MAXLAYER)
	REAL BCARB(MAXLAYER)
	REAL HCARB(MAXLAYER)
	REAL WRATEDM
	REAL TRATEM
	REAL PLANTUP(MAXLAYER)		! IN(RUN2_GIS_CROP)/OUT(TEST1_RES):Plant uptake from each soil layer
	REAl DDAYS

C
C Save current state of the soil 
C
      ISAVE=1 ! Save soil 
      CALL SAVE_GIS_SOIL(ISAVE,LU1,
     &                   SOILN,SOIL15,AMMN,AMMN15,
     &                   SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
     &                   DPMCARB0,DPMNIT0,DPMNLAB0,
     &                   RPMCARB0,RPMNIT0,RPMNLAB0,
     &                   BCARB0,BNIT0,BNLAB0,
     &                   HCARB0,HNIT0,HNLAB0,
     &                   IOM,PI_CEQ_MON)
C
C Calculate C at start of simulation
C
	STARTC=0
	MEASLAY=SOMDEPTH(ISERIES,LU1,NSOMLAY(ISERIES,LU1))
	MEASLAY=MEASLAY*MAXLAYER1/MAXDEPTH
	STARTPIC=TOTPIC(LU1)
	DO 600 IL=1,MEASLAY
        STARTC=STARTC+BCARB0(IL)+HCARB0(IL)+DPMCARB0(IL)+RPMCARB0(IL)
600   CONTINUE
C******************************
C Full ECOSSE simulation using
C 1. 1 x PIeq
C 2. 0.9 x PIeq
C
      MSS=1
400   CONTINUE
C
C For each growing season... 
C
      IYEAR=1
      DO WHILE (IYEAR.LE.NXYEARS)
        THISLU=LU1
C
C Initialise soil and crop in first year
C	
        IF(IYEAR.EQ.1)THEN
C
C ...initialise the soil for initial land use,
C
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
	    ENDIF
C
C Adjust plant input according to modifier
C
	  TOTPIC(THISLU)=MSS*STARTPIC
c          CALL PARTITION_PLANTINPUT(C_TS,RNIN,RNIN15,
c     &                                 PI_C,PI_N,PI_N15)
           CALL ADD_PLANT_DIST(C_TS,RNIN,RNIN15,PI_CEQ_MON,IK,
     &                         PI_C,PI_N,PI_N15)
C
C Initialise current crop
C 
        CALL INIT_GIS_CROP(IYEAR,THISLU,SECONDS,
     &                     TOTPIC(THISLU),PNREQ(THISLU),
     &                     N_STEPS,ISOWN,NSOW,IHARV,LHARV,IANTHES,
     &                     MEND,IEND,NFERT,IFERT,ILAB,FERT,TFERT,
     &                     JSTOP,CTOT)
C
C For each week till end of growing season...
C           
        DO 100 IK=1,MEND
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
	    RAIN=AVERAIN(IK)
	    AIRTEMP=AVETEMP(IK)
	    EVAPW=AVEPET(IK)
          DO 200 IL=1,MAXLAYER1
	      SOILTEMP(IL)=AIRTEMP
200       CONTINUE
C
C ...Calculate crop C and N returns and N offtake 
C
	    CALL RUN1_GIS_CROP(IYEAR,IK,THISLU,IEND,MEND,NSOW,
     &                          N_STEPS,IS_TS,TOTPIC(THISLU),
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
C          CALL PARTITION_PLANTINPUT(C_TS,RNIN,RNIN15,
C     &                                 PI_C,PI_N,PI_N15)
          CALL ADD_PLANT_DIST(C_TS,RNIN,RNIN15,PI_CEQ_MON,IK,
     &                        PI_C,PI_N,PI_N15) 
          IF(INMODEL.EQ.INPASSCN)
     &      CALL GET_CTON_FOEREID(DPMCTON,RPMCTON,THISLU)
	    CALL SET_DPMRPMRATIO(THISLU,DRRAT)
c	    IF(DOMSOILISIMP(ISERIES,THISLU))
c     &      CALL LIMIT_DRAIN(DOMSOILIMPDEPTH(ISERIES,THISLU),SOILW)
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
          CALL RUN_SUNDIAL_WATER_SUBMONTHLY(RAIN,WMAX,
     &                                      SOILW,WSAT,FLOWPROP,
     &                                      DRAINW,REFIL,EVAPW,
     &                                      SOILTEMP(1),ICOVER)
          CALL PHYSICAL_SUNDIAL_SOILCN(DOCMODEL,SECONDS,
     &                                   NSOIL,SOILTEMP,
     &                                   WMAX,SOILW,DRAINW,REFIL,
     &                                   SLEACH,SLEA15,
     &                                   WLEACH,WLEACH15,CONC,CONC15,
     &                                   SOILN,SOIL15,AMMN,AMMN15,
     &                                   TIN,TAM,MOBDOC,MOBDON,
     &                                   LEACHDOC,LEACHDON)
100     CONTINUE
        IYEAR=IYEAR+1
C
C Go back and set up the next crop
C
        END DO
C End of ECOSSE simulation
C******************************
C
C Save C at end and retrieve soil characteristics from start of simulation
C
	ENDC=0
      DO 500 IL=1,MEASLAY
        ENDC=ENDC+BCARB0(IL)+HCARB0(IL)+DPMCARB0(IL)+RPMCARB0(IL)
500   CONTINUE
      ISAVE=0 ! Retrieve soil 
      CALL SAVE_GIS_SOIL(ISAVE,LU1,
     &                   SOILN,SOIL15,AMMN,AMMN15,
     &                   SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
     &                   DPMCARB0,DPMNIT0,DPMNLAB0,
     &                   RPMCARB0,RPMNIT0,RPMNLAB0,
     &                   BCARB0,BNIT0,BNLAB0,
     &                   HCARB0,HNIT0,HNLAB0,
     &                   IOM,PI_CEQ_MON)
C
C After run using MSS=1, calculate change in C (DC1) and loop back 
C to redo calculation using modifier MSS=0.9
C
      IF(MSS.GT.0.95)THEN
        DC1=ENDC-STARTC
        MSS=0.9
	  GOTO 400
C
C After run using MSS=0.9, calculate change in C using modifier MSS=0.9 (DC2)
C and calculate adjustment to plant inputs according to PI = mss x PIeq
C where MSS = 1-(DC1 x (1-0.9))/(DC1-DC2)
C       DC1 = change in carbon using equilibrium run plant inputs
C   and DC2 = change in carbon using 0.9 x equilibrium run plant inputs
C
      ELSEIF(MSS.LT.0.95)THEN
        DC2=ENDC-STARTC
	  MSS=(1-(DC1*(1-0.9))/(DC1-DC2))
	ENDIF
C
C Adjust plant inputs using calculated multiplier
C
	  TOTPIC(THISLU)=MSS*STARTPIC
     	  CALL GET_PLANT_DIST(TOTPIC(THISLU),PI_CEQ_MON,THISLU)
C
C Leave GET_MODPI_UNLAYERED
C
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE GET_NLIM(AMMN,AMMN15,ATM,AVEPET,AVERAIN,AVETEMP,
     &                BCARB0,BNIT0,BNLAB0,CH4MODEL,CNMODEL,CTOT,
     &			    CLOSSX,CLOSSX15,CRIT,DMODEL,DOCMODEL,DPMCARB0,
     &			    DPMNIT0,DPMNLAB0,DRAINW,FERT,FIXN,FLOWPROP,
     &                FYMFERT,FYMFERT15,HCARB0,HNIT0,HNLAB0,IANTHES,
     &                IEND,IFERT,IHARV,ILAB,IMFUNC,ICMODEL,INMODEL,
     &                IOLAB,IOM,IORGM,ISERIES,ISOWN,ITFUNC,IVOL,
     &                JORGM,JSTOP,LHARV,LU1,MEND,
     &                NFERT,NORGM,N_STEPS,NSOIL,NSOW,ORGMA,
     &                PNREQ,REFIL,RPMCARB0,RPMNIT0,RPMNLAB0,SECONDS,
     &                SOIL15,SOILN,SOILW,TAM,TFERT,THISFERT,
     &                THISFERT15,TIN,TOC,TOTPIC,VOLAT,VOLAT15,WLEACH,
     &                WLEACH15,WSAT,WMAX,SOILPH,PI_CEQ_MON)		  
C
C Initialise SUNDIAL crop from limited information available at large scale
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
	INTEGER MAXITER1			! Maximum allowed number of iterations
      DATA MAXITER1 /1000/
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXORGM				! Max.no.of organic manure applications allowed
	PARAMETER (MAXORGM=52)
	INTEGER MAXSEQ				! Max.no of sequences included in GIS
	PARAMETER (MAXSEQ=MAXLU*MAXLU)						
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXSOIL				! Max.no.of soil types
	PARAMETER (MAXSOIL=50)

      REAL FNULL					
      REAL FNULLARRAY(MAXLAYER)					
	INTEGER IYEAR				! Current growing season number - set to 1
	INTEGER IK					! No.timesteps from prev.harv. to current
      INTEGER IL					! Local layer counter
	INTEGER INULL				
	INTEGER ISATEQ				! Code to indicate whether equilibrium 
								! has been found 0=No 1=Yes
	INTEGER ISTEST				! Output balance sheet files 
								! 0 = Do not output test files
	                            ! 1 = Output balance sheet files (BALANCE_N.OUT, BALANCE_C.OUT & SOILW.OUT)
	INTEGER NITER				! No.of iterations
      INTEGER TESTLU				! LU to be output

	INTEGER ICFIXED				! Initialisation of C is fixed
	INTEGER ICROTHCEQ			! Initialisation of C by RothC 
								!  equilibrium run
	DATA ICFIXED,ICROTHCEQ /1,2/
	INTEGER INSTABLECN			! Initisalisation of N by assuming 
								!  steady state (after Bradbury)
	INTEGER INPASSCN			! Initialisation of N by passing the
								!  C:N ratio of DPM and RPM
	DATA INSTABLECN,INPASSCN /1,2/
	INTEGER MEASLAY				! Layer that soil is measured to

C
C Variables passed to/from other subroutines
C
	REAL AIRTEMP				! IN: Air temperature (deg.C/timestep)	
	REAL ATM15					! IN:Atmospheric N15 input (kgN15/ha/time)
      REAL AVEPET(12)				! IN:Long term average PET (mm/month)
      REAL AVERAIN(12)			! IN:Long term average rainfall(mm/month)
      REAL AVETEMP(12)			! IN:Long term average monthly average 
      REAL BALANCE(20)			! IN:C Results
	REAL BULKDENS(MAXLAYER)				! Soil bulk density (g cm-3)
	REAL C_TS					! IN:Litter C input in this timestep (kgC/ha)
	REAL CACT					! IN:Crop N uptake (kgN/ha/timestep)
	REAL CACT15					! IN:Crop N15 uptake (kgN15/ha/timestep)
	REAL CACTOT					! IN:N taken up by crop (kgN/ha)
	REAL CATOT15				! IN:N15 taken up by crop (kgN15/ha)
	REAL CH4(MAXLAYER)			! IN:CH4 emitted (kgC/ha/layer) (Aitkenhead CH4 model)
	REAL CH4TOAIR				! OUT:CH4 released to atmosphere (kgC/ha) (Aitkenhead CH4 model)
	REAL CH4_PROD(MAXLAYER)  	! OUT:CH4 produced in each soil layer [kgC/ha/layer] (Richards CH4 model)
	REAL CH4_SOIL_OX			! OUT:CH4 oxidised from methane produced in the soil [kgC/ha] (Richards CH4 model)
      REAL CH4_ATMOS_OX			! OUT:CH4 oxidised from atmospheric CH4 [kgC/ha] (Richards CH4 model)
      REAL CH4_FLUX				! OUT:Net CH4 flux from soil column [kgC/ha] (Richards CH4 model)
	REAL CO2(MAXLAYER)			! IN:CO2 emitted (kgC/ha/layer)
	REAL CO2FROMDOC(MAXLAYER)	! IN:CO2 output from DOC in each layer(kgC/ha)
	REAL CONC					! IN:Conc.N leached (kgN/ha/mm/layer)
	REAL CONC15					! IN:Conc.N15 leached (kgN/ha/mm/layer)
	REAL DENITRIFN(MAXLAYER)	! IN/OUT: N lost by denitrification (kgN/ha/layer
	REAL DENIT					! IN:N lost by denitrification (kgN/ha)
	REAL DN15					! IN:N15 lost by denitrification (kgN15/ha)
	REAL DNIT(MAXLAYER)			! IN:Net Mineralised N (kgN/ha/lay/timestep)
	REAL DNIT15(MAXLAYER)		! IN:Net Mineralised N15 (kgN15/ha/lay/time)
	REAL ENDTOTC				! IN:Total organic C at end (kgC/ha)
	REAL EVAPW					! IN: Potential evap. (mm/timestep)
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
      INTEGER ICOVER				! IN:Crop cover 1=Covered 0=Bare
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
	INTEGER IS_TS				! IN:Crop counter
      INTEGER ISAVE				! OUT:Code to save or retrieve variables
								!		0=save; 1=retrieve
	INTEGER ISERIES				! IN:Counter for soil series
	REAL LTA_AWC(12,MAXLAYER)	! IN:Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)	! IN:Average air temp this month (deg.C)
      INTEGER LUARRAY(MAXSOIL)	! IN/OUT:Land use before equilibrium 
	REAL MOBDOC(MAXLAYER)		! IN:Mobile DOC in each layer (kgC/ha)
	REAL MOBDON(MAXLAYER)		! IN:Mobile DON in each layer (kgN/ha)
	REAL OCROPN					! IN:N uptake of crop (kgN/ha) 
	REAL ORGMA(MAXORGM)			! IN:Amount of manure applied (t/ha)
	REAL ORGMANF				! IN:Amount of manure applied (t/ha)
	REAL ORGN					! IN:Total org.N input (kgN/ha) 
	REAL PI_CEQ_MON(12,MAXLAYER)! IN/OUT: Equilibrium plant C input each month
								!         in each layer (kgC/ha/month/layer)
      REAL PI_C(MAXLAYER)			! IN: Plant input C to soil (kgC/ha/step)
      REAL PI_N(MAXLAYER)			! IN: Plant input N to soil (kgC/ha/step)
      REAL PI_N15(MAXLAYER)		! IN: Plant input N15 to soil (kgC/ha/step)
	REAL RAIN					! IN: Rainfall (mm/timestep)
	REAL RNIN					! IN:Litter N input in this timestep (kgN/ha)
	REAL RNIN15					! IN:Litter N15 input in timestep (kgN15/ha)
	REAL SEEDIN					! IN: Data used by MAGEC
	REAL SLEA15(MAXLAYER)		! IN:Nitrate-N15 leached (kgN15/ha)
	REAL SLEACH(MAXLAYER)		! IN:Nitrate-N leached (kgN/ha)
	REAL SOILTEMP(MAXLAYER)		! IN: Soil temperature (deg.C/timestep)
	REAL SRNIN					! IN:Total litter N input (kgN/ha)
	REAL STARTTOTC				! IN:Total organic C at start (kgC/ha)
	REAL SX						! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL SXORGN					! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL T15					! IN:Net Mineralised N (kgN/ha/layer)
	REAL TACTOT					! IN:Total crop N uptake (kgN/ha)
	REAL TC						! IN:Total org.N15 input (kgN/ha) DONT PASS?
	REAL TCINP					! IN:Total org.C input (kgC/ha) DONT PASS?
	REAL TDNIT					! IN:Net Mineralised N (kgN/ha/layer)
	REAL TFYMC					! IN:Total FYM C input (kgC/ha) DONT PASS?
	REAL TOC(MAXLAYER)			! IN:Measured total organic C (kgC/ha/layer)
	REAL TRNIN					! IN:Total litter N input (kgN/ha)
	REAL WLEACH					! IN:N leached (kgN/ha)
	REAL WLEACH15				! IN:N15 leached (kgN15/ha)
	REAL WR						! IN:Root requirement (kgN/ha)
	REAL XORGN					! IN:Total org.N input (kgN/ha) 

	REAL HY(MAXSOIL,MAXLAYER)	! IN/OUT:Stable N:C ratio of BIO & HUM pools
	REAL BPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL HPART(MAXSOIL,MAXLAYER)	! IN/OUT:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL BPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from biomass decompositn
	REAL HPROP(MAXSOIL,MAXLAYER)	! IN/OUT: BIO/HUM from humus decompositn
C
C Soil factors
C
	REAL ALPHA(MAXLAYER)			! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL BETA(MAXLAYER)				! IN/OUT: Prop.HUM produced on BIO decompn
	REAL GAMMA(MAXLAYER)			! IN/OUT: Prop.BIO produced on HUM decompn 
	REAL DELTA(MAXLAYER)			! IN/OUT: Prop.HUM produced on HUM decompn 
	REAL HZ1(MAXLAYER)				! IN/OUT: N:C ratio for steady state
	REAL BIOP(MAXSOIL,MAXLAYER)		! IN/OUT: BIO/TOC 
	REAL DPMRATE(MAXLAYER)			! IN/OUT: Rate constant for DPM decompn
	REAL DRRAT						! IN:DPM:RPM ratio
	REAL RPMRATE(MAXLAYER)			! IN/OUT: Rate constant for RPM decompn
	REAL BRATE(MAXLAYER)			! IN/OUT: Rate constant for HUM decompn
	REAL HRATE(MAXLAYER)			! IN/OUT: Rate constant for BIO decompn
      REAL CONVER_F					! IN/OUT:Conversion to correct timestep
      INTEGER NUMSOIL					! Number of soils defined
	REAL BIORATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for BIO decompn/yr
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for HUM decompn/yr
	REAL IOMARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium IOM (kgC/ha/layer)
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium TOC (kgC/ha/layer)
  	REAL ICFACTOR(MAXLAYER)			! IN:Adjustment in rate needed to	
									!    achieve measured NPP and TOC	
  	REAL ICNULL(MAXLAYER)			! IN:Adjustment in rate needed to	
									!    achieve measured NPP and TOC - not used	
	CHARACTER*40 SNAME(MAXSOIL)		! IN/OUT:Soil name
      REAL FPROPX(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
	REAL FYMCX(MAXORGM)			! IN/OUT:Proportion of C in FYM
	REAL FYMNX(MAXORGM)			! IN/OUT:Proportion of Organic N in FYM
	REAL FYMAX(MAXORGM)			! IN/OUT:Proportion of ammonium-N in FYM
	REAL FYMWX(MAXORGM)			! IN/OUT:Amount of water added in FYM
	REAL FYMLOSS(MAXORGM)		! IN/OUT:Prop.FYM lost by volatn/timestep
	REAL FYMPOS(MAXORGM)		! IN/OUT:Proportion of FYM added to top 25cm
      REAL FYMSTART(MAXORGM)		! IN/OUT:Amount of rainfall in 1 week, 
								!		 below which volatilisation will occur
	REAL BYRAIN					! IN/OUT:Rain lost by bypass flow mm/timestep	
	REAL ANIT					! Ammonium immobilised (kgN/ha)
	REAL ANIT15					! Ammonium N15 nitrified (kgN/ha)
	REAL FANIT					! Fertiliser N nitrified (kgN/ha)
	REAL FANIT15				! Fertiliser N15 nitrified (kgN/ha)
	REAL DFACT(MAXLAYER)		! IN/OUT: Denitrification factor
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: pH above which decomp.max.
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: % clay content in this layer
C
C Variables passed from calling routines
C
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
      REAL ATM					! IN:Atmospheric N input (kgN/ha/timestep)
	REAL BCARB0(MAXLAYER)		! IN:C in soil biomass (kgC/ha/layer)
      REAL BNIT0(MAXLAYER)		! IN:N in soil biomass (kgN/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN:N15 in soil humus (kgN15/ha/layer)
	INTEGER CH4MODEL			! OUT:Methane model (off, Richards or Aitkenhead model) 
	REAL CLOSSX					! IN:N lost by senescence (kgN/ha) 
	REAL CLOSSX15				! IN:N15 lost by senescence (kgN15/ha) 	
	INTEGER CNMODEL				! IN:Type of Model for C:N ratio of DPM and 
	                            !     RPM (MAGEC or Foereid)
	REAL CRIT(MAXSOIL,MAXLAYER)	! IN:Critical N content (kgN/ha/layer)
	REAL CTOT					! IN:Total litter C input (kgC/ha)
	REAL CUPTN					! IN:Crop N uptake (kgN/ha/timestep)
	INTEGER DMODEL				! IN:Type denitrification model chosen 
	                            ! (Bradbury or NEMIS)
	INTEGER DOCMODEL			! IN:DOC model (on or off)
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
      REAL DPMCTON(MAXLAYER)		! OUT:DPM C:N ratio
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL DRAINW(MAXLAYER)		! IN:Water drained from this layer (mm/layer)
	REAL FERT(MAXFERT)			! IN:Amount of fertiliser applied (kgN/ha)
	REAL FERTNF					! INOUT:Amount of fertiliser applied (kgN/ha)
	REAL FIXN					! IN:N fixed from atmosphere (kgN/ha)
	REAL FLOWPROP				! IN:Proportion of flow needed to achieve
							    !	     water table at depth WTABLE
	REAL FYMFERT				! IN:N input by fertiliser & FYM (kgN/ha)
	REAL FYMFERT15				! IN:N15 input by fertiliser & FYM (kgN15/ha)
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN:N in soil humus (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN:N15 in soil biomass (kgN15/ha/layer)
	INTEGER IANTHES				! IN:No.timesteps from 01/01 to anthesis 
	INTEGER ICMODEL				! OUT:Type of C initialisation model 
	INTEGER IEND				! IN:Timesteps from sowing date to harvest
	INTEGER IFERT(MAXFERT)		! IN:No.timesteps to split fert.application
	INTEGER IHARV				! IN:Timesteps from 01/01/01 to harvest date 
	INTEGER ILAB(MAXFERT)		! IN:Labelling on fertiliser? 0=No 1=Yes
	INTEGER ILABNF				! INOUT:Labelling on fertiliser? 0=No 1=Yes
	INTEGER IMFUNC				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
	INTEGER INMODEL				! OUT:Type of N initialisation model
	INTEGER IOLAB(MAXORGM)		! IN:Labelling on manure? 0=No 1=Yes
	INTEGER IOLABNF				! IN:Labelling on manure? 0=No 1=Yes
	INTEGER IORGM(MAXORGM)		! IN:No.timesteps to split manure application
	INTEGER ISOWN				! IN:Timesteps from 01/01/01 to sowing date 
	INTEGER ITFUNC				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IVOL(MAXFERT)		! IN:N volatilised from fertiliser (kgN/ha)
	INTEGER IVOLNF				! INOUT:N volatilised from fertiliser (kgN/ha)
	INTEGER JORGMNF				! IN:Type of manure application
	INTEGER JORGM(MAXORGM)		! IN:Type of manure application
	INTEGER JSTOP				! IN:Code to stop simulation 0=No 1=Yes
	REAL LEACHDOC(MAXLAYER)		! IN/OUT:DOC leached from the layer (kgC/ha/layer/timestep)
	REAL LEACHDON(MAXLAYER)		! IN/OUT:DON leached from the layer (kgN/ha/layer/timestep)
	INTEGER LHARV				! IN:No.timesteps from 01/01 to prev. harvest
      INTEGER LU1					! IN:Land use code
	INTEGER MEND				! IN:Number of timesteps from prev.harv.-harv.
	REAL NITRIFN(MAXLAYER)		! OUT:Nitrification from this layer (kg/ha/layer/timestep)
	INTEGER N_STEPS				! IN:No.timesteps in a year
	INTEGER NFERT				! IN:No.fertiliser applications to this crop
	INTEGER NORGM				! IN:No.manure applications to this crop
	INTEGER NSOIL				! IN:Soil code number
	INTEGER NSOW				! IN:Timesteps from prev.harv. to sowing date
	INTEGER NXYEARS				! IN:No.growing seasons simulated
	REAL PNREQ(MAXLU)			! IN:Amount of N required by plant each yr
	                            !        (kgN/ha/year)
	REAL REFIL(MAXLAYER)		! IN:Water needed to fill to FC (mm/lay)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
      REAL RPMCTON(MAXLAYER)		! OUT:RPM C:N ratio
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL SATWATCONT(MAXLAYER)	    ! IN:Total water content at saturation (mm/layer)
	REAL SECONDS				! IN:Number of seconds in one timestep
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOILPH(MAXLAYER)		! IN/OUT:pH of soil in this layer
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	REAL TAM(MAXLAYER)			! IN:15N/N in the ammonium in the layer
	REAL TFERT(MAXFERT,3)		! IN:Prop.NO3,NH4,urea in fertiliser
	REAL TFERTNF(3)				! INOUT:Prop.NO3,NH4,urea in fertiliser
	REAL THISFERT				! IN:N input by fertiliser (kgN/ha)
	REAL THISFERT15				! IN:N15 input by fertiliser (kgN15/ha)
	REAL TIN(MAXLAYER)			! IN:15N/N in the nitrate in the layer
	REAL TOTPIC(MAXLU)			! IN:Plant C input calculated using 
	REAL VOLAT					! IN:N lost by volatilisation (kgN/ha)
	REAL VOLAT15				! IN:N15 lost by volatilisation (kgN15/ha)
	REAL VOLATN(MAXLAYER)		! ON/OUT:Volatilisation from this layer (kg/ha/layer/timestep)
	REAL WILTPOINT(MAXLAYER)	! IN:Water content at wilting point (mm/layer)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
C
C Required for Monte Carlo - 
C
     	REAL DPMCARB(MAXLAYER)
	REAL RPMCARB(MAXLAYER)
	REAL BCARB(MAXLAYER)
	REAL HCARB(MAXLAYER)
	REAL WRATEDM
	REAL TRATEM
	REAL PLANTUP(MAXLAYER)		! IN(RUN2_GIS_CROP)/OUT(TEST1_RES):Plant uptake from each soil layer
	REAL DDAYS

C
C Set results recording
C
      NXYEARS=10
      ISTEST=0	! Set to 1 to output results in C_BALANCE.OUT and N_BALANCE.OUT
	TESTLU=1	! Set to code number of LU that you want to look at
	ISATEQ=0
C
C Save current state of the soil 
C
      ISAVE=1 ! Save soil 
      CALL SAVE_GIS_SOIL(ISAVE,LU1,
     &                   SOILN,SOIL15,AMMN,AMMN15,
     &                   SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
     &                   DPMCARB0,DPMNIT0,DPMNLAB0,
     &                   RPMCARB0,RPMNIT0,RPMNLAB0,
     &                   BCARB0,BNIT0,BNLAB0,
     &                   HCARB0,HNIT0,HNLAB0,
     &                   IOM,PI_CEQ_MON)
      CALL SAVE_GIS_RATES(ISAVE,LU1,ICFACTOR)
C
C Get the total organic C at start
C
      STARTTOTC=0
      DO 100 IL=1,MAXLAYER1
	  STARTTOTC=STARTTOTC+TOC(IL)-IOM(IL)
100   CONTINUE	 
C
C Loop to calculate new ICFACTOR according to N limitation etc
C
      NITER=0
101   CONTINUE
C
C Initialise current crop
C 
      DO 200 IYEAR=1,NXYEARS
        CALL INIT_GIS_CROP(IYEAR,LU1,SECONDS,
     &                   TOTPIC(LU1),PNREQ(LU1),
     &                   N_STEPS,ISOWN,NSOW,IHARV,LHARV,IANTHES,
     &                   MEND,IEND,NFERT,IFERT,ILAB,FERT,TFERT,
     &                   JSTOP,CTOT)
C
C For each timestep in year ...
C
             
        DO 300 IK=1,MEND
C
C Initialize THISFERT=current weeks fertilizer addition
C            WLEACH=current weeks water leaching
C      
          CALL SETWEEKVARS(JSTOP,CLOSSX,CLOSSX15,FYMFERT,FYMFERT15,
     &                   FIXN,THISFERT,THISFERT15,WLEACH,WLEACH15,
     &                   VOLAT,VOLAT15)
C
C Get weather data
C
C
C Set soil temperature to air temperature
C
          CALL GETWEATHER_GIS(IK,LHARV,AVERAIN,AVEPET,AVETEMP, 
     &                      RAIN,EVAPW,SOILTEMP,AIRTEMP)
C
C If current week is first in the growing season set results
C
          CALL SETRES(BALANCE,CO2,SX,SXORGN)
C
C ...Calculate crop C and N returns and N offtake 
C
	    CALL RUN1_GIS_CROP(IYEAR,IK,LU1,IEND,MEND,NSOW,
     &                     N_STEPS,IS_TS,TOTPIC(LU1),
     &                     PNREQ(LU1),CACTOT,CATOT15,
     &                     ICOVER,SXORGN,ORGN,
     &                     TRNIN,RNIN,RNIN15,C_TS,TCINP,WR)      
C
C ...Add stuff to soil
C
          CALL GETFYM(IK,NORGM,IORGM,ORGMA,JORGM,IOLAB,
     &              ORGMANF,JORGMNF,IOLABNF) 
          CALL GETFERT(IK,NFERT,IFERT,FERT,TFERT,ILAB,IVOL,
     &               FERTNF,TFERTNF,ILABNF,IVOLNF) 
          CALL ADD_SUNDIAL_SOILCN(
C INPUTS: time, environment, organic manure, fertiliser, atmospheric, 
     &                          SECONDS,
     &                          NSOIL,RAIN,
     &                          ORGMANF,JORGMNF,IOLABNF,      
     &  					      FERTNF,TFERTNF,ILABNF,IVOLNF,
     &                          ATM,								  
C OUTPUTS: organic manure, atmospheric
     &                          FYMFERT,FYMFERT15,TFYMC,			  
     &                          VOLAT,VOLAT15,					  
C INPUTS/OUTPUTS: Soil C&N
     &                          SOILN,SOIL15,AMMN,AMMN15,			  
     &                          DPMCARB0,DPMNIT0,DPMNLAB0,		  
     &                          RPMCARB0,RPMNIT0,RPMNLAB0,		  
     &                          HCARB0,HNIT0,HNLAB0,TIN,TAM)		  
C
C Run SUNDIAL Soil C and N routines
C
          CALL PARTITION_PLANTINPUT(C_TS,RNIN,RNIN15,
     &                            PI_C,PI_N,PI_N15)
          IF(INMODEL.EQ.INPASSCN)
     &      CALL GET_CTON_FOEREID(DPMCTON,RPMCTON,LU1)
	    CALL SET_DPMRPMRATIO(LU1,DRRAT)
	    MEASLAY=MAXLAYER
          CALL MICROBIAL_SUNDIAL_SOILCN(DMODEL,DOCMODEL,CNMODEL,
     &                                ITFUNC,IMFUNC,CH4MODEL,INMODEL,
     &                                SECONDS,ICOVER,
     &                                NSOIL,SOILTEMP,RAIN,
     & 							    PI_C,PI_N,PI_N15,DRRAT,                                                                     
     &                                WMAX,SOILW,WSAT,
     &                                WILTPOINT,SATWATCONT,
     &                                THISFERT,THISFERT15,
     &                                CO2,VOLAT,VOLAT15,VOLATN,
     &                                CH4,CH4_PROD,CH4_SOIL_OX,
     &							      CH4_ATMOS_OX,CH4_FLUX,
     &                                DENITRIFN,NITRIFN,
     &                                DENIT,DN15,GNN2O,GNNO,GPNN2O,	
     &                                G15NN2O,G15NNO,G15PNN2O,		
     &                                GDN2,GDN2O,G15DN2,G15DN2O,		
     &                                DNIT,DNIT15,TDNIT,T15,
     &                                WLEACH,WLEACH15,
     &                                SOILN,SOIL15,AMMN,AMMN15,
     &                                DPMCARB0,DPMNIT0,DPMNLAB0,
     &                                RPMCARB0,RPMNIT0,RPMNLAB0,
     &                                BCARB0,BNIT0,BNLAB0,BULKDENS,
     &                                TOC,CH4TOAIR,HCARB0,HNIT0,HNLAB0,
     &                                TIN,TAM,MOBDOC,MOBDON,CO2FROMDOC,
     &                                ICFACTOR,DPMCTON,RPMCTON,
C Required for Monte Carlo - 
     &					DPMCARB,RPMCARB,BCARB,HCARB,
     &					WRATEDM,TRATEM,SOILPH,MEASLAY)
C
C Run within season crop routines
C
          CALL RUN2_GIS_CROP(SOILTEMP(1),AMMN,AMMN15,ATM,ATM15,
     &                     CACT,CACT15,CACTOT,CATOT15,CLOSSX,
     &                     CLOSSX15,PNREQ(LU1),CRIT,CTOT,
     &                     CUPTN,IK,IS_TS,IYEAR,JSTOP,MEND,NSOIL,
     &                     NSOW,OCROPN,ORGN,RNIN,SOIL15,SOILN,SRNIN,
     &                     SXORGN,TACTOT,TAM,TC,TIN,TRNIN,VOLAT,
     &                     VOLAT15,XORGN,SEEDIN,PLANTUP,DDAYS)	
C
C Leaching routines
C
          CALL RUN_SUNDIAL_WATER_SUBMONTHLY(RAIN,WMAX,
     &                                      SOILW,WSAT,FLOWPROP,
     &                                      DRAINW,REFIL,EVAPW,
     &                                      SOILTEMP(1),ICOVER)
          CALL PHYSICAL_SUNDIAL_SOILCN(DOCMODEL,SECONDS,
     &                               NSOIL,SOILTEMP,
     &                               WMAX,SOILW,DRAINW,REFIL,
     &                               SLEACH,SLEA15,
     &                               WLEACH,WLEACH15,CONC,CONC15,
     &                               SOILN,SOIL15,AMMN,AMMN15,
     &                               TIN,TAM,MOBDOC,MOBDON,
     &                               LEACHDOC,LEACHDON)
300     CONTINUE
200   CONTINUE
C
C Get the total organic C at end
C
      CALL GET_TOTC(DPMCARB0,RPMCARB0,BCARB0,HCARB0,ENDTOTC)
C
C Calculate N limiting factor
C ...if equilibrium already found save soil
C
      IF(ISATEQ.EQ.1)THEN
        ISAVE=1 ! Save soil 
        CALL SAVE_GIS_SOIL(ISAVE,LU1,
     &                     SOILN,SOIL15,AMMN,AMMN15,
     &                     SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
     &                     DPMCARB0,DPMNIT0,DPMNLAB0,
     &                     RPMCARB0,RPMNIT0,RPMNLAB0,
     &                     BCARB0,BNIT0,BNLAB0,
     &                     HCARB0,HNIT0,HNLAB0,
     &                     IOM,PI_CEQ_MON)
C
C ...if carbon in soil has changed...
C
	ELSEIF((NITER.LT.MAXITER1).AND.
     &((NINT(STARTTOTC/10)).NE.(NINT(ENDTOTC/10))))THEN
	  NITER=NITER+1
C
C ......retrieve soil parameters
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
C ......calculate change in ICFACTOR to account for N limitation etc,
C
        DO 400 IL=1,MAXLAYER1 
	    IF(ICFACTOR(IL).LT.0.0001)THEN
	      ICFACTOR(IL)=ICFACTOR(IL)+((ENDTOTC-STARTTOTC)/STARTTOTC)
	    ELSE
	      ICFACTOR(IL)=ICFACTOR(IL)*(ENDTOTC/STARTTOTC)
	    ENDIF
          IF(ICFACTOR(IL).LT.0)ICFACTOR(IL)=0.1
400     CONTINUE
C
C ......put revised ICFACTOR into soil routines,
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
C ...and go back an redo calculation with new ICFACTOR
C
        GOTO 101
	ELSEIF(NITER.GE.MAXITER1)THEN
	  ISATEQ=-1
	  PRINT*,'                    Warning - equilibrium not found!'
      ENDIF												
C
C If equilibrium found,...
C
      IF(ISATEQ.EQ.0)THEN
C
C ......save as initial state of the soil
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
        CALL SAVE_GIS_RATES(ISAVE,LU1,ICFACTOR)
C ...and rerun once for duration of main simulation
 	  ISATEQ=1
	  GOTO 101
	ENDIF
C
C Leave GET_NLIM
C
      END
C
C-------------------------------------------------------------
C
	SUBROUTINE GET_SOIL_MULT(BCARB0,BNIT0,BNLAB0,
     &                         DPMCARB0,DPMNIT0,DPMNLAB0,
     &                         HCARB0,HNIT0,HNLAB0,                       
     &                         LU1,MULT,NLAYER,
     &                         RPMCARB0,RPMNIT0,RPMNLAB0)
C
C Subroutine to multiply soil pools by set amount
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)

      INTEGER IL					! Local layer counter
C
C Variables passed to / from this subroutine
C
	REAL BCARB0(MAXLAYER)		! IN(CALL):C in soil biomass (kgC/ha/layer)
      REAL BNIT0(MAXLAYER)		! IN(CALL):N in soil biomass (kgN/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN(CALL):N15 in soil humus (kgN15/ha/layer)
	REAL DPMCARB0(MAXLAYER)		! IN(CALL):C in decomposable PM (kgC/ha/layer)
	REAL DPMNIT0(MAXLAYER)		! IN(CALL):N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN(CALL):N15 in decomposable PM (kgN15/ha/layer)
	REAL HCARB0(MAXLAYER)		! IN(CALL):C in soil humus (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN(CALL):N in soil humus (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN(CALL):N15 in soil biomass (kgN15/ha/layer)
      INTEGER LU1					! IN(CALL):Land use code
	INTEGER NLAYER				! IN(CALL):Layer that soil is multiplied to
	REAL MULT					! OUT(GET_SOIL_PLUS):Multiplication factor for current soil pools
	REAL RPMCARB0(MAXLAYER)		! IN(CALL):C in resistant PM (kgC/ha/layer)
	REAL RPMNIT0(MAXLAYER)		! IN(CALL):N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN(CALL):N15 in resistant PM (kgN15/ha/layer)
C
C Increase by fraction PLUS
C
	DO 100 IL=1,NLAYER
	  DPMCARB0(IL)=DPMCARB0(IL)*MULT
	  RPMCARB0(IL)=RPMCARB0(IL)*MULT
	  BCARB0(IL)=BCARB0(IL)*MULT
	  HCARB0(IL)=HCARB0(IL)*MULT
	  DPMNIT0(IL)=DPMNIT0(IL)*MULT
	  RPMNIT0(IL)=RPMNIT0(IL)*MULT
	  BNIT0(IL)=BNIT0(IL)*MULT
	  HNIT0(IL)=HNIT0(IL)*MULT
	  DPMNLAB0(IL)=DPMNLAB0(IL)*MULT
	  RPMNLAB0(IL)=RPMNLAB0(IL)*MULT
	  BNLAB0(IL)=BNLAB0(IL)*MULT
	  HNLAB0(IL)=HNLAB0(IL)*MULT
100   CONTINUE
      END
c     
c
c
C
C-------------------------------------------------------------
C
      SUBROUTINE INIT_GIS_CROP_LIM(IYEAR,LUCODE,SECONDS,
     &                         PIC_ANN,PLANT_REQN,N_STEPS,
     &                         ISOWN,NSOW,IHARV,LHARV,IANTHES,MEND,IEND,
     &						 NFERT,IFERT,ILAB,FERT,TFERT,JSTOP,CTOT,miscticker)
C
C Initialise SUNDIAL crop from limited information available at large scale
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=10)
	INTEGER MAXFERT				! Max.no.of fertiliser applications allowed 
	PARAMETER(MAXFERT=5)

	INTEGER NF
	INTEGER ILU					! LU counter 
								!(1=arable,2=grassland,3=forestry,4=semi-natural)
C For arable assume C:N ratio for 7t/ha WWheat crop
C Debris C = 1346(1.23 + 1.4(1-exp(-0.24 x Yld)))= 3189 kg/ha
C Straw C = 0.4 x((850 x Yld) x ((1/0.45)-1)) = 2909 kg/ha
C Total C = Debris C + Straw C = 3189+2909 = 6098 kg/ha
C Below ground crop N = 63 x (1-exp(-0.63xYld)) = 62 kg N / ha
C Above ground crop N = 236 x (exp(0.071xYld)-1) = 152 kg N / ha
C Total N requirement = 214 kg N / ha
C C:N ratio of total crop = total C / total N = 6098/214 = 28.5
	REAL CTON(MAXLU),miscticker			! Standard C:N ratio of land use type
	DATA (CTON(ILU),ILU=1,MAXLU) /25,100,100,100,94,100,80,50,100,25/  ! 15/02/13 !MLR decreased SRC 21/02/2014

      INTEGER SOWMONTH(MAXLU)		! Set month of sowing 
	DATA (SOWMONTH(ILU),ILU=1,MAXLU) /2,1,1,1,1,1,2,2,1,2/ 
C**** TEMP CHANGE TO GET RID OF JUMP AT START OF ARABLE >>>
      INTEGER HARVMONTH(MAXLU)	! Set month of harvest 
c	DATA (HARVMONTH(ILU),ILU=1,MAXLU) /8,12,12,12/ CRASHES IF HAVE HARVMONTH 8! SORT OUT!!
	DATA (HARVMONTH(ILU),ILU=1,MAXLU) /7,12,12,12,12,12,9,9,12,7/ 
C**** <<< TEMP CHANGE TO GET RID OF JUMP AT START OF ARABLE >>>
      INTEGER FERTMONTH(MAXLU)	! Set month of fertiliser application month
   
	DATA (FERTMONTH(ILU),ILU=1,MAXLU) /4,4,0,0,5,4,4,4,4,4/ ! 15/02/2013
         
	INTEGER ANTHMONTH(MAXLU)
	DATA (ANTHMONTH(ILU),ILU=1,MAXLU) /7,11,11,11,11,11,11,11,11,7/
C Fraction N lost by senescence = 9 kg N / ha
C Seed N = 4 kg N / ha
C
C Assume plant C is given by 7t/ha WWheat parameters for arable crops, grass (also forestry and natural)
C Stubble C = 1446(1-0.94(exp(-0.175xYld))) = 1047 kg/ha
C Straw C = 0.4 x((850 x Yld) x ((1/0.45)-1)) = 2909 kg/ha
C Total C in tops at harvest = Straw C + Stubble C = 2909 + 1047 = 3956 kg / ha
C Debris C (inc.stubble) = 1346(1.23 + 1.4(1-exp(-0.24 x Yld)))= 3189 kg/ha
C Ratio of debris inc.before harv.to above ground C = (Debris C - Stubble C) / Tot.C at harvest = (3189 - 1047)/3956 = 0.54

c***** Temporary change 138 to correct concepts of limited data run
      REAL DEBTOPLANTC(MAXLU)		! Ratio of C that goes below ground as 
c								! debris before harvest to total C in 
c								! above group crop at harvest
c	DATA (DEBTOPLANTC(ILU),ILU=1,MAXLU) /0.54,0.54,1.00,1.00/		! Temp change 118
c	DATA (DEBTOPLANTC(ILU),ILU=1,MAXLU) /0.54,0.75,1.00,1.00/		! Temp change 144 to correct C change A->G and G->A	
c	DATA (DEBTOPLANTC(ILU),ILU=1,MAXLU) /0.54,1.0,1.0,1.0,0.75,1.0/	! Temp change 146
c      DATA (DEBTOPLANTC(ILU),ILU=1,MAXLU)/1.0,1.0,1.0,1.0,1.0,1.0,
c     &                                                1.0,1.0,1.0/	! Temp change 146 ! MLR: commented out
      ! Fraction of annual plant production that goes into the soil (at equilibrium)
      ! MLR: based on MIAMI_DYCE frac values (extended for new land use types)
      DATA (DEBTOPLANTC(ILU),ILU=1,MAXLU)/0.53, 0.71, 0.88, 0.95, 0.36,
     &                                    0.35, 0.25, 0.75, 0.4, 0.53/
c***** End of temporary change 138 to correct concepts of limited data run
C
	REAL TOPFERT				! Maximum allowed application rate of fert.
	DATA TOPFERT /500/

	
********************************************
* Need data on                             *
* C:N of LU types                          *
* Soil improvement factors                 *
* Sowing month of LU types                 *
* Harvest month of LU types                *
* Fertiliser application month of LU types *
* Timing of anthesis for different LUs     *
********************************************
C
C Variables passed to/from other subroutines
C
	REAL CTOT					! OUT:N uptake to this week (kgN/ha)
	                            !     (kgN/ha/year)
	INTEGER LUCODE				! IN:LU code
								!(1=arable,2=grassland,3=forestry,4=semi-natural)
	REAL PIC_ANN				! IN:Annual plant C input (kgC/ha/year)
	REAL PLANT_REQN				! OUT:Annual plant N requirement (kgN/ha/year)
	REAL SECONDS				! IN:Number of seconds in one timestep

C
C ...Timing factors
C
	REAL CONVER_F				! OUT:Conversion between this timestep & weeks
	INTEGER IANTHES				! OUT:No.timesteps from 01/01 to anthesis 
	INTEGER IEND				! IN:Timesteps from sowing date to harvest
	INTEGER IHARV				! OUT:Timesteps from 01/01/01 to harvest date 
      INTEGER ISAVE				! Code to save or retrieve variables   
	INTEGER ISOWN				! OUT:Timesteps from 01/01/01 to sowing date 
	INTEGER IYEAR				! IN:Current growing season number
	INTEGER JSTOP				! IN/OUT:Code to stop simulation 0=No 1=Yes
	INTEGER LHARV				! OUT:No.timesteps from 01/01 to prev.harv.
	INTEGER MCROP				! IN:Crop code number
	INTEGER MEND				! OUT:Number of timesteps from prev.harv.-harv.
	INTEGER NSOW				! OUT:Timesteps from prev.harv. to sowing date
	INTEGER N_STEPS				! OUT:No.timesteps in a year
      INTEGER N_REMAIN
	REAL TCNEED					! IN/OUT:Total N requirement (kgN/ha)
C
C ...Fertiliser applications
C
	INTEGER IFERT(MAXFERT)		! IN:No.timesteps to split fert.application
	INTEGER ILAB(MAXFERT)		! IN:Labelling on fertiliser? 0=No 1=Yes
	REAL FERT(MAXFERT)			! IN:Amount of fertiliser applied (kgN/ha)
	INTEGER NFERT				! IN:No.fertiliser applications to this crop
	REAL TFERT(MAXFERT,3)		! IN:Prop.NO3,NH4,urea in fertiliser
C
C Set time factors from SECONDS
C
      N_STEPS=(365.25*24*60*60)/SECONDS
	CONVER_F=SECONDS/(7*24*60*60)
	N_REMAIN=10/CONVER_F
C
C Assume N content can be set from a constant C:N ratio for each LU type
C
      MCROP=LUCODE
c***** Temporary change 138 to correct concepts of limited data run
c	PLANT_REQN=(PIC_ANN*(1+DEBTOPLANTC(LUCODE)))/CTON(LUCODE)	! Temp change 146 !MLR Temportary commented out
c	PLANT_REQN=PIC_ANN/CTON(LUCODE)								! Temp change 148 !MLR Temportary uncommented out
	PLANT_REQN=(PIC_ANN*(1/DEBTOPLANTC(LUCODE)))/CTON(LUCODE)	! MLR new versions using MIAMI_DYCE frac values
c***** End of temporary 138 change to correct concepts of limited data run
	CTOT=0
C
C Get start and end of growing season (sowing and harvest dates) from standard values
C
      IF(IYEAR.EQ.1)THEN
	  LHARV=0
      ELSE
	  LHARV=IHARV
	ENDIF
      ISOWN=SOWMONTH(LUCODE)+(IYEAR-1)*12
      NSOW=ISOWN-(IYEAR-1)*12 !LHARV
	IHARV=HARVMONTH(LUCODE)+(IYEAR-1)*12
      MEND=IHARV-(IYEAR-1)*12 !LHARV
      MEND=min(MEND,13)   		!new line to prevent crash with 8 months eoj
      IEND=MEND-NSOW
      IANTHES=ANTHMONTH(LUCODE)+(IYEAR-1)*12
	TCNEED=0
      lharv=miscticker
C
C Fertiliser applications - assume plant N requirement is applied as inorganic fertiliser
c (assume no FYM)
C
      ! Land use types: 1 ara, 2 gra, 3 for, 4 nat, 5 mis, 6 src, 7 sug, 8 osr, 9 srf, 10 whe)
      ! Miscanthus
      if (lucode == 5) then
          if (miscticker > 2) then
              fertmonth(lucode) = 5
              fert(1) = 30
          else
              fertmonth(lucode) = 0
          endif
      ! SRC      
      elseif (lucode == 6) then
          if (miscticker > 2) then
              fertmonth(lucode) = 4
              fert(1) = 60
          else
              fertmonth(lucode) = 0
          endif
      ! SRF
      elseif (lucode == 9) then
          if (miscticker > 2) then
              fertmonth(lucode) = 4
              fert(1) = 45
          else
              fertmonth(lucode) = 0
          endif
      ! All other land use types
      else
          fert(1) = plant_reqn
      endif
c >>> Commented out by MLR (replaced by above code)
c      if ((MOD(lharv,3)).eq.0.and.LUCODE.EQ.6)then !fert for SRC every 3 years
c          FERTMONTH(6)=4
c          !PLANT_REQN=PLANT_REQN*3  MLR commented out for elum
c      endif
c      if (miscticker.LT.3.and.lucode.eq.5)then 
c          FERTMONTH(5)=0
c      endif
c End comment <<<
      IF(FERTMONTH(LUCODE).GT.0)NFERT=1
      DO 100 NF=1,NFERT
          IFERT(NF)=(FERTMONTH(LUCODE))!+(IYEAR-1)*12)!-LHARV
      ILAB(NF)=0
      IF(FERT(NF).GT.TOPFERT)FERT(NF)=TOPFERT
      TFERT(NF,1)=50
      TFERT(NF,2)=50
      TFERT(NF,3)=0
100   CONTINUE
      JSTOP=0
C
C Save internal variables
C
      ISAVE=1
	CALL SAVE_GIS_CROP(ISAVE,N_REMAIN,CONVER_F,MCROP,TCNEED)
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE INIT_GIS_CROP(IYEAR,LUCODE,SECONDS,
     &                         PIC_ANN,PLANT_REQN,N_STEPS,
     &                         ISOWN,NSOW,IHARV,LHARV,IANTHES,MEND,IEND,
     &						 NFERT,IFERT,ILAB,FERT,TFERT,JSTOP,CTOT)
C
C Initialise SUNDIAL crop from limited information available at large scale
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXFERT				! Max.no.of fertiliser applications allowed 
	PARAMETER(MAXFERT=5)

	INTEGER NF
	INTEGER ILU					! LU counter 
								!(1=arable,2=grassland,3=forestry,4=semi-natural)
C For arable assume C:N ratio for 7t/ha WWheat crop
C Debris C = 1346(1.23 + 1.4(1-exp(-0.24 x Yld)))= 3189 kg/ha
C Straw C = 0.4 x((850 x Yld) x ((1/0.45)-1)) = 2909 kg/ha
C Total C = Debris C + Straw C = 3189+2909 = 6098 kg/ha
C Below ground crop N = 63 x (1-exp(-0.63xYld)) = 62 kg N / ha
C Above ground crop N = 236 x (exp(0.071xYld)-1) = 152 kg N / ha
C Total N requirement = 214 kg N / ha
C C:N ratio of total crop = total C / total N = 6098/214 = 28.5
	REAL CTON(MAXLU)			! Standard C:N ratio of land use type
	DATA (CTON(ILU),ILU=1,MAXLU) /25,100,188,100,94,188/  ! 15/02/13
C
      INTEGER SOWMONTH(MAXLU)		! Set month of sowing 
	DATA (SOWMONTH(ILU),ILU=1,MAXLU) /2,1,1,1,1,1/ 
C**** TEMP CHANGE TO GET RID OF JUMP AT START OF ARABLE >>>
      INTEGER HARVMONTH(MAXLU)	! Set month of harvest 
c	DATA (HARVMONTH(ILU),ILU=1,MAXLU) /8,12,12,12/ CRASHES IF HAVE HARVMONTH 8! SORT OUT!!
	DATA (HARVMONTH(ILU),ILU=1,MAXLU) /7,12,12,12,12,12/ 
C**** <<< TEMP CHANGE TO GET RID OF JUMP AT START OF ARABLE >>>
      INTEGER FERTMONTH(MAXLU)	! Set month of fertiliser application month
	DATA (FERTMONTH(ILU),ILU=1,MAXLU) /4,4,0,0,5,4/   ! 15/02/2013
	INTEGER ANTHMONTH(MAXLU)
	DATA (ANTHMONTH(ILU),ILU=1,MAXLU) /7,11,11,11,11,11/
C Fraction N lost by senescence = 9 kg N / ha
C Seed N = 4 kg N / ha
C
C Assume plant C is given by 7t/ha WWheat parameters for arable crops, grass (also forestry and natural)
C Stubble C = 1446(1-0.94(exp(-0.175xYld))) = 1047 kg/ha
C Straw C = 0.4 x((850 x Yld) x ((1/0.45)-1)) = 2909 kg/ha
C Total C in tops at harvest = Straw C + Stubble C = 2909 + 1047 = 3956 kg / ha
C Debris C (inc.stubble) = 1346(1.23 + 1.4(1-exp(-0.24 x Yld)))= 3189 kg/ha
C Ratio of debris inc.before harv.to above ground C = (Debris C - Stubble C) / Tot.C at harvest = (3189 - 1047)/3956 = 0.54

c***** Temporary change 138 to correct concepts of limited data run
      REAL DEBTOPLANTC(MAXLU)		! Ratio of C that goes below ground as 
c								! debris before harvest to total C in 
c								! above group crop at harvest
c	DATA (DEBTOPLANTC(ILU),ILU=1,MAXLU) /0.54,0.54,1.00,1.00/		! Temp change 118
c	DATA (DEBTOPLANTC(ILU),ILU=1,MAXLU) /0.54,0.75,1.00,1.00/		! Temp change 144 to correct C change A->G and G->A	
c	DATA (DEBTOPLANTC(ILU),ILU=1,MAXLU) /0.54,1.0,1.0,1.0,0.75,1.0/	! Temp change 146
	DATA (DEBTOPLANTC(ILU),ILU=1,MAXLU) /1.0,1.0,1.0,1.0,1.0,1.0/	! Temp change 146
c***** End of temporary change 138 to correct concepts of limited data run
C
	REAL TOPFERT				! Maximum allowed application rate of fert.
	DATA TOPFERT /500/

********************************************
* Need data on                             *
* C:N of LU types                          *
* Soil improvement factors                 *
* Sowing month of LU types                 *
* Harvest month of LU types                *
* Fertiliser application month of LU types *
* Timing of anthesis for different LUs     *
********************************************
      
C
C Variables passed to/from other subroutines
C
	REAL CTOT					! OUT:N uptake to this week (kgN/ha)
	                            !     (kgN/ha/year)
	INTEGER LUCODE				! IN:LU code
								!(1=arable,2=grassland,3=forestry,4=semi-natural)
	REAL PIC_ANN				! IN:Annual plant C input (kgC/ha/year)
	REAL PLANT_REQN				! OUT:Annual plant N requirement (kgN/ha/year)
	REAL SECONDS				! IN:Number of seconds in one timestep

C
C ...Timing factors
C
	REAL CONVER_F				! OUT:Conversion between this timestep & weeks
	INTEGER IANTHES				! OUT:No.timesteps from 01/01 to anthesis 
	INTEGER IEND				! IN:Timesteps from sowing date to harvest
	INTEGER IHARV				! OUT:Timesteps from 01/01/01 to harvest date 
      INTEGER ISAVE				! Code to save or retrieve variables   
	INTEGER ISOWN				! OUT:Timesteps from 01/01/01 to sowing date 
	INTEGER IYEAR				! IN:Current growing season number
	INTEGER JSTOP				! IN/OUT:Code to stop simulation 0=No 1=Yes
	INTEGER LHARV				! OUT:No.timesteps from 01/01 to prev.harv.
	INTEGER MCROP				! IN:Crop code number
	INTEGER MEND				! OUT:Number of timesteps from prev.harv.-harv.
	INTEGER NSOW				! OUT:Timesteps from prev.harv. to sowing date
	INTEGER N_STEPS				! OUT:No.timesteps in a year
      INTEGER N_REMAIN
	REAL TCNEED					! IN/OUT:Total N requirement (kgN/ha)
C
C ...Fertiliser applications
C
	INTEGER IFERT(MAXFERT)		! IN:No.timesteps to split fert.application
	INTEGER ILAB(MAXFERT)		! IN:Labelling on fertiliser? 0=No 1=Yes
	REAL FERT(MAXFERT)			! IN:Amount of fertiliser applied (kgN/ha)
	INTEGER NFERT				! IN:No.fertiliser applications to this crop
	REAL TFERT(MAXFERT,3)		! IN:Prop.NO3,NH4,urea in fertiliser
C
C Set time factors from SECONDS
C
      N_STEPS=(365.25*24*60*60)/SECONDS
	CONVER_F=SECONDS/(7*24*60*60)
	N_REMAIN=10/CONVER_F
C
C Assume N content can be set from a constant C:N ratio for each LU type
C
      MCROP=LUCODE
c***** Temporary change 138 to correct concepts of limited data run
	PLANT_REQN=(PIC_ANN*(1+DEBTOPLANTC(LUCODE)))/CTON(LUCODE)	! Temp change 146
c	PLANT_REQN=PIC_ANN/CTON(LUCODE)								! Temp change 148
c***** End of temporary 138 change to correct concepts of limited data run
	CTOT=0
C
C Get start and end of growing season (sowing and harvest dates) from standard values
C
      IF(IYEAR.EQ.1)THEN
	  LHARV=0
      ELSE
	  LHARV=IHARV
	ENDIF
      ISOWN=SOWMONTH(LUCODE)+(IYEAR-1)*12
      NSOW=ISOWN-(IYEAR-1)*12 !LHARV
	IHARV=HARVMONTH(LUCODE)+(IYEAR-1)*12
      MEND=IHARV-(IYEAR-1)*12 !LHARV
	  MEND=min(MEND,12)    					!new line to prevent crash with 8 months eoj
      IEND=MEND-NSOW+1
      IANTHES=ANTHMONTH(LUCODE)+(IYEAR-1)*12
	TCNEED=0
C
C Fertiliser applications - assume plant N requirement is applied as inorganic fertiliser
c (assume no FYM)
C
       nfert=1
	IF(FERTMONTH(LUCODE).GT.0)NFERT=1
      DO 100 NF=1,NFERT
        IFERT(NF)=(FERTMONTH(LUCODE)+(IYEAR-1)*12)-LHARV
	  ILAB(NF)=0
	  FERT(NF)=PLANT_REQN 
	  IF(FERT(NF).GT.TOPFERT)FERT(NF)=TOPFERT
	  TFERT(NF,1)=50	
	  TFERT(NF,2)=50	
	  TFERT(NF,3)=0		
100   CONTINUE
      JSTOP=0
C
C Save internal variables
C
      ISAVE=1
	CALL SAVE_GIS_CROP(ISAVE,N_REMAIN,CONVER_F,MCROP,TCNEED)
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE RUN1_GIS_CROP(IYEAR,IK,LUCODE,IEND,MEND,NSOW,
     &                         N_STEPS,IS_TS,PI_ANN,CREQN,
     &                         CACTOT,CATOT15,ICOVER,SXORGN,ORGN,
     &                         TRNIN,RNIN,RNIN15,C_TS,TCINP,WR)
C
C Calculate crop C and N returns and N offtake 
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXCROP				! Max.no.of crop types allows
      PARAMETER (MAXCROP=36)
      INTEGER I_TS				! Number of weeks since sowing
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=10)
	INTEGER ILU					! Local counter for LU
	INTEGER	INULL				! Unused variable
	INTEGER CHANGEYR			! Year when LU changed
	INTEGER LASTLU				! Last LU type
	REAL FRACIN					! Fraction of equilibrium plant input added
C Plants do not reach equilibrium plant input immediately after LU change
	INTEGER EQUIYEARS(MAXLU)	! Years to reach equilibrium plant input 
c	DATA (EQUIYEARS(ILU),ILU=1,MAXLU) /1,5,10,15/ ! Temp change 118
	DATA (EQUIYEARS(ILU),ILU=1,MAXLU) /1,1,20,1,6,5,1,1,15,1/ ! MLR: Updated for ELUM LUTs
	REAL NONEQPI(MAXLU)			! Fraction change to non-eq PI after LU change
c	DATA (NONEQPI(ILU),ILU=1,MAXLU) /0.0,0.25,1.0,1.0/ ! Temp change 118
	DATA (NONEQPI(ILU),ILU=1,MAXLU) /0.0,0.0,0.4,0.0,0.0,0.5,0.0,0.0,0.4,
     &                                 0.0/ !MLR updated for ELUM LUTs
C Assume straw not incorporated
	INTEGER INSTRAW				! IN:Straw incorporated? 0=No 1=Yes
	DATA INSTRAW /0/
C Fraction N lost by senescence = 9 kg N / ha
C Seed N = 4 kg N / ha
C
C Assume below ground plant N is given by 7t/ha WWheat parameters for arable crops, grass (also forestry and natural)
C Debris N = 63 x (1-exp(-0.63xYld)) = 62 kg N / ha
C Above ground crop N = 236 x (exp(0.071xYld)-1) = 152 kg N / ha
C Total N requirement = 214 kg N / ha
C Prop.of total N req.that goes below ground before harv = Debris N/Total N requirement
C Debris N/(Total N req) = 62/214 = 0.14
      REAL DEBTONREQ(MAXLU)	! Prop.of total N req.that goes below 
							! ground before harvest
c	DATA (DEBTONREQ(ILU),ILU=1,MAXLU) /0.29,1.00,1.00,1.00,0.67,1.00/		! Temp change 146
	DATA (DEBTONREQ(ILU),ILU=1,MAXLU) /1.00,1.00,1.00,1.00,1.0,1.0,
     &                                   1.0,1.0,1.0,1.0/		! Temp change 148 ! MLR
C Assume plant C is given by 7t/ha WWheat parameters for arable crops, grass (also forestry and natural)
C Stubble C = 1446(1-0.94(exp(-0.175xYld))) = 1047 kg/ha
C Straw C = 0.4 x((850 x Yld) x ((1/0.45)-1)) = 2909 kg/ha
C Total C in tops at harvest = Straw C + Stubble C = 2909 + 1047 = 3956 kg / ha
C Debris C (inc.stubble) = 1346(1.23 + 1.4(1-exp(-0.24 x Yld)))= 3189 kg/ha
C Ratio of debris inc.before harv.to above ground C = (Debris C - Stubble C) / Tot.C at harvest = (3189 - 1047)/3956 = 0.54
      REAL DEBTOPLANTC(MAXLU)		! Ratio of C that goes below ground as 
								! debris before harvest to total C in 
								! above group crop at harvest
c	DATA (DEBTOPLANTC(ILU),ILU=1,MAXLU) /0.54,0.54,1.00,1.00/	! Temp change 118
c	DATA (DEBTOPLANTC(ILU),ILU=1,MAXLU) /0.54,0.75,1.00,1.00/	! Temp change 144
c	DATA (DEBTOPLANTC(ILU),ILU=1,MAXLU) /0.54,1.0,1.0,1.0,0.75,1.0/	! Temp change 146
	DATA (DEBTOPLANTC(ILU),ILU=1,MAXLU) /1.00,1.00,1.00,1.00,1.00,1.00,
     &                                     1.00,1.00,1.00,1.0/	! Temp change 148 ! MLR
C Assume 7t/ha WWheat parameters for arable crops
C Assume debris that must be left in field (stubble for arable) is as for cereals
C Above ground crop N = 236 x (exp(0.071xYld)-1) = 152 kg N / ha
C Stubble N =0.12 x Above ground crop N =46 kg N / ha
C Debris N = 63 x (1-exp(-0.63xYld)) = 62 kg N / ha
C Crop N req = Above ground crop N + Debris N = 152+62 = 214 kg N / ha
C Prop.stubble in tot.crop N req.= Stubble N / Crop N req = 46/214 = 0.085
C For other land use types it is assumed to be 0
c***** Temporary change 138 to correct concepts of limited data run ! Temp change 146
	REAL PSTUBBLEN(MAXLU)		! Proportion stubble in total N requirement
c	DATA (PSTUBBLEN(ILU),ILU=1,MAXLU) /0.085,0.0,0.0,0.0,0.0,0.0/	
	DATA (PSTUBBLEN(ILU),ILU=1,MAXLU) /0.085,0.0,0.0,0.0,0.0,0.0,
     &                                   0.0,0.0,0.0,0.0/	!MLR
C Assume 7t/ha WWheat parameters for arable crops
C Stubble C = 1446(1-0.94(exp(-0.175xYld))) = 1047 kg/ha
C Straw C = 0.4 x((850 x Yld) x ((1/0.45)-1)) = 2909 kg/ha
C Total C in tops at harvest = Straw C + Stubble C = 2909 + 1047 = 3956 kg / ha
C Ratio of stubble to above ground plant C = Stubble C / (Stubble C + StrawC) = 1047/(2909+1047) = 0.26
	REAL STUBTOPLANTC(MAXLU)	! Ratio of C that goes below ground as 
								! stubble to total C in 
								! above group crop at harvest
c	DATA (STUBTOPLANTC(ILU),ILU=1,MAXLU) /0.0,0.0,0.0,0.0,0.0,0.0/	
	DATA (STUBTOPLANTC(ILU),ILU=1,MAXLU) /0.0,0.0,0.0,0.0,0.0,0.0,
     &                                      0.0,0.0,0.0,0.0/	!MLR
c***** End of temporary change 138 to correct concepts of limited data run ! Temp change 146
********************************************
* Need data on                             *
* Fraction of N lost by senescence per LU  *
* Proportion non-cartable debris per LU    *
* Proportion stubble per LU                *
* Prop.of total N that is below ground     *
* Proportion straw per LU                  *
* Years for plant input to come to =brium  *
********************************************
C
C Variables passed to/from other subroutines
C
	REAL CACTOT					! IN:N taken up by crop (kgN/ha)
	REAL CATOT15				! IN:N15 taken up by crop (kgN15/ha)
	REAL CONVER_F				! IN:Conversion between this timestep & weeks
	INTEGER IANTHES				! IN:No.timesteps from 01/01 to anthesis 
      INTEGER ICROP				! OUT:Current crop number
	INTEGER IEND				! IN:Timesteps from sowing date to harvest
	INTEGER IK					! IN:No.timesteps from prev.harv. to current
	INTEGER IS_TS				! IN:Crop counter
	INTEGER ISAVE				! Code to save (0) or retrieve (1) internal 
								! gis values
	INTEGER ISOWN				! IN:Timesteps from 01/01/01 to sowing date 
	INTEGER ISTART				! IN:Sowing month
	INTEGER IYEAR				! IN:Current growing season number
	INTEGER LUCODE				! IN/OUT:LU code for current land use
								!(1=arable,2=grassland,3=forestry,4=semi-natural)
	INTEGER MCROP				! IN:Crop code number
	INTEGER MEND				! OUT:Number of timesteps from prev.harv.-harv.
	INTEGER NSOW				! IN:Timesteps from prev.harv. to sowing date
	INTEGER N_STEPS				! OUT:No.timesteps in a year
      INTEGER N_REMAIN

	REAL ANINP					! IN:N added to soil as non-cartable deb. 
								!    in this timestep (kgN/ha/timestep)
C Assume rate of N deb.return follows the same pattern as for cereal crop
	REAL ANRATE(0:MAXCROP)		! OUT:Rate factor for C in non-cartable deb
c	DATA (ANRATE(ILU),ILU=1,MAXLU) /0.1,0.1,0.1,0.1,0.1,0.1/ 
	DATA (ANRATE(ILU),ILU=1,MAXLU) /0.1,0.1,0.1,0.1,0.1,0.1,
     &                                0.1,0.1,0.1,0.1/ !MLR
	REAL CINP					! IN:C added to soil as non-cartable deb.
								!    in this timestep (kgN/ha/timestep)
C Assume rate of C deb.return follows the same pattern as for cereal crop
	REAL CRATE(0:MAXCROP)		! OUT:Rate factor for C in non-cartable deb
	DATA (CRATE(ILU),ILU=1,MAXLU) /0.15,0.15,0.15,0.15,0.15,0.15,
     &                               0.15,0.15,0.15,0.15/  !MLR
	REAL CREQN					! IN:Amount of N required by plant each year
	                            !     (kgN/ha/year)
	REAL CTOT					! N uptake to this week (kgN/ha)
	REAL C_TS					! IN:Deb.C input in this timestep (kgC/ha)
      INTEGER ICOVER				! OUT:Crop cover 1=Covered 0=Bare
	INTEGER JSTOP				! IN/OUT:Code to stop simulation 0=No 1=Yes
	REAL ORGC					! OUT:Total org.C input (kgC/ha) 
	REAL ORGN					! OUT:Total org.N input (kgN/ha) 
	REAL PI_ANN					! IN:Annual plant input (kgC/ha/year)
	REAL RNIN					! IN/OUT:Deb.N input in this timestep (kgN/ha)
	REAL RNIN15					! IN/OUT:Deb.N15 input in timestep (kgN15/ha)
	REAL SORGC					! OUT:Stubble org.C input (kgC/ha) 
	REAL SORGN					! OUT:Stubble org.N input (kgN/ha) 
	REAL SXORGC					! OUT:Prev.stubble org.C input (kgC/ha) 
	REAL SXORGN					! OUT:Prev.stubble org.N input (kgN/ha) 
	REAL SXORGN15				! OUT:Prev.stubble org.N15 input (kgN15/ha) 
	REAL TCINP					! IN/OUT:Total org.C input (kgC/ha) 
	REAL TCNEED					! IN/OUT:Total N requirement (kgN/ha)
	REAL TRNIN					! IN/OUT:Total litter N input (kgN/ha)
	REAL WR						! IN/OUT:Root requirement (kgN/ha)	
	REAL XORGC					! OUT:Prev.total org.C input (kgC/ha) 
	REAL XORGN					! OUT:Prev.total org.N input (kgN/ha) 
C
C Save LASTLU and CHANGEYR
C
      SAVE
C
C Retrieve crop characteristics
C
      ISAVE=0
	CALL SAVE_GIS_CROP(ISAVE,N_REMAIN,CONVER_F,MCROP,TCNEED)
C
C Set total organic C and N input
C
      IF(IK.LE.NSOW)THEN
	  IF(IYEAR.GT.1)THEN 
	    XORGC=ORGC
	    XORGN=ORGN
	    SXORGC=SORGC
	    SXORGN=SORGN
        ENDIF
	  MCROP=LUCODE
	  ICROP=LUCODE
c***** Temporary change 138 to correct concepts of limited data run
	  SORGC=STUBTOPLANTC(LUCODE)*(PI_ANN)	! Temp change 146
        SORGN=PSTUBBLEN(LUCODE)*CREQN		! Temp change 146
c        SORGC=(1-DEBTOPLANTC(LUCODE))*PI_ANN	! Temp change 148
c        SORGN=(1-DEBTONREQ(LUCODE))*CREQN		! Temp change 148
c***** End of temporary change 138 to correct concepts of limited data run
	  IF(SORGC.LT.0.OR.SORGN.LT.0)THEN
	    SORGC=0
	    SORGN=0
        ENDIF
	  ORGC=DEBTOPLANTC(LUCODE)*PI_ANN
	  ORGN=DEBTONREQ(LUCODE)*CREQN
        
C Non-equilibrium plant input after land use change - ASSUMED LINEAR - SHOULD BE EXP
! >>> Commented out by MLR:
!        IF(IYEAR.EQ.1)THEN
!         CHANGEYR=0
!c***** Temporary change 146 to correct concepts of limited data run
!	   SXORGN=SORGN
!	   SXORGC=SORGC
!c***** End of temporary change 146 to correct concepts of limited data run
!	  ELSEIF(IYEAR.GT.1.AND.LUCODE.NE.LASTLU)THEN
!	    CHANGEYR=IYEAR
!	  ENDIF
!        LASTLU=LUCODE
!	  IF(CHANGEYR.GT.0)THEN
!	    IF((IYEAR-CHANGEYR+1).LT.EQUIYEARS(LUCODE))THEN
!		  FRACIN=1.-((IYEAR-CHANGEYR+1.)/EQUIYEARS(LUCODE))
!	      FRACIN=1.-(NONEQPI(LUCODE)*FRACIN)
!	    ELSE
!	      FRACIN=1
!	    ENDIF
!          IF(FRACIN.LT.0)FRACIN=0
!          ORGC=ORGC*FRACIN
!          ORGN=ORGN*FRACIN
!        ENDIF
! >>> End of change by MLR:
        
C Add stubble later
c***** Temporary change 139 to correct concepts of limited data run
	  ORGC=ORGC-SORGC					! Temp change 146
	  ORGN=ORGN-SORGN					! Temp change 146
	  IF(ORGC.LT.0.OR.ORGN.LT.0)THEN	! Temp change 146
	    ORGC=0							! Temp change 146
	    ORGN=0							! Temp change 146
	  ENDIF								! Temp change 146
c***** End of temporary change 139 to correct concepts of limited data run
	  IF(IYEAR.EQ.1)THEN 
   	    XORGC=0
	    XORGN=0
	    SXORGC=0
	    SXORGN=0
        ENDIF
	ENDIF
C
C Count number of timesteps since sowing
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
C Set the end of the growing season for the sown crop
C
      CALL SETEND(IYEAR,IEND,INULL,N_STEPS,ISOWN,MEND,NSOW)
C
C Initialize crop variables
C
      IF(IK.EQ.NSOW)THEN
        CALL SETVAR3(IS_TS,CACTOT,CATOT15,CTOT,TCNEED,JSTOP)
	ENDIF
C
C Cycle root N from total crop N back to RO pool
C  N and C is cycled back from dead roots to RO pool, starting
C in the week after crop uptake starts (i.e. IK.GT.NSOW)
C 
      IF(IYEAR.GT.0.AND.IK.GT.NSOW.AND.IK.LE.MEND)THEN
        CALL RCYCLE(MEND,IK,N_REMAIN,ICROP, 
     &              WR,CACTOT,CATOT15,RNIN,RNIN15,CONVER_F)
	  TRNIN=TRNIN+RNIN
      ENDIF
C
C Work out if there is crop cover or not
C
      IF(IYEAR.GT.0.AND.IK.GT.NSOW.AND.IK.LE.MEND)THEN !+ANINT(4/CONVER_F) removed eoj. So spinup Icover agrees with model run
	  ICOVER=1
      ELSE
        ICOVER=0
      ENDIF
C
C Save crop characteristics
C
      ISAVE=1
	CALL SAVE_GIS_CROP(ISAVE,N_REMAIN,CONVER_F,MCROP,TCNEED)
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE RUN2_GIS_CROP(AIRTEMP,AMMN,AMMN15,ATM,ATM15,
     &                         CACT,CACT15,CACTOT,CATOT15,CLOSSX,
     &                         CLOSSX15,CREQN,CRIT,CTOT,CUPTN,
     &                         IK,IS_TS,IYEAR,JSTOP,MEND,NSOIL,NSOW,
     &                         OCROPN,ORGN,RNIN,SOIL15,SOILN,SRNIN,
     &                         SXORGN,TACTOT,TAM,TC,TIN,TRNIN,VOLAT,
     &                         VOLAT15,XORGN,SEEDIN,PLANTUP,DDAYS)	
C
C Subroutine to take up N and return C and N during season
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXCROP				! Max.no.of crop types allows
      PARAMETER (MAXCROP=36)
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=10)
	INTEGER MAXSOIL				! Max.no.of soil types
	PARAMETER (MAXSOIL=50)
	INTEGER ILU					! Local LU counter 
C Assume arable parameters are equivalent to WWheat
	REAL C1(0:MAXCROP)			! Crop uptake parameter
	DATA (C1(ILU),ILU=1,MAXLU) /0.003,0.003,0.003,0.003,0.003,0.003,
     &                            0.003,0.003,0.003,0.003/
      INTEGER IROCKS(0:MAXCROP)	! Maximum rooting depth
	DATA (IROCKS(ILU),ILU=1,MAXLU) /150,50,300,50,50,300,200,150,300,150/  ! MLR changed for to 50 for temp fix
      INTEGER L_TSS(0:MAXCROP)	! Weeks from anthesis to end of growing season
      DATA (L_TSS(ILU),ILU=1,MAXLU) /5,20,15,15,20,15,5,20,15,5/
	INTEGER IL_TSS
	REAL RRG(0:MAXCROP)			! Rate of root growth (cm/week)
	DATA (RRG(ILU),ILU=1,MAXLU) /5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0/
      REAL RRX(0:MAXCROP)			! Root exploration factor
	DATA (RRX(ILU),ILU=1,MAXLU) /25.0, 25.0, 25.0, 25.0, 25.0, 25.0,
     &                            25.0, 25.0, 25.0, 25.0/
	REAL SEED(0:MAXCROP)		! Seed N (arable only- others LUs perenniel)
	DATA (SEED(ILU),ILU=1,MAXLU) /4,0,0,0,0,0,4,4,0,4/
      REAL T1(0:MAXCROP)			! Crop uptake parameter
	DATA (T1(ILU),ILU=1,MAXLU) /1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5/
C
C Variables passed to/from calling subroutine
C
	REAL AIRTEMP				! IN:Air temperature (deg.C)
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
      REAL ATM					! IN:Atmospheric N input (kgN/ha/timestep)
	REAL ATM15					! IN:Atmospheric N15 input (kgN15/ha/time)
	REAL CACT					! IN:Actual N uptake (kgN/ha)
	REAL CACT15					! IN:Actual N15 uptake (kgN15/ha)
	REAL CACTOT					! IN:N taken up by crop (kgN/ha)
	REAL CATOT15				! IN:N15 taken up by crop (kgN15/ha)
	REAL CLOSSX					! IN:N lost by senescence (kgN/ha) 
	REAL CLOSSX15				! IN:N15 lost by senescence (kgN15/ha) 	
	REAL CONVER_F				! IN:Conversion between this timestep & weeks
	REAL CREQN					! IN/OUT:Amount of N required by plant each yr
	                            !        (kgN/ha/year)
	REAL CRIT(MAXSOIL,MAXLAYER)	! IN:Critical N content (kgN/ha/layer)
	REAL CTOT					! IN:Total litter C input (kgC/ha)
      REAL CUPTN					! IN:Crop N requirement (kgN/ha)
      REAL DDAYS					! OUT:Degree days since sowing (deg.C)
	INTEGER IK					! IN:No.timesteps from prev.harv. to current
	INTEGER ISAVE				! OUT:Code to save (0) or retrieve (1) internal 
								! gis values
	INTEGER IS_TS				! IN:Crop counter
	INTEGER IYEAR				! IN:Current growing season number
	INTEGER JSTOP				! IN/OUT:Code to stop simulation 0=No 1=Yes
	INTEGER MCROP				! IN:Crop code number
	INTEGER MEND				! IN:Number of timesteps from prev.harv.-harv.
	INTEGER NSOIL				! IN:Soil code number
	INTEGER NSOW				! IN:Timesteps from prev.harv. to sowing date
	INTEGER N_STEPS				! OUT:No.timesteps in a year
      INTEGER N_REMAIN			! IN/OUT
	REAL OCROPN					! IN:N uptake of crop (kgN/ha) 
								! (0->calculate within model)
	REAL ORGN					! IN:Total org.N input (kgN/ha) 
	REAL PLANTUP(MAXLAYER)		! IN(CROPGROW_GIS)/OUT:Plant uptake from each soil layer
	REAL RNIN					! IN:Litter N input in this timestep (kgN/ha)
      REAL ROOTS					! IN:root depth at 5cm/week since sowing
	REAL RORGN					! OUT:Total org.N input (kgN/ha) 
	REAL SEEDIN					! IN: Data used by MAGEC
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SRNIN					! IN:Total litter N input (kgN/ha)
	REAL SXORGN					! IN:Stubble org.N15 input (kgN/ha) DONT PASS?
	REAL TACTOT					! IN:Total crop N uptake (kgN/ha)
	REAL TAM(MAXLAYER)			! IN:15N/N in the ammonium in the layer
	REAL TC						! IN/OUT:Total org.N15 input (kgN/ha) DONT PASS?
	REAL TCNEED					! IN/OUT:Total N requirement (kgN/ha)
	REAL TIN(MAXLAYER)			! IN:15N/N in the nitrate in the layer
	REAL TRNIN					! IN:Total litter N input (kgN/ha)
	REAL VOLAT					! IN:N lost by volatilisation (kgN/ha)
	REAL VOLAT15				! IN:N15 lost by volatilisation (kgN15/ha)
	REAL XORGN					! IN:Total org.N input (kgN/ha) 
C
C Retrieve crop characteristics
C
      ISAVE=0
	CALL SAVE_GIS_CROP(ISAVE,N_REMAIN,CONVER_F,MCROP,TCNEED)
C
C Calulate the number of degree days since sowing
C
	CALL DAYD(IYEAR,IK,NSOW,DDAYS,AIRTEMP,CONVER_F) 
C
C If freezing or not in growing season set crop uptake parameters to zero
C
      IF((AIRTEMP.LE.0.0).OR.
     &   (IYEAR.GT.0.AND.IK.LE.NSOW).OR.(JSTOP.EQ.1))THEN
        CALL ZEROUP(CUPTN,CACT,CACT15)
	ENDIF
C
C Otherwise calculate the N/C turnover associated with growing the crop
C
      IF(JSTOP.EQ.1)THEN
       CALL ZEROUP(CUPTN,CACT,CACT15)
	ENDIF
 	JSTOP=0 ! Added to avoid non-zero difference in N balance. 
			! Error due to notentering the condition in EXTRACT...
			!          IF(UPTN.GE.(SCRIT+ACRIT))THEN
			! Solve properly when time allows
      IF((AIRTEMP.GT.0.0).AND.
     &   (IYEAR.GT.0.AND.IK.GE.NSOW).AND.(JSTOP.NE.1))THEN
        CALL CROPGROW_GIS(MCROP,IYEAR,IK,NSOW,MEND,IS_TS,NSOIL, 
     &                CACT15,SOILN,TIN,SOIL15,AMMN15,AMMN,TAM,
     &                DDAYS,CUPTN,ROOTS,CONVER_F,ATM,ATM15,
     &                T1,C1,RRG,IROCKS,RRX,CRIT,CREQN,CACT,CTOT,TCNEED,
     &                PLANTUP)

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
      IF(JSTOP.EQ.0)CLOSSX=0
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
        RORGN=SXORGN+ORGN
        CALL SENES(IYEAR,MCROP,JSTOP,OCROPN,SRNIN,RORGN,TACTOT,
     &             CONVER_F,CLOSSX,L_TSS,XORGN)
      ENDIF
C
C If stopping take off losses due to senescence
C
      IF(JSTOP.EQ.1)THEN
        IF(CLOSSX.LT.0)CLOSSX=0
	  IF(CLOSSX.LT.CACTOT)THEN
         CLOSSX15=CLOSSX*CATOT15/CACTOT
          CACTOT=CACTOT-CLOSSX
          CATOT15=CATOT15-CLOSSX15
          VOLAT=VOLAT+CLOSSX
          VOLAT15=VOLAT15+CLOSSX15
	  ELSE
	    IF(CACTOT.GT.0)THEN
            CLOSSX15=CLOSSX*CATOT15/CACTOT
	    ELSE
	      CLOSSX15=0
	    ENDIF
	    CLOSSX=CACTOT
          CACTOT=0
          CATOT15=CATOT15-CLOSSX15
          VOLAT=VOLAT+CLOSSX
          VOLAT15=VOLAT15+CLOSSX15
	  ENDIF
      ENDIF
C
C Save crop characteristics
C
      ISAVE=1
	CALL SAVE_GIS_CROP(ISAVE,N_REMAIN,CONVER_F,MCROP,TCNEED)
      CALL GETCAC(CACTOT,CACT,CATOT15,CACT15)
C
C Leave RUN2_SUNDIAL_CROP
C
      END

C
C-------------------------------------------------------------
C
      SUBROUTINE SET_NLIM(THISLU,ICFACTOR,
     &                           SOILN,SOIL15,AMMN,AMMN15,SOILW,
     &                           DPMCARB0,DPMNIT0,DPMNLAB0,
     &                           RPMCARB0,RPMNIT0,RPMNLAB0,
     &                           BCARB0,BNIT0,BNLAB0,
     &                           HCARB0,HNIT0,HNLAB0,IOM,PI_CEQ_MON)		  
C
C Set soil state and parameters after N limited equilibrium run
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXORGM				! Max.no.of organic manure applications allowed
	PARAMETER (MAXORGM=52)
	INTEGER MAXSOIL				! Max.no.of soil types
	PARAMETER (MAXSOIL=50)
C
C Variables passed to/from other subroutines
C
	REAL ALPHA(MAXLAYER)			! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL AMMN(MAXLAYER)				! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)			! IN:Soil ammonium-N15 (kgN15/ha/layer)
	REAL ANIT						! Ammonium immobilised (kgN/ha)
	REAL ANIT15						! Ammonium N15 nitrified (kgN/ha)
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
	REAL BYRAIN						! IN/OUT:Rain lost by bypass flow mm/timestep	
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: % clay content in this layer
      REAL CONVER_F					! IN/OUT:Conversion to correct timestep
	REAL CRIT(MAXSOIL,MAXLAYER)		! IN:Critical N content (kgN/ha/layer)
	REAL DELTA(MAXLAYER)			! IN/OUT: Prop.HUM produced on HUM decompn 
	REAL DFACT(MAXLAYER)			! IN/OUT: Denitrification factor
	REAL DPMCARB0(MAXLAYER)			! IN:C in decomposable PM (kgC/ha/layer)
	REAL DPMNIT0(MAXLAYER)			! IN:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)			! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL DPMRATE(MAXLAYER)			! IN/OUT: Rate constant for DPM decompn
	REAL FANIT						! Fertiliser N nitrified (kgN/ha)
	REAL FANIT15					! Fertiliser N15 nitrified (kgN/ha)
	REAL FLOWPROP				! IN:Proportion of flow needed to achieve
							    !	     water table at depth WTABLE
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
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
	REAL IOMARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium IOM (kgC/ha/layer)
      INTEGER ISAVE					! OUT:Code to save or retrieve variables
									!		0=save; 1=retrieve
	REAL LTA_AWC(12,MAXLAYER)	! IN:Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)	! IN:Average air temp this month (deg.C)
      INTEGER LUARRAY(MAXSOIL)	! IN/OUT:Land use before equilibrium 
      INTEGER THISLU					! IN:Land use code
      INTEGER NUMSOIL					! Number of soils defined
	CHARACTER*20 PETARRAY(MAXSOIL)	! IN/OUT:File containing monthly PET for 
									! the period between equilibrium and start of simulation. 

	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)! IN/OUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)! IN/OUT: pH above which decomp.max.
	REAL PI_CEQ_MON(12,MAXLAYER)! IN/OUT: Equilibrium plant C input each month
								!         in each layer (kgC/ha/month/layer)
	CHARACTER*20 RAINARRAY(MAXSOIL)	! IN/OUT:File containing monthly rain for 
									! the period between equilibrium and start of simulation. 
	REAL RPMCARB0(MAXLAYER)			! IN:C in resistant PM (kgC/ha/layer)
	REAL RPMNIT0(MAXLAYER)			! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)			! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL RPMRATE(MAXLAYER)			! IN/OUT: Rate constant for RPM decompn
	REAL SECONDS					! IN:Number of seconds in one timestep
	CHARACTER*40 SNAME(MAXSOIL)		! IN/OUT:Soil name
	REAL SOIL15(MAXLAYER)			! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILN(MAXLAYER)			! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOILW(MAXLAYER)			! IN:Available water (mm/layer)
	CHARACTER*20 TEMPARRAY(MAXSOIL) ! IN/OUT:File containing monthly temp.for 
									! the period between equilibrium and start of simulation. 
	INTEGER TIMEARRAY(MAXSOIL)		! IN/OUT:Time between soil at equilibrium 
									! and start of simulation
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium TOC (kgC/ha/layer)
	REAL WMAX(MAXLAYER)			! IN:Available water at field capacity(mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
C
C Retrieve current state of the soil and save for C and N model
C
      ISAVE=0 ! Retrieve soil 
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
      CALL SAVE_GIS_SOIL(ISAVE,THISLU,
     &            SOILN,SOIL15,AMMN,AMMN15,
     &            SOILW,WSAT,WMAX,FLOWPROP,LTA_AWC,LTA_TEMP,
     &            DPMCARB0,DPMNIT0,DPMNLAB0,
     &            RPMCARB0,RPMNIT0,RPMNLAB0,
     &            BCARB0,BNIT0,BNLAB0,
     &            HCARB0,HNIT0,HNLAB0,
     &            IOM,PI_CEQ_MON)
      CALL SAVE_GIS_RATES(ISAVE,THISLU,ICFACTOR)
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
C Leave SET_NLIM
C
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE SET_NLIM_RATE(THISLU,ICFACTOR)		  
C
C Set soil state and parameters after N limited equilibrium run
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXORGM				! Max.no.of organic manure applications allowed
	PARAMETER (MAXORGM=52)
	INTEGER MAXSOIL				! Max.no.of soil types
	PARAMETER (MAXSOIL=50)
C
C Variables passed to/from other subroutines
C
	REAL ALPHA(MAXLAYER)			! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL AMMN(MAXLAYER)				! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)			! IN:Soil ammonium-N15 (kgN15/ha/layer)
	REAL ANIT						! Ammonium immobilised (kgN/ha)
	REAL ANIT15						! Ammonium N15 nitrified (kgN/ha)
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
	REAL CRIT(MAXSOIL,MAXLAYER)		! IN:Critical N content (kgN/ha/layer)
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
      INTEGER ISAVE					! OUT:Code to save or retrieve variables
									!		0=save; 1=retrieve
      INTEGER LUARRAY(MAXSOIL)	! IN/OUT:Land use before equilibrium 
      INTEGER THISLU					! IN:Land use code
      INTEGER NUMSOIL					! Number of soils defined

	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)! IN/OUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)! IN/OUT: pH above which decomp.max.
	REAL RPMRATE(MAXLAYER)			! IN/OUT: Rate constant for RPM decompn
	REAL SECONDS					! IN:Number of seconds in one timestep
	CHARACTER*40 SNAME(MAXSOIL)		! IN/OUT:Soil name
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium TOC (kgC/ha/layer)

C
C Retrieve current state of the soil and save for C and N model
C
      ISAVE=0 ! Retrieve soil 
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
      CALL SAVE_GIS_RATES(ISAVE,THISLU,ICFACTOR)
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
C Leave SET_NLIM_RATE
C
      END
C
C-----------------------------------------------------------------------
C
C INTERNAL SUBROUTINES
C 1. CROPGROW_GIS
C 2. GET_PET
C 3. GET_LD
C 4. GET_TOTC
C 5. SAVE_GIS_CROP
C 6. SAVE_GIS_SOIL
C
C-----------------------------------------------------------
C
      SUBROUTINE CROPGROW_GIS(MCROP,IYEAR,IK,NSOW,MEND,IS_TS,NSOIL, 
     &                    CACT15,SOILN,TIN,SOIL15,AMMN15,AMMN,TAM,
     &                    DDAYS,CUPTN,ROOTS,CONVER_F,ATM,ATM15,
     &                    T1,C1,RRG,IROCKS,RRX,CRIT,CREQN,CACT,CTOT,
     &                    TCNEED,PLANTUP)
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
      INTEGER M
      REAL FIXFRAC(0:MAXCROP),CNEED,CFIXNEED
	DATA (FIXFRAC(M),M=0,6) /0,0,1,1,1,1,1/  !MLR changed so that forest has fixfrac
C
C Variables passed to/from calling subroutine
C
	REAL AMMN(MAXLAYER)
	REAL AMMN15(MAXLAYER)
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
	REAL CTOT
	REAL CUPTN
	REAL CURRENTLAI			! Current LAI
	REAL CURRENTPLBIOM		! Current plant biomass [ kg ha-1]
	REAL CURRENTPLHEI		! Current plant height [cm]
	REAL DDAYS
	REAL FIXN
	INTEGER IK
      INTEGER IROCKS(0:MAXCROP)
	INTEGER IS_TS
	INTEGER IYEAR
	REAL LAI				! Plant LAI at harvest
	INTEGER MCROP
	INTEGER MEND
	INTEGER NSOW
	INTEGER NSOIL
	REAL PLANTBIOMASS		! Plant biomass at harvest [ kg ha-1]
	REAL PLANTHEIGHT		! Plant height at harvest [cm]
	REAL PLANTUP(MAXLAYER)	! IN(EXTRACT)/OUT:Plant uptake from each soil layer
	INTEGER ROOTLAYER		! Layer until which roots have grown
	REAL ROOTS
	REAL RRG(0:MAXCROP)
      REAL RRX(0:MAXCROP)
	REAL SOIL15(MAXLAYER)
      REAL SOILN(MAXLAYER)
      REAL T1(0:MAXCROP)
	REAL TIN(MAXLAYER)
	REAL TAM(MAXLAYER)
	REAL TCNEED
C
C In first year or in growing season set the crop N requirement for this week
C
      IF(IYEAR.EQ.1.OR.(IYEAR.GT.1.AND.IK.GE.NSOW.AND.
     &   IK.LE.MEND))CALL ADDUP(BREQN,CACT,CREQN,TCNEED)
C
C In first year or first week of the growing season
C set the inital rate of crop N uptake
C
      IF(IYEAR.EQ.1.AND.IK.EQ.1)CALL INITUP(CTOT,CREQN,BREQN,CACTOT)
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
C      CUPTN=(1-FIXFRAC(MCROP))*CUPTN !230909
      CALL EXTRACT(M,CACT,NSOIL,CUPTN,SOILN,CRIT,AMMN,PLANTUP)
      CALL EXTR15(CACT15,SOILN,TIN,SOIL15,AMMN15,AMMN,TAM)
      CFIXNEED=CNEED-CACT
	IF(CFIXNEED.GT.(FIXFRAC(MCROP)*CUPTN))THEN !230909
	  CFIXNEED=FIXFRAC(MCROP)*CUPTN !230909
	ENDIF
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
C-------------------------------------------------------------
C
	SUBROUTINE GET_PET(LAT,AVETEMP,AVEPET)
C
C Subroutine to calculate PET from oter weather data 
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      INTEGER IMON				! Month counter
      REAL*8 APAR					! APAR=(1.6/100)I + 0.5 
	REAL*8 HINDEX					! Annual heat index = HINDEX
C
C Variables passed to/from other subroutines
C
      REAL AVEPET(12)				! Long term average PET (mm/month)
      REAL AVETEMP(12)			! Long term average temp. (deg.C/month)
	REAL LAT					! Latitude
      REAL*8 RLD					! Daytime hours in units of 12 
								! (eg for 10 daylight hours, Ld=10/12 (daytime hours = DAYHRS)
C
C Calculate PET by Thornthwaite method
C PET = 1.6 x RLD x pow((10T/I),APAR)*10 (multiply by 10: cm/month -> mm/month)
C PET = potential evapotranspiration in cm/month
C T = average temperature for month in deg.C
C RLD = daytime hours in units of 12 (eg for 10 daylight hours, Ld=10/12 (daytime hours = DAYHRS)
C I = annual heat index = HINDEX
C I= sum(i(month)) 
C i(month)=power((T(month)/5),1.514)
C APAR=(1.6/100)I + 0.5 
C daylight hours calculated according to latitude
C
C Calculate annual heat index, HINDEX
C
	HINDEX=0
	DO 100 IMON=1,12
	  IF(AVETEMP(IMON)/5.GE.0)
     &    HINDEX=HINDEX+((AVETEMP(IMON)/5)**1.514)
100   CONTINUE
C
C Calculate parameter a, APAR
C
      APAR=((1.6/100)*HINDEX)+0.5
C
C Calculate daytime hours, RLD
C
	DO 200 IMON=1,12
 	  CALL GET_LD(LAT,IMON,RLD)
	  AVEPET(IMON)=1.6*RLD*(((10*AVETEMP(IMON))/HINDEX)**APAR)*10
        AVEPET(IMON)=(ANINT(10*AVEPET(IMON)))/10
	  IF(AVETEMP(IMON).LE.0.OR.HINDEX.LT.0)AVEPET(IMON)=0
200   CONTINUE
C
C Leave GET_PET
C
      END
C
C-------------------------------------------------------------
C
	SUBROUTINE GET_THIS_PET(LAT,AIRTEMP,EVAPW,THISMON)
C
C Subroutine to calculate PET from other weather data for this month
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      REAL*8 APAR					! APAR=(1.6/100)I + 0.5 
	REAL*8 HINDEX					! Annual heat index = HINDEX
C
C Variables passed to/from other subroutines
C
	REAL AIRTEMP				! IN: Air temperature (deg.C/timestep)	
	REAL EVAPW					! OUT: Potential evap. (mm/timestep)
	REAL LAT					! IN:Latitude
      REAL*8 RLD					! IN:Daytime hours in units of 12 
								! (eg for 10 daylight hours, Ld=10/12 (daytime hours = DAYHRS)
	INTEGER THISMON				! IN:Current month
C
C Calculate PET by Thornthwaite method
C PET = 1.6 x RLD x pow((10T/I),APAR)*10 (multiply by 10: cm/month -> mm/month)
C PET = potential evapotranspiration in cm/month
C T = average temperature for month in deg.C
C RLD = daytime hours in units of 12 (eg for 10 daylight hours, Ld=10/12 (daytime hours = DAYHRS)
C I = annual heat index = HINDEX
C I= sum(i(month)) 
C i(month)=power((T(month)/5),1.514)
C APAR=(1.6/100)I + 0.5 
C daylight hours calculated according to latitude
C
C Calculate annual heat index, HINDEX
C
	HINDEX=0
	IF(AIRTEMP.GE.0)HINDEX=HINDEX+((AIRTEMP/5)**1.514)
C
C Calculate parameter a, APAR
C
      APAR=((1.6/100)*HINDEX)+0.5
C
C Calculate daytime hours, RLD
C
 	CALL GET_LD(LAT,THISMON,RLD)
	EVAPW=1.6*RLD*(((10*AIRTEMP)/HINDEX)**APAR)*10
	IF(AIRTEMP.LE.0.OR.HINDEX.LT.0)EVAPW=0
C
C Leave GET_THIS_PET
C
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE GET_LD(LAT,IMON,RLD)
C Originally extracted from E. H. Wiser's 
C subprogram in DRAINMOD's climate portion. 
C modified etc, J. E. Parsons 
C INPUTS
C LAT = Latitude of location degrees and minutes 
C form -- ddmm, for example, Raleigh NC we would enter 3552
C OUTPUTS
C Thornewaite's daylight hours (RLD) by day
C RLD = Ld is the fraction of daylight hours 
C where Ld=1 corresponds to 12 hours.
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER IDAYS(12)
	INTEGER M
	DATA (IDAYS(M),M=1,12)/31,28,31,30,31,30,31,31,30,31,30,31/
      REAL*8 DAYS(2)
	INTEGER I,J
	REAL*8 XLAT1
	REAL*8 XLAT2
	REAL*8 RLAT
	REAL*8 SINLAT
	REAL*8 COSLAT
      REAL*8 XM
      REAL*8 XLAM
	REAL*8 YD
      REAL*8 XD
C
C Variables passed to/from this subroutine
C
	REAL LAT					! Latitude
      REAL*8 RLD					! Daytime hours in units of 12 
								! (eg for 10 daylight hours, Ld=10/12 (daytime hours = DAYHRS)
	REAL*8 D
      INTEGER IMON				! Month counter
C
C Extract degrees and minutes from entered latitude and calculate the sin and cos of latitude
C
      I=LAT
      J=(LAT-I)*60
      XLAT1=DFLOAT(I)
      XLAT2=DFLOAT(J)
      RLAT=(0.0174533*XLAT1)+(0.0002909*XLAT2)
      SINLAT=DSIN(RLAT)
      COSLAT=DCOS(RLAT)
C
C Calculate day of start (DAYS(1)) and end (DAYS(2)) of month
C
      DAYS(1)=1
      DO 100 I=1,IMON-1
        DAYS(1)=DAYS(1)+IDAYS(I)
100   CONTINUE
      DAYS(2)=DAYS(1)+IDAYS(IMON)-1
C
C Calculate RLD at start and end of month and take average
C
      RLD=0
      DO 200 I=1,2
        XM=0.0172264*(-0.6+DAYS(I))
        XLAM=4.874239+XM+0.0334762*DSIN(XM)+0.0003502*DSIN(XM+XM)
        YD=0.397900*DSIN(XLAM)
	  IF(1.-YD*YD.LT.0)THEN
	    RLD=1
          RETURN
	  ENDIF
        XD=DSQRT(1.-YD*YD)
        D=DATAN2(YD,XD)
        XD=(-0.0414544-(SINLAT*DSIN(D)))/(COSLAT*DCOS(D))
	  IF(1-XD*XD.LT.0)THEN
	    RLD=1
	    RETURN
	  ENDIF
        YD=DSQRT(1-XD*XD)
        RLD=RLD+0.0111111*DATAN2(YD,XD)*57.29578
200   CONTINUE
      RLD=RLD/2
C
	END

C
C-------------------------------------------------------------
C
      SUBROUTINE GET_TOTC(DPMCARB0,RPMCARB0,BCARB0,HCARB0,TOTC)
C
C Get total carbon
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/

      INTEGER IL					! Local layer counter
C
C Variables passed to/from other subroutines
C
	REAL BCARB0(MAXLAYER)		! IN:C in soil biomass (kgC/ha/layer)
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
	REAL TOTC					! OUT: Total organic C (kgC/ha)
C
C Initialise total C to zero before summing
C
	TOTC=0
C
C Add up C pools in different layers
C
      DO 100 IL=1,MAXLAYER1
	  TOTC=TOTC+BCARB0(IL)+HCARB0(IL)+RPMCARB0(IL)+DPMCARB0(IL)
100   CONTINUE
C
C Leave GET_TOTC
C
      END
C
C------------------------------------------------------------------
C
      SUBROUTINE SAVE_GIS_CROP(ISAVE,N_REMAIN,CONVER_F,MCROP,TCNEED)
C
C Subroutine to save and retrieve crop parameters 
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER N_REMAIN1
	REAL CONVER_F1				! Conversion between this timestep & weeks
	INTEGER MCROP1				! Crop code number
	REAL TCNEED1				! Total N requirement (kgN/ha)
C
C Variables passed to/from calling subroutine
C
      INTEGER ISAVE				! IN/OUT:Code to save (1) or retrieve (0) vars
	INTEGER N_REMAIN
	REAL CONVER_F				! IN/OUT:Conversion between this timestep &wks
	INTEGER MCROP				! IN/OUT:Crop code number
	REAL TCNEED					! IN/OUT:Total N requirement (kgN/ha)
C
C Save variable
C
	SAVE
C
C Save parameters
C
	IF(ISAVE.EQ.1)THEN
	  N_REMAIN1=N_REMAIN
	  CONVER_F1=CONVER_F
	  MCROP1=MCROP
	  TCNEED1=TCNEED
C
C Retrieve parameters
C
      ELSEIF(ISAVE.EQ.0)THEN
	  N_REMAIN=N_REMAIN1
	  CONVER_F=CONVER_F1
	  MCROP=MCROP1
	  TCNEED=TCNEED1
	ENDIF
C
C Leave SAVE_GIS
C
	END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SAVE_GIS_RATES(ISAVE,LU1,ICFACTOR)
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

  	REAL ICFACTOR1(MAXLAYER,MAXLU)	! IN:Adjustment in rate needed to	
									!    achieve measured NPP and TOC	
	INTEGER IL						! Local layer counter variable
C
C Variables passed to / from this routine
C 
  	REAL ICFACTOR(MAXLAYER)		! IN:Adjustment in rate needed to	
								!    achieve measured NPP and TOC	
      INTEGER ISAVE				! IN:Code to save or retrieve variables
      INTEGER LU1					! IN: Land use code
C
C Save soil 
C
      IF(ISAVE.EQ.1)THEN
	  DO 100 IL=1,MAXLAYER1
		ICFACTOR1(IL,LU1)=ICFACTOR(IL)
100     CONTINUE
C
C Retrieve soil
C
	ELSEIF(ISAVE.EQ.0)THEN
	  DO 200 IL=1,MAXLAYER1
		ICFACTOR(IL)=ICFACTOR1(IL,LU1)
200     CONTINUE
	ENDIF
C
C Leave SAVE_GIS_RATE
C
      END

C
C-----------------------------------------------------------------------
C
      SUBROUTINE SAVE_GIS_SOIL(ISAVE,LU1,
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
C Leave SAVE_GIS_SOIL
C
      END

C*************************************************************
C RESULTS ROUTINES
C*************************************************************
C
C EXTERNAL SUBROUTINES
C 1. CHECK_GIS_CELL
C 2. CHECK_GIS_GRID
C 3. CHECK_GIS_LU
C 4. CHECK_GIS_SOIL
C 5. CLOSECHAN2
C 6. GET_GIS_CHANGE
C 7. GET_GIS_EMISSIONS
C 8. GETWEATHER_GIS
C 9. INIT_GIS_RES
C 10. OPENCHAN_GIS
C 11. PUT_GIS_RES 
C 12. SETERROR
C 13. SETFILE_GIS
C 14. SET_GIS_FERT
C 15. SETNOERROR
C 16. TEST1_OPENCHAN
C 16. TEST2_OPENCHAN
C 17. TEST1_RES
C 18. TEST2_RES
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CHECK_GIS_CELL(GISOK,NCELL,AREA,DOMSOIL,PERSOIL,
     &                          WETCLASS,CELLNPP,CELLLAT,RAIN,MAIRTEMP,
     &                          LU1TOLU2,PI_SOURCE)
C
C Subroutine to check 1km2 cell data
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
      INTEGER MAXDEC				! Max.no.of decades included in calculation
	PARAMETER(MAXDEC=9)		
	INTEGER MAXDEC1				! Max.no.of decades included in calculation
	DATA MAXDEC1 /9/
	INTEGER MAXLU				! Max.no.of land use types
      PARAMETER (MAXLU=6)
	INTEGER MAXLU1				! Max.no.of land use types
	DATA MAXLU1 /6/
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXSERIES1
	DATA MAXSERIES1 /5/
      INTEGER MAXWC				! Max. no. wetness classes
	PARAMETER (MAXWC=6)

      INTEGER IDEC				! Local decade counter
	INTEGER ILU1				! Local land use counter
	INTEGER ILU2				! Local land use counter
      INTEGER IMON				! Local month counter
	INTEGER IS					! Local counter for soil series

C
C Variables passed to / from this routine
C 
      REAL MAIRTEMP(12)			! Long term ave. airtemp for 1km2 grid (deg.C)
	REAL AREA					! IN:Area of cell accounted for (m2)
	REAL CELLLAT				! Latitude of 1km2 grid cell
	REAL CELLNPP				! IN:NPP kgC/m2 x 1000 for this 1km2 grid cell
      INTEGER DOMSOIL(CELLSTOGRID,MAXSERIES)	! IN:5 major soil series in the  		
								! 1km2 cell x no. of cells in 20km2 grid
	INTEGER GISOK				! OUT:Check of GIS data
	INTEGER NCELL				! IN:Counter for cells within the 20km2 grid
	REAL LU1TOLU2(CELLSTOGRID,MAXDEC,MAXLU+2,MAXLU+2)	! Fracn of LU1 changed 
														! to LU2 in decade
	REAL PERSOIL(CELLSTOGRID,MAXSERIES)	! IN:% of each major soil type in the
								! 1km2 cell x no. of cells in 20km2 grid
	INTEGER PI_SOURCE			! Source of PI 1=input 2=calc from TOC
	INTEGER PI_INPUT			! Input PI and use for initialisation
	INTEGER PI_FROMTOC			! Calculate PI from TOC during initialisation
	DATA PI_INPUT,PI_FROMTOC /1,2/
      REAL RAIN(12)				! Long term ave. rain for 1km2 cell (mm/month)
      INTEGER WETCLASS(CELLSTOGRID,MAXSERIES) !OUT:Wetness class for 
								! each major soil type in the 1km2 square cell 
								! x 20km2 cell
C
C Check for non-valid entries
C
      GISOK=0
C
C Check area
C
      IF(AREA.LT.0.OR.AREA.GT.1000000)RETURN
C
C Check for one valid soil 
C
      DO 100 IS=1,MAXSERIES1
        IF((DOMSOIL(NCELL,IS).GT.0.AND.PERSOIL(NCELL,IS).GT.0).AND.
     &     (WETCLASS(NCELL,IS).GT.0.AND.WETCLASS(NCELL,IS).LE.MAXWC))
     &     GOTO 101 
100   CONTINUE
      RETURN
101   CONTINUE
C
C Check NPP and latitude
C
      IF(PI_SOURCE.EQ.PI_INPUT)THEN
        IF(CELLNPP.LT.0.OR.CELLNPP.GT.10000)RETURN
	ENDIF
	IF(CELLLAT.LT.0.OR.CELLLAT.GT.100)RETURN
C
C Check rain and air temperature
C
      DO 200 IMON=1,12
	  IF(RAIN(IMON).LT.0)RETURN
200   CONTINUE
      DO 300 IMON=1,12
	  IF(MAIRTEMP(IMON).LT.-5000)RETURN
300   CONTINUE
C
C Look for any LU change
C
	DO 400 IDEC=1,MAXDEC1
	  DO 500 ILU1=1,MAXLU1
	    DO 600 ILU2=1,MAXLU1
     	      IF(LU1TOLU2(NCELL,IDEC,ILU1,ILU2).LT.0)RETURN
600       CONTINUE
500     CONTINUE
400   CONTINUE
C
C All cell data is OK
C	 
     	GISOK=1
	END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CHECK_GIS_GRID(GISOK,TRAPERR,STYPES,TOTSEQ,
     &                          GRIDNPP,AVERAIN,AVEPET,AVETEMP)
C
C Subroutine to check soil characteristics
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXSEQ				! Max.no of sequences included in GIS
	PARAMETER (MAXSEQ=MAXLU*MAXLU)						
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						

	INTEGER IMON						! Local month counter
	INTEGER ISERIES						! Counter for soil series
	INTEGER NSEQ						! Counter for number of land use seq 
C
C Variables passed to / from this routine
C 
      REAL AVEPET(12)						! IN:Long term ave.PET (mm/month)
      REAL AVERAIN(12)					! IN:Long term ave.rainfall(mm/month)
      REAL AVETEMP(12)					! IN:Long term ave.monthly average 
	REAL GRIDNPP						! IN:Ave.NPP in the 1km cell (kgC/ha)
	INTEGER GISOK						! IN:Check of GIS data
	INTEGER STYPES						! IN:Number of soil types in 20km2 grid
	INTEGER TOTSEQ						! IN:Total number of sequences
      INTEGER TRAPERR(CELLSTOGRID*MAXSERIES/2,MAXSEQ) 
								! IN:Error trapping
								!		0=No error
								!		1=Error in NPP data
								!		2=Error in weather data
								!		3=Error in LU data
								!		4=Error in LU change data
								!		5=Soil error
								!		6=Equilibrium not found
	INTEGER NPPERR
	DATA NPPERR /1/						!		1=Error in NPP data
	INTEGER METERR
	DATA METERR /2/						!		2=Error in weather data
C
C Check for non-valid entries
C
      GISOK=0
C
C Check NPP
C
	IF(GRIDNPP.LT.0.OR.GRIDNPP.GT.10000)THEN
	  PRINT*,'Error in 20km2 grid cell NPP: ',GRIDNPP,' kgC/ha/yr'
	  CALL TRAPERROR(TRAPERR,NPPERR,1,STYPES,1,TOTSEQ)
	  RETURN
	ENDIF
C
C Check weather
C
      DO 200 IMON=1,12
        IF(AVERAIN(IMON).LT.0.OR.AVERAIN(IMON).GT.700)THEN
	    PRINT*,'Error in rainfall data: ',AVERAIN(IMON),'mm'
	    CALL TRAPERROR(TRAPERR,METERR,1,STYPES,1,TOTSEQ)
	    RETURN	  
	  ENDIF
	  IF(AVEPET(IMON).LT.0.OR.AVEPET(IMON).GT.200)THEN
	    PRINT*,'Error in calculated PET: ',AVEPET(IMON),'mm'
	    CALL TRAPERROR(TRAPERR,METERR,1,STYPES,1,TOTSEQ)
	    RETURN	  
	  ENDIF
        IF(AVETEMP(IMON).LT.-20.OR.AVETEMP(IMON).GT.40)THEN
	    PRINT*,'Error in temperature data: ',AVETEMP(IMON),'deg.C'
	    CALL TRAPERROR(TRAPERR,METERR,1,STYPES,1,TOTSEQ)
	    RETURN	  
	  ENDIF
200   CONTINUE
      DO 300 IMON=1,12
        IF(AVETEMP(IMON).GT.0)GOTO 303
300   CONTINUE
	PRINT*,'Error in temperature data (all 0 or below)'
	CALL TRAPERROR(TRAPERR,METERR,1,STYPES,1,TOTSEQ)
	RETURN	  
303   CONTINUE
	GISOK=1
	END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CHECK_GIS_LU(GISOK,LUSEQ,LU1TOLU2,
     &                        PERSOIL,DOMSOIL,SOILID,
     &                        NUMCELLSIN20KM2,TRAPERR,ISERIES,NSEQ)
C
C Subroutine to check soil characteristics
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
      INTEGER MAXDEC				! Max.no.of decades included in calc.
	PARAMETER(MAXDEC=9)		
      INTEGER MAXDEC1				! Max.no.of decades included in calc.
	DATA MAXDEC1 /9/		
	INTEGER MAXGROW				! Max.no.of growing seasons allowed
	PARAMETER (MAXGROW=300)
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)						
	INTEGER MAXLU1				! Max.no.of land use types
	DATA MAXLU1 /6/
	INTEGER MAXSEQ				! Max.no of sequences included in GIS
	PARAMETER (MAXSEQ=MAXLU*MAXLU)						
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						

      INTEGER ICELL				! Local cell counter
      INTEGER IDEC				! Local decade counter
	INTEGER IS					! Local counter for soil series
	INTEGER LU2					! Land use changed to
C
C Variables passed to / from this routine
C 
      INTEGER DOMSOIL(CELLSTOGRID,MAXSERIES)	! 5 major soil series in the 1km 		
								! square cell x no.cells in 20km2 grid
	INTEGER GISOK				! IN:Check of GIS data
	INTEGER ISERIES				! IN:Soil series
	REAL LU1TOLU2(CELLSTOGRID,MAXDEC,MAXLU+2,MAXLU+2) ! IN:Fracn of LU1 to LU2 
													  ! in dec 
      INTEGER LUSEQ(MAXGROW)		! IN:Land use sequence
	INTEGER NSEQ				! IN:Sequence number
	INTEGER NUMCELLSIN20KM2		! IN:Number of cells in 20km2 cell
	REAL PERSOIL(CELLSTOGRID,MAXSERIES)		! % of each major soil type in the
								!	1km square cell
      INTEGER SOILID(CELLSTOGRID*MAXSERIES/2)	! IN:Soil series integer codes
      INTEGER TRAPERR(CELLSTOGRID*MAXSERIES/2,MAXSEQ) 
								! IN:Error trapping
								!		0=No error
								!		1=Error in NPP data
								!		2=Error in weather data
								!		3=Error in LU data
								!		4=Error in LU change data
								!		5=Soil error
								!		6=Equilibrium not found
      INTEGER LUERR
	DATA LUERR /3/				!		3=Error in LU data
      INTEGER LUCERR
	DATA LUCERR /4/				!		4=Error in LU change data
C
C Check for non-valid entries
C
      GISOK=0
C
C Look for this LU change anywhere in the 20km2 grid cell
C
      DO 100 ICELL=1,NUMCELLSIN20KM2,1
C
C If this is not a LU change (ie LU1->LU1), run this LU if LU change from LU1 
C occurs in any decade...
C
        IF(LUSEQ(1).EQ.LUSEQ(6))THEN
          DO 200 IDEC=1,MAXDEC1,1
	      DO 300 LU2=1,MAXLU1,1
	        IF(LU1TOLU2(ICELL,IDEC,LUSEQ(1),LU2).GT.0)THEN
C
C ...if LU change found, look for this soil series in the same grid cell.
C
                DO 400 IS=1,MAXSERIES
C
C ......if soil series found, do not mark as an error
C
                  IF(DOMSOIL(ICELL,IS).EQ.SOILID(ISERIES).AND.
     &               PERSOIL(ICELL,IS).GT.0)THEN
C Note check range of wetness class
                     GOTO 101
	            ENDIF
400             CONTINUE              
              ENDIF
300         CONTINUE
200       CONTINUE
        ENDIF
C
C If this is a LU change, look for this LU change occuring in the same grid cell 
C in any decade
C
	  DO 500 IDEC=1,MAXDEC1,1
     	    IF(LU1TOLU2(ICELL,IDEC,LUSEQ(1),LUSEQ(6)).GT.0)THEN
C
C ...if LU change found, look for this soil series in the same grid cell 
C
            DO 600 IS=1,MAXSERIES
C
C ......if soil series found, do not mark as an error
C
              IF(DOMSOIL(ICELL,IS).EQ.SOILID(ISERIES).AND.
     &           PERSOIL(ICELL,IS).GT.0.AND.
     &           PERSOIL(ICELL,IS).GT.0)GOTO 101
C Note check range of wetness class
600         CONTINUE              
          ENDIF
500     CONTINUE
100   CONTINUE
C
C If LU change or soil not found, mark as an error and leave without marking as OK
C so the simulation won't be run
C
	CALL TRAPERROR(TRAPERR,LUERR,ISERIES,ISERIES,NSEQ,NSEQ)
	RETURN
101   CONTINUE
C
C If LU change found on this soil series, mark as OK so simulation will run
C
	GISOK=1
	END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CHECK_GIS_SOIL(GISOK,SOILN,SOIL15,AMMN,AMMN15,SOILW,
     &			  DPMCARB0,DPMNIT0,DPMNLAB0,RPMCARB0,RPMNIT0,RPMNLAB0,
     &			  BCARB0,BNIT0,BNLAB0,HCARB0,HNIT0,HNLAB0,
     &              TRAPERR,ISERIES,NSEQ)
C
C Subroutine to check soil characteristics
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
	INTEGER MAXLAYER					! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)						
	INTEGER MAXLAYER1					! No.of layers in the soil profile
	DATA MAXLAYER1 /60/						
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXSEQ				! Max.no of sequences included in GIS
	PARAMETER (MAXSEQ=MAXLU*MAXLU)						
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						

	INTEGER IL							! Local layer counter
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
	INTEGER GISOK				! OUT:Check of GIS data
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN:N in soil humus (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN:N15 in soil biomass (kgN15/ha/layer)
	INTEGER ISERIES				! IN:Current soil series
	INTEGER NSEQ				! IN:Current LU sequence
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
      INTEGER TRAPERR(CELLSTOGRID*MAXSERIES/2,MAXSEQ) 
								! IN:Error trapping
								!		0=No error
								!		1=Error in NPP data
								!		2=Error in weather data
								!		3=Error in LU data
								!		4=Error in LU change data
								!		5=Soil error
								!		6=Equilibrium not found
	INTEGER SOILERR				
	DATA SOILERR /5/			!		5=Soil error
C
C Check for non-valid entries
C
      GISOK=0
      DO 100 IL=1,MAXLAYER
	  IF(AMMN(IL).LT.0)GOTO 101
	  IF(AMMN15(IL).LT.0)GOTO 101
	  IF(BCARB0(IL).LT.0)GOTO 101
        IF(BNIT0(IL).LT.0)GOTO 101
	  IF(BNLAB0(IL).LT.0)GOTO 101
	  IF(DPMCARB0(IL).LT.0)GOTO 101
	  IF(DPMNIT0(IL).LT.0)GOTO 101
	  IF(DPMNLAB0(IL).LT.0)GOTO 101
	  IF(HCARB0(IL).LT.0)GOTO 101
	  IF(HNIT0(IL).LT.0)GOTO 101
        IF(HNLAB0(IL).LT.0)GOTO 101
	  IF(RPMCARB0(IL).LT.0)GOTO 101
	  IF(RPMNIT0(IL).LT.0)GOTO 101
	  IF(RPMNLAB0(IL).LT.0)GOTO 101
	  IF(SOIL15(IL).LT.0)GOTO 101
	  IF(SOILN(IL).LT.0)GOTO 101
	  IF(SOILW(IL).LT.0)GOTO 101
	  GOTO 102
101     CONTINUE
	  CALL TRAPERROR(TRAPERR,SOILERR,ISERIES,ISERIES,NSEQ,NSEQ)
	  RETURN
102     CONTINUE
100   CONTINUE
	GISOK=1
	END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CLOSECHAN2()
C
C Subroutine to close channels opened for site specific calculations
C
      IMPLICIT NONE
	CLOSE(30)
	CLOSE(31)
	CLOSE(32)
	CLOSE(33)
	CLOSE(34)
	CLOSE(35)
	CLOSE(36)
	CLOSE(37)  ! Summary results file
	CLOSE(61)
	CLOSE(62)
	CLOSE(63)
	CLOSE(64)
	CLOSE(65)
	CLOSE(66)
	CLOSE(67)
      END
C
C-----------------------------------------------------------
C
      SUBROUTINE GET_GIS_CHANGE(REPORTDEPTH,ISERIES,NSEQ,IWC,
     &                          CCHANGE10,CO2C10,CH4C10,LU1,LU2,IDEC,
     &                          DPMCARB0,RPMCARB0,BCARB0,HCARB0)
C
C Subroutine to get change in carbon content over the decade
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
      INTEGER MAXDEC				! Max.no.of decades included in calculation
	PARAMETER(MAXDEC=9)		
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXSEQ				! Max.no of sequences included in GIS
	PARAMETER (MAXSEQ=MAXLU*MAXLU)						
	INTEGER MAXWC				! Max.no of wetness classes
      PARAMETER (MAXWC=6)

	INTEGER IL					! Local layer counter
	INTEGER IS					! Local soil series
	REAL SOMC
	INTEGER REPORTIL			! Layer reported to
	REAL TEMP					! Temporary real used in calculation
C
C Variables passed from calling subroutine
C
	REAL BCARB0(MAXLAYER)		! IN:C in soil biomass (kgC/ha/layer)
	REAL CCHANGE10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ) ! OUT:Change 
								! in C content over the decade (kgC/ha/decade)
	REAL CH4C10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! IN:CH4-C 
								! emitted over the decade (kgC/ha/decade)
	REAL CO2C10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! IN:CO2-C  
								! emitted over the decade (kgC/ha/decade)
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	INTEGER IDEC				! IN:Decade counter
	INTEGER ISERIES				! IN:Counter for soil series
	INTEGER IWC					! IN: Code for this wetness classes
	INTEGER LU1					! Land use code of first land use type
	INTEGER LU2					! Land use code of land use type changed to
	INTEGER NSEQ				! Counter for number of land use sequences 
	REAL REPORTDEPTH			! IN:Depth of reporting (cm)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
C
C Save values at the start of the decade 
C
      REPORTIL=REPORTDEPTH*MAXLAYER1/MAXDEPTH
      IS=ISERIES
	SOMC=0
	DO 100 IL=1,REPORTIL
	  SOMC=SOMC+DPMCARB0(IL)+RPMCARB0(IL)+BCARB0(IL)+HCARB0(IL)
100   CONTINUE
      CCHANGE10(IDEC,IS,IWC,NSEQ)=SOMC
C
C Save changes at the end of the decade
C
	IF(IDEC.GT.1)THEN
        CCHANGE10(IDEC-1,IS,IWC,NSEQ)=SOMC-CCHANGE10(IDEC-1,IS,IWC,NSEQ)
	  TEMP=CH4C10(IDEC-1,IS,IWC,NSEQ)-CCHANGE10(IDEC-1,IS,IWC,NSEQ)
	  CO2C10(IDEC-1,IS,IWC,NSEQ)=TEMP
      ENDIF
C
C Leave GET_GIS_CHANGE
C
      END
C
C-----------------------------------------------------------
C
      SUBROUTINE GET_GIS_EMISSIONS(REPORTDEPTH,CH4,GNN2O,GPNN2O,GDN2O,
     &                             CH4C10,N2ON10,ISERIES,NSEQ,IDEC,IWC)

C
C Subroutine to get emissions over the decade
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
      INTEGER MAXDEC				! Max.no.of decades included in calculation
	PARAMETER(MAXDEC=9)		
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXSEQ				! Max.no of sequences included in GIS
	PARAMETER (MAXSEQ=MAXLU*MAXLU)				
	INTEGER MAXWC				! Max.no of wetness classes
      PARAMETER (MAXWC=6)
			
	REAL FTEMP					! Local temporary real number
	INTEGER IL					! Local layer counter
	INTEGER REPORTIL			! Reporting depth
C
C Variables passed from calling subroutine
C
	REAL CH4(MAXLAYER)			! IN:CH4 emitted (kgC/ha/layer)
	REAL CH4C10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! IN:CH4-C  
								! emitted over the decade (kgC/ha/decade)
	REAL GNN2O					! IN:N2O lost by nitrification (kgN/ha)
	REAL GPNN2O					! IN:N2O lost by part.nitrification (kgN/ha)
	REAL GDN2O					! IN:N2O lost by denitrification (kgN/ha)
	INTEGER IDEC				! IN:Decade counter
	INTEGER ISERIES				! IN:Counter for soil series
	INTEGER IWC					! IN: Counter for wetness classes
	REAL N2ON10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! IN/OUT:N2O-N 
								! emitted over the decade (kgN/ha/decade)
	INTEGER NSEQ				! IN:Counter for number of land use sequences 
	REAL REPORTDEPTH			! IN:Depth of reporting (cm)
C
C Add up all emissions
C
      REPORTIL=REPORTDEPTH*MAXLAYER1/MAXDEPTH
	FTEMP=N2ON10(IDEC,ISERIES,IWC,NSEQ)
      N2ON10(IDEC,ISERIES,IWC,NSEQ)=FTEMP+GNN2O+GPNN2O+GDN2O
	DO 100 IL=1,REPORTIL
	  FTEMP=CH4C10(IDEC,ISERIES,IWC,NSEQ)
	  CH4C10(IDEC,ISERIES,IWC,NSEQ)=FTEMP+CH4(IL)
100   CONTINUE
C
C Leave GET_GIS_EMISSIONS
C
      END

C
C-----------------------------------------------------------
C
	SUBROUTINE GETFUTWEATHER(IK,IYEAR,LHARV,
     &                         FUTRAIN,FUTPET,FUTTEMP,
     &                         RAIN,EVAPW,SOILTEMP,AIRTEMP)
C
C Subroutine to get this months weather for GIS run
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
      INTEGER MAXMET				! No of met.years
	PARAMETER (MAXMET=90)

	INTEGER IL					! Local counter for months
	INTEGER IMON				! Current month
	REAL*8 RMON					! Working for current month
C
C Variables passed from calling subroutine
C
	INTEGER IK					! IN:No.timesteps from prev.harv. to current
	INTEGER IYEAR				! IN: Year number in simulation (assume starts at 1990) 
	INTEGER LHARV				! IN:No.timesteps from 01/01 to prev. harvest
      REAL FUTRAIN(MAXMET,12)		! IN:Future monthly rainfall(mm/month)
      REAL FUTPET(MAXMET,12)		! IN:Future monthly PET (mm/month)
      REAL FUTTEMP(MAXMET,12)		! IN:Future monthly average 
								!	air temp (deg.C)
	REAL RAIN					! IN: Rainfall (mm/timestep)
	REAL EVAPW					! IN: Potential evap. (mm/timestep)
	REAL AIRTEMP				! IN: Air temperature (deg.C/timestep)	
	REAL SOILTEMP(MAXLAYER)		! IN: Soil temperature (deg.C/timestep)
C
C Get the current month
C
	IMON=AINT((IK+LHARV)/12.)
	RMON=(REAL(IK+LHARV)/12.)
	IMON=NINT(12.*(RMON-IMON))
	IF(IMON.EQ.0)IMON=12
C
C Get weather for current month
C
      RAIN=FUTRAIN(IYEAR,IMON)
	EVAPW=FUTPET(IYEAR,IMON)
	AIRTEMP=FUTTEMP(IYEAR,IMON)
C
C Set soil temperature from air temperature
C Assume soil temperature is same thoughout profile
C
      DO 100 IL=1,MAXLAYER1
	  SOILTEMP(IL)=AIRTEMP
100   CONTINUE
      END
C
C-----------------------------------------------------------
C
	SUBROUTINE GETFUTWEATHER_JULES(IK,ILU,IYEAR,LHARV,
     &                              FUTNPP,FUTRAIN,FUTPET,FUTTEMP,
     &                              GRIDNPP,RAIN,EVAPW,SOILTEMP,AIRTEMP)
C
C Subroutine to get this months weather for future GIS run and translate NPP into an adjusted PI
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
      INTEGER MAXMET				! No of met.years
	PARAMETER (MAXMET=90)
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER IMON				! Current month
	REAL*8 RMON					! Working for current month
	INTEGER IL					! Local counter for months
C
C Variables passed from calling subroutine
C
	INTEGER IK					! IN(CALL):No.timesteps from prev.harv. to current
	INTEGER ILU					! IN(CALL):Land use code
	INTEGER IYEAR				! IN(CALL): Year number in simulation (assume starts at 1990) 
	INTEGER LHARV				! IN(CALL):No.timesteps from 01/01 to prev. harvest
      REAL FUTNPP(MAXMET,12,MAXLU)! IN(CALL):Future monthly NPP (kgC/ha/yr)
      REAL FUTRAIN(MAXMET,12)		! IN(CALL):Future monthly rainfall(mm/month)
      REAL FUTPET(MAXMET,12)		! IN(CALL):Future monthly PET (mm/month)
      REAL FUTTEMP(MAXMET,12)		! IN(CALL):Future monthly average 
								!	air temp (deg.C)
	REAL GRIDNPP				! OUT(CALL):Average NPP in the 1km cell (kgC/ha/yr)
	REAL RAIN					! OUT(CALL): Rainfall (mm/timestep)
	REAL EVAPW					! OUT(CALL): Potential evap. (mm/timestep)
	REAL AIRTEMP				! OUT(CALL): Air temperature (deg.C/timestep)	
	REAL SOILTEMP(MAXLAYER)		! OUT(CALL): Soil temperature (deg.C/timestep)
C
C Get the current month
C
	IMON=AINT((IK+LHARV)/12.)
	RMON=(REAL(IK+LHARV)/12.)
	IMON=NINT(12.*(RMON-IMON))
	IF(IMON.EQ.0)IMON=12
C
C Get weather for current month
C
      GRIDNPP=FUTNPP(IYEAR,IMON,ILU)
      RAIN=FUTRAIN(IYEAR,IMON)
	EVAPW=FUTPET(IYEAR,IMON)
	AIRTEMP=FUTTEMP(IYEAR,IMON)
C
C Set soil temperature from air temperature
C Assume soil temperature is same thoughout profile
C
      DO 100 IL=1,MAXLAYER1
	  SOILTEMP(IL)=AIRTEMP
100   CONTINUE
      END
C
C-----------------------------------------------------------
C
      SUBROUTINE GETWEATHER_GIS(IK,LHARV,AVERAIN,AVEPET,AVETEMP, 
     &                           RAIN,EVAPW,SOILTEMP,AIRTEMP)
C
C Subroutine to get this months weather for GIS run
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER IMON				! Current month
	REAL*8 RMON					! Working for current month
	INTEGER IL					! Local counter for months
C
C Variables passed from calling subroutine
C
	INTEGER IK					! IN:No.timesteps from prev.harv. to current
	INTEGER LHARV				! IN:No.timesteps from 01/01 to prev. harvest
      REAL AVERAIN(12)			! Long term average rainfall(mm/month)
      REAL AVEPET(12)				! Long term average PET (mm/month)
      REAL AVETEMP(12)			! Long term average monthly average 
	REAL RAIN					! IN: Rainfall (mm/timestep)
	REAL EVAPW					! IN: Potential evap. (mm/timestep)
	REAL AIRTEMP				! IN: Air temperature (deg.C/timestep)	
	REAL SOILTEMP(MAXLAYER)		! IN: Soil temperature (deg.C/timestep)
C
C Get the current month
C
	IMON=AINT((IK+LHARV)/12.)
	RMON=(REAL(IK+LHARV)/12.)
	IMON=NINT(12.*(RMON-IMON))
	IF(IMON.EQ.0)IMON=12
C
C Get weather for current month
C
      RAIN=AVERAIN(IK)
	EVAPW=AVEPET(IK)
	AIRTEMP=AVETEMP(IK)
C
C Set soil temperature from air temperature
C Assume soil temperature is same thoughout profile
C
      DO 100 IL=1,MAXLAYER1
	  SOILTEMP(IL)=AIRTEMP
100   CONTINUE
      END
C
C-----------------------------------------------------------
C
	SUBROUTINE INIT_GIS_RES(CO2C10,CH4C10,N2ON10,ISERIES,NSEQ,IWC)
C
C Subroutine to get this months weather for GIS run
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
      INTEGER MAXDEC				! Max.no.of decades included in calculation
	PARAMETER(MAXDEC=9)		
	INTEGER MAXDEC1				! Max.no.of decades included in calculation
	DATA MAXDEC1 /9/
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXSEQ				! Max.no of sequences included in GIS
	PARAMETER (MAXSEQ=MAXLU*MAXLU)						
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXWC				! Max.no of wetness classes
      PARAMETER (MAXWC=6)

	INTEGER IDEC				! Local decade counter
C
C Variables passed from calling subroutine
C
	REAL CO2C10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! OUT:CO2-C 
								! emitted over the decade (kgC/ha/decade)
	REAL CH4C10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! OUT:CH4-C  
								! emitted over the decade (kgC/ha/decade)
	INTEGER ISERIES				! IN:Counter for soil series
	INTEGER IWC					! IN:Counter for wetness classes
	REAL N2ON10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! OUT:N2O-N 
								! emitted over the decade (kgN/ha/decade)
	INTEGER NSEQ				! Counter for number of land use sequences 
C
C Initialise to zero
C
      DO 100 IDEC=1,MAXDEC1
        CO2C10(IDEC,ISERIES,IWC,NSEQ)=0
	  CH4C10(IDEC,ISERIES,IWC,NSEQ)=0
	  N2ON10(IDEC,ISERIES,IWC,NSEQ)=0
100   CONTINUE
	END
C
C--------------------------------------------------------
C
	SUBROUTINE OPENCHAN_GIS(GIS_INDATA,REPORTDEPTH,RUNFUTMET)
C
C Subroutine to open channels
C
      IMPLICIT NONE
C
C Variables local to this routine
C
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXSOMLAY			! Max.no.of soil layers entered			
	PARAMETER (MAXSOMLAY=10)

	INTEGER CODECHAN(MAXLU)		! Channel for soil codes 
      INTEGER DATACHAN			! Channel for GIS data
      INTEGER FRESCHAN			! Channel for results of future met.data outputs
      INTEGER FUTCHAN				! Channel for file of future met.data inputs
	CHARACTER*80 FUTMETFILE		! File of future met.data
	CHARACTER*80 FUTNPPFILE		! File of future NPP data 
	CHARACTER*80 FUTRAINFILE	! File of future rain.data 
	CHARACTER*80 FUTTEMPFILE	! File of future temp.data 
	CHARACTER*80 FUTRESFILE		! File of results using future met.data
	INTEGER GISCHAN				! Channel for file specifying GIS inputs
      CHARACTER*80 GISDATA		! File for GIS data
	INTEGER GRIDCHAN			! Channel for file for 20km2 grid outputs
      CHARACTER*80 GRIDFILE		! File for 20km2 grid output
	INTEGER ICHAN				! Channel counter
	INTEGER IL					! Layer counter
	INTEGER ILU1				! Land use counter
	INTEGER ILU2				! Land use counter
	INTEGER ISERIES				! Soil series counter
	INTEGER ISOMLAY				! SOM soil layer counter
      INTEGER LUCHAN				! Channel for LU data
      CHARACTER*80 LUDATA			! File for LU data
	CHARACTER*80 METID			! File for location of met.stations
      INTEGER MITICHAN(MAXLU*MAXLU)		! Channel for results of files for calculation of mitigation impacts
	CHARACTER*80 MITIFILE(MAXLU*MAXLU)	! File for calculation of mitigation options
	INTEGER NF					! Counter
      INTEGER NPPCHAN				! Channel for file of future NPP data inputs
	INTEGER NUMFILES			! Number of soil code files used (non-SSKIB=1; SSKIB=4)
      INTEGER OFRESCHAN			! Channel for results of future met.data outputs organic soils
	CHARACTER*80 OFUTRESFILE	! File of results using future met.data organic soils
	INTEGER OGRIDCHAN			! Channel for file for 20km2 grid outputs - organic soil
      CHARACTER*80 OGRIDFILE		! File for 20km2 grid output - organic soil
	INTEGER OOUTCHAN			! Channel for GIS output - organic soil
	CHARACTER*80 OOUTFILE		! File for GIS output - organic soil
	CHARACTER*80 OUTFILE		! File for GIS output
	INTEGER OUTCHAN				! Channel for GIS output
	CHARACTER*80 PETFILE		! File for PET data
	CHARACTER*80 RAINFILE		! File for rainfall data
      INTEGER RAINCHAN			! Channel for file of future rain data inputs
	INTEGER STARTCHAN			! Channel for outputting start data
      INTEGER TEMPCHAN			! Channel for file of future temp.data inputs
	CHARACTER*10 TITLE1(MAXSERIES,MAXLU,MAXSOMLAY)	! Titles for depth columns
	CHARACTER*10 TITLE2(MAXSERIES,MAXLU,MAXSOMLAY)	! Titles for soil C columns
      CHARACTER*10 TITLE3(MAXSERIES)					! Titles for percentage soil columns
      CHARACTER*10 TITLE4(MAXLU)						! Titles for land use columns
	CHARACTER*80 SOILCODES(6)	! Files for soil codes 
	CHARACTER*100 TEMP
C
C Variables passed to/from this routine
C
      INTEGER GIS_INDATA			! IN(CALL): Type of data is used in this run
	 INTEGER GIS_NONSSKIB		! Scottish non-SSKIB
	 INTEGER GIS_SSKIB			! Scottish SSKIB
	 INTEGER GIS_JULES			! JULES data
	 DATA GIS_NONSSKIB,GIS_SSKIB,GIS_JULES /1,2,3/
	REAL REPORTDEPTH			! OUT(CALL):Depth of reporting (cm)
	INTEGER RUNFUTMET			! OUT(CALL):Integer code for future met.run
								!	  0=No, 1=Yes
C
C Channels
C
      DATACHAN=30
      LUCHAN=32
	OUTCHAN=37
	GRIDCHAN=38
	OOUTCHAN=47
	OGRIDCHAN=48
	GISCHAN=39
	NPPCHAN=55
	TEMPCHAN=56
	RAINCHAN=57
	CODECHAN(1)=61
	CODECHAN(2)=62
	CODECHAN(3)=63
	CODECHAN(4)=64
	FUTCHAN=65
	FRESCHAN=66
	OFRESCHAN=67
	STARTCHAN=68
	DO 500 ILU1=1,MAXLU
	  DO 600 ILU2=1,MAXLU
	   IF(ILU1.NE.ILU2)THEN
	     ICHAN=((ILU1-1)*MAXLU)+ILU2
	     MITICHAN(ICHAN)=OFRESCHAN+ICHAN
	   ENDIF
600     CONTINUE
500   CONTINUE
C
C Open input files
C
C File input
C Get names of input files from file GNAMES.DAT
C
	IF(GIS_INDATA.EQ.GIS_NONSSKIB)THEN
	  NUMFILES=1
	ELSEIF(GIS_INDATA.EQ.GIS_SSKIB)THEN
	  NUMFILES=MAXLU
	ELSEIF(GIS_INDATA.EQ.GIS_JULES)THEN ! JULES runs currently assumed to use NON-SSKIB format for soils data 
	  NUMFILES=1
	ENDIF
      OPEN(GISCHAN,FILE='GNAMES.DAT',STATUS='OLD',ERR=11)
	READ(GISCHAN,10)GISDATA
	READ(GISCHAN,10)LUDATA
	DO 100 NF=1,NUMFILES
	  READ(GISCHAN,10)SOILCODES(NF)
100   CONTINUE
	READ(GISCHAN,10)OUTFILE
	READ(GISCHAN,10)GRIDFILE
	READ(GISCHAN,10)OOUTFILE
	READ(GISCHAN,10)OGRIDFILE
	READ(GISCHAN,*)REPORTDEPTH
	IF(REPORTDEPTH.LE.0)REPORTDEPTH=300
C
C Type of simulation
C
      IF(GIS_INDATA.EQ.GIS_JULES)THEN
	  READ(GISCHAN,10,END=12)FUTRAINFILE
	  READ(GISCHAN,10,END=12)FUTTEMPFILE
	  READ(GISCHAN,10,END=12)FUTNPPFILE
	ELSE
	  READ(GISCHAN,10,END=12)FUTMETFILE
	ENDIF
	READ(GISCHAN,10,END=12)FUTRESFILE
	READ(GISCHAN,10,END=12)OFUTRESFILE
12    CONTINUE
	CLOSE(GISCHAN)
C
C Soils data
C
      OPEN(DATACHAN,FILE=GISDATA,STATUS='OLD',ERR=1111)
      GOTO 1102
1111  CONTINUE
      PRINT*,'Cannot open file ',GISDATA,' Please check GNAMES.DAT'
	PRINT*,'...Press any key to continue'
	READ(*,*)
	STOP
1102  CONTINUE
      READ(DATACHAN,*,ERR=1103,END=1103)TEMP
      GOTO 1104
1103  CONTINUE
      PRINT*,'Error in reading file ',GISDATA
	PRINT*,'Last line read = ',TEMP
	PRINT*,'Press any key to continue...'
	READ(*,*)
	STOP
1104  CONTINUE
C
C LU Data
C
      OPEN(LUCHAN,FILE=LUDATA,STATUS='OLD',ERR=1222)
      GOTO 1202
1222  CONTINUE
      PRINT*,'Cannot open file ',LUDATA,' Please check GNAMES.DAT'
      PRINT*,'...Press any key to continue'
	READ(*,*)
	STOP
1202  CONTINUE
      READ(LUCHAN,*)
C
C Soil codes 
C
      DO 2100 NF=1,NUMFILES
        OPEN(CODECHAN(NF),FILE=SOILCODES(NF),STATUS='OLD',ERR=1333)
        GOTO 1303
1333	  CONTINUE
        PRINT*,'Cannot open file ',SOILCODES(NF),
     &         ' Please check GNAMES.DAT'
        PRINT*,'...Press any key to continue'
	  READ(*,*)
	  STOP
1303    CONTINUE
        READ(CODECHAN(NF),*)
2100  CONTINUE
C
C Future met.data
C
      RUNFUTMET=1
      IF(GIS_INDATA.EQ.GIS_JULES)THEN
        OPEN(RAINCHAN,FILE=FUTRAINFILE,STATUS='OLD',ERR=1444)
        OPEN(TEMPCHAN,FILE=FUTTEMPFILE,STATUS='OLD',ERR=1444)
        OPEN(NPPCHAN,FILE=FUTNPPFILE,STATUS='OLD',ERR=1444)
	ELSE
        OPEN(FUTCHAN,FILE=FUTMETFILE,STATUS='OLD',ERR=1444)
	ENDIF
        GOTO 1404
1444	  CONTINUE
	  RUNFUTMET=0
1404    CONTINUE
C
C Open 1km2 output file
C
      OPEN(OUTCHAN,FILE=OUTFILE,STATUS='UNKNOWN',ERR=1888)
      GOTO 1802
1888  CONTINUE
      PRINT*,'Cannot open file ',OUTFILE,' Please check GNAMES.DAT'
      PRINT*,'...Press any key to continue'
	READ(*,*)
	STOP
1802  CONTINUE
C
C Open 20km2 output file
C
      OPEN(GRIDCHAN,FILE=GRIDFILE,STATUS='UNKNOWN',ERR=1999)
      GOTO 1902
1999  CONTINUE
      PRINT*,'Cannot open file ',GRIDFILE,' Please check GNAMES.DAT'
      PRINT*,'...Press any key to continue'
	READ(*,*)
	STOP
1902  CONTINUE
C
C Open 1km2 output file - organic soils
C
      OPEN(OOUTCHAN,FILE=OOUTFILE,STATUS='UNKNOWN',ERR=2001)
      GOTO 2002
2001  CONTINUE
      PRINT*,'Cannot open file ',OOUTFILE,' Please check GNAMES.DAT'
      PRINT*,'...Press any key to continue'
	READ(*,*)
	STOP
2002  CONTINUE
C
C Open 20km2 output file
C
      OPEN(OGRIDCHAN,FILE=OGRIDFILE,STATUS='UNKNOWN',ERR=2101)
      GOTO 2102
2101  CONTINUE
      PRINT*,'Cannot open file ',OGRIDFILE,' Please check GNAMES.DAT'
      PRINT*,'...Press any key to continue'
	READ(*,*)
	STOP
2102  CONTINUE
C
C Open future climate results file
C
      IF(RUNFUTMET.EQ.1)THEN
        OPEN(FRESCHAN,FILE=FUTRESFILE,STATUS='UNKNOWN',ERR=2222)
        OPEN(OFRESCHAN,FILE=OFUTRESFILE,STATUS='UNKNOWN',ERR=2222)
        GOTO 2202
2222    CONTINUE
        PRINT*,'Cannot open res.file ',FUTRESFILE,
     &         ' Please check GNAMES.DAT'
        PRINT*,'...Press any key to continue'
	  READ(*,*)
	  STOP
2202    CONTINUE
      ENDIF
C
C Open results files for calculating mitigation options
C
C LU1 = Arable
	MITIFILE(1)='MIT_A2A.OUT'
      MITIFILE(2)='MIT_A2P.OUT'
      MITIFILE(3)='MIT_A2W.OUT'
      MITIFILE(4)='MIT_A2N.OUT'
	MITIFILE(5)='MIT_A2M.OUT'
	MITIFILE(6)='MIT_A2S.OUT'
C LU1 = Grass
      MITIFILE(7)='MIT_P2A.OUT'
      MITIFILE(8)='MIT_P2P.OUT'
      MITIFILE(9)='MIT_P2W.OUT'
      MITIFILE(10)='MIT_P2N.OUT'
      MITIFILE(11)='MIT_P2M.OUT'
      MITIFILE(12)='MIT_P2S.OUT'
C LU1 = Forestry
      MITIFILE(13)='MIT_W2A.OUT'
      MITIFILE(14)='MIT_W2P.OUT'
      MITIFILE(15)='MIT_W2W.OUT'
      MITIFILE(16)='MIT_W2N.OUT'
      MITIFILE(17)='MIT_W2M.OUT'
      MITIFILE(18)='MIT_W2S.OUT'
C LU1 = Seminatural
      MITIFILE(19)='MIT_N2A.OUT'
      MITIFILE(20)='MIT_N2P.OUT'
      MITIFILE(21)='MIT_N2W.OUT'
      MITIFILE(22)='MIT_N2N.OUT'
      MITIFILE(23)='MIT_N2M.OUT'
      MITIFILE(24)='MIT_N2S.OUT'
C LU1 = Miscanthus
      MITIFILE(25)='MIT_M2A.OUT'
      MITIFILE(26)='MIT_M2P.OUT'
      MITIFILE(27)='MIT_M2W.OUT'
      MITIFILE(28)='MIT_M2N.OUT'
      MITIFILE(29)='MIT_M2M.OUT'
      MITIFILE(30)='MIT_M2S.OUT'
C LU1 = SRC
      MITIFILE(31)='MIT_S2A.OUT'
      MITIFILE(32)='MIT_S2P.OUT'
      MITIFILE(33)='MIT_S2W.OUT'
      MITIFILE(34)='MIT_S2N.OUT'
      MITIFILE(35)='MIT_S2M.OUT'
      MITIFILE(36)='MIT_S2S.OUT'
C      
	DO 2300 ILU1=1,MAXLU
	  DO 2400 ILU2=1,MAXLU
	   IF(ILU1.NE.ILU2)THEN
	     ICHAN=((ILU1-1)*MAXLU)+ILU2
	     OPEN(MITICHAN(ICHAN),FILE=MITIFILE(ICHAN),STATUS='UNKNOWN')
	   ENDIF
2400    CONTINUE
2300  CONTINUE
C
C Miss out interactive data input
C
      GOTO 22
C
C Interactive input
C GIS data
C
11    CONTINUE
      PRINT*,'GIS DATA'
	PRINT*,'========================================================='
	PRINT*,'Format:Line 1: Title line'
	PRINT*,'       Line 2: 20km Square ID (10 characters)' 
	PRINT*,'               1km Square ID (real number)' 
	PRINT*,'               National grid easting x 1000'
	PRINT*,'               National grid northing x 1000'
	PRINT*,'               Area'
	PRINT*,'               Code number of dominant soils series 1 '
	PRINT*,'               Percentage of soils series 1 in 1km cell'
	PRINT*,'               Code number of dominant soils series 2'
	PRINT*,'               Percentage of soils series 2 in 1km cell'
	PRINT*,'               Code number of dominant soils series 3'
	PRINT*,'               Percentage of soils series 3 in 1km cell'
	PRINT*,'               Code number of dominant soils series 4'
	PRINT*,'               Percentage of soils series 4 in 1km cell'
	PRINT*,'               Code number of dominant soils series 5'
	PRINT*,'               Percentage of soils series 5 in 1km cell'
	PRINT*,'               Percentage of other soils in the 1km cell'
	PRINT*,'               Percentage of organic soil in the 1km cell'
	PRINT*,'               Percentage of organomin.soil in 1km cell'
	PRINT*,'               NPP (UK) kg C / ha x 1000'
      PRINT*,'               Long term average rainfall (mm) /Jan'
      PRINT*,'               Long term average rainfall (mm) /Feb'
      PRINT*,'               Long term average rainfall (mm) /...'
      PRINT*,'               Long term average rainfall (mm) /Dec'
      PRINT*,'               Long term average temp. (deg.C) /Jan'
      PRINT*,'               Long term average temp. (deg.C) /Feb'
      PRINT*,'               Long term average temp. (deg.C) /...'
      PRINT*,'               Long term average temp. (deg.C) /Dec'
      PRINT*,'               Longitude'
      PRINT*,'               Latitude'
	PRINT*,'               ... (ONE LINE FOR EACH CELL ',
     &                           '- SPACE SEPARATED FREE FORMAT)'
	PRINT*,'========================================================='
	PRINT*,'Enter name of soils data file: '
101   CONTINUE
	READ(*,10)GISDATA
10    FORMAT(A40)
      OPEN(DATACHAN,FILE=GISDATA,STATUS='OLD',ERR=111)
      GOTO 102
111   CONTINUE
      PRINT*,'Cannot open file ',GISDATA,' Please enter other file: '
	GOTO 101
102   CONTINUE
C
C Miss out title line
C
      READ(DATACHAN,*)
	PRINT*,'========================================================='
C
C LU data
C
      PRINT*,'LU DATA'
	PRINT*,'========================================================='
	PRINT*,'Format:Line 1: Title line'
	PRINT*,'       Line 2: 20km Square ID (10 characters)' 
	PRINT*,'               1km Square ID (real number)' 
	PRINT*,'               National grid easting x 1000'
	PRINT*,'               National grid northing x 1000'
	PRINT*,'               Area'
	PRINT*,'               Fraction cell under arable at start'
	PRINT*,'               Dummy fraction (not used)'
	PRINT*,'               Fraction cell under seminat at start'
	PRINT*,'               Dummy fraction (not used)'
	PRINT*,'               Fraction cell under grass at start'
	PRINT*,'               Fraction cell under sea at start'
	PRINT*,'               Fraction cell under forestry at start'
	PRINT*,'               Fraction cell under urban at start'
	PRINT*,'               Fraction cell under water at start'
	PRINT*,'               Fraction cell under non-CORINE at start'
	PRINT*,'               Fraction cell under miscanthus at start'
	PRINT*,'               Fraction cell under SRC at start'
	PRINT*,'               % change in 1950s to forestry from'
	PRINT*,'											forestry'
	PRINT*,'											natural'
	PRINT*,'											grassland'
	PRINT*,'											arable'
	PRINT*,'											miscanthus'
	PRINT*,'											SRC'
	PRINT*,'               % change in 1950s to natural from'
	PRINT*,'											forestry...'
	PRINT*,'               % change in 1950s to grassland from'
	PRINT*,'											forestry...'
	PRINT*,'               % change in 1950s to arable from'
	PRINT*,'											forestry...'
	PRINT*,'               % change in 1950s to urban from'
	PRINT*,'											forestry...'
	PRINT*,'               % change in 1950s to other from'
	PRINT*,'											forestry...'
	PRINT*,'               % change in 1960s to forestry from'
	PRINT*,'											forestry'
	PRINT*,'               % change in 1960s ...'
	PRINT*,'               % change in 1970s to forestry from'
	PRINT*,'											forestry'
	PRINT*,'               % change in 1970s ...'
	PRINT*,'               % change in 1980s to forestry from'
	PRINT*,'											forestry'
	PRINT*,'               % change in 1980s ...'
	PRINT*,'               % change in 1990s to forestry from'
	PRINT*,'											forestry'
	PRINT*,'               % change in 1990s ...'
	PRINT*,'               % change in 2000s to forestry from'
	PRINT*,'											forestry'
	PRINT*,'               % change in 2000s ...'
	PRINT*,'               ... (ONE LINE FOR EACH CELL ',
     &                           '- SPACE SEPARATED FREE FORMAT)'
      OPEN(LUCHAN,FILE=LUDATA,STATUS='OLD',ERR=1222)
      GOTO 1202
222   CONTINUE
      PRINT*,'Cannot open file ',LUDATA,' Please check GNAMES.DAT'
      PRINT*,'...Press any key to continue'
	READ(*,*)
	STOP
202   CONTINUE
C
C Miss out title lines
C
      READ(LUCHAN,*)
C
C Soil codes
C Scottish non-SSKIB data
C
	IF(GIS_INDATA.EQ.GIS_NONSSKIB.OR.GIS_INDATA.EQ.GIS_JULES)THEN
      PRINT*,'SOIL CODES'
	PRINT*,'========================================================='
	PRINT*,'Format:Line 1: Title line'
	PRINT*,'       Line 2: Soil code number                          ' 
	PRINT*,'               Arable (0-30cm): % C by mass              '
	PRINT*,'			         Amount C in kg*10^6/km^2 '
	PRINT*,'			         % clay                   '
	PRINT*,'			         % silt                   '
	PRINT*,'			         % sand                   '
	PRINT*,'			         Bulk density (g/cm^3)    '
	PRINT*,'               Arable (30-100cm):% C by mass             '
	PRINT*,'			         Amount C in kg*10^6/km^2'
	PRINT*,'			         % clay                  '
	PRINT*,'			         % silt                  '
	PRINT*,'			         % sand                  '
	PRINT*,'			         Bulk density (g/cm^3)   '
	PRINT*,'               Grass (0-30cm):  % C by mass              '
	PRINT*,'			         Amount C in kg*10^6/km^2 '
	PRINT*,'			         % clay                   '
	PRINT*,'			         % silt                   '
	PRINT*,'			         % sand                   '
	PRINT*,'			         Bulk density (g/cm^3)    '
	PRINT*,'               Grass (30-100cm): % C by mass             '
	PRINT*,'			         Amount C in kg*10^6/km^2'
	PRINT*,'			         % clay                  '
	PRINT*,'			         % silt                  '
	PRINT*,'			         % sand                  '
	PRINT*,'			         Bulk density (g/cm^3)   '
	PRINT*,'               SemiNat(0-30cm): % C by mass              '
	PRINT*,'			         Amount C in kg*10^6/km^2 '
	PRINT*,'			         % clay                   '
	PRINT*,'			         % silt                   '
	PRINT*,'			         % sand                   '
	PRINT*,'			         Bulk density (g/cm^3)    '
	PRINT*,'               SemiNat(30-100cm):% C by mass             '
	PRINT*,'			         Amount C in kg*10^6/km^2'
	PRINT*,'			         % clay                  '
	PRINT*,'			         % silt                  '
	PRINT*,'			         % sand                  '
	PRINT*,'			         Bulk density (g/cm^3)   '
	PRINT*,'               Woodland(0-30cm):% C by mass              '
	PRINT*,'			         Amount C in kg*10^6/km^2 '
	PRINT*,'			         % clay                   '
	PRINT*,'			         % silt                   '
	PRINT*,'			         % sand                   '
	PRINT*,'			         Bulk density (g/cm^3)    '
	PRINT*,'               Wood.(30-100cm):  % C by mass             '
	PRINT*,'			         Amount C in kg*10^6/km^2'
	PRINT*,'			         % clay                  '
	PRINT*,'			         % silt                  '
	PRINT*,'			         % sand                  '
	PRINT*,'			         Bulk density (g/cm^3)   '
	PRINT*,'               Miscanthus(0-30cm):% C by mass            '
	PRINT*,'			         Amount C in kg*10^6/km^2 '
	PRINT*,'			         % clay                   '
	PRINT*,'			         % silt                   '
	PRINT*,'			         % sand                   '
	PRINT*,'			         Bulk density (g/cm^3)    '
	PRINT*,'               Misc.(30-100cm):  % C by mass             '
	PRINT*,'			         Amount C in kg*10^6/km^2'
	PRINT*,'			         % clay                  '
	PRINT*,'			         % silt                  '
	PRINT*,'			         % sand                  '
	PRINT*,'			         Bulk density (g/cm^3)   '
	PRINT*,'               SRC(0-30cm):% C by mass              '
	PRINT*,'			         Amount C in kg*10^6/km^2 '
	PRINT*,'			         % clay                   '
	PRINT*,'			         % silt                   '
	PRINT*,'			         % sand                   '
	PRINT*,'			         Bulk density (g/cm^3)    '
	PRINT*,'               SRC(30-100cm):  % C by mass             '
	PRINT*,'			         Amount C in kg*10^6/km^2'
	PRINT*,'			         % clay                  '
	PRINT*,'			         % silt                  '
	PRINT*,'			         % sand                  '
	PRINT*,'			         Bulk density (g/cm^3)   '
	PRINT*,'               (SPACE SEPARATED FREE FORMAT)'
	PRINT*,'       Line 3: Soil code number                          ' 
	PRINT*,'               ... (ONE LINE FOR EACH SOIL CODE)         ' 
	PRINT*,'========================================================='
	PRINT*,'Enter name of soil codes file: '
304   CONTINUE
	READ(*,10)SOILCODES(1)
      OPEN(CODECHAN(1),FILE=SOILCODES(1),STATUS='OLD',ERR=305)
      GOTO 306
305   CONTINUE
      PRINT*,'Cannot open file ',SOILCODES(1),' Please enter other file'
	GOTO 304
306   CONTINUE
C
C Scottish SSKIB data
C
      ELSEIF(GIS_INDATA.EQ.GIS_SSKIB)THEN 
      PRINT*,'SOIL CODES'
	PRINT*,'========================================================='
	PRINT*,'Format:Line 1: Title line'
	PRINT*,'       Line 2: Soil code number' 
	PRINT*,'               Layer 1'
	PRINT*,'                     Top depth (cm)'
	PRINT*,'                     Bottom depth (cm)'
	PRINT*,'                     Thickness (cm)'
	PRINT*,'                     pH measured in water'
	PRINT*,'                     % C by mass'
	PRINT*,'			         Amount C in kg / ha'
	PRINT*,'			         % clay'
	PRINT*,'			         % silt'
	PRINT*,'			         % sand'
	PRINT*,'			         Bulk density (g/cm^3)'
	PRINT*,'			         Percent stones'
	PRINT*,'               Layer 2'
	PRINT*,'                     Top depth (cm)...'
	PRINT*,'               Layer 3'
	PRINT*,'                     Top depth (cm)...'
	PRINT*,'               Layer 4'
	PRINT*,'                     Top depth (cm)...'
	PRINT*,'               Layer 5'
	PRINT*,'                     Top depth (cm)...'
	PRINT*,'               Layer 6'
	PRINT*,'                     Top depth (cm)...'
	PRINT*,'               Layer 7'
	PRINT*,'                     Top depth (cm)...'
	PRINT*,'               Layer 8'
	PRINT*,'                     Top depth (cm)...'
	PRINT*,'               Layer 9'
	PRINT*,'                     Top depth (cm)...'
	PRINT*,'       Line 3: Soil code number' 
	PRINT*,'               ... (ONE LINE FOR EACH SOIL CODE)' 
	PRINT*,'========================================================='
	DO 300 NF=1,NUMFILES
	  IF(NF.EQ.1)PRINT*,'Enter name of the arable soil codes file: '
	  IF(NF.EQ.2)PRINT*,'Enter name of the grassland codes file: '
	  IF(NF.EQ.3)PRINT*,'Enter name of the forestry soil codes file: '
	  IF(NF.EQ.4)PRINT*,'Enter name of the semi-nat.soil codes file: '
301     CONTINUE
	  READ(*,10)SOILCODES(NF)
        OPEN(CODECHAN(NF),FILE=SOILCODES(NF),STATUS='OLD',ERR=333)
        GOTO 303
333     CONTINUE
        PRINT*,'Cannot open ',SOILCODES(NF),'. Please enter other file'
	  GOTO 301
303     CONTINUE
300   CONTINUE
      ENDIF
C
C Miss out title lines
C
      DO 307 NF=1,NUMFILES
        READ(CODECHAN(NF),*)
307   CONTINUE
	PRINT*,'========================================================='
C
C Open output files
C
      PRINT*,'OUTPUT FILE'
	PRINT*,'========================================================='
	PRINT*,'Enter name of output file: '
801   CONTINUE
	READ(*,10)OUTFILE
      OPEN(OUTCHAN,FILE=OUTFILE,STATUS='UNKNOWN',ERR=888)
      GOTO 802
888   CONTINUE
      PRINT*,'Cannot open file ',OUTFILE,' Please enter other file: '
	GOTO 801
802   CONTINUE
	PRINT*,'========================================================='
C
C Open output files
C
      PRINT*,'OUTPUT FILE'
	PRINT*,'========================================================='
	PRINT*,'Enter name of output file: '
901   CONTINUE
	READ(*,10)GRIDFILE
      OPEN(GRIDCHAN, FILE=GRIDFILE,STATUS='UNKNOWN',ERR=999)
      GOTO 902
999   CONTINUE
      PRINT*,'Cannot open file ',GRIDFILE,' Please enter other file: '
	GOTO 901
902   CONTINUE
	PRINT*,'========================================================='
C
C Files have been opened
C
22    CONTINUE
C
C Title lines for output file
C   
      WRITE(OUTCHAN,440)
      WRITE(OOUTCHAN,440)
440   FORMAT('20KM2_GRID      SQUID     EAST    NORTH GLEY ',
     &   ' PlantC_Input  Land_Use  ',
     &   'Arable_C    Grassland_C    Forest_C    Semi-Natural_C ',
     &   'Arable_D    Grassland_D    Forest_D    Semi-Natural_D ',
     &   ' C_Dec1_A2A  C_Dec1_P2A  C_Dec1_W2A  C_Dec1_N2A ',
     &   ' C_Dec1_M2A  C_Dec1_S2A ',
     &   ' C_Dec1_A2P  C_Dec1_P2P  C_Dec1_W2P  C_Dec1_N2P ',
     &   ' C_Dec1_M2P  C_Dec1_S2P ',
     &   ' C_Dec1_A2W  C_Dec1_P2W  C_Dec1_W2W  C_Dec1_N2W ',
     &   ' C_Dec1_M2W  C_Dec1_S2W ',
     &   ' C_Dec1_A2N  C_Dec1_P2N  C_Dec1_W2N  C_Dec1_N2N ',
     &   ' C_Dec1_M2N  C_Dec1_S2N ',
     &   ' C_Dec1_A2M  C_Dec1_P2M  C_Dec1_W2M  C_Dec1_N2M ',
     &   ' C_Dec1_M2M  C_Dec1_S2M ',
     &   ' C_Dec1_A2S  C_Dec1_P2S  C_Dec1_W2S  C_Dec1_N2S ',
     &   ' C_Dec1_M2S  C_Dec1_S2S ',
     &   ' C_Dec2_A2A  C_Dec2_P2A  C_Dec2_W2A  C_Dec2_N2A ',
     &   ' C_Dec2_M2A  C_Dec2_S2A ',
     &   ' C_Dec2_A2P  C_Dec2_P2P  C_Dec2_W2P  C_Dec2_N2P ',
     &   ' C_Dec2_M2P  C_Dec2_S2P ',
     &   ' C_Dec2_A2W  C_Dec2_P2W  C_Dec2_W2W  C_Dec2_N2W ',
     &   ' C_Dec2_M2W  C_Dec2_S2W ',
     &   ' C_Dec2_A2N  C_Dec2_P2N  C_Dec2_W2N  C_Dec2_N2N ',
     &   ' C_Dec2_M2N  C_Dec2_S2N ',
     &   ' C_Dec2_A2M  C_Dec2_P2M  C_Dec2_W2M  C_Dec2_N2M ',
     &   ' C_Dec2_M2M  C_Dec2_S2M ',
     &   ' C_Dec2_A2S  C_Dec2_P2S  C_Dec2_W2S  C_Dec2_N2S ',
     &   ' C_Dec2_M2S  C_Dec2_S2S ',
     &   ' C_Dec3_A2A  C_Dec3_P2A  C_Dec3_W2A  C_Dec3_N2A ',
     &   ' C_Dec3_M2A  C_Dec3_S2A ',
     &   ' C_Dec3_A2P  C_Dec3_P2P  C_Dec3_W2P  C_Dec3_N2P ',
     &   ' C_Dec3_M2P  C_Dec3_S2P ',
     &   ' C_Dec3_A2W  C_Dec3_P2W  C_Dec3_W2W  C_Dec3_N2W ',
     &   ' C_Dec3_M2W  C_Dec3_S2W ',
     &   ' C_Dec3_A2N  C_Dec3_P2N  C_Dec3_W2N  C_Dec3_N2N ',
     &   ' C_Dec3_M2N  C_Dec3_S2N ',
     &   ' C_Dec3_A2M  C_Dec3_P2M  C_Dec3_W2M  C_Dec3_N2M ',
     &   ' C_Dec3_M2M  C_Dec3_S2M ',
     &   ' C_Dec3_A2S  C_Dec3_P2S  C_Dec3_W2S  C_Dec3_N2S ',
     &   ' C_Dec3_M2S  C_Dec3_S2S ',
     &   ' C_Dec4_A2A  C_Dec4_P2A  C_Dec4_W2A  C_Dec4_N2A ',
     &   ' C_Dec4_M2A  C_Dec4_S2A ',
     &   ' C_Dec4_A2P  C_Dec4_P2P  C_Dec4_W2P  C_Dec4_N2P ',
     &   ' C_Dec4_M2P  C_Dec4_S2P ',
     &   ' C_Dec4_A2W  C_Dec4_P2W  C_Dec4_W2W  C_Dec4_N2W ',
     &   ' C_Dec4_M2W  C_Dec4_S2W ',
     &   ' C_Dec4_A2N  C_Dec4_P2N  C_Dec4_W2N  C_Dec4_N2N ',
     &   ' C_Dec4_M2N  C_Dec4_S2N ',
     &   ' C_Dec4_A2M  C_Dec4_P2M  C_Dec4_W2M  C_Dec4_N2M ',
     &   ' C_Dec4_M2M  C_Dec4_S2M ',
     &   ' C_Dec4_A2S  C_Dec4_P2S  C_Dec4_W2S  C_Dec4_N2S ',
     &   ' C_Dec4_M2S  C_Dec4_S2S ',
     &   ' C_Dec5_A2A  C_Dec5_P2A  C_Dec5_W2A  C_Dec5_N2A ',
     &   ' C_Dec5_M2A  C_Dec5_S2A ',
     &   ' C_Dec5_A2P  C_Dec5_P2P  C_Dec5_W2P  C_Dec5_N2P ',
     &   ' C_Dec5_M2P  C_Dec5_S2P ',
     &   ' C_Dec5_A2W  C_Dec5_P2W  C_Dec5_W2W  C_Dec5_N2W ',
     &   ' C_Dec5_M2W  C_Dec5_S2W ',
     &   ' C_Dec5_A2N  C_Dec5_P2N  C_Dec5_W2N  C_Dec5_N2N ',
     &   ' C_Dec5_M2N  C_Dec5_S2N ',
     &   ' C_Dec5_A2M  C_Dec5_P2M  C_Dec5_W2M  C_Dec5_N2M ',
     &   ' C_Dec5_M2M  C_Dec5_S2M ',
     &   ' C_Dec5_A2S  C_Dec5_P2S  C_Dec5_W2S  C_Dec5_N2S ',
     &   ' C_Dec5_M2S  C_Dec5_S2S ',
     &   ' C_Dec6_A2A  C_Dec6_P2A  C_Dec6_W2A  C_Dec6_N2A ',
     &   ' C_Dec6_M2A  C_Dec6_S2A ',
     &   ' C_Dec6_A2P  C_Dec6_P2P  C_Dec6_W2P  C_Dec6_N2P ',
     &   ' C_Dec6_M2P  C_Dec6_S2P ',
     &   ' C_Dec6_A2W  C_Dec6_P2W  C_Dec6_W2W  C_Dec6_N2W ',
     &   ' C_Dec6_M2W  C_Dec6_S2W ',
     &   ' C_Dec6_A2N  C_Dec6_P2N  C_Dec6_W2N  C_Dec6_N2N ',
     &   ' C_Dec6_M2N  C_Dec6_S2N ',
     &   ' C_Dec6_A2M  C_Dec6_P2M  C_Dec6_W2M  C_Dec6_N2M ',
     &   ' C_Dec6_M2M  C_Dec6_S2M ',
     &   ' C_Dec6_A2S  C_Dec6_P2S  C_Dec6_W2S  C_Dec6_N2S ',
     &   ' C_Dec6_M2S  C_Dec6_S2S ',
     &   ' C_Dec7_A2A  C_Dec7_P2A  C_Dec7_W2A  C_Dec7_N2A ',
     &   ' C_Dec7_M2A  C_Dec7_S2A ',
     &   ' C_Dec7_A2P  C_Dec7_P2P  C_Dec7_W2P  C_Dec7_N2P ',
     &   ' C_Dec7_M2P  C_Dec7_S2P ',
     &   ' C_Dec7_A2W  C_Dec7_P2W  C_Dec7_W2W  C_Dec7_N2W ',
     &   ' C_Dec7_M2W  C_Dec7_S2W ',
     &   ' C_Dec7_A2N  C_Dec7_P2N  C_Dec7_W2N  C_Dec7_N2N ',
     &   ' C_Dec7_M2N  C_Dec7_S2N ',
     &   ' C_Dec7_A2M  C_Dec7_P2M  C_Dec7_W2M  C_Dec7_N2M ',
     &   ' C_Dec7_M2M  C_Dec7_S2M ',
     &   ' C_Dec7_A2S  C_Dec7_P2S  C_Dec7_W2S  C_Dec7_N2S ',
     &   ' C_Dec7_M2S  C_Dec7_S2S ',
     &   ' C_Dec8_A2A  C_Dec8_P2A  C_Dec8_W2A  C_Dec8_N2A ',
     &   ' C_Dec8_M2A  C_Dec8_S2A ',
     &   ' C_Dec8_A2P  C_Dec8_P2P  C_Dec8_W2P  C_Dec8_N2P ',
     &   ' C_Dec8_M2P  C_Dec8_S2P ',
     &   ' C_Dec8_A2W  C_Dec8_P2W  C_Dec8_W2W  C_Dec8_N2W ',
     &   ' C_Dec8_M2W  C_Dec8_S2W ',
     &   ' C_Dec8_A2N  C_Dec8_P2N  C_Dec8_W2N  C_Dec8_N2N ',
     &   ' C_Dec8_M2N  C_Dec8_S2N ',
     &   ' C_Dec8_A2M  C_Dec8_P2M  C_Dec8_W2M  C_Dec8_N2M ',
     &   ' C_Dec8_M2M  C_Dec8_S2M ',
     &   ' C_Dec8_A2S  C_Dec8_P2S  C_Dec8_W2S  C_Dec8_N2S ',
     &   ' C_Dec8_M2S  C_Dec8_S2S ',
     &   ' C_Dec9_A2A  C_Dec9_P2A  C_Dec9_W2A  C_Dec9_N2A ',
     &   ' C_Dec9_M2A  C_Dec9_S2A ',
     &   ' C_Dec9_A2P  C_Dec9_P2P  C_Dec9_W2P  C_Dec9_N2P ',
     &   ' C_Dec9_M2P  C_Dec9_S2P ',
     &   ' C_Dec9_A2W  C_Dec9_P2W  C_Dec9_W2W  C_Dec9_N2W ',
     &   ' C_Dec9_M2W  C_Dec9_S2W ',
     &   ' C_Dec9_A2N  C_Dec9_P2N  C_Dec9_W2N  C_Dec9_N2N ',
     &   ' C_Dec9_M2N  C_Dec9_S2N ',
     &   ' C_Dec9_A2M  C_Dec9_P2M  C_Dec9_W2M  C_Dec9_N2M ',
     &   ' C_Dec9_M2M  C_Dec9_S2M ',
     &   ' C_Dec9_A2S  C_Dec9_P2S  C_Dec9_W2S  C_Dec9_N2S ',
     &   ' C_Dec9_M2S  C_Dec9_S2S ',
     &   '  Soil_C_change       Soil_C_change       Soil_C_change     ',
     &   '  Soil_C_change       Soil_C_change       Soil_C_change     ',
     &   '  Soil_C_change       Soil_C_change       Soil_C_change     ',
     &   'OrgSoil_C_change    OrgSoil_C_change    OrgSoil_C_change    ',
     &   'OrgSoil_C_change    OrgSoil_C_change    OrgSoil_C_change    ',
     &   'OrgSoil_C_change    OrgSoil_C_change    OrgSoil_C_change    ',
     &   'CO2-C_emitted       CO2-C_emitted       CO2-C_emitted       ',
     &   'CO2-C_emitted       CO2-C_emitted       CO2-C_emitted       ',
     &   'CO2-C_emitted       CO2-C_emitted       CO2-C_emitted       ',
     &   'CH4-C_emitted       CH4-C_emitted       CH4-C_emitted       ',
     &   'CH4-C_emitted       CH4-C_emitted       CH4-C_emitted       ',
     &   'CH4-C_emitted       CH4-C_emitted       CH4-C_emitted       ',
     &   'N2O-N_emitted       N2O-N_emitted       N2O-N_emitted       ',
     &   'N2O-N_emitted       N2O-N_emitted       N2O-N_emitted       ',
     &   'N2O-N_emitted       N2O-N_emitted       N2O-N_emitted       ',
     &   'C-equivalents       C-equivalents       C-equivalents       ',
     &   'C-equivalents       C-equivalents       C-equivalents       ',
     &   'C-equivalents       C-equivalents       C-equivalents       ',
     &   'OrgC-equivalents    OrgC-equivalents    OrgC-equivalents    ',
     &   'OrgC-equivalents    OrgC-equivalents    OrgC-equivalents    ',
     &   'OrgC-equivalents    OrgC-equivalents    OrgC-equivalents    ',
     &   'LUCDec1_A2A LUCDec1_P2A LUCDec1_W2A LUCDec1_N2A ',
     &   'LUCDec1_M2A LUCDec1_S2A ',
     &   'LUCDec1_A2P LUCDec1_P2P LUCDec1_W2P LUCDec1_N2P ',
     &   'LUCDec1_M2P LUCDec1_S2P ',
     &   'LUCDec1_A2W LUCDec1_P2W LUCDec1_W2W LUCDec1_N2W ',
     &   'LUCDec1_M2W LUCDec1_S2W ',
     &   'LUCDec1_A2N LUCDec1_P2N LUCDec1_W2N LUCDec1_N2N ',
     &   'LUCDec1_M2N LUCDec1_S2N ',
     &   'LUCDec1_A2M LUCDec1_P2M LUCDec1_W2M LUCDec1_N2M ',
     &   'LUCDec1_M2M LUCDec1_S2M ',
     &   'LUCDec1_A2S LUCDec1_P2S LUCDec1_W2S LUCDec1_N2S ',
     &   'LUCDec1_M2S LUCDec1_S2S ',
     &   'LUCDec2_A2A LUCDec2_P2A LUCDec2_W2A LUCDec2_N2A ',
     &   'LUCDec2_M2A LUCDec2_S2A ',
     &   'LUCDec2_A2P LUCDec2_P2P LUCDec2_W2P LUCDec2_N2P ',
     &   'LUCDec2_M2P LUCDec2_S2P ',
     &   'LUCDec2_A2W LUCDec2_P2W LUCDec2_W2W LUCDec2_N2W ',
     &   'LUCDec2_M2W LUCDec2_S2W ',
     &   'LUCDec2_A2N LUCDec2_P2N LUCDec2_W2N LUCDec2_N2N ',
     &   'LUCDec2_M2N LUCDec2_S2N ',
     &   'LUCDec2_A2M LUCDec2_P2M LUCDec2_W2M LUCDec2_N2M ',
     &   'LUCDec2_M2M LUCDec2_S2M ',
     &   'LUCDec2_A2S LUCDec2_P2S LUCDec2_W2S LUCDec2_N2S ',
     &   'LUCDec2_M2S LUCDec2_S2S ',
     &   'LUCDec3_A2A LUCDec3_P2A LUCDec3_W2A LUCDec3_N2A ',
     &   'LUCDec3_M2A LUCDec3_S2A ',
     &   'LUCDec3_A2P LUCDec3_P2P LUCDec3_W2P LUCDec3_N2P ',
     &   'LUCDec3_M2P LUCDec3_S2P ',
     &   'LUCDec3_A2W LUCDec3_P2W LUCDec3_W2W LUCDec3_N2W ',
     &   'LUCDec3_M2W LUCDec3_S2W ',
     &   'LUCDec3_A2N LUCDec3_P2N LUCDec3_W2N LUCDec3_N2N ',
     &   'LUCDec3_M2N LUCDec3_S2N ',
     &   'LUCDec3_A2M LUCDec3_P2M LUCDec3_W2M LUCDec3_N2M ',
     &   'LUCDec3_M2M LUCDec3_S2M ',
     &   'LUCDec3_A2S LUCDec3_P2S LUCDec3_W2S LUCDec3_N2S ',
     &   'LUCDec3_M2S LUCDec3_S2S ',
     &   'LUCDec4_A2A LUCDec4_P2A LUCDec4_W2A LUCDec4_N2A ',
     &   'LUCDec4_M2A LUCDec4_S2A ',
     &   'LUCDec4_A2P LUCDec4_P2P LUCDec4_W2P LUCDec4_N2P ',
     &   'LUCDec4_M2P LUCDec4_S2P ',
     &   'LUCDec4_A2W LUCDec4_P2W LUCDec4_W2W LUCDec4_N2W ',
     &   'LUCDec4_M2W LUCDec4_S2W ',
     &   'LUCDec4_A2N LUCDec4_P2N LUCDec4_W2N LUCDec4_N2N ',
     &   'LUCDec4_M2N LUCDec4_S2N ',
     &   'LUCDec4_A2M LUCDec4_P2M LUCDec4_W2M LUCDec4_N2M ',
     &   'LUCDec4_M2M LUCDec4_S2M ',
     &   'LUCDec4_A2S LUCDec4_P2S LUCDec4_W2S LUCDec4_N2S ',
     &   'LUCDec4_M2S LUCDec4_S2S ',
     &   'LUCDec5_A2A LUCDec5_P2A LUCDec5_W2A LUCDec5_N2A ',
     &   'LUCDec5_M2A LUCDec5_S2A ',
     &   'LUCDec5_A2P LUCDec5_P2P LUCDec5_W2P LUCDec5_N2P ',
     &   'LUCDec5_M2P LUCDec5_S2P ',
     &   'LUCDec5_A2W LUCDec5_P2W LUCDec5_W2W LUCDec5_N2W ',
     &   'LUCDec5_M2W LUCDec5_S2W ',
     &   'LUCDec5_A2N LUCDec5_P2N LUCDec5_W2N LUCDec5_N2N ',
     &   'LUCDec5_M2N LUCDec5_S2N ',
     &   'LUCDec5_A2M LUCDec5_P2M LUCDec5_W2M LUCDec5_N2M ',
     &   'LUCDec5_M2M LUCDec5_S2M ',
     &   'LUCDec5_A2S LUCDec5_P2S LUCDec5_W2S LUCDec5_N2S ',
     &   'LUCDec5_M2S LUCDec5_S2S ',
     &   'LUCDec6_A2A LUCDec6_P2A LUCDec6_W2A LUCDec6_N2A ',
     &   'LUCDec6_M2A LUCDec6_S2A ',
     &   'LUCDec6_A2P LUCDec6_P2P LUCDec6_W2P LUCDec6_N2P ',
     &   'LUCDec6_M2P LUCDec6_S2P ',
     &   'LUCDec6_A2W LUCDec6_P2W LUCDec6_W2W LUCDec6_N2W ',
     &   'LUCDec6_M2W LUCDec6_S2W ',
     &   'LUCDec6_A2N LUCDec6_P2N LUCDec6_W2N LUCDec6_N2N ',
     &   'LUCDec6_M2N LUCDec6_S2N ',
     &   'LUCDec6_A2M LUCDec6_P2M LUCDec6_W2M LUCDec6_N2M ',
     &   'LUCDec6_M2M LUCDec6_S2M ',
     &   'LUCDec6_A2S LUCDec6_P2S LUCDec6_W2S LUCDec6_N2S ',
     &   'LUCDec6_M2S LUCDec6_S2S ',
     &   'LUCDec7_A2A LUCDec7_P2A LUCDec7_W2A LUCDec7_N2A ',
     &   'LUCDec7_M2A LUCDec7_S2A ',
     &   'LUCDec7_A2P LUCDec7_P2P LUCDec7_W2P LUCDec7_N2P ',
     &   'LUCDec7_M2P LUCDec7_S2P ',
     &   'LUCDec7_A2W LUCDec7_P2W LUCDec7_W2W LUCDec7_N2W ',
     &   'LUCDec7_M2W LUCDec7_S2W ',
     &   'LUCDec7_A2N LUCDec7_P2N LUCDec7_W2N LUCDec7_N2N ',
     &   'LUCDec7_M2N LUCDec7_S2N ',
     &   'LUCDec7_A2M LUCDec7_P2M LUCDec7_W2M LUCDec7_N2M ',
     &   'LUCDec7_M2M LUCDec7_S2M ',
     &   'LUCDec7_A2S LUCDec7_P2S LUCDec7_W2S LUCDec7_N2S ',
     &   'LUCDec7_M2S LUCDec7_S2S ',
     &   'LUCDec8_A2A LUCDec8_P2A LUCDec8_W2A LUCDec8_N2A ',
     &   'LUCDec8_M2A LUCDec8_S2A ',
     &   'LUCDec8_A2P LUCDec8_P2P LUCDec8_W2P LUCDec8_N2P ',
     &   'LUCDec8_M2P LUCDec8_S2P ',
     &   'LUCDec8_A2W LUCDec8_P2W LUCDec8_W2W LUCDec8_N2W ',
     &   'LUCDec8_M2W LUCDec8_S2W ',
     &   'LUCDec8_A2N LUCDec8_P2N LUCDec8_W2N LUCDec8_N2N ',
     &   'LUCDec8_M2N LUCDec8_S2N ',
     &   'LUCDec8_A2M LUCDec8_P2M LUCDec8_W2M LUCDec8_N2M ',
     &   'LUCDec8_M2M LUCDec8_S2M ',
     &   'LUCDec8_A2S LUCDec8_P2S LUCDec8_W2S LUCDec8_N2S ',
     &   'LUCDec8_M2S LUCDec8_S2S ',
     &   'LUCDec9_A2A LUCDec9_P2A LUCDec9_W2A LUCDec9_N2A ',
     &   'LUCDec9_M2A LUCDec9_S2A ',
     &   'LUCDec9_A2P LUCDec9_P2P LUCDec9_W2P LUCDec9_N2P ',
     &   'LUCDec9_M2P LUCDec9_S2P ',
     &   'LUCDec9_A2W LUCDec9_P2W LUCDec9_W2W LUCDec9_N2W ',
     &   'LUCDec9_M2W LUCDec9_S2W ',
     &   'LUCDec9_A2N LUCDec9_P2N LUCDec9_W2N LUCDec9_N2N ',
     &   'LUCDec9_M2N LUCDec9_S2N ',
     &   'LUCDec9_A2M LUCDec9_P2M LUCDec9_W2M LUCDec9_N2M ',
     &   'LUCDec9_M2M LUCDec9_S2M ',
     &   'LUCDec9_A2S LUCDec9_P2S LUCDec9_W2S LUCDec9_N2S ',
     &   'LUCDec9_M2S LUCDec9_S2S '/
     &   '20km 1km E N depth                                ',
     &   'kg/ha kg/ha kg/ha kg/ha ',
     &   'cm cm cm cm ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' kt/km2/10y  kt/km2/10y  kt/km2/10y  kt/km2/10y   ',
     &   ' kt/km2/10y  kt/km2/10y ',
     &   ' Dec1(kt/km2/10yr)   Dec2(kt/km2/10yr)   Dec3(kt/km2/10yr)  ',
     &   ' Dec4(kt/km2/10yr)   Dec5(kt/km2/10yr)   Dec6(kt/km2/10yr)  ',
     &   ' Dec7(kt/km2/10yr)   Dec8(kt/km2/10yr)   Dec9(kt/km2/10yr)  ',
     &   ' Dec1(kt/km2/10yr)   Dec2(kt/km2/10yr)   Dec3(kt/km2/10yr)  ',
     &   ' Dec4(kt/km2/10yr)   Dec5(kt/km2/10yr)   Dec6(kt/km2/10yr)  ',
     &   ' Dec7(kt/km2/10yr)   Dec8(kt/km2/10yr)   Dec9(kt/km2/10yr)  ',
     &   ' Dec1(kt/km2/10yr)   Dec2(kt/km2/10yr)   Dec3(kt/km2/10yr)  ',
     &   ' Dec4(kt/km2/10yr)   Dec5(kt/km2/10yr)   Dec6(kt/km2/10yr)  ',
     &   ' Dec7(kt/km2/10yr)   Dec8(kt/km2/10yr)   Dec9(kt/km2/10yr)  ',
     &   ' Dec1(kt/km2/10yr)   Dec2(kt/km2/10yr)   Dec3(kt/km2/10yr)  ',
     &   ' Dec4(kt/km2/10yr)   Dec5(kt/km2/10yr)   Dec6(kt/km2/10yr)  ',
     &   ' Dec7(kt/km2/10yr)   Dec8(kt/km2/10yr)   Dec9(kt/km2/10yr)  ',
     &   ' Dec1(kt/km2/10yr)   Dec2(kt/km2/10yr)   Dec3(kt/km2/10yr)  ',
     &   ' Dec4(kt/km2/10yr)   Dec5(kt/km2/10yr)   Dec6(kt/km2/10yr)  ',
     &   ' Dec7(kt/km2/10yr)   Dec8(kt/km2/10yr)   Dec9(kt/km2/10yr)  ',
     &   ' Dec1(kt/km2/10yr)   Dec2(kt/km2/10yr)   Dec3(kt/km2/10yr)  ',
     &   ' Dec4(kt/km2/10yr)   Dec5(kt/km2/10yr)   Dec6(kt/km2/10yr)  ',
     &   ' Dec7(kt/km2/10yr)   Dec8(kt/km2/10yr)   Dec9(kt/km2/10yr)  ',
     &   ' Dec1(kt/km2/10yr)   Dec2(kt/km2/10yr)   Dec3(kt/km2/10yr)  ',
     &   ' Dec4(kt/km2/10yr)   Dec5(kt/km2/10yr)   Dec6(kt/km2/10yr)  ',
     &   ' Dec7(kt/km2/10yr)   Dec8(kt/km2/10yr)   Dec9(kt/km2/10yr)  ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ',
     &   'km2/10y/km2 km2/10y/km2 km2/10y/km2 km2/10y/km2 ')
C
C Title lines for grid file
C   
      WRITE(GRIDCHAN,450)
      WRITE(OGRIDCHAN,450)
450   FORMAT('All results in t/20km2/10y'/
     &       '20KM2_GRID    EAST    NORTH ',
     &   ' C_Dec1_A2A  C_Dec1_P2A  C_Dec1_W2A  C_Dec1_N2A ',
     &   ' C_Dec1_M2A  C_Dec1_S2A ',
     &   ' C_Dec1_A2P  C_Dec1_P2P  C_Dec1_W2P  C_Dec1_N2P ',
     &   ' C_Dec1_M2P  C_Dec1_S2P ',
     &   ' C_Dec1_A2W  C_Dec1_P2W  C_Dec1_W2W  C_Dec1_N2W ',
     &   ' C_Dec1_M2W  C_Dec1_S2W ',
     &   ' C_Dec1_A2N  C_Dec1_P2N  C_Dec1_W2N  C_Dec1_N2N ',
     &   ' C_Dec1_M2N  C_Dec1_S2N ',
     &   ' C_Dec1_A2M  C_Dec1_P2M  C_Dec1_W2M  C_Dec1_N2M ',
     &   ' C_Dec1_M2M  C_Dec1_S2M ',
     &   ' C_Dec1_A2S  C_Dec1_P2S  C_Dec1_W2S  C_Dec1_N2S ',
     &   ' C_Dec1_M2S  C_Dec1_S2S ',
     &   ' C_Dec2_A2A  C_Dec2_P2A  C_Dec2_W2A  C_Dec2_N2A ',
     &   ' C_Dec2_M2A  C_Dec2_S2A ',
     &   ' C_Dec2_A2P  C_Dec2_P2P  C_Dec2_W2P  C_Dec2_N2P ',
     &   ' C_Dec2_M2P  C_Dec2_S2P ',
     &   ' C_Dec2_A2W  C_Dec2_P2W  C_Dec2_W2W  C_Dec2_N2W ',
     &   ' C_Dec2_M2W  C_Dec2_S2W ',
     &   ' C_Dec2_A2N  C_Dec2_P2N  C_Dec2_W2N  C_Dec2_N2N ',
     &   ' C_Dec2_M2N  C_Dec2_S2N ',
     &   ' C_Dec2_A2M  C_Dec2_P2M  C_Dec2_W2M  C_Dec2_N2M ',
     &   ' C_Dec2_M2M  C_Dec2_S2M ',
     &   ' C_Dec2_A2S  C_Dec2_P2S  C_Dec2_W2S  C_Dec2_N2S ',
     &   ' C_Dec2_M2S  C_Dec2_S2S ',
     &   ' C_Dec3_A2A  C_Dec3_P2A  C_Dec3_W2A  C_Dec3_N2A ',
     &   ' C_Dec3_M2A  C_Dec3_S2A ',
     &   ' C_Dec3_A2P  C_Dec3_P2P  C_Dec3_W2P  C_Dec3_N2P ',
     &   ' C_Dec3_M2P  C_Dec3_S2P ',
     &   ' C_Dec3_A2W  C_Dec3_P2W  C_Dec3_W2W  C_Dec3_N2W ',
     &   ' C_Dec3_M2W  C_Dec3_S2W ',
     &   ' C_Dec3_A2N  C_Dec3_P2N  C_Dec3_W2N  C_Dec3_N2N ',
     &   ' C_Dec3_M2N  C_Dec3_S2N ',
     &   ' C_Dec3_A2M  C_Dec3_P2M  C_Dec3_W2M  C_Dec3_N2M ',
     &   ' C_Dec3_M2M  C_Dec3_S2M ',
     &   ' C_Dec3_A2S  C_Dec3_P2S  C_Dec3_W2S  C_Dec3_N2S ',
     &   ' C_Dec3_M2S  C_Dec3_S2S ',
     &   ' C_Dec4_A2A  C_Dec4_P2A  C_Dec4_W2A  C_Dec4_N2A ',
     &   ' C_Dec4_M2A  C_Dec4_S2A ',
     &   ' C_Dec4_A2P  C_Dec4_P2P  C_Dec4_W2P  C_Dec4_N2P ',
     &   ' C_Dec4_M2P  C_Dec4_S2P ',
     &   ' C_Dec4_A2W  C_Dec4_P2W  C_Dec4_W2W  C_Dec4_N2W ',
     &   ' C_Dec4_M2W  C_Dec4_S2W ',
     &   ' C_Dec4_A2N  C_Dec4_P2N  C_Dec4_W2N  C_Dec4_N2N ',
     &   ' C_Dec4_M2N  C_Dec4_S2N ',
     &   ' C_Dec4_A2M  C_Dec4_P2M  C_Dec4_W2M  C_Dec4_N2M ',
     &   ' C_Dec4_M2M  C_Dec4_S2M ',
     &   ' C_Dec4_A2S  C_Dec4_P2S  C_Dec4_W2S  C_Dec4_N2S ',
     &   ' C_Dec4_M2S  C_Dec4_S2S ',
     &   ' C_Dec5_A2A  C_Dec5_P2A  C_Dec5_W2A  C_Dec5_N2A ',
     &   ' C_Dec5_M2A  C_Dec5_S2A ',
     &   ' C_Dec5_A2P  C_Dec5_P2P  C_Dec5_W2P  C_Dec5_N2P ',
     &   ' C_Dec5_M2P  C_Dec5_S2P ',
     &   ' C_Dec5_A2W  C_Dec5_P2W  C_Dec5_W2W  C_Dec5_N2W ',
     &   ' C_Dec5_M2W  C_Dec5_S2W ',
     &   ' C_Dec5_A2N  C_Dec5_P2N  C_Dec5_W2N  C_Dec5_N2N ',
     &   ' C_Dec5_M2N  C_Dec5_S2N ',
     &   ' C_Dec5_A2M  C_Dec5_P2M  C_Dec5_W2M  C_Dec5_N2M ',
     &   ' C_Dec5_M2M  C_Dec5_S2M ',
     &   ' C_Dec5_A2S  C_Dec5_P2S  C_Dec5_W2S  C_Dec5_N2S ',
     &   ' C_Dec5_M2S  C_Dec5_S2S ',
     &   ' C_Dec6_A2A  C_Dec6_P2A  C_Dec6_W2A  C_Dec6_N2A ',
     &   ' C_Dec6_M2A  C_Dec6_S2A ',
     &   ' C_Dec6_A2P  C_Dec6_P2P  C_Dec6_W2P  C_Dec6_N2P ',
     &   ' C_Dec6_M2P  C_Dec6_S2P ',
     &   ' C_Dec6_A2W  C_Dec6_P2W  C_Dec6_W2W  C_Dec6_N2W ',
     &   ' C_Dec6_M2W  C_Dec6_S2W ',
     &   ' C_Dec6_A2N  C_Dec6_P2N  C_Dec6_W2N  C_Dec6_N2N ',
     &   ' C_Dec6_M2N  C_Dec6_S2N ',
     &   ' C_Dec6_A2M  C_Dec6_P2M  C_Dec6_W2M  C_Dec6_N2M ',
     &   ' C_Dec6_M2M  C_Dec6_S2M ',
     &   ' C_Dec6_A2S  C_Dec6_P2S  C_Dec6_W2S  C_Dec6_N2S ',
     &   ' C_Dec6_M2S  C_Dec6_S2S ',
     &   ' C_Dec7_A2A  C_Dec7_P2A  C_Dec7_W2A  C_Dec7_N2A ',
     &   ' C_Dec7_M2A  C_Dec7_S2A ',
     &   ' C_Dec7_A2P  C_Dec7_P2P  C_Dec7_W2P  C_Dec7_N2P ',
     &   ' C_Dec7_M2P  C_Dec7_S2P ',
     &   ' C_Dec7_A2W  C_Dec7_P2W  C_Dec7_W2W  C_Dec7_N2W ',
     &   ' C_Dec7_M2W  C_Dec7_S2W ',
     &   ' C_Dec7_A2N  C_Dec7_P2N  C_Dec7_W2N  C_Dec7_N2N ',
     &   ' C_Dec7_M2N  C_Dec7_S2N ',
     &   ' C_Dec7_A2M  C_Dec7_P2M  C_Dec7_W2M  C_Dec7_N2M ',
     &   ' C_Dec7_M2M  C_Dec7_S2M ',
     &   ' C_Dec7_A2S  C_Dec7_P2S  C_Dec7_W2S  C_Dec7_N2S ',
     &   ' C_Dec7_M2S  C_Dec7_S2S ',
     &   ' C_Dec8_A2A  C_Dec8_P2A  C_Dec8_W2A  C_Dec8_N2A ',
     &   ' C_Dec8_M2A  C_Dec8_S2A ',
     &   ' C_Dec8_A2P  C_Dec8_P2P  C_Dec8_W2P  C_Dec8_N2P ',
     &   ' C_Dec8_M2P  C_Dec8_S2P ',
     &   ' C_Dec8_A2W  C_Dec8_P2W  C_Dec8_W2W  C_Dec8_N2W ',
     &   ' C_Dec8_M2W  C_Dec8_S2W ',
     &   ' C_Dec8_A2N  C_Dec8_P2N  C_Dec8_W2N  C_Dec8_N2N ',
     &   ' C_Dec8_M2N  C_Dec8_S2N ',
     &   ' C_Dec8_A2M  C_Dec8_P2M  C_Dec8_W2M  C_Dec8_N2M ',
     &   ' C_Dec8_M2M  C_Dec8_S2M ',
     &   ' C_Dec8_A2S  C_Dec8_P2S  C_Dec8_W2S  C_Dec8_N2S ',
     &   ' C_Dec8_M2S  C_Dec8_S2S ',
     &   ' C_Dec9_A2A  C_Dec9_P2A  C_Dec9_W2A  C_Dec9_N2A ',
     &   ' C_Dec9_M2A  C_Dec9_S2A ',
     &   ' C_Dec9_A2P  C_Dec9_P2P  C_Dec9_W2P  C_Dec9_N2P ',
     &   ' C_Dec9_M2P  C_Dec9_S2P ',
     &   ' C_Dec9_A2W  C_Dec9_P2W  C_Dec9_W2W  C_Dec9_N2W ',
     &   ' C_Dec9_M2W  C_Dec9_S2W ',
     &   ' C_Dec9_A2N  C_Dec9_P2N  C_Dec9_W2N  C_Dec9_N2N ',
     &   ' C_Dec9_M2N  C_Dec9_S2N ',
     &   ' C_Dec9_A2M  C_Dec9_P2M  C_Dec9_W2M  C_Dec9_N2M ',
     &   ' C_Dec9_M2M  C_Dec9_S2M ',
     &   ' C_Dec9_A2S  C_Dec9_P2S  C_Dec9_W2S  C_Dec9_N2S ',
     &   ' C_Dec9_M2S  C_Dec9_S2S       ',
     &   'LUCDec1_A2A LUCDec1_P2A LUCDec1_W2A LUCDec1_N2A ',
     &   'LUCDec1_M2A LUCDec1_S2A ',
     &   'LUCDec1_A2P LUCDec1_P2P LUCDec1_W2P LUCDec1_N2P ',
     &   'LUCDec1_M2P LUCDec1_S2P ',
     &   'LUCDec1_A2W LUCDec1_P2W LUCDec1_W2W LUCDec1_N2W ',
     &   'LUCDec1_M2W LUCDec1_S2W ',
     &   'LUCDec1_A2N LUCDec1_P2N LUCDec1_W2N LUCDec1_N2N ',
     &   'LUCDec1_M2N LUCDec1_S2N ',
     &   'LUCDec1_A2M LUCDec1_P2M LUCDec1_W2M LUCDec1_N2M ',
     &   'LUCDec1_M2M LUCDec1_S2M ',
     &   'LUCDec1_A2S LUCDec1_P2S LUCDec1_W2S LUCDec1_N2S ',
     &   'LUCDec1_M2S LUCDec1_S2S ',
     &   'LUCDec2_A2A LUCDec2_P2A LUCDec2_W2A LUCDec2_N2A ',
     &   'LUCDec2_M2A LUCDec2_S2A ',
     &   'LUCDec2_A2P LUCDec2_P2P LUCDec2_W2P LUCDec2_N2P ',
     &   'LUCDec2_M2P LUCDec2_S2P ',
     &   'LUCDec2_A2W LUCDec2_P2W LUCDec2_W2W LUCDec2_N2W ',
     &   'LUCDec2_M2W LUCDec2_S2W ',
     &   'LUCDec2_A2N LUCDec2_P2N LUCDec2_W2N LUCDec2_N2N ',
     &   'LUCDec2_M2N LUCDec2_S2N ',
     &   'LUCDec2_A2M LUCDec2_P2M LUCDec2_W2M LUCDec2_N2M ',
     &   'LUCDec2_M2M LUCDec2_S2M ',
     &   'LUCDec2_A2S LUCDec2_P2S LUCDec2_W2S LUCDec2_N2S ',
     &   'LUCDec2_M2S LUCDec2_S2S ',
     &   'LUCDec3_A2A LUCDec3_P2A LUCDec3_W2A LUCDec3_N2A ',
     &   'LUCDec3_M2A LUCDec3_S2A ',
     &   'LUCDec3_A2P LUCDec3_P2P LUCDec3_W2P LUCDec3_N2P ',
     &   'LUCDec3_M2P LUCDec3_S2P ',
     &   'LUCDec3_A2W LUCDec3_P2W LUCDec3_W2W LUCDec3_N2W ',
     &   'LUCDec3_M2W LUCDec3_S2W ',
     &   'LUCDec3_A2N LUCDec3_P2N LUCDec3_W2N LUCDec3_N2N ',
     &   'LUCDec3_M2N LUCDec3_S2N ',
     &   'LUCDec3_A2M LUCDec3_P2M LUCDec3_W2M LUCDec3_N2M ',
     &   'LUCDec3_M2M LUCDec3_S2M ',
     &   'LUCDec3_A2S LUCDec3_P2S LUCDec3_W2S LUCDec3_N2S ',
     &   'LUCDec3_M2S LUCDec3_S2S ',
     &   'LUCDec4_A2A LUCDec4_P2A LUCDec4_W2A LUCDec4_N2A ',
     &   'LUCDec4_M2A LUCDec4_S2A ',
     &   'LUCDec4_A2P LUCDec4_P2P LUCDec4_W2P LUCDec4_N2P ',
     &   'LUCDec4_M2P LUCDec4_S2P ',
     &   'LUCDec4_A2W LUCDec4_P2W LUCDec4_W2W LUCDec4_N2W ',
     &   'LUCDec4_M2W LUCDec4_S2W ',
     &   'LUCDec4_A2N LUCDec4_P2N LUCDec4_W2N LUCDec4_N2N ',
     &   'LUCDec4_M2N LUCDec4_S2N ',
     &   'LUCDec4_A2M LUCDec4_P2M LUCDec4_W2M LUCDec4_N2M ',
     &   'LUCDec4_M2M LUCDec4_S2M ',
     &   'LUCDec4_A2S LUCDec4_P2S LUCDec4_W2S LUCDec4_N2S ',
     &   'LUCDec4_M2S LUCDec4_S2S ',
     &   'LUCDec5_A2A LUCDec5_P2A LUCDec5_W2A LUCDec5_N2A ',
     &   'LUCDec5_M2A LUCDec5_S2A ',
     &   'LUCDec5_A2P LUCDec5_P2P LUCDec5_W2P LUCDec5_N2P ',
     &   'LUCDec5_M2P LUCDec5_S2P ',
     &   'LUCDec5_A2W LUCDec5_P2W LUCDec5_W2W LUCDec5_N2W ',
     &   'LUCDec5_M2W LUCDec5_S2W ',
     &   'LUCDec5_A2N LUCDec5_P2N LUCDec5_W2N LUCDec5_N2N ',
     &   'LUCDec5_M2N LUCDec5_S2N ',
     &   'LUCDec5_A2M LUCDec5_P2M LUCDec5_W2M LUCDec5_N2M ',
     &   'LUCDec5_M2M LUCDec5_S2M ',
     &   'LUCDec5_A2S LUCDec5_P2S LUCDec5_W2S LUCDec5_N2S ',
     &   'LUCDec5_M2S LUCDec5_S2S ',
     &   'LUCDec6_A2A LUCDec6_P2A LUCDec6_W2A LUCDec6_N2A ',
     &   'LUCDec6_M2A LUCDec6_S2A ',
     &   'LUCDec6_A2P LUCDec6_P2P LUCDec6_W2P LUCDec6_N2P ',
     &   'LUCDec6_M2P LUCDec6_S2P ',
     &   'LUCDec6_A2W LUCDec6_P2W LUCDec6_W2W LUCDec6_N2W ',
     &   'LUCDec6_M2W LUCDec6_S2W ',
     &   'LUCDec6_A2N LUCDec6_P2N LUCDec6_W2N LUCDec6_N2N ',
     &   'LUCDec6_M2N LUCDec6_S2N ',
     &   'LUCDec6_A2M LUCDec6_P2M LUCDec6_W2M LUCDec6_N2M ',
     &   'LUCDec6_M2M LUCDec6_S2M ',
     &   'LUCDec6_A2S LUCDec6_P2S LUCDec6_W2S LUCDec6_N2S ',
     &   'LUCDec6_M2S LUCDec6_S2S ',
     &   'LUCDec7_A2A LUCDec7_P2A LUCDec7_W2A LUCDec7_N2A ',
     &   'LUCDec7_M2A LUCDec7_S2A ',
     &   'LUCDec7_A2P LUCDec7_P2P LUCDec7_W2P LUCDec7_N2P ',
     &   'LUCDec7_M2P LUCDec7_S2P ',
     &   'LUCDec7_A2W LUCDec7_P2W LUCDec7_W2W LUCDec7_N2W ',
     &   'LUCDec7_M2W LUCDec7_S2W ',
     &   'LUCDec7_A2N LUCDec7_P2N LUCDec7_W2N LUCDec7_N2N ',
     &   'LUCDec7_M2N LUCDec7_S2N ',
     &   'LUCDec7_A2M LUCDec7_P2M LUCDec7_W2M LUCDec7_N2M ',
     &   'LUCDec7_M2M LUCDec7_S2M ',
     &   'LUCDec7_A2S LUCDec7_P2S LUCDec7_W2S LUCDec7_N2S ',
     &   'LUCDec7_M2S LUCDec7_S2S ',
     &   'LUCDec8_A2A LUCDec8_P2A LUCDec8_W2A LUCDec8_N2A ',
     &   'LUCDec8_M2A LUCDec8_S2A ',
     &   'LUCDec8_A2P LUCDec8_P2P LUCDec8_W2P LUCDec8_N2P ',
     &   'LUCDec8_M2P LUCDec8_S2P ',
     &   'LUCDec8_A2W LUCDec8_P2W LUCDec8_W2W LUCDec8_N2W ',
     &   'LUCDec8_M2W LUCDec8_S2W ',
     &   'LUCDec8_A2N LUCDec8_P2N LUCDec8_W2N LUCDec8_N2N ',
     &   'LUCDec8_M2N LUCDec8_S2N ',
     &   'LUCDec8_A2M LUCDec8_P2M LUCDec8_W2M LUCDec8_N2M ',
     &   'LUCDec8_M2M LUCDec8_S2M ',
     &   'LUCDec8_A2S LUCDec8_P2S LUCDec8_W2S LUCDec8_N2S ',
     &   'LUCDec8_M2S LUCDec8_S2S ',
     &   'LUCDec9_A2A LUCDec9_P2A LUCDec9_W2A LUCDec9_N2A ',
     &   'LUCDec9_M2A LUCDec9_S2A ',
     &   'LUCDec9_A2P LUCDec9_P2P LUCDec9_W2P LUCDec9_N2P ',
     &   'LUCDec9_M2P LUCDec9_S2P ',
     &   'LUCDec9_A2W LUCDec9_P2W LUCDec9_W2W LUCDec9_N2W ',
     &   'LUCDec9_M2W LUCDec9_S2W ',
     &   'LUCDec9_A2N LUCDec9_P2N LUCDec9_W2N LUCDec9_N2N ',
     &   'LUCDec9_M2N LUCDec9_S2N ',
     &   'LUCDec9_A2M LUCDec9_P2M LUCDec9_W2M LUCDec9_N2M ',
     &   'LUCDec9_M2M LUCDec9_S2M ',
     &   'LUCDec9_A2S LUCDec9_P2S LUCDec9_W2S LUCDec9_N2S ',
     &   'LUCDec9_M2S LUCDec9_S2S ',
     &     'CChge_Dec_1 CChge_Dec_2 CChge_Dec_3 CChge_Dec_4 ',
     &     'CChge_Dec_5 CChge_Dec_6 CChge_Dec_7 CChge_Dec_8 ',
     &     'CChge_Dec_9 ',
     &     'OCCh_Dec_1  OCCh_Dec_2  OCCh_Dec_3  OCCh_Dec_4  ',
     &     'OCCh_Dec_5  OCCh_Dec_6  OCCh_Dec_7  OCCh_Dec_8  ',
     &     'OCCh_Dec_9  ',
     &     ' C-eq_Dec_1  C-eq_Dec_2  C-eq_Dec_3  C-eq_Dec_4 ',
     &     ' C-eq_Dec_5  C-eq_Dec_6  C-eq_Dec_7  C-eq_Dec_8 ',
     &     ' C-eq_Dec_9 ',
     &     'OC-eq_Dec_1 OC-eq_Dec_2 OC-eq_Dec_3 OC-eq_Dec_4 ',
     &     'OC-eq_Dec_5 OC-eq_Dec_6 OC-eq_Dec_7 OC-eq_Dec_8 ',
     &     'OC-eq_Dec_9 ')
C
C Title lines for climate change output file
C   
      WRITE(FRESCHAN,460)
      WRITE(OFRESCHAN,460)
460   FORMAT('20KM2_GRID      SQUID     EAST    NORTH ',
     &     '  Dec_1_ARA   Dec_1_GRA   Dec_1_FOR   Dec_1_NAT ',
     &     '  Dec_1_MIS   Dec_1_SRC ',
     &     '  Dec_2_ARA   Dec_2_GRA   Dec_2_FOR   Dec_2_NAT ',
     &     '  Dec_2_MIS   Dec_2_SRC ',
     &     '  Dec_3_ARA   Dec_3_GRA   Dec_3_FOR   Dec_3_NAT ',
     &     '  Dec_3_MIS   Dec_3_SRC ',
     &     '  Dec_4_ARA   Dec_4_GRA   Dec_4_FOR   Dec_4_NAT ',
     &     '  Dec_4_MIS   Dec_4_SRC ',
     &     '  Dec_5_ARA   Dec_5_GRA   Dec_5_FOR   Dec_5_NAT ',
     &     '  Dec_5_MIS   Dec_5_SRC ',
     &     '  Dec_6_ARA   Dec_6_GRA   Dec_6_FOR   Dec_6_NAT ',
     &     '  Dec_6_MIS   Dec_6_SRC ',
     &     '  Dec_7_ARA   Dec_7_GRA   Dec_7_FOR   Dec_7_NAT ',
     &     '  Dec_7_MIS   Dec_7_SRC ',
     &     '  Dec_8_ARA   Dec_8_GRA   Dec_8_FOR   Dec_8_NAT ',
     &     '  Dec_8_MIS   Dec_8_SRC ',
     &     '  Dec_9_ARA   Dec_9_GRA   Dec_9_FOR   Dec_9_NAT ',
     &     '  Dec_9_MIS   Dec_9_SRC ')
C
C Title lines for mitigation output file
C   
	DO 700 ILU1=1,MAXLU
	  DO 800 ILU2=1,MAXLU
	    IF(ILU1.NE.ILU2)THEN
	      ICHAN=((ILU1-1)*MAXLU)+ILU2
	      WRITE(MITICHAN(ICHAN),810)
	    ENDIF
800     CONTINUE
700   CONTINUE
810   FORMAT('20KM2_GRID      SQUID     EAST    NORTH ',
     &     '   Series1    WetClass          %C       %Clay   ',
     &     '    LUChange  Dec1_CChange  Dec2_CChange  Dec3_CChange  ',
     &     'Dec4_CChange  Dec5_CChange  Dec6_CChange  Dec7_CChange  ',
     &     'Dec8_CChange  Dec9_CChange  ',
     &     '   Series2    WetClass          %C       %Clay  ',
     &     '    LUChange  Dec1_CChange  Dec2_CChange  Dec3_CChange  ',
     &     'Dec4_CChange  Dec5_CChange  Dec6_CChange  Dec7_CChange  ',
     &     'Dec8_CChange  Dec9_CChange  ',
     &     '   Series3    WetClass          %C       %Clay  ',
     &     '    LUChange  Dec1_CChange  Dec2_CChange  Dec3_CChange  ',
     &     'Dec4_CChange  Dec5_CChange  Dec6_CChange  Dec7_CChange  ',
     &     'Dec8_CChange  Dec9_CChange  ',
     &     '   Series4    WetClass          %C       %Clay  ',
     &     '    LUChange  Dec1_CChange  Dec2_CChange  Dec3_CChange  ',
     &     'Dec4_CChange  Dec5_CChange  Dec6_CChange  Dec7_CChange  ',
     &     'Dec8_CChange  Dec9_CChange  ',
     &     '   Series5    WetClass          %C       %Clay  ',
     &     '    LUChange  Dec1_CChange  Dec2_CChange  Dec3_CChange  ',
     &     'Dec4_CChange  Dec5_CChange  Dec6_CChange  Dec7_CChange  ',
     &     'Dec8_CChange  Dec9_CChange  ',
     &      /'                                       ',
     &     '                                                   ',
     &     '   (km2/km2)    (t/ha/dec)    (t/ha/dec)    (t/ha/dec)  ',
     &     '  (t/ha/dec)    (t/ha/dec)    (t/ha/dec)    (t/ha/dec)  ',
     &     '  (t/ha/dec)    (t/ha/dec)  ',
     &     '                                                   ',
     &     '   (km2/km2)    (t/ha/dec)   (t/km2/dec)    (t/ha/dec)  ',
     &     '  (t/ha/dec)   (t/km2/dec)    (t/ha/dec)    (t/ha/dec)  ',
     &     '  (t/ha/dec)    (t/ha/dec)  ',
     &     '                                                   ',
     &     '   (km2/km2)    (t/ha/dec)   (t/km2/dec)    (t/ha/dec)  ',
     &     '  (t/ha/dec)   (t/km2/dec)    (t/ha/dec)    (t/ha/dec)  ',
     &     '  (t/ha/dec)    (t/ha/dec)  ',
     &     '                                                   ',
     &     '   (km2/km2)    (t/ha/dec)   (t/km2/dec)    (t/ha/dec)  ',
     &     '  (t/ha/dec)   (t/km2/dec)    (t/ha/dec)    (t/ha/dec)  ',
     &     '  (t/ha/dec)    (t/ha/dec)  ',
     &     '                                                   ',
     &     '   (km2/km2)    (t/ha/dec)   (t/km2/dec)    (t/ha/dec)  ',
     &     '  (t/ha/dec)   (t/km2/dec)    (t/ha/dec)    (t/ha/dec)  ',
     &     '  (t/ha/dec)    (t/ha/dec)  ')
C
C Open channel for start data and write titles
C
	TITLE1(1,1,1)='S1LU1L1 '
	TITLE1(2,1,1)='S2LU1L1 '
	TITLE1(3,1,1)='S3LU1L1 '
	TITLE1(4,1,1)='S4LU1L1 '
	TITLE1(5,1,1)='S5LU1L1 '
	TITLE1(1,2,1)='S1LU2L1 '
	TITLE1(2,2,1)='S2LU2L1 '
	TITLE1(3,2,1)='S3LU2L1 '
	TITLE1(4,2,1)='S4LU2L1 '
	TITLE1(5,2,1)='S5LU2L1 '
	TITLE1(1,3,1)='S1LU3L1 '
	TITLE1(2,3,1)='S2LU3L1 '
	TITLE1(3,3,1)='S3LU3L1 '
	TITLE1(4,3,1)='S4LU3L1 '
	TITLE1(5,3,1)='S5LU3L1 '
	TITLE1(1,4,1)='S1LU4L1 '
	TITLE1(2,4,1)='S2LU4L1 '
	TITLE1(3,4,1)='S3LU4L1 '
	TITLE1(4,4,1)='S4LU4L1 '
	TITLE1(5,4,1)='S5LU4L1 '
	TITLE1(1,5,1)='S1LU5L1 '
	TITLE1(2,5,1)='S2LU5L1 '
	TITLE1(3,5,1)='S3LU5L1 '
	TITLE1(4,5,1)='S4LU5L1 '
	TITLE1(5,5,1)='S5LU5L1 '
	TITLE1(1,6,1)='S1LU6L1 '
	TITLE1(2,6,1)='S2LU6L1 '
	TITLE1(3,6,1)='S3LU6L1 '
	TITLE1(4,6,1)='S4LU6L1 '
	TITLE1(5,6,1)='S5LU6L1 '
	TITLE1(1,1,2)='S1LU1L2 '
	TITLE1(2,1,2)='S2LU1L2 '
	TITLE1(3,1,2)='S3LU1L2 '
	TITLE1(4,1,2)='S4LU1L2 '
	TITLE1(5,1,2)='S5LU1L2 '
	TITLE1(1,2,2)='S1LU2L2 '
	TITLE1(2,2,2)='S2LU2L2 '
	TITLE1(3,2,2)='S3LU2L2 '
	TITLE1(4,2,2)='S4LU2L2 '
	TITLE1(5,2,2)='S5LU2L2 '
	TITLE1(1,3,2)='S1LU3L2 '
	TITLE1(2,3,2)='S2LU3L2 '
	TITLE1(3,3,2)='S3LU3L2 '
	TITLE1(4,3,2)='S4LU3L2 '
	TITLE1(5,3,2)='S5LU3L2 '
	TITLE1(1,4,2)='S1LU4L2 '
	TITLE1(2,4,2)='S2LU4L2 '
	TITLE1(3,4,2)='S3LU4L2 '
	TITLE1(4,4,2)='S4LU4L2 '
	TITLE1(5,4,2)='S5LU4L2 '
	TITLE1(1,5,2)='S1LU5L2 '
	TITLE1(2,5,2)='S2LU5L2 '
	TITLE1(3,5,2)='S3LU5L2 '
	TITLE1(4,5,2)='S4LU5L2 '
	TITLE1(5,5,2)='S5LU5L2 '
	TITLE1(1,6,2)='S1LU6L2 '
	TITLE1(2,6,2)='S2LU6L2 '
	TITLE1(3,6,2)='S3LU6L2 '
	TITLE1(4,6,2)='S4LU6L2 '
	TITLE1(5,6,2)='S5LU6L2 '
	TITLE1(1,1,3)='S1LU1L3 '
	TITLE1(2,1,3)='S2LU1L3 '
	TITLE1(3,1,3)='S3LU1L3 '
	TITLE1(4,1,3)='S4LU1L3 '
	TITLE1(5,1,3)='S5LU1L3 '
	TITLE1(1,2,3)='S1LU2L3 '
	TITLE1(2,2,3)='S2LU2L3 '
	TITLE1(3,2,3)='S3LU2L3 '
	TITLE1(4,2,3)='S4LU2L3 '
	TITLE1(5,2,3)='S5LU2L3 '
	TITLE1(1,3,3)='S1LU3L3 '
	TITLE1(2,3,3)='S2LU3L3 '
	TITLE1(3,3,3)='S3LU3L3 '
	TITLE1(4,3,3)='S4LU3L3 '
	TITLE1(5,3,3)='S5LU3L3 '
	TITLE1(1,4,3)='S1LU4L3 '
	TITLE1(2,4,3)='S2LU4L3 '
	TITLE1(3,4,3)='S3LU4L3 '
	TITLE1(4,4,3)='S4LU4L3 '
	TITLE1(5,4,3)='S5LU4L3 '
	TITLE1(1,5,3)='S1LU5L3 '
	TITLE1(2,5,3)='S2LU5L3 '
	TITLE1(3,5,3)='S3LU5L3 '
	TITLE1(4,5,3)='S4LU5L3 '
	TITLE1(5,5,3)='S5LU5L3 '
	TITLE1(1,6,3)='S1LU6L3 '
	TITLE1(2,6,3)='S2LU6L3 '
	TITLE1(3,6,3)='S3LU6L3 '
	TITLE1(4,6,3)='S4LU6L3 '
	TITLE1(5,6,3)='S5LU6L3 '
	TITLE1(1,1,4)='S1LU1L4 '
	TITLE1(2,1,4)='S2LU1L4 '
	TITLE1(3,1,4)='S3LU1L4 '
	TITLE1(4,1,4)='S4LU1L4 '
	TITLE1(5,1,4)='S5LU1L4 '
	TITLE1(1,2,4)='S1LU2L4 '
	TITLE1(2,2,4)='S2LU2L4 '
	TITLE1(3,2,4)='S3LU2L4 '
	TITLE1(4,2,4)='S4LU2L4 '
	TITLE1(5,2,4)='S5LU2L4 '
	TITLE1(1,3,4)='S1LU3L4 '
	TITLE1(2,3,4)='S2LU3L4 '
	TITLE1(3,3,4)='S3LU3L4 '
	TITLE1(4,3,4)='S4LU3L4 '
	TITLE1(5,3,4)='S5LU3L4 '
	TITLE1(1,4,4)='S1LU4L4 '
	TITLE1(2,4,4)='S2LU4L4 '
	TITLE1(3,4,4)='S3LU4L4 '
	TITLE1(4,4,4)='S4LU4L4 '
	TITLE1(5,4,4)='S5LU4L4 '
	TITLE1(1,5,4)='S1LU5L4 '
	TITLE1(2,5,4)='S2LU5L4 '
	TITLE1(3,5,4)='S3LU5L4 '
	TITLE1(4,5,4)='S4LU5L4 '
	TITLE1(5,5,4)='S5LU5L4 '
	TITLE1(1,6,4)='S1LU6L4 '
	TITLE1(2,6,4)='S2LU6L4 '
	TITLE1(3,6,4)='S3LU6L4 '
	TITLE1(4,6,4)='S4LU6L4 '
	TITLE1(5,6,4)='S5LU6L5 '
	TITLE1(1,1,5)='S1LU1L5 '
	TITLE1(2,1,5)='S2LU1L5 '
	TITLE1(3,1,5)='S3LU1L5 '
	TITLE1(4,1,5)='S4LU1L5 '
	TITLE1(5,1,5)='S5LU1L5 '
	TITLE1(1,2,5)='S1LU2L5 '
	TITLE1(2,2,5)='S2LU2L5 '
	TITLE1(3,2,5)='S3LU2L5 '
	TITLE1(4,2,5)='S4LU2L5 '
	TITLE1(5,2,5)='S5LU2L5 '
	TITLE1(1,3,5)='S1LU3L5 '
	TITLE1(2,3,5)='S2LU3L5 '
	TITLE1(3,3,5)='S3LU3L5 '
	TITLE1(4,3,5)='S4LU3L5 '
	TITLE1(5,3,5)='S5LU3L5 '
	TITLE1(1,4,5)='S1LU4L5 '
	TITLE1(2,4,5)='S2LU4L5 '
	TITLE1(3,4,5)='S3LU4L5 '
	TITLE1(4,4,5)='S4LU4L5 '
	TITLE1(5,4,5)='S5LU4L5 '
	TITLE1(1,5,5)='S1LU5L5 '
	TITLE1(2,5,5)='S2LU5L5 '
	TITLE1(3,5,5)='S3LU5L5 '
	TITLE1(4,5,5)='S4LU5L5 '
	TITLE1(5,5,5)='S5LU5L5 '
	TITLE1(1,6,5)='S1LU6L5 '
	TITLE1(2,6,5)='S2LU6L5 '
	TITLE1(3,6,5)='S3LU6L5 '
	TITLE1(4,6,5)='S4LU6L5 '
	TITLE1(5,6,5)='S5LU6L5 '
	TITLE1(1,1,6)='S1LU1L6 '
	TITLE1(2,1,6)='S2LU1L6 '
	TITLE1(3,1,6)='S3LU1L6 '
	TITLE1(4,1,6)='S4LU1L6 '
	TITLE1(5,1,6)='S5LU1L6 '
	TITLE1(1,2,6)='S1LU2L6 '
	TITLE1(2,2,6)='S2LU2L6 '
	TITLE1(3,2,6)='S3LU2L6 '
	TITLE1(4,2,6)='S4LU2L6 '
	TITLE1(5,2,6)='S5LU2L6 '
	TITLE1(1,3,6)='S1LU3L6 '
	TITLE1(2,3,6)='S2LU3L6 '
	TITLE1(3,3,6)='S3LU3L6 '
	TITLE1(4,3,6)='S4LU3L6 '
	TITLE1(5,3,6)='S5LU3L6 '
	TITLE1(1,4,6)='S1LU4L6 '
	TITLE1(2,4,6)='S2LU4L6 '
	TITLE1(3,4,6)='S3LU4L6 '
	TITLE1(4,4,6)='S4LU4L6 '
	TITLE1(5,4,6)='S5LU4L6 '
	TITLE1(1,5,6)='S1LU5L6 '
	TITLE1(2,5,6)='S2LU5L6 '
	TITLE1(3,5,6)='S3LU5L6 '
	TITLE1(4,5,6)='S4LU5L6 '
	TITLE1(5,5,6)='S5LU5L6 '
	TITLE1(1,6,6)='S1LU6L6 '
	TITLE1(2,6,6)='S2LU6L6 '
	TITLE1(3,6,6)='S3LU6L6 '
	TITLE1(4,6,6)='S4LU6L6 '
	TITLE1(5,6,6)='S5LU6L6 '
	TITLE1(1,1,7)='S1LU1L7 '
	TITLE1(2,1,7)='S2LU1L7 '
	TITLE1(3,1,7)='S3LU1L7 '
	TITLE1(4,1,7)='S4LU1L7 '
	TITLE1(5,1,7)='S5LU1L7 '
	TITLE1(1,2,7)='S1LU2L7 '
	TITLE1(2,2,7)='S2LU2L7 '
	TITLE1(3,2,7)='S3LU2L7 '
	TITLE1(4,2,7)='S4LU2L7 '
	TITLE1(5,2,7)='S5LU2L7 '
	TITLE1(1,3,7)='S1LU3L7 '
	TITLE1(2,3,7)='S2LU3L7 '
	TITLE1(3,3,7)='S3LU3L7 '
	TITLE1(4,3,7)='S4LU3L7 '
	TITLE1(5,3,7)='S5LU3L7 '
	TITLE1(1,4,7)='S1LU4L7 '
	TITLE1(2,4,7)='S2LU4L7 '
	TITLE1(3,4,7)='S3LU4L7 '
	TITLE1(4,4,7)='S4LU4L7 '
	TITLE1(5,4,7)='S5LU4L7 '
	TITLE1(1,5,7)='S1LU5L7 '
	TITLE1(2,5,7)='S2LU5L7 '
	TITLE1(3,5,7)='S3LU5L7 '
	TITLE1(4,5,7)='S4LU5L7 '
	TITLE1(5,5,7)='S5LU5L7 '
	TITLE1(1,6,7)='S1LU6L7 '
	TITLE1(2,6,7)='S2LU6L7 '
	TITLE1(3,6,7)='S3LU6L7 '
	TITLE1(4,6,7)='S4LU6L7 '
	TITLE1(5,6,7)='S5LU6L7 '
	TITLE1(1,1,8)='S1LU1L8 '
	TITLE1(2,1,8)='S2LU1L8 '
	TITLE1(3,1,8)='S3LU1L8 '
	TITLE1(4,1,8)='S4LU1L8 '
	TITLE1(5,1,8)='S5LU1L8 '
	TITLE1(1,2,8)='S1LU2L8 '
	TITLE1(2,2,8)='S2LU2L8 '
	TITLE1(3,2,8)='S3LU2L8 '
	TITLE1(4,2,8)='S4LU2L8 '
	TITLE1(5,2,8)='S5LU2L8 '
	TITLE1(1,3,8)='S1LU3L8 '
	TITLE1(2,3,8)='S2LU3L8 '
	TITLE1(3,3,8)='S3LU3L8 '
	TITLE1(4,3,8)='S4LU3L8 '
	TITLE1(5,3,8)='S5LU3L8 '
	TITLE1(1,4,8)='S1LU4L8 '
	TITLE1(2,4,8)='S2LU4L8 '
	TITLE1(3,4,8)='S3LU4L8 '
	TITLE1(4,4,8)='S4LU4L8 '
	TITLE1(5,4,8)='S5LU4L8 '
	TITLE1(1,5,8)='S1LU5L8 '
	TITLE1(2,5,8)='S2LU5L8 '
	TITLE1(3,5,8)='S3LU5L8 '
	TITLE1(4,5,8)='S4LU5L8 '
	TITLE1(5,5,8)='S5LU5L8 '
	TITLE1(1,6,8)='S1LU6L8 '
	TITLE1(2,6,8)='S2LU6L8 '
	TITLE1(3,6,8)='S3LU6L8 '
	TITLE1(4,6,8)='S4LU6L8 '
	TITLE1(5,6,8)='S5LU6L8 '
	TITLE1(1,1,9)='S1LU1L9 '
	TITLE1(2,1,9)='S2LU1L9 '
	TITLE1(3,1,9)='S3LU1L9 '
	TITLE1(4,1,9)='S4LU1L9 '
	TITLE1(5,1,9)='S5LU1L9 '
	TITLE1(1,2,9)='S1LU2L9 '
	TITLE1(2,2,9)='S2LU2L9 '
	TITLE1(3,2,9)='S3LU2L9 '
	TITLE1(4,2,9)='S4LU2L9 '
	TITLE1(5,2,9)='S5LU2L9 '
	TITLE1(1,3,9)='S1LU3L9 '
	TITLE1(2,3,9)='S2LU3L9 '
	TITLE1(3,3,9)='S3LU3L9 '
	TITLE1(4,3,9)='S4LU3L9 '
	TITLE1(5,3,9)='S5LU3L9 '
	TITLE1(1,4,9)='S1LU4L9 '
	TITLE1(2,4,9)='S2LU4L9 '
	TITLE1(3,4,9)='S3LU4L9 '
	TITLE1(4,4,9)='S4LU4L9 '
	TITLE1(5,4,9)='S5LU4L9 '
	TITLE1(1,5,9)='S1LU5L9 '
	TITLE1(2,5,9)='S2LU5L9 '
	TITLE1(3,5,9)='S3LU5L9 '
	TITLE1(4,5,9)='S4LU5L9 '
	TITLE1(5,5,9)='S5LU5L9 '
	TITLE1(1,6,9)='S1LU6L9 '
	TITLE1(2,6,9)='S2LU6L9 '
	TITLE1(3,6,9)='S3LU6L9 '
	TITLE1(4,6,9)='S4LU6L9 '
	TITLE1(5,6,9)='S5LU6L9 '
	TITLE1(1,1,10)='S1LU1L10'
	TITLE1(2,1,10)='S2LU1L10'
	TITLE1(3,1,10)='S3LU1L10'
	TITLE1(4,1,10)='S4LU1L10'
	TITLE1(5,1,10)='S5LU1L10'
	TITLE1(1,2,10)='S1LU2L10'
	TITLE1(2,2,10)='S2LU2L10'
	TITLE1(3,2,10)='S3LU2L10'
	TITLE1(4,2,10)='S4LU2L10'
	TITLE1(5,2,10)='S5LU2L10'
	TITLE1(1,3,10)='S1LU3L10'
	TITLE1(2,3,10)='S2LU3L10'
	TITLE1(3,3,10)='S3LU3L10'
	TITLE1(4,3,10)='S4LU3L10'
	TITLE1(5,3,10)='S5LU3L10'
	TITLE1(1,4,10)='S1LU4L10'
	TITLE1(2,4,10)='S2LU4L10'
	TITLE1(3,4,10)='S3LU4L10'
	TITLE1(4,4,10)='S4LU4L10'
	TITLE1(5,4,10)='S5LU4L10'
	TITLE1(1,5,10)='S1LU5L10'
	TITLE1(2,5,10)='S2LU5L10'
	TITLE1(3,5,10)='S3LU5L10'
	TITLE1(4,5,10)='S4LU5L10'
	TITLE1(5,5,10)='S5LU5L10'
	TITLE1(1,6,10)='S1LU6L10'
	TITLE1(2,6,10)='S2LU6L10'
	TITLE1(3,6,10)='S3LU6L10'
	TITLE1(4,6,10)='S4LU6L10'
	TITLE1(5,6,10)='S5LU6L10'
	TITLE1(1,1,10)='S1LU1L10'
	TITLE1(2,1,10)='S2LU1L10'
	TITLE1(3,1,10)='S3LU1L10'
	TITLE1(4,1,10)='S4LU1L10'
	TITLE1(5,1,10)='S5LU1L10'
	TITLE1(1,2,10)='S1LU2L10'
	TITLE1(2,2,10)='S2LU2L10'
	TITLE1(3,2,10)='S3LU2L10'
	TITLE1(4,2,10)='S4LU2L10'
	TITLE1(5,2,10)='S5LU2L10'
	TITLE1(1,3,10)='S1LU3L10'
	TITLE1(2,3,10)='S2LU3L10'
	TITLE1(3,3,10)='S3LU3L10'
	TITLE1(4,3,10)='S4LU3L10'
	TITLE1(5,3,10)='S5LU3L10'
	TITLE1(1,4,10)='S1LU4L10'
	TITLE1(2,4,10)='S2LU4L10'
	TITLE1(3,4,10)='S3LU4L10'
	TITLE1(4,4,10)='S4LU4L10'
	TITLE1(5,4,10)='S5LU4L10'
	TITLE1(1,5,10)='S1LU5L10'
	TITLE1(2,5,10)='S2LU5L10'
	TITLE1(3,5,10)='S3LU5L10'
	TITLE1(4,5,10)='S4LU5L10'
	TITLE1(5,5,10)='S5LU5L10'
	TITLE1(1,6,10)='S1LU6L10'
	TITLE1(2,6,10)='S2LU6L10'
	TITLE1(3,6,10)='S3LU6L10'
	TITLE1(4,6,10)='S4LU6L10'
	TITLE1(5,6,10)='S5LU6L10'
	TITLE3(1)='    S1    '
	TITLE3(2)='    S2    '
	TITLE3(3)='    S3    '
	TITLE3(4)='    S4    '
	TITLE3(5)='    S4    '
	TITLE4(1)='   LU1    '
	TITLE4(2)='   LU2    '
	TITLE4(3)='   LU3    '
	TITLE4(4)='   LU4    '
	TITLE4(5)='   LU5    '
	TITLE4(6)='   LU6    '
      OPEN(STARTCHAN,FILE='STARTDATA.TXT',STATUS='UNKNOWN')
	WRITE(STARTCHAN,910)
910   FORMAT('20km2_ID    1km2_ID        Easting     Northing      ',
     &       '  Area    SoilC(kgC/ha) to 100cm')
C	WRITE(STARTCHAN,910)((((TITLE1(ISERIES,ILU1,IL),
C     &                        TITLE1(ISERIES,ILU1,IL)),
C     &                        IL=1,MAXSOMLAY),
C     &                        ILU1=1,MAXLU),
C     &                        ISERIES=1,MAXSERIES),
C     &                        (TITLE3(ISERIES),ISERIES=1,MAXSERIES),
C     &                        (TITLE4(ILU1),ILU1=1,MAXLU)
C910   FORMAT('Note: S=soil LU=land use L=layer'/
C     &       '  20km2 ID  1km2 ID   Easting    Northing      Area    ',
C     &       300('Depth(cm)  SoilC(kgC/ha)')/
C     &       55X,600(A10,2X),
C     &       5(A10,2X),6('  LU   ',A10,2X))
C
C Leave OPENCHAN_GIS
C
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE PUT_GIS_RES(SQUID,EAST,NORTH,NUMCELLSIN20KM2,
     &                       KM20GRIDID,DOMSOIL,PERSOIL,SOILID,STYPES,
     &                       WETCLASS,TOTSEQ,LUSEQ,SEQTYPE,
     &                       LU1TOLU2,CCHANGE10,CH4C10,CO2C10,N2ON10,
     &                       TRAPERR,DOMSOILC,DOMSOILBD,DOMSOILCLAY,
     &                       SOMDEPTH,gley,PIANN,lanu)	 
C
C Subroutine to record results of spatial simulation
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      INTEGER MAXDEC				! Max.no.of decades included in calculation
	PARAMETER(MAXDEC=9)		
	INTEGER MAXDEC1				! Max.no.of decades included in calculation
	DATA MAXDEC1 /9/
	INTEGER MAXGROW				! Max.no.of growing seasons allowed
	PARAMETER (MAXGROW=300)
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXLU1				! Max.no.of land use types
	DATA MAXLU1 /6/
	INTEGER MAXSEQ				! Max.no of sequences included in GIS
	PARAMETER (MAXSEQ=MAXLU*MAXLU)						
	INTEGER MAXSERIES			! Max.no of soil series
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXSOMLAY			! Max.no.of soil layers entered			
	PARAMETER (MAXSOMLAY=10)
	INTEGER MAXWC				! Max.no of wetness classes
      PARAMETER (MAXWC=6)

      REAL AWC					! Atomic weight of C g/mol
	REAL AWN					! Atomic weight of N g/mol
	DATA AWC,AWN /12,14/
	REAL CCHGE(MAXDEC)			! C Change for 1km2 cell t/km2/10y
	REAL CELLCCHANGE(MAXDEC)	! C Change for cell in t/20km2/10 yrs
	REAL CELLCEQUIV(MAXDEC)		! C equiv.emmision for cell in t/20km2/10 yrs
      CHARACTER*50 CELLERROR		! Error Message for 1km2 cell
	REAL CELLCH4C(MAXDEC)		! CH4C emissions for cell in t/20km2/10 yrs
	REAL CELLCO2C(MAXDEC)		! CO2C emissions for cell in t/20km2/10 yrs
	REAL CELLN2ON(MAXDEC)		! N2ON emissions for cell in t/20km2/10 yrs
	REAL CEQ(MAXDEC)			! C equiv.emmision for 1km2 cell in t/km2/10y
	REAL CLU1TOLU2(MAXDEC,MAXLU+2,MAXLU+2) ! IN:Carbon change for LU1 to LU2 
								! decade 1 to 6 in ktC/km2/10y
	REAL C20KM2LU1TOLU2(MAXDEC,MAXLU+2,MAXLU+2) ! IN:Carbon change LU1 to LU2 
								! decade 1 to 6 in tC/20km2/10y
	INTEGER DECERROR(MAXDEC)	! Error found in this decade
      INTEGER EQUIERR
	DATA EQUIERR /6/			!		2=Error in equilibrium run
	REAL FRACCELL(MAXSERIES,MAXLU,MAXLU) ! Fraction of cell in soil ISERIES and 
      INTEGER OFRESCHAN			! Channel for results of future met.data inputs
	DATA OFRESCHAN /67/
	INTEGER GRIDCHAN			! Channel for file for 20km2 grid outputs
	DATA GRIDCHAN /38/
	INTEGER OGRIDCHAN			! Channel for file for 20km2 grid outputs - organic soils
	DATA OGRIDCHAN /48/
      CHARACTER*50 GRIDERROR		! Error Message for 20km2 grid
      REAL GWPCO2					! Global warming potential of CO2 /mol
      REAL GWPCH4					! Global warming potential of CO2 /mol
      REAL GWPN2O					! Global warming potential of CO2 /mol
	DATA GWPCO2,GWPCH4,GWPN2O /1,23,296/
	INTEGER IBACK				! Local counter of number of decades previous 
								! to current when change occurred
	INTEGER ICELL				! Local cell counter
	INTEGER ICHAN				! Channel counter
      INTEGER IDEC				! Decade counter
	INTEGER ISORGANIC			! Marker for org.soil (0=mineral,1=organic)
	INTEGER IWC					! Local counter for wetness class
	INTEGER LUCERROR(MAXDEC,MAXLU+2,MAXLU+2) ! Error found in this LU change
	REAL LU1TOLU220KM2(MAXDEC,MAXLU+2,MAXLU+2) ! LU change km2/20km2
	INTEGER ISERIES				! Local counter for soil series
	INTEGER LU1					! Land use 1
	INTEGER LU2					! Land use 2
	INTEGER METERR
	DATA METERR /2/				!		2=Error in weather data
	REAL MWCH4					! Molecular weight of CH4 g/mol
	REAL MWCO2					! Molecular weight of CO2 g/mol
	REAL MWN2O					! Molecular weight of N2O g/mol
	DATA MWCO2,MWCH4,MWN2O /44,16,44/
	INTEGER NPPERR
	DATA NPPERR /1/				!		1=Error in NPP data
	INTEGER NSEQ				! Counter for number of land use sequences 
								! undergoing land use change given by NSEQ in given decade
	INTEGER OUTCHAN				! Channel for GIS output
      DATA OUTCHAN /37/
	INTEGER OOUTCHAN,lanu			! Channel for GIS output - organic soils
      DATA OOUTCHAN /47/
	REAL OCCHGE(MAXDEC)			! C Change on org.soil for 1km2 cell t/km2/10y
	REAL OCEQ(MAXDEC)			! C equiv.emmision on org.soil in t/km2/10y
	REAL OCLU1TOLU2(MAXDEC,MAXLU+2,MAXLU+2) ! IN:Carbon change on organic 
								! soils for LU1 to LU2 in decade 1 to 6 in ktC/km2/10y
	REAL OC20KM2LU1TOLU2(MAXDEC,MAXLU+2,MAXLU+2) ! IN:Carbon change LU1 to LU2 
								! on organic soils in decade 1 to 6 in tC/20km2/10y
	REAL OCELLCCHANGE(MAXDEC)	! C Change on org.soil for cell in t/20km2/10y
	REAL OCELLCH4C(MAXDEC)		! CH4C emissions on org.soil in t/20km2/10y
	REAL OCELLCO2C(MAXDEC)		! CO2C emissions on org.soil in t/20km2/10y
	REAL OCELLN2ON(MAXDEC)		! N2ON emissions on org.soil in t/20km2/10y
	REAL OCELLCEQUIV(MAXDEC)	! C equiv. on org.soil for cell in t/20km2/10y
	REAL OLU1TOLU2(CELLSTOGRID,MAXDEC,MAXLU+2,MAXLU+2) ! IN:Frac of LU1 km2/dec/km2
								! changed to LU2 in decade on organic soil
	REAL OLU1TOLU220KM2(MAXDEC,MAXLU+2,MAXLU+2) ! LU change km2/20km2 - organic soil
	REAL PC						! % C in the top soil layer
	REAL PCLAY(MAXSERIES,MAXLU)	! % Clay in this series and lu
      Real Gley(Maxlu)			! Gleying layer...
	REAL CCH(MAXDEC,MAXSERIES,MAXSEQ)	! CChange in this decade and series and lu change
	REAL PERCENTC(MAXSERIES,MAXLU)	! % C in the top soil layer for given series and LU
	REAL TEMP					! Temporary real for calculation
C
C Variables passed to/from calling subroutine
C
	REAL CCHANGE10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ) ! IN:Change 
								! in C content over the decade (kgC/ha/decade)
	REAL CH4C10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! IN:CH4-C 
								! emitted over the decade (kgC/ha/decade)
	REAL CO2C10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! IN:CO2-C 
								! emitted over the decade (kgC/ha/decade)
      INTEGER DOMSOIL(CELLSTOGRID,MAXSERIES)	! IN:5 major soil series in the  		
								! 1km2 cell x no. of cells in 20km2 grid
	REAL DOMSOILBD(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! Soil BD in 5 
												! major soil series under different LU 
												! in SOM layers (g/cm3)
	REAL DOMSOILC(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! Soil C in 5 
												! major soil series under different LU 
												! in SOM layers (kgC/ha)
	REAL DOMSOILCLAY(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! % Soil clay in 5 
												! major soil series under different LU 
	REAL EAST(CELLSTOGRID),PIANN		! IN:National grid easting x 1000
	CHARACTER*10 KM20GRIDID		! IN:20km sqare identifier
	REAL LU1TOLU2(CELLSTOGRID,MAXDEC,MAXLU+2,MAXLU+2) ! IN:Frac of LU1 km2/dec/km2
								! changed to LU2 in decade 
      INTEGER LUSEQ(MAXGROW)		! IN:Land use sequence
	REAL N2ON10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! IN:N2O-N 
								! emitted over the decade (kgN/ha/decade)
	INTEGER NGROW				! IN:No.growing seasons simulated
	REAL NORTH(CELLSTOGRID)		! IN:National grid northing x 1000
	INTEGER NUMCELLSIN20KM2		! IN:Number of cells in 20km2 cell
      INTEGER ORDER(MAXSERIES)	!IN:Order of soils calculations
	REAL PERSOIL(CELLSTOGRID,MAXSERIES)		! IN:Percentage of each major soil 
								!	type in the 1km2 square cell
	REAL SOMDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN:depth of soil organic matter layers(cm)
	REAL SQUID(CELLSTOGRID)		! IN:1km square identifier
	INTEGER SEQ1				! IN:Land use sequence of LU1->LU1
	INTEGER SEQTYPE				! IN:Type of sequence to be used
      INTEGER SOILID(CELLSTOGRID*MAXSERIES/2)	! IN:Soil series integer codes
	INTEGER STYPES				! IN:Number of soil types in 20km2 grid
	INTEGER TOTSEQ				! IN:Total number of sequences
      INTEGER TRAPERR(CELLSTOGRID*MAXSERIES/2,MAXSEQ) 
								! IN:Error trapping
								!		0=No error
								!		1=Error in NPP data
								!		2=Error in weather data
								!		3=Error in LU data
								!		4=Error in LU change data
								!		5=Soil error
								!		6=Equilibrium not found
      INTEGER WETCLASS(CELLSTOGRID,MAXSERIES) !OUT:Wetness class for 
								! each major soil type in the 1km2 square cell 
								! x 20km2 cell
C
C Initialise 20km2 results
C
      DO 100 IDEC=1,MAXDEC1
        DO 200 LU1=1,MAXLU1
	    DO 300 LU2=1,MAXLU1
            C20KM2LU1TOLU2(IDEC,LU1,LU2)=0
            OC20KM2LU1TOLU2(IDEC,LU1,LU2)=0
            LU1TOLU220KM2(IDEC,LU1,LU2)=0
            OLU1TOLU220KM2(IDEC,LU1,LU2)=0
300       CONTINUE	! LU2 loop
200     CONTINUE		! LU1 Loop
        CCHGE(IDEC)=0
	  OCCHGE(IDEC)=0
        CEQ(IDEC)=0
        OCEQ(IDEC)=0
100   CONTINUE		! IDEC loop
C
C Calculate overall changes for the cell
C
	GRIDERROR='        '
      DO 400 ICELL=1,NUMCELLSIN20KM2
        CELLERROR='         '
	  DO 500 IDEC=1,MAXDEC1
          CELLCCHANGE(IDEC)=0
          OCELLCCHANGE(IDEC)=0
	    OCELLCEQUIV(IDEC)=0
	    CELLCH4C(IDEC)=0
	    CELLCO2C(IDEC)=0
  	    CELLN2ON(IDEC)=0
	    OCELLCH4C(IDEC)=0
	    OCELLCO2C(IDEC)=0
  	    OCELLN2ON(IDEC)=0
	    DECERROR(IDEC)=0
	    DO 550 LU1=1,MAXLU1
	      DO 575 LU2=1,MAXLU1
              LUCERROR(IDEC,LU1,LU2)=0
              CLU1TOLU2(IDEC,LU1,LU2)=0
              OCLU1TOLU2(IDEC,LU1,LU2)=0
575         CONTINUE	! LU2 loop
550       CONTINUE	! LU1 loop
C
C For each soils series that exists in this cell...
C 
          DO 600 ISERIES=1,MAXSERIES
	      IF(PERSOIL(ICELL,ISERIES).GT.0)THEN
C
C ...get order of soils calculations
C
              CALL GET_SOIL_ORDER(DOMSOIL(ICELL,ISERIES),SOILID,
     &                            STYPES,ORDER(ISERIES))
C
C If this soil series is found...
C
              IF(ORDER(ISERIES).GT.0)THEN
C
C ...for each sequence in cell
C
	          DO 700 NSEQ=1,TOTSEQ
	            IF(TRAPERR(ORDER(ISERIES),NSEQ).EQ.EQUIERR)THEN
                    CELLERROR='Not initialised to equilibrium'
                    GRIDERROR='Not initialised to equilibrium'
			    ENDIF			! IF(TRAPERR(ORDER(ISERIES),NSEQ).EQ.EQUIERR)THEN
	            CALL GET_LUSEQ(LUSEQ,NGROW,NSEQ,SEQTYPE)
	            LU1=LUSEQ(1)
	            LU2=LUSEQ(6)
	            SEQ1=((LU1-1)*MAXLU1)+LU1
C
C ...check top soil layer for an organic soil
C
                  ISORGANIC=0
                  PC=DOMSOILC(ORDER(ISERIES),LU1,1)
				PC=PC/DOMSOILBD(ORDER(ISERIES),LU1,1)
				PC=PC/(SOMDEPTH(ORDER(ISERIES),LU1,1)*1000)
	            PERCENTC(ISERIES,LU1)=PC
                  IF(PC.GT.6)ISORGANIC=1
C 
C ...if there IS an error in results, record the error for this decade and LU change
C
	            IF(TRAPERR(ORDER(ISERIES),NSEQ).EQ.NPPERR.OR.
     &                   TRAPERR(ORDER(ISERIES),NSEQ).EQ.METERR)THEN
                    LUCERROR(IDEC,LU1,LU2)=1
	              DECERROR(IDEC)=1
	              IF(TRAPERR(ORDER(ISERIES),NSEQ).EQ.NPPERR)THEN
                      CELLERROR='Error in NPP'
                      GRIDERROR='Error in NPP'
	              ENDIF				! IF(TRAPERR(ORDER(ISERIES),NSEQ).EQ.NPPERR)THEN
	              IF(TRAPERR(ORDER(ISERIES),NSEQ).EQ.METERR)THEN
                      CELLERROR='Error in met.data'
                      GRIDERROR='Error in met.data'
	              ENDIF				! IF(TRAPERR(ORDER(ISERIES),NSEQ).EQ.METERR)THEN
C
C ...if there is no met or NPP error in results on this series and this sequence
C
	            ELSE				! ELSE to IF(TRAPERR(ORDER(ISERIES),NSEQ).EQ.NPPERR.OR.
									! TRAPERR(ORDER(ISERIES),NSEQ).EQ.METERR)THEN
C
C ......Add up land use change in this decade and sequence (cumulate over all cells in 20km2 grid)
C
                    TEMP=LU1TOLU220KM2(IDEC,LU1,LU2)
                    TEMP=TEMP+LU1TOLU2(ICELL,IDEC,LU1,LU2)
				  LU1TOLU220KM2(IDEC,LU1,LU2)=TEMP
	              IF(ISORGANIC.EQ.1)THEN
				    OLU1TOLU220KM2(IDEC,LU1,LU2)=TEMP
	        OLU1TOLU2(ICELL,IDEC,LU1,LU2)=LU1TOLU2(ICELL,IDEC,LU1,LU2)
	              ELSEIF(ISORGANIC.EQ.1)THEN
	                OLU1TOLU2(ICELL,IDEC,LU1,LU2)=0
	              ENDIF
C
C ...for each decade 
C
				  DO 750 IBACK=1,IDEC
C
C ...add in all the contributions to C change from previous decades
C ...and if LU change occurs...
C
	                IF(LU1TOLU2(ICELL,IBACK,LU1,LU2).GT.0)THEN
C
C
C ...get fraction of cell on this soil under this sequence and this LU change
C
	                  TEMP=LU1TOLU2(ICELL,IBACK,LU1,LU2) 
	                  TEMP=TEMP*(PERSOIL(ICELL,ISERIES)/100)
				      FRACCELL(ISERIES,LU1,LU2)=TEMP
C
C ...and cumulate results for this cell, for each decade and each LU change
C ......Carbon change in this decade and sequence (cumulate over all soils)
C .........On all soils
                        IWC=WETCLASS(ICELL,ISERIES)
c Check for erroneous C change values
      IF(CCHANGE10(IDEC-IBACK+1,ORDER(ISERIES),IWC,NSEQ).LT.-0.0001.OR.
     &   CCHANGE10(IDEC-IBACK+1,ORDER(ISERIES),IWC,NSEQ).GT.0.0001)THEN
	         TEMP=CCHANGE10(IDEC-IBACK+1,ORDER(ISERIES),IWC,NSEQ)
	         TEMP=TEMP-CCHANGE10(IDEC-IBACK+1,ORDER(ISERIES),IWC,SEQ1)
	         TEMP=TEMP*FRACCELL(ISERIES,LU1,LU2)
	         TEMP=TEMP+CLU1TOLU2(IDEC,LU1,LU2)
	         CLU1TOLU2(IDEC,LU1,LU2)=TEMP	 
C .........On organic soils
	                  IF(ISORGANIC.EQ.1)THEN
	         TEMP=CCHANGE10(IDEC-IBACK+1,ORDER(ISERIES),IWC,NSEQ)
		     TEMP=TEMP-CCHANGE10(IDEC-IBACK+1,ORDER(ISERIES),IWC,SEQ1)
	                    TEMP=TEMP*FRACCELL(ISERIES,LU1,LU2)
			            TEMP=TEMP+OCLU1TOLU2(IDEC,LU1,LU2)
	                    OCLU1TOLU2(IDEC,LU1,LU2)=TEMP
	                  ENDIF		! IF(ISORGANIC.EQ.1)THEN
C ......Methane emissions in this decade cumulated over all soils & sequences
C .........On all soils
			  	TEMP=CH4C10(IDEC-IBACK+1,ORDER(ISERIES),IWC,NSEQ)
				TEMP=TEMP-CH4C10(IDEC-IBACK+1,ORDER(ISERIES),IWC,SEQ1)
					  TEMP=TEMP*FRACCELL(ISERIES,LU1,LU2)
                        CELLCH4C(IDEC)=CELLCH4C(IDEC)+TEMP
C .........On organic soils
	                  IF(ISORGANIC.EQ.1)THEN
                          OCELLCH4C(IDEC)=OCELLCH4C(IDEC)+TEMP
	                  ENDIF		! IF(ISORGANIC.EQ.1)THEN
C ......Carbon dioxide emissions in this decade cumulated over all soils & sequences
C .........On all soils
				TEMP=CO2C10(IDEC-IBACK+1,ORDER(ISERIES),IWC,NSEQ)
			    TEMP=TEMP-CO2C10(IDEC-IBACK+1,ORDER(ISERIES),IWC,SEQ1)
					  TEMP=TEMP*FRACCELL(ISERIES,LU1,LU2)
	                  CELLCO2C(IDEC)=CELLCO2C(IDEC)+TEMP
C .........On organic soils
	                  IF(ISORGANIC.EQ.1)THEN
	                    OCELLCO2C(IDEC)=OCELLCO2C(IDEC)+TEMP
	                  ENDIF		! IF(ISORGANIC.EQ.1)THEN
C ......Nitrous oxide emissions in this decade cumulated cumulated over all soils & sequences
C .........On all soils
				TEMP=N2ON10(IDEC-IBACK+1,ORDER(ISERIES),IWC,NSEQ)
			    TEMP=TEMP-N2ON10(IDEC-IBACK+1,ORDER(ISERIES),IWC,SEQ1)
					  TEMP=TEMP*FRACCELL(ISERIES,LU1,LU2)
                        CELLN2ON(IDEC)=CELLN2ON(IDEC)+TEMP  
C .........On organic soils
	                  IF(ISORGANIC.EQ.1)THEN
                          OCELLN2ON(IDEC)=OCELLN2ON(IDEC)+TEMP  
					  ENDIF		! IF(ISORGANIC.EQ.1)THEN
	                  ENDIF		! IF(CCHANGE10(IDEC-IBACK+1,ORDER(ISERIES),IWC,NSEQ).LT.-0.0001).OR.(CCHANGE10(IDEC-IBACK+1,ORDER(ISERIES),IWC,NSEQ).GT.0.0001)THEN
				    ENDIF		! IF(LU1TOLU2(ICELL,IBACK,LU1,LU2).GT.0)THEN
750                 CONTINUE		! IBACK loop
	            ENDIF			! ELSE to IF(TRAPERR(ORDER(ISERIES),NSEQ).EQ.NPPERR.OR.
								! TRAPERR(ORDER,NSEQ).EQ.METERR)THEN
700		      CONTINUE			! NSEQ loop
              ENDIF				! IF(ORDER(ISERIES).GT.0)THEN
            ENDIF					! IF(PERSOIL(ICELL,ISERIES).GT.0)THEN
600       CONTINUE				! ISERIES loop
500     CONTINUE					! IDEC loop
C
C Recalculate in kt per 1km grid square
C kg/ha -> kt/km2: kg/1000*1000->kt; ha x 100 ->km2; kg/ha x 100/(1000*1000) -> kt/km2 =1/10000
C
        DO 800 IDEC=1,MAXDEC1
	    IF(DECERROR(IDEC).EQ.0)THEN
	      CELLCH4C(IDEC)=CELLCH4C(IDEC)/10000
	      CELLCO2C(IDEC)=CELLCO2C(IDEC)/10000
	      CELLN2ON(IDEC)=CELLN2ON(IDEC)/10000
	      CELLCEQUIV(IDEC)=CELLCEQUIV(IDEC)/10000
	      OCELLCH4C(IDEC)=OCELLCH4C(IDEC)/10000
	      OCELLCO2C(IDEC)=OCELLCO2C(IDEC)/10000
	      OCELLN2ON(IDEC)=OCELLN2ON(IDEC)/10000
	      OCELLCEQUIV(IDEC)=OCELLCEQUIV(IDEC)/10000
	      DO 820 LU1=1,MAXLU1
	        DO 840 LU2=1,MAXLU1
	          IF(LUCERROR(IDEC,LU1,LU2).EQ.0)THEN
      		    CLU1TOLU2(IDEC,LU1,LU2)=CLU1TOLU2(IDEC,LU1,LU2)/10000
      		   OCLU1TOLU2(IDEC,LU1,LU2)=OCLU1TOLU2(IDEC,LU1,LU2)/10000
C
C ...add up all the C changes over all sequences on mineral and organic soils
C
	            TEMP=CLU1TOLU2(IDEC,LU1,LU2)
	            CELLCCHANGE(IDEC)=CELLCCHANGE(IDEC)+TEMP
	            TEMP=OCLU1TOLU2(IDEC,LU1,LU2)
	            OCELLCCHANGE(IDEC)=OCELLCCHANGE(IDEC)+TEMP
	          ENDIF
840           CONTINUE		! LU2 loop
820         CONTINUE			! LU2 loop
C
C Set results for decades with errors in to zero
C
	    ELSE
	      CELLCCHANGE(IDEC)=0
	      OCELLCCHANGE(IDEC)=0
	      CELLCH4C(IDEC)=0
	      CELLCO2C(IDEC)=0
	      CELLN2ON(IDEC)=0
	      CELLCEQUIV(IDEC)=0
	      OCELLCH4C(IDEC)=0
	      OCELLCO2C(IDEC)=0
	      OCELLN2ON(IDEC)=0
	      OCELLCEQUIV(IDEC)=0
	      DO 860 LU1=1,MAXLU1
	        DO 880 LU2=1,MAXLU1
	          IF(LUCERROR(IDEC,LU1,LU2).EQ.1)THEN
                  CLU1TOLU2(IDEC,LU1,LU2)=0
                  OCLU1TOLU2(IDEC,LU1,LU2)=0
	          ENDIF
880           CONTINUE		! LU2 loop
860         CONTINUE			! LU1 loop
          ENDIF
C
C ... Work out C equivalents this decade from gaseous emissions 
C ......On all soils
          CELLCEQUIV(IDEC)=CELLCO2C(IDEC)*GWPCO2*MWCO2/AWC
	    TEMP=CELLCH4C(IDEC)*GWPCH4*MWCH4/AWC
	    CELLCEQUIV(IDEC)=CELLCEQUIV(IDEC)+TEMP
	    TEMP=CELLN2ON(IDEC)*GWPN2O*MWN2O/(2*AWN)
	    CELLCEQUIV(IDEC)=CELLCEQUIV(IDEC)+TEMP
	    CELLCEQUIV(IDEC)=CELLCEQUIV(IDEC)*(AWC/MWCO2)
C ......On organic soils
          OCELLCEQUIV(IDEC)=OCELLCO2C(IDEC)*GWPCO2*MWCO2/AWC
	    TEMP=OCELLCH4C(IDEC)*GWPCH4*MWCH4/AWC
	    OCELLCEQUIV(IDEC)=OCELLCEQUIV(IDEC)+TEMP
	    TEMP=OCELLN2ON(IDEC)*GWPN2O*MWN2O/(2*AWN)
	    OCELLCEQUIV(IDEC)=OCELLCEQUIV(IDEC)+TEMP
	    OCELLCEQUIV(IDEC)=OCELLCEQUIV(IDEC)*(AWC/MWCO2)
800     CONTINUE				! IDEC loop
C
C Output results
C .... Results for each 1km2 cell
C
        WRITE(OUTCHAN,10)KM20GRIDID,SQUID(ICELL),
     &     EAST(ICELL),NORTH(ICELL),gley(1),PIANN,lanu,
     &   sum(DOMSOILC(1,1,:),Mask=DOMSOILC(1,1,:)>0),
     &   sum(DOMSOILC(1,2,:),Mask=DOMSOILC(1,2,:)>0),
     &   sum(DOMSOILC(1,3,:),Mask=DOMSOILC(1,3,:)>0),
     &   sum(DOMSOILC(1,4,:),Mask=DOMSOILC(1,4,:)>0),
     &  maxval(SOMDEPTH(1,1,:)),
     &  maxval(SOMDEPTH(1,2,:)),
     &  maxval(SOMDEPTH(1,3,:)),
     &  maxval(SOMDEPTH(1,4,:)),
     &                (((CLU1TOLU2(IDEC,LU1,LU2),LU1=1,MAXLU1),
     &                         LU2=1,MAXLU1),IDEC=1,MAXDEC1),                               
     &                (CELLCCHANGE(IDEC),IDEC=1,MAXDEC1),
     &                (OCELLCCHANGE(IDEC),IDEC=1,MAXDEC1),
     &			    (CELLCO2C(IDEC),IDEC=1,MAXDEC1),
     &      			(CELLCH4C(IDEC),IDEC=1,MAXDEC1),
     &      			(CELLN2ON(IDEC),IDEC=1,MAXDEC1),
     &      			(CELLCEQUIV(IDEC),IDEC=1,MAXDEC1),
     &      			(OCELLCEQUIV(IDEC),IDEC=1,MAXDEC1),
     &                (((LU1TOLU2(ICELL,IDEC,LU1,LU2),LU1=1,MAXLU1),
     &                           LU2=1,MAXLU1),IDEC=1,MAXDEC1), 
     &                 CELLERROR
        WRITE(OOUTCHAN,10)KM20GRIDID,SQUID(ICELL),
     &                   EAST(ICELL),NORTH(ICELL),gley(1),PIANN,lanu,
     &   sum(DOMSOILC(1,1,:),Mask=DOMSOILC(1,1,:)>0),
     &   sum(DOMSOILC(1,2,:),Mask=DOMSOILC(1,2,:)>0),
     &   sum(DOMSOILC(1,3,:),Mask=DOMSOILC(1,3,:)>0),
     &   sum(DOMSOILC(1,4,:),Mask=DOMSOILC(1,4,:)>0),
     &  maxval(SOMDEPTH(1,1,:)),
     &  maxval(SOMDEPTH(1,2,:)),
     &  maxval(SOMDEPTH(1,3,:)),
     &  maxval(SOMDEPTH(1,4,:)),
     &                (((OCLU1TOLU2(IDEC,LU1,LU2),LU1=1,MAXLU1),
     &                           LU2=1,MAXLU1),IDEC=1,MAXDEC1),                               
     &                (CELLCCHANGE(IDEC),IDEC=1,MAXDEC1),
     &                (OCELLCCHANGE(IDEC),IDEC=1,MAXDEC1),
     &			    (CELLCO2C(IDEC),IDEC=1,MAXDEC1),
     &      			(CELLCH4C(IDEC),IDEC=1,MAXDEC1),
     &      			(CELLN2ON(IDEC),IDEC=1,MAXDEC1),
     &      			(CELLCEQUIV(IDEC),IDEC=1,MAXDEC1),
     &      			(OCELLCEQUIV(IDEC),IDEC=1,MAXDEC1),
     &                (((OLU1TOLU2(ICELL,IDEC,LU1,LU2),LU1=1,MAXLU1),
     &                           LU2=1,MAXLU1),IDEC=1,MAXDEC1), 
     &                 CELLERROR
10      FORMAT(A10,',',F11.0,',',3(F8.0,','),F8.0',',I8.0',',
     &         8(F15.0,','),6(6(9(F10.5,', '))),
     &         7(9(F18.5,', ')),
     &         6(6(9(F10.4,', '))),A50)	
C
C ...Sum 1km2 results over 20km2
C
        DO 900 IDEC=1,MAXDEC1
	    DO 1000 LU1=1,MAXLU1
	      DO 1100 LU2=1,MAXLU1
                TEMP=C20KM2LU1TOLU2(IDEC,LU1,LU2)
                TEMP=TEMP+CLU1TOLU2(IDEC,LU1,LU2)
	          C20KM2LU1TOLU2(IDEC,LU1,LU2)=TEMP
                TEMP=OC20KM2LU1TOLU2(IDEC,LU1,LU2)
                TEMP=TEMP+OCLU1TOLU2(IDEC,LU1,LU2)
	          OC20KM2LU1TOLU2(IDEC,LU1,LU2)=TEMP
1100        CONTINUE		! LU2 loop
1000      CONTINUE		! LU1 loop
          CCHGE(IDEC)=CCHGE(IDEC)+CELLCCHANGE(IDEC)
		OCCHGE(IDEC)=OCCHGE(IDEC)+OCELLCCHANGE(IDEC)
          CEQ(IDEC)=CEQ(IDEC)+CELLCEQUIV(IDEC)
          OCEQ(IDEC)=OCEQ(IDEC)+OCELLCEQUIV(IDEC)
900     CONTINUE			! IDEC loop
C
C ...Results for calculation of mitigation options
C
	  DO 1500 LU1=1,MAXLU1
	    DO 1600 LU2=1,MAXLU1
	      IF(LU1.NE.LU2)THEN
	        NSEQ=((LU1-1)*MAXLU1)+LU2
	        ICHAN=OFRESCHAN+NSEQ
	          DO 1700 ISERIES=1,MAXSERIES
	            DO 1800 IDEC=1,MAXDEC1
	              IF(ORDER(ISERIES).GT.0.AND.
     &                 ORDER(ISERIES).LE.CELLSTOGRID*MAXSERIES/2.AND.
     &                 WETCLASS(ICELL,ISERIES).GT.0.AND.
     &                 WETCLASS(ICELL,ISERIES).LE.MAXWC)THEN
	              PCLAY(ISERIES,LU1)=DOMSOILCLAY(ORDER(ISERIES),LU1,1)
	                 CCH(IDEC,ISERIES,NSEQ)=
     &                                    CCHANGE10(IDEC,ORDER(ISERIES),
     &                                    WETCLASS(ICELL,ISERIES),NSEQ)
     &                                    /1000
	              ELSEIF(ORDER(ISERIES).LE.0)THEN
	                PCLAY(ISERIES,LU1)=-9999.99
	                CCH(IDEC,ISERIES,NSEQ)=-9999.99
                    ENDIF
1800              CONTINUE
1700            CONTINUE
c	        WRITE(ICHAN,1510)KM20GRIDID,
c     &              SQUID(ICELL),EAST(ICELL),NORTH(ICELL),
c     &              DOMSOIL(ICELL,ISERIES),WETCLASS(ICELL,ISERIES),
c     &              PERCENTC(ISERIES,LU1),PCLAY(ISERIES,LU1),
c     &		      FRACCELL(ISERIES,LU1,LU2),
c     &              CCH(1,ISERIES,NSEQ),CCH(2,ISERIES,NSEQ),
c     &              CCH(3,ISERIES,NSEQ),CCH(4,ISERIES,NSEQ),
c     &              CCH(5,ISERIES,NSEQ),CCH(6,ISERIES,NSEQ),
c     &              CCH(7,ISERIES,NSEQ),ISERIES=1,MAXSERIES
	      ENDIF
1600      CONTINUE
1500    CONTINUE
1510      FORMAT(A10,',',F11.0,',',2(F8.0,','),
     &           5(2(I10,', '),2(F10.2,', '),8(F12.2,', ')))	! Output format for PC
400   CONTINUE			! ICELL loop
C
C Output results on a 20km2 grid cell
C
      WRITE(GRIDCHAN,20)KM20GRIDID,EAST(1),NORTH(1),
     &               (((C20KM2LU1TOLU2(IDEC,LU1,LU2),LU1=1,MAXLU1),
     &                               LU2=1,MAXLU1),IDEC=1,MAXDEC1),
     &               (((LU1TOLU220KM2(IDEC,LU1,LU2),LU1=1,MAXLU1),
     &                               LU2=1,MAXLU1),IDEC=1,MAXDEC1),
     &               (CCHGE(IDEC),IDEC=1,MAXDEC1),
     &               (OCCHGE(IDEC),IDEC=1,MAXDEC1),
     &               (CEQ(IDEC),IDEC=1,MAXDEC1),
     &               (OCEQ(IDEC),IDEC=1,MAXDEC1),
     &                GRIDERROR

      WRITE(OGRIDCHAN,20)KM20GRIDID,EAST(1),NORTH(1),
     &               (((OC20KM2LU1TOLU2(IDEC,LU1,LU2),LU1=1,MAXLU1),
     &                               LU2=1,MAXLU1),IDEC=1,MAXDEC1),
     &               (((OLU1TOLU220KM2(IDEC,LU1,LU2),LU1=1,MAXLU1),
     &                               LU2=1,MAXLU1),IDEC=1,MAXDEC1),
     &               (CCHGE(IDEC),IDEC=1,MAXDEC1),
     &               (OCCHGE(IDEC),IDEC=1,MAXDEC1),
     &               (CEQ(IDEC),IDEC=1,MAXDEC1),
     &               (OCEQ(IDEC),IDEC=1,MAXDEC1),
     &                GRIDERROR
20    FORMAT(A10,',',2(F8.0,','),
     &       6(6(9(F10.2,', '))),
     &       6(6(9(F12.2,', '))),
     &       4(9(F12.2,', ')),A50)		
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE PUT_GIS_FUTRES(SQUID,EAST,NORTH,NUMCELLSIN20KM2,
     &                       KM20GRIDID,DOMSOIL,PERSOIL,SOILID,STYPES,
     &                       WETCLASS,TOTSEQ,LUSEQ,SEQTYPE,FRACLU50,
     &                       CCHANGE10,CH4C10,CO2C10,N2ON10,
     &                       TRAPERR,DOMSOILC,DOMSOILBD,SOMDEPTH)	 
C
C Subroutine to record results of spatial simulation
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
      INTEGER MAXDEC				! Max.no.of decades included in calculation
	PARAMETER(MAXDEC=9)		
	INTEGER MAXDEC1				! Max.no.of decades included in calculation
	DATA MAXDEC1 /9/
	INTEGER MAXGROW				! Max.no.of growing seasons allowed
	PARAMETER (MAXGROW=300)
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXLU1				! Max.no.of land use types
	DATA MAXLU1 /6/
	INTEGER MAXSEQ				! Max.no of sequences included in GIS
	PARAMETER (MAXSEQ=MAXLU*MAXLU)						
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXSOMLAY			! Max.no.of soil layers entered			
	PARAMETER (MAXSOMLAY=10)
	INTEGER MAXWC				! Max.no of wetness classes
      PARAMETER (MAXWC=6)

      REAL CCLU1(MAXDEC,MAXLU)	! Carbon change for LU1 decade 1 to 6 in ktC/km2/10y
      CHARACTER*50 CELLERROR		! Error Message for 1km2 cell
      INTEGER EQUIERR
	DATA EQUIERR /6/			!		2=Error in equilibrium run
	REAL FRACCELL				! Fraction of cell in soil ISERIES and 
      INTEGER FRESCHAN			! Channel for results of future met.data inputs
	DATA FRESCHAN /66/
      CHARACTER*50 GRIDERROR		! Error Message for 20km2 grid
	INTEGER ICELL				! Local cell counter
      INTEGER IDEC				! Decade counter
	INTEGER ISERIES				! Local counter for soil series
	INTEGER ISORGANIC			! Marker for org.soil (0=mineral,1=organic)
	INTEGER IWC					! Local counter for wetness class
	INTEGER LUCERROR(MAXDEC,MAXLU+2,MAXLU+2) ! Error found in this LU change
	INTEGER LU1					! Land use 1
	INTEGER LU2					! Land use 2
	INTEGER METERR
	DATA METERR /2/				!		2=Error in weather data
	INTEGER NPPERR
	DATA NPPERR /1/				!		1=Error in NPP data
	INTEGER NSEQ				! Counter for number of land use sequences 
      REAL OCCLU1(MAXDEC,MAXLU)	! Carbon change for LU1 decade 1 to 6 from organic soils in ktC/km2/10y
      INTEGER OFRESCHAN			! Channel for results of future met.data inputs organic soils
	DATA OFRESCHAN /67/
	REAL PC						! % C in the top soil layer
	INTEGER SEQ1				! IN:Land use sequence of LU1->LU1
	REAL TEMP					! Temporary real for calculation

C
C Variables passed to/from calling subroutine
C
	REAL CCHANGE10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ) ! IN:Change 
								! in C content over the decade (kgC/ha/decade)
	REAL CH4C10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! IN:CH4-C 
								! emitted over the decade (kgC/ha/decade)
	REAL CO2C10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! IN:CO2-C 
								! emitted over the decade (kgC/ha/decade)
      INTEGER DOMSOIL(CELLSTOGRID,MAXSERIES)	! IN:5 major soil series in the  		
								! 1km2 cell x no. of cells in 20km2 grid
	REAL DOMSOILBD(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! Soil BD in 5 
												! major soil series under different LU 
												! in SOM layers (g/cm3)
	REAL DOMSOILC(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! Soil C in 5 
												! major soil series under different LU 
												! in SOM layers (kgC/ha)
	REAL EAST(CELLSTOGRID)		! IN:National grid easting x 1000
	REAL FRACLU50(CELLSTOGRID,MAXLU+2)	! IN:Fraction of 1km cell in LU1 in 1950
	CHARACTER*10 KM20GRIDID		! IN:20km sqare identifier
      INTEGER LUSEQ(MAXGROW)		! IN:Land use sequence
	REAL N2ON10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! IN:N2O-N 
								! emitted over the decade (kgN/ha/decade)
	INTEGER NGROW				! IN:No.growing seasons simulated
	REAL NORTH(CELLSTOGRID)		! IN:National grid northing x 1000
	INTEGER NUMCELLSIN20KM2		! IN:Number of cells in 20km2 cell
      INTEGER ORDER				! IN:Order of soils calculations
	REAL PERSOIL(CELLSTOGRID,MAXSERIES)		! IN:Percentage of each major soil 
								!	type in the 1km2 square cell
	INTEGER SEQTYPE				! IN:Type of sequence to be used
      INTEGER SOILID(CELLSTOGRID*MAXSERIES/2)	! IN:Soil series integer codes
	REAL SOMDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN:depth of soil organic matter layers(cm)
	REAL SQUID(CELLSTOGRID)		! IN:1km square identifier
	INTEGER STYPES				! IN:Number of soil types in 20km2 grid
	INTEGER TOTSEQ				! IN:Total number of sequences
      INTEGER TRAPERR(CELLSTOGRID*MAXSERIES/2,MAXSEQ) 
								! IN:Error trapping
								!		0=No error
								!		1=Error in NPP data
								!		2=Error in weather data
								!		3=Error in LU data
								!		4=Error in LU change data
								!		5=Soil error
								!		6=Equilibrium not found
      INTEGER WETCLASS(CELLSTOGRID,MAXSERIES) !IN:Wetness class for 
								! each major soil type in the 1km2 square cell 
								! x 20km2 cell
C
C Calculate overall changes for the cell
C
      DO 400 ICELL=1,NUMCELLSIN20KM2
        CELLERROR='         '
	  DO 500 IDEC=1,MAXDEC1
	    DO 550 LU1=1,MAXLU1
	      CCLU1(IDEC,LU1)=0
	      OCCLU1(IDEC,LU1)=0
550       CONTINUE	! LU1 loop
C
C For each soils series that exists in this cell...
C 
          DO 600 ISERIES=1,MAXSERIES
	      IF(PERSOIL(ICELL,ISERIES).GT.0)THEN
C
C ...get order of soils calculations
C
              CALL GET_SOIL_ORDER(DOMSOIL(ICELL,ISERIES),SOILID,
     &                            STYPES,ORDER)
C
C If this soil series is found...
C
              IF(ORDER.GT.0)THEN
C
C ...for each sequence in cell
C
	          DO 700 NSEQ=1,TOTSEQ
	            IF(TRAPERR(ORDER,NSEQ).EQ.EQUIERR)THEN
                    CELLERROR='Not initialised to equilibrium'
                    GRIDERROR='Not initialised to equilibrium'
			    ENDIF		! IF(TRAPERR(ORDER,NSEQ).EQ.EQUIERR)THEN
	            CALL GET_LUSEQ(LUSEQ,NGROW,NSEQ,SEQTYPE)
	            LU1=LUSEQ(1)
	            LU2=LUSEQ(6)
	            SEQ1=((LU1-1)*MAXLU1)+LU1
C
C ...check top soil layer for an organic soil
C
                  ISORGANIC=0
                  PC=DOMSOILC(ORDER,LU1,1)
				PC=PC/DOMSOILBD(ORDER,LU1,1)
				PC=PC/(SOMDEPTH(ORDER,LU1,1)*1000)
                  IF(PC.GT.6)ISORGANIC=1
C
C Record results for no land use change sequence only
C
                  IF(LU1.EQ.LU2)THEN
C 
C ...if there IS an error in results, record the error for this decade and LU change
C
	              IF(TRAPERR(ORDER,NSEQ).EQ.NPPERR.OR.
     &                   TRAPERR(ORDER,NSEQ).EQ.METERR)THEN
                      LUCERROR(IDEC,LU1,LU2)=1
	                IF(TRAPERR(ORDER,NSEQ).EQ.NPPERR)THEN
                        CELLERROR='Error in NPP'
	                ENDIF				! IF(TRAPERR(ORDER,NSEQ).EQ.NPPERR)THEN
	                IF(TRAPERR(ORDER,NSEQ).EQ.METERR)THEN
                        CELLERROR='Error in met.data'
	                ENDIF				! IF(TRAPERR(ORDER,NSEQ).EQ.METERR)THEN
C
C ...if there is no met or NPP error in results on this series and this sequence
C
	              ELSE				! ELSE to IF(TRAPERR(ORDER,NSEQ).EQ.NPPERR.OR.
									! TRAPERR(ORDER,NSEQ).EQ.METERR)THEN

C
C ...get fraction of cell on this soil under this sequence and this LU change
C and recalculate in t per 1km grid square
C kg/ha -> kt/km2: kg/(1000*1000)->t; ha x 100 ->km2; kg/ha x 100/1000000 -> kt/km2 =1/10000
C
	                TEMP=FRACLU50(ICELL,LU1)
	                TEMP=TEMP*(PERSOIL(ICELL,ISERIES)/100)
				    FRACCELL=TEMP
                      IWC=WETCLASS(ICELL,ISERIES)
                      TEMP=CCHANGE10(IDEC,ORDER,IWC,NSEQ)*FRACCELL/10000
                      CCLU1(IDEC,LU1)=CCLU1(IDEC,LU1)+TEMP
	                IF(ISORGANIC.EQ.1)THEN
                        OCCLU1(IDEC,LU1)=OCCLU1(IDEC,LU1)+TEMP
	                ENDIF
	              ENDIF
	            ENDIF
700		      CONTINUE			! NSEQ loop
              ENDIF				! IF(ORDER.GT.0)THEN
            ENDIF					! IF(PERSOIL(ICELL,ISERIES).GT.0)THEN
600       CONTINUE				! ISERIES loop
500     CONTINUE					! IDEC loop
C
C Output results
C      
        WRITE(FRESCHAN,10)KM20GRIDID,SQUID(ICELL),
     &                    EAST(ICELL),NORTH(ICELL),
     &                   ((CCLU1(IDEC,LU1),LU1=1,MAXLU1),IDEC=1,MAXDEC1)
     &                    ,CELLERROR
        WRITE(OFRESCHAN,10)KM20GRIDID,SQUID(ICELL),
     &                    EAST(ICELL),NORTH(ICELL),
     &                  ((OCCLU1(IDEC,LU1),LU1=1,MAXLU1),IDEC=1,MAXDEC1)
     &                    ,CELLERROR
10      FORMAT(A10,',',F11.0,',',2(F8.0,','),6(9(F10.2,', ')),A50)
400   CONTINUE
	END

C
C-------------------------------------------------------------
C
      SUBROUTINE SETERROR(FYMFERT,FYMFERT15,TFYMC,VOLAT,VOLAT15,					  
     &                      SOILN,SOIL15,AMMN,AMMN15,			  
     &                      DPMCARB0,DPMNIT0,DPMNLAB0,		  
     &                      RPMCARB0,RPMNIT0,RPMNLAB0,		  
     &                      HCARB0,HNIT0,HNLAB0,TIN,TAM,
     &                      CCHANGE10,CO2C10,CH4C10,N2ON10,
     &                      ISERIES,NSEQ,IWC)
C
C Subroutine to get filenames and simulation inputs
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID					! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
      INTEGER MAXDEC				! Max.no.of decades included in calculation
	PARAMETER(MAXDEC=9)		
      INTEGER MAXDEC1				! Max.no.of decades included in calculation
	DATA MAXDEC1 /9/		
	INTEGER MAXLAYER			! No.of layers in the soil profile
      PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXSEQ				! Max.no of sequences included in GIS
	PARAMETER (MAXSEQ=MAXLU*MAXLU)						
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXWC				! Max.no of wetness classes
      PARAMETER (MAXWC=6)

      INTEGER IDEC				! Local decade counter
	INTEGER IL					! Local layer counter
C
C Variables passed from calling subroutine
C
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
	REAL BCARB0(MAXLAYER)		! IN:C in soil biomass (kgC/ha/layer)
      REAL BNIT0(MAXLAYER)		! IN:N in soil biomass (kgN/ha/layer)
	REAL BNLAB0(MAXLAYER)		! IN:N15 in soil humus (kgN15/ha/layer)
	REAL CCHANGE10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ) ! IN/OUT:
								! Change in C content over the decade (kgC/ha/decade)
	REAL CH4C10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! IN/OUT:CH4-C
								!  emitted over the decade (kgC/ha/decade)
	REAL CO2C10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! IN/OUT:CO2-C 
								! emitted over the decade (kgC/ha/decade)
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
	REAL DPMNIT0(MAXLAYER)		! IN:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! IN:N15 in decomposable PM (kgN15/ha/layer)
	REAL FYMFERT				! IN:N input by fertiliser & FYM (kgN/ha)
	REAL FYMFERT15				! IN:N15 input by fertiliser & FYM (kgN15/ha)
	REAL HCARB0(MAXLAYER)		! IN:C in soil humus (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)		! IN:N in soil humus (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! IN:N15 in soil biomass (kgN15/ha/layer)
	INTEGER ISERIES				! IN:Counter for soil series
	INTEGER IWC					! IN:Wetness class
	REAL N2ON10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! IN/OUT:N2O-N 
								! emitted over the decade (kgN/ha/decade)
	INTEGER NSEQ				! Counter for number of land use sequences 
	REAL TAM(MAXLAYER)			! IN:15N/N in the ammonium in the layer
	REAL TIN(MAXLAYER)			! IN:15N/N in the nitrate in the layer
	REAL TFYMC					! IN:Total FYM C input (kgC/ha) DONT PASS?
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
	REAL RPMNIT0(MAXLAYER)		! IN:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! IN:N15 in resistant PM (kgN15/ha/layer)
	INTEGER SEQTYPE				! IN:Type of sequence to be used
	REAL SOILN(MAXLAYER)		! IN:Soil nitrate-N (kgN/ha/layer)
	REAL SOIL15(MAXLAYER)		! IN:Soil nitrate-N15 (kgN15/ha/layer)
	REAL VOLAT					! IN:N lost by volatilisation (kgN/ha)
	REAL VOLAT15				! IN:N15 lost by volatilisation (kgN15/ha)
C
C Set results to -999 to indicate error
C
      FYMFERT=-999
	FYMFERT15=-999
	TFYMC=-999
	VOLAT=-999
	VOLAT15=-999
	DO 100 IL=1,MAXLAYER1
	  SOILN(IL)=-999
	  SOIL15(IL)=-999
	  AMMN(IL)=-999
	  AMMN15(IL)=-999
	  BCARB0(IL)=-999
	  HCARB0(IL)=-999
	  DPMCARB0(IL)=-999
	  RPMCARB0(IL)=-999
        BNIT0(IL)=-999
	  HNIT0(IL)=-999
	  DPMNIT0(IL)=-999
	  RPMNIT0(IL)=-999
        HNLAB0(IL)=-999
	  BNLAB0(IL)=-999
	  DPMNLAB0(IL)=-999
	  RPMNLAB0(IL)=-999
	  TAM(IL)=-999
	  TIN(IL)=-999
100   CONTINUE
      DO 200 IDEC=1,MAXDEC1
        CCHANGE10(IDEC,ISERIES,IWC,NSEQ)=0
	  CH4C10(IDEC,ISERIES,IWC,NSEQ)=0
	  CO2C10(IDEC,ISERIES,IWC,NSEQ)=0
	  N2ON10(IDEC,ISERIES,IWC,NSEQ)=0
200   CONTINUE
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE SETNOERROR(TRAPERR)
C
C Subroutine to set the initial errors in the cell to none
C
	IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXLU1				! Max.no.of land use types
	DATA MAXLU1 /6/
	INTEGER MAXSEQ				! Max.no of sequences included in GIS
	PARAMETER (MAXSEQ=MAXLU*MAXLU)						
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						

	INTEGER ISERIES				! Soil series counter
	INTEGER NSEQ				! Sequence counter
C
C Variables passed to/from this subroutine
C
      INTEGER TRAPERR(CELLSTOGRID*MAXSERIES/2,MAXSEQ) 
								! IN:Error trapping
								!		0=No error
								!		1=Error in NPP data
								!		2=Error in weather data
								!		3=Error in LU data
								!		4=Error in LU change data
								!		5=Soil error
								!		6=Equilibrium not found
C
C Set TRAPERR to 0 (no error)
C
      DO 100 ISERIES=1,CELLSTOGRID*MAXSERIES/2
	  DO 200 NSEQ=1,MAXSEQ
	    TRAPERR(ISERIES,NSEQ)=0
200     CONTINUE
100   CONTINUE
C
C Leave SETNOERROR
C
      END      
C
C-------------------------------------------------------------
C
C
C------------------------------------------------------------
C
      SUBROUTINE PUT_START_DATA(ICELL,GRIDID,SQUID,EAST,NORTH,AREA,
     &                          SOMDEPTH,DOMSOIL,SOILID,STYPES,
     &                          DOMSOILC,DOMSOILCLAY,DOMSOILSILT,				
     &						  DOMSOILSAND,DOMSOILBD,DOMSOILPH,
     &						  DOMSOILISIMP,DOMSOILIMPDEPTH,
     &                          PERSOIL,FRACLU50,LU1TOLU2,NSOMLAY)
C
C Write out start soils and land use data
C
      IMPLICIT NONE
C
C Parameters
C
      INTEGER MAXDEC											! Max.no.of decades included in calculation
	PARAMETER(MAXDEC=7)		
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER MAXLU											! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXSERIES										! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXSOMLAY										! Max.no.of soil layers entered			
	PARAMETER (MAXSOMLAY=10)
C
C Variables local to this subroutine
C
	REAL BULKDENS(MAXLAYER)											! Soil bulk density (g cm-3)
	INTEGER CELLSTOGRID										! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
      REAL CLAY(MAXLAYER)										! Clay content of the layer (%)
	REAL FTEMP												! Temporary real number
	INTEGER IL												! Local layer counter
	INTEGER ILU												! Local land use counter
	REAL IOM(MAXLAYER)										! Inert organic C (kgC/ha/layer)
	INTEGER ISERIES											! Local counter for soil series
      INTEGER ISERROR											! Code for error in soil 1=error 0=noerror
	INTEGER LAYERMISSING									! Inidicates whether C data for layer is missing
	REAL SAND(MAXLAYER)										! Sand content of the layer (%)
	REAL SILT(MAXLAYER)										! Silt content of the layer (%)
	REAL SOILCTO100											! Soil carbon to 100 cm
	REAL SOILPH(MAXLAYER)									! pH of soil in this layer
	INTEGER STARTCHAN										! Channel for outputting start data
	REAL TOC(MAXLAYER)										! Total organic C (kgC/ha/layer)
	REAL WIDTH												! Width of layer (cm)
C
C Variables passed to/from calling subroutine
C
	REAL AREA(CELLSTOGRID)												! IN(CALL):Area of cell accounted for (m2)
      INTEGER DOMSOIL(CELLSTOGRID,MAXSERIES)					! IN:5 major soil series in the  		
															! 1km2 cell x no. of cells in 20km2 grid
	REAL DOMSOILBD(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! Soil BD in 5 
												! major soil series under different LU 
												! in SOM layers (g/cm3)
	REAL DOMSOILC(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN(CALL): Soil C in 5 
															!   major soil series under different LU 
															! in SOM layers (kgC/ha)
	REAL DOMSOILCLAY(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! % Soil clay in 5 
															! major soil series under different LU 
	REAL DOMSOILIMPDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN: Depth of impermeable layer if present (cm)
	INTEGER DOMSOILISIMP(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN: Flag for impermeable layer (0=No; 1=Yes)
	REAL DOMSOILPH(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! OUT:Soil PH in 5 
												! major soil series under different LU 
												! in SOM layers (g/cm3)
	REAL DOMSOILSAND(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! % Soil sand in 5 
												! major soil series under different LU 
	REAL DOMSOILSILT(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! % Soil silt in 5 
												! major soil series under different LU 
	REAL EAST(CELLSTOGRID)									! IN(CALL): National grid easting x 1000
	REAL FRACLU50(CELLSTOGRID,MAXLU+2)						! IN(CALL): Fraction of 1km cell in LU1 in dec.
      CHARACTER*10 GRIDID										! IN(CALL): 20km square identifier of last 1km cell
      INTEGER ICELL											! IN(CALL): Cell number
	REAL LU1TOLU2(CELLSTOGRID,MAXDEC,MAXLU+2,MAXLU+2)		! IN(CALL): Fracn of LU1 changed 
															!   to LU2 in dec
	REAL NORTH(CELLSTOGRID)									! IN(CALL): National grid northing x 1000
	INTEGER NSOMLAY(CELLSTOGRID*MAXSERIES/2,MAXLU)			! IN:Number of SOM layers
      INTEGER ORDER(MAXSERIES)								! IN(GET_SOIL_ORDER):Order of soils calculations
	REAL PERSOIL(CELLSTOGRID,MAXSERIES)						! IN(CALL): % of each major soil type in the
															!	1km square cell
      INTEGER SOILID(CELLSTOGRID*MAXSERIES/2)					! IN(CALL):Soil series integer codes
	REAL SOMDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN(CALL):depth of soil organic matter layers
	REAL SQUID(CELLSTOGRID)									! IN(CALL): 1km square identifier
	INTEGER STYPES											! IN(CALL):Number of soil types in 20km2 grid
C
C Output channel
C
      STARTCHAN=68
C
C Get order of soils calculations
C
      DO 100 ISERIES=1,MAXSERIES
        CALL GET_SOIL_ORDER(DOMSOIL(ICELL,ISERIES),SOILID,
     &                      STYPES,ORDER(ISERIES))
100   CONTINUE

C
C Write out starting conditions
C
c     WRITE(STARTCHAN,10)GRIDID,SQUID(ICELL),EAST(ICELL),NORTH(ICELL),
c    &                   AREA,(((SOMDEPTH(ORDER(ISERIES),ILU,IL),
c    &                          DOMSOILC(ORDER(ISERIES),ILU,IL),
c    &                          IL=1,MAXSOMLAY),
c    &                          ILU=1,MAXLU),
c    &                          ISERIES=1,MAXSERIES),
c    &                  (PERSOIL(ICELL,ISERIES),ISERIES=1,MAXSERIES),
c    &                  (ILU,FRACLU50(ICELL,ILU),ILU=1,MAXLU)
c10    FORMAT(A10,2X,I5,2X,3(F10.0,2X),
c    &      600(F10.0,2X),
c    &       5(F10.0,2X),I5,2X,6(F10.0,2X))
C
C Calculate soil C to 100 cm in each 1km2 grid
C
      SOILCTO100=0
      DO 200 ISERIES=1,MAXSERIES
        IF(ORDER(ISERIES).LT.1.OR.
     &     ORDER(ISERIES).GT.CELLSTOGRID*MAXSERIES/2)
     &    GOTO 303
	  DO 300 ILU=1,MAXLU
C
C Get soil data
C
c	    CALL GETSOIL(ISERIES,ILU,SOMDEPTH,NSOMLAY,							
c     &                 DOMSOILC,DOMSOILCLAY,DOMSOILSILT,				
c     &                 DOMSOILSAND,DOMSOILBD,DOMSOILPH,
c     &                 DOMSOILISIMP,DOMSOILIMPDEPTH,
c     &                 TOC,IOM,CLAY,SILT,SAND,BULKDENS,SOILPH,ISERROR)
C
C Add up the total C to 100cm
C
c          IF(ISERROR.EQ.0)THEN
c            DO 400 IL=1,MAXLAYER1
c              IF(IL*(MAXLAYER1/MAXDEPTH).LE.100)THEN
c	          FTEMP=FRACLU50(ICELL,ILU)*(PERSOIL(ICELL,ISERIES)/100)
c		      SOILCTO100=SOILCTO100+(TOC(IL)*FTEMP)
c	        ENDIF
c400         CONTINUE
c          ENDIF

c
c Only output initial soil C if all soil layers have available C data
c (i.e. are not 999)
	    LAYERMISSING = 0
		DO 500 IL=1,MAXSOMLAY
	      IF(DOMSOILC(ORDER(ISERIES),ILU,IL).GE.9990000)THEN
              LAYERMISSING=1
	      ENDIF
500		CONTINUE
	    IF(LAYERMISSING.EQ.0)THEN
			DO 400 IL=1,MAXSOMLAY
			  IF(FRACLU50(ICELL,ILU).GT.0.AND.
     &              FRACLU50(ICELL,ILU).LE.1)THEN
 				IF(SOMDEPTH(ORDER(ISERIES),ILU,IL).LT.100)THEN
 				  SOILCTO100=SOILCTO100+(FRACLU50(ICELL,ILU)
     &                     *DOMSOILC(ORDER(ISERIES),ILU,IL))
     &                     *(PERSOIL(ICELL,ISERIES)/100)
 				ELSEIF(SOMDEPTH(ORDER(ISERIES),ILU,IL).GE.100)THEN
 				  WIDTH=SOMDEPTH(ORDER(ISERIES),ILU,IL)
 				  IF(IL.GT.1)THEN
 					WIDTH=WIDTH-SOMDEPTH(ORDER(ISERIES),ILU,IL-1)
 				  ENDIF
 				  SOILCTO100=SOILCTO100+(FRACLU50(ICELL,ILU)
     &                     *DOMSOILC(ORDER(ISERIES),ILU,IL)
     &                     *(PERSOIL(ICELL,ISERIES)/100)
     &                     *(100-SOMDEPTH(ORDER(ISERIES),ILU,IL-1))
     &                     /WIDTH)
 				  GOTO 404
 				ENDIF
 			  ENDIF
400	        CONTINUE
		  ENDIF
404       CONTINUE
300     CONTINUE
303     CONTINUE
200   CONTINUE

      WRITE(STARTCHAN,10)GRIDID,SQUID(ICELL),EAST(ICELL),NORTH(ICELL),
     &                   AREA,SOILCTO100
10    FORMAT(A10,2X,I10,2X,4(F10.0,2X))
      END
C
C------------------------------------------------------------
C
      SUBROUTINE SET_GIS_FERT(FERTADD,FERTADD15)
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
      PARAMETER (MAXFERTSTEPS=30*60*24*7*3)	! Timesteps = 3 wks 
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
C--------------------------------------------------------
C
	SUBROUTINE TEST1_OPENCHAN(WATERMODEL,FULL_OUTPUT,SUMMARY_OUTPUT)
C
C Subroutine to open channels
C
      IMPLICIT NONE
C
C Arguments
C
	INTEGER WATERMODEL			! IN(READ_MODEL):Water model
	LOGICAL, INTENT(IN) :: FULL_OUTPUT     ! Flags whether to write detailed output files
	LOGICAL, INTENT(IN) :: SUMMARY_OUTPUT  ! Flags whether to write a summary output file
C
C Variables local to this subroutine
C
      INTEGER CHAN_TOTC			! Channel used to output TOTC and IOM
      INTEGER CHAN_BIOC			! Channel used to output BIOC
      INTEGER CHAN_HUMC			! Channel used to output HUMC
      INTEGER CHAN_RPMC			! Channel used to output RPMC
      INTEGER CHAN_DPMC			! Channel used to output DPMC
      INTEGER CHAN_TOTN			! Channel used to output TOTN 
      INTEGER CHAN_BION			! Channel used to output BION
      INTEGER CHAN_HUMN			! Channel used to output HUMN
      INTEGER CHAN_RPMN			! Channel used to output RPMN
      INTEGER CHAN_DPMN			! Channel used to output DPMN
      INTEGER CHAN_NH4N			! Channel used to output NH4N
      INTEGER CHAN_NO3N			! Channel used to output NO3N
      INTEGER CHAN_DENITN			! Channel used to output DENITN
      INTEGER CHAN_MINERN			! Channel used to output MINERN
      INTEGER CHAN_NITRIFN		! Channel used to output NITRIFN
      INTEGER CHAN_CROPN			! Channel used to output CROPN
      INTEGER CHAN_LEACHN			! Channel used to output LEACHN
      INTEGER CHAN_SENESN			! Channel used to output SENESN
      INTEGER CHAN_VOLN			! Channel used to output VOLN
      INTEGER CHAN_SOILW			! Channel used to outut SOILW
	INTEGER CHAN_SOILWC			! Channel used to outut SOILWC
	INTEGER CHAN_DOC			! Channel used to output DOC
	INTEGER CHAN_DON			! Channel used to output DON
	INTEGER CHAN_CH4			! Channel used to output CH4 from the Aitkenhead model
	INTEGER CHAN_CH4_RICHARDS	! Channel used to output CH4 from the Richards model
	INTEGER CHAN_CO2			! Channel used to output CO2
	INTEGER CHAN_PLANT_PARAM_SWAT	!Channel used to output PLANT_PARAM_SWAT
	INTEGER CHAN_EVAP_SWAT		! Channel used to output EVAP_SWAT
	INTEGER CHAN_BYRAIN			! Channel used to output BYRAIN
	INTEGER CHAN_EVAP_SUNDIAL	! Channel used to output EVAP_SUNDIAL
      INTEGER CHAN_SUMMARY      ! Channel used for summary output
	INTEGER SUNDIAL_WATER		!   SUNDIAL water model
	 INTEGER SWAT_WATER			!	SWAT water model
	 DATA SUNDIAL_WATER,SWAT_WATER /0,1/

      PARAMETER (CHAN_TOTN=124,CHAN_BION=125,CHAN_HUMN=126,
     &           CHAN_RPMN=127,CHAN_DPMN=128,CHAN_NH4N=129,
     &           CHAN_NO3N=130,CHAN_DENITN=131,CHAN_MINERN=132,
     &           CHAN_NITRIFN=133,CHAN_CROPN=134,CHAN_LEACHN=135,
     &           CHAN_SENESN=136,CHAN_VOLN=137,
     &           CHAN_CO2=138,CHAN_CH4=139, CHAN_CH4_RICHARDS=140)
      PARAMETER (CHAN_TOTC=24,CHAN_BIOC=25,CHAN_HUMC=26,CHAN_RPMC=27,
     &           CHAN_DPMC=28,CHAN_SOILW=20,CHAN_DOC=23,CHAN_DON=29)
	PARAMETER (CHAN_SOILWC=30,CHAN_PLANT_PARAM_SWAT=31,
     &			CHAN_EVAP_SWAT=32,CHAN_BYRAIN=33,
     &			CHAN_EVAP_SUNDIAL=34,CHAN_SUMMARY=37)
C
C Open Channels
C
      IF (SUMMARY_OUTPUT) THEN
        OPEN(CHAN_SUMMARY, FILE='SUMMARY.OUT',STATUS='UNKNOWN')
      ENDIF
      IF (FULL_OUTPUT) THEN
        OPEN(57, FILE='BALANCE_C.OUT',STATUS='UNKNOWN')
        OPEN(58, FILE='BALANCE_N.OUT',STATUS='UNKNOWN')
        OPEN(21, FILE='SOILN.OUT',  STATUS='UNKNOWN')
        OPEN(CHAN_SOILW, FILE='SOILW.OUT',  STATUS='UNKNOWN')
        OPEN(CHAN_DOC,  FILE='DOC.OUT', STATUS='UNKNOWN')
        OPEN(CHAN_DON,  FILE='DON.OUT', STATUS='UNKNOWN')
        OPEN(CHAN_TOTC, FILE='TOTC.OUT',STATUS='UNKNOWN')
        OPEN(CHAN_BIOC, FILE='BIOC.OUT',STATUS='UNKNOWN')
        OPEN(CHAN_HUMC, FILE='HUMC.OUT',STATUS='UNKNOWN')
        OPEN(CHAN_RPMC, FILE='RPMC.OUT',STATUS='UNKNOWN')
        OPEN(CHAN_DPMC, FILE='DPMC.OUT',STATUS='UNKNOWN')
        OPEN(CHAN_TOTN, FILE='TOTN.OUT',STATUS='UNKNOWN')
        OPEN(CHAN_BION, FILE='BION.OUT',STATUS='UNKNOWN')
        OPEN(CHAN_HUMN, FILE='HUMN.OUT',STATUS='UNKNOWN')
        OPEN(CHAN_DPMN, FILE='DPMN.OUT',STATUS='UNKNOWN')
        OPEN(CHAN_RPMN, FILE='RPMN.OUT',STATUS='UNKNOWN')
        OPEN(CHAN_NH4N, FILE='NH4N.OUT',STATUS='UNKNOWN')
        OPEN(CHAN_NO3N, FILE='NO3N.OUT',STATUS='UNKNOWN')
        OPEN(CHAN_DENITN, FILE='DENITN.OUT',STATUS='UNKNOWN')
        OPEN(CHAN_MINERN, FILE='MINERN.OUT',STATUS='UNKNOWN')
        OPEN(CHAN_NITRIFN, FILE='NITRIFN.OUT',STATUS='UNKNOWN')
        OPEN(CHAN_CROPN, FILE='CROPN.OUT',STATUS='UNKNOWN')
        OPEN(CHAN_LEACHN, FILE='LEACHN.OUT',STATUS='UNKNOWN')
C	OPEN(CHAN_SENESN, FILE='SENESN.OUT',STATUS='UNKNOWN')
        OPEN(CHAN_VOLN, FILE='VOLN.OUT',STATUS='UNKNOWN')
        OPEN(CHAN_CO2, FILE='CO2.OUT',STATUS='UNKNOWN')
        OPEN(CHAN_CH4, FILE='CH4.OUT',STATUS='UNKNOWN')
        OPEN(CHAN_CH4_RICHARDS, FILE='CH4_RICHARDS.OUT',
     &                STATUS='UNKNOWN')
        OPEN(CHAN_BYRAIN, FILE='BYRAIN.OUT',STATUS='UNKNOWN',
     &                ACTION='WRITE')
        OPEN(CHAN_EVAP_SUNDIAL, FILE='EVAP_SUNDIAL.OUT',
     &				STATUS='UNKNOWN', ACTION='WRITE')
        IF(WATERMODEL.EQ.SWAT_WATER)THEN
            OPEN(CHAN_SOILWC,FILE='SOILWC.OUT',STATUS='UNKNOWN',
     &				ACTION='WRITE')
            OPEN(CHAN_PLANT_PARAM_SWAT, FILE='PLANT_PARAM_SWAT.OUT',
     &				STATUS='UNKNOWN', ACTION='WRITE')

            OPEN(CHAN_EVAP_SWAT, FILE='EVAP_SWAT.OUT',STATUS='UNKNOWN',
     &				ACTION='WRITE')
        ENDIF
      ENDIF
C
C Titles in results files
C
      IF (SUMMARY_OUTPUT) THEN
        WRITE(CHAN_SUMMARY,259)
259     FORMAT("All units are kg/ha except avail_water ",
     &         "which is in mm")
        WRITE(CHAN_SUMMARY,260)
260     FORMAT(" year step growing_season ",
     &         "dpm_c rpm_c bio_c hum_c iom_c total_soc total_soc_0-30 ",
     &         "total_soc_30-100 co2_c ch4_c ",
     &         "dpm_n rpm_n bio_n hum_n total_son nh4_n no3_n ",
     &         "denitrified_n nitrification_n  volatilised_n ",
     &         "leached_n net_mineralised_n n2o_n no_n doc don ",
     &         "avail_water")
      ENDIF
      IF (FULL_OUTPUT) THEN
      WRITE(21,210)
210   FORMAT("Season Step",      
     &       "       NO3       NO3       NO3       NO3",
     &       "       NO3       NO3       NO3       NO3",
     &       "       NO3       NO3       NO3       NO3",
     &       "       NO3       NO3       NO3       NO3",
     &       "       NO3       NO3       NO3       NO3",
     &       "       NO3       NO3       NO3       NO3",
     &       "       NO3       NO3       NO3       NO3",
     &       "       NO3       NO3",
     &       "       NO3       NO3       NO3       NO3",
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
     &       "       NH4       NH4",
     &       "       NH4       NH4       NH4       NH4",
     &       "       NH4       NH4       NH4       NH4",
     &       "       NH4       NH4       NH4       NH4",
     &       "       NH4       NH4       NH4       NH4",
     &       "       NH4       NH4       NH4       NH4",
     &       "       NH4       NH4       NH4       NH4",
     &       "       NH4       NH4       NH4       NH4",
     &       "       NH4       NH4"/
     &       "                  0-5cm     5-10cm   10-15cm   15-20cm",
     &                    "   20-25cm    25-30cm   30-35cm   35-40cm",
     &                    "   40-45cm    45-50cm   50-55cm   55-60cm",
     &                    "   60-65cm    65-70cm   70-75cm   75-80cm",
     &                    "   80-85cm    85-90cm   90-95cm  95-100cm",
     &                    " 100-105cm 105-110cm 110-115cm  115-120cm",
     &                    " 120-125cm 125-130cm 130-135cm  135-140cm",
     &                    " 140-145cm 145-150cm 150-155cm, 155-160cm",
     &                    " 160-165cm 165-170cm 170-175cm, 175-180cm",
     &                    " 180-185cm 185-180cm 190-195cm, 195-200cm",
     &                    " 200-205cm 205-210cm 210-215cm, 215-220cm",
     &                    " 220-225cm 225-230cm 230-245cm, 245-250cm",
     &                    " 250-255cm 255-260cm 260-265cm, 265-270cm",
     &                    " 270-275cm 275-280cm 280-285cm, 285-290cm",
     &                    " 290-295cm 295-300cm",
     &                    "     0-5cm     5-10cm   10-15cm   15-20cm",
     &                    "   20-25cm    25-30cm   30-35cm   35-40cm",
     &                    "   40-45cm    45-50cm   50-55cm   55-60cm",
     &                    "   60-65cm    65-70cm   70-75cm   75-80cm",
     &                    "   80-85cm    85-90cm   90-95cm  95-100cm",
     &                    " 100-105cm 105-110cm 110-115cm  115-120cm",
     &                    " 120-125cm 125-130cm 130-135cm  135-140cm",
     &                    " 140-145cm 145-150cm 150-155cm, 155-160cm",
     &                    " 160-165cm 165-170cm 170-175cm, 175-180cm",
     &                    " 180-185cm 185-180cm 190-195cm, 195-200cm",
     &                    " 200-205cm 205-210cm 210-215cm, 215-220cm",
     &                    " 220-225cm 225-230cm 230-245cm, 245-250cm",
     &                    " 250-255cm 255-260cm 260-265cm, 265-270cm",
     &                    " 270-275cm 275-280cm 280-285cm, 285-290cm",
     &                    " 290-295cm 295-300cm")
c      WRITE(57,570)
570   FORMAT("Results output in order"/
     &       "DOY"/
     &       "H-H"/
     &       "YEAR"/
     &       "DPM-C in 5cm layers to depth of measurements (kg/ha)"/
     &       "RPM-C in 5cm layers to depth of measurements (kg/ha)"/
     &       "BIO-C in 5cm layers to depth of measurements (kg/ha)"/
     &       "HUM-C in 5cm layers to depth of measurements (kg/ha)"/
     &       "CO2-C in 5cm layers to depth of measurements (kg/ha)"/
     &       "CH4-C in 5cm layers to depth of measurements (kg/ha)"/
     &       "Tot.DPM-C summed to depth of measurements (kg/ha)"/
     &       "Tot.RPM-C summed to depth of measurements (kg/ha)"/
     &       "Tot.BIO-C summed to depth of measurements (kg/ha)"/
     &       "Tot.HUM-C summed to depth of measurements (kg/ha)"/
     &       "Tot.CO2-C summed to depth of measurements (kg/ha)"/
     &       "Tot.CH4-C summed to depth of measurements (kg/ha)"/
     &       "Tot.IOM-C summed to depth of measurements (kg/ha)"/
     &       "Ann.CO2-C emitted from depth of measurements (kg/ha)"/
     &       "Organic inputs (kg/ha)"/
     &       "Summed CO2-C emitted from depth of measurements (kg/ha)"/)
	WRITE(57,571) 
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
     &       "Tot.CH4    CH4.ATM    Tot.IOM  Tot.Input  CO2_for_year ",
     &       "Total_C    Plant_Input_C"      )

5719   FORMAT("    DOY  H-H YEAR      ",
     &       "DPM        DPM        DPM        DPM        DPM        ",
     &       "DPM        DPM        DPM        DPM        DPM        ",
     &       "DPM        DPM        DPM        DPM        DPM        ",
     &       "DPM        DPM        DPM        DPM        DPM        ",
     &       "DPM        DPM        DPM        DPM        DPM        ",
     &       "DPM        DPM        DPM        DPM        DPM        ",
     &       "DPM        DPM        ",
     &       "RPM        RPM        RPM        RPM        RPM        ",
     &       "RPM        RPM        RPM        RPM        RPM        ",
     &       "RPM        RPM        RPM        RPM        RPM        ",
     &       "RPM        RPM        RPM        RPM        RPM        ",
     &       "RPM        RPM        RPM        RPM        RPM        ",
     &       "RPM        RPM        RPM        RPM        RPM        ",
     &       "RPM        RPM        ",
     &       "BIO        BIO        BIO        BIO        BIO        ",
     &       "BIO        BIO        BIO        BIO        BIO        ",
     &       "BIO        BIO        BIO        BIO        BIO        ",
     &       "BIO        BIO        BIO        BIO        BIO        ",
     &       "BIO        BIO        BIO        BIO        BIO        ",
     &       "BIO        BIO        BIO        BIO        BIO        ",
     &       "BIO        BIO        ",
     &       "HUM        HUM        HUM        HUM        HUM        ",
     &       "HUM        HUM        HUM        HUM        HUM        ",
     &       "HUM        HUM        HUM        HUM        HUM        ",
     &       "HUM        HUM        HUM        HUM        HUM        ",
     &       "HUM        HUM        HUM        HUM        HUM        ",
     &       "HUM        HUM        HUM        HUM        HUM        ",
     &       "HUM        HUM        ",
     &       "CO2        CO2        CO2        CO2        CO2        ",
     &       "CO2        CO2        CO2        CO2        CO2        ",
     &       "CO2        CO2        CO2        CO2        CO2        ",
     &       "CO2        CO2        CO2        CO2        CO2        ",
     &       "CO2        CO2        CO2        CO2        CO2        ",
     &       "CO2        CO2        CO2        CO2        CO2        ",
     &       "CO2        CO2        ",
     &       "CH4        CH4        CH4        CH4        CH4        ",
     &       "CH4        CH4        CH4        CH4        CH4        ",
     &       "CH4        CH4        CH4        CH4        CH4        ",
     &       "CH4        CH4        CH4        CH4        CH4        ",
     &       "CH4        CH4        CH4        CH4        CH4        ",
     &       "CH4        CH4        CH4        CH4        CH4        ",
     &       "CH4        CH4    ",
     &       "Tot.DPM    Tot.RPM    Tot.BIO    Tot.HUM    Tot.CO2    ",
     &       "Tot.CH4    Tot.IOM  Tot.Input Annual CO2 (in kgC/ha)")
c	WRITE(57,572)
5729	FORMAT("                     0-5cm     5-10cm    10-15cm    ",
     &       "15-20cm    20-25cm    25-30cm    30-35cm    35-40cm    ",
     &       "40-45cm    45-50cm    ",
     &     "50-55cm    55-60cm    60-65cm    65-70cm    70-75cm    ",
     &     "80-85cm    80-85cm    85-90cm    90-95cm   95-100cm  ",
     &     "100-105cm  105-110cm  110-115cm  115-120cm  120-125cm  ",
     &     "125-130cm  130-135cm  135-140cm  140-145cm  145-150cm  ",
     &     "150-155cm  155-160cm  ",
     &     "    0-5cm     5-10cm    10-15cm    15-20cm    20-25cm    ",
     &       "25-30cm    30-35cm    35-40cm    40-45cm    45-50cm    ",
     &     "50-55cm    55-60cm    60-65cm    65-70cm    70-75cm    ",
     &     "80-85cm    80-85cm    85-90cm    90-95cm   95-100cm  ",
     &     "100-105cm  105-110cm  110-115cm  115-120cm  120-125cm  ",
     &     "125-130cm  130-135cm  135-140cm  140-145cm  145-150cm  ",
     &     "150-155cm  155-160cm  ",
     &     "    0-5cm     5-10cm    10-15cm    15-20cm    20-25cm    ",
     &       "25-30cm    30-35cm    35-40cm    40-45cm    45-50cm    ",
     &     "50-55cm    55-60cm    60-65cm    65-70cm    70-75cm    ",
     &     "80-85cm    80-85cm    85-90cm    90-95cm   95-100cm  ",
     &     "100-105cm  105-110cm  110-115cm  115-120cm  120-125cm  ",
     &     "125-130cm  130-135cm  135-140cm  140-145cm  145-150cm  ",
     &     "150-155cm  155-160cm  ",
     &     "    0-5cm     5-10cm    10-15cm    15-20cm    20-25cm    ",
     &       "25-30cm    30-35cm    35-40cm    40-45cm    45-50cm    ",
     &     "50-55cm    55-60cm    60-65cm    65-70cm    70-75cm    ",
     &     "80-85cm    80-85cm    85-90cm    90-95cm   95-100cm  ",
     &     "100-105cm  105-110cm  110-115cm  115-120cm  120-125cm  ",
     &     "125-130cm  130-135cm  135-140cm  140-145cm  145-150cm  ",
     &     "150-155cm  155-160cm  ",
     &     "    0-5cm     5-10cm    10-15cm    15-20cm    20-25cm    ",
     &       "25-30cm    30-35cm    35-40cm    40-45cm    45-50cm    ",
     &     "50-55cm    55-60cm    60-65cm    65-70cm    70-75cm    ",
     &     "80-85cm    80-85cm    85-90cm    90-95cm   95-100cm  ",
     &     "100-105cm  105-110cm  110-115cm  115-120cm  120-125cm  ",
     &     "125-130cm  130-135cm  135-140cm  140-145cm  145-150cm  ",
     &     "150-155cm  155-160cm  ",
     &     "    0-5cm     5-10cm    10-15cm    15-20cm    20-25cm    ",
     &       "25-30cm    30-35cm    35-40cm    40-45cm    45-50cm    ",
     &     "50-55cm    55-60cm    60-65cm    65-70cm    70-75cm    ",
     &     "80-85cm    80-85cm    85-90cm    90-95cm   95-100cm  ",
     &     "100-105cm  105-110cm  110-115cm  115-120cm  120-125cm  ",
     &     "125-130cm  130-135cm  135-140cm  140-145cm  145-150cm  ",
     &     "150-155cm  155-160cm  ")
572	FORMAT("                     0-5cm     5-10cm    10-15cm    ",
     &       "15-20cm    20-25cm    25-30cm    30-35cm    35-40cm    ",
     &       "40-45cm    45-50cm    ",
     &     "50-55cm    55-60cm    60-65cm    65-70cm    70-75cm    ",
     &     "80-85cm    80-85cm    85-90cm    90-95cm   95-100cm  ",
     &     "100-105cm  105-110cm  110-115cm  115-120cm  120-125cm  ",
     &     "125-130cm  130-135cm  135-140cm  140-145cm  145-150cm  ",
     &     "150-155cm  155-160cm  ",
     &     "    0-5cm     5-10cm    10-15cm    15-20cm    20-25cm    ",
     &       "25-30cm    30-35cm    35-40cm    40-45cm    45-50cm    ",
     &     "50-55cm    55-60cm    60-65cm    65-70cm    70-75cm    ",
     &     "80-85cm    80-85cm    85-90cm    90-95cm   95-100cm  ",
     &     "100-105cm  105-110cm  110-115cm  115-120cm  120-125cm  ",
     &     "125-130cm  130-135cm  135-140cm  140-145cm  145-150cm  ",
     &     "150-155cm  155-160cm  ",
     &     "    0-5cm     5-10cm    10-15cm    15-20cm    20-25cm    ",
     &       "25-30cm    30-35cm    35-40cm    40-45cm    45-50cm    ",
     &     "50-55cm    55-60cm    60-65cm    65-70cm    70-75cm    ",
     &     "80-85cm    80-85cm    85-90cm    90-95cm   95-100cm  ",
     &     "100-105cm  105-110cm  110-115cm  115-120cm  120-125cm  ",
     &     "125-130cm  130-135cm  135-140cm  140-145cm  145-150cm  ",
     &     "150-155cm  155-160cm  ",
     &     "    0-5cm     5-10cm    10-15cm    15-20cm    20-25cm    ",
     &       "25-30cm    30-35cm    35-40cm    40-45cm    45-50cm    ",
     &     "50-55cm    55-60cm    60-65cm    65-70cm    70-75cm    ",
     &     "80-85cm    80-85cm    85-90cm    90-95cm   95-100cm  ",
     &     "100-105cm  105-110cm  110-115cm  115-120cm  120-125cm  ",
     &     "125-130cm  130-135cm  135-140cm  140-145cm  145-150cm  ",
     &     "150-155cm  155-160cm  ",
     &     "    0-5cm     5-10cm    10-15cm    15-20cm    20-25cm    ",
     &       "25-30cm    30-35cm    35-40cm    40-45cm    45-50cm    ",
     &     "50-55cm    55-60cm    60-65cm    65-70cm    70-75cm    ",
     &     "80-85cm    80-85cm    85-90cm    90-95cm   95-100cm  ",
     &     "100-105cm  105-110cm  110-115cm  115-120cm  120-125cm  ",
     &     "125-130cm  130-135cm  135-140cm  140-145cm  145-150cm  ",
     &     "150-155cm  155-160cm  ",
     &     "    0-5cm     5-10cm    10-15cm    15-20cm    20-25cm    ",
     &       "25-30cm    30-35cm    35-40cm    40-45cm    45-50cm    ",
     &     "50-55cm    55-60cm    60-65cm    65-70cm    70-75cm    ",
     &     "80-85cm    80-85cm    85-90cm    90-95cm   95-100cm  ",
     &     "100-105cm  105-110cm  110-115cm  115-120cm  120-125cm  ",
     &     "125-130cm  130-135cm  135-140cm  140-145cm  145-150cm  ",
     &     "150-155cm  155-160cm  ")
	WRITE(58,22)"Timestep H-H Year BioHum    PM      DPM     RPM  ",
     &            "   NO3     NH4     Fert    Stub    ATM     Seed   M",
     &            "iner   Denit  NitN2O   NitNO PNitN2O Nit15N2O Nit15",
     &            "NO PNit15N2O DenN2 DenN2O Den15N2 Den15N2O Vol   Cr",
     &            "op Leach Senes  DON    Tot.In Tot.Out  Diff        "
22    FORMAT(5A51)
C
C Title lines for SOILW, TOTC, HUMC, BIOC, DPMC, RPMC.OUT 
C      
	WRITE(CHAN_SOILW,240)'Soil water - values in mm'
	WRITE(CHAN_TOTC,240)'Total soil C - values in kgC/ha/5cm layer'
	WRITE(CHAN_HUMC,240)'Humus C - values in kgC/ha/5cm layer'
	WRITE(CHAN_BIOC,240)'Biomass C - values in kgC/ha/5cm layer'
	WRITE(CHAN_DPMC,240)'DPM C - values in kgC/ha/5cm layer'
	WRITE(CHAN_RPMC,240)'RPM C - values in kgC/ha/5cm layer'
      WRITE(CHAN_TOTN,240)'Total soil N - values in kgN/ha/5cm layer'
	WRITE(CHAN_BION,240)'Biomass N - values in kgN/ha/5cm layer'
      WRITE(CHAN_HUMN,240)'Humus N - values in kgN/ha/5cm layer'
      WRITE(CHAN_RPMN,240)'RPM N - values in kgN/ha/5cm layer'
      WRITE(CHAN_DPMN,240)'DPM N - values in kgN/ha/5cm layer'
      WRITE(CHAN_NH4N,240)'Ammonium N - values in kgN/ha/5cm layer'
      WRITE(CHAN_NO3N,240)'Nitrate N - values in kgN/ha/5cm layer'
      WRITE(CHAN_DENITN,240)'Denitrified N - values in kgN/ha/5cm layer'
      WRITE(CHAN_MINERN,240)'Mineralised N - values in kgN/ha/5cm layer'
      WRITE(CHAN_NITRIFN,240)'Nitrified N - values in kgN/ha/5cm layer'
      WRITE(CHAN_CROPN,240)'Crop N - values in kgN/ha/5cm layer'
c      WRITE(CHAN_SENESN,240)'Crop N - values in kgN/ha/5cm layer'
      WRITE(CHAN_LEACHN,240)'Leached N - values in kgN/ha/5cm layer'
      WRITE(CHAN_VOLN,240)'Volatilised N - values in kgN/ha/5cm layer'
      WRITE(CHAN_DOC,240)'DOC - values in kgC/ha/5cm layer'
      WRITE(CHAN_DON,240)'DON - values in kgN/ha/5cm layer'
      WRITE(CHAN_CO2,240)'CO2 - values in kgC/ha/5cm layer'
      WRITE(CHAN_CH4,240)'CH4 - values in kgC/ha/5cm layer'
      WRITE(CHAN_CH4_RICHARDS,250)
     &    'CH4 production in kgC/ha/5cm layer, other values in kgC/ha'
240   FORMAT(A80/
     &       "   Year Step Season  0-5cm     5-10cm    10-15cm    ",
     &       "15-20cm    20-25cm    25-30cm    30-35cm    35-40cm    ",
     &       "40-45cm    45-50cm    50-55cm    55-60cm    60-65cm    ",
     &       "65-70cm    70-75cm    80-85cm    80-85cm    85-90cm    ",
     &       "90-95cm   95-100cm  100-105cm  105-110cm  110-115cm  ",
     &     "115-120cm  120-125cm  125-130cm  130-135cm  135-140cm  ",
     &     "140-145cm  145-150cm  150-155cm  155-160cm  160-165cm  ",
     &     "165-170cm  170-175cm  175-180cm  180-185cm  185-190cm  ",
     &     "190-195cm  195-200cm  200-205cm  205-210cm  210-215cm  ",
     &     "215-220cm  220-225cm  225-230cm  230-235cm  235-240cm  ",
     &     "240-245cm  245-250cm  250-255cm  255-260cm  260-265cm  ",
     &     "265-270cm  270-275cm  275-280cm  280-285cm  285-290cm  ",
     &     "290-295cm  295-300cm    0-300cm")
250   FORMAT(A80/
     &       "   Year Step Season  0-5cm     5-10cm    10-15cm    ",
     &       "15-20cm    20-25cm    25-30cm    30-35cm    35-40cm    ",
     &       "40-45cm    45-50cm    50-55cm    55-60cm    60-65cm    ",
     &       "65-70cm    70-75cm    80-85cm    80-85cm    85-90cm    ",
     &       "90-95cm   95-100cm  100-105cm  105-110cm  110-115cm  ",
     &     "115-120cm  120-125cm  125-130cm  130-135cm  135-140cm  ",
     &     "140-145cm  145-150cm  150-155cm  155-160cm  160-165cm  ",
     &     "165-170cm  170-175cm  175-180cm  180-185cm  185-190cm  ",
     &     "190-195cm  195-200cm  200-205cm  205-210cm  210-215cm  ",
     &     "215-220cm  220-225cm  225-230cm  230-235cm  235-240cm  ",
     &     "240-245cm  245-250cm  250-255cm  255-260cm  260-265cm  ",
     &     "265-270cm  270-275cm  275-280cm  280-285cm  285-290cm  ",
     &     "290-295cm  295-300cm    0-300cm	Soil_Oxidation	",
     &     "Atmos_Oxidation	Flux")
      ENDIF
      END

C
C------------------------------------------------------------------
C
      SUBROUTINE TEST1_RES(IYEAR,IRYEAR,NXYEARS,LHARV,SECONDS,
     &                   IK,N_TS,SUM_TS,FIXEND,NSOW,ISOWN,MEND,
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
     &                   DENITRIFN,NITRIFN,VOLATN,LEACHDOC,LEACHDON,
     &                   THISSOMDEPTH,PLANTUP,FULL_OUTPUT,
     &                   SUMMARY_OUTPUT)
C
C Subroutine to record weeks results
C

C
C Variables local to this subroutine
C
      INTEGER MAXCROP
	INTEGER MAXDEPTH
	DATA MAXDEPTH /300/ 
      PARAMETER (MAXCROP=360)
	INTEGER MAXLAYER
	PARAMETER (MAXLAYER=60)
      INTEGER MAXLAYER1
	DATA MAXLAYER1 /60/
      INTEGER MAXSOIL
	PARAMETER (MAXSOIL=50)
	INTEGER CH4_OFF				! CH4 model off
	INTEGER CH4_RICHARDS	    ! Richards CH4 model on		
	INTEGER CH4_AITKENHEAD      ! Aitkenhead CH4 model on
	DATA CH4_OFF,CH4_RICHARDS,CH4_AITKENHEAD /0,1,2/
 
	REAL ANNCO2
      REAL FTEMP
	REAL FTEMP15
      INTEGER IL
	INTEGER MEASLAY				! Measurement layer number
	REAL TOTBIO
      REAL TOTBION
	REAL TOTCH4
	REAL TOTCH4_RICHARDS
	REAL TOTCO2
	REAL TOTCROPN
      REAL TOTDENITN
	REAL TOTIOM
	REAL TOTDOC
	REAL TOTDON
	REAL TOTDPM
	REAL TOTDPMN
	REAL TOTHUM
	REAL TOTHUMN
	REAL TOTLEACHN
	REAL TOTMINERN
	REAL TOTNH4N
	REAL TOTNITRIFN
	REAL TOTNO3N
	REAL TOTRPM
	REAL TOTRPMN
	REAL TOTSENESN
	REAL TOTVOLN
	REAL TOTSOILW               ! Total available water in soil column 
C
C Variables passed to/from subroutine
C
	REAL AMMN(MAXLAYER)
	REAL BCARB0(MAXLAYER)
	REAL BNIT0(MAXLAYER)
	REAL C_TS					! IN:Litter C input in this timestep (kgC/ha)
	INTEGER CH4MODEL			! IN:Methane model (off, Richards or Aitkenhead CH4 model) 
	REAL CH4TOAIR				! IN:CH4 release to atmosphere (kgC/ha) (Aitkenhead model CH4 model)
	REAL CH4_PROD(MAXLAYER)  	! IN:CH4 produced in each soil layer [kgC/ha/layer] (Richards CH4 model)
	REAL CH4_SOIL_OX			! IN:CH4 oxidised from methane produced in the soil [kgC/ha] (Richards CH4 model)
      REAL CH4_ATMOS_OX			! IN:CH4 oxidised from atmospheric CH4 [kgC/ha] (Richards CH4 model)
      REAL CH4_FLUX				! IN:Net CH4 flux from soil column [kgC/ha] (Richards CH4 model)
	REAL CH4(MAXLAYER)			! IN:CH4 emitted (kgC/ha/layer) (Aitkenhead model only)
      INTEGER CHAN_TOTC			! Channel used to output TOTC and IOM
      INTEGER CHAN_BIOC			! Channel used to output BIOC
      INTEGER CHAN_BION			! Channel used to output BION
	INTEGER CHAN_CH4			! Channel used to output CH4
	INTEGER CHAN_CH4_RICHARDS	! Channel used to output CH4 from the Richards model
	INTEGER CHAN_CO2			! Channel used to output CO2
      INTEGER CHAN_CROPN			! Channel used to output CROPN
      INTEGER CHAN_DENITN			! Channel used to output DENITN
	INTEGER CHAN_DOC			! Channel used to output DOC
	INTEGER CHAN_DON			! Channel used to output DON
      INTEGER CHAN_DPMC			! Channel used to output DPMC
      INTEGER CHAN_DPMN			! Channel used to output DPMN
      INTEGER CHAN_HUMC			! Channel used to output HUMC
      INTEGER CHAN_HUMN			! Channel used to output HUMN
      INTEGER CHAN_LEACHN			! Channel used to output LEACHN
      INTEGER CHAN_MINERN			! Channel used to output MINERN
      INTEGER CHAN_NH4N			! Channel used to output NH4N
      INTEGER CHAN_NITRIFN		! Channel used to output NITRIFN
      INTEGER CHAN_NO3N			! Channel used to output NO3N
      INTEGER CHAN_RPMC			! Channel used to output RPMC
      INTEGER CHAN_RPMN			! Channel used to output RPMN
      INTEGER CHAN_SENESN			! Channel used to output SENESN
      INTEGER CHAN_SOILW			! Channel used to output SOILW
      INTEGER CHAN_SUMMARY      ! Channel used for summary output
      INTEGER CHAN_TOTN			! Channel used to output TOTN 
      INTEGER CHAN_VOLN			! Channel used to output VOLN
	REAL CO2(MAXLAYER)
	REAL CO2FROMDOC(MAXLAYER)
	REAL DENIT
	REAL DENITRIFN(MAXLAYER)	! IN: N lost by denitrification (kgN/ha/layer
	REAL DNIT(MAXLAYER)			! IN:Net Mineralised N (kgN/ha/lay/timestep)
	REAL DPMCARB0(MAXLAYER)
      REAL DPMNIT0(MAXLAYER)
	INTEGER FIXEND
	LOGICAL, INTENT(IN) :: FULL_OUTPUT     ! Flags whether detailed output files should be written
	REAL G15DN2
	REAL G15DN2O									
	REAL G15NN2O
	REAL G15NNO
	REAL G15PNN2O					
	REAL GDN2
	REAL GDN2O
	REAL GNN2O
	REAL GNNO
	REAL GPNN2O
	REAL HCARB0(MAXLAYER)
	REAL HNIT0(MAXLAYER)
	INTEGER I
	INTEGER ICROP
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
	INTEGER ISOWN
	REAL LEACHDOC(MAXLAYER)		! IN/OUT:DOC leached from the layer (kgC/ha/layer/timestep)
	REAL LEACHDON(MAXLAYER)		! IN/OUT:DON leached from the layer (kgN/ha/layer/timestep)
	INTEGER LHARV				! IN:No.timesteps from 01/01 to prev. harvest
	REAL MOBDOC(MAXLAYER)
	REAL MOBDON(MAXLAYER)
	INTEGER N_TS
	REAL NITRIFN(MAXLAYER)		! IN/OUT:Nitrification from this layer (kg/ha/layer/timestep)
	REAL PLANTUP(MAXLAYER)		! IN(CALL):Plant uptake from each soil layer
	REAL RNIN
	REAL RPMCARB0(MAXLAYER)
	REAL RPMNIT0(MAXLAYER)
      REAL SECONDS
	REAL SEEDIN
	REAL SLEACH(MAXLAYER)
	REAL SOILN(MAXLAYER)
	REAL SOILW(MAXLAYER)
	LOGICAL, INTENT(IN) :: SUMMARY_OUTPUT  ! Flags whether summary output file should be written
      INTEGER SUM_TS
	REAL TCINP
	REAL TFYMC
	REAL THISSOMDEPTH			! Measurement depth (cm)
	REAL TOTAST
	REAL TOTAST15
	REAL TOTNST
	REAL TOTNST15
        real TOTSOC0TO30
        real TOTSOC30TO100
	REAL VOLATN(MAXLAYER)		! IN/OUT:Volatilisation from this layer (kg/ha/layer/timestep)
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)

      PARAMETER (CHAN_TOTN=124,CHAN_BION=125,CHAN_HUMN=126,
     &           CHAN_RPMN=127,CHAN_DPMN=128,CHAN_NH4N=129,
     &           CHAN_NO3N=130,CHAN_DENITN=131,CHAN_MINERN=132,
     &           CHAN_NITRIFN=133,CHAN_CROPN=134,CHAN_LEACHN=135,
     &           CHAN_SENESN=136,CHAN_VOLN=137,
     &           CHAN_CO2=138,CHAN_CH4=139,CHAN_CH4_RICHARDS=140)
      PARAMETER (CHAN_TOTC=24,CHAN_BIOC=25,CHAN_HUMC=26,CHAN_RPMC=27,
     &           CHAN_DPMC=28,CHAN_SOILW=20,CHAN_DOC=23,CHAN_DON=29,
     &           CHAN_SUMMARY=37)				
C
C Common block
C
      INCLUDE 'CROPPAR.FOR'
      INCLUDE 'FIXER.FOR'
      INCLUDE 'LASTHARV.FOR'
      INCLUDE 'LEVEL.FOR'
      INCLUDE 'N15.FOR'
      INCLUDE 'RES.FOR'
      INCLUDE 'TOTVAR.FOR'
      INCLUDE 'TOTVAR15.FOR'
C
	SAVE
C
C Get measurement layer
C
 	MEASLAY=(INT(THISSOMDEPTH))*MAXLAYER1/MAXDEPTH
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
      TOTBION=0
	TOTHUMN=0
	TOTRPMN=0
	TOTDPMN=0
	TOTNH4N=0
	TOTNO3N=0
      TOTDENITN=0
	TOTNITRIFN=0
	TOTSENESN=0
	TOTLEACHN=0
	TOTVOLN=0
	TOTMINERN=0
	TOTDOC=0
	TOTDON=0
	TOTN20=0
	TOTSOILW=0
        TOTSOC0TO30=0
        TOTSOC30TO70=0
	DO 100 IL=1,MEASLAY
	  TOTDPM=TOTDPM+DPMCARB0(IL)
	  TOTRPM=TOTRPM+RPMCARB0(IL)
	  TOTBIO=TOTBIO+BCARB0(IL)
	  TOTHUM=TOTHUM+HCARB0(IL)
	  TOTCO2=TOTCO2+CO2(IL)
	  IF(CH4MODEL.EQ.CH4_RICHARDS)THEN
	    TOTCH4=TOTCH4+CH4_PROD(IL)
	  ELSEIF(CH4MODEL.EQ.CH4_AITKENHEAD) THEN
	    TOTCH4=TOTCH4+CH4(IL)
	  ENDIF
	  TOTIOM=TOTIOM+IOM(IL)
          if (IL==6) then
              TOTSOC0TO30=TOTDPM+TOTRPM+TOTBIO+TOTHUM+TOTIOM
          endif
          if (IL==20) then
              TOTSOC30TO100=(TOTDPM+TOTRPM+TOTBIO+TOTHUM+TOTIOM) - TOTSOC0TO30
          endif    
        TOTBION=TOTBION+BNIT0(IL)
	  TOTHUMN=TOTHUMN+HNIT0(IL)
	  TOTRPMN=TOTRPMN+RPMNIT0(IL)
	  TOTDPMN=TOTDPMN+DPMNIT0(IL)
	  TOTNH4N=TOTNH4N+AMMN(IL)
	  TOTNO3N=TOTNO3N+SOILN(IL)
        TOTDENITN=TOTDENITN+DENITRIFN(IL)
	  TOTNITRIFN=TOTNITRIFN+NITRIFN(IL)
	  TOTSENESN=0
	  TOTVOLN=TOTVOLN+VOLATN(IL)
	  TOTLEACHN=TOTLEACHN+SLEACH(IL)
	  TOTMINERN=TOTMINERN+DNIT(IL)
	  TOTDOC=TOTDOC+MOBDOC(IL)
	  TOTDON=TOTDON+MOBDON(IL)
	  TOTSOILW=TOTSOILW+SOILW(IL)
100   CONTINUE
      IF(IK.EQ.1)ANNCO2=0
	ANNCO2=ANNCO2+TOTCO2
C
C Write out summary results file
C
      IF (SUMMARY_OUTPUT) THEN
		WRITE(CHAN_SUMMARY,245)IRYEAR+1,1+(IK+LHARV-1)-IRYEAR*12,
     &       IYEAR,
     &	   TOTDPM,TOTRPM,TOTBIO,TOTHUM,TOTIOM,
     &       TOTDPM+TOTRPM+TOTBIO+TOTHUM+TOTIOM,
     &       TOTSOC0TO30, TOTSOC30TO100,
     &       TOTCO2,CH4_FLUX,
     &	   TOTDPMN,TOTRPMN,TOTBION,TOTHUMN,
     &       TOTDPMN+TOTRPMN+TOTBION+TOTHUMN,
     &       TOTNH4N,TOTNO3N,TOTDENITN,TOTNITRIFN,
     &       TOTVOLN,TOTLEACHN,TOTMINERN,
     &       GDN2O+GNN2O+GPNN2O,GNNO,
     &       TOTDOC,TOTDON,TOTSOILW
      ENDIF
245   FORMAT(3I5,8(F12.1),2(F12.5),F10.2,5(F10.1),
     &       5(F10.3),5(F12.5),F10.1)
      IF (FULL_OUTPUT) THEN
C Write out the C results
C
      FTEMP=-999
	IRYEAR=AINT((IK+LHARV-1)/12.)
C      WRITE(57,5701)IRYEAR+1,(IK+LHARV)-IRYEAR*12,IYEAR,
C     &             (DPMCARB0(IL),IL=1,MEASLAY),
C     &             (FTEMP,IL=1,MAXLAYER1-MEASLAY),
C     &             (RPMCARB0(IL),IL=1,MEASLAY),
C     &             (FTEMP,IL=1,MAXLAYER1-MEASLAY),
C     &             (BCARB0(IL),IL=1,MEASLAY),
C     &             (FTEMP,IL=1,MAXLAYER1-MEASLAY),
C     &             (HCARB0(IL),IL=1,MEASLAY),
C     &             (FTEMP,IL=1,MAXLAYER1-MEASLAY),
C     &             (CO2(IL),IL=1,MEASLAY),
C     &             (FTEMP,IL=1,MAXLAYER1-MEASLAY),
C     &             (CH4(IL),IL=1,MEASLAY),
C     &             (FTEMP,IL=1,MAXLAYER1-MEASLAY),
C     &              TOTDPM,TOTRPM,TOTBIO,TOTHUM,TOTCO2,CH4TOAIR,TOTIOM,
C     &              TCINP+TFYMC,ANNCO2
 5701 FORMAT(1X,3I5,369(F10.1,1X))
	TFYMC = 0.0

C
	IRYEAR=AINT((IK+LHARV-1)/12.)
	WRITE(21,10)IRYEAR+1,(IK+LHARV)-IRYEAR*12,IYEAR,
     &                     (SOILN(I),I=1,MEASLAY),
     &                     (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                     (AMMN(I),I=1,MEASLAY),
     &                     (FTEMP,IL=1,MAXLAYER1-MEASLAY)
10    FORMAT(I5,I5,120F10.2)
	IF(IYEAR.EQ.1.AND.IK.EQ.1)THEN
	  WRITE(CHAN_SOILW,*)'Available water at saturation (mm)'
	  WRITE(CHAN_SOILW,240)IRYEAR+1,(IK+LHARV)-IRYEAR*12,IYEAR,
     &                       (WSAT(I),I=1,MAXLAYER1)
	  WRITE(CHAN_SOILW,*)'Available water at field capacity(mm)'
	  WRITE(CHAN_SOILW,240)IRYEAR+1,(IK+LHARV)-IRYEAR*12,IYEAR,
     &                       (WMAX(I),I=1,MAXLAYER1)
	  WRITE(CHAN_SOILW,*)'Available water (mm)'
	ENDIF
	WRITE(CHAN_SOILW,240)IRYEAR+1,(IK+LHARV)-IRYEAR*12,IYEAR,
     &                     (SOILW(I),I=1,MAXLAYER1)
C
C Output of TOTC, HUMC, BIOC, DPMC, RPMC
C
      IF(IK.EQ.1.AND.IYEAR.EQ.1)THEN 
        WRITE(CHAN_TOTC,*)'Inert Organic matter in kgC/ha/5cm layer'
        WRITE(CHAN_TOTC,240)IRYEAR+1,(IK+LHARV)-IRYEAR*12,IYEAR,
     &                     (IOM(IL),IL=1,MEASLAY),
     &                     (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                      TOTIOM
        WRITE(CHAN_TOTC,*)'Total soil C in kgC/ha/5cm layer'
	ENDIF
	WRITE(CHAN_TOTC,240)IRYEAR+1,(IK+LHARV)-IRYEAR*12,IYEAR,
     &       ((IOM(IL)+HCARB0(IL)+BCARB0(IL)+DPMCARB0(IL)+RPMCARB0(IL)),
     &        IL=1,MEASLAY),
     &        (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &        TOTIOM+TOTHUM+TOTBIO+TOTDPM+TOTRPM
      WRITE(CHAN_BIOC,240)IRYEAR+1,1+(IK+LHARV-1)-IRYEAR*12,IYEAR,
     &                    (BCARB0(IL),IL=1,MEASLAY),
     &                    (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                     TOTBIO
      WRITE(CHAN_HUMC,240)IRYEAR+1,1+(IK+LHARV-1)-IRYEAR*12,IYEAR,
     &                    (HCARB0(IL),IL=1,MEASLAY),
     &                    (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                     TOTHUM
      WRITE(CHAN_DPMC,240)IRYEAR+1,1+(IK+LHARV-1)-IRYEAR*12,IYEAR,
     &                    (DPMCARB0(IL),IL=1,MEASLAY),
     &                    (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                     TOTDPM
      WRITE(CHAN_RPMC,240)IRYEAR+1,1+(IK+LHARV-1)-IRYEAR*12,IYEAR,
     &                    (RPMCARB0(IL),IL=1,MEASLAY),
     &                    (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                     TOTRPM
      WRITE(CHAN_TOTN,240)IRYEAR+1,1+(IK+LHARV-1)-IRYEAR*12,IYEAR,
     &       ((HNIT0(IL)+BNIT0(IL)+DPMNIT0(IL)+RPMNIT0(IL)),
     &        IL=1,MEASLAY),
     &        (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &        TOTHUMN+TOTBION+TOTDPMN+TOTRPMN
      WRITE(CHAN_BION,240)IRYEAR+1,1+(IK+LHARV-1)-IRYEAR*12,IYEAR,
     &                    (BNIT0(IL),IL=1,MEASLAY),
     &                    (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                    TOTBION
      WRITE(CHAN_HUMN,240)IRYEAR+1,1+(IK+LHARV-1)-IRYEAR*12,IYEAR,
     &                    (HNIT0(IL),IL=1,MEASLAY),
     &                    (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                     TOTHUMN
      WRITE(CHAN_RPMN,240)IRYEAR+1,1+(IK+LHARV-1)-IRYEAR*12,IYEAR,
     &                    (RPMNIT0(IL),IL=1,MEASLAY),
     &                    (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                     TOTRPMN
      WRITE(CHAN_DPMN,240)IRYEAR+1,1+(IK+LHARV-1)-IRYEAR*12,IYEAR,
     &                    (DPMNIT0(IL),IL=1,MEASLAY),
     &                    (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                     TOTDPMN
      WRITE(CHAN_NH4N,240)IRYEAR+1,1+(IK+LHARV-1)-IRYEAR*12,IYEAR,
     &                    (AMMN(IL),IL=1,MEASLAY),
     &                    (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                     TOTNH4N
      WRITE(CHAN_NO3N,240)IRYEAR+1,1+(IK+LHARV-1)-IRYEAR*12,IYEAR,
     &                    (SOILN(IL),IL=1,MEASLAY),
     &                    (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                     TOTNO3N
      WRITE(CHAN_DENITN,240)IRYEAR+1,1+(IK+LHARV-1)-IRYEAR*12,IYEAR,
     &                    (DENITRIFN(IL),IL=1,MEASLAY),
     &                    (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                     TOTDENITN
      WRITE(CHAN_MINERN,240)IRYEAR+1,1+(IK+LHARV-1)-IRYEAR*12,IYEAR,
     &                    (DNIT(IL),IL=1,MEASLAY),
     &                    (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                     TOTMINERN
      WRITE(CHAN_NITRIFN,240)IRYEAR+1,1+(IK+LHARV-1)-IRYEAR*12,IYEAR,
     &                    (NITRIFN(IL),IL=1,MEASLAY),
     &                    (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                     TOTNITRIFN
      WRITE(CHAN_CROPN,240)IRYEAR+1,1+(IK+LHARV-1)-IRYEAR*12,IYEAR,
     &                    (PLANTUP(IL),IL=1,MEASLAY),
     &                    (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                     TOTBION
      IF(MEASLAY.LT.MAXLAYER1)THEN 
        WRITE(CHAN_LEACHN,240)IRYEAR+1,1+(IK+LHARV-1)-IRYEAR*12,IYEAR,
     &                    (SLEACH(IL),IL=1,MEASLAY),
     &                    (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                     SLEACH(MEASLAY+1)
	ELSE
        WRITE(CHAN_LEACHN,240)IRYEAR+1,1+(IK+LHARV-1)-IRYEAR*12,IYEAR,
     &                    (SLEACH(IL),IL=1,MEASLAY),
     &                    (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                     WLEACH
      ENDIF
C      WRITE(CHAN_SENES,240)IRYEAR+1,1+(IK+LHARV-1)-IRYEAR*12,IYEAR,
C     &                    (BNIT0(IL),IL=1,MAXLAYER1),TOTBION
      WRITE(CHAN_VOLN,240)IRYEAR+1,1+(IK+LHARV-1)-IRYEAR*12,IYEAR,
     &                    (VOLATN(IL),IL=1,MEASLAY),
     &                    (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                     TOTVOLN
      WRITE(CHAN_CO2,240)IRYEAR+1,1+(IK+LHARV-1)-IRYEAR*12,IYEAR,
     &                    (CO2(IL),IL=1,MEASLAY),
     &                    (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                     TOTCO2
      IF(CH4MODEL.EQ.CH4_RICHARDS)THEN
		WRITE(CHAN_CH4_RICHARDS,250)IRYEAR+1,1+(IK+LHARV-1)-
     &                    IRYEAR*12,IYEAR,
     &                    (CH4_PROD(IL),IL=1,MEASLAY),
     &                    (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                     TOTCH4,CH4_SOIL_OX,CH4_ATMOS_OX,CH4_FLUX
	ELSEIF(CH4MODEL.EQ.CH4_AITKENHEAD)THEN
		WRITE(CHAN_CH4,240)IRYEAR+1,1+(IK+LHARV-1)-IRYEAR*12,
     &					IYEAR,(CH4(IL),IL=1,MEASLAY),
     &                    (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                     TOTCH4 
	ENDIF
	WRITE(CHAN_DOC,240)IRYEAR+1,(IK+LHARV)-IRYEAR*12,IYEAR,
     &                     (MOBDOC(I),I=1,MEASLAY),
     &                     (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                      TOTDOC
	WRITE(CHAN_DON,240)IRYEAR+1,(IK+LHARV)-IRYEAR*12,IYEAR,
     &                     (MOBDON(I),I=1,MEASLAY),
     &                     (FTEMP,IL=1,MAXLAYER1-MEASLAY),
     &                      TOTDON
240   FORMAT(1X,3I5,62(F10.1,1X))
250   FORMAT(1X,3I5,64(F10.5,1X))

C
C Initialize variables at sowing
C
      IF(IK.LT.NSOW)THEN
        CACL=CACTOT
      ELSE
        CACL=0
      ENDIF
C
C Record name of setup file
C
      IF(IRYEAR.EQ.0.AND.N_TS.EQ.IANTHES)THEN
       WRITE(22,5)
5      FORMAT(//'TotN Min '/
     &        'Year, Week, Day, 0-23cm, 23-50cm, 50-70cm, 70-100cm ') 
      ENDIF
C
      START=RSTART+BSTART+TOTNST+TOTAST
C
C Sum inputs during week
C
      WATM=ATM+FIXN
      SUBTOTIN=THISFERT+WATM+TOTMINERN+FYMFERT
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
      DO 400 IL=1,MEASLAY
        TOTNIT=TOTNIT+SOILN(IL)
400   CONTINUE
C
C Ammonium...
C
      DO 500 IL=1,MEASLAY
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
      DO 550 IL=1,MEASLAY
        ROOTN=ROOTN+DPMNIT0(IL)+RPMNIT0(IL)
	  ROOTDPM=ROOTDPM+DPMNIT0(IL)
	  ROOTRPM=ROOTRPM+RPMNIT0(IL)
        BIOHUM=BIOHUM+BNIT0(IL)+HNIT0(IL)
550   CONTINUE
      END=ROOTN+BIOHUM+TOTNIT+TOTAMM
C
C Sum outputs during week
C
C      SUBTOTOUT=DENIT+CACT+WLEACH+VOLAT-CLOSSX+LEACHDON(MAXLAYER1)
      SUBTOTOUT=DENIT+CACT+WLEACH+VOLAT-CLOSSX
C
C Sum N at end of week
C
      TOTOUT=END+SUBTOTOUT
C
C Write out Weekly N Balance sheet
C
      BALIN=BALIN+THISFERT+FYMFERT+TOTMINERN+WATM+SEEDIN
	BALOUT=BIOHUM+ROOTDPM+ROOTRPM+TOTNIT+TOTAMM
	BALOUT=BALOUT+TOTDENITN+TOTVOLN-CLOSSX+CACT+TOTLEACHN
	BALOUT=BALOUT+LEACHDON(MEASLAY)
	DIFF=BALIN-BALOUT
      WRITE(58,5801) N_TS,IK,IYEAR, BIOHUM,ROOTN,
     &    ROOTDPM,ROOTRPM,TOTNIT,TOTAMM,THISFERT+FYMFERT,RNIN,
     &    WATM,SEEDIN,TOTMINERN,TOTDENITN,
     &    GNN2O,GNNO,GPNN2O,G15NN2O,G15NNO,
     &    G15PNN2O,GDN2,GDN2O,G15DN2,G15DN2O,TOTVOLN-CLOSSX,CACT,
     &    TOTLEACHN,
     &    CLOSSX,LEACHDON(MEASLAY),BALIN,BALOUT,DIFF
	BALIN=BIOHUM+ROOTDPM+ROOTRPM+TOTNIT+TOTAMM

 5801 FORMAT (2X,3I5,1X,F10.1,5(F10.2,1X),16(F15.9,1X),5(F10.2,1X),
     &        3(F10.1,1X))

      ENDIF
      RETURN
      END
C
C--------------------------------------------------------
C
	SUBROUTINE TEST2_OPENCHAN()
C
C Subroutine to open channels
C
      IMPLICIT NONE
C
C Open Channels
C
      OPEN(41,FILE='CCHANGE.OUT',STATUS='UNKNOWN')
C
C Titles in results files
C
      WRITE(41,10)
10    FORMAT('20KM2_GRID       SOIL  LAND_USE1  LAND_USE2 '
     &       'LU1_to_LU2 ',
     &       'CCHANGE(t/ha/10y)1 CO2(tC/ha)1 CH4(tC/ha)1 ',
     &       'N2O(tN/ha)1 CEquiv(tC/ha)1 ',
     &       'CCHANGE(t/ha/10y)2 CO2(tC/ha)2 CH4(tC/ha)2 ',
     &       'N2O(tN/ha)2 CEquiv(tC/ha)2 ',
     &       'CCHANGE(t/ha/10y)3 CO2(tC/ha)3 CH4(tC/ha)3 ',
     &       'N2O(tN/ha)3 CEquiv(tC/ha)3 ',
     &       'CCHANGE(t/ha/10y)4 CO2(tC/ha)4 CH4(tC/ha)4 ',
     &       'N2O(tN/ha)4 CEquiv(tC/ha)4 ',
     &       'CCHANGE(t/ha/10y)5 CO2(tC/ha)5 CH4(tC/ha)5 ',
     &       'N2O(tN/ha)5 CEquiv(tC/ha)5 ',
     &       'CCHANGE(t/ha/10y)6 CO2(tC/ha)6 CH4(tC/ha)6 ',
     &       'N2O(tN/ha)6 CEquiv(tC/ha)6 ',
     &       'CCHANGE(t/ha/10y)7 CO2(tC/ha)7 CH4(tC/ha)7 ',
     &       'N2O(tN/ha)7 CEquiv(tC/ha)7 ',
     &       '  %C    %Clay   Error?'/)
      END

C
C------------------------------------------------------------------
C
      SUBROUTINE TEST2_RES(KM20GRIDID,SOILID,NSEQ,IWC,
     &                     SOMDEPTH,ISERIES,CCHANGE10,
     &                     CO2C10,CH4C10,N2ON10,
     &                     LU1TOLU2,
     &                     DOMSOILC,DOMSOILCLAY,DOMSOILBD,
     &                     NUMCELLSIN20KM2,
     &	   			     LU1,LU2,TRAPERR)
C Subroutine to record weeks results
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
      INTEGER MAXDEC				! Max.no.of decades included in calculation
	PARAMETER(MAXDEC=9)		
      INTEGER MAXDEC1				! Max.no.of decades included in calculation
	DATA MAXDEC1 /9/		
	INTEGER MAXGROW				! Max.no.of growing seasons allowed
	PARAMETER (MAXGROW=300)
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXLU1				! Max.no.of land use types
	DATA MAXLU1 /6/
	INTEGER MAXSEQ				! Max.no of sequences included in GIS
	PARAMETER (MAXSEQ=MAXLU*MAXLU)						
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXSOMLAY			! Max.no.of soil layers entered			
	PARAMETER (MAXSOMLAY=10)
	INTEGER MAXWC				! Max.no of wetness classes
      PARAMETER (MAXWC=6)

      REAL AWC					! Atomic weight of C g/mol
	REAL AWN					! Atomic weight of N g/mol
	DATA AWC,AWN /12,14/
      REAL GWPCO2					! Global warming potential of CO2 /mol
      REAL GWPCH4					! Global warming potential of CO2 /mol
      REAL GWPN2O					! Global warming potential of CO2 /mol
	DATA GWPCO2,GWPCH4,GWPN2O /1,23,296/
	INTEGER IDEC				! Local decade counter
	REAL MWCH4					! Molecular weight of CH4 g/mol
	REAL MWCO2					! Molecular weight of CO2 g/mol
	REAL MWN2O					! Molecular weight of N2O g/mol
	DATA MWCO2,MWCH4,MWN2O /44,16,44/

      REAL CEQUIV(MAXDEC)			! Total emissions in carbon equiv.(t/ha/10y)
      INTEGER ICELL				! Local cell counter
	INTEGER IS					! Local counter for soil series
	REAL PERCENTC				! % C in the soil layer
	REAL TEMP					! Local real used in calculation
C
C Variables passed to / from this subroutine
C
	REAL CCHANGE10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! IN:Change in C 
								! content over the decade (kgC/ha/decade)
	REAL CH4C10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! IN:CH4-C emitted 
								! over the decade (kgC/ha/decade)
	REAL CO2C10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! IN:CO2-C emitted 
								! over the decade (kgC/ha/decade)
	REAL DOMSOILC(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! Soil C in 5 
												! major soil series under different LU 
												! in SOM layers (kgC/ha)
	REAL DOMSOILCLAY(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! % Soil clay in 5 
												! major soil series under different LU 
	REAL DOMSOILBD(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! Soil BD in 5 
												! major soil series under different LU 
												! in SOM layers (g/cm3)
									! LU in SOM layers	// (1c) Del
	INTEGER ISERIES				! IN:Soil series code number
	INTEGER IWC					! IN: Counter for wetness classes
	CHARACTER*10 KM20GRIDID		! IN:20km sqare identifier
	INTEGER	LU1					! IN:Starting land use
	INTEGER	LU2					! IN:Land use after change
	REAL LU1TOLU2(CELLSTOGRID,MAXDEC,MAXLU+2,MAXLU+2)	! Fracn of LU1 changed 
	REAL N2ON10(MAXDEC+1,CELLSTOGRID*MAXSERIES/2,MAXWC,MAXSEQ)	! IN:N2O-N emitted 
								! over the decade (kgN/ha/decade)
	INTEGER NSEQ				! IN:Counter for number of land use sequences 
	INTEGER NUMCELLSIN20KM2		! IN:Number of cells in 20km2 cell
	INTEGER SEQ1				! IN:Land use sequence of LU1->LU1
	INTEGER SEQTYPE				! Type of sequence to be used
      INTEGER SOILID(CELLSTOGRID*MAXSERIES/2)	! IN:Soil series integer codes
	REAL SOMDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN:depth of soil organic matter layers
      INTEGER TRAPERR(CELLSTOGRID*MAXSERIES/2,MAXSEQ) 
								! IN:Error trapping
								!		0=No error
								!		1=Error in NPP data
								!		2=Error in weather data
								!		3=Error in LU data
								!		4=Error in LU change data
								!		5=Soil error
								!		6=Equilibrium not found
C
C Output C change
C Get C equivalents
C
      DO 100 IDEC=1,MAXDEC1
	  SEQ1=((LU1-1)*MAXLU1)+LU1
        TEMP=CO2C10(IDEC,ISERIES,IWC,NSEQ)-CO2C10(IDEC,ISERIES,IWC,SEQ1)
	  CEQUIV(IDEC)=TEMP*GWPCO2*MWCO2/AWC
	  TEMP=CH4C10(IDEC,ISERIES,IWC,NSEQ)-CH4C10(IDEC,ISERIES,IWC,SEQ1)
	  CEQUIV(IDEC)=CEQUIV(IDEC)+(TEMP*GWPCH4*MWCH4/AWC)
	  TEMP=N2ON10(IDEC,ISERIES,IWC,NSEQ)-N2ON10(IDEC,ISERIES,IWC,SEQ1)
	  CEQUIV(IDEC)=CEQUIV(IDEC)+(TEMP*GWPN2O*MWN2O/AWN)
	  CEQUIV(IDEC)=CEQUIV(IDEC)*AWC/(MWCO2*1000)
100   CONTINUE
C
C Get percent C
C
	PERCENTC=DOMSOILC(ISERIES,LU1,1)/DOMSOILBD(ISERIES,LU1,1)
	PERCENTC=PERCENTC/(1000*SOMDEPTH(ISERIES,LU1,1))
C
C Get C Change for LU change scenarios
C
c      IF(LU1.NE.LU2)THEN
c   	  IF(LU1TOLU2(1,5,LU1,LU2).GT.0)
c     &    WRITE(41,10)KM20GRIDID,SOILID(ISERIES),LU1,LU2,
c          WRITE(41,10)KM20GRIDID,SOILID(ISERIES),LU1,LU2,
c     &    LU1TOLU2(1,1,LU1,LU2),
c     &    (((CCHANGE10(IDEC,ISERIES,IWC,NSEQ)
c     &       -CCHANGE10(IDEC,ISERIES,IWC,SEQ1))/1000,
c     &    (CO2C10(IDEC,ISERIES,IWC,NSEQ)-CO2C10(IDEC,ISERIES,IWC,SEQ1))
c     &      /1000,
c     &    (CH4C10(IDEC,ISERIES,IWC,NSEQ)-CH4C10(IDEC,ISERIES,IWC,SEQ1))
c     &      /1000,
c     &    (N2ON10(IDEC,ISERIES,IWC,NSEQ)-N2ON10(IDEC,ISERIES,IWC,SEQ1))
c     &      /1000,
c     &    CEQUIV(IDEC)),IDEC=1,MAXDEC1),PERCENTC,
c     &    DOMSOILCLAY(ISERIES,LU1,1),TRAPERR(ISERIES,NSEQ)
10        FORMAT(A10,1X,3(I10,1X),F10.7,7X,5(7(F10.5,2X)),
     &           2(F7.2,2X),I7/)
C
C Get C Change for future climate change scenarios
C
c      ELSE
c       WRITE(41,10)KM20GRIDID,SOILID(ISERIES),LU1,LU2,
c    &    LU1TOLU2(1,5,LU1,LU2),
c    &    (((CCHANGE10(IDEC,ISERIES,IWC,NSEQ))/1000,
c    &      (CO2C10(IDEC,ISERIES,IWC,NSEQ))/1000,
c    &      (CH4C10(IDEC,ISERIES,IWC,NSEQ))/1000,
c    &      (N2ON10(IDEC,ISERIES,IWC,NSEQ))/1000,
c    &       CEQUIV(IDEC)),IDEC=1,MAXDEC1),PERCENTC,
c    &    DOMSOILCLAY(ISERIES,LU1,1),TRAPERR(ISERIES,NSEQ)
c     ENDIF
      RETURN
      END
C
C--------------------------------------------------------
C
	SUBROUTINE TEST3_OPENCHAN()
C
C Subroutine to open channels
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      INTEGER ILU					! Land use counter
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
      INTEGER SOLCHAN
C
C Set channels
C
      SOLCHAN=150
C
C Open Channels
C
      OPEN(SOLCHAN+1,FILE='LEACH_ARA.OUT',STATUS='UNKNOWN')
      OPEN(SOLCHAN+2,FILE='LEACH_GRA.OUT',STATUS='UNKNOWN')
      OPEN(SOLCHAN+3,FILE='LEACH_FOR.OUT',STATUS='UNKNOWN')
      OPEN(SOLCHAN+4,FILE='LEACH_NAT.OUT',STATUS='UNKNOWN')
      OPEN(SOLCHAN+5,FILE='LEACH_MIS.OUT',STATUS='UNKNOWN')
      OPEN(SOLCHAN+6,FILE='LEACH_SRC.OUT',STATUS='UNKNOWN')
C
C Titles in results files
C
      DO 100,ILU=1,MAXLU
        WRITE(SOLCHAN+ILU,10)
100   CONTINUE
10    FORMAT('GRIDID     MONTH YEAR Nitrate(mgN/l) DOC(mgC/l)')
      END

C
C-------------------------------------------------------------
C
      SUBROUTINE TEST3_RES(KM20GRIDID,IK,IYEAR,LU1,
     &                     SLEACH,LEACHDOC,DRAINW)
C
C Subroutine to get filenames and simulation inputs
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	REAL CONCNO3                ! Concentration of nitrate mg/l
	REAL CONCDOC				! Concentration of DOC (lg/l)
	INTEGER LAYER				! Layer that represents the bottom of the profile
	INTEGER MAXLAYER			! No.of layers in the soil profile
      PARAMETER (MAXLAYER=60)
	INTEGER SOLCHAN
C
C Variables passed from calling subroutine
C
	INTEGER IK					! IN(Call):No.timesteps from prev.harv. to current
	INTEGER IYEAR				! IN(Call):Current growing season number - set to 1
	CHARACTER*10 KM20GRIDID		! IN(Call):20km sqare identifier
	REAL LEACHDOC(MAXLAYER)		! IN(Call):DOC leached from the layer (kgC/ha/layer/timestep)
      INTEGER LU1					! IN(Call):Land use code
      REAL SLEACH(MAXLAYER)		! IN:Nitrate-N leached (kgN/ha)
	REAL DRAINW(MAXLAYER)		! IN(Call):Water drained from this layer (mm/layer)
C
C Set channel
C
      SOLCHAN=150
	LAYER=8
C
C Calculate concentrations
C
      CONCNO3=0
	CONCDOC=0
      IF(DRAINW(LAYER).GT.0)THEN
        CONCNO3=SLEACH(LAYER)*100/DRAINW(LAYER)
        CONCDOC=LEACHDOC(LAYER)*100/DRAINW(LAYER)
	ENDIF
C
C Record the results
C     
!      WRITE(SOLCHAN+LU1,10)KM20GRIDID,IK,IYEAR,CONCNO3,CONCDOC
!10    FORMAT(A10,2X,I3,3X,I3,2X,F10.2,6X,F10.2)
      END
C
C-----------------------------------------------------------------------
C
C INTERNAL SUBROUTINES
C 1a. FINDSERIES
C 1b. FIND_SSKIB_SERIES
C 2. GET_SOIL_ORDER
C 3. READINGRID
C 4. TRAPERROR
C
C-------------------------------------------------------------
C
      SUBROUTINE FINDSERIES(THISSERIES,PERC,CKG106PKM2,
     &                      CLAY,SILT,SAND,BD,SOIL_PH,SDEPTH,
     &                      IS_IMPLAYER,IMPDEPTH)
C
C Find soils series from series codes file 
C
C BD(MAXLU,MAXSOMLAY) = BD (g/cm3)
C CKG106PKM2(MAXLU,MAXSOMLAY) = C (kg*10^6/km2 0-30cm)
C CLAY(MAXLU,MAXSOMLAY) = % clay 
C CODECHAN() = channel for reading in soil codes
C ISL = Counter for SOM layers
C ILU = Counter for Land uses
C PERC(MAXLU,MAXSOMLAY) = % C 
C SAND(MAXLU,MAXSOMLAY) = % sand
C SERIES = Soil series code number 
C SILT(MAXLU,MAXSOMLAY) = % silt
C THISSOIL = Soil series being looked for
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      INTEGER GISSOMLAY		! Max number of SOM layers in GIS
	PARAMETER (GISSOMLAY=2)
	INTEGER MAXSOMLAY,MAXLU
	INTEGER MAXLU1
	PARAMETER (MAXSOMLAY=10,MAXLU=6)
	DATA MAXLU1 /6/
	INTEGER CODECHAN(4)		! Channel for soil codes 
      INTEGER SERIES
	INTEGER ISL,ILU
C
C Variables passed to / from this routine
C 
      INTEGER THISSERIES				! IN: Code number of this soil series
	REAL PERC(MAXLU,MAXSOMLAY)		! OUT: % carbon
	REAL CKG106PKM2(MAXLU,MAXSOMLAY)! OUT: Carbon (kg*10^6/km2 0-30cm)
	REAL BD(MAXLU,MAXSOMLAY)		! OUT: Bulk density (g/cm3)
	REAL CLAY(MAXLU,MAXSOMLAY)		! OUT: % clay 
	REAL IMPDEPTH(MAXLU)			! OUT: Depth of impermeable layer if present (cm)
	INTEGER IS_IMPLAYER(MAXLU)		! OUT: Flag for impermeable layer (0=No; 1=Yes)
	REAL SAND(MAXLU,MAXSOMLAY)		! OUT: % sand
	REAL SILT(MAXLU,MAXSOMLAY)		! OUT: % silt
	REAL SOIL_PH(MAXLU,MAXSOMLAY)	! OUT: pH of soil in this layer
      REAL SDEPTH(MAXLU,MAXSOMLAY)		! OUT: Depth of each SOM layer
C
C Read in data until this series is found
C
      CODECHAN(1)=61
101   CONTINUE
	READ(CODECHAN(1),*,END=111)SERIES,
     &  PERC(1,1),CKG106PKM2(1,1),CLAY(1,1),SILT(1,1),SAND(1,1),BD(1,1),
     &  PERC(1,2),CKG106PKM2(1,2),CLAY(1,2),SILT(1,2),SAND(1,2),BD(1,2),
     &  PERC(2,1),CKG106PKM2(2,1),CLAY(2,1),SILT(2,1),SAND(2,1),BD(2,1),
     &  PERC(2,2),CKG106PKM2(2,2),CLAY(2,2),SILT(2,2),SAND(2,2),BD(2,2),
     &  PERC(4,1),CKG106PKM2(4,1),CLAY(4,1),SILT(4,1),SAND(4,1),BD(4,1),
     &  PERC(4,2),CKG106PKM2(4,2),CLAY(4,2),SILT(4,2),SAND(4,2),BD(4,2),
     &  PERC(3,1),CKG106PKM2(3,1),CLAY(3,1),SILT(3,1),SAND(3,1),BD(3,1),
     &  PERC(3,2),CKG106PKM2(3,2),CLAY(3,2),SILT(3,2),SAND(3,2),BD(3,2),
     &  PERC(5,1),CKG106PKM2(5,1),CLAY(5,1),SILT(5,1),SAND(5,1),BD(5,1),
     &  PERC(5,2),CKG106PKM2(5,2),CLAY(5,2),SILT(5,2),SAND(5,2),BD(5,2),
     &  PERC(6,1),CKG106PKM2(6,1),CLAY(6,1),SILT(6,1),SAND(6,1),BD(6,1),
     &  PERC(6,2),CKG106PKM2(6,2),CLAY(6,2),SILT(6,2),SAND(6,2),BD(6,2)
      IF(SERIES.EQ.THISSERIES)THEN
	  DO 400 ILU=1,MAXLU1
          DO 300 ISL=1,GISSOMLAY
	      SOIL_PH(ILU,ISL)=7
300       CONTINUE
	    SDEPTH(ILU,1)=30
	    SDEPTH(ILU,2)=100
		IS_IMPLAYER(ILU)=0
		IMPDEPTH(ILU)=300
400       CONTINUE
	  REWIND(CODECHAN(1))
	  READ(CODECHAN(1),*)
	  RETURN
	ENDIF
      GOTO 101
C
C If the series has not been found, set data to 999
C
111   CONTINUE
	REWIND(CODECHAN(1))
	READ(CODECHAN(1),*)
      DO 100 ISL=1,GISSOMLAY
	  DO 200 ILU=1,MAXLU1
          PERC(ILU,ISL)=999
		CKG106PKM2(ILU,ISL)=999
		CLAY(ILU,ISL)=999
		SILT(ILU,ISL)=999
		SAND(ILU,ISL)=999
		BD(ILU,ISL)=999
200     CONTINUE
100   CONTINUE
	REWIND(CODECHAN(1))
	READ(CODECHAN(1),*)
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE FIND_SSKIB_SERIES(THISSERIES,PERC,CKG106PKM2,
     &                      CLAY,SILT,SAND,BD,SOIL_PH,SDEPTH,
     &                      IS_IMPLAYER,IMPDEPTH,GLEY)
C
C Find soils series from series codes file 
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      INTEGER GISSOMLAY				! Max number of SOM layers in GIS
	PARAMETER (GISSOMLAY=9)
	INTEGER MAXLU					! Max.no.of land use types
	PARAMETER (MAXLU=6)	
	INTEGER MAXSOMLAY				! Max.no.of soil layers entered			
	PARAMETER (MAXSOMLAY=10)
	INTEGER MAXLU1					! Max.no.of land use types
	DATA MAXLU1 /6/
	INTEGER CODECHAN(MAXLU)			! Channel for soil codes for LUs
      INTEGER SERIES					! Series code number
	INTEGER ISL						! Counter for SOM layers
	INTEGER ILU						! Counter for Land uses
	REAL DUMMY						! Dummy 
C
C Variables passed to / from this routine
C 
      INTEGER THISSERIES				! IN: Code number of this soil series
	REAL PERC(MAXLU,MAXSOMLAY)		! OUT: % carbon
	REAL CKG106PKM2(MAXLU,MAXSOMLAY)! OUT: Carbon (kg*10^6/km2 0-30cm)
	REAL BD(MAXLU,MAXSOMLAY)		! OUT: Bulk density (g/cm3)
	REAL CLAY(MAXLU,MAXSOMLAY)		! OUT: % clay 
	REAL IMPDEPTH(MAXLU)			! OUT: Depth of impermeable layer if present (cm)
	INTEGER IS_IMPLAYER(MAXLU)		! OUT: Flag for impermeable layer (0=No; 1=Yes)
	REAL SAND(MAXLU,MAXSOMLAY)		! OUT: % sand
	REAL SILT(MAXLU,MAXSOMLAY)		! OUT: % silt
	REAL SOIL_PH(MAXLU,MAXSOMLAY)	! OUT: pH of soil in this layer
      REAL SDEPTH(MAXLU,MAXSOMLAY)		! OUT: Depth of each SOM layer
      REAL GLEY(MAXLU)		! OUT: Gleying depth, analogous to saturated water level -> WTABLE
C
C Read in data until this series is found
C
	CODECHAN(1)=61
	CODECHAN(2)=62
	CODECHAN(3)=63
	CODECHAN(4)=64
C
C For each LU...
C  
      DO 100 ILU=1,MAXLU1
C
C Read through file until series is found
C
101     CONTINUE
c  	  READ(CODECHAN(ILU),*,END=111)SERIES,GLEY(ILU),
c     &        IS_IMPLAYER(ILU),IMPDEPTH(ILU),
c     &      ((DUMMY,SDEPTH(ILU,ISL),DUMMY,SOIL_PH(ILU,ISL),
c     &	    PERC(ILU,ISL),CKG106PKM2(ILU,ISL),
c     &        CLAY(ILU,ISL),SILT(ILU,ISL),SAND(ILU,ISL),
c     &        BD(ILU,ISL),DUMMY),ISL=1,MAXSOMLAY)
C
C If series found, find out number of layers and stop reading through file
C
        IF(SERIES.EQ.THISSERIES)GOTO 102
        GOTO 101
111     CONTINUE
C
C If at end of file and data not found, set to 999 to show not found
C
        DO 300 ISL=1,GISSOMLAY
          PERC(ILU,ISL)=-999
		CKG106PKM2(ILU,ISL)=-9990000
		CLAY(ILU,ISL)=-999
		SILT(ILU,ISL)=-999
		SAND(ILU,ISL)=-999
		BD(ILU,ISL)=-999
300     CONTINUE
C
C Series found
C
102     CONTINUE
C
C Find bottom of measured layer
C
        DO 400 ISL=1,MAXSOMLAY
          IF(PERC(ILU,ISL).LE.-999.OR.
     &       CKG106PKM2(ILU,ISL).LE.-9990000.OR.
     &       CLAY(ILU,ISL).LE.-999.OR.
     &       SILT(ILU,ISL).LE.-999.OR.
     &       SAND(ILU,ISL).LE.-999.OR.
     &       BD(ILU,ISL).LE.-999)THEN
	       GOTO 404
	    ENDIF
400     CONTINUE
404     CONTINUE
C
C Rewind file ready for next series
C 
	  REWIND(CODECHAN(ILU))
	  READ(CODECHAN(ILU),*)
C
C Convert from kg/ha to (kg x 10^6 / km^2)
C
	  DO 200 ISL=1,MAXSOMLAY
	    CKG106PKM2(ILU,ISL)=CKG106PKM2(ILU,ISL)/10000
200     CONTINUE
100   CONTINUE
	END
C
C-----------------------------------------------------------
C
	SUBROUTINE GET_SOIL_ORDER(THISSOIL,SOILID,STYPES,ORDER)
C
C Subroutine to get the soil order
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						

	INTEGER ISG					! Counter for soil series in 20km2 grid
C
C Variables passed from calling subroutine
C
      INTEGER ORDER					!OUT:Order of soils calculations
      INTEGER SOILID(CELLSTOGRID*MAXSERIES/2)	! IN:Soil series integer codes
 	INTEGER STYPES					! IN:Number of soil types in 20km2 grid
      INTEGER THISSOIL				! IN: Soil series code
C
C Search for soil
C
      ORDER=0
      DO 100 ISG=1,STYPES
	  IF(SOILID(ISG).EQ.THISSOIL)THEN
	    ORDER=ISG
	    GOTO 101
        ENDIF
100   CONTINUE
101   CONTINUE
C
C Leave GET_SOIL_ORDER
C
	END
C
C-------------------------------------------------------------
C
	SUBROUTINE READINFUTMET(NUMCELLSIN20KM2,GRIDLAT,
     &                        AVETEMP,AVERAIN,
     &                        FUTTEMP,FUTRAIN,FUTPET)
C
C Subroutine to read and interpolate in future weather data
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
      INTEGER MAXMET				! No. met.years
      PARAMETER (MAXMET=90)
      INTEGER FUTCHAN				! Channel for file of future met.data inputs
      CHARACTER*10 GRIDID			! 20km sqare identifier of last 1km cell
	INTEGER ICELL				! Counter for the cell number in this 20km2 grid
	INTEGER IMON				! Month counter
	INTEGER ITIME				! Timeslice counter (1=2020,2=2050,3=2080)
	INTEGER IYEAR				! Year counter
	REAL RCELL(3,12)			! Total monthly rainfall Jan-Dec 2020 for the 1km2 cell (mm/month)
	REAL RAIN(3,12)				! Total monthly rainfall Jan-Dec 2020 averaged across 20km2 cell (mm/month)
	REAL TCELL(3,12)			! Average monthly temperature Jan-Dec 2020 for the 1km2 cell (deg.C/month)
	REAL TEMP(3,12)				! Average monthly temperature Jan-Dec 2020 averaged across 20km2 cell (deg.C/month)
	REAL THISSQUID				! 1km sqare identifier
	REAL THISEAST				! National grid easting x 1000
	REAL THISNORTH				! National grid northing x 1000
	REAL THISPET(12)			! Temporary holder for PET in this year
	REAL THISTEMP(12)			! Temporary holder for temperature in this year
C
C Variables passed to/from calling subroutine
C
      REAL AVERAIN(12)			! IN:Long term average rainfall(mm/month)
      REAL AVETEMP(12)			! IN:Long term average monthly average 
								!	air temp (deg.C)
      REAL FUTRAIN(MAXMET,12)		! OUT:Future monthly rainfall(mm/month)
      REAL FUTPET(MAXMET,12)		! OUT:Future monthly PET (mm/month)
      REAL FUTTEMP(MAXMET,12)		! OUT:Future monthly average 
								!	air temp (deg.C)
      REAL GRIDLAT				! IN: Average latitude of this 20km2 grid cell
	INTEGER NUMCELLSIN20KM2		! IN:Number of cells in 20km2 cell
C
C Set channel for reading future met data
C
	FUTCHAN=65
C
C Initialise variables for averaging temp and rain across 20km2 cell
C
      DO 100 IMON=1,12
	  DO 200 ITIME=1,3
	    TEMP(ITIME,IMON)=0
200     CONTINUE
100   CONTINUE
C
C Read in future met data for all cells in this 20km2 grid
C
      DO 300 ICELL=1,NUMCELLSIN20KM2
        READ(FUTCHAN,*,END=333)
     &    GRIDID,THISSQUID,THISEAST,THISNORTH,
     &    (TCELL(1,IMON),IMON=1,12),(RCELL(1,IMON),IMON=1,12),
     &    (TCELL(1,IMON),IMON=2,12),(RCELL(2,IMON),IMON=1,12),
     &    (TCELL(1,IMON),IMON=3,12),(RCELL(3,IMON),IMON=1,12)
        DO 400 ITIME=1,3
	    DO 500 IMON=1,12
	      TEMP(ITIME,IMON)=TEMP(ITIME,IMON)+TCELL(ITIME,IMON)
	      RAIN(ITIME,IMON)=RAIN(ITIME,IMON)+RCELL(ITIME,IMON)
500       CONTINUE
400     CONTINUE
300   CONTINUE
      GOTO 303
333   CONTINUE
      PRINT*,'Error in future weather data file.'
      PRINT*,'     ....Press and key to continue'
	READ(*,*)
	STOP
303   CONTINUE
C
C Average data across 20km2 cell
C
      DO 600 ITIME=1,3
	  DO 700 IMON=1,12
	    TEMP(ITIME,IMON)=TEMP(ITIME,IMON)/NUMCELLSIN20KM2
	    RAIN(ITIME,IMON)=RAIN(ITIME,IMON)/NUMCELLSIN20KM2
700     CONTINUE
600   CONTINUE
C
C Interpolate weather data for each year
C
      DO 800 IYEAR=1,30
	  DO 900 IMON=1,12
C
C ... years 1990-2020
C
	    FUTRAIN(IYEAR,IMON)=IYEAR*(RAIN(1,IMON)-AVERAIN(IMON))/30
	    FUTRAIN(IYEAR,IMON)=FUTRAIN(IYEAR,IMON)+AVERAIN(IMON)
	    FUTTEMP(IYEAR,IMON)=IYEAR*(TEMP(1,IMON)-AVETEMP(IMON))/30
	    FUTTEMP(IYEAR,IMON)=FUTTEMP(IYEAR,IMON)+AVETEMP(IMON)
C
C ... years 2021-2050
C
	    FUTRAIN(IYEAR+30,IMON)=IYEAR*(RAIN(2,IMON)-RAIN(1,IMON))/30
	    FUTRAIN(IYEAR+30,IMON)=FUTRAIN(IYEAR+30,IMON)+RAIN(1,IMON)
	    FUTTEMP(IYEAR+30,IMON)=IYEAR*(TEMP(2,IMON)-TEMP(1,IMON))/30
	    FUTTEMP(IYEAR+30,IMON)=FUTTEMP(IYEAR+30,IMON)+TEMP(1,IMON)
C
C ... years 2051-2080
C
	    FUTRAIN(IYEAR+60,IMON)=IYEAR*(RAIN(3,IMON)-RAIN(2,IMON))/30
	    FUTRAIN(IYEAR+60,IMON)=FUTRAIN(IYEAR+60,IMON)+RAIN(2,IMON)
	    FUTTEMP(IYEAR+60,IMON)=IYEAR*(TEMP(3,IMON)-TEMP(2,IMON))/30
	    FUTTEMP(IYEAR+60,IMON)=FUTTEMP(IYEAR+60,IMON)+TEMP(2,IMON)
900     CONTINUE
800   CONTINUE     
C
C Get PET from Thornthwaite equation
C
      DO 1100 IYEAR=1,90
	  DO 1200 IMON=1,12
	    THISTEMP(IMON)=FUTTEMP(IYEAR,IMON)
1200    CONTINUE
	  CALL GET_PET(GRIDLAT,THISTEMP,THISPET)
	  DO 1300 IMON=1,12
         FUTPET(IYEAR,IMON)=THISPET(IMON)
1300    CONTINUE
1100  CONTINUE
      END
C
C-------------------------------------------------------------
C
	SUBROUTINE READINFUTMET_JULES(NUMCELLSIN20KM2,GRIDLAT,
     &                        AVETEMP,AVERAIN,
     &                        FUTNPP,FUTTEMP,FUTRAIN,FUTPET)
C
C Subroutine to read and interpolate in future weather data - ADD NPP INPUT
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
      INTEGER MAXMET				! Number of metyears
	PARAMETER(MAXMET=90)

      INTEGER FUTCHAN				! Channel for file of future met.data inputs
	REAL FNULL
      CHARACTER*10 GRIDID			! 20km sqare identifier of last 1km cell
	INTEGER ICELL				! Counter for the cell number in this 20km2 grid
	INTEGER ILU					! Counter for land uses
	INTEGER IMON				! Month counter
	INTEGER IYEAR				! Year counter
      REAL NCELL(MAXMET,12,MAXLU)	! Monthly NPP Jan-Dec for the 1km2 cell (kg C/ha/month)
      REAL NPP(MAXMET,12,MAXLU)	! Monthly NPP Jan-Dec averaged across 20km2 cell (kg C/ha/month) 
      INTEGER NPPCHAN				! Channel for file of future NPP data inputs
      INTEGER RAINCHAN			! Channel for file of future rain data inputs
      INTEGER TEMPCHAN			! Channel for file of future temp.data inputs

	REAL RCELL(MAXMET,12)		! Total monthly rainfall Jan-Dec for the 1km2 cell (mm/month)
	REAL RAIN(MAXMET,12)		! Total monthly rainfall Jan-Dec averaged across 20km2 cell (mm/month)
	REAL TCELL(MAXMET,12)		! Average monthly temperature Jan-Dec for the 1km2 cell (deg.C/month)
	REAL TEMP(MAXMET,12)		! Average monthly temperature Jan-Dec averaged across 20km2 cell (deg.C/month)
	REAL THISSQUID				! 1km sqare identifier
	REAL THISEAST				! National grid easting x 1000
	REAL THISNORTH				! National grid northing x 1000
	REAL THISPET(12)			! Temporary holder for PET in this year
	REAL THISTEMP(12)			! Temporary holder for temperature in this year
C
C Variables passed to/from calling subroutine
C
      REAL AVERAIN(12)			! IN(CALL):Long term average rainfall(mm/month)
      REAL AVETEMP(12)			! IN(CALL):Long term average monthly average 
								!	air temp (deg.C)
      REAL FUTNPP(MAXMET,12,MAXLU)! OUT(CALL):Future monthly NPP (kg C/ha/month)
      REAL FUTRAIN(MAXMET,12)		! OUT(CALL):Future monthly rainfall(mm/month)
      REAL FUTPET(MAXMET,12)		! OUT(CALL):Future monthly PET (mm/month)
      REAL FUTTEMP(MAXMET,12)		! OUT(CALL):Future monthly average 
								!	air temp (deg.C)
      REAL GRIDLAT				! IN(CALL): Average latitude of this 20km2 grid cell
	INTEGER NUMCELLSIN20KM2		! IN(CALL):Number of cells in 20km2 cell
C
C Set channel for reading future met data
C
	NPPCHAN=55
	TEMPCHAN=56
	RAINCHAN=57
C
C Initialise variables for averaging temp and rain across 20km2 cell
C
      DO 100 IYEAR=1,MAXMET
	  DO 200 IMON=1,12
	    FUTRAIN(IYEAR,IMON)=0
	    FUTTEMP(IYEAR,IMON)=0
	    DO 300 ILU=1,4
	      FUTNPP(IYEAR,IMON,ILU)=0
300       CONTINUE
200     CONTINUE
100   CONTINUE
C
C Read in future met data for all cells in this 20km2 grid
C
      DO 400 ICELL=1,NUMCELLSIN20KM2
        READ(NPPCHAN,*,END=333)
     &    GRIDID,THISSQUID,THISEAST,THISNORTH,
     &    (((NCELL(IYEAR,IMON,ILU),ILU=1,4),IMON=1,12),IYEAR=1,MAXMET)
        READ(TEMPCHAN,*,END=333)
     &    GRIDID,THISSQUID,THISEAST,THISNORTH,
     &    ((TCELL(IYEAR,IMON),IMON=1,12),IYEAR=1,MAXMET)
        READ(RAINCHAN,*,END=333)
     &    GRIDID,THISSQUID,THISEAST,THISNORTH,
     &    ((RCELL(IYEAR,IMON),IMON=1,12),IYEAR=1,MAXMET)
        DO 500 IYEAR=1,MAXMET
	    DO 600 IMON=1,12
            DO 700 ILU=1,4
	        FNULL=FUTNPP(IYEAR,IMON,ILU)+NCELL(IYEAR,IMON,ILU)
		    FUTNPP(IYEAR,IMON,ILU)=FNULL
700         CONTINUE
	      FUTTEMP(IYEAR,IMON)=FUTTEMP(IYEAR,IMON)+TCELL(IYEAR,IMON)
	      FUTRAIN(IYEAR,IMON)=FUTRAIN(IYEAR,IMON)+RCELL(IYEAR,IMON)
600       CONTINUE
500     CONTINUE
400   CONTINUE
      GOTO 303
333   CONTINUE
      PRINT*,'Error in future weather data file.'
      PRINT*,'     ....Press and key to continue'
	READ(*,*)
	STOP
303   CONTINUE
C
C Average data across 20km2 cell
C
      DO 800 IYEAR=1,MAXMET
	  DO 900 IMON=1,12
	    DO 1000 ILU=1,4
	     FUTNPP(IYEAR,IMON,ILU)=FUTNPP(IYEAR,IMON,ILU)/NUMCELLSIN20KM2
1000      CONTINUE
	    FUTTEMP(IYEAR,IMON)=FUTTEMP(IYEAR,IMON)/NUMCELLSIN20KM2
	    FUTRAIN(IYEAR,IMON)=FUTRAIN(IYEAR,IMON)/NUMCELLSIN20KM2
900     CONTINUE
800   CONTINUE
C
C Get PET from Thornthwaite equation
C
      DO 1100 IYEAR=1,MAXMET
	  DO 1200 IMON=1,12
	    THISTEMP(IMON)=FUTTEMP(IYEAR,IMON)
1200    CONTINUE
	  CALL GET_PET(GRIDLAT,THISTEMP,THISPET)
	  DO 1300 IMON=1,12
         FUTPET(IYEAR,IMON)=THISPET(IMON)
1300    CONTINUE
1100  CONTINUE
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE READINGRID(NUMCELLSIN20KM2,KM20GRIDID,GRIDNPP,GRIDLAT,
     &                      FAREA,AVERAIN,AVETEMP,SQUID,EAST,NORTH,
     &                      DOMSOIL,PERSOIL,WETCLASS,FRACLU50,LU1TOLU2,
     &                      ISWAIT,PI_SOURCE)
C
C Subroutine to read in all the cells in the grid square (note the files must be ordered accoridng to 20km2 grid)
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
      INTEGER MAXDEC				! Max.no.of decades included in calculation
	PARAMETER(MAXDEC=9)		
	INTEGER MAXDEC1				! Max.no.of decades included in calculation
	DATA MAXDEC1 /9/
	INTEGER MAXLU				! Max.no.of land use types
      PARAMETER (MAXLU=6)
	INTEGER MAXLU1				! Max.no.of land use types
	DATA MAXLU1 /6/
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXSERIES1
	DATA MAXSERIES1 /5/

	REAL AREA					! Area of cell accounted for (m2)
	REAL CELLLAT				! Latitude of 1km2 grid cell
	REAL CELLNPP				! NPP kgC/m2 x 1000 for this 1km2 grid cell
      INTEGER DATACHAN			! channel for reading in soil&climate data
	REAL DEL_RES				! Null
      CHARACTER*10 GRIDID			! 20km sqare identifier of last 1km cell
	INTEGER ICELL				! Counter for cells within the 20km2 grid
	INTEGER IDEC				! Decade counter
	INTEGER ILU1				! Counter for land use 
	INTEGER ILU2				! Counter for land use 
      INTEGER IMON				! Month counter
	INTEGER IS					! Counter for soil series
	INTEGER ISEND				! Marker indicating that end of file found
	INTEGER ISG					! Counter for soil series in 20km2 grid
	REAL LONG					! Longitude (not used)
      INTEGER LUCHAN				! channel for reading in LU data
      REAL MAIRTEMP(12)			! Long term ave. airtemp for 1km2 grid (deg.C)
	INTEGER NCELL				! Bumber of valid cells read in so far
	INTEGER OD(MAXLU+2)			! Order of land uses in file
C                                     F N G A M S
      DATA (OD(ILU1),ILU1=1,MAXLU+2) /3,4,2,1,5,6,7,8/
	REAL PEROTHER				! percentage of other soil series in the 
								!	1km square cell
      REAL RAIN(12)				! Long term ave. rain for 1km2 cell (mm/month)

	INTEGER STARTCELL			! Cell number in 20km2 grid cell to start read
	REAL TAREA					! Total area of 20km2 cell included (m2)
	REAL TEMP					! temporary real number
C
C Variables passed to/from calling subroutine
C
      REAL AVERAIN(12)			! OUT:Long term average rainfall(mm/month)
      REAL AVETEMP(12)			! OUT:Long term average monthly average 
								!	air temp (deg.C)
      INTEGER DOMSOIL(CELLSTOGRID,MAXSERIES)	! OUT:5 major soil series in the 1km 		
								! square cell x no.cells
	REAL EAST(CELLSTOGRID)		! OUT:National grid easting x 1000
	REAL FAREA(CELLSTOGRID)		! IN:Fraction of 20km2 in this 1km2 cell
	REAL FRACLU50(CELLSTOGRID,MAXLU+2)	! Fraction of 1km cell in LU1 in 1950
	INTEGER GISOK				! IN:Check GIS data 0=WRONG, 1=OK
      REAL GRIDLAT				! OUT:Average latitude of this 20km2 grid cell
	REAL GRIDNPP				! OUT:Average NPP in 20km2 cell (kgC/ha)
	INTEGER ISWAIT				! Code to wait for key press (1) or not (0)
	CHARACTER*10 KM20GRIDID		! OUT:20km sqare identifier
	REAL LU1TOLU2(CELLSTOGRID,MAXDEC,MAXLU+2,MAXLU+2)	! Fracn of LU1 changed 
														! to LU2 in decade
	REAL NORTH(CELLSTOGRID)		! OUT:National grid northing x 1000
	INTEGER NUMCELLSIN20KM2		! OUT:Number of cells in 20km2 cell
	REAL PERSOIL(CELLSTOGRID,MAXSERIES)	! OUT:% of each major soil type in the
								!	1km2 square cell x 20km2 cell
	INTEGER PI_SOURCE			! Source of PI 1=input 2=calc from TOC
	REAL SQUID(CELLSTOGRID)		! OUT:1km sqare identifier
      INTEGER WETCLASS(CELLSTOGRID,MAXSERIES) !OUT:Wetness class for 
								! each major soil type in the 1km2 square cell 
								! x 20km2 cell
C
C Save information for next 20km2 grid
C
      SAVE 
C
C Channels
C
      DATACHAN=30
	LUCHAN=32
C
C Read in all the cells in the grid square (note the files must be ordered accoridng to 20km2 grid)
C ... If first cell has not already been read in, read from cell 1
C
      IF(NUMCELLSIN20KM2.EQ.CELLSTOGRID)THEN
        STARTCELL=1
	  ISEND=0
	  GRIDNPP=0
	  GRIDLAT=0
	  FAREA(1)=0
	  DO 100 IMON=1,12
	    AVERAIN(IMON)=0
	    AVETEMP(IMON)=0
100     CONTINUE
        NUMCELLSIN20KM2=0
C
C ... If first cell in new 20km2 grid square has already been read in, copy data to cell 1
C     and read data from cell 2
C
      ELSEIF(NUMCELLSIN20KM2.LT.CELLSTOGRID)THEN
	  STARTCELL=2
        SQUID(1)=SQUID(NUMCELLSIN20KM2+1)
	  EAST(1)=EAST(NUMCELLSIN20KM2+1)
	  NORTH(1)=NORTH(NUMCELLSIN20KM2+1)
	  FAREA(1)=AREA
        GRIDNPP=CELLNPP
        GRIDLAT=CELLLAT
	  KM20GRIDID=GRIDID
	  DO 200 IS=1,MAXSERIES1
          DOMSOIL(1,IS)=DOMSOIL(NUMCELLSIN20KM2+1,IS)
		PERSOIL(1,IS)=PERSOIL(NUMCELLSIN20KM2+1,IS)
          WETCLASS(1,IS)=WETCLASS(NUMCELLSIN20KM2+1,IS)
200      CONTINUE
        DO 300 IMON=1,12
          AVERAIN(IMON)=RAIN(IMON)/100
	    AVETEMP(IMON)=MAIRTEMP(IMON)/100
300      CONTINUE
        DO 400 ILU1=1,MAXLU1
          FRACLU50(1,ILU1)=FRACLU50(NUMCELLSIN20KM2+1,ILU1)
		DO 500 ILU2=1,MAXLU1
	      DO 600 IDEC=1,MAXDEC1
	        TEMP=LU1TOLU2(NUMCELLSIN20KM2+1,IDEC,ILU1,ILU2)
	        LU1TOLU2(1,IDEC,ILU1,ILU2)=TEMP
600         CONTINUE
500       CONTINUE
400     CONTINUE
	NUMCELLSIN20KM2=1
      ENDIF      
C
C ... Read in up to 400 cells at a time, or until the 20km2 grid cell changes
C
      NCELL=STARTCELL
      DO 700 ICELL=STARTCELL,CELLSTOGRID
        READ(DATACHAN,*,END=777)
     &    GRIDID,SQUID(NCELL),EAST(NCELL),NORTH(NCELL),AREA,
     &    (DOMSOIL(NCELL,IS),PERSOIL(NCELL,IS),
     &     WETCLASS(NCELL,IS),IS=1,MAXSERIES),
     &    PEROTHER,
     &    CELLNPP,(RAIN(IMON),IMON=1,12),(MAIRTEMP(IMON),IMON=1,12),
     &    LONG,CELLLAT
	  IF(NUMCELLSIN20KM2.EQ.0)KM20GRIDID=GRIDID
	  READ(LUCHAN,*,END=777)
     &    GRIDID,SQUID(NCELL),EAST(NCELL),NORTH(NCELL),AREA,
     &    FRACLU50(NCELL,1),DEL_RES,FRACLU50(NCELL,4),DEL_RES,
     &    FRACLU50(NCELL,2),DEL_RES,FRACLU50(NCELL,3),DEL_RES,
     &    DEL_RES,DEL_RES,FRACLU50(NCELL,5),FRACLU50(NCELL,6),
     &    ((LU1TOLU2(NCELL,1,OD(ILU1),OD(ILU2)),ILU1=1,MAXLU1),
     &      ILU2=1,MAXLU1),
     &    ((LU1TOLU2(NCELL,2,OD(ILU1),OD(ILU2)),ILU1=1,MAXLU1),
     &      ILU2=1,MAXLU1),
     &    ((LU1TOLU2(NCELL,3,OD(ILU1),OD(ILU2)),ILU1=1,MAXLU1),
     &      ILU2=1,MAXLU1),
     &    ((LU1TOLU2(NCELL,4,OD(ILU1),OD(ILU2)),ILU1=1,MAXLU1),
     &      ILU2=1,MAXLU1),
     &    ((LU1TOLU2(NCELL,5,OD(ILU1),OD(ILU2)),ILU1=1,MAXLU1),
     &      ILU2=1,MAXLU1),
     &    ((LU1TOLU2(NCELL,6,OD(ILU1),OD(ILU2)),ILU1=1,MAXLU1),
     &      ILU2=1,MAXLU1),
     &    ((LU1TOLU2(NCELL,7,OD(ILU1),OD(ILU2)),ILU1=1,MAXLU1),
     &      ILU2=1,MAXLU1)
C
C Work out proportion of 1km2 that is ILU1 at start (km2/km2)
C
        DO 750 ILU1=1,MAXLU1
	    FRACLU50(NCELL,ILU1)=FRACLU50(NCELL,ILU1)*AREA/1000000
C
C Set LU change in decades 7-9 equivalent to LU change in decade 6
C
          DO 775 ILU2=1,MAXLU1
	!	 ! DO 774 IDEC=7, MAXDEC
		  DO 774 IDEC=8, MAXDEC
			LU1TOLU2(NCELL,IDEC,ILU1,ILU2)=LU1TOLU2(NCELL,7,ILU1,ILU2)
774		  CONTINUE
775       CONTINUE
750     CONTINUE
C
C End of soil cells
C
	  GOTO 708
777     CONTINUE
        IF(ISEND.EQ.0.AND.NUMCELLSIN20KM2.GT.0)THEN
          ISEND=1		! Do simulation for the last 20km2 cell in the file
	    GOTO 708
	  ENDIF
	  CALL CLOSECHAN2()
	  PRINT*,'Simulation of cells completed'
c	  IF(ISWAIT)THEN
c	    PRINT*,'    ...Press any key to continue'
c	    READ(*,*)
c	  ENDIF
	  STOP
708   CONTINUE
C
C Check read in data
C
        CALL CHECK_GIS_CELL(GISOK,NCELL,AREA,DOMSOIL,PERSOIL,WETCLASS,
     &                          CELLNPP,CELLLAT,RAIN,MAIRTEMP,
     &                          LU1TOLU2,PI_SOURCE)
	  IF(GISOK.EQ.0)GOTO 699
	    
        FAREA(NCELL)=AREA
C
C Check that this 1km2 cell is the same 20km2 cell, if not, save data for next run
C
        IF(NUMCELLSIN20KM2.GT.0.AND.GRIDID.NE.KM20GRIDID)GOTO 707	  
C
C Count the number of valid cells in this 20km2 grid
C
	  NUMCELLSIN20KM2=NUMCELLSIN20KM2+1
	  NCELL=NCELL+1
	  IF(NUMCELLSIN20KM2.GT.CELLSTOGRID)THEN
	    PRINT*,'Error in number of cells in 20km2 grid ',GRIDID
	    PRINT*,'                     ....Press any key to continue'
	    READ(*,*)
	    STOP
	  ENDIF
	  KM20GRIDID=GRIDID
C
C Get rainfall and temp (read in as mm x 100 and deg.C x 100)
C Sum values to calculate average rain and air temp for 20km2 grid
C
        DO 800 IMON=1,12
	    AVERAIN(IMON)=AVERAIN(IMON)+(RAIN(IMON)/100)
	    AVETEMP(IMON)=AVETEMP(IMON)+(MAIRTEMP(IMON)/100)
800    CONTINUE
C
C Sum cell NPPs and latitudes to calculate average NPP and lat over 20km grid and 
C
        GRIDNPP=GRIDNPP+CELLNPP
	  GRIDLAT=GRIDLAT+CELLLAT
699     CONTINUE
700   CONTINUE
C
C When all cells in this 20km2 grid read in...
C
707   CONTINUE
C
C ...proportion land use according to the amount of area in the cell
C
      TAREA=0
      DO 900 ICELL=1,NUMCELLSIN20KM2
	  TAREA=TAREA+FAREA(ICELL)
900   CONTINUE
      DO 1000 ICELL=1,NUMCELLSIN20KM2
	  FAREA(ICELL)=FAREA(ICELL)/TAREA
1000   CONTINUE
C
C Convert land use change into km2/km2
C Note original numbers (50s - 70s) are in ha change per decade per 20km2
C      original numbers (80s - 90s) are in kha change per year per 20km2
C We want km2/dec/1km2cell 
C
C For decades 1-3
C LU1TOLU2 ha/dec/20km2 -> LU1TOLU2 (ha / 100ha/km2) / dec / (20km2 / FAREA km2/20km2)
C where FAREA = fraction of the land area in the 20km2 grid cell in this 1km2 cell
C (LU1TOLU2 ha/dec/20km2 x FAREA/100) = LU1TOLU2 in km2/dec/km2
C
C For decades 4-6
C LU1TOLU2 kha/yr/20km2 -> LU1TOLU2 (kha x 10km2/kha) / (yr x 1dec/10yr) / (20km2 / FAREA km2/20km2)
C where FAREA = fraction of the land area in the 20km2 grid cell in this 1km2 cell
C (LU1TOLU2 kha/yr/20km2 x 100xFAREA) = LU1TOLU2 in km2/dec/km2
C
      DO 1100 ICELL=1,NUMCELLSIN20KM2
        DO 1200 IDEC=1,MAXDEC1
	    DO 1300 ILU1=1,MAXLU1+2
	      DO 1400 ILU2=1,MAXLU1+2
	        IF(ILU1.NE.ILU2)THEN
	          TEMP=LU1TOLU2(ICELL,IDEC,ILU1,ILU2)*FAREA(ICELL) ! value/20km2->value/km2
	          LU1TOLU2(ICELL,IDEC,ILU1,ILU2)=TEMP
	          IF(IDEC.LE.3)THEN
	            LU1TOLU2(ICELL,IDEC,ILU1,ILU2)=TEMP/100	! ha/dec->km2/dec
	          ELSEIF(IDEC.GT.3)THEN
	            LU1TOLU2(ICELL,IDEC,ILU1,ILU2)=TEMP*100	! kha/y->km2/dec
	          ENDIF
	        ENDIF
1400        CONTINUE
1300      CONTINUE	  
1200    CONTINUE
1100  CONTINUE  
C
C ...average NPP, latitude, rainfall, temp over all cells in the 20km2 grid
C
      GRIDNPP=GRIDNPP/NUMCELLSIN20KM2
      GRIDLAT=GRIDLAT/NUMCELLSIN20KM2
      DO 1600 IMON=1,12
	  AVERAIN(IMON)=AVERAIN(IMON)/NUMCELLSIN20KM2
	  AVETEMP(IMON)=AVETEMP(IMON)/NUMCELLSIN20KM2
1600  CONTINUE
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE TRAPERROR(TRAPERR,IERROR,NSERIES1,NSERIES2,NSEQ1,NSEQ2)
C
C Subroutine to trap errors in run
C
      IMPLICIT NONE    
C
C Variables local to this subroutine
C
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
	INTEGER MAXSEQ				! Max.no of sequences included in GIS
	PARAMETER (MAXSEQ=MAXLU*MAXLU)						
	INTEGER MAXSERIES			! Max.no of soil series
	PARAMETER (MAXSERIES=5)						

	INTEGER ISERIES						! Counter for soil series
	INTEGER NSEQ						! Counter for number of land use seq 
C
C Variables passed to/from this subroutine
C      
      INTEGER IERROR						! IN:Error type code
	INTEGER NSERIES1					! IN:First soil type with error
	INTEGER NSERIES2					! IN:Last soil type with error
	INTEGER NSEQ1						! IN:First sequences with error
	INTEGER NSEQ2						! IN:Last sequences with error
      INTEGER TRAPERR(CELLSTOGRID*MAXSERIES/2,MAXSEQ) 
								! IN:Error trapping
								!		0=No error
								!		1=Error in NPP data
								!		2=Error in weather data
								!		3=Error in LU data
								!		4=Error in LU change data
								!		5=Soil error
								!		6=Equilibrium not found
C
C Trap errors
C
	DO 100 ISERIES=NSERIES1,NSERIES2
	  DO 200 NSEQ=NSEQ1,NSEQ2
	    TRAPERR(ISERIES,NSEQ)=IERROR
200     CONTINUE
100   CONTINUE
C
C Leave TRAPERROR
C
      END

C*************************************************************
C SOIL C&N ROUTINES
C*************************************************************
C
C EXTERNAL SUBROUTINES
C 1. ADJUST_PI
C 2. CHECK_DOMSOIL
C 3. GET_DOMSOIL_LAYER
C 4. GET_GIS_PI_FROM_ROTHC
C 5. GETSOIL
C 6. INIT_GIS_SOILCN_NOOPT
C 7. INIT_GIS_SOILCN_NOPARS
C 8. USE_WETNESS_CLASS

C
C-------------------------------------------------------------
C
      SUBROUTINE ADJUST_PI_BYNPP(FUTNPP,IK,IYEAR,THISLU,PIANN)
C
C Subroutine to adjust PI according to measured NPP
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLU				! Max.no.of land use types
	PARAMETER (MAXLU=6)
      INTEGER MAXMET				! No.of met years
	PARAMETER (MAXMET=90)
	
	INTEGER IMON				! Month counter
      REAL LASTNPP				! Annual NPP last year
      REAL THISNPP				! Annual NPP this year
C
C Variables passed to/from this subroutine
C
      REAL FUTNPP(MAXMET,12,MAXLU)! IN(CALL):Future monthly NPP (kg/ha/month)
	INTEGER IK					! IN(CALL):No.timesteps from prev.harv. to current
	INTEGER IYEAR				! IN(CALL):Current growing season number - set to 1
      INTEGER THISLU				! IN(CALL):This LU
	REAL PIANN					! IN(CALL)/OUT(CALL): Total annual plant C input 
	                            !      (kg C / ha / yr) 
								!      - used as temprary value to allow adjustment of PI according 
								!        to external factors
C
C Adjust after year 1
C
      IF(IYEAR.GT.1)THEN
C
C Get annual NPP for this year
C
        THISNPP=0
	  LASTNPP=0
        DO 100 IMON=1,12
	    THISNPP=THISNPP+FUTNPP(IYEAR,IMON,THISLU)
	    LASTNPP=LASTNPP+FUTNPP(IYEAR-1,IMON,THISLU)
100     CONTINUE
C
C Adjust plant inputs by the ratio of annual NPP for this year and the previous year
C

        PIANN=PIANN*THISNPP/LASTNPP
	ENDIF
      END
C
C-------------------------------------------------------------
C
	
	SUBROUTINE MIAMI_DYCE(lc,temperature,precipitation,si)
! modification of the miami model by altering coefficients (no need to reparameterise exponent terms since model is effectively linear in climate range of the UK)
! multiply npp values according to land cover type:
! forest (3): *7/8 (Ecology, 89(8), 2008, pp. 2117-2126)
! semi-natural (4), grassland (2), arable (1): forest/2 (Ecology, 89(8), 2008, pp. 2117-2126)
! Miscanthus (5): multiply by 1.6 (from comparison of unadjusted Miami results with Miscanfor peak yield results)
! SRC (6): as forest (Biomass and Bioenergy 32(5), 2008, pp. 407-421 - peak yield is just over half that of Miscanthus)
! for plant inputs to soil, multiply npp values by different fractions according to land cover (Global Change Biology, 16, 2008, pp. 1451-1469; sum of net biome productivities divided by npp)
! forest: 0.15
! semi-natural, grassland: 0.08
! arable: 0.02
! Miscanthus: 0.3 (widely reported as losing around 1/3 of peak mass before harvest; small amount returns to rhizome)
! SRC: 0.15 (assumed as forest)
	implicit none
	real:: temperature		! mean annual temperature
	real:: precipitation	! mean total annual precipitation
	real:: nppt	! temperature-limited npp
	real:: nppp	! precipitation-limited npp
	real:: npp	! minimum of nppt and nppp
	real:: si	! soil input of vegetation
	real:: resc(10)	! npp rescaling factor for each land cover type
	real:: frac(10)	! soil input as fraction of npp for each land cover type
	integer:: lc	! land cover type
!	resc = [0.44,0.44,0.88,0.44,1.6,0.88]  !MLR: commented out
	! [ara, gra, for, nat, mis, src, sug, osr, srf, whe)
	resc = [0.73, 1.0, 0.88, 0.5, 1.6, 1.05, 1.2, 0.55, 1.0, 0.73]   !MLR: altered for ELUM
	frac = [0.53, 0.71, 0.88, 0.71, 0.36, 0.35, 0.25, 0.75, 0.4, 0.53] !MLR: altered for ELUM
      	!resc = [0.44,0.44,0.88,0.7,1.6,0.88] !changed semi-nat 0.44 to 0.7 from Del Grosso et al 2008
     	!resc = [0.8,0.8,0.88,0.6,1.6,0.88]
        !frac = [0.53,0.71,0.88,0.8,0.3,0.323]
        !frac = [0.77,0.77,0.77,0.77,0.3,0.323]

	nppt = 3000/(1+exp(1.315-0.119*temperature))
	nppp = 3000*(1-exp(-0.000664*precipitation))
	npp = 0.5*10*resc(lc)*minval([nppt,nppp],1)      !Times 10 for unit conversion (g/m^2 to Kg/ha) and .5 for C
	si = frac(lc)*npp

	return
	end subroutine miami_dyce

c
c----------------------------------------------------------------------------------------------------------------------
c
      SUBROUTINE ADJUST_PI(AVERAIN,AVETEMP,IK,LHARV,PIANN,RAIN,
     &                          SOILTEMP)
C
C Subroutine to adjust PI according to external factors
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	
	INTEGER IMON				! Current month (rounded to nearest integer)
	REAL RMON					! Current month (real)
	REAL NPPRAIN				! NPP limited by rain (t DM / ha / yr)
	REAL NPPTEMP				! NPP limited by temperature (t DM / ha / yr)
	REAL NPP_LTA				! NPP according to long term average weather data (t DM / ha / yr)
	REAL NPP_MET				! NPP according to current weather data (t DM / ha / yr)
C
C Variables passed to/from this subroutine
C
      REAL AVEPET(12)				! IN:Long term average PET (mm/month)
      REAL AVERAIN(12)			! IN:Long term average rainfall(mm/month)
      REAL AVETEMP(12)			! IN:Long term average monthly average 
	REAL EVAPW					! IN: Potential evap. (mm/timestep)
	INTEGER IK					! IN:No.timesteps from prev.harv. to current
	INTEGER LHARV				! IN:No.timesteps from 01/01 to prev. harvest
	REAL PIANN					! IN/OUT: Total annual plant C input 
	                            !      (kg C / ha / yr) 
								!      - used as temprary value to allow adjustment of PI according 
								!        to external factors
	REAL RAIN					! IN: Rainfall (mm/timestep)
	REAL SOILTEMP(MAXLAYER)		! IN: Soil temperature (deg.C/timestep)
C
C Get current month
C
	IMON=AINT((IK+LHARV)/12.)
	RMON=(REAL(IK+LHARV)/12.)
	IMON=NINT(12.*(RMON-IMON))
	IF(IMON.EQ.0)IMON=12
C
C Work out PI according to the MIAMI model using LTA weather data this month
C
	NPPRAIN=3000*(1-EXP(-0.000664*AVERAIN(IMON)))
	NPPTEMP=3000*1/(1+EXP(1.315-(0.119*AVETEMP(IMON))))
      IF(NPPRAIN.LT.NPPTEMP)THEN
	  NPP_LTA=NPPRAIN
	ELSE
	  NPP_LTA=NPPTEMP
      ENDIF
C
C Work out PI according to the MIAMI model using actual weather data this month
C
	NPPRAIN=3000*(1-EXP(-0.000664*RAIN))
	NPPTEMP=3000*1/(1+EXP(1.315-(0.119*SOILTEMP(1))))
      IF(NPPRAIN.LT.NPPTEMP)THEN
	  NPP_MET=NPPRAIN
	ELSE
	  NPP_MET=NPPTEMP
      ENDIF
C
C Adjust PI according to the ratio of NPP worked out using actual and long term average weather data
C      
      PIANN=PIANN*NPP_MET/NPP_LTA
C
C Leave ADJUST_PI
C
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE CHECK_DOMSOIL(ISERROR,ISMISS,SDEPTH1,SDEPTH0,
     &                         DSC,DSCLAY,DSSILT,DSSAND,DSBD,DSPH,
     &                         DSISIMP,DSIMPDEPTH)
C
C Subroutine to check data of dominant soil types
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/
C
C Variables passed to/from this subroutine
C
      INTEGER ISERROR				! OUT:Code for error in soil 1=error 0=noerror
      INTEGER ISMISS				! OUT:Code for missing values in soil 1=missing 0=nomissing
	REAL DSC					! IN:Soil C 
	REAL DSCLAY					! IN:% Soil clay 
	REAL DSIMPDEPTH				! IN: Depth of impermeable layer if present (cm)
	INTEGER DSISIMP				! IN: Flag for impermeable layer (0=No; 1=Yes)
	REAL DSSILT					! IN:% Soil silt 
	REAL DSSAND					! IN:% Soil sand 
	REAL DSBD					! IN:Soil BD  
	REAL DSPH					! IN:Soil PH 
	REAL SDEPTH1				! IN:Depth of SOM this layer
	REAL SDEPTH0				! IN:Depth of SOM previous layer
C
C Check for missing values
C
      IF(SDEPTH1.LT.-9998.AND.
     &   DSC.LT.-9998.AND.
     &   DSCLAY.LT.-9998.AND.
     &   DSCLAY.LT.-9998.AND.
     &   DSSILT.GT.-9998.AND.
     &   DSSAND.GT.-9998.AND.
     &   DSBD.GT.-9998.AND.
     &   DSPH.GT.-9998)THEN
	     ISMISS=1
	ELSE
	     ISMISS=0
	ENDIF
C
C Check for valid values
C
      IF(SDEPTH1.LT.0.OR.
     &   SDEPTH1.LT.SDEPTH0.OR.
     &   SDEPTH0.LT.0.OR.
     &   DSC.GE.9990000.OR.
     &   DSC.LT.0.OR.
     &   DSCLAY.GT.100.OR.
     &   DSCLAY.LT.0.OR.
     &   DSSILT.GT.100.OR.
     &   DSSILT.LT.0.OR.
     &   DSSAND.GT.100.OR.
     &   DSSAND.LT.0.OR.
     &   DSBD.GT.2.3.OR.
     &   DSBD.LE.0.OR.
     &   DSPH.GT.14.OR.
     &   DSPH.LT.0.OR.
     &   ((DSISIMP.EQ.1).AND.
     &    (DSIMPDEPTH.LT.0.OR.DSIMPDEPTH.GT.MAXDEPTH)))THEN
	     ISERROR=1
	ELSE
	     ISERROR=0
	ENDIF
	END

C
C-------------------------------------------------------------
C
      SUBROUTINE GET_DOMSOIL_LAYER(DSC1,DSC0,
     &                             DSCLAY1,DSCLAY0,
     &                             DSSILT1,DSSILT0,
     &                             DSSAND1,DSSAND0,
     &                             DSBD1,DSBD0,
     &                             DSPH1,DSPH0,
     *                             DEP1,DEP0)
C
C Get dominant soil characteristics from the layer above
C
      IMPLICIT NONE
C
C Variable local to this subroutine
C
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/
C
C Variables passed to/from this subroutine
C
	REAL DSC1					! OUT:Soil C 
	REAL DSC0					! IN:Soil C 
	REAL DSCLAY1				! OUT:% Soil clay 
	REAL DSCLAY0				! IN:% Soil clay 
	REAL DSSILT1				! OUT:% Soil silt 
	REAL DSSILT0				! IN:% Soil silt 
	REAL DSSAND1				! OUT:% Soil sand 
	REAL DSSAND0				! IN:% Soil sand 
	REAL DSBD1					! OUT:Soil BD  
	REAL DSBD0					! IN:Soil BD  
	REAL DSPH1					! OUT:Soil PH 
	REAL DSPH0					! IN:Soil PH 
	REAL DEP1					! IN: Depth 
	REAL DEP0					! IN: Depth 
C
C Get layer 1 from layer 0
C
      DSC1=DSC0*(MAXDEPTH-DEP1)/(DEP1-DEP0)
	DSCLAY1=DSCLAY0
	DSSILT1=DSSILT0				
	DSSAND1=DSSAND0				
	DSBD1=DSBD0					
	DSPH1=DSPH0					
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE GET_GIS_PI_FROM_ROTHC(GRIDNPP,FRACLU50,ISERIES,
     &                      SOMDEPTH,DOMSOILC,DOMSOILCLAY,DOMSOILSILT,
     &                      DOMSOILSAND,DOMSOILBD,DOMSOILPH,
     &                      IAWC,IROCK,NSOIL,
     &                      WTABLE,ROOT,AVEPET,AVETEMP,AVERAIN,
     &                      ITFUNC,IMFUNC,SPARMODEL,TOTPIC,
     &                      NUMCELLSIN20KM2,NSOMLAY)
C
C Get plant inputs to this soil
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
	INTEGER MAXLU1				! Max.no.of land use types
	DATA MAXLU1 /6/
	INTEGER MAXSOIL				! Max.no.of soil types
	PARAMETER (MAXSOIL=50)
	INTEGER MAXSERIES			! Max.no of soil series
	INTEGER CELLSTOGRID					! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
	PARAMETER (MAXSERIES=5)						
	INTEGER MAXSOMLAY			! Max.no.of soil layers entered			
	PARAMETER (MAXSOMLAY=10)
	INTEGER MAXORGM				! Max.no.of organic manure applications allowed
	PARAMETER (MAXORGM=52)
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/

	INTEGER ICELL				! Local cell counter
      INTEGER ILU					! Land use counter
	REAL DEPTH					! Depth
	INTEGER IL					! Layer counter
	INTEGER IMON				! Month counter
	REAL PIALL					! Sum of calculated PI over all LUs
	REAL CALCNPP				! NPP Calculated from proportioned PIs
	REAL FRACLU20KM2(MAXLU)		! Fraction of 20km2 grid in LU1 in 1950
C
C Variables passed to/from calling subroutine
C
      INTEGER ISAVE				! Code to save or retrieve variables
C
C ...Model factors
C
	INTEGER EQMODEL				! OUT:Type of equilibrium run 
	                            ! (NPP, TOC or both)
	INTEGER EQNPP				! Model initialised using plant inputs 
								!    from measured NPP
	INTEGER EQTOC				! Model initialised using measured TOC
	INTEGER EQNPPTOC			! Model initialised using both plant inputs
								!    from measured NPP and measured TOC
	INTEGER EQHILLIER,EQJONES			! Model initialised using TOC and
								!    the Hillier solver to get C pools and plant inputs
	DATA EQNPP,EQTOC,EQNPPTOC,EQHILLIER,EQJONES /1,2,3,4,5/
	INTEGER ITFUNC				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IMFUNC				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
      INTEGER ISERROR				! IN:Code for error in soil 1=error 0=noerror
C
C ...Timing factors
C
	REAL SECONDS				! IN:Number of seconds in one timestep
	REAL CONVER_F
C
C ...Weather factors
C
	REAL BYRAIN					! IN:Rain lost by bypass flow mm/timestep	
	REAL LTA_AWC(12,MAXLAYER)	! IN:Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)	! IN:Average air temp this month (deg.C)
C
C ...Plant factors
C
	REAL GRIDNPP					! IN:Average NPP in the 1km cell (kgC/ha)
	REAL FRACLU50(CELLSTOGRID,MAXLU+2)	! IN:Fraction of 1km cell in LU1 
	INTEGER NUMCELLSIN20KM2				! IN:Number of cells in 20km2 cell
	REAL PI_CEQ_MON(12,MAXLAYER)	! IN:Equilibrium plant C input each month
									!    in each layer (kgC/ha/month/layer)
	REAL TOTPIC(MAXLU)				! OUT:Plant input calculated using 
									!	  ROTHC equilibrium run for each LU (kgC/ha/yr)
C
C ...Soil factors
C
	REAL ALPHA(MAXLAYER)			! IN:Prop.BIO produced on BIO decompn 
	REAL ANIT						! IN:Ammonium immobilised (kgN/ha)
	REAL ANIT15						! IN:Ammonium N15 nitrified (kgN/ha)
      REAL AVEPET(12)					! IN:Long term average PET (mm/month)
      REAL AVERAIN(12)				! IN:Long term average rainfall(mm/month)
      REAL AVETEMP(12)				! IN:Long term average monthly average 
	REAL BETA(MAXLAYER)				! IN:Prop.HUM produced on BIO decompn
	REAL BIOP(MAXSOIL,MAXLAYER)		! IN:BIO/TOC 
	REAL BIORATE(MAXSOIL,MAXLAYER)	! IN:Rate constant for BIO decompn/yr
	REAL BORGM						! IN:C in soil biomass (kgC/ha/layer)
	REAL BPART(MAXSOIL,MAXLAYER)	! IN:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL BPROP(MAXSOIL,MAXLAYER)	! IN:BIO/HUM from biomass decompositn
	REAL BRATE(MAXLAYER)			! IN:Rate constant for HUM decompn
	REAL BULKDENS(MAXLAYER)		! IN:Bulk density of the layer (g/cm3)
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN:% clay content in this layer
      REAL CLAY(MAXLAYER)			! IN:Clay content of the layer (%)
	REAL CRIT(MAXSOIL,MAXLAYER)		! IN:Minimum nitrate level (Sozanska)
	REAL DELTA(MAXLAYER)			! IN:Prop.HUM produced on HUM decompn 
	REAL DFACT(MAXLAYER)			! IN:Denitrification factor
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
									! LU in SOM layers	// (1c) Del
	REAL DORGM						! IN:C in soil DPM (kgC/ha/layer)
      REAL DPM_RPM_RATIO				! OUT:Ratio of DPM:RPM. If set to 
									!     zero, this will be worked out using a call to SETLU
	REAL DPMRATE(MAXLAYER)			! IN:Rate constant for DPM decompn
	REAL FANIT						! IN:Fertiliser N nitrified (kgN/ha)
	REAL FANIT15					! IN:Fertiliser N15 nitrified (kgN/ha)
	REAL FIELDCAP(MAXLAYER)	        ! IN:Soil water content at field capacity (mm/layer)
	REAL FLOWPROP					! IN:Proportion of flow needed to achieve
									!	     water table at depth WTABLE
	REAL GAMMA(MAXLAYER)			! IN:Prop.BIO produced on HUM decompn 
	REAL HORGM						! IN:C in soil humus (kgC/ha/layer)
	REAL HPART(MAXSOIL,MAXLAYER)	! IN:Decomposition efficiency 
									! = (BIO + HUM) / (Total decomposition)
	REAL HPROP(MAXSOIL,MAXLAYER)	! IN:BIO/HUM from humus decompositn
	REAL HRATE(MAXLAYER)			! IN:Rate constant for BIO decompn
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN:Rate constant for HUM decompn/yr
	REAL HY(MAXSOIL,MAXLAYER)		! IN:Stable N:C ratio of BIO&HUM pools
	REAL HZ1(MAXLAYER)				! IN:N:C ratio for steady state
	INTEGER IAWC					! IN:Water movement code number
  	REAL ICFACTOR(MAXLAYER)			! IN:Adjustment in rate needed to	
									!    achieve measured NPP and TOC	
	REAL IOMARRAY(MAXSOIL,MAXLAYER)	! IN:Equilibrium IOM (kgC/ha/layer)
	REAL IOM(MAXLAYER)				! IN:Inert organic C (kgC/ha/layer)
	INTEGER IROCK					! IN:No.layers: 1=50cm, 2=100cm , 3=150cm
	INTEGER ISERIES				    ! IN:Counter for soil series
      INTEGER LUARRAY(MAXSOIL)		! IN/OUT:Land use before equilibrium 
	INTEGER NSOIL					! IN:Soil code number
	INTEGER NSOMLAY(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN:Number of SOM layers
      INTEGER NUMSOIL					! IN:Number of soils defined
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)! IN:pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)! IN:pH above which decomp.max.
	REAL ROOT						! IN:Rooting depth according to restriction (cm)
	REAL RORGM						! IN:C in soil RPM (kgC/ha/layer)
	REAL RPMRATE(MAXLAYER)			! IN:Rate constant for RPM decompn
									! the period between equilibrium and start of simulation. 
	REAL SAND(MAXLAYER)				! IN:Sand content of the layer (%)
	REAL SATWATCONT(MAXLAYER)	    ! IN:Total water content at saturation (mm/layer)
	REAL SILT(MAXLAYER)				! IN:Silt content of the layer (%)
	CHARACTER*40 SNAME(MAXSOIL)		! IN:Soil name
	REAL SOILPH(MAXLAYER)		! IN/OUT:pH of soil in this layer
	REAL SOILW(MAXLAYER)			! IN:Available water (mm/layer)
	REAL SOMDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)		! IN:depth of soil organic matter layers
	INTEGER SPARMODEL				! IN/OUT:Soil parameter model (from file or calc)
	CHARACTER*20 TEMP
									! the period between equilibrium and start of simulation. 
      REAL THISPH						! OUT:pH of soil in this layer
	REAL THISPHP1					! OUT:pH below which decomp.zero
	REAL THISPHP2					! OUT:pH above which decomp.max.
	REAL THISPI(12)				    ! IN/OUT:Equilibrium plant C input each 
								    !        month in each layer (tC/ha/month/layer)
	INTEGER THISTIME   				! OUT:Time between soil at equilibrium 
									! and start of simulation
	INTEGER TIMEARRAY(MAXSOIL)		! IN:Time between soil at equilibrium 
									! and start of simulation
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! IN:Equilibrium TOC (kgC/ha/layer)
	REAL TOC(MAXLAYER)				! IN:Total organic C (kgC/ha/layer)
	REAL WILTPOINT(MAXLAYER)	    ! IN:Water conent at wilting point (mm/layer)
	REAL WMAX(MAXLAYER)				! IN:Avail.water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)				! IN:Avail.water at saturation (mm/layer)
	REAL WTABLE						! IN:Water table depth in cm
C
C ...Manure factors
C
      REAL FPROPX(MAXORGM)			! IN:Prop.of FYM added to top 25cm
	REAL FYMCX(MAXORGM)				! IN:Prop.of C in FYM
	REAL FYMNX(MAXORGM)				! IN:Prop.of Organic N in FYM
	REAL FYMAX(MAXORGM)				! IN:Proportion of ammonium-N in FYM
	REAL FYMWX(MAXORGM)				! IN:Amount of water added in FYM
	REAL FYMLOSS(MAXORGM)			! IN:Prop.FYM lost by volatn/timestep
	REAL FYMPOS(MAXORGM)			! IN:Prop.of FYM added to top 25cm
	REAL FYMSTART(MAXORGM)			! IN:Amount of rainfall in 1 week, 
									!	 below which volatilisation will occur
      real crop(12)                                      
C
C Retrieve soil characteristics
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
!	EQMODEL=EQTOC
C
C For each land use type
C
      DO 100 ILU=1,MAXLU1
C
C Get plant distribution for this land use (assuming plant input is equal to cell NPP)
C
        CALL GET_PLANT_DIST(GRIDNPP,PI_CEQ_MON,ILU)
        CALL GETSOIL(ISERIES,ILU,SOMDEPTH,NSOMLAY,							
     &               DOMSOILC,DOMSOILCLAY,DOMSOILSILT,				
     &               DOMSOILSAND,DOMSOILBD,DOMSOILPH,
     &               DOMSOILISIMP,DOMSOILIMPDEPTH,
     &               TOC,IOM,CLAY,SILT,SAND,BULKDENS,SOILPH,ISERROR)
c	  IF(ISERROR)THEN
c	    TOTPIC(ILU)=0
c		GOTO 100
c	  ENDIF

        CALL GET_PLANT_DIST(TOTPIC(ILU),PI_CEQ_MON,ILU) 
C
C Initialise water
C              
        CALL INIT_GIS_WATER(BULKDENS,CLAY,SAND,FLOWPROP,IAWC,IROCK,
     &                      NSOIL,SILT,TOC,WMAX,WSAT,WTABLE,ROOT,
     &                      AVEPET,AVERAIN,AVETEMP,SOILW,PI_CEQ_MON,
     &                      WILTPOINT,FIELDCAP,SATWATCONT,SECONDS)
C
C Get weather data
C

	DO IMON=1,12
	IF(PI_CEQ_MON(IMON,ILU).GT.0)THEN
	CROP(IMON)=1
	ELSE
	CROP(IMON)=0
	ENDIF
	ENDDO
C
	  CALL GETLTA_LIM(LTA_AWC,LTA_TEMP,WMAX,WSAT,FLOWPROP,  !Gets long term average water, temp etc. uses crop to denote cover,  
     &                AVERAIN,AVEPET,AVETEMP,CROP)            !thus matches water in spinup and the runs.
C
C Get soil C and N parameters
C
	  CALL GET_SUNDIAL_SOILPARS(NSOIL,TOC,IOM,CRIT,ILU,
     &                              CLAY,BULKDENS,SPARMODEL)
C Set the amounts in initial organic matter
C of biomass(=BORGM) and humus (=HORGM)
C and partition down the soil profile
C
        THISPH=7
	  THISPHP1=5.5
        THISPHP2=2.5
	  THISTIME=0
	  TOTPIC(ILU)=0
        DO 200 IL=1,MAXLAYER1
C
C Get depth of this layer
C
          DEPTH=IL*MAXDEPTH/(MAXLAYER1*1.)
C
C Set the amounts in initial organic matter
C of biomass(=BORGM) and humus (=HORGM)
C
            DO 300 IMON=1,12
	        THISPI(IMON)=PI_CEQ_MON(IMON,IL)/1000
300         CONTINUE
            CALL SET_DPMRPMRATIO(ILU,DPM_RPM_RATIO)
            CALL GETEQOM_FROM_ROTHC(TOC(IL),IOM(IL),
     &             BORGM,HORGM,RORGM,DORGM,
     &             CLAY(IL),DEPTH,
     &             THISPH,THISPHP1,THISPHP2,
     &             EQMODEL,THISPI,ICFACTOR(IL),LTA_AWC,LTA_TEMP,
     &		     ITFUNC,IMFUNC,WMAX,WSAT,DPM_RPM_RATIO)
            DO 400 IMON=1,12
	        TOTPIC(ILU)=TOTPIC(ILU)+THISPI(IMON)
400         CONTINUE
200     CONTINUE
100   CONTINUE
C
C Partition GRIDNPP between different land uses
C
      PIALL=0
      DO 500 ILU=1,MAXLU1
        PIALL=PIALL+TOTPIC(ILU)	  
500   CONTINUE
      DO 600 ILU=1,MAXLU1
	  TOTPIC(ILU)=(GRIDNPP*MAXLU1*TOTPIC(ILU))/PIALL
C
C Get average fraction of LU in 20km grid cell
C
        FRACLU20KM2(ILU)=0
        DO 700 ICELL=1,NUMCELLSIN20KM2
          FRACLU20KM2(ILU)=FRACLU20KM2(ILU)+FRACLU50(ICELL,ILU)
700     CONTINUE
        FRACLU20KM2(ILU)=FRACLU20KM2(ILU)/NUMCELLSIN20KM2
600   CONTINUE
C
C Correct plant input so that it gives the measured average NPP across the whole cell
C
      CALCNPP=0
	DO 800 ILU=1,MAXLU1
	  CALCNPP=CALCNPP+FRACLU20KM2(ILU)*TOTPIC(ILU)
800   CONTINUE
      DO 900 ILU=1,MAXLU1
	  TOTPIC(ILU)=TOTPIC(ILU)*GRIDNPP/CALCNPP
900   CONTINUE
C
C Leave GET_PI
C
      END
C
C-------------------------------------------------------------
C
	SUBROUTINE GETSOIL(ISERIES,ILU,SOMDEPTH,NSOMLAY,							
     &                   DOMSOILC,DOMSOILCLAY,DOMSOILSILT,				
     &                   DOMSOILSAND,DOMSOILBD,DOMSOILPH,
     &                   DOMSOILISIMP,DOMSOILIMPDEPTH,
     &                   TOC,IOM,CLAY,SILT,SAND,BULKDENS,SOILPH,ISERROR)
C
C Subroutine to get soil data from 
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLU				! Max.no.of land use types
	INTEGER MAXSOMLAY			! Max.no.of soil layers entered			
	INTEGER MAXLAYER			! No.of layers in the soil profile
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	INTEGER MAXSERIES			! Maximum number of soil series
	PARAMETER (MAXLAYER=60)
	PARAMETER (MAXSOMLAY=10)
	PARAMETER (MAXLU=6)
	INTEGER CELLSTOGRID			! No.of 1km2 cells to each 20km2 grid 
	PARAMETER (CELLSTOGRID=400)
	PARAMETER (MAXSERIES=5)						
	DATA MAXLAYER1 /60/
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/
	INTEGER IL					! Local layer counter
	INTEGER ISOMLAY				! Current SOM layer
	REAL DEPTH					! Depth of current soil layer
	REAL SOMTHICK1				! Thickness of SOM layer (cm)
	REAL LAYRAT1				! Ratio of SOM layer to soil layer
	REAL SOMTHICK0				! Thickness of previous SOM layer (cm)
	REAL LAYRAT0				! Ratio of previous SOM layer to soil layer
	REAL LAYTHICK				! Thickness of soil layer (cm)
	REAL BD						! Bulk density of the layer (g/cm3)
	REAL PH						! pH of soil in this layer
	REAL FTEMP					! Temporary real number
C
C Variables passed to/from calling subroutine
C
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
      INTEGER ISERROR				! OUT:Code for error in soil 1=error 0=noerror
      INTEGER ISMISS				! IN:Code for missing values in soil 1=missing 0=nomissing
	INTEGER ISERIES				! IN:Counter for soil series
	INTEGER ILU					! IN:Counter for land use 
	INTEGER NSOMLAY(CELLSTOGRID*MAXSERIES/2,MAXLU)	! IN:Number of SOM layers
	REAL SOMDEPTH(CELLSTOGRID*MAXSERIES/2,MAXLU,MAXSOMLAY)	! IN:depth of soil organic matter layers
	REAL TOC(MAXLAYER)			! IN:Total organic C (kgC/ha/layer)
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
      REAL CLAY(MAXLAYER)			! IN:Clay content of the layer (%)
	REAL SILT(MAXLAYER)			! IN:Silt content of the layer (%)
	REAL SAND(MAXLAYER)			! IN:Sand content of the layer (%)
	REAL BULKDENS(MAXLAYER)		! IN:Bulk density of the layer (g/cm3)
	REAL SOILPH(MAXLAYER)		! IN/OUT:pH of soil in this layer
C
C Get layer thickness used in simulations
C
      LAYTHICK=(1.*MAXDEPTH/MAXLAYER1)
      ISOMLAY=1
	SOMTHICK1=0
	SOMTHICK0=0
	ISERROR=0
C
C If no layers defined, mark immediately as an error
C
      IF(NSOMLAY(ISERIES,ILU).LE.0)THEN
	  ISERROR=1
        DO 50 IL=1,MAXLAYER1
          TOC(IL)=-9999
	    IOM(IL)=-9999
	    CLAY(IL)=-9999
	    SILT(IL)=-9999
	    SAND(IL)=-9999
	    BULKDENS(IL)=-9999
	    SOILPH(IL)=-9999
50      CONTINUE
C
C If layers are defined, carry on checking
C
      ELSE
C
C For each 5cm layer down the soil profile
C
        DO 100 IL=1,MAXLAYER1
C
C Get depth of this layer
C 
  	    DEPTH=IL*LAYTHICK
C
C Fill bottom of profile, below lowest measurement with values from the previous layers
C
          IF(IL.GT.1.AND.
     &      DEPTH.GT.SOMDEPTH(ISERIES,ILU,NSOMLAY(ISERIES,ILU)))THEN
            TOC(IL)=TOC(IL-1)
            IOM(IL)=IOM(IL-1)
            CLAY(IL)=CLAY(IL-1)
            SILT(IL)=SILT(IL-1)
            SAND(IL)=SAND(IL-1)
            BULKDENS(IL)=BULKDENS(IL-1)
            SOILPH(IL)=SOILPH(IL-1)
            GOTO 101
          ENDIF
C
C Get SOM layer that contains the top of this soil layer
C
          IF(DEPTH.GT.SOMDEPTH(ISERIES,ILU,ISOMLAY))ISOMLAY=ISOMLAY+1
	    IF(ISOMLAY.GT.MAXSOMLAY)ISOMLAY=MAXSOMLAY
C
C Get thickness of this and the previous SOM layer
C
	    IF(ISOMLAY.GT.1)THEN
	      SOMTHICK1=SOMDEPTH(ISERIES,ILU,ISOMLAY)
		  SOMTHICK1=SOMTHICK1-SOMDEPTH(ISERIES,ILU,ISOMLAY-1)
	    ELSEIF(ISOMLAY.EQ.1)THEN
	      SOMTHICK1=SOMDEPTH(ISERIES,ILU,ISOMLAY)
	    ENDIF
	    IF(ISOMLAY.GT.2)THEN
	      SOMTHICK0=SOMDEPTH(ISERIES,ILU,ISOMLAY-1)
		  SOMTHICK0=SOMTHICK0-SOMDEPTH(ISERIES,ILU,ISOMLAY-2)
	    ELSEIF(ISOMLAY.EQ.2)THEN
	      SOMTHICK0=SOMDEPTH(ISERIES,ILU,ISOMLAY-1)
	    ELSEIF(ISOMLAY.EQ.1)THEN
	      SOMTHICK0=0
	    ENDIF
	    IF(SOMDEPTH(ISERIES,ILU,ISOMLAY).EQ.0)SOMTHICK1=0
	    IF(ISOMLAY.GT.1)THEN
		  IF(SOMDEPTH(ISERIES,ILU,ISOMLAY-1).EQ.0)THEN
              SOMTHICK0=0
	      ENDIF
	    ENDIF
C
C Work out the proportions of this and the previous SOM layer that are
C being included in this soil layer
C
C ...If previous SOM layer is below the bottom of this soil layer, 
C    the whole of the soil layer is in this SOM layer and the proportion 
C    of this layer included is calculated as a whole soil later (LAYTHICK)
C
	    IF(SOMDEPTH(ISERIES,ILU,ISOMLAY)-SOMTHICK1.LE.DEPTH-LAYTHICK)THEN
	      LAYRAT1=LAYTHICK
C
C ...otherwise, take off the amount of the soil layer that is in the 
C    previous SOM layer
C
	    ELSE
	      LAYRAT1=SOMDEPTH(ISERIES,ILU,ISOMLAY)-SOMTHICK1
		  LAYRAT1=LAYRAT1-((IL-1)*LAYTHICK)
	      LAYRAT1=(LAYTHICK-LAYRAT1)
	    ENDIF
C
C ...The proportion of the previous soil layer included is obtained from
C    the soil layer thickness minus the amount already included from the 
C    upper layer
C
	    IF(SOMTHICK0.GT.0)THEN
	      LAYRAT0=(LAYTHICK-LAYRAT1)/SOMTHICK0
	    ELSE
	      LAYRAT0=0
	    ENDIF
          IF(SOMTHICK1.GT.0)THEN
	      LAYRAT1=LAYRAT1/SOMTHICK1
	    ELSE
	      LAYRAT1=0
	    ENDIF
C
C Check for invalid soil characteristics - 
C -9999 used to indicate simulation should not be run for this 
C combination of series and LU type
C
	    IF(ISOMLAY.GT.1)THEN
	      CALL CHECK_DOMSOIL(ISERROR,ISMISS,
     &                       SOMDEPTH(ISERIES,ILU,ISOMLAY),
     &                       SOMDEPTH(ISERIES,ILU,ISOMLAY-1),
     &                       DOMSOILC(ISERIES,ILU,ISOMLAY),
     &                       DOMSOILCLAY(ISERIES,ILU,ISOMLAY),
     &                       DOMSOILSILT(ISERIES,ILU,ISOMLAY),
     &                       DOMSOILSAND(ISERIES,ILU,ISOMLAY),
     &                       DOMSOILBD(ISERIES,ILU,ISOMLAY),
     &                       DOMSOILPH(ISERIES,ILU,ISOMLAY),
     &                       DOMSOILISIMP(ISERIES,ILU),
     &                       DOMSOILIMPDEPTH(ISERIES,ILU))
	    ELSE
	      FTEMP=0
	      CALL CHECK_DOMSOIL(ISERROR,ISMISS,
     &                       SOMDEPTH(ISERIES,ILU,ISOMLAY),
     &                       FTEMP,
     &                       DOMSOILC(ISERIES,ILU,ISOMLAY),
     &                       DOMSOILCLAY(ISERIES,ILU,ISOMLAY),
     &                       DOMSOILSILT(ISERIES,ILU,ISOMLAY),
     &                       DOMSOILSAND(ISERIES,ILU,ISOMLAY),
     &                       DOMSOILBD(ISERIES,ILU,ISOMLAY),
     &                       DOMSOILPH(ISERIES,ILU,ISOMLAY),
     &                       DOMSOILISIMP(ISERIES,ILU),
     &                       DOMSOILIMPDEPTH(ISERIES,ILU))
          ENDIF
          IF(ISERROR.eq.1)THEN
             TOC(IL)=-9999
	       IOM(IL)=-9999
	       CLAY(IL)=-9999
	       SILT(IL)=-9999
	       SAND(IL)=-9999
	       BULKDENS(IL)=-9999
	       SOILPH(IL)=-9999
	       GOTO 101
	    ENDIF
          IF(ISOMLAY.GT.2)THEN
            CALL CHECK_DOMSOIL(ISERROR,ISMISS,
     &                       SOMDEPTH(ISERIES,ILU,ISOMLAY-1),
     &                       SOMDEPTH(ISERIES,ILU,ISOMLAY-2),
     &                       DOMSOILC(ISERIES,ILU,ISOMLAY-1),
     &                       DOMSOILCLAY(ISERIES,ILU,ISOMLAY-1),
     &                       DOMSOILSILT(ISERIES,ILU,ISOMLAY-1),
     &                       DOMSOILSAND(ISERIES,ILU,ISOMLAY-1),
     &                       DOMSOILBD(ISERIES,ILU,ISOMLAY-1),
     &                       DOMSOILPH(ISERIES,ILU,ISOMLAY-1),
     &                       DOMSOILISIMP(ISERIES,ILU),
     &                       DOMSOILIMPDEPTH(ISERIES,ILU))
	    ELSEIF(ISOMLAY.GT.1)THEN
	      FTEMP=0
	      CALL CHECK_DOMSOIL(ISERROR,ISMISS,
     &                       SOMDEPTH(ISERIES,ILU,ISOMLAY-1),
     &                       FTEMP,
     &                       DOMSOILC(ISERIES,ILU,ISOMLAY-1),
     &                       DOMSOILCLAY(ISERIES,ILU,ISOMLAY-1),
     &                       DOMSOILSILT(ISERIES,ILU,ISOMLAY-1),
     &                       DOMSOILSAND(ISERIES,ILU,ISOMLAY-1),
     &                       DOMSOILBD(ISERIES,ILU,ISOMLAY-1),
     &                       DOMSOILPH(ISERIES,ILU,ISOMLAY-1),
     &                       DOMSOILISIMP(ISERIES,ILU),
     &                       DOMSOILIMPDEPTH(ISERIES,ILU))
	    ENDIF
          IF(ISERROR.eq.1)THEN
            TOC(IL)=-9999
	      IOM(IL)=-9999
	      CLAY(IL)=-9999
	      SILT(IL)=-9999
	      SAND(IL)=-9999
	      BULKDENS(IL)=-9999
	      SOILPH(IL)=-9999
	      GOTO 101
          ENDIF
C 
C Get the amount of C in this soil layer from the sum of the proportions
C of C included from each SOM layer, and
C %clay, %silt, %sand, bulk density and pH from the averages weighted 
C according to proportion of the value in each SOM layer
C	    
          IF(LAYRAT0+LAYRAT1.GT.0)THEN
            TOC(IL)=LAYRAT1*DOMSOILC(ISERIES,ILU,ISOMLAY)
	      CLAY(IL)=LAYRAT1*DOMSOILCLAY(ISERIES,ILU,ISOMLAY)
	      SILT(IL)=LAYRAT1*DOMSOILSILT(ISERIES,ILU,ISOMLAY)
	      SAND(IL)=LAYRAT1*DOMSOILSAND(ISERIES,ILU,ISOMLAY)
	      BD=LAYRAT1*DOMSOILBD(ISERIES,ILU,ISOMLAY)
	      PH=LAYRAT1*DOMSOILPH(ISERIES,ILU,ISOMLAY)
	      IF(ISOMLAY.GT.1)THEN
	        TOC(IL)=TOC(IL)+(LAYRAT0*DOMSOILC(ISERIES,ILU,ISOMLAY-1))
	        FTEMP=(LAYRAT0*DOMSOILCLAY(ISERIES,ILU,ISOMLAY-1))
	        CLAY(IL)=CLAY(IL)+FTEMP
	        FTEMP=(LAYRAT0*DOMSOILSILT(ISERIES,ILU,ISOMLAY-1))
	        SILT(IL)=SILT(IL)+FTEMP
	        FTEMP=(LAYRAT0*DOMSOILSAND(ISERIES,ILU,ISOMLAY-1))
	        SAND(IL)=SAND(IL)+FTEMP
	        BD=BD+(LAYRAT0*DOMSOILBD(ISERIES,ILU,ISOMLAY-1))
	        PH=PH+(LAYRAT0*DOMSOILPH(ISERIES,ILU,ISOMLAY-1))
	      ENDIF
	      CLAY(IL)=CLAY(IL)/(LAYRAT0+LAYRAT1)
	      SILT(IL)=SILT(IL)/(LAYRAT0+LAYRAT1)
	      SAND(IL)=SAND(IL)/(LAYRAT0+LAYRAT1)
	      BULKDENS(IL)=BD/(LAYRAT0+LAYRAT1)
	      SOILPH(IL)=PH/(LAYRAT0+LAYRAT1)
C
C Derive IOM from Falloon equation
C
            CALL GET_IOM_FROM_FALLOON_EQN(TOC(IL),IOM(IL))
C
C If below the bottom defined SOM layer, set TOC and IOM to 0, and 
C all other parameters to the same value as the previous layer
C
          ELSE
            TOC(IL)=0
            IOM(IL)=0
            IF(IL.GT.1)THEN 
	        CLAY(IL)=CLAY(IL-1)
	        SILT(IL)=SILT(IL-1)
	        SAND(IL)=SAND(IL-1)
	        BULKDENS(IL)=BULKDENS(IL-1)
	        SOILPH(IL)=SOILPH(IL-1)
	      ENDIF
	    ENDIF
101       CONTINUE
100     CONTINUE
      ENDIF
               il=sum(TOC)
C
C Leave GETSOIL	
C
      END

C
C-------------------------------------------------------------
C
      SUBROUTINE INIT_GIS_SOILCN_NOOPT(
C INPUTS: Models, time, environment	
     &                               DOCMODEL,EQMODEL,
     &                               SECONDS,NSOIL,
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
     &                               WMAX,WSAT,ICFACTOR,PH_MODEL,SOILPH)
C
C Subroutine to initialise soil C and N without optimisation 
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
      PARAMETER (MAXFERTSTEPS=30*60*24*7*3)	! Timesteps = 3 wks 
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
	DATA EQNPP,EQTOC,EQNPPTOC,EQHILLIER,EQJONES/1,2,3,4,5/
	INTEGER PH_MODEL			! How is pH calculated?
	INTEGER PH_PARAM			! pH set to neutral or passed from parameter file
      INTEGER PH_STATIC			! pH is read in from input
								!   file and does not change
	INTEGER PH_FROM_VSD			! pH is calculated using VSD (a very 
								!   simple dynamic version of the MAGIC model
								!   by Ed Rowe & Chris Evans, CEH, Bangor)
	DATA PH_PARAM,PH_STATIC,PH_FROM_VSD /0,1,2/			 
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
      INTEGER LUARRAY(MAXSOIL)	! IN/OUT:Land use before equilibrium 
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
	REAL PI_CEQ_MON(12,MAXLAYER)	! IN/OUT: Equilibrium plant C input each 
									!         month to this layer (kgC/ha/month/layer)
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
	REAL AMMN(MAXLAYER)			! IN:Soil ammonium-N (kgN/ha/layer)
	REAL AMMN15(MAXLAYER)		! IN:Soil ammonium-N15 (kgN15/ha/layer)
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
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: % clay content in this layer
	INTEGER ILU						! IN:Counter for land use 
  	REAL ICFACTOR(MAXLAYER)			! IN:Adjustment in rate needed to	
	REAL SOILPH(MAXLAYER)		! IN/OUT:pH of soil in this layer
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
C If pH set as neutral or read from parameter file, pass soil pH out of routine
C
      IF(PH_MODEL.EQ.PH_PARAM)THEN
        DO 100 IL=1,MAXLAYER
	    SOILPH(IL)=PHARRAY(NSOIL,IL)
100     CONTINUE
C
C Otherwise, save pH to here
C
      ELSEIF(PH_MODEL.EQ.PH_STATIC.OR.PH_MODEL.EQ.PH_FROM_VSD)THEN
        DO 200 IL=1,MAXLAYER
	    PHARRAY(NSOIL,IL)=SOILPH(IL)
200     CONTINUE
      ENDIF
C
C Set time factors according to SECONDS
C
      CALL SETTIME_CN(SECONDS,CONVER_F)
C
C Set BIO and HUM decompostion parameters
C
c      CALL OMPROP(NSOIL,HY,HZ1,ALPHA,BETA,GAMMA,DELTA,
c     &                  BPART,HPART,BPROP,HPROP,PHARRAY)
C
C Set residual DOC content
C
      IF(DOCMODEL.EQ.DOC_ON)CALL SETDOC(NSOIL,MOBDOC,MOBDON)
C
C Set initial fertiliser additions
C
      CALL SET_GIS_FERT(FERTADD,FERTADD15)
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
      SUBROUTINE INIT_GIS_SOILCN_NOPARS(
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
     &                               MOBDOC,MOBDON,SATWATCONT,
     &                               CLAY,BULKDENS,PI_CEQ_MON,
     &                               LTA_AWC,LTA_TEMP,ITFUNC,IMFUNC,
     &                               WMAX,WSAT,ICFACTOR,wiltpoint,
     &                               DPM_RPM_RATIO,DPMCTON,RPMCTON,
     &                               PH_MODEL,SOILPH,
     &                               CACCUM,ANCH4,ANCO2,ANDOC)
C
C Subroutine to initialise soil C and N without reading in parameters
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
      PARAMETER (MAXFERTSTEPS=30*60*24*7*3)	! Timesteps = 3 wks 
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
	INTEGER PH_PARAM			! pH set to neutral or passed from parameter file
      INTEGER PH_STATIC			! pH is read in from input
								!   file and does not change
	INTEGER PH_FROM_VSD			! pH is calculated using VSD (a very 
								!   simple dynamic version of the MAGIC model
								!   by Ed Rowe & Chris Evans, CEH, Bangor)
	DATA PH_PARAM,PH_STATIC,PH_FROM_VSD /0,1,2/			 
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
	INTEGER ICMODEL				! IN:Type of C initialisation model 
	INTEGER ICFIXED				! Initialisation of C is fixed
	INTEGER ICROTHCEQ			! Initialisation of C by RothC 
								!  equilibrium run
	DATA ICFIXED,ICROTHCEQ /1,2/
	INTEGER INMODEL				! IN:Type of N initialisation model
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
	INTEGER PH_MODEL			! IN:How is pH calculated?
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
      REAL ANCH4					! OUT: Measured annual CH4 emissions
								!     kgC/ha/yr
	REAL ANCO2					! OUT: Measured annual CO2 emissions
								!     kgC/ha/yr
	REAL ANDOC					! OUT: Measured annual DOC loss
								!     kgC/ha/yr
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
      REAL CACCUM					! OUT: Measured annual C accumulation 
								!     in the soil kgC/ha/yr
      REAL CLAY(MAXLAYER)			! IN:Clay content of the layer (%)
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: % clay content in this layer
	REAL CRIT(MAXSOIL,MAXLAYER)		! IN/OUT: Minimum nitrate level (Sozanska)
	REAL DELTA(MAXLAYER)			! IN/OUT: Prop.HUM produced on HUM decompn 
	REAL DFACT(MAXLAYER)			! IN/OUT: Denitrification factor
      REAL DPM_RPM_RATIO			! IN/OUT:Ratio of DPM:RPM. If set to 
								!    zero, this will be worked out using a call to SET_DPMRPMRATIO
      REAL DPMCTON(MAXLAYER)		! IN:Ratio of C:N in DPM
	REAL DPMCARB0(MAXLAYER)		! IN:C in decomposable PM (kgC/ha/layer)
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
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for HUM decompn/yr
	REAL HRATE(MAXLAYER)			! IN/OUT: Rate constant for BIO decompn
	REAL HY(MAXSOIL,MAXLAYER)	! IN/OUT:Stable N:C ratio of BIO & HUM pools
	REAL HZ1(MAXLAYER)				! IN/OUT: N:C ratio for steady state
  	REAL ICFACTOR(MAXLAYER)			! IN:Adjustment in rate needed to	
									!    achieve measured NPP and TOC	
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
	REAL IOMARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: Equilibrium IOM (kgC/ha/layer)
      INTEGER LUARRAY(MAXSOIL)	! IN/OUT:Land use before equilibrium 
	INTEGER NSOIL				! IN:Soil code number
      INTEGER NUMSOIL				! IN/OUT:Number of soils defined
	REAL PHARRAY(MAXSOIL,MAXLAYER)	! IN/OUT:pH of soil in this layer
	REAL PHP1ARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: pH below which decomp.zero
	REAL PHP1(MAXLAYER)			! IN/OUT: pH below which decomp.zero
	REAL PHP2ARRAY(MAXSOIL,MAXLAYER)	! IN/OUT: pH above which decomp.max.
	REAL PHP2(MAXLAYER)			! IN/OUT: pH above which decomp.max.
	REAL PI_CEQ_MON(12,MAXLAYER)	! IN/OUT: Equilibrium plant C input each 
									!         month to this layer (kgC/ha/month/layer)
	REAL RPMCARB0(MAXLAYER)		! IN:C in resistant PM (kgC/ha/layer)
      REAL RPMCTON(MAXLAYER)		! IN:Ratio of C:N in RPM
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
      real satwatcont(MAXLAYER)
      real wiltpoint(maxlayer)
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
C If pH set as neutral or read from parameter file, pass soil pH out of routine
C
      IF(PH_MODEL.EQ.PH_PARAM)THEN
        DO 100 IL=1,MAXLAYER
	    SOILPH(IL)=PHARRAY(NSOIL,IL)
100     CONTINUE
C
C Otherwise, save pH to here
C
      ELSEIF(PH_MODEL.EQ.PH_STATIC.OR.PH_MODEL.EQ.PH_FROM_VSD)THEN
        DO 200 IL=1,MAXLAYER
	    PHARRAY(NSOIL,IL)=SOILPH(IL)
200     CONTINUE
      ENDIF
C
C Set time factors according to SECONDS
C
      CALL SETTIME_CN(SECONDS,CONVER_F)
C
C Set starting value of organic matter
C
 !     CALL STARTORG(NSOIL,TOC,TOCARRAY,IOM,IOMARRAY)
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
        CALL SETORGM_GIS_EQRUN(NSOIL,DPMRATE,RPMRATE,TOC,IOM,CONVER_F,
     &                   BRATE,BIORATE,HRATE,HUMRATE,HZ1,HY,
     &                   BCARB0,HCARB0,DPMCARB0,RPMCARB0,
     &                   BNIT0,HNIT0,DPMNIT0,RPMNIT0,
     &                   HNLAB0,BNLAB0,DPMNLAB0,RPMNLAB0,
     &                   CLARRAY,wiltpoint,
     &                   PHARRAY,PHP1ARRAY,PHP2ARRAY,
     &                   EQMODEL,PI_CEQ_MON,ICFACTOR,LTA_AWC,LTA_TEMP,
     &                   ITFUNC,IMFUNC,WMAX,WSAT,satwatcont,
     &                   DPM_RPM_RATIO,DPMCTON,RPMCTON,ALPHA,BETA,CLAY)	
      ELSEIF(ICMODEL.EQ.ICFIXED)THEN
        CALL SETORGM_GIS_FIXED(NSOIL,DPMRATE,RPMRATE,TOC,IOM,CONVER_F,
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
      CALL SET_GIS_FERT(FERTADD,FERTADD15)
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
C-----------------------------------------------------------------------
C
      SUBROUTINE USE_WETNESS_CLASS(WETCLASS)
C
C Subroutine to use the soil wetness class supplied
C
      IMPLICIT NONE
C
C Variables local to this routine
C
      
C
C Variables passed to/from this subroutine
C
      INTEGER WETCLASS				! Soil wetness class
C***********************************************
C Add code here to use soil wetness class
C***********************************************

C
C  Leave USE_WETNESS_CLASS
C
      END
C
C-----------------------------------------------------------------------
C
C INTERNAL SUBROUTINES
C 1. GET_CTON_GIS
C 2. SETORGM_GIS_ACCUM
C 3. SETORGM_GIS_EQRUN
C 4. SETORGM_GIS_FIXED
C
C------------------------------------------------------------
C
      SUBROUTINE GET_CTON_GIS(DPMCTON,RPMCTON,LUCODE)
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
	DATA (DCTON(ILU),ILU=1,MAXLU) /14.3,60.8,51.6,60.8,60.8,51.6/
	REAL RCTON(MAXLU)			! C:N ratio of RPM
C                                    Ara  Gra  For  Nat  Mis SRC
	DATA (RCTON(ILU),ILU=1,MAXLU) /41.9,97.9,79.7,97.9,97.9,79.7/
C
C Variables passed to/from this subroutine
C
      REAL DPMCTON(MAXLAYER)	! OUT:DPM C:N ratio
	INTEGER LUCODE		! IN:Land use code 1=Ara;2=Gra;3=For;4=Natural
      REAL RPMCTON(MAXLAYER)	! OUT:RPM C:N ratio
C
C Get C:N ratio for this LU from array
C
      DO 100 IL=1,MAXLAYER
        DPMCTON(IL)=DCTON(LUCODE)
	  RPMCTON(IL)=RCTON(LUCODE)
100   CONTINUE
C
C Leave GET_CTON_GIS
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE SETORGM_GIS_ACCUM(
C Inputs: Soil code, time conversion factor,... 
     &                   NSOIL,CONVER_F,
C      ...DPM:RPM, C:N of PI SOM pools,...
     &                   DPM_RPM_RATIO,DPMCTON,RPMCTON,HZ1,HY,
C      ...proportion of biomass (alpha) & humus (beta) produced on decomposition
     &                   ALPHA,BETA,
C      ...measured carbon, inert OM, crop cover, LU...
     &                   TOC,IOM,LUCODE,
C      ...annual C accumulation, CH4 emission, CO2 emission and DOC emission
     &                   CACCUM,ANCH4,ANCO2,ANDOC,
C      ...bulk density, pH & pH factors
     &                   BULKDENS,PH,PHP1,PHP2,
C      ...longterm average temp & rate mod.codes, soil water & rate mod.codes, 
C         field cap and sat
     &                   LTA_TEMP,ITFUNC,LTA_AWC,IMFUNC,WMAX,WSAT,
C Outputs: Decomposition rate constants,... 
     &                   DPMRATE,RPMRATE,BRATE,BIORATE,HRATE,HUMRATE,
     &                   ICFACTOR,
C      ...pool sizes (C, N and N15)
     &                   DPMCARB0,RPMCARB0,BCARB0,HCARB0,
     &                   DPMNIT0,RPMNIT0,BNIT0,HNIT0,
     &                   DPMNLAB0,RPMNLAB0,HNLAB0,BNLAB0) 
C
C Subroutine to set the residual N and C in the RO,BIO and HUM pools
C for accumulating peats
C NOTE: STILL UNDER DEVELOPMENT!
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	INTEGER MAXDEPTH			! Maximum depth of the soil profile
	PARAMETER (MAXLAYER=60)								
	DATA MAXDEPTH /300/
	INTEGER MAXSOIL				! Max.no.of soil types
      PARAMETER (MAXSOIL=50)

      REAL COMRAT					! C:OM ratio (assumed 0.58)
	DATA COMRAT /0.58/

      INTEGER ACT_LAY				! Active layers
	REAL AVEALPHA				! Average prop.bio produced on decomp in active peat layer
      REAL AVEBD					! Average bulk density of active peat layer (g/cm3)
	REAL AVEBETA				! Average prop.hum produced on decomp in active peat layer
	REAL AVEBRATE				! Ave.rate constant for BIO decomposition
	REAL AVEDRATE				! Ave.rate constant for DPM decomposition
	REAL AVEHRATE				! Ave.rate constant for HUM decomposition
	REAL AVERRATE				! Ave.rate constant for RPM decomposition
      REAL CGAIN_BIO				! C gain as BIO over time required to grow active peat layer (kg C / ha)
      REAL CGAIN_HUM				! C gain as HUM over time required to grow active peat layer (kg C / ha)
      REAL CLOSS_BIO				! C loss (as gas, bio and hum) from BIO over time required to grow active peat layer (kg C / ha)
      REAL CLOSS_BIOHUM			! C loss (as gas, bio and hum) from BIO & HUM over time required to grow active peat layer (kg C / ha)
      REAL CLOSS_DPM				! C loss (as gas, bio and hum) from DPM over time required to grow active peat layer (kg C / ha)
      REAL CLOSS_HUM   			! C loss (as gas, bio and hum) from HUM over time required to grow active peat layer (kg C / ha)
      REAL CLOSS_RPM				! C loss (as gas, bio and hum) from RPM over time required to grow active peat layer (kg C / ha)
      REAL CLOSS_RPMBIOHUM		! C loss (as gas, bio and hum) from RPM,BIO & HUM over time required to grow active peat layer (kg C / ha)
      REAL CLOSS_RPMHUM			! C loss (as gas, bio and hum) from RPM & HUM over time required to grow active peat layer (kg C / ha)
 	REAL DEPTHPI				! Maximum depth of plant input (cm)
								! (assumed to be the actively growing peat layer)
      REAL FTEMP					! Temprary real used in calculations
	REAL GASCLOSS				! Gaseous C loss over time required to grow active peat layer (kg C / ha)
	REAL GASCLOSS_DPM			! Gaseous C loss from DPM over time required to grow active peat layer (kg C / ha)
	REAL GASCLOSS_RPM			! Gaseous C loss from RPM over time required to grow active peat layer (kg C / ha)
      REAL HUM_TO_RPMHUM			! Ratio of losses of HUM to (RPM+HUM)
	INTEGER IL					! Local counter for layers
	INTEGER IMON				! Local counter for months
      REAL PI_ANN					! Annual plant input C (kgC/ha/year)
	REAL PI_DPM					! Total plant input to DPM during growth of the active layer (kg C / ha)
	REAL PI_GROW				! Plant input needed to grow the active layer (kg C / ha)
      REAL PI_LAY					! Total plant input in this layer (kgC/ha/layer)
	REAL PI_RPM					! Total plant input to RPM during growth of the active layer (kg C / ha)
	REAL RATEMOD				! Average rate modifier across year for active layer
      REAL RPM_TO_RPMHUM			! Ratio of losses of RPM to (RPM+HUM)
	REAL T_GROW					! Time required to grow the active peat layer (years)
C
C Variables passed to/from calling subroutine
C
	REAL ALPHA(MAXLAYER)			! IN/OUT: Prop.BIO produced on BIO decompn 
      REAL ANCH4					! IN: Measured annual CH4 emissions
								!     kgC/ha/yr
	REAL ANCO2					! IN: Measured annual CO2 emissions
								!     kgC/ha/yr
	REAL ANDOC					! IN: Measured annual DOC loss
								!     kgC/ha/yr
  	REAL BCARB0(MAXLAYER)		! OUT:C in soil biomass (kgC/ha/layer)
	REAL BETA(MAXLAYER)				! IN/OUT: Prop.HUM produced on BIO decompn
	REAL BIORATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for BIO decompn/yr
      REAL BNIT0(MAXLAYER)		! OUT:N in soil biomass (kgN/ha/layer)
	REAL BNLAB0(MAXLAYER)		! OUT:N15 in soil humus (kgN15/ha/layer)
      REAL CACCUM					! IN: Measured annual C accumulation 
								!     in the soil kgC/ha/yr
      REAL CONVER_F				! IN:Conversion to correct timestep
	INTEGER NSOIL				! IN:Soil code number
	REAL BRATE(MAXLAYER)		! IN/OUT: Rate constant for HUM decompn
	REAL BULKDENS(MAXLAYER)		! IN:Bulk density of the layer (g/cm3)
      REAL CRRATE					! IN:crop modifying factor (proportion)
      REAL DPM_RPM_RATIO			! IN:Ratio of DPM:RPM. If set to 
								!     zero, this will be worked out using a call to SET_DPMRPMRATIO
	REAL DPMCARB0(MAXLAYER)		! OUT:C in decomposable PM (kgC/ha/layer)
      REAL DPMCTON(MAXLAYER)		! IN:Ratio of C:N in DPM
	REAL DPMNIT0(MAXLAYER)		! OUT:N in decomposable PM (kgN/ha/layer)
	REAL DPMNLAB0(MAXLAYER)		! OUT:N15 in decomposable PM (kgN15/ha/layer)
	REAL DPMRATE(MAXLAYER)		! IN/OUT: Rate constant for DPM decompn
	REAL HCARB0(MAXLAYER)		! OUT:C in soil humus (kgC/ha/layer)
	REAL HNIT0(MAXLAYER)		! OUT:N in soil humus (kgN/ha/layer)
      REAL HNLAB0(MAXLAYER)		! OUT:N15 in soil biomass (kgN15/ha/layer)
	REAL HRATE(MAXLAYER)		! IN/OUT: Rate constant for BIO decompn
	REAL HUMRATE(MAXSOIL,MAXLAYER)	! IN/OUT: Rate constant for HUM decompn/yr
	REAL HY(MAXSOIL,MAXLAYER)	! IN:Stable N:C ratio of BIO & HUM pools
	REAL HZ1(MAXLAYER)			! IN: N:C ratio for steady state
	REAL ICFACTOR(MAXLAYER)		! IN:Adjustment in rate needed to 
								!    achieve measured NPP and TOC
      INTEGER ICOVER				! OUT:Code for crop cover: 1=covered; 2=bare
	INTEGER IMFUNC				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
	REAL IOM(MAXLAYER)			! IN:Inert organic C (kgC/ha/layer)
	INTEGER ITFUNC				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER LUCODE				! IN:Code for land use type
	REAL LTA_AWC(12,MAXLAYER)	! IN/OUT:Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)	! IN/OUT:Average air temp this month (deg.C)

	REAL PH(MAXLAYER)			! IN/OUT:pH of soil in this layer
	REAL PHP1(MAXLAYER)			! IN/OUT: pH below which decomp.zero
	REAL PHP2(MAXLAYER)			! IN/OUT: pH above which decomp.max.
      REAL PHRATE					! IN:pH rate modifier (proportion)
	REAL PI_CEQ_MON(12,MAXLAYER)	! IN:Equilibrium plant C input each 
									!         month in each layer (kgC/ha/month/layer)
	REAL RPMCARB0(MAXLAYER)		! OUT:C in resistant PM (kgC/ha/layer)
      REAL RPMCTON(MAXLAYER)		! IN:Ratio of C:N in RPM
	REAL RPMNIT0(MAXLAYER)		! OUT:N in resistant PM (kgN/ha/layer)
	REAL RPMNLAB0(MAXLAYER)		! OUT:N15 in resistant PM (kgN15/ha/layer)
	REAL RPMRATE(MAXLAYER)		! IN/OUT: Rate constant for RPM decompn
	REAL TOC(MAXLAYER)			! IN:Total organic C (kgC/ha/layer)
	REAL TRATE					! IN:Temperature rate modifier summed for layers
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WRATED					! IN:Moisture rate modifyer for DPM & BIO(prop)
	REAL WRATER					! IN:Moisture rate modifyer for RPM & HUM(prop)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
C
C Set Rate constants for decomposition of RO (=R), BIO (=B) and HUM (=H)
C
	CALL GET_DECOMP_RATECONSTS(CONVER_F,NSOIL,BIORATE,HUMRATE,
     &                           DPMRATE,RPMRATE,BRATE,HRATE)	  
C
C Work out depth of plant input (assumed to be the active peat layer),
C average bulk density, proportions of bio and hum produced on decomp., 
C and average rate constants
C
      PI_ANN=CACCUM+ANCH4+ANCO2+ANDOC
      CALL GET_PLANT_DIST(PI_ANN,PI_CEQ_MON,LUCODE)
	AVEBD=0
	AVEALPHA=0
	AVEBETA=0
	AVEBRATE=0
	AVEHRATE=0
	AVEDRATE=0
	AVERRATE=0
	RATEMOD=0
	DO 100 IL=1,MAXLAYER
	  PI_LAY=0
	  DO 200 IMON=1,12
	    PI_LAY=PI_LAY+PI_CEQ_MON(IMON,IL)
200     CONTINUE             
        IF(PI_LAY.LE.0)DEPTHPI=(IL-1)*MAXDEPTH/MAXLAYER
	  IF(PI_LAY.GT.0)DEPTHPI=IL*MAXDEPTH/MAXLAYER
	  AVEBD=AVEBD+BULKDENS(IL)
	  AVEALPHA=AVEALPHA+ALPHA(IL)
	  AVEBETA=AVEBETA+BETA(IL)
	  AVEBRATE=AVEBRATE+BRATE(IL)
	  AVEDRATE=AVEDRATE+DPMRATE(IL)
	  AVEHRATE=AVEHRATE+HRATE(IL)
	  AVERRATE=AVERRATE+RPMRATE(IL)
100   CONTINUE
      ACT_LAY=(PI_LAY*MAXLAYER/MAXDEPTH)
      IF(DEPTHPI.LT.MAXDEPTH)THEN
	  AVEBD=AVEBD-BULKDENS(MAXLAYER)
	  AVEALPHA=AVEALPHA-ALPHA(MAXLAYER)
	  AVEBETA=AVEBETA-BETA(MAXLAYER)
	  AVEBRATE=AVEBRATE-BRATE(MAXLAYER)
	  AVEHRATE=AVEHRATE-HRATE(MAXLAYER)
	  AVEDRATE=AVEDRATE-DPMRATE(MAXLAYER)
	  AVERRATE=AVERRATE-RPMRATE(MAXLAYER)	
      ENDIF
	AVEBD=AVEBD/ACT_LAY
	AVEALPHA=AVEALPHA/ACT_LAY
	AVEBETA=AVEBETA/ACT_LAY
	AVEBRATE=AVEBRATE/ACT_LAY
	AVEHRATE=AVEHRATE/ACT_LAY
	AVEDRATE=AVEDRATE/ACT_LAY
	AVERRATE=AVERRATE/ACT_LAY
C
C Work out average rate modifying factor across all months and over active layer
C
      ICOVER=1
      DO 300 IL=1,ACT_LAY
	  DO 400 IMON=1,12
	    CALL MODFACTS_MINER(LTA_AWC(IMON,IL),WMAX(IL),WSAT(IL),
     &                WRATED,WRATER,
     &                LTA_TEMP(IMON,IL),TRATE,								
     &			    PH(IL),PHP1(IL),PHP2(IL),PHRATE,
     &                ICOVER,CRRATE,ITFUNC,IMFUNC)									
	    RATEMOD=RATEMOD+((WRATED+WRATER)/2)*TRATE*PHRATE*CRRATE
400     CONTINUE
300   CONTINUE
	RATEMOD=RATEMOD/(12*ACT_LAY)
C
C Get time required to grow the active layer of peat from rate of C accumulation, 
C depth of the active layer and the average bulk density of the active layer
C
      T_GROW=(COMRAT*100000*DEPTHPI*AVEBD)/CACCUM
C
C Get plant input needed to grow the active peat layer
C 
      PI_GROW=T_GROW*PI_ANN
C
C Get inputs to DPM and RPM from DPM:RPM ratio
C
      PI_DPM=(DPM_RPM_RATIO*PI_ANN)/(1+DPM_RPM_RATIO)
	PI_RPM=PI_ANN/(1+DPM_RPM_RATIO)
C
C Get gaseous losses over time required to grow peat from annual emissions x time
C
      GASCLOSS=T_GROW*(ANCH4+ANCO2)
C
C Partition gaseous C losses into losses from DPM and RPM
C
      FTEMP=PI_DPM*EXP(AVEDRATE*RATEMOD)+PI_RPM*EXP(AVERRATE*RATEMOD)
      GASCLOSS_DPM=GASCLOSS*(PI_DPM*EXP(AVEDRATE*RATEMOD))/FTEMP
      GASCLOSS_RPM=GASCLOSS*(PI_RPM*EXP(AVERRATE*RATEMOD))/FTEMP
C
C Work total C losses from DPM from (CO2 loss from DPM)*(1+((BIO+HUM)/CO2))
C
      CLOSS_DPM=1+((AVEALPHA+AVEBETA)/(1-AVEALPHA-AVEBETA))
	CLOSS_DPM=GASCLOSS_DPM*CLOSS_DPM
C
C Limit DPM losses to plant input to DPM pool
C
      IF(CLOSS_DPM.GT.PI_DPM)THEN
	  CLOSS_DPM=PI_DPM
	ENDIF
C
C Work out losses from RPM
C
	GASCLOSS_DPM=(1-AVEALPHA-AVEBETA)*CLOSS_DPM
	CLOSS_RPMBIOHUM=(GASCLOSS-GASCLOSS_DPM)/(1-AVEALPHA-AVEBETA)
      CLOSS_RPM=CLOSS_RPMBIOHUM-((AVEALPHA+AVEBETA)*CLOSS_DPM)
	CLOSS_RPM=CLOSS_RPM/(1+AVEALPHA+AVEBETA)
C
C Limit RPM losses to plant input to RPM pool
C
101   CONTINUE
      IF(CLOSS_RPM.GT.PI_RPM)THEN
	  CLOSS_RPM=PI_RPM
	ENDIF
C
C Work out losses from BIO and HUM
C
	GASCLOSS_RPM=(1-AVEALPHA-AVEBETA)*CLOSS_RPM
      CLOSS_BIOHUM=(GASCLOSS-GASCLOSS_DPM-GASCLOSS_RPM)
	CLOSS_BIOHUM=CLOSS_BIOHUM/(1-AVEALPHA-AVEBETA)
	CLOSS_BIO=(AVEALPHA/(AVEALPHA+AVEBETA))*CLOSS_BIOHUM
	CLOSS_HUM=(AVEBETA/(AVEALPHA+AVEBETA))*CLOSS_BIOHUM
C
C Work out gains in biomass and humus 
C	
      CGAIN_BIO=AVEALPHA*(CLOSS_DPM+CLOSS_RPM+CLOSS_BIOHUM)        
      CGAIN_HUM=AVEBETA*(CLOSS_DPM+CLOSS_RPM+CLOSS_BIOHUM)        
C
C Limit losses from biomass to the amount gained
C	
      IF(CLOSS_BIO.GT.CGAIN_BIO)THEN
	  CLOSS_BIO=CGAIN_BIO
	  CLOSS_RPMHUM=CLOSS_RPMBIOHUM-CLOSS_BIO
	  RPM_TO_RPMHUM=CLOSS_RPM/(CLOSS_RPM+CLOSS_HUM)
	  HUM_TO_RPMHUM=CLOSS_RPM/(CLOSS_RPM+CLOSS_HUM)
	  CLOSS_RPM=RPM_TO_RPMHUM*CLOSS_RPMHUM
	  CLOSS_HUM=HUM_TO_RPMHUM*CLOSS_RPMHUM
	  GOTO 101
	ENDIF
C
C Work out pool sizes 
C
      DO 500 IL=1,MAXLAYER
	  DPMCARB0(IL)=(PI_DPM-CLOSS_DPM)/ACT_LAY
	  RPMCARB0(IL)=(PI_RPM-CLOSS_RPM)/ACT_LAY
	  BCARB0(IL)=(CGAIN_BIO-CLOSS_BIO)/ACT_LAY
	  HCARB0(IL)=(CGAIN_HUM-CLOSS_HUM)/ACT_LAY
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
500   CONTINUE
C
C Leave SETORGM_GIS_ACCUM
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE SETORGM_GIS_EQRUN(NSOIL,DPMRATE,RPMRATE,TOC,IOM,
     &                   CONVER_F,BRATE,BIORATE,HRATE,HUMRATE,HZ1,HY,
     &                   BCARB0,HCARB0,DPMCARB0,RPMCARB0,
     &                   BNIT0,HNIT0,DPMNIT0,RPMNIT0,
     &                   HNLAB0,BNLAB0,DPMNLAB0,RPMNLAB0,
     &                   CLARRAY,wiltpoint,
     &                   PHARRAY,PHP1ARRAY,PHP2ARRAY,
     &                   EQMODEL,PI_CEQ_MON,ICFACTOR,LTA_AWC,LTA_TEMP,
     &                   ITFUNC,IMFUNC,WMAX,WSAT,satwatcont,
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
	REAL DEPTH					! Depth of this layer
	INTEGER IL					! Local counter for layers
	INTEGER IMON				! Local counter for months
	REAL THISPI(12)			    ! Equilibrium plant C input each month 
								! to this layer (t C / ha / month)
C
C Variables passed to/from calling subroutine
C ...Model descriptors
C
	INTEGER EQMODEL				! IN:Type of equilibrium run (NPP,TOC or both)
	INTEGER EQNPP				! Model initialised using plant inputs 
								!    from measured NPP
	INTEGER EQTOC				! Model initialised using measured TOC
	INTEGER EQNPPTOC			! Model initialised using both plant inputs
								!    from measured NPP and measured TOC
	INTEGER EQHILLIER,EQJONES			! Model initialised using TOC and
								!    the Hillier solver to get C pools and plant inputs
	DATA EQNPP,EQTOC,EQNPPTOC,EQHILLIER,EQJONES /1,2,3,4,5/

	INTEGER ITFUNC				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IMFUNC				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)     
C
C ...Weather factors
C
	REAL AWC_IL(12)				! OUT:Long term ave.soil moist.def.(mm)
	REAL TEMP_IL(12)			! OUT:Long term ave.soil moist.def.(mm)   
	REAL LTA_AWC(12,MAXLAYER)	! IN:Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)	! IN:Average air temp this month (deg.C)
C
C ...Soil descriptors
C
	REAL ALPHA(MAXLAYER)		! IN: Prop.BIO produced on BIO decompn 
	REAL BETA(MAXLAYER)			! IN: Prop.HUM produced on BIO decompn
      REAL CLAY(MAXLAYER)			! IN:Clay content of the layer (%)
	REAL BORGM					! IN:C in soil biomass (kgC/ha/layer)
	REAL HORGM					! IN:C in soil humus (kgC/ha/layer)
	REAL RORGM					! IN:C in soil RPM (kgC/ha/layer)
	REAL DORGM					! IN:C in soil DPM (kgC/ha/layer)
	INTEGER NSOIL				! IN:Soil code number
      REAL DPM_RPM_RATIO			! IN/OUT:Ratio of DPM:RPM. If set to 
								!    zero, this will be worked out using a call to SET_DPMRPMRATIO
      REAL DPMCTON(MAXLAYER)		! IN:Ratio of C:N in DPM
      REAL RPMCTON(MAXLAYER)		! IN:Ratio of C:N in RPM
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
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
	real satwatcont(maxlayer)
      real wiltpoint(maxlayer)
	REAL PI_CEQ_MON(12,MAXLAYER)! IN/OUT: Equilibrium plant C input each month
	REAL PI_IL(12)              ! IN/OUT: Equilibrium plant C input each month
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
C
C Model initialised using an analytical solution of the RothC equilibrium run
C worked out by Jon Hillier ! Edit by EOJ to adjust for standard decomposition rate
C
	  ELSEIF(EQMODEL.EQ.EQHILLIER)THEN
          DO 250 IMON=1,12
		  AWC_IL(IMON)=LTA_AWC(IMON,IL)
	      TEMP_IL(IMON)=LTA_TEMP(IMON,IL)
		  PI_IL(IMON)=PI_CEQ_MON(IMON,IL)
250       CONTINUE
		DPMRATE=52*DPMRATE/CONVER_F
		RPMRATE=52*RPMRATE/CONVER_F
		BRATE=52*BRATE/CONVER_F
		HRATE=52*HRATE/CONVER_F
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
     &                   BORGM,DORGM,HORGM,RORGM,il)

		DPMRATE=DPMRATE*CONVER_F/52
		RPMRATE=RPMRATE*CONVER_F/52
		BRATE=BRATE*CONVER_F/52
		HRATE=HRATE*CONVER_F/52
		
		
          DO 260 IMON=1,12
	      THISPI(IMON)=PI_IL(IMON)/1000
260       CONTINUE
	
	ELSEIF(EQMODEL.EQ.EQJONES)THEN
      DO 251 IMON=1,12
		  AWC_IL(IMON)=LTA_AWC(IMON,IL)
	      TEMP_IL(IMON)=LTA_TEMP(IMON,IL)
		  PI_IL(IMON)=PI_CEQ_MON(IMON,IL)
251       CONTINUE


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
     &                   BORGM,DORGM,HORGM,RORGM,IL,
     &                    wiltpoint(il),satwatcont(il))

	DPMRATE=DPMRATE*CONVER_F/52 !convert back
	RPMRATE=RPMRATE*CONVER_F/52
	BRATE=BRATE*CONVER_F/52
	HRATE=HRATE*CONVER_F/52

          DO 261 IMON=1,12
	      THISPI(IMON)=PI_IL(IMON)/1000
261       CONTINUE
	  ENDIF
C
C Set the amounts in initial organic matter
C of biomass(=BORGM) and humus (=HORGM)
C
        DO 300 IMON=1,12
	    PI_CEQ_MON(IMON,IL)=THISPI(IMON)*1000 
300     CONTINUE
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
C Leave SETORGM_GIS_EQRUN
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE SETORGM_GIS_FIXED(NSOIL,DPMRATE,RPMRATE,TOC,IOM,
     &                   CONVER_F,BRATE,BIORATE,HRATE,HUMRATE,HZ1,HY,
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
	    DPM_PLUS_RPM=(DPM_OM_FIXED+RPM_OM_FIXED)*ACT_OM_FIXED
	    RPMCARB0(IL)=DPM_PLUS_RPM/(1+DPM_RPM_RATIO)
	    DPMCARB0(IL)=DPM_PLUS_RPM-RPMCARB0(IL)
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
C Leave SETORGM_FIXED
C
      RETURN
      END

C*************************************************************
C WATER ROUTINES
C*************************************************************
C
C EXTERNAL SUBROUTINES
C 1. INIT_GIS_WATER
C-------------------------------------------------------------
      SUBROUTINE INIT_GIS_WATER(BULKDENS,CLAY,SAND,FLOWPROP,IAWC,IROCK,
     &                          NSOIL,SILT,TOC,WMAX,WSAT,WTABLE,ROOT,
     &                          AVEPET,AVERAIN,AVETEMP,SOILW,PI_CEQ_MON,
     &                          WILTPOINT,FIELDCAP,SATWATCONT,SECONDS)
C
C Subroutine to initialise soil water 
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
	DATA MAXLAYER1 /60/
	INTEGER IL				! Local layer counter
	REAL WILTPOINT_LYR(MAXSOIL,MAXLAYER) ! Water conent at wilting point (mm/layer)
C
C Variables passed to/from calling subroutine
C ...Soil factors
C
	REAL BULKDENS(MAXLAYER)	! IN:Bulk density of the layer (g/cm3)
      REAL CLAY(MAXLAYER)		! IN:Clay content of the layer (%)
	REAL FLOWPROP			! IN/OUT:Proportion of flow needed to achieve
							!	     water table at depth WTABLE
	INTEGER IAWC			! IN:Water movement code number
	INTEGER IROCK			! IN:No.layers: 1=50cm, 2=100cm , 3=150cm
      INTEGER ISAVE			! OUT:Code to save or retrieve variables
	REAL MAXWAT(MAXSOIL,MAXLAYER)	! IN/OUT: Maximum water content (mm/layer)
      INTEGER NUMSOIL			! IN/OUT: Number of soils defined
	INTEGER NSOIL			! IN:Soil code number
	REAL SAND(MAXLAYER)		! IN:Sand content of the layer (%)
	REAL SILT(MAXLAYER)		! IN:Silt content of the layer (%)
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)
	CHARACTER*40 SNAME(MAXSOIL)		! IN/OUT: Soil name
	REAL TOC(MAXLAYER)		! IN:Total organic C (kgC/ha/layer)
	REAL WATSAT(MAXSOIL,MAXLAYER)	! IN/OUT:Avail.water at saturation (mm/layer)
	REAL WMAX(MAXLAYER)		! IN/OUT:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)		! IN/OUT:Available water at saturation (mm/layer)
	REAL WTABLE				! IN:Water table depth in cm
	REAL FIELDCAP(MAXLAYER)	! OUT:soil water content at field capacity (mm/layer)
	REAL WILTPOINT(MAXLAYER)	! OUT:Water conent at wilting point (mm/layer)
	REAL SATWATCONT(MAXLAYER)	! OUT:Total water content at saturation (mm/layer)
C
C ... Timing factors
C
	REAL SECONDS			! IN:Number of seconds in one timestep
C
C ... Crop factors
C
	REAL ROOT				! IN:Rooting depth according to restriction (cm)
C
C ...Weather factors
C
	REAL AVEPET(12)			! IN/OUT: Long term average PET (mm)
	REAL AVERAIN(12)		! IN/OUT: Long term average rainfall (mm)
	REAL AVETEMP(12)		! IN/OUT: Long term average temperature (deg.C)
      REAL PI_CEQ_MON(12,MAXLAYER)
C
C Set parameters and save for future use
C
	CALL MAKE_PAR_WAT(NSOIL,MAXWAT,TOC,CLAY,SAND,SILT,BULKDENS,WATSAT,
     &                  WILTPOINT_LYR)
	ISAVE=1
	CALL SAVE_WAT(NUMSOIL,MAXWAT,SNAME,ISAVE,WATSAT)
C
C Set the available water at field capacity and saturation
C
      CALL SWATER(IAWC,WMAX,MAXWAT,WSAT,WATSAT,FIELDCAP,
     &		    WILTPOINT_LYR,WILTPOINT,SATWATCONT)
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
C Set soil water to field capacity
C
C >>> Commented out by Mark Richards, 05/03/2012
C Commented out so that the initial soil water conditions are at the
C equilbrium achieved by the spin-up carried out in GETEQWATER()
C      DO 100,IL=1,MAXLAYER1
C        SOILW(IL)=WMAX(IL)
C100   CONTINUE
C <<< End of commented out code by Mark Richards, 05/03/2012 
C
C Leave INIT_WAT
C
      END
!
!--------------------------------------------------------------------------------

      subroutine adj_pi_by_lu_age(iyear, ik, lu_code, last_lu,
     &                            lu_change_yr, pi_ann)
      
      ! Adjust annual plant inputs due to vegetation not being at equilibrium 
      ! after a land use change event.
      !
      ! The plant inputs are currently assumed to increase linearly from the
      ! year of land use change until equilibrium is reached. This should be
      ! an exponential or sigmoidal curve.
      
      implicit none
      
	! Arguments with intent(in)
      integer, intent(in) :: ik
      integer, intent(in) :: iyear
      integer, intent(in) :: lu_code

      ! Arguments with intent(inout)      
      real, intent(inout) :: pi_ann
      integer, intent(inout) :: last_lu
      integer,intent(inout) :: lu_change_yr
      
      ! Local variables
      real :: pi_frac         ! Fraction of equilibrium plant input added
      integer :: yrs_since_change
      
      ! Local constants
      real, parameter, dimension(10) :: yrs_to_equil =          ! Years to reach equilibrium plant input 
     &                        (/1, 1, 20, 20, 6, 5, 1, 1, 15, 1/)
      real, parameter, dimension(10) :: pi_reduction =          ! Max reduction to non-equliibrium PI after LU change
     &              (/0.0, 0.0, 0.4, 0.2, 0.0, 0.5, 0.0, 0.0, 0.4, 0.0/)

c      if (ik > nsow) then
	if (iyear > 1 .and. lu_code /= last_lu) then
            lu_change_yr = iyear
      endif
      last_lu = lu_code
      yrs_since_change = iyear - lu_change_yr
          
      if (lu_change_yr > 0) then
	    if (yrs_since_change < yrs_to_equil(lu_code)) then
	        pi_frac = 1.0 - (yrs_since_change / 
     &                         yrs_to_equil(lu_code))
	        pi_frac = 1.0 - (pi_reduction(lu_code) * pi_frac)
	    else
	        pi_frac = 1
          endif
          if (pi_frac < 0) pi_frac = 0
          pi_ann = pi_ann * pi_frac  
      endif
c      endif
      
      end subroutine  ! adjust_pi_by_lu_age()

