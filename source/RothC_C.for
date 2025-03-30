C ****************************************************************************
C
C ROTHC-26.3 : A MODEL FOR THE TURNOVER OF ORGANIC CARBON IN THE SOIL
C July  1997
C K.W. Coleman, D.S. Jenkinson, L.C. Parry and J.H. Rayner
C IACR - Rothamsted, Harpenden, Herts. AL5 2JQ. U.K.
C Copyright IACR - Rothamsted
C Modified for application in ATEAM J.U. Smith - May 2003
C Modified for application in ECOSSE J.U. Smith - Feb 2005
C 
C ****************************************************************************
      SUBROUTINE BPCALC(CLAYPER,FDEC1,FPDEC1,FDEC2,FPDEC2)
C
C
C  To produce CO2 : (BIO + POM) ratio from percent clay value
C
C
      IMPLICIT NONE
C
C Variables local to this routine
C
      REAL CLAYF		! Clay factor
	REAL SOILF		! Soil factor
	REAL X			! Temp value used in calc
	REAL XP			! Temp value used in calc
	INTEGER I		! Pool counter
C
C Variables passed to/from this routine
C
	REAL CLAYPER	! % clay
	REAL FDEC1(5)	! Prop. of plant material, and humus that goes to 
					! pool I: I=1=DPM; I=2=RPM; I=3=BIO; I=4=BIO; I=5=HUM 
	REAL FPDEC1(5)	! Prop. of plant material, and humus that goes to 
					! pool I: I=1=DPM; I=2=RPM; I=3=BIO; I=4=BIO; I=5=HUM 
	REAL FDEC2(5)	! Prop. of humus (check) that goes to pool 
					! I: I=1=DPM; I=2=RPM; I=3=BIO; I=4=BIO; I=5=HUM 
					! after loss of CO2
	REAL FPDEC2(5)	! Prop. of decomposing material that goes to pool 
					! I: I=1=DPM; I=2=RPM; I=3=BIO; I=4=BIO; I=5=HUM 
					! after loss of CO2
C
      CLAYF=1.67
      SOILF=1.00
c NOTE: If BPCALC is changed also change BPART and HPART in MAKE_PAR_CN
      CALL GET_DECOMP_EFF(CLAYPER,X)
	X=(1/X)-1
C      X=CLAYF*(1.85+1.60*EXP(-0.0786*CLAYPER))
      XP=X/SOILF
      DO 10 I=1,5
        FDEC2(I)=0.0
        FPDEC2(I)=0.0
   10 CONTINUE
C
      FDEC2(3)=FDEC1(3)/(X+1)
      FDEC2(5)=FDEC1(5)/(X+1)
C
      FPDEC2(4)=FPDEC1(4)/(XP+1)
      FPDEC2(5)=FPDEC1(5)/(XP+1)
C
      RETURN
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE CALC1(SOIL,SOILB,DEC,DECOMP,MCEND,MCSTART)
C
C
C  Primary calculations on input data
C
C
      IMPLICIT NONE
C
C Variables local to this routine
C
	INTEGER I
C
C Variables passed to/from this routine
C
	REAL SOIL(0:6),SOILB(0:6),DEC(5),DECOMP(5)
	INTEGER MCEND,MCSTART
C
C  Calculating initial monthly decomposition rates
C
      DO 10 I=1,6
        SOILB(I)=SOIL(I)
   10 CONTINUE
C
      DO 20 I=1,5
        DEC(I)=DECOMP(I)/12.0
   20 CONTINUE
C
C Find 'last' month of the year
C
      MCEND=MCSTART-1
      IF(MCEND.LT.1)MCEND=12
C
      RETURN
      END

C
C--------------------------------------------------------------
C
      SUBROUTINE GETEQOM_FROM_ROTHC(THISTOC,THISIOM,
     &                   BORGM,HORGM,RORGM,DORGM,
     &                   CLAYPER,DEPTH,
     &				   PH,PHP1,PHP2,
     &                   EQMODEL,THISPI,ICFACT,LTA_AWC,LTA_TEMP,
     &				   ITFUNC,IMFUNC,WMAX,WSAT,DPM_RPM_RATIO)
C
C Subroutine to get OM pools at equilibrium
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
	INTEGER NITER				! Number of iterations
	INTEGER MAXITER				! Maximum allowed number of iterations
	INTEGER M					! Month counter
	INTEGER IMON				! Month counter
	INTEGER IPOOL				! Pool counter
	INTEGER MAXLU				! Maximum number of land use types
	PARAMETER (MAXLU=6)	
	INTEGER MAXSOMLAY			! Maximum number of SOM layers
	PARAMETER (MAXSOMLAY=10)
	REAL PL_MON					! Total standard plant input in a month
	REAL TOTPI					! Total plant input kgC/ha/yr
C
C Variables passed to/from this routine
C ... Model descriptors
C
	INTEGER EQMODEL				! OUT:Type of equil.run (NPP,TOC or both)
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
C ...Soil factors
C
      REAL DPM_RPM_RATIO			! IN:Ratio of DPM:RPM. If set to 
								!    zero, this will be worked out using a call to SET_DPMRPMRATIO
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
      REAL THISTOC				! IN:Total organic C (kgC/ha/layer)
	REAL THISIOM				! IN:Inert organic C (kgC/ha/layer)
	REAL BORGM					! IN:C in soil biomass (kgC/ha/layer)
	REAL HORGM					! IN:C in soil humus (kgC/ha/layer)
	REAL RORGM					! IN:C in soil RPM (kgC/ha/layer)
	REAL DORGM					! IN:C in soil DPM (kgC/ha/layer)
      REAL DRRAT					! Ratio of DPM:RPM
	REAL TPAR
	REAL TPPAR
	REAL TFPAR
      REAL SOIL(0:6)				! Soil C 1=DPM;2=RPM;3=BIOa;4=BIOz;5=HUM;6=IOM 
								! (t/ha)
      REAL SOIL1(0:6)				! Soil C 1=DPM;2=RPM;3=BIOa;4=BIOz;5=HUM;6=IOM 
								! (t/ha)
	REAL FDEC1(5)
	REAL FPDEC1(5)
	REAL DECOMP(5)				
	REAL CLAYPER				! % Clay
	REAL DEPTH					! Depth
	REAL RATEM(12)				! Rate modifier per month
	REAL TOTCSIM				! Simulated C input (t/ha)
	REAL TOTCMEAS				! Measured C input (t/ha)
	REAL DEC(5)
	REAL FDEC2(5)
	REAL FPDEC2(5)
	REAL PH						! IN:pH of this soil layer
	REAL PHP1					! IN: pH below which decomp.zero
	REAL PHP2					! IN: pH above which decomp.max.
C
C ...Plant factors
C
	REAL FPLANT(5)				! Proportion plant material going to each pool
	INTEGER ICROP(12)			! Crop cover (0=no;1-yes)
	REAL THISPI(12)				! IN/OUT: Equilibrium plant C input each month
								!         in each layer (tC/ha/month/layer)
  	REAL ICFACT					! IN:Adjustment in rate needed to	
								!    achieve measured NPP and TOC	
C
C ...Manure factors
C
	REAL FFYM(5)				! Proportion fym going to each pool
	REAL FYMADD(12)				! FYM input (t/ha/month)
C
C ...Weather factors
C
	REAL LTA_AWC(12,MAXLAYER)	! IN/OUT:Long term ave.soil moist.def.(mm)
	REAL LTA_TEMP(12,MAXLAYER)	! IN/OUT:Long term ave.soil moist.def.(mm)
C
C ...Timing factors
C
	INTEGER ISYEAR				! Starting year of ROTHC simulation
	INTEGER NCYEARS				! Number of years for ROTHC run
	INTEGER MCSTART				! Number of starting month for ROTHC run
C
C Get crop cover from plant input
C Assumption: If there is any plant input, crop cover =1; otherwise assume no crop cover
C
      DO 100 IMON=1,12
        IF(THISPI(IMON).GT.0)THEN
	    ICROP(IMON)=1
	  ELSEIF(THISPI(IMON).LE.0)THEN
	    ICROP(IMON)=0
        ENDIF
100   CONTINUE
C
C Set IOM
C
      SOIL(6)=THISIOM/1000
C
C Set Measured TOC
C
      TOTCMEAS=THISTOC/1000
	IF(TOTCMEAS.LE.0.0001.AND.EQMODEL.NE.EQNPP)THEN
	  BORGM=0
	  HORGM=0
        DORGM=0
	  RORGM=0
	  RETURN
	ENDIF
C
C Setup land use and soil characteristics
C
	DRRAT=DPM_RPM_RATIO
C
C Set up simulation for steady state phase calculation
C
	CALL SETEQ(ISYEAR,NCYEARS,MCSTART,SOIL,DRRAT,
     &                 FPLANT,FFYM,FDEC1,FPDEC1,DECOMP,
     &                 TPAR,TPPAR,TFPAR)
C
C Equilibrium run
C
      ICFACT=1														
      CALL RUNTOEQ(TOTCSIM,SOIL,MCSTART,CLAYPER,DEPTH,
     &                   ICROP,RATEM,DEC,DECOMP,
     &                   FDEC1,FPDEC1,FDEC2,FPDEC2,
     &                   FPLANT,THISPI,FFYM,FYMADD,NCYEARS,
     &                   PH,PHP1,PHP2,ICFACT,LTA_AWC,LTA_TEMP,
     &                   ITFUNC,IMFUNC,WMAX,WSAT)
C																		
C For equilibrium runs using TOC or NPP&TOC, optimise by comparison to TOC
C
      IF(EQMODEL.EQ.EQTOC.OR.EQMODEL.EQ.EQNPPTOC)THEN
C
C Match equilibrium phase of the calculation to measurement
C
        NITER=0
	  MAXITER=10000
        DO 300 IPOOL=0,6
	    SOIL1(IPOOL)=SOIL(IPOOL)
300     CONTINUE
202	  CONTINUE
        NITER=NITER+1
C
C If measured tot C is non-zero...
C and calculated tot C does not match measured tot C,...
C
C	  IF(TOTCMEAS.GE.0.0001.AND.CLAYPER.GE.0.0001.AND.
	  IF(TOTCMEAS.GE.0.0001.AND.
     &  (TOTCSIM.GT.(1.0001*TOTCMEAS).OR.TOTCSIM.LT.(0.9999*TOTCMEAS))
     &                                                             )THEN
	    DO 400 IPOOL=0,6
	      SOIL(IPOOL)=SOIL1(IPOOL)
400       CONTINUE
C
C ...adjust inputs and... 
C
          IF(EQMODEL.EQ.EQTOC)THEN
            CALL GETPI(THISPI,SOIL,TOTCMEAS,TOTCSIM)
C
C ...for equilibrium run using TOC&NPP, adjust internal consistency factor... 
C
          ELSEIF(EQMODEL.EQ.EQNPPTOC)THEN	
	      IF(ICFACT.LT.0.0001)THEN
  	        ICFACT=ICFACT+(TOTCSIM-SOIL(6))/(TOTCMEAS-SOIL(6))
	      ELSE
		    ICFACT=ICFACT*(TOTCSIM-SOIL(6))/(TOTCMEAS-SOIL(6))
	      ENDIF
		ENDIF													
C      
C ...repeat the calculation
C
	    DRRAT=DPM_RPM_RATIO
	    CALL RUNTOEQ(TOTCSIM,SOIL,MCSTART,CLAYPER,DEPTH,
     &                   ICROP,RATEM,DEC,DECOMP,
     &                   FDEC1,FPDEC1,FDEC2,FPDEC2,
     &                   FPLANT,THISPI,FFYM,FYMADD,NCYEARS,
     &                   PH,PHP1,PHP2,ICFACT,LTA_AWC,LTA_TEMP,
     &                   ITFUNC,IMFUNC,WMAX,WSAT)
C
C Check for zero plant inputs
C
          TOTPI=0
          DO 200 IMON=1,12
	      TOTPI=TOTPI+THISPI(IMON)
200       CONTINUE
C
C Loop back to test the fit with measured data
C
	    IF(NITER.LT.MAXITER.AND.TOTPI.GE.0)THEN
	      GOTO 202
	     ELSE
	      !TOTCSIM=-999
          !BORGM=-999
	      !HORGM=-999
	      !DORGM=-999
	      !RORGM=-999
	  print*,'Warning: Simulated C and Measured C differ:',
     & TOTCSIM/TOTCMEAS
            RETURN
	    ENDIF          	        
	  ENDIF																
      ENDIF
C
C Set BORGM and HORGM from equilibrium run
C
      BORGM=1000*(SOIL(3)+SOIL(4))
	HORGM=1000*(SOIL(5))
	DORGM=1000*(SOIL(1))
	RORGM=1000*(SOIL(2))
	END
C
C--------------------------------------------------------------
C
      SUBROUTINE GETMET(RAIN,PET,TEMP)														
C
C Set this years weather data
C
C
C Variables local to this routine
C
      INTEGER I
C
C Variables passed to / from this routine
C
      REAL RAIN(12),PET(12),TEMP(12)
C
C Read in weather data from RAINFILE, PETFILE and TEMPFILE
C
      READ(12,*,ERR=111)(RAIN(I),I=1,12)
      READ(13,*,ERR=111)(PET(I),I=1,12)
      READ(14,*,ERR=111)(TEMP(I),I=1,12)
	GOTO 101
 111  CONTINUE
      PRINT*,'Error in weather files. Press any key to continue'
	READ(*,*)
	STOP
 101  CONTINUE
      END        
C
C--------------------------------------------------------------
C
      SUBROUTINE GETPI(THISPI,SOIL,TOTCMEAS,TOTCSIM)
C
C Estimates plant inputs
C
      IMPLICIT NONE
C
C Variables local to this routine
C
      REAL TOTPI					! Total plant input over the year (kg C/ha/yr)
	REAL TOTNEED				! Total plant input required to achieve 
								! measured SOM (kg C/ha/yr
	INTEGER IMON				! Local month counter
C
C Variables passed to/from this routine
C
	REAL THISPI(12)				! IN/OUT: Equilibrium plant C input each month
								!         in each layer (tC/ha/month/layer)
      REAL SOIL(0:6)				! Soil C 1=DPM;2=RPM;3=BIOa;4=BIOz;5=HUM;
	                            ! 6=IOM (t/ha)
	REAL TOTCSIM				! Simulated C input (t/ha)
	REAL TOTCMEAS				! Measured C input (t/ha)
C
C Calculate annual inputs
C
      TOTPI=0
      DO 100 IMON=1,12
	  TOTPI=TOTPI+THISPI(IMON)
100   CONTINUE
C
C Correct annual inputs
C
      IF(TOTPI.EQ.0)THEN
       TOTNEED=TOTPI+(TOTCMEAS-SOIL(6))/(TOTCSIM-SOIL(6))
	ELSE
        TOTNEED=TOTPI*(TOTCMEAS-SOIL(6))/(TOTCSIM-SOIL(6))
	ENDIF
C
C Apportion corrected annual inputs by the months
C
      IF(TOTPI.GT.0)THEN 
        DO 200 IMON=1,12
	    THISPI(IMON)=TOTNEED*THISPI(IMON)/TOTPI
200     CONTINUE
      ELSE
	  DO 300 IMON=1,12
	    THISPI(IMON)=0
300     CONTINUE
      ENDIF
      END
C
C------------------------------------------------------------
C
      SUBROUTINE HILLIER_SOLVER(
C Input 
     &                   ALPHA_IL,AWC_IL,BETA_IL,CLAYPER,
     &                   DRRAT,ITFUNC,IMFUNC,KDPM,KRPM,KBIO,KHUM,
     &                   MEASTOC_IL,PH_IL,PHP1_IL,PHP2_IL,
     &                   TEMP_IL,WMAX_IL,WSAT_IL,
C Input / Output
     &                   PI_IL,ICFACTOR_IL,
C Output
     &                   BORGM,DORGM,HORGM,RORGM,IL)
C
C Initialisation using an analytical solution of the RothC equilibrium run
C worked out by Jon Hillier
C
C
C Variables local to this routine
C
      REAL ACTTOC					! Active organic C (kgC/ha/layer)
	REAL ANNPI					! Annual plant inputs (kg C / ha / yr)
     	INTEGER ICROP(12)			! Crop cover (0=no;1-yes)
	INTEGER IMON				! Month counter
	REAL PROPCO2				! Proportion of decomposition emitted as CO2
	REAL RATEMODD				! Average rate modifier (dpm & bio) over 1 year
	REAL RATEMODR				! Average rate modifier (rpm & hum) over 1 year
	REAL SDPM					! Amount leaving DPM in each timestep
	REAL SRPM					! Amount leaving RPM in each timestep
	REAL SBIO					! Amount leaving BIO in each timestep
	REAL SHUM					! Amount leaving HUM in each timestep
	REAL TOTPI					! Total annual plant input before recalculation (kg C / ha / yr)
C
C Variables passed to/from this routine
C
	REAL ALPHA_IL				! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL AWC_IL(12)				! IN:Long term ave.soil moist.def.(mm)
	REAL BETA_IL				! IN/OUT: Prop.HUM produced on BIO decompn
	REAL BORGM					! OUT:C in soil biomass (kgC/ha/layer)
	REAL CLAYPER				! IN:% Clay
	REAL DORGM					! OUT:C in soil DPM (kgC/ha/layer)
      REAL DRRAT					! IN:Ratio of DPM:RPM. If set to 
								!    zero, this will be worked out using a call to SET_DPMRPMRATIO
	REAL HORGM					! OUT:C in soil humus (kgC/ha/layer)
  	REAL ICFACTOR_IL			! IN/OUT:Adjustment in rate needed to	
								!    achieve measured NPP and TOC for this layer	
	REAL IOM_IL					! IN:Inert organic C (kgC/ha/layer)
	INTEGER ITFUNC,IL				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IMFUNC				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
	REAL KBIO					! IN: Rate constant for BIO decompn/yr
	REAL KDPM					! IN: Rate constant for DPM decompn/yr
	REAL KHUM					! IN: Rate constant for HUM decompn/yr
	REAL KRPM					! IN: Rate constant for RPM decompn/yr
      REAL MEASTOC_IL				! IN:Total measured organic C (kgC/ha/layer)
	REAL PH_IL					! IN:pH of this soil layer
	REAL PHP1_IL				! IN: pH below which decomp.zero
	REAL PHP2_IL				! IN: pH above which decomp.max.
	REAL PI_IL(12)				! IN/OUT: Equilibrium plant C input each month
								!         in each layer (tC/ha/month/layer)
	REAL RORGM					! OUT:C in soil RPM (kgC/ha/layer)
	REAL TEMP_IL(12)			! IN:Long term ave.soil moist.def.(mm)
	REAL WMAX_IL				! IN:Available water at field cap. (mm/layer)
	REAL WSAT_IL				! IN:Available water at saturation (mm/layer)
      REAL EDPM(12),EDPM1			! Rate mods for spinup
      REAL ERPM(12),ERPM1			!
      REAL EB(12),EB1				!
      REAL EH(12),EH1				!
C
C Set ICFACOR to 1 (no effect)
C
	icfactor_il=1
      NITER=0
C
C Get active OM = TOC - IOM
C
	CALL GET_IOM_FROM_FALLOON_EQN(MEASTOC_IL,IOM_IL)
	ACTTOC=MEASTOC_IL-IOM_IL      
C
C Work out the proportion decomposition emitted as CO2
C
      PROPCO2=1.-(ALPHA_IL+BETA_IL)
C
C Get crop cover from plant input
C Assumption: If there is any plant input, crop cover =1; otherwise assume no crop cover
C
      RATEMODD=0
      RATEMODR=0
	TOTPI=0
      DO 100 M=1,12
		IMON=M
	  TOTPI=TOTPI+PI_IL(IMON)
        IF(PI_IL(IMON).GT.0)THEN
	    ICROP(IMON)=1
	  ELSEIF(PI_IL(IMON).LE.0)THEN
	    ICROP(IMON)=0
        ENDIF
C
C Calculate rate modifiers (Source Bradbury et al, 1993)
C
        CALL MODFACTS_MINER(AWC_IL(IMON),WMAX_IL,WSAT_IL,
     &                WRATED,WRATER,
     &                TEMP_IL(IMON),TRATE,								
     &			    PH_IL,PHP1_IL,PHP2_IL,PHRATE,
     &                ICROP(IMON),CRRATE,ITFUNC,IMFUNC)									
	  
	  EDPM(M)=-KDPM*WRATED*TRATE*CRRATE*PHRATE/12
	  ERPM(M)=-KRPM*WRATER*TRATE*CRRATE*PHRATE/12
      EB(M)=-KBIO*WRATED*TRATE*CRRATE*PHRATE/12
      EH(M)=-KHUM*WRATER*TRATE*CRRATE*PHRATE/12
	  
	  RATEMODD=RATEMODD+(WRATED*TRATE*PHRATE*CRRATE)									
	  RATEMODR=RATEMODR+(WRATER*TRATE*PHRATE*CRRATE)
100   CONTINUE
      RATEMODD=RATEMODD/12
	RATEMODR=RATEMODR/12
C
C Work out the amount leaving each pool during the annual cycle
C
      SDPM=EXP(-1*(KDPM/12)*RATEMODD)*(DRRAT/(1+DRRAT))*TOTPI
     &/(1-EXP(-1*(KDPM/12)*RATEMODD))
      SRPM=EXP(-1*(KRPM/12)*RATEMODR)*(1./(1+DRRAT))*TOTPI
     &/(1-EXP(-1*(KRPM/12)*RATEMODR))
      
      SBIO=TOTPI*ALPHA_IL/(PROPCO2*(1-EXP(-1*(KBIO/12)*RATEMODD)))
      
      SHUM=TOTPI*BETA_IL/(PROPCO2*(1-EXP(-1*(KHUM/12)*RATEMODR)))
		ACTTOCDASH=(SDPM+SRPM+SBIO+SHUM)
      SDPM=SDPM*(20+ACTTOC)/ACTTOCDASH
      SRPM=SRPM*(200+ACTTOC)/ACTTOCDASH
      SBIO=SBIO*(20+ACTTOC)/ACTTOCDASH
      SHUM=SHUM*(200+ACTTOC)/ACTTOCDASH
      ACTTOCDASH=(SDPM+SRPM+SBIO+SHUM)

2022	CONTINUE 

         NITER=NITER+1
	
	IF(ACTTOC.GE.0.0001) THEN 


      PI_IL=PI_IL*(1+(ACTTOCDASH+ACTTOC)/(2*ACTTOCDASH))/2.0
  
      goto 222
	else
	goto 221
	endif
222   continue
           
       DO 101 M=1,12
	   SDPM1=SDPM
       SRPM1=SRPM
       SBIO1=SBIO
       SHUM1=SHUM

           									
      EDPM1=EXP(EDPM(M))
	  ERPM1=EXP(ERPM(M))
      EB1=EXP(EB(M))
      EH1=EXP(EH(M))
        
C
C Adjust C(=BCARB()) in BIO pool due to biological activity
C
        SBIO=(SBIO1*EB1)+
     &         ALPHA_IL*(SBIO1*(1-EB1)+SHUM1*(1-EH1)+         
     &         (PI_IL(M)*(DRRAT/(1+DRRAT))+SDPM1)*(1-EDPM1)+
     &         (PI_IL(M)*(1./(1+DRRAT))+SRPM1)*(1-ERPM1))
C
C Adjust C(=HCARB()) in HUM pool due to biological activity
C
        SHUM=(SHUM1*EH1)+
     &         (BETA_IL*SHUM1*(1-EH1))+
     &         (BETA_IL*SBIO1*(1-EB1))+
     &         (BETA_IL*(PI_IL(M)*(DRRAT/(1+DRRAT))+SDPM1)*(1-EDPM1))+
     &         (BETA_IL*(PI_IL(M)*(1./(1+DRRAT))+SRPM1)*(1-ERPM1))
        
        SDPM=(PI_IL(M)*(DRRAT/(1+DRRAT))+SDPM1)*EDPM1
        SRPM=(PI_IL(M)*(1./(1+DRRAT))+SRPM1)*ERPM1
        
101   continue    
      ACTTOCDASH=(SDPM+SRPM+SBIO+SHUM)

      IF(NITER.GT.10000)then !.OR.ACTTOCDASH.EQ.ACTTOC)THEN
	GOTO 221
      ELSE
      GOTO 2022
	ENDIF


C
C Recalculate plant inputs
C
  221 continue

C
C Get pools at equilibrium
C
          Niter=0
300    continue
          Niter=niter+1
          DO 105 M=1,12
	 SDPM1=SDPM
       SRPM1=SRPM
       SBIO1=SBIO
       SHUM1=SHUM
      									
        EDPM1=EXP(EDPM(M))
	  ERPM1=EXP(ERPM(M))
        EB1=EXP(EB(M))
        EH1=EXP(EH(M))
      
C
C Adjust C(=BCARB()) in BIO pool due to biological activity

C
        SBIO=(SBIO1*EB1)+
     &         ALPHA_IL*(SBIO1*(1-EB1)+
     &         SHUM1*(1-EH1)+
     &         (PI_IL(M)*(DRRAT/(1+DRRAT))+SDPM1)*(1-EDPM1)+
     &         (PI_IL(M)*(1./(1+DRRAT))+SRPM1)*(1-ERPM1))



C
C Adjust C(=HCARB()) in HUM pool due to biological activity

C
        SHUM=(SHUM1*EH1)+
     &         (BETA_IL*SHUM1*(1-EH1))+
     &         (BETA_IL*SBIO1*(1-EB1))+
     &         (BETA_IL*(PI_IL(M)*(DRRAT/(1+DRRAT))+SDPM1)*(1-EDPM1))+
     &         (BETA_IL*(PI_IL(M)*(1./(1+DRRAT))+SRPM1)*(1-ERPM1))
        
        SDPM=(PI_IL(M)*(DRRAT/(1+DRRAT))+SDPM1)*EDPM1
        SRPM=(PI_IL(M)*(1./(1+DRRAT))+SRPM1)*ERPM1
105       continue 
          
          if(niter.lt.100)goto 300

      DORGM=SDPM
      RORGM=SRPM
      BORGM=SBIO
      HORGM=SHUM
      do 111 hh=1,12
          
      CALL MODFACTS_MINER(AWC_IL(hh),WMAX_IL,WSAT_IL,
     &                WRATED,WRATER,
     &                TEMP_IL(hh),TRATE,								
     &			    PH_IL,PHP1_IL,PHP2_IL,PHRATE,
     &                ICROP(hh),CRRATE,ITFUNC,IMFUNC) 
111   continue      


C
C Leave ED_SOLVER
C
      END	

C-------------------------------------------------------------------

C
      SUBROUTINE ED_SOLVER_1(
C Input 
     &                   ALPHA_IL,AWC_IL,BETA_IL,CLAYPER,
     &                   DRRAT,ITFUNC,IMFUNC,KDPM,KRPM,KBIO,KHUM,
     &                   MEASTOC_IL,PH_IL,PHP1_IL,PHP2_IL,
     &                   TEMP_IL,WMAX_IL,WSAT_IL,
C Input / Output
     &                   PI_IL,ICFACTOR_IL,
C Output
     &                   BORGM,DORGM,HORGM,RORGM,IL,
     &                   wiltpoint_il,satwatcont_il)
C
C Initialisation using an analytical solution of the RothC equilibrium run
C worked out by Jon Hillier
C
C
C Variables local to this routine
C
      REAL ACTTOC					! Active organic C (kgC/ha/layer)
	REAL ANNPI					! Annual plant inputs (kg C / ha / yr)
     	INTEGER ICROP(12)			! Crop cover (0=no;1-yes)
	INTEGER IMON				! Month counter
       INTEGER IL
	REAL PROPCO2				! Proportion of decomposition emitted as CO2
	REAL RATEMODD				! Average rate modifier (dpm & bio) over 1 year
	REAL RATEMODR				! Average rate modifier (rpm & hum) over 1 year
	REAL SDPM,SDPM1					! Amount leaving DPM in each timestep
	REAL SRPM,SRPM1					! Amount leaving RPM in each timestep
	REAL SBIO,SBIO1,ACTTOCDASH		! Amount leaving BIO in each timestep
	REAL SHUM,SHUM1					! Amount leaving HUM in each timestep
	REAL TOTPI					! Total annual plant input before recalculation (kg C / ha / yr)
C
C Variables passed to/from this routine
C
	REAL ALPHA_IL				! IN/OUT: Prop.BIO produced on BIO decompn 
	REAL AWC_IL(12)				! IN:Long term ave.soil moist.def.(mm)
	REAL BETA_IL				! IN/OUT: Prop.HUM produced on BIO decompn
	REAL BORGM					! OUT:C in soil biomass (kgC/ha/layer)
	REAL CLAYPER				! IN:% Clay
	REAL DORGM					! OUT:C in soil DPM (kgC/ha/layer)
      REAL DRRAT					! IN:Ratio of DPM:RPM. If set to 
								!    zero, this will be worked out using a call to SET_DPMRPMRATIO
	REAL HORGM					! OUT:C in soil humus (kgC/ha/layer)
  	REAL ICFACTOR_IL			! IN/OUT:Adjustment in rate needed to	
								!    achieve measured NPP and TOC for this layer	
	REAL IOM_IL					! IN:Inert organic C (kgC/ha/layer)
	INTEGER ITFUNC				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IMFUNC,NITER				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
	REAL KBIO					! IN: Rate constant for BIO decompn/yr
	REAL KDPM					! IN: Rate constant for DPM decompn/yr
	REAL KHUM					! IN: Rate constant for HUM decompn/yr
	REAL KRPM					! IN: Rate constant for RPM decompn/yr
      REAL MEASTOC_IL				! IN:Total measured organic C (kgC/ha/layer)
	REAL PH_IL					! IN:pH of this soil layer
	REAL PHP1_IL				! IN: pH below which decomp.zero
	REAL PHP2_IL				! IN: pH above which decomp.max.
	REAL PI_IL(12)				! IN/OUT: Equilibrium plant C input each month
								!         in each layer (tC/ha/month/layer)
	REAL RORGM					! OUT:C in soil RPM (kgC/ha/layer)
	REAL TEMP_IL(12)			! IN:Long term ave.soil moist.def.(mm)
	REAL WMAX_IL				! IN:Available water at field cap. (mm/layer)
	REAL WSAT_IL				! IN:Available water at saturation (mm/layer)
      REAL EDPM(12),EDPM1
      REAL ERPM(12),ERPM1
      REAL EB(12),EB1
      REAL EH(12),EH1,layerdepth
      real wiltpoint_il,watsatcont_il
C
C Set ICFACOR to 1 (no effect)
C
	icfactor_il=1
      NITER=0
C
C Get active OM = TOC - IOM
C
	CALL GET_IOM_FROM_FALLOON_EQN(MEASTOC_IL,IOM_IL)
	ACTTOC=MEASTOC_IL-IOM_IL      
C
C Work out the proportion decomposition emitted as CO2
C
      PROPCO2=1.-(ALPHA_IL+BETA_IL)
C
C Get crop cover from plant input
C Assumption: If there is any plant input, crop cover =1; otherwise assume no crop cover
C
      RATEMODD=0
      RATEMODR=0
	TOTPI=0
      DO 100 M=1,12
          IMON=M
	  TOTPI=TOTPI+PI_IL(IMON) 
        IF(PI_IL(IMON).GT.0)THEN
	    ICROP(IMON)=1
          CRRATE=0.6
	  ELSEIF(PI_IL(IMON).LE.0)THEN
	    ICROP(IMON)=0
          Crrate=1
        ENDIF
        layerdepth=5
C
C Calculate rate modifiers (Source Bradbury et al, 1993)
C
       call modfacts_miner2(layerdepth, AWC_IL(IMON), WMAX_IL,
     &   WSAT_IL, TEMP_IL(IMON), wiltpoint_il, satwatcont_il,
     &   PH_IL,PHP1_IL,PHP2_IL,ITFUNC,IMFUNC, 
     &    WRATED, WRATER,TRATE, PHRATE)
       
    	 EDPM(M)=-KDPM*WRATED*TRATE*CRRATE*PHRATE/12
    	 ERPM(M)=-KRPM*WRATER*TRATE*CRRATE*PHRATE/12
         EB(M)=-KBIO*WRATED*TRATE*CRRATE*PHRATE/12
         EH(M)=-KHUM*WRATER*TRATE*CRRATE*PHRATE/12
     									
	  RATEMODD=RATEMODD+(WRATED*TRATE*PHRATE*CRRATE*ICFACTOR_IL)									
	  RATEMODR=RATEMODR+(WRATER*TRATE*PHRATE*CRRATE*ICFACTOR_IL)


100   CONTINUE
      RATEMODD=RATEMODD/12
	RATEMODR=RATEMODR/12
C
C Work out the amount leaving each pool during the annual cycle
C	
      SDPM=EXP(-1*(KDPM/12)*RATEMODD)*(DRRAT/(1+DRRAT))*TOTPI
     &/(1-EXP(-1*(KDPM/12)*RATEMODD))
      SRPM=EXP(-1*(KRPM/12)*RATEMODR)*(1./(1+DRRAT))*TOTPI
     &/(1-EXP(-1*(KRPM/12)*RATEMODR))
      
      SBIO=TOTPI*ALPHA_IL/(PROPCO2*(1-EXP(-1*(KBIO/12)*RATEMODD)))
      
      SHUM=TOTPI*BETA_IL/(PROPCO2*(1-EXP(-1*(KHUM/12)*RATEMODR)))
      
      ACTTOCDASH=(SDPM+SRPM+SBIO+SHUM)
      SDPM=SDPM*(20+ACTTOC)/ACTTOCDASH
      SRPM=SRPM*(200+ACTTOC)/ACTTOCDASH
      SBIO=SBIO*(20+ACTTOC)/ACTTOCDASH
      SHUM=SHUM*(200+ACTTOC)/ACTTOCDASH
      
      !print*, SDPM,SRPM,SBIO,SHUM

      ACTTOCDASH=(SDPM+SRPM+SBIO+SHUM)

2022	CONTINUE 

         NITER=NITER+1
	
	IF(ACTTOC.GE.0.0001) THEN 
 


      ICFACTOR_IL=ICFACTOR_IL*(1+2*ACTTOCDASH/(ACTTOCDASH+ACTTOC))/2.0
      goto 222
	else
	goto 221
	endif
222   continue
           
       DO 101 M=1,12
	 SDPM1=SDPM
       SRPM1=SRPM
       SBIO1=SBIO
       SHUM1=SHUM
          									
        EDPM1=EXP(EDPM(M)*ICFACTOR_IL)
	  ERPM1=EXP(ERPM(M)*ICFACTOR_IL)
        EB1=EXP(EB(M)*ICFACTOR_IL)
        EH1=EXP(EH(M)*ICFACTOR_IL)

        
C
C Adjust C(=BCARB()) in BIO pool due to biological activity
C
        SBIO=(SBIO1*EB1)+
     &         ALPHA_IL*(SBIO1*(1-EB1)+SHUM1*(1-EH1)+         
     &         (PI_IL(M)*(DRRAT/(1+DRRAT))+SDPM1)*(1-EDPM1)+
     &         (PI_IL(M)*(1./(1+DRRAT))+SRPM1)*(1-ERPM1))
C
C Adjust C(=HCARB()) in HUM pool due to biological activity
C
        SHUM=(SHUM1*EH1)+
     &         (BETA_IL*SHUM1*(1-EH1))+
     &         (BETA_IL*SBIO1*(1-EB1))+
     &         (BETA_IL*(PI_IL(M)*(DRRAT/(1+DRRAT))+SDPM1)*(1-EDPM1))+
     &         (BETA_IL*(PI_IL(M)*(1./(1+DRRAT))+SRPM1)*(1-ERPM1))
        
        SDPM=(PI_IL(M)*(DRRAT/(1+DRRAT))+SDPM1)*EDPM1
        SRPM=(PI_IL(M)*(1./(1+DRRAT))+SRPM1)*ERPM1
   !     write(1216,*)SDPM1,SRPM1,SBIO1,SHUM1
        
101   continue    
      ACTTOCDASH=(SDPM+SRPM+SBIO+SHUM)
      IF(NITER.GT.5000)then !.OR.ACTTOCDASH.EQ.ACTTOC)THEN
	GOTO 221
      ELSE
      GOTO 2022
	ENDIF
C
C Recalculate plant inputs
C
  221 continue
C
C Get pools at equilibrium
C
   !       print*,icfactor_il
          Niter=0
300    continue
          Niter=niter+1
          DO 105 M=1,12
	   SDPM1=SDPM
       SRPM1=SRPM
       SBIO1=SBIO
       SHUM1=SHUM
      									
        EDPM1=EXP(EDPM(M)*ICFACTOR_IL)
        ERPM1=EXP(ERPM(M)*ICFACTOR_IL)
        EB1=EXP(EB(M)*ICFACTOR_IL)
        EH1=EXP(EH(M)*ICFACTOR_IL)
        
C
C Adjust C(=BCARB()) in BIO pool due to biological activity
C
        SBIO=(SBIO1*EB1)+
     &         ALPHA_IL*(SBIO1*(1-EB1)+SHUM1*(1-EH1)+
     &         (PI_IL(M)*(DRRAT/(1+DRRAT))+SDPM1)*(1-EDPM1)+
     &         (PI_IL(M)*(1./(1+DRRAT))+SRPM1)*(1-ERPM1))
C
C Adjust C(=HCARB()) in HUM pool due to biological activity
C
        SHUM=(SHUM1*EH1)+
     &         (BETA_IL*SHUM1*(1-EH1))+
     &         (BETA_IL*SBIO1*(1-EB1))+
     &         (BETA_IL*(PI_IL(M)*(DRRAT/(1+DRRAT))+SDPM1)*(1-EDPM1))+
     &         (BETA_IL*(PI_IL(M)*(1./(1+DRRAT))+SRPM1)*(1-ERPM1))
        
        SDPM=(PI_IL(M)*(DRRAT/(1+DRRAT))+SDPM1)*EDPM1
        SRPM=(PI_IL(M)*(1./(1+DRRAT))+SRPM1)*ERPM1

105       continue 
          
          if(niter.lt.1000)goto 300
      DORGM=SDPM
      RORGM=SRPM
      BORGM=SBIO
      HORGM=SHUM
 
C
C Leave ED_SOLVER
C
      END	


C-----------------------------------------------------------
      SUBROUTINE MONAA(SOIL,CMONTH,MON)
C
C  To calculate amount in carbon compartments for this month
C
      IMPLICIT NONE
C
C Variables local to this routine
C
      REAL SOILB(0:6)
	REAL CMONTH(12,0:5,0:5)
      INTEGER MON
	INTEGER I,J
C
C Variables passed to / from this routine
C
      REAL SOIL(0:6)
C
C Save soil pools
C
      SOILB(0)=SOIL(0)
      DO 10 I=1,5
        SOILB(I)=SOIL(I)
   10 CONTINUE
C
C Adjust pools each month
C
      DO 30 I=1,5
        SOIL(I)=CMONTH(MON,I,0)*SOILB(0)
        DO 20 J=1,5
          SOIL(I)=SOIL(I)+CMONTH(MON,I,J)*SOILB(J)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE MULT1(A,B,C)
C
C Matrix Multiplication    A = B * C
C
      IMPLICIT NONE
C
C Variables local to this routine
C
	INTEGER I,J,K
	REAL ZERO2
C
C Variables passed to/from this routine
C
      REAL A(0:5,0:5),B(0:5,0:5),C(0:5,0:5)
C
C Initialise variables
C
      ZERO2=1E-10
C
C Multiply B and C
C
      DO 30 I=0,5
        DO 20 J=0,5
          A(I,J)=0.0
          DO 10 K=0,5
            IF(B(I,K).LT.ZERO2)B(I,K)=0.0
            IF(C(K,J).LT.ZERO2)C(K,J)=0.0
            A(I,J)=A(I,J)+(B(I,K)*C(K,J))
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
C
      RETURN
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE MULT3(A,B)
C
C Matrix Multiplication    A = B * A
C
C
      IMPLICIT NONE
C
C Variables local to this routine
C
      REAL C(0:5,0:5)
	INTEGER I,J,K
C
C Variables passed to/from this routine
C
      REAL A(0:5,0:5),B(0:5,0:5)
C
      DO 30 I=0,5
        DO 20 J=0,5
          C(I,J)=0.0
          DO 10 K=0,5
            C(I,J)=C(I,J)+(B(I,K)*A(K,J))
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
C
      DO 50 I=0,5
        DO 40 J=0,5
          A(I,J)=C(I,J)
   40   CONTINUE
   50 CONTINUE
C
      RETURN
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE PCAR(TOTCSIM,SOIL)
C
C  Totals carbon content
C
C
      IMPLICIT NONE
C
C Variables local to this routine
C
      INTEGER I
C
C Variables passed to/from this routine
C
      REAL TOTCSIM,SOIL(0:6)

      TOTCSIM=0.0
C
      DO 10 I=1,6
        TOTCSIM=TOTCSIM+SOIL(I)
   10 CONTINUE
C
      RETURN
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE POW10(A)
C
C Matrix operation    A = A ** 10
C
      IMPLICIT NONE
C
C Variables local to this routine
C
	REAL B(0:5,0:5),C(0:5,0:5),D(0:5,0:5)
C
C Variables passed to / from this routine
C
      REAL A(0:5,0:5)
C
C .. B = A ** 2
C
      CALL MULT1(B,A,A)
C
C .. C = A ** 4
C
      CALL MULT1(C,B,B)
C
C .. D = A ** 8
C
      CALL MULT1(D,C,C)
C
C .. A = A ** 10
C
      CALL MULT1(A,D,B)
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE RATEF(AVERAIN,AVEPET,AVETEMP,CLAYPER,DEPTH,
     &                 ICROP,RATEM,PH,PHP1,PHP2,ICFACT)
C
C  To produce monthly rate modifying factors from environmental
C  and crop data.
C
      IMPLICIT NONE
C
C Integers local to this routine
C
      INTEGER MAXLU,MAXSOMLAY
	PARAMETER (MAXLU=6,MAXSOMLAY=10)
      REAL CROPRF,TM
	REAL RMTEMP(12),RMSMD(12),RMCROP(12)
	REAL X1,X2,Y1,Y2,D1,D2,DF,EE
	INTEGER M,MM
	REAL PHRATE
C
C Integers passed to/from this routine
C
      REAL AVERAIN(12),AVEPET(12),AVETEMP(12)
	REAL CLAYPER,DEPTH
	INTEGER ICROP(12)
	REAL RATEM(12)
	REAL PH,PHP1,PHP2
  	REAL ICFACT					! IN:Adjustment in rate needed to	
								!    achieve measured NPP and TOC	
C
C Crop retainment factor
C
      CROPRF=0.6
C
C Calculate the pH modifying factor, PHRATE
C
      IF(PH.LE.PHP1)THEN
	  PHRATE=0
	ELSEIF(PH.GT.PHP1.AND.PH.LT.PHP2)THEN
	  PHRATE=(PH-PHP1)/(PHP2-PHP1)
	ELSEIF(PH.GE.PHP2)THEN
	  PHRATE=1
	ENDIF
C
C Temperature rate modifier
C
      DO 10 M=1,12
        TM=AVETEMP(M)
        IF(TM.LT.-5.0)THEN
          RMTEMP(M)=0.0
        ELSE
          RMTEMP(M)=47.91/(EXP(106.06/(TM+18.27))+1)
        ENDIF
   10 CONTINUE
C
C Find month where rainfall exceeds evaporation (starting from January)
C
      DO 30 M=1,12
        IF(AVERAIN(M).GT.AVEPET(M))GOTO 40
   30 CONTINUE
   40 CONTINUE
C
C Find first month where rainfall < evaporation after the
C month found above where rainfall > evaporation.
C
      DO 50 MM=1,12
        M=M+1
        IF(M.GT.12)M=1
        IF(AVERAIN(M).LT.AVEPET(M))GOTO 70
   50 CONTINUE
C
C If rainfall always exceeds evaporation so no moisture deficit
C occurs and rate modifier always = 1
C
      DO 60 MM=1,12
        RMSMD(MM)=1.0
   60 CONTINUE
      GOTO 130
C
   70 CONTINUE
C
C Rainfall is less than evaporation in month 'M'
C start calculating soil moisture DF its from month before 'M'
C
      M=M-1
      IF(M.LT.1)M=12
      X2=-(20+1.3*CLAYPER-0.01*(CLAYPER*CLAYPER))
      X2=X2*DEPTH/23
      Y2=0.2
      D1=0.556*X2
      D2=X2
      X1=0.444*X2
      Y1=1.0
      DF=0.0
C
C Cumulate moisture deficit (DF) in months following rainfall > PET
C
      DO 120 MM=1,12
        EE=AVEPET(M)
C
C ...If rain DOES EXCEED evapotranspiration, add moisture to the soil profile
C 
        IF(AVERAIN(M).GT.EE)THEN
          DF=DF+(AVERAIN(M)-EE)
          GOTO 90
        END IF
C
C ...If rain DOES NOT EXCEED evapotranspiration, crop IS NOT growing 
C    and moisture deficit is NOT below a critical minimum (D1) then
C    take moisture out of soil profile
C
        IF(ICROP(M).EQ.1)GOTO 80
        IF(DF.LT.D1)GOTO 100
        DF=DF+(AVERAIN(M)-EE)
        IF(DF.LT.D1)DF=D1
        GOTO 100
C
C ...If rain DOES NOT exceed evapotranspiration, crop IS growing 
C    take as moisture out of layer down to a critical minimum (D2)
C
   80   CONTINUE
        DF=DF+(AVERAIN(M)-EE)
        IF(DF.LT.D2)DF=D2
C
C If moisture deficit is positive, set to 0
C
   90   CONTINUE
        IF(DF.GT.0.0)DF=0.0
  100   CONTINUE
C
C Calculate moisture rate modifier from SMD
C...if SMD >= critical value X1, set to maximum rate (1)
C
        IF(DF.GE.X1)THEN
           RMSMD(M)=1.0
           GOTO 110
        END IF
C
C...if SMD >= critical value X2, calculate from equation
C
        IF(DF.GE.X2)THEN
          RMSMD(M)=Y2+(Y1-Y2)*((X2-DF)/(X2-X1))
        GOTO 110
        END IF
        STOP
C
  110   CONTINUE
        M=M+1
        IF(M.GT.12)M=1
  120 CONTINUE
  130 CONTINUE
C
C Obtain rate modifier from crop retaining factor
C
      DO 140 M=1,12
        IF(ICROP(M).EQ.1)THEN
          RMCROP(M)=CROPRF
        ELSE
          RMCROP(M)=1.0
        END IF
  140 CONTINUE
C
C Combine rate modifiers
C
      DO 150 M=1,12
        RATEM(M)=RMTEMP(M)*RMSMD(M)*RMCROP(M)*PHRATE*ICFACT
  150 CONTINUE
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE RATEF_LTA(LTA_AWC,LTA_TEMP,ITFUNC,IMFUNC,WMAX,WSAT,
     &                     DEPTH,RATEM,PH,PHP1,PHP2,ICFACT,
     &                     ICROP)
C
C  To produce monthly rate modifying factors from environmental
C  and crop data.
C
      IMPLICIT NONE
C
C Integers local to this routine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1			! No.of layers in the soil profile
	INTEGER MAXDEPTH			! Maximum depth of the soil profile
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
	INTEGER MAXLU				! Maximum number of land use types
	PARAMETER (MAXLU=6)	
	INTEGER MAXSOMLAY			! Maximum number of SOM layers
	PARAMETER (MAXSOMLAY=3)
	REAL TRATE					! Temperature rate modifier
	REAL WRATED					! Moisture rate modifier
	REAL WRATER					! Moisture rate modifier
	REAL CRRATE					! Crop rate modifier
	REAL PHRATE					! pH rate modifier
	INTEGER M					! Counter 1 for months
	INTEGER IL					! Layer counter
C
C Variables passed to/from this routine
C...Weather factors
C
	REAL LTA_AWC(12,MAXLAYER)	! IN:Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)	! IN:Average air temp this month (deg.C)
C
C ...Model descriptors
C
	INTEGER ITFUNC				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IMFUNC				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
C
C...Soil factors
C
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
	REAL DEPTH					! IN:Depth (cm)
	REAL RATEM(12)				! OUT:Rate modifier per month
	REAL PH						! IN:pH of this soil layer
	REAL PHP1					! IN: pH below which decomp.zero
	REAL PHP2					! IN: pH above which decomp.max.
  	REAL ICFACT					! IN:Adjustment in rate needed to	
								!    achieve measured NPP and TOC	
C
C...Plant factors
C
	INTEGER ICROP(12)			! Crop cover (0=no;1-yes)
C
C Calculate current SUNDIAL layer
C
      IL=DEPTH*MAXLAYER1/MAXDEPTH
C
C For each month...
C
      DO 10 M=1,12
C
C Calculate rate modifiers (Source Bradbury et al, 1993)
C
        CALL MODFACTS_MINER(LTA_AWC(M,IL),WMAX(IL),WSAT(IL),
     &                WRATED,WRATER,
     &                LTA_TEMP(M,IL),TRATE,								
     &			    PH,PHP1,PHP2,PHRATE,
     &                ICROP(M),CRRATE,ITFUNC,IMFUNC)									
C
C Overall rate modifier
C      
        RATEM(M)=((WRATED+WRATER)/2)*TRATE*PHRATE*CRRATE*ICFACT
10    CONTINUE
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE RUNC14(ISTARTYR,ISTOPYR,
     &                  RAINFILE,PETFILE,TEMPFILE,SOIL,
     &                  CLAYPER,DEPTH,ICROP,PH,PHP1,PHP2,
     &                  MCSTART,DEC,FDEC2,FPDEC2,FPLANT,THISPI,FFYM,
     &                  FYMADD,ICFACT)
C
C  Running Carbon model with varying radiocarbon activity
C
      IMPLICIT NONE
C
C Variables local to this routine
C
	INTEGER NYEARS,K,IYEAR,M,MAXLU,MAXSOMLAY
	PARAMETER (MAXLU=6,MAXSOMLAY=10)
C
C Variables passed to / from this routine
C
	INTEGER MON					! Month
	INTEGER ISTARTYR			! Starting year
	INTEGER ISTOPYR				! End year
	CHARACTER*20 RAINFILE		! IN/OUT:File of average rainfall (mm/month)
	CHARACTER*20 TEMPFILE		! IN/OUT:File of average monthly temp. (deg.C)
	CHARACTER*20 PETFILE		! IN/OUT:File of average PET (mm/month)
	REAL CLAYPER				! % Clay
	REAL DEPTH					! Depth
	INTEGER ICROP(12)			! Crop cover (0=no;1-yes)
	REAL PH						! Current pH of the soil
	REAL PHP1					! pH below which rate is zero
	REAL PHP2					! pH above which rate is optimum
      INTEGER MCSTART
	REAL DEC(5)
	REAL FDEC2(5)
	REAL FPDEC2(5)
	REAL FPLANT(5)				! Proportion plant material going to each pool
	REAL THISPI(12)				! IN/OUT: Equilibrium plant C input each month
								!         in each layer (tC/ha/month/layer)
	REAL FFYM(5)				! Proportion fym going to each pool
	REAL FYMADD(12)				! FYM input (t/ha/month)
	REAL CYEAREQ(0:5,0:5)
      REAL SOIL(0:6)				! Soil C 1=DPM;2=RPM;3=BIOa;4=BIOz;5=HUM;6=IOM 
								! (t/ha)
	REAL RAIN(12)				! Monthly rainfall (mm)
	REAL PET(12)				! Monthly potential evapotranspiration (mm)
	REAL TEMP(12)				! Monthly average air temperature (deg.C)
	REAL RATEM(12)				! Rate modifier per month
	REAL CMONTH(12,0:5,0:5)
  	REAL ICFACT					! IN:Adjustment in rate needed to	
								!    achieve measured NPP and TOC	
C
C Set start of simulation to January
C
      MON=1
	IYEAR=0
C
C Open weather data files
C
      OPEN(12,FILE=RAINFILE,STATUS='OLD',ERR=111)
      OPEN(13,FILE=PETFILE,STATUS='OLD',ERR=111)
      OPEN(14,FILE=TEMPFILE,STATUS='OLD',ERR=111)
	GOTO 101
111   CONTINUE
      PRINT*,'Cannot open weather files. Press any key to continue.'
	READ(*,*)
	STOP
101   CONTINUE
C
C Do calculation for NYEARS years
C
      NYEARS=ISTOPYR
	OPEN(55,FILE='TEMPC.OUT',STATUS='UNKNOWN',ERR=111)
      DO 20 K=ISTARTYR,ISTOPYR
        IYEAR=IYEAR+1
C
C Set this years weather data
C        
        CALL GETMET(RAIN,PET,TEMP)														
C
C calculate (1) monthly rate modifying factors (2) set transition matrices
C
	  CALL RATEF(RAIN,PET,TEMP,CLAYPER,DEPTH,
     &                 ICROP,RATEM,PH,PHP1,PHP2,ICFACT)														
	  CALL SETMAT(MCSTART,RATEM,DEC,FDEC2,FPDEC2,CYEAREQ,
     &                 FPLANT,THISPI,FFYM,FYMADD,CMONTH)
C
C Do calculation for 1 year
C
        DO 10 M=1,12
C
C Calculate amount in soil carbon compartments for this month
C
          CALL MONAA(SOIL,CMONTH,MON)
C
C Move on to next month
C
          MON=MON+1
          IF(MON.GT.12)THEN
            MON=1
          END IF
   10   CONTINUE
	WRITE(55,30)IYEAR,SOIL(1),SOIL(2),SOIL(3),SOIL(4),SOIL(5),SOIL(6)
   30 FORMAT(I4,6(F10.0,2X))
   20 CONTINUE
      CLOSE(12)
      CLOSE(13)
	CLOSE(14)
	CLOSE(55)
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE RUNCEQ(NCYEARS,SOIL,SOILB,CYEAREQ,TOTCSIM)
C
C
C  Running Equilibrium ( steady state ) Carbon model only
C
C
      IMPLICIT NONE
C
C Variables local to this routine
C
      INTEGER K,K1,K2,KS,ICYEAR
C
C Variables passed to/from this routine
C
	INTEGER NCYEARS
      REAL SOIL(0:6),SOILB(0:6),CYEAREQ(0:5,0:5)
	REAL TOTCSIM
C
C Output simple run results at year 50 and 100
C
      K1=1
      K2=10
      KS=1
	ICYEAR=0
C
   10 CONTINUE
C
      IF(NCYEARS.LT.K2)K2=NCYEARS
C
      DO 20 K=K1,K2,KS
C
C Set weather year (for access to correct weather data)
C
        ICYEAR=ICYEAR+KS
C
C Transform soil carbon
C
        CALL YAAEQ(SOIL,SOILB,CYEAREQ)
C
C Sum result
C
	  IF(K.EQ.NCYEARS)THEN
	    CALL PCAR(TOTCSIM,SOIL)
        END IF
   20 CONTINUE
C   
C Multiply number of years by power 10
C
      IF(K2.EQ.(KS*10))THEN												
        K1=K2*2															
        K2=K2*10														
        KS=KS*10														
        CALL POW10(CYEAREQ)												
        GOTO 10
      END IF
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE RUNTOEQ(TOTCSIM,SOIL,MCSTART,CLAYPER,DEPTH,
     &                   ICROP,RATEM,DEC,DECOMP,
     &                   FDEC1,FPDEC1,FDEC2,FPDEC2,
     &                   FPLANT,THISPI,FFYM,FYMADD,NCYEARS,
     &                   PH,PHP1,PHP2,ICFACT,LTA_AWC,LTA_TEMP,
     &                   ITFUNC,IMFUNC,WMAX,WSAT)
C
C Steady-state calculation
C
      IMPLICIT NONE
C
C Variables local to this routine
C
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLU				! Maximum number of land use types
	PARAMETER (MAXLU=6)	
	INTEGER MAXSOMLAY			! Maximum number of SOM layers
	PARAMETER (MAXSOMLAY=10)
C
C Variables passed to/from this routine
C ...Weather factors
C
	REAL LTA_AWC(12,MAXLAYER)	! IN:Long term available water (mm)
	REAL LTA_TEMP(12,MAXLAYER)		! IN:Average air temp this month (deg.C)
C
C ...Model descriptors
C
	INTEGER ITFUNC				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER IMFUNC				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
C
C...Soil factors
C
	REAL WMAX(MAXLAYER)			! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)			! IN:Available water at saturation (mm/layer)
	REAL TOTCSIM				! Simulated C input (t/ha)
      REAL SOIL(0:6)				! Soil C 1=DPM;2=RPM;3=BIOa;4=BIOz;5=HUM;
	                            ! 6=IOM (t/ha)
      REAL SOILB(0:6)				! Soil C at beginning 1=DPM;2=RPM;3=BIOa;
	                            ! 4=BIOz;5=HUM;6=IOM (t/ha)
	REAL DEC(5)
	REAL DECOMP(5)
	INTEGER MCEND
	INTEGER MCSTART
	REAL CLAYPER				! % Clay
	REAL DEPTH					! Depth (cm)
	INTEGER NCYEARS
	REAL RATEM(12)				! Rate modifier per month
	REAL CYEAREQ(0:5,0:5)
	REAL FDEC1(5)
	REAL FPDEC1(5)
	REAL FDEC2(5)
	REAL FPDEC2(5)
	REAL PH						! IN:pH of soil in this layer
	REAL PHP1					! IN: pH below which decomp.zero
	REAL PHP2					! IN: pH above which decomp.max.
	REAL CMONTH(12,0:5,0:5)
  	REAL ICFACT					! IN:Adjustment in rate needed to	
								!    achieve measured NPP and TOC	
C
C ...Plant factors
C
	INTEGER ICROP(12)			! Crop cover (0=no;1-yes)
	REAL FPLANT(5)				! Prop.plant material going to each pool
								! 1=arable;2=grass;3=forestry;4=semi-natural
	REAL THISPI(12)				! Plant input (tC/ha/month)
C
C ...Manure factors
C
	REAL FFYM(5)				! IN:Proportion fym going to each pool
	REAL FYMADD(12)				! IN:FYM input (t/ha/month)
C
c	REAL AVERAIN(12),AVETEMP(12),AVEPET(12)								! ROTHC rate modifiers
C
C Primary calculations of radiocarbon activity, monthly decomposition
C and last month of year
C
	CALL CALC1(SOIL,SOILB,DEC,DECOMP,MCEND,MCSTART)
C
C Calculate CO2:(BIO+POM) ratio from % clay value; 
C
	CALL BPCALC(CLAYPER,FDEC1,FPDEC1,FDEC2,FPDEC2)
c      CALL AVEMET(AVERAIN,AVEPET,AVETEMP)								! ROTHC rate modifiers
c      CALL RATEF(AVERAIN,AVEPET,AVETEMP,CLAYPER,DEPTH,					! ROTHC rate modifiers
c     &                 ICROP,RATEM,PH,PHP1,PHP2,ICFACT)					! ROTHC rate modifiers
      CALL RATEF_LTA(LTA_AWC,LTA_TEMP,ITFUNC,IMFUNC,WMAX,WSAT,			! ECOSSE rate modifiers
     &                     DEPTH,RATEM,PH,PHP1,PHP2,ICFACT,				! ECOSSE rate modifiers
     &                     ICROP)											! ECOSSE rate modifiers
	CALL PCAR(TOTCSIM,SOIL)
C																		
C Set monthly transition matrices										
C																		
 	CALL SETMAT(MCSTART,RATEM,DEC,FDEC2,FPDEC2,CYEAREQ,
     &            FPLANT,THISPI,FFYM,FYMADD,CMONTH)
C
C Run to equilibrium
C
      CALL RUNCEQ(NCYEARS,SOIL,SOILB,CYEAREQ,TOTCSIM)
	END
C
C------------------------------------------------------------
C
      SUBROUTINE RUNTONOW(TOTCSIM,SOIL,MCSTART,CLAYPER,TIMEFROMEQ,
     &                   DEPTH,RAINFILE,PETFILE,TEMPFILE,
     &                   ICROP,DEC,DECOMP,
     &                   FDEC1,FPDEC1,FDEC2,FPDEC2,
     &                   FPLANT,THISPI,FFYM,FYMADD,
     &                   PH,PHP1,PHP2,ICFACT,DRRAT)
C
C Run from equilibrium to present day
C
      IMPLICIT NONE
C
C Variables local to this routine
C
	INTEGER MAXLU				! Maximum number of land use types
	PARAMETER (MAXLU=6)	
	INTEGER MAXSOMLAY			! Maximum number of SOM layers
	PARAMETER (MAXSOMLAY=10)
C
C Variables passed to this routine
C...Timing factors
C
      INTEGER INYEAR
	INTEGER ISTARTYR
	INTEGER ISTOPYR
	INTEGER TIMEFROMEQ
	INTEGER MCSTART
C
C...Soil factors
C
	REAL TOTCSIM				! Simulated C input (t/ha)
      REAL SOIL(0:6)				! Soil C 1=DPM;2=RPM;3=BIOa;4=BIOz;5=HUM;
	                            ! 6=IOM (t/ha)
	REAL CLAYPER				! % Clay
	REAL DEPTH					! Depth
	CHARACTER*20 RAINFILE		! IN/OUT:File of average rainfall (mm/month)
	CHARACTER*20 TEMPFILE		! IN/OUT:File of average monthly temp. (deg.C)
	CHARACTER*20 PETFILE		! IN/OUT:File of average PET (mm/month)
	REAL FDEC1(5)
	REAL FPDEC1(5)
	REAL FDEC2(5)
	REAL FPDEC2(5)
	REAL PH							! IN:pH of soil in this layer
	REAL PHP1						! IN: pH below which decomp.zero
	REAL PHP2						! IN: pH above which decomp.max.
  	REAL ICFACT						! IN:Adjustment in rate needed to	
	REAL DECOMP(5)
      REAL DEC(5)
	REAL TPAR
	REAL TPPAR
	REAL TFPAR
C
C ...Plant factors
C
	INTEGER ICROP(12)				! Crop cover (0=no;1-yes)
	REAL FPLANT(5)					! Prop.plant material going to each pool
									! 1=arable;2=grass;3=forestry;4=semi-natural
	REAL THISPI(12)					! Plant input (tC/ha/month)				
	REAL DRRAT						! DPM:RPM ratio
C
C ...Manure factors
C
	REAL FFYM(5)					! IN:Proportion fym going to each pool
	REAL FYMADD(12)					! IN:FYM input (t/ha/month)
C
C Set up ROTH-C calculation for short-term phase
C
      INYEAR=1
	CALL PCAR(TOTCSIM,SOIL)
	IF(TOTCSIM.GT.0)THEN
	  CALL SETSIM(MCSTART,SOIL,DRRAT,
     &                 FPLANT,FFYM,FDEC1,FPDEC1,DECOMP,
     &                 TPAR,TPPAR,TFPAR)
C
C Calculate (1) CO2:(BIO+POM) ratio from % clay value; 
C (2) run carbon model for 100 years
C
	  CALL BPCALC(CLAYPER,FDEC1,FPDEC1,FDEC2,FPDEC2)
	  ISTARTYR=1
	  ISTOPYR=TIMEFROMEQ
        CALL RUNC14(ISTARTYR,ISTOPYR,
     &                  RAINFILE,PETFILE,TEMPFILE,SOIL,
     &                  CLAYPER,DEPTH,ICROP,PH,PHP1,PHP2,
     &                  MCSTART,DEC,FDEC2,FPDEC2,FPLANT,THISPI,FFYM,
     &                  FYMADD,ICFACT)											
      ENDIF
      TOTCSIM=SOIL(1)+SOIL(2)+SOIL(3)+SOIL(4)+SOIL(5)+SOIL(6)
	END
C
C------------------------------------------------------------
C
      SUBROUTINE SETEQ(ISYEAR,NCYEARS,MCSTART,SOIL,DRRAT,
     &                 FPLANT,FFYM,FDEC1,FPDEC1,DECOMP,
     &                 TPAR,TPPAR,TFPAR)
C
C Sets up simulation for steady sate phase
C
      IMPLICIT NONE
C
C Variables local to this routine
C
      INTEGER I,IPAR
C
C Variables passed to/from this routine
C
 	REAL SOIL(0:6)
      INTEGER ISYEAR,NCYEARS,MCSTART
	REAL FPLANT(5),FFYM(5),FDEC1(5),FPDEC1(5),DECOMP(5)
	REAL DRRAT,TPAR,TPPAR,TFPAR
C
C Set values
C Starting year of simulation 0 (equilibrium run)
      ISYEAR=0
C Run simulation for 10000 years
      NCYEARS=10000
C Start month = January
      MCSTART=1
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
C
C Set all active soil pools to zero ready for equilibrium run
C
      SOIL(0)=1
	DO 100 I=1,5
	  SOIL(I)=0
100   CONTINUE
      END
C
C------------------------------------------------------------
C
      SUBROUTINE SET_DPMRPMRATIO(LUCODE,DRRAT)
C
C Set up previous land use and DPM:RPM ratio
C
      IMPLICIT NONE
C
C Variables passed to/from this routine
C
	INTEGER LUCODE
	REAL DRRAT
C
C Set variables 
C
      IF(LUCODE.EQ.1)THEN
        DRRAT = 1.44
      ELSEIF(LUCODE.EQ.2)THEN
c        DRRAT = 0.41				! after Shirato & Yokozava (2006)
        DRRAT = 1.44				
      ELSEIF(LUCODE.EQ.3)THEN
	   DRRAT = 0.25
       ! DRRAT = 0.19
c        
		!DRRAT =0.67
      ELSEIF(LUCODE.EQ.4)THEN
        DRRAT = 1.44
      ELSEIF(LUCODE.EQ.5)THEN
        DRRAT = 1.44
      ELSEIF(LUCODE.EQ.6)THEN
        DRRAT = 1.0  ! MLR changed for ELUM
      ELSEIF(LUCODE.EQ.7)THEN
        DRRAT = 1.44  ! MLR changed for ELUM
      ELSEIF(LUCODE.EQ.8)THEN
        DRRAT = 1.44  ! MLR changed for ELUM
      ELSEIF(LUCODE.EQ.9)THEN
        DRRAT = 0.20  ! MLR changed for ELUM
      ELSEIF(LUCODE.EQ.10)THEN
        DRRAT = 1.44  ! MLR changed for ELUM
      ENDIF
	END
C
C------------------------------------------------------------
C
      SUBROUTINE SETMAT(MCSTART,RATEM,DEC,FDEC2,FPDEC2,CYEAREQ,
     &                  FPLANT,THISPI,FFYM,FYMADD,CMONTH)
C
C
C  To set up monthly Transition matrices
C  ( for carbon and nitrogen calculations )
C
C  and a yearly Transition matrix
C  ( for carbon calculations only )
C
C
      IMPLICIT NONE
C
C Variables local to this routine
C
      INTEGER MAXLU
	INTEGER MAXSOMLAY
	PARAMETER (MAXLU=6,MAXSOMLAY=10)
      REAL D(5)
	REAL E(5)
	REAL G(5)
	REAL GP(5)
      REAL CT(0:5,0:5)
	INTEGER M
	INTEGER MM
	INTEGER I
	INTEGER J
C
C Variables passed to/from this routine
C
	INTEGER MCSTART				! Number of starting month for ROTHC run
	REAL RATEM(12)
	REAL DEC(5)
	REAL FDEC2(5)
	REAL FPDEC2(5)
	REAL FPLANT(5)				! Proportion plant material going to each pool
	REAL THISPI(12)				! IN/OUT: Equilibrium plant C input each month
								!         in each layer (tC/ha/month/layer)
	REAL FFYM(5)				! Proportion fym going to each pool
	REAL FYMADD(12)				! FYM input (t/ha/month)
	REAL CYEAREQ(0:5,0:5)
      REAL CMONTH(12,0:5,0:5)
C
      M=MCSTART
C
C Clear yearly transition matrices
C
      CALL UNITM(CYEAREQ)
C
      DO 120 MM=1,12
C
C Clear this month's transition matrix
C
        CALL UNITM(CT)
C
        IF(RATEM(M).LT.0.0001)GOTO 50
C
        DO 10 I=1,5
          D(I)=DEC(I)*RATEM(M)
   10   CONTINUE
C
        DO 20 I=1,5
          E(I)=EXP(-D(I))
          G(I)=1-E(I)
          GP(I)=0.0
   20   CONTINUE
C
        GP(5)=G(5)
        G(5)=0.0
C
        DO 40 I=1,5
          DO 30 J=1,5
            CT(I,J)=FDEC2(I)*G(J)+FPDEC2(I)*GP(J)
   30     CONTINUE
          CT(I,I)=CT(I,I)+E(I)
   40   CONTINUE
C
   50   CONTINUE
        DO 60 I=1,5
          CT(I,0)=FPLANT(I)*THISPI(M)
		CT(I,0)=CT(I,0)+FFYM(I)*FYMADD(M)
   60   CONTINUE
C
C Save monthly carbon inputs
C
        DO 100 I=0,5
          DO 90 J=0,5
            CMONTH(M,I,J)=CT(I,J)
   90     CONTINUE
  100   CONTINUE
C
        CALL MULT3(CYEAREQ,CT)
C
        M=M+1
        IF(M.GT.12)M=1
  120 CONTINUE
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE UNITM(A)
C
C Makes matrix A a unit matrix
C
C
      IMPLICIT NONE
C
C Variables local to this routine
C
	INTEGER I,J
C
C Variables passed to/from this routine
C
      REAL A(0:5,0:5)
C
      DO 20 I=0,5
        DO 10 J=0,5
          A(I,J)=0.0
   10   CONTINUE
        A(I,I)=1.0
   20 CONTINUE
C
      RETURN
      END
C
C------------------------------------------------------------
C
      SUBROUTINE YAAEQ(SOIL,SOILB,CYEAREQ)
C
C
C  To calculate yearly (1,10,100 or 1000 yearly) amounts and ages
C  of carbon compartments  ( For equilibrium carbon model only )
C
      IMPLICIT NONE
C
C Variables local to this routine
C
	INTEGER I,J
C
C Variables passed to/from this routine
C
      REAL SOIL(0:6),SOILB(0:6),CYEAREQ(0:5,0:5)
C
      DO 10 I=0,5
        SOILB(I)=SOIL(I)
   10 CONTINUE
C
      DO 30 I=0,5
        SOIL(I)=0.0
        DO 20 J=0,5
          SOIL(I)=SOIL(I)+CYEAREQ(I,J)*SOILB(J)
   20   CONTINUE
   30 CONTINUE
C
      RETURN
      END

