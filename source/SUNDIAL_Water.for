C------------------------------------------------------------
C
C Rothamsted Carbon and Nitrogen Turnover Model
C Soil Water Routines
C
C Adapted from SUNDIAL (MAGEC)
C by Jo Smith & Kevin Coleman
C 02/03/05
C-------------------------------------------------------------
C
C EXTERNAL SUBROUTINES
C 1. RUN_SUNDIAL_WATER_SUBMONTHLY() (author: Mark Richards)
C 2. DRAIN_SUNDIAL_WATER2 (author: Mark Richards)
C 3. DRAIN_SUNDIAL_WATER
C 4. EVAP_SUNDIAL_WATER2 (author: Mark Richards) 
C 3. EVAP_SUNDIAL_WATER
C 4. INIT_SUNDIAL_WATER
C 5. LIMIT_DRAIN
C
C-------------------------------------------------------------------------------
C
      subroutine run_sundial_water_submonthly(month_rain, wmax, soilw,
     &                                        wsat, flowprop,
     &                                        month_drainw, 
     &                                        month_swd_below_fc,
     &                                        month_evap, airtemp,
     &                                        icover)
      ! Runs the water model (drainage and evaporation) for one month using a daily timestep.
      ! When the model runs on a monthly timestep the water model is run on a daily timestep using
      ! this subroutine. It was developed because if the water model is run on a monthly timestep
      ! when monthly rainfall is greater than monthly PET the soil can end up drying out too much.
      ! This is because the month's rainfall is added first, with a proportion of any excess water
      ! above field being lost due to flowrprop. The month's evaporation is then applied resulting
      ! in the upper profile drying out to field capacity. When rain and evaporation are applied
      ! daily the loss of excess water is reduced resulting in a wetter profile. Daily rainfall
      ! and evaporation are obtained by averaging the monthly figures out over 30 days. 
      !
      ! Developed by Mark Richards by merging DRAIN_SUNDIAL_WATER2() and EVAP_SUNDIAL_WATER()
      
      implicit none

      ! Constants
      integer, parameter :: maxlayer = 60  ! No. of layers in the soil profile
      integer, parameter :: maxdepth = 300 ! Depth of soil profile [cm]

      ! Arguments with intent in
      real, intent(in) :: airtemp         ! Mean air temperature [deg C]
      real, intent(in) :: month_evap      ! Potential evapotranspiration for the month [mm/month]
      real, intent(in) :: flowprop        ! Proportion of the excess water that 
                                          ! drains out of the profile each iteration
      real, intent(in) :: month_rain      ! Rainfall for the month [mm/month]
      real, intent(in) :: wmax(maxlayer)  ! Available water at field capacity [mm/layer]
      real, intent(in) :: wsat(maxlayer)  ! Available water at saturation [mm/layer]
      integer, intent(in) :: icover       ! 1 if land is covered with vegetation, otherwise 0

      ! Arguments with intent inout
      real, intent(inout) :: soilw(maxlayer)  ! Available water [mm/layer]

      ! Arguments with intent out
      real, intent(out) :: month_drainw(maxlayer) ! Net water draining through each layer
                                                  ! (used by leaching functions) [mm/layer]
      real, intent(out) :: month_swd_below_fc(maxlayer)  ! Soil water deficit below field capacity [mm/layer]

      ! Local variables
      integer :: d                      ! Loop counter
      real :: day_evap                  ! Daily evaporation [mm/day]
      real :: day_rain                  ! Daily rainfall [mm/day]
      real :: day_drainw(maxlayer)      ! Net water draining through each layer each day
                                        ! (used by leaching functions) [mm/layer/day]
      integer :: i                      ! Local layer counter
                                        ! is above field capacity [mm]  
      real :: day_swd_below_fc(maxlayer) ! Daily soil water deficit below field capacity [mm/layer]

      ! Initialise variables
      day_rain = month_rain / 30.0      ! Calculate daily rain from monthly rain
      day_evap = month_evap / 30.0      ! Calculate daily evap from monthly evap
      month_drainw = 0.0
      month_swd_below_fc = 0.0

      do d=1, 30  ! days
        call drain_sundial_water2(day_rain, wmax, soilw, wsat,
     &                            flowprop, day_drainw,
     &                            day_swd_below_fc)
        
      ! Add the day's drainw and swd_below_fc to the monthly totals
      do i=1, maxlayer
            month_drainw(i) = month_drainw(i) + day_drainw(i)
            month_swd_below_fc(i) = month_swd_below_fc(i) + 
     &                              day_swd_below_fc(i)
        enddo

        call evap_sundial_water2(day_evap, airtemp, icover, wmax, soilw)

      enddo  ! Days

      end subroutine run_sundial_water_submonthly
C
C-------------------------------------------------------------------------------
C
	subroutine drain_sundial_water2(rain, wmax, soilw, wsat, flowprop,
     &								drainw, swd_below_fc)

	! Simulates the percolation of infiltrating water through the soil profile
	! using the "tipping bucket" approach.

	! The subroutine operates as follows:
	! 1. The water above field capacity in each layer is removed and added to
	!    an "excess water" pool (called water_above_fc).
	! 2. Infiltration is added to each layer, starting at the top of the 
	!    profile, until each layer is filled to field capacity. Once a layer 
	!    is filled, the remaining infiltration "tips" into the next layer down 
	!    until it is also filled to field capacity and so on.
	! 3. Any infiltration left-over after (2) is added to the excess water pool
	! 4. The excess pool is partitioned between drainage (water lost from the
	!    profile) and refill according to the value of flowprop. The refill
      !    water is used to fill up layers to saturation from the bottom of the 
	!    profile upwards.
      !
      ! Note: The drainw array records the drainage though each layer for
      ! use by leaching subroutines. All downward movement of water through a 
      ! layer is added to drainw and all upward movement (i.e. from refill)
      ! is subtracted from drainw in order to record the net downward drainage
      ! of water.

	implicit none

	! Constants
	integer, parameter :: maxlayer = 60	 ! No.of layers in the soil profile

	! Arguments with intent in
	real, intent(in) :: flowprop        ! Proportion of the excess water that 
	                                    ! drains out of the profile each iteration
	real, intent(in) :: rain            ! Rainfall [mm/timestep]
	real, intent(in) :: wmax(maxlayer)  ! Available water at field capacity [mm/layer]
	real, intent(in) :: wsat(maxlayer)	! Available water at saturation [mm/layer]

	! Arguments with intent inout
	real, intent(inout) :: soilw(maxlayer)  ! Available water [mm/layer]

	! Arguments with intent out
	real, intent(out) :: drainw(maxlayer) ! Net water draining through each layer
	                                      ! (used by leaching functions) [mm/layer]
	real, intent(out) :: swd_below_fc(maxlayer)	! Soil water deficit below field capacity [mm/layer]

	! Local variables
	integer :: i			! Local layer counter
	real :: drainage		! Water draining out of bottom of profile [mm]
	real :: infiltration   	! Water infiltrating into soil [mm]
	real :: refill			! The excess water available for filling layers up 
	                        ! to saturation [mm]
	real :: swd_below_sat   ! Soil water deficit below saturation of a layer [mm]
	real :: water_above_fc  ! Total amount of (excess) water in profile that
	                        ! is above field capacity [mm]

      ! Initialise variables
      drainw = 0.0
      water_above_fc = 0.0
      infiltration = rain
      
    	! 1. Drain the water above each field capacity to the bottom of the profile
      do i=1, maxlayer
          water_above_fc = water_above_fc + max(soilw(i) - wmax(i), 0.0)
          drainw(i) = water_above_fc  ! Amount of saturated flow through the layer
          soilw(i) = min(soilw(i), wmax(i))
      enddo  ! i (layers)   
    
	! 2. Fill up the layers to field capacity using infiltration water
      do i=1, maxlayer 
	    swd_below_fc(i) = wmax(i) - soilw(i)  ! Calc soil water deficit of layer
	    ! If there is sufficient infiltration left to fill layer to FC...
	    if (infiltration >= swd_below_fc(i)) then	
	        soilw(i) = wmax(i)
	        infiltration = infiltration - swd_below_fc(i)
	    ! ...else there is insufficient infiltration left to fill layer to FC
	    else  
		    soilw(i) = soilw(i) + infiltration
	        infiltration = 0.0
	        exit  ! No infiltration left so quit loop early
	    endif
	    ! Record the amount of infiltration water passing through the layer
          ! (used in leaching functions)
	    drainw(i) = drainw(i) + infiltration  
      enddo  ! i (layers)
      
	! 3. Partition left-over water into refill and drainage
	water_above_fc = water_above_fc + infiltration  
	drainage = water_above_fc * flowprop  ! Water draining out of the profile
	refill = water_above_fc - drainage    ! Water available to fill layers to saturation 
     
	! 4. Use refill water to fill layers up to saturation from the bottom of 
	!    the profile
	if (refill > 0.0) then
          do i=maxlayer, 1, -1
              ! Subtract the refill passing into this layer from the original drainage
              ! in order to get the NET drainage through the layer
              drainw(i) = drainw(i) - refill
              swd_below_sat = wsat(i) - soilw(i)
	        ! If sufficient refill available to fill layer to saturation...
	        if (refill >= swd_below_sat) then 
	            soilw(i) = wsat(i)
                  refill = refill - swd_below_sat
		    ! ...else insufficient refill available to fill layer to saturation
		    else  
	            soilw(i) = soilw(i) + refill
	            refill = 0.0
	            exit  ! No refill left so quit loop early
	        endif
	    enddo  ! i (layers)
      endif
      
	end subroutine drain_sundial_water2
C
C-------------------------------------------------------------------------------
C
      SUBROUTINE DRAIN_SUNDIAL_WATER(RAIN,WMAX,SOILW,DRAINW,REFIL,
     &                               WSAT,FLOWPROP)
C
C Subroutine to calculate the water drainage
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXSOIL			! Max.no.of soil types
      PARAMETER(MAXSOIL=50)
	INTEGER MAXLAYER		! No.of layers in the soil profile
	INTEGER MAXLAYER1		! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/
	INTEGER IL				! Local layer counter
	REAL EXCESS				! Excess drainage to achieve water table (mm)
C
C Variables passed from calling subroutine
C
	REAL WMAX(MAXLAYER)		! IN:Available water at field cap. (mm/layer)
	REAL SOILW(MAXLAYER)	! IN:Available water (mm/layer)
	REAL RAIN				! IN: Rainfall (mm/timestep)
	REAL DRAINW(MAXLAYER)	! IN:Water drained from this layer (mm/layer)
	REAL REFIL(MAXLAYER)	! Water deficit in the layer mm/layer
	REAL FLOWPROP			! IN:Proportion of flow needed to achieve
							!	     water table at depth WTABLE
	REAL WSAT(MAXLAYER)		! IN:Available water at saturation (mm/layer)
C
C If no rainfall do not allow any leaching
C
      IF(RAIN.LT.0.0001)THEN
	  DO 100 IL=1,MAXLAYER1
          DRAINW(IL)=0.0
100     CONTINUE
        GOTO 200
      END IF
C
C For each layer in profile...
C
      DO 300 IL=1,MAXLAYER1
C
C Set water deficit
C
        REFIL(IL)=WMAX(IL)-SOILW(IL)
	  IF(REFIL(IL).LT.0)REFIL(IL)=0
C
C If rainfall is less than the water deficit
C add all rainfall to the current layer - no drainage
C
        IF(RAIN.LE.REFIL(IL))THEN
          SOILW(IL)=SOILW(IL)+RAIN
          RAIN=0
          DRAINW(IL)=0
C
C If rainfall exceeds the water deficit allow drainage
C
        ELSE
C
C ...Fill layer to available water
C
          SOILW(IL)=WMAX(IL)
C
C ...Take the amount of water left behind off the rainfall
C
	    RAIN=RAIN-REFIL(IL)
C
C ... The amount drained is what is left from the rainwater
C
		DRAINW(IL)=RAIN
        ENDIF
300     CONTINUE
C
C ...Restrict drainage and allow layers to fill up to saturation from bottom of the profile
C
      EXCESS=0
	DO 400 IL=1,MAXLAYER1
	  EXCESS=EXCESS+DRAINW(IL)*(1-FLOWPROP)
400   CONTINUE
      DO 500 IL=MAXLAYER1,1,-1
        DRAINW(IL)=DRAINW(IL)-EXCESS
	  IF(EXCESS.GT.(WSAT(IL)-SOILW(IL)))THEN
          DRAINW(IL)=DRAINW(IL)+EXCESS-(WSAT(IL)-SOILW(IL))
	    EXCESS=EXCESS-(WSAT(IL)-SOILW(IL))
		SOILW(IL)=WSAT(IL)
	  ELSEIF(EXCESS.LE.(WSAT(IL)-SOILW(IL)))THEN
		SOILW(IL)=SOILW(IL)+EXCESS
	    EXCESS=0
	  ENDIF
500   CONTINUE 
200   CONTINUE
      END
C
C-------------------------------------------------------------------------------
C
      subroutine evap_sundial_water2(evapw, airtemp, icover, wmax,
     &                               soilw)

    ! This routine is the same as evap_sundial_water() except that it has been
    ! optimised to run faster because it is now used on a daily timestep in
    ! limited data and spatial ecosse simulations

      ! Constants
      integer, parameter :: maxlayer = 60  ! No. of layers in the soil profile
      integer, parameter :: maxdepth = 300 ! Depth of soil profile [cm]

      ! Arguments with intent in
      real, intent(in) :: airtemp         ! Mean air temperature [deg C]
      real, intent(in) :: evapw           ! Potential evapotranspiration [mm]
      integer, intent(in) :: icover       ! 1 if land is covered with vegetation, otherwise 0
      real, intent(in) :: wmax(maxlayer)  ! Available water at field capacity [mm/layer]

      ! Arguments with intent inout
      real, intent(inout) :: soilw(maxlayer)  ! Available water [mm/layer]

      ! Local variables
      real :: availw                      ! The water available for evaporation in a layer [mm/layer]
      real :: evap                        ! Potential evapotranspiration [mm]
      integer :: i                        ! Layer counter
      integer :: idepth                   ! Depth of bottom of soil layer [cm]

      evap = evapw

      if (evap > 0) then  ! Protect against negative evap
        ! If not freezing and crop present calculate water lost by evapotranspiration
        if (airtemp >= 0.0 .and. icover == 1) then
            do i=1, maxlayer
                idepth = i * maxdepth / maxlayer  ! Calculte depth to bottom of layer

                ! In 0-50cm layers any available water can be lost by evapotranspiration
                if (idepth <= 50) then
                    if (evap < soilw(i)) then
                        soilw(i) = soilw(i) - evap
                        evap = 0
                        exit  ! No evaporation left so quit layer loop early
                    else
                        evap = evap - soilw(i)
                        soilw(i) = 0
                    endif
                ! If the soil moisture now exceeds half the capacity of the 50-100cm layer 
                ! the excess is lost by evapouration as far as possible
                elseif (idepth > 50 .and. idepth < 100) then
                    availw = wmax(i) * 0.50
                    if (soilw(i) > availw) then
                        if (evap > soilw(i) - availw) then
                            soilw(i) = availw
                            evap = evap - (soilw(i) - availw)
                        else
                            soilw(i) = soilw(i) - evap
                            evap = 0
                            exit  ! No evaporation left so quit layer loop early
                        endif
                    endif
                ! If the soil moisture now exceeds 75% the capacity of the 100-150cm layer the
                ! excess is lost by evapouration as far as possible
                else  ! idepth must be > 100)
                    availw = wmax(i) * 0.75
                    if (soilw(i) > availw) then
                        if (evap > soilw(i) - availw) then
                            soilw(i) = availw
                            evap = evap - (soilw(i) - availw)             
                        else
                            soilw(i) = soilw(i) - evap
                            evap = 0
                            exit  ! No evaporation left so quit layer loop early
                        endif
                    endif
                endif
            enddo  ! Layers

        ! If freezing or no crop any available water can be lost by evaporation from the top 5cm only
        else
            soilw(1) = max(soilw(1) - evap, 0.0)
            evap = 0.0
        endif
      endif

      end subroutine evap_sundial_water2
C
C-------------------------------------------------------------
C
   	SUBROUTINE EVAP_SUNDIAL_WATER(EVAPW,AIRTEMP,ICOVER,
     &                              WMAX,SOILW)
C
C Subroutine to calculate the water drainage
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      INTEGER MAXSOIL
      PARAMETER(MAXSOIL=50)
	INTEGER MAXLAYER,MAXDEPTH,IDEPTH,MAXLAYER1
	PARAMETER (MAXLAYER=60)
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
	INTEGER IBARE,IL
	REAL AVAILW,EVAP
C
C Variables passed from calling subroutine
C
      INTEGER ICOVER
	REAL AIRTEMP,WMAX(MAXLAYER),EVAPW,SOILW(MAXLAYER)
	REAL TOTAL_EVAPW
C
C Set evapotranspiration for calculations
C
      EVAP=EVAPW
C
C For each layer in profile...
C
      DO 30 IL=1,MAXLAYER1
C
C Protect against negative EVAP
C
        IF(EVAP.GT.0)THEN
C
C Work out depth
C
          IDEPTH=IL*MAXDEPTH/(MAXLAYER1*1.)
C
C ...If freezing or no crop set IBARE to 1
C
		IF(AIRTEMP.LT.0.0.OR.ICOVER.EQ.0)THEN
            IBARE=1
          ELSE
            IBARE=0
          END IF
C
C ...If not freezing or crop present calculate water lost by evapotranspiration
C
          IF(IBARE.EQ.0)THEN
C
C ......If the soil moisture now exceeds half the capacity of the 50-100cm layer 
C the excess is lost by evapouration as far as possible
C
            IF(IDEPTH.GT.50.AND.IDEPTH.LE.100)THEN
              AVAILW=WMAX(IL)*0.50
              IF(SOILW(IL).GT.AVAILW)THEN
                IF(EVAP.GT.SOILW(IL)-AVAILW)THEN
                  SOILW(IL)=AVAILW
	            EVAP=EVAP-(SOILW(IL)-AVAILW)
                ELSE
                  SOILW(IL)=SOILW(IL)-EVAP
                  EVAP=0
                END IF
              END IF
C
C ......If the soil moisture now exceeds 75% the capacity of the 100-150cm layer the
C excess is lost by evapouration as far as possible
C
            ELSE IF(IDEPTH.GT.100)THEN
              AVAILW=WMAX(IL)*0.75
              IF(SOILW(IL).GT.AVAILW)THEN
                IF(EVAP.GT.SOILW(IL)-AVAILW)THEN
                  SOILW(IL)=AVAILW
	            EVAP=EVAP-(SOILW(IL)-AVAILW)             
                ELSE
                  SOILW(IL)=SOILW(IL)-EVAP
                  EVAP=0
                END IF
              END IF
C 
C ......In 0-50cm layers any available water can be lost by evapotranspiration
C

            ELSE
              IF(EVAP.LT.SOILW(IL))THEN
                SOILW(IL)=SOILW(IL)-EVAP
                EVAP=0
              ELSE
                EVAP=EVAP-SOILW(IL)
                SOILW(IL)=0
              END IF
            END IF
C
C ...If freezing or no crop any available water can be lost by evaporation from the top 5cm
C
          ELSE
            IF(EVAP.GT.SOILW(IL))THEN
              SOILW(IL)=0
            ELSE
              SOILW(IL)=SOILW(IL)-EVAP
            END IF
            EVAP = 0.0
          END IF
	  ENDIF
   30 CONTINUE
C
C Sum of effectivly evapotranspirated water
	TOTAL_EVAPW=EVAPW-EVAP
C	WRITE(34,*) TOTAL_EVAPW	
C
C Leave EVAPOTRANSPIRATION
C
      RETURN
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE INIT_SUNDIAL_WATER(SOILW,ROOT,IROCK,SUM_TS,IDATEFC,
     &                              IANTHES,IRYEAR,IFILE,WMAX,K_TS,RAIN,
     &                              EVAP,AIRTEMP,RDD,TMMN,TMMX,VP,WN,
     &                              IAWC,NSOIL,MODTYPE,SPARMODEL,
     &                              TOC,CLAY,SAND,SILT,BULKDENS,WSAT,
     &                              WTABLE,
     &                              FLOWPROP,AVERAIN,AVEPET,AVETEMP,
     &                              IYEAR,ISTHARV, ISOWNJ,N_STEPS,
     &                              WILTPOINT, SATWATCONT, SECONDS)
C
C Subroutine to initialise soil water 
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXGROW			! Max.no.of growing seasons allowed
	PARAMETER (MAXGROW=300)
	INTEGER MAXSOIL			! Max.no.of soil types
	PARAMETER (MAXSOIL=50)
	INTEGER MAXLAYER		! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER SPARFILE		! IN:Soil parameters read in from file
	INTEGER SPARCALC		! IN:Soil parameters calculated from TOC etc
	DATA SPARFILE,SPARCALC /1,2/
	REAL WILTPOINT_LYR(MAXSOIL,MAXLAYER) ! Soil water content at wilting point (mm/layer)

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
	REAL SECONDS			! IN:Number of seconds in one timestep
	INTEGER SUM_TS			! IN:Total number of timesteps passed
	INTEGER N_STEPS			! IN:No.timesteps in a year
	INTEGER IYEAR			! IN:Current growing season number
	INTEGER ISTHARV			! IN:Timesteps from 01/01/01 to first harvest
	INTEGER ISOWNJ(0:MAXGROW) ! IN:Timesteps from 01/01/01 to sowing date 
C
C ...Soil factors
C
      INTEGER ISAVE			! OUT:Code to save or retrieve variables
	INTEGER IAWC			! IN:Water movement code number
	INTEGER IROCK			! IN:No.layers: 1=50cm, 2=100cm , 3=150cm
      INTEGER NUMSOIL			! IN/OUT: Number of soils defined
	INTEGER NSOIL			! IN:Soil code number
	REAL SOILW(MAXLAYER)	! IN/OUT:Available water (mm/layer)
	REAL WMAX(MAXLAYER)		! IN:Available water at field cap. (mm/layer)
	REAL MAXWAT(MAXSOIL,MAXLAYER)	! IN/OUT: Maximum water content (mm/layer)
	REAL WATSAT(MAXSOIL,MAXLAYER)	! IN:Avail.water at saturation (mm/layer)
	REAL WILTPOINT(MAXLAYER)  ! OUT:Soil water content at wilting point (mm/layer)
	REAL FIELDCAP(MAXLAYER)	  ! OUT:Soil water content at field capacity (mm/layer)
	REAL SATWATCONT(MAXLAYER) ! OUT:Soil water content at saturation (mm/layer)
	REAL WSAT(MAXLAYER)		! IN:Available water at saturation (mm/layer)
	REAL WTABLE				! IN:Water table depth in cm
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
C
C ...Weather factors
C
	INTEGER IFILE			! IN:Current weather file
	REAL RAIN				! IN: Rainfall (mm/timestep)
	REAL EVAP				! IN: Potential evap. (mm/timestep)
	REAL AIRTEMP			! IN: Air temperature (deg.C/timestep)
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
C Set parameters and save for future use
C
      IF(SPARMODEL.EQ.SPARFILE)THEN
        CALL PAR_WAT(NUMSOIL,MAXWAT,SNAME,WATSAT,WILTPOINT_LYR)
	ELSEIF(SPARMODEL.EQ.SPARCALC)THEN
	  CALL MAKE_PAR_WAT(NSOIL,MAXWAT,TOC,CLAY,SAND,SILT,BULKDENS,
     &       WATSAT,WILTPOINT_LYR)
	ENDIF
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
C Calculate soil moisture deficit in mid-July.
C
      CALL DEF(SUM_TS,SOILW,WMAX,IDATEFC,IANTHES,
     &               K_TS,RAIN,EVAP,AIRTEMP,RDD,TMMN,TMMX,VP,WN,
     &               MODTYPE,WSAT,FLOWPROP,IYEAR,ISTHARV, ISOWNJ,
     &				N_STEPS)
C
C Leave INIT_WAT
C
	END
C
C------------------------------------------------------------
C
      SUBROUTINE LIMIT_DRAIN(IMPDEPTH,SOILW)
C
C Subroutine to limit water movement to above the impermeable layer
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXDEPTH			! Maximum depth of the soil profile (cm)
	DATA MAXDEPTH /300/
	INTEGER MAXLAYER			! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)

      REAL DEPTH					! Depth of layer (cm)
	INTEGER IL					! Local layer counter
C
C Variables passed to/from calling subroutine
C
	REAL IMPDEPTH	! IN: Depth of impermeable layer if present (cm)
	REAL SOILW(MAXLAYER)		! IN:Available water (mm/layer)      
C
C Set water content of all layers below the impermeable layer to 0
C Note: fate of water following restricted drainage is dealt with through
C water table depth and FLOWPROP. Here water is simply removed from 
C lower layers to ensure decomposition stops, as the fate of the removed 
C water will depend on factors such as gradient etc. which are encompassed 
C by the value of FLOWPROP
C
      DO 100 IL=1,MAXLAYER
	  DEPTH=IL*MAXDEPTH/MAXLAYER
	  IF(DEPTH.GE.IMPDEPTH)SOILW(IL)=0
100   CONTINUE
C
C Leave LIMIT_DRAIN
C
      END
C
C------------------------------------------------------------
C
C INTERNAL SUBROUTINES
C
C------------------------------------------------------------
C
      SUBROUTINE DEF(SUM_TS,SOILW,WMAX,IDATEFC,IANTHES,
     &               K_TS,RAIN,EVAP,AIRTEMP,RDD,TMMN,TMMX,VP,WN,
     &               MODTYPE,WSAT,FLOWPROP,IYEAR,ISTHARV, ISOWNJ,
     &				N_STEPS)
C
C Subroutine to calculate soil moisture deficit in mid-July.
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
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
C
C Set the variable for summing weeks
C
      SUM_TS=0
C
C Set soil water to max
C
C >>> Commented out by Mark Richards, 31/05/2011
C Commented out so that the initial soil water conditions are at the
C equilbrium achieved by the spin-up carried out in GETEQWATER()
C      DO 1 IL=1,MAXLAYER1
C       SOILW(IL)=WMAX(IL)
C    1 CONTINUE
C End of commented out code by Mark Richards, 31/05/2011 <<<
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
     &   (K_TS.GT.(ISTHARV-(ISTHARV+ISOWNJ(1)-N_STEPS))))THEN
	    IBARE=0
        ELSE
          IBARE=1
        END IF
C
C Evapourate and drain
C
	  CALL DRAIN_SUNDIAL_WATER2(RAIN, WMAX, SOILW, WSAT, FLOWPROP,
     &                            DRAIN, REFIL)
C        CALL DRAIN_SUNDIAL_WATER(RAIN,WMAX,SOILW,DRAIN,REFIL,
C     &                           WSAT,FLOWPROP)
	  CALL EVAP_SUNDIAL_WATER(EVAP,AIRTEMP,IBARE,WMAX,
     &                          SOILW)
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
      SUBROUTINE GETEQWATER(AVERAIN,AVEPET,AVETEMP,
     &                      WMAX,WSAT,LTABLE,
     &                      AVAILWAT,DRAIN,FLOWPROP,

     &                      SECONDS,PI_CEQ_MON)
C
C Subtoutine to get water content at equilibrium
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER		! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1		! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER MAXCYCLES		! Maximum number of spin-up cycles before giving up trying to reach equilibrium!
	PARAMETER (MAXCYCLES=50)
	LOGICAL AT_EQ           ! Flag to indicate whether soil water has reached equilibrium
	INTEGER CYCLES          ! Spin-up cycle counter
	INTEGER LTABLE			! Layer at depth of water table
      INTEGER IL				! Local layer counter
	INTEGER IMON			! Local month counter
      REAL MONTH_SECONDS      ! Number of seconds in one month
	REAL DIFF		        ! Difference between current and previous soil water content [mm/layer]
	REAL RAINTHISMON		! Rainfall this month
	REAL REF(MAXLAYER)		! Refill space in the layer (mm/5cm)
	INTEGER ICOVER			! Crop cover - 1=yes 0=No
	REAL PREV_AVAILWAT(MAXLAYER)  ! Available water in the soil layer in the previous spin-up cycle [mm/layer]
C
C Variables passed to/from this subroutine
C
	REAL AVERAIN(12)		! Long term average rainfall (mm)
	REAL AVETEMP(12)		! Long term average temperature (deg.C)
	REAL AVEPET(12)			! Long term average PET (mm)
	REAL DRAIN(MAXLAYER)	! Drainage from the layer (mm/5cm)
	REAL AVAILWAT(MAXLAYER)	! Available water in the soil layer (mm)
	REAL FLOWPROP			! IN/OUT:Proportion of flow needed to achieve water table at depth WTABLE
      REAL SECONDS			! IN:Number of seconds in one timestep
	REAL WMAX(MAXLAYER)		! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)		! IN:Available water at saturation (mm/layer)
      REAL PI_CEQ_MON(12,MAXLAYER)

      MONTH_SECONDS = (365.25 / 12) * 60 * 60 * 24
C
C Calculate available water after long term average rainfall and evaporation
C Set soil water to field capacity above the water table and saturated below the water table
C
      DO 100 IL=1,LTABLE-1
       AVAILWAT(IL)=WMAX(IL)
100	CONTINUE
      DO 200 IL=LTABLE,MAXLAYER1
       AVAILWAT(IL)=WSAT(IL)
200	CONTINUE
C
C Calculate water deficit from soil 
C
	AT_EQ = .false.
	CYCLES = 0
	DO WHILE(.not.AT_EQ .and. CYCLES < MAXCYCLES)
        PREV_AVAILWAT = AVAILWAT
               DO IMON=1,12
            IF(PI_CEQ_MON(IMON,1).GT.0)THEN
                  ICOVER=1    ! Assume soil covered - results in soil being slightly drier 
            ELSE
                  ICOVER=0
            ENDIF
            
            RAINTHISMON=AVERAIN(IMON)
            CALL RUN_SUNDIAL_WATER_SUBMONTHLY(RAINTHISMON,WMAX,
     &                                        AVAILWAT,WSAT,FLOWPROP,
     &                                        DRAIN,REF,AVEPET(IMON),

     &                                        AVETEMP(IMON),ICOVER)
        ENDDO  ! IMON (months)

	! Compare the soil water in each layer at the end of the cycle to the soil
	! water at the end of the previous cycle to see if all the layers are at
	! equilibrium yet
	  AT_EQ = .TRUE.
	  DO IL=1, MAXLAYER
	    DIFF = ABS(AVAILWAT(IL) - PREV_AVAILWAT(IL))
	    IF (DIFF > 0.1) THEN
	      AT_EQ = .FALSE.
            EXIT
	    ENDIF
	  ENDDO  ! IL (layers)
	  CYCLES = CYCLES + 1
	ENDDO  ! Spin-up cycles
C
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE GETFLOWPROP(DRAIN,LTABLE,WSAT,AVAILWAT,FLOWPROP)
C
C Calculate the flow proportion 
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER		! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1		! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
      INTEGER IL				! Local layer counter	
	INTEGER IMON			! Local month counter
	REAL DRAIN_UNRESTRICTED	! Drainage below observed water table 
							!	assuming no restrictions (mm)
	REAL DRAIN_RESTRICTED	! Drainage below observed water table 
							!	assuming restrictions that result in water table depth (mm)
C
C Variables passed to/from this subroutine
C
	REAL DRAIN(MAXLAYER)	! Drainage from the layer (mm/5cm)
	REAL AVAILWAT(MAXLAYER)	! Available water in the soil layer (mm)
	INTEGER LTABLE			! Layer at depth of water table
	REAL FLOWPROP			! IN/OUT:Proportion of flow needed to achieve		
							!	     water table at depth WTABLE				
	REAL WSAT(MAXLAYER)		! IN:Available water at saturation (mm/layer)
C
C Work out drainage assuming restricted and unrestricted flow
C
      DRAIN_UNRESTRICTED=0
      DRAIN_RESTRICTED=0
      DO 100 IL=1,MAXLAYER1
         DRAIN_UNRESTRICTED=DRAIN_UNRESTRICTED+DRAIN(IL)
         DRAIN_RESTRICTED=DRAIN_RESTRICTED+DRAIN(IL)
	   IF(IL.GE.LTABLE)THEN
		 DRAIN_RESTRICTED=DRAIN_RESTRICTED-(WSAT(IL)-AVAILWAT(IL))
	   ENDIF
100	CONTINUE
C
C Work out flow proportion from the ratio of restricted and unrestricted drainage
C
      IF(DRAIN_UNRESTRICTED.GT.0)THEN
        FLOWPROP=DRAIN_RESTRICTED/DRAIN_UNRESTRICTED
	ELSEIF(DRAIN_UNRESTRICTED.LE.0)THEN
	  IF(DRAIN_RESTRICTED.GT.0)FLOWPROP=0
	  IF(DRAIN_RESTRICTED.LE.0)FLOWPROP=1
	ENDIF
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE GET_LAYERFC(TOC_IL,BULKDENS_IL,CLAY_IL,SILT_IL,
     &                       MAXWAT_IL,IL)
C
C Subroutine to calculate maximum water content
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	REAL FRACOC				! Fraction organic C in the soil
	REAL FTEMP				! Temporary number used in calculation
	INTEGER MAXLAYER1		! No.of layers in the soil profile
	INTEGER MAXDEPTH		! Depth of profile (cm)
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
C
C Variables passed to/from calling subroutine
C
      INTEGER IL				! IN:Layer number
	REAL TOC_IL				! IN:Total organic C (kgC/ha/layer)
      REAL CLAY_IL			! IN:Clay content of the layer (%)
	REAL SILT_IL			! IN:Silt content of the layer (%)
	REAL BULKDENS_IL		! IN:Bulk density of the layer (g/cm3)
	REAL MAXWAT_IL			! OUT: Maximum water content (mm/layer)
C
C Caclulate FC 
C MAXWAT_IL = Maximum available water in the layer (Hall et al, 1977; Soil Survey technical Monograph 9:32-42)
C
	FRACOC=TOC_IL/(BULKDENS_IL*1000000)
c	MAXWAT_IL=47+(0.25*CLAY_IL)+(0.1*SILT_IL)+(1.12*FRACOC)
	MAXWAT_IL=47+(0.25*CLAY_IL)+(0.1*SILT_IL)+(1.12*FRACOC*100)
	MAXWAT_IL=MAXWAT_IL-(16.52*BULKDENS_IL)
	FTEMP=(2.94+0.83*CLAY_IL-0.0054*CLAY_IL*CLAY_IL)
	MAXWAT_IL=MAXWAT_IL-FTEMP
	MAXWAT_IL=MAXWAT_IL*(MAXDEPTH/(MAXLAYER1*1.))/10.
	IF(MAXWAT_IL.LT.0)MAXWAT_IL=32/5	
	IF(IL*(MAXDEPTH/(MAXLAYER1*1.)).GE.80)				
     &  MAXWAT_IL=MAXWAT_IL*1.33			
C
C Temporary test of Lilly's pedotransfer funtion
C
c      MAXWAT_IL=(MAXDEPTH/MAXLAYER1)*(56.1-(7.709*BULKDENS_IL))/25
	END
C
C-------------------------------------------------------------
C
	SUBROUTINE GET_LAYERWSAT(TOC_IL,CLAY_IL,MAXWAT_IL,WATSAT_IL)
C
C Subroutine calculate saturated water content 
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER1		! No.of layers in the soil profile
	INTEGER MAXDEPTH		! Depth of profile (cm)
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
C
C Variables passed to/from calling subroutine
C
	REAL TOC_IL				! IN:Total organic C (kgC/ha/layer)
      REAL CLAY_IL			! IN:Clay content of the layer (%)
	REAL MAXWAT_IL			! IN: Maximum water content (mm/layer)
	REAL WATSAT_IL			! OUT:Available water at saturation (mm/layer)
C
C Calculate saturated water content using Aitkenhead Pedotransfer function
C
	WATSAT_IL=(0.0002*TOC_IL)+(0.25*CLAY_IL)
	WATSAT_IL=WATSAT_IL*(MAXDEPTH/MAXLAYER1)/5.0

	IF(WATSAT_IL.LT.MAXWAT_IL)
     &  WATSAT_IL=MAXWAT_IL/0.8
      END
C
C-------------------------------------------------------------
C
	SUBROUTINE GET_LAYERWSAT_OM(BULKDENS_IL,MAXWAT_IL,WATSAT_IL)
C
C Subroutine calculate saturated water content 
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER1		! No.of layers in the soil profile
	INTEGER MAXDEPTH		! Maximum depth of the soil profile
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
      REAL SOLID_DENSITY		! Density of solids (g/cm)
	DATA SOLID_DENSITY /1.5/	! Reported to vary between 1.4 and 
	                            ! 1.6 g cm-3 (Puustjarvi, 1970, 
								! Degree of humification. 
								! Peat Plant News 3,48-52.)
								! Average taken = 1.5 g/cm
C
C Variables local to this subroutine
C
      REAL WILTPOINT			! Volumetric water content at wilting 
							! point (pressure head =-15,000 cm H2O) (mm)
C
C Variables passed to/from calling subroutine
C
	REAL BULKDENS_IL		! IN:Bulk density of the layer (g/cm3)
	REAL MAXWAT_IL			! IN: Maximum water content (mm/layer)
	REAL WATSAT_IL			! OUT:Available water at saturation (mm/layer)
C
C Calculate saturated water content using water retention properties from
C Weiss et al, 2006. Simulation of water table level and peat temperatures
C in boreal peatlands. Ecological Modelling, 192,441-456.
C
	WATSAT_IL=10*(SOLID_DENSITY-BULKDENS_IL)/SOLID_DENSITY
	WATSAT_IL=WATSAT_IL*MAXDEPTH/MAXLAYER1
      WILTPOINT=WATSAT_IL*0.1 ! Temporary code
	WATSAT_IL=WATSAT_IL-WILTPOINT

C	IF(MAXWAT_IL.GT.WATSAT_IL)
C     &  MAXWAT_IL=WATSAT_IL*0.8
      MAXWAT_IL=WATSAT_IL*0.9
C
C Temporary test of Lilly's pedotransfer funtion
C
c      WATSAT_IL=22.1+(0.348*49.4)-(7.427*BULKDENS_IL)
c      WATSAT_IL=(MAXDEPTH/MAXLAYER1)*WATSAT_IL/25
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE MAKE_PAR_WAT(NSOIL,MAXWAT,TOC,CLAY,SAND,SILT,BULKDENS,
     &                        WATSAT,WILTPOINT_LYR)
C
C Subroutine to make water parameters 
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXSOIL			! Max.no.of soil types
	INTEGER MAXLAYER		! No.of layers in the soil profile
      PARAMETER (MAXSOIL=50,MAXLAYER=60)
	INTEGER MAXLAYER1		! No.of layers in the soil profile
	INTEGER MAXDEPTH		! Depth of profile (cm)
	DATA MAXLAYER1 /60/
	DATA MAXDEPTH /300/
	real :: c_frac          ! Carbon concentration of a layer [fraction]
	integer :: i			! Local counter variable
      logical :: is_peat(maxlayer) ! Flags whether soil is peat or not (1 = peat, 0 = non-peat)
	real :: lyr_depth       ! Depth/thickness of a layer [mm]
	real :: maxsat			! Maximum value for availabel water at saturation per layer [mm/layer]
	real :: orig_watsat     ! Temporary store of the watsat value before modification
      logical :: ptf_failed   ! Flag to indicate if PTF failed, .true. means it produced bad values
	logical :: topsoil      ! True=layer belongs to topsoil (top 30cm), False=subsoil
C
C Variables passed to/from calling subroutine
C
	REAL WATSAT(MAXSOIL,MAXLAYER) ! OUT:Available water at saturation (mm/layer)
	INTEGER NSOIL			! IN:Soil code number
	REAL TOC(MAXLAYER)		! IN:Total organic C (kgC/ha/layer)
      REAL CLAY(MAXLAYER)		! IN:Clay content of the layer (%)
	real SAND(MAXLAYER)     ! IN: Sand content of the layer [%]
	REAL SILT(MAXLAYER)		! IN:Silt content of the layer (%)
	REAL BULKDENS(MAXLAYER)	! IN:Bulk density of the layer (g/cm3)
	REAL MAXWAT(MAXSOIL,MAXLAYER)	! OUT: Maximum water content (mm/layer)
	REAL WILTPOINT_LYR(MAXSOIL,MAXLAYER) ! OUT: Soil water content at wilting point (mm/layer)

    	lyr_depth = maxdepth / maxlayer * 10

      ! Determine if each soil layer is peat or not. This infor is used when deciding which
      ! pedotransfer function to use.
      do i=1, maxlayer
          ! Layer is considered to be peat if C conc > 30% and bulk density < 0.25 g/cm3
		c_frac = toc(i) / (bulkdens(i) * 10000 * lyr_depth)
		if (c_frac > 0.30 .and. bulkdens(i) <= 0.25) then         
              is_peat(i) = .true.
          else
              is_peat(i) = .false.
          endif
      enddo
      
      ! Calculate wiltpoint, field capacity and water at saturation
	do i=1, maxlayer
		if (is_peat(i)) then
			call boelter_ptf(bulkdens(i), maxwat(nsoil,i),
     &                         watsat(nsoil,i), wiltpoint_lyr(nsoil,i),
     &                         lyr_depth)
		! ...for all other soils, use the British Soil Survey PTF
		else
			! If this layer in the topsoil (i.e. top 30 cm)...
			if (i < 7) then
				topsoil = .true.
			! ...else its part of the subsoil
			else
				topsoil = .false.
			endif
	        call bss_ptf(sand(i), silt(i), clay(i), c_frac*100.0, 
     &		             bulkdens(i), lyr_depth, topsoil,
     &					 wiltpoint_lyr(nsoil,i), maxwat(nsoil,i),
     &					 watsat(nsoil,i))
		endif
	end do  ! i (layers)	

	! TEMPORARY FIX TO PTF PROBLEM!
	! It is likely that the the combination of the Boelter PTF (for peats)
	! and the BSS PTF (for non-peats) are not sufficient to predict
	! sensible water retention properties for all soil types. Therefore
	! the following code has been put in to correct nonsensical values and
	! issue a warning to the user.
	maxsat = lyr_depth * 0.9
	do i=1, maxlayer
          ptf_failed = .true.
          if ((watsat(nsoil,i) + wiltpoint_lyr(nsoil,i)) > maxsat) then
			print *,'Warning: PTF predicted saturation > maxsat'
          elseif (maxwat(nsoil,i) > watsat(nsoil,i)) then
			print *,'Warning: PTF gave a field capacity > saturation'
          elseif (wiltpoint_lyr(nsoil,i) > maxwat(nsoil,i)) then
			print *,'Warning: PTF gave a wilt point > field capacity'
          else  ! PTF estimated values that were within "sensibile" bounds
              ptf_failed = .false.
          endif
          if (ptf_failed) then
              print *,'Warning: using default soil water properties'
              if (is_peat(i)) then
                  wiltpoint_lyr(nsoil,i) = lyr_depth * 0.15
                  maxwat(nsoil,i) = 0.65 * lyr_depth
                  watsat(nsoil,i) = 0.9 * lyr_depth
              else  ! Non-peat layer so use approx. vlaues for a loam soil
                  wiltpoint_lyr(nsoil,i) = lyr_depth * 0.1
                  maxwat(nsoil,i) = 0.25 * lyr_depth
                  watsat(nsoil,i) = 0.45 * lyr_depth
              endif
          endif
	enddo  ! i (layers)
	! END OF TEMPORARY FIX

      ! Convert the predicted actual water contents at FC and sat
	! to available water contents because SWATER() will add 
      ! wiltpoint on later
	do i=1, maxlayer
		watsat(nsoil,i) = watsat(nsoil,i) - wiltpoint_lyr(nsoil,i)
		maxwat(nsoil,i) = maxwat(nsoil,i) - wiltpoint_lyr(nsoil,i)
	enddo  ! i (layers)
			
C>>> Temp change JUS 28/02/08 >>>>
C      CALL GET_LAYERWSAT_OM(BULKDENS(IL),MAXWAT(NSOIL,IL),
C     &                     WATSAT(NSOIL,IL))
C<<< Temp change JUS 28/02/08 <<<<
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE PAR_WAT(NUMSOIL,MAXWAT,SNAME,WATSAT,WILTPOINT_LYR)
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
	INTEGER I,J,N,IL			! Local counter variables
	REAL FNULL
	REAL TOCARRAY(MAXSOIL,MAXLAYER)	! Equilibrium TOC (kgC/ha/layer)
	REAL CLARRAY(MAXSOIL,MAXLAYER)	! % clay content in this layer
	REAL WILTPOINT25(MAXSOIL,MAXLAYER)	! Soil water content at wilting point (mm/0-25cm)
C
C Variables passed to/from calling subroutine
C
	REAL WATSAT(MAXSOIL,MAXLAYER) ! OUT:Available water at saturation (mm/layer)
      INTEGER NUMSOIL				! OUT:Number of soils defined
	CHARACTER*40 SNAME(MAXSOIL)	! OUT:Soil name
	REAL MAXWAT(MAXSOIL,MAXLAYER)	! IN/OUT:Maximum water content (mm/layer)
	REAL WILTPOINT_LYR(MAXSOIL,MAXLAYER) ! OUT: Soil water content at wilting point (mm/layer)
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
	WILTPOINT_LYR(I,1)=WILTPOINT25(I,1)*(MAXDEPTH/(MAXLAYER1*1.))/25.
	DO 200 IL=2,MAXLAYER1
	  MAXWAT(I,IL)=MAXWAT(I,1)
	  WATSAT(I,IL)=WATSAT(I,1)
	  WILTPOINT_LYR(I,IL)=WILTPOINT_LYR(I,1)
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
	READ(46,20,ERR=222)MAXWAT(I,1),WATSAT(I,1)
	READ(46,*,ERR=222)WILTPOINT25(I,1)
	MAXWAT(I,1)=MAXWAT(I,1)*(MAXDEPTH/(MAXLAYER1*1.))/25.
	WATSAT(I,1)=WATSAT(I,1)*(MAXDEPTH/(MAXLAYER1*1.))/25.
	WILTPOINT_LYR(I,1)=WILTPOINT25(I,1)*(MAXDEPTH/(MAXLAYER1*1.))/25.
c20    FORMAT(F9.1/,16/,F9.1)
20    FORMAT(F9.1/,17/,F9.1)
	DO 100 IL=2,MAXLAYER1
	  MAXWAT(I,IL)=MAXWAT(I,1)
	  WATSAT(I,IL)=WATSAT(I,1)
	  WILTPOINT_LYR(I,IL)=WILTPOINT_LYR(I,1)
100   CONTINUE
      SNAME(I)=TEMP
      NUMSOIL=I
      REWIND(46)
	RETURN
      PRINT*,'SOIL: Water parameters from SOIL.DAT'
333   CONTINUE
222   CONTINUE
      WRITE(*,*)'Warning! Error in general soil parameters!'
      WRITE(*,*)'Check format of soil parameter file'
      WRITE(15,*)'Warning! Error in general soil parameters!'
      WRITE(15,*)'Check format of soil parameter file'
C
C Check value of WATSAT and correct using Aitkenhead Pedotransfer function
C
      DO 400 I=1,NUMSOIL
        DO 300 IL=1,MAXLAYER1
          IF(WATSAT(I,IL).LT.MAXWAT(I,IL))THEN
            CALL GET_LAYERWSAT(TOCARRAY(I,IL),CLARRAY(I,IL),MAXWAT(I,IL)
     &                        ,WATSAT(I,IL))
	    ENDIF
300     CONTINUE
400	CONTINUE
C
C Leave RESPAR
C
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE RESTRICT_DRAINAGE(FLOWPROP,WTABLE,WSAT,WMAX,PI_CEQ_MON,
     &                             AVERAIN,AVEPET,AVETEMP,AVAILWAT,
     &                             SECONDS)
C
C Calculate the proportion restriction occuring given the minimum water table depth
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXLAYER		! No.of layers in the soil profile
	PARAMETER (MAXLAYER=60)
	INTEGER MAXLAYER1		! No.of layers in the soil profile
	DATA MAXLAYER1 /60/
	INTEGER MAXDEPTH		! Maximum depth of the soil profile
	DATA MAXDEPTH /300/
	REAL AVERAIN(12)		! Long term average rainfall (mm)
	REAL AVETEMP(12)		! Long term average temperature (deg.C)
	REAL AVEPET(12)			! Long term average PET (mm)
	REAL FLOWPROP_HI     	! High flowprop value used in binary search algorithm
	REAL FLOWPROP_LO        ! Low flowprop value used in binary search algorithm
	REAL PREV_FLOWPROP   	! Previous flowprop value (used in binary search algorithm
	INTEGER LTABLE			! Layer at depth of water table
      INTEGER IL				! Local layer counter	
	INTEGER IMON			! Local month counter
	REAL DRAIN(MAXLAYER)	! Drainage from the layer (mm/5cm)
	REAL AVAILWAT(MAXLAYER)	! Available water in the soil layer (mm)
      INTEGER LTABLESIM		! Simulated layer containing the water table
	REAL MONTH_SECONDS      ! Number of seconds in one month
	REAL TSTEPS_IN_MONTH    ! Number of timesteps in one month
C
C Variables passed to/from this subroutine
C
	REAL FLOWPROP			! IN/OUT:Proportion of flow needed to achieve		
							!	     water table at depth WTABLE				
	REAL SECONDS			! IN:Number of seconds in one timestep
	REAL WMAX(MAXLAYER)		! IN:Available water at field cap. (mm/layer)
	REAL WSAT(MAXLAYER)		! IN:Available water at saturation (mm/layer)
	REAL WTABLE				! IN:Water table depth in cm
      REAL PI_CEQ_MON(12,MAXLAYER)

C
C Initialise flowprop variables for start of binary search algorithm
C
	FLOWPROP_HI = 1.0
	FLOWPROP_LO  = 0.0
C
C Work out layer containing the water table
C
      LTABLE=1+INT(WTABLE*MAXLAYER1/MAXDEPTH)
	IF(LTABLE.GT.MAXLAYER1)LTABLE=MAXLAYER1
      FLOWPROP = 0.5  ! Assign an initial "middle" value for FLOWPROP

202   CONTINUE
      CALL GETEQWATER(AVERAIN,AVEPET,AVETEMP,
     &                        WMAX,WSAT,LTABLE,
     &                        AVAILWAT,DRAIN,FLOWPROP,

     &                        SECONDS,PI_CEQ_MON)

      LTABLESIM=MAXLAYER1
      DO 200 IL=1,MAXLAYER1 
	  IF(LTABLESIM.EQ.MAXLAYER1.AND.AVAILWAT(IL).EQ.WSAT(IL))
     &    LTABLESIM=IL
	  IF(LTABLESIM.LT.IL.AND.AVAILWAT(IL).LT.WSAT(IL))
     &    LTABLESIM=MAXLAYER1
200   CONTINUE

C >>> Added by Mark Richards, 19/05/2011
C Binary search algorithm to find the appropriate value of FLOWPROP. This works by 
C iteratively checking whether the current value of FLOWPROP gives a water table (WT)
C higher or lower than the target WT. If the WT is higher (i.e. a lower soil layer 
C number) then the next value of FLOWPROP is half way between the previous upper
C estimate (FLOWPROP_HI) and the current FLOWPROP value. If the WT is lower then the
C next value of FLOWPROP will be half way between the previous lower estimate 
C (FLOWPROP_LO) and the current value of FLOWPROP. FLOWPROP_HI and FLOWPROP_LO are set
C to 1 and 0 respectively.
C The algorithm exits when either the target WT has been achieved or FLOWPROP
C changes by less than a given threshold value (in which case it is assumed the closest
C possible WT position to the target has been acheived).
C 
	IF (LTABLESIM < LTABLE) THEN  ! Simulated WT is too high so increase FLOWPROP
	  ! Increase flow prop
	  PREV_FLOWPROP = FLOWPROP
	  FLOWPROP_LO = FLOWPROP
	  FLOWPROP = (FLOWPROP_HI + FLOWPROP) / 2.0
	  IF (FLOWPROP - PREV_FLOWPROP < 0.0000001) THEN
	    GOTO 203
	  ENDIF
      ELSEIF (LTABLESIM > LTABLE) THEN  ! Simulated WT is too low so decrease FLOWPROP
	  PREV_FLOWPROP = FLOWPROP
	  FLOWPROP_HI = FLOWPROP
	  FLOWPROP = (FLOWPROP_LO + FLOWPROP) / 2.0
	  IF (PREV_FLOWPROP - FLOWPROP < 0.0000001) THEN
	    GOTO 203
	  ENDIF
	ELSEIF (LTABLESIM == LTABLE) THEN
	  GOTO 203
      ENDIF
      GOTO 202
203   CONTINUE
     
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE SAVE_WAT(NUMSOIL,MAXWAT,SNAME,ISAVE,WATSAT)
C
C Subroutine to save water parameters for future use
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER MAXSOIL			! Max.no.of soil types
	INTEGER MAXLAYER		! No.of layers in the soil profile
      PARAMETER (MAXSOIL=50)
	PARAMETER (MAXLAYER=60)
	INTEGER I,J				! Local counter variables
      INTEGER NUMSOIL1		! Number of soils defined
	CHARACTER*40 SNAME1(MAXSOIL) ! Soil name
	REAL MAXWAT1(MAXSOIL,MAXLAYER)	! Maximum water content (mm/layer)
      REAL WATSAT1(MAXSOIL,MAXLAYER)	! Avail.water at saturation (mm/layer)
C
C Variables passed to/from calling subroutine
C
      INTEGER NUMSOIL			! IN/OUT: Number of soils defined
      INTEGER ISAVE			! IN:Code to save or retrieve variables
	CHARACTER*40 SNAME(MAXSOIL)		! IN/OUT: Soil name
	REAL MAXWAT(MAXSOIL,MAXLAYER)	! IN/OUT: Maximum water content (mm/layer)
	REAL WATSAT(MAXSOIL,MAXLAYER)	! IN:Avail.water at saturation (mm/layer)
C
C Save parameters
C
      SAVE
	NUMSOIL=MAXSOIL
	IF(ISAVE.EQ.1)THEN
	  NUMSOIL1=NUMSOIL
	  DO 100 I=1,NUMSOIL1
	    SNAME1(I)=SNAME(I)
		DO 200 J=1,MAXLAYER
	      MAXWAT1(I,J)=MAXWAT(I,J)  
	      WATSAT1(I,J)=WATSAT(I,J)  
200       CONTINUE
100     CONTINUE
C
C Retrieve parameters
C
      ELSEIF(ISAVE.EQ.0)THEN
	  NUMSOIL=NUMSOIL1
	  DO 300 I=1,NUMSOIL
	    SNAME(I)=SNAME1(I)
	    DO 400 J=1,20
	      MAXWAT(I,J)=MAXWAT1(I,J)  
	      WATSAT(I,J)=WATSAT1(I,J)  
400       CONTINUE
300     CONTINUE
      ENDIF
C
C Leave SAVE_WAT
C
      END
C
C---------------------------------------------------------
C
      SUBROUTINE SWATER(IAWC,WMAX,MAXWAT,WSAT,WATSAT,FIELDCAP,
     &				  WILTPOINT_LYR,WILTPOINT,SATWATCONT)
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
	INTEGER IL				! Local counter variable
C
C Variables passed to/from calling subroutine
C
	REAL FIELDCAP(MAXLAYER)	! OUT:soil water content at field capacity (mm/layer)
	INTEGER IAWC			! IN:Water movement code number
	REAL WMAX(MAXLAYER)		! IN:Available water at field cap. (mm/layer)
	REAL MAXWAT(MAXSOIL,MAXLAYER)	! IN/OUT: Maximum water content (mm/layer)
	REAL WILTPOINT_LYR(MAXSOIL,MAXLAYER) ! IN:Water conent at wilting point (mm/layer)
	REAL WILTPOINT(MAXLAYER)	! OUT:Water conent at wilting point (mm/layer)
	REAL SATWATCONT(MAXLAYER)	! OUT:Total water content at saturation (mm/layer)
	REAL WATSAT(MAXSOIL,MAXLAYER)	! IN:Avail.water at saturation (mm/layer)
	REAL WSAT(MAXLAYER)		! IN:Available water at saturation (mm/layer)

C
C Set values of available water 
C
      DO 100 IL=1,MAXLAYER1
        WMAX(IL)=(MAXWAT(IAWC,IL))
        WSAT(IL)=(WATSAT(IAWC,IL))
	  WILTPOINT(IL)= WILTPOINT_LYR(IAWC,IL)
	  FIELDCAP(IL)=WMAX(IL)+WILTPOINT(IL)
	  SATWATCONT(IL)=WSAT(IL)+WILTPOINT(IL)
100   CONTINUE
C
C Leave SWATER
C
      RETURN
      END
C
C-------------------------------------------------------------
C
      SUBROUTINE WATFIX(IROCK,ROOT,NSOIL)
C
C Subroutine to set the available water at field capacity and
C the rooting depth
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      INTEGER MAXSOIL,IL
      PARAMETER (MAXSOIL=50)
C
C Variables passed to/from calling subroutine
C
      INTEGER NSOIL,IROCK
      REAL ROOT
C
C Set rooting depth according to rooting restriction
C
      IF(IROCK.EQ.1)THEN
        ROOT=50
      ELSE IF(IROCK.EQ.2)THEN
        ROOT=100
      ELSE
        ROOT=150
      END IF
C
C Leave WATFIX
C
      RETURN
      END
C
C-------------------------------------------------------------------------------
C
      subroutine boelter_ptf(bd, fc, sat, wp, dz)

	! Estimates the water retention properties of peat soils using the
	! pedotransfer functions of Boelter (1969).
	!
	! Boelter, D.H. (1969) Physical properties of peats as related to degree
	!     of decomposition. Soil sci. soc. amer. proc. vol 33, 1969

	implicit none

	! Arguments with intent(in)
	real, intent(in) :: bd   ! Bulk density of each layer [g/cm3]
	real, intent(in) :: dz   ! Depth (thickness) of layer [mm]

	! Arguments with intent(out)
	real, intent(out) :: sat ! Available water at saturation [mm/layer]
	real, intent(out) :: fc  ! Available water at field capacity [mm/layer]
	real, intent(out) :: wp  ! Soil water content at wilting point [mm/layer]

	! Calculate wilt point of layer (based on Boelter's equation for 15 bar)
	wp = 1.57 + (115.28 * bd) - (107.77 * bd**2)
	wp = wp / 100.0 * dz  ! Convert % vol to mm/lyr

	! Calculate field capacity of layer (based on Boelter's equation for 0.1 bar)
	fc = 2.06 + (719.35 * bd) - (1809.68 * bd**2)
	fc = fc / 100.0 * dz  ! Convert % vol to mm/lyr

      ! Calculate saturated water content of layer
	sat = 99.0 - (123.45 * bd) + (252.92 * bd**2)
	sat = sat / 100.0 * dz  ! Convert % vol to mm/lyr

	end subroutine  ! boelter_ptf()
!
!-------------------------------------------------------------------------------
!
	subroutine bss_ptf(sand,silt,clay,oc,bd,depth,top,wp,fc,sa)
	
	implicit none
	
	! For field soils, use equation of British Soil Service to find wilt point, field capacity and saturation water content
	! from sand, silt, clay, carbon and bulk density data [H92]
	! Equation seems to perform well [D04, G04]
	! Updated values obtained from LEACHM February 2011 revision, provided by John Hutson [H11]

	! eqn of the form: SWC_p = a + b*clay + c*silt + d*OC + e*BD
	! where SWC_p is soil water capacity at pressure p, a is a constant
	! b is coefficient for % clay, c for silt, d for organic carbon and e for bulk density
	! change of form for equation in subsoil at high pressures, as noted below

	! At pressures: [-1500., -200., -40., -10., -5.]kPa
	! For topsoil:
	! a = [0.0611,0.0938,0.2668,0.4030,0.4981]
	! b = [0.004,0.0047,0.0039,0.0034,0.0027]
	! c = [0.0005,0.0011,0.0013,0.0013,0.0011]
	! d = [0.005,0.0069,0.0046,0.004,0.003]
	! e = [0.,0.,-0.0764,-0.125,-0.1778]

	! For subsoil*:
	! a = [0.0125,0.0431,0.2205,0.3086,0.4216]
	! b = [0.0092,0.0108,0.0047,0.0040,0.0034]
	! c = [-0.000062,-0.000079,0.0020,0.0021,0.0018]
	! d = [0.,0.,0.0093,0.0126,0.0022]
	! e = [0.,0.,-0.0956,-0.1246,-0.1697]

	! *in subsoil at presssures -1500 and -200kPa, equation form changes to:
	! SWC_p = a + b*clay + c*clay^2

	! We are interested in p = -1500kPa and 0kPa for wilt point and saturation
	! and a range of pressures for field capacity [R44]:
	! -1kPa (>5% OC), -10kPa (loamy), -100kPa (heavy clay) and -5kPa for typical UK field soil
	! loamy and heavy clay properties defined by soil texture triangle [H12]

	! To obtain coefficients for equation at small p, simply extend line
	! Perform similar for coefficients between any other values.  Similar to procedure used by Donatelli et al. [D04]

	! For field soils saturation (0kPa), rather than use PTF (which doesn't apply at very low pressure), find total pore space per unit volume
	! The fraction f of volume occupied by solid matter is roughly (bulk density)/(particle density), hence saturation is given by 1-f
	! Particle density is assumed 2.65g/cm3 for mineral fractions, and 1.5g/cm3 for organic fractions [P70,R73]

	! References:
	! [D04] M. Donatelli, J.H.M. Wösten, G. Belocchi (2004). Methods to evaluate pedotransfer functions
	!		Developmennt of Pedotransfer Functions in Soil Hydrology, Developments in Soil Science 30, 357-411
	! [G04] Givi et al. (2004), Agricultural Water Management 70, 83-96
	! [H92] Hutson et al. (1992), LEACHM manual (http://www.apesimulator.org/help/utilities/pedotransfers/British_Soil_Service_(topsoil).html)
	! [H11] Hutson (2011), LEACHM manual update
	! [H12] HWSD Documentation 2012
	! [P70] Puustjarvi (1970), Peat Plant News 3,48-52
	! [R44] Richards and Weaver (1944), J Ag Res 215-235
	! [R73] Reid (1973), Area 5(1) 10-12

	real,intent(in):: sand,silt,clay,oc,bd,depth	! values for the current layer [%,%,%,%,g/cm3,mm]
	logical,intent(in):: top ! true if topsoil, else false
	real,intent(out):: wp,fc,sa	! wilt point, field capacity, saturation for layer [mm]

	integer:: i,j,k,N	! loop indices and number of grid points
	integer:: nullv=999	! null value
	integer:: np	! number of pressures of interest for British Soil Service eqns
	real:: wpfcsa(3)	! wilt point, field capacity, saturation [vol/vol, then converted to mm]
	real :: p(5)	! pressures of interest (wilt point, all possible field capacities, saturation)
	real :: coef(5,5)	! coefficients a-e (rows) at all pressure of interest (columns) for layer
	integer:: np0,nc0	! number of input pressures and coefficients
	real :: coef0(5,5)	! coefficients at original defined pressures in layer
	real :: p0(5)	! original defined pressures
	real :: m(5,5),c(5,5)	! coefficient and constant for straight lines between each pressure for each coefficient in layer
	integer:: vi(2)		! index of coefficient vectors for wp and fc

	! constants ----------------------------------------------------------------
	! field soils
	np0 = 5
	nc0 = 5
	p0 = [-1500., -200., -40., -10., -5.]
	if (top) then
		coef0(1,:) = [0.0611,0.0938,0.2668,0.4030,0.4981]
		coef0(2,:) = [0.004,0.0047,0.0039,0.0034,0.0027]
		coef0(3,:) = [0.0005,0.0011,0.0013,0.0013,0.0011]
		coef0(4,:) = [0.005,0.0069,0.0046,0.004,0.003]
		coef0(5,:) = [0.,0.,-0.0764,-0.125,-0.1778]
	else
		coef0(1,:) = [0.0125,0.0431,0.2205,0.3086,0.4216]
		coef0(2,:) = [0.0092,0.0108,0.0047,0.0040,0.0034]
		coef0(3,:) = [-0.000062,-0.000079,0.0020,0.0021,0.0018]
		coef0(4,:) = [0.,0.,0.0093,0.0126,0.0022]
		coef0(5,:) = [0.,0.,-0.0956,-0.1246,-0.1697]
	end if
	! define pressures of interest (n.b. if you change these, you must also change the method to define vi, which must correspond to indices)
	np = 5
	p = [-1500.,-100.,-10.,-5.,-1.]	! wilt point pressure first, then possible field capacities
	vi(1) = 1		! i.e. wp always -1500kPa
	! end constants ------------------------------------------------------------

	! find straight line between each pressure for each coefficient and layer
	do i = 1,nc0
		do j = 1,np0-1
			m(i,j) = (coef0(i,j+1)-coef0(i,j))/(p0(j+1)-p0(j))
			c(i,j) = coef0(i,j)-p0(j)*m(i,j)	! -ve as it's really +(0-p0(j))
		end do
	end do
	! account for change of form for equations in subsoil at -200kPa and below
	if (.not.top) then
		m(3,2) = m(3,3)	! assume equation form changes at exactly -200kPa, hence line from -200 to -40 is same as -40 to -10kPa
		c(3,2) = c(3,3)
	end if

	! calculate coefficients for pressures of interest
	do i = 1,np
		do j = 1,nc0
			do k = 1,np0-1	! find which straight line to use
				if ((p(i)>p0(k).and.p(i)<=p0(k+1)).or.(p(i)<=p0(1)
     &                       .and.k==1).or.(p(i)>p0(np0).and.
     &                        k==np0-1)) then
					coef(j,i) = m(j,k)*p(i)+c(j,k)	! use line between points, or extend line at either end if smaller/larger than quoted range
				end if
			end do
		end do
	end do

	! determine index for p at field capacity from soil properties
	vi(2) = 4	! default index for -5kPa
	if ((clay<=25 .and. clay>20 .and. silt<=50 .and. silt>=28).or.
     &    (clay<=20 .and. clay>=8 .and. silt<=50 .and.
     &     sand<=52)) vi(2) = 3 ! loamy -10kPa
	if (clay>=60) vi(2) = 2	! heavy clay -100kPa
	if (oc>5) vi(2) = 5	! organic soil -1kPa

	! calculate WP, FC, sat (eqns assume BD in g/cm3=kg/dm3, all others as % (not fractions); results as fractions (vol/vol))
	do j = 1,3
		if (j==3) then	! saturation found from bulk density and estimated particle density
			wpfcsa(j) = 1-BD/(2.65*(1-OC/100)+1.5*OC/100)
		else	! use British Soil Service PTF
			wpfcsa(j) = coef(1,vi(j)) + coef(2,vi(j))*clay + 
     &                    coef(3,vi(j))*silt + coef(4,vi(j))*OC + 
     &                    coef(5,vi(j))*BD
			if (.not.top .and. p(vi(j))<=-200) then	! low pressures in subsoil
				wpfcsa(j) = wpfcsa(j) + coef(3,vi(j))*((clay**2)-silt)	! replacing silt term by clay^2 term for -1500 to -200kPa in subsoil
			end if
		end if
	end do

	! convert from vol/vol to mm/ref_depth by multiplying by the layer depth
	do j = 1,3
		wpfcsa(j) = depth*wpfcsa(j)
	end do

	! set any nulls in either layer to zero
	do j = 1,3
		if (sand==nullv) then
			wpfcsa(j) = 0
		end if
	end do

	wp = wpfcsa(1)
	fc = wpfcsa(2)
	sa = wpfcsa(3)

	return
	end subroutine bss_ptf
