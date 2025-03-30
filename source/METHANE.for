!===============================================================================
! Soil Methane Subroutine
! Author: Mark Richards, University of Aberdeen
! Created: September, 2010
!
! External subroutines:
!     methane_richards() - calculates methane production, oxidation & net flux
!
! Internal subroutines:
!     andec() - calculates anaerobic decomposition and reduces CH4 production
!               under N-limiting conditions
!
! Description
! -----------
! Calculates the net methane flux during a timestep.
! Net flux = CH4 produced - CH4 oxidised
!
! The gross CH4 produced in each soil layer is calculated based on the amount of 
! CO2 produced from aerobic decomposition and the CH4:CO2 ratio. The CH4:CO2 is 
! weighted by the water filled pore space (WFPS) of the layer (as a proxy for 
! the degree of anoxia).
!
! The anaerobic decomposition (i.e. the change in C & N pool sizes) is
! calculated from the gross CH4 production. If there is insufficient N
! available to maintain the stable N:C ratio of the BIO and HUM pools then
! the CH4 production (and therefore anaerobic decomposition) is decreased 
! until the N-demand can be met.
!
! Methane from the soil and atmoshere are both subject to oxidation.
! Oxidatoin takes place in the "oxidation zone" which is defined the top 30 cm
! of the soil profile.
! Oxidation of soil derived CH4 is a function of the total CH4 produced during
! the time step and the mean soil temperature, mean WFPS and mean inorganic N
! concentration of the "oxidation zone". 
!
! Oxidation of methane fromt thr atmosphere is a function of the maximum
! atmospheric methane oxidation rate and the soil T, WFPS and inorganic N
! conentration and pH in the oxidation zone. The pH and soil T modifiers are not
! applied to the oxidation of soil derived methane since a pH modifier is already 
! applied to the aerobic decomposition rate that controls CO2 production, which in
! turn effects the CH4 production and hence the amount of CH4 oxidation taking place.
!
! Note, the net methane flux can be negative (i.e. the soil is a methane sink).
!
! References
! ----------
! Chan, A.S.K & Parkin, T.B. (2001) Methane oxidation and productivity activity 
!     soils from natural and agricultural ecosystems. J. Environ. Qual.
!     30:1896-1903
! Jang, I., Lee, S., Hong, J., Kang, H. (2006) Methane oxidation rates in 
!     forest soils and their controlling variables: a review and a case study
!     in Korea. Ecol. Res. 21:849-854
! Saari, A., Rinnan, R. & Martikainen, P.J. (2004) Methane oxidation in boreal
!     forest soils: kinetics and sensitivity to pH and ammonium. Soil Biology
!     & Chemisty 36:1037-1046
! Smith, K.A., Dobbie, K.E., Ball, B.C. et al (2000) Oxidation of atmospheric
!    methane in Northern European soils, comparison with other ecosystems, and
!    uncertainties in the global terrestial sink. Global Change Biology 6:791-803
!===============================================================================

      subroutine methane_richards(maxlayer, maxdepth, measlay, maxsoil, 
     &							ch4_prod, ch4_soil_ox, ch4_atmos_ox,
     &							ch4_flux, soilt, asw, sw_wilt,
     &                            asw_fc, asw_sat, sw_sat, bulkdens,
     &                            co2, co2_dpm, co2_rpm,
     &                            co2_hum, co2_bio, dpm_c, rpm_c,
     &                            bio_c, hum_c, dpm_n, rpm_n, bio_n,
     &                            hum_n, nh4_n, no3_n, nsoil,
     &                            biohum_n2c, crit,
     &                            alpha, beta, delta, gamma, ph,
     &                            seconds)
	
      ! Calculates CH4 production, oxidation and net methane flux for a timestep.
      
      implicit none

      ! Local constants
      real, parameter :: bd_max = 0.7                ! Bulk density at which rate modifier is at its maximum [g/cm3]
      real, parameter :: bd_min = 1.4                ! Bulk density at which rate modifier is at is minimum [g/cm3]
      real, parameter :: bd_rate_min = 0.05          ! Minimum value of the bulk density oxidation rate modifier [proportion]
      real, parameter :: ch4co2_ratio_fc = 0.0001     ! Maximum CH4:CO2 ratio at field capacity
      real, parameter :: ch4co2_ratio_sat = 1.67     ! Maximum CH4:CO2 ratio at saturation from Table 3 in Segers et al (1998)
      real, parameter :: daysecs = 86400             ! Number of seconds in a day
      real, parameter :: max_daily_atmos_ox = 0.0671 ! Maximum atmospheric methane oxidation rate [kgC/ha/day] from Jang et al, 2006
      real, parameter :: max_ox_frac = 0.95          ! Maximum proportion of CH4 produced that can be oxidised [prop]
      real, parameter :: n_rate_min = 0.05           ! Minimum value the inorganic N conc oxidation rate modifier can take [proportion]
      real, parameter :: ox_depth = 30               ! Depth (lower boundary) of the oxidation zone [cm]
      real, parameter :: ph_min1 = 2.5               ! pH at which minimum oxidation rate is reached at acidic end of pH scale
      real, parameter :: ph_max1 = 4.0               ! pH at which maximum oxidation rate is reached at acidic end of pH scale
      real, parameter :: ph_min2 = 8.5               ! pH at which minimum oxidation rate is reached at acidic end of pH scale
      real, parameter :: ph_max2 = 7.0               ! pH at which maximum oxidation rate is reached at acidic end of pH scale
      real, parameter :: ph_rate_min = 0.2           ! Minimum value the pH oxidation rate modifier can take [proportion]

	! Arguments with intent(in)
	integer, intent(in) :: maxdepth    ! Maximum depth of the soil profile [cm]
	integer, intent(in) :: maxlayer    ! No.of layers in the soil profile
	integer, intent(in) :: measlay     ! Layer number that soil is measured to 
	integer, intent(in) :: maxsoil     ! Maximum number of soil types
	integer, intent(in) :: nsoil       ! Soil code number
	real, intent(in) :: alpha(maxlayer)      ! Proportion of BIO produced on BIO/DPM/RPM decomposition
      real, intent(in) :: asw(maxlayer)          ! Available soil water content [mm/layer]
      real, intent(in) :: asw_fc(maxlayer)       ! Available soil water content at field capacity [mm/layer]
      real, intent(in) :: asw_sat(maxlayer)      ! Available soil water content at saturation [mm/layer]
	real, intent(in) :: beta(maxlayer)       ! Proportion of HUM produced on BIO/DPM/RPM decomposition
      real, intent(in) :: biohum_n2c(maxsoil, maxlayer)  ! Stable N:C of the BIO & HUM pools
	real, intent(in) :: bulkdens(maxlayer)	 ! Bulk density of the layer [g/cm3]
	real, intent(in) :: co2(maxlayer)        ! CO2 produced in layer [kgC/ha]
      real, intent(in) :: co2_dpm(maxlayer)    ! CO2 produced from DPM pool in layer [kgC/ha]
      real, intent(in) :: co2_rpm(maxlayer)    ! CO2 produced from RPM pool in layer [kgC/ha]
      real, intent(in) :: co2_bio(maxlayer)    ! CO2 produced from BIO pool in layer [kgC/ha]
      real, intent(in) :: co2_hum(maxlayer)    ! CO2 produced from HUM pool in layer [kgC/ha]
	real, intent(in) :: crit(maxsoil,maxlayer) ! Minimum NO3-N and NH4-N level [kgN/ha/50cm]
	real, intent(in) :: delta(maxlayer)      ! Proportion of HUM produced on HUM decomposition
	real, intent(in) :: gamma(maxlayer)      ! Proportion of BIO produced on HUM decomposition
	real, intent(in) :: ph(maxsoil, maxlayer)  ! pH of each layer
      real, intent(in) :: seconds              ! Number of seconds in a timestep [seconds]
	real, intent(in) :: soilt(maxlayer) 	 ! Soil temperature [degC]
      real, intent(in) :: sw_sat(maxlayer)       ! Soil water content at saturation [mm/layer]
      real, intent(in) :: sw_wilt(maxlayer)      ! Soil water content at wilting point [mm/layer]

	! Arguments with intent(inout)
      real, intent(inout) :: bio_c(maxlayer)     ! C in BIO pool before decomposition [kgC/ha]
      real, intent(inout) :: bio_n(maxlayer)     ! N in BIO pool before decomposition [kgN/ha]
      real, intent(inout) :: dpm_c(maxlayer)     ! C in DPM pool before decomposition [kgC/ha]
      real, intent(inout) :: dpm_n(maxlayer)     ! N in DPM pool before decomposition [kgN/ha]
      real, intent(inout) :: hum_c(maxlayer)     ! C in HUM pool before decomposition [kgC/ha]
      real, intent(inout) :: hum_n(maxlayer)     ! N in HUM pool before decomposition [kgN/ha]
      real, intent(inout) :: no3_n(maxlayer)     ! N in the nitrate pool [kgN/ha]
      real, intent(inout) :: nh4_n(maxlayer)     ! N in the ammonium pool [kgN/ha/layer]
      real, intent(inout) :: rpm_c(maxlayer)     ! C in RPM pool before decomposition [kgC/ha]
      real, intent(inout) :: rpm_n(maxlayer)     ! N in RPM pool before decomposition [kgN/ha]

	! Arguments with intent(out)
	real, intent(out) :: ch4_atmos_ox 	     ! CH4 oxidised from atmospheric CH4 [kgC/ha]
	real, intent(out) :: ch4_flux	         ! Net CH4 flux from soil column (can be neg) [kgC/ha]
	real, intent(out) :: ch4_prod(maxlayer)  ! Gross CH4 produced in layer [kgC/ha/layer]
	real, intent(out) :: ch4_soil_ox 	     ! CH4 oxidised from CH4 produced in the soil [kgC/ha]

	! Local variables
	real :: bd_rate         ! Bulk density oxidation rate modifier [proportion]
	real :: ch4co2_ratio    ! CH4:CO2 ratio [proportion]
	real :: ch4_totprod     ! Total gross CH4 prudction of soil column [kgC/ha]
	integer :: i            ! Soil layer loop counter
	real :: layer_depth     ! Depth of a soil layer [cm]
	real :: max_atmos_ox    ! Maximum atmospheric methane oxidation rate [kgC/ha/timestep] (from Jang et al, 2006)
      real :: n_rate          ! Inorganic nitrogen concentration oxidation rate modifier [proportion]
	real :: ox_bd           ! Mean bulk density in the oxidation zone [g/cm3]
      real :: ox_mass         ! Mass of soil in the oxidation zone [kg/ha]
      real :: ox_n            ! Inorganic nitrogen in the oxidation zone [mgN/ha]
      real :: ox_n_conc       ! Inorganic N concentration in the oxidation zone [mgN/kg]
	real :: ox_ph           ! Mean pH in the oxidation zone
      real :: ox_wfps         ! Mean water-filled pore space of the oxidation zone [prop]
	real :: ox_t            ! Mean temperature of the oxidation zone [degC]
	real :: ph_rate         ! pH oxidation rate modifier [proportion]
      real :: sw(maxlayer)    ! Soil water content of each layer [mm/layer]
      real :: sw_fc(maxlayer) ! Soil water content at field capacity [mm/layer]
	real :: t_rate          ! Soil temperature oxidation rate modifier [proportion]
	real :: tot_oxdepth     ! Tracks the depth of oxidation zone. Used in case 
					        ! soil depth  is less than the ox_depth [cm]
	real :: w_rate          ! Soil water oxidation rate modifier [proportion]
	real :: z1, z2          ! Upper and lower depth of layer respectively [cm]

      ! Intialise variables
      max_atmos_ox = max_daily_atmos_ox * (seconds / daysecs)  ! Convert daily rate to timestep rate
      layer_depth = maxdepth / maxlayer

      ! Precalculate some useful soil water variables
      do i=1, measlay
          sw(i) = asw(i) + sw_wilt(i)
          sw_fc(i) = asw_fc(i) + sw_wilt(i)
      enddo
      
	!====================================================================
	! Methane Production 
      !====================================================================
 
      !          CH4:CO2 ratio
      !
      ! max|                       /
	!    |                      /
	!    |                     /
	!    |                    /
	!    |                   /
	!    |                  /
	!    |_________________/
	!   0|________________________
	! Wiltpoint          FC    SAT
	!        Soil Water Content

      do i=1, measlay
         ! If soil water is above field capacity...
          if (asw(i) > asw_fc(i)) then
              ch4co2_ratio = (ch4co2_ratio_sat - ch4co2_ratio_fc) * 
     &                (1 - (sw_sat(i) - sw(i)) / (sw_sat(i) - sw_fc(i)))
     &                + ch4co2_ratio_fc
          ! ...else soil water is less than or equal to field capacity
          else
              ch4co2_ratio = ch4co2_ratio_fc * 
     &                      (1 - (sw_fc(i) - sw(i)) / sw_fc(i))
          endif
          ch4_prod(i) = co2(i) * ch4co2_ratio
		if (ch4_prod(i) > 0.0) then
			call anaerobic_deocmposition(
     &                   ch4_prod(i), co2(i), co2_dpm(i), co2_rpm(i),
     &                   co2_bio(i), co2_hum(i), dpm_c(i), rpm_c(i),
     &                   bio_c(i), hum_c(i), dpm_n(i), rpm_n(i), 
     &                   bio_n(i), hum_n(i), nh4_n(i), no3_n(i),
     &                   biohum_n2c(nsoil,i), 
     &                   alpha(i), beta(i), delta(i), gamma(i), 
     &                   crit(nsoil,i), maxdepth, maxlayer)                   
		endif
      enddo  ! i (layers)	
      ch4_totprod = sum(ch4_prod)

      !====================================================================
	!  Methane Oxidation
      !====================================================================
	! Calculate soil properties/conditions within the oxidation zone
	ox_t = 0
	ox_wfps = 0
      ox_n = 0
	ox_ph = 0
	tot_oxdepth = 0
	z2 = 0
	do i=1, measlay
		z1 = z2
		z2 = z2 + layer_depth
		if (z2 <= ox_depth) then	! Layer is completely within oxidation zone
			! Weight T and WFPS of the layer by the layer's depth
			ox_t = ox_t + (soilt(i) * layer_depth)
              ox_wfps = ox_wfps + sw(i) / sw_sat(i) * layer_depth
			tot_oxdepth = tot_oxdepth + layer_depth
              ox_n = ox_n + no3_n(i) + nh4_n(i)
              ox_mass = ox_mass + (layer_depth * bulkdens(i) * 1e5)      ! Convert g/cm3 --> kg/ha/layer
			ox_ph = ox_ph + (ph(nsoil,i) * layer_depth)
			ox_bd = ox_bd + (bulkdens(i) * layer_depth)
		elseif (z2 > ox_depth .and. z1 < ox_depth) then
			! Layer overlaps bottom of oxidation zone so weight the T and WFPS
			! by the part of the layer that is within the oxidation zone
			ox_t = ox_t + (soilt(i) * (ox_depth - z1))
              ox_wfps = ox_wfps + sw(i) / sw_sat(i) * (ox_depth - z1)
              ox_n = ox_n + (((ox_depth - z1) / layer_depth) * 
     &               (no3_n(i) + nh4_n(i)))
              ox_mass = ox_mass + (bulkdens(i) * 1e5 * (ox_depth - z1))  ! Convert g/cm3 --> kg/ha/layer
			ox_ph = ox_ph + (ph(nsoil,i) * (ox_depth - z1))
			ox_bd = ox_bd + (bulkdens(i) * (ox_depth - z1))
			tot_oxdepth = tot_oxdepth + (ox_depth - z1)
			exit
		endif
      enddo  ! i (layers)
      ox_n = ox_n * 1e6     ! Convert kg/ha --> mg/ha
      ox_n_conc = ox_n / ox_mass
      ox_t = ox_t / tot_oxdepth
	ox_ph = ox_ph / tot_oxdepth
	ox_bd = ox_bd / tot_oxdepth
	ox_wfps = ox_wfps / tot_oxdepth

      !--------------------------------------------------------------------
      ! Soil temperature rate modifier of CH4 oxidation
      !--------------------------------------------------------------------
      ! Follows Ridgewell et al (1999) which is then normalised between 0-1
      ! Note this modifier is only applied to atmospheric methane oxidation
      ! because temperature already modifies CH4 production (via effects on
      ! CO2 production) and hence on the amount CH4 oxidation
      ! 1|              ____
      !  |             /    \
      !  |            /      \
      !  |          /         \
      !  |        /            \
      !  |      /               \
      !  |   /                   \
      ! 0|/_______________________\_
      !  0               29       43.25
      !               T (degC)

      t_rate = exp((0.0693 * ox_t) - (8.56e-7 * ox_t**4))
      t_rate = (t_rate - 1) / 3.1224  ! Normalise to between 0 and 1 for 0 <= T <= 43.25
      t_rate = min(1.0, t_rate)
      t_rate = max(0.0, t_rate)
	! Alternative Michaelis-Menton T curve that reaches asymptote at 25 deg C
	! t_rate = 1.12 * ox_t / (ox_t + 3.0)

      !--------------------------------------------------------------------
      ! Soil water rate modifier of CH4 oxidation
      !--------------------------------------------------------------------
      ! A Weibull distribution that approximates the soil water-oxidation
      ! relationships in Castro et al (1995), Whalen & Reeburgh (1996),
      ! Gulledge & Schimel (1998), van den Pol-van Dasselaar et al (1998)
      ! 1|       ____
      !  |      /     \
      !  |     /         \
      !  |    /             \
      !  |   /                 \
      !  |  /                    \
      !  | /                        \
      ! 0|/____________________________\_
      !  0        0.3                  1
      !           WFPS (proportion)

      w_rate = (2.172 / 0.3926) * ((ox_wfps / 0.3926)**1.172) *
     &         exp(-((ox_wfps / 0.3926)**2.172)) * 0.4324

      !--------------------------------------------------------------------
      ! Inorganic N concentration rate modifier of CH4 oxidation
      !--------------------------------------------------------------------
      ! Approximates the relationship given in Chan & Parkin (2001)

      !     1| |
      !      |  |
      !      |  |
      !      |   | 
      !      |   |
      !      |    \
      !      |      \
      !      |        \
      ! 0.05 |          \_________________
      !     0|____________________________
      !      0 1   5   10   15   20   25
      !         Inorganic N conc [mgN/kg]

      n_rate = min(1.12 * exp(-0.24 * ox_n_conc), 1.0)
      n_rate = max(n_rate, n_rate_min)
      
      !--------------------------------------------------------------------
      ! pH rate modifier of oxidation of CH4 oxidation
      !--------------------------------------------------------------------
      ! Note this modifier is only applied to atmospheric methane oxidation
      ! because pH already modifies CH4 production (via effects on CO2 production)
      ! and hence on the amount CH4 oxidation
      ! This response approximates that given in Saari et al (2004)
      !
      !   1  |        ____________         
      !      |       /:          :\
      !      |      / :          : \
      !      |     /  :          :  \
      !      |    /   :          :   \
      ! 0.2  | __/    :          :    \_____
      !      |________:__________:___________
      !      0  2.5   4        7.0   8.5  9.5
      !         min1 max1      max2  min2
      !                     pH
      if (ox_ph > ph_max1 .and. ox_ph < ph_max2) then
          ph_rate = 1.0
      elseif (ox_ph > ph_min1 .and. ox_ph <= ph_max1) then
          ph_rate = ph_rate_min + (1 - ph_rate_min) * 
     &              (ox_ph - ph_min1) / (ph_max1 - ph_min1)
          !ph_rate = 0.6333 * ox_ph - 1.5333
      elseif (ox_ph >= ph_max2 .and. ox_ph < ph_min2) then
          ph_rate = ph_rate_min + (1 - ph_rate_min) * 
     &              (ox_ph - ph_min2) / (ph_max2 - ph_min2)
              !ph_rate = -0.475 * ox_ph + 4.0875
      else  ! ph <= ph_min1 or ph >= ph_min2
         ph_rate = ph_rate_min
      endif
	
      !--------------------------------------------------------------------
      ! Bulk density rate modifier of CH4 oxidation
      !--------------------------------------------------------------------
      ! Based on Borken & Brumme (1997), Smith et al (2000)

      !      1| ________
      !       |         \
      !       |           \
      !       |             \
      !       |               \
      !       |                 \
      ! 0.05  |                   \_____
      !      0|________________________
      !       0        0.7        1.4 
      !              bd_max      bd_min
      !          Bulk density [g/cm3]

      if (ox_bd > bd_max .and. ox_bd < bd_min) then
          bd_rate = bd_rate_min + (1 - bd_rate_min) * 
     &              (ox_bd - bd_min) / (bd_max - bd_min)
          !bd_rate = -1.3571428571 * bulkdens(i) + 1.95
          !bd_rate = max(bd_rate_min, bd_rate)
          !bd_rate = min(1.0, bd_rate)
      elseif (ox_bd <= bd_max) then
          bd_rate = 1.0
      else  ! Bulk density >= bd_min
          bd_rate = bd_rate_min
      endif
	
	!--------------------------------------------------------------------
	! Calculate oxidation of soil and atmospheric methane
      !--------------------------------------------------------------------
	ch4_atmos_ox = max_atmos_ox * w_rate * t_rate * n_rate * ph_rate * 
     &               bd_rate
	ch4_soil_ox = ch4_totprod * max_ox_frac * w_rate * n_rate  ! Note: t_rate & ph_rate are already accounted for in aerobic decomposition
	
	!====================================================================
	! Calculate net methane flux
	!====================================================================
	ch4_flux = ch4_totprod - ch4_soil_ox - ch4_atmos_ox
	
	end subroutine methane_richards
!
!-------------------------------------------------------------------------------
!	
	subroutine anaerobic_deocmposition(
     &                 ch4_prod, co2, co2_dpm, co2_rpm, co2_bio,
     &                 co2_hum, dpm_c, rpm_c, bio_c, hum_c, dpm_n, 
     &                 rpm_n, bio_n, hum_n, nh4_n, no3_n,
     &                 biohum_n2c, alpha, beta, delta, gamma, 
     &                 crit, maxdepth, maxlayer)
      
      ! Calculates the anaerobic decomposition from the amount of CH4 produced
      ! and decreases CH4 production/decomposition under N-limiting conditions
      
	implicit none

    	! Arguments with intent(in)
	real, intent(in) :: alpha        ! Proportion of BIO produced on BIO/DPM/RPM decomposition
	real, intent(in) :: beta         ! Proportion of HUM produced on BIO/DPM/RPM decomposition
      real, intent(in) :: biohum_n2c  ! Stable N:C ratio of the BIO & HUM pools
      real, intent(in) :: co2          ! CO2 produced from all pools in layer [kgC/ha]
      real, intent(in) :: co2_dpm      ! CO2 produced from the DPM pool in layer [kgC/ha]
      real, intent(in) :: co2_rpm      ! CO2 produced from the RPM pool in layer [kgC/ha]
      real, intent(in) :: co2_bio      ! CO2 produced from the BIO pool in layer [kgC/ha]
      real, intent(in) :: co2_hum      ! CO2 produced from the HUM pool in layer [kgC/ha]
	real, intent(in) :: crit         ! Minimum NO3-N and NH4-N level [kgN/ha/50cm]
	real, intent(in) :: delta        ! Proportion of HUM produced on HUM decomposition
	real, intent(in) :: gamma        ! Proportion of BIO produced on HUM decomposition
      integer, intent(in) :: maxdepth  ! Maximum depth of the soil profile [cm]
      integer, intent(in) :: maxlayer  ! Number of layers in the soil profile

	! Arguments with intent(inout)
      real, intent(inout) :: bio_c     ! C in BIO pool before decomposition [kgC/ha]
      real, intent(inout) :: bio_n     ! N in BIO pool before decomposition [kgN/ha]
	real, intent(inout) :: ch4_prod  ! Gross CH4 produced [kgC/ha]
      real, intent(inout) :: dpm_c     ! C in DPM pool before decomposition [kgC/ha]
      real, intent(inout) :: dpm_n     ! N in DPM pool before decomposition [kgN/ha]
      real, intent(inout) :: hum_c     ! C in HUM pool before decomposition [kgC/ha]
      real, intent(inout) :: hum_n     ! N in HUM pool before decomposition [kgN/ha]
      real, intent(inout) :: no3_n     ! N in the nitrate pool [kgN/ha]
      real, intent(inout) :: nh4_n     ! N in the ammonium pool [kgN/ha/layer]
      real, intent(inout) :: rpm_c     ! C in RPM pool before decomposition [kgC/ha]
      real, intent(inout) :: rpm_n     ! N in RPM pool before decomposition [kgN/ha]
	
	! Local variables
      real :: alphabeta     ! alpha + beta
      real :: bio_c_prod    ! BIO C produced during anaerobic decomposition [kgC/ha]
	real :: bio_delta_c   ! Change in C in BIO pool after anaerobic decomposition [kgC/ha]
      real :: ch4_from_dpm  ! CH4-C produced from the DPM pool [kgC/ha]
      real :: ch4_from_rpm  ! CH4-C produced from the RPM pool [kgC/ha]
      real :: ch4_from_bio  ! CH4-C produced from the BIO pool [kgC/ha]
      real :: ch4_from_hum  ! CH4-C produced from the HUM pool [kgC/ha]
      real :: deltagamma    ! delta + gamma
      real :: dpm_delta_c   ! Change in C in DPM pool after anaerobic decomposition [kgC/ha]
      real :: dpm_n2c       ! N:C ratio of the DPM pool
      real :: hum_c_prod    ! HUM C produced during anaerobic decomposition [kgC/ha]
      real :: hum_delta_c   ! Change in C in HUM pool after anaerobic decomposition [kgC/ha]
	real :: minimum_n     ! Minimum N in the NO3 and NH4 pools [kgN/ha/layer]
	real :: modifier      ! Modifier used to reduce CH4 production under 
	                      ! N-limited conditions [proportion]
      real :: net_mineralised_n  ! Net mineralised N [kgN/ha/layer]
      real :: new_dpm_c     ! Size of DPM C pool after anaerobic decomposition [kgC/ha]
      real :: new_dpm_n     ! Size of DPM N pool after anaerobic decomposition [kgN/ha]
      real :: new_bio_c     ! Size of BIO C pool after anaerobic decomposition [kgC/ha]
      real :: new_bio_n     ! Size of BIO N pool after anaerobic decomposition [kgN/ha]
      real :: new_hum_c     ! Size of HUM C pool after anaerobic decomposition [kgC/ha]
      real :: new_hum_n     ! Size of HUM N pool after anaerobic decomposition [kgN/ha]
      real :: new_rpm_c     ! Size of RPM C pool after anaerobic decomposition [kgC/ha]
      real :: new_rpm_n     ! Size of RPM N pool after anaerobic decomposition [kgN/ha]
	real :: nh4_n_avail   ! N available for immobilisation from the
	                      ! ammonium pool [kg NH4-N/ha]
	real :: no3_n_avail   ! N available for immobilisation from the
	                      ! nitrate pool [kg NO3-N/ha]
	real :: orig_ch4_prod ! Gross CH4 produced in layer before modification due to N-limitation [kgC/ha]
      real :: rpm_n2c       ! N:C ratio of the RPM pool
      real :: rpm_delta_c   ! Change in C in RPM pool after anaerobic decomposition [kgC/ha]
	real :: tot_n_avail   ! Total amount of N available for immobilisation [kgN/ha]

      ! Local constants
      real, parameter :: mod_step = 0.025  ! Size of decrement to apply to modifier
	
      ! Initialise variables
      modifier = 1.0
      orig_ch4_prod = ch4_prod
      alphabeta = alpha + beta
      deltagamma = delta + gamma
      
      ! Calculate NO3-N and NH4-N available for immobilisation
	minimum_n = crit * maxdepth / (50.0 * maxlayer) ! Minimum amount of NO3-N and NH4-N in the layer
	no3_n_avail = max(no3_n - minimum_n, 0.0)
	nh4_n_avail = max(nh4_n - minimum_n, 0.0)
      tot_n_avail = no3_n_avail + nh4_n_avail
      
      ! Loop until the N needed by decomposition can be met by the 
      ! available N in the NO3 and NH4 pools
	do  
	    ch4_prod = orig_ch4_prod * modifier
		if (ch4_prod == 0.0) exit
		 
		! Calculate how much CH4-C comes from each pool
          ! Assumption: CH4-C comes from each pool in the same proportions
          ! as the CO2-C
          ch4_from_dpm = co2_dpm / co2 * ch4_prod
          ch4_from_rpm = co2_rpm / co2 * ch4_prod
	    ch4_from_bio = co2_bio / co2 * ch4_prod
          ch4_from_hum = co2_hum / co2 * ch4_prod
               
          ! Calculate the amount of C decomposed from each pool from the amount 
          ! of CH4 produced from each pool
          dpm_delta_c = ch4_from_dpm + ch4_from_dpm * 
     &                  (alphabeta) / (1 - alphabeta)
          rpm_delta_c  = ch4_from_rpm + ch4_from_rpm *
     &                   (alphabeta) / (1 - alphabeta)
          bio_delta_c  = ch4_from_bio + ch4_from_bio *
     &                   (alphabeta) / (1 - alphabeta)
          hum_delta_c  = ch4_from_hum + ch4_from_hum *
     &                   (deltagamma) / (1 - deltagamma)

          ! Calculate amount of BIO and HUM C produced during anaerobic decomposition
          bio_c_prod = alpha * (dpm_delta_c + rpm_delta_c + bio_delta_c)  ! DPM, RPM & BIO --> BIO
          bio_c_prod = bio_c_prod + (gamma * hum_delta_c)                 ! HUM --> BIO
          hum_c_prod = beta * (dpm_delta_c + rpm_delta_c + bio_delta_c)   ! DPM, RPM & BIO --> HUM
          hum_c_prod = hum_c_prod + (delta * hum_delta_c)                 ! HUM --> HUM

		! Calculate the new carbon pool sizes
          new_dpm_c = dpm_c - dpm_delta_c
          new_rpm_c = rpm_c - rpm_delta_c
          new_bio_c = bio_c - bio_delta_c + bio_c_prod  
          new_hum_c = hum_c - hum_delta_c + hum_c_prod
   
          ! Calculate N:C of the DPM & RPM pools
          if (dpm_c <= 0.0) then
              dpm_n2c = 0
          else
              dpm_n2c = dpm_n / dpm_c
          endif
          if (rpm_c <= 0.0) then
              rpm_n2c = 0
          else
              rpm_n2c = rpm_n / rpm_c
          endif
       
          ! Calculate new nitrogen pool sizes after mineralisation
          new_dpm_n = dpm_n - dpm_delta_c * dpm_n2c
          new_rpm_n = rpm_n - rpm_delta_c * rpm_n2c
          new_bio_n = new_bio_c * biohum_n2c
          new_hum_n = new_hum_c * biohum_n2c
          
          ! Calculate net mineralised N
          net_mineralised_n = (dpm_n - new_dpm_n) +
     &                        (rpm_n - new_rpm_n) +
     &                        (bio_n - new_bio_n) +
     &                        (hum_n - new_hum_n)

          ! If the net mineralised N is negative then some mineral N needs to be
          ! immobilised to makeup the shortfall.If there is not enough mineral N
          ! available to meet the requirement then CH4 production becomes N-limited
          if (-net_mineralised_n <= tot_n_avail) then
              exit
          ! ...else decomposition is N-limited so decrement the decomposition modifier
          else
              modifier = modifier - mod_step
          endif
      enddo  ! until total N needed < total N available
     
      if (ch4_prod > 0.0) then
          ! Calculate the new NO3 and NH4 pool sizes
          ! If there is excess N, add it to the ammonium pool...
          if (net_mineralised_n >= 0) then
              nh4_n = nh4_n + net_mineralised_n  ! Remember, net_mineralised_n is negative if there is excess N!
          !...or if pools need more N, take it from the ammonium pool...
          elseif (-net_mineralised_n < nh4_n_avail) then
              nh4_n = nh4_n + net_mineralised_n
          !...then when that runs out, take it from the nitrate pool as well...
          elseif (-net_mineralised_n < tot_n_avail) then
              nh4_n = minimum_n
              no3_n = no3_n + net_mineralised_n + nh4_n_avail
          else  ! Should never get here
              print *,"ERROR: tot_n_needed > tot_n_avail"
          endif      

          ! Calculate the new carbon pool sizes
          dpm_c = dpm_c - dpm_delta_c
          rpm_c = rpm_c - rpm_delta_c
          bio_c = bio_c - bio_delta_c + bio_c_prod  ! TODO: Check that need to add bio_c_prod, may already be accounted for in CO2 production??
          hum_c = hum_c - hum_delta_c + hum_c_prod  ! TODO: Check that need to add hum_c_prod, may already be accounted for in CO2 production??

          ! Calculate the new nitrogen pool sizes
          dpm_n = max(dpm_n - (dpm_delta_c * dpm_n2c), 0.0)
          rpm_n = max(rpm_n - (rpm_delta_c * rpm_n2c), 0.0)
          bio_n = bio_c * biohum_n2c
          hum_n = hum_c * biohum_n2c
      endif
     
      end subroutine  ! anaerobic_decomposition()

