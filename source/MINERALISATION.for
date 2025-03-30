!===============================================================================
! Name   : MINERALISATION.FOR
! Purpose: SUNDIAL & RothC mineralisation subroutines
! Author : Mark Richards, University of Aberdeen 
! Notes  : Adapted from the MINER1_BRADBURY() routine originally written by
!          Nick Bradbury
!
! Subroutines
! --------------------
! miner1_richards()         - High level/public mineralisation subroutine
! add_pi_to_dpmrpm()        - Adds plant inputs to the DPM & RPM pools
! hadley_water_rate_mods()  - Calculates the Hadley Centre water rate modifer of decomposition
! mineralisation()          - Aerobic mineralisation routine
! modfacts_miner2()         - Calculates the rate modifiers of decomposition
! rothc_water_rate_mods()   - Calculates the RotchC water rate modifer of decomposition
!===============================================================================

      subroutine miner1_richards(maxlayer, maxdepth, measlay,
     &                 nsoil, maxsoil, cover,
     &                 dpm_rate, rpm_rate, bio_rate, hum_rate,
     &                 alpha, beta, delta, gamma, biohum_n2c,
     &                 t_func, w_func, ic_factor,
     &                 soilw, wmax, wsat, soilt, min_n_lvl,
     &                 wiltpoint, satwatcont, ph, ph_min, ph_max,
     &                 no3_n, nh4_n,
     &                 dpm_c, dpm_n, rpm_c, rpm_n,
     &                 bio_c, bio_n, hum_c, hum_n,
     &                 co2, co2_dpm, co2_rpm, co2_bio, co2_hum,
     &                 net_mineralised_n, tot_mineralised_n,
     &                 w_ratedm, t_ratem,PI_C)

      implicit none

      ! Scalars with intent (in)
      integer, intent(in) :: cover              ! Code for crop cover: 1=covered; 2=bare
      integer, intent(in) :: maxdepth           ! Maximum depth of the soil profile [cm]
      integer, intent(in) :: maxlayer           ! Number of layers in the soil profile
      integer, intent(in) :: maxsoil            ! Max no. of soil types
      integer, intent(in) :: measlay            ! Layer to which soil is measured
      real, intent(in) :: min_n_lvl             ! Minimum amount of N in the soil [kgN/ha/50cm]
      integer, intent(in) :: nsoil              ! Soil code number
      real, intent(in) :: ph(maxsoil,maxlayer)       ! Soil pH in layer
      real, intent(in) :: ph_min(maxsoil,maxlayer)   ! pH below which rate of decomposition is zero
      real, intent(in) :: ph_max (maxsoil,maxlayer)  ! pH above which rate of decomposition is maximum
      real, intent(in) :: biohum_n2c(maxsoil, maxlayer)  ! Stable N:C ratio of the BIO a& HUM pools [ratio]
      integer, intent(in) :: t_func             ! Choice of temperature rate modifier 1=RothC; 2=Hadley
      integer, intent(in) :: w_func             ! Choice of water rate modifier 1=RothC; 2=Hadley
      
      ! Arrays with intent (in)
      real, intent(in) :: alpha(maxlayer)       ! Proportion of BIO produced on BIO decomposition
      real, intent(in) :: beta(maxlayer)        ! Proportion of HUM produced on BIO decomposition
      real, intent(in) :: bio_rate(maxlayer)    ! Rate constant for decomposition of the BIO pool
      real, intent(in) :: delta(maxlayer)       ! Proportion of HUM produced on HUM decomposition
      real, intent(in) :: dpm_rate(maxlayer)    ! Rate constant for decomposition of the DPM pool
      real, intent(in) :: gamma(maxlayer)       ! Proportion of BIO produced on HUM decomposition
      real, intent(in) :: hum_rate(maxlayer)    ! Rate constant for decomposition of the HUM pool
      real, intent(in) :: ic_factor(maxlayer)   ! Adjustment in decomposition rate needed to achieve measured NPP and TOC
      real, intent(in) :: rpm_rate(maxlayer)    ! Rate constant for decomposition of the RPM pool
      real, intent(in) :: satwatcont(maxlayer)  ! Total water content at saturation [mm/layer]
      real, intent(in) :: soilt(maxlayer)       ! Soil temperature [degC]
      real, intent(in) :: soilw(maxlayer)       ! Available soil water content [mm/layer]
      real, intent(in) :: wiltpoint(maxlayer)   ! Water content at wilting point [mm/layer]
      real, intent(in) :: wmax(maxlayer)        ! Available soil water at field capacity [mm/layer]
      real, intent(in) :: wsat(maxlayer)        ! Available soil water at saturation [mm/layer]

      ! Arrays with intent (inout)
      real, intent(inout) :: bio_c(maxlayer)    ! C in Biomass pool before decomposition [kgC/ha/layer]
      real, intent(inout) :: bio_n(maxlayer)    ! N in Biomass pool before decomposition [kgN/ha/layer]
      real, intent(inout) :: dpm_c(maxlayer)    ! C in Decomposable Plant Matter pool before decomposition [kgC/ha/layer]
      real, intent(inout) :: dpm_n(maxlayer)    ! N in Decomposable Plant Matter pool before decomposition [kgN/ha/layer]
      real, intent(inout) :: hum_c(maxlayer)    ! C in Humus pool before decomposition [kgC/ha/layer]
      real, intent(inout) :: hum_n(maxlayer)    ! N in Humus pool before decomposition [kgN/ha/layer]
      real, intent(inout) :: nh4_n(maxlayer)    ! N in the ammonium pool [kgN/ha/layer]
      real, intent(inout) :: no3_n(maxlayer)    ! N in the nitrate pool [kgN/ha/layer]
      real, intent(inout) :: rpm_c(maxlayer)    ! C in Resistant Plant Matter pool before decomposition [kgC/ha/layer]
      real, intent(inout) :: rpm_n(maxlayer),PI_C(maxlayer)    ! N in Resistant Plant Matter pool before decomposition [kgN/ha/layer]
      real, intent(inout) :: net_mineralised_n(maxlayer)   ! Net mineralised N (kgN/ha/layer/timestep)

      ! Scalars with intent (out)
      real, intent(out) :: tot_mineralised_n    ! Total net mineralised N [kgN/ha]

      ! Arrays with intent (out)
      real, intent(out) :: co2(maxlayer)     ! CO2 produced during aerobic decomposition [kgCO2-C/ha/layer]
      real, intent(out) :: co2_dpm(maxlayer) ! CO2 emitted from the DPM pool [kgCO2-C/ha/layer]
      real, intent(out) :: co2_rpm(maxlayer) ! CO2 emitted from the RPM pool [kgCO2-C/ha/layer]
      real, intent(out) :: co2_bio(maxlayer) ! CO2 emitted from the BIO pool [kgCO2-C/ha/layer]
      real, intent(out) :: co2_hum(maxlayer) ! CO2 emitted from the HUM pool [kgCO2-C/ha/layer]
      
      ! Local scalar variables
      real :: adj_factor           ! Value used to adjust the eff. of decompisition in order to maintain the stable C:N ratio of BIO & HUM [proportion]
      real :: crop_rate            ! Crop rate modifier for decomposition [proportion]
      real*8 :: exp_bio            ! Exponential factor for decomposition of the BIO pool
      real*8 :: exp_dpm            ! Exponential factor for decomposition of the DPM pool
      real*8 :: exp_hum            ! Exponential factor for decomposition of the HUM pool
      real*8 :: exp_rpm            ! Exponential factor for decomposition of the RPM pool
      integer :: i                 ! Soil layer loop counter
      real :: layer_depth          ! Depth (thickness) of soil layer [cm]
      real :: minimum_n            ! Minimum level of N in a layer [kgN/ha/layer]
      real :: mod_alpha            ! Proportion of decomposed BIO that goes to BIO
      real :: mod_beta             ! Proportion of decomposed BIO that goes to HUM
      real :: mod_delta            ! Proportion of decomposed HUM that goes to HUM
      real :: mod_gamma            ! Proportion of decomposed HUM that goes to BIO
      real :: new_bio_c            ! C in Biomass pool after decomposition [kgC/ha/layer]
      real :: new_bio_n            ! N in Biomass pool after decomposition [kgN/ha/layer]
      real :: new_dpm_c            ! C in Decomposable Plant Matter pool after decomposition [kgC/ha/layer]
      real :: new_dpm_n            ! N in Decomposable Plant Matter pool after decomposition [kgN/ha/layer]
      real :: new_hum_c            ! C in Humus pool after decomposition [kgC/ha/layer]
      real :: new_hum_n            ! N in Humus pool after decomposition [kgN/ha/layer]
      real :: new_rpm_c            ! C in Resistant Plant Matter pool after decomposition [kgC/ha/layer]
      real :: new_rpm_n            ! N in Resistant Plant Matter pool after decomposition [kgN/ha/layer]
      real :: nh4_n_avail4immob    ! Amount of ammonium-N available for immobilisation [kgN/ha]
      real :: no3_n_avail4immob    ! Amount of nitrate-N available for immobilisation [kgN/ha]
      real :: ph_rate              ! pH rate modifier for decomposition [proportion]
      real*8 :: rate_mod           ! Product of the rate modifiers for decomposition [proportion]
      real :: t_rate               ! Water rate modifier for decomposition [proportion]
      real :: t_ratem
      real :: tot_n_avail4immob    ! Total N available for immobilisation [kgN/ha/layer]
      real :: w_ratedm
      real :: w_rate_dpmbio        ! Water rate modifier for DPM & BIO pools [proportion]
      real :: w_rate_rpmhum        ! Water rate modifier for RPM & HUM pools [proportion]
      
      ! Local logical variables
      logical :: adj_factor_found ! Flag to indicate whether the search for the appropriate eff. of decomp. value has finished
      
      ! Initialise variables
      w_ratedm = 0
      t_ratem = 0
      layer_depth = maxdepth / maxlayer
      minimum_n = min_n_lvl * (layer_depth / 50.0)

      ! Crop modification factor
      if (cover == 1) then
          crop_rate = 0.6
      else
          crop_rate = 1.0
      endif
     
      do i=1, measlay
          ! Calculate the water, temperature & pH rate modifiying factors for decomposition
          call modfacts_miner2(layer_depth, soilw(i), wmax(i),
     &             wsat(i), soilt(i), wiltpoint(i), satwatcont(i),
     &    ph(nsoil,i), ph_min(nsoil,i), ph_max(nsoil,i), t_func, w_func,
     &             w_rate_dpmbio,w_rate_rpmhum, t_rate, ph_rate)
          
      
          ! Sum the water and temperature rate modifiers over the first 10 layers
          ! (needed for Monte Carlo?)
          if (i < 11) then
              w_ratedm = w_ratedm + w_rate_dpmbio
              t_ratem = t_ratem + t_rate
          endif
                      
          ! Calculate each pool's exponential factor for decomposition
          rate_mod = w_rate_dpmbio * t_rate * crop_rate * ph_rate *
     &               ic_factor(i)
          exp_dpm = exp(-dpm_rate(i) * rate_mod)
          exp_rpm = exp(-rpm_rate(i) * rate_mod)
          exp_bio = exp(-bio_rate(i) * rate_mod)
          exp_hum = exp(-hum_rate(i) * rate_mod)
     
          mod_alpha = alpha(i)
          mod_beta = beta(i)
          mod_delta = delta(i)
          mod_gamma = gamma(i)
          
          call mineralisation(biohum_n2c(nsoil,i),
     &                    mod_alpha, mod_beta, mod_delta, mod_gamma, 
     &                    exp_dpm, exp_rpm, exp_bio, exp_hum,
     &                    dpm_c(i), dpm_n(i), rpm_c(i), rpm_n(i),
     &                    bio_c(i), bio_n(i), hum_c(i), hum_n(i),
     &                    new_dpm_c, new_dpm_n, new_rpm_c, new_rpm_n,
     &                    new_bio_c, new_bio_n, new_hum_c, new_hum_n,
     &                    net_mineralised_n(i),
     &                    minimum_n, no3_n(i), nh4_n(i),
     &                    no3_n_avail4immob, nh4_n_avail4immob, 
     &                    tot_n_avail4immob)

          ! If the N demand (i.e. negative net mineralised N) exceeds the N available for
          ! immobilisation then adjust the efficiency of decomposition (fraction of C released 
          ! from DPM and RPM retained in the soil) until enough N is available. This maintains 
          ! a stable BIO & HUM N:C ratio.
          if (-net_mineralised_n(i) > tot_n_avail4immob) then
              adj_factor = 0.9
              adj_factor_found = .false.
              do while (.not. adj_factor_found)
                  ! Adjust the efficiency of decomposition (alpha, beta etc.) 
                  ! by the adjustmnent factor
                  mod_alpha = alpha(i) * adj_factor
                  mod_beta = beta(i) * adj_factor
                  mod_delta = delta(i) * adj_factor
                  mod_gamma = gamma(i) * adj_factor

                  if (adj_factor < 0.01 .or. 
     &                          tot_n_avail4immob < 0.0001) then
                      ! Switch off decomposition
                      exp_dpm = 1
                      exp_rpm = 1
                      exp_bio = 1
                      exp_hum = 1
                      adj_factor_found = .true.
                  endif

                  call mineralisation(biohum_n2c(nsoil,i),
     &                    mod_alpha, mod_beta, mod_delta, mod_gamma, 
     &                    exp_dpm, exp_rpm, exp_bio, exp_hum,
     &                    dpm_c(i), dpm_n(i), rpm_c(i), rpm_n(i),
     &                    bio_c(i), bio_n(i), hum_c(i), hum_n(i),
     &                    new_dpm_c, new_dpm_n, new_rpm_c, new_rpm_n,
     &                    new_bio_c, new_bio_n, new_hum_c, new_hum_n,
     &                    net_mineralised_n(i),
     &                    minimum_n, no3_n(i), nh4_n(i), 
     &                    no3_n_avail4immob, nh4_n_avail4immob,
     &                    tot_n_avail4immob)
                  
                  ! If there's still not enough N available then decrease the 
                  ! efficiency of decomposition further and try again
                  if (-net_mineralised_n(i) > tot_n_avail4immob) then 
                      adj_factor = adj_factor * 0.9
                  else 
                      adj_factor_found = .true.
                  endif    
              enddo  ! while .not. adj_factor_found
                                
          ! If the ammount of N immobilised exceeds the N available as ammonium take the
          ! remaining mineralised N from the nitrate pool (note: net_mineralised_n(i) is negative)
          elseif (-net_mineralised_n(i) > nh4_n_avail4immob) then
              no3_n(i) = no3_n(i) + nh4_n_avail4immob + 
     &                   net_mineralised_n(i)
              nh4_n(i) = minimum_n
          ! If the amount of N immobilised does not exceed N available as ammonium take
          ! all the N from the ammonium pool (note: net_mineralised_n(i) is negative)
          else
              nh4_n(i) = nh4_n(i) + net_mineralised_n(i)
          endif
          
          ! Calculate the CO2 released due to biological activity
          co2_dpm(i) = 
     &         (1 - (mod_alpha + mod_beta)) * dpm_c(i) * (1 - exp_dpm)
          co2_rpm(i) = 
     &         (1 - (mod_alpha + mod_beta)) * rpm_c(i) * (1 - exp_rpm)
          co2_bio(i) =
     &         (1 - (mod_alpha + mod_beta)) * bio_c(i) * (1 - exp_bio)
          co2_hum(i) =
     &         (1 - (mod_gamma + mod_delta)) * hum_c(i) * (1 - exp_hum)
          co2(i) = co2_dpm(i) + co2_rpm(i) + co2_bio(i) + co2_hum(i)

         ! Set the pools to the size of the new pools (i.e. after decomposition)
          dpm_c(i) = new_dpm_c
          dpm_n(i) = new_dpm_n
          rpm_c(i) = new_rpm_c
          rpm_n(i) = new_rpm_n
          bio_c(i) = new_bio_c
          bio_n(i) = new_bio_n
          hum_c(i) = new_hum_c
          hum_n(i) = new_hum_n

          tot_mineralised_n = tot_mineralised_n + net_mineralised_n(i)
     &                          
      enddo  ! i (layers) 
 
      end subroutine miner1_richards

!---------------------------------------------------------------------------------------------
    
      subroutine mineralisation(stable_biohum_n2c,
     &                          alpha, beta, delta, gamma, 
     &                          exp_dpm, exp_rpm, exp_bio, exp_hum,
     &                          dpm_c, dpm_n, rpm_c, rpm_n,
     &                          bio_c, bio_n, hum_c, hum_n,
     &                          new_dpm_c, new_dpm_n,
     &                          new_rpm_c, new_rpm_n,
     &                          new_bio_c, new_bio_n,
     &                          new_hum_c, new_hum_n,
     &                          net_mineralised_n,
     &                          minimum_n, no3_n, nh4_n,
     &                          no3_n_avail4immob,
     &                          nh4_n_avail4immob, tot_n_avail4immob)     
      implicit none
      
      ! Scalar arguments with intent(in)
      real, intent(in) :: alpha       ! Proportion of BIO produced on BIO decomposition
      real, intent(in) :: beta        ! Proportion of HUM produced on BIO decomposition
      real, intent(in) :: bio_c       ! C in Biomass pool before decomposition [kgC/ha/layer]
      real, intent(in) :: bio_n       ! N in Biomass pool before decomposition [kgN/ha/layer]
      real, intent(in) :: delta       ! Proportion of HUM produced on HUM decomposition
      real, intent(in) :: dpm_c       ! C in Decomposable Plant Matter pool before decomposition [kgC/ha/layer]
      real, intent(in) :: dpm_n       ! N in Decomposable Plant Matter pool before decomposition [kgN/ha/layer] 
      real*8, intent(in) :: exp_bio   ! Exponential factor for decomposition of the BIO pool
      real*8, intent(in) :: exp_dpm   ! Exponential factor for decomposition of the DPM pool
      real*8, intent(in) :: exp_hum   ! Exponential factor for decomposition of the HUM pool
      real*8, intent(in) :: exp_rpm   ! Exponential factor for decomposition of the RPM pool
      real, intent(in) :: gamma       ! Proportion of BIO produced on HUM decomposition
      real, intent(in) :: hum_c       ! C in Humus pool before decomposition [kgC/ha/layer]
      real, intent(in) :: hum_n       ! N in Humus pool before decomposition [kgN/ha/layer]
      real, intent(in) :: minimum_n   ! Minimum level of N in layer [kgN/ha/layer]
      real, intent(in) :: no3_n       ! N in the nitrate pool [kgN/ha/layer]
      real, intent(in) :: nh4_n       ! N in the ammonium pool [kgN/ha/layer]
      real, intent(in) :: rpm_c       ! C in Resistant Plant Matter pool before decomposition [kgC/ha/layer]
      real, intent(in) :: rpm_n       ! N in Resistant Plant Matter pool before decomposition [kgN/ha/layer]
      real, intent(in) :: stable_biohum_n2c  ! Stable N:C ratio of the BIO a& HUM pools [ratio]

      ! Scalar arguments with intent(out)
      real, intent(out) :: net_mineralised_n  ! Net mineralised N [kgN/ha/layer/timestep]
      real, intent(out) :: new_bio_c          ! C in Biomass pool after decomposition [kgC/ha/layer]
      real, intent(out) :: new_bio_n          ! N in Biomass pool after decomposition [kgN/ha/layer]
      real, intent(out) :: new_dpm_c          ! C in Decomposable Plant Matter pool after decomposition [kgC/ha/layer]
      real, intent(out) :: new_dpm_n          ! N in Decomposable Plant Matter pool after decomposition [kgN/ha/layer]
      real, intent(out) :: new_hum_c          ! C in Humus pool after decomposition [kgC/ha/layer]
      real, intent(out) :: new_hum_n          ! N in Humus pool after decomposition [kgN/ha/layer]
      real, intent(out) :: new_rpm_c          ! C in Resistant Plant Matter pool after decomposition [kgC/ha/layer]
      real, intent(out) :: new_rpm_n          ! N in Resistant Plant Matter pool after decomposition [kgN/ha/layer]
      real, intent(out) :: nh4_n_avail4immob  ! Amount of ammonium-N available for immobilisation [kgN/ha]
      real, intent(out) :: no3_n_avail4immob  ! Amount of nitrate-N available for immobilisation [kgN/ha]
      real, intent(out) :: tot_n_avail4immob  ! Total N available for immobilisation [kgN/ha/layer]
     
      !-------------------------------------------------------------------------------------
      ! Calculate the new carbon pool size after mineralisation
      !-------------------------------------------------------------------------------------
      new_dpm_c = dpm_c * exp_dpm
      new_rpm_c = rpm_c * exp_rpm
      new_bio_c = (bio_c * exp_bio) +
     &            (alpha * bio_c * (1 - exp_bio)) +  ! BIO to BIO
     &            (gamma * hum_c * (1 - exp_hum)) +  ! HUM to BIO
     &            (alpha * dpm_c * (1 - exp_dpm)) +  ! DPM to BIO
     &            (alpha * rpm_c * (1 - exp_rpm))    ! RPM to BIO     
      new_hum_c = (hum_c * exp_hum) +
     &            (delta * hum_c * (1 - exp_hum)) +  ! HUM to BIO
     &            (beta * bio_c * (1 - exp_bio)) +   ! BIO to HUM
     &            (beta * dpm_c * (1 - exp_dpm)) +   ! DPM to HUM
     &            (beta * rpm_c * (1 - exp_rpm))     ! RPM to HUM

      !-------------------------------------------------------------------------------------
      ! Calculate the new nitrogen pool sizes after mineralisation
      !-------------------------------------------------------------------------------------     
      new_dpm_n = dpm_n * exp_dpm
      new_rpm_n = rpm_n * exp_rpm
      new_bio_n = new_bio_c * stable_biohum_n2c
      new_hum_n = new_hum_c * stable_biohum_n2c

      !-------------------------------------------------------------------------------------
      ! Calculate the net mineralised N
      !-------------------------------------------------------------------------------------
      net_mineralised_n = (dpm_n - new_dpm_n) +
     &                    (rpm_n - new_rpm_n) +
     &                    (bio_n - new_bio_n) +
     &                    (hum_n - new_hum_n)

      !-------------------------------------------------------------------------------------
      ! Calculate the amount of nitrate-N and ammomium-N available for immobilisation
      !-------------------------------------------------------------------------------------
      nh4_n_avail4immob = max(nh4_n - minimum_n, 0.0)
      no3_n_avail4immob = max(no3_n - minimum_n, 0.0)
      tot_n_avail4immob = nh4_n_avail4immob + no3_n_avail4immob

      end subroutine mineralisation

!----------------------------------------------------------------------------------------------

      subroutine modfacts_miner2(layer_depth,
     &                           soilw, wmax, wsat, soilt, 
     &                           wiltpoint, satwatcont,
     &                           ph, ph_min, ph_max, t_func, w_func,
     &                           w_rate_dpmbio, w_rate_rpmhum, t_rate,
     &                           ph_rate)
      
      ! Calculates the rate modifiying factors for mineralisation
   
      implicit none

      ! Scalar arguments with intent(in)
      integer, intent(in) :: t_func   ! Choice of temperature rate modifier
      integer, intent(in) :: w_func   ! Choice of water rate modifier
      real, intent(in) :: layer_depth ! Depth (thickness) of soil layer [cm]
      real, intent(in) :: ph          ! Soil pH in layer
      real, intent(in) :: ph_min      ! pH below which rate of decomposition is zero
      real, intent(in) :: ph_max      ! pH above which rate of decomposition is maximum
      real, intent(in) :: satwatcont  ! Total water content at saturation [mm/layer]
      real, intent(in) :: soilt       ! Soil temperature [degC]
      real, intent(in) :: soilw       ! Available soil water content [mm/layer]
      real, intent(in) :: wiltpoint   ! Soil water content at wilting point [mm/layer]
      real, intent(in) :: wmax        ! Available soil water at field capacity [mm/layer]
      real, intent(in) :: wsat        ! Available soil water at saturation [mm/layer]

      ! Scalar arguments with intent(out)
      real, intent(out) :: ph_rate         ! pH rate modifier [proportion]
      real, intent(out) :: t_rate          ! Temperature rate modifying factor [prop.]
      real, intent(out) :: w_rate_dpmbio   ! Water rate modifier for DPM and BIO
      real, intent(out) :: w_rate_rpmhum   ! Water rate modifier for RPM and HUM

      ! Local scalar variables
      real :: soilt_kelvin     ! Soil temperature [deg Kelvin]
      real :: swc_at_onebar    ! Soil water content per layer at -1 bar (or 20mm/25cm)

      ! Local constants
      integer, parameter :: rothc_wmod = 0   ! RothC soil water rate modifier function
      integer, parameter :: hadley_wmod = 1  ! Hadley soil water rate modifier function
      integer, parameter :: rothc_tmod = 0   ! RothC temperature rate modifier function
      integer, parameter :: hadley_tmod = 1  ! Hadley temperature rate modifier function
      real, parameter :: q10 = 2.0           ! Q10 value used in Hadley temperature rate modifier calculation
      real, parameter :: ph_rate_min = 0.2   ! Minimum value of pH rate modifier [proportion]
      
      swc_at_onebar = 20.0 / 25.0 * layer_depth  ! SWC per layer at -1 bar, i.e. 20mm/25cm

      !-----------------------------------------------------------------------------------
      ! Soil water rate modifiers
      !-----------------------------------------------------------------------------------
      if (w_func == rothc_wmod) then
          call rothc_water_rate_mods(layer_depth, soilw, wmax, wsat,
     &                               w_rate_dpmbio, w_rate_rpmhum)
      elseif (w_func == hadley_wmod) then
          call hadley_water_rate_mods(soilw, wmax, wsat, wiltpoint,
     &                                satwatcont, w_rate_dpmbio,
     &                                w_rate_rpmhum)
      endif
      
      !-----------------------------------------------------------------------------------
      ! Soil temperature rate modifier
      !-----------------------------------------------------------------------------------
      if (t_func == rothc_tmod) then
          if (soilt < -10) then
              t_rate = 0.0
          else
              t_rate = 47.91 / (1 + exp(106.06 / (soilt + 18.27)))
          endif
      elseif (t_func == hadley_tmod) then
          soilt_kelvin = soilt + 273.15
          t_rate = q10**(0.1 * (soilt_kelvin - 298.15))
      endif

      !-----------------------------------------------------------------------------------
      ! pH rate modifier
      !-----------------------------------------------------------------------------------
      if (ph <= ph_min) then
          ph_rate = ph_rate_min
      elseif (ph > ph_min .and. ph < ph_max) then
          ph_rate = ph_rate_min + (1 - ph_rate_min) * 
     &              (ph - ph_min) / (ph_max - ph_min)
      elseif (ph >= ph_max) then
          ph_rate = 1.0
      endif
      
      end subroutine modfacts_miner2

!---------------------------------------------------------------------------------------------

      subroutine rothc_water_rate_mods(layer_depth, soilw, wmax, wsat,
     &                                 w_rate_dpmbio, w_rate_rpmhum)
      
      ! Subroutine to calculate the RothC water rate modifiying factors for mineralisation
   
      implicit none

      ! Scalar arguments with intent(in)
      real, intent(in) :: layer_depth  ! Thickness of the soil layer [cm]
      real, intent(in) :: soilw    ! Available soil water content [mm/layer]
      real, intent(in) :: wmax     ! Available soil water at field capacity [mm/layer]
      real, intent(in) :: wsat     ! Available soil water at saturation [mm/layer]

      ! Scalar arguments with intent(out)
      real, intent(out) :: w_rate_dpmbio   ! Water rate modifier for DPM and BIO
      real, intent(out) :: w_rate_rpmhum   ! Water rate modifier for RPM and HUM

      ! Local scalar variables
      real :: swd_below_fc   ! Soil water deficit below field capacity [mm/layer]
      real :: swc_at_onebar  ! Soil water content at -1 bar, i.e. 20mm/25cm [mm/layer]

      ! Local constants
      real, parameter :: wmin = 0.2   ! Min soil water content as a fraction of swc at field capacity

      swc_at_onebar = 20.0 / 25.0 * layer_depth  ! SWC per layer at -1 bar, i.e. 20mm/25cm
      swd_below_fc = wmax - soilw
      
      if (swd_below_fc <= 0) then  
          ! Soil water content of layer is above field capacity
          w_rate_dpmbio = ((wmin - 1) * soilw) + wsat - (wmin * wmax)
          w_rate_dpmbio = w_rate_dpmbio / (wsat - wmax)
          w_rate_rpmhum = w_rate_dpmbio
      elseif (swd_below_fc <= swc_at_onebar .and. swd_below_fc > 0) then
          ! Soil water content of layer is below field capacity but above swc at -1 bar
          w_rate_dpmbio = 1.0
          w_rate_rpmhum = 1.0
      elseif (swd_below_fc > swc_at_onebar) then
          ! Soil water content of layer is less than swc at -1 bar
          w_rate_dpmbio = 1.0 - ((1-wmin) * 
     &                    (swd_below_fc - swc_at_onebar)) /
     &                    (wmax - swc_at_onebar)
          w_rate_rpmhum = w_rate_dpmbio          
      endif
      IF(W_RATE_DPMbio.LE.0.2)then
       w_rate_dpmbio=0.2
       w_rate_rpmhum = w_rate_dpmbio 
       endif
      end subroutine rothc_water_rate_mods

!---------------------------------------------------------------------------------------------

      subroutine hadley_water_rate_mods(soilw, wmax, wsat,
     &                                  wiltpoint, satwatcont,
     &                                  w_rate_dpmbio, w_rate_rpmhum)
      
      ! Subroutine to calculate the Hadley Centre water rate modifiying factors for mineralisation
   
      implicit none

      ! Scalar arguments with intent(in)
      real, intent(in) :: soilw           ! Available soil water content [mm/layer]
      real, intent(in) :: wmax            ! Available soil water at field capacity [mm/layer]
      real, intent(in) :: wsat            ! Available soil water at saturation [mm/layer]
      real, intent(in) :: wiltpoint       ! Soil water content at wilting point [mm/layer]
      real, intent(in) :: satwatcont      ! Soil water content at saturation [mm/layer]

      ! Scalar arguments with intent(out)
      real, intent(out) :: w_rate_dpmbio   ! Water rate modifier for DPM and BIO [proportion]
      real, intent(out) :: w_rate_rpmhum   ! Water rate modifier for RPM and HUM [proportion]

      ! Local scalar variables
      real :: swd_below_fc    ! Soil water deficit below field capacity [mm/layer]
      real :: wfps            ! Water filled pore space [porportion]
      real :: wfps_wiltpoint  ! Water filled pore space at wilting point [porportion]
      real :: wfps_opt        ! Water filled pore space at which respiration is maxmimum [proportion]

      wfps_wiltpoint = wiltpoint / satwatcont
      wfps = (soilw + wiltpoint) / satwatcont
      wfps_opt = 0.5 * (1 + wfps_wiltpoint)
      
      if (wfps <= wfps_wiltpoint) then
          w_rate_dpmbio = 0.2
          w_rate_rpmhum = 0.2
      elseif (wfps > wfps_wiltpoint .and. wfps <= wfps_opt) then
          w_rate_dpmbio = 0.2 + 0.8 * 
     &                    ((wfps - wfps_wiltpoint) /
     &                    (wfps_opt - wfps_wiltpoint))
          w_rate_rpmhum= w_rate_dpmbio  
      elseif (wfps > wfps_opt) then
          w_rate_dpmbio = 1 - 0.8 * (wfps - wfps_opt)
          w_rate_rpmhum = w_rate_dpmbio
      endif

      end subroutine hadley_water_rate_mods
