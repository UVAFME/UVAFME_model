module Model

!*******************************************************************************
  !
  ! This contains the five main subroutines of UVAFME:
  !
  ! BioGeoClimate: computes daily and yearly site- and plot-level weather and
  !                soil dynamics
  !
  ! Canopy:        computes the plot-level LAI and light availability
  !
  ! Growth:        computes annual tree growth and branch thinning
  !
  ! Mortality:     determines which trees die and adds their components to the
  !                soil
  !
  ! Renewal        updates the seed and seedling banks for each species and
  !                regenerates new trees
  !
!*******************************************************************************


    use Parameters
    use Constants
    use Soil
    use Site
    use Species
    use Tree
    use Random
    use Climate
    use Input

    implicit none


contains

    !:.........................................................................:

    subroutine BioGeoClimate(site, year)
        !
        !  Computes daily weather data and annaul sums of weather and soil
        !  characteristics
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong          Original Code
        !    01/01/12     K. Holcomb           Updated to OOP structure
        !    10/10/16     A. C. Foster         Updated for soil/plot overhaul
        !                                       and permafrost updates
        !    05/01/17     A. C. Foster         Updated for moss and nutrient
        !                                       updates
        !

        ! Data dictionary: constants
        real,    parameter    :: MIN_GROW_TEMP = 5.0   ! Minimum temperature for growing season (degC)
        real,    parameter    :: MAX_DRY_PARM = 1.0001 ! Threshold for below wilting point
        real,    parameter    :: DRY_THRESH = 0.8      ! Threshold for droughty conditions
        real,    parameter    :: MIN_FLOOD_PARM = 0.8  ! Threshold for flooded conditions
        integer, parameter    :: N_MEM = 3             ! Number of days to calculate temperature "memory"

        ! Last Julian Day in each month
        integer, dimension(12), parameter :: MODAYS = [31, 59, 90, 120, 151,   &
            181, 212, 243, 273, 304, 334, 365]


        ! Data dictionary: calling arguments
        integer,        intent(in)    :: year ! Year of simulation
        type(SiteData), intent(inout) :: site ! Site object

        ! Data dictionary: local variables
        real,    dimension(NTEMPS, 2)     :: tdd            ! Thawing degree-days (>0degC)
        real,    dimension(NTEMPS, 2)     :: fdd            ! Freezing degree-days (<0degC)
        real,    dimension(NTEMPS)        :: tmin           ! Monthly minimum temperature (degC)
        real,    dimension(NTEMPS)        :: tmax           ! Monthly maximum temperature (degC)
        real,    dimension(NTEMPS)        :: prcp           ! Monthly precipitation (cm)
        real,    dimension(NTEMPS)        :: tmean          ! Monthly average temperature (deg)
        real,    dimension(NTEMPS)        :: cld            ! Monthly cloudiness (tenths of sky covered)
        real,    dimension(NTEMPS)        :: rh             ! Monthly relative humidity (%)
        real,    dimension(NTEMPS)        :: wind           ! Monthly wind speed
        real,    dimension(NTEMPS)        :: strikes        ! Monthly lightning (strikes/km2/day)
        real,    dimension(NTEMPS)        :: tmptmin        ! Temporary variable for calculating actual tmin (degC)
        real,    dimension(NTEMPS)        :: tmptmax        ! Temporary variable for calculating actual tmax (degC)
        real,    dimension(NTEMPS)        :: tmpprec        ! Temporary variable for calculating actual prcp (cm)
        real,    dimension(NTEMPS)        :: tmpcld         ! Temporary variable for calculating actual cld (tenths of sky)
        real,    dimension(DAYS_PER_YEAR) :: daytemp        ! Daily temperature (degC)
        real,    dimension(DAYS_PER_YEAR) :: daytemp_min    ! Daily minimum temperature (degC)
        real,    dimension(DAYS_PER_YEAR) :: daytemp_max    ! Daily maximum temperature (degC)
        real,    dimension(DAYS_PER_YEAR) :: daycld         ! Daily cloud cover (tenths of sky covered)
        real,    dimension(DAYS_PER_YEAR) :: dayprecip      ! Daily precipitation (cm)
        real,    dimension(DAYS_PER_YEAR) :: sun            ! Surface solar radiation (cal/cm2/day)
        real,    dimension(DAYS_PER_YEAR) :: st             ! Horizontal surface solar radiation (cal/cm2/day)
        real,    dimension(DAYS_PER_YEAR) :: exrad          ! Top of atmosphere solar radiation (cal/cm2/day)
        real,    dimension(DAYS_PER_YEAR) :: pot_ev_day     ! Potential evapotranspiration (cm)
        character(len = MAX_CHAR)         :: message        ! Error message
        real                              :: rain           ! Annual precipitation (cm)
        real                              :: rain_n         ! Annual N deposition (tN)
        real                              :: temp_f         ! Factor for creating temperature randomness
        real                              :: prcp_f         ! Factor for creating precipitation randomness
        real                              :: cld_f          ! Factor for creating cloud cover randomness
        real                              :: temp_max       ! Maximum temperature of warmest month (degC)
        real                              :: temp_min       ! Mininum temperature of warmest month (degC)
        real                              :: daytemp_mem    ! Average temperature over last N_MEM days
        real                              :: tmean_max      ! Maximum average temperature - for finding warmest month
        real                              :: n_avail        ! Plant-available nitrogen (tN/ha)
        real                              :: pet            ! Potential evapotranspiration (cm)
        real                              :: e1             ! Saturation vapor presstion at tmin of warmest month
        real                              :: e2             ! Saturation vapor presstion at tmax of warmest month
        real                              :: aet            ! Annual actual evapotranspiration (cm)
        real                              :: aet_mm         ! Annual actual evapotranspiration (mm)
        real                              :: growdays       ! Growing season length (days)
        real                              :: soildays       ! Soil degree-days (>0degC)
        real                              :: flooddays      ! Proportion of growing season with flooded conditions
        real                              :: wpdays         ! Proportion of growing season below wilting point
        real                              :: drydays        ! Proportion of growing season with drought conditions
        real                              :: degday         ! Growing degree-days (>5degC)
        real                              :: outwater       ! Runoff (cm)
        real                              :: tot_sun        ! Annual surface solar radiation (cal/cm2/day)
        real                              :: tot_st         ! Annual horizontal surface solar radiation (cal/cm2/day)
        real                              :: cfs            ! Ratio of surface:horizontal surface solar radiation
        real                              :: act_ev_day     ! Actual evapotranspiration (cm)
        real                              :: tcum           ! Cumulative thawing degree-days (>0degC)
        real                              :: fcum           ! Cumulative freezing degree-days (<0degC)
        real                              :: amlt           ! Last year's active layer depth (m)
        real                              :: xmlt           ! Thaw depth (m)
        real                              :: xfrz           ! Freezing depth (m)
        real                              :: zh             ! Soil layer depth
        real                              :: alff           ! Available light on the forest floor (0-1)
        real                              :: pc_germ        ! Effect of temperature on germination (0-1)
        real                              :: aow0_ByMin     ! Organic layer moisture scaled by wilting point
        real                              :: saw0_ByFC      ! Mineral layer moisture scaled by field capacity
        real                              :: saw0_BySAT     ! Mineral layer moisture scaled by saturation capacity
        real                              :: saw0_ByWP      ! Mineral layer moisture scaled by wilting point
        real                              :: saw0_ByFC_sum  ! Sum of mineral layer moisture scaled by field capacity
        real                              :: aow0_ByMin_sum ! Sum of organic layer moisture scaled by wilting point
        real                              :: saw0_BySAT_sum ! Sum of mineral layer moisture scaled by saturation capacity
        real                              :: tmpstep1       ! Temporary variable for implementing linear climate change
        real                              :: tmpstep2       ! Temporary variable for implementing linear climate change
        real                              :: tmp            ! Temporary variable for implementing linear climate change
        integer                           :: gcm_year       ! Year of climate change simulation
        integer                           :: siteid         ! Site ID
        integer                           :: warmest_month  ! Warmest month
        integer                           :: hrise          ! Hour of sunrise
        integer                           :: i, j, m, ip    ! Looping indices
        integer                           :: l

        ! Initialize accumulators
        rain = 0.0
        rain_n = 0.0
        tmean_max = RNVALID

        ! Set site ID - in case we need to warn user
        siteid = site%site_id

        ! Check for and implement climate chnage
        ! The user is expected to input decr_by values as positive
        if (linear_cc) then
            ! Using linear climate change
            if (year .ge. site%gcm_year .and. year .le.                        &
                (site%gcm_year + gcm_duration)) then
                site%accum_tmin = site%accum_tmin + tmin_change
                site%accum_tmax = site%accum_tmax + tmax_change
                do m = 1, NTEMPS
                    tmpstep1 = site%precip(m) + site%accum_precip(m)
                    tmpstep2 = tmpstep1 * precip_change
                    site%accum_precip(m) = site%accum_precip(m) + tmpstep2
                end do
            endif
        else if (use_gcm) then

            ! Using climate change from input file - figure out which year
            ! we are in the file
            gcm_year = start_gcm + year - site%gcm_year
            if (gcm_year .ge. start_gcm .and. gcm_year .le. end_gcm) then

                ! Read in climate change data
                call read_gcm_climate(site%site_id, gcm_year, start_gcm, tmin, &
                    tmax, prcp)

                if (site%site_id .ne. INVALID) then
                    ! Update climate values
                    site%tmin = tmin
                    site%tmax = tmax
                    site%precip = prcp*MM_TO_CM

                    ! Readjust for altitude if needed
                    if (adjust_altitude) then
                        call adjustForAltitude(site)
                    end if
                else
                    ! Problem reading in climate data - tell user, but don't
                    ! update
                    write(message, '(A, I6, A)') "Bad climate data for site ", &
                        siteid, " - using historical climate."
                    call warning(message)
                end if
            endif
        endif

        ! Generate current year's weather from distributions of input climate
        ! data
        do i = 1, NTEMPS
            if (linear_cc) then
                ! Adjust for linear climate change
                tmptmin(i) = site%tmin(i) + site%accum_tmin
                tmptmax(i) = site%tmax(i) + site%accum_tmax
                tmpprec(i) = site%precip(i) + site%accum_precip(i)
                tmpcld(i) = site%cld(i)
            else
                tmptmin(i) = site%tmin(i)
                tmptmax(i) = site%tmax(i)
                tmpprec(i) = site%precip(i)
                tmpcld(i)  = site%cld(i)
            endif

            ! Calculate climate fluctuations
             temp_f = clim_nrand(0.0, 1.0)
             prcp_f = clim_nrand(0.0, 1.0)
             cld_f = clim_nrand(0.0, 1.0)

            ! Adjust
            prcp_f = max(-1.0, min(prcp_f, 1.0))
             temp_f = max(-1.0, min(temp_f, 1.0))
             cld_f = max(-1.0, min(cld_f, 1.0))

             ! Adjust monthly climate vars with std and random numbers
             if    (use_gcm .and. (year .ge. site%gcm_year .and. year .le.        &
                    (site%gcm_year + gcm_duration))) then

                ! Just use input values for that month for precip, temperature,
                ! relative humidity, and lightning
                tmin(i) = tmptmin(i)
                tmax(i) = tmptmax(i)
                prcp(i) = max(tmpprec(i), 0.0)

                ! Adjust for std
                cld(i) = max(tmpcld(i) + cld_f*site%cld_std(i), 0.0)


             else
                ! Adjust for std
                 tmin(i) = tmptmin(i) + temp_f*site%tmin_std(i)
                 tmax(i)= tmptmax(i) + temp_f*site%tmax_std(i)

                ! Can't be less than 0.0
                 cld(i) = max(tmpcld(i) + cld_f*site%cld_std(i), 0.0)
                 prcp(i) = max(tmpprec(i) + prcp_f*site%precip_std(i), 0.0)

             end if

             ! Accumulate precipitation and N deposition
             rain = rain + prcp(i)
             rain_n = rain_n + prcp(i)*PRCP_N

             ! Get mean monthly temperature for warmest month calculation
             tmean(i) = (site%tmin(i) + site%tmax(i))/2.0
         end do

         ! Find warmest month of this year
         do i = 1, NTEMPS
             tmean_max = max(tmean_max, tmean(i))
             if (tmean_max .eq. tmean(i)) warmest_month = i
         end do

        ! Get tmax and tmin of warmest month
         temp_max = site%tmax(warmest_month)
         temp_min = site%tmin(warmest_month)

         ! Calculate e2 and e1 (used for PET calculation)
         e1 = esat(temp_min)
         e2 = esat(temp_max)

         ! Convert monthly weather data into daily weather data
         call cov365_state(tmin, daytemp_min)
         call cov365_state(tmax, daytemp_max)
         call cov365_integr(prcp, dayprecip)
         call cov365_state(cld, daycld)

         ! Initialize accumulators
         pet = 0.0
         tot_sun = 0.0
         tot_st = 0.0
         m = 1

         ! Calculate mean daily temperature, solar radiation, and PET
         do i = 1, DAYS_PER_YEAR

            ! Mean daily temperature (degC)
            daytemp(i) = 0.5*(daytemp_min(i) + daytemp_max(i))

            ! Calculate solar radiation (cal/cm2/day)
            call solar_rad(i, site%latitude, site%slope, site%aspect,         &
                daycld(i), exrad(i), sun(i), st(i), hrise)

            ! Accumulate surface and horizontal surface radiation
            tot_sun = tot_sun + sun(i) ! Actual surface
            tot_st = tot_st + st(i)    ! Horizontal surface

            ! Calculate PET (cm)
            pot_ev_day(i) = pot_evap(daytemp(i), sun(i), site%altitude, e2, e1)

            ! Accumulate PET (cm)
            pet = pet + pot_ev_day(i)
        end do

        ! Calculate ratio of actual surface to horizontal surface radiation
        cfs = tot_sun/tot_st

         ! Calculate freezing and thawing degree days for permafrost subroutine
        tdd = 0.0
        fdd = 0.0
        m = 1
        do j = 1, DAYS_PER_YEAR
            if (j .gt. MODAYS(m)) m = m + 1
            if (tmean(m) > epsilon(1.0) .and. daytemp(j) > epsilon(1.0)) then
                 tdd(m, 1) = tdd(m, 1) + daytemp(j)
            end if
            if (tmean(m) <= epsilon(1.0) .and. daytemp(j) <= epsilon(1.0)) then
                 fdd(m, 1) = fdd(m, 1) + abs(daytemp(j))
            end if
        end do

        ! Calculate cumulative freezing and thawing degree days
        tcum = 0.0
        fcum = 0.0
        do m = 12, 1, -1
            if (fdd(m, 1) .gt. 0.0) fcum = fcum + fdd(m, 1)
            if (fdd(m, 1) .eq. 0.0) exit
        end do
        do m = 1, 12
            if (tdd(m, 1) .eq. 0.0) tcum = 0.0
            tcum = tcum + tdd(m, 1)
            tdd(m, 2) = tcum
            if (fdd(m, 1) .eq. 0.0) fcum = 0.0
            fcum = fcum + fdd(m, 1)
            fdd(m, 2) = fcum
         end do

        ! Loop through each plot to calculate soil dynamics
        do ip = 1, site%numplots

            ! Initialize accumulators
            aet = 0.0
            degday = 0.0
            growdays = 0.0
            soildays = 0.0
            flooddays = 0.0
            drydays = 0.0
            outwater = 0.0
            wpdays = 0.0
            aow0_ByMin_sum = 0.0
            saw0_ByFC_sum = 0.0
            saw0_BySAT_sum = 0.0

            ! Store depth of thaw from previous year (m)
            amlt = min(site%plots(ip)%soil%active, site%plots(ip)%soil%A_depth)
            site%plots(ip)%amlt = amlt

            ! Reset
            site%plots(ip)%soil%active = 0.0
            site%plots(ip)%soil%z_freeze = 0.0

            ! Initialize depths freeze and thaw
            xmlt = 0.0
            xfrz = site%plots(ip)%soil%M_depth +                               &
               site%plots(ip)%soil%O_depth + site%plots(ip)%soil%A_depth

            ! Calculate light on the forest floor
            alff = 1.0*exp(-0.25*site%plots(ip)%cla/plotsize)

            do l = 1, 2
                ! Calculate drainage conditions
                site%plots(ip)%soil%z_drain(l) =                               &
                    (site%plots(ip)%soil%sat(l)*(1.0 - amlt) +                &
                    site%plots(ip)%soil%fc(l)*(amlt - 0.32))/(1.0 - 0.32)

                ! Must be between field capacity and saturation capacity
                site%plots(ip)%soil%z_drain(l) =                               &
                    min(site%plots(ip)%soil%z_drain(l),                        &
                    site%plots(ip)%soil%sat(l))
                site%plots(ip)%soil%z_drain(l) =                               &
                    max(site%plots(ip)%soil%z_drain(l),                        &
                     site%plots(ip)%soil%fc(l))

                ! Set soil to fully saturated at onset of year
                if (l .eq. 1) then
                     site%plots(ip)%soil%wc(l) =                                &
                        site%plots(ip)%soil%z_drain(l)*H2O_D/                  &
                         site%plots(ip)%soil%O_bulk_dens
                    zh = site%plots(ip)%soil%O_depth +                         &
                        site%plots(ip)%soil%M_depth
                else
                    site%plots(ip)%soil%wc(l) =                                &
                        site%plots(ip)%soil%z_drain(l)*H2O_D/                  &
                         site%plots(ip)%soil%A_bulk_dens
                     zh = site%plots(ip)%soil%A_depth
                end if

                ! Soil is completely frozen at start of year
                site%plots(ip)%soil%H2Oice(l) =                                &
                    site%plots(ip)%soil%z_drain(l)*zh
                site%plots(ip)%soil%water(l) = 0.0
                site%plots(ip)%soil%d_melt(l) = 0.0
                site%plots(ip)%soil%d_freeze(l) = zh
                site%plots(ip)%soil%minWC = site%plots(ip)%soil%wc(2)

            end do

            ! Loop on days - thaw depths are calculated monthly so have to
            ! increment the months as well
            m = 1
            do j = 1, DAYS_PER_YEAR

                if (j .gt. MODAYS(m)) m = m + 1

                ! Calculate freeze/thaw depths (xfrz, xmlt) and maximum depths
                ! of freeze and thaw
                call permf(site%plots(ip)%soil, m, 1, alff, tdd, fdd, cfs, xfrz)
                call permf(site%plots(ip)%soil, m, 2, alff, tdd, fdd, cfs, xmlt)

                ! Update maximum depths (z_freeze and active)
                site%plots(ip)%soil%z_freeze = max((xfrz -                     &
                    site%plots(ip)%soil%M_depth -                              &
                    site%plots(ip)%soil%O_depth), site%plots(ip)%soil%z_freeze)
                site%plots(ip)%soil%active = max((xmlt -                       &
                    site%plots(ip)%soil%M_depth -                              &
                    site%plots(ip)%soil%O_depth), site%plots(ip)%soil%active)

                 ! Calculate soil water dynamics for the day
                call moist(site%plots(ip)%soil, site%site_id, ip, year, j,     &
                    daytemp(j), dayprecip(j), pot_ev_day(j),                   &
                    site%leaf_area_ind, site%slope, amlt, xmlt, xfrz, tdd, m,  &
                    act_ev_day, aow0_ByMin, saw0_ByFC, saw0_ByWP, saw0_BySAT)

                ! Update minimum water content for the year
                site%plots(ip)%soil%minWC = min(site%plots(ip)%soil%minWC,     &
                    site%plots(ip)%soil%wc(2))

                ! Accumulate variables
                outwater = outwater + site%plots(ip)%soil%runoff
                aet = act_ev_day + aet

                ! Compute degday, dry days, flood days, and growing season
                ! length (days)
                if (daytemp(j) .ge. MIN_GROW_TEMP) then

                    ! Growing degree-days
                    degday = degday + (daytemp(j) - MIN_GROW_TEMP)

                    ! Growing season length
                    growdays = growdays + 1.0

                    ! For averageing values
                    saw0_ByFC_sum = saw0_ByFC_sum + saw0_ByFC
                    saw0_BySAT_sum = saw0_BySAT_sum + saw0_BySAT
                    aow0_ByMin_sum = aow0_ByMin_sum + aow0_ByMin

                    if (saw0_ByFC .lt. DRY_THRESH) then
                        drydays = drydays + 1.0
                     end if

                    if (aow0_ByMin .lt. MAX_DRY_PARM) then
                        wpdays = wpdays + 1.0
                    end if

                    if (saw0_BySAT .gt. MIN_FLOOD_PARM) then
                       flooddays = flooddays + 1.0
                    endif

                end if

                ! Accumulate soil degree-days
                if (daytemp(j) .ge. 0.0) then
                    soildays = soildays + (daytemp(j) - 0.0)
                end if

            end do

            ! Convert drydays, flooddays, and wpdays to proportion of growing
            ! season
            if (growdays .eq. 0) then
                drydays = 0.0
                flooddays = 0.0
                wpdays = 0.0
            else
                tmp = max(min(rain/pet, 1.0), min(aet/pet, 1.0))
                drydays = ((drydays/growdays) + (1.0 - tmp))/2.0
                flooddays = flooddays/growdays
                wpdays = wpdays/growdays
            endif

            ! Convert aet to mm for decomposition
            aet_mm = aet*10.0

            call moss(site%plots(ip)%soil, alff, site%plots(ip)%cla,        &
                site%plots(ip)%soil%dec_fuel, drydays, site%site_id, ip, year)

            call soiln(site%plots(ip)%soil, aet_mm, site%plots(ip)%cla,        &
                soildays, flooddays, n_avail)

            ! Set fan to 0.0 (used last year's fan for this year's soiln
            ! calculation)
            site%plots(ip)%soil%fan = 0.0

            ! Set plot-level attributes for year-end values
            site%plots(ip)%soil%avail_N = n_avail + rain_n
            site%plots(ip)%saw0_ByFC = saw0_ByFC_sum/growdays
            site%plots(ip)%aow0_ByMin = aow0_ByMin_sum/growdays
            site%plots(ip)%saw0_BySAT = saw0_BySAT_sum/growdays
            site%plots(ip)%soil%runoff = outwater
            site%plots(ip)%act_evap_day = aet
            site%plots(ip)%flood_days = flooddays
            site%plots(ip)%dry_days = drydays
            site%plots(ip)%wilt_days = wpdays

        end do

        !calculate site's base aridity (calculated from first 10 years) and
        !yearly aridity
        if (year .le. 9) then
            if (year .eq. 0) then
                 site%aridity_base = min(rain/pet, 1.0)
            else if (year .gt. 0 .and. year .lt. 9) then
                 site%aridity_base = site%aridity_base + min(rain/pet, 1.0)
            else if (year .eq. 9) then
                 site%aridity_base = (site%aridity_base + min(rain/pet, 1.0))/10.0
            endif
        end if
        site%aridity = min(rain/pet, 1.0)

        ! Set site-level attributes to yearly sums of climate values
        site%deg_days = degday
        site%grow_days = growdays
        site%pot_evap_day = pet
        site%solar = tot_sun
        site%rain = rain

    end subroutine BioGeoClimate

    !:.........................................................................:

    subroutine Canopy(site)
        !
        !  Calculates plot-level leaf area, LAI, and light environment
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong          Original Code
        !    01/01/12     K. Holcomb           Updated to OOP structure

        ! Data dictionary: constants
        real, parameter :: XT = -0.40      ! Light extinction coefficient

        ! Data dictionary: calling arguments
        type(SiteData), intent(inout) :: site ! Site object

        ! Data dictionary: local variables
        real, dimension(maxheight) :: la_dec      ! Leaf area experienced by deciduous plants (m2)
        real, dimension(maxheight) :: la_con      ! Leaf area experienced by evergreen plants (m2)
        real, dimension(maxheight) :: cla_dec     ! Cumulative leaf area experienced by deciduous plants (m2)
        real, dimension(maxheight) :: cla_con     ! Cumulative leaf area experienced by evergreen plants (m2)
        real                       :: forht       ! Tree height (m)
        real                       :: canht       ! Clear branch bole height (m)
        real                       :: tla         ! Tree leaf area (m2)
        real                       :: tla_adj     ! Leaf area per 1-m height bin (m2)
        integer                    :: ntrees      ! Number of trees on plot
        integer                    :: iht         ! Tree height (m)
        integer                    :: cl          ! Canopy length (m)
        integer                    :: i, ip       ! Looping indices
        integer                    :: ih, it
        integer                    :: m

        ! Initialize accumulators
        site%leaf_area_ind = 0.0
        site%lai_array = 0.0

        ! Loop through plots to calculate leaf area of each tree and LAI of
        ! stand
        do ip = 1, site%numplots

            ! Get number of trees on plot
            ntrees = site%plots(ip)%numtrees

            if (ntrees .eq. 0) then

                ! Full light conditions and no nutrient pressure
                site%plots(ip)%con_light = 1.0
                site%plots(ip)%dec_light = 1.0
                site%plots(ip)%fc_nutr = 1.0

            else

                ! Initialize leaf area and canopy biomass arrays
                la_dec = 0.0
                la_con = 0.0
                cla_dec = 0.0
                cla_con = 0.0

                do it = 1, ntrees

                    ! Total tree height (m)
                    forht = max(site%plots(ip)%trees(it)%forska_ht, 2.0)

                    ! Integer of tree height (m)
                    iht = min(int(forht), maxheight)

                    ! Clear branch bole height (m)
                    canht = max(site%plots(ip)%trees(it)%canopy_ht, 1.0)

                    ! Calculate leaf area (m2)
                    tla = leaf_area(site%plots(ip)%trees(it))

                    ! Accumulate site leaf area
                    site%leaf_area_ind = site%leaf_area_ind + tla

                    !Calculate canopy depth and divide leaf area and biomass
                    ! into 1-m sections
                    cl = max(int(forht) - int(canht) + 1, 1)
                    tla_adj = tla/float(cl)

                    ! Fill temporary arrays with leaf area/biomass
                    if (site%plots(ip)%trees(it)%conifer) then
                        ! Leaf experienced by evergreens reduced for deciduous
                        ! plants. This accounts for part of each year without
                        ! decid. leaves
                        do ih = int(canht), int(forht)
                            la_dec(ih) = la_dec(ih) + tla_adj
                            la_con(ih) = la_con(ih) + tla_adj
                        end do
                    else
                        do ih = int(canht), int(forht)
                            la_dec(ih) = la_dec(ih) + tla_adj
                            la_con(ih) = la_con(ih) + tla_adj*0.8
                        end do
                    end if

                end do

                ! Calculate cumulative leaf area from top down
                cla_dec(maxheight) = la_dec(maxheight)
                cla_con(maxheight) = la_con(maxheight)

                do ih = 1, maxheight - 1
                    ! Reduced leaf area from deciduous plants means higher
                    ! light environment for evergreens
                    ! Leaf area for deciduous plants normal b/c evergreen and
                    ! decid both present when decid has leaves
                    cla_dec(maxheight - ih) = cla_dec(maxheight - ih + 1) +    &
                        la_dec(maxheight - ih)
                    cla_con(maxheight - ih) = cla_con(maxheight - ih + 1) +    &
                        la_con(maxheight - ih)
                end do

                ! Calculate light environment
                do ih = 1, maxheight - 1
                    site%plots(ip)%con_light(ih) = exp(XT*cla_con(ih + 1)/     &
                        plotsize)
                    site%plots(ip)%dec_light(ih) = exp(XT*cla_dec(ih + 1)/     &
                        plotsize)
                end do


            end if ! end if any trees

            ! Save plot attributes
            site%plots(ip)%cla = cla_dec(1)
            site%lai_array = site%lai_array + la_dec

        end do !end plot loop

        ! Get average LAI (m2/m2) for site
        site%leaf_area_ind = site%leaf_area_ind/float(site%numplots)/plotsize
        site%lai_array = site%lai_array/float(site%numplots)/plotsize

    end subroutine Canopy

    !:.........................................................................:

    subroutine Growth(site)
        !
        !  Calculates annual growth and branch thinning of each tree
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong          Original Code
        !    01/01/12     K. Holcomb           Updated to OOP structure
        !    10/10/16     A. C. Foster         Updated for soil/plot overhaul
        !                                       and permafrost updates
        !

        ! Data dictionary: constants
        integer, parameter :: MCOUNT = 2

        ! Data dictionary: calling arguments
        type(SiteData), intent(inout) :: site ! Site object

        ! Data dictionary: local variables
        real,    dimension(size(site%species)) :: recr_trees  ! Number of reproductively active trees
        real,    dimension(maxcells*maxcells)  :: diam        ! Tree diameter (dbh)
        real,    dimension(maxcells*maxcells)  :: shade       ! Shade stress at top of tree's canopy (0-1)
        real,    dimension(maxcells*maxcells)  :: can_shade   ! Shade stress at bottom of tree's canopy (0-1)
        real,    dimension(maxcells*maxcells)  :: biomC       ! Woody biomass (tC)
        real,    dimension(maxcells*maxcells)  :: bleaf       ! Leaf biomass (tC)
        real                                   :: ht          ! Tree height (m)
        real                                   :: canht       ! Clear branch bole height(m)
        real                                   :: envstress   ! Growth stress factor (0-1)
        real                                   :: totbiomC    ! Total biomass on plot (tC)
        real                                   :: NPP         ! Net primary production
        real                                   :: d_leafb     ! Change in leaf biomass (tC)
        real                                   :: N_used      ! Nitrogen used by plant growth (tN)
        real                                   :: N_req       ! Nitrogen required for plant growth (tN)
        real                                   :: Npavail     ! Percent N available of required N
        real                                   :: dt          ! Diameter increment (cm)
        real                                   :: mindt       ! Minimum diameter increment before "stressed" (cm)
        real                                   :: dbiomC      ! Change in woody biomass from previous year (tC)
        real                                   :: leafbm      ! Leaf biomass (tC)
        real                                   :: bct         ! Aboveground woody biomass (tC)
        real                                   :: bcr         ! Branch biomass (tC)
        real                                   :: bcs         ! Stem biomass (tC)
        real                                   :: bcbr        ! Total branch biomass (tC)
        real                                   :: d_bc        ! Woody biomass lost to branch thinning (tC)
        real                                   :: d_bcs       ! Stem biomass lost to branch thinning (tC)
        real                                   :: d_bcr       ! Root biomass lost to branch thinning (tC)
        real                                   :: d_bcbr      ! Total branch biomass lost to branch thinning (tC)
        real                                   :: d_bctw      ! Twig biomass lost to branch thinning (tC)
        real                                   :: d_bcsb      ! Small branch biomass lost to branch thinning (tC)
        real                                   :: d_bclb      ! Large branch biomass lost to branch thinning (tC)
        integer                                :: hc          ! Canopy height (m)
        integer                                :: h           ! Tree height (m)
        integer                                :: ntrees      ! Number of trees on plot
        integer                                :: num_species ! Number of species in site
        integer                                :: it, is, ip  ! Looping indices
        integer                                :: lc          ! Litter class

        ! Get number of species at site
        num_species = size(site%species)

        ! Initialize recruiting trees accumulator
        recr_trees = 0.0

        plot: do ip = 1, site%numplots

            ! Initialize accumulators
            N_used = 0.0
            N_req = 0.0
            NPP = 0.0
            totbiomC = 0.0
            site%plots(ip)%mature(:) = 0
            site%plots(ip)%avail_spec = 0.0

            ! Calculate species-level response to drought, over-saturation,
            ! temperature, and permafrost
            do is = 1, num_species
                call temp_rsp(site%species(is), site%deg_days,                 &
                    site%plots(ip)%fc_gdd(is))
                call drought_rsp(site%species(is), site%plots(ip)%dry_days,    &
                    site%plots(ip)%fc_drought(is))
                call perm_rsp(site%species(is)%perm_tol,                       &
                    site%plots(ip)%soil%active, site%plots(ip)%fc_perm(is))
            end do

            ! Get number of trees
            ntrees = site%plots(ip)%numtrees

            numtrees: if (ntrees .gt. 0) then

                stress: do it = 1, ntrees

                    ! Get species index and update tree
                    is = site%plots(ip)%trees(it)%species_index
                    !call update_tree(site%plots(ip)%trees(it), site%species(is))

                    ! Save diameter here
                    diam(it) = site%plots(ip)%trees(it)%diam_bht

                    ! Convenience variables to reduce table lookups
                    canht = site%plots(ip)%trees(it)%canopy_ht
                    ht = site%plots(ip)%trees(it)%forska_ht

                    ! Calculate if species is able to regenerate
                    site%plots(ip)%avail_spec(is) = max(kron(diam(it) -        &
                        site%species(is)%max_diam*site%species(is)%dbh_min),   &
                        site%plots(ip)%avail_spec(is))

                    if (site%plots(ip)%trees(it)%tree_age .ge.                 &
                        site%species(is)%recr_age .and.                        &
                        site%plots(ip)%trees(it)%diam_bht .gt. 10.0) then
                        site%plots(ip)%mature(is) = site%plots(ip)%mature(is) + 1
                    end if

                    ! Get canopy and total tree height as integers
                    h = max(int(ht), 1)
                    hc = max(int(canht), 1)

                    ! Get leaf biomass and maximum possible DBH growth
                    call leaf_biomass_c(site%plots(ip)%trees(it))
                    call max_growth(site%plots(ip)%trees(it))

                    ! Save current value of leaf biomass
                    bleaf(it) = site%plots(ip)%trees(it)%leaf_bm

                    ! Calculate shading effect on tree
                    if (site%plots(ip)%trees(it)%conifer) then
                        shade(it) = light_rsp(site%species(is),                &
                            site%plots(ip)%con_light(h))
                        can_shade(it) = light_rsp(site%species(is),            &
                            site%plots(ip)%con_light(hc))
                    else
                        shade(it) = light_rsp(site%species(is),                &
                            site%plots(ip)%dec_light(h))
                        can_shade(it) = light_rsp(site%species(is),            &
                            site%plots(ip)%dec_light(hc))
                    end if

                    ! Calculate environmental stress (excluding nutrients)
                    call env_stress(site%plots(ip)%trees(it), shade(it),       &
                        site%plots(ip)%fc_gdd(is),                             &
                        site%plots(ip)%fc_drought(is),                         &
                        site%plots(ip)%fc_perm(is), envstress)

                    ! Increment tree diameter using potential DBH growth
                    site%plots(ip)%trees(it)%diam_bht =                        &
                        site%plots(ip)%trees(it)%diam_bht +                    &
                        site%plots(ip)%trees(it)%diam_opt*envstress

                    ! Compute total height for new diameter
                    call forska_height(site%plots(ip)%trees(it))

                    ! Update clear branch bole height and leaf biomass with
                    ! new height
                    call stem_shape(site%plots(ip)%trees(it))
                    call leaf_biomass_c(site%plots(ip)%trees(it))

                    ! Calculate leaf and fine root growth N requirement
                    if (site%species(is)%conifer) then
                        N_req = N_req +                                        &
                            (site%plots(ip)%trees(it)%leaf_bm -                &
                            bleaf(it)*(1.0 - CON_LEAF_RATIO))/CON_LEAF_C_N
                    else
                        N_req = N_req +                                        &
                            site%plots(ip)%trees(it)%leaf_bm/DEC_LEAF_C_N
                    end if

                    ! Store old value of woody biomass
                    biomC(it) = site%plots(ip)%trees(it)%biomC +               &
                        site%plots(ip)%trees(it)%rootC

                    ! Compute new value
                    call biomass_c(site%plots(ip)%trees(it))
                    call biomass_n(site%plots(ip)%trees(it))

                    ! Calculate woody growth N requirement
                    N_req = N_req + ((site%plots(ip)%trees(it)%biomC +         &
                        site%plots(ip)%trees(it)%rootC) - biomC(it))/STEM_C_N

                end do stress

                ! Convert N_req tonnes N/ha
                N_req = max(N_req*HEC_TO_M2/plotsize, epsilon(1.0))

                ! Calculate percent available N
                Npavail = site%plots(ip)%soil%avail_N/N_req

                ! Calculate species-level response to available N
                do is = 1, num_species
                    site%plots(ip)%fc_nutr(is) = poor_soil_rsp(Npavail,        &
                        site%species(is)%lownutr_tol)
                end do

                ! Calculate actual DBH growth and N used
                grow: do it = 1, ntrees

                    ! Get species index and update tree
                    is = site%plots(ip)%trees(it)%species_index
                    !call update_tree(site%plots(ip)%trees(it), site%species(is))

                    ! Calculate environmental stress - including nutrients
                    call env_stress(site%plots(ip)%trees(it), shade(it),       &
                        site%plots(ip)%fc_gdd(is),                             &
                        site%plots(ip)%fc_drought(is),                         &
                        site%plots(ip)%fc_perm(is),                            &
                        envstress, site%plots(ip)%fc_nutr(is))

                    ! Calculate actual diameter increment growth
                    dt = site%plots(ip)%trees(it)%diam_opt*envstress

                    ! Increment old diameter
                    site%plots(ip)%trees(it)%diam_bht = diam(it) + dt

                    ! Get check value for age and growth-related mortality
                    mindt = min(site%species(is)%max_diam/                     &
                        site%species(is)%max_age*0.1, site%species(is)%dbh_min)

                    ! Check for possible mortality age/growth stress mortality
                    if (dt .le. site%species(is)%dbh_min) then

                        ! Diameter growth is below minimum level, increment
                        ! mortality counter
                        site%plots(ip)%trees(it)%mort_count =                  &
                            site%plots(ip)%trees(it)%mort_count + 1

                        if (site%plots(ip)%trees(it)%mort_count .ge. MCOUNT) then
                            ! Tree has been stressed for too many years,
                            ! turn on mortality flag
                            site%plots(ip)%trees(it)%mort_marker = .true.
                        else
                            ! Still possible to live
                            site%plots(ip)%trees(it)%mort_count = 0
                        endif
                    else
                        ! Rest mortality counter and set flag to false
                        site%plots(ip)%trees(it)%mort_count = 0
                        site%plots(ip)%trees(it)%mort_marker = .false.
                    endif

                    ! Compute actual new height and diameter
                    call forska_height(site%plots(ip)%trees(it))
                    call stem_shape(site%plots(ip)%trees(it))

                    ! Update biomass, saving leaf biomass into convenience
                    ! variable
                    call leaf_biomass_c(site%plots(ip)%trees(it))
                    call biomass_c(site%plots(ip)%trees(it))
                    call biomass_n(site%plots(ip)%trees(it))
                    leafbm = site%plots(ip)%trees(it)%leaf_bm

                    ! Calculate change in biomass from previous year
                    dbiomC = (site%plots(ip)%trees(it)%biomC +                 &
                        site%plots(ip)%trees(it)%rootC) - biomC(it)

                    ! Calculate C and N used
                    NPP = NPP + dbiomC
                    N_used = N_used + dbiomC/STEM_C_N

                    ! Update NPP and N_used from leaf biomass
                    if (site%plots(ip)%trees(it)%conifer) then

                        ! Conifers don't have to put all of their leaf biomass
                        ! back on
                        NPP = NPP + leafbm - bleaf(it)*(1.0 - CON_LEAF_RATIO)
                        N_used = N_used + (leafbm -                           &
                            bleaf(it)*(1.0 - CON_LEAF_RATIO))/CON_LEAF_C_N

                        ! Accumulate total biomass
                        totbiomC = totbiomC + site%plots(ip)%trees(it)%biomC + &
                            site%plots(ip)%trees(it)%rootC + leafbm
                    else

                        NPP = NPP + leafbm
                        N_used = N_used + leafbm/DEC_LEAF_C_N

                        ! Accumulate total woody biomass (no leaves)
                        totbiomC = totbiomC + site%plots(ip)%trees(it)%biomC + &
                            site%plots(ip)%trees(it)%rootC

                    end if

                    ! Calculate stand age as maximum tree age
                    site%plots(ip)%stand_age = max(site%plots(ip)%stand_age,   &
                        site%plots(ip)%trees(it)%tree_age)

                    ! Get updated tree height and clear branch bole height
                    ht = site%plots(ip)%trees(it)%forska_ht
                    canht = site%plots(ip)%trees(it)%canopy_ht
                    hc = max(int(canht), 1)

                    ! Check for lower branch thinning
                    ! This will increase clear branch bole height
                    branchfall: if (dt <= site%species(is)%dbh_min) then

                        ! Tree form and will drop some branches

                        ! Increment clear branch bole height
                        hc = hc + 1

                        htcheck: if (hc < int(ht)) then

                            ! Clear branch bole height is still less than tree
                            ! height
                            ! So we can update it without any issues
                            site%plots(ip)%trees(it)%canopy_ht = float(hc) +   &
                                0.01

                            ! Update diameter at clear branch bole height
                            call stem_shape(site%plots(ip)%trees(it))

                            ! Save old woody biomass
                            bct = site%plots(ip)%trees(it)%biomC +             &
                                site%plots(ip)%trees(it)%rootC

                            ! Save old root, twig, stem C biomass
                            bcr = site%plots(ip)%trees(it)%rootC
                            bcbr = site%plots(ip)%trees(it)%branchC
                            bcs = site%plots(ip)%trees(it)%stemC

                            ! Update biomass C and N given new clear branch
                            ! bole height
                            call biomass_c(site%plots(ip)%trees(it))
                            call biomass_n(site%plots(ip)%trees(it))

                            ! How much wood litter did we lose?
                            d_bc = bct - (site%plots(ip)%trees(it)%biomC +     &
                                site%plots(ip)%trees(it)%rootC)
                            d_bcr = bcr - site%plots(ip)%trees(it)%rootC
                            d_bcbr = bcbr - site%plots(ip)%trees(it)%branchC
                            d_bcs = bcs -  site%plots(ip)%trees(it)%stemC

                            ! Divide branch biomass into twigs, large branches,
                            ! and small branches (Thonicke et al. 2010)
                            d_bctw = d_bcbr*PERC_BRANCHES(1)
                            d_bcsb = d_bcbr*PERC_BRANCHES(2)
                            d_bclb = d_bcbr*PERC_BRANCHES(3)

                            ! Add litter loss to litter pools
                            ! Here we convert to dry biomass because soiln
                            ! subroutine calculates weight loss not C loss

                            ! Roots
                            site%plots(ip)%soil%litter(IROOT) =                &
                                site%plots(ip)%soil%litter(IROOT) + d_bcr/B_TO_C

                            ! Twigs
                            site%plots(ip)%soil%litter(ITW) =                   &
                                site%plots(ip)%soil%litter(ITW) + d_bctw/B_TO_C

                            ! Small branches
                            site%plots(ip)%soil%litter(ISBR) =                 &
                                site%plots(ip)%soil%litter(ISBR) + d_bcsb/B_TO_C

                            !large branches
                            site%plots(ip)%soil%litter(ILBR) =                 &
                                site%plots(ip)%soil%litter(ILBR) + d_bclb/B_TO_C

                            ! Tree DBH < 10 cm go into smallwood, otherwise into
                            ! large wood
                            if (site%plots(ip)%trees(it)%diam_bht  .gt.        &
                                10.0) then
                                site%plots(ip)%soil%litter(ILBL) =             &
                                    site%plots(ip)%soil%litter(ILBL) +         &
                                    d_bcs/B_TO_C
                            else
                                site%plots(ip)%soil%litter(ISBL) =             &
                                    site%plots(ip)%soil%litter(ISBL) +         &
                                    d_bcs/B_TO_C
                            end if

                            ! Save previous value of leaf biomass
                            leafbm = site%plots(ip)%trees(it)%leaf_bm

                            ! Update leaf bm and get difference
                            call leaf_biomass_c(site%plots(ip)%trees(it))
                            d_leafb = leafbm - site%plots(ip)%trees(it)%leaf_bm

                            ! Add that litter to correct leaf litter class
                            lc = site%species(is)%litter_class
                            site%plots(ip)%soil%litter(lc) =                   &
                                site%plots(ip)%soil%litter(lc) + d_leafb/B_TO_C

                        end if htcheck

                    end if branchfall

                end do grow

            end if numtrees

            ! Update plot-level soil characteristics
            N_used = N_used*HEC_TO_M2/plotsize ! tN/ha
            site%plots(ip)%NPP = NPP*HEC_TO_M2/plotsize ! tC/ha
            site%plots(ip)%soil%N_used = N_used
            site%plots(ip)%soil%avail_N = max(0.0,                             &
                site%plots(ip)%soil%avail_N - site%plots(ip)%soil%N_used)

        end do plot


    end subroutine Growth

    !:.........................................................................:

    subroutine Mortality(site, year)
        !
        !  Determines which trees die by age, stress, or disturbances, and adds
        !  their biomass to the appropriate soil litter pools
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    01/01/12     K. Holcomb          Updated to OOP structure
        !    10/10/16     A. C. Foster        Updated for soil/plot overhaul
        !                                       and permafrost updates
        !    03/01/17     A. C. Foster        Updated for fire updates
        !

        ! Data dictionary: constants
        real, parameter :: SBR_RT = 0.1 ! Rate of consumption of small branches
        real, parameter :: LBR_RT = 0.1 ! Rate of consumption of large branches
        real, parameter :: BL_RT = 0.05 ! Rate of consumption of boles

        ! Data dictionary: calling argumnets
        type(SiteData),        intent(inout) :: site ! Site object
        integer,               intent(in)    :: year ! Simulation year

        ! Data dictionary: local variables
        integer, dimension(:), allocatable :: dbh_ind          ! Index of sorted tree DBH array
        real,    dimension(:), allocatable :: dbh              ! Tree DBH (cm)
        real                               :: totbiomC         ! Plot-wide biomass (tC)
        real                               :: NPP_loss         ! Loss of NPP from mortality
        real                               :: fan              ! N volatilized by fires (tN)
        real                               :: consRoot         ! Proportion roots consumed by fire (0-1)
        real                               :: wind_prob        ! Random number for windthrow
        real                               :: fire_prob        ! Random number for fire
        real                               :: fire_p           ! Modified fire probability
        real                               :: biomC            ! Total tree biomass (tC)
        real                               :: leaf_bm          ! Leaf biomass (tC)
        real                               :: bcr              ! Root biomass (tC)
        real                               :: N_cons           ! Proportion of N consumed by fire (0-1)
        real                               :: burn             ! Amount of tree burned (tC)
        real                               :: bctw             ! Twig biomass (tC)
        real                               :: bcs              ! Stem biomass (tC)
        real                               :: bcbr             ! Total branch biomass (tC)
        real                               :: bcsb             ! Small branch biomass (tC)
        real                               :: bclb             ! Large branch biomass (tC)
        real                               :: av_fuel          ! Available fuel for fire consumption
        integer                            :: num_species      ! Number of species on site
        integer                            :: it, ip           ! Looping indices
        integer                            :: dt               ! Counter for dead trees
        integer                            :: lt               ! Counter for live trees
        integer                            :: is               ! Species index
        integer                            :: trow             ! Row of tree
        integer                            :: tcol             ! Column of tree
        logical                            :: age_survive      ! Does tree survive age check?
        logical                            :: growth_survive   ! Does tree survive stress check?
        logical                            :: fire_survive     ! Does tree survive fire?
        logical                            :: wind_survive     ! Does tree survive windthrow?
        integer                            :: snum             ! Tree growth stressor
        integer                            :: lc               ! Litter class

        ! Get number of species on site
        num_species = size(site%species)

        plot: do ip = 1, site%numplots

            ! Initialize accumulators
            site%plots(ip)%num_dead = 0
            totbiomC = 0.0
            NPP_loss = 0.0
            fan = 0.0
            site%plots(ip)%d_type = 0.0

            ! Set fire and wind to 0
            site%plots(ip)%fire = 0
            site%plots(ip)%wind = 0

            ! Get random number for fire and wind throw
            wind_prob = urand()
            fire_prob = urand()

            !increase or decrease fire probability based on aridity
            if (year >= 10) then
                if (site%aridity .lt. site%aridity_base) then
                    fire_p = site%fire_prob + ((site%aridity_base -            &
                        site%aridity)/site%aridity_base)*site%fire_prob
                else if (site%aridity .gt. site%aridity_base) then
                    fire_p = site%fire_prob - ((site%aridity -                 &
                        site%aridity_base)/site%aridity_base)*site%fire_prob
                else
                    fire_p = site%fire_prob
                end if
            else
                fire_p = site%fire_prob
            endif

            treatments: if (fire_prob < fire_p) then

                ! We have a fire on the plot

                fntrees: if (site%plots(ip)%numtrees > 0) then

                    site%plots(ip)%fire = 1

                    ! Calculate and consume available litter fuels
                    call forest_fuels(site%plots(ip)%soil,                     &
                        site%plots(ip)%dry_days, av_fuel, N_cons, consRoot)

                    ! Kill trees that died by fire, age, or low growth - only
                    ! copy surviving trees, rest go into soil or burn

                    ! Initialize dead and live tree counters
                    dt = 0
                    lt = 0

                    fireloop: do it = 1, site%plots(ip)%numtrees

                        ! Get species index
                        is = site%plots(ip)%trees(it)%species_index

                        ! Get leaf biomass
                        call leaf_biomass_c(site%plots(ip)%trees(it))
                        leaf_bm = site%plots(ip)%trees(it)%leaf_bm

                        ! Check for growth and age survival
                        call growth_survival(site%plots(ip)%trees(it),         &
                            growth_survive)
                        call age_survival(site%plots(ip)%trees(it),            &
                            age_survive)

                        ! Check for fire survival
                        call fire_survival(site%plots(ip)%trees(it), av_fuel,  &
                            fire_survive)

                        fdeathcheck: if (growth_survive .and.                   &
                            age_survive .and. fire_survive) then

                            ! Tree survives

                            lt = lt + 1 ! Increment live tree counter

                            ! Increment tree age
                            site%plots(ip)%trees(it)%tree_age =                &
                            site%plots(ip)%trees(it)%tree_age + 1

                            ! Copy tree to front of list
                            call copy_tree(site%plots(ip)%trees(lt),           &
                                site%plots(ip)%trees(it))

                            ! Calculate leaf litter
                            if (site%species(is)%conifer) then
                                ! Conifers only drop some of their needles

                                lc = site%species(is)%litter_class
                                site%plots(ip)%soil%litter(lc) =               &
                                    site%plots(ip)%soil%litter(lc) +           &
                                    leaf_bm*CON_LEAF_RATIO/B_TO_C

                                ! Accumulate NPP losses
                                NPP_loss = NPP_loss + leaf_bm*CON_LEAF_RATIO

                            else

                                lc = site%species(is)%litter_class
                                site%plots(ip)%soil%litter(lc) =               &
                                    site%plots(ip)%soil%litter(lc) +           &
                                    leaf_bm/B_TO_C

                                ! Acumulate NPP losses
                                NPP_loss = NPP_loss + leaf_bm

                            end if

                        else fdeathcheck

                            ! Tree dies from something (need to check which)

                            dt = dt + 1 ! Increment dead tree counter

                            ! Copy to dead tree list for output
                            call copy_tree(site%plots(ip)%deadtrees(dt),       &
                                site%plots(ip)%trees(it))

                            ! Set cells of that tree to unfilled
                            trow = site%plots(ip)%trees(it)%row
                            tcol = site%plots(ip)%trees(it)%col
                            site%plots(ip)%cells(trow, tcol) = 0

                            ! Get most limiting growth factor
                            snum = site%plots(ip)%trees(it)%stressor

                            firecheck: if (.not. fire_survive) then

                                ! Died by fire, fire consumes some of each
                                !litter class. also calculate N volatilized in
                                ! tree burning

                                ! Add biomass to array of dead biomass
                                site%plots(ip)%d_type(IFIRE) =                 &
                                    site%plots(ip)%d_type(IFIRE) +             &
                                    site%plots(ip)%trees(it)%biomC + leaf_bm

                                ! Update "stressor"
                                site%plots(ip)%deadtrees(dt)%stressor = IFIRE

                                ! Roots - fire consumption from Bonan (1989)
                                bcr = site%plots(ip)%trees(it)%rootC
                                burn = bcr*consRoot
                                bcr = bcr - burn
                                site%plots(ip)%soil%litter(IROOT) =            &
                                    site%plots(ip)%soil%litter(IROOT) +        &
                                    bcr/B_TO_C

                                ! Accumulate volatilized N
                                fan = fan + burn*litter_params(IROOT, 2)*      &
                                    (1.0 - N_cons)

                                ! Branch biomass
                                bcbr = site%plots(ip)%trees(it)%branchC

                                ! Convert branch litter into twigs, small
                                ! branches, and large branches (Thonicke et al. 2010)
                                bctw = bcbr*PERC_BRANCHES(1)
                                bcsb = bcbr*PERC_BRANCHES(2)
                                bclb = bcbr*PERC_BRANCHES(3)

                                ! Twigs
                                burn = bctw*(site%plots(ip)%trees(it)%CK)
                                bctw = bctw - burn
                                site%plots(ip)%soil%litter(ITW) =              &
                                    site%plots(ip)%soil%litter(ITW) +          &
                                    bctw/B_TO_C

                                ! Accumulate volatilized N
                                fan = fan + burn*litter_params(ITW, 2)*        &
                                    (1.0 - N_cons)

                                ! Small branches
                                burn = bcsb*(SBR_RT*site%plots(ip)%trees(it)%CK)
                                bcsb = bcsb - burn
                                site%plots(ip)%soil%litter(ISBR) =             &
                                    site%plots(ip)%soil%litter(ISBR) +         &
                                    bcsb/B_TO_C

                                ! Accumulate volatilized N
                                fan = fan + burn*litter_params(ISBR, 2)*       &
                                    (1.0 - N_cons)

                                ! Large branches
                                burn = bclb*(LBR_RT*site%plots(ip)%trees(it)%CK)
                                bclb = bclb - burn
                                site%plots(ip)%soil%litter(ILBR) =             &
                                    site%plots(ip)%soil%litter(ILBR) +         &
                                    bclb/B_TO_C

                                ! Accumulate volatilized N
                                fan = fan + burn*litter_params(ILBR, 2)*       &
                                    (1.0 - N_cons)

                                ! Stems
                                bcs = site%plots(ip)%trees(it)%stemC
                                burn = bcs*(BL_RT*site%plots(ip)%trees(it)%CK)
                                bcs = bcs - burn

                                ! Small boles (DBH < 10) vs. large boles
                                if (site%plots(ip)%trees(it)%diam_bht          &
                                    > 10.0) then
                                    site%plots(ip)%soil%litter(ILBL) =         &
                                        site%plots(ip)%soil%litter(ILBL) +     &
                                        bcs/B_TO_C

                                    ! Accumulate volatilized N
                                    fan = fan + burn*litter_params(ILBL, 2)*   &
                                        (1.0 - N_cons)
                                else
                                    site%plots(ip)%soil%litter(ISBL) =         &
                                        site%plots(ip)%soil%litter(ISBL) +     &
                                        bcs/B_TO_C

                                    ! Accumulate volatilized N
                                    fan = fan + burn*litter_params(ISBL, 2)*   &
                                        (1.0 - N_cons)
                                end if

                                ! Leaves
                                lc = site%species(is)%litter_class
                                burn = leaf_bm*(site%plots(ip)%trees(it)%CK)
                                leaf_bm = leaf_bm - burn

                                site%plots(ip)%soil%litter(lc) =               &
                                site%plots(ip)%soil%litter(lc) + leaf_bm/B_TO_C

                                ! Accumulate volatilized N
                                fan = fan + burn*litter_params(lc, 2)*         &
                                    (1.0 - N_cons)

                            else if (.not. growth_survive .or.                 &
                                .not. age_survive) then firecheck

                                ! Died from growth or age-related stress, all
                                ! litter goes into soil

                                ! Add biomass to array of dead biomass
                                site%plots(ip)%d_type(snum) =                  &
                                    site%plots(ip)%d_type(snum) +              &
                                    site%plots(ip)%trees(it)%biomC + leaf_bm

                                ! Add total tree litter components to all
                                ! litter categories

                                ! Roots
                                bcr = site%plots(ip)%trees(it)%rootC
                                site%plots(ip)%soil%litter(IROOT) =            &
                                    site%plots(ip)%soil%litter(IROOT) +        &
                                    bcr/B_TO_C

                                ! Branches
                                bcbr = site%plots(ip)%trees(it)%branchC
                                ! Convert branch litter into twigs, small
                                ! branches, and large branches
                                bctw = bcbr*PERC_BRANCHES(1)
                                bcsb = bcbr*PERC_BRANCHES(2)
                                bclb = bcbr*PERC_BRANCHES(3)

                                ! Twigs
                                site%plots(ip)%soil%litter(ITW) =              &
                                    site%plots(ip)%soil%litter(ITW) +          &
                                    bctw/B_TO_C

                                ! Small branches
                                site%plots(ip)%soil%litter(ISBR) =             &
                                    site%plots(ip)%soil%litter(ISBR) +         &
                                    bcsb/B_TO_C

                                ! Large branches
                                site%plots(ip)%soil%litter(ILBR) =             &
                                    site%plots(ip)%soil%litter(ILBR) +         &
                                    bclb/B_TO_C

                                ! Stems
                                bcs = site%plots(ip)%trees(it)%stemC

                                ! Small boles (DBH < 10) vs. large boles
                                if (site%plots(ip)%trees(it)%diam_bht          &
                                    > 10.0) then
                                    site%plots(ip)%soil%litter(ILBL) =         &
                                        site%plots(ip)%soil%litter(ILBL) +     &
                                        bcs/B_TO_C
                                else
                                    site%plots(ip)%soil%litter(ISBL) =         &
                                        site%plots(ip)%soil%litter(ISBL) +     &
                                        bcs/B_TO_C
                                end if

                                ! Leaves
                                lc = site%species(is)%litter_class
                                site%plots(ip)%soil%litter(lc) =               &
                                      site%plots(ip)%soil%litter(lc) +           &
                                    leaf_bm/B_TO_C

                            end if firecheck

                            ! Get biomass of tree
                            biomC = site%plots(ip)%trees(it)%biomC

                            ! Acumulate NPP losses
                            NPP_loss = NPP_loss + biomC + leaf_bm

                        end if  fdeathcheck

                    end do fireloop

                    ! Set number of live and dead trees
                    site%plots(ip)%numtrees = lt
                    site%plots(ip)%num_dead = dt

                end if fntrees

            else if (wind_prob < site%wind_prob) then treatments

                ! We have a windthrow event

                ! Set fire to 0 and wind to 1
                site%plots(ip)%fire = 0
                site%plots(ip)%wind = 1

                ! Set wind counter for regeneration
                site%plots(ip)%windCount = 3

                wntrees: if (site%plots(ip)%numtrees > 0) then

                    ! Initialize counters for live and dead trees
                    lt = 0
                    dt = 0

                    windloop: do it = 1, site%plots(ip)%numtrees

                        ! Get species index and update tree
                        is = site%plots(ip)%trees(it)%species_index
                        !call update_tree(site%plots(ip)%trees(it),             &
                        !    site%species(is))

                        ! Get leaf biomass
                        call leaf_biomass_c(site%plots(ip)%trees(it))
                        leaf_bm = site%plots(ip)%trees(it)%leaf_bm

                        ! Check for growth and age mortality
                        call growth_survival(site%plots(ip)%trees(it),         &
                            growth_survive)
                        call age_survival(site%plots(ip)%trees(it),            &
                            age_survive)

                        ! Check for wind survival
                        call wind_survival(site%plots(ip)%trees(it),           &
                            wind_survive)

                        wdeathcheck: if (growth_survive .and. age_survive      &
                            .and. wind_survive) then

                            ! Tree survives

                            lt = lt + 1 ! Increment live tree counter

                            ! Increment tree age
                            site%plots(ip)%trees(it)%tree_age =                &
                                site%plots(ip)%trees(it)%tree_age + 1

                            ! Copy tree object to top of list
                            call copy_tree(site%plots(ip)%trees(lt),           &
                                site%plots(ip)%trees(it))

                            ! Calculate leaf litter
                            if (site%species(is)%conifer) then

                                ! Evergreens only lose some of their leaves
                                lc = site%species(is)%litter_class
                                site%plots(ip)%soil%litter(lc) =               &
                                    site%plots(ip)%soil%litter(lc) +           &
                                    leaf_bm*CON_LEAF_RATIO/B_TO_C

                                ! Acumulate NPP losses
                                NPP_loss = NPP_loss + leaf_bm*CON_LEAF_RATIO

                            else

                                lc = site%species(is)%litter_class
                                site%plots(ip)%soil%litter(lc) =               &
                                    site%plots(ip)%soil%litter(lc) +           &
                                    leaf_bm/B_TO_C

                                ! Accumulate NPP losses
                                NPP_loss = NPP_loss + leaf_bm
                            end if

                        else wdeathcheck

                            ! Tree dies
                            dt = dt + 1 ! Increment dead tree counter

                            ! Copy to list of dead trees for output
                            call copy_tree(site%plots(ip)%deadtrees(dt),       &
                                site%plots(ip)%trees(it))

                            ! Set cells of that tree to unfilled
                            trow = site%plots(ip)%trees(it)%row
                            tcol = site%plots(ip)%trees(it)%col
                            site%plots(ip)%cells(trow, tcol) = 0

                            ! Get most limiting growth factor
                            snum = site%plots(ip)%trees(it)%stressor

                            windcheck: if (.not. growth_survive  .or.          &
                                .not. age_survive) then

                                ! Tree died from low growth or age

                                ! Add biomass to mortality array
                                site%plots(ip)%d_type(snum) =                  &
                                    site%plots(ip)%d_type(snum) +              &
                                    site%plots(ip)%trees(it)%biomC + leaf_bm

                            else windcheck

                                ! Died by windthrow

                                ! Add biomass to mortality array
                                site%plots(ip)%d_type(IWIND) =                 &
                                    site%plots(ip)%d_type(IWIND) +             &
                                    site%plots(ip)%trees(it)%biomC + leaf_bm

                                site%plots(ip)%deadtrees(dt)%stressor = IWIND

                            end if windcheck

                            ! Add total tree litter components to all
                            ! litter categories

                            ! Roots
                            bcr = site%plots(ip)%trees(it)%rootC
                            site%plots(ip)%soil%litter(IROOT) =                &
                                site%plots(ip)%soil%litter(IROOT) + bcr/B_TO_C

                            ! Branches
                            bcbr = site%plots(ip)%trees(it)%branchC
                            ! Convert branch litter into twigs, small
                            ! branches, and large branches
                            bctw = bcbr*PERC_BRANCHES(1)
                            bcsb = bcbr*PERC_BRANCHES(2)
                            bclb = bcbr*PERC_BRANCHES(3)

                            ! Twigs
                            site%plots(ip)%soil%litter(ITW) =                  &
                                site%plots(ip)%soil%litter(ITW) + bctw/B_TO_C

                            ! Small branches
                            site%plots(ip)%soil%litter(ISBR) =                 &
                                site%plots(ip)%soil%litter(ISBR) + bcsb/B_TO_C

                            ! Large branches
                            site%plots(ip)%soil%litter(ILBR) =                 &
                                site%plots(ip)%soil%litter(ILBR) + bclb/B_TO_C

                            ! Stems
                            bcs = site%plots(ip)%trees(it)%stemC

                            ! Small boles (DBH < 10) vs. large boles
                            if (site%plots(ip)%trees(it)%diam_bht > 10.0) then
                                site%plots(ip)%soil%litter(ILBL) =             &
                                    site%plots(ip)%soil%litter(ILBL) +         &
                                    bcs/B_TO_C
                            else
                                site%plots(ip)%soil%litter(ISBL) =             &
                                    site%plots(ip)%soil%litter(ISBL) +         &
                                    bcs/B_TO_C
                            end if

                            ! Leaves
                            lc = site%species(is)%litter_class
                            site%plots(ip)%soil%litter(lc) =                   &
                                site%plots(ip)%soil%litter(lc) + leaf_bm/B_TO_C

                            ! Acumulate NPP losses
                            NPP_loss = NPP_loss +                              &
                                site%plots(ip)%trees(it)%biomC  + leaf_bm

                        end if wdeathcheck
                    end do windloop

                        ! Set number of live and trees
                        site%plots(ip)%numtrees = lt
                        site%plots(ip)%num_dead = dt

                end if wntrees

            else treatments

                ! No disturbances - just check for age and growth stress
                ! mortality

                ! Set fire and wind to 0
                site%plots(ip)%fire = 0
                site%plots(ip)%wind = 0

                gntrees: if (site%plots(ip)%numtrees > 0) then

                    ! Initialize counters for live and dead trees
                    lt = 0
                    dt = 0

                    deathloop: do it = 1, site%plots(ip)%numtrees

                        ! Get species index and update tree
                        is = site%plots(ip)%trees(it)%species_index
                        !call update_tree(site%plots(ip)%trees(it),             &
                        !    site%species(is))

                        ! Get leaf biomass
                        call leaf_biomass_c(site%plots(ip)%trees(it))
                        leaf_bm = site%plots(ip)%trees(it)%leaf_bm

                        ! Check for age and growth survival
                        call growth_survival(site%plots(ip)%trees(it),         &
                            growth_survive)
                        call age_survival(site%plots(ip)%trees(it),            &
                            age_survive)

                        deathcheck: if (growth_survive .and. age_survive) then

                            ! Tree survives

                            lt = lt + 1 ! Increment live tree counter

                            ! Increment tree age
                            site%plots(ip)%trees(it)%tree_age=                 &
                                site%plots(ip)%trees(it)%tree_age + 1

                            ! Copy tree to top of list
                            call copy_tree(site%plots(ip)%trees(lt),           &
                                site%plots(ip)%trees(it))

                            ! Calculate litterfall
                            if (site%species(is)%conifer) then

                                ! Conifers only drop some of their leaves
                                lc = site%species(is)%litter_class
                                site%plots(ip)%soil%litter(lc) =               &
                                    site%plots(ip)%soil%litter(lc) +           &
                                    leaf_bm*CON_LEAF_RATIO/B_TO_C

                                ! Accumulate NPP losses
                                NPP_loss = NPP_loss + leaf_bm*CON_LEAF_RATIO

                            else

                                lc = site%species(is)%litter_class
                                site%plots(ip)%soil%litter(lc) =               &
                                    site%plots(ip)%soil%litter(lc) +           &
                                    leaf_bm/B_TO_C

                                ! Accumulate NPP losses
                                NPP_loss = NPP_loss + leaf_bm

                            end if

                        else deathcheck

                            ! Tree dies

                            dt = dt + 1

                            ! Copy to list of dead trees for output
                            call copy_tree(site%plots(ip)%deadtrees(dt),       &
                                site%plots(ip)%trees(it))

                            ! Set cells of that tree to unfilled
                            trow = site%plots(ip)%trees(it)%row
                            tcol = site%plots(ip)%trees(it)%col
                            site%plots(ip)%cells(trow, tcol) = 0

                            ! Get most limiting growth factor
                            snum = site%plots(ip)%trees(it)%stressor

                            ! Add biomass to mortality array
                            site%plots(ip)%d_type(snum) =                      &
                                site%plots(ip)%d_type(snum) +                  &
                                site%plots(ip)%trees(it)%biomC + leaf_bm

                            ! Add total tree litter components to all
                            ! litter categories

                            ! Roots
                            bcr = site%plots(ip)%trees(it)%rootC
                            site%plots(ip)%soil%litter(IROOT) =                &
                                site%plots(ip)%soil%litter(IROOT) + bcr/B_TO_C

                            ! Branches
                            bcbr = site%plots(ip)%trees(it)%branchC
                            ! Convert branch litter into twigs, small
                            ! branches, and large branches
                            bctw = bcbr*PERC_BRANCHES(1)
                            bcsb = bcbr*PERC_BRANCHES(2)
                            bclb = bcbr*PERC_BRANCHES(3)

                            ! Twigs
                            site%plots(ip)%soil%litter(ITW) =                  &
                                site%plots(ip)%soil%litter(ITW) + bctw/B_TO_C

                            ! Small branches
                            site%plots(ip)%soil%litter(ISBR) =                 &
                                site%plots(ip)%soil%litter(ISBR) + bcsb/B_TO_C

                            ! Large branches
                            site%plots(ip)%soil%litter(ILBR) =                 &
                                site%plots(ip)%soil%litter(ILBR) + bclb/B_TO_C

                            ! Stems
                            bcs = site%plots(ip)%trees(it)%stemC

                            ! Small boles (DBH < 10) vs. large boles
                            if (site%plots(ip)%trees(it)%diam_bht > 10.0) then
                                site%plots(ip)%soil%litter(ILBL) =             &
                                    site%plots(ip)%soil%litter(ILBL) +         &
                                    bcs/B_TO_C
                            else
                                site%plots(ip)%soil%litter(ISBL) =             &
                                    site%plots(ip)%soil%litter(ISBL) +         &
                                    bcs/B_TO_C
                            end if

                            ! Leaves
                            lc = site%species(is)%litter_class
                            site%plots(ip)%soil%litter(lc) =                   &
                                site%plots(ip)%soil%litter(lc) + leaf_bm/B_TO_C

                            ! Acumulate NPP losses
                            NPP_loss = NPP_loss +                              &
                                site%plots(ip)%trees(it)%biomC  + leaf_bm

                        end if deathcheck

                    end do deathloop

                    ! Set number of live and dead trees
                    site%plots(ip)%num_dead = dt
                    site%plots(ip)%numtrees = lt

                end if gntrees

            end if treatments

            ! Update N volatilized N (tN/ha)
            site%plots(ip)%soil%fan = fan/plotsize/M2_TO_HEC

            ! Update NPP (tC/ha)
            site%plots(ip)%NPP = site%plots(ip)%NPP -                          &
                NPP_loss/plotsize/M2_TO_HEC

        end do plot

    end subroutine Mortality

    !:.........................................................................:

    subroutine Renewal(site, year)
        !
        !  Calculates the seedling and seed banks of each species and
        !  regenerates new trees and shrubs
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    01/01/12     K. Holcomb          Updated to OOP structure
        !    10/10/16     A. C. Foster        Updated for soil/plot overhaul
        !                                       and permafrost updates
        !

        ! Data dictionary: calling arguments
        type(SiteData), intent(inout) :: site ! Site object
        integer,        intent(in)    :: year ! Year of simulation

        ! Data dictionary: constants
        real, parameter :: regmin = 0.01 ! Minimum threshold for regeneration

        ! Data dictionary: local variables
        real,    dimension(size(site%species)) :: regstress   ! Regeneration stress factor (0-1)
        real,    dimension(site%num_trees)     :: probt       ! Probability of regeneration (trees)
        integer, dimension(site%num_trees)     :: tree_sp     ! Id locations of tree species
        integer, dimension(maxcells*maxcells)  :: r_empty     ! Rows of empty cells
        integer, dimension(maxcells*maxcells)  :: c_empty     ! Columns of empty cells
        integer, dimension(:), allocatable     :: locs        ! Locations of empty cells
        real                                   :: NPP         ! Net primary productivity (tC/ha)
        real                                   :: N_used      ! N used in regeneration (tN/ha)
        real                                   :: tregmax     ! Maximum growth stress factor for trees
        real                                   :: fc_org      ! Effect of organic layer depth on regeneration
        real                                   :: germinants  ! Number of seedlings regenerating (#/m2)
        real                                   :: probtsum    ! Sum of probability of regeneration for trees
        real                                   :: rand        ! Random number for determining species (uniform)
        real                                   :: dbh         ! Tree dbh (cm)
        real                                   :: leafbm      ! Leaf biomass (tC)
        real                                   :: shade       ! Shade stress (0-1)
        real                                   :: envstress   ! Growth stress factor (0-1)
        integer                                :: num_species ! Number of species at site
        integer                                :: new_trees   ! Number of trees regenerated
        integer                                :: numtrees    ! Number of trees on plot
        integer                                :: max_trenew  ! Max number of trees to renew
        integer                                :: ntrenew     ! Number of trees to renew
        integer                                :: org_tol     ! Ability to regenerate on deep soil
        integer                                :: n_empty     ! Number of empty cells on plot
        integer                                :: ht          ! Tree height (m)
        integer                                :: lc          ! Litter class
        integer                                :: is, ip, it  ! Looping indices
        integer                                :: stp         ! Species counters
        integer                                :: t, r, c     ! Looping indices
        integer                                :: irenew, i   ! Looping indices

        ! Get number of species at site
        num_species = size(site%species)

        plot: do ip = 1, site%numplots

            ! Initialize accumulators
            new_trees = 0
            NPP = 0.0
            N_used = 0.0
            tregmax = 0.0

            numtrees = site%plots(ip)%numtrees

            navail: if (site%plots(ip)%soil%avail_N > epsilon(1.0)) then

                ! Check to make sure we aren't still waiting on wind counter
                windcheck: if (site%plots(ip)%windCount == 0) then

                    ! Calculate species-level environmental stressors for
                    ! the plot
                    do is = 1, num_species

                        ! First get minimum of gdd, nutrient, and drought
                        ! stress
                        regstress(is) = min(site%plots(ip)%fc_gdd(is),         &
                            site%plots(ip)%fc_nutr(is),                        &
                            site%plots(ip)%fc_drought(is))

                        if (site%plots(ip)%numtrees .eq. 0) then
                            ! We don't have to worry about shade stress
                            regstress(is) = regstress(is)*                     &
                                site%plots(ip)%fc_perm(is)
                        else
                            ! Need to consider shading
                            if (site%species(is)%conifer) then

                                regstress(is) = min(regstress(is),             &
                                    light_rsp(site%species(is),                &
                                    site%plots(ip)%con_light(1)))*             &
                                    site%plots(ip)%fc_perm(is)
                            else
                                regstress(is) = min(regstress(is),             &
                                    light_rsp(site%species(is),                &
                                    site%plots(ip)%dec_light(1)))*             &
                                    site%plots(ip)%fc_perm(is)

                            end if
                        end if

                        ! Check for enough mature trees
                        if (site%plots(ip)%mature(is) <= 5 .and.               &
                            site%species(is)%conifer) then
                            regstress(is) = regstress(is)*0.5
                        end if

                        ! Can't regenerate if below minimum
                        if (regstress(is) <= site%species(is)%dbh_min) then
                            regstress(is) = 0.0
                        end if

                        ! Calculate maximum regrowth capacity across all species
                        tregmax = max(tregmax, regstress(is))

                    end do

                    ! Compute the max renew number trees and shrubs
                    max_trenew = max(min(int(maxtrees*tregmax) - numtrees,     &
                        maxtrees), 0)

                    ! Compute actual number of renewable trees and shrubs
                    ntrenew = min(max_trenew, maxtrees - numtrees)

                    ! Update the seed bank size and seedling bank size (#/m2)
                    seedbank: do is = 1, num_species

                        ! Get species-level regeneration response to fire
                        call fire_rsp(site%species(is), site%plots(ip)%fire,   &
                            site%plots(ip)%fc_fire(is))

                        site%plots(ip)%seedbank(is) =                          &
                            site%plots(ip)%seedbank(is) +                      &
                            site%species(is)%invader +                         &
                            site%species(is)%seed_num*                         &
                            site%plots(ip)%avail_spec(is)*                     &
                            site%plots(ip)%fc_fire(is)

                        ! We don't allow seedling regeneration the first
                        ! year of a fire
                        firecheck: if (site%plots(ip)%fire == 0) then

                            ! Put seeds into seedling bank if envstress is
                            ! high enough
                            seedcheck: if (regstress(is) > site%species(is)%dbh_min) then

                                ! Ability to reproduce on moss-covered soil
                                org_tol = site%species(is)%org_tol
                                fc_org = exp(ORG_GF(org_tol)*                  &
                                    (site%plots(ip)%soil%O_depth +             &
                                    site%plots(ip)%soil%M_depth))

                                ! Calculate the number of new seedlings
                                germinants = site%plots(ip)%seedbank(is)*      &
                                    fc_org

                                ! Remove new seedlings from seedbank
                                site%plots(ip)%seedbank(is) =                  &
                                    site%plots(ip)%seedbank(is) - germinants

                                ! Add germinants to seedling bank
                                site%plots(ip)%seedling(is) =                  &
                                    site%plots(ip)%seedling(is) + germinants

                                ! Decrease seedbank for survival
                                site%plots(ip)%seedbank(is) =                  &
                                    site%plots(ip)%seedbank(is)*               &
                                    site%species(is)%seed_surv

                            else seedcheck

                                ! No new seedlings from seedbank

                                ! Decrease seedbank for survival
                                site%plots(ip)%seedbank(is) =                  &
                                    site%plots(ip)%seedbank(is)*               &
                                    site%species(is)%seed_surv

                            endif seedcheck

                            ! Add seedlings from layering
                            if (site%species(is)%layering) then
                                if ((site%plots(ip)%soil%M_depth +             &
                                    site%plots(ip)%soil%O_depth) > 0.05) then
                                    site%plots(ip)%seedling(is) =              &
                                        site%plots(ip)%seedling(is)*1.8
                                end if
                            end if

                            ! Add seedlings from sprouts
                            site%plots(ip)%seedling(is) =                      &
                                site%plots(ip)%seedling(is) +                  &
                                site%species(is)%sprout_num*                   &
                                site%plots(ip)%avail_spec(is)

                            ! Convert seedling bank to #/plot
                            site%plots(ip)%seedling(is) =                      &
                                site%plots(ip)%seedling(is)*plotsize

                        end if firecheck

                    end do seedbank

                    ! Calculate probability of regeneration
                    probtsum = 0.0
                    stp = 1
                    do is = 1, num_species
                        probt(stp) = site%plots(ip)%seedling(is)*regstress(is)
                        probtsum = probtsum + probt(stp)
                        tree_sp(stp) = is
                        stp = stp + 1
                    end do

                else windcheck

                    ! We are still waiting after a windthrow event
                    ! Decrease counter
                    site%plots(ip)%windCount = max(0,                          &
                        site%plots(ip)%windCount - 1)

                    ! Set probsum to 0
                    probtsum = 0.0

                end if windcheck

                ! After setting seed and seedling banks

                ! Calculate cumulative probability of regeneration for trees
                if (probtsum .gt. epsilon(1.0)) then
                    do is = 1, site%num_trees
                        probt(is) = probt(is)/probtsum
                    end do
                    do is = 2, site%num_trees
                        probt(is) = probt(is - 1) + probt(is)
                    end do
                else
                    ntrenew = 0
                end if

                ! Get current number of trees on plot
                it = site%plots(ip)%numtrees

                renew: if (ntrenew >= 1) then

                    ! Count number of unfilled cells
                    n_empty = count(site%plots(ip)%cells(:,:) == 0)

                    nempty: if (n_empty > 0) then

                        ! Some cells are unfilled - can place trees
                        allocate(locs(n_empty))

                        ! Loop through whole rows and columns and fill
                        ! r_empty and c_empty with empty cell indices -
                        ! keeping them together with the same 't' index
                        t = 1
                        do while (t <= n_empty)
                            do r = 1, maxcells
                                do c = 1, maxcells
                                    if (site%plots(ip)%cells(r, c) == 0) then
                                        r_empty(t) = r
                                        c_empty(t) = c
                                        locs(t) = t
                                        t = t + 1
                                    end if
                                end do
                            end do
                        end do

                        ! Shuffle locations array so we can randomly place
                        ! new trees
                        call shuffle(locs)

                        ntrees: if (ntrenew >= 1) then

                            ! Renew trees
                            tree_renew: do irenew = 1, ntrenew

                                ! Determine species of new tree
                                rand = urand()
                                is = 1
                                do while (rand .gt. probt(is))
                                    is = is + 1
                                    if (is .gt. site%num_trees) then
                                        is = 1 + int(urand(0.0,                &
                                            real(site%num_trees)))
                                        rand = urand()
                                    endif
                                end do

                                ! Increment new tree counter
                                new_trees = new_trees + 1
                                is = tree_sp(is)

                                ! Decrement seedling bank of species
                                site%plots(ip)%seedling(is) =                  &
                                    max(0.0, site%plots(ip)%seedling(is) - 1.0)

                                ! Increment number of plants and
                                ! initialize new tree
                                it = it + 1
                                call initialize_tree(site%plots(ip)%trees(it), &
                                    site%species(is), is)

                                ! Grab the r and c values of that index
                                r = r_empty(locs(irenew))
                                c = c_empty(locs(irenew))

                                ! Set tree location to that value
                                site%plots(ip)%trees(it)%row = r
                                site%plots(ip)%trees(it)%col = c

                                ! Set this to filled in cells array
                                site%plots(ip)%cells(r,c) = 1

                                ! Get dbh value of new tree
                                ! (must be between 0.5 and 2.5 cm)
                                dbh = 1.5 + nrand(0.0, 1.0)
                                if (dbh >= 2.5) dbh = 2.5
                                if (dbh <= 0.5) dbh = 0.5
                                site%plots(ip)%trees(it)%diam_bht = dbh

                                ! Set clear branch bole height
                                site%plots(ip)%trees(it)%canopy_ht = 1.0

                                ! Get other characteristics
                                call forska_height(site%plots(ip)%trees(it))
                                call stem_shape(site%plots(ip)%trees(it))
                                call biomass_c(site%plots(ip)%trees(it))
                                call biomass_n(site%plots(ip)%trees(it))
                                call leaf_biomass_c(site%plots(ip)%trees(it))

                                ! Leaf biomass
                                leafbm = site%plots(ip)%trees(it)%leaf_bm

                                ! Tree height
                                ht = max(int(site%plots(ip)%trees(it)%forska_ht), 1)

                                !calculate shading effect on tree
                                if (site%plots(ip)%trees(it)%conifer) then
                                    shade = light_rsp(site%species(is),        &
                                        site%plots(ip)%con_light(ht))
                                else
                                    shade = light_rsp(site%species(is),        &
                                        site%plots(ip)%dec_light(ht))
                                end if

                                !calculate environmental stressors
                                call env_stress(site%plots(ip)%trees(it),      &
                                    shade, site%plots(ip)%fc_gdd(is),          &
                                    site%plots(ip)%fc_drought(is),             &
                                    site%plots(ip)%fc_perm(is),                &
                                    envstress, site%plots(ip)%fc_nutr(is))

                                ! Add leaf litter to soil and update NPP and
                                ! N used
                                if (site%species(is)%conifer) then

                                    NPP = NPP + leafbm*                        &
                                    (1.0 - CON_LEAF_RATIO) +                   &
                                    site%plots(ip)%trees(it)%biomC +           &
                                    site%plots(ip)%trees(it)%rootC

                                    N_used = N_used + leafbm/CON_LEAF_C_N +    &
                                        site%plots(ip)%trees(it)%biomN

                                    lc = site%species(is)%litter_class
                                    site%plots(ip)%soil%litter(lc) =           &
                                        site%plots(ip)%soil%litter(lc) +       &
                                        leafbm*CON_LEAF_RATIO/B_TO_C
                                else
                                    NPP = NPP + leafbm +                       &
                                    site%plots(ip)%trees(it)%biomC +           &
                                    site%plots(ip)%trees(it)%rootC

                                    N_used = N_used +                          &
                                        site%plots(ip)%trees(it)%biomN +       &
                                        leafbm/DEC_LEAF_C_N

                                    lc = site%species(is)%litter_class
                                    site%plots(ip)%soil%litter(lc) =           &
                                        site%plots(ip)%soil%litter(lc) +       &
                                        leafbm/B_TO_C
                                end if

                            end do tree_renew

                        end if ntrees

                    end if nempty

                end if renew

                ! Set new number of trees
                site%plots(ip)%numtrees = it

                ! Decrease seedling bank for survivability
                ! Also convert to #/m2
                do is = 1, num_species
                    site%plots(ip)%seedling(is) = site%plots(ip)%seedling(is)* &
                        site%species(is)%seedling_surv/plotsize
                end do

            end if navail

            ! Update site and soil variables
            N_used = N_used*HEC_TO_M2/plotsize
            site%plots(ip)%NPP = site%plots(ip)%NPP + NPP*HEC_TO_M2/plotsize

            ! Convert to tonnes/ha
            do i = 1, 18
                site%plots(ip)%soil%litter(i) =                                &
                    site%plots(ip)%soil%litter(i)*HEC_TO_M2/plotsize
            end do

            ! Deallocate locs array
            if (allocated(locs)) deallocate(locs)

        end do plot

    end subroutine Renewal

    !:.........................................................................:

end module Model
