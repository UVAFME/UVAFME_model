module Site

!*******************************************************************************
  !
  ! The Site module contains the definition of the SiteData type, which
  ! holds attributes (physical and geographical) of a particular site,
  ! as well as procedures directly operating on site variables.
  !
!*******************************************************************************

    use Constants
    use Input
    use Species
    use Soil
    use Plot
    use Utilities
    use lists

    implicit none

    ! Define SiteData types
    type SiteData
        type(PlotData),    dimension(:),     allocatable :: plots           ! Array of plot objects
        type(SpeciesData), dimension(:),     allocatable :: species         ! Array of species objects
        real,              dimension(:),     allocatable :: lai_array       ! LAI with height (m2/m2)
        real,              dimension(NTEMPS)             :: temp_lapse_r    ! Temperature lapse rate (degC/km)
        real,              dimension(NTEMPS)             :: precip_lapse_r  ! Precipitation lapse rate (mm/km)
        real,              dimension(NTEMPS)             :: tmin            ! Mean monthly minimum temperature (degC)
        real,              dimension(NTEMPS)             :: tmax            ! Mean monthly maximum temperature (degC)
        real,              dimension(NTEMPS)             :: precip          ! Mean monthly precipitation (cm)
        real,              dimension(NTEMPS)             :: tmin_std        ! SD of monthly minimum temperature (degC)
        real,              dimension(NTEMPS)             :: tmax_std        ! SD of monthly maximum temperature (degC)
        real,              dimension(NTEMPS)             :: precip_std      ! SD of monthly precipitation (cm)
        real,              dimension(NTEMPS)             :: cld             ! Mean monthly cloudiness (tenths of sky covered)
        real,              dimension(NTEMPS)             :: cld_std         ! SD of monthly cloudiness (tenths of sky covered)
        real,              dimension(NTEMPS)             :: accum_precip    ! Accumulated precipitation for linear climate change (cm)
        character(len = MAX_NLEN)                        :: region          ! Site region
        character(len = MAX_NLEN)                        :: site_name       ! Site name
        integer                                          :: site_id         ! Site ID
        integer                                          :: runID           ! Run ID
        integer                                          :: numplots        ! Number of plots
        integer                                          :: gcm_year        ! Year to start climate change
        integer                                          :: num_trees       ! Number of tree species
        real                                             :: latitude        ! Site latitude (degrees)
        real                                             :: longitude       ! Site longitude (degrees)
        real                                             :: elevation       ! Site elevation (m)
        real                                             :: altitude        ! Site altitude (for altitude climate adjustments) (m)
        real                                             :: slope           ! Site slope (degrees)
        real                                             :: aspect          ! Site aspect (degrees)
        real                                             :: leaf_area_ind   ! Leaf area index (m2/m2)
        real                                             :: rain            ! Annual rainfall (cm)
        real                                             :: pot_evap_day    ! Annual potential evapotranspiration (cm)
        real                                             :: solar           ! Annual solar radiation (cal/cm2/day)
        real                                             :: grow_days       ! Growing season (>5degC) length (days)
        real                                             :: deg_days        ! Growing degree-days (>5degC)
        real                                             :: wind_prob       ! Windthrow probability (0-1)
        real                                             :: fire_prob       ! Fire probability (0-1)
        real                                             :: accum_tmax      ! Accumulated minimum temperature for linear climate change (degC)
        real                                             :: accum_tmin      ! Accumulated maximum temperature for linear climate change (degC)
        real                                             :: aridity         ! Current year's index
        real                                             :: aridity_base    ! "Base" aridity index (calculated from first 10 years of simulation)
        real                                             :: alff            ! Available light on the forst floor (0-1)
    end type SiteData

contains

    !:...........................................................................:

    subroutine initialize_site(self, site_vals, species_data, species_ids, sndx)
        !
        !  Inititalizes a site with input site parameters and starting values.
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb          Original Code
        !    01/26/21     A. C. Foster        Separated attach_species and
        !                                       initialize_plot calls and skip
        !                                       site if bad rangelist data
        !

        ! Data dictionary: calling arguments
        class(SiteData),                    intent(inout) :: self         ! Site object
        type(SpeciesData),  dimension(:),   intent(inout) :: species_data ! Array of species objects
        integer,            dimension(:,:), intent(inout) :: species_ids  ! Array of species present at each site
        real,               dimension(:),   intent(in)    :: site_vals    ! Site-level runtime parameters
        integer,                            intent(in)    :: sndx         ! Site index in site array

        ! Data dictionary: local variables
        character(len = 8), dimension(:),    allocatable :: range_species_ids ! Unique IDs of species at this site
        real,              dimension(NTEMPS)             :: tmin              ! Mean monthly minimum temperature (degC)
        real,              dimension(NTEMPS)             :: tmax              ! Mean monthly minimum temperature (degC)
        real,              dimension(NTEMPS)             :: prcp              ! Mean monthly precipitation (mm)
        real,              dimension(NTEMPS)             :: tmin_std          ! SD of monthly minimum temperature (degC)
        real,              dimension(NTEMPS)             :: tmax_std          ! SD of monthly maximum temperature (degC)
        real,              dimension(NTEMPS)             :: prcp_std          ! SD of monthly precipitation (mm)
        real,              dimension(NTEMPS)             :: cld               ! Mean monthly cloud cover (tenths of sky covered)
        real,              dimension(NTEMPS)             :: cld_std           ! SD monthly cloud cover (tenths of sky covered)
        character(len = MAX_NLEN)                        :: sitename          ! Name of site
        character(len = MAX_NLEN)                        :: siteregion        ! Site region
        character(len = MAX_CHAR)                        :: message           ! Warning message
        real                                             :: gcm_year          ! Year to start climate change
        real                                             :: lat               ! Latitude (degrees)
        real                                             :: long              ! Longitude (degrees)
        real                                             :: elevation         ! Elevation (m)
        real                                             :: slope             ! Slope (degrees)
        real                                             :: aspect            ! Aspect (degrees)
        real                                             :: altitude          ! Site altitude (for altitude climate adjustments) (m)
        real                                             :: wind_prob         ! Number of windthrow events in 1000 years
        real                                             :: fire_prob         ! Number of fires in 1000 years
        real                                             :: a_sat             ! A-layer saturation capacity (volumetric)
        real                                             :: a_fc              ! A-layer field capacity (volumetric)
        real                                             :: a_pwp             ! A-layer permanent wilting point (volumetric)
        real                                             :: hum_input         ! Initial humus amount (t/ha)
        real                                             :: o_sat             ! Organic layer saturation capacity (volumetric)
        real                                             :: o_fc              ! Organic layer field capacity (volumetric)
        real                                             :: o_pwp             ! Organic permanent wilting point (volumetric)
        real                                             :: o_bd              ! Organic layer bulk density (kg/m3)
        real                                             :: a_bd              ! A-layer bulk density (kg/m3)
        real                                             :: flow              ! Moisture input from overland flow (mm)
        integer                                          :: itxt              ! Soil texture (0: very coarse; 1: coarse; 2: fine)
        integer                                          :: ip                ! Looping index

        ! Initialize properties from the sitelist file
        self%site_id = int(site_vals(1))
        self%runID = int(site_vals(2))
        altitude = site_vals(3)

        ! Read in site file
        call read_site(self%site_id, sitename, siteregion, lat, long,          &
            elevation, slope, aspect, a_sat, a_fc, a_pwp, o_sat, o_fc, o_pwp,  &
            o_bd, a_bd, itxt, hum_input, fire_prob, wind_prob, gcm_year)

        ! Initialize values
        self%site_name = sitename
        self%region = siteregion
        self%latitude = lat
        self%longitude = long
        self%elevation = elevation
        self%slope = slope
        self%aspect = aspect
        self%wind_prob = wind_prob
        self%gcm_year = int(gcm_year)

        ! Standard adjustments and settings
        self%wind_prob = self%wind_prob/1000.0
        self%fire_prob = self%fire_prob/1000.0
        self%leaf_area_ind = 1.0
        allocate(self%lai_array(maxheight))
        self%lai_array = 0.0

        ! Read in climate data
        call read_climate(self%site_id, tmin, tmax, prcp, cld)

        ! Adjust climate data for correct units
        self%tmin = tmin
        self%tmax = tmax
        self%precip = prcp*MM_TO_CM
        self%cld = cld/10.0
        self%accum_tmax = 0.0
        self%accum_tmin = 0.0
        self%accum_precip = 0.0

        ! Adjust temperature and precipitation for altitude if necessary
        if (altitude .ne. RNVALID) then
            adjust_altitude = .true.
            self%altitude = altitude
            call adjustForAltitude(self)
        else
            adjust_altitude = .false.
            self%altitude = self%elevation
        endif

        ! Read in climate stddev data
        if (self%site_id .ne. INVALID) then
            if (use_climstd) then
                call read_climate_stds(self%site_id, tmin_std, tmax_std,       &
                    prcp_std, cld_std)
                self%tmin_std = tmin_std
                self%tmax_std = tmax_std
                self%precip_std = prcp_std*MM_TO_CM
                self%cld_std = cld_std/10.0
            end if
        endif

        ! Read in rangelist for site and add to site, also initialize plot-level
        ! soil information
        if (self%site_id .ne. INVALID) then
            if  (use_rangelist) then
                ! We are using a rangelist - read in
                call read_rangelist(self%site_id, range_species_ids)

                ! Add range information to site
                if (size(range_species_ids) /= 0) then
                    ! Read was successfull - attach some species
                    call attach_species(self, species_data, species_ids, sndx, &
                        range_species_ids)
                else
                    ! Read was not successful or something weird happened
                    ! Tell user
                    write(message, '(A,A,I9)') "Bad read of rangelist for ",   &
                        " site ", self%site_id
                    call warning(message)
                    self%site_id = INVALID
                endif
            else
                ! Not using rangelist
                ! All species present in this site
                call attach_species(self, species_data, species_ids, sndx)
            endif
        endif

        ! Now that we have the species info, initialize plots
        if (self%site_id .ne. INVALID) then
            self%numplots = numplots
            allocate(self%plots(self%numplots))
            do ip = 1, self%numplots
                call initialize_plot(self%plots(ip), size(self%species),       &
                    a_sat, a_fc, a_pwp, o_sat, o_fc, o_pwp, o_bd, a_bd, itxt,  &
                    hum_input)
            enddo
        end if

    end subroutine initialize_site

    !:.........................................................................:

    subroutine attach_species(self, species_data, species_ids, sndx,           &
         range_species_ids)
         !
         !  Attaches correct species list for current site, depending on the
         !  rangelist for the site.
         !
         !  Record of revisions:
         !      Date       Programmer          Description of change
         !      ====       ==========          =====================
         !    01/01/12     K. Holcomb          Original Code
         !    10/10/19     A. C. Foster        Updated for seed rain so that we can
         !                                     have 'native' and 'non-native'
         !                                     species
         !    01/26/21     A. C. Foster        Moved initialize plot outside of this
         !                                       subroutine
         !

        ! Data dictionary: calling arguments
        class(SiteData),                              intent(inout) :: self              ! Site object
        type(SpeciesData),            dimension(:),   intent(inout) :: species_data      ! Array of species objects
        character(len = *), optional, dimension(:),   intent(in)    :: range_species_ids ! Unique IDs for species present at this site
        integer,                      dimension(:,:), intent(inout) :: species_ids       ! Site IDs x species IDs
        integer,                                      intent(in)    :: sndx              ! Index of site in site array

        ! Data dictionary: local variables
        integer :: num_all_species   ! Number of total species in input file
        integer :: num_site_species  ! Number of total species in site
        integer :: num_range_species ! Number of total species in rangelist file
        integer :: n, nn, t          ! Looping indices

        if (present(range_species_ids)) then

            ! Total number of species in specieslist
            num_all_species = size(species_data)

            ! Number of species in rangelist file
            num_range_species = size(range_species_ids)

            ! Number of species native to this site
            num_site_species = count(range_species_ids .ne. 'NP')

            ! Loop through and add native species to species object
            if (num_site_species == 0) then
                allocate(self%species(0)) ! Nothing to allocate, no native species
            else
                ! Loop through and add species data for each native species
                ! Also add species id to species_id array for specific site
                t = 1
                do n = 1, num_all_species
                    do nn = 1, num_range_species
                        if (species_data(n)%unique_id .eq.                     &
                            range_species_ids(nn)) then
                            call append(self%species, species_data(n))
                            species_ids(sndx, t) = n
                            t = t + 1
                        endif
                    enddo
                enddo
            endif
        else
            t = 1
            ! No range list, all species in this site
            do n = 1, num_all_species
                call append(self%species, species_data(n))
                species_ids(sndx, t) = n
                t = t + 1
            enddo

        endif

        ! Count number of tree species
        self%num_trees = size(self%species)

    end subroutine attach_species

    !:.........................................................................:

    subroutine adjustForAltitude(self)
        !
        !  Adjusts precipitaiton and temperature data for altitude
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb          Original Code
        !

        ! Data dictionary: calling arguments
        class(SiteData), intent(inout) :: self ! Site object

        ! Data dictionary: local variables
        integer :: z ! Looping index

        ! Loop through and adjust each month
        ! Lapse rates are in degC/km and mm/km, elevation is in m
        if (self%altitude .ne. RNVALID .and. adjust_altitude) then
            do  z = 1, 12
                self%tmax(z) = self%tmax(z) -                                  &
                    (self%altitude - self%elevation)*self%temp_lapse_r(z)*0.01
                self%tmin(z) = self%tmin(z) -                                  &
                    (self%altitude - self%elevation)*self%temp_lapse_r(z)*0.01
                self%precip(z) = (max(self%precip(z) +                         &
                    (self%altitude - self%elevation)*self%precip_lapse_r(z)*   &
                    0.001, 0.0))
            end do
        endif

    end subroutine adjustForAltitude

    !:.........................................................................:

    subroutine delete_site(self)
        !
        !  Frees memory associated with a site
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb          Original Code
        !

        ! Data dictionary: calling arguments
        class(SiteData), intent(inout) :: self ! Site object

        ! Data dictionary: local variables
        integer :: ip ! Looping index

        ! Free memory from plots
        do ip = 1, numplots
            call delete_plot(self%plots(ip))
        enddo

        ! Deallocate site arrays
        if (allocated(self%plots)) deallocate(self%plots)
        if (allocated(self%species)) deallocate(self%species)
        if (allocated(self%lai_array)) deallocate(self%lai_array)

    end subroutine delete_site

    !:.........................................................................:

    subroutine write_site_csv(self, site_unit)
        !
        !  Helper function for writing climate & site data to the output climate
        !   file
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb          Original Code
        !
        use csv_file

        ! Data dictionary: calling arguments
        class(SiteData), intent(in) :: self      ! Site object
        integer,         intent(in) :: site_unit ! File unit for output climate

        ! Data dictionary: local variables
        real, dimension(numplots) :: aet           ! AET (cm)
        real, dimension(numplots) :: drydays       ! Drought index
        real, dimension(numplots) :: flooddays     ! Flooding index
        real, dimension(numplots) :: active        ! Active layer depth (cm)
        real, dimension(numplots) :: org           ! Organic layer depth (cm)
        real, dimension(numplots) :: availn        ! Plant-available N (kgN/ha)
        real, dimension(numplots) :: wilt_days     ! Wilting-point index
        real, dimension(numplots) :: saw0_ByFC     ! A-layer moisture scaled by field capacity
        real, dimension(numplots) :: aow0_ByMin    ! Organic-layer moisture scaled by wilting point
        real, dimension(numplots) :: saw0_BySAT    ! A-layer moisture scaled by saturation capacity
        real                      :: aet_mn        ! Average AET (cm)
        real                      :: aet_sd        ! SD of AET (cm)
        real                      :: drydays_mn    ! Average drought index
        real                      :: drydays_sd    ! SD of drought index
        real                      :: flooddays_mn  ! Average Flooding index
        real                      :: flooddays_sd  ! SD of flooding index
        real                      :: active_mn     ! Average active layer depth (cm)
        real                      :: active_sd     ! SD of active layer depth (cm)
        real                      :: org_mn        ! Average organic layer depth (cm)
        real                      :: org_sd        ! SD of organic layer depth (cm)
        real                      :: availn_mn     ! Average plant-available N (kgN/ha)
        real                      :: availn_sd     ! SD of plant-available N (tkgN/ha)
        real                      :: wilt_days_mn  ! Average wilting point index
        real                      :: wilt_days_sd  ! SD of wilting point index
        real                      :: saw0_ByFC_mn  ! Average A-layer moisture scaled by field capacity
        real                      :: saw0_ByFC_sd  ! SD of A-layer moisture scaled by field capacity
        real                      :: saw0_BySAT_mn ! Average A-layer moisture scaled by saturation capacity
        real                      :: saw0_BySAT_sd ! SD of A-layer moisture scaled by saturation capacity
        real                      :: aow0_ByMin_mn ! Average organic layer moisture scaled by wilting point
        real                      :: aow0_ByMin_sd ! SD of organic layer moisture scaled by wilting point
        integer                   :: ip            ! Looping index

        do ip = 1, numplots

            ! Read in relevant variables
            aet(ip) = self%plots(ip)%act_evap_day
            drydays(ip) = self%plots(ip)%dry_days
            flooddays(ip) = self%plots(ip)%flood_days
            active(ip) = self%plots(ip)%soil%active*M_TO_CM
            org(ip) = self%plots(ip)%soil%O_depth*M_TO_CM
            availn(ip) = self%plots(ip)%soil%avail_N*T_TO_KG
            wilt_days(ip) = self%plots(ip)%wilt_days
            saw0_ByFC(ip) = self%plots(ip)%saw0_ByFC
            saw0_BySAT(ip) = self%plots(ip)%saw0_BySAT
            aow0_ByMin(ip) = self%plots(ip)%aow0_ByMin
        end do

        ! Get mean and sd
        call stddev(aet, aet_mn, aet_sd, RNVALID)
        call stddev(drydays, drydays_mn, drydays_sd, RNVALID)
        call stddev(flooddays, flooddays_mn, flooddays_sd, RNVALID)
        call stddev(active, active_mn, active_sd, RNVALID)
        call stddev(org, org_mn, org_sd, RNVALID)
        call stddev(availn, availn_mn, availn_sd, RNVALID)
        call stddev(wilt_days, wilt_days_mn, wilt_days_sd, RNVALID)
        call stddev(saw0_ByFC, saw0_ByFC_mn, saw0_ByFC_sd, RNVALID)
        call stddev(saw0_BySAT, saw0_BySAT_mn, saw0_BySAT_sd, RNVALID)
        call stddev(aow0_ByMin, aow0_ByMin_mn, aow0_ByMin_sd, RNVALID)

        ! Write to file
        call csv_write(site_unit, self%rain, .false.)
        call csv_write(site_unit, self%pot_evap_day, .false.)
        call csv_write(site_unit, self%solar, .false.)
        call csv_write(site_unit, active_mn, .false.)
        call csv_write(site_unit, org_mn, .false.)
        call csv_write(site_unit, availn_mn, .false.)
        call csv_write(site_unit, aet_mn, .false.)
        call csv_write(site_unit, self%grow_days, .false.)
        call csv_write(site_unit, self%deg_days, .false.)
        call csv_write(site_unit, drydays_mn, .false.)
        call csv_write(site_unit, saw0_ByFC_mn, .false.)
        call csv_write(site_unit, saw0_BySAT_mn, .false.)
        call csv_write(site_unit, aow0_ByMin_mn, .false.)
        call csv_write(site_unit, wilt_days_mn, .false.)
        call csv_write(site_unit, flooddays_mn, .true.)

    end subroutine write_site_csv

    !:.........................................................................:

   end module Site
