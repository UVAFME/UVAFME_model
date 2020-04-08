module Site

  use Constants
  use Input
  use Species
  use Soil
  use Plot
  use Utilities
  use lists

  implicit none

!*******************************************************************************
  !
  !The Site module contains the definition of the SiteData type, which
  !holds attributes (physical and geographical) of a particular site,
  !as well as procedures directly operating on site variables
  !
!*******************************************************************************

	type SiteData
		type(PlotData),    dimension(:),   allocatable :: plots
		type(SpeciesData), dimension(:),   allocatable :: species
		real,              dimension(:),   allocatable :: fc_flood
		real,              dimension(:),   allocatable :: fc_drought
		real,              dimension(:),   allocatable :: fc_degday
		character(len = MAX_NLEN)                      :: region
		character(len = MAX_NLEN)                      :: site_name
		integer                                        :: site_id, runID
		integer                                        :: numplots
		integer                                        :: gcm_year
		real                                           :: latitude
		real                                           :: longitude
		real                                           :: elevation
		real                                           :: altitude, slope
		real                                           :: aspect
		real                                           :: heat_load
        real                                           :: pot_evap_day
		real                                           :: leaf_area_ind
		real, dimension(NTEMPS)                        :: temp_lapse_r
		real, dimension(NTEMPS)                        :: precip_lapse_r
		real, dimension(NTEMPS)                        :: tmin, tmax
		real, dimension(NTEMPS)                        :: precip
		real, dimension(NTEMPS)                        :: tmin_std
		real, dimension(NTEMPS)                        :: tmax_std
		real, dimension(NTEMPS)                        :: cld
		real, dimension(NTEMPS)                        :: cld_std
		real, dimension(NTEMPS)                        :: precip_std
		real, dimension(60)                            :: lai_array
		real                                           :: rain
		real                                           :: solar
		real                                           :: aridity, aridity_base
		real                                           :: grow_days
		real                                           :: e1, e2
		real                                           :: deg_days
		real                                           :: fire_prob
		real                                           :: wind_prob
		real                                           :: wd_ind
		real                                           :: alff, pmax
        real                                           :: pc_germ
	end type SiteData


contains

!-------------------------------------------------------------------------------
!  Methods
!-------------------------------------------------------------------------------

	!:...........................................................................:

	subroutine initialize_site(self, site_vals, species_data)
		!initializes site with input site variables
		!Author: Katherine Holcomb 2012, v. 1.0
		!Inputs/Outputs:
		!	self:          site instance
		!	species_data:  species data instance
		!Inputs:
		!	site_vals:     _sitelist.csv parameters

		class(SiteData),                   intent(inout) :: self
		type(SpeciesData),   dimension(:), intent(inout) :: species_data
		real,                dimension(:), intent(in)    :: site_vals

		character(len = 8),  dimension(:),   allocatable :: range_species_ids
		character(len = MAX_NLEN)                        :: sitename
		character(len = MAX_NLEN)                        :: siteregion
		real                                             :: gcm_year
		real                                             :: lat
		real                                             :: long
		real                                             :: elevation
		real                                             :: slope
		real                                             :: aspect
		real, dimension(NTEMPS)                          :: tmin, tmax
		real, dimension(NTEMPS)                          :: prcp
		real, dimension(NTEMPS)                          :: tmin_std
		real, dimension(NTEMPS)                          :: tmax_std
		real, dimension(NTEMPS)                          :: prcp_std
		real, dimension(NTEMPS)                          :: cld
		real, dimension(NTEMPS)                          :: cld_std
		real                                             :: altitude
		real                                             :: wind_level
		real                                             :: fire_level
		real                                             :: fire_prob
		real                                             :: wind_prob
		real                                             :: a_sat
		real                                             :: a_fc, wd
		real                                             :: o_depth, itxt
		integer                                          :: i_text

		!initialize properties from the site file
		self%site_id = int(site_vals(1))
		self%runID = int(site_vals(2))

		call read_site(self%site_id, sitename, siteregion, lat, long,          &
						elevation, slope, aspect, a_sat, a_fc, itxt,           &
						fire_prob, wind_prob, gcm_year, wd)

		self%site_name = sitename
		self%region    = siteregion
		self%latitude = lat
		self%longitude = long
		self%elevation = elevation
		self%slope = slope
		self%aspect = aspect
		self%fire_prob = fire_prob
		self%wind_prob = wind_prob
		self%gcm_year = int(gcm_year)
		i_text = int(itxt)
		self%wd_ind = wd

		!adjust according to global parameters set by user
		altitude = site_vals(3)
		fire_level = site_vals(4)
		wind_level = site_vals(5)
		o_depth = site_vals(6)

		if (fire_level .ne. rnvalid)  self%fire_prob = fire_level
		if (wind_level .ne. rnvalid)  self%wind_prob = wind_level

		!standard adjustments
		self%fire_prob = self%fire_prob/1000.0
		self%wind_prob = self%wind_prob/1000.0
		self%leaf_area_ind = 1.0
        self%pc_germ = .false.

		!read in climate data
		call read_climate(self%site_id, tmin, tmax, prcp)
		call read_extra_climate(self%site_id, cld)

		!adjust climate data for correct units
		self%tmin = tmin
		self%tmax = tmax
		self%precip = prcp*mm_to_cm
		self%cld = cld/10.0

		!adjust temperature and precipitation for altitude if necessary
		if (altitude .ne. rnvalid) then
			adjust_altitude = .true.
			self%altitude = altitude
			call adjustForAltitude(self)
		else
			self%altitude = self%elevation
		endif

		!read in climate stddev data
		if (self%site_id .ne. invalid) then
			if (use_climstd) then
				call read_climate_stds(self%site_id, tmin_std, tmax_std,       &
					prcp_std)
				self%tmin_std = tmin_std
				self%tmax_std = tmax_std
				self%precip_std = prcp_std*mm_to_cm
			endif

			if (use_cexstd) then
				call read_extra_climate_stds(self%site_id, cld_std)
				self%cld_std = cld_std/10.0
			end if
		endif

		!read in rangelist for site and add to site, also initialize plot-level
		  !soil information
		if (self%site_id .ne. invalid) then
			if  (use_rangelist) then
				call read_rangelist(self%site_id, range_species_ids)

				!add range information to site
				if (size(range_species_ids) /= 0) then
					call attach_species(self, species_data, a_sat, a_fc,       &
						o_depth, i_text, range_species_ids)
				else
					call attach_species(self, species_data, a_sat, a_fc,       &
						o_depth, i_text)
				endif

			else
				!all species present in this site and initialize plots
				call attach_species(self, species_data, a_sat, a_fc, o_depth,  &
					i_text)
			endif
		endif

	end subroutine initialize_site

	!:.........................................................................:

	subroutine attach_species(self, species_data, a_sat, a_fc, o_depth,        &
		i_text, range_species_ids)
		!attaches correct species list for current site, depending on
		  !rangelist, and initializes plots on site. Makes use of
		  !overloading = for Species type
		!Author: Katherine Holcomb 2012, v. 1.0
		!Inputs/Outputs:
		!	self:              site instance
		!	species_data:      species data instance
		!Inputs:
		!	a_sat:             mineral (A-layer) saturation capacity (volumetric)
		!	a_fc:              mineral (A-layer) field capacity (volumetric)
		!	o_depth:           organic soil layer depth (m)
		!	i_text: 	       soil texture (1: coarse-grained; 2: fine-grained)
		!	range_species_ids: species eleible for colonizaiton in site

		class(SiteData),                            intent(inout) :: self
		type(SpeciesData),            dimension(:), intent(inout) :: species_data
		character(len = *), optional, dimension(:), intent(in)    :: range_species_ids
		real,                                       intent(in)    :: a_sat, a_fc, o_depth
		integer,                                    intent(in)    :: i_text
		integer                                                   :: num_all_species
		integer                                                   :: num_site_species
		integer                                                   :: num_range_species
		integer                                                   :: n, nn

		if (present(range_species_ids)) then

			num_all_species = size(species_data)
			num_range_species = size(range_species_ids)
			num_site_species = count(range_species_ids .ne. 'NP')

			if (num_site_species == 0) then
				allocate(self%species(0))
			else
				do n = 1, num_all_species
					do nn = 1, num_range_species
						if (species_data(n)%unique_id .eq.                     &
							range_species_ids(nn)) then
							call append(self%species, species_data(n))
						endif
					enddo
				enddo
			endif

		else
			!no range list, all species in this site
			do n = 1, num_all_species
				call append(self%species, species_data(n))
			enddo

		endif

		!now that we have the species info, initialize plots
		self%numplots = numplots
		allocate(self%plots(self%numplots))
		do n = 1, self%numplots
			call initialize_plot(self%plots(n), self%species, maxcells,        &
				maxheight, a_sat, a_fc, o_depth, i_text)
		enddo

	end subroutine attach_species

	!:.........................................................................:

	subroutine adjustForAltitude(self)
		!adjusts precip and temperature data for altitude
		!Author: Katherine Holcomb 2012, v. 1.0
		!Inputs/Outputs:
		!	self:              site instance

		class(SiteData), intent(inout) ::  self

		integer                        :: z

		if (self%altitude .ne. rnvalid) then
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
		!deletes attributes from a site
		!Author: Katherine Holcomb 2012, v. 1.0
		!Inputs/Outputs:
		!	self:              site instance

		class(SiteData), intent(inout) :: self
		integer                        :: n

		do n = 1, numplots
			call delete_plot(self%plots(n))
		enddo

		if (allocated(self%plots)) deallocate(self%plots)
		if (allocated(self%species)) deallocate(self%species)
		if (allocated(self%fc_flood)) deallocate(self%fc_flood)
		if (allocated(self%fc_drought)) deallocate(self%fc_drought)
		if (allocated(self%fc_degday)) deallocate(self%fc_degday)

	end subroutine delete_site

	!:.........................................................................:

	subroutine write_site_csv(self, site_unit)
	use csv_file
		!helper function for writing climate data to the Climate.csv file
		!Author: Katherine Holcomb 2012, v. 1.0
		!Inputs/Outputs:
		!	self:       site instance
		!	site_unit:  Climate.csv file unit number


		class(SiteData), intent(in)  :: self
		integer,         intent(in)  :: site_unit

		real, dimension(numplots)    :: aet
		real, dimension(numplots)    :: drydays_up, drydays_b
		real, dimension(numplots)    :: flooddays, active, org, availn
		real, dimension(numplots)    :: wldays, wilt_days

		real                         :: aet_mn, aet_sd
		real                         :: drydays_up_mn
		real                         :: drydays_up_sd, drydays_b_mn
		real                         :: drydays_b_sd, flooddays_mn
		real                         :: flooddays_sd, active_mn
		real                         :: active_sd, org_mn, org_sd
		real                         :: availn_mn, availn_sd
		real                         :: wldays_mn, wldays_sd
        real                         :: wilt_days_mn, wilt_days_sd
		integer                      :: ip

		do ip = 1, numplots

			aet(ip) = self%plots(ip)%act_evap_day
			drydays_up(ip) = self%plots(ip)%dry_days_upper_layer
			drydays_b(ip) = self%plots(ip)%dry_days_base_layer
			flooddays(ip) = self%plots(ip)%flood_days
			active(ip) = self%plots(ip)%soil%active
			org(ip) = self%plots(ip)%soil%O_depth
			availn(ip) = self%plots(ip)%soil%avail_N
			wldays(ip) = self%plots(ip)%wl_days
            wilt_days(ip) = self%plots(ip)%wilt_days

		end do

		call stddev(aet, aet_mn, aet_sd, rnvalid)
		call stddev(drydays_up, drydays_up_mn, drydays_up_sd, rnvalid)
		call stddev(drydays_b, drydays_b_mn, drydays_b_sd, rnvalid)
		call stddev(flooddays, flooddays_mn, flooddays_sd, rnvalid)
		call stddev(active, active_mn, active_sd, rnvalid)
		call stddev(org, org_mn, org_sd, rnvalid)
		call stddev(availn, availn_mn, availn_sd, rnvalid)
		call stddev(wldays, wldays_mn, wldays_sd, rnvalid)
        call stddev(wilt_days, wilt_days_mn, wilt_days_sd, rnvalid)

		call csv_write(site_unit, self%rain, .false.)
		call csv_write(site_unit, self%pot_evap_day, .false.)
		call csv_write(site_unit, self%solar, .false.)
		call csv_write(site_unit, self%heat_load, .false.)
		call csv_write(site_unit, active_mn, .false.)
		call csv_write(site_unit, org_mn, .false.)
		call csv_write(site_unit, availn_mn, .false.)
		call csv_write(site_unit, aet_mn, .false.)
		call csv_write(site_unit, self%grow_days, .false.)
        call csv_write(site_unit, self%pc_germ, .false.)
		call csv_write(site_unit, self%deg_days, .false.)
		call csv_write(site_unit, drydays_up_mn, .false.)
		call csv_write(site_unit, drydays_b_mn, .false.)
        call csv_write(site_unit, wilt_days_mn, .false.)
		call csv_write(site_unit, flooddays_mn, .false.)
		call csv_write(site_unit, wldays_mn, .true.)

	end subroutine write_site_csv

	!:.........................................................................:

   end module Site
