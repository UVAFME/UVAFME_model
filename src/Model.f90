module Model

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

!*******************************************************************************
  !
  !This contains the five main subroutines of UVAFME:
  !
  !BioGeoClimate: computes daily and yearly site- and plot-level weather and
  !               soil dynamics
  !
  !Canopy:        computes the plot-level LAI and light availability
  !
  !Growth:        computes annual tree growth and branch thinning
  !
  !Mortality:     determines which trees die and adds their components to the
  !               soil
  !
  !Renewal        updates the seed and seedling banks for each species and
  !               regenerates new trees
  !
!*******************************************************************************

    real,    parameter                         :: growth_min = 0.03
    integer, parameter                         :: mcount = 2

contains

	!:.........................................................................:

	subroutine BioGeoClimate(site, year)
		!computes daily data and annual sums of weather and soil characteristics
		!Inputs/Outputs:
		!	site: site instance
		!Inputs:
		!	year: year of simulation

		integer,                         intent(in)    :: year
		type(SiteData),                  intent(inout) :: site

		integer                                        :: gcm_year
		real,    dimension(NTEMPS)                     :: tmin, tmax
		real,    dimension(NTEMPS)                     :: prcp
		real,    dimension(NTEMPS)                     :: tmean, cld
		real,    dimension(NTEMPS)                     :: tmptmin
		real,    dimension(NTEMPS)                     :: tmptmax
		real,    dimension(NTEMPS)                     :: tmpprec
		real,    dimension(NTEMPS)                     :: tmpcld
		real,    dimension(days_per_year)              :: daytemp
		real,    dimension(days_per_year)              :: daytemp_min
		real,    dimension(days_per_year)              :: daytemp_max
		real,    dimension(days_per_year)              :: daycld
		real,    dimension(days_per_year)              :: dayprecip
		real,  	 dimension(days_per_year)              :: sun, st
		real,    dimension(days_per_year)              :: exrad
		real,    dimension(days_per_year)              :: pot_ev_day
		real,    dimension(12, 2)                      :: tdd, fdd
		integer, dimension(12)                         :: modays

		real                                           :: rain, rain_n
		real                                           :: temp_f, prcp_f
		real                                           :: cld_f
		real                                           :: temp_max
		real                                           :: temp_min, daytemp_mem
		real                                           :: tmean_max
		real                                           :: n_avail
		real                                           :: pet
		real                                           :: aet, aet_mm
		real                                           :: growdays
		real                                           :: soildays, wpdays
		real                                           :: drydays_upper
		real                                           :: drydays_base
		real                                           :: flooddays, wl_days
		real                                           :: degday
		real                                           :: outwater
		real                                           :: tot_sun
		real                                           :: tot_st, cfs
		real                                           :: cum_tdd
		real                                           :: act_ev_day
		real                                           :: tcum, fcum
		real                                           :: amlt, xmlt
		real                                           :: xfrz, zh
		real                                           :: alff, pc_germ
		real                                           :: aow0_ScaledByMax
		real                                           :: aow0_ScaledByMin
		real                                           :: saw0_ScaledByFC
		real                                           :: saw0_ScaledByWP
		real                                           :: saw0_ScaledBySAT

		!used to temporarily hold accumulated climate variables
		real                                           :: tmpstep1
		real                                           :: tmpstep2
		real                                           :: tmp

		!iterator integers and other indexing variables
		integer                                        :: i, j, m, ip
		integer                                  	   :: l
		integer                                        :: warmest_month, hrise

		!parameters for degdays, flooddays, and drydays
		real, parameter                                :: min_grow_temp = 5.0
		real, parameter                                :: max_dry_parm = 1.0001
		real, parameter                                :: min_flood_parm = 2.0

		!last julian day in each month
		data modays /31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365/

		save

		rain = 0.0
		rain_n = 0.0
		tmean_max = 0.0

		!check for and implement climate chnage
		!the user is expected to input decr_by values as positive
		if (linear_cc) then
			if (year .ge. site%gcm_year .AND. year .le.                        &
				(site%gcm_year + gcm_duration)) then
				accumulated_tmin = accumulated_tmin + tmin_change
				accumulated_tmax = accumulated_tmax + tmax_change
				do m = 1, 12
					tmpstep1 = site%precip(m) + accumulated_precip(m)
					tmpstep2 = tmpstep1 * precip_change
					accumulated_precip(m) = accumulated_precip(m) + tmpstep2
				end do
			endif
		else if (use_gcm) then
			!begin_change_year is site-specific so that different sites
			!can start climate change forcing at different times
			gcm_year = start_gcm + year - site%gcm_year
			if ( gcm_year .ge. start_gcm .and. gcm_year .le. end_gcm ) then
				call read_gcm_climate(site%site_id, gcm_year, start_gcm, tmin, &
                    tmax, prcp)
				site%tmin = tmin
				site%tmax = tmax
				site%precip = prcp*mm_to_cm
				call adjustForAltitude(site)
			endif
		endif

		!generate current year's weather from distributions of input climate
        !data, adjust for linear cc if needed
 		do i = 1, NTEMPS
 		  if (linear_cc) then
 				tmptmin(i) = site%tmin(i)   + accumulated_tmin
 				tmptmax(i) = site%tmax(i)   + accumulated_tmax
 				tmpprec(i) = site%precip(i) + accumulated_precip(i)
 				tmpcld(i)  = site%cld(i)
 			else
 				tmptmin(i) = site%tmin(i)
 				tmptmax(i) = site%tmax(i)
 				tmpprec(i) = site%precip(i)
 				tmpcld(i)  = site%cld(i)
 			endif

			!calculate climate fluctuations
 			temp_f = clim_nrand(0.0, 1.0)
 			prcp_f = clim_nrand(0.0, 1.0)
 			cld_f = clim_nrand(0.0, 1.0)

 			temp_f = max(-1.0, min(temp_f, 1.0))
 			cld_f = max(-1.0, min(cld_f, 1.0))

 			!forest cover can increase rainfall by maximum 15%
 			prcp_f = max(-0.5, min(prcp_f, 0.5))

 			!adjust monthly climate vars with std and random numbers
 			if	(use_gcm .and. (year .ge. site%gcm_year .and. year .le.        &
					(site%gcm_year + gcm_duration))) then

				!just use input values for that month
				tmin(i) = tmptmin(i)
				tmax(i) = tmptmax(i)
				prcp(i) = max(tmpprec(i), 0.0)
				cld(i) = max(tmpcld(i) + cld_f*site%cld_std(i), 0.0)
 			else
 				tmin(i) = tmptmin(i) + temp_f*site%tmin_std(i)
 				tmax(i)= tmptmax(i) + temp_f*site%tmax_std(i)

				!prcp and cld can't be less than 0.0
 				cld(i) = max(tmpcld(i) + cld_f*site%cld_std(i), 0.0)
 				prcp(i) = max(tmpprec(i) + prcp_f*site%precip_std(i), 0.0)
 			end if

 			!accumulate rainfall and N deposition
 			rain = rain + prcp(i)
 			rain_n = rain_n + prcp(i)*prcp_n

 			!get mean monthly temperature for warmest month calculation
 			tmean(i) = (site%tmin(i) + site%tmax(i))/2
 		end do

 		!find warmest month of this year
 		do i = 1, NTEMPS
 			tmean_max = max(tmean_max, tmean(i))
 			if (tmean_max .eq. tmean(i)) warmest_month = i
 		end do

		!get tmax and tmin of warmest month
 		temp_max = site%tmax(warmest_month)
 		temp_min = site%tmin(warmest_month)

 		!calculate e2 and e1 (used for PET calculation)
 		site%e2 = 33.8639*((0.00738*temp_max + 0.8072)**8.0 -                  &
            0.000019*abs(1.8*temp_max + 48.0) + 0.001316)

 		site%e1 = 33.8639*((0.00738*temp_min + 0.8072)**8.0 -                  &
 			0.000019*abs(1.8*temp_min + 48.0) + 0.001316)

 		!convert monthly weather data into daily weather data
 		call cov365(tmin, daytemp_min)
 		call cov365(tmax, daytemp_max)
 		call cov365a(prcp, dayprecip)
 		call cov365(cld, daycld)

 		!initialize pet and solar radiation accumulators
 		pet = 0.0
 		tot_sun = 0.0
 		tot_st = 0.0

 		!calculate mean daily temperature, solar radiation, and PET
        site%pc_germ = 0.0
		do i = 1, days_per_year
 			!mean daily temperature (degC)
 			daytemp(i) = 0.5*(daytemp_min(i) + daytemp_max(i))

            if (daytemp(i) > 5.0 .and. i .le. 200) then

                daytemp_mem = (daytemp(i-2) + daytemp(i-1) + daytemp(i))/3.0

                if (daytemp_mem <= 15.0) then
                    pc_germ = min(1.0, max(0.0, daytemp_mem*0.1 - 0.5))
                else
                    pc_germ = 1.0
                endif

                site%pc_germ = max(site%pc_germ, pc_germ)
            end if


 			!solar radiation (cal/cm2/day)
 			call ex_rad(i, site%latitude, site%slope, site%aspect, daycld(i),  &
						exrad(i), sun(i), st(i), hrise)

 			!accumulate surface and horizontal surface radiation
 			tot_sun = tot_sun + sun(i) !actual surface
 			tot_st = tot_st + st(i) !horizontal surface

			!PET (cm)
 			pot_ev_day(i) = jensenhaise(daytemp(i), sun(i), site%altitude,     &
								site%e2, site%e1)

            !accumulate PET (cm)
			pet = pet + pot_ev_day(i)

 		end do

 		!calculate fraction sunlight attenuated through atm.
 		cfs = tot_sun/tot_st

 		!calculate freezing and thawing degree days for permafrost equation
 		tdd = 0.0
 		fdd = 0.0
 		m = 1
 		do j = 1, days_per_year
			if (j .gt. modays(m)) m = m + 1

			if ((tmean(m) .gt. 0.001) .and. (daytemp(j) .gt. 0.001)) then
 				tdd(m, 1) = tdd(m, 1) + daytemp(j)
 			end if
 			if ((tmean(m) .le. 0.001) .and. (daytemp(j) .le. 0.001)) then
 				fdd(m, 1) = fdd(m, 1) + abs(daytemp(j))
 			end if

 		end do

 		!calculate cumulative freezing and thawing degree days
 		tcum = 0.0
 		fcum = 0.0
 		do m = 12, 1, -1
 			if (fdd(m, 1) .gt. 0.0) fcum = fcum + fdd(m, 1)
 			if (fdd(m, 1) .eq. 0.0) GO TO 100
 		end do
 100     do m = 1, 12
 			if (tdd(m, 1) .eq. 0.0) tcum = 0.0
 			tcum = tcum + tdd(m, 1)
 			tdd(m, 2) = tcum
 			if (fdd(m, 1) .eq. 0.0) fcum = 0.0
 			fcum = fcum + fdd(m, 1)
 			fdd(m, 2) = fcum
 		end do
 		cum_tdd = maxval(tdd(:,2))

 		!loop through each plot to calculate soil dynamics
 		do ip = 1, site%numplots

 			!initialize accumulators
 			aet = 0.0
 			degday = 0.0
 			growdays = 0.0
 			soildays = 0.0
 			drydays_upper = 0.0
 			drydays_base = 0.0
 			flooddays = 0.0
 			outwater = 0.0
 			wl_days = 0.0
            wpdays = 0.0

 			!store depth of thaw from previous year (amlt)
 			amlt = min(site%plots(ip)%soil%active, site%plots(ip)%soil%A_depth)
 			site%plots(ip)%amlt = amlt
 			site%plots(ip)%soil%active = 0.0
 			site%plots(ip)%soil%z_freeze = 0.0

 			!initialize depths of seasonal freeze and thaw
 			xmlt = 0.0
 			xfrz = site%plots(ip)%soil%M_depth +                               &
                site%plots(ip)%soil%O_depth + site%plots(ip)%soil%A_depth

 			!calculate light on the forest floor
 			alff = 1.0*exp(-0.25*site%plots(ip)%cla(1)/plotsize)

 			!calculate drainage conditions and initialize gravimetric
 			!organic and mineral soil moisture contents
 			do l = 1, 2
 				site%plots(ip)%soil%z_drain(l) =                               &
                    (site%plots(ip)%soil%ASAT(l)*(1.0 - amlt) +                &
                    site%plots(ip)%soil%AFC(l)*(amlt - 0.32))/(1.0 - 0.32)

                !must be between field capacity and saturation capacity
 				site%plots(ip)%soil%z_drain(l) =                               &
                    min(site%plots(ip)%soil%z_drain(l),                        &
                    site%plots(ip)%soil%ASAT(l))
 				site%plots(ip)%soil%z_drain(l) =                               &
                    max(site%plots(ip)%soil%z_drain(l),                        &
 					site%plots(ip)%soil%AFC(l))

                !set soil to fully saturated at onset of year
 				if (l .eq. 1) then
 					site%plots(ip)%soil%wc(l) =                                &
                        site%plots(ip)%soil%z_drain(l)*1000.0/                 &
 						site%plots(ip)%soil%O_bulk_dens
 					zh = site%plots(ip)%soil%O_depth +                         &
                        site%plots(ip)%soil%M_depth
 				else
 					site%plots(ip)%soil%wc(l) =                                &
                        site%plots(ip)%soil%z_drain(l)*1000.0/                 &
 						site%plots(ip)%soil%A_bulk_dens
 					zh = site%plots(ip)%soil%A_depth
 				end if

                !soil is completely frozen at start of year
 				site%plots(ip)%soil%H2Oice(l) =                                &
                    site%plots(ip)%soil%z_drain(l)*zh
 				site%plots(ip)%soil%water(l) = 0.0
 				site%plots(ip)%soil%owater(l) = 0.0
 				site%plots(ip)%soil%d_melt(l) = 0.0
 				site%plots(ip)%soil%d_freeze(l) = zh
 				site%plots(ip)%soil%minWC = site%plots(ip)%soil%wc(2)

 			end do

 			!loop on days - thaw depths are calculated monthly so have to
				!increment the months as well
 			m = 1
 			do j = 1, days_per_year
 				if (j .gt. modays(m)) m = m + 1

 				!calculate freeze/thaw depths (xfrz, xmlt) and maximum depths of
                !freeze and thaw
 				call permf(site%plots(ip)%soil, m, 1, alff, tdd, fdd, cfs, xfrz)
 				call permf(site%plots(ip)%soil, m, 2, alff, tdd, fdd, cfs, xmlt)

 				!update maximum depths (z_freeze and active)
 				site%plots(ip)%soil%z_freeze = max((xfrz -                     &
                    site%plots(ip)%soil%M_depth -                              &
                    site%plots(ip)%soil%O_depth), site%plots(ip)%soil%z_freeze)
 				site%plots(ip)%soil%active = max((xmlt -                       &
 					site%plots(ip)%soil%M_depth -                              &
                    site%plots(ip)%soil%O_depth), site%plots(ip)%soil%active)

 				!calculate soil water dynamics for the day
 				call moist(site%plots(ip)%soil, daytemp(j), dayprecip(j),      &
					pot_ev_day(j), site%leaf_area_ind, site%slope, amlt, xmlt, &
					xfrz, tdd, m, site%wd_ind, act_ev_day, aow0_ScaledByMax,   &
					aow0_ScaledByMin, saw0_ScaledByFC, saw0_ScaledByWP,        &
					saw0_ScaledBySAT)

                !update minimum water content for year
 				site%plots(ip)%soil%minWC = min(site%plots(ip)%soil%minWC,     &
                    site%plots(ip)%soil%wc(2))

 				!accumulate variables
 				outwater = outwater + site%plots(ip)%soil%runoff
 				aet = act_ev_day + aet

 				!compute degday, dry days, flood days, and growing season
					!length (days)
 				if (daytemp(j) .ge. min_grow_temp) then
 					degday = degday + (daytemp(j) - min_grow_temp)
 					growdays = growdays + 1.0

 					if ((saw0_ScaledByFC .lt. 0.68)) then
                        drydays_upper = drydays_upper + 1.0
 						if (saw0_ScaledByWP .lt. 1.2) then
 							drydays_base = drydays_base + 1.0
 						end if
 					end if

                    if ((saw0_ScaledByWP .lt. 1.001)) then
                        wpdays = wpdays + 1.0
                    end if

 					if (saw0_ScaledByFC .gt. 1.5) then
 						flooddays = flooddays + 1.01

						if (saw0_ScaledBySAT .gt. 0.85) then
							wl_days = wl_days + 1.0
						end if
 					endif
 				end if

 				if (daytemp(j) .ge. 0.0) then
					soildays = soildays + (daytemp(j) - 0.0)
				end if

 			end do

 			!convert drydays and flooddays to proportion of growing season
 			if (growdays .eq. 0) then
 				drydays_upper = 0.0
 				drydays_base = 0.0
 				flooddays = 0.0
 				wl_days = 0.0
                wpdays = 0.0
 			else
 				tmp = max(min(rain/pet, 1.0), min(aet/pet, 1.0))
 				drydays_upper = ((drydays_upper/growdays) + (1.0 - tmp))/2.0
 				drydays_base = drydays_base/growdays
 				flooddays = flooddays/growdays
 				wl_days = wl_days/growdays
                wpdays = wpdays/growdays
 			endif

 			!convert aet to mm for decomposition
			aet_mm = aet*10.0

			call moss(site%plots(ip)%soil, alff, site%plots(ip)%cla(1),        &
						site%plots(ip)%soil%dec_fuel, drydays_upper)

 			call soiln(site%plots(ip)%soil, aet_mm, site%plots(ip)%cla(1),     &
						soildays, flooddays, n_avail)

 			!set fan to 0.0 (used last year's fan for this year's soiln
            !calculation)
 			site%plots(ip)%soil%fan = 0.0

  		    !set plot-level attributes for year-end values
 			site%plots(ip)%soil%avail_N = n_avail + rain_n
 			site%plots(ip)%soil%runoff = outwater
 		    site%plots(ip)%act_evap_day = aet
 		    site%grow_days = growdays
 			site%deg_days = degday
 			site%plots(ip)%flood_days = flooddays
 			site%plots(ip)%wl_days = wl_days
 			site%plots(ip)%dry_days_upper_layer = drydays_upper
 			site%plots(ip)%dry_days_base_layer = drydays_base
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
 				site%aridity_base = (site%ridity_base + min(rain/pet, 1.0))/10.0
 			endif
 		end if
 		site%aridity = min(rain/pet, 1.0)

 		!set site-level attributes to yearly sums of soil/climate values
 		site%pot_evap_day = pet
 		site%solar = tot_sun
 		site%rain = rain
 		site%heat_load = cum_tdd

	end subroutine BioGeoClimate

	!:.........................................................................:

	subroutine Canopy(site)
		!calculates plot-level LAI and light level
		!Inputs/Outputs:
		!	site:  site instance

		type(SiteData),              intent(inout) :: site

		real, dimension(maxheight)                 :: lvd_c1, lvd_c2
		real, dimension(maxheight)                 :: lvd_c3, lvd_c4
		real, dimension(maxheight)                 :: lai_matrix
		real                                       :: forht, canht
		real                                       :: tlai, lvd_adj, aream2
        real, parameter                            :: xt = -0.40
		integer                                    :: num_species
		integer                                    :: ntrees
		integer                                    :: i, ip, ih, it
		integer                                    :: k, jtmp, iht, m

		!set LAI to 0 and get number of species at site
		site%leaf_area_ind = 0.0
		site%lai_array = 0.0
		lai_matrix = 0.0

		num_species = size(site%species)

		do ip = 1, site%numplots

			!get number of trees on plot
			ntrees = site%plots(ip)%numtrees
			if (ntrees .eq. 0) then
				site%plots(ip)%con_light = 1.0
				site%plots(ip)%dec_light = 1.0
				site%plots(ip)%nutrient = 1.0
				site%plots(ip)%cla = 0.0
			else
				!initialize temporary leaf area arrays
				lvd_c1 = 0.0
				lvd_c2 = 0.0
				lvd_c3 = 0.0
				lvd_c4 = 0.0

				do i = 1, maxheight
					site%plots(ip)%cla(i) = 0.0
				end do

				do it = 1, ntrees

					!get species index
					k = site%plots(ip)%trees(it)%species_index

					!calculate total tree and clear branch bole height
					forht = site%plots(ip)%trees(it)%forska_ht
					iht = min(int(forht), maxheight)
					canht = site%plots(ip)%trees(it)%canopy_ht

					aream2 = 0.160694*site%plots(ip)%trees(it)%diam_bht**2.129
					site%plots(ip)%cla(iht) = site%plots(ip)%cla(iht) + aream2

					tlai = lai_biomass_c(site%plots(ip)%trees(it))

					!accumulate site leaf area
					site%leaf_area_ind = site%leaf_area_ind + tlai

					!calculate canopy depth and divide lai into 1-m sections
					jtmp = max(int(forht) - int(canht) + 1, 1)
					lvd_adj = tlai/float(jtmp)

					!fill temporary arrays with leaf area
					if (site%plots(ip)%trees(it)%conifer) then
						do ih = int(canht), int(forht)
							lai_matrix(ih) = lai_matrix(ih) + lvd_adj
							lvd_c1(ih) = lvd_c1(ih) + lvd_adj
							lvd_c2(ih) = lvd_c2(ih) + lvd_adj
						end do
					else
						do ih = int(canht), int(forht)
							lai_matrix(ih) = lai_matrix(ih) + lvd_adj
							lvd_c2(ih) = lvd_c2(ih) + lvd_adj*0.8
							lvd_c1(ih) = lvd_c1(ih) + lvd_adj
						end do
					end if

				end do

				!cumulative leaf area (Bonan)
				do m = maxheight-1, 1, -1
					site%plots(ip)%cla(m) = site%plots(ip)%cla(m) +            &
                        site%plots(ip)%cla(m+1)
				enddo

				!calculate cumulative leaf area from top down
				lvd_c3(maxheight) = lvd_c1(maxheight)
				lvd_c4(maxheight) = lvd_c2(maxheight)
				!site%plots(ip)%cla(maxheight) = lvd_c2(maxheight)
				do ih = 1, maxheight - 1
					lvd_c3(maxheight - ih) = lvd_c3(maxheight - ih + 1) +      &
						lvd_c1(maxheight - ih)
					lvd_c4(maxheight - ih) = lvd_c4(maxheight - ih + 1) +      &
					  lvd_c2(maxheight - ih)
				end do

				!calculate light level
				do ih = 1, maxheight - 1
					!lai experienced by evergreens reduced by decid.
					!this accounts for part of each year w/out decid.
					  !leaves
					!reduced lai means higher light for evergreens cumulative
					  !for year
					site%plots(ip)%con_light(ih) = exp(xt*lvd_c4(ih + 1)/      &
                        plotsize)
					!deciduous lai normal b/c evergreen and decid both present
					  !when decid has leaves
					site%plots(ip)%dec_light(ih) = exp(xt*lvd_c3(ih + 1)/      &
                        plotsize)
				end do
			end if
		end do

		!get average LAI (m2/m2) for site
		site%leaf_area_ind = site%leaf_area_ind/float(site%numplots)/plotsize
		site%lai_array = lai_matrix/float(site%numplots)/plotsize

	end subroutine Canopy

	!:.........................................................................:

	subroutine Growth(site)
		!calculates annual growth of each tree on each plot
		!Inputs/Outputs:
		!	site:  site instance

		type(SiteData),              intent(inout) :: site

		real,    dimension(maxcells*maxcells)      :: N_stress, perm_resp
		real,    dimension(maxcells*maxcells)      :: bleaf
		real,    dimension(maxcells*maxcells)      :: diam, biom_C
		real,    dimension(maxcells*maxcells)      :: forska_shade
		real,    dimension(maxcells*maxcells)      :: forht
		integer, dimension(maxcells*maxcells)      :: khc, kh
		real                                       :: canopy_shade
		real                                       :: leaf_b
		real                                       :: leafbm
		real                                       :: d_leafb
		real                                       :: canht
		real                                       :: fc_n, pp
		real                                       :: bct, d_bc
		real                                       :: dt, d_bioC
		real                                       :: N_used, N_req
		real                                       :: prim_prod
		real                                       :: net_prim_prodC
		real                                       :: net_prim_prodN
		real                                       :: biomc, biomn
		real                                       :: check
		real                                       :: uconvert
		real                                       :: min_stress
		real                                       :: bcr, bctw, bcs
		real                                       :: d_bcr, d_bctw
		real                                       :: d_bcs
		integer                                    :: stress_cat
		integer                                    :: ntrees, lc
		integer                                    :: num_species
		integer                                    :: it, is, ip, k

		!get number of species at site
		num_species = size(site%species)

		!calculate leaf_b used for calculating conifer leaf litter
		leaf_b = 1.0 + con_leaf_ratio

		do ip = 1, site%numplots

			N_used = 0.0
			net_prim_prodC = 0.0
			net_prim_prodN = 0.0
			biomc = 0.0
			biomn = 0.0
			site%plots(ip)%mature(:) = 0

			!calculate species-level response to drought, over-saturation, and
				!temperature
			do k = 1, num_species
				call temp_rsp(site%species(k), site%deg_days)
				call drought_rsp(site%species(k),                              &
                    site%plots(ip)%dry_days_upper_layer,                       &
					site%plots(ip)%dry_days_base_layer)
				call flood_rsp(site%species(k), site%plots(ip)%wl_days)
			end do

			!initialize 'available' species on plot and get number of trees
			site%plots(ip)%avail_spec = 0.0
			ntrees = site%plots(ip)%numtrees

			!initialize N required
			N_req = 0.0

			if (ntrees .gt. 0) then

				!initialize N required
				N_req = 0.0

				do it = 1, ntrees

					!get species index and update tree
					k = site%plots(ip)%trees(it)%species_index
					call update_tree(site%plots(ip)%trees(it), site%species(k))

					!convenience variables to reduce table lookups
					diam(it) = site%plots(ip)%trees(it)%diam_bht
					canht = site%plots(ip)%trees(it)%canopy_ht
					forht(it) = site%plots(ip)%trees(it)%forska_ht

					!calculate available species
					site%plots(ip)%avail_spec(k) = max(kron(diam(it) -         &
                        site%species(k)%max_diam*growth_thresh),               &
                        site%plots(ip)%avail_spec(k))

					!get canopy and total tree height
					khc(it) = int(canht)
					kh(it) = int(forht(it))

					!get leaf biomass and maximum possible DBH growth
					call leaf_biomass_c(site%plots(ip)%trees(it))
					call max_growth(site%plots(ip)%trees(it))

					!calculate shading effect on tree
					if (site%plots(ip)%trees(it)%conifer) then
						canopy_shade = light_rsp(site%species(k),              &
                            site%plots(ip)%con_light(kh(it)))
						forska_shade(it) = light_rsp(site%species(k),          &
							site%plots(ip)%con_light(khc(it)))
					else
						canopy_shade = light_rsp(site%species(k),              &
							site%plots(ip)%dec_light(kh(it)))
						forska_shade(it) = light_rsp(site%species(k),          &
							site%plots(ip)%dec_light(khc(it)))
					end if

					!calculate environmental stress (excluding N)
					call env_stress(site%plots(ip)%trees(it), canopy_shade,    &
                        min_stress, stress_cat)
					N_stress(it) = min_stress

                    !get stressor
					site%plots(ip)%trees(it)%stressor = stress_cat

                    !calculate permafrost effect
                    call perm_rsp(site%plots(ip)%trees(it)%perm_tol,           &
                        site%plots(ip)%soil%active, perm_resp(it))
                    !perm_resp(it) = 1.0

                    !update stress category if permafrost is more limiting
                    if (perm_resp(it) .lt. N_stress(it)) then
						site%plots(ip)%trees(it)%stressor = 4
					end if

                    site%plots(ip)%trees(it)%env_resp(4) = perm_resp(it)

                    if (site%species(k)%fc_flood .lt. perm_resp(it) .and.      &
						site%species(k)%fc_flood .lt. N_stress(it)) then
							site%plots(ip)%trees(it)%stressor = 5
					end if
					site%plots(ip)%trees(it)%env_resp(5) =                     &
						site%species(k)%fc_flood

					!increment tree diameter using actual DBH growth
					site%plots(ip)%trees(it)%diam_bht =                        &
						site%plots(ip)%trees(it)%diam_bht +                    &
                        site%plots(ip)%trees(it)%diam_max*N_stress(it)*        &
                        perm_resp(it)*site%species(k)%fc_flood

					!compute total height for this diameter
					call forska_height(site%plots(ip)%trees(it))

					!save current value of leaf biomass
					bleaf(it) = site%plots(ip)%trees(it)%leaf_bm

					!update canopy height and leaf biomass with new height
					call stem_shape(site%plots(ip)%trees(it))
					call leaf_biomass_c(site%plots(ip)%trees(it))

					!calculate leaf and fine root growth N requirement
					if (site%species(k)%conifer) then
						N_req = N_req + (leaf_b*                               &
							site%plots(ip)%trees(it)%leaf_bm -                 &
							bleaf(it))/con_leaf_c_n
					else
						N_req = N_req +                                        &
							site%plots(ip)%trees(it)%leaf_bm/dec_leaf_c_n
					end if

					!calculate wood growth N requirement and store old
                      !value of total biomiomass in biom_C
					biom_C(it) = site%plots(ip)%trees(it)%biomC

					!compute new value
					call biomass_c(site%plots(ip)%trees(it))
					call biomass_n(site%plots(ip)%trees(it))

					N_req = N_req + (site%plots(ip)%trees(it)%biomC -          &
						biom_C(it))/stem_c_n

				end do  ! first tree loop to calculate N required for growth

				!conver N_req tonnes N/ha and then convert to percent available
				N_req = max(N_req*hec_to_m2/plotsize, 0.00001)
				N_req = site%plots(ip)%soil%avail_N/N_req

				!calculate species-level response to available N
				do is = 1, num_species
					site%plots(ip)%nutrient(is) = poor_soil_rsp(N_req,         &
                        site%species(is)%lownutr_tol)
				end do

				!calculate actual DBH growth and N used
				do it=1,ntrees

					!get species index and update tree
					k = site%plots(ip)%trees(it)%species_index
					call update_tree(site%plots(ip)%trees(it), site%species(k))

					!calculate nitrogen effect on growth
			        fc_n = min(N_stress(it), site%plots(ip)%nutrient(k))

                    !update stressor if nutrient stress is the most stressful
					if (site%plots(ip)%nutrient(k) .lt. N_stress(it) .and.     &
							site%plots(ip)%nutrient(k) .lt. perm_resp(it) .and. &
							site%plots(ip)%nutrient(k) .lt.                    &
								site%species(k)%fc_flood) then
						site%plots(ip)%trees(it)%stressor = 6
					end if

					site%plots(ip)%trees(it)%env_resp(6) =                     &
						site%plots(ip)%nutrient(k)

					!calculate actual dbh increment
					fc_n = fc_n*perm_resp(it)*site%species(k)%fc_flood
					dt = fc_n*site%plots(ip)%trees(it)%diam_max

					!increment old diameter, store as current tree diameter
					site%plots(ip)%trees(it)%diam_bht = diam(it) + dt

					!get check value for age and growth-related mortality
					pp = min(site%species(k)%max_diam/                         &
						site%species(k)%max_age*0.1, growth_thresh)

                    !check for possible mortality age/growth stress mortality
					if ((dt .le. pp) .or. (fc_n .le. growth_thresh)) then
						site%plots(ip)%trees(it)%mort_count =                  &
							site%plots(ip)%trees(it)%mort_count + 1

						if (site%plots(ip)%trees(it)%mort_count .ge.           &
							mcount) then
							site%plots(ip)%trees(it)%mort_marker = .true.
						endif

					else
						site%plots(ip)%trees(it)%mort_count = 0
						site%plots(ip)%trees(it)%mort_marker = .false.
					endif

					!compute N_used and NPP

					!compute actual height and diameter w/o intermediate
					  !adjustments
					call forska_height(site%plots(ip)%trees(it))
					call stem_shape(site%plots(ip)%trees(it))

					!update biomass, saving leaf biomass into convenience
					  !variable
					call leaf_biomass_c(site%plots(ip)%trees(it))
					leafbm = site%plots(ip)%trees(it)%leaf_bm

					call biomass_c(site%plots(ip)%trees(it))
					call biomass_n(site%plots(ip)%trees(it))

					!calculate change in biomass from previous year
					d_bioC = site%plots(ip)%trees(it)%biomC - biom_C(it)

					!calculate C and N used
					net_prim_prodC = net_prim_prodC + d_bioC
					N_used = N_used + d_bioC/stem_c_n

					!update NPP and N_used from leaf biomass
					if (site%plots(ip)%trees(it)%conifer) then
						prim_prod = leaf_b*leafbm - bleaf(it)
						net_prim_prodC = net_prim_prodC + prim_prod
						N_used = N_used + prim_prod/con_leaf_c_n

						!accumulate total biomass
						biomc = biomc + site%plots(ip)%trees(it)%biomC + leafbm
						biomn = biomn + site%plots(ip)%trees(it)%biomN +       &
							leafbm/con_leaf_c_n
					else
						net_prim_prodC = net_prim_prodC + leafbm
						N_used = N_used + leafbm/dec_leaf_c_n

						!accumulate total biomass (no leaves)
						biomc = biomc + site%plots(ip)%trees(it)%biomC
						biomn = biomn + site%plots(ip)%trees(it)%biomN
					end if

                    !calculate stand age as maximum tree age
					site%plots(ip)%stand_age = max(site%plots(ip)%stand_age,   &
                        site%plots(ip)%trees(it)%tree_age)

				end do  ! second tree loop for calculating actual DBH growth

				!update clear branch bole height and wood/leaf litter fall
				do it = 1, ntrees

					!get species index and update tree
					k = site%plots(ip)%trees(it)%species_index
					call update_tree(site%plots(ip)%trees(it), site%species(k))

					! if (site%species(k)%unique_id .eq. 'PICEmari') then
					! 	if (site%plots(ip)%trees(it)%tree_age .gt. 30 .and.    &
                    !         site%plots(ip)%trees(it)%diam_bht .gt. 10.0) then
					! 		site%plots(ip)%mature(k) = site%plots(ip)%mature(k) + 1
					! 	end if
					! else
					! 	if (site%plots(ip)%trees(it)%diam_bht .gt. 10.0) then
					! 		site%plots(ip)%mature(k) = site%plots(ip)%mature(k) + 1
					! 	end if
					! end if

                    if (site%plots(ip)%trees(it)%tree_age .ge.                 &
                        site%species(k)%recr_age .and.                         &
                        site%plots(ip)%trees(it)%diam_bht .gt. 10.0) then
                        site%plots(ip)%mature(k) = site%plots(ip)%mature(k) + 1
                    end if

					!get total tree height
					forht(it) = site%plots(ip)%trees(it)%forska_ht

					!check for lower branch thinning
					check = min(site%species(k)%fc_degday,                     &
								site%species(k)%fc_drought,                    &
								forska_shade(it), site%plots(ip)%nutrient(k))* &
								perm_resp(it)*site%species(k)%fc_flood

					!check if lower branches fall (thus increasing clear branch
					  !bole height)
					if (check .le. growth_thresh) then

						!increment clear branch bole height
                        khc(it) = khc(it) + 1

                        if (khc(it) .lt.  int(forht(it))) then
							site%plots(ip)%trees(it)%canopy_ht =               &
								float(khc(it)) + 0.01

                            !update diameter at CH
                            call stem_shape(site%plots(ip)%trees(it))

							!save old biomass
							bct = site%plots(ip)%trees(it)%biomC

							!save old root, twig, stem C biomass
							bcr = site%plots(ip)%trees(it)%rootC
							bctw = site%plots(ip)%trees(it)%twigC
							bcs = site%plots(ip)%trees(it)%stemC

							!update biomass C and N
							call biomass_c(site%plots(ip)%trees(it))
							call biomass_n(site%plots(ip)%trees(it))

							!how much wood litter did we lose?
							d_bc = bct - site%plots(ip)%trees(it)%biomC
							d_bcr = bcr - site%plots(ip)%trees(it)%rootC
							d_bctw = bctw - site%plots(ip)%trees(it)%twigC
							d_bcs = bcs -  site%plots(ip)%trees(it)%stemC

                            !add litter loss to litter pools
                            !here we divide by 0.45 because soiln subroutine
                              !calculates weight loss not C loss

							!roots
							site%plots(ip)%soil%litter(13) =                   &
								site%plots(ip)%soil%litter(13) + d_bcr/0.45

							!twigs
							site%plots(ip)%soil%litter(16) =                   &
								site%plots(ip)%soil%litter(16) + d_bctw/0.45

							!tree DBH < 10 cm go into smallwood, otherwise into
							  !large wood
							if (site%plots(ip)%trees(it)%diam_bht  .gt.        &
								10.0) then
								site%plots(ip)%soil%litter(15) =               &
								site%plots(ip)%soil%litter(15) + d_bcs/0.45
							else
								site%plots(ip)%soil%litter(14) =               &
								site%plots(ip)%soil%litter(14) + d_bcs/0.45
							end if

							!save previous value of leaf biomass
							leafbm = site%plots(ip)%trees(it)%leaf_bm

							!update leaf bm and get difference
							call leaf_biomass_c(site%plots(ip)%trees(it))
							d_leafb = leafbm - site%plots(ip)%trees(it)%leaf_bm

							!add that litter to correct leaf litter class
							if (site%species(k)%conifer) then
								lc = site%species(k)%litter_class
                                site%plots(ip)%soil%litter(lc) =               &
                                    site%plots(ip)%soil%litter(lc) +           &
                                    d_leafb*leaf_b/0.45
							else
								lc = site%species(k)%litter_class
								site%plots(ip)%soil%litter(lc) =               &
									site%plots(ip)%soil%litter(lc) +           &
									d_leafb/0.45
							end if
						end if
					end if
				end do  !end third tree loop for lower branch thinning
			end if

			!update plot-level soil characteristics
			uconvert = hec_to_m2/plotsize
			N_used = N_used*uconvert
			site%plots(ip)%soil%net_prim_prodC = net_prim_prodC*uconvert
			site%plots(ip)%soil%N_used = N_used
			site%plots(ip)%soil%net_prim_prodN = N_used
			site%plots(ip)%soil%avail_N = max(0.0,                             &
				site%plots(ip)%soil%avail_N - site%plots(ip)%soil%N_used)

		end do

	end subroutine Growth

	!:.........................................................................:

	subroutine Mortality(site, year)
		!determines which trees die by age, stress, or disturbances, and adds
          !their biomass to soil litter pools
        !Inputs/Outputs:
		!	site:  site instance
		!Inputs:
		!	year: simulation year


		type(SiteData),        intent(inout) :: site
		integer,               intent(in)    :: year

		real                                 :: biomc, biomn
		real                                 :: leaf_b
		real                                 :: fire_prob, wind_prob
		real                                 :: wind_int, wind_std
		real                                 :: fire_p
		real                                 :: leaf_bm, bmc
		real                                 :: NPP_loss, NPPn_loss
		real                                 :: bcr, bctw, bcs
		real                                 :: av_fuel, N_cons
		real, parameter                      :: fol_rt = 0.5
		real, parameter                      :: tw_rt = 0.1
		real, parameter                      :: stem_rt = 0.05
		integer                              :: num_species
		integer                              :: jt, k
		integer                              :: ip, it, ht
		integer                              :: snum, lc

		logical                              :: age_survive
		logical                              :: growth_survive
		logical                              :: fire_survive
		logical                              :: wind_survive

		real, parameter                      :: tfc1 = 0.8
		real, parameter                      :: tfc2 = 0.2
		real, parameter                      :: bfc = 0.4
		real                                 :: consRoot, burn
		real                                 :: fan1, fan2, fan3, fan4
		integer                              :: trow, tcol

		!get number of species on site
		num_species = size(site%species)

		!set up conifer litter loss
		leaf_b = 1.0 + con_leaf_ratio

		do ip = 1, site%numplots

			ht = 0
			site%plots(ip)%num_dead = 0

			!initialize biomass and NPP losses
			biomc = 0.0
			biomn = 0.0
			NPP_loss = 0.0
			NPPn_loss = 0.0
			fan1 = 0.0
			fan2 = 0.0
			fan3 = 0.0
			fan4 = 0.0

			!set root fire consumption to 0.0
			consRoot = 0.0

			!set disturbance type to 0.0
			site%plots(ip)%d_type = 0.0

			!get random numbers for fire and wind throw
			fire_prob = urand()
			wind_prob = urand()

			!increase or decrease fire probability based on aridity
			if (year .ge. 10) then
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

   
			!check for fire or windthrow
			if (fire_prob < fire_p .or. wind_prob < site%wind_prob) then

				if (site%plots(ip)%numtrees > 0) then !must be some trees

					if (fire_prob < fire_p) then !fire occurs

						!fire disturbance occurs
						site%plots(ip)%fire = 1

						!calculate and consume available litter fuels
						call forest_fuels(site%plots(ip)%soil,                 &
                            site%plots(ip)%dry_days_upper_layer, av_fuel,      &
                            N_cons, consRoot)

						!kill trees that died by fire, age, or low
                          !growth - only copy surviving trees, rest go into soil

						!initialize new and dead tree counters
						ht = 0
						jt = 0

						do it = 1, site%plots(ip)%numtrees

							!get species index and update tree
							k = site%plots(ip)%trees(it)%species_index
							call update_tree(site%plots(ip)%trees(it),         &
                                site%species(k))

                            !get leaf biomass
							call leaf_biomass_c(site%plots(ip)%trees(it))
							leaf_bm = site%plots(ip)%trees(it)%leaf_bm

							!check for growth and age survival
							call growth_survival(site%plots(ip)%trees(it),     &
                                growth_survive)
							call age_survival(site%plots(ip)%trees(it),        &
                                age_survive)

							!check for fire survival
							call fire_survival(site%plots(ip)%trees(it),       &
                                av_fuel, fire_survive)

							if (growth_survive .and. age_survive .and.         &
                                fire_survive) then

								!tree survives, copy its attributes from old
								  !list to new

								jt = jt + 1
								site%plots(ip)%trees(it)%tree_age =            &
                                site%plots(ip)%trees(it)%tree_age + 1

								call copy_tree(site%plots(ip)%trees(jt),       &
                                    site%plots(ip)%trees(it))

								!calculate leaf litter
								if (site%species(k)%conifer) then
									lc = site%species(k)%litter_class
									site%plots(ip)%soil%litter(lc) =           &
										site%plots(ip)%soil%litter(lc) +       &
                                        leaf_bm*(leaf_b - 1.0)/0.45

									!accumulate NPP losses
									NPP_loss = NPP_loss + leaf_bm*(leaf_b - 1.0)
									NPPn_loss = NPPn_loss +                    &
                                        leaf_bm*(leaf_b - 1.0)/con_leaf_c_n
								else

									lc = site%species(k)%litter_class
									site%plots(ip)%soil%litter(lc) =           &
										site%plots(ip)%soil%litter(lc) +       &
										leaf_bm/0.45

									!acumulate NPP losses
									NPP_loss = NPP_loss + leaf_bm
									NPPn_loss = NPPn_loss + leaf_bm/dec_leaf_c_n
								end if

							else

								!tree dies, copy to list of dead trees

								ht = ht + 1
								call copy_tree(site%plots(ip)%deadtrees(ht),   &
                                    site%plots(ip)%trees(it))

								!set cells of that tree to unfilled
								trow = site%plots(ip)%trees(it)%row
								tcol = site%plots(ip)%trees(it)%col
								site%plots(ip)%cells(trow, tcol) = 0

								!get most limiting growth factor
								snum = site%plots(ip)%trees(it)%stressor

								if (fire_survive .eq. .false.) then

                                    !died by fire, fire consumes some of each
                                      !litter class. also calculate N
                                      !mineralized in tree burning
									site%plots(ip)%d_type(7) =                 &
                                        site%plots(ip)%d_type(7) +             &
										site%plots(ip)%trees(it)%biomC + leaf_bm
									site%plots(ip)%deadtrees(ht)%stressor = 7

									!roots - fire consumption from Bonan
									bcr = site%plots(ip)%trees(it)%rootC
									burn = bcr*consRoot
									fan1 = fan1 + burn*litter_params(13, 2)*   &
                                        (1 - N_cons)
									bcr = bcr - burn
									site%plots(ip)%soil%litter(13) =           &
                                        site%plots(ip)%soil%litter(13) +       &
                                        bcr/0.45

									!twigs
									bctw = site%plots(ip)%trees(it)%twigC
									burn = bctw*tw_rt
									bctw = bctw - burn
									fan4 = fan4 + burn*litter_params(16, 2)*   &
										(1 - N_cons)
									site%plots(ip)%soil%litter(16) =           &
                                        site%plots(ip)%soil%litter(16) +       &
                                        bctw/0.45

									!stems
									bcs = site%plots(ip)%trees(it)%stemC
									burn = bcs*stem_rt
									bcs = bcs - burn

									if (site%plots(ip)%trees(it)%diam_bht      &
										.gt. 10.0) then
										site%plots(ip)%soil%litter(15) =       &
											site%plots(ip)%soil%litter(15) +   &
											bcs/0.45
										fan3 = fan3 +                          &
											burn*litter_params(15, 2)*         &
											(1 - N_cons)
									else
										site%plots(ip)%soil%litter(14) =       &
											site%plots(ip)%soil%litter(14) +   &
											bcs/0.45
										fan3 = fan3 +                          &
											burn*litter_params(14, 2)*         &
												(1 - N_cons)
									end if

									!leaves
									lc = site%species(k)%litter_class

									burn = leaf_bm*fol_rt
									if (site%species(k)%conifer) then
										leaf_bm = leaf_bm*leaf_b - burn
									else
										leaf_bm = leaf_bm - burn
									end if
                                    site%plots(ip)%soil%litter(lc) =           &
                                    site%plots(ip)%soil%litter(lc) +           &
										leaf_bm/0.45
									fan2 = fan2 +                              &
										burn*litter_params(lc, 2)*(1 - N_cons)

								else if (growth_survive .eq. .false. .or.      &
                                    age_survive .eq. .false.) then

                                    !died from growth/age-related stress, all
                                      !litter goes into soil accumulate biomass
                                      !into correct stress category

									site%plots(ip)%d_type(snum) =              &
                                        site%plots(ip)%d_type(snum) +          &
									    site%plots(ip)%trees(it)%biomC + leaf_bm

									!add total tree litter components to all
									  !categories

									!roots
									bcr = site%plots(ip)%trees(it)%rootC
									site%plots(ip)%soil%litter(13) =           &
										site%plots(ip)%soil%litter(13) +       &
										bcr/0.45

									!twigs
									bctw = site%plots(ip)%trees(it)%twigC
									site%plots(ip)%soil%litter(16) =           &
									  site%plots(ip)%soil%litter(16) +         &
									  bctw/0.45

									!stems
									bcs = site%plots(ip)%trees(it)%stemC
									if (site%plots(ip)%trees(it)%diam_bht      &
										.gt. 10.0) then
										site%plots(ip)%soil%litter(15) =       &
											site%plots(ip)%soil%litter(15) +   &
											bcs/0.45
									else
										site%plots(ip)%soil%litter(14) =       &
											site%plots(ip)%soil%litter(14) +   &
											bcs/0.45
									end if

									!leaves
									lc = site%species(k)%litter_class
									if (site%species(k)%conifer) then
										site%plots(ip)%soil%litter(lc) =       &
											site%plots(ip)%soil%litter(lc) +   &
                                            leaf_bm*leaf_b/0.45
									else
										site%plots(ip)%soil%litter(lc) =       &
									  		site%plots(ip)%soil%litter(lc) +   &
                                            leaf_bm/0.45
									end if

								end if

								!get biomass of tree
								bmc = site%plots(ip)%trees(it)%biomC

								!acumulate NPP losses
								if (site%species(k)%conifer) then
									NPP_loss = NPP_loss + bmc + leaf_bm*leaf_b
									NPPn_loss = NPPn_loss + bmc/stem_c_n +     &
                                        leaf_bm/con_leaf_c_n*leaf_b
								else
									NPP_loss = NPP_loss + bmc + leaf_bm
									NPPn_loss = NPPn_loss + bmc/stem_c_n +     &
                                        leaf_bm/dec_leaf_c_n
								end if
							end if  !end mortality check
						end do !end tree loop for fire occurrence

						!set number of trees, wind to 0
						site%plots(ip)%numtrees = jt
						site%plots(ip)%wind = 0
						site%plots(ip)%num_dead = ht

					else !wind disturbance occurs

                        site%plots(ip)%fire = 0
						site%plots(ip)%wind = 1
						wind_std = 0.5
						wind_int = 0.2

						!calculate wind intensity
						site%plots(ip)%wind_cat = nrand(wind_int, wind_std)

						!make sure it is between 0.0 and 1.0
						do while (site%plots(ip)%wind_cat .lt. 0.01            &
                            .or. site%plots(ip)%wind_cat .gt. 2.0)
							site%plots(ip)%wind_cat=nrand(wind_int, wind_std)
						end do

						!if 1.0, regeneration has to wait
						if (site%plots(ip)%wind_cat .ge. 1.0) then
							site%plots(ip)%windCount = 3
						else
							site%plots(ip)%windCount = 0
						endif

						!initialize counters for live and dead trees
						jt = 0
						ht = 0

                        do it = 1, site%plots(ip)%numtrees

							!get species index and update tree
							k = site%plots(ip)%trees(it)%species_index
							call update_tree(site%plots(ip)%trees(it),         &
                                site%species(k))

							!get leaf biomass
							call leaf_biomass_c(site%plots(ip)%trees(it))
							leaf_bm = site%plots(ip)%trees(it)%leaf_bm

							!check for growth and age mortality
							call growth_survival(site%plots(ip)%trees(it),     &
                                growth_survive)
							call age_survival(site%plots(ip)%trees(it),        &
                                age_survive)

							!check for wind survival
							call wind_survival(site%plots(ip)%trees(it),       &
                                site%plots(ip)%wind_cat, wind_survive)

							if (growth_survive .and. age_survive .and.         &
                                wind_survive) then

                                !tree survives, copy from old list to new
								jt = jt + 1
								site%plots(ip)%trees(it)%tree_age =            &
									site%plots(ip)%trees(it)%tree_age + 1

								call copy_tree(site%plots(ip)%trees(jt),       &
                                    site%plots(ip)%trees(it))

                                !calculate leaf litter
								if (site%species(k)%conifer) then
									lc = site%species(k)%litter_class
									site%plots(ip)%soil%litter(lc) =           &
										site%plots(ip)%soil%litter(lc) +       &
                                        leaf_bm*(leaf_b - 1.0)/0.45

									!acumulate NPP losses
									NPP_loss = NPP_loss + leaf_bm*(leaf_b - 1.0)
									NPPn_loss = NPPn_loss +                    &
                                        leaf_bm*(leaf_b-1.0)/con_leaf_c_n
								else
									lc = site%species(k)%litter_class
									site%plots(ip)%soil%litter(lc) =           &
										site%plots(ip)%soil%litter(lc) +       &
                                        leaf_bm/0.45

									!accumulate NPP losses
									NPP_loss = NPP_loss + leaf_bm
									NPPn_loss = NPPn_loss + leaf_bm/dec_leaf_c_n
								end if

							else

								!tree dies, copy to list of dead trees
								bmc = site%plots(ip)%trees(it)%biomC
								ht = ht + 1
								call copy_tree(site%plots(ip)%deadtrees(ht),   &
                                    site%plots(ip)%trees(it))

                                !set cells of that tree to unfilled
								trow = site%plots(ip)%trees(it)%row
								tcol = site%plots(ip)%trees(it)%col
								site%plots(ip)%cells(trow, tcol) = 0

								!get most limiting growth factor
								snum = site%plots(ip)%trees(it)%stressor

								if (growth_survive .eq. .false. .or.           &
                                    age_survive .eq. .false.) then

									site%plots(ip)%d_type(snum) =              &
                                        site%plots(ip)%d_type(snum) +          &
									    site%plots(ip)%trees(it)%biomC + leaf_bm

									!add total tree litter components to all
									  !categories

									!roots
									bcr = site%plots(ip)%trees(it)%rootC
									site%plots(ip)%soil%litter(13) =           &
										site%plots(ip)%soil%litter(13) +       &
										bcr/0.45

									!twigs
									bctw = site%plots(ip)%trees(it)%twigC
									site%plots(ip)%soil%litter(16) =           &
										site%plots(ip)%soil%litter(16) +       &
										bctw/0.45

									!stems
									bcs = site%plots(ip)%trees(it)%stemC
									if (site%plots(ip)%trees(it)%diam_bht      &
										.gt. 10.0) then
										site%plots(ip)%soil%litter(15) =       &
											site%plots(ip)%soil%litter(15) +   &
											bcs/0.45
									else
										site%plots(ip)%soil%litter(14) =       &
											site%plots(ip)%soil%litter(14) +   &
											bcs/0.45
									end if

									!leaves
									lc = site%species(k)%litter_class
									if (site%species(k)%conifer) then
										site%plots(ip)%soil%litter(lc) =       &
											site%plots(ip)%soil%litter(lc) +   &
                                            leaf_bm*leaf_b/0.45
									else
										site%plots(ip)%soil%litter(lc) =       &
											site%plots(ip)%soil%litter(lc) +   &
											leaf_bm/0.45
									end if

								else if (wind_survive .eq. .false.) then

                                    !died by windthrow, add to windthrow mort
									  !markers
                                    site%plots(ip)%d_type(8) =                 &
                                        site%plots(ip)%d_type(8) +             &
										site%plots(ip)%trees(it)%biomC + leaf_bm

									site%plots(ip)%deadtrees(ht)%stressor = 8

									!add total tree litter components to all
									  !categories

									!roots
									bcr = site%plots(ip)%trees(it)%rootC
									site%plots(ip)%soil%litter(13) =           &
										site%plots(ip)%soil%litter(13) +       &
										bcr/0.45

									!twigs
									bctw = site%plots(ip)%trees(it)%twigC
									site%plots(ip)%soil%litter(16) =           &
										site%plots(ip)%soil%litter(16) +       &
										bctw/0.45

									!stems
									bcs = site%plots(ip)%trees(it)%stemC
									if (site%plots(ip)%trees(it)%diam_bht      &
										.gt. 10.0) then
										site%plots(ip)%soil%litter(15) =       &
											site%plots(ip)%soil%litter(15) +   &
											bcs/0.45
									else
										site%plots(ip)%soil%litter(14) =       &
											site%plots(ip)%soil%litter(14) +   &
											bcs/0.45
									end if

									!leaves
									lc = site%species(k)%litter_class

									if (site%species(k)%conifer) then
										site%plots(ip)%soil%litter(lc) =       &
									  	    site%plots(ip)%soil%litter(lc) +   &
                                            leaf_bm*leaf_b/0.45
									else
										site%plots(ip)%soil%litter(lc) =       &
									  	    site%plots(ip)%soil%litter(lc) +   &
                                            leaf_bm/0.45
									end if
								end if

								!acumulate NPP losses
								if (site%species(k)%conifer) then
									NPP_loss = NPP_loss + bmc + leaf_bm*leaf_b
									NPPn_loss = NPPn_loss + bmc/stem_c_n +     &
                                        leaf_bm/con_leaf_c_n*leaf_b
								else
									NPP_loss = NPP_loss + bmc + leaf_bm
									NPPn_loss = NPPn_loss + bmc/stem_c_n +     &
                                        leaf_bm/dec_leaf_c_n
								end if
							end if !end mortality check for windthrow
						end do ! end tree loop for wind disturbance

						!set number of trees
						site%plots(ip)%numtrees = jt
						site%plots(ip)%num_dead = ht

					end if !end if fire or wind throw
				end if !end if any trees

			else !no disturbances occur

				if (site%plots(ip)%numtrees > 0) then

					jt = 0
					ht = 0

					!set fire and wind to 0
					site%plots(ip)%fire = 0
					site%plots(ip)%wind = 0

					do it = 1, site%plots(ip)%numtrees

						!get species index and update tree
						k = site%plots(ip)%trees(it)%species_index
						call update_tree(site%plots(ip)%trees(it),             &
                            site%species(k))

						!get leaf biomass

						call leaf_biomass_c(site%plots(ip)%trees(it))
						leaf_bm = site%plots(ip)%trees(it)%leaf_bm

                        !check for age and growth survival
						call growth_survival(site%plots(ip)%trees(it),         &
                            growth_survive)
						call age_survival(site%plots(ip)%trees(it),            &
                            age_survive)

						if (growth_survive .and. age_survive) then

							!tree survives, copy to new list
							jt = jt + 1
							site%plots(ip)%trees(it)%tree_age=                 &
								site%plots(ip)%trees(it)%tree_age + 1

							call copy_tree(site%plots(ip)%trees(jt),           &
                                site%plots(ip)%trees(it))

							!calculate litterfall
							if (site%species(k)%conifer) then
								lc = site%species(k)%litter_class
                                site%plots(ip)%soil%litter(lc) =               &
									site%plots(ip)%soil%litter(lc) +           &
                                    leaf_bm*(leaf_b - 1.0)/0.45

								NPP_loss = NPP_loss + leaf_bm*(leaf_b - 1.0)
								NPPn_loss = NPPn_loss +                        &
                                    leaf_bm*(leaf_b-1.0)/con_leaf_c_n
							else
								lc = site%species(k)%litter_class
                                site%plots(ip)%soil%litter(lc) =               &
									site%plots(ip)%soil%litter(lc) +           &
									leaf_bm/0.45

								NPP_loss = NPP_loss + leaf_bm
								NPPn_loss = NPPn_loss + leaf_bm/dec_leaf_c_n
							end if

						else
							!tree dies, copy to list of dead trees
							ht = ht + 1
							call copy_tree(site%plots(ip)%deadtrees(ht),       &
                                site%plots(ip)%trees(it))

							!set cells of that tree to unfilled
							trow = site%plots(ip)%trees(it)%row
							tcol = site%plots(ip)%trees(it)%col
							site%plots(ip)%cells(trow, tcol) = 0

							!get most limiting growth factor
							snum = site%plots(ip)%trees(it)%stressor

							if (growth_survive .eq. .false. .or. age_survive   &
                                .eq. .false.) then

								site%plots(ip)%d_type(snum) =                  &
                                    site%plots(ip)%d_type(snum) +              &
								    site%plots(ip)%trees(it)%biomC + leaf_bm

								!add total tree litter components to all
								  !categories

								!roots
								bcr = site%plots(ip)%trees(it)%rootC
								site%plots(ip)%soil%litter(13) =               &
									site%plots(ip)%soil%litter(13) + bcr/0.45

								!twigs
								bctw = site%plots(ip)%trees(it)%twigC
								site%plots(ip)%soil%litter(16) =               &
									site%plots(ip)%soil%litter(16) + bctw/0.45

								!stems
								bcs = site%plots(ip)%trees(it)%stemC
								if (site%plots(ip)%trees(it)%diam_bht .gt.     &
									10.0) then
									site%plots(ip)%soil%litter(15) =           &
									    site%plots(ip)%soil%litter(15) +       &
									    bcs/0.45
								else
									site%plots(ip)%soil%litter(14) =           &
									    site%plots(ip)%soil%litter(14) +       &
									    bcs/0.45
								end if

								!leaves
								lc = site%species(k)%litter_class
								if (site%species(k)%conifer) then
									site%plots(ip)%soil%litter(lc) =           &
									  site%plots(ip)%soil%litter(lc) +         &
                                      leaf_bm*leaf_b/0.45
								else
									site%plots(ip)%soil%litter(lc) =           &
									  site%plots(ip)%soil%litter(lc) +         &
                                      leaf_bm/0.45
								end if

                            end if

                            !get biomass
						    bmc = site%plots(ip)%trees(it)%biomC

						    !decrease NPP
						    if (site%species(k)%conifer) then
							    NPP_loss = NPP_loss + bmc + leaf_bm*leaf_b
							    NPPn_loss = NPPn_loss + bmc/stem_c_n +         &
                                    leaf_bm/con_leaf_c_n*leaf_b
						    else
							    NPP_loss = NPP_loss + bmc + leaf_bm
							    NPPn_loss = NPPn_loss + bmc/stem_c_n +         &
                                    leaf_bm/dec_leaf_c_n
						    end if
						end if !end mortality check
					end do !end tree loop for no disturbances

					site%plots(ip)%num_dead = ht
					site%plots(ip)%numtrees = jt

				end if !end if any trees


			end if !end disturbances/no disturbances check

			!update N mineralized (tonnes N/ha)
			site%plots(ip)%soil%fan = (fan1 + fan2 + fan3 + fan4)/             &
				plotsize/0.0001

		end do
	return

	end subroutine Mortality

	!:.........................................................................:

	subroutine Renewal(site)
		!calculates the seedling and seed banks of each species and regenerates
          !new trees
        !Inputs/Outputs:
		!	site: site instance
		!	year: simulation year

		type(SiteData),              intent(inout) :: site

		real,    dimension(size(site%species))     :: regrowth
		real,    dimension(size(site%species))     :: prob
		real,    dimension(size(site%species))     :: perm_resp
		integer, dimension(maxcells*maxcells)      :: r_empty, c_empty
		integer, dimension(:), allocatable         :: t_vals
		real                                       :: net_prim_prod
		real                                       :: N_used
		real                                       :: leaf_b
		real                                       :: growmax, grow_cap
		real                                       :: probsum, q0
		real                                       :: z, zz
		real                                       :: uconvert
		real                                       :: canopy_shade
		real                                       :: min_stress
		real                                       :: org_resp, ilitter
        real                                       :: germinants
		integer                                    :: org_fact
		integer                                    :: num_species
		integer                                    :: max_renew
		integer                                    :: new, kh
		integer                                    :: stress_cat
		integer                                    :: ip, is, it, k
		integer                                    :: lc, lit, t, r, c
		integer                                    :: nrenew, irenew
		integer                                    :: n_empty

		!calculate conifer leaf ratio
		leaf_b = 1.0 + con_leaf_ratio

		!get number of species at site
		num_species = size(site%species)

		do ip = 1, site%numplots

			!initialize new trees, NPP, N_used
			new = 0
			net_prim_prod = 0.0
			N_used = 0.0

			if (site%plots(ip)%soil%avail_N .gt. 0.0) then

				!recalculate species-level response to soil moisture and
				  !temperature
                !have to do this each time we run through plot due to the way
                  !it is currently set up
				do k = 1, num_species
					call temp_rsp(site%species(k), site%deg_days)
					call drought_rsp(site%species(k),                          &
                        site%plots(ip)%dry_days_upper_layer,                   &
						site%plots(ip)%dry_days_base_layer)
					call flood_rsp(site%species(k), site%plots(ip)%wl_days)
					call perm_rsp(site%species(k)%perm_tol,                    &
                        site%plots(ip)%soil%active, perm_resp(k))
				end do

				!check to make sure we aren't still waiting on wind counter
				if (site%plots(ip)%windCount .eq. 0) then

                    !calculate species-level regrowth, and growcap and growmax
                      !for plot
					growmax = 0.0
					do is = 1, num_species
						grow_cap = min(site%species(is)%fc_degday,             &
                            !site%species(is)%fc_drought,                       &
                            site%plots(ip)%nutrient(is))*site%species(is)%fc_drought


						if (site%plots(ip)%numtrees .eq. 0) then
							regrowth(is) = grow_cap
						else
							if (site%species(is)%conifer) then
								regrowth(is) = min(grow_cap,                   &
                                    light_rsp(site%species(is),                &
                                    site%plots(ip)%con_light(1)))*             &
                                    site%species(is)%fc_flood*perm_resp(is)
							else
								regrowth(is) =  min(grow_cap,                  &
                                    light_rsp(site%species(is),                &
									site%plots(ip)%dec_light(1)))*             &
                                    site%species(is)%fc_flood*perm_resp(is)
							end if
						end if

						!check for enough mature trees
						if (site%plots(ip)%mature(is) .le. 5 .and.             &
                            site%species(is)%conifer) then
							regrowth(is) = regrowth(is)*0.5
						end if

						if (regrowth(is).le. growth_min) regrowth(is) = 0.0

                        !calculate maximum regrowth across all species
						growmax = max(growmax, regrowth(is))
					end do

					!compute the max renew number
					max_renew = min(int((maxcells*maxcells)*growmax) -         &
                        site%plots(ip)%numtrees, int(maxcells*maxcells))

					if (max_renew .lt. 0) max_renew = 0

					!compute actual number of renewable trees
					nrenew = min(max(max_renew, 3), int(maxcells*maxcells) -   &
						site%plots(ip)%numtrees)

                    !if (site%plots(ip)%fire) nrenew = 0

					!compute the seed bank size and seedling bank size (1/m^2)
					do is = 1, num_species
						site%plots(ip)%seedbank(is) =                          &
                            site%plots(ip)%seedbank(is) +                      &
                            site%species(is)%invader +                         &
                            site%species(is)%seed_num*                         &
                            site%plots(ip)%avail_spec(is)

						!get species-level regeneration response to fire
						call fire_rsp(site%species(is), site%plots(ip)%fire)

						!update seedling bank
						site%plots(ip)%seedbank(is) =                         &
                            site%plots(ip)%seedbank(is)*                      &
                            site%species(is)%fc_fire

                        if (.not. site%plots(ip)%fire) then

						    !put seeds into seedling bank if enough regrowth
						    if (regrowth(is) .ge. growth_min) then

                                if (site%species(is)%unique_id == 'PICEmari') then

                                    germinants = site%plots(ip)%seedbank(is)*  &
                                        site%pc_germ

        						    !ability to reproduce on moss-covered soil
        						    org_fact = site%species(is)%org_tol

            						org_resp = exp(org_gf(org_fact)*           &
            							(site%plots(ip)%soil%O_depth +         &
            							site%plots(ip)%soil%M_depth))

        						    germinants = germinants*org_resp

            						site%plots(ip)%seedbank(is) =              &
                                        site%plots(ip)%seedbank(is) - germinants

                                    site%plots(ip)%seedling(is) =              &
                                        site%plots(ip)%seedling(is) +          &
                                        germinants

                                    !decrement seedbank
                                    site%plots(ip)%seedbank(is) =              &
                                        site%plots(ip)%seedbank(is)*           &
                                        site%species(is)%seed_surv

                                else

                                    germinants = site%plots(ip)%seedbank(is)

                                    !ability to reproduce on moss-covered soil
                                    org_fact = site%species(is)%org_tol

                                    org_resp = exp(org_gf(org_fact)*           &
                                        (site%plots(ip)%soil%O_depth +         &
                                        site%plots(ip)%soil%M_depth))

                                    germinants = germinants*org_resp

                                    site%plots(ip)%seedbank(is) =              &
                                        site%plots(ip)%seedbank(is) - germinants

                                    site%plots(ip)%seedling(is) =              &
                                        site%plots(ip)%seedling(is) +          &
                                        germinants

                                    !decrement seedbank
                                    site%plots(ip)%seedbank(is) =              &
                                        site%plots(ip)%seedbank(is)*           &
                                        site%species(is)%seed_surv
                                endif
						    else
    							!decrement seedbank
    							site%plots(ip)%seedbank(is) =                  &
                                    site%plots(ip)%seedbank(is)*               &
    								site%species(is)%seed_surv
						    endif

                            !check for layering
                            if (site%species(is)%layering) then
                                if ((site%plots(ip)%soil%M_depth +             &
                                    site%plots(ip)%soil%O_depth) .gt. 0.05) then
    								site%plots(ip)%seedling(is) =              &
                                        site%plots(ip)%seedling(is)*1.8
    							end if
                            end if

                            !sprouts
                            site%plots(ip)%seedling(is) =                      &
    							site%plots(ip)%seedling(is) +                  &
    							site%species(is)%sprout_num*                   &
                                site%plots(ip)%avail_spec(is)

    						!calculate seedling number
    						site%plots(ip)%seedling_number =                   &
                                max(kron(site%plots(ip)%seedling(is)),         &
                                site%plots(ip)%seedling_number)

                            !convert seedling bank to per plot
    						site%plots(ip)%seedling(is) =                      &
                                site%plots(ip)%seedling(is)*plotsize
                        end if
					end do

					!calculate probability of regeneration
					probsum = 0.0
					do is = 1, num_species
						prob(is) = site%plots(ip)%seedling(is)*regrowth(is)
						probsum = probsum + prob(is)
					end do

                else !still waiting on windCount

					if (site%plots(ip)%windCount .eq. 1) then

						!final count year
						do is = 1, num_species

							!calculate grow cap
							grow_cap = min(site%species(is)%fc_degday,         &
                                site%species(is)%fc_drought,                   &
								site%plots(ip)%nutrient(is))*perm_resp(is)*    &
                                site%species(is)%fc_flood

							!calculate seedling and seedling number
							site%plots(ip)%seedling(is) =                      &
                                site%plots(ip)%seedling(is)*plotsize

							site%plots(ip)%seedling_number =                   &
								max(kron(site%plots(ip)%seedling(is)),         &
								site%plots(ip)%seedling_number)

                            !calculate probability of regeneration
                            probsum = 0.0
                            prob(is) = site%plots(ip)%seedling(is)*grow_cap
                            probsum  = probsum + prob(is)

						end do

						!reset counter
						site%plots(ip)%windCount = 0

					else

						!still counting down, decrement counter
						site%plots(ip)%windCount = max(0,                      &
							site%plots(ip)%windCount - 1)

						!set probsum to 0
						probsum=0.0

					end if
				end if

				!after setting seed and seedling banks

				!calculate cumulative probability of regeneration
				if (probsum .gt. epsilon(1.0)) then
					do is = 1, num_species
						prob(is) = prob(is)/probsum
					end do

					do is = 2, num_species
						prob(is) = prob(is - 1) + prob(is)
					end do
				else
					nrenew = 0
				end if

				!get number of trees on plot
				it = site%plots(ip)%numtrees

				if (nrenew .ge. 1) then

					!get new tree locations ---
					!count number of unfilled cells
					n_empty = count(site%plots(ip)%cells(:,:) .eq. 0)

					if (n_empty .gt. 0) then
						!some cells are unfilled - can place trees

						allocate(t_vals(n_empty))

						!loop through whole rows and columns and fill
						!r_empty and c_empty with empty cell indices -
						!keeping them together with the same 't' index

						t = 1
						do while (t .le. n_empty)
							do r = 1, maxcells
								do c = 1, maxcells
									if (site%plots(ip)%cells(r, c) .eq. 0) then
										r_empty(t) = r
										c_empty(t) = c
										t_vals(t) = t
										t = t + 1
									end if
								end do
							end do
						end do

						!shuffle t_vals array
						call shuffle(t_vals)

						!renew trees
						do irenew = 1, nrenew

							!determine species of new tree
							q0=urand()
							is = 1
							do while (q0 .gt. prob(is))
								is = is + 1
								if (is .gt. num_species) then
									is = 1 + int(urand(0.0, real(num_species)))
									q0	=	urand()
								endif
							end do

							!increment new tree counter
							new = new + 1

							!decrement seedling bank ofspecies
							site%plots(ip)%seedling(is) =                      &
								site%plots(ip)%seedling(is) - 1.0

							!increment number of trees and initialize new tree
							it = it + 1
							call initialize_tree(site%plots(ip)%trees(it),     &
								site%species(is), is)

							!get species index
							k = is

							!grab the r and c values of that index
							r = r_empty(t_vals(irenew))
							c = c_empty(t_vals(irenew))

							!set tree location to that value
							site%plots(ip)%trees(it)%row = r
							site%plots(ip)%trees(it)%col = c

							!set this to filled in cells array
							site%plots(ip)%cells(r,c) = 1


							!get dbh value of new tree (must be between
							!0.5 and 2.5 cm)
							z = 1.5 + nrand(0.0, 1.0)
							if (z .ge. 2.5) z = 2.5
							if (z .le. 0.5) z = 0.5
							site%plots(ip)%trees(it)%diam_bht = z

							!set canopy height
							site%plots(ip)%trees(it)%canopy_ht = 1.0

							!get height, canopy diameter, and biomass
							call forska_height(site%plots(ip)%trees(it))
							call stem_shape(site%plots(ip)%trees(it))
							call biomass_c(site%plots(ip)%trees(it))
							call biomass_n(site%plots(ip)%trees(it))
							call leaf_biomass_c(site%plots(ip)%trees(it))

							zz = lai_biomass_c(site%plots(ip)%trees(it))*      &
								site%species(k)%leafarea_c*2.0

							kh = int(site%plots(ip)%trees(it)%forska_ht)

							!calculate shading effect on tree
							if (site%plots(ip)%trees(it)%conifer) then
								canopy_shade = light_rsp(site%species(k),      &
									site%plots(ip)%con_light(kh))
							else
								canopy_shade = light_rsp(site%species(k),      &
									site%plots(ip)%dec_light(kh))
							end if

							!calculate environmental stressors
							call env_stress(site%plots(ip)%trees(it),          &
								canopy_shade, min_stress, stress_cat)
							site%plots(ip)%trees(it)%stressor = stress_cat

							if (site%plots(ip)%nutrient(k) .lt. min_stress) then
								min_stress = site%plots(ip)%nutrient(k)
								site%plots(ip)%trees(it)%stressor = 6
							end if

							if (perm_resp(k) .lt. min_stress) then
								site%plots(ip)%trees(it)%stressor = 4
								min_stress = perm_resp(k)
							end if

							if (site%species(k)%fc_flood .lt. min_stress) then
								site%plots(ip)%trees(it)%stressor = 5
							end if

							site%plots(ip)%trees(it)%env_resp(5) =             &
								site%species(k)%fc_flood
							site%plots(ip)%trees(it)%env_resp(6) =             &
								site%plots(ip)%nutrient(k)
							site%plots(ip)%trees(it)%env_resp(4) = perm_resp(k)

							!add leaf litter to soil and update N_used and NPP
							if (site%species(k)%conifer) then
								net_prim_prod = net_prim_prod + zz*leaf_b +    &
									site%plots(ip)%trees(it)%biomC

								N_used = N_used + zz/con_leaf_c_n +            &
                                    site%plots(ip)%trees(it)%biomN

								lc = site%species(k)%litter_class
								site%plots(ip)%soil%litter(lc) =               &
									site%plots(ip)%soil%litter(lc) +           &
                                    zz*(leaf_b - 1.0)/0.45
							else
								net_prim_prod = net_prim_prod +                &
                                    site%plots(ip)%trees(it)%biomC + zz

								N_used = N_used +                              &
                                    site%plots(ip)%trees(it)%biomN +           &
									zz/dec_leaf_c_n

								lc = site%species(k)%litter_class
								site%plots(ip)%soil%litter(lc) =               &
									site%plots(ip)%soil%litter(lc) + zz/0.45
							end if

						end do !end renewing trees

					end if !end if unfilled cells

				end if !end if renewing trees
				site%plots(ip)%numtrees = it

				!decrease seedling bank
				do is = 1, num_species
					site%plots(ip)%seedling(is) = site%plots(ip)%seedling(is)* &
						site%species(is)%seedling_lg/plotsize
				end do

			end if

			!update site and soil variables
			uconvert = hec_to_m2/plotsize
			N_used = N_used*uconvert
			net_prim_prod = net_prim_prod*uconvert
			site%plots(ip)%soil%net_prim_prodC =                               &
				site%plots(ip)%soil%net_prim_prodC + net_prim_prod
			site%plots(ip)%soil%net_prim_prodN =                               &
				site%plots(ip)%soil%net_prim_prodN + N_used

			!convert to tonnes/ha
			do lit = 1, 16
				if (lit .ne. 17) then
					ilitter = site%plots(ip)%soil%litter(lit)
					if (ilitter .gt. 0.00001) then
						site%plots(ip)%soil%litter(lit) =                      &
                            site%plots(ip)%soil%litter(lit)*uconvert
					else
						site%plots(ip)%soil%litter(lit) = 0.0
					end if
				end if
			end do

			!deallocate t_vals
			if (allocated(t_vals)) deallocate(t_vals)

		end do

	end subroutine Renewal

	!:.........................................................................:

end module Model
