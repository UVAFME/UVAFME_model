module Input
  use Constants
  use Parameters
  use IO
  use FileUtils
  use Species

  implicit none

!*******************************************************************************
  !
  !This module contains subroutines to open and read input files
  !and to initialize each site, plot, species, and tree list with
  !initial input parameters/variables.
  !
!*******************************************************************************

contains

	!:.........................................................................:

	subroutine initialize_inputFiles(filelist)
		!calls subroutines to open and read file_list.txt and runtime file,
		!initializes runtime parameters, and opens all other input files

		character(len = *)  :: filelist

		!read the file_list.txt file
		call read_file_list(filelist)

		!open runtime file
		call read_runFile

		!read runtime file and initialize runtime parameters
		call initialize_parameters

		!open all other input files
		call open_inputFiles

    !read and initialize litter parameters file
    call read_litterpars(litter_params)

	end subroutine initialize_inputFiles

	!:.........................................................................:

	subroutine initialize_parameters
		!initializes runtime parameters or sets to default values
		!Author: Katherine Holcomb, 2012 v. 1.0

		integer  ::   ios

		!list of runtime parameters
		namelist /uvafme/ numyears, numplots, plotsize, growth_thresh,         &
						maxtrees, maxcells, FRI, with_clim_change, linear_cc,  &
						use_gcm, start_gcm, end_gcm, plot_level_data,          &
						tree_level_data, incr_tmin_by, incr_tmax_by,           &
						year_print_interval, incr_precip_by, decr_tmin_by,     &
						decr_tmax_by, decr_precip_by, gcm_duration,            &
						incr_or_decr_prcp, incr_or_decr_temp, fixed_seed,      &
						debug

		!first initialize to default values
		debug = .false.
		year_print_interval = 10
		gcm_duration = 0
		use_gcm = .false.
		with_clim_change = .false.
		plot_level_data = .false.
		tree_level_data = .false.
		fixed_seed = .false.
		start_gcm = 0
		end_gcm = 100

		numyears = 1000
		numplots = 200
		maxtrees = 1000
		maxcells = 30
		plotsize = 500.0

		incr_tmin_by = 0.00
		incr_tmax_by = 0.00
		incr_precip_by = 0.0000
		decr_tmin_by = 0.00
		decr_tmax_by = 0.00
		decr_precip_by = 0.0000

		growth_thresh = 0.03

		incr_or_decr_prcp = 'decr'
		incr_or_decr_temp = 'incr'

		!read parameters namelist.  Any variables not set in this file will take
			!the default values above.
		read(rt_file, uvafme, iostat = ios)

		!create climate change variables
		if (with_clim_change) then
			if (.not. use_gcm) then

				linear_cc = .true.

				!check for linear cc inputs to make sure they make sense, if so,
					!create yearly changes in temperature/precipitation
				if (gcm_duration == 0) then
					stop "Inconsistent input for climate change: 0 years duration"
				else if (incr_tmin_by == 0.0 .and.                             &
					decr_tmax_by == 0.0 .and.                                  &
					decr_tmin_by == 0.0 .and.                                  &
					decr_tmax_by == 0.0 .and.                                  &
					incr_precip_by == 0.0 .and.                                &
					decr_precip_by == 0.0) then
					stop "Inconsistent climate change: no temp or precip changes"

				else if (incr_or_decr_temp == 'incr' .and.                     &
					incr_or_decr_prcp == 'incr') then

					if (incr_tmin_by <= 0.0 .and. incr_tmax_by <= 0.0          &
							.and. incr_precip_by <= 0.0) then
						stop "Inconsistent climate change: bad incr input"
					else
						tmin_change = incr_tmin_by/(gcm_duration + 1)
						tmax_change = incr_tmax_by/(gcm_duration + 1)
						precip_change = incr_precip_by/(gcm_duration + 1)
					end if

				else if (incr_or_decr_temp == 'decr' .and.                     &
						incr_or_decr_prcp == 'decr') then

					if (decr_tmin_by == 0.0 .and. decr_tmax_by == 0.0          &
							.and. decr_precip_by == 0.0) then
						stop "Inconsistent climate change: bad decr data"
					else

						!input for decr should be positive, assume user made
						  !mistake
						if (decr_tmin_by < 0.0) then
							write(*, *) "Assuming decrease intended"
							decr_tmin_by = -decr_tmin_by
						endif
						if (decr_tmax_by < 0.0) then
							write(*, *) "Assuming decrease intended"
							decr_tmax_by = -decr_tmax_by
						endif
						if (decr_precip_by < 0.0) then
							write(*,*) "Assuming decrease intended"
							decr_precip_by = -decr_precip_by
						endif

						tmin_change = -decr_tmin_by/(gcm_duration + 1)
						tmax_change = -decr_tmax_by/(gcm_duration + 1)
						precip_change = -decr_precip_by/(gcm_duration + 1)
					end if

				else if (incr_or_decr_temp == 'incr' .and.                     &
						incr_or_decr_prcp == 'decr') then
					if (incr_tmin_by <= 0.0 .and. incr_tmax_by <= 0.0          &
							.and. decr_precip_by == 0.0) then
						stop "Bad incr or decr data"
					else
						if (decr_precip_by < 0.0) then
							write(*,*) "Assuming decrease intended"
							decr_precip_by = -decr_precip_by
						end if
						tmin_change = incr_tmin_by/(gcm_duration + 1)
						tmax_change = incr_tmax_by/(gcm_duration + 1)
						precip_change = -decr_precip_by/(gcm_duration + 1)
					end if

				else if (incr_or_decr_temp == 'decr' .and.                     &
						incr_or_decr_prcp == 'incr') then
					if (decr_tmin_by == 0.0 .and. decr_tmax_by == 0.0          &
							.and. incr_precip_by == 0.0) then
						stop "Inconsistent climate change data"
					else
						if (decr_tmax_by < 0.0) then
							write(*,*) "Assuming decrease intended"
							decr_tmax_by = -decr_tmax_by
						end if
						if (decr_tmin_by < 0.0) then
							write(*,*) "Asuming decrease intended"
							decr_tmin_by = -decr_tmin_by
						end if
						tmin_change = -decr_tmin_by/(gcm_duration + 1)
						tmax_change = -decr_tmax_by/(gcm_duration + 1)
						precip_change = incr_precip_by/(gcm_duration + 1)
					end if
				end if

			else

				linear_cc = .false.
				tmin_change = 0.0
				tmax_change = 0.0
				precip_change = 0.0

				!for consistency, it needs to change one way or the other.
				!Note that with the GCM we will not add changes to the old data,
				!we will read in a new array for the tmin, tmax, and precip data

				if (gcm_duration == 0) then
					stop "Inconsistent input for climate change: 0 years duration"
				else
					!we start counting at 0
					end_gcm = start_gcm + gcm_duration - 1
				endif

			endif

		endif

		if (debug) then
			fixed_seed = .true.
			same_climate = .true.
		endif

		!close runtime file
		close(rt_file)

	end subroutine initialize_parameters

	!:.........................................................................:

	subroutine read_litterpars(litter_pars)
		!reads in parameters for litter decomposition
		!Author: Adrianna Foster, 2018 v. 1.0
		!Inputs/Outputs:
		!   litter_pars: litter parameters to be used in decomposition routine

		real, dimension(18, 10), intent(out) :: litter_pars

		character(len = MAX_NLEN)            :: litter_name, header
		integer                              :: lc, m, i
		integer                              :: ios

		litter_pars = rnvalid

		!we don't care about the header
		read(litterfile, '(a)'), header

		!read litter parameters and set to litter_pars array
		lc = 0
		do while (lc .le. 18)
			lc = lc + 1
			read(litterfile, *, iostat = ios, end = 10) litter_name,           &
				(litter_pars(lc, m), m = 2, 10)
		enddo

10  	continue

		!set name column to id
		do i = 1, 18
			litter_pars(i, 1) = i
		end do

		!check if any are invalid
		if (any(litter_pars .eq. rnvalid)) then
			write(logf, *) 'Error in litter parameters file, rnvalid'
		endif

		!close file
		close(litterfile)

	end subroutine read_litterpars

	!:.........................................................................:

	subroutine read_sitelist(site_vals)
		!reads sitelist file
		!Author: Katherine Holcomb, 2014 v. 1.0
		!Inputs/Outputs:
		!   site_vals: values from sitelist file

		real, allocatable, dimension(:,:), intent(inout)  :: site_vals
		integer                                           :: site_id
		integer                                           :: numoptions
		character(len = MAX_LINE)                         :: row
		character(len = :), allocatable, dimension(:)     :: fields
		character(len = :), allocatable, dimension(:)     :: rowvals

		integer                                           :: numsites
		integer                                           :: ival
		integer                                           :: ios
		integer                                           :: lc, l

		!count rows in sitelist file to get number of sites to run
		numsites = count_records(slist, 1)
		if (numsites <= 0) then
			call fatal_error("Error in site list file")
		endif

		!read header to obtain the number of site variables to modify
		read(slist, '(a)') row
		call split(row, fields, ',')
		numoptions = size(fields) - 1

		!allocate size of site_vals to number of sites x number of variables to
			!change
		allocate(site_vals(numsites, numoptions + 1))

		ios = 0
		lc = 0
		!do while doesn't really work quite as we'd like (we want repeat until)
		!also note that some compilers return nonzero iostat values for harmless
		  !errors, such as format conversions, so we can't rely on testing
		  !for ios=0

		!set all site_vals values (except first column) to -999
		site_vals(:,2:) = rnvalid
		do
			!read each row and increment line counter
			read(slist, '(a)', end = 10) row
			lc = lc + 1

			!split row by comma
			call split(row, rowvals, ',')

			!if nothing there, skip site
			if (size(rowvals) == 0) then
				cycle
			else if (len(rowvals(1)) == 0 .or. rowvals(1) == ' ') then
				cycle
			else
				!set first column to site_id and add to site_vals array
				read(rowvals(1), '(i9)') site_id
				site_vals(lc, 1) = site_id
				!print *, rowvals

				!grab all other site_vals data
				if (size(rowvals) > 1) then
					do l = 2, size(rowvals)
						if (len(rowvals(l)) == 0 .or. rowvals(l) == ' ') then
							cycle
						else
							if (scan(rowvals(l), '.') .eq. 0) then
								read(rowvals(l), '(i8)') ival
								site_vals(lc, l) = real(ival)
							else
								read(rowvals(l), '(f15.7)') site_vals(lc, l)
							endif
						endif
					enddo
				endif
			endif
		enddo

10  	continue

		write(logf, *) 'Site data initialized. Total read in: ', numsites

	end subroutine read_sitelist

	!:.........................................................................:

	subroutine read_rangelist(site_id,species_present)
		!reads in species rangelist for current site
		!Author: Katherine Holcomb, 2012 v. 1.0
		!Inputs/Outputs:
		!   site_id:         siteID of site, set to -999 if no data here
		!Outputs:
		!   species_present: species eligible for growth at site

		integer,                                       intent(inout) :: site_id
		character(len = 8), dimension(:), allocatable, intent(out)   :: species_present

		character(len = MAX_LONG_LINE)                               :: line
		character(len = :),        dimension(:),  allocatable        :: fields
		character(len = 8),        dimension(:),  allocatable        :: species_ids
		integer,                   dimension(:),  allocatable        :: species_pres
		real                                                         :: lat, long
		integer                                                      :: iwmo
		integer                                                      :: max_numspec
		integer                                                      :: lc, n
		integer                                                      :: ios

		!first line has the list of species.  Parse the line.
		read(rglist, '(a)', iostat = ios) line
		if (ios .eq. 0) then
			call split(line, fields, ',')
			max_numspec = size(fields) - 3
		else
			allocate(species_present(0))
			call warning("Problem in rangelist file.")
			return
		endif

		allocate(species_ids(max_numspec))
		allocate(species_pres(max_numspec))
		allocate(species_present(max_numspec))

		!set up the list of species
		do n = 1, max_numspec
			species_ids(n) = fields(n + 3)
		enddo

		species_present = 'NP'

		!now read the species lists, per site. Do while doesn't really work
		  !quite as we'd like (we want repeat until)
		lc = 0
		do
			read(rglist, *, iostat = ios, end = 10) iwmo, lat, long,           &
				species_pres
			lc = lc + 1
			if (site_id == iwmo ) then
				where (species_pres /= 0)
					species_present = species_ids
				endwhere
				exit
			endif
		end do
10  	continue
		rewind(rglist)

		if (all(species_present == 'NP')) then
			write(*, *) "No species present in site", site_id
			site_id = invalid
		else
			write(logf, *) 'Species data initialized for site ', site_id
		endif

	end subroutine read_rangelist

	!:.........................................................................:

	subroutine read_climate(site_id, tmin, tmax, prcp)
		!reads in monthly tmin, tmax, and precipitation data for site
		!Author: Katherine Holcomb, 2012 v. 1.0
		!Inputs/Outputs:
		!   site_id:  siteID of site, set to -999 if no data here
		!Outputs:
		!   tmin:     mean monthly minimum temperature (degrees C)
		!   tmax:     mean monthly maximum temperature (degrees C)
		!   prcp:     mean monthly precipitation (mm)

		integer,                  intent(inout) :: site_id
		real,  dimension(NTEMPS), intent(out)   :: tmin,tmax,prcp
		integer                                 :: siteid
		real                                    :: lat,long
		character(len = MAX_LONG_LINE)          :: header
		integer                                 :: ios
		integer                                 :: m

		tmin = rnvalid; tmax = rnvalid; prcp = rnvalid

		!we don't care about the header
		read(cfile, '(a)'), header

		do
			read(cfile, *, iostat = ios, end = 10) siteid, lat, long,          &
                        (tmin(m), m = 1, NTEMPS), (tmax(m), m = 1, NTEMPS),    &
                        (prcp(m), m = 1, NTEMPS)

			if (site_id .eq. siteid) then
				exit
			endif
		enddo

10  	continue
		rewind(cfile)

		if ((any(tmin .eq. rnvalid)) .or. (any(tmax .eq. rnvalid)) .or.        &
				(any(prcp .eq. rnvalid)) ) then
			write(logf, *) 'No climate data for site number ', site_id
			site_id = invalid
		endif

	end subroutine read_climate

	!:.........................................................................:

	subroutine read_climate_stds(site_id, tmin_std, tmax_std, prcp_std)
		!reads in monthly temp and precip standard deviations for site
		!Author: Katherine Holcomb, 2012 v. 1.0
		!Inputs/Outputs:
		!   site_id:  siteID of site, set to -999 if no data here
		!Outputs:
		!   tmin_std: sd of monthly minimum temperature (degrees C)
		!   tmax_std: sd of monthly maximum temperature (degrees C)
		!   prcp_std: sd of monthly precipitation (mm)

		integer,                 intent(inout) :: site_id
		real, dimension(NTEMPS), intent(out)   :: tmin_std, tmax_std
		real, dimension(NTEMPS), intent(out)   :: prcp_std

		integer                                :: siteid
		real                                   :: lat, long
		character(len = MAX_LONG_LINE)         :: header
		integer                                :: ios
		integer                                :: m

		!set standard deviations to invalid
		tmin_std = rnvalid; tmax_std = rnvalid; prcp_std = rnvalid

		!we don't care about the header
		read(cstdfile, '(a)'), header

		do
			read(cstdfile, *, iostat = ios, end = 10) siteid, lat, long,       &
				(tmin_std(m), m = 1, NTEMPS), (tmax_std(m), m = 1, NTEMPS),    &
				(prcp_std(m), m = 1, NTEMPS)

			if (site_id .eq. siteid) then
				exit
			endif
		enddo

10  	continue
		rewind(cstdfile)

		if ((any(tmin_std .eq. rnvalid)) .or.                                  &
			(any(tmax_std .eq. rnvalid)) .or.                                  &
			(any(prcp_std  .eq. rnvalid))) then
				write(logf, *) 'No climate stdev data for site number ', site_id
				site_id = invalid
		endif

	end subroutine read_climate_stds

	!:.........................................................................:

	subroutine read_extra_climate(site_id, cld)
		!reads cloudiness for site
		!Author: Adrianna Foster, 2018 v. 1.0
		!Inputs/Outputs:
		!   site_id:  siteID of site, set to -999 if no data here
		!Outputs:
		!   cld:     mean monthly cloudiness (%)


		integer,                 intent(inout) :: site_id
		real, dimension(NTEMPS), intent(out)   :: cld

		integer                                :: siteid
		real                                   :: lat, long
		character(len = MAX_LONG_LINE)         :: header
		integer                                :: ios
		integer                                :: m

		!set extra climate variables to -999
		cld = rnvalid

		!we don't care about the header
		read(cexfile, '(a)'), header

		do
			read(cexfile, *, iostat = ios,end = 10) siteid, lat, long,         &
								  (cld(m), m = 1, NTEMPS)
			if (site_id .eq. siteid) then
				exit
			end if
		end do
10  	continue
		rewind(cexfile)

		if(any(cld .eq. rnvalid)) then
			write(logf, *) 'No extra climate data for site number', site_id
			site_id = invalid
		end if
	end subroutine read_extra_climate

	!:.........................................................................:

	subroutine read_extra_climate_stds(site_id, cld_std)
		!reads in standard deviations for cloudiness
		!Author: Adrianna Foster, 2018 v. 1.0
		!Inputs/Outputs:
		!   site_id:  siteID of site, set to -999 if no data here
		!Outputs:
		!   cld_std:  sd of monthly cloudiness (%)

		integer,                 intent(inout) :: site_id
		real, dimension(NTEMPS), intent(out)   :: cld_std

		integer                                :: siteid
		real                                   :: lat, long
		character(len = MAX_LONG_LINE)         :: header
		integer                                :: ios
		integer                                :: m

		!set extra climate variables to -999
		cld_std = rnvalid

		!we don't care about the header
		read(cexstdfile, '(a)'), header

		do
			read(cexstdfile, *, iostat = ios, end = 10) siteid, lat, long,     &
				(cld_std(m), m = 1, NTEMPS)
			if (site_id .eq. siteid) then
				exit
			end if
		end do
10  	continue
		rewind(cexstdfile)

		if(any(cld_std  .eq. rnvalid)) then
			write(logf, *) 'No extra climate sttdev data for site', site_id
			site_id = invalid
		end if
	end subroutine read_extra_climate_stds

	!:.........................................................................:

	subroutine read_gcm_climate(site_id, ygcm, start_gcm, tmin, tmax, prcp)
		!reads in optional GCM climate change file
		!Author: Katherine Holcomb, 2012 v. 1.0
		!Inputs/Outputs:
		!   site_id:   siteID of site, set to -999 if no data here
		!   ygcm:        current 'gcm year' of simulation
		!   start_gcm: year to start gcm
		!
		!Outputs:
		!   tmin:      mean monthly minimum temperature (degrees C)
		!   tmax:      mean monthly maximum temperature (degrees C)
		!   prcp:      mean monthly precipitation (mm)

		integer,                 intent(inout) :: site_id
		integer,                 intent(in)    :: ygcm, start_gcm
		real, dimension(NTEMPS), intent(out)   :: tmin, tmax, prcp

		integer                                :: siteid
		real                                   :: lat, long, year
		character(len = MAX_LONG_LINE)         :: header
		integer                                :: m

		tmin = rnvalid; tmax = rnvalid; prcp = rnvalid

		if (ygcm .eq. start_gcm) then
			!new site.  We don't assume numerical ordering in the site
			!file.
			rewind(cgcmfile)

			!we don't care about the header
			read(cgcmfile, '(a)'), header

		endif

		do
			read(cgcmfile, *, end = 10) siteid, lat, long, year,               &
				(tmin(m), m = 1, NTEMPS), (tmax(m), m = 1, NTEMPS),            &
				(prcp(m), m = 1, NTEMPS)

			if (site_id .eq. siteid .and. int(year) .eq. ygcm) then
				exit
			endif
		enddo

10  	continue

		if (any(tmin .eq. rnvalid) .or. any(tmax .eq. rnvalid) .or.            &
				any(prcp .eq. rnvalid) ) then
				write(logf, *) 'No climate change data for site number ',      &
					site_id
			site_id = invalid
		endif

	end subroutine read_gcm_climate

	!:.........................................................................:

	subroutine read_site(site_id, sitename, siteregion, lat, long, elevation,  &
                        slope, aspect, asat, afc, itxt, fprob, wprob, gcm_year, &
                        wd)
        !reads in site and soil input data for site
        !Author: Katherine Holcomb, 2012 v. 1.0
		!Inputs/Outputs:
		!   site_id:   siteID of site, set to -999 if no data here
		!Outputs:
		!   sitename:   name of site (char)
		!   siteregion: region of site (char)
		!   lat:        latitude of site (degrees)
		!   long:       longitude of site (degrees)
		!   elevation:  elevation of site (m)
		!	slope:      slope of site (degrees)
		!   aspect:     aspect of site (degrees)
		!   asat:       mineral/A-layer saturation capacity (volumetric)
		!   afc:        mineral/A-layer field capacity (volumetric)
		!   itxt:       soil texture (1: coarse-grained; 2: fine-grained)
		!   fprob:      fires in 1000 years
		!   wprob:      windthrow events in 1000 years
		!   gcm_year:   year to start GCM file (usually based on stand age)

		character(len = MAX_NLEN), intent(out)   :: sitename, siteregion
		integer,                   intent(inout) :: site_id
		real,                      intent(out)   :: lat, long
		real,                      intent(out)   :: elevation, slope, aspect
		real,                      intent(out)   :: asat, afc, itxt
		real,                      intent(out)   :: fprob, wprob, gcm_year, wd

		integer                                  :: siteid
		character(MAX_LINE)                      :: header
		integer                                  :: ios
		real                                     :: wd1

		siteid = invalid

		!we don't care about the header
		read(sfile, '(a)'), header

		do
			read(sfile, *, iostat = ios, end = 10) siteid, lat, long,          &
				sitename, siteregion, elevation, slope, aspect, asat, afc,     &
				itxt, fprob, wprob, gcm_year, wd1, wd

			if (site_id .eq. siteid) then
				exit
			endif
		enddo

10  	continue
		rewind(sfile)

		if (siteid .eq. invalid) then
			call warning("Site not found")
			site_id = invalid
		endif

	end subroutine read_site

	!:.........................................................................:

	subroutine read_speciesdata(species_data)
		!reads in species parameters and initiliazes species list
		!Author: Katherine Holcomb, 2012 v. 1.0
		!Outputs:
		!   species_data:  species data type

		type(SpeciesData), dimension(:), allocatable, intent(out) :: species_data

		integer                                    :: ios, lc = 0
		character(len = MAX_LONG_LINE)             :: line
		character(len = MAX_NLEN)                  :: genus_name
		character(len = MAX_NLEN)                  :: taxonomic_name
		character(len = 8)                         :: unique_id
		character(len = MAX_NLEN)                  :: common_name
		integer                                    :: leaf_type
		integer                                    :: layer_type
		integer                                    :: genus_sort
		integer                                    :: genus_id
		integer                                    :: shade_tol
		integer                                    :: lownutr_tol
		integer                                    :: stress_tol
		integer                                    :: age_tol, dry_tol
		integer                                    :: flood_tol
		integer                                    :: perm_tol
		integer                                    :: org_tol, recr_age
		integer                                    :: fire_regen
		integer                                    :: litter_class
		real                                       :: max_age, max_diam
		real                                       :: max_ht
		real                                       :: wood_bulk_dens
		real                                       :: leafdiam_a
		real                                       :: leafarea_c
		real                                       :: deg_day_min
		real                                       :: deg_day_opt
		real                                       :: deg_day_max
		real                                       :: seed_surv
		real                                       :: seedling_lg
		real                                       :: invader
		real                                       :: bark_thick
		real                                       :: seed_num
		real                                       :: sprout_num
		real                                       :: arfa_0, g
		logical                                    :: conifer
		logical                                    :: layering
		integer                                    :: num_all_species

		!the species file has a header line
		num_all_species = count_records(splist, 1)
		if (num_all_species < 0) then
			call fatal_error('Error in species-list file')
		endif

		allocate(species_data(num_all_species))

		ios = 0
		!read and discard the first line
		read(splist, '(a500)', iostat = ios) line
		if (ios .ne. 0) then
			write(*,*) 'Unable to read the specieslist file'
			stop
		endif

		lc = 0
		!do while doesn't really work quite as we'd like
		!(we want repeat until)
		do while (lc .lt. num_all_species)
			lc = lc + 1
			read(splist, *, end = 100)                                         &
            genus_id,                                                          &
            genus_name,                                                        &
            genus_sort,                                                        &
            taxonomic_name,                                                    &
            common_name,                                                       &
            max_age,                                                           &
            max_diam,                                                          &
            max_ht,                                                            &
            arfa_0,                                                            &
            g,                                                                 &
            wood_bulk_dens,                                                    &
            leafdiam_a,                                                        &
            leafarea_c,                                                        &
            deg_day_min,                                                       &
            deg_day_opt,                                                       &
            deg_day_max,                                                       &
            shade_tol,                                                         &
            dry_tol,                                                           &
            flood_tol,                                                         &
            perm_tol,                                                          &
            lownutr_tol,                                                       &
            bark_thick,                                                        &
            fire_regen,                                                        &
            stress_tol,                                                        &
            age_tol,                                                           &
            leaf_type,                                                         &
            litter_class,                                                      &
            invader,                                                           &
            seed_num,                                                          &
            sprout_num,                                                        &
            layer_type,                                                        &
            org_tol,                                                           &
            recr_age,                                                          &
            seed_surv,                                                         &
            seedling_lg,                                                       &
            unique_id

			if (leaf_type == 0) then
				conifer = .false.
			else
				conifer = .true.
			endif

			if (layer_type == 0) then
				layering = .false.
			else
				layering = .true.
			end if

			call initialize_species(species_data(lc), lc, genus_name,          &
				taxonomic_name, unique_id, common_name, genus_id, shade_tol,   &
				lownutr_tol, stress_tol, age_tol, dry_tol, flood_tol,          &
				perm_tol, org_tol, bark_thick, fire_regen, max_age, max_diam,  &
				max_ht, wood_bulk_dens, rootdepth, leafdiam_a, leafarea_c,     &
				deg_day_min, deg_day_opt, deg_day_max, seedling_lg, invader,   &
				seed_num, sprout_num, layering, seed_surv, arfa_0, g, conifer, &
				litter_class, recr_age)

		enddo

100 	continue

		write(logf, *) 'Species data initialized. total read in: ',            &
			num_all_species

	end subroutine read_speciesdata

	!:.........................................................................:

end module Input
