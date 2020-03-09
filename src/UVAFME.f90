!:.............................................................................:

program UVAFME

  use Constants
  use Parameters
  use Input
  use Output
  use Species
  use Site
  use GenusGroups
  use Model


!*******************************************************************************
  !
  !This is the main UVAFME program which runs the model and prints out the main
  !banner and progress messages to the screen
  !Current version: UVAFME v. 3.0 October 2018
  !
!*******************************************************************************

!:.............................................................................:

	type(SiteData)                                   :: current_site
	type(Groups)                                     :: species_present
	type(SpeciesData), dimension(:),    allocatable  :: species_data
	real,              dimension(:, :), allocatable  :: all_site_vals
	real                                             :: start_time
	real                                             :: total_time
	integer                                          :: sndx, year
	integer                                          :: nargs
	integer                                          :: print_interval
	character(len = 80)                              :: filelist

	interface

     subroutine drawBanner(numsites, species_present)
		   use Parameters
		   use GenusGroups
		   implicit none
		   type(Groups)   ::  species_present
		   integer        ::  numsites
		 end subroutine drawBanner

		 subroutine showProgress(asite)
		   use Constants
		   use Parameters
		   use Site
		   implicit none
		   type(SiteData), intent(in) :: asite
		 end subroutine showProgress

	end interface

!:.............................................................................:

	!get command-line filelist name argument
	nargs = command_argument_count()
	if (nargs .ne. 1) then
		filelist = ''
	else
		call get_command_argument(1, filelist)
	endif

	!open and read file list & runtime file, initialize runtime parameters,
      !open all input files, and read and initialize litter parameters file
 	call initialize_inputFiles(filelist)

 	!read site list file and count how many sites to run
 	call read_sitelist(all_site_vals)
 	numsites = size(all_site_vals(:, 1))

	!read in species data file and initialize groups of genera and species
 	call read_speciesdata(species_data)
	call initialize_genus_groups(species_present, species_data)

 	!open output files and write headers
 	call initialize_outputFiles

 	!write runtime vars to screen
 	call drawBanner(numsites, species_present)

 	!start timing
 	call cpu_time(start_time)

 	do sndx = 1, numsites

 		!initialize current site
 		call initialize_site(current_site, all_site_vals(sndx, :), species_data)

        !make sure the site exists and has valid climate data
		if (current_site%site_id == invalid) then
 			write(*, *) '       No valid site or climate data for site ',      &
				current_site%site_name
 			write(*, *) '            Skipping site ', current_site%site_name
 			write(*, *) '                          '
			!free the memory used by this site
			call delete_site(current_site)
 			cycle
 		endif

 		!load site-specific random seed
 		call set_site_rng_seed(fixed_seed)

		!set accumulated climate to 0 for linear/GCM climate change
 		call set_site_climate

 		!skip this site if no species are present
		if ( size(current_site%species) .eq. 0 )  then
			write(*, *) '              No species present in site ',           &
				current_site%site_id
			write(*, *) '            Skipping site ', current_site%site_name
 			write(*, *) '             '
			!free the memory used by this site
			call delete_site(current_site)
 			cycle
 		endif

 		!print out current site and site vars to screen
 		call showProgress(current_site)

 		!run the model
 		do year = 0, numyears

 			!calculate current weather/site variables for the year
 			call BioGeoClimate(current_site, year)

 			!write climate and soil data
 			call write_site_data(current_site, year)
 			call write_soiln_data(current_site, year)

 			!calculate current LAI
 			call Canopy(current_site)

 			!run annual growth, mortality, and renewal
 			call Growth(current_site)

 			call Mortality(current_site, year)
 			call Renewal(current_site)

 			print_interval = year_print_interval

 			!print output
 			if ((mod(year, print_interval) .eq. 0) .or. year == numyears) then

 				!across-species attributes
 				call total_plot_values(current_site, maxcells, year)

 				!genus-and species-level attributes
				call write_genus_data(current_site, species_present, year)
 				call write_species_data(current_site, species_present, year)
 				call write_dead_species_data(current_site, species_present,    &
					year)
 				call write_dead_genus_data(current_site, species_present, year)

 				!tree-level attributes
 				if (tree_level_data) then
 					call write_tree_data(current_site, year)
 				endif
 			end if

 		end do

 		!free the memory used by this site
 		call delete_site(current_site)
 		write(sitef, *) 'Finished site ', current_site%site_id

 		!print current run time
 		call cpu_time(total_time)
 		write(*, *) '  Cumulative time : ', total_time - start_time
 		write(*, '(a80)')                                                      &
 '============================================================================='

 	end do

 	!close output files created
 	call close_outputFiles

end program UVAFME

!:.............................................................................:

subroutine drawBanner(numsites, species_present)
	!writes the UVAFME banner and runtime parameters currently in use
	!Inputs:
	!	numsites:         number of sites to run
	!	species_present:  species present in species file

	use Parameters
	use GenusGroups

	implicit none

		type(Groups), intent(in)  ::  species_present
		integer,      intent(in)  ::  numsites

		write(*, 500)                                                          &
'============================================================================='
		write(*, 500)                                                          &
'                       UVA Forest Model Enhanced                      '
		write(*, 500)                                                          &
'                           2018 Version 3.0                           '
		write(*, 500)                                                          &
'               Center For Regional Environmental Studies              '
		write(*, 500)                                                          &
'                        University of Virginia                        '
		write(*, 500)                                                          &
'                   Department of Environmental Sciences               '
		write(*, 500)                                                          &
'============================================================================='

		write(*, *) '  Running with parameters:'
		write(*, 400) 'Number of sites:', numsites
		write(*, 400) 'Number of years:', numyears
		write(*, 400) 'Number of plots:', numplots
		write(*, 400) 'Number of species:', species_present%numspecies
		write(*, 400) 'Maximum number of trees:', maxtrees
		write(*, 401) 'Plotsize:', plotsize

		if (with_clim_change) then

			write(*, *) 'Running with climate change'
			write(*, 400) 'Duration in years:', gcm_duration

			if (linear_cc) then

				write(*, *) 'Running with linear cc'

				if (incr_or_decr_temp .eq. 'incr' .and.                        &
					incr_or_decr_prcp .eq. 'incr') then
					write(*, 401) 'Total tmin increase', incr_tmin_by
					write(*, 401) 'Total tmax increase', incr_tmax_by
					write(*, 401) 'Total precip increase', incr_precip_by
				else if (incr_or_decr_temp .eq. 'decr' .and.                   &
					incr_or_decr_prcp .eq. 'decr') then
					write(*, 401) 'Total tmin decrease', decr_tmin_by
					write(*, 401) 'Total tmax decrease', decr_tmax_by
					write(*, 401) 'Total precip decrease', decr_precip_by
				else if (incr_or_decr_temp .eq. 'incr' .and.                   &
					incr_or_decr_prcp .eq. 'decr') then
					write(*, 401) 'Total tmin increase', incr_tmin_by
					write(*, 401) 'Total tmax increase', incr_tmax_by
					write(*, 401) 'Total precip decrease', decr_precip_by
				else if (incr_or_decr_temp .eq. 'decr' .and.                   &
					incr_or_decr_prcp .eq. 'incr') then
					write(*, 401) 'Total tmin decrease', decr_tmin_by
					write(*, 401) 'Total tmax decrease', decr_tmin_by
					write(*, 401) 'Total precip increase', incr_precip_by
				end if

			else if (use_gcm) then

				write(*, *) 'Using GCM data:'
				write(*, 400) 'GCM start year ', start_gcm
				write(*, 400) 'GCM end year ', end_gcm
			end if
		end if

		write(*, 400) 'Printing interval in years:', year_print_interval
		write(*, 500)                                                          &
'============================================================================='
		write(*, *)

400 	format(a30, i10)
401 	format(a30, f10.3)
402 	format(a30, a)
500 	format(a80)

end subroutine drawBanner

!:.............................................................................:

subroutine showProgress(asite)
	!prints out current site being run and site parameters
	!Inputs:
	!	asite: current site instance

	use Parameters
	use Site
	implicit none

		type(SiteData), intent(in)     ::  asite

		integer                        ::  num_site_species

		!get number of species at site
		num_site_species = size(asite%species)

		write(*, 500) 'Running for site ',  asite%site_id, asite%site_name
		write(*,*) '                    '
		write(*, 501) 'Number of species present: ', num_site_species

		!get altitude adjustment (if present), and print out
		if (asite%altitude .ne. rnvalid) then
			write(*, *) '             Site altitude adjustment ', asite%altitude
		endif

		write(*,*)

		!write all other site parameters
		write(*, 502) asite%elevation, asite%slope, asite%fire_prob,           &
		asite%wind_prob

500 	format(14x, a, i10, 4x, a)
501 	format(14x, a, i8)
502 	format(7x, 'Site parameters: elevation ', f9.3, '   slope     ', f7.3, &
        /  23x, ' fire/1000   ', f7.3, '   wind/1000 ', f7.3)

end subroutine showProgress

!:.............................................................................:
