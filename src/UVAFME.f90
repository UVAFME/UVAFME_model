program UVAFME

!*******************************************************************************
  !
  ! This is the main UVAFME program which runs the model and prints out the main
  ! banner and progress messages to the screen
  ! Current version: UVAFME v. 3.0 October 2018
  !
!*******************************************************************************

    use Constants
    use Parameters
    use Input
    use Output
    use Species
    use Site
    use GenusGroups
    use Model

    implicit none

    ! Data dictionary: constants
    character(len=8), parameter :: GENFIELD = 'genus'    ! Field name for genus output
    character(len=8), parameter :: SPECFIELD = 'species' ! Field name for species output

    ! Data dictonary: global variables
    type(SpeciesData), dimension(:),    allocatable  :: species_data    ! Array of species data objects
    real,              dimension(:, :), allocatable  :: all_site_vals   ! Array of site runtime parameters
    integer,           dimension(:, :), allocatable  :: species_ids     ! Array of species ids for each site
    type(SiteData)                                   :: current_site    ! Site object
    type(Groups)                                     :: species_present ! Species/genus names
    real                                             :: start_time      ! Start time of run
    real                                             :: total_time      ! Time for run
    integer                                          :: year            ! Year of simulation
    integer                                          :: nargs           ! Number of command-line arguments
    integer                                          :: sndx            ! Looping index
    integer                                          :: numsites        ! Number of sites to run
    character(len=80)                                :: filelist        ! File list file name
    interface

    subroutine drawBanner(numsites, species_present)
        !
        !  Writes the UVAFME banner and runtime parameters currently in use
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    07/26/12    K. Holcomb           Original Code
        !

        use Parameters
        use GenusGroups
        implicit none

        type(Groups) ::  species_present ! Species/genus names
        integer      ::  numsites        ! Number of sites to run

    end subroutine drawBanner

    subroutine showProgress(asite)
        !
        !  Prints out current site being run and some site parameters
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    07/26/12    K. Holcomb           Original Code
        !

        use Constants
        use Parameters
        use Site
        implicit none

        type(SiteData), intent(in) :: asite ! Site object

    end subroutine showProgress

    end interface

!:.............................................................................:

    ! Get command-line filelist name argument
    nargs = command_argument_count()
    if (nargs /= 1) then
        filelist = ''
    else
        call get_command_argument(1, filelist)
    endif

    ! Open and read file list & runtime file, initialize runtime parameters,
    ! open all other input files, and read and initialize litter parameters file
    call initialize_inputFiles(filelist)

    ! Read site list file and count how many sites to run
    call read_sitelist(all_site_vals)
    numsites = size(all_site_vals(:, 1))

    ! Read in species data file and initialize groups of genera and species
    call read_speciesdata(species_data)
    call initialize_genus_groups(species_present, species_data)

    ! Open output files and write headers
    call initialize_outputFiles

    ! Write runtime vars to screen
    call drawBanner(numsites, species_present)

    ! Start timing
    call cpu_time(start_time)

    ! Allocate species_ids to the number of site and species objects
    allocate(species_ids(numsites, size(species_data)))


    do sndx = 1, numsites

        ! Set the random number generator seed
        call set_site_rng_seed(fixed_seed)

        ! Initialize current site
        call initialize_site(current_site, all_site_vals(sndx, :),            &
            species_data, species_ids, sndx)

        ! Make sure the site exists and has valid climate data
        if (current_site%site_id == INVALID) then
            write(*, *) '       No valid site or climate data for site ',      &
                current_site%site_name
            write(*, *) '            Skipping site ', current_site%site_name
            write(*, *) '                          '

            ! Free the memory used by this site
            call delete_site(current_site)
            cycle
         endif

         ! Also skip this site if no species are present
        if ( size(current_site%species) .eq. 0 )  then
            write(*, *) '              No species present in site ',           &
                current_site%site_id
            write(*, *) '            Skipping site ', current_site%site_name
            write(*, *) '             '

            ! Free the memory used by this site
            call delete_site(current_site)
            cycle
         endif

        ! Print out current site and site vars to screen
        call showProgress(current_site)

        ! Run the model
        do year = 0, numyears

            ! Calculate current weather/site variables for the year
            call BioGeoClimate(current_site, year)

            ! Write climate and soil data
            if (mod(year, year_print_interval) == 0 .or. year == numyears) then
                call write_site_data(current_site, year)
                call write_soiln_data(current_site, year)
            end if

            ! Calculate current LAI and light environment
            call Canopy(current_site)

            ! Run annual growth, mortality, and renewal
            call Growth(current_site)
            call Mortality(current_site, year)
            call Renewal(current_site, year)

            ! Print output
            if (mod(year, year_print_interval) == 0 .or. year == numyears) then

                ! Across-species attributes
                call total_plot_values(current_site, year)

                ! Genus-and species-level attributes
                call write_genus_or_species_data(current_site,                 &
                    species_present, year, SPECFIELD,                          &
                    species_present%numspecies, biom_by_s, pl_biom_by_s)

                call write_genus_or_species_data(current_site,                 &
                    species_present, year, GENFIELD,                           &
                    species_present%numgenera, biom_by_g, pl_biom_by_g)

                call write_dead_genus_or_species_data(current_site,            &
                    species_present, year, SPECFIELD,                          &
                    species_present%numspecies, dead_s, dead_ps)

                call write_dead_genus_or_species_data(current_site,            &
                    species_present, year, GENFIELD,                           &
                    species_present%numgenera, dead_g, dead_pg)

                ! Tree level data
                if (tree_level_data) then
                    call write_tree_data(current_site, year)
                end if

            end if

        end do

        ! Free the memory used by this site and tell user we finished the
        ! site
        call delete_site(current_site)
        write(sitef, *) 'Finished site ', current_site%site_id

        ! Print current run time
        call cpu_time(total_time)
        write(*, *) '  Cumulative time : ', total_time - start_time
        write(*, '(A80)')                                                      &
 '============================================================================='

    end do

    ! Close output files created
    call close_outputFiles

end program UVAFME

!:.............................................................................:

subroutine drawBanner(numsites, species_present)
    !
    !  Writes the UVAFME banner and runtime parameters currently in use
    !
    !  Record of revisions:
    !      Date       Programmer          Description of change
    !      ====       ==========          =====================
    !    07/26/12    K. Holcomb           Original Code
    !

    use Parameters
    use GenusGroups
    implicit none

    ! Data dictionary: calling arguments
    type(Groups), intent(in) :: species_present ! Species/genus names
    integer,      intent(in) :: numsites        ! Number of sites to run


    write(*, 500)                                                              &
'============================================================================='
    write(*, 500)                                                              &
'                       UVA Forest Model Enhanced                      '
    write(*, 500)                                                              &
'                           2018 Version 3.0                           '
    write(*, 500)                                                              &
'============================================================================='

    write(*, *) '  Running with parameters:'
    write(*, 400) 'Number of sites:', numsites
    write(*, 400) 'Number of years:', numyears
    write(*, 400) 'Number of plots:', numplots
    write(*, 400) 'Number of species:', species_present%numspecies
    write(*, 400) 'Maximum number of trees:', maxtrees
    write(*, 400) 'Maximum number of cells:', maxcells*maxcells
    write(*, 401) 'Plotsize:', plotsize

    if (with_clim_change) then

        write(*, *) 'Running with climate change'
        write(*, 400) 'Duration in years:', gcm_duration

        if (linear_cc) then

            write(*, *) 'Running with linear cc'

            if (incr_or_decr_temp .eq. 'incr' .and.                            &
                incr_or_decr_prcp .eq. 'incr') then
                write(*, 401) 'Total tmin increase', incr_tmin_by
                write(*, 401) 'Total tmax increase', incr_tmax_by
                write(*, 401) 'Total precip increase', incr_precip_by
            else if (incr_or_decr_temp .eq. 'decr' .and.                       &
                incr_or_decr_prcp .eq. 'decr') then
                write(*, 401) 'Total tmin decrease', decr_tmin_by
                write(*, 401) 'Total tmax decrease', decr_tmax_by
                write(*, 401) 'Total precip decrease', decr_precip_by
            else if (incr_or_decr_temp .eq. 'incr' .and.                       &
                incr_or_decr_prcp .eq. 'decr') then
                write(*, 401) 'Total tmin increase', incr_tmin_by
                write(*, 401) 'Total tmax increase', incr_tmax_by
                write(*, 401) 'Total precip decrease', decr_precip_by
            else if (incr_or_decr_temp .eq. 'decr' .and.                       &
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
    write(*, 500)                                                              &
'============================================================================='
    write(*, *)

    400 format(A30, I10)
    401 format(A30, F10.3)
    402 format(A30, A)
    500 format(A80)

end subroutine drawBanner

!:.............................................................................:

subroutine showProgress(asite)
    !
    !  Prints out current site being run and some site parameters
    !
    !  Record of revisions:
    !      Date       Programmer          Description of change
    !      ====       ==========          =====================
    !    07/26/12    K. Holcomb           Original Code
    !

    use Parameters
    use Site
    implicit none

    ! Data dictionary: calling arguments
    type(SiteData), intent(in) :: asite ! Site object

    ! Data dictionary: local variables
    integer ::  num_site_species ! Number of species present at site

    ! Get number of species at site
    num_site_species = size(asite%species)

    write(*, 500) 'Running for site ',  asite%site_id, asite%site_name
    write(*,*) '                    '
    write(*, 501) 'Number of species present: ', num_site_species

    ! Get altitude adjustment (if present), and print out
    if (adjust_altitude .and. asite%altitude .ne. RNVALID) then
        write(*, *) '             Site altitude adjustment ', asite%altitude
    endif

    write(*,*)

    ! Write some other site parameters
    write(*, 502) asite%elevation, asite%slope, asite%aspect, asite%fire_prob

    500 format(14X, A, I10, 4X, a)
    501 format(14X, A, I8)
    502 format(7X, 'Site parameters: elevation ', F9.3, '   slope     ', F7.3, &
        /  23X, ' aspect      ', F7.3, '   fire/1000 ', F7.3)

end subroutine showProgress

!:.............................................................................:
