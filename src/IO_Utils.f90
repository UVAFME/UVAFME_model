module IO

!*******************************************************************************
  !
  ! This module contains utility functions to open and close I/O files.
  !
!*******************************************************************************

  use Constants
  use Parameters
  use FileUtils
  use csv_file

  implicit none


contains

  !:...........................................................................:

    subroutine read_file_list(filelist)
        !
        !  Read the names of the directories for the input and output files.
        !  If the filelist file is not present, we set them to default values.
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *), optional, intent(in) :: filelist ! File list filename

        ! Data dictonary: local variables
        character(len = MAX_FILE) :: input_directory    ! Input directory name
        character(len = MAX_FILE) :: output_directory   ! Output directory name
        character(len = MAX_DIR)  :: climate_directory  ! Climate input directory name
        character(len = MAX_DIR)  :: GCM_directory      ! Climate change input directory name
        character(len = MAX_DIR)  :: site_directory     ! Site input directory name
        character(len = MAX_DIR)  :: sitelist_directory ! Sitelist input directory name
        character(len = MAX_DIR)  :: rt_directory       ! Runtime file input directory name
        character(len = MAX_DIR)  :: speclist_directory ! Input species file directory name
        integer                   :: nunit              ! File unit for filelist file
        integer                   :: ios                ! I/O status
        logical                   :: file_exists        ! Does the file exist?
        character(len = MAX_CHAR) :: message            ! Error message

        ! Directory namelist
        namelist /filenames/ input_directory, output_directory,                &
                            climate_directory, site_directory,                 &
                            sitelist_directory, GCM_directory,                 &
                            speclist_directory, rt_directory

        ! Initialize to default values
        input_directory = 'input_data'
        output_directory = 'output_data'
        climate_directory = 'input_data'
        GCM_directory = 'input_data'
        site_directory = 'input_data'
        sitelist_directory = 'input_data'
        rt_directory = 'input_data'
        speclist_directory = 'input_data'

        ! If filelist present, check that it exists
        if (present(filelist)) then

            inquire(file = filelist, exist = file_exists)

            if (.not. file_exists) then
                write(message, '(A, A)') "Can't find specified file list: ",   &
                    filelist
               call fatal_error(message)
            else
                ! File list exists - get the unit number
                nunit = unit_number()
                open(nunit, file = filelist, status = 'unknown',               &
                    iostat = ios)
                if (ios == 0) then
                    read(nunit, filenames, iostat = ios)
                    if (ios /= 0) then
                        write(message, '(A, A)') "Error reading ", filelist
                        call fatal_error(message)
                    end if
                endif
                close(nunit)
            end if
        else
            write(message, '(A)') "File list not specified - using defaults."
            call warning(message)
        endif

        ! Trim directory names of whitespace and set to global values
        inputdir = trim(input_directory)
        outputdir = trim(output_directory)
        climatedir = trim(climate_directory)
        GCMdir = trim(GCM_directory)
        sitedir = trim(site_directory)
        slistdir = trim(sitelist_directory)
        rtdir = trim(rt_directory)
        splistdir = trim(speclist_directory)

    end subroutine read_file_list

    !:.........................................................................:

    subroutine read_runFile
        !
        !  Opens the runtime file if it exists
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: local variables
        character(len = MAX_FILE) :: filename ! Runtime file name
        character(len = MAX_PATH) :: pathname ! Full path to runtime file

        ! Get the current working directory
        call get_cwd(fqcwd)

        ! Set file prefix
        write(CODE_ID, '(a,a)') CODENAME, VERSION_ID

        ! Open runtime file
        call build_filename(CODE_ID, '_runtime.txt', filename)
        call build_pathname(rtdir, filename, pathname)
        rt_file = open_file(pathname, 'r')

        ! If can't find runtime file, set to defaults
        if (rt_file == INVALID) then
            call warning('Unable to open runtime (parameter) file.')
            call warning('Using all default values.')
        endif

    end subroutine read_runFile

    !:.........................................................................:

    subroutine open_inputFiles
        !
        !  Opens all other input files and checks to make sure they aren't empty
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !    01/25/21     A. C. Foster         Updated to use helper function to
        !                                       reduce copy pasta

        ! Data dictionary: local variables
        character(len = MAX_PATH) :: pathname ! Full path to file

        ! Get current working directory
        call get_cwd(fqcwd)

        ! Logfile (must be opened as an "input" file so other input can write
        ! to it)
        call build_pathname(outputdir, 'log.txt', pathname)
        logf = open_file(pathname)

        ! Site list
        call open_and_check('_sitelist.csv', slistdir, slist, 1, 'site-list')

        ! Species attributes file
        call open_and_check('_specieslist.csv', splistdir, splist, 1,          &
            'species-list')

        ! Site file
        call open_and_check('_site.csv', sitedir, sfile, 1, 'site data')

        ! Climate file
        call open_and_check('_climate.csv', climatedir, cfile, 1, 'climate')

        ! Extra climate file (cld, rh, wind)
        call open_and_check('_climate_ex.csv', climatedir, cexfile, 1,         &
            'extra climate')

        ! Lightning file
        call open_and_check('_lightning.csv', climatedir, lightfile, 1,        &
            'lightning')

        ! Litter parameters
        call open_and_check('_litterpars.csv', sitedir, litterfile, 1,         &
            'litter')


        if (use_climstd) then
            ! Climate standard deviation
            call open_and_check('_climate_stddev.csv', climatedir, cstdfile,   &
                1, 'climate stddev')

            ! Extra climate standard deviation
            call open_and_check('_climate_ex_stddev.csv', climatedir,          &
                cexstdfile, 1, 'extra climate stddev')
        end if

        ! Climate change data from input file
        if (use_gcm) then
            ! Climate input
            call open_and_check('_climate_GCM.csv', GCMdir, cgcmfile, 1, 'GCM')

            ! Lightning climate change
            call open_and_check('_lightning_GCM.csv', GCMdir, cgclightfile,    &
                1, 'lightning GCM')
        endif

        ! Rangelist file
        if (use_rangelist) then
            call open_and_check('_rangelist.csv', inputdir, rglist, 1, 'range')
        else
            call warning('No rangelist file specified.')
            call warning('Using all species in all sites.')
        endif

    end subroutine open_inputFiles

    !:.........................................................................:

    subroutine open_and_check(suffix, dir, fileunit, nheaders, name)
        !
        !  Opens a file and checks to make sure it isn't empty. Helper
        !   subroutine for open_inputFiles
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/25/21     A. C. Foster        Original Code
        !

        ! Data dictionary - calling arguments
        character(len = *),  intent(in)    :: suffix   ! File suffix
        character(len = *),  intent(in)    :: dir      ! File directory
        integer,             intent(inout) :: fileunit ! Unit number of file
        integer,             intent(in)    :: nheaders ! Number of headers to skip in counting rows
        character(len = *),  intent(in)    :: name     ! File name for error message


        ! Data dictonary: local variables
        character(len = MAX_FILE) :: filename ! File name
        character(len = MAX_PATH) :: pathname ! Full path to file
        character(len = MAX_CHAR) :: message  ! Error message

        ! Build file name and path name
        call build_filename(CODE_ID, suffix, filename)
        call build_pathname(dir, filename, pathname)

        ! Open file
        fileunit = open_file(pathname, 'r')
        if (fileunit == INVALID) then
            write(message, '(A,A)') "Unable to find ", name, " file."
            call fatal_error(message)
        endif

        ! Count records
        if (count_records(fileunit, nheaders) <= 0) then
            write(message, '(A,A)') "Error in ", name, " file."
            call fatal_error(message)
        endif

    end subroutine open_and_check

    !:.........................................................................:

    subroutine open_outputFiles
        !
        !  Opens all output files
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code

        ! Data dictionary: local variables
        character(len = MAX_PATH) :: pathname ! Full pathname for file

        ! Open files

        ! Site log
        call build_pathname(outputdir, 'site_log.txt', pathname)
        sitef = open_file(pathname)

        ! Soil decomposition
        call build_pathname(outputdir, 'SoilDecomp.csv', pathname)
        soildecomp = open_file(pathname)

        ! Climate
        call build_pathname(outputdir, 'Climate.csv', pathname)
        clim_unit = open_file(pathname)

        ! Genus-level data (live)
        call build_pathname(outputdir, 'Genus_Data.csv', pathname)
        biom_by_g = open_file(pathname)

        ! Species-level data (live)
        call build_pathname(outputdir, 'Species_Data.csv', pathname)
        biom_by_s = open_file(pathname)

        ! Species-level data (dead)
        call build_pathname(outputdir, 'Dead_Species_Data.csv', pathname)
        dead_s = open_file(pathname)

        ! Genus-level data (dead)
        call build_pathname(outputdir, 'Dead_Genus_Data.csv', pathname)
        dead_g = open_file(pathname)

        ! Across-species data
        call build_pathname(outputdir, 'Total_Plot_Values.csv', pathname)
        plotvals = open_file(pathname)

        ! Plotlevel data
        if (plot_level_data) then

            ! Genus-level data (live)
            call build_pathname(outputdir, 'Plot_Genus_Data.csv', pathname)
            pl_biom_by_g = open_file(pathname)

            ! Genus-level data (dead)
            call build_pathname(outputdir, 'Plot_Dead_Genus.csv', pathname)
            dead_pg = open_file(pathname)

            ! Species-level data (live)
            call build_pathname(outputdir, 'Plot_Species_Data.csv', pathname)
            pl_biom_by_s = open_file(pathname)

            ! Species-level data (dead)
            call build_pathname(outputdir, 'Plot_Dead_Species.csv', pathname)
            dead_ps = open_file(pathname)
        end if

        ! Treelevel data
        if (tree_level_data) then
            call build_pathname(outputdir, 'Plot_Tree_Data.csv', pathname)
            pl_tree = open_file(pathname)
        end if

        ! Testing files
        if (testing) then
            call build_pathname(outputdir, 'Fuel_Conds.csv', pathname)
            fuel_conds = open_file(pathname)

            call build_pathname(outputdir, 'Fire_Conds.csv', pathname)
            dayfire = open_file(pathname)

            call build_pathname(outputdir, 'Cons_Data.csv', pathname)
            cons_out = open_file(pathname)

        end if

        if (reg_testing) then

            call build_pathname(outputdir, 'Spec_Regen.csv', pathname)
            spec_regen = open_file(pathname)

        end if

    end subroutine open_outputFiles

    !:.........................................................................:

    subroutine write_headers()
        !
        !  Writes headers for all output files
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code

        ! Write file headers

        ! Soil decomposition
        call csv_write(soildecomp, 'siteID', .false.)
        call csv_write(soildecomp, 'runID', .false.)
        call csv_write(soildecomp, 'year', .false.)
        call csv_write(soildecomp, 'odepth', .false.)
        call csv_write(soildecomp, 'odepth_sd', .false.)
        call csv_write(soildecomp, 'mdepth', .false.)
        call csv_write(soildecomp, 'moss_biom', .false.)
        call csv_write(soildecomp, 'active', .false.)
        call csv_write(soildecomp, 'OM', .false.)
        call csv_write(soildecomp, 'OM_N', .false.)
        call csv_write(soildecomp, 'flit_cornus', .false.)
        call csv_write(soildecomp, 'flit_acerfrax', .false.)
        call csv_write(soildecomp, 'flit_prunus', .false.)
        call csv_write(soildecomp, 'flit_betula', .false.)
        call csv_write(soildecomp, 'flit_queralba', .false.)
        call csv_write(soildecomp, 'flit_tsugthuj', .false.)
        call csv_write(soildecomp, 'flit_populus', .false.)
        call csv_write(soildecomp, 'flit_fagus', .false.)
        call csv_write(soildecomp, 'flit_querrubr', .false.)
        call csv_write(soildecomp, 'flit_abies', .false.)
        call csv_write(soildecomp, 'flit_picea', .false.)
        call csv_write(soildecomp, 'flit_pinus', .false.)
        call csv_write(soildecomp, 'flit_roots', .false.)
        call csv_write(soildecomp, 'flit_smboles', .false.)
        call csv_write(soildecomp, 'flit_lboles', .false.)
        call csv_write(soildecomp, 'flit_twigs', .false.)
        call csv_write(soildecomp, 'flit_smbranch', .false.)
        call csv_write(soildecomp, 'flit_lbranch', .false.)
        call csv_write(soildecomp, 'flit_WDW', .false.)
        call csv_write(soildecomp, 'flit_moss', .false.)
        call csv_write(soildecomp, 'lit_cornus', .false.)
        call csv_write(soildecomp, 'lit_acerfrax', .false.)
        call csv_write(soildecomp, 'lit_prunus', .false.)
        call csv_write(soildecomp, 'lit_betula', .false.)
        call csv_write(soildecomp, 'lit_queralba', .false.)
        call csv_write(soildecomp, 'lit_tsugthuj', .false.)
        call csv_write(soildecomp, 'lit_populus', .false.)
        call csv_write(soildecomp, 'lit_fagus', .false.)
        call csv_write(soildecomp, 'lit_querrubr', .false.)
        call csv_write(soildecomp, 'lit_abies', .false.)
        call csv_write(soildecomp, 'lit_picea', .false.)
        call csv_write(soildecomp, 'lit_pinus', .false.)
        call csv_write(soildecomp, 'lit_roots', .false.)
        call csv_write(soildecomp, 'lit_smboles', .false.)
        call csv_write(soildecomp, 'lit_lboles', .false.)
        call csv_write(soildecomp, 'lit_twigs', .false.)
        call csv_write(soildecomp, 'lit_smbranch', .false.)
        call csv_write(soildecomp, 'lit_lbranch', .false.)
        call csv_write(soildecomp, 'lit_WDW', .false.)
        call csv_write(soildecomp, 'lit_moss', .false.)
        call csv_write(soildecomp, 'avail_n', .true.)

        ! Climate
        call csv_write(clim_unit, 'siteID', .false.)
        call csv_write(clim_unit, 'runID',  .false.)
        call csv_write(clim_unit, 'year', .false.)
        call csv_write(clim_unit, 'rain', .false.)
        call csv_write(clim_unit, 'pet', .false.)
        call csv_write(clim_unit, 'solar_rad', .false.)
        call csv_write(clim_unit, 'thaw_depth', .false.)
        call csv_write(clim_unit, 'organic_depth', .false.)
        call csv_write(clim_unit, 'avail_n', .false.)
        call csv_write(clim_unit, 'aet', .false.)
        call csv_write(clim_unit, 'grow', .false.)
        call csv_write(clim_unit, 'pc_germ', .false.)
        call csv_write(clim_unit, 'degd', .false.)
        call csv_write(clim_unit, 'drydays', .false.)
        call csv_write(clim_unit, 'saw0_ByFC', .false.)
        call csv_write(clim_unit, 'saw0_BySAT', .false.)
        call csv_write(clim_unit, 'aow0_ByMin', .false.)
        call csv_write(clim_unit, 'wilt_days', .false.)
        call csv_write(clim_unit, 'flood_d', .true.)

        ! Genus-level data (live)
        call csv_write(biom_by_g, 'siteID', .false.)
        call csv_write(biom_by_g, 'runID',  .false.)
        call csv_write(biom_by_g, 'year', .false.)
        call csv_write(biom_by_g, 'genus', .false.)
        call csv_write(biom_by_g, '0 to 5', .false.)
        call csv_write(biom_by_g, '5 to 10', .false.)
        call csv_write(biom_by_g, '10 to 20', .false.)
        call csv_write(biom_by_g, '20 to 30', .false.)
        call csv_write(biom_by_g, '30 to 40', .false.)
        call csv_write(biom_by_g, '40 to 50', .false.)
        call csv_write(biom_by_g, '50 to 60', .false.)
        call csv_write(biom_by_g, '60 to 70', .false.)
        call csv_write(biom_by_g, '70 to 80', .false.)
        call csv_write(biom_by_g, '80 to 90', .false.)
        call csv_write(biom_by_g, 'over 90', .false.)
        call csv_write(biom_by_g, '0 to 5 biom', .false.)
        call csv_write(biom_by_g, '5 to 10 biom', .false.)
        call csv_write(biom_by_g, '10 to 20 biom', .false.)
        call csv_write(biom_by_g, '20 to 30 biom', .false.)
        call csv_write(biom_by_g, '30 to 40 biom', .false.)
        call csv_write(biom_by_g, '40 to 50 biom', .false.)
        call csv_write(biom_by_g, '50 to 60 biom', .false.)
        call csv_write(biom_by_g, '60 to 70 biom', .false.)
        call csv_write(biom_by_g, '70 to 80 biom', .false.)
        call csv_write(biom_by_g, '80 to 90 biom', .false.)
        call csv_write(biom_by_g, 'over 90 biom', .false.)
        call csv_write(biom_by_g, 'degday_resp', .false.)
        call csv_write(biom_by_g, 'drought_resp', .false.)
        call csv_write(biom_by_g, 'shade_resp', .false.)
        call csv_write(biom_by_g, 'perm_resp', .false.)
        call csv_write(biom_by_g, 'flood_resp', .false.)
        call csv_write(biom_by_g, 'nutrient_resp', .false.)
        call csv_write(biom_by_g, 'max_diam', .false.)
        call csv_write(biom_by_g, 'mean_diam', .false.)
        call csv_write(biom_by_g, 'mean_age', .false.)
        call csv_write(biom_by_g, 'max_hgt', .false.)
        call csv_write(biom_by_g, 'leaf_area_ind', .false.)
        call csv_write(biom_by_g, 'basal_area', .false.)
        call csv_write(biom_by_g, 'basal_sd', .false.)
        call csv_write(biom_by_g, 'total_biomC', .false.)
        call csv_write(biom_by_g, 'total_biomC_sd', .false.)
        call csv_write(biom_by_g, 'biomC_lg', .false.)
        call csv_write(biom_by_g, 'biomC_std_lg', .false.)
        call csv_write(biom_by_g, 'biomC_sm', .false.)
        call csv_write(biom_by_g, 'biomC_std_sm', .false.)
        call csv_write(biom_by_g, 'basal_lg', .false.)
        call csv_write(biom_by_g, 'basal_std_lg', .false.)
        call csv_write(biom_by_g, 'basal_sm', .false.)
        call csv_write(biom_by_g, 'basal_std_sm', .false.)
        call csv_write(biom_by_g, 'dens_lg', .false.)
        call csv_write(biom_by_g, 'dens_std_lg', .false.)
        call csv_write(biom_by_g, 'dens_sm', .false.)
        call csv_write(biom_by_g, 'dens_std_sm', .false.)
        call csv_write(biom_by_g, 'dbh_lg', .false.)
        call csv_write(biom_by_g, 'dbh_std_lg', .false.)
        call csv_write(biom_by_g, 'dbh_sm', .false.)
        call csv_write(biom_by_g, 'dbh_std_sm', .true.)

        ! Species-level data (live)
        call csv_write(biom_by_s,'siteID', .false.)
        call csv_write(biom_by_s, 'runID', .false.)
        call csv_write(biom_by_s, 'year', .false.)
        call csv_write(biom_by_s, 'genus', .false.)
        call csv_write(biom_by_s, 'species', .false.)
        call csv_write(biom_by_s, '0 to 5', .false.)
        call csv_write(biom_by_s, '5 to 10', .false.)
        call csv_write(biom_by_s, '10 to 20', .false.)
        call csv_write(biom_by_s, '20 to 30', .false.)
        call csv_write(biom_by_s, '30 to 40', .false.)
        call csv_write(biom_by_s, '40 to 50', .false.)
        call csv_write(biom_by_s, '50 to 60', .false.)
        call csv_write(biom_by_s, '60 to 70', .false.)
        call csv_write(biom_by_s, '70 to 80', .false.)
        call csv_write(biom_by_s, '80 to 90', .false.)
        call csv_write(biom_by_s, 'over 90', .false.)
        call csv_write(biom_by_s, '0 to 5 biom', .false.)
        call csv_write(biom_by_s, '5 to 10 biom', .false.)
        call csv_write(biom_by_s, '10 to 20 biom', .false.)
        call csv_write(biom_by_s, '20 to 30 biom', .false.)
        call csv_write(biom_by_s, '30 to 40 biom', .false.)
        call csv_write(biom_by_s, '40 to 50 biom', .false.)
        call csv_write(biom_by_s, '50 to 60 biom', .false.)
        call csv_write(biom_by_s, '60 to 70 biom', .false.)
        call csv_write(biom_by_s, '70 to 80 biom', .false.)
        call csv_write(biom_by_s, '80 to 90 biom', .false.)
        call csv_write(biom_by_s, 'over 90 biom', .false.)
        call csv_write(biom_by_s, 'degday_resp', .false.)
        call csv_write(biom_by_s, 'drought_resp', .false.)
        call csv_write(biom_by_s, 'shade_resp', .false.)
        call csv_write(biom_by_s, 'perm_resp', .false.)
        call csv_write(biom_by_s, 'flood_resp', .false.)
        call csv_write(biom_by_s, 'nutrient_resp', .false.)
        call csv_write(biom_by_s, 'max_diam', .false.)
        call csv_write(biom_by_s, 'mean_diam', .false.)
        call csv_write(biom_by_s, 'mean_age', .false.)
        call csv_write(biom_by_s, 'max_hgt', .false.)
        call csv_write(biom_by_s, 'leaf_area_ind', .false.)
        call csv_write(biom_by_s, 'basal_area', .false.)
        call csv_write(biom_by_s, 'basal_sd', .false.)
        call csv_write(biom_by_s, 'total_biomC', .false.)
        call csv_write(biom_by_s, 'total_biomC_sd', .false.)
        call csv_write(biom_by_s, 'biomC_lg', .false.)
        call csv_write(biom_by_s, 'biomC_std_lg', .false.)
        call csv_write(biom_by_s, 'biomC_sm', .false.)
        call csv_write(biom_by_s, 'biomC_std_sm', .false.)
        call csv_write(biom_by_s, 'basal_lg', .false.)
        call csv_write(biom_by_s, 'basal_std_lg', .false.)
        call csv_write(biom_by_s, 'basal_sm', .false.)
        call csv_write(biom_by_s, 'basal_std_sm', .false.)
        call csv_write(biom_by_s, 'dens_lg', .false.)
        call csv_write(biom_by_s, 'dens_std_lg', .false.)
        call csv_write(biom_by_s, 'dens_sm', .false.)
        call csv_write(biom_by_s, 'dens_std_sm', .false.)
        call csv_write(biom_by_s, 'dbh_lg', .false.)
        call csv_write(biom_by_s, 'dbh_std_lg', .false.)
        call csv_write(biom_by_s, 'dbh_sm', .false.)
        call csv_write(biom_by_s, 'dbh_std_sm', .true.)

        ! Species-level data (dead)
        call csv_write(dead_s, 'siteID', .false.)
        call csv_write(dead_s, 'runID',  .false.)
        call csv_write(dead_s, 'year', .false.)
        call csv_write(dead_s, 'genus', .false.)
        call csv_write(dead_s, 'species', .false.)
        call csv_write(dead_s, 'degday_death', .false.)
        call csv_write(dead_s, 'drought_death', .false.)
        call csv_write(dead_s, 'shade_death', .false.)
        call csv_write(dead_s, 'perm_death', .false.)
        call csv_write(dead_s, 'flood_death', .false.)
        call csv_write(dead_s, 'nutrient_death', .false.)
        call csv_write(dead_s, 'fire_death', .false.)
        call csv_write(dead_s, 'harv_death', .false.)
        call csv_write(dead_s, 'wind_death', .false.)
        call csv_write(dead_s, 'mean_diam', .false.)
        call csv_write(dead_s, 'total_biomC', .false.)
        call csv_write(dead_s, 'total_biomC_sd', .true.)

        ! Genus-level data (dead)
        call csv_write(dead_g, 'siteID', .false.)
        call csv_write(dead_g, 'runID',  .false.)
        call csv_write(dead_g, 'year', .false.)
        call csv_write(dead_g, 'genus', .false.)
        call csv_write(dead_g, 'degday_death', .false.)
        call csv_write(dead_g, 'drought_death', .false.)
        call csv_write(dead_g, 'shade_death', .false.)
        call csv_write(dead_g, 'perm_death', .false.)
        call csv_write(dead_g, 'flood_death', .false.)
        call csv_write(dead_g, 'nutrient_death', .false.)
        call csv_write(dead_g, 'fire_death', .false.)
        call csv_write(dead_g, 'harv_death', .false.)
        call csv_write(dead_g, 'wind_death', .false.)
        call csv_write(dead_g, 'mean_diam', .false.)
        call csv_write(dead_g, 'total_biomC', .false.)
        call csv_write(dead_g, 'total_biomC_sd', .true.)

        ! Across-species data
        call csv_write(plotvals, 'siteID', .false.)
        call csv_write(plotvals, 'runID', .false.)
        call csv_write(plotvals, 'year', .false.)
        call csv_write(plotvals, 'gdd_death', .false.)
        call csv_write(plotvals, 'drought_death', .false.)
        call csv_write(plotvals, 'shade_death', .false.)
        call csv_write(plotvals, 'perm_death', .false.)
        call csv_write(plotvals, 'flood_death', .false.)
        call csv_write(plotvals, 'nutrient_death', .false.)
        call csv_write(plotvals, 'fire_death', .false.)
        call csv_write(plotvals, 'harvest_death', .false.)
        call csv_write(plotvals, 'wind_death', .false.)
        call csv_write(plotvals, 'gddresp_1', .false.)
        call csv_write(plotvals, 'gddresp_2', .false.)
        call csv_write(plotvals, 'gddresp_3', .false.)
        call csv_write(plotvals, 'gddresp_4', .false.)
        call csv_write(plotvals, 'droughtresp_1', .false.)
        call csv_write(plotvals, 'droughtresp_2', .false.)
        call csv_write(plotvals, 'droughtresp_3', .false.)
        call csv_write(plotvals, 'droughtresp_4', .false.)
        call csv_write(plotvals, 'shaderesp_1', .false.)
        call csv_write(plotvals, 'shaderesp_2', .false.)
        call csv_write(plotvals, 'shaderesp_3', .false.)
        call csv_write(plotvals, 'shaderesp_4', .false.)
        call csv_write(plotvals, 'permresp_1', .false.)
        call csv_write(plotvals, 'permresp_2', .false.)
        call csv_write(plotvals, 'permresp_3', .false.)
        call csv_write(plotvals, 'permresp_4', .false.)
        call csv_write(plotvals, 'floodresp_1', .false.)
        call csv_write(plotvals, 'floodresp_2', .false.)
        call csv_write(plotvals, 'floodresp_3', .false.)
        call csv_write(plotvals, 'floodresp_4', .false.)
        call csv_write(plotvals, 'nutrientresp_1', .false.)
        call csv_write(plotvals, 'nutrientresp_2', .false.)
        call csv_write(plotvals, 'nutrientresp_3', .false.)
        call csv_write(plotvals, 'nutrientresp_4', .false.)
        call csv_write(plotvals, 'Loreys_height', .false.)
        call csv_write(plotvals, 'Loreys_height_sd', .false.)
        call csv_write(plotvals, 'max_height', .false.)
        call csv_write(plotvals, 'max_height_sd', .false.)
        call csv_write(plotvals, 'total_biomC', .false.)
        call csv_write(plotvals, 'total_biomC_sd', .false.)
        call csv_write(plotvals, 'basal_area', .false.)
        call csv_write(plotvals, 'basal_area_sd', .false.)
        call csv_write(plotvals, 'total_stems', .false.)
        call csv_write(plotvals, 'total_stems_sd', .false.)
        call csv_write(plotvals, 'small_stems', .false.)
        call csv_write(plotvals, 'small_stems_sd', .false.)
        call csv_write(plotvals, 'med_stems', .false.)
        call csv_write(plotvals, 'med_stems_sd', .false.)
        call csv_write(plotvals, 'lg_stems', .false.)
        call csv_write(plotvals, 'lg_stems_sd', .false.)
        call csv_write(plotvals, 'stand_age', .false.)
        call csv_write(plotvals, 'stand_age_sd', .false.)
        call csv_write(plotvals, 'LAI_1', .false.)
        call csv_write(plotvals, 'LAI_2', .false.)
        call csv_write(plotvals, 'LAI_3', .false.)
        call csv_write(plotvals, 'LAI_4', .false.)
        call csv_write(plotvals, 'LAI_5', .false.)
        call csv_write(plotvals, 'LAI_6', .false.)
        call csv_write(plotvals, 'LAI_7', .false.)
        call csv_write(plotvals, 'LAI_8', .false.)
        call csv_write(plotvals, 'LAI_9', .false.)
        call csv_write(plotvals, 'LAI_10', .false.)
        call csv_write(plotvals, 'LAI_11', .false.)
        call csv_write(plotvals, 'LAI_12', .true.)

        ! Plot-level data
        if (plot_level_data) then

           ! Genus-level data (live)
           call csv_write(pl_biom_by_g, 'siteID', .false.)
           call csv_write(pl_biom_by_g, 'runID', .false.)
           call csv_write(pl_biom_by_g, 'year', .false.)
           call csv_write(pl_biom_by_g, 'plot', .false.)
           call csv_write(pl_biom_by_g, 'genus', .false.)
           call csv_write(pl_biom_by_g, '0 to 5', .false.)
           call csv_write(pl_biom_by_g, '5 to 10', .false.)
           call csv_write(pl_biom_by_g, '10 to 20', .false.)
           call csv_write(pl_biom_by_g, '20 to 30', .false.)
           call csv_write(pl_biom_by_g, '30 to 40', .false.)
           call csv_write(pl_biom_by_g, '40 to 50', .false.)
           call csv_write(pl_biom_by_g, '50 to 60', .false.)
           call csv_write(pl_biom_by_g, '60 to 70', .false.)
           call csv_write(pl_biom_by_g, '70 to 80', .false.)
           call csv_write(pl_biom_by_g, '80 to 90', .false.)
           call csv_write(pl_biom_by_g, 'over 90', .false.)
           call csv_write(pl_biom_by_g, '0 to 5 biom', .false.)
           call csv_write(pl_biom_by_g, '5 to 10 biom', .false.)
           call csv_write(pl_biom_by_g, '10 to 20 biom', .false.)
           call csv_write(pl_biom_by_g, '20 to 30 biom', .false.)
           call csv_write(pl_biom_by_g, '30 to 40 biom', .false.)
           call csv_write(pl_biom_by_g, '40 to 50 biom', .false.)
           call csv_write(pl_biom_by_g, '50 to 60 biom', .false.)
           call csv_write(pl_biom_by_g, '60 to 70 biom', .false.)
           call csv_write(pl_biom_by_g, '70 to 80 biom', .false.)
           call csv_write(pl_biom_by_g, '80 to 90 biom', .false.)
           call csv_write(pl_biom_by_g, 'over 90 biom', .false.)
           call csv_write(pl_biom_by_g, 'max_diam', .false.)
           call csv_write(pl_biom_by_g, 'mean_diam', .false.)
           call csv_write(pl_biom_by_g, 'max_hgt', .false.)
           call csv_write(pl_biom_by_g, 'basal_area', .false.)
           call csv_write(pl_biom_by_g, 'total_biomC', .true.)

           ! Species-level data (live)
           call csv_write(pl_biom_by_s, 'siteID', .false.)
           call csv_write(pl_biom_by_s, 'runID', .false.)
           call csv_write(pl_biom_by_s, 'year', .false.)
           call csv_write(pl_biom_by_s, 'plot', .false.)
           call csv_write(pl_biom_by_s, 'genus', .false.)
           call csv_write(pl_biom_by_s, 'species', .false.)
           call csv_write(pl_biom_by_s, '0 to 5', .false.)
           call csv_write(pl_biom_by_s, '5 to 10', .false.)
           call csv_write(pl_biom_by_s, '10 to 20', .false.)
           call csv_write(pl_biom_by_s, '20 to 30', .false.)
           call csv_write(pl_biom_by_s, '30 to 40', .false.)
           call csv_write(pl_biom_by_s, '40 to 50', .false.)
           call csv_write(pl_biom_by_s, '50 to 60', .false.)
           call csv_write(pl_biom_by_s, '60 to 70', .false.)
           call csv_write(pl_biom_by_s, '70 to 80', .false.)
           call csv_write(pl_biom_by_s, '80 to 90', .false.)
           call csv_write(pl_biom_by_s, 'over 90', .false.)
           call csv_write(pl_biom_by_s, '0 to 5 biom', .false.)
           call csv_write(pl_biom_by_s, '5 to 10 biom', .false.)
           call csv_write(pl_biom_by_s, '10 to 20 biom', .false.)
           call csv_write(pl_biom_by_s, '20 to 30 biom', .false.)
           call csv_write(pl_biom_by_s, '30 to 40 biom', .false.)
           call csv_write(pl_biom_by_s, '40 to 50 biom', .false.)
           call csv_write(pl_biom_by_s, '50 to 60 biom', .false.)
           call csv_write(pl_biom_by_s, '60 to 70 biom', .false.)
           call csv_write(pl_biom_by_s, '70 to 80 biom', .false.)
           call csv_write(pl_biom_by_s, '80 to 90 biom', .false.)
           call csv_write(pl_biom_by_s, 'over 90 biom', .false.)
           call csv_write(pl_biom_by_s, 'max_diam', .false.)
           call csv_write(pl_biom_by_s, 'mean_diam', .false.)
           call csv_write(pl_biom_by_s, 'max_hgt', .false.)
           call csv_write(pl_biom_by_s, 'basal_area', .false.)
           call csv_write(pl_biom_by_s, 'total_biomC', .true.)

           ! Species-level data (dead)
           call csv_write(dead_ps, 'siteID', .false.)
           call csv_write(dead_ps, 'runID', .false.)
           call csv_write(dead_ps, 'year', .false.)
           call csv_write(dead_ps, 'plot', .false.)
           call csv_write(dead_ps, 'genus', .false.)
           call csv_write(dead_ps, 'species', .false.)
           call csv_write(dead_ps, 'degday_death', .false.)
           call csv_write(dead_ps, 'drought_death', .false.)
           call csv_write(dead_ps, 'shade_death', .false.)
           call csv_write(dead_ps, 'perm_death', .false.)
           call csv_write(dead_ps, 'flood_death', .false.)
           call csv_write(dead_ps, 'nutrient_death', .false.)
           call csv_write(dead_ps, 'wind_death', .false.)
           call csv_write(dead_ps, 'harv_death', .false.)
           call csv_write(dead_ps, 'mean_diam', .false.)
           call csv_write(dead_ps, 'total_biomC', .true.)

           ! Genus-level data (dead)
           call csv_write(dead_pg, 'siteID', .false.)
           call csv_write(dead_pg, 'runID', .false.)
           call csv_write(dead_pg, 'year', .false.)
           call csv_write(dead_pg, 'plot', .false.)
           call csv_write(dead_pg, 'genus', .false.)
           call csv_write(dead_pg, 'degday_death', .false.)
           call csv_write(dead_pg, 'drought_death', .false.)
           call csv_write(dead_pg, 'shade_death', .false.)
           call csv_write(dead_pg, 'perm_death', .false.)
           call csv_write(dead_pg, 'flood_death', .false.)
           call csv_write(dead_pg, 'nutrient_death', .false.)
           call csv_write(dead_pg, 'fire_death', .false.)
           call csv_write(dead_pg, 'harv_death', .false.)
           call csv_write(dead_pg, 'wind_death', .false.)
           call csv_write(dead_pg, 'mean_diam', .false.)
           call csv_write(dead_pg, 'total_biomC', .true.)
        endif

        if (tree_level_data) then
            ! Plotlevel individual tree output
            call csv_write(pl_tree, 'siteID', .false.)
            call csv_write(pl_tree, 'runID', .false.)
            call csv_write(pl_tree, 'year', .false.)
            call csv_write(pl_tree, 'plot', .false.)
            call csv_write(pl_tree, 'genus', .false.)
            call csv_write(pl_tree, 'species', .false.)
            call csv_write(pl_tree, 'treeID', .false.)
            call csv_write(pl_tree, 'row', .false.)
            call csv_write(pl_tree, 'col', .false.)
            call csv_write(pl_tree, 'age', .false.)
            call csv_write(pl_tree, 'diam', .false.)
            call csv_write(pl_tree, 'dcbb', .false.)
            call csv_write(pl_tree, 'height', .false.)
            call csv_write(pl_tree, 'cbb_height', .false.)
            call csv_write(pl_tree, 'leaf_biomass', .false.)
            call csv_write(pl_tree, 'leaf_area', .false.)
            call csv_write(pl_tree, 'woody_biomC', .false.)
            call csv_write(pl_tree, 'degd_resp', .false.)
            call csv_write(pl_tree, 'drought_resp', .false.)
            call csv_write(pl_tree, 'shade_resp', .false.)
            call csv_write(pl_tree, 'perm_resp', .false.)
            call csv_write(pl_tree, 'flood_resp', .false.)
            call csv_write(pl_tree, 'nutrient_resp', .true.)
        endif



        ! Testing files
        if (reg_testing) then

            ! Species regeneration
            call csv_write(spec_regen, 'siteID', .false.)
            call csv_write(spec_regen, 'ip', .false.)
            call csv_write(spec_regen, 'year', .false.)
            call csv_write(spec_regen, 'probsum', .false.)
            call csv_write(spec_regen, 'probsum_shrub', .false.)
            call csv_write(spec_regen, 'species', .false.)
            call csv_write(spec_regen, 'seedling', .false.)
            call csv_write(spec_regen, 'seedbank', .false.)
            call csv_write(spec_regen, 'fc_gdd', .false.)
            call csv_write(spec_regen, 'fc_nutr', .false.)
            call csv_write(spec_regen, 'fc_drought', .false.)
            call csv_write(spec_regen, 'fc_flood', .false.)
            call csv_write(spec_regen, 'fc_perm', .false.)
            call csv_write(spec_regen, 'con_light', .false.)
            call csv_write(spec_regen, 'dec_light', .false.)
            call csv_write(spec_regen, 'regrowth', .true.)

        end if

        if (testing) then

            ! Fuel conditions
            call csv_write(fuel_conds, 'siteID', .false.)
            call csv_write(fuel_conds, 'plot', .false.)
            call csv_write(fuel_conds, 'year', .false.)
            call csv_write(fuel_conds, 'day', .false.)
            call csv_write(fuel_conds, 'fuel_sum', .false.)
            call csv_write(fuel_conds, 'fuel_dec', .false.)
            call csv_write(fuel_conds, 'fuel_con', .false.)
            call csv_write(fuel_conds, 'fuel_twig', .false.)
            call csv_write(fuel_conds, 'fuel_smbr', .false.)
            call csv_write(fuel_conds, 'fuel_lgbr', .false.)
            call csv_write(fuel_conds, 'fuel_bole', .false.)
            call csv_write(fuel_conds, 'fuel_moss', .false.)
            call csv_write(fuel_conds, 'fuel_dmoss', .false.)
            call csv_write(fuel_conds, 'fuel_root', .false.)
            call csv_write(fuel_conds, 'fuel_shrub', .false.)
            call csv_write(fuel_conds, 'BD_dec', .false.)
            call csv_write(fuel_conds, 'BD_con', .false.)
            call csv_write(fuel_conds, 'BD_twig', .false.)
            call csv_write(fuel_conds, 'BD_smbr', .false.)
            call csv_write(fuel_conds, 'BD_lgbr', .false.)
            call csv_write(fuel_conds, 'BD_bole', .false.)
            call csv_write(fuel_conds, 'BD_moss', .false.)
            call csv_write(fuel_conds, 'BD_dmoss', .false.)
            call csv_write(fuel_conds, 'BD_root', .false.)
            call csv_write(fuel_conds, 'BD_shrub', .false.)
            call csv_write(fuel_conds, 'SAV_dec', .false.)
            call csv_write(fuel_conds, 'SAV_con', .false.)
            call csv_write(fuel_conds, 'SAV_twig', .false.)
            call csv_write(fuel_conds, 'SAV_smbr', .false.)
            call csv_write(fuel_conds, 'SAV_lgbr', .false.)
            call csv_write(fuel_conds, 'SAV_bole', .false.)
            call csv_write(fuel_conds, 'SAV_moss', .false.)
            call csv_write(fuel_conds, 'SAV_dmoss', .false.)
            call csv_write(fuel_conds, 'SAV_root', .false.)
            call csv_write(fuel_conds, 'SAV_shrub', .false.)
            call csv_write(fuel_conds, 'moist_dec', .false.)
            call csv_write(fuel_conds, 'moist_con', .false.)
            call csv_write(fuel_conds, 'moist_twig', .false.)
            call csv_write(fuel_conds, 'moist_smbr', .false.)
            call csv_write(fuel_conds, 'moist_lgbr', .false.)
            call csv_write(fuel_conds, 'moist_bole', .false.)
            call csv_write(fuel_conds, 'moist_moss', .false.)
            call csv_write(fuel_conds, 'moist_dmoss', .false.)
            call csv_write(fuel_conds, 'moist_root', .false.)
            call csv_write(fuel_conds, 'moist_shrub', .false.)
            call csv_write(fuel_conds, 'mef_dec', .false.)
            call csv_write(fuel_conds, 'mef_con', .false.)
            call csv_write(fuel_conds, 'mef_twig', .false.)
            call csv_write(fuel_conds, 'mef_smbr', .false.)
            call csv_write(fuel_conds, 'mef_lgbr', .false.)
            call csv_write(fuel_conds, 'mef_bole', .false.)
            call csv_write(fuel_conds, 'mef_moss', .false.)
            call csv_write(fuel_conds, 'mef_dmoss', .false.)
            call csv_write(fuel_conds, 'mef_root', .false.)
            call csv_write(fuel_conds, 'mef_shrub', .false.)
            call csv_write(fuel_conds, 'sumlit_moist', .false.)
            call csv_write(fuel_conds, 'sumlit_SAV', .false.)
            call csv_write(fuel_conds, 'sumlit_BD', .false.)
            call csv_write(fuel_conds, 'MEF', .true.)

            call csv_write(dayfire, 'siteID', .false.)
            call csv_write(dayfire, 'plot', .false.)
            call csv_write(dayfire, 'year', .false.)
            call csv_write(dayfire, 'day', .false.)
            call csv_write(dayfire, 'FDI', .false.)
            call csv_write(dayfire, 'ffmc', .false.)
            call csv_write(dayfire, 'MEF', .false.)
            call csv_write(dayfire, 'fuel_moisture', .false.)
            call csv_write(dayfire, 'fuel_BD', .false.)
            call csv_write(dayfire, 'fuel_SAV', .false.)
            call csv_write(dayfire, 'fuel_dec', .false.)
            call csv_write(dayfire, 'fuel_con', .false.)
            call csv_write(dayfire, 'fuel_twig', .false.)
            call csv_write(dayfire, 'fuel_smbr', .false.)
            call csv_write(dayfire, 'fuel_lgbr', .false.)
            call csv_write(dayfire, 'fuel_bole', .false.)
            call csv_write(dayfire, 'fuel_moss', .false.)
            call csv_write(dayfire, 'fuel_dmoss', .false.)
            call csv_write(dayfire, 'fuel_root', .false.)
            call csv_write(dayfire, 'fuel_shrub', .false.)
            call csv_write(dayfire, 'Uf', .false.)
            call csv_write(dayfire, 'A', .false.)
            call csv_write(dayfire, 'beta', .false.)
            call csv_write(dayfire, 'beta_op', .false.)
            call csv_write(dayfire, 'beta_frac', .false.)
            call csv_write(dayfire, 'rVol_max', .false.)
            call csv_write(dayfire, 'rVol_opt', .false.)
            call csv_write(dayfire, 'n_M', .false.)
            call csv_write(dayfire, 'net_fuel', .false.)
            call csv_write(dayfire, 'I_r', .false.)
            call csv_write(dayfire, 'flux_rat', .false.)
            call csv_write(dayfire, 'phi_wind', .false.)
            call csv_write(dayfire, 'Qig', .false.)
            call csv_write(dayfire, 'rosf', .false.)
            call csv_write(dayfire, 'a_f', .false.)
            call csv_write(dayfire, 'burn_dec', .false.)
            call csv_write(dayfire, 'burn_con', .false.)
            call csv_write(dayfire, 'burn_twig', .false.)
            call csv_write(dayfire, 'burn_smbr', .false.)
            call csv_write(dayfire, 'burn_lgbr', .false.)
            call csv_write(dayfire, 'burn_bole', .false.)
            call csv_write(dayfire, 'burn_moss', .false.)
            call csv_write(dayfire, 'burn_dmoss', .false.)
            call csv_write(dayfire, 'burn_root', .false.)
            call csv_write(dayfire, 'burn_shrub', .false.)
            call csv_write(dayfire, 'I_surf', .false.)
            call csv_write(dayfire, 'tau_l', .true.)

            call csv_write(cons_out, 'siteID', .false.)
            call csv_write(cons_out, 'ip', .false.)
            call csv_write(cons_out, 'year', .false.)
            call csv_write(cons_out, 'dmc', .false.)
            call csv_write(cons_out, 'duff_moist', .false.)
            call csv_write(cons_out, 'rfs', .false.)
            call csv_write(cons_out, 'N_cons', .false.)
            call csv_write(cons_out, 'consRoot', .false.)
            call csv_write(cons_out, 'emis', .false.)
            call csv_write(cons_out, 't_r', .false.)
            call csv_write(cons_out, 'duff_cons', .false.)
            call csv_write(cons_out, 'hum_avail', .false.)
            call csv_write(cons_out, 'pre_depth', .false.)
            call csv_write(cons_out, 'O_depth', .false.)
            call csv_write(cons_out, 'hum_combust', .false.)
            call csv_write(cons_out, 'bg_combust', .false.)
            call csv_write(cons_out, 'agw_combust', .false.)
            call csv_write(cons_out, 'agw_prefire', .false.)
            call csv_write(cons_out, 'not_burn', .false.)
            call csv_write(cons_out, 'canopy_bd', .false.)
            call csv_write(cons_out, 'canopy_bh', .false.)
            call csv_write(cons_out, 'canopy_biom', .false.)
            call csv_write(cons_out, 'R_a', .false.)
            call csv_write(cons_out, 'rosf_active', .false.)
            call csv_write(cons_out, 'CFB', .false.)
            call csv_write(cons_out, 'R_final', .false.)
            call csv_write(cons_out, 'I_final', .false.)
            call csv_write(cons_out, 'abcombust', .false.)
            call csv_write(cons_out, 'bgr_combust', .true.)

        end if

    end subroutine write_headers

    !:.........................................................................:

    subroutine close_outputFiles
        !
        !  Closes all opened output files
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code

        ! Data dictionary: local variables
        logical :: isopen ! Is the file open?

        ! Close files if they are open
        inquire(clim_unit, opened = isopen)
        if (isopen) close(clim_unit)

        inquire(biom_by_s, opened = isopen)
        if (isopen) close(biom_by_s)

        inquire(biom_by_g, opened = isopen)
        if (isopen) close(biom_by_g)

        inquire(plotvals, opened = isopen)
        if (isopen) close(plotvals)

        inquire(pl_biom_by_g, opened = isopen)
        if (isopen) close(pl_biom_by_g)

        inquire(pl_biom_by_s, opened = isopen)
        if (isopen) close(pl_biom_by_s)

        inquire(dead_g, opened = isopen)
        if (isopen) close(dead_g)

        inquire(dead_s, opened = isopen)
        if (isopen) close(dead_s)

        inquire(dead_ps, opened = isopen)
        if (isopen) close(dead_ps)

        inquire(dead_pg, opened = isopen)
        if (isopen) close(dead_pg)

        inquire(pl_tree, opened = isopen)
        if (isopen) close(pl_tree)

        inquire(fuel_conds, opened = isopen)
        if (isopen) close(fuel_conds)

        inquire(dayfire, opened = isopen)
        if (isopen) close(dayfire)

        inquire(moist_out, opened = isopen)
        if (isopen) close(moist_out)

        inquire(spec_regen, opened = isopen)
        if (isopen) close(spec_regen)

        inquire(mclimfile, opened = isopen)
        if (isopen) close(mclimfile)

    end subroutine close_outputFiles

    !:.........................................................................:

end module IO
