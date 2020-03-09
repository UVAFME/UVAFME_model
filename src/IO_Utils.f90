module IO

  use Constants
  use Parameters
  use FileUtils
  use csv_file

  implicit none

!*******************************************************************************
  !
  !This module contains utility functions to open and close I/O files
  !
!*******************************************************************************


contains

  !:...........................................................................:

	subroutine read_file_list(filelist)
		!this reads the names of the directories for the input files.
		!if the filelist file is not present it sets it them to default values.
		!Author: Katherine Holcomb, 2012 v. 1.0
		!Inputs:
		!	filelist:  name of filelist

		character(len = *)           :: filelist
		character(len = MAX_FILE)    :: input_directory
		character(len = MAX_FILE)    :: output_directory
		character(len = MAX_DIR)     :: climate_directory
		character(len = MAX_DIR)     :: GCM_directory
		character(len = MAX_DIR)     :: site_directory
		character(len = MAX_DIR)     :: sitelist_directory
		integer                      :: nunit, ios
		logical                      :: file_exists

		!directory list
		namelist /filenames/ input_directory, output_directory,                &
							climate_directory, site_directory,                 &
							sitelist_directory, GCM_directory

		!initialize to default values
		input_directory = 'input_data'
		output_directory = 'output_data'
		climate_directory = 'input_data'
		GCM_directory = 'input_data'
		site_directory = 'input_data'
		sitelist_directory = 'input_data'

		!read file list if it exists
		inquire(file = filelist, exist = file_exists)

		if (.not. file_exists) then
			write(*, *) 'File list file not specified, using defaults'
		else
			nunit = unit_number()
			open(nunit, file = filelist, status = 'unknown', iostat = ios)
			if (ios .eq. 0) then
				read(nunit, filenames, iostat = ios)
				if (ios .ne. 0) then
					write(*, *) 'Error reading ', filelist, ': ', ios,         &
						' Using defaults'
				end if
			endif

			close(nunit)

		endif

		!trim directory names of whitespace
		inputdir = trim(input_directory)
		outputdir = trim(output_directory)
		climatedir = trim(climate_directory)
		GCMdir = trim(GCM_directory)
		sitedir = trim(site_directory)
		slistdir = trim(sitelist_directory)

	end subroutine read_file_list

	!:.........................................................................:

	subroutine read_runFile
		!opens the runtime.txt file if it exists
		!Author: Katherine Holcomb, 2012 v. 1.0

		character(len = MAX_FILE)  :: filename
		character(len = MAX_PATH)  :: pathname

		!find the current working directory
		call get_cwd(fqcwd)

		!set file prefix
		write(code_id, '(a,a)') codename, version_id

		!open runtime file
		call build_filename(code_id, '_runtime.txt', filename)
		call build_pathname(inputdir, filename, pathname)
		rt_file = open_file(pathname, 'r')

		!if can't find runtime file, set to defaults
		if (rt_file .eq. invalid) then
			write(*, *) 'Unable to open runtime (parameter) file.'
			write(*, *) 'Using all default values.'
		endif

	end subroutine read_runFile

	!:.........................................................................:

	subroutine open_inputFiles
		!opens all other input files
		!if mandatory files do not exist, program is halted
		!Author: Katherine Holcomb, 2012 v. 1.0

		character(len = MAX_FILE)  :: filename
		character(len = MAX_PATH)  :: pathname

		!get current working directory
		call get_cwd(fqcwd)

		!logfile (must be opened as an "input" file so other input can write
			!to it)
		call build_pathname(outputdir, 'log.txt', pathname)
		logf = open_file(pathname)

		!site log file (must be opened as "input" file so others can write to it)
		call build_pathname(outputdir, 'site_log.txt', pathname)
		sitef = open_file(pathname)

		!site list
		call build_filename(code_id, '_sitelist.csv', filename)
		call build_pathname(slistdir, filename,pathname)
		slist = open_file(pathname, 'r')
		if (slist .eq. invalid) then
			call fatal_error('Unable to find site list file.')
		endif

		!species attributes file
		call build_filename(code_id, '_specieslist.csv', filename)
		call build_pathname(inputdir, filename, pathname)
		splist = open_file(pathname, 'r')
		if (splist .eq. invalid) then
			call fatal_error('Unable to find species-list file.')
		endif

		!site file
		call build_filename(code_id, '_site.csv', filename)
		call build_pathname(sitedir, filename, pathname)
		sfile = open_file(pathname, 'r')
		if (sfile .eq. invalid) then
			call fatal_error('Unable to find site data file.')
		endif
		!the site file has a header line
		if (count_records(sfile, 1) <=0 ) then
			call fatal_error("Error in site data file")
		endif

		!climate file
		call build_filename(code_id, '_climate.csv', filename)
		call build_pathname(climatedir, filename, pathname)
		cfile = open_file(pathname, 'r')
		if (cfile .eq. invalid) then
			call fatal_error('Unable to find climate file.')
		endif
		if ( count_records(cfile, 1) <=0 ) then
			call fatal_error("Error in climate data file")
		endif

		!extra climate file (cld)
		call build_filename(code_id, '_climate_ex.csv', filename)
		call build_pathname(climatedir, filename, pathname)
		cexfile = open_file(pathname, 'r')
		if (cexfile .eq. invalid) then
			call fatal_error('Unable to find climte file')
		end if
		if (count_records(cexfile, 1) <= 0) then
			call fatal_error("Error in extra climate data file")
		end if

		!litter parameters
		call build_filename(code_id, '_litterpars.csv', filename)
		call build_pathname(sitedir, filename, pathname)
		litterfile = open_file(pathname, 'r')
		if (litterfile .eq. invalid) then
			call fatal_error('Unable to find litter parameter file')
		end if
		if (count_records(litterfile, 1) <= 0) then
			call fatal_error("Error in litter parameter file")
		end if

		!climate standard-deviation file
		call build_filename(code_id, '_climate_stddev.csv', filename)
		call build_pathname(climatedir, filename, pathname)
		cstdfile = open_file(pathname, 'r')
		if (cstdfile .eq. invalid) then
			write(*, *)                                                		   &
				'No climate stdev file specified or invalid file name.'
			use_climstd = .false.
		else
			if (count_records(cstdfile, 1) <=0 ) then
				call warning("Error in climatology std file, ignoring")
				use_climstd = .false.
			else
				use_climstd = .true.
			endif
		endif

		!extra climate data standard deviation file
		call build_filename(code_id, '_climate_ex_stddev.csv', filename)
		call build_pathname(climatedir, filename, pathname)
		cexstdfile = open_file(pathname, 'r')
		if (cexstdfile .eq. invalid) then
			write(*, *)                                                        &
				'No extra climate stdev file specified or invalid file name.'
			use_cexstd = .false.
		else
			if (count_records(cexstdfile, 1) <= 0) then
				call warning('Error in extra climate stddev file, ignoring')
				use_cexstd = .false.
			else
				use_cexstd = .true.
			end if
		end if

		!climate data from GCM
		if ( use_gcm ) then
			call build_filename(code_id, '_climate_GCM.csv', filename)
			call build_pathname(GCMdir, filename, pathname)
			cgcmfile = open_file(pathname, 'r')
		endif

		!rangelist file
		call build_filename(code_id, '_rangelist.csv', filename)
		call build_pathname(inputdir, filename, pathname)
		rglist = open_file(pathname, 'r')
		if (rglist .eq. invalid) then
			write(*, *) 'No rangelist file specified or invalid file name.'
			write(*, *) 'Using all species in all sites.'
			use_rangelist = .false.
		else
			use_rangelist = .true.
		endif

	end subroutine open_inputFiles

	!:.........................................................................:

	subroutine open_outputFiles
		!opens all output files
		!Author: Katherine Holcomb, 2012 v. 1.0

		character(len = MAX_PATH) :: pathname

		!open files

		!soil decomposition
		call build_pathname(outputdir, 'SoilDecomp.csv', pathname)
		soildecomp = open_file(pathname)

		!climate
		call build_pathname(outputdir, 'Climate.csv', pathname)
		clim_unit = open_file(pathname)

		!genus-level data (live)
		call build_pathname(outputdir, 'Genus_Data.csv', pathname)
		biom_by_g = open_file(pathname)

		!species-level data (live)
		call build_pathname(outputdir, 'Species_Data.csv', pathname)
		biom_by_s = open_file(pathname)

		!species-level data (dead)
		call build_pathname(outputdir, 'Dead_Species_Data.csv', pathname)
		dead_s = open_file(pathname)

		!genus-level data (dead)
		call build_pathname(outputdir, 'Dead_Genus_Data.csv', pathname)
		dead_g = open_file(pathname)

		!across-species data
		call build_pathname(outputdir, 'Total_Plot_Values.csv', pathname)
		plotvals = open_file(pathname)

		!plotlevel data
		if (plot_level_data) then
			call build_pathname(outputdir, 'Plot_Genus_Data.csv', pathname)
			pl_biom_by_g = open_file(pathname)

			call build_pathname(outputdir, 'Plot_Dead_Genus.csv', pathname)
			dead_pg = open_file(pathname)

			call build_pathname(outputdir, 'Plot_Species_Data.csv', pathname)
			pl_biom_by_s = open_file(pathname)

			call build_pathname(outputdir, 'Plot_Dead_Species.csv', pathname)
			dead_ps = open_file(pathname)
		end if

		!treelevel data
		if (tree_level_data) then
			call build_pathname(outputdir, 'Plot_Tree_Data.csv', pathname)
			tld = open_file(pathname)
		end if


	end subroutine open_outputFiles

	!:.........................................................................:

	subroutine write_headers()
		!writes headers for all output files
		!Author: Katherine Holcomb, 2012 v. 1.0

		!write file headers

		!soil
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
		call csv_write(soildecomp, 'lit_WDW', .false.)
		call csv_write(soildecomp, 'lit_moss', .false.)
		call csv_write(soildecomp, 'avail_n', .true.)

		!climate
		call csv_write(clim_unit, 'siteID', .false.)
		call csv_write(clim_unit, 'runID',  .false.)
		call csv_write(clim_unit, 'year', .false.)
		call csv_write(clim_unit, 'rain', .false.)
		call csv_write(clim_unit, 'pet', .false.)
		call csv_write(clim_unit, 'solar_rad', .false.)
		call csv_write(clim_unit, 'tdd', .false.)
		call csv_write(clim_unit, 'thaw_depth', .false.)
		call csv_write(clim_unit, 'organic_depth', .false.)
		call csv_write(clim_unit, 'avail_n', .false.)
		call csv_write(clim_unit, 'aet', .false.)
		call csv_write(clim_unit, 'grow', .false.)
        call csv_write(clim_unit, 'pc_germ', .false.)
		call csv_write(clim_unit, 'degd', .false.)
		call csv_write(clim_unit, 'dryd_upper', .false.)
		call csv_write(clim_unit, 'dryd_base', .false.)
        call csv_write(clim_unit, 'wilt_days', .false.)
		call csv_write(clim_unit, 'flood_d', .false.)
		call csv_write(clim_unit, 'waterlog_d', .true.)

		!genus-level data (live)
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
		call csv_write(biom_by_g, 'total_biomC_sd', .true.)


		!species-level data (live)
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
		call csv_write(biom_by_s, 'total_biomC_sd', .true.)


		!species-level data (dead)
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
		call csv_write(dead_s, 'wind_death', .false.)
		call csv_write(dead_s, 'mean_diam', .false.)
		call csv_write(dead_s, 'total_biomC', .false.)
		call csv_write(dead_s, 'total_biomC_sd', .true.)

		!genus-level data (dead)
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
		call csv_write(dead_g, 'wind_death', .false.)
		call csv_write(dead_g, 'mean_diam', .false.)
		call csv_write(dead_g, 'total_biomC', .false.)
		call csv_write(dead_g, 'total_biomC_sd', .true.)

		!across-species data
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
		call csv_write(plotvals, 'wind_death', .false.)
		call csv_write(plotvals, 'gddresp_1', .false.)
		call csv_write(plotvals, 'droughtresp_1', .false.)
		call csv_write(plotvals, 'shaderesp_1', .false.)
		call csv_write(plotvals, 'permresp_1', .false.)
		call csv_write(plotvals, 'floodresp_1', .false.)
		call csv_write(plotvals, 'nutrientresp_1', .false.)
		call csv_write(plotvals, 'gddresp_2', .false.)
		call csv_write(plotvals, 'droughtresp_2', .false.)
		call csv_write(plotvals, 'shaderesp_2', .false.)
		call csv_write(plotvals, 'permresp_2', .false.)
		call csv_write(plotvals, 'floodresp_2', .false.)
		call csv_write(plotvals, 'nutrientresp_2', .false.)
		call csv_write(plotvals, 'gddresp_3', .false.)
		call csv_write(plotvals, 'droughtresp_3', .false.)
		call csv_write(plotvals, 'shaderesp_3', .false.)
		call csv_write(plotvals, 'permresp_3', .false.)
		call csv_write(plotvals, 'floodresp_3', .false.)
		call csv_write(plotvals, 'nutrientresp_3', .false.)
		call csv_write(plotvals, 'gddlim_1', .false.)
		call csv_write(plotvals, 'droughtlim_1', .false.)
		call csv_write(plotvals, 'shadelim_1', .false.)
		call csv_write(plotvals, 'permlim_1', .false.)
		call csv_write(plotvals, 'floodlim_1', .false.)
		call csv_write(plotvals, 'nutrientlim_1', .false.)
		call csv_write(plotvals, 'gddlim_2', .false.)
		call csv_write(plotvals, 'droughtlim_2', .false.)
		call csv_write(plotvals, 'shadelim_2', .false.)
		call csv_write(plotvals, 'permlim_2', .false.)
		call csv_write(plotvals, 'floodlim_2', .false.)
		call csv_write(plotvals, 'nutrientlim_2', .false.)
		call csv_write(plotvals, 'gddlim_3', .false.)
		call csv_write(plotvals, 'droughtlim_3', .false.)
		call csv_write(plotvals, 'shadelim_3', .false.)
		call csv_write(plotvals, 'permlim_3', .false.)
		call csv_write(plotvals, 'floodlim_3', .false.)
		call csv_write(plotvals, 'nutrientlim_3', .false.)
		call csv_write(plotvals, 'Loreys_height', .false.)
		call csv_write(plotvals, 'Loreys_height_sd', .false.)
		call csv_write(plotvals, 'Max_height', .false.)
		call csv_write(plotvals, 'Max_height_sd', .false.)
		call csv_write(plotvals, 'Total_plot_biomc', .false.)
		call csv_write(plotvals, 'Total_plot_biomc_sd', .false.)
		call csv_write(plotvals, 'Total_basal_area', .false.)
		call csv_write(plotvals, 'Total_basal_area_sd', .false.)
		call csv_write(plotvals, 'Total_stems', .false.)
		call csv_write(plotvals, 'Total_stems_sd', .false.)
		call csv_write(plotvals, 'Stand_age', .false.)
		call csv_write(plotvals, 'Stand_age_sd', .false.)
		call csv_write(plotvals, 'Weighted_age', .false.)
		call csv_write(plotvals, 'Weighted_age_sd', .false.)
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

		if (plot_level_data) then

		   !genus-level data (live)
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
		   call csv_write(pl_biom_by_g, 'leaf_area_ind', .false.)
		   call csv_write(pl_biom_by_g, 'basal_area', .false.)
		   call csv_write(pl_biom_by_g, 'total_biomC', .true.)

		   !species-level data (live)
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
		   call csv_write(pl_biom_by_s, 'leaf_area_ind', .false.)
		   call csv_write(pl_biom_by_s, 'basal_area', .false.)
		   call csv_write(pl_biom_by_s, 'total_biomC', .true.)

		   !species-level data (dead)
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
		   call csv_write(dead_ps, 'fire_death', .false.)
		   call csv_write(dead_ps, 'wind_death', .false.)
		   call csv_write(dead_ps, 'mean_diam', .false.)
		   call csv_write(dead_ps, 'total_biomC', .true.)


		   !genus-level data (dead)
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
		   call csv_write(dead_pg, 'wind_death', .false.)
		   call csv_write(dead_pg, 'mean_diam', .false.)
		   call csv_write(dead_pg, 'total_biomC', .true.)

		endif

		if (tree_level_data) then
			!plotlevel individual tree output
			call csv_write(tld, 'siteID', .false.)
			call csv_write(tld, 'runID', .false.)
			call csv_write(tld, 'year', .false.)
			call csv_write(tld, 'plot', .false.)
			call csv_write(tld, 'treenum', .false.)
			call csv_write(tld, 'genus', .false.)
			call csv_write(tld, 'species', .false.)
			call csv_write(tld, 'row', .false.)
			call csv_write(tld, 'col', .false.)
            call csv_write(tld, 'age', .false.)
			call csv_write(tld, 'diam', .false.)
			call csv_write(tld, 'height', .false.)
			call csv_write(tld, 'cbb_height', .false.)
			call csv_write(tld, 'leaf_biomass', .false.)
			call csv_write(tld, 'stem_biomC', .false.)
            call csv_write(tld, 'degd_resp', .false.)
            call csv_write(tld, 'drought_resp', .false.)
            call csv_write(tld, 'shade_resp', .false.)
            call csv_write(tld, 'perm_resp', .false.)
            call csv_write(tld, 'flood_resp', .false.)
            call csv_write(tld, 'nutrient_resp', .true.)
		endif

	end subroutine write_headers

	!:.........................................................................:

	subroutine close_outputFiles
		!closes all opened output files
		!Author: Katherine Holcomb, 2012 v. 1.0
		logical                  :: isopen

		close(c_and_n)
		close(clim_unit)
		close(biom_by_s)
		close(biom_by_g)
		close(plotvals)

		!optional files
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

		inquire(tld, opened = isopen)
		if (isopen) close(pl_tree)

	end subroutine close_outputFiles

	!:.........................................................................:

end module IO
