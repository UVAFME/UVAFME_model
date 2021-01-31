module Output

  use Constants
  use Parameters
  use Tree
  use Plot
  use GenusGroups
  use IO
  use csv_file

  implicit none

!*******************************************************************************
  !
  !This module initializes and writes data to output files.
  !
!*******************************************************************************

	!plot-scale unit adjustments
	real  :: plotscale, plotadj, plotrenorm

contains

	!:.........................................................................:

	subroutine initialize_outputFiles()
		!opens output files and writes their headers

		call open_outputFiles
		call write_headers()

	end subroutine initialize_outputFiles

	!:.........................................................................:

	subroutine total_plot_values(site, maxcells, year)
		!writes across-species data to Total_Plot_Values.csv file
		!Authors: Adrianna Foster and Jacquelyn Shuman 2014, v. 1.0
		!Inputs:
		!	site:     site instance
		!	maxcells: maximum cells per plot (maxcells*maxcells = max trees)
		!	year:     simulation year

		type(SiteData),   intent(in)       :: site
		integer,          intent(in)       :: maxcells
		integer,          intent(in)       :: year

		real, dimension(:, :), allocatable :: tree_ht, tree_biom
		real, dimension(:, :), allocatable :: tree_basal_area, tree_age
		real, dimension(site%numplots)     :: LH, basal_ar, WA
		real, dimension(site%numplots)     :: plbiomc, stems, plstandage
		real, dimension(site%numplots)     :: plmaxheight, sm_stems, med_stems
		real, dimension(site%numplots)     :: n_1, n_2, n_3, lg_stems
		real, dimension(site%numplots, 6)  :: envresp_1, envresp_2
		real, dimension(site%numplots, 6)  :: envresp_3
		real, dimension(site%numplots, 6)  :: scounts_1, scounts_2
		real, dimension(site%numplots, 6)  :: scounts_3
		real, dimension(12)                :: site_lai
		real, dimension(8)                 :: mort_markers
		real, dimension(6)                 :: envresp_1mn, envresp_2mn
		real, dimension(6)                 :: envresp_3mn, envresp_1sd
		real, dimension(6)                 :: envresp_2sd, envresp_3sd
		real, dimension(6)                 :: scounts_1mn, scounts_2mn
		real, dimension(6)                 :: scounts_3mn, scounts_1sd
		real, dimension(6)                 :: scounts_2sd, scounts_3sd
		real                               :: LH_num, LH_denom
		real                               :: WA_num, WA_denom
		real                               :: LH_mean, LH_std
		real                               :: Totalbiomc, Totalbiomc_std
		real                               :: Totalbasal, Totalbasal_std
		real                               :: Totalstems, Totalstems_std
		real                               :: Maxheight, Maxheight_std
		real                               :: Avg_age, Avg_age_std
		real                               :: WA_mean, WA_std
        real                               :: lg_stems_mn, lg_stems_std
        real                               :: sm_stems_mn, sm_stems_std
        real                               :: med_stems_mn, med_stems_std
	    integer                            :: ip, it, m, stress

		!for converting to per ha
		plotscale = hec_to_m2/plotsize

		!for converting to per ha and averaging
		plotadj = plotscale/float(site%numplots)

		!for converting to per m2 and averaging
		plotrenorm = 1.0/plotsize/float(site%numplots)

		allocate(tree_ht(site%numplots, maxcells*maxcells))
		allocate(tree_basal_area(site%numplots, maxcells*maxcells))
		allocate(tree_biom(site%numplots, maxcells*maxcells))
		allocate(tree_age(site%numplots, maxcells*maxcells))

		plmaxheight = 0.0
		basal_ar = 0.0
		plbiomc = 0.0
		tree_ht = 0.0
		tree_basal_area = 0.0
		tree_biom = 0.0
		tree_age = 0.0
		site_lai = 0.0
		mort_markers = 0.0
		envresp_1 = 0.0
		envresp_2 = 0.0
		envresp_3 = 0.0
		scounts_1 = 0.0
		scounts_2 = 0.0
		scounts_3 = 0.0
		n_1 = 0.0
		n_2 = 0.0
		n_3 = 0.0
        lg_stems = 0.0
        sm_stems = 0.0
        med_stems = 0.0

		!bin LAI into 5-m sections
        site_lai(1) = sum(site%lai_array(1:5))
        site_lai(2) = sum(site%lai_array(6:10))
        site_lai(3) = sum(site%lai_array(11:15))
        site_lai(4) = sum(site%lai_array(16:20))
        site_lai(5) = sum(site%lai_array(21:25))
        site_lai(6) = sum(site%lai_array(26:30))
        site_lai(7) = sum(site%lai_array(31:35))
        site_lai(8) = sum(site%lai_array(36:40))
        site_lai(9) = sum(site%lai_array(41:45))
        site_lai(10) = sum(site%lai_array(46:50))
        site_lai(11) = sum(site%lai_array(51:55))
        site_lai(12) = sum(site%lai_array(56:60))


		do ip = 1, site%numplots

			do it = 1, site%plots(ip)%numtrees

                if (site%plots(ip)%trees(it)%diam_bht .le. 5.0) then
                    sm_stems(ip) = sm_stems(ip) + 1.0
                else if (site%plots(ip)%trees(it)%diam_bht .gt. 5.0 .and.      &
                    site%plots(ip)%trees(it)%diam_bht .le. 20.0) then
                    med_stems(ip) = med_stems(ip) + 1.0
                else if (site%plots(ip)%trees(it)%diam_bht .gt. 20.0) then
                    lg_stems(ip) = lg_stems(ip) + 1.0
                end if

				!grab individual tree height (m)
				tree_ht(ip, it) = site%plots(ip)%trees(it)%forska_ht

				!calculate individual tree basal area (m2)
				tree_basal_area(ip, it) =                                      &
					0.000025*pi*site%plots(ip)%trees(it)%diam_bht**2

				!sum up plot-level biomass (tC)
				plbiomc(ip) = plbiomc(ip) + (site%plots(ip)%trees(it)%biomC +  &
					site%plots(ip)%trees(it)%leaf_bm)

				!sum up plot-level basal area (cm2)
				basal_ar(ip) = basal_ar(ip) +                                  &
					0.25*pi*site%plots(ip)%trees(it)%diam_bht**2

				!grab individual tree biomass (tC)
				tree_biom(ip, it) = site%plots(ip)%trees(it)%biomC +           &
					site%plots(ip)%trees(it)%leaf_bm

				!grab individual tree age (years)
				tree_age(ip, it) = site%plots(ip)%trees(it)%tree_age

				!get maximum height on plot (m)
				plmaxheight(ip) = max(plmaxheight(ip),                         &
					site%plots(ip)%trees(it)%forska_ht)

				!sum up growth response to each stressor by height class
				if (site%plots(ip)%trees(it)%forska_ht .lt. 10.0) then

					envresp_1(ip, 1) = envresp_1(ip, 1) +                      &
						site%plots(ip)%trees(it)%env_resp(1)
					envresp_1(ip, 2) = envresp_1(ip, 2) +                      &
						site%plots(ip)%trees(it)%env_resp(2)
					envresp_1(ip, 3) = envresp_1(ip, 3) +                      &
							site%plots(ip)%trees(it)%env_resp(3)
					envresp_1(ip, 4) = envresp_1(ip, 4) +                      &
							site%plots(ip)%trees(it)%env_resp(4)
					envresp_1(ip, 5) = envresp_1(ip, 5) +                      &
							site%plots(ip)%trees(it)%env_resp(5)
					envresp_1(ip, 6) = envresp_1(ip, 6) +                      &
							site%plots(ip)%trees(it)%env_resp(6)

					!get number in this class
					n_1(ip) = n_1(ip) + 1.0

				else if ((site%plots(ip)%trees(it)%forska_ht .ge. 10.0)        &
					.and. (site%plots(ip)%trees(it)%forska_ht .lt. 20.0)) then

					envresp_2(ip, 1) = envresp_2(ip, 1) +                      &
						site%plots(ip)%trees(it)%env_resp(1)
					envresp_2(ip, 2) = envresp_2(ip, 2) +                      &
						site%plots(ip)%trees(it)%env_resp(2)
					envresp_2(ip, 3) = envresp_2(ip, 3) +                      &
						site%plots(ip)%trees(it)%env_resp(3)
					envresp_2(ip, 4) = envresp_2(ip, 4) +                      &
						site%plots(ip)%trees(it)%env_resp(4)
					envresp_2(ip, 5) = envresp_2(ip, 5) +                      &
						site%plots(ip)%trees(it)%env_resp(5)
					envresp_2(ip, 6) = envresp_2(ip, 6) +                      &
						site%plots(ip)%trees(it)%env_resp(6)

					!get number in this class
					n_2(ip) = n_2(ip) + 1.0

				else if (site%plots(ip)%trees(it)%forska_ht .ge. 20.0) then

					envresp_3(ip, 1) = envresp_3(ip, 1) +                      &
						site%plots(ip)%trees(it)%env_resp(1)
					envresp_3(ip, 2) = envresp_3(ip, 2) +                      &
						site%plots(ip)%trees(it)%env_resp(2)
					envresp_3(ip, 3) = envresp_3(ip, 3) +                      &
						site%plots(ip)%trees(it)%env_resp(3)
					envresp_3(ip, 4) = envresp_3(ip, 4) +                      &
						site%plots(ip)%trees(it)%env_resp(4)
					envresp_3(ip, 5) = envresp_3(ip, 5) +                      &
						site%plots(ip)%trees(it)%env_resp(5)
					envresp_3(ip, 6) = envresp_3(ip, 6) +                      &
						site%plots(ip)%trees(it)%env_resp(6)

					!get number in this class
					n_3(ip) = n_3(ip) + 1.0

				end if

				!get number of trees for each most limiting factor by height
				  !class
				if (site%plots(ip)%trees(it)%forska_ht .lt. 10.0) then

					stress = site%plots(ip)%trees(it)%stressor
					scounts_1(ip, stress) = scounts_1(ip, stress) + 1.0

				else if ((site%plots(ip)%trees(it)%forska_ht .ge. 10.0) &
					.and. (site%plots(ip)%trees(it)%forska_ht .lt. 20.0)) then

					stress = site%plots(ip)%trees(it)%stressor
					scounts_2(ip, stress) = scounts_2(ip, stress) + 1.0

				else if (site%plots(ip)%trees(it)%forska_ht .ge. 20.0) then

					stress = site%plots(ip)%trees(it)%stressor
					scounts_3(ip, stress) = scounts_3(ip, stress) + 1.0

				end if
			enddo

		    !Lorey's Height calculation
			LH_num = sum((tree_basal_area(ip, :))*(tree_ht(ip, :)))
			LH_denom = sum(tree_basal_area(ip, :))
			if (LH_denom .eq. 0.0) then
				LH(ip) = 0.0
			else
				LH(ip) = LH_num/LH_denom
			endif

			!weighted stand age calculation
			WA_num = sum((tree_biom(ip, :))*(tree_age(ip, :)))
			WA_denom = sum(tree_biom(ip, :))
			if (WA_denom .eq. 0.0) then
				WA(ip) = 0.0
			else
				WA(ip) = WA_num/WA_denom
			endif

			!plot level stems
			stems(ip) = site%plots(ip)%numtrees

			!plot level stand age (not weighted) - age of oldest tree
			plstandage(ip) = site%plots(ip)%stand_age

			!mortality markers - biomass killed per stressor
			mort_markers(1) = mort_markers(1) + site%plots(ip)%d_type(1)
			mort_markers(2) = mort_markers(2) + site%plots(ip)%d_type(2)
			mort_markers(3) = mort_markers(3) + site%plots(ip)%d_type(3)
			mort_markers(4) = mort_markers(4) + site%plots(ip)%d_type(4)
			mort_markers(5) = mort_markers(5) + site%plots(ip)%d_type(5)
			mort_markers(6) = mort_markers(6) + site%plots(ip)%d_type(6)
			mort_markers(7) = mort_markers(7) + site%plots(ip)%d_type(7)
			mort_markers(8) = mort_markers(8) + site%plots(ip)%d_type(8)


			!get average or proportion of each stressor type
			if (n_1(ip) .ge. 1.0) then
				envresp_1(ip, :) = envresp_1(ip, :)/n_1(ip)
				scounts_1(ip, :) = scounts_1(ip, :)/n_1(ip)
			else
				envresp_1(ip, :) = rnvalid
				scounts_1(ip, :) = rnvalid
			end if

			if (n_2(ip) .ge. 1.0) then
				envresp_2(ip, :) = envresp_2(ip, :)/n_2(ip)
				scounts_2(ip, :) = scounts_2(ip, :)/n_2(ip)
			else
				envresp_2(ip, :) = rnvalid
				scounts_2(ip, :) = rnvalid
			end if

			if (n_3(ip) .ge. 1.0) then
				envresp_3(ip, :) = envresp_3(ip, :)/n_3(ip)
				scounts_3(ip, :) = scounts_3(ip, :)/n_3(ip)
			else
				envresp_3(ip, :) = rnvalid
				scounts_3(ip, :) = rnvalid
			end if

		enddo

		!get average of mortality bins and convert to tC/ha
		mort_markers(:) = mort_markers(:)/float(site%numplots)/plotsize*       &
			hec_to_m2

		!compute means and stds
		call stddev(LH, LH_mean, LH_std, rnvalid)
		call stddev(plbiomc, Totalbiomc, Totalbiomc_std, rnvalid)
		call stddev(basal_ar, Totalbasal, Totalbasal_std, rnvalid)
		call stddev(stems, Totalstems, Totalstems_std, rnvalid)
		call stddev(plmaxheight, Maxheight, Maxheight_std, rnvalid)
		call stddev(plstandage, Avg_age, Avg_age_std, rnvalid)
		call stddev(WA, WA_mean, WA_std, rnvalid)
        call stddev(lg_stems, lg_stems_mn, lg_stems_std, rnvalid)
        call stddev(med_stems, med_stems_mn, med_stems_std, rnvalid)
        call stddev(sm_stems, sm_stems_mn, sm_stems_std, rnvalid)

		do m = 1, 6
			call stddev(envresp_1(:, m), envresp_1mn(m), envresp_1sd(m),       &
				rnvalid)
			call stddev(envresp_2(:, m), envresp_2mn(m), envresp_2sd(m),       &
				rnvalid)
			call stddev(envresp_3(:, m), envresp_3mn(m), envresp_3sd(m),       &
				rnvalid)
		end do

		do m = 1, 6
			call stddev(scounts_1(:, m), scounts_1mn(m), scounts_1sd(m),       &
				rnvalid)
			call stddev(scounts_2(:, m), scounts_2mn(m), scounts_2sd(m),       &
				rnvalid)
			call stddev(scounts_3(:, m), scounts_3mn(m), scounts_3sd(m),       &
				rnvalid)
		end do

		Totalbiomc = plotscale*Totalbiomc  !converted to tc/ha
		Totalbiomc_std = plotscale*Totalbiomc_std !tc/ha
		Totalbasal = Totalbasal/plotsize*hec_to_m2*0.0001 !m2/ha
		Totalbasal_std = Totalbasal_std/plotsize*hec_to_m2 !cm2/ha
		Totalstems = Totalstems*plotscale !converted to total stems/ha
		Totalstems_std = Totalstems*plotscale !converted to stems/ha
        lg_stems_mn = lg_stems_mn*plotscale !converted to total stems/ha
        lg_stems_std = lg_stems_std*plotscale !converted to stems/ha
        med_stems_mn = med_stems_mn*plotscale !converted to total stems/ha
        med_stems_std = med_stems_std*plotscale !converted to stems/ha
        sm_stems_mn = sm_stems_mn*plotscale !converted to total stems/ha
        sm_stems_std = sm_stems_std*plotscale !converted to stems/ha

		!write to csv
		call csv_write(plotvals, site%site_id, .false.)
		call csv_write(plotvals, site%runID, .false.)
		call csv_write(plotvals, year, .false.)
		call csv_write(plotvals, mort_markers(1), .false.)
		call csv_write(plotvals, mort_markers(2), .false.)
		call csv_write(plotvals, mort_markers(3), .false.)
		call csv_write(plotvals, mort_markers(4), .false.)
		call csv_write(plotvals, mort_markers(5), .false.)
		call csv_write(plotvals, mort_markers(6), .false.)
		call csv_write(plotvals, mort_markers(7), .false.)
		call csv_write(plotvals, mort_markers(8), .false.)
		call csv_write(plotvals, envresp_1mn(1), .false.)
		call csv_write(plotvals, envresp_1mn(2), .false.)
		call csv_write(plotvals, envresp_1mn(3), .false.)
		call csv_write(plotvals, envresp_1mn(4), .false.)
		call csv_write(plotvals, envresp_1mn(5), .false.)
		call csv_write(plotvals, envresp_1mn(6), .false.)
		call csv_write(plotvals, envresp_2mn(1), .false.)
		call csv_write(plotvals, envresp_2mn(2), .false.)
		call csv_write(plotvals, envresp_2mn(3), .false.)
		call csv_write(plotvals, envresp_2mn(4), .false.)
		call csv_write(plotvals, envresp_2mn(5), .false.)
		call csv_write(plotvals, envresp_2mn(6), .false.)
		call csv_write(plotvals, envresp_3mn(1), .false.)
		call csv_write(plotvals, envresp_3mn(2), .false.)
		call csv_write(plotvals, envresp_3mn(3), .false.)
		call csv_write(plotvals, envresp_3mn(4), .false.)
		call csv_write(plotvals, envresp_3mn(5), .false.)
		call csv_write(plotvals, envresp_3mn(6), .false.)
		call csv_write(plotvals, scounts_1mn(1), .false.)
		call csv_write(plotvals, scounts_1mn(2), .false.)
		call csv_write(plotvals, scounts_1mn(3), .false.)
		call csv_write(plotvals, scounts_1mn(4), .false.)
		call csv_write(plotvals, scounts_1mn(5), .false.)
		call csv_write(plotvals, scounts_1mn(6), .false.)
		call csv_write(plotvals, scounts_2mn(1), .false.)
		call csv_write(plotvals, scounts_2mn(2), .false.)
		call csv_write(plotvals, scounts_2mn(3), .false.)
		call csv_write(plotvals, scounts_2mn(4), .false.)
		call csv_write(plotvals, scounts_2mn(5), .false.)
		call csv_write(plotvals, scounts_2mn(6), .false.)
		call csv_write(plotvals, scounts_3mn(1), .false.)
		call csv_write(plotvals, scounts_3mn(2), .false.)
		call csv_write(plotvals, scounts_3mn(3), .false.)
		call csv_write(plotvals, scounts_3mn(4), .false.)
		call csv_write(plotvals, scounts_3mn(5), .false.)
		call csv_write(plotvals, scounts_3mn(6), .false.)
		call csv_write(plotvals, LH_mean, .false.)
		call csv_write(plotvals, LH_std, .false.)
		call csv_write(plotvals, Maxheight, .false.)
		call csv_write(plotvals, Maxheight_std, .false.)
		call csv_write(plotvals, Totalbiomc, .false.)
		call csv_write(plotvals, Totalbiomc_std, .false.)
		call csv_write(plotvals, Totalbasal, .false.)
		call csv_write(plotvals, Totalbasal_std, .false.)
		call csv_write(plotvals, Totalstems, .false.)
		call csv_write(plotvals, Totalstems_std, .false.)
        call csv_write(plotvals, sm_stems_mn, .false.)
        call csv_write(plotvals, sm_stems_std, .false.)
        call csv_write(plotvals, med_stems_mn, .false.)
        call csv_write(plotvals, med_stems_std, .false.)
        call csv_write(plotvals, lg_stems_mn, .false.)
        call csv_write(plotvals, lg_stems_std, .false.)
		call csv_write(plotvals, Avg_age, .false.)
		call csv_write(plotvals, Avg_age_std, .false.)
		call csv_write(plotvals, WA_mean, .false.)
		call csv_write(plotvals, WA_std, .false.)
		call csv_write(plotvals, site_lai(1), .false.)
		call csv_write(plotvals, site_lai(2), .false.)
		call csv_write(plotvals, site_lai(3), .false.)
		call csv_write(plotvals, site_lai(4), .false.)
		call csv_write(plotvals, site_lai(5), .false.)
		call csv_write(plotvals, site_lai(6), .false.)
		call csv_write(plotvals, site_lai(7), .false.)
		call csv_write(plotvals, site_lai(8), .false.)
		call csv_write(plotvals, site_lai(9), .false.)
		call csv_write(plotvals, site_lai(10), .false.)
		call csv_write(plotvals, site_lai(11), .false.)
		call csv_write(plotvals, site_lai(12), .true.)

	end subroutine total_plot_values

  !:...........................................................................:

	subroutine write_genus_data(site, species_pres, year)
		!writes output data to Genus_Data.csv and Plot_Genus_Data.csv
		  !file if writing plot-level data
		!Author: Katherine Holcomb 2012, v. 1.0
		!Inputs:
		!	site:         site instance
		!	species_pres: species present at site
		!	year:         simulation year

		type(SiteData),   intent(in)            :: site
		type(Groups),     intent(in)            :: species_pres
		integer,          intent(in)            :: year

		integer, dimension(site%numplots, species_pres%numgenera, NHC) :: diam_categories
		real,    dimension(site%numplots, species_pres%numgenera, NHC) :: biom_categories
		real,    dimension(site%numplots, species_pres%numgenera, 6) :: envi_resp
		real,    dimension(site%numplots, species_pres%numgenera) :: basal_area
		real,    dimension(site%numplots, species_pres%numgenera) :: biomC,    &
																	 biomN,    &
												                     leaf_bm,  &
												                     leaf_area, &
												                     max_ht,   &
												                     max_diam, &
												                     plotlevel_biomc, &
												                     plotlevel_biomn, &
												                     plotlevel_basal, &
												                     mean_dbh, &
												                     mean_year

		real, dimension(species_pres%numgenera, NHC) :: total_dm_cats
		real,    dimension(species_pres%numgenera, NHC) :: total_biom_cats
		real,    dimension(species_pres%numgenera, 6) :: total_envi_resp
		real,    dimension(species_pres%numgenera) :: total_biomc,             &
													  total_biomn,             &
													  total_basal,             &
                                                      total_laibm,             &
                                                      total_max_ht,            &
                                                      total_max_diam,          &
                                                      plot_mean_biomc,         &
                                                      plot_std_biomc,          &
                                                      plot_mean_biomn,         &
                                                      plot_std_biomn,          &
                                                      plot_mean_basal,         &
                                                      plot_std_basal,          &
                                                      plot_mean_dbh,           &
                                                      plot_std_dbh,            &
                                                      plot_mean_year,          &
                                                      plot_std_year,           &
                                                      plot_mean_Maxdbh,        &
                                                      plot_std_Maxdbh,         &
                                                      plot_mean_Maxht,         &
                                                      plot_std_Maxht
        real,   dimension(species_pres%numspecies, 6)    :: plot_env_resp_mn,  &
												            plot_env_resp_sd

		integer, dimension(species_pres%numspecies) :: spec_index
		logical, dimension(species_pres%numgenera)  :: in_site

		character(len = 8)                     :: field = 'genus'
		integer                                :: iwmo
		integer                                :: ip, is, ns, l

		plotscale = hec_to_m2/plotsize
		plotadj = plotscale/float(site%numplots)
		plotrenorm = 1.0/plotsize/float(site%numplots)

		iwmo = site%site_id

	    !check to see which genera are in the site
		in_site = .false.
		do is = 1, species_pres%numgenera
			do ip = 1, site%numplots
				do ns = 1, size(site%plots(ip)%species)
					if (species_pres%genusgroups(is) .eq.                      &
						site%plots(ip)%species(ns)%genus_name) then
						spec_index(is) = ns
						in_site(is) = .true.
						exit
					endif
				enddo
			enddo
		enddo

		!get individual plot data for each genus
		do ip = 1,site%numplots
			call sum_over_sg(site%plots(ip), species_pres%genusgroups, field,  &
							basal_area(ip, :), leaf_bm(ip, :), biomC(ip, :),   &
							biomN(ip, :), max_ht(ip, :), max_diam(ip, :),      &
							mean_dbh(ip, :), mean_year(ip, :),                 &
							biom_categories(ip, :, :), envi_resp(ip, :, :))

			call tree_dm_cats(site%plots(ip), species_pres%genusgroups,        &
							field, diam_categories(ip, :, :))
		enddo

		total_biomc = 0.0
		total_biomn = 0.0
		total_basal = 0.0
		total_laibm = 0.0
		total_max_ht = 0.0
		total_max_diam = 0.0
		plotlevel_biomc = rnvalid
		plotlevel_biomn = rnvalid
		plotlevel_basal = rnvalid
		total_dm_cats = 0
		total_biom_cats = 0.0
		total_envi_resp = 0.0

		!accumulate across-plot data
		do is = 1, species_pres%numgenera
			do ip = 1, site%numplots
				if (in_site(is)) then
					ns = spec_index(is)
					total_biomc(is) = total_biomc(is) + biomc(ip, is)
					total_biomn(is) = total_biomn(is) + biomn(ip, is)
					total_basal(is) = total_basal(is) + basal_area(ip, is)
					total_laibm(is) = total_laibm(is) + leaf_bm(ip, is)/       &
						(site%plots(ip)%species(ns)%leafarea_c*2.0)
					total_max_ht(is) = max(total_max_ht(is), max_ht(ip, is))
					total_max_diam(is) = max(total_max_diam(is),               &
						max_diam(ip, is))
					total_dm_cats(is, :) = total_dm_cats(is, :) +              &
						float(diam_categories(ip, is, :))
					total_biom_cats(is, :) = total_biom_cats(is, :) +          &
						biom_categories(ip, is, :)
					total_envi_resp(is, :) = total_envi_resp(is, :) +          &
						envi_resp(ip, is, :)
				endif
			enddo

			!get mean and standard deviation of data
			call stddev(biomc(:, is), plot_mean_biomc(is), plot_std_biomc(is), &
				rnvalid)
			call stddev(biomn(:, is), plot_mean_biomn(is), plot_std_biomn(is), &
				rnvalid)
			call stddev(mean_dbh(:, is), plot_mean_dbh(is), plot_std_dbh(is),  &
				rnvalid)
			call stddev(mean_year(:, is), plot_mean_year(is),                  &
				plot_std_year(is), rnvalid)
			call stddev(basal_area(:, is), plot_mean_basal(is),                &
				plot_std_basal(is), rnvalid)

			call stddev(max_diam(:, is), plot_mean_Maxdbh(is),                 &
				plot_std_Maxdbh(is), rnvalid)
			call stddev(max_ht(:, is), plot_mean_Maxht(is),                    &
				plot_std_Maxht(is), rnvalid)

			do l = 1, 6
				call stddev(envi_resp(:, is, l), plot_env_resp_mn(is, l),      &
					plot_env_resp_sd(is, l), rnvalid)
			end do
		enddo

		!convert values to per ha
		plotlevel_biomc = plotscale*biomc
		plotlevel_basal = plotscale*basal_area*0.0001 !m2/ha


		plot_mean_biomc = plotscale*plot_mean_biomc
		plot_std_biomc = plotscale*plot_std_biomc
		plotlevel_biomn = plotscale*biomn
		plot_mean_biomn = plotscale*plot_mean_biomn
		plot_std_biomn = plotscale*plot_std_biomn
		plot_mean_basal = plot_mean_basal/plotsize*hec_to_m2*0.0001
		plot_std_basal = plot_std_basal/plotsize*hec_to_m2*0.0001

		!write to csv
		do is = 1, species_pres%numgenera
			call csv_write(biom_by_g, iwmo, .false.)
			call csv_write(biom_by_g, site%runID, .false.)
			call csv_write(biom_by_g, year, .false.)
			call csv_write(biom_by_g,                                          &
				trim(adjustl(species_pres%genusgroups(is))), .false.)

			if (in_site(is)) then

				total_basal(is) = total_basal(is)*plotadj*0.0001
				total_laibm(is) = total_laibm(is)*plotrenorm
				total_biomc(is) = total_biomc(is)*plotadj
				total_biomn(is) = total_biomn(is)*plotadj*10.0
				total_dm_cats(is, :) = total_dm_cats(is, :)*plotadj
				total_biom_cats(is, :) = total_biom_cats(is, :)*plotadj
				total_envi_resp(is, :) = total_envi_resp(is, :)/numplots

			    call csv_write(biom_by_g, total_dm_cats(is, 1), .false.)
			    call csv_write(biom_by_g, total_dm_cats(is, 2), .false.)
			    call csv_write(biom_by_g, total_dm_cats(is, 3), .false.)
				call csv_write(biom_by_g, total_dm_cats(is, 4), .false.)
				call csv_write(biom_by_g, total_dm_cats(is, 5), .false.)
				call csv_write(biom_by_g, total_dm_cats(is, 6), .false.)
				call csv_write(biom_by_g, total_dm_cats(is, 7), .false.)
				call csv_write(biom_by_g, total_dm_cats(is, 8), .false.)
				call csv_write(biom_by_g, total_dm_cats(is, 9), .false.)
				call csv_write(biom_by_g, total_dm_cats(is, 10), .false.)
				call csv_write(biom_by_g, total_dm_cats(is, 11), .false.)

				call csv_write(biom_by_g, total_biom_cats(is, 1), .false.)
				call csv_write(biom_by_g, total_biom_cats(is, 2), .false.)
				call csv_write(biom_by_g, total_biom_cats(is, 3), .false.)
				call csv_write(biom_by_g, total_biom_cats(is, 4), .false.)
				call csv_write(biom_by_g, total_biom_cats(is, 5), .false.)
				call csv_write(biom_by_g, total_biom_cats(is, 6), .false.)
				call csv_write(biom_by_g, total_biom_cats(is, 7), .false.)
				call csv_write(biom_by_g, total_biom_cats(is, 8), .false.)
				call csv_write(biom_by_g, total_biom_cats(is, 9), .false.)
				call csv_write(biom_by_g, total_biom_cats(is, 10), .false.)
				call csv_write(biom_by_g, total_biom_cats(is, 11), .false.)

				call csv_write(biom_by_g, plot_env_resp_mn(is, 1), .false.)
				call csv_write(biom_by_g, plot_env_resp_mn(is, 2), .false.)
				call csv_write(biom_by_g, plot_env_resp_mn(is, 3), .false.)
				call csv_write(biom_by_g, plot_env_resp_mn(is, 4), .false.)
				call csv_write(biom_by_g, plot_env_resp_mn(is, 5), .false.)
				call csv_write(biom_by_g, plot_env_resp_mn(is, 6), .false.)

				call csv_write(biom_by_g, plot_mean_Maxdbh(is), .false.)
				call csv_write(biom_by_g, plot_mean_dbh(is), .false.)
				call csv_write(biom_by_g, plot_mean_year(is), .false.)
				call csv_write(biom_by_g, plot_mean_Maxht(is), .false.)
				call csv_write(biom_by_g, total_laibm(is), .false.)
				call csv_write(biom_by_g, total_basal(is), .false.)
				call csv_write(biom_by_g, plot_std_basal(is), .false.)
				call csv_write(biom_by_g, total_biomc(is),.false.)
				call csv_write(biom_by_g, plot_std_biomc(is), .true.)

			else

				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .false.)
				call csv_write(biom_by_g, rnvalid, .true.)

			endif

		enddo

		!write plot-level data if turned on
		if (plot_level_data) then

			do ip = 1, site%numplots
				do is = 1, species_pres%numgenera
					call csv_write(pl_biom_by_g, iwmo, .false.)
					call csv_write(pl_biom_by_g, site%runID, .false.)
					call csv_write(pl_biom_by_g, year, .false.)
					call csv_write(pl_biom_by_g, ip, .false.)
					call csv_write(pl_biom_by_g,                               &
						trim(adjustl(species_pres%genusgroups(is))), .false.)

					if (in_site(is)) then

						diam_categories(ip, is, :) =                           &
							diam_categories(ip, is, :)*plotscale

						biom_categories(ip, is, :) =                           &
							biom_categories(ip, is, :)*plotscale

						ns = spec_index(is)
						leaf_area(ip, is) = leaf_bm(ip, is)/                   &
							(site%plots(ip)%species(ns)%leafarea_c*2.0)/plotsize


						call csv_write(pl_biom_by_g,                           &
							diam_categories(ip, is, 1), .false.)
						call csv_write(pl_biom_by_g,                           &
							diam_categories(ip, is, 2), .false.)
						call csv_write(pl_biom_by_g,                           &
							diam_categories(ip, is, 3), .false.)
						call csv_write(pl_biom_by_g,                           &
							diam_categories(ip, is, 4), .false.)
						call csv_write(pl_biom_by_g,                           &
							diam_categories(ip, is, 5), .false.)
						call csv_write(pl_biom_by_g,                           &
							diam_categories(ip, is, 6), .false.)
						call csv_write(pl_biom_by_g,                           &
							diam_categories(ip, is, 7), .false.)
						call csv_write(pl_biom_by_g,                           &
							diam_categories(ip, is, 8), .false.)
						call csv_write(pl_biom_by_g,                           &
							diam_categories(ip, is, 9), .false.)
						call csv_write(pl_biom_by_g,                           &
							diam_categories(ip, is, 10), .false.)
						call csv_write(pl_biom_by_g,                           &
							diam_categories(ip, is, 11), .false.)

						call csv_write(pl_biom_by_g,                           &
							biom_categories(ip, is, 1), .false.)
						call csv_write(pl_biom_by_g,                           &
							biom_categories(ip, is, 2), .false.)
						call csv_write(pl_biom_by_g,                           &
							biom_categories(ip, is, 3), .false.)
						call csv_write(pl_biom_by_g,                           &
							biom_categories(ip, is, 4), .false.)
						call csv_write(pl_biom_by_g,                           &
							biom_categories(ip, is, 5), .false.)
						call csv_write(pl_biom_by_g,                           &
							biom_categories(ip, is, 6), .false.)
						call csv_write(pl_biom_by_g,                           &
							biom_categories(ip, is, 7), .false.)
						call csv_write(pl_biom_by_g,                           &
							biom_categories(ip, is, 8), .false.)
						call csv_write(pl_biom_by_g,                           &
							biom_categories(ip, is, 9), .false.)
						call csv_write(pl_biom_by_g,                           &
							biom_categories(ip, is, 10), .false.)
						call csv_write(pl_biom_by_g,                           &
							biom_categories(ip, is, 11), .false.)

						call csv_write(pl_biom_by_g, max_diam(ip, is), .false.)
						call csv_write(pl_biom_by_g, mean_dbh(ip, is), .false.)
						call csv_write(pl_biom_by_g, max_ht(ip, is), .false.)
						call csv_write(pl_biom_by_g, leaf_area(ip, is), .false.)
						call csv_write(pl_biom_by_g, plotlevel_basal(ip, is),  &
							.false.)
						call csv_write(pl_biom_by_g, plotlevel_biomc(ip, is),  &
							.true.)

					else
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .false.)
						call csv_write(pl_biom_by_g, rnvalid, .true.)
					endif
				enddo
			enddo
		endif

	end subroutine write_genus_data

	!:.........................................................................:

	subroutine write_species_data(site, species_pres, year)
		!writes species-level data to Species_Data.csv file and
		  !Plot_Species_Data.csv if plot-level output turned on
		!Author: Katherine Holcomb 2012, v. 1.0
		!Inputs:
		!	site:         site instance
		!	species_pres: species present at site
		!	year:         simulation year

		type(SiteData),   intent(in)            :: site
		type(Groups),     intent(in)            :: species_pres
		integer,          intent(in)            :: year

		integer, dimension(site%numplots, species_pres%numspecies, NHC) :: diam_categories
		real,    dimension(site%numplots, species_pres%numspecies, NHC) :: biom_categories
		real,    dimension(site%numplots, species_pres%numspecies, 6)   :: envi_resp
		real,    dimension(site%numplots, species_pres%numspecies) :: biomC,   &
																	  basal_area, &
                                                                      biomN,   &
                                                                      leaf_bm, &
                                                                      leaf_area, &
                                                                      max_ht,  &
                                                                      max_diam, &
                                                                      plotlevel_biomc, &
                                                                      plotlevel_biomn, &
                                                                      plotlevel_basal, &
                                                                      mean_dbh, &
                                                                      mean_year
		real,    dimension(species_pres%numspecies, NHC) :: total_dm_cats
		real,    dimension(species_pres%numspecies, NHC) :: total_biom_cats
		real,    dimension(species_pres%numspecies, 6)   :: total_envi_resp
		real,    dimension(species_pres%numspecies)      :: total_biomc,       &
                                                            total_biomn,       &
                                                            total_basal,       &
                                                            total_laibm,       &
                                                            total_max_ht,      &
                                                            total_max_diam,    &
                                                            plot_mean_biomc,   &
                                                            plot_std_biomc,    &
                                                            plot_mean_biomn,   &
                                                            plot_std_biomn,    &
                                                            plot_mean_basal,   &
                                                            plot_std_basal,    &
                                                            plot_mean_dbh,     &
                                                            plot_std_dbh,      &
                                                            plot_mean_year,    &
                                                            plot_std_year,     &
                                                            plot_mean_Maxdbh,  &
                                                            plot_std_Maxdbh,   &
                                                            plot_mean_Maxht,   &
                                                            plot_std_Maxht
        real,   dimension(species_pres%numspecies, 6)    :: plot_env_resp_mn,  &
												            plot_env_resp_sd

		integer, dimension(species_pres%numspecies)      :: spec_index
		logical, dimension(species_pres%numspecies)      :: in_site

		character(len = 8)                               :: field = 'species'
		integer                                          :: iwmo
		integer                                          :: ip, is, ns, l

		plotscale = hec_to_m2/plotsize
		plotadj = plotscale/float(site%numplots)
		plotrenorm = 1.0/plotsize/float(site%numplots)

		iwmo = site%site_id

		!check which species are in the site
		in_site = .false.
		do is = 1, species_pres%numspecies
			do ip = 1, site%numplots
				do ns = 1, size(site%plots(ip)%species)
					if (species_pres%spec_names(is, 2) .eq.                    &
							site%plots(ip)%species(ns)%unique_id) then
						in_site(is) = .true.
						spec_index(is) = ns
						exit
					endif
				enddo
			enddo
		enddo

		!get individual plot-level data
		do ip = 1, site%numplots
			call sum_over_sg(site%plots(ip), species_pres%spec_names(:, 2),    &
							field, basal_area(ip, :), leaf_bm(ip, :),          &
							biomC(ip, :), biomN(ip, :), max_ht(ip, :),         &
							max_diam(ip, :), mean_dbh(ip, :),                  &
							mean_year(ip, :), biom_categories(ip, :, :),       &
							envi_resp(ip, :, :))

			call tree_dm_cats(site%plots(ip), species_pres%spec_names(:, 2),   &
							field, diam_categories(ip, :, :))
		enddo

		total_biomc = 0.0
		total_biomn = 0.0
		total_basal = 0.0
		total_laibm = 0.0
		total_max_ht = 0.0
		total_max_diam = 0.0

		plotlevel_biomc = rnvalid
		plotlevel_biomn = rnvalid
		plotlevel_basal = rnvalid

		total_dm_cats = 0
		total_biom_cats = 0.0
		total_envi_resp = 0.0

		!accumulate data across plots
		do is = 1, species_pres%numspecies
			if (in_site(is)) then

				ns = spec_index(is)

				do ip = 1, site%numplots

					total_biomc(is) = total_biomc(is) + biomc(ip, is)
					total_biomn(is) = total_biomn(is) + biomn(ip, is)

					total_basal(is) = total_basal(is) +  basal_area(ip, is)
					total_laibm(is) = total_laibm(is) +  leaf_bm(ip, is)/      &
                         (site%plots(ip)%species(ns)%leafarea_c*2.0)
					total_max_ht(is) = max(total_max_ht(is), max_ht(ip, is))
					total_max_diam(is) = max(total_max_diam(is),               &
						max_diam(ip, is))
					total_dm_cats(is, :) = total_dm_cats(is, :) +              &
						float(diam_categories(ip, is, :))
					total_biom_cats(is, :) = total_biom_cats(is, :) +          &
						biom_categories(ip, is, :)
					total_envi_resp(is, :) = total_envi_resp(is, :) +          &
						envi_resp(ip, is, :)
				enddo
			endif

			!calculate standard deviations of data
			call stddev(basal_area(:, is), plot_mean_basal(is),                &
				plot_std_basal(is), rnvalid)
			call stddev(biomc(:, is), plot_mean_biomc(is),                     &
                plot_std_biomc(is), rnvalid)
			call stddev(biomn(:, is), plot_mean_biomn(is),                     &
                plot_std_biomn(is), rnvalid)
			call stddev(mean_dbh(:, is), plot_mean_dbh(is),                    &
				plot_std_dbh(is), rnvalid)
			call stddev(mean_year(:, is), plot_mean_year(is),                  &
				plot_std_year(is), rnvalid)
			call stddev(max_diam(:, is), plot_mean_Maxdbh(is),                 &
				plot_std_Maxdbh(is), rnvalid)
			call stddev(max_ht(:, is), plot_mean_Maxht(is),                    &
				plot_std_Maxht(is), rnvalid)

			do l = 1, 6
				call stddev(envi_resp(:, is, l), plot_env_resp_mn(is, l),      &
					plot_env_resp_sd(is, l), rnvalid)
			end do

		enddo

		!convert values to per ha
		plotlevel_biomc = plotscale*biomc !tC/ha
		plotlevel_basal = plotscale*basal_area*0.0001 !m2/ha
		plotlevel_biomn = plotscale*biomn !tN/ha

		plot_mean_biomc = plotscale*plot_mean_biomc !tC/ha
		plot_std_biomc = plotscale*plot_std_biomc !tC/ha


		plot_mean_biomn = plotscale*plot_mean_biomn !tN/ha
		plot_std_biomn = plotscale*plot_std_biomn !tN/ha

		plot_mean_basal = plot_mean_basal*plotscale*0.0001 !m2/ha
		plot_std_basal = plot_std_basal*plotscale*0.0001 !m2/ha

		!write to csv
		do is = 1, species_pres%numspecies

			call csv_write(biom_by_s, iwmo, .false.)
			call csv_write(biom_by_s, site%runID, .false.)
			call csv_write(biom_by_s, year, .false.)
			call csv_write(biom_by_s,                                          &
				trim(adjustl(species_pres%spec_names(is, 1))), .false.)
			call csv_write(biom_by_s,                                          &
				trim(adjustl(species_pres%spec_names(is, 2))), .false.)

			if (in_site(is)) then

				total_basal(is) = total_basal(is)*plotadj*0.0001
				total_laibm(is) = total_laibm(is)*plotrenorm
				total_biomc(is) = total_biomc(is)*plotadj
				total_biomn(is) = total_biomn(is)*plotadj*10.0

				total_dm_cats(is, :) = total_dm_cats(is, :)*plotadj
				total_biom_cats(is, :) = total_biom_cats(is, :)*plotadj
				total_envi_resp(is, :) = total_envi_resp(is, :)/numplots

				call csv_write(biom_by_s, total_dm_cats(is, 1), .false.)
				call csv_write(biom_by_s, total_dm_cats(is, 2), .false.)
				call csv_write(biom_by_s, total_dm_cats(is, 3), .false.)
				call csv_write(biom_by_s, total_dm_cats(is, 4), .false.)
				call csv_write(biom_by_s, total_dm_cats(is, 5), .false.)
				call csv_write(biom_by_s, total_dm_cats(is, 6), .false.)
				call csv_write(biom_by_s, total_dm_cats(is, 7), .false.)
				call csv_write(biom_by_s, total_dm_cats(is, 8), .false.)
				call csv_write(biom_by_s, total_dm_cats(is, 9), .false.)
				call csv_write(biom_by_s, total_dm_cats(is, 10), .false.)
				call csv_write(biom_by_s, total_dm_cats(is, 11), .false.)

				call csv_write(biom_by_s, total_biom_cats(is, 1), .false.)
				call csv_write(biom_by_s, total_biom_cats(is, 2), .false.)
				call csv_write(biom_by_s, total_biom_cats(is, 3), .false.)
				call csv_write(biom_by_s, total_biom_cats(is, 4), .false.)
				call csv_write(biom_by_s, total_biom_cats(is, 5), .false.)
				call csv_write(biom_by_s, total_biom_cats(is, 6), .false.)
				call csv_write(biom_by_s, total_biom_cats(is, 7), .false.)
				call csv_write(biom_by_s, total_biom_cats(is, 8), .false.)
				call csv_write(biom_by_s, total_biom_cats(is, 9), .false.)
				call csv_write(biom_by_s, total_biom_cats(is, 10), .false.)
				call csv_write(biom_by_s, total_biom_cats(is, 11), .false.)

				call csv_write(biom_by_s, plot_env_resp_mn(is, 1), .false.)
				call csv_write(biom_by_s, plot_env_resp_mn(is, 2), .false.)
				call csv_write(biom_by_s, plot_env_resp_mn(is, 3), .false.)
				call csv_write(biom_by_s, plot_env_resp_mn(is, 4), .false.)
				call csv_write(biom_by_s, plot_env_resp_mn(is, 5), .false.)
				call csv_write(biom_by_s, plot_env_resp_mn(is, 6), .false.)

				call csv_write(biom_by_s, plot_mean_Maxdbh(is), .false.)
				call csv_write(biom_by_s, plot_mean_dbh(is), .false.)
				call csv_write(biom_by_s, plot_mean_year(is), .false.)
				call csv_write(biom_by_s, plot_mean_Maxht(is), .false.)
				call csv_write(biom_by_s, total_laibm(is), .false.)
				call csv_write(biom_by_s, total_basal(is), .false.)
				call csv_write(biom_by_s, plot_std_basal(is), .false.)
				call csv_write(biom_by_s, total_biomc(is), .false.)
				call csv_write(biom_by_s, plot_std_biomc(is), .true.)

			else

			    call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .false.)
				call csv_write(biom_by_s, rnvalid, .true.)

			endif

		enddo

		!write plot-level data if turned on
		if (plot_level_data) then

			do ip = 1, site%numplots
				do is = 1, species_pres%numspecies
					call csv_write(pl_biom_by_s, iwmo, .false.)
					call csv_write(pl_biom_by_s, site%runID, .false.)
					call csv_write(pl_biom_by_s, year, .false.)
					call csv_write(pl_biom_by_s, ip, .false.)
					call csv_write(pl_biom_by_s,                               &
						trim(adjustl(species_pres%spec_names(is, 1))), .false.)
					call csv_write(pl_biom_by_s,                               &
						trim(adjustl(species_pres%spec_names(is, 2))), .false.)

					if (in_site(is)) then

						biom_categories(ip, is, :) =                           &
							biom_categories(ip, is, :)*plotscale

						diam_categories(ip, is, :) =                           &
							diam_categories(ip, is, :)*plotscale

						ns = spec_index(is)
						leaf_area(ip, is) = leaf_bm(ip, is)/                   &
							(site%plots(ip)%species(ns)%leafarea_c*2.0)/plotsize

						call csv_write(pl_biom_by_s,                           &
							diam_categories(ip, is, 1), .false.)
						call csv_write(pl_biom_by_s,                           &
							diam_categories(ip, is, 2), .false.)
						call csv_write(pl_biom_by_s,                           &
							diam_categories(ip, is, 3), .false.)
						call csv_write(pl_biom_by_s,                           &
							diam_categories(ip, is, 4), .false.)
						call csv_write(pl_biom_by_s,                           &
							diam_categories(ip, is, 5), .false.)
						call csv_write(pl_biom_by_s,                           &
							diam_categories(ip, is, 6), .false.)
						call csv_write(pl_biom_by_s,                           &
							diam_categories(ip, is, 7), .false.)
						call csv_write(pl_biom_by_s,                           &
							diam_categories(ip, is, 8), .false.)
						call csv_write(pl_biom_by_s,                           &
							diam_categories(ip, is, 9), .false.)
						call csv_write(pl_biom_by_s,                           &
							diam_categories(ip, is, 10), .false.)
						call csv_write(pl_biom_by_s,                           &
							diam_categories(ip, is, 11), .false.)

						call csv_write(pl_biom_by_s,                           &
							biom_categories(ip, is, 1), .false.)
						call csv_write(pl_biom_by_s,                           &
							biom_categories(ip, is, 2), .false.)
						call csv_write(pl_biom_by_s,                           &
							biom_categories(ip, is, 3), .false.)
						call csv_write(pl_biom_by_s,                           &
							biom_categories(ip, is, 4), .false.)
						call csv_write(pl_biom_by_s,                           &
							biom_categories(ip, is, 5), .false.)
						call csv_write(pl_biom_by_s,                           &
							biom_categories(ip, is, 6), .false.)
						call csv_write(pl_biom_by_s,                           &
							biom_categories(ip, is, 7), .false.)
						call csv_write(pl_biom_by_s,                           &
							biom_categories(ip, is, 8), .false.)
						call csv_write(pl_biom_by_s,                           &
							biom_categories(ip, is, 9), .false.)
						call csv_write(pl_biom_by_s,                           &
							biom_categories(ip, is, 10), .false.)
						call csv_write(pl_biom_by_s,                           &
							biom_categories(ip, is, 11), .false.)

						call csv_write(pl_biom_by_s, max_diam(ip, is), .false.)
						call csv_write(pl_biom_by_s, mean_dbh(ip, is), .false.)
						call csv_write(pl_biom_by_s, max_ht(ip, is), .false.)
						call csv_write(pl_biom_by_s, leaf_area(ip, is), .false.)
						call csv_write(pl_biom_by_s, plotlevel_basal(ip, is),  &
							.false.)
						call csv_write(pl_biom_by_s, plotlevel_biomc(ip, is),  &
							.true.)
					else
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .false.)
						call csv_write(pl_biom_by_s, rnvalid, .true.)
					endif

				enddo
			enddo

		endif

	end subroutine write_species_data

	!:.........................................................................:

	subroutine write_dead_genus_data(site, species_pres, year)
		!writes dead genus-level data to Dead_Genus_Data.csv and to
		  !Plotlevel_Dead_Genus.csv if turned on
		!Authors: Adrianna Foster and Jacquelyn Shuman 2015, v. 1.0
		!Inputs:
		!	site:         site instance
		!	species_pres: species present at site
		!	year:         simulation year

		type(SiteData),   intent(in)            :: site
		type(Groups),     intent(in)            :: species_pres
		integer,          intent(in)            :: year

		integer, dimension(site%numplots, species_pres%numgenera, NHC) :: diam_categories
		real,    dimension(site%numplots, species_pres%numgenera, 8)   :: d_markers
		real,    dimension(site%numplots, species_pres%numgenera) ::   biomC,  &
                                                                       plotlevel_biomc, &
                                                                       mean_dbh
		integer, dimension(species_pres%numgenera, NHC) :: total_dm_cats
		real,    dimension(species_pres%numgenera, 8)   :: total_mort_markers

		real,    dimension(species_pres%numgenera)      :: total_biomc,        &
													       plot_mean_biomc,    &
                                                           plot_std_biomc,     &
                                                           plot_mean_dbh,      &
                                                           plot_std_dbh

		integer, dimension(species_pres%numspecies)     :: spec_index
		logical, dimension(species_pres%numgenera)      :: in_site
		character(len = 8)                              :: field='genus'
		integer                                         :: iwmo
		integer                                         :: ip, is, ns

		plotscale = hec_to_m2/plotsize
		plotadj = plotscale/float(site%numplots)
		plotrenorm = 1.0/plotsize/float(site%numplots)

		iwmo = site%site_id

		!check if in site
		in_site = .false.
		do is = 1, species_pres%numgenera
			do ip = 1, site%numplots
				do ns = 1, size(site%plots(ip)%species)
					if (species_pres%genusgroups(is) .eq.                      &
							site%plots(ip)%species(ns)%genus_name) then
						spec_index(is) = ns
						in_site(is) = .true.
						exit
					endif
				enddo
			enddo
		enddo

		!get individual plot data
		do ip = 1, site%numplots
			call sum_over_deadsg(site%plots(ip), species_pres%genusgroups,     &
				field, biomC(ip, :), mean_dbh(ip, :), d_markers(ip, :, :))

			call tree_dm_cats(site%plots(ip), species_pres%genusgroups,        &
				field, diam_categories(ip, :, :))
		enddo

		total_biomc = 0.0
		plotlevel_biomc = rnvalid
		total_dm_cats = 0
		total_mort_markers = 0

		do is = 1, species_pres%numgenera
			do ip = 1, site%numplots
				if (in_site(is)) then
					ns = spec_index(is)
					total_biomc(is) = total_biomc(is) + biomc(ip,is)
					total_dm_cats(is, :) = total_dm_cats(is, :) +              &
						diam_categories(ip, is, :)
					total_mort_markers(is, :) = total_mort_markers(is, :) +    &
						d_markers(ip, is, :)
				endif
			enddo

			!calculate standard deviation of data
			call stddev(biomc(:, is), plot_mean_biomc(is), plot_std_biomc(is), &
				rnvalid)
			call stddev(mean_dbh(:, is), plot_mean_dbh(is), plot_std_dbh(is),  &
				rnvalid)
		enddo

		!convert values to per ha
		plotlevel_biomc = plotscale*biomc
		plot_mean_biomc = plotscale*plot_mean_biomc
		plot_std_biomc = plotscale*plot_std_biomc

		!write to csv
		do is = 1, species_pres%numgenera

			call csv_write(dead_g, iwmo, .false.)
			call csv_write(dead_g, site%runID, .false.)
			call csv_write(dead_g, year, .false.)
			call csv_write(dead_g,                                             &
				trim(adjustl(species_pres%genusgroups(is))), .false.)

			if (in_site(is)) then

				total_biomc(is) = total_biomc(is)*plotadj
				total_mort_markers(is, :) = total_mort_markers(is, :)*plotadj

				call csv_write(dead_g, total_mort_markers(is, 1), .false.)
				call csv_write(dead_g, total_mort_markers(is, 2), .false.)
				call csv_write(dead_g, total_mort_markers(is, 3), .false.)
				call csv_write(dead_g, total_mort_markers(is, 4), .false.)
				call csv_write(dead_g, total_mort_markers(is, 5), .false.)
				call csv_write(dead_g, total_mort_markers(is, 6), .false.)
				call csv_write(dead_g, total_mort_markers(is, 7), .false.)
				call csv_write(dead_g, total_mort_markers(is, 8), .false.)
				call csv_write(dead_g, plot_mean_dbh(is), .false.)
				call csv_write(dead_g, total_biomc(is), .false.)
				call csv_write(dead_g, plot_std_biomc(is), .true.)

			else

				call csv_write(dead_g, rnvalid, .false.)
				call csv_write(dead_g, rnvalid, .false.)
				call csv_write(dead_g, rnvalid, .false.)
				call csv_write(dead_g, rnvalid, .false.)
				call csv_write(dead_g, rnvalid, .false.)
				call csv_write(dead_g, rnvalid, .false.)
				call csv_write(dead_g, rnvalid, .false.)
				call csv_write(dead_g, rnvalid, .false.)
				call csv_write(dead_g, rnvalid, .false.)
                call csv_write(dead_g, rnvalid, .false.)
				call csv_write(dead_g, rnvalid, .true.)

			endif
		enddo

		!write plot-level data if turned on
		if (plot_level_data) then

			do ip = 1, site%numplots
				do is = 1, species_pres%numgenera
					call csv_write(dead_pg, iwmo,.false.)
					call csv_write(dead_pg, site%runID, .false.)
					call csv_write(dead_pg, year, .false.)
					call csv_write(dead_pg, ip, .false.)
					call csv_write(dead_pg,                                    &
						trim(adjustl(species_pres%genusgroups(is))), .false.)

					if (in_site(is)) then

						d_markers(ip, is, :) = d_markers(ip, is, :)*plotscale

						call csv_write(dead_pg, d_markers(ip, is, 1), .false.)
						call csv_write(dead_pg, d_markers(ip, is, 2), .false.)
						call csv_write(dead_pg, d_markers(ip, is, 3), .false.)
						call csv_write(dead_pg, d_markers(ip, is, 4), .false.)
						call csv_write(dead_pg, d_markers(ip, is, 5), .false.)
						call csv_write(dead_pg, d_markers(ip, is, 6), .false.)
						call csv_write(dead_pg, d_markers(ip, is, 7), .false.)
						call csv_write(dead_pg, d_markers(ip, is, 8), .false.)
						call csv_write(dead_pg, mean_dbh(ip, is), .false.)
						call csv_write(dead_pg, plotlevel_biomc(ip, is), .true.)

					else
						call csv_write(dead_pg, rnvalid, .false.)
						call csv_write(dead_pg, rnvalid, .false.)
						call csv_write(dead_pg, rnvalid, .false.)
						call csv_write(dead_pg, rnvalid, .false.)
						call csv_write(dead_pg, rnvalid, .false.)
						call csv_write(dead_pg, rnvalid, .false.)
						call csv_write(dead_pg, rnvalid, .false.)
						call csv_write(dead_pg, rnvalid, .false.)
						call csv_write(dead_pg, rnvalid, .false.)
						call csv_write(dead_pg, rnvalid, .true.)
					endif

				enddo
			enddo
		endif

	end subroutine write_dead_genus_data

	!:.........................................................................:

	subroutine write_dead_species_data(site, species_pres, year)
		!writes dead species-level data to Dead_Species_Data.csv and
		  !to Plotlevel_Dead_Species.csv if turned on
		!Authors: Adrianna Foster and Jacquelyn Shuman 2015, v. 1.0
		!Inputs:
		!	site:         site instance
		!	species_pres: species present at site
		!	year:         simulation year

		type(SiteData),   intent(in)                :: site
		type(Groups),     intent(in)                :: species_pres
		integer,          intent(in)                :: year

		real,    dimension(site%numplots, species_pres%numspecies, 8) :: d_markers

		real,    dimension(site%numplots, species_pres%numspecies) :: biomC,   &
                                                                      plotlevel_biomc, &
                                                                      mean_dbh


		real,    dimension(species_pres%numspecies, 8)   :: total_mort_markers
		real,    dimension(species_pres%numspecies)                            &
												         :: total_biomc,       &
                                                            plot_mean_biomc,   &
                                                            plot_std_biomc,    &
                                                            plot_mean_dbh,     &
                                                            plot_std_dbh

		integer, dimension(species_pres%numspecies)      :: spec_index
		logical, dimension(species_pres%numspecies)      :: in_site
		character(len = 8)                               :: field = 'species'
		integer                                          :: iwmo
		integer                                          :: ip, is, ns, il

		plotscale = hec_to_m2/plotsize
		plotadj = plotscale/float(site%numplots)
		plotrenorm = 1.0/plotsize/float(site%numplots)

		iwmo = site%site_id

		!check which species are in site
		in_site = .false.
		do is = 1, species_pres%numspecies
			do ip = 1, site%numplots
				do ns = 1, size(site%plots(ip)%species)
					if (species_pres%spec_names(is, 2) .eq.                    &
							site%plots(ip)%species(ns)%unique_id) then
						in_site(is) = .true.
						spec_index(is) = ns
						exit
					endif
				enddo
			enddo
		enddo

		!get individual plot data
		do ip = 1, site%numplots
			call sum_over_deadsg(site%plots(ip),                               &
				species_pres%spec_names(:, 2), field, biomC(ip, :),            &
				mean_dbh(ip, :), d_markers(ip, :, :))
		enddo

		total_biomc = 0.0
		plotlevel_biomc = 0.0
		total_mort_markers = 0.0

		!accumulate across plots
		do is = 1, species_pres%numspecies
			if (in_site(is)) then
				ns = spec_index(is)
				do ip = 1, site%numplots
					total_biomc(is) = total_biomc(is) + biomc(ip, is)

					do il = 1, 8
						total_mort_markers(is, il) =                           &
							total_mort_markers(is, il) + d_markers(ip, is, il)
					end do
				end do
			end if

			!calculate mean and standard deviations of data
			call stddev(biomc(:, is), plot_mean_biomc(is), plot_std_biomc(is), &
				rnvalid)
			call stddev(mean_dbh(:, is), plot_mean_dbh(is), plot_std_dbh(is),  &
				rnvalid)
		enddo

		!convert values to per ha
		plotlevel_biomc = plotscale*biomc
		plot_mean_biomc = plotscale*plot_mean_biomc
		plot_std_biomc = plotscale*plot_std_biomc

		!write to csv
		do is = 1, species_pres%numspecies
			call csv_write(dead_s, iwmo, .false.)
			call csv_write(dead_s, site%runID, .false.)
			call csv_write(dead_s, year, .false.)
			call csv_write(dead_s,                                             &
				trim(adjustl(species_pres%spec_names(is, 1))), .false.)
			call csv_write(dead_s,                                             &
				trim(adjustl(species_pres%spec_names(is, 2))), .false.)

			if (in_site(is)) then

				total_biomc(is) = total_biomc(is)*plotadj
				total_mort_markers(is, :) =  total_mort_markers(is, :)*plotadj

				call csv_write(dead_s, total_mort_markers(is, 1), .false.)
				call csv_write(dead_s, total_mort_markers(is, 2), .false.)
				call csv_write(dead_s, total_mort_markers(is, 3), .false.)
				call csv_write(dead_s, total_mort_markers(is, 4), .false.)
				call csv_write(dead_s, total_mort_markers(is, 5), .false.)
				call csv_write(dead_s, total_mort_markers(is, 6), .false.)
				call csv_write(dead_s, total_mort_markers(is, 7), .false.)
				call csv_write(dead_s, total_mort_markers(is, 8), .false.)

				call csv_write(dead_s, plot_mean_dbh(is), .false.)
				call csv_write(dead_s, total_biomc(is), .false.)
				call csv_write(dead_s, plot_std_biomc(is), .true.)
			else
				call csv_write(dead_s, rnvalid, .false.)
				call csv_write(dead_s, rnvalid, .false.)
				call csv_write(dead_s, rnvalid, .false.)
				call csv_write(dead_s, rnvalid, .false.)
				call csv_write(dead_s, rnvalid, .false.)
				call csv_write(dead_s, rnvalid, .false.)
				call csv_write(dead_s, rnvalid, .false.)
				call csv_write(dead_s, rnvalid, .false.)
				call csv_write(dead_s, rnvalid, .false.)
				call csv_write(dead_s, rnvalid, .false.)
				call csv_write(dead_s, rnvalid, .true.)
			endif

		enddo

		!write plot-level data if turned on
		if (plot_level_data) then

			do ip = 1, site%numplots
				do is = 1, species_pres%numspecies
					call csv_write(dead_ps, iwmo, .false.)
					call csv_write(dead_ps, site%runID, .false.)
					call csv_write(dead_ps, year, .false.)
					call csv_write(dead_ps, ip, .false.)
					call csv_write(dead_ps,                                    &
						trim(adjustl(species_pres%spec_names(is, 1))), .false.)
					call csv_write(dead_ps,                                    &
						trim(adjustl(species_pres%spec_names(is, 2))), .false.)

					if (in_site(is)) then

						d_markers(ip, is, :) = d_markers(ip, is, :)*plotscale

						call csv_write(dead_ps, d_markers(ip, is, 1), .false.)
						call csv_write(dead_ps, d_markers(ip, is, 2), .false.)
						call csv_write(dead_ps, d_markers(ip, is, 3), .false.)
						call csv_write(dead_ps, d_markers(ip, is, 4), .false.)
						call csv_write(dead_ps, d_markers(ip, is, 5), .false.)
						call csv_write(dead_ps, d_markers(ip, is, 6), .false.)
						call csv_write(dead_ps, d_markers(ip, is, 7), .false.)
						call csv_write(dead_ps, d_markers(ip, is, 8), .false.)

						call csv_write(dead_ps, mean_dbh(ip, is), .false.)
						call csv_write(dead_ps, plotlevel_biomc(ip, is), .true.)
					else
						call csv_write(dead_ps, rnvalid, .false.)
						call csv_write(dead_ps, rnvalid, .false.)
						call csv_write(dead_ps, rnvalid, .false.)
						call csv_write(dead_ps, rnvalid, .false.)
						call csv_write(dead_ps, rnvalid, .false.)
						call csv_write(dead_ps, rnvalid, .false.)
						call csv_write(dead_ps, rnvalid, .false.)
						call csv_write(dead_ps, rnvalid, .false.)
						call csv_write(dead_ps, rnvalid, .false.)
						call csv_write(dead_ps, rnvalid, .true.)

					endif

				enddo
			enddo
		endif

	end subroutine write_dead_species_data

	!:.........................................................................:

	subroutine write_site_data(site, year)
		!writes yearly climate data to Climate.csv
		!Author: Katherine Holcomb 2012, v. 1.0
		!Inputs:
		!	site: site instance
		!	year: simulation year

		class(SiteData), intent(in) :: site
		integer,         intent(in) :: year
		integer                     :: iwmo

		iwmo = site%site_id

		call csv_write(clim_unit, iwmo, .false.)
		call csv_write(clim_unit, site%runID, .false.)
		call csv_write(clim_unit, year, .false.)
		call write_site_csv(site, clim_unit)

	end subroutine write_site_data

	!:.........................................................................:

	subroutine write_soiln_data(site, year)
		!writes average yearly soil decomposition variables, averaged by
		  !plot, to SoilDecomp.csv
		!Author: Adrianna Foster 2018, v. 1.0
		!Inputs:
		!	site: site instance
		!	year: simulation year

		class(SiteData), intent(in) :: site
		integer,         intent(in) :: year

		real, dimension(numplots)     :: O_depth, M_depth, moss_biom
		real, dimension(numplots)     :: active, OM, OM_N, avail_n
		real, dimension(numplots, 18) :: litter_in, forest_litter
		real                          :: O_depth_mn, O_depth_sd
		real                          :: M_depth_mn, M_depth_sd
		real                          :: moss_biom_mn, moss_biom_sd
		real                          :: active_mn, active_sd, OM_mn
		real                          :: OM_sd, OM_N_mn, OM_N_sd
		real                          :: avail_n_mn, avail_n_sd
		real, dimension(18)           :: litter_in_mn, litter_in_sd
		real, dimension(18)           :: forlitter_mn, forlitter_sd
		integer                       :: ip, iwmo, il

		do ip = 1, numplots
			O_depth(ip) = site%plots(ip)%soil%O_depth
			M_depth(ip) = site%plots(ip)%soil%M_depth
			!convert from kg/plot to t/ha
			moss_biom(ip) = site%plots(ip)%soil%moss_biom/1000/plotsize*10000
			active(ip) = site%plots(ip)%soil%active
			OM(ip) = site%plots(ip)%soil%cohorts(1, 1)
			OM_N(ip) = site%plots(ip)%soil%cohorts(1, 2)
			avail_n(ip) = site%plots(ip)%soil%avail_N

			do il = 1, 18
				litter_in(ip, il) = site%plots(ip)%soil%litter(il)
				forest_litter(ip, il) =                                        &
					site%plots(ip)%soil%forest_litter(il, 2)
			end do

		end do

		call stddev(O_depth, O_depth_mn, O_depth_sd, rnvalid)
		call stddev(M_depth, M_depth_mn, M_depth_sd, rnvalid)
		call stddev(moss_biom, moss_biom_mn, moss_biom_sd, rnvalid)
		call stddev(active, active_mn, active_sd, rnvalid)
		call stddev(OM, OM_mn, OM_sd, rnvalid)
		call stddev(OM_N, OM_N_mn, OM_N_sd, rnvalid)
		call stddev(avail_n, avail_n_mn, avail_n_sd, rnvalid)

		do il = 1, 18

			call stddev(litter_in(:,il), litter_in_mn(il), litter_in_sd(il))
			call stddev(forest_litter(:,il), forlitter_mn(il), forlitter_sd(il))
		end do

		iwmo = site%site_id

		call csv_write(soildecomp, iwmo, .false.)
		call csv_write(soildecomp, site%runID, .false.)
		call csv_write(soildecomp, year, .false.)
		call csv_write(soildecomp, O_depth_mn, .false.)
		call csv_write(soildecomp, O_depth_sd, .false.)
		call csv_write(soildecomp, M_depth_mn, .false.)
		call csv_write(soildecomp, moss_biom_mn, .false.)
		call csv_write(soildecomp, active_mn, .false.)
		call csv_write(soildecomp, OM_mn, .false.)
		call csv_write(soildecomp, OM_N_mn, .false.)

		do il = 1, 18
			call csv_write(soildecomp, forlitter_mn(il), .false.)
		end do

		call csv_write(soildecomp, avail_n_mn, .true.)

	end subroutine write_soiln_data

	!:.........................................................................:

	subroutine write_tree_data(site, year)
		!writes tree-level data to TreeData.csv
		!Author: Katherine Holcomb 2012, v. 1.0
		!Inputs:
		!	site: site instance
		!	year: simulation year

		class(SiteData), intent(in) :: site
		integer,         intent(in) :: year
		integer                     :: iwmo
		integer                     :: ip, it

		iwmo = site%site_id

		do ip = 1, site%numplots
			do it = 1, site%plots(ip)%numtrees
				call csv_write(tld, iwmo, .false.)
				call csv_write(tld, site%runID, .false.)
				call csv_write(tld, year,.false.)
				call csv_write(tld, ip, .false.)
				call csv_write(tld, it, .false.)
				call write_tree_csv(site%plots(ip)%trees(it), tld)
			enddo
		enddo

	end subroutine write_tree_data

	!:.........................................................................:

end module Output
