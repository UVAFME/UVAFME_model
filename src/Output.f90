module Output

!*******************************************************************************
  !
  ! This module initializes and writes data to output files.
  !
!*******************************************************************************

  use Constants
  use Parameters
  use Tree
  use Plot
  use GenusGroups
  use IO
  use csv_file

  implicit none

contains

    !:.........................................................................:

    subroutine initialize_outputFiles()
        !
        !  Opens output files and writes their headers
        !

        call open_outputFiles
        call write_headers()

    end subroutine initialize_outputFiles

    !:.........................................................................:

    subroutine total_plot_values(site, year)
        !
        !  Writes across-species data to an output file
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/14     A. C. Foster
        !                 J. K. Shuman        Original Code

        ! Data dictionary: constants
        integer, parameter :: LAI_BINS = 12 ! Number of LAI bins

        ! Data dictionary: calling arguments
        type(SiteData), intent(in) :: site ! Site object
        integer,        intent(in) :: year ! Year of simulation

        ! Data dictionary: local variables
        real, dimension(:,:), allocatable      :: ht           ! Tree height (m)
        real, dimension(:,:), allocatable      :: basal_area   ! Tree basal area (m2)
        real, dimension(:,:), allocatable      :: age          ! Tree age (years)
        real, dimension(site%numplots, FC_NUM) :: envresp_1    ! Growth response for trees < 10 m height (0-1)
        real, dimension(site%numplots, FC_NUM) :: envresp_2    ! Growth response for trees >= 10 m and < 20 m height (0-1)
        real, dimension(site%numplots, FC_NUM) :: envresp_3    ! Growth response for trees >= 20 m height (0-1)
        real, dimension(site%numplots)         :: loreys_ht    ! Lorey's Height (m)
        real, dimension(site%numplots)         :: plbasal_area ! Plot-level basal area (m2)
        real, dimension(site%numplots)         :: plbiomC      ! Plot-level biomass (tC)
        real, dimension(site%numplots)         :: plstems      ! Plot-level stem density
        real, dimension(site%numplots)         :: plstandage   ! Plot stand age (years)
        real, dimension(site%numplots)         :: plmaxheight  ! Plot-level maximum height (m)
        real, dimension(site%numplots)         :: sm_stems     ! Plot-level stems <= 5 cm DBH (#)
        real, dimension(site%numplots)         :: med_stems    ! Plot-level stems > 5 cm and <= 20 cm DBH (#)
        real, dimension(site%numplots)         :: lg_stems     ! Plot-level stems > 20 cm DBH (#)
        real, dimension(site%numplots)         :: n_1          ! Trees < 10 m height (#)
        real, dimension(site%numplots)         :: n_2          ! Trees >= 10 m and < 20 m height (#)
        real, dimension(site%numplots)         :: n_3          ! Trees >= 20 m height (#)
        real, dimension(LAI_BINS)              :: site_lai     ! LAI binned by height (m2/m2)
        real, dimension(FC_NUM)                :: envresp_1_mn ! Average growth response for height bin 1 (0-1)
        real, dimension(FC_NUM)                :: envresp_2_mn ! Average growth response for height bin 2 (0-1)
        real, dimension(FC_NUM)                :: envresp_3_mn ! Average growth response for height bin 3 (0-1)
        real, dimension(FC_NUM)                :: envresp_1_sd ! Sd of growth response for height bin 1 (0-1)
        real, dimension(FC_NUM)                :: envresp_2_sd ! Sd of growth response for height bin 2 (0-1)
        real, dimension(FC_NUM)                :: envresp_3_sd ! Sd of growth response for height bin 3 (0-1)
        real, dimension(M_TYPES)               :: mort_markers ! Biomass killed by different stressors (tC/ha)
        real                                   :: loreys_ht_mn ! Average Lorey's Height (m)
        real                                   :: loreys_ht_sd ! SD of Lorey's Height (m)
        real                                   :: totbiomC     ! Average biomass (tC/ha)
        real                                   :: totbiomC_sd  ! SD of  biomass (tC/ha)
        real                                   :: totbasal     ! Average basal area (m2/ha)
        real                                   :: totbasal_sd  ! SD of basal area (m2/ha)
        real                                   :: totstems     ! Average stems (trees/ha)
        real                                   :: totstems_sd  ! SD of stems (trees/ha)
        real                                   :: age_mn       ! Average stand age (years)
        real                                   :: age_sd       ! SD of stand age (years)
        real                                   :: max_ht       ! Average max height (m)
        real                                   :: max_ht_sd    ! SD of max height (m)
        real                                   :: lg_stems_mn  ! Average stems <= 5 cm DBH (trees/ha)
        real                                   :: lg_stems_sd  ! SD of stems <= 5 cm DBH (trees/ha)
        real                                   :: med_stems_mn ! Average stems > 5 cm and <= 20 cm DBH (trees/ha)
        real                                   :: med_stems_sd ! SD of stems > 5 cm and <= 20 cm DBH (trees/ha)
        real                                   :: sm_stems_mn  ! Average stems > 20 cm DBH (trees/ha)
        real                                   :: sm_stems_sd  ! SD of stems > 20 cm DBH (trees/ha)
        integer                                :: laibin_size  ! Size of LAI bins (m)
        integer                                :: ip, it, i    ! Looping indices
        integer                                :: start        ! Start of LAI bin location
        integer                                :: fin          ! End of LAI bin location

        ! Initialize accumulators
        loreys_ht = 0.0
        plbasal_area = 0.0
        plbiomC = 0.0
        plstems = 0.0
        plstandage = 0.0
        plmaxheight = 0.0
        mort_markers = 0.0
        envresp_1 = 0.0
        envresp_2 = 0.0
        envresp_3 = 0.0
        lg_stems = 0.0
        sm_stems = 0.0
        med_stems = 0.0
        site_lai = 0.0
        totbiomC = 0.0
        basal_area = 0.0

        allocate(ht(site%numplots, maxcells*maxcells))
        allocate(basal_area(site%numplots, maxcells*maxcells))
        allocate(age(site%numplots, maxcells*maxcells))

        ! Bin LAI into height sections
        laibin_size = maxheight/LAI_BINS
        do i = 0, LAI_BINS-1
            start = i*laibin_size + 1
            fin = i*laibin_size + laibin_size
            site_lai(i + 1) = sum(site%lai_array(start:fin))
        end do

        do ip = 1, site%numplots

            do it = 1, site%plots(ip)%numtrees

                ! Get number of small, medium and large stems
                if (site%plots(ip)%trees(it)%diam_bht <= 5.0) then
                    sm_stems(ip) = sm_stems(ip) + 1.0
                else if (site%plots(ip)%trees(it)%diam_bht > 5.0 .and.         &
                    site%plots(ip)%trees(it)%diam_bht <= 20.0) then
                    med_stems(ip) = med_stems(ip) + 1.0
                else if (site%plots(ip)%trees(it)%diam_bht > 20.0) then
                    lg_stems(ip) = lg_stems(ip) + 1.0
                end if

                ! Grab individual tree height (m)
                ht(ip, it) = site%plots(ip)%trees(it)%forska_ht

                ! Calculate individual tree basal area (m2)
                basal_area(ip, it) = CM_TO_M*CM_TO_M*0.25*pi*                  &
                    site%plots(ip)%trees(it)%diam_bht**2

                ! Sum up plot-level aboveground biomass (tC)
                plbiomC(ip) = plbiomC(ip) + site%plots(ip)%trees(it)%branchC + &
                    site%plots(ip)%trees(it)%stemC +                           &
                    site%plots(ip)%trees(it)%leaf_bm

                ! Sum up plot-level basal area (m2)
                plbasal_area(ip) = plbasal_area(ip) + basal_area(ip, it)

                ! Grab individual tree age (years)
                age(ip, it) = site%plots(ip)%trees(it)%tree_age

                ! Get maximum height on plot (m)
                plmaxheight(ip) = max(plmaxheight(ip),                         &
                    site%plots(ip)%trees(it)%forska_ht)

                ! Sum up growth response to each stressor by height class
                if (site%plots(ip)%trees(it)%forska_ht < 10.0) then

                    do i = 1, FC_NUM
                        envresp_1(ip, i) = envresp_1(ip, 1) +                  &
                            site%plots(ip)%trees(it)%env_resp(i)
                    end do
                    n_1(ip) = n_1(ip) + 1.0

                else if ((site%plots(ip)%trees(it)%forska_ht >= 10.0)        &
                    .and. (site%plots(ip)%trees(it)%forska_ht < 20.0)) then

                    do i = 1, FC_NUM
                        envresp_2(ip, i) = envresp_2(ip, 1) +                  &
                            site%plots(ip)%trees(it)%env_resp(i)
                    end do
                    n_2(ip) = n_2(ip) + 1.0

                else if (site%plots(ip)%trees(it)%forska_ht .ge. 20.0) then

                    do i = 1, FC_NUM
                        envresp_3(ip, i) = envresp_3(ip, 1) +                  &
                            site%plots(ip)%trees(it)%env_resp(i)
                    end do
                    n_3(ip) = n_3(ip) + 1.0

                end if
            enddo

            ! Lorey's Height calculation
            if (site%plots(ip)%numtrees == 0.0) then
                loreys_ht(ip) = 0.0
            else
                loreys_ht(ip) = sum((basal_area(ip, :))*(ht(ip, :)))/          &
                    sum(basal_area(ip, :))
            endif

            ! Plot level stems
            plstems(ip) = site%plots(ip)%numtrees

            ! Plot level stand - age of oldest tree
            plstandage(ip) = site%plots(ip)%stand_age

            ! Mortality markers - biomass killed per stressor
            do i = 1, M_TYPES
                mort_markers(i) = mort_markers(i) + site%plots(ip)%d_type(i)
            end do

            ! Get average of each stressor type
            if (n_1(ip) >= 1.0) then
                envresp_1(ip, :) = envresp_1(ip, :)/n_1(ip)
            else
                envresp_1(ip, :) = RNVALID
            end if

            if (n_2(ip) >= 1.0) then
                envresp_2(ip, :) = envresp_2(ip, :)/n_2(ip)
            else
                envresp_2(ip, :) = RNVALID
            end if

            if (n_3(ip) >= 1.0) then
                envresp_3(ip, :) = envresp_3(ip, :)/n_3(ip)
            else
                envresp_3(ip, :) = RNVALID
            end if

        end do

        ! Get average of mortality bins and convert to tC/ha
        mort_markers(:) = mort_markers(:)/float(site%numplots)/plotsize*       &
            HEC_TO_M2

        ! Compute means and standard deviations
        call stddev(loreys_ht, loreys_ht_mn, loreys_ht_sd, RNVALID)
        call stddev(plbiomC, totbiomC, totbiomC_sd, RNVALID)
        call stddev(plbasal_area, totbasal, totbasal_sd, RNVALID)
        call stddev(plstems, totstems, totstems_sd, RNVALID)
        call stddev(plmaxheight, max_ht, max_ht_sd, RNVALID)
        call stddev(plstandage, age_mn, age_sd, RNVALID)
        call stddev(lg_stems, lg_stems_mn, lg_stems_sd, RNVALID)
        call stddev(med_stems, med_stems_mn, med_stems_sd, RNVALID)
        call stddev(sm_stems, sm_stems_mn, sm_stems_sd, RNVALID)

        do i = 1, FC_NUM
            call stddev(envresp_1(:, i), envresp_1_mn(i), envresp_1_sd(i),     &
                RNVALID)
            call stddev(envresp_2(:, i), envresp_2_mn(i), envresp_2_sd(i),     &
                RNVALID)
            call stddev(envresp_3(:, i), envresp_3_mn(i), envresp_3_sd(i),     &
                RNVALID)
        end do

        ! Convert units
        totbiomC = totbiomC*HEC_TO_M2/plotsize         ! tC/ha
        totbiomC_sd = totbiomC_sd*HEC_TO_M2/plotsize   ! tC/ha
        totbasal = totbasal*HEC_TO_M2/plotsize         ! m2/ha
        totbasal_sd = totbasal_sd*HEC_TO_M2/plotsize   ! m2/ha
        totstems = totstems*HEC_TO_M2/plotsize         ! stems/ha
        totstems_sd = totstems_sd*HEC_TO_M2/plotsize   ! stems/ha
        lg_stems_mn = lg_stems_mn*HEC_TO_M2/plotsize   ! stems/ha
        lg_stems_sd = lg_stems_sd*HEC_TO_M2/plotsize   ! stems/ha
        med_stems_mn = med_stems_mn*HEC_TO_M2/plotsize ! stems/ha
        med_stems_sd = med_stems_sd*HEC_TO_M2/plotsize ! stems/ha
        sm_stems_mn = sm_stems_mn*HEC_TO_M2/plotsize   ! stems/ha
        sm_stems_sd = sm_stems_sd*HEC_TO_M2/plotsize   ! stems/ha

        ! Write to csv
        call csv_write(plotvals, site%site_id, .false.)
        call csv_write(plotvals, site%runID, .false.)
        call csv_write(plotvals, year, .false.)
        do i = 1, M_TYPES
            call csv_write(plotvals, mort_markers(i), .false.)
        end do
        do i = 1, FC_NUM
            call csv_write(plotvals, envresp_1_mn(i), .false.)
            call csv_write(plotvals, envresp_2_mn(i), .false.)
            call csv_write(plotvals, envresp_3_mn(i), .false.)
        end do
        call csv_write(plotvals, loreys_ht_mn, .false.)
        call csv_write(plotvals, loreys_ht_sd, .false.)
        call csv_write(plotvals, max_ht, .false.)
        call csv_write(plotvals, max_ht_sd, .false.)
        call csv_write(plotvals, totbiomC, .false.)
        call csv_write(plotvals, totbiomC_sd, .false.)
        call csv_write(plotvals, totbasal, .false.)
        call csv_write(plotvals, totbasal_sd, .false.)
        call csv_write(plotvals, totstems, .false.)
        call csv_write(plotvals, totstems_sd, .false.)
        call csv_write(plotvals, sm_stems_mn, .false.)
        call csv_write(plotvals, sm_stems_sd, .false.)
        call csv_write(plotvals, med_stems_mn, .false.)
        call csv_write(plotvals, med_stems_sd, .false.)
        call csv_write(plotvals, lg_stems_mn, .false.)
        call csv_write(plotvals, lg_stems_sd, .false.)
        call csv_write(plotvals, age_mn, .false.)
        call csv_write(plotvals, age_sd, .false.)
        do i = 1, LAI_BINS-1
            call csv_write(plotvals, site_lai(i), .false.)
        end do
        call csv_write(plotvals, site_lai(LAI_BINS), .true.)

    end subroutine total_plot_values

  !:...........................................................................:

    subroutine write_genus_or_species_data(site, species_pres, year, field,    &
        num_types, funit, funit_p)
        !
        !  Writes dead species- or genus- level data to output
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/12     K. Holcomb          Original code
        !    01/29/21     A. C. Foster        Update to write either species
        !                                         or genus data

        ! Data dictionary: calling arguments
        type(SiteData),          intent(in)           :: site         ! Site object
        type(Groups),            intent(in)           :: species_pres ! List of genus and species names
        integer,                 intent(in)           :: year         ! Simulation year
        character(len = 8),      intent(in)           :: field        ! 'genus' or 'species'
        integer,                 intent(in)           :: num_types    ! Number of species or genera
        integer,                 intent(in)           :: funit        ! File unit for output file
        integer,                 intent(in)           :: funit_p      ! File unit for plot-level output file

        ! Data dictionary: local variables
        character(len=MAX_NLEN), dimension(num_types)        :: genera        ! Species/genus names
        integer, dimension(site%numplots, num_types, NHC)    :: diam_cats     ! Size structure (stems/dbh bin)
        real,    dimension(site%numplots, num_types, NHC)    :: biom_cats     ! Biomass by dbh bin (tC)
        real,    dimension(site%numplots, num_types, FC_NUM) :: env_resp      ! Growth response to stressors (0-1)
        real,    dimension(site%numplots, num_types)         :: basal_area    ! Basal area (m2)
        real,    dimension(site%numplots, num_types)         :: biomC         ! Aboveground biomass (tC)
        real,    dimension(site%numplots, num_types)         :: biomN         ! N content (tC)
        real,    dimension(site%numplots, num_types)         :: leaf_bm       ! Leaf biomass (tC)
        real,    dimension(site%numplots, num_types)         :: max_ht        ! Maximum height (m)
        real,    dimension(site%numplots, num_types)         :: max_diam      ! Maximum diameter (cm)
        real,    dimension(site%numplots, num_types)         :: mean_diam     ! Average diameter (cm)
        real,    dimension(site%numplots, num_types)         :: mean_age      ! Average age
        real,    dimension(site%numplots, num_types)         :: biomC_lg      ! Biomass > 9 cm DBH (tC)
        real,    dimension(site%numplots, num_types)         :: biomC_sm      ! Biomass < 9 cm DBH (tC)
        real,    dimension(site%numplots, num_types)         :: basal_lg      ! Basal area > 9 cm DBH (m2)
        real,    dimension(site%numplots, num_types)         :: basal_sm      ! Basal area < 9 cm DBH (m2)
        real,    dimension(site%numplots, num_types)         :: dens_lg       ! Stems > 9 cm DBH (#)
        real,    dimension(site%numplots, num_types)         :: dens_sm       ! Stems < 9 cm DBH (#)
        real,    dimension(site%numplots, num_types)         :: dbh_lg        ! Average dbh > 9 cm DBH (cm)
        real,    dimension(site%numplots, num_types)         :: dbh_sm        ! Average dbh < 9 cm DBH (cm)
        real,    dimension(num_types, FC_NUM)                :: env_resp_mn   ! Average growth response (0-1)
        real,    dimension(num_types, FC_NUM)                :: env_resp_sd   ! SD of growth response (0-1)
        real,    dimension(num_types, NHC)                   :: tot_diam_cats ! Average size structure (stems/ha)
        real,    dimension(num_types, NHC)                   :: tot_biom_cats ! Average biomass by bin (tC/ha)
        real,    dimension(num_types)                        :: tot_basal     ! Average basal area (m2/ha)
        real,    dimension(num_types)                        :: tot_basal_sd  ! SD of basal area (m2/ha)
        real,    dimension(num_types)                        :: tot_biomC     ! Average aboveground biomass (tC/ha)
        real,    dimension(num_types)                        :: tot_biomC_sd  ! SD aboveground biomass (tC/ha)
        real,    dimension(num_types)                        :: tot_biomN     ! Average N content (tN/ha)
        real,    dimension(num_types)                        :: tot_biomN_sd  ! SD of N content (tN/ha)
        real,    dimension(num_types)                        :: tot_leafbm    ! Average leaf biomass (tC)
        real,    dimension(num_types)                        :: tot_leafbm_sd ! SD of leaf biomass (tC)
        real,    dimension(num_types)                        :: max_ht_mn     ! Average maximum height (m)
        real,    dimension(num_types)                        :: max_ht_sd     ! SD of maximum height (m)
        real,    dimension(num_types)                        :: max_diam_mn   ! Average maximum diameter (cm)
        real,    dimension(num_types)                        :: max_diam_sd   ! SD of maximum diameter (cm)
        real,    dimension(num_types)                        :: mean_diam_mn  ! Average diameter (cm)
        real,    dimension(num_types)                        :: mean_diam_sd  ! SD of diameter (cm)
        real,    dimension(num_types)                        :: age_mn        ! Average age (years)
        real,    dimension(num_types)                        :: age_sd        ! SD of age (years)
        real,    dimension(num_types)                        :: biomC_lg_mn   ! Average biomass > 9 cm DBH (tC/ha)
        real,    dimension(num_types)                        :: biomC_lg_sd   ! SD biomass > 9 cm DBH (tC/ha)
        real,    dimension(num_types)                        :: biomC_sm_mn   ! Average basal < 9 cm DBH (tC/ha)
        real,    dimension(num_types)                        :: biomC_sm_sd   ! SD biomass < 9 cm DBH (tC/ha)
        real,    dimension(num_types)                        :: basal_lg_mn   ! Average basal area > 9 cm DBH (m2/ha)
        real,    dimension(num_types)                        :: basal_lg_sd   ! SD basal_area > 9 cm DBH (m2/ha)
        real,    dimension(num_types)                        :: basal_sm_mn   ! Average basal area < 9 cm DBH (m2/ha)
        real,    dimension(num_types)                        :: basal_sm_sd   ! SD basal area < 9 cm DBH (m2/ha)
        real,    dimension(num_types)                        :: dens_lg_mn    ! Average stems > 9 cm DBH (stems/ha)
        real,    dimension(num_types)                        :: dens_lg_sd    ! SD stems > 9 cm DBH (stems/ha)
        real,    dimension(num_types)                        :: dens_sm_mn    ! Average stems < 9 cm DBH (stems/ha)
        real,    dimension(num_types)                        :: dens_sm_sd    ! SD stems > 9 cm DBH (stems/ha)
        real,    dimension(num_types)                        :: dbh_lg_mn     ! Average diameter > 9 cm DBH (cm)
        real,    dimension(num_types)                        :: dbh_lg_sd     ! SD of diameter > 9 cm DBH (cm)
        real,    dimension(num_types)                        :: dbh_sm_mn     ! Average diameter < 9 cm DBH (cm)
        real,    dimension(num_types)                        :: dbh_sm_sd     ! SD of diameter < 9 cm DBH (cm)
        integer, dimension(num_types)                        :: spec_index    ! Location of present species/genera
        logical, dimension(num_types)                        :: in_site       ! Is the species/genus in the site?
        character(len=MAX_NLEN)                              :: comp          ! Species/genus name
        character(len=MAX_NLEN)                              :: genname       ! Genus name
        integer                                              :: ip, is, ns, l ! Looping indices

        ! Get either genera or species
        if (field == 'genus') then
            genera = species_pres%genusgroups
        else if (field == 'species') then
            genera = species_pres%spec_names(:, 2)
        end if

        ! Check to see which species/genera are in the site
        in_site = .false.
        do is = 1, num_types
            do ip = 1, site%numplots
                do ns = 1, size(site%species)
                    if (field == 'genus') then
                        comp = site%species(ns)%genus_name
                    else
                        comp = site%species(ns)%unique_id
                    end if
                    if (genera(is) == comp) then
                        spec_index(is) = ns
                        in_site(is) = .true.
                        exit
                    endif
                enddo
            enddo
        enddo

        ! Get individual plot data for each genus
        do ip = 1, site%numplots
            call sum_over_sg(site%plots(ip), genera, field, basal_area(ip, :), &
                leaf_bm(ip, :), biomC(ip, :), biomN(ip, :), max_ht(ip, :),     &
                max_diam(ip, :), mean_diam(ip, :), mean_age(ip, :),            &
                biom_cats(ip, :, :), env_resp(ip, :, :), basal_lg(ip, :),      &
                basal_sm(ip, :), dbh_lg(ip, :), dbh_sm(ip, :), dens_lg(ip, :), &
                dens_sm(ip, :), biomC_lg(ip, :), biomC_sm(ip, :))

            call tree_dm_cats(site%plots(ip), genera, field,                   &
                diam_cats(ip, :, :))
        enddo


        tot_biomC = RNVALID
        tot_biomN = RNVALID
        tot_basal = RNVALID
        env_resp_mn = RNVALID
        tot_diam_cats = 0
        tot_biom_cats = 0.0

        ! Accumulate across-plot data
        do is = 1, num_types
            do ip = 1, site%numplots
                if (in_site(is)) then
                    tot_diam_cats(is, :) = tot_diam_cats(is, :) +              &
                        float(diam_cats(ip, is, :))
                    tot_biom_cats(is, :) = tot_biom_cats(is, :) +              &
                        biom_cats(ip, is, :)
                endif
            enddo

            ! Get mean and standard deviation of data
            call stddev(basal_area(:, is), tot_basal(is), tot_basal_sd(is),    &
                RNVALID)
            call stddev(biomC(:, is), tot_biomC(is), tot_biomC_sd(is), RNVALID)
            call stddev(biomN(:, is), tot_biomN(is), tot_biomN_sd(is), &
                RNVALID)
            call stddev(leaf_bm(:, is), tot_leafbm(is), tot_leafbm_sd(is),     &
                RNVALID)
            call stddev(max_ht(:, is), max_ht_mn(is), max_ht_sd(is), RNVALID)
            call stddev(max_diam(:, is), max_diam_mn(is), max_diam_sd(is),     &
                RNVALID)
            call stddev(mean_diam(:, is), mean_diam_mn(is), mean_diam_sd(is),  &
                RNVALID)
            call stddev(mean_age(:, is), age_mn(is), age_sd(is), RNVALID)
            call stddev(biomC_lg(:, is), biomC_lg_mn(is), biomC_lg_sd(is),     &
                RNVALID)
            call stddev(biomC_sm(:, is), biomC_sm_mn(is), biomC_sm_sd(is),     &
                RNVALID)
            call stddev(basal_lg(:, is), basal_lg_mn(is), basal_lg_sd(is),     &
                RNVALID)
            call stddev(basal_sm(:, is), basal_sm_mn(is), basal_sm_sd(is),     &
                RNVALID)
            call stddev(dbh_lg(:, is), dbh_lg_mn(is), dbh_lg_sd(is), RNVALID)
            call stddev(dbh_sm(:, is), dbh_sm_mn(is), dbh_sm_sd(is), RNVALID)
            call stddev(dens_lg(:, is), dens_lg_mn(is), dens_lg_sd(is), RNVALID)
            call stddev(dens_sm(:, is), dens_sm_mn(is), dens_sm_sd(is), RNVALID)

            do l = 1, FC_NUM
                call stddev(env_resp(:, is, l), env_resp_mn(is, l),            &
                    env_resp_sd(is, l), RNVALID)
            end do
        enddo

        ! Convert to correct units
        tot_basal = tot_basal*HEC_TO_M2/plotsize         ! m2/ha
        tot_basal_sd = tot_basal_sd*HEC_TO_M2/plotsize   ! m2/ha
        tot_biomC = tot_biomC*HEC_TO_M2/plotsize         ! tC/ha
        tot_biomC_sd = tot_biomC_sd*HEC_TO_M2/plotsize   ! tC/ha
        tot_biomN = tot_biomN*HEC_TO_M2/plotsize         ! tN/ha
        tot_biomN_sd = tot_biomN_sd*HEC_TO_M2/plotsize   ! tN/ha
        tot_leafbm = tot_leafbm/plotsize                 ! tC/m2
        tot_leafbm_sd = tot_leafbm_sd/plotsize           ! tC/m2
        biomC_lg_mn = biomC_lg_mn*HEC_TO_M2/plotsize     ! tC/ha
        biomC_lg_sd = biomC_lg_sd*HEC_TO_M2/plotsize     ! tC/ha
        biomC_sm_mn = biomC_sm_mn*HEC_TO_M2/plotsize     ! tC/ha
        biomC_sm_sd = biomC_sm_sd*HEC_TO_M2/plotsize     ! tC/ha
        basal_lg_mn = basal_lg_mn*HEC_TO_M2/plotsize     ! m2/ha
        basal_lg_sd = basal_lg_sd*HEC_TO_M2/plotsize     ! m2/ha
        basal_sm_mn = basal_sm_mn*HEC_TO_M2/plotsize     ! m2/ha
        basal_sm_sd = basal_sm_sd*HEC_TO_M2/plotsize     ! m2/ha
        dens_lg_mn = dens_lg_mn*HEC_TO_M2/plotsize       ! trees/ha
        dens_lg_sd = dens_lg_sd*HEC_TO_M2/plotsize       ! trees/ha
        dens_sm_mn = dens_sm_mn*HEC_TO_M2/plotsize       ! trees/ha
        dens_sm_sd = dens_sm_sd*HEC_TO_M2/plotsize       ! trees/ha

        !trees/ha and tC/ha
        tot_diam_cats = tot_diam_cats*HEC_TO_M2/plotsize/float(site%numplots)
        tot_biom_cats = tot_biom_cats*HEC_TO_M2/plotsize/float(site%numplots)

        ! Write to csv
        do is = 1, num_types

            ! We write lines out for every species/genus, but write -999.0s
            ! to ones not 'present'
            call csv_write(funit, site%site_id, .false.)
            call csv_write(funit, site%runID, .false.)
            call csv_write(funit, year, .false.)

            ! Get species/genus name
            if (field == 'genus') then
                comp = trim(adjustl(species_pres%genusgroups(is)))
                call csv_write(funit, comp, .false.)
            else
                genname = trim(adjustl(species_pres%spec_names(is,1)))
                comp = trim(adjustl(species_pres%spec_names(is, 2)))
                call csv_write(funit, genname, .false.)
                call csv_write(funit, comp, .false.)
            end if

            if (in_site(is)) then

                ns = spec_index(is)

                ! Convert to LAI
                tot_leafbm(is) = tot_leafbm(is)/                               &
                    (site%species(ns)%leafarea_c*  2.0)
                tot_leafbm_sd(is) = tot_leafbm_sd(is)/                         &
                    (site%species(ns)%leafarea_c*2.0)

                ! Write out the data
                do l = 1, NHC
                    call csv_write(funit, tot_diam_cats(is, l), .false.)
                end do
                do l = 1, NHC
                    call csv_write(funit, tot_biom_cats(is, l), .false.)
                end do
                do l = 1, FC_NUM
                    call csv_write(funit, env_resp_mn(is, l), .false.)
                end do
                call csv_write(funit, max_diam_mn(is), .false.)
                call csv_write(funit, mean_diam_mn(is), .false.)
                call csv_write(funit, age_mn(is), .false.)
                call csv_write(funit, max_ht_mn(is), .false.)
                call csv_write(funit, tot_leafbm(is), .false.)
                call csv_write(funit, tot_basal(is), .false.)
                call csv_write(funit, tot_basal_sd(is), .false.)
                call csv_write(funit, tot_biomC(is), .false.)
                call csv_write(funit, tot_biomC_sd(is), .false.)
                call csv_write(funit, biomC_lg_mn(is), .false.)
                call csv_write(funit, biomC_lg_sd(is), .false.)
                call csv_write(funit, biomC_sm_mn(is), .false.)
                call csv_write(funit, biomC_sm_sd(is), .false.)
                call csv_write(funit, basal_lg_mn(is), .false.)
                call csv_write(funit, basal_lg_sd(is), .false.)
                call csv_write(funit, basal_sm_mn(is), .false.)
                call csv_write(funit, basal_sm_sd(is), .false.)
                call csv_write(funit, dens_lg_mn(is), .false.)
                call csv_write(funit, dens_lg_sd(is), .false.)
                call csv_write(funit, dens_sm_mn(is), .false.)
                call csv_write(funit, dens_sm_sd(is), .false.)
                call csv_write(funit, dbh_lg_mn(is), .false.)
                call csv_write(funit, dbh_lg_sd(is), .false.)
                call csv_write(funit, dbh_sm_mn(is), .false.)
                call csv_write(funit, dbh_sm_sd(is), .true.)
            else
                ! Write -999.0s
                do l = 1, NHC + NHC + FC_NUM + 24
                    call csv_write(funit, RNVALID, .false.)
                end do
                call csv_write(funit, RNVALID, .true.)
            endif
        enddo

        ! Write plot-level data if turned on
        if (plot_level_data) then

            do ip = 1, site%numplots
                do is = 1, num_types

                    ! We write lines out for every species/genus, but write
                    ! -999.0s to ones not 'present'
                    call csv_write(funit_p, site%site_id, .false.)
                    call csv_write(funit_p, site%runID, .false.)
                    call csv_write(funit_p, year, .false.)
                    call csv_write(funit_p, ip, .false.)

                    ! Get species/genus name
                    if (field == 'genus') then
                        comp = trim(adjustl(species_pres%genusgroups(is)))
                        call csv_write(funit_p, comp, .false.)
                    else
                        genname = trim(adjustl(species_pres%spec_names(is,1)))
                        comp = trim(adjustl(species_pres%spec_names(is, 2)))
                        call csv_write(funit_p, genname, .false.)
                        call csv_write(funit_p, comp, .false.)
                    end if

                    if (in_site(is)) then

                        ! Convert some units
                        diam_cats(ip, is,:) = diam_cats(ip,is,:)*              &
                            HEC_TO_M2/plotsize
                        biom_cats(ip,is,:) =                                   &
                            biom_cats(ip,is,:)*HEC_TO_M2/plotsize
                        basal_area(ip, is) = basal_area(ip, is )*              &
                            HEC_TO_M2/plotsize
                        biomC(ip, is) = biomC(ip, is )*HEC_TO_M2/plotsize

                        do l = 1, NHC
                            call csv_write(funit_p, diam_cats(ip, is, l),      &
                                .false.)
                        end do
                        do l = 1, NHC
                            call csv_write(funit_p, biom_cats(ip, is, l),      &
                                .false.)
                        end do
                        call csv_write(funit_p, max_diam(ip, is), .false.)
                        call csv_write(funit_p, mean_diam(ip, is), .false.)
                        call csv_write(funit_p, max_ht(ip, is), .false.)
                        call csv_write(funit_p, basal_area(ip, is), .false.)
                        call csv_write(funit_p, biomC(ip, is), .true.)

                    else
                        ! Write -999.0s
                        do l = 1, NHC + NHC + 4
                            call csv_write(funit_p, RNVALID, .false.)
                        end do
                        call csv_write(funit_p, RNVALID, .true.)
                    endif
                enddo
            enddo
        endif

    end subroutine write_genus_or_species_data

    !:.........................................................................:

    subroutine write_dead_genus_or_species_data(site, species_pres, year,      &
        field, num_types, funit, funit_p)
        !
        !  Writes dead species- or genus- level data to output
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/15     A. C. Foster        Original code
        !                 J. K. Shuman
        !    01/30/21     A. C. Foster        Update to write either species
        !                                         or genus data

        ! Data dictionary: calling arguments
        type(SiteData),     intent(in) :: site         ! Site object
        type(Groups),       intent(in) :: species_pres ! List of genus and species names
        integer,            intent(in) :: year         ! Simulation year
        character(len = 8), intent(in) :: field        ! 'genus' or 'species'
        integer,            intent(in) :: num_types    ! Number of species or genera
        integer,            intent(in) :: funit        ! File unit for output file
        integer,            intent(in) :: funit_p      ! File unit for plot-level output file

        ! Data dictionary: local variables
        character(len=MAX_NLEN), dimension(num_types)         :: genera        ! Species/genus names
        real,    dimension(site%numplots, num_types, M_TYPES) :: d_markers     ! Causes of mortality (tC)
        real,    dimension(site%numplots, num_types)          :: biomC         ! Aboveground biomass (tC)
        real,    dimension(site%numplots, num_types)          :: mean_diam     ! Average diameter (cm)
        real,    dimension(num_types, M_TYPES)                :: tot_d_markers ! Average causes of mortality (tC/ha)
        real,    dimension(num_types)                         :: tot_biomC     ! Average aboveground biomass (tC/ha)
        real,    dimension(num_types)                         :: tot_biomC_sd  ! SD aboveground biomass (tC/ha)
        real,    dimension(num_types)                         :: mean_diam_mn  ! Average diameter (cm)
        real,    dimension(num_types)                         :: mean_diam_sd  ! SD of diameter (cm)
        integer, dimension(num_types)                         :: spec_index    ! Location of present species/genera
        logical, dimension(num_types)                         :: in_site       ! Is the species/genus in the site?
        character(len=MAX_NLEN)                               :: comp          ! Species/genus name
        character(len=MAX_NLEN)                               :: genname       ! Genus name
        integer                                               :: ip, is, ns, l ! Looping indices

        ! Get either genera or species
        if (field == 'genus') then
            genera = species_pres%genusgroups
        else if (field == 'species') then
            genera = species_pres%spec_names(:, 2)
        end if

        ! Check to see which species/genera are in the site
        in_site = .false.
        do is = 1, num_types
            do ip = 1, site%numplots
                do ns = 1, size(site%species)
                    if (field == 'genus') then
                        comp = site%species(ns)%genus_name
                    else
                        comp = site%species(ns)%unique_id
                    end if
                    if (genera(is) == comp) then
                        spec_index(is) = ns
                        in_site(is) = .true.
                        exit
                    endif
                enddo
            enddo
        enddo

        ! Get individual plot data for each genus
        do ip = 1, site%numplots
            call sum_over_deadsg(site%plots(ip), genera, field, biomC(ip, :),  &
                mean_diam(ip, :), d_markers(ip, :, :))
        enddo

        tot_biomC = RNVALID
        tot_d_markers = 0.0

        ! Accumulate across-plot data
        do is = 1, num_types
            do ip = 1, site%numplots
                if (in_site(is)) then
                    tot_d_markers(is, :) = tot_d_markers(is, :) +              &
                        d_markers(ip, is, :)
                endif
            enddo

            ! Get mean and standard deviation of data
            call stddev(biomC(:, is), tot_biomC(is), tot_biomC_sd(is), RNVALID)
            call stddev(mean_diam(:, is), mean_diam_mn(is), mean_diam_sd(is),  &
                RNVALID)
        enddo

        ! Convert to correct units
        tot_biomC = tot_biomC*HEC_TO_M2/plotsize       ! tC/ha
        tot_biomC_sd = tot_biomC_sd*HEC_TO_M2/plotsize ! tC/ha

        ! tC/ha
        tot_d_markers = tot_d_markers*HEC_TO_M2/plotsize/float(site%numplots)

        ! Write to csv
        do is = 1, num_types

            ! We write lines out for every species/genus, but write -999.0s
            ! to ones not 'present'
            call csv_write(funit, site%site_id, .false.)
            call csv_write(funit, site%runID, .false.)
            call csv_write(funit, year, .false.)

            ! Get species/genus name
            if (field == 'genus') then
                comp = trim(adjustl(species_pres%genusgroups(is)))
                call csv_write(funit, comp, .false.)
            else
                genname = trim(adjustl(species_pres%spec_names(is,1)))
                comp = trim(adjustl(species_pres%spec_names(is, 2)))
                call csv_write(funit, genname, .false.)
                call csv_write(funit, comp, .false.)
            end if

            if (in_site(is)) then

                ! Write out the data
                do l = 1, M_TYPES
                    call csv_write(funit, tot_d_markers(is, l), .false.)
                end do
                call csv_write(funit, mean_diam_mn(is), .false.)
                call csv_write(funit, tot_biomC(is), .false.)
                call csv_write(funit, tot_biomC_sd(is), .true.)
            else
                ! Write -999.0s
                do l = 1, M_TYPES + 2
                    call csv_write(funit, RNVALID, .false.)
                end do
                call csv_write(funit, RNVALID, .true.)
            endif
        enddo

        !write plot-level data if turned on
        if (plot_level_data) then

            do ip = 1, site%numplots
                do is = 1, num_types

                    ! We write lines out for every species/genus, but write
                    ! -999.0s to ones not 'present'
                    call csv_write(funit_p, site%site_id, .false.)
                    call csv_write(funit_p, site%runID, .false.)
                    call csv_write(funit_p, year, .false.)
                    call csv_write(funit_p, ip, .false.)

                    ! Get species/genus name
                    if (field == 'genus') then
                        comp = trim(adjustl(species_pres%genusgroups(is)))
                        call csv_write(funit_p, comp, .false.)
                    else
                        genname = trim(adjustl(species_pres%spec_names(is,1)))
                        comp = trim(adjustl(species_pres%spec_names(is, 2)))
                        call csv_write(funit_p, genname, .false.)
                        call csv_write(funit_p, comp, .false.)
                    end if

                    if (in_site(is)) then

                        ! Convert some units
                        d_markers(ip,is,:) = d_markers(ip,is,:)*               &
                            HEC_TO_M2/plotsize
                        biomC(ip, is) = biomC(ip, is )*HEC_TO_M2/plotsize

                        do l = 1, M_TYPES
                            call csv_write(funit_p, d_markers(ip, is, l),      &
                                .false.)
                        end do
                        call csv_write(funit_p, mean_diam(ip, is), .false.)
                        call csv_write(funit_p, biomC(ip, is), .true.)

                    else
                        ! Write -999.0s
                        do l = 1, M_TYPES + 1
                            call csv_write(funit_p, RNVALID, .false.)
                        end do
                        call csv_write(funit_p, RNVALID, .true.)
                    endif
                enddo
            enddo
        endif

    end subroutine write_dead_genus_or_species_data


    !:.........................................................................:

    subroutine write_site_data(site, year)
        !
        !  Writes site/climate output
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/12     K. Holcomb          Original code

        ! Data dictionary: calling arguments
        class(SiteData), intent(in) :: site ! Site object
        integer,         intent(in) :: year ! Simulation year

        call csv_write(clim_unit, site%site_id, .false.)
        call csv_write(clim_unit, site%runID, .false.)
        call csv_write(clim_unit, year, .false.)
        call write_site_csv(site, clim_unit)

    end subroutine write_site_data

    !:.........................................................................:

    subroutine write_soiln_data(site, year)
        !
        !  Writes soil decomposition variables, averaged by plot
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    10/26/18     A. C. Foster        Original code

        ! Data dictionary: calling arguments
        class(SiteData), intent(in) :: site ! Site object
        integer,         intent(in) :: year ! Simulation year

        ! Data dictionary: local variables
        real, dimension(site%numplots, LIT_LEVS) :: forest_litter ! Total litter (t/ha)
        real, dimension(site%numplots)           :: O_depth       ! Organic layer depth (cm)
        real, dimension(site%numplots)           :: M_depth       ! Moss layer depth (cm)
        real, dimension(site%numplots)           :: moss_biom     ! Moss biomass (kg/ha)
        real, dimension(site%numplots)           :: active        ! Active layer thickness (cm)
        real, dimension(site%numplots)           :: OM            ! Humus content (t/ha)
        real, dimension(site%numplots)           :: OM_N          ! Humus N content (tN/ha)
        real, dimension(site%numplots)           :: avail_n       ! Plant-available N (kgN/ha)
        real, dimension(LIT_LEVS)                :: fresh_lit_mn  ! Average fresh litter (t/ha)
        real, dimension(LIT_LEVS)                :: fresh_lit_sd  ! SD of fresh litter (t/ha)
        real, dimension(LIT_LEVS)                :: forlitter_mn  ! Average litter (t/ha)
        real, dimension(LIT_LEVS)                :: forlitter_sd  ! SD of litter (t/ha)
        real                                     :: O_depth_mn    ! Average organic layer depth (cm)
        real                                     :: O_depth_sd    ! SD of organic layer depth (cm)
        real                                     :: M_depth_mn    ! Average moss depth (cm)
        real                                     :: M_depth_sd    ! SD of moss depth (cm)
        real                                     :: moss_biom_mn  ! Average moss biomass (kg)
        real                                     :: moss_biom_sd  ! SD of moss biomass (kg)
        real                                     :: active_mn     ! Average active layer depth (cm)
        real                                     :: active_sd     ! SD of active layer depth (cm)
        real                                     :: OM_mn         ! Average humus content (t/ha)
        real                                     :: OM_sd         ! SD of humus content (t/ha)
        real                                     :: OM_N_mn       ! Average humus N content (tN/ha)
        real                                     :: OM_N_sd       ! SD of humus N content (tN/ha)
        real                                     :: avail_n_mn    ! Average plant-available N (kgN/ha)
        real                                     :: avail_n_sd    ! SD of plant-available N (kgN/ha)
        integer                                  :: ip, il        ! Looping indices

        ! Pull in data from plots
        do ip = 1, site%numplots

            O_depth(ip) = site%plots(ip)%soil%O_depth*M_TO_CM
            M_depth(ip) = site%plots(ip)%soil%M_depth*M_TO_CM
            moss_biom(ip) = site%plots(ip)%soil%moss_biom/plotsize/M2_TO_HEC
            active(ip) = site%plots(ip)%soil%active*M_TO_CM
            OM(ip) = site%plots(ip)%soil%cohorts(1, 1)
            OM_N(ip) = site%plots(ip)%soil%cohorts(1, 2)
            avail_n(ip) = site%plots(ip)%soil%avail_N*T_TO_KG

            do il = 1, LIT_LEVS
                forest_litter(ip, il) =                                        &
                    site%plots(ip)%soil%forest_litter(il, 1)
            end do

        end do

        ! Get means and standard deviations
        call stddev(O_depth, O_depth_mn, O_depth_sd, RNVALID)
        call stddev(M_depth, M_depth_mn, M_depth_sd, RNVALID)
        call stddev(moss_biom, moss_biom_mn, moss_biom_sd, RNVALID)
        call stddev(active, active_mn, active_sd, RNVALID)
        call stddev(OM, OM_mn, OM_sd, RNVALID)
        call stddev(OM_N, OM_N_mn, OM_N_sd, RNVALID)
        call stddev(avail_n, avail_n_mn, avail_n_sd, RNVALID)
        do il = 1, LIT_LEVS
            call stddev(forest_litter(:,il), forlitter_mn(il), forlitter_sd(il))
        end do

        ! Write to file
        call csv_write(soildecomp, site%site_id, .false.)
        call csv_write(soildecomp, site%runID, .false.)
        call csv_write(soildecomp, year, .false.)
        call csv_write(soildecomp, O_depth_mn, .false.)
        call csv_write(soildecomp, O_depth_sd, .false.)
        call csv_write(soildecomp, M_depth_mn, .false.)
        call csv_write(soildecomp, moss_biom_mn, .false.)
        call csv_write(soildecomp, active_mn, .false.)
        call csv_write(soildecomp, OM_mn, .false.)
        call csv_write(soildecomp, OM_N_mn, .false.)
        do il = 1, LIT_LEVS
            call csv_write(soildecomp, forlitter_mn(il), .false.)
        end do
        call csv_write(soildecomp, avail_n_mn, .true.)

    end subroutine write_soiln_data

    !:.........................................................................:

    subroutine write_tree_data(site, year)
        !
        !  Writes tree level data
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/12     K. Holcomb          Original code

        ! Data dictionary: calling arguments
        class(SiteData), intent(in) :: site ! Site object
        integer,         intent(in) :: year ! Simulation year

        ! Data dictionary: local variables
        integer :: ip, it ! Looping indices

        do ip = 1, site%numplots
            do it = 1, site%plots(ip)%numtrees
                call csv_write(pl_tree, site%site_id, .false.)
                call csv_write(pl_tree, site%runID, .false.)
                call csv_write(pl_tree, year, .false.)
                call csv_write(pl_tree, ip, .false.)
                call write_tree_csv(site%plots(ip)%trees(it), pl_tree)
            enddo
        enddo

    end subroutine write_tree_data

    !:.........................................................................:

end module Output
