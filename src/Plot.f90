module Plot

!*******************************************************************************
  !
  ! The Plot module contains the definition of the PlotData type, which holds
  ! attributes of an individual plot, as well as procedures directly operating
  ! on plot-level variables.
  !
!*******************************************************************************

use Constants
use Input
use Species
use Tree
use Soil

implicit none

    ! Define PlotData type
    type PlotData
        type(TreeData),    dimension(:),   allocatable :: trees        ! Array of tree objects
        type(TreeData),    dimension(:),   allocatable :: deadtrees    ! Array of dead tree objects
        type(SoilData)                                 :: soil         ! Soil object
        real,              dimension(:,:), allocatable :: cells        ! 2D array of plant locations
        real,              dimension(:),   allocatable :: avail_spec   ! Species "available" for regeneration
        real,              dimension(:),   allocatable :: seedbank     ! Species-specific seedbank (seeds/m2)
        real,              dimension(:),   allocatable :: seedling     ! Species-specific seedling bank (seedlings/m2)
        real,              dimension(:),   allocatable :: con_light    ! Conifer light environment (0-1)
        real,              dimension(:),   allocatable :: dec_light    ! Deciduous light environment (0-1)
        real,              dimension(:),   allocatable :: fc_nutr      ! Impact of nutrients on species-specific growth (0-1)
        real,              dimension(:),   allocatable :: fc_gdd       ! Impact of growing degree-days on species-specific growth (0-1)
        real,              dimension(:),   allocatable :: fc_drought   ! Impact of drought on species-specific growth (0-1)
        real,              dimension(:),   allocatable :: fc_perm      ! Impact of permafrost on species-specific growth (0-1)
        real,              dimension(:),   allocatable :: fc_fire      ! Impact of fire on species-specific regeneration
        integer,           dimension(:),   allocatable :: mature       ! "Mature" species - that can reproduce
        real,              dimension(M_TYPES)          :: d_type       ! Biomass killed by different types of mortality (tC/ha)
        real                                           :: wind_cat     ! Windthrow intensity category
        real                                           :: act_evap_day ! Actual evapotranspiration (cm)
        real                                           :: wilt_days    ! Proportion of growing season below wilting point
        real                                           :: flood_days   ! Proportion of growing season with flooded conditions
        real                                           :: dry_days     ! Proportion of growing season with drought conditions
        real                                           :: saw0_ByFC    ! A-layer moisture scaled by field capacity
        real                                           :: aow0_ByMin   ! Organic layer moisture scaled by wilting point
        real                                           :: saw0_BySAT   ! A-layer moisture scaled by saturation capacity
        real                                           :: amlt         ! Seasonal maximum depth of thaw (m)
        real                                           :: cla          ! Cumulative leaf area on the forest floor (m2)
        real                                           :: NPP          ! Net primary production (tC/ha)
        integer                                        :: numtrees     ! Number of live trees on plot
        integer                                        :: num_dead     ! Number of dead trees on plot
        integer                                        :: fire         ! Just had a fire (1) or no (0)
        integer                                        :: wind         ! Just had a windthrow event (1) or no (0)
        integer                                        :: windCount    ! Number of years since last windthrow event
        integer                                        :: stand_age    ! Stand age of site (years)
    end type PlotData

contains

    !:.........................................................................:

    subroutine initialize_plot(self, numspecies, a_sat, a_fc, a_pwp, o_sat,   &
            o_fc, o_pwp, o_bd, a_bd, itxt, hum_input)
        !
        !  Initializes a plot at a site with starting/input values.
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: constants
        real, parameter :: INIT_LAIW = 0.15  ! Initial canopy moisture (cm)
        real, parameter :: INIT_BIOWD = 42.5 ! Initial litter amount (t/ha)
        real, parameter :: BIOFF = 0.01617   ! Initial humus (kg/ha)
        real, parameter :: HCON = 0.658      ! Humus consumption amount (proportion)
        real, parameter :: PCTN = 0.008      ! Initial N % in humus

        ! Data dictionary: calling arguments
        class(PlotData), intent(inout) :: self       ! Plot object
        integer,         intent(in)    :: numspecies ! Number of species on plot
        real,            intent(in)    :: a_sat      ! A-layer saturation capacity (volumetric)
        real,            intent(in)    :: a_fc       ! A-layer field capacity (volumetric)
        real,            intent(in)    :: a_pwp      ! A-layer permanent wilting point (volumetric)
        real,            intent(in)    :: o_sat      ! Organic layer saturation capacity (volumetric)
        real,            intent(in)    :: o_fc       ! Organic layer field capacity (volumetric)
        real,            intent(in)    :: o_pwp      ! Organic layer permanent wilting point (volumetric)
        real,            intent(in)    :: o_bd       ! Organic layer bulk density (kg/m3)
        real,            intent(in)    :: a_bd       ! A-layer layer bulk density (kg/m3)
        real,            intent(in)    :: hum_input  ! Initial humus amount (t/ha)
        integer,         intent(in)    :: itxt       ! Soil texture (0: very coarse; 1: coarse; 2: fine)

        ! Data dictionary: local variables
        real            :: duff ! Duff content (kg)
        integer         :: n, i ! Looping indices

        if (maxcells .ne. 0) then
            ! Allocate tree and cells arrays
            allocate(self%trees(maxcells*maxcells))
            allocate(self%deadtrees(maxcells*maxcells))
            allocate(self%cells(maxcells, maxcells))
        else
            stop "Must allow at least a few trees"
        endif

        if (maxheight .ne. 0) then
            ! Allocate light arrays
            allocate(self%con_light(maxheight))
            allocate(self%dec_light(maxheight))
        else
            stop "Must have a nonzero maximum height"
        endif

        ! Allocate species-specific arrays
        allocate(self%seedling(numspecies))
        allocate(self%seedbank(numspecies))
        allocate(self%avail_spec(numspecies))
        allocate(self%fc_nutr(numspecies))
        allocate(self%fc_perm(numspecies))
        allocate(self%fc_gdd(numspecies))
        allocate(self%fc_drought(numspecies))
        allocate(self%fc_fire(numspecies))
        allocate(self%mature(numspecies))

        ! Set to starting values
        self%cells = .false.
        self%numtrees = 0
        self%num_dead = 0
        self%fire = 0
        self%wind = 0
        self%windCount = 0
        self%stand_age = 0
        self%wind_cat = 0.0
        self%cla = 0.0
        self%act_evap_day = 0.0
        self%flood_days = 0.0
        self%wilt_days = 0.0
        self%dry_days = 0.0
        self%aow0_ByMin = 0.0
        self%saw0_ByFC = 0.0
        self%saw0_BySAT = 0.0
        self%amlt = 0.0
        self%seedbank = 0.0
        self%seedling = 0.0
        self%mature = 0
        self%avail_spec = 0.0
        self%dec_light = 0.0
        self%con_light =0.0
        self%d_type = 0.0
        self%fc_nutr = 1.0
        self%fc_gdd = 1.0
        self%fc_drought = 1.0
        self%fc_fire = 1.0
        self%fc_perm = 1.0
        self%soil%lai_w0 = INIT_LAIW
        self%soil%twig_fuel = 0.0
        self%soil%smbl_fuel = 0.0
        self%soil%lrbl_fuel = 0.0
        self%soil%avail_fuel = 0.0

        ! Initialize soil properties

        ! Snowpack initialized to 0.0
        self%soil%snowpack = 0.0
        self%soil%swe = 0.0

        ! Litter initialization
        do i = 1, LIT_LEVS
            self%soil%litter(i) = 0.0
        end do
        self%soil%cohorts = 0.0

        ! Initializing from bare ground - assume secondary succession conditions
        ! wood litter layer (fresh wood)
        self%soil%litter(15) = INIT_BIOWD ! small boles (t/ha)
        self%soil%litter(14) = INIT_BIOWD ! large boles (t/ha)

        ! Humus layer
        self%soil%cohorts(1, 1) = hum_input
        self%soil%cohorts(1, 2) = self%soil%cohorts(1, 1)*PCTN ! tN/ha
        self%soil%cohorts(1, 5) = 21.0
        self%soil%ncohort = 1

        ! N from recent fires
        self%soil%fan = 0.0

        ! Calculate organic layer depth
        duff = self%soil%cohorts(1, 1)/HEC_TO_M2*plotsize*T_TO_KG ! kg/plot
        self%soil%O_depth = (1.0/plotsize)*(duff/BULK_DUFF)

        ! Initialize A-layer, active, and moss-layer depths
        self%soil%A_depth = 1.0
        self%soil%M_depth = 0.0
        self%soil%moss_biom = 0.0
        self%soil%active = self%soil%A_depth

        ! Initialize values from site file
        self%soil%O_bulk_dens = o_bd
        self%soil%A_bulk_dens = a_bd
        self%soil%sat(1) = o_sat
        self%soil%fc(1) = o_fc
        self%soil%pwp(1) = a_pwp
        self%soil%pwp(2) = o_pwp
        self%soil%sat(2) = a_sat
        self%soil%fc(2) = a_fc
        self%soil%itxt = itxt

    end subroutine initialize_plot

    !:.........................................................................:

    subroutine delete_plot(self)
        !
        !  Frees the memory associated with a plot
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        class(PlotData), intent(inout) :: self ! Plot object

        ! Data dictionary: local variables
        integer :: it ! Looping index


        ! Deallocate arrays
        do it = 1, self%numtrees
            call delete_tree(self%trees(it))
        end do
        do it = 1, self%num_dead
            call delete_tree(self%deadtrees(it))
        end do
        if (allocated(self%trees)) deallocate(self%trees)
        if (allocated(self%deadtrees)) deallocate(self%deadtrees)
        if (allocated(self%cells)) deallocate(self%cells)
        if (allocated(self%avail_spec)) deallocate(self%avail_spec)
        if (allocated(self%seedbank)) deallocate(self%seedbank)
        if (allocated(self%seedling)) deallocate(self%seedling)
        if (allocated(self%con_light)) deallocate(self%con_light)
        if (allocated(self%dec_light)) deallocate(self%dec_light)
        if (allocated(self%fc_nutr)) deallocate(self%fc_nutr)
        if (allocated(self%fc_drought)) deallocate(self%fc_drought)
        if (allocated(self%fc_perm)) deallocate(self%fc_perm)
        if (allocated(self%fc_gdd)) deallocate(self%fc_gdd)
        if (allocated(self%fc_fire)) deallocate(self%fc_fire)
        if (allocated(self%mature)) deallocate(self%mature)

    end subroutine delete_plot

    !:.........................................................................:

    subroutine tree_dm_cats(self, genera, field, diam_categories)
        !
        !  Bins the number of trees on a plot into DBH size classes, for
        !    writing to output files
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        class(PlotData),                         intent(in)  :: self            ! Plot object
        character(len = MAX_NLEN), dimension(:), intent(in)  :: genera          ! Names of genera/species
        character(len = 8),                      intent(in)  :: field           ! 'genus' or 'species'
        integer, dimension(size(genera), NHC),   intent(out) :: diam_categories ! Array of binned tree density

        ! Data dictionary: local variables
        character(len = MAX_NLEN) :: comp
        integer, dimension(NHC)   :: tree_diams
        integer                   :: numitems   ! Number of species/genera
        integer                   :: ig, n      ! Looping indices

        ! Get number of species/genera
        numitems = size(genera)

        ! Initialize accumulation of trees
        diam_categories = 0

        ! Loop through genera and trees
        do ig = 1, numitems
            do n = 1, self%numtrees

                ! Get correct identifier - genus or species ID
                if (field .eq. 'genus') then
                    comp = self%trees(n)%spec_ptr%genus_name
                else if (field .eq. 'species') then
                    comp = self%trees(n)%spec_ptr%unique_id
                endif

                ! If we are on the right one - get the DBH categories
                if (comp == genera(ig)) then
                    call get_diam_category(self%trees(n), tree_diams)
                    diam_categories(ig,:) = diam_categories(ig,:) + tree_diams
                endif
            enddo
        enddo

    end subroutine tree_dm_cats

    !:.........................................................................:

    subroutine sum_over_sg(self, genera, field, basal_area, leaf_bm, biomC,    &
        biomN, max_ht, max_diam, mean_dbh, mean_year, biom_cats, envresp_mn,   &
        basal_lg, basal_sm, dbh_lg, dbh_sm, n_lg, n_sm, biomC_lg, biomC_sm)
        !
        !  Aggregates live species- or genus-level data for plot
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !
        class(PlotData),                         intent(in)  :: self       ! Plot object
        character(len = MAX_NLEN), dimension(:), intent(in)  :: genera     ! Species or genus names
        character(len = 8),                      intent(in)  :: field      ! 'species' or 'genus'
        real, dimension(size(genera)),           intent(out) :: basal_area ! Basal area sum (cm2)
        real, dimension(size(genera)),           intent(out) :: leaf_bm    ! Leaf biomass sum (tC)
        real, dimension(size(genera)),           intent(out) :: mean_dbh   ! Mean DBH (cm)
        real, dimension(size(genera)),           intent(out) :: mean_year  ! Average tree age (years)
        real, dimension(size(genera)),           intent(out) :: basal_lg   ! Basal area sum > 9 cm DBH (cm2)
        real, dimension(size(genera)),           intent(out) :: basal_sm   ! Basal area sum < 9 cm DBH (cm2)
        real, dimension(size(genera)),           intent(out) :: dbh_lg     ! Mean DBH > 9 cm DBH (cm)
        real, dimension(size(genera)),           intent(out) :: dbh_sm     ! Mean DBH < 9 cm DBH (cm)
        real, dimension(size(genera)),           intent(out) :: n_lg       ! Total stems > 9 cm DBH (#)
        real, dimension(size(genera)),           intent(out) :: n_sm       ! Total stems < 9 cm DBH (#)
        real, dimension(size(genera)),           intent(out) :: biomC      ! Biomass sum (tC)
        real, dimension(size(genera)),           intent(out) :: biomN      ! N content sum (tN)
        real, dimension(size(genera)),           intent(out) :: biomC_lg   ! Biomass sum > 9 cm DBH (tC)
        real, dimension(size(genera)),           intent(out) :: biomC_sm   ! Biomass sum < 9 cm DBH (tC)
        real, dimension(size(genera)),           intent(out) :: max_ht     ! Maximum height (m)
        real, dimension(size(genera)),           intent(out) :: max_diam   ! Maximum diameter (cm)
        real, dimension(size(genera), NHC),      intent(out) :: biom_cats  ! Biomass sum binned by DBH (tC)
        real, dimension(size(genera), FC_NUM),   intent(out) :: envresp_mn ! Mean growth response (0-z1)

        ! Data dictionary: local variables
        real, dimension(size(genera), FC_NUM) :: envresp      ! Sum of growth responses (0-1)
        real, dimension(size(genera), FC_NUM) :: n_er         ! Number of trees check for each growth response
        real, dimension(size(genera))         :: n            ! Total stems (#)
        real, dimension(size(genera))         :: sum_dbh      ! Sum of diameters (cm)
        real, dimension(size(genera))         :: sumdbh_lg    ! Sum of diameters > 9.0 cm
        real, dimension(size(genera))         :: sumdbh_sm    ! Sum of diameters < 9.0 cm
        real, dimension(size(genera))         :: sum_year     ! Sum of tree ages (year)
        character(len=MAX_NLEN)               :: comp         ! Tree identifier (species or genus)
        real                                  :: lf_biom      ! Local tree leaf biomass (tC)
        real                                  :: l_cn         ! Tree C:N ratio
        real                                  :: tot_tree_bc  ! Total tree biomass (tC)
        real                                  :: tot_tree_bn  ! Total tree N content (tN)
        real                                  :: dm           ! Local tree diameter (cm)
        integer                               :: numitems     ! Number of species/genera
        integer                               :: it, is, k    ! Looping indices

        ! Initialize accumulators
        basal_area = 0.0
        biomC = 0.0
        biomN = 0.0
        leaf_bm = 0.0
        max_ht = 0.0
        max_diam = 0.0
        n = 0.0
        n_er = 0.0
        sum_dbh = 0.0
        sum_year = 0.0
        biom_cats = 0.0
        envresp = 0.0
        basal_lg = 0.0
        basal_sm = 0.0
        biomC_lg = 0.0
        biomC_sm = 0.0
        sumdbh_lg = 0.0
        sumdbh_sm = 0.0
        n_lg = 0.0
        n_sm = 0.0
        dbh_lg = 0.0
        dbh_sm = 0.0

        ! Get size of species/genera
        numitems = size(genera)

        do is = 1, numitems
            do it = 1, self%numtrees

                ! Get identifier (genus or unique ID)
                if (field .eq. 'genus') then
                    comp = self%trees(it)%spec_ptr%genus_name
                else if (field .eq. 'species') then
                    comp = self%trees(it)%spec_ptr%unique_id
                endif

                if (comp == genera(is)) then
                    ! If we have the correct species/genus, sum up values

                    ! Grab diameter
                    dm = self%trees(it)%diam_bht

                    if (self%trees(it)%diam_bht .ge. 9.0) then

                        ! Calculate and sum basal area (m2)
                        basal_lg(is) = basal_lg(is) +                          &
                            CM_TO_M*CM_TO_M*0.25*pi*dm**2

                        ! Calculate and sum biomass (tC)
                        lf_biom = self%trees(it)%leaf_bm
                        tot_tree_bc = self%trees(it)%stemC +                   &
                            self%trees(it)%branchC + lf_biom
                        biomC_lg(is) = biomC_lg(is) + tot_tree_bc

                        ! Sum up number and diameter
                        sumdbh_lg(is) = sumdbh_lg(is) + dm
                        n_lg(is) = n_lg(is) + 1.0

                    else if (self%trees(it)%diam_bht .lt. 9.0) then
                        if (self%trees(it)%forska_ht .ge. 1.3) then

                            ! Calculate and sum basal area (m2)
                            basal_sm(is) = basal_sm(is) +                      &
                                CM_TO_M*CM_TO_M*0.25*pi*dm**2

                            ! Calculate and sum biomass (tC)
                            lf_biom = self%trees(it)%leaf_bm
                            tot_tree_bc = self%trees(it)%branchC +             &
                                self%trees(it)%stemC + lf_biom
                            biomC_sm(is) = biomC_sm(is) + tot_tree_bc

                            ! Sum up number and diameter
                            sumdbh_sm(is) = sumdbh_sm(is) + dm
                            n_sm(is) = n_sm(is) + 1.0
                        end if
                    end if

                    ! Number of trees in that genera
                    n(is) = n(is) + 1.0

                    ! Get C:N ratio
                    if (self%trees(it)%conifer) then
                        l_cn = CON_LEAF_C_N
                    else
                        l_cn = DEC_LEAF_C_N
                    endif

                    ! Calculate and sum up basal area (m2)
                    basal_area(is) = basal_area(is) +                          &
                        CM_TO_M*CM_TO_M*0.25*pi*dm**2

                    ! Calculate and sum biomass (tC) and N content
                    lf_biom = self%trees(it)%leaf_bm
                    leaf_bm(is) = leaf_bm(is) + lf_biom

                    tot_tree_bc = self%trees(it)%branchC +                     &
                        self%trees(it)%stemC + lf_biom
                    tot_tree_bn = self%trees(it)%biomN + lf_biom/l_cn

                    biomC(is) = biomC(is) + tot_tree_bc
                    biomN(is) = biomN(is) + tot_tree_bn

                    ! Get biomass in DBH bins
                    do k = 1, NHC - 1
                        if (dm > DBC_MIN(k) .and. dm <= DBC_MAX(k)) then
                            biom_cats(is, k) = biom_cats(is, k) + tot_tree_bc
                        end if
                    end do
                    if (dm > DBC_MIN(NHC)) then
                        biom_cats(is, NHC) = biom_cats(is, NHC) + tot_tree_bc
                    end if

                    ! Get maximum height and DBH
                    max_ht(is) = max(max_ht(is), self%trees(it)%forska_ht)
                    max_diam(is) = max(max_diam(is), self%trees(it)%diam_bht)

                    ! Sum up DBH and age
                    sum_dbh(is) = sum_dbh(is) + self%trees(it)%diam_bht
                    sum_year(is) = sum_year(is) + float(self%trees(it)%tree_age)

                    ! Sum up the growth responses for each factor
                    do k = 1, FC_NUM
                        if (self%trees(it)%env_resp(k) >= epsilon(1.0)) then
                            n_er(is, k) = n_er(is, k) + 1
                            envresp(is, k) = envresp(is, k) +                  &
                            self%trees(it)%env_resp(k)
                        end if
                    end do

                endif
            enddo
        enddo

        ! Get average of DBH, age, and environmental responses
        do is = 1, numitems
            if (n(is) .eq. 0.0) then
                mean_dbh(is) = 0.0
                mean_year(is) = 0.0
                envresp_mn(is, :) = RNVALID
            else
                mean_dbh(is) = sum_dbh(is)/n(is)
                mean_year(is) = sum_year(is)/n(is)
                do k = 1, FC_NUM
                    if (n_er(is, k) > 0.0) then
                        envresp_mn(is, k) = envresp(is, k)/n_er(is, k)
                    else
                        envresp_mn(is, k) = RNVALID
                    end if
                end do
            endif

            if (n_lg(is) .eq. 0.0) then
                dbh_lg(is) = 0.0
            else
                dbh_lg(is) = sumdbh_lg(is)/n_lg(is)
            end if
            if (n_sm(is) .eq. 0.0) then
                dbh_sm(is) = 0.0
            else
                dbh_sm(is) = sumdbh_sm(is)/n_sm(is)
            end if
        enddo

    end subroutine sum_over_sg

    !:.........................................................................:

    subroutine sum_over_deadsg(self, genera, field, biomC, mean_dbh,           &
        death_markers)
        !
        !  Aggregates dead species- or genus-level data for plot
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/15     A. C. Foster          Original Code
        !

        ! Data dictionary: calling arguments
        class(PlotData),                         intent(in)  :: self          ! Plot object
        character(len = MAX_NLEN), dimension(:), intent(in)  :: genera        ! Species or genus names
        character(len = 8),                      intent(in)  :: field         ! 'species' or 'genus'
        real, dimension(size(genera)),           intent(out) :: mean_dbh      ! Mean DBH (cm)
        real, dimension(size(genera)),           intent(out) :: biomC         ! Biomass sum (tC)
        real, dimension(size(genera), M_TYPES),  intent(out) :: death_markers ! Stems lost to each mortality type (#)

        ! Data dictionary: local variables
        integer, dimension(size(genera)) :: n           ! Total stems (#)
        real,    dimension(size(genera)) :: sum_dbh     ! Sum of DBHs (cm)
        character(len=MAX_NLEN)          :: comp        ! Tree identifier (species/genus)
        real                             :: lf_biom     ! Tree leaf biomass (tC)
        real                             :: l_cn        ! Tree C:N ratio
        real                             :: tot_tree_bc ! Total tree biomass (tC)
        real                             :: dm          ! Tree diameter (cm)
        integer                          :: numitems    ! Total number of species/genera
        integer                          :: it, is      ! Looping indices
        integer                          :: ks          ! Killing factor for tree

        ! Get total number of species/genera
        numitems = size(genera)

        ! Initialize accumulators
        biomC = 0.0
        n = 0
        sum_dbh = 0.0
        death_markers = 0.0

        do is = 1, numitems
            do it = 1, self%num_dead

                ! Get identifier (genus or unique ID)
                if (field .eq. 'genus') then
                    comp = self%deadtrees(it)%spec_ptr%genus_name
                else if (field .eq. 'species') then
                    comp = self%deadtrees(it)%spec_ptr%unique_id
                endif

                ! If we have the correct species/genus
                if (comp == genera(is)) then

                    ! Get number of trees in that genera
                    n(is) = n(is) + 1

                    ! Get C:N ratio
                    if (self%deadtrees(it)%conifer) then
                        l_cn = CON_LEAF_C_N
                    else
                        l_cn = DEC_LEAF_C_N
                    endif

                    ! What was the killing factor
                    ks = self%deadtrees(it)%stressor

                    ! Get diameter
                    dm = self%deadtrees(it)%diam_bht

                    ! Get and sum up biomass C
                    lf_biom = self%deadtrees(it)%leaf_bm
                    tot_tree_bc = self%deadtrees(it)%branchC +                 &
                        self%deadtrees(it)%stemC + lf_biom
                    biomC(is) = biomC(is) + tot_tree_bc

                    ! Bin biomass into the limiting factors/disturbances
                    death_markers(is, ks) = death_markers(is, ks) + tot_tree_bc

                    ! Sum up DBH
                    sum_dbh(is) = sum_dbh(is) + self%deadtrees(it)%diam_bht

                endif
            end do
        end do

        ! Get average DBH
        do is = 1, numitems
            if (n(is) .eq. 0) then
                mean_dbh(is) = 0.0
            else
                mean_dbh(is) = sum_dbh(is)/n(is)
            endif
        enddo

    end subroutine sum_over_deadsg

    !:.........................................................................:

end module Plot
