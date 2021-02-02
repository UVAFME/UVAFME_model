module Tree
!*******************************************************************************
  !
  ! The Tree module contains attributes and procedures relevant to tree growth
  ! and mortality
  !
!*******************************************************************************

    use Constants
    use Species
    use Random

    implicit none

    ! Define tree data type
    type                            :: TreeData
        class(SpeciesData), pointer :: spec_ptr => null() ! Pointer to species data
        real, dimension(FC_NUM)     :: env_resp           ! Species growth response to stressors (0-1)
        real                        :: diam_opt           ! Optimum diameter increment growth (cm)
        real                        :: diam_bht           ! Diameter at breast height (cm)
        real                        :: diam_canht         ! Diameter at clear branch bole height (cm)
        real                        :: canopy_ht          ! Clear branch bole height (m)
        real                        :: forska_ht          ! Tree height (m)
        real                        :: leaf_bm            ! Leaf biomass (tC)
        real                        :: biomC              ! Aboveground woody biomass (tC)
        real                        :: stemC              ! Stem biomass (tC)
        real                        :: branchC            ! Branch biomass (tC)
        real                        :: rootC              ! Root biomass (tC)
        real                        :: biomN              ! N content (tC)
        real                        :: CK                 ! Proportion of crown scorched by fire (%)
        real                        :: row                ! Row location
        real                        :: col                ! Column location
        integer                     :: mort_count         ! How many years tree has experiences low growth
        integer                     :: tree_age           ! Tree age (years)
        integer                     :: stressor           ! Which factor is most stressful
        integer                     :: conifer            ! 1: conifer; 2: deciduous
        integer                     :: species_index      ! Species index - points to SpeciesData array location
        logical                     :: mort_marker        ! Is tree marked for death?
    end type TreeData

contains


    !:.........................................................................:

    subroutine initialize_tree(self, tree_species, is)
        !
        !  Initializes a tree object with initial conditions
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        class(TreeData),           intent(inout) :: self         ! Tree object
        type(SpeciesData), target, intent(in)    :: tree_species ! Species object
        integer,                   intent(in)    :: is           ! Integer for species array

        ! Constructor
        self%env_resp = RNVALID
        self%diam_opt = 0.0
        self%diam_bht = 0.0
        self%diam_canht = 0.0
        self%canopy_ht = STD_HT
        self%forska_ht = STD_HT
        self%leaf_bm = 0.0
        self%biomC = 0.0
        self%stemC = 0.0
        self%branchC = 0.0
        self%rootC = 0.0
        self%biomN = 0.0
        self%CK = 0.0
        self%row = 0
        self%col = 0
        self%mort_count = 0
        self%tree_age = 0
        self%stressor = 0
        self%mort_marker = .false.
        self%species_index = is

        ! Attach species data and species index
        if (associated(self%spec_ptr)) nullify(self%spec_ptr)
        allocate(self%spec_ptr)
        self%spec_ptr => tree_species

        ! Grab this from species data for easier use
        self%conifer = self%spec_ptr%conifer

    end subroutine initialize_tree

    !:.........................................................................:

    subroutine copy_tree(self, tree)
        !
        !  Copies attributes from one tree object to another
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        class(TreeData), intent(inout) :: self
        class(TreeData), intent(in)    :: tree

        ! Copy attributes
        self%env_resp       = tree%env_resp
        self%spec_ptr       => tree%spec_ptr
        self%species_index  = tree%species_index
        self%diam_opt       = tree%diam_opt
        self%diam_bht       = tree%diam_bht
        self%diam_canht     = tree%diam_canht
        self%canopy_ht      = tree%canopy_ht
        self%forska_ht      = tree%forska_ht
        self%leaf_bm        = tree%leaf_bm
        self%biomC          = tree%biomC
        self%stemC          = tree%stemC
        self%branchC        = tree%branchC
        self%rootC          = tree%rootC
        self%biomN          = tree%biomN
        self%CK             = tree%CK
        self%row            = tree%row
        self%col            = tree%col
        self%mort_count     = tree%mort_count
        self%tree_age       = tree%tree_age
        self%stressor       = tree%stressor
        self%mort_marker    = tree%mort_marker
        self%conifer        = tree%conifer

    end subroutine copy_tree


    !:.........................................................................:

    subroutine env_stress(self, shade, fc_gdd, fc_drought, fc_perm, envstress, &
        fc_nutr)
        !
        !  Calculates the environmental stress on a tree using Liebig's Law of
        !   the Minimum, and sets what the most stressful factor is
        !
        !    stressor:  which factor was most stressful
        !                1: temperature
        !                2: drought
        !                3: shade
        !                4: permafrost
        !                5: nutrients
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/16     A. C. Foster        Original Code
        !    01/27/21     A. C. Foster        Updated to read in effects rather
        !                                       than use attached species
        !                                       attributes

        ! Data dictionary: calling arguments
        class(TreeData), intent(inout)        :: self       ! Tree object
        real,            intent(in)           :: shade      ! Growth response to shading (0-1)
        real,            intent(in)           :: fc_gdd     ! Growth response to growing degree-days (0-1)
        real,            intent(in)           :: fc_drought ! Growth response to drought (0-1)
        real,            intent(in)           :: fc_perm    ! Growth response to permafrost (0-1)
        real,            intent(out)          :: envstress  ! Growth response to stress
        real,            intent(in), optional :: fc_nutr

        ! Data dictionary: local variables
        real    :: minstress  ! Minimum stress from gdd, drought, shade, and (optionally) nutrients
        integer :: stressor   ! Which factor had the highest impact

        ! Get the minimum value
        if (present(fc_nutr)) then

            ! Calculate using nutrient stress
            minstress = min(fc_gdd, fc_drought, shade, fc_nutr)

            ! Find the most limiting factor
            if (minstress .eq. fc_gdd) then
                stressor = 1
            else if (minstress .eq. fc_drought) then
                stressor = 2
            else if (minstress .eq. shade) then
                stressor = 3
            else if (minstress .eq. fc_nutr) then
                stressor = 5
            end if

            if (minstress .gt. fc_perm) then
                stressor = 4
            endif

            ! Output is minimum of shade, gdd, nutrient, and drought times flood
            ! and permafrost
            envstress = minstress*fc_perm

        else
            ! Don't use nutrient stress
            ! This is so we can calculate the potential growth and then
            ! decrease it based on available N
            minstress = min(fc_gdd, fc_drought, shade)

            ! Find the most limiting factor
            if (minstress .eq. fc_gdd) then
                stressor = 1
            else if (minstress .eq. fc_drought) then
                stressor = 2
            else if (minstress .eq. shade) then
                stressor = 3
            end if

            if (minstress .gt. fc_perm) then
                stressor = 4
            endif

            ! Output is minimum of shade, gdd, and drought times flood and
            ! permafrost
            envstress = minstress*fc_perm

        end if

        ! Set the environmental stressors
        self%stressor = stressor
        self%env_resp(1) = fc_gdd
        self%env_resp(2) = fc_drought
        self%env_resp(3) = shade
        self%env_resp(4) = fc_perm
        if(present(fc_nutr)) self%env_resp(5) = fc_nutr

    end subroutine env_stress

    !:.........................................................................:

    subroutine stem_shape(tree)
        !
        !  Calculates diameter at clear branch bole height (cm)
        !   Notes from original code:
        !       basal diameter = dhshape(0.0, h, dbh)
        !       diameter @ ccb = dhshape(hc, h, sdbh)
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure
        !    10/10/20     A. C. Foster        Update for shrub forms
        !

        ! Data dictionary: calling arguments
        class(TreeData), intent(inout) :: tree ! Tree object

        ! Data dictionary: local variables
        real :: hcbb ! Clear branch bole height (m)
        real :: ht   ! Tree height (m)
        real :: dbh  ! Diameter at breast height (cm)
        real :: dcbb ! Diameter at clear branch bole height (cm)
        real :: beta ! Stem shape parameter

        hcbb = tree%canopy_ht
        ht = tree%forska_ht
        dbh = tree%diam_bht
        beta = tree%spec_ptr%beta

        if (ht .le. hcbb .or. ht .le. STD_HT) then
            ! Just set to dbh
            dcbb = dbh
        else
            ! Calculate
            dcbb = ((ht - hcbb)/(ht - STD_HT))**(1.0/beta)*dbh
        end if
        ! Set
        tree%diam_canht = dcbb


    end subroutine stem_shape

    !:.........................................................................:

    subroutine forska_height(tree)
        !
        !  Calculates total tree height (m)
        !  Equation adapted from FORSKA model Leemans & Prentice 1989
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure
        !    10/10/20     A. C. Foster        Update for shrub forms
        !

        ! Data dictionary: calling arguments
        class(TreeData), intent(inout) :: tree ! Tree object

        ! Data dictionary: local variables
        real :: dbh      ! Diameter at breast height (cm)
        real :: s        ! Growth parameter
        real :: hmax     ! Average maximum height (m)
        real :: delta_ht ! Temporary variable

        dbh = tree%diam_bht
        hmax = tree%spec_ptr%max_ht
        s = tree%spec_ptr%s

        ! Tree or tree-like shrub form
        delta_ht = hmax - STD_HT
        tree%forska_ht = STD_HT + delta_ht*(1.0 - exp(-(s*dbh/delta_ht)))

    end subroutine forska_height

    !:.........................................................................:

    subroutine biomass_c(tree)
        !
        !  Calculates aboveground woody biomass (tC)
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure
        !    10/10/20     A. C. Foster        Update for shrub forms
        !

        ! Data dictionary: calling arguments
        class(TreeData), intent(inout) :: tree ! Tree object

        ! Data dictionary: local varibles
        real            :: dbh       ! Diameter at breask height (cm)
        real            :: ht        ! Tree height (m)
        real            :: hcbb      ! Clear branch bole height (m)
        real            :: hrt       ! Rooting depth (m)
        real            :: bulk      ! Wood bulk density (t/m3)
        real            :: stembc    ! Stem biomass (tC)
        real            :: branchbc  ! Branch biomass (tC)
        real            :: abovegr_c ! Aboveground woody biomass (tC)
        real            :: root_c    ! Root biomass (tC)

        ht = tree%forska_ht
        hrt = tree%spec_ptr%rootdepth

        ! Calculate diameter at cbb
        call stem_shape(tree)

        ! Calculate stem biomass (tC)
        stembc   = stem_biomass_c(tree)

        ! Calculate branch biomass (tC)
        branchbc = branch_biomass_c(tree)

        ! Aboveground woody biomass is branch + stem
        abovegr_c = stembc + branchbc

        ! Calculate root biomass (tC)
        root_c    = stembc*hrt/ht + branchbc

        ! Set values
        tree%biomC = abovegr_c + root_c
        tree%stemC = stembc
        tree%branchC = branchbc
        tree%rootC = root_c


    end subroutine biomass_c

    !:.........................................................................:

    subroutine biomass_n(tree)
        !
        !  Calculates stem and branch N biomass (tN)
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure
        !

        ! Data dictionary: calling arguments
        class(TreeData), intent(inout) :: tree ! Tree object

        ! Use stem C:N ratio
        tree%biomN = tree%biomC/STEM_C_N

    end subroutine biomass_n

    !:.........................................................................:

    subroutine leaf_biomass_c(tree)
        !
        !  Calculates leaf biomass (tC)
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure

        ! Data dictionary: calling arguments
        class(TreeData), intent(inout) :: tree ! Tree object

        tree%leaf_bm = leaf_area(tree)*tree%spec_ptr%leafarea_c*2.0

    end subroutine leaf_biomass_c

    !:.........................................................................:

    real function leaf_area(tree)
        !
        !  Calculates leaf area of a tree (m2)
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure

        ! Data dictionary: calling arguments
        class(TreeData), intent(in) :: tree ! Tree object

        ! Data dictionary: local variables
        real :: dcbb    ! Diameter at clear branch bole height
        real :: leafd_a ! Scalar parameter for DCBB:leaf area relationship

        dcbb = tree%diam_canht
        leafd_a = tree%spec_ptr%leafdiam_a

        leaf_area = (dcbb*dcbb)*leafd_a

    end function leaf_area

    !:.........................................................................:

    subroutine age_survival(tree, a_survive)
        !
        !  Random death check. Trees can reach maximum age as a probability of
        !  1% ('often seen'), 0.1% ('some survive'), or 0.01% ('few seen'), as
        !  set by input age_tol parameter (1-3)
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure

        ! Data dictionary: constants
        ! Age-tolerance related parameters
        ! 1% ---> 4.605  0.1% ---> 6.908  0.01% ---> 11.5129
        real, dimension(3), parameter :: CHECK = [4.605, 6.908, 11.51]

        ! Data dictionary: calling arguments
        class(TreeData), intent(in)  :: tree      ! Tree object
        logical,         intent(out) :: a_survive ! Does tree survive?

        ! Data dictionary: local variables
        integer :: k        ! Probability of surviving to max age (1-3; 3 = least likely)
        real    :: agemax   ! Average maximum age (years)
        real    :: rand_val ! Random number (uniform)

        k = tree%spec_ptr%age_tol
        agemax = tree%spec_ptr%max_age

        ! Get random uniform between 0 and 1
        rand_val = urand()

        ! Check against parameter
        if (rand_val .lt. CHECK(k)/agemax) then
            a_survive = .false.
        else
            a_survive = .true.
        end if

    end subroutine age_survival

    !:.........................................................................:

    subroutine growth_survival(tree, g_survive)
        !
        !  Death check by low growth stress.
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure

        ! Data dictionary: constants
        ! Values to check against (depends on stress_tol parameter)
        real, dimension(5), parameter :: CHECK = [0.21, 0.24, 0.27, 0.30, 0.33]

        ! Data dictionary: calling arguments
        class(TreeData), intent(in) :: tree      ! Tree object
        logical, intent(out)        :: g_survive ! Does tree survive?

        ! Data dictionary: local variables
        integer :: k        ! Relative stress tolerance (1-5; 5 = least tolerant)
        real    :: rand_val ! Random value (uniform)

        k = tree%spec_ptr%stress_tol

        ! Get random uniform between 0 and 1.0
        rand_val = urand()

        ! Check against parameter
        if (tree%mort_marker .and. rand_val .lt. CHECK(k)) then
            g_survive = .false.
        else
            g_survive = .true.
        end if

    end subroutine growth_survival

    subroutine fire_survival(tree, av_fuel, f_survive)
        !
        !  Calculates survival of a tree from a fire event
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    10/10/18     A. C. Foster         Original Code

        ! Data dictionary: constants
        real, parameter :: CK1 = 0.21111
        real, parameter :: CK2 = -0.00445

        ! Data dictionary: calling arguments
        class(TreeData), intent(inout) :: tree      ! Tree object
        real,            intent(in)    :: av_fuel   ! Available fuel for burning (t/ha)
        logical,         intent(out)   :: f_survive ! Did tree survive?

        ! Data dictionary: local variables
        real :: pfire    ! Probability of fire mortality
        real :: dbh_eff  ! Effective DBH (cm)
        real :: rand_val ! Random uniform

        ! Get effective DBH (cm)
        if (tree%diam_bht >= 40.0) then
            dbh_eff = 40.0
        else
            dbh_eff = tree%diam_bht
        end if

        ! Calculate crown scorch (%)
        tree%CK = min(100.0, 100.0*(CK1 + CK2*dbh_eff)*av_fuel)
        if (tree%CK <= epsilon(1.0)) tree%CK = 0.0

        ! Fire probability
        pfire = 1/(1 + exp(-1.466 +                                            &
                        1.91*(tree%spec_ptr%bark_thick*tree%diam_bht) -        &
                        0.1775*((tree%spec_ptr%bark_thick*tree%diam_bht)**2) - &
                        0.000541*tree%CK**2))

        ! Check for fire mortality
        rand_val = urand()
        if (rand_val .lt. Pfire) then
            f_survive = .false.
        else
            f_survive = .true.
        endif

    end subroutine fire_survival

    !:.........................................................................:

    subroutine wind_survival(tree, w_survive)
        !
        !  Calculates survival from a windthrow event
        !  Adapted from Rich et al. 2007 Journal of Ecology 95:1261-1273
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    08/01/15     A. C. Foster        Original Code

        ! Data dictionary: calling arguments
        class(TreeData), intent(in)  :: tree      ! Tree object
        logical,         intent(out) :: w_survive ! Does tree survive?

        ! Data dictionary: local variables
        real :: logit    ! Temporary variable
        real :: pwind    ! Probability of mortality from windthrow
        real :: rand_val ! Random number (uniform)

        ! Probability of mortality from windthrow
        ! Equation from Rich et al. (2007)
        logit = 0.75*log(tree%diam_bht)
        pwind = 1.0/(1.0+exp(-logit))

        ! Check for survival
        rand_val=urand()
        if (rand_val .lt. pwind) then
            w_survive = .false.
        else
            w_survive = .true.
        endif


    end subroutine wind_survival

    !:.........................................................................:

    subroutine max_growth(tree)
        !
        !  Calculates maximum DBH increment growth given "optimal" environmental
        !   conditions
        !  Based on equation from Botkin et al. 1972 Journal of Ecology 60(3):849
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure
        !    10/10/20     A. C. Foster        Update for shrub forms

        ! Data dictionary: calling arguments
        class(TreeData), intent(inout) :: tree ! Tree object

        real :: dbh  ! Diameter at breast height (cm)
        real :: ht   ! Tree height (m)
        real :: dmax ! Average maximum dbh (cm)
        real :: hmax ! Average maximum height (m)
        real :: g    ! Growth parameter
        real :: s    ! Growth parameter
        real :: temp ! Temporary variable

        dbh = tree%diam_bht
        ht = tree%forska_ht
        dmax = tree%spec_ptr%max_diam
        hmax = tree%spec_ptr%max_ht
        g = tree%spec_ptr%g
        s = tree%spec_ptr%s

        temp = s*exp(-s*dbh/(hmax - STD_HT))*dbh
        tree%diam_opt = g*dbh*(1.0 - dbh*ht/dmax/hmax)/(2.0*ht + temp)

    end subroutine max_growth

    !:.........................................................................:

    real function stem_biomass_c(tree)
        !
        !  Calculates stem biomass (tC)
        !   Notes from previous model version:
        !       bstemC = bc(bd, 1.0, bulk_density)
        !       Here bd was "basal diameter" - diameter at clear branch bole
        !       height, but for a clear branch bole height of 0.0, so it was
        !       equal to dbh
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure
        !    01/28/21     A. C. Foster        Update to use dbh of tree, without
        !                                       setting cbb to 0.0
        !

        ! Data dictionary: calling arguments
        class(TreeData), intent(in) :: tree ! Tree object

        real :: dbh          ! Diameter at breast height (cm)
        real :: ht           ! Tree height (m)
        real :: bulk_density ! Wood bulk density (t/m3)
        real :: beta         ! Stem shape parameter
        real :: temp         ! Temporary variable

        dbh = tree%diam_bht
        ht = tree%forska_ht
        bulk_density = tree%spec_ptr%wood_bulk_dens
        beta = tree%spec_ptr%beta

        temp = 1E-4*PI*0.25*B_TO_C*bulk_density*beta/(beta + 2.0)
        stem_biomass_c = temp*dbh*dbh*ht

    end function stem_biomass_c

    !:.........................................................................:

    real function branch_biomass_c(tree)
        !
        !  Calculates branch biomass of a tree (tC)
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure
        !

        ! Data dictionary: calling arguments
        class(TreeData), intent(in) :: tree ! Tree object

        ! Data dictionary: local variables
        real :: dcbb         ! Diameter at clear branch bole height (cm)
        real :: hcbb         ! Clear branch bole height (m)
        real :: ht           ! Height (m)
        real :: bulk_density ! Wood bulk density (t/m3)
        real :: beta         ! Stem shape parameter
        real :: temp         ! Temporary variable

        dcbb = tree%diam_canht
        hcbb = tree%canopy_ht
        ht = tree%forska_ht
        bulk_density = tree%spec_ptr%wood_bulk_dens
        beta = tree%spec_ptr%beta

        temp = 1E-4*PI*0.25*B_TO_C*bulk_density*(2.0/(beta + 2.0)-0.33)
        branch_biomass_c = temp*dcbb*dcbb*(ht-hcbb)

    end function branch_biomass_c

    !:.........................................................................:

    subroutine get_diam_category(self, diam_category)
        !
        !  Gets the diameter class of a tree
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    07/27/12     K. Holcomb          Original code
        !    01/28/21     A. C. Foster        Updated to use parameters for DBH
        !                                       bins
        !

        ! Data dictionary: calling arguments
        class(TreeData),         intent(in)  :: self          ! Tree object
        integer, dimension(NHC), intent(out) :: diam_category ! List with 1 in correct bin

        ! Data dictionary: local variables
        real    :: dbh ! Diameter at breast height (cm)
        integer :: i   ! Looping index

        dbh = self%diam_bht

        ! Reset to 0
        diam_category = 0

        ! Find correct bin
        do i = 1, NHC - 1
            if (dbh > DBC_MIN(i) .and. dbh <= DBC_MAX(i)) then
                diam_category(i) = 1
            end if
        end do
        if (dbh > DBC_MIN(NHC)) then
            diam_category(NHC) = 1
        end if

    end subroutine get_diam_category

    !:.........................................................................:

    subroutine write_tree_csv(self, tree_unit)
        !
        !  Helper subroutine for writing tree-level data to the output Tree file
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/16     A. C. Foster        Original code
        !

        use csv_file

        ! Data dictionary: calling arguments
        class(TreeData), intent(in) :: self      ! Tree object
        integer,         intent(in) :: tree_unit ! Unit number for output tree file

        ! Data dictionary: local variables
        real    :: tla ! Leaf area (m2)
        integer :: i   ! Looping index

        ! Get leaf area
        tla = leaf_area(self)

        ! Write out attributes
        call csv_write(tree_unit, trim(adjustl(self%spec_ptr%genus_name)),     &
            .false.)
        call csv_write(tree_unit, self%spec_ptr%unique_id, .false.)
        call csv_write(tree_unit, self%row, .false.)
        call csv_write(tree_unit, self%col, .false.)
        call csv_write(tree_unit, self%tree_age, .false.)
        call csv_write(tree_unit, self%diam_bht, .false.)
        call csv_write(tree_unit, self%diam_canht, .false.)
        call csv_write(tree_unit, self%forska_ht, .false.)
        call csv_write(tree_unit, self%canopy_ht, .false.)
        call csv_write(tree_unit, self%leaf_bm, .false.)
        call csv_write(tree_unit, tla, .false.)
        call csv_write(tree_unit, self%biomC, .false.)
        do i = 1, FC_NUM-1
            call csv_write(tree_unit, self%env_resp(i), .false.)
        end do
        call csv_write(tree_unit, self%env_resp(FC_NUM), .true.)

    end subroutine write_tree_csv

    !:.........................................................................:

end module Tree
