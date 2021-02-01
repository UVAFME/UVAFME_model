module Species

!*******************************************************************************
  !
  ! The Species module contains the definition of the SpeciesData type,
  ! which holds attributes pertinent to all members of a species, and
  ! associated procedures.
  !
  !
!*******************************************************************************

    use Constants
    implicit none

    ! Define SpeciesData type
    type SpeciesData
        character(len = MAX_NLEN)  :: genus_name     ! Genus name
        character(len = MAX_NLEN)  :: taxonomic_name ! Taxonomic name
        character(len = 8)         :: unique_id      ! 8-character unique species id
        character(len = MAX_NLEN)  :: common_name    ! Common name
        integer                    :: shade_tol      ! Relative shade tolerance (1-5; 5 = least tolerant)
        integer                    :: lownutr_tol    ! Relative low nutrient tolerance (1-3; 3 = least tolerant)
        integer                    :: stress_tol     ! Relative stress tolerance (1-5; 5 = least tolerant)
        integer                    :: age_tol        ! Probability of surviving to max age (1-3; 3 = least likely)
        integer                    :: drought_tol    ! Relative drought tolerance (1-6; 6 = least tolerant)
        integer                    :: flood_tol      ! Relative flood tolernace (1-7; 7 = least tolerant)
        integer                    :: perm_tol       ! Relative permafrost tolerane (1-2; 2 = least tolerant)
        integer                    :: fire_regen     ! Relative ability to reproduce following fire
                                                     ! (1-2: fire beneficial; 3: 4-6 fire detrimental)
        integer                    :: litter_class   ! Litter class
        integer                    :: org_tol        ! Relative ability to germinate on deep soil (1-3; 3 = least able)
        integer                    :: recr_age       ! Minimum age for reproduction (years)
        real                       :: max_age        ! Average maximum age (years)
        real                       :: max_diam       ! Average maximum diameter (cm)
        real                       :: max_ht         ! Average maximum height (m)
        real                       :: wood_bulk_dens ! Wood bulk density (t/m3)
        real                       :: rootdepth      ! Average rooting depth (m)
        real                       :: leafdiam_a     ! Scalar for leaf area-diameter relationship
        real                       :: leafarea_c     ! Scalar for leaf area:biomass relationship
        real                       :: deg_day_min    ! Minimum growing degree-days (>5degC)
        real                       :: deg_day_opt    ! Optimum growing degree-days (>5degC)
        real                       :: deg_day_max    ! Maximum growing degree-days (>5degC)
        real                       :: seed_surv      ! Proportion of seebank lost annually (0-1)
        real                       :: seedling_surv  ! Proportion of seedling bank lost annually (0-1)
        real                       :: invader        ! Seed rain from outside the plot (seeds/m2)
        real                       :: bark_thick     ! Bark thickness (cm bark/cm DBH)
        real                       :: seed_num       ! Seed rain from within the plot (seeds/m2)
        real                       :: sprout_num     ! Regeneration from sprouting (sprouts/m2)
        real                       :: s, g, beta     ! Growth parameters
        real                       :: dbh_min        ! Minimum diameter increment before "stressed" (cm)
        logical                    :: conifer        ! Is species a conifer?
        logical                    :: layering       ! Can species reproduce by layering?
    end type SpeciesData

    interface assignment(=)
        module procedure copy_species
    end interface

    private copy_species

contains

	!:.........................................................................:

	subroutine initialize_species(self, genus_name, taxonomic_name, unique_id, &
        common_name, shade_tol, lownutr_tol, stress_tol, age_tol,              &
        drought_tol, flood_tol, perm_tol, org_tol, bark_thick,                 &
        fire_regen, max_age, max_diam, max_ht, wood_bulk_dens, rootdepth,      &
        leafdiam_a, leafarea_c, deg_day_min, deg_day_opt, deg_day_max,         &
        seedling_surv, invader, seed_num, sprout_num, layering, seed_surv, s,  &
        g, beta, conifer, litter_class, recr_age, dbh_min)
        !
        !  Initializes a species object using input parameters and initial
        !  conditions.
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    07/27/12     K. Holcomb           Original Code
        !

        ! Data dictionary: constants
        ! Scalar for modifying g by shade tolerance
        real, dimension(5), parameter :: SS = [1.1, 1.15, 1.2, 1.23, 1.25]
        ! Scalar for modifying leafarea_c by shade tolerance
        real, dimension(5), parameter :: ADJUST = [1.5, 1.55, 1.6, 1.65, 1.7]

        ! Data dictionary: calling argumnets
		class(SpeciesData),   	   intent(inout)  :: self           ! Species object
		character(len = MAX_NLEN), intent(in)     :: genus_name     ! Genus name
		character(len = MAX_NLEN), intent(in)     :: taxonomic_name ! Taxonomic name
		character(len = 8),        intent(in)     :: unique_id      ! 8-character unique species id
		character(len = MAX_NLEN), intent(in)     :: common_name    ! Common name
		integer,                   intent(in)     :: shade_tol      ! Relative shade tolerance (1-5; 5 = least tolerant)
		integer,                   intent(in)     :: lownutr_tol    ! Relative low nutrient tolerance (0-3; 3 = least tolerant)
		integer,                   intent(in)     :: stress_tol     ! Relative stress tolerance (1-5; 5 = least tolerant)
		integer,                   intent(in)     :: age_tol        ! Probability of surviving to max age (1-3; 3 = least likely)
        integer,                   intent(in)     :: fire_regen     ! Relative ability to reproduce following fire
                                                                    ! (1-2: fire beneficial; 3: 4-6 fire detrimental)
		integer,                   intent(in)     :: drought_tol    ! Relative drought tolerance (1-6; 6 = least tolerant)
		integer,                   intent(in)     :: flood_tol      ! Relative flood tolernace (1-7; 7 = least tolerant)
		integer,                   intent(in)     :: perm_tol       ! Relative permafrost tolerane (1-2; 2 = least tolerant)
		integer,                   intent(in)     :: org_tol        ! Relative ability to germinate on deep soil (1-3; 3 = least able)
        integer,                   intent(in)     :: recr_age       ! Minimum age for reproduction (years)
		integer,                   intent(in)     :: litter_class   ! Litter class
		real,                      intent(in)     :: max_age        ! Average maximum age (years)
        real,                      intent(in)     :: max_diam       ! Average maximum diameter (cm)
		real,                      intent(in)     :: max_ht         ! Average maximum height (m)
        real,                      intent(in)     :: rootdepth      ! Average rooting depth (m)
		real,                      intent(in)     :: wood_bulk_dens ! Wood bulk density (t/m3)
		real,                      intent(in)     :: leafdiam_a     ! Scalar for leaf area-diameter relationship
		real,                      intent(in)     :: leafarea_c     ! Scalar for leaf area:biomass relationship
		real,                      intent(in)     :: deg_day_min    ! Minimum growing degree-days (>5degC)
		real,                      intent(in)     :: deg_day_opt    ! Optimum growing degree-days (>5degC)
		real,                      intent(in)     :: deg_day_max    ! Maximum growing degree-days (>5degC)
		real,                      intent(in)     :: seed_surv      ! Proportion of seebank lost annually (0-1)
		real,                      intent(in)     :: seedling_surv  ! Proportion of seedling bank lost annually (0-1)
		real,                      intent(in)     :: invader        ! Seed rain from outside the plot (seeds/m2)
        real,                      intent(in)     :: bark_thick     ! Bark thickness (cm bark/cm DBH)
		real,                      intent(in)     :: seed_num       ! Seed rain from within the plot (seeds/m2)
		real,                      intent(in)     :: sprout_num     ! Regeneration from sprouting (sprouts/m2)
		real,                      intent(in)     :: s, g, beta     ! Growth parameters
        real,                      intent(in)     :: dbh_min        ! Minimum diameter increment before "stressed" (cm)
		logical,                   intent(in)     :: conifer        ! Can species reproduce by layering?
        logical,                   intent(in)     :: layering       ! Can species reproduce by layering?




        ! Set instance variables from input file
        self%genus_name     = genus_name
        self%taxonomic_name = taxonomic_name
        self%common_name    = common_name
        self%unique_id      = unique_id
        self%max_age        = max_age
        self%max_diam       = max_diam
        self%max_ht         = max_ht
        self%rootdepth      = rootdepth
        self%wood_bulk_dens = wood_bulk_dens
        self%deg_day_min    = deg_day_min
        self%deg_day_max    = deg_day_max
        self%deg_day_opt    = deg_day_opt
        self%shade_tol      = shade_tol
        self%lownutr_tol    = lownutr_tol
        self%drought_tol    = drought_tol
        self%bark_thick     = bark_thick
        self%fire_regen     = fire_regen
        self%flood_tol      = flood_tol
        self%perm_tol       = perm_tol
        self%org_tol        = org_tol
        self%stress_tol     = stress_tol
        self%age_tol        = age_tol
        self%conifer        = conifer
        self%layering       = layering
        self%invader        = invader
        self%seed_num       = seed_num
        self%sprout_num     = sprout_num
        self%seed_surv      = seed_surv
        self%seedling_surv  = seedling_surv
        self%s              = s
        self%beta           = beta
        self%litter_class   = litter_class
        self%recr_age       = recr_age
        self%dbh_min        = dbh_min

		! Convert to tC/ha
        self%leafarea_c = leafarea_c/HEC_TO_M2

        ! Adjust for shade tolerance
        self%leafdiam_a = leafdiam_a*ADJUST(shade_tol)
        self%g = g*SS(shade_tol)

	end subroutine initialize_species

	!:.........................................................................:

	subroutine copy_species(self, species_data)
        !
        !  Copies a instance variables from one object to another
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    07/27/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
		class(SpeciesData),  intent(out)  :: self         ! New species object
		class(SpeciesData),  intent(in)   :: species_data ! Old species object

        ! Copy attributes
        self%genus_name     = species_data%genus_name
        self%taxonomic_name = species_data%taxonomic_name
        self%common_name    = species_data%common_name
        self%unique_id      = species_data%unique_id
        self%max_ht         = species_data%max_ht
        self%max_age        = species_data%max_age
        self%max_diam       = species_data%max_diam
        self%rootdepth      = species_data%rootdepth
        self%wood_bulk_dens = species_data%wood_bulk_dens
        self%leafarea_c     = species_data%leafarea_c
        self%leafdiam_a     = species_data%leafdiam_a
        self%deg_day_min    = species_data%deg_day_min
        self%deg_day_max    = species_data%deg_day_max
        self%deg_day_opt    = species_data%deg_day_opt
        self%shade_tol      = species_data%shade_tol
        self%lownutr_tol    = species_data%lownutr_tol
        self%drought_tol    = species_data%drought_tol
        self%bark_thick     = species_data%bark_thick
        self%fire_regen     = species_data%fire_regen
        self%flood_tol      = species_data%flood_tol
        self%perm_tol       = species_data%perm_tol
        self%org_tol        = species_data%org_tol
        self%stress_tol     = species_data%stress_tol
        self%age_tol        = species_data%age_tol
        self%conifer        = species_data%conifer
        self%layering       = species_data%layering
        self%litter_class   = species_data%litter_class
        self%invader        = species_data%invader
        self%seed_num       = species_data%seed_num
        self%sprout_num     = species_data%sprout_num
        self%seed_surv      = species_data%seed_surv
        self%seedling_surv  = species_data%seedling_surv
        self%s              = species_data%s
        self%g              = species_data%g
        self%beta           = species_data%beta
        self%recr_age       = species_data%recr_age
        self%dbh_min        = species_data%dbh_min

	end subroutine copy_species

	!:.........................................................................:

	real function light_rsp(self, al)
        !
        !  Computes growth response to available light
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure
        !

        ! Data dictionary: constants
        ! Light response parameters - based on shade tolerance
        real, dimension(5), parameter :: LIGHT_C1 = [1.0, 1.31, 1.62, 1.93,    &
            2.24]
        real, dimension(5), parameter :: LIGHT_C2 = [4.64, 3.76, 2.888, 2.012, &
            1.136]
        real, dimension(5), parameter :: LIGHT_C3 = [0.05, 0.0575, 0.065,      &
            0.0725, 0.08]

        ! Data dictionary: calling arguments
		class(SpeciesData), intent(in) :: self ! Species object
		real,               intent(in) :: al   ! Available light (0-1)

        ! Data dictionary: local variables
		real    :: flight ! Growth response to available light (0-1)
		integer :: kt     ! Relative shade tolerance (1-5; 5 = least tolerant)

        ! Get shade tolerance
		kt = self%shade_tol

        ! Calculate growth response
		flight = LIGHT_C1(kt)*(1.0 - exp(-LIGHT_C2(kt)*(al - LIGHT_C3(kt))))
		if (flight .lt. 0.0) flight = 0.0
		if (flight .gt. 1.0) flight = 1.0

		light_rsp = flight

	end function light_rsp

	!:.........................................................................:

	real function poor_soil_rsp(n_avail, n_tol)
        !
        !  Computes quadratic poor nutrient response growth factors by poor
        !  nutrient response class and relative N availability (available N
        !  relative to N demand by all plants)
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure
        !

        ! Data dictionary: constants
        ! Nutrient response parameters - based on nutrient tolerance
        real, dimension(3), parameter :: FERT_C1 = [-0.6274, -0.2352, 0.2133]
        real, dimension(3), parameter :: FERT_C2 = [3.600, 2.771, 1.789]
        real, dimension(3), parameter :: FERT_C3 = [-1.994, -1.550, -1.014]

        ! Data dictionary: calling arguments
		real,    intent(in) :: n_avail ! Relative N availability (0-1)
		integer, intent(in) :: n_tol   ! Relative low nutrient tolerance (0-3; 3 = least tolerant)

        ! Data dictionary: local variables
	    integer :: k      ! Relative low nutrient tolerance (0-3; 3 = most tolerant)
	    real    :: navail ! Relative N availability
		real    :: fpoor  ! Growth response to low nutrients (0-1)

        ! Here the arrays are backwards so 1 is intolerant and 3 is tolerant
		! (in input csv 3 is intolerant, 1 is tolerant), so must switch
        k = 4 - n_tol

		! Max of 1.0
		navail = min(n_avail, 1.0)

        ! Calculate impact of low nutrients
		fpoor = FERT_C1(k) + FERT_C2(k)*navail + FERT_C3(k)*navail**2
		if (fpoor.le. 0.0) fpoor = 0.0
		if (fpoor.ge. 1.0) fpoor = 1.0

		poor_soil_rsp = fpoor*navail


	end function poor_soil_rsp

	!:.........................................................................:

	subroutine fire_rsp(self, fire, fc_fire)
        !
        !  Calculates species-level effect of fire on regeneration
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure
        !    01/27/21     A. C. Foster        Updated to have output be a plot
        !                                       attribute
        !

        ! Data dictionary: constants:
        ! Regeneration response parameters - based on fire regeneration tolerance
        real, dimension(6), parameter :: GAMA = [100.0, 10.0, 1.00, 0.1,       &
            0.01, 0.001]

        ! Data dictionary: calling arguments
		class(SpeciesData), intent(in)  :: self    ! Species object
		integer,            intent(in)  :: fire    ! Did a fire occur (1: yes; 0: no)
        real,               intent(out) :: fc_fire ! Regeneration response to fire

        ! Data dictionary: local variables
		real    :: resp ! Regeneration response to fire
		integer :: k    ! Relative ability to reproduce following fire
                        ! (1-2: fire beneficial; 3: 4-6 fire detrimental)

        ! Get fire regeneration tolerance
		k = self%fire_regen

		if (fire == 1) then
            ! Fire occurred - calculate response
			resp = GAMA(k)
		else
            ! No fire occurred
			resp = 1.0
		end if

		fc_fire = resp

	end subroutine fire_rsp

	!:.........................................................................:

	subroutine temp_rsp(self, x, fc_gdd)
        !
        !  Calculates species-level effect of growing degree-days on growth
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure
        !    10/10/14     A. C. Foster        Asymptotic growth response
        !    01/27/21     A. C. Foster        Updated to have output be a plot
        !                 J. K. Shuman           attribute
        !

        ! Data dictionary: calling arguments
		class(SpeciesData), intent(in)  :: self   ! Species object
		real,               intent(in)  :: x      ! Growing degree-days (>5degC)
        real,               intent(out) :: fc_gdd ! Growth response to growing degree-days (0-1)

        ! Data dictionary: local variables
		real :: ftemp ! Growth response to growing degree-days (0-1)
		real :: ddmin ! Minimum growing degree-days (>5degC)
        real :: ddopt ! Optimum growing degree-days (>5degC)
        real :: ddmax ! Maximum growing degree-days (>5degC)
		real :: a, b  ! Exponents for growing degree-day effect
        real :: tmp   ! Temporary variable

        ! Reducing OOP lookups
		ddmin = self%deg_day_min
        ddmax = self%deg_day_max
		ddopt = self%deg_day_opt

        ! Calculate exponents
		a = (ddopt - ddmin)/(ddmax - ddmin)
		b = (ddmax - ddopt)/(ddmax - ddmin)

		! Asymptotic GDD
		if (x .le. ddmin) then
            ! Below minimum - can't grow
			ftemp = 0.0
		elseif (x .ge. ddopt) then
            ! Above optimum - no effect
			ftemp = 1.0
		else
            ! Calculate effect
			tmp = ((x - ddmin)/(ddopt - ddmin))**a
			ftemp = tmp*((ddmax - x)/(ddmax - ddopt))**b
		end if

		fc_gdd = ftemp

	end subroutine temp_rsp

	!:.........................................................................:

	subroutine drought_rsp(self, drydays, fc_drought)
        !
        !  Calculates species-level effect low soil moisture on growth
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure
        !    01/27/21     A. C. Foster        Updated to have output be a plot
        !                                       attribute
        !

        ! Data dictionary: calling arguments
		class(SpeciesData), intent(in)  :: self       ! Species object
		real,               intent(in)  :: drydays    ! Drought index (0-1)
        real,               intent(out) :: fc_drought ! Growth response to drought (0-1)


		if (self%drought_tol .eq. 1) then

            ! Very drought tolerant - give extra boost
			if (self%conifer) then
				fc_drought = fdry(drydays, 1)*1.1
			else
				fc_drought = fdry(drydays, 1)*1.05
			endif
		else

			fc_drought = fdry(drydays, self%drought_tol)

		end if

        fc_drought = min(max(fc_drought, 0.0), 1.0)

	end subroutine drought_rsp

	!:.........................................................................:


	subroutine perm_rsp(p_tol, amlt, fc_perm)
        !
        !  Calculates species-level effect of permafrost on growth
        !  Adapted from Bonan 1989 Ecological Modelling 45:275-306
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    09/08/19     A. C. Foster        Original Code
        !

        ! Data dictionary: calling arguments
		integer, intent(in)  :: p_tol    ! Relative permafrost tolerane (1-2; 2 = least tolerant)
		real,    intent(in)  :: amlt     ! Active layer depth (m)
		real,    intent(out) :: fc_perm  ! Growth response to permafrost (0-1)

        ! Data dictionary: local variables
		real :: y ! Growth response to permafrost (0-1)

		if (amlt .le. 0.6) then
            ! Permafrost fairly impactful

            if (p_tol .eq. 1) then
				y = 0.7*amlt
			else
				y = 0.494*amlt
			end if
		else if (amlt .gt. 0.6 .and. amlt .le. 1.0) then
            ! Permafrost somewhat impactful
			if (p_tol .eq. 1) then
				y = 0.9*amlt
			else
				y = 0.8*amlt
			end if
		else if (amlt .gt. 1.0) then
            ! Permafrost not impactful
			y = 1.0
		end if

        ! Must be between 0.0 and 1.0
		if (y .ge. 1.0) y = 1.0
		if (y .le. 0.0) y = 0.0

		fc_perm = y

	end subroutine perm_rsp

	!:.........................................................................:

	real function fdry(dryday, k)
        !
        !  Helper function for drought_rsp subroutine
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !

        ! Data dictionary: constants
        ! Drought response parameters - based on drought tolerance
        real, dimension(6), parameter :: GAMA = [0.533, 0.4364, 0.3398,        &
            0.2432, 0.1466, 0.050]

        ! Data dictionary: calling arguments
		real,    intent(in)  :: dryday ! Drought index (0-1)
		integer, intent(in)  :: k      ! Relative drought tolerance (1-6; 6 = least tolerant)

        ! Data dictionary: local variables
		real :: tmp ! Temporary variable

        ! Calculate growth response
		tmp = max(GAMA(k) - dryday, 0.0)
		fdry = (tmp/GAMA(k))**0.5

	end function fdry

	!:.........................................................................:

end module Species
