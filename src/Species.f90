module Species
  use Constants

  implicit none

!*******************************************************************************
  !
  !The Species module contains the definition of the SpeciesData type,
  !which holds attributes pertinent to all members of a species, and
  !associated procedures.
  !
  !
  !Methods:
  !  initialize_species
  !  copy_species: overloads assignments
  !  light_rsp(self, al)
  !  temp_rsp(self, x)
  ! 	verified against original 07/27/2012
  !  drought_rsp(self, x)
  !     verified against original 07/27/2012
  !  poor_soil_rsp(self, x)
  !     verified against original 07/27/2012
  !  fire_rsp(self, x)
  !  flood_rsp(self, x)
  !
!*******************************************************************************

  type SpeciesData
     character(len = MAX_NLEN)  :: genus_name
     character(len = MAX_NLEN)  :: taxonomic_name
     character(len = 8)         :: unique_id
     character(len = MAX_NLEN)  :: common_name
     integer                    :: genus_id, species_id
     integer                    :: shade_tol, lownutr_tol, stress_tol
     integer                    :: age_tol, drought_tol, flood_tol
     integer                    :: perm_tol
     integer                    :: fire_tol, fire_regen
     integer                    :: litter_class, seed_moist
     integer                    :: org_tol, recr_age
     real                       :: max_age, max_diam, max_ht
     real                       :: wood_bulk_dens
     real                       :: rootdepth
     real                       :: leafdiam_a, leafarea_c
     real                       :: deg_day_min, deg_day_opt, deg_day_max
     real                       :: seed_surv, seedling_lg
     real                       :: invader, bark_thick
     real                       :: seed_num, sprout_num
     real                       :: arfa_0, g
     real                       :: fc_fire, fc_wind
     real                       :: fc_degday, fc_drought, fc_flood
     logical                    :: conifer, fire_kill
     logical                    :: layering
  end type SpeciesData

  interface assignment(=)
     module procedure copy_species
  end interface

  private copy_species

contains

!-------------------------------------------------------------------------------
! Methods
!-------------------------------------------------------------------------------

	!:.........................................................................:

	subroutine initialize_species(self, species_id, genus_name,                &
								taxonomic_name, unique_id, common_name,        &
								genus_id, shade_tol, lownutr_tol, stress_tol,  &
								age_tol, drought_tol, flood_tol, perm_tol,     &
								org_tol, bark_thick, fire_regen, max_age,      &
								max_diam, max_ht, wood_bulk_dens, rootdepth,   &
								leafdiam_a, leafarea_c, deg_day_min,           &
								deg_day_opt, deg_day_max, seedling_lg,         &
								invader, seed_num, sprout_num, layering,       &
								seed_surv, arfa_0, g, conifer, litter_class,   &
								recr_age)
		!initializes species instances
		!Author: Katherine Holcomb, 2012 v. 1.0
		!Inputs/Outputs:
		!   self:           species data type
		!Inputs:
		!   genus_name:     genus of species
		!   taxonomic_name: latin name of species
		!   unique_id:      8 character ID built from genus and species names
		!   common_name:    english commmon name of species
		!   genus_id:       genus ID number
		!	species_id:     species ID number (within genus)
		!   shade_tol:      relative shade tolerance (1-5; 5 = least tolerant)
		!   lownutr_tol:    relative nutrient tolerance (1-3; 3 = least tolerant)
		!   stress_tol:     relative stress tolerance (1-3; 3 = least tolerant)
		!   age_tol:        probability of reaching max_age (1-3; 3 = least likely)
		!   fire_regen:     relative ability to regenerate after fire (1-6;
		!						1 = benefit from fire; 3 = no effect;
		!						6 = fire detrimental to regeneration
		!   drought_tol:    relative drought tolerance (1-6; 6 = least tolerant)
		!   flood_tol:      relative flood tolerance (1-5; 5 = least tolerant)
		!	perm_tol:       permafrost tolerance (1 = tolerant; 2 = intolerant)
		!   org_tol:        ability to regenerate on organic layer (1-3;
		!						3 = least tolerant)
		!	litter_class:   litter class of species
		!	max_age:        average maximum age (years)
		! 	max_diam:       average maximum DBH (cm)
		!	max_ht:         average maximum height (m)
		!	rootdepth:      rootdepth of species (m)
		!	wood_bulk_dens: bulk density of wood (tonnes/m3)
		!	leafdiam_a:     scalar parameter for leaf LA-Dcbb relationship
		!   leafarea_c:     leaf area ratio
		!   deg_day_min:    minimum growing degree-days
		!   deg_day_opt:    optimum growing degree-days
		!	deg_day_max:    maximum growing degree-days
		! 	seed_surv:      proportion seedbank lost annually (0 to 1)
		!   seedling_lg:    proportion seedling bank lost annually (0 to 1)
		!   invader:        seed rain from outside plot (seeds/m2)
		!	seed_num:       seed rain within plot (seeds/m2)
		!   sprout_num:     average sprouts per individual
		!   arfa_0:         growth parameter (i.e. s) (m/cm)
		!   g:              growth parameter
		!   conifer:        evergreen (1) or deciduous (0)
		!   layering:       ability of species to reproduce by layer (0 or 1)

		class(SpeciesData),   	   intent(inout)  :: self
		character(len = MAX_NLEN), intent(in)     :: genus_name
		character(len = MAX_NLEN), intent(in)     :: taxonomic_name
		character(len = 8),        intent(in)     :: unique_id
		character(len = MAX_NLEN), intent(in)     :: common_name
		integer,                   intent(in)     :: genus_id
		integer,                   intent(in)     :: species_id
		integer,                   intent(in)     :: shade_tol
		integer,                   intent(in)     :: lownutr_tol
		integer,                   intent(in)     :: stress_tol
		integer,                   intent(in)     :: age_tol, fire_regen
		integer,                   intent(in)     :: drought_tol
		integer,                   intent(in)     :: flood_tol
		integer,                   intent(in)     :: perm_tol
		integer,                   intent(in)     :: org_tol, recr_age
		integer,                   intent(in)     :: litter_class
		real,                      intent(in)     :: max_age, max_diam
		real,                      intent(in)     :: max_ht, rootdepth
		real,                      intent(in)     :: wood_bulk_dens
		real,                      intent(in)     :: leafdiam_a
		real,                      intent(in)     :: leafarea_c
		real,                      intent(in)     :: deg_day_min
		real,                      intent(in)     :: deg_day_opt
		real,                      intent(in)     :: deg_day_max
		real,                      intent(in)     :: seed_surv
		real,                      intent(in)     :: seedling_lg
		real,                      intent(in)     :: invader, bark_thick
		real,                      intent(in)     :: seed_num
		real,                      intent(in)     :: sprout_num
		real,                      intent(in)     :: arfa_0, g
		logical,                   intent(in)     :: conifer, layering

		real, dimension(5)                        :: ss, adjust
		data ss/1.1, 1.15, 1.2, 1.23, 1.25/
		data adjust/1.5, 1.55, 1.6, 1.65, 1.7/

        self%species_id     = species_id
        self%genus_name     = genus_name
        self%taxonomic_name = taxonomic_name
        self%common_name    = common_name
        self%unique_id      = unique_id
        self%genus_id       = genus_id
        self%max_age        = max_age
        self%max_diam       = max_diam
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
        self%seedling_lg    = seedling_lg
        self%arfa_0         = arfa_0
        self%litter_class   = litter_class
        self%recr_age       = recr_age

		!adjustments
        self%leafarea_c = leafarea_c/hec_to_m2
        self%max_ht = min(max_ht, rootdepth*80.0/(1 + rootdepth))
        self%leafdiam_a = leafdiam_a*adjust(shade_tol)
        self%g = g*ss(shade_tol)

		!initialized to constants
        self%fc_fire = 0.0
        self%fc_wind = 0.0
        self%fc_degday = 0.0
        self%fc_drought = 0.0
        self%fc_flood = 0.0
        self%fire_kill = .false.

	end subroutine initialize_species

	!:.........................................................................:

	subroutine copy_species(self, species_data)
		!copies species list from old list to new
		!Author: Katherine Holcomb, 2012 v. 1.0
		!Inputs/Outputs:
		!   self:           species data type instance (new)
		!Inputs:
		!	species_data:   species data type instance (old)

		class(SpeciesData),  intent(out)  :: self
		class(SpeciesData),  intent(in)   :: species_data

        self%species_id     = species_data%species_id
        self%genus_name     = species_data%genus_name
        self%taxonomic_name = species_data%taxonomic_name
        self%common_name    = species_data%common_name
        self%unique_id      = species_data%unique_id
        self%genus_id       = species_data%genus_id
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
        self%fc_degday      = species_data%fc_degday
        self%fc_fire        = species_data%fc_fire
        self%fc_drought     = species_data%fc_drought
        self%fc_flood       = species_data%fc_flood
        self%fc_wind        = species_data%fc_wind
        self%conifer        = species_data%conifer
        self%layering       = species_data%layering
        self%litter_class   = species_data%litter_class
        self%fire_kill      = species_data%fire_kill
        self%invader        = species_data%invader
        self%seed_num       = species_data%seed_num
        self%sprout_num     = species_data%sprout_num
        self%seed_surv      = species_data%seed_surv
        self%seedling_lg    = species_data%seedling_lg
        self%arfa_0         = species_data%arfa_0
        self%g              = species_data%g
        self%recr_age       = species_data%recr_age

	end subroutine copy_species

	!:.........................................................................:


	function light_rsp(self, al)
		!computes growth response to available light
		!Inputs:
		!	self:      species data type instance
		!	al:        available light (0-1)
		!Outputs:
		!	light_rsp: growth response to available light (0-1)

		real                            :: light_rsp
		class(SpeciesData), intent(in)  :: self
		real,               intent(in)  :: al

		real                            :: flight
		integer                         :: kt

		real, dimension(5)              :: light_c1, light_c2, light_c3

		data light_c1 /1.01, 1.02, 1.11, 1.24, 1.49/
		data light_c2 /4.62, 4, 2.52, 1.78, 1.23/
		data light_c3 /0.05, 0.055, 0.07, 0.08, 0.09/

		kt = self%shade_tol

		flight = light_c1(kt)*(1.0 - exp(-light_c2(kt)*(al - light_c3(kt))))

		if (flight .lt. 0.0) flight = 0.0
		if (flight .gt. 1.0) flight = 1.0
		light_rsp = flight

	end function light_rsp

	!:.........................................................................:


	function poor_soil_rsp(n_av, n_tol)
		!computes quadratic poor nutrient response growth factors by poor
		  !nutrient response class (1 = intol, 3 = tol) and soil fertility
		!Inputs:
		!	n_req:         relative N availability (0-1)
		!	n_tol:         species' nutrient tolerance (1-3; 3 = least tolerant)
		!	species:       species instance
		!Outputs:
		!	poor_soil_rsp: growth response to nutrients (0-1)

		real                               :: poor_soil_rsp
		real,               intent(in)     :: n_av
		integer,            intent(in)     :: n_tol

	    integer                            :: k
	    real                               :: Nav
		real                               :: fpoor
		real, dimension(3)                 :: fert_c1, fert_c2, fert_c3
		data  fert_c1 /-0.6274, -0.2352, 0.2133/
		data  fert_c2 /3.600, 2.771, 1.789/
		data  fert_c3 /-1.994, -1.550, -1.014/

		!here the arrays are backwards so 1 is intolerant and 3 is tolerant
		  !(in input csv 3 is intolerant, 1 is tolerant), so must switch
		k = 4 - n_tol

		!N req is max of 1.0
		Nav = min(n_av, 1.0)

		fpoor = fert_c1(k) + fert_c2(k)*Nav + fert_c3(k)*Nav**2
		if (fpoor.le. 0.0) fpoor = 0.0
		if (fpoor.ge. 1.0) fpoor = 1.0

		poor_soil_rsp = fpoor*Nav

	end function poor_soil_rsp

	!:.........................................................................:

	subroutine fire_rsp(self, fire)
		!calculates species-level effect of fire on regeneration
		!Inputs/Outputs:
		!	self:     species instance
		!Inputs:
		!	fire:     did a fire occur? (0 or 1)
		!Outputs:
		!	fire_rsp: regeneration response to fire (0.001 - 1000)

		class(SpeciesData), intent(inout)  :: self
		integer,            intent(in)     :: fire

		real                               :: resp
		integer                            :: k
		real, dimension(6)                 :: gama
		data gama/100.0, 10.0, 1.00, 0.1, 0.01, 0.001/

		k = self%fire_regen

		if (fire == 1) then
			resp = gama(k)
		else
			resp = 1.0
		end if

		self%fc_fire = resp

	end subroutine fire_rsp

	!:.........................................................................:

	subroutine temp_rsp(self, x)
		!calculates species-level effect of GDD on growth
		!Author: Yan Xiaodong 2005 v. 1.0 - parabolic growth response
		!Changelog: Adrianna Foster & Jacquelyn Shuman 2014 v. 2.0 - updated to
		  !asymptotic growth response
		!Inputs/Outputs:
		!	self:     species instance
		!Inputs:
		!	x:        growing degree-days
		!Outputs:
		!	temp_rsp: growth response to temperature (0-1)

		class(SpeciesData), intent(inout)  :: self
		real,               intent(in)     :: x

		real                               :: ftemp
		real                               :: ddmin, ddopt, ddmax
		real                               :: a, b, tmp

		ddmin = self%deg_day_min; ddmax = self%deg_day_max
		ddopt = self%deg_day_opt

		a = (ddopt - ddmin)/(ddmax - ddmin)
		b = (ddmax - ddopt)/(ddmax - ddmin)

		!asymptotic GDD
		if (x .le. ddmin) then
			ftemp = 0.0
		elseif (x .ge. ddopt) then
			ftemp = 1.0
		else
			tmp = ((x - ddmin)/(ddopt - ddmin))**a
			ftemp = tmp*((ddmax - x)/(ddmax - ddopt))**b
		end if

		self%fc_degday = ftemp

	end subroutine temp_rsp

	!:.........................................................................:

	subroutine drought_rsp(self, drydays, drydays_base)
		!calculates species-level growth response to low soil moisture
		!Inputs/Outputs:
		!	self:         species instance
		!Inputs:
		!	drydays:      drought index (upper) (0-1)
		!	drydays_base: drought index (lower) (0-1)


		class(SpeciesData), intent(inout) :: self
		real,               intent(in)    :: drydays, drydays_base

		real                              :: fcdry1, fcdry2

		if (self%drought_tol .eq. 1) then

			if (self%conifer) then
				fcdry1 = fdry(drydays_base, 1)*0.33
			else
				fcdry1 = fdry(drydays_base, 1)*0.2
			endif

			fcdry2 = fdry(drydays, 1)
			self%fc_drought = max(fcdry1, fcdry2)

		else

			self%fc_drought = fdry(drydays, self%drought_tol)

		end if

	end subroutine drought_rsp

	!:.........................................................................:

	subroutine flood_rsp(self, floodday)
		!calculates species-level response to high soil moisture
		!currently inactive
		!Inputs/Outputs:
		!	self:       species instance
		!Inputs:
		!	flooddays:  flooding index

		class(SpeciesData), intent(inout)  :: self
		real,               intent(in)     :: floodday

		real                               :: fflood
		integer                            :: k

		real, dimension(6)                 :: gama
		data gama/0.85, 0.8, 0.75, 0.7, 0.65, 0.6/

		k = self%flood_tol

		if (floodday .le. gama(k)) then
			fflood = 1.0
		else
			fflood = 1.0 - ((floodday - gama(k)**2)/((floodday - gama(k))**2 + &
				gama(k)**2))
		end if

		if (fflood .ge. 1.0) fflood = 1.0
		if (fflood .le. 0.0) fflood = 0.0

		self%fc_flood = fflood

	end subroutine flood_rsp

	!:.........................................................................:

	subroutine perm_rsp(p_tol, amlt, resp)
		!calculates species-level response to active layer thickness
		  !modified from Bonan (1989)
		!Author: Adrianna Foster 2018 v. 1.0
		!Inputs:
		!	p_tol:  species permafrost tolerance (1: tolerant; 2: intolerant)
		!	amlt:   active layer depth (m)
		!Inputs:
		!	resp:   growth response to permafrost (0-1)

		integer, intent(in)  :: p_tol
		real,    intent(in)  :: amlt
		real,    intent(out) :: resp

		real                 :: y


		if (amlt .le. 0.6) then
			if (p_tol .eq. 1) then
				!y = 1.28*amlt
				y = 0.9*amlt
			else
				y = 0.494*amlt
			end if
		else if (amlt .gt. 0.6 .and. amlt .le. 1.0) then
			if (p_tol .eq. 1) then
				y = 1.0
			else
				y = 0.8*amlt
			end if
		else if (amlt .gt. 1.0) then
			y = 1.0
		end if

		if (y .ge. 1.0) y = 1.0
		if (y .le. 0.0) y = 0.0

		resp = y

	end subroutine perm_rsp

	!:.........................................................................:

	function fdry(dryday, k)
		!helper function for subroutine drought_rsp
		!Inputs:
		!	dryday:  drought index (0-1)
		!	k:       species drought tolerance (1-6; 6 = least tolerant)
		!Outputs:
		!	fdry:    response to drought

		real                 :: fdry
		real,    intent(in)  :: dryday
		integer, intent(in)  :: k

		real                 :: tmp
		real, dimension(6)   :: gama
		data gama/0.5, 0.45, 0.4, 0.3, 0.085, 0.025/

		tmp = max(gama(k) - dryday, 0.0)
		fdry = (tmp/gama(k))**0.5

	end function fdry

	!:.........................................................................:

end module Species
