module Tree

   use Constants
   use Species
   use Random

   implicit none

!*******************************************************************************
  !
  !The Tree module contains attributes and procedures relevant to tree growth
  !and mortality
  !
!*******************************************************************************

	real, parameter              :: tC = 3.92699e-5
	real, parameter              :: beta = 1.0

	type, extends(SpeciesData)   :: TreeData
		real                     :: diam_max
		real                     :: diam_bht
		real                     :: diam_canht
		real                     :: canopy_ht
		real                     :: foret_ht
		real                     :: forska_ht
		real                     :: leaf_bm
		real                     :: biomC
		real                     :: stemC
		real                     :: twigC
		real                     :: rootC
		real                     :: biomN
		real                     :: CK
		real                     :: row, col
		real, dimension(6)       :: env_resp
		integer                  :: species_index
		integer                  :: mort_count
		integer                  :: tree_age
		integer                  :: stressor
		logical                  :: mort_marker
	end type TreeData

contains

!===============================================================================
! Methods
!===============================================================================

	!:.........................................................................:

	subroutine initialize_tree(self, tree_species, si)
		!initializes tree with initial conditions
		!Inputs/Outputs:
		!	self:         tree instance
		!Inputs:
		!	tree_species: species instance
		!	si:	          species index

		class(TreeData),             intent(inout) :: self
		type(SpeciesData), optional, intent(in)    :: tree_species
		integer,           optional, intent(in)    :: si

		!constructor
		self%diam_bht = 0.0
		self%diam_canht = 0.0
		self%canopy_ht = std_ht
		self%foret_ht = 1.0
		self%forska_ht = 1.0
		self%biomC = 0.0
		self%twigC = 0.0
		self%stemC = 0.0
		self%rootC = 0.0
		self%biomN = 0.0
		self%mort_count = 0
		self%CK = 0.0
		self%tree_age = 0
		self%stressor = 0
		self%mort_marker = .false.
		self%env_resp = -999.0
		self%row = 0
		self%col = 0

		!attach species data and species index
		if (present(tree_species)) then
			self%SpeciesData = tree_species
		endif

		if ( present(si) ) then
			self%species_index = si
		else
			self%species_index = 0
		endif

	end subroutine initialize_tree

	!:.........................................................................:

	subroutine copy_tree(self, tree)
		!copies tree from old list to new
		!Inputs/Outputs:
		!	self: tree instance (new)
		!Inputs:
		!	tree: tree instance (old)

		class(TreeData), intent(inout) :: self
		class(TreeData), intent(in)    :: tree

		self%diam_max       = tree%diam_max
		self%diam_bht       = tree%diam_bht
		self%diam_canht     = tree%diam_canht
		self%canopy_ht      = tree%canopy_ht
		self%foret_ht       = tree%foret_ht
		self%forska_ht      = tree%forska_ht
		self%leaf_bm        = tree%leaf_bm
		self%biomC          = tree%biomC
		self%twigC          = tree%twigC
		self%stemC          = tree%stemC
		self%rootC          = tree%rootC
		self%biomN          = tree%biomN
		self%species_index  = tree%species_index
		self%mort_count     = tree%mort_count
		self%tree_age       = tree%tree_age
		self%CK             = tree%CK
		self%mort_marker    = tree%mort_marker
		self%stressor       = tree%stressor
		self%env_resp       = tree%env_resp
		self%row            = tree%row
		self%col            = tree%col
		self%SpeciesData    = tree%SpeciesData

	end subroutine copy_tree

	!:.........................................................................:

	subroutine update_tree(self, tree_species)
		!a helper function required to have more flexible object structure
		!Inputs/Outputs:
		!	self:         tree instance
		!Inputs:
		!	tree_species: tree species

		class(TreeData),   intent(inout) :: self
		type(SpeciesData), intent(in)    :: tree_species

		self%SpeciesData = tree_species

	end subroutine update_tree

	!:.........................................................................:

    subroutine env_stress(self, shade, minstress, stressor)
		!calculates the environmental stress on a tree (using Liebig's Law of
		  !the Minimum, and sets what the most stressful factor is
		!Author: Adrianna Foster 2016 v. 1.0
		!Inputs/Outputs:
		!	self:      tree instance
		!Inputs:
		!	shade:     tree growth response to shade (0-1)
		!Outputs:
		!	minstress: environmental stress (0-1)
		!	stressor:  which factor was most stressful
		!				1: temperature
		!				2: drought
		!				3: shade
		!				4: permafrost (set later)
		!				5: nutrients (set later)

		class(TreeData), intent(inout) :: self
		real,            intent(in)    :: shade
		real,            intent(out)   :: minstress
		integer,         intent(out)   :: stressor
		real                           :: envstress

		envstress = min(self%fc_degday, self%fc_drought, shade)
		minstress = envstress

		if (minstress .eq. self%fc_degday) then
			stressor = 1
		else if (minstress .eq. self%fc_drought) then
			stressor = 2
		else if (minstress .eq. shade) then
			stressor = 3
		endif

		self%env_resp(1) = self%fc_degday
		self%env_resp(2) = self%fc_drought
		self%env_resp(3) = shade

	end subroutine env_stress

	!:.........................................................................:

	subroutine stem_shape(tree)
		!calculates diameter at clear branch bole height (m)
		!basal diameter = dhshape(0.0, h, dbh)
		!diameter at canopy height = dhshape(hc, h, dbh)
		!Note: mentioned above "beta" can be 1, 1.5 ,2 , 3, but right now beta
		  !is set to 1 for all these routines.
		!Inputs/Outputs:
		!	self:      tree instance

		class(TreeData), intent(inout) :: tree

		real                           :: hc, h, dbh, dcbb

		hc = tree%canopy_ht; h = tree%forska_ht; dbh = tree%diam_bht

		if (h .le. hc .or. h .le. std_ht) then
			dcbb = dbh
		else
			dcbb = ((h - hc)/(h - std_ht))**(1.0/beta)*dbh
		end if

		tree%diam_canht = dcbb

	end subroutine stem_shape

	!:.........................................................................:

	subroutine forska_height(tree)
		!calculates total tree height using height-dbh relation equation from
		  !Forska model
		!Author: Yan Xiaodong & Hank Shugart 2005 v. 1.0
		!Inputs/Outputs:
		!	tree: tree instance

		class(TreeData), intent(inout) :: tree

		real                           :: d, par, hmax
		real                           :: delta_ht

		d = tree%diam_bht
		hmax = tree%max_ht; par = tree%arfa_0

		delta_ht = hmax - std_ht
		tree%forska_ht = std_ht + delta_ht*(1.0 - exp(-(par*d/delta_ht)))

	end subroutine forska_height

	!:.........................................................................:

	subroutine biomass_c(tree)
		!calculates trunk and branch biomass C (tonnes C) of tree
		!Author: Yan Xiaodong & Hank Shugart 2005 v. 1.0
		!Inputs/Outputs:
		!	tree: tree instance

		class(TreeData), intent(inout) :: tree

		type (TreeData)                :: tree_0
		real                           :: d, h, hc, hrt, bulk
		real                           :: stembc, twigbc
		real                           :: abovegr_c, root_c

		d = tree%diam_bht; h = tree%forska_ht; hc = tree%canopy_ht
		hrt = tree%rootdepth; bulk = tree%wood_bulk_dens

		call copy_tree(tree_0, tree)
		tree_0%canopy_ht = 0.0

		call stem_shape(tree)
		call stem_shape(tree_0)

		stembc = stem_biomass_c(tree_0)
		twigbc = twig_biomass_c(tree)

		abovegr_c = stembc + twigbc
		root_c    = stembc*hrt/h + twigbc/2.0

		tree%biomC = (abovegr_c + root_c)
		tree%stemC = stembc
		tree%twigC = twigbc
		tree%rootC = root_c

	end subroutine biomass_c

	!:.........................................................................:

	subroutine biomass_n(tree)
		!calculates stem and branch N biomass
		!Author: Yan Xiaodong & Hank Shugart 2005 v. 1.0
		!Inputs/Outputs:
		!	tree: tree instance

		class(TreeData), intent(inout) :: tree

		tree%biomN = tree%biomC/stem_c_n

	end subroutine biomass_n

	!:.........................................................................:

	subroutine leaf_biomass_c(tree)
		!calculates leaf biomass of tree (tonnes C)
		!Author: Yan Xiaodong & Hank Shugart 2005 v. 1.0
		!Inputs/Outputs:
		!	tree: tree instance

		class(TreeData), intent(inout) :: tree

		tree%leaf_bm = lai_biomass_c(tree)*tree%leafarea_c*2.0

	end subroutine leaf_biomass_c

	!:.........................................................................:

	function lai_biomass_c(tree)
		!calculates leaf area of tree
		!Author: Yan Xiaodong & Hank Shugart 2005 v. 1.0
		!Inputs:
		!	tree: tree instance

		class(TreeData), intent(in) :: tree

		real                        :: lai_biomass_c
		real                        :: dc, leafd_a

		dc = tree%diam_canht; leafd_a = tree%leafdiam_a
		lai_biomass_c = dc*dc*leafd_a

	end function lai_biomass_c

	!:.........................................................................:

	subroutine age_survival(tree, a_survive)
		!death check by age. Trees can reach life sapn (max age) as a
		  !probability of 1% ('often seen'), 0.1% ('some survive'), or 0.01%
		  !('few seen'), as set up by input parameter
		!Author: Yan Xiaodong & Hank Shugart 2005 v. 1.0
		!Inputs:
		!	tree:      tree instance
		!Outputs:
		!	a_survive: does tree survive this check?

		class(TreeData), intent(in)  :: tree
		logical,         intent(out) :: a_survive

		integer                      :: k
		real                         :: agemax
		real                         :: rand_val
		real, dimension(3)           :: check

		data check /4.605, 6.908, 11.51/
		!1% ---> 4.605  0.1% ---> 6.908   0.01% ---> 11.5129

		k = tree%age_tol
		agemax = tree%max_age

		rand_val = urand()
		if (rand_val .lt. (check(k)/agemax)) then
			a_survive = .false.
		else
			a_survive = .true.
		end if

	end subroutine age_survival

	!:.........................................................................:

	subroutine growth_survival(tree, g_survive)
		!death check by growth insufficiency
		!Author: Yan Xiaodong & Hank Shugart 2005 v. 1.0
		!Inputs:
		!	tree:       tree instance
		!Outputs:
		!	g_survive:  does tree survive growth check?

		class(TreeData), intent(in) :: tree
		logical, intent(out)        :: g_survive

		integer                     :: k
		real                        :: rand_val
		real, dimension(5)          :: check

		data check /0.21, 0.24, 0.27, 0.30, 0.33/

		k = tree%stress_tol
		rand_val = urand()

		if (tree%mort_marker .and. rand_val .lt. check(k)) then
			g_survive=.false.
		else
			g_survive=.true.
		end if

	end subroutine growth_survival

	!:.........................................................................:

	subroutine fire_survival(tree, av_fuel, f_survive)
		!calculates survival of tree from fire event
		!Authors: Adrianna Foster & Jacquelyn Shuman 2015 v. 1.0
		!Changelog:
		!	Adrianna Foster updated to 1.1, updated to equations from
		!		Shumacher et al. 2007
		!Inputs/Outputs:
		!	tree:       tree instance
		!Inputs:
		!	av_fuel:    available fuel for burning (t/ha)
		!Outputs:
		!	f_survive:  does tree survive fire?


		class(TreeData), intent(inout) :: tree
		real,            intent(in)    :: av_fuel
		logical,         intent(out)   :: f_survive

		logical                        :: fkill
		real                           :: Pfire, dbh_eff
		real                           :: rand_val
		real, parameter                :: BT = 0.063
		real, parameter                :: ck1 = 0.21111
		real, parameter                :: ck2 = -0.00445

		if (tree%diam_bht .ge. 40.0) then
			dbh_eff = 40.0
		else
			dbh_eff = tree%diam_bht
		end if

		tree%CK = min(100.0, 100.0*(ck1 + ck2*dbh_eff)*av_fuel)
		if (tree%CK .le. 0.0001) tree%CK = 0.0


		Pfire = 1/(1 + exp(-1.466 +                                            &
		                1.91*(tree%bark_thick*tree%diam_bht) -                 &
		                0.1775*((tree%bark_thick*tree%diam_bht)**2) -          &
		                0.000541*tree%CK**2))

		rand_val = urand()

		if (rand_val .lt. Pfire) then
			fkill = .true.
		else
			fkill = .false.
		endif

		if (fkill .eq. .true.) then
			f_survive = .false.
		else
			f_survive = .true.
		endif

	end subroutine fire_survival

	!:.........................................................................:

	subroutine wind_survival(tree, wcat, w_survive)
		!calculates mortality of tree from a windthrow event
		!Author: Adrianna Foster 2015 v. 1.0
		!Inputs:
		!	tree:       tree instance
		!	wcat:       wind intensity category
		!Outputs:
		!	w_survive:  does tree survive windthrow?


		class(TreeData), intent(in)  :: tree
		real,            intent(in)  :: wcat
		logical,         intent(out) :: w_survive

		real                         :: logit, pwind, rand_val


		!from Rich et al 2007 J. of Ecol.
		if (wcat .ge. 0.1) then
			logit = 0.75*log(tree%diam_bht)
			pwind = 1.0/(1.0+exp(-logit))
			rand_val=urand()
			if (rand_val .lt. pwind) then
				w_survive = .false.
			else
				w_survive = .true.
			endif
		else
			w_survive = .true.
		end if

	end subroutine wind_survival

	!:.........................................................................:

	subroutine max_growth(tree)
		!calculates optimal DBH increment growth of tree
		!!Author: Yan Xiaodong & Hank Shugart 2005 v. 1.0
		!Inputs/Outputs:
		!	tree: tree instance

		class(TreeData), intent(inout) :: tree

		real                           :: d, dc, h, dm, hm, g, arfa0
		real                           :: beginning

		dc = tree%diam_canht; d = tree%diam_bht ;h = tree%forska_ht
		dm = tree%max_diam; hm = tree%max_ht; g = tree%g
		arfa0 = tree%arfa_0

		beginning = arfa0*exp(-arfa0*d/(hm - std_ht))*d
		tree%diam_max = g*d*(1.0 - d*h/dm/hm)/(2.0*h + beginning)

	end subroutine max_growth

	!:.........................................................................:

	function stem_biomass_c(tree)
		!calculates stem biomass C (tonnes C) of tree
		!note: root stem biomass C = bstemc(bd,1.0,bulk_density)
		!Author: Yan Xiaodong & Hank Shugart 2005 v. 1.0
		!Inputs:
		!	tree: tree instance

		class(TreeData), intent(in) :: tree

		real                        :: stem_biomass_c
		real                        :: bd, h, bulk_density
		real                        :: yxd

		!bd is basal diameter but when we call this routine we use an object
		  !with a "canopy height" of zero.
		bd = tree%diam_canht; h = tree%forska_ht
		bulk_density = tree%wood_bulk_dens

		yxd = tC*bulk_density*beta/(beta + 2.0)
		stem_biomass_c = yxd*bd*bd*h*0.90

	end function stem_biomass_c

	!:.........................................................................:

	function twig_biomass_c(tree)
		!calculates twig C biomass of tree
		!Author: Yan Xiaodong & Hank Shugart 2005 v. 1.0
		!Inputs:
		!	tree: tree instance

		class(TreeData), intent(in) :: tree

		real                        :: twig_biomass_c
		real                        :: dc, hc, h, bulk_density
		real                        :: yxd

		dc = tree%diam_canht; hc = tree%canopy_ht; h = tree%forska_ht
		bulk_density = tree%wood_bulk_dens

		yxd = tC*bulk_density*(2.0/(beta + 2.0)-0.33)
		twig_biomass_c = yxd*dc*dc*(h-hc)

	end function twig_biomass_c

	!:.........................................................................:

	subroutine get_diam_category(self, diam_category)
		!gets the diameter class of the tree
		!Author: Katherine Holcomb 2012 v. 1.0
		!Inputs:
		!	self:           tree instance
		!Outputs:
		!	diam_category:  list with 1 in the correct DBH category

		class(TreeData),         intent(in)  :: self
		integer, dimension(NHC), intent(out) :: diam_category

		real                                 :: dm

		dm = self%diam_bht
		diam_category = 0

		if (dm >= 0.5 .and. dm <= 5.0) then
			diam_category(1) = 1
		else if (dm > 5.0 .and. dm <= 10.0) then
			diam_category(2) = 1
		else if (dm > 10.0 .and. dm <= 20.0) then
			diam_category(3) = 1
		else if (dm > 20.0 .and. dm <= 30.0) then
			diam_category(4) = 1
		else if (dm > 30.0 .and. dm <= 40.0) then
			diam_category(5) = 1
		else if (dm > 40.0 .and. dm <= 50.0) then
			diam_category(6) = 1
		else if (dm > 50.0 .and. dm <= 60.0) then
			diam_category(7) = 1
		else if (dm > 60.0 .and. dm <=70.0) then
			diam_category(8) = 1
		else if (dm > 70.0 .and. dm <= 80.0) then
			diam_category(9) = 1
		else if (dm > 80.0 .and. dm <= 90.0) then
			diam_category(10) = 1
		else if (dm > 90.0) then
			diam_category(11) = 1
		endif

	end subroutine get_diam_category

	!:.........................................................................:

	subroutine write_tree_csv(self, tree_unit)
		use csv_file
		!heper function for writing tree-level data to Tree_Data.csv
		!Author: Adrianna Foster 2016 v. 1.0
		!Inputs:
		!	self:       tree instance
		!	tree_unit:  file unit number for Tree_Data.csv

		class(TreeData), intent(in) :: self
		integer,         intent(in) :: tree_unit

		call csv_write(tree_unit, trim(adjustl(self%genus_name)), .false.)
		call csv_write(tree_unit, self%unique_id, .false.)
		call csv_write(tree_unit, self%row, .false.)
		call csv_write(tree_unit, self%col, .false.)
        call csv_write(tree_unit, self%tree_age, .false.)
		call csv_write(tree_unit, self%diam_bht, .false.)
		call csv_write(tree_unit, self%forska_ht, .false.)
		call csv_write(tree_unit, self%canopy_ht, .false.)
		call csv_write(tree_unit, self%leaf_bm, .false.)
		call csv_write(tree_unit, self%biomC, .false.)
        call csv_write(tree_unit, self%env_resp(1), .false.)
        call csv_write(tree_unit, self%env_resp(2), .false.)
        call csv_write(tree_unit, self%env_resp(3), .false.)
        call csv_write(tree_unit, self%env_resp(4), .false.)
        call csv_write(tree_unit, self%env_resp(5), .false.)
        call csv_write(tree_unit, self%env_resp(6), .true.)

	end subroutine write_tree_csv

	!:.........................................................................:

end module Tree
