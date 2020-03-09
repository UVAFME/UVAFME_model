
module Plot
use Constants
use Input
use Species
use Tree
use Soil

implicit none

!*******************************************************************************
  !
  !The Plot module contains the definition of the PlotData type, which holds
	!attributes of an individual plot, as well as procedures directly operating
	!on plot-level variables.
  !
!*******************************************************************************

	type PlotData
		type(TreeData),    dimension(:),   allocatable :: trees
		type(TreeData),    dimension(:),   allocatable :: deadtrees
		type(SpeciesData), dimension(:),   allocatable :: species
		type(SoilData)                                 :: soil
		real,              dimension(:),   allocatable :: avail_spec
		real,              dimension(:),   allocatable :: seedbank
		real,              dimension(:),   allocatable :: seedling
		real,              dimension(:),   allocatable :: con_light
		real,              dimension(:),   allocatable :: dec_light
		real,              dimension(:),   allocatable :: cla
		real,              dimension(:),   allocatable :: nutrient
		real,              dimension(:,:), allocatable :: cells
		integer,           dimension(:),   allocatable :: mature
		real,              dimension(8)                :: d_type
		real                                           :: seedling_number
		real                                           :: wind_cat
		real                                           :: act_evap_day
		real                                           :: wilt_days
		real                                           :: flood_days, wl_days
		real                                           :: dry_days_upper_layer
		real                                           :: dry_days_base_layer
		real                                           :: amlt, leaf_area_w0
		real                                           :: aridity
		real                                           :: aridity_base
		integer                                        :: numspecies
		integer                                        :: numtrees, num_dead
		integer                                        :: fire, wind
		integer                                        :: windCount
		integer                                        :: stand_age
	end type PlotData

contains

	!:.........................................................................:

	subroutine initialize_plot(self, species, maxcells, maxheight, a_sat,      &
								a_fc, o_depth, i_text)
		!initializes invidual plots at a site
		!Author: Katherine Holcomb 2012, v. 1.0
		!Inputs/Outputs:
		!	self:       plot instance
		!	species:    species type instance
		!Inputs:
		!	maxcells:   maximum cells on plot (where maxtrees = maxcells*maxcells)
		!	maxheight:  maximum tree height (m)
		!	a_sat:      mineral (A-layer) saturation capacity (volumetric)
		!	a_fc: 		mineral (A-layer) field capacity (volumetric)
		!	o_depth: 	soil organic layer depth (m)
		!	i_text: 	soil texture (1: coarse-grained; 2: fine-grained)


		class(PlotData),                 intent(inout) :: self
		type(SpeciesData), dimension(:), intent(inout) :: species
		real,                            intent(in)    :: a_sat, a_fc
		real,                            intent(in)    :: o_depth
		integer,                         intent(in)    :: i_text
		integer,                         intent(in)    :: maxcells
		integer,                         intent(in)    :: maxheight

		real                                           :: burn
		real                                           :: duff
		real, parameter                                :: biowd = 0.0107
		real, parameter                                :: bioff = 0.01617
		real, parameter                                :: pctn = 0.008
		real, parameter                                :: bcon = 0.50
		real, parameter                                :: hcon = 0.658
		integer                                        :: n, i

		if (maxcells .ne. 0) then
			allocate(self%trees(maxcells*maxcells))
			allocate(self%deadtrees(maxcells*maxcells))
			allocate(self%cells(maxcells, maxcells))
		else
			stop "Must allow at least a few trees"
		endif

		if (maxheight .ne. 0) then
			allocate(self%con_light(maxheight))
			allocate(self%dec_light(maxheight))
			allocate(self%cla(maxheight))
		else
			stop "Must have a nonzero maximum height"
		endif

		!for now the possible species in each plot is assumed to be all in
		  !rangelist for the site. But they may not all be able to grow
		self%numspecies = size(species)
		self%seedling_number = 0.0

		allocate(self%seedling(self%numspecies))
		allocate(self%seedbank(self%numspecies))
		allocate(self%avail_spec(self%numspecies))
		allocate(self%nutrient(self%numspecies))
		allocate(self%species(self%numspecies))
		allocate(self%mature(self%numspecies))

		!copy the species list for the plot from the site's list.
		do n = 1, self%numspecies
			self%species(n) = species(n)
		enddo

		self%cells(:,:) = .false.
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
		self%wl_days = 0.0
		self%wilt_days = 0.0
		self%dry_days_upper_layer = 0.0
		self%dry_days_base_layer = 0.0
		self%amlt = 0.0

		self%seedbank = 0.0
		self%seedling = 0.0
		self%seedling_number = 0
		self%mature = 0

		self%avail_spec = 0.0
		self%dec_light = 0.0
		self%con_light =0.0
		self%d_type = 0.0
		self%nutrient = 1.0
		self%soil%lai_w0 = 1.0*0.15

		!initialize soil properties

		!snowpack initialized to 0.0
		self%soil%snowpack = 0.0
		self%soil%swe = 0.0

		!permafrost and litter initialization
		do i = 1, 20
			self%soil%litter(i) = 0.0
		end do

		self%soil%cohorts = 0.0
		self%soil%ncohort = 0


		!initializing from bare ground - assume secondary succession conditions
		  !wood litter layer (fresh wood)

		burn = (biowd/2.0)*bcon !amount burned in tonnes

		self%soil%litter(15) = (biowd/2.0 - burn)/0.0001 !small boles (t/ha)
		self%soil%litter(14) = (biowd/2.0 - burn)/0.0001 !large boles (t/ha)

		!humus layer
		self%soil%cohorts(1, 1) = (bioff*(1 - hcon))/0.0001 !t/ha
		self%soil%cohorts(1, 2) = self%soil%cohorts(1, 1)*pctn !tN/ha
		self%soil%cohorts(1, 5) = 20.0
		self%soil%ncohort = 1

		!initialize forest fuels to 0
		self%soil%avail_fuel = 0.0
		self%soil%fol_fuel = 0.0
		self%soil%twig_fuel = 0.0
		self%soil%smbl_fuel = 0.0
		self%soil%lrbl_fuel = 0.0
		self%soil%fan = 0.0

		!calculate organic layer depth
		duff = self%soil%cohorts(1, 1)/10000*plotsize*1000 !kg/plot
		self%soil%O_depth = (1.0/plotsize)*(duff/hum_bulk)

		!initialize A-layer, active, and moss-layer depths
		self%soil%A_depth = 1.0
		self%soil%M_depth = 0.0
		self%soil%moss_biom = 0.0
		self%soil%active = self%soil%A_depth

		!default values
		self%soil%O_bulk_dens = 84.0
		self%soil%A_bulk_dens = 1250.0

		self%soil%ASAT(1) = 0.39
		self%soil%AFC(1) = 0.39
		self%soil%APWP(1) = 0.039
		self%soil%APWP(2) = 0.06

		!values from site file
		self%soil%ASAT(2) = a_sat
		self%soil%AFC(2) = a_fc
		self%soil%itxt = i_text

		!from sitelist file
		if (o_depth .ne. rnvalid) self%soil%O_depth = o_depth

	end subroutine initialize_plot

	!:.........................................................................:

	subroutine delete_plot(self)
		!deletes a plot and its attributes
		!Author: Katherine Holcomb 2012, v. 1.0
		!Inputs/Outputs:
		!	self:       plot instance

		class(PlotData), intent(inout) :: self

		if (allocated(self%trees)) deallocate(self%trees)
		if (allocated(self%deadtrees)) deallocate(self%deadtrees)
		if (allocated(self%cells)) deallocate(self%cells)
		if (allocated(self%species)) deallocate(self%species)
		if (allocated(self%avail_spec)) deallocate(self%avail_spec)
		if (allocated(self%seedbank)) deallocate(self%seedbank)
		if (allocated(self%seedling)) deallocate(self%seedling)
		if (allocated(self%con_light)) deallocate(self%con_light)
		if (allocated(self%dec_light)) deallocate(self%dec_light)
		if (allocated(self%cla)) deallocate(self%cla)
		if (allocated(self%nutrient)) deallocate(self%nutrient)
		if (allocated(self%mature)) deallocate(self%mature)

	end subroutine delete_plot

	!:.........................................................................:

	subroutine tree_dm_cats(self, genera, field, diam_categories)
		!bins number of trees on a plot into DBH size classes, for writing to
		  !output files
		!Author: Katherine Holcomb 2012, v. 1.0
		!Inputs:
		!	self:           plot instance
		!	genera:         species or genus names
		!	field:          'species' or 'genus'
		!Outputs:
		!	diam_categories: number of stems in each DBH bin

		class(PlotData),                         intent(in)  :: self
		character(len = MAX_NLEN), dimension(:), intent(in)  :: genera
		character(len = 8),                      intent(in)  :: field
		integer, dimension(size(genera), NHC),   intent(out) :: diam_categories

		character(len = MAX_NLEN)                            :: comp
		integer, dimension(NHC)                              :: tree_diams
		integer                                              :: numitems
		integer                                              :: ig, n

		numitems = size(genera)
		diam_categories = 0

		do ig = 1, numitems
			do n = 1, self%numtrees
				if (field .eq. 'genus') then
					comp = self%trees(n)%genus_name
				else if (field .eq. 'species') then
					comp = self%trees(n)%unique_id
				endif

				if (comp == genera(ig)) then
					call get_diam_category(self%trees(n), tree_diams)
					diam_categories(ig, :) = diam_categories(ig, : ) +         &
						tree_diams
				endif
			enddo
		enddo

	end subroutine tree_dm_cats

	!:.........................................................................:

	subroutine sum_over_sg(self, genera, field, basal_area, leaf_bm, biomC,    &
							biomN, max_ht, max_diam, mean_dbh, mean_year,      &
							biom_categories, envresp_mn)
		!aggregates live species- or genus-level data for plot
		!Author: Katherine Holcomb 2012, v. 1.0
		!Inputs:
		!	self:            plot instance
		!	genera:          species or genus names
		!	field:           'species' or 'genus'
		!Outputs:
		!	basal_area:      basal area on plot (cm2)
		!   leaf_bm:         leaf biomass (tC)
		!   mean_dbh:        average DBH (cm)
		!   mean_year:       average tree age (years)
		!   bioMC:           biomass (tC)
		!   biomN:           N content (tN)
		!   max_ht:          maximum height (m)
		!   max_diam:        maximum DBH (cm)
		!   biom_categories: biomass binned by DBH (tC)
		!   envresp_mn:      mean growth response for each limiting factor (0-1)

		class(PlotData),                         intent(in)  :: self
		character(len = MAX_NLEN), dimension(:), intent(in)  :: genera
		character(len = 8),                      intent(in)  :: field
		real, dimension(size(genera)),           intent(out) :: basal_area,    &
													            leaf_bm,       &
													            mean_dbh,      &
													            mean_year
		real, dimension(size(genera)),          intent(out)  :: biomC
		real, dimension(size(genera)),          intent(out)  :: biomN
		real, dimension(size(genera)),          intent(out)  :: max_ht
		real, dimension(size(genera)),          intent(out)  :: max_diam
		real, dimension(size(genera), NHC),     intent(out)  :: biom_categories
		real, dimension(size(genera), 6),       intent(out)  :: envresp_mn

		integer, dimension(size(genera))                     :: n
		real,    dimension(size(genera))                     :: sum_dbh
		real,    dimension(size(genera))                     :: sum_year
		real,    dimension(size(genera), 6)                  :: envresp, n_er
		character(len=MAX_NLEN)                              :: comp
		real                                                 :: lf_biom
		real                                                 :: l_cn
		real                                                 :: tot_tree_bc
		real                                                 :: tot_tree_bn
		real                                                 :: dm
		integer                                              :: numitems
		integer                                              :: it, is

		numitems = size(genera)
		basal_area = 0.0
		biomC = 0.0
		biomN = 0.0
		leaf_bm = 0.0
		max_ht = 0.0
		max_diam = 0.0
		n = 0
		n_er = 0
		sum_dbh = 0.0
		sum_year = 0.0
		biom_categories = 0.0
		envresp = 0.0

		do is = 1, numitems
			do it = 1, self%numtrees

				if (field .eq. 'genus') then
					comp = self%trees(it)%genus_name
				else if (field .eq. 'species') then
					comp = self%trees(it)%unique_id
				endif

				if (comp == genera(is)) then

					!number of trees in that genera
					n(is) = n(is) + 1

					!get C:N ratio
					if (self%trees(it)%conifer) then
						l_cn = con_leaf_c_n
					else
						l_cn = dec_leaf_c_n
					endif

					!grab diameter
					dm = self%trees(it)%diam_bht

					!calculate and sum up basal area
					basal_area(is) = basal_area(is) + 0.25*pi*dm**2

					!calculate and sum biomass (tC) and N content
					lf_biom = self%trees(it)%leaf_bm
					leaf_bm(is) = leaf_bm(is) + lf_biom

					tot_tree_bc = self%trees(it)%biomC + lf_biom
					tot_tree_bn = self%trees(it)%biomN + lf_biom/l_cn

					biomC(is) = biomC(is) + tot_tree_bc
					biomN(is) = biomN(is) + tot_tree_bn

					!get biomass in DBH bins
					if (dm >= 0.5 .and. dm <= 5.0) then
						biom_categories(is, 1) =                               &
							biom_categories(is, 1) + tot_tree_bc
					else if (dm > 5.0 .and. dm <= 10.0) then
						biom_categories(is, 2) =                               &
							biom_categories(is, 2) + tot_tree_bc
					else if (dm > 10.0 .and. dm <= 20.0) then
						biom_categories(is, 3) =                               &
							biom_categories(is, 3) + tot_tree_bc
					else if (dm > 20.0 .and. dm <= 30.0) then
						biom_categories(is, 4) =                               &
							biom_categories(is, 4) + tot_tree_bc
					else if (dm > 30.0 .and. dm <= 40.0) then
						biom_categories(is, 5) =                               &
							biom_categories(is, 5) + tot_tree_bc
					else if (dm > 40.0 .and. dm <= 50.0) then
						biom_categories(is, 6) =                               &
							biom_categories(is, 6) + tot_tree_bc
					else if (dm > 50.0 .and. dm <= 60.0) then
						biom_categories(is, 7) =                               &
							biom_categories(is, 7) + tot_tree_bc
					else if (dm > 60.0 .and. dm <= 70.0) then
						biom_categories(is, 8) =                               &
							biom_categories(is, 8) + tot_tree_bc
					else if (dm > 70.0 .and. dm <= 80.0) then
						biom_categories(is, 9) =                               &
							biom_categories(is, 9) + tot_tree_bc
					else if (dm > 80.0 .and. dm <= 90.0) then
						biom_categories(is, 10) =                              &
							biom_categories(is, 10) + tot_tree_bc
					else if (dm > 90.0) then
						biom_categories(is, 11) =                              &
							biom_categories(is, 11) + tot_tree_bc
					endif

					!get maximum height and DBH
					max_ht(is) = max(max_ht(is), self%trees(it)%forska_ht)
					max_diam(is) = max(max_diam(is), self%trees(it)%diam_bht)

					!sum up DBH and age
					sum_dbh(is) = sum_dbh(is) + self%trees(it)%diam_bht
					sum_year(is) = sum_year(is) + float(self%trees(it)%tree_age)

					!sum up the growth responses for each factor
					if (self%trees(it)%env_resp(1) .ge. 0.0) then
						n_er(is, 1) = n_er(is, 1) + 1
						envresp(is, 1) = envresp(is, 1) +                      &
							self%trees(it)%env_resp(1)
					end if

					if (self%trees(it)%env_resp(2) .ge. 0.0) then
						n_er(is, 2) = n_er(is, 2) + 1
						envresp(is, 2) = envresp(is, 2) +                      &
							self%trees(it)%env_resp(2)
					end if
					if (self%trees(it)%env_resp(3) .ge. 0.0) then
						n_er(is, 3) = n_er(is, 3) + 1
						envresp(is, 3) = envresp(is, 3) +                      &
							self%trees(it)%env_resp(3)
					end if
					if (self%trees(it)%env_resp(4) .ge. 0.0) then
						n_er(is, 4) = n_er(is, 4) + 1
						envresp(is, 4) = envresp(is, 4) +                      &
							self%trees(it)%env_resp(4)
					end if
					if (self%trees(it)%env_resp(5) .ge. 0.0) then
						n_er(is, 5) = n_er(is, 5) + 1
						envresp(is, 5) = envresp(is, 5) +                      &
							self%trees(it)%env_resp(5)
					end if
					if (self%trees(it)%env_resp(6) .ge. 0.0) then
						n_er(is, 6) = n_er(is, 6) + 1
						envresp(is, 6) = envresp(is, 6) +                      &
							self%trees(it)%env_resp(6)
					end if
				endif
			enddo
		enddo

		!get average of DBH, age, and environmental responses
		do is = 1, numitems
			if (n(is) .eq. 0) then
				mean_dbh(is) = 0.0
				mean_year(is) = 0.0
				envresp_mn(is, :) = rnvalid
			else
				mean_dbh(is) = sum_dbh(is)/n(is)
				mean_year(is) = sum_year(is)/n(is)
				envresp_mn(is, 1) = envresp(is, 1)/n_er(is, 1)
				envresp_mn(is, 2) = envresp(is, 2)/n_er(is, 2)
				envresp_mn(is, 3) = envresp(is, 3)/n_er(is, 3)
				envresp_mn(is, 4) = envresp(is, 4)/n_er(is, 4)
				envresp_mn(is, 5) = envresp(is, 5)/n_er(is, 5)
				envresp_mn(is, 6) = envresp(is, 6)/n_er(is, 6)

			endif
		enddo

	end subroutine sum_over_sg

	!:.........................................................................:

	subroutine sum_over_deadsg(self, genera, field, biomC, mean_dbh,           &
								death_markers)
		!aggregates dead species- or genus-level data for plot
		!Authors: Adrianna Foster and Jacquelyn Shuman 2015, v. 1.0
		!Inputs:
		!	self:          plot instance
		!	genera:        species or genus names
		!	field:         'species' or 'genus'
		!Outputs:
		!	biomC:         biomass (tC)
		!   mean_dbh:      average DBH (cm)
		!   death_markers: biomass mortality from each limiting factor or
		!					disturbance (tC)


		class(PlotData),                         intent(in)  :: self
		character(len = MAX_NLEN), dimension(:), intent(in)  :: genera
		character(len = 8),                      intent(in)  :: field
		real, dimension(size(genera)),           intent(out) :: mean_dbh
		real, dimension(size(genera)),           intent(out) :: biomC
		real, dimension(size(genera), 8),        intent(out) :: death_markers

		integer, dimension(size(genera))                     :: n
		real,    dimension(size(genera))                     :: sum_dbh
		character(len=MAX_NLEN)                              :: comp
		real                                                 :: lf_biom
		real                                                 :: l_cn
		real                                                 :: tot_tree_bc
		real                                                 :: dm
		integer                                              :: numitems
		integer                                              :: it, is
		integer                                              :: ks

		numitems = size(genera)
		biomC = 0.0
		n = 0
		sum_dbh = 0.0
		death_markers = 0.0

		do is = 1, numitems
			do it = 1, self%num_dead

				if (field .eq. 'genus') then
					comp = self%deadtrees(it)%genus_name
				else if (field .eq. 'species') then
					comp = self%deadtrees(it)%unique_id
				endif

				if (comp == genera(is)) then

					!get number of trees in that genera
					n(is) = n(is) + 1

					!get C:N ratio
					if (self%deadtrees(it)%conifer) then
						l_cn = con_leaf_c_n
					else
						l_cn = dec_leaf_c_n
					endif

					!what was the killing factor
					ks = self%deadtrees(it)%stressor

					!get diameter
					dm = self%deadtrees(it)%diam_bht

					!get and sum up biomass C

					lf_biom = self%deadtrees(it)%leaf_bm
					tot_tree_bc = self%deadtrees(it)%biomC + lf_biom
					biomC(is) = biomC(is) + tot_tree_bc

					!bin biomass into the limiting factors/disturbances
					death_markers(is, ks) = death_markers(is, ks) + tot_tree_bc

					!sum up DBH
					sum_dbh(is) = sum_dbh(is) + self%deadtrees(it)%diam_bht

				endif
			end do
		end do

		!get average DBH
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
