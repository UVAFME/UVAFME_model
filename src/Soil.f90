module Soil

use Constants
use Parameters
use Random

implicit none

!*******************************************************************************
  !
  !The Soil module contains attributes and procedures relevant to soil
  !  properties.
  !Author: Adrianna Foster, 2018 v. 1.0
  !
  !Methods:
  ! permf:       adapted from Bonan (1989), added by AC Foster 5/15/2017
  ! moist:       adapted from Bonan (1989), added by AC Foster 5/16/2017
  ! soiln:       adapted from Bonan (1990) and Pastor & Post (1985), added by
  !					AC Foster 10/10/2017
  ! moss:        adapted from Bonan & Korzukin (1989), added by AC Foster
  !  				10/10/2017
  ! forest_fuels: adapted from Bonan (1990) and Schumacher et al. (2007), added
  !					by AC Foster 10/10/2017
  !
!*******************************************************************************

	real       :: lai_min = 0.01
	real       :: lai_max = 0.15

	type SoilData
		real                      :: lai_w0
		real                      :: N_used
		real                      :: avail_N
		real                      :: net_prim_prodN
		real                      :: net_prim_prodC
		real                      :: runoff
		real                      :: total_C_rsp

		!permafrost and water
		real, dimension(2)        :: z_drain
		real, dimension(2)        :: ASAT, AFC, APWP
		real, dimension(2)        :: wc, H2Oice, water, owater
		real, dimension(2)        :: d_melt, d_freeze, od_melt, od_freeze
		real, dimension(2)        :: csm, pwp, sat, fc, wt, wf, xs
		real                      :: active, z_freeze
		real                      :: A_depth, O_depth
		real                      :: O_bulk_dens, A_bulk_dens
		real                      :: snowpack, swe
		real                      :: minWC

		!litter and moss
		real, dimension(20)       :: litter
		real, dimension(20, 3)    :: forest_litter
		real, dimension(1500, 15) :: cohorts
		real                      :: leaf_litterN
		real                      :: moss_biom, M_depth, dec_fuel

		!fuels tracking
		real                      :: avail_fuel, fol_fuel, twig_fuel
		real                      :: smbl_fuel, lrbl_fuel, fan

		integer                   :: itxt, ncohort

	end type SoilData


  contains

	!:.........................................................................:

	subroutine soiln(soil, aet_mm, cla, soildays, flooddays, avail_n)
		!calculates annual decomposition of litter cohorts and calculates
		  !palnt-available N
		!Author: Adrianna Foster 2018, v. 1.0
		!Inputs/Outputs:
		!	soil:   soil instance
		!Inputs:
		!	aet_mm:    actual evapotranspiration (mm/year)
		!	cla:       cumulative leaf area
		!	drI:       drought index
		!	soildays:  degree-days above 0degC
		!	flooddays: flood days
		!Outputs:
		!	avail_n:   plant-available N (tN/ha)

		class(SoilData),   intent(inout) :: soil
		real,              intent(in)    :: aet_mm, cla
		real,              intent(in)    :: soildays, flooddays
		real,              intent(out)   :: avail_n

		real, dimension(1500, 15)        :: C
		real, dimension(19)              :: pwl
		real                             :: litCO2, humCO2, totCO2
		real                             :: WDW_new
		real                             :: tot_Nimob, smult
		real                             :: lit_Nmin, hum_Nmin, tot_Nmin
		real                             :: xaet, aetm
		real                             :: leaf_litterC, can_prod
		real                             :: fc, dry, dec_lit
		real                             :: lit_decaym, perc_wtloss, lt
		real                             :: cla_decaym
		real                             :: tfallNmin
		real                             :: perc_rem, wtloss
		real                             :: delta_N, crit_wtloss
		real                             :: hum_Nnew, hum_OMnew, humCN
		real                             :: floor_litter
		real                             :: moss_litter
		real                             :: con_fuel, dec_fuel, twig_fuel
		real                             :: duff, min_Nmin, perc_wtloss2
		real                             :: ddmult
		real                             :: UPN, mdepth
		integer                          :: i, j, nc, ix, ilt

		!initialize accumulators
		litCO2 = 0.0
		humCO2 = 0.0
		totCO2 = 0.0
		WDW_new = 0.0
		tot_Nimob = 0.0
		avail_n = 0.0
		lit_Nmin = 0.0
		hum_Nmin = 0.0
		tot_Nmin = 0.0
		soil%leaf_litterN = 0.0
		leaf_litterC = 0.0

		do i = 1, 20
			do j = 1, 3
				soil%forest_litter(i, j) = 0.0
			end do
		end do

		!reduce OOP table lookups
		do i = 1, 1500
			do j = 1, 15
				C(i, j) = soil%cohorts(i, j)
			end do
		end do
		nc = soil%ncohort
		fc = (soil%fc(1) + soil%fc(2))*100.0
		dry = (soil%pwp(1) + soil%pwp(2))*100.0

		!calculate leaf litter N
		do i = 1, 12
			soil%leaf_litterN = soil%leaf_litterN +                            &
				(soil%litter(i)*litter_params(i, 2))
		end do

		!calculate AET multiplier - from Pastor & Post (1985)
		xaet = aet_mm
		if (xaet .gt. 600.0) xaet = 600.0
		aetm = (-1.0*xaet)/(-1200.0 + xaet)

		ddmult = 2**(0.005*(soildays - 1950))

		mdepth = (soil%M_depth +                                               &
			((soil%forest_litter(18,2)/10000*plotsize*1000)*(1/plotsize)/      &
			moss_bulk))*100.0

		smult = (max((1.6 - sqrt(mdepth*3)), 0.0) +                            &
			max(3*(1.0-flooddays/0.8)**2, 0.2))/2

		!populate cohort array with new cohorts and cohort properties
		do i = 1, 16
			if (soil%litter(i) .gt. epsilon(1.0)) then
				nc = nc + 1
				if (nc .gt. 1500) print *, 'ERROR: ncohort > 1500'
				C(nc, 1) = soil%litter(i)*litter_params(i, 10) !current weight (t/ha)
				C(nc, 2) = soil%litter(i)*litter_params(i, 2) !current N (tN/ha)
				do j = 3, 9
					C(nc, j) = litter_params(i, j) !litter parameters (from input file)
				end do
				C(nc, 10) = soil%litter(i)*litter_params(i, 10) !initial weight (t/ha)
				C(nc, 11) = litter_params(i, 2) !current N %
				C(nc, 12) = litter_params(i, 7)*1.7039 + 0.0955 !critical weight (t/ha)
				if (C(nc, 5) .eq. 14.0) C(nc, 12) = 0.30
				if (C(nc, 5) .eq. 15.0) C(nc, 12) = 0.30
				if (C(nc, 5) .eq. 16.0) C(nc, 12) = 0.30
				C(nc, 13) = 0.0 !year
			end if
		end do

		!deal with moss litter
		if (soil%litter(18) .gt. epsilon(1.0)) then
			nc = nc + 1
			C(nc, 1) = soil%litter(18)*litter_params(18, 10)
			C(nc, 2) = soil%litter(18)*litter_params(18, 2)
			do j = 3, 9
				C(nc, j) = litter_params(18, j)
			end do
			C(nc, 10) = soil%litter(18)*litter_params(18, 10)
			C(nc, 11) = litter_params(18, 2)
			C(nc, 12) = 0.0
			C(nc, 13) = 0.0
		end if

		!calculate litter decay multiplier, simulate effects of gaps on decay
		  !from Pastor & Post (1985)
		do i = 1, 12
			leaf_litterC = leaf_litterC + soil%litter(i)
		end do
		can_prod = 1.54 + 0.0457*(fc - dry)
		if (leaf_litterC .gt. can_prod) leaf_litterC = can_prod
		lit_decaym = 1.0 + (-0.5 + 0.075*(fc - dry))*(1.0 - leaf_litterC/can_prod)

		!calculate bonan's decay multiplier - based on LAI
		  !Bonan (1990)
		if (cla/plotsize .ge. 2.5) then
			cla_decaym = 1.0
		else
			cla_decaym = 1.0 + 1.5*sqrt(1.0 - (cla/plotsize)/2.5)
		end if

		!bypass litter cohort calculations if there is nothing to
		  !decay
		if (nc .gt. 1) then

			!loop to calculate litter decay, N immobilization, lignin
			  !decay, and litter CO2 evolution
			do i = 2, nc
				!calculate percent weight loss based on AET and lignin:N ratio
				!note: can't use this equation for moss litter decay
				perc_wtloss = (0.9804 + 0.09352*aet_mm) -                      &
					((-0.4956 + 0.00193*aet_mm)*(C(i, 7)/C(i, 11)))
				perc_wtloss = (perc_wtloss*lit_decaym)/100.0
				if (perc_wtloss .gt. 0.99) perc_wtloss = 0.99

				!bonan's equation uses N:C ratio and active layer thickness
				  !note: CAN use this for moss litter decay
				perc_wtloss2 = (-0.0052 + 2.08*(C(i,2)/C(i, 1)))*              &
					exp(0.898*soil%active)
				perc_wtloss2 = max(perc_wtloss2*cla_decaym, 0.0)
				if (perc_wtloss2 .gt. 0.99) perc_wtloss2 = 0.99

				if (soil%active .le. 1.5) then
					!average these two percent weight losses when
					  !permafrost present
					perc_wtloss = ((perc_wtloss2 + perc_wtloss)/2.0)*ddmult*   &
						smult
					perc_wtloss = min(max(perc_wtloss, 0.01), 0.99)
				else
					perc_wtloss = min(max(perc_wtloss*ddmult*smult, 0.01), 0.99)
				endif

				!litter cohort type
				lt = C(i, 5)

				!weight loss of moss litter can only use Bonan's eq.
				if (lt .eq. 18.0) perc_wtloss = min(max(perc_wtloss2*ddmult*   &
					smult, 0.01), 0.99)

				!max weight loss of large wood (DBH > 10cm) is 3%
				if (lt .eq. 15.0) perc_wtloss = 0.03

				!max weight loss of small wood is 10%
				if (lt .eq. 14.0) perc_wtloss = 0.1

				!maximum weight loss of well-decayed wood is 5%
				if (lt .eq. 17.0) perc_wtloss = 0.05

				!weight loss of twigs is less than 20%
				if (lt .eq. 16.0) perc_wtloss = min(perc_wtloss, 0.2)

				pwl(lt) = perc_wtloss

				!calculate actual weight loss (tonnes/ha)
				wtloss = perc_wtloss*C(i, 1)

				!calculate weight loss for litter other than moss
				if (lt .ne. 18.0) then

					!calculate fraction of organic matter remaining following
					  !this weight loss
					perc_rem = (C(i, 1) - wtloss)/(C(i, 10))

					!calculate new N concentration
					C(i, 11) = min(max(C(i, 3) - C(i, 4)*perc_rem, 0.0001),    &
						0.99)

					!retain cohort another year if fraction remaining is
					  !greater than fraction which will become humus/WDW

					if (perc_rem .le. C(i, 12)) then !gets transferred

						!recalculate weight loss and N concentration
						wtloss = max(C(i, 1) - C(i, 12)*C(i, 10), 0.0)
						C(i, 11) = min(max(C(i, 3)- C(i, 4)*C(i, 12),          &
							0.0001), 0.99)

						!calculate absolute change in N content
						delta_N = C(i, 2) - C(i, 11)*(C(i, 1) - wtloss)

						!if positive - immobilize N
						!if negative - mineralize N
						if (delta_N .lt. epsilon(1.0)) then
							tot_Nimob = tot_Nimob  - delta_N
						else
							lit_Nmin = lit_Nmin + delta_N
						end if

						!tranfer cohorts to humus or well-decayed wood
						if (C(i, 6) .eq. 1.0) then
							C(1, 1) = C(1, 1) + C(i, 1) - wtloss
							C(1, 2) = C(1, 2) + C(i, 11)*(C(i, 1) - wtloss)
							C(i, 1) = 0.0
						else
							WDW_new = WDW_new + C(i, 1) - wtloss
							C(i, 1) = 0.0
						end if

					end if

				else if (lt .eq. 18.0) then
					!decay moss litter

					crit_wtloss = (C(i, 4)*C(i, 1) - C(i, 2))/                 &
						(C(i, 3) + C(i, 4))

					if (wtloss .ge. crit_wtloss) then
						!transfer cohort to humus, weight loss proceeds
						 !at critical weight loss rate

						wtloss = crit_wtloss

						!transfer cohort
						UPN = wtloss*C(i,3)
						C(i,11) = C(i, 4)
						C(1, 1) = C(1, 1) + C(i, 1) - wtloss
						C(1, 2) = C(1, 2) + UPN

						!calculate absolute change in N content
						delta_N = C(i, 2) - C(i, 11)*(C(i, 1) - wtloss)
						if (delta_N .lt. epsilon(1.0)) then
							tot_Nimob = tot_Nimob - delta_N
						else
							lit_Nmin = lit_Nmin + delta_N
						end if

						C(i, 1) = 0.0
					else
						!cohort is not transferred
						UPN = wtloss*C(i,3)
						C(i,2) = C(i,2) + UPN
						C(i,11) = C(i,2)/C(i,1)
					end if
				end if

				!update weight and N content of cohorts that didn't get
				  !transferred to humus/WDW
				if (C(i, 1) .gt. epsilon(1.0)) then
					C(i, 1) = C(i, 1) - wtloss

					if (C(i, 5) .ne. 18.0) then
						C(i, 2) = C(i, 1)*C(i, 11)
						C(i, 7) = C(i, 8) - C(i, 9)*(C(i, 1)/C(i, 10))
					else
						C(i, 7) = 0.0
					end if
					C(i, 13) = C(i, 13) + 1.0
				end if

				!calculate litter cohort CO2 evolution
				litCO2 = litCO2 + (wtloss*0.48)

			end do !end litter cohort loop

			!throughfall is 16% of leaf litter N
			!tot_Nimob = tot_Nimob - 0.16*soil%leaf_litterN
			tfallNmin = 0.16*soil%leaf_litterN

		end if !end if ncohort > 1

		!calculate humus N mineralization
		perc_wtloss = (-0.0052 + 2.08*(C(1, 2)/C(1, 1)))*exp(0.898*soil%active)
		perc_wtloss = perc_wtloss*cla_decaym*ddmult*smult
		perc_wtloss = min(perc_wtloss, 1.0)
		perc_wtloss = max(perc_wtloss, 0.0)

		if (soil%active .le. 1.5) then !0.9
			hum_Nmin = C(1, 2)*perc_wtloss
			pwl(19) = perc_wtloss
		else
			perc_wtloss = max(0.0001, min(0.035*lit_decaym*ddmult*smult, 0.99))
			hum_Nmin = C(1, 2)*perc_wtloss
			pwl(19) = hum_Nmin/C(1, 2)
		end if

		!subtract mineralized N from humus N and calculate humus CO2
		hum_Nnew = C(1, 2) - hum_Nmin
		hum_OMnew = C(1,1)*(hum_Nnew/C(1,2))

		humCO2 = (C(1,1) - hum_OMnew)*0.48
		C(1, 1) = hum_OMnew
		C(1, 2) = hum_Nnew

		!humus C:N ratio
		humCN = (0.48*C(1, 1))/C(1, 2)

		!add humus N mineralization to cohort N mineralization to get
		  !total N mineralization
		tot_Nmin = lit_Nmin + hum_Nmin

		min_Nmin = 10.0*min(soil%active, soil%A_depth)*0.001

		!subtract immobilization from total mineralization to get
		  !plant-available N
		!add throughfall N mineralization, mineral N, and N from fires
		avail_n = max(0.0, tot_Nmin - tot_Nimob + soil%fan + tfallNmin +       &
			min_Nmin)

		!calculate total soil respiration
		totCO2 = litCO2 + humCO2

		!remove transferred cohorts
		ix = 0
		do 20 i = 1, nc
			if (C(i, 1) .eq. 0.0) GO TO 16
				do 12 j = 1, 13
12                C(i-ix,j) = C(i, j)
			GO TO 20
16          ix = ix + 1
20       CONTINUE
		nc = nc - ix

		!create new well-decayed wood cohort
		if (WDW_new .gt. epsilon(1.0)) then
			nc = nc + 1
			if (nc .gt. 1500) print *, 'ERROR, ncohort > 1500'
			C(nc, 1) = WDW_new
			C(nc, 2) = WDW_new*litter_params(17, 2)
			do j = 3, 9
				C(nc, j) = litter_params(17, j)
			end do
			C(nc, 10) = WDW_new
			C(nc, 11) = litter_params(17, 2)
			C(nc, 12) = 0.50
		end if

		!calculate total weight and N content by forest floor
		!compartment
		do i = 1, nc
			ilt = int(C(i, 5))
			soil%forest_litter(ilt, 1) = C(i, 5)
			soil%forest_litter(ilt, 2) = soil%forest_litter(ilt, 2) + C(i, 1)
			soil%forest_litter(ilt, 3) = soil%forest_litter(ilt, 3) + C(i, 2)
		end do

		!calculate organic layer depth
		con_fuel = 0.0
		dec_fuel = 0.0
		dec_lit = 0.0
		do i = 1, 12
			if (i .eq. 6 .or. i .ge. 10) then
				con_fuel = con_fuel + soil%forest_litter(i, 2)
			else
				dec_fuel = dec_fuel + soil%forest_litter(i, 2)
			end if
		end do

		do i = 1, 12
			if (i .ne. 6 .or. i .lt. 10) then
				dec_lit = dec_lit + soil%litter(i)
			end if
		end do

		soil%dec_fuel = dec_lit !save for moss subroutine


		con_fuel = con_fuel/10000*plotsize*1000 !kg/plot
		dec_fuel = dec_fuel/10000*plotsize*1000 !kg/plot
		twig_fuel = soil%forest_litter(16, 2)/10000*plotsize*1000
		moss_litter = soil%forest_litter(18, 2)/10000*plotsize*1000
		duff = C(1, 1)/10000*plotsize*1000

		floor_litter = 0.0
		do i = 1,12
			floor_litter = floor_litter + soil%forest_litter(i, 2)
		end do
		floor_litter = (floor_litter + soil%forest_litter(16, 2) +             &
						soil%forest_litter(18, 2))/10000*plotsize*1000
		duff = C(1, 1)/10000*plotsize*1000


		soil%O_depth = (1.0/plotsize)*((con_fuel/con_bulk) +                   &
			(dec_fuel/dec_bulk) + (twig_fuel/bulk_l) + (duff/bulk_duff) +      &
			(moss_litter/moss_bulk))

		!reassign attributes
		do i = 1, 1500
			do j = 1, 15
				soil%cohorts(i, j) = C(i, j)
			end do
		end do
		soil%ncohort = nc


		do i = 1, 18
			soil%litter(i) = 0.0
		end do

		return

	end subroutine soiln

	!:.........................................................................:

	subroutine forest_fuels(soil, drI, avail_fuel, consN, consRoot)
		!calculates available fuels for forest fires and burns up litter
		  !cohorts accordingly
		!Author: Adrianna Foster 2018 v. 1.0
		!Inputs/Outputs:
		!	soil:       soil instance
		!Inputs:
		!	drI:        drought index
		!Outputs:
		!	avail_fuel: available fuel for fire (t/ha)
		!   consN:      portion of N consumed by fire
		!   consRoot:   portion of live roots consumed by fire
		!
		class(SoilData),          intent(inout) :: soil
		real,                     intent(in)    :: drI
		real,                     intent(out)   :: avail_fuel, consN
		real,                     intent(out)   :: consRoot

		real, dimension(1500, 15)               :: C
		real                                    :: fol_fuel
		real                                    :: lt
		real                                    :: rfs, vwc, root_kill
		real                                    :: WDW_kill, hum_kill
		real                                    :: N_cons, M_loss, BO
		real                                    :: sfan
		real                                    :: oldN
		real                                    :: volN, humloss, oldHN
		real                                    :: twigloss, stemloss
		real                                    :: WDWloss, rootloss
		real                                    :: con_fuel, dec_fuel
		real                                    :: moss_fuel
		real                                    :: twig_fuel, duff
		real, parameter                         :: chcoal = 0.25
		integer                                 :: i, j, nc, ix, ilt

		!initialize accumulators
		sfan = 0.0

		!reduce table lookups
		nc = soil%ncohort
		do i = 1, 1500
			do j = 1, 15
				C(i, j) = soil%cohorts(i, j)
			end do
		end do

		!calculate available fuels

		!foliage - all leaf litter is available for burning
		fol_fuel = 0.0
		do i = 1, 12
			fol_fuel = fol_fuel + soil%forest_litter(i, 2)
		end do
		soil%fol_fuel = fol_fuel

		!twigs, small boles, and large boles
		  !equations for fuel amount are from Schumacher et al. (2006),
		  !Landscape Ecology
		soil%twig_fuel = (tfc1 + tfc2*drI)*(soil%forest_litter(16, 2))
		soil%smbl_fuel = (tfc1 + tfc2*drI)*(soil%forest_litter(14, 2))
		soil%lrbl_fuel = (bfc*drI)*soil%forest_litter(15, 2)

		soil%avail_fuel = fol_fuel + soil%twig_fuel + soil%smbl_fuel +         &
			soil%lrbl_fuel

		avail_fuel = soil%avail_fuel

		!calculate fire severity for humus, root, and well-decayed wood
		  !consumption
		vwc = soil%minWC*soil%A_bulk_dens/1000.0

		rfs = max(0.0, (soil%ASAT(2) - vwc)/(soil%ASAT(2) - soil%APWP(2)))

		root_kill = max(0.0, min((0.302 + 0.597*rfs +                          &
			3.34*(soil%M_depth + soil%O_depth)), 0.9))
		consRoot = root_kill

		WDW_kill = max(0.0, min((0.098 + 0.597*rfs +                           &
			3.34*(soil%M_depth + soil%O_depth)), 0.5))

		hum_kill = max(0.0, min((0.079 + 0.5744*rfs +                          &
			3.34*(soil%M_depth + soil%O_depth)), 0.658))

		N_cons = max(0.0, min((0.6426*rfs +                                    &
			3.34*(soil%M_depth + soil%O_depth)), 0.7))
		consN = N_cons

		!set % nitrogen in forest floor not volatalized
		volN = 1.0 - N_cons

		!save old value of moss_biom
		BO = soil%moss_biom

		!remove moss layer
		M_loss = min(soil%M_depth, hum_kill*soil%M_depth)
		soil%M_depth = max(soil%M_depth - M_loss, 0.0)
		soil%moss_biom = soil%M_depth*28.6*plotsize
		sfan = sfan + (BO - soil%moss_biom)*litter_params(18, 2)*volN

		!remove humus layer
		humloss = min(hum_kill*C(1, 1), C(1, 1))
		oldHN = C(1, 2)
		C(1, 2) = (1.0 - humloss/C(1,1))*C(1, 2)
		sfan = sfan + (oldHN - C(1, 2))*volN
		C(1, 1) = max(C(1, 1) - hum_kill*C(1, 1), 0.0)

		!get rid of cohorts
		do i = 2, nc

			lt = C(i, 5)

			!all leaf litter fuels decomposed
			if (lt .le. 12.0) then
				C(i, 1) = 0.0
				sfan = sfan + C(i, 2)*volN
			end if

			!roots
			if (lt .eq. 13.0) then
				rootloss = min(root_kill*C(i, 1), C(i, 1))
				C(i, 1) = max(C(i, 1) - root_kill*C(i, 1), 0.0)
				oldN = C(i, 2)
				C(i, 2) = C(i, 2)*(1.0 - root_kill)
				C(i, 11) = C(i, 2)/C(i,1)
				sfan = sfan + (oldN - C(i, 2))*volN
			end if

			!twigs and small trees
			if (lt .eq. 14.0 .or. lt .eq. 16.0) then
				twigloss = min((tfc1 + tfc2*drI)*C(i, 1), C(i, 1))
				C(i, 1) = max(C(i, 1) - (tfc1 + tfc2*drI)*C(i, 1), 0.0)
				oldN = C(i, 2)
				C(i, 2) = C(i, 2)*(1.0 - (tfc1 + tfc2*drI))
				C(i, 11) = C(i, 2)/C(i,1)
				sfan = sfan + (oldN - C(i, 2))*volN
			end if

			!boles
			if (lt .eq. 15.0) then
				stemloss = min((bfc*drI)*C(i, 1), C(i, 1))
				C(i, 1) = max(C(i, 1) - (bfc*drI)*C(i, 1), 0.0)
				oldN = C(i, 2)
				C(i, 2) = C(i, 2)*(1.0 - (bfc*drI))
				C(i, 11) = C(i, 2)/C(i,1)
				sfan = sfan + (oldN - C(i, 2))*volN
			end if

			!well-decayed wood
			if (lt .eq. 17.0) then
				WDWloss = min(WDW_kill*C(i, 1), C(i, 1))
				C(i, 1) = max(C(i, 1) - WDW_kill*C(i, 1), 0.0)
				oldN = C(i, 2)
				C(i, 2) = C(i, 2)*(1.0 - WDW_kill)
				C(i, 11) = C(i, 2)/C(i,1)
				sfan = sfan + (oldN - C(i, 2))*volN
			end if

			!all moss litter consumed
			if (lt .eq. 18.0) then
				C(i, 1) = 0.0
			end if

		end do

		!remove completely burned cohorts
		ix = 0
		do i = 1, nc
			if (C(i, 1) .ge. epsilon(1.0)) then
				do j = 1, 12
					C(i-ix, j) = C(i,j)
				end do
			end if
			ix = ix + 1
		end do

		!calculate total weight and N content by forest floor
		!compartment
		do i = 1, nc
			ilt = int(C(i, 5))
			soil%forest_litter(ilt, 1) = C(i, 5)
			soil%forest_litter(ilt, 2) = soil%forest_litter(ilt, 2) +  C(i, 1)
			soil%forest_litter(ilt, 3) = soil%forest_litter(ilt, 3) +  C(i, 2)
		end do

		soil%forest_litter(19, 1) = 19.0
		do ilt = 1, 12
			soil%forest_litter(19, 2) = soil%forest_litter(19, 2) +            &
				soil%forest_litter(ilt, 2)
			soil%forest_litter(19, 3) = soil%forest_litter(19, 3) +            &
				soil%forest_litter(ilt, 3)
		end do

		soil%forest_litter(19, 2) = soil%forest_litter(19, 2) +                &
			soil%forest_litter(18, 2) + soil%forest_litter(13, 2)
		soil%forest_litter(19, 3) = soil%forest_litter(19, 3) +                &
			soil%forest_litter(18, 3) + soil%forest_litter(13, 3)

		!reassign attributes
		do i = 1, 1500
			do j = 1, 15
				soil%cohorts(i, j) = C(i, j)
			end do
		end do
		soil%ncohort = nc
		soil%fan = sfan

		!calculate organic layer depth
		con_fuel = 0.0
		dec_fuel = 0.0
		do i = 1, 12
			if (i .eq. 6 .or. i .ge. 10) then
				con_fuel = con_fuel + soil%forest_litter(i, 2)
			else
				dec_fuel = dec_fuel + soil%forest_litter(i, 2)
			end if
		end do

		soil%dec_fuel = dec_fuel !save for moss subroutine


		con_fuel = con_fuel/10000*plotsize*1000
		dec_fuel = dec_fuel/10000*plotsize*1000
		twig_fuel = soil%forest_litter(16, 2)/10000*plotsize*1000
		moss_fuel = soil%forest_litter(18, 2)/10000*plotsize*1000
		duff = C(1, 1)/10000*plotsize*1000

		soil%O_depth = (1.0/plotsize)*((con_fuel/con_bulk) +                   &
			(dec_fuel/dec_bulk) + (twig_fuel/bulk_l) + (duff/bulk_duff) +      &
			(moss_fuel/con_bulk))

	    return

	end subroutine forest_fuels

	!:.........................................................................:

	subroutine permf(soil, month, knd, alff, tdd,  fdd, cfs, zdepth)
		!calculates depth of freeze/thaw in the soil profile (zdepth, m)
		!Author: Adrianna Foster 2018 v. 1.0
		!Inputs/Outputs:
		!	soil:   soil instance
		!Inputs:
		!	month:  simulation month
		!	alff:   available light at forest floor (0-1)
		!   tdd:    thawing degree-days
		!	fdd:    freezing degree-days
		!	cfs: 	fraction sunlight attenuated through atm.
		!Outputs
		!   zdepth: depth of freeze/thaw (m)

		class(SoilData),        intent(inout) :: soil
		real, dimension(12, 2), intent(in)    :: tdd, fdd
		real,                   intent(in)    :: alff, cfs
		integer,                intent(in)    :: month, knd
		real,                   intent(out)   :: zdepth

		real                                  :: tyr, cft, cff, cadd, dl
		real                                  :: vwc, tku, tkf, gwc, bde
		real                                  :: tk, ql, res, degd, depth
		real                                  :: ax, bx, cx
		real                                  :: Ad, Od, ABD, OBD
		integer                               :: l

		!convenience variables to reduce table lookups
		Ad = soil%A_depth
		Od = soil%O_depth
		ABD = soil%A_bulk_dens
		OBD = soil%O_bulk_dens

		!initialize zdepth to 0.0 and tyr to snowpack
		zdepth = 0.0
		tyr = 0.0

		if (soil%M_depth .gt. 0) then
			tyr = soil%M_depth/(0.291*0.86)
		end if

		!adjust freeze/thaw days for surface conditions
		if (alff .gt. 0.75) then
			cft = 0.92
			cff = 0.36
		else
			if (alff .le. 0.50) then
				cft = 0.62
				cff = 0.38
			else
				cft = 0.77
				cff = 0.37
			end if
		end if

		if (knd .eq. 1) then
			cadd = fdd(month, 2)*cff*(2.0 - cfs)
		else
			cadd = tdd(month, 2)*cft*cfs
		end if

		!calculate freeze/thaw soil layers
		do l = 1, 2
			if (l .eq. 1) then
				dl = Od
			else
				dl = Ad  + 2.0
			end if

			if (l .eq. 1) then
				!moss/humus layer conductivities
				vwc = soil%wc(l)*OBD/1000.0 !volumetric moisture content
				tku = (0.5*(soil%APWP(l) - vwc) + 0.08*(vwc - soil%AFC(l)))/   &
					(soil%APWP(l) - soil%AFC(l))
				tkf = (2.0*tku*(soil%APWP(l) - vwc) + tku*(vwc -               &
					soil%AFC(l)))/(soil%APWP(l) - soil%AFC(l))

			else

				!mineral soil thermal conductivities
				gwc = soil%wc(2)
				bde = ABD*0.06243
				if (soil%itxt .eq. 1) then !granular
					tku = (0.7*log10(gwc*100.0) + 0.4)*10.0**(0.01*bde)
					tkf = 0.076*10.0**(0.013*bde) + (gwc*100.0)*0.032*         &
						10.0**(0.0146*bde)
				else !fine-textured
					tku = (0.9*log10(gwc*100.0) - 0.2)*10.0**(0.01*bde)
					tkf = 0.01*10.0**(0.022*bde) + (gwc*100.0)*0.085*          &
						10.0**(0.008*bde)
				end if

				!convert from btu-in/ft^2/hr/f to kcal/m/hr/c
				tkf = tkf/8.0645
				tku = tku/8.0645

			end if

			!determine to use freezing or thawing conductivity
			if (knd .eq. 1) then
				tk = tkf
			else
				tk = tku
			end if

			!calculate depth of soil freezing or thawing (m)
			if (l .eq. 1) then
				ql = 80.0*(soil%wc(l))*OBD
			else
				ql = 80.0*(soil%wc(l))*ABD
			end if

			res = dl/tk
			degd = ql*dl/24.0*(tyr + res/2.0)

			if (degd .le. cadd) then
				depth = dl
				cadd = cadd - degd
			else
				ax = 0.5*ql/tk
				bx = ql*tyr
				cx = -24.0*cadd
				depth = (-bx + sqrt(bx*bx - 4.0*ax*cx))/(2.0*ax)
				depth = min(depth, dl)
				depth = max(depth, 0.0)
				cadd = 0.0

			end if

			tyr = tyr + res

			zdepth = zdepth + depth

		end do

		return

	end subroutine permf

	!:.........................................................................:

	subroutine moist(soil, ta, rain, pot_ev_day, nlai, slope, amlt, xmlt,      &
					xfrz, tdd, k, wd, aet, aow0_ScaledByMax, aow0_ScaledByMin, &
					saw0_ScaledByFC, saw0_ScaledByWP, saw0_ScaledBySAT)
		!calculates daily soil moisture dynamics
		!Author: Adrianna Foster 2018 v. 1.0
		!Inputs/Outputs:
		!	soil:   soil instance
		!Inputs:
		!	ta:               average temperature (degC)
		!	rain:             precipitation (cm)
		!   pot_ev_day:       potential evapotranspiration (cm)
		!	nlai:             leaf area index (m2/m2)
		!	slope: 	          slope of site (degrees)
		!	amlt:             active layer depth (m)
		!	xmlt: 		      depth of thaw (m)
		!	xfrz: 		      depth of freeze (m)
		!	tdd:		      thawing degree-days
		!	k:			      simulation month
		!Outputs
		!   aet:              actual evapotranspiration (cm)
		!	aow0_ScaledByMax: organic layer water scaled by field capacity
		!	aow0_ScaledByMin: organic layer water scaled by wilting point
		!	saw0_ScaledByFC:  A-layer water scaled by field capacity
		!	saw0_ScaledByWP:  organic layer water scaled by wilting point


		class(SoilData),        intent(inout) :: soil
		real, dimension(12, 2), intent(in)    :: tdd
		real,                   intent(in)    :: ta, rain, pot_ev_day
		real,                   intent(in)    :: amlt, xmlt, xfrz, nlai
		real,                   intent(in)    :: slope, wd
		integer,                intent(in)    :: k
		real,                   intent(out)   :: aet
		real,                   intent(out)   :: aow0_ScaledByMax
		real,                   intent(out)   :: aow0_ScaledByMin
		real,                   intent(out)   :: saw0_ScaledByFC
		real,                   intent(out)   :: saw0_ScaledByWP
		real,                   intent(out)   :: saw0_ScaledBySAT

		real, dimension(2)                    :: zh, wcf, gwl
		real, dimension(2)                    :: odmelt, dmelt, aet_loss
		real, dimension(2)                    :: odfreeze, dfreeze
		real, dimension(2)                    :: csm, pwp, fc, sat
		real, dimension(2)                    :: wt, water, owater, wf, wc
		real, dimension(2)                    :: zdrain, H2Oice, wc_w, wc_f
		real, dimension(2)                    :: exs, water_pet
		real, dimension(2)                    :: standing, water_exs
		real                                  :: canopy_int, lai, tfall
		real                                  :: laiw, laiw_min
		real                                  :: laiw_max, laiw0
		real                                  :: slope_fact, loss_slp
		real                                  :: can_evap
		real                                  :: xabove
		real                                  :: xcsm
		real                                  :: standing_b, laiw_int, laiw_evap
		real                                  :: precip, pet, pwl_can
		real                                  :: pwl, rz, pwll, pwg
		real                                  :: dthaw, B
		real                                  :: Od, Ad, OBD, ABD, pflood, rflood
		real                                  :: ps, pw, pr, thaw, pwl_init
		real, parameter                       :: meltfact = 4.0, tc = 0.0
		real, parameter                       :: tsmax = 3.3, tsmin = -1.1
		integer                               :: l


		!reducing table lookups
		dmelt = soil%d_melt
		dfreeze = soil%d_freeze
		Od = soil%O_depth
		Ad = soil%A_depth
		OBD = soil%O_bulk_dens
		ABD = soil%A_bulk_dens
		zdrain = soil%z_drain
		water = soil%water
		owater = soil%owater
		H2Oice = soil%H2Oice
		wc = soil%wc
		can_evap = 0.0
		canopy_int = 0.0
		laiw_int = 0.0
		laiw_evap = 0.0
		gwl = 0.0

		!calculate LAI metrics
		lai = max(nlai, 1.0)
		laiw_min = (lai*lai_min)/100.0
		laiw_max = (lai*lai_max)/100.0
		laiw0 = soil%lai_w0/100.0

		!initialize xabove, xcsm, and aet to 0.0
		xabove = 0.0
		xcsm = 0.0
		aet = 0.0
		exs(1) = 0.0
		exs(2) = 0.0

		!set previous day's depth of freeze/thaw
		odmelt(1) = dmelt(1)
		odmelt(2) = dmelt(2)
		odfreeze(1) = dfreeze(1)
		odfreeze(2) = dfreeze(2)

		!initialize csm, pwp to 0.0
		!set zh to appropriate soil layer depths
		csm(1) = 0.0
		csm(2) = 0.0
		pwp(1) = 0.0
		pwp(2) = 0.0
		zh(1) = soil%M_depth + Od
		zh(2) = Ad
		aet_loss = 0.0

		pflood = wd*0.01
		rflood = urand()

		!if (rflood .le. pflood) then

			!precip = (rain/100.0)*200.0*wd
		!else
			!calculate this day's water input (precip, m)
			!precip = rain/100.0 + wd*0.25
		!endif

		if (wd .gt. 0.0) then
			precip = rain/100 + wd*0.20
		else

			precip = rain/100.0
		endif

		!calculate daily PET (m)
		pet = pot_ev_day/100.0

		!partition precipitation between rain and snow
		if (ta .ge. tsmax) then
			pr = precip
			ps = 0.0

			!calculate canopy interception and throughfall
			canopy_int = min(max((laiw_max - laiw0), 0.0), pr)
			tfall = max(pr - canopy_int, 0.0)
			laiw = laiw0 + canopy_int
			laiw_int = laiw
			pw = tfall

		else if (ta .lt. tsmax .and. ta .gt. tsmin) then
			ps = (tsmax - ta)/(tsmax - tsmin)*precip
			pr = precip - ps

			!accumulate snowpack
			soil%swe = soil%swe + ps

			!calculate canopy interception and throughfall
			canopy_int = min(max((laiw_max - laiw0), 0.0), pr)
			tfall = max(pr - canopy_int, 0.0)
			laiw = laiw0 + canopy_int
			laiw_int = laiw
			pw = tfall

		else if (ta .le. tsmin) then
			pr = 0.0
			ps = precip

			!accumulate snowpack
			canopy_int = 0.0
			soil%swe = soil%swe + ps
			laiw = laiw0
			laiw_int = 0.0

			!no throughfall
			pw = 0.0
		end if

		!melt snowpack if it exists and if ta is above tc
		if (soil%swe .ge. 0.0001) then
			if (ta .gt. tc) then

				!thaw the snowpack
				thaw = (meltfact*(ta - tc))*0.001

				if (soil%swe .ge. thaw) then
					!decrement the snowpack
					soil%swe = soil%swe - thaw
				else
					!thaw whole snowpack - none left
					thaw = soil%swe
					soil%swe = 0.0
					soil%snowpack = 0.0
				end if
			else
				!temperature too cold to thaw
				thaw = 0.0
			end if
		else
			!no snowpack to thaw
			thaw = 0.0
		end if

		!calculate snow density
		soil%snowpack = soil%swe/(snow_dens/1000.0)

		tfall = pw + thaw

		!calculate losses due to slope runoff
		slope_fact = (slope/90.0)**2
		loss_slp = slope_fact*tfall

		!calculate potential water loss or gain this day
		pwl = tfall - loss_slp - pet

		!relative root density in forest floor organic layer
		if (zh(1) .gt. 0.0) then
			rz = min((zh(1) + amlt), 1.0)
			rz = 2.0*zh(1)/rz*(1.0 - zh(1)/(2.0*rz))
		else
			rz = 0.0
		end if

		!water released in soil thawing
		do l = 1, 2
			dthaw = max(0.0, (xmlt - xabove))
			dmelt(l) = min(zh(l), dthaw)
			dmelt(l) = max(dmelt(l), odmelt(l))
			wt(l) = zdrain(l)*(dmelt(l) - odmelt(l))

			!check to make sure enough ice for thawing amount
			if ((H2Oice(l) - wt(l)) .ge. 0.0001) then
				water(l) = water(l) + wt(l)
				owater(l) = water(l)
				H2Oice(l) = H2Oice(l) - wt(l)
			else
				water(l) = water(l) + H2Oice(l)
				owater(l) = water(l)
				H2Oice(l) = 0.0
			endif

			!set moisture conditions
			sat(l) = soil%ASAT(l)*dmelt(l)
			fc(l) = zdrain(l)*dmelt(l)
			pwp(l) = soil%APWP(l)*dmelt(l)

			xabove = xabove + zh(l)
		end do

		!add positive water balance to layer, adjusting for excess
		  !water from above layers
		do l = 1, 2
			if (l .eq. 1) then
				if (pwl .ge. 0.0) then
					!more water available than pet, aet = pet
					aet = pet
					!water gain is water from throughfall and snowmelt
					pwg = pwl
				else
					!no excess water from rainfall
					pwg = 0.0
				end if
			else
				if (pwl .ge. 0.0) then
					!excess water is from above layer
					pwg = exs(1)
				else
					!no excess water from rainfall/groundwater flow
					pwg = 0.0
				endif
			end if

			!add throughall/snowmelt/groundwater flow to layers
			water(l) = water(l) + pwg
			!drain off excess water
			exs(l) = max((water(l) - fc(l)), 0.0) !calculate excess water
			standing(l) = max((water(l) - sat(l)), 0.0)
			if (exs(l) .gt. 0.0) then
				if (soil%itxt .eq. 1) then
					gwl(l) = 2.0*((exs(l)**2)/(pet + exs(l)))*(1-(fc(l)-pwp(l)))
				else
					gwl(l) = 0.6*((exs(l)**2)/(pet + exs(l)))*(1-(fc(l)-pwp(l)))
				end if
				if (gwl(l) .le. 0.0) gwl(l) = 0.0
				if (gwl(l) .ge. 1.0) gwl(l) = 1.0
				exs(l) = gwl(l)*exs(l)
			end if
			water(l) = water(l) - exs(l) !subtract excess water
			water_exs(l) = water(l)
			csm(l) = water(l) - owater(l)
			xcsm = xcsm + csm(l)
		end do

		standing_b = standing(2)
		soil%runoff = exs(2) + loss_slp

		pwl_init = pwl
		if (pwl .lt. 0.0) then
			aet = 0.0
			!calculate canopy evaporation
			can_evap = min(-pwl, max(laiw - laiw_min, 0.0))
			laiw = max(laiw - can_evap, laiw_min)
			laiw_evap = laiw
			pwl = min(pwl + can_evap, 0.0)
			pwl_can = pwl
			aet = aet + can_evap
		end if

		!partition negative water balance between layers based on root density
		do l = 1, 2
			if (l .eq. 1) then
				pwll = min(0.0, pwl*rz) !either negative or 0.0
			else
				pwll = min(0.0, pwl*(1.0-rz)) !either negative or 0.0
			end if

			if (dmelt(l) .gt. 0.0) then !if any unfrozen soil
				B = 0.461 - 1.10559/(zdrain(l)*dmelt(l))
				aet_loss(l) = min(water(l) - water(l)*exp(B*abs(pwll)), -pwll)
				water(l) = max(water(l) - aet_loss(l), 0.0)
				aet = aet + aet_loss(l)
			else
				aet_loss(l) = 0.0
			end if
			water_pet(l) = water(l)
		end do

		!water content of soil layer when able to freeze
		do l = 1, 2
			if (tdd(k,1) .gt. 0.0) then
				if (dmelt(l) .gt. 0.0) then
					wcf(l) = water(l)/dmelt(l)
				else
					wcf(l) = 0.0
				end if
				odfreeze(l) = 0.0
			end if
		end do

		!water frozen as ice in soil freezing
		xabove = 0.0
		do l = 1, 2
			dfreeze(l) = min(zh(l), xfrz - xabove)
			dfreeze(l) = max(dfreeze(l), odfreeze(l))
			wf(l) = wcf(l)*(dfreeze(l) - odfreeze(l))

			!check to make sure enough liquid water for freezing amount
			if ((water(l) - wf(l)) .ge. 0.0001) then
				water(l) = water(l) - wf(l)
				H2Oice(l) = H2Oice(l) + wf(l)
			else
				H2Oice(l) = H2Oice(l) + water(l)
				water(l) = 0.0
			end if

			!gravimetric water content
			if (zh(l) .gt. 0.0) then
				wc(l) = (water(l)+H2Oice(l))/zh(l) !volumetric
				wc_w(l) = water(l)/zdrain(l)
				wc_f(l) = H2Oice(l)/zdrain(l)
			else
				wc(l) = 0.0
				wc_w(l) = 0.0
				wc_f(l) = 0.0
			end if
			if (l .eq. 1) then
				wc(l) = wc(l)*1000.0/OBD !volumetric to gravimetric
			else
				wc(l) = wc(l)*1000.0/ABD
			end if

			!adjust for above layers
			xabove = xabove + zh(l)
		end do

		soil%od_melt = odmelt
		soil%od_freeze = odfreeze
		soil%d_freeze = dfreeze
		soil%d_melt = dmelt
		soil%csm = csm
		soil%pwp = pwp
		soil%sat = sat
		soil%wt = wt
		soil%water = water
		soil%owater = owater
		soil%H2Oice = H2Oice
		soil%wf = wf
		soil%wc = wc
		soil%fc = fc
		soil%lai_w0 = laiw*100.0
		soil%runoff = soil%runoff*100.0


		if (soil%fc(1) .gt. 0.0) then
			aow0_ScaledByMax = ((soil%water(1) + soil%H2Oice(1))/zh(1))/(soil%fc(1))
		else
			aow0_ScaledByMax = 0.0
		end if

		if (soil%pwp(1) .gt. 0.0) then
			aow0_ScaledByMin = ((soil%water(1) + soil%H2Oice(1))/zh(1))/(soil%pwp(1))
		else
			aow0_ScaledByMin = 0.0
		end if

		if (soil%fc(2) .gt. 0.0) then
			saw0_ScaledByFC = (soil%water(2)/soil%fc(2))
		else
			saw0_ScaledByFC = 1.0
		end if

		!saw0_ScaledByFC = (soil%water(2) + soil%H2Oice(2))/soil%AFC(2)

		if (soil%pwp(2) .gt. 0.0) then
			saw0_ScaledByWP = soil%water(2)/(soil%pwp(2))
		else
			saw0_ScaledByWP = 2.0
		end if

		if (soil%sat(2) .gt. 0.0) then
			saw0_ScaledBySAT = (soil%water(2)/soil%sat(2))
		else
			saw0_ScaledBySAT = 0.0
		end if

		aet = aet*100.0

		return

	end subroutine moist

	!:.........................................................................:

	subroutine moss(soil, alff, cla, decLit, drydays)
		!calculates daily soil moisture dynamics
		!Author: Adrianna Foster 2018 v. 1.0
		!Inputs/Outputs:
		!	soil:   soil instance
		!Inputs:
		!	alff:    available light at forest floor
		!	cla:     cumulative leaf area
		!   decLit:  deciduous litter (t/ha)
		!	drydays: drought index

		class (SoilData), intent(inout) :: soil
		real,             intent(in)    :: alff, cla, decLit
		real,             intent(in)    :: drydays

		real                            :: biokg, biokg_1
		real                            :: al, algf, fcgf, dlgf, ddgf
		real                            :: assim, binp, binc1, binc
		real                            :: xx, mossbiom_tha

		real, parameter                 :: a1 = 3.41, a2 = 2.14
		real, parameter                 :: a3 = 0.08, q = 0.12
		real, parameter                 :: b1 = 0.136, ext = 3.5
		real, parameter                 :: q1 = 1.7, spores = 0.001
		real, parameter                 :: pmax = 0.2
		real, parameter                 :: cf = 5.2

		!convert moss biomass in kg to kg/m^2
		biokg = soil%moss_biom/plotsize

		!save old value of biomass (kg)
		biokg_1 = biokg

		!competition growth multiplier
		al = exp(-ext*q1*biokg/2.0)
		algf = a1*(al - a3)/(1.0 + a2*al)
		algf = max(0.0, algf)
		algf = min(1.0, algf)

		!available light growth multiplier (alff > 0.5)
		if (alff .gt. 0.5) then
			fcgf = 1.25 - alff**2.0
		else
			fcgf = 1.0
		end if

		!deciduous leaf litter growth multiplier
		if (cla .gt. 0.0) then
			if (decLit .gt. 0.0) then
				dlgf = 1.0 - 0.0095*(decLit)**2
				if (dlgf .le. 0.0) dlgf = 0.0
				if (dlgf .ge. 1.0) dlgf = 1.0
			else
				dlgf = 1.0
			end if
		else
			dlgf = 1.0
		end if

		!soil moisture growth multiplier
		if (drydays .gt. 0.10) then
			ddgf = 0.0
		else
			ddgf = 1.0
		end if

		!grow moss
		if (alff .le. 0.5) then
			assim = pmax*(1.27 + 0.3*sqrt(alff))
		else
			assim = pmax
		end if

		assim = assim*algf*ddgf*fcgf*dlgf

		binp = spores*ddgf*dlgf
		binc1 = q1*biokg*assim - q1*biokg*(q + b1) + binp

		if ((biokg + binc1) .lt. 0.0) then
			binc = -biokg
		else
			binc = binc1
		end if

		soil%moss_biom = (biokg + binc)*plotsize !kg/plot

		!calculate litter
		xx = (q1*biokg*(assim - q) + binp - binc)*plotsize !kg/plot
		if (xx .lt. 0.0) then
			soil%litter(18) = 0.0
		else
			soil%litter(18) = xx
		end if
		
		!convert to tonnes/ha from kg/plot
		soil%litter(18) = soil%litter(18)/plotsize/0.0001/1000

		!moss in t/ha
		mossbiom_tha = soil%moss_biom/plotsize/0.0001/1000

		!thickness of moss layer (m)
		soil%M_depth = soil%moss_biom/plotsize/cf


	end subroutine moss

	!:.........................................................................:

end module Soil
