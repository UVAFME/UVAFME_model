module Random

	implicit none

!*******************************************************************************
  !
  !This module contains various random number generators.
  !
!*******************************************************************************


!variables for the climate-specific RNG
integer   :: clim_seed = 0
integer   :: rng_seed
logical   :: reset_clim = .true.

private      clim_seed, rng_seed, reset_clim

contains

	!:.........................................................................:

	function urand(lb, ub, seed)
		!returns a uniformly-distributed random number in the range lb to ub,
		  !default 0.0 to 1.0
		!Author: Katherine Holcomb 2012, v. 1.0
		!Inputs/Outputs:
		!	lb:    lower bound
		!	ub:    upper bound
		!	seed:  random seed
		!Outputs:
		!	urand: uniform random number

		real                       :: urand
		real, optional, intent(in) :: lb, ub
		real, optional, intent(in) :: seed
		integer                    :: iseed
		real                       :: rnd
		real                       :: lower, upper

		if (present(seed)) then
			iseed = int(seed)
			call set_random_seed(iseed)
		endif

		if (present(lb)) then
			lower = lb
		else
			lower = 0.0
		endif

		if (present(ub)) then
			upper = ub
		else
			upper = 1.0
		endif

		call random_number(rnd)
		urand = lower + (upper - lower)*rnd

		return

	end function urand

	!:.........................................................................:

	function nrand(mean, std, seed)
		!returns a normally-distributed random number with given mean
		  !and standard deviation, default to mean = 0.0, std = 1.0
		!Author: Katherine Holcomb 2012, v. 1.0
		!Inputs/Outputs:
		!	mean:   mean
		!	std:    standard deviation
		!   seed:   random seed
		!Outputs:
		!	nrand:  normally-distributed random numer


		real                       :: nrand
		real, optional, intent(in) :: mean, std
		real, optional, intent(in) :: seed
		integer                    :: iseed
		real                       :: mn, st
		real                       :: x1, y1, x2, y2, w

		if (present(mean)) then
			mn = mean
		else
			mn = 0.0
		endif

		if (present(std)) then
			st = std
		else
			st = 1.0
		endif

		!Box-Muller polar method
1 		continue

		if (present(seed)) then
			iseed = int(seed)
			x1 = urand(-1.0, 1.0, seed)
			x2 = urand(-1.0, 1.0, seed)
		else
			x1 = urand(-1.0, 1.0)
			x2 = urand(-1.0, 1.0)
		endif

		w = x1**2 + x2**2
		if (w .eq. 0.0 .or. w .gt. 1.0) go to 1

		w = sqrt((-2.0 * log( w ))/w)
		y1 = x1*w
		y2 = x2*w

		!pick one, adjust its mean and std
		nrand = y1*st + mn

		return

	end function nrand

	!:.........................................................................:

	function clim_urand(lb, ub)
		!used in order to have a separate RNG for climate data. It is
		  !simple-minded but decent linear congruence with a shuffle, from the
		  !Numerical Recipes First Edition. Passes "Runs" test but probably not
		  !DIEHARD.
		!Author: Katherine Holcomb 2012, v. 1.0
		!Inputs/Outputs:
		!	lb:         lower bound
		!	ub:         upper bound
		!Outputs:
		!	clim_urand: uniform random number

		real                       :: clim_urand
		real, optional, intent(in) :: lb, ub
		real                       :: lower, upper
		real                       :: ran1
		real                       :: rm1, rm2
		real, dimension(97)        :: r
		integer, parameter         :: m1 = 259200, ia1 = 7141
		integer, parameter         :: ic1 = 54773
		integer, parameter         :: m2 = 134456, ia2 = 8121
		integer, parameter         :: ic2 = 28411
		integer, parameter         :: m3 = 243000, ia3 = 4561
		integer, parameter         :: ic3 = 51349
		integer                    :: ix1, ix2, ix3
		integer                    :: j
		logical                    :: first = .true.
		save

		if (present(lb)) then
			lower = lb
		else
			lower = 0.0
		endif

		if (present(ub)) then
			upper = ub
		else
			upper = 1.0
		endif

		rm1 = 1.0/m1
		rm2 = 1.0/m2

		if (first .or. reset_clim) then
			ix1 = mod(abs(ic1 - clim_seed), m1)
			ix1 = mod(ia1*ix1 + ic1, m1)
			ix2 = mod(ix1, m2)
			ix1 = mod(ia1*ix1 + ic1, m1)
			ix3 = mod(ix1, m3)
			do j = 1, 97
				ix1 = mod(ia1*ix1 + ic1, m1)
				ix2 = mod(ia2*ix2 + ic2, m2)
				r(j) = (real(ix1) + real(ix2)*rm2)*rm1
			enddo
			first = .false.
		  reset_clim = .false.
		endif

		ix1 = mod(ia1*ix1 + ic1, m1)
		ix2 = mod(ia2*ix2 + ic2, m2)
		ix3 = mod(ia3*ix3 + ic3, m3)

		j = 1 + (97*ix3)/m3
		if (j .gt. 97 .or. j .lt. 1) then
			j = 97
			write(*, *) 'Error in climate RNG'
		else
			ran1 = r(j)
			r(j) = (real(ix1) + real(ix2)*rm2)*rm1
		endif

		clim_urand = lower + (upper - lower)*ran1

		return

	end function clim_urand

	!:.........................................................................:

	function clim_nrand(mean, std)
		!returns a normally-distributed random number with given mean and standard
		  !deviation, default is mean = 0.0 and std = 1.0
		!Author: Katherine Holcomb 2012, v. 1.0
		!Inputs/Outputs:
		!	mean:       mean
		!	std:        standard deviation
		!Outputs:
		!	clim_nrand: normally distributed random number

		real           :: clim_nrand
		real, optional :: mean, std
		real           :: mn, st
		real           :: x1, y1, x2, y2, w

		if (present(mean)) then
			mn = mean
		else
			mn = 0.0
		endif

		if (present(std)) then
			st = std
		else
			st = 1.0
		endif

		!Box-Muller polar method
1 		continue

		x1 = clim_urand(-1.0, 1.0)
		x2 = clim_urand(-1.0, 1.0)

		w = x1**2 + x2**2
		if (w .eq. 0.0 .or. w .gt. 1.0) go to 1

		w = sqrt((-2.0 * log(w))/w)
		y1 = x1*w
		y2 = x2*w

		!pick one, adjust its mean and std
		clim_nrand = y1*st + mn

		return

	end function clim_nrand

	!:.........................................................................:

	subroutine set_site_rng_seed(fixed_seed, seed)
		!selets the site's random number generator seed. Fixed_seed is for
		  !debugging. If true, set the seed to be the default seed each time
		  !it is called. If false, get a random seed (or optionally use a
		  !specified seed) on the first call
		!Author: Katherine Holcomb 2012, v. 1.0
		!Inputs/Outputs:
		!	fixed_seed:  should a fixed seed be used (.true./.false.)
		!	siteID:      unique site iD
		!	seed:        random seed

		logical,        intent(in) :: fixed_seed
		real, optional, intent(in) :: seed
		integer                    :: idate(8)
		integer                    :: iseed
		integer, parameter         :: default_seed = 2345678
		logical                    :: first = .true.
		real                       :: rnd

		!If fixed_seed is true, reset the seed for each site. This is done
		  !because if we ever parallelize the code, the only way to ensure
		  !reproducibility is to start each site from the same state. It's not
		  !so important for serial runs but it doesn't take all that much time.
		  !Even for serial runs, however, it means we can run the sites in any
		  !order.

		if (.not. fixed_seed) then
			if (first) then
				if (.not. present(seed)) then
					!create seed value
					!get date and time values
					call date_and_time(values = idate)

					!get a random seed
					call get_random_seed(iseed)

					!calculate actual seed value
					if (iseed .ne. 0) then
						iseed = iseed * (idate(8)) !idate(8) = milliseconds
					else
						iseed = default_seed * (idate(8))
					endif

				else
					!iseed is equal to input seed
					iseed = seed
				endif
				first = .false.
				rng_seed = iseed

				!set the random seed to the calculated or input seed for any
				  !future random number calls
				call set_random_seed(rng_seed)
			endif
		else
			!fixed seed - use same seed for all sites. Seed is set to default
		      !seeds unless seed is present in call
			if (.not. present(seed)) then
				iseed = default_seed
			else
				iseed = seed
			endif
			rng_seed = iseed

			!set the random number seed and call a random number
			call set_random_seed(rng_seed)
			call random_number(rnd)

		endif

	end subroutine set_site_rng_seed

	!:.........................................................................:

	subroutine get_random_seed(seed)
		!gets a random seed
		!Author: Katherine Holcomb 2012, v. 1.0
		!Outputs:
		!	seed:  random seed

		integer, intent(out)               :: seed
		integer                            :: isize
		integer, dimension(:), allocatable :: iseed

		!get random seed size and allocate
		call random_seed(size = isize)
		if (.not. allocated(iseed)) allocate(iseed(isize))

		!get a random seed of size isize
		iseed = seed
		call random_seed(get = iseed)
		seed = iseed(1)

	end subroutine get_random_seed

	!:.........................................................................:

	subroutine set_random_seed(seed)
		!sets the random seed for random_number calls
		!Author: Katherine Holcomb 2012, v. 1.0
		!Inputs
		!	seed:  random seed

		integer, intent(in)                :: seed
		integer                            :: isize
		integer, dimension(:), allocatable :: iseed

		!get size of random seed
		call random_seed(size = isize)
		if (.not. allocated(iseed)) allocate(iseed(isize))

		!put the seed as the random seed for any future random number calls
		iseed = seed
		call random_seed(put = iseed)

	end subroutine set_random_seed

	!:.........................................................................:
	subroutine shuffle(a)
		!shuffles the array a
		!Author: Katherine Holcomb 2012, v. 1.0
		!Inputs:
		!	a: array to shuffle

		integer, dimension(:), intent(inout) :: a
		integer                              :: i
		integer                              :: randpos,tempval
		real                                 :: rnd

		do i = size(a), 1, -1
			call random_number(rnd)
			randpos = int(rnd*i) + 1
			tempval = a(randpos)
			a(randpos) = a(i)
			a(i) = tempval
		enddo

	end subroutine shuffle

	!:.........................................................................:

end module Random
