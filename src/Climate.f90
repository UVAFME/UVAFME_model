module Climate

  use Constants
  use Random

  implicit none

!*******************************************************************************
  !
  !This module contains subroutines to compute climate-related variables
  !
!*******************************************************************************

  real                              :: accumulated_tmin
  real                              :: accumulated_tmax
  real, dimension(12)               :: accumulated_precip

contains

	!:.........................................................................:

	subroutine set_site_climate()
		!sets accumulated tmin, tmax, and precip to 0 for potential climate
		!change

		!Note - ACF 4/30/19: I'm not sure why this is here as an individual
		 !subroutine

		accumulated_tmin = 0.0
		accumulated_tmax = 0.0
		accumulated_precip = 0.0

	end subroutine set_site_climate

	!:.........................................................................:

	subroutine cov365(month_vals, day_vals)
		!converts monthly state data into daily data
		!Author: Yan Xiaodong, 1996 v. 1.0
		!Inputs:
		!   month_vals: vector of monthly state data (temperature/cloud cover)
		!Outputs:
		!   day_vals:   vector of daily weather data

		real, dimension(:), intent(in)    :: month_vals
		real, dimension(:), intent(inout) :: day_vals

		real, dimension(13)               :: monthly_vals
		real, dimension(381)              :: dayly_vals
		real                              :: mean_val
		integer                           :: middays(13), k, md
		data middays/16, 45, 75, 105, 136, 166, 196, 227, 258, 288, 319, 349,  &
						381/

		monthly_vals(13) = month_vals(1)
		do k = 1, 12
			monthly_vals(k) = month_vals(k)
		end do

		do k = 1, 12
			mean_val = (monthly_vals(k + 1) - monthly_vals(k))/                &
				float(middays(k + 1) - middays(k))
			do md = middays(k), middays(k + 1)
				dayly_vals(md) = monthly_vals(k) + mean_val*float(md -         &
					middays(k))
			end do
		end do

		do md = 16, 365
			day_vals(md) = dayly_vals(md)
		end do

		do md = 1, 15
			day_vals(md) = dayly_vals(365 + md)
		end do

		return

	end subroutine cov365

	!:.........................................................................:

	subroutine cov365a(month_vals, day_vals)
		!converts monthly integrated data into daily data randomly
		!Author: Yan Xiaodong, 1996 v. 1.0
		!Inputs:
		!   month_vals: vector of monthly integrated data (precipitation)
		!Outputs:
		!   day_vals:   vector of daily weather data

		real, dimension(:), intent(in)    :: month_vals
		real, dimension(:), intent(inout) :: day_vals

		real                              :: raindays, r_day, p_rday, r_rand
		integer                           :: modays(12), md, i, k, iraindays
		integer                           :: inum
		data modays/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

		!number of rain days is function of monthly rainfall amount (cm)
			!basic feature is: 100 cm = 25 days, 1 cm rain = 1 day

		md = 0
		do k = 1, 12
			!calculate how many days it rains this month
			raindays = min1(25.0, month_vals(k)/4.0 + 1.0)
			iraindays = int(raindays)

			!calculate how much rain per rain day
			r_day = month_vals(k)/float(iraindays)

			!calculate what percentage of days it rains this month
			p_rday = raindays/modays(k)

			!for each day in month, add rain amount to each day randomly such
				!that it adds up to amount of rain that month
			inum = iraindays
			do i = 1, modays(k)
				md = md + 1
				!if any rain left to add
				if (inum .gt. 0) then
					r_rand = clim_urand()
					if(r_rand .le. p_rday) then
						!add amount of rain each day to that day's rain
						day_vals(md) = r_day
						!decrement number of rain days
						inum = inum - 1
					else
						day_vals(md) = 0.0
					end if
				else
					day_vals(md) = 0.0
				end if
			end do

			!if any rain days left over add rest of rain to middle of month
			if (inum .gt. 0) then
				day_vals(md - 15) = float(inum)*r_day
			end if

		end do

		return

	end subroutine cov365a

	!:.........................................................................:

	subroutine ex_rad(jd, latit, slope, aspect, cld, erad, sun, st, sunrise_hr)
		!calculates daily TOA and surface solar radiation
		!adapted from Bonan 1989 Ecological Modelling 45:275-306
		!Author: Adrianna Foster, 2018 v. 1.0
		!Inputs:
		!  jd:         Julian Day
		!  latit:      latitude of site (degrees)
		!  slope:      slope of site (degrees)
		!  aspect:     aspect of site (degrees)
		!  cld:        cloud cover (tenths of sky covered)
		!Outputs:
		!  erad:       top-of-atmosphere solar radiation (MJ/m2/day)
		!  sun:        surface solar radiation (cal/cm2/day)
		!  st:         horizontal surface solar radiation (cal/cm2/day)
		!  sunrise_hr: hour of sunrise

		integer,  intent(in)   :: jd
		real,     intent(in)   :: latit, slope, aspect, cld
		real,     intent(out)  :: erad, sun, st
		integer,  intent(out)  :: sunrise_hr

		real                   :: temp_omega, rlat, zs, azw, dist, dec
		real                   :: omega, erad1, erad2, sol_hr, sol_alt
		real                   :: azs, ang1, ang2, xx, yy, st1, df
		real                   :: dr, angsum, altsum, akt, Rb, Rd
		real                   :: A
		integer                :: hr
		integer, dimension(24) :: sun_hr

		angsum = 0.0
		altsum = 0.0

		!convert latitude to radians
		rlat = deg2rad*latit

		!convert slope degrees to radians
		zs = deg2rad*slope

		!convert aspect from degrees to radians and calculate wall azimuth angle
		azw = (180.0 - aspect)*pi/180.0

		!calculate earth-sun distance
		dist = 1.0 + Ac*cos(2*pi*float(jd)/365.0)

		!calculate solar declination (radians)
		A = (284.0 + float(jd))/365.0*2.0*pi
		dec = 23.45*sin(A)*pi/180.0

		!initial value for sunset hour angle (radians)
		temp_omega = -tan(rlat)*tan(dec)

		!modify sunset hour angle based on sunset/sunrise occurence
		if (temp_omega .ge. 1.0) then !sun never rises
		   omega = 0.0
		else if (temp_omega .le. -1.0) then !sun never sets
		   omega = pi
		else !all other times
		   omega = acos(temp_omega)
		end if

		!calculate extraterrestrial radiation (MJ/m2/day)
		erad1 = Amp*dist*cos(rlat)*cos(dec)*(sin(omega) - omega*cos(omega))
		erad = max(erad1, 0.0)

		!convert to cal/cm2/day
		erad2 = erad*239.006*1000/10000

		!calculate hourly solar incidence and altitude angles
		do hr = 1, 24
			sol_hr = 15.0*(12.0 - float(2*hr - 1)/2.0)*pi/180
			sol_alt = asin(sin(rlat)*sin(dec) + cos(rlat)*cos(dec)*cos(sol_hr))

			!azimuth angle of sun from south
			azs = asin(cos(dec)*sin(sol_hr)/cos(sol_alt))

			!solar incidence angle
			ang1 = sin(dec)*sin(rlat)*cos(zs) -                                &
				   sin(dec)*cos(rlat)*sin(zs)*cos(azw) +                       &
				   cos(dec)*cos(sol_hr)*cos(rlat)*cos(zs) +                    &
				   cos(dec)*cos(sol_hr)*sin(rlat)*sin(zs)*cos(azw) +           &
				   cos(dec)*sin(zs)*sin(azw)*sin(sol_hr)

			ang2 = acos(ang1)

			!sunrise/sunset occurs when ang2 = pi/2 or sol_alt = 0
			if (ang2 .ge. (pi/2.0) .or. sol_alt .le. 0.0) then
				xx = 0.0
				yy = 0.0
				sun_hr(hr) = 0
			else
				xx = cos(ang2)
				yy = sin(sol_alt)
				sun_hr(hr) = 1
			end if

			angsum = angsum + xx
			altsum = altsum + yy

		end do

		sunrise_hr = 12
		!find sunrise hour
		do hr = 1, 24
			if (sun_hr(hr) > 0) then
				sunrise_hr = hr
				exit
			endif

		end do

		!attenuate solar radiation through atmosphere
		st1 = -7.13 + 0.812*erad2-0.44*erad2*(cld/10.0)
		st = max(st1, 0.0)

		if (erad2 .gt. 0.0) then
			akt = st/erad2
		else
			akt = 0.0
		end if

		df = (1.0045 + 0.04349*akt - 3.5227*akt**2.0 + 2.6313*akt**3.0)*st
		if (akt .gt. (0.75)) df = 0.166*st
		if (df .lt. 0.0) df = 0.0
		dr = st - df

		!calculate direct beam tilt factor
		if (altsum .ne. 0.0) then
			Rb = angsum/altsum
		else
			Rb = 0.0
		end if

		!calculate diffuse radiation tilt factor
		Rd = (cos(zs/2.0))**2.0

		!calculate mean solar radiation at surface
		sun = Rb*dr + Rd*df

		return

	end subroutine ex_rad

	!:.........................................................................:

	function jensenhaise(ta, sun, elev, e2, e1)
		!calculation of potential evapotranspiration (cm/day)
		!adapted from Bonan 1989 Ecological Modelling 45:275-306
		!Author: Adrianna Foster, 2018 v. 1.0
		!Inputs:
		!   ta:    mean daily air temperature (degrees C)
		!   sun:   mean solar radiation (cal/cm2/day)
		!   elev:  elevation of site (m)
		!   e2:    saturation vapor pressure (mbar) at tmax of warmest
		!           month of year
		!   e1:    saturation vapor pressure (mbar) at tmin of warmest
		!           month of year
		!Outputs:
		!   jensenhaise: potential evapotranspiration (cm/day)

		real                :: jensenhaise
		real,    intent(in) :: ta, sun, elev, e2, e1

		real                :: a, b, hvap, evap

		!calculate latent heat of vaporization
		hvap = 597.391 - 0.5680*ta !kcal/kg

		!calculate coefficients
		a = 1.0/(38.0 - (2.0*elev/305.0) + 380.0/(e2 - e1))
		b = -2.5  - 0.14*(e2 - e1) - elev/550.0

		!calculate PET (cm)
		if (ta .gt. 0.0) then
			evap = a*(ta - b)*sun/hvap
			evap = max(evap, 0.0)
		else
			evap = 0.0
		end if

		jensenhaise = evap

	end function jensenhaise

	!:.........................................................................:

	function hrsabove(tmin, tmax, hrise, thresh)
		!calculates cumulative hours above a threshold for the day
		!Author: Adrianna Foster, 2016 v. 1.0
		!Changelog: ACF - 4/30/10, v. 2.0
		!Inputs:
		!   tmin:     minimum temperature (degrees C)
		!   tmax:     maximum temperature (degrees C)
		!   hrise:    hour of sunrise
		!   thresh:   threshold temperature (degrees C)
		!Outputs:
		!   hrsabove: hours above threshold temperature


		real                    :: hrsabove
		real,    intent(in)     :: tmin, tmax, thresh
		integer, intent(in)     :: hrise

		real, dimension(24)     :: td
		real                    :: tmean
		integer                 :: i, couni
		integer                 :: hmax = 14
		integer, dimension(24)  :: h

		data h/0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,   &
				18, 19, 20, 21, 22, 23/

		couni = 0

		tmean = (tmax + tmin)/2.0

		!calculate hourly temperature for the day
		do i = 1, 24
			if ((h(i) .ge. 0) .and. (h(i) .lt. hrise)) then
				td(i) = tmean + ((tmax - tmin)/2.0)*                           &
					cos((pi*(float(h(i)) + 10.0))/(10.0 + hrise))
		    else if ((h(i) .ge. hrise) .and. (h(i) .le. hmax)) then
				td(i) = tmean - ((tmax - tmin)/2.0)*                           &
					cos((pi*(float(h(i)) - hrise))/(hmax - hrise))
			else if ((h(i) .gt. hmax) .and. (h(i) .lt. 24)) then
				td(i) = tmean + ((tmax - tmin)/2.0)*                           &
					cos((pi*(float(h(i)) + hmax))/(10.0 + hrise))
			end if
		end do

		!accumulate hours above the threshold temperature
		do i = 1, 24
			if (td(i) .gt. thresh) then
				couni = couni + 1
			endif
		enddo

		hrsabove = couni

	end function hrsabove

	!:.........................................................................:

end module Climate
