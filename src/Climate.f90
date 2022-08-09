module Climate

!*******************************************************************************
  !
  ! This module contains subroutines to compute climate-related variables
  !
!*******************************************************************************

  use Constants
  use Random

  implicit none

contains

    !:.........................................................................:

    subroutine cov365_state(month_vals, day_vals)
        !
        !  Converts monthly state data into daily data
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/96     Y. Xiaodong         Original Code
        !

        ! Data dictionary: constants:
        ! Julian day for middle of each month, plus one for wrapping around
        integer, parameter, dimension(13) :: MIDDAYS = [16, 45, 75, 105, 136,  &
            166, 196, 227, 258, 288, 319, 349, 381]

        ! Data dictionary: calling argumnets
        real, dimension(:), intent(in)    :: month_vals ! Monthly state data
        real, dimension(:), intent(inout) :: day_vals   ! Daily data

        ! Data dictionary: local variables
        real, dimension(13)  :: monthly_vals ! Local array for monthly values
        real, dimension(381) :: daily_vals   ! Local array for daily values
        real                 :: mean_val     ! Mean value
        integer              :: k, md        ! Looping indices

        ! Set 13th index to January so we can wrap around
        monthly_vals(13) = month_vals(1)

        ! Fill the rest with input values
        do k = 1, 12
            monthly_vals(k) = month_vals(k)
        end do

        do k = 1, 12
            ! Average value between next month's and this month, scaled by
            ! days inbetween the midpoints
            mean_val = (monthly_vals(k + 1) - monthly_vals(k))/                &
                float(MIDDAYS(k + 1) - MIDDAYS(k))

            ! Average plus scaling towards the midpoint
            do md = MIDDAYS(k), MIDDAYS(k + 1)
                daily_vals(md) = monthly_vals(k) + mean_val*float(md -         &
                    MIDDAYS(k))
            end do
        end do

        ! Fill with daily values, accounting for wrapping around
        do md = 16, 365
            day_vals(md) = daily_vals(md)
        end do
        do md = 1, 15
            day_vals(md) = daily_vals(365 + md)
        end do

        return

    end subroutine cov365_state

    !:.........................................................................:

    subroutine cov365_integr(month_vals, day_vals)
        !
        !  Converts monthly integrated data (i.e. rainfall) into daily data
        !  Number of rain days is a function of monthly rainfall amount (cm)
        !  Basic feature is: 100 cm = 25 days, 1 cm rain = 1 day
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/96     Y. Xiaodong         Original Code
        !

        ! Data dictionary: constants
        ! Number of days in each month
        integer, dimension(12), parameter :: MODAYS = [31, 28, 31, 30, 31, 30, &
            31, 31, 30, 31, 30, 31]

        ! Data dictionary: calling arguments
        real, dimension(:), intent(in)    :: month_vals ! Monthly integrated data
        real, dimension(:), intent(inout) :: day_vals   ! Daily data

        ! Data dictionary: local variables
        real     :: raindays  ! Number of days with rainfall in a month
        real     :: r_day     ! How much it rains per day (cm)
        real     :: p_rday    ! Proportion of rainy days in a month (0-1)
        real     :: r_rand    ! Random number (uniform)
        integer  :: iraindays ! Number of days with rainfall in a month
        integer  :: md, i, k  ! Looping indices
        integer  :: inum      ! Number of rain days left to distribute

        md = 0
        do k = 1, 12

            ! Calculate how many days it rains each month
            raindays = min1(25.0, month_vals(k)/4.0 + 1.0)
            iraindays = int(raindays)

            ! Calculate how much rain per rain day
            r_day = month_vals(k)/float(iraindays)

            ! Calculate what percentage of days it rains this month
            p_rday = raindays/MODAYS(k)

            ! For each day in month, add rain amount to each day randomly such
            ! that it adds up to amount of rain that month
            inum = iraindays
            do i = 1, MODAYS(k)
                md = md + 1
                ! If any rain left to add
                if (inum > 0) then
                    ! Get random number
                    r_rand = clim_urand()
                    if(r_rand <= p_rday) then
                        ! Add amount of rain each day to that day's rain
                        day_vals(md) = r_day
                        ! Decrement number of rain days left to add
                        inum = inum - 1
                    else
                        day_vals(md) = 0.0
                    end if
                else
                    day_vals(md) = 0.0
                end if
            end do

            ! If any rain days left over add rest of rain to middle of month
            if (inum > 0) then
                day_vals(md - 15) = float(inum)*r_day
            end if

        end do

        return

    end subroutine cov365_integr

    !:.........................................................................:

    subroutine solar_rad(jd, latit, slope, aspect, cld, erad, sun, st,         &
            sunrise_hr)
        !
        !  Calculates daily top-of-atmosphere and surface solar radiation.
        !  Also returns hour of sunrise.
        !
        !  Adapted from Bonan 1989 Ecological Modelling 45:275-306
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    10/10/18     A. C. Foster         Original Code
        !

        ! Data dictionary: constants:
        real, parameter :: S = 2809.44 ! Solar constant (cal/cm2/day)
        ! Parameters for attenuating solar radiation through atmosphere (Bonan 1989)
        ! These are for North America
        real, dimension(3), parameter :: AT = [-7.130, 0.812, 0.440]

        ! Data dictionary: calling arguments
        integer,  intent(in)  :: jd         ! Julian day
        real,     intent(in)  :: latit      ! Latitude (degrees)
        real,     intent(in)  :: slope      ! Slope (degrees)
        real,     intent(in)  :: aspect     ! Aspect (degrees)
        real,     intent(in)  :: cld        ! Cloud cover (tenths of sky covered)
        real,     intent(out) :: erad       ! Top of atmosphere radiation (cal/cm2/day)
        real,     intent(out) :: sun        ! Surface solar radiation (cal/cm2/day)
        real,     intent(out) :: st         ! Horizontal surface solar radiation (cal/cm2/day)
        integer,  intent(out) :: sunrise_hr ! Hour of sunrise

        ! Data dictionary: local variables
        integer, dimension(24) :: sun_hr     ! Hours where we have sun
        real                   :: rlat       ! Latitude (radians)
        real                   :: rslope     ! Slope (radians)
        real                   :: az_wall    ! Wall azimuth angle (radians)
        real                   :: r_dist     ! Relative Earth-Sun distance (AU; 1 AU = 1.485E11 m)
        real                   :: A          ! Temporary variable for calculating solar declination
        real                   :: dec        ! Solar declination (radians)
        real                   :: temp_omega ! Initial value for sunset hour angle (radians)
        real                   :: omega      ! Sunset/sunrise hour angle (radians)
        real                   :: sol_hr     ! Solar hour angle (radians)
        real                   :: sol_alt    ! Solar altitude angle (radians)
        real                   :: azs        ! Azimuth angle of sun from south (radians)
        real                   :: ang1       ! Temporary variable for calculating incidence angle
        real                   :: sol_inc    ! Solar incidence angle (radians)
        real                   :: xx, yy     ! Temporary variables for summing altitude & incidence angles
        real                   :: angsum     ! Sum of hourly incidence angles (radians)
        real                   :: altsum     ! Sum of hourly altitude angles (radians)
        real                   :: akt        ! Fraction of TOA radiation attenuated by atmosphere
        real                   :: df         ! Horizontal surface diffuse radiation (cal/cm2/day)
        real                   :: dr         ! Horizontal surface direct beam radiation (cal/cm2/day)
        real                   :: Rb         ! Direct beam tilt factor
        real                   :: Rd         ! Diffuse tilt factor
        integer                :: hr         ! Looping index


        ! Initialize accumulators
        angsum = 0.0
        altsum = 0.0

        ! Convert latitude to radians
        rlat = DEG2RAD*latit

        ! Convert slope degrees to radians
        rslope = DEG2RAD*slope

        ! Calculate wall azimuth angle
        az_wall = (180.0 - aspect)*pi/180.0

        ! Calculate realtive Earth-Sun distance (AU)
        ! Brock 1981 Ecological Modelling 14:1-19
        ! Nicholls & Child 1979 Solar Energy 22(1)
        !r_dist = 1.0 + 0.033*cos(2*pi*float(jd)/365.0)
        r_dist = 1.0/(1.0 + (0.033*cos(360.0*float(jd)/365.0)))**(1/2)

        ! Calculate solar declination (radians)
        ! Cooper 1969
        A = ((284.0 + float(jd)))/365.0*2.0*pi
        dec = 23.45*sin(A)*pi/180.0

        ! Initial value for sunset hour angle (radians)
        ! Brock 1981 Ecological Modelling 14:1-19
        temp_omega = -tan(rlat)*tan(dec)

        ! Modify sunset hour angle based on sunset/sunrise occurence
        ! Keith & Kreider 1978 Principles of Solar Engineering
        if (temp_omega >= 1.0) then
           ! Sun never rises
            omega = 0.0
        else if (temp_omega <= -1.0) then
            ! Sun never sets
            omega = pi
        else
            ! Calculate angle
            omega = acos(temp_omega)
        end if

        ! Calculate top of atmosphere radiation (cal/cm2/day)
         erad = S/pi*(1.0/r_dist**2)*cos(rlat)*cos(dec)*(sin(omega) -          &
            omega*cos(omega))
         erad = max(erad, 0.0)

        ! Calculate hourly solar incidence and altitude angles
        ! Bonan 1989
        do hr = 1, 24

            ! Solar hour angle (radians)
            sol_hr = 15.0*(12.0 - float(2*hr - 1)/2.0)*pi/180.0

            ! Solar altitude angle (radians)
            sol_alt = asin(sin(rlat)*sin(dec) + cos(rlat)*cos(dec)*cos(sol_hr))

            ! Azimuth angle of sun from south (radians)
            azs = asin(cos(dec)*sin(sol_hr)/cos(sol_alt))

            ! Solar incidence angle (radians)
            ang1 = sin(dec)*sin(rlat)*cos(rslope) -                            &
                   sin(dec)*cos(rlat)*sin(rslope)*cos(az_wall) +               &
                   cos(dec)*cos(sol_hr)*cos(rlat)*cos(rslope) +                &
                   cos(dec)*cos(sol_hr)*sin(rlat)*sin(rslope)*cos(az_wall) +   &
                   cos(dec)*sin(rslope)*sin(az_wall)*sin(sol_hr)
            sol_inc = acos(ang1)

            ! Sum up incidence and altitude angles
            if (sol_inc >= (pi/2.0) .or. sol_alt <= 0.0) then
                ! Sun is not above horizon
                xx = 0.0
                yy = 0.0
                sun_hr(hr) = 0
            else
                xx = cos(sol_inc)
                yy = sin(sol_alt)
                sun_hr(hr) = 1
            end if

            angsum = angsum + xx
            altsum = altsum + yy

        end do

        sunrise_hr = 12
        ! Find sunrise hour
        do hr = 1, 24
            if (sun_hr(hr) > 0) then
                sunrise_hr = hr
                exit
            endif
        end do

        ! Calculate horizontal surface radiation (cal/cm2/day)
        ! Bonan 1989
        st = max(AT(1) + AT(2)*erad - AT(3)*erad*(cld/10.0), 0.0)

        ! Calculate fraction of TOA radiation attenuated by atmosphere
        if (erad > 0.0) then
            akt = st/erad
        else
            akt = 0.0
        end if

        ! Calculate horizontal surface diffuse radiation (cal/cm2/day)
        ! Keith & Kreider 1978 Principles of Solar Engineering
        df = (1.0045 + 0.04349*akt - 3.5227*akt**2 + 2.6313*akt**3)*st
        if (akt > 0.75) df = 0.166*st
        if (df < 0.0) df = 0.0

        ! Calculate horizontal surface direct beam radiation (cal/cm2/day)
        dr = st - df

        ! Calculate direct beam tilt factor
        if (altsum /= 0.0) then
            Rb = angsum/altsum
        else
            Rb = 0.0
        end if

        ! Calculate diffuse radiation tilt factor
        Rd = (cos(rslope/2.0))**2

        ! Calculate surface solar radiation (cal/cm2/day)
        sun = Rb*dr + Rd*df

    end subroutine solar_rad

    !:.........................................................................:

    real function pot_evap(ta, sun, elev, e2, e1)
        !
        !  Calculates daily potential evapotranspiration (cm)
        !  Uses modified Priestley-Taylor equation
        !  Adapted from Bonan 1989 Ecological Modelling 45:275-306
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    10/10/17     A. C. Foster         Original Code
        !

        ! Data dictionary: calling arguments
        real, intent(in) :: ta   ! Air temperature (degC)
        real, intent(in) :: sun  ! Surface solar radiation (cal/cm2/day)
        real, intent(in) :: elev ! Elevation (m)
        real, intent(in) :: e1   ! Saturation vapor pressure (mbar) at tmin of warmest month
        real, intent(in) :: e2   ! Saturation vapor pressure (mbar) at tmax of warmest month

        ! Data dictionary: local variables
        real :: a, b ! Coefficients
        real :: hvap ! Latent heat of vaporization (kca/kg)
        real :: evap ! Local variable for evaporation (cm)

        ! Calculate latent heat of vaporization (kcal/kg)
        hvap = 597.391 - 0.5680*ta

        ! Calculate coefficients
        ! Campbell 1977 An Introduction to Environmental Biophysics
        a = 1.0/(38.0 - (2.0*elev/305.0) + 380.0/(e2 - e1))
        b = -2.5  - 0.14*(e2 - e1) - elev/550.0

        ! Calculate PET (cm)
        if (ta > 0.0) then
            evap = a*(ta - b)*sun/hvap
            evap = max(evap, 0.0)
        else
            ! Assume 0.0 PET if ta < 0.0
            evap = 0.0
        end if

        pot_evap = evap

    end function pot_evap

    !:.........................................................................:

    subroutine calc_ffmc(temp, RH, W_ms, precip_cm, Fo, ffmc)
        !
        !  Calculates the daily Canadian Forest Fire Rating System (CFRS) fine
        !  fuel moisture code (ffmc).
        !
        !  Adapted from Wang 2015 NRC Information Report NOR-X-424
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/15/20     A. C. Foster         Original Code
        !

        ! Data dictionary: calling arguments
        real,  intent(in)  :: temp      ! Noon daily air temperature (degC)
        real,  intent(in)  :: RH        ! Noon relative humidity (%)
        real,  intent(in)  :: W_ms      ! Noon wind speed (m/s)
        real,  intent(in)  :: precip_cm ! Precipitation over last 24 hrs (cm)
        real,  intent(in)  :: Fo        ! Previous day's ffmc
        real,  intent(out) :: ffmc      ! Fine fuel moisture code

        ! Data dictionary: local variables
        real :: precip ! Precipitation (mm)
        real :: ws     ! Wind speed (km/hr)
        real :: wmo    ! Previous day's fine fuel moisture (%)
        real :: ra     ! Effective rainfall (mm)
        real :: ed     ! Equilibrium moisture content from drying (%)
        real :: ew     ! Equilibrium moisture content from wetting (%)
        real :: z      ! Log drying/wetting rate at 21.1 degC
        real :: x      ! Actual drying/wetting rate
        real :: wm     ! Fine fuel moisture (%)

        precip = precip_cm*10.0    ! Convert rainfall to mm
        ws = W_ms/1000.0*60.0*60.0 ! Convert wind to km/hr

        ! Fine fuel moisture content from previous day
        wmo = (147.2*(101.0 - Fo))/(59.5 + Fo)

        ! Rain reduction to allow for loss in overhead canopy
        if (precip > 0.5) then
            ra = precip - 0.5
        else
            ra = 0.5
        end if

        ! Moisture content from wetting
        if (precip > 0.5) then
            if (wmo > 150.0) then
                wmo = wmo + 0.0015*(wmo - 150.0)*(wmo - 150.0)*                &
                    sqrt(ra) + 42.5*ra*exp(-100.0/(251.0 - wmo))*              &
                    (1.0 - exp(-6.93/ra))
            else
                wmo = wmo + 42.5*ra*exp(-100.0/(251.0 - wmo))*                 &
                    (1.0 - exp(-6.93)/ra)
            end if
        end if

        ! Cap wmo at 250%
        if (wmo > 250.0) wmo = 250.0

        ! Equilibrium moisture content from drying
          ed = 0.942*(RH**0.679) + (11.0 * exp((RH - 100.0)/10.0)) + 0.18*     &
            (21.1 - temp)*(1.0 - 1.0/exp(RH*0.115))

        ! Equilibrium moisture content from wetting
        ew = 0.618*(RH**0.753) + (10.0*exp((RH - 100.0)/10)) + 0.18*           &
            (21.1 - temp)*(1.0 - 1.0/exp(RH*0.115))

        ! Log drying rate at normal temperature (21.1degC)
        if (wmo < ed .and. wmo < ew) then
            z = 0.424*(1.0 - (((100.0 - RH)/100.0)**1.7)) + 0.0694*sqrt(ws)*   &
                (1.0 - ((100.0 - RH)/100.0)**8)
        else
            z = 0.0
        end if

        ! Effect of temperature on drying rate
        x = z*0.581*exp(0.0365*temp)

        ! Calculate moisture
        if (wmo < ed .and. wmo < ew) then
            wm = ew - (ew - wmo)/(10.0**x)
        else
            wm = wmo
        end if

        ! Log of wetting rate at normal temperature of 21.1degC
        if (wmo > ed) then
            z = 0.424*(1.0 - (RH/100.0)**1.7) + 0.0694*sqrt(ws)*               &
                (1.0 - (RH/100.0)**8.0)
        end if

        ! Effect of temperature on  wetting rate
        x = z*0.581*exp(0.0365*temp)

        ! Calculate moisture
        if (wmo > ed) then
            wm = ed + (wmo - ed)/(10.0**x)
        end if

        ! Calculate ffmc and correct for outside bounds
        ffmc = (59.5*(250.0 - wm))/(147.2 + wm)

        if (ffmc > 101) ffmc = 101
        if (ffmc <= 0.0) ffmc = 0.0

    end subroutine calc_ffmc

    !:.........................................................................:

    subroutine calc_dmc(temp, RH, precip_cm, Po, i, latitude, dmc)
        !
        !  Calculates the Canadian Forest Fire Rating System (CFRS) duff
        !  moisture code (dmc)
        !
        !  Adapted from Wang 2015 NRC Information Report NOR-X-424
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/15/20     A. C. Foster         Original Code
        !

        ! Data dictionary: constants

        ! Day length adjustments
        ! For latitude near equation (-10, 10 degrees), use a factor of 9
        ! for all months
        real, dimension(12), parameter :: ELL01 = [6.5, 7.5, 9.0, 12.8, 13.9,  &
            13.9, 12.4, 10.9, 9.4, 8.0, 7.0, 6.0] ! latititude >= 30 N
        real, dimension(12), parameter :: ELL02 = [7.9, 8.4, 8.9, 9.5, 9.9,    &
            10.2, 10.1, 9.7, 9.1, 8.6, 8.1, 7.8] ! 30 > latitude >= 10
        real, dimension(12), parameter :: ELL03 = [10.1, 9.6, 9.1, 8.5, 8.1,   &
            7.8, 7.9, 8.3, 8.9, 9.4, 9.9, 10.2] ! -10 > latitude >= -30
        real, dimension(12), parameter :: ELL04 = [11.5, 10.5, 9.2, 7.9, 6.8,  &
            6.2, 6.5, 7.4, 8.7, 10.0, 11.2, 11.8] ! latitude < -30

        ! Data dictionary: calling arguments
        real,    intent(in)  :: temp      ! Noon air temperature (degC)
        real,    intent(in)  :: RH        ! Noon elative humidity (%)
        real,    intent(in)  :: precip_cm ! Precipitation over last 24 hrs (cm)
        real,    intent(in)  :: Po        ! Previous day's dmc
        real,    intent(in)  :: latitude  ! Site latitude (degrees)
        integer, intent(in)  :: i         ! Month of simulation
        real,    intent(out) :: dmc       ! Duff moisture code

        ! Data dictionary: local variables
        real, dimension(12) :: ELL    ! Day length adjustment
        real                :: precip ! Precipitation (mm)
        real                :: rw     ! Net rainfall (mm)
        real                :: b      ! Temporary variable to calculate moisture content after rain
        real                :: temp0  ! Corrected temperature (degC)
        real                :: wmi    ! Previous day's duff moisture (%)
        real                :: wmr    ! Moisture content after rain (%)
        real                :: pr     ! Corrected rainfall (mm)
        real                :: rk     ! Log drying rate

        precip = precip_cm*10.0 ! Convert rain to mm

        ! Constrain low end of temperature
        if (temp < -1.1) then
            temp0 = -1.1
        else
            temp0 = temp
        end if

        ! Determine day length adjustment based on latitude
        if (latitude > 30.0) then
            ELL = ELL01
        else if (latitude <= 30.0 .and. latitude > 10.0) then
            ELL = ELL02
        else if (latitude <= 10.0 .and. latitude > -10.0) then
            ELL = 9.0
        else if (latitude <= -10.0 .and. latitude > -30.0) then
            ELL = ELL03
        else if (latitude <= -30.0) then
            ELL = ELL04
        end if

        ! Log drying rate
        rk = 1.894*(temp0 + 1.1)*(100.0 - RH)*(ELL(i)*0.0001)

        ! Net rainfall (mm)
        rw = 0.92*precip - 1.27

        ! Alteration to Eq. 12 to calculate more accurately
        wmi = 20.0 + 280.0/exp(0.023*Po)

        ! Eqs 13a-c
        if (Po <= 33.0) then
            b = 100.0/(0.5 + 0.3*Po)
        else if (Po > 33.0 .and. Po <= 65.0) then
            b = 14.0 - 1.3*log(Po)
        else
            b = 6.2*log(Po) - 17.2
        end if

        ! Moisture content after rain
        wmr = wmi + 1000.0*rw/(48.77 + b*rw)

        ! Constrain p
        if (precip <= 1.5) then
            pr = Po
        else
            pr = 43.43*(5.6348 - log(wmr - 20.0))
        end if

        if (pr < 0.0) pr = 0.0

        ! Calculate dmc
        dmc = pr + rk
        if (dmc < 0.0) dmc = 0.0

    end subroutine calc_dmc

    !:.........................................................................:

    subroutine calc_dc(temp, precip_cm, Dd, latitude, i, dc)
        !
        !  Calculates the Canadian Forest Fire Rating System (CFRS) drought
        !  code (dc)
        !  Adapted from Wang 2015 NRC Information Report NOR-X-424
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/15/20     A. C. Foster         Original Code
        !

        ! Data dictionary: constants
        ! Day length adjustment
        ! Near equator, just use 1.4 for all months
        real, dimension(12), parameter :: FL01 = [-1.6, -1.6, -1.6, 0.9, 3.8,  &
            5.8, 6.4, 5.0, 2.4, 0.4, -1.6, -1.6] ! latitude > 20 N
        real, dimension(12), parameter :: FL02 = [6.4, 5.0, 2.4, 0.4, -1.6,    &
            -1.6, -1.6, -1.6, -1.6, 0.9, 3.8, 5.8] ! latitude < -20

        ! Data dictionary: calling arguments
        real,    intent(in)  :: temp      ! Max air temperature (degC)
        real,    intent(in)  :: precip_cm ! Precipitation (cm)
        real,    intent(in)  :: Dd        ! Previous day's drought code
        real,    intent(in)  :: latitude  ! Site latitude (degrees)
        integer, intent(in)  :: i         ! Month of simulation
        real,    intent(out) :: dc        ! Drought code

        ! Data dictionary: local variables
        real :: precip ! Precipitation (mm)
        real :: temp0  ! Corrected temperature (degC)
        real :: pe     ! Potential evapotranspiration (mm)
        real :: ra     ! Precipitation (mm)
        real :: rw     ! Effective rainfall (mm)
        real :: smi    ! Temporary variable for calculating dc
        real :: dr0    ! Temporary variable for calculating drying rate
        real :: dr     ! Drying rate

        precip = precip_cm*10 ! Convert rainfall to mm

        ! Constrain temperature
        if (temp < -2.8) then
            temp0 = -2.8
        else
            temp0 = temp
        end if

        ! Potential evapotranspiration (mm)
        if (latitude > 20.0) then
            pe = (0.36*(temp + 2.8) + FL01(i))/2.0
        else if (latitude <= -20.0) then
            pe = (0.36*(temp + 2.8) + FL02(i))/2.0
        else
            pe = (0.36*(temp + 2.8) + 1.4)/2.0
        end if

        ! Cap pe at 0 for negative winter DC values
        if (pe < 0.0) pe = 0.0
        ra = precip

        ! Effective rainfall
        rw = 0.83*ra - 1.27

        ! Eq. 19
        smi = 800.0*exp(-1.0*Dd/400.0)

        ! Alteration to Eq 21
        dr0 = Dd - 400.0*log(1.0 + 3.937*rw/smi)
        if (dr0 < 0.0) dr0 = 0.0

        ! Drying rate, use yesterday's DC if precip < 2.8
        if (precip <= 2.8) then
            dr = Dd
        else
            dr = dr0
        end if

        ! Calculate dc
        dc = dr + pe

        if (dc < 0.0) dc = 0.0

    end subroutine calc_dc

    !:.........................................................................:

    real function hrsabove(tmin, tmax, hrise, thresh)
        !
        !  Calculates cumulative hours above a threshold for the day
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/01/16     A. C. Foster         Original Code
        !    04/30/10     A. C. Foster         Added input hour of sunrise
        !

        ! Data dictionary: constants
        integer, parameter :: HMAX = 14 ! Hour of maximum temperature

        ! Data dictionary: calling arguments
        real,    intent(in) :: tmin   ! Minimum temperature (degC)
        real,    intent(in) :: tmax   ! Maximum temperature (degC)
        real,    intent(in) :: thresh ! Threshold temperature (degC)
        integer, intent(in) :: hrise  ! Hour of sunrise

        ! Data dictionary: local varuables
        real,    dimension(24) :: th    ! Hourly temperature (degC)
        real                   :: tmean ! Mean daily temperature (degC)
        integer                :: i     ! Looping index
        integer                :: couni ! Local count for hours above

        ! Calculate average temperature
        tmean = (tmax + tmin)/2.0

        ! Calculate hourly temperature for the day
        do i = 0, 23
            if ((i >= 0) .and. (i < hrise)) then
                th(i+1) = tmean + ((tmax - tmin)/2.0)*                         &
                    cos((pi*(float(i) + 10.0))/(10.0 + hrise))
            else if ((i >= hrise) .and. (i <= HMAX)) then
                th(i+1) = tmean - ((tmax - tmin)/2.0)*                         &
                    cos((pi*(float(i) - hrise))/(HMAX - hrise))
            else if ((i > HMAX) .and. (i < 24)) then
                th(i+1) = tmean + ((tmax - tmin)/2.0)*                         &
                    cos((pi*(float(i) + HMAX))/(10.0 + hrise))
            end if
        end do

        ! Accumulate hours above the threshold temperature
        couni = 0
        do i = 1, 24
            if (th(i) > thresh) then
                couni = couni + 1
            endif
        enddo

        hrsabove = float(couni)

    end function hrsabove

    !:.........................................................................:

    real function esat(ta)
        !
        !  Calculates saturation vapor pressure (mbar)
        !  Adapted from Campbell 1977 An Introduction to Biophysics
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !     10/10/17     A. C. Foster         Original Code
        !

        ! Data dictionary: calling arguments
        real,    intent(in) :: ta ! Air temperature (degC)

        esat = 33.8639*((0.00738*ta + 0.8072)**8.0 -                           &
            0.000019*abs(1.8*ta + 48.0) + 0.001316)

    end function esat

    !:.........................................................................:

end module Climate
