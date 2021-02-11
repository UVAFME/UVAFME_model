module Soil

!*******************************************************************************
  !
  ! The Soil module contains attributes and procedures relevant to soil
  !  properties.
  !
!*******************************************************************************

    use Constants
    use Parameters
    use Random
    use FileUtils
    use csv_file

    implicit none

    ! Data dictionary: global constants
    real,    parameter :: LAI_MIN = 0.01   ! Scalar for calculating min canopy moisture
    real,    parameter :: LAI_MAX = 0.15   ! Scalar for calculating max canopy moisture
    real,    parameter :: TFC1 = 0.8
    real,    parameter :: TFC2 = 0.2
    real,    parameter :: BFC = 0.4
    integer, parameter :: NCOH_MAX = 1500  ! Maximum number of decaying cohorts
    integer, parameter :: NCOH_CHR = 13    ! Number of cohort characteristics
    integer, parameter :: IMOSS = 20       ! Index location of moss in soil%litter array
    integer, parameter :: IWDW = 19        ! Index location of well-decayed wood in soil%litter array

    ! Define soil type
    type SoilData
        real, dimension(NCOH_MAX, NCOH_CHR) :: cohorts       ! Array of decaying cohorts of litter
        real, dimension(LIT_LEVS)           :: litter        ! Fresh litter content (t/ha)
        real, dimension(LIT_LEVS+1, 2)      :: forest_litter ! Total forest litter (t/ha)
        real, dimension(2)                  :: z_drain       ! Drainage capacity (volumetric)
        real, dimension(2)                  :: sat           ! Saturation capacity (volumetric)
        real, dimension(2)                  :: fc            ! Field capacity (volumetric)
        real, dimension(2)                  :: pwp           ! Permanent wilting point (volumetric)
        real, dimension(2)                  :: wc            ! Gravimetric water content
        real, dimension(2)                  :: H2Oice        ! Ice water content (m)
        real, dimension(2)                  :: water         ! Liquid water content (m)
        real, dimension(2)                  :: d_melt        ! Depth of thaw (m)
        real, dimension(2)                  :: d_freeze      ! Depth of freezing (m)
        real                                :: lai_w0        ! Canopy moisture (cm)
        real                                :: N_used        ! N used by plant growth (tN/ha)
        real                                :: avail_N       ! Plant-available N (tN/ha)
        real                                :: runoff        ! Runoff (cm)
        real                                :: active        ! Seasonal maximum thaw depth (m)
        real                                :: z_freeze      ! Seasonal maximum freezing depth (m)
        real                                :: A_depth       ! A-layer depth (m)
        real                                :: O_depth       ! Organic layer depth (m)
        real                                :: O_bulk_dens   ! Organic layer bulk density (kg/m3)
        real                                :: A_bulk_dens   ! A layer bulk density (kg/m3)
        real                                :: snowpack      ! Snowpack depth (m)
        real                                :: swe           ! Snowpack snow-water equivalent (m)
        real                                :: fan           ! N volatilized from frecent fires
        real                                :: moss_biom     ! Moss biomass (kg)
        real                                :: M_depth       ! Moss depth (m)
        real                                :: dec_fuel      ! Fresh deciduous leaf litter (t/ha)
        real                                :: fol_fuel      ! Leaf litter fuel (t/ha)
        real                                :: twig_fuel     ! Branch litter fuel (t/ha)
        real                                :: smbl_fuel     ! Small bole fuel (t/ha)
        real                                :: lrbl_fuel     ! Large bole fuel (t/ha)
        real                                :: avail_fuel    ! Available fuel for burning (t/ha)
        real                                :: minWC         ! Minimum gravimetric water content for the year
        integer                             :: itxt          ! Soil texture (0: very coarse; 1: coarse; 2: fine)
        integer                             :: ncohort       ! Number of decaying cohorts
    end type SoilData


  contains

    !:.........................................................................:

    subroutine permf(soil, month, knd, alff, tdd, fdd, cfs, zdepth)
        !
        !  Calculates depth of freeze/thaw in the soil profile (zdepth, m)
        !  Adapted from Bonan 1989 Ecological Modelling 45:275-306
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    5/15/17     A. C. Foster        Original Code
        !

        ! Data dictionary: calling arguments
        class(SoilData),            intent(inout) :: soil   ! Soil object
        real, dimension(NTEMPS, 2), intent(in)    :: tdd    ! Thawing degree-days (>0degC)
        real, dimension(NTEMPS, 2), intent(in)    :: fdd    ! Freezing degree-days (<0degC)
        real,                       intent(in)    :: alff   ! Available light on the forest floor (0-1)
        real,                       intent(in)    :: cfs    ! Ratio of actual:horizontal surface radiation (0-1)
        integer,                    intent(in)    :: month  ! Month of simulation
        integer,                    intent(in)    :: knd    ! 1:freezing; 2: thawing
        real,                       intent(out)   :: zdepth ! Total depth of freeze/thaw (m)

        ! Data dictionary: local variables
        real    :: tyr   ! Sum of thermal resistance of above layers (m2degChr/kcal)
        real    :: cft   ! Correction factor for thawing degree-days
        real    :: cff   ! Correction factor for freezing deagree-days
        real    :: cadd  ! Avaiable thawing/freezing degree-days
        real    :: dl    ! Soil depth (m)
        real    :: vwc   ! Volumetric moisture content
        real    :: tku   ! Unfrozen thermal conductivity (kcal/m2/degC/hr)
        real    :: tkf   ! Frozen thermal conductivity (kcal/m2/degC/hr)
        real    :: gwc   ! Gravimetric moisture content
        real    :: bde   ! Soil bulk density (lb/ft3)
        real    :: tk    ! Thermal conductivity (kcal/m2/degC/hr)
        real    :: ql    ! Latent heat of fusion (kcal/m3)
        real    :: res   ! Thermal resistance (m2degChr/kcal)
        real    :: degd  ! Required degree-days to thaw entire soil layer
        real    :: depth ! Freeze/thaw depth (m)
        real    :: ax    ! Temporary variable for calculating partial freeze/thaw depth
        real    :: bx    ! Temporary variable for calculating partial freeze/thaw depth
        real    :: cx    ! Temporary variable for calculating partial freeze/thaw depth
        real    :: Ad    ! Mineral layer depth (m)
        real    :: Od    ! Organic layer depth (m)
        real    :: ABD   ! Mineral layer bulk density (kg/m3)
        real    :: OBD   ! Organic layer bulk density (kg/m3)
        integer :: l     ! Looping index

        ! Convenience variables to reduce table lookups
        Ad = soil%A_depth
        Od = soil%O_depth
        ABD = soil%A_bulk_dens
        OBD = soil%O_bulk_dens

        ! Initialize zdepth
        zdepth = 0.0

        ! Tyr is initially just moss
        if (soil%M_depth .gt. 0) then
            tyr = soil%M_depth/(0.291*0.86)
        else
            tyr = 0.0
        end if

        ! Adjust freeze/thaw days for surface conditions
        if (alff <= 0.50) then
            cft = 0.62
            cff = 0.38
        else if (alff > 0.50 .and. alff <= 0.75) then
            cft = 0.77
            cff = 0.37
        else if (alff > 0.75) then
            cft = 0.92
              cff = 0.36
        end if

        ! Calculate actual freezing/thawing degree-days
        if (knd .eq. 1) then
            cadd = fdd(month, 2)*cff*(2.0 - cfs)
        else
            cadd = tdd(month, 2)*cft*cfs
        end if

        ! Calculate freeze/thaw of soil layers
        do l = 1, 2
            if (l .eq. 1) then
                dl = Od
            else
                dl = Ad  + 2.0
            end if

            if (l .eq. 1) then
                ! Moss/humus thermal layer conductivities
                vwc = soil%wc(l)*OBD/H2O_D ! Volumetric moisture content
                tku = (0.5*(soil%pwp(l) - vwc) + 0.08*(vwc - soil%fc(l)))/     &
                    (soil%pwp(l) - soil%fc(l))
                tkf = (2.0*tku*(soil%pwp(l) - vwc) + tku*(vwc -                &
                    soil%fc(l)))/(soil%pwp(l) - soil%fc(l))

            else
                ! Mineral soil thermal conductivities
                gwc = soil%wc(2) ! Gravimetric moisture content
                bde = ABD*0.06243 ! Bulk density in lb/ft3
                if (soil%itxt .le. 1 ) then ! Granular
                    tku = (0.7*log10(gwc*100.0) + 0.4)*10.0**(0.01*bde)
                    tkf = 0.076*10.0**(0.013*bde) + (gwc*100.0)*0.032*         &
                        10.0**(0.0146*bde)
                else ! Fine-textured
                    tku = (0.9*log10(gwc*100.0) - 0.2)*10.0**(0.01*bde)
                    tkf = 0.01*10.0**(0.022*bde) + (gwc*100.0)*0.085*          &
                        10.0**(0.008*bde)
                end if

                ! Convert from btu-in/ft2/hr/f to kcal/m2/hr/c
                tkf = tkf/8.0645
                tku = tku/8.0645

            end if

            ! Determine to use freezing or thawing conductivity
            if (knd .eq. 1) then
                tk = tkf
            else
                tk = tku
            end if

            ! Latent heat of fusion (kcal/m3)
            if (l .eq. 1) then
                ql = 80.0*(soil%wc(l))*OBD
            else
                ql = 80.0*(soil%wc(l))*ABD
            end if

            ! Calculate thermal resistance
            res = dl/tk

            ! Calculate degree-days required to freeze entire layer
            degd = ql*dl/24.0*(tyr + res/2.0)

            if (degd .le. cadd) then
                ! We have plenty - freeze/thaw whole layer
                depth = dl
                cadd = cadd - degd
            else
                ! Not enough - calculate actual freeze/thaw depth
                ax = 0.5*ql/tk
                bx = ql*tyr
                cx = -24.0*cadd
                depth = (-bx + sqrt(bx*bx - 4.0*ax*cx))/(2.0*ax)
                depth = min(depth, dl)
                depth = max(depth, 0.0)
                cadd = 0.0
            end if

            ! Update resistance
            tyr = tyr + res

            ! Update depth
            zdepth = zdepth + depth

        end do

    end subroutine permf

    !:.........................................................................:

    subroutine moist(soil, siteID, ip, year, day, ta, rain, pot_ev_day, lai ,  &
        slope, amlt, xmlt, xfrz, tdd, k, aet, aow0_ByMin, saw0_ByFC,           &
        saw0_ByWP, saw0_BySAT)
        !
        !  Calculates daily soil moisture dynamics
        !  Adapted from Bonan 1989 Ecological Modelling 45:275-306
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    5/16/17     A. C. Foster        Original Code
        !

        ! Data dictionary: constants
        real, parameter :: MELTFACT = 4.0 ! Scalar for melting of snowpack
        real, parameter :: TC = 0.0       ! Minimum temperature for snow melting (degC)
        real, parameter :: TSMAX = 3.3    ! Minimum temperature for rainfall (degC)
        real, parameter :: TSMIN = -1.1   ! Maximum temperature for snowfall (degC)

        ! Scalar for how quickly moisture drains from field capacity
        ! Depends on soil texture
        real, parameter, dimension(3) :: EXS_P = [4.0, 1.0, 0.6]

        ! Data dictionary:, calling arguments
        class(SoilData),            intent(inout) :: soil       ! Soil object
        real, dimension(NTEMPS, 2), intent(in)    :: tdd        ! Thawing degree-days (>0degC)
        integer,                    intent(in)    :: siteID     ! Site ID
        integer,                    intent(in)    :: ip         ! Plot number
        integer,                    intent(in)    :: year       ! Year
        integer,                    intent(in)    :: k          ! Month of simulation
        integer,                    intent(in)    :: day        ! Julian day
        real,                       intent(in)    :: ta         ! Air temperature (degC)
        real,                       intent(in)    :: rain       ! Precipitation (cm)
        real,                       intent(in)    :: pot_ev_day ! Potential evapotranspiration (cm)
        real,                       intent(in)    :: amlt       ! Previous year's seasonal maximum depth of thaw (m)
        real,                       intent(in)    :: xmlt       ! Current depth of thaw (m)
        real,                       intent(in)    :: xfrz       ! Current freezing depth
        real,                       intent(in)    :: lai        ! Leaf area index (m2/m2)
        real,                       intent(in)    :: slope      ! Slope (degrees)
        real,                       intent(out)   :: aet        ! Actual evapotranspiration (cm)
        real,                       intent(out)   :: aow0_ByMin ! Organic layer moisture scaled by wilting point
        real,                       intent(out)   :: saw0_ByFC  ! Mineral layer moisture scaled by field capacity
        real,                       intent(out)   :: saw0_ByWP  ! Mineral layer moisture scaled by wilting point
        real,                       intent(out)   :: saw0_BySAT ! Mineral layer moisture scaled by saturation capacity

        ! Data dictionary: local variables
        real, dimension(2) :: zh         ! Soil layer depth (m)
        real, dimension(2) :: dmelt      ! Current thaw depth (m)
        real, dimension(2) :: odmelt     ! Previous day's thaw depth (m)
        real, dimension(2) :: dfreeze    ! Depth of freezing (m)
        real, dimension(2) :: odfreeze   ! Previous day's freezing depth (m)
        real, dimension(2) :: zdrain     ! Drainage capacity (volumetric)
        real, dimension(2) :: water      ! Liquid moisture content (m)
        real, dimension(2) :: H2Oice     ! Frozen moisture content (m)
        real, dimension(2) :: wc         ! Gravimetric moisture content
        real, dimension(2) :: wcf        ! Liquid water available for freezing (m)
        real, dimension(2) :: aet_loss   ! Evaporative loss from soil (m)
        real, dimension(2) :: gwl        ! Factor for how quickly water drains from field capacity
        real, dimension(2) :: wt         ! Thawing water content (m)
        real, dimension(2) :: pwp        ! Permanent wilting point scaled by thaw depth (m)
        real, dimension(2) :: fc         ! Field capacity scaled by thaw depth (m)
        real, dimension(2) :: sat        ! Saturation capacity scaled by thaw depth (m)
        real, dimension(2) :: wf         ! Freezing water content (m)
        real, dimension(2) :: exs        ! Excess moisture above field capacity (m)
        real, dimension(2) :: standing   ! Excess moisture above saturation capacity (m)
        real               :: Od         ! Organic layer depth (m)
        real               :: Ad         ! Mineral layer depth (m)
        real               :: OBD        ! Organic layer bulk density (kg/m3)
        real               :: ABD        ! A layer bulk density (kg/m3)
        real               :: can_evap   ! Canopy evaporation (m)
         real               :: canopy_int ! Canopy interception (m)
        real               :: laic       ! LAI (m2/m2)
        real               :: laiw_min   ! Minimum canopy moisture content (m)
        real               :: laiw_max   ! Maximum canopy moisture content (m)
        real               :: laiw       ! Canopy moisture content (m)
        real               :: xabove     ! Depth of above soil layers (m)
        real               :: tfall      ! Throughfall (m)
        real               :: slope_fact ! Proportion of surface water than runs off
        real               :: loss_slp   ! Surface water lost to slope runoff (m)
        real               :: precip     ! Precipitation (m)
        real               :: pet        ! Potential evapotranspiration (m)
        real               :: pr         ! Liquid precipitation (m)
        real               :: ps         ! Snowfall water equivalent (m)
        real               :: melt       ! Snowpack melt (m)
        real               :: pwl        ! Potential water losses/gains (m)
        real               :: rz         ! Relative rooting density in organic layer
        real               :: dthaw      ! Depth of thaw in current layer (m)
        real               :: pwll       ! Atmospheric demand for current layer (m)
        real               :: pwg        ! Water gain from throughfall, snowmelt (m)
        real               :: B          ! Temporary variable for calcualating aet loss from soil
        integer            :: l          ! Looping index

        ! Convenience varaibles
        dmelt = soil%d_melt
        dfreeze = soil%d_freeze
        odmelt = soil%d_melt
        odfreeze = soil%d_freeze
        Od = soil%O_depth
        Ad = soil%A_depth
        OBD = soil%O_bulk_dens
        ABD = soil%A_bulk_dens
        zdrain = soil%z_drain
        water = soil%water
        H2Oice = soil%H2Oice
        wc = soil%wc

        ! Initialize variables
        can_evap = 0.0
        canopy_int = 0.0
        xabove = 0.0
        aet = 0.0
        aet_loss = 0.0
        exs = 0.0

        ! Calculate LAI metrics
        laic = max(lai, 1.0)
        laiw_min = (laic*LAI_MIN)*CM_TO_M
        laiw_max = (laic*LAI_MAX)*CM_TO_M
        laiw = soil%lai_w0*CM_TO_M

        ! Set zh to appropriate soil layer depths
        zh(1) = soil%M_depth + Od
        zh(2) = Ad

        ! Convert precipitation to m
        precip = rain*CM_TO_M

        ! Convert PET to m
        pet = pot_ev_day*CM_TO_M

        ! Partition precipitation between rain and snow
        if (ta .ge. TSMAX) then

            ! All liquid rainfall
            pr = precip
            ps = 0.0

            ! Calculate canopy interception and throughfall
            canopy_int = min(max((laiw_max - laiw), 0.0), pr)
            laiw = laiw + canopy_int
            tfall = max(pr - canopy_int, 0.0)

        else if (ta .lt. TSMAX .and. ta .gt. TSMIN) then

            ! Mixture of snow and rain
            ps = (TSMAX - ta)/(TSMAX - TSMIN)*precip
            pr = precip - ps

            ! Accumulate snowpack
            soil%swe = soil%swe + ps

            ! Calculate canopy interception and throughfall
            ! For now we only intercept liquid rainfall
            canopy_int = min(max((laiw_max - laiw), 0.0), pr)
            laiw = laiw + canopy_int
            tfall = max(pr - canopy_int, 0.0)

        else if (ta .le. TSMIN) then

            ! Only snow - no canopy interception or throughfall
            pr = 0.0
            ps = precip
            canopy_int = 0.0
            tfall = 0.0

            ! Accumulate snowpack
            soil%swe = soil%swe + ps

        end if

        ! Melt snowpack if it exists and if ta is above tc
        if (soil%swe .ge. epsilon(1.0)) then
            if (ta .gt. TC) then

                ! Melt the snowpack
                melt = (MELTFACT*(ta - TC))*0.001

                if (soil%swe .ge. melt) then

                    ! Melt some snowpack
                    soil%swe = max(soil%swe - melt, 0.0)
                else

                    ! Melt whole snowpack - none left
                    melt = soil%swe
                    soil%swe = 0.0

                end if
            else
                ! Temperature too cold to melt
                melt = 0.0
            end if
        else
            ! No snowpack to melt
            melt = 0.0
        end if

        ! Calculate snow density
        soil%snowpack = soil%swe/(SNOW_DENS/H2O_D)

        ! Add up throughfall, snowpack melt
        tfall = tfall + melt

        ! Calculate losses due to slope runoff
        slope_fact = (slope/90.0)**2
        loss_slp = slope_fact*tfall

        ! Calculate potential water loss or gain this day
        pwl = tfall - loss_slp - pet

        ! Relative root density in forest floor organic layer
        if (zh(1) .gt. 0.0) then
            rz = min((zh(1) + amlt), 1.0)
            rz = 2.0*zh(1)/rz*(1.0 - zh(1)/(2.0*rz))
        else
            rz = 0.0
        end if

        ! Water released in soil thawing
        do l = 1, 2

            ! Depth of thaw in current layer
            dthaw = max(0.0, (xmlt - xabove))

            ! Cap to depth of layer
            dmelt(l) = min(zh(l), dthaw)

            ! Cap to max of today's and previous day's
            dmelt(l) = max(dmelt(l), odmelt(l))
            wt(l) = zdrain(l)*max(dmelt(l) - odmelt(l), 0.0)

            ! Check to make sure enough frozen water for thawing amount
            if ((H2Oice(l) - wt(l)) .ge. epsilon(1.0)) then

                ! Plenty of ice for melting - add thawed content to liquid
                ! water pool, and remove from ice content
                water(l) = water(l) + wt(l)
                H2Oice(l) = H2Oice(l) - wt(l)
            else

                ! Not enough ice - add whole amount to liquid content and 0.0
                ! the ice content
                water(l) = water(l) + H2Oice(l)
                H2Oice(l) = 0.0
            endif

            ! Set moisture conditions
            sat(l) = soil%sat(l)*dmelt(l)
            fc(l) = zdrain(l)*dmelt(l)
            pwp(l) = soil%pwp(l)*dmelt(l)

            ! Update xabove
            xabove = xabove + zh(l)
        end do

        ! Add positive water balance to layer, adjusting for excess
          ! water from above layers
        do l = 1, 2
            if (l .eq. 1) then
                if (pwl .ge. 0.0) then
                    ! More water available than pet, aet = pet
                    aet = pet

                    ! Water gain is water from throughfall and snowmelt
                    pwg = pwl
                else
                    ! No excess water from rainfall (aet < pet)
                    pwg = 0.0
                end if
            else
                if (pwl .ge. 0.0) then
                    ! Excess water from above layer
                    pwg = exs(1)
                else
                    ! No excess water from rainfall/groundwater flow (aet < pet)
                    pwg = 0.0
                endif
            end if

            ! Add throughall/snowmelt/groundwater flow to layers
            water(l) = water(l) + pwg

            ! Calculate excess water
            exs(l) = max((water(l) - fc(l)), 0.0)
            standing(l) = max((water(l) - sat(l)), 0.0)

            ! Drain off excess water
            if (exs(l) .gt. 0.0) then
                gwl(l) = min(1.0, max(0.0, EXS_P(soil%itxt)*                   &
                    ((exs(l)**2)/(pet + exs(l))*(1 - sat(l))))/exs(l))
            end if
            exs(l) = standing(l) + (exs(l) - standing(l))*gwl(l)
            water(l) = max(water(l) - exs(l), 0.0)
        end do

        ! Calculate runoff
        soil%runoff = exs(2) + loss_slp

        if (pwl .lt. 0.0) then
            ! We have atmospheric demand - need to evaporate from canopy and
            ! soil layers

            ! Calculate canopy evaporation
            can_evap = min(-1.0*pwl, max(laiw - laiw_min, 0.0))
            laiw = max(laiw - can_evap, laiw_min)

            ! Update atmosperic demand and aet
            pwl = min(pwl + can_evap, 0.0)
            aet = aet + can_evap
        end if

        ! Partition negative water balance between layers based on root density
        do l = 1, 2
            if (l .eq. 1) then
                pwll = min(0.0, pwl*rz) ! Either negative or 0.0
            else
                pwll = min(0.0, pwl*(1.0-rz)) ! Either negative or 0.0
            end if

            if (dmelt(l) .gt. 0.0) then
                ! We have some unfrozen soil to take water from

                ! Calculate aet loss from soil
                B = 0.461 - 1.10559/(water(l)/zh(l))
                aet_loss(l) = min(water(l) - water(l)*exp(B*abs(pwll)),        &
                    -1.0*pwll)
                water(l) = max(water(l) - aet_loss(l), 0.0)

                ! Update aet
                aet = aet + aet_loss(l)
            else
                ! No unfrozen soil - nothing to evaporate
                aet_loss(l) = 0.0
            end if
        end do

        ! Water content of soil layer when able to freeze
        do l = 1, 2
            if (tdd(k,1) .gt. 0.0) then
                if (dmelt(l) .gt. 0.0) then

                    ! Some unfrozen water to freeze
                    wcf(l) = water(l)/dmelt(l)
                else
                    wcf(l) = 0.0
                end if
                odfreeze(l) = 0.0
            end if
        end do

        ! Water frozen as ice in soil freezing
        xabove = 0.0
        do l = 1, 2
            ! Calcualte actual depth of freezing in layer
            dfreeze(l) = max(min(zh(l), xfrz - xabove), 0.0)

            ! Calculate water available for frezing
            dfreeze(l) = max(dfreeze(l), odfreeze(l))
            wf(l) = wcf(l)*max(dfreeze(l) - odfreeze(l), 0.0)

            ! Check to make sure enough liquid water for freezing amount
            if ((water(l) - wf(l)) .ge. epsilon(1.0)) then

                ! Plenty of water to freeze - remove from liquid water pool and
                ! add to frozen
                water(l) = water(l) - wf(l)
                H2Oice(l) = H2Oice(l) + wf(l)
            else
                ! Not enough - add whole liquid pool to frozen and 0.0 liquid
                H2Oice(l) = H2Oice(l) + water(l)
                water(l) = 0.0
            end if

            ! Calculate gravimetric water content
            if (zh(l) .gt. 0.0) then
                if (l .eq. 1) then
                    wc(l) = ((water(l)+H2Oice(l))/zh(l)*H2O_D/OBD)
                else
                    wc(l) = ((water(l)+H2Oice(l))/zh(l)*H2O_D/ABD)
                end if
            else
                wc(l) = 0.0
            end if

            ! Adjust for above layers
            xabove = xabove + zh(l)
        end do

        ! Calculate scaled water content
        if (pwp(1) .gt. 0.0) then
            aow0_ByMin = water(1)/pwp(1)
        else
            aow0_ByMin = 0.0
        end if

        if (fc(2) .gt. 0.0) then
            saw0_ByFC = water(2)/fc(2)
        else
            saw0_ByFC = 1.0
        end if

        if (pwp(2) .gt. 0.0) then
            saw0_ByWP = water(2)/pwp(2)
        else
            saw0_ByWP = 2.0
        end if

        if (sat(2) .gt. 0.0) then
            saw0_BySAT = water(2)/sat(2)
        else
            saw0_BySAT = 0.0
        end if

        ! Update instance variables and outputs
        soil%d_freeze = dfreeze
        soil%d_melt = dmelt
        soil%water = water
        soil%H2Oice = H2Oice
        soil%wc = wc
        soil%lai_w0 = laiw*M_TO_CM
        soil%runoff = soil%runoff*M_TO_CM
        aet = aet*M_TO_CM

    end subroutine moist

    !:.........................................................................:

    subroutine moss(soil, alff, cla, decLit, drydays, siteID, plot, year)
        !
        !  Calculates annual moss growth and mortality
        !  Adapted from Bonan and Korzukhin 1989 Vegetatio 84:31-44
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    10/10/17     A. C. Foster        Original Code
        !

        ! Data dictionary: constants - from Bonan & Korzukhin (1989)
        real, parameter :: A1 = 3.41      ! Light curve parameter
        real, parameter :: A2 = 2.14      ! Light curve parameter
        real, parameter :: A3 = 0.08      ! Light compensation point
        real, parameter :: Q = 0.12       ! Mortality parameter
        real, parameter :: B1 = 0.136     ! Mortality parameter
        real, parameter :: EXT = 3.5      ! Light extinction coefficient
        real, parameter :: Q1 = 1.0       ! Specific leaf area
        real, parameter :: SPORES = 0.001 ! Proportion of biomass spent on reproduction
        real, parameter :: PMAX = 0.3     ! Maximum moss production (kg/m2)

        ! Data dictionary: calling arguments
        class (SoilData), intent(inout) :: soil    ! Soil object
        real,             intent(in)    :: alff    ! Available light on the forest floor (0-1)
        real,             intent(in)    :: cla     ! Cumulative leaf area on forest floor (m2)
        real,             intent(in)    :: decLit  ! Fresh deciduous leaf litter (t/ha)
        real,             intent(in)    :: drydays ! Drought index (0-1)
        integer,          intent(in)    :: siteID  ! Site ID
        integer,          intent(in)    :: plot    ! Plot number
        integer,          intent(in)    :: year    ! Year

        ! Data dictionary: local variables
        real :: biokg   ! Moss biomass (kg/m2)
        real :: biokg_1 ! Previous year's moss biomass (kg/m2)
        real :: al      ! Available light
        real :: algf    ! Available light growth factor
        real :: fcgf    ! Forest cover growth factor
        real :: dlgf    ! Deciduous leaf litter growth factor
        real :: ddgf    ! Moisture growth factor
        real :: assim   ! Moss assimilation rate (kg/kg)
        real :: binp    ! Respiration rate (kg/kg)
        real :: binc1   ! Potential moss growth/mortality (kg/m2)
        real :: binc    ! Actual moss growth/mortality (kg/m2)
        real :: xx      ! Moss litter (kg)
        real :: bnpp    ! Moss NPP (kg)

        ! Convert moss biomass in kg to kg/m2
        biokg = soil%moss_biom/plotsize

        ! Save old value of biomass (kg/m2)
        biokg_1 = biokg

        ! Light growth multiplier
        al = exp(-1.0*EXT*Q1*biokg/2.0)
        algf = A1*(al - A3)/(1.0 + A2*al)
        algf = max(0.0, algf)
        algf = min(1.0, algf)

        ! Forest cover growth multiplier (alff > 0.5)
        if (alff .gt. 0.5) then
            fcgf = 1.25 - alff**2.0
        else
            fcgf = 1.0
        end if

        ! Deciduous leaf litter growth multiplier
        if (decLit .gt. 0.0) then
            dlgf = exp(-0.45*decLit)
            if (dlgf .le. 0.0) dlgf = 0.0
            if (dlgf .ge. 1.0) dlgf = 1.0
        else
            dlgf = 1.0
        end if

        ! Soil moisture growth multiplier
        if (drydays .gt. 0.10) then
            ddgf = 0.0
        else
            ddgf = 1.0
        end if

        ! Calculate assimilation
        if (alff .le. 0.5) then
            assim = PMAX*(1.27 + 0.3*sqrt(alff))*algf*ddgf*fcgf*dlgf
        else
            assim = PMAX*algf*ddgf*fcgf*dlgf
        end if

        ! Calculate respiration rate (kg/kg)
        binp = SPORES*ddgf*dlgf

        ! Calculate moss growth/mortality
        binc1 = Q1*biokg*assim - Q1*biokg*(Q + B1) + binp

        if (biokg + binc1 .lt. 0.0) then
            ! Not enough moss to account for mortality/respiration
            ! Set moss biomass loss to all of current moss biomass
            binc = -1.0*biokg
        else
            ! Enough to account - set loss/gain to temporary value
            binc = binc1
        end if

        ! update moss biomass (kg)
        soil%moss_biom = (biokg + binc)*plotsize

        ! Calculate litter (kg)
        xx = (Q1*biokg*(assim - Q) + binp - binc)*plotsize
        if (xx .lt. 0.0) then
            soil%litter(IMOSS) = 0.0
        else
            soil%litter(IMOSS) = xx
        end if

        ! Calculate moss NPP (kg)
        bnpp = binc*plotsize + soil%litter(IMOSS)

        ! Convert to tonnes/ha from kg/plot
        soil%litter(IMOSS) = soil%litter(IMOSS)/plotsize*HEC_TO_M2*KG_TO_T

        ! Thickness of moss layer (m)
        soil%M_depth = soil%moss_biom/plotsize/BULK_MOSS


    end subroutine moss

    !:.........................................................................:

    subroutine soiln(soil, aet_mm, cla, soildays, flooddays, avail_n)
        !
        !  Calculates annual soil decomposition and plant-available N
        !  Adapted from Bonan 1990 Biogeochemistry 10:1-28 &
        !  Pastor & Post 1985 ORNL/TM-9519
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    10/10/17    A. C. Foster         Original Code
        !    05/01/20    A. C. Foster         Broke branch litter into twigs,
        !                                        small branches, large branches
        !

        ! Data dictionary: calling arguments
        class(SoilData),   intent(inout) :: soil      ! Soil object
        real,              intent(in)    :: aet_mm    ! Actual evapotranspirtion (mm)
        real,              intent(in)    :: cla       ! Cumulative leaf area at forest floor (m2)
        real,              intent(in)    :: soildays  ! Soil degree-days (>0degC)
        real,              intent(in)    :: flooddays ! Flooding index (0-1)
        real,              intent(out)   :: avail_n   ! Plant-available nitrogen (tN/ha)

        ! Data dictionary: local variables
        real, dimension(NCOH_MAX, NCOH_CHR) :: C            ! Local array of decaying cohorts
        real                                :: litCO2       ! CO2 release from litter decay
        real                                :: humCO2       ! CO2 release from humus decay
        real                                :: totCO2       ! Total CO2 release from decay
        real                                :: WDW_new      ! New cohort of well-decayed wood (t/ha)
        real                                :: perc_rem     ! Percent remaining of original cohort (0-1)
        real                                :: tot_Nimob    ! N immobilization (tN/ha)
        real                                :: smult        ! Decay factor due to moisture
        real                                :: lit_Nmin     ! N mineralization from litter decay (tN/ha)
        real                                :: hum_Nmin     ! N mineralization from humus decay (tN/ha)
        real                                :: tot_Nmin     ! Total N mineralization (tN/ha)
        real                                :: dec_lit      ! Fresh deciduous leaf litter (t/ha)
        real                                :: perc_wtloss  ! Percent weight loss from decay (0-1)
        real                                :: perc_wtloss2 ! Temporary variable for percent weight loss
        real                                :: cla_decaym   ! Decay factor due to canopy cover
        real                                :: tfallNmin    ! N mineralization from throughfall (tN/ha)
        real                                :: wtloss       ! Weight loss from decay (t/ha)
        real                                :: crit_perc    ! Critical percent remaining before transfer to humus (0-1)
        real                                :: delta_N      ! Absolute change in N content from decay (tN/ha)
        real                                :: crit_wtloss  ! Critical weight loss before transfer to humus (t/ha)
        real                                :: fuel_frac    ! Fraction of total fuel type
        real                                :: hum_Nnew     ! New value for humus N content (tN/ha)
        real                                :: hum_OMnew    ! New value for humus content (t/ha)
        real                                :: humCN        ! Humus C:N ratio
        real                                :: moss_litter  ! Moss litter (kg)
        real                                :: con_fuel     ! Conifer needle litter (kg)
        real                                :: dec_fuel     ! Deciduous leaf litter (kg)
        real                                :: twig_fuel    ! Twig litter (kg)
        real                                :: duff         ! Duff content (kg)
        real                                :: min_Nmin     ! N mineralization from the A-layer (tN/ha)
        real                                :: ddmult       ! Decay factor due to soil degree-days
        real                                :: mdepth       ! Depth of moss (live and dead) layer (cm)
        real                                :: leaf_litterN ! Leaf litter N content (tN/ha)
        integer                             :: lt           ! Litter class
        integer                             :: nc           ! Total number of decaying cohorts
        integer                             :: i, j, ix     ! Looping indices

        ! Initialize accumulators
        litCO2 = 0.0
        humCO2 = 0.0
        totCO2 = 0.0
        WDW_new = 0.0
        tot_Nimob = 0.0
        avail_n = 0.0
        lit_Nmin = 0.0
        hum_Nmin = 0.0
        tot_Nmin = 0.0
        soil%forest_litter = 0.0
        dec_lit = 0.0
        leaf_litterN = 0.0

        ! Reduce OOP table lookups
        do i = 1, NCOH_MAX
            do j = 1, NCOH_CHR
                C(i, j) = soil%cohorts(i, j)
            end do
        end do

        ! Get total number of cohorts
        nc = soil%ncohort

        ! Calculate leaf litter N and fresh deciduous leaf litter
        do i = 1, 12
            leaf_litterN = leaf_litterN + soil%litter(i)*litter_params(i, 2)
            if (i .ne. 6 .or. i .lt. 10) then
                dec_lit = dec_lit + soil%litter(i)
            end if
        end do
        soil%dec_fuel = dec_lit ! Save for moss subroutine (t/ha)

        ! Calculate effect of soil degree-days (>0degC) on decay rate
        ddmult = 1.8**(0.005*(soildays - 1900.0))
        if (ddmult .le. epsilon(1.0)) ddmult = 0.0

        ! Calculae moss (live and dead) depth (cm)
        mdepth = soil%M_depth*M_TO_CM +                                        &
            soil%litter(IMOSS)/HEC_TO_M2*plotsize*T_TO_KG/BULK_MOSS*M_TO_CM

        ! Calculate effect of moss and moisture on decay rate
        smult = (max((1.0 - sqrt(mdepth*3.0)), 0.0) +                          &
            max(1.0*(1.0-1.25*flooddays)**2.0, 0.1))/2.0

        ! Populate cohort array with new cohorts and cohort properties
        do i = 1, LIT_LEVS
            if (i .ne. IWDW) then
                ! Add in parameters for everything but well-decayed wood
                if (soil%litter(i) .gt. epsilon(1.0)) then
                    nc = nc + 1
                    ! Current weight (t/ha)
                    C(nc, 1) = soil%litter(i)*litter_params(i, 10)
                    ! Current N (tN/ha)
                    C(nc, 2) = soil%litter(i)*litter_params(i, 2)
                    ! Litter parameters (from input file)
                    do j = 3, 9
                        C(nc, j) = litter_params(i, j)
                    end do
                    ! Initial weight (t/ha)
                    C(nc, 10) = soil%litter(i)*litter_params(i, 10)
                    ! Current N %
                    C(nc, 11) = litter_params(i, 2)
                    ! Critical % remaining
                    C(nc, 12) = litter_params(i, 7)*1.7039 + 0.0955
                    ! Critical % remaining of freshwood and branches is 30%
                    if (C(nc, 5) >= 14.0 .and. C(nc, 5) <= 18.0) C(nc, 12) = 0.30
                    ! Age (years)
                    C(nc, 13) = 0.0
                end if
            end if
        end do

        ! Calculate decay multiplier based on LAI
        if (cla/plotsize .ge. 3.88) then
            cla_decaym = 1.0
        else
            cla_decaym = 1.0 + 1.75*sqrt(1.0 - (cla/plotsize)/3.88)
        end if

        ! Bypass litter cohort calculations if there is nothing to
          ! decay
        if (nc .gt. 1) then
            ! Loop to calculate litter decay, N immobilization, lignin
            ! decay, and litter CO2 evolution
            do i = 2, nc

                ! Calculate percent weight loss based on AET and lignin:N ratio
                ! Pastor & Post (1984)
                ! Note: can't use this equation for moss litter decay
                perc_wtloss = (0.9404 + 0.09352*aet_mm) -                      &
                    ((-0.4956 + 0.00193*aet_mm)*(C(i, 7)/C(i, 11)))
                perc_wtloss = (perc_wtloss*cla_decaym)/100.0
                if (perc_wtloss .gt. 1.0) perc_wtloss = 1.0

                ! Bonan's equation uses N:C ratio and active layer thickness
                ! Note: CAN use this for moss litter decay
                perc_wtloss2 = (-0.0052 + 2.08*(C(i,2)/C(i, 1)))*              &
                    exp(0.898*soil%active)
                perc_wtloss2 = max(perc_wtloss2*cla_decaym, 0.0)
                if (perc_wtloss2 .gt. 1.0) perc_wtloss2 = 1.0

                if (soil%active .le. 1.5) then
                    ! Average these two percent weight losses when
                    ! permafrost present
                    perc_wtloss = ((perc_wtloss2 + perc_wtloss)/2.0)*ddmult*smult
                    perc_wtloss = min(max(perc_wtloss, 0.00), 1.0)
                else
                    perc_wtloss = min(max(perc_wtloss*ddmult*smult, 0.0), 1.0)
                endif

                ! Litter cohort type
                lt = int(C(i, 5))

                ! Weight loss of moss litter can only use Bonan's eq.
                if (lt .eq. IMOSS) perc_wtloss = min(max(perc_wtloss2*ddmult*  &
                    smult, 0.0), 1.0)

                ! Max weight loss of large wood (DBH > 10cm) is 3%
                if (lt .eq. 15) perc_wtloss = min(perc_wtloss, 0.03)

                ! Max weight loss of small wood is 10%
                if (lt .eq. 14) perc_wtloss = min(perc_wtloss, 0.1)

                ! Maximum weight loss of well-decayed wood is 5%
                if (lt .eq. 19) perc_wtloss = min(perc_wtloss, 0.05)

                ! Weight loss of branches is less than 20%
                if (lt >= 16 .and. lt <= 18) then
                    perc_wtloss = min(perc_wtloss, 0.2)
                end if

                ! Calculate actual weight loss (t/ha)
                wtloss = perc_wtloss*C(i, 1)

                ! Calculate weight loss for litter other than moss
                if (lt .ne. IMOSS) then

                    ! Calculate fraction of organic matter remaining following
                    ! this weight loss
                    perc_rem = (C(i, 1) - wtloss)/(C(i, 10))

                    ! Calculate new N concentration
                    C(i, 11) = min(max(C(i, 3) - C(i, 4)*perc_rem,             &
                        epsilon(1.)), 1.0)

                    ! Retain cohort another year if fraction remaining is
                    ! greater than fraction which will become humus/WDW
                    if (perc_rem .le. C(i, 12)) then

                        ! Transfer cohort

                        ! Calculate actual weight loss and N concentration
                        wtloss = max(C(i, 1) - C(i, 12)*C(i, 10), 0.0)

                        C(i, 11) = min(max(C(i, 3)- C(i, 4)*C(i, 12),          &
                            epsilon(1.0)), 1.0)

                        ! Calculate absolute change in N content (tN/ha)
                        delta_N = C(i, 2) - C(i, 11)*(C(i, 1) - wtloss)

                        if (delta_N .lt. epsilon(1.0)) then
                            tot_Nimob = tot_Nimob  - delta_N
                            ! Negative - immobilize N
                        else
                            lit_Nmin = lit_Nmin + delta_N
                            ! Positive - mineralize N
                        end if

                        ! Tranfer cohorts to humus or well-decayed wood
                        if (C(i, 6) .eq. 1.0) then
                            ! Humus
                            C(1, 1) = C(1, 1) + C(i, 1) - wtloss
                            C(1, 2) = C(1, 2) + C(i, 11)*(C(i, 1) - wtloss)
                            C(i, 1) = 0.0
                        else
                            ! Well-decayed wood
                            WDW_new = WDW_new + C(i, 1) - wtloss
                            C(i, 1) = 0.0
                        end if
                    end if

                else if (lt .eq. IMOSS ) then

                    ! Decay moss litter
                    crit_wtloss = (C(i, 4)*C(i, 1) - C(i, 2))/                 &
                        (C(i, 3) + C(i, 4))

                    if (wtloss .ge. crit_wtloss) then
                        ! Transfer cohort to humus, weight loss proceeds
                        ! at critical weight loss rate
                        wtloss = crit_wtloss

                        ! Transfer cohort
                        C(i, 11) = C(i, 4)
                        C(1, 1) = C(1, 1) + C(i, 1) - wtloss
                        C(1, 2) = C(1, 2) + wtloss*C(i,3)
                        C(i, 1) = 0.0

                        ! Calculate absolute change in N content
                        delta_N = C(i, 2) - C(i, 11)*(C(i, 1) - wtloss)

                        if (delta_N .lt. epsilon(1.0)) then
                            ! Negative - immobilize N
                            tot_Nimob = tot_Nimob - delta_N
                        else
                            ! Posivite - mineralize N
                            lit_Nmin = lit_Nmin + delta_N
                        end if
                    end if
                end if

                ! Update weight and N content of cohorts that didn't get
                ! transferred to humus/WDW
                if (C(i, 1) .gt. epsilon(1.0)) then
                    ! Remove decay loss
                    C(i, 1) = C(i, 1) - wtloss

                    if (int(C(i,5)) .ne. IMOSS .and. int(C(i,5)) .ne. IWDW) then
                        ! Not moss or WDW
                        C(i, 2) = C(i, 1)*C(i, 11)
                        C(i, 7) = C(i, 8) - C(i, 9)*(C(i, 1)/C(i, 10))
                    else if (int(C(i, 5)) .eq. IMOSS) then
                        ! Moss
                        C(i, 2) = C(i, 2) + wtloss*C(i,3)
                        C(i, 11) = C(i,2)/C(i,1)
                    else
                        ! WDW wood
                        C(i, 2) = C(i, 1)*C(i, 11)
                        C(i, 7) = 0.0
                    end if
                    ! Update year
                    C(i, 13) = C(i, 13) + 1.0
                end if

                ! Update litter cohort CO2 release
                litCO2 = litCO2 + (wtloss*B_TO_C)

            end do ! End litter cohort loop

        end if !end if ncohort > 1

        ! Calculate humus percent weight loss
        perc_wtloss = (-0.0052 + 2.08*(C(1, 2)/C(1, 1)))*exp(0.898*soil%active)
        perc_wtloss = perc_wtloss*ddmult*smult
        perc_wtloss = min(perc_wtloss, 1.0)
        perc_wtloss = max(perc_wtloss, 0.0)

        ! Calculate humus N mineralization
        if (soil%active .le. 1.5) then
            ! Use Bonan (1990) equation if permafrost present
            hum_Nmin = C(1, 2)*perc_wtloss
        else
            ! Otherwise use Pastor & Post (1984)
            perc_wtloss = max(0.0, min(0.025*ddmult*smult, 1.0))
            hum_Nmin = C(1, 2)*perc_wtloss
        end if

        ! Subtract mineralized N from humus N and calculate humus CO2 release
        hum_Nnew = C(1, 2) - hum_Nmin
        hum_OMnew = C(1,1)*(hum_Nnew/C(1,2))
        humCO2 = (C(1,1) - hum_OMnew)*B_TO_C

        ! Update cohorts array
        C(1, 1) = hum_OMnew
        C(1, 2) = hum_Nnew

        ! Calculate humus C:N ratio
        humCN = (B_TO_C*C(1, 1))/C(1, 2)

        ! Add humus N mineralization to cohort N mineralization to get
        !total N mineralization
        tot_Nmin = lit_Nmin + hum_Nmin

        ! Throughfall is 16% of leaf litter N
        tfallNmin = 0.16*leaf_litterN

        ! Mineralization from A-layer
        min_Nmin = 10.0*min(soil%active, soil%A_depth)*0.001

        ! Subtract immobilization from total mineralization to get
        ! plant-available N
        ! Add throughfall N mineralization, mineral N, and N from fires
        avail_n = max(0.0, tot_Nmin - tot_Nimob + tfallNmin + soil%fan +       &
            min_Nmin)

        ! Calculate total soil respiration
        totCO2 = litCO2 + humCO2

        !remove transferred cohorts
        ix = 0
        do 20 i = 2, nc
            if (C(i, 1) .eq. 0.0) GO TO 16
                do 12 j = 1, NCOH_CHR
12                C(i-ix,j) = C(i, j)
            GO TO 20
16          ix = ix + 1
20       CONTINUE
        nc = nc - ix

        ! Create new well-decayed wood cohort
        if (WDW_new .gt. epsilon(1.0)) then
            nc = nc + 1
            ! Initital weight (t/ha)
            C(nc, 1) = WDW_new
            ! Current N content (tN/ha)
            C(nc, 2) = WDW_new*litter_params(IWDW, 2)
            ! Litter parameters (from input file)
            do j = 3, 9
                C(nc, j) = litter_params(IWDW, j)
            end do
            ! Current weight (t/ha)
            C(nc, 10) = WDW_new
            ! Current N %
            C(nc, 11) = litter_params(IWDW, 2)
            ! Critical percent remaining
            C(nc, 12) = 0.50
            ! Age (years)
            C(nc, 13) = 0.0
        end if

        ! Calculate total weight and N content by forest floor compartment and
        ! fuel conditions
         do i = 2, nc

            ! litter type
            lt = int(C(i, 5))

            ! Weight (t/ha)
            soil%forest_litter(lt, 1) = soil%forest_litter(lt, 1) + C(i, 1)

            ! N content (tN/ha)
            soil%forest_litter(lt, 2) = soil%forest_litter(lt, 2) + C(i, 2)

         end do

        ! Calculate organic layer depth

        ! Sum up leaf litter (t/ha)
        con_fuel = 0.0
        dec_fuel = 0.0
        do i = 1, 12
            if (i .eq. 6 .or. i .ge. 10) then
                con_fuel = con_fuel + soil%forest_litter(i, 1)
            else
                dec_fuel = dec_fuel + soil%forest_litter(i, 1)
            end if
        end do

        ! Convert to kg
        con_fuel = con_fuel/HEC_TO_M2*plotsize*T_TO_KG
        dec_fuel = dec_fuel/HEC_TO_M2*plotsize*T_TO_KG
        twig_fuel = (soil%forest_litter(16, 1)+                                &
            soil%forest_litter(17, 1) +                                        &
            soil%forest_litter(18, 1))/HEC_TO_M2*plotsize*T_TO_KG
        moss_litter = soil%forest_litter(20, 1)/HEC_TO_M2*plotsize*T_TO_KG
        duff = C(1, 1)/HEC_TO_M2*plotsize*T_TO_KG

        ! Calculate organic layer depth
        soil%O_depth = (1.0/plotsize)*((con_fuel/BULK_CON) +                   &
            (dec_fuel/BULK_DEC) + (twig_fuel/BULK_LITTER) + (duff/BULK_DUFF) + &
            (moss_litter/BULK_MOSS))

        ! Reassign incidence variables
        do i = 1, NCOH_MAX
            do j = 1, NCOH_CHR
                soil%cohorts(i, j) = C(i, j)
            end do
        end do
        soil%ncohort = nc

        ! Set fresh litter to 0.0
        do i = 1, LIT_LEVS
            soil%litter(i) = 0.0
        end do

    end subroutine soiln

    !:.........................................................................:

    subroutine forest_fuels(soil, drI, avail_fuel, consN, consRoot)
        !
        !  Calculates available fuels for forest fires and burns up litter
        !    cohorts accordingly
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    10/10/18    A. C. Foster         Original Code
        !

        ! Data dictionary: calling arguments
        class(SoilData), intent(inout) :: soil       ! Soil object
        real,            intent(in)    :: drI        ! Drought index (0-1)
        real,            intent(out)   :: avail_fuel ! Available fuel for forest fires (t/ha)
        real,            intent(out)   :: consN      ! Proportion of litter N consumed by fire
        real,            intent(out)   :: consRoot   ! Proportion of live roots consumed by fire

        ! Data dictionary: local variables
        real, dimension(NCOH_MAX, NCOH_CHR) :: C         ! Decaying cohorts
        real                                :: fol_fuel  ! Leaf litter available for burning (t/ha)
        real                                :: lt        ! Litter class
        real                                :: rfs       ! Moisture index (0-1)
        real                                :: vwc       ! Mineral layer volumetric moisture content
        real                                :: WDW_kill  ! Proportion well-decayed wood burned (0-1)
        real                                :: hum_kill  ! Proportion well-decayed wood burned (0-1)
        real                                :: M_loss    ! Depth of moss burned (m)
        real                                :: BO        ! Temporary value of moss biomass (kg)
        real                                :: sfan      ! Total N content volatilized by burning (tN/ha)
        real                                :: oldN      ! Temporary value for litter N content (tN/ha)
        real                                :: volN      ! N volitilized by burning (tN/ha)
        real                                :: humloss   ! Amount of humus burned (t/ha)
        real                                :: oldHN     ! Old value of humus N content
        real                                :: twigloss  ! Twig litter burned (t/ha)
        real                                :: stemloss  ! Bole litter burned (t/ha)
        real                                :: WDWloss   ! Well-decayed wood litter burned (t/ha)
        real                                :: rootloss  ! Root litter burned (t/ha)
        real                                :: con_fuel  ! Coniferous needle litter (kg)
        real                                :: dec_fuel  ! Deciduous leaf litter (kg)
        real                                :: moss_fuel ! Moss litter (kg)
        real                                :: twig_fuel ! Twig litter (kg)
        real                                :: duff      ! Duff content (kg)
        integer                             :: i, j, ix  ! Looping indices
        integer                             :: nc        ! Number of decaying cohorts
        integer                             :: ilt       ! Litter class

        ! Initialize accumulator
        sfan = 0.0

        ! Reduce table lookups
        nc = soil%ncohort
        do i = 1, NCOH_MAX
            do j = 1, NCOH_CHR
                C(i, j) = soil%cohorts(i, j)
            end do
        end do

        ! Calculate available fuels

        ! Foliage - all leaf litter is available for burning
        fol_fuel = 0.0
        do i = 1, 12
            fol_fuel = fol_fuel + soil%forest_litter(i, 1)
        end do
        soil%fol_fuel = fol_fuel

        ! Branches, small boles, and large boles
          ! Equations for fuel amount are from Schumacher et al. (2006),
        ! Landscape Ecology
        soil%twig_fuel = (TFC1 + TFC2*drI)*(soil%forest_litter(ITW, 1) +       &
            soil%forest_litter(ISBR, 1) + soil%forest_litter(ILBR, 1))
        soil%smbl_fuel = (TFC1 + TFC2*drI)*(soil%forest_litter(ISBL, 1))
        soil%lrbl_fuel = (BFC*drI)*soil%forest_litter(ILBL, 1)

        avail_fuel = fol_fuel + soil%twig_fuel + soil%smbl_fuel +         &
            soil%lrbl_fuel
        soil%avail_fuel = avail_fuel

        ! Calculate fire severity for humus, root, and well-decayed wood
        ! consumption
        vwc = soil%minWC*soil%A_bulk_dens/H2O_D
        rfs = max(0.0, (soil%sat(2) - vwc)/(soil%sat(2) - soil%pwp(2)))

        consRoot = max(0.0, min((0.302 + 0.597*rfs +                          &
            3.34*(soil%M_depth + soil%O_depth)), 0.9))

        WDW_kill = max(0.0, min((0.098 + 0.597*rfs +                           &
            3.34*(soil%M_depth + soil%O_depth)), 0.5))

        hum_kill = max(0.0, min((0.079 + 0.5744*rfs +                          &
            3.34*(soil%M_depth + soil%O_depth)), 0.658))

        consN = max(0.0, min((0.6426*rfs +                                     &
            3.34*(soil%M_depth + soil%O_depth)), 0.7))

        ! Set % nitrogen in forest floor not volatalized
        volN = 1.0 - consN

        ! Save old value of moss_biom
        BO = soil%moss_biom

        ! Remove moss layer
        M_loss = min(soil%M_depth, hum_kill*soil%M_depth)
        soil%M_depth = max(soil%M_depth - M_loss, 0.0)
        soil%moss_biom = soil%M_depth*28.6*plotsize
        sfan = sfan + (BO - soil%moss_biom)*litter_params(IMOSS, 2)*volN

        ! Remove humus layer
        humloss = min(hum_kill*C(1, 1), C(1, 1))
        oldHN = C(1, 2)
        C(1, 2) = (1.0 - humloss/C(1,1))*C(1, 2)
        sfan = sfan + (oldHN - C(1, 2))*volN
        C(1, 1) = max(C(1, 1) - hum_kill*C(1, 1), 0.0)

        ! Get rid of cohorts
        do i = 2, nc

            lt = C(i, 5)

            ! All leaf litter fuels decomposed
            if (lt .le. 12.0) then
                C(i, 1) = 0.0
                sfan = sfan + C(i, 2)*volN
            end if

            !roots
            if (lt == float(IROOT)) then
                rootloss = min(consRoot*C(i, 1), C(i, 1))
                C(i, 1) = max(C(i, 1) - consRoot*C(i, 1), 0.0)
                oldN = C(i, 2)
                C(i, 2) = C(i, 2)*(1.0 - consRoot)
                C(i, 11) = C(i, 2)/C(i,1)
                sfan = sfan + (oldN - C(i, 2))*volN
            end if

            ! Twigs and small trees
            if (lt >= 14.0 .or. lt <= 16.0) then
                if (lt /=15.0) then
                    twigloss = min((TFC1 + TFC2*drI)*C(i, 1), C(i, 1))
                    C(i, 1) = max(C(i, 1) - (TFC1 + TFC2*drI)*C(i, 1), 0.0)
                    oldN = C(i, 2)
                    C(i, 2) = C(i, 2)*(1.0 - (TFC1 + TFC2*drI))
                    C(i, 11) = C(i, 2)/C(i,1)
                    sfan = sfan + (oldN - C(i, 2))*volN
                end if
            end if

            ! Boles
            if (lt .eq. 15.0) then
                stemloss = min((BFC*drI)*C(i, 1), C(i, 1))
                C(i, 1) = max(C(i, 1) - (BFC*drI)*C(i, 1), 0.0)
                oldN = C(i, 2)
                C(i, 2) = C(i, 2)*(1.0 - (BFC*drI))
                C(i, 11) = C(i, 2)/C(i,1)
                sfan = sfan + (oldN - C(i, 2))*volN
            end if

            ! Well-decayed wood
            if (lt == float(IWDW)) then
                WDWloss = min(WDW_kill*C(i, 1), C(i, 1))
                C(i, 1) = max(C(i, 1) - WDW_kill*C(i, 1), 0.0)
                oldN = C(i, 2)
                C(i, 2) = C(i, 2)*(1.0 - WDW_kill)
                C(i, 11) = C(i, 2)/C(i,1)
                sfan = sfan + (oldN - C(i, 2))*volN
            end if

            ! All moss litter consumed
            if (lt == float(IMOSS)) then
                C(i, 1) = 0.0
            end if

        end do

        ! Remove completely burned cohorts
        ix = 0
        do i = 1, nc
            if (C(i, 1) .ge. epsilon(1.0)) then
                do j = 1, NCOH_CHR
                    C(i-ix, j) = C(i,j)
                end do
            end if
            ix = ix + 1
        end do

        ! Calculate total weight and N content by forest floor compartment and
        ! fuel conditions
         do i = 2, nc

            ! litter type
            ilt = int(C(i, 5))

            ! Weight (t/ha)
            soil%forest_litter(ilt, 1) = soil%forest_litter(ilt, 1) + C(i, 1)

            ! N content (tN/ha)
            soil%forest_litter(ilt, 2) = soil%forest_litter(ilt, 2) + C(i, 2)

         end do

        ! Calculate organic layer depth

        ! Sum up leaf litter (t/ha)
        con_fuel = 0.0
        dec_fuel = 0.0
        do i = 1, 12
            if (i == 6 .or. i >= 10) then
                con_fuel = con_fuel + soil%forest_litter(i, 1)
            else
                dec_fuel = dec_fuel + soil%forest_litter(i, 1)
            end if
        end do

        ! Convert to kg
        con_fuel = con_fuel/HEC_TO_M2*plotsize*T_TO_KG
        dec_fuel = dec_fuel/HEC_TO_M2*plotsize*T_TO_KG
        twig_fuel = (soil%forest_litter(16, 1)+                                &
            soil%forest_litter(17, 1) +                                        &
            soil%forest_litter(18, 1))/HEC_TO_M2*plotsize*T_TO_KG
        moss_fuel = soil%forest_litter(20, 1)/HEC_TO_M2*plotsize*T_TO_KG
        duff = C(1, 1)/HEC_TO_M2*plotsize*T_TO_KG

        ! Calculate organic layer depth
        soil%O_depth = (1.0/plotsize)*((con_fuel/BULK_CON) +                   &
            (dec_fuel/BULK_DEC) + (twig_fuel/BULK_LITTER) + (duff/BULK_DUFF) + &
            (moss_fuel/BULK_MOSS))

        ! Reassign incidence variables
        do i = 1, NCOH_MAX
            do j = 1, NCOH_CHR
                soil%cohorts(i, j) = C(i, j)
            end do
        end do
        soil%ncohort = nc


    end subroutine forest_fuels

    !:.........................................................................:

end module Soil
