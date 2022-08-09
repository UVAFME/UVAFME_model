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
    real,    parameter :: LAI_MIN = 0.01    ! Scalar for calculating min canopy moisture
    real,    parameter :: LAI_MAX = 0.15    ! Scalar for calculating max canopy moisture
    real,    parameter :: SAV_DUFF = 20.0   ! Duff SAV (/cm)
    real,    parameter :: MOSS_SAV = 70.0   ! Moss SAV (/cm)
    real,    parameter :: SHRUB_SAV = 80.0  ! Live shrub SAV (/cm)
    real,    parameter :: BULK_SHRUB = 20.0 ! Live shrub bulk density (kg/m3)
    integer, parameter :: NCOH_MAX = 1500   ! Maximum number of decaying cohorts
    integer, parameter :: NCOH_CHR = 15     ! Number of cohort characteristics
    integer, parameter :: IMOSS = 20        ! Index location of moss in soil%litter array
    integer, parameter :: IWDW = 19         ! Index location of well-decayed wood in soil%litter array
    real,    parameter :: H = 18000.0       ! Heat content of fuel (kJ/kg)

    ! Define soil type
    type SoilData
        real, dimension(NCOH_MAX, NCOH_CHR) :: cohorts       ! Array of decaying cohorts of litter
        real, dimension(LIT_LEVS)           :: litter        ! Fresh litter content (t/ha)
        real, dimension(LIT_LEVS+1, 2)      :: forest_litter ! Total forest litter (t/ha)
        real, dimension(FL_LEVS)            :: fuels         ! Fuel loading (kg/m2)
        real, dimension(FL_LEVS)            :: litBD         ! Fuel bulk density (kgC/m3)
        real, dimension(FL_LEVS)            :: litSAV        ! Fuel SAV (/cm)
        real, dimension(FL_LEVS)            :: litter_moist  ! Fuel moisture (volumetric)
        real, dimension(FL_LEVS)            :: litter_mef    ! Fuel moisture of extinction (vol)
        real, dimension(FL_LEVS)            :: frac_burnt    ! Fraction burnt
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
        real                                :: shrubLitter   ! Shrub live foliage and fine twigs (kg/m2)
        real                                :: fuel_sum      ! Total fuel loading (kg/m2)
        real                                :: fuel_SAV      ! Average fuel SAV (/cm)
        real                                :: fuel_BD       ! Average fuel bulk density (kg/m3)
        real                                :: fuel_moisture ! Average fuel moisture (volumetric)
        real                                :: MEF           ! Average moisture of extinction (volumetric)
        real                                :: I_surf        ! Surface fire intensity (kW/m)
        real                                :: rosf          ! Surface rate of spread (m/min)
        real                                :: rosf_phi      ! Surface rate of spread without phi (m/min)
        real                                :: fuel_consumed ! Fuel consumed by fire (kg/m2)
        real                                :: tau_l         ! Residence time of fire (min)
        real                                :: dmc_fire      ! Duff moisture code for day of fire
        real                                :: ffmc_fire     ! Fine fuel moisture code for day of fire
        real                                :: wind_fire     ! Wind speed for day of fire
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
        real,                      intent(out)    :: zdepth ! Total depth of freeze/thaw (m)

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
        if (soil%M_depth > 0) then
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
        if (knd == 1) then
            cadd = fdd(month, 2)*cff*(2.0 - cfs)
        else
            cadd = tdd(month, 2)*cft*cfs
        end if

        ! Calculate freeze/thaw of soil layers
        do l = 1, 2
            if (l == 1) then
                dl = Od
            else
                dl = Ad  + 2.0
            end if

            if (l == 1) then
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
                if (soil%itxt <= 1) then ! Granular
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
            if (knd == 1) then
                tk = tkf
            else
                tk = tku
            end if

            ! Latent heat of fusion (kcal/m3)
            if (l == 1) then
                ql = 80.0*(soil%wc(l))*OBD
            else
                ql = 80.0*(soil%wc(l))*ABD
            end if

            ! Calculate thermal resistance
            res = dl/tk

            ! Calculate degree-days required to freeze entire layer
            degd = ql*dl/24.0*(tyr + res/2.0)

            if (degd <= cadd) then
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
        slope, amlt, xmlt, xfrz, tdd, k, aet, flow, aow0_ByMin, saw0_ByFC,     &
        saw0_ByWP, saw0_BySAT, soilW, soilI)
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
        real, parameter, dimension(4) :: EXS_P = [4.0, 1.0, 0.4, 0.1]

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
        real,                       intent(in)    :: flow       ! Moisture input from overland flow (mm)
        real,                       intent(out)   :: aet        ! Actual evapotranspiration (cm)
        real,                       intent(out)   :: aow0_ByMin ! Organic layer moisture scaled by wilting point
        real,                       intent(out)   :: saw0_ByFC  ! Mineral layer moisture scaled by field capacity
        real,                       intent(out)   :: saw0_ByWP  ! Mineral layer moisture scaled by wilting point
        real,                       intent(out)   :: saw0_BySAT ! Mineral layer moisture scaled by saturation capacity
        real,                intent(in), optional :: soilW      ! Moisture from monthly climate input (volumetric)
        real,                intent(in), optional :: soilI      ! Ice content from monthly climate input (volumetric)

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
        if (ta >= TSMAX) then

            ! All liquid rainfall
            pr = precip
            ps = 0.0

            ! Calculate canopy interception and throughfall
            canopy_int = min(max((laiw_max - laiw), 0.0), pr)
            laiw = laiw + canopy_int
            tfall = max(pr - canopy_int, 0.0)

        else if (ta < TSMAX .and. ta > TSMIN) then

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

        else if (ta <= TSMIN) then

            ! Only snow - no canopy interception or throughfall
            pr = 0.0
            ps = precip
            canopy_int = 0.0
            tfall = 0.0

            ! Accumulate snowpack
            soil%swe = soil%swe + ps

        end if

        ! Melt snowpack if it exists and if ta is above tc
        if (soil%swe >= epsilon(1.0)) then
            if (ta > TC) then

                ! Melt the snowpack
                melt = (MELTFACT*(ta - TC))*0.001

                if (soil%swe >= melt) then

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
        if (zh(1) > 0.0) then
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
            if ((H2Oice(l) - wt(l)) >= epsilon(1.0)) then

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
            if (l == 1) then
                if (pwl >= 0.0) then
                    ! More water available than pet, aet = pet
                    aet = pet

                    ! Water gain is water from throughfall and snowmelt
                    pwg = pwl
                else
                    ! No excess water from rainfall (aet < pet)
                    pwg = 0.0
                end if
            else
                if (pwl >= 0.0) then
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
            if (exs(l) > 0.0) then
                gwl(l) = ((EXS_P(soil%itxt + 1)*exs(l)**2)*                    &
                    (water(l)))/((pet + exs(l))*(fc(l)))
                gwl(l) = max(gwl(l), 0.0)
                gwl(l) = min(gwl(l), 1.0)
            end if
            exs(l) = exs(l)*gwl(l)
            water(l) = max(water(l) - exs(l), 0.0)

        end do

        ! Calculate runoff
        soil%runoff = exs(2) + loss_slp

        if (pwl < 0.0) then
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
            if (l == 1) then
                pwll = min(0.0, pwl*rz) ! Either negative or 0.0
            else
                pwll = min(0.0, pwl*(1.0-rz)) ! Either negative or 0.0
            end if

            if (dmelt(l) > 0.0) then
                ! We have some unfrozen soil to take water from

                ! Calculate aet loss from soil
                B = 0.461 - 1.10559/(zdrain(l)/zh(l))
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
            if (tdd(k,1) > 0.0) then
                if (dmelt(l) > 0.0) then

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
            if ((water(l) - wf(l)) >= epsilon(1.0)) then

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
            if (zh(l) > 0.0) then
                if (l == 1) then
                    wc(l) = ((water(l) + H2Oice(l))/zh(l)*H2O_D/OBD)
                else
                    wc(l) = ((water(l) + H2Oice(l))/zh(l)*H2O_D/ABD)
                end if
            else
                wc(l) = 0.0
            end if

            ! Adjust for above layers
            xabove = xabove + zh(l)
        end do

        ! Calculate scaled water content
        if (pwp(1) > 0.0) then
            aow0_ByMin = water(1)/pwp(1)
        else
            aow0_ByMin = 0.0
        end if

        if (fc(2) > 0.0) then
            saw0_ByFC = water(2)/fc(2)
        else
            saw0_ByFC = 1.0
        end if

        if (pwp(2) > 0.0) then
            saw0_ByWP = water(2)/pwp(2)
        else
            saw0_ByWP = 2.0
        end if

        if (sat(2) > 0.0) then
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
        real, parameter :: Q = 0.12       ! Mortality parameter
        real, parameter :: B = 0.136      ! Mortality parameter
        real, parameter :: EXT = 0.5      ! Light extinction coefficient
        real, parameter :: SLA = 1.0       ! Specific leaf area (m2/kg)
        real, parameter :: SPORES = 0.001 ! Proportion of biomass spent on reproduction
        real, parameter :: PMAX = 0.35    ! Maximum moss production (kg/m2/yr)
        real, parameter :: LRMIN = 0.01   ! Light compensation point
        real, parameter :: LRMAX = 0.05   ! Light compensation point

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
        real :: biokg     ! Moss biomass (kg/m2)
        real :: al        ! Available light
        real :: algf      ! Available light growth factor
        real :: fcgf      ! Forest cover growth factor
        real :: dlgf      ! Deciduous leaf litter growth factor
        real :: ddgf      ! Moisture growth factor
        real :: assim     ! Moss assimilation rate (kg/m2)
        real :: assim_eff ! Effective assimilation (kg/kg)
        real :: prod      ! Moss production (kg/m2)
        real :: litter    ! Moss litter (kg)

        ! Convert moss biomass in kg to kg/m2
        biokg = soil%moss_biom/plotsize

        ! Light growth multiplier
        al = exp(-1.0*EXT*(cla/plotsize + biokg*SLA))
        algf = (al - LRMIN)/(LRMAX - LRMIN)
        algf = max(0.0, algf)
        algf = min(1.0, algf)

        ! Forest cover growth multiplier (alff > 0.75)
        if (alff > 0.75) then
           fcgf = 1.5625 - alff**2
        else
            fcgf = 1.0
        end if
        if (fcgf > 1.0) fcgf = 1.0
        if (fcgf < 0.0) fcgf = 0.0

        ! Deciduous leaf litter growth multiplier
        if (decLit > 0.0) then
            dlgf = exp(-0.2932*decLit)
            if (dlgf <= 0.0) dlgf = 0.0
            if (dlgf >= 1.0) dlgf = 1.0
        else
            dlgf = 1.0
        end if

        ! Soil moisture growth multiplier
        if (drydays > 0.025) then
            ddgf = 0.0
        else
            ddgf = 1.0
        end if

        ! Moss assimilation
        assim = PMAX*algf*fcgf*dlgf*ddgf

        ! Effective assimilation
        assim_eff = SLA*assim*(1.0 - SPORES*dlgf*ddgf)
        prod = assim_eff*biokg - biokg*Q - biokg*B + SPORES*dlgf*ddgf

        if (biokg + prod < 0.0) then
            ! Not enough moss to account for mortality/respiration
            ! Set moss loss to all of current biomass
            prod = -1.0*biokg
        end if

        soil%moss_biom = (biokg + prod)*plotsize

        ! Calculate litter (kg)
        litter = (assim_eff*biokg + SPORES*dlgf*ddgf - prod)*plotsize
        if (litter < 0.0) then
            soil%litter(IMOSS) = 0.0
        else
            soil%litter(IMOSS) = litter
        end if

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
        !                                        small branches, large branches,
        !                                        and added calculation of BD and
        !                                        SAV of each cohort
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
        real, dimension(FL_LEVS)            :: fuels        ! Fuel loading (kg/m2)
        real, dimension(FL_LEVS)            :: litBD        ! Fuel bulk density (kg/m3)
        real, dimension(FL_LEVS)            :: litSAV       ! Fuel SAV (/cm)
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
        real                                :: fuel_sum     ! Total fuel (kg/m2)
        real                                :: sav          ! Cohort SAV (/cm)
        real                                :: bd           ! Cohort bulk density (kg/m3)
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
        real                                :: bd_slope     ! Temporary variable for calculating bulk density
        real                                :: sav_slope    ! Temporary variable for calculating SAV
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
        fuels = 0.0
        litBD = 0.0
        litSAV = 0.0
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
            if (i /= 6 .or. i < 10) then
                dec_lit = dec_lit + soil%litter(i)
            end if
        end do
        soil%dec_fuel = dec_lit ! Save for moss subroutine (t/ha)

        ! Calculate effect of soil degree-days (>0degC) on decay rate
        ddmult = 1.8**(0.0025*(soildays - 1900.0))
        if (ddmult <= epsilon(1.0)) ddmult = 0.0

        ! Calculae moss (live and dead) depth (m)
        mdepth = soil%M_depth +                                                &
            soil%litter(IMOSS)/HEC_TO_M2*plotsize*T_TO_KG/BULK_MOSS

        ! Calculate effect of moss and moisture on decay rate
        !smult = (1.0*(1.0 - 1.2*flooddays)**2)*(1.0 - sqrt(mdepth)*0.3)
        smult = ((1.0 - flooddays)**2)*(1.0 - sqrt(mdepth)*2.0)
       !  smult = min(1.0, smult)
        smult = max(0.05, smult) ! add a random comment

        ! Populate cohort array with new cohorts and cohort properties
        do i = 1, LIT_LEVS
            if (i /= IWDW) then
                ! Add in parameters for everything but well-decayed wood
                if (soil%litter(i) > epsilon(1.0)) then
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
                    ! Current BD (kg/m3)
                    C(nc, 14) = litter_params(i, 11)
                    ! Current SAV (/cm)
                    C(nc, 15) = litter_params(i, 12)
                end if
            end if
        end do

        ! Calculate decay multiplier based on LAI
        if (cla/plotsize >= 2.5) then
            cla_decaym = 1.0
        else
            cla_decaym = 1.0 + 1.5*sqrt(1.0 - (cla/plotsize)/2.5)
        end if

        ! Bypass litter cohort calculations if there is nothing to
          ! decay
        if (nc > 1) then
            ! Loop to calculate litter decay, N immobilization, lignin
            ! decay, and litter CO2 evolution
            do i = 2, nc

                ! Calculate percent weight loss based on AET and lignin:N ratio
                ! Pastor & Post (1984)
                ! Note: can't use this equation for moss litter decay
                perc_wtloss = (0.9804 + 0.09352*aet_mm) -                      &
                    ((-0.4956 + 0.00193*aet_mm)*(C(i, 7)/C(i, 11)))
                perc_wtloss = (perc_wtloss*cla_decaym)/100.0
                if (perc_wtloss > 1.0) perc_wtloss = 1.0

                ! Bonan's equation uses N:C ratio and active layer thickness
                ! Note: CAN use this for moss litter decay
                perc_wtloss2 = (-0.0052 + 2.08*(C(i, 2)/C(i, 1)))*             &
                    exp(0.898*soil%active)
                perc_wtloss2 = max(perc_wtloss2*cla_decaym, 0.0)
                if (perc_wtloss2 > 1.0) perc_wtloss2 = 1.0

                if (soil%active <= 2.0) then
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
                if (lt == IMOSS) perc_wtloss = min(max(perc_wtloss2*ddmult*  &
                    smult, 0.0), 1.0)

                ! Max weight loss of large wood (DBH > 10cm) is 3%
                if (lt == 15) perc_wtloss = 0.03

                ! Max weight loss of small wood is 10%
                if (lt == 14) perc_wtloss = 0.10

                ! Maximum weight loss of well-decayed wood is 5%
                if (lt == 19) perc_wtloss = 0.05

                ! Weight loss of branches is less than 20%
                if (lt >= 16 .and. lt <= 18) then
                    perc_wtloss = min(perc_wtloss, 0.2)
                end if

                ! Calculate actual weight loss (t/ha)
                wtloss = perc_wtloss*C(i, 1)

                ! Calculate weight loss for litter other than moss
                if (lt /= IMOSS) then

                    ! Calculate fraction of organic matter remaining following
                    ! this weight loss
                    perc_rem = max((C(i, 1) - wtloss)/C(i, 10), 0.0)

                    ! Calculate new N concentration
                    C(i, 11) = min(max(C(i, 3) - C(i, 4)*perc_rem,             &
                        epsilon(1.0)), 1.0)

                    ! Retain cohort another year if fraction remaining is
                    ! greater than fraction which will become humus/WDW
                    if (perc_rem <= C(i, 12)) then

                        ! Transfer cohort

                        ! Calculate actual weight loss and N concentration
                        wtloss = max(C(i, 1) - C(i, 12)*C(i, 10), 0.0)

                        C(i, 11) = min(max(C(i, 3)- C(i, 4)*C(i, 12),          &
                            epsilon(1.0)), 1.0)

                        ! Calculate absolute change in N content (tN/ha)
                        delta_N = C(i, 2) - C(i, 11)*(C(i, 1) - wtloss)

                        if (delta_N < epsilon(1.0)) then
                            tot_Nimob = tot_Nimob  - delta_N
                            ! Negative - immobilize N
                        else
                            lit_Nmin = lit_Nmin + delta_N
                            ! Positive - mineralize N
                        end if

                        ! Tranfer cohorts to humus or well-decayed wood
                        if (C(i, 6) == 1.0) then
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

                else if (lt == IMOSS) then

                    ! Decay moss litter
                    crit_wtloss = (C(i, 4)*C(i, 1) - C(i, 2))/                 &
                        (C(i, 3) + C(i, 4))

                    if (wtloss >= crit_wtloss) then
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

                        if (delta_N < epsilon(1.0)) then
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
                if (C(i, 1) > epsilon(1.0)) then

                    ! Remove decay loss
                    C(i, 1) = max(C(i, 1) - wtloss, 0.0)

                    ! Calculate critical percent left before transfer to humus
                    perc_rem = C(i, 1)/C(i, 10)

                    if (lt == IMOSS) then
                         crit_wtloss = (C(i, 4)*C(i, 1) - C(i, 2))/             &
                             (C(i, 3) + C(i, 4))
                         crit_perc = (C(i, 1) - crit_wtloss)/C(i,10)
                     else
                         crit_perc = C(i, 12)
                     end if

                    ! Update bulk density and SAV
                     sav_slope = (litter_params(lt, 12) - SAV_DUFF)/           &
                        (1.0 - crit_perc)
                     sav = (sav_slope*perc_rem) + (-1.0*sav_slope) +           &
                        litter_params(lt, 12)
                    C(i, 15) = sav

                    bd_slope = (litter_params(lt, 11) - BULK_LITTER)/          &
                        (1.0 - crit_perc)
                    bd = (bd_slope*perc_rem) + (-1.0*bd_slope) +               &
                        litter_params(lt, 11)
                    C(i, 14) = bd

                    if (C(i, 1)/C(i, 10) > 0.7) then
                    ! Divide into compartments (kg/m2)
                    if (lt <= 12) then
                        if (lt /= 6 .and. lt < 10) then
                            ! Deciduous leaf litter
                            fuels(FLDEC) = fuels(FLDEC) + C(i, 1)*T_TO_KG/     &
                                HEC_TO_M2
                        else
                            ! Conifer needle litter
                            fuels(FLCON) = fuels(FLCON) + C(i, 1)*T_TO_KG/     &
                                HEC_TO_M2
                        end if
                    else if (lt == 16) then
                        ! Twig litter
                        fuels(FLTW) = fuels(FLTW) + C(i, 1)*T_TO_KG/HEC_TO_M2
                    else if (lt == 17) then
                        ! Small branch litter
                        fuels(FLSBR) = fuels(FLSBR) + C(i, 1)*T_TO_KG/HEC_TO_M2
                    else if (lt == 18) then
                        ! Large branch litter
                        fuels(FLLBR) = fuels(FLLBR) + C(i, 1)*T_TO_KG/HEC_TO_M2
                    else if (lt == 14 .or. lt == 15) then
                        ! Bole litter
                        fuels(FLBL) = fuels(FLBL) + C(i, 1)*T_TO_KG/HEC_TO_M2
                    else if (lt == 13) then
                        fuels(FLRT) = fuels(FLRT) + C(i, 1)*T_TO_KG/HEC_TO_M2
                    else if (lt == IMOSS) then
                        fuels(FLDM) = fuels(FLDM) + C(i, 1)*T_TO_KG/HEC_TO_M2
                    end if
                    end if

                    if (int(C(i,5)) /= IMOSS .and. int(C(i,5)) /= IWDW) then
                        ! Not moss or WDW
                        C(i, 2) = C(i, 1)*C(i, 11)
                        C(i, 7) = C(i, 8) - C(i, 9)*(C(i, 1)/C(i, 10))
                    else if (int(C(i, 5)) == IMOSS) then
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

        ! Calculate humus N mineralization and percent weight loss
        if (soil%active <= 2.0) then
            ! Use Bonan (1990) equation if permafrost present
            perc_wtloss = (-0.0052 + 2.08*(C(1, 2)/C(1, 1)))*exp(0.898*soil%active)
            perc_wtloss = perc_wtloss*ddmult*smult*cla_decaym
            perc_wtloss = min(perc_wtloss, 1.0)
            perc_wtloss = max(perc_wtloss, 0.0)
            hum_Nmin = C(1, 2)*perc_wtloss
        else
            ! Otherwise use Pastor & Post (1984)
            perc_wtloss = max(0.0, min(0.035*ddmult*smult*cla_decaym, 1.0))
            hum_Nmin = C(1, 2)*perc_wtloss
        end if

        ! Subtract mineralized N from humus N and calculate humus CO2 release
        hum_Nnew = C(1, 2) - hum_Nmin
        hum_OMnew = max(C(1,1) - C(1,1)*perc_wtloss, 0.0)
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
        avail_n = max(0.0, tot_Nmin - tot_Nimob +  tfallNmin + soil%fan +      &
            min_Nmin)

        ! Calculate total soil respiration
        totCO2 = litCO2 + humCO2

        ! Remove transferred cohorts
        ix = 0
        do i = 2, nc
            if (C(i, 1) <= epsilon(1.0)) then
                ix = ix + 1
            else
                do j = 1, NCOH_CHR
                    C(i - ix, j) = C(i, j)
                end do
            end if
        end do
        nc = nc - ix

        ! Create new well-decayed wood cohort
        if (WDW_new > epsilon(1.0)) then
            nc = nc + 1
            if (nc > NCOH_MAX) then
                print *, 'no'
            end if
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
            ! Current bulk density (kg/m3)
            C(nc, 14) = litter_params(IWDW, 11)
            ! Current SAV (/cm)
            C(nc, 15) = litter_params(IWDW, 12)
        end if

        ! Live moss
        fuels(FLLM) = soil%moss_biom/plotsize

        ! Shrubs
        fuels(FLSH) = soil%shrubLitter

        ! Sum up fuels (no boles or roots)
        fuel_sum = sum(fuels) - fuels(FLBL) - fuels(FLRT)

        ! Calculate total weight and N content by forest floor compartment and
        ! fuel conditions
         do i = 2, nc

            ! litter type
            lt = int(C(i, 5))

            ! Weight (t/ha)
            soil%forest_litter(lt, 1) = soil%forest_litter(lt, 1) + C(i, 1)

            ! N content (tN/ha)
            soil%forest_litter(lt, 2) = soil%forest_litter(lt, 2) + C(i, 2)

            ! Divide into compartments
            if (C(i, 1)/C(i, 10) > 0.7) then
            if (lt <= 12) then
                if (lt /= 6 .and. lt < 10) then
                    ! Deciduous leaf litter
                    fuel_frac = (C(i, 1)*T_TO_KG/HEC_TO_M2)/fuels(FLDEC)
                    litBD(FLDEC) = litBD(FLDEC) + C(i, 14)*fuel_frac
                    litSAV(FLDEC) = litSAV(FLDEC) + C(i, 15)*fuel_frac
                else
                    ! Coniferous leaf litter
                    fuel_frac = (C(i, 1)*T_TO_KG/HEC_TO_M2)/fuels(FLCON)
                    litBD(FLCON) = litBD(FLCON) + C(i, 14)*fuel_frac
                    litSAV(FLCON) = litSAV(FLCON) + C(i, 15)*fuel_frac
                end if
            else if (lt == 16) then
                ! Twig litter
                fuel_frac = (C(i, 1)*T_TO_KG/HEC_TO_M2)/fuels(FLTW)
                litBD(FLTW) = litBD(FLTW) + C(i, 14)*fuel_frac
                litSAV(FLTW) = litSAV(FLTW) + C(i, 15)*fuel_frac

            else if (lt == 17) then
                ! Small branch litter
                fuel_frac = (C(i, 1)*T_TO_KG/HEC_TO_M2)/fuels(FLSBR)
                litBD(FLSBR) = litBD(FLSBR) + C(i, 14)*fuel_frac
                litSAV(FLSBR) = litSAV(FLSBR) + C(i, 15)*fuel_frac
            else if (lt == 18) then
                ! Large branch litter
                fuel_frac = (C(i, 1)*T_TO_KG/HEC_TO_M2)/fuels(FLLBR)
                litBD(FLLBR) = litBD(FLLBR) + C(i, 14)*fuel_frac
                litSAV(FLLBR) = litSAV(FLLBR) + C(i, 15)*fuel_frac
            else if (lt == 14 .or. lt == 15) then
                ! Bole litter
                fuel_frac = (C(i, 1)*T_TO_KG/HEC_TO_M2)/fuels(FLBL)
                litBD(FLBL) = litBD(FLBL) + C(i, 14)*fuel_frac
                litSAV(FLBL) = litSAV(FLBL) + C(i, 15)*fuel_frac

            else if (lt == IMOSS) then
                ! Moss litter
                fuel_frac = (C(i, 1)*T_TO_KG/HEC_TO_M2)/fuels(FLDM)
                litBD(FLDM) = litBD(FLDM) + C(i, 14)*fuel_frac
                litSAV(FLDM) = litSAV(FLDM) + C(i, 15)*fuel_frac
            end if
            end if
         end do

         ! Live moss
         litBD(FLLM) = BULK_MOSS
         litSAV(FLLM) = MOSS_SAV

         ! Dead moss
         !litBD(FLDM) = BULK_MOSS
         !litSAV(FLDM) = MOSS_SAV

         ! Live shrubs
         litBD(FLSH) = BULK_SHRUB
         litSAV(FLSH) = SHRUB_SAV

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
        twig_fuel = soil%forest_litter(ITW, 1)/HEC_TO_M2*plotsize*T_TO_KG
        moss_litter = soil%forest_litter(IMOSS, 1)/HEC_TO_M2*plotsize*T_TO_KG
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
        soil%fuel_sum = fuel_sum
        soil%fuels = fuels
        soil%litBD = litBD
        soil%litSAV = litSAV

        ! Set fresh litter to 0.0
        do i = 1, LIT_LEVS
            soil%litter(i) = 0.0
        end do


    end subroutine soiln

    !:.........................................................................:

    subroutine active_passive_check(soil, canopy_bh, canopy_bd, canopy_biom,   &
        wind_ms, ffmc, R_a, rosf_active, CFB, R_final, I_final,                &
        passive_crowning, active_crowning)
        !
        !  Checks whether a fire is and active or passive crown fire
        !  From Van Wagner 1977
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    02/03/21    A. C. Foster         Original Code
        !    06/25/21    A. C. Foster         Updated to use Scott & Reinhardt
        !                                     2001 equations
        !

        ! Data dictionary: constants
        real, parameter :: S_0 = 0.05  ! Critical mass flow rate (kg/m2/s)
        real, parameter :: HFL = 12700 ! Heat content of live fuel (kJ/kg)
        real, parameter :: FMC = 100.0 ! Default foliar moisture content (%)

        ! Data dictionary: calling arguments
        class(SoilData), intent(inout) :: soil             ! Soil object
        real,            intent(in)    :: canopy_bh        ! Canopy base height (m)
        real,            intent(in)    :: canopy_bd        ! Canopy bulk density (kg/m3)
        real,            intent(in)    :: canopy_biom      ! Canopy biomass (kg/m2)
        real,            intent(in)    :: wind_ms          ! Wind speed (m/s)
        real,            intent(in)    :: ffmc             ! Fine fuel moisture code
        real,            intent(out)   :: R_a              ! Critical rate of spread for active crown fire (m/min)
        real,            intent(out)   :: rosf_active      ! Rate of spread for active crown fire (m/min)
        real,            intent(out)   :: CFB              ! Crown fraction burnt (0-1)
        real,            intent(out)   :: R_final          ! Final rate of spread (surface + crown) (m/min)
        real,            intent(out)   :: I_final          ! Final fire intensity (surface + crown) (kW/m)
        logical,         intent(out)   :: passive_crowning ! Do we have a passive crown fire?
        logical,         intent(out)   :: active_crowning  ! Do we have an active crown fire?

        ! Data dictionary: local variables
        real, dimension(FL_LEVS) :: tau_b       ! Temporary variable for calculating tau_l
        real                     :: I_p         ! Critical surface fire intensity for passive crowning (kW/m)
        real                     :: R_p         ! Critical ROS for passive crowning (m/min)
        real                     :: C_p         ! Check for passive crowning
        real                     :: C_a         ! Check for active crowning
        real                     :: rosf_SA     ! Rate of spread at which R_active = R'_active
        real                     :: tau_l       ! Residence time of fire (min)
        integer                  :: ic          ! Looping index

        ! Initialize to false
        passive_crowning = .false.
        active_crowning = .false.

        ! Check for passive crowning
        if (canopy_biom > 0.0) then
            I_p = ((canopy_bh*(460.0 + 25.9*FMC))/(100.0))**1.5
            C_p = soil%I_surf/I_p
        else
            I_p = 0.0
            C_p = 0.0
        end if

        ! Rate of spread of passive crowning
        R_p = 60.0*I_p/(H*soil%fuel_consumed)

        if (C_p >= 1.0) then

            passive_crowning = .true.

            ! Now check for active crowning

            ! Critical rate of spread for active crowning (m/min)
            R_a = S_0/canopy_bd*60.0

            ! Calculate active rate of spread (m/min)
            call active_rate_of_spread(soil, wind_ms, ffmc, R_a, rosf_active,  &
                rosf_SA)

            C_a = rosf_active/R_a

            if (C_a >= 1.0) then

                active_crowning = .true.

                ! Update residence time to include canopy fuels
                do ic = 1, FL_LEVS
                    tau_b(ic) = 39.4*(soil%fuels(ic)/10.0)*                    &
                    (1.0 - ((1.0 - soil%frac_burnt(ic))**0.50))
                end do
                tau_b(FLBL) = 0.0 ! Boles not included
                tau_b(FLRT) = 0.0 ! Roots not included

                tau_l = sum(tau_b) + 39.4*((canopy_biom/10.0)*CFB)

                soil%tau_l = tau_l

                ! Calculate crown fraction burned
                CFB = (soil%rosf - R_p)/(rosf_SA - R_p)
                if (CFB > 1.0) CFB = 1.0
                if (CFB < 0.0) CFB = 0.0
            else
                CFB = 0.0
            end if
        else
            R_a = RNVALID
            rosf_active = 0.0
            rosf_SA = RNVALID
            CFB = 0.0
        end if

        if (active_crowning) then
            R_final = soil%rosf + CFB*(rosf_active - soil%rosf)
            I_final = ((H*soil%fuel_consumed + (canopy_biom*H*CFB))*R_final)/60.0
        else
            R_final = soil%rosf
            I_final = soil%I_surf
        end if



    end subroutine active_passive_check

    !:.........................................................................:

    subroutine update_fuels(soil)
        !
        !  Update fuel conditions from fresh litter
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    02/15/21    A. C. Foster         Original Code
        !

        ! Data dicitonary: calling arguments
        class(SoilData), intent(inout) :: soil ! Soil object

        ! Data dictionary: local variables
        real, dimension(FL_LEVS) :: fuel_amnt   ! Litter amount in each type (kg/m2)
        real                     :: fuel_frac   ! Fraction of fuel
        real                     :: litter_amnt ! Amount of fuel (kg/m2)
        integer                  :: lc          ! Looping index

        ! Start with just existing
        fuel_amnt = soil%fuels

        ! Update BD and SAV from fresh litter
        do lc = 1, LIT_LEVS

            ! Litter amount (kg/m2)
            litter_amnt = soil%litter(lc)*T_TO_KG/HEC_TO_M2
            if (litter_amnt > epsilon(1.0)) then

                if (lc /= 6 .and. lc < 10) then

                    ! Deciduous leaf litter
                    fuel_frac = litter_amnt/(fuel_amnt(FLDEC) + litter_amnt)
                    soil%litBD(FLDEC) = soil%litBD(FLDEC)*(1 - fuel_frac) +    &
                        litter_params(lc, 11)*fuel_frac
                    soil%litSAV(FLDEC) = soil%litSAV(FLDEC)*(1 - fuel_frac) +  &
                        litter_params(lc, 12)*fuel_frac
                    fuel_amnt(FLDEC) = fuel_amnt(FLDEC) + litter_amnt

                else if (lc == 6 .or. (lc >= 10 .and. lc <= 12)) then

                    ! Coniferous needle litter
                    fuel_frac = litter_amnt/(fuel_amnt(FLCON) + litter_amnt)
                    soil%litBD(FLCON) = soil%litBD(FLCON)*(1 - fuel_frac) +    &
                        litter_params(lc, 11)*fuel_frac
                    soil%litSAV(FLCON) = soil%litSAV(FLCON)*(1 - fuel_frac) +  &
                        litter_params(lc, 12)*fuel_frac
                    fuel_amnt(FLCON) = fuel_amnt(FLCON) + litter_amnt

                else if (lc == IROOT) then

                    ! Root litter
                    fuel_frac = litter_amnt/(fuel_amnt(FLRT) + litter_amnt)
                    soil%litBD(FLRT) = soil%litBD(FLRT)*(1 - fuel_frac) +      &
                        litter_params(lc, 11)*fuel_frac
                    soil%litSAV(FLRT) = soil%litSAV(FLRT)*(1 - fuel_frac) +    &
                        litter_params(lc, 12)*fuel_frac
                    fuel_amnt(FLRT) = fuel_amnt(FLRT) + litter_amnt

                else if (lc == ISBL .or. lc == ISBR) then

                    ! Bole litter
                    fuel_frac = litter_amnt/(fuel_amnt(FLBL) + litter_amnt)
                    soil%litBD(FLBL) = soil%litBD(FLBL)*(1 - fuel_frac) +      &
                        litter_params(lc, 11)*fuel_frac
                    soil%litSAV(FLBL) = soil%litSAV(FLBL)*(1 - fuel_frac) +    &
                        litter_params(lc, 12)*fuel_frac
                    fuel_amnt(FLBL) = fuel_amnt(FLBL) + litter_amnt

                else if (lc == ITW) then

                    ! Twig litter
                    fuel_frac = litter_amnt/(fuel_amnt(FLTW) + litter_amnt)
                    soil%litBD(FLTW) = soil%litBD(FLTW)*(1 - fuel_frac) +      &
                        litter_params(lc, 11)*fuel_frac
                    soil%litSAV(FLTW) = soil%litSAV(FLTW)*(1 - fuel_frac) +    &
                        litter_params(lc, 12)*fuel_frac
                    fuel_amnt(FLTW) = fuel_amnt(FLTW) + litter_amnt

                else if (lc == ISBR) then

                    ! Small branch litter
                    fuel_frac = litter_amnt/(fuel_amnt(FLSBR) + litter_amnt)
                    soil%litBD(FLSBR) = soil%litBD(FLSBR)*(1 - fuel_frac) +    &
                        litter_params(lc, 11)*fuel_frac
                    soil%litSAV(FLSBR) = soil%litSAV(FLSBR)*(1 - fuel_frac) +  &
                        litter_params(lc, 12)*fuel_frac
                    fuel_amnt(FLSBR) = fuel_amnt(FLSBR) + litter_amnt

                else if (lc == ILBR) then

                    ! Large branch litter
                    fuel_frac = litter_amnt/(fuel_amnt(FLLBR) + litter_amnt)
                    soil%litBD(FLLBR) = soil%litBD(FLLBR)*(1 - fuel_frac) +    &
                        litter_params(lc, 11)*fuel_frac
                    soil%litSAV(FLLBR) = soil%litSAV(FLLBR)*(1 - fuel_frac) +  &
                        litter_params(lc, 12)*fuel_frac
                    fuel_amnt(FLLBR) = fuel_amnt(FLLBR) + litter_amnt

                else if (lc == IMOSS) then
                    ! Moss litter
                    fuel_frac = litter_amnt/(fuel_amnt(FLDM) + litter_amnt)
                    soil%litBD(FLDM) = soil%litBD(FLDM)*(1 - fuel_frac) +      &
                        litter_params(lc, 11)*fuel_frac
                    soil%litSAV(FLDM) = soil%litSAV(FLDM)*(1 - fuel_frac) +    &
                        litter_params(lc, 12)*fuel_frac
                    fuel_amnt(FLDM) = fuel_amnt(FLDM) + litter_amnt

                end if
            end if
        end do

        soil%fuels = fuel_amnt

        ! Update total litter - don't include boles or root litter
        soil%fuel_sum = sum(soil%fuels) - soil%fuels(FLBL) - soil%fuels(FLRT)

    end subroutine update_fuels

    !:.........................................................................:

    subroutine fuel_conditions(soil, ffmc, dmc, FDI, siteID, plot, year, day)
        !
        !  Calculates daily fuel conditions for calculation of reate of spread
        !  and fire intensity as in Thonicke et al. (2010) and the Canadian
        !  Fire Rating System
        !
        ! Relative total litter moisture is calculated from weighted average of
        ! relative moisture content of 1-h, 10-h, and 100-h fuels

        ! In Thonicke et al. (2010):
        !  1-h fuels: leaves, twigs, live herbaceous leaves
        !  10-h fuels: small branches
        !  100-h fuels: large branches
        !
        ! For our purposes we will try to co-opt these classes
        !  1-h fuels: leaf litter, moss litter, live moss, live shrub foliage
        !             13.6% of branch litter and live shrub twigs
        !  10-h fuels: 22.3% of branch litter
        !  100-h fuels: 63.4% of branch litter
        !
        ! Fuel array:
        ! 1: deciduous leaf litter
        ! 2: conifer needle litter
        ! 3: twig litter
        ! 4: small branch litter
        ! 5: large branch litter
        ! 6: bole litter (not included in ROS calculations)
        ! 7: live moss
        ! 8: moss litter
        ! 9: root litter (not included in ROS calculations)
        ! 10: live shrub foliage and fine twigs
        !
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    03/10/19    A. C. Foster         Original Code
        !

        ! Data dictionary: constants
        real, parameter :: MOSS_MEF = 0.76 ! Moss moisture of extinction

        ! Data dictionary: calling arguments
        class(SoilData), intent(inout) :: soil    ! Soil object
        real,            intent(in)    :: ffmc    ! Fine fuel moisture code
        real,            intent(in)    :: dmc     ! Duff moisture code
        real,            intent(out)   :: FDI     ! Fire danger index (0-1)
        integer,         intent(in)    :: siteID  ! Site ID
        integer,         intent(in)    :: plot    ! Plot number
        integer,         intent(in)    :: year    ! Simulation year
        integer,         intent(in)    :: day     ! Simulation day

        ! Data dictionary: local variables
        real, dimension(FL_LEVS) :: fuel_frac      ! Fraction of total fuel loading
        real, dimension(FL_LEVS) :: alpha_FMC      ! Drying rate
        real, dimension(FL_LEVS) :: m_e            ! Moisture of extinction (volumetric)
        real, dimension(FL_LEVS) :: fuel_moist     ! Fuel moisture (volumetric)
        real, dimension(FL_LEVS) :: fuel_eff_moist ! Relative fuel moisture
        real                     :: MEF            ! Average moisture of extinction (volumetric)
        real                     :: sumlit_SAV     ! Average SAV (/cm)
        real                     :: sumlit_BD      ! Average bulk density (kg/m3)
        real                     :: sumlit_moist   ! Average moisture (volumetric)
        integer                  :: lc             ! Looping index

        ! Initialize accumulators
        MEF = 0.0
        sumlit_moist = 0.0
        sumlit_SAV = 0.0
        sumlit_BD = 0.0

        if (soil%fuel_sum > epsilon(1.0)) then
            do lc = 1, FL_LEVS

                ! Don't include boles or roots
                if (lc /= FLBL .and. lc /= FLRT .and. lc /= FLLM .and.         &
                    lc /= FLSH) then

                    ! Calculate fuel fraction, drying rate, moisture of
                    ! extinction, and moisture
                    if (soil%fuels(lc) > epsilon(1.0)) then
                        fuel_frac(lc) = soil%fuels(lc)/soil%fuel_sum
                        alpha_FMC(lc) = exp(-1.0*(soil%litSAV(lc)/250.0))
                        m_e(lc) = 0.524 - 0.066*log(soil%litSAV(lc))
                        fuel_moist(lc) = (147.2*(101.0 - ffmc))/               &
                            (59.5 + ffmc)/100.0*alpha_FMC(lc)
                    else
                        fuel_frac(lc) = 0.0
                        alpha_FMC(lc) = 0.0
                        m_e(lc) = 0.0
                        fuel_moist(lc) = 0.0
                    end if
                else if (lc == FLBL .or. lc == FLRT) then
                    ! Boles and roots
                    if (soil%fuels(lc) > epsilon(1.0)) then
                        fuel_frac(lc) = soil%fuels(lc)/soil%fuel_sum
                        alpha_FMC(lc) = exp(-1.0*(soil%litSAV(lc)/250.0))
                        m_e(lc) = 0.524 - 0.066*log(soil%litSAV(lc))
                        fuel_moist(lc) = (20.0 + 100.0/exp(0.023*dmc))/100.0*  &
                            alpha_FMC(lc)
                    else
                        fuel_frac(lc) = 0.0
                        alpha_FMC(lc) = 0.0
                        m_e(lc) = 0.0
                        fuel_moist(lc) = 0.0
                    end if
                else if (lc == FLLM) then
                    ! Live moss
                    if (soil%fuels(lc) > epsilon(1.0)) then
                        fuel_frac(lc) = soil%fuels(lc)/soil%fuel_sum
                        alpha_FMC(lc) = 1.0
                        fuel_moist(lc) = ((199.2*(101.0 - ffmc))/              &
                            ( 59.5 + ffmc) + 12.0)/100.0
                    else
                        fuel_frac(lc) = 0.0
                        alpha_FMC(lc) = 0.0
                        fuel_moist(lc) = 0.0
                    end if
                else if (lc == FLSH) then
                    ! Live shrubs
                    if (soil%fuels(lc) > epsilon(1.0)) then
                        fuel_frac(lc) = soil%fuels(lc)/soil%fuel_sum
                        alpha_FMC(lc) = 1.0
                        m_e(lc) = 100.0
                        fuel_moist(lc) = (90.0 + 30.0/exp(0.023*dmc))/100.0
                    else
                        fuel_frac(lc) = 0.0
                        alpha_FMC(lc) = 0.0
                        m_e(lc) = 0.0
                        fuel_moist(lc) = 0.0
                    end if
                end if
            end do

            ! Moss
            m_e(FLDM) = MOSS_MEF
            m_e(FLLM) = MOSS_MEF

            ! Calculate weighted average fuel conditions
            do lc = 1, FL_LEVS
                if (lc /= FLBL .and. lc /= FLRT .and.                          &
                    soil%fuels(lc) >= epsilon(1.0)) then
                    MEF = MEF + fuel_frac(lc)*m_e(lc)
                    sumlit_moist = sumlit_moist + fuel_frac(lc)*fuel_moist(lc)
                    sumlit_SAV = sumlit_SAV + fuel_frac(lc)*soil%litSAV(lc)
                    sumlit_BD = sumlit_BD + fuel_frac(lc)*soil%litBD(lc)
                end if

                ! Calculate relative moisture content
                if (m_e(lc) > 0.0 .and. soil%fuels(lc) >= epsilon(1.0)) then
                    fuel_eff_moist(lc) = fuel_moist(lc)/m_e(lc)
                else
                    fuel_eff_moist(lc) = 0.0
                end if
            end do

            ! Calculate fire danger index
            FDI = max(0.0, (1.0 - (1.0/MEF)*sumlit_moist))

        else
            ! No fuels
            fuel_frac = 0.0
            alpha_FMC = 0.0
            m_e = 0.0
            fuel_moist = 0.0
            fuel_eff_moist = 0.0
        end if

        ! Save instance variables
        soil%MEF = MEF
        soil%fuel_moisture = sumlit_moist
        soil%fuel_BD = sumlit_BD
        soil%fuel_SAV = sumlit_SAV
        soil%litter_moist = fuel_moist
        soil%litter_mef = m_e

        if (conds_testing) then
            call csv_write(fuel_conds, siteID, .false.)
            call csv_write(fuel_conds, plot, .false.)
            call csv_write(fuel_conds, year, .false.)
            call csv_write(fuel_conds, day, .false.)
            call csv_write(fuel_conds, soil%fuel_sum, .false.)
            do lc = 1, FL_LEVS
                call csv_write(fuel_conds, soil%fuels(lc), .false.)
            end do
            do lc = 1, FL_LEVS
                call csv_write(fuel_conds, soil%litBD(lc), .false.)
            end do
            do lc = 1, FL_LEVS
                call csv_write(fuel_conds, soil%litSAV(lc), .false.)
            end do
            do lc = 1, FL_LEVS
                call csv_write(fuel_conds, soil%litter_moist(lc), .false.)
            end do
            do lc = 1, FL_LEVS
                call csv_write(fuel_conds, soil%litter_mef(lc), .false.)
            end do
            call csv_write(fuel_conds, sumlit_moist, .false.)
            call csv_write(fuel_conds, sumlit_SAV, .false.)
            call csv_write(fuel_conds, sumlit_BD, .false.)
            call csv_write(fuel_conds, MEF, .true.)
        end if

    end subroutine fuel_conditions

    !:.........................................................................:

    subroutine rate_of_spread(soil, wind_ms, ffmc, FDI, rosf, siteID, plot,    &
            year, day)
        !
        ! Calculates daily potential rate of spread as in Thonicke et al. 2010
        ! and Rothermel 1972
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    03/10/19    A. C. Foster         Original Code
        !

        ! Data dictionary: constants
        real, parameter :: MINER_TOTAL = 0.055 ! Total mineral fraction of litter
        real, parameter :: PART_DENS = 513.0   ! Oven-dry particle density (kg/m3)
        real, parameter :: N_S = 0.41739       ! Mineral dampening coefficient

        ! Data dictionary: calling arguments
        class(SoilData), intent(inout) :: soil    ! Soil object
        real,            intent(in)    :: wind_ms ! Wind speed (m/s)
        real,            intent(in)    :: ffmc    ! Fine fuel moisture content
        real,            intent(in)    :: FDI     ! Fire Danger Index
        integer,         intent(in)    :: siteID  ! Site ID
        integer,         intent(in)    :: plot    ! Plot number
        integer,         intent(in)    :: year    ! Simulation year
        integer,         intent(in)    :: day     ! Simulation day
        real,            intent(out)   :: rosf    ! Rate of forward spread (m/min)

        ! Data dictionary: local variables
        real :: wind      ! Wind speed (m/min)
        real :: net_fuel  ! Net fuel loading (kg/m2)
        real :: beta      ! Packing ratio
        real :: beta_op   ! Optimum packing ratio
        real :: beta_frac ! Fraction of beta:beta_op
        real :: rVol_max  ! Maximum reaction velocity (/min)
        real :: rVol_opt  ! Optimum reaction velocity (/min)
        real :: A         ! Temporary variable
        real :: mw_weight ! Relative moisture content
        real :: n_M       ! Moisture dampening coefficient
        real :: I_r       ! Reaction intensity (kJ/kg/m2)
        real :: flux_rat  ! Propagating flux ratio
        real :: B         ! Temporary variable for phi_wind
        real :: C         ! Temporary variable for phi_wind
        real :: E         ! Temporary variable for phi_wind
        real :: Uf        ! Effective wind speed (m/min)
        real :: Uf_max    ! Maximum wind speed (m/min)
        real :: phi_wind  ! Wind coefficient
        real :: eps       ! Effective heating number
        real :: Qig       ! Heat of pre-ignition (kJ/kg)
        real :: rosb      ! Rate of backwards spread (m/min)
        real :: Lb        ! Length-to-breadth ratio of fire shape (ellipse)
        real :: t_fire    ! Fire duration (min)
        real :: Dt        ! Length of major axis of fire (m)
        real :: Db        ! Length of major axis of fire from backwards spread (m)
        real :: Df        ! Length of major axis of fire from forward spread (m)
        real :: a_f       ! Fire area (ha)
        integer :: lc     ! Looping index

        ! Calculate wind speed in m/min (input is m/s)
        !wind = min(wind_ms, 14.0)*MIN_TO_SEC
        wind = wind_ms*MIN_TO_SEC

        ! Net fuel loading
        net_fuel = soil%fuel_sum*(1.0 - MINER_TOTAL)

        ! Calculate packing ratio
        beta = soil%fuel_BD/PART_DENS

        if (soil%fuel_sum <= epsilon(1.0)) then
            ! No fuel
            A = 0.0
            rVol_opt = 0.0
            n_M = 0.0
            rVol_max = 0.0
            beta_op = 0.0
            beta_frac = 0.0
        else
            ! Optimum packing ratio
            beta_op = 0.200395*(soil%fuel_SAV**(-0.8189))
            beta_frac = beta/beta_op

            ! Maximum reaction velocity (/min)
            rVol_max = 1.0/(0.0591 + 2.926*(soil%fuel_SAV**(-1.5)))

            ! Optimum reaction velocity (/min)
            A = 8.9033*(soil%fuel_SAV**(-0.7913))
            rVol_opt = rVol_max*(beta_frac**A)*exp(A*(1.0 - beta_frac))

            ! Calculate moisture dampening coefficient
            mw_weight = soil%fuel_moisture/soil%MEF
            n_M = max(0.0, 1.0 - 2.59*mw_weight + 5.11*(mw_weight**2) -        &
                3.52*(mw_weight**3))
        end if

        ! Reaction intensity (kJ/m2/min)
        I_r = rVol_opt*net_fuel*H*n_M*N_S

        ! Propagating flux ratio
        flux_rat = (exp((0.792 + 3.7597*sqrt(soil%fuel_SAV))*(beta + 0.1)))/   &
            (192.0 + 7.9095*soil%fuel_SAV)

        ! phi_wind parameters
        B = 0.15988*(soil%fuel_SAV**0.54)
        C = 7.47*exp(-0.8711*(soil%fuel_SAV**0.55))
        E = 0.715*exp(-0.01094*soil%fuel_SAV)

        ! Effective wind speed (m/min)
        Uf = wind*0.3

        ! Maximum effective wind speed
        Uf_max = 13.12629*I_r**(1.0/3.0)

        ! Cap Uf maximum
        if (Uf > Uf_max) Uf = Uf_max

        if (soil%fuel_sum <= epsilon(1.0)) then
            phi_wind = 0.0
        else
            ! Wind effect on flux ratio
            phi_wind = C*((3.281*Uf)**B)*(beta_frac**(-1.0*E))
        end if

        ! Effective heating number
        eps = exp(-4.528/soil%fuel_SAV)

        ! Heat of pre-ignition (kJ/kg)
        Qig = 581.0 + 2594.0*soil%fuel_moisture

        if (soil%fuel_BD <= 0.0 .or. eps <= 0.0 .or. Qig <= 0.0) then
            rosf = 0.0
            rosb = 0.0
            soil%rosf_phi = 0.0
        else
            ! Rate of forward spread (m/min)
            rosf = (I_r*flux_rat*(1.0 + phi_wind))/(soil%fuel_BD*eps*Qig)
            soil%rosf_phi = (I_r*flux_rat)/(soil%fuel_BD*eps*Qig)

            ! Rate of backwards spread (m/min)
            rosb = rosf*exp(-0.012*Uf)
        end if

        ! Length-to-breadth ratio of ellipse
        Lb = 1.0 + 8.729*((1.0 - exp(-0.03*Uf))**2.115)

        ! Fire duration (min)
        t_fire = 241.0/(1 + 240.0*exp(-11.06*FDI))

        ! Length of major axis (m):
        Df = rosf*t_fire
        Db = rosb*t_fire
        Dt = Df + Db

        ! Area of fire (ha)
        a_f = ((PI/4.0*Lb)*(Dt**2))/HEC_TO_M2

        soil%rosf = rosf


         if (testing) then
             call csv_write(dayfire, siteID, .false.)
             call csv_write(dayfire, plot, .false.)
             call csv_write(dayfire, year, .false.)
             call csv_write(dayfire, day, .false.)
             call csv_write(dayfire, FDI, .false.)
             call csv_write(dayfire, ffmc, .false.)
             call csv_write(dayfire, soil%MEF, .false.)
             call csv_write(dayfire, soil%fuel_moisture, .false.)
             call csv_write(dayfire, soil%fuel_BD, .false.)
             call csv_write(dayfire, soil%fuel_SAV, .false.)
             do lc = 1, FL_LEVS
                 call csv_write(dayfire, soil%fuels(lc), .false.)
             end do
             call csv_write(dayfire, Uf, .false.)
             call csv_write(dayfire, A, .false.)
             call csv_write(dayfire, beta, .false.)
             call csv_write(dayfire, beta_op, .false.)
             call csv_write(dayfire, beta_frac, .false.)
             call csv_write(dayfire, rVol_max, .false.)
             call csv_write(dayfire, rVol_opt, .false.)
             call csv_write(dayfire, n_M, .false.)
             call csv_write(dayfire, net_fuel, .false.)
             call csv_write(dayfire, I_r, .false.)
             call csv_write(dayfire, flux_rat, .false.)
             call csv_write(dayfire, phi_wind, .false.)
             call csv_write(dayfire, Qig, .false.)
             call csv_write(dayfire, rosf, .false.)
             call csv_write(dayfire, a_f, .false.)
         end if

    end subroutine rate_of_spread

    !:.........................................................................:

    subroutine active_rate_of_spread(soil, wind_ms, ffmc, R_a, rosf_active,    &
            rosf_SA)
        !
        ! Calculates active rate of spread based on Rothermel 1991
        !
        ! This uses the standard Rothermel rate of spread equation using
        ! fuel parameters from Fire Behavior Model (FM) 10 - Anderson 1982 and
        ! a mid-flame wind-speed 40% of the measured wind speed.
        !
        ! FM 10 loading, MEF, and depth from Anderson 1982, "Aids to determining
        ! fuel models for estimating fire behavior"
        !
        ! SAV values from the BEHAVE model (Burgan & Rothermel 1984)
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/24/21    A. C. Foster         Original Code
        !

        ! Data dictionary: constants
        real, parameter :: MINER_TOTAL = 0.055 ! Total mineral fraction of litter
        real, parameter :: PART_DENS = 513.0   ! Oven-dry particle density (kg/m3)
        real, parameter :: N_S = 0.41739       ! Mineral dampening coefficient
        real, parameter :: FUEL_1HR = 3.01     ! FM 10 1-hr fuel loading (US tons/acre)
        real, parameter :: FUEL_10HR = 2.0     ! FM 10 10-hr fuel loading (US tons/acre)
        real, parameter :: FUEL_100HR = 5.01   ! FM 10 100-hr fuel loading (US tons/acre)
        real, parameter :: FUEL_LIVE = 2.0     ! FM 10 live fuel loading (US tons/acre)
        real, parameter :: FUEL_MEF = 0.25     ! FM 10 moisture of extinction (volumetric)
        real, parameter :: FUEL_DEPTH = 1.0    ! FM 10 fuel depth (ft)
        real, parameter :: SAV_1HR = 2000.0    ! FM 10 1-hr SAV (ft2/ft3)
        real, parameter :: SAV_10HR = 109.0    ! FM 10 10-hr SAV (ft2/ft3)
        real, parameter :: SAV_100HR = 30.0    ! FM 10 100-hr SAV (ft2/ft3)
        real, parameter :: SAV_LIVE = 1650.0   ! FM 10 live SAV (ft2/ft3)

        ! Data dictionary: calling arguments
        class(SoilData), intent(inout) :: soil         ! Soil object
        real,            intent(in)    :: wind_ms      ! Wind speed (m/s)
        real,            intent(in)    :: ffmc         ! Fine fuel moisture code
        real,            intent(in)    :: R_a          ! Critical active rate of spread
        real,            intent(out)   :: rosf_active  ! Rate of forward spread (m/min)
        real,            intent(out)   :: rosf_SA      ! Rate of spread at conditions where R_active = R'_active

        ! Data dictionary: local variables
        real :: total_fuel ! Total fuel loading (US tons/acre)
        real :: net_fuel   ! Net fuel loading (kg/m2)
        real :: fuel_BD    ! Fuel bulk density (kg/m3)
        real :: fuel_SAV   ! Fuel SAV (/cm)
        real :: beta       ! Packing ratio
        real :: beta_op    ! Optimum packing ratio
        real :: beta_frac  ! Fraction of beta:beta_op
        real :: rVol_max   ! Maximum reaction velocity (/min)
        real :: rVol_opt   ! Optimum reaction velocity (/min)
        real :: A          ! Temporary variable
        real :: moisture   ! Fuel moisture (volumetric)
        real :: mw_weight  ! Relative moisture content
        real :: n_M        ! Moisture dampening coefficient
        real :: I_r        ! Reaction intensity (kJ/kg/m2)
        real :: flux_rat   ! Propagating flux ratio
        real :: B          ! Temporary variable for phi_wind
        real :: C          ! Temporary variable for phi_wind
        real :: E          ! Temporary variable for phi_wind
        real :: Uf         ! Wind speed (m/min)
        real :: phi_wind   ! Wind coefficient
        real :: eps        ! Effective heating number
        real :: Qig        ! Heat of pre-ignition (kJ/kg)
        real :: Uf_check   ! Wind speed at which R_active = R'_active

        ! Calculate wind speed in m/min (input is m/s)
        ! 40% of measured
        Uf = (wind_ms*MIN_TO_SEC)*0.4

        ! Net fuel loading (kg/m2) for FM 10
        total_fuel = FUEL_1HR + FUEL_10HR + FUEL_100HR + FUEL_LIVE
        net_fuel = (total_fuel*907.185/4046.86)*(1.0 - MINER_TOTAL)

        ! FM 10 bulk density (kg/m3)
        fuel_BD = (total_fuel*907.185/4046.86)/(FUEL_DEPTH/3.28084)

        ! FM 10 SAV (/cm)
        fuel_SAV = (FUEL_1HR/total_fuel*SAV_1HR +                              &
            FUEL_10HR/total_fuel*SAV_10HR + FUEL_100HR/total_fuel*SAV_100HR +  &
            FUEL_LIVE/total_fuel*SAV_LIVE)*929.03/28316.8

        ! Calculate packing ratio
        beta = fuel_BD/PART_DENS

        ! Optimum packing ratio
        beta_op = 0.200395*(fuel_SAV**(-0.8189))
        beta_frac = beta/beta_op

        ! Maximum reaction velocity (/min)
        rVol_max = 1.0/(0.0591 + 2.926*(fuel_SAV**(-1.5)))

        ! Optimum reaction velocity (/min)
        A = 8.9033*(fuel_SAV**(-0.7913))
        rVol_opt = rVol_max*(beta_frac**A)*exp(A*(1.0 - beta_frac))

        ! Calculate moisture dampening coefficient
        moisture = (147.2*(101.0 - ffmc))/(59.5 + ffmc)/100.0
        mw_weight = moisture/FUEL_MEF
        n_M = max(0.0, 1.0 - 2.59*mw_weight + 5.11*(mw_weight**2) -            &
            3.52*(mw_weight**3))

        ! Reaction intensity (kJ/m2/min)
        I_r = rVol_opt*net_fuel*H*n_M*N_S

        ! Propagating flux ratio
        flux_rat = (exp((0.792 + 3.7597*sqrt(fuel_SAV))*(beta + 0.1)))/        &
            (192.0 + 7.9095*fuel_SAV)

        ! phi_wind coefficients
        B = 0.15988*(fuel_SAV**0.54)
        C = 7.47*exp(-0.8711*(fuel_SAV**0.55))
        E = 0.715*exp(-0.01094*fuel_SAV)

        ! Wind effect
        phi_wind = C*((3.281*Uf)**B)*(beta_frac**(-1.0*E))

        ! Effective heating number
        eps = exp(-4.528/fuel_SAV)

        ! Heat of pre-ignition (kJ/kg)
        Qig = 581.0 + 2594.0*moisture

        ! Active rate spread (m/min)
        rosf_active = (I_r*flux_rat*(1.0 + phi_wind))/(fuel_BD*eps*Qig)

        Uf_check = ((((R_a*fuel_BD*eps*Qig)/(I_r*flux_rat) - 1.0)/             &
            (C*(beta_frac**(-1.0*E))))**(1.0/B))/3.281

        ! Recalculate phi_wind with actual fuel characteristics
        B = 0.15988*(soil%fuel_SAV**0.54)
        C = 7.47*exp(-0.8711*(soil%fuel_SAV**0.55))
        E = 0.715*exp(-0.01094*soil%fuel_SAV)
        beta = soil%fuel_BD/PART_DENS
        beta_op = 0.200395*(soil%fuel_SAV**(-0.8189))
        beta_frac = beta/beta_op
        phi_wind = C*((3.281*Uf_check)**B)*(beta_frac**(-1.0*E))

        !Surface rate of spread under wind conditions where R_active = R'_active
        rosf_SA = soil%rosf_phi*(1.0 + phi_wind)

    end subroutine active_rate_of_spread

    !:.........................................................................:

    subroutine fire_intensity(soil, rosf, I_surf)
        !
        ! Calculates daily fire intensity as in  Thonicke et al. 2010
        ! and Rothermel 1972
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    03/10/19    A. C. Foster         Original Code
        !

        ! Data dictionary: constants
        real,                     parameter :: MINER_TOTAL = 0.055
        real, dimension(FL_LEVS), parameter :: MIN_MOISTURE = [0.24, 0.24,     &
            0.18, 0.12, 0.0, 0.0, 0.45, 0.45, 0.0, 0.0]
        real, dimension(FL_LEVS), parameter :: PARA = [-1.36, -1.36, -1.32,    &
            -0.83, -0.17, 0.06, -1.72, -1.72, 0.06, -1.72]
        real, dimension(FL_LEVS), parameter :: PARB = [0.54, 0.54, 0.41,       &
            -0.22, -0.82, -0.87, 0.57, 0.57, -0.87, 0.57]
        real, dimension(FL_LEVS), parameter :: PARC = [0.95, 0.95, 0.96, 1.02, &
            0.98, 0.55, 0.98, 0.98, 0.82, 0.98]

        ! Data dictionary: calling arguments
        class(SoilData), intent(inout) :: soil   ! Soil object
        real,            intent(in)    :: rosf   ! Rate of forward spread of fire (m/min)
        real,            intent(out)   :: I_surf ! Surface fire intensity (kW/m)

        ! Data dictionary: local variables
        real,   dimension(FL_LEVS) :: tau_b         ! Temporary variable
        real,   dimension(FL_LEVS) :: burned        ! Amount consumed (kg/m2)
        real                       :: fuel_consumed ! Total amount consumed (kg/m2)
        real                       :: lmoist        ! Relative fuel moisture
        real                       :: tau_l         ! Residence time of fire (min)
        integer                    :: ic            ! Looping index

        ! 0 values
        soil%frac_burnt = 0.0

        ! Calculate potential fuel consumption from surface fuels

        do ic = 1, FL_LEVS
            if (soil%fuels(ic) > epsilon(1.0)) then

                lmoist = soil%litter_moist(ic)/soil%litter_mef(ic)

                if (lmoist <= MIN_MOISTURE(ic)) then
                    ! Dry litter
                    soil%frac_burnt(ic) = 1.0
                    burned(ic) = soil%fuels(ic)

                else if (lmoist > MIN_MOISTURE(ic) .and. lmoist <= 1.0) then

                    soil%frac_burnt(ic) = max(0.0, min(1.0, PARA(ic)*lmoist**2 +   &
                        PARB(ic)*lmoist + PARC(ic)))
                    burned(ic) = soil%frac_burnt(ic)*soil%fuels(ic)

                else if (lmoist > 1.0) then
                    ! High moisture
                    soil%frac_burnt(ic) = 0.0
                    burned(ic) = 0.0
                end if
            end if
        end do

        soil%frac_burnt(:) = soil%frac_burnt(:)*(1.0 - MINER_TOTAL)
        fuel_consumed = (sum(burned) - burned(FLBL) - burned(FLRT))

        ! Intensity of surface fire (kW/m)
        I_surf = H*fuel_consumed*(rosf/60.0)

        ! Fire residence time (min)
        do ic = 1, FL_LEVS
            tau_b(ic) = 39.4*(soil%fuels(ic)/10.0)*                            &
               (1.0 - ((1.0 - soil%frac_burnt(ic))**0.50))
        end do
        tau_b(FLBL) = 0.0 ! Don't include boles
        tau_b(FLRT) = 0.0 ! Don't include roots
        tau_l = sum(tau_b)

        if (tau_l > 10.0) tau_l = 10.0

        soil%tau_l = tau_l
        soil%I_surf = I_surf
        soil%fuel_consumed = fuel_consumed

        if (testing) then
             do ic = 1, FL_LEVS
                 call csv_write(dayfire, soil%frac_burnt(ic), .false.)
             end do
             call csv_write(dayfire, I_surf, .false.)
             call csv_write(dayfire, tau_l, .true.)
        end if

    end subroutine fire_intensity

    !:.........................................................................:

    subroutine fuel_consumption(siteID, ip, year, soil, consRoot, N_cons)
        !
        ! Calculates consumption of litter and humus layers from fire
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    03/10/19    A. C. Foster         Original Code
        !

        ! Data dictionary: calling arguments
        class(SoilData), intent(inout) :: soil     ! Soil object
        integer,         intent(in)    :: siteID   ! Site ID
        integer,         intent(in)    :: ip       ! Plot number
        integer,         intent(in)    :: year     ! Simulation year
        real,            intent(out)   :: consRoot ! Proportion roots consumed by fire (0-1)
        real,            intent(out)   :: N_cons   ! Proportion of N consumed by fire (0-1)

        ! Data dictionary: local variables
        real, dimension(NCOH_MAX, NCOH_CHR) :: C           ! Array of cohorts
        real                                :: sfan        ! Volatilized N from fires (tN/ha)
        real                                :: pre_depth   ! Pre-fire organic layer depth (m)
        real                                :: duff_moist  ! Duff moisture content (volumetric)
        real                                :: rfs         ! Relative duff moisture content
        real                                :: volN        ! Proportion N volatilized by fires
        real                                :: emis        ! Duff emissivity
        real                                :: duff_cons   ! Duff consumption through smoldering (kg/m2)
        real                                :: hum_avail   ! Humus content (kg/m2)
        real                                :: root_avail  ! Root content (kg/m2)
        real                                :: root_cons   ! Root consumption
        real                                :: humloss     ! Proportion humus consumed by fire
        real                                :: Nloss       ! N content consumed by fire (tN/ha)
        real                                :: frac_loss   ! Fraction of litter consumed by fire
        real                                :: WDWcons     ! Fraction of well-decayed wood consumed by fire
        real                                :: weightloss  ! Litter amount consumed by fire (t/ha)
        real                                :: m_loss      ! Live moss consumed by fire (kg)
        real                                :: con_fuel    ! Conifer needle litter amount (kg)
        real                                :: dec_fuel    ! Deciduous leaf litter amount (kg)
        real                                :: twig_fuel   ! Twig litter amount (kg)
        real                                :: moss_fuel   ! Moss litter amount (kg)
        real                                :: duff        ! Duff amount (kg)
        real                                :: bg_combust  ! Belowground combustion (kg/m2)
        real                                :: agw_combust ! Aboveground woody combustion (kg/m2)
        real                                :: agw_prefire ! Aboveground woody prefire (kg/m2)
        integer                             :: nc          ! Number of decaying cohorts
        integer                             :: i, j, ix    ! Looping indices
        integer                             :: lt          ! Litter class

        ! Initialize accumulators
        sfan = 0.0
        soil%forest_litter = 0.0
        bg_combust = 0.0
        agw_combust = 0.0
        agw_prefire = 0.0

        ! Reduce table lookups
        nc = soil%ncohort
        do i = 1, nc
            do j = 1, NCOH_CHR
                C(i, j) = soil%cohorts(i, j)
            end do
        end do

        ! Pre-fire organic layer depth (m)
        pre_depth = soil%O_depth

        ! Get fuel moisture conditions
        duff_moist = (exp(-1.0*((soil%dmc_fire - 244.72)/43.43)) + 20.0)/100.0
        rfs = max(0.0, (3.0 - duff_moist)/(3.0 - soil%pwp(1)))

        ! Duff emmisivity
        emis = 0.770 - 0.004*(duff_moist*100.0)
        if (emis <= 0.1) emis = 0.1

        ! Set % nitrogen in forest floor consumed
        N_cons = max(0.0, min(0.6426*rfs + 3.34*(soil%M_depth +               &
            soil%O_depth), 0.7))
        volN = 1.0 - N_cons

        ! Calculate root consumption
        !consRoot = max(0.0, min((0.302 + 0.597*rfs + 3.34*(soil%M_depth +      &
        !   soil%O_depth)), 0.3))
        root_avail = soil%forest_litter(IROOT, 1)*T_TO_KG/HEC_TO_M2

        if (root_avail < epsilon(1.0)) then
            root_cons = 0.0
            consRoot = 0.0
        else
            root_cons = min(max(((60.0*soil%tau_l)*7.473*emis)/                &
                (110.0 + 6.2*(duff_moist*100.0)), 0.0), root_avail)
            consRoot = max(0.0, min(1.0, root_cons/root_avail))
        end if

        ! Well-decayed wood consumption
        WDWcons = max(0.0, min((0.098 + 0.597*rfs +                            &
            3.34*(soil%M_depth + soil%O_depth)), 0.5))

        ! Available duff = humus (kg/m2)
        hum_avail = C(1, 1)*T_TO_KG/HEC_TO_M2

        ! Duff consumption through smoldering (kg/m2)
        duff_cons = min(max(((60.0*soil%tau_l)*7.473*emis)/                    &
            (110.0 + 6.2*(duff_moist*100.0)), 0.0), hum_avail)

        ! Calculate humus loss (proportion) and update pool
        if (hum_avail < epsilon(1.0)) then
            humloss = 0.0
        else
            humloss = max(0.0, min(1.0, duff_cons/hum_avail))
            C(1, 1) = max(C(1, 1) - duff_cons/T_TO_KG/M2_TO_HEC, 0.0)
        end if

        ! Calculate N loss from humus and update pool
        Nloss = C(1, 2)*humloss
        C(1, 2) = max(C(1, 2) - Nloss, 0.0)
        sfan = sfan + Nloss*volN ! Some is volatilized

        WDWcons = humloss

        ! Get rid of litter cohorts
        do i = 2, nc

            ! Litter type
            lt = int(C(i, 5))

            if (lt <= 12) then
                ! Leaf litter
                if (lt /= 6 .and. lt < 10) then !deciduous
                    ! Deciduous leaf litter
                    if (C(i, 1)/C(i, 10) > 0.7) then
                        frac_loss = soil%frac_burnt(FLDEC)
                    else
                        frac_loss = WDWcons
                    end if
                    ! Calculate total weight loss
                    weightloss = min(C(i, 1)*frac_loss, C(i, 1))
                    bg_combust = bg_combust + weightloss*T_TO_KG/HEC_TO_M2
                else
                    ! Conifer needle litter
                    if (C(i, 1)/C(i, 10) > 0.7) then
                        frac_loss = soil%frac_burnt(FLCON)
                    else
                        frac_loss = WDWcons
                    end if
                    ! Calculate total weight loss
                    weightloss = min(C(i, 1)*frac_loss, C(i, 1))
                    bg_combust = bg_combust + weightloss*T_TO_KG/HEC_TO_M2
                end if
            else if (lt == IROOT) then

                ! Root litter
                frac_loss = consRoot

                ! Calculate total weight loss
                weightloss = min(C(i, 1)*frac_loss, C(i, 1))
                bg_combust = bg_combust + weightloss*T_TO_KG/HEC_TO_M2

            else if (lt == ISBL .or. lt == ILBL) then
                ! Bole litter
                if (C(i, 1)/C(i, 10) > 0.7) then
                    frac_loss = soil%frac_burnt(FLBL)
                    ! Calculate total weight loss
                    weightloss = min(C(i, 1)*frac_loss, C(i, 1))
                    agw_combust = agw_combust + weightloss*T_TO_KG/HEC_TO_M2
                    agw_prefire = agw_prefire + (C(i, 1) - weightloss)*T_TO_KG/HEC_TO_M2
                else
                    frac_loss = WDWcons
                    weightloss = min(C(i, 1)*frac_loss, C(i, 1))
                    bg_combust = bg_combust + weightloss*T_TO_KG/HEC_TO_M2
                end if
            else if (lt == ITW) then
                ! Twig litter
                if (C(i, 1)/C(i, 10) > 0.7) then
                    frac_loss = soil%frac_burnt(FLTW)
                else
                    frac_loss = WDWcons
                end if
                ! Calculate total weight loss
                weightloss = min(C(i, 1)*frac_loss, C(i, 1))
                bg_combust = bg_combust + weightloss*T_TO_KG/HEC_TO_M2
            else if (lt == ISBR) then
                ! Small branch litter
                if (C(i, 1)/C(i, 10) > 0.7) then
                    frac_loss = soil%frac_burnt(FLSBR)
                    ! Calculate total weight loss
                    weightloss = min(C(i, 1)*frac_loss, C(i, 1))
                    agw_combust = agw_combust + weightloss*T_TO_KG/HEC_TO_M2
                    agw_prefire = agw_prefire + (C(i, 1) - weightloss)*T_TO_KG/HEC_TO_M2
                else
                    frac_loss = WDWcons
                    weightloss = min(C(i, 1)*frac_loss, C(i, 1))
                    bg_combust = bg_combust + weightloss*T_TO_KG/HEC_TO_M2
                end if
            else if (lt == ILBR) then
                ! Large branch litter
                if (C(i, 1)/C(i, 10) > 0.7) then
                    frac_loss = soil%frac_burnt(FLLBR)
                    ! Calculate total weight loss
                    weightloss = min(C(i, 1)*frac_loss, C(i, 1))
                    agw_combust = agw_combust + weightloss*T_TO_KG/HEC_TO_M2
                    agw_prefire = agw_prefire + (C(i, 1) - weightloss)*T_TO_KG/HEC_TO_M2
                else
                    frac_loss = WDWcons
                    weightloss = min(C(i, 1)*frac_loss, C(i, 1))
                    bg_combust = bg_combust + weightloss*T_TO_KG/HEC_TO_M2
                end if
            else if (lt == IWDW) then
                frac_loss = WDWcons
                ! Calculate total weight loss
                weightloss = min(C(i, 1)*frac_loss, C(i, 1))
                bg_combust = bg_combust + weightloss*T_TO_KG/HEC_TO_M2
            else if (lt == IMOSS) then
                if (C(i, 1)/C(i, 10) > 0.7) then
                    frac_loss = soil%frac_burnt(FLDM)
                else
                    frac_loss = WDWcons
                end if
                weightloss = min(C(i, 1)*frac_loss, C(i, 1))
                bg_combust = bg_combust + weightloss*T_TO_KG/HEC_TO_M2
            end if

            ! Subtract weight loss
            C(i, 1) = max(C(i, 1) - weightloss, 0.0)

            ! Subtract N loss
            Nloss = C(i, 2)*frac_loss
            C(i, 2) = max(C(i, 2) - Nloss, 0.0)

            ! Update percent N
            C(i, 11) = C(i, 2)/C(i, 1)

            ! Add volatilized N to N pool
            sfan = sfan + Nloss*volN

        end do

        ! Burn live moss
        m_loss = min(soil%moss_biom*soil%frac_burnt(FLLM), soil%moss_biom)
        soil%moss_biom = max(soil%moss_biom - m_loss, 0.0)
        bg_combust = bg_combust + m_loss/plotsize
        Nloss = (soil%moss_biom*0.0046)*soil%frac_burnt(FLLM)
        sfan = sfan + Nloss*volN
        soil%M_depth = soil%moss_biom/plotsize/BULK_MOSS

        ! Remove completely removed cohorts
        ix = 0
        do i = 2, nc
            if (C(i, 1) == 0.0) then
                ix = ix + 1
            else
                do j = 1, NCOH_CHR
                    C(i - ix, j) = C(i, j)
                end do
            end if
        end do
        nc = nc - ix

        ! Calculate total weight and N content by forest floor compartment
        do i = 1, nc
            lt = int(C(i, 5))
            soil%forest_litter(lt, 1) = soil%forest_litter(lt, 1) +  C(i, 1)
            soil%forest_litter(lt, 2) = soil%forest_litter(lt, 2) +  C(i, 2)
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
        soil%dec_fuel = dec_fuel

        ! Convert to kg
        con_fuel = con_fuel/HEC_TO_M2*plotsize*T_TO_KG
        dec_fuel = dec_fuel/HEC_TO_M2*plotsize*T_TO_KG
        twig_fuel = soil%forest_litter(ITW,1)/HEC_TO_M2*plotsize*T_TO_KG
        moss_fuel = soil%forest_litter(IMOSS, 1)/HEC_TO_M2*plotsize*T_TO_KG
        duff = C(1, 1)/HEC_TO_M2*plotsize*T_TO_KG

        ! Calculate organic layer depth
        soil%O_depth = (1.0/plotsize)*((con_fuel/BULK_CON) +                   &
            (dec_fuel/BULK_DEC) + (twig_fuel/BULK_LITTER) + (duff/BULK_DUFF) + &
            (moss_fuel/BULK_MOSS))

        ! Also get rid of this year's litter cohorts
        do i = 1, FL_LEVS
            if (lt <= 12) then
                ! Leaf litter
                if (lt /= 6 .and. lt < 10) then !deciduous
                    ! Deciduous leaf litter
                    frac_loss = soil%frac_burnt(FLDEC)
                    ! Get rid of cohort
                    soil%litter(i) = max(soil%litter(i) - soil%litter(i)*      &
                        frac_loss, 0.0)
                    bg_combust = bg_combust + (soil%litter(i)*frac_loss)*      &
                        T_TO_KG/HEC_TO_M2
                else
                    ! Conifer needle litter
                    frac_loss = soil%frac_burnt(FLCON)
                    ! Get rid of cohort
                    soil%litter(i) = max(soil%litter(i) - soil%litter(i)*      &
                        frac_loss, 0.0)
                    bg_combust = bg_combust + (soil%litter(i)*frac_loss)*      &
                        T_TO_KG/HEC_TO_M2
                end if
            else if (lt == IROOT) then
                ! Root litter
                frac_loss = consRoot
                ! Get rid of cohort
                soil%litter(i) = max(soil%litter(i) - soil%litter(i)*          &
                    frac_loss, 0.0)
                bg_combust = bg_combust + (soil%litter(i)*frac_loss)*          &
                    T_TO_KG/HEC_TO_M2
            else if (lt == ISBL .or. lt == ILBL) then
                ! Bole litter
                frac_loss = soil%frac_burnt(FLBL)
                ! Get rid of cohort
                soil%litter(i) = max(soil%litter(i) - soil%litter(i)*          &
                    frac_loss, 0.0)
                agw_combust = agw_combust + (soil%litter(i)*frac_loss)*        &
                    T_TO_KG/HEC_TO_M2
                agw_prefire = agw_prefire + (soil%litter(i)*(1.0 - frac_loss))*  &
                    T_TO_KG/HEC_TO_M2
            else if (lt == ITW) then
                ! Twig litter
                frac_loss = soil%frac_burnt(FLTW)
                ! Get rid of cohort
                soil%litter(i) = max(soil%litter(i) - soil%litter(i)*          &
                    frac_loss, 0.0)
                bg_combust = bg_combust + (soil%litter(i)*frac_loss)*          &
                    T_TO_KG/HEC_TO_M2
            else if (lt == ISBR) then
                ! Small branch litter
                frac_loss = soil%frac_burnt(FLSBR)
                ! Get rid of cohort
                soil%litter(i) = max(soil%litter(i) - soil%litter(i)*          &
                    frac_loss, 0.0)
                agw_combust = agw_combust + (soil%litter(i)*frac_loss)*        &
                    T_TO_KG/HEC_TO_M2
                agw_prefire = agw_prefire + (soil%litter(i)*(1.0 - frac_loss))*  &
                    T_TO_KG/HEC_TO_M2

            else if (lt == ILBR) then
                ! Large branch litter
                frac_loss = soil%frac_burnt(FLLBR)
                ! Get rid of cohort
                soil%litter(i) = max(soil%litter(i) - soil%litter(i)*          &
                    frac_loss, 0.0)
                agw_combust = agw_combust + (soil%litter(i)*frac_loss)*        &
                    T_TO_KG/HEC_TO_M2
                agw_prefire = agw_prefire + (soil%litter(i)*(1.0 - frac_loss))*  &
                    T_TO_KG/HEC_TO_M2
            else if (lt == IMOSS) then
                frac_loss = soil%frac_burnt(FLDM)
                ! Get rid of cohort
                soil%litter(i) = max(soil%litter(i) - soil%litter(i)*          &
                    frac_loss, 0.0)
                bg_combust = bg_combust + (soil%litter(i)*frac_loss)*          &
                    T_TO_KG/HEC_TO_M2
            end if

            ! Volatilize N
            Nloss = soil%litter(i)*litter_params(i, 2)*frac_loss
            sfan = sfan + Nloss*volN

        end do

        ! Reassign attributes
        do i = 1, NCOH_MAX
            do j = 1, NCOH_CHR
                soil%cohorts(i, j) = C(i, j)
            end do
        end do
        soil%ncohort = nc

        soil%fan = sfan

        if (testing) then
             call csv_write(cons_out, siteID, .false.)
             call csv_write(cons_out, ip, .false.)
             call csv_write(cons_out, year, .false.)
             call csv_write(cons_out, soil%dmc_fire, .false.)
             call csv_write(cons_out, duff_moist, .false.)
             call csv_write(cons_out, rfs, .false.)
             call csv_write(cons_out, N_cons, .false.)
             call csv_write(cons_out, consRoot, .false.)
             call csv_write(cons_out, emis, .false.)
             call csv_write(cons_out, soil%tau_l, .false.)
             call csv_write(cons_out, duff_cons, .false.)
             call csv_write(cons_out, hum_avail, .false.)
             call csv_write(cons_out, pre_depth, .false.)
             call csv_write(cons_out, soil%O_depth, .false.)
             call csv_write(cons_out, humloss, .false.)
             call csv_write(cons_out, bg_combust, .false.)
             call csv_write(cons_out, agw_combust, .false.)
             call csv_write(cons_out, agw_prefire, .false.)
         end if

    end subroutine fuel_consumption

    !:.........................................................................:

end module Soil
