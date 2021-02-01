module Constants

!*******************************************************************************
  !
  ! This module sets up constants used throughout UVAFME
  !
!*******************************************************************************


    ! Code identifiers
    character(len = 6)  :: CODENAME = 'UVAFME' ! Model name
    character(len = 4)  :: VERSION_ID = '2018' ! Model version
    character(len = 10) :: CODE_ID             ! Model name + version

    ! Values for invalid/missing data
    real,    parameter :: RNVALID = -999.0 ! Real invalid
    integer, parameter :: INVALID = -999   ! Integer invalid

    ! Unit conversions
    real, parameter :: PI = 4.0*atan(1.0)   ! Pi (3.1415..)
    real, parameter :: DEG2RAD = PI/180.0   ! Convert from degrees to radians
    real, parameter :: M_TO_CM = 100.0      ! Convert from m to cm
    real, parameter :: CM_TO_M = 0.01       ! Convert from cm to m
    real, parameter :: MM_TO_CM = 0.1       ! Convert from mm to cm
    real, parameter :: HEC_TO_M2 = 10000.0  ! Convert from hectares to m2
    real, parameter :: M2_TO_HEC = 0.0001   ! Convert from m2 to hectares
    real, parameter :: T_TO_KG = 1000.0     ! Convert from tonnes to kg
    real, parameter :: KG_TO_T = 0.001      ! Convert from kg to tonnes
    real, parameter :: B_TO_C = 0.45        ! Convert from dry biomass to Carbon
    real, parameter :: H2O_D = 1000.0       ! Density of water

    ! Date-related constants
    integer, parameter ::  NTEMPS = 12         ! Months in a year
    integer, parameter ::  DAYS_PER_YEAR = 365 ! Days in a year

    ! Characters and files
    integer, parameter :: MAX_NLEN = 30        ! Maximum length for names
    integer, parameter :: MAX_FILE = 132       ! Maximum file name length
    integer, parameter :: MAX_LINE = 256       ! Maximum line for file
    integer, parameter :: MAX_LONG_LINE = 1000 ! Maximum length of line (headers)
    integer, parameter :: MAX_DIR = 132        ! Maximum directory name length
    integer, parameter :: MAX_CHAR = 80        ! Maximum length for messages
    integer, parameter :: MAX_PATH = 256       ! Maximum path length

    ! Array-related constants
    integer, parameter :: LIT_LEVS = 20 ! Number of litter classes
    integer, parameter :: M_TYPES = 7   ! Number of mortality types
    integer, parameter :: FC_NUM = 5    ! Number of environmental stressors
    integer, parameter :: IROOT = 13    ! Array location of roots in litter array
    integer, parameter :: ITW = 16      ! Array location of twigs in litter array
    integer, parameter :: ISBR = 17     ! Array location of small branches in litter array
    integer, parameter :: ILBR = 18     ! Array location of large branches in litter array
    integer, parameter :: ISBL = 14     ! Array location of small boles in litter array
    integer, parameter :: ILBL = 15     ! Array location of large boles in litter array
    integer, parameter :: IFIRE = 7     ! Array location of burned trees
    integer, parameter :: IWIND = 8     ! Array location of wind-killed trees

    ! Tree-related constants
    real,    parameter :: STD_HT = 1.3         ! Standard height for DBH measurements (m)
    integer, parameter :: NHC = 11             ! Number of DBH categories
    real               :: CON_LEAF_C_N = 60.0  ! Conifer leaf C:N ratio
    real               :: DEC_LEAF_C_N = 40.0  ! Deciduous leaf C:N ratio
    real               :: STEM_C_N = 450.0     ! Stem C:N ratio
    real               :: CON_LEAF_RATIO = 0.3 ! Conifer to deciduous leaf area ratio

    ! DBH bins (minimum)
    real, dimension(NHC),   parameter :: DBC_MIN = [0.5, 5.0, 10.0, 20.0,      &
        30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0]

    ! DBH bins (maximum)
    real, dimension(NHC-1), parameter :: DBC_MAX = [5.0, 10.0, 20.0, 30.0,     &
        40.0, 50.0, 60.0, 70.0, 80.0, 90.0]

    ! Other constants
    real, parameter :: PRCP_N = 0.00002     ! N deposition from rainfall (tN/cmrain)
    real, parameter :: SNOW_DENS = 217.0    ! Snow density (kg/m3) (Sturm et al. 2010)
    real, parameter :: BULK_LITTER = 44.1   ! Litter bulk density (kg/m3)
    real, parameter :: BULK_DUFF = 55.3     ! Duff bulk density (kg/m3)
    real, parameter :: BULK_MOSS = 5.2      ! Moss bulk density (kg/m3)
    real, parameter :: BULK_CON = 28.6      ! Conifer needles bulk density (kg/m3)
    real, parameter :: BULK_DEC = 100.0     ! Deciduous leaves bulk density (kg/m3)

    ! Percent of crown woody biomass in twigs (1), small branches (2),
    ! and large branches (3) -- Thonicke et al. (2010)
    real, dimension(3), parameter :: PERC_BRANCHES = [0.136, 0.223, 0.634]

    ! Parameter for effect of organic layer depth on regeneration
    ! Depends on org_tol
    real, dimension(3), parameter :: ORG_GF = [-2.4, -25.4, -50.4]

end module Constants
