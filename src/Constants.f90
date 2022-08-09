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
    real, parameter :: KM_TO_M = 1000.0     ! Convert from km to m
    real, parameter :: MIN_TO_SEC = 60.0    ! Convert from minutes to seconds
    real, parameter :: HR_TO_MIN = 60.0     ! Convert from hours to minutes

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
    integer, parameter :: M_TYPES = 9   ! Number of mortality types
    integer, parameter :: FC_NUM = 6    ! Number of environmental stressors
    integer, parameter :: FL_LEVS = 10  ! Number of fuel categories
    integer, parameter :: FLDEC = 1     ! Array location of deciduous leaf litter in fuel array
    integer, parameter :: FLCON = 2     ! Array location of conifer needle litter in fuel array
    integer, parameter :: FLTW = 3      ! Array location of twig litter in fuel array
    integer, parameter :: FLSBR = 4     ! Array location of small branch litter in fuel array
    integer, parameter :: FLLBR = 5     ! Array location of large branch litter in fuel array
    integer, parameter :: FLBL = 6      ! Array location of bole litter in fuel array
    integer, parameter :: FLSH = 10     ! Aray location of live shrub fuel in fuel array
    integer, parameter :: FLRT = 9      ! Array location of root litter in fuel array
    integer, parameter :: FLDM = 8      ! Array location of moss litter in fuel array
    integer, parameter :: FLLM = 7      ! Array location of live moss in fuel array
    integer, parameter :: IROOT = 13    ! Array location of roots in litter array
    integer, parameter :: ITW = 16      ! Array location of twigs in litter array
    integer, parameter :: ISBR = 17     ! Array location of small branches in litter array
    integer, parameter :: ILBR = 18     ! Array location of large branches in litter array
    integer, parameter :: ISBL = 14     ! Array location of small boles in litter array
    integer, parameter :: ILBL = 15     ! Array location of large boles in litter array
    integer, parameter :: IFIRE = 7     ! Array location of burned trees
    integer, parameter :: IHARV = 8     ! Array location of harvested trees
    integer, parameter :: IWIND = 9     ! Array location of wind-killed trees

    ! Tree-related constants
    real,    parameter :: STD_HT = 1.3         ! Standard height for DBH measurements (m)
    integer, parameter :: NHC = 11             ! Number of DBH categories
    real               :: CON_LEAF_C_N = 60.0  ! Conifer leaf C:N ratio
    real               :: DEC_LEAF_C_N = 40.0  ! Deciduous leaf C:N ratio
    real               :: STEM_C_N = 450.0     ! Stem C:N ratio
    real               :: CON_LEAF_RATIO = 0.3 ! Conifer to deciduous leaf area ratio

    ! DBH bins (minimum)
    real, dimension(NHC),   parameter :: DBC_MIN = [0.1, 5.0, 10.0, 20.0,      &
        30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0]

    ! DBH bins (maximum)
    real, dimension(NHC-1), parameter :: DBC_MAX = [5.0, 10.0, 20.0, 30.0,     &
        40.0, 50.0, 60.0, 70.0, 80.0, 90.0]

    ! Other constants
    real, parameter :: PRCP_N = 0.00002     ! N deposition from rainfall (tN/cmrain)
    real, parameter :: SNOW_DENS = 217.0    ! Snow density (kg/m3) (Sturm et al. 2010)
    real, parameter :: BULK_LITTER = 44.1   ! Litter bulk density (kg/m3)
    real, parameter :: BULK_DUFF = 55.3     ! Duff bulk density (kg/m3)
    real, parameter :: BULK_MOSS = 18       ! Moss bulk density (kg/m3)
    real, parameter :: BULK_CON = 28.6      ! Conifer needles bulk density (kg/m3)
    real, parameter :: BULK_DEC = 100.0     ! Deciduous leaves bulk density (kg/m3)
    real, parameter :: DRYING_RAT = 66000.0 ! Drying ratio for fuel (Thonicke et al. 2010)

    ! Percent of crown woody biomass in twigs (1), small branches (2),
    ! and large branches (3) -- Thonicke et al. (2010)
    real, dimension(3), parameter :: PERC_BRANCHES = [0.136, 0.223, 0.634]

    ! Parameter for effect of organic layer depth on regeneration
    ! Depends on org_tol
    real, dimension(3), parameter :: ORG_GF = [-0.1, -10.4, -52.4]

end module Constants
