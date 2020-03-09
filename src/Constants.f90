module Constants

!*******************************************************************************
  !
  !This module sets up constants used throughout UVAFME
  !
!*******************************************************************************


  !code identifiers
  character(len = 6)   :: codename = 'UVAFME'
  character(len = 4)   :: version_id = '2018'
  character(len = 10)  :: code_id

  !global constants
  real, parameter      :: pi = 4.0*atan(1.0)
  real, parameter      :: deg2rad = pi/180.0

  !characters and files
  integer, parameter   :: MAX_NLEN = 30
  integer, parameter   :: MAX_FILE = 132
  integer, parameter   :: MAX_LINE = 256
  integer, parameter   :: MAX_LONG_LINE = 1000
  integer, parameter   :: MAX_DIR = 132
  integer, parameter   :: MAX_CHAR = 80
  integer, parameter   :: MAX_FIELDS = 100

  !unit conversions
  real, parameter      :: m_to_cm = 100.0
  real, parameter      :: hec_to_m2 = 10000
  real, parameter      :: m2_to_hec = 0.0001
  real, parameter      :: mm_to_cm = 0.1

  !radiation constants

  !solar radiation model parameters
  real, parameter        :: b = 0.017214
  real, parameter        :: As = 0.409
  real, parameter        :: Ac = 0.034
  real, parameter        :: phase = -1.39
  !radiation latitude dependencies
  real, parameter        :: Amp = 38.38176
  real, parameter        :: dl_omega = 7.639437
  real, parameter        :: exrad_coef = 0.0820

  !Hargreaves evaporation constants
  real, parameter        :: H_coeff = 0.000093876
  real, parameter        :: H_addon = 17.8

  !climate-related constants
  integer, parameter     ::  NTEMPS = 12
  integer, parameter     ::  max_days_per_year = 366
  integer, parameter     ::  days_per_year = 365

  !precipitation nitrogen content
  real, parameter        :: prcp_n = 0.00002

  !global tree attributes

  !standard for DBH measurements
  real, parameter        :: std_ht = 1.3
  !number of DBH categories
  integer, parameter     :: NHC = 11
  integer, parameter     :: n_indbh = 22
  !conifer and deciduous leaf C:N ratio
  real                   :: con_leaf_c_n = 60.0
  real                   :: dec_leaf_c_n = 40.0
  !stem C:N ratio
  real                   :: stem_c_n = 450.0
  !conifer to deciduous leaf area ratio
  real                   :: con_leaf_ratio = 0.3

  !snow parameters
  real, parameter        :: snow_dens = 217.0 !Sturm et al. 2010

  !forest fuels and litter parameters
  real, parameter        :: tfc1 = 0.8
  real, parameter        :: tfc2 = 0.2
  real, parameter        :: bfc = 0.4
  real, parameter        :: bulk_l = 44.1
  real, parameter        :: bulk_duff = 76.94
  real, parameter        :: moss_bulk = 60.0
  real, parameter        :: con_bulk = 28.6
  real, parameter        :: dec_bulk = 100.0
  real, parameter        :: hum_bulk = 55.3

  !organic layer effects on regeneration
  real, dimension(3)     :: org_gf

  !data org_gf /-7.4, -22.4, -52.4/
  data org_gf /-7.4, -18.4, -32.4/


end module Constants
