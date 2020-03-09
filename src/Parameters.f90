module Parameters
  use Constants

  implicit none

!*******************************************************************************
  !
  !This module sets parameters used in climate change, runtime parameters, and
  !I/O routines
  !
!*******************************************************************************

!:.............................................................................:

  integer             :: clim_counter = 0, rand_counter = 0

  !basic parameters
  integer             :: numyears
  integer             :: numplots
  integer             :: maxtrees, maxcells
  integer             :: maxheight = 60.0
  real                :: growth_thresh

  !values for invalid/missing data
  real                :: rnvalid = -999.0
  integer             :: invalid = -999

  !variables determining whether to use explicit seeds for RGNs
  logical             :: fixed_seed
  logical             :: same_climate = .true.
  logical             :: debug = .false.

  !adjustments
  logical             :: adjust_altitude

  !optional files
  logical             :: use_rangelist
  logical             :: use_climstd
  logical             :: use_cexstd
  logical             :: use_gcm

  !litter parameters
  real, dimension(18, 10) :: litter_params

!:.............................................................................:

  !climate change variables
  real                :: incr_tmin_by
  real                :: incr_tmax_by
  real                :: incr_precip_by
  real                :: decr_tmin_by
  real                :: decr_tmax_by
  real                :: decr_precip_by
  real                :: tmin_change, tmax_change, precip_change
  integer             :: begin_change_year
  integer             :: start_gcm
  integer             :: end_gcm
  integer             :: gcm_duration
  integer             :: FRI
  character(len = 4)  :: incr_or_decr_prcp, incr_or_decr_temp
  integer             :: year_print_interval
  logical             :: linear_cc
  logical             :: with_clim_change
  logical             :: plot_level_data
  logical             :: tree_level_data
  logical             :: daily_soil_data
  
!:.............................................................................:

  !plot parameters
  real                :: plotsize
  real                :: rootdepth = 0.8

end module Parameters
