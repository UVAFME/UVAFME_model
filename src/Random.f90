module Random

!*******************************************************************************
  !
  !This module contains various random number generators.
  !
!*******************************************************************************

    use FileUtils
    implicit none

    ! Data dictionary: global variables
    integer   :: clim_seed = 0       ! Climate RNG seed
    integer   :: rng_seed            ! Regular RNG seed
    logical   :: reset_clim = .true. ! Reset the climate seed?

    private clim_seed, rng_seed, reset_clim

contains

    !:.........................................................................:

    real function urand(lb, ub, seed)
        !
        !  Returns a uniformly-distributed random number in the range (lb, ub],
        !    default to (0.0, 1.0]
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        real, optional, intent(in) :: lb   ! Upper bound
        real, optional, intent(in) :: ub   ! Lower bound
        real, optional, intent(in) :: seed ! Optional RNG seed

        ! Data dictionary: local variables
        integer :: iseed ! Local seed
        real    :: rnd   ! Random number
        real    :: lower ! Local lower bound
        real    :: upper ! Local upper bound

        ! Set the seed if it is present
        if (present(seed)) then
            iseed = int(seed)
            call set_random_seed(iseed)
        endif

        ! Set the lower bound if present, otherwise default to 0.0
        if (present(lb)) then
            lower = lb
        else
            lower = 0.0
        endif

        ! Set the upper bound if present, otherwise default to 1.0
        if (present(ub)) then
            upper = ub
        else
            upper = 1.0
        endif

        ! Get a random number with intrinsic
        call random_number(rnd)

        ! Rescale to the correct bounds
        urand = lower + (upper - lower)*rnd

    end function urand

    !:.........................................................................:

    real function nrand(mean, std, seed)
        !
        !  Returns a normally-distributed random number with given mean and
        !    standard deviation, default to mean = 0.0, sd = 1.0
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        real, optional, intent(in) :: mean ! Input mean
        real, optional, intent(in) :: std  ! Optional standard deviation
        real, optional, intent(in) :: seed ! Optional input RNG seed

        ! Data dictionary: local variables
        integer :: iseed ! Local seed
        real    :: mn    ! Local mean
        real    :: st    ! Local standard deviation
        real    :: x1    ! Random uniform between (-1.0, 1.0]
        real    :: x2    ! Random uniform between (-1.0, 1.0]
        real    :: w     ! Temporary variable
        real    :: y1    ! Temporary variable to calculate nrand
        real    :: y2    ! Temporary variable to calculate nrand

        ! Set the mean and sd if present, otherwise use defaults
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

        ! Box-Muller polar method
        do
            ! Set the seed if present, and generate 2 uniform random numbers
            ! between -1.0 and 1.0
            if (present(seed)) then
                iseed = int(seed)
                x1 = urand(-1.0, 1.0, seed)
                x2 = urand(-1.0, 1.0, seed)
            else
                x1 = urand(-1.0, 1.0)
                x2 = urand(-1.0, 1.0)
            endif

            ! Calculate sum of squared values
            ! Check to make sure not 0.0 or > 1.0
            w = x1**2 + x2**2
            if (w /= 0.0 .and. w < 1.0) exit
        end do

        ! Calculate y1 and y2 from w and x1, x2
        w = sqrt((-2.0*log(w))/w)
        y1 = x1*w
        y2 = x2*w

        ! Pick one, adjust its mean and std
        nrand = y1*st + mn

    end function nrand

    !:.........................................................................:

    real function clim_urand(lb, ub)
        !
        !  Returns a uniformly distributed random number (used for climate
        !    data) between a lower and upper bound (default to 0.0 and 1.0).
        !  Adapted from the Numerical Recipes First Edition
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictonary: constants
        ! Parameters for creating some randomness
        integer, parameter :: M1 = 259200
        integer, parameter :: M2 = 134456
        integer, parameter :: M3 = 243000
        integer, parameter :: IA1 = 7141
        integer, parameter :: IA2 = 8121
        integer, parameter :: IA3 = 4561
        integer, parameter :: IC1 = 54773
        integer, parameter :: IC2 = 28411
        integer, parameter :: IC3 = 51349

        ! Data dictionary: calling arguments
        real, optional, intent(in) :: lb ! Optional lower bound
        real, optional, intent(in) :: ub ! Optional upper bound

        ! Data dictionary: local variables
        real, dimension(97) :: r              ! Array of numbers to pull from
        real                :: lower          ! Local lower bound
        real                :: upper          ! Local upper bound
        real                :: ran1           ! Number from r
        real                :: rm1            ! Temp. variable
        real                :: rm2            ! Temp. variable
        integer             :: ix1            ! Temp. variable
        integer             :: ix2            ! Temp. variable
        integer             :: ix3            ! Temp. variable
        integer             :: j              ! Indexing variable
        logical             :: first = .true. ! First time we have called this?
        save

        ! Set the lower bound if present, otherwise default to 0.0
        if (present(lb)) then
            lower = lb
        else
            lower = 0.0
        endif

        ! Set the upper bound if present, otherwise default to 1.0
        if (present(ub)) then
            upper = ub
        else
            upper = 1.0
        endif

        ! Generate two numbers from division of large integers
        rm1 = 1.0/M1
        rm2 = 1.0/M1

        if (first .or. reset_clim) then

            ! First time we are calling this
            ! Generate some random numbers by doing some integer math
            ix1 = mod(abs(IC1 - clim_seed), M1)
            ix1 = mod(IA1*ix1 + IC1, m1)
            ix2 = mod(ix1, M2)
            ix1 = mod(IA1*ix1 + IC1, M1)
            ix3 = mod(ix1, M3)
            do j = 1, 97
                ix1 = mod(IA1*ix1 + IC1, M1)
                ix2 = mod(IA2*ix2 + IC2, M2)
                r(j) = (real(ix1) + real(ix2)*rm2)*rm1
            enddo
            first = .false.
            reset_clim = .false.
        endif

        ! Do some more math - we will keep doing just this on future calls
        ix1 = mod(IA1*ix1 + IC1, M1)
        ix2 = mod(IA2*ix2 + IC2, M2)
        ix3 = mod(IA3*ix3 + IC3, M3)

        ! Check to make sure j is in the correct range
        j = 1 + (97*ix3)/M3
        if (j > 97 .or. j < 1) then
            j = 97
            write(logf, *) 'Error in climate RNG'
        else
            ran1 = r(j)
            r(j) = (real(ix1) + real(ix2)*rm2)*rm1
        endif

        ! Rescale to correct bounds
        clim_urand = lower + (upper - lower)*ran1

    end function clim_urand

    !:.........................................................................:

    real function clim_nrand(mean, std)
        !
        !  Returns a normally distributed random number (used for climate
        !    data) with a given mean and standard deviation
        !    (default to mean = 0.0 and sd = 1.0)
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        real, optional :: mean ! Input mean
        real, optional :: std  ! Input standard deviation

        ! Data dictionary: local variables
        real :: mn ! Local mean
        real :: st ! Local standard deviation
        real :: x1 ! Random uniform between (-1.0, 1.0]
        real :: x2 ! Random uniform between (-1.0, 1.0]
        real :: w  ! Temporary variable
        real :: y1 ! Temporary variable to calculate nrand
        real :: y2 ! Temporary variable to calculate nrand

        ! Set mean and sd if present, otherwise use defaults
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

        ! Box-Muller polar method
        do
            ! Generate 2 uniform random numbers between -1.0 and 1.0
            x1 = clim_urand(-1.0, 1.0)
            x2 = clim_urand(-1.0, 1.0)

            ! Check value of w and exit if possible
            w = x1**2 + x2**2
            if (w /= 0.0 .and. w < 1.0) exit
        end do

        ! Calculate y1 and y2 from w, x1, and x2
        w = sqrt((-2.0*log(w))/w)
        y1 = x1*w
        y2 = x2*w

        ! Pick one, adjust its mean and std
        clim_nrand = y1*st + mn

    end function clim_nrand

    !:.........................................................................:

    subroutine set_site_rng_seed(fixed_seed, seed)
        !
        !  Sets the site's random number generator seed. Fixed_seed is for
        !  debugging. If true, we set the seed to be the default seed each time
        !  it is called - this allows for duplication of runs.
        !  If false - we get a random seed, or optionally use an input seed
        !  on the first call.
        !  This also means we can run the sites in any order for serial runs
        !  with fixed_seed on.
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: constants
        integer, parameter :: DEFAULT_SEED = 2345678 ! Default RNG seed

        ! Data dictionary: calling arguments
        logical,        intent(in) :: fixed_seed ! Fixed_seed true?
        real, optional, intent(in) :: seed       ! Optional RNG seed

        ! Data dictionary: local variables
        integer :: idate(8)       ! Array of outputs from date_and_time
        integer :: iseed          ! Local seed
        logical :: first = .true. ! First time we are calling this?

        save

        if (.not. fixed_seed) then

            ! We are not using a fixed seed - need to generate one on first call
            if (first) then

                ! First time we are calling - need to generate a seed

                if (.not. present(seed)) then
                    ! Create seed value
                    ! Get date and time values
                    call date_and_time(values = idate)

                    ! Get a random seed
                    call get_random_seed(iseed)

                    ! Calculate actual seed value
                    if (iseed /= 0) then
                        iseed = iseed*idate(8) ! idate(8) = milliseconds
                    else
                        ! If it's 0 somehow, use the default seed
                        iseed = DEFAULT_SEED*idate(8)
                    endif
                else
                    ! use the input seed
                    iseed = seed
                endif

                ! Set first to false, and set rng_seed to generated seed
                first = .false.
                rng_seed = iseed

                ! Set the random seed to the calculated or input seed for any
                ! future random number calls
                call set_random_seed(rng_seed)
            endif
        else
            ! fixed seed - use same seed for all sites. Seed is set to default
            ! seeds unless seed is present in call
            if (.not. present(seed)) then
                ! Use default seed
                iseed = DEFAULT_SEED
            else
                ! Use input seed
                iseed = seed
            endif

            ! Set the seed
            rng_seed = iseed

            ! Set the random number seed
            call set_random_seed(rng_seed)

        endif

    end subroutine set_site_rng_seed

    !:.........................................................................:

    subroutine get_random_seed(seed)
        !
        !  Gets an RGN seed
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        integer, intent(out) :: seed ! RNG seed

        ! Data dictionary: local variables
        integer, dimension(:), allocatable :: iseed ! Array of RNG seeds
        integer                            :: isize ! Max size of seed array


        ! Get max size of random seed
        call random_seed(size = isize)

        ! Allocate iseed to this size
        if (.not. allocated(iseed)) allocate(iseed(isize))

        ! Get a random seed
        call random_seed(get = iseed)
        seed = iseed(1)

    end subroutine get_random_seed

    !:.........................................................................:

    subroutine set_random_seed(seed)
        !
        !  Sets the random seed for random_number calls
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        integer, intent(in) :: seed ! Input RNG seed

        ! Data dictionary: local variables
        integer, dimension(:), allocatable :: iseed ! Array of RNG seeds
        integer                            :: isize ! Max size of RNG seed array

        ! Get size of random seed array and allocate
        call random_seed(size = isize)
        if (.not. allocated(iseed)) allocate(iseed(isize))

        ! Put the seed as the random seed for any future random number calls
        iseed = seed
        call random_seed(put = iseed)

    end subroutine set_random_seed

    !:.........................................................................:

    subroutine shuffle(a)
        !
        !  Shuffles an input integer array
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        integer, dimension(:), intent(inout) :: a ! Array to shuffle

        ! Data dictionary: local variables
        integer :: i       ! Looping index
        integer :: randpos ! Random location in array
        integer :: tempval ! Temporary variable for swapping
        real    :: rnd     ! Random number

        ! Loop backwards through array
        do i = size(a), 1, -1

            ! Call random number
            call random_number(rnd)

            ! Get a random position to swap with
            randpos = int(rnd*i) + 1

            ! Swap values
            tempval = a(randpos)
            a(randpos) = a(i)
            a(i) = tempval
        enddo

    end subroutine shuffle

    !:.........................................................................:

end module Random
