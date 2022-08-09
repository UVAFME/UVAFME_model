module Utilities

!***********************************************************************
  !
  ! This  moudule contains helper functions/subroutines used throughout the
  ! program
  !
!***********************************************************************

  use Constants
  implicit none

  ! Interface for swap procedure
  interface swap
      module procedure swap_int
      module procedure swap_real
  end interface

contains

    !:.........................................................................:

    real function kron(x)
        !
        !  Calculates the Kronecker delta
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        real, intent(in) :: x ! Input number

        if (x > 0.0) then
            kron = 1.0
        else
            kron = 0.0
        endif

    end function kron

    !:.........................................................................:

    real function roundtoN(x, N)
        !
        !  Rounds to N decimal places
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        real,    intent(in) :: x ! Input number
        integer, intent(in) :: N ! Decimal places to round to

        ! Data dictionary: local variables
        real :: factor ! Factor for rounding

        ! Get factor
        factor = 10.0**N

        ! Round to N
        roundtoN = float(int(x* factor + 0.5))/factor

    end function roundtoN

    !:.........................................................................:

    subroutine stddev(array, arith_mean, std_deviation, missing)
        !
        !  Computs mean and standard deviation of an array, ignoring values
        !   equal to missing
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        real, dimension(:), intent(in)  :: array         ! Input array
        real,               intent(out) :: arith_mean    ! Arithmetic mean
        real,               intent(out) :: std_deviation ! Standard deviation
        real, optional,     intent(in)  :: missing       ! Flag for missing values

        ! Data dictionary: local variables
        real    :: sqr_sum ! Squared sum
        integer :: nelems  ! Number of elements in array

        if (present(missing)) then

            ! Count how many items in the array aren't equal to missing
            nelems = count(array /= missing)

            ! Get mean, masking missing values
            arith_mean = sum(array, mask = array /= missing)/float(nelems)

            ! Get standard deviation, maasking missing values
            sqr_sum = sum(array**2, mask = array /= missing)
        else

            ! Use whol array
            nelems = size(array)
            arith_mean = sum(array)/float(nelems)
            sqr_sum = sum(array**2)
        endif

        ! Check for 0.0
        if (nelems == 0) then
            std_deviation = 0.0
            return
        endif

        ! Calculate sample standard deviation
        std_deviation = sqrt((sqr_sum - nelems*arith_mean**2)/(nelems - 1))

    end subroutine stddev

    !:.........................................................................:

    subroutine sort(array, length)
        !
        !  Sorts a character arry in alphabetical order.
        !  Sorts using Shell sort with lexigraphical sorting
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *), dimension(length), intent(inout) :: array  ! Input array
        integer,                               intent(in)    :: length ! Length of array

        ! Data dictionary: local variables
        character(len = MAX_NLEN) :: temp  ! Temporary character
        integer                   :: inc   ! Gap size
        integer                   :: lower ! Lower bound of array
        integer                   :: upper ! Upper bound of array
        integer                   :: i, j  ! Looping indices

        lower = lbound(array, 1)
        upper = ubound(array, 1)

        inc = length/2

        ! We continue until the stride is the lower bound
        do while (inc > 0)
            do i = inc + lower, upper
                temp = adjustl(array(i))
                j = i
                ! Cannot use j>inc .and. array(j-inc)>temp if it doesn't
                ! short-circuit
                do while (j > inc)
                    if (array(j - inc) > temp) then
                        array(j) = array(j - inc)
                        j = j - inc
                    else
                        exit
                    endif
                enddo
                array(j) = temp
            enddo

            if (inc == 2) then
                inc = 1
            else
                inc = int(real(inc)/2.2)
            endif
        enddo

    end subroutine sort

    !:.........................................................................:

    subroutine merge_index(A, Ai)
        !
        !  Sorts input real array A and also returns array with un-sorted array
        !   locations of sorted array.
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/01/20     A. C. Foster        Original Code
        !

        ! Data dictionary: calling arguments
        real,    dimension(:), intent(inout) :: A  ! Input array
        integer, dimension(:), intent(out)   :: Ai ! Array of sorted indices

        ! Data dictionary: local variables
        real,    dimension((size(A) + 1)/2) :: work  ! Working array
        integer, dimension((size(A) + 1)/2) :: worki ! Working array - indices
        integer                             :: n     ! Size of A
        integer                             :: i     ! Looping index

        ! Create the index array
        n = size(A)
        Ai = (/(i, i = 1,n)/)

        ! Sort the array
        call MergeSort(A, Ai, work, worki)

    end subroutine merge_index

    recursive subroutine MergeSort(A, Ai, work, worki)
        !
        !  Sorts input array A and also returns array with un-sorted array
        !   locations of sorted array.
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/01/20     A. C. Foster        Original Code
        !    01/25/21     A. C. Foster        Updated to a use a recursive
        !                                        subroutine
        !

        implicit none

        ! Data dictionary: calling arguments
        real,    dimension(:), intent(inout) :: A     ! Array to sort
        real,    dimension(:), intent(inout) :: work  ! Working array
        integer, dimension(:), intent(inout) :: Ai    ! Array indices
        integer, dimension(:), intent(inout) :: worki ! Working array - indices

        ! Data dictionary: local variables
        integer :: half ! Half the size of the array

        ! Get half the size of A
        half = (size(A) + 1)/2

        if (size(A) < 2) then
            ! We are done
            continue
        else if (size(A) == 2) then

            ! Compare sizes
            if (A(1) > A(2)) then
                ! Swap so that they are sorted
                call swap(A(1), A(2))
                call swap(Ai(1), Ai(2))
            end if
        else
            ! Still too big - split in two
            call MergeSort(A(:half), Ai(:half), work, worki)
            call MergeSort(A(half+1:), Ai(half+1:), work, worki)
            if (A(half) > A(half + 1)) then
                work(1:half) = A(1:half)
                worki(1:half) = Ai(1:half)
                call merge(work(1:half), A(half+1:), A, worki(1:half),         &
                    Ai(half+1:), Ai)
            end if
        end if
    end subroutine MergeSort

    !:.........................................................................:

    subroutine swap_real(x, y)
        !
        !  Helper subroutine to swap x and y
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/25/21     A. C. Foster        Original Code
        !

        implicit none

        ! Data dictionary: calling arguments
        real, intent(inout) :: x, y ! Values to swap

        ! Data dictionary: local variables
        real :: tmp ! Temporary variable

        tmp = x
        x = y
        y = tmp

    end subroutine swap_real

    !:.........................................................................:

    subroutine swap_int(x, y)
        !
        !  Helper subroutine to swap x and y
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/25/21     A. C. Foster        Original Code
        !

        implicit none

        ! Data dictionary: calling arguments
        integer, intent(inout) :: x, y ! Values to swap

        ! Data dictionary: local variables
        integer :: tmp ! Temporary variable

        tmp = x
        x = y
        y = tmp

    end subroutine swap_int

    !:.........................................................................:

    subroutine merge(A, B, C, Ai, Bi, Ci)
        !
        !  Merges arrays A and B, along with their accompanying array index
        !   locations
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/01/20     A. C. Foster        Original Code
        !    01/25/21     A. C. Foster        Updated to a use a recursive
        !                                       call

        implicit none

        ! Data dictionary: calling arguments
        real,    dimension(:), intent(in)    :: A   ! 1st array to merge
        real,    dimension(:), intent(in)    :: B   ! 2nd array to merge
        integer, dimension(:), intent(in)    :: Ai  ! Indices of array A
        integer, dimension(:), intent(in)    :: Bi  ! Indices of array B
        real,    dimension(:), intent(inout) :: C   ! Merged array
        integer, dimension(:), intent(inout) :: Ci  ! Merged index array

        ! Data dictionary: local variables
        integer :: i ! Looping index (A)
        integer :: j ! Looping index (B)
        integer :: k ! Looping index (C)

        if (size(A) + size(B) > size(C)) stop

        i = 1
        j = 1
        do k = 1, size(C)
            if (i <= size(A) .and. j <= size(B)) then
                ! We are still looking through arrays

                if (A(i) <= B(j)) then
                    ! Next value is from array A
                    Ci(k) = Ai(i)
                    C(k) = A(i)
                    i = i + 1
                else
                    ! Next value is from array B
                    Ci(k) = Bi(j)
                    C(k) = B(j)
                    j = j + 1
                end if

            else if (i <= size(A)) then
                ! We've stepped through all of B
                Ci(k) = Ai(i)
                C(k) = A(i)
                i = i + 1
            else if (j <= size(B)) then
                ! We've stepped through all of A
                Ci(k) = Bi(j)
                C(k) = B(j)
                j = j + 1
            end if
        end do
    end subroutine merge

end module Utilities
