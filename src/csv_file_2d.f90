!*******************************************************************************
  !
  ! Include file for csv_file.f90:
  !  Contains the body of the two-dimensional version of the
  !  writing routines
  !
  !  Record of revisions:
  !
  !       Date            Programmer                   Change
  !    ==========       ==============            =================
  !     03/26/06         Arjen Markus              Original Code
!*******************************************************************************

    ! Data dictionary: calling arguments
    integer, intent(in)  :: lun ! Unit number of csv file

    ! Data dictionary: local variables
    logical :: adv ! Advance - yes or no?
    integer :: i   ! Looping index

    ! We are advancing here
    adv = .true.

    ! Write line
    do i = 1, size(array, 2)
          call csv_write(lun, array(:, i), adv)
    enddo
