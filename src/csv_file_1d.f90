!*******************************************************************************
  !
  ! Include file for csv_file.f90:
  !  contains the body of the one-dimensional version of the
  !  writing routines
  !
  !  Record of revisions:
  !
  !       Date            Programmer                   Change
  !    ==========       ==============            =================
  !     03/26/06         Arjen Markus              Original Code
!*******************************************************************************

    ! Data dictionary: calling arguments
    integer, intent(in)           :: lun     ! Unit number of csv file
    logical, intent(in), optional :: advance ! Next line or no?

    ! Data dictionary: local variables
    logical :: adv ! Advance or no?
    integer :: i   ! Looping index

    ! Set advance to true as default and change if present
    adv = .true.
    if (present(advance)) adv = advance

    ! Write the line
    do i = 1, size(array) - 1
          call csv_write(lun, array(i), .false.)
    end do
    call csv_write(lun, array(size(array)), adv)
