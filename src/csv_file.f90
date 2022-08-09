module csv_file
!*******************************************************************************
  !
  ! This module was written to facilitate writing csv files
  !
  ! For convenience, the generic name "csv_write" can be used
  !  instead of the individual routines.
  !
  ! The file to write to must already be opened as a LU-number
  !   is passed.
  !
  ! Layout of the CSV-file:
  !   - single items are written to the end of the current record
  !   - one-dimensional items are also written to the end of the current record
  !   - two-dimensional items are written to separate records, one for each row
  !   - except for the two-dimensional versions, all routines allow you to
  !       suppress advancing to the next record:
  !   - for single items you must indicate whether to advance or not
  !   - for one-dimensional items, the argument is optional. Default is to
  !       advance.
  !
  ! Note on the format:
  !  CSV-files apparently come in different guises (Kernighan and Pike, The
  !  practice of Programming, Addison-Wesley, 1999).
  !  This module uses the following rules:
  !     - items are always separated by a single comma (,)
  !     - string items are delimited by double quotes (")
  !     - embedded double quotes are treated by doubling the quote
  !     - trailing blanks are considered irrelevant
  !
!*******************************************************************************

    implicit none

    ! Interface for procedures
    interface csv_write
        module procedure csv_write_integer
        module procedure csv_write_integer_1d
        module procedure csv_write_integer_2d
        module procedure csv_write_char
        module procedure csv_write_char_1d
        module procedure csv_write_char_2d
        module procedure csv_write_real
        module procedure csv_write_real_1d
        module procedure csv_write_real_2d
        module procedure csv_write_dble
        module procedure csv_write_dble_1d
        module procedure csv_write_dble_2d
    end interface

contains

    !:.........................................................................:

    subroutine csv_next_record(lun)
        !
        !  Helper subroutine to go to the next record in a file.
        !  Note: it may result in a superfluous comma at the end of the previous
        !  record.
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    03/26/06    Arjen Markus         Original Code
        !

        ! Data dictionary - calling arguments
        integer, intent(in) :: lun ! Unit number of csv file

        ! Write the end, the current record is closed so the next write will be
        !   to the new record
        write(lun, '(a)') ''

    end subroutine csv_next_record

    !:.........................................................................:

    subroutine csv_write_integer(lun, value, advance)
        !
        !  Writes a single integer to the csv file
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    03/26/06    Arjen Markus         Original Code
        !

        ! Data dictionary - calling arguments
        integer, intent(in) :: lun     ! Unit number of csv file
        integer, intent(in) :: value   ! Value to write
        logical, intent(in) :: advance ! Advance (T or F), so that more items can be written to same record

        ! Data dictionary - local variables
        character(len = 40) :: buffer ! Character of value so we can write

        ! Write the value to a character variable first and trim
        write(buffer, '(I10)') value
        buffer = adjustl(buffer)

        ! If we are advancing, just write
        if (advance) then
            write(lun,'(a)') trim(buffer)
        else
            ! Most likely: write the comma only when needed
            ! - depends on other actions
            write(lun, '(a,a)', advance = 'no') trim(buffer), ','
        endif

    end subroutine csv_write_integer

    !:.........................................................................:

    subroutine csv_write_real(lun, value, advance)
        !
        !  Writes a single real to the csv file
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    03/26/06    Arjen Markus         Original Code
        !

        ! Data dictionary - calling arguments
        integer, intent(in) :: lun     ! Unit number of file
        real, intent(in)    :: value   ! Value to write
        logical, intent(in) :: advance ! Advance (T or F), so that more items can be written to same record

        ! Data dictionary - local variables
        character(len = 40) :: buffer ! Character of value so we can write

        ! Write the value to a character variable first and trim
        write(buffer, '(G14.6)') value
        buffer = adjustl(buffer)

         ! If we are advancing, just write
        if (advance) then
            write(lun, '(a)') trim(buffer)
        else
           ! Most likely: write the comma only when needed
           ! - depends on other actions
           write(lun, '(a,a)', advance = 'no') trim(buffer), ','
        endif

    end subroutine csv_write_real

    !:.........................................................................:

    subroutine csv_write_dble(lun, value, advance)
        !
        !  Writes a single double-precision real to the csv file
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    03/26/06    Arjen Markus         Original Code
        !

        ! Data dictionary - calling arguments
        integer, intent(in)                  :: lun     ! Unit number of file
        real(kind = kind(1.0d0)), intent(in) :: value   ! Value to write
        logical, intent(in)                  :: advance ! Advance (T or F), so that more items can be written to same record

        ! Data dictionary - local variables
        character(len = 40) :: buffer ! Character of value so we can write

        ! Write the value to a character value so we can write
        write(buffer, '(G20.12)') value
        buffer = adjustl(buffer)

        ! If we are advancing, just write
        if (advance) then
            write(lun, '(a)') trim(buffer)
        else
            ! Most likely: write the comma only when needed
            ! - depends on other actions
            write(lun, '(a,a)', advance = 'no') trim(buffer), ','
        endif

    end subroutine csv_write_dble

    !:.........................................................................:

    subroutine csv_write_char(lun, value, advance)
        !
        !  Writes a single character string to the csv file
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    03/26/06    Arjen Markus         Original Code
        !

        ! Data dictionary - calling arguments
        integer, intent(in)            :: lun     ! Unit number of file
        character(len = *), intent(in) :: value   ! Value to write
        logical, intent(in)            :: advance ! Advance (T or F), so that more items can be written to same record

        ! Data dictionary - local variables
        integer                       :: k      ! Looping index
        integer                       :: pos    ! Looping index
        integer                       :: posb   ! Looping index
        character(len = 2*len(value)) :: buffer ! Character of value so we can write

        ! Set buffer to input value
        buffer = value

        ! Check for nasty characters (")
        k = index(value, '"')
        pos = 1
        posb = 1

        ! Get rid of quotes
        do while (k >= 1)
            buffer(posb:) = value(pos:)
            buffer(posb+k:) = '"' // value(pos+k:)
            pos = pos + k + 1
            posb = posb + k + 2
            k = index(value(pos:), '"')
        enddo

        ! If we are advancing, just write
        if (advance) then
            write(lun,'(a)') trim(buffer)
        else
            ! Most likely: write the comma only when needed
            ! - depends on other actions
            write(lun, '(a,a)', advance = 'no') trim(buffer), ','
        endif

    end subroutine csv_write_char

    !:.........................................................................:
    ! csv_write_integer/real/dble_1d --
    !     Write a one-dimensional array of items to the CSV-file
    ! Inputs:
    !    lun:      LU-number of the CSV-file
    !    array:    array to write
    !    advance:  advance (.true.) or not, so that more items can be written
    !               to the same record
    ! Result:
    !    The array is written to the current record of the CSV-file
    ! Note:
    !    Because the four routines of this type differ only in the data type,
    !     we use an include file for the body.
    !
    !:.........................................................................:

    subroutine csv_write_integer_1d(lun, array, advance)
        integer, dimension(:), intent(in)  :: array

        include 'csv_file_1d.f90'

    end subroutine csv_write_integer_1d

    !:.........................................................................:

    subroutine csv_write_real_1d(lun, array, advance)
        real, dimension(:), intent(in)  :: array

        include 'csv_file_1d.f90'

    end subroutine csv_write_real_1d

    !:.........................................................................:

    subroutine csv_write_dble_1d(lun, array, advance)
        real(kind = kind(1.0d0)), dimension(:), intent(in)  :: array

        include 'csv_file_1d.f90'

    end subroutine csv_write_dble_1d

    !:.........................................................................:

    subroutine csv_write_char_1d(lun, array, advance)
        character(len = *), dimension(:), intent(in)  :: array

        include 'csv_file_1d.f90'

    end subroutine csv_write_char_1d

    !:.........................................................................:
    ! csv_write_integer/real/dble_2d --
    !     Write a two-dimensional array of items to the CSV-file
    ! Inputs:
    !    lun:    LU-number of the CSV-file
    !    array:  Array to write
    ! Result:
    !    The array is written to the current record of the CSV-file
    ! Note:
    !    Because the four routines of this type differ only in the data type,
    !     we use an include file for the body.
    !:.........................................................................:

    subroutine csv_write_integer_2d(lun, array)
        integer, dimension(:,:), intent(in)  :: array

        include 'csv_file_2d.f90'

    end subroutine csv_write_integer_2d

    !:.........................................................................:

    subroutine csv_write_real_2d(lun, array)
        real, dimension(:,:), intent(in)  :: array

        include 'csv_file_2d.f90'

    end subroutine csv_write_real_2d

    !:.........................................................................:

    subroutine csv_write_dble_2d(lun, array)
        real(kind = kind(1.0d0)), dimension(:,:), intent(in)  :: array

        include 'csv_file_2d.f90'

    end subroutine csv_write_dble_2d

    !:.........................................................................:

    subroutine csv_write_char_2d(lun, array)
        character(len=*), dimension(:,:), intent(in)  :: array

        include 'csv_file_2d.f90'

    end subroutine csv_write_char_2d

    !:.........................................................................:

end module csv_file
