module FileUtils

!*******************************************************************************
  !
  ! This module manages variables related to files and their paths.
  !
!*******************************************************************************

    use Constants
    use Parameters
    use string

    implicit none

    ! Data dictionary: global constants
    character(len = 1), parameter :: SEPARATOR = '/' ! File path separater
    integer,            parameter :: BASE_UNIT = 10  ! Base unit for files the first time unit_number is called

    ! Data dictionary: global variables
    character(len = MAX_PATH) :: fqcwd        ! Current working directory
    character(len = MAX_PATH) :: inputdir     ! Input directory
    character(len = MAX_PATH) :: outputdir    ! Output directory
    character(len = MAX_PATH) :: climatedir   ! Input climate directory
    character(len = MAX_PATH) :: GCMdir       ! Input climate change directory
    character(len = MAX_PATH) :: sitedir      ! Input site directory
    character(len = MAX_PATH) :: slistdir     ! Input sitelist directory
    character(len = MAX_PATH) :: rtdir        ! Input runtime file directory
    character(len = MAX_PATH) :: splistdir    ! Input specieslist file directory
    integer                   :: rt_file      ! Unit number for input runtime file
    integer                   :: slist        ! Unit number for input sitelist file
    integer                   :: sfile        ! Unit number for input site file
    integer                   :: cfile        ! Unit number for input climate file
    integer                   :: cstdfile     ! Unit number for input climate stddev file
    integer                   :: dailyclim    ! Unit number for input daily climate file
    integer                   :: cgcmfile     ! Unit number for input climate file
    integer                   :: splist       ! Unit number for input species file
    integer                   :: rglist       ! Unit number for input rangelist file
    integer                   :: litterfile   ! Unit number for input litter parameters file
    integer                   :: cexfile      ! Unit number for input extra climate file
    integer                   :: cexstdfile   ! Unit number for input extra cliamte stddev file
    integer                   :: logf         ! Unit number for output log file
    integer                   :: sitef        ! Unit number for output site/year log file
    integer                   :: clim_unit    ! Unit number for output climate file
    integer                   :: plotvals     ! Unit number for output across-species values file
    integer                   :: biom_by_s    ! Unit number for species-specific output file
    integer                   :: biom_by_g    ! Unit number for genus-specific output file
    integer                   :: soildecomp   ! Unit number for soil output file
    integer                   :: pl_biom_by_g ! Unit number for plot-level genus-specific output file
    integer                   :: pl_biom_by_s ! Unit number for plot-level species-specific output file
    integer                   :: pl_tree      ! Unit number for tree-level species-specific output file
    integer                   :: dead_g       ! Unit number for dead genus-specific output file
    integer                   :: dead_s       ! Unit number for dead species-specific output file
    integer                   :: dead_pg      ! Unit number for plot-level dead genus-specific output file
    integer                   :: dead_ps      ! Unit number for plot-level dead species-specific output file
    integer                   :: moist_out    ! Unit number for soil moisture output file
    integer                   :: soiln_day    ! Unit number for daily soil output file

    contains

    !:.........................................................................:

    subroutine get_cwd(cwd)
        !
        !  Gets current working directory (may vary by compiler)
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *), intent(out) :: cwd ! Working directory

        ! Use intrinsic to get working directory
        call getcwd(cwd)

    end subroutine get_cwd

    !:.........................................................................:

    integer function unit_number()
        !
        !  Generates a unit number to be used in opening files.
        !  The first time the function is called, it returns BASE_UNIT (10),
        !  this should be the file list file, then increments by 1
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictonary: local variables
        integer :: iunit          ! File unit (increments after first call)
        logical :: first = .true. ! First time this has been called?
        save

        if (first) then
            ! Set first to false and iunit to base unit on first call
            iunit = BASE_UNIT
            first = .false.
        else
            ! Otherwise, increment
            iunit = iunit + 1
        endif

        ! Set to output
        unit_number = iunit

    end function unit_number

    !:.........................................................................:

    integer function open_file(filename, mode)
        !
        !  Opens the file filename if it can, returns a unit number for it.
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !    01/25/21     A. C. Foster         Added warning statements for file
        !                                       not existing or bad open
        !

        ! Data dictionary: calling arguments
        character(len = *), intent(in)           :: filename ! Name of file to open
        character(len = *), intent(in), optional :: mode     ! Optional mode ('r', 'w', 'rw')

        ! Data dictionary: local variables
        character(len = 9)           :: fmode       ! File open mode
        logical                      :: file_exists ! Does the file exist?
        character(len = MAX_PATH)    :: fname       ! Local filename (trimmed)
        character(len = MAX_CHAR)    :: message     ! Error message
        integer                      :: i           ! Looping index
        integer                      :: ios         ! I/O status
        integer                      :: iunit       ! File unit number
        integer, dimension(MAX_PATH) :: farray      ! Array of characters of file name

        ! Get mode of open (read, write, or read/write)
        ! Defaults to read/write
        if (present(mode)) then
            if (mode == 'r' .or. mode == 'R') then
                fmode = 'read'
            else if (mode == 'w' .or. mode == 'W') then
                fmode = 'write'
            else if (mode == 'rw' .or. mode == 'RW'                            &
                       .or. mode == 'wr' .or. mode == 'WR') then
                fmode = 'readwrite'
            else
                fmode = 'readwrite'
            endif
        else
            fmode = 'readwrite'
        endif

        ! Trim filename of whitespace
        fname = trim(adjustl(filename))

        if (fmode == 'read' .or. fmode == 'readwrite') then

            ! Check for valid name of file - only ASCII allowed
            farray = 0
            do i = 1, len_trim(fname)
                farray(i) = ichar(fname(i:i))
            enddo
            if (any(farray > 127)) then
                call fatal_error("INVALID filename")
            endif
        endif

        ! Does the file exist?
        inquire(file = fname, exist = file_exists)

        ! Open file if conditions are correct
        if (file_exists .and. fmode == 'write') then
            write(message, '(A,A,A)') "File ", fname(1:len_trim(fname)),       &
                " exists. Cannot open write only."
            call warning(message)
            open_file = INVALID
        else if (.not. file_exists .and. fmode == 'read') then
            write(message, '(A,A,A)') "File ", fname(1:len_trim(fname)),       &
                " does not exist. Can't read."
            call warning(message)
            iunit = INVALID
        else
            iunit = unit_number()
            open(iunit, file = fname, action = fmode, iostat = ios)
            if (ios .ne. 0) then
                write(message, '(A,A,A,I6)') "Problem opening",                &
                    fname(1:len_trim(fname)), " ios: ", ios
                call warning(message)
                iunit = INVALID
            endif
        endif

        open_file = iunit

    end function open_file

    !:.........................................................................:

    integer function count_records(funit, nheaders)
        !
        !  Counts the number of rows in a file, excluding first nheaders rows
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !    01/25/21     A. C. Foster         Added unit number to warning
        !                                        message
        !

        ! Data dictionary: calling arguments
        integer, intent(in) :: funit    ! Unit number of file
        integer, intent(in) :: nheaders ! Number of lines to skip

        ! Data dictionary: local variables
        integer                   :: ios     ! I/O status
        integer                   :: lpcount ! Line counter
        character(len = MAX_LINE) :: line    ! Line read in
        character(len = MAX_CHAR) :: message ! Error message

        lpcount = 0

        ! Read each line in file until read end, incrementing line counter
         do
            read(funit, *, end = 10, iostat = ios) line
            if (ios == 0) then
                ! Success, increment counter
                lpcount = lpcount + 1
             else
                 write(message, '(A,I4,A,I6)') "Problem reading file unit",    &
                    funit, " ios: ", ios
                call fatal_error(message)
             endif
         end do

    10  continue

        ! Return number of total lines minus headers
        count_records = lpcount - nheaders

        ! Rewind the file
        rewind(funit)

    end function count_records

    !:.........................................................................:

    subroutine build_pathname(subdir, filename, pathname)
        !
        !  Returns the full path to a file
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *), intent(in)  ::  subdir   ! Subdirectory
        character(len = *), intent(in)  ::  filename ! File name
        character(len = *), intent(out) ::  pathname ! Full path

        ! Data dictionary: local variables
        character(len = MAX_PATH) :: pname   ! Pathname
        character(len = MAX_FILE) :: tmpsub  ! Trimmed subdirectory name
        character(len = MAX_FILE) :: tmpfile ! Trimmed filename
        character(len = 15)       :: fmtstr  ! Format string
        integer                   :: flen    ! File name length
        integer                   :: slen    ! Subdirectory name length

        ! Reset pathname
        pathname = ''

        ! Create temporary directory and filename and remove whitespace
        tmpsub = subdir
        tmpfile = filename
        call remove(tmpsub)
        call remove(tmpfile)

        ! Get length of directory and filename
        slen = len_trim(adjustl(tmpsub))
        flen = len_trim(adjustl(tmpfile))

        ! Construct format string
        write(fmtstr, '(a,i2,a,a,i2,a)'), '(a', slen, ',a1,' ,'a', flen, ')'

        ! Construct path
        write(pname, fmtstr) tmpsub(1:slen), SEPARATOR, tmpfile(1:flen)

        ! Remove any whitespace
        call remove(pname)

        pathname = pname

    end subroutine build_pathname

    !:.........................................................................:

    subroutine build_filename(str, suffix, filename)
        !
        !  Returns a filename from a specified string and a suffix
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *), intent(in)  :: suffix   ! Suffix
        character(len = *), intent(in)  :: str      ! Input string
        character(len = *), intent(out) :: filename ! Output filename

        ! Data dictionary: local variables
        character(len = MAX_FILE/2) :: tmpstr ! Trimmed string
        character(len = MAX_FILE/2) :: tmpsuf ! Trimmed suffix
        character(len = MAX_FILE)   :: fname  ! Output filename
        character(len = 15)         :: fmtstr ! Format string
        integer                     :: strlen ! Length of string
        integer                     :: suflen ! Length of suffix

        ! Reset filename
        filename = ''

        ! Create temporary string and suffix and remove whitespace
        tmpstr = str
        tmpsuf = suffix
        call remove(tmpstr)
        call remove(tmpsuf)

        ! Get length of string and suffix
        strlen = len_trim(adjustl(tmpstr))
        suflen = len_trim(adjustl(tmpsuf))

        ! Construct format string
        write(fmtstr, '(a,i2,a,i2,a)'), '(a', strlen, ',a', suflen, ')'

        ! Construct filename and remove whitespace
        write(fname, fmtstr), tmpstr(1:strlen), tmpsuf(1:suflen)
        call remove(fname)

        filename = fname

    end subroutine build_filename

    !:.........................................................................:

    subroutine warning(message)
        !
        !  Prints a warning message to the screen but does not stop the program
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *), intent(in) :: message ! Error message

        ! Print the message
        write(*, *) message

    end subroutine warning

    !:.........................................................................:

    subroutine fatal_error(message)
        !
        !  Stops the program and prints an error message to the screen
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *), intent(in) :: message ! Error message

        ! Write error message
        write(*, *) message

        ! Stop the program
        stop 10

    end subroutine fatal_error

    !:.........................................................................:

end module FileUtils
