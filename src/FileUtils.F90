module FileUtils
  use Constants
  use Parameters
  use string

  implicit none

!*******************************************************************************
  !
  ! This module manages variables related to files and their paths.
  !
!*******************************************************************************

  !set global parameters related to operating system
#ifdef PC
  integer, parameter   :: MAX_PATH = 260
  character(len = 1)   :: separator = '\'
#else
  integer, parameter   :: MAX_PATH = 256
  character(len = 1)   :: separator = '/'
#endif

  !names of directories for input/output
  character(len = MAX_PATH) :: fqcwd
  character(len = MAX_PATH) :: inputdir, outputdir, climatedir, GCMdir
  character(len = MAX_PATH) :: sitedir, slistdir

  integer, parameter        :: base_unit = 10

  !names of unit numbers so we can more easily tell which file is which
  !input
  integer                 :: rt_file = 0, slist = 0, sfile = 0
  integer                 :: cfile = 0, cstdfile = 0
  integer                 :: cgcmfile = 0, splist = 0, rglist = 0
  integer                 :: altlist = 0, harvfile = 0, litterfile = 0
  integer                 :: initdbh = 0, initspec = 0, initgliht = 0
  integer                 :: cexfile = 0, cexstdfile = 0

  !output
  integer                 :: logf, sitef
  integer                 :: c_and_n, clim_unit, tld, plotvals
  integer                 :: biom_by_s, biom_by_g, soildecomp
  integer                 :: pl_biom_by_g, pl_biom_by_s
  integer                 :: pl_tree, dead_g, dead_s, dead_pg, dead_ps
  integer                 :: perm_day, soil_day, carb_day, clim_day
  integer                 :: moss_unit, soiln_day, regen, fuels
  integer                 :: testRand, testCC


contains

	!:.........................................................................:

  subroutine get_cwd(fqcwd)
		!gets current working directory, may vary by compiler
		!Outputs:
		!	fqcwd:  working directory

    character(len = *), intent(out) :: fqcwd

		call getcwd(fqcwd)

	end subroutine get_cwd

	!:.........................................................................:

	function unit_number()
		!generates a unit number, to be used in opening files. 
		!The first time the function is called, it returns base_unit (10), this 
		 !should be file_list.txt file, then increments by 1
		!Outputs:
		!	unit_number:  unit number of file

		integer          :: unit_number
		integer          :: base_unit = 10
		integer          :: iunit
		logical          :: first = .true.
		save

		if (first) then
			iunit = base_unit
			first = .false.
		else
			iunit = iunit + 1
		endif

		unit_number = iunit

	end function unit_number

	!:.........................................................................:

	function open_file(filename, mode)
		!opens the file filename if it can, returns a unit number for it.
		!Author: Katherine Holcomb, 2012 v. 1.0
		!Inputs:
		!   filename: name of file
		!   mode:     optional mode (r/w/rw)
		!Results:
		!   the file is opened

		integer                                   :: open_file
		character(len = *), intent(in)            :: filename
		character(len = *), intent(in), optional  :: mode
		character(len = 9)                        :: fmode
		logical                                   :: file_exists
		character(len = MAX_PATH)                 :: fname
		character(len = MAX_CHAR)                 :: message
		integer                                   :: i, ios, iunit
		integer, dimension(MAX_PATH)              :: farray

		!get mode of open (read, write, or read/write)
		!defaults to read/write
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

		!trim filename of whitespace
		fname = trim(adjustl(filename))

		if (fmode == 'read' .or. fmode == 'readwrite') then
			!file needs to exist
			farray = 0
			do i = 1, len_trim(fname)
				farray(i) = ichar(fname(i:i))
			enddo

			if (any(farray > 127)) then
				call fatal_error("Invalid filename")
			endif

		endif

		!does the file exist?
		inquire(file = fname, exist = file_exists)

		!open file if conditions are correct
		if (file_exists .and. fmode == 'write') then
			write(message, '(a,a,a)') "File ", fname(1:len_trim(fname)),       &
				" exists. Cannot open write only."
			call warning(message)
			open_file = invalid
		else if (.not. file_exists .and. fmode == 'read') then
			iunit = invalid
		else
			iunit = unit_number()
			open(iunit, file = fname, action = fmode, iostat = ios)
			if (ios .ne. 0) then
				iunit = invalid
			endif
		endif

		open_file = iunit

	end function open_file

	!:.........................................................................:

	function count_records(funit, nheaders)
		!counts the number of rows in a file not including the first nheaders 
		  !rows.
		!Author: Katherine Holcomb, 2012 v. 1.0
		!Inputs:
		!   funit:     unit number of file
		!   nheaders:  number of headers
		!Outputs:
		!	count_records: the number of rows in a file
		integer                    :: count_records
		integer, intent(in)        :: funit
		integer, intent(in)        :: nheaders
		integer                    :: ios, lpcount
		character(len = MAX_LINE)  :: line

		lpcount = 0

		!read each line in file until read end, incrementing line counter
		 do
			read(funit, *, end = 10, iostat = ios) line
			if (ios == 0) then
				!success, increment counter
				lpcount = lpcount + 1
			 else
				call fatal_error("Unable to read file")
			 endif
		 end do

	10  continue

		count_records = lpcount - nheaders

		rewind(funit)

	end function count_records

	!:.........................................................................:

	subroutine build_pathname(subdir, filename, pathname)
		!returns the full path to a file
		!Author: Katherine Holcomb, 2012 v. 1.0
		!Inputs:
		!   subdir:    subdirectory
		!   filename:  name of file
		!Outputs:
		!	pathname: full path to file

		character(len = *), intent(in)   ::  subdir
		character(len = *), intent(in)   ::  filename
		character(len = *), intent(out)  ::  pathname

		character(len = MAX_PATH)        ::  pname
		character(len = MAX_FILE)        ::  tmpsub, tmpfile
		character(len = 15)              ::  fmtstr
		integer                          ::  flen, slen

		!reset pathname
		pathname = ''

		!create temporary directory and filename and remove whitespace
		tmpsub = subdir
		tmpfile = filename
		call remove(tmpsub)
		call remove(tmpfile)

        !get length of directory and filename
		slen = len_trim(adjustl(tmpsub))
		flen = len_trim(adjustl(tmpfile))

        !construct path
		write(fmtstr, '(a,i2,a,a,i2,a)'), '(a', slen, ',a1,' ,'a', flen, ')'
		write(pname, fmtstr) tmpsub(1:slen), separator, tmpfile(1:flen)
		call remove(pname)

		pathname = pname

	end subroutine build_pathname

	!:.........................................................................:

	subroutine build_filename(string, suffix, filename)
		!returns a filename from a specified string and a suffix
		!Author: Katherine Holcomb, 2012 v. 1.0
		!Inputs:
		!   string:   subdirectory
		!   suffix:   name of file
		!Outputs:
		!	filename: filename

		character(len = *), intent(in)  ::  suffix
		character(len = *), intent(in)  ::  string
		character(len = *), intent(out) ::  filename

		character(len = MAX_FILE/2)     ::  tmpstr, tmpsuf
		character(len = MAX_FILE)       ::  fname
		character(len = 15)             ::  fmtstr
		integer                         ::  strlen, suflen

		!reset filename
		filename = ''

        !create temporary string and suffix and remove whitespace
		tmpstr = string
		tmpsuf = suffix
		call remove(tmpstr)
		call remove(tmpsuf)

        !get length of string and suffix
		strlen = len_trim(adjustl(tmpstr))
		suflen = len_trim(adjustl(tmpsuf))

		!construct name of file
		write(fmtstr, '(a,i2,a,i2,a)'), '(a', strlen, ',a', suflen, ')'
		write(fname, fmtstr), tmpstr(1:strlen), tmpsuf(1:suflen)
		call remove(fname)

		filename = fname

	end subroutine build_filename

	!:.........................................................................:

	subroutine warning(message)
		!prints warning message to screen but does not stop program
		!Author: Katherine Holcomb, 2012 v. 1.0
		!Inputs:
		!   message:  message to print
		!Result:
		!	the message is printed to the screen

		 character(len = *), intent(in)  :: message

		 write(*, *) message

	end subroutine warning

	!:.........................................................................:

	subroutine fatal_error(message)
		!stops the program and prints error message to screen
		!Author: Katherine Holcomb, 2012 v. 1.0
		!Inputs:
		!   message:  message to print
		!Result:
		!	the message is printed to the screen and the program is stopped

		 character(len = *), intent(in)  :: message

		 write(*, *) message

		 stop

	end subroutine fatal_error

	!:.........................................................................:

end module FileUtils
