module string

implicit none

!*******************************************************************************
  !
  !This module implements several handy Python-like string functions. It works
  !with ordinary Fortran character variables but requires deferred-size
  !characters from F2003.
  !Author:  Katherine Holcomb
  !Date:    2014-06-26
  !
!*******************************************************************************

   character(len = 26)  :: uppercase = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
   character(len = 26)  :: lowercase = "abcdefghijklmnopqrstuvwxyz"
   character(len = 28)  :: punct = "#$%&()*+,-./:;<=>?@[]^_`{|}~"
   character(len = 2)   :: ws = ' \t'
   character            :: single_quote = "'"
   character            :: double_quote = '"'

contains

	!:.........................................................................:

	subroutine lower(a)
		!converts all uppercase letters in input string to lowercase
		!Inputs/Outputs:
		!	a:  character to be switched to lower case

		character(len = *), intent(inout)  :: a
		integer                            :: n, numchars
		character                          :: c

		!check if string has any characters
		numchars = len(a)
		if (numchars == 0) then
			return
		endif

		!cycle through characters in string and change any uppercase letters to
		  !lowercase
		do n = 1, numchars
			c = a(n:n)
			if (is_alpha(c)) then
				if (is_lower(c)) then
					cycle
				else
				a(n:n) = char(ichar(c) + 32)
				endif
			endif
		enddo

	end subroutine lower

	!:.........................................................................:

	subroutine upper(a)
		!converts all lowercase letters in input string to uppercase
		!Inputs/Outputs:
		!	a:  character to be switched to upper case

		character(len = *), intent(inout)  :: a
		integer                            :: n, numchars
		character                          :: c

		!check if string has any characters
		numchars = len(a)
		if (numchars == 0) then
			return
		endif

		!cycle through characters in string and change any lowercase letters to
		  !uppercase
		do n = 1, numchars
			c = a(n:n)
			if (is_alpha(c)) then
				if (is_upper(c)) then
					cycle
				else
					a(n:n) = char(ichar(c) - 32)
				endif
			endif
		enddo

	end subroutine upper

	!:.........................................................................:

	logical function is_upper(a)
		!checks if a string is all uppercase
		!Inputs:
		!	a:        character to be checked for uppercase
		!Outputs:
		!	is_upper: all upper case?

		character(len = *), intent(in)   :: a
		integer                          :: n, numchars

		!if no characters, can't be uppercase
		numchars = len(a)
		if (numchars == 0) then
			is_upper = .false.
			return
		endif

		!cycle through characters anc check for uppercase letters
		do n = 1, numchars
			!just space or punctuation?
			if (is_space(a(n:n)) .or. is_punctuation(a(n:n))) then
				if (numchars > 1) cycle
			else
				is_upper = .false.
			endif

			!check for all uppercase
			if (verify(a(n:n), uppercase) == 0) then
				is_upper = .true.
			else
				is_upper = .false.
				return
			endif
		enddo

	end function

	!:.........................................................................:

	logical function is_lower(a)
		!checks to see if string is all lowercase
		!Inputs:
		!	a:        character to be checked for lowercase
		!Outputs:
		!	is_lower: all lower case?

		character(len = *), intent(in)  :: a
		integer                         :: n, numchars

		!check to see if string has any characters
		numchars = len(a)
		if (numchars == 0) then
			is_lower = .false.
			return
		endif

		!cycle through characters and check for lowercase letters
		do n = 1, numchars
				!just space or punctuation?
			if (is_space(a(n:n)) .or. is_punctuation(a(n:n))) then
				if ( numchars > 1) cycle
			else
				is_lower = .false.
			endif

			!check for all lowercase
			if (verify(a(n:n), lowercase) == 0) then
				is_lower = .true.
			else
				is_lower = .false.
				return
			endif
		enddo

	end function

	!:.........................................................................:

	logical function is_space(a)
		!checks to see if string is comprised of just space(s)
		!Inputs:
		!	a:        character to be checked for spaces
		!Outputs:
		!	is_space: all spaces?

		character(len = *), intent(in)  :: a
		integer                         ::n, numchars

		!see if string has any characters
		numchars = len(a)
		if (numchars == 0) then
			is_space = .false.
			return
		endif

		!cycle through characters and check for all spaces
		do n = 1, numchars
			if (a(n:n) == ' ' .or. a(n:n) == '\t') then
				is_space = .true.
			else
				is_space = .false.
				return
			endif
		enddo

	end function is_space

	!:.........................................................................:

	logical function is_punctuation(a)
		!checks to see if string is comprised of all punctuation characters
		!Inputs:
		!	a:               character to be checked for punctuation
		!Outputs:
		!	is_punctutation: all punctuation?

		character(len=*), intent(in)  :: a
		integer                       :: n, numchars

		!check to see if string has any characters
		numchars = len(a)
		if (numchars == 0) then
			is_punctuation = .false.
			return
		endif

		!cycle through characters and check for all punctuation
		do n = 1, numchars
			if (verify(a(n:n), punct) == 0 .or. a(n:n) == single_quote .or.    &
						a(n:n) == double_quote)   then
				is_punctuation = .true.
			else
				is_punctuation = .false.
				return
			endif
		enddo

	end function is_punctuation

	!:.........................................................................:

	logical function is_alpha(a)
		!checks to see if string is all alpha characters
		!Inputs:
		!	a:        character to be checked for alpha
		!Outputs:
		!	is_alpha: all alpha characters?

		character(len = *), intent(in)  :: a
		integer                         :: n, numchars

		!check to see if string has any characters
		numchars = len(a)
		if (numchars == 0) return

		!cycle through characters and check for alpha
		do n = 1, numchars
			if (is_upper(a(n:n)) .or. is_lower(a(n:n))) then
				is_alpha=.true.
			else
				is_alpha=.false.
				return
			endif
		enddo

	end function

	!:.........................................................................:

	logical function is_digit(a)
		!checks if string is numbers only
		!Inputs:
		!	a:        character to be checked for numbers
		!Outputs:
		!	is_digit: all numbers?

		character(len=*), intent(in)  :: a
		integer                       :: n, numchars

		!check if string has any characters
		numchars = len(a)
		if (numchars == 0) return

		!cycle through characters and check for all numeric
		do n = 1, numchars
			if ((ichar(a(n:n)) >= 48 .and. ichar(a(n:n)) <= 57)) then
				is_digit = .true.
			else
				is_digit = .false.
				return
			endif
		enddo

	end function

	!:.........................................................................:

	logical function is_alnum(a)
		!check to see if string is only alphanumeric
		!Inputs:
		!	a:        character to be checked for alphanumeric
		!Outputs:
		!	is_alnum: all alphanumeric?

		character(len = *), intent(in)  :: a
		integer                         :: n, numchars

		!check if string has any characters
		numchars = len(a)
		if (numchars == 0) return

		!cycle through characters and check for all alphanumeric
		do n = 1, numchars
			if (is_alpha(a(n:n)) .or. is_digit(a(n:n))) then
				is_alnum = .true.
			else
				is_alnum = .false.
				return
			endif
		enddo

	end function

	!:.........................................................................:

	subroutine split(string, fields, delimiter)
		!splits a string (e.g. row of a file) by a delimiter and returns fields
		  !of values
		!Inputs/Outputs:
		!	string:    string to split
		!Inputs:
		!	delimiter: delimiter by which to split string
		!Outputs:
		!	fields:    array with fields from string split by delimiter

		character(len = *),              intent(inout)  :: string
		character(len = *), optional,    intent(in)     :: delimiter
		character(len = :), allocatable, intent(out)    :: fields(:)

		integer,            allocatable, dimension(:)   :: locs, strlens
		character(len = :), allocatable                 :: delimit
		character(len = 1)                              :: tab = '\t'
		integer                                         :: n, numchars
		integer                                         :: len_delim
		integer                                         :: maxlen
		integer                                         :: num_delims
		integer                                         :: numfields
		integer                                         :: counter

        !deallocate fields if already allocated
		if (allocated(fields)) deallocate(fields)

        !get length of string - check to make sure anything in string
		numchars = len(string)
		if (numchars == 0) then
			allocate(character(1) :: fields(0))
			return
		endif

        !if delimiter present, set to delimit, otherwise set to space
		if (present(delimiter)) then
			allocate(delimit, source = delimiter)
		else
			allocate(character(len = 1) :: delimit)
			delimit = ' '

			!need to remove tabs and compress white space from string if
			  !delimiter is a space
			call remove(string, tab)
			call compress_ws(string)
		endif

        !find locations where delimiter is in string
		call findlocs(locs, string, delimit)
		num_delims = size(locs)
		if (num_delims == 0) then
			!delimiter not found, return empty list
			allocate(character(1) :: fields(0))
			return
		endif

        !get length of delimiter
		len_delim = len(delimit)

		!check to see how many fields there should be
		if (locs(num_delims) + len_delim >= len_trim(string)) then
			!nothing after last delimiter
			numfields = num_delims
		else
			!value after last delimiter
			numfields = num_delims + 1
		endif

		!allocate strlens to number of fields
		allocate(strlens(numfields))

		!get the length of each string within fields
		if (size(locs) == 1) then
			if (locs(1) > 1) then
				strlens(1) = locs(1) - 1
			else
				strlens(1) = 0
			endif
		else
			strlens(1) = locs(1) - 1
			do n = 2, numfields - 1
				strlens(n) = locs(n) - (locs(n-1) + len_delim)
			enddo
			if (numfields > num_delims) then
				strlens(numfields) = len_trim(string(locs(num_delims):))
			else
				strlens(numfields) = locs(num_delims) - (locs(num_delims-1) +  &
					len_delim)
			endif
		endif

		!must be values within the delimiters
		maxlen = maxval(strlens)
		if (maxlen > 0) then
			allocate(character(maxlen) :: fields(numfields))
		else
			allocate(character(1) :: fields(0))
			return
		endif

		!populate fields with strings within delimiter
		counter = 1
		do n = 1, numfields - 1
			fields(n) = string(counter:locs(n) - 1)
			counter = locs(n) + len_delim
		enddo
		if (numfields > size(locs)) then
			fields(numfields) = string(counter:)
		else
			fields(numfields) = string(counter:locs(n) - 1)
		endif

	end subroutine split

	!:.........................................................................:

	subroutine join(string, fields, delimiter)
		!joins a string together to form a delimeter-separated row of data
		!Inputs:
		!	fields:    array of values to join
		!	delimiter: delimiter by which to split string
		!Outputs:
		!	string: string with values from fields separated by delimiter

		character(len = *), dimension(:),  intent(in)  :: fields
		character(len = *), optional,      intent(in)  :: delimiter
		character(len = :), allocatable,   intent(out) :: string

		character(len = :), allocatable                :: delimit
		integer                                        :: n, numchars
		integer                                        :: numfields
		integer                                        :: len_delim
		integer                                        :: maxlen
		integer                                        :: counter
		integer                                        :: nextfield_len

		!use space if no delimiter present
		if (present(delimiter)) then
			allocate(delimit, source = delimiter)
		else
			allocate(character(1) :: delimit)
			delimit = ' '
		endif
		len_delim = len(delimit)

		!get number of fields to join
		numfields = size(fields)

		!must be some values to join
		if (numfields == 0) then
			allocate(character(len = 1) :: string)
			string = ''
			return
		endif

		!maximum length within fields
		maxlen = maxval(len_trim(fields(:)))

		!number of characters to join
		numchars = sum(len_trim(fields(:)))

		allocate(character(numchars + (numfields - 1)*len_delim) :: string)

		!add fields to string and split by delimiter
		counter = 1
		do n = 1, numfields - 1
			nextfield_len = len_trim(fields(n)) + len_delim
			string(counter:counter + nextfield_len) = trim(fields(n))//delimit
			counter = counter + nextfield_len
		enddo
		string(counter:) = trim(fields(numfields))

	end subroutine join

	!:.........................................................................:

	subroutine replace(string, old, new)
		!replaces character(s) in string with new character(s). string will have
          !to be long enough to accommodate changes in this version
        !Inputs/Outputs:
		!	string:  input string to change
		!Inputs:
		!	old:     character to replace
		!	new:     replacement character

		character(len = *), intent(inout)  :: string
		character(len = *), intent(in)     :: old, new

		character(len = len(string))       :: tmpstring
		integer, allocatable, dimension(:) :: locs
		integer                            :: len_old, len_new
		integer                            :: j, m, n, numchars
		logical                            :: insert

		!get length string, new, and old substrings
		numchars = len(string)
		len_old = len(old); len_new = len(new)

		!check to be sure strings have characters
		if (numchars == 0) return
		if (len_old == 0 .or. len_new == 0) return

		!find locations of old string within input string
		call findlocs(locs, string, old)
		if (size(locs) == 0) return

		tmpstring = ''

		!replace old with new
		n = 1
		m = 1
		do while (m <= numchars)
			do j = 1, size(locs)
				insert = .false.
				if (m == locs(j)) then
				   insert = .true.
				   exit
				endif
			enddo
			if (insert) then
				tmpstring(n:n + len_new - 1) = new
				n = n + len_new
				m = m + len_old
			else
				tmpstring(n:n) = string(m:m)
				n = n + 1
				m = m + 1
			endif
		enddo

		string = trim(adjustl(tmpstring))

	end subroutine replace

	!:.........................................................................:

	subroutine strip(string, substr)
		!strips leading and trailing ends of string of substring
		!Inputs/Outputs:
		!	string:  input string to change
		!Inputs:
		!	substr:  character(s) to get rid of

		character(len = *), intent(inout)  :: string
		character(len = *), intent(in)     :: substr

		integer                            :: sub_len
		integer                            :: numchars

		!get length of string and substring
		numchars = len(string)
		sub_len = len(substr)

		if (numchars == 0) return
		if (sub_len == 0 ) return

		!strip leading and trailing ends
		call lstrip(string, substr)
		call rstrip(string, substr)

	end subroutine strip

	!:.........................................................................:

	subroutine lstrip(string, substr)
		!strips left side of string of substring
		!Inputs/Outputs:
		!	string:  input string to change
		!Inputs:
		!	substr:  character(s) to get rid of

		character(len = *), intent(inout)  :: string
		character(len = *), intent(in)     :: substr

		character(len = len(string))       :: tmpstring
		integer, allocatable, dimension(:) :: locs
		integer                            :: sub_len
		integer                            :: numchars

		!get length of string and substring
		numchars = len(string)
		sub_len = len(substr)

		if (numchars == 0) return
		if (sub_len == 0) return

		!find locations of substring within string
		call findlocs(locs, string, substr)
		if (size(locs) == 0) return

		!string leading edge of string of substring
		if (locs(1) == 1) then
			tmpstring = string(sub_len + 1:)
		else
			tmpstring = string
		endif

		string = tmpstring

	end subroutine lstrip

	!:.........................................................................:

	subroutine rstrip(string, substr)
		!strips right side of string of substring
		!Inputs/Outputs:
		!	string:  input string to change
		!Inputs:
		!	substr:  character(s) to get rid of

		character(len = *), intent(inout)  :: string
		character(len = *), intent(in)     :: substr

		character(len = len(string))       :: tmpstring
		integer, allocatable, dimension(:) :: locs
		integer                            :: sub_len
		integer                            :: numchars

		!get length of string and substring
		numchars = len_trim(string)
		sub_len = len(substr)

		if (numchars == 0) return
		if (sub_len == 0) return


		!find locations of substring within string
		call findlocs(locs, string, substr)
		if (size(locs) == 0) return

		!string trailing edge of string of substring
		if (locs(size(locs)) == (numchars - sub_len + 1)) then
			tmpstring = string(:numchars-sub_len)
		else
			tmpstring = string
		endif

		string = tmpstring

	end subroutine rstrip

	!:.........................................................................:

	subroutine remove(string, set)
		!removes set from an input string
		!Inputs/Outputs:
		!	string:  input string to change
		!Inputs:
		!	set:     character(s) to get rid of

		character(len = *),           intent(inout)  :: string
		character(len = *), optional, intent(in)     :: set

		character(len = len(string))                 :: tmpstring
		character(len = :), allocatable              :: strset
		character                                    :: c
		integer                                      :: i, n, numchars
		integer                                      :: len_set

		!get length of input string, ingoring trailing blanks
		numchars = len_trim(string)

		!exit if no characters in string
		if (numchars == 0) return

		!if set is present, allocate strset to that, otherwise allocate to
         !whitespace
		if (present(set)) then
			allocate(strset, source = set)
		else
			allocate(strset, source = ws)
		endif

        !get length of strset and exit if 0
		len_set = len(strset)
		if (len_set == 0) return

        !loop through characters in input string. Check each character to see if
          !it is in the set to be removed. If not, add it to a new string.
		tmpstring = ''
		n = 0
		do i = 1, numchars
			c = string(i:i)
			if (is_in(c, strset)) then
				cycle
			else
				n = n + 1
				tmpstring(n:n) = c
			endif
		enddo

        !return new string without values in 'set', trimmed of whitespace and
          !adjusted left
		string = trim(adjustl(tmpstring))

	end subroutine remove

	!:.........................................................................:

	integer function find(string, substr, start, last)
		!finds first integer location of substring in string
		!may search within set location in string by using optional start, last
		!Inputs:
		!	string: input string
		!	substr: string to find within substring
		!Optional Inputs:
		!	start:  starting index to search within string
		!	last:	ending index to search within string
		!Outputs:
		!	find: integer location of substring in string

		character(len = *),  intent(in)  :: string, substr
		integer, optional,   intent(in)  :: start, last

		integer                          :: st, lst, numchars
		integer                          :: i, sub_len

		!get length of string and substring
		numchars = len(string)
		sub_len = len(substr)

		!must have length > 0
		if (numchars== 0) return
		if (sub_len == 0) return

		!set starting and ending indices if present
		if (present(start)) then
			st = start
		else
			st = 1
		endif

		if (present(last)) then
			lst = last
		else
			lst = numchars
		endif

		!find first instance of substr within string using set st/lst indices
		find = 0
		do i = st, lst
			if (string(i:i+sub_len-1) == substr) then
				find = i
				exit
			endif
		enddo

	end function find

	!:.........................................................................:

	logical function starts_with(string, substr)
		!checks whether input string starts with substring
		!Inputs:
		!	string:       input string to check
		!	substr:       substring to check against
		!Outputs:
		!	starts_with:  does the string start with substr?

		character(len = *), intent(in)  :: string
		character(len = *), intent(in)  :: substr

		if (find(string, substr) == 1) then
			starts_with = .true.
		else
			starts_with = .false.
		endif

   end function starts_with

	!:.........................................................................:

	logical function ends_with(string, substr)
		!checks whether input string ends with substring
		!Inputs:
		!	string:     input string to check
		!	substr:     substring to check against
		!Outputs:
		!	ends_with:  does the string start with substr?

		character(len = *), intent(in)  :: string
		character(len = *), intent(in)  :: substr

		if (find(string, substr) == (len(string) - len(substr) + 1)) then
			ends_with = .true.
		else
			ends_with = .false.
		endif

	end function ends_with

	!:.........................................................................:

	logical function is_in(string, set)
		!checks whether input string is in set of characters
		!Inputs:
		!	string:  input string to check
		!	set:     set of characters to check against
		!Outputs:
		!	is_in:   is string in set of characters?

		character(len = *), intent(in)  :: string
		character(len = *), intent(in)  :: set

		if (scan(string, set) == 0) then
			is_in = .false.
		else
			is_in = .true.
		endif

	end function is_in

	!:.........................................................................:

	subroutine compress_ws(string)
		!compresses white space in string
		!Inputs/Outputs:
		!	string:  input string to change

		character(len = *), intent(inout) :: string

		character(len = len(string))      :: tmpstring
		character                         :: ci
		integer                           :: i, n, numchars

		!adjust left and get length of string, ignoring trailing blanks
		string = adjustl(string)
		numchars = len_trim(string)

		!must have some characters
		if (numchars == 0) return

		tmpstring = ''
		n = 0

		!loop through number of characters, checking if each individual
		  !character and the following is whitespace
		do i = 1, numchars - 1
			ci = string(i:i)
			if (verify(ci, ws) == 0 .and.                                      &
				verify(next_char(string, i), ws) == 0) then
				cycle
			else
				n = n + 1
				if (ci == '\t' ) ci = ' '
				tmpstring(n:n) = ci
			endif
		enddo

		tmpstring(n+1:n+1) = string(numchars:numchars)

		string = trim(adjustl(tmpstring))

	end subroutine compress_ws

!:.............................................................................:

	subroutine findlocs(locs, string, substr)
		!finds locations of substring within string
		!Inputs:
		!	string: input string
		!	substr: substring to find within string
		!Outputs:
		!	locs:   index locations within string with substring

		integer, dimension(:), allocatable, intent(out) :: locs
		character(len = *),                 intent(in)  :: string
		character(len = *),                 intent(in)  :: substr

		integer                                         :: counter
		integer                                         :: pos, sub_len
		integer                                         :: numchars
		integer, dimension(:), allocatable              :: templocs

		numchars = len(string)
		sub_len  = len(substr)

		!must find substring somewhere within string
		if (find(string, substr) == 0) then
			allocate(locs(0))
			return
		endif

		!allocate locs to 1 to start
		if (.not. allocated(locs)) allocate(locs(1))
		!find first instance of substr
		pos = find(string, substr)

		!set the first value to first instance
		locs(1) = pos

		!move forward in string by length of substring
		counter = pos + sub_len

		!find all instances
		do while (pos > 0)
			pos = find(string, substr, counter)
			if (pos > 0) then
				allocate(templocs(size(locs)))
				templocs = locs(:)
				deallocate(locs)
				allocate(locs(size(templocs) + 1))
				locs(:size(templocs)) = templocs(:)
				deallocate(templocs)
				locs(size(locs)) = pos
				counter = pos + sub_len
			endif
		enddo

	end subroutine findlocs

	!:.........................................................................:

	character function next_char(string, l)
		!gets character after l in string
		!Inputs:
		!	string:    input string
		!	l:         index to get following character
		!Outputs:
		!	next_char:  character after l in string

		character(len = *), intent(in)  :: string
		integer,            intent(in)  :: l

		integer                         :: next

		next = l + 1
		if (next < len(string)) then
			next_char = string(next:next)
		else if (next < 0) then
           next_char = ''
        else
           next_char = ''
        endif

	end function

	!:.........................................................................:

end module string
