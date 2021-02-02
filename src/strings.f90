module string

!*******************************************************************************
  !
  !This module implements several handy Python-like string functions. It works
  !with ordinary Fortran character variables but requires deferred-size
  !characters from F2003.
  !
!*******************************************************************************

    implicit none

    ! Data dictionary, constants:
    character(len = 26), parameter :: UPPERCASE = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" ! Upper case letters
    character(len = 26), parameter :: LOWERCASE = "abcdefghijklmnopqrstuvwxyz" ! Lower case letters
    character(len = 28), parameter :: PUNCT = "#$%&()*+,-./:;<=>?@[]^_`{|}~"   ! Punctuation
    character(len = 2),  parameter :: WS = ' \t'                               ! White space
    character,           parameter :: SINGLE_QUOTE = "'"                       ! Single quote
    character,           parameter :: DOUBLE_QUOTE = '"'                       ! Double quote
    character,           parameter :: TAB = '\t'                               ! Tab

contains

    !:.........................................................................:

    subroutine lower(a)
        !
        !  Converts all uppercase letters in the input string to lowercase
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/14     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *), intent(inout) :: a ! Input string

        ! Data dictionary: local variables
        integer   :: n        ! Looping index
        integer   :: numchars ! Length of string
        character :: c        ! Temporary string

        ! Check length of string
        numchars = len(a)
        if (numchars == 0) then
            return
        endif

        ! Cycle through characters in string and change any uppercase letters to
        ! lowercase
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
        !
        !  Converts all lowercase letters in the input string to uppercase
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/14     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *), intent(inout) :: a ! Input string

        ! Data dictionary: local variables
        integer   :: n        ! Looping index
        integer   :: numchars ! Length of string
        character :: c        ! Temporary string

        ! Check length of string
        numchars = len(a)
        if (numchars == 0) then
            return
        endif

        ! Cycle through characters in string and change any lowercase letters to
        ! uppercase
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
        !
        !  Checks if string is all uppercase
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/14     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *), intent(in) :: a ! Input string

        ! Data dictionary: local variables
        integer :: n        ! Looping index
        integer :: numchars ! Length of string

        ! Check length of string
        numchars = len(a)
        if (numchars == 0) then
            is_upper = .false.
            return
        endif

        ! Cycle through characters and check for uppercase letters
        do n = 1, numchars

            ! Space or punctuation?
            if (is_space(a(n:n)) .or. is_punctuation(a(n:n))) then
                if (numchars > 1) cycle
            else
                is_upper = .false.
            endif

            ! Check for all uppercase
            if (verify(a(n:n), UPPERCASE) == 0) then
                is_upper = .true.
            else
                is_upper = .false.
                return
            endif
        enddo

    end function

    !:.........................................................................:

    logical function is_lower(a)
        !
        !  Checks if string is all lowercase
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/14     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *), intent(in) :: a ! Input string

        ! Data dictionary: local variables
        integer :: n        ! Looping index
        integer :: numchars ! Length of string

        ! Check length of string
        numchars = len(a)
        if (numchars == 0) then
            is_lower = .false.
            return
        endif

        ! Cycle through characters and check for lowercase letters
        do n = 1, numchars
            ! Space or punctuation?
            if (is_space(a(n:n)) .or. is_punctuation(a(n:n))) then
                if ( numchars > 1) cycle
            else
                is_lower = .false.
            endif

            ! Check for all lowercase
            if (verify(a(n:n), LOWERCASE) == 0) then
                is_lower = .true.
            else
                is_lower = .false.
                return
            endif
        enddo

    end function

    !:.........................................................................:

    logical function is_space(a)
        !
        !  Checks if string is comprised of just space(s)
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/14     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *), intent(in) :: a ! Input string

        ! Data dictionary: local variables
        integer :: n        ! Looping index
        integer :: numchars ! Length of string

        ! Check length of string
        numchars = len(a)
        if (numchars == 0) then
            is_space = .false.
            return
        endif

        ! Cycle through characters and check for spaces
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
        !
        !  Checks if string is comprised of just punctuation characters
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/14     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len=*), intent(in) :: a ! Input string

        ! Data dictionary: local variables
        integer :: n        ! Looping index
        integer :: numchars ! Length of string

        ! Check length of string
        numchars = len(a)
        if (numchars == 0) then
            is_punctuation = .false.
            return
        endif

        ! Cycle through characters and check for punctuation
        do n = 1, numchars
            if (verify(a(n:n), PUNCT) == 0 .or. a(n:n) == SINGLE_QUOTE .or.    &
                        a(n:n) == DOUBLE_QUOTE)   then
                is_punctuation = .true.
            else
                is_punctuation = .false.
                return
            endif
        enddo

    end function is_punctuation

    !:.........................................................................:

    logical function is_alpha(a)
        !
        !  Checks if string is alpha only
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/14     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *), intent(in) :: a ! Input string

        ! Data dictionary: local variables
        integer :: n        ! Looping index
        integer :: numchars ! Length of string

        ! Check length of string
        numchars = len(a)
        if (numchars == 0) return

        ! Cycle through characters and check for alpha
        do n = 1, numchars
            if (is_upper(a(n:n)) .or. is_lower(a(n:n))) then
                is_alpha = .true.
            else
                is_alpha = .false.
                return
            endif
        enddo

    end function

    !:.........................................................................:

    logical function is_digit(a)
        !
        !  Checks if string is numeric only
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/14     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len=*), intent(in) :: a ! Input string

        ! Data dictionary: local variables
        integer :: n        ! Looping index
        integer :: numchars ! Length of string

        ! Check length of string
        numchars = len(a)
        if (numchars == 0) return

        ! Cycle through characters and check for numeric
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
        !
        !  Checks if string is alphanumeric
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/14     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *), intent(in) :: a ! Input string

        ! Data dictionary: local variables
        integer :: n        ! Looping index
        integer :: numchars ! Length of string

        ! Check length of string
        numchars = len(a)
        if (numchars == 0) return

        ! Cycle through characters and check for alphanumeric
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
        !
        !  Splits a string (e.g. row of a file) by a delimiter and returns
        !   an array of values.
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/14     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *),              intent(inout) :: string    ! Input string
        character(len = *), optional,    intent(in)    :: delimiter ! Delimiter
        character(len = :), allocatable, intent(out)   :: fields(:) ! Array with fields from split string

        ! Data dictionary: local variables
        integer,            allocatable, dimension(:) :: locs       ! Locations of the deliminter in input string
        integer,            allocatable, dimension(:) :: strlens    ! Length of each string in the output array
        character(len = :), allocatable               :: delimit    ! Delimiter
        integer                                       :: n          ! Looping index
        integer                                       :: numchars   ! Length of input string
        integer                                       :: len_delim  ! Length of delimiter
        integer                                       :: maxlen     ! Maximum length of strings in output array
        integer                                       :: num_delims ! Number of delimiters found in input string
        integer                                       :: numfields  ! Number of fields in output array
        integer                                       :: counter    ! Looping index

        ! Deallocate fields if already allocated
        if (allocated(fields)) deallocate(fields)

        ! Get length of string - check to make sure anything in string
        numchars = len(string)
        if (numchars == 0) then
            allocate(character(1) :: fields(0))
            return
        endif

        ! Ff delimiter present, set to delimit, otherwise set to space
        if (present(delimiter)) then
            allocate(delimit, source = delimiter)
        else
            allocate(character(len = 1) :: delimit)
            delimit = ' '

            ! Need to remove tabs and compress white space from string if
            ! delimiter is a space
            call remove(string, TAB)
            call compress_ws(string)
        endif

        ! Find locations where delimiter is in string
        call findlocs(locs, string, delimit)
        num_delims = size(locs)
        if (num_delims == 0) then
            ! Delimiter not found, return empty list
            allocate(character(1) :: fields(0))
            return
        endif

        ! Get length of delimiter
        len_delim = len(delimit)

        ! Check to see how many fields there should be
        if (locs(num_delims) + len_delim >= len_trim(string)) then
            ! Nothing after last delimiter
            numfields = num_delims
        else
            ! Value after last delimiter
            numfields = num_delims + 1
        endif

        ! Allocate strlens to number of fields
        allocate(strlens(numfields))

        ! Get the length of each string within fields
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

        ! Must be values within the delimiters
        maxlen = maxval(strlens)
        if (maxlen > 0) then
            allocate(character(maxlen) :: fields(numfields))
        else
            allocate(character(1) :: fields(0))
            return
        endif

        ! Populate fields with strings within delimiter
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
        !
        !  Joins an array of strings together to form a delimeter-separated row
        !   of data
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/14     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *), dimension(:), intent(in)  :: fields    ! Array to join
        character(len = *), optional,     intent(in)  :: delimiter ! Delimiter
        character(len = :), allocatable,  intent(out) :: string    ! Output string

        ! Data dictionary: local variables
        character(len = :), allocatable :: delimit       ! Delimiter
        integer                         :: n             ! Looping index
        integer                         :: numchars      ! Number of characters to join
        integer                         :: numfields     ! Number of fields to join
        integer                         :: len_delim     ! Length of delimiter
        integer                         :: maxlen        ! Maximum length within fields
        integer                         :: counter       ! Looping index
        integer                         :: nextfield_len ! Length of next field

        ! Use space if no delimiter present
        if (present(delimiter)) then
            allocate(delimit, source = delimiter)
        else
            allocate(character(1) :: delimit)
            delimit = ' '
        endif
        len_delim = len(delimit)

        ! Get number of fields to join
        numfields = size(fields)

        ! Must be some values to join
        if (numfields == 0) then
            allocate(character(len = 1) :: string)
            string = ''
            return
        endif

        ! Maximum length within fields
        maxlen = maxval(len_trim(fields(:)))

        ! Number of characters to join
        numchars = sum(len_trim(fields(:)))

        ! Allocate correct length to output string
        allocate(character(numchars + (numfields - 1)*len_delim) :: string)

        ! Add fields to string and split by delimiter
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
        !
        !  Replaces character(s) in a string with new character(s).
        !  String will have to be long enough to accomodate changes.
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/14     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *), intent(inout) :: string ! String to modify
        character(len = *), intent(in)    :: old    ! Characters to be replaced
        character(len=*),   intent(in)    :: new    ! Charcters to replace with

        ! Data dictionary: local variables
        character(len = len(string))       :: tmpstring ! Temporary string
        integer, allocatable, dimension(:) :: locs      ! Locations of old string within input string
        integer                            :: len_old   ! Length of old substring
        integer                            :: len_new   ! Length of new substring
        integer                            :: j, m, n   ! Looping indices
        integer                            :: numchars  ! Length of input string
        logical                            :: insert    ! Insert or no?

        ! Get length string, new, and old substrings
        numchars = len(string)
        len_old = len(old)
        len_new = len(new)

        ! Check to be sure strings have length
        if (numchars == 0) return
        if (len_old == 0 .or. len_new == 0) return

        ! Find locations of old string within input string
        call findlocs(locs, string, old)
        if (size(locs) == 0) return

        tmpstring = ''

        ! Replace old with new
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

        ! Trim whitespace
        string = trim(adjustl(tmpstring))

    end subroutine replace

    !:.........................................................................:

    subroutine strip(string, substr)
        !
        !  Strips leading and trailing ends of string of an input substring
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/14     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *), intent(inout) :: string ! Input string
        character(len = *), intent(in)    :: substr ! Input substring

        ! Data dictionary: local variables
        integer :: sub_len  ! Length of input substring
        integer :: numchars ! Length of input string

        ! Check length of string and substring
        numchars = len(string)
        sub_len = len(substr)
        if (numchars == 0) return
        if (sub_len == 0) return

        ! Strip leading and trailing ends
        call lstrip(string, substr)
        call rstrip(string, substr)

    end subroutine strip

    !:.........................................................................:

    subroutine lstrip(string, substr)
        !
        !  Strips left side of string of an input substring
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/14     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *), intent(inout) :: string ! Input string
        character(len = *), intent(in)    :: substr ! Input substring

        ! Data dictionary: local variables
        character(len = len(string))       :: tmpstring ! Temporary string
        integer, allocatable, dimension(:) :: locs      ! Locations of substring in string
        integer                            :: sub_len   ! Length of substring
        integer                            :: numchars  ! Length of string

        ! Check length of string and substring
        numchars = len(string)
        sub_len = len(substr)
        if (numchars == 0) return
        if (sub_len == 0) return

        ! Find locations of substring within string
        call findlocs(locs, string, substr)
        if (size(locs) == 0) return

        ! Strip leading edge of string of substring
        if (locs(1) == 1) then
            tmpstring = string(sub_len + 1:)
        else
            tmpstring = string
        endif

        string = tmpstring

    end subroutine lstrip

    !:.........................................................................:

    subroutine rstrip(string, substr)
        !
        !  Strips right side of string of an input substring
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/14     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *), intent(inout) :: string ! Input string
        character(len = *), intent(in)    :: substr ! Input substring

        ! Data dictionary: local variables
        character(len = len(string))       :: tmpstring ! Temporary string
        integer, allocatable, dimension(:) :: locs      ! Locations of substring within string
        integer                            :: sub_len   ! Length of substring
        integer                            :: numchars  ! Length of string

        ! Check length of string and substring
        numchars = len_trim(string)
        sub_len = len(substr)
        if (numchars == 0) return
        if (sub_len == 0) return


        ! Find locations of substring within string
        call findlocs(locs, string, substr)
        if (size(locs) == 0) return

        ! Strip trailing edge of string of substring
        if (locs(size(locs)) == (numchars - sub_len + 1)) then
            tmpstring = string(:numchars-sub_len)
        else
            tmpstring = string
        endif

        string = tmpstring

    end subroutine rstrip

    !:.........................................................................:

    subroutine remove(string, set)
        !
        !  Removes an input substring from an input string
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/14     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *),           intent(inout) :: string ! Input string
        character(len = *), optional, intent(in)    :: set    ! Input substring (default to whitespace)

        ! Data dictionary: local variables
        character(len = len(string))    :: tmpstring ! Temporary string
        character(len = :), allocatable :: strset    ! String to remove
        character                       :: c         ! Temporary string
        integer                         :: i, n      ! Looping indices
        integer                         :: numchars  ! Length of input string
        integer                         :: len_set   ! Length of substring

        ! Get length of input string, ingoring trailing blanks
        numchars = len_trim(string)

        ! Exit if no characters in string
        if (numchars == 0) return

        ! If set is present, allocate strset to that, otherwise allocate to
        ! whitespace
        if (present(set)) then
            allocate(strset, source = set)
        else
            allocate(strset, source = WS)
        endif

        ! Get length of strset and exit if 0
        len_set = len(strset)
        if (len_set == 0) return

        ! Loop through characters in input string. Check each character to see if
        ! it is in the set to be removed. If not, add it to a new string.
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

        ! Return new string without values in 'set', trimmed of whitespace
        string = trim(adjustl(tmpstring))

    end subroutine remove

    !:.........................................................................:

    integer function find(string, substr, start, last)
        !
        !  Finds the first location of a substring in string
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/14     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *), intent(in) :: string ! Input string
        character(len = *), intent(in) :: substr ! Input substring
        integer, optional,  intent(in) :: start  ! Optional starting location to search
        integer, optional,  intent(in) :: last   ! Optional ending location to search

        ! Data dictionary: local variables
        integer :: st       ! Starting location to search
        integer :: lst      ! Ending location to search
        integer :: numchars ! Length of string
        integer :: i        ! Looping index
        integer :: sub_len  ! Length of substring

        ! Check length of string and substring
        numchars = len(string)
        sub_len = len(substr)
        if (numchars== 0) return
        if (sub_len == 0) return

        ! Set starting and ending indices if present, otherwise search whole
        ! string
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

        ! Find first instance of substr within string using set st/lst indices
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
        !
        !  Checks whether input string starts with substring
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/14     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *), intent(in) :: string ! Input string
        character(len = *), intent(in) :: substr ! Input substring

        ! Check if it's the first location
        if (find(string, substr) == 1) then
            starts_with = .true.
        else
            starts_with = .false.
        endif

   end function starts_with

    !:.........................................................................:

    logical function ends_with(string, substr)
        !
        !  Checks whether input string ends with substring
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/14     K. Holcomb           Original Code
        !

        !checks whether input string ends with substring
        !Inputs:
        !    string:     input string to check
        !    substr:     substring to check against
        !Outputs:
        !    ends_with:  does the string start with substr?

        ! Data dictionary: calling arguments
        character(len = *), intent(in) :: string ! Input string
        character(len = *), intent(in) :: substr ! Input substring

        ! Check if it ends with that string
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
        !    string:  input string to check
        !    set:     set of characters to check against
        !Outputs:
        !    is_in:   is string in set of characters?

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
        !
        !  Deletes white space in a string
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/14     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        character(len = *), intent(inout) :: string ! Input string

        ! Data dictionary: local variables
        character(len = len(string)) :: tmpstring ! Temporary string
        character                    :: ci        ! Temporary string (sub of input)
        integer                      :: i, n      ! Looping indices
        integer                      :: numchars  ! Length of string

        ! Adjust left and get length of string, ignoring trailing blanks
        string = adjustl(string)
        numchars = len_trim(string)

        ! Must have some characters
        if (numchars == 0) return

        tmpstring = ''
        n = 0
        ! Loop through number of characters, checking if each individual
        ! character and the following is whitespace
        do i = 1, numchars - 1
            ci = string(i:i)
            if (verify(ci, WS) == 0 .and.                                      &
                verify(next_char(string, i), WS) == 0) then
                cycle
            else
                n = n + 1
                if (ci == '\t' ) ci = ' '
                tmpstring(n:n) = ci
            endif
        enddo

        ! Set to new string and trim
        tmpstring(n+1:n+1) = string(numchars:numchars)
        string = trim(adjustl(tmpstring))

    end subroutine compress_ws

!:.............................................................................:

    subroutine findlocs(locs, string, substr)
        !
        !  Finds locations of substring within string
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/14     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        integer, dimension(:), allocatable, intent(out) :: locs   ! Array of location of substring
        character(len = *),                 intent(in)  :: string ! Input string
        character(len = *),                 intent(in)  :: substr ! Input substring

        ! Data dictionary: local variables
        integer                            :: counter  ! Looping index
        integer                            :: pos      ! Current position in string
        integer                            :: sub_len  ! Length of substring
        integer                            :: numchars ! Length of input string
        integer, dimension(:), allocatable :: templocs ! Temporary array of locations

        ! Get length of string and substring
        numchars = len(string)
        sub_len  = len(substr)

        ! Must find substring somewhere within string
        if (find(string, substr) == 0) then
            allocate(locs(0))
            return
        endif

        ! Allocate locs to 1 to start
        if (.not. allocated(locs)) allocate(locs(1))

        ! Find first instance of substr
        pos = find(string, substr)

        ! Set the first value to first instance
        locs(1) = pos

        ! Move forward in string by length of substring
        counter = pos + sub_len

        ! Now find all instances
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
        !
        !  Gets the character after input index (l) in string
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    06/26/14     K. Holcomb           Original Code
        !

        ! Data dictonary: calling arguments
        character(len = *), intent(in) :: string ! Input string
        integer,            intent(in) :: l      ! Index to check

        ! Data dictionary: local variables
        integer :: next ! Looping index

        ! Check the next string
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
