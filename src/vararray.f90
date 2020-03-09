module DATA_MODULE
  use Constants
  
!*******************************************************************************
  !
  !This module defines the CHAR_DATA type to be used in the vararray module.
  !
!*******************************************************************************

	type CHAR_DATA
		character(len = MAX_NLEN)  ::  value
	end type CHAR_DATA

	type(CHAR_DATA), parameter :: empty_data = CHAR_DATA('')

end module DATA_MODULE


module vararray
  use DATA_MODULE, LIST_DATA => CHAR_DATA, empty_list_data => empty_data

!*******************************************************************************
  !
  !This module implements dynamically growing arrays
  !
  !It is essentially Arjen Markus's vector.f90 module but we have renamed it to 
  !vararray, because most scientists don't care about fine distinctions between 
  !'linked lists' and C++-type 'vectors'.
  !This object behaves like a Python/Perl/etc. list.
  
  !The module is straightforward: it defines a suitable data structure.
  !Data can be added to the list and you can retrieve data from it
  
  !Note: 
  !  For the function list_at() we need a parameter that represents the 'empty 
  !	 list data' value.
  !
  !Id: vectors.f90, v1.4 12/06/2008 relaxmike Exp
  !Written by Arjen Markis
  !
!*******************************************************************************

type LIST
    private
    integer                                  :: not_used
    type(LIST_DATA), dimension(:), pointer   :: data    => null()
end type list

private
public :: LIST
public :: LIST_DATA
public :: list_create
public :: list_append
public :: list_at
public :: list_size
public :: list_put
public :: list_delete_elements
public :: list_destroy
public :: list_insert_empty

real, parameter :: growth_rate = 1.1

contains

	!:.........................................................................:
	
	subroutine list_create(lst, capacity)
		!creates a new list
		!Author: Katherine Holcomb 2012 v. 1
		!Inputs:
		!  lst:       variable that should hold the list
		!  capacity:  initial capacity 
		!Note:
	    !  !the fields of the list data structure are set
		
		type(LIST)         :: lst
		integer, optional  :: capacity

		integer            :: cap

		!check that the list does not have any data left
		if (associated(lst%data)) then
			call list_destroy(lst)
		endif

		if (present(capacity)) then
			cap = max(1, capacity)
		else
			cap = 10
		endif
		
		allocate(lst%data(1:cap))
		lst%not_used = 0
	end subroutine list_create
	
	!:.........................................................................:

	subroutine list_destroy(lst)
		!destroys a list
		!Author: Katherine Holcomb 2012 v. 1
		!Inputs:
		!  lst: input list
		
		type(LIST) :: lst
		
		!check that the list does not have any data left
		if (associated(lst%data)) then
			deallocate(lst%data)
		endif
		lst%not_used = 0
	end subroutine list_destroy
	
	!:.........................................................................:

	integer function list_size(lst)
		!returns the number of elements in use
		!Author: Katherine Holcomb 2012 v. 1
		!Inputs:
		!  lst:         input list
		!Outputs:
		!	list_size:  size of list
	
		type(LIST) :: lst

		list_size = lst%not_used
	end function list_size
	
	!:.........................................................................:

	type(LIST_DATA) function list_at(lst, n)
		!gets the value of the nth element of the list
		!Author: Katherine Holcomb 2012 v. 1
		!Inputs:
		!	lst:      input list
		!	n:        index to get
		!Outputs:
		!	list_at:  value of nth element
	
		type(LIST) :: lst
		integer    :: n

		if (n .lt. 1 .or. n .gt. lst%not_used) then
			list_at = empty_list_data
		else
			list_at = lst%data(n)
		endif
	end function list_at
	
	!:.........................................................................:

	subroutine list_insert_empty(lst, pos, number)
		!inserts one or more empty elements
		!Author: Katherine Holcomb 2012 v. 1
		!Inputs:
		!   lst:     list in question
		!   pos:     position to insert the emtpy elements
		!   number:  number of empty elements
	
		type(LIST)          :: lst
		integer, intent(in) :: pos
		integer, intent(in) :: number
		
		integer             :: i

		if (number .lt. 1 .or. pos .lt. 1 .or. pos .gt. lst%not_used) then
			return
		endif

		if (lst%not_used + number .ge. size(lst%data)) then
			call list_increase_capacity(lst, lst%not_used + number)
		endif

		do i = lst%not_used, pos, -1
			lst%data(i + number) = lst%data(i)
		enddo
		
		do i = 1, number
			lst%data(pos + i - 1) = empty_list_data
		enddo
		
		lst%not_used = lst%not_used + number
		
	end subroutine list_insert_empty
	
	!:.........................................................................:

	subroutine list_delete_elements(lst, pos, number)
		!deletes one or more elements
		!Author: Katherine Holcomb 2012 v. 1
		!Inputs:
		!	lst:     list in question
		!	pos:     position to start deletion
		! 	number:  number of elements

		type(LIST)          :: lst
		integer, intent(in) :: pos
		integer, intent(in) :: number
		
		integer             :: i

		if (number .lt. 1 .or. pos .lt. 1 .or. pos .gt. lst%not_used) then
			return
		endif

		do i = pos, lst%not_used - number
			lst%data(i) = lst%data(i + number)
		enddo
		
		lst%not_used = lst%not_used - number
		
	end subroutine list_delete_elements
	
	!:.........................................................................:

	subroutine list_append(lst, data)
		!appends data to a list
		!Author: Katherine Holcomb 2012 v. 1
		!Inputs:
		!	lst:     list in question
		!	pos:     data to append

		type(LIST)      :: lst
		type(LIST_DATA) :: data

		if (lst%not_used .ge. size(lst%data)) then
			call list_increase_capacity(lst, lst%not_used + 1)
		endif
		
		lst%not_used = lst%not_used + 1
		lst%data(lst%not_used) = data
		
	end subroutine list_append
	
	!:.........................................................................:

	subroutine list_put(lst, n, data)
		!puts a value at a specific element of the list
		  !(it does not need to exist yet)
		!Author: Katherine Holcomb 2012 v. 1.0
		!Inputs:
		!	lst:    list in question
		!	n:      index of element
		!  	data:   data to be put in the list

		type(LIST)      :: lst
		integer         :: n
		type(LIST_DATA) :: data

		if (n .lt. 1) then
			return
		endif
		
		if (n .gt. size(lst%data)) then
			call list_increase_capacity(lst, n)
		endif

		lst%not_used = max(lst%not_used, n)
		lst%data(n) = data
		
	end subroutine list_put
	
	!:.........................................................................:

	subroutine list_increase_capacity(lst, capacity)
		!expands the array holding the data
		!Author: Katherine Holcomb 2012 v. 1.0
		!Inputs:
		!	lst:       list in question
		!	capacity:  minimum capacity

		type(LIST)                             :: lst
		integer                                :: capacity
		
		integer                                :: new_cap
		type(LIST_DATA), dimension(:), pointer :: new_data

		new_cap = max(capacity, nint(growth_rate*size(lst%data)))

		if (new_cap .gt. size(lst%data)) then
			allocate(new_data(1:new_cap))
			new_data(1:lst%not_used) = lst%data(1:lst%not_used)
			new_data(lst%not_used + 1:new_cap) = empty_list_data

			deallocate(lst%data)
			lst%data => new_data
		endif
		
	end subroutine list_increase_capacity
	
	!:.........................................................................:

end module vararray
