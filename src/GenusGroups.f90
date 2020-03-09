module Genusgroups
  use Constants
  use Site
  use vararray
  implicit none

!*******************************************************************************
  !
  !This module gets the genus and species names of the input Species file.
  !
!*******************************************************************************


	type Groups
		character(len = MAX_NLEN), allocatable, dimension(:)   :: genusgroups
		character(len = MAX_NLEN), allocatable, dimension(:,:) :: spec_names
		integer                                                :: numgenera
		integer                                                :: numspecies
	end type Groups

contains

	!:.........................................................................:

	subroutine initialize_genus_groups(group, species_data)
		!initializes the list of genus and species names
		!Author: Katherine Holcomb, 2012 v. 1.0
		!Inputs:
		!   species_data:  species type
		!Outputs:
		!   species_data:  species type
		!   group:         group type

		class(Groups),                   intent(out)         :: group
		type(SpeciesData), dimension(:), intent(inout)       :: species_data

		character(len = MAX_NLEN), dimension(:), allocatable :: species_names
		integer                                              :: num_all_species
		integer                                              :: ns
		integer                                              :: na

        !get how many species there are at site
		num_all_species = size(species_data)

		call get_unique_items(species_data(:)%unique_id, species_names)
		call get_unique_items(species_data(:)%genus_name, group%genusgroups)

		group%numspecies = size(species_names)
		group%numgenera = size(group%genusgroups)

		allocate(group%spec_names(group%numspecies, 2))

		do na = 1, num_all_species
			do ns = 1, group%numspecies
				if (species_names(ns) .eq.                                     &
								species_data(na)%unique_id) then
					group%spec_names(ns, 1) =                                  &
										species_data(na)%genus_name
					group%spec_names(ns, 2) =                                  &
										species_data(na)%unique_id
					cycle
				endif
			enddo
		enddo

	end subroutine initialize_genus_groups

	!:.........................................................................:

	subroutine get_unique_items(array, unique_array)
		!gets unique items from an array
		!Author: Katherine Holcomb, 2012 v. 1.0
		!Inputs/Outputs:
		!   array:         input array
		!   unique_array:  unique items from input array

		character(len = *),     dimension(:), intent(inout)                 :: array
		character(len = MAX_NLEN), dimension(:), allocatable, intent(out)   :: unique_array
		character(len = MAX_NLEN), dimension(size(array))                   :: temp_array
		integer                                                             :: n, ncount
		integer                                                             :: nsize

		nsize = size(array)
        temp_array = array

		!if the array has only one element then we are done
		if (nsize .eq. 1) then
			allocate(unique_array(1))
			unique_array(1) = array(1)
			return
		endif

		!sort the array
		call sort(temp_array, nsize)

		!first count the number of unique items (only required since we want the
			!array returned to be the correct length)
		ncount = 1
		do n = 2, nsize
			if (temp_array(n - 1) .ne. temp_array(n)) then
				ncount = ncount + 1
			endif
		enddo

		allocate(unique_array(ncount))

		unique_array(1) = temp_array(1)
		ncount = 2
		do n = 2, nsize
			if (temp_array(n - 1) .ne. temp_array(n)) then
				unique_array(ncount) = trim(temp_array(n))
				ncount = ncount + 1
			endif
		enddo

	end subroutine get_unique_items

	!:.........................................................................:

end module GenusGroups
