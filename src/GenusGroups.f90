module Genusgroups

!*******************************************************************************
  !
  ! This defines a type to hold the genus and species names of SpeciesObject
  ! arrays.
  !
!*******************************************************************************

  use Constants
  use Site
  implicit none

    ! Define Groups type
    type Groups
        character(len = MAX_NLEN), allocatable, dimension(:)   :: genusgroups ! Array of genus names
        character(len = MAX_NLEN), allocatable, dimension(:,:) :: spec_names  ! Array of genus x species names
        integer                                                :: numgenera   ! Number of genera
        integer                                                :: numspecies  ! Number of species
    end type Groups

contains

    !:.........................................................................:

    subroutine initialize_genus_groups(group, species_data)
        !
        !  Initializes a list of genus and species names from an input array
        !   of Species objects
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        class(Groups),                   intent(out)   :: group        ! Group object
        type(SpeciesData), dimension(:), intent(inout) :: species_data ! Species opjects array

        ! Data dictonary: calling arguments
        character(len = MAX_NLEN), dimension(:), allocatable :: species_names   ! Array of species names
        integer                                              :: num_all_species ! Number of species
        integer                                              :: ns, na          ! Looping indices

        ! Get how many species there are at site
        num_all_species = size(species_data)

        ! Get unique arrays of the unique id and the genera
        call get_unique_items(species_data(:)%unique_id, species_names)
        call get_unique_items(species_data(:)%genus_name, group%genusgroups)

        ! Set number of species and genera
        group%numspecies = size(species_names)
        group%numgenera = size(group%genusgroups)

        ! Allocate size of spec_names array
        allocate(group%spec_names(group%numspecies, 2))

        ! Loop through species names and genera and add to instance array
        do na = 1, num_all_species
            do ns = 1, group%numspecies
                if (species_names(ns) ==                                       &
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
        !
        !  Gets the unique character items from an array (sorted)
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments:
        character(len = *),        dimension(:),              intent(inout) :: array        ! Input array
        character(len = MAX_NLEN), dimension(:), allocatable, intent(out)   :: unique_array ! Unique items from array

        ! Data dictionary: local variables
        character(len = MAX_NLEN), dimension(size(array)) :: temp_array ! Temporary array
        integer                                           :: n          ! Looping index
        integer                                           :: ncount     ! Number of unique items
        integer                                           :: nsize      ! Size of input array

        nsize = size(array)
        temp_array = array

        ! If the array has only one element then we are done
        if (nsize == 1) then
            allocate(unique_array(1))
            unique_array(1) = array(1)
            return
        endif

        ! Sort the array
        call sort(temp_array, nsize)

        ! First count the number of unique items (only required since we want the
        ! array returned to be the correct length)
        ncount = 1
        do n = 2, nsize
            if (temp_array(n - 1) /= temp_array(n)) then
                ncount = ncount + 1
            endif
        enddo

        ! Allocate output array
        allocate(unique_array(ncount))

        ! Loop through and fill with values, skipping duplicates
        unique_array(1) = temp_array(1)
        ncount = 2
        do n = 2, nsize
            if (temp_array(n - 1) /= temp_array(n)) then
                unique_array(ncount) = trim(temp_array(n))
                ncount = ncount + 1
            endif
        enddo

    end subroutine get_unique_items

    !:.........................................................................:

end module GenusGroups
