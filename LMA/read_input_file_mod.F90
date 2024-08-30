module read_input_file_mod
  use constants_mod, only: r_def, i_def
  use dino_mod, only : dino_type

contains

  subroutine read_dinodump_file(ncell, ncell_3d, ncolours, nlayers, ncells_per_colour, cmap, ndf1, undf1, &
                              map1, ndf2, undf2, map2, data1, data2, op_data, answer)
    implicit none

    integer(kind=i_def), intent(out) :: ncell, ncell_3d, ncolours, nlayers, ndf1, undf1, ndf2, undf2
    integer(kind=i_def), allocatable, intent(out) :: ncells_per_colour(:), cmap(:,:), map1(:,:), map2(:,:)
    real(kind=r_def), allocatable, intent(out) :: op_data(:,:,:), data1(:), data2(:), answer(:)

    type(dino_type) :: dino

    ! make the reader
    dino = dino_type()    

    !ingest the data
    call dino%input_scalar(ncell)
    call dino%input_scalar(ncell_3d)
    call dino%input_scalar(ncolours)
    call dino%input_scalar(nlayers)

    ! allocate the colour arrays
    allocate(ncells_per_colour(ncolours))
    call dino%input_array(ncells_per_colour,ncolours)
    allocate( cmap(ncolours,maxval(ncells_per_colour)) )
    call dino%input_array(cmap,ncolours,maxval(ncells_per_colour))

    ! read space sizes and allocate arrays
    call dino%input_scalar(ndf1)
    call dino%input_scalar(undf1)
    allocate(map1(ndf1,ncell))
    call dino%input_array(map1,ndf1, ncell)
  
    call dino%input_scalar(ndf2)
    call dino%input_scalar(undf2)
    allocate(map2(ndf2,ncell))
    call dino%input_array(map2,ndf2, ncell)

    ! allocate the floating point data arrays
    allocate( data1(undf1), data2(undf2) )
    allocate( op_data(ndf1,ndf2,ncell_3d) )
    allocate( answer(undf1) )

    ! read the floating point data
    call dino%input_array(op_data, ndf1, ndf2, ncell_3d)  
    call dino%input_array(data1, undf1)
    call dino%input_array(data2, undf2)
    call dino%input_array(answer, undf1)

    call dino%io_close()

  end subroutine read_dinodump_file

  subroutine read_binary_file(ncell, ncell_3d, ncolours, nlayers, ncells_per_colour, cmap, ndf1, undf1, &
                              map1, ndf2, undf2, map2, data1, data2, op_data, answer)
    implicit none

    integer(kind=i_def), intent(out) :: ncell, ncell_3d, ncolours, nlayers, ndf1, undf1, ndf2, undf2
    integer(kind=i_def), allocatable, intent(out) :: ncells_per_colour(:), cmap(:,:), map1(:,:), map2(:,:)
    real(kind=r_def), allocatable, intent(out) :: op_data(:,:,:), data1(:), data2(:), answer(:)

    integer :: funit

    funit = 888

    open(unit=funit, file='dinodump_binary.dat', status='old', action='read', form='unformatted')
    read(funit) ncell
    read(funit) ncell_3d
    read(funit) ncolours
    read(funit) nlayers

    allocate(ncells_per_colour(ncolours))
    read(funit) ncells_per_colour(:)

    allocate(cmap(ncolours,maxval(ncells_per_colour)))
    read(funit) cmap(:,:)

    read(funit) ndf1
    read(funit) undf1

    allocate(map1(ndf1,ncell))
    read(funit) map1(:,:)
    read(funit) ndf2
    read(funit) undf2

    allocate(map2(ndf2,ncell))
    read(funit) map2(:,:)

    allocate(op_data(ndf1,ndf2,ncell_3d))
    read(funit) op_data(:,:,:)

    allocate(data1(undf1))
    read(funit) data1(:)

    allocate(data2(undf2))
    read(funit) data2(:)

    allocate(answer(undf1))
    read(funit) answer(:)

    close(funit)

  end subroutine read_binary_file

end module read_input_file_mod
