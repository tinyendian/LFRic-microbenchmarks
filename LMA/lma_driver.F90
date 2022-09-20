!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

program lma_driver
  use constants_mod, only: r_def, i_def, l_def
  use matrix_vector_kernel_mod, only: matrix_vector_code
  use dino_mod, only : dino_type
  use compare_mod, only : compare
  implicit none
  
  ! mesh sizes
  integer(kind=i_def) :: ncell, ncell_3d, nlayers
  
  ! colouring sizes and arrays
  integer(kind=i_def)                              :: ncolours
  integer(kind=i_def), allocatable, dimension(:)   :: ncells_per_colour
  integer(kind=i_def), allocatable, dimension(:,:) :: cmap
  
  ! dof-maps for space 1
  integer(kind=i_def) :: ndf1, undf1
  integer(kind=i_def), allocatable, dimension(:,:) :: map1
  
  ! dof-maps for space 2
  integer(kind=i_def) :: ndf2, undf2
  integer(kind=i_def), allocatable, dimension(:,:) :: map2
  
  ! the data
  real(kind=r_def), allocatable, dimension(:)     :: data1
  real(kind=r_def), allocatable, dimension(:)     :: data2
  real(kind=r_def), allocatable, dimension(:)     :: answer
  real(kind=r_def), allocatable, dimension(:,:,:) :: op_data    

  real(kind=r_def), allocatable, dimension(:)     :: data1_snapshot

  type(dino_type) :: dino

  ! loop counters
  integer(kind=i_def) :: colour, cell

  integer(kind=i_def) :: count
  
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
  allocate( data1_snapshot(undf1) )

  ! read the floating point data
  call dino%input_array(op_data, ndf1, ndf2, ncell_3d)  
  call dino%input_array(data1, undf1)
  call dino%input_array(data2, undf2)
  call dino%input_array(answer, undf1)

  call dino%io_close()

  !  write(*,'(A)') "lma_driver:ingested dinodump",ndf1,ndf2
  write(*,*) "lma_driver:ingested dinodump",ndf1,ndf2  

  ! Compute new colour map
  call compute_colour_map(1, 1, 48, cmap, ncells_per_colour, .true.)

  ! call reorder_data(ncell, ndf1, nlayers, undf1, reshape(transpose(cmap), (/ncell/)), map1, data1)
  ! call reorder_data(ncell, ndf2, nlayers, undf2, reshape(transpose(cmap), (/ncell/)), map2, data2)

  ! call check_colours(48, ncell, ndf1, map1, ncolours, maxval(ncells_per_colour), cmap)

  !
  ! Repeat the work 1000 times to hide the cost of reading the data.
  do count =1 , 1000 
  do colour=1, ncolours
  !$omp parallel default(shared), private(cell)
  !$omp do schedule(static)
     do cell=1,ncells_per_colour(colour)
        call matrix_vector_code(cmap(colour,cell), nlayers, data1, data2, &
             ncell_3d, op_data, ndf1, undf1, map1(:,cmap(colour,cell)), &
             ndf2, undf2, map2(:,cmap(colour,cell)) )
     end do
  !$omp end do
  !$omp end parallel
  end do
  if (count.eq.1) data1_snapshot = data1
  end do

  write(*,'(A)') "lma_driver:Kernel run, checking answer ..."
  !check the answer
  count = compare(data1_snapshot, answer, undf1, .false.)
  write(*,'(A,I6,A,I6,A)') "lma_driver:checked ",undf1," answers, found ",count, " errors" 
  
  ! deallocate the arrays
  deallocate(ncells_per_colour, cmap)
  deallocate(map1, map2)
  deallocate(data1, data2)
  deallocate(op_data, answer)

contains

  subroutine compute_colour_map(tile_x, tile_y, ncells_per_dimension, colour_map, ncells_per_colour, write_maps)
    implicit none
    integer(kind=i_def), intent(in) :: tile_x, tile_y, ncells_per_dimension
    integer(kind=i_def), intent(out), allocatable :: colour_map(:,:)
    integer(kind=i_def), intent(out) :: ncells_per_colour(4)
    logical(kind=l_def), intent(in) :: write_maps
    ! Local variables
    integer(kind=i_def) :: cell_colour_order_per_panel(2,2,6)
    integer(kind=i_def) :: cell, ncells, ncells_per_panel, panel_number, tile_number
    integer(kind=i_def) :: tile_id_in_panel, ntiles_per_panel
    integer(kind=i_def) :: ntiles_per_row, ntiles_per_column, row, column
    integer(kind=i_def) :: colour, counter(4)
    integer(kind=i_def) :: cell_properties(6*ncells_per_dimension*ncells_per_dimension,3)

    if (mod(ncells_per_dimension, tile_x) .ne. 0 .or. tile_x .gt. ncells_per_dimension) then
       print *, 'Invalid tile_x'
       stop
    end if
    if (mod(ncells_per_dimension, tile_y) .ne. 0 .or. tile_y .gt. ncells_per_dimension) then
       print *, 'Invalid tile_y'
       stop
    end if

    ! Colour order for 4 colours
    cell_colour_order_per_panel(:,:,1) = reshape((/1,2,3,4/),(/2,2/))
    cell_colour_order_per_panel(:,:,2) = reshape((/3,2,1,4/),(/2,2/))
    cell_colour_order_per_panel(:,:,3) = reshape((/3,4,1,2/),(/2,2/))
    cell_colour_order_per_panel(:,:,4) = reshape((/1,4,3,2/),(/2,2/))
    cell_colour_order_per_panel(:,:,5) = reshape((/3,2,4,1/),(/2,2/))
    cell_colour_order_per_panel(:,:,6) = reshape((/2,3,1,4/),(/2,2/))

    ncells_per_panel = ncells_per_dimension*ncells_per_dimension
    ncells = 6*ncells_per_panel
    ntiles_per_row = ncells_per_dimension/tile_x
    ntiles_per_column = ncells_per_dimension/tile_y
    ntiles_per_panel = ntiles_per_row*ntiles_per_column

    counter = 0
    ! Without domain decomposition, global cell ID = local cell ID
    do cell = 1, ncells

       ! Cubed-sphere panel number, 1...6
       panel_number = (cell - 1) / ncells_per_panel + 1

       if (panel_number .lt. 1 .or. panel_number .gt. 6) then
          print *, 'Invalid panel number', panel_number
          stop
       end if

       ! Compute global tile number for this cell
       tile_number = compute_tile_number(cell, ncells_per_dimension, tile_x, tile_y)

       ! Compute tile number relative to panel (starts from 1 in each panel)
       tile_id_in_panel = tile_number - (ntiles_per_panel * (panel_number-1))

       if (tile_id_in_panel .lt. 1 .or. tile_id_in_panel .gt. ntiles_per_panel) then
          print *, 'Invalid tile_id_in_panel', tile_id_in_panel
          stop
       end if

       ! Compute column and row for this tile
       column = mod(tile_id_in_panel-1, ntiles_per_row) + 1
       row = (tile_id_in_panel-1)/ntiles_per_column + 1

       ! Look up colour
       colour = cell_colour_order_per_panel(2-mod(column,2), mod(row,2)+1, panel_number)

       counter(colour) = counter(colour) + 1
       cell_properties(cell,:) = (/panel_number, tile_number, colour/)

    end do

    if (allocated(colour_map)) then
       deallocate(colour_map)
    end if

    allocate(colour_map(4, maxval(counter)))

    ! Reorder so that cells appear in tile order
    colour_map = -1
    counter = 0
    do tile_number = 1, 6*ntiles_per_panel
       do cell = 1, ncells
          if (cell_properties(cell, 2) .eq. tile_number) then
             colour = cell_properties(cell, 3)
             counter(colour) = counter(colour) + 1
             colour_map(colour, counter(colour)) = cell
          end if
       end do
    end do
    ncells_per_colour = counter

    if (write_maps) then
       print *, 'Writing cell property maps...'
       open(12345, file='lma_driver_maps.bin', form='unformatted', status='unknown', &
            & access='direct', recl=ncells)
       write(12345, rec=1) int(cell_properties(:,1), kind=4)
       write(12345, rec=2) int(cell_properties(:,2), kind=4)
       write(12345, rec=3) int(cell_properties(:,3), kind=4)
       close(12345)

    end if

  end subroutine compute_colour_map

  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  integer(kind=i_def) function compute_tile_number(cellid, ncells_per_dimension, &
       &                                           tile_x, tile_y) result(tile_number)
    implicit none
    integer(kind=i_def), intent(in) :: cellid, ncells_per_dimension, tile_x, tile_y

    integer(kind=i_def) :: panel_number, cellid_in_panel, column_in_panel, row_in_panel
    integer(kind=i_def) :: tile_column, tile_row, ntiles_per_panel, ntiles_per_row, ncells_per_panel

    ncells_per_panel = ncells_per_dimension*ncells_per_dimension

    ! Cubed-sphere panel number, 1...6
    panel_number = (cellid - 1) / ncells_per_panel + 1

    if (panel_number .lt. 1 .or. panel_number .gt. 6) then
       print *, 'Invalid panel number', panel_number
       stop
    end if

    ! Compute ID relative to panel (starts from 1 in each panel)
    cellid_in_panel = cellid - (ncells_per_panel * (panel_number - 1))

    if (cellid_in_panel .lt. 1 .or. cellid_in_panel .gt. ncells_per_panel) then
       print *, 'Invalid cellid in panel', cellid_in_panel
       stop
    end if

    ! Compute column and row IDs for this cell
    column_in_panel = mod(cellid_in_panel-1, ncells_per_dimension) + 1
    row_in_panel    = (cellid_in_panel-1)/ncells_per_dimension + 1

    if (column_in_panel .lt. 1 .or. column_in_panel .gt. ncells_per_dimension) then
       print *, 'Invalid column_in_panel', column_in_panel
       stop
    else if (row_in_panel .lt. 1 .or. row_in_panel .gt. ncells_per_dimension) then
       print *, 'Invalid row_in_panel', row_in_panel
       stop
    end if

    ! Compute column and row numbers of tiles in panel
    tile_column = (column_in_panel-1)/tile_x + 1
    tile_row = (row_in_panel-1)/tile_y + 1

    if (tile_column .lt. 1 .or. tile_column .gt. ncells_per_dimension/tile_x) then
       print *, 'Invalid tile_column', tile_column
       stop
    else if (tile_row .lt. 1 .or. tile_row .gt. ncells_per_dimension/tile_y) then
       print *, 'Invalid tile_row', tile_row
       stop
    end if

    ntiles_per_panel = ncells_per_panel/(tile_x*tile_y)
    ntiles_per_row = ncells_per_dimension/tile_x

    ! Compute tile number for this cell
    tile_number = (panel_number-1)*ntiles_per_panel + &
         & (tile_row-1)*ntiles_per_row + &
         & tile_column

    if (tile_number .lt. 1 .or. tile_number .gt. 6*ntiles_per_panel) then
       print *, 'Invalid tile_number', tile_number
       stop
    end if

  end function compute_tile_number

  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine reorder_data(ncells, ndofs_per_cell, nlayers, ndofs, cell_list, dofmap, data)
    implicit none
    integer(kind=i_def), intent(in) :: ncells, ndofs_per_cell, nlayers, ndofs, cell_list(ncells)
    integer(kind=i_def), intent(inout) :: dofmap(ndofs_per_cell, ncells)
    real(kind=r_def), intent(inout) :: data(ndofs)
    ! Local variables
    integer(kind=i_def) :: i, cell, dof, base_idx, data_idx
    real(kind=r_def) :: buffer(nlayers)
    integer(kind=l_def) :: moved(ndofs/nlayers)

    moved = 0
    base_idx = 1

    ! Walk through cells in the given order and reorder data array
    do i = 1, ncells
       cell = cell_list(i)

       ! Unused cells in list are marked by a negative value
       if (cell > 0) then
          do dof = 1, ndofs_per_cell

             data_idx = dofmap(dof, cell)

             ! Check if data has been moved already
             if ( moved((data_idx-1)/nlayers+1) .eq. 0) then

                if (base_idx .gt. ndofs-nlayers+1) then
                   print *, 'This should not happen'
                end if

                moved((data_idx-1)/nlayers+1) = 1

                ! Swap data
                buffer = data(base_idx:base_idx+nlayers-1)
                data(base_idx:base_idx+nlayers-1) = data(data_idx:data_idx+nlayers-1)
                data(data_idx:data_idx+nlayers-1) = data(base_idx:base_idx+nlayers-1)

                ! Update all pointers to this location
                where (dofmap .eq. data_idx)
                   dofmap = base_idx
                end where

                base_idx = base_idx + nlayers

             end if
          end do
       end if
    end do

  end subroutine reorder_data

  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine check_colours(ncells_per_dimension, ncells, ndofs_per_cell, dofmap, ncolours, ncellmax, colour_map)
    implicit none
    integer(kind=i_def), intent(in) :: ncells_per_dimension, ncells, ndofs_per_cell, ncolours, ncellmax
    integer(kind=i_def), intent(in) :: dofmap(ndofs_per_cell,ncells)
    integer(kind=i_def), intent(in) :: colour_map(ncolours,ncellmax)
    ! Local variables
    integer(kind=i_def) :: cell_colours(ncell), colour, dof, i, j, k, cell, data_idx
    integer(kind=i_def) :: ncells_per_panel, panel_number_1, panel_number_2

    print *, 'Checking colours...'

    ncells_per_panel = ncells_per_dimension*ncells_per_dimension

    cell_colours = -1
    do colour = 1, ncolours
       do i = 1, ncellmax
          cell = colour_map(colour, i)
          if (cell .gt. 0) then
             cell_colours(cell) = colour
          end if
       end do
    end do

    if (minval(cell_colours) .lt. 1) then
       print *, 'Some cells are uncoloured'
       stop
    end if

    ! Loop over all cells
    do cell = 1, ncells

       colour = cell_colours(cell)

       ! Check colours of all horizontal neighbour cells
       do dof = 1, ndofs_per_cell

          ! Identify neighbours by DOF storage location
          data_idx = dofmap(dof, cell)

          do j = 1, ncells
             do k = 1, ndofs_per_cell
                if (dofmap(k,j) .eq. data_idx) then
                   ! Compare colours
                   if (j .ne. cell .and. cell_colours(j) .eq. colour) then

                      ! Cubed-sphere panel number, 1...6
                      panel_number_1 = (cell - 1) / ncells_per_panel + 1
                      panel_number_2 = (j - 1) / ncells_per_panel + 1
                      print *, 'The following neighbour cells have the same colour:', cell, j, ' panels ', panel_number_1, panel_number_2
                   end if
                end if
             end do
          end do
       end do
    end do

    print *, 'Cell colouring is correct.'

  end subroutine check_colours

end program lma_driver
