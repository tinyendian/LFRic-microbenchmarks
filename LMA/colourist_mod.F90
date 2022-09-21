module colourist_mod
  use constants_mod, only: i_def, r_def, l_def
  implicit none

contains

  subroutine compute_colour_map(tile_x, tile_y, ncells_per_dimension, colour_map, ncells_per_colour, write_maps)
    implicit none
    integer(kind=i_def), intent(in) :: tile_x, tile_y, ncells_per_dimension
    integer(kind=i_def), intent(out), allocatable :: colour_map(:,:)
    integer(kind=i_def), intent(out) :: ncells_per_colour(24)
    logical(kind=l_def), intent(in) :: write_maps
    ! Local variables
    integer(kind=i_def) :: cell_colour_order_per_panel(2,2,6)
    integer(kind=i_def) :: cell, ncells, ncells_per_panel, panel_number, tile_number
    integer(kind=i_def) :: tile_id_in_panel, ntiles_per_panel
    integer(kind=i_def) :: ntiles_per_row, ntiles_per_column, row, column
    integer(kind=i_def) :: colour, counter(24)
    integer(kind=i_def) :: cell_properties(6*ncells_per_dimension*ncells_per_dimension,3)

    if (mod(ncells_per_dimension, tile_x) .ne. 0 .or. tile_x .gt. ncells_per_dimension) then
       print *, 'Invalid tile_x'
       stop
    end if
    if (mod(ncells_per_dimension, tile_y) .ne. 0 .or. tile_y .gt. ncells_per_dimension) then
       print *, 'Invalid tile_y'
       stop
    end if

    ! Colour order for 4 colours - modified so that each panel uses a different set of colours
    cell_colour_order_per_panel(:,:,1) = reshape((/1,2,3,4/),(/2,2/))
    cell_colour_order_per_panel(:,:,2) = reshape((/3,2,1,4/),(/2,2/)) + 4
    cell_colour_order_per_panel(:,:,3) = reshape((/3,4,1,2/),(/2,2/)) + 8
    cell_colour_order_per_panel(:,:,4) = reshape((/1,4,3,2/),(/2,2/)) + 12
    cell_colour_order_per_panel(:,:,5) = reshape((/3,2,4,1/),(/2,2/)) + 16
    cell_colour_order_per_panel(:,:,6) = reshape((/2,3,1,4/),(/2,2/)) + 20

    ncells_per_panel = ncells_per_dimension*ncells_per_dimension
    ncells = 6*ncells_per_panel
    ntiles_per_row = ncells_per_dimension/tile_x
    ntiles_per_column = ncells_per_dimension/tile_y
    ntiles_per_panel = ntiles_per_row*ntiles_per_column

    counter = 0
    ! Colouring across cubed-sphere panel edges is not yet correct for non-square tiles
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
       row = (tile_id_in_panel-1)/ntiles_per_row + 1

       ! Look up colour
       colour = cell_colour_order_per_panel(2-mod(column,2), mod(row,2)+1, panel_number)

       counter(colour) = counter(colour) + 1
       cell_properties(cell,:) = (/panel_number, tile_number, colour/)

    end do

    write(*, '(A)') 'Number of cells per colour:'
    do colour = 1, 24
       write(*, '(A,X,I,X,A,X,I)') 'Colour:', colour, 'Cells:', counter(colour)
    end do

    if (allocated(colour_map)) then
       deallocate(colour_map)
    end if

    allocate(colour_map(24, maxval(counter)))

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
       write(*,'(A)') 'Writing cell property maps...'
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

    write(*,'(A)') 'Reordering data...'

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

  subroutine check_colours(tile_x, tile_y, ncells_per_dimension, ncells, ndofs_per_cell, dofmap, ncolours, ncellmax, colour_map)
    implicit none
    integer(kind=i_def), intent(in) :: tile_x, tile_y
    integer(kind=i_def), intent(in) :: ncells_per_dimension, ncells, ndofs_per_cell, ncolours, ncellmax
    integer(kind=i_def), intent(in) :: dofmap(ndofs_per_cell,ncells)
    integer(kind=i_def), intent(in) :: colour_map(ncolours, ncellmax)
    ! Local variables
    integer(kind=i_def) :: colour, cell, dof, nb_cell, nb_dof, i
    integer(kind=i_def) :: ncells_per_panel, panel_number_1, panel_number_2, num_neighbours
    integer(kind=i_def) :: cell_colours(ncells), tile_number(ncells)

    write(*, '(A)') 'Checking colours...'

    ncells_per_panel = ncells_per_dimension*ncells_per_dimension

    ! Reconstruct colour table for each cell
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
       write(*, '(A)') 'Some cells are uncoloured'
       stop
    end if

    ! Compute tile numbers
    do cell = 1, ncells
       tile_number(cell) = compute_tile_number(cell, ncells_per_dimension, tile_x, tile_y)
    end do

    ! Loop over all DOFs in all cells, identify neighbours with shared DOFs, and compare colours
    ! cells are not in the same tile
    num_neighbours = 0
    do cell = 1, ncells
       do nb_cell = 1, ncells

          if (nb_cell .ne. cell) then

             do dof = 1, ndofs_per_cell
                do nb_dof = 1, ndofs_per_cell

                   if ( dofmap(nb_dof, nb_cell) .eq. dofmap(dof, cell)) then

                      num_neighbours = num_neighbours + 1
                      if (tile_number(nb_cell) .ne. tile_number(cell) .and. &
                           & cell_colours(nb_cell) .eq. cell_colours(cell)) then
                         ! Cubed-sphere panel number, 1...6
                         panel_number_1 = (cell - 1) / ncells_per_panel + 1
                         panel_number_2 = (nb_cell - 1) / ncells_per_panel + 1
                         write(*, '(2(A,X,I,X,I,X))') 'Neighbour cells with same colour:', cell, nb_cell, &
                              & ' panels ', panel_number_1, panel_number_2
                      end if

                   end if
                end do
             end do

          end if

       end do
    end do

    write(*, '(A,X,F4.1,X,A)') 'Cells share ', real(num_neighbours)/real(ncells), ' DOFs on average'
    write(*,'(A)') 'Cell colouring check done.'

  end subroutine check_colours

  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine check_cell_order(tile_x, tile_y, ncells_per_dimension, ncolours, ncellmax, colour_map, ncells_per_colour)
    implicit none
    integer(kind=i_def), intent(in) :: tile_x, tile_y, ncells_per_dimension, ncolours, ncellmax
    integer(kind=i_def), intent(in) :: colour_map(ncolours, ncellmax), ncells_per_colour(ncolours)
    ! Local variables
    integer(kind=i_def) :: colour, cell, previous_cell, i, tile_number
    integer(kind=i_def) :: ncells_per_tile

    write(*,'(A)') 'Checking cell order in colour map...'

    ncells_per_tile = tile_x*tile_y

    ! Loop over all cells in the colour map and check that they belong to the same tile,
    ! and that they are ordered continuously within tile rows
    do colour = 1, ncolours
       previous_cell = 0
       tile_number = 0
       do i = 1, ncells_per_colour(colour)

          cell = colour_map(colour, i)

          ! First cell in a new tile - reset tile number
          if (mod(i-1, ncells_per_tile) .eq. 0) then
             tile_number = compute_tile_number(cell, ncells_per_dimension, tile_x, tile_y)
          else if (compute_tile_number(cell, ncells_per_dimension, tile_x, tile_y) .ne. tile_number ) then
             write(*, '(A, 3(X,I))') 'Cell does not belong to current tile - cell, tile number, current tile:', &
                  & cell, compute_tile_number(cell, ncells_per_dimension, tile_x, tile_y), tile_number
          end if

          if (mod(i-1, tile_x) .ne. 0 .and. cell .ne. previous_cell + 1) then
             write(*, '(A, 2(X,I))') 'Cell order is not continous - cell, previous cell:', &
                  & cell, previous_cell
          end if

          previous_cell = cell

       end do
    end do

    write(*,'(A)') 'Cell order check done.'

  end subroutine check_cell_order

end module colourist_mod
