module colourist_mod
  use constants_mod, only: i_def, r_def, l_def
  implicit none

contains

  subroutine compute_tiled_colour_map(tile_x, tile_y, ncells_per_dimension, colour_map, ncolours, &
       & ntiles_per_colour, colouring, write_maps)
    implicit none
    integer(kind=i_def), intent(in) :: tile_x, tile_y, ncells_per_dimension
    integer(kind=i_def), intent(out), allocatable :: colour_map(:,:,:)
    integer(kind=i_def), intent(out) :: ncolours, ntiles_per_colour
    logical(kind=l_def), intent(in) :: colouring, write_maps
    ! Local variables
    integer(kind=i_def) :: cell, ncells, ncells_per_panel, panel_number, tile_number, tile
    integer(kind=i_def) :: tile_id_in_panel, ntiles_per_panel
    integer(kind=i_def) :: ntiles_per_row, ntiles_per_column, row, column, colour
    integer(kind=i_def), allocatable :: counter(:)
    integer(kind=i_def) :: cell_properties(6*ncells_per_dimension*ncells_per_dimension,3)

    ! Panel dimensions must be integer multiples of tile dimension
    if (mod(ncells_per_dimension, tile_x) .ne. 0 .or. tile_x .gt. ncells_per_dimension) then
       print *, 'Invalid tile_x'
       stop
    end if
    if (mod(ncells_per_dimension, tile_y) .ne. 0 .or. tile_y .gt. ncells_per_dimension) then
       print *, 'Invalid tile_y'
       stop
    end if

    ! Compute tiling configuration
    ncells_per_panel = ncells_per_dimension*ncells_per_dimension
    ncells = 6*ncells_per_panel
    ntiles_per_row = ncells_per_dimension/tile_x
    ntiles_per_column = ncells_per_dimension/tile_y
    ntiles_per_panel = ntiles_per_row*ntiles_per_column

    ! Set the number of colours according to what is required for the requested tiling configuration
    ! Colour entire panels if no colouring is requested, to prevent full parallelisation over the entire
    ! mesh, which could lead to an unfair comparison with the other setups
    if ((.not. colouring) .or. (ntiles_per_row .eq. 1 .and. ntiles_per_column .eq. 1)) then
       ncolours = 6
    else if (ntiles_per_row .eq. 1 .or. ntiles_per_column .eq. 1) then
       ncolours = 12
    else
       ncolours = 24
    end if
    allocate(counter(ncolours))

    write(*,'(A)') 'Tiling configuration:'
    write(*,'(A,X,I7)') 'Cells:', ncells
    write(*,'(A,X,I7)') 'Cells per panel:', ncells_per_panel
    write(*,'(A,2(X,I4))') 'Tile size:', tile_x, tile_y
    write(*,'(A,X,I7)') 'Tiles per panel:', ntiles_per_panel
    write(*,'(A,X,I7)') 'Tiles per panel row:', ntiles_per_row
    write(*,'(A,X,I7)') 'Tiles per panel column:', ntiles_per_column
    write(*,'(A,X,I7)') 'Colours:', ncolours
    if (colouring) then
       write(*,'(A)') 'Running with colouring'
    else
       write(*,'(A)') 'Running WITHOUT colouring (single colour per panel)'
    end if
    write(*,'(A)') 'Computing tiled colour map...'

    ! Loop over cells and assign tile numbers to each cell and colours to each tile
    counter = 0
    do cell = 1, ncells

       ! Cubed-sphere panel number, 1...6
       panel_number = (cell - 1) / ncells_per_panel + 1

       if (panel_number .lt. 1 .or. panel_number .gt. 6) then
          write(*,'(A,X,I3)') 'Invalid panel number', panel_number
          stop
       end if

       ! Compute global tile number for this cell
       tile_number = compute_tile_number(cell, ncells_per_dimension, tile_x, tile_y)

       ! Compute tile number relative to panel (starts from 1 in each panel)
       tile_id_in_panel = tile_number - (ntiles_per_panel * (panel_number-1))

       if (tile_id_in_panel .lt. 1 .or. tile_id_in_panel .gt. ntiles_per_panel) then
          write(*,'(A,X,I9)') 'Invalid tile_id_in_panel', tile_id_in_panel
          stop
       end if

       ! Compute column and row for this tile
       column = mod(tile_id_in_panel-1, ntiles_per_row) + 1
       row = (tile_id_in_panel-1)/ntiles_per_row + 1

       ! Assign tile colour, each panel uses a different set of colours to avoid race conditions
       if ((.not. colouring) .or. (ntiles_per_row .eq. 1 .and. ntiles_per_column .eq. 1)) then
          colour = panel_number ! No tiling or no (in-panel) colouring - single colour per panel
       else if (ntiles_per_row .eq. 1) then
          colour = (panel_number-1)*2 + 1 + 1-mod(row,2) ! Row tiling - two colours per panel
       else if (ntiles_per_column .eq. 1) then
          colour = (panel_number-1)*2 + 1 + 1-mod(column,2) ! Column tiling - two colours per panel
       else
          colour = (panel_number-1)*4 + 1 + 1-mod(column,2)+ 2*(1-mod(row,2)) ! Tiling - four colours per panel
       end if

       ! Sanity check
       if (colour .lt. 1 .or. colour .gt. ncolours) then
          write(*,'(A)') 'Invalid colour'
          stop
       end if

       counter(colour) = counter(colour) + 1
       cell_properties(cell,:) = (/panel_number, tile_number, colour/)

    end do

    ! Sanity checks
    if (minval(counter) .ne. maxval(counter)) then
       write(*,'(A)') 'Unequal number of cells per colour'
       stop
    end if
    if (sum(counter) .ne. ncells) then
       write(*,'(A)') 'Unexpected number of cells'
       stop
    end if

    ! Number of tiles of the same colour
    ntiles_per_colour = counter(1)/(tile_x*tile_y)

    if (allocated(colour_map)) deallocate(colour_map)
    allocate(colour_map(ncolours, ntiles_per_colour, tile_x*tile_y))

    ! Loop over tiles and arrange cells contiguously per tile to improve cache utilisation
    colour_map = -1
    counter = 0
    do tile_number = 1, 6*ntiles_per_panel
       do cell = 1, ncells
          if (cell_properties(cell, 2) .eq. tile_number) then
             colour = cell_properties(cell, 3)
             ! Fill up tiles of the same colour
             tile = counter(colour)/(tile_x*tile_y) + 1
             counter(colour) = counter(colour) + 1
             colour_map(colour, tile, counter(colour) - (tile-1)*tile_x*tile_y) = cell
          end if
       end do
    end do

    ! Sanity check
    if (minval(colour_map) .lt. 1) then
       write(*,'(A)') 'Colour map incomplete'
       stop
    end if

    if (write_maps) then
       write(*,'(A)') 'Writing cell property maps...'
       open(12345, file='lma_driver_maps.bin', form='unformatted', status='unknown', &
            & access='direct', recl=ncells)
       write(12345, rec=1) int(cell_properties(:,1), kind=4)
       write(12345, rec=2) int(cell_properties(:,2), kind=4)
       write(12345, rec=3) int(cell_properties(:,3), kind=4)
       close(12345)

    end if

  end subroutine compute_tiled_colour_map

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

  subroutine reorder_data(ncells, ndofs_per_cell, nlayers, ndofs_total, cell_list, dofmap, data, data2)
    implicit none
    integer(kind=i_def), intent(in) :: ncells, ndofs_per_cell, nlayers, ndofs_total, cell_list(:)
    integer(kind=i_def), intent(inout) :: dofmap(ndofs_per_cell, ncells)
    real(kind=r_def), intent(inout) :: data(ndofs_total)
    real(kind=r_def), intent(inout), optional :: data2(ndofs_total)
    ! Local variables
    integer(kind=i_def) :: i, cell, dof, base_idx, data_idx, j, k
    real(kind=r_def) :: reordered_data(ndofs_total), reordered_data2(ndofs_total)
    integer(kind=i_def) :: moved(ndofs_total), col_length, reordered_dofmap(ndofs_per_cell, ncells)

    write(*,'(A)') 'Reordering data...'

    ! W2 function spaces:
    ! - 6 DOFs per cell, 4 horizontal DOFs and 2 vertical ones
    ! - Horizontal DOFs live in individual columns of nlayers length
    ! - Vertical DOFs share a column of nlayers+1 length

    if (ndofs_per_cell .ne. 6) then
       write(*,'(A)') 'reorder_data currently only works for W2 fields'
       stop
    end if

    moved = 0
    base_idx = 1
    reordered_dofmap = -1

    ! Walk through cells in the given order and reorder data array
    do i = 1, size(cell_list)
       cell = cell_list(i)

       ! Skip the last DOF, it shares the same column as the second-to-last DOF
       do dof = 1, ndofs_per_cell-1

          data_idx = dofmap(dof, cell)

          ! Check if data has been moved already (DOFs are shared)
          if ( moved(data_idx) .eq. 0 ) then

             moved(data_idx) = 1

             ! The first 4 dofs are horizontal ones where column length is nlayers,
             ! the 5th dof is vertical, where column length is nlayers+1
             if ( dof .le. 4 ) then
                col_length = nlayers
             else
                col_length = nlayers+1
             end if

             reordered_data(base_idx:base_idx+col_length-1) = data(data_idx:data_idx+col_length-1)

             if (present(data2)) then
                reordered_data2(base_idx:base_idx+col_length-1) = data2(data_idx:data_idx+col_length-1)
             end if

             ! Update all pointers to this location
             do j = 1, ncells
                do k = 1, ndofs_per_cell-1
                   if (dofmap(k, j) .eq. data_idx) then
                      reordered_dofmap(k, j) = base_idx
                      ! Also update 6th DOF
                      if (k .eq. ndofs_per_cell-1) reordered_dofmap(k+1, j) = base_idx+1
                   end if
                end do
             end do

             base_idx = base_idx + col_length

          end if
       end do

    end do
    if ( base_idx .ne. (ndofs_total+1) ) then
       print *, 'Reordered array is not complete!'
       stop
    end if
    if ( minval(reordered_dofmap) .eq. -1 ) then
       print *, 'Reordered dofmap is not complete!'
       stop
    end if

    dofmap(:,:) = reordered_dofmap(:,:)
    data(:) = reordered_data(:)

    if (present(data2)) then
       data2(:) = reordered_data2(:)
    end if

  end subroutine reorder_data

  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine check_tile_colouring(tile_x, tile_y, ncells_per_dimension, ncells, ndofs_per_cell, dofmap, ncolours, &
       &                          ntiles_per_colour, colour_map)
    implicit none
    integer(kind=i_def), intent(in) :: tile_x, tile_y
    integer(kind=i_def), intent(in) :: ncells_per_dimension, ncells, ndofs_per_cell, ncolours, ntiles_per_colour
    integer(kind=i_def), intent(in) :: dofmap(ndofs_per_cell,ncells)
    integer(kind=i_def), intent(in) :: colour_map(ncolours, ntiles_per_colour, tile_x*tile_y)
    ! Local variables
    integer(kind=i_def) :: colour, cell, dof, nb_cell, nb_dof, cell_in_tile, tile
    integer(kind=i_def) :: ncells_per_panel, ncells_per_tile, panel_number_1, panel_number_2, num_neighbours
    integer(kind=i_def) :: cell_colours(ncells), tile_number(ncells)

    write(*, '(A)') 'Checking tile colouring...'

    ncells_per_panel = ncells_per_dimension*ncells_per_dimension
    ncells_per_tile = tile_x*tile_y

    ! Look up colour for each cell
    cell_colours = -1
    do colour = 1, ncolours
       do tile = 1, ntiles_per_colour
          do cell_in_tile = 1, ncells_per_tile
             cell_colours(colour_map(colour, tile, cell_in_tile)) = colour
          end do
       end do
    end do

    if (minval(cell_colours) .lt. 1) then
       write(*, '(A)') 'Some cells are uncoloured'
       stop
    end if

    ! Compute tile number for each cell
    do cell = 1, ncells
       tile_number(cell) = compute_tile_number(cell, ncells_per_dimension, tile_x, tile_y)
    end do

    ! Loop over all DOFs in all cells, identify neighbours with shared DOFs, and compare colours
    ! cells are not in the same tile
    num_neighbours = 0
    !$omp parallel do default(shared) private(cell, nb_cell, dof, nb_dof, panel_number_1, panel_number_2) &
    !$omp & reduction(+:num_neighbours)
    do cell = 1, ncells
       do nb_cell = 1, ncells

          if (nb_cell .ne. cell) then

             ! If cells belong to the same tile, their colours must be the same
             if (tile_number(nb_cell) .eq. tile_number(cell) .and. cell_colours(nb_cell) .ne. cell_colours(cell)) then
                write(*,'(A)') 'Inconsistent colours in the same tile'
                stop
             end if

             do dof = 1, ndofs_per_cell
                do nb_dof = 1, ndofs_per_cell

                   ! If cells share the same DOF, they are neighbours
                   if ( dofmap(nb_dof, nb_cell) .eq. dofmap(dof, cell)) then

                      num_neighbours = num_neighbours + 1

                      ! Flag neighbour cells that have the same colour, but don't belong to the same tile
                      if (tile_number(nb_cell) .ne. tile_number(cell) .and. &
                           & cell_colours(nb_cell) .eq. cell_colours(cell)) then
                         ! Cubed-sphere panel number, 1...6
                         panel_number_1 = (cell - 1) / ncells_per_panel + 1
                         panel_number_2 = (nb_cell - 1) / ncells_per_panel + 1
                         write(*, '(2(A,X,I6,X,I6,X))') 'Neighbour cells with same colour:', cell, nb_cell, &
                              & ' panels ', panel_number_1, panel_number_2
                      end if

                   end if
                end do
             end do

          end if

       end do
    end do
    !$omp end parallel do

    write(*, '(A,X,F4.1,X,A)') 'Cells share ', real(num_neighbours)/real(ncells), ' DOFs on average'
    write(*,'(A)') 'Cell colouring check done.'

  end subroutine check_tile_colouring

  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine check_cell_order(tile_x, tile_y, ncells_per_dimension, ncolours, ntiles_per_colour, colour_map)
    implicit none
    integer(kind=i_def), intent(in) :: tile_x, tile_y, ncells_per_dimension, ncolours, ntiles_per_colour
    integer(kind=i_def), intent(in) :: colour_map(ncolours, ntiles_per_colour, tile_x*tile_y)
    ! Local variables
    integer(kind=i_def) :: colour, cell, previous_cell, cell_in_tile, tile_number, tile
    integer(kind=i_def) :: ncells_per_tile

    write(*,'(A)') 'Checking cell order in colour map...'

    ncells_per_tile = tile_x*tile_y

    ! Loop over all cells in the colour map, check that are all valid or invalid, that they belong to the same tile,
    ! and that they are ordered contiguously within tile rows
    do colour = 1, ncolours
       do tile = 1, ntiles_per_colour
          tile_number = 0
          previous_cell = 0
          do cell_in_tile = 1, ncells_per_tile

             cell = colour_map(colour, tile, cell_in_tile)

             ! First cell in a new tile - reset tile number, check for consistency otherwise
             if (cell_in_tile .eq. 1) then
                tile_number = compute_tile_number(cell, ncells_per_dimension, tile_x, tile_y)
             else if (compute_tile_number(cell, ncells_per_dimension, tile_x, tile_y) .ne. tile_number) then
                write(*, '(A, 3(X,I6))') 'Cell does not belong to current tile - cell, tile number, current tile:', &
                     & cell, compute_tile_number(cell, ncells_per_dimension, tile_x, tile_y), tile_number
             end if

             ! Check if cells are contiguous within tile rows
             if (mod(cell_in_tile-1, tile_x) .ne. 0 .and. cell .ne. previous_cell + 1) then
                write(*, '(A, 2(X,I6))') 'Cell order is not continuous - cell, previous cell:', &
                     & cell, previous_cell
             end if

             previous_cell = cell

          end do
       end do
    end do

    write(*,'(A)') 'Cell order check done.'

  end subroutine check_cell_order

end module colourist_mod
