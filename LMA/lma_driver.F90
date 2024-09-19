!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

program lma_driver
  use constants_mod, only: r_def, i_def, l_def
  use read_input_file_mod, only: read_dinodump_file, read_binary_file
  use compare_mod, only : compare
  use colourist_mod, only: compute_tiled_colour_map, reorder_data, &
       check_cell_order, check_tile_colouring
  use compute_loop_mod, only: &
       invoke_matrix_vector_kernel, &
       invoke_matrix_vector_kernel_coloured, &
       invoke_matrix_vector_kernel_coloured_tiled, &
       invoke_matrix_vector_kernel_coloured_tiled_vert
  implicit none

  ! mesh sizes
  integer(kind=i_def) :: ncell, ncell_3d, nlayers

  ! colouring sizes and arrays
  integer(kind=i_def)                              :: ncolours, ntiles_per_colour
  integer(kind=i_def), allocatable, dimension(:)   :: ncells_per_colour, cellmap
  integer(kind=i_def), allocatable, dimension(:,:) :: cmap
  integer(kind=i_def), allocatable, dimension(:,:,:) :: tmap

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
  real(kind=r_def), allocatable, dimension(:,:,:) :: op_data, op_data_transposed

  ! Copy of output for comparison with KGO
  real(kind=r_def), allocatable, dimension(:)     :: data1_snapshot

  ! loop counters
  integer(kind=i_def) :: i, j

  integer(kind=i_def) :: count

  ! Tiling
  integer(kind=i_def) :: tile_x = 1, tile_y = 1
  integer(kind=i_def), allocatable :: tile(:), kstart(:), kstop(:)
  integer(kind=i_def) :: ntiles, nksections, nblocks, ksectionlength, ksection

  ! Timing
  integer(kind=8) :: startclock, stopclock, clockrate
  real(kind=8) :: timing

  ! CLI
  logical(kind=l_def) :: write_cell_props = .false.
  logical(kind=l_def) :: reorder_fields = .false.
  logical(kind=l_def) :: tiling = .false.
  logical(kind=l_def) :: colouring = .false.
  logical(kind=l_def) :: binary = .false.
  logical(kind=l_def) :: check = .false.
  logical(kind=l_def) :: inline = .false.
  logical(kind=l_def) :: vertical_tiling = .false.
  integer(kind=i_def) :: ntimes = 1000
  character(len=256) :: arg

  ! CLI
  do i = 1, command_argument_count()
     call get_command_argument(i, arg)
     select case (trim(arg))
        case ('-t', '--tile')
           tiling = .true.
           call get_command_argument(i+1, arg)
           read(arg,*) tile_x
           call get_command_argument(i+2, arg)
           read(arg,*) tile_y
        case ('-w', '--write-cell-props')
           write_cell_props = .true.
        case ('-r', '--reorder')
           reorder_fields = .true.
        case ('-c', '--colouring')
           colouring = .true.
        case ('-b', '--binary')
           binary = .true.
        case ('-C', '--check')
           check = .true.
        case ('-i', '--inline')
           inline = .true.
        case ('-n', '--ntimes')
           call get_command_argument(i+1, arg)
           read(arg,*) ntimes
        case ('-v', '--vertical')
           vertical_tiling = .true.
     end select
  end do
  if (vertical_tiling .and. .not. tiling) then
    write(*, '(A)') 'Vertical tiling requires (horizontal) tiling'
    stop
  end if

  ! Read data file
  if (binary) then
    write(*,'(A)') 'Loading input data from binary file...'

    call read_binary_file(ncell, ncell_3d, ncolours, nlayers, ncells_per_colour, cmap, ndf1, undf1, &
         map1, ndf2, undf2, map2, data1, data2, op_data, answer)

  else
    write(*,'(A)') 'Loading input data from dinodump file...'

    call read_dinodump_file(ncell, ncell_3d, ncolours, nlayers, ncells_per_colour, cmap, ndf1, undf1, &
         map1, ndf2, undf2, map2, data1, data2, op_data, answer)
  end if

  ! Reorder matrix to k-first storage order (Ticket 3811)
  allocate(op_data_transposed, source=reshape(transpose(reshape(op_data, (/ndf1*ndf2, ncell_3d/))), &
                                              (/ncell_3d, ndf1, ndf2/)))

  deallocate(op_data)
  allocate(op_data, source=op_data_transposed)
  deallocate(op_data_transposed)

  write(*,'(A,6(X,A,I10))') "lma_driver: ingested dinodump", "ndf1=", ndf1, "ndf2=", ndf2, &
       "undf1=", undf1, "undf2=", undf2, "nlayers=", nlayers, "ncell=", ncell

  if (colouring) then
    write(*,'(A)') 'Using colouring'
  else
    write(*,'(A)') 'Not using colouring, resetting colour map'

    ! Reset colour map to a single colour, maintaining computation order of cells
    ! to ensure that results are numerically comparable
    allocate(cellmap(ncell))
    j = 1
    do i = 1, ncolours
      cellmap(j:j+ncells_per_colour(i)-1) = cmap(i,1:ncells_per_colour(i))
      j = j + ncells_per_colour(i)
    end do
    deallocate(cmap, ncells_per_colour)

    ncolours = 1
    allocate(ncells_per_colour(1), source=ncell)
    allocate(cmap(1, ncell), source=reshape(cellmap,(/1,ncell/)))
    deallocate(cellmap)
  end if

  ! Compute new tiled colour map
  if (tiling) then
    write(*,'(A)') 'Using tiling'
    call compute_tiled_colour_map(tile_x, tile_y, ncell, tmap, ncolours, ntiles_per_colour, &
         & colouring, write_cell_props)
    if (check) then
      call check_tile_colouring(tile_x, tile_y, int(sqrt(ncell/6.0)), ncell, ndf1, map1, ncolours, &
           ntiles_per_colour, tmap)
    end if
  else
    write(*,'(A)') 'Not using tiling'
  end if

  ! Use cmap to reorder field data (need to transpose to make cell dimension fastest-moving)
  ! map1 - lhs
  ! Also need to reorder answers, otherwise comparison won't work
  ! if (reorder_fields) then
  !    call reorder_data(ncell, ndf1, nlayers, undf1, &
  !         & reshape(transpose(reshape(cmap_tiled, (/ncolours, ntiles_per_colour*tile_x*tile_y/))), (/size(cmap_tiled)/)), &
  !         & map1, data1, answer)
  ! end if
  ! call check_cell_order(tile_x, tile_y, cells_per_dimension, ncolours, ntiles_per_colour, cmap_tiled)

  allocate( data1_snapshot(undf1) )

  ! Call the compute loops and measure timings, with or without tiling
  ! Use explicit OpenACC data transfer directives here to ensure that timings only include GPU compute
  if (vertical_tiling) then

    ! Set up vertical tiling, dividing the vertical axis into the same number of sections as there are cells per tile
    ! to keep the number of total number of work items = tiles*sections constant
    ntiles = ntiles_per_colour/6
    nksections = tile_x*tile_y
    nblocks = ntiles*nksections
    ksectionlength = nlayers/nksections
    if (ksectionlength < 1) then
      write(*,'(A)') 'ERROR - too many vertical sections'
      stop
    end if
    write(*,'(3(A,X,I3,X),A)') 'Vertical tiling: using', nksections, 'vertical sections of >=', ksectionlength, &
         'layers for a total of', nblocks, '3D blocks'

    ! Set up look-up arrays
    allocate(tile(nblocks))
    allocate(kstart(nblocks))
    allocate(kstop(nblocks))

    ! Work out horizontal tile number and vertical k section from block number
    do i = 1, nblocks
      tile(i) = mod(i-1, ntiles)

      ksection = (i-1)/ntiles
      kstart(i) = ksection*ksectionlength
      kstop(i) = kstart(i) + ksectionlength-1
      if (ksection == (nksections-1)) kstop(i) = nlayers-1
    end do

    !$acc data copyin(ntiles_per_colour, tmap, tile, kstart, kstop, data1, data2, op_data, map1, map2)
    call system_clock(startclock, clockrate)

    call invoke_matrix_vector_kernel_coloured_tiled_vert(ncolours,nlayers, ndf1, ndf2, &
         ntiles_per_colour, tile_x*tile_y, tmap, tile, kstart, kstop, nblocks, map1, map2, data1, data2, op_data, data1_snapshot, &
         inline, ntimes)

    call system_clock(stopclock, clockrate)
    !$acc end data

    deallocate(tile, kstart, kstop)

  else if (tiling) then
    !$acc data copyin(ntiles_per_colour, tmap, data1, data2, op_data, map1, map2)
    call system_clock(startclock, clockrate)

    call invoke_matrix_vector_kernel_coloured_tiled(ncolours, ncell_3d, nlayers, ndf1, undf1, ndf2, undf2, &
         ntiles_per_colour, tile_x*tile_y, tmap, map1, map2, data1, data2, op_data, data1_snapshot, &
         inline, ntimes)

    call system_clock(stopclock, clockrate)
    !$acc end data
  else if (colouring) then
    !$acc data copyin(ncells_per_colour, cmap, data1, data2, op_data, map1, map2)
    call system_clock(startclock, clockrate)

    call invoke_matrix_vector_kernel_coloured(ncolours, ncell_3d, nlayers, ndf1, undf1, ndf2, undf2, &
         ncells_per_colour, cmap, map1, map2, data1, data2, op_data, data1_snapshot, inline, ntimes)

    call system_clock(stopclock, clockrate)
    !$acc end data
  else
    !$acc data copyin(data1, data2, op_data, map1, map2)
    call system_clock(startclock, clockrate)

    call invoke_matrix_vector_kernel(ncell_3d, nlayers, ndf1, undf1, ndf2, undf2, &
         map1, map2, data1, data2, op_data, data1_snapshot, inline, ntimes)

    call system_clock(stopclock, clockrate)
    !$acc end data
  end if

  timing = dble(stopclock - startclock)/dble(clockrate)
  write(*,'(A, 1X, F7.2)') 'Timing:', timing

  write(*,'(A)') "lma_driver: Kernel run, checking answer ..."
  count = compare(data1_snapshot, answer, undf1, .true., 1.0d-5)
  write(*,'(A,I8,A,I8,A)') "lma_driver:checked ",undf1," answers, found ",count, " errors" 

  ! deallocate the arrays
  deallocate(ncells_per_colour, cmap)
  deallocate(map1, map2)
  deallocate(data1, data2)
  deallocate(op_data, answer)
  if (allocated(tmap)) deallocate(tmap)

end program lma_driver
