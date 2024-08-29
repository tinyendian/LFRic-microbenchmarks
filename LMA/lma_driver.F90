!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

program lma_driver
  use constants_mod, only: r_def, i_def, l_def
  use dino_mod, only : dino_type
  use compare_mod, only : compare
  use colourist_mod, only: compute_tiled_colour_map, reorder_data, &
       check_cell_order, check_tile_colouring
  use compute_loop_mod, only: &
       invoke_matrix_vector_kernel_coloured, &
       invoke_matrix_vector_kernel_coloured_tiled
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

  type(dino_type) :: dino

  ! loop counters
  integer(kind=i_def) :: i, j

  integer(kind=i_def) :: count

  ! Tiling
  integer(kind=i_def) :: tile_x = 1, tile_y = 1

  ! Timing
  integer(kind=8) :: startclock, stopclock, clockrate
  real(kind=8) :: timing

  ! CLI
  logical(kind=l_def) :: write_cell_props = .false., reorder_fields = .false.
  logical(kind=l_def) :: tiling = .false.
  logical(kind=l_def) :: colouring = .false.
  character(len=256) :: arg

  ! Read data file - file format is simple sequential ASCII

  write(*,'(A)') 'Loading input data...'

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

  ! Reorder matrix to k-first storage order (Ticket 3811)
  allocate(op_data_transposed, source=reshape(transpose(reshape(op_data, (/ndf1*ndf2, ncell_3d/))), &
                                              (/ncell_3d, ndf1, ndf2/)))

  deallocate(op_data)
  allocate(op_data, source=op_data_transposed)
  deallocate(op_data_transposed)

  write(*,'(A,5(X,A,I10))') "lma_driver: ingested dinodump", "ndf1=", ndf1, "ndf2=", ndf2, &
       "undf1=", undf1, "undf2=", undf2, "nlayers=", nlayers

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
     end select
  end do

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
    call check_tile_colouring(tile_x, tile_y, int(sqrt(ncell/6.0)), ncell, ndf1, map1, ncolours, &
         ntiles_per_colour, tmap)
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

  ! Call the compute loops and measure timings, with or without tiling
  ! Use explicit OpenACC data transfer directives here to ensure that timings only include GPU compute
  if (tiling) then
    !$acc data copyin(ntiles_per_colour, tmap, data1, data2, op_data, map1, map2)
    call system_clock(startclock, clockrate)

    call invoke_matrix_vector_kernel_coloured_tiled(ncolours, ncell_3d, nlayers, ndf1, undf1, ndf2, undf2, &
         & ntiles_per_colour, tile_x*tile_y, tmap, map1, map2, data1, data2, op_data, data1_snapshot)

    call system_clock(stopclock, clockrate)
    !$acc end data
  else
    !$acc data copyin(ncells_per_colour, cmap, data1, data2, op_data, map1, map2)
    call system_clock(startclock, clockrate)

    call invoke_matrix_vector_kernel_coloured(ncolours, ncell_3d, nlayers, ndf1, undf1, ndf2, undf2, &
         & ncells_per_colour, cmap, map1, map2, data1, data2, op_data, data1_snapshot)

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
