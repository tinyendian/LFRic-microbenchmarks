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
       & check_cell_order, check_tile_colouring
  use compute_loop_mod, only: compute_loop_tiled, compute_loop_inlined
  implicit none
  
  ! mesh sizes
  integer(kind=i_def) :: ncell, ncell_3d, nlayers
  
  ! colouring sizes and arrays
  integer(kind=i_def)                              :: ncolours, ntiles_per_colour
  integer(kind=i_def), allocatable, dimension(:)   :: ncells_per_colour
  integer(kind=i_def), allocatable, dimension(:,:) :: cmap
  integer(kind=i_def), allocatable, dimension(:,:,:) :: cmap_tiled
  
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
  integer(kind=i_def) :: i

  integer(kind=i_def) :: count

  ! Tiles
  integer(kind=i_def) :: cells_per_dimension = 48, tile_x = 1, tile_y = 1

  ! Timing
  integer(kind=8) :: startclock, stopclock, clockrate
  real(kind=8) :: timing

  ! CLI
  logical(kind=l_def) :: write_cell_props = .false., reorder_fields = .false.
  character(len=256) :: arg

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
  write(*,*) "lma_driver:ingested dinodump",ndf1,ndf2,undf1,undf2,nlayers

  do i = 1, command_argument_count()
     call get_command_argument(i, arg)
     select case (trim(arg))
        case ('-t', '--tile')
           call get_command_argument(i+1, arg)
           read(arg,*) tile_x
           call get_command_argument(i+2, arg)
           read(arg,*) tile_y
        case ('-w', '--write-cell-props')
           write_cell_props = .true.
        case ('-r', '--reorder')
           reorder_fields = .true.
     end select
  end do

  ! Compute new tiled colour map
  call compute_tiled_colour_map(tile_x, tile_y, cells_per_dimension, cmap_tiled, ncolours, ntiles_per_colour, write_cell_props)

  ! Use cmap to reorder field data (need to transpose to make cell dimension fastest-moving)
  ! map1 - lhs
  ! Also need to reorder answers, otherwise comparison won't work
  if (reorder_fields) then
     call reorder_data(ncell, ndf1, nlayers, undf1, &
          & reshape(transpose(reshape(cmap_tiled, (/ncolours, ntiles_per_colour*tile_x*tile_y/))), (/size(cmap_tiled)/)), &
          & map1, data1, answer)
  end if

  call check_tile_colouring(tile_x, tile_y, cells_per_dimension, ncell, ndf1, map1, ncolours, ntiles_per_colour, cmap_tiled)
  call check_cell_order(tile_x, tile_y, cells_per_dimension, ncolours, ntiles_per_colour, cmap_tiled)

  !$acc data copyin(ncells_per_colour, cmap, data1, data2, op_data, map1, map2)
  call system_clock(startclock, clockrate)

  call compute_loop_tiled(ncolours, tile_x, tile_y, ncell_3d, nlayers, ndf1, undf1, ndf2, undf2, &
       & ntiles_per_colour, cmap_tiled, map1, map2, data1, data2, op_data, data1_snapshot)

  call system_clock(stopclock, clockrate)
  timing = dble(stopclock - startclock)/dble(clockrate)
  write(*,'(A, 1X, F7.2)') 'Timing:', timing

  !$acc end data

  write(*,'(A)') "lma_driver:Kernel run, checking answer ..."
  !check the answer
  count = compare(data1_snapshot, answer, undf1, .false.)
  write(*,'(A,I8,A,I8,A)') "lma_driver:checked ",undf1," answers, found ",count, " errors" 

  ! deallocate the arrays
  deallocate(ncells_per_colour, cmap)
  deallocate(map1, map2)
  deallocate(data1, data2)
  deallocate(op_data, answer)

end program lma_driver
