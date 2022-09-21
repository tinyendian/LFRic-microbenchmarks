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
  use colourist_mod, only: compute_colour_map, reorder_data, &
       & check_cell_order, check_colours
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
  integer(kind=i_def) :: colour, cell, i

  integer(kind=i_def) :: count

  ! Tiles
  integer(kind=i_def) :: cells_per_dimension = 48, tile_x = 1, tile_y = 1

  ! Timing
  integer(kind=8) :: startclock, stopclock, clockrate
  real(kind=8) :: timing

  ! CLI
  logical(kind=l_def) :: write_cell_props = .false.
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
  write(*,*) "lma_driver:ingested dinodump",ndf1,ndf2  

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
     end select
  end do
  write(*, '(A,2(X,I))') ' Using tiling configuration', tile_x, tile_y

  ! Compute new colour map
  deallocate(ncells_per_colour)
  ncolours = 24
  allocate( ncells_per_colour(ncolours) )
  call compute_colour_map(tile_x, tile_y, cells_per_dimension, cmap, ncells_per_colour, write_cell_props)

  ! call reorder_data(ncell, ndf1, nlayers, undf1, reshape(transpose(cmap), (/ncell/)), map1, data1)
  ! call reorder_data(ncell, ndf2, nlayers, undf2, reshape(transpose(cmap), (/ncell/)), map2, data2)

  call check_colours(tile_x, tile_y, cells_per_dimension, ncell, ndf1, map1, ncolours, maxval(ncells_per_colour), cmap)
  call check_cell_order(tile_x, tile_y, cells_per_dimension, ncolours, maxval(ncells_per_colour), cmap, ncells_per_colour)
  !
  ! Repeat the work 1000 times to hide the cost of reading the data.
  call system_clock(startclock, clockrate)
  do count =1 , 100
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
  call system_clock(stopclock, clockrate)
  timing = dble(stopclock - startclock)/dble(clockrate)

  write(*,'(A, 1X, F6.1)') 'Timing:', timing

  write(*,'(A)') "lma_driver:Kernel run, checking answer ..."
  !check the answer
  count = compare(data1_snapshot, answer, undf1, .false.)
  write(*,'(A,I6,A,I6,A)') "lma_driver:checked ",undf1," answers, found ",count, " errors" 
  
  ! deallocate the arrays
  deallocate(ncells_per_colour, cmap)
  deallocate(map1, map2)
  deallocate(data1, data2)
  deallocate(op_data, answer)

end program lma_driver
