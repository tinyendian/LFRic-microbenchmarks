module compute_loop_mod

  use constants_mod, only: r_def, i_def

  implicit none

contains

  subroutine invoke_matrix_vector_kernel_coloured(ncolours, ncell_3d, nlayers, ndf1, undf1, ndf2, undf2, &
       ncells_per_colour, cmap, map1, map2, data1, data2, op_data, data1_snapshot)

    use matrix_vector_kernel_mod, only: matrix_vector_code
    use omp_lib, only: omp_get_num_threads, omp_get_thread_num

    implicit none

    integer(kind=i_def), intent(in) :: ncolours, ncell_3d, nlayers, ndf1, undf1, ndf2, undf2
    integer(kind=i_def), intent(in) :: ncells_per_colour(:), cmap(:,:), map1(:,:), map2(:,:)
    real(kind=r_def), intent(in) :: op_data(:,:,:), data2(:)
    real(kind=r_def), intent(inout) :: data1(:)
    real(kind=r_def), intent(out) :: data1_snapshot(:)

    integer(kind=i_def) :: count, panel, colour, cell, first_cell, last_cell

    write(*,'(A)') 'Running version invoke_matrix_vector_kernel_coloured'

    ! Repeat loop many times to hide the cost of reading the data
    do count = 1, 1000

      ! Compute each panel separately, to simulate partitioning into rectangles
      do panel = 1, 6
        do colour = 1, ncolours
          first_cell = (panel-1)*ncells_per_colour(colour)/6+1
          last_cell = panel*ncells_per_colour(colour)/6

          !$omp parallel default(none) private(cell) shared(count, colour, ncells_per_colour) &
          !$             shared(cmap, nlayers, data1, data2, ncell_3d, op_data, ndf1, undf1)  &
          !$             shared(map1, ndf2, undf2, map2, panel, ncolours, first_cell, last_cell)

          ! Sanity checks (OpenMP parallelisation relies on colouring)
          if (omp_get_num_threads() > 1 .and. ncolours == 1 .and. count == 1 .and. panel == 1 &
               .and. omp_get_thread_num() == 1) then
            write(*,'(A)') 'WARNING: RUNNING WITH MULTIPLE THREADS BUT WITHOUT COLOURING!'
          else if (omp_get_num_threads() > ncells_per_colour(colour) .and. count == 1 .and. &
               & colour == 1 .and. panel == 1 .and. omp_get_thread_num() == 1) then
            write(*,'(A)') 'WARNING: There are more threads than cells per colour'
          end if

          !$acc parallel present(ncells_per_colour, cmap, map1, map2, data1, data2, op_data) &
          !$acc          firstprivate(panel, colour, nlayers, ncell_3d, ndf1, undf1, ndf2) &
          !$acc          firstprivate(undf2, first_cell, last_cell) default(none)

          !$omp do schedule(static)
          !$acc loop
          do cell = first_cell, last_cell
            call matrix_vector_code(cmap(colour, cell), nlayers, data1, data2, &
                 ncell_3d, op_data, ndf1, undf1, map1(:, cmap(colour, cell)), &
                 ndf2, undf2, map2(:, cmap(colour, cell)))
          end do
          !$omp end do
          !$acc end parallel
          !$omp end parallel

        end do
      end do

      ! Array data1 will be overwritten in each iteration, keep a snapshot for KGO checking
      if (count.eq.1) then
        !$acc update host(data1)
        data1_snapshot = data1
      end if
    end do

  end subroutine invoke_matrix_vector_kernel_coloured


  subroutine invoke_matrix_vector_kernel_coloured_tiled(ncolours, ncell_3d, nlayers, ndf1, undf1, ndf2, undf2, &
       ntiles_per_colour, ncells_per_tile, tmap, map1, map2, data1, data2, op_data, data1_snapshot)

    use matrix_vector_kernel_mod, only: matrix_vector_code
    use omp_lib, only: omp_get_num_threads, omp_get_thread_num

    implicit none

    integer(kind=i_def), intent(in) :: ncolours, ncell_3d, nlayers, ndf1, undf1, ndf2, undf2
    integer(kind=i_def), intent(in) :: ntiles_per_colour, ncells_per_tile, tmap(:,:,:), map1(:,:), map2(:,:)
    real(kind=r_def), intent(in) :: op_data(:,:,:), data2(:)
    real(kind=r_def), intent(inout) :: data1(:)
    real(kind=r_def), intent(out) :: data1_snapshot(:)

    integer(kind=i_def) :: count, panel, colour, tile, cell, first_tile, last_tile

    write(*,'(A)') 'Running version invoke_matrix_vector_kernel_coloured_tiled'

    do count = 1, 1000
      do panel = 1, 6
          first_tile = (panel-1)*ntiles_per_colour/6+1
          last_tile = panel*ntiles_per_colour/6

        do colour = 1, ncolours

          !$omp parallel default(none) private(tile, cell) shared(count, colour, ntiles_per_colour) &
          !$             shared(ncells_per_tile, tmap, nlayers, data1, data2, ncell_3d, op_data)  &
          !$             shared(ndf1, undf1, map1, ndf2, undf2, map2, panel, ncolours, first_tile) &
          !$             shared(last_tile)
          if (omp_get_num_threads() > 1 .and. ncolours == 1 .and. count == 1 .and. panel == 1 &
               .and. omp_get_thread_num() == 1) then
            write(*,'(A)') 'WARNING: RUNNING WITH MULTIPLE THREADS BUT WITHOUT COLOURING!'
          else if (omp_get_num_threads() > ntiles_per_colour .and. count == 1 .and. &
               & colour == 1 .and. omp_get_thread_num() == 1) then
            write(*,'(A)') 'WARNING: MORE THREADS AVAILABLE THAN TILES PER COLOUR'
          end if


          !$acc parallel present(ntiles_per_colour, tmap, map1, map2, data1, data2, op_data) &
          !$acc          firstprivate(panel, colour, nlayers, ncell_3d, ndf1, undf1, ndf2) &
          !$acc          firstprivate(undf2, first_tile, last_tile) default(none) private(cell)

          ! Tile numbers are global across all 6 cubed-sphere panels
          !$omp do schedule(static)
          do tile = first_tile, last_tile
            !$acc loop
            do cell = 1, ncells_per_tile
              call matrix_vector_code(tmap(colour, tile, cell), nlayers, data1, data2, &
                   ncell_3d, op_data, ndf1, undf1, map1(:, tmap(colour, tile, cell)), &
                   ndf2, undf2, map2(:, tmap(colour, tile, cell)))
            end do
          end do
          !$omp end do
          !$acc end parallel
          !$omp end parallel

        end do
      end do

      if (count.eq.1) then
        !$acc update host(data1)
        data1_snapshot = data1
      end if
    end do

  end subroutine invoke_matrix_vector_kernel_coloured_tiled

end module compute_loop_mod
