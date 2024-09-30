module compute_loop_mod

  use constants_mod, only: r_def, i_def, l_def

  implicit none

contains

  subroutine invoke_matrix_vector_kernel(ncell_3d, nlayers, ndf1, undf1, ndf2, undf2, &
       map1, map2, data1, data2, op_data, data1_snapshot, inline, ntimes)

    use matrix_vector_kernel_mod, only: matrix_vector_code

    implicit none

    integer(kind=i_def), intent(in) :: ncell_3d, nlayers, ndf1, undf1, ndf2, undf2
    integer(kind=i_def), intent(in) :: map1(:,:), map2(:,:)
    real(kind=r_def), intent(in) :: op_data(:,:,:), data2(:)
    real(kind=r_def), intent(inout) :: data1(:)
    real(kind=r_def), intent(out) :: data1_snapshot(:)
    logical(kind=l_def), intent(in) :: inline
    integer(kind=i_def), intent(in) :: ntimes

    integer(kind=i_def) :: count, panel, cell, first_cell, last_cell

    ! Inlined kernel
    integer(kind=i_def) :: df, k, ik, df2

    if (inline) then
      write(*,'(A,X,I4,X,A)') 'Running invoke_matrix_vector_kernel', ntimes, 'times - inlined'
    else
      write(*,'(A,X,I4,X,A)') 'Running invoke_matrix_vector_kernel', ntimes, 'times - not inlined'
    end if

    ! Repeat loop many times to hide the cost of reading the data
    do count = 1, ntimes

      ! Compute each panel separately, to simulate partitioning into rectangles
      do panel = 1, 6
        first_cell = (panel-1)*(ncell_3d/nlayers)/6+1
        last_cell = panel*(ncell_3d/nlayers)/6

        if (inline) then
          !$acc parallel present(map1, map2, data1, data2, op_data) &
          !$acc default(none) firstprivate(first_cell, last_cell, df, ndf1, df2, ndf2, k, nlayers, ik)
          !$acc loop
          do cell = first_cell, last_cell
            do df2 = 1, ndf2
              do df = 1,ndf1
                do k = 0, nlayers-1
                  ik = (cell-1)*nlayers+k+1
                  !$acc atomic update
                  data1(map1(df, cell)+k) = data1(map1(df, cell)+k) + op_data(ik,df,df2)*data2(map2(df2, cell)+k)
                end do
              end do
            end do
          end do
          !$acc end parallel
        else
          !$acc parallel present(map1, map2, data1, data2, op_data) vector_length(128) &
          !$acc default(none) firstprivate(first_cell, last_cell, nlayers, ncell_3d, ndf1, undf1, ndf2, undf2)
          !$acc loop
          do cell = first_cell, last_cell
            call matrix_vector_code(cell, nlayers, data1, data2, &
                 ncell_3d, op_data, ndf1, undf1, map1(:, cell), &
                 ndf2, undf2, map2(:, cell))
          end do
          !$acc end parallel
        end if
      end do

      ! Array data1 will be overwritten in each iteration, keep a snapshot for KGO checking
      if (count.eq.1) then
        !$acc update host(data1)
        data1_snapshot = data1
      end if
    end do

  end subroutine invoke_matrix_vector_kernel

  subroutine invoke_matrix_vector_kernel_coloured(ncolours, ncell_3d, nlayers, ndf1, undf1, ndf2, undf2, &
       ncells_per_colour, cmap, map1, map2, data1, data2, op_data, data1_snapshot, inline, ntimes)

    use matrix_vector_kernel_mod, only: matrix_vector_code
    use omp_lib, only: omp_get_num_threads, omp_get_thread_num

    implicit none

    integer(kind=i_def), intent(in) :: ncolours, ncell_3d, nlayers, ndf1, undf1, ndf2, undf2
    integer(kind=i_def), intent(in) :: ncells_per_colour(:), cmap(:,:), map1(:,:), map2(:,:)
    real(kind=r_def), intent(in) :: op_data(:,:,:), data2(:)
    real(kind=r_def), intent(inout) :: data1(:)
    real(kind=r_def), intent(out) :: data1_snapshot(:)
    logical(kind=l_def), intent(in) :: inline
    integer(kind=i_def), intent(in) :: ntimes

    integer(kind=i_def) :: count, panel, colour, cell, first_cell, last_cell

    ! Inlined kernel
    integer(kind=i_def) :: cell_number, df, k, ik, df2

    if (inline) then
      write(*,'(A,X,I4,X,A)') 'Running invoke_matrix_vector_kernel_coloured', ntimes, 'times - inlined'
    else
      write(*,'(A,X,I4,X,A)') 'Running invoke_matrix_vector_kernel_coloured', ntimes, 'times - not inlined'
    end if

    ! Repeat loop many times to hide the cost of reading the data
    do count = 1, ntimes

      ! Compute each panel separately, to simulate partitioning into rectangles
      do panel = 1, 6
        do colour = 1, ncolours
          first_cell = (panel-1)*ncells_per_colour(colour)/6+1
          last_cell = panel*ncells_per_colour(colour)/6

          if (inline) then
            !$acc parallel present(cmap, map1, map2, data1, data2, op_data) &
            !$acc default(none) firstprivate(first_cell, last_cell, df, ndf1, df2, ndf2, k, nlayers, ik, colour, cell)
            !$acc loop
            do cell_number = first_cell, last_cell
              cell = cmap(colour, cell_number)
              do df2 = 1, ndf2
                do df = 1,ndf1
                  do k = 0, nlayers-1
                    ik = (cell-1)*nlayers+k+1
                    !$acc atomic update
                    data1(map1(df, cell)+k) = data1(map1(df, cell)+k) + op_data(ik,df,df2)*data2(map2(df2, cell)+k)
                  end do
                end do
              end do
            end do
            !$acc end parallel
          else
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
          end if
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
       ntiles_per_colour, ncells_per_tile, tmap, map1, map2, data1, data2, op_data, data1_snapshot, inline, ntimes)

    use matrix_vector_kernel_mod, only: matrix_vector_code
    use omp_lib, only: omp_get_num_threads, omp_get_thread_num

    implicit none

    integer(kind=i_def), intent(in) :: ncolours, ncell_3d, nlayers, ndf1, undf1, ndf2, undf2
    integer(kind=i_def), intent(in) :: ntiles_per_colour, ncells_per_tile, tmap(:,:,:), map1(:,:), map2(:,:)
    real(kind=r_def), intent(in) :: op_data(:,:,:), data2(:)
    real(kind=r_def), intent(inout) :: data1(:)
    real(kind=r_def), intent(out) :: data1_snapshot(:)
    logical(kind=l_def), intent(in) :: inline
    integer(kind=i_def), intent(in) :: ntimes

    integer(kind=i_def) :: count, panel, colour, tile, cell, first_tile, last_tile

    ! Inlined kernel
    integer(kind=i_def) :: cell_number, df, k, ik, df2

    if (inline) then
      write(*,'(A,X,I4,X,A)') 'Running invoke_matrix_vector_kernel_coloured_tiled', ntimes, 'times - inlined'
    else
      write(*,'(A,X,I4,X,A)') 'Running invoke_matrix_vector_kernel_coloured_tiled', ntimes, 'times - not inlined'
    end if

    do count = 1, ntimes
      do panel = 1, 6
        first_tile = (panel-1)*ntiles_per_colour/6+1
        last_tile = panel*ntiles_per_colour/6

        do colour = 1, ncolours

          if (inline) then
            !$acc parallel present(tmap, map1, map2, data1, data2, op_data) &
            !$acc default(none) firstprivate(last_tile, first_tile, cell_number, ncells_per_tile, df) &
            !$acc firstprivate(ndf1, df2, ndf2, k, nlayers, ik, colour, cell)
            !$acc loop
            do tile = first_tile, last_tile
              do cell_number = 1, ncells_per_tile
                cell = tmap(colour, tile, cell_number)
                do df2 = 1, ndf2
                  do df = 1,ndf1
                    do k = 0, nlayers-1
                      ik = (cell-1)*nlayers+k+1
                      !$acc atomic update
                      data1(map1(df, cell)+k) = data1(map1(df, cell)+k) + op_data(ik,df,df2)*data2(map2(df2, cell)+k)
                    end do
                  end do
                end do
              end do
            end do
            !$acc end parallel
          else
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
            !$acc loop
            do tile = first_tile, last_tile
              do cell = 1, ncells_per_tile
                call matrix_vector_code(tmap(colour, tile, cell), nlayers, data1, data2, &
                     ncell_3d, op_data, ndf1, undf1, map1(:, tmap(colour, tile, cell)), &
                     ndf2, undf2, map2(:, tmap(colour, tile, cell)))
              end do
            end do
            !$omp end do

            !$acc end parallel
            !$omp end parallel
          end if
        end do
      end do

      if (count.eq.1) then
        !$acc update host(data1)
        data1_snapshot = data1
      end if
    end do

  end subroutine invoke_matrix_vector_kernel_coloured_tiled

  subroutine invoke_matrix_vector_kernel_coloured_tiled_vert(ncolours, nlayers, ndf1, ndf2, &
       ntiles_per_colour, ncells_per_tile, tmap, tile, kstart, kstop, nblocks, map1, map2, data1, data2, op_data, data1_snapshot, inline, ntimes)

    implicit none

    integer(kind=i_def), intent(in) :: ncolours, nlayers, ndf1, ndf2
    integer(kind=i_def), intent(in) :: ntiles_per_colour, ncells_per_tile, tmap(:,:,:), map1(:,:), map2(:,:)
    integer(kind=i_def), intent(in) :: tile(:), kstart(:), kstop(:), nblocks
    real(kind=r_def), intent(in) :: op_data(:,:,:), data2(:)
    real(kind=r_def), intent(inout) :: data1(:)
    real(kind=r_def), intent(out) :: data1_snapshot(:)
    logical(kind=l_def), intent(in) :: inline
    integer(kind=i_def), intent(in) :: ntimes

    integer(kind=i_def) :: count, panel, colour, cell, first_tile

    ! Inlined kernel
    integer(kind=i_def) :: cell_number, df, k, ik, df2

    ! Vertical tiling
    integer(kind=i_def) :: iblock

    if (inline) then
      write(*,'(A,X,I4,X,A)') 'Running invoke_matrix_vector_kernel_coloured_tiled_vert', ntimes, 'times - inlined'
    else
      write(*,'(A,X,I4,X,A)') 'invoke_matrix_vector_kernel_coloured_tiled_vert: only inlined variant is available'
      stop
    end if

    do count = 1, ntimes
      do panel = 1, 6
        first_tile = (panel-1)*ntiles_per_colour/6+1
        do colour = 1, ncolours

          !$acc parallel present(tile, kstart, kstop, tmap, map1, map2, data1, data2, op_data) &
          !$acc default(none) firstprivate(first_tile, cell_number, ncells_per_tile, df) &
          !$acc firstprivate(ndf1, df2, ndf2, k, nlayers, ik, colour, cell, nblocks) &
          !$acc vector_length(128)

          ! Outer loop over 3D blocks (2D tiles + 1D sections)
          !$acc loop gang
          do iblock = 1, nblocks
            do cell_number = 1, ncells_per_tile
              cell = tmap(colour, first_tile + tile(iblock), cell_number)
              do df2 = 1, ndf2
                do df = 1,ndf1
                  do k = kstart(iblock), kstop(iblock)
                    ik = (cell-1)*nlayers+k+1
                    !$acc atomic update
                    data1(map1(df, cell)+k) = data1(map1(df, cell)+k) + op_data(ik,df,df2)*data2(map2(df2, cell)+k)
                  end do
                end do
              end do
            end do
          end do
          !$acc end parallel
        end do
      end do

      if (count.eq.1) then
        !$acc update host(data1)
        data1_snapshot = data1
      end if
    end do

  end subroutine invoke_matrix_vector_kernel_coloured_tiled_vert

end module compute_loop_mod
