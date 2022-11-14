module compute_loop_mod

  use constants_mod, only: r_def, i_def, l_def

  implicit none

contains

  subroutine compute_loop_tiled(ncolours, tile_x, tile_y, ncell_3d, nlayers, ndf1, undf1, ndf2, undf2, &
       & ntiles_per_colour, cmap, map1, map2, data1, data2, op_data, data1_snapshot, colouring)

    use matrix_vector_kernel_mod, only: matrix_vector_code
    use omp_lib, only: omp_get_num_threads, omp_get_thread_num, omp_get_max_threads

    implicit none

    integer(kind=i_def), intent(in) :: ncolours, tile_x, tile_y, ncell_3d, nlayers, ndf1, undf1, ndf2, undf2
    integer(kind=i_def), intent(in) :: ntiles_per_colour, cmap(:,:,:), map1(:,:), map2(:,:)
    real(kind=r_def), intent(in) :: op_data(:,:,:), data2(:)
    real(kind=r_def), intent(inout) :: data1(:)
    real(kind=r_def), intent(out) :: data1_snapshot(:)
    logical(kind=l_def), intent(in) :: colouring

    integer(kind=i_def) :: count, colour, ncells_per_tile, tile, cell
    logical(kind=l_def) :: atomic

    write(*,'(A)') 'Running version compute_loop_tiled'

    ncells_per_tile = tile_x*tile_y

    ! Use atomic updates if running without colouring
    atomic = .not. colouring
    if (atomic) then
       write(*,'(A)') 'Running with atomic updates'
    else
       write(*,'(A)') 'Running without atomic updates'
    end if

    write(*,'(A,2(X,I6))') 'Number of tiles for parallel computation vs number of threads:', &
         & ntiles_per_colour, omp_get_max_threads()

    ! Repeat the work 1000 times to hide the cost of reading the data.
    do count = 1, 1000
       do colour = 1, ncolours

          !$omp parallel default(shared), private(tile, cell)
          if (omp_get_num_threads() .gt. ntiles_per_colour .and. count .eq. 1 .and. &
               & colour .eq. 1 .and. omp_get_thread_num() .eq. 1) then
             write(*,'(A)') 'WARNING: MORE THREADS AVAILABLE THAN TILES'
          end if
          !$omp do schedule(static)
          do tile = 1, ntiles_per_colour

             do cell = 1, ncells_per_tile
                call matrix_vector_code(cmap(colour,tile, cell), nlayers, data1, data2, &
                     ncell_3d, op_data, ndf1, undf1, map1(:,cmap(colour,tile,cell)), &
                     ndf2, undf2, map2(:,cmap(colour,tile,cell)), atomic)
             end do
          end do
          !$omp end do
          !$omp end parallel

       end do
       if (count.eq.1) data1_snapshot = data1
    end do

  end subroutine compute_loop_tiled


  subroutine compute_loop_inlined(ncolours, tile_x, tile_y, ncell_3d, nlayers, ndf1, undf1, ndf2, undf2, &
       & ntiles_per_colour, cmap, map1, map2, data1, data2, op_data, data1_snapshot)

    implicit none

    integer(kind=i_def), intent(in) :: ncolours, tile_x, tile_y, ncell_3d, nlayers, ndf1, undf1, ndf2, undf2
    integer(kind=i_def), intent(in) :: ntiles_per_colour, cmap(:,:,:), map1(:,:), map2(:,:)
    real(kind=r_def), intent(in) :: op_data(:,:,:), data2(:)
    real(kind=r_def), intent(inout) :: data1(:)
    real(kind=r_def), intent(out) :: data1_snapshot(:)

    integer(kind=i_def) :: count, colour, ncells_per_tile, tile, cell
    integer(kind=i_def) :: k, ik, df, df2, m1, m2

    write(*,*) 'Running version compute_loop_inlined'

    ncells_per_tile = tile_x*tile_y

    ! Repeat the work 1000 times to hide the cost of reading the data.
    do count = 1, 1000
       do colour = 1, ncolours

          !$omp parallel default(shared), private(tile, cell, k, ik, df2, m2, df, m1)
          !$omp do schedule(static)
          !$acc parallel loop gang num_gangs(ntiles_per_colour) num_workers(ncells_per_tile) &
          !$acc present(cmap, map1, map2, data1, data2, op_data)
          do tile = 1, ntiles_per_colour

             !$acc loop worker
             do cell = 1, ncells_per_tile
                !$omp simd
                !$acc loop vector
                do k = 0, nlayers-1
                   ik = (cmap(colour,tile,cell)-1)*nlayers
                   !$acc loop seq
                   do df2 = 1, ndf2
                      m2 = map2(df2, cmap(colour,tile,cell))
                      !$acc loop seq
                      do df = 1,ndf1
                         m1 = map1(df, cmap(colour,tile,cell))
                         !$acc atomic update
                         data1(m1+k) = data1(m1+k) + op_data(df,df2,ik+k+1)*data2(m2+k)
                      end do
                   end do
                end do

             end do
          end do
          !$acc end parallel
          !$omp end do
          !$omp end parallel
       end do
       if (count.eq.1) then
          !$acc update host(data1)
          data1_snapshot = data1
       end if
    end do

  end subroutine compute_loop_inlined


  ! subroutine compute_loop_nocolour_notile(ncolours, tile_x, tile_y, ncell_3d, nlayers, ndf1, undf1, ndf2, undf2, &
  !      & ntiles_per_colour, cmap, map1, map2, data1, data2, op_data, data1_snapshot)

  !   implicit none

  !   integer(kind=i_def), intent(in) :: ncolours, tile_x, tile_y, ncell_3d, nlayers, ndf1, undf1, ndf2, undf2
  !   integer(kind=i_def), intent(in) :: ntiles_per_colour, cmap(:,:,:), map1(:,:), map2(:,:)
  !   real(kind=r_def), intent(in) :: op_data(:,:,:), data2(:)
  !   real(kind=r_def), intent(inout) :: data1(:)
  !   real(kind=r_def), intent(out) :: data1_snapshot(:)

  !   integer(kind=i_def) :: count, colour, cell
  !   integer(kind=i_def) :: k, ik, df, df2, m1, m2
  !   integer(kind=i_def) :: ncolours_nocolour = 6
  !   integer(kind=i_def), allocatable :: ncells_per_colour_nocolour(:), cmap_nocolour(:,:)

  !   write(*,*) 'Running version compute_loop_nocolour_notile'

  !   ! Reset colouring to have one colour per cubed-sphere panel, for a consistent comparison
  !   ! with coloured loops that use distinct colours for each panel
  !   count = sum(ncells_per_colour)
  !   allocate(ncells_per_colour_nocolour(ncolours_nocolour))
  !   ncells_per_colour_nocolour(:) = count/6

  !   ! Rewrite colour map with contiguous cell IDs for each panel
  !   allocate(cmap_nocolour(ncolours_nocolour, maxval(ncells_per_colour_nocolour)))
  !   count = 1
  !   do colour = 1, ncolours_nocolour
  !      do cell = 1, ncells_per_colour_nocolour(colour)
  !         cmap_nocolour(colour, cell) = count
  !         count = count + 1
  !      end do
  !   end do

  !   !$acc data copyin(ncells_per_colour_nocolour, cmap_nocolour)

  !   ! Repeat the work 1000 times to hide the cost of reading the data.
  !   do count = 1, 1000
  !      do colour = 1, ncolours_nocolour

  !         !$omp parallel default(shared), private(cell, k, ik, df2, m2, df, m1)
  !         !$omp do schedule(static)
  !         !$acc parallel loop &
  !         !$acc present(ncells_per_colour_nocolour, cmap_nocolour, map1, map2, data1, data2, op_data)
  !         do cell = 1, ncells_per_colour_nocolour(colour)
  !            !$omp simd
  !            do k = 0, nlayers-1
  !               ik = (cmap_nocolour(colour,cell)-1)*nlayers
  !               !$acc loop seq
  !               do df2 = 1, ndf2
  !                  m2 = map2(df2, cmap_nocolour(colour,cell))
  !                  !$acc loop seq
  !                  do df = 1,ndf1
  !                     m1 = map1(df, cmap_nocolour(colour,cell))
  !                     !$acc atomic update
  !                     data1(m1+k) = data1(m1+k) + op_data(df,df2,ik+k+1)*data2(m2+k)
  !                  end do
  !               end do
  !            end do
  !         end do
  !         !$acc end parallel
  !         !$omp end do
  !         !$omp end parallel
  !      end do
  !      if (count.eq.1) then
  !         !$acc update host(data1)
  !         data1_snapshot = data1
  !      end if
  !   end do

  !   !$acc end data

  ! end subroutine compute_loop_nocolour_notile

end module compute_loop_mod
