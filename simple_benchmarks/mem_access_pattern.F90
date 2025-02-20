program mem_access_pattern
  implicit none
  integer(kind=4) :: nxsize, nysize, veclength, ntimes
  integer(kind=4) :: k, stride, nxsize_delta
  real(kind=4) :: timing

  ! Test performance of different memory access patterns on GPUs:
  ! - Match/mismatch between loop vector length and contiguous array dimension
  ! - Strided access
  ! The different settings will require different numbers of memory
  ! transactions for a saxpy/daxpy calculation, which may affect throughput
  ! Each tests reports runtime per array element to measure efficiency
  ! IMPORTANT: Total array size must be larger than last-level cache size,
  !            to ensure that the cache does not compensate for non-optimal
  !            memory access patterns

  ! Set parameters
  ! Contiguous array dimension
  nxsize = 32
  ! Second array dimension
  nysize = 100000
  ! Inner loop vector length
  veclength = 32
  ! Repeat benchmark for robust measurements
  ntimes = 1

  call initialise_random()

  ! Vary contiguous array size
  do nxsize_delta = -30, 30, 5
    timing = 0.0
    do k = 1, ntimes
      timing = timing + saxpy_2d(nxsize+nxsize_delta, nysize, veclength)
    end do
    print '(3(A,X,I6,X),A,X,E12.6)', 'Saxpy nx=', nxsize+nxsize_delta, 'ny=', nysize, 'veclength=', veclength, 'timing [s]:', timing/dble(ntimes)
  end do

  ! Vary inner loop stride
  do stride = 1, 8
    timing = 0.0
    do k = 1, ntimes
      timing = timing + saxpy_2d_strided(nxsize, nysize, stride, veclength)
    end do
    print '(3(A,X,I6,X),A,X,E12.6)', 'Saxpy strided nx=', nxsize, 'ny=', nysize, 'veclength=', veclength, 'timing [s]:', timing/dble(ntimes)
  end do

  ! FP64 variant
  timing = 0.0
  do k = 1, ntimes
    timing = timing + daxpy_2d(nxsize, nysize, veclength)
  end do
  print '(3(A,X,I6,X),A,X,E12.6)', 'Daxpy nx=', nxsize, 'ny=', nysize, 'veclength=', veclength, 'timing [s]:', timing/dble(ntimes)

contains

  subroutine initialise_random()
    implicit none
    integer(kind=4) :: seed_size
    integer(kind=8) :: startclock, clockrate
    integer(kind=8), allocatable :: seedarr(:)

    ! Set random seed using system_clock
    ! Call system_clock twice, as startclock may be zero on first call
    call system_clock(startclock, clockrate)
    call random_seed(size=seed_size)
    call system_clock(startclock, clockrate)
    allocate(seedarr(seed_size), source=startclock)
    call random_seed(put=seedarr)
    deallocate(seedarr)
  end subroutine initialise_random

  ! ======================================================================================================

  function saxpy_2d(nx, ny, acc_vector_length) result(timing)
    implicit none
    integer, intent(in) :: nx, ny, acc_vector_length
    real(kind=8) :: timing
    integer(kind=4) :: i, j
    real(kind=4) :: dataarr(nx, ny), kgoarr(nx, ny)
    integer(kind=8) :: startclock, stopclock, clockrate

    call random_number(dataarr)
    kgoarr(:,:) = dataarr(:,:)

    !$acc data copy(dataarr)

    call system_clock(startclock, clockrate)

    !$acc parallel vector_length(acc_vector_length) default(none) firstprivate(nx,ny) present(dataarr)

    !$acc loop gang
    do j = 1, ny
      !$acc loop vector
      do i = 1, nx
        dataarr(i,j) = 2.345_4*dataarr(i,j) + 0.432_4
      end do
    end do

    !$acc end parallel

    call system_clock(stopclock, clockrate)

    !$acc end data

    ! Compare with CPU result
    kgoarr = 2.345_4*kgoarr + 0.432_4
    if (maxval(abs((dataarr-kgoarr)/kgoarr)) > 1.0e-7) print *, 'WARNING - GPU and CPU results do not match'

    timing = dble(stopclock - startclock)/dble(clockrate*nx*ny)

  end function saxpy_2d

  ! ======================================================================================================

  function saxpy_2d_strided(nx, ny, dx, acc_vector_length) result(timing)
    implicit none
    integer, intent(in) :: nx, ny, dx, acc_vector_length
    real(kind=8) :: timing
    integer(kind=4) :: i, j
    real(kind=4) :: dataarr(nx*dx, ny), kgoarr(nx*dx, ny)
    integer(kind=8) :: startclock, stopclock, clockrate

    call random_number(dataarr)
    kgoarr(:,:) = dataarr(:,:)

    !$acc data copy(dataarr)

    call system_clock(startclock, clockrate)

    !$acc parallel vector_length(acc_vector_length) default(none) firstprivate(nx,ny,dx) present(dataarr)
    !$acc loop gang
    do j = 1, ny
      !$acc loop vector
      do i = 1, nx*dx, dx
        dataarr(i,j) = 2.345_4*dataarr(i,j) + 0.432_4
      end do
    end do
    !$acc end parallel

    call system_clock(stopclock, clockrate)

    !$acc end data

    kgoarr(1:nx*dx:dx,:) = 2.345_4*kgoarr(1:nx*dx:dx,:) + 0.432_4
    if (maxval(abs((dataarr-kgoarr)/kgoarr)) > 1.0e-7) print *, 'WARNING - GPU and CPU results do not match'

    timing = dble(stopclock - startclock)/dble(clockrate*nx*ny)

  end function saxpy_2d_strided

  ! ======================================================================================================

  function daxpy_2d(nx, ny, acc_vector_length) result(timing)
    implicit none
    integer, intent(in) :: nx, ny, acc_vector_length
    real(kind=8) :: timing
    integer(kind=4) :: i, j
    real(kind=8) :: dataarr(nx, ny), kgoarr(nx, ny)
    integer(kind=8) :: startclock, stopclock, clockrate

    call random_number(dataarr)
    kgoarr(:,:) = dataarr(:,:)

    !$acc data copy(dataarr)

    call system_clock(startclock, clockrate)

    !$acc parallel vector_length(acc_vector_length) default(none) firstprivate(nx,ny) present(dataarr)
    !$acc loop gang
    do j = 1, ny
      !$acc loop vector
      do i = 1, nx
        dataarr(i,j) = 2.345_8*dataarr(i,j) + 0.432_8
      end do
    end do
    !$acc end parallel

    call system_clock(stopclock, clockrate)

    !$acc end data

    ! Compare with CPU result
    kgoarr = 2.345_8*kgoarr + 0.432_8
    if (maxval(abs((dataarr-kgoarr)/kgoarr)) > 1.0e-7) print *, 'WARNING - GPU and CPU results do not match'

    timing = dble(stopclock - startclock)/dble(clockrate*nx*ny)

  end function daxpy_2d

end program mem_access_pattern
