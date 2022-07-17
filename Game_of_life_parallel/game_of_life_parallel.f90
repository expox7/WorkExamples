program game_of_life
    use mpi_f08
    implicit none
    logical, dimension(:, :), pointer :: old_world, new_world, tmp_world
    type(MPI_Comm) :: comm
    type(MPI_Status) :: status
    integer(kind=MPI_ADDRESS_KIND) :: lb, real_extent
    type(MPI_Datatype) :: a_col, a_row, a_row2, a_tmp_row
    integer :: n_ranks, my_rank, root, north_rank, south_rank, west_rank, east_rank
    integer :: n_rows, n_cols, row, col, ib, ie, jb, je, rb, re, cb, ce
    integer :: height, width, i
    integer :: max_gen, gen
    logical :: still

    call MPI_Init()
    comm = MPI_COMM_WORLD
    call MPI_Comm_size(comm, n_ranks)
    call MPI_Comm_rank(comm, my_rank)

    root = 0

    if (my_rank == root) read *, height, width, max_gen, n_rows, n_cols
    call MPI_Bcast(height, 1, MPI_INTEGER, root, comm)
    call MPI_Bcast(width, 1, MPI_INTEGER, root, comm)
    call MPI_Bcast(max_gen, 1, MPI_INTEGER, root, comm)
    call MPI_Bcast(n_rows, 1, MPI_INTEGER, root, comm)
    call MPI_Bcast(n_cols, 1, MPI_INTEGER, root, comm)

    ! Create a 2D topology.
    if (my_rank == root .and. n_rows * n_cols /= n_ranks) then
        print "(a, 1X, i0, a)", "Run this code with `-np", n_rows * n_cols, "`."
        call MPI_Abort(comm, MPI_ERR_TOPOLOGY)
    end if

    call get_coords(my_rank, n_rows, n_cols, row, col)

    ! Define neighbor ranks with periodic boundaries.
    north_rank = get_rank(row - 1, col,     n_rows, n_cols)
    south_rank = get_rank(row + 1, col,     n_rows, n_cols)
    west_rank  = get_rank(row,     col - 1, n_rows, n_cols)
    east_rank  = get_rank(row,     col + 1, n_rows, n_cols)

    ! Partition the global 2D array into a Cartisian product of two 1D array
    ! partitions for each of the two dimensions.
    call partition(row, n_rows, height, ib, ie)
    call partition(col, n_cols, width, jb, je)

    ! Create sub-arrays for the grid.  Include two outer cells on
    ! both sides for ghost columns and rows.
    allocate(old_world(ib - 1:ie + 1, jb - 1:je + 1))
    allocate(new_world(ib - 1:ie + 1, jb - 1:je + 1))

    ! One column covers the internal range, [ib:ie].
    call MPI_Type_contiguous(ie - ib + 1, MPI_LOGICAL, a_col)
    call MPI_Type_commit(a_col)

    ! One row covers the internal range and two outer ghost cells, [jb-1:je+1].
    call MPI_Type_get_extent(MPI_LOGICAL, lb, real_extent)
    call MPI_Type_vector(je - jb + 3, 1, ie - ib + 3, MPI_LOGICAL, a_tmp_row)
    call MPI_Type_create_resized(a_tmp_row, lb, real_extent, a_row)
    call MPI_Type_commit(a_row)

    ! One row covers the internal range, [jb:je].
    call MPI_Type_get_extent(MPI_LOGICAL, lb, real_extent)
    call MPI_Type_vector(je - jb + 1, 1, ie - ib + 3, MPI_LOGICAL, a_tmp_row)
    call MPI_Type_create_resized(a_tmp_row, lb, real_extent, a_row2)
    call MPI_Type_commit(a_row2)

    call read_map()
    call update_borders(old_world, ib, jb, ie, je)

    do gen = 1, max_gen
        if (my_rank == root) print "(a, i0)", "Generation ", gen
        call print_map()
        call next_gen(old_world, new_world, ib, ie, jb, je)
        call update_borders(new_world, ib, jb, ie, je)

        !call wait_cls(100) !-> Activate to print the game in terminal.

        ! Convergence test
        call MPI_Allreduce(world_is_still(old_world, new_world), still, 1, MPI_LOGICAL, MPI_LAND, comm)
        if (still) exit

        ! Swap maps
        tmp_world => old_world; old_world => new_world; new_world => tmp_world
    end do
    
    if (associated(old_world)) deallocate(old_world)
    if (associated(new_world)) deallocate(new_world)

    call MPI_Type_free(a_col)
    call MPI_Type_free(a_row)
    call MPI_Type_free(a_row2)

    call MPI_Finalize()

contains

    logical function world_is_still(old_map, new_map)
        logical, dimension(:, :), pointer, intent(in) :: old_map, new_map
        world_is_still = all(old_map .eqv. new_map)
        return
    end function world_is_still

    !---------------------------------------------------------------------------

    integer function get_rank(row, col, n_rows, n_cols)
        integer, intent(in) :: row, col, n_rows, n_cols
        integer             :: t_row, t_col

        if (0 <= col) then
            t_col = merge(0, col, col == n_cols)
        else
            t_col = n_cols - 1
        end if

        if (0 <= row) then
            t_row = merge(0, row, row == n_rows)
        else
            t_row = n_rows - 1
        end if

        get_rank = t_row + t_col * n_rows
        return
    end function get_rank

    subroutine get_coords(rank, n_rows, n_cols, row, col)
        integer, intent(in)    :: rank, n_rows, n_cols
        integer, intent(inout) :: row, col
        
        row = modulo(rank, n_rows)
        col = (rank - row) / n_rows
        if (0 <= col .and. col < n_cols) then
            return
        else
            print "(a, 2(i0, a))", "get_coords: rank ", rank, &
                " is outside the column range [0, ", n_cols, ")."
            call MPI_Abort(comm, MPI_ERR_TOPOLOGY)
        end if
    end subroutine get_coords

    subroutine partition(id, n_ids, size, b, e)
        integer, intent(in)    :: id, n_ids, size
        integer, intent(inout) :: b, e
        integer :: remainder, quotient

        remainder = modulo(size, n_ids)
        quotient  = (size - remainder) / n_ids
        b = 1 + quotient * (id    ) + min(remainder, id    )
        e =     quotient * (id + 1) + min(remainder, id + 1)
        return
    end subroutine partition

    !---------------------------------------------------------------------------

    subroutine update_borders(map, ib, jb, ie, je)
        logical, dimension(:, :), pointer, intent(inout) :: map
        integer, intent(in) :: ib, jb, ie, je

        ! Synchronize columns.
        call MPI_Sendrecv(map(ib, jb    ), 1, a_col, west_rank, 1, &
                          map(ib, je + 1), 1, a_col, east_rank, 1, comm, status)
        call MPI_Sendrecv(map(ib, je    ), 1, a_col, east_rank, 2, &
                          map(ib, jb - 1), 1, a_col, west_rank, 2, comm, status)

        ! Synchronize rows.
        call MPI_Sendrecv(map(ib,     jb - 1), 1, a_row, north_rank, 3, &
                          map(ie + 1, jb - 1), 1, a_row, south_rank, 3, comm, status)
        call MPI_Sendrecv(map(ie,     jb - 1), 1, a_row, south_rank, 4, &
                          map(ib - 1, jb - 1), 1, a_row, north_rank, 4, comm, status)
        return
    end subroutine update_borders

    subroutine next_gen(old_map, new_map, ib, ie, jb, je)
        logical, dimension(:, :), pointer, intent(inout) :: old_map, new_map
        integer, intent(in) :: ib, ie, jb, je
        integer :: i, j
        integer :: c ! number of live neighbors

        do j = jb, je
            do i = ib, ie
                c = count(old_map(i - 1:i + 1, j - 1:j + 1))
                if (old_map(i, j)) then ! cell is live
                    new_map(i, j) = merge(.true., .false., 3 <= c .and. c <= 4)
                else ! cell is dead
                    new_map(i, j) = merge(.true., .false., c == 3)
                end if
            end do
        end do
        return
    end subroutine next_gen

    !---------------------------------------------------------------------------

    subroutine read_line(logic_line, w)
        logical, dimension(:), pointer, intent(inout) :: logic_line
        integer, intent(in) :: w
        character(len=:), allocatable :: line
        integer :: i
        
        allocate(character(len=w) :: line)
        read *, line
        do i = 1, w
            select case (line(i:i))
            case ('X')
                logic_line(i:i) = .true.
            case ('.')
                logic_line(i:i) = .false.
            case default
                stop "read_map: wrong input character `" // line(i:i) // "`"
            end select
        end do
        if (allocated(line)) deallocate(line)
        return
    end subroutine read_line

    subroutine read_map()
        if (my_rank == root) then
            block
                logical, dimension(:), pointer :: line
                integer :: dst
                allocate(line(1:width))
                line = .false.
                do row = 0, n_rows - 1
                    call partition(row, n_rows, height, rb, re)
                    do i = rb, re
                        call read_line(line,width)
                        do col = 0, n_cols - 1
                            call partition(col, n_cols, width, cb, ce)
                            dst = get_rank(row, col, n_rows, n_cols)
                            if (dst == root) then
                                old_world(i, cb:ce) = line(cb:ce)
                            else
                                call MPI_Send(line(cb), ce - cb + 1, MPI_LOGICAL, dst, 0, comm)
                            end if
                        end do
                    end do
                end do
                deallocate(line)
            end block
        else
            do i = ib, ie
                call MPI_Recv(old_world(i, jb), 1, a_row2, root, 0, comm, status)
            end do
        end if
        return
    end subroutine read_map

    !---------------------------------------------------------------------------

    subroutine print_line(logic_line, w)
        logical, dimension(:), pointer, intent(in) :: logic_line
        integer, intent(in) :: w
        character(len=:), allocatable :: line
        integer :: i

        allocate(character(len=w) :: line)
        do i = 1, w
            line(i:i) = merge('X', '.', logic_line(i))
        end do
        print "(a)", line
        if (allocated(line)) deallocate(line)
        return
    end subroutine print_line

    subroutine print_map()
        if (my_rank == root) then
            block
                logical, dimension(:), pointer :: line
                integer :: sender
                allocate(line(1:width))
                do row = 0, n_rows - 1
                    call partition(row, n_rows, height, rb, re)
                    do i = rb, re
                        do col = 0, n_cols - 1
                            call partition(col, n_cols, width, cb, ce)
                            sender = get_rank(row, col, n_rows, n_cols)
                            if (sender == root) then
                                line(cb:ce) = old_world(i, cb:ce)
                            else
                                call MPI_Recv(line(cb), ce - cb + 1, MPI_LOGICAL, sender, 0, comm, status)
                            end if
                        end do
                        call print_line(line,width)
                    end do
                end do
                print *
                deallocate(line)
            end block
        else
            do i = ib, ie
                call MPI_Send(old_world(i, jb), 1, a_row2, root, 0, comm)
            end do
        end if
        return
    end subroutine print_map

    !---------------------------------------------------------------------------

    ! Wait specified number of ms and then clear the terminal screen.
    subroutine wait_cls(ms)
        integer, intent(in) :: ms
        integer :: tick, tack
        real :: rate
        call system_clock(count=tick, count_rate=rate)
        do
            call system_clock(count=tack)
            if (real(tack - tick) / rate >= ms * 1e-3) exit
        end do
        ! Clear the terminal screen using console escape codes.
        print "(2a)", achar(27), '[2J'
        return
    end subroutine wait_cls

end program game_of_life