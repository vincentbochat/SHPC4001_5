program Ass5_Q2_b

    use CreateHerm
    implicit none
    complex(8), allocatable :: L(:, :), vec(:,:), WORK(:)
    real(8), allocatable :: ev(:), RWORK(:)
    Integer :: N = 4
    Integer :: i, k, info   !defined in Herm
    complex(8) :: k1, k2    !Defined in Herm

    write(*,*) "What is the dimension of your square vector? - N"
    read(*,*) N

    allocate(ev(N), vec(N, N), WORK(N), RWORK(3*N - 2))

    call Herm(N = N, L = L)

    call ZHBEV('V', 'U', N, 1, L, 2, ev, vec, N, WORK, RWORK, info)

    if (info .ne. 0) then
        write(*,*) "Error"
    endif

    write(*,*) ev

end program Ass5_Q2_b


