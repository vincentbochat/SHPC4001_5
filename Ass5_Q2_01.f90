program Ass5_Q2_01

    use CreateHerm
    implicit none
    complex(8), allocatable :: L(:, :)
    Integer :: N, j
    Integer :: i, k
    complex(8) :: k1, k2

    N = 4

    call HermSq(N = 4, L = L)

    do j = 1, N
        write(*,*) L(j, :)
    end do
 
end program Ass5_Q2_01


