
module CreateHerm
    use, intrinsic :: omp_lib

    implicit none
    
contains

    subroutine HermSq(N, L)
        complex(8), allocatable, intent(inout) :: L(:, :)
        Integer, intent(in) :: N
        Integer :: i
        complex(8) :: k1, k2
    
        k1 = 2.0    !main diagonal
        k2 = -1.0   !second diagonal
    
        Allocate(L(N, N))
    
        L = 0.0 !Assign 0 to entire matrix
   
        do i = 1, N
            L(i, i) = k1
            if (i .lt. N) then
                L(i, i + 1) = k2
                L(i + 1, i) = k2
            endif
        enddo
            
    end subroutine HermSq
       
end module CreateHerm