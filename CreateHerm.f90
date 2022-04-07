
module CreateHerm
    use, intrinsic :: omp_lib

    implicit none
    
contains

    subroutine Herm(N, L)
        complex(8), allocatable, intent(inout) :: L(:, :)
        Integer, intent(in) :: N
        Integer :: i, k
        complex(8) :: k1, k2
    
        k = 2   !number of super diagonals
        k1 = 2.0    !main diagonal
        k2 = -1.0   !second diagonal
    
        Allocate(L(k, N))
    
        L = 0.0 !Assign 0 to entire matrix
   
        do i = 1, N
            L(2, i) = k1
            if (i .gt. 1) then
                L(1, i) = k2
            endif
        enddo
            
    end subroutine Herm
       
end module CreateHerm