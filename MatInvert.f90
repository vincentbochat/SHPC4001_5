module MatInvert
    use, intrinsic :: omp_lib
    implicit none
    
contains

    function Matinverter(Mat) result(inverse)
        implicit none
        complex(8), intent(in) :: Mat(3, 3)
        integer(8) :: i, j
        complex(8):: inverse(3,3) ,detone, dettwo, det, current = 1.0


        !calculate det
        detone = Mat(1,1)*Mat(2,2)*Mat(3,3) + Mat(1,2)*Mat(2,3)*Mat(3,1) + Mat(1,3)*Mat(2,1)*Mat(3,2)
        dettwo = Mat(1,1)*Mat(2,3)*Mat(3,2) + Mat(1,2)*Mat(2,1)*Mat(3,3) + Mat(1,3)*Mat(2,2)*Mat(3,1)
        det = detone - dettwo
        
        !check if Mat is invertible
        if(abs(det) .eq. 0) then
            print*, 'Error: Matrix not invertable'
        else
            !set the corresponding Small determinents
    !$OMP parallel do
            do i=1, 3
                do j=1, 3
                    inverse(j, i) = current*smalldet(i, j, Mat)
                    current  = current*(-1.0)
                end do
            end do
    !$OMP end parallel do

        end if  
        
        !normalize by determinant
        inverse = inverse/det
    end function Matinverter

    function smalldet(i, j, Mat) result(small_det)
        implicit none
        integer(8) :: i, j, k, p, Fourcount
        complex(8) :: Mat(3,3), small_det
        complex(8), allocatable :: Mattemp(:,:), Small(:, :)
        allocate(Small(2, 2), Mattemp(3, 3))
        
        Mattemp = Mat + 0.0
        Fourcount = 1
        Small = 0.0
    
      do k= 1, 3
            do p = 1, 3      
                !ignore the ith row and jth column
                if (k .eq. i .or. p .eq. j) then
                else
                    ! update the Small Mat values depending on success count value
                    if(Fourcount .eq. 1) then
                            Small(1,1) = Mattemp(k,p)
                    else if(Fourcount .eq. 2) then
                            Small(1,2) = Mattemp(k,p)
                    else if(Fourcount .eq. 3) then
                            Small(2,1) = Mattemp(k,p)
                    else if(Fourcount .eq. 4) then
                            Small(2,2) = Mattemp(k,p)
                    end if
                    !update count
                    Fourcount = Fourcount + 1
                end if
            end do
        end do
    
        small_det = Small(1,1)*Small(2,2) - Small(1,2)*Small(2,1)
        ! clear 
        deallocate(Mattemp, Small)
        
        end function smalldet

end module MatInvert