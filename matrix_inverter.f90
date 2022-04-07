module matrix_inverter
    use, intrinsic :: omp_lib
    implicit none
    
contains

function calculate_minor_det(i,j,matrix) result(minor_det)
    implicit none
    integer(8) :: i,j,k,m,succ_count
    complex(8) :: matrix(3,3),minor_det
    complex(8), allocatable :: matrix_copy(:,:),minor(:,:)
    allocate(minor(2,2), matrix_copy(3,3))
    

    matrix_copy = matrix + 0.0
    succ_count = 1
    minor =minor*0

  do k=1,3
        do m = 1,3
          
            !ignore the ith row and jth column
            if (k == i .or. m == j) then
            else
                ! update the minor matrix values depending on success count value
                if(succ_count == 1) then
                        minor(1,1) = matrix_copy(k,m)
                else if(succ_count == 2) then
                        minor(1,2) = matrix_copy(k,m)
                else if(succ_count == 3) then
                        minor(2,1) = matrix_copy(k,m)
                else if(succ_count == 4) then
                        minor(2,2) = matrix_copy(k,m)
                end if

                !update success count
                succ_count = succ_count + 1

            end if

          

       
        end do

    end do

    minor_det = minor(1,1)*minor(2,2) - minor(1,2)*minor(2,1)
    ! clear memory
    deallocate(matrix_copy,minor)



    
    end function calculate_minor_det
    function mat_3_inverter(matrix) result(inverse)
        implicit none
        complex(8), intent(in) :: matrix(3,3)
        integer(8) :: i ,j
        complex(8):: inverse(3,3),det_pos,det_neg,det,current_identifier=1.0


        !calculate det
        det_pos = matrix(1,1)*matrix(2,2)*matrix(3,3) + matrix(1,2)*matrix(2,3)*matrix(3,1) + matrix(1,3)*matrix(2,1)*matrix(3,2)
        det_neg = matrix(1,1)*matrix(2,3)*matrix(3,2)+matrix(1,2)*matrix(2,1)*matrix(3,3)+matrix(1,3)*matrix(2,2)*matrix(3,1)
        det = det_pos-det_neg
        
        !check if matrix is invertible
        if(abs(det) == 0) then
            print*, 'ERROR: MATRIX NOT INVERTIBLE'
        else
            !set the corresponding minor determinents
    !$OMP parallel do
            do i=1,3
                do j=1,3
                    inverse(j,i) = current_identifier*calculate_minor_det(i,j,matrix)
                    current_identifier = current_identifier*(-1.0)
                end do
            end do
    !$OMP end parallel do


        end if  
        
        !normalize by determinant
        inverse = inverse/det
    end function mat_3_inverter

end module matrix_inverter