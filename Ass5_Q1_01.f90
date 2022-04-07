program Ass5_Q1_01
    use Matinvert
    implicit none
    complex(8) :: Mat(3,3), inverted_Mat(3,3), EigVal
    interface 
    subroutine PowerItSys(Mat)
        implicit none
        complex(8) :: Mat(3,3), EigVal, prev = 10, vectormult(3), initial = 10
        end subroutine
    end interface

    Mat = reshape((/ (4,0), (0,1), (2,0), (0,-1), (2,0), (2,-7), (2,0), (2,7), (-2,0) /),shape(Mat))

    inverted_Mat = Matinverter(Mat)  

   call PowerItSys(inverted_Mat)
    
end program Ass5_Q1_01


function rayleigh_quotient(vectormult,Mat) result(quot)
    implicit none
    complex(8), intent(in) :: Mat(3, 3), vectormult(3)
    complex(8), allocatable :: quot

    !calculate eigen value through rayleigh method
    quot = dot_product(vectormult, matmul(Mat, vectormult))
    quot = quot / dot_product(vectormult, vectormult)
    
end function rayleigh_quotient

subroutine PowerItSys(Mat)
    implicit none
    complex(8) :: Mat(3,3), EigVal, initial=10, prev=200, vectormult(3)
    integer(8) :: count
    
    interface   
    function rayleigh_quotient(vectormult,Mat) result(quot)
    implicit none
    complex(8), intent(in) :: Mat(3,3),vectormult(3)
    complex(8), allocatable :: quot
    end function
    end interface

    !set up vectormult Mat for iteration method
    vectormult = reshape((/ 1.0,1.0,1.0 /), shape(vectormult))
    open(1, file='EVQ1.csv')

    do while(abs(prev-initial) .gt. 1e-7)   !set desired precision
        prev = initial
        vectormult = MATMUL(Mat, vectormult)
        vectormult = vectormult/norm2(abs(vectormult))
        !check for eigen value
        initial = rayleigh_quotient(vectormult, Mat)
        
        !update count
        count = count + 1
        write(1,*) real(initial), count !save the eigenvalues to file

    end do

    close(1)

    print*,'Dominant eigenvector =', vectormult
    print*,'Dominant eigenvalue =', initial

end subroutine PowerItSys


