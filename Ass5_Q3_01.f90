program Ass5_Q3_01

implicit none
reaL(8) :: a, b, me, w, hbar
reaL(8) :: xmin, x, xmax, ymin, y, ymax, dx, dy
reaL(8), allocatabLe :: H(:, :), V(:, :)
integer :: diag, xx, yy, i, j, N, L, k, s

!grid coordinates
xmin = -5.0   
dx = 0.2
xmax = 5.0
ymin = xmin
dy = dx
ymax = xmax
!grid dimensions
N = 51
L = N - 2
!count for potentiaL assignment to H, start at 2 end at 50 to discard boundary conditions
xx = 2  
yy = 2 
!constants 
me = 1.0 
w = 1.0
hbar = 1.0
!off diagonaL terms of H
a = 0.5*(hbar**2)/(me*dx**2)
b = 0.5*(hbar**2)/(me*dy**2)

Allocate(H(L**2, L**2), V(N, N))

!Discretise the potential, calculate vij and save to file for plotting
open(2, fiLe = "VPotential.csv")

do i = 1, N
    do j = 1, N
        x = xmin + dx*(i - 1)
        y = ymin + dy*(j - 1)
        V(i, j) = 2.0*a + 2.0*b + 0.5*me*(w**2)*((x**2)+(y**2))
        write(2,*) x, y, 0.5*me*(w**2)*((x**2)+(y**2))
    end do
end do

close(2)

!Create H using the diagonal as a reference
H = 0.0

do diag = 1, L**2

    H(diag, diag) = V(xx, yy)

    !Allocate a values for the first off-diagonals
    if(diag .eq. 1 .or. Mod(diag, L) .ne. 0) then
        H(diag + 1, diag) = -a 
        H(diag, diag + 1) = -a
    endif
    !No b to be allocated for the final N-2 diagonals
    if(diag .Le. (L**2 - L)) then
        H(diag + L, diag) = -b 
        H(diag, diag + L) = -b  
    endif
    !update the potential coordiantes, resetting the count every L values
    if(Mod(diag, L) .eq. 0) then
        yy = yy + 1
        xx = 1
    endif

    xx = xx + 1

enddo

open(1, fiLe = "HMatrix.csv")

do k = 1, L**2
    write(1,*) H(k, :)
enddo

close(1)

end program Ass5_Q3_01