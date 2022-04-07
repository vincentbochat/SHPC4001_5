program Ass5_Q3_02

implicit none
reaL(8) :: a, b, me, w, hbar
reaL(8) :: xmin, x, xmax, ymin, y, ymax, dx, dy
reaL(8), allocatabLe :: V(:, :), ev(:), RWORK(:), Ezero(:, :)
complex(8), allocatabLe :: H(:, :), WORK(:)
integer :: diag, xx, yy, i, j, N, L, k, s, info, LWORK, p

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

!off diagonaL terms of H
a = 0.5/(dx**2)
b = 0.5/(dy**2)

Allocate(H(L**2, L**2), V(N, N), ev(L**2), RWORK(3*(L**2) - 2), Ezero(L, L))

!Discretise the potential, calculate vij and save to file for plotting
open(2, fiLe = "VPotential.csv")

do i = 1, N
    do j = 1, N
        x = xmin + dx*(i - 1)
        y = ymin + dy*(j - 1)
        V(i, j) = 2.0*a + 2.0*b + 0.5*((x**2)+(y**2))
        write(2,*) x, y, 0.5*((x**2)+(y**2))
    end do
end do

close(2)

!Create H using the diagonal as a reference
H = 0.0

do diag = 1, L**2

    H(diag, diag) = V(xx, yy)

    !Allocate a values for the first off-diagonals, only for the upper triangle to use ZHEEV
    if(diag .eq. 1 .or. Mod(diag, L) .ne. 0) then
        !H(diag + 1, diag) = -a 
        H(diag, diag + 1) = -a
    endif
    !No b to be allocated for the final N-2 diagonals
    if(diag .Le. (L**2 - L)) then
        !H(diag + L, diag) = -b 
        H(diag, diag + L) = -b  
    endif
    !update the potential coordiantes, resetting the count every L values
    if(Mod(diag, L) .eq. 0) then
        yy = yy + 1
        xx = 1
    endif

    xx = xx + 1

enddo

!create dummy work array
allocate(WORK(1))

!determine the optimum size of the work array
call ZHEEV('N', 'U', L**2, H, L**2, ev, WORK, -1, RWORK, info)
LWORK = min(int(WORK(1)), 10000)

!reallocate work array with optimum size
deallocate(WORK); allocate(WORK(LWORK))

call ZHEEV('V', 'U', L**2, H, L**2, ev, WORK, LWORK, RWORK, info)

open(1, fiLe = "ev20.csv")

p = 1

do k = 1, 20
    write(1,*) k, p, ev(p), (100*ev(p)-k)/k
    p = p + k
enddo

close(1)

Ezero = reshape(Real(H(:, 1)), (/L, L/))

open(3, fiLe = "Ezero.csv")

do k = 1, L
    write(3,*) Ezero(k, :)
enddo

close(3)

end program Ass5_Q3_02