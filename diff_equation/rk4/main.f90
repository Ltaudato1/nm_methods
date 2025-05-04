program rk4_adaptive_study
    implicit none
    integer, parameter :: p = 4
    integer :: j, i, n, iter, max_points
    real(8), parameter :: a = 1.0d0, b = 3.0d0, h_min = 1d-14
    real(8), dimension(13) :: eps_list, h_list, err_list, iter_list
    real(8), allocatable :: x_vals(:), y_vals(:)
    real(8) :: x, y, h, eps, max_err, y_full, y_half, err, exact

    eps_list = [1d-1, 1d-2, 1d-3, 1d-4, 1d-5, 1d-6, 1d-7, 1d-8, 1d-9, 1d-10, 1d-11, 1d-12, 1d-13]

    write(*,*) ' eps         max_error       h_last        iters'
    do j = 1, size(eps_list)
        eps = eps_list(j)
        h = 0.1d0
        x = a
        y = f_exact(a)  ! Начальное условие
        max_points = 10000
        allocate(x_vals(0:max_points), y_vals(0:max_points))
        x_vals(0) = x
        y_vals(0) = y
        iter = 0
        n = 0

        do while (x < b .and. h >= h_min)
            iter = iter + 1

            call rk4_single(x, y, h, y_full)
            call rk4_double_step(x, y, h, y_half)

            err = abs(y_half - y_full) / (2.0d0**p - 1.0d0)

            if (err > eps) then
                h = h / 2.0d0
            else
                x = x + h
                y = y_half
                n = n + 1
                x_vals(n) = x
                y_vals(n) = y
            end if
        end do

        max_err = 0.0d0
        do i = 0, n
            exact = f_exact(x_vals(n))
            err = abs(exact - y_vals(n))
            if (err > max_err) max_err = err
        end do

        h_list(j) = h
        err_list(j) = max_err
        iter_list(j) = iter

        write(*,'(1PE10.1, 3X, 1PE10.1, 3X, 1PE10.2, 3X, I6)') eps, max_err, h, iter
        deallocate(x_vals, y_vals)
    end do

contains
    subroutine rk4_double_step(x, y, h, y_out)
        real(8), intent(in) :: x, y, h
        real(8), intent(out) :: y_out
        real(8) :: y_temp
        call rk4_single(x, y, h/2.0d0, y_temp)
        call rk4_single(x + h/2.0d0, y_temp, h/2.0d0, y_out)
    end subroutine rk4_double_step

    function f(x, y) result(res)
        real(8), intent(in) :: x, y
        real(8) :: res
        res = exp(x) / x + y 
    end function f

    real(8) function f_exact(x) result(res)
        real(8), intent(in) :: x
        res = exp(x) * (log(abs(x)) + 1.0d0)
    end function f_exact

    subroutine rk4_single(x, y, h, y_out)
        real(8), intent(in) :: x, y, h
        real(8), intent(out) :: y_out
        real(8) :: k1, k2, k3, k4
        k1 = f(x, y)
        k2 = f(x + h/3.d0, y + h*k1/3.d0)
        k3 = f(x + 2.d0*h/3.d0, y - h*k1/3.d0 + h*k2)
        k4 = f(x + h, y + h*k1 - h*k2 + h*k3)
        y_out = y + h * (k1 + 3.d0*k2 + 3.d0*k3 + k4) / 8.d0  ! Метод РК4
    end subroutine rk4_single

end program rk4_adaptive_study