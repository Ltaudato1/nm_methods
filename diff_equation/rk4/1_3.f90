program rk4_method
    implicit none
    integer, parameter :: n_max = 10000 
    integer, parameter :: m = 10        
    integer, parameter :: Nx = 20       
    real(8), parameter :: a = 1.0d0, b = 3.0d0
    real(8) :: h0, h, eps
    real(8) :: y_h, y_h2, err_runge, y0, y_prev, error, cur_error
    real(8), dimension(Nx) :: x_list
    integer :: i, j, iter, k

    h0 = (b - a) / real(Nx - 1)
    do k = 1, Nx
        x_list(k) = a + (k - 1) * h0
    end do
    y0 = f_exact(a)
    open(unit=10, file='results.csv', status='replace', action='write')
    write(10, '(A)') 'j,eps,h_final,error,iters'
    iter = 1
    do j = 1, m
        eps = 10.0d0**(-j) 
        y_prev = y0
        error = 0
        do i = 1, Nx - 1
            do
                h = h0 / (2.0d0**(iter-1))
                y_h  = integrate_rk4(x_list(i), y_prev, h, x_list(i + 1))
                y_h2 = integrate_rk4(x_list(i), y_prev, h/2.0d0, x_list(i + 1))
                err_runge = abs(y_h2 - y_h) / (2.0d0**4 - 1.0d0)
                if (err_runge <= eps) exit
                if (iter == n_max) exit
                iter = iter + 1
            end do
            y_prev = y_h
            cur_error = abs(y_h2 - f_exact(x_list(i + 1)))
            if (error < cur_error) error = cur_error
        end do
    write(10,'(I0,",",1P,E10.3,",",1P,E10.3,",",1P,E10.3,",",I0)') &
             j, eps, h, error, iter
    end do
    close(10)

contains

    function f(x, y) result(res)
        implicit none
        real(8), intent(in) :: x, y
        real(8) :: res
        res = dexp(x)/x + y
    end function f

    function f_exact(x) result(res)
        implicit none
        real(8), intent(in) :: x
        real(8) :: res
        res = dexp(x)*(log(abs(x)) + 1.d0)
    end function f_exact

    function rk4_step(x, y, h) result(y_next)
        implicit none
        real(8), intent(in) :: x, y, h
        real(8) :: y_next, k1, k2, k3, k4
        k1 = f(x, y)
        k2 = f(x + h/3.d0, y + h*k1/3.d0)
        k3 = f(x + 2.d0*h/3.d0, y - h*k1/3.d0 + h*k2)
        k4 = f(x + h, y + h*k1 - h*k2 + h*k3)
        y_next = y + h*(k1 + 3.d0*k2 + 3.d0*k3 + k4)/8.0d0
    end function rk4_step

    function integrate_rk4(x0, y0, h, x_end) result(y_end)
        implicit none
        real(8), intent(in) :: x0, y0, h, x_end
        real(8) :: y_end, x_loc
        integer :: m_steps, kk
        m_steps = max(1, int((x_end - x0)/h))
        y_end = y0
        x_loc = x0
        do kk = 1, m_steps
            y_end = rk4_step(x_loc, y_end, h)
            x_loc = x_loc + h
        end do
        if (x_loc < x_end) y_end = rk4_step(x_loc, y_end, x_end - x_loc)
    end function integrate_rk4

end program rk4_method