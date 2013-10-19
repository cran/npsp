
!-----------------------------------------------------------------------
    module grid_module
!-----------------------------------------------------------------------
!   Modulo clases rejilla (grid) y rejilla binning (grid_bin)
!       type(grid) :: g
!           g%ndim          = N� de dimensiones
!           g%n(1:g%ndim)   = N� de nodos en cada dimensi�n
!           g%ngrid         = N� total de nodos (= PRODUCT(n) rejilla est�ndar)
!           g%min(1:g%ndim) = M�nimo de la rejilla binning en cada dimensi�n
!           g%max(1:g%ndim) = M�ximo      ""
!           g%lag(1:g%ndim) = Espaciado   ""
!           g%i             = indice unidimensional
!           g%ii(1:g%ndim)  = indice multidimensional
!           procedure :: set => set_grid, end => end_grid, ind, set_ind, incii
!       type(grid_den) :: bin
!           bin%w(1:bin%ngrid)  = Peso/frecuencia nodo binning (unidim.)
!           bin%ny              = N� de datos (suma de los pesos binning)
!           procedure :: set_bin_den, end_bin_den
!       type(grid_bin) :: bin
!           bin%y(1:bin%ngrid)  = Valor nodo binning (indice unidimensional)
!           bin%med             = Media ponderada de bin%y
!                                 (media de los datos & binning; para estimaci�n aditiva)
!           procedure :: set_bin => set_grid_bin, end_bin => end_grid_bin
!       set_grid1d
!
!   Interfaces con R:
!       binning             Rejilla binning ("binning" en "bin.data.R")
!       interp_data_grid    Interpolacion lineal de una rejilla ("interp.grid.par" en "interp.R")
!
!   PENDENTE:
!       - Crear clase grid_den intermedia entre grid y grid_bin
!         (renombrar %ny %nw? %nx?)  
!       - Implementar ndim como par�metro LEN de type grid
!       - Crear clase grid_data (intermedia entre grid y grid_bin?)
!       - procedure :: end type-bound / final
!
!   Autor: (c) Ruben Fernandez-Casal    Ultima revision: Oct 2012  Nov 2013
!-----------------------------------------------------------------------
        implicit none

!       ----------------------------------------------------------------
        type grid
!       ----------------------------------------------------------------
!       Indice de una rejilla multidimensional
!       ----------------------------------------------------------------
            integer ndim, ngrid, igrid
            integer, allocatable :: n(:), ii(:)
!           n(g%ndim) = NBM(1..NDimBM)= N� nodos de la rejilla binning
            real*8, allocatable :: min(:), max(:), lag(:)
        contains
            procedure :: set => set_grid
!           final :: end_grid
            procedure :: end => end_grid
            procedure :: ind
            procedure :: set_ind
            procedure :: incii
        end type

!       ----------------------------------------------------------------
        type, extends(grid) :: grid_den
!       ----------------------------------------------------------------
!       Rejilla binning multidimensional densidad
!           bin%w(1:bin%ngrid) = Peso/frecuencia nodo binning
!       ----------------------------------------------------------------
            integer ny    !NBMc
            real*8, allocatable :: w(:)
        contains
            procedure :: set_bin_den
!           final :: end_bin_den
            procedure :: end_bin_den
        end type


!       ----------------------------------------------------------------
        type, extends(grid_den) :: grid_bin
!       ----------------------------------------------------------------
!       Rejilla binning multidimensional datos
!           bin%y(1:bin%ngrid) = Valor nodo binning
!           bin%med = Media ponderada de BMy (media de los datos binning)
!       ----------------------------------------------------------------
            real*8 :: med
            real*8, allocatable :: y(:)
        contains
            procedure :: set_bin => set_grid_bin
!           final :: end_grid_bin
            procedure :: end_bin => end_grid_bin
        end type


!   --------------------------------------------------------------------
    contains
!   --------------------------------------------------------------------


!       ----------------------------------------------------------------
        subroutine set_grid(g, ndim, n, min, max)
!       ----------------------------------------------------------------
!       Establece la rejilla
!       ----------------------------------------------------------------
        implicit none
        class(grid) :: g
        integer ndim, n(ndim)
        real*8 min(ndim), max(ndim)
!       ----------------------------------------------------------------
            g%ndim = ndim
            allocate(g%n(ndim), g%ii(ndim), g%min(ndim), g%max(ndim), g%lag(ndim))
            g%n = n
            g%ngrid = PRODUCT(n)
            g%min = min
            g%max = max
            g%lag = (max - min)/(n - 1.0d0)
        return
        end subroutine set_grid


!       ----------------------------------------------------------------
        subroutine set_grid1d(g, n, min, max)
!       Establece una rejilla unidimensional
!       Avoid rank mismatch in argument(rank-1 and scalar)
!       ----------------------------------------------------------------
        implicit none
        class(grid) :: g
        integer n
        real*8 min, max
!       ----------------------------------------------------------------
            g%ndim = 1
            allocate(g%n(1), g%ii(1), g%min(1), g%max(1), g%lag(1))
            g%n(1) = n
            g%ngrid = n
            g%min(1) = min
            g%max(1) = max
            g%lag(1) = (max - min)/(n - 1.0d0)
        return
        end subroutine set_grid1d


!       ----------------------------------------------------------------
        subroutine end_grid(g)
!       Libera memoria
!       ----------------------------------------------------------------
        implicit none
        class(grid) :: g
!       ----------------------------------------------------------------
            deallocate(g%n, g%ii, g%min, g%max, g%lag)
        return
        end subroutine end_grid


!       ----------------------------------------------------------------
        integer function ind(g, ii)
!       Devuelve el indice unidimensional equivalente a uno
!       multidimensional. Por ejemplo, en el caso tridimensional:
!       iibm  = ix + (iy-1)*nbmx + (iz-1)*nbmx*nbmy
!             = ix + ((iy-1) + (iz-1)*nbmy)*nbmx
!       (se evalua empleando una regla tipo Horner)
!       ----------------------------------------------------------------
        implicit none
        class(grid) :: g
        integer ii(g%ndim)
        integer i, k
!       ----------------------------------------------------------------
            i = 0
            do k = g%ndim, 2, -1 
                i = g%n(k-1)*(i + ii(k) - 1)
            end do
            ind = i + ii(1)
        return
        end function ind


!       ----------------------------------------------------------------
        subroutine set_ind(g, ii)
!       Establece el indice unidimensional y multidimensional, a partir
!       de un indice multidimensional
!       PENDIENTE: VERIFICAR RANGO INDICES
!       ----------------------------------------------------------------
        implicit none
        class(grid) :: g
        integer ii(g%ndim)
            g%ii = ii
            g%igrid = ind(g, ii)
        return
        end subroutine set_ind

!       ----------------------------------------------------------------
        subroutine incii(g)
!       Incrementa el indice (uni y multi dimensional)
!       NO VERIFICA RANGOS
!       ----------------------------------------------------------------
        implicit none
        class(grid) :: g
        integer j
!           ------------------------------------------------------------
            g%igrid = g%igrid+1
            do j = 1, g%ndim
                g%ii(j) = g%ii(j)+1
                if (g%ii(j) <= g%n(j)) return
                g%ii(j) = 1
            end do
            return
        end subroutine incii


!       ----------------------------------------------------------------
        subroutine set_bin_den(g, nd, nbin, x, ny)
!       Establece la rejilla binning (lineal) para densidad
!       ----------------------------------------------------------------
        implicit none
        class(grid_den) :: g
        integer nd, nbin(nd), ny
        real*8  x(nd,ny)
!
        integer ii(nd), iib(nd), ib, i, j, k
        real*8  minx(nd), maxx(nd), w(2,nd), tmp
        integer niinc, iinc(nd, 2**nd)
!       real*8, external :: dmach
!           ------------------------------------------------------------
!           Calcular dimensiones rejilla binning
            minx = x(1:nd, 1)
            maxx = minx
            do i = 2, ny
                do j = 1, nd
                    if (minx(j)>x(j, i)) then
                        minx(j) = x(j, i)
                    else if (maxx(j)<x(j, i)) then
                        maxx(j) = x(j, i)
                    end if
                end do
            end do
!           Expandir un poco y establecer rejilla
            tmp = 1.0d2*epsilon(1.0d0)
            minx = minx-tmp
            maxx = maxx+tmp
            call set_grid(g, nd, nbin, minx, maxx)
!           Asignar memoria rejilla binning
            allocate(g%w(g%ngrid))
            g%ny = ny
            g%w = 0.0d0
!           Indice rejilla actualizaci�n
            niinc = 2**nd
            ii = 0
            do i = 1, niinc
                do j = 1, nd-1
                    if (ii(j) > 1) then
                        ii(j) = 0
                        ii(j+1) = ii(j+1)+1
                    else
                        exit
                    end if
                end do
                iinc(1:nd, i) = ii(1:nd)
                ii(1) = ii(1) + 1
            end do
!           Recorrer datos
            do i = 1, ny
                do j = 1, nd
                    iib(j) = 1 + int((x(j, i) - g%min(j)) / g%lag(j))
!                   calculo de los pesos
                    w(2, j) = (x(j, i) - g%min(j) - (iib(j)-1)*g%lag(j)) / g%lag(j)
                    w(1, j) = 1.0d0 - w(2, j)
                end do
!               Actualizar valores
                do k = 1, niinc
                    tmp = 1.0d0
                    do j = 1, nd
                        ii(j) = iib(j) + iinc(j, k)
                        tmp = tmp * w(iinc(j, k) + 1, j)
                    end do
                    ! if (tmp < epsilon) cycle
                    ib = g%ind(ii)
                    g%w(ib) = g%w(ib) + tmp
                end do
            end do
        return
        end subroutine set_bin_den


!       ----------------------------------------------------------------
        subroutine end_bin_den(g)
!       Libera memoria
!       ----------------------------------------------------------------
        implicit none
        class(grid_den) :: g
            call g%end
            deallocate(g%w)
        return
        end subroutine end_bin_den

        
!       ----------------------------------------------------------------
        subroutine set_grid_bin(g, nd, nbin, x, ny, y)
!       Establece la rejilla binning (lineal)
!       ----------------------------------------------------------------
        implicit none
        class(grid_bin) :: g
        integer nd, nbin(nd), ny
        real*8  x(nd,ny), y(ny)
!
        integer ii(nd), iib(nd), ib, i, j, k
        real*8  minx(nd), maxx(nd), w(2,nd), tmp
        integer niinc, iinc(nd, 2**nd)
!       real*8, external :: dmach
!           ------------------------------------------------------------
!           Calcular dimensiones rejilla binning
            minx = x(1:nd, 1)
            maxx = minx
            do i = 2, ny
                do j = 1, nd
                    if (minx(j) > x(j, i)) then
                        minx(j) = x(j, i)
                    else if (maxx(j) < x(j, i)) then
                        maxx(j) = x(j, i)
                    end if
                end do
            end do
!           Expandir un poco y establecer rejilla
            tmp = 1.0d2*epsilon(1.0d0)
            minx = minx - tmp
            maxx = maxx + tmp
            call set_grid(g, nd, nbin, minx, maxx)
!           Asignar memoria rejilla binning
            allocate(g%y(g%ngrid), g%w(g%ngrid))
            g%ny = ny
            g%y = 0.0d0
            g%w = 0.0d0
!           Indice rejilla actualizaci�n
            niinc = 2**nd
            ii = 0
            do i = 1, niinc
                do j = 1, nd-1
                    if (ii(j) > 1) then
                        ii(j) = 0
                        ii(j+1) = ii(j+1)+1
                    else
                        exit
                    end if
                end do
                iinc(1:nd, i) = ii(1:nd)
                ii(1) = ii(1) + 1
            end do
!           Recorrer datos
            do i = 1, ny
                do j = 1, nd
                    iib(j) = 1 + int((x(j, i) - g%min(j)) / g%lag(j))
!                   calculo de los pesos
                    w(2, j) = (x(j, i) - g%min(j) - (iib(j)-1)*g%lag(j)) / g%lag(j)
                    w(1, j) = 1.0d0 - w(2, j)
                end do
!               Actualizar valores
                do k = 1, niinc
                    tmp = 1.0d0
                    do j = 1, nd
                        ii(j) = iib(j) + iinc(j, k)
                        tmp = tmp * w(iinc(j, k) + 1, j)
                    end do
                    ! if (tmp < epsilon) cycle
                    ib = g%ind(ii)
                    g%y(ib) = g%y(ib) + tmp * y(i)
                    g%w(ib) = g%w(ib) + tmp
                end do
            end do
!           promediar y calcular valor medio
            g%med = 0.0d0
            do i = 1, g%ngrid
                if (g%w(i) > 0.d0) then
                    g%med = g%med + g%y(i)/dble(g%ny)
                    g%y(i) = g%y(i) / g%w(i)
                end if
            end do
        return
        end subroutine set_grid_bin


!       ----------------------------------------------------------------
        subroutine end_grid_bin(g)
!       Libera memoria
!       ----------------------------------------------------------------
        implicit none
        class(grid_bin) :: g
            call g%end
            deallocate(g%y, g%w)
        return
        end subroutine end_grid_bin
        
!   --------------------------------------------------------------------
    end module grid_module
!   --------------------------------------------------------------------



!   --------------------------------------------------------------------
    subroutine binning(nd, nbin, x, ny, y, bin_min, bin_max, bin_med, bin_y, bin_w)
!   --------------------------------------------------------------------
!       Interfaz para la rutina de R "binning"
!       Devuelve la rejilla binning type(grid_bin)%set_bin(nd, nbin, x, ny, y)
!
!   Autor: (c) Ruben Fernandez-Casal    Ultima revision: Jun 2012
!   --------------------------------------------------------------------
    use grid_module
    implicit none
    integer nd, nbin(nd), ny
    real*8  x(nd,ny), y(ny)
    real*8  bin_min(nd), bin_max(nd), bin_med, bin_y(*), bin_w(*)
    type(grid_bin) :: bin
        call bin%set_bin(nd, nbin, x, ny, y) ! Establece la rejilla binning (lineal)
        bin_min = bin%min
        bin_max = bin%max
        bin_med = bin%med
        bin_y(1:bin%ngrid) = bin%y
        bin_w(1:bin%ngrid) = bin%w
        call bin%end_bin
    return
    end subroutine binning


!   --------------------------------------------------------------------
    subroutine bin_den(nd, nbin, x, ny, bin_min, bin_max, bin_w)
!   --------------------------------------------------------------------
!       Interfaz para la rutina de R "bin.den"
!       Devuelve la rejilla binning type(grid_den)%set_bin_den(nd, nbin, x, ny, y)
!
!   Autor: (c) Ruben Fernandez-Casal    Ultima revision: Nov 2013
!   --------------------------------------------------------------------
    use grid_module
    implicit none
    integer nd, nbin(nd), ny
    real*8  x(nd,ny)
    real*8  bin_min(nd), bin_max(nd), bin_w(*)
    type(grid_den) :: bin
        call bin%set_bin_den(nd, nbin, x, ny) ! Establece la rejilla binning (lineal)
        bin_min = bin%min
        bin_max = bin%max
        bin_w(1:bin%ngrid) = bin%w
        call bin%end_bin_den
    return
    end subroutine bin_den

!   --------------------------------------------------------------------
    subroutine interp_data_grid( nd, nbin, bin_min, bin_max,     &
   &                    ngrid, gy, x, ny, y)
!   --------------------------------------------------------------------
!       Interfaz para la rutina de R "interp.data.grid"
!       Interpolacion lineal de una rejilla
!       Establece type(grid) :: g a partir de par�metros
!   --------------------------------------------------------------------
!       PENDENTE: FALLA CON DATOS MISSING
!   ----------------------------------------------------------------
    use grid_module
    implicit none
    integer nd, nbin(nd), ngrid, ny
    real*8  bin_min(nd), bin_max(nd)
    type(grid) :: g
    real*8  x(nd,ny), gy(ngrid), y(ny)
!       ----------------------------------------------------------------
!       Establecer la rejilla
        call g%set(nd, nbin, bin_min, bin_max)
!       Interpolaci�n
        call interp_grid(g, gy, x, ny, y)
        call g%end
    return
    end subroutine interp_data_grid



!   ----------------------------------------------------------------
    subroutine interp_grid(g, gy, x, ny, y)
!       Interpola linealmente valores en una rejilla
!   ----------------------------------------------------------------
!       PENDENTE: FALLA CON DATOS MISSING
!   ----------------------------------------------------------------
    use grid_module
    implicit none
    type(grid) :: g
    integer ny
    real*8  x(g%ndim, ny), gy(g%ngrid), y(ny)
!
    integer nd, ii(g%ndim), iib(g%ndim), ib, i, j, k
    real*8  w(2,g%ndim), tmp
    integer niinc, iinc(g%ndim, 2**g%ndim)
!       ------------------------------------------------------------
        nd = g%ndim
        niinc = 2**nd
!       Indice rejilla actualizaci�n
        ii = 0
        do i = 1, niinc
            do j = 1, nd-1
                if (ii(j) > 1) then
                    ii(j) = 0
                    ii(j+1) = ii(j+1)+1
                else
                    exit
                end if
            end do
            iinc(1:nd, i) = ii(1:nd)
            ii(1) = ii(1) + 1
        end do
!       Recorrer datos
        y = 0.0d0
        do i = 1, ny
            do j = 1, nd
                iib(j) = 1 + int((x(j, i) - g%min(j)) / g%lag(j))
!               Extrapolaci�n:
                if (iib(j) < 1) iib(j) = 1
                if (iib(j) >= g%n(j)) iib(j) = g%n(j) - 1
!               Calculo de los pesos
                w(2, j) = (x(j, i) - g%min(j) - (iib(j)-1)*g%lag(j)) / g%lag(j)
                w(1, j) = 1.0d0 - w(2, j)
            end do
!           Actualizar valores
            do k = 1, niinc
                tmp = 1.0d0
                do j = 1, nd
                    ii(j) = iib(j) + iinc(j, k)
                    tmp = tmp * w(iinc(j, k) + 1, j)
                end do
                ib = g%ind(ii)
                y(i) = y(i) + tmp * gy(ib)
            end do
        end do
    return
    end subroutine interp_grid





!   --------------------------------------------------------------------
!   NOTAS SOBRE VARIABLES PARA CONVERSI�N DE C�DIGO SKUET
!       type(grid_bin) :: bin
!           bin%min             !   MinBM(1..NDimBM)= M�nimo de la rejilla binning
!                                                       en cada dimensi�n
!           bin%max             !   MaxBM(1..NDimBM)= M�ximo ""
!           bin%lag             !   LBM(1..NDimBM)= Espaciado ""
!           bin%y(1:bin%ngrid)  !   BMy(1..NTotBM)= Valor nodo binning
!                                                       (indice unidimensional)
!           bin%w(1:bin%ngrid)  !   BMc(1..NTotBM)= Peso/frecuencia nodo binning
!                                                       (indice unidimensional)
!           bin%med             !   MedBMy= Media ponderada de BMy
!                                                       (media de los datos binning)
!       iBM(1:ND,i)= indice multidimensional correspondiente al indice unidimensional i
!   --------------------------------------------------------------------
!           print *, g%ndim
!           print *, g%igrid
!           print *, g%ii
!           print *, g%n
!           print *, g%ngrid


