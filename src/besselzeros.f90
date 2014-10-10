!-----------------------------------------------------------------------
!   [besselzeros.f90]   Utilidades para discretización de la distribución
!                       espectral de un variograma multidimensional
!
!   Interfaces con R:
!       disc_sbv      ptos de discretización para un modelo de variograma de 
!                     Shapiro-Botha (R "disc_sbv")
!
!   Autor: (c) Ruben Fernandez-Casal                Creacion: Abr 2002
!   Revisiones: Mar 2013
!-----------------------------------------------------------------------

!     ------------------------------------------------------------------
!     [disc_sbv]  Obtiene los ptos de discretización de la función de 
!                 distribución espectral para un modelo de Shapiro-Botha 
!                 extendido. Basado en el artículo:
!                 Gorsich y Genton (2001) "On the discretization of 
!                 nonparametric covariogram estimators"
!
!     PARÁMETROS:
!         nx = nº de nodos                                            (I)
!         x(nx) = nodos                                               (O)
!         dim = dimensión correspondiente                             (I)
!         ((dim-2.0)/2.0 orden de la función de Bessel)
!         rango = máximo salto                                        (I)
!     ------------------------------------------------------------------
      SUBROUTINE disc_sbv(nx, x, dim, rango)
      IMPLICIT NONE
!     Parámetros
      INTEGER nx, dim, i
      REAL*8 x(nx), rango, a
!     En el caso infinito se toman equidistantes
      IF (dim.LE.0) THEN
          DO i = 1, nx
              x(i) = i*0.3
          END DO
      ELSE
          a =(dim-2.0)/2.0
          CALL besselzeros(nx, a, x)
      END IF
!     Se reescala por máximo salto
      DO i = 1,nx
          x(i) = x(i)/rango
      END DO
      RETURN
      END SUBROUTINE disc_sbv


!     ------------------------------------------------------------------
!     [Besselzeros]   Calcula los ceros de una función de Bessel de orden
!                     real utilizando el algoritmo propuesto por J.S. Ball
!                     (2000) "Automatic computation of zeros of Bessel
!                     functions and other special functions", 1458-1464,
!                     J. Sci. Comput.
!
!     PARÁMETROS:
!         nt = nº de ceros                                            (I)
!         a = orden de la función de Bessel                           (I)
!         c(nt) = ceros (ordenados)                                   (O)
!     ------------------------------------------------------------------
      SUBROUTINE besselzeros(nt, a, c)
      IMPLICIT NONE
      INTEGER nt
      REAL*8 a, c(nt)
!     Variables locales
      INTEGER i, j, nmax, ierr       
      REAL*8 fl, a1, aa
      REAL*8, ALLOCATABLE :: e(:), d(:), z(:,:)
!     Asignar memoria a variables locales
      nmax = 2*nt
      ALLOCATE (e(nmax+1), d(nmax), z(nmax, nmax), STAT=i)
      IF (i.NE.0) CALL Error(i)
!     Valores iniciales
      z=0.0D0
      do 10 i=1,nmax
          fl=dfloat(i)+.5d0*a+.5d0
          d(i)=.125d0/fl/(fl-1.d0)
          e(i+1)=.125d0/fl/dsqrt(4.d0*fl**2-1.d0)
          z(i,i)=1.0D0
   10 continue
!     input e(i)=alpha, d(i)= beta; tql2 returns eigenvalues
!     of symmetric tridiagonal matrix in d(i); unsorted
      call tql2(nmax,nmax,d,e,z,ierr)
!     form zeros by inverse of g(x)
      do 20 i=1,nmax
          d(i)=1.d0/dsqrt(d(i))
   20 continue
!     sort zeros
      do 30 i=1,nmax
          a1=d(i)
          do 25 j=i+1,nmax
              if(d(j).lt.a1)then
              aa=d(j)
              d(j)=a1
              a1=aa
              endif
   25     continue
          d(i)=a1
   30 continue
      do 40 i=1,nt
          c(i) = d(i)
   40 continue
!     Liberar memoria de variables locales
      DEALLOCATE (e,d,z,STAT=i)
      IF (i.NE.0) CALL Error(i)
      RETURN
      END SUBROUTINE Besselzeros


!     ------------------------------------------------------------------
!     [tql2]  Subrutina del paquete EISPACK (descargada de internet)
!             para el cálculo de autovalores (y autovectores) de una 
!             matriz tridiagonal simétrica.
!
!     NOTA: Se eliminó la ordenación de autovalores y autovectores.
!     ------------------------------------------------------------------
      subroutine tql2(nm,n,d,e,z,ierr)
      integer n,nm,ierr
      double precision d(n),e(n),z(nm,n)
!
!     this subroutine is a translation of the algol procedure tql2,
!     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
!     wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a symmetric tridiagonal matrix by the ql method.
!     the eigenvectors of a full symmetric matrix can also
!     be found if  tred2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        z contains the transformation matrix produced in the
!          reduction by  tred2, if performed.  if the eigenvectors
!          of the tridiagonal matrix are desired, z must contain
!          the identity matrix.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1,2,...,ierr-1.
!
!        e has been destroyed.
!
!        z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  if an error exit is made,
!          z contains the eigenvectors associated with the stored
!          eigenvalues.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     calls pythag for  dsqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
      integer i,j,k,l,m,ii,l1,l2,mml
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
!
      ierr = 0
      if (n .eq. 1) go to 1001
!
      do 100 i = 2, n
  100 e(i-1) = e(i)
!
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
!
      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
!     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
!     .......... e(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
  110    continue
!
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
!     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
!
         do 140 i = l2, n
  140    d(i) = d(i) - h
!
  145    f = f + h
!     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
!     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
!
  200    continue
!
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
!     .......... order eigenvalues and eigenvectors ..........
!      do 300 ii = 2, n
!         i = ii - 1
!         k = i
!         p = d(i)
!
!         do 260 j = ii, n
!            if (d(j) .ge. p) go to 260
!            k = j
!            p = d(j)
!  260    continue
!
!         if (k .eq. i) go to 300
!         d(k) = d(i)
!         d(i) = p
!
!         do 280 j = 1, n
!            p = z(j,i)
!            z(j,i) = z(j,k)
!            z(j,k) = p
!  280    continue
!
!  300 continue
!
      go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end



!     ------------------------------------------------------------------
!     [pythag]    Función utilizada por la subrutina tql2 
!                 (paquete EISPACK).
!     ------------------------------------------------------------------
!     finds dsqrt(a**2+b**2) without overflow or destructive underflow
!     ------------------------------------------------------------------
      double precision function pythag(a,b)
      double precision a,b
!
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
