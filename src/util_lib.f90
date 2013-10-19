!     ------------------------------------------------------------------
!     [DSYTRFI]   Calcula el determinante y la inversa de una matriz
!                 simétrica a partir de la factorización UDU'.
!     ------------------------------------------------------------------
      SUBROUTINE DSYTRFI(N, A, AI, DET)
      IMPLICIT NONE
      INTEGER N
      REAL*8  A(N,N), AI(N,N), DET
!     Variables locales
      INTEGER IPIV(N), LWORK, INFO, i
      REAL*8 tmp
      REAL*8, ALLOCATABLE ::  WORK(:)
!     ------------------------------------------------------------------
      AI = A
!     Factorizar matriz
!     DSYTRF computes the factorization of a real symmetric matrix A using
!           the Bunch-Kaufman diagonal pivoting method.  The form of the
!           factorization is A = U*D*U**T  or  A = L*D*L**T
      LWORK = -1    ! Determine the block size
      CALL DSYTRF( 'U', N, AI, N, IPIV, tmp, LWORK, INFO )
      LWORK = NINT(tmp)
      ALLOCATE(WORK(LWORK),STAT=i)
      IF (i.NE.0) CALL Error(i)
      CALL DSYTRF( 'U', N, AI, N, IPIV, WORK, LWORK, INFO )
      IF (INFO.NE.0) CALL Error(INFO)
!     Calcular el determinante
!     If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!           interchanged and D(k,k) is a 1-by-1 diagonal block.
!     If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!           columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!           is a 2-by-2 diagonal block.
      tmp = 1.0d0
      DO i = 1, N
         IF (IPIV(i).GT.0) THEN
            tmp = tmp * AI(i,i)
         ELSE IF((i.GT.1).AND.(IPIV(i).LT.0).AND.(IPIV(i).EQ.IPIV(i-1))) THEN
            tmp = tmp * ( AI(i,i) * AI(i-1,i-1) - AI(i-1,i) * AI(i,i-1) )
         END IF
      END DO
      DET = tmp
!     Invertir la matriz
!     DSYTRI computes the inverse of a real symmetric indefinite matrix
!           A using the factorization computed by DSYTRF.
      CALL DSYTRI( 'U', N, AI, N, IPIV, WORK, INFO )
      IF (INFO.NE.0) CALL Error(INFO)
!     Liberar memoria
      DEALLOCATE (WORK,STAT=i)
      IF (i.NE.0) CALL Error(i)
      RETURN
      END SUBROUTINE DSYTRFI


!     ------------------------------------------------------------------
!     [DMXYX] Calcula, dependiendo de TRANSY:
!                         R=XYX' ó XY'X' si TRANSX='N'
!                         R=X'YX ó X'Y'X si TRANSX='T'
!             M y N se establecen de forma que el resultado tiene
!             dimensión MxM
!             Si Y no es necesaria, R puede almacenarse en Y.
!     PENDENTE: PASAR A F90/2003
!     ------------------------------------------------------------------
      SUBROUTINE DMXYX(TRANSX,TRANSY,M,N,X,LDX,Y,LDY,R,LDR)
      IMPLICIT NONE
      INTEGER M,N,LDX,LDY,LDR
      REAL*8  X(LDX,*),Y(LDY,*),R(LDR,*)
      CHARACTER*1  TRANSX,TRANSY
!     Variables locales
      INTEGER i
      CHARACTER*1  TRANSX2
      REAL*8, ALLOCATABLE :: XY(:,:)
!     ------------------------------------------------------------------
!     Asignar memoria
      ALLOCATE (XY(M,N),STAT=i)
      IF (i.NE.0) CALL Error(i)
!     Determinar operación
      IF ((TRANSX.EQ.'t').OR.(TRANSX.EQ.'T')) THEN
          TRANSX2='N'
      ELSE
          TRANSX2='T'
      END IF
      CALL DGEMM (TRANSX,TRANSY,M,N,N,1.0D0,X,LDX,Y,LDY,0.0D0,XY,M) ! XY
      CALL DGEMM ('N',TRANSX2,M,M,N,1.0D0,XY,M,X,LDX,0.0D0,R,LDR)   ! XYX'
!     Liberar memoria
      DEALLOCATE (XY,STAT=i)
      IF (i.NE.0) CALL Error(i)
      RETURN
      END SUBROUTINE DMXYX


!     ------------------------------------------------------------------
!     [DMURRV]    Multiplies a real rectangular matrix by a vector.
!           A — Real NRA by NCA rectangular matrix.   (Input)
!           X — Real vector of length NX.   (Input)
!           Y — Real vector of length NY containing the product
!               A * X if IPATH is equal to 1 and the product
!               trans(A) * X if IPATH is equal to 2.   (Output)
!     PENDENTE: PASAR A F90/2003
!     ------------------------------------------------------------------
      SUBROUTINE DMURRV(NRA,NCA,A,LDA,NX,X,IPATH,NY,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(LDA,*),X(NX),Y(NY)
!     ------------------------------------------------------------------
      IF (IPATH.EQ.1) THEN
         DO 10 I=1,NRA
            SUM=0D0
            DO 20 J=1,NCA
               SUM=SUM+A(I,J)*X(J)
 20         CONTINUE
            Y(I)=SUM
 10      CONTINUE
      ENDIF
      IF (IPATH.EQ.2) THEN
         DO 30 J=1,NCA
            SUM=0D0
            DO 40 I=1,NRA
               SUM=SUM+A(I,J)*X(I)
 40         CONTINUE
            Y(J)=SUM
 30      CONTINUE
      ENDIF
      RETURN
      END
