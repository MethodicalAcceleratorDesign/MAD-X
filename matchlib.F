      SUBROUTINE DLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B, LDB, BX, LDBX,  &
     &PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM,                     &
     &POLES, DIFL, DIFR, Z, K, C, S, WORK, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            GIVPTR, ICOMPQ, INFO, K, LDB, LDBX, LDGCOL,    &
     &LDGNUM, NL, NR, NRHS, SQRE
      DOUBLE PRECISION   C, S
!     ..
!     .. Array Arguments ..
      INTEGER            GIVCOL( LDGCOL, * ), PERM( * )
      DOUBLE PRECISION   B( LDB, * ), BX( LDBX, * ), DIFL( * ),         &
     &DIFR( LDGNUM, * ), GIVNUM( LDGNUM, * ),                           &
     &POLES( LDGNUM, * ), WORK( * ), Z( * )
!     ..
!
!  Purpose
!  =======
!
!  DLALS0 applies back the multiplying factors of either the left or the
!  right singular vector matrix of a diagonal matrix appended by a row
!  to the right hand side matrix B in solving the least squares problem
!  using the divide-and-conquer SVD approach.
!
!  For the left singular vector matrix, three types of orthogonal
!  matrices are involved:
!
!  (1L) Givens rotations: the number of such rotations is GIVPTR; the
!       pairs of columns/rows they were applied to are stored in GIVCOL;
!       and the C- and S-values of these rotations are stored in GIVNUM.
!
!  (2L) Permutation. The (NL+1)-st row of B is to be moved to the first
!       row, and for J=2:N, PERM(J)-th row of B is to be moved to the
!       J-th row.
!
!  (3L) The left singular vector matrix of the remaining matrix.
!
!  For the right singular vector matrix, four types of orthogonal
!  matrices are involved:
!
!  (1R) The right singular vector matrix of the remaining matrix.
!
!  (2R) If SQRE = 1, one extra Givens rotation to generate the right
!       null space.
!
!  (3R) The inverse transformation of (2L).
!
!  (4R) The inverse transformation of (1L).
!
!  Arguments
!  =========
!
!  ICOMPQ (input) INTEGER
!         Specifies whether singular vectors are to be computed in
!         factored form:
!         = 0: Left singular vector matrix.
!         = 1: Right singular vector matrix.
!
!  NL     (input) INTEGER
!         The row dimension of the upper block. NL >= 1.
!
!  NR     (input) INTEGER
!         The row dimension of the lower block. NR >= 1.
!
!  SQRE   (input) INTEGER
!         = 0: the lower block is an NR-by-NR square matrix.
!         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.
!
!         The bidiagonal matrix has row dimension N = NL + NR + 1,
!         and column dimension M = N + SQRE.
!
!  NRHS   (input) INTEGER
!         The number of columns of B and BX. NRHS must be at least 1.
!
!  B      (input/output) DOUBLE PRECISION array, dimension ( LDB, NRHS )
!         On input, B contains the right hand sides of the least
!         squares problem in rows 1 through M. On output, B contains
!         the solution X in rows 1 through N.
!
!  LDB    (input) INTEGER
!         The leading dimension of B. LDB must be at least
!         max(1,MAX( M, N ) ).
!
!  BX     (workspace) DOUBLE PRECISION array, dimension ( LDBX, NRHS )
!
!  LDBX   (input) INTEGER
!         The leading dimension of BX.
!
!  PERM   (input) INTEGER array, dimension ( N )
!         The permutations (from deflation and sorting) applied
!         to the two blocks.
!
!  GIVPTR (input) INTEGER
!         The number of Givens rotations which took place in this
!         subproblem.
!
!  GIVCOL (input) INTEGER array, dimension ( LDGCOL, 2 )
!         Each pair of numbers indicates a pair of rows/columns
!         involved in a Givens rotation.
!
!  LDGCOL (input) INTEGER
!         The leading dimension of GIVCOL, must be at least N.
!
!  GIVNUM (input) DOUBLE PRECISION array, dimension ( LDGNUM, 2 )
!         Each number indicates the C or S value used in the
!         corresponding Givens rotation.
!
!  LDGNUM (input) INTEGER
!         The leading dimension of arrays DIFR, POLES and
!         GIVNUM, must be at least K.
!
!  POLES  (input) DOUBLE PRECISION array, dimension ( LDGNUM, 2 )
!         On entry, POLES(1:K, 1) contains the new singular
!         values obtained from solving the secular equation, and
!         POLES(1:K, 2) is an array containing the poles in the secular
!         equation.
!
!  DIFL   (input) DOUBLE PRECISION array, dimension ( K ).
!         On entry, DIFL(I) is the distance between I-th updated
!         (undeflated) singular value and the I-th (undeflated) old
!         singular value.
!
!  DIFR   (input) DOUBLE PRECISION array, dimension ( LDGNUM, 2 ).
!         On entry, DIFR(I, 1) contains the distances between I-th
!         updated (undeflated) singular value and the I+1-th
!         (undeflated) old singular value. And DIFR(I, 2) is the
!         normalizing factor for the I-th right singular vector.
!
!  Z      (input) DOUBLE PRECISION array, dimension ( K )
!         Contain the components of the deflation-adjusted updating row
!         vector.
!
!  K      (input) INTEGER
!         Contains the dimension of the non-deflated matrix,
!         This is the order of the related secular equation. 1 <= K <=N.
!
!  C      (input) DOUBLE PRECISION
!         C contains garbage if SQRE =0 and the C-value of a Givens
!         rotation related to the right null space if SQRE = 1.
!
!  S      (input) DOUBLE PRECISION
!         S contains garbage if SQRE =0 and the S-value of a Givens
!         rotation related to the right null space if SQRE = 1.
!
!  WORK   (workspace) DOUBLE PRECISION array, dimension ( K )
!
!  INFO   (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ming Gu and Ren-Cang Li, Computer Science Division, University of
!       California at Berkeley, USA
!     Osni Marques, LBNL/NERSC, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO, NEGONE
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0, NEGONE = -1.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, M, N, NLP1
      DOUBLE PRECISION   DIFLJ, DIFRJ, DJ, DSIGJ, DSIGJP, TEMP
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPL, DGEMV, DLACPY, DLASCL, DROT, DSCAL,     &
     &XERBLA
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3, DNRM2
      EXTERNAL           DLAMC3, DNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) THEN
        INFO = -1
      ELSE IF( NL.LT.1 ) THEN
        INFO = -2
      ELSE IF( NR.LT.1 ) THEN
        INFO = -3
      ELSE IF( ( SQRE.LT.0 ) .OR. ( SQRE.GT.1 ) ) THEN
        INFO = -4
      END IF
!
      N = NL + NR + 1
!
      IF( NRHS.LT.1 ) THEN
        INFO = -5
      ELSE IF( LDB.LT.N ) THEN
        INFO = -7
      ELSE IF( LDBX.LT.N ) THEN
        INFO = -9
      ELSE IF( GIVPTR.LT.0 ) THEN
        INFO = -11
      ELSE IF( LDGCOL.LT.N ) THEN
        INFO = -13
      ELSE IF( LDGNUM.LT.N ) THEN
        INFO = -15
      ELSE IF( K.LT.1 ) THEN
        INFO = -20
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DLALS0', -INFO )
        RETURN
      END IF
!
      M = N + SQRE
      NLP1 = NL + 1
!
      IF( ICOMPQ.EQ.0 ) THEN
!
!        Apply back orthogonal transformations from the left.
!
!        Step (1L): apply back the Givens rotations performed.
!
        DO I = 1, GIVPTR
          CALL DROT( NRHS, B( GIVCOL( I, 2 ), 1 ), LDB,                 &
     &B( GIVCOL( I, 1 ), 1 ), LDB, GIVNUM( I, 2 ),                      &
     &GIVNUM( I, 1 ) )
        enddo
!
!        Step (2L): permute rows of B.
!
        CALL DCOPL( NRHS, B( NLP1, 1 ), LDB, BX( 1, 1 ), LDBX )
        DO I = 2, N
          CALL DCOPL( NRHS, B( PERM( I ), 1 ), LDB, BX( I, 1 ), LDBX )
        enddo
!
!        Step (3L): apply the inverse of the left singular vector
!        matrix to BX.
!
        IF( K.EQ.1 ) THEN
          CALL DCOPL( NRHS, BX, LDBX, B, LDB )
          IF( Z( 1 ).LT.ZERO ) THEN
            CALL DSCAL( NRHS, NEGONE, B, LDB )
          END IF
        ELSE
          DO J = 1, K
            DIFLJ = DIFL( J )
            DJ = POLES( J, 1 )
            DSIGJ = -POLES( J, 2 )
            IF( J.LT.K ) THEN
              DIFRJ = -DIFR( J, 1 )
              DSIGJP = -POLES( J+1, 2 )
            END IF
            IF( ( Z( J ).EQ.ZERO ) .OR. ( POLES( J, 2 ).EQ.ZERO ) )     &
     &THEN
              WORK( J ) = ZERO
            ELSE
              WORK( J ) = -POLES( J, 2 )*Z( J ) / DIFLJ /               &
     &( POLES( J, 2 )+DJ )
            END IF
            DO I = 1, J - 1
              IF( ( Z( I ).EQ.ZERO ) .OR.                               &
     &( POLES( I, 2 ).EQ.ZERO ) ) THEN
                WORK( I ) = ZERO
              ELSE
                WORK( I ) = POLES( I, 2 )*Z( I ) /                      &
     &( DLAMC3( POLES( I, 2 ), DSIGJ )-                                 &
     &DIFLJ ) / ( POLES( I, 2 )+DJ )
              END IF
            enddo
            DO I = J + 1, K
              IF( ( Z( I ).EQ.ZERO ) .OR.                               &
     &( POLES( I, 2 ).EQ.ZERO ) ) THEN
                WORK( I ) = ZERO
              ELSE
                WORK( I ) = POLES( I, 2 )*Z( I ) /                      &
     &( DLAMC3( POLES( I, 2 ), DSIGJP )+                                &
     &DIFRJ ) / ( POLES( I, 2 )+DJ )
              END IF
            enddo
            WORK( 1 ) = NEGONE
            TEMP = DNRM2( K, WORK, 1 )
            CALL DGEMV( 'T', K, NRHS, ONE, BX, LDBX, WORK, 1, ZERO,     &
     &B( J, 1 ), LDB )
            CALL DLASCL( 'G', 0, 0, TEMP, ONE, 1, NRHS, B( J, 1 ),      &
     &LDB, INFO )
          enddo
        END IF
!
!        Move the deflated rows of BX to B also.
!
        IF( K.LT.MAX( M, N ) )                                          &
     &CALL DLACPY( 'A', N-K, NRHS, BX( K+1, 1 ), LDBX,                  &
     &B( K+1, 1 ), LDB )
      ELSE
!
!        Apply back the right orthogonal transformations.
!
!        Step (1R): apply back the new right singular vector matrix
!        to B.
!
        IF( K.EQ.1 ) THEN
          CALL DCOPL( NRHS, B, LDB, BX, LDBX )
        ELSE
          DO J = 1, K
            DSIGJ = POLES( J, 2 )
            IF( Z( J ).EQ.ZERO ) THEN
              WORK( J ) = ZERO
            ELSE
              WORK( J ) = -Z( J ) / DIFL( J ) /                         &
     &( DSIGJ+POLES( J, 1 ) ) / DIFR( J, 2 )
            END IF
            DO I = 1, J - 1
              IF( Z( J ).EQ.ZERO ) THEN
                WORK( I ) = ZERO
              ELSE
                WORK( I ) = Z( J ) / ( DLAMC3( DSIGJ, -POLES( I+1,      &
     &2 ) )-DIFR( I, 1 ) ) /                                            &
     &( DSIGJ+POLES( I, 1 ) ) / DIFR( I, 2 )
              END IF
            enddo
            DO I = J + 1, K
              IF( Z( J ).EQ.ZERO ) THEN
                WORK( I ) = ZERO
              ELSE
                WORK( I ) = Z( J ) / ( DLAMC3( DSIGJ, -POLES( I,        &
     &2 ) )-DIFL( I ) ) /                                               &
     &( DSIGJ+POLES( I, 1 ) ) / DIFR( I, 2 )
              END IF
            enddo
            CALL DGEMV( 'T', K, NRHS, ONE, B, LDB, WORK, 1, ZERO,       &
     &BX( J, 1 ), LDBX )
          enddo
        END IF
!
!        Step (2R): if SQRE = 1, apply back the rotation that is
!        related to the right null space of the subproblem.
!
        IF( SQRE.EQ.1 ) THEN
          CALL DCOPL( NRHS, B( M, 1 ), LDB, BX( M, 1 ), LDBX )
          CALL DROT( NRHS, BX( 1, 1 ), LDBX, BX( M, 1 ), LDBX, C, S )
        END IF
        IF( K.LT.MAX( M, N ) )                                          &
     &CALL DLACPY( 'A', N-K, NRHS, B( K+1, 1 ), LDB, BX( K+1, 1 ),      &
     &LDBX )
!
!        Step (3R): permute rows of B.
!
        CALL DCOPL( NRHS, BX( 1, 1 ), LDBX, B( NLP1, 1 ), LDB )
        IF( SQRE.EQ.1 ) THEN
          CALL DCOPL( NRHS, BX( M, 1 ), LDBX, B( M, 1 ), LDB )
        END IF
        DO I = 2, N
          CALL DCOPL( NRHS, BX( I, 1 ), LDBX, B( PERM( I ), 1 ), LDB )
        enddo
!
!        Step (4R): apply back the Givens rotations performed.
!
        DO I = GIVPTR, 1, -1
          CALL DROT( NRHS, B( GIVCOL( I, 2 ), 1 ), LDB,                 &
     &B( GIVCOL( I, 1 ), 1 ), LDB, GIVNUM( I, 2 ),                      &
     &-GIVNUM( I, 1 ) )
        enddo
      END IF
!
      RETURN
!
!     End of DLALS0
!
      END
      SUBROUTINE DLALSA( ICOMPQ, SMLSIZ, N, NRHS, B, LDB, BX, LDBX, U,  &
     &LDU, VT, K, DIFL, DIFR, Z, POLES, GIVPTR,                         &
     &GIVCOL, LDGCOL, PERM, GIVNUM, C, S, WORK,                         &
     &IWORK, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            ICOMPQ, INFO, LDB, LDBX, LDGCOL, LDU, N, NRHS, &
     &SMLSIZ
!     ..
!     .. Array Arguments ..
      INTEGER            GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ),  &
     &K( * ), PERM( LDGCOL, * )
      DOUBLE PRECISION   B( LDB, * ), BX( LDBX, * ), C( * ),            &
     &DIFL( LDU, * ), DIFR( LDU, * ),                                   &
     &GIVNUM( LDU, * ), POLES( LDU, * ), S( * ),                        &
     &U( LDU, * ), VT( LDU, * ), WORK( * ),                             &
     &Z( LDU, * )
!     ..
!
!  Purpose
!  =======
!
!  DLALSA is an itermediate step in solving the least squares problem
!  by computing the SVD of the coefficient matrix in compact form (The
!  singular vectors are computed as products of simple orthorgonal
!  matrices.).
!
!  If ICOMPQ = 0, DLALSA applies the inverse of the left singular vector
!  matrix of an upper bidiagonal matrix to the right hand side; and if
!  ICOMPQ = 1, DLALSA applies the right singular vector matrix to the
!  right hand side. The singular vector matrices were generated in
!  compact form by DLALSA.
!
!  Arguments
!  =========
!
!
!  ICOMPQ (input) INTEGER
!         Specifies whether the left or the right singular vector
!         matrix is involved.
!         = 0: Left singular vector matrix
!         = 1: Right singular vector matrix
!
!  SMLSIZ (input) INTEGER
!         The maximum size of the subproblems at the bottom of the
!         computation tree.
!
!  N      (input) INTEGER
!         The row and column dimensions of the upper bidiagonal matrix.
!
!  NRHS   (input) INTEGER
!         The number of columns of B and BX. NRHS must be at least 1.
!
!  B      (input/output) DOUBLE PRECISION array, dimension ( LDB, NRHS )
!         On input, B contains the right hand sides of the least
!         squares problem in rows 1 through M.
!         On output, B contains the solution X in rows 1 through N.
!
!  LDB    (input) INTEGER
!         The leading dimension of B in the calling subprogram.
!         LDB must be at least max(1,MAX( M, N ) ).
!
!  BX     (output) DOUBLE PRECISION array, dimension ( LDBX, NRHS )
!         On exit, the result of applying the left or right singular
!         vector matrix to B.
!
!  LDBX   (input) INTEGER
!         The leading dimension of BX.
!
!  U      (input) DOUBLE PRECISION array, dimension ( LDU, SMLSIZ ).
!         On entry, U contains the left singular vector matrices of all
!         subproblems at the bottom level.
!
!  LDU    (input) INTEGER, LDU = > N.
!         The leading dimension of arrays U, VT, DIFL, DIFR,
!         POLES, GIVNUM, and Z.
!
!  VT     (input) DOUBLE PRECISION array, dimension ( LDU, SMLSIZ+1 ).
!         On entry, VT' contains the right singular vector matrices of
!         all subproblems at the bottom level.
!
!  K      (input) INTEGER array, dimension ( N ).
!
!  DIFL   (input) DOUBLE PRECISION array, dimension ( LDU, NLVL ).
!         where NLVL = INT(log_2 (N/(SMLSIZ+1))) + 1.
!
!  DIFR   (input) DOUBLE PRECISION array, dimension ( LDU, 2 * NLVL ).
!         On entry, DIFL(*, I) and DIFR(*, 2 * I -1) record
!         distances between singular values on the I-th level and
!         singular values on the (I -1)-th level, and DIFR(*, 2 * I)
!         record the normalizing factors of the right singular vectors
!         matrices of subproblems on I-th level.
!
!  Z      (input) DOUBLE PRECISION array, dimension ( LDU, NLVL ).
!         On entry, Z(1, I) contains the components of the deflation-
!         adjusted updating row vector for subproblems on the I-th
!         level.
!
!  POLES  (input) DOUBLE PRECISION array, dimension ( LDU, 2 * NLVL ).
!         On entry, POLES(*, 2 * I -1: 2 * I) contains the new and old
!         singular values involved in the secular equations on the I-th
!         level.
!
!  GIVPTR (input) INTEGER array, dimension ( N ).
!         On entry, GIVPTR( I ) records the number of Givens
!         rotations performed on the I-th problem on the computation
!         tree.
!
!  GIVCOL (input) INTEGER array, dimension ( LDGCOL, 2 * NLVL ).
!         On entry, for each I, GIVCOL(*, 2 * I - 1: 2 * I) records the
!         locations of Givens rotations performed on the I-th level on
!         the computation tree.
!
!  LDGCOL (input) INTEGER, LDGCOL = > N.
!         The leading dimension of arrays GIVCOL and PERM.
!
!  PERM   (input) INTEGER array, dimension ( LDGCOL, NLVL ).
!         On entry, PERM(*, I) records permutations done on the I-th
!         level of the computation tree.
!
!  GIVNUM (input) DOUBLE PRECISION array, dimension ( LDU, 2 * NLVL ).
!         On entry, GIVNUM(*, 2 *I -1 : 2 * I) records the C- and S-
!         values of Givens rotations performed on the I-th level on the
!         computation tree.
!
!  C      (input) DOUBLE PRECISION array, dimension ( N ).
!         On entry, if the I-th subproblem is not square,
!         C( I ) contains the C-value of a Givens rotation related to
!         the right null space of the I-th subproblem.
!
!  S      (input) DOUBLE PRECISION array, dimension ( N ).
!         On entry, if the I-th subproblem is not square,
!         S( I ) contains the S-value of a Givens rotation related to
!         the right null space of the I-th subproblem.
!
!  WORK   (workspace) DOUBLE PRECISION array.
!         The dimension must be at least N.
!
!  IWORK  (workspace) INTEGER array.
!         The dimension must be at least 3 * N
!
!  INFO   (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ming Gu and Ren-Cang Li, Computer Science Division, University of
!       California at Berkeley, USA
!     Osni Marques, LBNL/NERSC, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I1, IC, IM1, INODE, J, LF, LL, LVL, LVL2,   &
     &ND, NDB1, NDIML, NDIMR, NL, NLF, NLP1, NLVL,                      &
     &NR, NRF, NRP1, SQRE
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPL, DGEMM, DLALS0, DLASDT, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) THEN
        INFO = -1
      ELSE IF( SMLSIZ.LT.3 ) THEN
        INFO = -2
      ELSE IF( N.LT.SMLSIZ ) THEN
        INFO = -3
      ELSE IF( NRHS.LT.1 ) THEN
        INFO = -4
      ELSE IF( LDB.LT.N ) THEN
        INFO = -6
      ELSE IF( LDBX.LT.N ) THEN
        INFO = -8
      ELSE IF( LDU.LT.N ) THEN
        INFO = -10
      ELSE IF( LDGCOL.LT.N ) THEN
        INFO = -19
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DLALSA', -INFO )
        RETURN
      END IF
!
!     Book-keeping and  setting up the computation tree.
!
      INODE = 1
      NDIML = INODE + N
      NDIMR = NDIML + N
!
      CALL DLASDT( N, NLVL, ND, IWORK( INODE ), IWORK( NDIML ),         &
     &IWORK( NDIMR ), SMLSIZ )
!
!     The following code applies back the left singular vector factors.
!     For applying back the right singular vector factors, go to 50.
!
      IF( ICOMPQ.EQ.1 ) THEN
        GO TO 50
      END IF
!
!     The nodes on the bottom level of the tree were solved
!     by DLASDQ. The corresponding left and right singular vector
!     matrices are in explicit form. First apply back the left
!     singular vector matrices.
!
      NDB1 = ( ND+1 ) / 2
      DO I = NDB1, ND
!
!        IC : center row of each node
!        NL : number of rows of left  subproblem
!        NR : number of rows of right subproblem
!        NLF: starting row of the left   subproblem
!        NRF: starting row of the right  subproblem
!
        I1 = I - 1
        IC = IWORK( INODE+I1 )
        NL = IWORK( NDIML+I1 )
        NR = IWORK( NDIMR+I1 )
        NLF = IC - NL
        NRF = IC + 1
        CALL DGEMM( 'T', 'N', NL, NRHS, NL, ONE, U( NLF, 1 ), LDU,      &
     &B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX )
        CALL DGEMM( 'T', 'N', NR, NRHS, NR, ONE, U( NRF, 1 ), LDU,      &
     &B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX )
      enddo
!
!     Next copy the rows of B that correspond to unchanged rows
!     in the bidiagonal matrix to BX.
!
      DO I = 1, ND
        IC = IWORK( INODE+I-1 )
        CALL DCOPL( NRHS, B( IC, 1 ), LDB, BX( IC, 1 ), LDBX )
      enddo
!
!     Finally go through the left singular vector matrices of all
!     the other subproblems bottom-up on the tree.
!
      J = 2**NLVL
      SQRE = 0
!
      DO LVL = NLVL, 1, -1
        LVL2 = 2*LVL - 1
!
!        find the first node LF and last node LL on
!        the current level LVL
!
        IF( LVL.EQ.1 ) THEN
          LF = 1
          LL = 1
        ELSE
          LF = 2**( LVL-1 )
          LL = 2*LF - 1
        END IF
        DO I = LF, LL
          IM1 = I - 1
          IC = IWORK( INODE+IM1 )
          NL = IWORK( NDIML+IM1 )
          NR = IWORK( NDIMR+IM1 )
          NLF = IC - NL
          NRF = IC + 1
          J = J - 1
          CALL DLALS0( ICOMPQ, NL, NR, SQRE, NRHS, BX( NLF, 1 ), LDBX,  &
     &B( NLF, 1 ), LDB, PERM( NLF, LVL ),                               &
     &GIVPTR( J ), GIVCOL( NLF, LVL2 ), LDGCOL,                         &
     &GIVNUM( NLF, LVL2 ), LDU, POLES( NLF, LVL2 ),                     &
     &DIFL( NLF, LVL ), DIFR( NLF, LVL2 ),                              &
     &Z( NLF, LVL ), K( J ), C( J ), S( J ), WORK,                      &
     &INFO )
        enddo
      enddo
      GO TO 90
!
!     ICOMPQ = 1: applying back the right singular vector factors.
!
   50 CONTINUE
!
!     First now go through the right singular vector matrices of all
!     the tree nodes top-down.
!
      J = 0
      DO LVL = 1, NLVL
        LVL2 = 2*LVL - 1
!
!        Find the first node LF and last node LL on
!        the current level LVL.
!
        IF( LVL.EQ.1 ) THEN
          LF = 1
          LL = 1
        ELSE
          LF = 2**( LVL-1 )
          LL = 2*LF - 1
        END IF
        DO I = LL, LF, -1
          IM1 = I - 1
          IC = IWORK( INODE+IM1 )
          NL = IWORK( NDIML+IM1 )
          NR = IWORK( NDIMR+IM1 )
          NLF = IC - NL
          NRF = IC + 1
          IF( I.EQ.LL ) THEN
            SQRE = 0
          ELSE
            SQRE = 1
          END IF
          J = J + 1
          CALL DLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B( NLF, 1 ), LDB,    &
     &BX( NLF, 1 ), LDBX, PERM( NLF, LVL ),                             &
     &GIVPTR( J ), GIVCOL( NLF, LVL2 ), LDGCOL,                         &
     &GIVNUM( NLF, LVL2 ), LDU, POLES( NLF, LVL2 ),                     &
     &DIFL( NLF, LVL ), DIFR( NLF, LVL2 ),                              &
     &Z( NLF, LVL ), K( J ), C( J ), S( J ), WORK,                      &
     &INFO )
        enddo
      enddo
!
!     The nodes on the bottom level of the tree were solved
!     by DLASDQ. The corresponding right singular vector
!     matrices are in explicit form. Apply them back.
!
      NDB1 = ( ND+1 ) / 2
      DO I = NDB1, ND
        I1 = I - 1
        IC = IWORK( INODE+I1 )
        NL = IWORK( NDIML+I1 )
        NR = IWORK( NDIMR+I1 )
        NLP1 = NL + 1
        IF( I.EQ.ND ) THEN
          NRP1 = NR
        ELSE
          NRP1 = NR + 1
        END IF
        NLF = IC - NL
        NRF = IC + 1
        CALL DGEMM( 'T', 'N', NLP1, NRHS, NLP1, ONE, VT( NLF, 1 ), LDU, &
     &B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX )
        CALL DGEMM( 'T', 'N', NRP1, NRHS, NRP1, ONE, VT( NRF, 1 ), LDU, &
     &B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX )
      enddo
!
   90 CONTINUE
!
      RETURN
!
!     End of DLALSA
!
      END
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
!     .. Scalar Arguments ..
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
!     ..
!
!  Purpose
!  =======
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      DOUBLE PRECISION DMAX
      INTEGER I,IX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DABS
!     ..
      IDAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IDAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) GO TO 20
!
!        code for increment not equal to 1
!
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO I = 2,N
        IF (DABS(DX(IX)).LE.DMAX) GO TO 5
        IDAMAX = I
        DMAX = DABS(DX(IX))
    5   IX = IX + INCX
      enddo
      RETURN
!
!        code for increment equal to 1
!
 20   continue
      DMAX = DABS(DX(1))
      DO I = 2,N
        IF (DABS(DX(I)).LE.DMAX) GO TO 30
        IDAMAX = I
        DMAX = DABS(DX(I))
 30     CONTINUE
      ENDDO
      RETURN
      END
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
      implicit none
!     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  Purpose
!  =======
!
!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      DDOT = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
      enddo
      DDOT = DTEMP
      RETURN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
      enddo
      IF (N.LT.5) GO TO 60
   40 MP1 = M + 1
      DO I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) +                 &
     &DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
      enddo
   60 DDOT = DTEMP
      RETURN
      END
      SUBROUTINE DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,    &
     &WORK, LWORK, IWORK, INFO )
      implicit none
!
!  -- LAPACK driver routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
      DOUBLE PRECISION   RCOND
!     ..
!     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGELSD computes the minimum-norm solution to a real linear least
!  squares problem:
!      minimize 2-norm(| b - A*x |)
!  using the singular value decomposition (SVD) of A. A is an M-by-N
!  matrix which may be rank-deficient.
!
!  Several right hand side vectors b and solution vectors x can be
!  handled in a single call; they are stored as the columns of the
!  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
!  matrix X.
!
!  The problem is solved in three steps:
!  (1) Reduce the coefficient matrix A to bidiagonal form with
!      Householder transformations, reducing the original problem
!      into a "bidiagonal least squares problem" (BLS)
!  (2) Solve the BLS using a divide and conquer approach.
!  (3) Apply back all the Householder tranformations to solve
!      the original least squares problem.
!
!  The effective rank of A is determined by treating as zero those
!  singular values which are less than RCOND times the largest singular
!  value.
!
!  The divide and conquer algorithm makes very mild assumptions about
!  floating point arithmetic. It will work on machines with a guard
!  digit in add/subtract, or on those binary machines without guard
!  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
!  Cray-2. It could conceivably fail on hexadecimal or decimal machines
!  without guard digits, but we know of none.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of A. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of A. N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrices B and X. NRHS >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit, A has been destroyed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the M-by-NRHS right hand side matrix B.
!          On exit, B is overwritten by the N-by-NRHS solution
!          matrix X.  If m >= n and RANK = n, the residual
!          sum-of-squares for the solution in the i-th column is given
!          by the sum of squares of elements n+1:m in that column.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B. LDB >= max(1,max(M,N)).
!
!  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
!          The singular values of A in decreasing order.
!          The condition number of A in the 2-norm = S(1)/S(min(m,n)).
!
!  RCOND   (input) DOUBLE PRECISION
!          RCOND is used to determine the effective rank of A.
!          Singular values S(i) <= RCOND*S(1) are treated as zero.
!          If RCOND < 0, machine precision is used instead.
!
!  RANK    (output) INTEGER
!          The effective rank of A, i.e., the number of singular values
!          which are greater than RCOND*S(1).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK must be at least 1.
!          The exact minimum amount of workspace needed depends on M,
!          N and NRHS. As long as LWORK is at least
!              12*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2,
!          if M is greater than or equal to N or
!              12*M + 2*M*SMLSIZ + 8*M*NLVL + M*NRHS + (SMLSIZ+1)**2,
!          if M is less than N, the code will execute correctly.
!          SMLSIZ is returned by ILAENV and is equal to the maximum
!          size of the subproblems at the bottom of the computation
!          tree (usually about 25), and
!             NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 )
!          For good performance, LWORK should generally be larger.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  IWORK   (workspace) INTEGER array, dimension (MAX(1,LIWORK))
!          LIWORK >= 3 * MINMN * NLVL + 11 * MINMN,
!          where MINMN = MIN( M,N ).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  the algorithm for computing the SVD failed to converge;
!                if INFO = i, i off-diagonal elements of an intermediate
!                bidiagonal form did not converge to zero.
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ming Gu and Ren-Cang Li, Computer Science Division, University of
!       California at Berkeley, USA
!     Osni Marques, LBNL/NERSC, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            IASCL, IBSCL, IE, IL, ITAU, ITAUP, ITAUQ,      &
     &LDWORK, MAXMN, MAXWRK, MINMN, MINWRK, MM,                         &
     &MNTHR, NLVL, NWORK, SMLSIZ, WLALSD
      DOUBLE PRECISION   ANRM, BIGNUM, BNRM, EPS, SFMIN, SMLNUM
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEBRD, DGELQF, DGEQRF, DLABAD, DLACPY, DLALSD,&
     &DLASCL, DLASET, DORMBR, DORMLQ, DORMQR, XERBLA
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           ILAENV, DLAMCH, DLANGE
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, LOG, MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
!
      INFO = 0
      MINMN = MIN( M, N )
      MAXMN = MAX( M, N )
      MNTHR = ILAENV( 6, 'DGELSD', ' ', M, N, NRHS, -1 )
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
        INFO = -1
      ELSE IF( N.LT.0 ) THEN
        INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
        INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -5
      ELSE IF( LDB.LT.MAX( 1, MAXMN ) ) THEN
        INFO = -7
      END IF
!
      SMLSIZ = ILAENV( 9, 'DGELSD', ' ', 0, 0, 0, 0 )
!
!     Compute workspace.
!     (Note: Comments in the code beginning "Workspace:" describe the
!     minimal amount of workspace needed at that point in the code,
!     as well as the preferred amount for good performance.
!     NB refers to the optimal block size for the immediately
!     following subroutine, as returned by ILAENV.)
!
      MINWRK = 1
      MINMN = MAX( 1, MINMN )
      NLVL = MAX( INT( LOG( DBLE( MINMN ) / DBLE( SMLSIZ+1 ) ) /        &
     &LOG( TWO ) ) + 1, 0 )
!
      IF( INFO.EQ.0 ) THEN
        MAXWRK = 0
        MM = M
        IF( M.GE.N .AND. M.GE.MNTHR ) THEN
!
!           Path 1a - overdetermined, with many more rows than columns.
!
          MM = N
          MAXWRK = MAX( MAXWRK, N+N*ILAENV( 1, 'DGEQRF', ' ', M, N,     &
     &-1, -1 ) )
          MAXWRK = MAX( MAXWRK, N+NRHS*                                 &
     &ILAENV( 1, 'DORMQR', 'LT', M, NRHS, N, -1 ) )
        END IF
        IF( M.GE.N ) THEN
!
!           Path 1 - overdetermined or exactly determined.
!
          MAXWRK = MAX( MAXWRK, 3*N+( MM+N )*                           &
     &ILAENV( 1, 'DGEBRD', ' ', MM, N, -1, -1 ) )
          MAXWRK = MAX( MAXWRK, 3*N+NRHS*                               &
     &ILAENV( 1, 'DORMBR', 'QLT', MM, NRHS, N, -1 ) )
          MAXWRK = MAX( MAXWRK, 3*N+( N-1 )*                            &
     &ILAENV( 1, 'DORMBR', 'PLN', N, NRHS, N, -1 ) )
          WLALSD = 9*N+2*N*SMLSIZ+8*N*NLVL+N*NRHS+(SMLSIZ+1)**2
          MAXWRK = MAX( MAXWRK, 3*N+WLALSD )
          MINWRK = MAX( 3*N+MM, 3*N+NRHS, 3*N+WLALSD )
        END IF
        IF( N.GT.M ) THEN
          WLALSD = 9*M+2*M*SMLSIZ+8*M*NLVL+M*NRHS+(SMLSIZ+1)**2
          IF( N.GE.MNTHR ) THEN
!
!              Path 2a - underdetermined, with many more columns
!              than rows.
!
            MAXWRK = M + M*ILAENV( 1, 'DGELQF', ' ', M, N, -1, -1 )
            MAXWRK = MAX( MAXWRK, M*M+4*M+2*M*                          &
     &ILAENV( 1, 'DGEBRD', ' ', M, M, -1, -1 ) )
            MAXWRK = MAX( MAXWRK, M*M+4*M+NRHS*                         &
     &ILAENV( 1, 'DORMBR', 'QLT', M, NRHS, M, -1 ) )
            MAXWRK = MAX( MAXWRK, M*M+4*M+( M-1 )*                      &
     &ILAENV( 1, 'DORMBR', 'PLN', M, NRHS, M, -1 ) )
            IF( NRHS.GT.1 ) THEN
              MAXWRK = MAX( MAXWRK, M*M+M+M*NRHS )
            ELSE
              MAXWRK = MAX( MAXWRK, M*M+2*M )
            END IF
            MAXWRK = MAX( MAXWRK, M+NRHS*                               &
     &ILAENV( 1, 'DORMLQ', 'LT', N, NRHS, M, -1 ) )
            MAXWRK = MAX( MAXWRK, M*M+4*M+WLALSD )
          ELSE
!
!              Path 2 - remaining underdetermined cases.
!
            MAXWRK = 3*M + ( N+M )*ILAENV( 1, 'DGEBRD', ' ', M, N,      &
     &-1, -1 )
            MAXWRK = MAX( MAXWRK, 3*M+NRHS*                             &
     &ILAENV( 1, 'DORMBR', 'QLT', M, NRHS, N, -1 ) )
            MAXWRK = MAX( MAXWRK, 3*M+M*                                &
     &ILAENV( 1, 'DORMBR', 'PLN', N, NRHS, M, -1 ) )
            MAXWRK = MAX( MAXWRK, 3*M+WLALSD )
          END IF
          MINWRK = MAX( 3*M+NRHS, 3*M+M, 3*M+WLALSD )
        END IF
        MINWRK = MIN( MINWRK, MAXWRK )
        WORK( 1 ) = MAXWRK
        IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
          INFO = -12
        END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DGELSD', -INFO )
        RETURN
      ELSE IF( LQUERY ) THEN
        GO TO 10
      END IF
!
!     Quick return if possible.
!
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
        RANK = 0
        RETURN
      END IF
!
!     Get machine parameters.
!
      EPS = DLAMCH( 'P' )
      SFMIN = DLAMCH( 'S' )
      SMLNUM = SFMIN / EPS
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
!
!     Scale A if max entry outside range [SMLNUM,BIGNUM].
!
      ANRM = DLANGE( 'M', M, N, A, LDA, WORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
!
!        Scale matrix norm up to SMLNUM.
!
        CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
        IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN
!
!        Scale matrix norm down to BIGNUM.
!
        CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
        IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN
!
!        Matrix all zero. Return zero solution.
!
        CALL DLASET( 'F', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
        CALL DLASET( 'F', MINMN, 1, ZERO, ZERO, S, 1 )
        RANK = 0
        GO TO 10
      END IF
!
!     Scale B if max entry outside range [SMLNUM,BIGNUM].
!
      BNRM = DLANGE( 'M', M, NRHS, B, LDB, WORK )
      IBSCL = 0
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
!
!        Scale matrix norm up to SMLNUM.
!
        CALL DLASCL( 'G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO )
        IBSCL = 1
      ELSE IF( BNRM.GT.BIGNUM ) THEN
!
!        Scale matrix norm down to BIGNUM.
!
        CALL DLASCL( 'G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO )
        IBSCL = 2
      END IF
!
!     If M < N make sure certain entries of B are zero.
!
      IF( M.LT.N )                                                      &
     &CALL DLASET( 'F', N-M, NRHS, ZERO, ZERO, B( M+1, 1 ), LDB )
!
!     Overdetermined case.
!
      IF( M.GE.N ) THEN
!
!        Path 1 - overdetermined or exactly determined.
!
        MM = M
        IF( M.GE.MNTHR ) THEN
!
!           Path 1a - overdetermined, with many more rows than columns.
!
          MM = N
          ITAU = 1
          NWORK = ITAU + N
!
!           Compute A=Q*R.
!           (Workspace: need 2*N, prefer N+N*NB)
!
          CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ),       &
     &LWORK-NWORK+1, INFO )
!
!           Multiply B by transpose(Q).
!           (Workspace: need N+NRHS, prefer N+NRHS*NB)
!
          CALL DORMQR( 'L', 'T', M, NRHS, N, A, LDA, WORK( ITAU ), B,   &
     &LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
!
!           Zero out below R.
!
          IF( N.GT.1 ) THEN
            CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA )
          END IF
        END IF
!
        IE = 1
        ITAUQ = IE + N
        ITAUP = ITAUQ + N
        NWORK = ITAUP + N
!
!        Bidiagonalize R in A.
!        (Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB)
!
        CALL DGEBRD( MM, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ),       &
     &WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1,                      &
     &INFO )
!
!        Multiply B by transpose of left bidiagonalizing vectors of R.
!        (Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB)
!
        CALL DORMBR( 'Q', 'L', 'T', MM, NRHS, N, A, LDA, WORK( ITAUQ ), &
     &B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
!
!        Solve the bidiagonal least squares problem.
!
        CALL DLALSD( 'U', SMLSIZ, N, NRHS, S, WORK( IE ), B, LDB,       &
     &RCOND, RANK, WORK( NWORK ), IWORK, INFO )
        IF( INFO.NE.0 ) THEN
          GO TO 10
        END IF
!
!        Multiply B by right bidiagonalizing vectors of R.
!
        CALL DORMBR( 'P', 'L', 'N', N, NRHS, N, A, LDA, WORK( ITAUP ),  &
     &B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
!
      ELSE IF( N.GE.MNTHR .AND. LWORK.GE.4*M+M*M+                       &
     &MAX( M, 2*M-4, NRHS, N-3*M, WLALSD ) ) THEN
!
!        Path 2a - underdetermined, with many more columns than rows
!        and sufficient workspace for an efficient algorithm.
!
        LDWORK = M
        IF( LWORK.GE.MAX( 4*M+M*LDA+MAX( M, 2*M-4, NRHS, N-3*M ),       &
     &M*LDA+M+M*NRHS, 4*M+M*LDA+WLALSD ) )LDWORK = LDA
        ITAU = 1
        NWORK = M + 1
!
!        Compute A=L*Q.
!        (Workspace: need 2*M, prefer M+M*NB)
!
        CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ),         &
     &LWORK-NWORK+1, INFO )
        IL = NWORK
!
!        Copy L to WORK(IL), zeroing out above its diagonal.
!
        CALL DLACPY( 'L', M, M, A, LDA, WORK( IL ), LDWORK )
        CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, WORK( IL+LDWORK ),      &
     &LDWORK )
        IE = IL + LDWORK*M
        ITAUQ = IE + M
        ITAUP = ITAUQ + M
        NWORK = ITAUP + M
!
!        Bidiagonalize L in WORK(IL).
!        (Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB)
!
        CALL DGEBRD( M, M, WORK( IL ), LDWORK, S, WORK( IE ),           &
     &WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ),                      &
     &LWORK-NWORK+1, INFO )
!
!        Multiply B by transpose of left bidiagonalizing vectors of L.
!        (Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)
!
        CALL DORMBR( 'Q', 'L', 'T', M, NRHS, M, WORK( IL ), LDWORK,     &
     &WORK( ITAUQ ), B, LDB, WORK( NWORK ),                             &
     &LWORK-NWORK+1, INFO )
!
!        Solve the bidiagonal least squares problem.
!
        CALL DLALSD( 'U', SMLSIZ, M, NRHS, S, WORK( IE ), B, LDB,       &
     &RCOND, RANK, WORK( NWORK ), IWORK, INFO )
        IF( INFO.NE.0 ) THEN
          GO TO 10
        END IF
!
!        Multiply B by right bidiagonalizing vectors of L.
!
        CALL DORMBR( 'P', 'L', 'N', M, NRHS, M, WORK( IL ), LDWORK,     &
     &WORK( ITAUP ), B, LDB, WORK( NWORK ),                             &
     &LWORK-NWORK+1, INFO )
!
!        Zero out below first M rows of B.
!
        CALL DLASET( 'F', N-M, NRHS, ZERO, ZERO, B( M+1, 1 ), LDB )
        NWORK = ITAU + M
!
!        Multiply transpose(Q) by B.
!        (Workspace: need M+NRHS, prefer M+NRHS*NB)
!
        CALL DORMLQ( 'L', 'T', N, NRHS, M, A, LDA, WORK( ITAU ), B,     &
     &LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
!
      ELSE
!
!        Path 2 - remaining underdetermined cases.
!
        IE = 1
        ITAUQ = IE + M
        ITAUP = ITAUQ + M
        NWORK = ITAUP + M
!
!        Bidiagonalize A.
!        (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)
!
        CALL DGEBRD( M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ),        &
     &WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1,                      &
     &INFO )
!
!        Multiply B by transpose of left bidiagonalizing vectors.
!        (Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB)
!
        CALL DORMBR( 'Q', 'L', 'T', M, NRHS, N, A, LDA, WORK( ITAUQ ),  &
     &B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
!
!        Solve the bidiagonal least squares problem.
!
        CALL DLALSD( 'L', SMLSIZ, M, NRHS, S, WORK( IE ), B, LDB,       &
     &RCOND, RANK, WORK( NWORK ), IWORK, INFO )
        IF( INFO.NE.0 ) THEN
          GO TO 10
        END IF
!
!        Multiply B by right bidiagonalizing vectors of A.
!
        CALL DORMBR( 'P', 'L', 'N', N, NRHS, M, A, LDA, WORK( ITAUP ),  &
     &B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
!
      END IF
!
!     Undo scaling.
!
      IF( IASCL.EQ.1 ) THEN
        CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO )
        CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN,       &
     &INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
        CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO )
        CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN,       &
     &INFO )
      END IF
      IF( IBSCL.EQ.1 ) THEN
        CALL DLASCL( 'G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO )
      ELSE IF( IBSCL.EQ.2 ) THEN
        CALL DLASCL( 'G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO )
      END IF
!
   10 CONTINUE
      WORK( 1 ) = MAXWRK
      RETURN
!
!     End of DGELSD
!
      END
#ifdef _G95
      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
      implicit none
!...
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
!...
      BETA = RADIX(0.0D0)
      T = DIGITS(0.0D0)
      RND = .TRUE.
      IEEE1 = .TRUE.

      END
#endif
      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
!     ..
!
!  Purpose
!  =======
!
!  DLAMC2 determines the machine parameters specified in its argument
!  list.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  EPS     (output) DOUBLE PRECISION
!          The smallest positive number such that
!
!             fl( 1.0 - EPS ) .LT. 1.0,
!
!          where fl denotes the computed value.
!
!  EMIN    (output) INTEGER
!          The minimum exponent before (gradual) underflow occurs.
!
!  RMIN    (output) DOUBLE PRECISION
!          The smallest normalized number for the machine, given by
!          BASE**( EMIN - 1 ), where  BASE  is the floating point value
!          of BETA.
!
!  EMAX    (output) INTEGER
!          The maximum exponent before overflow occurs.
!
!  RMAX    (output) DOUBLE PRECISION
!          The largest positive number for the machine, given by
!          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
!          value of BETA.
!
!  Further Details
!  ===============
!
!  The computation of  EPS  is based on a routine PARANOIA by
!  W. Kahan of the University of California at Berkeley.
!
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,      &
     &NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE, &
     &SIXTH, SMALL, THIRD, TWO, ZERO
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAMC1, DLAMC4, DLAMC5
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,&
     &LRMIN, LT
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
        ZERO = 0
        ONE = 1
        TWO = 2
!
!        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
!        BETA, T, RND, EPS, EMIN and RMIN.
!
!        Throughout this routine  we use the function  DLAMC3  to ensure
!        that relevant values are stored  and not held in registers,  or
!        are not affected by optimizers.
!
!        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
!
        CALL DLAMC1( LBETA, LT, LRND, LIEEE1 )
!
!        Start to find EPS.
!
        B = LBETA
        A = B**( -LT )
        LEPS = A
!
!        Try some tricks to see whether or not this is the correct  EPS.
!
        B = TWO / 3
        HALF = ONE / 2
        SIXTH = DLAMC3( B, -HALF )
        THIRD = DLAMC3( SIXTH, SIXTH )
        B = DLAMC3( THIRD, -HALF )
        B = DLAMC3( B, SIXTH )
        B = ABS( B )
        IF( B.LT.LEPS )                                                 &
     &B = LEPS
!
        LEPS = 1
!
!+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10   CONTINUE
        IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
          LEPS = B
          C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
          C = DLAMC3( HALF, -C )
          B = DLAMC3( HALF, C )
          C = DLAMC3( HALF, -B )
          B = DLAMC3( HALF, C )
          GO TO 10
        END IF
!+       END WHILE
!
        IF( A.LT.LEPS )                                                 &
     &LEPS = A
!
!        Computation of EPS complete.
!
!        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
!        Keep dividing  A by BETA until (gradual) underflow occurs. This
!        is detected when we cannot recover the previous A.
!
        RBASE = ONE / LBETA
        SMALL = ONE
        DO I = 1, 3
          SMALL = DLAMC3( SMALL*RBASE, ZERO )
        enddo
        A = DLAMC3( ONE, SMALL )
        CALL DLAMC4( NGPMIN, ONE, LBETA )
        CALL DLAMC4( NGNMIN, -ONE, LBETA )
        CALL DLAMC4( GPMIN, A, LBETA )
        CALL DLAMC4( GNMIN, -A, LBETA )
        IEEE = .FALSE.
!
        IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
          IF( NGPMIN.EQ.GPMIN ) THEN
            LEMIN = NGPMIN
!            ( Non twos-complement machines, no gradual underflow;
!              e.g.,  VAX )
          ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
            LEMIN = NGPMIN - 1 + LT
            IEEE = .TRUE.
!            ( Non twos-complement machines, with gradual underflow;
!              e.g., IEEE standard followers )
          ELSE
            LEMIN = MIN( NGPMIN, GPMIN )
!            ( A guess; no known machine )
            IWARN = .TRUE.
          END IF
!
        ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
          IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
            LEMIN = MAX( NGPMIN, NGNMIN )
!            ( Twos-complement machines, no gradual underflow;
!              e.g., CYBER 205 )
          ELSE
            LEMIN = MIN( NGPMIN, NGNMIN )
!            ( A guess; no known machine )
            IWARN = .TRUE.
          END IF
!
        ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.                    &
     &( GPMIN.EQ.GNMIN ) ) THEN
          IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
            LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
!            ( Twos-complement machines with gradual underflow;
!              no known machine )
          ELSE
            LEMIN = MIN( NGPMIN, NGNMIN )
!            ( A guess; no known machine )
            IWARN = .TRUE.
          END IF
!
        ELSE
          LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
!         ( A guess; no known machine )
          IWARN = .TRUE.
        END IF
        FIRST = .FALSE.
!**
! Comment out this if block if EMIN is ok
        IF( IWARN ) THEN
          FIRST = .TRUE.
          WRITE( 6, FMT = 9999 )LEMIN
        END IF
!**
!
!        Assume IEEE arithmetic if we found denormalised  numbers above,
!        or if arithmetic seems to round in the  IEEE style,  determined
!        in routine DLAMC1. A true IEEE machine should have both  things
!        true; however, faulty machines may have one or the other.
!
        IEEE = IEEE .OR. LIEEE1
!
!        Compute  RMIN by successive division by  BETA. We could compute
!        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
!        this computation.
!
        LRMIN = 1
        DO I = 1, 1 - LEMIN
          LRMIN = DLAMC3( LRMIN*RBASE, ZERO )
        enddo
!
!        Finally, call DLAMC5 to compute EMAX and RMAX.
!
        CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
!
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
!
      RETURN
!
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',        &
     &'  EMIN = ', I8, /                                                &
     &' If, after inspection, the value EMIN looks',                    &
     &' acceptable please comment out ',                                &
     &/ ' the IF block as marked within the code of routine',           &
     &' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )
!
!     End of DLAMC2
!
      END
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
!     ..
!
!  Purpose
!  =======
!
!  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
!  the addition of  A  and  B ,  for use in situations where optimizers
!  might hold one of these in a register.
!
!  Arguments
!  =========
!
!  A       (input) DOUBLE PRECISION
!  B       (input) DOUBLE PRECISION
!          The values A and B.
!
! =====================================================================
!
!     .. Executable Statements ..
!
      DLAMC3 = A + B
!
      RETURN
!
!     End of DLAMC3
!
      END
      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          CMACH
!     ..
!
!  Purpose
!  =======
!
!  DLAMCH determines double precision machine parameters.
!
!  Arguments
!  =========
!
!  CMACH   (input) CHARACTER*1
!          Specifies the value to be returned by DLAMCH:
!          = 'E' or 'e',   DLAMCH := eps
!          = 'S' or 's ,   DLAMCH := sfmin
!          = 'B' or 'b',   DLAMCH := base
!          = 'P' or 'p',   DLAMCH := eps*base
!          = 'N' or 'n',   DLAMCH := t
!          = 'R' or 'r',   DLAMCH := rnd
!          = 'M' or 'm',   DLAMCH := emin
!          = 'U' or 'u',   DLAMCH := rmin
!          = 'L' or 'l',   DLAMCH := emax
!          = 'O' or 'o',   DLAMCH := rmax
!
!          where
!
!          eps   = relative machine precision
!          sfmin = safe minimum, such that 1/sfmin does not overflow
!          base  = base of the machine
!          prec  = eps*base
!          t     = number of (base) digits in the mantissa
!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!          emin  = minimum exponent before (gradual) underflow
!          rmin  = underflow threshold - base**(emin-1)
!          emax  = largest exponent before overflow
!          rmax  = overflow threshold  - (base**emax)*(1-eps)
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      DOUBLE PRECISION   BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,&
     &RND, SFMIN, SMALL, T
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAMC2
!     ..
!     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,   &
     &EMAX, RMAX, PREC
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
        CALL DLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
        BASE = BETA
        T = IT
        IF( LRND ) THEN
          RND = ONE
          EPS = ( BASE**( 1-IT ) ) / 2
        ELSE
          RND = ZERO
          EPS = BASE**( 1-IT )
        END IF
        PREC = EPS*BASE
        EMIN = IMIN
        EMAX = IMAX
        SFMIN = RMIN
        SMALL = ONE / RMAX
        IF( SMALL.GE.SFMIN ) THEN
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
          SFMIN = SMALL*( ONE+EPS )
        END IF
      END IF
!
      IF( LSAME( CMACH, 'E' ) ) THEN
        RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
        RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
        RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
        RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
        RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
        RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
        RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
        RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
        RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
        RMACH = RMAX
      END IF
!
      DLAMCH = RMACH
      FIRST  = .FALSE.
      RETURN
!
!     End of DLAMCH
!
      END
      DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DLANGE  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  real matrix A.
!
!  Description
!  ===========
!
!  DLANGE returns the value
!
!     DLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in DLANGE as described
!          above.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.  When M = 0,
!          DLANGE is set to zero.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.  When N = 0,
!          DLANGE is set to zero.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The m by n matrix A.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(M,1).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
!          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
!          referenced.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   SCALE, SUM, VALUE
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASSQ
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      IF( MIN( M, N ).EQ.0 ) THEN
        VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
        VALUE = ZERO
        DO J = 1, N
          DO I = 1, M
            VALUE = MAX( VALUE, ABS( A( I, J ) ) )
          enddo
        enddo
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
!
!        Find norm1(A).
!
        VALUE = ZERO
        DO J = 1, N
          SUM = ZERO
          DO I = 1, M
            SUM = SUM + ABS( A( I, J ) )
          enddo
          VALUE = MAX( VALUE, SUM )
        enddo
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
!
!        Find normI(A).
!
        DO I = 1, M
          WORK( I ) = ZERO
        enddo
        DO J = 1, N
          DO I = 1, M
            WORK( I ) = WORK( I ) + ABS( A( I, J ) )
          enddo
        enddo
        VALUE = ZERO
        DO I = 1, M
          VALUE = MAX( VALUE, WORK( I ) )
        enddo
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
        SCALE = ZERO
        SUM = ONE
        DO J = 1, N
          CALL DLASSQ( M, A( 1, J ), 1, SCALE, SUM )
        enddo
        VALUE = SCALE*SQRT( SUM )
      END IF
!
      DLANGE = VALUE
      RETURN
!
!     End of DLANGE
!
      END
      DOUBLE PRECISION FUNCTION DLANST( NORM, N, D, E )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
!     ..
!
!  Purpose
!  =======
!
!  DLANST  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  real symmetric tridiagonal matrix A.
!
!  Description
!  ===========
!
!  DLANST returns the value
!
!     DLANST = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in DLANST as described
!          above.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, DLANST is
!          set to zero.
!
!  D       (input) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of A.
!
!  E       (input) DOUBLE PRECISION array, dimension (N-1)
!          The (n-1) sub-diagonal or super-diagonal elements of A.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   ANORM, SCALE, SUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASSQ
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
      IF( N.LE.0 ) THEN
        ANORM = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
        ANORM = ABS( D( N ) )
        DO I = 1, N - 1
          ANORM = MAX( ANORM, ABS( D( I ) ) )
          ANORM = MAX( ANORM, ABS( E( I ) ) )
        enddo
      ELSE IF( LSAME( NORM, 'O' ) .OR. NORM.EQ.'1' .OR.                 &
     &LSAME( NORM, 'I' ) ) THEN
!
!        Find norm1(A).
!
        IF( N.EQ.1 ) THEN
          ANORM = ABS( D( 1 ) )
        ELSE
          ANORM = MAX( ABS( D( 1 ) )+ABS( E( 1 ) ),                     &
     &ABS( E( N-1 ) )+ABS( D( N ) ) )
          DO I = 2, N - 1
            ANORM = MAX( ANORM, ABS( D( I ) )+ABS( E( I ) )+            &
     &ABS( E( I-1 ) ) )
          enddo
        END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
        SCALE = ZERO
        SUM = ONE
        IF( N.GT.1 ) THEN
          CALL DLASSQ( N-1, E, 1, SCALE, SUM )
          SUM = 2*SUM
        END IF
        CALL DLASSQ( N, D, 1, SCALE, SUM )
        ANORM = SCALE*SQRT( SUM )
      END IF
!
      DLANST = ANORM
      RETURN
!
!     End of DLANST
!
      END
      DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y
!     ..
!
!  Purpose
!  =======
!
!  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
!  overflow.
!
!  Arguments
!  =========
!
!  X       (input) DOUBLE PRECISION
!  Y       (input) DOUBLE PRECISION
!          X and Y specify the values x and y.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, Z
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      XABS = ABS( X )
      YABS = ABS( Y )
      W = MAX( XABS, YABS )
      Z = MIN( XABS, YABS )
      IF( Z.EQ.ZERO ) THEN
        DLAPY2 = W
      ELSE
        DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
      END IF
      RETURN
!
!     End of DLAPY2
!
      END
      DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
      implicit none
!     .. Scalar Arguments ..
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION X(*)
!     ..
!
!  Purpose
!  =======
!
!  DNRM2 returns the euclidean norm of a vector via the function
!  name, so that
!
!     DNRM2 := sqrt( x'*x )
!
!
!  -- This version written on 25-October-1982.
!     Modified on 14-October-1993 to inline the call to DLASSQ.
!     Sven Hammarling, Nag Ltd.
!
!
!     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION ABSXI,NORM,SCALE,SSQ
      INTEGER IX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
!     ..
      IF (N.LT.1 .OR. INCX.LT.1) THEN
        NORM = ZERO
      ELSE IF (N.EQ.1) THEN
        NORM = ABS(X(1))
      ELSE
        SCALE = ZERO
        SSQ = ONE
!        The following loop is equivalent to this call to the LAPACK
!        auxiliary routine:
!        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
!
        DO IX = 1,1 + (N-1)*INCX,INCX
          IF (X(IX).NE.ZERO) THEN
            ABSXI = ABS(X(IX))
            IF (SCALE.LT.ABSXI) THEN
              SSQ = ONE + SSQ* (SCALE/ABSXI)**2
              SCALE = ABSXI
            ELSE
              SSQ = SSQ + (ABSXI/SCALE)**2
            END IF
          END IF
        enddo
        NORM = SCALE*SQRT(SSQ)
      END IF
!
      DNRM2 = NORM
      RETURN
!
!     End of DNRM2.
!
      END
      SUBROUTINE DLALSD( UPLO, SMLSIZ, N, NRHS, D, E, B, LDB, RCOND,    &
     &RANK, WORK, IWORK, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDB, N, NRHS, RANK, SMLSIZ
      DOUBLE PRECISION   RCOND
!     ..
!     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   B( LDB, * ), D( * ), E( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DLALSD uses the singular value decomposition of A to solve the least
!  squares problem of finding X to minimize the Euclidean norm of each
!  column of A*X-B, where A is N-by-N upper bidiagonal, and X and B
!  are N-by-NRHS. The solution X overwrites B.
!
!  The singular values of A smaller than RCOND times the largest
!  singular value are treated as zero in solving the least squares
!  problem; in this case a minimum norm solution is returned.
!  The actual singular values are returned in D in ascending order.
!
!  This code makes very mild assumptions about floating point
!  arithmetic. It will work on machines with a guard digit in
!  add/subtract, or on those binary machines without guard digits
!  which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2.
!  It could conceivably fail on hexadecimal or decimal machines
!  without guard digits, but we know of none.
!
!  Arguments
!  =========
!
!  UPLO   (input) CHARACTER*1
!         = 'U': D and E define an upper bidiagonal matrix.
!         = 'L': D and E define a  lower bidiagonal matrix.
!
!  SMLSIZ (input) INTEGER
!         The maximum size of the subproblems at the bottom of the
!         computation tree.
!
!  N      (input) INTEGER
!         The dimension of the  bidiagonal matrix.  N >= 0.
!
!  NRHS   (input) INTEGER
!         The number of columns of B. NRHS must be at least 1.
!
!  D      (input/output) DOUBLE PRECISION array, dimension (N)
!         On entry D contains the main diagonal of the bidiagonal
!         matrix. On exit, if INFO = 0, D contains its singular values.
!
!  E      (input/output) DOUBLE PRECISION array, dimension (N-1)
!         Contains the super-diagonal entries of the bidiagonal matrix.
!         On exit, E has been destroyed.
!
!  B      (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!         On input, B contains the right hand sides of the least
!         squares problem. On output, B contains the solution X.
!
!  LDB    (input) INTEGER
!         The leading dimension of B in the calling subprogram.
!         LDB must be at least max(1,N).
!
!  RCOND  (input) DOUBLE PRECISION
!         The singular values of A less than or equal to RCOND times
!         the largest singular value are treated as zero in solving
!         the least squares problem. If RCOND is negative,
!         machine precision is used instead.
!         For example, if diag(S)*X=B were the least squares problem,
!         where diag(S) is a diagonal matrix of singular values, the
!         solution would be X(i) = B(i) / S(i) if S(i) is greater than
!         RCOND*max(S), and X(i) = 0 if S(i) is less than or equal to
!         RCOND*max(S).
!
!  RANK   (output) INTEGER
!         The number of singular values of A greater than RCOND times
!         the largest singular value.
!
!  WORK   (workspace) DOUBLE PRECISION array, dimension at least
!         (9*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2),
!         where NLVL = max(0, INT(log_2 (N/(SMLSIZ+1))) + 1).
!
!  IWORK  (workspace) INTEGER array, dimension at least
!         (3*N*NLVL + 11*N)
!
!  INFO   (output) INTEGER
!         = 0:  successful exit.
!         < 0:  if INFO = -i, the i-th argument had an illegal value.
!         > 0:  The algorithm failed to compute an singular value while
!               working on the submatrix lying in rows and columns
!               INFO/(N+1) through MOD(INFO,N+1).
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ming Gu and Ren-Cang Li, Computer Science Division, University of
!       California at Berkeley, USA
!     Osni Marques, LBNL/NERSC, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            BX, BXST, C, DIFL, DIFR, GIVCOL, GIVNUM,       &
     &GIVPTR, I, ICMPQ1, ICMPQ2, IWK, J, K, NLVL,                       &
     &NM1, NSIZE, NSUB, NWORK, PERM, POLES, S, SIZEI,                   &
     &SMLSZP, SQRE, ST, ST1, U, VT, Z
      DOUBLE PRECISION   CS, EPS, ORGNRM, R, RCND, SN, TOL
!     ..
!     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH, DLANST
      EXTERNAL           IDAMAX, DLAMCH, DLANST
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPL, DGEMM, DLACPY, DLALSA, DLARTG, DLASCL,  &
     &DLASDA, DLASDQ, DLASET, DLASRT, DROT, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, LOG, SIGN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( N.LT.0 ) THEN
        INFO = -3
      ELSE IF( NRHS.LT.1 ) THEN
        INFO = -4
      ELSE IF( ( LDB.LT.1 ) .OR. ( LDB.LT.N ) ) THEN
        INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DLALSD', -INFO )
        RETURN
      END IF
!
      EPS = DLAMCH( 'Epsilon' )
!
!     Set up the tolerance.
!
      IF( ( RCOND.LE.ZERO ) .OR. ( RCOND.GE.ONE ) ) THEN
        RCND = EPS
      ELSE
        RCND = RCOND
      END IF
!
      RANK = 0
!
!     Quick return if possible.
!
      IF( N.EQ.0 ) THEN
        RETURN
      ELSE IF( N.EQ.1 ) THEN
        IF( D( 1 ).EQ.ZERO ) THEN
          CALL DLASET( 'A', 1, NRHS, ZERO, ZERO, B, LDB )
        ELSE
          RANK = 1
          CALL DLASCL( 'G', 0, 0, D( 1 ), ONE, 1, NRHS, B, LDB, INFO )
          D( 1 ) = ABS( D( 1 ) )
        END IF
        RETURN
      END IF
!
!     Rotate the matrix if it is lower bidiagonal.
!
      IF( UPLO.EQ.'L' ) THEN
        DO I = 1, N - 1
          CALL DLARTG( D( I ), E( I ), CS, SN, R )
          D( I ) = R
          E( I ) = SN*D( I+1 )
          D( I+1 ) = CS*D( I+1 )
          IF( NRHS.EQ.1 ) THEN
            CALL DROT( 1, B( I, 1 ), 1, B( I+1, 1 ), 1, CS, SN )
          ELSE
            WORK( I*2-1 ) = CS
            WORK( I*2 ) = SN
          END IF
        enddo
        IF( NRHS.GT.1 ) THEN
          DO I = 1, NRHS
            DO J = 1, N - 1
              CS = WORK( J*2-1 )
              SN = WORK( J*2 )
              CALL DROT( 1, B( J, I ), 1, B( J+1, I ), 1, CS, SN )
            enddo
          enddo
        END IF
      END IF
!
!     Scale.
!
      NM1 = N - 1
      ORGNRM = DLANST( 'M', N, D, E )
      IF( ORGNRM.EQ.ZERO ) THEN
        CALL DLASET( 'A', N, NRHS, ZERO, ZERO, B, LDB )
        RETURN
      END IF
!
      CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, N, 1, D, N, INFO )
      CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, NM1, 1, E, NM1, INFO )
!
!     If N is smaller than the minimum divide size SMLSIZ, then solve
!     the problem with another solver.
!
      IF( N.LE.SMLSIZ ) THEN
        NWORK = 1 + N*N
        CALL DLASET( 'A', N, N, ZERO, ONE, WORK, N )
        CALL DLASDQ( 'U', 0, N, N, 0, NRHS, D, E, WORK, N, WORK, N, B,  &
     &LDB, WORK( NWORK ), INFO )
        IF( INFO.NE.0 ) THEN
          RETURN
        END IF
        TOL = RCND*ABS( D( IDAMAX( N, D, 1 ) ) )
        DO I = 1, N
          IF( D( I ).LE.TOL ) THEN
            CALL DLASET( 'A', 1, NRHS, ZERO, ZERO, B( I, 1 ), LDB )
          ELSE
            CALL DLASCL( 'G', 0, 0, D( I ), ONE, 1, NRHS, B( I, 1 ),    &
     &LDB, INFO )
            RANK = RANK + 1
          END IF
        enddo
        CALL DGEMM( 'T', 'N', N, NRHS, N, ONE, WORK, N, B, LDB, ZERO,   &
     &WORK( NWORK ), N )
        CALL DLACPY( 'A', N, NRHS, WORK( NWORK ), N, B, LDB )
!
!        Unscale.
!
        CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO )
        CALL DLASRT( 'D', N, D, INFO )
        CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, N, NRHS, B, LDB, INFO )
!
        RETURN
      END IF
!
!     Book-keeping and setting up some constants.
!
      NLVL = INT( LOG( DBLE( N ) / DBLE( SMLSIZ+1 ) ) / LOG( TWO ) ) + 1
!
      SMLSZP = SMLSIZ + 1
!
      U = 1
      VT = 1 + SMLSIZ*N
      DIFL = VT + SMLSZP*N
      DIFR = DIFL + NLVL*N
      Z = DIFR + NLVL*N*2
      C = Z + NLVL*N
      S = C + N
      POLES = S + N
      GIVNUM = POLES + 2*NLVL*N
      BX = GIVNUM + 2*NLVL*N
      NWORK = BX + N*NRHS
!
      SIZEI = 1 + N
      K = SIZEI + N
      GIVPTR = K + N
      PERM = GIVPTR + N
      GIVCOL = PERM + NLVL*N
      IWK = GIVCOL + NLVL*N*2
!
      ST = 1
      SQRE = 0
      ICMPQ1 = 1
      ICMPQ2 = 0
      NSUB = 0
!
      DO I = 1, N
        IF( ABS( D( I ) ).LT.EPS ) THEN
          D( I ) = SIGN( EPS, D( I ) )
        END IF
      enddo
!
      DO I = 1, NM1
        IF( ( ABS( E( I ) ).LT.EPS ) .OR. ( I.EQ.NM1 ) ) THEN
          NSUB = NSUB + 1
          IWORK( NSUB ) = ST
!
!           Subproblem found. First determine its size and then
!           apply divide and conquer on it.
!
          IF( I.LT.NM1 ) THEN
!
!              A subproblem with E(I) small for I < NM1.
!
            NSIZE = I - ST + 1
            IWORK( SIZEI+NSUB-1 ) = NSIZE
          ELSE IF( ABS( E( I ) ).GE.EPS ) THEN
!
!              A subproblem with E(NM1) not too small but I = NM1.
!
            NSIZE = N - ST + 1
            IWORK( SIZEI+NSUB-1 ) = NSIZE
          ELSE
!
!              A subproblem with E(NM1) small. This implies an
!              1-by-1 subproblem at D(N), which is not solved
!              explicitly.
!
            NSIZE = I - ST + 1
            IWORK( SIZEI+NSUB-1 ) = NSIZE
            NSUB = NSUB + 1
            IWORK( NSUB ) = N
            IWORK( SIZEI+NSUB-1 ) = 1
            CALL DCOPL( NRHS, B( N, 1 ), LDB, WORK( BX+NM1 ), N )
          END IF
          ST1 = ST - 1
          IF( NSIZE.EQ.1 ) THEN
!
!              This is a 1-by-1 subproblem and is not solved
!              explicitly.
!
            CALL DCOPL( NRHS, B( ST, 1 ), LDB, WORK( BX+ST1 ), N )
          ELSE IF( NSIZE.LE.SMLSIZ ) THEN
!
!              This is a small subproblem and is solved by DLASDQ.
!
            CALL DLASET( 'A', NSIZE, NSIZE, ZERO, ONE,                  &
     &WORK( VT+ST1 ), N )
            CALL DLASDQ( 'U', 0, NSIZE, NSIZE, 0, NRHS, D( ST ),        &
     &E( ST ), WORK( VT+ST1 ), N, WORK( NWORK ),                        &
     &N, B( ST, 1 ), LDB, WORK( NWORK ), INFO )
            IF( INFO.NE.0 ) THEN
              RETURN
            END IF
            CALL DLACPY( 'A', NSIZE, NRHS, B( ST, 1 ), LDB,             &
     &WORK( BX+ST1 ), N )
          ELSE
!
!              A large problem. Solve it using divide and conquer.
!
            CALL DLASDA( ICMPQ1, SMLSIZ, NSIZE, SQRE, D( ST ),          &
     &E( ST ), WORK( U+ST1 ), N, WORK( VT+ST1 ),                        &
     &IWORK( K+ST1 ), WORK( DIFL+ST1 ),                                 &
     &WORK( DIFR+ST1 ), WORK( Z+ST1 ),                                  &
     &WORK( POLES+ST1 ), IWORK( GIVPTR+ST1 ),                           &
     &IWORK( GIVCOL+ST1 ), N, IWORK( PERM+ST1 ),                        &
     &WORK( GIVNUM+ST1 ), WORK( C+ST1 ),                                &
     &WORK( S+ST1 ), WORK( NWORK ), IWORK( IWK ),                       &
     &INFO )
            IF( INFO.NE.0 ) THEN
              RETURN
            END IF
            BXST = BX + ST1
            CALL DLALSA( ICMPQ2, SMLSIZ, NSIZE, NRHS, B( ST, 1 ),       &
     &LDB, WORK( BXST ), N, WORK( U+ST1 ), N,                           &
     &WORK( VT+ST1 ), IWORK( K+ST1 ),                                   &
     &WORK( DIFL+ST1 ), WORK( DIFR+ST1 ),                               &
     &WORK( Z+ST1 ), WORK( POLES+ST1 ),                                 &
     &IWORK( GIVPTR+ST1 ), IWORK( GIVCOL+ST1 ), N,                      &
     &IWORK( PERM+ST1 ), WORK( GIVNUM+ST1 ),                            &
     &WORK( C+ST1 ), WORK( S+ST1 ), WORK( NWORK ),                      &
     &IWORK( IWK ), INFO )
            IF( INFO.NE.0 ) THEN
              RETURN
            END IF
          END IF
          ST = I + 1
        END IF
      enddo
!
!     Apply the singular values and treat the tiny ones as zero.
!
      TOL = RCND*ABS( D( IDAMAX( N, D, 1 ) ) )
!
      DO I = 1, N
!
!        Some of the elements in D can be negative because 1-by-1
!        subproblems were not solved explicitly.
!
        IF( ABS( D( I ) ).LE.TOL ) THEN
          CALL DLASET( 'A', 1, NRHS, ZERO, ZERO, WORK( BX+I-1 ), N )
        ELSE
          RANK = RANK + 1
          CALL DLASCL( 'G', 0, 0, D( I ), ONE, 1, NRHS,                 &
     &WORK( BX+I-1 ), N, INFO )
        END IF
        D( I ) = ABS( D( I ) )
      enddo
!
!     Now apply back the right singular vectors.
!
      ICMPQ2 = 1
      DO I = 1, NSUB
        ST = IWORK( I )
        ST1 = ST - 1
        NSIZE = IWORK( SIZEI+I-1 )
        BXST = BX + ST1
        IF( NSIZE.EQ.1 ) THEN
          CALL DCOPL( NRHS, WORK( BXST ), N, B( ST, 1 ), LDB )
        ELSE IF( NSIZE.LE.SMLSIZ ) THEN
          CALL DGEMM( 'T', 'N', NSIZE, NRHS, NSIZE, ONE,                &
     &WORK( VT+ST1 ), N, WORK( BXST ), N, ZERO,                         &
     &B( ST, 1 ), LDB )
        ELSE
          CALL DLALSA( ICMPQ2, SMLSIZ, NSIZE, NRHS, WORK( BXST ), N,    &
     &B( ST, 1 ), LDB, WORK( U+ST1 ), N,                                &
     &WORK( VT+ST1 ), IWORK( K+ST1 ),                                   &
     &WORK( DIFL+ST1 ), WORK( DIFR+ST1 ),                               &
     &WORK( Z+ST1 ), WORK( POLES+ST1 ),                                 &
     &IWORK( GIVPTR+ST1 ), IWORK( GIVCOL+ST1 ), N,                      &
     &IWORK( PERM+ST1 ), WORK( GIVNUM+ST1 ),                            &
     &WORK( C+ST1 ), WORK( S+ST1 ), WORK( NWORK ),                      &
     &IWORK( IWK ), INFO )
          IF( INFO.NE.0 ) THEN
            RETURN
          END IF
        END IF
      enddo
!
!     Unscale and sort the singular values.
!
      CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO )
      CALL DLASRT( 'D', N, D, INFO )
      CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, N, NRHS, B, LDB, INFO )
!
      RETURN
!
!     End of DLALSD
!
      END
      SUBROUTINE DLAMC4( EMIN, START, BASE )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
!     ..
!
!  Purpose
!  =======
!
!  DLAMC4 is a service routine for DLAMC2.
!
!  Arguments
!  =========
!
!  EMIN    (output) INTEGER
!          The minimum exponent before (gradual) underflow, computed by
!          setting A = START and dividing by BASE until the previous A
!          can not be recovered.
!
!  START   (input) DOUBLE PRECISION
!          The starting point for determining EMIN.
!
!  BASE    (input) INTEGER
!          The base of the machine.
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
!     ..
!     .. Executable Statements ..
!
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
!+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
!    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.         &
     &( D2.EQ.A ) ) THEN
        EMIN = EMIN - 1
        A = B1
        B1 = DLAMC3( A / BASE, ZERO )
        C1 = DLAMC3( B1*BASE, ZERO )
        D1 = ZERO
        DO I = 1, BASE
          D1 = D1 + B1
        enddo
        B2 = DLAMC3( A*RBASE, ZERO )
        C2 = DLAMC3( B2 / RBASE, ZERO )
        D2 = ZERO
        DO I = 1, BASE
          D2 = D2 + B2
        enddo
        GO TO 10
      END IF
!+    END WHILE
!
      RETURN
!
!     End of DLAMC4
!
      END
      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
!     ..
!
!  Purpose
!  =======
!
!  DLAMC5 attempts to compute RMAX, the largest machine floating-point
!  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
!  approximately to a power of 2.  It will fail on machines where this
!  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
!  EMAX = 28718).  It will also fail if the value supplied for EMIN is
!  too large (i.e. too close to zero), probably with overflow.
!
!  Arguments
!  =========
!
!  BETA    (input) INTEGER
!          The base of floating-point arithmetic.
!
!  P       (input) INTEGER
!          The number of base BETA digits in the mantissa of a
!          floating-point value.
!
!  EMIN    (input) INTEGER
!          The minimum exponent before (gradual) underflow.
!
!  IEEE    (input) LOGICAL
!          A logical flag specifying whether or not the arithmetic
!          system is thought to comply with the IEEE standard.
!
!  EMAX    (output) INTEGER
!          The largest exponent before overflow
!
!  RMAX    (output) DOUBLE PRECISION
!          The largest machine floating-point number.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MOD
!     ..
!     .. Executable Statements ..
!
!     First compute LEXP and UEXP, two powers of 2 that bound
!     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
!     approximately to the bound that is closest to abs(EMIN).
!     (EMAX is the exponent of the required number RMAX).
!
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
        LEXP = TRY
        EXBITS = EXBITS + 1
        GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
        UEXP = LEXP
      ELSE
        UEXP = TRY
        EXBITS = EXBITS + 1
      END IF
!
!     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
!     than or equal to EMIN. EXBITS is the number of bits needed to
!     store the exponent.
!
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
        EXPSUM = 2*LEXP
      ELSE
        EXPSUM = 2*UEXP
      END IF
!
!     EXPSUM is the exponent range, approximately equal to
!     EMAX - EMIN + 1 .
!
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
!
!     NBITS is the total number of bits needed to store a
!     floating-point number.
!
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
!
!        Either there are an odd number of bits used to store a
!        floating-point number, which is unlikely, or some bits are
!        not used in the representation of numbers, which is possible,
!        (e.g. Cray machines) or the mantissa has an implicit bit,
!        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
!        most likely. We have to assume the last alternative.
!        If this is true, then we need to reduce EMAX by one because
!        there must be some way of representing zero in an implicit-bit
!        system. On machines like Cray, we are reducing EMAX by one
!        unnecessarily.
!
        EMAX = EMAX - 1
      END IF
!
      IF( IEEE ) THEN
!
!        Assume we are on an IEEE machine which reserves one exponent
!        for infinity and NaN.
!
        EMAX = EMAX - 1
      END IF
!
!     Now create RMAX, the largest machine number, which should
!     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
!
!     First compute 1.0 - BETA**(-P), being careful that the
!     result is less than 1.0 .
!
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO I = 1, P
        Z = Z*RECBAS
        IF( Y.LT.ONE )                                                  &
     &OLDY = Y
        Y = DLAMC3( Y, Z )
      enddo
      IF( Y.GE.ONE )                                                    &
     &Y = OLDY
!
!     Now multiply by BETA**EMAX to get RMAX.
!
      DO I = 1, EMAX
        Y = DLAMC3( Y*BETA, ZERO )
      enddo
!
      RMAX = Y
      RETURN
!
!     End of DLAMC5
!
      END
      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
!     ..
!
!  Purpose
!  =======
!
!  ILAENV is called from the LAPACK routines to choose problem-dependent
!  parameters for the local environment.  See ISPEC for a description of
!  the parameters.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines (DEPRECATED)
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR method
!               for nonsymmetric eigenvalue problems (DEPRECATED)
!          = 9: maximum size of the subproblems at the bottom of the
!               computation tree in the divide-and-conquer algorithm
!               (used by xGELSD and xGESDD)
!          =10: ieee NaN arithmetic can be trusted not to trap
!          =11: infinity arithmetic can be trusted not to trap
!          12 <= ISPEC <= 16:
!               xHSEQR or one of its subroutines,
!               see IPARMQ for detailed explanation
!
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
! (ILAENV) (output) INTEGER
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB.LE.1 ) NB = MAX( 1, N )
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IC, IZ, NB, NBMIN, NX
      LOGICAL            CNAME, SNAME
      CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*6
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
!     ..
!     .. External Functions ..
      INTEGER            IEEECK, IPARMQ
      EXTERNAL           IEEECK, IPARMQ
!     ..
!     .. Executable Statements ..
!
      GO TO ( 10, 10, 10, 80, 90, 100, 110, 120,                        &
     &130, 140, 150, 160, 160, 160, 160, 160 )ISPEC
!
!     Invalid value for ISPEC
!
      ILAENV = -1
      RETURN
!
   10 CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1: 1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
!
!        ASCII character set
!
        IF( IC.GE.97 .AND. IC.LE.122 ) THEN
          SUBNAM( 1: 1 ) = CHAR( IC-32 )
          DO I = 2, 6
            IC = ICHAR( SUBNAM( I: I ) )
            IF( IC.GE.97 .AND. IC.LE.122 )                              &
     &SUBNAM( I: I ) = CHAR( IC-32 )
          enddo
        END IF
!
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
!
!        EBCDIC character set
!
        IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.                          &
     &( IC.GE.145 .AND. IC.LE.153 ) .OR.                                &
     &( IC.GE.162 .AND. IC.LE.169 ) ) THEN
          SUBNAM( 1: 1 ) = CHAR( IC+64 )
          DO I = 2, 6
            IC = ICHAR( SUBNAM( I: I ) )
            IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.                      &
     &( IC.GE.145 .AND. IC.LE.153 ) .OR.                                &
     &( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I:                         &
     &I ) = CHAR( IC+64 )
          enddo
        END IF
!
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
!
!        Prime machines:  ASCII+128
!
        IF( IC.GE.225 .AND. IC.LE.250 ) THEN
          SUBNAM( 1: 1 ) = CHAR( IC-32 )
          DO I = 2, 6
            IC = ICHAR( SUBNAM( I: I ) )
            IF( IC.GE.225 .AND. IC.LE.250 )                             &
     &SUBNAM( I: I ) = CHAR( IC-32 )
          enddo
        END IF
      END IF
!
      C1 = SUBNAM( 1: 1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )                                   &
     &RETURN
      C2 = SUBNAM( 2: 3 )
      C3 = SUBNAM( 4: 6 )
      C4 = C3( 2: 3 )
!
      GO TO ( 50, 60, 70 )ISPEC
!
 50   CONTINUE
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
      NB = 1
!
      IF( C2.EQ.'GE' ) THEN
        IF( C3.EQ.'TRF' ) THEN
          IF( SNAME ) THEN
            NB = 64
          ELSE
            NB = 64
          END IF
        ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.     &
     &C3.EQ.'QLF' ) THEN
          IF( SNAME ) THEN
            NB = 32
          ELSE
            NB = 32
          END IF
        ELSE IF( C3.EQ.'HRD' ) THEN
          IF( SNAME ) THEN
            NB = 32
          ELSE
            NB = 32
          END IF
        ELSE IF( C3.EQ.'BRD' ) THEN
          IF( SNAME ) THEN
            NB = 32
          ELSE
            NB = 32
          END IF
        ELSE IF( C3.EQ.'TRI' ) THEN
          IF( SNAME ) THEN
            NB = 64
          ELSE
            NB = 64
          END IF
        END IF
      ELSE IF( C2.EQ.'PO' ) THEN
        IF( C3.EQ.'TRF' ) THEN
          IF( SNAME ) THEN
            NB = 64
          ELSE
            NB = 64
          END IF
        END IF
      ELSE IF( C2.EQ.'SY' ) THEN
        IF( C3.EQ.'TRF' ) THEN
          IF( SNAME ) THEN
            NB = 64
          ELSE
            NB = 64
          END IF
        ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
          NB = 32
        ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
          NB = 64
        END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
        IF( C3.EQ.'TRF' ) THEN
          NB = 64
        ELSE IF( C3.EQ.'TRD' ) THEN
          NB = 32
        ELSE IF( C3.EQ.'GST' ) THEN
          NB = 64
        END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
        IF( C3( 1: 1 ).EQ.'G' ) THEN
          IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.    &
     &'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )            &
     &THEN
            NB = 32
          END IF
        ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
          IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.    &
     &'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )            &
     &THEN
            NB = 32
          END IF
        END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
        IF( C3( 1: 1 ).EQ.'G' ) THEN
          IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.    &
     &'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )            &
     &THEN
            NB = 32
          END IF
        ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
          IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.    &
     &'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )            &
     &THEN
            NB = 32
          END IF
        END IF
      ELSE IF( C2.EQ.'GB' ) THEN
        IF( C3.EQ.'TRF' ) THEN
          IF( SNAME ) THEN
            IF( N4.LE.64 ) THEN
              NB = 1
            ELSE
              NB = 32
            END IF
          ELSE
            IF( N4.LE.64 ) THEN
              NB = 1
            ELSE
              NB = 32
            END IF
          END IF
        END IF
      ELSE IF( C2.EQ.'PB' ) THEN
        IF( C3.EQ.'TRF' ) THEN
          IF( SNAME ) THEN
            IF( N2.LE.64 ) THEN
              NB = 1
            ELSE
              NB = 32
            END IF
          ELSE
            IF( N2.LE.64 ) THEN
              NB = 1
            ELSE
              NB = 32
            END IF
          END IF
        END IF
      ELSE IF( C2.EQ.'TR' ) THEN
        IF( C3.EQ.'TRI' ) THEN
          IF( SNAME ) THEN
            NB = 64
          ELSE
            NB = 64
          END IF
        END IF
      ELSE IF( C2.EQ.'LA' ) THEN
        IF( C3.EQ.'UUM' ) THEN
          IF( SNAME ) THEN
            NB = 64
          ELSE
            NB = 64
          END IF
        END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
        IF( C3.EQ.'EBZ' ) THEN
          NB = 1
        END IF
      END IF
      ILAENV = NB
      RETURN
!
 60   CONTINUE
!
!     ISPEC = 2:  minimum block size
!
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
        IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.   &
     &'QLF' ) THEN
          IF( SNAME ) THEN
            NBMIN = 2
          ELSE
            NBMIN = 2
          END IF
        ELSE IF( C3.EQ.'HRD' ) THEN
          IF( SNAME ) THEN
            NBMIN = 2
          ELSE
            NBMIN = 2
          END IF
        ELSE IF( C3.EQ.'BRD' ) THEN
          IF( SNAME ) THEN
            NBMIN = 2
          ELSE
            NBMIN = 2
          END IF
        ELSE IF( C3.EQ.'TRI' ) THEN
          IF( SNAME ) THEN
            NBMIN = 2
          ELSE
            NBMIN = 2
          END IF
        END IF
      ELSE IF( C2.EQ.'SY' ) THEN
        IF( C3.EQ.'TRF' ) THEN
          IF( SNAME ) THEN
            NBMIN = 8
          ELSE
            NBMIN = 8
          END IF
        ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
          NBMIN = 2
        END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
        IF( C3.EQ.'TRD' ) THEN
          NBMIN = 2
        END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
        IF( C3( 1: 1 ).EQ.'G' ) THEN
          IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.    &
     &'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )            &
     &THEN
            NBMIN = 2
          END IF
        ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
          IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.    &
     &'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )            &
     &THEN
            NBMIN = 2
          END IF
        END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
        IF( C3( 1: 1 ).EQ.'G' ) THEN
          IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.    &
     &'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )            &
     &THEN
            NBMIN = 2
          END IF
        ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
          IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.    &
     &'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )            &
     &THEN
            NBMIN = 2
          END IF
        END IF
      END IF
      ILAENV = NBMIN
      RETURN
!
 70   CONTINUE
!
!     ISPEC = 3:  crossover point
!
      NX = 0
      IF( C2.EQ.'GE' ) THEN
        IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.   &
     &'QLF' ) THEN
          IF( SNAME ) THEN
            NX = 128
          ELSE
            NX = 128
          END IF
        ELSE IF( C3.EQ.'HRD' ) THEN
          IF( SNAME ) THEN
            NX = 128
          ELSE
            NX = 128
          END IF
        ELSE IF( C3.EQ.'BRD' ) THEN
          IF( SNAME ) THEN
            NX = 128
          ELSE
            NX = 128
          END IF
        END IF
      ELSE IF( C2.EQ.'SY' ) THEN
        IF( SNAME .AND. C3.EQ.'TRD' ) THEN
          NX = 32
        END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
        IF( C3.EQ.'TRD' ) THEN
          NX = 32
        END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
        IF( C3( 1: 1 ).EQ.'G' ) THEN
          IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.    &
     &'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )            &
     &THEN
            NX = 128
          END IF
        END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
        IF( C3( 1: 1 ).EQ.'G' ) THEN
          IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.    &
     &'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )            &
     &THEN
            NX = 128
          END IF
        END IF
      END IF
      ILAENV = NX
      RETURN
!
 80   CONTINUE
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
      ILAENV = 6
      RETURN
!
 90   CONTINUE
!
!     ISPEC = 5:  minimum column dimension (not used)
!
      ILAENV = 2
      RETURN
!     
 100  CONTINUE
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
!
 110  CONTINUE
!
!     ISPEC = 7:  number of processors (not used)
!
      ILAENV = 1
      RETURN
!
  120 CONTINUE
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
      ILAENV = 50
      RETURN
!
 130  CONTINUE
!
!     ISPEC = 9:  maximum size of the subproblems at the bottom of the
!                 computation tree in the divide-and-conquer algorithm
!                 (used by xGELSD and xGESDD)
!
      ILAENV = 25
      RETURN
!
 140  CONTINUE
!
!     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
!
!     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
        ILAENV = IEEECK( 0, 0.0, 1.0 )
      END IF
      RETURN
!
 150  CONTINUE
!
!     ISPEC = 11: infinity arithmetic can be trusted not to trap
!
!     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
        ILAENV = IEEECK( 1, 0.0, 1.0 )
      END IF
      RETURN
!
 160  CONTINUE
!
!     12 <= ISPEC <= 16: xHSEQR or one of its subroutines.
!
      ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
      RETURN
!
!     End of ILAENV
!
      END
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
      IMPLICIT NONE
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
!     ..
!
!  Purpose
!  =======
!
!  IEEECK is called from the ILAENV to verify that Infinity and
!  possibly NaN arithmetic is safe (i.e. will not trap).
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies whether to test just for inifinity arithmetic
!          or whether to test for infinity and NaN arithmetic.
!          = 0: Verify infinity arithmetic only.
!          = 1: Verify infinity and NaN arithmetic.
!
!  ZERO    (input) REAL
!          Must contain the value 0.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  ONE     (input) REAL
!          Must contain the value 1.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  RETURN VALUE:  INTEGER
!          = 0:  Arithmetic failed to produce the correct answers
!          = 1:  Arithmetic produced the correct answers
!
!     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,    &
     &NEGZRO, NEWZRO, POSINF
!     ..
!     .. Executable Statements ..
      IEEECK = 1
!
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
        IEEECK = 0
        RETURN
      END IF
!
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
        IEEECK = 0
        RETURN
      END IF
!
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
        IEEECK = 0
        RETURN
      END IF
!
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
        IEEECK = 0
        RETURN
      END IF
!
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
        IEEECK = 0
        RETURN
      END IF
!
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
        IEEECK = 0
        RETURN
      END IF
!
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
        IEEECK = 0
        RETURN
      END IF
!
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
        IEEECK = 0
        RETURN
      END IF
!
!
!
!
!     Return if we were only asked to check infinity arithmetic
!
      IF( ISPEC.EQ.0 )                                                  &
     &RETURN
!
      NAN1 = POSINF + NEGINF
!
      NAN2 = POSINF / NEGINF
!
      NAN3 = POSINF / POSINF
!
      NAN4 = POSINF*ZERO
!
      NAN5 = NEGINF*NEGZRO
!
      NAN6 = NAN5*0.0
!
      IF( NAN1.EQ.NAN1 ) THEN
        IEEECK = 0
        RETURN
      END IF
!
      IF( NAN2.EQ.NAN2 ) THEN
        IEEECK = 0
        RETURN
      END IF
!
      IF( NAN3.EQ.NAN3 ) THEN
        IEEECK = 0
        RETURN
      END IF
!
      IF( NAN4.EQ.NAN4 ) THEN
        IEEECK = 0
        RETURN
      END IF
!
      IF( NAN5.EQ.NAN5 ) THEN
        IEEECK = 0
        RETURN
      END IF
!
      IF( NAN6.EQ.NAN6 ) THEN
        IEEECK = 0
        RETURN
      END IF
!
      RETURN
      END
      INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, ISPEC, LWORK, N
      CHARACTER          NAME*( * ), OPTS*( * )
!
!  Purpose
!  =======
!
!       This program sets problem and machine dependent parameters
!       useful for xHSEQR and its subroutines. It is called whenever
!       ILAENV is called with 12 <= ISPEC <= 16
!
!  Arguments
!  =========
!
!       ISPEC  (input) integer scalar
!              ISPEC specifies which tunable parameter IPARMQ should
!              return.
!
!              ISPEC=12: (INMIN)  Matrices of order nmin or less
!                        are sent directly to xLAHQR, the implicit
!                        double shift QR algorithm.  NMIN must be
!                        at least 11.
!
!              ISPEC=13: (INWIN)  Size of the deflation window.
!                        This is best set greater than or equal to
!                        the number of simultaneous shifts NS.
!                        Larger matrices benefit from larger deflation
!                        windows.
!
!              ISPEC=14: (INIBL) Determines when to stop nibbling and
!                        invest in an (expensive) multi-shift QR sweep.
!                        If the aggressive early deflation subroutine
!                        finds LD converged eigenvalues from an order
!                        NW deflation window and LD.GT.(NW*NIBBLE)/100,
!                        then the next QR sweep is skipped and early
!                        deflation is applied immediately to the
!                        remaining active diagonal block.  Setting
!                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a
!                        multi-shift QR sweep whenever early deflation
!                        finds a converged eigenvalue.  Setting
!                        IPARMQ(ISPEC=14) greater than or equal to 100
!                        prevents TTQRE from skipping a multi-shift
!                        QR sweep.
!
!              ISPEC=15: (NSHFTS) The number of simultaneous shifts in
!                        a multi-shift QR iteration.
!
!              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the
!                        following meanings.
!                        0:  During the multi-shift QR sweep,
!                            xLAQR5 does not accumulate reflections and
!                            does not use matrix-matrix multiply to
!                            update the far-from-diagonal matrix
!                            entries.
!                        1:  During the multi-shift QR sweep,
!                            xLAQR5 and/or xLAQRaccumulates reflections and uses
!                            matrix-matrix multiply to update the
!                            far-from-diagonal matrix entries.
!                        2:  During the multi-shift QR sweep.
!                            xLAQR5 accumulates reflections and takes
!                            advantage of 2-by-2 block structure during
!                            matrix-matrix multiplies.
!                        (If xTRMM is slower than xGEMM, then
!                        IPARMQ(ISPEC=16)=1 may be more efficient than
!                        IPARMQ(ISPEC=16)=2 despite the greater level of
!                        arithmetic work implied by the latter choice.)
!
!       NAME    (input) character string
!               Name of the calling subroutine
!
!       OPTS    (input) character string
!               This is a concatenation of the string arguments to
!               TTQRE.
!
!       N       (input) integer scalar
!               N is the order of the Hessenberg matrix H.
!
!       ILO     (input) INTEGER
!       IHI     (input) INTEGER
!               It is assumed that H is already upper triangular
!               in rows and columns 1:ILO-1 and IHI+1:N.
!
!       LWORK   (input) integer scalar
!               The amount of workspace available.
!
!  Further Details
!  ===============
!
!       Little is known about how best to choose these parameters.
!       It is possible to use different values of the parameters
!       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.
!
!       It is probably best to choose different parameters for
!       different matrices and different parameters at different
!       times during the iteration, but this has not been
!       implemented --- yet.
!
!
!       The best choices of most of the parameters depend
!       in an ill-understood way on the relative execution
!       rate of xLAQR3 and xLAQR5 and on the nature of each
!       particular eigenvalue problem.  Experiment may be the
!       only practical way to determine which choices are most
!       effective.
!
!       Following is a list of default values supplied by IPARMQ.
!       These defaults may be adjusted in order to attain better
!       performance in any particular computational environment.
!
!       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.
!                        Default: 75. (Must be at least 11.)
!
!       IPARMQ(ISPEC=13) Recommended deflation window size.
!                        This depends on ILO, IHI and NS, the
!                        number of simultaneous shifts returned
!                        by IPARMQ(ISPEC=15).  The default for
!                        (IHI-ILO+1).LE.500 is NS.  The default
!                        for (IHI-ILO+1).GT.500 is 3*NS/2.
!
!       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.
!
!       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.
!                        a multi-shift QR iteration.
!
!                        If IHI-ILO+1 is ...
!
!                        greater than      ...but less    ... the
!                        or equal to ...      than        default is
!
!                                0               30       NS =   2+
!                               30               60       NS =   4+
!                               60              150       NS =  10
!                              150              590       NS =  **
!                              590             3000       NS =  64
!                             3000             6000       NS = 128
!                             6000             infinity   NS = 256
!
!                    (+)  By default matrices of this order are
!                         passed to the implicit double shift routine
!                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These
!                         values of NS are used only in case of a rare
!                         xLAHQR failure.
!
!                    (**) The asterisks (**) indicate an ad-hoc
!                         function increasing from 10 to 64.
!
!       IPARMQ(ISPEC=16) Select structured matrix multiply.
!                        (See ISPEC=16 above for details.)
!                        Default: 3.
!
!     ================================================================
!     .. Parameters ..
      INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22
      PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14,          &
     &ISHFTS = 15, IACC22 = 16 )
      INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
      PARAMETER          ( NMIN = 75, K22MIN = 14, KACMIN = 14,         &
     &NIBBLE = 14, KNWSWP = 500 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0 )
!     ..
!     .. Local Scalars ..
      INTEGER            NH, NS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          LOG, MAX, MOD, NINT, REAL
!     ..
!     .. Executable Statements ..
      IF( ( ISPEC.EQ.ISHFTS ) .OR. ( ISPEC.EQ.INWIN ) .OR.              &
     &( ISPEC.EQ.IACC22 ) ) THEN
!
!        ==== Set the number simultaneous shifts ====
!
      NH = IHI - ILO + 1
      NS = 2
      IF( NH.GE.30 )                                                    &
     &NS = 4
      IF( NH.GE.60 )                                                    &
     &NS = 10
      IF( NH.GE.150 )                                                   &
     &NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )
      IF( NH.GE.590 )                                                   &
     &NS = 64
      IF( NH.GE.3000 )                                                  &
     &NS = 128
      IF( NH.GE.6000 )                                                  &
     &NS = 256
      NS = MAX( 2, NS-MOD( NS, 2 ) )
      END IF
!
      IF( ISPEC.EQ.INMIN ) THEN
!
!
!        ===== Matrices of order smaller than NMIN get sent
!        .     to xLAHQR, the classic double shift algorithm.
!        .     This must be at least 11. ====
!
      IPARMQ = NMIN
!
      ELSE IF( ISPEC.EQ.INIBL ) THEN
!
!        ==== INIBL: skip a multi-shift qr iteration and
!        .    whenever aggressive early deflation finds
!        .    at least (NIBBLE*(window size)/100) deflations. ====
!
      IPARMQ = NIBBLE
!
      ELSE IF( ISPEC.EQ.ISHFTS ) THEN
!
!        ==== NSHFTS: The number of simultaneous shifts =====
!
      IPARMQ = NS
!
      ELSE IF( ISPEC.EQ.INWIN ) THEN
!
!        ==== NW: deflation window size.  ====
!
      IF( NH.LE.KNWSWP ) THEN
      IPARMQ = NS
      ELSE
      IPARMQ = 3*NS / 2
      END IF
!
      ELSE IF( ISPEC.EQ.IACC22 ) THEN
!
!        ==== IACC22: Whether to accumulate reflections
!        .     before updating the far-from-diagonal elements
!        .     and whether to use 2-by-2 block structure while
!        .     doing it.  A small amount of work could be saved
!        .     by making this choice dependent also upon the
!        .     NH=IHI-ILO+1.
!
      IPARMQ = 0
      IF( NS.GE.KACMIN )                                                &
     &IPARMQ = 1
      IF( NS.GE.K22MIN )                                                &
     &IPARMQ = 2
!
      ELSE
!        ===== invalid value of ispec =====
      IPARMQ = -1
!
      END IF
!
!     ==== End of IPARMQ ====
!
      END
      LOGICAL          FUNCTION LSAME( CA, CB )
      IMPLICIT NONE
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          CA, CB
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
!     ..
!     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
      LSAME = CA.EQ.CB
      IF( LSAME )                                                       &
     &RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR( 'Z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
!
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
      IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
      IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
!
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
      IF( INTA.GE.129 .AND. INTA.LE.137 .OR.                            &
     &INTA.GE.145 .AND. INTA.LE.153 .OR.                                &
     &INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
      IF( INTB.GE.129 .AND. INTB.LE.137 .OR.                            &
     &INTB.GE.145 .AND. INTB.LE.153 .OR.                                &
     &INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
      IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
      IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
!
!     RETURN
!
!     End of LSAME
!
      END
      SUBROUTINE DBDSDC( UPLO, COMPQ, N, D, E, U, LDU, VT, LDVT, Q, IQ, &
     &WORK, IWORK, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          COMPQ, UPLO
      INTEGER            INFO, LDU, LDVT, N
!     ..
!     .. Array Arguments ..
      INTEGER            IQ( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), Q( * ), U( LDU, * ),           &
     &VT( LDVT, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DBDSDC computes the singular value decomposition (SVD) of a real
!  N-by-N (upper or lower) bidiagonal matrix B:  B = U * S * VT,
!  using a divide and conquer method, where S is a diagonal matrix
!  with non-negative diagonal elements (the singular values of B), and
!  U and VT are orthogonal matrices of left and right singular vectors,
!  respectively. DBDSDC can be used to compute all singular values,
!  and optionally, singular vectors or singular vectors in compact form.
!
!  This code makes very mild assumptions about floating point
!  arithmetic. It will work on machines with a guard digit in
!  add/subtract, or on those binary machines without guard digits
!  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
!  It could conceivably fail on hexadecimal or decimal machines
!  without guard digits, but we know of none.  See DLASD3 for details.
!
!  The code currently calls DLASDQ if singular values only are desired.
!  However, it can be slightly modified to compute singular values
!  using the divide and conquer method.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  B is upper bidiagonal.
!          = 'L':  B is lower bidiagonal.
!
!  COMPQ   (input) CHARACTER*1
!          Specifies whether singular vectors are to be computed
!          as follows:
!          = 'N':  Compute singular values only;
!          = 'P':  Compute singular values and compute singular
!                  vectors in compact form;
!          = 'I':  Compute singular values and singular vectors.
!
!  N       (input) INTEGER
!          The order of the matrix B.  N >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the n diagonal elements of the bidiagonal matrix B.
!          On exit, if INFO=0, the singular values of B.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, the elements of E contain the offdiagonal
!          elements of the bidiagonal matrix whose SVD is desired.
!          On exit, E has been destroyed.
!
!  U       (output) DOUBLE PRECISION array, dimension (LDU,N)
!          If  COMPQ = 'I', then:
!             On exit, if INFO = 0, U contains the left singular vectors
!             of the bidiagonal matrix.
!          For other values of COMPQ, U is not referenced.
!
!  LDU     (input) INTEGER
!          The leading dimension of the array U.  LDU >= 1.
!          If singular vectors are desired, then LDU >= max( 1, N ).
!
!  VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)
!          If  COMPQ = 'I', then:
!             On exit, if INFO = 0, VT' contains the right singular
!             vectors of the bidiagonal matrix.
!          For other values of COMPQ, VT is not referenced.
!
!  LDVT    (input) INTEGER
!          The leading dimension of the array VT.  LDVT >= 1.
!          If singular vectors are desired, then LDVT >= max( 1, N ).
!
!  Q       (output) DOUBLE PRECISION array, dimension (LDQ)
!          If  COMPQ = 'P', then:
!             On exit, if INFO = 0, Q and IQ contain the left
!             and right singular vectors in a compact form,
!             requiring O(N log N) space instead of 2*N**2.
!             In particular, Q contains all the DOUBLE PRECISION data in
!             LDQ >= N*(11 + 2*SMLSIZ + 8*INT(LOG_2(N/(SMLSIZ+1))))
!             words of memory, where SMLSIZ is returned by ILAENV and
!             is equal to the maximum size of the subproblems at the
!             bottom of the computation tree (usually about 25).
!          For other values of COMPQ, Q is not referenced.
!
!  IQ      (output) INTEGER array, dimension (LDIQ)
!          If  COMPQ = 'P', then:
!             On exit, if INFO = 0, Q and IQ contain the left
!             and right singular vectors in a compact form,
!             requiring O(N log N) space instead of 2*N**2.
!             In particular, IQ contains all INTEGER data in
!             LDIQ >= N*(3 + 3*INT(LOG_2(N/(SMLSIZ+1))))
!             words of memory, where SMLSIZ is returned by ILAENV and
!             is equal to the maximum size of the subproblems at the
!             bottom of the computation tree (usually about 25).
!          For other values of COMPQ, IQ is not referenced.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          If COMPQ = 'N' then LWORK >= (4 * N).
!          If COMPQ = 'P' then LWORK >= (6 * N).
!          If COMPQ = 'I' then LWORK >= (3 * N**2 + 4 * N).
!
!  IWORK   (workspace) INTEGER array, dimension (8*N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  The algorithm failed to compute an singular value.
!                The update process of divide and conquer failed.
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ming Gu and Huan Ren, Computer Science Division, University of
!     California at Berkeley, USA
!
!  =====================================================================
!  Changed dimension statement in comment describing E from (N) to
!  (N-1).  Sven, 17 Feb 05.
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            DIFL, DIFR, GIVCOL, GIVNUM, GIVPTR, I, IC,     &
     &ICOMPQ, IERR, II, IS, IU, IUPLO, IVT, J, K, KK,                   &
     &MLVL, NM1, NSIZE, PERM, POLES, QSTART, SMLSIZ,                    &
     &SMLSZP, SQRE, START, WSTART, Z
      DOUBLE PRECISION   CS, EPS, ORGNRM, P, R, SN
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANST
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANST
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPL, DLARTG, DLASCL, DLASD0, DLASDA, DLASDQ, &
     &DLASET, DLASR, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, LOG, SIGN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IUPLO = 0
      IF( LSAME( UPLO, 'U' ) )                                          &
     &IUPLO = 1
      IF( LSAME( UPLO, 'L' ) )                                          &
     &IUPLO = 2
      IF( LSAME( COMPQ, 'N' ) ) THEN
        ICOMPQ = 0
      ELSE IF( LSAME( COMPQ, 'P' ) ) THEN
        ICOMPQ = 1
      ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
        ICOMPQ = 2
      ELSE
        ICOMPQ = -1
      END IF
      IF( IUPLO.EQ.0 ) THEN
        INFO = -1
      ELSE IF( ICOMPQ.LT.0 ) THEN
        INFO = -2
      ELSE IF( N.LT.0 ) THEN
        INFO = -3
      ELSE IF( ( LDU.LT.1 ) .OR. ( ( ICOMPQ.EQ.2 ) .AND. ( LDU.LT.      &
     &N ) ) ) THEN
        INFO = -7
      ELSE IF( ( LDVT.LT.1 ) .OR. ( ( ICOMPQ.EQ.2 ) .AND. ( LDVT.LT.    &
     &N ) ) ) THEN
        INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DBDSDC', -INFO )
        RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 )                                                      &
     &RETURN
      SMLSIZ = ILAENV( 9, 'DBDSDC', ' ', 0, 0, 0, 0 )
      IF( N.EQ.1 ) THEN
        IF( ICOMPQ.EQ.1 ) THEN
          Q( 1 ) = SIGN( ONE, D( 1 ) )
          Q( 1+SMLSIZ*N ) = ONE
        ELSE IF( ICOMPQ.EQ.2 ) THEN
          U( 1, 1 ) = SIGN( ONE, D( 1 ) )
          VT( 1, 1 ) = ONE
        END IF
        D( 1 ) = ABS( D( 1 ) )
        RETURN
      END IF
      NM1 = N - 1
!
!     If matrix lower bidiagonal, rotate to be upper bidiagonal
!     by applying Givens rotations on the left
!
      WSTART = 1
      QSTART = 3
      IF( ICOMPQ.EQ.1 ) THEN
        CALL DCOPL( N, D, 1, Q( 1 ), 1 )
        CALL DCOPL( N-1, E, 1, Q( N+1 ), 1 )
      END IF
      IF( IUPLO.EQ.2 ) THEN
        QSTART = 5
        WSTART = 2*N - 1
        DO I = 1, N - 1
          CALL DLARTG( D( I ), E( I ), CS, SN, R )
          D( I ) = R
          E( I ) = SN*D( I+1 )
          D( I+1 ) = CS*D( I+1 )
          IF( ICOMPQ.EQ.1 ) THEN
            Q( I+2*N ) = CS
            Q( I+3*N ) = SN
          ELSE IF( ICOMPQ.EQ.2 ) THEN
            WORK( I ) = CS
            WORK( NM1+I ) = -SN
          END IF
        enddo
      END IF
!
!     If ICOMPQ = 0, use DLASDQ to compute the singular values.
!
      IF( ICOMPQ.EQ.0 ) THEN
        CALL DLASDQ( 'U', 0, N, 0, 0, 0, D, E, VT, LDVT, U, LDU, U,     &
     &LDU, WORK( WSTART ), INFO )
        GO TO 40
      END IF
!
!     If N is smaller than the minimum divide size SMLSIZ, then solve
!     the problem with another solver.
!
      IF( N.LE.SMLSIZ ) THEN
        IF( ICOMPQ.EQ.2 ) THEN
          CALL DLASET( 'A', N, N, ZERO, ONE, U, LDU )
          CALL DLASET( 'A', N, N, ZERO, ONE, VT, LDVT )
          CALL DLASDQ( 'U', 0, N, N, N, 0, D, E, VT, LDVT, U, LDU, U,   &
     &LDU, WORK( WSTART ), INFO )
        ELSE IF( ICOMPQ.EQ.1 ) THEN
          IU = 1
          IVT = IU + N
          CALL DLASET( 'A', N, N, ZERO, ONE, Q( IU+( QSTART-1 )*N ),    &
     &N )
          CALL DLASET( 'A', N, N, ZERO, ONE, Q( IVT+( QSTART-1 )*N ),   &
     &N )
          CALL DLASDQ( 'U', 0, N, N, N, 0, D, E,                        &
     &Q( IVT+( QSTART-1 )*N ), N,                                       &
     &Q( IU+( QSTART-1 )*N ), N,                                        &
     &Q( IU+( QSTART-1 )*N ), N, WORK( WSTART ),                        &
     &INFO )
        END IF
        GO TO 40
      END IF
!
      IF( ICOMPQ.EQ.2 ) THEN
        CALL DLASET( 'A', N, N, ZERO, ONE, U, LDU )
        CALL DLASET( 'A', N, N, ZERO, ONE, VT, LDVT )
      END IF
!
!     Scale.
!
      ORGNRM = DLANST( 'M', N, D, E )
      IF( ORGNRM.EQ.ZERO )                                              &
     &RETURN
      CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, N, 1, D, N, IERR )
      CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, NM1, 1, E, NM1, IERR )
!
      EPS = DLAMCH( 'Epsilon' )
!
      MLVL = INT( LOG( DBLE( N ) / DBLE( SMLSIZ+1 ) ) / LOG( TWO ) ) + 1
      SMLSZP = SMLSIZ + 1
!
      IF( ICOMPQ.EQ.1 ) THEN
        IU = 1
        IVT = 1 + SMLSIZ
        DIFL = IVT + SMLSZP
        DIFR = DIFL + MLVL
        Z = DIFR + MLVL*2
        IC = Z + MLVL
        IS = IC + 1
        POLES = IS + 1
        GIVNUM = POLES + 2*MLVL
!
        K = 1
        GIVPTR = 2
        PERM = 3
        GIVCOL = PERM + MLVL
      END IF
!
      DO I = 1, N
        IF( ABS( D( I ) ).LT.EPS ) THEN
          D( I ) = SIGN( EPS, D( I ) )
        END IF
      enddo
!
      START = 1
      SQRE = 0
!
      DO I = 1, NM1
        IF( ( ABS( E( I ) ).LT.EPS ) .OR. ( I.EQ.NM1 ) ) THEN
!
!        Subproblem found. First determine its size and then
!        apply divide and conquer on it.
!
          IF( I.LT.NM1 ) THEN
!
!        A subproblem with E(I) small for I < NM1.
!
            NSIZE = I - START + 1
          ELSE IF( ABS( E( I ) ).GE.EPS ) THEN
!
!        A subproblem with E(NM1) not too small but I = NM1.
!
            NSIZE = N - START + 1
          ELSE
!
!        A subproblem with E(NM1) small. This implies an
!        1-by-1 subproblem at D(N). Solve this 1-by-1 problem
!        first.
!
            NSIZE = I - START + 1
            IF( ICOMPQ.EQ.2 ) THEN
              U( N, N ) = SIGN( ONE, D( N ) )
              VT( N, N ) = ONE
            ELSE IF( ICOMPQ.EQ.1 ) THEN
              Q( N+( QSTART-1 )*N ) = SIGN( ONE, D( N ) )
              Q( N+( SMLSIZ+QSTART-1 )*N ) = ONE
            END IF
            D( N ) = ABS( D( N ) )
          END IF
          IF( ICOMPQ.EQ.2 ) THEN
            CALL DLASD0( NSIZE, SQRE, D( START ), E( START ),           &
     &U( START, START ), LDU, VT( START, START ),                       &
     &LDVT, SMLSIZ, IWORK, WORK( WSTART ), INFO )
          ELSE
            CALL DLASDA( ICOMPQ, SMLSIZ, NSIZE, SQRE, D( START ),       &
     &E( START ), Q( START+( IU+QSTART-2 )*N ), N,                      &
     &Q( START+( IVT+QSTART-2 )*N ),                                    &
     &IQ( START+K*N ), Q( START+( DIFL+QSTART-2 )*                      &
     &N ), Q( START+( DIFR+QSTART-2 )*N ),                              &
     &Q( START+( Z+QSTART-2 )*N ),                                      &
     &Q( START+( POLES+QSTART-2 )*N ),                                  &
     &IQ( START+GIVPTR*N ), IQ( START+GIVCOL*N ),                       &
     &N, IQ( START+PERM*N ),                                            &
     &Q( START+( GIVNUM+QSTART-2 )*N ),                                 &
     &Q( START+( IC+QSTART-2 )*N ),                                     &
     &Q( START+( IS+QSTART-2 )*N ),                                     &
     &WORK( WSTART ), IWORK, INFO )
            IF( INFO.NE.0 ) THEN
              RETURN
            END IF
          END IF
          START = I + 1
        END IF
      enddo
!
!     Unscale
!
      CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, N, 1, D, N, IERR )
   40 CONTINUE
!
!     Use Selection Sort to minimize swaps of singular vectors
!
      DO II = 2, N
        I = II - 1
        KK = I
        P = D( I )
        DO J = II, N
          IF( D( J ).GT.P ) THEN
            KK = J
            P = D( J )
          END IF
        enddo
        IF( KK.NE.I ) THEN
          D( KK ) = D( I )
          D( I ) = P
          IF( ICOMPQ.EQ.1 ) THEN
            IQ( I ) = KK
          ELSE IF( ICOMPQ.EQ.2 ) THEN
            CALL DSWAP( N, U( 1, I ), 1, U( 1, KK ), 1 )
            CALL DSWAP( N, VT( I, 1 ), LDVT, VT( KK, 1 ), LDVT )
          END IF
        ELSE IF( ICOMPQ.EQ.1 ) THEN
          IQ( I ) = I
        END IF
      enddo
!
!     If ICOMPQ = 1, use IQ(N,1) as the indicator for UPLO
!
      IF( ICOMPQ.EQ.1 ) THEN
        IF( IUPLO.EQ.1 ) THEN
          IQ( N ) = 1
        ELSE
          IQ( N ) = 0
        END IF
      END IF
!
!     If B is lower bidiagonal, update U by those Givens rotations
!     which rotated B to be upper bidiagonal
!
      IF( ( IUPLO.EQ.2 ) .AND. ( ICOMPQ.EQ.2 ) )                        &
     &CALL DLASR( 'L', 'V', 'B', N, N, WORK( 1 ), WORK( N ), U, LDU )
!
      RETURN
!
!     End of DBDSDC
!
      END
      SUBROUTINE DBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U,    &
     &LDU, C, LDC, WORK, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), D( * ), E( * ), U( LDU, * ),      &
     &VT( LDVT, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DBDSQR computes the singular values and, optionally, the right and/or
!  left singular vectors from the singular value decomposition (SVD) of
!  a real N-by-N (upper or lower) bidiagonal matrix B using the implicit
!  zero-shift QR algorithm.  The SVD of B has the form
!
!     B = Q * S * P**T
!
!  where S is the diagonal matrix of singular values, Q is an orthogonal
!  matrix of left singular vectors, and P is an orthogonal matrix of
!  right singular vectors.  If left singular vectors are requested, this
!  subroutine actually returns U*Q instead of Q, and, if right singular
!  vectors are requested, this subroutine returns P**T*VT instead of
!  P**T, for given real input matrices U and VT.  When U and VT are the
!  orthogonal matrices that reduce a general matrix A to bidiagonal
!  form:  A = U*B*VT, as computed by DGEBRD, then
!
!     A = (U*Q) * S * (P**T*VT)
!
!  is the SVD of A.  Optionally, the subroutine may also compute Q**T*C
!  for a given real input matrix C.
!
!  See "Computing  Small Singular Values of Bidiagonal Matrices With
!  Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
!  LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,
!  no. 5, pp. 873-912, Sept 1990) and
!  "Accurate singular values and differential qd algorithms," by
!  B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics
!  Department, University of California at Berkeley, July 1992
!  for a detailed description of the algorithm.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  B is upper bidiagonal;
!          = 'L':  B is lower bidiagonal.
!
!  N       (input) INTEGER
!          The order of the matrix B.  N >= 0.
!
!  NCVT    (input) INTEGER
!          The number of columns of the matrix VT. NCVT >= 0.
!
!  NRU     (input) INTEGER
!          The number of rows of the matrix U. NRU >= 0.
!
!  NCC     (input) INTEGER
!          The number of columns of the matrix C. NCC >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the n diagonal elements of the bidiagonal matrix B.
!          On exit, if INFO=0, the singular values of B in decreasing
!          order.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, the N-1 offdiagonal elements of the bidiagonal
!          matrix B.
!          On exit, if INFO = 0, E is destroyed; if INFO > 0, D and E
!          will contain the diagonal and superdiagonal elements of a
!          bidiagonal matrix orthogonally equivalent to the one given
!          as input.
!
!  VT      (input/output) DOUBLE PRECISION array, dimension (LDVT, NCVT)
!          On entry, an N-by-NCVT matrix VT.
!          On exit, VT is overwritten by P**T * VT.
!          Not referenced if NCVT = 0.
!
!  LDVT    (input) INTEGER
!          The leading dimension of the array VT.
!          LDVT >= max(1,N) if NCVT > 0; LDVT >= 1 if NCVT = 0.
!
!  U       (input/output) DOUBLE PRECISION array, dimension (LDU, N)
!          On entry, an NRU-by-N matrix U.
!          On exit, U is overwritten by U * Q.
!          Not referenced if NRU = 0.
!
!  LDU     (input) INTEGER
!          The leading dimension of the array U.  LDU >= max(1,NRU).
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC, NCC)
!          On entry, an N-by-NCC matrix C.
!          On exit, C is overwritten by Q**T * C.
!          Not referenced if NCC = 0.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C.
!          LDC >= max(1,N) if NCC > 0; LDC >=1 if NCC = 0.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)
!          if NCVT = NRU = NCC = 0, (max(1, 4*N-4)) otherwise
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  If INFO = -i, the i-th argument had an illegal value
!          > 0:  the algorithm did not converge; D and E contain the
!                elements of a bidiagonal matrix which is orthogonally
!                similar to the input matrix B;  if INFO = i, i
!                elements of E have not converged to zero.
!
!  Internal Parameters
!  ===================
!
!  TOLMUL  DOUBLE PRECISION, default = max(10,min(100,EPS**(-1/8)))
!          TOLMUL controls the convergence criterion of the QR loop.
!          If it is positive, TOLMUL*EPS is the desired relative
!             precision in the computed singular values.
!          If it is negative, abs(TOLMUL*EPS*sigma_max) is the
!             desired absolute accuracy in the computed singular
!             values (corresponds to relative accuracy
!             abs(TOLMUL*EPS) in the largest singular value.
!          abs(TOLMUL) should be between 1 and 1/EPS, and preferably
!             between 10 (for fast convergence) and .1/EPS
!             (for there to be some accuracy in the results).
!          Default is to lose at either one eighth or 2 of the
!             available decimal digits in each computed singular value
!             (whichever is smaller).
!
!  MAXITR  INTEGER, default = 6
!          MAXITR controls the maximum number of passes of the
!          algorithm through its inner loop. The algorithms stops
!          (and so fails to converge) if the number of passes
!          through the inner loop exceeds MAXITR*N**2.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   NEGONE
      PARAMETER          ( NEGONE = -1.0D0 )
      DOUBLE PRECISION   HNDRTH
      PARAMETER          ( HNDRTH = 0.01D0 )
      DOUBLE PRECISION   TEN
      PARAMETER          ( TEN = 10.0D0 )
      DOUBLE PRECISION   HNDRD
      PARAMETER          ( HNDRD = 100.0D0 )
      DOUBLE PRECISION   MEIGTH
      PARAMETER          ( MEIGTH = -0.125D0 )
      INTEGER            MAXITR
      PARAMETER          ( MAXITR = 6 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LOWER, ROTATE
      INTEGER            I, IDIR, ISUB, ITER, J, LL, LLL, M, MAXIT, NM1,&
     &NM12, NM13, OLDLL, OLDM
      DOUBLE PRECISION   ABSE, ABSS, COSL, COSR, CS, EPS, F, G, H, MU,  &
     &OLDCS, OLDSN, R, SHIFT, SIGMN, SIGMX, SINL,                       &
     &SINR, SLL, SMAX, SMIN, SMINL, SMINOA,                             &
     &SN, THRESH, TOL, TOLMUL, UNFL
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARTG, DLAS2, DLASQ1, DLASR, DLASV2, DROT,    &
     &DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      LOWER = LSAME( UPLO, 'L' )
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LOWER ) THEN
        INFO = -1
      ELSE IF( N.LT.0 ) THEN
        INFO = -2
      ELSE IF( NCVT.LT.0 ) THEN
        INFO = -3
      ELSE IF( NRU.LT.0 ) THEN
        INFO = -4
      ELSE IF( NCC.LT.0 ) THEN
        INFO = -5
      ELSE IF( ( NCVT.EQ.0 .AND. LDVT.LT.1 ) .OR.                       &
     &( NCVT.GT.0 .AND. LDVT.LT.MAX( 1, N ) ) ) THEN
        INFO = -9
      ELSE IF( LDU.LT.MAX( 1, NRU ) ) THEN
        INFO = -11
      ELSE IF( ( NCC.EQ.0 .AND. LDC.LT.1 ) .OR.                         &
     &( NCC.GT.0 .AND. LDC.LT.MAX( 1, N ) ) ) THEN
        INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DBDSQR', -INFO )
        RETURN
      END IF
      IF( N.EQ.0 )                                                      &
     &RETURN
      IF( N.EQ.1 )                                                      &
     &GO TO 160
!
!     ROTATE is true if any singular vectors desired, false otherwise
!
      ROTATE = ( NCVT.GT.0 ) .OR. ( NRU.GT.0 ) .OR. ( NCC.GT.0 )
!
!     If no singular vectors desired, use qd algorithm
!
      IF( .NOT.ROTATE ) THEN
        CALL DLASQ1( N, D, E, WORK, INFO )
        RETURN
      END IF
!
      NM1 = N - 1
      NM12 = NM1 + NM1
      NM13 = NM12 + NM1
      IDIR = 0
!
!     Get machine constants
!
      EPS = DLAMCH( 'Epsilon' )
      UNFL = DLAMCH( 'Safe minimum' )
!
!     If matrix lower bidiagonal, rotate to be upper bidiagonal
!     by applying Givens rotations on the left
!
      IF( LOWER ) THEN
        DO I = 1, N - 1
          CALL DLARTG( D( I ), E( I ), CS, SN, R )
          D( I ) = R
          E( I ) = SN*D( I+1 )
          D( I+1 ) = CS*D( I+1 )
          WORK( I ) = CS
          WORK( NM1+I ) = SN
        enddo
!
!        Update singular vectors if desired
!
        IF( NRU.GT.0 )                                                  &
     &CALL DLASR( 'R', 'V', 'F', NRU, N, WORK( 1 ), WORK( N ), U,       &
     &LDU )
        IF( NCC.GT.0 )                                                  &
     &CALL DLASR( 'L', 'V', 'F', N, NCC, WORK( 1 ), WORK( N ), C,       &
     &LDC )
      END IF
!
!     Compute singular values to relative accuracy TOL
!     (By setting TOL to be negative, algorithm will compute
!     singular values to absolute accuracy ABS(TOL)*norm(input matrix))
!
      TOLMUL = MAX( TEN, MIN( HNDRD, EPS**MEIGTH ) )
      TOL = TOLMUL*EPS
!
!     Compute approximate maximum, minimum singular values
!
      SMAX = ZERO
      DO I = 1, N
        SMAX = MAX( SMAX, ABS( D( I ) ) )
      enddo
      DO I = 1, N - 1
        SMAX = MAX( SMAX, ABS( E( I ) ) )
      enddo
      SMINL = ZERO
      IF( TOL.GE.ZERO ) THEN
!
!        Relative accuracy desired
!
        SMINOA = ABS( D( 1 ) )
        IF( SMINOA.EQ.ZERO )                                            &
     &GO TO 50
        MU = SMINOA
        DO I = 2, N
          MU = ABS( D( I ) )*( MU / ( MU+ABS( E( I-1 ) ) ) )
          SMINOA = MIN( SMINOA, MU )
          IF( SMINOA.EQ.ZERO )                                          &
     &GO TO 50
        enddo
   50   CONTINUE
        SMINOA = SMINOA / SQRT( DBLE( N ) )
        THRESH = MAX( TOL*SMINOA, MAXITR*N*N*UNFL )
      ELSE
!
!        Absolute accuracy desired
!
        THRESH = MAX( ABS( TOL )*SMAX, MAXITR*N*N*UNFL )
      END IF
!
!     Prepare for main iteration loop for the singular values
!     (MAXIT is the maximum number of passes through the inner
!     loop permitted before nonconvergence signalled.)
!
      MAXIT = MAXITR*N*N
      ITER = 0
      OLDLL = -1
      OLDM = -1
!
!     M points to last element of unconverged part of matrix
!
      M = N
!
!     Begin main iteration loop
!
   60 CONTINUE
!
!     Check for convergence or exceeding iteration count
!
      IF( M.LE.1 )                                                      &
     &GO TO 160
      IF( ITER.GT.MAXIT )                                               &
     &GO TO 200
!
!     Find diagonal block of matrix to work on
!
      IF( TOL.LT.ZERO .AND. ABS( D( M ) ).LE.THRESH )                   &
     &D( M ) = ZERO
      SMAX = ABS( D( M ) )
      SMIN = SMAX
      DO LLL = 1, M - 1
        LL = M - LLL
        ABSS = ABS( D( LL ) )
        ABSE = ABS( E( LL ) )
        IF( TOL.LT.ZERO .AND. ABSS.LE.THRESH )                          &
     &D( LL ) = ZERO
        IF( ABSE.LE.THRESH )                                            &
     &GO TO 80
        SMIN = MIN( SMIN, ABSS )
        SMAX = MAX( SMAX, ABSS, ABSE )
      enddo
      LL = 0
      GO TO 90
   80 CONTINUE
      E( LL ) = ZERO
!
!     Matrix splits since E(LL) = 0
!
      IF( LL.EQ.M-1 ) THEN
!
!        Convergence of bottom singular value, return to top of loop
!
        M = M - 1
        GO TO 60
      END IF
   90 CONTINUE
      LL = LL + 1
!
!     E(LL) through E(M-1) are nonzero, E(LL-1) is zero
!
      IF( LL.EQ.M-1 ) THEN
!
!        2 by 2 block, handle separately
!
        CALL DLASV2( D( M-1 ), E( M-1 ), D( M ), SIGMN, SIGMX, SINR,    &
     &COSR, SINL, COSL )
        D( M-1 ) = SIGMX
        E( M-1 ) = ZERO
        D( M ) = SIGMN
!
!        Compute singular vectors, if desired
!
        IF( NCVT.GT.0 )                                                 &
     &CALL DROT( NCVT, VT( M-1, 1 ), LDVT, VT( M, 1 ), LDVT, COSR,      &
     &SINR )
        IF( NRU.GT.0 )                                                  &
     &CALL DROT( NRU, U( 1, M-1 ), 1, U( 1, M ), 1, COSL, SINL )
        IF( NCC.GT.0 )                                                  &
     &CALL DROT( NCC, C( M-1, 1 ), LDC, C( M, 1 ), LDC, COSL,           &
     &SINL )
        M = M - 2
        GO TO 60
      END IF
!
!     If working on new submatrix, choose shift direction
!     (from larger end diagonal element towards smaller)
!
      IF( LL.GT.OLDM .OR. M.LT.OLDLL ) THEN
        IF( ABS( D( LL ) ).GE.ABS( D( M ) ) ) THEN
!
!           Chase bulge from top (big end) to bottom (small end)
!
          IDIR = 1
        ELSE
!
!           Chase bulge from bottom (big end) to top (small end)
!
          IDIR = 2
        END IF
      END IF
!
!     Apply convergence tests
!
      IF( IDIR.EQ.1 ) THEN
!
!        Run convergence test in forward direction
!        First apply standard test to bottom of matrix
!
        IF( ABS( E( M-1 ) ).LE.ABS( TOL )*ABS( D( M ) ) .OR.            &
     &( TOL.LT.ZERO .AND. ABS( E( M-1 ) ).LE.THRESH ) ) THEN
          E( M-1 ) = ZERO
          GO TO 60
        END IF
!
        IF( TOL.GE.ZERO ) THEN
!
!           If relative accuracy desired,
!           apply convergence criterion forward
!
          MU = ABS( D( LL ) )
          SMINL = MU
          DO LLL = LL, M - 1
            IF( ABS( E( LLL ) ).LE.TOL*MU ) THEN
              E( LLL ) = ZERO
              GO TO 60
            END IF
            MU = ABS( D( LLL+1 ) )*( MU / ( MU+ABS( E( LLL ) ) ) )
            SMINL = MIN( SMINL, MU )
          enddo
        END IF
!
      ELSE
!
!        Run convergence test in backward direction
!        First apply standard test to top of matrix
!
        IF( ABS( E( LL ) ).LE.ABS( TOL )*ABS( D( LL ) ) .OR.            &
     &( TOL.LT.ZERO .AND. ABS( E( LL ) ).LE.THRESH ) ) THEN
          E( LL ) = ZERO
          GO TO 60
        END IF
!
        IF( TOL.GE.ZERO ) THEN
!
!           If relative accuracy desired,
!           apply convergence criterion backward
!
          MU = ABS( D( M ) )
          SMINL = MU
          DO LLL = M - 1, LL, -1
            IF( ABS( E( LLL ) ).LE.TOL*MU ) THEN
              E( LLL ) = ZERO
              GO TO 60
            END IF
            MU = ABS( D( LLL ) )*( MU / ( MU+ABS( E( LLL ) ) ) )
            SMINL = MIN( SMINL, MU )
          enddo
        END IF
      END IF
      OLDLL = LL
      OLDM = M
!
!     Compute shift.  First, test if shifting would ruin relative
!     accuracy, and if so set the shift to zero.
!
      IF( TOL.GE.ZERO .AND. N*TOL*( SMINL / SMAX ).LE.                  &
     &MAX( EPS, HNDRTH*TOL ) ) THEN
!
!        Use a zero shift to avoid loss of relative accuracy
!
        SHIFT = ZERO
      ELSE
!
!        Compute the shift from 2-by-2 block at end of matrix
!
        IF( IDIR.EQ.1 ) THEN
          SLL = ABS( D( LL ) )
          CALL DLAS2( D( M-1 ), E( M-1 ), D( M ), SHIFT, R )
        ELSE
          SLL = ABS( D( M ) )
          CALL DLAS2( D( LL ), E( LL ), D( LL+1 ), SHIFT, R )
        END IF
!
!        Test if shift negligible, and if so set to zero
!
        IF( SLL.GT.ZERO ) THEN
          IF( ( SHIFT / SLL )**2.LT.EPS )                               &
     &SHIFT = ZERO
        END IF
      END IF
!
!     Increment iteration count
!
      ITER = ITER + M - LL
!
!     If SHIFT = 0, do simplified QR iteration
!
      IF( SHIFT.EQ.ZERO ) THEN
        IF( IDIR.EQ.1 ) THEN
!
!           Chase bulge from top to bottom
!           Save cosines and sines for later singular vector updates
!
          CS = ONE
          OLDCS = ONE
          DO I = LL, M - 1
            CALL DLARTG( D( I )*CS, E( I ), CS, SN, R )
            IF( I.GT.LL )                                               &
     &E( I-1 ) = OLDSN*R
            CALL DLARTG( OLDCS*R, D( I+1 )*SN, OLDCS, OLDSN, D( I ) )
            WORK( I-LL+1 ) = CS
            WORK( I-LL+1+NM1 ) = SN
            WORK( I-LL+1+NM12 ) = OLDCS
            WORK( I-LL+1+NM13 ) = OLDSN
          enddo
          H = D( M )*CS
          D( M ) = H*OLDCS
          E( M-1 ) = H*OLDSN
!
!           Update singular vectors
!
          IF( NCVT.GT.0 )                                               &
     &CALL DLASR( 'L', 'V', 'F', M-LL+1, NCVT, WORK( 1 ),               &
     &WORK( N ), VT( LL, 1 ), LDVT )
          IF( NRU.GT.0 )                                                &
     &CALL DLASR( 'R', 'V', 'F', NRU, M-LL+1, WORK( NM12+1 ),           &
     &WORK( NM13+1 ), U( 1, LL ), LDU )
          IF( NCC.GT.0 )                                                &
     &CALL DLASR( 'L', 'V', 'F', M-LL+1, NCC, WORK( NM12+1 ),           &
     &WORK( NM13+1 ), C( LL, 1 ), LDC )
!
!           Test convergence
!
          IF( ABS( E( M-1 ) ).LE.THRESH )                               &
     &E( M-1 ) = ZERO
!
        ELSE
!
!           Chase bulge from bottom to top
!           Save cosines and sines for later singular vector updates
!
          CS = ONE
          OLDCS = ONE
          DO I = M, LL + 1, -1
            CALL DLARTG( D( I )*CS, E( I-1 ), CS, SN, R )
            IF( I.LT.M )                                                &
     &E( I ) = OLDSN*R
            CALL DLARTG( OLDCS*R, D( I-1 )*SN, OLDCS, OLDSN, D( I ) )
            WORK( I-LL ) = CS
            WORK( I-LL+NM1 ) = -SN
            WORK( I-LL+NM12 ) = OLDCS
            WORK( I-LL+NM13 ) = -OLDSN
          enddo
          H = D( LL )*CS
          D( LL ) = H*OLDCS
          E( LL ) = H*OLDSN
!
!           Update singular vectors
!
          IF( NCVT.GT.0 )                                               &
     &CALL DLASR( 'L', 'V', 'B', M-LL+1, NCVT, WORK( NM12+1 ),          &
     &WORK( NM13+1 ), VT( LL, 1 ), LDVT )
          IF( NRU.GT.0 )                                                &
     &CALL DLASR( 'R', 'V', 'B', NRU, M-LL+1, WORK( 1 ),                &
     &WORK( N ), U( 1, LL ), LDU )
          IF( NCC.GT.0 )                                                &
     &CALL DLASR( 'L', 'V', 'B', M-LL+1, NCC, WORK( 1 ),                &
     &WORK( N ), C( LL, 1 ), LDC )
!
!           Test convergence
!
          IF( ABS( E( LL ) ).LE.THRESH )                                &
     &E( LL ) = ZERO
        END IF
      ELSE
!
!        Use nonzero shift
!
        IF( IDIR.EQ.1 ) THEN
!
!           Chase bulge from top to bottom
!           Save cosines and sines for later singular vector updates
!
          F = ( ABS( D( LL ) )-SHIFT )*                                 &
     &( SIGN( ONE, D( LL ) )+SHIFT / D( LL ) )
          G = E( LL )
          DO I = LL, M - 1
            CALL DLARTG( F, G, COSR, SINR, R )
            IF( I.GT.LL )                                               &
     &E( I-1 ) = R
            F = COSR*D( I ) + SINR*E( I )
            E( I ) = COSR*E( I ) - SINR*D( I )
            G = SINR*D( I+1 )
            D( I+1 ) = COSR*D( I+1 )
            CALL DLARTG( F, G, COSL, SINL, R )
            D( I ) = R
            F = COSL*E( I ) + SINL*D( I+1 )
            D( I+1 ) = COSL*D( I+1 ) - SINL*E( I )
            IF( I.LT.M-1 ) THEN
              G = SINL*E( I+1 )
              E( I+1 ) = COSL*E( I+1 )
            END IF
            WORK( I-LL+1 ) = COSR
            WORK( I-LL+1+NM1 ) = SINR
            WORK( I-LL+1+NM12 ) = COSL
            WORK( I-LL+1+NM13 ) = SINL
          enddo
          E( M-1 ) = F
!
!           Update singular vectors
!
          IF( NCVT.GT.0 )                                               &
     &CALL DLASR( 'L', 'V', 'F', M-LL+1, NCVT, WORK( 1 ),               &
     &WORK( N ), VT( LL, 1 ), LDVT )
          IF( NRU.GT.0 )                                                &
     &CALL DLASR( 'R', 'V', 'F', NRU, M-LL+1, WORK( NM12+1 ),           &
     &WORK( NM13+1 ), U( 1, LL ), LDU )
          IF( NCC.GT.0 )                                                &
     &CALL DLASR( 'L', 'V', 'F', M-LL+1, NCC, WORK( NM12+1 ),           &
     &WORK( NM13+1 ), C( LL, 1 ), LDC )
!
!           Test convergence
!
          IF( ABS( E( M-1 ) ).LE.THRESH )                               &
     &E( M-1 ) = ZERO
!
        ELSE
!
!           Chase bulge from bottom to top
!           Save cosines and sines for later singular vector updates
!
          F = ( ABS( D( M ) )-SHIFT )*( SIGN( ONE, D( M ) )+SHIFT /     &
     &D( M ) )
          G = E( M-1 )
          DO I = M, LL + 1, -1
            CALL DLARTG( F, G, COSR, SINR, R )
            IF( I.LT.M )                                                &
     &E( I ) = R
            F = COSR*D( I ) + SINR*E( I-1 )
            E( I-1 ) = COSR*E( I-1 ) - SINR*D( I )
            G = SINR*D( I-1 )
            D( I-1 ) = COSR*D( I-1 )
            CALL DLARTG( F, G, COSL, SINL, R )
            D( I ) = R
            F = COSL*E( I-1 ) + SINL*D( I-1 )
            D( I-1 ) = COSL*D( I-1 ) - SINL*E( I-1 )
            IF( I.GT.LL+1 ) THEN
              G = SINL*E( I-2 )
              E( I-2 ) = COSL*E( I-2 )
            END IF
            WORK( I-LL ) = COSR
            WORK( I-LL+NM1 ) = -SINR
            WORK( I-LL+NM12 ) = COSL
            WORK( I-LL+NM13 ) = -SINL
          enddo
          E( LL ) = F
!
!           Test convergence
!
          IF( ABS( E( LL ) ).LE.THRESH )                                &
     &E( LL ) = ZERO
!
!           Update singular vectors if desired
!
          IF( NCVT.GT.0 )                                               &
     &CALL DLASR( 'L', 'V', 'B', M-LL+1, NCVT, WORK( NM12+1 ),          &
     &WORK( NM13+1 ), VT( LL, 1 ), LDVT )
          IF( NRU.GT.0 )                                                &
     &CALL DLASR( 'R', 'V', 'B', NRU, M-LL+1, WORK( 1 ),                &
     &WORK( N ), U( 1, LL ), LDU )
          IF( NCC.GT.0 )                                                &
     &CALL DLASR( 'L', 'V', 'B', M-LL+1, NCC, WORK( 1 ),                &
     &WORK( N ), C( LL, 1 ), LDC )
        END IF
      END IF
!
!     QR iteration finished, go back and check convergence
!
      GO TO 60
!
!     All singular values converged, so make them positive
!
  160 CONTINUE
      DO I = 1, N
        IF( D( I ).LT.ZERO ) THEN
          D( I ) = -D( I )
!
!           Change sign of singular vectors, if desired
!
          IF( NCVT.GT.0 )                                               &
     &CALL DSCAL( NCVT, NEGONE, VT( I, 1 ), LDVT )
        END IF
      enddo
!
!     Sort the singular values into decreasing order (insertion sort on
!     singular values, but only one transposition per singular vector)
!
      DO I = 1, N - 1
!
!        Scan for smallest D(I)
!
        ISUB = 1
        SMIN = D( 1 )
        DO J = 2, N + 1 - I
          IF( D( J ).LE.SMIN ) THEN
            ISUB = J
            SMIN = D( J )
          END IF
        enddo
        IF( ISUB.NE.N+1-I ) THEN
!
!           Swap singular values and vectors
!
          D( ISUB ) = D( N+1-I )
          D( N+1-I ) = SMIN
          IF( NCVT.GT.0 )                                               &
     &CALL DSWAP( NCVT, VT( ISUB, 1 ), LDVT, VT( N+1-I, 1 ),            &
     &LDVT )
          IF( NRU.GT.0 )                                                &
     &CALL DSWAP( NRU, U( 1, ISUB ), 1, U( 1, N+1-I ), 1 )
          IF( NCC.GT.0 )                                                &
     &CALL DSWAP( NCC, C( ISUB, 1 ), LDC, C( N+1-I, 1 ), LDC )
        END IF
      enddo
      GO TO 220
!
!     Maximum number of iterations exceeded, failure to converge
!
  200 CONTINUE
      INFO = 0
      DO I = 1, N - 1
        IF( E( I ).NE.ZERO )                                            &
     &INFO = INFO + 1
      enddo
  220 CONTINUE
      RETURN
!
!     End of DBDSQR
!
      END
      SUBROUTINE DCOPL(N,DX,INCX,DY,INCY)
      implicit none
!     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  Purpose
!  =======
!
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
      enddo
      RETURN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 M = MOD(N,7)
      IF (M.EQ.0) GO TO 40
      DO I = 1,M
        DY(I) = DX(I)
      enddo
      IF (N.LT.7) RETURN
   40 MP1 = M + 1
      DO I = MP1,N,7
        DY(I) = DX(I)
        DY(I+1) = DX(I+1)
        DY(I+2) = DX(I+2)
        DY(I+3) = DX(I+3)
        DY(I+4) = DX(I+4)
        DY(I+5) = DX(I+5)
        DY(I+6) = DX(I+6)
      enddo
      RETURN
      END
      SUBROUTINE DGEBD2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAUP( * ),        &
     &TAUQ( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGEBD2 reduces a real general m by n matrix A to upper or lower
!  bidiagonal form B by an orthogonal transformation: Q' * A * P = B.
!
!  If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows in the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns in the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the m by n general matrix to be reduced.
!          On exit,
!          if m >= n, the diagonal and the first superdiagonal are
!            overwritten with the upper bidiagonal matrix B; the
!            elements below the diagonal, with the array TAUQ, represent
!            the orthogonal matrix Q as a product of elementary
!            reflectors, and the elements above the first superdiagonal,
!            with the array TAUP, represent the orthogonal matrix P as
!            a product of elementary reflectors;
!          if m < n, the diagonal and the first subdiagonal are
!            overwritten with the lower bidiagonal matrix B; the
!            elements below the first subdiagonal, with the array TAUQ,
!            represent the orthogonal matrix Q as a product of
!            elementary reflectors, and the elements above the diagonal,
!            with the array TAUP, represent the orthogonal matrix P as
!            a product of elementary reflectors.
!          See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  D       (output) DOUBLE PRECISION array, dimension (min(M,N))
!          The diagonal elements of the bidiagonal matrix B:
!          D(i) = A(i,i).
!
!  E       (output) DOUBLE PRECISION array, dimension (min(M,N)-1)
!          The off-diagonal elements of the bidiagonal matrix B:
!          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
!          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
!
!  TAUQ    (output) DOUBLE PRECISION array dimension (min(M,N))
!          The scalar factors of the elementary reflectors which
!          represent the orthogonal matrix Q. See Further Details.
!
!  TAUP    (output) DOUBLE PRECISION array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors which
!          represent the orthogonal matrix P. See Further Details.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (max(M,N))
!
!  INFO    (output) INTEGER
!          = 0: successful exit.
!          < 0: if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The matrices Q and P are represented as products of elementary
!  reflectors:
!
!  If m >= n,
!
!     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)
!
!  Each H(i) and G(i) has the form:
!
!     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
!
!  where tauq and taup are real scalars, and v and u are real vectors;
!  v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i);
!  u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
!  tauq is stored in TAUQ(i) and taup in TAUP(i).
!
!  If m < n,
!
!     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)
!
!  Each H(i) and G(i) has the form:
!
!     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
!
!  where tauq and taup are real scalars, and v and u are real vectors;
!  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);
!  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
!  tauq is stored in TAUQ(i) and taup in TAUP(i).
!
!  The contents of A on exit are illustrated by the following examples:
!
!  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
!
!    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
!    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
!    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
!    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
!    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
!    (  v1  v2  v3  v4  v5 )
!
!  where d and e denote diagonal and off-diagonal elements of B, vi
!  denotes an element of the vector defining H(i), and ui an element of
!  the vector defining G(i).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARF, DLARFG, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      IF( M.LT.0 ) THEN
        INFO = -1
      ELSE IF( N.LT.0 ) THEN
        INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -4
      END IF
      IF( INFO.LT.0 ) THEN
        CALL XERBLA( 'DGEBD2', -INFO )
        RETURN
      END IF
!
      IF( M.GE.N ) THEN
!
!        Reduce to upper bidiagonal form
!
        DO I = 1, N
!
!           Generate elementary reflector H(i) to annihilate A(i+1:m,i)
!
          CALL DLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1,      &
     &TAUQ( I ) )
          D( I ) = A( I, I )
          A( I, I ) = ONE
!
!           Apply H(i) to A(i:m,i+1:n) from the left
!
          IF( I.LT.N )                                                  &
     &CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAUQ( I ),          &
     &A( I, I+1 ), LDA, WORK )
          A( I, I ) = D( I )
!
          IF( I.LT.N ) THEN
!
!              Generate elementary reflector G(i) to annihilate
!              A(i,i+2:n)
!
            CALL DLARFG( N-I, A( I, I+1 ), A( I, MIN( I+2, N ) ),       &
     &LDA, TAUP( I ) )
            E( I ) = A( I, I+1 )
            A( I, I+1 ) = ONE
!
!              Apply G(i) to A(i+1:m,i+1:n) from the right
!
            CALL DLARF( 'Right', M-I, N-I, A( I, I+1 ), LDA,            &
     &TAUP( I ), A( I+1, I+1 ), LDA, WORK )
            A( I, I+1 ) = E( I )
          ELSE
            TAUP( I ) = ZERO
          END IF
        enddo
      ELSE
!
!        Reduce to lower bidiagonal form
!
        DO I = 1, M
!
!           Generate elementary reflector G(i) to annihilate A(i,i+1:n)
!
          CALL DLARFG( N-I+1, A( I, I ), A( I, MIN( I+1, N ) ), LDA,    &
     &TAUP( I ) )
          D( I ) = A( I, I )
          A( I, I ) = ONE
!
!           Apply G(i) to A(i+1:m,i:n) from the right
!
          IF( I.LT.M )                                                  &
     &CALL DLARF( 'Right', M-I, N-I+1, A( I, I ), LDA,                  &
     &TAUP( I ), A( I+1, I ), LDA, WORK )
          A( I, I ) = D( I )
!
          IF( I.LT.M ) THEN
!
!              Generate elementary reflector H(i) to annihilate
!              A(i+2:m,i)
!
            CALL DLARFG( M-I, A( I+1, I ), A( MIN( I+2, M ), I ), 1,    &
     &TAUQ( I ) )
            E( I ) = A( I+1, I )
            A( I+1, I ) = ONE
!
!              Apply H(i) to A(i+1:m,i+1:n) from the left
!
            CALL DLARF( 'Left', M-I, N-I, A( I+1, I ), 1, TAUQ( I ),    &
     &A( I+1, I+1 ), LDA, WORK )
            A( I+1, I ) = E( I )
          ELSE
            TAUQ( I ) = ZERO
          END IF
        enddo
      END IF
      RETURN
!
!     End of DGEBD2
!
      END
      SUBROUTINE DGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK,   &
     &INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAUP( * ),        &
     &TAUQ( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGEBRD reduces a general real M-by-N matrix A to upper or lower
!  bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
!
!  If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows in the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns in the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N general matrix to be reduced.
!          On exit,
!          if m >= n, the diagonal and the first superdiagonal are
!            overwritten with the upper bidiagonal matrix B; the
!            elements below the diagonal, with the array TAUQ, represent
!            the orthogonal matrix Q as a product of elementary
!            reflectors, and the elements above the first superdiagonal,
!            with the array TAUP, represent the orthogonal matrix P as
!            a product of elementary reflectors;
!          if m < n, the diagonal and the first subdiagonal are
!            overwritten with the lower bidiagonal matrix B; the
!            elements below the first subdiagonal, with the array TAUQ,
!            represent the orthogonal matrix Q as a product of
!            elementary reflectors, and the elements above the diagonal,
!            with the array TAUP, represent the orthogonal matrix P as
!            a product of elementary reflectors.
!          See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  D       (output) DOUBLE PRECISION array, dimension (min(M,N))
!          The diagonal elements of the bidiagonal matrix B:
!          D(i) = A(i,i).
!
!  E       (output) DOUBLE PRECISION array, dimension (min(M,N)-1)
!          The off-diagonal elements of the bidiagonal matrix B:
!          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
!          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
!
!  TAUQ    (output) DOUBLE PRECISION array dimension (min(M,N))
!          The scalar factors of the elementary reflectors which
!          represent the orthogonal matrix Q. See Further Details.
!
!  TAUP    (output) DOUBLE PRECISION array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors which
!          represent the orthogonal matrix P. See Further Details.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The length of the array WORK.  LWORK >= max(1,M,N).
!          For optimum performance LWORK >= (M+N)*NB, where NB
!          is the optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The matrices Q and P are represented as products of elementary
!  reflectors:
!
!  If m >= n,
!
!     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)
!
!  Each H(i) and G(i) has the form:
!
!     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
!
!  where tauq and taup are real scalars, and v and u are real vectors;
!  v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i);
!  u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
!  tauq is stored in TAUQ(i) and taup in TAUP(i).
!
!  If m < n,
!
!     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)
!
!  Each H(i) and G(i) has the form:
!
!     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
!
!  where tauq and taup are real scalars, and v and u are real vectors;
!  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);
!  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
!  tauq is stored in TAUQ(i) and taup in TAUP(i).
!
!  The contents of A on exit are illustrated by the following examples:
!
!  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
!
!    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
!    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
!    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
!    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
!    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
!    (  v1  v2  v3  v4  v5 )
!
!  where d and e denote diagonal and off-diagonal elements of B, vi
!  denotes an element of the vector defining H(i), and ui an element of
!  the vector defining G(i).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IINFO, J, LDWRKX, LDWRKY, LWKOPT, MINMN, NB,&
     &NBMIN, NX
      DOUBLE PRECISION   WS
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEBD2, DGEMM, DLABRD, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      NB = MAX( 1, ILAENV( 1, 'DGEBRD', ' ', M, N, -1, -1 ) )
      LWKOPT = ( M+N )*NB
      WORK( 1 ) = DBLE( LWKOPT )
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
        INFO = -1
      ELSE IF( N.LT.0 ) THEN
        INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, M, N ) .AND. .NOT.LQUERY ) THEN
        INFO = -10
      END IF
      IF( INFO.LT.0 ) THEN
        CALL XERBLA( 'DGEBRD', -INFO )
        RETURN
      ELSE IF( LQUERY ) THEN
        RETURN
      END IF
!
!     Quick return if possible
!
      MINMN = MIN( M, N )
      IF( MINMN.EQ.0 ) THEN
        WORK( 1 ) = 1
        RETURN
      END IF
!
      WS = MAX( M, N )
      LDWRKX = M
      LDWRKY = N
!
      IF( NB.GT.1 .AND. NB.LT.MINMN ) THEN
!
!        Set the crossover point NX.
!
        NX = MAX( NB, ILAENV( 3, 'DGEBRD', ' ', M, N, -1, -1 ) )
!
!        Determine when to switch from blocked to unblocked code.
!
        IF( NX.LT.MINMN ) THEN
          WS = ( M+N )*NB
          IF( LWORK.LT.WS ) THEN
!
!              Not enough work space for the optimal NB, consider using
!              a smaller block size.
!
            NBMIN = ILAENV( 2, 'DGEBRD', ' ', M, N, -1, -1 )
            IF( LWORK.GE.( M+N )*NBMIN ) THEN
              NB = LWORK / ( M+N )
            ELSE
              NB = 1
              NX = MINMN
            END IF
          END IF
        END IF
      ELSE
        NX = MINMN
      END IF
!
      DO I = 1, MINMN - NX, NB
!
!        Reduce rows and columns i:i+nb-1 to bidiagonal form and return
!        the matrices X and Y which are needed to update the unreduced
!        part of the matrix
!
        CALL DLABRD( M-I+1, N-I+1, NB, A( I, I ), LDA, D( I ), E( I ),  &
     &TAUQ( I ), TAUP( I ), WORK, LDWRKX,                               &
     &WORK( LDWRKX*NB+1 ), LDWRKY )
!
!        Update the trailing submatrix A(i+nb:m,i+nb:n), using an update
!        of the form  A := A - V*Y' - X*U'
!
        CALL DGEMM( 'No transpose', 'Transpose', M-I-NB+1, N-I-NB+1,    &
     &NB, -ONE, A( I+NB, I ), LDA,                                      &
     &WORK( LDWRKX*NB+NB+1 ), LDWRKY, ONE,                              &
     &A( I+NB, I+NB ), LDA )
        CALL DGEMM( 'No transpose', 'No transpose', M-I-NB+1, N-I-NB+1, &
     &NB, -ONE, WORK( NB+1 ), LDWRKX, A( I, I+NB ), LDA,                &
     &ONE, A( I+NB, I+NB ), LDA )
!
!        Copy diagonal and off-diagonal elements of B back into A
!
        IF( M.GE.N ) THEN
          DO J = I, I + NB - 1
            A( J, J ) = D( J )
            A( J, J+1 ) = E( J )
          enddo
        ELSE
          DO J = I, I + NB - 1
            A( J, J ) = D( J )
            A( J+1, J ) = E( J )
          enddo
        END IF
      enddo
!
!     Use unblocked code to reduce the remainder of the matrix
!
      CALL DGEBD2( M-I+1, N-I+1, A( I, I ), LDA, D( I ), E( I ),        &
     &TAUQ( I ), TAUP( I ), WORK, IINFO )
      WORK( 1 ) = WS
      RETURN
!
!     End of DGEBRD
!
      END
      SUBROUTINE DGELQ2( M, N, A, LDA, TAU, WORK, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGELQ2 computes an LQ factorization of a real m by n matrix A:
!  A = L * Q.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the m by n matrix A.
!          On exit, the elements on and below the diagonal of the array
!          contain the m by min(m,n) lower trapezoidal matrix L (L is
!          lower triangular if m <= n); the elements above the diagonal,
!          with the array TAU, represent the orthogonal matrix Q as a
!          product of elementary reflectors (see Further Details).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (M)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  The matrix Q is represented as a product of elementary reflectors
!
!     Q = H(k) . . . H(2) H(1), where k = min(m,n).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n),
!  and tau in TAU(i).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, K
      DOUBLE PRECISION   AII
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARF, DLARFG, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M.LT.0 ) THEN
        INFO = -1
      ELSE IF( N.LT.0 ) THEN
        INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DGELQ2', -INFO )
        RETURN
      END IF
!
      K = MIN( M, N )
!
      DO I = 1, K
!
!        Generate elementary reflector H(i) to annihilate A(i,i+1:n)
!
        CALL DLARFG( N-I+1, A( I, I ), A( I, MIN( I+1, N ) ), LDA,      &
     &TAU( I ) )
        IF( I.LT.M ) THEN
!
!           Apply H(i) to A(i+1:m,i:n) from the right
!
          AII = A( I, I )
          A( I, I ) = ONE
          CALL DLARF( 'Right', M-I, N-I+1, A( I, I ), LDA, TAU( I ),    &
     &A( I+1, I ), LDA, WORK )
          A( I, I ) = AII
        END IF
      enddo
      RETURN
!
!     End of DGELQ2
!
      END
      SUBROUTINE DGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGELQF computes an LQ factorization of a real M-by-N matrix A:
!  A = L * Q.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit, the elements on and below the diagonal of the array
!          contain the m-by-min(m,n) lower trapezoidal matrix L (L is
!          lower triangular if m <= n); the elements above the diagonal,
!          with the array TAU, represent the orthogonal matrix Q as a
!          product of elementary reflectors (see Further Details).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= max(1,M).
!          For optimum performance LWORK >= M*NB, where NB is the
!          optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  The matrix Q is represented as a product of elementary reflectors
!
!     Q = H(k) . . . H(2) H(1), where k = min(m,n).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n),
!  and tau in TAU(i).
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB,      &
     &NBMIN, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGELQ2, DLARFB, DLARFT, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      NB = ILAENV( 1, 'DGELQF', ' ', M, N, -1, -1 )
      LWKOPT = M*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
        INFO = -1
      ELSE IF( N.LT.0 ) THEN
        INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, M ) .AND. .NOT.LQUERY ) THEN
        INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DGELQF', -INFO )
        RETURN
      ELSE IF( LQUERY ) THEN
        RETURN
      END IF
!
!     Quick return if possible
!
      K = MIN( M, N )
      IF( K.EQ.0 ) THEN
        WORK( 1 ) = 1
        RETURN
      END IF
!
      NBMIN = 2
      NX = 0
      IWS = M
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
        NX = MAX( 0, ILAENV( 3, 'DGELQF', ' ', M, N, -1, -1 ) )
        IF( NX.LT.K ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
          LDWORK = M
          IWS = LDWORK*NB
          IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DGELQF', ' ', M, N, -1,         &
     &-1 ) )
          END IF
        END IF
      END IF
!
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!
!        Use blocked code initially
!
        DO I = 1, K - NX, NB
          IB = MIN( K-I+1, NB )
!
!           Compute the LQ factorization of the current block
!           A(i:i+ib-1,i:n)
!
          CALL DGELQ2( IB, N-I+1, A( I, I ), LDA, TAU( I ), WORK,       &
     &IINFO )
          IF( I+IB.LE.M ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
            CALL DLARFT( 'Forward', 'Rowwise', N-I+1, IB, A( I, I ),    &
     &LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H to A(i+ib:m,i:n) from the right
!
            CALL DLARFB( 'Right', 'No transpose', 'Forward',            &
     &'Rowwise', M-I-IB+1, N-I+1, IB, A( I, I ),                        &
     &LDA, WORK, LDWORK, A( I+IB, I ), LDA,                             &
     &WORK( IB+1 ), LDWORK )
          END IF
        enddo
      ELSE
        I = 1
      END IF
!
!     Use unblocked code to factor the last or only block.
!
      IF( I.LE.K )                                                      &
     &CALL DGELQ2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK,        &
     &IINFO )
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of DGELQF
!
      END
      SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, &
     &INFO )
      implicit none
!
!  -- LAPACK driver routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGELS solves overdetermined or underdetermined real linear systems
!  involving an M-by-N matrix A, or its transpose, using a QR or LQ
!  factorization of A.  It is assumed that A has full rank.
!
!  The following options are provided:
!
!  1. If TRANS = 'N' and m >= n:  find the least squares solution of
!     an overdetermined system, i.e., solve the least squares problem
!                  minimize || B - A*X ||.
!
!  2. If TRANS = 'N' and m < n:  find the minimum norm solution of
!     an underdetermined system A * X = B.
!
!  3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
!     an undetermined system A**T * X = B.
!
!  4. If TRANS = 'T' and m < n:  find the least squares solution of
!     an overdetermined system, i.e., solve the least squares problem
!                  minimize || B - A**T * X ||.
!
!  Several right hand side vectors b and solution vectors x can be
!  handled in a single call; they are stored as the columns of the
!  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
!  matrix X.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          = 'N': the linear system involves A;
!          = 'T': the linear system involves A**T.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of
!          columns of the matrices B and X. NRHS >=0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit,
!            if M >= N, A is overwritten by details of its QR
!                       factorization as returned by DGEQRF;
!            if M <  N, A is overwritten by details of its LQ
!                       factorization as returned by DGELQF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the matrix B of right hand side vectors, stored
!          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
!          if TRANS = 'T'.
!          On exit, if INFO = 0, B is overwritten by the solution
!          vectors, stored columnwise:
!          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least
!          squares solution vectors; the residual sum of squares for the
!          solution in each column is given by the sum of squares of
!          elements N+1 to M in that column;
!          if TRANS = 'N' and m < n, rows 1 to N of B contain the
!          minimum norm solution vectors;
!          if TRANS = 'T' and m >= n, rows 1 to M of B contain the
!          minimum norm solution vectors;
!          if TRANS = 'T' and m < n, rows 1 to M of B contain the
!          least squares solution vectors; the residual sum of squares
!          for the solution in each column is given by the sum of
!          squares of elements M+1 to N in that column.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B. LDB >= MAX(1,M,N).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          LWORK >= max( 1, MN + max( MN, NRHS ) ).
!          For optimal performance,
!          LWORK >= max( 1, MN + max( MN, NRHS )*NB ).
!          where MN = min(M,N) and NB is the optimum block size.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO =  i, the i-th diagonal element of the
!                triangular factor of A is zero, so that A does not have
!                full rank; the least squares solution could not be
!                computed.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, TPSD
      INTEGER            BROW, I, IASCL, IBSCL, J, MN, NB, SCLLEN, WSIZE
      DOUBLE PRECISION   ANRM, BIGNUM, BNRM, SMLNUM
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   RWORK( 1 )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           LSAME, ILAENV, DLABAD, DLAMCH, DLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGELQF, DGEQRF, DLASCL, DLASET, DORMLQ, DORMQR,&
     &DTRTRS, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
!
      INFO = 0
      MN = MIN( M, N )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.( LSAME( TRANS, 'N' ) .OR. LSAME( TRANS, 'T' ) ) ) THEN
        INFO = -1
      ELSE IF( M.LT.0 ) THEN
        INFO = -2
      ELSE IF( N.LT.0 ) THEN
        INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
        INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -6
      ELSE IF( LDB.LT.MAX( 1, M, N ) ) THEN
        INFO = -8
      ELSE IF( LWORK.LT.MAX( 1, MN+MAX( MN, NRHS ) ) .AND. .NOT.LQUERY )&
     &THEN
        INFO = -10
      END IF
!
!     Figure out optimal block size
!
      IF( INFO.EQ.0 .OR. INFO.EQ.-10 ) THEN
!
        TPSD = .TRUE.
        IF( LSAME( TRANS, 'N' ) )                                       &
     &TPSD = .FALSE.
!
        IF( M.GE.N ) THEN
          NB = ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )
          IF( TPSD ) THEN
            NB = MAX( NB, ILAENV( 1, 'DORMQR', 'LN', M, NRHS, N,        &
     &-1 ) )
          ELSE
            NB = MAX( NB, ILAENV( 1, 'DORMQR', 'LT', M, NRHS, N,        &
     &-1 ) )
          END IF
        ELSE
          NB = ILAENV( 1, 'DGELQF', ' ', M, N, -1, -1 )
          IF( TPSD ) THEN
            NB = MAX( NB, ILAENV( 1, 'DORMLQ', 'LT', N, NRHS, M,        &
     &-1 ) )
          ELSE
            NB = MAX( NB, ILAENV( 1, 'DORMLQ', 'LN', N, NRHS, M,        &
     &-1 ) )
          END IF
        END IF
!
        WSIZE = MAX( 1, MN+MAX( MN, NRHS )*NB )
        WORK( 1 ) = DBLE( WSIZE )
!
      END IF
!
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DGELS ', -INFO )
        RETURN
      ELSE IF( LQUERY ) THEN
        RETURN
      END IF
!
!     Quick return if possible
!
      IF( MIN( M, N, NRHS ).EQ.0 ) THEN
        CALL DLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
        RETURN
      END IF
!
!     Get machine parameters
!
      SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
!
!     Scale A, B if max element outside range [SMLNUM,BIGNUM]
!
      ANRM = DLANGE( 'M', M, N, A, LDA, RWORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
!
!        Scale matrix norm up to SMLNUM
!
        CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
        IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN
!
!        Scale matrix norm down to BIGNUM
!
        CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
        IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN
!
!        Matrix all zero. Return zero solution.
!
        CALL DLASET( 'F', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
        GO TO 50
      END IF
!
      BROW = M
      IF( TPSD )                                                        &
     &BROW = N
      BNRM = DLANGE( 'M', BROW, NRHS, B, LDB, RWORK )
      IBSCL = 0
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
!
!        Scale matrix norm up to SMLNUM
!
        CALL DLASCL( 'G', 0, 0, BNRM, SMLNUM, BROW, NRHS, B, LDB,       &
     &INFO )
        IBSCL = 1
      ELSE IF( BNRM.GT.BIGNUM ) THEN
!
!        Scale matrix norm down to BIGNUM
!
        CALL DLASCL( 'G', 0, 0, BNRM, BIGNUM, BROW, NRHS, B, LDB,       &
     &INFO )
        IBSCL = 2
      END IF
!
      IF( M.GE.N ) THEN
!
!        compute QR factorization of A
!
        CALL DGEQRF( M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN,   &
     &INFO )
!
!        workspace at least N, optimally N*NB
!
        IF( .NOT.TPSD ) THEN
!
!           Least-Squares Problem min || A * X - B ||
!
!           B(1:M,1:NRHS) := Q' * B(1:M,1:NRHS)
!
          CALL DORMQR( 'Left', 'Transpose', M, NRHS, N, A, LDA,         &
     &WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN,                        &
     &INFO )
!
!           workspace at least NRHS, optimally NRHS*NB
!
!           B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)
!
          CALL DTRTRS( 'Upper', 'No transpose', 'Non-unit', N, NRHS,    &
     &A, LDA, B, LDB, INFO )
!
          IF( INFO.GT.0 ) THEN
            RETURN
          END IF
!
          SCLLEN = N
!
        ELSE
!
!           Overdetermined system of equations A' * X = B
!
!           B(1:N,1:NRHS) := inv(R') * B(1:N,1:NRHS)
!
          CALL DTRTRS( 'Upper', 'Transpose', 'Non-unit', N, NRHS,       &
     &A, LDA, B, LDB, INFO )
!
          IF( INFO.GT.0 ) THEN
            RETURN
          END IF
!
!           B(N+1:M,1:NRHS) = ZERO
!
          DO J = 1, NRHS
            DO I = N + 1, M
              B( I, J ) = ZERO
            enddo
          enddo
!
!           B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)
!
          CALL DORMQR( 'Left', 'No transpose', M, NRHS, N, A, LDA,      &
     &WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN,                        &
     &INFO )
!
!           workspace at least NRHS, optimally NRHS*NB
!
          SCLLEN = M
!
        END IF
!
      ELSE
!
!        Compute LQ factorization of A
!
        CALL DGELQF( M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN,   &
     &INFO )
!
!        workspace at least M, optimally M*NB.
!
        IF( .NOT.TPSD ) THEN
!
!           underdetermined system of equations A * X = B
!
!           B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)
!
          CALL DTRTRS( 'Lower', 'No transpose', 'Non-unit', M, NRHS,    &
     &A, LDA, B, LDB, INFO )
!
          IF( INFO.GT.0 ) THEN
            RETURN
          END IF
!
!           B(M+1:N,1:NRHS) = 0
!
          DO J = 1, NRHS
            DO I = M + 1, N
              B( I, J ) = ZERO
            enddo
          enddo
!
!           B(1:N,1:NRHS) := Q(1:N,:)' * B(1:M,1:NRHS)
!
          CALL DORMLQ( 'Left', 'Transpose', N, NRHS, M, A, LDA,         &
     &WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN,                        &
     &INFO )
!
!           workspace at least NRHS, optimally NRHS*NB
!
          SCLLEN = N
!
        ELSE
!
!           overdetermined system min || A' * X - B ||
!
!           B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)
!
          CALL DORMLQ( 'Left', 'No transpose', N, NRHS, M, A, LDA,      &
     &WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN,                        &
     &INFO )
!
!           workspace at least NRHS, optimally NRHS*NB
!
!           B(1:M,1:NRHS) := inv(L') * B(1:M,1:NRHS)
!
          CALL DTRTRS( 'Lower', 'Transpose', 'Non-unit', M, NRHS,       &
     &A, LDA, B, LDB, INFO )
!
          IF( INFO.GT.0 ) THEN
            RETURN
          END IF
!
          SCLLEN = M
!
        END IF
!
      END IF
!
!     Undo scaling
!
      IF( IASCL.EQ.1 ) THEN
        CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, SCLLEN, NRHS, B, LDB,     &
     &INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
        CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, SCLLEN, NRHS, B, LDB,     &
     &INFO )
      END IF
      IF( IBSCL.EQ.1 ) THEN
        CALL DLASCL( 'G', 0, 0, SMLNUM, BNRM, SCLLEN, NRHS, B, LDB,     &
     &INFO )
      ELSE IF( IBSCL.EQ.2 ) THEN
        CALL DLASCL( 'G', 0, 0, BIGNUM, BNRM, SCLLEN, NRHS, B, LDB,     &
     &INFO )
      END IF
!
   50 CONTINUE
      WORK( 1 ) = DBLE( WSIZE )
!
      RETURN
!
!     End of DGELS
!
      END
      SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      implicit none
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
!     ..
!
!  Purpose
!  =======
!
!  DGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Arguments
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
      LOGICAL NOTA,NOTB
!     ..
!     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      IF (NOTA) THEN
        NROWA = M
        NCOLA = K
      ELSE
        NROWA = K
        NCOLA = M
      END IF
      IF (NOTB) THEN
        NROWB = K
      ELSE
        NROWB = N
      END IF
!
!     Test the input parameters.
!
      INFO = 0
      IF ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND.              &
     &(.NOT.LSAME(TRANSA,'T'))) THEN
        INFO = 1
      ELSE IF ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND.         &
     &(.NOT.LSAME(TRANSB,'T'))) THEN
        INFO = 2
      ELSE IF (M.LT.0) THEN
        INFO = 3
      ELSE IF (N.LT.0) THEN
        INFO = 4
      ELSE IF (K.LT.0) THEN
        INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
        INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
        INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
        INFO = 13
      END IF
      IF (INFO.NE.0) THEN
        CALL XERBLA('DGEMM ',INFO)
        RETURN
      END IF
!
!     Quick return if possible.
!
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.                                   &
     &(((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
!
!     And if  alpha.eq.zero.
!
      IF (ALPHA.EQ.ZERO) THEN
        IF (BETA.EQ.ZERO) THEN
          DO J = 1,N
            DO I = 1,M
              C(I,J) = ZERO
            enddo
          enddo
        ELSE
          DO J = 1,N
            DO I = 1,M
              C(I,J) = BETA*C(I,J)
            enddo
          enddo
        END IF
        RETURN
      END IF
!
!     Start the operations.
!
      IF (NOTB) THEN
        IF (NOTA) THEN
!
!           Form  C := alpha*A*B + beta*C.
!
          DO J = 1,N
            IF (BETA.EQ.ZERO) THEN
              DO I = 1,M
                C(I,J) = ZERO
              enddo
            ELSE IF (BETA.NE.ONE) THEN
              DO I = 1,M
                C(I,J) = BETA*C(I,J)
              enddo
            END IF
            DO L = 1,K
              IF (B(L,J).NE.ZERO) THEN
                TEMP = ALPHA*B(L,J)
                DO I = 1,M
                  C(I,J) = C(I,J) + TEMP*A(I,L)
                enddo
              END IF
            enddo
          enddo
        ELSE
!
!           Form  C := alpha*A'*B + beta*C
!
          DO J = 1,N
            DO I = 1,M
              TEMP = ZERO
              DO L = 1,K
                TEMP = TEMP + A(L,I)*B(L,J)
              enddo
              IF (BETA.EQ.ZERO) THEN
                C(I,J) = ALPHA*TEMP
              ELSE
                C(I,J) = ALPHA*TEMP + BETA*C(I,J)
              END IF
            enddo
          enddo
        END IF
      ELSE
        IF (NOTA) THEN
!
!           Form  C := alpha*A*B' + beta*C
!
          DO J = 1,N
            IF (BETA.EQ.ZERO) THEN
              DO I = 1,M
                C(I,J) = ZERO
              enddo
            ELSE IF (BETA.NE.ONE) THEN
              DO I = 1,M
                C(I,J) = BETA*C(I,J)
              enddo
            END IF
            DO L = 1,K
              IF (B(J,L).NE.ZERO) THEN
                TEMP = ALPHA*B(J,L)
                DO I = 1,M
                  C(I,J) = C(I,J) + TEMP*A(I,L)
                enddo
              END IF
            enddo
          enddo
        ELSE
!
!           Form  C := alpha*A'*B' + beta*C
!
          DO J = 1,N
            DO I = 1,M
              TEMP = ZERO
              DO L = 1,K
                TEMP = TEMP + A(L,I)*B(J,L)
              enddo
              IF (BETA.EQ.ZERO) THEN
                C(I,J) = ALPHA*TEMP
              ELSE
                C(I,J) = ALPHA*TEMP + BETA*C(I,J)
              END IF
            enddo
          enddo
        END IF
      END IF
!
      RETURN
!
!     End of DGEMM .
!
      END
      SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      implicit none
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!     ..
!
!  Purpose
!  =======
!
!  DGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Arguments
!  ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.       &
     &.NOT.LSAME(TRANS,'C')) THEN
        INFO = 1
      ELSE IF (M.LT.0) THEN
        INFO = 2
      ELSE IF (N.LT.0) THEN
        INFO = 3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
        INFO = 6
      ELSE IF (INCX.EQ.0) THEN
        INFO = 8
      ELSE IF (INCY.EQ.0) THEN
        INFO = 11
      END IF
      IF (INFO.NE.0) THEN
        CALL XERBLA('DGEMV ',INFO)
        RETURN
      END IF
!
!     Quick return if possible.
!
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.                                   &
     &((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      IF (LSAME(TRANS,'N')) THEN
        LENX = N
        LENY = M
      ELSE
        LENX = M
        LENY = N
      END IF
      IF (INCX.GT.0) THEN
        KX = 1
      ELSE
        KX = 1 - (LENX-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
        KY = 1
      ELSE
        KY = 1 - (LENY-1)*INCY
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
      IF (BETA.NE.ONE) THEN
        IF (INCY.EQ.1) THEN
          IF (BETA.EQ.ZERO) THEN
            DO I = 1,LENY
              Y(I) = ZERO
            enddo
          ELSE
            DO I = 1,LENY
              Y(I) = BETA*Y(I)
            enddo
          END IF
        ELSE
          IY = KY
          IF (BETA.EQ.ZERO) THEN
            DO I = 1,LENY
              Y(IY) = ZERO
              IY = IY + INCY
            enddo
          ELSE
            DO I = 1,LENY
              Y(IY) = BETA*Y(IY)
              IY = IY + INCY
            enddo
          END IF
        END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(TRANS,'N')) THEN
!
!        Form  y := alpha*A*x + y.
!
        JX = KX
        IF (INCY.EQ.1) THEN
          DO J = 1,N
            IF (X(JX).NE.ZERO) THEN
              TEMP = ALPHA*X(JX)
              DO I = 1,M
                Y(I) = Y(I) + TEMP*A(I,J)
              enddo
            END IF
            JX = JX + INCX
          enddo
        ELSE
          DO J = 1,N
            IF (X(JX).NE.ZERO) THEN
              TEMP = ALPHA*X(JX)
              IY = KY
              DO I = 1,M
                Y(IY) = Y(IY) + TEMP*A(I,J)
                IY = IY + INCY
              enddo
            END IF
            JX = JX + INCX
          enddo
        END IF
      ELSE
!
!        Form  y := alpha*A'*x + y.
!
        JY = KY
        IF (INCX.EQ.1) THEN
          DO J = 1,N
            TEMP = ZERO
            DO I = 1,M
              TEMP = TEMP + A(I,J)*X(I)
            enddo
            Y(JY) = Y(JY) + ALPHA*TEMP
            JY = JY + INCY
          enddo
        ELSE
          DO J = 1,N
            TEMP = ZERO
            IX = KX
            DO I = 1,M
              TEMP = TEMP + A(I,J)*X(IX)
              IX = IX + INCX
            enddo
            Y(JY) = Y(JY) + ALPHA*TEMP
            JY = JY + INCY
          enddo
        END IF
      END IF
!
      RETURN
!
!     End of DGEMV .
!
      END
      SUBROUTINE DGEQR2( M, N, A, LDA, TAU, WORK, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGEQR2 computes a QR factorization of a real m by n matrix A:
!  A = Q * R.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the m by n matrix A.
!          On exit, the elements on and above the diagonal of the array
!          contain the min(m,n) by n upper trapezoidal matrix R (R is
!          upper triangular if m >= n); the elements below the diagonal,
!          with the array TAU, represent the orthogonal matrix Q as a
!          product of elementary reflectors (see Further Details).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  The matrix Q is represented as a product of elementary reflectors
!
!     Q = H(1) H(2) . . . H(k), where k = min(m,n).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
!  and tau in TAU(i).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, K
      DOUBLE PRECISION   AII
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARF, DLARFG, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M.LT.0 ) THEN
        INFO = -1
      ELSE IF( N.LT.0 ) THEN
        INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DGEQR2', -INFO )
        RETURN
      END IF
!
      K = MIN( M, N )
!
      DO I = 1, K
!
!        Generate elementary reflector H(i) to annihilate A(i+1:m,i)
!
        CALL DLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1,        &
     &TAU( I ) )
        IF( I.LT.N ) THEN
!
!           Apply H(i) to A(i:m,i+1:n) from the left
!
          AII = A( I, I )
          A( I, I ) = ONE
          CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),       &
     &A( I, I+1 ), LDA, WORK )
          A( I, I ) = AII
        END IF
      enddo
      RETURN
!
!     End of DGEQR2
!
      END
      SUBROUTINE DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGEQRF computes a QR factorization of a real M-by-N matrix A:
!  A = Q * R.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit, the elements on and above the diagonal of the array
!          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
!          upper triangular if m >= n); the elements below the diagonal,
!          with the array TAU, represent the orthogonal matrix Q as a
!          product of min(m,n) elementary reflectors (see Further
!          Details).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= max(1,N).
!          For optimum performance LWORK >= N*NB, where NB is
!          the optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  The matrix Q is represented as a product of elementary reflectors
!
!     Q = H(1) H(2) . . . H(k), where k = min(m,n).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
!  and tau in TAU(i).
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB,      &
     &NBMIN, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEQR2, DLARFB, DLARFT, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      NB = ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
        INFO = -1
      ELSE IF( N.LT.0 ) THEN
        INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
        INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DGEQRF', -INFO )
        RETURN
      ELSE IF( LQUERY ) THEN
        RETURN
      END IF
!
!     Quick return if possible
!
      K = MIN( M, N )
      IF( K.EQ.0 ) THEN
        WORK( 1 ) = 1
        RETURN
      END IF
!
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
        NX = MAX( 0, ILAENV( 3, 'DGEQRF', ' ', M, N, -1, -1 ) )
        IF( NX.LT.K ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
          LDWORK = N
          IWS = LDWORK*NB
          IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DGEQRF', ' ', M, N, -1,         &
     &-1 ) )
          END IF
        END IF
      END IF
!
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!
!        Use blocked code initially
!
        DO I = 1, K - NX, NB
          IB = MIN( K-I+1, NB )
!
!           Compute the QR factorization of the current block
!           A(i:m,i:i+ib-1)
!
          CALL DGEQR2( M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK,       &
     &IINFO )
          IF( I+IB.LE.N ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
            CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB,            &
     &A( I, I ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H' to A(i:m,i+ib:n) from the left
!
            CALL DLARFB( 'Left', 'Transpose', 'Forward',                &
     &'Columnwise', M-I+1, N-I-IB+1, IB,                                &
     &A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),                       &
     &LDA, WORK( IB+1 ), LDWORK )
          END IF
        enddo
      ELSE
        I = 1
      END IF
!
!     Use unblocked code to factor the last or only block.
!
      IF( I.LE.K )                                                      &
     &CALL DGEQR2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK,        &
     &IINFO )
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of DGEQRF
!
      END
      SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      implicit none
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER INCX,INCY,LDA,M,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!     ..
!
!  Purpose
!  =======
!
!  DGER   performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Arguments
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JY,KX
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      IF (M.LT.0) THEN
        INFO = 1
      ELSE IF (N.LT.0) THEN
        INFO = 2
      ELSE IF (INCX.EQ.0) THEN
        INFO = 5
      ELSE IF (INCY.EQ.0) THEN
        INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
        INFO = 9
      END IF
      IF (INFO.NE.0) THEN
        CALL XERBLA('DGER  ',INFO)
        RETURN
      END IF
!
!     Quick return if possible.
!
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF (INCY.GT.0) THEN
        JY = 1
      ELSE
        JY = 1 - (N-1)*INCY
      END IF
      IF (INCX.EQ.1) THEN
        DO J = 1,N
          IF (Y(JY).NE.ZERO) THEN
            TEMP = ALPHA*Y(JY)
            DO I = 1,M
              A(I,J) = A(I,J) + X(I)*TEMP
            enddo
          END IF
          JY = JY + INCY
        enddo
      ELSE
        IF (INCX.GT.0) THEN
          KX = 1
        ELSE
          KX = 1 - (M-1)*INCX
        END IF
        DO J = 1,N
          IF (Y(JY).NE.ZERO) THEN
            TEMP = ALPHA*Y(JY)
            IX = KX
            DO I = 1,M
              A(I,J) = A(I,J) + X(IX)*TEMP
              IX = IX + INCX
            enddo
          END IF
          JY = JY + INCY
        enddo
      END IF
!
      RETURN
!
!     End of DGER  .
!
      END
      SUBROUTINE DGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, &
     &LWORK, IWORK, INFO )
      implicit none
!
!  -- LAPACK driver routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ
      INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ),              &
     &VT( LDVT, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGESDD computes the singular value decomposition (SVD) of a real
!  M-by-N matrix A, optionally computing the left and right singular
!  vectors.  If singular vectors are desired, it uses a
!  divide-and-conquer algorithm.
!
!  The SVD is written
!
!       A = U * SIGMA * transpose(V)
!
!  where SIGMA is an M-by-N matrix which is zero except for its
!  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
!  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
!  are the singular values of A; they are real and non-negative, and
!  are returned in descending order.  The first min(m,n) columns of
!  U and V are the left and right singular vectors of A.
!
!  Note that the routine returns VT = V**T, not V.
!
!  The divide and conquer algorithm makes very mild assumptions about
!  floating point arithmetic. It will work on machines with a guard
!  digit in add/subtract, or on those binary machines without guard
!  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
!  Cray-2. It could conceivably fail on hexadecimal or decimal machines
!  without guard digits, but we know of none.
!
!  Arguments
!  =========
!
!  JOBZ    (input) CHARACTER*1
!          Specifies options for computing all or part of the matrix U:
!          = 'A':  all M columns of U and all N rows of V**T are
!                  returned in the arrays U and VT;
!          = 'S':  the first min(M,N) columns of U and the first
!                  min(M,N) rows of V**T are returned in the arrays U
!                  and VT;
!          = 'O':  If M >= N, the first N columns of U are overwritten
!                  on the array A and all rows of V**T are returned in
!                  the array VT;
!                  otherwise, all columns of U are returned in the
!                  array U and the first M rows of V**T are overwritten
!                  in the array A;
!          = 'N':  no columns of U or rows of V**T are computed.
!
!  M       (input) INTEGER
!          The number of rows of the input matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the input matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit,
!          if JOBZ = 'O',  A is overwritten with the first N columns
!                          of U (the left singular vectors, stored
!                          columnwise) if M >= N;
!                          A is overwritten with the first M rows
!                          of V**T (the right singular vectors, stored
!                          rowwise) otherwise.
!          if JOBZ .ne. 'O', the contents of A are destroyed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
!          The singular values of A, sorted so that S(i) >= S(i+1).
!
!  U       (output) DOUBLE PRECISION array, dimension (LDU,UCOL)
!          UCOL = M if JOBZ = 'A' or JOBZ = 'O' and M < N;
!          UCOL = min(M,N) if JOBZ = 'S'.
!          If JOBZ = 'A' or JOBZ = 'O' and M < N, U contains the M-by-M
!          orthogonal matrix U;
!          if JOBZ = 'S', U contains the first min(M,N) columns of U
!          (the left singular vectors, stored columnwise);
!          if JOBZ = 'O' and M >= N, or JOBZ = 'N', U is not referenced.
!
!  LDU     (input) INTEGER
!          The leading dimension of the array U.  LDU >= 1; if
!          JOBZ = 'S' or 'A' or JOBZ = 'O' and M < N, LDU >= M.
!
!  VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)
!          If JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the
!          N-by-N orthogonal matrix V**T;
!          if JOBZ = 'S', VT contains the first min(M,N) rows of
!          V**T (the right singular vectors, stored rowwise);
!          if JOBZ = 'O' and M < N, or JOBZ = 'N', VT is not referenced.
!
!  LDVT    (input) INTEGER
!          The leading dimension of the array VT.  LDVT >= 1; if
!          JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N;
!          if JOBZ = 'S', LDVT >= min(M,N).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= 1.
!          If JOBZ = 'N',
!            LWORK >= 3*min(M,N) + max(max(M,N),7*min(M,N)).
!          If JOBZ = 'O',
!            LWORK >= 3*min(M,N)*min(M,N) +
!                     max(max(M,N),5*min(M,N)*min(M,N)+4*min(M,N)).
!          If JOBZ = 'S' or 'A'
!            LWORK >= 3*min(M,N)*min(M,N) +
!                     max(max(M,N),4*min(M,N)*min(M,N)+4*min(M,N)).
!          For good performance, LWORK should generally be larger.
!          If LWORK = -1 but other input arguments are legal, WORK(1)
!          returns the optimal LWORK.
!
!  IWORK   (workspace) INTEGER array, dimension (8*min(M,N))
!
!  INFO    (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  DBDSDC did not converge, updating process failed.
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ming Gu and Huan Ren, Computer Science Division, University of
!     California at Berkeley, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, WNTQA, WNTQAS, WNTQN, WNTQO, WNTQS
      INTEGER            BDSPAC, BLK, CHUNK, I, IE, IERR, IL,           &
     &IR, ISCL, ITAU, ITAUP, ITAUQ, IU, IVT, LDWKVT,                    &
     &LDWRKL, LDWRKR, LDWRKU, MAXWRK, MINMN, MINWRK,                    &
     &MNTHR, NWORK, WRKBL
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, SMLNUM
!     ..
!     .. Local Arrays ..
      INTEGER            IDUM( 1 )
      DOUBLE PRECISION   DUM( 1 )
!     ..
!     .. External Subroutines ..
      EXTERNAL           DBDSDC, DGEBRD, DGELQF, DGEMM, DGEQRF, DLACPY, &
     &DLASCL, DLASET, DORGBR, DORGLQ, DORGQR, DORMBR,                   &
     &XERBLA
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           DLAMCH, DLANGE, ILAENV, LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      MINMN = MIN( M, N )
      WNTQA = LSAME( JOBZ, 'A' )
      WNTQS = LSAME( JOBZ, 'S' )
      WNTQAS = WNTQA .OR. WNTQS
      WNTQO = LSAME( JOBZ, 'O' )
      WNTQN = LSAME( JOBZ, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
!
      IF( .NOT.( WNTQA .OR. WNTQS .OR. WNTQO .OR. WNTQN ) ) THEN
        INFO = -1
      ELSE IF( M.LT.0 ) THEN
        INFO = -2
      ELSE IF( N.LT.0 ) THEN
        INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -5
      ELSE IF( LDU.LT.1 .OR. ( WNTQAS .AND. LDU.LT.M ) .OR.             &
     &( WNTQO .AND. M.LT.N .AND. LDU.LT.M ) ) THEN
        INFO = -8
      ELSE IF( LDVT.LT.1 .OR. ( WNTQA .AND. LDVT.LT.N ) .OR.            &
     &( WNTQS .AND. LDVT.LT.MINMN ) .OR.                                &
     &( WNTQO .AND. M.GE.N .AND. LDVT.LT.N ) ) THEN
        INFO = -10
      END IF
!
!     Compute workspace
!      (Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of workspace needed at that point in the code,
!       as well as the preferred amount for good performance.
!       NB refers to the optimal block size for the immediately
!       following subroutine, as returned by ILAENV.)
!
      IF( INFO.EQ.0 ) THEN
        MINWRK = 1
        MAXWRK = 1
        IF( M.GE.N .AND. MINMN.GT.0 ) THEN
!
!           Compute space needed for DBDSDC
!
          MNTHR = INT( MINMN*11.0D0 / 6.0D0 )
          IF( WNTQN ) THEN
            BDSPAC = 7*N
          ELSE
            BDSPAC = 3*N*N + 4*N
          END IF
          IF( M.GE.MNTHR ) THEN
            IF( WNTQN ) THEN
!
!                 Path 1 (M much larger than N, JOBZ='N')
!
              WRKBL = N + N*ILAENV( 1, 'DGEQRF', ' ', M, N, -1,         &
     &-1 )
              WRKBL = MAX( WRKBL, 3*N+2*N*                              &
     &ILAENV( 1, 'DGEBRD', ' ', N, N, -1, -1 ) )
              MAXWRK = MAX( WRKBL, BDSPAC+N )
              MINWRK = BDSPAC + N
            ELSE IF( WNTQO ) THEN
!
!                 Path 2 (M much larger than N, JOBZ='O')
!
              WRKBL = N + N*ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )
              WRKBL = MAX( WRKBL, N+N*ILAENV( 1, 'DORGQR', ' ', M,      &
     &N, N, -1 ) )
              WRKBL = MAX( WRKBL, 3*N+2*N*                              &
     &ILAENV( 1, 'DGEBRD', ' ', N, N, -1, -1 ) )
              WRKBL = MAX( WRKBL, 3*N+N*                                &
     &ILAENV( 1, 'DORMBR', 'QLN', N, N, N, -1 ) )
              WRKBL = MAX( WRKBL, 3*N+N*                                &
     &ILAENV( 1, 'DORMBR', 'PRT', N, N, N, -1 ) )
              WRKBL = MAX( WRKBL, BDSPAC+3*N )
              MAXWRK = WRKBL + 2*N*N
              MINWRK = BDSPAC + 2*N*N + 3*N
            ELSE IF( WNTQS ) THEN
!
!                 Path 3 (M much larger than N, JOBZ='S')
!
              WRKBL = N + N*ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )
              WRKBL = MAX( WRKBL, N+N*ILAENV( 1, 'DORGQR', ' ', M,      &
     &N, N, -1 ) )
              WRKBL = MAX( WRKBL, 3*N+2*N*                              &
     &ILAENV( 1, 'DGEBRD', ' ', N, N, -1, -1 ) )
              WRKBL = MAX( WRKBL, 3*N+N*                                &
     &ILAENV( 1, 'DORMBR', 'QLN', N, N, N, -1 ) )
              WRKBL = MAX( WRKBL, 3*N+N*                                &
     &ILAENV( 1, 'DORMBR', 'PRT', N, N, N, -1 ) )
              WRKBL = MAX( WRKBL, BDSPAC+3*N )
              MAXWRK = WRKBL + N*N
              MINWRK = BDSPAC + N*N + 3*N
            ELSE IF( WNTQA ) THEN
!
!                 Path 4 (M much larger than N, JOBZ='A')
!
              WRKBL = N + N*ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )
              WRKBL = MAX( WRKBL, N+M*ILAENV( 1, 'DORGQR', ' ', M,      &
     &M, N, -1 ) )
              WRKBL = MAX( WRKBL, 3*N+2*N*                              &
     &ILAENV( 1, 'DGEBRD', ' ', N, N, -1, -1 ) )
              WRKBL = MAX( WRKBL, 3*N+N*                                &
     &ILAENV( 1, 'DORMBR', 'QLN', N, N, N, -1 ) )
              WRKBL = MAX( WRKBL, 3*N+N*                                &
     &ILAENV( 1, 'DORMBR', 'PRT', N, N, N, -1 ) )
              WRKBL = MAX( WRKBL, BDSPAC+3*N )
              MAXWRK = WRKBL + N*N
              MINWRK = BDSPAC + N*N + 3*N
            END IF
          ELSE
!
!              Path 5 (M at least N, but not much larger)
!
            WRKBL = 3*N + ( M+N )*ILAENV( 1, 'DGEBRD', ' ', M, N, -1,   &
     &-1 )
            IF( WNTQN ) THEN
              MAXWRK = MAX( WRKBL, BDSPAC+3*N )
              MINWRK = 3*N + MAX( M, BDSPAC )
            ELSE IF( WNTQO ) THEN
              WRKBL = MAX( WRKBL, 3*N+N*                                &
     &ILAENV( 1, 'DORMBR', 'QLN', M, N, N, -1 ) )
              WRKBL = MAX( WRKBL, 3*N+N*                                &
     &ILAENV( 1, 'DORMBR', 'PRT', N, N, N, -1 ) )
              WRKBL = MAX( WRKBL, BDSPAC+3*N )
              MAXWRK = WRKBL + M*N
              MINWRK = 3*N + MAX( M, N*N+BDSPAC )
            ELSE IF( WNTQS ) THEN
              WRKBL = MAX( WRKBL, 3*N+N*                                &
     &ILAENV( 1, 'DORMBR', 'QLN', M, N, N, -1 ) )
              WRKBL = MAX( WRKBL, 3*N+N*                                &
     &ILAENV( 1, 'DORMBR', 'PRT', N, N, N, -1 ) )
              MAXWRK = MAX( WRKBL, BDSPAC+3*N )
              MINWRK = 3*N + MAX( M, BDSPAC )
            ELSE IF( WNTQA ) THEN
              WRKBL = MAX( WRKBL, 3*N+M*                                &
     &ILAENV( 1, 'DORMBR', 'QLN', M, M, N, -1 ) )
              WRKBL = MAX( WRKBL, 3*N+N*                                &
     &ILAENV( 1, 'DORMBR', 'PRT', N, N, N, -1 ) )
              MAXWRK = MAX( MAXWRK, BDSPAC+3*N )
              MINWRK = 3*N + MAX( M, BDSPAC )
            END IF
          END IF
        ELSE IF( MINMN.GT.0 ) THEN
!
!           Compute space needed for DBDSDC
!
          MNTHR = INT( MINMN*11.0D0 / 6.0D0 )
          IF( WNTQN ) THEN
            BDSPAC = 7*M
          ELSE
            BDSPAC = 3*M*M + 4*M
          END IF
          IF( N.GE.MNTHR ) THEN
            IF( WNTQN ) THEN
!
!                 Path 1t (N much larger than M, JOBZ='N')
!
              WRKBL = M + M*ILAENV( 1, 'DGELQF', ' ', M, N, -1,         &
     &-1 )
              WRKBL = MAX( WRKBL, 3*M+2*M*                              &
     &ILAENV( 1, 'DGEBRD', ' ', M, M, -1, -1 ) )
              MAXWRK = MAX( WRKBL, BDSPAC+M )
              MINWRK = BDSPAC + M
            ELSE IF( WNTQO ) THEN
!
!                 Path 2t (N much larger than M, JOBZ='O')
!
              WRKBL = M + M*ILAENV( 1, 'DGELQF', ' ', M, N, -1, -1 )
              WRKBL = MAX( WRKBL, M+M*ILAENV( 1, 'DORGLQ', ' ', M,      &
     &N, M, -1 ) )
              WRKBL = MAX( WRKBL, 3*M+2*M*                              &
     &ILAENV( 1, 'DGEBRD', ' ', M, M, -1, -1 ) )
              WRKBL = MAX( WRKBL, 3*M+M*                                &
     &ILAENV( 1, 'DORMBR', 'QLN', M, M, M, -1 ) )
              WRKBL = MAX( WRKBL, 3*M+M*                                &
     &ILAENV( 1, 'DORMBR', 'PRT', M, M, M, -1 ) )
              WRKBL = MAX( WRKBL, BDSPAC+3*M )
              MAXWRK = WRKBL + 2*M*M
              MINWRK = BDSPAC + 2*M*M + 3*M
            ELSE IF( WNTQS ) THEN
!
!                 Path 3t (N much larger than M, JOBZ='S')
!
              WRKBL = M + M*ILAENV( 1, 'DGELQF', ' ', M, N, -1, -1 )
              WRKBL = MAX( WRKBL, M+M*ILAENV( 1, 'DORGLQ', ' ', M,      &
     &N, M, -1 ) )
              WRKBL = MAX( WRKBL, 3*M+2*M*                              &
     &ILAENV( 1, 'DGEBRD', ' ', M, M, -1, -1 ) )
              WRKBL = MAX( WRKBL, 3*M+M*                                &
     &ILAENV( 1, 'DORMBR', 'QLN', M, M, M, -1 ) )
              WRKBL = MAX( WRKBL, 3*M+M*                                &
     &ILAENV( 1, 'DORMBR', 'PRT', M, M, M, -1 ) )
              WRKBL = MAX( WRKBL, BDSPAC+3*M )
              MAXWRK = WRKBL + M*M
              MINWRK = BDSPAC + M*M + 3*M
            ELSE IF( WNTQA ) THEN
!
!                 Path 4t (N much larger than M, JOBZ='A')
!
              WRKBL = M + M*ILAENV( 1, 'DGELQF', ' ', M, N, -1, -1 )
              WRKBL = MAX( WRKBL, M+N*ILAENV( 1, 'DORGLQ', ' ', N,      &
     &N, M, -1 ) )
              WRKBL = MAX( WRKBL, 3*M+2*M*                              &
     &ILAENV( 1, 'DGEBRD', ' ', M, M, -1, -1 ) )
              WRKBL = MAX( WRKBL, 3*M+M*                                &
     &ILAENV( 1, 'DORMBR', 'QLN', M, M, M, -1 ) )
              WRKBL = MAX( WRKBL, 3*M+M*                                &
     &ILAENV( 1, 'DORMBR', 'PRT', M, M, M, -1 ) )
              WRKBL = MAX( WRKBL, BDSPAC+3*M )
              MAXWRK = WRKBL + M*M
              MINWRK = BDSPAC + M*M + 3*M
            END IF
          ELSE
!
!              Path 5t (N greater than M, but not much larger)
!
            WRKBL = 3*M + ( M+N )*ILAENV( 1, 'DGEBRD', ' ', M, N, -1,   &
     &-1 )
            IF( WNTQN ) THEN
              MAXWRK = MAX( WRKBL, BDSPAC+3*M )
              MINWRK = 3*M + MAX( N, BDSPAC )
            ELSE IF( WNTQO ) THEN
              WRKBL = MAX( WRKBL, 3*M+M*                                &
     &ILAENV( 1, 'DORMBR', 'QLN', M, M, N, -1 ) )
              WRKBL = MAX( WRKBL, 3*M+M*                                &
     &ILAENV( 1, 'DORMBR', 'PRT', M, N, M, -1 ) )
              WRKBL = MAX( WRKBL, BDSPAC+3*M )
              MAXWRK = WRKBL + M*N
              MINWRK = 3*M + MAX( N, M*M+BDSPAC )
            ELSE IF( WNTQS ) THEN
              WRKBL = MAX( WRKBL, 3*M+M*                                &
     &ILAENV( 1, 'DORMBR', 'QLN', M, M, N, -1 ) )
              WRKBL = MAX( WRKBL, 3*M+M*                                &
     &ILAENV( 1, 'DORMBR', 'PRT', M, N, M, -1 ) )
              MAXWRK = MAX( WRKBL, BDSPAC+3*M )
              MINWRK = 3*M + MAX( N, BDSPAC )
            ELSE IF( WNTQA ) THEN
              WRKBL = MAX( WRKBL, 3*M+M*                                &
     &ILAENV( 1, 'DORMBR', 'QLN', M, M, N, -1 ) )
              WRKBL = MAX( WRKBL, 3*M+M*                                &
     &ILAENV( 1, 'DORMBR', 'PRT', N, N, M, -1 ) )
              MAXWRK = MAX( WRKBL, BDSPAC+3*M )
              MINWRK = 3*M + MAX( N, BDSPAC )
            END IF
          END IF
        END IF
        MAXWRK = MAX( MAXWRK, MINWRK )
        WORK( 1 ) = MAXWRK
!
        IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
          INFO = -12
        END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DGESDD', -INFO )
        RETURN
      ELSE IF( LQUERY ) THEN
        RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
        RETURN
      END IF
!
!     Get machine constants
!
      EPS = DLAMCH( 'P' )
      SMLNUM = SQRT( DLAMCH( 'S' ) ) / EPS
      BIGNUM = ONE / SMLNUM
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
      ANRM = DLANGE( 'M', M, N, A, LDA, DUM )
      ISCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
        ISCL = 1
        CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, IERR )
      ELSE IF( ANRM.GT.BIGNUM ) THEN
        ISCL = 1
        CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, IERR )
      END IF
!
      IF( M.GE.N ) THEN
!
!        A has at least as many rows as columns. If A has sufficiently
!        more rows than columns, first reduce using the QR
!        decomposition (if sufficient workspace available)
!
        IF( M.GE.MNTHR ) THEN
!
          IF( WNTQN ) THEN
!
!              Path 1 (M much larger than N, JOBZ='N')
!              No singular vectors to be computed
!
            ITAU = 1
            NWORK = ITAU + N
!
!              Compute A=Q*R
!              (Workspace: need 2*N, prefer N+N*NB)
!
            CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ),     &
     &LWORK-NWORK+1, IERR )
!
!              Zero out below R
!
            CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA )
            IE = 1
            ITAUQ = IE + N
            ITAUP = ITAUQ + N
            NWORK = ITAUP + N
!
!              Bidiagonalize R in A
!              (Workspace: need 4*N, prefer 3*N+2*N*NB)
!
            CALL DGEBRD( N, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ),    &
     &WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1,                      &
     &IERR )
            NWORK = IE + N
!
!              Perform bidiagonal SVD, computing singular values only
!              (Workspace: need N+BDSPAC)
!
            CALL DBDSDC( 'U', 'N', N, S, WORK( IE ), DUM, 1, DUM, 1,    &
     &DUM, IDUM, WORK( NWORK ), IWORK, INFO )
!
          ELSE IF( WNTQO ) THEN
!
!              Path 2 (M much larger than N, JOBZ = 'O')
!              N left singular vectors to be overwritten on A and
!              N right singular vectors to be computed in VT
!
            IR = 1
!
!              WORK(IR) is LDWRKR by N
!
            IF( LWORK.GE.LDA*N+N*N+3*N+BDSPAC ) THEN
              LDWRKR = LDA
            ELSE
              LDWRKR = ( LWORK-N*N-3*N-BDSPAC ) / N
            END IF
            ITAU = IR + LDWRKR*N
            NWORK = ITAU + N
!
!              Compute A=Q*R
!              (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
!
            CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ),     &
     &LWORK-NWORK+1, IERR )
!
!              Copy R to WORK(IR), zeroing out below it
!
            CALL DLACPY( 'U', N, N, A, LDA, WORK( IR ), LDWRKR )
            CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, WORK( IR+1 ),       &
     &LDWRKR )
!
!              Generate Q in A
!              (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
!
            CALL DORGQR( M, N, N, A, LDA, WORK( ITAU ),                 &
     &WORK( NWORK ), LWORK-NWORK+1, IERR )
            IE = ITAU
            ITAUQ = IE + N
            ITAUP = ITAUQ + N
            NWORK = ITAUP + N
!
!              Bidiagonalize R in VT, copying result to WORK(IR)
!              (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
!
            CALL DGEBRD( N, N, WORK( IR ), LDWRKR, S, WORK( IE ),       &
     &WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ),                      &
     &LWORK-NWORK+1, IERR )
!
!              WORK(IU) is N by N
!
            IU = NWORK
            NWORK = IU + N*N
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in WORK(IU) and computing right
!              singular vectors of bidiagonal matrix in VT
!              (Workspace: need N+N*N+BDSPAC)
!
            CALL DBDSDC( 'U', 'I', N, S, WORK( IE ), WORK( IU ), N,     &
     &VT, LDVT, DUM, IDUM, WORK( NWORK ), IWORK,                        &
     &INFO )
!
!              Overwrite WORK(IU) by left singular vectors of R
!              and VT by right singular vectors of R
!              (Workspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB)
!
            CALL DORMBR( 'Q', 'L', 'N', N, N, N, WORK( IR ), LDWRKR,    &
     &WORK( ITAUQ ), WORK( IU ), N, WORK( NWORK ),                      &
     &LWORK-NWORK+1, IERR )
            CALL DORMBR( 'P', 'R', 'T', N, N, N, WORK( IR ), LDWRKR,    &
     &WORK( ITAUP ), VT, LDVT, WORK( NWORK ),                           &
     &LWORK-NWORK+1, IERR )
!
!              Multiply Q in A by left singular vectors of R in
!              WORK(IU), storing result in WORK(IR) and copying to A
!              (Workspace: need 2*N*N, prefer N*N+M*N)
!
            DO I = 1, M, LDWRKR
              CHUNK = MIN( M-I+1, LDWRKR )
              CALL DGEMM( 'N', 'N', CHUNK, N, N, ONE, A( I, 1 ),        &
     &LDA, WORK( IU ), N, ZERO, WORK( IR ),                             &
     &LDWRKR )
              CALL DLACPY( 'F', CHUNK, N, WORK( IR ), LDWRKR,           &
     &A( I, 1 ), LDA )
            enddo
!
          ELSE IF( WNTQS ) THEN
!
!              Path 3 (M much larger than N, JOBZ='S')
!              N left singular vectors to be computed in U and
!              N right singular vectors to be computed in VT
!
            IR = 1
!
!              WORK(IR) is N by N
!
            LDWRKR = N
            ITAU = IR + LDWRKR*N
            NWORK = ITAU + N
!
!              Compute A=Q*R
!              (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
!
            CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ),     &
     &LWORK-NWORK+1, IERR )
!
!              Copy R to WORK(IR), zeroing out below it
!
            CALL DLACPY( 'U', N, N, A, LDA, WORK( IR ), LDWRKR )
            CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, WORK( IR+1 ),       &
     &LDWRKR )
!
!              Generate Q in A
!              (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
!
            CALL DORGQR( M, N, N, A, LDA, WORK( ITAU ),                 &
     &WORK( NWORK ), LWORK-NWORK+1, IERR )
            IE = ITAU
            ITAUQ = IE + N
            ITAUP = ITAUQ + N
            NWORK = ITAUP + N
!
!              Bidiagonalize R in WORK(IR)
!              (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
!
            CALL DGEBRD( N, N, WORK( IR ), LDWRKR, S, WORK( IE ),       &
     &WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ),                      &
     &LWORK-NWORK+1, IERR )
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagoal matrix in U and computing right singular
!              vectors of bidiagonal matrix in VT
!              (Workspace: need N+BDSPAC)
!
            CALL DBDSDC( 'U', 'I', N, S, WORK( IE ), U, LDU, VT,        &
     &LDVT, DUM, IDUM, WORK( NWORK ), IWORK,                            &
     &INFO )
!
!              Overwrite U by left singular vectors of R and VT
!              by right singular vectors of R
!              (Workspace: need N*N+3*N, prefer N*N+2*N+N*NB)
!
            CALL DORMBR( 'Q', 'L', 'N', N, N, N, WORK( IR ), LDWRKR,    &
     &WORK( ITAUQ ), U, LDU, WORK( NWORK ),                             &
     &LWORK-NWORK+1, IERR )
!
            CALL DORMBR( 'P', 'R', 'T', N, N, N, WORK( IR ), LDWRKR,    &
     &WORK( ITAUP ), VT, LDVT, WORK( NWORK ),                           &
     &LWORK-NWORK+1, IERR )
!
!              Multiply Q in A by left singular vectors of R in
!              WORK(IR), storing result in U
!              (Workspace: need N*N)
!
            CALL DLACPY( 'F', N, N, U, LDU, WORK( IR ), LDWRKR )
            CALL DGEMM( 'N', 'N', M, N, N, ONE, A, LDA, WORK( IR ),     &
     &LDWRKR, ZERO, U, LDU )
!
          ELSE IF( WNTQA ) THEN
!
!              Path 4 (M much larger than N, JOBZ='A')
!              M left singular vectors to be computed in U and
!              N right singular vectors to be computed in VT
!
            IU = 1
!
!              WORK(IU) is N by N
!
            LDWRKU = N
            ITAU = IU + LDWRKU*N
            NWORK = ITAU + N
!
!              Compute A=Q*R, copying result to U
!              (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
!
            CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ),     &
     &LWORK-NWORK+1, IERR )
            CALL DLACPY( 'L', M, N, A, LDA, U, LDU )
!
!              Generate Q in U
!              (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
            CALL DORGQR( M, M, N, U, LDU, WORK( ITAU ),                 &
     &WORK( NWORK ), LWORK-NWORK+1, IERR )
!
!              Produce R in A, zeroing out other entries
!
            CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA )
            IE = ITAU
            ITAUQ = IE + N
            ITAUP = ITAUQ + N
            NWORK = ITAUP + N
!
!              Bidiagonalize R in A
!              (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
!
            CALL DGEBRD( N, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ),    &
     &WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1,                      &
     &IERR )
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in WORK(IU) and computing right
!              singular vectors of bidiagonal matrix in VT
!              (Workspace: need N+N*N+BDSPAC)
!
            CALL DBDSDC( 'U', 'I', N, S, WORK( IE ), WORK( IU ), N,     &
     &VT, LDVT, DUM, IDUM, WORK( NWORK ), IWORK,                        &
     &INFO )
!
!              Overwrite WORK(IU) by left singular vectors of R and VT
!              by right singular vectors of R
!              (Workspace: need N*N+3*N, prefer N*N+2*N+N*NB)
!
            CALL DORMBR( 'Q', 'L', 'N', N, N, N, A, LDA,                &
     &WORK( ITAUQ ), WORK( IU ), LDWRKU,                                &
     &WORK( NWORK ), LWORK-NWORK+1, IERR )
            CALL DORMBR( 'P', 'R', 'T', N, N, N, A, LDA,                &
     &WORK( ITAUP ), VT, LDVT, WORK( NWORK ),                           &
     &LWORK-NWORK+1, IERR )
!
!              Multiply Q in U by left singular vectors of R in
!              WORK(IU), storing result in A
!              (Workspace: need N*N)
!
            CALL DGEMM( 'N', 'N', M, N, N, ONE, U, LDU, WORK( IU ),     &
     &LDWRKU, ZERO, A, LDA )
!
!              Copy left singular vectors of A from A to U
!
            CALL DLACPY( 'F', M, N, A, LDA, U, LDU )
!
          END IF
!
        ELSE
!
!           M .LT. MNTHR
!
!           Path 5 (M at least N, but not much larger)
!           Reduce to bidiagonal form without QR decomposition
!
          IE = 1
          ITAUQ = IE + N
          ITAUP = ITAUQ + N
          NWORK = ITAUP + N
!
!           Bidiagonalize A
!           (Workspace: need 3*N+M, prefer 3*N+(M+N)*NB)
!
          CALL DGEBRD( M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ),      &
     &WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1,                      &
     &IERR )
          IF( WNTQN ) THEN
!
!              Perform bidiagonal SVD, only computing singular values
!              (Workspace: need N+BDSPAC)
!
            CALL DBDSDC( 'U', 'N', N, S, WORK( IE ), DUM, 1, DUM, 1,    &
     &DUM, IDUM, WORK( NWORK ), IWORK, INFO )
          ELSE IF( WNTQO ) THEN
            IU = NWORK
            IF( LWORK.GE.M*N+3*N+BDSPAC ) THEN
!
!                 WORK( IU ) is M by N
!
              LDWRKU = M
              NWORK = IU + LDWRKU*N
              CALL DLASET( 'F', M, N, ZERO, ZERO, WORK( IU ),           &
     &LDWRKU )
            ELSE
!
!                 WORK( IU ) is N by N
!
              LDWRKU = N
              NWORK = IU + LDWRKU*N
!
!                 WORK(IR) is LDWRKR by N
!
              IR = NWORK
              LDWRKR = ( LWORK-N*N-3*N ) / N
            END IF
            NWORK = IU + LDWRKU*N
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in WORK(IU) and computing right
!              singular vectors of bidiagonal matrix in VT
!              (Workspace: need N+N*N+BDSPAC)
!
            CALL DBDSDC( 'U', 'I', N, S, WORK( IE ), WORK( IU ),        &
     &LDWRKU, VT, LDVT, DUM, IDUM, WORK( NWORK ),                       &
     &IWORK, INFO )
!
!              Overwrite VT by right singular vectors of A
!              (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
!
            CALL DORMBR( 'P', 'R', 'T', N, N, N, A, LDA,                &
     &WORK( ITAUP ), VT, LDVT, WORK( NWORK ),                           &
     &LWORK-NWORK+1, IERR )
!
            IF( LWORK.GE.M*N+3*N+BDSPAC ) THEN
!
!                 Overwrite WORK(IU) by left singular vectors of A
!                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
!
              CALL DORMBR( 'Q', 'L', 'N', M, N, N, A, LDA,              &
     &WORK( ITAUQ ), WORK( IU ), LDWRKU,                                &
     &WORK( NWORK ), LWORK-NWORK+1, IERR )
!
!                 Copy left singular vectors of A from WORK(IU) to A
!
              CALL DLACPY( 'F', M, N, WORK( IU ), LDWRKU, A, LDA )
            ELSE
!
!                 Generate Q in A
!                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
!
              CALL DORGBR( 'Q', M, N, N, A, LDA, WORK( ITAUQ ),         &
     &WORK( NWORK ), LWORK-NWORK+1, IERR )
!
!                 Multiply Q in A by left singular vectors of
!                 bidiagonal matrix in WORK(IU), storing result in
!                 WORK(IR) and copying to A
!                 (Workspace: need 2*N*N, prefer N*N+M*N)
!
              DO I = 1, M, LDWRKR
                CHUNK = MIN( M-I+1, LDWRKR )
                CALL DGEMM( 'N', 'N', CHUNK, N, N, ONE, A( I, 1 ),      &
     &LDA, WORK( IU ), LDWRKU, ZERO,                                    &
     &WORK( IR ), LDWRKR )
                CALL DLACPY( 'F', CHUNK, N, WORK( IR ), LDWRKR,         &
     &A( I, 1 ), LDA )
              enddo
            END IF
!
          ELSE IF( WNTQS ) THEN
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in U and computing right singular
!              vectors of bidiagonal matrix in VT
!              (Workspace: need N+BDSPAC)
!
            CALL DLASET( 'F', M, N, ZERO, ZERO, U, LDU )
            CALL DBDSDC( 'U', 'I', N, S, WORK( IE ), U, LDU, VT,        &
     &LDVT, DUM, IDUM, WORK( NWORK ), IWORK,                            &
     &INFO )
!
!              Overwrite U by left singular vectors of A and VT
!              by right singular vectors of A
!              (Workspace: need 3*N, prefer 2*N+N*NB)
!
            CALL DORMBR( 'Q', 'L', 'N', M, N, N, A, LDA,                &
     &WORK( ITAUQ ), U, LDU, WORK( NWORK ),                             &
     &LWORK-NWORK+1, IERR )
            CALL DORMBR( 'P', 'R', 'T', N, N, N, A, LDA,                &
     &WORK( ITAUP ), VT, LDVT, WORK( NWORK ),                           &
     &LWORK-NWORK+1, IERR )
          ELSE IF( WNTQA ) THEN
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in U and computing right singular
!              vectors of bidiagonal matrix in VT
!              (Workspace: need N+BDSPAC)
!
            CALL DLASET( 'F', M, M, ZERO, ZERO, U, LDU )
            CALL DBDSDC( 'U', 'I', N, S, WORK( IE ), U, LDU, VT,        &
     &LDVT, DUM, IDUM, WORK( NWORK ), IWORK,                            &
     &INFO )
!
!              Set the right corner of U to identity matrix
!
            IF( M.GT.N ) THEN
              CALL DLASET( 'F', M-N, M-N, ZERO, ONE, U( N+1, N+1 ),     &
     &LDU )
            END IF
!
!              Overwrite U by left singular vectors of A and VT
!              by right singular vectors of A
!              (Workspace: need N*N+2*N+M, prefer N*N+2*N+M*NB)
!
            CALL DORMBR( 'Q', 'L', 'N', M, M, N, A, LDA,                &
     &WORK( ITAUQ ), U, LDU, WORK( NWORK ),                             &
     &LWORK-NWORK+1, IERR )
            CALL DORMBR( 'P', 'R', 'T', N, N, M, A, LDA,                &
     &WORK( ITAUP ), VT, LDVT, WORK( NWORK ),                           &
     &LWORK-NWORK+1, IERR )
          END IF
!
        END IF
!
      ELSE
!
!        A has more columns than rows. If A has sufficiently more
!        columns than rows, first reduce using the LQ decomposition (if
!        sufficient workspace available)
!
        IF( N.GE.MNTHR ) THEN
!
          IF( WNTQN ) THEN
!
!              Path 1t (N much larger than M, JOBZ='N')
!              No singular vectors to be computed
!
            ITAU = 1
            NWORK = ITAU + M
!
!              Compute A=L*Q
!              (Workspace: need 2*M, prefer M+M*NB)
!
            CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ),     &
     &LWORK-NWORK+1, IERR )
!
!              Zero out above L
!
            CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, A( 1, 2 ), LDA )
            IE = 1
            ITAUQ = IE + M
            ITAUP = ITAUQ + M
            NWORK = ITAUP + M
!
!              Bidiagonalize L in A
!              (Workspace: need 4*M, prefer 3*M+2*M*NB)
!
            CALL DGEBRD( M, M, A, LDA, S, WORK( IE ), WORK( ITAUQ ),    &
     &WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1,                      &
     &IERR )
            NWORK = IE + M
!
!              Perform bidiagonal SVD, computing singular values only
!              (Workspace: need M+BDSPAC)
!
            CALL DBDSDC( 'U', 'N', M, S, WORK( IE ), DUM, 1, DUM, 1,    &
     &DUM, IDUM, WORK( NWORK ), IWORK, INFO )
!
          ELSE IF( WNTQO ) THEN
!
!              Path 2t (N much larger than M, JOBZ='O')
!              M right singular vectors to be overwritten on A and
!              M left singular vectors to be computed in U
!
            IVT = 1
!
!              IVT is M by M
!
            IL = IVT + M*M
            IF( LWORK.GE.M*N+M*M+3*M+BDSPAC ) THEN
!
!                 WORK(IL) is M by N
!
              LDWRKL = M
              CHUNK = N
            ELSE
              LDWRKL = M
              CHUNK = ( LWORK-M*M ) / M
            END IF
            ITAU = IL + LDWRKL*M
            NWORK = ITAU + M
!
!              Compute A=L*Q
!              (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
!
            CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ),     &
     &LWORK-NWORK+1, IERR )
!
!              Copy L to WORK(IL), zeroing about above it
!
            CALL DLACPY( 'L', M, M, A, LDA, WORK( IL ), LDWRKL )
            CALL DLASET( 'U', M-1, M-1, ZERO, ZERO,                     &
     &WORK( IL+LDWRKL ), LDWRKL )
!
!              Generate Q in A
!              (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
!
            CALL DORGLQ( M, N, M, A, LDA, WORK( ITAU ),                 &
     &WORK( NWORK ), LWORK-NWORK+1, IERR )
            IE = ITAU
            ITAUQ = IE + M
            ITAUP = ITAUQ + M
            NWORK = ITAUP + M
!
!              Bidiagonalize L in WORK(IL)
!              (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
!
            CALL DGEBRD( M, M, WORK( IL ), LDWRKL, S, WORK( IE ),       &
     &WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ),                      &
     &LWORK-NWORK+1, IERR )
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in U, and computing right singular
!              vectors of bidiagonal matrix in WORK(IVT)
!              (Workspace: need M+M*M+BDSPAC)
!
            CALL DBDSDC( 'U', 'I', M, S, WORK( IE ), U, LDU,            &
     &WORK( IVT ), M, DUM, IDUM, WORK( NWORK ),                         &
     &IWORK, INFO )
!
!              Overwrite U by left singular vectors of L and WORK(IVT)
!              by right singular vectors of L
!              (Workspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB)
!
            CALL DORMBR( 'Q', 'L', 'N', M, M, M, WORK( IL ), LDWRKL,    &
     &WORK( ITAUQ ), U, LDU, WORK( NWORK ),                             &
     &LWORK-NWORK+1, IERR )
            CALL DORMBR( 'P', 'R', 'T', M, M, M, WORK( IL ), LDWRKL,    &
     &WORK( ITAUP ), WORK( IVT ), M,                                    &
     &WORK( NWORK ), LWORK-NWORK+1, IERR )
!
!              Multiply right singular vectors of L in WORK(IVT) by Q
!              in A, storing result in WORK(IL) and copying to A
!              (Workspace: need 2*M*M, prefer M*M+M*N)
!
            DO I = 1, N, CHUNK
              BLK = MIN( N-I+1, CHUNK )
              CALL DGEMM( 'N', 'N', M, BLK, M, ONE, WORK( IVT ), M,     &
     &A( 1, I ), LDA, ZERO, WORK( IL ), LDWRKL )
              CALL DLACPY( 'F', M, BLK, WORK( IL ), LDWRKL,             &
     &A( 1, I ), LDA )
            enddo
!
          ELSE IF( WNTQS ) THEN
!
!              Path 3t (N much larger than M, JOBZ='S')
!              M right singular vectors to be computed in VT and
!              M left singular vectors to be computed in U
!
            IL = 1
!
!              WORK(IL) is M by M
!
            LDWRKL = M
            ITAU = IL + LDWRKL*M
            NWORK = ITAU + M
!
!              Compute A=L*Q
!              (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
!
            CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ),     &
     &LWORK-NWORK+1, IERR )
!
!              Copy L to WORK(IL), zeroing out above it
!
            CALL DLACPY( 'L', M, M, A, LDA, WORK( IL ), LDWRKL )
            CALL DLASET( 'U', M-1, M-1, ZERO, ZERO,                     &
     &WORK( IL+LDWRKL ), LDWRKL )
!
!              Generate Q in A
!              (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
!
            CALL DORGLQ( M, N, M, A, LDA, WORK( ITAU ),                 &
     &WORK( NWORK ), LWORK-NWORK+1, IERR )
            IE = ITAU
            ITAUQ = IE + M
            ITAUP = ITAUQ + M
            NWORK = ITAUP + M
!
!              Bidiagonalize L in WORK(IU), copying result to U
!              (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
!
            CALL DGEBRD( M, M, WORK( IL ), LDWRKL, S, WORK( IE ),       &
     &WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ),                      &
     &LWORK-NWORK+1, IERR )
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in U and computing right singular
!              vectors of bidiagonal matrix in VT
!              (Workspace: need M+BDSPAC)
!
            CALL DBDSDC( 'U', 'I', M, S, WORK( IE ), U, LDU, VT,        &
     &LDVT, DUM, IDUM, WORK( NWORK ), IWORK,                            &
     &INFO )
!
!              Overwrite U by left singular vectors of L and VT
!              by right singular vectors of L
!              (Workspace: need M*M+3*M, prefer M*M+2*M+M*NB)
!
            CALL DORMBR( 'Q', 'L', 'N', M, M, M, WORK( IL ), LDWRKL,    &
     &WORK( ITAUQ ), U, LDU, WORK( NWORK ),                             &
     &LWORK-NWORK+1, IERR )
            CALL DORMBR( 'P', 'R', 'T', M, M, M, WORK( IL ), LDWRKL,    &
     &WORK( ITAUP ), VT, LDVT, WORK( NWORK ),                           &
     &LWORK-NWORK+1, IERR )
!
!              Multiply right singular vectors of L in WORK(IL) by
!              Q in A, storing result in VT
!              (Workspace: need M*M)
!
            CALL DLACPY( 'F', M, M, VT, LDVT, WORK( IL ), LDWRKL )
            CALL DGEMM( 'N', 'N', M, N, M, ONE, WORK( IL ), LDWRKL,     &
     &A, LDA, ZERO, VT, LDVT )
!
          ELSE IF( WNTQA ) THEN
!
!              Path 4t (N much larger than M, JOBZ='A')
!              N right singular vectors to be computed in VT and
!              M left singular vectors to be computed in U
!
            IVT = 1
!
!              WORK(IVT) is M by M
!
            LDWKVT = M
            ITAU = IVT + LDWKVT*M
            NWORK = ITAU + M
!
!              Compute A=L*Q, copying result to VT
!              (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
!
            CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ),     &
     &LWORK-NWORK+1, IERR )
            CALL DLACPY( 'U', M, N, A, LDA, VT, LDVT )
!
!              Generate Q in VT
!              (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
!
            CALL DORGLQ( N, N, M, VT, LDVT, WORK( ITAU ),               &
     &WORK( NWORK ), LWORK-NWORK+1, IERR )
!
!              Produce L in A, zeroing out other entries
!
            CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, A( 1, 2 ), LDA )
            IE = ITAU
            ITAUQ = IE + M
            ITAUP = ITAUQ + M
            NWORK = ITAUP + M
!
!              Bidiagonalize L in A
!              (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
!
            CALL DGEBRD( M, M, A, LDA, S, WORK( IE ), WORK( ITAUQ ),    &
     &WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1,                      &
     &IERR )
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in U and computing right singular
!              vectors of bidiagonal matrix in WORK(IVT)
!              (Workspace: need M+M*M+BDSPAC)
!
            CALL DBDSDC( 'U', 'I', M, S, WORK( IE ), U, LDU,            &
     &WORK( IVT ), LDWKVT, DUM, IDUM,                                   &
     &WORK( NWORK ), IWORK, INFO )
!
!              Overwrite U by left singular vectors of L and WORK(IVT)
!              by right singular vectors of L
!              (Workspace: need M*M+3*M, prefer M*M+2*M+M*NB)
!
            CALL DORMBR( 'Q', 'L', 'N', M, M, M, A, LDA,                &
     &WORK( ITAUQ ), U, LDU, WORK( NWORK ),                             &
     &LWORK-NWORK+1, IERR )
            CALL DORMBR( 'P', 'R', 'T', M, M, M, A, LDA,                &
     &WORK( ITAUP ), WORK( IVT ), LDWKVT,                               &
     &WORK( NWORK ), LWORK-NWORK+1, IERR )
!
!              Multiply right singular vectors of L in WORK(IVT) by
!              Q in VT, storing result in A
!              (Workspace: need M*M)
!
            CALL DGEMM( 'N', 'N', M, N, M, ONE, WORK( IVT ), LDWKVT,    &
     &VT, LDVT, ZERO, A, LDA )
!
!              Copy right singular vectors of A from A to VT
!
            CALL DLACPY( 'F', M, N, A, LDA, VT, LDVT )
!
          END IF
!
        ELSE
!
!           N .LT. MNTHR
!
!           Path 5t (N greater than M, but not much larger)
!           Reduce to bidiagonal form without LQ decomposition
!
          IE = 1
          ITAUQ = IE + M
          ITAUP = ITAUQ + M
          NWORK = ITAUP + M
!
!           Bidiagonalize A
!           (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)
!
          CALL DGEBRD( M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ),      &
     &WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1,                      &
     &IERR )
          IF( WNTQN ) THEN
!
!              Perform bidiagonal SVD, only computing singular values
!              (Workspace: need M+BDSPAC)
!
            CALL DBDSDC( 'L', 'N', M, S, WORK( IE ), DUM, 1, DUM, 1,    &
     &DUM, IDUM, WORK( NWORK ), IWORK, INFO )
          ELSE IF( WNTQO ) THEN
            LDWKVT = M
            IVT = NWORK
            IF( LWORK.GE.M*N+3*M+BDSPAC ) THEN
!
!                 WORK( IVT ) is M by N
!
              CALL DLASET( 'F', M, N, ZERO, ZERO, WORK( IVT ),          &
     &LDWKVT )
              NWORK = IVT + LDWKVT*N
            ELSE
!
!                 WORK( IVT ) is M by M
!
              NWORK = IVT + LDWKVT*M
              IL = NWORK
!
!                 WORK(IL) is M by CHUNK
!
              CHUNK = ( LWORK-M*M-3*M ) / M
            END IF
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in U and computing right singular
!              vectors of bidiagonal matrix in WORK(IVT)
!              (Workspace: need M*M+BDSPAC)
!
            CALL DBDSDC( 'L', 'I', M, S, WORK( IE ), U, LDU,            &
     &WORK( IVT ), LDWKVT, DUM, IDUM,                                   &
     &WORK( NWORK ), IWORK, INFO )
!
!              Overwrite U by left singular vectors of A
!              (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
!
            CALL DORMBR( 'Q', 'L', 'N', M, M, N, A, LDA,                &
     &WORK( ITAUQ ), U, LDU, WORK( NWORK ),                             &
     &LWORK-NWORK+1, IERR )
!
            IF( LWORK.GE.M*N+3*M+BDSPAC ) THEN
!
!                 Overwrite WORK(IVT) by left singular vectors of A
!                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
!
              CALL DORMBR( 'P', 'R', 'T', M, N, M, A, LDA,              &
     &WORK( ITAUP ), WORK( IVT ), LDWKVT,                               &
     &WORK( NWORK ), LWORK-NWORK+1, IERR )
!
!                 Copy right singular vectors of A from WORK(IVT) to A
!
              CALL DLACPY( 'F', M, N, WORK( IVT ), LDWKVT, A, LDA )
            ELSE
!
!                 Generate P**T in A
!                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
!
              CALL DORGBR( 'P', M, N, M, A, LDA, WORK( ITAUP ),         &
     &WORK( NWORK ), LWORK-NWORK+1, IERR )
!
!                 Multiply Q in A by right singular vectors of
!                 bidiagonal matrix in WORK(IVT), storing result in
!                 WORK(IL) and copying to A
!                 (Workspace: need 2*M*M, prefer M*M+M*N)
!
              DO I = 1, N, CHUNK
                BLK = MIN( N-I+1, CHUNK )
                CALL DGEMM( 'N', 'N', M, BLK, M, ONE, WORK( IVT ),      &
     &LDWKVT, A( 1, I ), LDA, ZERO,                                     &
     &WORK( IL ), M )
                CALL DLACPY( 'F', M, BLK, WORK( IL ), M, A( 1, I ),     &
     &LDA )
              enddo
            END IF
          ELSE IF( WNTQS ) THEN
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in U and computing right singular
!              vectors of bidiagonal matrix in VT
!              (Workspace: need M+BDSPAC)
!
            CALL DLASET( 'F', M, N, ZERO, ZERO, VT, LDVT )
            CALL DBDSDC( 'L', 'I', M, S, WORK( IE ), U, LDU, VT,        &
     &LDVT, DUM, IDUM, WORK( NWORK ), IWORK,                            &
     &INFO )
!
!              Overwrite U by left singular vectors of A and VT
!              by right singular vectors of A
!              (Workspace: need 3*M, prefer 2*M+M*NB)
!
            CALL DORMBR( 'Q', 'L', 'N', M, M, N, A, LDA,                &
     &WORK( ITAUQ ), U, LDU, WORK( NWORK ),                             &
     &LWORK-NWORK+1, IERR )
            CALL DORMBR( 'P', 'R', 'T', M, N, M, A, LDA,                &
     &WORK( ITAUP ), VT, LDVT, WORK( NWORK ),                           &
     &LWORK-NWORK+1, IERR )
          ELSE IF( WNTQA ) THEN
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in U and computing right singular
!              vectors of bidiagonal matrix in VT
!              (Workspace: need M+BDSPAC)
!
            CALL DLASET( 'F', N, N, ZERO, ZERO, VT, LDVT )
            CALL DBDSDC( 'L', 'I', M, S, WORK( IE ), U, LDU, VT,        &
     &LDVT, DUM, IDUM, WORK( NWORK ), IWORK,                            &
     &INFO )
!
!              Set the right corner of VT to identity matrix
!
            IF( N.GT.M ) THEN
              CALL DLASET( 'F', N-M, N-M, ZERO, ONE, VT( M+1, M+1 ),    &
     &LDVT )
            END IF
!
!              Overwrite U by left singular vectors of A and VT
!              by right singular vectors of A
!              (Workspace: need 2*M+N, prefer 2*M+N*NB)
!
            CALL DORMBR( 'Q', 'L', 'N', M, M, N, A, LDA,                &
     &WORK( ITAUQ ), U, LDU, WORK( NWORK ),                             &
     &LWORK-NWORK+1, IERR )
            CALL DORMBR( 'P', 'R', 'T', N, N, M, A, LDA,                &
     &WORK( ITAUP ), VT, LDVT, WORK( NWORK ),                           &
     &LWORK-NWORK+1, IERR )
          END IF
!
        END IF
!
      END IF
!
!     Undo scaling if necessary
!
      IF( ISCL.EQ.1 ) THEN
        IF( ANRM.GT.BIGNUM )                                            &
     &CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN,         &
     &IERR )
        IF( ANRM.LT.SMLNUM )                                            &
     &CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN,         &
     &IERR )
      END IF
!
!     Return optimal workspace in WORK(1)
!
      WORK( 1 ) = MAXWRK
!
      RETURN
!
!     End of DGESDD
!
      END
      SUBROUTINE DLABAD( SMALL, LARGE )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   LARGE, SMALL
!     ..
!
!  Purpose
!  =======
!
!  DLABAD takes as input the values computed by DLAMCH for underflow and
!  overflow, and returns the square root of each of these values if the
!  log of LARGE is sufficiently large.  This subroutine is intended to
!  identify machines with a large exponent range, such as the Crays, and
!  redefine the underflow and overflow limits to be the square roots of
!  the values computed by DLAMCH.  This subroutine is needed because
!  DLAMCH does not compensate for poor arithmetic in the upper half of
!  the exponent range, as is found on a Cray.
!
!  Arguments
!  =========
!
!  SMALL   (input/output) DOUBLE PRECISION
!          On entry, the underflow threshold as computed by DLAMCH.
!          On exit, if LOG10(LARGE) is sufficiently large, the square
!          root of SMALL, otherwise unchanged.
!
!  LARGE   (input/output) DOUBLE PRECISION
!          On entry, the overflow threshold as computed by DLAMCH.
!          On exit, if LOG10(LARGE) is sufficiently large, the square
!          root of LARGE, otherwise unchanged.
!
!  =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          LOG10, SQRT
!     ..
!     .. Executable Statements ..
!
!     If it looks like we're on a Cray, take the square root of
!     SMALL and LARGE to avoid overflow and underflow problems.
!
      IF( LOG10( LARGE ).GT.2000.D0 ) THEN
        SMALL = SQRT( SMALL )
        LARGE = SQRT( LARGE )
      END IF
!
      RETURN
!
!     End of DLABAD
!
      END
      SUBROUTINE DLABRD( M, N, NB, A, LDA, D, E, TAUQ, TAUP, X, LDX, Y, &
     &LDY )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            LDA, LDX, LDY, M, N, NB
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAUP( * ),        &
     &TAUQ( * ), X( LDX, * ), Y( LDY, * )
!     ..
!
!  Purpose
!  =======
!
!  DLABRD reduces the first NB rows and columns of a real general
!  m by n matrix A to upper or lower bidiagonal form by an orthogonal
!  transformation Q' * A * P, and returns the matrices X and Y which
!  are needed to apply the transformation to the unreduced part of A.
!
!  If m >= n, A is reduced to upper bidiagonal form; if m < n, to lower
!  bidiagonal form.
!
!  This is an auxiliary routine called by DGEBRD
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows in the matrix A.
!
!  N       (input) INTEGER
!          The number of columns in the matrix A.
!
!  NB      (input) INTEGER
!          The number of leading rows and columns of A to be reduced.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the m by n general matrix to be reduced.
!          On exit, the first NB rows and columns of the matrix are
!          overwritten; the rest of the array is unchanged.
!          If m >= n, elements on and below the diagonal in the first NB
!            columns, with the array TAUQ, represent the orthogonal
!            matrix Q as a product of elementary reflectors; and
!            elements above the diagonal in the first NB rows, with the
!            array TAUP, represent the orthogonal matrix P as a product
!            of elementary reflectors.
!          If m < n, elements below the diagonal in the first NB
!            columns, with the array TAUQ, represent the orthogonal
!            matrix Q as a product of elementary reflectors, and
!            elements on and above the diagonal in the first NB rows,
!            with the array TAUP, represent the orthogonal matrix P as
!            a product of elementary reflectors.
!          See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  D       (output) DOUBLE PRECISION array, dimension (NB)
!          The diagonal elements of the first NB rows and columns of
!          the reduced matrix.  D(i) = A(i,i).
!
!  E       (output) DOUBLE PRECISION array, dimension (NB)
!          The off-diagonal elements of the first NB rows and columns of
!          the reduced matrix.
!
!  TAUQ    (output) DOUBLE PRECISION array dimension (NB)
!          The scalar factors of the elementary reflectors which
!          represent the orthogonal matrix Q. See Further Details.
!
!  TAUP    (output) DOUBLE PRECISION array, dimension (NB)
!          The scalar factors of the elementary reflectors which
!          represent the orthogonal matrix P. See Further Details.
!
!  X       (output) DOUBLE PRECISION array, dimension (LDX,NB)
!          The m-by-nb matrix X required to update the unreduced part
!          of A.
!
!  LDX     (input) INTEGER
!          The leading dimension of the array X. LDX >= M.
!
!  Y       (output) DOUBLE PRECISION array, dimension (LDY,NB)
!          The n-by-nb matrix Y required to update the unreduced part
!          of A.
!
!  LDY     (input) INTEGER
!          The leading dimension of the array Y. LDY >= N.
!
!  Further Details
!  ===============
!
!  The matrices Q and P are represented as products of elementary
!  reflectors:
!
!     Q = H(1) H(2) . . . H(nb)  and  P = G(1) G(2) . . . G(nb)
!
!  Each H(i) and G(i) has the form:
!
!     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
!
!  where tauq and taup are real scalars, and v and u are real vectors.
!
!  If m >= n, v(1:i-1) = 0, v(i) = 1, and v(i:m) is stored on exit in
!  A(i:m,i); u(1:i) = 0, u(i+1) = 1, and u(i+1:n) is stored on exit in
!  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).
!
!  If m < n, v(1:i) = 0, v(i+1) = 1, and v(i+1:m) is stored on exit in
!  A(i+2:m,i); u(1:i-1) = 0, u(i) = 1, and u(i:n) is stored on exit in
!  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).
!
!  The elements of the vectors v and u together form the m-by-nb matrix
!  V and the nb-by-n matrix U' which are needed, with X and Y, to apply
!  the transformation to the unreduced part of the matrix, using a block
!  update of the form:  A := A - V*Y' - X*U'.
!
!  The contents of A on exit are illustrated by the following examples
!  with nb = 2:
!
!  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
!
!    (  1   1   u1  u1  u1 )           (  1   u1  u1  u1  u1  u1 )
!    (  v1  1   1   u2  u2 )           (  1   1   u2  u2  u2  u2 )
!    (  v1  v2  a   a   a  )           (  v1  1   a   a   a   a  )
!    (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )
!    (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )
!    (  v1  v2  a   a   a  )
!
!  where a denotes an element of the original matrix which is unchanged,
!  vi denotes an element of the vector defining H(i), and ui an element
!  of the vector defining G(i).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMV, DLARFG, DSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( M.LE.0 .OR. N.LE.0 )                                          &
     &RETURN
!
      IF( M.GE.N ) THEN
!
!        Reduce to upper bidiagonal form
!
        DO I = 1, NB
!
!           Update A(i:m,i)
!
          CALL DGEMV( 'No transpose', M-I+1, I-1, -ONE, A( I, 1 ),      &
     &LDA, Y( I, 1 ), LDY, ONE, A( I, I ), 1 )
          CALL DGEMV( 'No transpose', M-I+1, I-1, -ONE, X( I, 1 ),      &
     &LDX, A( 1, I ), 1, ONE, A( I, I ), 1 )
!
!           Generate reflection Q(i) to annihilate A(i+1:m,i)
!
          CALL DLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1,      &
     &TAUQ( I ) )
          D( I ) = A( I, I )
          IF( I.LT.N ) THEN
            A( I, I ) = ONE
!
!              Compute Y(i+1:n,i)
!
            CALL DGEMV( 'Transpose', M-I+1, N-I, ONE, A( I, I+1 ),      &
     &LDA, A( I, I ), 1, ZERO, Y( I+1, I ), 1 )
            CALL DGEMV( 'Transpose', M-I+1, I-1, ONE, A( I, 1 ), LDA,   &
     &A( I, I ), 1, ZERO, Y( 1, I ), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, -ONE, Y( I+1, 1 ),    &
     &LDY, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 )
            CALL DGEMV( 'Transpose', M-I+1, I-1, ONE, X( I, 1 ), LDX,   &
     &A( I, I ), 1, ZERO, Y( 1, I ), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, -ONE, A( 1, I+1 ),       &
     &LDA, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 )
            CALL DSCAL( N-I, TAUQ( I ), Y( I+1, I ), 1 )
!
!              Update A(i,i+1:n)
!
            CALL DGEMV( 'No transpose', N-I, I, -ONE, Y( I+1, 1 ),      &
     &LDY, A( I, 1 ), LDA, ONE, A( I, I+1 ), LDA )
            CALL DGEMV( 'Transpose', I-1, N-I, -ONE, A( 1, I+1 ),       &
     &LDA, X( I, 1 ), LDX, ONE, A( I, I+1 ), LDA )
!
!              Generate reflection P(i) to annihilate A(i,i+2:n)
!
            CALL DLARFG( N-I, A( I, I+1 ), A( I, MIN( I+2, N ) ),       &
     &LDA, TAUP( I ) )
            E( I ) = A( I, I+1 )
            A( I, I+1 ) = ONE
!
!              Compute X(i+1:m,i)
!
            CALL DGEMV( 'No transpose', M-I, N-I, ONE, A( I+1, I+1 ),   &
     &LDA, A( I, I+1 ), LDA, ZERO, X( I+1, I ), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, Y( I+1, 1 ), LDY,     &
     &A( I, I+1 ), LDA, ZERO, X( 1, I ), 1 )
            CALL DGEMV( 'No transpose', M-I, I, -ONE, A( I+1, 1 ),      &
     &LDA, X( 1, I ), 1, ONE, X( I+1, I ), 1 )
            CALL DGEMV( 'No transpose', I-1, N-I, ONE, A( 1, I+1 ),     &
     &LDA, A( I, I+1 ), LDA, ZERO, X( 1, I ), 1 )
            CALL DGEMV( 'No transpose', M-I, I-1, -ONE, X( I+1, 1 ),    &
     &LDX, X( 1, I ), 1, ONE, X( I+1, I ), 1 )
            CALL DSCAL( M-I, TAUP( I ), X( I+1, I ), 1 )
          END IF
        enddo
      ELSE
!
!        Reduce to lower bidiagonal form
!
        DO I = 1, NB
!
!           Update A(i,i:n)
!
          CALL DGEMV( 'No transpose', N-I+1, I-1, -ONE, Y( I, 1 ),      &
     &LDY, A( I, 1 ), LDA, ONE, A( I, I ), LDA )
          CALL DGEMV( 'Transpose', I-1, N-I+1, -ONE, A( 1, I ), LDA,    &
     &X( I, 1 ), LDX, ONE, A( I, I ), LDA )
!
!           Generate reflection P(i) to annihilate A(i,i+1:n)
!
          CALL DLARFG( N-I+1, A( I, I ), A( I, MIN( I+1, N ) ), LDA,    &
     &TAUP( I ) )
          D( I ) = A( I, I )
          IF( I.LT.M ) THEN
            A( I, I ) = ONE
!
!              Compute X(i+1:m,i)
!
            CALL DGEMV( 'No transpose', M-I, N-I+1, ONE, A( I+1, I ),   &
     &LDA, A( I, I ), LDA, ZERO, X( I+1, I ), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, Y( I, 1 ), LDY,   &
     &A( I, I ), LDA, ZERO, X( 1, I ), 1 )
            CALL DGEMV( 'No transpose', M-I, I-1, -ONE, A( I+1, 1 ),    &
     &LDA, X( 1, I ), 1, ONE, X( I+1, I ), 1 )
            CALL DGEMV( 'No transpose', I-1, N-I+1, ONE, A( 1, I ),     &
     &LDA, A( I, I ), LDA, ZERO, X( 1, I ), 1 )
            CALL DGEMV( 'No transpose', M-I, I-1, -ONE, X( I+1, 1 ),    &
     &LDX, X( 1, I ), 1, ONE, X( I+1, I ), 1 )
            CALL DSCAL( M-I, TAUP( I ), X( I+1, I ), 1 )
!
!              Update A(i+1:m,i)
!
            CALL DGEMV( 'No transpose', M-I, I-1, -ONE, A( I+1, 1 ),    &
     &LDA, Y( I, 1 ), LDY, ONE, A( I+1, I ), 1 )
            CALL DGEMV( 'No transpose', M-I, I, -ONE, X( I+1, 1 ),      &
     &LDX, A( 1, I ), 1, ONE, A( I+1, I ), 1 )
!
!              Generate reflection Q(i) to annihilate A(i+2:m,i)
!
            CALL DLARFG( M-I, A( I+1, I ), A( MIN( I+2, M ), I ), 1,    &
     &TAUQ( I ) )
            E( I ) = A( I+1, I )
            A( I+1, I ) = ONE
!
!              Compute Y(i+1:n,i)
!
            CALL DGEMV( 'Transpose', M-I, N-I, ONE, A( I+1, I+1 ),      &
     &LDA, A( I+1, I ), 1, ZERO, Y( I+1, I ), 1 )
            CALL DGEMV( 'Transpose', M-I, I-1, ONE, A( I+1, 1 ), LDA,   &
     &A( I+1, I ), 1, ZERO, Y( 1, I ), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, -ONE, Y( I+1, 1 ),    &
     &LDY, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 )
            CALL DGEMV( 'Transpose', M-I, I, ONE, X( I+1, 1 ), LDX,     &
     &A( I+1, I ), 1, ZERO, Y( 1, I ), 1 )
            CALL DGEMV( 'Transpose', I, N-I, -ONE, A( 1, I+1 ), LDA,    &
     &Y( 1, I ), 1, ONE, Y( I+1, I ), 1 )
            CALL DSCAL( N-I, TAUQ( I ), Y( I+1, I ), 1 )
          END IF
        enddo
      END IF
      RETURN
!
!     End of DLABRD
!
      END
      SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DLACPY copies all or part of a two-dimensional matrix A to another
!  matrix B.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies the part of the matrix A to be copied to B.
!          = 'U':      Upper triangular part
!          = 'L':      Lower triangular part
!          Otherwise:  All of the matrix A
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The m by n matrix A.  If UPLO = 'U', only the upper triangle
!          or trapezoid is accessed; if UPLO = 'L', only the lower
!          triangle or trapezoid is accessed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  B       (output) DOUBLE PRECISION array, dimension (LDB,N)
!          On exit, B = A in the locations specified by UPLO.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,M).
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( UPLO, 'U' ) ) THEN
        DO J = 1, N
          DO I = 1, MIN( J, M )
            B( I, J ) = A( I, J )
          enddo
        enddo
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
        DO J = 1, N
          DO I = J, M
            B( I, J ) = A( I, J )
          enddo
        enddo
      ELSE
        DO J = 1, N
          DO I = 1, M
            B( I, J ) = A( I, J )
          enddo
        enddo
      END IF
      RETURN
!
!     End of DLACPY
!
      END
      SUBROUTINE DLAED6( KNITER, ORGATI, RHO, D, Z, FINIT, TAU, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      LOGICAL            ORGATI
      INTEGER            INFO, KNITER
      DOUBLE PRECISION   FINIT, RHO, TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( 3 ), Z( 3 )
!     ..
!
!  Purpose
!  =======
!
!  DLAED6 computes the positive or negative root (closest to the origin)
!  of
!                   z(1)        z(2)        z(3)
!  f(x) =   rho + --------- + ---------- + ---------
!                  d(1)-x      d(2)-x      d(3)-x
!
!  It is assumed that
!
!        if ORGATI = .true. the root is between d(2) and d(3);
!        otherwise it is between d(1) and d(2)
!
!  This routine will be called by DLAED4 when necessary. In most cases,
!  the root sought is the smallest in magnitude, though it might not be
!  in some extremely rare situations.
!
!  Arguments
!  =========
!
!  KNITER       (input) INTEGER
!               Refer to DLAED4 for its significance.
!
!  ORGATI       (input) LOGICAL
!               If ORGATI is true, the needed root is between d(2) and
!               d(3); otherwise it is between d(1) and d(2).  See
!               DLAED4 for further details.
!
!  RHO          (input) DOUBLE PRECISION
!               Refer to the equation f(x) above.
!
!  D            (input) DOUBLE PRECISION array, dimension (3)
!               D satisfies d(1) < d(2) < d(3).
!
!  Z            (input) DOUBLE PRECISION array, dimension (3)
!               Each of the elements in z must be positive.
!
!  FINIT        (input) DOUBLE PRECISION
!               The value of f at 0. It is more accurate than the one
!               evaluated inside this routine (if someone wants to do
!               so).
!
!  TAU          (output) DOUBLE PRECISION
!               The root of the equation f(x).
!
!  INFO         (output) INTEGER
!               = 0: successful exit
!               > 0: if INFO = 1, failure to converge
!
!  Further Details
!  ===============
!
!  30/06/99: Based on contributions by
!     Ren-Cang Li, Computer Science Division, University of California
!     at Berkeley, USA
!
!  10/02/03: This version has a few statements commented out for thread safety
!     (machine parameters are computed on each entry). SJH.
!
!  05/10/06: Modified from a new version of Ren-Cang Li, use
!     Gragg-Thornton-Warner cubic convergent scheme for better stability.
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 40 )
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE, FOUR, EIGHT
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,      &
     &THREE = 3.0D0, FOUR = 4.0D0, EIGHT = 8.0D0 )
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   DSCALE( 3 ), ZSCALE( 3 )
!     ..
!     .. Local Scalars ..
      LOGICAL            SCALE
      INTEGER            I, ITER, NITER
      DOUBLE PRECISION   A, B, BASE, C, DDF, DF, EPS, ERRETM, ETA, F,   &
     &FC, SCLFAC, SCLINV, SMALL1, SMALL2, SMINV1,                       &
     &SMINV2, TEMP, TEMP1, TEMP2, TEMP3, TEMP4,                         &
     &LBD, UBD
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      INFO = 0
!
      IF( ORGATI ) THEN
        LBD = D(2)
        UBD = D(3)
      ELSE
        LBD = D(1)
        UBD = D(2)
      END IF
      IF( FINIT .LT. ZERO )THEN
        LBD = ZERO
      ELSE
        UBD = ZERO
      END IF
!
      NITER = 1
      TAU = ZERO
      IF( KNITER.EQ.2 ) THEN
        IF( ORGATI ) THEN
          TEMP = ( D( 3 )-D( 2 ) ) / TWO
          C = RHO + Z( 1 ) / ( ( D( 1 )-D( 2 ) )-TEMP )
          A = C*( D( 2 )+D( 3 ) ) + Z( 2 ) + Z( 3 )
          B = C*D( 2 )*D( 3 ) + Z( 2 )*D( 3 ) + Z( 3 )*D( 2 )
        ELSE
          TEMP = ( D( 1 )-D( 2 ) ) / TWO
          C = RHO + Z( 3 ) / ( ( D( 3 )-D( 2 ) )-TEMP )
          A = C*( D( 1 )+D( 2 ) ) + Z( 1 ) + Z( 2 )
          B = C*D( 1 )*D( 2 ) + Z( 1 )*D( 2 ) + Z( 2 )*D( 1 )
        END IF
        TEMP = MAX( ABS( A ), ABS( B ), ABS( C ) )
        A = A / TEMP
        B = B / TEMP
        C = C / TEMP
        IF( C.EQ.ZERO ) THEN
          TAU = B / A
        ELSE IF( A.LE.ZERO ) THEN
          TAU = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
        ELSE
          TAU = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
        END IF
        IF( TAU .LT. LBD .OR. TAU .GT. UBD )                            &
     &TAU = ( LBD+UBD )/TWO
        TEMP = FINIT + TAU*Z(1)/( D(1)*( D( 1 )-TAU ) ) +               &
     &TAU*Z(2)/( D(2)*( D( 2 )-TAU ) ) +                                &
     &TAU*Z(3)/( D(3)*( D( 3 )-TAU ) )
        IF( TEMP .LE. ZERO )THEN
          LBD = TAU
        ELSE
          UBD = TAU
        END IF
        IF( ABS( FINIT ).LE.ABS( TEMP ) )                               &
     &TAU = ZERO
      END IF
!
!     get machine parameters for possible scaling to avoid overflow
!
!     modified by Sven: parameters SMALL1, SMINV1, SMALL2,
!     SMINV2, EPS are not SAVEd anymore between one call to the
!     others but recomputed at each call
!
      EPS = DLAMCH( 'Epsilon' )
      BASE = DLAMCH( 'Base' )
      SMALL1 = BASE**( INT( LOG( DLAMCH( 'SafMin' ) ) / LOG( BASE ) /   &
     &THREE ) )
      SMINV1 = ONE / SMALL1
      SMALL2 = SMALL1*SMALL1
      SMINV2 = SMINV1*SMINV1
!
!     Determine if scaling of inputs necessary to avoid overflow
!     when computing 1/TEMP**3
!
      IF( ORGATI ) THEN
        TEMP = MIN( ABS( D( 2 )-TAU ), ABS( D( 3 )-TAU ) )
      ELSE
        TEMP = MIN( ABS( D( 1 )-TAU ), ABS( D( 2 )-TAU ) )
      END IF
      SCALE = .FALSE.
      IF( TEMP.LE.SMALL1 ) THEN
        SCALE = .TRUE.
        IF( TEMP.LE.SMALL2 ) THEN
!
!        Scale up by power of radix nearest 1/SAFMIN**(2/3)
!
          SCLFAC = SMINV2
          SCLINV = SMALL2
        ELSE
!
!        Scale up by power of radix nearest 1/SAFMIN**(1/3)
!
          SCLFAC = SMINV1
          SCLINV = SMALL1
        END IF
!
!        Scaling up safe because D, Z, TAU scaled elsewhere to be O(1)
!
        DO I = 1, 3
          DSCALE( I ) = D( I )*SCLFAC
          ZSCALE( I ) = Z( I )*SCLFAC
        enddo
        TAU = TAU*SCLFAC
        LBD = LBD*SCLFAC
        UBD = UBD*SCLFAC
      ELSE
!
!        Copy D and Z to DSCALE and ZSCALE
!
        DO I = 1, 3
          DSCALE( I ) = D( I )
          ZSCALE( I ) = Z( I )
        enddo
      END IF
!
      FC = ZERO
      DF = ZERO
      DDF = ZERO
      DO I = 1, 3
        TEMP = ONE / ( DSCALE( I )-TAU )
        TEMP1 = ZSCALE( I )*TEMP
        TEMP2 = TEMP1*TEMP
        TEMP3 = TEMP2*TEMP
        FC = FC + TEMP1 / DSCALE( I )
        DF = DF + TEMP2
        DDF = DDF + TEMP3
      enddo
      F = FINIT + TAU*FC
!
      IF( ABS( F ).LE.ZERO )                                            &
     &GO TO 60
      IF( F .LE. ZERO )THEN
        LBD = TAU
      ELSE
        UBD = TAU
      END IF
!
!        Iteration begins -- Use Gragg-Thornton-Warner cubic convergent
!                            scheme
!
!     It is not hard to see that
!
!           1) Iterations will go up monotonically
!              if FINIT < 0;
!
!           2) Iterations will go down monotonically
!              if FINIT > 0.
!
      ITER = NITER + 1
!
      DO NITER = ITER, MAXIT
!
        IF( ORGATI ) THEN
          TEMP1 = DSCALE( 2 ) - TAU
          TEMP2 = DSCALE( 3 ) - TAU
        ELSE
          TEMP1 = DSCALE( 1 ) - TAU
          TEMP2 = DSCALE( 2 ) - TAU
        END IF
        A = ( TEMP1+TEMP2 )*F - TEMP1*TEMP2*DF
        B = TEMP1*TEMP2*F
        C = F - ( TEMP1+TEMP2 )*DF + TEMP1*TEMP2*DDF
        TEMP = MAX( ABS( A ), ABS( B ), ABS( C ) )
        A = A / TEMP
        B = B / TEMP
        C = C / TEMP
        IF( C.EQ.ZERO ) THEN
          ETA = B / A
        ELSE IF( A.LE.ZERO ) THEN
          ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
        ELSE
          ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
        END IF
        IF( F*ETA.GE.ZERO ) THEN
          ETA = -F / DF
        END IF
!
        TAU = TAU + ETA
        IF( TAU .LT. LBD .OR. TAU .GT. UBD )                            &
     &TAU = ( LBD + UBD )/TWO
!
        FC = ZERO
        ERRETM = ZERO
        DF = ZERO
        DDF = ZERO
        DO I = 1, 3
          TEMP = ONE / ( DSCALE( I )-TAU )
          TEMP1 = ZSCALE( I )*TEMP
          TEMP2 = TEMP1*TEMP
          TEMP3 = TEMP2*TEMP
          TEMP4 = TEMP1 / DSCALE( I )
          FC = FC + TEMP4
          ERRETM = ERRETM + ABS( TEMP4 )
          DF = DF + TEMP2
          DDF = DDF + TEMP3
        enddo
        F = FINIT + TAU*FC
        ERRETM = EIGHT*( ABS( FINIT )+ABS( TAU )*ERRETM ) +             &
     &ABS( TAU )*DF
        IF( ABS( F ).LE.EPS*ERRETM )                                    &
     &GO TO 60
        IF( F .LE. ZERO )THEN
          LBD = TAU
        ELSE
          UBD = TAU
        END IF
      enddo
      INFO = 1
   60 CONTINUE
!
!     Undo scaling
!
      IF( SCALE )                                                       &
     &TAU = TAU*SCLINV
      RETURN
!
!     End of DLAED6
!
      END
      SUBROUTINE DLAMRG( N1, N2, A, DTRD1, DTRD2, INDEX )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            DTRD1, DTRD2, N1, N2
!     ..
!     .. Array Arguments ..
      INTEGER            INDEX( * )
      DOUBLE PRECISION   A( * )
!     ..
!
!  Purpose
!  =======
!
!  DLAMRG will create a permutation list which will merge the elements
!  of A (which is composed of two independently sorted sets) into a
!  single set which is sorted in ascending order.
!
!  Arguments
!  =========
!
!  N1     (input) INTEGER
!  N2     (input) INTEGER
!         These arguements contain the respective lengths of the two
!         sorted lists to be merged.
!
!  A      (input) DOUBLE PRECISION array, dimension (N1+N2)
!         The first N1 elements of A contain a list of numbers which
!         are sorted in either ascending or descending order.  Likewise
!         for the final N2 elements.
!
!  DTRD1  (input) INTEGER
!  DTRD2  (input) INTEGER
!         These are the strides to be taken through the array A.
!         Allowable strides are 1 and -1.  They indicate whether a
!         subset of A is sorted in ascending (DTRDx = 1) or descending
!         (DTRDx = -1) order.
!
!  INDEX  (output) INTEGER array, dimension (N1+N2)
!         On exit this array will contain a permutation such that
!         if B( I ) = A( INDEX( I ) ) for I=1,N1+N2, then B will be
!         sorted in ascending order.
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IND1, IND2, N1SV, N2SV
!     ..
!     .. Executable Statements ..
!
      N1SV = N1
      N2SV = N2
      IF( DTRD1.GT.0 ) THEN
        IND1 = 1
      ELSE
        IND1 = N1
      END IF
      IF( DTRD2.GT.0 ) THEN
        IND2 = 1 + N1
      ELSE
        IND2 = N1 + N2
      END IF
      I = 1
!     while ( (N1SV > 0) & (N2SV > 0) )
   10 CONTINUE
      IF( N1SV.GT.0 .AND. N2SV.GT.0 ) THEN
        IF( A( IND1 ).LE.A( IND2 ) ) THEN
          INDEX( I ) = IND1
          I = I + 1
          IND1 = IND1 + DTRD1
          N1SV = N1SV - 1
        ELSE
          INDEX( I ) = IND2
          I = I + 1
          IND2 = IND2 + DTRD2
          N2SV = N2SV - 1
        END IF
        GO TO 10
      END IF
!     end while
      IF( N1SV.EQ.0 ) THEN
        DO N1SV = 1, N2SV
          INDEX( I ) = IND2
          I = I + 1
          IND2 = IND2 + DTRD2
        enddo
      ELSE
!     N2SV .EQ. 0
        DO N2SV = 1, N1SV
          INDEX( I ) = IND1
          I = I + 1
          IND1 = IND1 + DTRD1
        enddo
      END IF
!
      RETURN
!
!     End of DLAMRG
!
      END
      SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      DOUBLE PRECISION   TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DLARF applies a real elementary reflector H to a real m by n matrix
!  C, from either the left or the right. H is represented in the form
!
!        H = I - tau * v * v'
!
!  where tau is a real scalar and v is a real vector.
!
!  If tau = 0, then H is taken to be the unit matrix.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': form  H * C
!          = 'R': form  C * H
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  V       (input) DOUBLE PRECISION array, dimension
!                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!          The vector v in the representation of H. V is not used if
!          TAU = 0.
!
!  INCV    (input) INTEGER
!          The increment between elements of v. INCV <> 0.
!
!  TAU     (input) DOUBLE PRECISION
!          The value tau in the representation of H.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!          or C * H if SIDE = 'R'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension
!                         (N) if SIDE = 'L'
!                      or (M) if SIDE = 'R'
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( SIDE, 'L' ) ) THEN
!
!        Form  H * C
!
        IF( TAU.NE.ZERO ) THEN
!
!           w := C' * v
!
          CALL DGEMV( 'Transpose', M, N, ONE, C, LDC, V, INCV, ZERO,    &
     &WORK, 1 )
!
!           C := C - v * w'
!
          CALL DGER( M, N, -TAU, V, INCV, WORK, 1, C, LDC )
        END IF
      ELSE
!
!        Form  C * H
!
        IF( TAU.NE.ZERO ) THEN
!
!           w := C * v
!
          CALL DGEMV( 'No transpose', M, N, ONE, C, LDC, V, INCV,       &
     &ZERO, WORK, 1 )
!
!           C := C - w * v'
!
          CALL DGER( M, N, -TAU, WORK, 1, V, INCV, C, LDC )
        END IF
      END IF
      RETURN
!
!     End of DLARF
!
      END
      SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,  &
     &T, LDT, C, LDC, WORK, LDWORK )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ),         &
     &WORK( LDWORK, * )
!     ..
!
!  Purpose
!  =======
!
!  DLARFB applies a real block reflector H or its transpose H' to a
!  real m by n matrix C, from either the left or the right.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply H or H' from the Left
!          = 'R': apply H or H' from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply H (No transpose)
!          = 'T': apply H' (Transpose)
!
!  DIRECT  (input) CHARACTER*1
!          Indicates how H is formed from a product of elementary
!          reflectors
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!
!  STOREV  (input) CHARACTER*1
!          Indicates how the vectors which define the elementary
!          reflectors are stored:
!          = 'C': Columnwise
!          = 'R': Rowwise
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  K       (input) INTEGER
!          The order of the matrix T (= the number of elementary
!          reflectors whose product defines the block reflector).
!
!  V       (input) DOUBLE PRECISION array, dimension
!                                (LDV,K) if STOREV = 'C'
!                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
!                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
!          The matrix V. See further details.
!
!  LDV     (input) INTEGER
!          The leading dimension of the array V.
!          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
!          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
!          if STOREV = 'R', LDV >= K.
!
!  T       (input) DOUBLE PRECISION array, dimension (LDT,K)
!          The triangular k by k matrix T in the representation of the
!          block reflector.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= K.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDA >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,K)
!
!  LDWORK  (input) INTEGER
!          The leading dimension of the array WORK.
!          If SIDE = 'L', LDWORK >= max(1,N);
!          if SIDE = 'R', LDWORK >= max(1,M).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPL, DGEMM, DTRMM
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( M.LE.0 .OR. N.LE.0 )                                          &
     &RETURN
!
      IF( LSAME( TRANS, 'N' ) ) THEN
        TRANST = 'T'
      ELSE
        TRANST = 'N'
      END IF
!
      IF( LSAME( STOREV, 'C' ) ) THEN
!
        IF( LSAME( DIRECT, 'F' ) ) THEN
!
!           Let  V =  ( V1 )    (first K rows)
!                     ( V2 )
!           where  V1  is unit lower triangular.
!
          IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
!
!              W := C1'
!
            DO J = 1, K
              CALL DCOPL( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
            enddo
!
!              W := W * V1
!
            CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,    &
     &K, ONE, V, LDV, WORK, LDWORK )
            IF( M.GT.K ) THEN
!
!                 W := W + C2'*V2
!
              CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K,       &
     &ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV,                          &
     &ONE, WORK, LDWORK )
            END IF
!
!              W := W * T'  or  W * T
!
            CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,     &
     &ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W'
!
            IF( M.GT.K ) THEN
!
!                 C2 := C2 - V2 * W'
!
              CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K,       &
     &-ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE,                        &
     &C( K+1, 1 ), LDC )
            END IF
!
!              W := W * V1'
!
            CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K,    &
     &ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W'
!
            DO J = 1, K
              DO I = 1, N
                C( J, I ) = C( J, I ) - WORK( I, J )
              enddo
            enddo
!
          ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C1
!
            DO J = 1, K
              CALL DCOPL( M, C( 1, J ), 1, WORK( 1, J ), 1 )
            enddo
!
!              W := W * V1
!
            CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,    &
     &K, ONE, V, LDV, WORK, LDWORK )
            IF( N.GT.K ) THEN
!
!                 W := W + C2 * V2
!
              CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K,    &
     &ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV,                          &
     &ONE, WORK, LDWORK )
            END IF
!
!              W := W * T  or  W * T'
!
            CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,      &
     &ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V'
!
            IF( N.GT.K ) THEN
!
!                 C2 := C2 - W * V2'
!
              CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K,       &
     &-ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE,                        &
     &C( 1, K+1 ), LDC )
            END IF
!
!              W := W * V1'
!
            CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K,    &
     &ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
            DO J = 1, K
              DO I = 1, M
                C( I, J ) = C( I, J ) - WORK( I, J )
              enddo
            enddo
          END IF
!
        ELSE
!
!           Let  V =  ( V1 )
!                     ( V2 )    (last K rows)
!           where  V2  is unit upper triangular.
!
          IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
!
!              W := C2'
!
            DO J = 1, K
              CALL DCOPL( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
            enddo
!
!              W := W * V2
!
            CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,    &
     &K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
            IF( M.GT.K ) THEN
!
!                 W := W + C1'*V1
!
              CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K,       &
     &ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
            END IF
!
!              W := W * T'  or  W * T
!
            CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,     &
     &ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W'
!
            IF( M.GT.K ) THEN
!
!                 C1 := C1 - V1 * W'
!
              CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K,       &
     &-ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
            END IF
!
!              W := W * V2'
!
            CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K,    &
     &ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W'
!
            DO J = 1, K
              DO I = 1, N
                C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
              enddo
            enddo
!
          ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C2
!
            DO J = 1, K
              CALL DCOPL( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
            enddo
!
!              W := W * V2
!
            CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,    &
     &K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
            IF( N.GT.K ) THEN
!
!                 W := W + C1 * V1
!
              CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K,    &
     &ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
            END IF
!
!              W := W * T  or  W * T'
!
            CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,      &
     &ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V'
!
            IF( N.GT.K ) THEN
!
!                 C1 := C1 - W * V1'
!
              CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K,       &
     &-ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
            END IF
!
!              W := W * V2'
!
            CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K,    &
     &ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W
!
            DO J = 1, K
              DO I = 1, M
                C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
              enddo
            enddo
          END IF
        END IF
!
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
!
        IF( LSAME( DIRECT, 'F' ) ) THEN
!
!           Let  V =  ( V1  V2 )    (V1: first K columns)
!           where  V1  is unit upper triangular.
!
          IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
!
!              W := C1'
!
            DO J = 1, K
              CALL DCOPL( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
            enddo
!
!              W := W * V1'
!
            CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K,    &
     &ONE, V, LDV, WORK, LDWORK )
            IF( M.GT.K ) THEN
!
!                 W := W + C2'*V2'
!
              CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,     &
     &C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE,                          &
     &WORK, LDWORK )
            END IF
!
!              W := W * T'  or  W * T
!
            CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,     &
     &ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V' * W'
!
            IF( M.GT.K ) THEN
!
!                 C2 := C2 - V2' * W'
!
              CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE,    &
     &V( 1, K+1 ), LDV, WORK, LDWORK, ONE,                              &
     &C( K+1, 1 ), LDC )
            END IF
!
!              W := W * V1
!
            CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,    &
     &K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W'
!
            DO J = 1, K
              DO I = 1, N
                C( J, I ) = C( J, I ) - WORK( I, J )
              enddo
            enddo
!
          ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
!
!              W := C1
!
            DO J = 1, K
              CALL DCOPL( M, C( 1, J ), 1, WORK( 1, J ), 1 )
            enddo
!
!              W := W * V1'
!
            CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K,    &
     &ONE, V, LDV, WORK, LDWORK )
            IF( N.GT.K ) THEN
!
!                 W := W + C2 * V2'
!
              CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K,       &
     &ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV,                          &
     &ONE, WORK, LDWORK )
            END IF
!
!              W := W * T  or  W * T'
!
            CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,      &
     &ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
            IF( N.GT.K ) THEN
!
!                 C2 := C2 - W * V2
!
              CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K,    &
     &-ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE,                        &
     &C( 1, K+1 ), LDC )
            END IF
!
!              W := W * V1
!
            CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,    &
     &K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
            DO J = 1, K
              DO I = 1, M
                C( I, J ) = C( I, J ) - WORK( I, J )
              enddo
            enddo
!
          END IF
!
        ELSE
!
!           Let  V =  ( V1  V2 )    (V2: last K columns)
!           where  V2  is unit lower triangular.
!
          IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
!
!              W := C2'
!
            DO J = 1, K
              CALL DCOPL( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
            enddo
!
!              W := W * V2'
!
            CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K,    &
     &ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
            IF( M.GT.K ) THEN
!
!                 W := W + C1'*V1'
!
              CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,     &
     &C, LDC, V, LDV, ONE, WORK, LDWORK )
            END IF
!
!              W := W * T'  or  W * T
!
            CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,     &
     &ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V' * W'
!
            IF( M.GT.K ) THEN
!
!                 C1 := C1 - V1' * W'
!
              CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE,    &
     &V, LDV, WORK, LDWORK, ONE, C, LDC )
            END IF
!
!              W := W * V2
!
            CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,    &
     &K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W'
!
            DO J = 1, K
              DO I = 1, N
                C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
              enddo
            enddo
!
          ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
!
!              W := C2
!
            DO J = 1, K
              CALL DCOPL( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
            enddo
!
!              W := W * V2'
!
            CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K,    &
     &ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
            IF( N.GT.K ) THEN
!
!                 W := W + C1 * V1'
!
              CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K,       &
     &ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
            END IF
!
!              W := W * T  or  W * T'
!
            CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,      &
     &ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
            IF( N.GT.K ) THEN
!
!                 C1 := C1 - W * V1
!
              CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K,    &
     &-ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
            END IF
!
!              W := W * V2
!
            CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,    &
     &K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
            DO J = 1, K
              DO I = 1, M
                C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
              enddo
            enddo
!
          END IF
!
        END IF
      END IF
!
      RETURN
!
!     End of DLARFB
!
      END
      SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   ALPHA, TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
!     ..
!
!  Purpose
!  =======
!
!  DLARFG generates a real elementary reflector H of order n, such
!  that
!
!        H * ( alpha ) = ( beta ),   H' * H = I.
!            (   x   )   (   0  )
!
!  where alpha and beta are scalars, and x is an (n-1)-element real
!  vector. H is represented in the form
!
!        H = I - tau * ( 1 ) * ( 1 v' ) ,
!                      ( v )
!
!  where tau is a real scalar and v is a real (n-1)-element
!  vector.
!
!  If the elements of x are all zero, then tau = 0 and H is taken to be
!  the unit matrix.
!
!  Otherwise  1 <= tau <= 2.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the elementary reflector.
!
!  ALPHA   (input/output) DOUBLE PRECISION
!          On entry, the value alpha.
!          On exit, it is overwritten with the value beta.
!
!  X       (input/output) DOUBLE PRECISION array, dimension
!                         (1+(N-2)*abs(INCX))
!          On entry, the vector x.
!          On exit, it is overwritten with the vector v.
!
!  INCX    (input) INTEGER
!          The increment between elements of X. INCX > 0.
!
!  TAU     (output) DOUBLE PRECISION
!          The value tau.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   BETA, RSAFMN, SAFMIN, XNORM
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2, DNRM2
      EXTERNAL           DLAMCH, DLAPY2, DNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
!     ..
!     .. External Subroutines ..
      EXTERNAL           DSCAL
!     ..
!     .. Executable Statements ..
!
      IF( N.LE.1 ) THEN
        TAU = ZERO
        RETURN
      END IF
!
      XNORM = DNRM2( N-1, X, INCX )
!
      IF( XNORM.EQ.ZERO ) THEN
!
!        H  =  I
!
        TAU = ZERO
      ELSE
!
!        general case
!
        BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
        SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
        IF( ABS( BETA ).LT.SAFMIN ) THEN
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
          RSAFMN = ONE / SAFMIN
          KNT = 0
   10     CONTINUE
          KNT = KNT + 1
          CALL DSCAL( N-1, RSAFMN, X, INCX )
          BETA = BETA*RSAFMN
          ALPHA = ALPHA*RSAFMN
          IF( ABS( BETA ).LT.SAFMIN )                                   &
     &GO TO 10
!
!           New BETA is at most 1, at least SAFMIN
!
          XNORM = DNRM2( N-1, X, INCX )
          BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
          TAU = ( BETA-ALPHA ) / BETA
          CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
!
!           If ALPHA is subnormal, it may lose relative accuracy
!
          ALPHA = BETA
          DO J = 1, KNT
            ALPHA = ALPHA*SAFMIN
          enddo
        ELSE
          TAU = ( BETA-ALPHA ) / BETA
          CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
          ALPHA = BETA
        END IF
      END IF
!
      RETURN
!
!     End of DLARFG
!
      END
      SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   T( LDT, * ), TAU( * ), V( LDV, * )
!     ..
!
!  Purpose
!  =======
!
!  DLARFT forms the triangular factor T of a real block reflector H
!  of order n, which is defined as a product of k elementary reflectors.
!
!  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
!
!  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
!
!  If STOREV = 'C', the vector which defines the elementary reflector
!  H(i) is stored in the i-th column of the array V, and
!
!     H  =  I - V * T * V'
!
!  If STOREV = 'R', the vector which defines the elementary reflector
!  H(i) is stored in the i-th row of the array V, and
!
!     H  =  I - V' * T * V
!
!  Arguments
!  =========
!
!  DIRECT  (input) CHARACTER*1
!          Specifies the order in which the elementary reflectors are
!          multiplied to form the block reflector:
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!
!  STOREV  (input) CHARACTER*1
!          Specifies how the vectors which define the elementary
!          reflectors are stored (see also Further Details):
!          = 'C': columnwise
!          = 'R': rowwise
!
!  N       (input) INTEGER
!          The order of the block reflector H. N >= 0.
!
!  K       (input) INTEGER
!          The order of the triangular factor T (= the number of
!          elementary reflectors). K >= 1.
!
!  V       (input/output) DOUBLE PRECISION array, dimension
!                               (LDV,K) if STOREV = 'C'
!                               (LDV,N) if STOREV = 'R'
!          The matrix V. See further details.
!
!  LDV     (input) INTEGER
!          The leading dimension of the array V.
!          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i).
!
!  T       (output) DOUBLE PRECISION array, dimension (LDT,K)
!          The k by k triangular factor T of the block reflector.
!          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
!          lower triangular. The rest of the array is not used.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= K.
!
!  Further Details
!  ===============
!
!  The shape of the matrix V and the storage of the vectors which define
!  the H(i) is best illustrated by the following example with n = 5 and
!  k = 3. The elements equal to 1 are not stored; the corresponding
!  array elements are modified but restored on exit. The rest of the
!  array is not used.
!
!  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
!
!               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!                   ( v1  1    )                     (     1 v2 v2 v2 )
!                   ( v1 v2  1 )                     (        1 v3 v3 )
!                   ( v1 v2 v3 )
!                   ( v1 v2 v3 )
!
!  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
!
!               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!                   (     1 v3 )
!                   (        1 )
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   VII
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMV, DTRMV
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.EQ.0 )                                                      &
     &RETURN
!
      IF( LSAME( DIRECT, 'F' ) ) THEN
        DO I = 1, K
          IF( TAU( I ).EQ.ZERO ) THEN
!
!              H(i)  =  I
!
            DO J = 1, I
              T( J, I ) = ZERO
            enddo
          ELSE
!
!              general case
!
            VII = V( I, I )
            V( I, I ) = ONE
            IF( LSAME( STOREV, 'C' ) ) THEN
!
!                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i)
!
              CALL DGEMV( 'Transpose', N-I+1, I-1, -TAU( I ),           &
     &V( I, 1 ), LDV, V( I, I ), 1, ZERO,                               &
     &T( 1, I ), 1 )
            ELSE
!
!                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) * V(i,i:n)'
!
              CALL DGEMV( 'No transpose', I-1, N-I+1, -TAU( I ),        &
     &V( 1, I ), LDV, V( I, I ), LDV, ZERO,                             &
     &T( 1, I ), 1 )
            END IF
            V( I, I ) = VII
!
!              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
!
            CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T,    &
     &LDT, T( 1, I ), 1 )
            T( I, I ) = TAU( I )
          END IF
        enddo
      ELSE
        DO I = K, 1, -1
          IF( TAU( I ).EQ.ZERO ) THEN
!
!              H(i)  =  I
!
            DO J = I, K
              T( J, I ) = ZERO
            enddo
          ELSE
!
!              general case
!
            IF( I.LT.K ) THEN
              IF( LSAME( STOREV, 'C' ) ) THEN
                VII = V( N-K+I, I )
                V( N-K+I, I ) = ONE
!
!                    T(i+1:k,i) :=
!                            - tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i)
!
                CALL DGEMV( 'Transpose', N-K+I, K-I, -TAU( I ),         &
     &V( 1, I+1 ), LDV, V( 1, I ), 1, ZERO,                             &
     &T( I+1, I ), 1 )
                V( N-K+I, I ) = VII
              ELSE
                VII = V( I, N-K+I )
                V( I, N-K+I ) = ONE
!
!                    T(i+1:k,i) :=
!                            - tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)'
!
                CALL DGEMV( 'No transpose', K-I, N-K+I, -TAU( I ),      &
     &V( I+1, 1 ), LDV, V( I, 1 ), LDV, ZERO,                           &
     &T( I+1, I ), 1 )
                V( I, N-K+I ) = VII
              END IF
!
!                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
!
              CALL DTRMV( 'Lower', 'No transpose', 'Non-unit', K-I,     &
     &T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
            END IF
            T( I, I ) = TAU( I )
          END IF
        enddo
      END IF
      RETURN
!
!     End of DLARFT
!
      END
      SUBROUTINE DLARTG( F, G, CS, SN, R )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   CS, F, G, R, SN
!     ..
!
!  Purpose
!  =======
!
!  DLARTG generate a plane rotation so that
!
!     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
!     [ -SN  CS  ]     [ G ]     [ 0 ]
!
!  This is a slower, more accurate version of the BLAS1 routine DROTG,
!  with the following other differences:
!     F and G are unchanged on return.
!     If G=0, then CS=1 and SN=0.
!     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
!        floating point operations (saves work in DBDSQR when
!        there are zeros on the diagonal).
!
!  If F exceeds G in magnitude, CS will be positive.
!
!  Arguments
!  =========
!
!  F       (input) DOUBLE PRECISION
!          The first component of vector to be rotated.
!
!  G       (input) DOUBLE PRECISION
!          The second component of vector to be rotated.
!
!  CS      (output) DOUBLE PRECISION
!          The cosine of the rotation.
!
!  SN      (output) DOUBLE PRECISION
!          The sine of the rotation.
!
!  R       (output) DOUBLE PRECISION
!          The nonzero component of the rotated vector.
!
!  This version has a few statements commented out for thread safety
!  (machine parameters are computed on each entry). 10 feb 03, SJH.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
!     ..
!     .. Local Scalars ..
!     LOGICAL            FIRST
      INTEGER            COUNT, I
      DOUBLE PRECISION   EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, SQRT
!     ..
!     .. Save statement ..
!     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
!     ..
!     .. Data statements ..
!     DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
!     IF( FIRST ) THEN
      SAFMIN = DLAMCH( 'S' )
      EPS = DLAMCH( 'E' )
      SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) /                &
     &LOG( DLAMCH( 'B' ) ) / TWO )
      SAFMX2 = ONE / SAFMN2
!        FIRST = .FALSE.
!     END IF
      IF( G.EQ.ZERO ) THEN
        CS = ONE
        SN = ZERO
        R = F
      ELSE IF( F.EQ.ZERO ) THEN
        CS = ZERO
        SN = ONE
        R = G
      ELSE
        F1 = F
        G1 = G
        SCALE = MAX( ABS( F1 ), ABS( G1 ) )
        IF( SCALE.GE.SAFMX2 ) THEN
          COUNT = 0
   10     CONTINUE
          COUNT = COUNT + 1
          F1 = F1*SAFMN2
          G1 = G1*SAFMN2
          SCALE = MAX( ABS( F1 ), ABS( G1 ) )
          IF( SCALE.GE.SAFMX2 )                                         &
     &GO TO 10
          R = SQRT( F1**2+G1**2 )
          CS = F1 / R
          SN = G1 / R
          DO I = 1, COUNT
            R = R*SAFMX2
          enddo
        ELSE IF( SCALE.LE.SAFMN2 ) THEN
          COUNT = 0
   30     CONTINUE
          COUNT = COUNT + 1
          F1 = F1*SAFMX2
          G1 = G1*SAFMX2
          SCALE = MAX( ABS( F1 ), ABS( G1 ) )
          IF( SCALE.LE.SAFMN2 )                                         &
     &GO TO 30
          R = SQRT( F1**2+G1**2 )
          CS = F1 / R
          SN = G1 / R
          DO I = 1, COUNT
            R = R*SAFMN2
          enddo
        ELSE
          R = SQRT( F1**2+G1**2 )
          CS = F1 / R
          SN = G1 / R
        END IF
        IF( ABS( F ).GT.ABS( G ) .AND. CS.LT.ZERO ) THEN
          CS = -CS
          SN = -SN
          R = -R
        END IF
      END IF
      RETURN
!
!     End of DLARTG
!
      END
      SUBROUTINE DLAS2( F, G, H, SSMIN, SSMAX )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   F, G, H, SSMAX, SSMIN
!     ..
!
!  Purpose
!  =======
!
!  DLAS2  computes the singular values of the 2-by-2 matrix
!     [  F   G  ]
!     [  0   H  ].
!  On return, SSMIN is the smaller singular value and SSMAX is the
!  larger singular value.
!
!  Arguments
!  =========
!
!  F       (input) DOUBLE PRECISION
!          The (1,1) element of the 2-by-2 matrix.
!
!  G       (input) DOUBLE PRECISION
!          The (1,2) element of the 2-by-2 matrix.
!
!  H       (input) DOUBLE PRECISION
!          The (2,2) element of the 2-by-2 matrix.
!
!  SSMIN   (output) DOUBLE PRECISION
!          The smaller singular value.
!
!  SSMAX   (output) DOUBLE PRECISION
!          The larger singular value.
!
!  Further Details
!  ===============
!
!  Barring over/underflow, all output quantities are correct to within
!  a few units in the last place (ulps), even in the absence of a guard
!  digit in addition/subtraction.
!
!  In IEEE arithmetic, the code works correctly if one matrix element is
!  infinite.
!
!  Overflow will not occur unless the largest singular value itself
!  overflows, or is within a few ulps of overflow. (On machines with
!  partial overflow, like the Cray, overflow may occur if the largest
!  singular value is within a factor of 2 of overflow.)
!
!  Underflow is harmless if underflow is gradual. Otherwise, results
!  may correspond to a matrix modified by perturbations of size near
!  the underflow threshold.
!
!  ====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   AS, AT, AU, C, FA, FHMN, FHMX, GA, HA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      FA = ABS( F )
      GA = ABS( G )
      HA = ABS( H )
      FHMN = MIN( FA, HA )
      FHMX = MAX( FA, HA )
      IF( FHMN.EQ.ZERO ) THEN
        SSMIN = ZERO
        IF( FHMX.EQ.ZERO ) THEN
          SSMAX = GA
        ELSE
          SSMAX = MAX( FHMX, GA )*SQRT( ONE+                            &
     &( MIN( FHMX, GA ) / MAX( FHMX, GA ) )**2 )
        END IF
      ELSE
        IF( GA.LT.FHMX ) THEN
          AS = ONE + FHMN / FHMX
          AT = ( FHMX-FHMN ) / FHMX
          AU = ( GA / FHMX )**2
          C = TWO / ( SQRT( AS*AS+AU )+SQRT( AT*AT+AU ) )
          SSMIN = FHMN*C
          SSMAX = FHMX / C
        ELSE
          AU = FHMX / GA
          IF( AU.EQ.ZERO ) THEN
!
!              Avoid possible harmful underflow if exponent range
!              asymmetric (true SSMIN may not underflow even if
!              AU underflows)
!
            SSMIN = ( FHMN*FHMX ) / GA
            SSMAX = GA
          ELSE
            AS = ONE + FHMN / FHMX
            AT = ( FHMX-FHMN ) / FHMX
            C = ONE / ( SQRT( ONE+( AS*AU )**2 )+                       &
     &SQRT( ONE+( AT*AU )**2 ) )
            SSMIN = ( FHMN*C )*AU
            SSMIN = SSMIN + SSMIN
            SSMAX = GA / ( C+C )
          END IF
        END IF
      END IF
      RETURN
!
!     End of DLAS2
!
      END
      SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      DOUBLE PRECISION   CFROM, CTO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DLASCL multiplies the M by N real matrix A by the real scalar
!  CTO/CFROM.  This is done without over/underflow as long as the final
!  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
!  A may be full, upper triangular, lower triangular, upper Hessenberg,
!  or banded.
!
!  Arguments
!  =========
!
!  TYPE    (input) CHARACTER*1
!          TYPE indices the storage type of the input matrix.
!          = 'G':  A is a full matrix.
!          = 'L':  A is a lower triangular matrix.
!          = 'U':  A is an upper triangular matrix.
!          = 'H':  A is an upper Hessenberg matrix.
!          = 'B':  A is a symmetric band matrix with lower bandwidth KL
!                  and upper bandwidth KU and with the only the lower
!                  half stored.
!          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
!                  and upper bandwidth KU and with the only the upper
!                  half stored.
!          = 'Z':  A is a band matrix with lower bandwidth KL and upper
!                  bandwidth KU.
!
!  KL      (input) INTEGER
!          The lower bandwidth of A.  Referenced only if TYPE = 'B',
!          'Q' or 'Z'.
!
!  KU      (input) INTEGER
!          The upper bandwidth of A.  Referenced only if TYPE = 'B',
!          'Q' or 'Z'.
!
!  CFROM   (input) DOUBLE PRECISION
!  CTO     (input) DOUBLE PRECISION
!          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
!          without over/underflow if the final result CTO*A(I,J)/CFROM
!          can be represented without over/underflow.  CFROM must be
!          nonzero.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
!          storage type.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  INFO    (output) INTEGER
!          0  - successful exit
!          <0 - if INFO = -i, the i-th argument had an illegal value.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
!
      IF( LSAME( TYPE, 'G' ) ) THEN
        ITYPE = 0
      ELSE IF( LSAME( TYPE, 'L' ) ) THEN
        ITYPE = 1
      ELSE IF( LSAME( TYPE, 'U' ) ) THEN
        ITYPE = 2
      ELSE IF( LSAME( TYPE, 'H' ) ) THEN
        ITYPE = 3
      ELSE IF( LSAME( TYPE, 'B' ) ) THEN
        ITYPE = 4
      ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
        ITYPE = 5
      ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
        ITYPE = 6
      ELSE
        ITYPE = -1
      END IF
!
      IF( ITYPE.EQ.-1 ) THEN
        INFO = -1
      ELSE IF( CFROM.EQ.ZERO ) THEN
        INFO = -4
      ELSE IF( M.LT.0 ) THEN
        INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR.             &
     &( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
        INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
        INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
        IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
          INFO = -2
        ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR.                  &
     &( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) )               &
     &THEN
          INFO = -3
        ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR.                  &
     &( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR.                             &
     &( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
          INFO = -9
        END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DLASCL', -INFO )
        RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. M.EQ.0 )                                          &
     &RETURN
!
!     Get machine parameters
!
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
!
      CFROMC = CFROM
      CTOC = CTO
!
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      CTO1 = CTOC / BIGNUM
      IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
        MUL = SMLNUM
        DONE = .FALSE.
        CFROMC = CFROM1
      ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
        MUL = BIGNUM
        DONE = .FALSE.
        CTOC = CTO1
      ELSE
        MUL = CTOC / CFROMC
        DONE = .TRUE.
      END IF
!
      IF( ITYPE.EQ.0 ) THEN
!
!        Full matrix
!
        DO J = 1, N
          DO I = 1, M
            A( I, J ) = A( I, J )*MUL
          enddo
        enddo
!
      ELSE IF( ITYPE.EQ.1 ) THEN
!
!        Lower triangular matrix
!
        DO J = 1, N
          DO I = J, M
            A( I, J ) = A( I, J )*MUL
          enddo
        enddo
!
      ELSE IF( ITYPE.EQ.2 ) THEN
!
!        Upper triangular matrix
!
        DO J = 1, N
          DO I = 1, MIN( J, M )
            A( I, J ) = A( I, J )*MUL
          enddo
        enddo
!
      ELSE IF( ITYPE.EQ.3 ) THEN
!
!        Upper Hessenberg matrix
!
        DO J = 1, N
          DO I = 1, MIN( J+1, M )
            A( I, J ) = A( I, J )*MUL
          enddo
        enddo
!
      ELSE IF( ITYPE.EQ.4 ) THEN
!
!        Lower half of a symmetric band matrix
!
        K3 = KL + 1
        K4 = N + 1
        DO J = 1, N
          DO I = 1, MIN( K3, K4-J )
            A( I, J ) = A( I, J )*MUL
          enddo
        enddo
!
      ELSE IF( ITYPE.EQ.5 ) THEN
!
!        Upper half of a symmetric band matrix
!
        K1 = KU + 2
        K3 = KU + 1
        DO J = 1, N
          DO I = MAX( K1-J, 1 ), K3
            A( I, J ) = A( I, J )*MUL
          enddo
        enddo
!
      ELSE IF( ITYPE.EQ.6 ) THEN
!
!        Band matrix
!
        K1 = KL + KU + 2
        K2 = KL + 1
        K3 = 2*KL + KU + 1
        K4 = KL + KU + 1 + M
        DO J = 1, N
          DO I = MAX( K1-J, K2 ), MIN( K3, K4-J )
            A( I, J ) = A( I, J )*MUL
          enddo
        enddo
!
      END IF
!
      IF( .NOT.DONE )                                                   &
     &GO TO 10
!
      RETURN
!
!     End of DLASCL
!
      END
      SUBROUTINE DLASD0( N, SQRE, D, E, U, LDU, VT, LDVT, SMLSIZ, IWORK,&
     &WORK, INFO )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDU, LDVT, N, SMLSIZ, SQRE
!     ..
!     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), U( LDU, * ), VT( LDVT, * ),    &
     &WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  Using a divide and conquer approach, DLASD0 computes the singular
!  value decomposition (SVD) of a real upper bidiagonal N-by-M
!  matrix B with diagonal D and offdiagonal E, where M = N + SQRE.
!  The algorithm computes orthogonal matrices U and VT such that
!  B = U * S * VT. The singular values S are overwritten on D.
!
!  A related subroutine, DLASDA, computes only the singular values,
!  and optionally, the singular vectors in compact form.
!
!  Arguments
!  =========
!
!  N      (input) INTEGER
!         On entry, the row dimension of the upper bidiagonal matrix.
!         This is also the dimension of the main diagonal array D.
!
!  SQRE   (input) INTEGER
!         Specifies the column dimension of the bidiagonal matrix.
!         = 0: The bidiagonal matrix has column dimension M = N;
!         = 1: The bidiagonal matrix has column dimension M = N+1;
!
!  D      (input/output) DOUBLE PRECISION array, dimension (N)
!         On entry D contains the main diagonal of the bidiagonal
!         matrix.
!         On exit D, if INFO = 0, contains its singular values.
!
!  E      (input) DOUBLE PRECISION array, dimension (M-1)
!         Contains the subdiagonal entries of the bidiagonal matrix.
!         On exit, E has been destroyed.
!
!  U      (output) DOUBLE PRECISION array, dimension at least (LDQ, N)
!         On exit, U contains the left singular vectors.
!
!  LDU    (input) INTEGER
!         On entry, leading dimension of U.
!
!  VT     (output) DOUBLE PRECISION array, dimension at least (LDVT, M)
!         On exit, VT' contains the right singular vectors.
!
!  LDVT   (input) INTEGER
!         On entry, leading dimension of VT.
!
!  SMLSIZ (input) INTEGER
!         On entry, maximum size of the subproblems at the
!         bottom of the computation tree.
!
!  IWORK  (workspace) INTEGER work array.
!         Dimension must be at least (8 * N)
!
!  WORK   (workspace) DOUBLE PRECISION work array.
!         Dimension must be at least (3 * M**2 + 2 * M)
!
!  INFO   (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = 1, an singular value did not converge
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ming Gu and Huan Ren, Computer Science Division, University of
!     California at Berkeley, USA
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, I1, IC, IDXQ, IDXQC, IM1, INODE, ITEMP, IWK,&
     &J, LF, LL, LVL, M, NCC, ND, NDB1, NDIML, NDIMR,                   &
     &NL, NLF, NLP1, NLVL, NR, NRF, NRP1, SQREI
      DOUBLE PRECISION   ALPHA, BETA
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASD1, DLASDQ, DLASDT, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( N.LT.0 ) THEN
        INFO = -1
      ELSE IF( ( SQRE.LT.0 ) .OR. ( SQRE.GT.1 ) ) THEN
        INFO = -2
      END IF
!
      M = N + SQRE
!
      IF( LDU.LT.N ) THEN
        INFO = -6
      ELSE IF( LDVT.LT.M ) THEN
        INFO = -8
      ELSE IF( SMLSIZ.LT.3 ) THEN
        INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DLASD0', -INFO )
        RETURN
      END IF
!
!     If the input matrix is too small, call DLASDQ to find the SVD.
!
      IF( N.LE.SMLSIZ ) THEN
        CALL DLASDQ( 'U', SQRE, N, M, N, 0, D, E, VT, LDVT, U, LDU, U,  &
     &LDU, WORK, INFO )
        RETURN
      END IF
!
!     Set up the computation tree.
!
      INODE = 1
      NDIML = INODE + N
      NDIMR = NDIML + N
      IDXQ = NDIMR + N
      IWK = IDXQ + N
      CALL DLASDT( N, NLVL, ND, IWORK( INODE ), IWORK( NDIML ),         &
     &IWORK( NDIMR ), SMLSIZ )
!
!     For the nodes on bottom level of the tree, solve
!     their subproblems by DLASDQ.
!
      NDB1 = ( ND+1 ) / 2
      NCC = 0
      DO I = NDB1, ND
!
!     IC : center row of each node
!     NL : number of rows of left  subproblem
!     NR : number of rows of right subproblem
!     NLF: starting row of the left   subproblem
!     NRF: starting row of the right  subproblem
!
        I1 = I - 1
        IC = IWORK( INODE+I1 )
        NL = IWORK( NDIML+I1 )
        NLP1 = NL + 1
        NR = IWORK( NDIMR+I1 )
        NRP1 = NR + 1
        NLF = IC - NL
        NRF = IC + 1
        SQREI = 1
        CALL DLASDQ( 'U', SQREI, NL, NLP1, NL, NCC, D( NLF ), E( NLF ), &
     &VT( NLF, NLF ), LDVT, U( NLF, NLF ), LDU,                         &
     &U( NLF, NLF ), LDU, WORK, INFO )
        IF( INFO.NE.0 ) THEN
          RETURN
        END IF
        ITEMP = IDXQ + NLF - 2
        DO J = 1, NL
          IWORK( ITEMP+J ) = J
        enddo
        IF( I.EQ.ND ) THEN
          SQREI = SQRE
        ELSE
          SQREI = 1
        END IF
        NRP1 = NR + SQREI
        CALL DLASDQ( 'U', SQREI, NR, NRP1, NR, NCC, D( NRF ), E( NRF ), &
     &VT( NRF, NRF ), LDVT, U( NRF, NRF ), LDU,                         &
     &U( NRF, NRF ), LDU, WORK, INFO )
        IF( INFO.NE.0 ) THEN
          RETURN
        END IF
        ITEMP = IDXQ + IC
        DO J = 1, NR
          IWORK( ITEMP+J-1 ) = J
        enddo
      enddo
!
!     Now conquer each subproblem bottom-up.
!
      DO LVL = NLVL, 1, -1
!
!        Find the first node LF and last node LL on the
!        current level LVL.
!
        IF( LVL.EQ.1 ) THEN
          LF = 1
          LL = 1
        ELSE
          LF = 2**( LVL-1 )
          LL = 2*LF - 1
        END IF
        DO I = LF, LL
          IM1 = I - 1
          IC = IWORK( INODE+IM1 )
          NL = IWORK( NDIML+IM1 )
          NR = IWORK( NDIMR+IM1 )
          NLF = IC - NL
          IF( ( SQRE.EQ.0 ) .AND. ( I.EQ.LL ) ) THEN
            SQREI = SQRE
          ELSE
            SQREI = 1
          END IF
          IDXQC = IDXQ + NLF - 1
          ALPHA = D( IC )
          BETA = E( IC )
          CALL DLASD1( NL, NR, SQREI, D( NLF ), ALPHA, BETA,            &
     &U( NLF, NLF ), LDU, VT( NLF, NLF ), LDVT,                         &
     &IWORK( IDXQC ), IWORK( IWK ), WORK, INFO )
          IF( INFO.NE.0 ) THEN
            RETURN
          END IF
        enddo
      enddo
!
      RETURN
!
!     End of DLASD0
!
      END
      SUBROUTINE DLASD1( NL, NR, SQRE, D, ALPHA, BETA, U, LDU, VT, LDVT,&
     &IDXQ, IWORK, WORK, INFO )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDU, LDVT, NL, NR, SQRE
      DOUBLE PRECISION   ALPHA, BETA
!     ..
!     .. Array Arguments ..
      INTEGER            IDXQ( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), U( LDU, * ), VT( LDVT, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASD1 computes the SVD of an upper bidiagonal N-by-M matrix B,
!  where N = NL + NR + 1 and M = N + SQRE. DLASD1 is called from DLASD0.
!
!  A related subroutine DLASD7 handles the case in which the singular
!  values (and the singular vectors in factored form) are desired.
!
!  DLASD1 computes the SVD as follows:
!
!                ( D1(in)  0    0     0 )
!    B = U(in) * (   Z1'   a   Z2'    b ) * VT(in)
!                (   0     0   D2(in) 0 )
!
!      = U(out) * ( D(out) 0) * VT(out)
!
!  where Z' = (Z1' a Z2' b) = u' VT', and u is a vector of dimension M
!  with ALPHA and BETA in the NL+1 and NL+2 th entries and zeros
!  elsewhere; and the entry b is empty if SQRE = 0.
!
!  The left singular vectors of the original matrix are stored in U, and
!  the transpose of the right singular vectors are stored in VT, and the
!  singular values are in D.  The algorithm consists of three stages:
!
!     The first stage consists of deflating the size of the problem
!     when there are multiple singular values or when there are zeros in
!     the Z vector.  For each such occurence the dimension of the
!     secular equation problem is reduced by one.  This stage is
!     performed by the routine DLASD2.
!
!     The second stage consists of calculating the updated
!     singular values. This is done by finding the square roots of the
!     roots of the secular equation via the routine DLASD4 (as called
!     by DLASD3). This routine also calculates the singular vectors of
!     the current problem.
!
!     The final stage consists of computing the updated singular vectors
!     directly using the updated singular values.  The singular vectors
!     for the current problem are multiplied with the singular vectors
!     from the overall problem.
!
!  Arguments
!  =========
!
!  NL     (input) INTEGER
!         The row dimension of the upper block.  NL >= 1.
!
!  NR     (input) INTEGER
!         The row dimension of the lower block.  NR >= 1.
!
!  SQRE   (input) INTEGER
!         = 0: the lower block is an NR-by-NR square matrix.
!         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.
!
!         The bidiagonal matrix has row dimension N = NL + NR + 1,
!         and column dimension M = N + SQRE.
!
!  D      (input/output) DOUBLE PRECISION array,
!                        dimension (N = NL+NR+1).
!         On entry D(1:NL,1:NL) contains the singular values of the
!         upper block; and D(NL+2:N) contains the singular values of
!         the lower block. On exit D(1:N) contains the singular values
!         of the modified matrix.
!
!  ALPHA  (input/output) DOUBLE PRECISION
!         Contains the diagonal element associated with the added row.
!
!  BETA   (input/output) DOUBLE PRECISION
!         Contains the off-diagonal element associated with the added
!         row.
!
!  U      (input/output) DOUBLE PRECISION array, dimension(LDU,N)
!         On entry U(1:NL, 1:NL) contains the left singular vectors of
!         the upper block; U(NL+2:N, NL+2:N) contains the left singular
!         vectors of the lower block. On exit U contains the left
!         singular vectors of the bidiagonal matrix.
!
!  LDU    (input) INTEGER
!         The leading dimension of the array U.  LDU >= max( 1, N ).
!
!  VT     (input/output) DOUBLE PRECISION array, dimension(LDVT,M)
!         where M = N + SQRE.
!         On entry VT(1:NL+1, 1:NL+1)' contains the right singular
!         vectors of the upper block; VT(NL+2:M, NL+2:M)' contains
!         the right singular vectors of the lower block. On exit
!         VT' contains the right singular vectors of the
!         bidiagonal matrix.
!
!  LDVT   (input) INTEGER
!         The leading dimension of the array VT.  LDVT >= max( 1, M ).
!
!  IDXQ  (output) INTEGER array, dimension(N)
!         This contains the permutation which will reintegrate the
!         subproblem just solved back into sorted order, i.e.
!         D( IDXQ( I = 1, N ) ) will be in ascending order.
!
!  IWORK  (workspace) INTEGER array, dimension( 4 * N )
!
!  WORK   (workspace) DOUBLE PRECISION array, dimension( 3*M**2 + 2*M )
!
!  INFO   (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = 1, an singular value did not converge
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ming Gu and Huan Ren, Computer Science Division, University of
!     California at Berkeley, USA
!
!  =====================================================================
!
!     .. Parameters ..
!
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            COLTYP, I, IDX, IDXC, IDXP, IQ, ISIGMA, IU2,   &
     &IVT2, IZ, K, LDQ, LDU2, LDVT2, M, N, N1, N2
      DOUBLE PRECISION   ORGNRM
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAMRG, DLASCL, DLASD2, DLASD3, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( NL.LT.1 ) THEN
        INFO = -1
      ELSE IF( NR.LT.1 ) THEN
        INFO = -2
      ELSE IF( ( SQRE.LT.0 ) .OR. ( SQRE.GT.1 ) ) THEN
        INFO = -3
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DLASD1', -INFO )
        RETURN
      END IF
!
      N = NL + NR + 1
      M = N + SQRE
!
!     The following values are for bookkeeping purposes only.  They are
!     integer pointers which indicate the portion of the workspace
!     used by a particular array in DLASD2 and DLASD3.
!
      LDU2 = N
      LDVT2 = M
!
      IZ = 1
      ISIGMA = IZ + M
      IU2 = ISIGMA + N
      IVT2 = IU2 + LDU2*N
      IQ = IVT2 + LDVT2*M
!
      IDX = 1
      IDXC = IDX + N
      COLTYP = IDXC + N
      IDXP = COLTYP + N
!
!     Scale.
!
      ORGNRM = MAX( ABS( ALPHA ), ABS( BETA ) )
      D( NL+1 ) = ZERO
      DO I = 1, N
        IF( ABS( D( I ) ).GT.ORGNRM ) THEN
          ORGNRM = ABS( D( I ) )
        END IF
      enddo
      CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, N, 1, D, N, INFO )
      ALPHA = ALPHA / ORGNRM
      BETA = BETA / ORGNRM
!
!     Deflate singular values.
!
      CALL DLASD2( NL, NR, SQRE, K, D, WORK( IZ ), ALPHA, BETA, U, LDU, &
     &VT, LDVT, WORK( ISIGMA ), WORK( IU2 ), LDU2,                      &
     &WORK( IVT2 ), LDVT2, IWORK( IDXP ), IWORK( IDX ),                 &
     &IWORK( IDXC ), IDXQ, IWORK( COLTYP ), INFO )
!
!     Solve Secular Equation and update singular vectors.
!
      LDQ = K
      CALL DLASD3( NL, NR, SQRE, K, D, WORK( IQ ), LDQ, WORK( ISIGMA ), &
     &U, LDU, WORK( IU2 ), LDU2, VT, LDVT, WORK( IVT2 ),                &
     &LDVT2, IWORK( IDXC ), IWORK( COLTYP ), WORK( IZ ),                &
     &INFO )
      IF( INFO.NE.0 ) THEN
        RETURN
      END IF
!
!     Unscale.
!
      CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO )
!
!     Prepare the IDXQ sorting permutation.
!
      N1 = K
      N2 = N - K
      CALL DLAMRG( N1, N2, D, 1, -1, IDXQ )
!
      RETURN
!
!     End of DLASD1
!
      END
      SUBROUTINE DLASD2( NL, NR, SQRE, K, D, Z, ALPHA, BETA, U, LDU, VT,&
     &LDVT, DSIGMA, U2, LDU2, VT2, LDVT2, IDXP, IDX,                    &
     &IDXC, IDXQ, COLTYP, INFO )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDU, LDU2, LDVT, LDVT2, NL, NR, SQRE
      DOUBLE PRECISION   ALPHA, BETA
!     ..
!     .. Array Arguments ..
      INTEGER            COLTYP( * ), IDX( * ), IDXC( * ), IDXP( * ),   &
     &IDXQ( * )
      DOUBLE PRECISION   D( * ), DSIGMA( * ), U( LDU, * ),              &
     &U2( LDU2, * ), VT( LDVT, * ), VT2( LDVT2, * ),                    &
     &Z( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASD2 merges the two sets of singular values together into a single
!  sorted set.  Then it tries to deflate the size of the problem.
!  There are two ways in which deflation can occur:  when two or more
!  singular values are close together or if there is a tiny entry in the
!  Z vector.  For each such occurrence the order of the related secular
!  equation problem is reduced by one.
!
!  DLASD2 is called from DLASD1.
!
!  Arguments
!  =========
!
!  NL     (input) INTEGER
!         The row dimension of the upper block.  NL >= 1.
!
!  NR     (input) INTEGER
!         The row dimension of the lower block.  NR >= 1.
!
!  SQRE   (input) INTEGER
!         = 0: the lower block is an NR-by-NR square matrix.
!         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.
!
!         The bidiagonal matrix has N = NL + NR + 1 rows and
!         M = N + SQRE >= N columns.
!
!  K      (output) INTEGER
!         Contains the dimension of the non-deflated matrix,
!         This is the order of the related secular equation. 1 <= K <=N.
!
!  D      (input/output) DOUBLE PRECISION array, dimension(N)
!         On entry D contains the singular values of the two submatrices
!         to be combined.  On exit D contains the trailing (N-K) updated
!         singular values (those which were deflated) sorted into
!         increasing order.
!
!  Z      (output) DOUBLE PRECISION array, dimension(N)
!         On exit Z contains the updating row vector in the secular
!         equation.
!
!  ALPHA  (input) DOUBLE PRECISION
!         Contains the diagonal element associated with the added row.
!
!  BETA   (input) DOUBLE PRECISION
!         Contains the off-diagonal element associated with the added
!         row.
!
!  U      (input/output) DOUBLE PRECISION array, dimension(LDU,N)
!         On entry U contains the left singular vectors of two
!         submatrices in the two square blocks with corners at (1,1),
!         (NL, NL), and (NL+2, NL+2), (N,N).
!         On exit U contains the trailing (N-K) updated left singular
!         vectors (those which were deflated) in its last N-K columns.
!
!  LDU    (input) INTEGER
!         The leading dimension of the array U.  LDU >= N.
!
!  VT     (input/output) DOUBLE PRECISION array, dimension(LDVT,M)
!         On entry VT' contains the right singular vectors of two
!         submatrices in the two square blocks with corners at (1,1),
!         (NL+1, NL+1), and (NL+2, NL+2), (M,M).
!         On exit VT' contains the trailing (N-K) updated right singular
!         vectors (those which were deflated) in its last N-K columns.
!         In case SQRE =1, the last row of VT spans the right null
!         space.
!
!  LDVT   (input) INTEGER
!         The leading dimension of the array VT.  LDVT >= M.
!
!  DSIGMA (output) DOUBLE PRECISION array, dimension (N)
!         Contains a copy of the diagonal elements (K-1 singular values
!         and one zero) in the secular equation.
!
!  U2     (output) DOUBLE PRECISION array, dimension(LDU2,N)
!         Contains a copy of the first K-1 left singular vectors which
!         will be used by DLASD3 in a matrix multiply (DGEMM) to solve
!         for the new left singular vectors. U2 is arranged into four
!         blocks. The first block contains a column with 1 at NL+1 and
!         zero everywhere else; the second block contains non-zero
!         entries only at and above NL; the third contains non-zero
!         entries only below NL+1; and the fourth is dense.
!
!  LDU2   (input) INTEGER
!         The leading dimension of the array U2.  LDU2 >= N.
!
!  VT2    (output) DOUBLE PRECISION array, dimension(LDVT2,N)
!         VT2' contains a copy of the first K right singular vectors
!         which will be used by DLASD3 in a matrix multiply (DGEMM) to
!         solve for the new right singular vectors. VT2 is arranged into
!         three blocks. The first block contains a row that corresponds
!         to the special 0 diagonal element in SIGMA; the second block
!         contains non-zeros only at and before NL +1; the third block
!         contains non-zeros only at and after  NL +2.
!
!  LDVT2  (input) INTEGER
!         The leading dimension of the array VT2.  LDVT2 >= M.
!
!  IDXP   (workspace) INTEGER array dimension(N)
!         This will contain the permutation used to place deflated
!         values of D at the end of the array. On output IDXP(2:K)
!         points to the nondeflated D-values and IDXP(K+1:N)
!         points to the deflated singular values.
!
!  IDX    (workspace) INTEGER array dimension(N)
!         This will contain the permutation used to sort the contents of
!         D into ascending order.
!
!  IDXC   (output) INTEGER array dimension(N)
!         This will contain the permutation used to arrange the columns
!         of the deflated U matrix into three groups:  the first group
!         contains non-zero entries only at and above NL, the second
!         contains non-zero entries only below NL+2, and the third is
!         dense.
!
!  IDXQ   (input/output) INTEGER array dimension(N)
!         This contains the permutation which separately sorts the two
!         sub-problems in D into ascending order.  Note that entries in
!         the first hlaf of this permutation must first be moved one
!         position backward; and entries in the second half
!         must first have NL+1 added to their values.
!
!  COLTYP (workspace/output) INTEGER array dimension(N)
!         As workspace, this will contain a label which will indicate
!         which of the following types a column in the U2 matrix or a
!         row in the VT2 matrix is:
!         1 : non-zero in the upper half only
!         2 : non-zero in the lower half only
!         3 : dense
!         4 : deflated
!
!         On exit, it is an array of dimension 4, with COLTYP(I) being
!         the dimension of the I-th type columns.
!
!  INFO   (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ming Gu and Huan Ren, Computer Science Division, University of
!     California at Berkeley, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, EIGHT
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0,   &
     &EIGHT = 8.0D+0 )
!     ..
!     .. Local Arrays ..
      INTEGER            CTOT( 4 ), PSM( 4 )
!     ..
!     .. Local Scalars ..
      INTEGER            CT, I, IDXI, IDXJ, IDXJP, J, JP, JPREV, K2, M, &
     &N, NLP1, NLP2
      DOUBLE PRECISION   C, EPS, HLFTOL, S, TAU, TOL, Z1
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2
      EXTERNAL           DLAMCH, DLAPY2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPL, DLACPY, DLAMRG, DLASET, DROT, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( NL.LT.1 ) THEN
        INFO = -1
      ELSE IF( NR.LT.1 ) THEN
        INFO = -2
      ELSE IF( ( SQRE.NE.1 ) .AND. ( SQRE.NE.0 ) ) THEN
        INFO = -3
      END IF
!
      N = NL + NR + 1
      M = N + SQRE
!
      IF( LDU.LT.N ) THEN
        INFO = -10
      ELSE IF( LDVT.LT.M ) THEN
        INFO = -12
      ELSE IF( LDU2.LT.N ) THEN
        INFO = -15
      ELSE IF( LDVT2.LT.M ) THEN
        INFO = -17
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DLASD2', -INFO )
        RETURN
      END IF
!
      NLP1 = NL + 1
      NLP2 = NL + 2
!
!     Generate the first part of the vector Z; and move the singular
!     values in the first part of D one position backward.
!
      Z1 = ALPHA*VT( NLP1, NLP1 )
      Z( 1 ) = Z1
      DO I = NL, 1, -1
        Z( I+1 ) = ALPHA*VT( I, NLP1 )
        D( I+1 ) = D( I )
        IDXQ( I+1 ) = IDXQ( I ) + 1
      enddo
!
!     Generate the second part of the vector Z.
!
      DO I = NLP2, M
        Z( I ) = BETA*VT( I, NLP2 )
      enddo
!
!     Initialize some reference arrays.
!
      DO I = 2, NLP1
        COLTYP( I ) = 1
      enddo
      DO I = NLP2, N
        COLTYP( I ) = 2
      enddo
!
!     Sort the singular values into increasing order
!
      DO I = NLP2, N
        IDXQ( I ) = IDXQ( I ) + NLP1
      enddo
!
!     DSIGMA, IDXC, IDXC, and the first column of U2
!     are used as storage space.
!
      DO I = 2, N
        DSIGMA( I ) = D( IDXQ( I ) )
        U2( I, 1 ) = Z( IDXQ( I ) )
        IDXC( I ) = COLTYP( IDXQ( I ) )
      enddo
!
      CALL DLAMRG( NL, NR, DSIGMA( 2 ), 1, 1, IDX( 2 ) )
!
      DO I = 2, N
        IDXI = 1 + IDX( I )
        D( I ) = DSIGMA( IDXI )
        Z( I ) = U2( IDXI, 1 )
        COLTYP( I ) = IDXC( IDXI )
      enddo
!
!     Calculate the allowable deflation tolerance
!
      EPS = DLAMCH( 'Epsilon' )
      TOL = MAX( ABS( ALPHA ), ABS( BETA ) )
      TOL = EIGHT*EPS*MAX( ABS( D( N ) ), TOL )
!
!     There are 2 kinds of deflation -- first a value in the z-vector
!     is small, second two (or more) singular values are very close
!     together (their difference is small).
!
!     If the value in the z-vector is small, we simply permute the
!     array so that the corresponding singular value is moved to the
!     end.
!
!     If two values in the D-vector are close, we perform a two-sided
!     rotation designed to make one of the corresponding z-vector
!     entries zero, and then permute the array so that the deflated
!     singular value is moved to the end.
!
!     If there are multiple singular values then the problem deflates.
!     Here the number of equal singular values are found.  As each equal
!     singular value is found, an elementary reflector is computed to
!     rotate the corresponding singular subspace so that the
!     corresponding components of Z are zero in this new basis.
!
      K = 1
      K2 = N + 1
      DO J = 2, N
        IF( ABS( Z( J ) ).LE.TOL ) THEN
!
!           Deflate due to small z component.
!
          K2 = K2 - 1
          IDXP( K2 ) = J
          COLTYP( J ) = 4
          IF( J.EQ.N )                                                  &
     &GO TO 120
        ELSE
          JPREV = J
          GO TO 90
        END IF
      enddo
   90 CONTINUE
      J = JPREV
  100 CONTINUE
      J = J + 1
      IF( J.GT.N )                                                      &
     &GO TO 110
      IF( ABS( Z( J ) ).LE.TOL ) THEN
!
!        Deflate due to small z component.
!
        K2 = K2 - 1
        IDXP( K2 ) = J
        COLTYP( J ) = 4
      ELSE
!
!        Check if singular values are close enough to allow deflation.
!
        IF( ABS( D( J )-D( JPREV ) ).LE.TOL ) THEN
!
!           Deflation is possible.
!
          S = Z( JPREV )
          C = Z( J )
!
!           Find sqrt(a**2+b**2) without overflow or
!           destructive underflow.
!
          TAU = DLAPY2( C, S )
          C = C / TAU
          S = -S / TAU
          Z( J ) = TAU
          Z( JPREV ) = ZERO
!
!           Apply back the Givens rotation to the left and right
!           singular vector matrices.
!
          IDXJP = IDXQ( IDX( JPREV )+1 )
          IDXJ = IDXQ( IDX( J )+1 )
          IF( IDXJP.LE.NLP1 ) THEN
            IDXJP = IDXJP - 1
          END IF
          IF( IDXJ.LE.NLP1 ) THEN
            IDXJ = IDXJ - 1
          END IF
          CALL DROT( N, U( 1, IDXJP ), 1, U( 1, IDXJ ), 1, C, S )
          CALL DROT( M, VT( IDXJP, 1 ), LDVT, VT( IDXJ, 1 ), LDVT, C,   &
     &S )
          IF( COLTYP( J ).NE.COLTYP( JPREV ) ) THEN
            COLTYP( J ) = 3
          END IF
          COLTYP( JPREV ) = 4
          K2 = K2 - 1
          IDXP( K2 ) = JPREV
          JPREV = J
        ELSE
          K = K + 1
          U2( K, 1 ) = Z( JPREV )
          DSIGMA( K ) = D( JPREV )
          IDXP( K ) = JPREV
          JPREV = J
        END IF
      END IF
      GO TO 100
  110 CONTINUE
!
!     Record the last singular value.
!
      K = K + 1
      U2( K, 1 ) = Z( JPREV )
      DSIGMA( K ) = D( JPREV )
      IDXP( K ) = JPREV
!
  120 CONTINUE
!
!     Count up the total number of the various types of columns, then
!     form a permutation which positions the four column types into
!     four groups of uniform structure (although one or more of these
!     groups may be empty).
!
      DO J = 1, 4
        CTOT( J ) = 0
      enddo
      DO J = 2, N
        CT = COLTYP( J )
        CTOT( CT ) = CTOT( CT ) + 1
      enddo
!
!     PSM(*) = Position in SubMatrix (of types 1 through 4)
!
      PSM( 1 ) = 2
      PSM( 2 ) = 2 + CTOT( 1 )
      PSM( 3 ) = PSM( 2 ) + CTOT( 2 )
      PSM( 4 ) = PSM( 3 ) + CTOT( 3 )
!
!     Fill out the IDXC array so that the permutation which it induces
!     will place all type-1 columns first, all type-2 columns next,
!     then all type-3's, and finally all type-4's, starting from the
!     second column. This applies similarly to the rows of VT.
!
      DO J = 2, N
        JP = IDXP( J )
        CT = COLTYP( JP )
        IDXC( PSM( CT ) ) = J
        PSM( CT ) = PSM( CT ) + 1
      enddo
!
!     Sort the singular values and corresponding singular vectors into
!     DSIGMA, U2, and VT2 respectively.  The singular values/vectors
!     which were not deflated go into the first K slots of DSIGMA, U2,
!     and VT2 respectively, while those which were deflated go into the
!     last N - K slots, except that the first column/row will be treated
!     separately.
!
      DO J = 2, N
        JP = IDXP( J )
        DSIGMA( J ) = D( JP )
        IDXJ = IDXQ( IDX( IDXP( IDXC( J ) ) )+1 )
        IF( IDXJ.LE.NLP1 ) THEN
          IDXJ = IDXJ - 1
        END IF
        CALL DCOPL( N, U( 1, IDXJ ), 1, U2( 1, J ), 1 )
        CALL DCOPL( M, VT( IDXJ, 1 ), LDVT, VT2( J, 1 ), LDVT2 )
      enddo
!
!     Determine DSIGMA(1), DSIGMA(2) and Z(1)
!
      DSIGMA( 1 ) = ZERO
      HLFTOL = TOL / TWO
      IF( ABS( DSIGMA( 2 ) ).LE.HLFTOL )                                &
     &DSIGMA( 2 ) = HLFTOL
      IF( M.GT.N ) THEN
        Z( 1 ) = DLAPY2( Z1, Z( M ) )
        IF( Z( 1 ).LE.TOL ) THEN
          C = ONE
          S = ZERO
          Z( 1 ) = TOL
        ELSE
          C = Z1 / Z( 1 )
          S = Z( M ) / Z( 1 )
        END IF
      ELSE
        IF( ABS( Z1 ).LE.TOL ) THEN
          Z( 1 ) = TOL
        ELSE
          Z( 1 ) = Z1
        END IF
      END IF
!
!     Move the rest of the updating row to Z.
!
      CALL DCOPL( K-1, U2( 2, 1 ), 1, Z( 2 ), 1 )
!
!     Determine the first column of U2, the first row of VT2 and the
!     last row of VT.
!
      CALL DLASET( 'A', N, 1, ZERO, ZERO, U2, LDU2 )
      U2( NLP1, 1 ) = ONE
      IF( M.GT.N ) THEN
        DO I = 1, NLP1
          VT( M, I ) = -S*VT( NLP1, I )
          VT2( 1, I ) = C*VT( NLP1, I )
        enddo
        DO I = NLP2, M
          VT2( 1, I ) = S*VT( M, I )
          VT( M, I ) = C*VT( M, I )
        enddo
      ELSE
        CALL DCOPL( M, VT( NLP1, 1 ), LDVT, VT2( 1, 1 ), LDVT2 )
      END IF
      IF( M.GT.N ) THEN
        CALL DCOPL( M, VT( M, 1 ), LDVT, VT2( M, 1 ), LDVT2 )
      END IF
!
!     The deflated singular values and their corresponding vectors go
!     into the back of D, U, and V respectively.
!
      IF( N.GT.K ) THEN
        CALL DCOPL( N-K, DSIGMA( K+1 ), 1, D( K+1 ), 1 )
        CALL DLACPY( 'A', N, N-K, U2( 1, K+1 ), LDU2, U( 1, K+1 ),      &
     &LDU )
        CALL DLACPY( 'A', N-K, M, VT2( K+1, 1 ), LDVT2, VT( K+1, 1 ),   &
     &LDVT )
      END IF
!
!     Copy CTOT into COLTYP for referencing in DLASD3.
!
      DO J = 1, 4
        COLTYP( J ) = CTOT( J )
      enddo
!
      RETURN
!
!     End of DLASD2
!
      END
      SUBROUTINE DLASD3( NL, NR, SQRE, K, D, Q, LDQ, DSIGMA, U, LDU, U2,&
     &LDU2, VT, LDVT, VT2, LDVT2, IDXC, CTOT, Z,                        &
     &INFO )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDQ, LDU, LDU2, LDVT, LDVT2, NL, NR,  &
     &SQRE
!     ..
!     .. Array Arguments ..
      INTEGER            CTOT( * ), IDXC( * )
      DOUBLE PRECISION   D( * ), DSIGMA( * ), Q( LDQ, * ), U( LDU, * ), &
     &U2( LDU2, * ), VT( LDVT, * ), VT2( LDVT2, * ),                    &
     &Z( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASD3 finds all the square roots of the roots of the secular
!  equation, as defined by the values in D and Z.  It makes the
!  appropriate calls to DLASD4 and then updates the singular
!  vectors by matrix multiplication.
!
!  This code makes very mild assumptions about floating point
!  arithmetic. It will work on machines with a guard digit in
!  add/subtract, or on those binary machines without guard digits
!  which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2.
!  It could conceivably fail on hexadecimal or decimal machines
!  without guard digits, but we know of none.
!
!  DLASD3 is called from DLASD1.
!
!  Arguments
!  =========
!
!  NL     (input) INTEGER
!         The row dimension of the upper block.  NL >= 1.
!
!  NR     (input) INTEGER
!         The row dimension of the lower block.  NR >= 1.
!
!  SQRE   (input) INTEGER
!         = 0: the lower block is an NR-by-NR square matrix.
!         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.
!
!         The bidiagonal matrix has N = NL + NR + 1 rows and
!         M = N + SQRE >= N columns.
!
!  K      (input) INTEGER
!         The size of the secular equation, 1 =< K = < N.
!
!  D      (output) DOUBLE PRECISION array, dimension(K)
!         On exit the square roots of the roots of the secular equation,
!         in ascending order.
!
!  Q      (workspace) DOUBLE PRECISION array,
!                     dimension at least (LDQ,K).
!
!  LDQ    (input) INTEGER
!         The leading dimension of the array Q.  LDQ >= K.
!
!  DSIGMA (input) DOUBLE PRECISION array, dimension(K)
!         The first K elements of this array contain the old roots
!         of the deflated updating problem.  These are the poles
!         of the secular equation.
!
!  U      (output) DOUBLE PRECISION array, dimension (LDU, N)
!         The last N - K columns of this matrix contain the deflated
!         left singular vectors.
!
!  LDU    (input) INTEGER
!         The leading dimension of the array U.  LDU >= N.
!
!  U2     (input/output) DOUBLE PRECISION array, dimension (LDU2, N)
!         The first K columns of this matrix contain the non-deflated
!         left singular vectors for the split problem.
!
!  LDU2   (input) INTEGER
!         The leading dimension of the array U2.  LDU2 >= N.
!
!  VT     (output) DOUBLE PRECISION array, dimension (LDVT, M)
!         The last M - K columns of VT' contain the deflated
!         right singular vectors.
!
!  LDVT   (input) INTEGER
!         The leading dimension of the array VT.  LDVT >= N.
!
!  VT2    (input/output) DOUBLE PRECISION array, dimension (LDVT2, N)
!         The first K columns of VT2' contain the non-deflated
!         right singular vectors for the split problem.
!
!  LDVT2  (input) INTEGER
!         The leading dimension of the array VT2.  LDVT2 >= N.
!
!  IDXC   (input) INTEGER array, dimension ( N )
!         The permutation used to arrange the columns of U (and rows of
!         VT) into three groups:  the first group contains non-zero
!         entries only at and above (or before) NL +1; the second
!         contains non-zero entries only at and below (or after) NL+2;
!         and the third is dense. The first column of U and the row of
!         VT are treated separately, however.
!
!         The rows of the singular vectors found by DLASD4
!         must be likewise permuted before the matrix multiplies can
!         take place.
!
!  CTOT   (input) INTEGER array, dimension ( 4 )
!         A count of the total number of the various types of columns
!         in U (or rows in VT), as described in IDXC. The fourth column
!         type is any column which has been deflated.
!
!  Z      (input) DOUBLE PRECISION array, dimension (K)
!         The first K elements of this array contain the components
!         of the deflation-adjusted updating row vector.
!
!  INFO   (output) INTEGER
!         = 0:  successful exit.
!         < 0:  if INFO = -i, the i-th argument had an illegal value.
!         > 0:  if INFO = 1, an singular value did not converge
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ming Gu and Huan Ren, Computer Science Division, University of
!     California at Berkeley, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO, NEGONE
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0,                 &
     &NEGONE = -1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            CTEMP, I, J, JC, KTEMP, M, N, NLP1, NLP2, NRP1
      DOUBLE PRECISION   RHO, TEMP
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3, DNRM2
      EXTERNAL           DLAMC3, DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPL, DGEMM, DLACPY, DLASCL, DLASD4, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( NL.LT.1 ) THEN
        INFO = -1
      ELSE IF( NR.LT.1 ) THEN
        INFO = -2
      ELSE IF( ( SQRE.NE.1 ) .AND. ( SQRE.NE.0 ) ) THEN
        INFO = -3
      END IF
!
      N = NL + NR + 1
      M = N + SQRE
      NLP1 = NL + 1
      NLP2 = NL + 2
!
      IF( ( K.LT.1 ) .OR. ( K.GT.N ) ) THEN
        INFO = -4
      ELSE IF( LDQ.LT.K ) THEN
        INFO = -7
      ELSE IF( LDU.LT.N ) THEN
        INFO = -10
      ELSE IF( LDU2.LT.N ) THEN
        INFO = -12
      ELSE IF( LDVT.LT.M ) THEN
        INFO = -14
      ELSE IF( LDVT2.LT.M ) THEN
        INFO = -16
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DLASD3', -INFO )
        RETURN
      END IF
!
!     Quick return if possible
!
      IF( K.EQ.1 ) THEN
        D( 1 ) = ABS( Z( 1 ) )
        CALL DCOPL( M, VT2( 1, 1 ), LDVT2, VT( 1, 1 ), LDVT )
        IF( Z( 1 ).GT.ZERO ) THEN
          CALL DCOPL( N, U2( 1, 1 ), 1, U( 1, 1 ), 1 )
        ELSE
          DO I = 1, N
            U( I, 1 ) = -U2( I, 1 )
          enddo
        END IF
        RETURN
      END IF
!
!     Modify values DSIGMA(i) to make sure all DSIGMA(i)-DSIGMA(j) can
!     be computed with high relative accuracy (barring over/underflow).
!     This is a problem on machines without a guard digit in
!     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
!     The following code replaces DSIGMA(I) by 2*DSIGMA(I)-DSIGMA(I),
!     which on any of these machines zeros out the bottommost
!     bit of DSIGMA(I) if it is 1; this makes the subsequent
!     subtractions DSIGMA(I)-DSIGMA(J) unproblematic when cancellation
!     occurs. On binary machines with a guard digit (almost all
!     machines) it does not change DSIGMA(I) at all. On hexadecimal
!     and decimal machines with a guard digit, it slightly
!     changes the bottommost bits of DSIGMA(I). It does not account
!     for hexadecimal or decimal machines without guard digits
!     (we know of none). We use a subroutine call to compute
!     2*DSIGMA(I) to prevent optimizing compilers from eliminating
!     this code.
!
      DO I = 1, K
        DSIGMA( I ) = DLAMC3( DSIGMA( I ), DSIGMA( I ) ) - DSIGMA( I )
      enddo
!
!     Keep a copy of Z.
!
      CALL DCOPL( K, Z, 1, Q, 1 )
!
!     Normalize Z.
!
      RHO = DNRM2( K, Z, 1 )
      CALL DLASCL( 'G', 0, 0, RHO, ONE, K, 1, Z, K, INFO )
      RHO = RHO*RHO
!
!     Find the new singular values.
!
      DO J = 1, K
        CALL DLASD4( K, J, DSIGMA, Z, U( 1, J ), RHO, D( J ),           &
     &VT( 1, J ), INFO )
!
!        If the zero finder fails, the computation is terminated.
!
        IF( INFO.NE.0 ) THEN
          RETURN
        END IF
      enddo
!
!     Compute updated Z.
!
      DO I = 1, K
        Z( I ) = U( I, K )*VT( I, K )
        DO J = 1, I - 1
          Z( I ) = Z( I )*( U( I, J )*VT( I, J ) /                      &
     &( DSIGMA( I )-DSIGMA( J ) ) /                                     &
     &( DSIGMA( I )+DSIGMA( J ) ) )
        enddo
        DO J = I, K - 1
          Z( I ) = Z( I )*( U( I, J )*VT( I, J ) /                      &
     &( DSIGMA( I )-DSIGMA( J+1 ) ) /                                   &
     &( DSIGMA( I )+DSIGMA( J+1 ) ) )
        enddo
        Z( I ) = SIGN( SQRT( ABS( Z( I ) ) ), Q( I, 1 ) )
      enddo
!
!     Compute left singular vectors of the modified diagonal matrix,
!     and store related information for the right singular vectors.
!
      DO I = 1, K
        VT( 1, I ) = Z( 1 ) / U( 1, I ) / VT( 1, I )
        U( 1, I ) = NEGONE
        DO J = 2, K
          VT( J, I ) = Z( J ) / U( J, I ) / VT( J, I )
          U( J, I ) = DSIGMA( J )*VT( J, I )
        enddo
        TEMP = DNRM2( K, U( 1, I ), 1 )
        Q( 1, I ) = U( 1, I ) / TEMP
        DO J = 2, K
          JC = IDXC( J )
          Q( J, I ) = U( JC, I ) / TEMP
        enddo
      enddo
!
!     Update the left singular vector matrix.
!
      IF( K.EQ.2 ) THEN
        CALL DGEMM( 'N', 'N', N, K, K, ONE, U2, LDU2, Q, LDQ, ZERO, U,  &
     &LDU )
        GO TO 100
      END IF
      IF( CTOT( 1 ).GT.0 ) THEN
        CALL DGEMM( 'N', 'N', NL, K, CTOT( 1 ), ONE, U2( 1, 2 ), LDU2,  &
     &Q( 2, 1 ), LDQ, ZERO, U( 1, 1 ), LDU )
        IF( CTOT( 3 ).GT.0 ) THEN
          KTEMP = 2 + CTOT( 1 ) + CTOT( 2 )
          CALL DGEMM( 'N', 'N', NL, K, CTOT( 3 ), ONE, U2( 1, KTEMP ),  &
     &LDU2, Q( KTEMP, 1 ), LDQ, ONE, U( 1, 1 ), LDU )
        END IF
      ELSE IF( CTOT( 3 ).GT.0 ) THEN
        KTEMP = 2 + CTOT( 1 ) + CTOT( 2 )
        CALL DGEMM( 'N', 'N', NL, K, CTOT( 3 ), ONE, U2( 1, KTEMP ),    &
     &LDU2, Q( KTEMP, 1 ), LDQ, ZERO, U( 1, 1 ), LDU )
      ELSE
        CALL DLACPY( 'F', NL, K, U2, LDU2, U, LDU )
      END IF
      CALL DCOPL( K, Q( 1, 1 ), LDQ, U( NLP1, 1 ), LDU )
      KTEMP = 2 + CTOT( 1 )
      CTEMP = CTOT( 2 ) + CTOT( 3 )
      CALL DGEMM( 'N', 'N', NR, K, CTEMP, ONE, U2( NLP2, KTEMP ), LDU2, &
     &Q( KTEMP, 1 ), LDQ, ZERO, U( NLP2, 1 ), LDU )
!
!     Generate the right singular vectors.
!
  100 CONTINUE
      DO I = 1, K
        TEMP = DNRM2( K, VT( 1, I ), 1 )
        Q( I, 1 ) = VT( 1, I ) / TEMP
        DO J = 2, K
          JC = IDXC( J )
          Q( I, J ) = VT( JC, I ) / TEMP
        enddo
      enddo
!
!     Update the right singular vector matrix.
!
      IF( K.EQ.2 ) THEN
        CALL DGEMM( 'N', 'N', K, M, K, ONE, Q, LDQ, VT2, LDVT2, ZERO,   &
     &VT, LDVT )
        RETURN
      END IF
      KTEMP = 1 + CTOT( 1 )
      CALL DGEMM( 'N', 'N', K, NLP1, KTEMP, ONE, Q( 1, 1 ), LDQ,        &
     &VT2( 1, 1 ), LDVT2, ZERO, VT( 1, 1 ), LDVT )
      KTEMP = 2 + CTOT( 1 ) + CTOT( 2 )
      IF( KTEMP.LE.LDVT2 )                                              &
     &CALL DGEMM( 'N', 'N', K, NLP1, CTOT( 3 ), ONE, Q( 1, KTEMP ),     &
     &LDQ, VT2( KTEMP, 1 ), LDVT2, ONE, VT( 1, 1 ),                     &
     &LDVT )
!
      KTEMP = CTOT( 1 ) + 1
      NRP1 = NR + SQRE
      IF( KTEMP.GT.1 ) THEN
        DO I = 1, K
          Q( I, KTEMP ) = Q( I, 1 )
        enddo
        DO I = NLP2, M
          VT2( KTEMP, I ) = VT2( 1, I )
        enddo
      END IF
      CTEMP = 1 + CTOT( 2 ) + CTOT( 3 )
      CALL DGEMM( 'N', 'N', K, NRP1, CTEMP, ONE, Q( 1, KTEMP ), LDQ,    &
     &VT2( KTEMP, NLP2 ), LDVT2, ZERO, VT( 1, NLP2 ), LDVT )
!
      RETURN
!
!     End of DLASD3
!
      END
      SUBROUTINE DLASD4( N, I, D, Z, DELTA, RHO, SIGMA, WORK, INFO )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            I, INFO, N
      DOUBLE PRECISION   RHO, SIGMA
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), DELTA( * ), WORK( * ), Z( * )
!     ..
!
!  Purpose
!  =======
!
!  This subroutine computes the square root of the I-th updated
!  eigenvalue of a positive symmetric rank-one modification to
!  a positive diagonal matrix whose entries are given as the squares
!  of the corresponding entries in the array d, and that
!
!         0 <= D(i) < D(j)  for  i < j
!
!  and that RHO > 0. This is arranged by the calling routine, and is
!  no loss in generality.  The rank-one modified system is thus
!
!         diag( D ) * diag( D ) +  RHO *  Z * Z_transpose.
!
!  where we assume the Euclidean norm of Z is 1.
!
!  The method consists of approximating the rational functions in the
!  secular equation by simpler interpolating rational functions.
!
!  Arguments
!  =========
!
!  N      (input) INTEGER
!         The length of all arrays.
!
!  I      (input) INTEGER
!         The index of the eigenvalue to be computed.  1 <= I <= N.
!
!  D      (input) DOUBLE PRECISION array, dimension ( N )
!         The original eigenvalues.  It is assumed that they are in
!         order, 0 <= D(I) < D(J)  for I < J.
!
!  Z      (input) DOUBLE PRECISION array, dimension ( N )
!         The components of the updating vector.
!
!  DELTA  (output) DOUBLE PRECISION array, dimension ( N )
!         If N .ne. 1, DELTA contains (D(j) - sigma_I) in its  j-th
!         component.  If N = 1, then DELTA(1) = 1.  The vector DELTA
!         contains the information necessary to construct the
!         (singular) eigenvectors.
!
!  RHO    (input) DOUBLE PRECISION
!         The scalar in the symmetric updating formula.
!
!  SIGMA  (output) DOUBLE PRECISION
!         The computed sigma_I, the I-th updated eigenvalue.
!
!  WORK   (workspace) DOUBLE PRECISION array, dimension ( N )
!         If N .ne. 1, WORK contains (D(j) + sigma_I) in its  j-th
!         component.  If N = 1, then WORK( 1 ) = 1.
!
!  INFO   (output) INTEGER
!         = 0:  successful exit
!         > 0:  if INFO = 1, the updating process failed.
!
!  Internal Parameters
!  ===================
!
!  Logical variable ORGATI (origin-at-i?) is used for distinguishing
!  whether D(i) or D(i+1) is treated as the origin.
!
!            ORGATI = .true.    origin at i
!            ORGATI = .false.   origin at i+1
!
!  Logical variable SWTCH3 (switch-for-3-poles?) is for noting
!  if we are working with THREE poles!
!
!  MAXIT is the maximum number of iterations allowed for each
!  eigenvalue.
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ren-Cang Li, Computer Science Division, University of California
!     at Berkeley, USA
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 20 )
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE, FOUR, EIGHT, TEN
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0,   &
     &THREE = 3.0D+0, FOUR = 4.0D+0, EIGHT = 8.0D+0,                    &
     &TEN = 10.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            ORGATI, SWTCH, SWTCH3
      INTEGER            II, IIM1, IIP1, IP1, ITER, J, NITER
      DOUBLE PRECISION   A, B, C, DELSQ, DELSQ2, DPHI, DPSI, DTIIM,     &
     &DTIIP, DTIPSQ, DTISQ, DTNSQ, DTNSQ1, DW, EPS,                     &
     &ERRETM, ETA, PHI, PREW, PSI, RHOINV, SG2LB,                       &
     &SG2UB, TAU, TEMP, TEMP1, TEMP2, W
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   DD( 3 ), ZZ( 3 )
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAED6, DLASD5
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Since this routine is called in an inner loop, we do no argument
!     checking.
!
!     Quick return for N=1 and 2.
!
      INFO = 0
      IF( N.EQ.1 ) THEN
!
!        Presumably, I=1 upon entry
!
        SIGMA = SQRT( D( 1 )*D( 1 )+RHO*Z( 1 )*Z( 1 ) )
        DELTA( 1 ) = ONE
        WORK( 1 ) = ONE
        RETURN
      END IF
      IF( N.EQ.2 ) THEN
        CALL DLASD5( I, D, Z, DELTA, RHO, SIGMA, WORK )
        RETURN
      END IF
!
!     Compute machine epsilon
!
      EPS = DLAMCH( 'Epsilon' )
      RHOINV = ONE / RHO
!
!     The case I = N
!
      IF( I.EQ.N ) THEN
!
!        Initialize some basic variables
!
        II = N - 1
        NITER = 1
!
!        Calculate initial guess
!
        TEMP = RHO / TWO
!
!        If ||Z||_2 is not one, then TEMP should be set to
!        RHO * ||Z||_2^2 / TWO
!
        TEMP1 = TEMP / ( D( N )+SQRT( D( N )*D( N )+TEMP ) )
        DO J = 1, N
          WORK( J ) = D( J ) + D( N ) + TEMP1
          DELTA( J ) = ( D( J )-D( N ) ) - TEMP1
        enddo
!
        PSI = ZERO
        DO J = 1, N - 2
          PSI = PSI + Z( J )*Z( J ) / ( DELTA( J )*WORK( J ) )
        enddo
!
        C = RHOINV + PSI
        W = C + Z( II )*Z( II ) / ( DELTA( II )*WORK( II ) ) +          &
     &Z( N )*Z( N ) / ( DELTA( N )*WORK( N ) )
!
        IF( W.LE.ZERO ) THEN
          TEMP1 = SQRT( D( N )*D( N )+RHO )
          TEMP = Z( N-1 )*Z( N-1 ) / ( ( D( N-1 )+TEMP1 )*              &
     &( D( N )-D( N-1 )+RHO / ( D( N )+TEMP1 ) ) ) +                    &
     &Z( N )*Z( N ) / RHO
!
!           The following TAU is to approximate
!           SIGMA_n^2 - D( N )*D( N )
!
          IF( C.LE.TEMP ) THEN
            TAU = RHO
          ELSE
            DELSQ = ( D( N )-D( N-1 ) )*( D( N )+D( N-1 ) )
            A = -C*DELSQ + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
            B = Z( N )*Z( N )*DELSQ
            IF( A.LT.ZERO ) THEN
              TAU = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
            ELSE
              TAU = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
            END IF
          END IF
!
!           It can be proved that
!               D(N)^2+RHO/2 <= SIGMA_n^2 < D(N)^2+TAU <= D(N)^2+RHO
!
        ELSE
          DELSQ = ( D( N )-D( N-1 ) )*( D( N )+D( N-1 ) )
          A = -C*DELSQ + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
          B = Z( N )*Z( N )*DELSQ
!
!           The following TAU is to approximate
!           SIGMA_n^2 - D( N )*D( N )
!
          IF( A.LT.ZERO ) THEN
            TAU = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
          ELSE
            TAU = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
          END IF
!
!           It can be proved that
!           D(N)^2 < D(N)^2+TAU < SIGMA(N)^2 < D(N)^2+RHO/2
!
        END IF
!
!        The following ETA is to approximate SIGMA_n - D( N )
!
        ETA = TAU / ( D( N )+SQRT( D( N )*D( N )+TAU ) )
!
        SIGMA = D( N ) + ETA
        DO J = 1, N
          DELTA( J ) = ( D( J )-D( I ) ) - ETA
          WORK( J ) = D( J ) + D( I ) + ETA
        enddo
!
!        Evaluate PSI and the derivative DPSI
!
        DPSI = ZERO
        PSI = ZERO
        ERRETM = ZERO
        DO J = 1, II
          TEMP = Z( J ) / ( DELTA( J )*WORK( J ) )
          PSI = PSI + Z( J )*TEMP
          DPSI = DPSI + TEMP*TEMP
          ERRETM = ERRETM + PSI
        enddo
        ERRETM = ABS( ERRETM )
!
!        Evaluate PHI and the derivative DPHI
!
        TEMP = Z( N ) / ( DELTA( N )*WORK( N ) )
        PHI = Z( N )*TEMP
        DPHI = TEMP*TEMP
        ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV +           &
     &ABS( TAU )*( DPSI+DPHI )
!
        W = RHOINV + PHI + PSI
!
!        Test for convergence
!
        IF( ABS( W ).LE.EPS*ERRETM ) THEN
          GO TO 240
        END IF
!
!        Calculate the new step
!
        NITER = NITER + 1
        DTNSQ1 = WORK( N-1 )*DELTA( N-1 )
        DTNSQ = WORK( N )*DELTA( N )
        C = W - DTNSQ1*DPSI - DTNSQ*DPHI
        A = ( DTNSQ+DTNSQ1 )*W - DTNSQ*DTNSQ1*( DPSI+DPHI )
        B = DTNSQ*DTNSQ1*W
        IF( C.LT.ZERO )                                                 &
     &C = ABS( C )
        IF( C.EQ.ZERO ) THEN
          ETA = RHO - SIGMA*SIGMA
        ELSE IF( A.GE.ZERO ) THEN
          ETA = ( A+SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
        ELSE
          ETA = TWO*B / ( A-SQRT( ABS( A*A-FOUR*B*C ) ) )
        END IF
!
!        Note, eta should be positive if w is negative, and
!        eta should be negative otherwise. However,
!        if for some reason caused by roundoff, eta*w > 0,
!        we simply use one Newton step instead. This way
!        will guarantee eta*w < 0.
!
        IF( W*ETA.GT.ZERO )                                             &
     &ETA = -W / ( DPSI+DPHI )
        TEMP = ETA - DTNSQ
        IF( TEMP.GT.RHO )                                               &
     &ETA = RHO + DTNSQ
!
        TAU = TAU + ETA
        ETA = ETA / ( SIGMA+SQRT( ETA+SIGMA*SIGMA ) )
        DO J = 1, N
          DELTA( J ) = DELTA( J ) - ETA
          WORK( J ) = WORK( J ) + ETA
        enddo
!
        SIGMA = SIGMA + ETA
!
!        Evaluate PSI and the derivative DPSI
!
        DPSI = ZERO
        PSI = ZERO
        ERRETM = ZERO
        DO J = 1, II
          TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
          PSI = PSI + Z( J )*TEMP
          DPSI = DPSI + TEMP*TEMP
          ERRETM = ERRETM + PSI
        enddo
        ERRETM = ABS( ERRETM )
!
!        Evaluate PHI and the derivative DPHI
!
        TEMP = Z( N ) / ( WORK( N )*DELTA( N ) )
        PHI = Z( N )*TEMP
        DPHI = TEMP*TEMP
        ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV +           &
     &ABS( TAU )*( DPSI+DPHI )
!
        W = RHOINV + PHI + PSI
!
!        Main loop to update the values of the array   DELTA
!
        ITER = NITER + 1
!
        DO NITER = ITER, MAXIT
!
!           Test for convergence
!
          IF( ABS( W ).LE.EPS*ERRETM ) THEN
            GO TO 240
          END IF
!
!           Calculate the new step
!
          DTNSQ1 = WORK( N-1 )*DELTA( N-1 )
          DTNSQ = WORK( N )*DELTA( N )
          C = W - DTNSQ1*DPSI - DTNSQ*DPHI
          A = ( DTNSQ+DTNSQ1 )*W - DTNSQ1*DTNSQ*( DPSI+DPHI )
          B = DTNSQ1*DTNSQ*W
          IF( A.GE.ZERO ) THEN
            ETA = ( A+SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
          ELSE
            ETA = TWO*B / ( A-SQRT( ABS( A*A-FOUR*B*C ) ) )
          END IF
!
!           Note, eta should be positive if w is negative, and
!           eta should be negative otherwise. However,
!           if for some reason caused by roundoff, eta*w > 0,
!           we simply use one Newton step instead. This way
!           will guarantee eta*w < 0.
!
          IF( W*ETA.GT.ZERO )                                           &
     &ETA = -W / ( DPSI+DPHI )
          TEMP = ETA - DTNSQ
          IF( TEMP.LE.ZERO )                                            &
     &ETA = ETA / TWO
!
          TAU = TAU + ETA
          ETA = ETA / ( SIGMA+SQRT( ETA+SIGMA*SIGMA ) )
          DO J = 1, N
            DELTA( J ) = DELTA( J ) - ETA
            WORK( J ) = WORK( J ) + ETA
          enddo
!
          SIGMA = SIGMA + ETA
!
!           Evaluate PSI and the derivative DPSI
!
          DPSI = ZERO
          PSI = ZERO
          ERRETM = ZERO
          DO J = 1, II
            TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
          enddo
          ERRETM = ABS( ERRETM )
!
!           Evaluate PHI and the derivative DPHI
!
          TEMP = Z( N ) / ( WORK( N )*DELTA( N ) )
          PHI = Z( N )*TEMP
          DPHI = TEMP*TEMP
          ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV +         &
     &ABS( TAU )*( DPSI+DPHI )
!
          W = RHOINV + PHI + PSI
        enddo
!
!        Return with INFO = 1, NITER = MAXIT and not converged
!
        INFO = 1
        GO TO 240
!
!        End for the case I = N
!
      ELSE
!
!        The case for I < N
!
        NITER = 1
        IP1 = I + 1
!
!        Calculate initial guess
!
        DELSQ = ( D( IP1 )-D( I ) )*( D( IP1 )+D( I ) )
        DELSQ2 = DELSQ / TWO
        TEMP = DELSQ2 / ( D( I )+SQRT( D( I )*D( I )+DELSQ2 ) )
        DO J = 1, N
          WORK( J ) = D( J ) + D( I ) + TEMP
          DELTA( J ) = ( D( J )-D( I ) ) - TEMP
        enddo
!
        PSI = ZERO
        DO J = 1, I - 1
          PSI = PSI + Z( J )*Z( J ) / ( WORK( J )*DELTA( J ) )
        enddo
!
        PHI = ZERO
        DO J = N, I + 2, -1
          PHI = PHI + Z( J )*Z( J ) / ( WORK( J )*DELTA( J ) )
        enddo
        C = RHOINV + PSI + PHI
        W = C + Z( I )*Z( I ) / ( WORK( I )*DELTA( I ) ) +              &
     &Z( IP1 )*Z( IP1 ) / ( WORK( IP1 )*DELTA( IP1 ) )
!
        IF( W.GT.ZERO ) THEN
!
!           d(i)^2 < the ith sigma^2 < (d(i)^2+d(i+1)^2)/2
!
!           We choose d(i) as origin.
!
          ORGATI = .TRUE.
          SG2LB = ZERO
          SG2UB = DELSQ2
          A = C*DELSQ + Z( I )*Z( I ) + Z( IP1 )*Z( IP1 )
          B = Z( I )*Z( I )*DELSQ
          IF( A.GT.ZERO ) THEN
            TAU = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
          ELSE
            TAU = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
          END IF
!
!           TAU now is an estimation of SIGMA^2 - D( I )^2. The
!           following, however, is the corresponding estimation of
!           SIGMA - D( I ).
!
          ETA = TAU / ( D( I )+SQRT( D( I )*D( I )+TAU ) )
        ELSE
!
!           (d(i)^2+d(i+1)^2)/2 <= the ith sigma^2 < d(i+1)^2/2
!
!           We choose d(i+1) as origin.
!
          ORGATI = .FALSE.
          SG2LB = -DELSQ2
          SG2UB = ZERO
          A = C*DELSQ - Z( I )*Z( I ) - Z( IP1 )*Z( IP1 )
          B = Z( IP1 )*Z( IP1 )*DELSQ
          IF( A.LT.ZERO ) THEN
            TAU = TWO*B / ( A-SQRT( ABS( A*A+FOUR*B*C ) ) )
          ELSE
            TAU = -( A+SQRT( ABS( A*A+FOUR*B*C ) ) ) / ( TWO*C )
          END IF
!
!           TAU now is an estimation of SIGMA^2 - D( IP1 )^2. The
!           following, however, is the corresponding estimation of
!           SIGMA - D( IP1 ).
!
          ETA = TAU / ( D( IP1 )+SQRT( ABS( D( IP1 )*D( IP1 )+          &
     &TAU ) ) )
        END IF
!
        IF( ORGATI ) THEN
          II = I
          SIGMA = D( I ) + ETA
          DO J = 1, N
            WORK( J ) = D( J ) + D( I ) + ETA
            DELTA( J ) = ( D( J )-D( I ) ) - ETA
          enddo
        ELSE
          II = I + 1
          SIGMA = D( IP1 ) + ETA
          DO J = 1, N
            WORK( J ) = D( J ) + D( IP1 ) + ETA
            DELTA( J ) = ( D( J )-D( IP1 ) ) - ETA
          enddo
        END IF
        IIM1 = II - 1
        IIP1 = II + 1
!
!        Evaluate PSI and the derivative DPSI
!
        DPSI = ZERO
        PSI = ZERO
        ERRETM = ZERO
        DO J = 1, IIM1
          TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
          PSI = PSI + Z( J )*TEMP
          DPSI = DPSI + TEMP*TEMP
          ERRETM = ERRETM + PSI
        enddo
        ERRETM = ABS( ERRETM )
!
!        Evaluate PHI and the derivative DPHI
!
        DPHI = ZERO
        PHI = ZERO
        DO J = N, IIP1, -1
          TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
          PHI = PHI + Z( J )*TEMP
          DPHI = DPHI + TEMP*TEMP
          ERRETM = ERRETM + PHI
        enddo
!
        W = RHOINV + PHI + PSI
!
!        W is the value of the secular function with
!        its ii-th element removed.
!
        SWTCH3 = .FALSE.
        IF( ORGATI ) THEN
          IF( W.LT.ZERO )                                               &
     &SWTCH3 = .TRUE.
        ELSE
          IF( W.GT.ZERO )                                               &
     &SWTCH3 = .TRUE.
        END IF
        IF( II.EQ.1 .OR. II.EQ.N )                                      &
     &SWTCH3 = .FALSE.
!
        TEMP = Z( II ) / ( WORK( II )*DELTA( II ) )
        DW = DPSI + DPHI + TEMP*TEMP
        TEMP = Z( II )*TEMP
        W = W + TEMP
        ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV +              &
     &THREE*ABS( TEMP ) + ABS( TAU )*DW
!
!        Test for convergence
!
        IF( ABS( W ).LE.EPS*ERRETM ) THEN
          GO TO 240
        END IF
!
        IF( W.LE.ZERO ) THEN
          SG2LB = MAX( SG2LB, TAU )
        ELSE
          SG2UB = MIN( SG2UB, TAU )
        END IF
!
!        Calculate the new step
!
        NITER = NITER + 1
        IF( .NOT.SWTCH3 ) THEN
          DTIPSQ = WORK( IP1 )*DELTA( IP1 )
          DTISQ = WORK( I )*DELTA( I )
          IF( ORGATI ) THEN
            C = W - DTIPSQ*DW + DELSQ*( Z( I ) / DTISQ )**2
          ELSE
            C = W - DTISQ*DW - DELSQ*( Z( IP1 ) / DTIPSQ )**2
          END IF
          A = ( DTIPSQ+DTISQ )*W - DTIPSQ*DTISQ*DW
          B = DTIPSQ*DTISQ*W
          IF( C.EQ.ZERO ) THEN
            IF( A.EQ.ZERO ) THEN
              IF( ORGATI ) THEN
                A = Z( I )*Z( I ) + DTIPSQ*DTIPSQ*( DPSI+DPHI )
              ELSE
                A = Z( IP1 )*Z( IP1 ) + DTISQ*DTISQ*( DPSI+DPHI )
              END IF
            END IF
            ETA = B / A
          ELSE IF( A.LE.ZERO ) THEN
            ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
          ELSE
            ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
          END IF
        ELSE
!
!           Interpolation using THREE most relevant poles
!
          DTIIM = WORK( IIM1 )*DELTA( IIM1 )
          DTIIP = WORK( IIP1 )*DELTA( IIP1 )
          TEMP = RHOINV + PSI + PHI
          IF( ORGATI ) THEN
            TEMP1 = Z( IIM1 ) / DTIIM
            TEMP1 = TEMP1*TEMP1
            C = ( TEMP - DTIIP*( DPSI+DPHI ) ) -                        &
     &( D( IIM1 )-D( IIP1 ) )*( D( IIM1 )+D( IIP1 ) )*TEMP1
            ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
            IF( DPSI.LT.TEMP1 ) THEN
              ZZ( 3 ) = DTIIP*DTIIP*DPHI
            ELSE
              ZZ( 3 ) = DTIIP*DTIIP*( ( DPSI-TEMP1 )+DPHI )
            END IF
          ELSE
            TEMP1 = Z( IIP1 ) / DTIIP
            TEMP1 = TEMP1*TEMP1
            C = ( TEMP - DTIIM*( DPSI+DPHI ) ) -                        &
     &( D( IIP1 )-D( IIM1 ) )*( D( IIM1 )+D( IIP1 ) )*TEMP1
            IF( DPHI.LT.TEMP1 ) THEN
              ZZ( 1 ) = DTIIM*DTIIM*DPSI
            ELSE
              ZZ( 1 ) = DTIIM*DTIIM*( DPSI+( DPHI-TEMP1 ) )
            END IF
            ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
          END IF
          ZZ( 2 ) = Z( II )*Z( II )
          DD( 1 ) = DTIIM
          DD( 2 ) = DELTA( II )*WORK( II )
          DD( 3 ) = DTIIP
          CALL DLAED6( NITER, ORGATI, C, DD, ZZ, W, ETA, INFO )
          IF( INFO.NE.0 )                                               &
     &GO TO 240
        END IF
!
!        Note, eta should be positive if w is negative, and
!        eta should be negative otherwise. However,
!        if for some reason caused by roundoff, eta*w > 0,
!        we simply use one Newton step instead. This way
!        will guarantee eta*w < 0.
!
        IF( W*ETA.GE.ZERO )                                             &
     &ETA = -W / DW
        IF( ORGATI ) THEN
          TEMP1 = WORK( I )*DELTA( I )
          TEMP = ETA - TEMP1
        ELSE
          TEMP1 = WORK( IP1 )*DELTA( IP1 )
          TEMP = ETA - TEMP1
        END IF
        IF( TEMP.GT.SG2UB .OR. TEMP.LT.SG2LB ) THEN
          IF( W.LT.ZERO ) THEN
            ETA = ( SG2UB-TAU ) / TWO
          ELSE
            ETA = ( SG2LB-TAU ) / TWO
          END IF
        END IF
!
        TAU = TAU + ETA
        ETA = ETA / ( SIGMA+SQRT( SIGMA*SIGMA+ETA ) )
!
        PREW = W
!
        SIGMA = SIGMA + ETA
        DO J = 1, N
          WORK( J ) = WORK( J ) + ETA
          DELTA( J ) = DELTA( J ) - ETA
        enddo
!
!        Evaluate PSI and the derivative DPSI
!
        DPSI = ZERO
        PSI = ZERO
        ERRETM = ZERO
        DO J = 1, IIM1
          TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
          PSI = PSI + Z( J )*TEMP
          DPSI = DPSI + TEMP*TEMP
          ERRETM = ERRETM + PSI
        enddo
        ERRETM = ABS( ERRETM )
!
!        Evaluate PHI and the derivative DPHI
!
        DPHI = ZERO
        PHI = ZERO
        DO J = N, IIP1, -1
          TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
          PHI = PHI + Z( J )*TEMP
          DPHI = DPHI + TEMP*TEMP
          ERRETM = ERRETM + PHI
        enddo
!
        TEMP = Z( II ) / ( WORK( II )*DELTA( II ) )
        DW = DPSI + DPHI + TEMP*TEMP
        TEMP = Z( II )*TEMP
        W = RHOINV + PHI + PSI + TEMP
        ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV +              &
     &THREE*ABS( TEMP ) + ABS( TAU )*DW
!
        IF( W.LE.ZERO ) THEN
          SG2LB = MAX( SG2LB, TAU )
        ELSE
          SG2UB = MIN( SG2UB, TAU )
        END IF
!
        SWTCH = .FALSE.
        IF( ORGATI ) THEN
          IF( -W.GT.ABS( PREW ) / TEN )                                 &
     &SWTCH = .TRUE.
        ELSE
          IF( W.GT.ABS( PREW ) / TEN )                                  &
     &SWTCH = .TRUE.
        END IF
!
!        Main loop to update the values of the array   DELTA and WORK
!
        ITER = NITER + 1
!
        DO NITER = ITER, MAXIT
!
!           Test for convergence
!
          IF( ABS( W ).LE.EPS*ERRETM ) THEN
            GO TO 240
          END IF
!
!           Calculate the new step
!
          IF( .NOT.SWTCH3 ) THEN
            DTIPSQ = WORK( IP1 )*DELTA( IP1 )
            DTISQ = WORK( I )*DELTA( I )
            IF( .NOT.SWTCH ) THEN
              IF( ORGATI ) THEN
                C = W - DTIPSQ*DW + DELSQ*( Z( I ) / DTISQ )**2
              ELSE
                C = W - DTISQ*DW - DELSQ*( Z( IP1 ) / DTIPSQ )**2
              END IF
            ELSE
              TEMP = Z( II ) / ( WORK( II )*DELTA( II ) )
              IF( ORGATI ) THEN
                DPSI = DPSI + TEMP*TEMP
              ELSE
                DPHI = DPHI + TEMP*TEMP
              END IF
              C = W - DTISQ*DPSI - DTIPSQ*DPHI
            END IF
            A = ( DTIPSQ+DTISQ )*W - DTIPSQ*DTISQ*DW
            B = DTIPSQ*DTISQ*W
            IF( C.EQ.ZERO ) THEN
              IF( A.EQ.ZERO ) THEN
                IF( .NOT.SWTCH ) THEN
                  IF( ORGATI ) THEN
                    A = Z( I )*Z( I ) + DTIPSQ*DTIPSQ*                  &
     &( DPSI+DPHI )
                  ELSE
                    A = Z( IP1 )*Z( IP1 ) +                             &
     &DTISQ*DTISQ*( DPSI+DPHI )
                  END IF
                ELSE
                  A = DTISQ*DTISQ*DPSI + DTIPSQ*DTIPSQ*DPHI
                END IF
              END IF
              ETA = B / A
            ELSE IF( A.LE.ZERO ) THEN
              ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            ELSE
              ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
            END IF
          ELSE
!
!              Interpolation using THREE most relevant poles
!
            DTIIM = WORK( IIM1 )*DELTA( IIM1 )
            DTIIP = WORK( IIP1 )*DELTA( IIP1 )
            TEMP = RHOINV + PSI + PHI
            IF( SWTCH ) THEN
              C = TEMP - DTIIM*DPSI - DTIIP*DPHI
              ZZ( 1 ) = DTIIM*DTIIM*DPSI
              ZZ( 3 ) = DTIIP*DTIIP*DPHI
            ELSE
              IF( ORGATI ) THEN
                TEMP1 = Z( IIM1 ) / DTIIM
                TEMP1 = TEMP1*TEMP1
                TEMP2 = ( D( IIM1 )-D( IIP1 ) )*                        &
     &( D( IIM1 )+D( IIP1 ) )*TEMP1
                C = TEMP - DTIIP*( DPSI+DPHI ) - TEMP2
                ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
                IF( DPSI.LT.TEMP1 ) THEN
                  ZZ( 3 ) = DTIIP*DTIIP*DPHI
                ELSE
                  ZZ( 3 ) = DTIIP*DTIIP*( ( DPSI-TEMP1 )+DPHI )
                END IF
              ELSE
                TEMP1 = Z( IIP1 ) / DTIIP
                TEMP1 = TEMP1*TEMP1
                TEMP2 = ( D( IIP1 )-D( IIM1 ) )*                        &
     &( D( IIM1 )+D( IIP1 ) )*TEMP1
                C = TEMP - DTIIM*( DPSI+DPHI ) - TEMP2
                IF( DPHI.LT.TEMP1 ) THEN
                  ZZ( 1 ) = DTIIM*DTIIM*DPSI
                ELSE
                  ZZ( 1 ) = DTIIM*DTIIM*( DPSI+( DPHI-TEMP1 ) )
                END IF
                ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
              END IF
            END IF
            DD( 1 ) = DTIIM
            DD( 2 ) = DELTA( II )*WORK( II )
            DD( 3 ) = DTIIP
            CALL DLAED6( NITER, ORGATI, C, DD, ZZ, W, ETA, INFO )
            IF( INFO.NE.0 )                                             &
     &GO TO 240
          END IF
!
!           Note, eta should be positive if w is negative, and
!           eta should be negative otherwise. However,
!           if for some reason caused by roundoff, eta*w > 0,
!           we simply use one Newton step instead. This way
!           will guarantee eta*w < 0.
!
          IF( W*ETA.GE.ZERO )                                           &
     &ETA = -W / DW
          IF( ORGATI ) THEN
            TEMP1 = WORK( I )*DELTA( I )
            TEMP = ETA - TEMP1
          ELSE
            TEMP1 = WORK( IP1 )*DELTA( IP1 )
            TEMP = ETA - TEMP1
          END IF
          IF( TEMP.GT.SG2UB .OR. TEMP.LT.SG2LB ) THEN
            IF( W.LT.ZERO ) THEN
              ETA = ( SG2UB-TAU ) / TWO
            ELSE
              ETA = ( SG2LB-TAU ) / TWO
            END IF
          END IF
!
          TAU = TAU + ETA
          ETA = ETA / ( SIGMA+SQRT( SIGMA*SIGMA+ETA ) )
!
          SIGMA = SIGMA + ETA
          DO J = 1, N
            WORK( J ) = WORK( J ) + ETA
            DELTA( J ) = DELTA( J ) - ETA
          enddo
!
          PREW = W
!
!           Evaluate PSI and the derivative DPSI
!
          DPSI = ZERO
          PSI = ZERO
          ERRETM = ZERO
          DO J = 1, IIM1
            TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
          enddo
          ERRETM = ABS( ERRETM )
!
!           Evaluate PHI and the derivative DPHI
!
          DPHI = ZERO
          PHI = ZERO
          DO J = N, IIP1, -1
            TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
            PHI = PHI + Z( J )*TEMP
            DPHI = DPHI + TEMP*TEMP
            ERRETM = ERRETM + PHI
          enddo
!
          TEMP = Z( II ) / ( WORK( II )*DELTA( II ) )
          DW = DPSI + DPHI + TEMP*TEMP
          TEMP = Z( II )*TEMP
          W = RHOINV + PHI + PSI + TEMP
          ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV +            &
     &THREE*ABS( TEMP ) + ABS( TAU )*DW
          IF( W*PREW.GT.ZERO .AND. ABS( W ).GT.ABS( PREW ) / TEN )      &
     &SWTCH = .NOT.SWTCH
!
          IF( W.LE.ZERO ) THEN
            SG2LB = MAX( SG2LB, TAU )
          ELSE
            SG2UB = MIN( SG2UB, TAU )
          END IF
!
        enddo
!
!        Return with INFO = 1, NITER = MAXIT and not converged
!
        INFO = 1
!
      END IF
!
  240 CONTINUE
      RETURN
!
!     End of DLASD4
!
      END
      SUBROUTINE DLASD5( I, D, Z, DELTA, RHO, DSIGMA, WORK )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            I
      DOUBLE PRECISION   DSIGMA, RHO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( 2 ), DELTA( 2 ), WORK( 2 ), Z( 2 )
!     ..
!
!  Purpose
!  =======
!
!  This subroutine computes the square root of the I-th eigenvalue
!  of a positive symmetric rank-one modification of a 2-by-2 diagonal
!  matrix
!
!             diag( D ) * diag( D ) +  RHO *  Z * transpose(Z) .
!
!  The diagonal entries in the array D are assumed to satisfy
!
!             0 <= D(i) < D(j)  for  i < j .
!
!  We also assume RHO > 0 and that the Euclidean norm of the vector
!  Z is one.
!
!  Arguments
!  =========
!
!  I      (input) INTEGER
!         The index of the eigenvalue to be computed.  I = 1 or I = 2.
!
!  D      (input) DOUBLE PRECISION array, dimension ( 2 )
!         The original eigenvalues.  We assume 0 <= D(1) < D(2).
!
!  Z      (input) DOUBLE PRECISION array, dimension ( 2 )
!         The components of the updating vector.
!
!  DELTA  (output) DOUBLE PRECISION array, dimension ( 2 )
!         Contains (D(j) - sigma_I) in its  j-th component.
!         The vector DELTA contains the information necessary
!         to construct the eigenvectors.
!
!  RHO    (input) DOUBLE PRECISION
!         The scalar in the symmetric updating formula.
!
!  DSIGMA (output) DOUBLE PRECISION
!         The computed sigma_I, the I-th updated eigenvalue.
!
!  WORK   (workspace) DOUBLE PRECISION array, dimension ( 2 )
!         WORK contains (D(j) + sigma_I) in its  j-th component.
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ren-Cang Li, Computer Science Division, University of California
!     at Berkeley, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE, FOUR
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0,   &
     &THREE = 3.0D+0, FOUR = 4.0D+0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   B, C, DEL, DELSQ, TAU, W
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
!     ..
!     .. Executable Statements ..
!
      DEL = D( 2 ) - D( 1 )
      DELSQ = DEL*( D( 2 )+D( 1 ) )
      IF( I.EQ.1 ) THEN
        W = ONE + FOUR*RHO*( Z( 2 )*Z( 2 ) / ( D( 1 )+THREE*D( 2 ) )-   &
     &Z( 1 )*Z( 1 ) / ( THREE*D( 1 )+D( 2 ) ) ) / DEL
        IF( W.GT.ZERO ) THEN
          B = DELSQ + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
          C = RHO*Z( 1 )*Z( 1 )*DELSQ
!
!           B > ZERO, always
!
!           The following TAU is DSIGMA * DSIGMA - D( 1 ) * D( 1 )
!
          TAU = TWO*C / ( B+SQRT( ABS( B*B-FOUR*C ) ) )
!
!           The following TAU is DSIGMA - D( 1 )
!
          TAU = TAU / ( D( 1 )+SQRT( D( 1 )*D( 1 )+TAU ) )
          DSIGMA = D( 1 ) + TAU
          DELTA( 1 ) = -TAU
          DELTA( 2 ) = DEL - TAU
          WORK( 1 ) = TWO*D( 1 ) + TAU
          WORK( 2 ) = ( D( 1 )+TAU ) + D( 2 )
!           DELTA( 1 ) = -Z( 1 ) / TAU
!           DELTA( 2 ) = Z( 2 ) / ( DEL-TAU )
        ELSE
          B = -DELSQ + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
          C = RHO*Z( 2 )*Z( 2 )*DELSQ
!
!           The following TAU is DSIGMA * DSIGMA - D( 2 ) * D( 2 )
!
          IF( B.GT.ZERO ) THEN
            TAU = -TWO*C / ( B+SQRT( B*B+FOUR*C ) )
          ELSE
            TAU = ( B-SQRT( B*B+FOUR*C ) ) / TWO
          END IF
!
!           The following TAU is DSIGMA - D( 2 )
!
          TAU = TAU / ( D( 2 )+SQRT( ABS( D( 2 )*D( 2 )+TAU ) ) )
          DSIGMA = D( 2 ) + TAU
          DELTA( 1 ) = -( DEL+TAU )
          DELTA( 2 ) = -TAU
          WORK( 1 ) = D( 1 ) + TAU + D( 2 )
          WORK( 2 ) = TWO*D( 2 ) + TAU
!           DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
!           DELTA( 2 ) = -Z( 2 ) / TAU
        END IF
!        TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
!        DELTA( 1 ) = DELTA( 1 ) / TEMP
!        DELTA( 2 ) = DELTA( 2 ) / TEMP
      ELSE
!
!        Now I=2
!
        B = -DELSQ + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
        C = RHO*Z( 2 )*Z( 2 )*DELSQ
!
!        The following TAU is DSIGMA * DSIGMA - D( 2 ) * D( 2 )
!
        IF( B.GT.ZERO ) THEN
          TAU = ( B+SQRT( B*B+FOUR*C ) ) / TWO
        ELSE
          TAU = TWO*C / ( -B+SQRT( B*B+FOUR*C ) )
        END IF
!
!        The following TAU is DSIGMA - D( 2 )
!
        TAU = TAU / ( D( 2 )+SQRT( D( 2 )*D( 2 )+TAU ) )
        DSIGMA = D( 2 ) + TAU
        DELTA( 1 ) = -( DEL+TAU )
        DELTA( 2 ) = -TAU
        WORK( 1 ) = D( 1 ) + TAU + D( 2 )
        WORK( 2 ) = TWO*D( 2 ) + TAU
!        DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
!        DELTA( 2 ) = -Z( 2 ) / TAU
!        TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
!        DELTA( 1 ) = DELTA( 1 ) / TEMP
!        DELTA( 2 ) = DELTA( 2 ) / TEMP
      END IF
      RETURN
!
!     End of DLASD5
!
      END
      SUBROUTINE DLASD6( ICOMPQ, NL, NR, SQRE, D, VF, VL, ALPHA, BETA,  &
     &IDXQ, PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM,                       &
     &LDGNUM, POLES, DIFL, DIFR, Z, K, C, S, WORK,                      &
     &IWORK, INFO )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            GIVPTR, ICOMPQ, INFO, K, LDGCOL, LDGNUM, NL,   &
     &NR, SQRE
      DOUBLE PRECISION   ALPHA, BETA, C, S
!     ..
!     .. Array Arguments ..
      INTEGER            GIVCOL( LDGCOL, * ), IDXQ( * ), IWORK( * ),    &
     &PERM( * )
      DOUBLE PRECISION   D( * ), DIFL( * ), DIFR( * ),                  &
     &GIVNUM( LDGNUM, * ), POLES( LDGNUM, * ),                          &
     &VF( * ), VL( * ), WORK( * ), Z( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASD6 computes the SVD of an updated upper bidiagonal matrix B
!  obtained by merging two smaller ones by appending a row. This
!  routine is used only for the problem which requires all singular
!  values and optionally singular vector matrices in factored form.
!  B is an N-by-M matrix with N = NL + NR + 1 and M = N + SQRE.
!  A related subroutine, DLASD1, handles the case in which all singular
!  values and singular vectors of the bidiagonal matrix are desired.
!
!  DLASD6 computes the SVD as follows:
!
!                ( D1(in)  0    0     0 )
!    B = U(in) * (   Z1'   a   Z2'    b ) * VT(in)
!                (   0     0   D2(in) 0 )
!
!      = U(out) * ( D(out) 0) * VT(out)
!
!  where Z' = (Z1' a Z2' b) = u' VT', and u is a vector of dimension M
!  with ALPHA and BETA in the NL+1 and NL+2 th entries and zeros
!  elsewhere; and the entry b is empty if SQRE = 0.
!
!  The singular values of B can be computed using D1, D2, the first
!  components of all the right singular vectors of the lower block, and
!  the last components of all the right singular vectors of the upper
!  block. These components are stored and updated in VF and VL,
!  respectively, in DLASD6. Hence U and VT are not explicitly
!  referenced.
!
!  The singular values are stored in D. The algorithm consists of two
!  stages:
!
!        The first stage consists of deflating the size of the problem
!        when there are multiple singular values or if there is a zero
!        in the Z vector. For each such occurence the dimension of the
!        secular equation problem is reduced by one. This stage is
!        performed by the routine DLASD7.
!
!        The second stage consists of calculating the updated
!        singular values. This is done by finding the roots of the
!        secular equation via the routine DLASD4 (as called by DLASD8).
!        This routine also updates VF and VL and computes the distances
!        between the updated singular values and the old singular
!        values.
!
!  DLASD6 is called from DLASDA.
!
!  Arguments
!  =========
!
!  ICOMPQ (input) INTEGER
!         Specifies whether singular vectors are to be computed in
!         factored form:
!         = 0: Compute singular values only.
!         = 1: Compute singular vectors in factored form as well.
!
!  NL     (input) INTEGER
!         The row dimension of the upper block.  NL >= 1.
!
!  NR     (input) INTEGER
!         The row dimension of the lower block.  NR >= 1.
!
!  SQRE   (input) INTEGER
!         = 0: the lower block is an NR-by-NR square matrix.
!         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.
!
!         The bidiagonal matrix has row dimension N = NL + NR + 1,
!         and column dimension M = N + SQRE.
!
!  D      (input/output) DOUBLE PRECISION array, dimension ( NL+NR+1 ).
!         On entry D(1:NL,1:NL) contains the singular values of the
!         upper block, and D(NL+2:N) contains the singular values
!         of the lower block. On exit D(1:N) contains the singular
!         values of the modified matrix.
!
!  VF     (input/output) DOUBLE PRECISION array, dimension ( M )
!         On entry, VF(1:NL+1) contains the first components of all
!         right singular vectors of the upper block; and VF(NL+2:M)
!         contains the first components of all right singular vectors
!         of the lower block. On exit, VF contains the first components
!         of all right singular vectors of the bidiagonal matrix.
!
!  VL     (input/output) DOUBLE PRECISION array, dimension ( M )
!         On entry, VL(1:NL+1) contains the  last components of all
!         right singular vectors of the upper block; and VL(NL+2:M)
!         contains the last components of all right singular vectors of
!         the lower block. On exit, VL contains the last components of
!         all right singular vectors of the bidiagonal matrix.
!
!  ALPHA  (input/output) DOUBLE PRECISION
!         Contains the diagonal element associated with the added row.
!
!  BETA   (input/output) DOUBLE PRECISION
!         Contains the off-diagonal element associated with the added
!         row.
!
!  IDXQ   (output) INTEGER array, dimension ( N )
!         This contains the permutation which will reintegrate the
!         subproblem just solved back into sorted order, i.e.
!         D( IDXQ( I = 1, N ) ) will be in ascending order.
!
!  PERM   (output) INTEGER array, dimension ( N )
!         The permutations (from deflation and sorting) to be applied
!         to each block. Not referenced if ICOMPQ = 0.
!
!  GIVPTR (output) INTEGER
!         The number of Givens rotations which took place in this
!         subproblem. Not referenced if ICOMPQ = 0.
!
!  GIVCOL (output) INTEGER array, dimension ( LDGCOL, 2 )
!         Each pair of numbers indicates a pair of columns to take place
!         in a Givens rotation. Not referenced if ICOMPQ = 0.
!
!  LDGCOL (input) INTEGER
!         leading dimension of GIVCOL, must be at least N.
!
!  GIVNUM (output) DOUBLE PRECISION array, dimension ( LDGNUM, 2 )
!         Each number indicates the C or S value to be used in the
!         corresponding Givens rotation. Not referenced if ICOMPQ = 0.
!
!  LDGNUM (input) INTEGER
!         The leading dimension of GIVNUM and POLES, must be at least N.
!
!  POLES  (output) DOUBLE PRECISION array, dimension ( LDGNUM, 2 )
!         On exit, POLES(1,*) is an array containing the new singular
!         values obtained from solving the secular equation, and
!         POLES(2,*) is an array containing the poles in the secular
!         equation. Not referenced if ICOMPQ = 0.
!
!  DIFL   (output) DOUBLE PRECISION array, dimension ( N )
!         On exit, DIFL(I) is the distance between I-th updated
!         (undeflated) singular value and the I-th (undeflated) old
!         singular value.
!
!  DIFR   (output) DOUBLE PRECISION array,
!                  dimension ( LDGNUM, 2 ) if ICOMPQ = 1 and
!                  dimension ( N ) if ICOMPQ = 0.
!         On exit, DIFR(I, 1) is the distance between I-th updated
!         (undeflated) singular value and the I+1-th (undeflated) old
!         singular value.
!
!         If ICOMPQ = 1, DIFR(1:K,2) is an array containing the
!         normalizing factors for the right singular vector matrix.
!
!         See DLASD8 for details on DIFL and DIFR.
!
!  Z      (output) DOUBLE PRECISION array, dimension ( M )
!         The first elements of this array contain the components
!         of the deflation-adjusted updating row vector.
!
!  K      (output) INTEGER
!         Contains the dimension of the non-deflated matrix,
!         This is the order of the related secular equation. 1 <= K <=N.
!
!  C      (output) DOUBLE PRECISION
!         C contains garbage if SQRE =0 and the C-value of a Givens
!         rotation related to the right null space if SQRE = 1.
!
!  S      (output) DOUBLE PRECISION
!         S contains garbage if SQRE =0 and the S-value of a Givens
!         rotation related to the right null space if SQRE = 1.
!
!  WORK   (workspace) DOUBLE PRECISION array, dimension ( 4 * M )
!
!  IWORK  (workspace) INTEGER array, dimension ( 3 * N )
!
!  INFO   (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = 1, an singular value did not converge
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ming Gu and Huan Ren, Computer Science Division, University of
!     California at Berkeley, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IDX, IDXC, IDXP, ISIGMA, IVFW, IVLW, IW, M, &
     &N, N1, N2
      DOUBLE PRECISION   ORGNRM
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPL, DLAMRG, DLASCL, DLASD7, DLASD8, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      N = NL + NR + 1
      M = N + SQRE
!
      IF( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) THEN
        INFO = -1
      ELSE IF( NL.LT.1 ) THEN
        INFO = -2
      ELSE IF( NR.LT.1 ) THEN
        INFO = -3
      ELSE IF( ( SQRE.LT.0 ) .OR. ( SQRE.GT.1 ) ) THEN
        INFO = -4
      ELSE IF( LDGCOL.LT.N ) THEN
        INFO = -14
      ELSE IF( LDGNUM.LT.N ) THEN
        INFO = -16
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DLASD6', -INFO )
        RETURN
      END IF
!
!     The following values are for bookkeeping purposes only.  They are
!     integer pointers which indicate the portion of the workspace
!     used by a particular array in DLASD7 and DLASD8.
!
      ISIGMA = 1
      IW = ISIGMA + N
      IVFW = IW + M
      IVLW = IVFW + M
!
      IDX = 1
      IDXC = IDX + N
      IDXP = IDXC + N
!
!     Scale.
!
      ORGNRM = MAX( ABS( ALPHA ), ABS( BETA ) )
      D( NL+1 ) = ZERO
      DO I = 1, N
        IF( ABS( D( I ) ).GT.ORGNRM ) THEN
          ORGNRM = ABS( D( I ) )
        END IF
      enddo
      CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, N, 1, D, N, INFO )
      ALPHA = ALPHA / ORGNRM
      BETA = BETA / ORGNRM
!
!     Sort and Deflate singular values.
!
      CALL DLASD7( ICOMPQ, NL, NR, SQRE, K, D, Z, WORK( IW ), VF,       &
     &WORK( IVFW ), VL, WORK( IVLW ), ALPHA, BETA,                      &
     &WORK( ISIGMA ), IWORK( IDX ), IWORK( IDXP ), IDXQ,                &
     &PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, C, S,               &
     &INFO )
!
!     Solve Secular Equation, compute DIFL, DIFR, and update VF, VL.
!
      CALL DLASD8( ICOMPQ, K, D, Z, VF, VL, DIFL, DIFR, LDGNUM,         &
     &WORK( ISIGMA ), WORK( IW ), INFO )
!
!     Save the poles if ICOMPQ = 1.
!
      IF( ICOMPQ.EQ.1 ) THEN
        CALL DCOPL( K, D, 1, POLES( 1, 1 ), 1 )
        CALL DCOPL( K, WORK( ISIGMA ), 1, POLES( 1, 2 ), 1 )
      END IF
!
!     Unscale.
!
      CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO )
!
!     Prepare the IDXQ sorting permutation.
!
      N1 = K
      N2 = N - K
      CALL DLAMRG( N1, N2, D, 1, -1, IDXQ )
!
      RETURN
!
!     End of DLASD6
!
      END
      SUBROUTINE DLASD7( ICOMPQ, NL, NR, SQRE, K, D, Z, ZW, VF, VFW, VL,&
     &VLW, ALPHA, BETA, DSIGMA, IDX, IDXP, IDXQ,                        &
     &PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM,                     &
     &C, S, INFO )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            GIVPTR, ICOMPQ, INFO, K, LDGCOL, LDGNUM, NL,   &
     &NR, SQRE
      DOUBLE PRECISION   ALPHA, BETA, C, S
!     ..
!     .. Array Arguments ..
      INTEGER            GIVCOL( LDGCOL, * ), IDX( * ), IDXP( * ),      &
     &IDXQ( * ), PERM( * )
      DOUBLE PRECISION   D( * ), DSIGMA( * ), GIVNUM( LDGNUM, * ),      &
     &VF( * ), VFW( * ), VL( * ), VLW( * ), Z( * ),                     &
     &ZW( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASD7 merges the two sets of singular values together into a single
!  sorted set. Then it tries to deflate the size of the problem. There
!  are two ways in which deflation can occur:  when two or more singular
!  values are close together or if there is a tiny entry in the Z
!  vector. For each such occurrence the order of the related
!  secular equation problem is reduced by one.
!
!  DLASD7 is called from DLASD6.
!
!  Arguments
!  =========
!
!  ICOMPQ  (input) INTEGER
!          Specifies whether singular vectors are to be computed
!          in compact form, as follows:
!          = 0: Compute singular values only.
!          = 1: Compute singular vectors of upper
!               bidiagonal matrix in compact form.
!
!  NL     (input) INTEGER
!         The row dimension of the upper block. NL >= 1.
!
!  NR     (input) INTEGER
!         The row dimension of the lower block. NR >= 1.
!
!  SQRE   (input) INTEGER
!         = 0: the lower block is an NR-by-NR square matrix.
!         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.
!
!         The bidiagonal matrix has
!         N = NL + NR + 1 rows and
!         M = N + SQRE >= N columns.
!
!  K      (output) INTEGER
!         Contains the dimension of the non-deflated matrix, this is
!         the order of the related secular equation. 1 <= K <=N.
!
!  D      (input/output) DOUBLE PRECISION array, dimension ( N )
!         On entry D contains the singular values of the two submatrices
!         to be combined. On exit D contains the trailing (N-K) updated
!         singular values (those which were deflated) sorted into
!         increasing order.
!
!  Z      (output) DOUBLE PRECISION array, dimension ( M )
!         On exit Z contains the updating row vector in the secular
!         equation.
!
!  ZW     (workspace) DOUBLE PRECISION array, dimension ( M )
!         Workspace for Z.
!
!  VF     (input/output) DOUBLE PRECISION array, dimension ( M )
!         On entry, VF(1:NL+1) contains the first components of all
!         right singular vectors of the upper block; and VF(NL+2:M)
!         contains the first components of all right singular vectors
!         of the lower block. On exit, VF contains the first components
!         of all right singular vectors of the bidiagonal matrix.
!
!  VFW    (workspace) DOUBLE PRECISION array, dimension ( M )
!         Workspace for VF.
!
!  VL     (input/output) DOUBLE PRECISION array, dimension ( M )
!         On entry, VL(1:NL+1) contains the  last components of all
!         right singular vectors of the upper block; and VL(NL+2:M)
!         contains the last components of all right singular vectors
!         of the lower block. On exit, VL contains the last components
!         of all right singular vectors of the bidiagonal matrix.
!
!  VLW    (workspace) DOUBLE PRECISION array, dimension ( M )
!         Workspace for VL.
!
!  ALPHA  (input) DOUBLE PRECISION
!         Contains the diagonal element associated with the added row.
!
!  BETA   (input) DOUBLE PRECISION
!         Contains the off-diagonal element associated with the added
!         row.
!
!  DSIGMA (output) DOUBLE PRECISION array, dimension ( N )
!         Contains a copy of the diagonal elements (K-1 singular values
!         and one zero) in the secular equation.
!
!  IDX    (workspace) INTEGER array, dimension ( N )
!         This will contain the permutation used to sort the contents of
!         D into ascending order.
!
!  IDXP   (workspace) INTEGER array, dimension ( N )
!         This will contain the permutation used to place deflated
!         values of D at the end of the array. On output IDXP(2:K)
!         points to the nondeflated D-values and IDXP(K+1:N)
!         points to the deflated singular values.
!
!  IDXQ   (input) INTEGER array, dimension ( N )
!         This contains the permutation which separately sorts the two
!         sub-problems in D into ascending order.  Note that entries in
!         the first half of this permutation must first be moved one
!         position backward; and entries in the second half
!         must first have NL+1 added to their values.
!
!  PERM   (output) INTEGER array, dimension ( N )
!         The permutations (from deflation and sorting) to be applied
!         to each singular block. Not referenced if ICOMPQ = 0.
!
!  GIVPTR (output) INTEGER
!         The number of Givens rotations which took place in this
!         subproblem. Not referenced if ICOMPQ = 0.
!
!  GIVCOL (output) INTEGER array, dimension ( LDGCOL, 2 )
!         Each pair of numbers indicates a pair of columns to take place
!         in a Givens rotation. Not referenced if ICOMPQ = 0.
!
!  LDGCOL (input) INTEGER
!         The leading dimension of GIVCOL, must be at least N.
!
!  GIVNUM (output) DOUBLE PRECISION array, dimension ( LDGNUM, 2 )
!         Each number indicates the C or S value to be used in the
!         corresponding Givens rotation. Not referenced if ICOMPQ = 0.
!
!  LDGNUM (input) INTEGER
!         The leading dimension of GIVNUM, must be at least N.
!
!  C      (output) DOUBLE PRECISION
!         C contains garbage if SQRE =0 and the C-value of a Givens
!         rotation related to the right null space if SQRE = 1.
!
!  S      (output) DOUBLE PRECISION
!         S contains garbage if SQRE =0 and the S-value of a Givens
!         rotation related to the right null space if SQRE = 1.
!
!  INFO   (output) INTEGER
!         = 0:  successful exit.
!         < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ming Gu and Huan Ren, Computer Science Division, University of
!     California at Berkeley, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, EIGHT
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0,   &
     &EIGHT = 8.0D+0 )
!     ..
!     .. Local Scalars ..
!
      INTEGER            I, IDXI, IDXJ, IDXJP, J, JP, JPREV, K2, M, N,  &
     &NLP1, NLP2
      DOUBLE PRECISION   EPS, HLFTOL, TAU, TOL, Z1
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPL, DLAMRG, DROT, XERBLA
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2
      EXTERNAL           DLAMCH, DLAPY2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      N = NL + NR + 1
      M = N + SQRE
!
      IF( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) THEN
        INFO = -1
      ELSE IF( NL.LT.1 ) THEN
        INFO = -2
      ELSE IF( NR.LT.1 ) THEN
        INFO = -3
      ELSE IF( ( SQRE.LT.0 ) .OR. ( SQRE.GT.1 ) ) THEN
        INFO = -4
      ELSE IF( LDGCOL.LT.N ) THEN
        INFO = -22
      ELSE IF( LDGNUM.LT.N ) THEN
        INFO = -24
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DLASD7', -INFO )
        RETURN
      END IF
!
      NLP1 = NL + 1
      NLP2 = NL + 2
      IF( ICOMPQ.EQ.1 ) THEN
        GIVPTR = 0
      END IF
!
!     Generate the first part of the vector Z and move the singular
!     values in the first part of D one position backward.
!
      Z1 = ALPHA*VL( NLP1 )
      VL( NLP1 ) = ZERO
      TAU = VF( NLP1 )
      DO I = NL, 1, -1
        Z( I+1 ) = ALPHA*VL( I )
        VL( I ) = ZERO
        VF( I+1 ) = VF( I )
        D( I+1 ) = D( I )
        IDXQ( I+1 ) = IDXQ( I ) + 1
      enddo
      VF( 1 ) = TAU
!
!     Generate the second part of the vector Z.
!
      DO I = NLP2, M
        Z( I ) = BETA*VF( I )
        VF( I ) = ZERO
      enddo
!
!     Sort the singular values into increasing order
!
      DO I = NLP2, N
        IDXQ( I ) = IDXQ( I ) + NLP1
      enddo
!
!     DSIGMA, IDXC, IDXC, and ZW are used as storage space.
!
      DO I = 2, N
        DSIGMA( I ) = D( IDXQ( I ) )
        ZW( I ) = Z( IDXQ( I ) )
        VFW( I ) = VF( IDXQ( I ) )
        VLW( I ) = VL( IDXQ( I ) )
      enddo
!
      CALL DLAMRG( NL, NR, DSIGMA( 2 ), 1, 1, IDX( 2 ) )
!
      DO I = 2, N
        IDXI = 1 + IDX( I )
        D( I ) = DSIGMA( IDXI )
        Z( I ) = ZW( IDXI )
        VF( I ) = VFW( IDXI )
        VL( I ) = VLW( IDXI )
      enddo
!
!     Calculate the allowable deflation tolerence
!
      EPS = DLAMCH( 'Epsilon' )
      TOL = MAX( ABS( ALPHA ), ABS( BETA ) )
      TOL = EIGHT*EIGHT*EPS*MAX( ABS( D( N ) ), TOL )
!
!     There are 2 kinds of deflation -- first a value in the z-vector
!     is small, second two (or more) singular values are very close
!     together (their difference is small).
!
!     If the value in the z-vector is small, we simply permute the
!     array so that the corresponding singular value is moved to the
!     end.
!
!     If two values in the D-vector are close, we perform a two-sided
!     rotation designed to make one of the corresponding z-vector
!     entries zero, and then permute the array so that the deflated
!     singular value is moved to the end.
!
!     If there are multiple singular values then the problem deflates.
!     Here the number of equal singular values are found.  As each equal
!     singular value is found, an elementary reflector is computed to
!     rotate the corresponding singular subspace so that the
!     corresponding components of Z are zero in this new basis.
!
      K = 1
      K2 = N + 1
      DO J = 2, N
        IF( ABS( Z( J ) ).LE.TOL ) THEN
!
!           Deflate due to small z component.
!
          K2 = K2 - 1
          IDXP( K2 ) = J
          IF( J.EQ.N )                                                  &
     &GO TO 100
        ELSE
          JPREV = J
          GO TO 70
        END IF
      enddo
   70 CONTINUE
      J = JPREV
   80 CONTINUE
      J = J + 1
      IF( J.GT.N )                                                      &
     &GO TO 90
      IF( ABS( Z( J ) ).LE.TOL ) THEN
!
!        Deflate due to small z component.
!
        K2 = K2 - 1
        IDXP( K2 ) = J
      ELSE
!
!        Check if singular values are close enough to allow deflation.
!
        IF( ABS( D( J )-D( JPREV ) ).LE.TOL ) THEN
!
!           Deflation is possible.
!
          S = Z( JPREV )
          C = Z( J )
!
!           Find sqrt(a**2+b**2) without overflow or
!           destructive underflow.
!
          TAU = DLAPY2( C, S )
          Z( J ) = TAU
          Z( JPREV ) = ZERO
          C = C / TAU
          S = -S / TAU
!
!           Record the appropriate Givens rotation
!
          IF( ICOMPQ.EQ.1 ) THEN
            GIVPTR = GIVPTR + 1
            IDXJP = IDXQ( IDX( JPREV )+1 )
            IDXJ = IDXQ( IDX( J )+1 )
            IF( IDXJP.LE.NLP1 ) THEN
              IDXJP = IDXJP - 1
            END IF
            IF( IDXJ.LE.NLP1 ) THEN
              IDXJ = IDXJ - 1
            END IF
            GIVCOL( GIVPTR, 2 ) = IDXJP
            GIVCOL( GIVPTR, 1 ) = IDXJ
            GIVNUM( GIVPTR, 2 ) = C
            GIVNUM( GIVPTR, 1 ) = S
          END IF
          CALL DROT( 1, VF( JPREV ), 1, VF( J ), 1, C, S )
          CALL DROT( 1, VL( JPREV ), 1, VL( J ), 1, C, S )
          K2 = K2 - 1
          IDXP( K2 ) = JPREV
          JPREV = J
        ELSE
          K = K + 1
          ZW( K ) = Z( JPREV )
          DSIGMA( K ) = D( JPREV )
          IDXP( K ) = JPREV
          JPREV = J
        END IF
      END IF
      GO TO 80
   90 CONTINUE
!
!     Record the last singular value.
!
      K = K + 1
      ZW( K ) = Z( JPREV )
      DSIGMA( K ) = D( JPREV )
      IDXP( K ) = JPREV
!
  100 CONTINUE
!
!     Sort the singular values into DSIGMA. The singular values which
!     were not deflated go into the first K slots of DSIGMA, except
!     that DSIGMA(1) is treated separately.
!
      DO J = 2, N
        JP = IDXP( J )
        DSIGMA( J ) = D( JP )
        VFW( J ) = VF( JP )
        VLW( J ) = VL( JP )
      enddo
      IF( ICOMPQ.EQ.1 ) THEN
        DO J = 2, N
          JP = IDXP( J )
          PERM( J ) = IDXQ( IDX( JP )+1 )
          IF( PERM( J ).LE.NLP1 ) THEN
            PERM( J ) = PERM( J ) - 1
          END IF
        enddo
      END IF
!
!     The deflated singular values go back into the last N - K slots of
!     D.
!
      CALL DCOPL( N-K, DSIGMA( K+1 ), 1, D( K+1 ), 1 )
!
!     Determine DSIGMA(1), DSIGMA(2), Z(1), VF(1), VL(1), VF(M), and
!     VL(M).
!
      DSIGMA( 1 ) = ZERO
      HLFTOL = TOL / TWO
      IF( ABS( DSIGMA( 2 ) ).LE.HLFTOL )                                &
     &DSIGMA( 2 ) = HLFTOL
      IF( M.GT.N ) THEN
        Z( 1 ) = DLAPY2( Z1, Z( M ) )
        IF( Z( 1 ).LE.TOL ) THEN
          C = ONE
          S = ZERO
          Z( 1 ) = TOL
        ELSE
          C = Z1 / Z( 1 )
          S = -Z( M ) / Z( 1 )
        END IF
        CALL DROT( 1, VF( M ), 1, VF( 1 ), 1, C, S )
        CALL DROT( 1, VL( M ), 1, VL( 1 ), 1, C, S )
      ELSE
        IF( ABS( Z1 ).LE.TOL ) THEN
          Z( 1 ) = TOL
        ELSE
          Z( 1 ) = Z1
        END IF
      END IF
!
!     Restore Z, VF, and VL.
!
      CALL DCOPL( K-1, ZW( 2 ), 1, Z( 2 ), 1 )
      CALL DCOPL( N-1, VFW( 2 ), 1, VF( 2 ), 1 )
      CALL DCOPL( N-1, VLW( 2 ), 1, VL( 2 ), 1 )
!
      RETURN
!
!     End of DLASD7
!
      END
      SUBROUTINE DLASD8( ICOMPQ, K, D, Z, VF, VL, DIFL, DIFR, LDDIFR,   &
     &DSIGMA, WORK, INFO )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            ICOMPQ, INFO, K, LDDIFR
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), DIFL( * ), DIFR( LDDIFR, * ),          &
     &DSIGMA( * ), VF( * ), VL( * ), WORK( * ),                         &
     &Z( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASD8 finds the square roots of the roots of the secular equation,
!  as defined by the values in DSIGMA and Z. It makes the appropriate
!  calls to DLASD4, and stores, for each  element in D, the distance
!  to its two nearest poles (elements in DSIGMA). It also updates
!  the arrays VF and VL, the first and last components of all the
!  right singular vectors of the original bidiagonal matrix.
!
!  DLASD8 is called from DLASD6.
!
!  Arguments
!  =========
!
!  ICOMPQ  (input) INTEGER
!          Specifies whether singular vectors are to be computed in
!          factored form in the calling routine:
!          = 0: Compute singular values only.
!          = 1: Compute singular vectors in factored form as well.
!
!  K       (input) INTEGER
!          The number of terms in the rational function to be solved
!          by DLASD4.  K >= 1.
!
!  D       (output) DOUBLE PRECISION array, dimension ( K )
!          On output, D contains the updated singular values.
!
!  Z       (input) DOUBLE PRECISION array, dimension ( K )
!          The first K elements of this array contain the components
!          of the deflation-adjusted updating row vector.
!
!  VF      (input/output) DOUBLE PRECISION array, dimension ( K )
!          On entry, VF contains  information passed through DBEDE8.
!          On exit, VF contains the first K components of the first
!          components of all right singular vectors of the bidiagonal
!          matrix.
!
!  VL      (input/output) DOUBLE PRECISION array, dimension ( K )
!          On entry, VL contains  information passed through DBEDE8.
!          On exit, VL contains the first K components of the last
!          components of all right singular vectors of the bidiagonal
!          matrix.
!
!  DIFL    (output) DOUBLE PRECISION array, dimension ( K )
!          On exit, DIFL(I) = D(I) - DSIGMA(I).
!
!  DIFR    (output) DOUBLE PRECISION array,
!                   dimension ( LDDIFR, 2 ) if ICOMPQ = 1 and
!                   dimension ( K ) if ICOMPQ = 0.
!          On exit, DIFR(I,1) = D(I) - DSIGMA(I+1), DIFR(K,1) is not
!          defined and will not be referenced.
!
!          If ICOMPQ = 1, DIFR(1:K,2) is an array containing the
!          normalizing factors for the right singular vector matrix.
!
!  LDDIFR  (input) INTEGER
!          The leading dimension of DIFR, must be at least K.
!
!  DSIGMA  (input) DOUBLE PRECISION array, dimension ( K )
!          The first K elements of this array contain the old roots
!          of the deflated updating problem.  These are the poles
!          of the secular equation.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension at least 3 * K
!
!  INFO    (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = 1, an singular value did not converge
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ming Gu and Huan Ren, Computer Science Division, University of
!     California at Berkeley, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IWK1, IWK2, IWK2I, IWK3, IWK3I, J
      DOUBLE PRECISION   DIFLJ, DIFRJ, DJ, DSIGJ, DSIGJP, RHO, TEMP
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPL, DLASCL, DLASD4, DLASET, XERBLA
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DDOT, DLAMC3, DNRM2
      EXTERNAL           DDOT, DLAMC3, DNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) THEN
        INFO = -1
      ELSE IF( K.LT.1 ) THEN
        INFO = -2
      ELSE IF( LDDIFR.LT.K ) THEN
        INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DLASD8', -INFO )
        RETURN
      END IF
!
!     Quick return if possible
!
      IF( K.EQ.1 ) THEN
        D( 1 ) = ABS( Z( 1 ) )
        DIFL( 1 ) = D( 1 )
        IF( ICOMPQ.EQ.1 ) THEN
          DIFL( 2 ) = ONE
          DIFR( 1, 2 ) = ONE
        END IF
        RETURN
      END IF
!
!     Modify values DSIGMA(i) to make sure all DSIGMA(i)-DSIGMA(j) can
!     be computed with high relative accuracy (barring over/underflow).
!     This is a problem on machines without a guard digit in
!     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
!     The following code replaces DSIGMA(I) by 2*DSIGMA(I)-DSIGMA(I),
!     which on any of these machines zeros out the bottommost
!     bit of DSIGMA(I) if it is 1; this makes the subsequent
!     subtractions DSIGMA(I)-DSIGMA(J) unproblematic when cancellation
!     occurs. On binary machines with a guard digit (almost all
!     machines) it does not change DSIGMA(I) at all. On hexadecimal
!     and decimal machines with a guard digit, it slightly
!     changes the bottommost bits of DSIGMA(I). It does not account
!     for hexadecimal or decimal machines without guard digits
!     (we know of none). We use a subroutine call to compute
!     2*DSIGMA(I) to prevent optimizing compilers from eliminating
!     this code.
!
      DO I = 1, K
        DSIGMA( I ) = DLAMC3( DSIGMA( I ), DSIGMA( I ) ) - DSIGMA( I )
      enddo
!
!     Book keeping.
!
      IWK1 = 1
      IWK2 = IWK1 + K
      IWK3 = IWK2 + K
      IWK2I = IWK2 - 1
      IWK3I = IWK3 - 1
!
!     Normalize Z.
!
      RHO = DNRM2( K, Z, 1 )
      CALL DLASCL( 'G', 0, 0, RHO, ONE, K, 1, Z, K, INFO )
      RHO = RHO*RHO
!
!     Initialize WORK(IWK3).
!
      CALL DLASET( 'A', K, 1, ONE, ONE, WORK( IWK3 ), K )
!
!     Compute the updated singular values, the arrays DIFL, DIFR,
!     and the updated Z.
!
      DO J = 1, K
        CALL DLASD4( K, J, DSIGMA, Z, WORK( IWK1 ), RHO, D( J ),        &
     &WORK( IWK2 ), INFO )
!
!        If the root finder fails, the computation is terminated.
!
        IF( INFO.NE.0 ) THEN
          RETURN
        END IF
        WORK( IWK3I+J ) = WORK( IWK3I+J )*WORK( J )*WORK( IWK2I+J )
        DIFL( J ) = -WORK( J )
        DIFR( J, 1 ) = -WORK( J+1 )
        DO I = 1, J - 1
          WORK( IWK3I+I ) = WORK( IWK3I+I )*WORK( I )*                  &
     &WORK( IWK2I+I ) / ( DSIGMA( I )-                                  &
     &DSIGMA( J ) ) / ( DSIGMA( I )+                                    &
     &DSIGMA( J ) )
        enddo
        DO I = J + 1, K
          WORK( IWK3I+I ) = WORK( IWK3I+I )*WORK( I )*                  &
     &WORK( IWK2I+I ) / ( DSIGMA( I )-                                  &
     &DSIGMA( J ) ) / ( DSIGMA( I )+                                    &
     &DSIGMA( J ) )
        enddo
      enddo
!
!     Compute updated Z.
!
      DO I = 1, K
        Z( I ) = SIGN( SQRT( ABS( WORK( IWK3I+I ) ) ), Z( I ) )
      enddo
!
!     Update VF and VL.
!
      DO J = 1, K
        DIFLJ = DIFL( J )
        DJ = D( J )
        DSIGJ = -DSIGMA( J )
        IF( J.LT.K ) THEN
          DIFRJ = -DIFR( J, 1 )
          DSIGJP = -DSIGMA( J+1 )
        END IF
        WORK( J ) = -Z( J ) / DIFLJ / ( DSIGMA( J )+DJ )
        DO I = 1, J - 1
          WORK( I ) = Z( I ) / ( DLAMC3( DSIGMA( I ), DSIGJ )-DIFLJ )   &
     &/ ( DSIGMA( I )+DJ )
        enddo
        DO I = J + 1, K
          WORK( I ) = Z( I ) / ( DLAMC3( DSIGMA( I ), DSIGJP )+DIFRJ )  &
     &/ ( DSIGMA( I )+DJ )
        enddo
        TEMP = DNRM2( K, WORK, 1 )
        WORK( IWK2I+J ) = DDOT( K, WORK, 1, VF, 1 ) / TEMP
        WORK( IWK3I+J ) = DDOT( K, WORK, 1, VL, 1 ) / TEMP
        IF( ICOMPQ.EQ.1 ) THEN
          DIFR( J, 2 ) = TEMP
        END IF
      enddo
!
      CALL DCOPL( K, WORK( IWK2 ), 1, VF, 1 )
      CALL DCOPL( K, WORK( IWK3 ), 1, VL, 1 )
!
      RETURN
!
!     End of DLASD8
!
      END
      SUBROUTINE DLASDA( ICOMPQ, SMLSIZ, N, SQRE, D, E, U, LDU, VT, K,  &
     &DIFL, DIFR, Z, POLES, GIVPTR, GIVCOL, LDGCOL,                     &
     &PERM, GIVNUM, C, S, WORK, IWORK, INFO )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            ICOMPQ, INFO, LDGCOL, LDU, N, SMLSIZ, SQRE
!     ..
!     .. Array Arguments ..
      INTEGER            GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ),  &
     &K( * ), PERM( LDGCOL, * )
      DOUBLE PRECISION   C( * ), D( * ), DIFL( LDU, * ), DIFR( LDU, * ),&
     &E( * ), GIVNUM( LDU, * ), POLES( LDU, * ),                        &
     &S( * ), U( LDU, * ), VT( LDU, * ), WORK( * ),                     &
     &Z( LDU, * )
!     ..
!
!  Purpose
!  =======
!
!  Using a divide and conquer approach, DLASDA computes the singular
!  value decomposition (SVD) of a real upper bidiagonal N-by-M matrix
!  B with diagonal D and offdiagonal E, where M = N + SQRE. The
!  algorithm computes the singular values in the SVD B = U * S * VT.
!  The orthogonal matrices U and VT are optionally computed in
!  compact form.
!
!  A related subroutine, DLASD0, computes the singular values and
!  the singular vectors in explicit form.
!
!  Arguments
!  =========
!
!  ICOMPQ (input) INTEGER
!         Specifies whether singular vectors are to be computed
!         in compact form, as follows
!         = 0: Compute singular values only.
!         = 1: Compute singular vectors of upper bidiagonal
!              matrix in compact form.
!
!  SMLSIZ (input) INTEGER
!         The maximum size of the subproblems at the bottom of the
!         computation tree.
!
!  N      (input) INTEGER
!         The row dimension of the upper bidiagonal matrix. This is
!         also the dimension of the main diagonal array D.
!
!  SQRE   (input) INTEGER
!         Specifies the column dimension of the bidiagonal matrix.
!         = 0: The bidiagonal matrix has column dimension M = N;
!         = 1: The bidiagonal matrix has column dimension M = N + 1.
!
!  D      (input/output) DOUBLE PRECISION array, dimension ( N )
!         On entry D contains the main diagonal of the bidiagonal
!         matrix. On exit D, if INFO = 0, contains its singular values.
!
!  E      (input) DOUBLE PRECISION array, dimension ( M-1 )
!         Contains the subdiagonal entries of the bidiagonal matrix.
!         On exit, E has been destroyed.
!
!  U      (output) DOUBLE PRECISION array,
!         dimension ( LDU, SMLSIZ ) if ICOMPQ = 1, and not referenced
!         if ICOMPQ = 0. If ICOMPQ = 1, on exit, U contains the left
!         singular vector matrices of all subproblems at the bottom
!         level.
!
!  LDU    (input) INTEGER, LDU = > N.
!         The leading dimension of arrays U, VT, DIFL, DIFR, POLES,
!         GIVNUM, and Z.
!
!  VT     (output) DOUBLE PRECISION array,
!         dimension ( LDU, SMLSIZ+1 ) if ICOMPQ = 1, and not referenced
!         if ICOMPQ = 0. If ICOMPQ = 1, on exit, VT' contains the right
!         singular vector matrices of all subproblems at the bottom
!         level.
!
!  K      (output) INTEGER array,
!         dimension ( N ) if ICOMPQ = 1 and dimension 1 if ICOMPQ = 0.
!         If ICOMPQ = 1, on exit, K(I) is the dimension of the I-th
!         secular equation on the computation tree.
!
!  DIFL   (output) DOUBLE PRECISION array, dimension ( LDU, NLVL ),
!         where NLVL = floor(log_2 (N/SMLSIZ))).
!
!  DIFR   (output) DOUBLE PRECISION array,
!                  dimension ( LDU, 2 * NLVL ) if ICOMPQ = 1 and
!                  dimension ( N ) if ICOMPQ = 0.
!         If ICOMPQ = 1, on exit, DIFL(1:N, I) and DIFR(1:N, 2 * I - 1)
!         record distances between singular values on the I-th
!         level and singular values on the (I -1)-th level, and
!         DIFR(1:N, 2 * I ) contains the normalizing factors for
!         the right singular vector matrix. See DLASD8 for details.
!
!  Z      (output) DOUBLE PRECISION array,
!                  dimension ( LDU, NLVL ) if ICOMPQ = 1 and
!                  dimension ( N ) if ICOMPQ = 0.
!         The first K elements of Z(1, I) contain the components of
!         the deflation-adjusted updating row vector for subproblems
!         on the I-th level.
!
!  POLES  (output) DOUBLE PRECISION array,
!         dimension ( LDU, 2 * NLVL ) if ICOMPQ = 1, and not referenced
!         if ICOMPQ = 0. If ICOMPQ = 1, on exit, POLES(1, 2*I - 1) and
!         POLES(1, 2*I) contain  the new and old singular values
!         involved in the secular equations on the I-th level.
!
!  GIVPTR (output) INTEGER array,
!         dimension ( N ) if ICOMPQ = 1, and not referenced if
!         ICOMPQ = 0. If ICOMPQ = 1, on exit, GIVPTR( I ) records
!         the number of Givens rotations performed on the I-th
!         problem on the computation tree.
!
!  GIVCOL (output) INTEGER array,
!         dimension ( LDGCOL, 2 * NLVL ) if ICOMPQ = 1, and not
!         referenced if ICOMPQ = 0. If ICOMPQ = 1, on exit, for each I,
!         GIVCOL(1, 2 *I - 1) and GIVCOL(1, 2 *I) record the locations
!         of Givens rotations performed on the I-th level on the
!         computation tree.
!
!  LDGCOL (input) INTEGER, LDGCOL = > N.
!         The leading dimension of arrays GIVCOL and PERM.
!
!  PERM   (output) INTEGER array,
!         dimension ( LDGCOL, NLVL ) if ICOMPQ = 1, and not referenced
!         if ICOMPQ = 0. If ICOMPQ = 1, on exit, PERM(1, I) records
!         permutations done on the I-th level of the computation tree.
!
!  GIVNUM (output) DOUBLE PRECISION array,
!         dimension ( LDU,  2 * NLVL ) if ICOMPQ = 1, and not
!         referenced if ICOMPQ = 0. If ICOMPQ = 1, on exit, for each I,
!         GIVNUM(1, 2 *I - 1) and GIVNUM(1, 2 *I) record the C- and S-
!         values of Givens rotations performed on the I-th level on
!         the computation tree.
!
!  C      (output) DOUBLE PRECISION array,
!         dimension ( N ) if ICOMPQ = 1, and dimension 1 if ICOMPQ = 0.
!         If ICOMPQ = 1 and the I-th subproblem is not square, on exit,
!         C( I ) contains the C-value of a Givens rotation related to
!         the right null space of the I-th subproblem.
!
!  S      (output) DOUBLE PRECISION array, dimension ( N ) if
!         ICOMPQ = 1, and dimension 1 if ICOMPQ = 0. If ICOMPQ = 1
!         and the I-th subproblem is not square, on exit, S( I )
!         contains the S-value of a Givens rotation related to
!         the right null space of the I-th subproblem.
!
!  WORK   (workspace) DOUBLE PRECISION array, dimension
!         (6 * N + (SMLSIZ + 1)*(SMLSIZ + 1)).
!
!  IWORK  (workspace) INTEGER array.
!         Dimension must be at least (7 * N).
!
!  INFO   (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = 1, an singular value did not converge
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ming Gu and Huan Ren, Computer Science Division, University of
!     California at Berkeley, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I1, IC, IDXQ, IDXQI, IM1, INODE, ITEMP, IWK,&
     &J, LF, LL, LVL, LVL2, M, NCC, ND, NDB1, NDIML,                    &
     &NDIMR, NL, NLF, NLP1, NLVL, NR, NRF, NRP1, NRU,                   &
     &NWORK1, NWORK2, SMLSZP, SQREI, VF, VFI, VL, VLI
      DOUBLE PRECISION   ALPHA, BETA
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPL, DLASD6, DLASDQ, DLASDT, DLASET, XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) THEN
        INFO = -1
      ELSE IF( SMLSIZ.LT.3 ) THEN
        INFO = -2
      ELSE IF( N.LT.0 ) THEN
        INFO = -3
      ELSE IF( ( SQRE.LT.0 ) .OR. ( SQRE.GT.1 ) ) THEN
        INFO = -4
      ELSE IF( LDU.LT.( N+SQRE ) ) THEN
        INFO = -8
      ELSE IF( LDGCOL.LT.N ) THEN
        INFO = -17
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DLASDA', -INFO )
        RETURN
      END IF
!
      M = N + SQRE
!
!     If the input matrix is too small, call DLASDQ to find the SVD.
!
      IF( N.LE.SMLSIZ ) THEN
        IF( ICOMPQ.EQ.0 ) THEN
          CALL DLASDQ( 'U', SQRE, N, 0, 0, 0, D, E, VT, LDU, U, LDU,    &
     &U, LDU, WORK, INFO )
        ELSE
          CALL DLASDQ( 'U', SQRE, N, M, N, 0, D, E, VT, LDU, U, LDU,    &
     &U, LDU, WORK, INFO )
        END IF
        RETURN
      END IF
!
!     Book-keeping and  set up the computation tree.
!
      INODE = 1
      NDIML = INODE + N
      NDIMR = NDIML + N
      IDXQ = NDIMR + N
      IWK = IDXQ + N
!
      NCC = 0
      NRU = 0
!
      SMLSZP = SMLSIZ + 1
      VF = 1
      VL = VF + M
      NWORK1 = VL + M
      NWORK2 = NWORK1 + SMLSZP*SMLSZP
!
      CALL DLASDT( N, NLVL, ND, IWORK( INODE ), IWORK( NDIML ),         &
     &IWORK( NDIMR ), SMLSIZ )
!
!     for the nodes on bottom level of the tree, solve
!     their subproblems by DLASDQ.
!
      NDB1 = ( ND+1 ) / 2
      DO I = NDB1, ND
!
!        IC : center row of each node
!        NL : number of rows of left  subproblem
!        NR : number of rows of right subproblem
!        NLF: starting row of the left   subproblem
!        NRF: starting row of the right  subproblem
!
        I1 = I - 1
        IC = IWORK( INODE+I1 )
        NL = IWORK( NDIML+I1 )
        NLP1 = NL + 1
        NR = IWORK( NDIMR+I1 )
        NLF = IC - NL
        NRF = IC + 1
        IDXQI = IDXQ + NLF - 2
        VFI = VF + NLF - 1
        VLI = VL + NLF - 1
        SQREI = 1
        IF( ICOMPQ.EQ.0 ) THEN
          CALL DLASET( 'A', NLP1, NLP1, ZERO, ONE, WORK( NWORK1 ),      &
     &SMLSZP )
          CALL DLASDQ( 'U', SQREI, NL, NLP1, NRU, NCC, D( NLF ),        &
     &E( NLF ), WORK( NWORK1 ), SMLSZP,                                 &
     &WORK( NWORK2 ), NL, WORK( NWORK2 ), NL,                           &
     &WORK( NWORK2 ), INFO )
          ITEMP = NWORK1 + NL*SMLSZP
          CALL DCOPL( NLP1, WORK( NWORK1 ), 1, WORK( VFI ), 1 )
          CALL DCOPL( NLP1, WORK( ITEMP ), 1, WORK( VLI ), 1 )
        ELSE
          CALL DLASET( 'A', NL, NL, ZERO, ONE, U( NLF, 1 ), LDU )
          CALL DLASET( 'A', NLP1, NLP1, ZERO, ONE, VT( NLF, 1 ), LDU )
          CALL DLASDQ( 'U', SQREI, NL, NLP1, NL, NCC, D( NLF ),         &
     &E( NLF ), VT( NLF, 1 ), LDU, U( NLF, 1 ), LDU,                    &
     &U( NLF, 1 ), LDU, WORK( NWORK1 ), INFO )
          CALL DCOPL( NLP1, VT( NLF, 1 ), 1, WORK( VFI ), 1 )
          CALL DCOPL( NLP1, VT( NLF, NLP1 ), 1, WORK( VLI ), 1 )
        END IF
        IF( INFO.NE.0 ) THEN
          RETURN
        END IF
        DO J = 1, NL
          IWORK( IDXQI+J ) = J
        enddo
        IF( ( I.EQ.ND ) .AND. ( SQRE.EQ.0 ) ) THEN
          SQREI = 0
        ELSE
          SQREI = 1
        END IF
        IDXQI = IDXQI + NLP1
        VFI = VFI + NLP1
        VLI = VLI + NLP1
        NRP1 = NR + SQREI
        IF( ICOMPQ.EQ.0 ) THEN
          CALL DLASET( 'A', NRP1, NRP1, ZERO, ONE, WORK( NWORK1 ),      &
     &SMLSZP )
          CALL DLASDQ( 'U', SQREI, NR, NRP1, NRU, NCC, D( NRF ),        &
     &E( NRF ), WORK( NWORK1 ), SMLSZP,                                 &
     &WORK( NWORK2 ), NR, WORK( NWORK2 ), NR,                           &
     &WORK( NWORK2 ), INFO )
          ITEMP = NWORK1 + ( NRP1-1 )*SMLSZP
          CALL DCOPL( NRP1, WORK( NWORK1 ), 1, WORK( VFI ), 1 )
          CALL DCOPL( NRP1, WORK( ITEMP ), 1, WORK( VLI ), 1 )
        ELSE
          CALL DLASET( 'A', NR, NR, ZERO, ONE, U( NRF, 1 ), LDU )
          CALL DLASET( 'A', NRP1, NRP1, ZERO, ONE, VT( NRF, 1 ), LDU )
          CALL DLASDQ( 'U', SQREI, NR, NRP1, NR, NCC, D( NRF ),         &
     &E( NRF ), VT( NRF, 1 ), LDU, U( NRF, 1 ), LDU,                    &
     &U( NRF, 1 ), LDU, WORK( NWORK1 ), INFO )
          CALL DCOPL( NRP1, VT( NRF, 1 ), 1, WORK( VFI ), 1 )
          CALL DCOPL( NRP1, VT( NRF, NRP1 ), 1, WORK( VLI ), 1 )
        END IF
        IF( INFO.NE.0 ) THEN
          RETURN
        END IF
        DO J = 1, NR
          IWORK( IDXQI+J ) = J
        enddo
      enddo
!
!     Now conquer each subproblem bottom-up.
!
      J = 2**NLVL
      DO LVL = NLVL, 1, -1
        LVL2 = LVL*2 - 1
!
!        Find the first node LF and last node LL on
!        the current level LVL.
!
        IF( LVL.EQ.1 ) THEN
          LF = 1
          LL = 1
        ELSE
          LF = 2**( LVL-1 )
          LL = 2*LF - 1
        END IF
        DO I = LF, LL
          IM1 = I - 1
          IC = IWORK( INODE+IM1 )
          NL = IWORK( NDIML+IM1 )
          NR = IWORK( NDIMR+IM1 )
          NLF = IC - NL
          NRF = IC + 1
          IF( I.EQ.LL ) THEN
            SQREI = SQRE
          ELSE
            SQREI = 1
          END IF
          VFI = VF + NLF - 1
          VLI = VL + NLF - 1
          IDXQI = IDXQ + NLF - 1
          ALPHA = D( IC )
          BETA = E( IC )
          IF( ICOMPQ.EQ.0 ) THEN
            CALL DLASD6( ICOMPQ, NL, NR, SQREI, D( NLF ),               &
     &WORK( VFI ), WORK( VLI ), ALPHA, BETA,                            &
     &IWORK( IDXQI ), PERM, GIVPTR( 1 ), GIVCOL,                        &
     &LDGCOL, GIVNUM, LDU, POLES, DIFL, DIFR, Z,                        &
     &K( 1 ), C( 1 ), S( 1 ), WORK( NWORK1 ),                           &
     &IWORK( IWK ), INFO )
          ELSE
            J = J - 1
            CALL DLASD6( ICOMPQ, NL, NR, SQREI, D( NLF ),               &
     &WORK( VFI ), WORK( VLI ), ALPHA, BETA,                            &
     &IWORK( IDXQI ), PERM( NLF, LVL ),                                 &
     &GIVPTR( J ), GIVCOL( NLF, LVL2 ), LDGCOL,                         &
     &GIVNUM( NLF, LVL2 ), LDU,                                         &
     &POLES( NLF, LVL2 ), DIFL( NLF, LVL ),                             &
     &DIFR( NLF, LVL2 ), Z( NLF, LVL ), K( J ),                         &
     &C( J ), S( J ), WORK( NWORK1 ),                                   &
     &IWORK( IWK ), INFO )
          END IF
          IF( INFO.NE.0 ) THEN
            RETURN
          END IF
        enddo
      enddo
!
      RETURN
!
!     End of DLASDA
!
      END
      SUBROUTINE DLASDQ( UPLO, SQRE, N, NCVT, NRU, NCC, D, E, VT, LDVT, &
     &U, LDU, C, LDC, WORK, INFO )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU, SQRE
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), D( * ), E( * ), U( LDU, * ),      &
     &VT( LDVT, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASDQ computes the singular value decomposition (SVD) of a real
!  (upper or lower) bidiagonal matrix with diagonal D and offdiagonal
!  E, accumulating the transformations if desired. Letting B denote
!  the input bidiagonal matrix, the algorithm computes orthogonal
!  matrices Q and P such that B = Q * S * P' (P' denotes the transpose
!  of P). The singular values S are overwritten on D.
!
!  The input matrix U  is changed to U  * Q  if desired.
!  The input matrix VT is changed to P' * VT if desired.
!  The input matrix C  is changed to Q' * C  if desired.
!
!  See "Computing  Small Singular Values of Bidiagonal Matrices With
!  Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
!  LAPACK Working Note #3, for a detailed description of the algorithm.
!
!  Arguments
!  =========
!
!  UPLO  (input) CHARACTER*1
!        On entry, UPLO specifies whether the input bidiagonal matrix
!        is upper or lower bidiagonal, and wether it is square are
!        not.
!           UPLO = 'U' or 'u'   B is upper bidiagonal.
!           UPLO = 'L' or 'l'   B is lower bidiagonal.
!
!  SQRE  (input) INTEGER
!        = 0: then the input matrix is N-by-N.
!        = 1: then the input matrix is N-by-(N+1) if UPLU = 'U' and
!             (N+1)-by-N if UPLU = 'L'.
!
!        The bidiagonal matrix has
!        N = NL + NR + 1 rows and
!        M = N + SQRE >= N columns.
!
!  N     (input) INTEGER
!        On entry, N specifies the number of rows and columns
!        in the matrix. N must be at least 0.
!
!  NCVT  (input) INTEGER
!        On entry, NCVT specifies the number of columns of
!        the matrix VT. NCVT must be at least 0.
!
!  NRU   (input) INTEGER
!        On entry, NRU specifies the number of rows of
!        the matrix U. NRU must be at least 0.
!
!  NCC   (input) INTEGER
!        On entry, NCC specifies the number of columns of
!        the matrix C. NCC must be at least 0.
!
!  D     (input/output) DOUBLE PRECISION array, dimension (N)
!        On entry, D contains the diagonal entries of the
!        bidiagonal matrix whose SVD is desired. On normal exit,
!        D contains the singular values in ascending order.
!
!  E     (input/output) DOUBLE PRECISION array.
!        dimension is (N-1) if SQRE = 0 and N if SQRE = 1.
!        On entry, the entries of E contain the offdiagonal entries
!        of the bidiagonal matrix whose SVD is desired. On normal
!        exit, E will contain 0. If the algorithm does not converge,
!        D and E will contain the diagonal and superdiagonal entries
!        of a bidiagonal matrix orthogonally equivalent to the one
!        given as input.
!
!  VT    (input/output) DOUBLE PRECISION array, dimension (LDVT, NCVT)
!        On entry, contains a matrix which on exit has been
!        premultiplied by P', dimension N-by-NCVT if SQRE = 0
!        and (N+1)-by-NCVT if SQRE = 1 (not referenced if NCVT=0).
!
!  LDVT  (input) INTEGER
!        On entry, LDVT specifies the leading dimension of VT as
!        declared in the calling (sub) program. LDVT must be at
!        least 1. If NCVT is nonzero LDVT must also be at least N.
!
!  U     (input/output) DOUBLE PRECISION array, dimension (LDU, N)
!        On entry, contains a  matrix which on exit has been
!        postmultiplied by Q, dimension NRU-by-N if SQRE = 0
!        and NRU-by-(N+1) if SQRE = 1 (not referenced if NRU=0).
!
!  LDU   (input) INTEGER
!        On entry, LDU  specifies the leading dimension of U as
!        declared in the calling (sub) program. LDU must be at
!        least max( 1, NRU ) .
!
!  C     (input/output) DOUBLE PRECISION array, dimension (LDC, NCC)
!        On entry, contains an N-by-NCC matrix which on exit
!        has been premultiplied by Q'  dimension N-by-NCC if SQRE = 0
!        and (N+1)-by-NCC if SQRE = 1 (not referenced if NCC=0).
!
!  LDC   (input) INTEGER
!        On entry, LDC  specifies the leading dimension of C as
!        declared in the calling (sub) program. LDC must be at
!        least 1. If NCC is nonzero, LDC must also be at least N.
!
!  WORK  (workspace) DOUBLE PRECISION array, dimension (4*N)
!        Workspace. Only referenced if one of NCVT, NRU, or NCC is
!        nonzero, and if N is at least 2.
!
!  INFO  (output) INTEGER
!        On exit, a value of 0 indicates a successful exit.
!        If INFO < 0, argument number -INFO is illegal.
!        If INFO > 0, the algorithm did not converge, and INFO
!        specifies how many superdiagonals did not converge.
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ming Gu and Huan Ren, Computer Science Division, University of
!     California at Berkeley, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            ROTATE
      INTEGER            I, ISUB, IUPLO, J, NP1, SQRE1
      DOUBLE PRECISION   CS, R, SMIN, SN
!     ..
!     .. External Subroutines ..
      EXTERNAL           DBDSQR, DLARTG, DLASR, DSWAP, XERBLA
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IUPLO = 0
      IF( LSAME( UPLO, 'U' ) )                                          &
     &IUPLO = 1
      IF( LSAME( UPLO, 'L' ) )                                          &
     &IUPLO = 2
      IF( IUPLO.EQ.0 ) THEN
        INFO = -1
      ELSE IF( ( SQRE.LT.0 ) .OR. ( SQRE.GT.1 ) ) THEN
        INFO = -2
      ELSE IF( N.LT.0 ) THEN
        INFO = -3
      ELSE IF( NCVT.LT.0 ) THEN
        INFO = -4
      ELSE IF( NRU.LT.0 ) THEN
        INFO = -5
      ELSE IF( NCC.LT.0 ) THEN
        INFO = -6
      ELSE IF( ( NCVT.EQ.0 .AND. LDVT.LT.1 ) .OR.                       &
     &( NCVT.GT.0 .AND. LDVT.LT.MAX( 1, N ) ) ) THEN
        INFO = -10
      ELSE IF( LDU.LT.MAX( 1, NRU ) ) THEN
        INFO = -12
      ELSE IF( ( NCC.EQ.0 .AND. LDC.LT.1 ) .OR.                         &
     &( NCC.GT.0 .AND. LDC.LT.MAX( 1, N ) ) ) THEN
        INFO = -14
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DLASDQ', -INFO )
        RETURN
      END IF
      IF( N.EQ.0 )                                                      &
     &RETURN
!
!     ROTATE is true if any singular vectors desired, false otherwise
!
      ROTATE = ( NCVT.GT.0 ) .OR. ( NRU.GT.0 ) .OR. ( NCC.GT.0 )
      NP1 = N + 1
      SQRE1 = SQRE
!
!     If matrix non-square upper bidiagonal, rotate to be lower
!     bidiagonal.  The rotations are on the right.
!
      IF( ( IUPLO.EQ.1 ) .AND. ( SQRE1.EQ.1 ) ) THEN
        DO I = 1, N - 1
          CALL DLARTG( D( I ), E( I ), CS, SN, R )
          D( I ) = R
          E( I ) = SN*D( I+1 )
          D( I+1 ) = CS*D( I+1 )
          IF( ROTATE ) THEN
            WORK( I ) = CS
            WORK( N+I ) = SN
          END IF
        enddo
        CALL DLARTG( D( N ), E( N ), CS, SN, R )
        D( N ) = R
        E( N ) = ZERO
        IF( ROTATE ) THEN
          WORK( N ) = CS
          WORK( N+N ) = SN
        END IF
        IUPLO = 2
        SQRE1 = 0
!
!        Update singular vectors if desired.
!
        IF( NCVT.GT.0 )                                                 &
     &CALL DLASR( 'L', 'V', 'F', NP1, NCVT, WORK( 1 ),                  &
     &WORK( NP1 ), VT, LDVT )
      END IF
!
!     If matrix lower bidiagonal, rotate to be upper bidiagonal
!     by applying Givens rotations on the left.
!
      IF( IUPLO.EQ.2 ) THEN
        DO I = 1, N - 1
          CALL DLARTG( D( I ), E( I ), CS, SN, R )
          D( I ) = R
          E( I ) = SN*D( I+1 )
          D( I+1 ) = CS*D( I+1 )
          IF( ROTATE ) THEN
            WORK( I ) = CS
            WORK( N+I ) = SN
          END IF
        enddo
!
!        If matrix (N+1)-by-N lower bidiagonal, one additional
!        rotation is needed.
!
        IF( SQRE1.EQ.1 ) THEN
          CALL DLARTG( D( N ), E( N ), CS, SN, R )
          D( N ) = R
          IF( ROTATE ) THEN
            WORK( N ) = CS
            WORK( N+N ) = SN
          END IF
        END IF
!
!        Update singular vectors if desired.
!
        IF( NRU.GT.0 ) THEN
          IF( SQRE1.EQ.0 ) THEN
            CALL DLASR( 'R', 'V', 'F', NRU, N, WORK( 1 ),               &
     &WORK( NP1 ), U, LDU )
          ELSE
            CALL DLASR( 'R', 'V', 'F', NRU, NP1, WORK( 1 ),             &
     &WORK( NP1 ), U, LDU )
          END IF
        END IF
        IF( NCC.GT.0 ) THEN
          IF( SQRE1.EQ.0 ) THEN
            CALL DLASR( 'L', 'V', 'F', N, NCC, WORK( 1 ),               &
     &WORK( NP1 ), C, LDC )
          ELSE
            CALL DLASR( 'L', 'V', 'F', NP1, NCC, WORK( 1 ),             &
     &WORK( NP1 ), C, LDC )
          END IF
        END IF
      END IF
!
!     Call DBDSQR to compute the SVD of the reduced real
!     N-by-N upper bidiagonal matrix.
!
      CALL DBDSQR( 'U', N, NCVT, NRU, NCC, D, E, VT, LDVT, U, LDU, C,   &
     &LDC, WORK, INFO )
!
!     Sort the singular values into ascending order (insertion sort on
!     singular values, but only one transposition per singular vector)
!
      DO I = 1, N
!
!        Scan for smallest D(I).
!
        ISUB = I
        SMIN = D( I )
        DO J = I + 1, N
          IF( D( J ).LT.SMIN ) THEN
            ISUB = J
            SMIN = D( J )
          END IF
        enddo
        IF( ISUB.NE.I ) THEN
!
!           Swap singular values and vectors.
!
          D( ISUB ) = D( I )
          D( I ) = SMIN
          IF( NCVT.GT.0 )                                               &
     &CALL DSWAP( NCVT, VT( ISUB, 1 ), LDVT, VT( I, 1 ), LDVT )
          IF( NRU.GT.0 )                                                &
     &CALL DSWAP( NRU, U( 1, ISUB ), 1, U( 1, I ), 1 )
          IF( NCC.GT.0 )                                                &
     &CALL DSWAP( NCC, C( ISUB, 1 ), LDC, C( I, 1 ), LDC )
        END IF
      enddo
!
      RETURN
!
!     End of DLASDQ
!
      END
      SUBROUTINE DLASDT( N, LVL, ND, INODE, NDIML, NDIMR, MSUB )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            LVL, MSUB, N, ND
!     ..
!     .. Array Arguments ..
      INTEGER            INODE( * ), NDIML( * ), NDIMR( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASDT creates a tree of subproblems for bidiagonal divide and
!  conquer.
!
!  Arguments
!  =========
!
!   N      (input) INTEGER
!          On entry, the number of diagonal elements of the
!          bidiagonal matrix.
!
!   LVL    (output) INTEGER
!          On exit, the number of levels on the computation tree.
!
!   ND     (output) INTEGER
!          On exit, the number of nodes on the tree.
!
!   INODE  (output) INTEGER array, dimension ( N )
!          On exit, centers of subproblems.
!
!   NDIML  (output) INTEGER array, dimension ( N )
!          On exit, row dimensions of left children.
!
!   NDIMR  (output) INTEGER array, dimension ( N )
!          On exit, row dimensions of right children.
!
!   MSUB   (input) INTEGER.
!          On entry, the maximum row dimension each subproblem at the
!          bottom of the tree can be of.
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ming Gu and Huan Ren, Computer Science Division, University of
!     California at Berkeley, USA
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IL, IR, LLST, MAXN, NCRNT, NLVL
      DOUBLE PRECISION   TEMP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, LOG, MAX
!     ..
!     .. Executable Statements ..
!
!     Find the number of levels on the tree.
!
      MAXN = MAX( 1, N )
      TEMP = LOG( DBLE( MAXN ) / DBLE( MSUB+1 ) ) / LOG( TWO )
      LVL = INT( TEMP ) + 1
!
      I = N / 2
      INODE( 1 ) = I + 1
      NDIML( 1 ) = I
      NDIMR( 1 ) = N - I - 1
      IL = 0
      IR = 1
      LLST = 1
      DO NLVL = 1, LVL - 1
!
!        Constructing the tree at (NLVL+1)-st level. The number of
!        nodes created on this level is LLST * 2.
!
        DO I = 0, LLST - 1
          IL = IL + 2
          IR = IR + 2
          NCRNT = LLST + I
          NDIML( IL ) = NDIML( NCRNT ) / 2
          NDIMR( IL ) = NDIML( NCRNT ) - NDIML( IL ) - 1
          INODE( IL ) = INODE( NCRNT ) - NDIMR( IL ) - 1
          NDIML( IR ) = NDIMR( NCRNT ) / 2
          NDIMR( IR ) = NDIMR( NCRNT ) - NDIML( IR ) - 1
          INODE( IR ) = INODE( NCRNT ) + NDIML( IR ) + 1
        enddo
        LLST = LLST*2
      enddo
      ND = LLST*2 - 1
!
      RETURN
!
!     End of DLASDT
!
      END
      SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      DOUBLE PRECISION   ALPHA, BETA
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DLASET initializes an m-by-n matrix A to BETA on the diagonal and
!  ALPHA on the offdiagonals.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies the part of the matrix A to be set.
!          = 'U':      Upper triangular part is set; the strictly lower
!                      triangular part of A is not changed.
!          = 'L':      Lower triangular part is set; the strictly upper
!                      triangular part of A is not changed.
!          Otherwise:  All of the matrix A is set.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  ALPHA   (input) DOUBLE PRECISION
!          The constant to which the offdiagonal elements are to be set.
!
!  BETA    (input) DOUBLE PRECISION
!          The constant to which the diagonal elements are to be set.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On exit, the leading m-by-n submatrix of A is set as follows:
!
!          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
!          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
!          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
!
!          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Set the strictly upper triangular or trapezoidal part of the
!        array to ALPHA.
!
        DO J = 2, N
          DO I = 1, MIN( J-1, M )
            A( I, J ) = ALPHA
          enddo
        enddo
!
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
!
!        Set the strictly lower triangular or trapezoidal part of the
!        array to ALPHA.
!
        DO J = 1, MIN( M, N )
          DO I = J + 1, M
            A( I, J ) = ALPHA
          enddo
        enddo
!
      ELSE
!
!        Set the leading m-by-n submatrix to ALPHA.
!
        DO J = 1, N
          DO I = 1, M
            A( I, J ) = ALPHA
          enddo
        enddo
      END IF
!
!     Set the first min(M,N) diagonal elements to BETA.
!
      DO I = 1, MIN( M, N )
        A( I, I ) = BETA
      enddo
!
      RETURN
!
!     End of DLASET
!
      END
      SUBROUTINE DLASQ1( N, D, E, WORK, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASQ1 computes the singular values of a real N-by-N bidiagonal
!  matrix with diagonal D and off-diagonal E. The singular values
!  are computed to high relative accuracy, in the absence of
!  denormalization, underflow and overflow. The algorithm was first
!  presented in
!
!  "Accurate singular values and differential qd algorithms" by K. V.
!  Fernando and B. N. Parlett, Numer. Math., Vol-67, No. 2, pp. 191-230,
!  1994,
!
!  and the present implementation is described in "An implementation of
!  the dqds Algorithm (Positive Case)", LAPACK Working Note.
!
!  Arguments
!  =========
!
!  N     (input) INTEGER
!        The number of rows and columns in the matrix. N >= 0.
!
!  D     (input/output) DOUBLE PRECISION array, dimension (N)
!        On entry, D contains the diagonal elements of the
!        bidiagonal matrix whose SVD is desired. On normal exit,
!        D contains the singular values in decreasing order.
!
!  E     (input/output) DOUBLE PRECISION array, dimension (N)
!        On entry, elements E(1:N-1) contain the off-diagonal elements
!        of the bidiagonal matrix whose SVD is desired.
!        On exit, E is overwritten.
!
!  WORK  (workspace) DOUBLE PRECISION array, dimension (4*N)
!
!  INFO  (output) INTEGER
!        = 0: successful exit
!        < 0: if INFO = -i, the i-th argument had an illegal value
!        > 0: the algorithm failed
!             = 1, a split was marked by a positive value in E
!             = 2, current block of Z not diagonalized after 30*N
!                  iterations (in inner while loop)
!             = 3, termination criterion of outer while loop not met
!                  (program created more than N unreduced blocks)
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IINFO
      DOUBLE PRECISION   EPS, SCALE, SAFMIN, SIGMN, SIGMX
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPL, DLAS2, DLASCL, DLASQ2, DLASRT, XERBLA
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
      INFO = 0
      IF( N.LT.0 ) THEN
        INFO = -2
        CALL XERBLA( 'DLASQ1', -INFO )
        RETURN
      ELSE IF( N.EQ.0 ) THEN
        RETURN
      ELSE IF( N.EQ.1 ) THEN
        D( 1 ) = ABS( D( 1 ) )
        RETURN
      ELSE IF( N.EQ.2 ) THEN
        CALL DLAS2( D( 1 ), E( 1 ), D( 2 ), SIGMN, SIGMX )
        D( 1 ) = SIGMX
        D( 2 ) = SIGMN
        RETURN
      END IF
!
!     Estimate the largest singular value.
!
      SIGMX = ZERO
      DO I = 1, N - 1
        D( I ) = ABS( D( I ) )
        SIGMX = MAX( SIGMX, ABS( E( I ) ) )
      enddo
      D( N ) = ABS( D( N ) )
!
!     Early return if SIGMX is zero (matrix is already diagonal).
!
      IF( SIGMX.EQ.ZERO ) THEN
        CALL DLASRT( 'D', N, D, IINFO )
        RETURN
      END IF
!
      DO I = 1, N
        SIGMX = MAX( SIGMX, D( I ) )
      enddo
!
!     Copy D and E into WORK (in the Z format) and scale (squaring the
!     input data makes scaling by a power of the radix pointless).
!
      EPS = DLAMCH( 'Precision' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      SCALE = SQRT( EPS / SAFMIN )
      CALL DCOPL( N, D, 1, WORK( 1 ), 2 )
      CALL DCOPL( N-1, E, 1, WORK( 2 ), 2 )
      CALL DLASCL( 'G', 0, 0, SIGMX, SCALE, 2*N-1, 1, WORK, 2*N-1,      &
     &IINFO )
!
!     Compute the q's and e's.
!
      DO I = 1, 2*N - 1
        WORK( I ) = WORK( I )**2
      enddo
      WORK( 2*N ) = ZERO
!
      CALL DLASQ2( N, WORK, INFO )
!
      IF( INFO.EQ.0 ) THEN
        DO I = 1, N
          D( I ) = SQRT( WORK( I ) )
        enddo
        CALL DLASCL( 'G', 0, 0, SCALE, SIGMX, N, 1, D, N, IINFO )
      END IF
!
      RETURN
!
!     End of DLASQ1
!
      END
      SUBROUTINE DLASQ2( N, Z, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     Modified to call DLAZQ3 in place of DLASQ3, 13 Feb 03, SJH.
!
!     .. Scalar Arguments ..
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   Z( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASQ2 computes all the eigenvalues of the symmetric positive
!  definite tridiagonal matrix associated with the qd array Z to high
!  relative accuracy are computed to high relative accuracy, in the
!  absence of denormalization, underflow and overflow.
!
!  To see the relation of Z to the tridiagonal matrix, let L be a
!  unit lower bidiagonal matrix with subdiagonals Z(2,4,6,,..) and
!  let U be an upper bidiagonal matrix with 1's above and diagonal
!  Z(1,3,5,,..). The tridiagonal is L*U or, if you prefer, the
!  symmetric tridiagonal to which it is similar.
!
!  Note : DLASQ2 defines a logical variable, IEEE, which is true
!  on machines which follow ieee-754 floating-point standard in their
!  handling of infinities and NaNs, and false otherwise. This variable
!  is passed to DLAZQ3.
!
!  Arguments
!  =========
!
!  N     (input) INTEGER
!        The number of rows and columns in the matrix. N >= 0.
!
!  Z     (workspace) DOUBLE PRECISION array, dimension ( 4*N )
!        On entry Z holds the qd array. On exit, entries 1 to N hold
!        the eigenvalues in decreasing order, Z( 2*N+1 ) holds the
!        trace, and Z( 2*N+2 ) holds the sum of the eigenvalues. If
!        N > 2, then Z( 2*N+3 ) holds the iteration count, Z( 2*N+4 )
!        holds NDIVS/NIN^2, and Z( 2*N+5 ) holds the percentage of
!        shifts that failed.
!
!  INFO  (output) INTEGER
!        = 0: successful exit
!        < 0: if the i-th argument is a scalar and had an illegal
!             value, then INFO = -i, if the i-th argument is an
!             array and the j-entry had an illegal value, then
!             INFO = -(i*100+j)
!        > 0: the algorithm failed
!              = 1, a split was marked by a positive value in E
!              = 2, current block of Z not diagonalized after 30*N
!                   iterations (in inner while loop)
!              = 3, termination criterion of outer while loop not met
!                   (program created more than N unreduced blocks)
!
!  Further Details
!  ===============
!  Local Variables: I0:N0 defines a current unreduced segment of Z.
!  The shifts are accumulated in SIGMA. Iteration count is in ITER.
!  Ping-pong is controlled by PP (alternates between 0 and 1).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   CBIAS
      PARAMETER          ( CBIAS = 1.50D0 )
      DOUBLE PRECISION   ZERO, HALF, ONE, TWO, FOUR, HUNDRD
      PARAMETER          ( ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0,     &
     &TWO = 2.0D0, FOUR = 4.0D0, HUNDRD = 100.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            IEEE
      INTEGER            I0, I4, IINFO, IPN4, ITER, IWHILA, IWHILB, K,  &
     &N0, NBIG, NDIV, NFAIL, PP, SPLT, TTYPE
      DOUBLE PRECISION   D, DESIG, DMIN, DMIN1, DMIN2, DN, DN1, DN2, E, &
     &EMAX, EMIN, EPS, OLDEMN, QMAX, QMIN, S, SAFMIN,                   &
     &SIGMA, T, TAU, TEMP, TOL, TOL2, TRACE, ZMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAZQ3, DLASRT, XERBLA
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH, ILAENV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
!     (in case DLASQ2 is not called by DLASQ1)
!
      INFO = 0
      EPS = DLAMCH( 'Precision' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      TOL = EPS*HUNDRD
      TOL2 = TOL**2
!
      IF( N.LT.0 ) THEN
        INFO = -1
        CALL XERBLA( 'DLASQ2', 1 )
        RETURN
      ELSE IF( N.EQ.0 ) THEN
        RETURN
      ELSE IF( N.EQ.1 ) THEN
!
!        1-by-1 case.
!
        IF( Z( 1 ).LT.ZERO ) THEN
          INFO = -201
          CALL XERBLA( 'DLASQ2', 2 )
        END IF
        RETURN
      ELSE IF( N.EQ.2 ) THEN
!
!        2-by-2 case.
!
        IF( Z( 2 ).LT.ZERO .OR. Z( 3 ).LT.ZERO ) THEN
          INFO = -2
          CALL XERBLA( 'DLASQ2', 2 )
          RETURN
        ELSE IF( Z( 3 ).GT.Z( 1 ) ) THEN
          D = Z( 3 )
          Z( 3 ) = Z( 1 )
          Z( 1 ) = D
        END IF
        Z( 5 ) = Z( 1 ) + Z( 2 ) + Z( 3 )
        IF( Z( 2 ).GT.Z( 3 )*TOL2 ) THEN
          T = HALF*( ( Z( 1 )-Z( 3 ) )+Z( 2 ) )
          S = Z( 3 )*( Z( 2 ) / T )
          IF( S.LE.T ) THEN
            S = Z( 3 )*( Z( 2 ) / ( T*( ONE+SQRT( ONE+S / T ) ) ) )
          ELSE
            S = Z( 3 )*( Z( 2 ) / ( T+SQRT( T )*SQRT( T+S ) ) )
          END IF
          T = Z( 1 ) + ( S+Z( 2 ) )
          Z( 3 ) = Z( 3 )*( Z( 1 ) / T )
          Z( 1 ) = T
        END IF
        Z( 2 ) = Z( 3 )
        Z( 6 ) = Z( 2 ) + Z( 1 )
        RETURN
      END IF
!
!     Check for negative data and compute sums of q's and e's.
!
      Z( 2*N ) = ZERO
      EMIN = Z( 2 )
      QMAX = ZERO
      ZMAX = ZERO
      D = ZERO
      E = ZERO
!
      DO K = 1, 2*( N-1 ), 2
        IF( Z( K ).LT.ZERO ) THEN
          INFO = -( 200+K )
          CALL XERBLA( 'DLASQ2', 2 )
          RETURN
        ELSE IF( Z( K+1 ).LT.ZERO ) THEN
          INFO = -( 200+K+1 )
          CALL XERBLA( 'DLASQ2', 2 )
          RETURN
        END IF
        D = D + Z( K )
        E = E + Z( K+1 )
        QMAX = MAX( QMAX, Z( K ) )
        EMIN = MIN( EMIN, Z( K+1 ) )
        ZMAX = MAX( QMAX, ZMAX, Z( K+1 ) )
      enddo
      IF( Z( 2*N-1 ).LT.ZERO ) THEN
        INFO = -( 200+2*N-1 )
        CALL XERBLA( 'DLASQ2', 2 )
        RETURN
      END IF
      D = D + Z( 2*N-1 )
      QMAX = MAX( QMAX, Z( 2*N-1 ) )
      ZMAX = MAX( QMAX, ZMAX )
!
!     Check for diagonality.
!
      IF( E.EQ.ZERO ) THEN
        DO K = 2, N
          Z( K ) = Z( 2*K-1 )
        enddo
        CALL DLASRT( 'D', N, Z, IINFO )
        Z( 2*N-1 ) = D
        RETURN
      END IF
!
      TRACE = D + E
!
!     Check for zero data.
!
      IF( TRACE.EQ.ZERO ) THEN
        Z( 2*N-1 ) = ZERO
        RETURN
      END IF
!
!     Check whether the machine is IEEE conformable.
!
      IEEE = ILAENV( 10, 'DLASQ2', 'N', 1, 2, 3, 4 ).EQ.1 .AND.         &
     &ILAENV( 11, 'DLASQ2', 'N', 1, 2, 3, 4 ).EQ.1
!
!     Rearrange data for locality: Z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...).
!
      DO K = 2*N, 2, -2
        Z( 2*K ) = ZERO
        Z( 2*K-1 ) = Z( K )
        Z( 2*K-2 ) = ZERO
        Z( 2*K-3 ) = Z( K-1 )
      enddo
!
      I0 = 1
      N0 = N
!
!     Reverse the qd-array, if warranted.
!
      IF( CBIAS*Z( 4*I0-3 ).LT.Z( 4*N0-3 ) ) THEN
        IPN4 = 4*( I0+N0 )
        DO I4 = 4*I0, 2*( I0+N0-1 ), 4
          TEMP = Z( I4-3 )
          Z( I4-3 ) = Z( IPN4-I4-3 )
          Z( IPN4-I4-3 ) = TEMP
          TEMP = Z( I4-1 )
          Z( I4-1 ) = Z( IPN4-I4-5 )
          Z( IPN4-I4-5 ) = TEMP
        enddo
      END IF
!
!     Initial split checking via dqd and Li's test.
!
      PP = 0
!
      DO K = 1, 2
!
        D = Z( 4*N0+PP-3 )
        DO I4 = 4*( N0-1 ) + PP, 4*I0 + PP, -4
          IF( Z( I4-1 ).LE.TOL2*D ) THEN
            Z( I4-1 ) = -ZERO
            D = Z( I4-3 )
          ELSE
            D = Z( I4-3 )*( D / ( D+Z( I4-1 ) ) )
          END IF
        enddo
!
!        dqd maps Z to ZZ plus Li's test.
!
        EMIN = Z( 4*I0+PP+1 )
        D = Z( 4*I0+PP-3 )
        DO I4 = 4*I0 + PP, 4*( N0-1 ) + PP, 4
          Z( I4-2*PP-2 ) = D + Z( I4-1 )
          IF( Z( I4-1 ).LE.TOL2*D ) THEN
            Z( I4-1 ) = -ZERO
            Z( I4-2*PP-2 ) = D
            Z( I4-2*PP ) = ZERO
            D = Z( I4+1 )
          ELSE IF( SAFMIN*Z( I4+1 ).LT.Z( I4-2*PP-2 ) .AND.             &
     &SAFMIN*Z( I4-2*PP-2 ).LT.Z( I4+1 ) ) THEN
            TEMP = Z( I4+1 ) / Z( I4-2*PP-2 )
            Z( I4-2*PP ) = Z( I4-1 )*TEMP
            D = D*TEMP
          ELSE
            Z( I4-2*PP ) = Z( I4+1 )*( Z( I4-1 ) / Z( I4-2*PP-2 ) )
            D = Z( I4+1 )*( D / Z( I4-2*PP-2 ) )
          END IF
          EMIN = MIN( EMIN, Z( I4-2*PP ) )
        enddo
        Z( 4*N0-PP-2 ) = D
!
!        Now find qmax.
!
        QMAX = Z( 4*I0-PP-2 )
        DO I4 = 4*I0 - PP + 2, 4*N0 - PP - 2, 4
          QMAX = MAX( QMAX, Z( I4 ) )
        enddo
!
!        Prepare for the next iteration on K.
!
        PP = 1 - PP
      enddo
!
!     Initialise variables to pass to DLAZQ3
!
      TTYPE = 0
      DMIN1 = ZERO
      DMIN2 = ZERO
      DN    = ZERO
      DN1   = ZERO
      DN2   = ZERO
      TAU   = ZERO
!
      ITER = 2
      NFAIL = 0
      NDIV = 2*( N0-I0 )
!
      DO IWHILA = 1, N + 1
        IF( N0.LT.1 )                                                   &
     &GO TO 150
!
!        While array unfinished do
!
!        E(N0) holds the value of SIGMA when submatrix in I0:N0
!        splits from the rest of the array, but is negated.
!
        DESIG = ZERO
        IF( N0.EQ.N ) THEN
          SIGMA = ZERO
        ELSE
          SIGMA = -Z( 4*N0-1 )
        END IF
        IF( SIGMA.LT.ZERO ) THEN
          INFO = 1
          RETURN
        END IF
!
!        Find last unreduced submatrix's top index I0, find QMAX and
!        EMIN. Find Gershgorin-type bound if Q's much greater than E's.
!
        EMAX = ZERO
        IF( N0.GT.I0 ) THEN
          EMIN = ABS( Z( 4*N0-5 ) )
        ELSE
          EMIN = ZERO
        END IF
        QMIN = Z( 4*N0-3 )
        QMAX = QMIN
        DO I4 = 4*N0, 8, -4
          IF( Z( I4-5 ).LE.ZERO )                                       &
     &GO TO 100
          IF( QMIN.GE.FOUR*EMAX ) THEN
            QMIN = MIN( QMIN, Z( I4-3 ) )
            EMAX = MAX( EMAX, Z( I4-5 ) )
          END IF
          QMAX = MAX( QMAX, Z( I4-7 )+Z( I4-5 ) )
          EMIN = MIN( EMIN, Z( I4-5 ) )
        enddo
        I4 = 4
!
  100   CONTINUE
        I0 = I4 / 4
!
!        Store EMIN for passing to DLAZQ3.
!
        Z( 4*N0-1 ) = EMIN
!
!        Put -(initial shift) into DMIN.
!
        DMIN = -MAX( ZERO, QMIN-TWO*SQRT( QMIN )*SQRT( EMAX ) )
!
!        Now I0:N0 is unreduced. PP = 0 for ping, PP = 1 for pong.
!
        PP = 0
!
        NBIG = 30*( N0-I0+1 )
        DO IWHILB = 1, NBIG
          IF( I0.GT.N0 )                                                &
     &GO TO 130
!
!           While submatrix unfinished take a good dqds step.
!
          CALL DLAZQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL,  &
     &ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1,                   &
     &DN2, TAU )
!
          PP = 1 - PP
!
!           When EMIN is very small check for splits.
!
          IF( PP.EQ.0 .AND. N0-I0.GE.3 ) THEN
            IF( Z( 4*N0 ).LE.TOL2*QMAX .OR.                             &
     &Z( 4*N0-1 ).LE.TOL2*SIGMA ) THEN
              SPLT = I0 - 1
              QMAX = Z( 4*I0-3 )
              EMIN = Z( 4*I0-1 )
              OLDEMN = Z( 4*I0 )
              DO I4 = 4*I0, 4*( N0-3 ), 4
                IF( Z( I4 ).LE.TOL2*Z( I4-3 ) .OR.                      &
     &Z( I4-1 ).LE.TOL2*SIGMA ) THEN
                  Z( I4-1 ) = -SIGMA
                  SPLT = I4 / 4
                  QMAX = ZERO
                  EMIN = Z( I4+3 )
                  OLDEMN = Z( I4+4 )
                ELSE
                  QMAX = MAX( QMAX, Z( I4+1 ) )
                  EMIN = MIN( EMIN, Z( I4-1 ) )
                  OLDEMN = MIN( OLDEMN, Z( I4 ) )
                END IF
              enddo
              Z( 4*N0-1 ) = EMIN
              Z( 4*N0 ) = OLDEMN
              I0 = SPLT + 1
            END IF
          END IF
!
        enddo
!
        INFO = 2
        RETURN
!
!        end IWHILB
!
  130   CONTINUE
!
      enddo
!
      INFO = 3
      RETURN
!
!     end IWHILA
!
  150 CONTINUE
!
!     Move q's to the front.
!
      DO K = 2, N
        Z( K ) = Z( 4*K-3 )
      enddo
!
!     Sort and compute sum of eigenvalues.
!
      CALL DLASRT( 'D', N, Z, IINFO )
!
      E = ZERO
      DO K = N, 1, -1
        E = E + Z( K )
      enddo
!
!     Store trace, sum(eigenvalues) and information on performance.
!
      Z( 2*N+1 ) = TRACE
      Z( 2*N+2 ) = E
      Z( 2*N+3 ) = DBLE( ITER )
      Z( 2*N+4 ) = DBLE( NDIV ) / DBLE( N**2 )
      Z( 2*N+5 ) = HUNDRD*NFAIL / DBLE( ITER )
      RETURN
!
!     End of DLASQ2
!
      END
      SUBROUTINE DLASQ5( I0, N0, Z, PP, TAU, DMIN, DMIN1, DMIN2, DN,    &
     &DNM1, DNM2, IEEE )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            I0, N0, PP
      DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DNM1, DNM2, TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   Z( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASQ5 computes one dqds transform in ping-pong form, one
!  version for IEEE machines another for non IEEE machines.
!
!  Arguments
!  =========
!
!  I0    (input) INTEGER
!        First index.
!
!  N0    (input) INTEGER
!        Last index.
!
!  Z     (input) DOUBLE PRECISION array, dimension ( 4*N )
!        Z holds the qd array. EMIN is stored in Z(4*N0) to avoid
!        an extra argument.
!
!  PP    (input) INTEGER
!        PP=0 for ping, PP=1 for pong.
!
!  TAU   (input) DOUBLE PRECISION
!        This is the shift.
!
!  DMIN  (output) DOUBLE PRECISION
!        Minimum value of d.
!
!  DMIN1 (output) DOUBLE PRECISION
!        Minimum value of d, excluding D( N0 ).
!
!  DMIN2 (output) DOUBLE PRECISION
!        Minimum value of d, excluding D( N0 ) and D( N0-1 ).
!
!  DN    (output) DOUBLE PRECISION
!        d(N0), the last value of d.
!
!  DNM1  (output) DOUBLE PRECISION
!        d(N0-1).
!
!  DNM2  (output) DOUBLE PRECISION
!        d(N0-2).
!
!  IEEE  (input) LOGICAL
!        Flag for IEEE or non IEEE arithmetic.
!
!  =====================================================================
!
!     .. Parameter ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            J4, J4P2
      DOUBLE PRECISION   D, EMIN, TEMP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
      IF( ( N0-I0-1 ).LE.0 )                                            &
     &RETURN
!
      J4 = 4*I0 + PP - 3
      EMIN = Z( J4+4 )
      D = Z( J4 ) - TAU
      DMIN = D
      DMIN1 = -Z( J4 )
!
      IF( IEEE ) THEN
!
!        Code for IEEE arithmetic.
!
        IF( PP.EQ.0 ) THEN
          DO J4 = 4*I0, 4*( N0-3 ), 4
            Z( J4-2 ) = D + Z( J4-1 )
            TEMP = Z( J4+1 ) / Z( J4-2 )
            D = D*TEMP - TAU
            DMIN = MIN( DMIN, D )
            Z( J4 ) = Z( J4-1 )*TEMP
            EMIN = MIN( Z( J4 ), EMIN )
          enddo
        ELSE
          DO J4 = 4*I0, 4*( N0-3 ), 4
            Z( J4-3 ) = D + Z( J4 )
            TEMP = Z( J4+2 ) / Z( J4-3 )
            D = D*TEMP - TAU
            DMIN = MIN( DMIN, D )
            Z( J4-1 ) = Z( J4 )*TEMP
            EMIN = MIN( Z( J4-1 ), EMIN )
          enddo
        END IF
!
!        Unroll last two steps.
!
        DNM2 = D
        DMIN2 = DMIN
        J4 = 4*( N0-2 ) - PP
        J4P2 = J4 + 2*PP - 1
        Z( J4-2 ) = DNM2 + Z( J4P2 )
        Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
        DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU
        DMIN = MIN( DMIN, DNM1 )
!
        DMIN1 = DMIN
        J4 = J4 + 4
        J4P2 = J4 + 2*PP - 1
        Z( J4-2 ) = DNM1 + Z( J4P2 )
        Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
        DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU
        DMIN = MIN( DMIN, DN )
!
      ELSE
!
!        Code for non IEEE arithmetic.
!
        IF( PP.EQ.0 ) THEN
          DO J4 = 4*I0, 4*( N0-3 ), 4
            Z( J4-2 ) = D + Z( J4-1 )
            IF( D.LT.ZERO ) THEN
              RETURN
            ELSE
              Z( J4 ) = Z( J4+1 )*( Z( J4-1 ) / Z( J4-2 ) )
              D = Z( J4+1 )*( D / Z( J4-2 ) ) - TAU
            END IF
            DMIN = MIN( DMIN, D )
            EMIN = MIN( EMIN, Z( J4 ) )
          enddo
        ELSE
          DO J4 = 4*I0, 4*( N0-3 ), 4
            Z( J4-3 ) = D + Z( J4 )
            IF( D.LT.ZERO ) THEN
              RETURN
            ELSE
              Z( J4-1 ) = Z( J4+2 )*( Z( J4 ) / Z( J4-3 ) )
              D = Z( J4+2 )*( D / Z( J4-3 ) ) - TAU
            END IF
            DMIN = MIN( DMIN, D )
            EMIN = MIN( EMIN, Z( J4-1 ) )
          enddo
        END IF
!
!        Unroll last two steps.
!
        DNM2 = D
        DMIN2 = DMIN
        J4 = 4*( N0-2 ) - PP
        J4P2 = J4 + 2*PP - 1
        Z( J4-2 ) = DNM2 + Z( J4P2 )
        IF( DNM2.LT.ZERO ) THEN
          RETURN
        ELSE
          Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
          DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU
        END IF
        DMIN = MIN( DMIN, DNM1 )
!
        DMIN1 = DMIN
        J4 = J4 + 4
        J4P2 = J4 + 2*PP - 1
        Z( J4-2 ) = DNM1 + Z( J4P2 )
        IF( DNM1.LT.ZERO ) THEN
          RETURN
        ELSE
          Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
          DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU
        END IF
        DMIN = MIN( DMIN, DN )
!
      END IF
!
      Z( J4+2 ) = DN
      Z( 4*N0-PP ) = EMIN
      RETURN
!
!     End of DLASQ5
!
      END
      SUBROUTINE DLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN,         &
     &DNM1, DNM2 )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            I0, N0, PP
      DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DNM1, DNM2
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   Z( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASQ6 computes one dqd (shift equal to zero) transform in
!  ping-pong form, with protection against underflow and overflow.
!
!  Arguments
!  =========
!
!  I0    (input) INTEGER
!        First index.
!
!  N0    (input) INTEGER
!        Last index.
!
!  Z     (input) DOUBLE PRECISION array, dimension ( 4*N )
!        Z holds the qd array. EMIN is stored in Z(4*N0) to avoid
!        an extra argument.
!
!  PP    (input) INTEGER
!        PP=0 for ping, PP=1 for pong.
!
!  DMIN  (output) DOUBLE PRECISION
!        Minimum value of d.
!
!  DMIN1 (output) DOUBLE PRECISION
!        Minimum value of d, excluding D( N0 ).
!
!  DMIN2 (output) DOUBLE PRECISION
!        Minimum value of d, excluding D( N0 ) and D( N0-1 ).
!
!  DN    (output) DOUBLE PRECISION
!        d(N0), the last value of d.
!
!  DNM1  (output) DOUBLE PRECISION
!        d(N0-1).
!
!  DNM2  (output) DOUBLE PRECISION
!        d(N0-2).
!
!  =====================================================================
!
!     .. Parameter ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            J4, J4P2
      DOUBLE PRECISION   D, EMIN, SAFMIN, TEMP
!     ..
!     .. External Function ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
      IF( ( N0-I0-1 ).LE.0 )                                            &
     &RETURN
!
      SAFMIN = DLAMCH( 'Safe minimum' )
      J4 = 4*I0 + PP - 3
      EMIN = Z( J4+4 )
      D = Z( J4 )
      DMIN = D
!
      IF( PP.EQ.0 ) THEN
        DO J4 = 4*I0, 4*( N0-3 ), 4
          Z( J4-2 ) = D + Z( J4-1 )
          IF( Z( J4-2 ).EQ.ZERO ) THEN
            Z( J4 ) = ZERO
            D = Z( J4+1 )
            DMIN = D
            EMIN = ZERO
          ELSE IF( SAFMIN*Z( J4+1 ).LT.Z( J4-2 ) .AND.                  &
     &SAFMIN*Z( J4-2 ).LT.Z( J4+1 ) ) THEN
            TEMP = Z( J4+1 ) / Z( J4-2 )
            Z( J4 ) = Z( J4-1 )*TEMP
            D = D*TEMP
          ELSE
            Z( J4 ) = Z( J4+1 )*( Z( J4-1 ) / Z( J4-2 ) )
            D = Z( J4+1 )*( D / Z( J4-2 ) )
          END IF
          DMIN = MIN( DMIN, D )
          EMIN = MIN( EMIN, Z( J4 ) )
        enddo
      ELSE
        DO J4 = 4*I0, 4*( N0-3 ), 4
          Z( J4-3 ) = D + Z( J4 )
          IF( Z( J4-3 ).EQ.ZERO ) THEN
            Z( J4-1 ) = ZERO
            D = Z( J4+2 )
            DMIN = D
            EMIN = ZERO
          ELSE IF( SAFMIN*Z( J4+2 ).LT.Z( J4-3 ) .AND.                  &
     &SAFMIN*Z( J4-3 ).LT.Z( J4+2 ) ) THEN
            TEMP = Z( J4+2 ) / Z( J4-3 )
            Z( J4-1 ) = Z( J4 )*TEMP
            D = D*TEMP
          ELSE
            Z( J4-1 ) = Z( J4+2 )*( Z( J4 ) / Z( J4-3 ) )
            D = Z( J4+2 )*( D / Z( J4-3 ) )
          END IF
          DMIN = MIN( DMIN, D )
          EMIN = MIN( EMIN, Z( J4-1 ) )
        enddo
      END IF
!
!     Unroll last two steps.
!
      DNM2 = D
      DMIN2 = DMIN
      J4 = 4*( N0-2 ) - PP
      J4P2 = J4 + 2*PP - 1
      Z( J4-2 ) = DNM2 + Z( J4P2 )
      IF( Z( J4-2 ).EQ.ZERO ) THEN
        Z( J4 ) = ZERO
        DNM1 = Z( J4P2+2 )
        DMIN = DNM1
        EMIN = ZERO
      ELSE IF( SAFMIN*Z( J4P2+2 ).LT.Z( J4-2 ) .AND.                    &
     &SAFMIN*Z( J4-2 ).LT.Z( J4P2+2 ) ) THEN
        TEMP = Z( J4P2+2 ) / Z( J4-2 )
        Z( J4 ) = Z( J4P2 )*TEMP
        DNM1 = DNM2*TEMP
      ELSE
        Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
        DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) )
      END IF
      DMIN = MIN( DMIN, DNM1 )
!
      DMIN1 = DMIN
      J4 = J4 + 4
      J4P2 = J4 + 2*PP - 1
      Z( J4-2 ) = DNM1 + Z( J4P2 )
      IF( Z( J4-2 ).EQ.ZERO ) THEN
        Z( J4 ) = ZERO
        DN = Z( J4P2+2 )
        DMIN = DN
        EMIN = ZERO
      ELSE IF( SAFMIN*Z( J4P2+2 ).LT.Z( J4-2 ) .AND.                    &
     &SAFMIN*Z( J4-2 ).LT.Z( J4P2+2 ) ) THEN
        TEMP = Z( J4P2+2 ) / Z( J4-2 )
        Z( J4 ) = Z( J4P2 )*TEMP
        DN = DNM1*TEMP
      ELSE
        Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
        DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) )
      END IF
      DMIN = MIN( DMIN, DN )
!
      Z( J4+2 ) = DN
      Z( 4*N0-PP ) = EMIN
      RETURN
!
!     End of DLASQ6
!
      END
      SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASR applies a sequence of plane rotations to a real matrix A,
!  from either the left or the right.
!
!  When SIDE = 'L', the transformation takes the form
!
!     A := P*A
!
!  and when SIDE = 'R', the transformation takes the form
!
!     A := A*P**T
!
!  where P is an orthogonal matrix consisting of a sequence of z plane
!  rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',
!  and P**T is the transpose of P.
!
!  When DIRECT = 'F' (Forward sequence), then
!
!     P = P(z-1) * ... * P(2) * P(1)
!
!  and when DIRECT = 'B' (Backward sequence), then
!
!     P = P(1) * P(2) * ... * P(z-1)
!
!  where P(k) is a plane rotation matrix defined by the 2-by-2 rotation
!
!     R(k) = (  c(k)  s(k) )
!          = ( -s(k)  c(k) ).
!
!  When PIVOT = 'V' (Variable pivot), the rotation is performed
!  for the plane (k,k+1), i.e., P(k) has the form
!
!     P(k) = (  1                                            )
!            (       ...                                     )
!            (              1                                )
!            (                   c(k)  s(k)                  )
!            (                  -s(k)  c(k)                  )
!            (                                1              )
!            (                                     ...       )
!            (                                            1  )
!
!  where R(k) appears as a rank-2 modification to the identity matrix in
!  rows and columns k and k+1.
!
!  When PIVOT = 'T' (Top pivot), the rotation is performed for the
!  plane (1,k+1), so P(k) has the form
!
!     P(k) = (  c(k)                    s(k)                 )
!            (         1                                     )
!            (              ...                              )
!            (                     1                         )
!            ( -s(k)                    c(k)                 )
!            (                                 1             )
!            (                                      ...      )
!            (                                             1 )
!
!  where R(k) appears in rows and columns 1 and k+1.
!
!  Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is
!  performed for the plane (k,z), giving P(k) the form
!
!     P(k) = ( 1                                             )
!            (      ...                                      )
!            (             1                                 )
!            (                  c(k)                    s(k) )
!            (                         1                     )
!            (                              ...              )
!            (                                     1         )
!            (                 -s(k)                    c(k) )
!
!  where R(k) appears in rows and columns k and z.  The rotations are
!  performed without ever forming P(k) explicitly.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          Specifies whether the plane rotation matrix P is applied to
!          A on the left or the right.
!          = 'L':  Left, compute A := P*A
!          = 'R':  Right, compute A:= A*P**T
!
!  PIVOT   (input) CHARACTER*1
!          Specifies the plane for which P(k) is a plane rotation
!          matrix.
!          = 'V':  Variable pivot, the plane (k,k+1)
!          = 'T':  Top pivot, the plane (1,k+1)
!          = 'B':  Bottom pivot, the plane (k,z)
!
!  DIRECT  (input) CHARACTER*1
!          Specifies whether P is a forward or backward sequence of
!          plane rotations.
!          = 'F':  Forward, P = P(z-1)*...*P(2)*P(1)
!          = 'B':  Backward, P = P(1)*P(2)*...*P(z-1)
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  If m <= 1, an immediate
!          return is effected.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  If n <= 1, an
!          immediate return is effected.
!
!  C       (input) DOUBLE PRECISION array, dimension
!                  (M-1) if SIDE = 'L'
!                  (N-1) if SIDE = 'R'
!          The cosines c(k) of the plane rotations.
!
!  S       (input) DOUBLE PRECISION array, dimension
!                  (M-1) if SIDE = 'L'
!                  (N-1) if SIDE = 'R'
!          The sines s(k) of the plane rotations.  The 2-by-2 plane
!          rotation part of the matrix P(k), R(k), has the form
!          R(k) = (  c(k)  s(k) )
!                 ( -s(k)  c(k) ).
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          The M-by-N matrix A.  On exit, A is overwritten by P*A if
!          SIDE = 'R' or by A*P**T if SIDE = 'L'.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, INFO, J
      DOUBLE PRECISION   CTEMP, STEMP, TEMP
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
        INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT,            &
     &'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
        INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) )&
     &THEN
        INFO = 3
      ELSE IF( M.LT.0 ) THEN
        INFO = 4
      ELSE IF( N.LT.0 ) THEN
        INFO = 5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = 9
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DLASR ', INFO )
        RETURN
      END IF
!
!     Quick return if possible
!
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )                                  &
     &RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
!
!        Form  P * A
!
        IF( LSAME( PIVOT, 'V' ) ) THEN
          IF( LSAME( DIRECT, 'F' ) ) THEN
            DO J = 1, M - 1
              CTEMP = C( J )
              STEMP = S( J )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                DO I = 1, N
                  TEMP = A( J+1, I )
                  A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                  A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
                enddo
              END IF
            enddo
          ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
            DO J = M - 1, 1, -1
              CTEMP = C( J )
              STEMP = S( J )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                DO I = 1, N
                  TEMP = A( J+1, I )
                  A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                  A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
                enddo
              END IF
            enddo
          END IF
        ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
          IF( LSAME( DIRECT, 'F' ) ) THEN
            DO J = 2, M
              CTEMP = C( J-1 )
              STEMP = S( J-1 )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                DO I = 1, N
                  TEMP = A( J, I )
                  A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                  A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
                enddo
              END IF
            enddo
          ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
            DO J = M, 2, -1
              CTEMP = C( J-1 )
              STEMP = S( J-1 )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                DO I = 1, N
                  TEMP = A( J, I )
                  A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                  A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
                enddo
              END IF
            enddo
          END IF
        ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
          IF( LSAME( DIRECT, 'F' ) ) THEN
            DO J = 1, M - 1
              CTEMP = C( J )
              STEMP = S( J )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                DO I = 1, N
                  TEMP = A( J, I )
                  A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                  A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
                enddo
              END IF
            enddo
          ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
            DO J = M - 1, 1, -1
              CTEMP = C( J )
              STEMP = S( J )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                DO I = 1, N
                  TEMP = A( J, I )
                  A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                  A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
                enddo
              END IF
            enddo
          END IF
        END IF
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!        Form A * P'
!
        IF( LSAME( PIVOT, 'V' ) ) THEN
          IF( LSAME( DIRECT, 'F' ) ) THEN
            DO J = 1, N - 1
              CTEMP = C( J )
              STEMP = S( J )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                DO I = 1, M
                  TEMP = A( I, J+1 )
                  A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
                enddo
              END IF
            enddo
          ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
            DO J = N - 1, 1, -1
              CTEMP = C( J )
              STEMP = S( J )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                DO I = 1, M
                  TEMP = A( I, J+1 )
                  A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
                enddo
              END IF
            enddo
          END IF
        ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
          IF( LSAME( DIRECT, 'F' ) ) THEN
            DO J = 2, N
              CTEMP = C( J-1 )
              STEMP = S( J-1 )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                DO I = 1, M
                  TEMP = A( I, J )
                  A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                  A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
                enddo
              END IF
            enddo
          ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
            DO J = N, 2, -1
              CTEMP = C( J-1 )
              STEMP = S( J-1 )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                DO I = 1, M
                  TEMP = A( I, J )
                  A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                  A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
                enddo
              END IF
            enddo
          END IF
        ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
          IF( LSAME( DIRECT, 'F' ) ) THEN
            DO J = 1, N - 1
              CTEMP = C( J )
              STEMP = S( J )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                DO I = 1, M
                  TEMP = A( I, J )
                  A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                  A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
                enddo
              END IF
            enddo
          ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
            DO J = N - 1, 1, -1
              CTEMP = C( J )
              STEMP = S( J )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                DO I = 1, M
                  TEMP = A( I, J )
                  A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                  A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
                enddo
              END IF
            enddo
          END IF
        END IF
      END IF
!
      RETURN
!
!     End of DLASR
!
      END
      SUBROUTINE DLASRT( ID, N, D, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * )
!     ..
!
!  Purpose
!  =======
!
!  Sort the numbers in D in increasing order (if ID = 'I') or
!  in decreasing order (if ID = 'D' ).
!
!  Use Quick Sort, reverting to Insertion sort on arrays of
!  size <= 20. Dimension of STACK limits N to about 2**32.
!
!  Arguments
!  =========
!
!  ID      (input) CHARACTER*1
!          = 'I': sort D in increasing order;
!          = 'D': sort D in decreasing order.
!
!  N       (input) INTEGER
!          The length of the array D.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the array to be sorted.
!          On exit, D has been sorted into increasing order
!          (D(1) <= ... <= D(N) ) or into decreasing order
!          (D(1) >= ... >= D(N) ), depending on ID.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            SELECT
      PARAMETER          ( SELECT = 20 )
!     ..
!     .. Local Scalars ..
      INTEGER            DIR, ENDD, I, J, START, STKPNT
      DOUBLE PRECISION   D1, D2, D3, DMNMX, TMP
!     ..
!     .. Local Arrays ..
      INTEGER            STACK( 2, 32 )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input paramters.
!
      INFO = 0
      DIR = -1
      IF( LSAME( ID, 'D' ) ) THEN
        DIR = 0
      ELSE IF( LSAME( ID, 'I' ) ) THEN
        DIR = 1
      END IF
      IF( DIR.EQ.-1 ) THEN
        INFO = -1
      ELSE IF( N.LT.0 ) THEN
        INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DLASRT', -INFO )
        RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.1 )                                                      &
     &RETURN
!
      STKPNT = 1
      STACK( 1, 1 ) = 1
      STACK( 2, 1 ) = N
   10 CONTINUE
      START = STACK( 1, STKPNT )
      ENDD = STACK( 2, STKPNT )
      STKPNT = STKPNT - 1
      IF( ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) THEN
!
!        Do Insertion sort on D( START:ENDD )
!
        IF( DIR.EQ.0 ) THEN
!
!           Sort into decreasing order
!
          DO I = START + 1, ENDD
            DO J = I, START + 1, -1
              IF( D( J ).GT.D( J-1 ) ) THEN
                DMNMX = D( J )
                D( J ) = D( J-1 )
                D( J-1 ) = DMNMX
              ELSE
                GO TO 30
              END IF
            enddo
   30       CONTINUE
          enddo
!
        ELSE
!
!           Sort into increasing order
!
          DO I = START + 1, ENDD
            DO J = I, START + 1, -1
              IF( D( J ).LT.D( J-1 ) ) THEN
                DMNMX = D( J )
                D( J ) = D( J-1 )
                D( J-1 ) = DMNMX
              ELSE
                GO TO 50
              END IF
            enddo
   50       CONTINUE
          enddo
!
        END IF
!
      ELSE IF( ENDD-START.GT.SELECT ) THEN
!
!        Partition D( START:ENDD ) and stack parts, largest one first
!
!        Choose partition entry as median of 3
!
        D1 = D( START )
        D2 = D( ENDD )
        I = ( START+ENDD ) / 2
        D3 = D( I )
        IF( D1.LT.D2 ) THEN
          IF( D3.LT.D1 ) THEN
            DMNMX = D1
          ELSE IF( D3.LT.D2 ) THEN
            DMNMX = D3
          ELSE
            DMNMX = D2
          END IF
        ELSE
          IF( D3.LT.D2 ) THEN
            DMNMX = D2
          ELSE IF( D3.LT.D1 ) THEN
            DMNMX = D3
          ELSE
            DMNMX = D1
          END IF
        END IF
!
        IF( DIR.EQ.0 ) THEN
!
!           Sort into decreasing order
!
          I = START - 1
          J = ENDD + 1
   60     CONTINUE
   70     CONTINUE
          J = J - 1
          IF( D( J ).LT.DMNMX )                                         &
     &GO TO 70
   80     CONTINUE
          I = I + 1
          IF( D( I ).GT.DMNMX )                                         &
     &GO TO 80
          IF( I.LT.J ) THEN
            TMP = D( I )
            D( I ) = D( J )
            D( J ) = TMP
            GO TO 60
          END IF
          IF( J-START.GT.ENDD-J-1 ) THEN
            STKPNT = STKPNT + 1
            STACK( 1, STKPNT ) = START
            STACK( 2, STKPNT ) = J
            STKPNT = STKPNT + 1
            STACK( 1, STKPNT ) = J + 1
            STACK( 2, STKPNT ) = ENDD
          ELSE
            STKPNT = STKPNT + 1
            STACK( 1, STKPNT ) = J + 1
            STACK( 2, STKPNT ) = ENDD
            STKPNT = STKPNT + 1
            STACK( 1, STKPNT ) = START
            STACK( 2, STKPNT ) = J
          END IF
        ELSE
!
!           Sort into increasing order
!
          I = START - 1
          J = ENDD + 1
   90     CONTINUE
  100     CONTINUE
          J = J - 1
          IF( D( J ).GT.DMNMX )                                         &
     &GO TO 100
  110     CONTINUE
          I = I + 1
          IF( D( I ).LT.DMNMX )                                         &
     &GO TO 110
          IF( I.LT.J ) THEN
            TMP = D( I )
            D( I ) = D( J )
            D( J ) = TMP
            GO TO 90
          END IF
          IF( J-START.GT.ENDD-J-1 ) THEN
            STKPNT = STKPNT + 1
            STACK( 1, STKPNT ) = START
            STACK( 2, STKPNT ) = J
            STKPNT = STKPNT + 1
            STACK( 1, STKPNT ) = J + 1
            STACK( 2, STKPNT ) = ENDD
          ELSE
            STKPNT = STKPNT + 1
            STACK( 1, STKPNT ) = J + 1
            STACK( 2, STKPNT ) = ENDD
            STKPNT = STKPNT + 1
            STACK( 1, STKPNT ) = START
            STACK( 2, STKPNT ) = J
          END IF
        END IF
      END IF
      IF( STKPNT.GT.0 )                                                 &
     &GO TO 10
      RETURN
!
!     End of DLASRT
!
      END
      SUBROUTINE DLASV2( F, G, H, SSMIN, SSMAX, SNR, CSR, SNL, CSL )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   CSL, CSR, F, G, H, SNL, SNR, SSMAX, SSMIN
!     ..
!
!  Purpose
!  =======
!
!  DLASV2 computes the singular value decomposition of a 2-by-2
!  triangular matrix
!     [  F   G  ]
!     [  0   H  ].
!  On return, abs(SSMAX) is the larger singular value, abs(SSMIN) is the
!  smaller singular value, and (CSL,SNL) and (CSR,SNR) are the left and
!  right singular vectors for abs(SSMAX), giving the decomposition
!
!     [ CSL  SNL ] [  F   G  ] [ CSR -SNR ]  =  [ SSMAX   0   ]
!     [-SNL  CSL ] [  0   H  ] [ SNR  CSR ]     [  0    SSMIN ].
!
!  Arguments
!  =========
!
!  F       (input) DOUBLE PRECISION
!          The (1,1) element of the 2-by-2 matrix.
!
!  G       (input) DOUBLE PRECISION
!          The (1,2) element of the 2-by-2 matrix.
!
!  H       (input) DOUBLE PRECISION
!          The (2,2) element of the 2-by-2 matrix.
!
!  SSMIN   (output) DOUBLE PRECISION
!          abs(SSMIN) is the smaller singular value.
!
!  SSMAX   (output) DOUBLE PRECISION
!          abs(SSMAX) is the larger singular value.
!
!  SNL     (output) DOUBLE PRECISION
!  CSL     (output) DOUBLE PRECISION
!          The vector (CSL, SNL) is a unit left singular vector for the
!          singular value abs(SSMAX).
!
!  SNR     (output) DOUBLE PRECISION
!  CSR     (output) DOUBLE PRECISION
!          The vector (CSR, SNR) is a unit right singular vector for the
!          singular value abs(SSMAX).
!
!  Further Details
!  ===============
!
!  Any input parameter may be aliased with any output parameter.
!
!  Barring over/underflow and assuming a guard digit in subtraction, all
!  output quantities are correct to within a few units in the last
!  place (ulps).
!
!  In IEEE arithmetic, the code works correctly if one matrix element is
!  infinite.
!
!  Overflow will not occur unless the largest singular value itself
!  overflows or is within a few ulps of overflow. (On machines with
!  partial overflow, like the Cray, overflow may occur if the largest
!  singular value is within a factor of 2 of overflow.)
!
!  Underflow is harmless if underflow is gradual. Otherwise, results
!  may correspond to a matrix modified by perturbations of size near
!  the underflow threshold.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   FOUR
      PARAMETER          ( FOUR = 4.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            GASMAL, SWAP
      INTEGER            PMAX
      DOUBLE PRECISION   A, CLT, CRT, D, FA, FT, GA, GT, HA, HT, L, M,  &
     &MM, R, S, SLT, SRT, T, TEMP, TSIGN, TT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. Executable Statements ..
!
      FT = F
      FA = ABS( FT )
      HT = H
      HA = ABS( H )
!
!     PMAX points to the maximum absolute element of matrix
!       PMAX = 1 if F largest in absolute values
!       PMAX = 2 if G largest in absolute values
!       PMAX = 3 if H largest in absolute values
!
      PMAX = 1
      SWAP = ( HA.GT.FA )
      IF( SWAP ) THEN
        PMAX = 3
        TEMP = FT
        FT = HT
        HT = TEMP
        TEMP = FA
        FA = HA
        HA = TEMP
!
!        Now FA .ge. HA
!
      END IF
      GT = G
      GA = ABS( GT )
      IF( GA.EQ.ZERO ) THEN
!
!        Diagonal matrix
!
        SSMIN = HA
        SSMAX = FA
        CLT = ONE
        CRT = ONE
        SLT = ZERO
        SRT = ZERO
      ELSE
        GASMAL = .TRUE.
        IF( GA.GT.FA ) THEN
          PMAX = 2
          IF( ( FA / GA ).LT.DLAMCH( 'EPS' ) ) THEN
!
!              Case of very large GA
!
            GASMAL = .FALSE.
            SSMAX = GA
            IF( HA.GT.ONE ) THEN
              SSMIN = FA / ( GA / HA )
            ELSE
              SSMIN = ( FA / GA )*HA
            END IF
            CLT = ONE
            SLT = HT / GT
            SRT = ONE
            CRT = FT / GT
          END IF
        END IF
        IF( GASMAL ) THEN
!
!           Normal case
!
          D = FA - HA
          IF( D.EQ.FA ) THEN
!
!              Copes with infinite F or H
!
            L = ONE
          ELSE
            L = D / FA
          END IF
!
!           Note that 0 .le. L .le. 1
!
          M = GT / FT
!
!           Note that abs(M) .le. 1/macheps
!
          T = TWO - L
!
!           Note that T .ge. 1
!
          MM = M*M
          TT = T*T
          S = SQRT( TT+MM )
!
!           Note that 1 .le. S .le. 1 + 1/macheps
!
          IF( L.EQ.ZERO ) THEN
            R = ABS( M )
          ELSE
            R = SQRT( L*L+MM )
          END IF
!
!           Note that 0 .le. R .le. 1 + 1/macheps
!
          A = HALF*( S+R )
!
!           Note that 1 .le. A .le. 1 + abs(M)
!
          SSMIN = HA / A
          SSMAX = FA*A
          IF( MM.EQ.ZERO ) THEN
!
!              Note that M is very tiny
!
            IF( L.EQ.ZERO ) THEN
              T = SIGN( TWO, FT )*SIGN( ONE, GT )
            ELSE
              T = GT / SIGN( D, FT ) + M / T
            END IF
          ELSE
            T = ( M / ( S+T )+M / ( R+L ) )*( ONE+A )
          END IF
          L = SQRT( T*T+FOUR )
          CRT = TWO / L
          SRT = T / L
          CLT = ( CRT+SRT*M ) / A
          SLT = ( HT / FT )*SRT / A
        END IF
      END IF
      IF( SWAP ) THEN
        CSL = SRT
        SNL = CRT
        CSR = SLT
        SNR = CLT
      ELSE
        CSL = CLT
        SNL = SLT
        CSR = CRT
        SNR = SRT
      END IF
!
!     Correct signs of SSMAX and SSMIN
!
      IF( PMAX.EQ.1 )                                                   &
     &TSIGN = SIGN( ONE, CSR )*SIGN( ONE, CSL )*SIGN( ONE, F )
      IF( PMAX.EQ.2 )                                                   &
     &TSIGN = SIGN( ONE, SNR )*SIGN( ONE, CSL )*SIGN( ONE, G )
      IF( PMAX.EQ.3 )                                                   &
     &TSIGN = SIGN( ONE, SNR )*SIGN( ONE, SNL )*SIGN( ONE, H )
      SSMAX = SIGN( SSMAX, TSIGN )
      SSMIN = SIGN( SSMIN, TSIGN*SIGN( ONE, F )*SIGN( ONE, H ) )
      RETURN
!
!     End of DLASV2
!
      END
      SUBROUTINE DLAZQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL,&
     &ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1,                   &
     &DN2, TAU )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            I0, ITER, N0, NDIV, NFAIL, PP, TTYPE
      DOUBLE PRECISION   DESIG, DMIN, DMIN1, DMIN2, DN, DN1, DN2, QMAX, &
     &SIGMA, TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   Z( * )
!     ..
!
!  Purpose
!  =======
!
!  DLAZQ3 checks for deflation, computes a shift (TAU) and calls dqds.
!  In case of failure it changes shifts, and tries again until output
!  is positive.
!
!  Arguments
!  =========
!
!  I0     (input) INTEGER
!         First index.
!
!  N0     (input) INTEGER
!         Last index.
!
!  Z      (input) DOUBLE PRECISION array, dimension ( 4*N )
!         Z holds the qd array.
!
!  PP     (input) INTEGER
!         PP=0 for ping, PP=1 for pong.
!
!  DMIN   (output) DOUBLE PRECISION
!         Minimum value of d.
!
!  SIGMA  (output) DOUBLE PRECISION
!         Sum of shifts used in current segment.
!
!  DESIG  (input/output) DOUBLE PRECISION
!         Lower order part of SIGMA
!
!  QMAX   (input) DOUBLE PRECISION
!         Maximum value of q.
!
!  NFAIL  (output) INTEGER
!         Number of times shift was too big.
!
!  ITER   (output) INTEGER
!         Number of iterations.
!
!  NDIV   (output) INTEGER
!         Number of divisions.
!
!  IEEE   (input) LOGICAL
!         Flag for IEEE or non IEEE arithmetic (passed to DLASQ5).
!
!  TTYPE  (input/output) INTEGER
!         Shift type.  TTYPE is passed as an argument in order to save
!         its value between calls to DLAZQ3
!
!  DMIN1  (input/output) REAL
!  DMIN2  (input/output) REAL
!  DN     (input/output) REAL
!  DN1    (input/output) REAL
!  DN2    (input/output) REAL
!  TAU    (input/output) REAL
!         These are passed as arguments in order to save their values
!         between calls to DLAZQ3
!
!  This is a thread safe version of DLASQ3, which passes TTYPE, DMIN1,
!  DMIN2, DN, DN1. DN2 and TAU through the argument list in place of
!  declaring them in a SAVE statment.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   CBIAS
      PARAMETER          ( CBIAS = 1.50D0 )
      DOUBLE PRECISION   ZERO, QURTR, HALF, ONE, TWO, HUNDRD
      PARAMETER          ( ZERO = 0.0D0, QURTR = 0.250D0, HALF = 0.5D0, &
     &ONE = 1.0D0, TWO = 2.0D0, HUNDRD = 100.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            IPN4, J4, N0IN, NN
      DOUBLE PRECISION   EPS, G, S, SAFMIN, T, TEMP, TOL, TOL2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASQ5, DLASQ6, DLAZQ4
!     ..
!     .. External Function ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      N0IN   = N0
      EPS    = DLAMCH( 'Precision' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      TOL    = EPS*HUNDRD
      TOL2   = TOL**2
      G      = ZERO
!
!     Check for deflation.
!
   10 CONTINUE
!
      IF( N0.LT.I0 )                                                    &
     &RETURN
      IF( N0.EQ.I0 )                                                    &
     &GO TO 20
      NN = 4*N0 + PP
      IF( N0.EQ.( I0+1 ) )                                              &
     &GO TO 40
!
!     Check whether E(N0-1) is negligible, 1 eigenvalue.
!
      IF( Z( NN-5 ).GT.TOL2*( SIGMA+Z( NN-3 ) ) .AND.                   &
     &Z( NN-2*PP-4 ).GT.TOL2*Z( NN-7 ) )                                &
     &GO TO 30
!
   20 CONTINUE
!
      Z( 4*N0-3 ) = Z( 4*N0+PP-3 ) + SIGMA
      N0 = N0 - 1
      GO TO 10
!
!     Check  whether E(N0-2) is negligible, 2 eigenvalues.
!
   30 CONTINUE
!
      IF( Z( NN-9 ).GT.TOL2*SIGMA .AND.                                 &
     &Z( NN-2*PP-8 ).GT.TOL2*Z( NN-11 ) )                               &
     &GO TO 50
!
   40 CONTINUE
!
      IF( Z( NN-3 ).GT.Z( NN-7 ) ) THEN
        S = Z( NN-3 )
        Z( NN-3 ) = Z( NN-7 )
        Z( NN-7 ) = S
      END IF
      IF( Z( NN-5 ).GT.Z( NN-3 )*TOL2 ) THEN
        T = HALF*( ( Z( NN-7 )-Z( NN-3 ) )+Z( NN-5 ) )
        S = Z( NN-3 )*( Z( NN-5 ) / T )
        IF( S.LE.T ) THEN
          S = Z( NN-3 )*( Z( NN-5 ) /                                   &
     &( T*( ONE+SQRT( ONE+S / T ) ) ) )
        ELSE
          S = Z( NN-3 )*( Z( NN-5 ) / ( T+SQRT( T )*SQRT( T+S ) ) )
        END IF
        T = Z( NN-7 ) + ( S+Z( NN-5 ) )
        Z( NN-3 ) = Z( NN-3 )*( Z( NN-7 ) / T )
        Z( NN-7 ) = T
      END IF
      Z( 4*N0-7 ) = Z( NN-7 ) + SIGMA
      Z( 4*N0-3 ) = Z( NN-3 ) + SIGMA
      N0 = N0 - 2
      GO TO 10
!
   50 CONTINUE
!
!     Reverse the qd-array, if warranted.
!
      IF( DMIN.LE.ZERO .OR. N0.LT.N0IN ) THEN
        IF( CBIAS*Z( 4*I0+PP-3 ).LT.Z( 4*N0+PP-3 ) ) THEN
          IPN4 = 4*( I0+N0 )
          DO J4 = 4*I0, 2*( I0+N0-1 ), 4
            TEMP = Z( J4-3 )
            Z( J4-3 ) = Z( IPN4-J4-3 )
            Z( IPN4-J4-3 ) = TEMP
            TEMP = Z( J4-2 )
            Z( J4-2 ) = Z( IPN4-J4-2 )
            Z( IPN4-J4-2 ) = TEMP
            TEMP = Z( J4-1 )
            Z( J4-1 ) = Z( IPN4-J4-5 )
            Z( IPN4-J4-5 ) = TEMP
            TEMP = Z( J4 )
            Z( J4 ) = Z( IPN4-J4-4 )
            Z( IPN4-J4-4 ) = TEMP
          enddo
          IF( N0-I0.LE.4 ) THEN
            Z( 4*N0+PP-1 ) = Z( 4*I0+PP-1 )
            Z( 4*N0-PP ) = Z( 4*I0-PP )
          END IF
          DMIN2 = MIN( DMIN2, Z( 4*N0+PP-1 ) )
          Z( 4*N0+PP-1 ) = MIN( Z( 4*N0+PP-1 ), Z( 4*I0+PP-1 ),         &
     &Z( 4*I0+PP+3 ) )
          Z( 4*N0-PP ) = MIN( Z( 4*N0-PP ), Z( 4*I0-PP ),               &
     &Z( 4*I0-PP+4 ) )
          QMAX = MAX( QMAX, Z( 4*I0+PP-3 ), Z( 4*I0+PP+1 ) )
          DMIN = -ZERO
        END IF
      END IF
!
      IF( DMIN.LT.ZERO .OR. SAFMIN*QMAX.LT.MIN( Z( 4*N0+PP-1 ),         &
     &Z( 4*N0+PP-9 ), DMIN2+Z( 4*N0-PP ) ) ) THEN
!
!        Choose a shift.
!
        CALL DLAZQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, DN1,  &
     &DN2, TAU, TTYPE, G )
!
!        Call dqds until DMIN > 0.
!
   80   CONTINUE
!
        CALL DLASQ5( I0, N0, Z, PP, TAU, DMIN, DMIN1, DMIN2, DN,        &
     &DN1, DN2, IEEE )
!
        NDIV = NDIV + ( N0-I0+2 )
        ITER = ITER + 1
!
!        Check status.
!
        IF( DMIN.GE.ZERO .AND. DMIN1.GT.ZERO ) THEN
!
!           Success.
!
          GO TO 100
!
        ELSE IF( DMIN.LT.ZERO .AND. DMIN1.GT.ZERO .AND.                 &
     &Z( 4*( N0-1 )-PP ).LT.TOL*( SIGMA+DN1 ) .AND.                     &
     &ABS( DN ).LT.TOL*SIGMA ) THEN
!
!           Convergence hidden by negative DN.
!
          Z( 4*( N0-1 )-PP+2 ) = ZERO
          DMIN = ZERO
          GO TO 100
        ELSE IF( DMIN.LT.ZERO ) THEN
!
!           TAU too big. Select new TAU and try again.
!
          NFAIL = NFAIL + 1
          IF( TTYPE.LT.-22 ) THEN
!
!              Failed twice. Play it safe.
!
            TAU = ZERO
          ELSE IF( DMIN1.GT.ZERO ) THEN
!
!              Late failure. Gives excellent shift.
!
            TAU = ( TAU+DMIN )*( ONE-TWO*EPS )
            TTYPE = TTYPE - 11
          ELSE
!
!              Early failure. Divide by 4.
!
            TAU = QURTR*TAU
            TTYPE = TTYPE - 12
          END IF
          GO TO 80
        ELSE IF( DMIN.NE.DMIN ) THEN
!
!           NaN.
!
          TAU = ZERO
          GO TO 80
        ELSE
!
!           Possible underflow. Play it safe.
!
          GO TO 90
        END IF
      END IF
!
!     Risk of underflow.
!
   90 CONTINUE
      CALL DLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, DN1, DN2 )
      NDIV = NDIV + ( N0-I0+2 )
      ITER = ITER + 1
      TAU = ZERO
!
  100 CONTINUE
      IF( TAU.LT.SIGMA ) THEN
        DESIG = DESIG + TAU
        T = SIGMA + DESIG
        DESIG = DESIG - ( T-SIGMA )
      ELSE
        T = SIGMA + TAU
        DESIG = SIGMA - ( T-TAU ) + DESIG
      END IF
      SIGMA = T
!
      RETURN
!
!     End of DLAZQ3
!
      END
      SUBROUTINE DLAZQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN,   &
     &DN1, DN2, TAU, TTYPE, G )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            I0, N0, N0IN, PP, TTYPE
      DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   Z( * )
!     ..
!
!  Purpose
!  =======
!
!  DLAZQ4 computes an approximation TAU to the smallest eigenvalue
!  using values of d from the previous transform.
!
!  I0    (input) INTEGER
!        First index.
!
!  N0    (input) INTEGER
!        Last index.
!
!  Z     (input) DOUBLE PRECISION array, dimension ( 4*N )
!        Z holds the qd array.
!
!  PP    (input) INTEGER
!        PP=0 for ping, PP=1 for pong.
!
!  N0IN  (input) INTEGER
!        The value of N0 at start of EIGTEST.
!
!  DMIN  (input) DOUBLE PRECISION
!        Minimum value of d.
!
!  DMIN1 (input) DOUBLE PRECISION
!        Minimum value of d, excluding D( N0 ).
!
!  DMIN2 (input) DOUBLE PRECISION
!        Minimum value of d, excluding D( N0 ) and D( N0-1 ).
!
!  DN    (input) DOUBLE PRECISION
!        d(N)
!
!  DN1   (input) DOUBLE PRECISION
!        d(N-1)
!
!  DN2   (input) DOUBLE PRECISION
!        d(N-2)
!
!  TAU   (output) DOUBLE PRECISION
!        This is the shift.
!
!  TTYPE (output) INTEGER
!        Shift type.
!
!  G     (input/output) DOUBLE PRECISION
!        G is passed as an argument in order to save its value between
!        calls to DLAZQ4
!
!  Further Details
!  ===============
!  CNST1 = 9/16
!
!  This is a thread safe version of DLASQ4, which passes G through the
!  argument list in place of declaring G in a SAVE statment.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   CNST1, CNST2, CNST3
      PARAMETER          ( CNST1 = 0.5630D0, CNST2 = 1.010D0,           &
     &CNST3 = 1.050D0 )
      DOUBLE PRECISION   QURTR, THIRD, HALF, ZERO, ONE, TWO, HUNDRD
      PARAMETER          ( QURTR = 0.250D0, THIRD = 0.3330D0,           &
     &HALF = 0.50D0, ZERO = 0.0D0, ONE = 1.0D0,                         &
     &TWO = 2.0D0, HUNDRD = 100.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I4, NN, NP
      DOUBLE PRECISION   A2, B1, B2, GAM, GAP1, GAP2, S
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     A negative DMIN forces the shift to take that absolute value
!     TTYPE records the type of shift.
!
      IF( DMIN.LE.ZERO ) THEN
        TAU = -DMIN
        TTYPE = -1
        RETURN
      END IF
!
      NN = 4*N0 + PP
      IF( N0IN.EQ.N0 ) THEN
!
!        No eigenvalues deflated.
!
        IF( DMIN.EQ.DN .OR. DMIN.EQ.DN1 ) THEN
!
          B1 = SQRT( Z( NN-3 ) )*SQRT( Z( NN-5 ) )
          B2 = SQRT( Z( NN-7 ) )*SQRT( Z( NN-9 ) )
          A2 = Z( NN-7 ) + Z( NN-5 )
!
!           Cases 2 and 3.
!
          IF( DMIN.EQ.DN .AND. DMIN1.EQ.DN1 ) THEN
            GAP2 = DMIN2 - A2 - DMIN2*QURTR
            IF( GAP2.GT.ZERO .AND. GAP2.GT.B2 ) THEN
              GAP1 = A2 - DN - ( B2 / GAP2 )*B2
            ELSE
              GAP1 = A2 - DN - ( B1+B2 )
            END IF
            IF( GAP1.GT.ZERO .AND. GAP1.GT.B1 ) THEN
              S = MAX( DN-( B1 / GAP1 )*B1, HALF*DMIN )
              TTYPE = -2
            ELSE
              S = ZERO
              IF( DN.GT.B1 )                                            &
     &S = DN - B1
              IF( A2.GT.( B1+B2 ) )                                     &
     &S = MIN( S, A2-( B1+B2 ) )
              S = MAX( S, THIRD*DMIN )
              TTYPE = -3
            END IF
          ELSE
!
!              Case 4.
!
            TTYPE = -4
            S = QURTR*DMIN
            IF( DMIN.EQ.DN ) THEN
              GAM = DN
              A2 = ZERO
              IF( Z( NN-5 ) .GT. Z( NN-7 ) )                            &
     &RETURN
              B2 = Z( NN-5 ) / Z( NN-7 )
              NP = NN - 9
            ELSE
              NP = NN - 2*PP
              B2 = Z( NP-2 )
              GAM = DN1
              IF( Z( NP-4 ) .GT. Z( NP-2 ) )                            &
     &RETURN
              A2 = Z( NP-4 ) / Z( NP-2 )
              IF( Z( NN-9 ) .GT. Z( NN-11 ) )                           &
     &RETURN
              B2 = Z( NN-9 ) / Z( NN-11 )
              NP = NN - 13
            END IF
!
!              Approximate contribution to norm squared from I < NN-1.
!
            A2 = A2 + B2
            DO I4 = NP, 4*I0 - 1 + PP, -4
              IF( B2.EQ.ZERO )                                          &
     &GO TO 20
              B1 = B2
              IF( Z( I4 ) .GT. Z( I4-2 ) )                              &
     &RETURN
              B2 = B2*( Z( I4 ) / Z( I4-2 ) )
              A2 = A2 + B2
              IF( HUNDRD*MAX( B2, B1 ).LT.A2 .OR. CNST1.LT.A2 )         &
     &GO TO 20
            enddo
   20       CONTINUE
            A2 = CNST3*A2
!
!              Rayleigh quotient residual bound.
!
            IF( A2.LT.CNST1 )                                           &
     &S = GAM*( ONE-SQRT( A2 ) ) / ( ONE+A2 )
          END IF
        ELSE IF( DMIN.EQ.DN2 ) THEN
!
!           Case 5.
!
          TTYPE = -5
          S = QURTR*DMIN
!
!           Compute contribution to norm squared from I > NN-2.
!
          NP = NN - 2*PP
          B1 = Z( NP-2 )
          B2 = Z( NP-6 )
          GAM = DN2
          IF( Z( NP-8 ).GT.B2 .OR. Z( NP-4 ).GT.B1 )                    &
     &RETURN
          A2 = ( Z( NP-8 ) / B2 )*( ONE+Z( NP-4 ) / B1 )
!
!           Approximate contribution to norm squared from I < NN-2.
!
          IF( N0-I0.GT.2 ) THEN
            B2 = Z( NN-13 ) / Z( NN-15 )
            A2 = A2 + B2
            DO I4 = NN - 17, 4*I0 - 1 + PP, -4
              IF( B2.EQ.ZERO )                                          &
     &GO TO 40
              B1 = B2
              IF( Z( I4 ) .GT. Z( I4-2 ) )                              &
     &RETURN
              B2 = B2*( Z( I4 ) / Z( I4-2 ) )
              A2 = A2 + B2
              IF( HUNDRD*MAX( B2, B1 ).LT.A2 .OR. CNST1.LT.A2 )         &
     &GO TO 40
            enddo
   40       CONTINUE
            A2 = CNST3*A2
          END IF
!
          IF( A2.LT.CNST1 )                                             &
     &S = GAM*( ONE-SQRT( A2 ) ) / ( ONE+A2 )
        ELSE
!
!           Case 6, no information to guide us.
!
          IF( TTYPE.EQ.-6 ) THEN
            G = G + THIRD*( ONE-G )
          ELSE IF( TTYPE.EQ.-18 ) THEN
            G = QURTR*THIRD
          ELSE
            G = QURTR
          END IF
          S = G*DMIN
          TTYPE = -6
        END IF
!
      ELSE IF( N0IN.EQ.( N0+1 ) ) THEN
!
!        One eigenvalue just deflated. Use DMIN1, DN1 for DMIN and DN.
!
        IF( DMIN1.EQ.DN1 .AND. DMIN2.EQ.DN2 ) THEN
!
!           Cases 7 and 8.
!
          TTYPE = -7
          S = THIRD*DMIN1
          IF( Z( NN-5 ).GT.Z( NN-7 ) )                                  &
     &RETURN
          B1 = Z( NN-5 ) / Z( NN-7 )
          B2 = B1
          IF( B2.EQ.ZERO )                                              &
     &GO TO 60
          DO I4 = 4*N0 - 9 + PP, 4*I0 - 1 + PP, -4
            A2 = B1
            IF( Z( I4 ).GT.Z( I4-2 ) )                                  &
     &RETURN
            B1 = B1*( Z( I4 ) / Z( I4-2 ) )
            B2 = B2 + B1
            IF( HUNDRD*MAX( B1, A2 ).LT.B2 )                            &
     &GO TO 60
          enddo
   60     CONTINUE
          B2 = SQRT( CNST3*B2 )
          A2 = DMIN1 / ( ONE+B2**2 )
          GAP2 = HALF*DMIN2 - A2
          IF( GAP2.GT.ZERO .AND. GAP2.GT.B2*A2 ) THEN
            S = MAX( S, A2*( ONE-CNST2*A2*( B2 / GAP2 )*B2 ) )
          ELSE
            S = MAX( S, A2*( ONE-CNST2*B2 ) )
            TTYPE = -8
          END IF
        ELSE
!
!           Case 9.
!
          S = QURTR*DMIN1
          IF( DMIN1.EQ.DN1 )                                            &
     &S = HALF*DMIN1
          TTYPE = -9
        END IF
!
      ELSE IF( N0IN.EQ.( N0+2 ) ) THEN
!
!        Two eigenvalues deflated. Use DMIN2, DN2 for DMIN and DN.
!
!        Cases 10 and 11.
!
        IF( DMIN2.EQ.DN2 .AND. TWO*Z( NN-5 ).LT.Z( NN-7 ) ) THEN
          TTYPE = -10
          S = THIRD*DMIN2
          IF( Z( NN-5 ).GT.Z( NN-7 ) )                                  &
     &RETURN
          B1 = Z( NN-5 ) / Z( NN-7 )
          B2 = B1
          IF( B2.EQ.ZERO )                                              &
     &GO TO 80
          DO I4 = 4*N0 - 9 + PP, 4*I0 - 1 + PP, -4
            IF( Z( I4 ).GT.Z( I4-2 ) )                                  &
     &RETURN
            B1 = B1*( Z( I4 ) / Z( I4-2 ) )
            B2 = B2 + B1
            IF( HUNDRD*B1.LT.B2 )                                       &
     &GO TO 80
          enddo
   80     CONTINUE
          B2 = SQRT( CNST3*B2 )
          A2 = DMIN2 / ( ONE+B2**2 )
          GAP2 = Z( NN-7 ) + Z( NN-9 ) -                                &
     &SQRT( Z( NN-11 ) )*SQRT( Z( NN-9 ) ) - A2
          IF( GAP2.GT.ZERO .AND. GAP2.GT.B2*A2 ) THEN
            S = MAX( S, A2*( ONE-CNST2*A2*( B2 / GAP2 )*B2 ) )
          ELSE
            S = MAX( S, A2*( ONE-CNST2*B2 ) )
          END IF
        ELSE
          S = QURTR*DMIN2
          TTYPE = -11
        END IF
      ELSE IF( N0IN.GT.( N0+2 ) ) THEN
!
!        Case 12, more than two eigenvalues deflated. No information.
!
        S = ZERO
        TTYPE = -12
      END IF
!
      TAU = S
      RETURN
!
!     End of DLAZQ4
!
      END
      SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORG2R generates an m by n real matrix Q with orthonormal columns,
!  which is defined as the first n columns of a product of k elementary
!  reflectors of order m
!
!        Q  =  H(1) H(2) . . . H(k)
!
!  as returned by DGEQRF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the i-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQRF in the first k columns of its array
!          argument A.
!          On exit, the m-by-n matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, L
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARF, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M.LT.0 ) THEN
        INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
        INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
        INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DORG2R', -INFO )
        RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 )                                                      &
     &RETURN
!
!     Initialise columns k+1:n to columns of the unit matrix
!
      DO J = K + 1, N
        DO L = 1, M
          A( L, J ) = ZERO
        enddo
        A( J, J ) = ONE
      enddo
!
      DO I = K, 1, -1
!
!        Apply H(i) to A(i:m,i:n) from the left
!
        IF( I.LT.N ) THEN
          A( I, I ) = ONE
          CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),       &
     &A( I, I+1 ), LDA, WORK )
        END IF
        IF( I.LT.M )                                                    &
     &CALL DSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
        A( I, I ) = ONE - TAU( I )
!
!        Set A(1:i-1,i) to zero
!
        DO L = 1, I - 1
          A( L, I ) = ZERO
        enddo
      enddo
      RETURN
!
!     End of DORG2R
!
      END
      SUBROUTINE DORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          VECT
      INTEGER            INFO, K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORGBR generates one of the real orthogonal matrices Q or P**T
!  determined by DGEBRD when reducing a real matrix A to bidiagonal
!  form: A = Q * B * P**T.  Q and P**T are defined as products of
!  elementary reflectors H(i) or G(i) respectively.
!
!  If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q
!  is of order M:
!  if m >= k, Q = H(1) H(2) . . . H(k) and DORGBR returns the first n
!  columns of Q, where m >= n >= k;
!  if m < k, Q = H(1) H(2) . . . H(m-1) and DORGBR returns Q as an
!  M-by-M matrix.
!
!  If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**T
!  is of order N:
!  if k < n, P**T = G(k) . . . G(2) G(1) and DORGBR returns the first m
!  rows of P**T, where n >= m >= k;
!  if k >= n, P**T = G(n-1) . . . G(2) G(1) and DORGBR returns P**T as
!  an N-by-N matrix.
!
!  Arguments
!  =========
!
!  VECT    (input) CHARACTER*1
!          Specifies whether the matrix Q or the matrix P**T is
!          required, as defined in the transformation applied by DGEBRD:
!          = 'Q':  generate Q;
!          = 'P':  generate P**T.
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q or P**T to be returned.
!          M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q or P**T to be returned.
!          N >= 0.
!          If VECT = 'Q', M >= N >= min(M,K);
!          if VECT = 'P', N >= M >= min(N,K).
!
!  K       (input) INTEGER
!          If VECT = 'Q', the number of columns in the original M-by-K
!          matrix reduced by DGEBRD.
!          If VECT = 'P', the number of rows in the original K-by-N
!          matrix reduced by DGEBRD.
!          K >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the vectors which define the elementary reflectors,
!          as returned by DGEBRD.
!          On exit, the M-by-N matrix Q or P**T.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) DOUBLE PRECISION array, dimension
!                                (min(M,K)) if VECT = 'Q'
!                                (min(N,K)) if VECT = 'P'
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i) or G(i), which determines Q or P**T, as
!          returned by DGEBRD in its array argument TAUQ or TAUP.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= max(1,min(M,N)).
!          For optimum performance LWORK >= min(M,N)*NB, where NB
!          is the optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, WANTQ
      INTEGER            I, IINFO, J, LWKOPT, MN, NB
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DORGLQ, DORGQR, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      WANTQ = LSAME( VECT, 'Q' )
      MN = MIN( M, N )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.WANTQ .AND. .NOT.LSAME( VECT, 'P' ) ) THEN
        INFO = -1
      ELSE IF( M.LT.0 ) THEN
        INFO = -2
      ELSE IF( N.LT.0 .OR. ( WANTQ .AND. ( N.GT.M .OR. N.LT.MIN( M,     &
     &K ) ) ) .OR. ( .NOT.WANTQ .AND. ( M.GT.N .OR. M.LT.               &
     &MIN( N, K ) ) ) ) THEN
        INFO = -3
      ELSE IF( K.LT.0 ) THEN
        INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -6
      ELSE IF( LWORK.LT.MAX( 1, MN ) .AND. .NOT.LQUERY ) THEN
        INFO = -9
      END IF
!
      IF( INFO.EQ.0 ) THEN
        IF( WANTQ ) THEN
          NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
        ELSE
          NB = ILAENV( 1, 'DORGLQ', ' ', M, N, K, -1 )
        END IF
        LWKOPT = MAX( 1, MN )*NB
        WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DORGBR', -INFO )
        RETURN
      ELSE IF( LQUERY ) THEN
        RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
        WORK( 1 ) = 1
        RETURN
      END IF
!
      IF( WANTQ ) THEN
!
!        Form Q, determined by a call to DGEBRD to reduce an m-by-k
!        matrix
!
        IF( M.GE.K ) THEN
!
!           If m >= k, assume m >= n >= k
!
          CALL DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, IINFO )
!
        ELSE
!
!           If m < k, assume m = n
!
!           Shift the vectors which define the elementary reflectors one
!           column to the right, and set the first row and column of Q
!           to those of the unit matrix
!
          DO J = M, 2, -1
            A( 1, J ) = ZERO
            DO I = J + 1, M
              A( I, J ) = A( I, J-1 )
            enddo
          enddo
          A( 1, 1 ) = ONE
          DO I = 2, M
            A( I, 1 ) = ZERO
          enddo
          IF( M.GT.1 ) THEN
!
!              Form Q(2:m,2:m)
!
            CALL DORGQR( M-1, M-1, M-1, A( 2, 2 ), LDA, TAU, WORK,      &
     &LWORK, IINFO )
          END IF
        END IF
      ELSE
!
!        Form P', determined by a call to DGEBRD to reduce a k-by-n
!        matrix
!
        IF( K.LT.N ) THEN
!
!           If k < n, assume k <= m <= n
!
          CALL DORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, IINFO )
!
        ELSE
!
!           If k >= n, assume m = n
!
!           Shift the vectors which define the elementary reflectors one
!           row downward, and set the first row and column of P' to
!           those of the unit matrix
!
          A( 1, 1 ) = ONE
          DO I = 2, N
            A( I, 1 ) = ZERO
          enddo
          DO J = 2, N
            DO I = J - 1, 2, -1
              A( I, J ) = A( I-1, J )
            enddo
            A( 1, J ) = ZERO
          enddo
          IF( N.GT.1 ) THEN
!
!              Form P'(2:n,2:n)
!
            CALL DORGLQ( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK,      &
     &LWORK, IINFO )
          END IF
        END IF
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of DORGBR
!
      END
      SUBROUTINE DORGL2( M, N, K, A, LDA, TAU, WORK, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORGL2 generates an m by n real matrix Q with orthonormal rows,
!  which is defined as the first m rows of a product of k elementary
!  reflectors of order n
!
!        Q  =  H(k) . . . H(2) H(1)
!
!  as returned by DGELQF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. N >= M.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. M >= K >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the i-th row must contain the vector which defines
!          the elementary reflector H(i), for i = 1,2,...,k, as returned
!          by DGELQF in the first k rows of its array argument A.
!          On exit, the m-by-n matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGELQF.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (M)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, L
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARF, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M.LT.0 ) THEN
        INFO = -1
      ELSE IF( N.LT.M ) THEN
        INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.M ) THEN
        INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DORGL2', -INFO )
        RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.LE.0 )                                                      &
     &RETURN
!
      IF( K.LT.M ) THEN
!
!        Initialise rows k+1:m to rows of the unit matrix
!
        DO J = 1, N
          DO L = K + 1, M
            A( L, J ) = ZERO
          enddo
          IF( J.GT.K .AND. J.LE.M )                                     &
     &A( J, J ) = ONE
        enddo
      END IF
!
      DO I = K, 1, -1
!
!        Apply H(i) to A(i:m,i:n) from the right
!
        IF( I.LT.N ) THEN
          IF( I.LT.M ) THEN
            A( I, I ) = ONE
            CALL DLARF( 'Right', M-I, N-I+1, A( I, I ), LDA,            &
     &TAU( I ), A( I+1, I ), LDA, WORK )
          END IF
          CALL DSCAL( N-I, -TAU( I ), A( I, I+1 ), LDA )
        END IF
        A( I, I ) = ONE - TAU( I )
!
!        Set A(i,1:i-1) to zero
!
        DO L = 1, I - 1
          A( I, L ) = ZERO
        enddo
      enddo
      RETURN
!
!     End of DORGL2
!
      END
      SUBROUTINE DORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORGLQ generates an M-by-N real matrix Q with orthonormal rows,
!  which is defined as the first M rows of a product of K elementary
!  reflectors of order N
!
!        Q  =  H(k) . . . H(2) H(1)
!
!  as returned by DGELQF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. N >= M.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. M >= K >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the i-th row must contain the vector which defines
!          the elementary reflector H(i), for i = 1,2,...,k, as returned
!          by DGELQF in the first k rows of its array argument A.
!          On exit, the M-by-N matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGELQF.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= max(1,M).
!          For optimum performance LWORK >= M*NB, where NB is
!          the optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK,       &
     &LWKOPT, NB, NBMIN, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORGL2, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      NB = ILAENV( 1, 'DORGLQ', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, M )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
        INFO = -1
      ELSE IF( N.LT.M ) THEN
        INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.M ) THEN
        INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, M ) .AND. .NOT.LQUERY ) THEN
        INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DORGLQ', -INFO )
        RETURN
      ELSE IF( LQUERY ) THEN
        RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.LE.0 ) THEN
        WORK( 1 ) = 1
        RETURN
      END IF
!
      NBMIN = 2
      NX = 0
      IWS = M
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
        NX = MAX( 0, ILAENV( 3, 'DORGLQ', ' ', M, N, K, -1 ) )
        IF( NX.LT.K ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
          LDWORK = M
          IWS = LDWORK*NB
          IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DORGLQ', ' ', M, N, K, -1 ) )
          END IF
        END IF
      END IF
!
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!
!        Use blocked code after the last block.
!        The first kk rows are handled by the block method.
!
        KI = ( ( K-NX-1 ) / NB )*NB
        KK = MIN( K, KI+NB )
!
!        Set A(kk+1:m,1:kk) to zero.
!
        DO J = 1, KK
          DO I = KK + 1, M
            A( I, J ) = ZERO
          enddo
        enddo
      ELSE
        KK = 0
      END IF
!
!     Use unblocked code for the last or only block.
!
      IF( KK.LT.M )                                                     &
     &CALL DORGL2( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA,              &
     &TAU( KK+1 ), WORK, IINFO )
!
      IF( KK.GT.0 ) THEN
!
!        Use blocked code
!
        DO I = KI + 1, 1, -NB
          IB = MIN( NB, K-I+1 )
          IF( I+IB.LE.M ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
            CALL DLARFT( 'Forward', 'Rowwise', N-I+1, IB, A( I, I ),    &
     &LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H' to A(i+ib:m,i:n) from the right
!
            CALL DLARFB( 'Right', 'Transpose', 'Forward', 'Rowwise',    &
     &M-I-IB+1, N-I+1, IB, A( I, I ), LDA, WORK,                        &
     &LDWORK, A( I+IB, I ), LDA, WORK( IB+1 ),                          &
     &LDWORK )
          END IF
!
!           Apply H' to columns i:n of current block
!
          CALL DORGL2( IB, N-I+1, IB, A( I, I ), LDA, TAU( I ), WORK,   &
     &IINFO )
!
!           Set columns 1:i-1 of current block to zero
!
          DO J = 1, I - 1
            DO L = I, I + IB - 1
              A( L, J ) = ZERO
            enddo
          enddo
        enddo
      END IF
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of DORGLQ
!
      END
      SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORGQR generates an M-by-N real matrix Q with orthonormal columns,
!  which is defined as the first N columns of a product of K elementary
!  reflectors of order M
!
!        Q  =  H(1) H(2) . . . H(k)
!
!  as returned by DGEQRF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the i-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQRF in the first k columns of its array
!          argument A.
!          On exit, the M-by-N matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= max(1,N).
!          For optimum performance LWORK >= N*NB, where NB is the
!          optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK,       &
     &LWKOPT, NB, NBMIN, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORG2R, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, N )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
        INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
        INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
        INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
        INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DORGQR', -INFO )
        RETURN
      ELSE IF( LQUERY ) THEN
        RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) THEN
        WORK( 1 ) = 1
        RETURN
      END IF
!
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
        NX = MAX( 0, ILAENV( 3, 'DORGQR', ' ', M, N, K, -1 ) )
        IF( NX.LT.K ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
          LDWORK = N
          IWS = LDWORK*NB
          IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DORGQR', ' ', M, N, K, -1 ) )
          END IF
        END IF
      END IF
!
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!
!        Use blocked code after the last block.
!        The first kk columns are handled by the block method.
!
        KI = ( ( K-NX-1 ) / NB )*NB
        KK = MIN( K, KI+NB )
!
!        Set A(1:kk,kk+1:n) to zero.
!
        DO J = KK + 1, N
          DO I = 1, KK
            A( I, J ) = ZERO
          enddo
        enddo
      ELSE
        KK = 0
      END IF
!
!     Use unblocked code for the last or only block.
!
      IF( KK.LT.N )                                                     &
     &CALL DORG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA,              &
     &TAU( KK+1 ), WORK, IINFO )
!
      IF( KK.GT.0 ) THEN
!
!        Use blocked code
!
        DO I = KI + 1, 1, -NB
          IB = MIN( NB, K-I+1 )
          IF( I+IB.LE.N ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
            CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB,            &
     &A( I, I ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H to A(i:m,i+ib:n) from the left
!
            CALL DLARFB( 'Left', 'No transpose', 'Forward',             &
     &'Columnwise', M-I+1, N-I-IB+1, IB,                                &
     &A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),                       &
     &LDA, WORK( IB+1 ), LDWORK )
          END IF
!
!           Apply H to rows i:m of current block
!
          CALL DORG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK,   &
     &IINFO )
!
!           Set rows 1:i-1 of current block to zero
!
          DO J = I, I + IB - 1
            DO L = 1, I - 1
              A( L, J ) = ZERO
            enddo
          enddo
        enddo
      END IF
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of DORGQR
!
      END
      SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,     &
     &WORK, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORM2R overwrites the general real m by n matrix C with
!
!        Q * C  if SIDE = 'L' and TRANS = 'N', or
!
!        Q'* C  if SIDE = 'L' and TRANS = 'T', or
!
!        C * Q  if SIDE = 'R' and TRANS = 'N', or
!
!        C * Q' if SIDE = 'R' and TRANS = 'T',
!
!  where Q is a real orthogonal matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(1) H(2) . . . H(k)
!
!  as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q' from the Left
!          = 'R': apply Q or Q' from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply Q  (No transpose)
!          = 'T': apply Q' (Transpose)
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
!          The i-th column must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          DGEQRF in the first k columns of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          If SIDE = 'L', LDA >= max(1,M);
!          if SIDE = 'R', LDA >= max(1,N).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension
!                                   (N) if SIDE = 'L',
!                                   (M) if SIDE = 'R'
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      DOUBLE PRECISION   AII
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARF, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
!
!     NQ is the order of Q
!
      IF( LEFT ) THEN
        NQ = M
      ELSE
        NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
        INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
        INFO = -2
      ELSE IF( M.LT.0 ) THEN
        INFO = -3
      ELSE IF( N.LT.0 ) THEN
        INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
        INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
        INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
        INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DORM2R', -INFO )
        RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )                              &
     &RETURN
!
      IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) )  &
     &THEN
        I1 = 1
        I2 = K
        I3 = 1
      ELSE
        I1 = K
        I2 = 1
        I3 = -1
      END IF
!
      IF( LEFT ) THEN
        NI = N
        JC = 1
      ELSE
        MI = M
        IC = 1
      END IF
!
      DO I = I1, I2, I3
        IF( LEFT ) THEN
!
!           H(i) is applied to C(i:m,1:n)
!
          MI = M - I + 1
          IC = I
        ELSE
!
!           H(i) is applied to C(1:m,i:n)
!
          NI = N - I + 1
          JC = I
        END IF
!
!        Apply H(i)
!
        AII = A( I, I )
        A( I, I ) = ONE
        CALL DLARF( SIDE, MI, NI, A( I, I ), 1, TAU( I ), C( IC, JC ),  &
     &LDC, WORK )
        A( I, I ) = AII
      enddo
      RETURN
!
!     End of DORM2R
!
      END
      SUBROUTINE DORMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C,    &
     &LDC, WORK, LWORK, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS, VECT
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  If VECT = 'Q', DORMBR overwrites the general real M-by-N matrix C
!  with
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      Q * C          C * Q
!  TRANS = 'T':      Q**T * C       C * Q**T
!
!  If VECT = 'P', DORMBR overwrites the general real M-by-N matrix C
!  with
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      P * C          C * P
!  TRANS = 'T':      P**T * C       C * P**T
!
!  Here Q and P**T are the orthogonal matrices determined by DGEBRD when
!  reducing a real matrix A to bidiagonal form: A = Q * B * P**T. Q and
!  P**T are defined as products of elementary reflectors H(i) and G(i)
!  respectively.
!
!  Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the
!  order of the orthogonal matrix Q or P**T that is applied.
!
!  If VECT = 'Q', A is assumed to have been an NQ-by-K matrix:
!  if nq >= k, Q = H(1) H(2) . . . H(k);
!  if nq < k, Q = H(1) H(2) . . . H(nq-1).
!
!  If VECT = 'P', A is assumed to have been a K-by-NQ matrix:
!  if k < nq, P = G(1) G(2) . . . G(k);
!  if k >= nq, P = G(1) G(2) . . . G(nq-1).
!
!  Arguments
!  =========
!
!  VECT    (input) CHARACTER*1
!          = 'Q': apply Q or Q**T;
!          = 'P': apply P or P**T.
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q, Q**T, P or P**T from the Left;
!          = 'R': apply Q, Q**T, P or P**T from the Right.
!
!  TRANS   (input) CHARACTER*1
!          = 'N':  No transpose, apply Q  or P;
!          = 'T':  Transpose, apply Q**T or P**T.
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          If VECT = 'Q', the number of columns in the original
!          matrix reduced by DGEBRD.
!          If VECT = 'P', the number of rows in the original
!          matrix reduced by DGEBRD.
!          K >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension
!                                (LDA,min(nq,K)) if VECT = 'Q'
!                                (LDA,nq)        if VECT = 'P'
!          The vectors which define the elementary reflectors H(i) and
!          G(i), whose products determine the matrices Q and P, as
!          returned by DGEBRD.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          If VECT = 'Q', LDA >= max(1,nq);
!          if VECT = 'P', LDA >= max(1,min(nq,K)).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (min(nq,K))
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i) or G(i) which determines Q or P, as returned
!          by DGEBRD in the array argument TAUQ or TAUP.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q
!          or P*C or P**T*C or C*P or C*P**T.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          If SIDE = 'L', LWORK >= max(1,N);
!          if SIDE = 'R', LWORK >= max(1,M).
!          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
!          blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            APPLYQ, LEFT, LQUERY, NOTRAN
      CHARACTER          TRANST
      INTEGER            I1, I2, IINFO, LWKOPT, MI, NB, NI, NQ, NW
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DORMLQ, DORMQR, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      APPLYQ = LSAME( VECT, 'Q' )
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
!
!     NQ is the order of Q or P and NW is the minimum dimension of WORK
!
      IF( LEFT ) THEN
        NQ = M
        NW = N
      ELSE
        NQ = N
        NW = M
      END IF
      IF( .NOT.APPLYQ .AND. .NOT.LSAME( VECT, 'P' ) ) THEN
        INFO = -1
      ELSE IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
        INFO = -2
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
        INFO = -3
      ELSE IF( M.LT.0 ) THEN
        INFO = -4
      ELSE IF( N.LT.0 ) THEN
        INFO = -5
      ELSE IF( K.LT.0 ) THEN
        INFO = -6
      ELSE IF( ( APPLYQ .AND. LDA.LT.MAX( 1, NQ ) ) .OR.                &
     &( .NOT.APPLYQ .AND. LDA.LT.MAX( 1, MIN( NQ, K ) ) ) )             &
     &THEN
        INFO = -8
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
        INFO = -11
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
        INFO = -13
      END IF
!
      IF( INFO.EQ.0 ) THEN
        IF( APPLYQ ) THEN
          IF( LEFT ) THEN
            NB = ILAENV( 1, 'DORMQR', SIDE // TRANS, M-1, N, M-1,       &
     &-1 )
          ELSE
            NB = ILAENV( 1, 'DORMQR', SIDE // TRANS, M, N-1, N-1,       &
     &-1 )
          END IF
        ELSE
          IF( LEFT ) THEN
            NB = ILAENV( 1, 'DORMLQ', SIDE // TRANS, M-1, N, M-1,       &
     &-1 )
          ELSE
            NB = ILAENV( 1, 'DORMLQ', SIDE // TRANS, M, N-1, N-1,       &
     &-1 )
          END IF
        END IF
        LWKOPT = MAX( 1, NW )*NB
        WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DORMBR', -INFO )
        RETURN
      ELSE IF( LQUERY ) THEN
        RETURN
      END IF
!
!     Quick return if possible
!
      WORK( 1 ) = 1
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &RETURN
!
      IF( APPLYQ ) THEN
!
!        Apply Q
!
        IF( NQ.GE.K ) THEN
!
!           Q was determined by a call to DGEBRD with nq >= k
!
          CALL DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,       &
     &WORK, LWORK, IINFO )
        ELSE IF( NQ.GT.1 ) THEN
!
!           Q was determined by a call to DGEBRD with nq < k
!
          IF( LEFT ) THEN
            MI = M - 1
            NI = N
            I1 = 2
            I2 = 1
          ELSE
            MI = M
            NI = N - 1
            I1 = 1
            I2 = 2
          END IF
          CALL DORMQR( SIDE, TRANS, MI, NI, NQ-1, A( 2, 1 ), LDA, TAU,  &
     &C( I1, I2 ), LDC, WORK, LWORK, IINFO )
        END IF
      ELSE
!
!        Apply P
!
        IF( NOTRAN ) THEN
          TRANST = 'T'
        ELSE
          TRANST = 'N'
        END IF
        IF( NQ.GT.K ) THEN
!
!           P was determined by a call to DGEBRD with nq > k
!
          CALL DORMLQ( SIDE, TRANST, M, N, K, A, LDA, TAU, C, LDC,      &
     &WORK, LWORK, IINFO )
        ELSE IF( NQ.GT.1 ) THEN
!
!           P was determined by a call to DGEBRD with nq <= k
!
          IF( LEFT ) THEN
            MI = M - 1
            NI = N
            I1 = 2
            I2 = 1
          ELSE
            MI = M
            NI = N - 1
            I1 = 1
            I2 = 2
          END IF
          CALL DORMLQ( SIDE, TRANST, MI, NI, NQ-1, A( 1, 2 ), LDA,      &
     &TAU, C( I1, I2 ), LDC, WORK, LWORK, IINFO )
        END IF
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of DORMBR
!
      END
      SUBROUTINE DORML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,     &
     &WORK, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORML2 overwrites the general real m by n matrix C with
!
!        Q * C  if SIDE = 'L' and TRANS = 'N', or
!
!        Q'* C  if SIDE = 'L' and TRANS = 'T', or
!
!        C * Q  if SIDE = 'R' and TRANS = 'N', or
!
!        C * Q' if SIDE = 'R' and TRANS = 'T',
!
!  where Q is a real orthogonal matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(k) . . . H(2) H(1)
!
!  as returned by DGELQF. Q is of order m if SIDE = 'L' and of order n
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q' from the Left
!          = 'R': apply Q or Q' from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply Q  (No transpose)
!          = 'T': apply Q' (Transpose)
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension
!                               (LDA,M) if SIDE = 'L',
!                               (LDA,N) if SIDE = 'R'
!          The i-th row must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          DGELQF in the first k rows of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1,K).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGELQF.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension
!                                   (N) if SIDE = 'L',
!                                   (M) if SIDE = 'R'
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      DOUBLE PRECISION   AII
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARF, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
!
!     NQ is the order of Q
!
      IF( LEFT ) THEN
        NQ = M
      ELSE
        NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
        INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
        INFO = -2
      ELSE IF( M.LT.0 ) THEN
        INFO = -3
      ELSE IF( N.LT.0 ) THEN
        INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
        INFO = -5
      ELSE IF( LDA.LT.MAX( 1, K ) ) THEN
        INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
        INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DORML2', -INFO )
        RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )                              &
     &RETURN
!
      IF( ( LEFT .AND. NOTRAN ) .OR. ( .NOT.LEFT .AND. .NOT.NOTRAN ) )  &
     &THEN
        I1 = 1
        I2 = K
        I3 = 1
      ELSE
        I1 = K
        I2 = 1
        I3 = -1
      END IF
!
      IF( LEFT ) THEN
        NI = N
        JC = 1
      ELSE
        MI = M
        IC = 1
      END IF
!
      DO I = I1, I2, I3
        IF( LEFT ) THEN
!
!           H(i) is applied to C(i:m,1:n)
!
          MI = M - I + 1
          IC = I
        ELSE
!
!           H(i) is applied to C(1:m,i:n)
!
          NI = N - I + 1
          JC = I
        END IF
!
!        Apply H(i)
!
        AII = A( I, I )
        A( I, I ) = ONE
        CALL DLARF( SIDE, MI, NI, A( I, I ), LDA, TAU( I ),             &
     &C( IC, JC ), LDC, WORK )
        A( I, I ) = AII
      enddo
      RETURN
!
!     End of DORML2
!
      END
      SUBROUTINE DORMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,     &
     &WORK, LWORK, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORMLQ overwrites the general real M-by-N matrix C with
!
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      Q * C          C * Q
!  TRANS = 'T':      Q**T * C       C * Q**T
!
!  where Q is a real orthogonal matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(k) . . . H(2) H(1)
!
!  as returned by DGELQF. Q is of order M if SIDE = 'L' and of order N
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q**T from the Left;
!          = 'R': apply Q or Q**T from the Right.
!
!  TRANS   (input) CHARACTER*1
!          = 'N':  No transpose, apply Q;
!          = 'T':  Transpose, apply Q**T.
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension
!                               (LDA,M) if SIDE = 'L',
!                               (LDA,N) if SIDE = 'R'
!          The i-th row must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          DGELQF in the first k rows of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1,K).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGELQF.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          If SIDE = 'L', LWORK >= max(1,N);
!          if SIDE = 'R', LWORK >= max(1,M).
!          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
!          blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      CHARACTER          TRANST
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWS, JC, LDWORK, &
     &LWKOPT, MI, NB, NBMIN, NI, NQ, NW
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   T( LDT, NBMAX )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORML2, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF( LEFT ) THEN
        NQ = M
        NW = N
      ELSE
        NQ = N
        NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
        INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
        INFO = -2
      ELSE IF( M.LT.0 ) THEN
        INFO = -3
      ELSE IF( N.LT.0 ) THEN
        INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
        INFO = -5
      ELSE IF( LDA.LT.MAX( 1, K ) ) THEN
        INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
        INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
        INFO = -12
      END IF
!
      IF( INFO.EQ.0 ) THEN
!
!        Determine the block size.  NB may be at most NBMAX, where NBMAX
!        is used to define the local array T.
!
        NB = MIN( NBMAX, ILAENV( 1, 'DORMLQ', SIDE // TRANS, M, N, K,   &
     &-1 ) )
        LWKOPT = MAX( 1, NW )*NB
        WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DORMLQ', -INFO )
        RETURN
      ELSE IF( LQUERY ) THEN
        RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
        WORK( 1 ) = 1
        RETURN
      END IF
!
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
        IWS = NW*NB
        IF( LWORK.LT.IWS ) THEN
          NB = LWORK / LDWORK
          NBMIN = MAX( 2, ILAENV( 2, 'DORMLQ', SIDE // TRANS, M, N, K,  &
     &-1 ) )
        END IF
      ELSE
        IWS = NW
      END IF
!
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
!
!        Use unblocked code
!
        CALL DORML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK,   &
     &IINFO )
      ELSE
!
!        Use blocked code
!
        IF( ( LEFT .AND. NOTRAN ) .OR.                                  &
     &( .NOT.LEFT .AND. .NOT.NOTRAN ) ) THEN
          I1 = 1
          I2 = K
          I3 = NB
        ELSE
          I1 = ( ( K-1 ) / NB )*NB + 1
          I2 = 1
          I3 = -NB
        END IF
!
        IF( LEFT ) THEN
          NI = N
          JC = 1
        ELSE
          MI = M
          IC = 1
        END IF
!
        IF( NOTRAN ) THEN
          TRANST = 'T'
        ELSE
          TRANST = 'N'
        END IF
!
        DO I = I1, I2, I3
          IB = MIN( NB, K-I+1 )
!
!           Form the triangular factor of the block reflector
!           H = H(i) H(i+1) . . . H(i+ib-1)
!
          CALL DLARFT( 'Forward', 'Rowwise', NQ-I+1, IB, A( I, I ),     &
     &LDA, TAU( I ), T, LDT )
          IF( LEFT ) THEN
!
!              H or H' is applied to C(i:m,1:n)
!
            MI = M - I + 1
            IC = I
          ELSE
!
!              H or H' is applied to C(1:m,i:n)
!
            NI = N - I + 1
            JC = I
          END IF
!
!           Apply H or H'
!
          CALL DLARFB( SIDE, TRANST, 'Forward', 'Rowwise', MI, NI, IB,  &
     &A( I, I ), LDA, T, LDT, C( IC, JC ), LDC, WORK,                   &
     &LDWORK )
        enddo
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of DORMLQ
!
      END
      SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,     &
     &WORK, LWORK, INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORMQR overwrites the general real M-by-N matrix C with
!
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      Q * C          C * Q
!  TRANS = 'T':      Q**T * C       C * Q**T
!
!  where Q is a real orthogonal matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(1) H(2) . . . H(k)
!
!  as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q**T from the Left;
!          = 'R': apply Q or Q**T from the Right.
!
!  TRANS   (input) CHARACTER*1
!          = 'N':  No transpose, apply Q;
!          = 'T':  Transpose, apply Q**T.
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
!          The i-th column must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          DGEQRF in the first k columns of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          If SIDE = 'L', LDA >= max(1,M);
!          if SIDE = 'R', LDA >= max(1,N).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          If SIDE = 'L', LWORK >= max(1,N);
!          if SIDE = 'R', LWORK >= max(1,M).
!          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
!          blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWS, JC, LDWORK, &
     &LWKOPT, MI, NB, NBMIN, NI, NQ, NW
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   T( LDT, NBMAX )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORM2R, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF( LEFT ) THEN
        NQ = M
        NW = N
      ELSE
        NQ = N
        NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
        INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
        INFO = -2
      ELSE IF( M.LT.0 ) THEN
        INFO = -3
      ELSE IF( N.LT.0 ) THEN
        INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
        INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
        INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
        INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
        INFO = -12
      END IF
!
      IF( INFO.EQ.0 ) THEN
!
!        Determine the block size.  NB may be at most NBMAX, where NBMAX
!        is used to define the local array T.
!
        NB = MIN( NBMAX, ILAENV( 1, 'DORMQR', SIDE // TRANS, M, N, K,   &
     &-1 ) )
        LWKOPT = MAX( 1, NW )*NB
        WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DORMQR', -INFO )
        RETURN
      ELSE IF( LQUERY ) THEN
        RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
        WORK( 1 ) = 1
        RETURN
      END IF
!
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
        IWS = NW*NB
        IF( LWORK.LT.IWS ) THEN
          NB = LWORK / LDWORK
          NBMIN = MAX( 2, ILAENV( 2, 'DORMQR', SIDE // TRANS, M, N, K,  &
     &-1 ) )
        END IF
      ELSE
        IWS = NW
      END IF
!
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
!
!        Use unblocked code
!
        CALL DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK,   &
     &IINFO )
      ELSE
!
!        Use blocked code
!
        IF( ( LEFT .AND. .NOT.NOTRAN ) .OR.                             &
     &( .NOT.LEFT .AND. NOTRAN ) ) THEN
          I1 = 1
          I2 = K
          I3 = NB
        ELSE
          I1 = ( ( K-1 ) / NB )*NB + 1
          I2 = 1
          I3 = -NB
        END IF
!
        IF( LEFT ) THEN
          NI = N
          JC = 1
        ELSE
          MI = M
          IC = 1
        END IF
!
        DO I = I1, I2, I3
          IB = MIN( NB, K-I+1 )
!
!           Form the triangular factor of the block reflector
!           H = H(i) H(i+1) . . . H(i+ib-1)
!
          CALL DLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I, I ),  &
     &LDA, TAU( I ), T, LDT )
          IF( LEFT ) THEN
!
!              H or H' is applied to C(i:m,1:n)
!
            MI = M - I + 1
            IC = I
          ELSE
!
!              H or H' is applied to C(1:m,i:n)
!
            NI = N - I + 1
            JC = I
          END IF
!
!           Apply H or H'
!
          CALL DLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI,    &
     &IB, A( I, I ), LDA, T, LDT, C( IC, JC ), LDC,                     &
     &WORK, LDWORK )
        enddo
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of DORMQR
!
      END
      SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
      implicit none
!     .. Scalar Arguments ..
      DOUBLE PRECISION C,S
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  Purpose
!  =======
!
!     applies a plane rotation.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY
!     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO I = 1,N
        DTEMP = C*DX(IX) + S*DY(IY)
        DY(IY) = C*DY(IY) - S*DX(IX)
        DX(IX) = DTEMP
        IX = IX + INCX
        IY = IY + INCY
      enddo
      RETURN
!
!       code for both increments equal to 1
!
   20 DO I = 1,N
        DTEMP = C*DX(I) + S*DY(I)
        DY(I) = C*DY(I) - S*DX(I)
        DX(I) = DTEMP
      enddo
      RETURN
      END
      SUBROUTINE DSCAL(N,DA,DX,INCX)
      implicit none
!     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
!     ..
!
!  Purpose
!  =======
!*
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      INTEGER I,M,MP1,NINCX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
!
!        code for increment not equal to 1
!
      NINCX = N*INCX
      DO I = 1,NINCX,INCX
        DX(I) = DA*DX(I)
      enddo
      RETURN
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO I = 1,M
        DX(I) = DA*DX(I)
      enddo
      IF (N.LT.5) RETURN
   40 MP1 = M + 1
      DO I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I+1) = DA*DX(I+1)
        DX(I+2) = DA*DX(I+2)
        DX(I+3) = DA*DX(I+3)
        DX(I+4) = DA*DX(I+4)
      enddo
      RETURN
      END
      SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
      implicit none
!     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  Purpose
!  =======
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
!     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO I = 1,N
        DTEMP = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP
        IX = IX + INCX
        IY = IY + INCY
      enddo
      RETURN
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
   20 M = MOD(N,3)
      IF (M.EQ.0) GO TO 40
      DO I = 1,M
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
      enddo
      IF (N.LT.3) RETURN
   40 MP1 = M + 1
      DO I = MP1,N,3
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
        DTEMP = DX(I+1)
        DX(I+1) = DY(I+1)
        DY(I+1) = DTEMP
        DTEMP = DX(I+2)
        DX(I+2) = DY(I+2)
        DY(I+2) = DTEMP
      enddo
      RETURN
      END
      SUBROUTINE DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      implicit none
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*)
!     ..
!
!  Purpose
!  =======
!
!  DTRMM  performs one of the matrix-matrix operations
!
!     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
!
!  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  Arguments
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry,  SIDE specifies whether  op( A ) multiplies B from
!           the left or right as follows:
!
!              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
!
!              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = A'.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain the matrix  B,  and  on exit  is overwritten  by the
!           transformed matrix.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOUNIT,UPPER
!     ..
!     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!
!     Test the input parameters.
!
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
        NROWA = M
      ELSE
        NROWA = N
      END IF
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
!
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
        INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
        INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.                           &
     &(.NOT.LSAME(TRANSA,'T')) .AND.                                    &
     &(.NOT.LSAME(TRANSA,'C'))) THEN
        INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
        INFO = 4
      ELSE IF (M.LT.0) THEN
        INFO = 5
      ELSE IF (N.LT.0) THEN
        INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
        INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
        INFO = 11
      END IF
      IF (INFO.NE.0) THEN
        CALL XERBLA('DTRMM ',INFO)
        RETURN
      END IF
!
!     Quick return if possible.
!
      IF (N.EQ.0) RETURN
!
!     And when  alpha.eq.zero.
!
      IF (ALPHA.EQ.ZERO) THEN
        DO J = 1,N
          DO I = 1,M
            B(I,J) = ZERO
          enddo
        enddo
        RETURN
      END IF
!
!     Start the operations.
!
      IF (LSIDE) THEN
        IF (LSAME(TRANSA,'N')) THEN
!
!           Form  B := alpha*A*B.
!
          IF (UPPER) THEN
            DO J = 1,N
              DO K = 1,M
                IF (B(K,J).NE.ZERO) THEN
                  TEMP = ALPHA*B(K,J)
                  DO I = 1,K - 1
                    B(I,J) = B(I,J) + TEMP*A(I,K)
                  enddo
                  IF (NOUNIT) TEMP = TEMP*A(K,K)
                  B(K,J) = TEMP
                END IF
              enddo
            enddo
          ELSE
            DO J = 1,N
              DO K = M,1,-1
                IF (B(K,J).NE.ZERO) THEN
                  TEMP = ALPHA*B(K,J)
                  B(K,J) = TEMP
                  IF (NOUNIT) B(K,J) = B(K,J)*A(K,K)
                  DO I = K + 1,M
                    B(I,J) = B(I,J) + TEMP*A(I,K)
                  enddo
                END IF
              enddo
            enddo
          END IF
        ELSE
!
!           Form  B := alpha*A'*B.
!
          IF (UPPER) THEN
            DO J = 1,N
              DO I = M,1,-1
                TEMP = B(I,J)
                IF (NOUNIT) TEMP = TEMP*A(I,I)
                DO K = 1,I - 1
                  TEMP = TEMP + A(K,I)*B(K,J)
                enddo
                B(I,J) = ALPHA*TEMP
              enddo
            enddo
          ELSE
            DO J = 1,N
              DO I = 1,M
                TEMP = B(I,J)
                IF (NOUNIT) TEMP = TEMP*A(I,I)
                DO K = I + 1,M
                  TEMP = TEMP + A(K,I)*B(K,J)
                enddo
                B(I,J) = ALPHA*TEMP
              enddo
            enddo
          END IF
        END IF
      ELSE
        IF (LSAME(TRANSA,'N')) THEN
!
!           Form  B := alpha*B*A.
!
          IF (UPPER) THEN
            DO J = N,1,-1
              TEMP = ALPHA
              IF (NOUNIT) TEMP = TEMP*A(J,J)
              DO I = 1,M
                B(I,J) = TEMP*B(I,J)
              enddo
              DO K = 1,J - 1
                IF (A(K,J).NE.ZERO) THEN
                  TEMP = ALPHA*A(K,J)
                  DO I = 1,M
                    B(I,J) = B(I,J) + TEMP*B(I,K)
                  enddo
                END IF
              enddo
            enddo
          ELSE
            DO J = 1,N
              TEMP = ALPHA
              IF (NOUNIT) TEMP = TEMP*A(J,J)
              DO I = 1,M
                B(I,J) = TEMP*B(I,J)
              enddo
              DO K = J + 1,N
                IF (A(K,J).NE.ZERO) THEN
                  TEMP = ALPHA*A(K,J)
                  DO I = 1,M
                    B(I,J) = B(I,J) + TEMP*B(I,K)
                  enddo
                END IF
              enddo
            enddo
          END IF
        ELSE
!
!           Form  B := alpha*B*A'.
!
          IF (UPPER) THEN
            DO K = 1,N
              DO J = 1,K - 1
                IF (A(J,K).NE.ZERO) THEN
                  TEMP = ALPHA*A(J,K)
                  DO I = 1,M
                    B(I,J) = B(I,J) + TEMP*B(I,K)
                  enddo
                END IF
              enddo
              TEMP = ALPHA
              IF (NOUNIT) TEMP = TEMP*A(K,K)
              IF (TEMP.NE.ONE) THEN
                DO I = 1,M
                  B(I,K) = TEMP*B(I,K)
                enddo
              END IF
            enddo
          ELSE
            DO K = N,1,-1
              DO J = K + 1,N
                IF (A(J,K).NE.ZERO) THEN
                  TEMP = ALPHA*A(J,K)
                  DO I = 1,M
                    B(I,J) = B(I,J) + TEMP*B(I,K)
                  enddo
                END IF
              enddo
              TEMP = ALPHA
              IF (NOUNIT) TEMP = TEMP*A(K,K)
              IF (TEMP.NE.ONE) THEN
                DO I = 1,M
                  B(I,K) = TEMP*B(I,K)
                enddo
              END IF
            enddo
          END IF
        END IF
      END IF
!
      RETURN
!
!     End of DTRMM .
!
      END
      SUBROUTINE DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
      implicit none
!     .. Scalar Arguments ..
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*)
!     ..
!
!  Purpose
!  =======
!
!  DTRMV  performs one of the matrix-vector operations
!
!     x := A*x,   or   x := A'*x,
!
!  where x is an n element vector and  A is an n by n unit, or non-unit,
!  upper or lower triangular matrix.
!
!  Arguments
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   x := A*x.
!
!              TRANS = 'T' or 't'   x := A'*x.
!
!              TRANS = 'C' or 'c'   x := A'*x.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular matrix and the strictly lower triangular part of
!           A is not referenced.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular matrix and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!           A are not referenced either, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x. On exit, X is overwritten with the
!           tranformed vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOUNIT
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
        INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.  &
     &.NOT.LSAME(TRANS,'C')) THEN
        INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
        INFO = 3
      ELSE IF (N.LT.0) THEN
        INFO = 4
      ELSE IF (LDA.LT.MAX(1,N)) THEN
        INFO = 6
      ELSE IF (INCX.EQ.0) THEN
        INFO = 8
      END IF
      IF (INFO.NE.0) THEN
        CALL XERBLA('DTRMV ',INFO)
        RETURN
      END IF
!
!     Quick return if possible.
!
      IF (N.EQ.0) RETURN
!
      NOUNIT = LSAME(DIAG,'N')
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
      IF (INCX.LE.0) THEN
        KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
        KX = 1
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF (LSAME(TRANS,'N')) THEN
!
!        Form  x := A*x.
!
        IF (LSAME(UPLO,'U')) THEN
          IF (INCX.EQ.1) THEN
            DO J = 1,N
              IF (X(J).NE.ZERO) THEN
                TEMP = X(J)
                DO I = 1,J - 1
                  X(I) = X(I) + TEMP*A(I,J)
                enddo
                IF (NOUNIT) X(J) = X(J)*A(J,J)
              END IF
            enddo
          ELSE
            JX = KX
            DO J = 1,N
              IF (X(JX).NE.ZERO) THEN
                TEMP = X(JX)
                IX = KX
                DO I = 1,J - 1
                  X(IX) = X(IX) + TEMP*A(I,J)
                  IX = IX + INCX
                enddo
                IF (NOUNIT) X(JX) = X(JX)*A(J,J)
              END IF
              JX = JX + INCX
            enddo
          END IF
        ELSE
          IF (INCX.EQ.1) THEN
            DO J = N,1,-1
              IF (X(J).NE.ZERO) THEN
                TEMP = X(J)
                DO I = N,J + 1,-1
                  X(I) = X(I) + TEMP*A(I,J)
                enddo
                IF (NOUNIT) X(J) = X(J)*A(J,J)
              END IF
            enddo
          ELSE
            KX = KX + (N-1)*INCX
            JX = KX
            DO J = N,1,-1
              IF (X(JX).NE.ZERO) THEN
                TEMP = X(JX)
                IX = KX
                DO I = N,J + 1,-1
                  X(IX) = X(IX) + TEMP*A(I,J)
                  IX = IX - INCX
                enddo
                IF (NOUNIT) X(JX) = X(JX)*A(J,J)
              END IF
              JX = JX - INCX
            enddo
          END IF
        END IF
      ELSE
!
!        Form  x := A'*x.
!
        IF (LSAME(UPLO,'U')) THEN
          IF (INCX.EQ.1) THEN
            DO J = N,1,-1
              TEMP = X(J)
              IF (NOUNIT) TEMP = TEMP*A(J,J)
              DO I = J - 1,1,-1
                TEMP = TEMP + A(I,J)*X(I)
              enddo
              X(J) = TEMP
            enddo
          ELSE
            JX = KX + (N-1)*INCX
            DO J = N,1,-1
              TEMP = X(JX)
              IX = JX
              IF (NOUNIT) TEMP = TEMP*A(J,J)
              DO I = J - 1,1,-1
                IX = IX - INCX
                TEMP = TEMP + A(I,J)*X(IX)
              enddo
              X(JX) = TEMP
              JX = JX - INCX
            enddo
          END IF
        ELSE
          IF (INCX.EQ.1) THEN
            DO J = 1,N
              TEMP = X(J)
              IF (NOUNIT) TEMP = TEMP*A(J,J)
              DO I = J + 1,N
                TEMP = TEMP + A(I,J)*X(I)
              enddo
              X(J) = TEMP
            enddo
          ELSE
            JX = KX
            DO J = 1,N
              TEMP = X(JX)
              IX = JX
              IF (NOUNIT) TEMP = TEMP*A(J,J)
              DO I = J + 1,N
                IX = IX + INCX
                TEMP = TEMP + A(I,J)*X(IX)
              enddo
              X(JX) = TEMP
              JX = JX + INCX
            enddo
          END IF
        END IF
      END IF
!
      RETURN
!
!     End of DTRMV .
!
      END
      SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      implicit none
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*)
!     ..
!
!  Purpose
!  =======
!
!  DTRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  The matrix X is overwritten on B.
!
!  Arguments
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = A'.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOUNIT,UPPER
!     ..
!     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!
!     Test the input parameters.
!
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
        NROWA = M
      ELSE
        NROWA = N
      END IF
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
!
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
        INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
        INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.                           &
     &(.NOT.LSAME(TRANSA,'T')) .AND.                                    &
     &(.NOT.LSAME(TRANSA,'C'))) THEN
        INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
        INFO = 4
      ELSE IF (M.LT.0) THEN
        INFO = 5
      ELSE IF (N.LT.0) THEN
        INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
        INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
        INFO = 11
      END IF
      IF (INFO.NE.0) THEN
        CALL XERBLA('DTRSM ',INFO)
        RETURN
      END IF
!
!     Quick return if possible.
!
      IF (N.EQ.0) RETURN
!
!     And when  alpha.eq.zero.
!
      IF (ALPHA.EQ.ZERO) THEN
        DO J = 1,N
          DO I = 1,M
            B(I,J) = ZERO
          enddo
        enddo
        RETURN
      END IF
!
!     Start the operations.
!
      IF (LSIDE) THEN
        IF (LSAME(TRANSA,'N')) THEN
!
!           Form  B := alpha*inv( A )*B.
!
          IF (UPPER) THEN
            DO J = 1,N
              IF (ALPHA.NE.ONE) THEN
                DO I = 1,M
                  B(I,J) = ALPHA*B(I,J)
                enddo
              END IF
              DO K = M,1,-1
                IF (B(K,J).NE.ZERO) THEN
                  IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                  DO I = 1,K - 1
                    B(I,J) = B(I,J) - B(K,J)*A(I,K)
                  enddo
                END IF
              enddo
            enddo
          ELSE
            DO J = 1,N
              IF (ALPHA.NE.ONE) THEN
                DO I = 1,M
                  B(I,J) = ALPHA*B(I,J)
                enddo
              END IF
              DO K = 1,M
                IF (B(K,J).NE.ZERO) THEN
                  IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                  DO I = K + 1,M
                    B(I,J) = B(I,J) - B(K,J)*A(I,K)
                  enddo
                END IF
              enddo
            enddo
          END IF
        ELSE
!
!           Form  B := alpha*inv( A' )*B.
!
          IF (UPPER) THEN
            DO J = 1,N
              DO I = 1,M
                TEMP = ALPHA*B(I,J)
                DO K = 1,I - 1
                  TEMP = TEMP - A(K,I)*B(K,J)
                enddo
                IF (NOUNIT) TEMP = TEMP/A(I,I)
                B(I,J) = TEMP
              enddo
            enddo
          ELSE
            DO J = 1,N
              DO I = M,1,-1
                TEMP = ALPHA*B(I,J)
                DO K = I + 1,M
                  TEMP = TEMP - A(K,I)*B(K,J)
                enddo
                IF (NOUNIT) TEMP = TEMP/A(I,I)
                B(I,J) = TEMP
              enddo
            enddo
          END IF
        END IF
      ELSE
        IF (LSAME(TRANSA,'N')) THEN
!
!           Form  B := alpha*B*inv( A ).
!
          IF (UPPER) THEN
            DO J = 1,N
              IF (ALPHA.NE.ONE) THEN
                DO I = 1,M
                  B(I,J) = ALPHA*B(I,J)
                enddo
              END IF
              DO K = 1,J - 1
                IF (A(K,J).NE.ZERO) THEN
                  DO I = 1,M
                    B(I,J) = B(I,J) - A(K,J)*B(I,K)
                  enddo
                END IF
              enddo
              IF (NOUNIT) THEN
                TEMP = ONE/A(J,J)
                DO I = 1,M
                  B(I,J) = TEMP*B(I,J)
                enddo
              END IF
            enddo
          ELSE
            DO J = N,1,-1
              IF (ALPHA.NE.ONE) THEN
                DO I = 1,M
                  B(I,J) = ALPHA*B(I,J)
                enddo
              END IF
              DO K = J + 1,N
                IF (A(K,J).NE.ZERO) THEN
                  DO I = 1,M
                    B(I,J) = B(I,J) - A(K,J)*B(I,K)
                  enddo
                END IF
              enddo
              IF (NOUNIT) THEN
                TEMP = ONE/A(J,J)
                DO I = 1,M
                  B(I,J) = TEMP*B(I,J)
                enddo
              END IF
            enddo
          END IF
        ELSE
!
!           Form  B := alpha*B*inv( A' ).
!
          IF (UPPER) THEN
            DO K = N,1,-1
              IF (NOUNIT) THEN
                TEMP = ONE/A(K,K)
                DO I = 1,M
                  B(I,K) = TEMP*B(I,K)
                enddo
              END IF
              DO J = 1,K - 1
                IF (A(J,K).NE.ZERO) THEN
                  TEMP = A(J,K)
                  DO I = 1,M
                    B(I,J) = B(I,J) - TEMP*B(I,K)
                  enddo
                END IF
              enddo
              IF (ALPHA.NE.ONE) THEN
                DO I = 1,M
                  B(I,K) = ALPHA*B(I,K)
                enddo
              END IF
            enddo
          ELSE
            DO K = 1,N
              IF (NOUNIT) THEN
                TEMP = ONE/A(K,K)
                DO I = 1,M
                  B(I,K) = TEMP*B(I,K)
                enddo
              END IF
              DO J = K + 1,N
                IF (A(J,K).NE.ZERO) THEN
                  TEMP = A(J,K)
                  DO I = 1,M
                    B(I,J) = B(I,J) - TEMP*B(I,K)
                  enddo
                END IF
              enddo
              IF (ALPHA.NE.ONE) THEN
                DO I = 1,M
                  B(I,K) = ALPHA*B(I,K)
                enddo
              END IF
            enddo
          END IF
        END IF
      END IF
!
      RETURN
!
!     End of DTRSM .
!
      END
      SUBROUTINE DTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB,    &
     &INFO )
      implicit none
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DTRTRS solves a triangular system of the form
!
!     A * X = B  or  A**T * X = B,
!
!  where A is a triangular matrix of order N, and B is an N-by-NRHS
!  matrix.  A check is made to verify that A is nonsingular.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  A is upper triangular;
!          = 'L':  A is lower triangular.
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A**T * X = B  (Transpose)
!          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
!
!  DIAG    (input) CHARACTER*1
!          = 'N':  A is non-unit triangular;
!          = 'U':  A is unit triangular.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The triangular matrix A.  If UPLO = 'U', the leading N-by-N
!          upper triangular part of the array A contains the upper
!          triangular matrix, and the strictly lower triangular part of
!          A is not referenced.  If UPLO = 'L', the leading N-by-N lower
!          triangular part of the array A contains the lower triangular
!          matrix, and the strictly upper triangular part of A is not
!          referenced.  If DIAG = 'U', the diagonal elements of A are
!          also not referenced and are assumed to be 1.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, if INFO = 0, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = i, the i-th diagonal element of A is zero,
!               indicating that the matrix is singular and the solutions
!               X have not been computed.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOUNIT
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
        INFO = -1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.                     &
     &LSAME( TRANS, 'T' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
        INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
        INFO = -3
      ELSE IF( N.LT.0 ) THEN
        INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
        INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
        INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
        INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DTRTRS', -INFO )
        RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 )                                                      &
     &RETURN
!
!     Check for singularity.
!
      IF( NOUNIT ) THEN
        DO INFO = 1, N
          IF( A( INFO, INFO ).EQ.ZERO )                                 &
     &RETURN
        enddo
      END IF
      INFO = 0
!
!     Solve A * x = b  or  A' * x = b.
!
      CALL DTRSM( 'Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B,   &
     &LDB )
!
      RETURN
!
!     End of DTRTRS
!
      END
      SUBROUTINE XERBLA( SRNAME, INFO )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
!     ..
!
!  Purpose
!  =======
!
!  XERBLA  is an error handler for the LAPACK routines.
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!  =========
!
!  SRNAME  (input) CHARACTER*6
!          The name of the routine which called XERBLA.
!
!  INFO    (input) INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
! =====================================================================
!
!     .. Executable Statements ..
!
      WRITE( *, FMT = 9999 )SRNAME, INFO
!
      STOP
!
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',&
     &'an illegal value' )
!
!     End of XERBLA
!
      END
      SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
      implicit none
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASSQ  returns the values  scl  and  smsq  such that
!
!     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
!
!  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
!  assumed to be non-negative and  scl  returns the value
!
!     scl = max( scale, abs( x( i ) ) ).
!
!  scale and sumsq must be supplied in SCALE and SUMSQ and
!  scl and smsq are overwritten on SCALE and SUMSQ respectively.
!
!  The routine makes only one pass through the vector x.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements to be used from the vector X.
!
!  X       (input) DOUBLE PRECISION array, dimension (N)
!          The vector for which a scaled sum of squares is computed.
!             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
!
!  INCX    (input) INTEGER
!          The increment between successive values of the vector X.
!          INCX > 0.
!
!  SCALE   (input/output) DOUBLE PRECISION
!          On entry, the value  scale  in the equation above.
!          On exit, SCALE is overwritten with  scl , the scaling factor
!          for the sum of squares.
!
!  SUMSQ   (input/output) DOUBLE PRECISION
!          On entry, the value  sumsq  in the equation above.
!          On exit, SUMSQ is overwritten with  smsq , the basic sum of
!          squares from which  scl  has been factored out.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   ABSXI
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
      IF( N.GT.0 ) THEN
        DO IX = 1, 1 + ( N-1 )*INCX, INCX
          IF( X( IX ).NE.ZERO ) THEN
            ABSXI = ABS( X( IX ) )
            IF( SCALE.LT.ABSXI ) THEN
              SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
              SCALE = ABSXI
            ELSE
              SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
            END IF
          END IF
        enddo
      END IF
      RETURN
!
!     End of DLASSQ
!
      END
