module bcg4denseMod
contains


  ! ====================================================== !
  ! === bicgstab4dense                                 === !
  ! ====================================================== !
  subroutine bicgstab4dense( Amat, bvec, xvec, LI, LJ )
    implicit none
    integer         , intent(in)    :: LI, LJ
    double precision, intent(inout) :: Amat(LI,LJ), bvec(LI), xvec(LJ)
    integer                         :: iter
    double precision                :: alpha, beta, wk
    double precision, allocatable   :: Adp(:), Ads(:)
    double precision, allocatable   :: pvec(:), rvec(:), svec(:), r0vec(:)
    integer                         :: iterMax
    double precision, parameter     :: p1 = +1.d0
    double precision, parameter     :: m1 = -1.d0
    integer         , parameter     :: incX = 1, incY = 1
    integer         , parameter     :: iterMax_Factor = 2

    ! -- BLAS Reference                                                -- !
    ! call DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
    ! call DAXPY ( n, alpha, x, incx, y, incy )
    ! see... http://www.netlib.org/blas/                               -- !
    ! -- algorithm  ( Bi-CGStab )                                      -- !
    ! see... http://www.jicfus.jp/wiki/index.php?Bi-CGSTAB%20%E6%B3%95 -- !
    ! ------------------------------------------------------------------- !
    
    ! ------------------------------------------------------ !
    ! --- [1] preparation                                --- !
    ! ------------------------------------------------------ !

    iterMax = int( max( LI, LJ ) * iterMax_Factor )
    allocate( Adp (LI), Ads (LI) )
    allocate( pvec(LI), rvec(LI), svec(LI), r0vec(LI) )
    
    ! ------------------------------------------------------ !
    ! --- [2] initial settings                           --- !
    ! ------------------------------------------------------ !
    !  -- (1) x = 0       -- !
    xvec(:)  = 0.0
    !  -- (2) r = b - Ax  -- !
    rvec(:)  = bvec(:)
    call dgemv(  "N", LI, LJ, -1.d0, Amat, LI, xvec, incX, +1.d0, rvec, incY )
    !  -- (3) r0 = r      -- !
    r0vec(:) = rvec(:)
    !  -- (4) p  = r      -- !
    pvec(:)  = rvec(:)

    ! ------------------------------------------------------ !
    ! --- [3] Main Loop                                  --- !
    ! ------------------------------------------------------ !
    !  - rvec :: rk
    !  - svec :: sk, r_k+1
    !  - pvec :: pk
    !  ----------------------------------------------------- !
    
    do iter=1, iterMax

       ! -- (6-1) calculate A.p                         -- !
       Adp      = 0.d0
       call dgemv( "N", LI, LJ, +1.d0, Amat, LI, pvec, incX, +0.d0, Adp, incY )
       
       ! -- (6-2) alpha = ( r0, r ) / ( r0, A.p )       -- !
       alpha    = inner_product( r0vec, rvec, LI ) / inner_product( r0vec, Adp , LI )

       ! -- (7)   s = r - alpha*A.p                     -- !
       call dcopy( LI, rvec, incX, svec, incY )
       call daxpy( LI, -alpha, Adp, incX, svec, incY )
       
       ! -- (8-1) calculate A.s                         -- !
       Ads = 0.d0
       call dgemv( "N", LI, LJ, +1.d0, Amat, LI, svec, incX, +0.d0, Ads, incY )
       ! -- (8-2) wk = ( A.s, s ) / ( A.s, A.s )        -- !
       wk       = inner_product( Ads, svec, LI ) / inner_product( Ads, Ads , LI )
       ! -- (9-1) x_k+1 = x_k +  alpha * p  +  wk * s   -- !
       call daxpy( LI, +alpha, pvec, incX, xvec, incY )
       call daxpy( LI, +wk   , svec, incX, xvec, incY )

       ! -- (10)  r_k+1 = s - wk * Ads                  -- !
       call daxpy( LI, -wk, Ads, incX, svec, incY )
       ! -- (11) beta = alpha/wk*(r0,r_k+1)/(r0,r_k)    -- !
       beta     = alpha / wk * inner_product( r0vec, svec, LI ) / &
            &                  inner_product( r0vec, rvec, LI )
       ! -- (12) p = r_k+1 + beta * ( p - wk * A.p )    -- !
       ! pvec(:)  = rvec(:) + beta*( pvec(:) - wk*Adp(:) )
       call dcopy( LI, svec, incX, rvec, incY )
       call daxpy( LI, -wk , Adp , incX, pvec, incY )
       call daxpy( LI, beta, pvec, incX, svec, incY )
       call dcopy( LI, svec, incX, pvec, incY )
       
    enddo
    
    return
  end subroutine bicgstab4dense

  
  ! ====================================================== !
  ! === Function :: inner product :: x.y               === !
  ! ====================================================== !
  Function inner_product( xvec, yvec, Npt )
    implicit none
    integer         , intent(in) :: Npt
    double precision, intent(in) :: xvec(Npt), yvec(Npt)
    integer                      :: i
    double precision             :: inner_product
    
    inner_product = 0.d0
    do i=1, Npt
       inner_product = inner_product + xvec(i) * yvec(i)
    enddo
    
    return
  end Function inner_product
  
end module bcg4denseMod
