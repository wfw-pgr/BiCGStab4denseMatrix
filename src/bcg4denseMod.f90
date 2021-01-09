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
    double precision                :: alpha, beta, wk, criterion, rNorm
    double precision, allocatable   :: Adp(:), Ads(:)
    double precision, allocatable   :: pvec(:), rvec(:), svec(:), r0vec(:)
    integer                         :: iterMax
    integer         , parameter     :: incX = 1, incY = 1
    integer         , parameter     :: iterMax_Factor = 2
    double precision, parameter     :: convergence    = 1.d-12

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
    allocate( pvec(LJ), rvec(LJ), svec(LJ), r0vec(LJ) )

    ! ------------------------------------------------------ !
    ! --- [2] initial settings                           --- !
    ! ------------------------------------------------------ !
    !  -- (1) x = 0       -- !
    xvec(:)   = 0.0
    !  -- (2) r = b - Ax  -- !
    rvec(:)   = bvec(:)
    call dgemv(  "N", LI, LJ, -1.d0, Amat, LI, xvec, incX, +1.d0, rvec, incY )
    !  -- (3) r0=r, p=r   -- !
    r0vec(:)  = rvec(:)
    pvec(:)   = rvec(:)
    criterion = inner_product( bvec, bvec, LI ) * convergence

    ! ------------------------------------------------------ !
    ! --- [3] Main Loop                                  --- !
    ! ------------------------------------------------------ !
    !  - rvec :: rk
    !  - svec :: sk, r_k+1
    !  - pvec :: pk
    !  ----------------------------------------------------- !
    
    do iter=1, iterMax

       ! -- convergence check :: Loop End                 -- !
       rNorm = inner_product( rvec, rvec, LJ )
       if ( rNorm.lt.criterion ) then
          write(6,*)
          write(6,*) "[bcg4denseMod.f90] BiCGStab Method reached convergence.  "
          write(6,*) "[bcg4denseMod.f90] residual        :: ", rNorm
          write(6,*) "[bcg4denseMod.f90] convergence     :: ", convergence
          write(6,*)
          exit
       endif
       
       ! -- (4) calculate A.p                             -- !
       Adp      = 0.d0
       call dgemv( "N", LI, LJ, +1.d0, Amat, LI, pvec, incX, +0.d0, Adp, incY )
       
       ! -- (5) alpha = ( r0, r ) / ( r0, A.p )           -- !
       alpha    = inner_product( r0vec, rvec, LJ ) / inner_product( r0vec, Adp , LJ )

       ! -- (6)   s = r - alpha*A.p                       -- !
       call dcopy( LI, rvec, incX, svec, incY )
       call daxpy( LI, -alpha, Adp, incX, svec, incY )
       
       ! -- (7) calculate A.s                             -- !
       Ads = 0.d0
       call dgemv( "N", LI, LJ, +1.d0, Amat, LI, svec, incX, +0.d0, Ads, incY )
       ! -- (8) wk = ( A.s, s ) / ( A.s, A.s )            -- !
       wk       = inner_product( Ads, svec, LJ ) / inner_product( Ads, Ads , LJ )
       ! -- (9) x_k+1 = x_k +  alpha * p  +  wk * s       -- !
       call daxpy( LI, +alpha, pvec, incX, xvec, incY )
       call daxpy( LI, +wk   , svec, incX, xvec, incY )

       ! -- (10)  r_k+1 = s - wk * Ads                    -- !
       call daxpy( LI, -wk, Ads, incX, svec, incY )
       ! -- (11) beta = alpha/wk*(r0,r_k+1)/(r0,r_k)      -- !
       beta     = alpha / wk * inner_product( r0vec, svec, LJ ) / &
            &                  inner_product( r0vec, rvec, LJ )
       ! -- (12) p = r_k+1 + beta * ( p - wk * A.p )      -- !
       ! pvec(:)  = rvec(:) + beta*( pvec(:) - wk*Adp(:) )
       call dcopy( LJ, svec, incX, rvec, incY )
       call daxpy( LJ, -wk , Adp , incX, pvec, incY )
       call daxpy( LJ, beta, pvec, incX, svec, incY )
       call dcopy( LJ, svec, incX, pvec, incY )
       
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
