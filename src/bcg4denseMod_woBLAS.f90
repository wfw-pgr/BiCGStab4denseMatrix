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
    double precision, allocatable   :: Adp(:), Adx(:), Ads(:)
    double precision, allocatable   :: pvec(:), rvec(:), svec(:), r0vec(:)
    integer                         :: iterMax
    double precision, parameter     :: p1 = +1.d0
    double precision, parameter     :: m1 = -1.d0
    integer         , parameter     :: incX = 1, incY = 1
    integer         , parameter     :: iterMax_Factor = 2

    ! -- algorithm  ( Bi-CGStab )                                      -- !
    ! see... http://www.jicfus.jp/wiki/index.php?Bi-CGSTAB%20%E6%B3%95 -- !
    ! ------------------------------------------------------------------- !
    
    ! ------------------------------------------------------ !
    ! --- [1] preparation                                --- !
    ! ------------------------------------------------------ !

    iterMax = int( max( LI, LJ ) * iterMax_Factor )
    allocate( Adx (LI), Adp (LI), Ads (LI) )
    allocate( pvec(LI), rvec(LI), svec(LI), r0vec(LI) )
    
    ! ------------------------------------------------------ !
    ! --- [2] initial settings                           --- !
    ! ------------------------------------------------------ !
    !  -- (1) x = 0       -- !
    xvec(:)  = 0.0
    !  -- (2) r = b - Ax  -- !
    rvec(:)  = bvec(:) - matmul( Amat, xvec )
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
       Adp      = matmul( Amat, pvec )
       
       ! -- (6-2) alpha = ( r0, r ) / ( r0, A.p )       -- !
       alpha    = dot_product( r0vec, rvec ) / dot_product( r0vec, Adp  )
       
       ! -- (7)   s = r - alpha*A.p                     -- !
       svec(:)  = rvec(:) - alpha * Adp(:)
       
       ! -- (8-1) calculate A.s                         -- !
       Ads(:)   = matmul( Amat, svec )
       
       ! -- (8-2) wk = ( A.s, s ) / ( A.s, A.s )        -- !
       wk       = dot_product( Ads, svec ) / dot_product( Ads, Ads )
       ! -- (9-1) x_k+1 = x_k +  alpha * p  +  wk * s   -- !
       xvec(:)  = xvec(:) + alpha*pvec(:) + wk*svec(:)
       
       ! -- (10)  r_k+1 = s - wk * Ads                  -- !
       svec(:)  = svec(:) - wk * Ads(:)
       
       ! -- (11) beta = alpha/wk*(r0,r_k+1)/(r0,r_k)    -- !
       beta     = alpha / wk * dot_product( r0vec, svec ) / &
            &                  dot_product( r0vec, rvec )
       ! -- (12) p = r_k+1 + beta * ( p - wk * A.p )    -- !
       rvec(:)  = svec(:)
       pvec(:)  = rvec(:) + beta*( pvec(:) - wk*Adp(:) )
       
    enddo
    
    return
  end subroutine bicgstab4dense



end module bcg4denseMod
