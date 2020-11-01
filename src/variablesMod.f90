module variablesMod
  
  implicit none
  integer         , parameter   :: cLen     = 300
  integer         , parameter   :: lun      = 50
  double precision, parameter   :: eps      = 1.e-8
  character(cLen) , parameter   :: AmatFile = "dat/Amat.dat"
  character(cLen) , parameter   :: bvecFile = "dat/bvec.dat"
  character(cLen) , parameter   :: xvecFile = "dat/xvec.dat"  
  integer                       :: LI, LJ
  double precision, allocatable :: Amat(:,:), bvec(:), xvec(:)
  
end module variablesMod
