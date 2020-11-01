program main

  use variablesMod
  use ioUtilityMod
  use bcg4denseMod
  implicit none

  call load__MatrixVector
  call bicgstab4dense( Amat, bvec, xvec, LI, LJ )
  call save__xvector
  

end program main
  
