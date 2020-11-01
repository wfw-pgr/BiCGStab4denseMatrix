module ioUtilityMod
contains

  
  ! ====================================================== !
  ! === load__MatrixVector                             === !
  ! ====================================================== !
  subroutine load__MatrixVector
    use variablesMod
    implicit none
    integer :: i, j, LI_

    ! ------------------------------------------------------ !
    ! --- [1] load Amatrix                               --- !
    ! ------------------------------------------------------ !
    
    open(lun,file=trim(AmatFile),status="old")
    read(lun,*) LI, LJ
    allocate( Amat(LI,LJ), xvec(LJ) )
    xvec(:) = 0.d0
    do j=1, LJ
       read(lun,*) Amat(:,j)
    enddo
    close(lun)
    write(6,*) "[load__MatrixVector] AmatFile is loaded...."

    ! ------------------------------------------------------ !
    ! --- [2] load bvector                               --- !
    ! ------------------------------------------------------ !
    open(lun,file=trim(bvecFile),status="old")
    read(lun,*) LI_
    if ( LI_.ne.LI ) then
       write(6,*) "[load__bvector] LI_ != LI :: incompatible size... [ERROR]"
       stop
    endif
    allocate( bvec(LI) )
    do i=1, LI
       read(lun,*) bvec(i)
    enddo
    close(lun)
    write(6,*) "[load__MatrixVector] bvecFile is loaded...."

    ! ------------------------------------------------------ !
    ! --- [3] display matrix and bvector                 --- !
    ! ------------------------------------------------------ !
    write(6,*) "[load__MatrixVector] loaded matrix and bvector..."
    do i=1, LI
       write(6,*) ( Amat(i,j), j=1, LJ ), " | ", bvec(i)
    enddo
    write(6,*)
    
    return
  end subroutine load__MatrixVector
  
  

  ! ====================================================== !
  ! === save xvector in a file                         === !
  ! ====================================================== !
  subroutine save__xvector
    use variablesMod
    implicit none
    integer :: i

    ! ------------------------------------------------------ !
    ! --- [1] save xvector                               --- !
    ! ------------------------------------------------------ !
    open(lun,file=trim(xvecFile),status="replace")
    do i=1, LI
       write(lun,*) xvec(i)
       write(6  ,*) xvec(i)
    enddo
    close(lun)
    write(6,*) "[save__xvector] xvecFile is saved...."
    

    
    return
  end subroutine save__xvector
  
end module ioUtilityMod
