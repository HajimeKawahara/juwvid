module ionufft
  implicit none
contains

  subroutine ionufft1d2(nj,ms,iflag,eps,mx,xj,fk,cj)
    implicit none
    integer, intent(in) :: ms,nj,iflag
    integer, intent(in) :: mx !=10000
    complex*16,intent(in) :: fk(mx)
    real*8,intent(in) :: xj(mx)
    real*8,intent(in) :: eps
    complex*16,intent(out) :: cj(mx)
    integer :: ier
    
    call nufft1d2f90(nj,xj,cj,iflag,eps,ms,fk,ier)
    
  end subroutine ionufft1d2


end module ionufft
