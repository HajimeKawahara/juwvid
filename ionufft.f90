module ionufft
  implicit none
contains

  subroutine ionufft1d2(nj,ms,mx,iflag,eps,xj,fk,cj)
    implicit none
    integer, intent(in) :: ms,nj,mx,iflag
    complex*16,intent(in) :: fk(ms)
    real*8,intent(in) :: xj(nj)
    real*8,intent(in) :: eps
    complex*16,intent(out) :: cj(mx)
    integer :: ier
    
    call nufft1d2f90(nj,xj,cj,iflag,eps,ms,fk,ier)
    
  end subroutine ionufft1d2

  subroutine ionufft1d3(nj,nk,iflag,eps,xj,sk,cj,fk)
    implicit none
    integer, intent(in) :: nj,nk,iflag
    real*8,intent(in) :: xj(nj)
    real*8,intent(in) :: sk(nk)
    real*8,intent(in) :: eps
    complex*16,intent(out) :: fk(nk)
    complex*16,intent(in) :: cj(nj)
    integer :: ier
    
    call nufft1d3f90(nj,xj,cj,iflag,eps,nk,sk,fk,ier)
    
  end subroutine ionufft1d3

end module ionufft
