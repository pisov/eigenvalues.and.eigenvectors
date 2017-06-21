!**************************************************
!* Hermite polynomial coefficients evaluation by  *
!* means of recursion relation. The order of the  *
!* polynomial is n. The coefficients are returned *
!* in A(i).                                       *
!**************************************************
subroutine hermite_coeff(n,A,B)  
  integer, intent(in) :: n
  double precision, dimension(:), intent(inout)   :: A
  double precision, dimension(:,:), intent(inout) :: B
  
  integer :: i, j
  
  !Establish l0 and l1 coefficients
  B(0,0)=1.d0 ; B(1,0)=0.d0 ; B(1,1)=2.d0
  !Return if order is less than two
  if (n>1) then
    do i = 2, n
      B(i,0)=-2.d0*(i-1)*B(i-2,0)
      do j = 1, i
        !Basic recursion relation
        B(i,j)=2.d0*B(i-1,j-1)-2.d0*(i-1)*B(i-2,j)
      end do
    end do
    do i = 0, n
      A(i)=B(n,i)
    end do
  end if
  return
end
