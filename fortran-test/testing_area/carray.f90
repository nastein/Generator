subroutine fortfunc(array)
      implicit none
      complex*16 :: array(4,4)
      integer :: i,j,k
      k = 1.0d0
      do i = 1,4
        do j = 1,4
          array(i,j) = k
          k=k+1
          !write(6,*) array(i,j)
                end do 
          end do

      !Array before shifting
      write(6,*) array(:,1)
      write(6,*) array(:,2)
      write(6,*) array(:,3)
      write(6,*) array(:,4)

      array = CSHIFT(array, 1, DIM=1)
      array = CSHIFT(array, 1, DIM=2)

          array = TRANSPOSE(array)

          !do j=1,3
          !     write(6,*) array(:,j)
          !end do

      return
end
