! ---------------------------------------------- !
!       セル中央→格子点への物理量変換関数       !
!   ただし，配列数がまちまちなので注意が必要     !
! ---------------------------------------------- !
module interpolation

	use common_hh

	implicit none

  contains

	subroutine center2grid_1(f_cen, f_grid)
		implicit none
		double precision,dimension(0:nx,0:ny),intent(out) :: f_grid
		double precision,dimension(1:nx,1:ny),intent( in) :: f_cen
		
		integer :: i,j
		
			!	領域の角の場合

		f_grid( 0, 0) = f_cen( 1, 1)
		f_grid(nx, 0) = f_cen(nx, 1)
		f_grid( 0,ny) = f_cen( 1,ny)
		f_grid(nx,ny) = f_cen(nx,ny)

			!	領域の辺の場合

!$omp do
		do i=1,nx-1
			f_grid(i, 0) = (f_cen(i, 1)+f_cen(i+1, 1))*0.5d0
			f_grid(i,ny) = (f_cen(i,ny)+f_cen(i+1,ny))*0.5d0
		end do

!$omp do
		do j=1,ny-1
			f_grid( 0,j) = (f_cen( 1,j)+f_cen( 1,j+1))*0.5d0
			f_grid(nx,j) = (f_cen(nx,j)+f_cen(nx,j+1))*0.5d0
		end do

			!	領域内部の場合

!$omp do
		do j=1,ny-1
			do i=1,nx-1
				f_grid(i,j)=(f_cen(i,j)+f_cen(i+1,j)+f_cen(i,j+1)+f_cen(i+1,j+1))*.25d0
			end do
		end do
		
	end subroutine center2grid_1

	subroutine center2grid_2(f_cen, f_grid)
		implicit none
		double precision,dimension(0:nx,0:ny),intent(out) :: f_grid
		double precision,dimension(nk,1:nx,1:ny),intent( in) :: f_cen
		
		integer :: i,j,k
		
!$omp do
		do j=0,ny
			do i=0,nx
				f_grid(i,j) = 0.d0
			end do
		end do
		
!$omp do private(k)
		do k=1,nk
			f_grid( 0, 0) = f_grid( 0, 0)+f_cen(k, 1, 1)
			f_grid(nx, 0) = f_grid(nx, 0)+f_cen(k,nx, 1)
			f_grid( 0,ny) = f_grid( 0,ny)+f_cen(k, 1,ny)
			f_grid(nx,ny) = f_grid(nx,ny)+f_cen(k,nx,ny)
		end do

				!	領域の辺の場合

!$omp do
		do i=1,nx-1
			do k=1,nk
				f_grid(i, 0) = f_grid(i, 0)+(f_cen(k,i, 1)+f_cen(k,i+1, 1))*0.5d0
				f_grid(i,ny) = f_grid(i,ny)+(f_cen(k,i,ny)+f_cen(k,i+1,ny))*0.5d0
			end do
		end do

!$omp do
		do j=1,ny-1
			do k=1,nk
				f_grid( 0,j) = f_grid( 0,j)+(f_cen(k, 1,j)+f_cen(k, 1,j+1))*0.5d0
				f_grid(nx,j) = f_grid(nx,j)+(f_cen(k,nx,j)+f_cen(k,nx,j+1))*0.5d0
			end do
		end do

				!	領域内部の場合

!$omp do
		do j=1,ny-1
			do i=1,nx-1
				do k=1,nk
					f_grid(i,j)=f_grid(i,j)+(f_cen(k,i,j)+f_cen(k,i+1,j)	&
								+f_cen(k,i,j+1)+f_cen(k,i+1,j+1))*.25d0
				end do
			end do
		end do
		
	end subroutine center2grid_2

	subroutine center2grid_3(f_cen, f_grid)
		implicit none
		double precision,dimension(0:nx,0:ny),intent(out) :: f_grid
		double precision,dimension(0:nx,0:ny+1),intent( in) :: f_cen
		
		integer :: i,j
		
			!	領域の角の場合

		f_grid( 0, 0) = f_cen( 1, 1)
		f_grid(nx, 0) = f_cen(nx, 1)
		f_grid( 0,ny) = f_cen( 1,ny)
		f_grid(nx,ny) = f_cen(nx,ny)

			!	領域の辺の場合

!$omp do
		do i=1,nx-1
			f_grid(i, 0) = (f_cen(i, 1)+f_cen(i+1, 1))*0.5d0
			f_grid(i,ny) = (f_cen(i,ny)+f_cen(i+1,ny))*0.5d0
		end do

!$omp do
		do j=1,ny-1
			f_grid( 0,j) = (f_cen( 1,j)+f_cen( 1,j+1))*0.5d0
			f_grid(nx,j) = (f_cen(nx,j)+f_cen(nx,j+1))*0.5d0
		end do

			!	領域内部の場合

!$omp do
		do j=1,ny-1
			do i=1,nx-1
				f_grid(i,j)=(f_cen(i,j)+f_cen(i+1,j)+f_cen(i,j+1)+f_cen(i+1,j+1))*.25d0
			end do
		end do
		
	end subroutine center2grid_3

	subroutine center2grid_4(f_cen, f_grid)
		implicit none
		double precision,dimension(nk,0:nx,0:ny),intent(out) :: f_grid
		double precision,dimension(nk,0:nx,0:ny+1),intent( in) :: f_cen
		
		integer :: i,j,k

!$omp single
		do k=1,nk
			f_grid(k, 0, 0) = f_cen(k, 1, 1)
			f_grid(k,nx, 0) = f_cen(k,nx, 1)
			f_grid(k, 0,ny) = f_cen(k, 1,ny)
			f_grid(k,nx,ny) = f_cen(k,nx,ny)
		end do
!$omp end single

					!	領域の辺の場合

!$omp do private(i,k)
		do i=1,nx-1
			do k=1,nk
				f_grid(k,i, 0) = (f_cen(k,i, 1)+f_cen(k,i+1, 1))*0.5d0
				f_grid(k,i,ny) = (f_cen(k,i,ny)+f_cen(k,i+1,ny))*0.5d0
			end do
		end do

!$omp do private(j,k)
		do j=1,ny-1
			do k=1,nk
				f_grid(k, 0,j) = (f_cen(k, 1,j)+f_cen(k, 1,j+1))*0.5d0
				f_grid(k,nx,j) = (f_cen(k,nx,j)+f_cen(k,nx,j+1))*0.5d0
			end do
		end do

					!	領域内部の場合

!$omp do private(i,j,k)
		do j=1,ny-1
			do i=1,nx-1
				do k=1,nk
					f_grid(k,i,j)=(f_cen(k,i,j)+f_cen(k,i+1,j)+f_cen(k,i,j+1)+f_cen(k,i+1,j+1))*.25d0
				end do
			end do
		end do
		
	end subroutine center2grid_4


	subroutine center2grid_5(f_cen, f_grid)
		implicit none
		double precision,dimension(nk,0:nx,0:ny),intent(out) :: f_grid
		double precision,dimension(nk,1:nx,1:ny),intent( in) :: f_cen
		
		integer :: i,j,k

!$omp single
		do k=1,nk
			f_grid(k, 0, 0) = f_cen(k, 1, 1)
			f_grid(k,nx, 0) = f_cen(k,nx, 1)
			f_grid(k, 0,ny) = f_cen(k, 1,ny)
			f_grid(k,nx,ny) = f_cen(k,nx,ny)
		end do
!$omp end single

					!	領域の辺の場合

!$omp do private(i,k)
		do i=1,nx-1
			do k=1,nk
				f_grid(k,i, 0) = (f_cen(k,i, 1)+f_cen(k,i+1, 1))*0.5d0
				f_grid(k,i,ny) = (f_cen(k,i,ny)+f_cen(k,i+1,ny))*0.5d0
			end do
		end do

!$omp do private(j,k)
		do j=1,ny-1
			do k=1,nk
				f_grid(k, 0,j) = (f_cen(k, 1,j)+f_cen(k, 1,j+1))*0.5d0
				f_grid(k,nx,j) = (f_cen(k,nx,j)+f_cen(k,nx,j+1))*0.5d0
			end do
		end do

					!	領域内部の場合

!$omp do private(i,j,k)
		do j=1,ny-1
			do i=1,nx-1
				do k=1,nk
					f_grid(k,i,j)=(f_cen(k,i,j)+f_cen(k,i+1,j)+f_cen(k,i,j+1)+f_cen(k,i+1,j+1))*.25d0
				end do
			end do
		end do
		
	end subroutine center2grid_5

! --------------------------------------------------------------------
	subroutine cell2grid(cell,grid)
		implicit none
		double precision,dimension(1:nx,0:ny+1),intent(in) :: cell
		double precision,dimension(0:nx,0:ny),intent(out) :: grid
		
		integer :: i,j
!
		do j=0,ny
			do i=0,nx
				if(i==0 .and. j==0) then
					grid(i,j) = cell(i+1,j+1)
				else if(i==0 .and. j==ny) then
					grid(i,j) = cell(i+1,j)
				else if(i==nx .and. j==0) then
					grid(i,j) = cell(i,j+1)
				else if(i==nx .and. j==ny) then
					grid(i,j) = cell(i,j)
				else if(i==0) then
					grid(i,j) = (cell(i+1,j)+cell(i+1,j+1))*.5
				else if(i==nx) then
					grid(i,j) = (cell(i,j)+cell(i,j+1))*.5
				else if(j==0) then
					grid(i,j) = (cell(i,j+1)+cell(i+1,j+1))*.5
				else if(j==ny) then
					grid(i,j) = (cell(i,j)+cell(i+1,j))*.5
				else
					grid(i,j) = (cell(i,j)+cell(i+1,j)+cell(i,j+1)+cell(i+1,j+1))*.25
				end if
			end do
		end do
!
	end subroutine cell2grid
!

end module interpolation