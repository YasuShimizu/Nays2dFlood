! --------------------------------------- !
!      —¬‘¬C…[‚Ì•ÏŠ·ƒ‚ƒWƒ…[ƒ‹         !
! --------------------------------------- !

module transform

	use common_hh
	
	implicit none
	
	double precision,dimension(:,:),allocatable :: up0,vp0
	
  contains

	! --------------------------------------------------------------------
	subroutine uvpcal(u,v,up,vp,hs)
		implicit none

		double precision,dimension(-1:nx+1,-1:ny+1),intent(in) :: u,v
		double precision,dimension(1:nx,0:ny+1) :: hs
		double precision,dimension(1:nx,0:ny+1),intent(out) :: up,vp

		integer :: i,j
!
!$omp do
		do j=1-j_side1,ny+j_side2
			do i=1,nx
				if( hs(i,j)<=hmin ) then
					up(i,j) = 0.
					vp(i,j) = 0.
				else
					up(i,j) = (u(i-1,j)+u(i,j))*.5d0
					vp(i,j) = (v(i,j-1)+v(i,j))*.5d0
				end if
			end do
		end do
!
	end subroutine uvpcal

       
! --------------------------------------------------------------------
	subroutine alloc_uxxyycal_temp_variables
		implicit none
          
		allocate( up0(0:nx,0:ny), vp0(0:nx,0:ny) )
          
		up0 = 0.d0;	vp0 = 0.d0
        
	end subroutine alloc_uxxyycal_temp_variables
!
! ---------------------------------------------------------------------
	subroutine uxxyycal(u0,v0,ux0,uy0)
		use variables
		implicit none
		double precision,dimension(-1:nx+1,-1:ny+1),intent(in) :: u0,v0
		double precision,dimension(0:nx,0:ny),intent(out) :: ux0,uy0
		integer(4) :: i,j,n

		if( j_side1==0 ) then
!$omp do
			do i=0,nx
				up0(i,0) = u0(i,1)
				vp0(i,0) = 0.d0
			end do
		else
!$omp do
			do i=0,nx
				up0(i,0) = u0(i,1)
				vp0(i,0) = (v0(i,0)+v0(i+1,0))*.5d0
			end do
		end if

		if( j_side2==0 ) then
!$omp do
			do i=0,nx
				up0(i,ny) = u0(i,ny)
				vp0(i,ny) = 0.d0
			end do
		else
!$omp do
			do i=0,nx
				up0(i,ny) = u0(i,ny)
				vp0(i,ny) = (v0(i,ny)+v0(i+1,ny))*.5d0
			end do
		end if
!
!$omp do
		do j=1,ny-1
			do i=0,nx
				up0(i,j) = (u0(i,j)+u0(i,j+1))*.5d0
				vp0(i,j) = (v0(i,j)+v0(i+1,j))*.5d0
			end do
		end do
!
!$omp single

		if( j_side_j1==2 ) then
		
			do i=1,nx
				vp0(i,0) = 0.d0
			end do
			
			do n=1,jsin1
				do i=j1s(n),j1e(n)
					vp0(i,0) = (v0(i,0)+v0(i+1,0))*.5
				end do
			end do
		end if
!
		if( j_side_j2==2 ) then
		
			do i=1,nx
				vp0(i,ny) = 0.d0
			end do
			
			do n=1,jsin2
				do i=j2s(n),j2e(n)
					vp0(i,ny) = (v0(i,ny)+v0(i+1,ny))*.5
				end do
			end do
		end if

!$omp end single
!
!$omp do
		do j=0,ny
			do i=0,nx
				ux0(i,j) = ( et_y0(i,j)*up0(i,j)-xi_y0(i,j)*vp0(i,j))*r_sj0(i,j)
				uy0(i,j) = (-et_x0(i,j)*up0(i,j)+xi_x0(i,j)*vp0(i,j))*r_sj0(i,j)
			end do
		end do
!
	end subroutine uxxyycal

! ---------------------------------------------------------------------
	subroutine uxuycal(up0,vp0,ux0,uy0)
		use variables
		implicit none

		double precision,dimension(1:nx,0:ny+1),intent(in) :: up0,vp0
		double precision,dimension(1:nx,0:ny+1),intent(out) :: ux0,uy0

		integer :: i,j
!
!$omp do
		do j=1-j_side1,ny+j_side2
			do i=1,nx
				ux0(i,j) = ( et_y(i,j)*up0(i,j)-xi_y(i,j)*vp0(i,j))*r_sj(i,j)
				uy0(i,j) = (-et_x(i,j)*up0(i,j)+xi_x(i,j)*vp0(i,j))*r_sj(i,j)
			end do
		end do
!
	end subroutine uxuycal

! ---------------------------------------------------------------------
	subroutine hsxxcal(eta,z,hs,hsxx,rain,rainxx, ibc0, ibcg0)
		implicit none
        double precision,dimension(1:nx,0:ny+1),intent(in) :: eta,hs
		double precision,dimension(0:nx,0:ny),intent(out) :: z,hsxx,rainxx
		double precision,dimension(1:nx,1:ny),intent(in) :: rain
        integer,dimension(1:nx,0:ny+1),intent(in) :: ibc0
        integer,dimension(0:nx,0:ny),intent(out) :: ibcg0
		
		integer :: i,j
!
			!	—Ìˆæ‚ÌŠp‚Ìê‡

		hsxx( 0, 0) = hs( 1, 1)
		hsxx(nx, 0) = hs(nx, 1)
		hsxx( 0,ny) = hs( 1,ny)
		hsxx(nx,ny) = hs(nx,ny)

		z( 0, 0) = eta( 1, 1)
		z(nx, 0) = eta(nx, 1)
		z( 0,ny) = eta( 1,ny)
		z(nx,ny) = eta(nx,ny)

		rainxx( 0, 0) = rain( 1, 1)
		rainxx(nx, 0) = rain(nx, 1)
		rainxx( 0,ny) = rain( 1,ny)
		rainxx(nx,ny) = rain(nx,ny)
        
        ibcg0( 0, 0) = ibc0( 1, 1)
		ibcg0(nx, 0) = ibc0(nx, 1)
		ibcg0( 0,ny) = ibc0( 1,ny)
		ibcg0(nx,ny) = ibc0(nx,ny)

			!	—Ìˆæ‚Ì•Ó‚Ìê‡

!$omp do
		do i=1,nx-1
			hsxx(i, 0) = (hs(i, 1)+hs(i+1, 1))*0.5d0
			hsxx(i,ny) = (hs(i,ny)+hs(i+1,ny))*0.5d0
			z(i, 0) = (eta(i, 1)+eta(i+1, 1))*0.5d0
			z(i,ny) = (eta(i,ny)+eta(i+1,ny))*0.5d0
			rainxx(i, 0) = (rain(i, 1)+rain(i+1, 1))*0.5d0
			rainxx(i,ny) = (rain(i,ny)+rain(i+1,ny))*0.5d0
            
            ibcg0(i, 0) = (ibc0(i, 1)+ibc0(i+1, 1))*0.5d0
			ibcg0(i,ny) = (ibc0(i,ny)+ibc0(i+1,ny))*0.5d0
		end do

!$omp do
		do j=1,ny-1
			hsxx( 0,j) = (hs( 1,j)+hs( 1,j+1))*0.5d0
			hsxx(nx,j) = (hs(nx,j)+hs(nx,j+1))*0.5d0
			z( 0,j) = (eta( 1,j)+eta( 1,j+1))*0.5d0
			z(nx,j) = (eta(nx,j)+eta(nx,j+1))*0.5d0
			rainxx( 0,j) = (rain( 1,j)+rain( 1,j+1))*0.5d0
			rainxx(nx,j) = (rain(nx,j)+rain(nx,j+1))*0.5d0
            
            ibcg0( 0,j) = (ibc0( 1,j)+ibc0( 1,j+1))*0.5d0
			ibcg0(nx,j) = (ibc0(nx,j)+ibc0(nx,j+1))*0.5d0
		end do

			!	—Ìˆæ“à•”‚Ìê‡

!$omp do
		do j=1,ny-1
			do i=1,nx-1
				hsxx(i,j)=(hs(i,j)+hs(i+1,j)+hs(i,j+1)+hs(i+1,j+1))*.25d0
				z(i,j)=(eta(i,j)+eta(i+1,j)+eta(i,j+1)+eta(i+1,j+1))*.25d0
				rainxx(i,j)=(rain(i,j)+rain(i+1,j)+rain(i,j+1)+rain(i+1,j+1))*.25d0
                
                ibcg0(i,j)=(ibc0(i,j)+ibc0(i+1,j)+ibc0(i,j+1)+ibc0(i+1,j+1))*.25d0
			end do
		end do
!
	end subroutine hsxxcal

end module transform