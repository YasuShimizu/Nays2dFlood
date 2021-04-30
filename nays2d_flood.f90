! ------------------------------------------------------ !
!         Nays2D Flood for parallel computation          !
!            modified by Toshiki Iwasaki from 2012-10-3  !
! ------------------------------------------------------ !

! ---- module ---
!
module avgeo_m			! --------------------------------------------------------------------
	use common_hh
	use variables
	implicit none

	double precision, parameter :: etmax = 9999.d0

	contains
	!
	! ---------------------------------------------- !
	subroutine avgeo
		implicit none
		integer(4) :: i,j,net,nnp,mv
		double precision :: eett, ee00, dnx0
!
		do j=1,ny
			do i=1,nx
				net=0
				eett=0.

				if( z(i,j)<etmax ) then
					net = net+1
					eett = eett+z(i,j)
				end if

				if( z(i-1,j)<etmax ) then
					net = net+1
					eett = eett+z(i-1,j)
				end if

				if(z(i,j-1)<etmax) then
					net=net+1
					eett=eett+z(i,j-1)
				end if

				if(z(i-1,j-1)<etmax) then
					net=net+1
					eett=eett+z(i-1,j-1)
				end if

				if(net>0) then
					eta(i,j) = eett/float(net)
				else
					eta(i,j) = etmax
				end if

				eta0(i,j) = eta(i,j)
			end do
		end do
!
		do i=1,nx
			eave(i) = 0.
			emin(i) = 9999.
			emax(i) = -9999.
			nnp=0
			
			do j=1,ny
				if( ijo_in(i,j)/=1 ) then
					nnp = nnp+1
					eave(i) = eave(i)+eta(i,j)
					emin(i) = min(emin(i),eta(i,j))
					emax(i) = max(emax(i),eta(i,j))
				end if
			end do
			eave(i) = eave(i)/dble(nnp)
		end do
!
		width = 0.
		
		do i=0,nx
			chb(i)=0.
			do j=1,ny
				if( ijobst(i,j)*ijobst(i,j-1)==0 ) then
					dnx0 = sqrt((x(i,j)-x(i,j-1))**2d0+(y(i,j)-y(i,j-1))**2d0)
					chb(i) = chb(i)+dnx0
				end if
			end do
			width = width+chb(i)
		end do
		
		width = width/dble(nx+1)

	end subroutine avgeo

end module avgeo_m		! -----------------------------------------------------------------


module gcoefs_m			! -----------------------------------------------------------------

	use common_hh
	use variables
	implicit none

	integer, parameter :: ltime = 0
	double precision, parameter :: ds00 = 0.001
	double precision, parameter :: dn00 = 0.001
	double precision,dimension(:,:),allocatable :: x_xi, x_et, y_xi, y_et
	double precision,dimension(:,:),allocatable :: btmp

	contains

! -------------------------------------------------------------------

	subroutine alloc_gcoefs_temp_variables
		implicit none
		integer :: i,j

		allocate( x_xi(-1:nx+1,-1:ny+1), x_et(-1:nx+1,-1:ny+1) )
		allocate( y_xi(-1:nx+1,-1:ny+1), y_et(-1:nx+1,-1:ny+1) )
		allocate( btmp(-1:nx+1,-1:ny+1) )

!$omp do
		do j=-1,ny+1
			do i=-1,nx+1
				x_xi(i,j) = 0.d0
				x_et(i,j) = 0.d0
				y_xi(i,j) = 0.d0
				y_et(i,j) = 0.d0
				btmp(i,j) = 0.d0
			end do
		end do

	end subroutine alloc_gcoefs_temp_variables
!
! ---------------------------------------------------------------------
	subroutine gcoefs(iout)
		implicit none
		integer :: i,j,iout, m, l
 		double precision :: dx, dy, theta 			&
			, x1, y1, x2, y2, x3, y3, x4, y4 		&
			, dx1, dy1, dx2, dy2, ds1, ds2, ds12 	&
			, x_xi0, y_xi0, x_et0, y_et0 			&
			, x_xixi, y_xixi, x_xiet, y_xiet, x_etet, y_etet	&
			, s1, s2, area1, area2

		dx=0.;	dy=0.;	theta=0.
!
		do j=0,ny
			do i=1,nx
				ds(i,j) = dsqrt((x(i,j)-x(i-1,j))**2d0+(y(i,j)-y(i-1,j))**2d0)
				ds(i,j) = max(ds(i,j),ds00)
				xi_r(i,j) = dxi/ds(i,j)
			end do
		end do
!
		do i=1,nx
			ds(i,-1) = ds(i,0)
			ds(i,ny+1) = ds(i,ny)
			xi_r(i,-1) = xi_r(i,0)
			xi_r(i,ny+1) = xi_r(i,ny)
		end do
!
		do j=1-j_side1,ny+j_side2
			do i=1,nx-1
				xi_r_up(i,j) = (xi_r(i,j)+xi_r(i,j-1)+xi_r(i+1,j)+xi_r(i+1,j-1))*.25d0
			end do
		end do

		do j=1,ny
			xi_r_up( 0,j) = (xi_r( 1,j)+xi_r( 1,j-1))*.5d0
			xi_r_up(nx,j) = (xi_r(nx,j)+xi_r(nx,j-1))*.5d0
		end do
!
		do j=1,ny
			do i=1,nx-1
				dsy(i,j) = (ds(i,j)+ds(i+1,j)+ds(i,j-1)+ds(i+1,j-1))*0.25d0
			end do
			dsy(0,j) = dsy(1,j)*0.5d0
		end do

		if( j_wl==3 ) then
			do j=1,ny
				dsy(nx,j) = dsy(nx-1,j)
			end do
		else
			do j=1,ny
				dsy(nx,j) = dsy(nx-1,j)*0.5d0
			end do
		end if
!
		do j=1,ny
			do i=0,nx
				dn(i,j) = dsqrt((x(i,j)-x(i,j-1))**2d0+(y(i,j)-y(i,j-1))**2d0)
				dn(i,j) = max(dn(i,j),dn00)
				et_r(i,j) = det/dn(i,j)
			end do
		end do

		do i=0,nx
			dn(i,0) = dn(i,1)
			dn(i,ny+1) = dn(i,ny)
			et_r(i,0) = et_r(i,1)
			et_r(i,ny+1) = et_r(i,ny)
		end do
!
		do j=1-j_side1,ny-1+j_side2
			do i=1,nx
				et_r_vp(i,j) = (et_r(i,j)+et_r(i-1,j)+et_r(i,j+1)+et_r(i-1,j+1))*.25d0
			end do
		end do
!
		if( j_side1==0 ) then
			do i=1,nx
				et_r_vp(i,0) = (et_r(i,1)+et_r(i-1,1))*.5d0
			end do
		else
			do i=1,nx
				et_r_vp(i,-1) = et_r_vp(i,0)
			end do
		end if

		if( j_side2==0 ) then
			do i=1,nx
				et_r_vp(i,ny) = (et_r(i,ny)+et_r(i-1,ny))*.5d0
			end do
		else
			do i=1,nx
				et_r_vp(i,ny+1) = et_r_vp(i,ny)
			end do
		end if
!
		do j=1,ny-1
			do i=1,nx
				dnx(i,j) = (dn(i,j)+dn(i,j+1)+dn(i-1,j)+dn(i-1,j+1))*.25d0
			end do
		end do

		if( j_side1==1 ) then
			do i=1,nx
				dnx(i,0) = dnx(i,1)
			end do
		else
			do i=1,nx
				dnx(i,0) = dnx(i,1)*0.5d0
			end do
		end if

		if( j_side2==1 ) then
			do i=1,nx
				dnx(i,ny) = dnx(i,ny-1)
			end do
		else
			do i=1,nx
				dnx(i,ny) = dnx(i,ny-1)*0.5d0
			end do
		end if
		
		do j=1,ny
			do i=1,nx
				dsn(i,j) = dsqrt( (x(i-1,j)-x(i,j-1))**2.+(y(i-1,j)-y(i,j-1))**2. )
			end do
		end do
		
		do j=1,ny
			do i=1,nx
				s1 = 0.5d0*(dn(i-1,j)+ds(i,j-1)+dsn(i,j))
				s2 = 0.5d0*(ds(i,j)+dn(i,j)+dsn(i,j))
				area1 = dsqrt(s1*(s1-dn(i-1,j))*(s1-ds(i,j-1))*(s1-dsn(i,j)))
				area2 = dsqrt(s2*(s2-dn(i,j))*(s2-ds(i,j))*(s2-dsn(i,j)))
				area(i,j) = area1+area2
			end do
		end do
		
!
		do j=1,ny
			do i=1,nx
				x_xi(i,j) = (x(i,j)+x(i,j-1)-x(i-1,j)-x(i-1,j-1))*0.5d0*r_dxi
				x_et(i,j) = (x(i,j)+x(i-1,j)-x(i,j-1)-x(i-1,j-1))*0.5d0*r_det
				y_xi(i,j) = (y(i,j)+y(i,j-1)-y(i-1,j)-y(i-1,j-1))*0.5d0*r_dxi
				y_et(i,j) = (y(i,j)+y(i-1,j)-y(i,j-1)-y(i-1,j-1))*0.5d0*r_det
				sj(i,j) = 1.d0/(x_xi(i,j)*y_et(i,j)-x_et(i,j)*y_xi(i,j))
				xi_x(i,j) = sj(i,j)*y_et(i,j)
				xi_y(i,j) = -sj(i,j)*x_et(i,j)
				et_x(i,j) = -sj(i,j)*y_xi(i,j)
				et_y(i,j) = sj(i,j)*x_xi(i,j)
			end do
		end do
!
		do j=1,ny
			sj(0,j)=sj(1,j)
		end do

		do i=0,nx
			xi_x(i,0) = xi_x(i,1)
			xi_x(i,ny+1) = xi_x(i,ny)
			xi_y(i,0) = xi_y(i,1)
			xi_y(i,ny+1) = xi_y(i,ny)
			et_x(i,0) = et_x(i,1)
			et_x(i,ny+1) = et_x(i,ny)
			et_y(i,0) = et_y(i,1)
			et_y(i,ny+1) = et_y(i,ny)
			sj(i,0) = sj(i,1)
			sj(i,ny+1) = sj(i,ny)

			x_xi(i,0) = x_xi(i,1)
			x_xi(i,ny+1) = x_xi(i,ny)
			y_xi(i,0) = y_xi(i,1)
			y_xi(i,ny+1) = y_xi(i,ny)
			x_et(i,0) = x_et(i,1)
			x_et(i,ny+1) = x_et(i,ny)
			y_et(i,0) = y_et(i,1)
			y_et(i,ny+1) = y_et(i,ny)
		end do
		
		do j=0,ny+1
			do i=0,nx
				r_sj(i,j) = 1.d0/sj(i,j)
			end do
		end do
!
		do j=0,ny
			do i=0,nx
				if( i==0 ) then
					x_xi0 = (x(i+1,j)-x(i,j))/dxi
					y_xi0 = (y(i+1,j)-y(i,j))/dxi
				else if( i==nx ) then
					x_xi0 = (x(i,j)-x(i-1,j))/dxi
					y_xi0 = (y(i,j)-y(i-1,j))/dxi
				else
					x_xi0 = (x(i+1,j)-x(i-1,j))/(2.d0*dxi)
					y_xi0 = (y(i+1,j)-y(i-1,j))/(2.d0*dxi)
				end if
				if( j==0 ) then
					x_et0 = (x(i,j+1)-x(i,j))/det
					y_et0 = (y(i,j+1)-y(i,j))/det
				else if( j==ny ) then
					x_et0 = (x(i,j)-x(i,j-1))/det
					y_et0 = (y(i,j)-y(i,j-1))/det
				else
					x_et0 = (x(i,j+1)-x(i,j-1))/(2.d0*det)
					y_et0 = (y(i,j+1)-y(i,j-1))/(2.d0*det)
				end if
				sj0(i,j) = 1.d0/(x_xi0*y_et0-x_et0*y_xi0)
				xi_x0(i,j) =  sj0(i,j)*y_et0
				xi_y0(i,j) = -sj0(i,j)*x_et0
				et_x0(i,j) = -sj0(i,j)*y_xi0
				et_y0(i,j) =  sj0(i,j)*x_xi0
			end do
		end do
		
		do j=0,ny
			do i=0,nx
				r_sj0(i,j) = 1.d0/sj0(i,j)
			end do
		end do
!
! ---- coefs at u grid point ------
!
		do i=1,nx-1
			do j=1,ny
				xi_x_up(i,j) = (xi_x(i,j)+xi_x(i+1,j))*.5d0
				et_x_up(i,j) = (et_x(i,j)+et_x(i+1,j))*.5d0
				xi_y_up(i,j) = (xi_y(i,j)+xi_y(i+1,j))*.5d0
				et_y_up(i,j) = (et_y(i,j)+et_y(i+1,j))*.5d0
				beta(i,j,1) = xi_x_up(i,j)**2d0+xi_y_up(i,j)**2d0
				beta(i,j,2) = xi_x_up(i,j)*et_x_up(i,j)+xi_y_up(i,j)*et_y_up(i,j)
				x_xixi = (x_xi(i+1,j)-x_xi(i,j))/dxi
				y_xixi = (y_xi(i+1,j)-y_xi(i,j))/dxi
				x_xiet = (x_et(i+1,j)-x_et(i,j))/dxi
				y_xiet = (y_et(i+1,j)-y_et(i,j))/dxi
				
				if( j==1 ) then
					x_etet = (x_et(i+1,j+1)+x_et(i,j+1)-x_et(i+1,j)-x_et(i,j))/(2.d0*det)
					y_etet = (y_et(i+1,j+1)+y_et(i,j+1)-y_et(i+1,j)-y_et(i,j))/(2.d0*det)
				else if( j==ny ) then
					x_etet = (x_et(i+1,j)+x_et(i,j)-x_et(i+1,j-1)-x_et(i,j-1))/(2.d0*det)
					y_etet = (y_et(i+1,j)+y_et(i,j)-y_et(i+1,j-1)-y_et(i,j-1))/(2.d0*det)
				else
					x_etet = (x_et(i+1,j+1)+x_et(i,j+1)-x_et(i+1,j-1)-x_et(i,j-1))/(4.d0*det)
					y_etet = (y_et(i+1,j+1)+y_et(i,j+1)-y_et(i+1,j-1)-y_et(i,j-1))/(4.d0*det)
				end if

				alpha(i,j,1) = xi_x_up(i,j)*x_xixi+xi_y_up(i,j)*y_xixi
				alpha(i,j,2) = 2.d0*(xi_x_up(i,j)*x_xiet+xi_y_up(i,j)*y_xiet)
				alpha(i,j,3) = xi_x_up(i,j)*x_etet+xi_y_up(i,j)*y_etet
			end do
		end do


! ---- coefs at v grid point ------
!
		do j=1,ny-1
			do i=1,nx
				xi_x_vp(i,j) = (xi_x(i,j)+xi_x(i,j+1))*.5d0
				et_x_vp(i,j) = (et_x(i,j)+et_x(i,j+1))*.5d0
				xi_y_vp(i,j) = (xi_y(i,j)+xi_y(i,j+1))*.5d0
				et_y_vp(i,j) = (et_y(i,j)+et_y(i,j+1))*.5d0
				beta(i,j,3) = et_x_vp(i,j)*xi_x_vp(i,j)+et_y_vp(i,j)*xi_y_vp(i,j)
				beta(i,j,4) = et_y_vp(i,j)**2d0+et_x_vp(i,j)**2d0
				
				if( i==1 ) then
					x_xixi = (x_xi(i+1,j+1)+x_xi(i+1,j)-x_xi(i,j+1)-x_xi(i,j))/(2.d0*dxi)
					y_xixi = (y_xi(i+1,j+1)+y_xi(i+1,j)-y_xi(i,j+1)-y_xi(i,j))/(2.d0*dxi)
				else if( i==nx ) then
					x_xixi = (x_xi(i,j+1)+x_xi(i,j)-x_xi(i-1,j+1)-x_xi(i-1,j))/(2.d0*dxi)
					y_xixi = (y_xi(i,j+1)+y_xi(i,j)-y_xi(i-1,j+1)-y_xi(i-1,j))/(2.d0*dxi)
				else
					x_xixi = (x_xi(i+1,j+1)+x_xi(i+1,j)-x_xi(i-1,j+1)-x_xi(i-1,j))/(4.d0*dxi)
					y_xixi = (y_xi(i+1,j+1)+y_xi(i+1,j)-y_xi(i-1,j+1)-y_xi(i-1,j))/(4.d0*dxi)
				end if
				
				x_xiet = (x_xi(i,j+1)-x_xi(i,j))/det
				y_xiet = (y_xi(i,j+1)-y_xi(i,j))/det
				x_etet = (x_et(i,j+1)-x_et(i,j))/det
				y_etet = (y_et(i,j+1)-y_et(i,j))/det
				alpha(i,j,4) = et_x_vp(i,j)*x_xixi+et_y_vp(i,j)*y_xixi
				alpha(i,j,5) = 2.d0*(et_x_vp(i,j)*x_xiet+et_y_vp(i,j)*y_xiet)
				alpha(i,j,6) = et_x_vp(i,j)*x_etet+et_y_vp(i,j)*y_etet
			end do
		end do

		do i=1,nx
			alpha(i,0,4) = alpha(i,1,4)
			alpha(i,0,5) = alpha(i,1,5)
			alpha(i,0,6) = alpha(i,1,6)
			alpha(i,ny,4) = alpha(i,ny-1,4)
			alpha(i,ny,5) = alpha(i,ny-1,5)
			alpha(i,ny,6) = alpha(i,ny-1,6)
			beta(i,0,3) = beta(i,1,3)
			beta(i,0,4) = beta(i,1,4)
			beta(i,ny,3) = beta(i,ny-1,3)
			beta(i,ny,4) = beta(i,ny-1,4)
		end do

!
! ------------------ cos(theta) ----
!
		do j=1,ny
			do i=1,nx
				x1 = ( x(i-1,j  ) + x(i-1,j-1) ) * 0.5d0
				y1 = ( y(i-1,j  ) + y(i-1,j-1) ) * 0.5d0
				x2 = (x(i,j)+x(i,j-1))*.5d0
				y2 = (y(i,j)+y(i,j-1))*.5d0
				x3 = (x(i,j)+x(i-1,j))*.5d0
				y3 = (y(i,j)+y(i-1,j))*.5d0
				x4 = ( x(i  ,j-1) + x(i-1,j-1) ) * 0.5d0
				y4 = ( y(i  ,j-1) + y(i-1,j-1) ) * 0.5d0
				dx1 = x2-x1
				dy1 = y2-y1
				dx2 = x3 - x4
				dy2 = y3 - y4
				ds1 = dsqrt((x2-x1)**2d0+(y2-y1)**2d0)
				ds2 = dsqrt((x3-x4)**2d0+(y3-y4)**2d0)
				ds12 = ds1*ds2
				cos_t(i,j) = (dx1*dx2+dy1*dy2)/ds12
			end do
		end do
!
! ------------------ sin(theta) ----
!
		do j=1,ny
			do i=1,nx-1
				x1 = (x(i,j)+x(i,j-1)+x(i-1,j)+x(i-1,j-1))*.25
				x2 = (x(i,j)+x(i,j-1)+x(i+1,j)+x(i+1,j-1))*.25
				y1 = (y(i,j)+y(i,j-1)+y(i-1,j)+y(i-1,j-1))*.25
				y2 = (y(i,j)+y(i,j-1)+y(i+1,j)+y(i+1,j-1))*.25
				x3 = x(i,j-1)
				x4 = x(i,j)
				y3 = y(i,j-1)
				y4 = y(i,j)
				dx1 = x2-x1
				dy1 = y2-y1
				dx2 = x4-x3
				dy2 = y4-y3
				ds1 = sqrt((x2-x1)**2d0+(y2-y1)**2d0)
				ds2 = sqrt((x4-x3)**2d0+(y4-y3)**2d0)
				ds12 = ds1*ds2
				sin_t(i,j) = (dx1*dy2-dx2*dy1)/ds12
			end do
		end do

	end subroutine gcoefs
! --------------------------------------------------------------------

end module gcoefs_m		! -----------------------------------------------------------------

module initl_m		! --------------------------------------------------------------------
	use common_hh
	use variables
	implicit none

	double precision,dimension(:),allocatable :: h_a,s_center,a_a,al_a,sie
	double precision,dimension(:,:),allocatable :: uti

	contains

! -------------------------------------------------------------------
	subroutine alloc_initl_temp_variables
		implicit none

		allocate( h_a(0:nx), s_center(0:nx) )
		allocate( a_a(0:nx), al_a(0:nx), sie(0:nx) )
		allocate( uti(0:nx,0:ny) )

		h_a=0.; uti=0.; s_center=0.; a_a=0.; al_a=0.; sie=0.
		uti = 0.d0

	end subroutine alloc_initl_temp_variables
! --------------------------------------------------------------------
	subroutine initl(qp,hnx,snu_0,i_flow,h_slope,x_bk,h_slope_1,h_slope_2)
		implicit none

		integer :: i, j, i_flow, nc

		double precision :: qp, hnx, snu_0 &
					, h_slope, x_bk, h_slope_1, h_slope_2	&
					, qp_ttl, xfdse, dxf, hup, hdw, hwt, qd
		double precision :: qpp, hs00, u_hp, bs, qu00, hmin01, b1, b2, dx	&
					, h_a1, dy, hsup, hsvp, qdiff, s_dummy
!
		s_center(0) = 0.d0

		do i=1,nx
			s_center(i) = s_center(i-1)+ds(i,nym)
		end do
!
		h_a(nx) = hnx
!
		if( i_flow<=1 ) then
			xfdse = ds(nx,nym)*.5d0
			do i=nx-1,1,-1
				dxf = (ds(i+1,nym)+ds(i,nym))*.5d0
				xfdse = xfdse+dxf
				if( i_flow==0 ) then
					h_a(i) = h_a(i+1)+h_slope*dxf
				else if( xfdse<=x_bk ) then
					h_a(i) = xfdse*h_slope_1+hnx
				else
					h_a(i) = (xfdse-x_bk)*h_slope_2+x_bk*h_slope_1+hnx
				end if
			end do
		end if
!
		if( i_flow==4 ) then
			do j=1,ny
				do i=1,nx
					hs(i,j) = hmin
					h(i,j) = hs(i,j)+eta(i,j)
					hn(i,j) = h(i,j)
				end do
			end do
				
			if( j_side1==1 ) then
				do i=1,nx
					s_dummy = (eta(i,2)-eta(i,1))/dnx(i,1)
					eta(i,0) = eta(i,1)-s_dummy*dnx(i,0)
					hs(i,0) = hmin
					h(i,0) = hs(i,0)+eta(i,0)
					hn(i,0) = h(i,0)
				end do
			end if
				
			if( j_side2==1 ) then
				do i=1,nx
					s_dummy = (eta(i,ny)-eta(i,ny-1))/dnx(i,ny-1)
					eta(i,ny+1) = eta(i,ny)+s_dummy*dnx(i,ny)
					hs(i,ny+1) = hmin
					h(i,ny+1) = hs(i,ny+1)+eta(i,ny+1)
					hn(i,ny+1) = h(i,ny+1)
				end do
			end if
		else if( i_flow==5 ) then
			do j=1,ny
				do i=1,nx
					h(i,j) = hs(i,j)+eta(i,j)
					hn(i,j) = h(i,j)
				end do
			end do
				
			if( j_side1==1 ) then
				do i=1,nx
					s_dummy = (eta(i,2)-eta(i,1))/dnx(i,1)
					eta(i,0) = eta(i,1)-s_dummy*dnx(i,0)
					hs(i,0) = hmin
					h(i,0) = hs(i,0)+eta(i,0)
					hn(i,0) = h(i,0)
				end do
			end if
				
			if( j_side2==1 ) then
				do i=1,nx
					s_dummy = (eta(i,ny)-eta(i,ny-1))/dnx(i,ny-1)
					eta(i,ny+1) = eta(i,ny)+s_dummy*dnx(i,ny)
					hs(i,ny+1) = hmin
					h(i,ny+1) = hs(i,ny+1)+eta(i,ny+1)
					hn(i,ny+1) = h(i,ny+1)
				end do
			end if
		else
			do j=1,ny
				do i=1,nx
					hs(i,j) = h_a(i)-eta(i,j)
					if( hs(i,j)<=hmin ) hs(i,j) = hmin
					h(i,j) = eta(i,j)+hs(i,j)
					hn(i,j) = h(i,j)
				end do
			end do
!i
!i 初期水面形が一定勾配も利用できるように。
!i
			if( j_side1==1 ) then
				do i=1,nx
					s_dummy = (eta(i,2)-eta(i,1))/dnx(i,1)
					eta(i,0) = eta(i,1)-s_dummy*dnx(i,0)
					hs(i,0) = hs(i,1)
					h(i,0) = hs(i,0)+eta(i,0)
					hn(i,0) = h(i,0)
				end do
			end if
				
			if( j_side2==1 ) then
				do i=1,nx
					s_dummy = (eta(i,ny)-eta(i,ny-1))/dnx(i,ny-1)
					eta(i,ny+1) = eta(i,ny)+s_dummy*dnx(i,ny)
					hs(i,ny+1) = hs(i,ny)
					h(i,ny+1) = hs(i,ny+1)+eta(i,ny+1)
					hn(i,ny+1) = h(i,ny+1)
				end do
			end if
!i
		end if
!
			! ---- 氾濫計算ではslopeが求められないので無効にする ----

!		do i=1,nx-1
!			qc(i)=0.d0
!			do j=1,ny
!				hsup = (hs(i,j)+hs(i+1,j))*.5d0
!				uti(i,j) = hsup**(2.d0/3.d0)*sqrt(slope)/sn_up(i,j)
!				yu(i,j) = uti(i,j)*xi_r_up(i,j)
				
!				if( (hs(i,j)<=hmin.and.yu(i,j)>0.).or.	&
!					(hs(i+1,j)<=hmin.and.yu(i,j)<0.) ) then
!					uti(i,j) = 0.d0
!					yu(i,j) = 0.d0
!				end if
				
!				if( ijobst(i,j)==1.and.ijobst(i,j-1)==1 ) then
!					uti(i,j) = 0.d0
!					yu(i,j) = 0.d0
!				end if
				
!				q_xi(i,j) = yu(i,j)*(hs(i,j)+hs(i+1,j))/(sj(i,j)+sj(i+1,j))
!				qu(i,j) = yu(i,j)/xi_r_up(i,j)*(hs(i,j)+hs(i+1,j))*.5d0*dn(i,j)
!				qc(i) = qc(i)+qu(i,j)
!			end do
!		end do
!
			! ---- 氾濫計算では初期流速は0 ----

		do j=1-j_side1,ny+j_side2
			do i=1,nx-1
				yu(i,j) = 0.d0
			end do
		end do

		do j=1,ny
			do i=1,nx
				q_xi(i,j) = 0.d0
			end do
		end do

		do j=1,ny
			qu(0,j) = qu(1,j)
			q_xi(0,j) = q_xi(1,j)
			yu(0,j) = yu(1,j)
			qu(nx,j) = qu(nx-1,j)
			q_xi(nx,j) = q_xi(nx-1,j)
			yu(nx,j) = yu(nx-1,j)
		end do
!
		do j=1,ny
			do i=0,nx
				yun(i,j) = yu(i,j)
			end do
		end do
!
		do j=1,ny
			do i=1,nx
				usta(i,j) = 0
				snu(i,j) = snu_0
			end do
		end do

		do j=1,ny-1
			do i=0,nx
				snu_x(i,j) = snu_0
			end do
		end do
!
		yv=0.
		yvn=0.
		q_et=0.
!		yc=0.
!		ycn=0.
!
	end subroutine initl
! --------------------------------------------------------------------

end module initl_m		! ----------------------------------------------------------------


! --------------------------------------------------------------------
	subroutine taustacal(snu00)
		use common_hh
		use variables
		implicit none
		
		double precision,intent(in) :: snu00
		
		integer :: i,j
		double precision :: us_2
!
!$omp do private(us_2)
		do j=1,ny
			do i=1,nx
				vti(i,j) = dsqrt(ux(i,j)**2d0+uy(i,j)**2d0)

				if( abs(vti(i,j))<1e-8.or.hs(i,j)<hmin ) then
					usta(i,j) = 0.
				else
					cf(i,j) = g*snmm(i,j)**2d0/hs(i,j)**(0.3333333d0)
					us_2 = cf(i,j)*vti(i,j)**2d0
					usta(i,j) = sqrt(us_2)
					re(i,j) = vti(i,j)*hs(i,j)/snu00
				end if
			end do
		end do
!
	end subroutine taustacal
!
! --------------------------------------------------------------------
	subroutine hshift(hn,h)
		use common_hh
		implicit none

		double precision,dimension(1:nx,0:ny+1),intent(in) :: hn
		double precision,dimension(1:nx,0:ny+1),intent(out) :: h

		integer :: i,j
!
!$omp do
		do j=1-j_side1,ny+j_side2
			do i=1,nx-1
				h(i,j)=hn(i,j)
			end do
		end do

	end subroutine hshift
!
! --------------------------------------------------------------------
	subroutine snucal(snu00,a_snu,b_snu)	!inoue 渦動粘性係数の補正
		use common_hh
		use variables
		implicit none
		
		double precision,intent(in) :: snu00,a_snu,b_snu
		
		integer :: i,j
		
!		if(time.eq.0.) write(*,*) a_snu, b_snu
!
!$omp do
		do j=1,ny
			do i=1,nx
				if( ijo_in(i,j)==1.or.hs(i,j)<hmin ) then
					snu(i,j) = 0.
				else if( re(i,j)<=2000. ) then
					snu(i,j) = snu00
				else
					snu(i,j) = 0.4/6.*usta(i,j)*hs(i,j)*a_snu+b_snu
!					snu(i,j) = 0.2*usta(i,j)*hs(i,j)
				end if
			end do
		end do

!$omp do
		do i=1,nx
			snu(i,   0) = snu(i, 1)
			snu(i,ny+1) = snu(i,ny)
		end do

!$omp do
		do j=1-j_side1,ny-1+j_side2
			do i=1,nx-1
				snu_x(i,j)=(snu(i,j)+snu(i+1,j)+snu(i,j+1)+snu(i+1,j+1))*.25d0*(1.d0-ijobst(i,j))
			end do
		end do

	end subroutine snucal
!
! --------------------------------------------------------------------
	subroutine hqtcal_inn(nq,h_down)
		use common_hh
		use variables
		implicit none
		
		integer,intent(in) :: nq
		double precision,intent(in) :: h_down
		
		integer :: i,j,n,ii
		double precision :: hmax1,hmin1,hh,qcc,hs1,qs,as,u0
		double precision,dimension(j_in) :: qmax,b_ups,hsmax,emin2
		double precision,dimension(j_in,0:nq) :: qmin
!
		do i=1,j_in
			qmax(i) = -999.
		end do
!
		do i=1,j_in
			do n=0,nq
				qmax(i) = max(qmax(i),q_ups_in(i,n))
				qmin(i,n) = q_ups_in(i,n)/1000.
			end do
		end do
!
		emin2 = 99999.
!
		do i=1,j_in
			do j=js(i),je(i)
				emin2(i) = min(emin2(i),eta(1,j))
			end do
		end do
!
		do i=1,j_in
			b_ups(i) = 0.
			do j=js(i),je(i)
				if( ijobst(0,j-1)+ijobst(0,j)==0 ) then
					b_ups(i) = b_ups(i)+dn(0,j)
				else if( ijobst(0,j-1)+ijobst(0,j)==1 ) then
					b_ups(i) = b_ups(i)+dn(0,j)*.5
				end if
			end do
		end do
!
		do i=1,j_in
			hsmax(i) = (snmm(1,(js(i)+je(i))*0.5)*qmax(i)		&
					/(b_ups(i)*sqrt(upv_slope_in(i))))**(3./5.)*100.
		end do
!
! 上流端水位の計算
!
		do ii=1,j_in
			do n=0,nq
!
				if( q_ups_in(ii,n)<=1e-6 ) then
					h_ups_in(ii,n)=emin2(ii)
					goto 103
				end if

				hmax1 = eta(1,(js(ii)+je(ii))*0.5)+hsmax(ii)*10.
				hmin1 = emin2(ii)
			120	hh = (hmax1+hmin1)*.5
				qcc = 0.

				do j=js(ii),je(ii)
					hs1 = hh-eta(1,j)
					if( hs1>0. .and. ijo_in(1,j)==0 ) then
						as = hs1*dn(0,j)
						u0 = 1./snmm(1,j)*hs1**(2./3.)*dsqrt(upv_slope_in(ii))
						qs = as*u0
						qcc = qcc+qs
					end if
				end do

				if( abs(qcc-q_ups_in(ii,n))<qmin(ii,n) ) goto 100
				if( qcc>q_ups_in(ii,n) ) then
					hmax1 = hh
				else
					hmin1 = hh
				end if
				goto 120
			100	continue

				h_ups_in(ii,n) = hh
!
			103	continue
!
			end do
!
		end do
!i
		if( j_wl==0 ) then
			do n=0,nq
				h_dse(n) = h_down
			end do
		end if
!
	end subroutine hqtcal_inn
!
! --------------------------------------------------------------------
	subroutine hqtcal_j1(nq)
		use common_hh
		use variables
		implicit none
		
		integer,intent(in) :: nq
		
		integer :: i,j,n,ii
		double precision :: hmax1,hmin1,hh,qcc,hs1,qs,as,u0
		double precision,dimension(jsin1) :: qmax,b_ups,hsmax,emin2
		double precision,dimension(jsin1,0:nq) :: qmin
!
		do i=1,jsin1
			qmax(i)=-999.
		end do
!
		do i=1,jsin1
			do n=0,nq
				qmax(i) = max(qmax(i),q1_ups_in(i,n))
				qmin(i,n) = q1_ups_in(i,n)/1000.
			end do
		end do
!
		emin2 = 99999.
!
		do i=1,jsin1
			do j=j1s(i),j1e(i)
				emin2(i) = min(emin2(i),eta(j,1))
			end do
		end do
!
		do i=1,jsin1
			b_ups(i) = 0.
			do j=j1s(i),j1e(i)
				if( ijobst(j-1,0)+ijobst(j,0)==0 ) then
					b_ups(i) = b_ups(i)+ds(j,0)
				else if( ijobst(j-1,0)+ijobst(j,0)==1 ) then
					b_ups(i) = b_ups(i)+ds(j,0)*.5
				end if
			end do
		end do
!
		do i=1,jsin1
			hsmax(i) = (snmm((j1s(i)+j1e(i))*0.5,1)*qmax(i)		&
					/(b_ups(i)*sqrt(vpv1_slope_in(i))))**(3./5.)*100.
		end do
!
! 上流端水位の計算
!
		do ii=1,jsin1
			do n=0,nq
!
				if( q1_ups_in(ii,n)<=1e-6 ) then
					h1_ups_in(ii,n) = emin2(ii)
					goto 103
				end if

				hmax1 = eta((j1s(ii)+j1e(ii))*0.5,1)+hsmax(ii)*10.
				hmin1 = emin2(ii)
			120	hh = (hmax1+hmin1)*.5
				qcc = 0.

				do j=j1s(ii),j1e(ii)
					hs1 = hh-eta(j,1)
					if( hs1>0. .and. ijo_in(j,1)==0 ) then
						as = hs1*ds(j,0)
						u0 = 1./snmm(j,1)*hs1**(2./3.)*sqrt(vpv1_slope_in(ii))
						qs = as*u0
						qcc = qcc+qs
					end if
				end do

				if( abs(qcc-q1_ups_in(ii,n))<qmin(ii,n) ) goto 100
				if( qcc>q1_ups_in(ii,n) ) then
					hmax1 = hh
				else
					hmin1 = hh
				end if

				goto 120
			100	continue

				h1_ups_in(ii,n) = hh
!
			103	continue
!
			end do
!
		end do

	end subroutine hqtcal_j1
!
! --------------------------------------------------------------------
	subroutine hqtcal_j2(nq)
		use common_hh
		use variables
		implicit none
		
		integer,intent(in) :: nq
		
		integer :: i,j,n,ii
		double precision :: hmax1,hmin1,hh,qcc,hs1,qs,as,u0
		double precision,dimension(jsin2) :: qmax,b_ups,hsmax,emin2
		double precision,dimension(jsin2,0:nq) :: qmin
!
		do i=1,jsin2
			qmax(i) = -999.
		end do
!
		do i=1,jsin2
			do n=0,nq
				qmax(i) = max(qmax(i),q2_ups_in(i,n))
				qmin(i,n) = q2_ups_in(i,n)/1000.
			end do
		end do
!
		emin2 = 99999.
!
		do i=1,jsin2
			do j=j2s(i),j2e(i)
				emin2(i) = min(emin2(i),eta(j,ny))
			end do
		end do
!
		do i=1,jsin2
			b_ups(i) = 0.
			do j=j2s(i),j2e(i)
				if( ijobst(j-1,ny)+ijobst(j,ny)==0 ) then
					b_ups(i) = b_ups(i)+ds(j,ny)
				else if( ijobst(j-1,ny)+ijobst(j,ny)==1 ) then
					b_ups(i) = b_ups(i)+ds(j,ny)*.5
				end if
			end do
		end do
!
		do i=1,jsin2
			hsmax(i) = (snmm((j2s(i)+j2e(i))*0.5,ny)*qmax(i)		&
					/(b_ups(i)*sqrt(vpv2_slope_in(i))))**(3./5.)*100.
		end do
!
! 上流端水位の計算
!
		do ii=1,jsin2
			do n=0,nq
!
				if( q2_ups_in(ii,n)<=1e-6 ) then
					h2_ups_in(ii,n) = emin2(ii)
					goto 103
				end if
				
				hmax1 = eta((j2s(ii)+j2e(ii))*0.5,ny)+hsmax(ii)*10.
				hmin1 = emin2(ii)
			120	hh = (hmax1+hmin1)*.5
				qcc = 0.

				do j=j2s(ii),j2e(ii)
					hs1 = hh-eta(j,ny)
					if( hs1>0. .and. ijo_in(j,ny)==0 ) then
						as = hs1*ds(j,ny)
						u0 = 1./snmm(j,ny)*hs1**(2./3.)*sqrt(vpv2_slope_in(ii))
						qs = as*u0
						qcc = qcc+qs
					end if
				end do

				if( abs(qcc-q2_ups_in(ii,n))<qmin(ii,n) ) goto 100
				if( qcc>q2_ups_in(ii,n) ) then
					hmax1 = hh
				else
					hmin1 = hh
				end if

				goto 120
			100	continue

				h2_ups_in(ii,n) = hh
!
			103	continue
!
			end do
!
		end do

	end subroutine hqtcal_j2

! --------------------------------------------------------------------
!
module upstream_m		! ----------------------------------------------------------------
	use common_hh
	use variables
	implicit none

	double precision,dimension(:),allocatable :: uinti,hsin

	double precision,dimension(:,:),allocatable :: uinti_inn,hsin_inn
	double precision,dimension(:),allocatable :: qc0_inn

	double precision,dimension(:,:),allocatable :: uinti_j1,hsin_j1
	double precision,dimension(:),allocatable :: qc0_j1

	double precision,dimension(:,:),allocatable :: uinti_j2,hsin_j2
 	double precision,dimension(:),allocatable :: qc0_j2

	contains

! --------------------------------------------------------------------
	subroutine alloc_upstream_temp_variables
		implicit none

		allocate( uinti(0:ny),hsin(0:ny) )
		allocate( uinti_inn(5,0:ny),hsin_inn(5,0:ny) )
		allocate( qc0_inn(5) )
		allocate( uinti_j1(5,0:nx),hsin_j1(5,0:nx) )
		allocate( qc0_j1(5) )
		allocate( uinti_j2(5,0:nx),hsin_j2(5,0:nx) )
		allocate( qc0_j2(5) )

		uinti=0.; hsin=0.
		uinti_inn=0.; hsin_inn=0.
		uinti_j1=0.; hsin_j1=0.
		uinti_j2=0.; hsin_j2=0.

	end subroutine alloc_upstream_temp_variables

! ---------------------------------------------------------------------
	subroutine upstream_inn()
		implicit none
		integer(4) :: i,j
		double precision :: as,qdiff
!
		do j=1,ny
			qu(0,j) = 0.
			yu(0,j) = 0.
			yun(0,j) = 0.
			q_xi(0,j) = 0.
		end do
!
		do i=1,j_in
			if(q_input_in(i) > 0.) then
				qc0_inn(i)=0.
				do j=js(i),je(i)
					hsin_inn(i,j)=h_input_in(i)-eta(1,j)
					if(hsin_inn(i,j) < 0.) hsin_inn(i,j)=0.
					uinti_inn(i,j)=1./snmm(1,j)*hsin_inn(i,j)**(2./3.)		&
									*sqrt(upv_slope_in(i))
					if(ijo_in(1,j) == 0) then
						as=hsin_inn(i,j)*dn(0,j)
					else
						as=0.
					end if
					
					qu(0,j)=uinti_inn(i,j)*as
					qc0_inn(i)=qc0_inn(i)+qu(0,j)
				end do

				qdiff=qc0_inn(i)/q_input_in(i)
				qc(0)=0.

				do j=js(i),je(i)
					if(ijo_in(1,j) == 0) then
						yu(0,j)=uinti_inn(i,j)/qdiff*xi_r_up(0,j)
					else
						yu(0,j)=0.
					end if

					qu(0,j)=qu(0,j)/qdiff
					q_xi(0,j)=yu(0,j)*hsin_inn(i,j)/sj(1,j)
					yun(0,j)=yu(0,j)
					qc(0)=qc(0)+qu(0,j)
				end do
			end if
!
		end do

	end subroutine upstream_inn
!
! ---------------------------------------------------------------------
	subroutine upstream_j1()
		implicit none
		integer(4) :: i,j
		double precision :: as,qdiff
!
		do i=1,nx
			qv(i,0)=0.
			yv(i,0)=0.
			yvn(i,0)=0.
			q_et(i,0)=0.
		end do

		do i=1,jsin1
			if(q_in_j1(i) > 0.) then
				qc0_j1(i)=0.
				do j=j1s(i),j1e(i)
					hsin_j1(i,j)=h_in_j1(i)-eta(j,1)
					if(hsin_j1(i,j) < 0.) hsin_j1(i,j)=0.
					uinti_j1(i,j)=1./snmm(j,1)*hsin_j1(i,j)**(2./3.)	&
									*sqrt(vpv1_slope_in(i))
									
					if(ijo_in(j,1) == 0) then
						as=hsin_j1(i,j)*ds(j,0)
					else
						as=0.
					end if
					
					qv(j,0)=uinti_j1(i,j)*as
					qc0_j1(i)=qc0_j1(i)+qv(j,0)
				end do

				qdiff=qc0_j1(i)/q_in_j1(i)
				qc(0)=0.

				do j=j1s(i),j1e(i)
					if(ijo_in(j,1) == 0) then
						yv(j,0)=uinti_j1(i,j)/qdiff*et_r_vp(j,0)
					else
						yv(j,0)=0.
					end if
					
					qv(j,0)=qv(j,0)/qdiff
					q_et(j,0)=yv(j,0)*hsin_j1(i,j)/sj(j,1)
					yvn(j,0)=yv(j,0)
					qc(0)=qc(0)+qv(j,0)
				end do
			end if
!
		end do

	end subroutine upstream_j1

! ---------------------------------------------------------------------
	subroutine upstream_j2()
		implicit none
		integer(4) :: i,j
		double precision :: as,qdiff
!
		do i=1,nx
			qv(i,ny)=0.
			yv(i,ny)=0.
			yvn(i,ny)=0.
			q_et(i,ny)=0.
		end do

		do i=1,jsin2
			if(q_in_j2(i) > 0.) then
				qc0_j2(i)=0.
				do j=j2s(i),j2e(i)
					hsin_j2(i,j)=h_in_j2(i)-eta(j,ny)
					if(hsin_j2(i,j) < 0.) hsin_j2(i,j)=0.
					uinti_j2(i,j)=1./snmm(j,ny)*hsin_j2(i,j)**(2./3.)	&
									*sqrt(vpv2_slope_in(i))

					if(ijo_in(j,ny) == 0) then
						as=hsin_j2(i,j)*ds(j,ny)
					else
						as=0.
					end if

					qv(j,ny)=uinti_j2(i,j)*as
					qc0_j2(i)=qc0_j2(i)+qv(j,ny)
				end do

				qdiff=qc0_j2(i)/q_in_j2(i)
				qc(0)=0.

				do j=j2s(i),j2e(i)
					if(ijo_in(j,ny) == 0) then
						yv(j,ny)=-uinti_j2(i,j)/qdiff*et_r_vp(j,ny)
					else
						yv(j,ny)=0.
					end if

					qv(j,ny)=qv(j,ny)/qdiff
					q_et(j,ny)=yv(j,ny)*hsin_j2(i,j)/sj(j,ny)
					yvn(j,ny)=yv(j,ny)
					qc(0)=qc(0)+qv(j,ny)
				end do
			end if
!
		end do

	end subroutine upstream_j2
! --------------------------------------------------------------------

end module upstream_m		! ------------------------------------------------------------

! --------------------------------------------------------------------
	subroutine downstream(hnx)
		use common_hh
		use variables
		implicit none
		
		double precision,intent(in) :: hnx
		
		integer :: j

		do j=1,ny
			if( ijo_in(nx,j)==0 ) then
				hs(nx,j) = hnx-eta(nx,j)
				if( hs(nx,j)<hmin ) hs(nx,j) = hmin
				h(nx,j) = eta(nx,j)+hs(nx,j)
				hn(nx,j) = h(nx,j)
			end if
		end do

	end subroutine downstream
!
! --------------------------------------------------------------------
	subroutine wfcal(s,d,snu,wf,g)
		implicit real*8 (a-h,o-z)

		if(d>0.001) then
			wf=32.8*sqrt(d*100.)*0.01
		else
			wf=sqrt(2./3.*s*g*d+(6.*snu/d)**2)-6.*snu/d
		end if

	end subroutine wfcal

! ----------------------------------------------------------------------------------------
!                           end module and subroutine
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
!                                   main program
! ----------------------------------------------------------------------------------------
!
program nays2d_flood_parallel

	use common_hh
	use variables
	use output_cgn
	use GridCoord
	use GridCond
	use avgeo_m
	use gcoefs_m
	use initl_m
	use bound
	use non_advec_1
	use cal_gradient
	use diffusion_m
	use advection_m
	use upstream_m
	use allocate_variables
	use interpolation
	use transform
!docon add start#####################################################################################
	use box_gate_pump   !140602
!docon add end#######################################################################################
	implicit none

	include "cgnslib_f.h" 
	include "iriclib_f.h"

! --------------------------------------------------------------------
	integer,dimension(:,:),allocatable :: nb_io
	double precision   , dimension(:,:),allocatable :: h_max,v_max
	double precision         :: time8, v_temp

! --------------------------------------------------------------------

	integer(4), dimension(:,:), allocatable :: obst4, vege4, fm4, ob, q_cell, ibc0
	double precision   , dimension(:,:), allocatable :: chcond4
	double precision   , dimension(:,:), allocatable :: x8, y8, z8, hs8
	integer(4) ni4,nj4,iobst4
	double precision hmin8

	integer(4) :: ioutflag

!  for cngs i/0
!
	INTEGER :: FID, IER
	INTEGER :: index = 1
	double precision :: qptemp
	real, DIMENSION(:), ALLOCATABLE :: xtmp, ytmp
	double precision, DIMENSION(:), ALLOCATABLE :: xtmp_d, ytmp_d
	INTEGER :: tmpint

	integer :: i, j, k, n
	integer :: icount, istatus, ierr, i_out
	integer :: j_slope, j_upv, j_upv_slope, i_flow, j_cip	&
				, j_snu, j_inn, indexmax, jsizemax, i_sec_hour, nq, kend	&
				, kmod, ndry, lcount, iofrg, ndeposit, n_parallel,total_bomb, qtmpin, nqtcell
	double precision :: skp, snu00, h_down, bh_slope, upv_slope, h_slope	&
						, x_bk, h_slope_1, h_slope_2, e_thic, tantc	&
						, t_out_start, alh, slambda, rho, a_snu, b_snu	&	!inoue 渦動粘性係数の補正
						, snst, r_tantc, t_xx, thstart, qp, etime	&
						, d10, d50, d90, dm0, sn_g, cw, sigma_c, calculated_slope	&
						, slope, hs_dse	&
						, ts0, errmax, usci, theta_b, tan_tb, err, q_inptotal	&
						, dm_ave, cos_tb, gamma, gamma_m, rsgd, dtanmax, snu_0, hnx	&
						, wf, sst, sum_area
	double precision,dimension(:),allocatable :: xqtmp,yqtmp
	double precision,dimension(:),allocatable :: t_qc_dis
	double precision,dimension(:),allocatable :: q_tmp_dis, ni_dis_cell
	double precision,dimension(:,:),allocatable :: q_cell_dis

!inoue 建物阻害率 hamaki ver
	double precision,dimension(:,:),allocatable :: sh4
	integer :: j_gam
!inoue

	character(50) :: dis_label, mix_label, cm
	
	! For Hot Stat
  integer :: is, ii, iii, i_re_flag_i, i_re_flag_o, n_rest, i_tmp_count
  real(8) :: opt_tmp(0:9)
  real(8) :: ttt
  character(len = strMax) :: tmp_file_o(0:9), tmp_caption(0:9) &
       ,tmp_file_i, tmp_dummy, tmp_pass
		 
  i_tmp_count = 0
  do ii=0,9
     tmp_dummy = "tmp0.d"
     write(tmp_dummy(4:4),'(i1)') ii
     tmp_file_o(ii) = tmp_dummy
     tmp_dummy = "outtime_tmp0"
     write(tmp_dummy(12:12),'(i1)') ii
     tmp_caption(ii) = tmp_dummy
  end do
	
	
	write(*,*) 'Nays2d_flood Solver Version 5.0.0000 Last updated 2014/5/14'
	write(*,*) 'Copyright(C) by Yasuyuki Shimizu, Hokkaido Univ., Japan'
	write(*,*) 'Modified by Ichiro Kimura, Toshiki Iwasaki, Satomi Kawamura, Takuya Inoue , Michihiro Hamaki , Takeshi Takemura'

	call system_clock(cal_t1)	!h160105 計算開始時時刻

	icount = nargs() 	!for intelfortran
!	icount = iargc()	!for gfortran
	
	i_out = 1

	if ( icount==2 ) then
		CALL GETARG (1, condFile, istatus)	!for intelfortran
!	if ( icount==1 ) then					!for gfortran
!		CALL GETARG (1, condFile)			!for gfortran
	else
		write( 6, * ) 'input filename: ?'
		read( 5, *) condFile
	endif

	g=9.8
	skp=0.4
	snu00=1e-6
	pi=3.14159265

	ioutflag = 0

! ---- open CGNS file and initialize iRIC -----

	CALL cg_open_f(trim(condFile), CG_MODE_MODIFY, FID, IER)

	CALL CG_IRIC_INIT_F(FID, IER)
	IF(IER .NE. 0) THEN
		call cg_error_print_f()
	ENDIF

   !guiにcgnsファイルを読込みであることを知らせるファイルを生成
   call iric_initoption_f(IRIC_OPTION_CANCEL, ier)

! ---- Read computational condition ----


	! ---- Parameters for Downstream Water Surface Elevation -----

	CALL CG_IRIC_READ_INTEGER_F('j_wl', j_wl, ier) !i
!
!			j_wl = 0 ...下流端水位一定値を与える(h_down)
!			j_wl = 1 ...下流端水位は等流計算で求める⇒!iはん濫モデルでは等流起算なし
!			j_wl = 2 ...下流端水位はファイルから読み込む
!			j_wl = 3 ...下流端の水深は，前のグリットのコピーにする


	CALL CG_IRIC_READ_REAL_F('h_down', h_down, ier)
!
!			h_dwown = 下流端水位の値(上記j_wl=0の時のみ有効)
!

	CALL CG_IRIC_READ_INTEGER_F('j_rain', j_rain, ier) !i
!
!			j_rain =1	降雨なし
!			j_rain =2	降雨あり(全領域同じ時系列データ)
!			j_rain =3	降雨あり(Xrain対応)
!

	j_slope=0  		!iはん濫モデルでは等流起算なし

	bh_slope=0.  	!iはん濫モデルでは等流起算なし
!
! ------ Parameters for Upstream Boundary ------

	j_upv = 1   	 !iはん濫モデルでは使わない。ダミーデータ
!
!			j_upv =1 上流端の流速を等流計算で与える
!			j_upv =2 上流端の流速を、上流端の水深を使って流量から逆算する
!
	j_upv_slope=1      !iはん濫モデルでは使わない。ダミーデータ
!
!			上記j_upv=1のときの等流計算に使用する勾配の与え方
!			j_upv_slope=0 .... 河床データから自動的に計算
!			j_upv_slope=1 .... 値を与える→この場合は次の項目のuvp_slope
!
	upv_slope=0.001      !iはん濫モデルでは使わない。ダミーデータ

! ---- Parameters for Initial Water Surface Profile-----
!
	CALL CG_IRIC_READ_INTEGER_F('i_flow', i_flow, ier)
!
!			i_flow=0 初期水面形は直線(一定勾配)
!			i_flow=1 初期水面形は折線(１折点と２直線) !i　はん濫計算モデルでは使用不可
!			i_flow=2 初期水面形は等流計算             !i　はん濫計算モデルでは使用不可
!			i_flow=3 初期水面形は不等流計算           !i　はん濫計算モデルでは使用不可
!			i_flow=4 初期水面形は0

	CALL CG_IRIC_READ_REAL_F('h_slope', h_slope, ier)
!
!			上記i_flow=0のときの初期水面勾配
!
!i			CALL CG_IRIC_READ_REALSINGLE_F('x_bk', x_bk, ier)

	x_bk = 5  !i　はん濫計算モデルでは使用不可　ダミー。
!
!			上記i_flow=1のときの勾配変化点の下流からの距離 x_bk
!
!		 CALL CG_IRIC_READ_REALSINGLE_F('h_slope_1', h_slope_1, ier)
!		CALL CG_IRIC_READ_REALSINGLE_F('h_slope_2', h_slope_2, ier)

	h_slope_1 =0.001  !i　はん濫計算モデルでは使用不可　ダミー。
	h_slope_2 =0.001  !i　はん濫計算モデルでは使用不可　ダミー。

!
!			上記i_flow=1のときの初期水面勾配(下流側)h_slope_1
!			上記i_flow=1のときの初期水面勾配(上流側)h_slope_2

!
!	側方自由流出パラメータ
!
! ----- Parameters for innundation with some inlet ----- !
!

	i_in=1      !iはん濫計算だけを行うモデル。

	if(i_in == 1) then
		j_upv = 1
		j_upv_slope = 1

		allocate(upv_slope_in(50))
		allocate(js(50),je(50))
		allocate(q_input_in(50),h_input_in(50))
	end if

	CALL CG_IRIC_READ_INTEGER_F('j_side_j1', j_side_j1, ier)

!			j_side_j1 : flag for in/out flow at j=1

	if(j_side_j1 <=1) then
		j_side1=j_side_j1
	else
		j_side1=0
	end if

	allocate(vpv1_slope_in(50))
	allocate(j1s(50),j1e(50))
	allocate(q_in_j1(50),h_in_j1(50))

	CALL CG_IRIC_READ_INTEGER_F('j_side_j2', j_side_j2, ier)

!			jside_j2 : flag for in/out flow at j=ny

	if(j_side_j2 <=1) then
		j_side2=j_side_j2
	else
		j_side2=0
	end if

	allocate(vpv2_slope_in(50))
	allocate(j2s(50),j2e(50))
	allocate(q_in_j2(50),h_in_j2(50))
!
! -----Parameters on Bed Roughness -----
!
!			CALL CG_IRIC_READ_INTEGER_F('j_drg', j_drg, ier)
!
!			j_drg = 0 .... 粗度は河床材料(diam)から自動的に計算される
!			j_drg = 1 .... 粗度は与えられる（以下、具体的粗度の値)

	j_drg = 1

! ---- Parameters on Time Setting -----
!
	CALL CG_IRIC_READ_REAL_F('tuk', tuk, ier)
	CALL CG_IRIC_READ_REAL_F('dt', dt, ier)
	CALL CG_IRIC_READ_REAL_F('t_out_start', t_out_start, ier)
	CALL CG_IRIC_READ_REAL_F('bomber_time', bomber_time, ier)
	
	CALL CG_IRIC_READ_REAL_F('ster', ster, ier)
	ster = 0.

	CALL CG_IRIC_READ_INTEGER_F('j_cip', j_cip, ier)

!				CALL CG_IRIC_READ_INTEGER_F('j_snu', j_snu, ier)
	j_snu=1  		!i Zero equation model

 !
   ! --- Parameters for Hot Start ---
   !
     CALL CG_IRIC_READ_INTEGER_F('write_flag', i_re_flag_o, ier)
     CALL CG_IRIC_READ_INTEGER_F('read_flag', i_re_flag_i, ier)
     CALL CG_IRIC_READ_INTEGER_F('n_tempfile', n_rest, ier)
     CALL CG_IRIC_READ_STRING_F('tmp_readfile', tmp_file_i, ier)
     CALL CG_IRIC_READ_STRING_F('tmp_pass', tmp_pass, ier)
     
     !暫定処理
     i_re_flag_o = 0
     i_re_flag_i = 0
     n_rest = 1

     do ii=0,9
        opt_tmp(ii) = -9999
        CALL CG_IRIC_READ_REAL_F(tmp_caption(ii), opt_tmp(ii), ier)
     end do
     !
     if(n_rest > 1) then
         do iii=1,n_rest
            do ii=0,n_rest-1
               if(opt_tmp(ii) /= opt_tmp(ii+1) &
                    .or.opt_tmp(ii+1) < opt_tmp(ii)+dt) then
                  if(opt_tmp(ii) > opt_tmp(ii+1)) then
                     ttt = opt_tmp(ii)
                     opt_tmp(ii) = opt_tmp(ii+1)
                     opt_tmp(ii+1) = ttt
                  end if
               else
                  opt_tmp(ii+1) = opt_tmp(ii+1)+dt
               end if
            end do
         end do
    end if
!
! --- Other Parameters and Constants -----
!
	CALL CG_IRIC_READ_INTEGER_F('lmax', lmax, ier)
	CALL CG_IRIC_READ_REAL_F('alh', alh, ier)
	CALL CG_IRIC_READ_REAL_F('hmin', hmin, ier)
	CALL CG_IRIC_READ_REAL_F('a_snu', a_snu, ier)
	CALL CG_IRIC_READ_REAL_F('b_snu', b_snu, ier)

! ---- Parameter for parallel computation ----

	CALL CG_IRIC_READ_INTEGER_F('n_parallel', n_parallel, ier)
!docon add start######################################################################################
	if(n_parallel<1) n_parallel=1
!docon add start######################################################################################

! ---- Parameter for Method of Share to gamma ---- !h131029

	CALL CG_IRIC_READ_INTEGER_F('j_qbl', j_qbl, ier)
	CALL CG_IRIC_READ_INTEGER_F('j_gam', j_gam, ier)
	CALL CG_IRIC_READ_REAL_F('cdd', cdd, ier)

! ---- Parameter for discharge cells ----

	CALL CG_IRIC_READ_INTEGER_F('n_q_cell', n_q_cell, ier)

	if( n_q_cell>0 ) then
		
		CALL CG_IRIC_READ_FUNCTIONALSIZE_F('t_qcells',qtmpin,ier)
		
		nqtcell = qtmpin-1
		
		allocate( xqtmp(qtmpin), yqtmp(qtmpin) )
		allocate( q_tmp_dis(0:n_q_cell) )
		allocate( t_qc_dis(0:nqtcell), q_cell_dis(0:nqtcell,n_q_cell)	&
				, ni_dis_cell(0:n_q_cell) )
		
		q_tmp_dis = 0.d0	! 配列0は流入がないとき用
		ni_dis_cell = 1.d0
		
		CALL CG_IRIC_READ_FUNCTIONALWITHNAME_F('t_qcells','time',xqtmp,ier)
		
		if( i_sec_hour==1 ) then
			do k=0,nqtcell
				t_qc_dis(k) = xqtmp(k+1)
			end do
		else
			do k=0,nqtcell
				t_qc_dis(k) = xqtmp(k+1)*3600.d0
			end do
		end if
	
		do k=1,n_q_cell
			write(cm,'(i1)') k
			dis_label = 'discharge'//trim(cm)
			
			CALL CG_IRIC_READ_FUNCTIONALWITHNAME_F('t_qcells',dis_label,yqtmp,ier)
			
			do n=0,nqtcell
				q_cell_dis(n,k) = yqtmp(n+1)
			end do
			
		end do
		
		deallocate( xqtmp, yqtmp )
		
	end if

! ------ Read initial bed and cell condition -----

	CALL CGNS_READ_GRIDCOORD(FID, ni4, nj4, x8, y8, IER)

	CALL CGNS_Read_GridCondition(FID, ni4, nj4, z8, hs8		&
								, obst4, vege4, chcond4, fm4, ob, q_cell, sh4, ibc0, IER)

!	CALL CG_IRIC_GOTOCC_F(FID, IER)
!

!i3------------------------------------------------------------
!i3 境界条件設定機能への対応
!i3
!
	j_in = 0
	jsin1 = 0
	jsin2 = 0

	call CGNS_Read_BoundaryCondition		&
			( j_inn,j_inlen,indices, indexmax, j_size, jsizemax, f_param, f_value, slopevalue,flowname )

	do i = 1, j_inn
		do j = 1, j_inlen(i)-1
			if(indices(i,1,1)==1 .and. indices(i,1,j_inlen(i))==1) then
				if(j==1) j_in=j_in+1
				upv_slope_in(j_in)=slopevalue(i)
				js(j_in) = indices(i,2,1)
				je(j_in) = indices(i,2,j_inlen(i)-1)
				write(*,*) 0,i,js(j_in),je(j_in)
			else if(indices(i,2,j)==1) then
				if(j_side_j1 == 2) then
					if(j==1) jsin1=jsin1+1
					vpv1_slope_in(jsin1)=slopevalue(i)
					j1s(jsin1) = indices(i,1,1)
					j1e(jsin1) = indices(i,1,j_inlen(i)-1)
				end if
				write(*,*) 1,i,j1s(jsin1),j1e(jsin1)
			else if(indices(i,2,j)==nj4) then
				if(j_side_j2 == 2) then
					if(j==1) jsin2=jsin2+1
					vpv2_slope_in(jsin2)=slopevalue(i)
					j2s(jsin2) = indices(i,1,1)
					j2e(jsin2) = indices(i,1,j_inlen(i)-1)
				end if
				write(*,*) 2,i,j2s(jsin2),j2e(jsin2)
			else
				write(*,*) 'Inflow point is not on boundary grid'
				stop
			end if
		end do
	end do
!
	if(j_in==0) then
		write(*,*) 'Inflow point is zero'
		stop
	end if

	if(j_in>50) then
		write(*,*) 'Inflow point is up to fifty on each side'
		stop
	end if
!
	if(j_side_j1 == 2) then
		if(jsin1>50) then
			write(*,*) 'Inflow point is up to fifty on each side'
			stop
		end if
	end if
!
	if(j_side_j2 == 2) then
		if(jsin2>50) then
			write(*,*) 'Inflow point is up to fifty on each side'
			stop
		end if
	end if
!
!i3--------- re-allocate array for subroutine -----------
!
	nx = ni4-1
	ny = nj4-1

	call alloc_var1
	call alloc_gcoefs_temp_variables
	call alloc_initl_temp_variables
	call alloc_diffusion_temp_variables
	call alloc_advection_temp_variables
	call alloc_upstream_temp_variables
	call alloc_uxxyycal_temp_variables
!
!-------------------------------------------------------------
!
! ---- set coordinate system -----
!
	do j=1,nj4
		do i=1,ni4
			x(i-1,j-1) = x8(i,j)
			y(i-1,j-1) = y8(i,j)
			z(i-1,j-1) = z8(i,j)
			z0(i-1,j-1) = z8(i,j)
		end do
	end do

	if( i_flow==5 ) then
		do j=1,ny
			do i=1,nx
				hs(i,j) = (hs8(i,j)+hs8(i+1,j)+hs8(i,j+1)+hs8(i+1,j+1))*0.25d0
				if( ijo_in(i,j)==1 ) hs(i,j) = hmin
			end do
		end do
	end if
!
! ------ set cell status -----
!
	do i=1,nx
		do j=1-j_side1,ny+j_side2

				! ----- set obstacle cells -----

			if( j==0 ) then
				ijo_in(i,j) = obst4(i,1)
                ibc(i,j) = 1
			elseif( j==ny+1 ) then
				ijo_in(i,j) = obst4(i,ny)
                ibc(i,j) = 1
			else
				ijo_in(i,j) = obst4(i,j)
				ibc(i,j) = ibc0(i,j)
				if( ob(i,j)==1 ) then
					ijo_in(i,j) = 1
				end if
			end if

				! ----- set fixed bed cells -----

			if( j==0 ) then
				ij_ero(i,j) = fm4(i,1)
			else if( j==ny+1 ) then
				ij_ero(i,j) = fm4(i,ny)
			else
				ij_ero(i,j) = fm4(i,j)
			end if

				! ----- set vegetation condition -----

!			if(    vege4(i,j) == 0 ) then
!				cd_veg(i,j) = 0.
!			elseif(    vege4(i,j) == 1 ) then
!				cd_veg(i,j) = veg_lamb_1 * c_tree * 0.5
!			elseif(    vege4(i,j) == 2 ) then
!				cd_veg(i,j) = veg_lamb_2 * c_tree * 0.5
!			elseif(    vege4(i,j) == 3 ) then
!				cd_veg(i,j) = veg_lamb_3 * c_tree * 0.5
!			end if

			cd_veg(i,j) = 0. 		!i　はん濫計算で樹木なし 

            
		end do
	end do
!
	ijobst=0

	do j=1,ny
		do i=1,nx
			if( ijo_in(i,j)==1 ) then
				lcb(i,j) = 0
				ijobst(i,j) = 1
				ijobst(i-1,j) = 1
				ijobst(i,j-1) = 1
				ijobst(i-1,j-1) = 1
			end if
		end do
	end do
	
		! --- 流速定義点の構造物判定 ---

	do j=1,ny
		do i=0,nx
			if( ijobst(i,j)==1 .and. ijobst(i,j-1)==1 ) then
				ijobst_u(i,j) = 0
			else
				ijobst_u(i,j) = 1
			end if
		end do
	end do

	do j=0,ny
		do i=1,nx
			if( ijobst(i,j)==1 .and. ijobst(i-1,j)==1 ) then
				ijobst_v(i,j) = 0
			else
				ijobst_v(i,j) = 1
			end if
		end do
	end do
	
	if( n_q_cell>0.and.n_q_cell<5 ) then
		do k=n_q_cell+1,5
			do j=1,ny
				do i=1,nx
					if( q_cell(i,j)==k ) then
						q_cell(i,j) = 0
					end if
				end do
			end do
		end do
	end if
	
! -------- Setup for output file -------
!
	allocate( h_max(0:nx,0:ny),v_max(0:nx,0:ny) )
!
	do j=0,ny
		do i=0,nx
			h_max(i,j) = -999.d0
			v_max(i,j) = -999.d0
		end do
	end do
!
! ----------------------------------------

	jt1 = 0
	jt2 = ny
	nym = ny*0.5

	dxi = 1.d0/dble(nx)
	det = 1.d0/dble(ny)
	r_dxi = 1.d0/dxi
	r_det = 1.d0/det

	dxi1 = dxi
	dxi2 = dxi1*dxi1
	dxi3 = dxi2*dxi1
	det1 = det
	det2 = det1*det1
	det3 = det2*det1

! ---------- Read Discharge Data -----------

	CALL CG_IRIC_READ_INTEGER_F('i_sec_hour', i_sec_hour, ier)

	nq = jsizemax-1

	allocate( t_hyd(0:nq),q_ups(0:nq),h_dse(0:nq),h_ups(0:nq),rain(0:nq) )

	if( i_in/=0 ) then
		allocate( q_ups_in(j_in,0:nq),h_ups_in(j_in,0:nq) )
		q_ups_in = 0.;  h_ups_in = 0.
		allocate(fname_in(j_in))
	end if

	t_hyd=0.; q_ups=0.; h_dse=0.; h_ups=0.;rain=0.

	if( j_side_j1==2 ) then
		allocate( q1_ups_in(jsin1,0:nq),h1_ups_in(jsin1,0:nq) )
		q1_ups_in = 0.;  h1_ups_in = 0.
		allocate(fname_1(jsin1))
	end if

	if( j_side_j2==2 ) then
		allocate( q2_ups_in(jsin2,0:nq),h2_ups_in(jsin2,0:nq) )
		q2_ups_in = 0.;  h2_ups_in=0.
		allocate(fname_2(jsin2))
	end if

	j_in = 0
	jsin1 = 0
	jsin2 = 0

	do i = 1, j_inn
		if(indices(i,1,1)==1 .and. indices(i,1,j_inlen(i))==1) then
			j_in=j_in+1
		
			fname_in(j_in)=flowname(i)
			write(*,*) 'inflow(i=1)', trim(fname_in(j_in))
			
			do k=1,jsizemax
				if(f_param(i,k).ne.f_param(1,k)) then
					write(*,*) 'error---Time of Inflow'
				end if

				t_xx=f_param(i,k)
				q_ups_in(j_in,k-1)=f_value(i,k)

				if(i_sec_hour==1) then
					t_hyd(k-1)=t_xx
				else
					t_hyd(k-1)=t_xx*3600.
				end if
			end do

		else if(indices(i,2,1)==1) then
			if(j_side_j1 == 2) then
				jsin1=jsin1+1

				fname_1(jsin1)=flowname(i)
				write(*,*) 'inflow(j=1)',  trim(fname_1(jsin1))

				do k=1,jsizemax
					if(f_param(i,k).ne.f_param(1,k)) then
						write(*,*) 'error---Time of Inflow'
					end if

					q1_ups_in(jsin1,k-1)=f_value(i,k)
				end do
			end if
		else if(indices(i,2,1)==nj4) then
			if(j_side_j2 == 2) then
				jsin2=jsin2+1
				
				fname_2(jsin2)=flowname(i)
				write(*,*) 'inflow(j=ny)', trim(fname_2(jsin2))

				do k=1,jsizemax
					if(f_param(i,k).ne.f_param(1,k)) then
						write(*,*) 'error---Time of Inflow'
					end if

					q2_ups_in(jsin2,k-1)=f_value(i,k)
				end do
			end if
		end if
	end do
	
	thstart = t_hyd(0)

	do n=0,nq
		t_hyd(n) = t_hyd(n)-thstart
	end do
	
	if( n_q_cell>0 ) then
		do n=0,nqtcell
			t_qc_dis(n) = t_qc_dis(n)-thstart
		end do
	end if

	qp = q_ups_in(1,0)
!
	if( j_wl==2 ) then

		CALL CG_IRIC_READ_FUNCTIONALSIZE_F('stage_at_downstream',tmpint,ier)

		allocate( xtmp(tmpint),ytmp(tmpint) )

		CALL CG_IRIC_READ_FUNCTIONAL_REALSINGLE_F('stage_at_downstream',xtmp,ytmp,ier)

		DO I= 0,nq
			h_dse(i) = ytmp(i+1)
		ENDDO

		DEALLOCATE(xtmp, STAT = ier)
		DEALLOCATE(ytmp, STAT = ier)
	end if

	if( j_rain==2 ) then

		CALL CG_IRIC_READ_FUNCTIONALSIZE_F('rain_time',tmpint,ier)

		allocate(xtmp(tmpint),ytmp(tmpint))
		CALL CG_IRIC_READ_FUNCTIONAL_REALSINGLE_F('rain_time',xtmp,ytmp,ier)

		DO I= 0,nq
			rain(i) = ytmp(i+1)
		ENDDO

		DEALLOCATE(xtmp, STAT = ier)
		DEALLOCATE(ytmp, STAT = ier)
	end if

!h160104 Xrain対応
	if( j_rain==3 ) then
		CALL CG_IRIC_READ_GRID_FUNCTIONALTIMESIZE_F('Xrain',tmpint,ier)
		if(ier .ne. 0) then
			write(*,*) 'Xrain Mapping Error!'
			stop
		end if
		allocate(xtmp_d(tmpint),t_rain(tmpint))
		CALL CG_IRIC_READ_GRID_FUNCTIONALTIME_F('Xrain',xtmp_d,ier)
		
		nr=tmpint
		do i=1,nr
			write(*,*) 'Xrain time',xtmp_d(i)
			t_rain(i)=xtmp_d(i)-xtmp_d(1)
		end do
		
		DEALLOCATE(xtmp_d, STAT = ier)
	
	end if
!
	etime = t_hyd(nq)

	if( ster<0. ) ster = etime
	
	total_bomb = 0
	
	do j=1,ny
		do i=1,nx
			total_bomb = total_bomb+ob(i,j)
		end do
	end do
	
	if( total_bomb==0 ) bomber_time = etime
	
	i_bomb = 0

! ---- 河床粗度に関する設定 -----

	do j=1,ny
		do i=1,nx
			snmm(i,j) = chcond4(i,j)

			if(snmm(i,j)<=0.d0) then
				write(*,'(a30,i5,a1,i5,a4,f10.3)') &
					'Manning roughness coefficient(',i,',',j,') is',snmm(i,j)
				write(*,*) 'This coefficient should be larger than 0'
				stop
			end if
		end do
	end do
	
	do j=1,ny
		do i=1,nx-1
			sn_up(i,j) = (snmm(i,j)+snmm(i+1,j))*.5d0
		end do
	end do
	
	do j=1,ny
		sn_up(0,j) = snmm(1,j)
		sn_up(nx,j) = snmm(nx,j)
	end do

	do j=1,ny-1
		do i=1,nx
			sn_vp(i,j) = (snmm(i,j)+snmm(i,j+1))*.5d0
		end do
	end do
	
	do i=1,nx
		sn_vp(i,0) = snmm(i,1)
		sn_vp(i,ny) = snmm(i,ny)
	end do

!inoue ----	建物阻害率 hamaki ver ----
	do j=1,ny
		do i=1,nx
			share(i,j) = sh4(i,j)
			gam_v(i,j) = 1.0d0 - share(i,j)
			if(j_cip == 2 .and. gam_v(i,j) /= 1.0d0) then	!h150113 建物阻害率ありでCIPならエラーメッセージ
				write(*,*) 'When considering buildings occupy, please select the upwind scheme.'
				stop
			end if
			if(gam_v(i,j)<=0.05) then
				write(*,'(a30,i5,a1,i5,a4,f10.3)') &
					'Share of Buildings(',i,',',j,') is',share(i,j)
				write(*,*) 'This coefficient should be smaller than 0.95'
!inoue				if(j_cdata == 0)	then
!					write(*,*) 'Please reset Share of Buildings'
!					stop
!				else
					write(*,*) 'Correct Share of Buildings to 0.95'
					gam_v(i,j)=0.05
!				endif
			end if
			gam_e(i,j) = 1.0d0 - sqrt(1.0d0 - gam_v(i,j))
			if(j_gam == 1)	gam_e(i,j) = gam_v(i,j)
		end do
	end do
	
	do j=1,ny
		do i=1,nx-1
			gam_e_up(i,j) = max(gam_e(i,j),gam_e(i+1,j))	!gam_eは最大値
			gam_v_up(i,j) = (gam_v(i,j)+gam_v(i+1,j))*.5d0	!gam_vは平均値
		end do
	end do
	
	do j=1,ny
		gam_e_up(0,j) = gam_e(1,j)
		gam_e_up(nx,j) = gam_e(nx,j)
		gam_v_up(0,j) = gam_v(1,j)
		gam_v_up(nx,j) = gam_v(nx,j)
	end do

	do j=1,ny-1
		do i=1,nx
			gam_e_vp(i,j) = max(gam_e(i,j),gam_e(i,j+1))	!gam_eは最大値
			gam_v_vp(i,j) = (gam_v(i,j)+gam_v(i,j+1))*.5d0	!gam_vは平均値
		end do
	end do
	
	do i=1,nx
		gam_e_vp(i,0) = gam_e(i,1)
		gam_e_vp(i,ny) = gam_e(i,ny)
		gam_v_vp(i,0) = gam_v(i,1)
		gam_v_vp(i,ny) = gam_v(i,ny)
	end do
! ------------------------------------------

	cw = 0.0d0
	sigma_c = 1.d0

	kend = int(etime/dt+.5)
	kmod = int(tuk/dt+.5)

! -------------------------------------------

	call avgeo
!
	call cell2grid(eta,z0)

	if(j_wl==0.and.h_down<-10.) then
		write(*,*) 'Downstrem Water Suarfece Value Wrong!!!'
		stop
	end if

	if(j_wl.ne.0) h_down=-999.
!
	call gcoefs(0)
	
	if( n_q_cell>0 ) then
		do k=1,n_q_cell
			sum_area = 0
			do j=1,ny
				do i=1,nx
					if( q_cell(i,j)==k ) then
						sum_area = sum_area+area(i,j)
					end if
				end do
			end do
			
			ni_dis_cell(k) = sum_area
		end do
	
		ni_dis_cell(0) = 1.d0
	end if

!  ------------------------------------------

	errmax=hmin*0.01
	hmin2=hmin*2.

	snu_0 = snu00

!
! ------ Initial Condition -----
!
	if( j_wl==0 ) then
		hnx  =h_down
	else if( j_wl==2 ) then
		hnx = h_dse(0)
	else if( j_wl==3 ) then
		hnx = emin(nx)
	end if

	call initl(qp,hnx,snu_0,i_flow,h_slope,x_bk,h_slope_1,h_slope_2)
	call bound_u(yu,ijobst_u)
	call bound_v(yv,ijobst_v)
	call bound_u(yun,ijobst_u)
	call bound_v(yvn,ijobst_v)

!   等流計算による境界水位時系列の計算

	call hqtcal_inn(nq,h_down)

	if( j_side_j1==2 ) call hqtcal_j1(nq)
	if( j_side_j2==2 ) call hqtcal_j2(nq)

!  ------------------------------------------

	time = 0.
	icount = 0
	iofrg = 0
	ndeposit = 0
	!
   ! ----- Read tempfile for hot start ----
   !
   if(i_re_flag_i == 1) then
      !
        call read_hotstart_conditions(icount, tmp_file_i)
        
        !online
        time = 0.
		call gcoefs(1)
		!
		is = -1
		do ii=0,n_rest-1
		if(opt_tmp(ii) < time) is=ii
		end do
		!
		i_tmp_count = is+1
		!
		
		!Pythonからの連続計算の場合
		n_rest = 1
		opt_tmp(0) = time + opt_tmp(0) 
		
		
   end if
	
! ------------- Time loop -------------

!$	call omp_set_num_threads(n_parallel)

!$omp parallel

!	2000 continue
	
	do			!  Time loop

!$omp single

!	icount = icount+1

	if( time>bomber_time.and.i_bomb==0 ) then

		do j=1,ny
			do i=1,nx
				if( ob(i,j)==1 ) then
					ijo_in(i,j) = 0
					ijobst(i,j) = 0
					ijobst(i-1,j) = 0
					ijobst(i,j-1) = 0
					ijobst(i-1,j-1) = 0
				end if
			end do
		end do
			
		do j=1,ny
			do i=1,nx
				if( ijo_in(i,j)==1 ) then
					lcb(i,j) = 0
					ijobst(i,j) = 1
					ijobst(i-1,j) = 1
					ijobst(i,j-1) = 1
					ijobst(i-1,j-1) = 1
				end if
			end do
		end do

		do j=1,ny
			do i=0,nx
				if( ijobst(i,j)==1 .and. ijobst(i,j-1)==1 ) then
					ijobst_u(i,j) = 0
				else
					ijobst_u(i,j) = 1
				end if
			end do
		end do

		do j=0,ny
			do i=1,nx
				if( ijobst(i,j)==1 .and. ijobst(i-1,j)==1 ) then
					ijobst_v(i,j) = 0
				else
					ijobst_v(i,j) = 1
				end if
			end do
		end do
			
		i_bomb = 1

	end if

	do i=1,j_in
		if(time<=0.) then
			h_input_in(i)=h_ups_in(i,0)
			q_input_in(i)=q_ups_in(i,0)
		else if(time>t_hyd(nq)) then
			h_input_in(i)=h_ups_in(i,nq)
			q_input_in(i)=q_ups_in(i,nq)
		else
			do n=1,nq
				if(time>=t_hyd(n-1).and.time<=t_hyd(n)) then
					sst=(time-t_hyd(n-1))/(t_hyd(n)-t_hyd(n-1))
					h_input_in(i)=h_ups_in(i,n-1)+(h_ups_in(i,n)-h_ups_in(i,n-1))*sst
					q_input_in(i)=q_ups_in(i,n-1)+(q_ups_in(i,n)-q_ups_in(i,n-1))*sst
				end if
			end do
		end if
	end do

	qp=q_input_in(1)

	if(j_side_j1 == 2) then
		do i=1,jsin1
			if(time<=0.) then
				h_in_j1(i)=h1_ups_in(i,0)
				q_in_j1(i)=q1_ups_in(i,0)
			else if(time>t_hyd(nq)) then
				h_in_j1(i)=h1_ups_in(i,nq)
				q_in_j1(i)=q1_ups_in(i,nq)
			else
				do n=1,nq
					if(time>=t_hyd(n-1).and.time<=t_hyd(n)) then
						sst=(time-t_hyd(n-1))/(t_hyd(n)-t_hyd(n-1))
						h_in_j1(i)=h1_ups_in(i,n-1)+(h1_ups_in(i,n)-h1_ups_in(i,n-1))*sst
						q_in_j1(i)=q1_ups_in(i,n-1)+(q1_ups_in(i,n)-q1_ups_in(i,n-1))*sst
					end if
				end do
			end if
		end do
	end if

	if(j_side_j2 == 2) then
		do i=1,jsin2
			if(time<=0.) then
				h_in_j2(i)=h2_ups_in(i,0)
				q_in_j2(i)=q2_ups_in(i,0)
			else if(time>t_hyd(nq)) then
				h_in_j2(i)=h2_ups_in(i,nq)
				q_in_j2(i)=q2_ups_in(i,nq)
			else
				do n=1,nq
					if(time>=t_hyd(n-1).and.time<=t_hyd(n)) then
						sst=(time-t_hyd(n-1))/(t_hyd(n)-t_hyd(n-1))
						h_in_j2(i)=h2_ups_in(i,n-1)+(h2_ups_in(i,n)-h2_ups_in(i,n-1))*sst
						q_in_j2(i)=q2_ups_in(i,n-1)+(q2_ups_in(i,n)-q2_ups_in(i,n-1))*sst
					end if
				end do
			end if
		end do
	end if

! ------------------------------------------

	if(j_rain==2) then
		if(time<=0.) then
			rain_t=rain(0)
		else if(time>t_hyd(nq)) then
			rain_t=rain(nq)
		else
			do n=1,nq
				if(time>=t_hyd(n-1).and.time<=t_hyd(n)) then
					sst=(time-t_hyd(n-1))/(t_hyd(n)-t_hyd(n-1))
					rain_t=rain(n-1)+(rain(n)-rain(n-1))*sst
				end if
			end do
		end if

		rain_t=rain_t*0.001/3600.
		rain_t2=rain_t		!h160105 全領域一定の雨の場合

	else if(j_rain==3) then	!h160105 Xrain対応
		if(time<=0.) then
			tmpint=1
		else if(time>t_rain(nr)) then
			tmpint=nr
		else if(nr>2) then
			do n=2,nr
				if( time>=t_rain(n-1) .and. time<t_rain(n) ) then
					tmpint=n-1
				end if
			end do
		end if
		
		CALL CG_IRIC_READ_GRID_FUNCTIONAL_REAL_CELL_F('Xrain',tmpint,rain_t2,ier)
		rain_t2=rain_t2*0.001/3600.
		rain_t2=rain_t2*0.1		!hスケールファクター

	end if

! -------------------------------------------

	if( n_q_cell>0 ) then
		do k=1,n_q_cell
			if(time<=0.) then
				q_tmp_dis(k) = q_cell_dis(0,k)
			else if(time>t_qc_dis(nqtcell)) then
				q_tmp_dis(k) = q_cell_dis(nqtcell,k)
			else
				do n=1,nqtcell
					if( time>=t_qc_dis(n-1).and.time<=t_qc_dis(n) ) then
						sst=(time-t_qc_dis(n-1))/(t_qc_dis(n)-t_qc_dis(n-1))
						q_tmp_dis(k)=q_cell_dis(n-1,k)+(q_cell_dis(n,k)-q_cell_dis(n-1,k))*sst
					end if
				end do
			end if
		end do
		
		do j=1,ny
			do i=1,nx
				cell_discharge(i,j) = q_tmp_dis(q_cell(i,j))/ni_dis_cell(q_cell(i,j))
			end do
		end do
		
	end if

! -------------------------------------------

	if(h_down<-100.) then
		if(time<=0.) then
			hnx=h_dse(0)
		else if(time>t_hyd(nq)) then
			hnx=h_dse(nq)
		else
			do n=1,nq
				if(time>=t_hyd(n-1).and.time<=t_hyd(n)) then
					sst=(time-t_hyd(n-1))/(t_hyd(n)-t_hyd(n-1))
					hnx=h_dse(n-1)+(h_dse(n)-h_dse(n-1))*sst
				end if
			end do
		end if
	else
		hnx=h_down
	end if

! ----- 水位・流量の境界条件 -----

	call upstream_inn()

	if( j_side_j1==2 ) call upstream_j1()
	if( j_side_j2==2 ) call upstream_j2()

	if( j_wl/=3 ) then
		call downstream(hnx)
	end if

		! ユーザがGUI上で "STOP" ボタンを押して実行をキャンセルしたか確認
		CALL IRIC_CHECK_CANCEL_F(istatus)
		if(istatus == 1) then
			write(*,*) "Solver is stopped because the STOP button was clicked."
			
			call system_clock(cal_t2, t_rate, t_max)	!h160105 計算終了時時刻
			if ( cal_t2 < cal_t1 ) then
				t_diff = (t_max - cal_t1) + cal_t2 + 1
			else
				t_diff = cal_t2 - cal_t1
			endif
			write(*,*) "Calcuration time",real(t_diff/t_rate),"sec."
			write(*,*) "Calcuration time",real(t_diff/t_rate)/60.,"min."
			write(*,*) "Calcuration time",real(t_diff/t_rate)/3600.,"hour."
			stop
		end if

!$omp end single


! ------------ 計算結果のアウトプット -------------

!	if ( icount == 1.or.mod(icount-1,kmod) == 0 ) then	!h time=0も出力
	if( icount==0 .or. mod(icount,kmod)==0 ) then
!
		if(iofrg==0) then
			call hsxxcal(eta0,z,hs,hsxx,rain_t2,rain_xx, ibc, ibcg)
			iofrg = 1
		else
			call hsxxcal(eta,z,hs,hsxx,rain_t2,rain_xx, ibc, ibcg)
		end if

		call uxxyycal(yu,yv,uxx,uyy)

!$omp do private(i,j,v_temp)
		do j=0,ny
			do i=0,nx
				v_temp = sqrt(uxx(i,j)**2+uyy(i,j)**2)
				h_max(i,j) = max(h_max(i,j),hsxx(i,j))
				v_max(i,j) = max(v_max(i,j),v_temp)
			end do
		end do

!$omp single

			if( time>=t_out_start ) then
				qptemp = qp

				! guiがcgnsファイルを読込中か否かを判定
				do
					call iric_check_lock_f(condFile, istatus)
					if(istatus == 1) then
						call sleep(1)
					elseif(istatus == 0)then  !読込中でなければdoループを抜ける
						call iric_write_sol_start_f(condFile, ier)   
						exit
					end if
				end do

				CALL Write_CGNS		&
					(condFile,Index,time,q_input_in,q_in_j1,q_in_j2		&
							 ,fname_in,fname_1,fname_2,j_in,jsin1,jsin2      &
							 ,j_side_j1,j_side_j2,nx,ny,x,y,uxx,uyy,hsxx,z,z0	&
					,iofrg,h_max,v_max,j_rain,rain_xx, ibcg)

				call cg_iric_flush_f(condFile, fid, ier)
				call iric_write_sol_end_f(condFile, ier)

				Index = Index+1
			end if

			ioutflag=ioutflag+1

			if(ioutflag==1) then
				write(*,'(a10,a10,a10)') 'time','q_input','h_down'
			end if

			q_inptotal=0.

			do i=1,j_in
				q_inptotal=q_inptotal+q_input_in(i)
			end do

			if(j_side_j1 == 2) then
				do i=1,jsin1
					q_inptotal=q_inptotal+q_in_j1(i)
				end do
			endif

			if(j_side_j2 == 2) then
				do i=1,jsin2
					q_inptotal=q_inptotal+q_in_j2(i)
				end do
			endif

			if(time<t_out_start) then
				write(*,'(f10.3,3f10.4)') time,q_inptotal,hnx,rain_t*1000.*3600.
			else
				write(*,'(f10.3,3f10.4,a4)') time,q_inptotal,hnx,rain_t*1000.*3600.,'out'
			end if
		
!$omp end single
		
	end if


	
!$omp single
	!Output for hotstart
	if(i_re_flag_o == 1 .and. time > opt_tmp(i_tmp_count)) then
	
		tmp_file_o(i_tmp_count)=trim(tmp_pass)//"/"//tmp_file_o(i_tmp_count)
		
		call output_hotstart(icount, tmp_file_o(i_tmp_count))
		
		i_tmp_count = i_tmp_count + 1
      if(i_tmp_count >= n_rest) i_re_flag_o = 0
		
	end if
!$omp end single

	
	
! ---------  flow calculation  ---------

	call hcal(errmax,err,lcount,alh)

!docon add start#####################################################################################
!$omp single
	q_swap(:,:)=0.0d0
	bc_qswap(:)=0.0d0 ; gt_qswap(:)=0.0d0 ; p_qswap(:)=0.0d0    !140602
	call culvert
	call gate
	call pump
!$omp end single

!$omp do
	do j=1,ny
		do i=1,nx-1
			if( ijo_in(i,j)==1 ) goto 4000
			if(q_swap(i,j)/=0.0d0) then
				hn(i,j) = hn(i,j)+q_swap(i,j)*dt/area(i,j)
			end if
			hs(i,j) = hn(i,j)-eta(i,j)
				
			if( hs(i,j)<=hmin ) then
				lcb(i,j) = 0
				hs(i,j) = hmin
				hn(i,j) = hs(i,j)+eta(i,j)
				if( yun(i-1,  j)<0. ) yun(i-1,  j) = 0.
				if( yun(i  ,  j)>0. ) yun(i  ,  j) = 0.
				if( yvn(i  ,j-1)<0. ) yvn(i  ,j-1) = 0.
				if( yvn(i  ,  j)>0. ) yvn(i  ,  j) = 0.
			else
				lcb(i,j) = 1
			end if
			4000 continue
		end do
	end do
!docon add end#######################################################################################

	call bound_h(hn,hs,eta)
	call bound_u(yun,ijobst_u)
	call bound_v(yvn,ijobst_v)

	call diffusion(cw)
	call bound_u(yun,ijobst_u)
	call newgrd_u
	call bound_v(yvn,ijobst_v)
	call newgrd_v

	if( j_cip==1 ) then
!ino		call upwind2d_u(yun,gux,guy)
	else
!ino		call dcip2d_u(yun,gux,guy)
	end if

	call dryck_u(yun,hs,gux,guy)
	call bound_u(yun,ijobst_u)

	if( j_cip==1 ) then
!ino		call upwind2d_v(yvn,gvx,gvy)
	else
!ino		call dcip2d_v(yvn,gvx,gvy)
	end if

	call dryck_v(yvn,hs,gvx,gvy)
	call bound_v(yvn,ijobst_v)

!$omp do
	do j=-1,ny+1
		do i=-1,nx+1
			yu(i,j) = yun(i,j)
			yv(i,j) = yvn(i,j)
		end do
	end do

	call uvpcal(yun,yvn,up,vp,hs)
	call uxuycal(up,vp,ux,uy)

	call taustacal(snu00)

	call snucal(snu00,a_snu,b_snu)

	call hshift(hn,h)

	call bound_h(h,hs,eta)

!$omp do
	do j=1,ny
		do i=1,nx
			if(isnan(h(i,j))) then
				write(*,*) 'Calculation is failure!'
				stop
			end if
		end do
	end do

!	if(icount-1>kend) goto 3000
	if( icount>kend ) exit

!$omp barrier

!$omp single
	icount = icount+1
	time = dble(icount)*dt
!$omp end single

!	goto 2000

!	3000 continue

	end do

! ------------ end time loop ------------

!$omp end parallel

	call cg_close_f(fid,ier)

	write(*,*) "Finish", ier

	call system_clock(cal_t2, t_rate, t_max)	!h160105 計算終了時時刻
	if ( cal_t2 < cal_t1 ) then
		t_diff = (t_max - cal_t1) + cal_t2 + 1
	else
		t_diff = cal_t2 - cal_t1
	endif
	write(*,*) "Calcuration time",real(t_diff/t_rate),"sec."
	write(*,*) "Calcuration time",real(t_diff/t_rate)/60.,"min."
	write(*,*) "Calcuration time",real(t_diff/t_rate)/3600.,"hour."

end program nays2d_flood_parallel

!docon add start#####################################################################################
subroutine culvert
	use variables
	use common_hh
	use box_gate_pump
	implicit none
	integer :: i,j,muki
	double precision :: h1,h2,qlim,q_bc,bc_wlin,bc_wlot

	if(bc_num>0) then
L1:	do i=1,bc_num
			q_bc=0.0d0  !140602念のためゼロクリア
			if(bc_inout(i)==1) cycle L1 !出口セルの場合
L2:		do j=1,bc_num
				if(i==j) cycle L2
				if(bc_inout(j)==0) cycle L2 !入口セルの場合
				if(bc_couple_num(i)/=bc_couple_num(j)) cycle L2 !BC番号の一致の確認

				!140602bc_wlin=h(bc_indx(i,1,1),bc_indx(i,2,1))        !内水位(入口側水位)
				!140602bc_wlot=h(bc_indx(j,1,1),bc_indx(j,2,1))        !外水位(出口側水位)
				bc_wlin=hn(bc_indx(i,1,1),bc_indx(i,2,1))        !内水位(入口側水位)
				bc_wlot=hn(bc_indx(j,1,1),bc_indx(j,2,1))        !外水位(出口側水位)

				if(bc_wlin>bc_wlot) then
					!ボックスカルバート入口地点の水位が敷高以下の場合何もしない
					if(bc_wlin<=bc_base(i)) cycle L1
					!ボックスカルバート入口地点の水深がhmin以下(演算誤差考慮)の場合何もしない
					if(bc_wlin<=eta(bc_indx(i,1,1),bc_indx(i,2,1))+hmin+1.0e-6) cycle L1
					h1=bc_wlin-bc_base(i)
					h2=bc_wlot-bc_base(i)
					!カルバート敷高を超えるボリューム分移動可能
					qlim=(dmin1(bc_wlin-eta(bc_indx(i,1,1),bc_indx(i,2,1))-hmin,bc_wlin-bc_base(i))) &
						*area(bc_indx(i,1,1),bc_indx(i,2,1))/dt
					muki=1
				else
					!ボックスカルバート出口地点の水位が敷高以下の場合何もしない
					if(bc_wlot<=bc_base(i)) cycle L1
					!ボックスカルバート出口地点の水深がhmin以下(演算誤差考慮)の場合何もしない
					if(bc_wlot<=eta(bc_indx(j,1,1),bc_indx(j,2,1))+hmin+1.0e-6) cycle L1
					h1=bc_wlot-bc_base(i)
					h2=bc_wlin-bc_base(i)
					!カルバート敷高を超えるボリューム分移動可能
					qlim=(dmin1(bc_wlot-eta(bc_indx(j,1,1),bc_indx(j,2,1))-hmin,bc_wlot-bc_base(i))) &
						*area(bc_indx(j,1,1),bc_indx(j,2,1))/dt
					muki=-1
				end if
				
				if(h1<0.0d0) then
					h1=0.0d0 ; cycle L1
				end if
				if(h2<0.0d0) h2=0.0d0

				!カルバート公式
				if(h2>=bc_height(i)) then
					q_bc=bc_c1(i)*bc_width(i)*bc_height(i)*sqrt(g*2.0d0*(h1-h2))*bc_ren(i)
				elseif(bc_c2_yesno(i)==1.and.h1>=1.5d0*bc_height(i)) then
					q_bc=bc_c2(i)*bc_width(i)*bc_height(i)*sqrt(g*2.0d0*h1)*bc_ren(i)
				else
					if(h1>=1.5d0*h2) h2=2.0d0*h1/3.0d0
					q_bc=bc_c3(i)*bc_width(i)*h2*sqrt(g*2.0d0*(h1-h2))*bc_ren(i)
				end if

				!移動可能量で制御
				if(qlim<q_bc) q_bc=qlim
							
				if(muki==1) then
					!入口セル(流出)
					q_swap(bc_indx(i,1,1),bc_indx(i,2,1))=q_swap(bc_indx(i,1,1),bc_indx(i,2,1))-q_bc
					!出口セル(流入)
					q_swap(bc_indx(j,1,1),bc_indx(j,2,1))=q_swap(bc_indx(j,1,1),bc_indx(j,2,1))+q_bc
					!出力用140602
					bc_qswap(i)=-q_bc
					bc_qswap(j)=q_bc
				else
					!出口セル(流出)
					q_swap(bc_indx(j,1,1),bc_indx(j,2,1))=q_swap(bc_indx(j,1,1),bc_indx(j,2,1))-q_bc
					!入口セル(流入)
					q_swap(bc_indx(i,1,1),bc_indx(i,2,1))=q_swap(bc_indx(i,1,1),bc_indx(i,2,1))+q_bc
					!出力用140602
					bc_qswap(i)=q_bc
					bc_qswap(j)=-q_bc
				end if
			end do L2
		end do L1
	end if
end subroutine

subroutine gate
	use variables
	use common_hh
	use box_gate_pump
	implicit none
	integer :: i,j,k,muki,gt_flg
	double precision :: h1,h2,sst,qlim,q_gt,gt_wlin,gt_wlot

	if(gt_num>0) then
L1:	do i=1,gt_num
			q_gt=0.0d0  !140602念のためゼロクリア
			if(gt_inout(i)==1) cycle L1 !出口セルの場合

			gt_flg=0    !樋門が開の場合
			if(gt_method(i)==1) then    !樋門開閉状態をデータで設定する場合
				!gt_stateは開=0、閉=1
				if(time<=0.0d0) then
					gt_flg=gt_state(i,1)
				else if(time>=gt_time(i,gt_size(i))) then
					gt_flg=gt_state(i,gt_size(i))
				else
					do k=1,gt_size(i)-1
						if(time>=gt_time(i,k).and.time<gt_time(i,k+1)) gt_flg=gt_state(i,k)
					end do
				end if
			end if

			if(gt_flg==0) then  !樋門が開の場合
L2:			do j=1,gt_num
					if(i==j) cycle L2
					if(gt_inout(j)==0) cycle L2 !入口セルの場合
					if(gt_couple_num(i)/=gt_couple_num(j)) cycle L2 !Gate番号の一致の確認

					!140602gt_wlin=h(gt_indx(i,1,1),gt_indx(i,2,1))        !内水位(入口側水位)
					gt_wlin=hn(gt_indx(i,1,1),gt_indx(i,2,1))        !内水位(入口側水位)
					if(gt_outwl(j)==0) then                         !外水位の指定がない場合
						!140602gt_wlot=h(gt_indx(j,1,1),gt_indx(j,2,1))    !外水位
						gt_wlot=hn(gt_indx(j,1,1),gt_indx(j,2,1))    !外水位
					else                                            !外水位の指定がある場合
						if(time<=0.0d0) then
							gt_wlot=gt_wl(j,1)
						else if(time>=gt_time2(j,gt_size2(j))) then
							gt_wlot=gt_wl(j,gt_size2(j))
						else
							do k=1,gt_size2(j)-1
								if(time>=gt_time2(j,k).and.time<gt_time2(j,k+1)) then
									sst=(time-gt_time2(j,k))/(gt_time2(j,k+1)-gt_time2(j,k))
									gt_wlot=gt_wl(j,k)+(gt_wl(j,k+1)-gt_wl(j,k))*sst
								end if
							end do
						end if
					end if

					!開閉自動で内水位≦外水位の場合何もしない
					if(gt_method(i)==0.and.gt_wlin<=gt_wlot) cycle L1
					!開閉自動で水深がhmin以下(演算誤差考慮)の場合何もしない
					if(gt_method(i)==0.and.gt_wlin<=eta(gt_indx(i,1,1),gt_indx(i,2,1))+hmin+1.0e-6) cycle L1
					!開閉手動で水深がhmin以下(演算誤差考慮)の場合何もしない
					if(gt_method(i)==1) then
						if(gt_wlin>gt_wlot) then
							if(gt_wlin-eta(gt_indx(i,1,1),gt_indx(i,2,1))<=hmin+1.0e-6) cycle L1
						else
							if(gt_wlot-eta(gt_indx(j,1,1),gt_indx(j,2,1))<=hmin+1.0e-6) cycle L1
						end if
					end if

					if(gt_wlin>gt_wlot) then
						!樋門入口/出口地点の水位が敷高以下の場合何もしない
						if(gt_wlin<=gt_base(i)) cycle L1
						h1=gt_wlin-gt_base(i)
						h2=gt_wlot-gt_base(i)
						!樋門敷高を超えるボリューム分移動可能
						qlim=(dmin1(gt_wlin-eta(gt_indx(i,1,1),gt_indx(i,2,1))-hmin,gt_wlin-gt_base(i))) &
							*area(gt_indx(i,1,1),gt_indx(i,2,1))/dt
						muki=-1
					else
						!樋門入口/出口地点の水位が敷高以下の場合何もしない
						if(gt_wlot<=gt_base(i)) cycle L1
						h1=gt_wlot-gt_base(i)
						h2=gt_wlin-gt_base(i)
						!樋門敷高を超えるボリューム分移動可能
						qlim=(dmin1(gt_wlot-eta(gt_indx(j,1,1),gt_indx(j,2,1))-hmin,gt_wlot-gt_base(i))) &
							*area(gt_indx(j,1,1),gt_indx(j,2,1))/dt
						muki=1
					end if

					if(h1<0.0d0) then
						h1=0.0d0 ; cycle L1
					end if
					if(h2<0.0d0) h2=0.0d0

					!樋門公式
					if(h2>=gt_height(i)) then
						q_gt=gt_c1(i)*gt_width(i)*gt_height(i)*sqrt(g*2.0d0*(h1-h2))*gt_ren(i)
					elseif(gt_c2_yesno(i)==1.and.h1>=1.5d0*gt_height(i)) then
						q_gt=gt_c2(i)*gt_width(i)*gt_height(i)*sqrt(g*2.0d0*h1)*gt_ren(i)
					else
						if(h1>=1.5d0*h2) h2=2.0d0*h1/3.0d0
						q_gt=gt_c3(i)*gt_width(i)*h2*sqrt(g*2.0d0*(h1-h2))*gt_ren(i)
					end if

					!移動可能量で制御
					if(qlim<q_gt) q_gt=qlim
					if(muki==1) then
						q_swap(gt_indx(i,1,1),gt_indx(i,2,1))=q_swap(gt_indx(i,1,1),gt_indx(i,2,1))+q_gt	!入口セル(流入)
						if(gt_outwl(j)==0) q_swap(gt_indx(j,1,1),gt_indx(j,2,1))= &
							q_swap(gt_indx(j,1,1),gt_indx(j,2,1))-q_gt	!出口セル(流出)
						!出力用140602
						gt_qswap(i)=q_gt
						if(gt_outwl(j)==0) gt_qswap(j)=-q_gt
					else
						q_swap(gt_indx(i,1,1),gt_indx(i,2,1))=q_swap(gt_indx(i,1,1),gt_indx(i,2,1))-q_gt	!入口セル(流出)
						if(gt_outwl(j)==0) q_swap(gt_indx(j,1,1),gt_indx(j,2,1))= &
							q_swap(gt_indx(j,1,1),gt_indx(j,2,1))+q_gt	!出口セル(流入)
						!出力用140602
						gt_qswap(i)=-q_gt
						if(gt_outwl(j)==0) gt_qswap(j)=q_gt
					end if
				end do L2
			end if
		end do L1
	end if
end subroutine

subroutine pump
	use variables
	use common_hh
	use box_gate_pump
	implicit none
    integer :: i,j,k,l
	double precision :: h1,h2,sst,qlim,q_p,p_wlin
    double precision,allocatable,dimension(:) :: q_pin	!take160401

	if(p_num>0) then
        allocate(q_pin(p_num)) !take160401
L1:	do i=1,p_num
            q_pin(i)=0.0d0 !take160401

!take160401            q_p=0.0d0  !140602念のためゼロクリア
            if(p_inout(i)==1) then
              if(p_mxindx(i)>1) then
                write(*,*) 'Pump-cell is not allowed multiple selection with pump outlet'
                stop
              end if
            end if

			if(p_inout(i)==1) cycle L1  !出口セルの場合

			!140602p_wlin=h(p_indx(i,1,1),p_indx(i,2,1))   !内水位(入口側水位)
            !p_wlin=hn(p_indx(i,1,1),p_indx(i,2,1))   !内水位(入口側水位)

            
L5:         do l=1,p_mxindx(i)
            q_p=0.0d0  !140602念のためゼロクリア	!take160401

			!ポンプ入口地点の水深がhmin以下(演算誤差考慮)の場合何もしない
            p_wlin=hn(p_indx(i,1,l),p_indx(i,2,l))   !内水位(入口側水位)

            !if(p_wlin<=eta(p_indx(i,1,1),p_indx(i,2,1))+hmin+1.0e-6) cycle L1
            if(p_wlin>=eta(p_indx(i,1,l),p_indx(i,2,l))+hmin+1.0e-6) then

			if(p_method(i)==0) then     !ポンプ排水量が自動の場合
				if(p_wlin>=p_startope(i).and.p_wlin<p_stopope(i)) then
					q_p=p_qmax(i)
                    qlim=(p_wlin-p_startope(i))*area(p_indx(i,1,l),p_indx(i,2,l))/dt
					!140602qlim=dmin1(qlim,(h(p_indx(i,1,1),p_indx(i,2,1))-eta(p_indx(i,1,1),p_indx(i,2,1))-hmin)* &
					!140602    area(p_indx(i,1,1),p_indx(i,2,1))/dt)
                    qlim=dmin1(qlim,(hn(p_indx(i,1,l),p_indx(i,2,l))-eta(p_indx(i,1,l),p_indx(i,2,l))-hmin)* &
                        area(p_indx(i,1,l),p_indx(i,2,l))/dt)
				else
					q_p=0.0d0
				end if
			elseif(p_method(i)==1) then !ポンプ排水量が手動の場合
				if(time<=0.0d0) then
					q_p=p_qout(i,1)
				else if(time>=p_time(i,p_size(i))) then
					q_p=p_qout(i,p_size(i))
				else
					do k=1,p_size(i)-1
						if(time>=p_time(i,k).and.time<p_time(i,k+1)) then
							sst=(time-p_time(i,k))/(p_time(i,k+1)-p_time(i,k))
							q_p=p_qout(i,k)+(p_qout(i,k+1)-p_qout(i,k))*sst
						end if
					end do
				end if
			    qlim=(p_wlin-eta(p_indx(i,1,l),p_indx(i,2,l))-hmin)*area(p_indx(i,1,l),p_indx(i,2,l))/dt
                            !if (l.eq.2) 
                            !write(*,*) p_indx(i,1,l), p_indx(i,1,l), qlim,q_p
			end if

			!移動可能量で制御
			if(qlim<q_p) q_p=qlim
            q_swap(p_indx(i,1,l),p_indx(i,2,l))=q_swap(p_indx(i,1,l),p_indx(i,2,l))-q_p !入口セル(流出)
            q_pin(i)=q_pin(i)+q_p !take160401
			!出力用140602
			p_qswap(i)=-q_p

            end if

            end do L5  

			!出口セルが存在する場合
L2:		do j=1,p_num
				if(i==j) cycle L2
				if(p_inout(j)==0) cycle L2  !入口セルの場合
				if(p_couple_num(i)/=p_couple_num(j)) cycle L2   !Gate番号の一致の確認
!take160401                q_swap(p_indx(j,1,1),p_indx(j,2,1))=q_swap(p_indx(j,1,1),p_indx(j,2,1))+q_p !出口セル(流入)
                q_swap(p_indx(j,1,1),p_indx(j,2,1))=q_swap(p_indx(j,1,1),p_indx(j,2,1))+q_pin(i) !出口セル(流入)
				!出力用140602
!take160401                p_qswap(j)=q_p
                p_qswap(j)=q_pin(i)
			end do L2
		end do L1
        deallocate(q_pin) !take160401

		!出口側にポンプ排水量をデータで与える場合
L3:	do i=1,p_num
			q_p=0.0d0  !140602念のためゼロクリア
			if(p_inout(i)==0) cycle L3  !入口セルの場合
L4:		do j=1,p_num
				if(i==j) cycle L4
				if(p_couple_num(i)==p_couple_num(j)) cycle L3 !PUMP番号の一致の確認
			end do L4

			!PUMPが出口のみの場合
			if(time<=0.0d0) then
				q_p=p_qout(i,1)
			else if(time>=p_time(i,p_size(i))) then
				q_p=p_qout(i,p_size(i))
			else
				do k=1,p_size(i)-1
					if(time>=p_time(i,k).and.time<p_time(i,k+1)) then
						sst=(time-p_time(i,k))/(p_time(i,k+1)-p_time(i,k))
						q_p=p_qout(i,k)+(p_qout(i,k+1)-p_qout(i,k))*sst
					end if
				end do
			end if

			q_swap(p_indx(i,1,1),p_indx(i,2,1))=q_swap(p_indx(i,1,1),p_indx(i,2,1))+q_p !出口セル(流入)
			!出力用140602
			p_qswap(j)=q_p

		end do L3

	end if

end subroutine
!docon add end#######################################################################################
