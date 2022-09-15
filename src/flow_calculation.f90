! ---------------------------------------- !
!      —¬‚ê‚ÌŒvŽZ‚ÉŠÖ‚·‚éƒ‚ƒWƒ…[ƒ‹ŒQ      !
! ---------------------------------------- !

module bound

	use common_hh
	
	implicit none

  contains

! --------------------------------------------------------------------
	subroutine bound_h(h0,hs0,eta_o)
		use variables
		implicit none

		double precision,dimension(1:nx,0:ny+1),intent(in) :: eta_o
		double precision,dimension(1:nx,0:ny+1),intent(inout) :: h0,hs0

		integer i,j
		double precision :: dhdxi,dhdet
!
		if(j_wl==3) then
!			do i=nx,nx+1
!				do j=1,ny
!					hs0(i,j)=hs0(i-1,j)
!					h0(i,j)=eta_o(i,j)+hs0(i,j)
!				end do
!			end do

!$omp do private(dhdxi)
			do j=1,ny

!				if(hs0(nx,j)>hmin.and.h0(nx,j)<=h0(nx-1,j)) then
!				if(hs0(nx-1,j)>hmin2.and.h0(nx,j)<=h0(nx-1,j)) then
!					h0(nx,j)=h0(nx-1,j)
!					hs0(nx,j)=h0(nx,j)-eta_o(nx,j)
!				end if

				if( ijo_in(nx,j)==0 .or. hs0(nx,j)>=hmin2 ) then
!!					if( eta_o(nx,j)<eta_o(nx-1,j) ) then
!!						hs0(nx,j) = hs0(nx-1,j)
!!						h0(nx,j) = hs0(nx,j)+eta_o(nx,j)
!!					else
!						dhdxi = (h0(nx-1,j)-h0(nx-2,j))/dsy(nx-1,j)
!						if( dhdxi>0.d0 ) then
!							h0(nx,j) = hs0(nx-1,j)+eta_o(nx,j)
!						else
!							h0(nx,j) = h0(nx-1,j)+dhdxi*dsy(nx,j)
!						end if
!!!						h0(nx,j) = h0(nx-1,j)
!						hs0(nx,j) = h0(nx,j)-eta_o(nx,j)
!						if( hs0(nx,j)<=hmin ) hs0(nx,j) = hmin
!!					end if


					if( eta_o(nx,j)>h0(nx-1,j) ) then
						hs0(nx,j) = hmin
						h0(nx,j) = hs0(nx,j)+eta_o(nx,j)
					else
						if( eta_o(nx,j)<eta_o(nx-1,j) ) then
							hs0(nx,j) = hs0(nx-1,j)
							h0(nx,j) = hs0(nx,j)+eta_o(nx,j)
						else
							dhdxi = (h0(nx-1,j)-h0(nx-2,j))/dsy(nx-2,j)
							if( dhdxi>0.d0 ) dhdxi = 0.d0
							h0(nx,j) = h0(nx-1,j)+dhdxi*dsy(nx-1,j)
							hs0(nx,j) = h0(nx,j)-eta_o(nx,j)
							if( hs0(nx,j)<=hmin ) hs0(nx,j) = hmin
						end if
					end if

					
				else
					hs0(nx,j) = hmin
					h0(nx,j) = eta_o(nx,j)+hmin
				end if
				
			end do
		end if

!		do i=nx,nx+1
!			do j=1,ny
!				h(i,j)=h(nx-1,j)
!				hs(i,j)=h(i,j)-eta(i,j)
!			end do
!		end do
!
		if( j_side_j1==1 ) then
!$omp do private(dhdet)
			do i=1,nx
			
				if( eta_o(i,0)>h0(i,1) ) then
					hs0(i,0) = hmin
					h0(i,0) = hs0(i,0)+eta_o(i,0)
				else
					if( eta_o(i,0)<eta_o(i,1) ) then
						hs0(i,0) = hs0(i,1)
						h0(i,0) = hs0(i,0)+eta_o(i,0)
					else
						dhdet = (h0(i,2)-h0(i,1))/dnx(i,1)
						if( dhdet<0. ) dhdet = 0.d0
						h0(i,0) = h0(i,1)-dhdet*dnx(i,0)
						hs0(i,0) = h0(i,0)-eta_o(i,0)
						if( hs0(i,0)<=hmin ) hs0(i,0) = hmin
					end if
				end if
			
			end do
		end if
		
		if( j_side_j2==1 ) then
!$omp do private(dhdet)
			do i=1,nx
				
				if( eta_o(i,ny+1)>h0(i,ny) ) then
					hs0(i,ny+1) = hmin
					h0(i,ny+1) = hs0(i,ny+1)+eta_o(i,ny+1)
				else
					if( eta_o(i,ny+1)<eta_o(i,ny) ) then
						hs0(i,ny+1) = hs0(i,ny)
						h0(i,ny+1) = hs0(i,ny+1)+eta_o(i,ny+1)
					else
						dhdet = (h0(i,ny)-h0(i,ny-1))/dnx(i,ny-1)
						if( dhdet>0.d0 ) dhdet = 0.d0
						h0(i,ny+1) = h0(i,ny)+dhdet*dnx(i,ny)
						hs0(i,ny+1) = h0(i,ny+1)-eta_o(i,ny+1)
						if( hs0(i,ny+1)<=hmin ) hs0(i,ny+1) = hmin
					end if
				end if
	
			end do
		end if
!
		do j=1,ny
			do i=1,nx
				hs0(i,j) = hs0(i,j)*(1.d0-ijo_in(i,j))+hmin*ijo_in(i,j)
			end do
		end do

!		do j=1,ny
!			if(i<nx.and.ijo_in(i+1,j)==1) then
!				h0(i+1,j)=h0(i,j)
!				hs0(i+1,j)=h0(i+1,j)-eta_o(i+1,j)
!			else if(i>1.and.ijo_in(i-1,j)==1) then
!				h0(i-1,j)=h0(i,j)
!				hs0(i-1,j)=h0(i-1,j)-eta_o(i-1,j)
!			else if(j<ny.and.ijo_in(i,j+1)==1) then
!				h0(i,j+1)=h0(i,j)
!				hs0(i,j+1)=h0(i,j+1)-eta_o(i,j+1)
!			else if(j>1.and.ijo_in(i,j-1)==1) then
!				h0(i,j-1)=h0(i,j)
!				hs0(i,j-1)=h0(i,j-1)-eta_o(i,j-1)
!			end if
!		end do
!		end do

	end subroutine bound_h

! --------------------------------------------------------------------
	subroutine bound_u(u,ijobst_u)
		implicit none

		integer,dimension(0:nx,-1:ny+1),intent(in)::ijobst_u
		double precision,dimension(-1:nx+1,-1:ny+1),intent(inout)::u
		
		integer :: i,j
!
!$omp do
		do j=1,ny
			u(nx,j) = u(nx-1,j)
	!		if( u(nx,j)<0. ) u(nx,j) = 0.
		end do
!

		! ---- ‚·‚×‚Ä‚Ì‘¤•û‹«ŠEðŒ‚ÅƒXƒŠƒbƒv‚ð‰¼’è ---- (newgrd—p)

!$omp do
		do i=0,nx
			u(i,0) = u(i,1)
			u(i,ny+1) = u(i,ny)
		end do
!
!$omp do
		do j=1,ny
			do i=0,nx
					u(i,j) = u(i,j)*ijobst_u(i,j)
			end do
		end do

	end subroutine bound_u

! --------------------------------------------------------------------
	subroutine bound_qu
		use variables
		implicit none
		
		integer :: i,j,dzdxi,uu
		
!
!!$omp do private(j,dzdxi,uu)
!		do j=1,ny
!			if( ijo_in(nx,j)==0 ) then
!				dzdxi = (eta(nx-1,j)-eta(nx,j))/dsy(nx-1,j)
!				if( dzdxi<0.d0 ) dzdxi = 0.d0
!				uu = hs(nx,j)**(0.666667d0)*dzdxi**0.5d0/snmm(nx,j)
!				q_xi(nx,j) = uu*xi_r_up(nx,j)*hs(nx,j)/sj(nx,j)
!			end if
!		end do

!$omp do private(j)
		do j=1,ny
			q_xi(nx,j) = q_xi(nx-1,j)
		end do

!
!$omp do
		do j=1,ny
			do i=0,nx
					q_xi(i,j) = q_xi(i,j)*ijobst_u(i,j)
			end do
		end do

	end subroutine bound_qu

! --------------------------------------------------------------------
	subroutine bound_v(v,ijobst_v)
		use common_hh
		implicit none

		integer,dimension(0:nx,-1:ny+1),intent(in)::ijobst_v
		double precision,dimension(-1:nx+1,-1:ny+1),intent(inout)::v
		
		integer :: i,j

!$omp do
		do j=0,ny
			v(0,j) = v(1,j)
			v(nx,j) = v(nx-1,j)
			v(nx+1,j) = v(nx,j)
		end do

		if( j_side_j1==0 ) then
!$omp do
			do i=0,nx+1
				v(i,0) = 0.d0
			end do
		else if( j_side_j1==1 ) then
!$omp do
			do i=0,nx+1
				v(i,-1) = v(i,0)
			end do
		end if

		if( j_side_j2==0 ) then
!$omp do
			do i=0,nx+1
				v(i,ny) = 0.d0
			end do
		else if( j_side_j2==1 ) then
!$omp do
			do i=0,nx+1
				v(i,ny+1) = v(i,ny)
			end do
		end if

!$omp do
		do j=0,ny
			do i=1,nx
				v(i,j) = v(i,j)*ijobst_v(i,j)
			end do
		end do

	end subroutine bound_v

! --------------------------------------------------------------------
	subroutine bound_qv(v,ijobst_v)
		use common_hh
		implicit none

		integer,dimension(0:nx,-1:ny+1),intent(in)::ijobst_v
		double precision,dimension(0:nx,0:ny),intent(inout)::v
		
		integer :: i,j

		if( j_side_j1==0 ) then
!$omp do
			do i=0,nx
				v(i,0) = 0.d0
			end do
		end if

		if( j_side_j2==0 ) then
!$omp do
			do i=0,nx
				v(i,ny) = 0.d0
			end do
		end if

!$omp do
		do j=0,ny
			do i=1,nx
				v(i,j) = v(i,j)*ijobst_v(i,j)
			end do
		end do

	end subroutine bound_qv

! --------------------------------------------------------------------
	subroutine dryck_u(u,hs0,gx,gy)
		use common_hh
		implicit none
		double precision,dimension(1:nx,0:ny+1),intent(in) :: hs0
		double precision,dimension(-1:nx+1,-1:ny+1),intent(inout) :: u,gx,gy
		
		integer :: i,j
!
!$omp do
		do j=1,ny
			do i=1,nx-1
				if( hs0(i,j)<=hmin.and.hs0(i+1,j)<=hmin ) then
					u(i,j) = 0.
					gx(i,j) = 0.
					gy(i,j) = 0.
				else if( hs0(i,j)<=hmin.and.u(i,j)>0. ) then
					u(i,j) = 0.
					gx(i,j) = 0.
					gy(i,j) = 0.
				else if( hs0(i+1,j)<=hmin.and.u(i,j)<0. ) then
					u(i,j) = 0.
					gx(i,j) = 0.
					gy(i,j) = 0.
				end if
			end do
		end do

	end subroutine dryck_u

! --------------------------------------------------------------------
	subroutine dryck_v(v,hs0,gx,gy)
		use common_hh
		implicit none
		double precision,dimension(1:nx,0:ny+1),intent(in) :: hs0
		double precision,dimension(-1:nx+1,-1:ny+1),intent(inout) :: v,gx,gy
		
		integer :: i,j
!
!$omp do
		do j=1,ny-1
			do i=1,nx
				if( hs0(i,j)<=hmin.and.hs0(i,j+1)<=hmin ) then
					v(i,j) = 0.
					gx(i,j) = 0.
					gy(i,j) = 0.
				else if( hs0(i,j)<=hmin.and.v(i,j)>0. ) then
					v(i,j) = 0.
					gx(i,j) = 0.
					gy(i,j) = 0.
				else if( hs0(i,j+1)<=hmin.and.v(i,j)<0. ) then
					v(i,j) = 0.
					gx(i,j) = 0.
					gy(i,j) = 0.
				end if
			end do
		end do

	end subroutine dryck_v

end module bound

module non_advec_1
	use common_hh
	use transform
	use bound
	
	implicit none

  contains

! --------------------------------------------------------------------
	subroutine hcal(errmax,err,lcount,alh)
		use variables
		implicit none
		
		integer,intent(out) :: lcount
		double precision,intent(in) :: errmax, alh
		double precision,intent(out) :: err
		
		integer :: i,j,l
		double precision :: hs_up, v_up, ux_up, uy_up, vv_up	&
							, c_xi, c_xi_shear, c_veg, f_xi, dhdxi, dhdet, p_xi
		double precision :: hs_vp, u_vp, ux_vp, uy_vp, vv_vp	&
							, c_et, c_et_shear, f_et, p_et
		double precision :: div, hsta, serr
		
!$omp do
		do j=1,ny
			do i=1,nx
				whs(i,j) = hs(i,j)
			end do
		end do
! 
		do l=1,lmax

			call uvpcal(yun,yvn,up,vp,hs)
			call uxuycal(up,vp,ux,uy)
			
!$omp do private(hs_up,v_up,ux_up,uy_up,vv_up,c_xi,c_xi_shear,c_veg,f_xi,dhdxi,dhdet,p_xi)
			do j=1,ny
				do i=1,nx-1
					hs_up = (hs(i,j)+hs(i+1,j))*.5

					if( hs_up<=hmin ) then
						yun(i,j) = 0.d0
						q_xi(i,j) = 0.d0
					else if( hs(i+1,j)<=hmin2 .and. hn(i+1,j)>=hn(i,j) ) then
						yun(i,j) = 0.d0
						q_xi(i,j) = 0.d0
					else if( hs(i,j)<=hmin2 .and. hn(i,j)>=hn(i+1,j) ) then
						yun(i,j) = 0.d0
						q_xi(i,j) = 0.d0
					else
						v_up = (yvn(i,j-1)+yvn(i+1,j-1)+yvn(i,j)+yvn(i+1,j))*.25d0
						ux_up = (ux(i,j)+ux(i+1,j))*.5d0
						uy_up = (uy(i,j)+uy(i+1,j))*.5d0
						vv_up = dsqrt(ux_up**2d0+uy_up**2d0)
						c_xi = -(alpha(i,j,1)*yun(i,j)**2d0+alpha(i,j,2)*yun(i,j)*v_up+alpha(i,j,3)*v_up**2d0)
                        c_xi = c_xi*gam_e_up(i,j)	!Œš•¨‘jŠQ hamaki
						c_xi_shear = g*sn_up(i,j)**2d0/hs_up**1.33333
                        c_xi_shear = c_xi_shear*gam_v_up(i,j)	!Œš•¨‘jŠQ—¦@hamaki
!h						c_xi_shear = c_xi_shear + 0.5*0.383*(1-gam_v_up(i,j))	!Œš•¨‘jŠQ—¦@hamaki
						c_xi_shear = c_xi_shear + 0.5*cdd*(1-gam_v_up(i,j))	!Œš•¨‘jŠQ—¦@hamaki
						c_veg = (cd_veg(i,j)+cd_veg(i+1,j))*.5
						c_xi_shear = c_xi_shear+c_veg

						f_xi = - c_xi_shear*vv_up

						if( (eta(i,j)>eta(i+1,j).and.		&
								hs(i,j)<=hmin2.and.hn(i+1,j)<eta(i,j)) .or.	&
								(eta(i,j)<eta(i+1,j).and.	&
								hs(i+1,j)<=hmin2.and.hn(i,j)<eta(i+1,j)) ) then
							dhdxi = 0.
						else
							dhdxi = (-hn(i,j)+hn(i+1,j))*r_dxi*ijobst_u(i,j)
						end if
						
						if( j==1-j_side1 ) then
							if( hs(i+1,j+1)<=hmin.or.hs(i,j+1)<=hmin.or.		&
									hs(i+1,j)<=hmin.or.hs(i,j)<=hmin ) then
								dhdet = 0.
							else
								dhdet = (hn(i+1,j+1)+hn(i,j+1)-hn(i+1,j)-hn(i,j))*0.5d0*r_det
							end if
						else if( j==ny+j_side2 ) then
							if( hs(i+1,j)<=hmin.or.hs(i,j)<=hmin.or.		&
									hs(i+1,j-1)<=hmin.or.hs(i,j-1)<=hmin ) then
								dhdet = 0.
							else
								dhdet = (hn(i+1,j)+hn(i,j)-hn(i+1,j-1)-hn(i,j-1))*0.5d0*r_det
							end if
						else
							if( hs(i+1,j+1)<=hmin.or.hs(i,j+1)<=hmin.or.		&
									hs(i+1,j-1)<=hmin.or.hs(i,j-1)<=hmin ) then
								dhdet = 0.
							else
								dhdet = (hn(i+1,j+1)+hn(i,j+1)-hn(i+1,j-1)-hn(i,j-1))*0.25d0*r_det
							end if
						end if
						
						p_xi = -g*(beta(i,j,1)*dhdxi+beta(i,j,2)*dhdet)
						p_xi = p_xi*gam_v_up(i,j)	!Œš•¨‘jŠQ—¦ hamaki

						yun(i,j) = (yu(i,j)+(c_xi+p_xi)*dt)/(1.d0-f_xi*dt)

						if( hs(i,j)<=hmin2.and.yun(i,j)>0. ) yun(i,j) = 0.
						if( hs(i+1,j)<=hmin2.and.yun(i,j)<0. ) yun(i,j) = 0.
					end if
					
					q_xi(i,j) = yun(i,j)*hs_up*2.d0/(sj(i,j)+sj(i+1,j))
					q_xi(i,j) = q_xi(i,j)*gam_e_up(i,j)		!Œš•¨‘jŠQ—¦ hamaki
                end do
			end do
!
!$omp do private(hs_vp,u_vp,ux_vp,uy_vp,vv_vp,c_et,c_et_shear,c_veg,f_et,dhdxi,dhdet,p_et)
			do j=1-j_side1,ny-1+j_side2
				do i=1,nx-1
					hs_vp = (hs(i,j)+hs(i,j+1))*.5
					
					if( hs_vp<=hmin ) then
						yvn(i,j) = 0.
						q_et(i,j) = 0.
					else if( hs(i,j+1)<=hmin2.and.hn(i,j+1)>=hn(i,j) ) then
						yvn(i,j) = 0.
						q_et(i,j) = 0.
					else if( hs(i,j)<=hmin2.and.hn(i,j)>=hn(i,j+1) ) then
						yvn(i,j) = 0.
						q_et(i,j) = 0.
					else
						u_vp = (yun(i,j)+yun(i-1,j)+yun(i,j+1)+yun(i-1,j+1))*.25
						ux_vp = (ux(i,j)+ux(i,j+1))*.5
						uy_vp = (uy(i,j)+uy(i,j+1))*.5
						vv_vp = dsqrt(ux_vp**2d0+uy_vp**2d0)
						c_et = -(alpha(i,j,4)*u_vp**2d0+alpha(i,j,5)*u_vp*yvn(i,j)+alpha(i,j,6)*yvn(i,j)**2d0)
						c_et = c_et*gam_e_vp(i,j)	!Œš•¨‘jŠQ—¦ hamaki
                        c_et_shear = g*sn_vp(i,j)**2d0/hs_vp**1.33333
						c_et_shear = c_et_shear*gam_v_vp(i,j)	!Œš•¨‘jŠQ—¦ hamaki
!						c_et_shear = c_et_shear + 0.5*0.383*(1-gam_v_vp(i,j))	!Œš•¨‘jŠQ—¦ hamaki
						c_et_shear = c_et_shear + 0.5*cdd*(1-gam_v_vp(i,j))	!Œš•¨‘jŠQ—¦ hamaki
						c_veg = (cd_veg(i,j)+cd_veg(i,j+1))*.5
						c_et_shear = c_et_shear+c_veg

						f_et = -c_et_shear*vv_vp
!
						if( (eta(i,j)>eta(i,j+1).and.		&
								hs(i,j)<=hmin2.and.hn(i,j+1)<eta(i,j)) .or.		&
								(eta(i,j)<eta(i,j+1).and.		&
								hs(i,j+1)<=hmin2.and.hn(i,j)<eta(i,j+1)) ) then
							dhdet = 0.
						else
							dhdet = (-hn(i,j)+hn(i,j+1))*r_det*ijobst_v(i,j)
						end if
				
						if( i==1 ) then
							if( hs(i+1,j+1)<=hmin.or.hs(i+1,j)<=hmin.or.		&
									hs(i,j+1)<=hmin.or.hs(i,j)<=hmin ) then
								dhdxi = 0.
							else
								dhdxi = (hn(i+1,j+1)+hn(i+1,j)-hn(i,j+1)-hn(i,j))*0.5d0*r_dxi
							end if
						else if( i==nx ) then
							if( hs(i,j+1)<=hmin.or.hs(i,j)<=hmin.or.		&
									hs(i-1,j+1)<=hmin.or.hs(i-1,j)<=hmin ) then
								dhdxi = 0.
							else
								dhdxi = (hn(i,j+1)+hn(i,j)-hn(i-1,j+1)-hn(i-1,j))*0.5d0*r_dxi
							end if
						else
							if( hs(i+1,j+1)<=hmin.or.hs(i+1,j)<=hmin.or.		&
									hs(i-1,j+1)<=hmin.or.hs(i-1,j)<=hmin ) then
								dhdxi = 0.
							else
								dhdxi = (hn(i+1,j+1)+hn(i+1,j)-hn(i-1,j+1)-hn(i-1,j))*0.25d0*r_dxi
							end if
						end if
				
						p_et = -g*(beta(i,j,3)*dhdxi+beta(i,j,4)*dhdet)
						p_et = p_et*gam_v_vp(i,j)	!Œš•¨‘jŠQ—¦ hamaki
                        
						yvn(i,j) = (yv(i,j)+(c_et+p_et)*dt)/(1.d0-f_et*dt)
						
						if( hs(i,j)<=hmin2.and.yvn(i,j)>0. ) yvn(i,j) = 0.
						if( hs(i,j+1)<=hmin2.and.yvn(i,j)<0. ) yvn(i,j) = 0.
					end if
			
					q_et(i,j) = yvn(i,j)*hs_vp*2.d0/(sj(i,j)+sj(i,j+1))
					q_et(i,j) = q_et(i,j)*gam_e_vp(i,j)		!Œš•¨‘jŠQ—¦ hamaki
				end do
			end do
!
			call bound_u(yun,ijobst_u)
			call bound_qu
			call bound_v(yvn,ijobst_v)
			call bound_qv(q_et,ijobst_v)

!$omp single
			err=0.d0
!$omp end single

!$omp do reduction(+:err) private(div,hsta,serr)
			do j=1,ny
				do i=1,nx-1
					if( ijo_in(i,j)==1 ) goto 201
						div = (-q_xi(i-1,j)+q_xi(i,j))*r_dxi+(-q_et(i,j-1)+q_et(i,j))*r_det	
						if(j_qbl == 0) then	!hŒš•¨‚Ö‚ÌZ…‚ðl—¶‚µ‚È‚¢ê‡
							hsta                 = h(i,j)-div*dt*sj(i,j)/gam_v(i,j)+cell_discharge(i,j)*dt
							if( j_rain>=2 ) then
								if( hs(i,j)<=hmin .and. rain_t2(i,j)<=0.)	rain_t2(i,j)=0.	!h151225 NegativeRain‘Î‰ž‚Ì”O‚Ì‚½‚ßˆ—
								hsta             = h(i,j)-div*dt*sj(i,j)/gam_v(i,j)+cell_discharge(i,j)*dt+rain_t2(i,j)*dt
							end if
						else				!hŒš•¨‚Ö‚ÌZ…‚ðl—¶‚·‚éê‡
							hsta                 = h(i,j)-div*dt*sj(i,j)+cell_discharge(i,j)*dt
							if( j_rain>=2 ) then
								if( hs(i,j)<=hmin .and. rain_t2(i,j)<=0.)	rain_t2(i,j)=0.	!h151225 NegativeRain‘Î‰ž‚Ì”O‚Ì‚½‚ßˆ—
								hsta             = h(i,j)-div*dt*sj(i,j)+cell_discharge(i,j)*dt+rain_t2(i,j)*dt
							end if
						end if
						serr = dabs(hsta-hn(i,j))
						if( hs(i,j)>hmin ) err = err+serr
						hn(i,j) = hsta*alh+hn(i,j)*(1.-alh)
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
					201	continue
				end do
			end do
!
			call bound_h(hn,hs,eta)
!
			if( err<errmax ) exit
		end do

		lcount=l
!
	end	subroutine hcal
! --------------------------------------------------------------------

end module non_advec_1

module diffusion_m		! ----------------------------------------------------------------
	use common_hh
	use variables
	implicit none
        
	double precision,dimension(:,:),allocatable :: uvis_x, uvis_y, vvis_x, vvis_y
	double precision,dimension(:,:),allocatable :: cvis_x,cvis_y
     
  contains
      
! --------------------------------------------------------------------
	subroutine alloc_diffusion_temp_variables
		implicit none

		allocate( uvis_x(0:nx,0:ny+1), uvis_y(0:nx,0:ny+1) )
		allocate( vvis_x(0:nx,0:ny+1), vvis_y(0:nx,0:ny+1) )

		uvis_x = 0.; uvis_y=0.; vvis_x=0.; vvis_y=0.

	end subroutine alloc_diffusion_temp_variables
  
! --------------------------------------------------------------------
	subroutine diffusion(cw)
		implicit none
		double precision,intent(in) :: cw
		integer(4) :: i,j
		double precision :: hhx1,hhx2,uvis,vvis
!
! ------- uvix_x(i=1,nx,j=1,ny) ------
!
!$omp do
		do j=1,ny
			do i=1,nx
				uvis_x(i,j)=snu(i,j)*(xi_r(i,j-1)+xi_r(i,j))*.5d0*(-yun(i-1,j)+yun(i,j))*r_dxi
            end do
		end do
!
! ------- uvix_y(i=1,nx-1,j=0,ny) ------
!
!$omp do private(hhx1,hhx2)
		do j=1,ny-1
			do i=1,nx
				if( i==nx ) then
					hhx1 = hs(i,j)
					hhx2 = hs(i,j+1)
				else
					hhx1=(hs(i,j)+hs(i+1,j))*.5
					hhx2=(hs(i,j+1)+hs(i+1,j+1))*.5
				end if
				
				if(hhx1<=hmin.and.hhx2<=hmin) then
					uvis_y(i,j)=0.
				else
					uvis_y(i,j)=snu_x(i,j)*(-yun(i,j)+yun(i,j+1))*r_det*((et_r(i,j)+et_r(i,j+1))*.5)
				end if
			end do
		end do

!$omp do 
		do i=1,nx
			uvis_y(i, 0) = 0.d0
			uvis_y(i,ny) = 0.d0
		end do
!
			!	---- \‘¢•¨‚Ì–€ŽC‚Íl‚¦‚È‚¢ ----

!		do i=1,nx
!			do j=1,ny-1
!				hhx1=(hs(i,j)+hs(i+1,j))*.5
!				hhx2=(hs(i,j+1)+hs(i+1,j+1))*.5
				
!				if(hhx1<=hmin.or.(ijobst(i,j)==1.and.ijobst(i,j+1)==0))then
!					uvis_y(i,j)=cw*yun(i,j+1)*abs(yun(i,j+1))/(xi_r(i,j+1)+xi_r(i+1,j+1))*2.
!				end if
				
!				if(hhx2<=hmin.or.(ijobst(i,j)==1.and.ijobst(i,j-1)==0)) then
!					uvis_y(i,j)=-cw*yun(i,j)*abs(yun(i,j))/(xi_r(i,j-1)+xi_r(i+1,j-1))*2.
!				end if
!			end do
!
!			uvis_y(i,ny)=-cw*yun(i,ny)*abs(yun(i,ny))/(xi_r(i,ny-1)+xi_r(i+1,ny-1))*2.
!			uvis_y(i,0)=cw*yun(i,1)*abs(yun(i,1))/(xi_r(i,1)+xi_r(i+1,1))*2.
!		end do

!
!------- yun(i=1,nx-1 j=1,ny)------
!
!$omp do private(uvis)
		do j=1,ny
			do i=1,nx-1
				uvis=(-uvis_x(i,j)+uvis_x(i+1,j))*r_dxi*xi_r_up(i,j)	&
						+(-uvis_y(i,j-1)+uvis_y(i,j))*r_det*et_r(i,j)
				uvis=uvis*gam_e_up(i,j)		!Œš•¨‘jŠQ—¦ hamaki
                yun(i,j)=yun(i,j)+uvis*dt
			end do
		end do
!
! ------ vvis_x(i=1,nx-1 j=1,ny-1)
!
!$omp do private(hhx1,hhx2)
		do j=1-j_side1,ny-1+j_side2
			do i=1,nx-1
				hhx1=(hs(i,j)+hs(i,j+1))*.5
				hhx2=(hs(i+1,j)+hs(i+1,j+1))*.5
				
				if(hhx1<=hmin.and.hhx2<=hmin) then
					vvis_x(i,j)=0.
				else
					vvis_x(i,j)=snu_x(i,j)*(-yvn(i,j)+yvn(i+1,j))*r_dxi*(xi_r(i,j)+xi_r(i+1,j))*.5
				end if
			end do
		end do

			!	---- \‘¢•¨‚Ì–€ŽC‚Íl‚¦‚È‚¢ ----

!		do i=1,nx-1
!			do j=1-j_side1,ny-1+j_side2
!				hhx1=(hs(i,j)+hs(i,j+1))*.5
!				hhx2=(hs(i+1,j)+hs(i+1,j+1))*.5
!				
!				if(hhx1<=hmin.and.hhx2<=hmin) then
!					vvis_x(i,j)=0.
!				else if(hhx1<=hmin.or.(ijobst(i,j)==1.and.ijobst(i+1,j)==0)) then
!					vvis_x(i,j)=cw*yvn(i+1,j)*abs(yvn(i+1,j))/(et_r(i+1,j)+et_r(i+1,j+1))*2.
!				else if(hhx2<=hmin.or.(ijobst(i,j)==1.and.ijobst(i-1,j)==0)) then
!					vvis_x(i,j)=-cw*yvn(i,j)*abs(yvn(i,j))/(et_r(i-1,j)+et_r(i-1,j+1))*2.
!				else
!					vvis_x(i,j)=snu_x(i,j)*(yvn(i+1,j)-yvn(i,j))/dxi*(xi_r(i,j)+xi_r(i+1,j))*.5
!				end if
!			end do
!		end do
!
! ------ vvis_y(i=2,nx-1 j=1,ny)
!
!$omp do
		do j=1-j_side1,ny+j_side2
			do i=2,nx-1
				vvis_y(i,j)=snu(i,j)*(-yvn(i,j-1)+yvn(i,j))*r_det*(et_r(i-1,j)+et_r(i,j))*.5
			end do
		end do
!
! ------ yvn(i=2,nx-1 j=1,ny-1)
!
!$omp do private(vvis)
		do j=1-j_side1,ny-1+j_side2
			do i=2,nx-1
				vvis=(-vvis_x(i-1,j)+vvis_x(i,j))*r_dxi*xi_r(i,j)	&
					  +(-vvis_y(i,j)+vvis_y(i,j+1))*r_det*et_r_vp(i,j)
				vvis=vvis*gam_e_vp(i,j)		!Œš•¨‘jŠQ—¦ hamaki
                yvn(i,j)=yvn(i,j)+vvis*dt
			end do
		end do
!
	end subroutine diffusion
! --------------------------------------------------------------------

end module diffusion_m		! ------------------------------------------------------------

module cal_gradient

	use common_hh
	use variables

	implicit none
	
  contains

! --------------------------------------------------------------------
	subroutine newgrd_u
		implicit none
		
		integer :: i,j
!
!$omp do
		do j=1,ny
			do i=1,nx-1
				gux(i,j) = gux(i,j)+(-yun(i-1,j)+yun(i+1,j)+yu(i-1,j)-yu(i+1,j))*0.5d0*r_dxi
			end do
		end do

!$omp do
		do i=1,nx-1
			gux(i,0) = 0.d0
			gux(i,ny+1) = 0.d0
		end do

!$omp do
		do j=1-j_side1,ny+j_side2
			gux(0,j) = 0.d0
			gux(nx,j) = 0.d0
		end do
			
!$omp do
		do j=1,ny
			do i=0,nx
				guy(i,j)=(guy(i,j)+(-yun(i,j-1)+yun(i,j+1)	&
								+yu(i,j-1)-yu(i,j+1))*0.5d0*r_det)*ijobst_u(i,j)
			end do
		end do

!$omp do
		do i=0,nx
			guy(i,0) = 0.d0
			guy(i,ny+1) = 0.d0
		end do


	end subroutine newgrd_u

! --------------------------------------------------------------------
	subroutine newgrd_v
		implicit none

		integer :: i,j
!
!$omp do
		do j=1-j_side1,ny-1+j_side2
			do i=1,nx-1
					gvx(i,j)=(gvx(i,j)+(-yvn(i-1,j)+yvn(i+1,j)	&
									+yv(i-1,j)-yv(i+1,j))*0.5*r_dxi)*ijobst_v(i,j)
			end do
		end do

!$omp do
		do j=1-j_side1,ny-1+j_side2
			gvx(0,j) = 0.d0
			gvx(nx,j) = 0.d0
		end do

!$omp do
		do j=1-j_side1,ny-1+j_side2
			do i=1,nx
				gvy(i,j)=gvy(i,j)+(-yvn(i,j-1)+yvn(i,j+1)+yv(i,j-1)-yv(i,j+1))*0.5*r_det
			end do
		end do

		if( j_side_j1==0 ) then
!$omp do
			do i=1,nx
		!		gvy(i,0) = gvy(i,0)+(yvn(i,1)-yvn(i,0)-yv(i,1)+yv(i,0))*r_det
				gvy(i,0) = 0.d0
			end do
		else if( j_side_j1==1 ) then
!$omp do
			do i=1,nx
				gvy(i,-1) = 0.d0
			end do
		else if( j_side_j1==2 ) then
!$omp do
			do i=1,nx
				gvy(i,0) = 0.d0
			end do
		end if

		if( j_side_j2==0 ) then
!$omp do
			do i=1,nx
		!		gvy(i,ny) = gvy(i,ny)+(yvn(i,ny)-yvn(i,ny-1)-yv(i,ny)+yv(i,ny-1))*r_det
				gvy(i,ny) = 0.d0
			end do
		else if( j_side_j2==1 ) then
!$omp do
			do i=1,nx
				gvy(i,ny+1) = 0.d0
			end do
		else if( j_side_j2==2 ) then
!$omp do
			do i=1,nx
				gvy(i,ny) = 0.d0
			end do
		end if

	end subroutine newgrd_v

end module cal_gradient

module advection_m		! ----------------------------------------------------------------
	use common_hh
	use variables
	implicit none
        
    integer, parameter :: n_cip2up = 0
	double precision,dimension(:,:),allocatable :: fn,gxn,gyn,u,v
	
  contains
       
! --------------------------------------------------------------------
	subroutine alloc_advection_temp_variables
		implicit none

		allocate( fn(-1:nx+1,-1:ny+1) )
		allocate( gxn(-1:nx+1,-1:ny+1), gyn(-1:nx+1,-1:ny+1) )
		allocate( u(-1:nx+1,-1:ny+1), v(-1:nx+1,-1:ny+1) )

		fn = 0.d0;	gxn = 0.d0;	gyn = 0.d0
		u = 0.d0;	v = 0.d0

	end subroutine alloc_advection_temp_variables
!
! --------------------------------------------------------------------
	subroutine dcip2d_u(f,gx,gy)
		implicit none
		double precision,dimension(-1:nx+1,-1:ny+1),intent(inout) :: f,gx,gy
		integer :: i,j,isn,jsn,im1,jm1
		double precision ::  xx, yy 	&
				, fis, fjs, a1, b1, c1, d1, e1, f1, g1, gxo, gyo, tmp, tmq, u_dfdx,v_dfdy
!
!$omp do
		do j=1,ny
			do i=1,nx-1
				u(i,j) = yu(i,j)
				v(i,j) = 0.25d0*(yv(i,j-1)+yv(i+1,j-1)+yv(i,j)+yv(i+1,j))
			end do
		end do
!
!$omp do
		do j=1,ny
			u(0,j) = yu(0,j)
			v(0,j) = (yv(1,j)+yv(1,j-1))*0.5d0
			u(nx,j) = yu(nx,j)
			v(nx,j) = (yv(nx,j)+yv(nx,j-1))*0.5d0
		end do

!$omp do
		do i=0,nx
			u(i,0) = u(i,1)
			u(i,ny+1) = u(i,ny)
			v(i,0) = v(i,1)
			v(i,ny+1) = v(i,ny)
		end do
!
!$omp do private(xx,yy,isn,jsn,fis,fjs,im1,jm1,a1,b1,c1,d1,e1,f1,g1,tmp,tmq)
		do j=1,ny
			do i=1,nx-1-n_cip2up
				xx = -u(i,j)*dt
				yy = -v(i,j)*dt
				isn = sign(1.0,u(i,j)) 
				jsn = sign(1.0,v(i,j)) 
				fis = dble(isn) 
				fjs = dble(jsn)
				im1 = i-isn            
				jm1 = j-jsn            
				a1 = ((gx(im1,j)+gx(i,j))*dxi1*fis-2.0*(f(i,j)-f(im1,j)))/(dxi3*fis)
				e1 = (3.0*(f(im1,j)-f(i,j))+(gx(im1,j)+2.0*gx(i,j))*dxi1*fis)/dxi2
				b1 = ((gy(i,jm1)+gy(i,j))*det1*fjs-2.0*(f(i,j)-f(i,jm1)))/(det3*fjs)
				f1 = (3.0*(f(i,jm1)-f(i,j))+(gy(i,jm1)+2.0*gy(i,j))*det1*fjs)/det2
				tmp = f(i,j)-f(i,jm1)-f(im1,j)+f(im1,jm1)
				tmq = gy(im1,j)-gy(i,j)
				d1 = (-tmp -tmq*det1*fjs)/(dxi1*det2*fis)
				c1 = (-tmp-(gx(i,jm1)-gx(i,j))*dxi1*fis)/(dxi2*det1*fjs)
				g1 = (-tmq+c1*dxi2)/(dxi1*fis)

				fn(i,j) = ((a1*xx+c1*yy+e1)*xx+g1*yy+gx(i,j))*xx+((b1*yy+d1*xx+f1)*yy+gy(i,j))*yy+f(i,j)
				gxn(i,j) = (3.0*a1*xx+2.0*(c1*yy+e1))*xx+(d1*yy+g1)*yy+gx(i,j)
				gyn(i,j) = (3.0*b1*yy+2.0*(d1*xx+f1))*yy+(c1*xx+g1)*xx+gy(i,j)
			end do
		end do

!!$omp do private(u_dfdx,v_dfdy)
!		do j=1,ny
!			do i=nx-n_cip2up,nx-1
!				u_dfdx = ((u(i,j)+dabs(u(i,j)))*(-f(i-1,j)+f(i,j))		&
!							+(u(i,j)-dabs(u(i,j)))*(-f(i,j)+f(i+1,j)))*0.5d0*r_dxi
!				v_dfdy = ((v(i,j)+dabs(v(i,j)))*(-f(i,j-1)+f(i,j))		&
!							+(v(i,j)-dabs(v(i,j)))*(-f(i,j)+f(i,j+1)))*0.5d0*r_det
!				fn(i,j) = f(i,j)-(u_dfdx+v_dfdy)*dt
!			end do
!		end do


!$omp do private(gxo,gyo)
		do j=1,ny
			do i=1,nx-1
				f(i,j) = fn(i,j)
				gxo = (-yu(i-1,j)+yu(i+1,j))*.5*r_dxi
				gyo = (-yu(i,j-1)+yu(i,j+1))*.5*r_det
				gx(i,j)=gxn(i,j)-(gxo*(-u(i-1,j)+u(i+1,j))+gyo*(-v(i-1,j)+v(i+1,j)))*0.5*dt*r_dxi
				gy(i,j)=gyn(i,j)-(gxo*(-u(i,j-1)+u(i,j+1))+gyo*(-v(i,j-1)+v(i,j+1)))*0.5*dt*r_det
			end do
		end do
!
	end subroutine dcip2d_u
!
! --------------------------------------------------------------------
	subroutine upwind2d_u(f,gx,gy)
		implicit none
		double precision,dimension(-1:nx+1,-1:ny+1),intent(inout) :: f,gx,gy
		integer(4) :: i,j
		double precision :: u_dfdx, v_dfdy
!
!$omp do
		do j=1,ny
			do i=1,nx-1
				u(i,j) = yu(i,j)
				v(i,j) = 0.25d0*(yv(i,j-1)+yv(i+1,j-1)+yv(i,j)+yv(i+1,j))
			end do
		end do
!
!$omp do private(u_dfdx,v_dfdy)
		do j=1,ny
			do i=1,nx-1
				u_dfdx = ((u(i,j)+dabs(u(i,j)))*(-f(i-1,j)+f(i,j))		&
							+(u(i,j)-dabs(u(i,j)))*(-f(i,j)+f(i+1,j)))*0.5d0*r_dxi
				v_dfdy = ((v(i,j)+dabs(v(i,j)))*(-f(i,j-1)+f(i,j))		&
							+(v(i,j)-dabs(v(i,j)))*(-f(i,j)+f(i,j+1)))*0.5d0*r_det
!hamaki				fn(i,j) = f(i,j)-(u_dfdx+v_dfdy)*dt
				fn(i,j) = f(i,j)-(u_dfdx+v_dfdy)*dt*gam_e_up(i,j)/gam_v_up(i,j)
            end do
		end do

!$omp do 
		do j=1,ny
			do i=1,nx-1
				f(i,j) = fn(i,j)
			end do
		end do

	end subroutine upwind2d_u
!
! --------------------------------------------------------------------
	subroutine dcip2d_v(f,gx,gy)
		implicit none
		double precision,dimension(-1:nx+1,-1:ny+1),intent(inout) :: f,gx,gy
		integer(4) :: i,j,isn,jsn,im1,jm1
		double precision ::  xx, yy 	&
				, fis, fjs, a1, b1, c1, d1, e1, f1, g1, gxo, gyo, tmp, tmq, u_dfdx,v_dfdy
!
!$omp do
		do j=1-j_side1,ny-1+j_side2
			do i=0,nx+1
				v(i,j) = yv(i,j)
			end do
		end do

!$omp do
		do j=1-j_side1,ny-1+j_side2
			do i=1,nx
				u(i,j) = 0.25*(yu(i-1,j)+yu(i,j)+yu(i-1,j+1)+yu(i,j+1)) 
			end do
		end do

!$omp do
		do j=1-j_side1,ny-1+j_side2
			u(0,j) = (yu(0,j)+yu(0,j+1))*0.5d0
			u(nx+1,j) = (yu(nx,j)+yu(nx,j+1))*0.5d0
		end do

		if( j_side_j1==1 ) then
!$omp do
			do i=0,nx
				u(i,-1) = u(i,0)
				v(i,-1) = v(i,0)
			end do
		else
!$omp do
			do i=0,nx
				u(i,0) = (yu(i,1)+yu(i-1,1))*0.5d0
				v(i,0) = yv(i,0)
			end do
		end if

		if( j_side_j2==1 ) then
!$omp do
			do i=0,nx
				u(i,ny+1) = u(i,ny)
				v(i,ny+1) = v(i,ny)
			end do
		else
!$omp do
			do i=0,nx
				u(i,ny) = (yu(i,ny)+yu(i-1,ny))*0.5d0
				v(i,ny) = yv(i,ny)
			end do
		end if
!
!$omp do private(xx,yy,isn,jsn,fis,fjs,im1,jm1,a1,b1,c1,d1,e1,f1,g1,tmp,tmq)
		do j=1-j_side1,ny-1+j_side2
			do i=1,nx-1-n_cip2up
				xx = -u(i,j)*dt
				yy = -v(i,j)*dt
				isn = sign(1.0,u(i,j)) 
				jsn = sign(1.0,v(i,j)) 
				fis = dble(isn) 
				fjs = dble(jsn)
				im1 = i-isn            
				jm1 = j-jsn            
				a1 = ((gx(im1,j)+gx(i,j))*dxi1*fis		&
						-2.0*(f(i,j)-f(im1,j)))/(dxi3*fis)
				e1 = (3.0*(f(im1,j)-f(i,j))			&
						+(gx(im1,j)+2.0*gx(i,j))*dxi1*fis)/dxi2
				b1 = ((gy(i,jm1)+gy(i,j))*det1*fjs		&
						-2.0*(f(i,j)-f(i,jm1)))/(det3*fjs)
				f1 = (3.0*(f(i,jm1)-f(i,j))			&
						+(gy(i,jm1)+2.0*gy(i,j))*det1*fjs)/det2
				tmp = f(i,j)-f(i,jm1)-f(im1,j)+f(im1,jm1)
				tmq = gy(im1,j)-gy(i,j)
				d1 = (-tmp-tmq*det1*fjs)/(dxi1*det2*fis)
				c1 = (-tmp-(gx(i,jm1)-gx(i,j))*dxi1*fis)/(dxi2*det1*fjs)
				g1 = (-tmq+c1*dxi2)/(dxi1*fis)

				fn(i,j) = ((a1*xx+c1*yy+e1)*xx+g1*yy+gx(i,j))*xx		&
							+((b1*yy+d1*xx+f1)*yy+gy(i,j))*yy+f(i,j)
				gxn(i,j) = (3.0*a1*xx+2.0*(c1*yy+e1))*xx				&
							+(d1*yy+g1)*yy+gx(i,j)
				gyn(i,j) = (3.0*b1*yy+2.0*(d1*xx+f1))*yy				&
							+(c1*xx+g1)*xx+gy(i,j)
			end do
		end do
!
!!$omp do private(u_dfdx,v_dfdy)
!		do j=1-j_side1,ny-1+j_side2
!			do i=nx-n_cip2up,nx-1
!				u_dfdx = ((u(i,j)+dabs(u(i,j)))*(-f(i-1,j)+f(i,j))		&
!							+(u(i,j)-dabs(u(i,j)))*(-f(i,j)+f(i+1,j)))*0.5d0*r_dxi
!				v_dfdy = ((v(i,j)+dabs(v(i,j)))*(-f(i,j-1)+f(i,j))		&
!							+(v(i,j)-dabs(v(i,j)))*(-f(i,j)+f(i,j+1)))*0.5d0*r_det
!				fn(i,j) = f(i,j)-(u_dfdx+v_dfdy)*dt
!			end do
!		end do

!$omp do private(gxo,gyo)
		do j=1,ny-1
			do i=1,nx-1
				f(i,j) = fn(i,j)
				gxo = (-yv(i-1,j)+yv(i+1,j))*.5*r_dxi
				gyo = (-yv(i,j-1)+yv(i,j+1))*.5*r_det
				gx(i,j) = gxn(i,j)-(gxo*(-u(i-1,j)+u(i+1,j))		&
							+gyo*(-v(i-1,j)+v(i+1,j)))*0.5*dt*r_dxi
				gy(i,j) = gyn(i,j)-(gxo*(-u(i,j-1)+u(i,j+1))		&
							+gyo*(-v(i,j-1)+v(i,j+1)))*0.5*dt*r_det
			end do
		end do

	end subroutine dcip2d_v
!
! --------------------------------------------------------------------
	subroutine upwind2d_v(f,gx,gy)
		implicit none
		double precision,dimension(-1:nx+1,-1:ny+1),intent(inout) :: f,gx,gy
		integer(4) :: i,j
		double precision :: u_dfdx, v_dfdy
!
!$omp do 
		do j=1-j_side1,ny-1+j_side2
			do i=1,nx
				u(i,j) = 0.25*(yu(i-1,j)+yu(i,j)+yu(i-1,j+1)+yu(i,j+1)) 
				v(i,j) = yv(i,j)
			end do
		end do
!
!$omp do private(u_dfdx,v_dfdy)
		do j=1-j_side1,ny-1+j_side2
			do i=1,nx-1
				u_dfdx = ((u(i,j)+dabs(u(i,j)))*(-f(i-1,j)+f(i,j))		&
							+(u(i,j)-dabs(u(i,j)))*(-f(i,j)+f(i+1,j)))*0.5d0*r_dxi
				v_dfdy = ((v(i,j)+dabs(v(i,j)))*(-f(i,j-1)+f(i,j))		&
							+(v(i,j)-dabs(v(i,j)))*(-f(i,j)+f(i,j+1)))*0.5d0*r_det
!hamaki				fn(i,j) = f(i,j)-(u_dfdx+v_dfdy)*dt
				fn(i,j) = f(i,j)-(u_dfdx+v_dfdy)*dt*gam_e_vp(i,j)/gam_v_vp(i,j)
            end do
		end do

!$omp do
		do j=1-j_side1,ny-1+j_side2
			do i=1,nx-1
				f(i,j) = fn(i,j)
			end do
		end do

	end subroutine upwind2d_v
!
! --------------------------------------------------------------------
       
end module advection_m		! ------------------------------------------------------------
!
