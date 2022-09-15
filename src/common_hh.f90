module common_hh		! ---------------------------------------------------------------
!
	double precision :: time
	integer :: nx,ny,nym,nobst,lmax							&
				,jt1,jt2,nk,nm,iscd,j_in,i_in			&
				,j_side_j1,j_side_j2,jsin1,jsin2,j_side1,j_side2	&
				,j_rain,j_wl,j_drg,n_q_cell,i_bomb
	double precision :: dxi,det,dt,g,hmin,hmin2	&
				,tuk,ster,topi,pi,rain_t,diam,chl,width		&
				,dxi1,dxi2,dxi3,det1,det2,det3,bomber_time

	double precision :: r_dxi, r_det
	integer,parameter :: strMax=250
	character(len=strMax) :: condFile,gridFile,outFile		&
							,qhFile,mixFile,outFile_ini_mix,outFile_mix	&
							,j_inputFile,js1inputFile,js2inputFile		&
							,qj1file,qj2file
	integer :: cal_t1,cal_t2,t_rate,t_max,t_diff	!h160105ŒvŽZŽžŠÔŒv‘ª
!
end module common_hh	! ----------------------------------------------------------------

module variables		! ----------------------------------------------------------------
!
	double precision,dimension(:,:),allocatable :: yu,yun,yv,yvn
	double precision,dimension(:,:),allocatable :: yc,ycn,ycb
	double precision,dimension(:,:),allocatable :: up,vp,ux,uy
	double precision,dimension(:),allocatable :: qc
    integer,dimension(:,:),allocatable :: ibc, ibcg
	double precision,dimension(:,:),allocatable :: h,hn,hs,whs
	double precision,dimension(:,:),allocatable :: gux,guy,gvx,gvy
	double precision,dimension(:,:),allocatable :: x,y,z,zz,eta,ds,dn,sj,r_sj	&
								,xi_x,xi_y,et_x,et_y,xi_r,et_r,xi_r_up,et_r_vp,z0	&
								,dsn,area
!docon add start#####################################################################################
	double precision,dimension(:,:),allocatable :: q_swap
!docon add end#######################################################################################
	double precision,dimension(:,:),allocatable :: xi_x_up,xi_y_up			&
								,et_x_up,et_y_up,xi_x_vp,xi_y_vp,et_x_vp,et_y_vp
	double precision,dimension(:,:),allocatable :: vti,tausta,qsu,usta
	double precision,dimension(:,:),allocatable :: uxx,uyy
	double precision,dimension(:),allocatable :: eave,chb,emin,emax
	double precision,dimension(:,:),allocatable :: eta0
	double precision,dimension(:,:),allocatable :: cos_t,sin_t
	integer,dimension(:,:),allocatable :: ijobst,ijo_in,ij_ero,ijobst_u,ijobst_v
	double precision,dimension(:,:),allocatable :: snmm,sn_up,sn_vp
	double precision,dimension(:,:),allocatable :: xi_x0,xi_y0		&
								,et_x0,et_y0,sj0,hsxx,r_sj0
	double precision,dimension(:),allocatable :: t_hyd,q_ups,h_ups,h_dse,rain
	double precision,dimension(:,:),allocatable :: cf,re,cd_veg
	double precision,dimension(:,:),allocatable :: snu,snu_x
	double precision,dimension(:,:,:),allocatable :: alpha
	double precision,dimension(:,:,:),allocatable :: beta
	double precision,dimension(:,:),allocatable :: dnx, dsy
	double precision,dimension(:,:),allocatable :: qu,qv
	double precision,dimension(:,:),allocatable :: q_xi,q_et
	double precision,dimension(:,:),allocatable :: dex
	double precision,dimension(:,:),allocatable :: q_ups_in, h_ups_in
	double precision,dimension(:),allocatable :: q_input_in, h_input_in
	double precision,dimension(:),allocatable :: upv_slope_in
	integer,dimension(:),allocatable :: js,je
	double precision,dimension(:),allocatable :: q_in_j1,h_in_j1
	double precision,dimension(:),allocatable :: vpv1_slope_in
	integer,dimension(:),allocatable :: j1s,j1e
	double precision,dimension(:,:),allocatable :: q1_ups_in, h1_ups_in
	double precision,dimension(:),allocatable :: q_in_j2,h_in_j2
	double precision,dimension(:),allocatable :: vpv2_slope_in
	integer,dimension(:),allocatable :: j2s,j2e
	double precision,dimension(:,:),allocatable :: q2_ups_in, h2_ups_in
	integer,dimension(:),allocatable :: j_inlen, j_size, k_size
	integer, dimension(:,:,:), allocatable:: indices
	double precision,dimension(:),allocatable :: slopevalue
	double precision, dimension(:,:), allocatable:: f_param
	double precision, dimension(:,:), allocatable:: f_value
	double precision,dimension(:,:),allocatable :: cell_discharge
	
	integer,dimension(:,:),allocatable :: lcb

!inoue Œš•¨hamaki ver
   	double precision,dimension(:,:),allocatable :: gam_v,gam_e,share
  	double precision,dimension(:,:),allocatable :: gam_v_up,gam_v_vp,gam_e_up,gam_e_vp
	integer :: j_qbl
	double precision :: cdd
!h160105 Xrain‘Î‰ž
	double precision,dimension(:,:),allocatable :: rain_t2,rain_xx
	double precision,dimension(:),allocatable :: t_rain	
	integer :: nr
	!
    character(250),dimension(:),allocatable :: flowname, fname_in, fname_1, fname_2 
    
end module variables	! --------------------------------------------------------------------

module allocate_variables
	implicit none
	
  contains
!
! --------------------------------------------------------------------
	subroutine alloc_var1
		use common_hh
		use variables
		implicit none
!
		allocate( qc(0:nx),eave(0:nx),chb(0:nx),emin(0:nx),emax(0:nx) )
		
			qc = 0.d0;	eave = 0.d0;	chb = 0.d0;	emin = 0.d0;	emax = 0.d0
!
		allocate( alpha(1:nx,0:ny+1,6),beta(1:nx,0:ny+1,4) )
		
			alpha = 0.d0;	beta = 0.d0
		
		allocate( dnx(0:nx,0:ny),dsy(0:nx,0:ny) )
		
		allocate( qu(0:nx,0:ny),qv(0:nx,0:ny)	&
				,q_xi(0:nx,0:ny), q_et(0:nx,0:ny) )

		allocate( x(0:nx,0:ny),y(0:nx,0:ny),z(0:nx,0:ny)	&
					,zz(0:nx,0:ny),z0(0:nx,0:ny)	&
					,uxx(0:nx,0:ny),uyy(0:nx,0:ny)	&
					,xi_x0(0:nx,0:ny),xi_y0(0:nx,0:ny)	&
					,et_x0(0:nx,0:ny),et_y0(0:nx,0:ny)	&
					,sj0(0:nx,0:ny),r_sj0(0:nx,0:ny)	&
					,hsxx(0:nx,0:ny),ijobst(0:nx,0:ny),ibcg(0:nx,0:ny) )
					
			dnx = 0.d0;	dsy = 0.d0;	qu = 0.d0;	qv = 0.d0;	q_xi = 0.d0;
			q_et = 0.d0;	x = 0.d0;	y = 0.d0;	z = 0.d0;	zz = 0.d0;	z0 = 0.d0;
			uxx = 0.d0;	uyy = 0.d0;	xi_x0 = 0.d0;	xi_y0 = 0.d0;
			et_x0 = 0.d0;	et_y0 = 0.d0;	sj0 = 0.d0;	r_sj0 = 0.d0;
			hsxx = 0.d0;
!		
		allocate( yu(-1:nx+1,-1:ny+1),yun(-1:nx+1,-1:ny+1)	&
					,yv(-1:nx+1,-1:ny+1),yvn(-1:nx+1,-1:ny+1) )

		allocate( gux(-1:nx+1,-1:ny+1),guy(-1:nx+1,-1:ny+1)	&
					,gvx(-1:nx+1,-1:ny+1),gvy(-1:nx+1,-1:ny+1) )
					
			yu = 0.d0;	yun = 0.d0;	yv = 0.d0;	yvn = 0.d0;
			gux = 0.d0;	guy = 0.d0;	gvx = 0.d0;	gvy = 0.d0;
					
		allocate( up(1:nx,0:ny+1),vp(1:nx,0:ny+1)	&
					,ux(1:nx,0:ny+1),uy(1:nx,0:ny+1) )

        
		allocate( hn(1:nx,0:ny+1),h(1:nx,0:ny+1)	&
					,hs(1:nx,0:ny+1),whs(1:nx,0:ny+1)	&
					,eta(1:nx,0:ny+1),eta0(1:nx,0:ny+1)	&
					,ijo_in(1:nx,0:ny+1),ij_ero(1:nx,0:ny+1)	&
					,cd_veg(1:nx,0:ny+1),snmm(1:nx,0:ny+1),snu(1:nx,0:ny+1), ibc(1:nx,0:ny+1) )
!docon add start#####################################################################################
		allocate( q_swap(1:nx,0:ny+1))
!docon add end#######################################################################################
					
			up = 0.d0;	vp = 0.d0;	ux = 0.d0;	uy = 0.d0;
			hn = 0.d0;	h = 0.d0;	hs = 0.d0;	whs = 0.d0;
			eta = 0.d0;	eta0 = 0.d0;	cd_veg = 0.d0;
			snmm = 0.d0;	snu = 0.d0;
			ijo_in = 0;	
					
		allocate( ds(0:nx,-1:ny+1),dn(0:nx,-1:ny+1)	&
					,sj(0:nx,-1:ny+1),r_sj(0:nx,-1:ny+1)	&
					,xi_x(0:nx,-1:ny+1),xi_y(0:nx,-1:ny+1)	&
					,et_x(0:nx,-1:ny+1),et_y(0:nx,-1:ny+1)	&
					,xi_r(0:nx,-1:ny+1),et_r(0:nx,-1:ny+1)	&
					,xi_r_up(0:nx,-1:ny+1),et_r_vp(0:nx,-1:ny+1)	&
					,xi_x_up(0:nx,-1:ny+1),xi_y_up(0:nx,-1:ny+1)	&
					,et_x_up(0:nx,-1:ny+1),et_y_up(0:nx,-1:ny+1)	&
					,xi_x_vp(0:nx,-1:ny+1),xi_y_vp(0:nx,-1:ny+1)	&
					,et_x_vp(0:nx,-1:ny+1),et_y_vp(0:nx,-1:ny+1)	&
					,ijobst_u(0:nx,-1:ny+1),ijobst_v(0:nx,-1:ny+1)	&
					,sn_up(0:nx,-1:ny+1),sn_vp(0:nx,-1:ny+1),snu_x(0:nx,-1:ny+1) )
					
			ds = 0.d0;	dn = 0.d0;	sj = 0.d0;	r_sj = 0.d0;	xi_x = 0.d0;
			xi_y = 0.d0;	et_x = 0.d0;	et_y = 0.d0;	xi_r = 0.d0;	
			et_r = 0.d0;	xi_r_up = 0.d0;	et_r_vp = 0.d0;	xi_x_up = 0.d0;	
			xi_y_up = 0.d0;	et_x_up = 0.d0;	et_y_up = 0.d0;	xi_x_vp = 0.d0;	
			xi_y_vp = 0.d0;	et_x_vp = 0.d0;	et_y_vp = 0.d0;	sn_up = 0.d0;	
			sn_vp = 0.d0;	snu_x = 0.d0;	ijobst_u = 0;	ijobst_v = 0;
!
		allocate( vti(1:nx,1:ny),usta(1:nx,1:ny),cos_t(1:nx,1:ny)	&
					,sin_t(1:nx,1:ny),cf(1:nx,1:ny),re(1:nx,1:ny) )
					
			vti = 0.d0;	usta = 0.d0;
			cos_t = 0.d0;	sin_t = 0.d0;	cf = 0.d0;	re = 0.d0;

		allocate( cell_discharge(1:nx,1:ny) )
		
			cell_discharge = 0.d0
			
		allocate( dsn(1:nx,1:ny), area(1:nx,1:ny) )
		
		allocate( lcb(1:nx,1:ny) )
		
		lcb = 1

!inoue Œš•¨hamaki ver
        allocate( gam_v(1:nx,0:ny+1),gam_e(1:nx,0:ny+1),share(1:nx,0:ny+1) )
        gam_v = 0.d0; gam_e = 0.0d0; share = 0.0d0;
        allocate( gam_v_up(0:nx,-1:ny+1),gam_v_vp(0:nx,-1:ny+1) &
                    ,gam_e_up(0:nx,-1:ny+1),gam_e_vp(0:nx,-1:ny+1) )
        gam_v_up = 0.d0; gam_v_vp = 0.0d0;
		gam_e_up = 0.d0; gam_e_vp = 0.0d0;
!h160105 Xrain‘Î‰ž
        allocate( rain_t2(1:nx,1:ny),rain_xx(0:nx,0:ny) )
        rain_t2 = 0.0d0; rain_xx = 0.0d0;
!
	end subroutine alloc_var1

end module allocate_variables

!docon add start#####################################################################################
module box_gate_pump
	implicit none
	integer :: bc_num,gt_num,p_num,gt_tmax,gt_tmax2,p_tmax
	integer,dimension(:),allocatable :: bc_couple_num,gt_couple_num,p_couple_num
	integer,dimension(:),allocatable :: bc_inout,bc_c2_yesno,gt_inout,gt_c2_yesno,gt_method,gt_outwl
	integer,dimension(:),allocatable :: p_inout,p_method,gt_size,gt_size2,p_size
	integer,dimension(:),allocatable :: bc_mxindx,gt_mxindx,p_mxindx
	integer,dimension(:,:,:),allocatable :: bc_indx,gt_indx,p_indx
	double precision,dimension(:),allocatable :: bc_width,bc_height,bc_base,bc_ren,bc_c1,bc_c2,bc_c3
	double precision,dimension(:),allocatable :: gt_width,gt_height,gt_base,gt_ren,gt_c1,gt_c2,gt_c3
	double precision,dimension(:),allocatable :: p_qmax,p_startope,p_stopope,p_qt
	double precision,dimension(:,:),allocatable :: gt_time,gt_state,gt_time2,gt_wl,p_time,p_qout
	double precision,dimension(:),allocatable :: bc_qswap,gt_qswap,p_qswap  !140602
	character(250),dimension(:),allocatable :: BCname,Gtname,Ppname !140602
end module box_gate_pump
!docon add end#######################################################################################

