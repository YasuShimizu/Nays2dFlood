!======================================================================
!------- output temporary file for hot start -----------
!======================================================================
subroutine output_hotstart(icount, output_file)
	use common_hh
	use variables
	!use output_cgn
	use GridCoord
	use GridCond
	!use avgeo_m
	!use gcoefs_m
	!use initl_m
	!use bound
	!use non_advec_1
	!use cal_gradient
	!use diffusion_m
	!use advection_m
	!use upstream_m
	!use allocate_variables
	!use interpolation
	!use transform
	!use box_gate_pump
	
	implicit none	
	integer, intent(in) :: icount
	character(len = strMax), intent(in) :: output_file
	
	integer :: i, j, ier
	
!======================================================================
   
	open(10,file=trim(output_file) ,form='unformatted', iostat=ier, err=100)
	!
	write(10) time,icount,dt
	write(10) nx,ny
	!
	write(10) ((x(i,j),i=0,nx),j=0,ny)
	write(10) ((y(i,j),i=0,nx),j=0,ny)
	write(10) ((hs(i,j),i=1,nx),j=1,ny)
	write(10) ((eta(i,j),i=1,nx),j=1,ny)
	write(10) ((h(i,j),i=1,nx),j=1,ny)
	write(10) ((yun(i,j),i=0,nx+1),j=0,ny+1)
	write(10) ((yvn(i,j),i=0,nx+1),j=0,ny+1)
	write(10) ((gux(i,j),i=0,nx+1),j=0,ny+1)
	write(10) ((guy(i,j),i=0,nx+1),j=0,ny+1)
	write(10) ((gvx(i,j),i=0,nx+1),j=0,ny+1)
	write(10) ((gvy(i,j),i=0,nx+1),j=0,ny+1)
	write(10) ((snu(i,j),i=1,nx),j=0,ny+1)
	write(10) ((snu_x(i,j),i=0,nx),j=-1,ny+1)
	close(10)
	
	return

100 write(*,*) "iostat = ", ier
   stop


end subroutine