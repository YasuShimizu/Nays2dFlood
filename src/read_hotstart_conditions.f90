subroutine read_hotstart_conditions(icount, tmp_file_i)
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
	integer, intent(inout) :: icount
	character(len = strMax), intent(in) :: tmp_file_i

	integer :: i, j
	integer :: is, ii
	integer :: nx2, ny2
	real(8) :: dt2
	
	integer :: ierr_tmp

!======================================================================	
	
	open(11,file=tmp_file_i,status='old',iostat = ierr_tmp,form='unformatted')
	if(ierr_tmp /= 0) then
		write(6,*) 'Input file error!'
		write(6,*) 'Temporary file for Hot Start does not exist'
		pause
		stop
	end if
	!
	read(11) time,icount,dt2
	!
	time=(icount-1)*dt2
	icount=time/dt
	!
	read(11) nx2, ny2
	if(nx /= nx2.or.ny /= ny2) then
		write(6,*) 'Number of grid is different between grid file and temporary file!'
		pause
		stop
	end if
	!
	read(11) ((x(i,j),i=0,nx2),j=0,ny2)
	read(11) ((y(i,j),i=0,nx2),j=0,ny2)
	read(11) ((hs(i,j),i=1,nx2),j=1,ny2)
	read(11) ((eta(i,j),i=1,nx2),j=1,ny2)
	read(11) ((h(i,j),i=1,nx2),j=1,ny2)
	read(11) ((yun(i,j),i=0,nx2+1),j=0,ny2+1)
	read(11) ((yvn(i,j),i=0,nx2+1),j=0,ny2+1)
	read(11) ((gux(i,j),i=0,nx2+1),j=0,ny2+1)
	read(11) ((guy(i,j),i=0,nx2+1),j=0,ny2+1)
	read(11) ((gvx(i,j),i=0,nx2+1),j=0,ny2+1)
	read(11) ((gvy(i,j),i=0,nx2+1),j=0,ny2+1)
	read(11) ((snu(i,j),i=1,nx2),j=0,ny2+1)
	read(11) ((snu_x(i,j),i=0,nx2),j=-1,ny2+1)

	close(11)
	
	!
	do j=1,ny2
		do i=1,nx2
			hn(i,j) = h(i,j)
		end do
	end do
	!
	do j=0,ny2+1
		do i=0,nx2+1
			yu(i,j) = yun(i,j)
			yv(i,j) = yvn(i,j)
		end do
	end do

	return
	
end subroutine
