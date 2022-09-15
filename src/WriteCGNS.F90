      module output_cgn
!=======================================================================
      contains

	SUBROUTINE Write_CGNS(InputFile,Index,time,disch1,disch2,disch3,fname_in,fname_1,fname_2,j_in,jsin1,jsin2,j_side_j1,j_side_j2,im,jm	&
							,x,y,u,v,hs,z,z0,iofrg,h_max,v_max,j_rain,rain, ibcg0)
!docon add start#####################################################################################
        use box_gate_pump   !140602
!docon add end#######################################################################################
		IMPLICIT NONE
		INCLUDE "cgnslib_f.h"

		CHARACTER(*), INTENT(IN) :: InputFile
		double precision, INTENT(IN) :: time
		double precision, INTENT(IN) :: disch1(1:j_in)
		double precision, INTENT(IN) :: disch2(1:jsin1)
		double precision, INTENT(IN) :: disch3(1:jsin2)
        character(250), INTENT(IN) :: fname_in(1:j_in)
        character(250), INTENT(IN) :: fname_1(1:jsin1)
        character(250), INTENT(IN) :: fname_2(1:jsin2)
		INTEGER, INTENT(IN) :: Index
		INTEGER, INTENT(IN) :: iofrg
		INTEGER, INTENT(IN) :: im, jm
		INTEGER, INTENT(IN) :: j_side_j1, j_side_j2,j_in,jsin1,jsin2
        integer, dimension(0:im,0:jm),intent(in) :: ibcg0
		double precision,dimension(0:im,0:jm),intent(in) :: u,v,hs,z,z0
		double precision,dimension(0:im,0:jm),intent(in) :: h_max,v_max
		double precision,dimension(0:im,0:jm),intent(in) :: x, y
!h160123 Xrain‘Î‰ž
		INTEGER, INTENT(IN) :: j_rain
		double precision,dimension(0:im,0:jm),intent(in) :: rain

		INTEGER :: NX, NY
		INTEGER :: FID, BID, ZID, IER, iret
		INTEGER :: I,J
		double precision, ALLOCATABLE, DIMENSION(:,:) :: Vdata1, Udata1, Hdata1, Zbdata1, WSE
		double precision, ALLOCATABLE, DIMENSION(:,:) :: Vdata2, Hdata2
		double precision, ALLOCATABLE, DIMENSION(:,:) :: xx, yy
		double precision, ALLOCATABLE, DIMENSION(:,:) :: Rdata	!h160123 Xrain‘Î‰ž
        integer, ALLOCATABLE, DIMENSION(:,:) :: ibc0	!h160123 Xrain‘Î‰ž

		nx = im+1
		ny = jm+1

		CALL CG_IRIC_WRITE_SOL_TIME_F(time, IER)

		do i = 1, j_in
			CALL CG_IRIC_WRITE_SOL_BASEITERATIVE_REAL_F(trim(fname_in(i)), disch1(i), IER)
        end do
        
        if(j_side_j1.eq.2) then
		do i = 1, jsin1
		    CALL CG_IRIC_WRITE_SOL_BASEITERATIVE_REAL_F(trim(fname_1(i)), disch2(i), IER)
		end do
        end if
		
        if(j_side_j2.eq.2) then
		do i = 1, jsin2
!docon edit start#####################################################################################
!			CALL CG_IRIC_WRITE_SOL_BASEITERATIVE_REAL_F(trim(fname_2(i)), disch3(1), IER)
			CALL CG_IRIC_WRITE_SOL_BASEITERATIVE_REAL_F(trim(fname_2(i)), disch3(i), IER)   !140602
!docon edit end#######################################################################################
		end do
		end if

!docon add start#####################################################################################
        if(bc_num.ge.1) then    !140602
		do i = 1, bc_num
			CALL CG_IRIC_WRITE_SOL_BASEITERATIVE_REAL_F(trim(Bcname(i)), bc_qswap(i), IER)
		end do
		end if

        if(gt_num.ge.1) then
		do i = 1, gt_num
			CALL CG_IRIC_WRITE_SOL_BASEITERATIVE_REAL_F(trim(Gtname(i)), gt_qswap(i), IER)
		end do
		end if

        if(p_num.ge.1) then
		do i = 1, p_num
			CALL CG_IRIC_WRITE_SOL_BASEITERATIVE_REAL_F(trim(Ppname(i)), p_qswap(i), IER)
		end do
		end if
!docon add end#######################################################################################

		ALLOCATE(xx(nx, ny), STAT=ier)
		ALLOCATE(yy(nx, ny), STAT=ier)
		ALLOCATE(Udata1(nx, ny), STAT=ier)
		ALLOCATE(Vdata1(nx, ny), STAT=ier)
		ALLOCATE(Vdata2(nx, ny), STAT=ier)
		ALLOCATE(Hdata1(nx, ny), STAT=ier)
		ALLOCATE(Hdata2(nx, ny), STAT=ier)
		ALLOCATE(Zbdata1(nx, ny), STAT=ier)
		ALLOCATE(WSE(nx, ny), STAT=ier)
		ALLOCATE(Rdata(nx, ny), STAT=ier)	!h160123 Xrain‘Î‰ž
        ALLOCATE(ibc0(nx, ny), STAT=ier)

		DO j=1,NY
			DO i=1,NX
				xx(i,j) = x(i-1,j-1)
				yy(i,j) = y(i-1,j-1)
				Zbdata1(i,j) = z(i-1,j-1)
				WSE(i,j) = z(i-1, j-1) + hs(i-1,j-1)
				Hdata1(i,j) = hs(i-1, j-1)
				Hdata2(i,j) = h_max(i-1, j-1)
				Vdata2(i,j) = v_max(i-1, j-1)
				Udata1(i,j) = u(i-1,j-1)
				Vdata1(i,j) = v(i-1,j-1)
				Rdata(i,j)  = rain(i-1,j-1)/0.001*3600.
                ibc0(i,j)   = ibcg0(i-1, j-1) 
			ENDDO
		ENDDO

		CALL CG_IRIC_WRITE_SOL_GRIDCOORD2D_F(xx,yy,IER)
		CALL CG_IRIC_WRITE_SOL_REAL_F("Velocity(ms-1)X",UData1,IER)
		CALL CG_IRIC_WRITE_SOL_REAL_F("Velocity(ms-1)Y",VData1,IER)
		CALL CG_IRIC_WRITE_SOL_REAL_F("Depth(Max)",HData2,IER)
		CALL CG_IRIC_WRITE_SOL_REAL_F("Depth",HData1,IER)
		CALL CG_IRIC_WRITE_SOL_REAL_F("Elevation",Zbdata1,IER)
		CALL CG_IRIC_WRITE_SOL_REAL_F("WaterSurfaceElevation",WSE,IER)
		if(j_rain>=2) CALL CG_IRIC_WRITE_SOL_REAL_F("Rain(mmh-1)",Rdata,IER)
		CALL CG_IRIC_WRITE_SOL_REAL_F("Velocity (magnitude Max)",VData2,IER)
        !CALL CG_IRIC_WRITE_SOL_integer_F("IBC",ibc0,IER)
		
		DEALLOCATE(xx, STAT=ier)
		DEALLOCATE(yy, STAT=ier)
		DEALLOCATE(Udata1, STAT=ier)
		DEALLOCATE(Vdata1, STAT=ier)
		DEALLOCATE(Hdata1, STAT=ier)
		DEALLOCATE(Zbdata1, STAT=ier)

	END SUBROUTINE

	END MODULE
