MODULE GridCond

	USE GridCoord

  CONTAINS

	SUBROUTINE CGNS_Read_GridCondition(FID,nx, ny, z1, z2, iobst, vege, ch, fm, ob, q_cell, sh, ibc0, IER)

		IMPLICIT NONE

		INTEGER, INTENT(IN) :: FID
		INTEGER, INTENT(INOUT) :: nx, ny
		INTEGER, INTENT(OUT) :: IER
        integer, dimension(:,:), allocatable :: ibc0
		double precision, DIMENSION(:,:), allocatable :: z1,z2, ch, sh
		INTEGER(4), DIMENSION(:,:), allocatable :: iobst, vege, fm, ob, q_cell
		INTEGER :: status
	
		ALLOCATE(z1(1:NX,1:NY), STAT=status)
		ALLOCATE(z2(1:NX,1:NY), STAT=status)
		ALLOCATE( iobst(1:NX-1,1:NY-1), STAT=status)
        ALLOCATE( ibc0(1:NX-1,1:NY-1), STAT=status)
		ALLOCATE(  vege(1:NX-1,1:NY-1), STAT=status)
		ALLOCATE(    ch(1:NX-1,1:NY-1), STAT=status)
		ALLOCATE(    fm(1:NX-1,1:NY-1), STAT=status)
		ALLOCATE(    ob(1:NX-1,1:NY-1), STAT=status)
		ALLOCATE(q_cell(1:NX-1,1:NY-1), STAT=status)
		ALLOCATE(    sh(1:NX-1,1:NY-1), STAT=status)
   
		CALL CG_IRIC_READ_GRID_REAL_NODE_F('Elevation', z1, IER);
		CALL CG_IRIC_READ_GRID_REAL_NODE_F('water_depth', z2, IER);
		
	!	z2 = z1
    
		CALL CG_IRIC_READ_GRID_INTEGER_CELL_F('Obstacle',iobst,IER)
		CALL CG_IRIC_READ_GRID_INTEGER_CELL_F('Fixed_movable',fm,IER)
		CALL CG_IRIC_READ_GRID_INTEGER_CELL_F('Bomb',ob,IER)
		CALL CG_IRIC_READ_GRID_REAL_CELL_F('channelcondition',ch,IER)
		CALL CG_IRIC_READ_GRID_INTEGER_CELL_F('q_cell',q_cell,IER)
!		CALL CG_IRIC_READ_GRID_INTEGER_CELL_F('vegetation',vege,IER)
		CALL CG_IRIC_READ_GRID_REAL_CELL_F('Share',sh,IER)
        
        CALL CG_IRIC_READ_GRID_INTEGER_CELL_F('ibc',ibc0,IER)

	END SUBROUTINE CGNS_Read_GridCondition

	SUBROUTINE CGNS_Read_BoundaryCondition(j_in,j_inlen,indices, indexmax, j_size, jsizemax,f_param,f_value,slopevalue,flowname)

!docon add start#####################################################################################
        use box_gate_pump
!docon add end#######################################################################################

		IMPLICIT NONE

		integer :: j_in, index, indexmax, i, ier, jsizemax
!docon add start#####################################################################################
        integer :: j
	    integer :: bc_mxindx2,gt_mxindx2,p_mxindx2
!docon add end#######################################################################################
		integer,dimension(:),allocatable :: j_inlen, j_size
		integer, dimension(:,:,:), allocatable:: indices
		double precision, dimension(:,:), allocatable:: f_param
		double precision, dimension(:,:), allocatable:: f_value
		double precision,dimension(:),allocatable :: slopevalue
        character(250),dimension(:),allocatable :: flowname
	
		call cg_iric_read_bc_count_f('inflow', j_in)

		allocate(j_inlen(j_in))
		allocate(j_size(j_in))
		allocate(slopevalue(j_in))
      allocate(flowname(j_in))

		indexmax = 0
		jsizemax = 0

		do i = 1, j_in
			call cg_iric_read_bc_indicessize_f('inflow', i, j_inlen(i), ier)

			if(indexmax.lt. j_inlen(i)) indexmax= j_inlen(i)

            call cg_iric_read_bc_string_f('inflow',i,'_caption',flowname(i), ier)
			call cg_iric_read_bc_functionalsize_f('inflow', i, 'qt1', j_size(i), ier);

			if (jsizemax < j_size(i)) then
				jsizemax = j_size(i)
			end if
			
		end do

		allocate(indices(j_in, 2, indexmax))
		allocate(f_param(j_in, jsizemax), f_value(j_in, jsizemax))

		do i = 1, j_in
			call cg_iric_read_bc_indices_f('inflow', i, indices(i:i,:,:), ier)
			call cg_iric_read_bc_real_f('inflow', i, 'slope1', slopevalue(i:i), ier)
			call cg_iric_read_bc_functional_f('inflow',i, 'qt1', f_param(i:i,:), f_value(i:i,:), ier)
		end do

!docon add start#####################################################################################
        !以下、ボックスカルバート・樋門・ポンプの境界条件追加
		call cg_iric_read_bc_count_f('BoxCulvert', bc_num)
		call cg_iric_read_bc_count_f('Gate', gt_num)
		call cg_iric_read_bc_count_f('Pump', p_num)

        !動的配列宣言-----------------------------------------------------------------------------
        allocate(bc_couple_num(bc_num),bc_inout(bc_num),bc_c2_yesno(bc_num),bc_width(bc_num),bc_height(bc_num))
        allocate(bc_base(bc_num),bc_ren(bc_num),bc_c1(bc_num),bc_c2(bc_num),bc_c3(bc_num),bc_mxindx(bc_num))
        allocate(gt_couple_num(gt_num),gt_inout(gt_num),gt_c2_yesno(gt_num),gt_method(gt_num))
        allocate(gt_width(gt_num),gt_height(gt_num),gt_base(gt_num),gt_ren(gt_num),gt_c1(gt_num))
        allocate(gt_c2(gt_num),gt_c3(gt_num),gt_outwl(gt_num))
        allocate(gt_size(gt_num),gt_size2(gt_num),gt_mxindx(gt_num))
        allocate(p_couple_num(p_num),p_method(p_num),p_inout(p_num),p_qmax(p_num))
        allocate(p_startope(p_num),p_stopope(p_num),p_qt(p_num),p_size(p_num),p_mxindx(p_num))
        allocate(BCname(bc_num),Gtname(gt_num),Ppname(p_num))   !140602
        allocate(bc_qswap(bc_num),gt_qswap(gt_num),p_qswap(p_num))  !140602

        !ボックスカルバート関連データ入力---------------------------------------------------------
        bc_mxindx2=0
        do i=1,bc_num   !データ数分繰り返し
			call cg_iric_read_bc_indicessize_f('BoxCulvert', i, bc_mxindx(i), ier)  !セル数取得
            call cg_iric_read_bc_string_f('BoxCulvert',i,'_caption',BCname(i), ier) !140602

            !複数メッシュの選択は不可
            if(bc_mxindx(i)>1) then
                write(*,*) 'BoxCulvert-cell is not allowed multiple selection'
                stop
            end if

			if(bc_mxindx2<bc_mxindx(i)) bc_mxindx2= bc_mxindx(i)
            call cg_iric_read_bc_INTEGER_F('BoxCulvert', i, 'bc_couple_num', bc_couple_num(i), ier)
            call cg_iric_read_bc_INTEGER_F('BoxCulvert', i, 'bc_inout', bc_inout(i), ier)   !入口=0,出口=1
            if(bc_inout(i)==0) then         !入口側の場合
                call cg_iric_read_bc_INTEGER_F('BoxCulvert', i, 'bc_c2_yesno', bc_c2_yesno(i), ier)
                call cg_iric_read_bc_REAL_F('BoxCulvert', i, 'bc_width', bc_width(i), ier)
                call cg_iric_read_bc_REAL_F('BoxCulvert', i, 'bc_height', bc_height(i), ier)
                call cg_iric_read_bc_REAL_F('BoxCulvert', i, 'bc_base', bc_base(i), ier)
                call cg_iric_read_bc_REAL_F('BoxCulvert', i, 'bc_ren', bc_ren(i), ier)
                call cg_iric_read_bc_REAL_F('BoxCulvert', i, 'bc_c1', bc_c1(i), ier)
                call cg_iric_read_bc_REAL_F('BoxCulvert', i, 'bc_c3', bc_c3(i), ier)
                if(bc_c2_yesno(i)==1) then  !中間流出考慮の場合
                    call cg_iric_read_bc_REAL_F('BoxCulvert', i, 'bc_c2', bc_c2(i), ier)
                end if
            end if
        end do

        allocate(bc_indx(bc_num,2,bc_mxindx2))
        do i=1,bc_num   !データ数分繰り返し
			call cg_iric_read_bc_indices_f('BoxCulvert', i, bc_indx(i,:,:), ier)
			if(i==1) write(*,*) '------------------------BoxCulvert------------------------'
			if(i==1) write(*,'(a,i6)') 'BoxCulvert data= ',bc_num
			do j=1,bc_mxindx(i)
                write(*,'(4(a,i6))') 'Num= ',i,',  Index= ',j,',  I= ',bc_indx(i,1,j),',  J= ',bc_indx(i,2,j)
            end do
        end do

        !樋門関連データ入力-----------------------------------------------------------------------
        gt_tmax=0
        gt_tmax2=0
        gt_mxindx2=0
        do i=1,gt_num   !データ数分繰り返し
			call cg_iric_read_bc_indicessize_f('Gate', i, gt_mxindx(i), ier)    !セル数取得
            call cg_iric_read_bc_string_f('Gate',i,'_caption',Gtname(i), ier)   !140602

            !複数メッシュの選択は不可
            if(gt_mxindx(i)>1) then
                write(*,*) 'Gate-cell is not allowed multiple selection'
                stop
            end if

			if(gt_mxindx2<gt_mxindx(i)) gt_mxindx2= gt_mxindx(i)
	        call cg_iric_read_bc_INTEGER_F('Gate', i, 'gt_couple_num', gt_couple_num(i), ier)
	        call cg_iric_read_bc_INTEGER_F('Gate', i, 'gt_inout', gt_inout(i), ier)
            if(gt_inout(i)==0) then         !入口側の場合
    	        call cg_iric_read_bc_INTEGER_F('Gate', i, 'gt_method', gt_method(i), ier)
	            call cg_iric_read_bc_INTEGER_F('Gate', i, 'gt_c2_yesno', gt_c2_yesno(i), ier)
	            call cg_iric_read_bc_REAL_F('Gate', i, 'gt_width', gt_width(i), ier)
	            call cg_iric_read_bc_REAL_F('Gate', i, 'gt_height', gt_height(i), ier)
	            call cg_iric_read_bc_REAL_F('Gate', i, 'gt_base', gt_base(i), ier)
	            call cg_iric_read_bc_REAL_F('Gate', i, 'gt_ren', gt_ren(i), ier)
	            call cg_iric_read_bc_REAL_F('Gate', i, 'gt_c1', gt_c1(i), ier)
	            call cg_iric_read_bc_REAL_F('Gate', i, 'gt_c3', gt_c3(i), ier)
                if(gt_c2_yesno(i)==1) then  !中間流出考慮の場合
                    call cg_iric_read_bc_REAL_F('Gate', i, 'gt_c2', gt_c2(i), ier)
                end if
		        if(gt_method(i)==1) then
		            call cg_iric_read_bc_functionalsize_f('Gate', i, 'gt_openclose', gt_size(i), ier);
		            if (gt_tmax < gt_size(i)) gt_tmax = gt_size(i)
		        end if
		    end if
	        if(gt_inout(i)==1) then         !出口側の場合
	            call cg_iric_read_bc_INTEGER_F('Gate', i, 'gt_outwl', gt_outwl(i), ier)
	            if(gt_outwl(i)==1) then     !外水位が外読みファイルで与えられているとき
	                call cg_iric_read_bc_functionalsize_f('Gate', i, 'wl1', gt_size2(i), ier)
		            if (gt_tmax2 < gt_size2(i)) gt_tmax2 = gt_size2(i)
                end if
            end if
        end do
        allocate(gt_time(gt_num,gt_tmax),gt_time2(gt_num,gt_tmax2))
        allocate(gt_state(gt_num,gt_tmax),gt_wl(gt_num,gt_tmax2))
        allocate(gt_indx(gt_num,2,gt_mxindx2))

        do i=1,gt_num   !データ数分繰り返し
			call cg_iric_read_bc_indices_f('Gate', i, gt_indx(i,:,:), ier)
			if(i==1) write(*,*) '------------------------Gate------------------------'
			if(i==1) write(*,'(a,i6)') 'Gate data= ',gt_num
			do j=1,gt_mxindx(i)
                write(*,'(4(a,i6))') 'Num= ',i,',  Index= ',j,',  I= ',gt_indx(i,1,j),',  J= ',gt_indx(i,2,j)
            end do
            if(gt_inout(i)==0) then     !入口側の場合
    			if(gt_method(i)==1) &   !樋門の開閉操作が外読みファイルで与えられているとき
    			call cg_iric_read_bc_functional_f('Gate', i, 'gt_openclose', gt_time(i,:), gt_state(i,:), ier)
            end if
	        if(gt_inout(i)==1) then     !出口側の場合
	            call cg_iric_read_bc_INTEGER_F('Gate', i, 'gt_outwl', gt_outwl(i), ier)
	            if(gt_outwl(i)==1) &    !外水位が外読みファイルで与えられているとき
	            call cg_iric_read_bc_functional_f('Gate', i, 'wl1', gt_time2(i,:), gt_wl(i,:), ier)
            end if
        end do

        !ポンプ関連データ入力---------------------------------------------------------------------
        p_tmax=0
        p_mxindx2=0
        do i=1,p_num    !データ数分繰り返し
			call cg_iric_read_bc_indicessize_f('Pump', i, p_mxindx(i), ier) !セル数取得
            call cg_iric_read_bc_string_f('Pump',i,'_caption',Ppname(i), ier)   !140602
            !複数メッシュの選択は不可
       !     if(p_mxindx(i)>1) then
       !         write(*,*) 'Pump-cell is not allowed multiple selection'
       !         stop
       !     end if

			if(p_mxindx2<p_mxindx(i)) p_mxindx2= p_mxindx(i)
	        call cg_iric_read_bc_INTEGER_F('Pump', i, 'p_couple_num', p_couple_num(i), ier)
	        call cg_iric_read_bc_INTEGER_F('Pump', i, 'p_inout', p_inout(i), ier)
            if(p_inout(i)==0) then          !入口側の場合
	            call cg_iric_read_bc_INTEGER_F('Pump', i, 'p_method', p_method(i), ier)
	            if(p_method(i)==1) then     !ポンプ排水量が外読みファイルで与えられているとき
		            call cg_iric_read_bc_functionalsize_f('Pump', i, 'qt2', p_size(i), ier);
                end if
	            call cg_iric_read_bc_REAL_F('Pump', i, 'p_qmax', p_qmax(i), ier)
	            call cg_iric_read_bc_REAL_F('Pump', i, 'p_startope', p_startope(i), ier)
	            call cg_iric_read_bc_REAL_F('Pump', i, 'p_stopope', p_stopope(i), ier)
            end if
            if(p_inout(i)==1) then          !出口側の場合
	            call cg_iric_read_bc_functionalsize_f('Pump', i, 'qt2', p_size(i), ier);
            end if
		    if (p_tmax < p_size(i)) p_tmax = p_size(i)
        end do
        allocate(p_time(p_num,p_tmax),p_qout(p_num,p_tmax))
        allocate(p_indx(p_num,2,p_mxindx2))

        do i=1,p_num   !データ数分繰り返し
			call cg_iric_read_bc_indices_f('Pump', i, p_indx(i,:,:), ier)
			if(i==1) write(*,*) '------------------------Pump------------------------'
			if(i==1) write(*,'(a,i6)') 'Pump data= ',p_num
			do j=1,p_mxindx(i)
                write(*,'(4(a,i6))') 'Num= ',i,',  Index= ',j,',  I= ',p_indx(i,1,j),',  J= ',p_indx(i,2,j)
            end do
            if(p_inout(i)==0) then          !入口側の場合
	            if(p_method(i)==1) then     !ポンプ排水量が外読みファイルで与えられているとき
	                call cg_iric_read_bc_functional_f('Pump', i, 'qt2', p_time(i,:), p_qout(i,:), ier)
                end if
            end if
            if(p_inout(i)==1) then          !出口側の場合
                call cg_iric_read_bc_functional_f('Pump', i, 'qt2', p_time(i,:), p_qout(i,:), ier)
            end if
        end do
!docon add end#######################################################################################

	END SUBROUTINE CGNS_Read_BoundaryCondition

END MODULE GridCond
