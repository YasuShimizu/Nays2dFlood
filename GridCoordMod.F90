Module GridCoord

	IMPLICIT NONE

  CONTAINS

	SUBROUTINE CGNS_READ_GRIDCOORD(FID, nx, ny, xx, yy, IER)
		IMPLICIT NONE

		INTEGER, INTENT(IN) :: FID
		INTEGER, INTENT(OUT) :: nx, ny
		double precision, DIMENSION(:,:), allocatable :: xx, yy
		INTEGER, INTENT(OUT) :: IER
		INTEGER :: status, i, j, count, countji, ierror
	
		CALL CG_IRIC_GOTOGRIDCOORD2D_F(nx, ny, IER)

		ALLOCATE(xx(1:NX, 1:NY), STAT = status)
		ALLOCATE(yy(1:NX, 1:NY), STAT = status)
	
		CALL CG_IRIC_GETGRIDCOORD2D_F(xx,yy,IER)

	END SUBROUTINE CGNS_READ_GRIDCOORD

END MODULE GridCoord