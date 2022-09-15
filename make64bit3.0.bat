ifort 	nays2d_flood\common_hh.f90 ^
	nays2d_flood\GridCoordMod.F90 ^
	nays2d_flood\GridCondMod.F90 ^
	nays2d_flood\interpolation.F90 ^
	nays2d_flood\transform.F90 ^
	nays2d_flood\WriteCGNS.F90 ^
	nays2d_flood\output_hotstart.f90 ^
	nays2d_flood\read_hotstart_conditions.f90 ^
	nays2d_flood\flow_calculation.f90 ^
	nays2d_flood\libs\cgnsdll.lib ^
	nays2d_flood\libs\iriclib.lib ^
	nays2d_flood\nays2d_flood.f90 /MD /Qopenmp /nostandard-realloc-lhs -o nays2dflood_v5_64bit\nays2d_flood.exe

del *.obj 
del *.mod 

