module mpas_atm_dimensions

#ifdef CONST_INNER_DIMS

#ifndef CONST_NVERTLEVELS
#error "Defining CONST_INNER_DIMS requires CONST_NVERTLEVELS, CONST_MAXEDGES, CONST_MAXEDGES2, and CONST_NUM_SCALARS to be defined as well."
#endif
#ifndef CONST_MAXEDGES
#error "Defining CONST_INNER_DIMS requires CONST_NVERTLEVELS, CONST_MAXEDGES, CONST_MAXEDGES2, and CONST_NUM_SCALARS to be defined as well."
#endif
#ifndef CONST_MAXEDGES2
#error "Defining CONST_INNER_DIMS requires CONST_NVERTLEVELS, CONST_MAXEDGES, CONST_MAXEDGES2, and CONST_NUM_SCALARS to be defined as well."
#endif
#ifndef CONST_NUM_SCALARS
#error "Defining CONST_INNER_DIMS requires CONST_NVERTLEVELS, CONST_MAXEDGES, CONST_MAXEDGES2, and CONST_NUM_SCALARS to be defined as well."
#endif

    integer, parameter :: nVertLevels = CONST_NVERTLEVELS
    integer, parameter :: maxEdges = CONST_MAXEDGES
    integer, parameter :: maxEdges2 = CONST_MAXEDGES2
    integer, parameter :: num_scalars = CONST_NUM_SCALARS
#else
    integer :: nVertLevels
    integer :: maxEdges
    integer :: maxEdges2
    integer :: num_scalars
#endif


    contains


    subroutine mpas_atm_set_dims(nVertLevels_val, maxEdges_val, maxEdges2_val, num_scalars_val)

        implicit none

        integer, intent(in) :: nVertLevels_val
        integer, intent(in) :: maxEdges_val
        integer, intent(in) :: maxEdges2_val
        integer, intent(in) :: num_scalars_val

#ifdef CONST_INNER_DIMS
        if (nVertLevels /= nVertLevels_val) then
            write(0,*) '********************************************************************************'
            write(0,*) 'ERROR: Dimension nVertLevels read from input file does not match the value used'
            write(0,*) '       when compiling MPAS-Atmosphere (CONST_NVERTLEVELS=', CONST_NVERTLEVELS ,')'
            call mpas_dmpar_global_abort('********************************************************************************')
        end if

        if (maxEdges /= maxEdges_val) then
            write(0,*) '********************************************************************************'
            write(0,*) 'ERROR: Dimension maxEdges read from input file does not match the value used'
            write(0,*) '       when compiling MPAS-Atmosphere (CONST_MAXEDGES=', CONST_MAXEDGES ,')'
            call mpas_dmpar_global_abort('********************************************************************************')
        end if

        if (maxEdges2 /= maxEdges2_val) then
            write(0,*) '********************************************************************************'
            write(0,*) 'ERROR: Dimension maxEdges2 read from input file does not match the value used'
            write(0,*) '       when compiling MPAS-Atmosphere (CONST_MAXEDGES2=', CONST_MAXEDGES2 ,')'
            call mpas_dmpar_global_abort('********************************************************************************')
        end if

        if (num_scalars /= num_scalars_val) then
            write(0,*) '********************************************************************************'
            write(0,*) 'ERROR: Number of scalars read from input file does not match the value used'
            write(0,*) '       when compiling MPAS-Atmosphere (CONST_NUM_SCALARS=', CONST_NUM_SCALARS ,')'
            call mpas_dmpar_global_abort('********************************************************************************')
        end if
#else
        nVertLevels = nVertLevels_val
        maxEdges = maxEdges_val
        maxEdges2 = maxEdges2_val
        num_scalars = num_scalars_val
#endif

    end subroutine mpas_atm_set_dims

end module mpas_atm_dimensions
