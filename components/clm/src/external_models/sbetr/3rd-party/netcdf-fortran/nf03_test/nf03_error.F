#if 0
    Copyright 1996, UCAR/Unidata
    See netcdf/COPYRIGHT file for copying and redistribution conditions.
    $Id: nf_error.F,v 1.4 1997/06/03 22:26:09 steve Exp $
#endif


C
C Use for logging error messages
C
        subroutine error(msg)
        use tests
        implicit        none
        character*(*)   msg

        nfails = nfails + 1
        if (nfails .le. max_nmpt) print *, msg
        end


C
C Use for logging error conditions
C
        subroutine errori(msg, i)
        use tests
        implicit        none
        character*(*)   msg
        integer         i

        nfails = nfails + 1
        if (nfails .le. max_nmpt) print *, msg, i
        end


C
C Use for logging error conditions
C
        subroutine errord(msg, d)
        use tests
        implicit        none
        character*(*)   msg
        doubleprecision d

        nfails = nfails + 1
        if (nfails .le. max_nmpt) print *, msg, d
        end


C
C Use for logging error conditions
C
        subroutine errorc(msg, string)
        use tests
        implicit        none
        character*(*)   msg
        character*(*)   string

        nfails = nfails + 1
        if (nfails .le. max_nmpt) print *, msg, 
     +      string(1:len_trim(string))
        end


C
C Use for logging error conditions
C
        subroutine errore(msg, err)
        use tests
        implicit        none
        character*(*)   msg
        integer         err

        nfails = nfails + 1
        call errorc(msg, nf_strerror(err))
        end
