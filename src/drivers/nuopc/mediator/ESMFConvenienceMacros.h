#if 0
// ----------- ERROR handling macros ------------------------------------------
#endif

#define ESMF_ERR_ABORT(rc) if (ESMF_LogFoundError(rc, msg="Aborting NEMS", line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

#define ESMF_ERR_RETURN(rc,rcOut) if (ESMF_LogFoundError(rc, msg="Breaking out of subroutine", line=__LINE__, file=__FILE__, rcToReturn=rcOut)) return
