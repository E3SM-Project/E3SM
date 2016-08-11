#ifdef FORTRAN_SAME
#define FC_FUNC(name,NAME) name
#elif FORTRAN_UNDERSCORE_
#define FC_FUNC(name,NAME) name ##_
#elif FORTRAN_DOUBLE_UNDERSCORE_
#define FC_FUNC(name,NAME)  name ##__
#endif
