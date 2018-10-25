/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_DEBUG_HPP
#define HOMMEXX_DEBUG_HPP

#ifndef NDEBUG

// Note: This does not seem to do much more than what the standard
//       MPI_ERRORS_ARE_FATAL handling does. However, by calling
//       MPI_Abort, one can easily stop the execution in gdb by
//       setting a breakpoint at MPI_Abort:
//
//         (gdb) b MPI_Abort
//
//       Furthermore, we display additional details about the location
//       of the error (so that one knows where the problem is, even
//       when running outside of a debugger or without a breakpoint).

#define HOMMEXX_MPI_CHECK_ERROR(X,comm)                       \
  {                                                           \
    int err_code = X;                                         \
    if (err_code!=MPI_SUCCESS) {                              \
      char err_str[MPI_MAX_ERROR_STRING];                     \
      int resultlen;                                          \
      MPI_Error_string(err_code,err_str,&resultlen);          \
      printf("Hommexx mpi error: %s\n",err_str);              \
      printf("   at line %d of file %s\n",__LINE__,__FILE__); \
      MPI_Abort(comm,err_code);                               \
    }                                                         \
  }

#else
#define HOMMEXX_MPI_CHECK_ERROR(X,comm) X
#endif

#endif // HOMMEXX_DEBUG_HPP
