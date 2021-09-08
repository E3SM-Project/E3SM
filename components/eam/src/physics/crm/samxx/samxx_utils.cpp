
#include "samxx_utils.h"

// scream session initialize/finalize
extern "C" {
 void scream_session_init()
 {
   scream::initialize_scream_session();
 }

 void scream_session_finalize()
 {
  scream::finalize_scream_session();
 }
}

