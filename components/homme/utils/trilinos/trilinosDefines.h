// Select which preconditioner 
// NOTE: When adding a new preconditioner, both trilinosNoxSolver and 
//       trilinosModelEvaluator need to be updated with the new #define
// ==============================================================================

// Identity preconditioner 
#define IDENT_PC

// Block preconditioner using spectral element operators
//#define BLOCK_PC_OP

// Block preconditioner using local dense Jacobian block
//#define BLOCK_PC_DENSE

// Block preconditioner using a global sparse Jacobian matrix and Belos
//#define BLOCK_PC_SPARSE

// Block preconditioner using a global sparse Jacobian matrix and ML
//#define BLOCK_PC_SPARSE_ML

// Block preconditioner using a global sparse Jacobian matrix and Belos
//#define GMRES_PC_OP

// Block preconditioner using a global sparse Jacobian matrix and Belos
//#define GMRES_PC_SPARSE

// Multilevel preconditioner using ML package
//#define ML_PC

// Compare analytic Jacobian calculation against finite difference
//#define CHECK_JAC


