//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                                             
//   trilinosModelEvaluator.cpp - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
//                                                              
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   Copyright (C) 2005-2013
//   Glimmer-CISM contributors - see AUTHORS file for list of contributors
//
//   This file is part of Glimmer-CISM.
//
//   Glimmer-CISM is free software: you can redistribute it and/or modify it
//   under the terms of the Lesser GNU General Public License as published
//   by the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   Glimmer-CISM is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   Lesser GNU General Public License for more details.
//
//   You should have received a copy of the Lesser GNU General Public License
//   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "trilinosModelEvaluator.hpp"
#include "Teuchos_StandardCatchMacros.hpp"


extern "C" {
  void calc_F(double* x, double* f, int N, void* bb, int ispert);
  void apply_precond_nox(double* x, double* y, int n, void* bb);
  void reset_effstrmin(const double* esm);
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

trilinosModelEvaluator::trilinosModelEvaluator (
                            int N_, double* statevector,
                            const Epetra_Comm& comm_, void* blackbox_res_)
			     : N(N_), comm(comm_), blackbox_res(blackbox_res_)
{
  bool succeeded=true;
  try {
    xMap = Teuchos::rcp(new Epetra_Map(-1, N, 0, comm));
    xVec = Teuchos::rcp(new Epetra_Vector(Copy, *xMap, statevector));

    precOp = Teuchos::rcp(new trilinosPreconditioner(N, xVec, xMap, blackbox_res));

    pMap = Teuchos::rcp(new Epetra_LocalMap(1, 0, comm));
    pVec = Teuchos::rcp(new Epetra_Vector(*pMap));
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, succeeded);
  if (!succeeded) exit(1);
}

/*******************************************************************************/
// Return solution vector map
Teuchos::RCP<const Epetra_Map> trilinosModelEvaluator::get_x_map() const{
  return xMap;
}

// Return residual vector map
Teuchos::RCP<const Epetra_Map> trilinosModelEvaluator::get_f_map() const{
  return xMap;
}

// Return initial solution and x_dot init
Teuchos::RCP<const Epetra_Vector> trilinosModelEvaluator::get_x_init() const{
  return xVec;
}

Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner>
trilinosModelEvaluator::create_WPrec() const
{
  // bool is answer to: "Prec is already inverted?"
  return Teuchos::rcp(new EpetraExt::ModelEvaluator::Preconditioner(precOp,true));
}

Teuchos::RCP<const Epetra_Map> trilinosModelEvaluator::get_p_map(int l) const{
  return pMap;
}
Teuchos::RCP<const Epetra_Vector> trilinosModelEvaluator::get_p_init(int l) const{
  return pVec;
}
Teuchos::RCP<const  Teuchos::Array<std::string> >  trilinosModelEvaluator::get_p_names(int l) const{
    RCP<Teuchos::Array<std::string> > p_names =
      rcp(new Teuchos::Array<std::string>(1) );
    (*p_names)[0] = "Effstrmin Factor";

  return p_names;
}

/*******************************************************************************/
// Create InArgs
EpetraExt::ModelEvaluator::InArgs trilinosModelEvaluator::createInArgs() const{
  InArgsSetup inArgs;

  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(IN_ARG_x,true);
  inArgs.set_Np(1);
  return inArgs;
}

/*******************************************************************************/
// Create OutArgs
EpetraExt::ModelEvaluator::OutArgs trilinosModelEvaluator::createOutArgs() const{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(1, 0);
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_WPrec, true);

  return outArgs;
}

/*******************************************************************************/
// Evaluate model on InArgs
void trilinosModelEvaluator::evalModel(const InArgs& inArgs, const OutArgs& outArgs) const{

  // Get the solution vector x from inArgs and residual vector from outArgs
  RCP<const Epetra_Vector> x = inArgs.get_x();
  EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector> f = outArgs.get_f();
  
  if (x == Teuchos::null) throw "trilinosModelEvaluator::evalModel: x was NOT specified!";

  // Check if a "Effminstr Factor" parameter is being set by LOCA
  Teuchos::RCP<const Epetra_Vector> p_in = inArgs.get_p(0);
  if (p_in.get()) reset_effstrmin(&(*p_in)[0]);

  // Save the current solution, which makes it initial guess for next nonlienar solve
  *xVec = *x;

  if (f != Teuchos::null) {
    // Check if this is a perturbed eval. Glimmer only saves off matrices for unperturbed case.
    int ispert =0;
    if  (f.getType() == EpetraExt::ModelEvaluator::EVAL_TYPE_APPROX_DERIV) ispert=1;

    f->PutScalar(0.0);
    calc_F(x->Values(), f->Values(), N, blackbox_res, ispert);
  }

  RCP<Epetra_Operator> WPrec = outArgs.get_WPrec();
  if (WPrec != Teuchos::null) {
     //cout << "evalModel called for WPrec -- doing nothing " <<  endl;
  }
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
trilinosPreconditioner::trilinosPreconditioner (
       int N_, RCP<Epetra_Vector> xVec_, RCP<Epetra_Map> xMap_, void* blackbox_res_)
       : N(N_), xVec(xVec_), xMap(xMap_), blackbox_res(blackbox_res_)
{
}

int trilinosPreconditioner::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  bool succeeded=true;
  try {
    apply_precond_nox(Y(0)->Values(), X(0)->Values(), N, blackbox_res);
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, succeeded);
  if (!succeeded) exit(1);

  return 0;
}


