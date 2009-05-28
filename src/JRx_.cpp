//  JRx_.cpp
// 
//  Copyright (C) 2008 Laboratoire Statistique & Genome
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or (at
//  your option) any later version.
// 
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 

//! \file    JRx_.cpp
//! \author  Julien Chiquet, Gilles Grasseau
//! \brief   See the 'SIMoNe' R package
//! \version 0.1
//! \date    June 2008

#include <math.h>

#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

#include <cmath>

#include "matrix.h"

//! \cond
extern "C" {
  void JRx_C( 
	                       // Input
	     int &n,           //  Node number
	     int &Q,           //  Class number
	     double *tau_ ,    //  Tau(n,Q) proba to belong to class q
	     double *x_,       //  X(n,n) adjacency matrix
	     double *mu_,      //  mu(Q,Q) mean matrix
	     double *lambda_,  //  lambda(Q,Q) dispersion matrix
	     double *alpha_,   //   alpha(Q) class ratio
	                       // Ouput
	     double &J         //  J relurn value
	    ); 
}

Matrix lnFLaplace( Matrix &X, const double mu, const double lambda) {

  Matrix Tmp(X.getRow(), X.getCol());

  double log_eps  = log(DBL_MIN) ;
  double log_huge = log(DBL_MAX);   
  double eps      = DBL_MIN;

  if( lambda > eps ) {
    double log_2lambda = log(2*lambda);
    Tmp = - abs( X - mu ) / lambda - log_2lambda;  
  } else {
    Tmp = Tmp.Compare("==", mu, 0.0, log_eps );  
  }
  Tmp = max( Tmp, log_eps  ); 
  Tmp = min( Tmp, log_huge ); 
  return (Tmp);
}


void JRx_C( int &n, int &Q, 
	    double *tau_ , double *x_, 
	    double *mu_, double *lambda_, double *alpha_,
	    double &J ) {

  // Argument copies
  Matrix Tau   (n,Q, tau_);
  Matrix X     (n,n, x_);
  Matrix mu    (Q,Q, mu_);
  Matrix lambda(Q,Q, lambda_);
  Matrix alpha (Q,1, alpha_);

  // Define the smallest double
  double eps     = DBL_MIN; 
  double log_eps = log( eps );
  
  Matrix logTau(n,Q); 
  Matrix logalpha(Q,1);
  
  try {
    logTau   = max( log( Tau ),   log_eps ); 
    logalpha = max( log( alpha ), log_eps );
    
    J =  - sum( Tau | logTau ) +  sum ( Tau * logalpha ); 
    
    Matrix lnf(Q,Q);
    Matrix tau(Q,Q);
    
    for (int q=0; q < Q; q++) {
      for (int l=0; l < Q; l++) {
	if ( lambda(q,l) > eps ) {
	  lnf = lnFLaplace(X, mu(q,l), lambda(q,l)) ;
	  tau = Tau.getColVector(q) * Tau.getColVector(l).t();
	  J = J + sum( tau.getUpperTriang() |  lnf.getUpperTriang() );
	}
	else {
	  // Not computed here
	}
      }
    }
  } catch ( std::exception &e) {
    LOGEXCEPT( e );
  }
  LOGFLUSH;
}

//! \endcond
