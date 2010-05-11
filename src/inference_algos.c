// Compile with :  R CMD SHLIB inference_algos.c
#include "inference_utils.h"

/** ______________________________________________________________ **/
/** ______________________________________________________________ **/
/**                                                                **/
/** 			THE MAIN FUNCTION GLASSO                   **/
/** ______________________________________________________________ **/
/** ______________________________________________________________ **/
/**                                                                **/
/** RÉSOLUTION EXACTE DE LA VRAISEMBLANCE PÉNALISÉE                **/
/** Problème de Banerjee et al, Friedman et al                     **/
/** ______________________________________________________________ **/
/**                                                                **/
/**
 ** in :
 **	p   = colonne en cours
 **	S   = matrice S
 **	Rho = matrice de pénalités Rho
 **	eps 	= seuil de convergence pour la matrice W (n°.int.) si approx=0
 **	maxIt   = nombre maximum d'iterations
 **
 ** out :
 **	W       = matrice de covaraince estimée
 ** 	finalIt = nombre d'itérations requises
 **/
void GLasso(int *p,
	    double S[*p][*p],
	    double Rho[*p][*p],
	    double *eps,
	    int *maxIt,
	    double W[*p][*p],
	    double Beta[*p][*p],
	    int *finalIt) {
  
  /** INITIALISATION**/
  int i,k,active, Sp[*p-1];
  double W11[*p-1][*p-1], wk[*p];
  double rho12[*p-1],s12[*p-1],beta[*p-1];
  double threshold1 = *eps * NormMat(*p,W) / *p;
  double threshold2 = *eps / *p;
  double err = 2* threshold1;
  
  while( (*finalIt < *maxIt) && (err > threshold1) ) {
    err = 0.0;
    for(k=0; k<*p; k++) { // let's work on the kth column of W

      /** Setup prépare les vecteurs pour la valeur courante de k **/
      /** W11, rho12, s12, beta et wk                             **/
      Setup_GLasso(k, *p, S, Rho, W, Beta, s12, rho12, wk, beta, W11, Sp, &active);
      
      /** Lasso via pathwise coordinate **/
      Pathwise_Coordinate(*p-1, W11, s12, rho12, threshold2, beta, Sp, &active);

      /** Reconstruction of the full matrix W **/
      Cleanup(k, *p, W, W11, beta, Beta, Sp, &active);

      /** Checking for convergence between old W[,k] and new W[,k]... **/
      for(i=0; i<*p; i++) {
	err += fabs(W[k][i]-wk[i]);
      }
      
      R_CheckUserInterrupt();
    }
    *finalIt = *finalIt + 1;
  }
} // end Glasso

/** ______________________________________________________________ **/
/** ______________________________________________________________ **/
/**                                                                **/
/** 			THE MAIN FUNCTION Lasso                    **/
/** ______________________________________________________________ **/
/** ______________________________________________________________ **/
/**                                                                **/
/** RÉSOLUTION APPROCHÉE DE LA VRAISEMBLANCE PÉNALISÉE             **/
/** Problème de Meinshausen & Bulhmann                             **/
/** ______________________________________________________________ **/
/**
 ** in :
 **	p    = current column
 **     rule = boolean. If zero, AND rule is performed; OR rule is
 **            performed for any other value
 **	S    = matrix S
 **	Rho  = penalty mask Rho
 ** 	eps  = seuil de convergence pour le Lasso
 **
 ** out :
 **	Theta   = matrice inferrée
 **/

void Lasso( int *p,
	    double S[*p][*p],
	    double Rho[*p][*p],
	    double *eps,
	    double Beta[*p][*p] ) {
  
  /** INITIALISATION**/
  int i,k, active, Sp[*p-1];
  double S11[*p-1][*p-1];
  double rho12[*p-1], s12[*p-1], beta[*p-1];
  
  for(k=0; k < *p; k++) {
    /** Setup_Lasso prépare les vecteurs pour la valeur courante de k  **/
    /** S11, rho12, s12, beta **/
    Setup_Lasso(k, *p ,S, Rho,  Beta, s12, rho12,  beta, S11, Sp, &active);

    /** Lasso via pathwise coordinate **/
    Pathwise_Coordinate(*p-1, S11, s12, rho12, *eps / *p, beta, Sp, &active);
    
    /** Reconstruction of the kth column of Beta **/
    for(i=0; i<k; i++) {
      Beta[k][i] = beta[i];
    }
    for(i=k+1; i<*p; i++) {
      Beta[k][i] = beta[i-1];
    }
  }
}

/** ______________________________________________________________ **/
/** ______________________________________________________________ **/
/**                                                                **/
/** 			THE MAIN FUNCTION ARLasso                  **/
/** ______________________________________________________________ **/
/** ______________________________________________________________ **/
/**                                                                **/
/** RÉSOLUTION DE LA VRAISEMBLANCE PÉNALISÉE MODÈLE AUTORÉGRESSIF  **/
/** ______________________________________________________________ **/
/**
 ** in :
 **	p    = current column
 **	S    = matrix S
 **	V    = matrix V
 **	Rho  = penalty mask Rho
 ** 	eps  = seuil de convergence pour le Lasso
 **
 ** out :
 **	A   = matrice inferrée
 **/

void ARLasso( int *p,
	      double S[*p][*p],
	      double V[*p][*p],
	      double Rho[*p][*p],
	      double *eps,
	      double A[*p][*p] ) {
  
  int i,k, active, Sp[*p];
  double Rhok[*p], Vk[*p], beta[*p];

  for(k=0; k<*p; k++) {
    /** Prépare les vecteurs S, Rhok et Vk pour la colonne courante **/
    Setup_ARLasso(k, *p ,V, Rho, A, Vk, Rhok, beta, Sp, &active);
    
    /** Lasso via pathwise coordinate **/
    Pathwise_Coordinate(*p, S, Vk, Rhok, *eps / *p, beta, Sp, &active);
    
    /** Reconstruction of the kth column of A **/
    for(i=0; i<*p; i++) {
      A[k][i] = beta[i];
    }
  }
}
