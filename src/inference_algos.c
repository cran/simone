// Compile with :  R CMD SHLIB glasso.c
#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <R_ext/Utils.h>

/** ______________________________________________________________ **/
/**                                                                **/
/** INITIALISATION DES VECTEURS POUR UNE COLONNE FIXE m            **/
/** Création de W11, Rho12, S12                                    **/
/**
 ** in :
 **	m   = numéro de la colonne considérée
 **	p   = nombre de noeuds
 **	S   = matrice de covariances
 **	Rho = matrice de pénalités
 ** 	W   = Matrice à optimiser
 **
 ** out :
 **	s   = vecteur s12
 **	r   = vecteur rho12
 **	V   = matrice W11
 **/
void Setup(int m,
	   int p,
	   double S[p][p],
	   double Rho[p][p],
	   double W[p][p],
	   double s[p-1],
	   double r[p-1],
	   double V[p-1][p-1]) {

  int i,j;
  
  // i<m
  for(i=0; i < m; i++){
      r[i] = Rho[i][m];	
      s[i] = S[i][m];	

      // j<m
      for(j=0; j < m; j++) {
        V[i][j] = W[i][j]/2; 
      }
      // j>m
      for(j=m+1; j < p; j++) {
        V[i][j-1] = W[i][j]/2; 
      }
  } 

  // i>m
  for(i=m+1; i < p; i++){
      r[i-1] = Rho[i][m];
      s[i-1] = S[i][m];	

      // j<m
      for(j=0; j < m; j++) {
        V[i-1][j] = W[i][j]/2; 
      }
      // j>m
      for(j=m+1; j < p; j++) {
        V[i-1][j-1] = W[i][j]/2; 
      }
  } 

}

/** ______________________________________________________________ **/
/**                                                                **/
/** S(x,y) = sign(x)(|x|-y)                                        **/
/** Soft Threshold Operator                                        **/
/**								   **/
double Sfunc(double x,double y) {
  if (y > 0) {		// si pénalité <0, alors 
			//    pénalité -> +Inf => renvoi de 0
    if (y < fabs(x)) {	// si pénalité >= |x|, alors renvoi de 0
      if (x > 0) {
	return(x-y);
      } else {
	return(x+y);
      }
    } 
  } else if (y==0) {	// pénalité nulle => solution des moindres carrés
    return(x);
  }
  return(0);
}



/** ______________________________________________________________ **/
/**                                                                **/
/** PRODUIT SCALAIRE partiel DE DEUX VECTEURS                      **/
/**								   **/
double scalar_product(	int n,
		      	double vecA[n],
		      	double vecB[n]) {
  int i;
  double result;
  
  result=0;
  for(i=0; i<n; i++) {
    result += vecA[i]*vecB[i];
  }
  return(result);
}

/** ______________________________________________________________ **/
/**                                                                **/
/** MAXIMUM DE DEUX SCALAIRES			                   **/
/**								   **/
double theMax(double a,double b) {
  if(a > b) {return(a);} else {return(b);}
}

/** ______________________________________________________________ **/
/**                                                                **/
/** AND RULE TO SYMMETRIZE			                   **/
/**								   **/
double andRule(double a,double b) {
  if(fabs(a) < fabs(b)) {return(a);} else {return(b);}
}

/** ______________________________________________________________ **/
/**                                                                **/
/** OR RULE TO SYMMETRIZE			                   **/
/**								   **/
double orRule(double a,double b) {
  if(fabs(a) > fabs(b)) {return(a);} else {return(b);}
}

/** ______________________________________________________________ **/
/**                                                                **/
/** LASSO ALGORITHM IN THE PATHWISE COORDINATE WAY                 **/
/**								   **/
/**
 ** in :
 **	p   = nombre de noeuds - 1
 **	r   = vecteur de pénalités rho12
 **	V   = matrice W11
 **	thr = threshold de convergence
 **	M   = vecteur d'indices d'éléments à ne pas skipper
 **	m   = indice d'élément à skipper 
 **
 ** in&out :
 **	s   = vecteur s12 de taille p
 **	beta= vecteur beta
 **/
void LassoPath(int p,
	       double s[p],
	       double r[p],
	       double V[p][p],
	       double thr,
	       double beta[p],
	       int M[p],
	       int m) {

  const int maxIt=1000;
  int j,k,it;
  double dlx,del,betaOld,Sx[p];

  /** calcul de produit scalaire **/  
  /*  Sx[indi] = s[indi] - dot_product(V[m[1:j]][indi],x[M[1:j]]); */
  for(j=0; j<p; j++){
    Sx[j] = s[j];
    for (k=0; k < m; k++) {
      if (j != M[k]) {
	Sx[j] -= V[j][M[k]] * beta[M[k]];
      }
    }
  }


  /** _____________________________________________________**/
  /**                                                      **/
  /** PATHWISE COORDINATE ALGORITHM                        **/
  dlx = thr+1;
  //  printf("\n %f\n",thr);
  it = 0;
  while( (dlx >= thr) && (it < maxIt) ) {
    dlx = 0;
    for(j=0; j<p; j++) {

      /** copy of old value of beta **/
      betaOld = beta[j];

      /** the S update function**/
      beta[j] = Sfunc(Sx[j],r[j])/V[j][j]; // DIV/0 ??????????

      /** compute the difference to check convergence **/
      if(beta[j] != betaOld) {
      	del = beta[j]-betaOld;
	dlx = theMax(dlx,fabs(del));
	for(k=0;k < p; k++) {
      	  if (k != j) {
	    Sx[k] -= del * V[j][k];
	  }
	}
      }
    }
    
    R_CheckUserInterrupt();
    it++;
  }
  /*if (it == maxIt) {
    printf("\n !! convergence problem in pathwise coordinate !!");
  }*/
}

/** ______________________________________________________________ **/
/**                                                                **/
/** CLEAN UP FUNCTION                                              **/
/** Reconstruit W et S à partir de leurs sous-parties		   **/
/**								   **/
/**
 ** in :
 **	k   = colonne en cours
 **	p   = nombre de noeuds
 **	V   = matrice W11
 **	beta= vecteur de coefficients
 **	M   = vecteur d'indices d'éléments à ne pas skipper (astuce numérique)
 **	m   = indice d'élément à skipper (astuce numérique)
 **
 ** out :
 **	w   = matrice W achevée
 **	s   = vecteur s12
 **/
void Cleanup(int k,
	     int p,
	     double w[p][p],
	     double V[p-1][p-1],
	     double beta[p-1],
	     int M[p-1],
	     int m) {

  int i,j;
  double tmp[p-1]; 

  /** calcul de produit scalaire indicé **/
  //  s[i] = V[i][m[1:j]] %*% beta[M[1:j]]
  for(i=0; i<p-1; i++){

    tmp[i] = 0;
    for (j=0; j < p-1; j++) {
      tmp[i] += V[j][i] * beta[j];
    }

  }

  /** Mise à jour de la matrice estimée W **/
  for (i=0; i<k; i++) {
    w[k][i] = tmp[i];		// save p values of W[k,]
    w[i][k] = w[k][i];		// save p values of W[,k] (symmetry)
  } 
  for (i=k+1; i<p; i++) {
    w[k][i] = tmp[i-1];		// save p values of W[k,] 
    w[i][k] = w[k][i];		// save p values of W[,k] (symmetry)
  }  

}

/** ______________________________________________________________ **/
/**                                                                **/
/** COMPUTE THETA BY BLOCKWISE INVERTION OF W                      **/
/**								   **/
/**
 ** in :
 **	p   = colonne en cours
 **	W   = matrice W
 **	beta  = vecteur de coefficients
 **
 ** out :
 **	Theta = matrice W inversée
 **/
void Inverse(int *p,
	     double W[*p][*p],
	     double beta[*p][*p-1],
	     double Theta[*p][*p]) {
  
  int i,j,k;
  double beta12[*p-1];
  double w12[*p-1];

  /** CONSTRUCTION DE LA MATRICE THETA **/
  for(j=0;j<*p;j++) {

    for(i=0;i<(*p-1);i++) {
      beta12[i] = -0.5 * beta[j][i];
    }

    k=0;
    for(i=0;i<*p;i++) {
      if (i != j) {w12[k] = W[j][i];}
      k++;
    }

    Theta[j][j] = 1/(W[j][j]+scalar_product(*p-1,beta12,w12));
    
    k=0;
    for(i=0; i<(*p-1);i++) {
      if (i < j) {

	Theta[j][i] = Theta[j][j] * beta12[i];
	/** Force la symmétrie des valeurs (erreur d'arrondi par invertion) **/
	Theta[i][j] = orRule(Theta[j][i],Theta[i][j]);
	Theta[j][i] = Theta[i][j];

      } else {

	Theta[j][i+1] = Theta[j][j] * beta12[i];
	/** Force la symmétrie des valeurs (erreur d'arrondi par invertion) **/
	Theta[i+1][j] = orRule(Theta[j][i+1],Theta[i+1][j]);
	Theta[j][i+1] = Theta[i+1][j];

      }
    }
  }
}


/** ______________________________________________________________ **/
/**                                                                **/
/** NormMat			                                   **/
/** Calculer la norme 1 d'une matrice				   **/
/**								   **/
double NormMat(int p, double S[p][p], int diag) {
  
  int i,j;
  double res;

  res = 0;
  if (diag==0) {  
    for(i=0; i<p; i++) {
      for(j=i+1; j<p; j++) {
        res = res + 2 * fabs(S[j][i]);
      }
    }

  } else {

    for(i=0; i<p; i++) {
      for(j=i+1; j<p; j++) {
        res = res + 2 * fabs(S[j][i]);
      }
      res = res + fabs(S[i][i]); 
    }

  }

  return(res);
}

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
	    double betaS[*p][*p-1],
	    int *finalIt) {
  
  /** INITIALISATION**/
  int i,k,m,indSparse[*p-1];
  double W11[*p-1][*p-1],WS[*p];
  double rho12[*p-1],s12[*p-1],beta[*p-1];
  double threshold, dlx, diff;

  threshold = *eps * NormMat(*p,W,0) / *p; // -> seuil = fn°(matrice)
  
  dlx = 2* threshold; // pour rentrer au moins une fois dans le while...
  
  while( (*finalIt < *maxIt) && (dlx > threshold) ) {
    dlx = 0.0;
    for(k=0; k<*p; k++) { // on va travailler sur la colonne k de W
      
      /** _____________________________________________________**/
      /**                                                      **/
      /** PRETRAITEMENTS                                       **/

      /** Saving current values of the kth column of W and beta...**/
      m=0;
      for(i=0; i<*p-1; i++){
	WS[i]   = W[k][i];	// -> save p-1 first values of W[,k]
        beta[i] = betaS[k][i];  // -> save all p-1 values of betaS[k,]

	/** Prépare le terme du Lasso en multiplication "sparse" **/
	if(fabs(beta[i]) > 0) {
	  indSparse[m] = i;
	  m++;
	}
      }
      WS[*p-1] = W[k][*p-1];	// -> save last value of W[,k]
      
      /** Setup prépare les vecteurs pour la valeur courante de k **/
      /** Crée W11, rho12 et s12 **/
      Setup(k, *p, S, Rho, W, s12, rho12, W11);
      
      /** Lasso via pathwise coordinate **/
      LassoPath(*p-1, s12, rho12, W11, 
		threshold/NormMat(*p-1,W11,1), 
		beta, indSparse, m);

      /** Reconstruction of the full matrix W **/
      Cleanup(k, *p, W, W11, beta, indSparse, m);

      /** Checking for convergence twixt old W[,k] and new W[,k]... **/
      diff = 0;
      for(i=0; i<*p-1; i++) {
	diff += fabs(W[k][i]-WS[i]); 	// p-1 premiers termes
	betaS[k][i] = beta[i];			// p-1 termes
      }
      diff += fabs(W[k][*p-1]-WS[*p-1]); 	// dernier terme
      dlx = theMax(dlx,diff);

      R_CheckUserInterrupt();

    }
    *finalIt = *finalIt + 1;

  }

} // end Glasso


/** ______________________________________________________________ **/
/** ______________________________________________________________ **/
/**                                                                **/
/** 			THE MAIN FUNCTION regLasso                 **/
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

void regLasso( int *p,
	      int *rule,
	      double S[*p][*p],
	      double Rho[*p][*p],
	      double *eps,
	      double Theta[*p][*p] ) {
  
  /** INITIALISATION**/
  int i,j,k,m,indSparse[*p-1];
  double W11[*p-1][*p-1];
  double rho12[*p-1],s12[*p-1],beta[*p-1];
    
  for(k=0; k < *p; k++) {
    

    /** Le vecteur beta contient la valeur de Theta12 **/
    j=0;
    m=0;
    for(i=0; i<*p; i++){
      if(i != k) {
	beta[j] = Theta[k][i];
	if(fabs(beta[i]) > 0) {
	  indSparse[m] = i;
	  m++;
	}
	j++;
      }
    }
    
    /** Setup prépare les vecteurs pour la valeur courante de k **/
    /** Crée W11, rho12 et s12 **/
    Setup(k,*p,S,Rho,S,s12,rho12,W11);
    
    /** Lasso via pathwise coordinate **/
    LassoPath(*p-1,s12,rho12,W11,*eps/NormMat(*p-1,W11,1),beta,indSparse,m);
    
    /** Reconstruction of the kth column of Theta **/
    for(i=0; i<k; i++) {
      Theta[k][i] = beta[i];
    }
    for(i=k+1; i<*p; i++) {
      Theta[k][i] = beta[i-1];
    }
    
  }

  // Post symmetrization of Theta...
  for (k=0; k<*p; k++) {
    for (i=0; i<*p; i++) {
      if (*rule == 0) {
	Theta[k][i] = andRule(Theta[i][k],Theta[k][i]);	
      } else {
	Theta[k][i] = orRule(Theta[i][k],Theta[k][i]);
      }
      Theta[i][k] = Theta[k][i];
    }
  }

}
