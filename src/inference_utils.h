#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <R_ext/Utils.h>

/** ______________________________________________________________ **/
/**                                                                **/
/** INITIALISATION DES VECTEURS POUR UNE COLONNE FIXE m            **/
/** Création de W11, rho12, s12, beta + une copie de la            **/
/** colonne courante de W                                          **/
/**
 ** in :
 **	k     = numéro de la colonne considérée
 **	p     = nombre de noeuds
 **	S     = matrice de covariances
 **	Rho   = matrice de pénalités
 ** 	W     = Matrice à optimiser
 ** 	Beta = Valeur courante des beta
 **
 ** out :
 **	Smkk  = kieme colonne de S privé de la kieme ligne
 **	Rmkk  = kieme colonne de Rho privé de la kieme ligne
 **	Wmkk  = copie de la  kieme colonne de W
 **	beta  = valeur de départ pour beta dans le pathwise
 **	Wmkmk = matrice W privée de la kieme ligne et kieme colonne
 **	Sp    = vecteur d'indices des éléments non nul du beta courant
 **	active = nombre d'éléments non nul dans le beta courant
 **/
void Setup_GLasso(int k,
		  int p,
		  double S[p][p],
		  double Rho[p][p],
		  double W[p][p], 
		  double Beta[p][p],
		  double Smkk[p-1],
		  double Rmkk[p-1],
		  double Wmkk[p], 
		  double beta[p-1],
		  double Wmkmk[p-1][p-1],
		  int    Sp[p-1],
		  int    *active) {
  
  int i,j;
  *active = 0;

  // i<k
  for(i=0; i < k; i++){
    Rmkk[i] = Rho[k][i];
    Smkk[i] = S[k][i];
    Wmkk[i] = W[k][i];
    beta[i] = Beta[k][i];
    if (fabs(beta[i]) > 0) {
      Sp[*active] = i;
      *active=*active+1;
    }
        
    // j<k
    for(j=0; j < k; j++) {
      Wmkmk[j][i] = 0.5 * W[j][i];
    }
    // j>k
    for(j=k+1; j < p; j++) {
      Wmkmk[j-1][i] = 0.5 * W[j][i]; 
    }
  }

  Wmkk[k] = W[k][k];

  // i>k
  for(i=k+1; i < p; i++){
      Rmkk[i-1] = Rho[k][i];
      Smkk[i-1] = S[k][i];	
      Wmkk[i]   = W[k][i];
      beta[i-1] = Beta[k][i];
      if (fabs(beta[i-1]) > 0) {
	Sp[*active] = i-1;
	*active=*active+1;
      }

      // j<k
      for(j=0; j < k; j++) {
        Wmkmk[j][i-1] = 0.5 * W[j][i]; 
      }
      // j>k
      for(j=k+1; j < p; j++) {
        Wmkmk[j-1][i-1] = 0.5 * W[j][i]; 
      }
  }
}

/** ______________________________________________________________ **/
/**                                                                **/
/** INITIALISATION DES VECTEURS POUR UNE COLONNE FIXE k            **/
/**
 ** in :
 **	k     = numéro de la colonne considérée
 **	p     = nombre de noeuds
 **	S     = matrice de covariances
 **	Rho   = matrice de pénalités
 ** 	Beta  = Valeur courante des beta
 **
 ** out :
 **	Smkk   = kieme colonne de S privé de la kieme ligne
 **	Rmkk   = kieme colonne de Rho privé de la kieme ligne
 **	beta   = valeur de départ pour beta dans le pathwise
 **	Smkmk  = matrice S privée de la kieme ligne et kieme colonne
 **	Sp     = vecteur d'indices des éléments non nul du beta courant
 **	active = nombre d'éléments non nul dans le beta courant
 **/
void Setup_Lasso(int k,
		 int p,
		 double S[p][p],
		 double Rho[p][p],
		 double Beta[p][p],
		 double Smkk[p-1],
		 double Rmkk[p-1],
		 double beta[p-1],
		 double Smkmk[p-1][p-1],
		 int    Sp[p-1],
		 int    *active) {
  
  int i,j;
  *active = 0;

  // i<k
  for(i=0; i < k; i++){
    Rmkk[i] = Rho[k][i];
    Smkk[i] = S[k][i];
    beta[i] = Beta[k][i];
    if (fabs(beta[i]) > 0) {
      Sp[*active] = i;
      *active=*active+1;
    }
        
    // j<k
    for(j=0; j < k; j++) {
      Smkmk[j][i] = S[j][i];
    }
    // j>k
    for(j=k+1; j < p; j++) {
      Smkmk[j-1][i] = S[j][i]; 
    }
  }

  // i>k
  for(i=k+1; i < p; i++){
      Rmkk[i-1] = Rho[k][i];
      Smkk[i-1] = S[k][i];	
      beta[i-1] = Beta[k][i];
      if (fabs(beta[i-1]) > 0) {
	Sp[*active] = i-1;
	*active=*active+1;
      }

      // j<k
      for(j=0; j < k; j++) {
        Smkmk[j][i-1] = S[j][i]; 
      }
      // j>k
      for(j=k+1; j < p; j++) {
        Smkmk[j-1][i-1] = S[j][i]; 
      }
  }
}

/** ______________________________________________________________ **/
/**                                                                **/
/** INITIALISATION DES VECTEURS POUR UNE COLONNE FIXE k            **/
/**
 ** in :
 **	k     = numéro de la colonne considérée
 **	p     = nombre de noeuds
 **	V     = matrice de covariances temporelles
 **	Rho   = matrice de pénalités
 ** 	A     = Matrice à optimiser
 **
 ** out :
 **	Vk   = kieme colonne de V
 **	Rk   = kieme colonne de Rho
 **	beta = valeur de départ pour beta dans le pathwise
 **	Sp   = vecteur d'indices des éléments non nul du beta courant
 **	active = nombre d'éléments non nul dans le beta courant
 **/
void Setup_ARLasso(int k,
		   int p,
		   double V[p][p],
		   double Rho[p][p],
		   double A[p][p],
		   double Vk[p],
		   double Rk[p],
		   double beta[p],
		   int    Sp[p],
		   int    *active) {

  int i;
  
  *active = 0;
  for(i=0; i<p; i++){
    Rk[i] = Rho[k][i];
    Vk[i] = V[k][i];
    beta[i] = A[k][i];
    if (fabs(beta[i]) > 0) {
      Sp[*active] = i;
      *active = *active + 1;
    }
  }
  
}

/** ______________________________________________________________ **/
/**                                                                **/
/** S(x,y) = sign(x)(|x|-y)                                        **/
/** Soft Threshold Operator                                        **/
/**								   **/
double Sfunc(double x,double y) {
  if (y > 0) {		// check if penalty > 0
    if (y < fabs(x)) {	
      if (x > 0) {
	return(x-y);
      } else {
	return(x+y);
      }
    } 
  } else if (y==0) {	// penalty == 0 => least squares solution
    return(x);
  }
  return(0); // if penalty < 0, consider by convention infinite penalty : return 0
             // if penalty >= |x|, return 0
}

/** ______________________________________________________________ **/
/**                                                                **/
/** LASSO ALGORITHM IN THE PATHWISE COORDINATE WAY                 **/
/**								   **/
/**
 ** 
 ** résout le problème de pénalisation l1
 ** 
 ** 1/2 beta' X beta - beta' y + rho |beta|_1
 ** 
 ** in :
 **	p   = nombre de composant de beta
 **	X   = matrice de taille p x p
 **	y   = vecteur de taille p
 **	lambda = vecteur de pénalités
 **	thr = threshold de convergence
 **	Sp  = vecteur d'indices des éléments non nul du beta courant
 **	active = nombre d'éléments non nul dans le beta courant
 **
 ** in&out :
 **	beta= vecteur beta
 **/
void Pathwise_Coordinate(int p,
			 double X[p][p],
			 double y[p],
			 double lambda[p],
			 double thr,
			 double beta[p],
			 int not_null[p],
			 int *act) {

  const int maxIt=1e4;
  int i,k,it;
  double dlx,del,betaOld,Sx[p];

  /**  Sx[i] = y[i] - dot_product(X[i,:],beta); **/
  for(i=0; i < p; i++){
    Sx[i] = y[i];
    for (k=0; k < *act; k++) {
      if (i != not_null[k]) {
	Sx[i] -= X[i][not_null[k]] * beta[not_null[k]];
      }
    }
  }

  /** _____________________________________________________**/
  /**                                                      **/
  /** PATHWISE COORDINATE ALGORITHM                        **/
  dlx = thr+1;
  it = 0;
  while( (dlx >= thr) && (it < maxIt) ) {
    dlx = 0;
    for(i=0; i<p; i++) {
      
      /** copy of old value of beta **/
      betaOld = beta[i];
      
      /** the S update function**/
      beta[i] = Sfunc(Sx[i],lambda[i])/X[i][i];

      /** compute the difference to check convergence **/
      if(beta[i] != betaOld) {
      	del = beta[i]-betaOld;
	dlx = fmax(dlx,fabs(del));
	for(k=0; k< p; k++) {
      	  if(k != i) {
	    Sx[k] -= del * X[k][i];
	  }
	}
      }
    }
    
    R_CheckUserInterrupt();
    it++;
  }
  if (it >= maxIt) {
    //printf("\n Shooting did not converge after %d iterations !!",maxIt);
  }
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
 **
 ** out :
 **	w   = matrice W achevée
 **	s   = vecteur s12
 **/
void Cleanup(int k,
	     int p,
	     double w[p][p],
	     double Wmkmk[p-1][p-1],
	     double beta[p-1],
	     double Beta[p][p],
	     int not_null[p-1],
	     int *active) {

  int i,j;
  double tmp[p-1]; 

  /** calcul de produit scalaire indicé **/
  //  s[i] = V[i][:] %*% beta[:]
  for(i=0; i<p-1; i++){
    tmp[i] = 0;
    for (j=0; j < p-1; j++) {
	tmp[i] += Wmkmk[j][i] * beta[j] ;
    }

  }

  /** Mise à jour de la matrice estimée W **/
  for (i=0; i<k; i++) {
    w[k][i] = tmp[i];       // save p values of W[k,]
    w[i][k] = w[k][i]; 	    // save p values of W[,k] (symmetry)
    Beta[k][i] = beta[i];   // same for Beta */
  }
  for (i=k+1; i<p; i++) {
    w[k][i] = tmp[i-1];	    // save p values of W[k,]
    w[i][k] = w[k][i];	    // save p values of W[,k] (symmetry)
    Beta[k][i] = beta[i-1]; // same for Beta 
  }
  
}

/** ______________________________________________________________ **/
/**                                                                **/
/** NormMat			                                   **/
/** Calculer la norme 1 d'une matrice				   **/
/**								   **/
double NormMat(int p, double A[p][p]) {
  
  int i,j;
  double res;

  res = 0;
  for(i=0; i<p; i++) {
    for(j=i+1; j<p; j++) {
      res = res + 2 * fabs(A[j][i]);
    }
  }
  
  return(res);
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
