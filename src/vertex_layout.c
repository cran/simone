/*
######################################################################
#
# vertex_layout.c
#
# This file contain the routine to compute of vertex positions.
# Extracted from R/sna package (routine initially written by 
# Carter T. Butts <buttsc@uci.edu> (2004)).
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 12/19/04
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
#
######################################################################
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>

void vertex_coord_C(int *d, double *pn, int *pniter, 
		    double *pmaxdelta, double *pvolume, 
		    double *pcoolexp, double *prepulserad,
		    double *x, double *y)
/*
Calculate a two-dimensional Fruchterman-Reingold layout for (symmetrized) 
adjacency matrix d.  Positions (stored in (x,y)) should be initialized
prior to calling this routine.
*/
{
  double frk,maxdelta,volume,coolexp,repulserad,t,ded,xd,yd,*dx,*dy;
  double rf,af;
  long int n,j,k;
  int niter,i;

  /*Define various things*/
  n          = (long int)*pn;
  niter      = *pniter;
  maxdelta   = *pmaxdelta;
  volume     = *pvolume;
  coolexp    = *pcoolexp;
  repulserad = *prepulserad;
  frk        = sqrt(volume/(double)n); /*Define the F-R constant*/

  /*Allocate memory for transient structures*/
  dx=(double *)R_alloc(n,sizeof(double));
  dy=(double *)R_alloc(n,sizeof(double));

  /*Run the annealing loop*/
  for( i=niter; i>=0; i--) {

    /*Set the temperature (maximum move/iteration)*/
    t = maxdelta*pow(i/(double)niter, coolexp);

    /*Clear the deltas*/
    for(j=0;j<n;j++){
      dx[j]=0.0;
      dy[j]=0.0;
    }

    /*Increment deltas for each undirected pair*/
    for(j=0;j<n;j++)
      for(k=j+1;k<n;k++){
        /*Obtain difference vector*/
        xd=x[j]-x[k];
        yd=y[j]-y[k];
        ded=sqrt(xd*xd+yd*yd);  /*Get dyadic euclidean distance*/
        xd/=ded;                /*Rescale differences to length 1*/
        yd/=ded;
        /*Calculate repulsive "force"*/
        rf=frk*frk*(1.0/ded-ded*ded/repulserad);
        dx[j]+=xd*rf;        /*Add to the position change vector*/
        dx[k]-=xd*rf;
        dy[j]+=yd*rf;
        dy[k]-=yd*rf;
        /*Calculate the attractive "force"*/
        if(d[j+k*n]||d[k+j*n]){
          af=ded*ded/frk;
          dx[j]-=xd*af;        /*Add to the position change vector*/
          dx[k]+=xd*af;
          dy[j]-=yd*af;
          dy[k]+=yd*af;
        }
      }
    /*Dampen motion, if needed, and move the points*/
    for(j=0;j<n;j++){
      ded=sqrt(dx[j]*dx[j]+dy[j]*dy[j]);
      if(ded>t){                 /*Dampen to t*/
        ded=t/ded;
        dx[j]*=ded;
        dy[j]*=ded;
      }
      x[j]+=dx[j];               /*Update positions*/
      y[j]+=dy[j];
    }
  }
}

