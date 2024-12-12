// main.cpp
//

#include "mex.h"

#include "iagm/iagm_ves.h"
#include<iostream>


using namespace std;

// Input: Configuration and configurational velocity
// Output: New configurational velocity
//
void getCollision(double n[], size_t m, double x1[], double x2[], double min_sep, size_t nbr, double x_final[], double vol[], size_t np, double ids[], double vols[], double tol, size_t nv, size_t nb, size_t Npv, size_t Npb, size_t nexten, double cellSize)
{

  //std::cout << "Vesicles: " << m << std::endl;
  
  Iagm_ves *simulation = new Iagm_ves(n,m,x1,x2,min_sep,nbr,np,tol,nv,nb,Npv,Npb,nexten,cellSize);
  //simulation->updateQ();
  //std::cout<<"simulation finish"<<std::endl;
  simulation->getVolume();
  size_t nCount = 0;
  for(size_t i=0, currPt=0;i<m;i++)
  {
    size_t nn = (size_t)n[i];
    for(size_t j=0;j<nn;j++,currPt++)
    {
      x_final[2*nCount+j] = simulation->m_volumeGra[2*currPt];
      x_final[2*nCount+j+nn] = simulation->m_volumeGra[2*currPt+1];
      
      ids[2*nCount+j] = simulation->m_ids[2*currPt];
      ids[2*nCount+j+nn] = simulation->m_ids[2*currPt+1];
      
      vols[2*nCount+j] = simulation->m_vols[2*currPt];
      vols[2*nCount+j+nn] = simulation->m_vols[2*currPt+1];

    }
    nCount = nCount + nn;

  } 
  vol[0] = simulation->m_volume;
  delete simulation;
  //std::cout<<"exiting getCollision"<<std::endl;
}

void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
  double *x1, *x2, *x_final,min_sep,*n,cellSize;
  size_t m,nbr,np;
  size_t nv,nb,Npv,Npb,nexten;
  double *ivolume;
  double tol;

  double *ids, *vols;
/*
  if (nrhs != )
	  mexErrMsgTxt("Three inputs required.");
  
  mwSize n = mxGetM(prhs[0]); // Number of DOFs
  plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
  
  mwSize m = mxGetN(prhs[2]); // Number of vesicles

  x = mxGetPr(prhs[0]);
  v = mxGetPr(prhs[1]);
  e = mxGetPr(prhs[2]);

  vNew = mxGetPr(plhs[0]);
*/
  n = mxGetPr(prhs[0]);	//number of points on each surface
  m = (size_t)mxGetScalar(prhs[1]);	//number of surfaces 
  x1 = mxGetPr(prhs[2]);
  x2 = mxGetPr(prhs[3]);
  min_sep = mxGetScalar(prhs[4]);
  nbr = (size_t)mxGetScalar(prhs[5]);
  np = (size_t)mxGetScalar(prhs[6]);	//total number of points on all surface
  tol = mxGetScalar(prhs[7]);
  nv = (size_t)mxGetScalar(prhs[8]);
  nb = (size_t)mxGetScalar(prhs[9]);
  Npv = (size_t)mxGetScalar(prhs[10]);
  Npb = (size_t)mxGetScalar(prhs[11]);
  nexten = (size_t)mxGetScalar(prhs[12]);
  cellSize = mxGetScalar(prhs[13]);
  //n=malloc(m);
  //for(int i=0;i<m;i++){
  //	n[i] = (size_t)nn[i];
	  //std::cout<<n[i]<<std::endl;
  //}

  // Store the total gradient
  plhs[0] = mxCreateDoubleMatrix(2*np,1,mxREAL);
  // Store the total volume
  plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
  x_final = mxGetPr(plhs[0]);
  ivolume = mxGetPr(plhs[1]);
  
  // The index vector which tells which volume does this gradient element belongs
  plhs[2] = mxCreateDoubleMatrix(2*np,1,mxREAL);
  ids = mxGetPr(plhs[2]);
  // The volume vector which tells the volume of the gradient element belongs
  plhs[3] = mxCreateDoubleMatrix(2*np,1,mxREAL);
  vols = mxGetPr(plhs[3]);
  
  getCollision(n, m, x1, x2, min_sep, nbr, x_final, ivolume, np, ids, vols, tol, nv, nb, Npv, Npb,nexten,cellSize);
  //free(n);
  //std::cout<<"exitting mex"<<std::endl;

}

