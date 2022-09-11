/*******************************************************************
  Copyright (C) 2001-2012 Leo Breiman, Adele Cutler and Merck & Co., Inc.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
*******************************************************************/
  
  #include <R.h>
  #include "rf.h"
  #include "regTree.c"


  void predictRegTree(double *x, int nsample, int mdim,
                      int *lDaughter, int *rDaughter, int *nodestatus,
                      double *ypred, double *split, double *nodepred,
                      int *splitVar, int treeSize, int *cat, int maxcat,
                      int *nodex) {
    int i, j, k, m, *cbestsplit;
    double dpack;
    
    /* decode the categorical splits */
    if (maxcat > 1) {
      cbestsplit = (int *) Calloc(maxcat * treeSize, int);
      zeroInt(cbestsplit, maxcat * treeSize);
      for (i = 0; i < treeSize; ++i) {
        if (nodestatus[i] != NODE_TERMINAL && cat[splitVar[i] - 1] > 1) {
          dpack = split[i];
          /* unpack `npack' into bits */
          /* unpack(dpack, maxcat, cbestsplit + i * maxcat); */
          for (j = 0; j < cat[splitVar[i] - 1]; ++j) {
            cbestsplit[j + i*maxcat] = ((unsigned long) dpack & 1) ? 1 : 0;
            dpack = dpack / 2.0 ;
            /* cbestsplit[j + i*maxcat] = npack & 1; */
          }
        }
      }
    }
    
    for (i = 0; i < nsample; ++i) {
      k = 0;
      while (nodestatus[k] != NODE_TERMINAL) { /* go down the tree */
            m = splitVar[k] - 1;
        if (cat[m] == 1) {
          k = (x[m + i*mdim] <= split[k]) ?
          lDaughter[k] - 1 : rDaughter[k] - 1;
        } else {
          /* Split by a categorical predictor */
          k = cbestsplit[(int) x[m + i * mdim] - 1 + k * maxcat] ?
          lDaughter[k] - 1 : rDaughter[k] - 1;
        }
      }
      /* terminal node: assign prediction and move on to next */
      ypred[i] = nodepred[k];
      nodex[i] = k + 1;
    }
    if (maxcat > 1) Free(cbestsplit);
  }
  
  
  
  /*
   Train a regression randomforest model
   */
  void regRF_MSE(double *x, double *y, int *xdim, int *sampsize,
             int *nthsize, int *nrnodes, int *nTree, int *mtry, int *imp,
             int *cat, int *maxcat, int *jprint, int *doProx, int *oobprox,
             int *biasCorr, double *yptr, double *errimp, double *impmat,
             double *impSD, double *prox, int *treeSize, int *nodestatus,
             int *lDaughter, int *rDaughter, double *avnode, int *mbest,
             double *upper, double *mse, int *keepf, int *replace,
             int *testdat, double *xts, int *nts, double *yts, int *labelts,
             double *yTestPred, double *proxts, double *msets, double *coef,
             int *nout, int *inbag, double *Sw) {
    /*************************************************************************
     Input:
     mdim=number of variables in data set
     nsample=number of cases
     nthsize=number of cases in a node below which the tree will not split,
     setting nthsize=5 generally gives good results.
     nTree=number of trees in run.  200-500 gives pretty good results
     mtry=number of variables to pick to split on at each node.  mdim/3
     seems to give genrally good performance, but it can be
     altered up or down
     imp=1 turns on variable importance.  This is computed for the
     mth variable as the percent rise in the test set mean sum-of-
     squared errors when the mth variable is randomly permuted.
     *************************************************************************/
    
    double errts = 0.0, averrb, meanY, meanYts, varY, varYts, r, xrand,
      errb = 0.0, resid=0.0, ooberr, ooberrperm, delta, *resOOB;
    
    double *yb, *xtmp, *xb, *ytr, *ytree, *tgini, *sw;
    
    int k, m, mr, n, nOOB, j, jout, idx, ntest, last, ktmp, nPerm,
    nsample, mdim, keepF, keepInbag;
    int *oobpair, varImp, localImp, *varUsed;
    
    int *in, *nind, *nodex, *nodexts;
    
    nsample = xdim[0];
    mdim = xdim[1];
    ntest = *nts;
    varImp = imp[0];
    localImp = imp[1];
    nPerm = imp[2];
    keepF = keepf[0];
    keepInbag = keepf[1];
    
    if (*jprint == 0) *jprint = *nTree + 1;
    
    sw         = (double *) S_alloc(mdim, sizeof(double));
    yb         = (double *) S_alloc(*sampsize, sizeof(double));
    xb         = (double *) S_alloc(mdim * *sampsize, sizeof(double));
    ytr        = (double *) S_alloc(nsample, sizeof(double));
    xtmp       = (double *) S_alloc(nsample, sizeof(double));
    resOOB     = (double *) S_alloc(nsample, sizeof(double));
    
    in        = (int *) S_alloc(nsample, sizeof(int));
    nodex      = (int *) S_alloc(nsample, sizeof(int));
    varUsed    = (int *) S_alloc(mdim, sizeof(int));
    nind = *replace ? NULL : (int *) S_alloc(nsample, sizeof(int));
    
    for (j = 0; j < mdim; ++j) {
      sw[j]=Sw[j];
    }
    
    if (*testdat) {
      ytree      = (double *) S_alloc(ntest, sizeof(double));
      nodexts    = (int *) S_alloc(ntest, sizeof(int));
    }
    oobpair = (*doProx && *oobprox) ?
    (int *) S_alloc(nsample * nsample, sizeof(int)) : NULL;
    
    /* If variable importance is requested, tgini points to the second
     "column" of errimp, otherwise it's just the same as errimp. */
    tgini = varImp ? errimp + mdim : errimp;
    
    averrb = 0.0;
    meanY = 0.0;
    varY = 0.0;
    
    zeroDouble(yptr, nsample);
    zeroInt(nout, nsample);
    for (n = 0; n < nsample; ++n) {
      varY += n * (y[n] - meanY)*(y[n] - meanY) / (n + 1);
      meanY = (n * meanY + y[n]) / (n + 1);
    }
    varY /= nsample;
    
    varYts = 0.0;
    meanYts = 0.0;
    if (*testdat) {
      for (n = 0; n < ntest; ++n) {
        varYts += n * (yts[n] - meanYts)*(yts[n] - meanYts) / (n + 1);
        meanYts = (n * meanYts + yts[n]) / (n + 1);
      }
      varYts /= ntest;
    }
    
    if (*doProx) {
      zeroDouble(prox, nsample * nsample);
      if (*testdat) zeroDouble(proxts, ntest * (nsample + ntest));
    }
    
    if (varImp) {
      zeroDouble(errimp, mdim * 2);
      if (localImp) zeroDouble(impmat, nsample * mdim);
    } else {
      zeroDouble(errimp, mdim);
    }
    if (*labelts) zeroDouble(yTestPred, ntest);
    
    /* print header for running output */
    if (*jprint <= *nTree) {
      Rprintf("     |      Out-of-bag   ");
      if (*testdat) Rprintf("|       Test set    ");
      Rprintf("|\n");
      Rprintf("Tree |      MSE  %%Var(y) ");
      if (*testdat) Rprintf("|      MSE  %%Var(y) ");
      Rprintf("|\n");
    }
    GetRNGstate();
    /*************************************
     * Start the loop over trees.
     *************************************/
   
    for (j = 0; j < *nTree; ++j) {
      idx = keepF ? j * *nrnodes : 0;
      zeroInt(in, nsample);
      zeroInt(varUsed, mdim);
      /* Draw a random sample for growing a tree. */
      /*Rprintf("The useweights flag was set to %d", *useweights);*/
      if (*replace) { /* sampling with replacement */
      for (n = 0; n < *sampsize; ++n) {
        xrand = unif_rand();
        k = xrand * nsample;
        in[k] = 1;
        yb[n] = y[k];
        for(m = 0; m < mdim; ++m) {
          xb[m + n * mdim] = x[m + k * mdim];
        }
      }
      } else { /* sampling w/o replacement */
      for (n = 0; n < nsample; ++n) nind[n] = n;
        last = nsample - 1;
        for (n = 0; n < *sampsize; ++n) {
          ktmp = (int) (unif_rand() * (last+1));
          k = nind[ktmp];
          swapInt(nind[ktmp], nind[last]);
          last--;
          in[k] = 1;
          yb[n] = y[k];
          for(m = 0; m < mdim; ++m) {
            xb[m + n * mdim] = x[m + k * mdim];
          }
        }
      }
      if (keepInbag) {
        for (n = 0; n < nsample; ++n) inbag[n + j * nsample] = in[n];
      }
      /* grow the regression tree */
      regTree(xb, yb, mdim, *sampsize, lDaughter + idx, rDaughter + idx,
              upper + idx, avnode + idx, nodestatus + idx, *nrnodes,
              treeSize + j, *nthsize, *mtry, mbest + idx, cat, tgini,
              varUsed, sw);
      /* predict the OOB data with the current tree */
      /* ytr is the prediction on OOB data by the current tree */
      predictRegTree(x, nsample, mdim, lDaughter + idx,
                     rDaughter + idx, nodestatus + idx, ytr, upper + idx,
                     avnode + idx, mbest + idx, treeSize[j], cat, *maxcat,
                     nodex);
      /* yptr is the aggregated prediction by all trees grown so far */
      errb = 0.0;
      ooberr = 0.0;
      jout = 0; /* jout is the number of cases that has been OOB so far */
      nOOB = 0; /* nOOB is the number of OOB samples for this tree */
      for (n = 0; n < nsample; ++n) {
        if (in[n] == 0) {
          nout[n]++;
          nOOB++;
          yptr[n] = ((nout[n]-1) * yptr[n] + ytr[n]) / nout[n];
          resOOB[n] = ytr[n] - y[n];
          ooberr += resOOB[n] * resOOB[n];
        }
        if (nout[n]) {
          jout++;
          errb += (y[n] - yptr[n]) * (y[n] - yptr[n]);
        }
      }
      errb /= jout;
      /* Do simple linear regression of y on yhat for bias correction. */
      
      /* predict testset data with the current tree */
      if (*testdat) {
        predictRegTree(xts, ntest, mdim, lDaughter + idx,
                       rDaughter + idx, nodestatus + idx, ytree,
                       upper + idx, avnode + idx,
                       mbest + idx, treeSize[j], cat, *maxcat, nodexts);
        /* ytree is the prediction for test data by the current tree */
        /* yTestPred is the average prediction by all trees grown so far */
        errts = 0.0;
        for (n = 0; n < ntest; ++n) {
          yTestPred[n] = (j * yTestPred[n] + ytree[n]) / (j + 1);
        }
        /* compute testset MSE */
        if (*labelts) {
          for (n = 0; n < ntest; ++n) {
            resid = *biasCorr ?
            yts[n] - (coef[0] + coef[1]*yTestPred[n]) :
            yts[n] - yTestPred[n];
            errts += resid * resid;
          }
          errts /= ntest;
        }
      }
      /* Print running output. */
      if ((j + 1) % *jprint == 0) {
        Rprintf("%4d |", j + 1);
        Rprintf(" %8.4g %8.2f ", errb, 100 * errb / varY);
        if(*labelts == 1) Rprintf("| %8.4g %8.2f ",
           errts, 100.0 * errts / varYts);
        Rprintf("|\n");
      }
      mse[j] = errb;
      if (*labelts) msets[j] = errts;
      

      
      /* Variable importance */
      if (varImp) {
        for (mr = 0; mr < mdim; ++mr) {
          if (varUsed[mr]) { /* Go ahead if the variable is used */
      /* make a copy of the m-th variable into xtmp */
      for (n = 0; n < nsample; ++n)
        xtmp[n] = x[mr + n * mdim];
            ooberrperm = 0.0;
            for (k = 0; k < nPerm; ++k) {
              permuteOOB(mr, x, in, nsample, mdim);
              predictRegTree(x, nsample, mdim, lDaughter + idx,
                             rDaughter + idx, nodestatus + idx, ytr,
                             upper + idx, avnode + idx, mbest + idx,
                             treeSize[j], cat, *maxcat, nodex);
              for (n = 0; n < nsample; ++n) {
                if (in[n] == 0) {
                  r = ytr[n] - y[n];
                  ooberrperm += r * r;
                  if (localImp) {
                    impmat[mr + n * mdim] +=
                      (r*r - resOOB[n]*resOOB[n]) / nPerm;
                  }
                }
              }
            }
            delta = (ooberrperm / nPerm - ooberr) / nOOB;
            errimp[mr] += delta;
            impSD[mr] += delta * delta;
            /* copy original data back */
            for (n = 0; n < nsample; ++n)
              x[mr + n * mdim] = xtmp[n];
          }
        }
      }
    }
    PutRNGstate();
    /* end of tree iterations=======================================*/
    
    if (*biasCorr) {  /* bias correction for predicted values */
    for (n = 0; n < nsample; ++n) {
      if (nout[n]) yptr[n] = coef[0] + coef[1] * yptr[n];
    }
    if (*testdat) {
      for (n = 0; n < ntest; ++n) {
        yTestPred[n] = coef[0] + coef[1] * yTestPred[n];
      }
    }
    }
    
    if (*doProx) {
      for (n = 0; n < nsample; ++n) {
        for (k = n + 1; k < nsample; ++k) {
          prox[nsample*k + n] /= *oobprox ?
          (oobpair[nsample*k + n] > 0 ? oobpair[nsample*k + n] : 1) :
          *nTree;
          prox[nsample * n + k] = prox[nsample * k + n];
        }
        prox[nsample * n + n] = 1.0;
      }
      if (*testdat) {
        for (n = 0; n < ntest; ++n)
          for (k = 0; k < ntest + nsample; ++k)
            proxts[ntest*k + n] /= *nTree;
      }
    }
    
    if (varImp) {
      for (m = 0; m < mdim; ++m) {
        errimp[m] = errimp[m] / *nTree;
        impSD[m] = sqrt( ((impSD[m] / *nTree) -
          (errimp[m] * errimp[m])) / *nTree );
        if (localImp) {
          for (n = 0; n < nsample; ++n) {
            impmat[m + n * mdim] /= nout[n];
          }
        }
      }
    }
    for (m = 0; m < mdim; ++m) tgini[m] /= *nTree;
  }
  