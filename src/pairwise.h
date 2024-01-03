/*
 * =====================================================================================
 *
 *       Filename:  pairwise.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/16/2023 13:00:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
/**
 * @author      : greg (greg@$HOSTNAME)
 * @file        : pairwise
 * @created     : Thursday Nov 16, 2023 13:00:26 EST
 */

#ifndef __PAIRWISE_H__
#define __PAIRWISE_H__

#include "bayes.h"

typedef struct pairwisedists 
    {
    PolyTree   *tree;
    MrBFlt     *dists;
/* int        *doubletCounts;   */   
    int        *clusters;
    int        nTaxa;
    int        nPairs;
    } PairwiseDists;

void CountFreqs(void);
int  CountDoublets(int nPairs);
void CountTriples(void);
int  triplePos(int i, int j, int k);

int DoEstQPairwise(void);
int DoPairwiseLogLike(void);
int DoTripletLogLike(void);

int DoPwSetParm(char *parmName, char *tkn); // param set -- see param list at end of command.c
void PrintPairwiseDists(PairwiseDists *pd);

#endif /* end of include guard PAIRWISE_H */

