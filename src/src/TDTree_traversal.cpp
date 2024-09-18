#ifdef CILK
    #include <cilk/cilk.h>
#endif
#include <omp.h>
#include <bits/stdc++.h>
#include "TDTree_traversal.h"

using namespace std;

// Follow the TDTree forward method to execute the given function in parallel tasks
void TDTree_traversal_forward (void (*userSeqFctPtr)  (char **, TDTreeArgs *),
                     void (*userVecFctPtr)  (char **, TDTreeArgs *),
                     char **userArgs, TDTreeNode *treePtr)
{
    // If current node is a leaf, call the user function
    if (treePtr->isLeaf) 
    {
        // Initialize the TDTree arguments
        TDTreeArgs treeArgs;
        
		treeArgs.firstCell  = treePtr->firstCell;
		treeArgs.lastCell   = treePtr->lastCell;
        treeArgs.firstFace  = treePtr->firstFace;
        treeArgs.lastFace   = treePtr->lastFace;
        treeArgs.firstNode  = treePtr->firstNode;
        treeArgs.lastNode   = treePtr->lastNode;
        treeArgs.boundryCell = &treePtr->boundryCell;
        treeArgs.boundryNode = &treePtr->boundryNode;
        treeArgs.boundryNbCell = &treePtr->boundryNbCell;
        treeArgs.boundryNbNode = &treePtr->boundryNbNode;      

        treeArgs.boundryCellNo = &treePtr->boundryCellNo;
        treeArgs.boundryN = &treePtr->boundryN;

		// Call user sequential function
		userSeqFctPtr (userArgs, &treeArgs); 
		// Call user vectorized function if neccessary
        if (userVecFctPtr != NULL)        
            userVecFctPtr (userArgs, &treeArgs);      
    }
    else 
    {
        // Forward method: first execute sons and then execute isolator
        int nbParts = treePtr->nbParts;
#ifndef FORKJOIN
#ifdef OMP
#pragma omp taskloop default(shared)
		for (int i = 0; i < nbParts; i++)
#elif CILK
        cilk_for(int i = 0; i < nbParts; i++)
#endif
#else
        #pragma omp parallel for 
		for (int i = 0; i < nbParts; i++)
#endif
            if (treePtr->son[i] != nullptr)
                TDTree_traversal_forward (userSeqFctPtr, userVecFctPtr, userArgs, treePtr->son[i]);
#ifndef FORKJOIN
        // Synchronization
#ifdef OMP
        #pragma omp taskwait
#elif CILK
        cilk_sync;
#endif
#endif
        if (treePtr->iso != nullptr) 
            TDTree_traversal_forward (userSeqFctPtr, userVecFctPtr, userArgs, treePtr->iso);
    }
}

// Follow the TDTree backward method to execute the given function in parallel tasks
void TDTree_traversal_backward (void (*userSeqFctPtr)  (char **, TDTreeArgs *),
                     void (*userVecFctPtr)  (char **, TDTreeArgs *),
                     char **userArgs, TDTreeNode *treePtr)
{
    // If current node is a leaf, call the appropriate function
    if (treePtr->isLeaf) 
    {

        // Initialize the TDTree arguments
        TDTreeArgs treeArgs;
        
		treeArgs.firstCell  = treePtr->firstCell;
		treeArgs.lastCell   = treePtr->lastCell;
        treeArgs.firstFace  = treePtr->firstFace;
        treeArgs.lastFace   = treePtr->lastFace;
        treeArgs.firstNode  = treePtr->firstNode;
        treeArgs.lastNode   = treePtr->lastNode;
        treeArgs.boundryCell = &treePtr->boundryCell;
        treeArgs.boundryNode = &treePtr->boundryNode;
        treeArgs.boundryNbCell = &treePtr->boundryNbCell;
        treeArgs.boundryNbNode = &treePtr->boundryNbNode;

		// Call user sequential function
		userSeqFctPtr (userArgs, &treeArgs); 
		// Call user vectorized function if neccessary
        if (userVecFctPtr != NULL)        
            userVecFctPtr (userArgs, &treeArgs);        
    }
    else 
    {
        int nbParts = treePtr->nbParts;
        // Backward method: first execute isolator and then execute sons
        if (treePtr->iso != nullptr)
            TDTree_traversal_backward (userSeqFctPtr, userVecFctPtr, userArgs, treePtr->iso);

#ifndef FORKJOIN
        // Synchronization
#ifdef OMP
        #pragma omp taskwait
#elif CILK
        cilk_sync;
#endif
#endif

#ifndef FORKJOIN    
#ifdef OMP
#pragma omp taskloop default(shared)
		for (int i = 0; i < nbParts; i++)
#elif CILK
        cilk_for(int i = 0; i < nbParts; i++)
#endif
#else
        #pragma omp parallel for default(shared)
		for (int i = 0; i < nbParts; i++)
#endif
            if (treePtr->son[i] != nullptr)
                TDTree_traversal_backward (userSeqFctPtr, userVecFctPtr, userArgs, treePtr->son[i]);
    }
}

// Follow the TDTree forward method to execute the given function in parallel tasks
void TDTree_traversal_noDependence (void (*userSeqFctPtr)  (char **, TDTreeArgs *),
                     void (*userVecFctPtr)  (char **, TDTreeArgs *),
                     char **userArgs, TDTreeNode *treePtr)
{
    // If current node is a leaf, call the appropriate function
    if (treePtr->isLeaf) 
    {

        // Initialize the TDTree arguments
        TDTreeArgs treeArgs;
        
		treeArgs.firstCell  = treePtr->firstCell;
		treeArgs.lastCell   = treePtr->lastCell;
        treeArgs.firstFace  = treePtr->firstFace;
        treeArgs.lastFace   = treePtr->lastFace;
        treeArgs.firstNode  = treePtr->firstNode;
        treeArgs.lastNode   = treePtr->lastNode;

		// Call user sequential function
		userSeqFctPtr (userArgs, &treeArgs);  

		// Call user vectorized function if neccessary
        if (userVecFctPtr != NULL)        
            userVecFctPtr (userArgs, &treeArgs);      
    }
    else {
        
        int nbParts = treePtr->nbParts;
        // noDependence method: execute sons and isolator at the same time
#ifndef FORKJOIN
#ifdef OMP
#pragma omp taskloop default(shared)
		for (int i = 0; i < nbParts; i++)
#elif CILK
        cilk_for(int i = 0; i < nbParts; i++)
#endif
#else
        #pragma omp parallel for 
		for (int i = 0; i < nbParts; i++)
#endif
            if (treePtr->son[i] != nullptr)
                TDTree_traversal_noDependence (userSeqFctPtr, userVecFctPtr, userArgs, treePtr->son[i]);
        
        if (treePtr->iso != nullptr) 
            TDTree_traversal_noDependence (userSeqFctPtr, userVecFctPtr, userArgs, treePtr->iso);
    }
}

// Wrapper used to parallel user function
void TDTree::task_traversal (void (*userSeqFctPtr)  (char **, TDTreeArgs *), 
                                void (*userVecFctPtr)  (char **, TDTreeArgs *), 
                                char **userArgs, int traversal_type)
{
#ifndef FORKJOIN
#ifdef OMP
        #pragma omp parallel
        #pragma omp single nowait
#endif
#endif
    if (traversal_type == Forward)
        TDTree_traversal_forward (userSeqFctPtr, userVecFctPtr, userArgs, treeRoot);
    else 
        if (traversal_type == Backward)
            TDTree_traversal_backward (userSeqFctPtr, userVecFctPtr, userArgs, treeRoot);
        else
            TDTree_traversal_noDependence (userSeqFctPtr, userVecFctPtr, userArgs, treeRoot);         
}
