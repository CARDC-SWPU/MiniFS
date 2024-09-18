#include <bits/stdc++.h>
#include <omp.h>
#include <metis.h>
#ifdef CILK
#include <cilk/cilk.h>
#endif

#include "TDTree.h"
#include "TDTree_permutations.h"
#include "TDTree_creation.h"
#include "TDTree_partitioning.h"

// Create a nodal graph from a item mesh for METIS
void create_nodal_graph (int *graphIndex, int **graphValue, int *i2n, int localNbItem, int dimItem, int localNbNode)
{
    unordered_map<int, bool> *n2n = new unordered_map<int, bool>[localNbNode]();

    for (int i = 0; i < localNbItem; i++)
    {
        for (int j = 0; j < dimItem; j++)
        {
            if (j && i2n[i*dimItem+j] == i2n[i*dimItem]) break;
            for (int k = 0; k < dimItem; k++)
            {
                if (k && i2n[i*dimItem+k] == i2n[i*dimItem]) break;
                n2n[i2n[i*dimItem+j]][i2n[i*dimItem+k]] = 1;
                n2n[i2n[i*dimItem+k]][i2n[i*dimItem+j]] = 1; 
            }
        }
    }

    graphIndex[0] = 0;
    for (int i = 0; i < localNbNode; i++) 
        graphIndex[i+1] = graphIndex[i] + n2n[i].size();
    (*graphValue) = new int[graphIndex[localNbNode]];
    for(int i = 0, j = 0; i < localNbNode; i++)
    {
        for (auto tmp:n2n[i])
            (*graphValue)[j++] = tmp.first;
        n2n[i].clear();
    }
    delete[] n2n;
}

// Create local itemToNode array containing elements indexed contiguously from 0 to
// localNbItem and return the number of nodes
int create_local_i2n (int *local_i2n, int *i2n, int firstItem, int lastItem, int dimItem)
{
    int localNbNode = 0;
    unordered_map<int, int> newNode;

    for (int i = firstItem*dimItem, j = 0; i < (lastItem+1)*dimItem; i++, j++) {
        int oldNode = i2n[i];
        if (!newNode.count(oldNode)) {
            newNode[oldNode] = localNbNode++;
        }
        local_i2n[j] = newNode[oldNode];
    }
    newNode.clear();
    return localNbNode;
}

// partitioning of isolators with more than PARTSIZE item
int TDTree::TDTree_partitioning (int *i2n, int *local_i2n, int dimItem, int firstItem, int lastItem, int *nodePart)
{
    int localNbItem = lastItem - firstItem + 1;
    int localNbNode = create_local_i2n (local_i2n, i2n, firstItem, lastItem, dimItem);

    int constraint = 1, objVal;
    int *graphIndex = new int [localNbNode + 1](), *graphValue;
    
    create_nodal_graph (graphIndex, &graphValue, local_i2n, localNbItem, dimItem, localNbNode);
      
    int parts = nbParts;

#ifdef OMP
    #pragma omp critical
#endif
    METIS_PartGraphRecursive (&localNbNode, &constraint, graphIndex, graphValue,
                        NULL, NULL, NULL, &parts, NULL, NULL,
                        NULL, &objVal, nodePart);

    delete[] graphValue, delete[] graphIndex;

    return parts;
}