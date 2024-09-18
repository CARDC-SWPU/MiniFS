#ifdef CILK
#include <cilk/cilk.h>
#endif
#include <bits/stdc++.h>

#include "TDTree.h"
#include "TDTree_permutations.h"
#include "TDTree_partitioning.h"
#include "TDTree_creation.h"
#include "omp.h"

// Initialize the content of TDTree nodes
void TDTree_node_init (TDTreeNode *treePtr, int firstCell, int lastCell, int firstFace, int lastFace, int firstNode, int lastNode, 
                            int nbIsoItem, int *nbPartItem, bool isIso, bool isLeaf, int nbParts)
{
    treePtr->firstCell    = firstCell;
    treePtr->lastCell     = lastCell;
    treePtr->firstFace    = firstFace;
    treePtr->lastFace     = lastFace;
    treePtr->firstNode    = firstNode;
    treePtr->lastNode     = lastNode;
    treePtr->isIso        = isIso;
    treePtr->iso          = nullptr;
	treePtr->son		  = nullptr;
	treePtr->isLeaf		  = isLeaf;
    treePtr->nbParts      = nbParts;

    if (isLeaf == false) {
		
		treePtr->son = new TDTreeNode*[nbParts];
        
		for (int i = 0; i < nbParts; i++)
			if (nbPartItem[i] > 0)
                treePtr->son[i] = new TDTreeNode();
			else 
				treePtr->son[i] = nullptr;

        if (nbIsoItem > 0) 
            treePtr->iso = new TDTreeNode();
    }
}

// Create TDTree cellPart/facePart
void TDTree_create_itemPart (int *itemPart, int *nodePart, int *i2n, int localNbItem, int dimItem, int *nbPartItem, int *nbIsoItem)
{
	for (int i = 0; i < localNbItem; i++) {
		int node, colorA = -1, colorB = -1;
		for (int j = 0; j < dimItem; j++) {
			node = i2n[i*dimItem+j];
			if (colorA == -1) 
                colorA = nodePart[node];
            else
                if (colorA != nodePart[node]) 
                    colorB = nodePart[node];
		}
		if (colorB == -1) {
			nbPartItem[colorA]++;
			itemPart[i] = colorA;
		}
		else {
			(*nbIsoItem) ++;
			itemPart[i] = 1e9;
		}
	}
}

// Create TDTree itemPart according to cellPart/facePart
void TDTree_create_itemPart_accordingly(int *itemPerm, int *ItemPart, int *itemPart, int *I2i, int dimItem, 
                                int firstItem, int lastItem, int firstitem, int lastitem, int *nbPartitem, 
                                int *nbIsoitem)
{
    int localNbItem = lastItem - firstItem + 1;
    int localNbitem = lastitem - firstitem + 1;

    for(int i = 0; i < localNbitem; i++)
        itemPart[i] = 1e9;

    for(int i = 0; i < localNbItem; i++)
    {
        if (ItemPart[i] == 1e9) continue;
        for (int j = 0; j < dimItem; j++)
        {
            int item = itemPerm[I2i[(i+firstItem)*dimItem+j]];
            if (item < firstitem || item > lastitem) continue;
            item -= firstitem;
            itemPart[item] = ItemPart[i];
        }
    }

    for(int i = 0; i < localNbItem; i++)
    {
        if (ItemPart[i] != 1e9) continue;
        for (int j = 0; j < dimItem; j++)
        {
            int item = itemPerm[I2i[(i+firstItem)*dimItem+j]];
            if (item < firstitem || item > lastitem) continue;
            item -= firstitem;
            itemPart[item] = 1e9;
        }
    }
    
    for(int i = 0; i < localNbitem; i++)
        if (itemPart[i] == 1e9)
            (*nbIsoitem) ++;
        else
            nbPartitem[itemPart[i]] ++;
}

// Create TDTree for boundary faces
void TDTree::TDTree_create_boundary  (TDTreeNode *tree, int *f2n, int *f2c, int globalNbFace, int dimFace, 
                                                int firstFace, int lastFace, int globalNbNode,
                                                int firstNode, int lastNode, bool isIso)
{
    int localNbFace = lastFace - firstFace + 1;
    int localNbNode = lastNode - firstNode + 1;
  

	if (localNbFace < PARTSIZE) {
        TDTree_node_init (tree, 0, -1, firstFace, lastFace, firstNode, lastNode, 0, 0, isIso, true, 0);
        return;
    }


    int *nodePart = new int[localNbFace*dimFace];
    int *local_f2n = new int[localNbFace*dimFace];
    int parts = TDTree_partitioning (f2n, local_f2n, dimFace, firstFace, lastFace, nodePart);


    int nbIsoFace = 0;
    int *nbPartFace = new int [parts]();

    int *facePart = new int [localNbFace];
    TDTree_create_itemPart (facePart, nodePart, local_f2n, localNbFace, dimFace, nbPartFace, &nbIsoFace);
    delete[] local_f2n;
    delete[] nodePart;

    //Correct some faces that share the same cell
    //Remain one of them unmoved, others become isolator.
    unordered_map<int,int> mp;
    for(int i = 0; i < localNbFace; i++)
    {
        if (facePart[i] == 1e9)
            continue;
        int face = i + firstFace;
        if (!mp.count(f2c[face+face])||mp[f2c[face+face]]==facePart[i])
            mp[f2c[face+face]]=facePart[i];
        else
        {
            nbPartFace[facePart[i]]--;
            facePart[i] = 1e9;
            nbIsoFace++;
        }
    }
    mp.clear();

    int nbIsoNode = 0;
    int *nbPartNode = new int [parts]();
    
    if (localNbNode > 0)
    {
        nodePart = new int[localNbNode];
        TDTree_create_itemPart_accordingly (nodePerm, facePart, nodePart, f2n, dimFace, firstFace, lastFace, firstNode, lastNode, nbPartNode, &nbIsoNode);
        
        int *localNodePerm = new int [localNbNode];
        TDTree_create_permutation (localNodePerm, nodePart, localNbNode);
        delete[] nodePart;

        merge_permutations (nodePerm, localNodePerm, globalNbNode, localNbNode, firstNode, lastNode);
        delete[] localNodePerm;
    }

    int *localFacePerm = new int [localNbFace];
    TDTree_create_permutation (localFacePerm, facePart, localNbFace);
    delete[] facePart;


    merge_permutations (facePerm, localFacePerm, globalNbFace, localNbFace, firstFace, lastFace);
   

    TDTree_permute (f2n, localFacePerm, localNbFace, dimFace, firstFace);
    TDTree_permute (f2c, localFacePerm, localNbFace, 2, firstFace);
    delete[] localFacePerm;


    int nbIsoItem = max(nbIsoFace, nbIsoNode);
    int *nbPartItem = new int[parts]();
    for (int i = 0; i < parts; i++)
        nbPartItem[i] = max(nbPartFace[i], nbPartNode[i]);

    TDTree_node_init (tree, 0, -1, firstFace, lastFace, firstNode, lastNode, nbIsoItem, nbPartItem, isIso, false, parts);

    tree->isLeaf = true;
    int stFace[parts], edFace[parts], st2 = firstFace, ed2;
    int stNode[parts], edNode[parts], st3 = firstNode, ed3;
    for(int i = 0; i < parts; i++, st2 = ed2 + 1, st3 = ed3 + 1)
    {
        stFace[i] = st2;
        ed2 = edFace[i] = st2 + nbPartFace[i]-1;
        stNode[i] = st3;
        ed3 = edNode[i] = st3 + nbPartNode[i]-1;
    }

#ifndef FORKJOIN
#ifdef OMP
    #pragma omp taskloop default(shared)
    for (int i = 0; i < parts; i++) {
#elif CILK
    cilk_for (int i = 0; i < parts; i++) {
#endif
#else
    #pragma omp parallel for
    for (int i = 0; i < parts; i++) {
#endif
        if (nbPartItem[i] > 0)
        {
            tree->isLeaf = false;
#ifndef FORKJOIN
            TDTree_create_boundary (tree->son[i], f2n, f2c, globalNbFace, dimFace, stFace[i], edFace[i], globalNbNode, stNode[i], edNode[i], isIso);
#else
            TDTree_node_init (tree->son[i], 0, -1, stFace[i], edFace[i], stNode[i], edNode[i], 0, 0, isIso, true, parts);
#endif
        }
    }

#ifndef FORKJOIN
    // Synchronization
#ifdef OMP
        #pragma omp taskwait
#elif CILK
        cilk_sync;
#endif        
#endif

    delete[] nbPartFace;
    delete[] nbPartNode;
    delete[] nbPartItem;
    

    if (nbIsoItem > 0 && !tree->isLeaf) 
		TDTree_create_boundary (tree->iso, f2n, f2c, globalNbFace, dimFace, lastFace-nbIsoFace+1, lastFace, 
                                        globalNbNode, lastNode-nbIsoNode+1, lastNode, true);

}           

//Create TDTree for non-boundary cell
void TDTree::TDTree_create_normal  (TDTreeNode *tree, int *c2c, int *c2f, int *c2n, int globalNbCell, 
                                        int dimCell1, int dimCell2, int dimCell3, int firstCell, int lastCell, int globalNbFace, 
                                        int firstFace, int lastFace, int globalNbNode, int firstNode, int lastNode, bool isIso)
{
    int localNbCell = lastCell - firstCell + 1;
    int localNbFace = lastFace - firstFace + 1;
    int localNbNode = lastNode - firstNode + 1;

	if (localNbCell < PARTSIZE) {
        TDTree_node_init (tree, firstCell, lastCell, firstFace, lastFace, firstNode, lastNode, 0, 0, isIso, true, 0);
        return;
    }

    int *tmpCellPart = new int[localNbCell*dimCell3];
    int *local_c2n = new int[localNbCell*dimCell3];
    int parts = TDTree_partitioning (c2n, local_c2n, dimCell3, firstCell, lastCell, tmpCellPart);

    int nbIsoCell = 0;
    int *nbPartCell = new int [parts]();

    int *cellPart = new int [localNbCell];
    TDTree_create_itemPart (cellPart, tmpCellPart, local_c2n, localNbCell, dimCell3, nbPartCell, &nbIsoCell);
    delete[] local_c2n;
    delete[] tmpCellPart;

    unordered_map<int,int> mp;
    for(int i=firstCell;i<=lastCell;i++)
    {
        if (cellPart[i-firstCell]==1e9) continue;
        for(int j=0;j<dimCell1;j++)
            if (mp.count(c2c[i*dimCell1+j])&&mp[c2c[i*dimCell1+j]]!=cellPart[i-firstCell])
            {
                nbPartCell[cellPart[i-firstCell]]--;
                nbIsoCell++;
                cellPart[i-firstCell]=1e9;
                break;
            }
        if (cellPart[i-firstCell]==1e9) continue;
        for(int j=0;j<dimCell1;j++)
            mp[c2c[i*dimCell1+j]]=cellPart[i-firstCell];
    }

    //Correct some cells that share the same node
    //Remain one of them unmoved, others become isolator.
    int nbIsoFace = 0;
    int *nbPartFace = new int [parts]();
    if (localNbFace > 0)
    {
        int *facePart = new int[localNbFace];
        TDTree_create_itemPart_accordingly (facePerm, cellPart, facePart, c2f, dimCell2, firstCell, lastCell, firstFace, lastFace, nbPartFace, &nbIsoFace);
        
        int *localFacePerm = new int [localNbFace];
        TDTree_create_permutation (localFacePerm, facePart, localNbFace);
        delete[] facePart;

        merge_permutations (facePerm, localFacePerm, globalNbFace, localNbFace, firstFace, lastFace);
        delete[] localFacePerm;
    }

    int nbIsoNode = 0;
    int *nbPartNode = new int [parts]();
    if (localNbNode > 0)
    {
        int *nodePart = new int[localNbNode];
        TDTree_create_itemPart_accordingly (nodePerm, cellPart, nodePart, c2n, dimCell3, firstCell, lastCell, firstNode, lastNode, nbPartNode, &nbIsoNode);
        
        int *localNodePerm = new int [localNbNode];
        TDTree_create_permutation (localNodePerm, nodePart, localNbNode);
        delete[] nodePart;

        merge_permutations (nodePerm, localNodePerm, globalNbNode, localNbNode, firstNode, lastNode);
        delete[] localNodePerm;       
    }

    // Create local element permutation
    int *localCellPerm = new int [localNbCell];
    TDTree_create_permutation (localCellPerm, cellPart, localNbCell);
    delete[] cellPart;


    merge_permutations (cellPerm, localCellPerm, globalNbCell, localNbCell, firstCell, lastCell);
   
    // cout << "tree_creation 4 " << endl << flush;

    TDTree_permute (c2c, localCellPerm, localNbCell, dimCell1, firstCell);
    TDTree_permute (c2f, localCellPerm, localNbCell, dimCell2, firstCell);
    TDTree_permute (c2n, localCellPerm, localNbCell, dimCell3, firstCell);
    delete[] localCellPerm;


    int nbIsoItem = max(nbIsoCell, max(nbIsoFace, nbIsoNode));
    int *nbPartItem = new int[parts]();
    for (int i = 0; i < parts; i++)
        nbPartItem[i] = max(nbPartCell[i], max(nbPartFace[i], nbPartNode[i]));

    TDTree_node_init (tree, firstCell, lastCell, firstFace, lastFace, firstNode, lastNode, nbIsoItem, nbPartItem, isIso, false, parts);

    tree->isLeaf = true;
    int stCell[parts], edCell[parts], st1 = firstCell, ed1;
    int stFace[parts], edFace[parts], st2 = firstFace, ed2;
    int stNode[parts], edNode[parts], st3 = firstNode, ed3;
    for(int i = 0; i < parts; i++, st1 = ed1+1, st2 = ed2 + 1, st3 = ed3 + 1)
    {
        stCell[i] = st1;
        ed1 = edCell[i] = st1 + nbPartCell[i]-1;
        stFace[i] = st2;
        ed2 = edFace[i] = st2 + nbPartFace[i]-1;
        stNode[i] = st3;
        ed3 = edNode[i] = st3 + nbPartNode[i]-1;
    }

#ifndef FORKJOIN
#ifdef OMP
    #pragma omp taskloop default(shared)
    for (int i = 0; i < parts; i++) {
#elif CILK
    cilk_for (int i = 0; i < parts; i++) {
#endif
#else
    #pragma omp parallel for
    for (int i = 0; i < parts; i++) {
#endif
        if (nbPartItem[i] > 0)
        {
            tree->isLeaf = false;
#ifndef FORKJOIN
            TDTree_create_normal (tree->son[i], c2c, c2f, c2n, globalNbCell, dimCell1, dimCell2, dimCell3, stCell[i], edCell[i],
                                globalNbFace, stFace[i], edFace[i], globalNbNode, stNode[i], edNode[i], isIso);
#else
            TDTree_node_init (tree->son[i], stCell[i], edCell[i], stFace[i], edFace[i], stNode[i], edNode[i], 0, 0, isIso, true, parts);
#endif
        }
    }
 
#ifndef FORKJOIN
    // Synchronization
    #ifdef OMP
        #pragma omp taskwait
    #elif CILK
        cilk_sync;
    #endif        
#endif

    delete[] nbPartCell;
    delete[] nbPartFace;
    delete[] nbPartNode;
    delete[] nbPartItem;
    
    // Partitioning of Isoarator elements
    if (nbIsoItem > 0 && !tree->isLeaf)
		TDTree_create_normal (tree->iso, c2c, c2f, c2n, globalNbCell, dimCell1, dimCell2, dimCell3, lastCell-nbIsoCell+1, lastCell, 
                                                                        globalNbFace, lastFace-nbIsoFace+1, lastFace, 
                                                                        globalNbNode, lastNode-nbIsoNode+1, lastNode, 
                                                                        true);
    else
        tree->nbParts = 0;
}

// Create the TDTree
//dimCell1: dimension of c2c; dimCell2: dimension of c2f; dimCell3: dimension of c2n; dimFace: dimension of f2n;
//nPFace: physical face; nBFace: boundary face and physical face;
void TDTree::TDTree_creation (int *c2c, int *c2f, int *c2n, int *f2n, int *f2c, int globalNbCell, int dimCell1, int dimCell2, int dimCell3, 
                                    int globalNbFace, int dimFace, int globalNbNode, int nPFace, int nBFace)
{
    double start = omp_get_wtime();

#ifdef OMP
    #pragma omp parallel for
    for (int i = 0; i < globalNbCell; i++) 
        cellPerm[i] = i;

    #pragma omp parallel for
    for (int i = 0; i < globalNbFace; i++) 
        facePerm[i] = i;

    #pragma omp parallel for
    for (int i = 0; i < globalNbNode; i++) 
        nodePerm[i] = i;

#elif CILK

    cilk_for (int i = 0; i < globalNbCell; i++) 
        cellPerm[i] = i;
        
    cilk_for (int i = 0; i < globalNbFace; i++) 
        facePerm[i] = i;
        
    cilk_for (int i = 0; i < globalNbNode; i++) 
        nodePerm[i] = i;

#endif

    int nbIsoFace = globalNbFace - nPFace;
    int *nbPartFace = new int [1]();
    int *facePart = new int [globalNbFace]();
    
    //first level: physical face as submesh, others as isolator
    nbPartFace[0] = nPFace;
    TDTree_node_init (treeRoot, 0, globalNbCell-1, 0, globalNbFace-1, 0, globalNbNode-1, nbIsoFace, nbPartFace, false, false, 1);

    if (nPFace > 0)
        TDTree_create_boundary (treeRoot->son[0], f2n, f2c, globalNbFace, dimFace, 0, nPFace-1, globalNbNode, 0, -1, false);

    nbIsoFace = globalNbFace - nBFace;
    facePart = new int [globalNbFace-nPFace]();

    //second level: boundary face as submesh, other cell as isolator
    nbPartFace[0] = nBFace - nPFace;
    TDTree_node_init (treeRoot->iso, 0, globalNbCell-1, nPFace, globalNbFace-1, 0, globalNbNode-1, nbIsoFace, nbPartFace, false, false, 1);

    if (nBFace - nPFace > 0)
        TDTree_create_boundary (treeRoot->iso->son[0], f2n, f2c, globalNbFace, dimFace, nPFace, nBFace-1, globalNbNode, 0, -1, false);

    //now begin cell decomposition
#ifndef FORKJOIN
#ifdef OMP
        #pragma omp parallel
        #pragma omp single nowait
#endif
#endif
    TDTree_create_normal (treeRoot->iso->iso, c2c, c2f, c2n, globalNbCell, dimCell1, dimCell2, dimCell3, 0, globalNbCell-1, globalNbFace, nBFace, globalNbFace-1,
                                                            globalNbNode, 0, globalNbNode-1, false);

    //check for correctness

    // int *t = new int[globalNbCell]();
    // for (int i = 0; i < globalNbCell; i++) 
    //     t[cellPerm[i]]++;
    // for (int i = 0; i < globalNbCell; i++) 
    //     assert(t[i]==1);
    // delete[] t;

    // t = new int[globalNbFace]();
    // for (int i = 0; i < globalNbFace; i++) 
    //     t[facePerm[i]]++;
    // for (int i = 0; i < globalNbFace; i++) 
    //     assert(t[i]==1);
    // delete[] t;

    // t = new int[globalNbNode]();
    // for (int i = 0; i < globalNbNode; i++) 
    //     t[nodePerm[i]]++;
    // for (int i = 0; i < globalNbNode; i++) 
    //     assert(t[i]==1);
    // delete[] t;

    //get ordering
#ifdef OMP
        #pragma omp parallel for
        for (int i = 0; i < globalNbCell; i++) 
            cellRev[cellPerm[i]] = i;
            
        #pragma omp parallel for
        for (int i = 0; i < globalNbFace; i++) 
            faceRev[facePerm[i]] = i;
            
        #pragma omp parallel for
        for (int i = 0; i < globalNbNode; i++) 
            nodeRev[nodePerm[i]] = i;

#elif CILK
        cilk_for (int i = 0; i < globalNbCell; i++) 
            cellRev[cellPerm[i]] = i;
            
        cilk_for (int i = 0; i < globalNbFace; i++) 
            faceRev[facePerm[i]] = i;
            
        cilk_for (int i = 0; i < globalNbNode; i++) 
            nodeRev[nodePerm[i]] = i;
#endif

    double end = omp_get_wtime();
    cout << "Build Tree Time: " << end-start << endl << flush;
}

int* TDTree::TDTree_get_cellPerm()
{
    return cellPerm;
}

int* TDTree::TDTree_get_cellRev()
{
    return cellRev;
}
int* TDTree::TDTree_get_facePerm()
{
    return facePerm;
}

int* TDTree::TDTree_get_faceRev()
{
    return faceRev;
}

int* TDTree::TDTree_get_nodePerm()
{
    return nodePerm;
}

int* TDTree::TDTree_get_nodeRev()
{
    return nodeRev;
}
