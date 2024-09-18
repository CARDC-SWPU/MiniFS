#include <bits/stdc++.h>

#ifndef TDTREE_H
#define TDTREE_H

#define Forward 0
#define Backward 1
#define noDependence 2
#define CELL 3
#define FACE 4 

// TDTree structure
struct TDTreeNode 
{
	int nbParts;
	bool isIso, isLeaf;
	int firstCell, lastCell, firstFace, lastFace, firstNode, lastNode, vecOffset;
    struct TDTreeNode **son, *iso;
	int *boundryCell, *boundryNode;
	int boundryNbCell, boundryNbNode;	
	int *boundryCellNo, *boundryN;
};

// TDTree traversal arguments structure
struct TDTreeArgs 
{
    int firstCell, lastCell, firstFace, lastFace, firstNode, lastNode;
	int **boundryCell, **boundryNode;
	int *boundryNbCell, *boundryNbNode;	
	int **boundryCellNo, **boundryN;
};

typedef struct 
{
    int *index, *value;
} index_t;

using namespace std;

class TDTree
{
	private:
		TDTreeNode *treeRoot;
		int *cellPerm, *cellRev;
		int *facePerm, *faceRev;
		int *nodePerm, *nodeRev;
		int nbParts, PARTSIZE;
	
		int TDTree_partitioning (int *i2n, int *local_i2n, int dimItem, int firstItem, int lastItem, int *nodePart);

		void TDTree_create_boundary  (TDTreeNode *tree, int *f2n, int *f2c, int globalNbFace, int dimFace, 
                                                int firstFace, int lastFace, int globalNbNode,
                                                int firstNode, int lastNode, bool isIso);

		void TDTree_create_normal  (TDTreeNode *tree, int *c2c, int *c2f, int *c2n, int globalNbCell, 
                                        int dimCell1, int dimCell2, int dimCell3, int firstCell, int lastCell, int globalNbFace, 
                                        int firstFace, int lastFace, int globalNbNode, int firstNode, int lastNode, 
                                        bool isIso);
		
		void TDTree_permute (int *tab, int *perm, int nbItem, int dimItem, int offset);
									  
		void merge_permutations (int *perm, int *localPerm, int globalNbItem, int localNbItem, int firstItem, int lastItem);

		void TDTree_face_coloring (TDTreeNode treePtr, index_t &n2f, int *f2n, int globalNbFace, int dimFace);

		void TDTree_cell_coloring (TDTreeNode treePtr, index_t &n2c, int *c2n, int globalNbCell, int dimCell);
	public:
	
		void TDTree_creation (int *c2c, int *c2f, int *c2n, int *f2n, int *f2c, int globalNbCell, int dimCell1, int dimCell2, 
                                    int dimCell3, int globalNbFace, int dimFace, int globalNbNode, int nPFace, int nBFace);

		void task_traversal (void (*userSeqFctPtr)  (char **, TDTreeArgs *), void (*userVecFctPtr)  (char **, TDTreeArgs *), char **userArgs, int traversal_type);

		int* TDTree_get_cellPerm();

		int* TDTree_get_cellRev();

		int* TDTree_get_facePerm();

		int* TDTree_get_faceRev();
		
		int* TDTree_get_nodePerm();

		int* TDTree_get_nodeRev();
							
		TDTree(int globalNbCell, int globalNbFace, int globalNbNode, int nbparts, int partSize) 
		{
			treeRoot = new TDTreeNode();
			
			cellPerm = new int [globalNbCell];
			cellRev = new int [globalNbCell];
			facePerm = new int[globalNbFace];
			faceRev = new int [globalNbFace];
			nodePerm = new int [globalNbNode];
			nodeRev = new int [globalNbNode];

			nbParts = nbparts;
			PARTSIZE = partSize;
		}

	
};

void TDTree_create_permutation (int *perm, int *part, int size);
#endif