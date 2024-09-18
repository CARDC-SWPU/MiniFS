#ifndef TDTree_creation_h
#define TDTree_creation_h

#include "TDTree.h"

// Initialize the content of D&C tree nodes
void TDTree_node_init (TDTreeNode *treePtr, int firstCell, int lastCell, int firstFace, int lastFace, int firstNode, int lastNode, 
                            int nbIsoItem, int *nbPartItem, bool isIso, bool isLeaf, int nbParts);
				   
void TDTree_create_itemPart (int *itemPart, int *nodePart, int *i2n, int localNbItem, int dimItem, int *nbPartItem, int *nbIsoItem);

void TDTree_create_itemPart_accordingly(int *itemPerm, int *ItemPart, int *itemPart, int *I2i, int dimItem, 
                                int firstItem, int lastItem, int firstitem, int lastitem, int *nbPartitem, 
                                int *nbIsoitem, int *i2I);

#endif
