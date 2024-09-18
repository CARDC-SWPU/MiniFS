#include "TDTree.h"

void TDTree_traversal_forward (void (*userSeqFctPtr)  (char **, TDTreeArgs *),
                     void (*userVecFctPtr)  (char **, TDTreeArgs *),
                     char **userArgs, TDTreeNode *treePtr, int nbParts);

void TDTree_traversal_backward (void (*userSeqFctPtr)  (char **, TDTreeArgs *),
                     void (*userVecFctPtr)  (char **, TDTreeArgs *),
                     char **userArgs, TDTreeNode *tree, int nbParts);

void TDTree_traversal_noDependence (void (*userSeqFctPtr)  (char **, TDTreeArgs *),
                     void (*userVecFctPtr)  (char **, TDTreeArgs *),
                     char **userArgs, TDTreeNode *treePtr, int nbParts);
