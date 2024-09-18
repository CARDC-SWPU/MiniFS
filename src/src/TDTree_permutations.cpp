#include "TDTree_permutations.h"
#include <bits/stdc++.h>
#include "TDTree.h"

using namespace std;

void TDTree::TDTree_permute(int *tab, int *perm, int nbItem, int dimItem, int offset)
{
    bool *vis = new bool [nbItem] ();
    int  *tmpSrc    = new int  [dimItem];
    int  *tmpDst    = new int  [dimItem];

    for (int i = 0; i < nbItem; i++) 
    {
        if (vis[i] == 1) continue;

        int init = i, src = i, dst;
        for (int j = 0; j < dimItem; j++) 
            tmpSrc[j] = tab[(i+offset)*dimItem+j];
        do 
        {
            dst = perm[src];
            for (int j = 0; j < dimItem; j++) 
            {
                tmpDst[j] = tab[(dst+offset)*dimItem+j];
                tab[(dst+offset)*dimItem+j] = tmpSrc[j];
                tmpSrc[j] = tmpDst[j];
            }
            src = dst;
            vis[src] = 1;
        }
        while (src != init);
    }
    delete[] tmpDst, delete[] tmpSrc, delete[] vis;
}

void TDTree::merge_permutations (int *perm, int *localItemPerm, int globalNbItem, int localNbItem, int firstItem, int lastItem)
{
    int ptr = 0;
    for (int i = 0; i < globalNbItem; i++) 
    {
        int dst = perm[i];
        if (dst >= firstItem && dst <= lastItem) 
        {
            perm[i] = localItemPerm[dst-firstItem] + firstItem;
            ptr++;
        }
        if (ptr == localNbItem)	break;
    }
}

// Create permutation array from partition array
void TDTree_create_permutation (int *perm, int *part, int size)
{
    pair<int,int> *a;
    a = new pair<int,int>[size];
    
    for (int i = 0; i < size; i++)
        a[i] = make_pair(part[i], i);
    
    sort(a, a + size);
    
    for (int i = 0; i < size; i++)
        perm[a[i].second] = i;
    
    delete[] a;
}
