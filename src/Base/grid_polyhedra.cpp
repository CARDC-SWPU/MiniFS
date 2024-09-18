/*!
 * @file        grid_polyhedra.cpp
 * @brief       A class for unstructured polyhedral grid
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @version     1.0
 * @date        2023-05-21
 * @copyright   Copyright (c) 2023  PERSONAL
 *
 * @par 修改日志:
 * <table>
 * <tr><th> Date        <th> Version  <th> Author  <th> Description
 * <tr><td> 2023-05-21  <td> 1.0      <td> Wisces  <td> 〔内容〕
 * </table>
 */

// direct head file
#include "grid_polyhedra.h"

// build-in head files
// #include <cstring>
// #include <fstream>
#include <iostream>
// #include <sstream>
// #include <string>
// #include <cstdlib>
// #include <cmath>
// #include <cassert>
// #include <deque>
// #include <queue>
// #include <list>
// #include <set>
// #include <map>
// #include <iomanip>
using namespace std;

// user defined head files
#include "number_type.h"
#include "memory_util.h"
#include "para_field_global.h"
#include "constant.h"
// #include "zone.h"
// #include "solver_ns.h"
// #include "utility_functions.h"
// #include "algm.h"
// #include "io_log.h"
// #include "io_base_format.h"（直接 copy 到下面）
// define separation of a blank between two words

// #include "parallel_base_functions.h"
// #include "system_base_functions.h"
// #include "grid_patch_type.h"（集成到 constant.h）

// this header file is copied from cart

// #ifdef MF_MPICH
// #include "mpi.h"
// #endif
// #ifdef MF_OPENMP
// #include "grid_decoupling.h"
// #endif

//!< head file relying on condition-compiling
#ifdef MF_MPICH
#include <mpi.h>
#endif

#ifdef TDTREE
#include "TDTree.h"
#endif

#ifdef MF_MPICH
extern int myid;
extern int numprocs;
extern int myZone;        // = myid+1
extern MPI_Comm GridComm; ///< for each grid
#endif

#ifdef TDTREE
/*!
 * @brief
 * @param       grid
 * @param       fgrid
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-12
 */
void TDTree_CellReordering(PolyGrid *grid, PolyGrid *fgrid)
{
    TDTree *TDTreeRoot = grid->GetTDTree();
    IntType *index = TDTreeRoot->TDTree_get_cellPerm();

    IntType i;
    IntType nTCell = grid->GetNTCell();
    IntType nTFace = grid->GetNTFace();
    IntType *f2c = grid->Getf2c();

    IntType *t = new int[nTCell]();
    for (i = 0; i < nTCell; i++)
        t[index[i]]++;

    // for (i = 0; i < nTCell; i++)
    //     assert(t[i] == 1);

    for (i = 0; i < nTFace; i++)
    {
        if (f2c[i << 1] >= 0 && f2c[i << 1] < nTCell)
            f2c[i << 1] = index[f2c[i << 1]];
        if (f2c[i << 1 | 1] >= 0 && f2c[i << 1 | 1] < nTCell)
            f2c[i << 1 | 1] = index[f2c[i << 1 | 1]];
    }

    if (fgrid != NULL)
    {
        IntType *c2cc = fgrid->Getc2cc();
        IntType *c2cc_backup = NULL;
        IntType fnTCell = fgrid->GetNTCell();
        snew_array_1D(c2cc_backup, fnTCell);
        for (IntType i = 0; i < fnTCell; i++)
            c2cc_backup[i] = c2cc[i];
        for (IntType i = 0; i < fnTCell; i++)
            c2cc[i] = index[c2cc_backup[i]];
        sdel_array_1D(c2cc_backup);
    }

    cout << "TDTree Cell Reordering successfully!" << endl;
}

/*!
 * @brief       use DC Tree index to do face reorder
 * @param       grid
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-12
 */
void TDTree_FaceReordering(PolyGrid *grid)
{
    IntType nBFace = grid->GetNBFace();
    IntType nTFace = grid->GetNTFace();
    IntType nIFace = grid->GetNIFace();
    IntType pfacenum = nBFace - nIFace;
    IntType *f2c = grid->Getf2c();
    IntType *f2n = grid->Getf2n();
    IntType *nNPF = grid->GetnNPF();
    IntType *nbz = grid->GetnbZ();
    IntType *nbf = grid->GetnbBF();

    TDTree *TDTreeRoot = grid->GetTDTree();
    IntType *index = TDTreeRoot->TDTree_get_faceRev();
    IntType *Index = TDTreeRoot->TDTree_get_facePerm();

    // IntType level = grid->GetLevel();
    // if (level != 0)
    // {
    //     BCRecord **bcr = grid->Getbcr();
    //     BCRecord **bcr_backup = NULL;
    //     snew_array_1D(bcr_backup, nBFace);
    //     for (IntType i = 0; i < nBFace; i++)
    //     {
    //         bcr_backup[i] = bcr[i];
    //     }
    //     for (IntType i = 0; i < nBFace; i++)
    //     {
    //         assert(index[i] < nBFace);
    //         bcr[i] = bcr_backup[index[i]];
    //     }
    //     sdel_array_1D(bcr_backup);
    // }

    IntType *f2c_backup = NULL;
    IntType *f2n_backup = NULL;
    IntType *nNPF_backup = NULL;

    snew_array_1D(f2c_backup, nTFace * 2);
    snew_array_1D(nNPF_backup, nTFace);

    IntType n = 0;
    for (IntType i = 0; i < nTFace; ++i)
    {
        f2c_backup[i + i] = f2c[i + i];
        f2c_backup[i + i + 1] = f2c[i + i + 1];
        nNPF_backup[i] = nNPF[i];
        n += nNPF_backup[i];
    }

    snew_array_1D(f2n_backup, n);
    for (IntType i = 0; i < n; ++i)
        f2n_backup[i] = f2n[i];

    IntType **F2N_backup = NULL;
    snew_array_1D(F2N_backup, nTFace);
    F2N_backup[0] = f2n_backup;
    for (IntType i = 1; i < nTFace; ++i)
    {
        F2N_backup[i] = &(F2N_backup[i - 1][nNPF_backup[i - 1]]);
    }

    n = 0;

    for (IntType i = 0; i < nTFace; ++i)
    {
        f2c[i + i] = f2c_backup[index[i] + index[i]];
        f2c[i + i + 1] = f2c_backup[index[i] + index[i] + 1];
        nNPF[i] = nNPF_backup[index[i]];
        for (IntType j = 0; j < nNPF[i]; ++j)
            f2n[n++] = F2N_backup[index[i]][j];
    }

#ifdef MF_MPICH

    IntType *nbz_backup = NULL;
    IntType *nbf_backup = NULL;
    snew_array_1D(nbz_backup, nIFace);
    snew_array_1D(nbf_backup, nIFace);

    for (IntType i = 0; i < nIFace; ++i)
    {
        nbz_backup[i] = nbz[i];
        nbf_backup[i] = nbf[i];
    }

    IntType size, id;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    IntType *nbFace = new int[size];
    IntType **facePerm = new int *[size];
    for (IntType i = 0; i < size; i++)
    {
        if (i == id)
        {
            for (IntType j = 0; j < size; j++)
                if (i != j)
                {
                    MPI_Send(&nBFace, 1, MPI_INT, j, j, MPI_COMM_WORLD);
                    MPI_Send(Index, nBFace, MPI_INT, j, j + size, MPI_COMM_WORLD);
                }
        }
        else
        {
            MPI_Recv(&nbFace[i], 1, MPI_INT, i, id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            facePerm[i] = new int[nbFace[i]];
            MPI_Recv(facePerm[i], nbFace[i], MPI_INT, i, id + size, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    for (IntType i = 0; i < nIFace; ++i)
    {
        nbz[i] = nbz_backup[index[i + pfacenum] - pfacenum];
        nbf[i] = facePerm[nbz[i]][nbf_backup[index[i + pfacenum] - pfacenum]];
    }

    for (IntType i = 0; i < size; i++)
        if (i != id)
            delete[] facePerm[i];
    delete[] facePerm;
    delete[] nbFace;

    sdel_array_1D(nbz_backup);
    sdel_array_1D(nbf_backup);

#endif

    cout << "TDTree Face Reordering successfully!" << endl;

    sdel_array_1D(f2c_backup);
    sdel_array_1D(f2n_backup);
    sdel_array_1D(nNPF_backup);
    sdel_array_1D(F2N_backup);
}

/*!
 * @brief
 * @param       grid
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-12
 */
void TDTree_NodeReordering(PolyGrid *grid)
{
    IntType nTFace = grid->GetNTFace();
    IntType nTNode = grid->GetNTNode();
    IntType nINode = grid->GetNINode();
    IntType *f2n = grid->Getf2n();
    IntType *nNPF = grid->GetnNPF();
    RealGeom *x = grid->GetX();
    RealGeom *y = grid->GetY();
    RealGeom *z = grid->GetZ();

    TDTree *TDTreeRoot = grid->GetTDTree();
    IntType *index = TDTreeRoot->TDTree_get_nodePerm();
    IntType *Index = TDTreeRoot->TDTree_get_nodeRev();

    IntType n = 0;
    for (IntType i = 0; i < nTFace; ++i)
        n += nNPF[i];

    for (IntType i = 0; i < n; ++i)
        f2n[i] = index[f2n[i]];

    n = 0;

    RealGeom *x_backup = NULL;
    RealGeom *y_backup = NULL;
    RealGeom *z_backup = NULL;

    snew_array_1D(x_backup, nTNode);
    snew_array_1D(y_backup, nTNode);
    snew_array_1D(z_backup, nTNode);

    for (IntType i = 0; i < nTNode; i++)
    {
        x_backup[i] = x[i];
        y_backup[i] = y[i];
        z_backup[i] = z[i];
    }

    for (IntType i = 0; i < nTNode; i++)
    {
        x[i] = x_backup[Index[i]];
        y[i] = y_backup[Index[i]];
        z[i] = z_backup[Index[i]];
    }

#ifdef MF_MPICH

    IntType *nbsN = grid->GetnbSN();
    IntType *nbzN = grid->GetnbZN();
    IntType *nbrN = grid->GetnbRN();

    IntType size, id;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    IntType *ntNode = new int[size];
    IntType **nodePerm = new int *[size];
    for (IntType i = 0; i < size; i++)
    {
        if (i == id)
        {
            for (IntType j = 0; j < size; j++)
                if (i != j)
                {
                    MPI_Send(&nTNode, 1, MPI_INT, j, j, MPI_COMM_WORLD);
                    MPI_Send(index, nTNode, MPI_INT, j, j + size, MPI_COMM_WORLD);
                }
        }
        else
        {
            MPI_Recv(&ntNode[i], 1, MPI_INT, i, id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            nodePerm[i] = new int[ntNode[i]];
            MPI_Recv(nodePerm[i], ntNode[i], MPI_INT, i, id + size, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    for (IntType i = 0; i < nINode; i++)
    {
        nbsN[i] = index[nbsN[i]];
        nbrN[i] = nodePerm[nbzN[i]][nbrN[i]];
    }

#endif

    cout << "TDTree Node Reordering successfully!" << endl;

    sdel_array_1D(x_backup);
    sdel_array_1D(y_backup);
    sdel_array_1D(z_backup);
}

/*!
 * @brief       create DC tree of face for DC task parallelism
 * @param       grid
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-12
 */
void CreateTree(PolyGrid *grid)
{
    int nbParts;
#ifndef FORKJOIN
    nbParts = 4;
#else
    nbParts = 24 * 8;
    // cout << "nbParts = " << nbParts << endl;
#endif

    IntType i, j, n = 0;
    IntType dimCell1 = 0, dimCell2 = 0;
    IntType dimCell3 = 0, dimFace = 0;
    IntType nTFace = grid->GetNTFace();
    IntType nTCell = grid->GetNTCell();
    IntType nTNode = grid->GetNTNode();

    IntType nBFace = grid->GetNBFace();
    IntType nIFace = grid->GetNIFace();
    IntType pfacenum = nBFace - nIFace;

    IntType *f2c = grid->Getf2c();
    IntType *f2n = grid->Getf2n();
    IntType *nNPF = grid->GetnNPF();

    IntType *f2n_backup = NULL;
    IntType *f2c_backup = NULL;

    // IntType level = grid->GetLevel();
    int partSize = 512;
    // for (int i = 0; i < level; i++)
    //     partSize /= 2;

    for (i = 0; i < nTFace; ++i)
        dimFace = std::max(dimFace, nNPF[i]);

    snew_array_1D(f2n_backup, nTFace * dimFace);
    snew_array_1D(f2c_backup, nTFace * 2);

    for (i = 0; i < nTFace; i++)
    {
        for (j = 0; j < nNPF[i]; j++)
            f2n_backup[i * dimFace + j] = f2n[n++];
        for (; j < dimFace; j++)
            f2n_backup[i * dimFace + j] = f2n_backup[i * dimFace];
        f2c_backup[i << 1] = f2c[i << 1];
        f2c_backup[i << 1 | 1] = f2c[i << 1 | 1];
    }

    IntType *c2c = NULL;
    IntType *c2f = NULL;
    IntType *c2n = NULL;

    unordered_map<int, bool> *C2C = new unordered_map<int, bool>[nTCell];
    unordered_map<int, bool> *C2F = new unordered_map<int, bool>[nTCell];
    unordered_map<int, bool> *C2N = new unordered_map<int, bool>[nTCell];

    n = 0;
    for (i = 0; i < nTFace; i++)
    {
        for (j = 0; j < nNPF[i]; j++, n++)
        {
            C2N[f2c[i + i]][f2n[n]] = 1;
            dimCell3 = std::max(dimCell3, (int)C2N[f2c[i + i]].size());
            if (i >= nBFace)
            {
                C2N[f2c[i + i + 1]][f2n[n]] = 1;
                dimCell3 = std::max(dimCell3, (int)C2N[f2c[i + i + 1]].size());
            }
        }

        C2F[f2c[i + i]][i] = 1;
        dimCell2 = std::max(dimCell2, (int)C2F[f2c[i + i]].size());
        if (i >= nBFace)
        {
            C2F[f2c[i + i + 1]][i] = 1;
            dimCell2 = std::max(dimCell2, (int)C2F[f2c[i + i + 1]].size());

            C2C[f2c[i + i]][f2c[i + i]] = 1;
            C2C[f2c[i + i + 1]][f2c[i + i + 1]] = 1;
            C2C[f2c[i + i]][f2c[i + i + 1]] = 1;
            C2C[f2c[i + i + 1]][f2c[i + i]] = 1;
            dimCell1 = std::max(dimCell1, (int)std::max(C2C[f2c[i + i]].size(), C2C[f2c[i + i + 1]].size()));
        }
    }

    c2n = new int[dimCell3 * nTCell];
    c2f = new int[dimCell2 * nTCell];
    c2c = new int[dimCell1 * nTCell];

    // printGrid(grid, C2F);

    for (i = 0; i < nTCell; i++)
    {
        j = 0;
        for (auto cell : C2N[i])
        {
            c2n[i * dimCell3 + j] = cell.first;
            j++;
        }
        for (; j < dimCell3; j++)
            c2n[i * dimCell3 + j] = c2n[i * dimCell3];
        C2N[i].clear();

        j = 0;
        for (auto cell : C2F[i])
        {
            c2f[i * dimCell2 + j] = cell.first;
            j++;
        }
        for (; j < dimCell2; j++)
            c2f[i * dimCell2 + j] = c2f[i * dimCell2];
        C2F[i].clear();

        j = 0;
        for (auto cell : C2C[i])
        {
            c2c[i * dimCell1 + j] = cell.first;
            j++;
        }
        for (; j < dimCell1; j++)
            c2c[i * dimCell1 + j] = c2c[i * dimCell1];
        C2C[i].clear();
    }

    delete[] C2N;
    delete[] C2F;
    delete[] C2C;

    cout << endl
         << "TDTree begin to create!\n"
         << flush;
    TDTree *TDTreeRoot = new TDTree(nTCell, nTFace, nTNode, nbParts, partSize);
    TDTreeRoot->TDTree_creation(c2c, c2f, c2n, f2n_backup, f2c_backup, nTCell, dimCell1, dimCell2, dimCell3, nTFace, dimFace, nTNode, pfacenum, nBFace);
    grid->SetTDTree(TDTreeRoot);
    cout << "TDTree created successfully!\n"
         << flush;
    // #ifdef GPU
    //     IntType **stCell = new int *[50], **edCell = new int *[50];
    //     int groups = TDTreeRoot->TDTree_get_GPU_Groups(&stCell, &edCell);
    //     grid->SetStCell(stCell);
    //     grid->SetEdCell(edCell);
    //     grid->SetGroups(groups);
    // #endif
    sdel_array_1D(f2n_backup);
    sdel_array_1D(f2c_backup);
    delete[] c2n;
    delete[] c2f;
    delete[] c2c;
}
#endif

/**
PolyGrid::PolyGrid(IntType zin, IntType lin) : Grid(zin), level(lin), nIFace(0), nNeighbor(0), nINode(0), nNeighborN(0), VolAvg(0.0)
{
    // 指针
    cGrid = NULL;
    fGrid = NULL;
    c2cc = NULL;
    nCCPN = NULL;
    nVPN = NULL;
    weight_prol = NULL;
    nNPF = NULL;
    f2c = NULL;
    nCPC = NULL;
    nNPC = NULL;
    nFPC = NULL;
    nbZ = NULL;
    nbBF = NULL;
    nZIFace = NULL;
    nZINode = NULL;
    nbN = NULL;
    nbSN = NULL;
    nbZN = NULL;
    nbRN = NULL;
    WeightNodeDist = NULL;
    Nmark = NULL;
    xcc = NULL;
    ycc = NULL;
    zcc = NULL;
    xfn = NULL;
    yfn = NULL;
    zfn = NULL;
    vol = NULL;
    area = NULL;
    xfc = NULL;
    yfc = NULL;
    zfc = NULL;
    vgn = NULL;
    BFacevgx = NULL;
    BFacevgy = NULL;
    BFacevgz = NULL;
    vccx = NULL;
    vccy = NULL;
    vccz = NULL;
    node_act = NULL;
    edge_act = NULL;
    vgn = NULL;
    nChCellS = NULL;
    nChCellR = NULL;
    nChRNoComp = NULL;
    Cell2Zone = NULL;
    facecentroidskewness = NULL;
    faceangleskewness = NULL;
    cellcentroidskewness = NULL;
    cellvolsmoothness = NULL;
    faceangle = NULL;
    cellwallnumber = NULL;
    WeightNodeProl = NULL;
    nCPN = NULL;
    N2C = NULL;
    f2n = NULL;
    nb = NULL;

    // 指针的指针
    N2CC = NULL;
    N2V = NULL;
    c2c = NULL;
    C2F = NULL;
    F2N = NULL;
    C2N = NULL;
    bcr = NULL;
    bNSNo = NULL;
    bNRNo = NULL;
    Seta_Center = NULL;
    nChSNo = NULL;
    nChRNo = NULL;
    nChxcc = NULL;
    nChycc = NULL;
    nChzcc = NULL;
    WeightNodeC2N = NULL;
    WeightNodeN2C = NULL;
    bCNo = NULL;
    nbg = NULL;
    bFNo = NULL;
    // #ifdef USING_PETSC
    //     ghost2global = NULL;
    // #endif
    // data for line implicit
    nLinesForLI = 0;
    nCellsInLineForLI = 0;
    CellsOfLinesForLI = 0;
    FacesOfLinesForLI = 0;
    LineIndexOfCellsForLI = 0;
    // 物面出发的线原有数目和线长,(在连接之前),用于CFL数放大使用
    n_lines_wall_start = 0;
    n_lines_wall = 0;
    n_cells_in_line_wall = 0;
    span_len_of_2d = -1.0;
}

PolyGrid::PolyGrid(IntType zin) : Grid(zin), level(0), nIFace(0), nNeighbor(0), nINode(0), nNeighborN(0), VolAvg(0.0)
{
    // 指针
    cGrid = NULL;
    fGrid = NULL;
    c2cc = NULL;
    nCCPN = NULL;
    nVPN = NULL;
    weight_prol = NULL;
    nNPF = NULL;
    f2c = NULL;
    nCPC = NULL;
    nNPC = NULL;
    nFPC = NULL;
    nbZ = NULL;
    nbBF = NULL;
    nZIFace = NULL;
    nZINode = NULL;
    nbN = NULL;
    nbSN = NULL;
    nbZN = NULL;
    nbRN = NULL;
    WeightNodeDist = NULL;
    Nmark = NULL;
    xcc = NULL;
    ycc = NULL;
    zcc = NULL;
    xfn = NULL;
    yfn = NULL;
    zfn = NULL;
    vol = NULL;
    area = NULL;
    xfc = NULL;
    yfc = NULL;
    zfc = NULL;
    vgn = NULL;
    BFacevgx = NULL;
    BFacevgy = NULL;
    BFacevgz = NULL;
    vccx = NULL;
    vccy = NULL;
    vccz = NULL;
    node_act = NULL;
    edge_act = NULL;
    vgn = NULL;
    nChCellS = NULL;
    nChCellR = NULL;
    nChRNoComp = NULL;
    Cell2Zone = NULL;
    facecentroidskewness = NULL;
    faceangleskewness = NULL;
    cellcentroidskewness = NULL;
    cellvolsmoothness = NULL;
    faceangle = NULL;
    cellwallnumber = NULL;
    WeightNodeProl = NULL;
    nCPN = NULL;
    f2n = NULL;
    nb = NULL;

    // 指针的指针;
    N2CC = NULL;
    N2V = NULL;
    c2c = NULL;
    C2F = NULL;
    F2N = NULL;
    C2N = NULL;
    bcr = NULL;
    bNSNo = NULL;
    bNRNo = NULL;
    Seta_Center = NULL;
    nChSNo = NULL;
    nChRNo = NULL;
    nChxcc = NULL;
    nChycc = NULL;
    nChzcc = NULL;
    WeightNodeC2N = NULL;
    N2C = NULL;
    WeightNodeN2C = NULL;
    bCNo = NULL;
    nbg = NULL;
    bFNo = NULL;
    // #ifdef USING_PETSC
    //     ghost2global = NULL;
    // #endif
    // data for line implicit
    nLinesForLI = 0;
    nCellsInLineForLI = 0;
    CellsOfLinesForLI = 0;
    FacesOfLinesForLI = 0;
    LineIndexOfCellsForLI = 0;
    // 物面出发的线原有数目和线长,(在连接之前),用于CFL数放大使用
    n_lines_wall_start = 0;
    n_lines_wall = 0;
    n_cells_in_line_wall = 0;
    span_len_of_2d = -1.0;
}
*/

PolyGrid::PolyGrid() : Grid(), level(0), nTFace(0), nTCell(0), nBFace(0), nIFace(0), nINode(0), nNeighbor(0), nNeighborN(0) //, VolAvg(0.0)
{
#ifdef MF_DEBUG
    cout << "Constructor [PolyGrid()] is called!" << endl;
#endif

    //!< flow variables, dt_timestep, residuals, DQ - gField
    //!< add by Wisces
    //!< 网格重排序
    LUSGSLayer = NULL;
    LUSGSCellOrder = NULL;
    LUSGScellsPerlayer = NULL;

    //!< flow variables
    rho = NULL;
    u = NULL;
    v = NULL;
    w = NULL;
    p = NULL;

    //!< gradient
    dqdx = NULL;
    dqdy = NULL;
    dqdz = NULL;

    //!< time step
    dt = NULL; ///< dt_timestep

    //!< residuals
    res = NULL;

    //!< DQ
    DQ = NULL;

    //!< Some parameters for the Gradient computations
    //!< add by Wisces
    // IntType GaussLayer;
    // IntType *CellLayerNo;
    nodesymm = NULL;
    xfn_n_symm = NULL;
    yfn_n_symm = NULL;
    zfn_n_symm = NULL;

    // 指针
    // cGrid = NULL;
    // fGrid = NULL;
    c2cc = NULL;
    // nCCPN = NULL;
    // nVPN = NULL;
    // weight_prol = NULL;
    nNPF = NULL;
    f2n = NULL;
    F2N = NULL; ///< **
    f2c = NULL;

    xcc = NULL;
    ycc = NULL;
    zcc = NULL;

    xfc = NULL;
    yfc = NULL;
    zfc = NULL;

    xfn = NULL;
    yfn = NULL;
    zfn = NULL;

    vol = NULL;
    area = NULL;

    bcr = NULL; ///< **

    nCPC = NULL;
    c2c = NULL; ///< **
    nNPC = NULL;
    C2N = NULL; ///< **
    nFPC = NULL;
    C2F = NULL; ///< **
    nCPN = NULL;
    N2C = NULL; ///< **

    nb = NULL;
    nbN = NULL;

    nbZ = NULL;
    nbBF = NULL;
    nbSN = NULL;
    nbZN = NULL;
    nbRN = NULL;

    nZIFace = NULL;
    bCNo = NULL; ///< **
    bFNo = NULL; ///< **
    nZINode = NULL;
    bNSNo = NULL; ///< **
    bNRNo = NULL; ///< **

    Nmark = NULL;
    WeightNodeDist = NULL;
    WeightNodeC2N = NULL;     ///< **
    WeightNodeN2C = NULL;     ///< **
    WeightNodeBFace2C = NULL; ///< **

    vgn = NULL;
    BFacevgx = NULL;
    BFacevgy = NULL;
    BFacevgz = NULL;
    vccx = NULL;
    vccy = NULL;
    vccz = NULL;
    // node_act = NULL;
    // edge_act = NULL;
    // nChCellS = NULL;
    // nChCellR = NULL;
    // nChRNoComp = NULL;
    // Cell2Zone = NULL;
    // facecentroidskewness = NULL;
    // faceangleskewness = NULL;
    // cellcentroidskewness = NULL;
    // cellvolsmoothness = NULL;
    // faceangle = NULL;
    // cellwallnumber = NULL;
    // WeightNodeProl = NULL;

    // 指针的指针;
    // N2CC = NULL;
    // N2V = NULL;

    // Seta_Center = NULL;
    // nChSNo = NULL;
    // nChRNo = NULL;
    // nChxcc = NULL;
    // nChycc = NULL;
    // nChzcc = NULL;

    // nbg = NULL;

    // #ifdef USING_PETSC
    //     ghost2global = NULL;
    // #endif
    //     // data for line implicit
    //     nLinesForLI = 0;
    //     nCellsInLineForLI = 0;
    //     CellsOfLinesForLI = 0;
    //     FacesOfLinesForLI = 0;
    //     LineIndexOfCellsForLI = 0;
    //     // 物面出发的线原有数目和线长,(在连接之前),用于CFL数放大使用
    //     n_lines_wall_start = 0;
    //     n_lines_wall = 0;
    //     n_cells_in_line_wall = 0;
    //     span_len_of_2d = -1.0;
#ifdef OMP_DIVREP
    DivRepSuccess = false;
    idx_pthreads_bface = NULL;
    id_division_bface = NULL;
    idx_pthreads_iface = NULL;
    id_division_iface = NULL;
    // #ifdef BoundedColoring
    //         endIndex_bFace_vec = NULL;
    //         endIndex_iFace_vec = NULL;
    // #endif ///< ~BoundedColoring
#endif ///< ~OMP_DIVREP
}

PolyGrid::~PolyGrid()
{
    //!< flow variables, dt_timestep, residuals, DQ - gField
    //!< add by Wisces
    sdel_array_1D(LUSGSLayer);
    sdel_array_1D(LUSGSCellOrder);
    sdel_array_1D(LUSGScellsPerlayer);

    sdel_array_1D(rho);
    sdel_array_1D(u);
    sdel_array_1D(v);
    sdel_array_1D(w);
    sdel_array_1D(p);

    sdel_array_1D(dqdx);
    sdel_array_1D(dqdy);
    sdel_array_1D(dqdz);

    sdel_array_1D(dt);
    sdel_array_1D(res);
    sdel_array_1D(DQ);

    //!< Some parameters for the Gradient computations
    //!< add by Wisces
    // IntType GaussLayer;
    // IntType *CellLayerNo;
    sdel_array_1D(nodesymm);
    sdel_array_1D(xfn_n_symm);
    sdel_array_1D(yfn_n_symm);
    sdel_array_1D(zfn_n_symm);

    sdel_array_1D(c2cc);
    // sdel_array_1D(nCCPN);
    // sdel_array_1D(nVPN);
    // sdel_array_1D(weight_prol);

    sdel_array_1D(nNPF);
    sdel_array_1D(f2n);
    sdel_array_1D(F2N); // F2N is the reference to f2n, so use 1D array operator, tangj
    sdel_array_1D(f2c);

    sdel_array_1D(xcc);
    sdel_array_1D(ycc);
    sdel_array_1D(zcc);

    sdel_array_1D(xfc);
    sdel_array_1D(yfc);
    sdel_array_1D(zfc);

    sdel_array_1D(xfn);
    sdel_array_1D(yfn);
    sdel_array_1D(zfn);

    sdel_array_1D(vol);
    sdel_array_1D(area);

    sdel_array_1D(bcr);

    sdel_array_1D(nCPC);
    sdel_array_2D(c2c);
    sdel_array_1D(nNPC);
    sdel_array_2D(C2N);
    sdel_array_1D(nFPC);
    sdel_array_2D(C2F);
    sdel_array_1D(nCPN);
    IntType nTNode = this->GetNTNode();
    sdel_array_2D(N2C, nTNode, false);

    sdel_array_1D(nb);
    sdel_array_1D(nbN);
    // sdel_array_2D(N2CC);
    // sdel_array_2D(N2V);

    sdel_array_1D(nbZ);
    sdel_array_1D(nbBF);
    sdel_array_1D(nbSN); // nbSN可以在initzone后删除,在transgrid OutParGrid后可删去
    sdel_array_1D(nbZN); // nbSN可以在initzone后删除,在transgrid OutParGrid后可删去
    sdel_array_1D(nbRN); // nbSN可以在initzone后删除,在transgrid OutParGrid后可删去

    sdel_array_1D(nZIFace);

#ifdef MF_MPICH
    sdel_array_2D(bCNo, nNeighbor, false);
    sdel_array_2D(bFNo, nNeighbor, false);
#endif
    sdel_array_1D(nZINode);
    sdel_array_2D(bNSNo, nNeighborN, false);
    sdel_array_2D(bNRNo, nNeighborN, false);

    sdel_array_1D(Nmark);
    sdel_array_1D(WeightNodeDist); // 程序已删除
    sdel_array_2D(WeightNodeC2N);
    sdel_array_2D(WeightNodeN2C);
    sdel_array_2D(WeightNodeBFace2C);

    if (vgn)
    {
        sdel_array_1D(vgn);
        sdel_array_1D(BFacevgx);
        sdel_array_1D(BFacevgy);
        sdel_array_1D(BFacevgz);
        sdel_array_1D(vccx);
        sdel_array_1D(vccy);
        sdel_array_1D(vccz);
    }
    // sdel_array_1D(node_act);
    // sdel_array_1D(edge_act); // 程序已删除
    // sdel_array_1D(nChCellS);
    // sdel_array_1D(nChCellR);
    // sdel_array_1D(nChRNoComp);
    // sdel_array_1D(WeightNodeProl);

    // sdel_array_2D(Seta_Center);
    // sdel_array_2D(nChSNo);
    // sdel_array_2D(nChRNo);
    // sdel_array_2D(nChxcc);
    // sdel_array_2D(nChycc);
    // sdel_array_2D(nChzcc);

// // add by dingxin
// #ifdef REORDER
//     sdel_array_1D(order_cell_oTon);
//     sdel_array_1D(order_cell_nToo);
// #endif // REORDER

// #ifdef OMP_GroupColor
//     if (GroupColorSuccess)
//     {
//         sdel_array_1D(id_bface);
//         sdel_array_1D(id_iface);
//     }
// #endif
#ifdef OMP_DIVREP
    if (DivRepSuccess)
    {
        sdel_array_1D(idx_pthreads_bface);
        sdel_array_1D(id_division_bface);
        sdel_array_1D(idx_pthreads_iface);
        sdel_array_1D(id_division_iface);
        // #ifdef BoundedColoring
        //         sdel_array_1D(endIndex_bFace_vec);
        //         sdel_array_1D(endIndex_iFace_vec);
        // #endif ///< ~BoundedColoring
    }
#endif ///< ~OMP_DIVREP

    // #ifdef OMP_DIVCON
    //     tree_free(this->treeHead);
    //     sdel_array_1D(this->treeHead);
    // #endif // ~OMP_DIVCON

    // if (facecentroidskewness)
    // {
    //     sdel_array_1D(facecentroidskewness);
    //     sdel_array_1D(faceangleskewness); // 在EquiangleSkewnessSummary后不再调用
    //     sdel_array_1D(cellcentroidskewness);
    //     sdel_array_1D(cellvolsmoothness);
    //     sdel_array_1D(faceangle); // 在FaceAngleSummary后不再调用
    //     sdel_array_1D(cellwallnumber);
    // }
    // sdel_array_1D(Cell2Zone);
    // sdel_array_1D(cellwallnumber);
    // if (nNeighbor > 0)
    // {
    //     // sdel_array_1D(nbg);
    // }

    // #ifdef USING_PETSC
    //     sdel_array_1D(ghost2global);
    // #endif

    // // data for line implicit
    // sdel_array_1D(nCellsInLineForLI);
    // sdel_array_2D(CellsOfLinesForLI, nLinesForLI);
    // sdel_array_2D(FacesOfLinesForLI, nLinesForLI);
    // sdel_array_1D(LineIndexOfCellsForLI);

    // // 物面出发的线原有数目和线长,(在连接之前),用于CFL数放大使用
    // sdel_array_1D(n_cells_in_line_wall);
}

void PolyGrid::ComputeMetrics()
{
    IntType n = nTCell + nBFace; // interior+ghost cells

    // if (vol == 0)
    snew_array_1D(vol, n);

    // if (xcc == 0)
    snew_array_1D(xcc, n);
    // if (ycc == 0)
    snew_array_1D(ycc, n);
    // if (zcc == 0)
    snew_array_1D(zcc, n);

    // if (area == 0)
    snew_array_1D(area, nTFace);

    // if (xfn == 0)
    snew_array_1D(xfn, nTFace);
    // if (yfn == 0)
    snew_array_1D(yfn, nTFace);
    // if (zfn == 0)
    snew_array_1D(zfn, nTFace);

    // if (xfc == 0)
    snew_array_1D(xfc, nTFace);
    // if (yfc == 0)
    snew_array_1D(yfc, nTFace);
    // if (zfc == 0)
    snew_array_1D(zfc, nTFace);

    FaceCellCenterbyAverage(this, xfc, yfc, zfc, xcc, ycc, zcc);
    FaceAreaNormalCentroid_cycle(this, area, xfn, yfn, zfn, xfc, yfc, zfc);
    CellVolCentroid(this, vol, xcc, ycc, zcc);
    CorrectFaceNormal(this, xfn, yfn, zfn);
    CorrectCellCentroid(this, xcc, ycc, zcc, xfc, yfc, zfc, xfn, yfn, zfn);

#ifdef MF_MPICH
    // 体心，体积的并行传值
    CommInterfaceDataMPI(xcc);
    CommInterfaceDataMPI(ycc);
    CommInterfaceDataMPI(zcc);
    CommInterfaceDataMPI(vol);
#endif

    // ClosureCheck(this, xfn, area);
    // ClosureCheck(this, yfn, area);
    // ClosureCheck(this, zfn, area);
}

/*!
 * @brief       无网格重排序，此函数是为了编程统一
 *
 * @author
 * @date
 */
void PolyGrid::ReorderCellforLUSGS_0()
{
    IntType i;
    IntType n = nBFace + nTCell;

    snew_array_1D(LUSGSLayer, n);
    snew_array_1D(LUSGSCellOrder, nTCell);
    for (i = 0; i < n; i++)
    {
        LUSGSLayer[i] = i;
    }
    for (i = 0; i < nTCell; i++)
    {
        LUSGSCellOrder[i] = i;
    }
}

//////////////////////////////////////////////////////////////////////////
/**
struct xcc_idx
{
    double val;
    int idx;
};
*/

/**
template <typename data_type>
void swap(data_type &a, data_type &b)
{
    data_type t = a;
    a = b;
    b = t;
}

// < min_heap, after sorted: 321 top_k max
// > max_heap, after sorted: 123 top_k min
// #define COMPARE( a, b ) ( (a) < (b) )
// #define COMPARE2( a, b ) ( (a) < (b) )
#define COMPARE(a, b) ((a.val) < (b.val))

// k/2  2*k  2*k+1
#define PARENT(k) ((k) >> 1)
#define LEFT__(k) ((k) << 1)
#define RIGHT_(k) (((k) << 1) | 0x1)

template <typename data_type>
class heap_inplace
{
    // heap_max root >= child for k 2*k 2*k+1
public:
    // typedef size_t size_type;
    typedef int size_type;

private:
    data_type *_data;
    size_type _N;

public:
    heap_inplace(data_type *data, size_type N) : _data(data - 1), _N(N) {}

    // private:
    void swim(size_type k) // 上浮 自下而上的堆有序化
    {
        while ((k > 1) && COMPARE(_data[k], _data[PARENT(k)])) // 第 k 个元素不是根节点 且 比其父节点大
        {
            swap(_data[PARENT(k)], _data[k]); // 交换 当前节点 与 其父节点，将 当前节点 上浮一层
            k = PARENT(k);
        }
    }

    void sink(size_type k, size_type N) // 下沉 自上而下的堆有序化
    {
        while (LEFT__(k) <= N) // 该节点存在(左)子节点  第 k 个元素不是根节点 且 比其父节点大
        {
            int largest_child_idx = LEFT__(k); // 该节点下 较大 子节点

            // 该节点存在 右子节点 且较大
            if (largest_child_idx < N && COMPARE(_data[largest_child_idx + 1], _data[largest_child_idx]))
                ++largest_child_idx;

            // 当前节点 不小于 子节点，则停止下沉
            if (!COMPARE(_data[largest_child_idx], _data[k]))
                break;

            swap(_data[k], _data[largest_child_idx]); // 交换 当前节点 与 较大子节点，将 当前节点 下沉一层
            k = largest_child_idx;
        }
    }

public:
    void sort()
    {
        size_type N = _N; // get_size();
        for (size_type i = N / 2; i > 0; --i)
        {
            sink(i, N);
        }

        while (N > 1)
        {
            swap(_data[1], _data[N]);
            --N;
            sink(1, N);
        }
    }

    size_type get_size() const { return _N; }

    void add(const data_type &data) { _data[++_N] = data; }
    void set(size_type i, const data_type &data) { _data[i] = data; }
    const data_type &get(size_type i) { return _data[i]; }
};

#undef PARENT
#undef LEFT__
#undef RIGHT_

// data : to be parse, N size
// result, k size
template <typename data_type, typename size_type, size_t k>
void top_k(const data_type *data, size_type N, data_type *result)
{
    heap_inplace<data_type> max_heap(result, 0);

    for (size_type i = 0; i < N; ++i)
    {
        size_type size = max_heap.get_size();

        if (size < k) // max_heap 未满
        {
            max_heap.add(data[i]);
            max_heap.swim(size + 1);
        }
        else
        {
            const data_type &max_element = max_heap.get(1);

            if (COMPARE(max_element, data[i])) // 新元素比堆中最大的元素小 //if( data[i] < max_element )
            {
                max_heap.set(1, data[i]); // 替换根节点元素

                max_heap.sink(1, size);
            }
        }
    }
}

template <typename data_type, typename size_type>
void top_k(const data_type *data, size_type N, data_type *result, size_t k)
{
    heap_inplace<data_type> max_heap(result, 0);

    for (size_type i = 0; i < N; ++i)
    {
        size_type size = max_heap.get_size();

        if (size < k) // max_heap 未满
        {
            max_heap.add(data[i]);
            max_heap.swim(size + 1);
        }
        else
        {
            const data_type &max_element = max_heap.get(1);

            if (COMPARE(max_element, data[i])) // 新元素比堆中最大的元素小 //if( data[i] < max_element )
            {
                max_heap.set(1, data[i]); // 替换根节点元素

                max_heap.sink(1, size);
            }
        }
    }

    max_heap.sort();
}

#undef COMPARE
//////////////////////////////////////////////////////////////////////////

RealGeom aspect_ratio(RealGeom dist, RealGeom area, RealGeom span_len, IntType nnpf)
{
    const RealGeom SQRT3 = 1.732050807568877;

    RealGeom length = 1.0;

    if (span_len > TINY)
    { // 2d
        if (nnpf == 4)
        { // quadrilateral
            length = area / span_len;
        }
        else
        {
            cout << endl
                 << "Need new code in function aspect_ratio!" << endl;
        }
    }
    else
    { // 3d
        if (nnpf == 3)
        { // triangle
            length = sqrt(area * 4.0 / SQRT3);
        }
        else if (nnpf == 4)
        { // quadrilateral
            length = sqrt(area);
        }
        else
        {
            cout << endl
                 << "Need new code in function aspect_ratio!" << endl;
        }
    }

    return length / (dist + TINY);
}
*/

/**
void PolyGrid::ReorderCellforLUSGS_101()
{
    IntType nTCell = GetNTCell();
    IntType nBFace = GetNBFace();
    IntType n = nTCell + nBFace;
    IntType nTNode = GetNTNode();
    IntType *f2c = Getf2c();
    IntType *nNPF = GetnNPF();
    IntType **F2N = CalF2N(this);
    IntType *nFPC = CalnFPC(this);
    IntType **C2F = CalC2F(this);
    IntType *nNPC = CalnNPC(this);

    RealGeom *xfn = GetXfn();
    RealGeom *yfn = GetYfn();
    RealGeom *zfn = GetZfn();

    RealGeom *xfc = GetXfc();
    RealGeom *yfc = GetYfc();
    RealGeom *zfc = GetZfc();

    IntType *LineOfCell = NULL; // line index of each cell
    snew_array_1D(LineOfCell, n);

    bool *mark_cell = NULL;
    snew_array_1D(mark_cell, nTCell);
    for (IntType i = 0; i < nTCell; i++)
    {
        mark_cell[i] = true;
    }

    vector<IntType> Line_tmp;
    vector<IntType> CellOfLine_tmp; // cell index of each line's cells
    vector<IntType> FaceOfLine_tmp; // face index of each line's faces
    IntType iLine = 0;

    RealFlow AR_cutoff = 1.0;
    GetData(&AR_cutoff, REAL_FLOW, 1, "AR_cutoff", 0);
    IntType max_length_of_line = 0xffff; // 20;//
    GetData(&max_length_of_line, INT, 1, "max_length_of_line", 0);
    IntType min_length_of_line = 20; // 3,4, 8
    GetData(&min_length_of_line, INT, 1, "min_length_of_line", 0);

    //////////////////////////////////////////////////////////////////////////
    IntType is_search_concave_region = 1;
    GetData(&is_search_concave_region, INT, 1, "is_search_concave_region", 0);

    vector<IntType> wall_face_tmp;

    for (IntType face = 0; face < nBFace; ++face)
    {
        if (bcr[face]->GetType() == WALL)
            wall_face_tmp.push_back(face);
    }

    IntType nTFace = GetNTFace();
    RealGeom *norm_dist_c2c = (RealGeom *)GetDataPtr(REAL_GEOM, nTFace, "norm_dist_c2c");
    RealGeom *area = GetFaceArea();

    {
        xcc_idx *data = new xcc_idx[wall_face_tmp.size()];

        for (IntType i = 0; i < wall_face_tmp.size(); ++i)
        {
            IntType face = wall_face_tmp[i];
            data[i].val = aspect_ratio(norm_dist_c2c[face], area[face], span_len_of_2d, nNPF[face]);
            data[i].idx = face;
        }

        heap_inplace<xcc_idx> tmp(data, wall_face_tmp.size());
        tmp.sort();

        for (IntType i = 0; i < wall_face_tmp.size(); ++i)
        {
            wall_face_tmp[i] = data[i].idx;
        }

        delete[] data;
    }

    if (is_search_concave_region) // new method consider concave region and AR_cutoff
    {
        // 寻找物面附近的线(一般为半结构网格， 线太短时强制推min_length_of_line层
        for (IntType i = 0; i < wall_face_tmp.size(); ++i)
        {
            IntType face1 = wall_face_tmp[i];
            IntType cell = f2c[face1 + face1];

            RealFlow nx = -xfn[face1];
            RealFlow ny = -yfn[face1];
            RealFlow nz = -zfn[face1];

            RealGeom AR_prev = BIG;

            IntType count = 0;

            while (true)
            {
                if ((cell >= nTCell) || (!mark_cell[cell]) || (count >= max_length_of_line))
                { // ghost cell or marked cell
                    if (count != 0)
                    {
                        Line_tmp.push_back(count);
                        ++iLine;
                    }
                    break;
                }

                ++count;
                LineOfCell[cell] = iLine;
                CellOfLine_tmp.push_back(cell);
                FaceOfLine_tmp.push_back(face1);
                mark_cell[cell] = false;

                if (count < min_length_of_line)
                { // 物面出发的线太短，即使不是三棱柱和六面体也强制找线，规则：面积最大，AR不能太小。

                    set<IntType> face1_nodes;
                    for (IntType k = 0; k < nNPF[face1]; ++k)
                        face1_nodes.insert(F2N[face1][k]);

                    bool found_flag = false;

                    for (IntType j = 0; j < nFPC[cell]; ++j)
                    { // 在非粗化的网格上，仅有三棱柱和六面体才会找到“对面”
                        IntType face2 = C2F[cell][j];

                        if (face2 == face1)
                            continue;

                        set<IntType> tot_face_nodes(face1_nodes);
                        for (IntType k = 0; k < nNPF[face2]; ++k)
                            tot_face_nodes.insert(F2N[face2][k]);

                        if (tot_face_nodes.size() == nNPF[face1] + nNPF[face2])
                        { // found opposite face2
                            face1 = face2;
                            cell = f2c[face1 + face1] + f2c[face1 + face1 + 1] - cell;
                            found_flag = true;
                            break;
                        }
                    } // end face of cell

                    if (!found_flag)
                    { // not found opposite face

                        RealGeom max_area = -1.0;
                        IntType face_of_max = -1;

                        for (IntType j = 0; j < nFPC[cell]; ++j)
                        {
                            IntType face2 = C2F[cell][j];
                            if (face2 == face1)
                                continue;

                            IntType cell2 = f2c[face2 + face2] + f2c[face2 + face2 + 1] - cell;

                            if ((cell2 >= nTCell) || (!mark_cell[cell2]))
                                continue; // ghost cell or marked cell

                            RealGeom t_area = area[face2];

                            RealFlow tmp = xfn[face2] * nx + yfn[face2] * ny + zfn[face2] * nz;

                            if (f2c[face2 + face2] != cell)
                            {
                                tmp = -tmp;
                            }

                            if (tmp > TINY)
                            {
                                if (t_area > max_area)
                                {
                                    max_area = t_area;
                                    face_of_max = face2;
                                }
                            }
                        } // end face of cell

                        if (max_area > TINY)
                        {
                            RealGeom dist = (xfc[face_of_max] - xfc[face1]) * xfn[face_of_max] + (yfc[face_of_max] - yfc[face1]) * yfn[face_of_max] + (zfc[face_of_max] - zfc[face1]) * zfn[face_of_max];
                            // dist = abs(dist);
                            dist = fabs(dist);

                            RealGeom t_AR = aspect_ratio(dist, max_area, span_len_of_2d, nNPF[face_of_max]);

                            if (t_AR >= 3.0 * AR_cutoff)
                            {
                                face1 = face_of_max;
                                cell = f2c[face1 + face1] + f2c[face1 + face1 + 1] - cell;
                            }
                        }

                    } // end not found
                }
                else
                { // 边界层网格增长率一般工程上取1.15-1.25 AR减小取

                    set<IntType> face1_nodes;
                    for (IntType k = 0; k < nNPF[face1]; ++k)
                        face1_nodes.insert(F2N[face1][k]);

                    for (IntType j = 0; j < nFPC[cell]; ++j)
                    {
                        IntType face2 = C2F[cell][j];

                        if (face2 == face1)
                            continue;

                        set<IntType> tot_face_nodes(face1_nodes);
                        for (IntType k = 0; k < nNPF[face2]; ++k)
                            tot_face_nodes.insert(F2N[face2][k]);

                        if (tot_face_nodes.size() == nNPF[face1] + nNPF[face2])
                        { // found opposite face2

                            RealGeom dist = (xfc[face2] - xfc[face1]) * xfn[face2] + (yfc[face2] - yfc[face1]) * yfn[face2] + (zfc[face2] - zfc[face1]) * zfn[face2];
                            // dist = abs(dist);
                            dist = fabs(dist);

                            RealGeom t_AR = aspect_ratio(dist, area[face2], span_len_of_2d, nNPF[face2]); // norm_dist_c2c[face2]

                            if (t_AR >= AR_cutoff)
                            {
                                if (t_AR >= 3.0 * AR_cutoff)
                                {
                                    face1 = face2;
                                    cell = f2c[face1 + face1] + f2c[face1 + face1 + 1] - cell;
                                    AR_prev = t_AR;
                                }
                                else
                                {
                                    if (t_AR < AR_prev * 1.001)
                                    {
                                        face1 = face2;
                                        cell = f2c[face1 + face1] + f2c[face1 + face1 + 1] - cell;
                                        AR_prev = t_AR;
                                    }
                                }
                            }

                            break;
                        }

                    } // end of face of cell

                } // if(count >= min_length_of_line )
            }     // end of search for this line

        } // end of nbface
    }
    else
    {

        // 寻找物面附近的线(半结构网格
        for (IntType i = 0; i < wall_face_tmp.size(); ++i)
        {
            IntType face1 = wall_face_tmp[i];
            IntType cell = f2c[face1 + face1];
            /// for(IntType i=0; i<nBFace; ++i) {
            ///     IntType type = bcr[i]->GetType();
            ///     if (type != WALL) continue;

            /// IntType cell = f2c[i+i];
            /// IntType face1 = i;
            IntType count = 0;

            while (true)
            {
                if ((cell >= nTCell) || (!mark_cell[cell]))
                { // ghost cell or marked cell
                    if (count != 0)
                    {
                        Line_tmp.push_back(count);
                        ++iLine;
                    }
                    break;
                }

                if (count >= max_length_of_line)
                {
                    Line_tmp.push_back(count);
                    ++iLine;
                    break;
                }

                ++count;
                LineOfCell[cell] = iLine;
                CellOfLine_tmp.push_back(cell);
                FaceOfLine_tmp.push_back(face1);
                mark_cell[cell] = false;

                if (nNPC[cell] != 6 && nNPC[cell] != 8)
                { // only prism and hexa can continue search, so stop 实际上对于找”对面“的算法，不用此判断也可以
                    Line_tmp.push_back(count);
                    ++iLine;
                    break;
                }

                for (IntType j = 0; j < nFPC[cell]; ++j)
                {
                    IntType face2 = C2F[cell][j];

                    bool key = false;
                    for (IntType k = 0; k < nNPF[face1]; ++k)
                    {

                        IntType p1 = F2N[face1][k];
                        for (IntType l = 0; l < nNPF[face2]; ++l)
                        {
                            IntType p2 = F2N[face2][l];
                            if (p2 == p1)
                            {
                                key = true;
                                break;
                            }
                        }

                        if (key)
                            break;
                    }

                    if (!key)
                    { // found
                        // if(count >= min_length_of_line ) {
                        //     RealGeom dist = (xfc[face2]-xfc[face1])*xfn[face1] + (yfc[face2]-yfc[face1])*yfn[face1] + (zfc[face2]-zfc[face1])*zfn[face1];
                        //     dist = abs( dist );
                        //     RealGeom len_of_face1 = 1; //face1的参考长度，怎么计算？二维网格如何处理？
                        //     if( len_of_face1/(dist+TINY) < AR_cutoff ) {
                        //         break;
                        //     }
                        // }

                        // if(count < min_length_of_line ) {
                        face1 = face2;
                        cell = f2c[face1 + face1] + f2c[face1 + face1 + 1] - cell;
                        break;
                        //}
                    }
                } // end of face of cell
            }     // end of search for this line
        }         // end of nbface
    }

    IntType nFaceOnWallLines = iLine;
    IntType nFaceOnWallCells = CellOfLine_tmp.size();
    // #ifdef MF_MPICH
    //     {
    //         IntType nFaceOnWallLines_glb;
    //         MPI_Allreduce(&nFaceOnWallLines, &nFaceOnWallLines_glb, 1, MPIIntType, MPI_SUM, MPI_COMM_WORLD);
    //         IntType nFaceOnWallCells_glb;
    //         MPI_Allreduce(&nFaceOnWallCells, &nFaceOnWallCells_glb, 1, MPIIntType, MPI_SUM, MPI_COMM_WORLD);
    //         if (nFaceOnWallLines_glb)
    //         {
    //             mflog::log.set_one_processor_out();
    //             mflog::log << endl
    //                        << "\nnumber of Lines (from wall face) " << nFaceOnWallLines_glb << ", average length " << nFaceOnWallCells_glb / nFaceOnWallLines_glb << endl;
    //         }
    //     }
    // #else
    if (nFaceOnWallLines)
    {
        cout << "\nnumber of Lines (from wall face) " << nFaceOnWallLines << ", average length " << nFaceOnWallCells / nFaceOnWallLines << endl;
    }
    // #endif

    IntType is_search_wall_pointandedge = 1;
    GetData(&is_search_wall_pointandedge, INT, 1, "is_search_wall_pointandedge", 0);

    if (is_search_wall_pointandedge)
    { // 尖后缘找线

        bool *node_on_wall_mark = NULL;
        snew_array_1D(node_on_wall_mark, nTNode);
        for (IntType i = 0; i < nTNode; ++i)
        {
            node_on_wall_mark[i] = false;
        }

        for (IntType i = 0; i < wall_face_tmp.size(); ++i)
        { // 找所有物面边界上的点（不重复）
            IntType face = wall_face_tmp[i];

            for (IntType j = 0; j < nNPF[face]; ++j)
                node_on_wall_mark[F2N[face][j]] = true;
        }

        IntType **C2N = CalC2N(this);

        vector<IntType> wall_cells;
        for (IntType cell = 0; cell < nTCell; ++cell)
        {
            if (mark_cell[cell])
            {
                for (IntType i = 0; i < nNPC[cell]; ++i)
                {
                    if (node_on_wall_mark[C2N[cell][i]])
                    {
                        wall_cells.push_back(cell);
                        break;
                    }
                }
            }
        }

        for (IntType i = 0; i < wall_cells.size(); ++i)
        { // 和物面共点，棱的单元(共面的单元已经由前一步找到并标记)
            IntType cell1 = wall_cells[i];
            if (!mark_cell[cell1])
            { // marked cell1
                // cout << "error in search line of bound nodes\n"; //测试中出现了此种情况，实际上此单元被前面的线（和物面共点、棱的线）用了
                continue;
            }

            IntType count = 1;
            // -1 cell1 face12 cell2 ...
            LineOfCell[cell1] = iLine;
            CellOfLine_tmp.push_back(cell1);
            FaceOfLine_tmp.push_back(-1);
            mark_cell[cell1] = false;

            IntType cell = cell1;
            IntType face1 = -1;

            {
                vector<IntType> candidate_faces;
                vector<IntType> candidate_cells;

                for (IntType j = 0; j < nFPC[cell1]; ++j)
                {
                    IntType face2 = C2F[cell1][j];

                    bool key = true;
                    for (IntType k = 0; k < nNPF[face2]; ++k)
                    {
                        if (node_on_wall_mark[F2N[face2][k]])
                        { // face2有点在物面上，不考虑
                            key = false;
                            break;
                        }
                    }

                    if (key)
                    { // found
                        IntType cell_tmp = f2c[face2 + face2] + f2c[face2 + face2 + 1] - cell1;
                        if ((cell_tmp < nTCell) && (mark_cell[cell_tmp]))
                        {
                            candidate_faces.push_back(face2);
                            candidate_cells.push_back(cell_tmp);
                        }
                    }
                } // end of face of cell1

                if (candidate_faces.size() == 0)
                {
                    ;
                }
                else if (candidate_faces.size() == 1)
                {
                    face1 = candidate_faces[0];
                    cell = candidate_cells[0];
                }
                else
                { // 有多个可选，一般为2，比如平板前面那个六面体单元。此时选择一个 “更内部的单元” or 选择 AR大的 找线 if ( candidate_faces.size() > 1 )
                    /// cout << "number of candidate cells (from point and edge on wall) " << candidate_cells.size() << endl;

                    IntType max_inner_cells = -1;
                    IntType cell_idx = 0;

                    for (IntType j = 0; j < candidate_cells.size(); ++j)
                    {
                        IntType cell_tmp = candidate_cells[j];

                        IntType n_inner_cells_tmp = 0;

                        for (IntType k = 0; k < nFPC[cell_tmp]; ++k)
                        {
                            IntType face_tmp = C2F[cell_tmp][k];

                            /// if(( f2c[face_tmp+face_tmp] + f2c[face_tmp+face_tmp+1] - cell_tmp ) < nTCell) ++n_inner_cells_tmp;
                            if (face_tmp >= nBFace)
                                ++n_inner_cells_tmp;
                        }

                        if (n_inner_cells_tmp > max_inner_cells)
                        {
                            max_inner_cells = n_inner_cells_tmp;
                            cell_idx = j;
                        }
                    }

                    face1 = candidate_faces[cell_idx];
                    cell = candidate_cells[cell_idx];
                }
            }

            while (true)
            {
                if ((cell >= nTCell) || (!mark_cell[cell]))
                { // ghost cell or marked cell
                    if (count != 0)
                    {
                        Line_tmp.push_back(count);
                        ++iLine;
                    }
                    break;
                }

                if (count >= max_length_of_line)
                {
                    Line_tmp.push_back(count);
                    ++iLine;
                    break;
                }

                ++count;
                LineOfCell[cell] = iLine;
                CellOfLine_tmp.push_back(cell);
                FaceOfLine_tmp.push_back(face1);
                mark_cell[cell] = false;

                set<IntType> face1_nodes;
                for (IntType k = 0; k < nNPF[face1]; ++k)
                    face1_nodes.insert(F2N[face1][k]);

                for (IntType j = 0; j < nFPC[cell]; ++j)
                {
                    IntType face2 = C2F[cell][j];

                    if (face2 == face1)
                        continue;

                    set<IntType> tot_face_nodes(face1_nodes);
                    for (IntType k = 0; k < nNPF[face2]; ++k)
                        tot_face_nodes.insert(F2N[face2][k]);

                    if (tot_face_nodes.size() == nNPF[face1] + nNPF[face2])
                    { // found opposite face2

                        if (count < min_length_of_line)
                        {
                            face1 = face2;
                            cell = f2c[face1 + face1] + f2c[face1 + face1 + 1] - cell;
                            break;
                        }
                        else
                        {

                            RealGeom dist = (xfc[face2] - xfc[face1]) * xfn[face2] + (yfc[face2] - yfc[face1]) * yfn[face2] + (zfc[face2] - zfc[face1]) * zfn[face2];
                            // dist = abs(dist);
                            dist = fabs(dist);

                            RealGeom t_AR = aspect_ratio(dist, area[face2], span_len_of_2d, nNPF[face2]); // norm_dist_c2c[face2]

                            if (t_AR >= AR_cutoff)
                            {
                                face1 = face2;
                                cell = f2c[face1 + face1] + f2c[face1 + face1 + 1] - cell;
                                break;
                            }
                        }
                    }
                } // end of face of cell

            } // end of search for this line while(true)

        } // end of wall_cells cell with wall node/edge

        sdel_array_1D(node_on_wall_mark);

        {
            IntType pointwallline = iLine - nFaceOnWallLines;
            // #ifdef MF_MPICH
            //             IntType pointwallline_glb;
            //             MPI_Allreduce(&pointwallline, &pointwallline_glb, 1, MPIIntType, MPI_SUM, MPI_COMM_WORLD);
            //             if (pointwallline_glb)
            //             {
            //                 mflog::log.set_one_processor_out();
            //                 mflog::log << endl
            //                            << "number of Lines (from point and/or edge on wall) " << pointwallline_glb << endl;
            //             }
            // #else
            if (pointwallline)
            {
                cout << "number of Lines (from point and/or edge on wall) " << pointwallline << endl;
            }
            // #endif
        }
    } // end if() 尖后缘找线

    RealFlow rho00, u00, v00, w00;
    GetData(&rho00, REAL_FLOW, 1, "rho");
    GetData(&u00, REAL_FLOW, 1, "u");
    GetData(&v00, REAL_FLOW, 1, "v");
    GetData(&w00, REAL_FLOW, 1, "w");

    IntType is_sort_lines_from_wall = 1;
    GetData(&is_sort_lines_from_wall, INT, 1, "is_sort_lines_from_wall", 0);

    if (is_sort_lines_from_wall && iLine > 0) // 物面线排序
    {
        RealGeom *xcc = GetXcc();
        RealGeom *ycc = GetYcc();
        RealGeom *zcc = GetZcc();

        xcc_idx *data = new xcc_idx[iLine];

        IntType count = 0;
        for (IntType i = 0; i < iLine; ++i)
        {
            /// data[i].val = -xcc[ CellOfLine_tmp[count] ];
            data[i].val = -(xcc[CellOfLine_tmp[count]] * u00 + ycc[CellOfLine_tmp[count]] * v00 + zcc[CellOfLine_tmp[count]] * w00);
            data[i].idx = i;

            count += Line_tmp[i];
        }

        heap_inplace<xcc_idx> tmp(data, iLine);
        tmp.sort();

        // IntType *LineOfCell2 = NULL;
        // snew_array_1D(LineOfCell2, n);

        IntType *old_elem_start_idx = new IntType[iLine];

        old_elem_start_idx[0] = 0;
        for (IntType i = 1; i < iLine; ++i)
        {
            old_elem_start_idx[i] = old_elem_start_idx[i - 1] + Line_tmp[i - 1];
        }

        vector<IntType> Line_tmp_old = Line_tmp;
        vector<IntType> CellOfLine_tmp_old = CellOfLine_tmp;
        vector<IntType> FaceOfLine_tmp_old = FaceOfLine_tmp;

        count = 0;
        for (IntType i = 0; i < iLine; ++i)
        {
            IntType old_idx = data[i].idx;
            Line_tmp[i] = Line_tmp_old[old_idx];

            for (IntType j = 0; j < Line_tmp[i]; j++)
            {
                CellOfLine_tmp[count] = CellOfLine_tmp_old[old_elem_start_idx[old_idx] + j];
                LineOfCell[CellOfLine_tmp[count]] = i;
                FaceOfLine_tmp[count] = FaceOfLine_tmp_old[old_elem_start_idx[old_idx] + j];

                ++count;
            }
        }

        // sdel_array_1D(LineOfCell);
        // LineOfCell = LineOfCell2;

        delete[] old_elem_start_idx;
        delete[] data;
    }

    // 存储物面原有线的数量和长度，用于附面层CFL放大处理
    const IntType nLinesInWall_start = 0;
    SetLI_nLine_wall_start(nLinesInWall_start);

    const IntType nLinesInWall = iLine - nLinesInWall_start;
    SetLI_nLine_wall(nLinesInWall);

    IntType *nCellsInLinesWall = NULL;
    snew_array_1D(nCellsInLinesWall, nLinesInWall);
    for (IntType i = 0; i < nLinesInWall; ++i)
    {
        nCellsInLinesWall[i] = Line_tmp[i];
    }
    SetLI_nCellinLine_wall(nCellsInLinesWall);

    RealFlow norm_of_rq = rho00 * sqrt(u00 * u00 + v00 * v00 + w00 * w00);

    IntType is_search_lines_faraway = 1;

    RealFlow mach;
    GetData(&mach, REAL_FLOW, 1, "mach");

    GetData(&is_search_lines_faraway, INT, 1, "is_search_lines_faraway", 0);

    if (is_search_lines_faraway)
    { // 外围找线

        queue<IntType> seeds_of_outer_line;

        // Get flow variables
        RealFlow *q[5];
        q[0] = (RealFlow *)GetDataPtr(REAL_FLOW, n, "rho");
        q[1] = (RealFlow *)GetDataPtr(REAL_FLOW, n, "u");
        q[2] = (RealFlow *)GetDataPtr(REAL_FLOW, n, "v");
        q[3] = (RealFlow *)GetDataPtr(REAL_FLOW, n, "w");
        // q[4] = (RealFlow *)GetDataPtr(REAL_FLOW, n, "p");

        RealGeom *area = GetFaceArea();
        // RealFlow norm_of_rq = rho00 * sqrt( u00*u00 + v00*v00 + w00*w00 );

        RealGeom *xcc = GetXcc();
        RealGeom *ycc = GetYcc();
        RealGeom *zcc = GetZcc();

        // 入口边界找种子
        {
            vector<xcc_idx> seeds_before_sorted;

            for (IntType i = 0; i < nBFace; i++)
            {
                IntType type = bcr[i]->GetType();
                if (type == WALL)
                    continue;

                IntType cell = f2c[i + i];

                RealFlow rho = rho00, u = u00, v = v00, w = w00;
                if (q[0] != NULL)
                {
                    rho = q[0][cell];
                    u = q[1][cell];
                    v = q[2][cell];
                    w = q[3][cell];
                }

                RealFlow tmp = rho * (xfn[i] * u + yfn[i] * v + zfn[i] * w);

                if (tmp / norm_of_rq < -1.0e-14)
                {
                    xcc_idx tmp2;
                    tmp2.val = area[i] * tmp * tmp;
                    tmp2.idx = cell;
                    seeds_before_sorted.push_back(tmp2);
                }
            } // end of bdface

            heap_inplace<xcc_idx> tmp(&seeds_before_sorted[0], seeds_before_sorted.size());
            tmp.sort();

            for (IntType i = 0; i < seeds_before_sorted.size(); ++i)
            {
                seeds_of_outer_line.push(seeds_before_sorted[i].idx);
                /// cout << seeds_before_sorted[i].idx << " " << xcc[seeds_before_sorted[i].idx] << " "<< ycc[seeds_before_sorted[i].idx] << " "<< zcc[seeds_before_sorted[i].idx] << endl;
            }
        }

        // 其余内部单元也暂时作为种子
        for (IntType i_cell = 0; i_cell < nTCell; ++i_cell)
        { // out of boundary condition
            if (mark_cell[i_cell])
            {
                seeds_of_outer_line.push(i_cell); // 部分重复，先不管了
            }
        }

        while (!seeds_of_outer_line.empty())
        {
            IntType cell = seeds_of_outer_line.front();
            seeds_of_outer_line.pop();

            if ((cell >= nTCell) || (!mark_cell[cell]))
            { // ghost cell or marked cell
                continue;
            }

            vector<IntType> CellOfLine_forward;
            vector<IntType> FaceOfLine_forward;
            vector<IntType> CellOfLine_backward;
            vector<IntType> FaceOfLine_backward;

            // IntType face1 = -1;
            IntType count = 0;

            ++count;
            LineOfCell[cell] = iLine;
            // CellOfLine_forward.push_back(cell);
            // FaceOfLine_forward.push_back(face1);
            mark_cell[cell] = false;

            IntType cell_tmp = cell;
            // forward
            while (true)
            {
                vector<xcc_idx> faces_may_in_line;

                for (IntType j = 0; j < nFPC[cell_tmp]; ++j)
                {
                    IntType face2 = C2F[cell_tmp][j];
                    IntType cell2 = f2c[face2 + face2] + f2c[face2 + face2 + 1] - cell_tmp;

                    if ((cell2 >= nTCell) || (!mark_cell[cell2]))
                    { // ghost cell or marked cell
                        continue;
                    }

                    RealFlow ru = rho00 * u00, rv = rho00 * v00, rw = rho00 * w00;
                    if (q[0] != NULL)
                    {
                        ru = 0.5 * (q[0][cell_tmp] * q[1][cell_tmp] + q[0][cell2] * q[1][cell2]);
                        rv = 0.5 * (q[0][cell_tmp] * q[2][cell_tmp] + q[0][cell2] * q[2][cell2]);
                        rw = 0.5 * (q[0][cell_tmp] * q[3][cell_tmp] + q[0][cell2] * q[3][cell2]);
                    }

                    RealFlow tmp = xfn[face2] * ru + yfn[face2] * rv + zfn[face2] * rw;

                    if (f2c[face2 + face2] != cell_tmp)
                    {
                        tmp = -tmp;
                    }

                    if (tmp / norm_of_rq > 1.0e-3)
                    {
                        xcc_idx tmp2;
                        tmp2.val = area[face2] * tmp * tmp;
                        tmp2.idx = face2;
                        faces_may_in_line.push_back(tmp2);
                    }
                } // end of face of cell

                if (faces_may_in_line.size() == 0)
                { // not found
                    break;
                }
                else
                {
                    heap_inplace<xcc_idx> tmp3(&faces_may_in_line[0], faces_may_in_line.size());
                    tmp3.sort();

                    IntType face2 = faces_may_in_line[0].idx;
                    cell_tmp = f2c[face2 + face2] + f2c[face2 + face2 + 1] - cell_tmp;

                    ++count;
                    LineOfCell[cell_tmp] = iLine;
                    CellOfLine_forward.push_back(cell_tmp);
                    FaceOfLine_forward.push_back(face2);
                    mark_cell[cell_tmp] = false;
                }
            } // end of search forward

            cell_tmp = cell;
            // backward
            while (true)
            {
                vector<xcc_idx> faces_may_in_line;

                for (IntType j = 0; j < nFPC[cell_tmp]; ++j)
                {
                    IntType face2 = C2F[cell_tmp][j];
                    IntType cell2 = f2c[face2 + face2] + f2c[face2 + face2 + 1] - cell_tmp;

                    if ((cell2 >= nTCell) || (!mark_cell[cell2]))
                    { // ghost cell or marked cell
                        continue;
                    }

                    RealFlow ru = rho00 * u00, rv = rho00 * v00, rw = rho00 * w00;
                    if (q[0] != NULL)
                    {
                        ru = 0.5 * (q[0][cell_tmp] * q[1][cell_tmp] + q[0][cell2] * q[1][cell2]);
                        rv = 0.5 * (q[0][cell_tmp] * q[2][cell_tmp] + q[0][cell2] * q[2][cell2]);
                        rw = 0.5 * (q[0][cell_tmp] * q[3][cell_tmp] + q[0][cell2] * q[3][cell2]);
                    }

                    RealFlow tmp = xfn[face2] * ru + yfn[face2] * rv + zfn[face2] * rw;

                    if (f2c[face2 + face2] != cell_tmp)
                    {
                        tmp = -tmp;
                    }

                    if (tmp / norm_of_rq < -1.0e-3)
                    {
                        xcc_idx tmp2;
                        tmp2.val = area[face2] * tmp * tmp;
                        tmp2.idx = face2;
                        faces_may_in_line.push_back(tmp2);
                    }
                } // end of face of cell

                if (faces_may_in_line.size() == 0)
                { // not found
                    break;
                }
                else
                {
                    heap_inplace<xcc_idx> tmp3(&faces_may_in_line[0], faces_may_in_line.size());
                    tmp3.sort();

                    IntType face2 = faces_may_in_line[0].idx;
                    cell_tmp = f2c[face2 + face2] + f2c[face2 + face2 + 1] - cell_tmp;

                    ++count;
                    LineOfCell[cell_tmp] = iLine;
                    CellOfLine_backward.push_back(cell_tmp);
                    FaceOfLine_backward.push_back(face2);
                    mark_cell[cell_tmp] = false;
                }
            } // end of search backward

            Line_tmp.push_back(count);
            ++iLine;
            FaceOfLine_tmp.push_back(-1);

            if (CellOfLine_backward.size() > 0)
            {
                for (int j = CellOfLine_backward.size() - 1; j >= 0; --j)
                {
                    CellOfLine_tmp.push_back(CellOfLine_backward[j]);
                    FaceOfLine_tmp.push_back(FaceOfLine_backward[j]);
                }
            }

            CellOfLine_tmp.push_back(cell);

            for (int j = 0; j < CellOfLine_forward.size(); ++j)
            {
                CellOfLine_tmp.push_back(CellOfLine_forward[j]);
                FaceOfLine_tmp.push_back(FaceOfLine_forward[j]);
            }
        } // end of seeds
    }
    else
    {
        // if(0) //已经没有必要，因为所有单元都标记了 （外围线）
        for (IntType i_cell = 0; i_cell < nTCell; ++i_cell)
        { // out of boundary condition
            if (mark_cell[i_cell])
            {
                Line_tmp.push_back(1);
                LineOfCell[i_cell] = iLine++;
                CellOfLine_tmp.push_back(i_cell);
                FaceOfLine_tmp.push_back(-1);
                /// mark_cell[i_cell] = false;
            }
        }
    }

    sdel_array_1D(mark_cell);

    /// Wisces: 先注释掉

    // IntType is_debug_for_lineimplicit = 0;
    // GetData(&is_debug_for_lineimplicit, INT, 1, "is_debug_for_lineimplicit", 0);

    // if (is_debug_for_lineimplicit)
    //     do
    //     { // 按照线长 排序调试
    //         if (mflog::log.rank_id() == 1)
    //         {
    //             xcc_idx *data = new xcc_idx[iLine];

    //             for (IntType i = 0; i < iLine; ++i)
    //             {
    //                 data[i].val = Line_tmp[i];
    //                 data[i].idx = Line_tmp[i]; // i
    //             }

    //             heap_inplace<xcc_idx> tmp(data, iLine);
    //             tmp.sort();

    //             for (IntType i = 0; i < iLine; ++i)
    //             {
    //                 if (fabs(data[i].val - 1.0) < 0.001)
    //                 {
    //                     cout << "\nline length 1 * " << iLine - i << endl;
    //                     break;
    //                 }
    //                 cout << (IntType)data[i].val << " ";
    //                 if ((i % 16) == 15)
    //                     cout << "\n";
    //             }
    //             cout << "\n";
    //             delete[] data;
    //         }

    //         std::string file_name = "line_info_zone";
    //         // #ifdef MF_MPICH
    //         //             /// file_name += std::to_string( mflog::log.rank_id() );
    //         //             file_name += int2str(mflog::log.rank_id());
    //         // #endif
    //         file_name += ".dat";

    //         FILE *sgrid = fopen(file_name.c_str(), "w");
    //         // Open for read (will fail if file does not exist)
    //         if (sgrid == 0)
    //         {
    //             printf("The file nodetopo.dat was not opened\n");
    //             break;
    //         }

    //         fprintf(sgrid, "TITLE = \"3D Unstructured Hybrid Mesh Line Info\"\n");
    //         fprintf(sgrid, "FILETYPE = FULL\n");
    //         fprintf(sgrid, "VARIABLES = x, y, z\n");

    //         IntType count = 0;
    //         for (IntType i = 0; i < iLine; ++i)
    //         {

    //             fprintf(sgrid, "ZONE T = line%d\n", i);
    //             fprintf(sgrid, "NODES = %d, ELEMENTS = %d, ZONETYPE = FELINESEG\n", Line_tmp[i], max(Line_tmp[i] - 1, static_cast<IntType>(1)));
    //             fprintf(sgrid, "DATAPACKING = POINT\n");

    //             for (IntType j = 0; j < Line_tmp[i]; ++j)
    //             {
    //                 IntType cell_tmp = CellOfLine_tmp[count + j];

    //                 fprintf(sgrid, "%20.16e %20.16e %20.16e\n", xcc[cell_tmp], ycc[cell_tmp], zcc[cell_tmp]);
    //             }

    //             if (Line_tmp[i] == 1)
    //             {
    //                 fprintf(sgrid, "% 4d % 4d\n", 1, 1);
    //             }
    //             else
    //             {
    //                 for (IntType j = 0; j < Line_tmp[i] - 1; ++j)
    //                 {
    //                     fprintf(sgrid, "% 4d % 4d\n", j + 1, j + 2);
    //                 }
    //             }

    //             count += Line_tmp[i];
    //         }

    //         fclose(sgrid);
    //     } while (0);

    IntType is_joint_lines = 1;
    GetData(&is_joint_lines, INT, 1, "is_joint_lines", 0);

    if (is_search_lines_faraway)
        if (is_joint_lines)
        { ////joint line nLinesInWall //物面线仅仅往后连
            // Get flow variables
            RealFlow *q[5];
            q[0] = (RealFlow *)GetDataPtr(REAL_FLOW, n, "rho");
            q[1] = (RealFlow *)GetDataPtr(REAL_FLOW, n, "u");
            q[2] = (RealFlow *)GetDataPtr(REAL_FLOW, n, "v");
            q[3] = (RealFlow *)GetDataPtr(REAL_FLOW, n, "w");
            // q[4] = (RealFlow *)GetDataPtr(REAL_FLOW, n, "p");

            IntType *old_elem_start_idx = new IntType[iLine];
            bool *old_line_mark = new bool[iLine];

            old_elem_start_idx[0] = 0;
            for (IntType i = 1; i < iLine; ++i)
            {
                old_elem_start_idx[i] = old_elem_start_idx[i - 1] + Line_tmp[i - 1];
            }

            for (IntType i = 0; i < iLine; ++i)
            {
                old_line_mark[i] = false;
            }

            vector<IntType> Line_tmp_old = Line_tmp;
            vector<IntType> CellOfLine_tmp_old = CellOfLine_tmp;
            vector<IntType> FaceOfLine_tmp_old = FaceOfLine_tmp;
            IntType *LineOfCell_old = LineOfCell;

            Line_tmp.resize(0);
            CellOfLine_tmp.resize(0);
            FaceOfLine_tmp.resize(0);
            LineOfCell = NULL;
            snew_array_1D(LineOfCell, n);
            for (IntType i_cell = 0; i_cell < nTCell; ++i_cell)
            { // out of boundary condition //此时实际的线的数目可能< iLine ，但不影响
                LineOfCell[i_cell] = -1;
            }

            IntType count = 0;
            for (IntType i = 0; i < iLine; ++i)
            {
                if (old_line_mark[i])
                { // i marked 已经被别人合并
                    continue;
                }

                /// old_line_mark[i] = true;

                vector<IntType> cells_of_this_line;
                vector<IntType> faces_of_this_line;

                for (IntType j = 0; j < Line_tmp_old[i]; ++j)
                {
                    cells_of_this_line.push_back(CellOfLine_tmp_old[old_elem_start_idx[i] + j]);
                    faces_of_this_line.push_back(FaceOfLine_tmp_old[old_elem_start_idx[i] + j]);
                }

                IntType end_cell = CellOfLine_tmp_old[old_elem_start_idx[i] + Line_tmp_old[i] - 1];
                IntType line_end_cell = LineOfCell_old[end_cell];

                while (true)
                {
                    vector<xcc_idx> cells_may_be_jointed;

                    for (IntType j = 0; j < nFPC[end_cell]; ++j)
                    {
                        IntType face2 = C2F[end_cell][j];
                        IntType cell2 = f2c[face2 + face2] + f2c[face2 + face2 + 1] - end_cell;

                        if (cell2 >= nTCell)
                        { // ghost cell
                            continue;
                        }

                        IntType line2 = LineOfCell_old[cell2];

                        if (line2 == line_end_cell)
                        {
                            continue;
                        }

                        if ((cell2 != CellOfLine_tmp_old[old_elem_start_idx[line2]]) || (old_line_mark[line2]))
                        { // cell2 不是line2的第一个单元 或者 line2 marked 已经被别人合并
                            continue;
                        }

                        RealFlow ru = rho00 * u00, rv = rho00 * v00, rw = rho00 * w00;
                        if (q[0] != NULL)
                        {
                            ru = 0.5 * (q[0][end_cell] * q[1][end_cell] + q[0][cell2] * q[1][cell2]);
                            rv = 0.5 * (q[0][end_cell] * q[2][end_cell] + q[0][cell2] * q[2][cell2]);
                            rw = 0.5 * (q[0][end_cell] * q[3][end_cell] + q[0][cell2] * q[3][cell2]);
                        }

                        RealFlow tmp = xfn[face2] * ru + yfn[face2] * rv + zfn[face2] * rw;

                        if (f2c[face2 + face2] != end_cell)
                        {
                            tmp = -tmp;
                        }

                        if (tmp / norm_of_rq > 1.0e-14)
                        {
                            xcc_idx tmp2;
                            tmp2.val = area[face2] * tmp * tmp;
                            tmp2.idx = face2;
                            cells_may_be_jointed.push_back(tmp2);
                        }
                    } // end of face of end_cell

                    if (cells_may_be_jointed.size() == 0)
                    { // not found

                        // old_line_mark[line_end_cell] = true;

                        break;
                    }
                    else
                    {
                        {
                            heap_inplace<xcc_idx> tmp3(&cells_may_be_jointed[0], cells_may_be_jointed.size());
                            tmp3.sort();
                        }

                        IntType face2 = cells_may_be_jointed[0].idx;
                        IntType cell2 = f2c[face2 + face2] + f2c[face2 + face2 + 1] - end_cell; // 找到一个可能连接的单元，但是还要判断 end_cell 是 cell2 的最大。
                        IntType line2 = LineOfCell_old[cell2];

                        //////////////////////////////////////////////////////////////////////////
                        vector<xcc_idx> cells_may_be_jointed_cell2;

                        for (IntType j = 0; j < nFPC[cell2]; ++j)
                        {
                            IntType face3 = C2F[cell2][j];
                            IntType cell3 = f2c[face3 + face3] + f2c[face3 + face3 + 1] - cell2;

                            if (cell3 >= nTCell)
                            { // ghost cell
                                continue;
                            }

                            IntType line3 = LineOfCell_old[cell3];

                            if ((cell3 != CellOfLine_tmp_old[old_elem_start_idx[line3] + Line_tmp_old[line3] - 1]) || (old_line_mark[line3]))
                            { // cell3 不是line3的最后一个单元 或者 line3 marked 已经被别人合并
                                continue;
                            }

                            RealFlow ru = rho00 * u00, rv = rho00 * v00, rw = rho00 * w00;
                            if (q[0] != NULL)
                            {
                                ru = 0.5 * (q[0][cell3] * q[1][cell3] + q[0][cell2] * q[1][cell2]);
                                rv = 0.5 * (q[0][cell3] * q[2][cell3] + q[0][cell2] * q[2][cell2]);
                                rw = 0.5 * (q[0][cell3] * q[3][cell3] + q[0][cell2] * q[3][cell2]);
                            }

                            RealFlow tmp = xfn[face3] * ru + yfn[face3] * rv + zfn[face3] * rw;

                            if (f2c[face3 + face3] != cell2)
                            {
                                tmp = -tmp;
                            }

                            if (tmp / norm_of_rq < -1.0e-14)
                            {
                                xcc_idx tmp2;
                                tmp2.val = -area[face3] * tmp * tmp;
                                tmp2.idx = face3;
                                cells_may_be_jointed_cell2.push_back(tmp2);
                            }
                        } // end of face of cell

                        if (cells_may_be_jointed_cell2.size() == 0)
                        { // not found //应该一定会出现！！！
                            // ok

                            break;
                        }
                        else
                        {
                            {
                                heap_inplace<xcc_idx> tmp3(&cells_may_be_jointed_cell2[0], cells_may_be_jointed_cell2.size());
                                tmp3.sort();
                            }

                            IntType face3 = cells_may_be_jointed_cell2[0].idx;
                            IntType cell3 = f2c[face3 + face3] + f2c[face3 + face3 + 1] - cell2; // 找到一个可能连接的单元，但是还要判断 end_cell 是 cell2 的最大。
                            // IntType line3 = LineOfCell_old[cell3];

                            if (cell3 != end_cell)
                            { // 不可连接
                                break;
                            }
                            else // 可连接
                            {
                                if (line2 < nLinesInWall) // 不可能？？
                                {
                                    /// int illal = 1;
                                    break;
                                }

                                for (IntType j = 0; j < Line_tmp_old[line2]; ++j)
                                {
                                    cells_of_this_line.push_back(CellOfLine_tmp_old[old_elem_start_idx[line2] + j]);
                                }

                                faces_of_this_line.push_back(face2);
                                for (IntType j = 1; j < Line_tmp_old[line2]; ++j)
                                {
                                    faces_of_this_line.push_back(FaceOfLine_tmp_old[old_elem_start_idx[line2] + j]);
                                }

                                old_line_mark[line_end_cell] = true;
                                end_cell = CellOfLine_tmp_old[old_elem_start_idx[line2] + Line_tmp_old[line2] - 1];
                                line_end_cell = line2;
                            }
                        }
                    } // end of end_cell's face
                }     // end of while

                Line_tmp.push_back(cells_of_this_line.size());

                for (size_t j = 0; j < cells_of_this_line.size(); ++j)
                {
                    CellOfLine_tmp.push_back(cells_of_this_line[j]);
                    FaceOfLine_tmp.push_back(faces_of_this_line[j]);
                    LineOfCell[cells_of_this_line[j]] = Line_tmp.size() - 1;
                }
                old_line_mark[line_end_cell] = true;
            } // end all old line

            sdel_array_1D(LineOfCell_old);
            delete[] old_line_mark;
            delete[] old_elem_start_idx;
        } // joint lines

    /// Wisces: 先注释掉
    // if (is_debug_for_lineimplicit)
    //     do
    //     { // 按照线长 排序调试
    //         if (mflog::log.rank_id() == 1)
    //         {
    //             xcc_idx *data = new xcc_idx[Line_tmp.size()];

    //             for (IntType i = 0; i < Line_tmp.size(); ++i)
    //             {
    //                 data[i].val = Line_tmp[i];
    //                 data[i].idx = i;
    //             }

    //             heap_inplace<xcc_idx> tmp(data, Line_tmp.size());
    //             tmp.sort();

    //             for (IntType i = 0; i < Line_tmp.size(); ++i)
    //             {
    //                 if (fabs(data[i].val - 1.0) < 0.001)
    //                 {
    //                     cout << "\nline length 1 * " << Line_tmp.size() - i << endl;
    //                     break;
    //                 }

    //                 cout << (IntType)data[i].val << " ";
    //                 if ((i % 16) == 15)
    //                     cout << "\n";
    //             }
    //             cout << "\n";
    //             delete[] data;
    //         }

    //         std::string file_name = "line_info_zone_after_jointed";
    //         // #ifdef MF_MPICH
    //         //             /// file_name += std::to_string( mflog::log.rank_id() );
    //         //             file_name += int2str(mflog::log.rank_id());
    //         // #endif
    //         file_name += ".dat";

    //         FILE *sgrid = fopen(file_name.c_str(), "w");
    //         // Open for read (will fail if file does not exist)
    //         if (sgrid == 0)
    //         {
    //             printf("The file nodetopo.dat was not opened\n");
    //             break;
    //         }

    //         fprintf(sgrid, "TITLE = \"3D Unstructured Hybrid Mesh Line Info\"\n");
    //         fprintf(sgrid, "FILETYPE = FULL\n");
    //         fprintf(sgrid, "VARIABLES = x, y, z\n");

    //         IntType count = 0;
    //         for (IntType i = 0; i < Line_tmp.size(); ++i)
    //         {

    //             fprintf(sgrid, "ZONE T = line%d\n", i);
    //             fprintf(sgrid, "NODES = %d, ELEMENTS = %d, ZONETYPE = FELINESEG\n", Line_tmp[i], max(Line_tmp[i] - 1, static_cast<IntType>(1)));
    //             fprintf(sgrid, "DATAPACKING = POINT\n");

    //             for (IntType j = 0; j < Line_tmp[i]; ++j)
    //             {
    //                 IntType cell_tmp = CellOfLine_tmp[count + j];

    //                 fprintf(sgrid, "%20.16e %20.16e %20.16e\n", xcc[cell_tmp], ycc[cell_tmp], zcc[cell_tmp]);
    //             }

    //             if (Line_tmp[i] == 1)
    //             {
    //                 fprintf(sgrid, "% 4d % 4d\n", 1, 1);
    //             }
    //             else
    //             {
    //                 for (IntType j = 0; j < Line_tmp[i] - 1; ++j)
    //                 {
    //                     fprintf(sgrid, "% 4d % 4d\n", j + 1, j + 2);
    //                 }
    //             }

    //             count += Line_tmp[i];
    //         }

    //         fclose(sgrid);
    //     } while (0);

    // #ifdef MF_MPICH
    //     {
    //         IntType n_line = Line_tmp.size();
    //         IntType n_line_glb;
    //         MPI_Allreduce(&n_line, &n_line_glb, 1, MPIIntType, MPI_SUM, MPI_COMM_WORLD);
    //         IntType n_cell = CellOfLine_tmp.size();
    //         IntType n_cell_glb;
    //         MPI_Allreduce(&n_cell, &n_cell_glb, 1, MPIIntType, MPI_SUM, MPI_COMM_WORLD);
    //         if (n_line_glb)
    //         {
    //             mflog::log.set_one_processor_out();
    //             mflog::log << endl
    //                        << "number of final Lines " << n_line_glb << ", average length " << n_cell_glb / n_line_glb << endl;
    //         }
    //     }
    // #else
    if (Line_tmp.size())
    {
        cout << "number of final Lines " << Line_tmp.size() << ", average length " << CellOfLine_tmp.size() / Line_tmp.size() << endl;
    }
    // #endif

    //////////////////////////////////////////////////////////////////////////
    for (IntType i_cell = nTCell; i_cell < n; ++i_cell)
    {                                         // out of boundary condition //此时实际的线的数目可能< iLine ，但不影响
        LineOfCell[i_cell] = Line_tmp.size(); // iLine
    }

    IntType nLine = (IntType)Line_tmp.size();
    IntType *nCellsInLines = NULL;
    snew_array_1D(nCellsInLines, nLine);
    for (IntType i = 0; i < nLine; ++i)
    {
        nCellsInLines[i] = Line_tmp[i];
    }

    IntType maxLine = 0;
    for (IntType i = 0; i < nLine; ++i)
    {
        maxLine = max(maxLine, nCellsInLines[i]);
    }

    // #ifdef MF_MPICH
    //     IntType maxLine_glb;
    //     MPI_Allreduce(&maxLine, &maxLine_glb, 1, MPIIntType, MPI_MAX, MPI_COMM_WORLD);
    //     maxLine = maxLine_glb;
    // #endif
    // mflog::log.set_one_processor_out();
    // mflog::log << endl
    //            << "line implicit, max cell number in a line is: " << maxLine << endl;
    cout << endl
         << "line implicit, max cell number in a line is: " << maxLine << endl;

    IntType **CellOfLine = NULL;
    snew_array_2D(CellOfLine, nLine, nCellsInLines, true);
    IntType count = 0;

    IntType **FaceOfLine = NULL;
    snew_array_2D(FaceOfLine, nLine, nCellsInLines, true);
    count = 0;
    for (IntType i = 0; i < nLine; ++i)
    {
        for (IntType j = 0; j < nCellsInLines[i]; j++)
        {
            CellOfLine[i][j] = CellOfLine_tmp[count];
            FaceOfLine[i][j] = FaceOfLine_tmp[count];
            ++count;
        }
    }

    SetLI_nLine(nLine);
    SetLI_nCellsInLine(nCellsInLines);
    SetLI_CellOfLine(CellOfLine);
    SetLI_FaceOfLine(FaceOfLine);
    SetLI_LineOfCell(LineOfCell);

    // sdel_array_1D(mark_cell);

    // order 0 in order above
    IntType *Layer = NULL;
    IntType *LUSGSCellOrder = NULL;
    snew_array_1D(Layer, n);
    snew_array_1D(LUSGSCellOrder, nTCell);

    for (IntType i = 0; i < nTCell; i++)
    {
        LUSGSCellOrder[i] = CellOfLine_tmp[i];
        Layer[LUSGSCellOrder[i]] = i;
    }

    for (IntType i = nTCell; i < n; i++)
    {
        Layer[i] = nTCell;
    }

    this->UpdateDataPtr(LUSGSCellOrder, INT, nTCell, "LUSGSCellOrder");
    this->UpdateDataPtr(Layer, INT, n, "LUSGSLayer");

    /// Wisces: 先注释掉
    // // tecplot file of line info
    // if (is_debug_for_lineimplicit >= 2)
    // {
    //     FILE *sgrid = fopen("line_info.dat", "w");
    //     // Open for read (will fail if file does not exist)
    //     if (sgrid == 0)
    //     {
    //         printf("The file nodetopo.dat was not opened\n");
    //         return;
    //     }

    //     for (IntType i = 0; i < nLine; ++i)
    //     {
    //         int ce22 = CellOfLine[i][0];
    //         // cout << ce22 << " " << xcc[ce22] << " "<< ycc[ce22] << " "<< zcc[ce22] << endl;
    //         fprintf(sgrid, "% d %20.16e %20.16e %20.16e\n", ce22, xcc[ce22], ycc[ce22], zcc[ce22]);
    //     }

    //     fprintf(sgrid, "TITLE = \"3D Unstructured Hybrid Mesh Line Info\"\n");
    //     fprintf(sgrid, "FILETYPE = FULL\n");
    //     fprintf(sgrid, "VARIABLES = x, y, z, line_idx\n");
    //     fprintf(sgrid, "ZONE T = \"Topo Zone\"\n");
    //     fprintf(sgrid, "NODES = %d, ELEMENTS = %d, ZONETYPE = FEBRICK\n", nTNode, nTCell);
    //     fprintf(sgrid, "DATAPACKING = BLOCK\n");
    //     fprintf(sgrid, "VARLOCATION = ([1-3]=NODAL,[4]=CELLCENTERED)\n");

    //     RealGeom *xxx = GetX();
    //     RealGeom *yyy = GetY();
    //     RealGeom *zzz = GetZ();

    //     for (IntType n = 0; n < nTNode; n++)
    //     {
    //         fprintf(sgrid, "%f ", xxx[n]);
    //         if ((n % 16) == 15)
    //             fprintf(sgrid, "\n");
    //     }
    //     fprintf(sgrid, "\n");
    //     for (IntType n = 0; n < nTNode; n++)
    //     {
    //         fprintf(sgrid, "%f ", yyy[n]);
    //         if ((n % 16) == 15)
    //             fprintf(sgrid, "\n");
    //     }
    //     fprintf(sgrid, "\n");
    //     for (IntType n = 0; n < nTNode; n++)
    //     {
    //         fprintf(sgrid, "%f ", zzz[n]);
    //         if ((n % 16) == 15)
    //             fprintf(sgrid, "\n");
    //     }
    //     fprintf(sgrid, "\n");

    //     for (IntType e = 0; e < nTCell; e++)
    //     {
    //         fprintf(sgrid, "%d ", LineOfCell[e]);
    //         if ((e % 32) == 31)
    //             fprintf(sgrid, "\n");
    //     }
    //     fprintf(sgrid, "\n");

    //     IntType *nNPC = CalnNPC(this);
    //     IntType **C2N = CalC2N(this);

    //     IntType node_max_idx[8];
    //     for (IntType e = 0; e < nTCell; e++)
    //     {
    //         for (IntType j = 0; j < nNPC[e]; j++)
    //         {
    //             node_max_idx[j] = C2N[e][j];
    //         }

    //         for (IntType j = nNPC[e]; j < 8; j++)
    //         {
    //             node_max_idx[j] = C2N[e][nNPC[e] - 1];
    //         }
    //         if (nNPC[e] == 4)
    //         { // tetra
    //             fprintf(sgrid, "% 4d % 4d % 4d % 4d % 4d % 4d % 4d % 4d\n", node_max_idx[0] + 1, node_max_idx[1] + 1, node_max_idx[2] + 1,
    //                     node_max_idx[2] + 1, node_max_idx[3] + 1, node_max_idx[3] + 1, node_max_idx[3] + 1, node_max_idx[3] + 1);
    //         }
    //         else if (nNPC[e] == 5)
    //         { // pyramid
    //             fprintf(sgrid, "% 4d % 4d % 4d % 4d % 4d % 4d % 4d % 4d\n", node_max_idx[0] + 1, node_max_idx[1] + 1, node_max_idx[2] + 1,
    //                     node_max_idx[3] + 1, node_max_idx[4] + 1, node_max_idx[4] + 1, node_max_idx[4] + 1, node_max_idx[4] + 1);
    //         }
    //         else if (nNPC[e] == 6)
    //         { // prism
    //             fprintf(sgrid, "% 4d % 4d % 4d % 4d % 4d % 4d % 4d % 4d\n", node_max_idx[0] + 1, node_max_idx[1] + 1, node_max_idx[2] + 1,
    //                     node_max_idx[2] + 1, node_max_idx[3] + 1, node_max_idx[4] + 1, node_max_idx[5] + 1, node_max_idx[5] + 1);
    //         }
    //         else if (nNPC[e] == 8)
    //         { // hexa
    //             fprintf(sgrid, "% 4d % 4d % 4d % 4d % 4d % 4d % 4d % 4d\n", node_max_idx[0] + 1, node_max_idx[1] + 1, node_max_idx[2] + 1,
    //                     node_max_idx[3] + 1, node_max_idx[4] + 1, node_max_idx[5] + 1, node_max_idx[6] + 1, node_max_idx[7] + 1);
    //         }
    //     }
    //     fclose(sgrid);
    // }
}
*/

/**
//!< Communicate data from the current grid to neighbor grid in zone
void PolyGrid::CommInterfaceData(IntType nbZone, PolyGrid *grid, const char *name)
{
    IntType i, n, c1, c2;
    RealFlow *cq, *nq;
    BCRecord **bcr = Getbcr();
    IntType *f2c_n = grid->Getf2c();
    IntType nBFace = GetNBFace();
    IntType nIFace = 0;

    n = GetNTCell() + nBFace;
    // cq = (RealFlow *)GetDataPtr(REAL_FLOW, n, name);
    if (!cq)
    {
        printf("Variable %s to be communicated not found\n", name);
        return;
    }

    n = grid->GetNTCell() + grid->GetNBFace();
    // nq = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, name);
    if (!nq)
    {
        printf("Variable %s to be communicated not found\n", name);
        return;
    }

    n = 0;
    IntType ncnb = grid->GetNTCell();
    for (i = 0; i < nBFace; i++)
    {
        if (bcr[i]->GetType() == INTERFACE)
        {
            if (nbZ[nIFace] == nbZone)
            {
                c1 = f2c[i * 2];
                c2 = nbBF[nIFace] + ncnb;
                nq[c2] = cq[c1];

                c2 = f2c[i * 2 + 1];
                c1 = f2c_n[nbBF[nIFace] * 2];
                cq[c2] = nq[c1];

                n++;
            }
            nIFace++;
        }
    }
    assert(nIFace == GetNIFace());
}

//!< Communicate data from the current grid to neighbor grid in zone
void PolyGrid::CommCellCenterData(IntType nbZone, PolyGrid *grid)
{
    IntType i, n, c1, c2;
    BCRecord **bcr = Getbcr();
    IntType *f2c_n = grid->Getf2c();
    IntType nBFace = GetNBFace();
    IntType nIFace = 0;

    RealGeom *nxcc = grid->GetXcc();
    RealGeom *nycc = grid->GetYcc();
    RealGeom *nzcc = grid->GetZcc();

    n = 0;
    IntType ncnb = grid->GetNTCell();

    for (i = 0; i < nBFace; i++)
    {
        if (bcr[i]->GetType() == INTERFACE)
        {
            if (nbZ[nIFace] == nbZone)
            {
                c1 = f2c[i * 2];
                c2 = nbBF[nIFace] + ncnb;
                nxcc[c2] = xcc[c1];
                nycc[c2] = ycc[c1];
                nzcc[c2] = zcc[c1];

                c2 = f2c[i * 2 + 1];
                c1 = f2c_n[nbBF[nIFace] * 2];
                xcc[c2] = nxcc[c1];
                ycc[c2] = nycc[c1];
                zcc[c2] = nzcc[c1];

                n++;
            }
            nIFace++;
        }
    }
    assert(nIFace == GetNIFace());
}

//!< Purpose: To compute the least distance from point to box
RealGeom FindRminbox(RealGeom xp, RealGeom yp, RealGeom zp, RealGeom x1, RealGeom x2,
                     RealGeom y1, RealGeom y2, RealGeom z1, RealGeom z2)
{
    RealGeom rr, rx, ry, rz;
    if (xp >= x1 && xp <= x2)
        rx = 0;
    else
        rx = (xp < x1) ? x1 - xp : xp - x2;
    if (yp >= y1 && yp <= y2)
        ry = 0;
    else
        ry = (yp < y1) ? y1 - yp : yp - y2;
    if (zp >= z1 && zp <= z2)
        rz = 0;
    else
        rz = (zp < z1) ? z1 - zp : zp - z2;

    rr = rx * rx + ry * ry + rz * rz;

    return (rr);
}

//!< Purpose: To compute the distance from point to point
RealGeom FindRp2p(RealGeom x1, RealGeom y1, RealGeom z1, RealGeom x2, RealGeom y2, RealGeom z2)
{
    RealGeom dx, dy, dz;
    dx = x1 - x2;
    dy = y1 - y2;
    dz = z1 - z2;
    return (dx * dx + dy * dy + dz * dz);
}

///!<urpose:  To find the closest distance from a field point to the actual surface (i.e. not simply the closest discrete surface point), using local triangulation of the surface.
void FindRp2tri(RealGeom &dist, RealGeom xp, RealGeom yp, RealGeom zp,
                RealGeom xa, RealGeom ya, RealGeom za, RealGeom xb, RealGeom yb, RealGeom zb,
                RealGeom xc, RealGeom yc, RealGeom zc)
{
    RealGeom pp[3], aa[3], bb[3], cc[3];
    pp[0] = xp;
    pp[1] = yp;
    pp[2] = zp;
    aa[0] = xa;
    aa[1] = ya;
    aa[2] = za;
    bb[0] = xb;
    bb[1] = yb;
    bb[2] = zb;
    cc[0] = xc;
    cc[1] = yc;
    cc[2] = zc;

    RealGeom p[3], a[3], b[3], r[3], rr;
    RealGeom daa = 0., dbb = 0., dab = 0., den, dap = 0., dbp = 0., s, t, dsq;
    for (IntType i = 0; i < 3; i++)
    {
        p[i] = pp[i] - aa[i];
        a[i] = bb[i] - aa[i];
        b[i] = cc[i] - aa[i];
        daa += a[i] * a[i];
        dbb += b[i] * b[i];
        dab += a[i] * b[i];
    }
    den = dab * dab - daa * dbb;

    /// Wisces: { return (x > -TINY) && (x < TINY); }
    // if (EqualZero(den))
    if ((den > -TINY) && (den < TINY))
        ; // zhyb: 面积不为零则den不为零
    else
    {
        for (IntType i = 0; i < 3; i++)
        {
            dap += a[i] * p[i];
            dbp += b[i] * p[i];
        }
        s = (dab * dbp - dbb * dap) / den;
        t = (dab * dap - daa * dbp) / den;
        if (s < 0. || t < 0. || (t + s) > 1.)
            ; // zhyb: 这三种情况垂足落在三角形外边
        else
        {
            for (IntType i = 0; i < 3; i++)
                r[i] = p[i] - s * a[i] - t * b[i];
            rr = 0.;
            for (IntType i = 0; i < 3; i++)
                rr += r[i] * r[i];
            if (rr < dist)
                dist = rr;
            return;
        }
    }

    dsq = dist;
    // if (EqualZero(daa))
    if ((daa > -TINY) && (daa < TINY))
        ; // zhyb: bb点和aa点重合
    else
    {
        dap = 0.;
        for (IntType i = 0; i < 3; i++)
            dap += a[i] * p[i];
        t = dap / daa;
        if (t < 0. || t > 1.)
            ; // zhyb: 这两种情况下垂足落在线段aabb外
        else
        {
            for (IntType i = 0; i < 3; i++)
                r[i] = p[i] - t * a[i];
            rr = 0.;
            for (IntType i = 0; i < 3; i++)
                rr += r[i] * r[i];
            if (rr < dsq)
                dsq = rr;
        }
    }

    // if (EqualZero(dbb))
    if ((dbb > -TINY) && (dbb < TINY))
        ; // zhyb: cc点和aa点重合
    else
    {
        dbp = 0.;
        for (IntType i = 0; i < 3; i++)
            dbp += b[i] * p[i];
        t = dbp / dbb;
        if (t < 0. || t > 1.)
            ; // zhyb: 这两种情况下垂足落在线段aacc外
        else
        {
            for (IntType i = 0; i < 3; i++)
                r[i] = p[i] - t * b[i];
            rr = 0.;
            for (IntType i = 0; i < 3; i++)
                rr += r[i] * r[i];
            if (rr < dsq)
                dsq = rr;
        }
    }

    daa = 0.;
    for (IntType i = 0; i < 3; i++)
    {
        p[i] = pp[i] - bb[i];
        a[i] = cc[i] - bb[i];
        daa += a[i] * a[i];
    }

    // if (EqualZero(daa))
    if ((daa > -TINY) && (daa < TINY))
        ; // zhyb: cc点和bb点重合
    else
    {
        dap = 0;
        for (IntType i = 0; i < 3; i++)
            dap += a[i] * p[i];
        t = dap / daa;
        if (t < 0. || t > 1.)
            ; // zhyb: 这两种情况下垂足落在线段bbcc外
        else
        {
            for (IntType i = 0; i < 3; i++)
                r[i] = p[i] - t * a[i];
            rr = 0.;
            for (IntType i = 0; i < 3; i++)
                rr += r[i] * r[i];
            if (rr < dsq)
                dsq = rr;
        }
    }

    if (dsq < dist)
        dist = dsq;
}

//!< Purpose:  To interchange v[i] and v[j]
void xswap(IntType *indx, RealGeom *v, IntType i, IntType j)
{
    RealGeom temp;
    temp = v[i];
    v[i] = v[j];
    v[j] = temp;
    IntType itmp;
    itmp = indx[i];
    indx[i] = indx[j];
    indx[j] = itmp;
}

//!< Purpose:  To sort a list of points
void quicksort(IntType s, IntType e, IntType *indx, RealGeom *v)
{
    IntType i, last;

    if (e - s <= 0) // nothing to do
        return;

    last = s - 1;
    for (i = s; i < e; i++)
    {
        if (v[i] < v[s - 1])
            xswap(indx, v, ++last, i);
    }

    xswap(indx, v, s - 1, last);
    quicksort(s, last, indx, v);
    quicksort(last + 2, e, indx, v);
}

//   Compute the distance to wall(triangles) about the cell center
//   Author: zm 20080810
//   Modify: zm 20150423
//           zhyb 20180914 考虑顶点全在物面上，而所有面均不是物面的单元
void PolyGrid::ComputeDist2WallTriang(RealGeom *dist2wall_cell, IntType mark)
{
    IntType i, j, k, count;
    FILE *fp;
    IntType nTNode = GetNTNode();
    RealGeom *x = GetX(), *y = GetY(), *z = GetZ();
    String filename, filename_tmp;

    // mflog::log.set_one_processor_out();
    // mflog::log << "Now compute the distance to wall!" << endl;
    cout << "Now compute the distance to wall!" << endl;

    // GetData(filename, STRING, 1, "griddir");
    sprintf(filename_tmp, "GeominfoDist.dat");
    strcat(filename, filename_tmp);

    // mflog::log.set_each_grid_out();
    fp = fopen(filename, "rb");
    if (!fp)
    {
        // mflog::log << "Failed to open file " << filename << " in function ComputeDist2WallTriang!" << endl;
        // mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
        cout << "Failed to open file " << filename << " in function ComputeDist2WallTriang!" << endl;
        exit(1);
    }
    IntType nSP, nSF, *nSfP, *SfP, *nPntS, *PntS, *SP;
    IntType nBox, *nPBox;
    RealGeom **BBox;

    fread(&nSP, sizeof(IntType), 1, fp);

    RealGeom *xSf = NULL; // the coordinate
    RealGeom *ySf = NULL;
    RealGeom *zSf = NULL;
    snew_array_1D(xSf, nSP);
    snew_array_1D(ySf, nSP);
    snew_array_1D(zSf, nSP);
    fread(xSf, sizeof(RealGeom), nSP, fp);
    fread(ySf, sizeof(RealGeom), nSP, fp);
    fread(zSf, sizeof(RealGeom), nSP, fp);

    fread(&nSF, sizeof(IntType), 1, fp); // Number of the Solid-Face

    nSfP = NULL;
    snew_array_1D(nSfP, nSP + 1);
    fread(nSfP, sizeof(IntType), nSP + 1, fp); // Number of the Solid-Face that connected Local Point
    SfP = NULL;
    snew_array_1D(SfP, nSfP[nSP]);
    fread(SfP, sizeof(IntType), nSfP[nSP], fp); // these Solid-Faces that connected Local Point

    nPntS = NULL;
    snew_array_1D(nPntS, nSF + 1);
    fread(nPntS, sizeof(IntType), nSF + 1, fp);
    PntS = NULL;
    snew_array_1D(PntS, nPntS[nSF]);
    fread(PntS, sizeof(IntType), nPntS[nSF], fp);

    SP = NULL;
    snew_array_1D(SP, nSP);
    fread(SP, sizeof(IntType), nSP, fp); // Number of the boxs

    fread(&nBox, sizeof(IntType), 1, fp);
    nPBox = NULL; // Number of the Points in the box
    snew_array_1D(nPBox, nBox + 1);
    fread(nPBox, sizeof(IntType), nBox + 1, fp);

    BBox = NULL; // Six Bound-Core of the Box
    snew_array_2D(BBox, nBox, 6, true);
    for (i = 0; i < nBox; i++)
        fread(BBox[i], sizeof(RealGeom), 6, fp);

    fclose(fp);

    // mflog::log.set_one_processor_out();
    // mflog::log << "Finished to Read the coordinate of all nodes on the solid surface" << endl;
    cout << "Finished to Read the coordinate of all nodes on the solid surface" << endl;

    RealGeom *distP = NULL;
    snew_array_1D(distP, nTNode);
    for (i = 0; i < nTNode; i++)
        distP[i] = BIG;

    for (i = 0; i < nTNode; i++)
    {
        IntType pntmin = -1;

        RealGeom *distB = NULL;
        snew_array_1D(distB, nBox);
        for (j = 0; j < nBox; j++)
        {
            distB[j] = FindRminbox(x[i], y[i], z[i], BBox[j][0], BBox[j][1], BBox[j][2],
                                   BBox[j][3], BBox[j][4], BBox[j][5]);
        }

        IntType *Bsort = NULL;
        snew_array_1D(Bsort, nBox);
        for (j = 0; j < nBox; j++)
            Bsort[j] = j;
        quicksort(1, nBox, Bsort, distB);

        RealFlow dd, distP2P = BIG;
        IntType B, pnt2;
        for (j = 0; j < nBox; j++)
        {
            if (distP[i] < distB[j])
                break;

            B = Bsort[j];
            for (k = nPBox[B]; k < nPBox[B + 1]; k++)
            {
                pnt2 = SP[k];
                dd = FindRp2p(x[i], y[i], z[i], xSf[pnt2], ySf[pnt2], zSf[pnt2]);
                if (distP2P > dd)
                {
                    distP2P = dd;
                    pntmin = pnt2;
                }
            }
        }

        for (k = nSfP[pntmin]; k < nSfP[pntmin + 1]; k++)
        {
            IntType sface = SfP[k];
            IntType pnt[4];
            for (IntType jj = nPntS[sface]; jj < nPntS[sface + 1]; jj++)
                pnt[jj - nPntS[sface]] = PntS[jj];
            if (nPntS[sface + 1] - nPntS[sface] == 4)
            {
                if (pntmin == pnt[0])
                    pnt[2] = pnt[3];
                else if (pntmin == pnt[2])
                    pnt[0] = pnt[3];
                else if (pntmin == pnt[3])
                    pnt[1] = pnt[3];
            }
            FindRp2tri(distP2P, x[i], y[i], z[i], xSf[pnt[0]], ySf[pnt[0]], zSf[pnt[0]],
                       xSf[pnt[1]], ySf[pnt[1]], zSf[pnt[1]], xSf[pnt[2]], ySf[pnt[2]], zSf[pnt[2]]);
        }
        if (distP2P < distP[i])
            distP[i] = distP2P;
        distP[i] = sqrt(distP[i]);
        sdel_array_1D(distB);
        sdel_array_1D(Bsort);
    }
    sdel_array_2D(BBox);
    sdel_array_1D(nPBox);

    RealGeom *distC = NULL;
    snew_array_1D(distC, nTCell);
    IntType *nNPC = CalnNPC(this);
    IntType **C2N = CalC2N(this);

    for (i = 0; i < nTCell; i++)
    {
        distC[i] = 0.;
        for (j = 0; j < nNPC[i]; j++)
        {
            distC[i] += distP[C2N[i][j]];
        }
        distC[i] /= nNPC[i];
    }

    IntType type, c1, c2;
    RealGeom d;
    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        if (type == WALL)
        {
            c1 = f2c[i + i];
            distC[c1] = BIG;
        }
    }
    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        if (type == WALL)
        {
            c1 = f2c[i + i];
            c2 = f2c[i + i + 1];

            d = FindRp2p(xcc[c2], ycc[c2], zcc[c2], xcc[c1], ycc[c1], zcc[c1]);
            d = 0.5 * sqrt(d);
            distC[c1] = MIN(distC[c1], d);
        }
    }

    // 考虑顶点全在物面上，而所有面均不是物面的单元
    IntType *mark_dist0 = NULL;
    snew_array_1D(mark_dist0, nTCell);
    count = 0;
    for (i = 0; i < nTCell; i++)
    {
        if (distC[i] < TINY)
        {
            mark_dist0[i] = 1;
            count++;
        }
        else
        {
            mark_dist0[i] = 0;
        }
    }

    // #ifdef MF_MPICH
    //     Parallel::parallel_sum(count, GridComm);
    // #endif

    if (count)
    { // 可能存在
        // mflog::log.set_one_processor_out();
        // mflog::log << endl
        //            << "There is " << count << " cell's distC=0! Now correcting!" << endl;
        cout << endl
             << "There is " << count << " cell's distC=0! Now correcting!" << endl;

        IntType p1, f1;
        IntType *nNPF = GetnNPF();
        IntType **F2N = CalF2N(this);
        IntType *nFPC = CalnFPC(this);
        IntType **C2F = CalC2F(this);

        IntType *mark_wallnode = NULL; // 物面点
        snew_array_1D(mark_wallnode, nTNode);
        for (i = 0; i < nTNode; i++)
        {
            mark_wallnode[i] = 0;
        }
        for (i = 0; i < nBFace; i++)
        {
            type = bcr[i]->GetType();
            if (type != WALL)
                continue;

            for (j = 0; j < nNPF[i]; j++)
            {
                p1 = F2N[i][j];
                mark_wallnode[p1] = 1;
            }
        }
        // #ifdef MF_MPICH
        //         CommInternodeDataMPISUM(mark_wallnode);
        // #endif

        for (i = 0; i < nTCell; i++)
        {
            if (!mark_dist0[i])
                continue; // 不是距离为0单元，跳出

            mark = 0;
            for (j = 0; j < nNPC[i]; j++)
            {
                p1 = C2N[i][j];
                if (!mark_wallnode[p1])
                { // 这个点不是物面点
                    mark = 1;
                }
            }
            if (mark)
                continue; // 不是所有点均在物面上，跳出

            mark = 0;
            for (j = 0; j < nFPC[i]; j++)
            {
                f1 = C2F[i][j];
                if (f1 < nBFace)
                {
                    type = bcr[f1]->GetType();
                    if (type == WALL)
                    {
                        mark = 1;
                    }
                }
            }
            if (mark)
                continue; // 有物面边界，跳出

            // 求体心到单元顶点的最小距离，作为其到物面的距离
            distC[i] = BIG;
            for (j = 0; j < nNPC[i]; j++)
            {
                p1 = C2N[i][j];

                d = FindRp2p(xcc[i], ycc[i], zcc[i], x[p1], y[p1], z[p1]);
                distC[i] = MIN(distC[i], d);
            }
            distC[i] = sqrt(distC[i]);
        }
        sdel_array_1D(mark_wallnode);
    }
    sdel_array_1D(mark_dist0);

    // mflog::log.set_all_processors_out();

    count = 0;
    for (i = 0; i < nTCell; i++)
    {
        if (distC[i] < TINY)
        {
            // #ifdef MF_MPICH
            //             mflog::log << endl
            //                        << "Zone " << myZone << "  Cell " << i << " dist2wall<TINY" << endl;
            // #else
            // mflog::log << endl
            //            << "Cell " << i << " dist2wall<TINY" << endl;
            cout << endl
                 << "Cell " << i << " dist2wall<TINY" << endl;
            // #endif
            count++;
        }
    }
    // #ifdef MF_MPICH
    //     Parallel::parallel_sum(count, GridComm);
    // #endif
    if (count != 0)
        // mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
        exit(1);

    for (i = 0; i < nTCell; i++)
        dist2wall_cell[i] = MIN(dist2wall_cell[i], distC[i]);
    //  DumpPlyhedra3D_FieldView(this,0);

    sdel_array_1D(distP);
    sdel_array_1D(distC);
    sdel_array_1D(xSf);
    sdel_array_1D(ySf);
    sdel_array_1D(zSf);
    sdel_array_1D(nSfP);
    sdel_array_1D(SfP);
    sdel_array_1D(nPntS);
    sdel_array_1D(PntS);
    sdel_array_1D(SP);

    // mflog::log.set_one_processor_out();
    // mflog::log << "Distance is OK !!" << endl;
    cout << "Distance is OK !!" << endl;
}

//!< Write the information for Computing the distance to wall about the cell center in the MPI.
void PolyGrid::WriteInfoDist()
{
    IntType i, j, k, count, type;
    IntType nTNode = GetNTNode();
    RealGeom *x = GetX(), *y = GetY(), *z = GetZ();
    FILE *fp;
    String filename;

    IntType nSF = 0; // 物面数

    IntType *nFN = NULL;
    ;
    snew_array_1D(nFN, nBFace + 1);

    nFN[0] = 0;
    for (i = 1; i <= nBFace; i++)
        nFN[i] = nFN[i - 1] + nNPF[i - 1];

    IntType *mrk = NULL;
    snew_array_1D(mrk, nTNode);
    for (i = 0; i < nTNode; i++)
        mrk[i] = 0;
    IntType nSP = 0; // nSP 物面上点的数量
    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        if (type == WALL)
        {
            nSF++;
            for (k = nFN[i]; k < nFN[i + 1]; k++)
            {
                if (mrk[f2n[k]] == 0)
                    nSP++;
                mrk[f2n[k]]++;
            }
        }
    }

    IntType *PtNew = NULL;
    RealGeom *xSf = NULL; // 物面上点的坐标
    RealGeom *ySf = NULL;
    RealGeom *zSf = NULL;
    snew_array_1D(PtNew, nTNode);
    snew_array_1D(xSf, nSP);
    snew_array_1D(ySf, nSP);
    snew_array_1D(zSf, nSP);
    count = 0;
    for (i = 0; i < nTNode; i++)
    {
        if (mrk[i] > 0)
        {
            PtNew[i] = count;
            xSf[count] = x[i];
            ySf[count] = y[i];
            zSf[count] = z[i];
            mrk[count] = mrk[i];
            count++;
        }
        else
        {
            PtNew[i] = -1;
        }
    }

    // mflog::log.set_all_processors_out();

    IntType *nPntS = NULL; // zhyb: 每个物面面单元对应的点数，第i个物面面单元对应的点数为nPntS[i+1]-nPntS[i]
    IntType *nSfP = NULL;  // zhyb: 每个物面点对应的面数，第i个物面点对应的面数为nSfP[i+1]-nSfP[i]
    snew_array_1D(nPntS, nSF + 1);
    snew_array_1D(nSfP, nSP + 1);
    nPntS[0] = 0;
    nSfP[0] = 0;
    for (i = 0; i < nSP; i++)
        nSfP[i + 1] = nSfP[i] + mrk[i];
    sdel_array_1D(mrk);
    // mflog::log << std::endl
    //            << " nPntS = " << nSfP[nSP] << std::endl;
    cout << std::endl
         << " nPntS = " << nSfP[nSP] << std::endl;

    IntType *ntmp = NULL;
    IntType *SfP = NULL;  // 物面点的相关面   //zhyb: 每个物面点对应的面号
    IntType *PntS = NULL; // zhyb: 每个物面面对应的点号
    snew_array_1D(ntmp, nSP);
    snew_array_1D(SfP, nSfP[nSP]);
    snew_array_1D(PntS, nSfP[nSP]);

    for (i = 0; i < nSP; i++)
        ntmp[i] = 0;
    nSF = 0;
    count = 0;
    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        if (type == WALL)
        {
            for (k = nFN[i]; k < nFN[i + 1]; k++)
            {
                IntType pnt = PtNew[f2n[k]];
                PntS[count++] = pnt;
                SfP[nSfP[pnt] + ntmp[pnt]] = nSF;
                ntmp[pnt]++;
            }
            nPntS[++nSF] = count;
        }
    }
    sdel_array_1D(ntmp);
    sdel_array_1D(nFN);
    sdel_array_1D(PtNew);
    // mflog::log << endl
    //            << " nPntS = " << nPntS[0] << IOS_SEP << nPntS[nSF] << endl;
    cout << endl
         << " nPntS = " << nPntS[0] << IOS_SEP << nPntS[nSF] << endl;

    IntType *SP = NULL; // 物面点的序号
    RealGeom *xtmp = NULL;
    snew_array_1D(SP, nSP);
    snew_array_1D(xtmp, nSP);
    for (i = 0; i < nSP; i++)
    {
        SP[i] = i;
        xtmp[i] = xSf[i];
    }
    quicksort(1, nSP, SP, xtmp);
    sdel_array_1D(xtmp);

    IntType nBox; // zhyb: 等于物面点数的开方+1
    count = IntType(sqrt(1. * nSP));
    nBox = static_cast<IntType>(nSP / (count + TINY) + 1);

    IntType *nPBox = NULL;
    snew_array_1D(nPBox, nBox + 1);
    for (i = 0; i < nBox; i++)
    {
        nPBox[i] = i * count;
        if (nPBox[i] > nSP)
        {
            nBox = i;
            break;
        }
    }
    nPBox[i] = nSP;
    // mflog::log << "nBox = " << nBox << " nSP = " << nSP << IOS_SEP << count << endl;
    cout << "nBox = " << nBox << " nSP = " << nSP << IOS_SEP << count << endl;

    RealGeom **BBox = NULL; // sort: xmin, xmax, ymin, ymax, zmin, zmax
    snew_array_2D(BBox, nBox, 6, true);
    for (i = 0; i < nBox; i++)
    {
        BBox[i][0] = BIG;
        BBox[i][1] = -BIG;
        BBox[i][2] = BIG;
        BBox[i][3] = -BIG;
        BBox[i][4] = BIG;
        BBox[i][5] = -BIG;
    }
    for (i = 0; i < nBox; i++)
    {
        IntType pnt;
        for (j = nPBox[i]; j < nPBox[i + 1]; j++)
        {
            pnt = SP[j];
            if (xSf[pnt] < BBox[i][0])
                BBox[i][0] = xSf[pnt];
            if (xSf[pnt] > BBox[i][1])
                BBox[i][1] = xSf[pnt];
            if (ySf[pnt] < BBox[i][2])
                BBox[i][2] = ySf[pnt];
            if (ySf[pnt] > BBox[i][3])
                BBox[i][3] = ySf[pnt];
            if (zSf[pnt] < BBox[i][4])
                BBox[i][4] = zSf[pnt];
            if (zSf[pnt] > BBox[i][5])
                BBox[i][5] = zSf[pnt];
        }
    }
    // mflog::log << "Finished to compute BBox" << endl;
    cout << "Finished to compute BBox" << endl;

    String filetmp;
    // GetData(filename, STRING, 1, "griddir");
    sprintf(filetmp, "GeominfoDist.dat");
    strcat(filename, filetmp);

    fp = fopen(filename, "wb");
    // fp=fopen(filetmp,"wb");
    if (!fp)
    {
        // mflog::log << "Failed to open file " << filename << " in function WriteInfoDist!" << endl;
        // mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
        cout << "Failed to open file " << filename << " in function WriteInfoDist!" << endl;
        exit(1);
    }
    // write the coorade of all nodes
    fwrite(&nSP, sizeof(IntType), 1, fp);
    fwrite(xSf, sizeof(RealGeom), nSP, fp);
    fwrite(ySf, sizeof(RealGeom), nSP, fp);
    fwrite(zSf, sizeof(RealGeom), nSP, fp);
    fwrite(&nSF, sizeof(IntType), 1, fp);
    fwrite(nSfP, sizeof(IntType), nSP + 1, fp);
    fwrite(SfP, sizeof(IntType), nSfP[nSP], fp);
    fwrite(nPntS, sizeof(IntType), nSF + 1, fp);
    fwrite(PntS, sizeof(IntType), nPntS[nSF], fp);
    fwrite(SP, sizeof(IntType), nSP, fp);

    fwrite(&nBox, sizeof(IntType), 1, fp);
    fwrite(nPBox, sizeof(IntType), nBox + 1, fp);
    for (i = 0; i < nBox; i++)
        fwrite(BBox[i], sizeof(RealGeom), 6, fp);

    fclose(fp);
    sdel_array_1D(xSf);
    sdel_array_1D(ySf);
    sdel_array_1D(zSf);
    sdel_array_1D(nSfP);
    sdel_array_1D(SfP);
    sdel_array_1D(nPntS);
    sdel_array_1D(PntS);
    sdel_array_1D(SP);
    sdel_array_1D(nPBox);
    sdel_array_2D(BBox);
}

//!< 计算网格单元体心到壁面的距离in the MPI.
void PolyGrid::ComputeCellDist()
{
#ifndef MF_MPICH
    WriteInfoDist();
#endif

    RealGeom *dist2wall_cell = 0;
    // dist2wall_cell = (RealGeom *)GetDataPtr(REAL_GEOM, nTCell, "dist2wall_cell");
    if (!dist2wall_cell)
    {
        snew_array_1D(dist2wall_cell, nTCell);
        // UpdateDataPtr(dist2wall_cell, REAL_GEOM, nTCell, "dist2wall_cell");
    }
    for (IntType j = 0; j < nTCell; j++)
        dist2wall_cell[j] = BIG;

    IntType mark = 1;
    ComputeDist2WallTriang(dist2wall_cell, mark);
}

void PolyGrid::Set_RecvSend(RealFlow ***bqs, RealFlow ***bqr, IntType nvar)
{
    IntType i, j, k, temp_nZIFace = 0;

    for (i = 0; i < nNeighbor; i++)
        temp_nZIFace += nZIFace[i];
    bqs[0] = NULL;
    bqr[0] = NULL;
    snew_array_1D(bqs[0], nvar * nNeighbor);
    snew_array_1D(bqr[0], nvar * nNeighbor);

    for (i = 1; i < nNeighbor; i++)
    {
        bqs[i] = &bqs[i - 1][nvar];
        bqr[i] = &bqr[i - 1][nvar];
    }
    bqs[0][0] = NULL;
    bqr[0][0] = NULL;
    snew_array_1D(bqs[0][0], nNeighbor * nvar * temp_nZIFace);
    snew_array_1D(bqr[0][0], nNeighbor * nvar * temp_nZIFace);
    for (i = 1; i < nNeighbor; i++)
    {
        bqs[i][0] = &bqs[i - 1][0][nvar * nZIFace[i - 1]];
        bqr[i][0] = &bqr[i - 1][0][nvar * nZIFace[i - 1]];
    }
    for (i = 0; i < nNeighbor; i++)
    {
        for (j = 1; j < nvar; j++)
        {
            bqs[i][j] = &bqs[i][j - 1][nZIFace[i]];
            bqr[i][j] = &bqr[i][j - 1][nZIFace[i]];
        }
    }

    for (i = 0; i < nNeighbor; i++)
    {
        for (j = 0; j < nvar; j++)
        {
            for (k = 0; k < nZIFace[i]; k++)
            {
                bqs[i][j][k] = 0.0;
                bqr[i][j][k] = 0.0;
            }
        }
    }
}

void PolyGrid::Add_RecvSend(RealFlow ***bqs, RealFlow *q, IntType num_var)
{
    IntType i, j;
    for (i = 0; i < nNeighbor; i++)
    {
        for (j = 0; j < nZIFace[i]; j++)
        {
            bqs[i][num_var][j] = q[bCNo[i][j]];
        }
    }
}

void PolyGrid::Read_RecvSend(RealFlow ***bqr, RealFlow *q, IntType num_var)
{
    IntType i, j, ghost;
    for (i = 0; i < nNeighbor; i++)
    {
        for (j = 0; j < nZIFace[i]; j++)
        {
            ghost = nTCell + bFNo[i][j];
            q[ghost] = bqr[i][num_var][j];
        }
    }
}

//!< 初始化非定常计算网格单元面心的速度和体心的速度
void PolyGrid::InitialVgn()
{
    if (vgn == 0)
    {
        IntType i;

        RealGeom BFacevgx_init = 0.0, BFacevgy_init = 0.0, BFacevgz_init = 0.0;

        // zhyb: 保存上一步的边界面速度，用于非定常流场中流体加速度对物面边界压力的影响
        // zhyb: 初始化为第一步的速度
        RealGeom *BFacevgx_last = NULL;
        RealGeom *BFacevgy_last = NULL;
        RealGeom *BFacevgz_last = NULL;
        snew_array_1D(BFacevgx_last, nBFace);
        snew_array_1D(BFacevgy_last, nBFace);
        snew_array_1D(BFacevgz_last, nBFace);
        for (i = 0; i < nBFace; i++)
        {
            BFacevgx_last[i] = BFacevgx_init;
            BFacevgy_last[i] = BFacevgy_init;
            BFacevgz_last[i] = BFacevgz_init;
        }
        // UpdateDataPtr(BFacevgx_last, REAL_GEOM, nBFace, "BFacevgx_last");
        // UpdateDataPtr(BFacevgy_last, REAL_GEOM, nBFace, "BFacevgy_last");
        // UpdateDataPtr(BFacevgz_last, REAL_GEOM, nBFace, "BFacevgz_last");

        vgn = NULL;
        BFacevgx = NULL;
        BFacevgy = NULL;
        BFacevgz = NULL;
        snew_array_1D(vgn, nTFace);
        snew_array_1D(BFacevgx, nBFace);
        snew_array_1D(BFacevgy, nBFace);
        snew_array_1D(BFacevgz, nBFace);
        for (i = 0; i < nBFace; i++)
        {
            BFacevgx[i] = BFacevgx_init;
            BFacevgy[i] = BFacevgy_init;
            BFacevgz[i] = BFacevgz_init;
        }
        for (i = 0; i < nTFace; i++)
        {
            vgn[i] = BFacevgx_init * xfn[i] + BFacevgy_init * yfn[i] + BFacevgz_init * zfn[i];
        }
    }
    else
    {
        // mflog::log.set_one_processor_out();
        // mflog::log << "Error! vgn!=0 in function InitialVgn!" << endl;
        // mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
        cout << "Error! vgn!=0 in function InitialVgn!" << endl;
        exit(1);
    }
}
*/

/*!
 * @brief       quicksort 排序法的原始子程序
 * @param       a
 * @param       is
 * @param       ie
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-06
 */
void PolyGrid::quick_sort(IntType *a, IntType is, IntType ie)
{
    IntType i, last;
    if (ie - is <= 1)
        return;
    swap(a, is, (ie + is) / 2);
    last = is;
    for (i = is + 1; i <= ie; i++)
    {
        if (a[i] < a[is])
            swap(a, ++last, i);
    }
    swap(a, is, last);
    quick_sort(a, is, last);
    quick_sort(a, last + 1, ie);
}

void PolyGrid::swap(IntType *a, IntType i, IntType j)
{
    IntType temp;
    temp = a[i];
    a[i] = a[j];
    a[j] = temp;
}
/**
//!< check grid quality
void PolyGrid::CheckGridQuality()
{

    if (level == 0)
        CheckGridScale();

    IntType n = nTCell + nBFace;
    if (facecentroidskewness == 0)
        snew_array_1D(facecentroidskewness, nTFace);
    if (faceangleskewness == 0)
        snew_array_1D(faceangleskewness, nTFace);
    if (cellcentroidskewness == 0)
        snew_array_1D(cellcentroidskewness, n);
    if (cellvolsmoothness == 0)
        snew_array_1D(cellvolsmoothness, n);
    if (cellwallnumber == 0)
        snew_array_1D(cellwallnumber, nTCell);
    if (faceangle == 0)
        snew_array_1D(faceangle, nTFace);

    SkewnessSummary();
    EquiangleSkewnessSummary();
    SmoothnessSummary();
    FindIllWallCell();
    CheckSymmetryFace();
    FaceAngleSummary();

    // deal bad grid for robust
    DealBadGrid();
}

//!< summary of skewness of a given grid
void PolyGrid::SkewnessSummary()
{
    //     IntType i, c1, c2, type, count;
    //     IntType angleC[18];
    //     IntType n = nTCell + nBFace;
    //     RealGeom dotp1, x1, y1, z1, dis1, dotp2, x2, y2, z2, dis2, angle, angle1, angle2;
    //     RealGeom minAng, maxAng;

    //     if (nNPC == 0)
    //         nNPC = CalnNPC(this);
    //     if (C2N == 0)
    //         C2N = CalC2N(this);

    //     // face skew check
    //     for (i = 0; i < 18; i++)
    //     {
    //         angleC[i] = 0;
    //     }

    //     // cell skew initialized
    //     for (i = 0; i < n; i++)
    //     {
    //         cellcentroidskewness[i] = 90.0;
    //     }

    //     minAng = 90.0;
    //     maxAng = -90.0;
    //     for (i = 0; i < nTFace; i++)
    //     {
    //         c1 = f2c[i + i];
    //         c2 = f2c[i + i + 1];

    //         x1 = xfc[i] - xcc[c1];
    //         y1 = yfc[i] - ycc[c1];
    //         z1 = zfc[i] - zcc[c1];
    //         dis1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
    //         dotp1 = (xfn[i] * x1 + yfn[i] * y1 + zfn[i] * z1) / (dis1 + TINY);
    //         dotp1 = MIN(dotp1, 1.0);
    //         dotp1 = MAX(dotp1, -1.0);
    //         angle1 = asin(dotp1) * 180 / PI;
    //         cellcentroidskewness[c1] = MIN(cellcentroidskewness[c1], angle1);

    //         x2 = xcc[c2] - xfc[i];
    //         y2 = ycc[c2] - yfc[i];
    //         z2 = zcc[c2] - zfc[i];
    //         dis2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
    //         dotp2 = (xfn[i] * x2 + yfn[i] * y2 + zfn[i] * z2) / (dis2 + TINY);
    //         dotp2 = MIN(dotp2, 1.0);
    //         dotp2 = MAX(dotp2, -1.0);
    //         angle2 = asin(dotp2) * 180 / PI;
    //         if (i >= nBFace)
    //         {
    //             cellcentroidskewness[c2] = MIN(cellcentroidskewness[c2], angle2);
    //         }

    //         angle = MIN(angle1, angle2);
    //         facecentroidskewness[i] = angle;

    //         minAng = MIN(minAng, angle);
    //         maxAng = MAX(maxAng, angle);

    //         angle = MIN(angle, 89.9);
    //         angle = MAX(angle, -89.9);
    //         if (angle < 0)
    //             angle -= 10;
    //         angleC[(IntType)(angle / 10) + 9]++;
    //     }

    //     mflog::log.set_one_processor_out();

    //     IntType nTFace_glb = nTFace;
    //     // #ifdef MF_MPICH
    //     //     Parallel::parallel_sum(nTFace_glb, MPI_COMM_WORLD);
    //     //     Parallel::parallel_sum(angleC, 18, MPI_COMM_WORLD);
    //     //     Parallel::parallel_min_max(minAng, maxAng, MPI_COMM_WORLD);
    //     // #endif

    //     mflog::log << SEP_LINE << endl;
    //     mflog::log << " Face Skewness Summary (face angle of 90 degrees being the best) " << endl;
    //     mflog::log << SEP_LINE << endl;
    //     mflog::log << "Total faces number: " << nTFace_glb << endl;
    //     mflog::log << "Min skewness angle: " << IOS_EP(6) << minAng << endl;
    //     mflog::log << "Max skewness angle: " << IOS_EP(6) << maxAng << endl;
    //     for (i = 0; i < 18; i++)
    //     {
    //         mflog::log << "Face angle from " << (int)(-90 + i * 10) << " to " << (int)(-80 + i * 10) << " is "
    //                    << IOS_FWP(5, 2) << angleC[i] / ((float)nTFace_glb) * 100 << " percent, " << (long)angleC[i]
    //                    << std::endl;
    //     }
    //     mflog::log << SEP_LINE << endl;

    //     if (level == 0)
    //     {
    //         RealGeom BadFaceAngle = -1.0;
    //         GetData(&BadFaceAngle, REAL_GEOM, 1, "BadFaceAngle", 0);

    //         count = 0;
    //         for (i = 0; i < nTFace; i++)
    //         {
    //             if (facecentroidskewness[i] < BadFaceAngle)
    //             {
    //                 count++;
    //             }
    //         }

    //         IntType count_glb = count;
    //         // #ifdef MF_MPICH
    //         //         Parallel::parallel_sum(count_glb, MPI_COMM_WORLD);
    //         // #endif
    //         mflog::log << "Bad Face Angle = " << BadFaceAngle << endl
    //                    << "The number of face centroid skewness angle less than BadFaceAngle is: "
    //                    << count_glb << endl;
    //     }

    //     // cell skew check
    //     // ghost cell's value
    //     for (i = 0; i < nBFace; i++)
    //     {
    //         c1 = f2c[i + i];
    //         c2 = f2c[i + i + 1];
    //         type = bcr[i]->GetType();
    //         if (type == INTERFACE)
    //             continue;

    //         cellcentroidskewness[c2] = cellcentroidskewness[c1];
    //     }
    //     // #ifdef MF_MPICH
    //     //     CommInterfaceDataMPI(cellcentroidskewness);
    //     // #endif

    //     for (i = 0; i < 18; i++)
    //         angleC[i] = 0;
    //     minAng = 90.0;
    //     maxAng = -90.0;
    //     for (i = 0; i < nTCell; i++)
    //     {
    //         minAng = MIN(minAng, cellcentroidskewness[i]);
    //         maxAng = MAX(maxAng, cellcentroidskewness[i]);

    //         angle = MIN(cellcentroidskewness[i], 89.9);
    //         angle = MAX(angle, -89.9);
    //         if (angle < 0)
    //             angle -= 10;
    //         angleC[(IntType)(angle / 10) + 9]++;
    //     }

    //     IntType nTCell_glb = nTCell;
    //     // #ifdef MF_MPICH
    //     //     Parallel::parallel_sum(nTCell_glb, MPI_COMM_WORLD);
    //     //     Parallel::parallel_sum(angleC, 18, MPI_COMM_WORLD);
    //     //     Parallel::parallel_min_max(minAng, maxAng, MPI_COMM_WORLD);
    //     // #endif

    //     mflog::log << SEP_LINE << endl;
    //     mflog::log << " Cell Skewness Summary(cell angle of 90 degrees being the best) " << endl;
    //     mflog::log << SEP_LINE << endl;
    //     mflog::log << "Total cell number: " << nTCell_glb << endl;
    //     mflog::log << "Min cell centroid skewness value is " << IOS_EP(6) << minAng << endl;
    //     mflog::log << "Max cell centroid skewness value is " << IOS_EP(6) << maxAng << endl;
    //     for (i = 0; i < 18; i++)
    //     {
    //         mflog::log << "Cell centroid skewness value from " << (int)(-90 + i * 10) << " to "
    //                    << (int)(-80 + i * 10) << " is " << IOS_FWP(5, 2) << angleC[i] / ((float)nTCell_glb) * 100
    //                    << " percent, " << (long)angleC[i] << std::endl;
    //     }
    //     mflog::log << SEP_LINE << endl;

    //     if (level == 0)
    //     {
    //         RealGeom BadCellAngle = 0.0;
    //         GetData(&BadCellAngle, REAL_GEOM, 1, "BadCellAngle", 0);

    //         count = 0;
    //         for (i = 0; i < nTCell; i++)
    //         {
    //             if (cellcentroidskewness[i] < BadCellAngle)
    //             {
    //                 count++;
    //             }
    //         }
    //         IntType count_glb = count;
    //         // #ifdef MF_MPICH
    //         //         Parallel::parallel_sum(count_glb, MPI_COMM_WORLD);
    //         // #endif
    //         mflog::log << "Bad Cell Angle = " << BadCellAngle << endl
    //                    << "The number of cell centroid skewness angle less than BadCellAngle is: "
    //                    << count_glb << endl;
    //     }
}

//!< summary of equi-angle skewness of a given grid
void PolyGrid::EquiangleSkewnessSummary()
{
    //     IntType i, j, p1, p2, p3;
    //     IntType face_skew[11];
    //     RealGeom xx1, yy1, zz1, l1, xx2, yy2, zz2, l2, dot;
    //     RealGeom angle, min_angle, max_angle, e_angle, min_faceskew, max_faceskew;

    //     if (F2N == 0)
    //     {
    //         F2N = CalF2N(this);
    //     }
    //     assert(F2N);

    //     RealGeom *x = GetX();
    //     RealGeom *y = GetY();
    //     RealGeom *z = GetZ();

    //     for (i = 0; i < 11; i++)
    //         face_skew[i] = 0;

    //     min_faceskew = 1.0;
    //     max_faceskew = 0.0;
    //     for (i = 0; i < nTFace; i++)
    //     {
    //         min_angle = 180.0;
    //         max_angle = 0.0;
    //         for (j = 0; j < nNPF[i]; j++)
    //         {
    //             p1 = F2N[i][j];
    //             if (j == 0)
    //             {
    //                 p2 = F2N[i][nNPF[i] - 1];
    //                 p3 = F2N[i][j + 1];
    //             }
    //             else if (j == nNPF[i] - 1)
    //             {
    //                 p2 = F2N[i][j - 1];
    //                 p3 = F2N[i][0];
    //             }
    //             else
    //             {
    //                 p2 = F2N[i][j - 1];
    //                 p3 = F2N[i][j + 1];
    //             }

    //             xx1 = x[p2] - x[p1];
    //             yy1 = y[p2] - y[p1];
    //             zz1 = z[p2] - z[p1];
    //             xx2 = x[p3] - x[p1];
    //             yy2 = y[p3] - y[p1];
    //             zz2 = z[p3] - z[p1];
    //             l1 = sqrt(xx1 * xx1 + yy1 * yy1 + zz1 * zz1);
    //             l2 = sqrt(xx2 * xx2 + yy2 * yy2 + zz2 * zz2);
    //             dot = xx1 * xx2 + yy1 * yy2 + zz1 * zz2;
    //             dot /= (l1 * l2);
    //             dot = MIN(dot, 1.0);
    //             dot = MAX(dot, -1.0);
    //             angle = acos(dot) * 180.0 / PI;
    //             max_angle = MAX(max_angle, angle);
    //             min_angle = MIN(min_angle, angle);
    //         }
    //         e_angle = 180.0 * (1.0 - 2.0 / (RealGeom)nNPF[i]);
    //         faceangleskewness[i] = MAX((max_angle - e_angle) / (180.0 - e_angle), (e_angle - min_angle) / e_angle);

    //         min_faceskew = MIN(min_faceskew, faceangleskewness[i]);
    //         max_faceskew = MAX(max_faceskew, faceangleskewness[i]);
    //         face_skew[(IntType)(faceangleskewness[i] * 10)]++;
    //     }
    //     // 是否删除faceangleskewness;lihuan-2018-11-26
    //     // mfmesh::sdel_array_1D(faceangleskewness);

    //     mflog::log.set_one_processor_out();

    //     IntType nTFace_glb = nTFace;
    //     // #ifdef MF_MPICH
    //     //     Parallel::parallel_sum(nTFace_glb, MPI_COMM_WORLD);
    //     //     Parallel::parallel_sum(face_skew, 10, MPI_COMM_WORLD);
    //     //     Parallel::parallel_min_max(min_faceskew, max_faceskew, MPI_COMM_WORLD);
    //     // #endif

    //     mflog::log << SEP_LINE << endl;
    //     mflog::log << "          Face Equiangle Skewness Summary" << endl;
    //     mflog::log << "(0.0 being the best, 1.0 is the worse, value below 0.9 is acceptable)" << endl;
    //     mflog::log << SEP_LINE << endl;
    //     mflog::log << "Min face equiangle skewness value is " << IOS_EP(6) << min_faceskew << endl;
    //     mflog::log << "Max face equiangle skewness value is " << IOS_EP(6) << max_faceskew << endl;
    //     for (i = 0; i < 10; i++)
    //     {
    //         mflog::log << "Face equiangle skewness value from " << IOS_FWP(3, 1) << i * 0.1 << " to "
    //                    << IOS_FWP(3, 1) << (i + 1) * 0.1 << " is " << IOS_FWP(5, 2) << face_skew[i] / ((float)nTFace_glb) * 100
    //                    << " percent, " << (long)face_skew[i] << std::endl;
    //     }
    //     mflog::log << SEP_LINE << endl;
}

//!< Summary of cell's smoothness of a given grid
void PolyGrid::SmoothnessSummary()
{
    //     IntType i, j, c1, c2, cell, type;
    //     IntType cell_skew[10];
    //     RealGeom min_vol, max_vol, min_cellskew, max_cellskew;

    //     if (nCPC == 0)
    //         nCPC = CalnCPC(this);
    //     if (c2c == 0)
    //         c2c = CalC2C(this);

    //     for (i = 0; i < 10; i++)
    //         cell_skew[i] = 0;

    //     min_vol = BIG;
    //     max_vol = 0.0;
    //     for (i = 0; i < nTCell; i++)
    //     {
    //         min_vol = MIN(min_vol, vol[i]);
    //         max_vol = MAX(max_vol, vol[i]);

    //         cellvolsmoothness[i] = BIG;
    //     }

    //     for (i = 0; i < nTCell; i++)
    //     {
    //         for (j = 0; j < nCPC[i]; j++)
    //         {
    //             cell = c2c[i][j];
    //             cellvolsmoothness[i] = MIN(cellvolsmoothness[i], MIN(vol[cell], vol[i]) / MAX(vol[cell], vol[i]));
    //         }
    //         cellvolsmoothness[i] = 1.0 - cellvolsmoothness[i];
    //     }

    //     for (i = 0; i < nBFace; i++)
    //     {
    //         c1 = f2c[i + i];
    //         c2 = f2c[i + i + 1];
    //         type = bcr[i]->GetType();
    //         if (type == INTERFACE)
    //             continue;

    //         cellvolsmoothness[c2] = cellvolsmoothness[c1];
    //     }
    //     // #ifdef MF_MPICH
    //     //     CommInterfaceDataMPI(cellvolsmoothness);
    //     // #endif

    //     min_cellskew = 1.0;
    //     max_cellskew = 0.0;
    //     for (i = 0; i < nTCell; i++)
    //     {
    //         min_cellskew = MIN(min_cellskew, cellvolsmoothness[i]);
    //         max_cellskew = MAX(max_cellskew, cellvolsmoothness[i]);
    //         cell_skew[(IntType)(cellvolsmoothness[i] * 10)]++;
    //     }

    //     mflog::log.set_one_processor_out();

    //     IntType nTCell_glb = nTCell;
    //     // #ifdef MF_MPICH
    //     //     Parallel::parallel_sum(nTCell_glb, MPI_COMM_WORLD);
    //     //     Parallel::parallel_sum(cell_skew, 10, MPI_COMM_WORLD);
    //     //     Parallel::parallel_min_max(min_cellskew, max_cellskew, MPI_COMM_WORLD);
    //     //     Parallel::parallel_min_max(min_vol, max_vol, MPI_COMM_WORLD);
    //     // #endif

    //     mflog::log << SEP_LINE << endl;
    //     mflog::log << "             Cell volume Smoothness Summary" << endl;
    //     mflog::log << "(0.0 being the best, 1.0 is the worse, value below 0.8 is acceptable)" << endl;
    //     mflog::log << SEP_LINE << endl;
    //     mflog::log << "Min cell volume value is " << IOS_EP(6) << min_vol << endl;
    //     mflog::log << "Max cell volume value is " << IOS_EP(6) << max_vol << endl;
    //     mflog::log << "Min cell volume smoothness value is " << IOS_EP(6) << min_cellskew << endl;
    //     mflog::log << "Max cell volume smoothness value is " << IOS_EP(6) << max_cellskew << endl;

    //     for (i = 0; i < 10; i++)
    //     {
    //         mflog::log << "Cell volume smoothness value from " << IOS_FWP(3, 1) << i * 0.1 << " to "
    //                    << IOS_FWP(3, 1) << (i + 1) * 0.1 << " is " << IOS_FWP(5, 2) << cell_skew[i] / ((float)nTCell_glb) * 100
    //                    << " percent, " << (long)cell_skew[i] << std::endl;
    //     }
    //     mflog::log << SEP_LINE << endl;
}

//!< find ill wall cell, viz. the cell that has two or more wall faces
void PolyGrid::FindIllWallCell()
{
    IntType i, type, c1;
    IntType *nNPC = CalnNPC(this);

    IntType vis_mode = LAMINAR;
    GetData(&vis_mode, INT, 1, "vis_mode");

    for (i = 0; i < nTCell; i++)
        cellwallnumber[i] = 0;

    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        c1 = f2c[i + i];
        if (type != WALL)
            continue;
        cellwallnumber[c1]++;
    }

    IntType num_tetra = 0, num_pyra = 0, num_wall2 = 0;
    for (i = 0; i < nTCell; i++)
    {
        if (cellwallnumber[i] == 0)
            continue;

        if (vis_mode != INVISCID)
        {
            if (nNPC[i] == 4)
            {
                num_tetra++;
            }
            else if (nNPC[i] == 5)
            {
                num_pyra++;
            }
        }

        if (cellwallnumber[i] > 1)
        {
            num_wall2++;
        }
    }

    mflog::log.set_one_processor_out();

    // #ifdef MF_MPICH
    //     Parallel::parallel_sum(num_wall2, MPI_COMM_WORLD);
    // #endif
    if (num_wall2 > 0)
    {
        mflog::log << endl
                   << "total have " << num_wall2 << " cell has two or more face on wall!!!" << endl;
    }

    if (vis_mode != INVISCID)
    {
        IntType tet_pyr_total[2] = {num_tetra, num_pyra};
        // #ifdef MF_MPICH
        //         Parallel::parallel_sum(tet_pyr_total, 2, MPI_COMM_WORLD);
        // #endif
        mflog::log << endl
                   << "total have " << tet_pyr_total[0] << " tetrahedron cell on wall!" << endl;
        mflog::log << endl
                   << "total have " << tet_pyr_total[1] << " pyramid cell on wall!" << endl;
    }
}

//!< check symmetry face, find the max distant of symmetry boundary node from symmetry plane
void PolyGrid::CheckSymmetryFace()
{
    IntType i, j, type, count, p1, nsymm;
    RealGeom dmax, maxx = 0, maxy = 0, maxz = 0;

    RealGeom *x = GetX();
    RealGeom *y = GetY();
    RealGeom *z = GetZ();

    count = 0;
    nsymm = 0;
    dmax = 0.0;
    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        if (type == SYMM)
        {
            nsymm++;
            for (j = 0; j < nNPF[i]; j++)
            {
                p1 = f2n[count++];
                if (dmax < fabs(y[p1]))
                {
                    dmax = fabs(y[p1]);
                    maxx = x[p1];
                    maxy = y[p1];
                    maxz = z[p1];
                }
            }
        }
        else
        {
            count += nNPF[i];
        }
    }

    mflog::log.set_all_processors_out();

    // #ifdef MF_MPICH
    //     struct
    //     {
    //         RealGeom dmax;
    //         IntType rank;
    //     } in, out;
    //     in.dmax = dmax;
    //     in.rank = myZone;
    // #ifdef SINGLE_PRECISION
    //     MPI_Allreduce(&in, &out, 1, MPI_FLOAT_INT, MPI_MINLOC, MPI_COMM_WORLD);
    // #else
    //     MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    // #endif
    //     IntType nsymm_total;
    //     MPI_Allreduce(&nsymm, &nsymm_total, 1, MPIIntType, MPI_SUM, MPI_COMM_WORLD);

    //     if (myZone == out.rank && nsymm_total != 0)
    //     {
    //         mflog::log << endl
    //                    << "There have " << nsymm_total << "symmetry boundary face."
    //                    << endl
    //                    << "The max distant of symmetry boundary node from symmetry plane is: " << IOS_EP(8) << maxy
    //                    << endl
    //                    << "The coordinate is: " << maxx << ", " << maxy << ", " << maxz << endl
    //                    << endl;
    //     }
    // #else
    if (nsymm != 0)
    {
        mflog::log << endl
                   << "There have " << nsymm << "symmetry boundary face."
                   << endl
                   << "The max distant of symmetry boundary node from symmetry plane is: " << IOS_EP(8) << maxy
                   << endl
                   << "The coordinate is: " << maxx << ", " << maxy << ", " << maxz << endl
                   << endl;
    }
    // #endif
}

//!< check grid scale, estimate precision enough or not
void PolyGrid::CheckGridScale()
{
    IntType i, j, k, type, p1, p2, c1, count;
    RealGeom len, lenscal, maxx, minx, maxy, miny, maxz, minz;

    RealGeom ratio;
// #ifdef SINGLE_PRECISION
//     ratio = 1.0e-7;
#else // DOUBLE_PRECISION
    ratio = 1.0e-14;
    // #endif
    if (nNPC == 0)
        nNPC = CalnNPC(this);
    if (C2N == 0)
        C2N = CalC2N(this);
    if (F2N == 0)
        F2N = CalF2N(this);

    RealGeom *x = GetX();
    RealGeom *y = GetY();
    RealGeom *z = GetZ();

    maxx = -BIG;
    minx = BIG;
    maxy = -BIG;
    miny = BIG;
    maxz = -BIG;
    minz = BIG;
    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        if (type != WALL)
            continue;

        for (j = 0; j < nNPF[i]; j++)
        {
            p1 = F2N[i][j];
            maxx = MAX(maxx, x[p1]);
            minx = MIN(minx, x[p1]);
            maxy = MAX(maxy, y[p1]);
            miny = MIN(miny, y[p1]);
            maxz = MAX(maxz, z[p1]);
            minz = MIN(minz, z[p1]);
        }
    }

    mflog::log.set_one_processor_out();

    // #ifdef MF_MPICH
    //     Parallel::parallel_min_max(minx, maxx, MPI_COMM_WORLD);
    //     Parallel::parallel_min_max(miny, maxy, MPI_COMM_WORLD);
    //     Parallel::parallel_min_max(minz, maxz, MPI_COMM_WORLD);
    // #endif
    lenscal = MAX(MAX(maxx - minx, maxy - miny), maxz - minz);
    mflog::log << endl
               << "length scale is: " << IOS_EP(4) << lenscal << endl
               << endl;

    count = 0;
    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        if (type != WALL)
            continue;

        c1 = f2c[i + i];
        for (j = 0; j < nNPC[c1] - 1; j++)
        {
            p1 = C2N[c1][j];
            for (k = j + 1; k < nNPC[c1]; k++)
            {
                p2 = C2N[c1][k];
                len = sqrt((x[p2] - x[p1]) * (x[p2] - x[p1]) +
                           (y[p2] - y[p1]) * (y[p2] - y[p1]) +
                           (z[p2] - z[p1]) * (z[p2] - z[p1]));
                if (len / lenscal < ratio)
                {
                    mflog::log.set_all_processors_out();
                    mflog::log << endl
                               << "p1 and p2 is too close!" << endl;
                    mflog::log << "p1 " << p1 << IOS_EP(10) << x[p1] << IOS_EP(10) << y[p1] << IOS_EP(10) << z[p1] << endl;
                    mflog::log << "p2 " << p2 << IOS_EP(10) << x[p2] << IOS_EP(10) << y[p2] << IOS_EP(10) << z[p2] << endl;
                    count++;
                }
            }
        }
    }
    // #ifdef MF_MPICH
    //     Parallel::parallel_sum(count, MPI_COMM_WORLD);
    // #endif
    if (count != 0)
    {
        mflog::log.set_one_processor_out();
        mflog::log << endl
                   << "There has " << count << " node too fine!" << endl
                   << "You should use higher precision or make grid's first layer being coarser!" << endl;
        mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
    }
}

//!< Find Cell Layer No. from wall, No.0,1,2,3...
void PolyGrid::FindCellLayerNo()
{
    IntType i, j, c1, c2, type, mark, count;

    IntType n = nBFace + nTCell;
    IntType nPFace = nBFace - nIFace;
    IntType *nCPC = CalnCPC(this); // 注意：只带并行边界的虚拟网格
    IntType **C2C = CalC2C(this);

    IntType *CellLayerNo = NULL;
    snew_array_1D(CellLayerNo, n);

    // UpdateDataPtr(CellLayerNo, INT, n, "CellLayerNo");
    for (i = 0; i < n; i++)
        CellLayerNo[i] = -1;

    for (i = 0; i < nPFace; i++)
    {
        type = bcr[i]->GetType();
        if (type != WALL)
            continue;
        c1 = f2c[i + i];
        CellLayerNo[c1] = 0;
    }

#ifdef MF_MPICH
    CommInterfaceDataMPI(CellLayerNo);
#endif

    count = 0;
    mark = 1;
    while (mark)
    {
        mark = 0;
        count++;

        for (i = nPFace; i < nBFace; i++)
        {
            c1 = f2c[i + i];
            c2 = f2c[i + i + 1];

            if (CellLayerNo[c2] != count - 1)
                continue;

            if (CellLayerNo[c1] == -1)
            {
                CellLayerNo[c1] = count;
                mark = 1;
            }
        }

        for (i = 0; i < nTCell; i++)
        {
            if (CellLayerNo[i] != count - 1)
                continue;

            for (j = 0; j < nCPC[i]; j++)
            {
                c1 = C2C[i][j];

                if (c1 > nTCell)
                    continue;

                if (CellLayerNo[c1] == -1)
                {
                    CellLayerNo[c1] = count;
                    mark = 1;
                }
            }
        }
        // #ifdef MF_MPICH
        //         CommInterfaceDataMPI(CellLayerNo);

        //         // 只要还有更新的，就继续循环
        //         IntType mark_glb;
        //         MPI_Allreduce(&mark, &mark_glb, 1, MPIIntType, MPI_MAX, MPI_COMM_WORLD);
        //         mark = mark_glb;
        // #endif
    }

#ifdef MF_DEBUG
// mflog::log.set_one_processor_out();
#ifdef MF_MPICH
    if (myid == 0)
    {
#endif
        cout << endl
             << "Find cell layer no. Ok!" << endl;
#ifdef MF_MPICH
    }
#endif
#endif
}

//!< Find  face parallel to wall face used for Roe's entropy fix
void PolyGrid::FindNormalFace()
{
    // IntType CellNum = 5;
    IntType CellNum = 20;
    IntType i, j, face1, face2, mark, type;

    IntType *nFPC = CalnFPC(this);
    IntType **C2F = CalC2F(this);
    IntType n = nTCell + nBFace;

    IntType *CellLayerNo;
    // CellLayerNo = (IntType *)this->GetDataPtr(INT, n, "CellLayerNo");

    IntType *IsNormalFace = NULL;
    snew_array_1D(IsNormalFace, nTFace);
    // this->UpdateDataPtr(IsNormalFace, INT, nTFace, "IsNormalFace");
    for (i = 0; i < nTFace; i++)
        IsNormalFace[i] = 0;

    // 首先物面边界面都是
    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        if (type != WALL)
            continue;

        IsNormalFace[i] = 1;
    }

    // 附面层里面积最大的两个面一般就是垂直于物面法向的面
    for (i = 0; i < nTCell; i++)
    {
        if ((CellLayerNo[i] == -1) || (CellLayerNo[i] >= CellNum))
            continue;
        if (nNPC[i] < 6)
            continue; // 排除四面体和金字塔

        // 面积最大的面
        face1 = C2F[i][0];
        mark = 0;
        for (j = 1; j < nFPC[i]; j++)
        {
            face2 = C2F[i][j];
            if (area[face1] < area[face2])
            {
                face1 = face2;
                mark = j;
            }
        }
        IsNormalFace[face1] = 1;

        // 面积第二大的面
        if (mark == 0)
        {
            face1 = C2F[i][1];
        }
        else
        {
            face1 = C2F[i][0];
        }
        for (j = 0; j < nFPC[i]; j++)
        {
            if (j == mark)
                continue;
            face2 = C2F[i][j];
            if (area[face1] < area[face2])
            {
                face1 = face2;
                mark = j;
            }
        }
        IsNormalFace[face1] = 1;
    }
#ifdef MF_MPICH
    IntType *mpitmp = NULL;
    snew_array_1D(mpitmp, n);
    RealGeom *mpifc[3];
    for (i = 0; i < 3; i++)
    {
        mpifc[i] = NULL;
        snew_array_1D(mpifc[i], n);
    }
    for (i = 0; i < n; i++)
    {
        mpitmp[i] = 0;
        for (j = 0; j < 3; j++)
        {
            mpifc[j][i] = 0.0;
        }
    }

    for (i = 0; i < nBFace; i++)
    {
        IntType c1 = f2c[i + i];
        mpitmp[c1] = IsNormalFace[i];
        mpifc[0][c1] = xfc[i];
        mpifc[1][c1] = yfc[i];
        mpifc[2][c1] = zfc[i];
    }

    this->CommInterfaceDataMPI(mpitmp);
    for (j = 0; j < 3; j++)
    {
        this->CommInterfaceDataMPI(mpifc[j]);
    }

    for (i = 0; i < nBFace; i++)
    {
        IntType c2 = f2c[i + i + 1];
        if (mpitmp[c2] == 1)
        {
            if (fabs(mpifc[0][c2] - xfc[i]) < TINY &&
                fabs(mpifc[1][c2] - yfc[i]) < TINY &&
                fabs(mpifc[2][c2] - zfc[i]) < TINY)
            {
                IsNormalFace[i] = 1;
            }
        }
    }
    sdel_array_1D(mpitmp);
    for (i = 0; i < 3; i++)
    {
        sdel_array_1D(mpifc[i]);
    }
#endif
}

//!< summary of angle of lines between face center with two center center
void PolyGrid::FaceAngleSummary()
{
    IntType i, c1, c2;
    IntType angleC[18];
    RealGeom x1, y1, z1, x2, y2, z2, dot, dis1, dis2;
    RealGeom minang, maxang, angle;

    for (i = 0; i < 18; i++)
        angleC[i] = 0;

    minang = 180.0;
    maxang = 0.0;
    for (i = 0; i < nTFace; i++)
    {
        c1 = f2c[i + i];
        c2 = f2c[i + i + 1];

        x1 = xfc[i] - xcc[c1];
        y1 = yfc[i] - ycc[c1];
        z1 = zfc[i] - zcc[c1];
        x2 = xfc[i] - xcc[c2];
        y2 = yfc[i] - ycc[c2];
        z2 = zfc[i] - zcc[c2];

        dot = x1 * x2 + y1 * y2 + z1 * z2;
        dis1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
        dis2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
        faceangle[i] = dot / (dis1 + TINY) / (dis2 + TINY);
        faceangle[i] = MIN(faceangle[i], 1.0);
        faceangle[i] = MAX(faceangle[i], -1.0);
        faceangle[i] = acos(faceangle[i]) * 180.0 / PI;

        minang = MIN(minang, faceangle[i]);
        maxang = MAX(maxang, faceangle[i]);

        angle = MIN(179.0, faceangle[i]);
        angleC[(IntType)(angle / 10)]++;
    }
    // 是否删除faceangle;lihuan-2018-11-26
    // mfmesh::sdel_array_1D(faceangle);

    mflog::log.set_one_processor_out();

    IntType nTFace_glb = nTFace;
#ifdef MF_MPICH
    Parallel::parallel_sum(nTFace_glb, MPI_COMM_WORLD);
    Parallel::parallel_sum(angleC, 18, MPI_COMM_WORLD);
    Parallel::parallel_min_max(minang, maxang, MPI_COMM_WORLD);
#endif

    mflog::log << SEP_LINE << endl;
    mflog::log << " Face Angle Summary (angle of 180 degrees being the best) " << endl;
    mflog::log << SEP_LINE << endl;
    mflog::log << "Total number of faces " << nTFace << endl;
    mflog::log << "Min face angle is " << IOS_EP(6) << minang << endl;
    mflog::log << "Max face angle is " << IOS_EP(6) << maxang << endl;
    for (i = 0; i < 18; i++)
    {

        mflog::log << "Face angle from " << (int)(i * 10) << " to " << (int)(i * 10 + 10)
                   << " is " << IOS_FWP(5, 2) << angleC[i] / ((float)nTFace) * 100 << " percent, "
                   << (long)angleC[i] << std::endl;
    }
    mflog::log << SEP_LINE << endl;
}

//!< 根据网格质量判断在做重构时是否降阶
void PolyGrid::DealBadGrid()
{
    IntType i;

    IntType *IfOneOrder = NULL;
    // IfOneOrder = (IntType *)this->GetDataPtr(INT, nTFace, "IfOneOrder");
    if (!IfOneOrder)
    {
        snew_array_1D(IfOneOrder, nTFace);
        // UpdateDataPtr(IfOneOrder, INT, nTFace, "IfOneOrder");
    }

    // 面心与左右两侧体心连线的夹角>100度
    for (i = 0; i < nTFace; i++)
    {
        IfOneOrder[i] = 1;
    }
}

*/

/*!
 * @brief       compute some additional info for geomety for efficiency
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-15
 */
/**
void PolyGrid::AdditionalInfoForGeometry()
{

    CalculateVolumnAverage(); // Calculate volume average of cell

    CalNormalDistanceOfC2C(); // calculate normal distance of two cell
}

//!< Calculate volume average of cell used for volume reference value in venkatakrishnan's limiter
RealGeom PolyGrid::CalculateVolumnAverage()
{
    IntType nTCell = this->GetNTCell();

    IntType i;

    RealGeom volume_average = 0.0;
    for (i = 0; i < nTCell; i++)
        volume_average += vol[i];
#ifdef MF_MPICH
    RealGeom vol_glb;
    IntType nTCell_glb;
    MPI_Allreduce(&volume_average, &vol_glb, 1, MPIReal, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nTCell, &nTCell_glb, 1, MPIIntType, MPI_SUM, MPI_COMM_WORLD);
    volume_average = vol_glb / nTCell_glb;
#else
    volume_average /= nTCell;
#endif

    SetVolAvg(volume_average);
    return volume_average;
}

//!< calculate normal distance of two cell, used for viscous spectral radii in LUSGS
RealGeom *PolyGrid::CalNormalDistanceOfC2C()
{
    IntType nTFace = GetNTFace();
    IntType *f2c = Getf2c();
    RealGeom *xfn = GetXfn();
    RealGeom *yfn = GetYfn();
    RealGeom *zfn = GetZfn();
    RealGeom *xcc = GetXcc();
    RealGeom *ycc = GetYcc();
    RealGeom *zcc = GetZcc();

    IntType i, c1, c2;

    // RealGeom *norm_dist_c2c = NULL;
    // norm_dist_c2c = (RealGeom *)this->GetDataPtr(REAL_GEOM, nTFace, "norm_dist_c2c");
    if (!norm_dist_c2c)
    {
        snew_array_1D(norm_dist_c2c, nTFace);
        // UpdateDataPtr(norm_dist_c2c, REAL_GEOM, nTFace, "norm_dist_c2c");
    }

    for (i = 0; i < nTFace; i++)
    {
        c1 = f2c[i + i];
        c2 = f2c[i + i + 1];

        norm_dist_c2c[i] = fabs(xfn[i] * (xcc[c2] - xcc[c1]) + yfn[i] * (ycc[c2] - ycc[c1]) + zfn[i] * (zcc[c2] - zcc[c1]));
    }

    return norm_dist_c2c;
}
*/

/**
void PolyGrid::PartitionGrids(IntType *CellToZone, IntType n_zone)
{
    idx_t *xadj = NULL;
    idx_t *adjncy = NULL;
    idx_t *adjwgt = NULL;
    snew_array_1D(xadj, nTCell + 1);
    snew_array_1D(adjncy, 2 * (nTFace - nBFace));
    snew_array_1D(adjwgt, 2 * (nTFace - nBFace));

    Getxadjadjncy(xadj, adjncy, adjwgt);

    SerialMetis(xadj, adjncy, adjwgt, n_zone, CellToZone);

    sdel_array_1D(xadj);
    sdel_array_1D(adjncy);
    sdel_array_1D(adjwgt);
}

//!< get xadj, adjncy and adjwgt
void PolyGrid::Getxadjadjncy(idx_t *xadj, idx_t *adjncy, idx_t *adjwgt)
{
    register IntType i, j, k;
    IntType *nCPC = CalnCPC(this);
    IntType **c2c = CalC2C(this);

    xadj[0] = 0;
    k = 0;
    for (i = 0; i < nTCell; i++)
    {
        xadj[i + 1] = xadj[i] + nCPC[i];

        for (j = 0; j < nCPC[i]; j++)
            adjncy[k++] = c2c[i][j];
    }

    // Now get adjwgt
    if (level == 0)
    { // fine grid
        IntType c1, c2;
        RealGeom dx, dy, dz, d;
        RealGeom maxlen = -1.0;

        RealGeom *len = NULL;
        snew_array_1D(len, k);

        k = 0;
        for (i = 0; i < nTCell; i++)
        {
            c1 = i;

            for (j = 0; j < nCPC[i]; j++)
            {
                c2 = c2c[i][j];

                dx = xcc[c2] - xcc[c1];
                dy = ycc[c2] - ycc[c1];
                dz = zcc[c2] - zcc[c1];
                d = sqrt(dx * dx + dy * dy + dz * dz);

                len[k++] = d;
                maxlen = max(maxlen, d);
            }
        }
        maxlen += TINY;

        k = 0;
        for (i = 0; i < nTCell; i++)
        {
            for (j = 0; j < nCPC[i]; j++)
            {
                adjwgt[k] = (idx_t)(sqrt(maxlen / len[k]));
                ++k;
            }
        }

        sdel_array_1D(len);
    }
    else
    { // coarse grid
        k = 0;
        for (i = 0; i < nTCell; i++)
        {
            for (j = 0; j < nCPC[i]; j++)
            {
                adjwgt[k++] = 1;
            }
        }
    }

    SetnCPC(NULL);
    Setc2c(NULL);
    c2c = NULL;
    nCPC = NULL;
}

//!< serial metis
void PolyGrid::SerialMetis(idx_t *xadj, idx_t *adjncy, idx_t *adjwgt, IntType n_zone, IntType *CellToZone)
{
    IntType i, status;
    idx_t nvtxs, ncon, nparts, objval;
    idx_t *part;

    mflog::log.set_one_processor_out();
    mflog::log << endl
               << "Now beginning partition graph!" << endl;

    nvtxs = (idx_t)nTCell;
    ncon = 1;
    nparts = (idx_t)n_zone;
    part = NULL;
    snew_array_1D(part, static_cast<size_t>(nvtxs));
    for (i = 0; i < nvtxs; i++)
    {
        part[i] = 0;
    }

    if (n_zone > 1)
    { // zhyb:n_zone==1时不划分网格，只输出！
        // choose Recursive bisection for nparts<64, k-way for nparts>=64, refer to EDGE
        if (n_zone < 64)
        {
            mflog::log << "Using Recursive Partitioning!" << endl;
            status = METIS_PartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, NULL, NULL, adjwgt,
                                              &nparts, NULL, NULL, NULL, &objval, part);
        }
        else
        {
            mflog::log << "Using k-way Partitioning!" << endl;
            status = METIS_PartGraphKway(&nvtxs, &ncon, xadj, adjncy, NULL, NULL, adjwgt,
                                         &nparts, NULL, NULL, NULL, &objval, part);
        }

        if (status != METIS_OK)
        {
            mflog::log << "Metis partition error! status code: " << status << endl;
            mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
        }
    }

    for (i = 0; i < nTCell; i++)
    {
        CellToZone[i] = (IntType)part[i];
    }

    mflog::log << "The interface number: " << objval << endl;
    mflog::log << "Serial Metis is OK!" << endl;
    sdel_array_1D(part);
}
*/

/**
#ifdef OMP_CellColor
//     Add by DZ 2021-8-20
//    Reference：AIAA94-0645
void PolyGrid::ReorderCellforLUSGS_1()
{
    // IntType nTCell = grid->GetNTCell();
    // IntType nBFace = grid->GetNBFace();

    IntType i, c1, c2, type, count, status1, status2, status3, status4;

    IntType n = nBFace + nTCell;
    IntType BigLayer = 1000 * nTCell;
    IntType NowLayer = 0;
    IntType MaxLayer = 0;

    IntType *Layer = NULL;
    IntType *LUSGSCellOrder = NULL;
    snew_array_1D(Layer, n);
    for (i = 0; i < n; i++)
        Layer[i] = BigLayer;
    snew_array_1D(LUSGSCellOrder, nTCell);
    for (i = 0; i < nTCell; i++)
        LUSGSCellOrder[i] = -1;
    // 找任意一个物面单元作为起始层
    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        if (type == WALL)
        {
            c1 = f2c[i + i];
            Layer[c1] = NowLayer++;
            break;
        }
    }
    // 如果没有物面，则找第一个边界单元作为起始层
    if (NowLayer == 0)
    {
        c1 = f2c[0 + 0];
        Layer[c1] = NowLayer++;
    }

    mflog::log.set_all_processors_out();

    while (1)
    {
        status1 = 0;
        count = 2 * nBFace;
        for (i = nBFace; i < nTFace; i++)
        {
            c1 = f2c[count++];
            c2 = f2c[count++];

            if (Layer[c1] < NowLayer && Layer[c2] < NowLayer)
            {
                continue;
            }
            else if (Layer[c1] > NowLayer && Layer[c2] > NowLayer)
            {
                continue;
            }
            else if (Layer[c1] == NowLayer || Layer[c2] == NowLayer)
            {
                continue;
            }
            else if (Layer[c1] == BigLayer)
            {
                Layer[c1] = NowLayer;
                status1++;
            }
            else if (Layer[c2] == BigLayer)
            {
                Layer[c2] = NowLayer;
                status1++;
            }
        }

        // 检查是否还有同层相邻的单元
        status2 = 0;
        for (i = nBFace; i < nTFace; i++)
        {
            c1 = f2c[i + i];
            c2 = f2c[i + i + 1];

            if (Layer[c1] == NowLayer && Layer[c1] == Layer[c2])
            {
                Layer[c2]++;
                status2++;
            }
        }
        // 检查是否有遗漏的单元
        status3 = 0;
        if (status1 == 0 && status2 == 0)
        {
            for (i = 0; i < nTCell; i++)
            {
                if (Layer[i] == BigLayer)
                {
                    status3++;
                    break;
                }
            }
        }

        // 如果没有，就跳出循环
        if (status1 == 0 && status2 == 0 && status3 == 0)
            break;

        // 该块存在多区不联通，需要再次启动
        if (status1 == 0 && status2 == 0 && status3 != 0)
        {
            status4 = 0;
            // 找任意一个物面单元作为起始层
            for (i = 0; i < nBFace; i++)
            {
                type = bcr[i]->GetType();
                if (type == WALL)
                {
                    c1 = f2c[i + i];
                    if (Layer[c1] == BigLayer)
                    {
                        Layer[c1] = 0;
                        status4++;
                        break;
                    }
                }
            }
            // 如果没有物面，则找任意边界单元作为起始层
            if (status4 == 0)
            {
                for (i = 0; i < nBFace; i++)
                {
                    c1 = f2c[i + i];
                    if (Layer[c1] == BigLayer)
                    {
                        Layer[c1] = 0;
                        status4++;
                        break;
                    }
                }
            }

            NowLayer = 0;
            // 肯定有边界单元存在，否则报错退出
            if (status4 == 0)
            {
#ifdef MF_MPICH
                mflog::log << endl
                           << "Error! Zone " << myZone << " is a multi-block, but I do not find"
                           << " the boundary cell for layer!" << endl;
#else
                mflog::log << endl
                           << "Error! This is a multi-block, but I do not find"
                           << " the boundary cell for layer!" << endl;
#endif
                mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
            }
        }

        NowLayer++;
    }

    for (i = 0; i < nTCell; i++)
    {
        MaxLayer = MAX(MaxLayer, Layer[i]);
    }
    MaxLayer++;

    for (i = 0; i < nTCell; i++)
    {
        if (Layer[i] > MaxLayer)
        {
#ifdef MF_MPICH
            mflog::log << endl
                       << "Wrong! Zone" << myZone << "Layer[" << i << "]=" << Layer[i] << endl;
#else
            mflog::log << endl
                       << "Wrong! Layer[" << i << "]=" << Layer[i] << endl;
#endif
        }
    }

    // 虚拟单元排在最后一层
    for (i = 0; i < nBFace; i++)
    {
        c2 = f2c[i + i + 1];
        Layer[c2] = MaxLayer;
    }

    IntType *cellsPerlayer = NULL;
    snew_array_1D(cellsPerlayer, nTCell);
    cellsPerlayer[0] = MaxLayer;
    cellsPerlayer[1] = 0;

    IntType temp;
    // 迭代顺序从低层到高层
    count = 0;
    for (NowLayer = 0; NowLayer < MaxLayer; NowLayer++)
    {
        temp = 0;
        for (i = 0; i < nTCell; i++)
        {
            if (Layer[i] == NowLayer)
            {
                LUSGSCellOrder[count++] = i;
                temp++;
            }
        }
        cellsPerlayer[NowLayer + 2] = count;
    }
    this->UpdateDataPtr(cellsPerlayer, INT, nTCell, "LUSGScellsPerlayer");

    for (i = 0; i < nTCell; i++)
    {
        Layer[LUSGSCellOrder[i]] = i;
    }
    for (i = nTCell; i < n; i++)
    {
        Layer[i] = i;
    }
    printf("LUSGS reordered by hyperplane approach 1!\n");
    this->UpdateDataPtr(LUSGSCellOrder, INT, nTCell, "LUSGSCellOrder");
    this->UpdateDataPtr(Layer, INT, n, "LUSGSLayer");
}

//   网格重排序
//   参考文献：航空计算技术35卷第3期，基于非结构网格流场计算的网格重排序，方法1
void PolyGrid::ReorderCellforLUSGS_2()
{
    // IntType nTCell = grid->GetNTCell();
    // IntType nBFace = grid->GetNBFace();

    IntType i, c1, c2, type, count, status1, status2, status3;
    IntType n = nBFace + nTCell;
    IntType BigLayer = 1000 * nTCell, NowLayer = 0, LastLayer, MaxLayer = 0;
    IntType *Layer = NULL;
    IntType *LUSGSCellOrder = NULL;
    snew_array_1D(Layer, n);
    for (i = 0; i < n; i++)
        Layer[i] = BigLayer;
    snew_array_1D(LUSGSCellOrder, nTCell);
    for (i = 0; i < nTCell; i++)
        LUSGSCellOrder[i] = -1;
    // 找第一个物面单元作为起始层
    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        if (type == WALL)
        {
            c1 = f2c[i + i];
            Layer[c1] = NowLayer++;
            break;
        }
    }
    // 如果没有物面，则找第一个非内部边界单元作为起始层
    if (NowLayer == 0)
    {
        for (i = 0; i < nBFace; i++)
        {
            type = bcr[i]->GetType();
            if (type != INTERFACE)
            {
                c1 = f2c[i + i];
                Layer[c1] = NowLayer++;
                break;
            }
        }
    }
    // 如果全是内部边界，则找第一个边界单元作为起始层
    if (NowLayer == 0)
    {
        c1 = f2c[0 + 0];
        Layer[c1] = NowLayer++;
    }

    mflog::log.set_all_processors_out();

    while (1)
    {
        LastLayer = NowLayer - 1;
        status1 = 0;
        count = 2 * nBFace;
        for (i = nBFace; i < nTFace; i++)
        {
            c1 = f2c[count++];
            c2 = f2c[count++];

            if (Layer[c1] == LastLayer && Layer[c2] == BigLayer)
            {
                Layer[c2] = NowLayer;
                status1++;
            }
            else if (Layer[c2] == LastLayer && Layer[c1] == BigLayer)
            {
                Layer[c1] = NowLayer;
                status1++;
            }
        }

        // 检查是否有遗漏的单元
        status2 = 0;
        if (status1 == 0)
        {
            for (i = 0; i < nTCell; i++)
            {
                if (Layer[i] == BigLayer)
                {
                    status2++;
                    break;
                }
            }
        }

        // 如果没有，就跳出循环
        if (status1 == 0 && status2 == 0)
            break;

        // 该块存在多区不联通，需要再次启动
        if (status1 == 0 && status2 != 0)
        {
            status3 = 0;
            // 找第一个物面单元作为起始层
            for (i = 0; i < nBFace; i++)
            {
                type = bcr[i]->GetType();
                if (type == WALL)
                {
                    c1 = f2c[i + i];
                    if (Layer[c1] == BigLayer)
                    {
                        Layer[c1] = 0;
                        status3++;
                        break;
                    }
                }
            }
            // 如果没有物面，则找任意非内部边界单元作为起始层
            if (status3 == 0)
            {
                for (i = 0; i < nBFace; i++)
                {
                    type = bcr[i]->GetType();
                    if (type != INTERFACE)
                    {
                        c1 = f2c[i + i];
                        if (Layer[c1] == BigLayer)
                        {
                            Layer[c1] = 0;
                            status3++;
                            break;
                        }
                    }
                }
            }
            // 如果全是内部边界，则找第一个边界单元作为起始层
            if (status3 == 0)
            {
                for (i = 0; i < nBFace; i++)
                {
                    c1 = f2c[i + i];
                    if (Layer[c1] == BigLayer)
                    {
                        Layer[c1] = 0;
                        status3++;
                        break;
                    }
                }
            }

            NowLayer = 0;
            // 肯定有边界单元存在，否则报错退出
            if (status3 == 0)
            {
#ifdef MF_MPICH
                mflog::log << endl
                           << "Error! Zone " << myZone << " is a multi-block, but I do not find"
                           << " the boundary cell for layer!" << endl;
#else
                mflog::log << endl
                           << "Error! This is a multi-block, but I do not find"
                           << " the boundary cell for layer!" << endl;
#endif
                mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
            }
        }

        NowLayer++;
    }

    for (i = 0; i < nTCell; i++)
    {
        MaxLayer = MAX(MaxLayer, Layer[i]);
    }
    MaxLayer++;

    // 迭代顺序从低层到高层
    count = 0;
    for (NowLayer = 0; NowLayer < MaxLayer; NowLayer++)
    {
        for (i = 0; i < nTCell; i++)
        {
            if (Layer[i] == NowLayer)
            {
                LUSGSCellOrder[count++] = i;
            }
        }
    }
    for (i = 0; i < nTCell; i++)
    {
        Layer[LUSGSCellOrder[i]] = i;
    }
    for (i = nTCell; i < n; i++)
    {
        Layer[i] = i;
    }
    printf("LUSGS reordered by hyperplane approach 2!\n");
    this->UpdateDataPtr(LUSGSCellOrder, INT, nTCell, "LUSGSCellOrder");
    this->UpdateDataPtr(Layer, INT, n, "LUSGSLayer");
}

//   网格重排序
//   参考文献：航空计算技术35卷第3期，基于非结构网格流场计算的网格重排序，方法2
//                                 Reordering of 3-D Hybrid Grids for Vectorized LU-SGS NS Computations
void PolyGrid::ReorderCellforLUSGS_3()
{
    IntType i, j, c1, c2, type, count, status1, status2, status3;
    IntType n = nBFace + nTCell;
    IntType BigLayer = 1000 * nTCell, NowLayer = 0, LastLayer, MaxLayer = 0;

    IntType *Layer = NULL;
    IntType *LUSGSCellOrder = NULL;
    snew_array_1D(Layer, n);
    for (i = 0; i < n; i++)
        Layer[i] = BigLayer;
    snew_array_1D(LUSGSCellOrder, nTCell);
    for (i = 0; i < nTCell; i++)
        LUSGSCellOrder[i] = -1;
    // 找任意一个物面单元作为起始层
    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        if (type == WALL)
        {
            c1 = f2c[i + i];
            Layer[c1] = NowLayer++;
            break;
        }
    }
    // 如果没有物面，则找第一个边界单元作为起始层
    if (NowLayer == 0)
    {
        c1 = f2c[0 + 0];
        Layer[c1] = NowLayer++;
    }

    mflog::log.set_all_processors_out();

    while (1)
    {
        LastLayer = NowLayer - 1;
        status1 = 0;
        count = 2 * nBFace;
        for (i = nBFace; i < nTFace; i++)
        {
            c1 = f2c[count++];
            c2 = f2c[count++];

            if (Layer[c1] == LastLayer && Layer[c2] == BigLayer)
            {
                Layer[c2] = NowLayer;
                status1++;
            }
            else if (Layer[c2] == LastLayer && Layer[c1] == BigLayer)
            {
                Layer[c1] = NowLayer;
                status1++;
            }
        }

        // 检查是否有遗漏的单元
        status2 = 0;
        if (status1 == 0)
        {
            for (i = 0; i < nTCell; i++)
            {
                if (Layer[i] == BigLayer)
                {
                    status2++;
                    break;
                }
            }
        }

        // 如果没有，就跳出循环
        if (status1 == 0 && status2 == 0)
            break;

        // 该块存在多区不联通，需要再次启动
        if (status1 == 0 && status2 != 0)
        {
            status3 = 0;
            // 找任意一个物面单元作为起始层
            for (i = 0; i < nBFace; i++)
            {
                type = bcr[i]->GetType();
                if (type == WALL)
                {
                    c1 = f2c[i + i];
                    if (Layer[c1] == BigLayer)
                    {
                        Layer[c1] = 0;
                        status3++;
                        break;
                    }
                }
            }
            // 如果没有物面，则找任意边界单元作为起始层
            if (status3 == 0)
            {
                for (i = 0; i < nBFace; i++)
                {
                    c1 = f2c[i + i];
                    if (Layer[c1] == BigLayer)
                    {
                        Layer[c1] = 0;
                        status3++;
                        break;
                    }
                }
            }

            NowLayer = 0;
            // 肯定有边界单元存在，否则报错退出
            if (status3 == 0)
            {
#ifdef MF_MPICH
                mflog::log << endl
                           << "Error! Zone " << myZone << " is a multi-block, but I do not find"
                           << " the boundary cell for layer!" << endl;
#else
                mflog::log << endl
                           << "Error! This is a multi-block, but I do not find"
                           << " the boundary cell for layer!" << endl;
#endif
                mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
            }
        }

        NowLayer++;
    }
    for (i = 0; i < nTCell; i++)
    {
        MaxLayer = MAX(MaxLayer, Layer[i]);
    }
    MaxLayer++;

    // 排列子层
    IntType *SubLayerNum = NULL;
    IntType *SubLayer = NULL;
    snew_array_1D(SubLayerNum, MaxLayer);
    snew_array_1D(SubLayer, nTCell);
    for (i = 0; i < nTCell; i++)
    {
        SubLayer[i] = 0;
    }

    for (NowLayer = 0; NowLayer < MaxLayer; NowLayer++)
    {
        SubLayerNum[NowLayer] = 1;
        while (1)
        {
            status1 = 0;
            for (i = nBFace; i < nTFace; i++)
            {
                c1 = f2c[i + i];
                c2 = f2c[i + i + 1];

                if (Layer[c1] == NowLayer && Layer[c1] == Layer[c2])
                {
                    if (SubLayer[c1] == SubLayer[c2] && SubLayer[c1] == SubLayerNum[NowLayer] - 1)
                    {
                        SubLayer[c2] = SubLayerNum[NowLayer];
                        status1++;
                    }
                }
            }
            if (status1 == 0)
                break;
            SubLayerNum[NowLayer]++;
        }
    }

    count = 0;
    for (NowLayer = 0; NowLayer < MaxLayer; NowLayer++)
    {
        count += SubLayerNum[NowLayer];
    }
    count--;
    for (NowLayer = MaxLayer - 1; NowLayer >= 0; NowLayer--)
    {
        for (i = SubLayerNum[NowLayer] - 1; i >= 0; i--)
        {
            for (j = 0; j < nTCell; j++)
            {
                if (Layer[j] == NowLayer && SubLayer[j] == i)
                {
                    Layer[j] = count;
                }
            }
            count--;
        }
    }
    sdel_array_1D(SubLayerNum);
    sdel_array_1D(SubLayer);
    for (i = 0; i < nTCell; i++)
    {
        MaxLayer = MAX(MaxLayer, Layer[i]);
    }
    MaxLayer++;

    // 虚拟单元排在最后一层
    for (i = 0; i < nBFace; i++)
    {
        c2 = f2c[i + i + 1];
        Layer[c2] = MaxLayer;
    }

    IntType *cellsPerlayer = NULL;
    snew_array_1D(cellsPerlayer, nTCell);
    cellsPerlayer[0] = MaxLayer;
    cellsPerlayer[1] = 0;

    IntType temp = 0;
    // 迭代顺序从低层到高层
    count = 0;
    for (NowLayer = 0; NowLayer < MaxLayer; NowLayer++)
    {
        for (i = 0; i < nTCell; i++)
        {
            if (Layer[i] == NowLayer)
            {
                LUSGSCellOrder[count++] = i;
            }
        }
        cellsPerlayer[NowLayer + 2] = count;
    }
    for (i = 0; i < nTCell; i++)
    {
        Layer[LUSGSCellOrder[i]] = i;
    }
    for (i = nTCell; i < n; i++)
    {
        Layer[i] = i;
    }
    printf("LUSGS reordered by hyperplane approach 3!\n");
    this->UpdateDataPtr(cellsPerlayer, INT, nTCell, "LUSGScellsPerlayer");
    this->UpdateDataPtr(LUSGSCellOrder, INT, nTCell, "LUSGSCellOrder");
    this->UpdateDataPtr(Layer, INT, n, "LUSGSLayer");
}
#endif
*/

#ifdef OMP_CellColor
/*!
 * @brief       网格重排序
 * @param       colort      0 or 1 to be chosen where 0 is imbalance greedy algorithm, 1 is balance algorithm
 * @return      LUSGSLayer, LUSGSCellOrder, LUSGScellsPerlayer
 * @note        参考文献： Reordering of 3-D Hybrid Grids for Vectorized LU-SGS NS Computations
 * @note        改进了子层排序算法以提升子层均衡性，提高 OpenMP 的算法均衡性
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-06
 */
void PolyGrid::LUSGSGridColor(IntType colort)
{
    IntType i, nf, j, maxcolor, colored, tempcolor, tempcolorNum, c, c1, c2;
    IntType count;
    IntType n = nTCell + nBFace;
    snew_array_1D(LUSGSLayer, n);              ///< the array to judge the upper or lower cell in LUSGS
    snew_array_1D(LUSGSCellOrder, nTCell);     ///< the coloring information array
    snew_array_1D(LUSGScellsPerlayer, nTCell); ///< the number of colors and the indce number in each color

    for (i = 0; i < n; i++)
    {
        LUSGSLayer[i] = -1;
    }
    maxcolor = 1;
    tempcolorNum = 1;
    tempcolor = 0;
    if (nCPC == 0)
        nCPC = CalnCPC(this);
    if (c2c == 0)
        c2c = CalC2C(this);

    for (i = 0; i < nTCell; i++)
    { // Get grid's maxcolor, usually the degree of the grid, N or N+1.
        maxcolor = max(maxcolor, nCPC[i]);
    }
    maxcolor++;
    // maxcolor = 10;
    LUSGScellsPerlayer[0] = maxcolor;
    for (i = 1; i < maxcolor + 3; i++)
    {
        LUSGScellsPerlayer[i] = 0;
    }

    for (i = 0; i < nBFace; i++)
    { // 对包含壁面和边界面体的着色
        c1 = f2c[i + i];
        if (LUSGSLayer[c1] >= 0)
        {
            continue;
        }
        colored = (1 << maxcolor) - 1;
        for (j = 0; j < nCPC[c1]; j++)
        {
            c = c2c[c1][j];
            if ((c < nTCell) && (LUSGSLayer[c] >= 0))
            {
                colored = colored & (~(1 << LUSGSLayer[c]));
            }
        }
        if (colored == 0)
        {
            colored = 1 << maxcolor;
            maxcolor++;
        }
        tempcolorNum = nTCell;
        if (colort == 0)
        { // imbalanced greedy algorithm for face coloring
            for (j = 0; j < maxcolor; j++)
            {
                if (colored & (1 << j))
                { // with satisfied color number
                    tempcolor = j;
                    break;
                }
            }
        }
        else
        { // balanced algorithm for face coloring
            for (j = 0; j < maxcolor; j++)
            {
                if (colored & (1 << j))
                { // with satisfied color number
                    if (LUSGScellsPerlayer[j + 2] < tempcolorNum)
                    {
                        tempcolorNum = LUSGScellsPerlayer[j + 2];
                        tempcolor = j;
                    }
                }
            }
        }
        LUSGScellsPerlayer[tempcolor + 2]++;
        LUSGSLayer[c1] = tempcolor;
    }
    for (c = 0; c < nTCell; c++)
    { // 对内部体进行着色
        if (LUSGSLayer[c] >= 0)
        {
            continue;
        }
        colored = (1 << maxcolor) - 1;
        for (j = 0; j < nCPC[c]; j++)
        {
            c1 = c2c[c][j];
            if ((c1 < nTCell) && (LUSGSLayer[c1] >= 0))
            {
                colored = colored & (~(1 << LUSGSLayer[c1]));
            }
        }
        if (colored == 0)
        {
            colored = 1 << maxcolor;
            maxcolor++;
        }
        tempcolorNum = nTCell;
        if (colort == 0)
        { // imbalanced greedy algorithm for face coloring
            for (j = 0; j < maxcolor; j++)
            {
                if (colored & (1 << j))
                { // with satisfied color number
                    tempcolor = j;
                    break;
                }
            }
        }
        else
        { // balanced algorithm for face coloring
            for (j = 0; j < maxcolor; j++)
            {
                if (colored & (1 << j))
                { // with satisfied color number
                    if (LUSGScellsPerlayer[j + 2] < tempcolorNum)
                    {
                        tempcolorNum = LUSGScellsPerlayer[j + 2];
                        tempcolor = j;
                    }
                }
            }
        }
        LUSGScellsPerlayer[tempcolor + 2]++;
        LUSGSLayer[c] = tempcolor;
    }

    // 虚拟单元排在最后一层
    maxcolor++;
    for (i = 0; i < nBFace; i++)
    {
        c2 = f2c[i + i + 1];
        LUSGSLayer[c2] = maxcolor;
    }
    count = 0;

    for (i = 0; i < maxcolor; i++)
    {
        for (j = 0; j < nTCell; j++)
        {
            if (LUSGSLayer[j] == i)
            {
                LUSGSCellOrder[count++] = j;
            }
        }
        LUSGScellsPerlayer[i + 2] += LUSGScellsPerlayer[i + 1];
    }

    for (i = 0; i < nTCell; i++)
    {
        LUSGSLayer[LUSGSCellOrder[i]] = i;
    }
    for (i = nTCell; i < n; i++)
    {
        LUSGSLayer[i] = i;
    }
#ifdef MF_DEBUG
    cout << endl
         << "LUSGS reordered by cell color approach 4!" << endl;
#endif
    // this->UpdateDataPtr(cellsPerlayer, INT, nTCell, "LUSGScellsPerlayer");
    // this->UpdateDataPtr(LUSGSCellOrder, INT, nTCell, "LUSGSCellOrder");
    // this->UpdateDataPtr(Layer, INT, n, "LUSGSLayer");
}
#endif ///< ~OMP_CellColor

/**
#ifdef USING_PETSC
void PolyGrid::CalNonZeroInfo(IntType *&d_nnz, IntType *&o_nnz)
{
    snew_array_1D(d_nnz, nTCell);
    snew_array_1D(o_nnz, nTCell);
    for (IntType iCell = 0; iCell < nTCell; iCell++)
    {
        d_nnz[iCell] = 1; // include itself
        o_nnz[iCell] = 0;
    }
    for (IntType iFace = nBFace - nIFace; iFace < nBFace; iFace++)
    {
        IntType c1 = f2c[iFace + iFace];
        o_nnz[c1]++;
    }
    // diag elements all belong to interior grids
    for (IntType iFace = nBFace; iFace < nTFace; iFace++)
    {
        IntType c1 = f2c[iFace + iFace];
        IntType c2 = f2c[iFace + iFace + 1];
        d_nnz[c1]++;
        d_nnz[c2]++;
    }
}

IntType *PolyGrid::CalGhost2Global(IntType Bstart)
{
    if (ghost2global)
        return ghost2global;
#ifdef MF_MPICH
    snew_array_1D(ghost2global, nIFace);
    IntType **bqs = NULL;
    IntType **bqr = NULL;

    MPI_Request *req_send = NULL;
    MPI_Request *req_recv = NULL;
    MPI_Status *status_array = NULL;
    snew_array_1D(status_array, nNeighbor);
    snew_array_1D(req_send, nNeighbor);
    snew_array_1D(req_recv, nNeighbor);
    snew_array_2D(bqr, nNeighbor, nZIFace, false);
    snew_array_2D(bqs, nNeighbor, nZIFace, false);

    for (IntType i = 0; i < nNeighbor; i++)
    {
        for (IntType j = 0; j < nZIFace[i]; j++)
        {
            // Global_index = Local_index + Block_start
            // bqs[i][j] = bCNo[i][j] + Bstart;
            bqs[i][j] = bCNo[i][j] + Bstart;
        }
    }
    for (IntType np = 1; np <= numprocs; np++)
    {
        if (myZone == np)
        {
            //       send   to other
            for (IntType g = 0; g < nNeighbor; g++)
            {
                IntType nbZone = nb[g];
                MPI_Isend(bqs[g], nZIFace[g], MPI_INT, nbZone, level, MPI_COMM_WORLD, &req_send[g]);
            }
        }
        else
        {
            //       receive   from np
            for (IntType g = 0; g < nNeighbor; g++)
            {
                IntType nbZone = nb[g];
                if (nbZone == (np - 1))
                {
                    MPI_Irecv(bqr[g], nZIFace[g], MPI_INT, nbZone, level, MPI_COMM_WORLD, &req_recv[g]);
                }
            }
        }
    } // np
    MPI_Waitall(nNeighbor, req_recv, status_array);
    IntType ifStart = nBFace - nIFace;
    for (IntType g = 0; g < nNeighbor; g++)
    {
        for (IntType i = 0; i < nZIFace[g]; i++)
        {
            IntType idx = bFNo[g][i] - ifStart;
            ghost2global[idx] = bqr[g][i];
        }
    }
    MPI_Waitall(nNeighbor, req_send, status_array);
    sdel_array_2D(bqr, nNeighbor, false);
    sdel_array_2D(bqs, nNeighbor, false);
    sdel_array_1D(req_recv);
    sdel_array_1D(req_send);
    sdel_array_1D(status_array);
    return ghost2global;
#endif
}

#endif

*/
