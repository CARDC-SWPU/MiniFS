/*!
 * @file        grid_polyhedra.h
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
 * <tr><td> 2023-06-23  <td> 2.0      <td> Wisces  <td> 添加了 MPI 并行（使用条件编译）
 * <tr><td> 2023-07-03  <td> 3.0      <td> Wisces  <td> 添加了 OpenMP 并行（使用条件编译）
 * </table>
 */

#ifndef MF_GRID_POLYHEDRA_H
#define MF_GRID_POLYHEDRA_H

//!< C++ build-in head files
#include <iostream>
#include <vector>
// #include <cstdlib>

//!< user defined head files
#include "number_type.h"
#include "constant.h"
#include "boundary_condition.h"
#include "grid_base.h"
#include "memory_util.h"

//!< head file relying on condition-compiling
#ifdef MF_MPICH
#include <mpi.h>
#endif

#ifdef GASPI
#include "GASPI.h"
#include "omp.h"
#include "GASPI_utility.h"
#endif

#ifdef OMP_DIVREP
#include "metis.h"
#endif

#ifdef TDTREE
#include "TDTree.h"
#include <unordered_map>
#endif

// #ifdef SCOTCH
// #include <stdint.h>
// #include <scotch.h>
// #endif // scotch

using namespace std;

/*!
 * @brief       Class PolyGrid
 *
 */
class PolyGrid : public Grid
{
public:
        //!< for face/group colouring
        // #ifdef OMP_FaceColor
#if (defined OMP_FaceColor) || (defined OMP_GroupColor)
        vector<IntType> bfacegroup; ///< 对于物理边界面，每组颜色（除第一组直接是从 0 开始外）中第一个 face 的编号
        vector<IntType> ifacegroup; ///< 对于内部面，每组颜色（除第一组直接是从 0 开始外）中第一个 face 的编号
#endif

        //!< for group coloring
#ifdef OMP_GroupColor
        bool GroupColorSuccess;
        const static IntType groupSize = 128;
        const static IntType MaxColors = 70; ///< 36 - for m6wing_12w
                                             ///< 63 - for m6wing_25w
                                             ///< 70 - for CHNT1_600w(min 69)
        // IntType *id_bface; ///< new to old
        // IntType *id_iface; ///< new to old
#endif

        //!< for Divide&Replicate
#ifdef OMP_DIVREP
        IntType *idx_pthreads_bface;
        IntType *id_division_bface;
        IntType *idx_pthreads_iface;
        IntType *id_division_iface;
        IntType threads;
        bool DivRepSuccess;
        // IntType *endIndex_bFace_vec; ///< Used for BoundedColoring
        // IntType *endIndex_iFace_vec; ///< Used for BoundedColoring
#endif

        // #ifdef TDTREE
        // #ifndef FORKJOIN
        //         const static IntType nbParts = 4;
        // #else
        //         const static IntType nbParts = 24 * 8;
        // #endif
        //         const static IntType partSize = 512;
        // #endif

        /**
        #if (defined FS_SIMD) && (defined FS_SIMD_AVX) && (defined OMP_FaceColor)
                // for OMP_FaceColor SIMD, based SSE, ruitian, 2021.12.27
                RealGeom *xfntile, *yfntile, *zfntile;
                RealGeom *areatile, *qsumtile;
                RealGeom *dqdxtile, *dqdytile, *dqdztile;
                IntType *f2c1, *f2c2;
        #elif (defined FS_SIMD) && (defined Tile)
                // for Tile SIMD, based SSE, ruitian, 2021.12.27
                // for boundary faces, not including nIFace
                IntType *bfacerow;
                IntType *bfacecol;
                IntType *bfaceval;
                IntType bfacennz;
                // for interior faces
                IntType *ifacerow;
                IntType *ifacecol;
                IntType *ifaceval;
                IntType ifacennz;
                // tile for SIMD, ruitian, 2021.12.19, Multi-thread and Multi-lane tiling
                // for boundary faces, not including nIFace, Multi-thread and Multi-lane execution
                IntType *bSIMDrow;
                IntType *bSIMDcol;
                IntType *bSIMDval;
                IntType bSIMDnnz;
                IntType *bfacezero;
                vector<vector<int>> *boffsets;
                // for interior faces, nTface-nBFace
                IntType *iSIMDrow;
                IntType *iSIMDcol;
                IntType *iSIMDval;
                IntType iSIMDnnz;
                IntType *ifacezero;
                std::vector<vector<int>> *ioffsets;

                RealGeom *qsumtile;
                RealGeom *dqdxtile, *dqdytile, *dqdztile;
                RealGeom *xfnt, *yfnt, *zfnt;
                RealGeom *areat, *qsumt;
                IntType *nNPFt;
                // add by ruitian
        #endif
        */

private:
#ifdef TDTREE
        TDTree *TDTreeRoot = NULL;
#endif
#ifdef GASPI
    //vector<pair<int,int>> *f2b; 

    pair<int,int> **C2B, **N2B;
    IntType *nC2B, *nN2B, maxLength; 
    IntType **Ghost = NULL;

    IntType *GASPI_mapping, *GASPI_mappingN;
    IntType GASPI_blockSize, GASPI_nbBlock, GASPI_notifications, *GASPI_notifications_node, *GASPI_recvNotifications;

    gaspi_segment_id_t GASPI_sendDataSegmentId, GASPI_recvDataSegmentId, GASPI_sendOffsetSegmentId, GASPI_recvOffsetSegmentId;
    gaspi_pointer_t GASPI_sendDataPtr, GASPI_recvDataPtr, GASPI_sendOffsetPtr, GASPI_recvOffsetPtr;
    gaspi_rank_t GASPI_iProc, GASPI_nProc;
    IntType *GASPI_offset, *GASPI_commed;
    IntType *GASPI_need_Comm;

    RealFlow *GASPI_sendData, *GASPI_recvData;
    IntType *GASPI_sendOffset, *GASPI_recvOffset;

    omp_lock_t *GASPI_lock;

#endif
        //!< flow variables, dt_timestep, residuals, DQ - gField
        //!< add by Wisces
        //!< 网格重排序
        IntType *LUSGSLayer;
        IntType *LUSGSCellOrder;
        IntType *LUSGScellsPerlayer;

        //!< flow variables
        RealFlow *rho;
        RealFlow *u;
        RealFlow *v;
        RealFlow *w;
        RealFlow *p;

        //!< gradient
        RealFlow *dqdx;
        RealFlow *dqdy;
        RealFlow *dqdz;

        //!< time step
        RealFlow *dt; ///< dt_timestep

        //!< residuals
        RealFlow *res;

        //!< DQ
        RealFlow *DQ;

        //!< Some parameters for the Gradient computations
        //!< add by Wisces
        // IntType GaussLayer;
        // IntType *CellLayerNo;
        IntType *nodesymm;
        RealGeom *xfn_n_symm;
        RealGeom *yfn_n_symm;
        RealGeom *zfn_n_symm;

        /*!
         * @brief       kernel grid information for face based FVM method
         *
         */
        IntType level;          // multi-grid level, level 0 being the finest
        IntType nTFace, nTCell; // @readfile no. of total faces, no of total cells
        IntType nBFace;         // @readfile boundary faces which include interfaces
        IntType nIFace;         // parallel interface faces(zero if serial program)
        IntType nINode;         // parallel interface nodes(zero if serial program)

        IntType *nNPF, *f2n, **F2N; // @readfile no. of node per face and face->node connectivity
        IntType *f2c;               // @readfile face->cell connectivity

        RealGeom *xcc, *ycc, *zcc; // cell center data, including ghosts
        RealGeom *xfc, *yfc, *zfc; // face center data;
        RealGeom *xfn, *yfn, *zfn; // face unit normal;
        RealGeom *area, *vol;      // face area, cell volume excluding ghosts;

        BCRecord **bcr; // associate each boundary. face with a record

        // Grid *cGrid, *fGrid; // coarse or finer grid
        IntType *c2cc; ///< Used for TDTree.
        // IntType *nCCPN, **N2CC;
        // IntType *nVPN, **N2V;
        // RealGeom *weight_prol;

        // Auxiliary connectivities
        IntType *nCPC, **c2c; ///< Used for CellColoring. number of neighbor cells per cell, and cell->cell.
        IntType *nNPC, **C2N; ///< Used for gradient calculation / Used for Divide&Conquer
        IntType *nFPC, **C2F; ///< number of faces in each real cell, real cell to face connection
        IntType *nCPN, **N2C; ///< Used for gradient calculation. number of cells in each node(not include all ghost cell), and node->cell

        // Neighbor information for zones
        IntType nNeighbor;  // no of zonal neighbors linked by face
        IntType nNeighborN; ///< Used for gradient calculation. no of zonal neighbors linked by node

        // neighbor information of Zones
        IntType *nb;  // the neighbor zone/rank number for face
        IntType *nbN; ///< Used for gradient calculation. the neighbor zone/rank number for node
        // PolyGrid **nbg; // ??? not used

        // parallel information for each interface (saved in grid file of 'mmgrid' type)
        IntType *nbZ, *nbBF;         // parallel information for interfaces: zone No. and face No. in that zone. Size is nIFace.
        IntType *nbSN, *nbZN, *nbRN; // parallel information for inter-nodes: nodes' id on current zone, zone No. and nodes' id on other zone. Size is nINode.
        // Actual parallel information used to communication through MPI.
        IntType *nZIFace; // No. of interfaces for each neighbor
        IntType **bCNo;   // cell number in the current zone for each neighbor, each interface, ordered
                          // according to the current zone. Used to send data on cell center to other processor.
        IntType **bFNo;   // boundary face number in the current zone for each neighbor, each interface, ordered
                          // according to the neighbor zone. Used to receive data on cell center from other processor.
        IntType *nZINode; ///< Used for gradient calculation. No. of inter-nodes for each neighbor
        IntType **bNSNo;  ///< Used for gradient calculation. node number for sending data on nodes for each neighbor, each node.
        IntType **bNRNo;  ///< Used for gradient calculation. node number for receiving data on nodes for each neighbor, each node.

        // // Active and interpolation information for overlap grid
        // IntType *node_act, *edge_act;          // the active mark of nodes, cells, faces and edges.
        // IntType *nChCellS;                     // No. of cells sending to other processors whose grid is different.
        // IntType *nChCellR;                     // No. of cells receiving from other processors whose grid is different.
        // IntType **nChSNo, **nChRNo;            // The cells id to send and receive
        // RealGeom **nChxcc, **nChycc, **nChzcc; // Centroid of cells to interpolation on other processors
        // IntType *nChRNoComp;                   // Cells id which have no donor-cells on current processor

        // // Zone id for all cells used to initialize with different value for different zone.
        // IntType *Cell2Zone;

        // // grid quality
        // RealGeom *facecentroidskewness;
        // RealGeom *faceangleskewness;
        // RealGeom *cellcentroidskewness;
        // RealGeom *cellvolsmoothness;
        // RealGeom *faceangle;
        // IntType *cellwallnumber;

        // Used for gradient calculation. flow reconstruction from cell to node
        IntType *Nmark;                                             // Node type
        RealGeom *WeightNodeDist, **WeightNodeC2N, **WeightNodeN2C; // weight of node-distance
        RealGeom **WeightNodeBFace2C;                               // weight of node-distance for BFace to related Cells

        // // weight for prolongation from coarse to fine grid
        // RealGeom *WeightNodeProl;

        // // Central Scheme using the face order, 0 for the first cell, 1 for the next cell
        // RealGeom **Seta_Center;

        // moving grid
        RealGeom *vgn;                            // normal velocity for all faces
        RealGeom *BFacevgx, *BFacevgy, *BFacevgz; // 3 components of velocity for boundary faces.
        RealGeom *vccx, *vccy, *vccz;             // velocity at the cell center, for preconditioning of unsteady flow.

        // additional grid information for efficiency
        // RealGeom VolAvg; // cell's average volume, using for venkatakrishnan's limiter

        // #ifdef USING_PETSC
        //         IntType *ghost2global;
        // #endif

        // // data for line implicit
        // IntType nLinesForLI;            // line number
        // IntType *nCellsInLineForLI;     // cell number of each line
        // IntType **CellsOfLinesForLI;    // cell index of each line's cells
        // IntType **FacesOfLinesForLI;    // face index of each line's faces
        // IntType *LineIndexOfCellsForLI; // line index of each cell

        // // çŠé˘ĺşĺççşżĺćć°çŽĺçşżé(ĺ¨čżćĽäšĺ,ç¨äşCFLć°ćžĺ¤§ä˝żç
        // IntType n_lines_wall_start, n_lines_wall;
        // IntType *n_cells_in_line_wall;

        // RealGeom span_len_of_2d;

        // IntType *order_cell_oTon; // dingxin:the order of cell,old->new
        // IntType *order_cell_nToo; // dingxin:the order of cell,new->old

        // // add by dingxin
        // IntType *pfaceOrder_reorder; // new->old
        // IntType *pfaceOrder_color;

public:
#ifdef TDTREE
        void SetTDTree(TDTree *treeRoot);
        TDTree *GetTDTree() const;
#endif
#ifdef GASPI
    void checkMessage(IntType nVar, RealFlow **q);

    void GASPI_wait_for_queue_half_full (gaspi_queue_id_t queueID);

    void GASPI_wait_for_flush_queues ();
    
    void GASPIUnzip(IntType nVar, RealFlow **q, IntType neighbor, IntType ns, IntType ne);

    void GASPIWaitAll(IntType nVar, RealFlow **q);

    void GASPIComm(RealFlow **q, int *boundryCell, int boundryNbCell);

    void GASPIComm(RealFlow **q, int *n, int *Dcell);

    void GASPIComm(IntType nVar, RealFlow **q, IntType neighbor, int block);

    void GASPIComm(IntType nVar, RealFlow **q, IntType neighbor, int block, int No);
    
    void GASPIUnzip_Node(IntType nVar, RealFlow **q, IntType neighbor, IntType ns, IntType ne);

    void GASPIWaitAll_Node(IntType nVar, RealFlow **q);

    void GASPIComm_Node(IntType nVar, RealFlow **q, IntType neighbor, int block);

    void GASPIPrintf();
    
    void GASPIRecvSendVarNeighbor_Togeth(IntType nvar, RealFlow **q);

    void GASPIRecvSendVarNeighbor_Over(IntType nvar);

    inline pair<int,int> **GetC2B();

    inline IntType *GetnC2B();

    inline pair<int,int> **GetN2B();    
    
    inline IntType *GetnN2B();    

    //inline vector<pair<int,int>> *Getf2b();

    inline IntType GetblockSize();

    inline IntType **GetGhost();
#endif

        //!< flow variables, dt_timestep, residuals, DQ - gField
        //!< add by Wisces
        void SetLUSGSLayer(IntType *layer);
        void SetLUSGSCellOrder(IntType *luorder);
        void SetLUSGScellsPerlayer(IntType *cellsPerlayer);
        void SetRho(RealFlow *q);
        void SetU(RealFlow *q);
        void SetV(RealFlow *q);
        void SetW(RealFlow *q);
        void SetP(RealFlow *q);
        void SetDqdx(RealFlow *dq);
        void SetDqdy(RealFlow *dq);
        void SetDqdz(RealFlow *dq);
        void SetDt(RealFlow *q);
        void SetRes(RealFlow *q);
        void SetDQ(RealFlow *q);

        IntType *GetLUSGSLayer() const;
        IntType *GetLUSGSCellOrder() const;
        IntType *GetLUSGScellsPerlayer() const;
        RealFlow *GetRho() const;
        RealFlow *GetU() const;
        RealFlow *GetV() const;
        RealFlow *GetW() const;
        RealFlow *GetP() const;
        RealFlow *GetDqdx() const;
        RealFlow *GetDqdy() const;
        RealFlow *GetDqdz() const;
        RealFlow *GetDt() const;
        RealFlow *GetRes() const;
        RealFlow *GetDQ() const;

        //!< Some parameters for the Gradient computations
        //!< add by Wisces
        void SetNodeSymm(IntType *symm);
        void SetXfnNSymm(RealGeom *symm);
        void SetYfnNSymm(RealGeom *symm);
        void SetZfnNSymm(RealGeom *symm);

        IntType *GetNodeSymm() const;
        RealGeom *GetXfnNSymm() const;
        RealGeom *GetYfnNSymm() const;
        RealGeom *GetZfnNSymm() const;

        void SetNTFace(const IntType in);
        void SetNTCell(const IntType in);
        void SetNBFace(const IntType in);
        void SetNIFace(const IntType in);
        void SetNINode(const IntType in);

        IntType GetNTFace() const;
        IntType GetNTCell() const;
        IntType GetNBFace() const;
        IntType GetNIFace() const;
        IntType GetNINode() const; ///< Used for gradient calculation

        void SetnNPF(IntType *in);
        void SetnFPC(IntType *in);
        void SetnNPC(IntType *in); ///< Used for gradient calculation / Used for Divide&Conquer
        void Setf2n(IntType *in);
        void Setf2c(IntType *in);

        // Get and Set the cell center information
        RealGeom *GetXcc() const;
        RealGeom *GetYcc() const;
        RealGeom *GetZcc() const;

        // Get and set cell volume
        RealGeom *GetCellVol() const;

        // Get and Set the face center information
        RealGeom *GetXfc() const;
        RealGeom *GetYfc() const;
        RealGeom *GetZfc() const;

        // Get and Set the face normal information
        RealGeom *GetXfn() const;
        RealGeom *GetYfn() const;
        RealGeom *GetZfn() const;

        // Get and set face area
        RealGeom *GetFaceArea() const;

        // ĺ¤éç˝ć ź
        // void SetcGrid(PolyGrid *grid);
        // void SetfGrid(PolyGrid *grid);
        // Grid *GetcGrid() const;
        // Grid *GetfGrid() const;

        void Setc2cc(IntType *in);
        IntType *Getc2cc() const; ///< Used for TDTree.
        // void SetnCCPN(IntType *in);
        // IntType *GetnCCPN() const;
        // void SetN2CC(IntType **in);
        // IntType **GetN2CC() const;
        // void SetnVPN(IntType *in);
        // IntType *GetnVPN() const;
        // void SetN2V(IntType **in);
        // IntType **GetN2V() const;
        // void SetWeightNodeProl(RealGeom *weight_prol);
        // RealGeom *GetWeightNodeProl(void) const;

        void SetnCPC(IntType *in); ///< Used for CellColoring
        void SetnCPN(IntType *in); ///< Used for gradient calculation

        void Setc2c(IntType **in); ///< Used for CellColoring
        void SetC2N(IntType **in); ///< Used for gradient calculation / Used for Divide&Conquer
        void SetC2F(IntType **in);
        void SetN2C(IntType **in); ///< Used for gradient calculation
        // F2N is special and just a reference to f2n, so use 1D array operator, tangj
        void SetF2N(IntType **in);

        IntType *GetnNPF() const;
        IntType *GetnFPC() const;
        IntType *GetnNPC() const; ///< Used for gradient calculation / Used for Divide&Conquer
        IntType *Getf2n() const;
        IntType *Getf2c() const;

        IntType *GetnCPC() const; ///< Used for CellColoring
        IntType **Getc2c() const; ///< Used for CellColoring
        IntType **GetC2N() const; ///< Used for gradient calculation / Used for Divide&Conquer
        IntType **GetC2F() const;
        IntType **GetF2N() const;
        IntType **GetN2C() const; ///< Used for gradient calculation
        IntType *GetnCPN() const; ///< Used for gradient calculation

        void SetNumberOfFaceNeighbors(const IntType n);
        IntType GetNumberOfFaceNeighbors() const;
        void SetNumberOfNodeNeighbors(const IntType n); ///< Used for gradient calculation
        IntType GetNumberOfNodeNeighbors() const;       ///< Used for gradient calculation

        void SetFaceNeighborZones(IntType *fnz);
        IntType *GetFaceNeighborZones() const;
        void SetNodeNeighborZones(IntType *nnz); ///< Used for gradient calculation
        // IntType *GetNodeNeighborZones() const;
        // void SetNeighborGrids(PolyGrid **grids);

        void SetnbZ(IntType *in);
        IntType *GetnbZ() const;
        void SetnbBF(IntType *in);
        IntType *GetnbBF() const;
        void SetnbSN(IntType *in);
        IntType *GetnbSN() const;
        void SetnbZN(IntType *in);
        IntType *GetnbZN() const; ///< Used for gradient calculation
        void SetnbRN(IntType *in);
        IntType *GetnbRN() const;
        // IntType *Getface_act() const;

        // void Setcoe_trans(RealGeom *in);
        // RealGeom *Getcoe_trans() const;

        // IntType GetLevel() const;
        // void SetLevel(IntType lin);

        void Setbcr(BCRecord **bcrs);
        BCRecord **Getbcr() const;

        // void SetVolAvg(RealGeom in);
        // RealGeom GetVolAvg() const;

        // RealGeom *GetGridQualityFaceCentroidSkewness(void) const;

        // IntType *GetGridQualityCellWallNumber(void) const;

        //!< Used for gradient calculation. flow reconstruction from cell to node
        void SetNodeType(IntType *node_type);                    ///< Used for gradient calculation
        IntType *GetNodeType(void) const;                        ///< Used for gradient calculation
        void SetWeightNodeDist(RealGeom *weight);                ///< Used for gradient calculation
        RealGeom *GetWeightNodeDist(void) const;                 ///< Used for gradient calculation
        void SetWeightNodeC2N(RealGeom **weight_c2n);            ///< Used for gradient calculation
        RealGeom **GetWeightNodeC2N(void) const;                 ///< Used for gradient calculation
        void SetWeightNodeN2C(RealGeom **weight_n2c);            ///< Used for gradient calculation
        RealGeom **GetWeightNodeN2C(void) const;                 ///< Used for gradient calculation
        void SetWeightNodeBFace2C(RealGeom **WeightNodebFace2c); ///< Used for gradient calculation
        RealGeom **GetWeightNodeBFace2C(void) const;             ///< Used for gradient calculation

        // moving grid
        RealGeom *GetFaceNormalVelocity(void) const;
        RealGeom *GetBoundaryFaceVelocityX(void) const;
        RealGeom *GetBoundaryFaceVelocityY(void) const;
        RealGeom *GetBoundaryFaceVelocityZ(void) const;

        void
        ComputeMetrics();
        // void InitialVgn();
        // void ComputeDist2WallTriang(RealGeom *dist2wall_cell, IntType mark);
        // void WriteInfoDist();
        // void ComputeCellDist();

        void SetUpComm();
        void SetUpComm_Node(); ///< Used for gradient calculation
        // void CommInterfaceData(IntType zone, PolyGrid *grid, const char *name);
        // void CommCellCenterData(IntType zone, PolyGrid *grid);

        // void FindCellLayerNo();
        // void FindNormalFace();

        void quick_sort(IntType *a, IntType is, IntType ie); ///< Used for gradient calculation / Used for Divide&Conquer
        void swap(IntType *a, IntType i, IntType j);         ///< Used for gradient calculation / Used for Divide&Conquer

        //!< Reorder cell for LUSGS
        void ReorderCellforLUSGS_0();
        // /// Reroder cell for line-implicit block lusgs
        // void ReorderCellforLUSGS_101();

        // // for line implicit
        // void SetLI_nLine(const IntType in);
        // IntType GetLI_nLine() const;
        // void SetLI_nCellsInLine(IntType *in);
        // IntType *GetLI_nCellsInLine() const;
        // void SetLI_CellOfLine(IntType **in);
        // IntType **GetLI_CellOfLine() const;
        // void SetLI_FaceOfLine(IntType **in);
        // IntType **GetLI_FaceOfLine() const;
        // void SetLI_LineOfCell(IntType *in);
        // IntType *GetLI_LineOfCell() const;

        // void SetLI_nLine_wall_start(const IntType in);
        // IntType GetLI_nLine_wall_start() const;
        // void SetLI_nLine_wall(const IntType in);
        // IntType GetLI_nLine_wall() const;
        // void SetLI_nCellinLine_wall(IntType *in);
        // IntType *GetLI_nCellinLine_wall() const;

        // RealGeom &span_length_of_2d();
        // RealGeom span_length_of_2d() const;

        // // check grid quality
        // void CheckGridQuality();
        // void SkewnessSummary();
        // void EquiangleSkewnessSummary();
        // void SmoothnessSummary();
        // void FindIllWallCell();
        // void CheckSymmetryFace();
        // void CheckGridScale();
        // void FaceAngleSummary();

        // // deal bad grid for robust
        // void DealBadGrid();

        // // compute some additional info for geomety for efficiency
        // void AdditionalInfoForGeometry();
        // RealGeom CalculateVolumnAverage();
        // RealGeom *CalNormalDistanceOfC2C();

        /// Wisces: 暂时先注释掉
        // void PartitionGrids(IntType *CellToZone, IntType n_zone);
        // void Getxadjadjncy(idx_t *xadj, idx_t *adjncy, idx_t *adjwgt);                                    //.h metis.h
        // void SerialMetis(idx_t *xadj, idx_t *adjncy, idx_t *adjwgt, IntType n_zone, IntType *CellToZone); //.h metis.h
        // IntType *GetCell2Zone() const
        // {
        //         return Cell2Zone;
        // }

#ifdef OMP_CellColor
        // void ReorderCellforLUSGS_1();
        // void ReorderCellforLUSGS_2();
        // void ReorderCellforLUSGS_3();
        void LUSGSGridColor(IntType colort);
#endif

// #ifdef USING_PETSC
//         /// calcute the non-zero info for the petsc block matrix
//         void CalNonZeroInfo(IntType *&d_nnz, IntType *&o_nnz);
//         IntType *CalGhost2Global(IntType Bstart);
// #endif
#ifdef MF_MPICH
        // void CommInterfaceDataMPI(IntType *q);
        void CommInterfaceDataMPI(RealFlow *q);
        // void CommInternodeDataMPI(IntType *q);    // MPI synchronize node overlap attribute
        void CommInternodeDataMPI2(IntType *q); // Used for gradient calculation
        // void CommInternodeDataMPISUM(IntType *q); // MPI Sum node value
        void CommInternodeDataMPI(RealFlow *q); ///< Used for gradient calculation
        // void CommInternodeDataMPI(RealFlow *q, IntType key);
        void RecvSendVarNeighbor(RealFlow *q);
        // void RecvSendVarNeighbor(IntType *q);

        // // communicate active attribute of nodes
        // void RecvSendVarNeighbor_Node(IntType *q);

        //!< communicate node variable and synchronize to the minimum value with some constrains
        void RecvSendVarNeighbor_Node2(IntType *q); ///< Used for gradient calculation

        // // communicate and accumulate node variable
        // void RecvSendVarNeighbor_NodeSUM(IntType *q);
        void RecvSendVarNeighbor_Node(RealFlow *q); ///< Used for gradient calculation
        // void RecvSendVarNeighbor_Node(RealFlow *q, IntType key);

        //!< 并行传值（合并版）
        void RecvSendVarNeighbor_Togeth(IntType nvar, RealFlow **q);
        void RecvSendVarNeighbor_Over(RealFlow ***bqs, RealFlow ***bqr, MPI_Request *req_send, MPI_Request *req_recv, MPI_Status *status_array, IntType nvar);
        void Set_RecvSend(RealFlow ***qs, RealFlow ***qr, IntType nvar);
        void Add_RecvSend(RealFlow ***qs, RealFlow *q, IntType num_var);
        void Read_RecvSend(RealFlow ***qr, RealFlow *q, IntType num_var);

        // /// communicate and update the vector's ghost data, the original working vector store the inner cells variables by nTCell*nVar layout
        // void UpdateVectorGhostVar(const RealFlow *vec, RealFlow *ghosts, IntType nVar);

#endif

        // /*add by dingxin*/
        // #ifdef REORDER
        //         void Update_f2c();
        //         void Update_c2cc(PolyGrid *fgrid);
        //         void CellReordering_CMK();
        //         void CellReordering_morton();
        //         void CellReordering_metis();
        //         void CellReordering_scotch();
        //         void FaceReordering();
        //         IntType *GetNewOrder(void) const
        //         {
        //                 return order_cell_oTon;
        //         }
        // #endif // REORDER
        //         void UpdateBcr_pface();
        //         void SetPfaceOrder_reorder(IntType *in);
        //         void SetPfaceOrder_color(IntType *in);

        // explicit PolyGrid(IntType zin);
        // PolyGrid(IntType zin, IntType lin);
        PolyGrid();
        ~PolyGrid();
}; ///< ~PolyGrid

// inline functions of class PolyGrid
#include "grid_polyhedra.inl"

//=============================================================================
//                       Auxiliary Functions for PolyGrid
//=============================================================================

//!< for face colouring
#ifdef OMP_FaceColor
void FaceColoring(PolyGrid *grid, IntType bgroupsize, IntType igroupsize);
#endif

//!< for group coloring
#ifdef OMP_GroupColor
bool GroupColoring(PolyGrid *grid, bool balanceColors = false);
void quicksort_vecint(IntType s, IntType e, vector<IntType> &order, vector<IntType> &vecint);
#endif

//!< for Divide&Replicate
#ifdef OMP_DIVREP
bool DivRep(PolyGrid *grid);
#endif

//!< for TDTree
#ifdef TDTREE
void CreateTree(PolyGrid *grid);
void TDTree_CellReordering(PolyGrid *grid, PolyGrid *fgrid = NULL);
void TDTree_FaceReordering(PolyGrid *grid);
void TDTree_NodeReordering(PolyGrid *grid);
#endif

//!< Update boundary-face/interior-face conectivities after reordering.
void UpdateFaceData(PolyGrid *grid, IntType *index_face, IntType flag);

// Compute face and cell centroid by simple average all node value
void FaceCellCenterbyAverage(PolyGrid *grid, RealGeom *xfc, RealGeom *yfc, RealGeom *zfc, RealGeom *xcc, RealGeom *ycc, RealGeom *zcc);

// Compute cell center and cell volume
void CellVolCentroid(PolyGrid *grid, RealGeom *vol, RealGeom *xcc, RealGeom *ycc, RealGeom *zcc);

// Compute the normal vector, the area and face center on each cell face in 3D
void FaceAreaNormalCentroid_cycle(PolyGrid *grid, RealGeom *area, RealGeom *xfn, RealGeom *yfn, RealGeom *zfn, RealGeom *xfc, RealGeom *yfc, RealGeom *zfc);

// Correct face normal vector for face of TINY AREA
void CorrectFaceNormal(PolyGrid *grid, RealGeom *xfn, RealGeom *yfn, RealGeom *zfn);

// Correct cell centroid if it is too close to wall
void CorrectCellCentroid(PolyGrid *grid, RealGeom *xcc, RealGeom *ycc, RealGeom *zcc,
                         RealGeom *xfc, RealGeom *yfc, RealGeom *zfc,
                         RealGeom *xfn, RealGeom *yfn, RealGeom *zfn);

// Return an real number whose value is the max absolute value of two real number and
// whose sign is equal to the first param 'a'
RealGeom AbsMaxSignFirst(RealGeom a, RealGeom b);

// Compute the weight for grid nodes based on distance
void ComputeWeight3D_Node(PolyGrid *grid); ///< Used for gradient calculation

// // Check the closure of grid
// void ClosureCheck(PolyGrid *grid, RealGeom *xfn, RealGeom *area);

// Calculate face to nodes connectivity
// NOTE: F2N is a reference to f2n so as to reduce memory use.
IntType **CalF2N(PolyGrid *grid);

// Calculate face number of each cell
IntType *CalnFPC(PolyGrid *grid);

// Calculate cell to faces connectivity
IntType **CalC2F(PolyGrid *grid);

// Calculate number of neighbor cells sharing face
IntType *CalnCPC(PolyGrid *grid); ///< Used for CellColoring

// Calculate cell to cells connectivity
IntType **CalC2C(PolyGrid *grid); ///< Used for CellColoring

// Calculate node number of each cell
IntType *CalnNPC(PolyGrid *grid); ///< Used for gradient calculation / Used for Divide&Conquer

// Calculate cell to nodes connectivity
IntType **CalC2N(PolyGrid *grid); ///< Used for gradient calculation / Used for Divide&Conquer

// //
// void CalCNNCF(PolyGrid *grid);

//!<  Calculate the information of nodes on symmetry plane, including type flag(1) and normal vector.
void FindNodeSYMM(PolyGrid *grid); ///< Used for gradient calculation

// // Sort the sub-array of size npt
// void ReorderPnts(IntType *pnt, IntType npt);

// void BreakFaceLoop(PolyGrid *grid);

// /// \brief  čŽĄçŽçťç˝ć źčçšç¸éťçç˛ç˝ć źć źĺżä¸Şć
// IntType *CalnCCPN(PolyGrid *grid);
// /// \brief  čŽĄçŽçťç˝ć źčçšç¸éťçç˛ç˝ć źć źĺżĺşĺ
// IntType **CalN2CC(PolyGrid *grid);

// // Calculate the volume number of node excluding ghost cell
// IntType *CalnVPN(PolyGrid *grid);
// // Calculate node to volume connectivity excluding ghost cell
// IntType **CalN2V(PolyGrid *grid);

// Calculate node to cells connectivity
IntType **CalN2C(PolyGrid *grid); ///< Used for gradient calculation
// Calculate cell number of each node
IntType *CalnCPN(PolyGrid *grid); ///< Used for gradient calculation

#endif //~MF_GRID_POLYHEDRA_H
