/*!
 * @file        load_grid.cpp
 * @brief       加载网格的主要函数，包含对网格数据的预处理
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
 * <tr><td> 2023-06-13  <td> 1.1      <td> Wisces  <td> 添加了 InitParameter() 函数，用于初始化相关参数
 * <tr><td> 2023-07-07  <td> 1.2      <td> Wisces  <td> 添加了针对梯度计算的相关函数
 * <tr><td> 2023-06-15  <td> 2.0      <td> Wisces  <td> 删除了 InitParameter() 函数，添加头文件 para_field_global.h，采用全局变量形式存储数据
 * <tr><td> 2023-06-19  <td> 3.0      <td> Wisces  <td> 添加了 MPI 并行（使用条件编译）
 * -# 添加了 SpecifyNeighbors()、CommGraph()、CommGraph_node() 函数，用来统计进程的邻居
 * -# 添加了 SpecifyBC() 函数，用来设置边界
 * -# 添加了 SetUpComm()、SetUpComm_Node() 函数，用于建立边界面 / 边界点的传值机制
 * -# 添加了 PassInterfaceData() 和 SetGhostVariables() 函数，用来进行并行传值和设计边界虚网格
 * <tr><td> 2023-06-28  <td> 3.1      <td> Wisces  <td> 精简代码
 * -# 删除了 SetUpComm_Node() 函数
 * -# 删除了 CommGraph_node() 函数
 * -# 删除了 AdditionalInfoForGeometry() 函数
 * <tr><td> 2023-07-03  <td> 4.0      <td> Wisces  <td> 添加了 OpenMP 并行（使用条件编译）
 * -# 添加了 Decoupling() 函数，用来实现四种共享内存并行编程模型的算法
 * -# 添加了 ReorderCellforLUSGS() 函数，用于针对 LUSGS 进行网格排序（内含 CellColoring）
 * <tr><td> 2023-07-12  <td> 5.0      <td> Wisces  <td> 添加了 TDTree 线程级并行（使用条件编译）
 * -# 在 Decoupling() 函数添加了构建任务树等函数
 * </table>
 */

//!< direct head file
#include "load_grid.h"

//!< C/C++ head files
// #include <string>
// #include <cassert>
#include <cmath>     ///< sqrt() sin() cos() pow()
#include <stdio.h>   ///< printf()
#include <iostream>  ///< std::cout
#include <algorithm> ///< std::swap(), std::reverse()
#include <cmath>     ///< std::pow()
#include <sstream>   ///< ostringstream
#include <fstream>   ///< ifstream
#include <unistd.h>  ///< get_current_dir_name()
#include <cstdlib>   ///< malloc() & free()

//!< user defined head files
#include "grid_polyhedra.h"
#include "memory_util.h"
#include "para_field_global.h"
#include "cal_gradient.h"

//!< head file relying on condition-compiling
#ifdef MF_MPICH
#include <mpi.h>
#endif

#ifdef MF_OPENMP
#include <omp.h>
#endif

#ifdef MF_MPICH
extern int myid;
extern int numprocs;
extern int myZone;        ///< = myid+1
extern MPI_Comm GridComm; ///< for each grid
#endif

using namespace std;

/*!
 * @brief       Read grid file through file name, which mainly contains related geometic data
 * @param       grid
 * @remarks     modify according to the fun [void Simulation::LoadGrid(vector<string> &grid_dir_of_zone)]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-05-21
 */
void LoadGrid(PolyGrid *grid, BCond *bc)
{
    //!< determine the directory where grids locate
    // string grid_dir = "/public5/home/sch8427/zhouzheng/Miniflux/MiniFlux_gaspi/grid-4z/grid_chnt-1_xf_v3/";
    string grid_dir = "/public5/home/sch8427/zhouzheng/Miniflux/MiniFlux_gaspi/grid-4z/";
    // string grid_dir = "/public5/home/sch8427/daidai/flowstar_lite_GMG_original3/flowstar_lite_GMG/flowstar_lite_GMG/grid_c_m/16/";
    // string grid_dir = "/public5/home/sch8427/lf/c_grid/";
    //  char *buff = get_current_dir_name();
    //  string grid_dir = "/grid-4z/";
    //  grid_dir = buff + grid_dir;
    //  free(buff);
    //!< extra grid data 边界条件（用于设置名称）
    // GridIO::ExtraGridData extra_grid_data;

    //!< Read grid files
    string grid_file;

#ifdef MF_MPICH
    grid_file = grid_dir + "mmgrid" + int2str(myZone) + ".in";
    if (numprocs == 1)
        grid_file = grid_dir + "serial_grid.mfl";
#else
    grid_file = grid_dir + "serial_grid.mfl";
#endif

    //!< const IntType mg_levels = 1; 多重网格相关
    ReadMFlowGrid_binary(grid, grid_file /*, extra_grid_data*/);

    ///!< 读入包含网格边界相关参数的文件 bc.para
    string bc_file = grid_dir + "bc.para";
    ReadMFlowBoundaryCondition(bc, bc_file.data()); // bc_file.c_str()

    //!< Bind grids to zone and pass patch names into zone
    //!< @todo
}

/*!
 * @brief       Read grid of MFlow binary format generated by MFlow preprocessing. Each data of parallel grid saved in one file, the filename is different for different processor.
 * @param       grids
 * @param       filename
 * @note        本项目中若使用多重网格，则其他粗网格并不是直接通过文件读取，而是通过方法调用进行计算的
 * @note        modify according to the fun [void ReadMFlowGrid_binary(PolyGrid **grids, const IntType mg,
 *                                              const string &filename, ExtraGridData &extra_grid_data)]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date
 */
void ReadMFlowGrid_binary(PolyGrid *grid, const string &filename /*, ExtraGridData &extra_grid_data*/)
{
    //!< make sure grid is a real object
    // assert(grid != NULL);

    FILE *fp = NULL;
    if ((fp = fopen(filename.c_str(), "rb")) == NULL)
    {
        cerr << endl
             << "Could not open file " << filename << endl;
        exit(EXIT_FAILURE);
    }

#ifdef MF_MPICH
    cout << endl
         << "Reading binary grid of FlowStar format - " << filename << " - proccessor " << myid << endl;
#else
    cout << endl
         << "Reading binary grid of FlowStar format - " << filename << endl;
#endif

    // grid->SetLevel(0);

    IntType nTNode, nTFace, nTCell;
    fread(&nTNode, sizeof(IntType), 1, fp);
    fread(&nTFace, sizeof(IntType), 1, fp);
    fread(&nTCell, sizeof(IntType), 1, fp);
    grid->SetNTNode(nTNode);
    grid->SetNTFace(nTFace);
    grid->SetNTCell(nTCell);

#ifdef MF_DEBUG
#ifdef MF_MPICH
    cout << "Grid has " << (long)nTNode << " points, " << (long)nTFace << " faces and " << (long)nTCell << " cells - proccessor " << myid << endl;
#else
    cout << "Grid has " << (long)nTNode << " points, " << (long)nTFace << " faces and " << (long)nTCell << " cells" << endl;
#endif

#endif

    RealGeom *x = NULL, *y = NULL, *z = NULL;
    snew_array_1D(x, nTNode);
    snew_array_1D(y, nTNode);
    snew_array_1D(z, nTNode);

    grid->SetX(x);
    grid->SetY(y);
    grid->SetZ(z);

    fread(x, sizeof(RealGeom), nTNode, fp);
    fread(y, sizeof(RealGeom), nTNode, fp);
    fread(z, sizeof(RealGeom), nTNode, fp);

    IntType *facnod = NULL, *nNPF = NULL;
    snew_array_1D(facnod, nTFace + 1);
    snew_array_1D(nNPF, nTFace);
    fread(nNPF, sizeof(IntType), nTFace, fp);
    grid->SetnNPF(nNPF);

    facnod[0] = 0;
    IntType n = 0;
    for (IntType i = 0; i < nTFace; ++i)
    {
        n += nNPF[i];
        facnod[i + 1] = n;
    }
    IntType *f2n = NULL;
    snew_array_1D(f2n, n);
    grid->Setf2n(f2n);
    fread(f2n, sizeof(IntType), n, fp);

    n = 2 * nTFace;
    IntType *f2c = NULL;
    snew_array_1D(f2c, n);
    grid->Setf2c(f2c);
    fread(f2c, sizeof(IntType), n, fp);

    IntType nBFace = 0;
    IntType count = 0;

    for (IntType i = 0; i < nTFace; ++i)
    {
        IntType c1 = count++;
        IntType c2 = count++;

        if (f2c[c1] < 0) ///< 面的左边是物理边界
        {
            //!< need to reverse the node ordering, 保证点的顺序为右手系
            //!< Note: the first point of face will not move.
            std::reverse(&(f2n[facnod[i] + 1]), &(f2n[facnod[i + 1]]));

            //!< exchange c1 and c2
            std::swap(f2c[c1], f2c[c2]);

            ++nBFace;
        }
        else if (f2c[c2] < 0) ///< 面的右边是物理边界
        {
            ++nBFace;
        }
    }
    sdel_array_1D(facnod);

    grid->SetNBFace(nBFace);

    //!< see if any interfaces
    if (!feof(fp)) ///< exist interfaces
    {
        IntType nIFace, nINode;
        IntType *nbz = NULL, *nbf = NULL, *nbsN = NULL, *nbzN = NULL, *nbrN = NULL;

        fread(&nIFace, sizeof(IntType), 1, fp);
        grid->SetNIFace(nIFace);
        snew_array_1D(nbz, nIFace); ///< 对应的分区号
        snew_array_1D(nbf, nIFace); ///< 对应的面号
        fread(nbz, sizeof(IntType), nIFace, fp);
        fread(nbf, sizeof(IntType), nIFace, fp);

        grid->SetnbZ(nbz);
        grid->SetnbBF(nbf);

        fread(&nINode, sizeof(IntType), 1, fp);
        grid->SetNINode(nINode);
        snew_array_1D(nbsN, nINode); ///< 本块中的并行点序号
        snew_array_1D(nbzN, nINode); ///< 对应点的分区号
        snew_array_1D(nbrN, nINode); ///< 对应点在对应分区中的序号
        fread(nbsN, sizeof(IntType), nINode, fp);
        fread(nbzN, sizeof(IntType), nINode, fp);
        fread(nbrN, sizeof(IntType), nINode, fp);

        grid->SetnbSN(nbsN);
        grid->SetnbZN(nbzN);
        grid->SetnbRN(nbrN);

        //!< ??????
        // IntType nIShadowFace = 0;
        // fread(&nIShadowFace, sizeof(IntType), 1, fp);
        // IntType nShadowFace = 0;
        // fread(&nShadowFace, sizeof(IntType), 1, fp);
    }

    fclose(fp);

    //!< 边界条件相关
}

/*!
 * @brief       读取存放网格对应边界条件的文件
 * @param       filename
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-19
 */
void ReadMFlowBoundaryCondition(BCond *bc, const char *filename)
{
    ifstream in(filename);
    if (!in.is_open())
    {
        cerr << endl
             << "Failed to open " << filename << endl;
        exit(EXIT_FAILURE);
    }

#ifdef MF_MPICH
    if (myid == 0)
#endif
        cout << "Reading boundary condition of Grid - " << filename << endl;

    int n_patch;
    int patch_id, type;
    in >> n_patch;
    // cout << "n_patch = " << n_patch << endl;

    for (int i = 0; i < n_patch; i++)
    {
        in >> patch_id >> type;
        // cout << patch_id << " - " << type << endl;
        BCRecord *bcr = new BCRecord(patch_id, type);
        bc->AddBCRecord(bcr);
    }

    in.clear();
    in.close();

    // cout << " Done loading " << filename << endl;
}

/*!
 * @brief       工具函数：int 转化为 string 类型函数
 * @details     利用字符串流 ostringstream 将整型转化为 string 类型
 * @param       source
 * @return      string
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-18
 */
string int2str(const IntType source)
{
    ostringstream ss("");
    ss << source;
    return ss.str();
}

/*!
 * @brief       设计实现了基于循环级与任务级两种共享内存并行编程模型的总共四种并行算法
 * @param       grid
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-03
 */
void Decoupling(PolyGrid *grid)
{
    //!< 面着色技术（循环并行类）
#ifdef OMP_FaceColor
    IntType igroupsize, bgroupsize;
    igroupsize = 208000;
    bgroupsize = 5400;
    FaceColoring(grid, bgroupsize, igroupsize);
#endif

    //!< 分组着色技术（循环并行类）
#if (defined MF_OPENMP) && (defined OMP_GroupColor)
    grid->GroupColorSuccess = false;
    grid->GroupColorSuccess = GroupColoring(grid, true);
    if (!grid->GroupColorSuccess)
    {
#ifdef MF_MPICH
        cout << "Rank - " << myid << " skip group coloring." << endl;
#else
        cout << "Skip group coloring." << endl;
#endif
    }
#endif

    //!< 剖分复制策略（任务并行类）
#if (defined MF_OPENMP) && (defined OMP_DIVREP)
    grid->threads = omp_get_max_threads();
    grid->DivRepSuccess = DivRep(grid);
    if (!grid->DivRepSuccess)
    {
        cerr << "Error in func[Decoupling()] of file[load_grid.cpp]!" << endl;
        exit(1);
    }
// #if (defined FS_SIMD) && (defined BoundedColoring)
//     grid->endIndex_bFace_vec = NULL;
//     grid->endIndex_iFace_vec = NULL;
//     snew_array_1D(grid->endIndex_bFace_vec, grid->threads);
//     snew_array_1D(grid->endIndex_iFace_vec, grid->threads);
//     colorAfterDivRep(grid);
// #endif
#endif

    //!< --------------- deleted ---------------
    //     //!< 分治法（任务并行类）
    // #if (defined MF_OPENMP) && (defined OMP_DIVCON)
    //     DC_create_tree(grid);
    // #endif
    //!< --------------- deleted ---------------

    //!< TDTree
#ifdef TDTREE
    CreateTree(grid);
    TDTree_CellReordering(grid);
    TDTree_NodeReordering(grid);
    TDTree_FaceReordering(grid);
#endif
}

/*!
 * @brief       Do some preprocessing on grid data for subsequent solution calculations
 * @param       grid
 * @remarks     modify according to the fun [void Zone::InitZone()] and [void NSSolver::Init()], which is called in [simu->Start() - Init() - InitZone()/InitSolver()]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-05
 */
void PreprocessGrid(PolyGrid *grid, BCond *bc)
{
    //!< TEST
#ifdef MF_DEBUG
#ifdef MF_MPICH
    if (myid == 0)
#endif
    {
        cout << endl
             << "The incoming flowing condition: " << endl;
        cout << "================================================================================" << endl;
        cout << "      T              = " << T << endl;
        cout << "      re             = " << re << endl;
        cout << "      uinf(u_s)      = " << u_s << endl;
        cout << "      vinf(v_s)      = " << v_s << endl;
        cout << "      winf(w_s)      = " << w_s << endl;
        cout << "      p(p_s)         = " << p_s << endl;
        cout << "      rhoinf(rho_s)  = " << rho_s << endl;
        cout << "      amu            = " << amu << endl;
        cout << "      ainf           = " << ainf << endl;
        cout << "      p_bar          = " << p_bar << endl;
        // cout << "      amu120   = " << amu120 << endl;
        cout << "      p_stag         = " << p_stag << endl;
        // cout << "      rho_stag = " << rho_stag << endl;
        cout << "      e_stag         = " << e_stag << endl;
        // cout << "      t_stag   = " << t_stag << endl;
        cout << "      p_min          = " << p_min << endl;
        cout << "      p_max          = " << p_max << endl;
        cout << "      rho_min        = " << rho_min << endl;
        cout << "      rho_max        = " << rho_max << endl;
        cout << "      p_break        = " << p_break << endl;
        cout << "      e_stag_max     = " << e_stag_max << endl;
        cout << "      the machine zero = " << -iexp << endl;
        cout << "================================================================================" << endl;
    }
#endif

    //!< put neighbor information into each grid (非多重网格时，只有一个 grid)
    SpecifyNeighbors(grid);

    //!< Assign boundary condition for each boundary face of grid.（对于多重网格另处理）
    SpecifyBC(grid, bc);

#ifdef MF_DEBUG
#ifdef MF_MPICH
    if (myid == 0)
#endif
        cout << "Exit SpecifyNeighbors/SpecifyBC!" << endl;
#endif

        //!< 建立边界面 / 边界点的传值机制
#ifdef MF_MPICH
    //!< 如果考虑多重网格，则 for (j = 0; j < nGrids; j++)
    grid->SetUpComm();
    grid->SetUpComm_Node(); ///< 或许可以省略
#ifdef MF_DEBUG
#ifdef MF_MPICH
    if (myid == 0)
#endif
        cout << "Exit SetUpComm/SetUpComm_Node!" << endl;
#endif
#endif

    //!< 如果考虑多重网格，则 for (j = 0; j < nGrids; j++)
    //!< 计算网格的几何量：面心、面积、面单位法向矢量、体心、体积
    grid->ComputeMetrics();
    // grid->AdditionalInfoForGeometry();
#ifdef MF_DEBUG
#ifdef MF_MPICH
    if (myid == 0)
#endif
        cout << "Exit ComputeMetrics!" << endl;
#endif

    ///< Wisces: 暂时不考虑湍流计算
    /**
    // Compute the distance to wall for turbulence model
    if (i >= 0 && vis_mode == S_A_MODEL)
    { // i==0
        grid->ComputeCellDist();
    }
    */

    // grid->FindCellLayerNo(); ///< Find Cell Layer No. from wall（为计算梯度做准备）
    //!< 重排序 内含 CellColoring
    ReorderCellforLUSGS(grid);
#ifdef MF_DEBUG
#ifdef MF_MPICH
    if (myid == 0)
#endif
        cout << "Exit ReorderCellforLUSGS!" << endl;
#endif
    //!< allocate memory for flow field
    AllocAndInitFlowfield(grid);

    // UpdateInterfaceData();
    // SetGhostVariables(pgrid);
    // ComputeVis_l(pgrid);
    // InitVis_t(pgrid);

    //!< 初始化时，并行传值、设置边界值、计算梯度、预处理计算，计算粘性，包括相关内存的分配
    //!< 并行传值
    PassInterfaceData(grid);
    //!< 计算边界虚网格
    SetGhostVariables(grid);

#ifdef MF_DEBUG
#ifdef MF_MPICH
    if (myid == 0)
#endif
        cout << "Exit AllocAndInitFlowfield/PassInterfaceData/SetGhostVariables!" << endl;
#endif

    //!< 计算距离分之一权，为计算梯度做准备
    ComputeWeight3D_Node(grid);

    //!< 寻找对称边界，为计算梯度做准备
    FindNodeSYMM(grid);
#ifdef MF_DEBUG
#ifdef MF_MPICH
    if (myid == 0)
#endif
        cout << "Exit ComputeWeight3D_Node/FindNodeSYMM!" << endl;
#endif

    //!< 为梯度分配内存并计算
    AllocAndCalQuantityGradient(grid);
#ifdef MF_DEBUG
#ifdef MF_MPICH
    if (myid == 0)
#endif
        cout << "Exit AllocAndCalQuantityGradient!" << endl;
#endif

    //!< 计算层流动力粘性
}

/*!
 * @brief       针对 LUSGS，进行网格排序（内含 CellColoring）
 * @param       grid
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-06
 */
void ReorderCellforLUSGS(PolyGrid *grid)
{
    register IntType i;
#if (defined MF_OPENMP) && (defined OMP_CellColor)
    // ((PolyGrid *)grids[i])->ReorderCellforLUSGS_1();
    // ((PolyGrid *)grids[i])->ReorderCellforLUSGS_2();
    // ((PolyGrid *)grids[i])->ReorderCellforLUSGS_3();
    //!< find cell color for lusgs, 0 or 1 to be chosen in param where 0 is imbalance greedy algorithm, 1 is balance algorithm
    grid->LUSGSGridColor(1); ///< cell coloring in LUSGS
#else
    grid->ReorderCellforLUSGS_0(); ///< no reorder
#endif
}

/*!
 * @brief       Allocate memory and Initialize for flow field
 * @param       grid
 * @remarks     modify according to the fun [void NSSolver::AllocateFlowfieldMemory(PolyGrid *grid)] and [void NSSolver::InitGridVar(PolyGrid *grid)]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-05
 */
void AllocAndInitFlowfield(PolyGrid *grid)
{
    RealFlow *rho = NULL, *u = NULL, *v = NULL, *w = NULL, *p = NULL;
    // RealFlow rhoP, uP, vP, wP, pP; // for the parameters
    IntType nTCell = grid->GetNTCell();
    IntType n = nTCell + grid->GetNBFace();

    snew_array_1D(rho, n);
    snew_array_1D(u, n);
    snew_array_1D(v, n);
    snew_array_1D(w, n);
    snew_array_1D(p, n);

    //!< 为全流场赋来流值
    for (int i = 0; i < n; i++)
    {
        rho[i] = rho_s;
        u[i] = u_s;
        v[i] = v_s;
        w[i] = w_s;
        p[i] = p_s;
    }

    grid->SetRho(rho);
    grid->SetU(u);
    grid->SetV(v);
    grid->SetW(w);
    grid->SetP(p);

    ///< Wisces: 不考虑非定常的情况
}

/*!
 * @brief       Assign neighbors for each grid
 * @param       grid
 * @note        modify according to the fun [void Zone::SpecifyNeighbors()]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-19
 */
void SpecifyNeighbors(PolyGrid *grid)
{
    IntType nNeighbor = 0, *nb = NULL;   ///< num. of neighbors linked by face, the neighbor rank no. for face
    IntType nNeighborN = 0, *nbN = NULL; ///< num. of neighbors linked by node, the neighbor rank no. for node

    IntType i;

    CommGraph(grid, nNeighbor, nb);
    CommGraph_node(grid, nNeighborN, nbN);

    //!< put neighbor information into each grid
    //!< 不保留对多重网格的处理，否则 for(i = 0; i < nGrids; i++)
    grid->SetNumberOfFaceNeighbors(nNeighbor);
    grid->SetNumberOfNodeNeighbors(nNeighborN);
    if (nNeighbor > 0)
        grid->SetFaceNeighborZones(nb);
    if (nNeighborN > 0)
        grid->SetNodeNeighborZones(nbN);

    //<! ifndef MF_MPICH，nNeighbor == 0
}

/*!
 * @brief       Generate the communication graph (faces)
 * @param       grid
 * @param       nNeighbor
 * @param       nb
 * @note        modify according to the fun [void Simulation::CommGraph()]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-19
 */
void CommGraph(PolyGrid *grid, IntType &nNeighbor, IntType *&nb)
{
    IntType *nbZ, nIFace;
    IntType maxz = 0;
    IntType *count;

    IntType i, j;

    nbZ = grid->GetnbZ(); ///< （对应的分区号）parallel information for interfaces: zone(process) No. in that zone(process). Size is nIFace.
    nIFace = grid->GetNIFace();
    for (i = 0; i < nIFace; i++)
        maxz = MAX(maxz, nbZ[i]);
    maxz++;
    //!< 省略判断
    // #ifndef MF_MPICH
    //     assert(maxz == 1);
    // #endif
    count = NULL;
    snew_array_1D(count, maxz);

    for (j = 0; j < maxz; j++)
        count[j] = 0;

    for (i = 0; i < nIFace; i++)
        count[nbZ[i]] = 1;

    nNeighbor = 0;
    for (j = 0; j < maxz; j++)
    {
        if (count[j])
            nNeighbor++;
    }

    if (nNeighbor > 0)
    {
        nb = NULL;
        snew_array_1D(nb, maxz);
        nNeighbor = 0;
        for (j = 0; j < maxz; j++)
        {
            if (count[j])
                nb[nNeighbor++] = j;
        }

        //!< reorder the nb so that nb[j] in the increasing order
        IntType nb_tmp, k;
        for (j = 0; j < nNeighbor; j++)
        {
            for (k = j + 1; k < nNeighbor; k++)
            {
                if (nb[j] > nb[k])
                {
                    nb_tmp = nb[j];
                    nb[j] = nb[k];
                    nb[k] = nb_tmp;
                }
            }
        }
    }
    sdel_array_1D(count);
#ifdef MF_DEBUG
#ifdef MF_MPICH
    cout << endl
         << "Parallel face linked neighbor number for process - " << myid << " is " << nNeighbor << endl;
#else
    cout << endl
         << "Face linked neighbor number is " << nNeighbor << endl;
#endif
#endif
}

/*!
 * @brief       Generate the communication graph (nodes)
 * @param       grid
 * @param       nNeighborN
 * @param       nbN
 * @note        modify according to the fun [void Simulation::CommGraph_node()]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-19
 */

void CommGraph_node(PolyGrid *grid, IntType &nNeighborN, IntType *&nbN)
{
    IntType *nbZN, nINode;
    IntType maxz = 0;
    IntType *count;

    IntType j;

    nbZN = grid->GetnbZN(); ///< （对应点的分区号）parallel information for inter-nodes: zone No. on other zone. Size is nINode.
    nINode = grid->GetNINode();
    for (j = 0; j < nINode; j++)
        maxz = MAX(maxz, nbZN[j]);
    maxz++;
    //!< 省略判断
    // #ifndef MF_MPICH
    //     assert(maxz == 1);
    // #endif
    count = NULL;
    snew_array_1D(count, maxz);

    for (j = 0; j < maxz; j++)
        count[j] = 0;

    for (j = 0; j < nINode; j++)
        count[nbZN[j]] = 1;

    nNeighborN = 0;
    for (j = 0; j < maxz; j++)
    {
        if (count[j])
            nNeighborN++;
    }

    if (nNeighborN > 0)
    {
        nbN = NULL;
        snew_array_1D(nbN, maxz);
        nNeighborN = 0;
        for (j = 0; j < maxz; j++)
        {
            if (count[j])
                nbN[nNeighborN++] = j;
        }

        //!< reorder the nbN so that nbN[j] in the increasing order
        IntType nbN_tmp, k;
        for (j = 0; j < nNeighborN; j++)
        {
            for (k = j + 1; k < nNeighborN; k++)
            {
                if (nbN[j] > nbN[k])
                {
                    nbN_tmp = nbN[j];
                    nbN[j] = nbN[k];
                    nbN[k] = nbN_tmp;
                }
            }
        }
    }
    sdel_array_1D(count);

#ifdef MF_DEBUG
#ifdef MF_MPICH
    cout << "Parallel node linked neighbor number for process - " << myid << " is " << nNeighborN << endl;
#else
    cout << "Node linked neighbor number is " << nNeighborN << endl;
#endif
#endif
}

/*!
 * @brief       Assign boundary conditions
 * @param       grid
 * @note        Before this, c2 is negative, indicating the patch number; After this, c2 = b_face_no + n_cells
 * @note        modify according to the fun [void Zone::SpecifyBC()]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-21
 */
void SpecifyBC(PolyGrid *grid, BCond *bc)
{
    IntType j, k, c2;
    IntType n_patch, nNeighbor;
    IntType *f2c, nBFace, nTCell;

    //!< physical boundary conditions only. Grid related bc is stored in Grid object
    // BCond bc;
    //!< for m6wing_12w
    // bc->AddBCRecord(bcr1);
    // bc->AddBCRecord(bcr2);
    // bc->AddBCRecord(bcr3);

    //!< for m6wing_25w
    // bc->AddBCRecord(bcr1);
    // bc->AddBCRecord(bcr2);
    // bc->AddBCRecord(bcr3);
    // bc->AddBCRecord(bcr4);
    // bc->AddBCRecord(bcr5);

    //!< for CHNT1_600w
    // bc->AddBCRecord(bcr1);
    // bc->AddBCRecord(bcr2);
    // bc->AddBCRecord(bcr3);
    // bc->AddBCRecord(bcr4);
    // bc->AddBCRecord(bcr5);
    // bc->AddBCRecord(bcr6);
    // bc->AddBCRecord(bcr7);
    // bc->AddBCRecord(bcr8);

    n_patch = bc->GetNoBCR();
    nNeighbor = grid->GetNumberOfFaceNeighbors();
    if (nNeighbor > 0)
    {
        IntType bcr_int;
        //!< there must be some interface bcs, find the BCRecord
        BCRecord *inter = NULL, *bcr = NULL;
        for (j = 0; j < n_patch; j++)
        {
            //!< find the interface bcr
            bcr = bc->GetBCRecord(j);
            if (bcr->GetType() == INTERFACE)
            {
                inter = bcr;
                bcr_int = inter->GetPatchID();
                break;
            }
        }
        if (inter == 0)
        {
            //!< if not have ,then add it
            bcr_int = ++n_patch;
            inter = new BCRecord(bcr_int, INTERFACE);
            // inter = new BCRecord();
            // inter->SetType(INTERFACE);
            // inter->SetTypeSymbol("Interface");
            // inter->SetPatchID(bcr_int);
            bc->AddBCRecord(inter);
        }
        IntType nIFace;
        //!< 如果考虑多重网格，则 for (j = 0; j < nGrids; j++)
        f2c = grid->Getf2c();
        nBFace = grid->GetNBFace();
        nIFace = grid->GetNIFace();
        for (k = nBFace - nIFace; k < nBFace; k++)
        {
            f2c[k + k + 1] = -bcr_int;
        }
    }

    //!< 如果考虑多重网格，则 for (j = 0; j < nGrids; j++)
    nBFace = grid->GetNBFace();
    nTCell = grid->GetNTCell();
    f2c = grid->Getf2c();

    BCRecord **bcrs = NULL;
    snew_array_1D(bcrs, nBFace);
    grid->Setbcr(bcrs); ///< associate each boundary. face with a record

    for (k = 0; k < nBFace; k++)
    {
        c2 = -f2c[k * 2 + 1];
        // if (c2 <= 0)
        // {
        //     assert(0);
        // }
        // if (c2 > n_patch)
        // {
        //     assert(0);
        //     bcrs[k] = bc->GetBCRecord(n_patch - 1);
        // }
        // else
        // {
        //     bcrs[k] = bc->GetBCRecord(c2 - 1);
        // }
        bcrs[k] = bc->GetBCRecord(c2 - 1);
        f2c[k * 2 + 1] = k + nTCell;
    }
}

/*!
 * @brief       初始化或得到新的流场变量后，进行并行传值
 * @param       grid
 * @note        modify according to the fun [void NSSolver::PassInterfaceData(PolyGrid *grid)]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-27
 */
void PassInterfaceData(PolyGrid *grid)
{
    // RealFlow *rho = grid->GetRho();
    // RealFlow *u = grid->GetU();
    // RealFlow *v = grid->GetV();
    // RealFlow *w = grid->GetW();
    // RealFlow *p = grid->GetP();

    // IntType kNVar = 5;
    RealFlow **q = NULL;
    snew_array_1D(q, 5);
    // q[0] = static_cast<RealFlow *>(rho);
    // q[1] = static_cast<RealFlow *>(u);
    // q[2] = static_cast<RealFlow *>(v);
    // q[3] = static_cast<RealFlow *>(w);
    // q[4] = static_cast<RealFlow *>(p);
    q[0] = static_cast<RealFlow *>(grid->GetRho());
    q[1] = static_cast<RealFlow *>(grid->GetU());
    q[2] = static_cast<RealFlow *>(grid->GetV());
    q[3] = static_cast<RealFlow *>(grid->GetW());
    q[4] = static_cast<RealFlow *>(grid->GetP());
#ifdef GASPI
    // for (IntType i = 0; i < 5; ++i)
    // {
    //     grid->CommInterfaceDataMPI(q[i]); ///< lite
    // }
    grid->GASPIRecvSendVarNeighbor_Togeth(5, q);
#elif MF_MPICH
    grid->RecvSendVarNeighbor_Togeth(5, q);
#endif
    sdel_array_1D(q);
}

/*!
 * @brief       Set the Ghost Variables object, namely set flow variables at ghost cells
 * @param       grid
 * @note        modify according to the fun [void SetGhostVariables(PolyGrid *grid) in solver_ns.cpp]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-27
 */
void SetGhostVariables(PolyGrid *grid)
{
    IntType i, c1, c2, count, /*steady, vis_mode,*/ type, wmark;
    IntType nBFace = grid->GetNBFace();
    IntType nTCell = grid->GetNTCell();
    IntType *f2c = grid->Getf2c();
    IntType n = nTCell + nBFace;

    // RealFlow *rho = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "rho");
    // RealFlow *u = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "u");
    // RealFlow *v = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "v");
    // RealFlow *w = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "w");
    // RealFlow *p = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "p");
    RealFlow *rho = grid->GetRho();
    RealFlow *u = grid->GetU();
    RealFlow *v = grid->GetV();
    RealFlow *w = grid->GetW();
    RealFlow *p = grid->GetP();

    RealGeom *xfn = grid->GetXfn();
    RealGeom *yfn = grid->GetYfn();
    RealGeom *zfn = grid->GetZfn();
    RealGeom *vgn = grid->GetFaceNormalVelocity();
    RealGeom *BFacevgx = grid->GetBoundaryFaceVelocityX();
    RealGeom *BFacevgy = grid->GetBoundaryFaceVelocityY();
    RealGeom *BFacevgz = grid->GetBoundaryFaceVelocityZ();
    BCRecord **bcr = grid->Getbcr();

    RealFlow vn;
    RealFlow rho00, u00, v00, w00, p00, vtx, vty, vtz, riemp, riemm;
    RealFlow vnp, vnm, cp, cm, /*gam, p_bar,*/ p_t, yc, entr, vnb, cb;
    // RealFlow gascon;
    RealFlow rhow;

    RealFlow rhom, um, vm, wm, pm, rhop, up, vp, wp, pp;
    RealFlow tw;

    rho00 = rho_s;
    u00 = u_s;
    v00 = v_s;
    w00 = w_s;
    p00 = p_s;

    // grid->GetData(&rho00, REAL_FLOW, 1, "rho");
    // grid->GetData(&u00, REAL_FLOW, 1, "u");
    // grid->GetData(&v00, REAL_FLOW, 1, "v");
    // grid->GetData(&w00, REAL_FLOW, 1, "w");
    // grid->GetData(&p00, REAL_FLOW, 1, "p");
    // grid->GetData(&p_bar, REAL_FLOW, 1, "p_bar");

    RealFlow norm_of_uvw = sqrt(u00 * u00 + v00 * v00 + w00 * w00);
    RealFlow eps_of_farfield_vn = 0.0;
    // grid->GetData(&eps_of_farfield_vn, REAL_FLOW, 1, "eps_of_farfield_vn", 0); ///< 数据库中不存在，同时也没有写入更新

    // RealFlow rho_min, rho_max, p_min, p_max;
    // grid->GetData(&rho_min, REAL_FLOW, 1, "rho_min");
    // grid->GetData(&rho_max, REAL_FLOW, 1, "rho_max");
    // grid->GetData(&p_min, REAL_FLOW, 1, "p_min");
    // grid->GetData(&p_max, REAL_FLOW, 1, "p_max");

    // grid->GetData(&steady, INT, 1, "steady");
    // grid->GetData(&vis_mode, INT, 1, "vis_mode");
    // grid->GetData(&gam, REAL_FLOW, 1, "gam");
    // grid->GetData(&gascon, REAL_FLOW, 1, "gascon");

    count = 0;
    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        c1 = f2c[count++];
        c2 = f2c[count++];
        wmark = 0;

        // Do nothing for interfaces.
        if (type == INTERFACE)
            continue;

        // Assign the variable values for each ghost cell whose index is c2.
        switch (type)
        {
        case WALL:
            p[c2] = p[c1];
            rho[c2] = rho[c1];

            if (vis_mode == INVISCID)
            {
                vn = 2. * (xfn[i] * u[c1] + yfn[i] * v[c1] + zfn[i] * w[c1]);
                // if (!steady)
                // {
                //     vn -= 2 * vgn[i];
                // }
                u[c2] = u[c1] - vn * xfn[i];
                v[c2] = v[c1] - vn * yfn[i];
                w[c2] = w[c1] - vn * zfn[i];
            }
            /**
            else
            {
                if (steady)
                {
                    u[c2] = -u[c1];
                    v[c2] = -v[c1];
                    w[c2] = -w[c1];
                }
                else
                {
                    u[c2] = -u[c1] + 2. * BFacevgx[i];
                    v[c2] = -v[c1] + 2. * BFacevgy[i];
                    w[c2] = -w[c1] + 2. * BFacevgz[i];
                }

                // viscous adiabatic wall
                // nothing!!!

                // viscous iso-thermal wall
                tw = -1.0;
                // bcr[i]->GetBCVar(&tw, REAL_FLOW, "tw", 0); ///< 数据库中不存在，同时也没有写入更新
                if (tw > 0.0)
                {
                    rhow = (p[c2] + p_bar) / gascon / tw;
                    rho[c2] = 2.0 * rhow - rho[c1];
                    if (rho[c2] < 0.0)
                    {
                        rho[c2] = rhow;
                    }
                }
            }
            */
            break;

        case SYMM:
            rho[c2] = rho[c1];
            p[c2] = p[c1];
            vn = 2. * (xfn[i] * u[c1] + yfn[i] * v[c1] + zfn[i] * w[c1]);
            // if (!steady)
            // {                     // zhyb:对称面 vgn 为0，此处本来可以不考虑。但是在粘性计算时，有时可能会采用对称边界条件表示无粘的物面，
            //     vn -= 2 * vgn[i]; // 因此在此需要加上非定常的情况
            // }
            u[c2] = u[c1] - vn * xfn[i];
            v[c2] = v[c1] - vn * yfn[i];
            w[c2] = w[c1] - vn * zfn[i];
            break;

        case FAR_FIELD:
            um = u00;
            vm = v00;
            wm = w00;
            up = u[c1];
            vp = v[c1];
            wp = w[c1];
            // if (!steady)
            // {
            //     um -= BFacevgx[i];
            //     vm -= BFacevgy[i];
            //     wm -= BFacevgz[i];
            //     up -= BFacevgx[i];
            //     vp -= BFacevgy[i];
            //     wp -= BFacevgz[i];
            // }
            rhom = rho00;
            pm = p00 + p_bar;
            rhop = rho[c1];
            pp = p[c1] + p_bar;

            vnm = xfn[i] * um + yfn[i] * vm + zfn[i] * wm;
            vnp = xfn[i] * up + yfn[i] * vp + zfn[i] * wp;
            cm = sqrt(gam * pm / rhom);
            cp = sqrt(gam * pp / rhop);
            riemm = vnm - 2. * cm / (gam - 1.);
            riemp = vnp + 2. * cp / (gam - 1.);

            vnb = 0.5 * (riemp + riemm);
            cb = 0.25 * (riemp - riemm) * (gam - 1.);

            if (fabs(vnb / cb) > 1.)
            { // supersonic
                if (vnb <= 0.0)
                { // inlet
                    rho[c2] = rhom;
                    u[c2] = um;
                    v[c2] = vm;
                    w[c2] = wm;
                    p[c2] = pm;
                }
                else
                { // exit
                    rho[c2] = rhop;
                    u[c2] = up;
                    v[c2] = vp;
                    w[c2] = wp;
                    p[c2] = pp;
                }
            }
            else
            { // subsonic
                RealFlow rela_vnb = vnb / norm_of_uvw;
                if (rela_vnb <= -eps_of_farfield_vn)
                { // inlet
                    vtx = um - vnm * xfn[i];
                    vty = vm - vnm * yfn[i];
                    vtz = wm - vnm * zfn[i];
                    entr = pm / pow(rhom, gam);

                    rho[c2] = pow((cb * cb / (entr * gam)), RealFlow(1. / (gam - 1.)));
                    u[c2] = vtx + vnb * xfn[i];
                    v[c2] = vty + vnb * yfn[i];
                    w[c2] = vtz + vnb * zfn[i];
                    p[c2] = cb * cb * rho[c2] / gam;
                }
                else if (rela_vnb > eps_of_farfield_vn)
                { // exit
                    vtx = up - vnp * xfn[i];
                    vty = vp - vnp * yfn[i];
                    vtz = wp - vnp * zfn[i];
                    entr = pp / pow(rhop, gam);

                    rho[c2] = pow((cb * cb / (entr * gam)), RealFlow(1. / (gam - 1.)));
                    u[c2] = vtx + vnb * xfn[i];
                    v[c2] = vty + vnb * yfn[i];
                    w[c2] = vtz + vnb * zfn[i];
                    p[c2] = cb * cb * rho[c2] / gam;
                }
                else
                {

                    rho[c2] = 0.5 * (rhop + rhom);
                    u[c2] = 0.5 * (up + um);
                    v[c2] = 0.5 * (vp + vm);
                    w[c2] = 0.5 * (wp + wm);
                    p[c2] = 0.5 * (pp + pm);
                }
            }

            rho[c2] = 2 * rho[c2] - rhop;
            u[c2] = 2 * u[c2] - up;
            v[c2] = 2 * v[c2] - vp;
            w[c2] = 2 * w[c2] - wp;
            p[c2] = 2 * p[c2] - pp;
            p[c2] -= p_bar;

            // if (!steady)
            // {
            //     u[c2] += BFacevgx[i];
            //     v[c2] += BFacevgy[i];
            //     w[c2] += BFacevgz[i];
            // }
            break;

        default:
            cout << "Error in SetGhostVariables 001" << endl;
            break;
        }

        // ZHYB:对 c2 单元的 rho 和 p 做限制，不能为负，不能大于 10 倍的驻点值
        rho[c2] = std::max(rho[c2], rho_min);
        rho[c2] = std::min(rho[c2], rho_max);
        p[c2] = std::max(p[c2], p_min);
        p[c2] = std::min(p[c2], p_max);
    }
}
