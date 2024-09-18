//****************************************************************************************************\
//*                                   National Numerical Windtunnel                                  *
//*          MiniFS -- Mini-application for Fluxes-calculation and LUSGS-solver of FlowStar          *
//*   Institute for Quantum Information & State Key Laboratory of High Performance Computing(HPCL)   *
//*                          Parallel Computing Application Technology Lab.                          *
//*                    National University of Defense Technology, Changsha, CHINA                    *
//*                                      C.All rights reserved.                                      *
//****************************************************************************************************/

/**
 * @mainpage MiniFS 代理应用程序
 * <table>
 * <tr><th> Project     <td> MiniFS
 * <tr><th> Author      <td> Wisces
 * <tr><th> Source      <td> /public5/home/sch8427/wwwangqs/PROJECT/MiniFS/
 * </table>
 *
 * @section 项目描述
 * 针对基于非结构网格技术的国产大型通用 CFD 软件 FlowStar，开发其 N-S 求解器中关于通量计算和 LUSGS 求解器部分的代
 * 理应用，以此开展非结构 CFD 软件代理应用的开发与验证、代理应用的并行实现、代理应用与目标应用的性能建模等研究工作。
 *
 * @section 项目目标
 * 形成一套通用的针对非结构 CFD 软件代理应用的开发验证方法和性能建模方法，用于指导其他大型国产 CFD 软件甚至其他领
 * 域的大型科学应用软件的代理应用研究。
 *
 * @section 更新日志
 * <table>
 * <tr><th> Date        <th> Version  <th> Author  <th> Description  </tr>
 * <tr><td> 2023-05-18  <td> 1.0      <td> Wisces  <td> 创建初始版本
 * -# 简单从源代码中摘取最基础的代码块，未考虑任何优化方法，暂时不支持 MPI 并行
 * -# 主要的代码块包括加载网格、梯度计算、无粘通量计算、LUSGS 求解计算等
 * </tr>
 * <tr><td> 2023-06-15  <td> 2.0      <td> Wisces  <td> 修改了数据存储模式
 * -# 删除了文件 data_pool.h 和 data_pool.cpp，即摒弃由 DataSafe 和 DataStore 两个类支撑的“项目数据库”
 * -# 增加了文件 para_field_global.h 和 para_field_global.cpp，将相关参数和流场数据设置为全局变量
 * </tr>
 * <tr><td> 2023-06-16  <td> 3.0      <td> Wisces  <td> 添加了 MPI 并行（使用条件编译）
 * -# 增加了文件 parallel_bash.h / parallel_bash.cpp / parallel_mpi.cpp，为项目代码提供了可选择性的 MPI 并行模式
 * -# 读取网格文件时读入事先划分好的并行网格文件
 * -# 事先建立边界面的传值机制
 * -# 在具体的迭代过程中实现全局规约、并行传值等必要操作
 * -# 删除了部分非必要的过程，进一步精简了代码
 * </tr>
 * <tr><td> 2023-07-01  <td> 4.0      <td> Wisces  <td> 添加了 OpenMP 并行（使用条件编译）
 * -# 添加了文件 parallel_omp.h / parallel_omp.cpp，为项目代码提供了可选择性的 OpenMP 并行模式
 * -# 分别实现了基于循环级共享内存并行编程模型的规约策略、面着色技术、分组着色技术，以及基于任务级共享内存并行编程模型
 *    的剖分复制策略、分治法总共四种并行算法
 * -# 规约策略：的基本原理是将面循环计算的 gather 和 scatter 过程分离，即面循环计算得到物理量后不更新到体上而是先存起
 *    来，然后再进行体循环，完成物理量的规约（缺点是耗费额外的内存，数据局部性较差，在读取数据和写入数据时不能很好地重
 *    用数据）
 * -# 面着色技术 FaceColoring：对相邻面使用不同颜色（color）进行标记，属于同一颜色的面之间不会产生任何数据冲突，从而
 *    保证同一颜色组内的计算可以并行开展（缺点是在不连续的面上循环计算，数据局部性较差，甚至低于规约策略）
 * -# 分组着色技术 GroupColoring：与面着色处理单个面单元不同的是，分组着色处理的基本单元是拥有固定面数的分组。着色后
 *    同一种颜色的面中 含有若干个 groupsize 的面分组，组与组之间不存在数据竞争，需要设置调度模式 schedule(static,
 *    groupsize) 表明每个线程会处理 groupsize 大小的面分组单元。当 groupsiz e取 1 时分组着色就退化成了普通的着色算
 *    法（改善面着色算法的数据局部性，缺点是 groupsize 的大小选取经验性强）
 * -# 剖分复制策略（Divide & Republication）：在为 MPI 并行计算完成区域分解后，在进程内部根据线程数再进行一次分区（可
 *    继续使用图分区工具 METIS）。该策略的核心是通过工作复制来保证线程安全，没有直接在循环上进行并行，如果两个线程共享
 *    一个面，那这两个进程将复制这个面的工作，在写回操作时，通过额外的逻辑判断保证每个线程只负责自己分区内的体单元数据
 *    更新，从而避免数据冲突（尽量让单个线程处理一组内存位置相近的单元来提升线程内部的数据局部性，优点是每个线程都使用
 *    其专用缓冲区，避免了工作线程之间进行全局同步的开销，而且不用改变原始的变量存储数据结构）
 * -# 分治法（Divide & Conquer）：对 MPI 进程内的网格进一步剖分，递归地构建任务树，并基于任务树遍历实现 task 并行算
 *    法。在每个递归级别上，都会创建三个子分区，包括两个独立的左右子分区和一个分隔分区，其中分隔分区由分区交界面上的网
 *    格实体所构成，左右子分区不具备连通性，从而无数据依赖可并行执行。所有子分区（包括分隔分区）将继续递归分下去，终止
 *    条件是叶子节点的网格单元数小于某个阈值，阈值大小一般选择 Cache 大小，以发挥最大性能。执行计算时，采用 task 并行
 *    编程模型对树节点进行递归遍历执行：左分区和右分区并行执行，一旦它们同步，就执行分隔区，从而实现基于任务级的共享内
 *    存并行
 * -# 单元着色技术 CellColoring
 * </tr>
 * <tr><td> 2023-07-12  <td> 5.0      <td> Wisces  <td> 添加了 TDTree 线程级并行（使用条件编译）
 * -# 删除了效率偏低的分治法，改用 TDTree 实现基于任务级的共享内存并行
 * -# 添加了文件 TDTree.h，改写无粘通量计算、梯度计算、时间步计算、LUSGS 求解计算等主要计算 Kernel
 * -#
 * </tr>
 * </table>
 * ******************************************************************************************************
 */

/*!
 * @file        main.cpp
 * @brief       项目主函数文件
 * @details     主要包含 main 函数入口
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @version     1.0
 * @date        2023-05-18
 * @copyright   Copyright (c) 2023  PERSONAL
 *
 * @par 修改日志:
 * <table>
 * <tr><th> Date        <th> Version  <th> Author  <th> Description
 * <tr><td> 2023-05-27  <td> 1.0      <td> Wisces  <td> 加载网格部分跑通
 * <tr><td> 2023-06-12  <td> 1.1      <td> Wisces  <td> 通量计算部分跑通
 * <tr><td> 2023-06-13  <td> 1.2      <td> Wisces  <td> 增加了计算时间步和 LUSGS 求解部分，并进行循环迭代 @line ../version_iteration/version1.0
 * <tr><td> 2023-07-07  <td> 1.3      <td> Wisces  <td> 增加了计算梯度部分
 * <tr><td> 2023-06-15  <td> 2.0      <td> Wisces  <td> 添加头文件 para_field_global.h，采用全局变量形式存储数据
 * <tr><td> 2023-06-18  <td> 3.0      <td> Wisces  <td> 添加头文件 parallel_bash.h，以及并行传值和计算边界虚网格的函数，支持 MPI 并行
 * <tr><td> 2023-07-03  <td> 4.0      <td> Wisces  <td> 添加实现了四种共享内存并行编程模型的算法，支持 OpenMP 并行
 * <tr><td> 2023-07-12  <td> 5.0      <td> Wisces  <td> 添加了 TDTree 线程级并行（使用条件编译）
 * </table>
 */

//!< C++ build-in head files
#include <iostream>
#include <sys/time.h> ///< gettimeofday()
#include <cmath>      ///< sqrt()

//!< user defined head files
#include "grid_polyhedra.h"
#include "para_field_global.h"
#include "parallel_base.h"

#include "load_grid.h"
#include "cal_flux.h"
#include "compute_timestep.h"
#include "forward_lusgs.h"
#include "cal_gradient.h"

//!< head file relying on condition-compiling
#ifdef MF_MPICH
#include <mpi.h>
#endif

#if (defined MF_OPENMP) || (defined TDTREE)
#include <omp.h>
#endif

#ifdef MF_MPICH
extern int myid;
extern int numprocs;
extern int myZone;        // = myid+1
extern MPI_Comm GridComm; ///< for each grid
extern RealFlow lusgs_comm;
extern RealFlow lusgs_cal;
#endif

using namespace std;

int main(int argc, char *argv[])
{
#ifdef MF_MPICH
        //!< 初始化 MPI 并行环境
        Parallel::InitMPI(argc, argv);

        if (myid == 0)
#endif
        {
                cout << endl;
                cout << "$#==================================================================================================#$" << endl;
                cout << "$#                                   National Numerical Windtunnel                                  #$" << endl;
                cout << "$#          MiniFS -- Mini-application for Fluxes-calculation and LUSGS-solver of FlowStar          #$" << endl;
                cout << "$#   Institute for Quantum Information & State Key Laboratory of High Performance Computing(HPCL)   #$" << endl;
                cout << "$#                          Parallel Computing Application Technology Lab.                          #$" << endl;
                cout << "$#                    National University of Defense Technology, Changsha, CHINA                    #$" << endl;
                cout << "$#                                      C.All rights reserved.                                      #$" << endl;
                cout << "$#==================================================================================================#$" << endl;
                cout << endl;
        }

#if (defined MF_OPENMP) || (defined TDTREE)
        cout << endl
             << "num_threads = " << omp_get_max_threads() << endl;
#endif

        //!< 加载网格
        PolyGrid *grid = new PolyGrid();
        BCond *bc = new BCond();
        LoadGrid(grid, bc);       ///< 读取网格文件
        Decoupling(grid);         ///< OpenMP 的算法实现
        PreprocessGrid(grid, bc); ///< 网格数据预处理

        //!< 无粘通量计算 和 LUSGS 求解
        /*!
         * @note steady = 1 定常气动力计算
         * @note vis_mode = INVISCID 暂时只关注无粘
         * @note 迭代步数
         *              总的物理迭代步 Unst_steps = 1
         *              第一个物理迭代步的子迭代步数 n_steps = 100（没有第二个及以后的物理迭代步数的子迭代步数 n_steps_unst）
         *              省略最外层关于 物理迭代步 的外循环（n_unst: 0），只关注 子迭代不 的内循环（inner_iter: 0-100）
         *
         */

#ifdef MF_TIMING
        double time_tmp = 0.0, time_tmp_total = 0.0, total_time = 0.0;
        int cvt = 1000000;
#ifdef MF_MPICH
        cvt = 1;

        MPI_Barrier(MPI_COMM_WORLD);
        time_tmp_total = -MPI_Wtime();
#else
        struct timeval starttimeTemRoe, endtimeTemRoe, starttimeTemRoe_Total, endtimeTemRoe_Total;
        gettimeofday(&starttimeTemRoe_Total, 0);

#endif
#endif

        IntType n_steps_inner = n_steps;

        //!< 相当于是内循环
        for (IntType inner_iter = 0; inner_iter < n_steps_inner; ++inner_iter)
        {
                #ifdef MF_MPICH
                                if (myid == 0)
                #endif
                                        cout << endl
                                             << "This is the " << inner_iter << "th step!" << endl;

                RealFlow **limit = NULL;
                IntType level = 0;

                //!< 残差置零（或在第一次迭代时初始化）
                ZeroResiduals(grid);

#ifdef MF_TIMING
#ifdef MF_MPICH
                MPI_Barrier(MPI_COMM_WORLD);
                time_tmp = -MPI_Wtime();
#else
                gettimeofday(&starttimeTemRoe, 0);
#endif
#endif
                //!< 计算无粘通量（原始代码中是更新残差 UpdateResiduals(grid, level)，这里忽略粘性通量的计算，只考虑无粘通量）
#ifdef TDTREE
                InviscidFlux_TDTree(grid, limit, level);
#else
                InviscidFlux(grid, limit, level);
#endif
#ifdef MF_TIMING
#ifdef MF_MPICH
                total_t[0] += time_tmp + MPI_Wtime();
#else
                gettimeofday(&endtimeTemRoe, 0);
                time_tmp = (RealGeom)1000000 * (endtimeTemRoe.tv_sec - starttimeTemRoe.tv_sec) + endtimeTemRoe.tv_usec - starttimeTemRoe.tv_usec;
                total_t[0] += time_tmp;
#endif
#endif

#ifdef MF_TIMING
#ifdef MF_MPICH
                MPI_Barrier(MPI_COMM_WORLD);
                time_tmp = -MPI_Wtime();
#else
                gettimeofday(&starttimeTemRoe, 0);
#endif
#endif
                //!< 计算时间步，用于后面 LUSGS 求解
                ComputeTimeStep(grid);
#ifdef MF_TIMING
#ifdef MF_MPICH
                total_t[3] += time_tmp + MPI_Wtime();
#else
                gettimeofday(&endtimeTemRoe, 0);
                time_tmp = (RealGeom)1000000 * (endtimeTemRoe.tv_sec - starttimeTemRoe.tv_sec) + endtimeTemRoe.tv_usec - starttimeTemRoe.tv_usec;
                total_t[3] += time_tmp;
#endif
#endif

#ifdef MF_TIMING
#ifdef MF_MPICH
                MPI_Barrier(MPI_COMM_WORLD);
                time_tmp = -MPI_Wtime();
#else
                gettimeofday(&starttimeTemRoe, 0);
#endif
#endif
                //!< 使用默认方法进行 LUSGS 求解
                ForwardLUSGS(grid, level);
#ifdef MF_TIMING
#ifdef MF_MPICH
                total_t[1] += time_tmp + MPI_Wtime();
#else
                gettimeofday(&endtimeTemRoe, 0);
                time_tmp = (RealGeom)1000000 * (endtimeTemRoe.tv_sec - starttimeTemRoe.tv_sec) + endtimeTemRoe.tv_usec - starttimeTemRoe.tv_usec;
                total_t[1] += time_tmp;
#endif
#endif
                //!< 得到新的流场变量后，并行传值、设置边界值、计算梯度、预处理计算，计算粘性
                //!< 并行传值
                PassInterfaceData(grid);

                //!< 计算边界虚网格
                // SetGhostVariables(grid);

                //!< 计算梯度
                CalculateGradient(grid);

                //!< 计算层流动力粘性

                iter_done++;
        }

        //!< 输出残差
        IntType nTCell = grid->GetNTCell();
        RealFlow *res[5];
        res[0] = grid->GetRes();
        // RealFlow *_res[5];
        // _res[0] = (RealFlow *)res;
        for (IntType i = 1; i < 5; i++)
                res[i] = &res[i - 1][nTCell];
        RealFlow norm[5] = {0., 0., 0., 0., 0.};
        for (IntType j = 0; j < 5; j++)
        {
                for (IntType i = 0; i < nTCell; i++)
                {
                        norm[j] += res[j][i] * res[j][i];
                }
        }
#ifdef MF_MPICH
        RealFlow norm_t[5];
        // MPI_Allreduce(norm, norm_t, 5, MPIReal, MPI_SUM, MPI_COMM_WORLD);
        MPI_Reduce(norm, norm_t, 5, MPIReal, MPI_SUM, 0, MPI_COMM_WORLD);
        for (IntType j = 0; j < 5; j++)
                norm[j] = norm_t[j];
#endif
        for (IntType j = 0; j < 5; j++)
                norm[j] = sqrt(norm[j]);

        String char_tmp;
        sprintf(char_tmp, "%5d %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n",
                (int)n_steps, norm[0], norm[1], norm[2], norm[3], norm[4], dt_max, dt_min);

#ifdef MF_MPICH
        if (myid == 0)
#endif
        {
                printf("\n#iter   rho_res      mx_res      my_res      mz_res     et_res     dt_max    dt_min\n");
                printf(char_tmp);
        }

        //!< 统计时间
#ifdef MF_TIMING
#ifdef MF_MPICH
// RealFlow  lusgs_time[2];
        struct 
        {
                RealFlow time;
                IntType id;
        } in[1],out[1];
        in[0].time = lusgs_comm + lusgs_cal;
        in[0].id = myid;
        MPI_Allreduce(in, out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        Parallel::parallel_max(total_t, 5, MPI_COMM_WORLD);
// Parallel::parallel_max(&lusgs_comm, 1, MPI_COMM_WORLD);
        // Parallel::parallel_max(&lusgs_cal, 1, MPI_COMM_WORLD);
        // if (myid == 0) cout << "LUSGS time = "<< lusgs_cal <<endl;
        if (myid == out[0].id) cout << "LUSGS COMM time = "<< lusgs_comm <<endl <<"LUSGS calculate time = " <<lusgs_cal <<endl;
        time_tmp_total = time_tmp_total + MPI_Wtime();
        MPI_Reduce(&time_tmp_total, &total_time, 1, MPIReal, MPI_MAX, 0, MPI_COMM_WORLD);

        if (myid == 0)
#else
        gettimeofday(&endtimeTemRoe_Total, 0);
        total_time = (RealGeom)1000000 * (endtimeTemRoe_Total.tv_sec - starttimeTemRoe_Total.tv_sec) + endtimeTemRoe_Total.tv_usec - starttimeTemRoe_Total.tv_usec;
#endif
        {
                cout << endl;
                cout << total_t[0] / cvt << endl;
                cout << total_t[1] / cvt << endl;
                cout << total_t[2] / cvt << endl;
                cout << "The total_t_inviscid = " << total_t[0] / cvt << endl;
                cout << "The total_t_lusgs = " << total_t[1] / cvt << endl;
                cout << "The total_t_gradient = " << total_t[2] / cvt << endl;
                cout << "The total_t_timestep = " << total_t[3] / cvt << endl;
                cout << "The total_t_comm = " << total_t[4] / cvt << endl;
                cout << "The total_t_123 = " << (total_t[0] + total_t[1] + total_t[2]) / cvt << endl;
                cout << "The total time = " << total_time / cvt << endl;
        }
#endif

        //!< 释放内存
        // sdel_array_1D(nodesymm);
        // sdel_array_1D(xfn_n_symm);
        // sdel_array_1D(yfn_n_symm);
        // sdel_array_1D(zfn_n_symm);
        // sdel_array_1D(norm_dist_c2c);
        // sdel_array_1D(LUSGSLayer);
        // sdel_array_1D(LUSGSCellOrder);
        // sdel_array_1D(LUSGScellsPerlayer);
        // sdel_array_1D(rho);
        // sdel_array_1D(u);
        // sdel_array_1D(v);
        // sdel_array_1D(w);
        // sdel_array_1D(p);
        // sdel_array_1D(dqdx);
        // sdel_array_1D(dqdy);
        // sdel_array_1D(dqdz);
        // sdel_array_1D(dt);
        // sdel_array_1D(res);
        // sdel_array_1D(DQ);
        sdel_object(grid);
        sdel_object(bc);

#ifdef MF_MPICH
        Parallel::FinalMPI();
#endif

        cout << "SUCCESSFUL!" << endl;
        return 0;
}
