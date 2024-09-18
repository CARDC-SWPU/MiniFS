/*!
 * @file        forward_lusgs.cpp
 * @brief       LUSGS 求解的主要函数
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @version     1.0
 * @date        2023-06-13
 * @copyright   Copyright (c) 2023  PERSONAL
 *
 * @par 修改日志:
 * <table>
 * <tr><th> Date        <th> Version  <th> Author  <th> Description
 * <tr><td> 2023-06-13  <td> 1.0      <td> Wisces  <td> 〔内容〕
 * <tr><td> 2023-06-15  <td> 2.0      <td> Wisces  <td> 添加头文件 para_field_global.h，采用全局变量形式存储数据
 * <tr><td> 2023-06-27  <td> 3.0      <td> Wisces  <td> 添加了 MPI 并行（使用条件编译）
 * <tr><td> 2023-07-06  <td> 4.0      <td> Wisces  <td> 添加了 OpenMP 并行（使用条件编译）
 * <tr><td> 2023-07-12  <td> 5.0      <td> Wisces  <td> 添加了 TDTree 线程级并行（使用条件编译）
 * </table>
 */

//!< C/C++ head files
#include <iostream>
// #include <cassert>
#include <cmath>

//!< direct head file
#include "forward_lusgs.h"

//!< user defined head files
#include "grid_polyhedra.h"
#include "memory_util.h"
#include "para_field_global.h"
#include "parallel_base.h"

//!< head file relying on condition-compiling
#if (defined MF_OPENMP) || (defined TDTREE)
#include <omp.h> ///< 可省略
#endif

#ifdef TDTREE
#include "TDTree.h"
#endif

#ifdef MF_MPICH
extern int myid;
extern int numprocs;
extern int myZone;        // = myid+1
extern MPI_Comm GridComm; ///< for each grid
RealFlow lusgs_cal = 0;
RealFlow lusgs_comm = 0 ;
#endif

using namespace std;

/*!
 * @brief       Forward the solution one time step using LU-SGS method.
 * @param       grid
 * @param       level       (= 0)
 * @remarks     modify according to the fun [void NSSolver::ForwardLUSGS(PolyGrid *grid, IntType level)]
 * @note        包含LU-SGS的改进型，单重网格单sweep的预估-校正型LU-SGS比原始的LU-SGS计算效率明显要高，别的情况下，效率提高不明显，甚至效率下降。
 *              参考文献：赵信文，预估-校正LU-SGS的隐式算法，航空计算技术第42卷第4期，2012年
 *              方法：在原始LU-SGS的基础上，增加一个校正步，将原始LU-SGS省略的高阶项L*(D-1)*U*(DQ)加进来。目前的经验是，多步LUSGS迭代收敛速度最快，在无多重网格时，建议sweeps取4，有多重网格时，建议sweeps取3收敛最快，且残差下降最好。epsilon建议取0.01~0.1之间，一般情况下取0.05收敛更快。
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-13
 */
void ForwardLUSGS(PolyGrid *grid, IntType level)
{
    IntType nTCell = grid->GetNTCell();
    IntType nBFace = grid->GetNBFace();
    IntType n = nTCell + nBFace;
    RealGeom *vol = grid->GetCellVol();
    RealFlow *res = grid->GetRes();
    RealFlow *dt = grid->GetDt();
    // RealFlow *DQ = grid->GetDQ();

    // RealFlow epsilon = 0.1;
    // grid->GetData(&epsilon, REAL_FLOW, 1, "epsilon");
    // if (epsilon < TINY)
    //     epsilon = 0.1;

    // Get number of faces for each cell
    IntType *nFPC = CalnFPC(grid);
    // Get cell to face conections
    IntType **C2F = CalC2F(grid);

    IntType i, j, ntemp;
    // Now diagonal term in LU-SGS, here we need information of time steps
    RealFlow *Diag = NULL;
    snew_array_1D(Diag, nTCell);
    // assert(Diag != 0);
    // 未修改 overlap

    for (i = 0; i < nTCell; i++)
    {
        Diag[i] = vol[i] / dt[i];
        // if(i > nTCell - 100)
        //     cout << Diag[i] << "\t";
    }

    // Note: As it has been shown, Diag = CFL/2*Vol/Dt.
    //       If function CalDiagLUSGS is not called, make sure CFL <= 2.
    // 未修改 overlap
#ifdef TDTREE
    CalDiagLUSGS_TDTree(grid, Diag, level);
#else
    CalDiagLUSGS(grid, Diag, level);
#endif

    // Allocate memories for RHS or DQ
    RealFlow *DQ[5];
    DQ[0] = grid->GetDQ();
    if (!DQ[0])
    {
        snew_array_1D(DQ[0], 5 * n);
        grid->SetDQ(DQ[0]);
    }
    // RealFlow *_DQ[5];
    // _DQ[0] = (RealFlow *)DQ; // grid->GetDataPtr(REAL_FLOW, 5 * n, "DQ");
    // if (!_DQ[0])
    // {
    //     snew_array_1D(_DQ[0], 5 * n);
    //     // grid->UpdateDataPtr(DQ[0], REAL_FLOW, 5 * n, "DQ");
    // }
    // assert(_DQ[0] != 0);
    for (i = 1; i < 5; i++)
        DQ[i] = &DQ[i - 1][n];
    for (j = 0; j < 5 * n; j++)
        DQ[0][j] = 0.;

    ///< Wisces: 我们只关注 sweeps == 1 的情况
    /**
    if (sweeps == -1) {}
    else if (sweeps == -2) {}
    else if (sweeps == 1)
    */
    {
        //!< 单步
        //!< Copy the residual to DQ
        ntemp = 0;
        for (i = 0; i < 5; i++)
        {
            for (j = 0; j < nTCell; j++)
            {
                DQ[i][j] = res[ntemp++];
            }
        }
        // Now the LU-SGS part
#ifdef TDTREE
        SolveLUSGS3D_TDTree(grid, Diag, DQ, nFPC, C2F, level);
#else
        SolveLUSGS3D(grid, Diag, DQ, nFPC, C2F, level);
#endif
    }
    /**
    else {}
    */
    //!< Update flow field
    UpdateFlowField3D_CFL3d(grid, DQ);

    //!< delete temporary memories
    sdel_array_1D(Diag);
}

/*!
 * @brief       Calculate diagonal term in LU-SGS.
 * @param       grid
 * @param       Diag
 * @param       level
 * @remarks     modify according to the fun [void CalDiagLUSGS(PolyGrid *grid, RealFlow *Diag, IntType level)]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-13
 */
void CalDiagLUSGS(PolyGrid *grid, RealFlow *Diag, IntType level)
{
    IntType nTCell = grid->GetNTCell();
    IntType nBFace = grid->GetNBFace();
    IntType nTFace = grid->GetNTFace();
    IntType n = nTCell + nBFace;
    IntType *f2c = grid->Getf2c();
    RealGeom *xfn = grid->GetXfn();
    RealGeom *yfn = grid->GetYfn();
    RealGeom *zfn = grid->GetZfn();
    RealGeom *area = grid->GetFaceArea();
    RealGeom *vgn = grid->GetFaceNormalVelocity();
    RealGeom *vol = grid->GetCellVol();

    RealFlow *rho = grid->GetRho();
    RealFlow *u = grid->GetU();
    RealFlow *v = grid->GetV();
    RealFlow *w = grid->GetW();
    RealFlow *p = grid->GetP();

    // RealFlow lhs_omga;
    // lhs_omga = 1.0000000e+00;
    // grid->GetData(&lhs_omga, REAL_FLOW, 1, "lhs_omga");

    IntType i, c1, c2, count;
    RealFlow vn_1, vn_2, ss_1, ss_2;
    RealFlow eig, temp;

    // assert(norm_dist_c2c); // must exist

#if (defined MF_OPENMP) && (defined OMP_GroupColor)
    if (grid->GroupColorSuccess)
    {
        IntType pfacenum = nBFace - grid->GetNIFace();
        IntType groupSize = grid->groupSize;
        IntType bfacegroup_num, ifacegroup_num;
        IntType startFace, endFace, count;
        bfacegroup_num = grid->bfacegroup.size();
        ifacegroup_num = grid->ifacegroup.size();

        //!< Physical faces:
        for (IntType fcolor = 0; fcolor < bfacegroup_num; fcolor++)
        {
            if (!fcolor)
            {
                startFace = 0;
            }
            else
            {
                startFace = grid->bfacegroup[fcolor - 1];
            }
            endFace = grid->bfacegroup[fcolor];
#pragma omp parallel for private(i, c1, vn_1, ss_1, eig) schedule(static, groupSize)
            for (i = startFace; i < endFace; i++)
            {
                c1 = f2c[2 * i];
                vn_1 = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
                // if (!steady)
                //     vn_1 -= vgn[i];
                vn_1 = fabs(vn_1);
                ss_1 = gam * (p[c1] + p_bar) / rho[c1];
                eig = vn_1 + sqrt(ss_1);

                Diag[c1] += 0.5 * area[i] * eig * lhs_omga; ///< 对角线的表达式
            }
        }

        //!< Iterface
#ifdef MF_MPICH
        count = 2 * pfacenum;
        for (i = pfacenum; i < nBFace; i++)
        {
            c1 = f2c[2 * i];
            vn_1 = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
            // if (!steady)
            //     vn_1 -= vgn[i];
            vn_1 = fabs(vn_1);
            ss_1 = gam * (p[c1] + p_bar) / rho[c1];
            eig = vn_1 + sqrt(ss_1);

            Diag[c1] += 0.5 * area[i] * eig * lhs_omga; ///< 对角线的表达式
        }
#endif

        //!< Interior faces
        for (IntType fcolor = 0; fcolor < ifacegroup_num; fcolor++)
        {
            if (!fcolor)
            {
                startFace = nBFace;
            }
            else
            {
                startFace = grid->ifacegroup[fcolor - 1];
            }
            endFace = grid->ifacegroup[fcolor];
#pragma omp parallel for private(i, count, c1, c2, vn_1, ss_1, eig, vn_2, ss_2) schedule(static, groupSize)
            for (i = startFace; i < endFace; i++)
            {
                count = 2 * i;
                c1 = f2c[count++];
                c2 = f2c[count++];

                //!< Cell c1
                vn_1 = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
                // if (!steady)
                //     vn_1 -= vgn[i];
                vn_1 = fabs(vn_1);
                ss_1 = gam * (p[c1] + p_bar) / rho[c1];
                eig = vn_1 + sqrt(ss_1);

                Diag[c1] += 0.5 * area[i] * eig * lhs_omga; ///< 对角线的表达式

                //!< Cell c2
                vn_2 = u[c2] * xfn[i] + v[c2] * yfn[i] + w[c2] * zfn[i];
                // if (!steady)
                //     vn_2 -= vgn[i];
                vn_2 = fabs(vn_2);
                ss_2 = gam * (p[c2] + p_bar) / rho[c2];
                eig = vn_2 + sqrt(ss_2);

                Diag[c2] += 0.5 * area[i] * eig * lhs_omga; ///< 对角线的表达式
            }
        }
    }
    else
    {
        //!< Boundary faces first
        for (i = 0; i < nBFace; i++)
        {
            c1 = f2c[2 * i];
            vn_1 = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
            // if (!steady)
            //     vn_1 -= vgn[i];
            vn_1 = fabs(vn_1);
            ss_1 = gam * (p[c1] + p_bar) / rho[c1];
            eig = vn_1 + sqrt(ss_1);

            Diag[c1] += 0.5 * area[i] * eig * lhs_omga; ///< 对角线的表达式
        }
        //!< Interior faces
        count = 2 * nBFace;
        for (i = nBFace; i < nTFace; i++)
        {
            c1 = f2c[count++];
            c2 = f2c[count++];

            //!< Cell c1
            vn_1 = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
            // if (!steady)
            //     vn_1 -= vgn[i];
            vn_1 = fabs(vn_1);
            ss_1 = gam * (p[c1] + p_bar) / rho[c1];
            eig = vn_1 + sqrt(ss_1);

            Diag[c1] += 0.5 * area[i] * eig * lhs_omga; ///< 对角线的表达式

            //!< Cell c2
            vn_2 = u[c2] * xfn[i] + v[c2] * yfn[i] + w[c2] * zfn[i];
            // if (!steady)
            //     vn_2 -= vgn[i];
            vn_2 = fabs(vn_2);
            ss_2 = gam * (p[c2] + p_bar) / rho[c2];
            eig = vn_2 + sqrt(ss_2);

            Diag[c2] += 0.5 * area[i] * eig * lhs_omga; ///< 对角线的表达式
        }
    }

#elif (defined MF_OPENMP) && (defined OMP_FaceColor)
    IntType nIFace = grid->GetNIFace();
    IntType pfacenum = nBFace - nIFace;
    IntType bfacegroup_num, ifacegroup_num;
    ifacegroup_num = (*grid).ifacegroup.size(); ///< 内部面的颜色数
    bfacegroup_num = (*grid).bfacegroup.size(); ///< 物理边界面的颜色数

    //!< Boundary faces:
    for (IntType fcolor = 0; fcolor < bfacegroup_num; fcolor++)
    {
        IntType startFace, endFace;
        if (fcolor == 0)
        {
            startFace = 0;
        }
        else
        {
            startFace = (*grid).bfacegroup[fcolor - 1];
        }
        endFace = (*grid).bfacegroup[fcolor];

#pragma omp parallel for private(i, c1, vn_1, ss_1, eig)
        for (i = startFace; i < endFace; i++)
        {
            c1 = f2c[2 * i];
            vn_1 = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
            // if (!steady)
            //     vn_1 -= vgn[i];
            vn_1 = fabs(vn_1);
            ss_1 = gam * (p[c1] + p_bar) / rho[c1];
            eig = vn_1 + sqrt(ss_1);

            Diag[c1] += 0.5 * area[i] * eig * lhs_omga; ///< 对角线的表达式
        }
    }

    //!< Interface:
#ifdef MF_MPICH
    for (IntType i = pfacenum; i < nBFace; i++)
    {
        count = 2 * pfacenum;
        for (i = pfacenum; i < nBFace; i++)
        {
            c1 = f2c[2 * i];
            vn_1 = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
            // if (!steady)
            //     vn_1 -= vgn[i];
            vn_1 = fabs(vn_1);
            ss_1 = gam * (p[c1] + p_bar) / rho[c1];
            eig = vn_1 + sqrt(ss_1);

            Diag[c1] += 0.5 * area[i] * eig * lhs_omga; ///< 对角线的表达式
        }
    }
#endif

    //!< Interior faces:
    for (IntType fcolor = 0; fcolor < ifacegroup_num; fcolor++)
    {
        IntType startFace, endFace;
        if (fcolor == 0)
        {
            startFace = nBFace;
        }
        else
        {
            startFace = (*grid).ifacegroup[fcolor - 1];
        }
        endFace = (*grid).ifacegroup[fcolor];

#pragma omp parallel for private(i, count, c1, c2, vn_1, ss_1, eig, vn_2, ss_2)
        for (i = startFace; i < endFace; i++)
        {
            count = 2 * i;
            c1 = f2c[count++];
            c2 = f2c[count++];

            //!< Cell c1
            vn_1 = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
            // if (!steady)
            //     vn_1 -= vgn[i];
            vn_1 = fabs(vn_1);
            ss_1 = gam * (p[c1] + p_bar) / rho[c1];
            eig = vn_1 + sqrt(ss_1);

            Diag[c1] += 0.5 * area[i] * eig * lhs_omga; ///< 对角线的表达式

            //!< Cell c2
            vn_2 = u[c2] * xfn[i] + v[c2] * yfn[i] + w[c2] * zfn[i];
            // if (!steady)
            //     vn_2 -= vgn[i];
            vn_2 = fabs(vn_2);
            ss_2 = gam * (p[c2] + p_bar) / rho[c2];
            eig = vn_2 + sqrt(ss_2);

            Diag[c2] += 0.5 * area[i] * eig * lhs_omga; ///< 对角线的表达式
        }
    }

#elif (defined MF_OPENMP) && (defined OMP_Reduction)
    IntType *nFPC = CalnFPC(grid);
    IntType **C2F = CalC2F(grid);
    IntType j, face;
#pragma omp parallel for private(j, face, vn_1, ss_1, eig)
    for (i = 0; i < nTCell; i++)
    {
        for (j = 0; j < nFPC[i]; j++)
        {
            face = C2F[i][j];
            vn_1 = u[i] * xfn[face] + v[i] * yfn[face] + w[i] * zfn[face];
            // if (!steady)
            //     vn_1 -= vgn[face];
            vn_1 = fabs(vn_1);
            ss_1 = gam * (p[i] + p_bar) / rho[i];
            eig = vn_1 + sqrt(ss_1);
            Diag[i] += 0.5 * area[face] * eig * lhs_omga;
        }
    }
#elif (defined MF_OPENMP) && (defined OMP_DIVREP)
    IntType threads = grid->threads;
    IntType startFace, endFace, t, k, face;
    if (grid->DivRepSuccess)
    {
#pragma omp parallel for private(i, k, startFace, endFace, face, c1, c2, count, vn_1, ss_1, eig, vn_2, ss_2)
        for (t = 0; t < threads; t++)
        {
            //!< Boundary faces
            startFace = grid->idx_pthreads_bface[t];
            endFace = grid->idx_pthreads_bface[t + 1];
            for (i = startFace; i < endFace; i++)
            {
                face = grid->id_division_bface[i];
                c1 = f2c[2 * face];
                vn_1 = u[c1] * xfn[face] + v[c1] * yfn[face] + w[c1] * zfn[face];
                // if (!steady)
                //     vn_1 -= vgn[face];
                vn_1 = fabs(vn_1);
                ss_1 = gam * (p[c1] + p_bar) / rho[c1];
                eig = vn_1 + sqrt(ss_1);

                Diag[c1] += 0.5 * area[face] * eig * lhs_omga;
            }
            //!< Interior faces
            startFace = grid->idx_pthreads_iface[t];
            endFace = grid->idx_pthreads_iface[t + 1];
            for (i = startFace; i < endFace; i++)
            {
                k = grid->id_division_iface[i];
                if (abs(k) < nTFace)
                    face = k;
                else
                    face = abs(k) - nTFace;
                count = 2 * face;
                c1 = f2c[count];
                c2 = f2c[count + 1];

                if (abs(k) < nTFace)
                {
                    //!< write back to c1 & c2
                    //!< Cell c1
                    vn_1 = u[c1] * xfn[face] + v[c1] * yfn[face] + w[c1] * zfn[face];
                    if (!steady)
                        vn_1 -= vgn[face];
                    vn_1 = fabs(vn_1);
                    ss_1 = gam * (p[c1] + p_bar) / rho[c1];
                    eig = vn_1 + sqrt(ss_1);

                    Diag[c1] += 0.5 * area[face] * eig * lhs_omga; //!< 对角线的表达式

                    //!< Cell c2
                    vn_2 = u[c2] * xfn[face] + v[c2] * yfn[face] + w[c2] * zfn[face];
                    // if (!steady)
                    //     vn_2 -= vgn[face];
                    vn_2 = fabs(vn_2);
                    ss_2 = gam * (p[c2] + p_bar) / rho[c2];
                    eig = vn_2 + sqrt(ss_2);

                    Diag[c2] += 0.5 * area[face] * eig * lhs_omga; //!< 对角线的表达式
                }
                else
                {
                    if (k > 0)
                    {
                        //!< just write back to c1
                        vn_1 = u[c1] * xfn[face] + v[c1] * yfn[face] + w[c1] * zfn[face];
                        // if (!steady)
                        //     vn_1 -= vgn[face];
                        vn_1 = fabs(vn_1);
                        ss_1 = gam * (p[c1] + p_bar) / rho[c1];
                        eig = vn_1 + sqrt(ss_1);

                        Diag[c1] += 0.5 * area[face] * eig * lhs_omga; ///< 对角线的表达式
                    }
                    else
                    {
                        //!< just write back to c2
                        vn_2 = u[c2] * xfn[face] + v[c2] * yfn[face] + w[c2] * zfn[face];
                        // if (!steady)
                        //     vn_2 -= vgn[face];
                        vn_2 = fabs(vn_2);
                        ss_2 = gam * (p[c2] + p_bar) / rho[c2];
                        eig = vn_2 + sqrt(ss_2);

                        Diag[c2] += 0.5 * area[face] * eig * lhs_omga; ///< 对角线的表达式
                    }
                }
            }
        }
    }

#else
    //!< Boundary faces first
    for (i = 0; i < nBFace; i++)
    {
        c1 = f2c[2 * i];
        vn_1 = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
        // if (!steady)
        //     vn_1 -= vgn[i];
        vn_1 = fabs(vn_1);
        ss_1 = gam * (p[c1] + p_bar) / rho[c1];
        eig = vn_1 + sqrt(ss_1);

        Diag[c1] += 0.5 * area[i] * eig * lhs_omga; ///< 对角线的表达式
    }
    //!< Interior faces
    count = 2 * nBFace;
    for (i = nBFace; i < nTFace; i++)
    {
        c1 = f2c[count++];
        c2 = f2c[count++];

        //!< Cell c1
        vn_1 = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
        // if (!steady)
        //     vn_1 -= vgn[i];
        vn_1 = fabs(vn_1);
        ss_1 = gam * (p[c1] + p_bar) / rho[c1];
        eig = vn_1 + sqrt(ss_1);

        Diag[c1] += 0.5 * area[i] * eig * lhs_omga; ///< 对角线的表达式

        //!< Cell c2
        vn_2 = u[c2] * xfn[i] + v[c2] * yfn[i] + w[c2] * zfn[i];
        // if (!steady)
        //     vn_2 -= vgn[i];
        vn_2 = fabs(vn_2);
        ss_2 = gam * (p[c2] + p_bar) / rho[c2];
        eig = vn_2 + sqrt(ss_2);

        Diag[c2] += 0.5 * area[i] * eig * lhs_omga; ///< 对角线的表达式
    }
#endif ///< ~MF_OPENMP

    //!< Wisces: 只考虑定常，steady == 1
    /**
    if (!steady)
    {
#ifdef MF_OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < nTCell; i++)
            Diag[i] += (1.0 + time_accuracy) * vol[i] / real_dt;
    }
    */

    // If flow is viscous, need to count the contribution from viscosity
    // 该程序的粘性增加预处理的需要修改,然后进行测试
    // IntType vis_mode;
    // IntType vis_run = 0;
    // vis_mode = INVISCID;
    // grid->GetData(&vis_mode, INT, 1, "vis_mode");
    //!< Wisces: 暂时只考虑无粘，vis_mode == INVISCID
    /**
    if (vis_mode != INVISCID)
    {
        /// Wisces: vis_mode == INVISCID, 未处理，直接注释掉
        vis_run = 1;
        // if coarse grid doesn't want to run the viscous flux, turn it off
        if (level != 0)
        {
            IntType cg_vis = 1;
            grid->GetData(&cg_vis, INT, 1, "cg_vis");
            if (cg_vis == 0)
                vis_run = 0;
        }
    }

    if (vis_run)
    {
        /// Wisces: vis_run == 0, 未处理，直接注释掉

        RealFlow dot;
        RealFlow *vis_l = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "vis_l");
        RealFlow *vis_t = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "vis_t");
        #ifdef MF_OPENMP
                RealFlow Diag_v;
        #if !(defined OMP_Reduction)
                IntType *nFPC = CalnFPC(grid);
                IntType **C2F = CalC2F(grid);
                IntType j, face;
        #endif
                RealFlow *tmp_diag = NULL;
                snew_array_1D(tmp_diag, nTFace);
        #pragma omp parallel for private(dot)
                for (i = 0; i < nTFace; i++)
                {
                    dot = norm_dist_c2c[i];
                    tmp_diag[i] = area[i] / (dot + TINY);
                }
        #pragma omp parallel for private(j, face, Diag_v)
                for (i = 0; i < nTCell; i++)
                {
                    Diag_v = 0;
                    for (j = 0; j < nFPC[i]; j++)
                    {
                        face = C2F[i][j];
                        Diag_v += tmp_diag[face];
                    }
                    Diag[i] += (vis_l[i] + vis_t[i]) * Diag_v / rho[i];
                }
                sdel_array_1D(tmp_diag);
        #else

        RealFlow *Diag_v = NULL;
        snew_array_1D(Diag_v, nTCell);
        assert(Diag_v != 0);
        for (i = 0; i < nTCell; i++)
            Diag_v[i] = 0.;

        count = 0;
        for (i = 0; i < nTFace; i++)
        {
            c1 = f2c[count++];
            c2 = f2c[count++];
            dot = norm_dist_c2c[i];
            temp = area[i] / (dot + TINY);
            // Cell c1
            Diag_v[c1] += temp;
            // Cell c2
            if (c2 < nTCell)
                Diag_v[c2] += temp;
        }

        for (i = 0; i < nTCell; i++)
            Diag[i] += (vis_l[i] + vis_t[i]) * Diag_v[i] / rho[i];
        sdel_array_1D(Diag_v);
        #endif

    }
    */
}

/*!
 * @brief       Solve linear systems using the LU-SGS in 3D ~~~ORIGINAL LU-SGS ONE SWEEP~~~
 * @param       grid
 * @param       Diag
 * @param       DQ
 * @param       nFPC
 * @param       C2F
 * @param       level
 * @remarks     modify according to the fun [void SolveLUSGS3D(PolyGrid *grid, RealFlow *Diag, RealFlow *DQ[5], IntType *nFPC, IntType **C2F, IntType level)]
 * @note        增加网格排序，在单重网格时有明显的加速效果，但是多重网格时加速不明显，甚至收敛变慢
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-13
 */
void SolveLUSGS3D(PolyGrid *grid, RealFlow *Diag, RealFlow *DQ[5], IntType *nFPC, IntType **C2F, IntType level)
{

    IntType nTCell = grid->GetNTCell();
    IntType nBFace = grid->GetNBFace();
    IntType n = nTCell + nBFace;
    IntType *f2c = grid->Getf2c();
    // Get grid metrics
    RealGeom *xfn = grid->GetXfn();
    RealGeom *yfn = grid->GetYfn();
    RealGeom *zfn = grid->GetZfn();
    RealGeom *area = grid->GetFaceArea();
    RealGeom *vgn = grid->GetFaceNormalVelocity();
    // Get flow variables
    RealFlow *q[5];
    // q[0] = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "rho");
    // q[1] = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "u");
    // q[2] = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "v");
    // q[3] = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "w");
    // q[4] = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "p");
    q[0] = (RealFlow *)grid->GetRho();
    q[1] = (RealFlow *)grid->GetU();
    q[2] = (RealFlow *)grid->GetV();
    q[3] = (RealFlow *)grid->GetW();
    q[4] = (RealFlow *)grid->GetP();

    // RealFlow lhs_omga;
    // lhs_omga = 1.0000000e+00;
    // grid->GetData(&lhs_omga, REAL_FLOW, 1, "lhs_omga");

    RealFlow rho00 = rho_s, u00 = u_s, v00 = v_s, w00 = w_s /*, e_stag*/;
    // grid->GetData(&rho00, REAL_FLOW, 1, "rho");
    // grid->GetData(&u00, REAL_FLOW, 1, "u");
    // grid->GetData(&v00, REAL_FLOW, 1, "v");
    // grid->GetData(&w00, REAL_FLOW, 1, "w");

    // RealFlow rho_min, rho_max, p_min, p_max;
    // RealFlow e_stag_max;

    // IntType DQ_limit = 1;
    // DQ_limit = 2;
    // grid->GetData(&DQ_limit, INT, 1, "DQ_limit");

    // IntType vis_mode;
    // IntType vis_run = 0;
    // vis_mode = INVISCID;
    // grid->GetData(&vis_mode, INT, 1, "vis_mode");
    /**
     if (vis_mode != INVISCID)
     {
         /// Wisces: vis_mode == INVISCID, 未处理，直接注释掉
         vis_run = 1;

         // if coarse grid doesn't want to run the viscous flux, turn it off
         if (level != 0)
         {
             IntType cg_vis = 1;
             grid->GetData(&cg_vis, INT, 1, "cg_vis");
             if (cg_vis == 0)
                 vis_run = 0;
         }
     }
    RealFlow *vis_l = NULL, *vis_t = NULL;
    if (vis_run)
    {
        /// Wisces: vis_run == 0, 未处理
        // vis_l = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "vis_l");
        // vis_t = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "vis_t");
    }
    */
    // Some temporary variables
    IntType i, j, ilu, face, cell, c1, c2, c_tmp, count;
    RealFlow flux_s[5], flux[5], q_loc[5], DQ_loc[5], visc, tmp, vgn_tmp;
    RealGeom face_n[3], dist;

    IntType *luorder = grid->GetLUSGSCellOrder();           // (IntType *) grid->GetDataPtr(INT, nTCell, "LUSGSCellOrder"); // 为了对角占优，网格排序后的序号
    IntType *layer = grid->GetLUSGSLayer();                 // (IntType *) grid->GetDataPtr(INT, n, "LUSGSLayer"); // LUSGS迭代的层号，层号小为下三角，层号大为上三角
    IntType *cellsPerlayer = grid->GetLUSGScellsPerlayer(); // (IntType *)grid->GetDataPtr(INT, nTCell, "LUSGScellsPerlayer");

    IntType nTFace = grid->GetNTFace();
    // RealGeom *norm_dist_c2c = NULL;
    // norm_dist_c2c = (RealGeom *)grid->GetDataPtr(REAL_GEOM, nTFace, "norm_dist_c2c");
    // assert(norm_dist_c2c); // must exist
#ifdef MF_TIMING
#ifdef MF_MPICH
MPI_Barrier(MPI_COMM_WORLD);
    double time_tmp = -MPI_Wtime();
#endif
#endif 
    /*!
     * @brief       Forward Sweep
     *
     */
    for (i = 0; i < 5; i++)
        DQ[i][luorder[0]] /= Diag[luorder[0]];
#if (defined MF_OPENMP) && (defined OMP_CellColor)
    //!< not containing SIMD
    IntType laynum;
    IntType start, end;
    for (laynum = 0; laynum < cellsPerlayer[0]; laynum++)
    {
        start = cellsPerlayer[laynum + 1];
        end = cellsPerlayer[laynum + 2];
        if (laynum == 0)
        {
            start++;
        }
#pragma omp parallel for private(ilu)
        for (ilu = start; ilu < end; ilu++)
        {
            IntType cell;
            cell = luorder[ilu];
#else
    for (ilu = 1; ilu < nTCell; ilu++)
    {
        cell = luorder[ilu];
        // cell = ilu;
#endif
            for (IntType j = 0; j < nFPC[cell]; j++)
            {
                IntType face, c1, c2, c_tmp, count;
                RealFlow flux_s[5], flux[5], q_loc[5], DQ_loc[5], visc, tmp, vgn_tmp;
                RealGeom face_n[3], dist;
                face = C2F[cell][j];
                count = face + face;
                c1 = f2c[count++];
                c2 = f2c[count];
                // One of c1 and c2 must be cell itself.
                if (layer[c1] > layer[cell] || layer[c2] > layer[cell])
                    // if (c1 > cell || c2 > cell)
                    continue;

                // Now its neighboring cell belongs to lower triangular
                face_n[0] = xfn[face];
                face_n[1] = yfn[face];
                face_n[2] = zfn[face];
                // if (!steady)
                //     vgn_tmp = vgn[face];
                if (c2 == cell)
                {
                    c_tmp = c1;
                    c1 = c2;
                    c2 = c_tmp;
                    face_n[0] = -face_n[0];
                    face_n[1] = -face_n[1];
                    face_n[2] = -face_n[2];
                    // if (!steady)
                    //     vgn_tmp = -vgn[face];
                }
                // assert(c1 == cell);

                for (IntType i = 0; i < 5; i++)
                {
                    q_loc[i] = q[i][c2];
                    DQ_loc[i] = DQ[i][c2];
                }

                // cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<"ilu:"<<ilu<<"    i:"<<i<<"DQ[i][cell]  "<<DQ[i][cell]<<endl;

                // Calculate everything (I call it Flux) in lower triangular
                //!< Wisces: 只考虑定常，steady ==1
                FluxLUSGS3D(flux, q_loc, DQ_loc, face_n, gam, p_bar, lhs_omga);
                /**
                if (steady)
                {
                    FluxLUSGS3D(flux, q_loc, DQ_loc, face_n, gam, p_bar, lhs_omga);
                }
                else
                {
                    FluxLUSGS3D_unsteady(flux, q_loc, DQ_loc, face_n, gam, p_bar, lhs_omga, vgn_tmp);
                }
                */

                /**
                if (vis_run)
                {
                    dist = norm_dist_c2c[face];
                    visc = vis_l[c2] + vis_t[c2];
                    tmp = 2.0 * visc / (q_loc[0] * dist + TINY);
                    for (IntType i = 0; i < 5; i++)
                        flux[i] -= tmp * DQ_loc[i];
                }
                */

                // if(ilu<1000)
                // cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<"ilu="<<ilu<<"............"<<"flux="<<flux[0]<<endl;
                // Add Flux together
                tmp = 0.5 * area[face];
                for (IntType i = 0; i < 5; i++)
                {
                    DQ[i][cell] -= tmp * flux[i];
                    // cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<"ilu:"<<ilu<<"    i:"<<i<<"DQ[i][cell]  "<<DQ[i][cell]<<endl;
                }
            }
            for (IntType i = 0; i < 5; i++)
                DQ[i][cell] /= Diag[cell];

            // mflog::log.set_all_processors_out();
            //         int rank_id = 0;
            // #ifdef MF_MPICH
            //         rank_id = myid;
            // #endif
            //         if (fabs(DQ[0][cell]) > 1.0e3 * rho00)
            //         {
            //             cout << "Forward sweep: drho>1.0e3*rho00!  " << DQ[0][cell] << " "
            //                  << cell << " " << rank_id << endl;
            //         }
            //         if (fabs(DQ[4][cell]) > 1.0e5 * e_stag)
            //         {
            //             cout << "Forward sweep: de>1.0e5*e_stag!  " << DQ[4][cell] << " "
            //                  << cell << " " << rank_id << endl;
            //         }
            //         if (fabs(DQ[0][cell]) > rho_max || DQ[4][cell] > e_stag_max)
            //         {
            //             cout << "Error!\n Maybe CFL too big or entropy correction coefficient too small!" << endl;
            //             // mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
            //             exit(1);
            //         }

            // limit for rho>0
            //!< Wisces: DQ_limit == 2
            /**
            if (DQ_limit == 1)
            {
                // do nothing!
            }
            else if (DQ_limit == 2)
            {
            */
            RealFlow dp, vv;
            vv = q[1][cell] * q[1][cell] + q[2][cell] * q[2][cell] + q[3][cell] * q[3][cell];
            dp = DQ[4][cell] + 0.5 * DQ[0][cell] * vv - (DQ[1][cell] * q[1][cell] + DQ[2][cell] * q[2][cell] + DQ[3][cell] * q[3][cell]);
            dp *= gamm1;
            if ((q[0][cell] + DQ[0][cell]) < rho_min || (q[0][cell] + DQ[0][cell]) > rho_max ||
                (q[4][cell] + dp) < p_min || (q[4][cell] + dp) > p_max)
            {
                DQ[0][cell] *= 0.1;
                DQ[1][cell] *= 0.1;
                DQ[2][cell] *= 0.1;
                DQ[3][cell] *= 0.1;
                DQ[4][cell] *= 0.1;
            }
            dp = DQ[4][cell] + 0.5 * DQ[0][cell] * vv - (DQ[1][cell] * q[1][cell] + DQ[2][cell] * q[2][cell] + DQ[3][cell] * q[3][cell]);
            dp *= gamm1;
            if ((q[0][cell] + DQ[0][cell]) < rho_min || (q[0][cell] + DQ[0][cell]) > rho_max ||
                (q[4][cell] + dp) < p_min || (q[4][cell] + dp) > p_max)
            {
                DQ[0][cell] *= 0.1;
                DQ[1][cell] *= 0.1;
                DQ[2][cell] *= 0.1;
                DQ[3][cell] *= 0.1;
                DQ[4][cell] *= 0.1;
            }
            dp = DQ[4][cell] + 0.5 * DQ[0][cell] * vv - (DQ[1][cell] * q[1][cell] + DQ[2][cell] * q[2][cell] + DQ[3][cell] * q[3][cell]);
            dp *= gamm1;
            if ((q[0][cell] + DQ[0][cell]) < rho_min || (q[0][cell] + DQ[0][cell]) > rho_max ||
                (q[4][cell] + dp) < p_min || (q[4][cell] + dp) > p_max)
            {
                DQ[0][cell] = 0.0;
                DQ[1][cell] = 0.0;
                DQ[2][cell] = 0.0;
                DQ[3][cell] = 0.0;
                DQ[4][cell] = 0.0;
            }
            /**
            }

            else if (DQ_limit == 3)
            {
                DQ[0][cell] = MAX(DQ[0][cell], rho_min - q[0][cell]);
                DQ[0][cell] = MIN(DQ[0][cell], rho_max - q[0][cell]);
            }
            else if (DQ_limit == 4)
            {
                RealFlow alph, alph_rho, alph_rhoe, alph_p, dp, vv, rhoe;
                vv = q[1][cell] * q[1][cell] + q[2][cell] * q[2][cell] + q[3][cell] * q[3][cell];
                rhoe = 0.5 * q[0][cell] * vv + (q[4][cell] + p_bar) / (gam - 1.0);
                dp = DQ[4][cell] + 0.5 * DQ[0][cell] * vv - (DQ[1][cell] * q[1][cell] + DQ[2][cell] * q[2][cell] + DQ[3][cell] * q[3][cell]);
                dp *= gamm1;

                alph_rho = q[0][cell] / (MAX(q[0][cell], 0.05 * rho00) + MAX(0.0, -DQ[0][cell]));
                alph_rhoe = rhoe / (MAX(rhoe, 0.05 * e_stag) + MAX(0.0, -DQ[4][cell]));
                alph_p = (q[4][cell] + p_bar) / (MAX((q[4][cell] + p_bar), 0.05 * p_bar) + MAX(0.0, -dp));
                alph = MIN(alph_rho, alph_rhoe);
                alph = MIN(alph, alph_p);
                for (IntType i = 0; i < 5; i++)
                    DQ[i][cell] *= alph;
            }
            else
            {
                // mflog::log.set_one_processor_out();
                // mflog::log << endl
                //            << "DQ_limit is greater to 4! Now only have 4 methods." << endl;
                // mflog::log << "Then we will use the first method, i.e. do nothing!" << endl;
                cout << endl
                     << "DQ_limit is greater to 4! Now only have 4 methods." << endl;
                cout << "Then we will use the first method, i.e. do nothing!" << endl;
            }
            */
#if (defined MF_OPENMP) && (defined OMP_CellColor)
        }
    }
#else
    }
#endif
#ifdef MF_TIMING
#ifdef MF_MPICH
MPI_Barrier(MPI_COMM_WORLD);
    lusgs_cal += time_tmp + MPI_Wtime();
    struct 
    {
            RealFlow time;
            IntType id;
    } in[1],out[1];
    in[0].time = lusgs_cal;
    in[0].id = myid;
    MPI_Allreduce(in, out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD); 
    lusgs_cal = out[0].time;   
#endif
#endif
#ifdef MF_MPICH
    IntType nvar = 5;
    RealFlow *q_mpi[5];
    for (IntType j = 0; j < 5; j++)
        q_mpi[j] = DQ[j];
MPI_Barrier(MPI_COMM_WORLD);
    time_tmp = -MPI_Wtime();
#ifdef GASPI
    grid->GASPIRecvSendVarNeighbor_Togeth(nvar, q_mpi);
#else
    grid->RecvSendVarNeighbor_Togeth(nvar, q_mpi);
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    lusgs_comm += time_tmp + MPI_Wtime();
    in[0].time = lusgs_comm;
    in[0].id = myid;
    MPI_Allreduce(in, out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD); 
    lusgs_comm = out[0].time;   
    // time_tmp = -MPI_Wtime();
    // for (IntType i = 0; i < nvar; ++i)
    // {
    //     grid->CommInterfaceDataMPI(q_mpi[i]);
    // }
#endif

    /*!
     * @brief       Backward Sweep
     *
     */

#if (defined MF_OPENMP) && (defined OMP_CellColor)
    //!< not containing SIMD
    for (laynum = cellsPerlayer[0] - 1; laynum >= 0; laynum--)
    {
        start = cellsPerlayer[laynum + 2];
        end = cellsPerlayer[laynum + 1];
        IntType ilu;
#pragma omp parallel for private(ilu, cell)
        for (ilu = start - 1; ilu >= end; ilu--)
        {
            cell = luorder[ilu];
            IntType face, c1, c2, c_tmp, count;
            RealFlow flux_s[5], flux[5], q_loc[5], DQ_loc[5], visc, tmp, vgn_tmp;
            RealGeom face_n[3], dist;

#else
    for (ilu = nTCell - 1; ilu >= 0; ilu--)
    {
        cell = luorder[ilu];
        // cell = ilu;
#endif
            for (IntType i = 0; i < 5; i++)
                flux_s[i] = 0.;
            for (IntType j = 0; j < nFPC[cell]; j++)
            {
                face = C2F[cell][j];
                count = face + face;
                c1 = f2c[count++];
                c2 = f2c[count];
                // One of c1 and c2 must be cell itself.
                if (layer[c1] < layer[cell] || layer[c2] < layer[cell])
                    // if (c1 < cell || c2 < cell)
                    continue;

                // Now its neighboring cell belongs to upper triangular
                face_n[0] = xfn[face];
                face_n[1] = yfn[face];
                face_n[2] = zfn[face];
                // if (!steady)
                //     vgn_tmp = vgn[face];
                if (c2 == cell)
                {
                    c_tmp = c1;
                    c1 = c2;
                    c2 = c_tmp;
                    face_n[0] = -face_n[0];
                    face_n[1] = -face_n[1];
                    face_n[2] = -face_n[2];
                    // if (!steady)
                    //     vgn_tmp = -vgn[face];
                }
                // assert(c1 == cell);
                for (IntType i = 0; i < 5; i++)
                {
                    q_loc[i] = q[i][c2];
                    DQ_loc[i] = DQ[i][c2];
                }
                // Calculate everything (I call it Flux) in upper triangular
                //!< Wisces: 只考虑定常，steady == 1
                FluxLUSGS3D(flux, q_loc, DQ_loc, face_n, gam, p_bar, lhs_omga);
                /**
                if (steady)
                {
                    FluxLUSGS3D(flux, q_loc, DQ_loc, face_n, gam, p_bar, lhs_omga);
                }
                else
                {
                    FluxLUSGS3D_unsteady(flux, q_loc, DQ_loc, face_n, gam, p_bar, lhs_omga, vgn_tmp);
                }
                */
                /**
                if (vis_run)
                {
                    dist = norm_dist_c2c[face];
                    visc = vis_l[c2] + vis_t[c2];
                    tmp = 2.0 * visc / (q_loc[0] * dist + TINY);
                    for (IntType i = 0; i < 5; i++)
                        flux[i] -= tmp * DQ_loc[i];
                }
                */

                // Add Flux together
                tmp = area[face];
                for (IntType i = 0; i < 5; i++)
                    flux_s[i] += tmp * flux[i];
            }
            tmp = 2.0 * Diag[cell];
            for (IntType i = 0; i < 5; i++)
                DQ[i][cell] -= flux_s[i] / tmp;

            // mflog::log.set_all_processors_out();
            //         int rank_id = 0;
            // #ifdef MF_MPICH
            //         rank_id = myid;
            // #endif
            //         if (fabs(DQ[0][cell]) > 1.0e3 * rho00)
            //         {
            //             cout << "Backward sweep: drho>1.0e3*rho00! " << DQ[0][cell] << " "
            //                  << cell << " " << rank_id << endl;
            //         }
            //         if (fabs(DQ[4][cell]) > 1.0e5 * e_stag)
            //         {
            //             cout << "Backward sweep: de>1.0e5*e_stag! " << DQ[4][cell] << " "
            //                  << cell << " " << rank_id << endl;
            //         }
            //         // if(fabs(DQ[0][cell])>1.0e4*rho00 || DQ[4][cell]>1.0e8*e_stag){
            //         if (fabs(DQ[0][cell]) > rho_max || DQ[4][cell] > e_stag_max)
            //         {
            //             std::cerr << "Error!" << endl
            //                       << "Maybe CFL too big or entropy correction coefficient too small!" << endl;
            //             // mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
            //             exit(1);
            //         }

            // limit for rho>0
            //!< Wisces: DQ_limit == 2
            /**
            if (DQ_limit == 1)
            {
                // do nothing!
            }
            else if (DQ_limit == 2)
            {
            */
            RealFlow dp, vv;
            vv = q[1][cell] * q[1][cell] + q[2][cell] * q[2][cell] + q[3][cell] * q[3][cell];
            dp = DQ[4][cell] + 0.5 * DQ[0][cell] * vv - (DQ[1][cell] * q[1][cell] + DQ[2][cell] * q[2][cell] + DQ[3][cell] * q[3][cell]);
            dp *= gamm1;
            if ((q[0][cell] + DQ[0][cell]) < rho_min || (q[0][cell] + DQ[0][cell]) > rho_max ||
                (q[4][cell] + dp) < p_min || (q[4][cell] + dp) > p_max)
            {
                DQ[0][cell] *= 0.1;
                DQ[1][cell] *= 0.1;
                DQ[2][cell] *= 0.1;
                DQ[3][cell] *= 0.1;
                DQ[4][cell] *= 0.1;
            }
            dp = DQ[4][cell] + 0.5 * DQ[0][cell] * vv - (DQ[1][cell] * q[1][cell] + DQ[2][cell] * q[2][cell] + DQ[3][cell] * q[3][cell]);
            dp *= gamm1;
            if ((q[0][cell] + DQ[0][cell]) < rho_min || (q[0][cell] + DQ[0][cell]) > rho_max ||
                (q[4][cell] + dp) < p_min || (q[4][cell] + dp) > p_max)
            {
                DQ[0][cell] *= 0.1;
                DQ[1][cell] *= 0.1;
                DQ[2][cell] *= 0.1;
                DQ[3][cell] *= 0.1;
                DQ[4][cell] *= 0.1;
            }
            dp = DQ[4][cell] + 0.5 * DQ[0][cell] * vv - (DQ[1][cell] * q[1][cell] + DQ[2][cell] * q[2][cell] + DQ[3][cell] * q[3][cell]);
            dp *= gamm1;
            if ((q[0][cell] + DQ[0][cell]) < rho_min || (q[0][cell] + DQ[0][cell]) > rho_max ||
                (q[4][cell] + dp) < p_min || (q[4][cell] + dp) > p_max)
            {
                DQ[0][cell] = 0.0;
                DQ[1][cell] = 0.0;
                DQ[2][cell] = 0.0;
                DQ[3][cell] = 0.0;
                DQ[4][cell] = 0.0;
            }
            /**
            }
            else if (DQ_limit == 3)
            {
                DQ[0][cell] = MAX(DQ[0][cell], rho_min - q[0][cell]);
                DQ[0][cell] = MIN(DQ[0][cell], rho_max - q[0][cell]);
            }
            else if (DQ_limit == 4)
            {
                RealFlow alph, alph_rho, alph_rhoe, alph_p, dp, vv, rhoe;
                vv = q[1][cell] * q[1][cell] + q[2][cell] * q[2][cell] + q[3][cell] * q[3][cell];
                rhoe = 0.5 * q[0][cell] * vv + (q[4][cell] + p_bar) / (gam - 1.0);
                dp = DQ[4][cell] + 0.5 * DQ[0][cell] * vv - (DQ[1][cell] * q[1][cell] + DQ[2][cell] * q[2][cell] + DQ[3][cell] * q[3][cell]);
                dp *= gamm1;

                alph_rho = q[0][cell] / (MAX(q[0][cell], 0.05 * rho00) + MAX(0.0, -DQ[0][cell]));
                alph_rhoe = rhoe / (MAX(rhoe, 0.05 * e_stag) + MAX(0.0, -DQ[4][cell]));
                alph_p = (q[4][cell] + p_bar) / (MAX((q[4][cell] + p_bar), 0.05 * p_bar) + MAX(0.0, -dp));
                alph = MIN(alph_rho, alph_rhoe);
                alph = MIN(alph, alph_p);
                for (IntType i = 0; i < 5; i++)
                    DQ[i][cell] *= alph;
            }
            */

#if (defined MF_OPENMP) && (defined OMP_CellColor)
        }
    }
#else
    }
#endif
// #ifdef MF_TIMING
// #ifdef MF_MPICH
//     lusgs_cal += time_tmp + MPI_Wtime();
// #endif
// #endif
}

/*!
 * @brief       Calculate the Flux in LUSGS in 3D
 * @param       flux
 * @param       q
 * @param       DQ
 * @param       fa_n
 * @param       gam
 * @param       p_bar
 * @param       lhs_omga
 * @remarks     modify according to the fun [void FluxLUSGS3D(RealFlow flux[5], RealFlow q[5], RealFlow DQ[5], RealGeom fa_n[3], RealFlow gam, RealFlow p_bar, RealFlow lhs_omga)]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-13
 */
void FluxLUSGS3D(RealFlow flux[5], RealFlow q[5], RealFlow DQ[5], RealGeom fa_n[3], RealFlow gam, RealFlow p_bar, RealFlow lhs_omga)
{
    IntType i;
    RealFlow Q[5], rv2, v_n, p, peff, c2, eig, gam1 = gam - 1.;
    RealGeom nx, ny, nz;

    nx = fa_n[0];
    ny = fa_n[1];
    nz = fa_n[2];

    Q[0] = q[0];
    Q[1] = q[0] * q[1];
    Q[2] = q[0] * q[2];
    Q[3] = q[0] * q[3];
    rv2 = 0.5 * q[0] * (q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    p = q[4];
    Q[4] = p / gam1 + rv2;

    // Normal velocity and Eigenvalues
    v_n = q[1] * nx + q[2] * ny + q[3] * nz;
    c2 = gam * (p + p_bar) / q[0];
    eig = fabs(v_n) + sqrt(c2);
    eig *= lhs_omga;

    // Need to find out the fluxes on level n
    peff = gam * p_bar / gam1;
    flux[1] = -Q[1] * v_n - p * nx;
    flux[2] = -Q[2] * v_n - p * ny;
    flux[3] = -Q[3] * v_n - p * nz;
    flux[4] = -(Q[4] + p + peff) * v_n;

    // Conservative variable on level n+1
    for (i = 0; i < 5; i++)
        Q[i] += DQ[i];
    rv2 = 0.5 * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]) / Q[0];
    p = gam1 * (Q[4] - rv2);

    // Now the flux difference due to DQ
    flux[0] = DQ[1] * nx + DQ[2] * ny + DQ[3] * nz;
    v_n *= q[0];
    v_n += flux[0];
    v_n /= Q[0];
    flux[1] += Q[1] * v_n + p * nx;
    flux[2] += Q[2] * v_n + p * ny;
    flux[3] += Q[3] * v_n + p * nz;
    flux[4] += (Q[4] + p + peff) * v_n;

    // Subtract eigenvalue terms from the flux difference
    for (i = 0; i < 5; i++)
        flux[i] -= eig * DQ[i];

    // Note: We do not check Q[0] and p here, because they are used to
    // calculate sound speed in LUSGS. Check them later in UpdateFlowField
}

/*!
 * @brief       Update flow field in 3D
 * @param       grid
 * @param       DQ
 * @remarks     modify according to the fun [void UpdateFlowField3D_CFL3d(PolyGrid *grid, RealFlow *DQ[5])]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-13
 */
void UpdateFlowField3D_CFL3d(PolyGrid *grid, RealFlow *DQ[5])
{
    IntType nTCell = grid->GetNTCell();
    IntType nBFace = grid->GetNBFace();
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

    IntType i, count[2];
    RealFlow alpq, phiq, betq;
    RealFlow /*gam,*/ gam1, rho00, p00, rhot, rhotr, ru, rv, rw, re;
    RealFlow rho_del, p_del, rho_rat, p_rat, ptmp;

    // gam = 1.4;
    rho00 = rho_s;
    p00 = p_bar;

    // grid->GetData(&gam, REAL_FLOW, 1, "gam");
    // grid->GetData(&rho00, REAL_FLOW, 1, "rho");
    // grid->GetData(&p00, REAL_FLOW, 1, "p_bar");
    gam1 = gam - 1.;

    // RealFlow rho_min, rho_max, p_min, p_max;
    // rho_min = 1.1362543e-04;
    // rho_max = 1.5799525e+01;
    // p_min = -9.3975106e+04;
    // p_max = 1.3970652e+06;
    // grid->GetData(&rho_min, REAL_FLOW, 1, "rho_min");
    // grid->GetData(&rho_max, REAL_FLOW, 1, "rho_max");
    // grid->GetData(&p_min, REAL_FLOW, 1, "p_min");
    // grid->GetData(&p_max, REAL_FLOW, 1, "p_max");

    // the const value is came from CFL3D
    alpq = -0.2;
    phiq = 1. / 0.5;
    betq = 1.0 + alpq * phiq;

    count[0] = 0;
    count[1] = 0;

#if (defined MF_OPENMP) || (defined TDTREE)
#pragma omp parallel for private(rhot, rhotr, ru, rv, rw, re, rho_del, p_del, rho_rat, p_rat, ptmp)
#endif ///<  end of MF_OPENMP or TDTREE
    for (i = 0; i < nTCell; i++)
    {
        //!< Convert q to conservative variable
        rhot = rho[i];
        ru = rhot * u[i];
        rv = rhot * v[i];
        rw = rhot * w[i];
        re = p[i] / gam1 + 0.5 * rhot * (u[i] * u[i] + v[i] * v[i] + w[i] * w[i]);

        rhot += DQ[0][i];
        ru += DQ[1][i];
        rv += DQ[2][i];
        rw += DQ[3][i];
        re += DQ[4][i];

        rhotr = 1. / (rhot + TINY);
        u[i] = ru * rhotr;
        v[i] = rv * rhotr;
        w[i] = rw * rhotr;

        rho_del = DQ[0][i];
        rho_rat = rho_del / rho[i];
        if (rho_rat < alpq)
        {
            rho_del /= betq + fabs(rho_rat) * phiq;
            // count[0]++; ///< for DEBUG
        }
        rho[i] += rho_del;
        rho[i] = MAX(rho[i], rho_min);
        rho[i] = MIN(rho[i], rho_max);

        ptmp = gam1 * (re - 0.5 * (u[i] * u[i] + v[i] * v[i] + w[i] * w[i]) * rho[i]);
        p_del = ptmp - p[i];
        p_rat = p_del / (p[i] + p00);
        if (p_rat < alpq)
        {
            p_del /= betq + fabs(p_rat) * phiq;
            // count[1]++; ///< for DEBUG
        }
        p[i] += p_del;
        p[i] = MAX(p[i], p_min);
        p[i] = MIN(p[i], p_max);
    }
    // #else
    // for (i = 0; i < nTCell; i++)
    // {
    //     // Convert q to conservative variable
    //     rhot = rho[i];
    //     ru = rhot * u[i];
    //     rv = rhot * v[i];
    //     rw = rhot * w[i];
    //     re = p[i] / gam1 + 0.5 * rhot * (u[i] * u[i] + v[i] * v[i] + w[i] * w[i]);

    //     rhot += DQ[0][i];
    //     ru += DQ[1][i];
    //     rv += DQ[2][i];
    //     rw += DQ[3][i];
    //     re += DQ[4][i];

    //     rhotr = 1. / (rhot + TINY);
    //     u[i] = ru * rhotr;
    //     v[i] = rv * rhotr;
    //     w[i] = rw * rhotr;

    //     rho_del = DQ[0][i];
    //     rho_rat = rho_del / rho[i];
    //     if (rho_rat < alpq)
    //     {
    //         rho_del /= betq + fabs(rho_rat) * phiq;
    //         count[0]++;
    //     }
    //     rho[i] += rho_del;
    //     rho[i] = MAX(rho[i], rho_min);
    //     rho[i] = MIN(rho[i], rho_max);

    //     ptmp = gam1 * (re - 0.5 * (u[i] * u[i] + v[i] * v[i] + w[i] * w[i]) * rho[i]);
    //     p_del = ptmp - p[i];
    //     p_rat = p_del / (p[i] + p00);
    //     if (p_rat < alpq)
    //     {
    //         p_del /= betq + fabs(p_rat) * phiq;
    //         count[1]++;
    //     }
    //     p[i] += p_del;
    //     p[i] = MAX(p[i], p_min);
    //     p[i] = MIN(p[i], p_max);
    // }

    // #ifdef MF_DEBUG
    // #ifdef MF_MPICH
    //     Parallel::parallel_sum(count, 2, MPI_COMM_WORLD);
    // #endif
    //     // mflog::log.set_one_processor_out();
    //     if (count[0] != 0)
    //     {
    //         cout << "Warning: " << count[0] << "Cells have been modify for rho in the UpdateFlowField3D_CFL3d!" << endl;
    //     }
    //     if (count[1] != 0)
    //     {
    //         cout << "Warning: " << count[1] << "Cells have been modify for p   in the UpdateFlowField3D_CFL3d!" << endl;
    //     }
    // #endif
}

#ifdef TDTREE

/*!
 * @brief
 * @param       grid
 * @param       Diag
 * @param       level
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-12
 */
void CalDiagLUSGS_TDTree(PolyGrid *grid, RealFlow *Diag, IntType level)
{
    IntType nTCell = grid->GetNTCell();
    IntType nBFace = grid->GetNBFace();
    IntType nTFace = grid->GetNTFace();
    IntType n = nTCell + nBFace;
    IntType *f2c = grid->Getf2c();
    RealGeom *xfn = grid->GetXfn();
    RealGeom *yfn = grid->GetYfn();
    RealGeom *zfn = grid->GetZfn();
    RealGeom *area = grid->GetFaceArea();
    RealGeom *vgn = grid->GetFaceNormalVelocity();
    RealGeom *vol = grid->GetCellVol();
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

    // IntType steady;
    // RealFlow gam, p_bar, time_accuracy, real_dt, lhs_omga;
    // grid->GetData(&steady, INT, 1, "steady");
    // grid->GetData(&gam, REAL_FLOW, 1, "gam");
    // grid->GetData(&p_bar, REAL_FLOW, 1, "p_bar");
    // grid->GetData(&time_accuracy, REAL_FLOW, 1, "time_accuracy");
    // grid->GetData(&real_dt, REAL_FLOW, 1, "real_dt");
    // grid->GetData(&lhs_omga, REAL_FLOW, 1, "lhs_omga");

    IntType i, c1, c2, count;
    RealFlow vn_1, vn_2, ss_1, ss_2;
    RealFlow eig, temp;

    // RealGeom *norm_dist_c2c = NULL;
    // norm_dist_c2c = (RealGeom *)grid->GetDataPtr(REAL_GEOM, nTFace, "norm_dist_c2c");
    // assert(norm_dist_c2c); // must exist
    // IntType vis_mode, vis_run = 0;
    // grid->GetData(&vis_mode, INT, 1, "vis_mode");
    // if (vis_mode != INVISCID)
    // {
    //     vis_run = 1;
    //     // if coarse grid doesn't want to run the viscous flux, turn it off
    //     if (level != 0)
    //     {
    //         IntType cg_vis = 1;
    //         grid->GetData(&cg_vis, INT, 1, "cg_vis");
    //         if (cg_vis == 0)
    //             vis_run = 0;
    //     }
    // }
    // RealFlow *Diag_v = NULL;
    // if (vis_run)
    // {
    //     // RealFlow dot;
    //     // RealFlow *vis_l= (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "vis_l");
    //     // RealFlow *vis_t = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "vis_t");
    //     mfmem::snew_array_1D(Diag_v, nTCell, dmrfl);
    //     assert(Diag_v != 0);
    //     for (i = 0; i < nTCell; i++)
    //         Diag_v[i] = 0.;
    // }

    // char *userArgs[17] = {(char *)grid, (char *)Diag, (char *)&level, (char *)&vis_run, (char *)Diag_v,
    //                       (char *)rho, (char *)u, (char *)v, (char *)w, (char *)p, (char *)&steady,
    //                       (char *)&gam, (char *)&p_bar, (char *)&time_accuracy, (char *)&real_dt,
    //                       (char *)&lhs_omga, (char *)norm_dist_c2c};
    char *userArgs[17] = {(char *)grid, (char *)Diag, (char *)&level, NULL, NULL,
                          (char *)rho, (char *)u, (char *)v, (char *)w, (char *)p, (char *)&steady,
                          (char *)&gam, (char *)&p_bar, (char *)&time_accuracy, (char *)&real_dt,
                          (char *)&lhs_omga, NULL};

    TDTree *TDTreeRoot = grid->GetTDTree();
    TDTreeRoot->task_traversal(CalDiagLUSGS_Kernel, NULL, userArgs, Forward);

    //     if (vis_run)
    //     {
    //         // RealFlow dot;
    //         RealFlow *vis_l = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "vis_l");
    //         RealFlow *vis_t = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "vis_t");
    // #pragma omp parallel for
    //         for (i = 0; i < nTCell; i++)
    //             Diag[i] += (vis_l[i] + vis_t[i]) * Diag_v[i] / rho[i];
    //         mfmem::sdel_array_1D(Diag_v);
    //     }
}

/*!
 * @brief
 * @param       userArgs
 * @param       treeArgs
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-12
 */
void CalDiagLUSGS_Kernel(char **userArgs, TDTreeArgs *treeArgs)
{
    // #pragma omp critical
    // {
    PolyGrid *grid = (PolyGrid *)userArgs[0];
    RealFlow *Diag = (RealFlow *)userArgs[1];
    IntType level = *(IntType *)userArgs[2];
    // IntType vis_run = *(IntType *)userArgs[3];
    // RealFlow *Diag_v = (RealFlow *)userArgs[4];
    RealFlow *rho = (RealFlow *)userArgs[5];
    RealFlow *u = (RealFlow *)userArgs[6];
    RealFlow *v = (RealFlow *)userArgs[7];
    RealFlow *w = (RealFlow *)userArgs[8];
    RealFlow *p = (RealFlow *)userArgs[9];
    IntType steady = *(IntType *)userArgs[10];
    RealFlow gam = *(RealFlow *)userArgs[11];
    RealFlow p_bar = *(RealFlow *)userArgs[12];
    RealFlow time_accuracy = *(RealFlow *)userArgs[13];
    RealFlow real_dt = *(RealFlow *)userArgs[14];
    RealFlow lhs_omga = *(RealFlow *)userArgs[15];
    // RealGeom *norm_dist_c2c = (RealGeom *)userArgs[16];

    IntType nTCell = grid->GetNTCell();
    IntType nBFace = grid->GetNBFace();
    IntType nTFace = grid->GetNTFace();
    IntType n = nTCell + nBFace;
    IntType *f2c = grid->Getf2c();
    RealGeom *xfn = grid->GetXfn();
    RealGeom *yfn = grid->GetYfn();
    RealGeom *zfn = grid->GetZfn();
    RealGeom *area = grid->GetFaceArea();
    RealGeom *vgn = grid->GetFaceNormalVelocity();
    RealGeom *vol = grid->GetCellVol();

    IntType i, c1, c2, count;
    RealFlow vn_1, vn_2, ss_1, ss_2;
    RealFlow eig, temp;

    IntType nsFace = treeArgs->firstFace;
    IntType neFace = treeArgs->lastFace + 1;
    IntType nsCell = treeArgs->firstCell;
    IntType neCell = treeArgs->lastCell + 1;

    if (neFace <= nBFace)
    {
        for (i = nsFace; i < neFace; i++)
        {
            c1 = f2c[2 * i];
            vn_1 = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
            // if (!steady)
            //     vn_1 -= vgn[i];
            vn_1 = fabs(vn_1);
            ss_1 = gam * (p[c1] + p_bar) / rho[c1];
            eig = vn_1 + sqrt(ss_1);

            Diag[c1] += 0.5 * area[i] * eig * lhs_omga; // 对角线的表达式
        }
    }
    else
    {
        for (i = nsFace; i < neFace; i++)
        {
            c1 = f2c[i + i];
            c2 = f2c[i + i + 1];

            // Cell c1
            vn_1 = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
            // if (!steady)
            //     vn_1 -= vgn[i];
            vn_1 = fabs(vn_1);
            ss_1 = gam * (p[c1] + p_bar) / rho[c1];
            eig = vn_1 + sqrt(ss_1);

            Diag[c1] += 0.5 * area[i] * eig * lhs_omga; // 对角线的表达式

            // Cell c2
            vn_2 = u[c2] * xfn[i] + v[c2] * yfn[i] + w[c2] * zfn[i];
            // if (!steady)
            //     vn_2 -= vgn[i];
            vn_2 = fabs(vn_2);
            ss_2 = gam * (p[c2] + p_bar) / rho[c2];
            eig = vn_2 + sqrt(ss_2);

            Diag[c2] += 0.5 * area[i] * eig * lhs_omga; // 对角线的表达式
        }
    }
    // if (!steady)
    // {
    //     for (i = nsCell; i < neCell; i++)
    //         Diag[i] += (1.0 + time_accuracy) * vol[i] / real_dt;
    // }
    // if (vis_run)
    // {
    //     RealFlow dot;
    //     for (i = nsFace; i < neFace; i++)
    //     {
    //         c1 = f2c[i + i];
    //         c2 = f2c[i + i + 1];
    //         dot = norm_dist_c2c[i];
    //         temp = area[i] / (dot + TINY);
    //         // Cell c1
    //         Diag_v[c1] += temp;
    //         // Cell c2
    //         if (c2 < nTCell)
    //             Diag_v[c2] += temp;
    //     }
    //     // for(i=0; i<nTCell; i++) Diag[i] += (vis_l[i] + vis_t[i])*Diag_v[i]/rho[i];
    // }
    // }
}

#ifdef GASPI
void SolveLUSGS3D_Comm_kernel(char ** userArgs, TDTreeArgs *treeArgs)
{
    PolyGrid *grid = (PolyGrid *)userArgs[0];
    RealFlow **DQ = (RealFlow **)userArgs[2];

    IntType ns = treeArgs->firstCell;
    IntType ne = treeArgs->lastCell + 1;
    if (ns >= ne) return;

    // int *boundryCell = (*treeArgs->boundryCell);
    // int boundryNbCell = (*treeArgs->boundryNbCell);
    // if (boundryNbCell)
    //     grid->GASPIComm(DQ, boundryCell, boundryNbCell);

    // pair<int,int> **C2B = grid->GetC2B();
    // IntType *nC2B = grid->GetnC2B();
    // IntType **ghost = grid->GetGhost();
    // IntType nNo , no;

    // IntType ilu, cell, i;
    // for(ilu=0;ilu<boundryNbCell;ilu++){
    //     cell = boundryCell[ilu];
    //     for(i=0;i<nC2B[cell];i++)
    //     {
    //         int neighbor = C2B[cell][i].first;
    //         no = C2B[cell][i].second;
    //         nNo = ghost[neighbor][no];
    //         grid->GASPIComm(5, DQ, neighbor, no, nNo);
    //         // grid->GASPIComm(5, DQ, neighbor, no);
    //     }
    // }

    int boundryNbCell = (*treeArgs->boundryNbCell);
    int *Dcell = (*treeArgs->boundryCellNo);
    int *n = (*treeArgs->boundryN);
    if (boundryNbCell)    
        grid->GASPIComm(DQ, n, Dcell);
}
#endif

/*!
 * @brief
 * @param       grid
 * @param       Diag
 * @param       DQ
 * @param       nFPC
 * @param       C2F
 * @param       level
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-12
 */
void SolveLUSGS3D_TDTree(PolyGrid *grid, RealFlow *Diag, RealFlow *DQ[5], IntType *nFPC, IntType **C2F, IntType level)
{
#ifdef MF_TIMING
#ifdef MF_MPICH
    double time_tmp ,lusgs_c ,lusgs_ca;
#endif
#endif    
    IntType nTCell = grid->GetNTCell();
    IntType nBFace = grid->GetNBFace();
    IntType n = nTCell + nBFace;
    IntType *f2c = grid->Getf2c();
    // Get grid metrics
    RealGeom *xfn = grid->GetXfn();
    RealGeom *yfn = grid->GetYfn();
    RealGeom *zfn = grid->GetZfn();
    RealGeom *area = grid->GetFaceArea();
    RealGeom *vgn = grid->GetFaceNormalVelocity();
    // Get flow variables
    RealFlow *q[5];
    // q[0] = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "rho");
    // q[1] = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "u");
    // q[2] = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "v");
    // q[3] = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "w");
    // q[4] = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "p");
    q[0] = (RealFlow *)grid->GetRho();
    q[1] = (RealFlow *)grid->GetU();
    q[2] = (RealFlow *)grid->GetV();
    q[3] = (RealFlow *)grid->GetW();
    q[4] = (RealFlow *)grid->GetP();

    // RealFlow gam, gamm1, p_bar, lhs_omga;
    // grid->GetData(&gam, REAL_FLOW, 1, "gam");
    // gamm1 = gam - 1.0;
    // grid->GetData(&p_bar, REAL_FLOW, 1, "p_bar");
    // grid->GetData(&lhs_omga, REAL_FLOW, 1, "lhs_omga");
    // IntType steady = 1;
    // grid->GetData(&steady, INT, 1, "steady");

    // RealFlow rho00, u00, v00, w00, e_stag;
    // grid->GetData(&rho00, REAL_FLOW, 1, "rho");
    // grid->GetData(&u00, REAL_FLOW, 1, "u");
    // grid->GetData(&v00, REAL_FLOW, 1, "v");
    // grid->GetData(&w00, REAL_FLOW, 1, "w");
    // grid->GetData(&e_stag, REAL_FLOW, 1, "e_stag");
    RealFlow rho00 = rho_s, u00 = u_s, v00 = v_s, w00 = w_s;

    // RealFlow rho_min, rho_max, p_min, p_max, e_stag_max;
    // grid->GetData(&rho_min, REAL_FLOW, 1, "rho_min");
    // grid->GetData(&rho_max, REAL_FLOW, 1, "rho_max");
    // grid->GetData(&p_min, REAL_FLOW, 1, "p_min");
    // grid->GetData(&p_max, REAL_FLOW, 1, "p_max");
    // grid->GetData(&e_stag_max, REAL_FLOW, 1, "e_stag_max");

    // IntType DQ_limit = 1;
    // grid->GetData(&DQ_limit, INT, 1, "DQ_limit");

    // IntType vis_mode, vis_run = 0;
    // grid->GetData(&vis_mode, INT, 1, "vis_mode");

    // if (vis_mode != INVISCID)
    // {
    //     vis_run = 1;
    //     // if coarse grid doesn't want to run the viscous flux, turn it off
    //     if (level != 0)
    //     {
    //         IntType cg_vis = 1;
    //         grid->GetData(&cg_vis, INT, 1, "cg_vis");
    //         if (cg_vis == 0)
    //             vis_run = 0;
    //     }
    // }
    // RealFlow *vis_l = NULL, *vis_t = NULL;
    // if (vis_run)
    // {
    //     vis_l = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "vis_l");
    //     vis_t = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "vis_t");
    // }

    // Some temporary variables
    IntType i, j, ilu, face, cell, c1, c2, c_tmp, count;
    RealFlow flux_s[5], flux[5], q_loc[5], DQ_loc[5], visc, tmp, vgn_tmp;
    RealGeom face_n[3], dist;

    IntType *luorder = grid->GetLUSGSCellOrder();           // (IntType *) grid->GetDataPtr(INT, nTCell, "LUSGSCellOrder"); // 为了对角占优，网格排序后的序号
    IntType *layer = grid->GetLUSGSLayer();                 // (IntType *) grid->GetDataPtr(INT, n, "LUSGSLayer"); // LUSGS迭代的层号，层号小为下三角，层号大为上三角
    IntType *cellsPerlayer = grid->GetLUSGScellsPerlayer(); // (IntType *)grid->GetDataPtr(INT, nTCell, "LUSGScellsPerlayer");

    IntType nTFace = grid->GetNTFace();
    // RealGeom *norm_dist_c2c = NULL;
    // norm_dist_c2c = (RealGeom *)grid->GetDataPtr(REAL_GEOM, nTFace, "norm_dist_c2c");
    // assert(norm_dist_c2c); // must exist

    // Now the Forward Sweep
    for (i = 0; i < 5; i++)
        DQ[i][luorder[0]] /= Diag[luorder[0]];

    // char *userArgs[29] = {(char *)grid, (char *)Diag, (char *)DQ, (char *)nFPC, (char *)C2F,
    //                       (char *)&level, (char *)q, (char *)&gam, (char *)&p_bar, (char *)&lhs_omga,
    //                       (char *)&steady, (char *)&rho00, (char *)&u00, (char *)&v00, (char *)&w00,
    //                       (char *)&e_stag, (char *)&rho_min, (char *)&rho_max, (char *)&p_min,
    //                       (char *)&p_max, (char *)&e_stag_max, (char *)&DQ_limit, (char *)&vis_run,
    //                       (char *)vis_l, (char *)vis_t, (char *)luorder, (char *)layer, (char *)cellsPerlayer,
    //                       (char *)norm_dist_c2c};
    char *userArgs[29] = {(char *)grid, (char *)Diag, (char *)DQ, (char *)nFPC, (char *)C2F,
                          (char *)&level, (char *)q, (char *)&gam, (char *)&p_bar, (char *)&lhs_omga,
                          (char *)&steady, (char *)&rho00, (char *)&u00, (char *)&v00, (char *)&w00,
                          (char *)&e_stag, (char *)&rho_min, (char *)&rho_max, (char *)&p_min,
                          (char *)&p_max, (char *)&e_stag_max, (char *)&DQ_limit, NULL,
                          NULL, NULL, (char *)luorder, (char *)layer, (char *)cellsPerlayer,
                          NULL};

    TDTree *TDTreeRoot = grid->GetTDTree();
#ifdef MF_TIMING
#ifdef MF_MPICH
    MPI_Barrier(MPI_COMM_WORLD);
    time_tmp = -MPI_Wtime();
#endif
#endif  
#ifdef GASPI
    // cout << "LUSGS Begin!!!!!" << endl << flush;
    TDTreeRoot->task_traversal(SolveLUSGS3D_forward_Kernel, SolveLUSGS3D_Comm_kernel, userArgs, Forward); 
#ifdef MF_TIMING
#ifdef MF_MPICH
    MPI_Barrier(MPI_COMM_WORLD);
    lusgs_ca = time_tmp + MPI_Wtime();
    lusgs_cal += lusgs_ca;
    MPI_Barrier(MPI_COMM_WORLD);
    time_tmp = -MPI_Wtime(); 
#endif
#endif    
    grid->GASPIWaitAll(5, DQ);
    // cout << "LUSGS Finish!!!!" << endl << flush;
#ifdef MF_TIMING
#ifdef MF_MPICH
    MPI_Barrier(MPI_COMM_WORLD);
    lusgs_c = time_tmp + MPI_Wtime();
    lusgs_comm += lusgs_c;
    time_tmp = lusgs_ca + lusgs_c;
    Parallel::parallel_max(&time_tmp, 1, MPI_COMM_WORLD);
    if (myid == 0) cout << "LUSGS forword sweep time = "<< lusgs_ca << endl << "GASPI wait time = "<< lusgs_c <<endl << flush;
    // Parallel::parallel_max(&lusgs_c, 1, MPI_COMM_WORLD);
    // if (myid == 0) cout << "LUSGS forward sweep + gaspi comm time = "<< lusgs_c <<endl;
#endif
#endif
    // grid->GASPIPrintf();    
#else
    TDTreeRoot->task_traversal(SolveLUSGS3D_forward_Kernel, NULL, userArgs, Forward);
#ifdef MF_TIMING
#ifdef MF_MPICH
    MPI_Barrier(MPI_COMM_WORLD);
    lusgs_ca = time_tmp + MPI_Wtime();
    lusgs_cal += lusgs_ca;
#endif
#endif
#ifdef MF_MPICH
    IntType nvar = 5;
    RealFlow *q_mpi[5];
    for (IntType j = 0; j < 5; j++)
        q_mpi[j] = DQ[j];
    MPI_Barrier(MPI_COMM_WORLD);
    time_tmp = -MPI_Wtime();        
#ifdef GASPI
    grid->GASPIRecvSendVarNeighbor_Togeth(nvar, q_mpi);
#else
    grid->RecvSendVarNeighbor_Togeth(nvar, q_mpi);
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    lusgs_c = time_tmp + MPI_Wtime();
    lusgs_comm += lusgs_c;
    time_tmp = lusgs_ca + lusgs_c;
    Parallel::parallel_max(&time_tmp, 1, MPI_COMM_WORLD);
    if (myid == 0) cout << "LUSGS forword sweep time = "<< lusgs_ca << endl << "LUSGS comm time = "<< lusgs_c <<endl << flush;
    // time_tmp = -MPI_Wtime();    
    // for (IntType i = 0; i < nvar; ++i)
    // {
    //     grid->CommInterfaceDataMPI(q_mpi[i]);
    // }
#endif
#endif
// #ifdef MF_TIMING
// #ifdef MF_MPICH
//     lusgs_comm = time_tmp + MPI_Wtime();
//     lusgs_cal += lusgs_comm;
//     Parallel::parallel_max(&lusgs_comm, 1, MPI_COMM_WORLD);
//     if (myid == 0) cout << "LUSGS time = "<< lusgs_comm <<endl;
// #endif
// #endif
    TDTreeRoot->task_traversal(SolveLUSGS3D_backward_Kernel, NULL, userArgs, Backward);
// #ifdef MF_TIMING
// #ifdef MF_MPICH
//     lusgs_comm = time_tmp + MPI_Wtime();
//     lusgs_cal += lusgs_comm;
//     Parallel::parallel_max(&lusgs_comm, 1, MPI_COMM_WORLD);
//     if (myid == 0) cout << "LUSGS time = "<< lusgs_comm <<endl;
// #endif
// #endif
}

/*!
 * @brief
 * @param       userArgs
 * @param       treeArgs
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-12
 */
void SolveLUSGS3D_forward_Kernel(char **userArgs, TDTreeArgs *treeArgs)
{
    // #pragma omp critical
    // {
    PolyGrid *grid = (PolyGrid *)userArgs[0];
    RealFlow *Diag = (RealFlow *)userArgs[1];
    RealFlow **DQ = (RealFlow **)userArgs[2];
    IntType *nFPC = (IntType *)userArgs[3];
    IntType **C2F = (IntType **)userArgs[4];
    IntType level = *(IntType *)userArgs[5];
    RealFlow **q = (RealFlow **)userArgs[6];
    RealFlow gam = *(RealFlow *)userArgs[7];
    RealFlow p_bar = *(RealFlow *)userArgs[8];
    RealFlow lhs_omga = *(RealFlow *)userArgs[9];
    IntType steady = *(IntType *)userArgs[10];
    RealFlow rho00 = *(RealFlow *)userArgs[11];
    RealFlow u00 = *(RealFlow *)userArgs[12];
    RealFlow v00 = *(RealFlow *)userArgs[13];
    RealFlow w00 = *(RealFlow *)userArgs[14];
    RealFlow e_stag = *(RealFlow *)userArgs[15];
    RealFlow rho_min = *(RealFlow *)userArgs[16];
    RealFlow rho_max = *(RealFlow *)userArgs[17];
    RealFlow p_min = *(RealFlow *)userArgs[18];
    RealFlow p_max = *(RealFlow *)userArgs[19];
    RealFlow e_stag_max = *(RealFlow *)userArgs[20];
    IntType DQ_limit = *(IntType *)userArgs[21];
    // IntType vis_run = *(IntType *)userArgs[22];
    // RealFlow *vis_l = (RealFlow *)userArgs[23];
    // RealFlow *vis_t = (RealFlow *)userArgs[24];
    IntType *luorder = (IntType *)userArgs[25];
    IntType *layer = (IntType *)userArgs[26];
    IntType *cellsPerlayer = (IntType *)userArgs[27];
    // RealGeom *norm_dist_c2c = (RealGeom *)userArgs[28];

    IntType nsCell = treeArgs->firstCell;
    IntType neCell = treeArgs->lastCell + 1;

    if (nsCell == 0)
        nsCell++;

    IntType nTCell = grid->GetNTCell();
    IntType nBFace = grid->GetNBFace();
    IntType n = nTCell + nBFace;
    IntType *f2c = grid->Getf2c();
    // Get grid metrics
    RealGeom *xfn = grid->GetXfn();
    RealGeom *yfn = grid->GetYfn();
    RealGeom *zfn = grid->GetZfn();
    RealGeom *area = grid->GetFaceArea();
    RealGeom *vgn = grid->GetFaceNormalVelocity();
    // RealFlow gamm1 = gam - 1.0;

    IntType i, j, ilu, face, cell, c1, c2, c_tmp, count;
    RealFlow flux_s[5], flux[5], q_loc[5], DQ_loc[5], visc, tmp, vgn_tmp;
    RealGeom face_n[3], dist;

    for (ilu = nsCell; ilu < neCell; ilu++)
    {
        cell = luorder[ilu];
        for (IntType j = 0; j < nFPC[cell]; j++)
        {
            face = C2F[cell][j];
            count = face + face;
            c1 = f2c[count++];
            c2 = f2c[count];
            // One of c1 and c2 must be cell itself.
            if (layer[c1] > layer[cell] || layer[c2] > layer[cell])
                continue;

            // Now its neighboring cell belongs to lower triangular
            face_n[0] = xfn[face];
            face_n[1] = yfn[face];
            face_n[2] = zfn[face];
            // if (!steady)
            //     vgn_tmp = vgn[face];
            if (c2 == cell)
            {
                c_tmp = c1;
                c1 = c2;
                c2 = c_tmp;
                face_n[0] = -face_n[0];
                face_n[1] = -face_n[1];
                face_n[2] = -face_n[2];
                // if (!steady)
                //     vgn_tmp = -vgn[face];
            }
            // assert(c1 == cell);

            for (IntType i = 0; i < 5; i++)
            {
                q_loc[i] = q[i][c2];
                DQ_loc[i] = DQ[i][c2];
            }
            // Calculate everything (I call it Flux) in lower triangular
            FluxLUSGS3D(flux, q_loc, DQ_loc, face_n, gam, p_bar, lhs_omga);
            // if (steady)
            // {
            //     FluxLUSGS3D(flux, q_loc, DQ_loc, face_n, gam, p_bar, lhs_omga);
            // }
            // else
            // {
            //     FluxLUSGS3D_unsteady(flux, q_loc, DQ_loc, face_n, gam, p_bar, lhs_omga, vgn_tmp);
            // }
            // if (vis_run)
            // {
            //     dist = norm_dist_c2c[face];
            //     visc = vis_l[c2] + vis_t[c2];
            //     tmp = 2.0 * visc / (q_loc[0] * dist + TINY);
            //     for (IntType i = 0; i < 5; i++)
            //         flux[i] -= tmp * DQ_loc[i];
            // }

            // Add Flux together
            tmp = 0.5 * area[face];
            for (IntType i = 0; i < 5; i++)
                DQ[i][cell] -= tmp * flux[i];
        }
        for (IntType i = 0; i < 5; i++)
            DQ[i][cell] /= Diag[cell];

        // limit for rho>0
        // if (DQ_limit == 1)
        // {
        //     // do nothing!
        // }
        // else if (DQ_limit == 2)
        // {
        RealFlow dp, vv;
        vv = q[1][cell] * q[1][cell] + q[2][cell] * q[2][cell] + q[3][cell] * q[3][cell];
        dp = DQ[4][cell] + 0.5 * DQ[0][cell] * vv - (DQ[1][cell] * q[1][cell] + DQ[2][cell] * q[2][cell] + DQ[3][cell] * q[3][cell]);
        dp *= gamm1;
        if ((q[0][cell] + DQ[0][cell]) < rho_min || (q[0][cell] + DQ[0][cell]) > rho_max ||
            (q[4][cell] + dp) < p_min || (q[4][cell] + dp) > p_max)
        {
            DQ[0][cell] *= 0.1;
            DQ[1][cell] *= 0.1;
            DQ[2][cell] *= 0.1;
            DQ[3][cell] *= 0.1;
            DQ[4][cell] *= 0.1;
        }
        dp = DQ[4][cell] + 0.5 * DQ[0][cell] * vv - (DQ[1][cell] * q[1][cell] + DQ[2][cell] * q[2][cell] + DQ[3][cell] * q[3][cell]);
        dp *= gamm1;
        if ((q[0][cell] + DQ[0][cell]) < rho_min || (q[0][cell] + DQ[0][cell]) > rho_max ||
            (q[4][cell] + dp) < p_min || (q[4][cell] + dp) > p_max)
        {
            DQ[0][cell] *= 0.1;
            DQ[1][cell] *= 0.1;
            DQ[2][cell] *= 0.1;
            DQ[3][cell] *= 0.1;
            DQ[4][cell] *= 0.1;
        }
        dp = DQ[4][cell] + 0.5 * DQ[0][cell] * vv - (DQ[1][cell] * q[1][cell] + DQ[2][cell] * q[2][cell] + DQ[3][cell] * q[3][cell]);
        dp *= gamm1;
        if ((q[0][cell] + DQ[0][cell]) < rho_min || (q[0][cell] + DQ[0][cell]) > rho_max ||
            (q[4][cell] + dp) < p_min || (q[4][cell] + dp) > p_max)
        {
            DQ[0][cell] = 0.0;
            DQ[1][cell] = 0.0;
            DQ[2][cell] = 0.0;
            DQ[3][cell] = 0.0;
            DQ[4][cell] = 0.0;
        }
        // }
        // else if (DQ_limit == 3)
        // {
        //     DQ[0][cell] = MAX(DQ[0][cell], rho_min - q[0][cell]);
        //     DQ[0][cell] = MIN(DQ[0][cell], rho_max - q[0][cell]);
        // }
        // else if (DQ_limit == 4)
        // {
        //     RealFlow alph, alph_rho, alph_rhoe, alph_p, dp, vv, rhoe;
        //     vv = q[1][cell] * q[1][cell] + q[2][cell] * q[2][cell] + q[3][cell] * q[3][cell];
        //     rhoe = 0.5 * q[0][cell] * vv + (q[4][cell] + p_bar) / (gam - 1.0);
        //     dp = DQ[4][cell] + 0.5 * DQ[0][cell] * vv - (DQ[1][cell] * q[1][cell] + DQ[2][cell] * q[2][cell] + DQ[3][cell] * q[3][cell]);
        //     dp *= gamm1;

        //     alph_rho = q[0][cell] / (MAX(q[0][cell], 0.05 * rho00) + MAX(0.0, -DQ[0][cell]));
        //     alph_rhoe = rhoe / (MAX(rhoe, 0.05 * e_stag) + MAX(0.0, -DQ[4][cell]));
        //     alph_p = (q[4][cell] + p_bar) / (MAX((q[4][cell] + p_bar), 0.05 * p_bar) + MAX(0.0, -dp));
        //     alph = MIN(alph_rho, alph_rhoe);
        //     alph = MIN(alph, alph_p);
        //     for (IntType i = 0; i < 5; i++)
        //         DQ[i][cell] *= alph;
        // }
    }
    // }
}

/*!
 * @brief
 * @param       userArgs
 * @param       treeArgs
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-12
 */
void SolveLUSGS3D_backward_Kernel(char **userArgs, TDTreeArgs *treeArgs)
{
    // #pragma omp critical
    // {
    PolyGrid *grid = (PolyGrid *)userArgs[0];
    RealFlow *Diag = (RealFlow *)userArgs[1];
    RealFlow **DQ = (RealFlow **)userArgs[2];
    IntType *nFPC = (IntType *)userArgs[3];
    IntType **C2F = (IntType **)userArgs[4];
    IntType level = *(IntType *)userArgs[5];
    RealFlow **q = (RealFlow **)userArgs[6];
    RealFlow gam = *(RealFlow *)userArgs[7];
    RealFlow p_bar = *(RealFlow *)userArgs[8];
    RealFlow lhs_omga = *(RealFlow *)userArgs[9];
    IntType steady = *(IntType *)userArgs[10];
    RealFlow rho00 = *(RealFlow *)userArgs[11];
    RealFlow u00 = *(RealFlow *)userArgs[12];
    RealFlow v00 = *(RealFlow *)userArgs[13];
    RealFlow w00 = *(RealFlow *)userArgs[14];
    RealFlow e_stag = *(RealFlow *)userArgs[15];
    RealFlow rho_min = *(RealFlow *)userArgs[16];
    RealFlow rho_max = *(RealFlow *)userArgs[17];
    RealFlow p_min = *(RealFlow *)userArgs[18];
    RealFlow p_max = *(RealFlow *)userArgs[19];
    RealFlow e_stag_max = *(RealFlow *)userArgs[20];
    IntType DQ_limit = *(IntType *)userArgs[21];
    // IntType vis_run = *(IntType *)userArgs[22];
    // RealFlow *vis_l = (RealFlow *)userArgs[23];
    // RealFlow *vis_t = (RealFlow *)userArgs[24];
    IntType *luorder = (IntType *)userArgs[25];
    IntType *layer = (IntType *)userArgs[26];
    IntType *cellsPerlayer = (IntType *)userArgs[27];
    // RealGeom *norm_dist_c2c = (RealGeom *)userArgs[28];

    IntType nsCell = treeArgs->firstCell;
    IntType neCell = treeArgs->lastCell + 1;

    IntType nTCell = grid->GetNTCell();
    IntType nBFace = grid->GetNBFace();
    IntType n = nTCell + nBFace;
    IntType *f2c = grid->Getf2c();
    // Get grid metrics
    RealGeom *xfn = grid->GetXfn();
    RealGeom *yfn = grid->GetYfn();
    RealGeom *zfn = grid->GetZfn();
    RealGeom *area = grid->GetFaceArea();
    RealGeom *vgn = grid->GetFaceNormalVelocity();
    // RealFlow gamm1 = gam - 1.0;

    IntType i, j, ilu, face, cell, c1, c2, c_tmp, count;
    RealFlow flux_s[5], flux[5], q_loc[5], DQ_loc[5], visc, tmp, vgn_tmp;
    RealGeom face_n[3], dist;

    // Backward Sweep
    for (ilu = neCell - 1; ilu >= nsCell; ilu--)
    {
        cell = luorder[ilu];
        for (IntType i = 0; i < 5; i++)
            flux_s[i] = 0.;
        for (IntType j = 0; j < nFPC[cell]; j++)
        {
            face = C2F[cell][j];
            count = face + face;
            c1 = f2c[count++];
            c2 = f2c[count];
            // One of c1 and c2 must be cell itself.
            if (layer[c1] < layer[cell] || layer[c2] < layer[cell])
                continue;

            // Now its neighboring cell belongs to upper triangular
            face_n[0] = xfn[face];
            face_n[1] = yfn[face];
            face_n[2] = zfn[face];
            // if (!steady)
            //     vgn_tmp = vgn[face];
            if (c2 == cell)
            {
                c_tmp = c1;
                c1 = c2;
                c2 = c_tmp;
                face_n[0] = -face_n[0];
                face_n[1] = -face_n[1];
                face_n[2] = -face_n[2];
                // if (!steady)
                //     vgn_tmp = -vgn[face];
            }
            // assert(c1 == cell);
            for (IntType i = 0; i < 5; i++)
            {
                q_loc[i] = q[i][c2];
                DQ_loc[i] = DQ[i][c2];
            }
            // Calculate everything (I call it Flux) in upper triangular
            FluxLUSGS3D(flux, q_loc, DQ_loc, face_n, gam, p_bar, lhs_omga);
            // if (steady)
            // {
            //     FluxLUSGS3D(flux, q_loc, DQ_loc, face_n, gam, p_bar, lhs_omga);
            // }
            // else
            // {
            //     FluxLUSGS3D_unsteady(flux, q_loc, DQ_loc, face_n, gam, p_bar, lhs_omga, vgn_tmp);
            // }
            // if (vis_run)
            // {
            //     dist = norm_dist_c2c[face];
            //     visc = vis_l[c2] + vis_t[c2];
            //     tmp = 2.0 * visc / (q_loc[0] * dist + TINY);
            //     for (IntType i = 0; i < 5; i++)
            //         flux[i] -= tmp * DQ_loc[i];
            // }

            // Add Flux together
            tmp = area[face];
            for (IntType i = 0; i < 5; i++)
                flux_s[i] += tmp * flux[i];
        }
        tmp = 2.0 * Diag[cell];
        for (IntType i = 0; i < 5; i++)
            DQ[i][cell] -= flux_s[i] / tmp;

        // limit for rho>0
        // if (DQ_limit == 1)
        // {
        //     // do nothing!
        // }
        // else if (DQ_limit == 2)
        // {
        RealFlow dp, vv;
        vv = q[1][cell] * q[1][cell] + q[2][cell] * q[2][cell] + q[3][cell] * q[3][cell];
        dp = DQ[4][cell] + 0.5 * DQ[0][cell] * vv - (DQ[1][cell] * q[1][cell] + DQ[2][cell] * q[2][cell] + DQ[3][cell] * q[3][cell]);
        dp *= gamm1;
        if ((q[0][cell] + DQ[0][cell]) < rho_min || (q[0][cell] + DQ[0][cell]) > rho_max ||
            (q[4][cell] + dp) < p_min || (q[4][cell] + dp) > p_max)
        {
            DQ[0][cell] *= 0.1;
            DQ[1][cell] *= 0.1;
            DQ[2][cell] *= 0.1;
            DQ[3][cell] *= 0.1;
            DQ[4][cell] *= 0.1;
        }
        dp = DQ[4][cell] + 0.5 * DQ[0][cell] * vv - (DQ[1][cell] * q[1][cell] + DQ[2][cell] * q[2][cell] + DQ[3][cell] * q[3][cell]);
        dp *= gamm1;
        if ((q[0][cell] + DQ[0][cell]) < rho_min || (q[0][cell] + DQ[0][cell]) > rho_max ||
            (q[4][cell] + dp) < p_min || (q[4][cell] + dp) > p_max)
        {
            DQ[0][cell] *= 0.1;
            DQ[1][cell] *= 0.1;
            DQ[2][cell] *= 0.1;
            DQ[3][cell] *= 0.1;
            DQ[4][cell] *= 0.1;
        }
        dp = DQ[4][cell] + 0.5 * DQ[0][cell] * vv - (DQ[1][cell] * q[1][cell] + DQ[2][cell] * q[2][cell] + DQ[3][cell] * q[3][cell]);
        dp *= gamm1;
        if ((q[0][cell] + DQ[0][cell]) < rho_min || (q[0][cell] + DQ[0][cell]) > rho_max ||
            (q[4][cell] + dp) < p_min || (q[4][cell] + dp) > p_max)
        {
            DQ[0][cell] = 0.0;
            DQ[1][cell] = 0.0;
            DQ[2][cell] = 0.0;
            DQ[3][cell] = 0.0;
            DQ[4][cell] = 0.0;
        }
        // }
        // else if (DQ_limit == 3)
        // {
        //     DQ[0][cell] = MAX(DQ[0][cell], rho_min - q[0][cell]);
        //     DQ[0][cell] = MIN(DQ[0][cell], rho_max - q[0][cell]);
        // }
        // else if (DQ_limit == 4)
        // {
        //     RealFlow alph, alph_rho, alph_rhoe, alph_p, dp, vv, rhoe;
        //     vv = q[1][cell] * q[1][cell] + q[2][cell] * q[2][cell] + q[3][cell] * q[3][cell];
        //     rhoe = 0.5 * q[0][cell] * vv + (q[4][cell] + p_bar) / (gam - 1.0);
        //     dp = DQ[4][cell] + 0.5 * DQ[0][cell] * vv - (DQ[1][cell] * q[1][cell] + DQ[2][cell] * q[2][cell] + DQ[3][cell] * q[3][cell]);
        //     dp *= gamm1;

        //     alph_rho = q[0][cell] / (MAX(q[0][cell], 0.05 * rho00) + MAX(0.0, -DQ[0][cell]));
        //     alph_rhoe = rhoe / (MAX(rhoe, 0.05 * e_stag) + MAX(0.0, -DQ[4][cell]));
        //     alph_p = (q[4][cell] + p_bar) / (MAX((q[4][cell] + p_bar), 0.05 * p_bar) + MAX(0.0, -dp));
        //     alph = MIN(alph_rho, alph_rhoe);
        //     alph = MIN(alph, alph_p);
        //     for (IntType i = 0; i < 5; i++)
        //         DQ[i][cell] *= alph;
        // }
    }
    // }
}

#endif