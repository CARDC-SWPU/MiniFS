/*!
 * @file        compute_timestep.cpp
 * @brief       计算时间步的函数
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @version     1.0
 * @date        2023-05-31
 * @copyright   Copyright (c) 2023  PERSONAL
 *
 * @par 修改日志:
 * <table>
 * <tr><th> Date        <th> Version  <th> Author  <th> Description
 * <tr><td> 2023-05-31  <td> 1.0      <td> Wisces  <td> 〔内容〕
 * <tr><td> 2023-06-15  <td> 2.0      <td> Wisces  <td> 添加头文件 para_field_global.h，采用全局变量形式存储数据
 * <tr><td> 2023-06-27  <td> 3.0      <td> Wisces  <td> 添加了 MPI 并行（使用条件编译）
 * <tr><td> 2023-07-03  <td> 4.0      <td> Wisces  <td> 添加了 OpenMP 并行（使用条件编译）
 * <tr><td> 2023-07-12  <td> 5.0      <td> Wisces  <td> 添加了 TDTree 线程级并行（使用条件编译）
 * </table>
 */

//!< C/C++ head files
#include <iostream>
// #include <cassert>
#include <cmath>   ///< sqrt()
#include <cstdlib> ///< exit()

//!< direct head file
#include "compute_timestep.h"

//!< user defined head files
#include "grid_polyhedra.h"
#include "memory_util.h"
#include "para_field_global.h"
#include "parallel_base.h"

//!< head file relying on condition-compiling
#ifdef MF_MPICH
#include <mpi.h>
#endif

#if (defined MF_OPENMP) || (defined TDTREE)
#include <omp.h>
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

using namespace std;

/*!
 * @brief       Compute time step
 * @param       grid
 * @remarks     modify according to the fun [void ComputeTimeStep(PolyGrid *grid)]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-05-31
 */
void ComputeTimeStep(PolyGrid *grid)
{
    IntType nTCell = grid->GetNTCell();
    // IntType level = grid->GetLevel();

    // IntType vis_mode;
    // grid->GetData(&vis_mode, INT, 1, "vis_mode");
    // vis_mode = INVISCID; // S_A_MODEL

    // Allocate memories for time step
    // RealFlow *dt = (RealFlow *)grid->GetDataPtr(REAL_FLOW, nTCell, "dt_timestep");
    RealFlow *dt = grid->GetDt();
    if (!dt)
    {
        snew_array_1D(dt, nTCell);
        grid->SetDt(dt);
        // assert(dt != 0);
        // grid->UpdateDataPtr(dt, REAL_FLOW, nTCell, "dt_timestep");
    }

    // If count viscous or not
    // IntType vis_run; // 0--inviscid   1--laminar   2--turbulence
    // vis_run = 0;
    //!< Wisces: 暂时只考虑无粘，vis_mode == INVISCID

    /**
    if (vis_mode == INVISCID)
    {
        vis_run = 0;
    }
    else if (vis_mode == LAMINAR)
    {
        vis_run = 1;
    }
    else
    {
        vis_run = 2;
    }
    if ((level != 0) && (vis_mode != INVISCID))
    {
        // if coarse grid doesn't want to run the viscous flux, turn it off
        IntType cg_vis = 1;
        grid->GetData(&cg_vis, INT, 1, "cg_vis");
        if (cg_vis == 0)
            vis_run = 0;
    }
    */

#ifdef TDTREE
    TimeStepNormal_new_TDTree(grid, dt, 0);
#else
    // TimeStepNormal_new(grid, dt, vis_run);
    TimeStepNormal_new(grid, dt, 0);
#endif

    LimitTimeStep(grid, dt); // note: cfl number in this function
}

/*!
 * @brief
 * @param       grid
 * @param       dt
 * @param       vis_run
 * @remarks     modify according to the fun [void TimeStepNormal_new(PolyGrid *grid, RealFlow *dt, IntType vis_run)]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-13
 */
void TimeStepNormal_new(PolyGrid *grid, RealFlow *dt, IntType vis_run)
{
    IntType nTCell = grid->GetNTCell();
    IntType nBFace = grid->GetNBFace();
    IntType nTFace = grid->GetNTFace();
    IntType n = nTCell + nBFace;
    IntType *f2c = grid->Getf2c();
    RealGeom *xfn = grid->GetXfn();
    RealGeom *yfn = grid->GetYfn();
    RealGeom *zfn = grid->GetZfn();
    RealGeom *xfc = grid->GetXfc();
    RealGeom *yfc = grid->GetYfc();
    RealGeom *zfc = grid->GetZfc();
    RealGeom *xcc = grid->GetXcc();
    RealGeom *ycc = grid->GetYcc();
    RealGeom *zcc = grid->GetZcc();
    RealGeom *vgn = grid->GetFaceNormalVelocity();
    RealGeom *area = grid->GetFaceArea();
    RealGeom *vol = grid->GetCellVol();

    RealFlow *rho = grid->GetRho();
    RealFlow *u = grid->GetU();
    RealFlow *v = grid->GetV();
    RealFlow *w = grid->GetW();
    RealFlow *p = grid->GetP();

    IntType i, c1, c2;
    RealFlow eigv, dn, vn, c2tmp, gam_tmp;

    RealFlow C = 4.0;
    RealFlow *vis_l, *vis_t;
    RealFlow prl, prt, muoopr;

//!< 暂时只考虑无粘，vis_run == 0
/**
if (vis_run)
{
    vis_l = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "vis_l");
    grid->GetData(&prl, REAL_FLOW, 1, "prl");
    vis_t = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "vis_t");
    grid->GetData(&prt, REAL_FLOW, 1, "prt");
}
*/

//!< Set dt to BIG
#ifdef MF_OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nTCell; i++)
        dt[i] = BIG;

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
#pragma omp parallel for private(i, count, c1, eigv, dn, vn, c2tmp, gam_tmp, muoopr) schedule(static, groupSize)
            for (i = startFace; i < endFace; i++)
            {
                count = 2 * i;
                c1 = f2c[count];
                c2tmp = gam * (p[c1] + p_bar) / rho[c1];
                dn = fabs((xfc[i] - xcc[c1]) * xfn[i] + (yfc[i] - ycc[c1]) * yfn[i] + (zfc[i] - zcc[c1]) * zfn[i]);

                vn = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
                // if (!steady)
                //     vn -= vgn[i];
                vn = fabs(vn);
                eigv = vn + sqrt(c2tmp);

                // if (vis_run)
                // {
                //     muoopr = vis_l[c1] / prl + vis_t[c1] / prt;
                //     gam_tmp = gam;
                //     // eigv += C*gam_tmp/rho[c1]*muoopr/(dn+TINY);
                //     eigv += C * gam_tmp / rho[c1] * muoopr * area[i] / vol[c1];
                // }
                dt[c1] = std::min(dt[c1], dn / eigv);
            }
        }

        //!< Iterface
#ifdef MF_MPICH
        count = 2 * pfacenum;
        for (i = pfacenum; i < nBFace; i++)
        {
            c1 = f2c[count];
            count += 2;
            c2tmp = gam * (p[c1] + p_bar) / rho[c1];
            dn = fabs((xfc[i] - xcc[c1]) * xfn[i] + (yfc[i] - ycc[c1]) * yfn[i] + (zfc[i] - zcc[c1]) * zfn[i]);

            vn = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
            // if (!steady)
            //     vn -= vgn[i];
            vn = fabs(vn);
            eigv = vn + sqrt(c2tmp);

            // if (vis_run)
            // {
            //     muoopr = vis_l[c1] / prl + vis_t[c1] / prt;
            //     gam_tmp = gam;
            //     // eigv += C*gam_tmp/rho[c1]*muoopr/(dn+TINY);
            //     eigv += C * gam_tmp / rho[c1] * muoopr * area[i] / vol[c1];
            // }
            dt[c1] = std::min(dt[c1], dn / eigv);
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
#pragma omp parallel for private(i, count, c1, c2, eigv, dn, vn, c2tmp, gam_tmp, muoopr) schedule(static, groupSize)
            for (i = startFace; i < endFace; i++)
            {
                count = 2 * i;
                c1 = f2c[count];
                c2 = f2c[count + 1];
                c2tmp = gam * (p[c1] + p_bar) / rho[c1];
                dn = fabs((xfc[i] - xcc[c1]) * xfn[i] + (yfc[i] - ycc[c1]) * yfn[i] + (zfc[i] - zcc[c1]) * zfn[i]);

                vn = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
                // if (!steady)
                //     vn -= vgn[i];
                vn = fabs(vn);
                eigv = vn + sqrt(c2tmp);

                // if (vis_run)
                // {
                //     muoopr = vis_l[c1] / prl + vis_t[c1] / prt;
                //     gam_tmp = gam;
                //     // eigv += C*gam_tmp/rho[c1]*muoopr/dn;
                //     eigv += C * gam_tmp / rho[c1] * muoopr * area[i] / vol[c1];
                // }
                dt[c1] = std::min(dt[c1], dn / eigv);

                c2tmp = gam * (p[c2] + p_bar) / rho[c2];
                dn = fabs((xfc[i] - xcc[c2]) * xfn[i] + (yfc[i] - ycc[c2]) * yfn[i] + (zfc[i] - zcc[c2]) * zfn[i]);

                vn = u[c2] * xfn[i] + v[c2] * yfn[i] + w[c2] * zfn[i];
                // if (!steady)
                //     vn -= vgn[i];
                vn = fabs(vn);
                eigv = vn + sqrt(c2tmp);

                // if (vis_run)
                // {
                //     muoopr = vis_l[c2] / prl + vis_t[c2] / prt;
                //     gam_tmp = gam;
                //     // eigv += C*gam_tmp/rho[c2]*muoopr/dn;
                //     eigv += C * gam_tmp / rho[c2] * muoopr * area[i] / vol[c2];
                // }
                dt[c2] = std::min(dt[c2], dn / eigv);
            }
        }
    }
    else
    {
        for (i = 0; i < nBFace; i++)
        {
            c1 = f2c[i + i];

            c2tmp = gam * (p[c1] + p_bar) / rho[c1];
            dn = fabs((xfc[i] - xcc[c1]) * xfn[i] + (yfc[i] - ycc[c1]) * yfn[i] + (zfc[i] - zcc[c1]) * zfn[i]);

            vn = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
            // if (!steady)
            //     vn -= vgn[i];
            vn = fabs(vn);
            eigv = vn + sqrt(c2tmp);

            // if (vis_run)
            // {
            //     muoopr = vis_l[c1] / prl + vis_t[c1] / prt;
            //     gam_tmp = gam;
            //     // eigv += C*gam_tmp/rho[c1]*muoopr/(dn+TINY);
            //     eigv += C * gam_tmp / rho[c1] * muoopr * area[i] / vol[c1];
            // }
            dt[c1] = std::min(dt[c1], dn / eigv);
        }

        //!< For interior faces
        for (i = nBFace; i < nTFace; i++)
        {
            c1 = f2c[i + i];
            c2 = f2c[i + i + 1];

            c2tmp = gam * (p[c1] + p_bar) / rho[c1];
            dn = fabs((xfc[i] - xcc[c1]) * xfn[i] + (yfc[i] - ycc[c1]) * yfn[i] + (zfc[i] - zcc[c1]) * zfn[i]);

            vn = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
            // if (!steady)
            //     vn -= vgn[i];
            vn = fabs(vn);
            eigv = vn + sqrt(c2tmp);

            // if (vis_run)
            // {
            //     muoopr = vis_l[c1] / prl + vis_t[c1] / prt;
            //     gam_tmp = gam;
            //     // eigv += C*gam_tmp/rho[c1]*muoopr/dn;
            //     eigv += C * gam_tmp / rho[c1] * muoopr * area[i] / vol[c1];
            // }
            dt[c1] = std::min(dt[c1], dn / eigv);

            c2tmp = gam * (p[c2] + p_bar) / rho[c2];
            dn = fabs((xfc[i] - xcc[c2]) * xfn[i] + (yfc[i] - ycc[c2]) * yfn[i] + (zfc[i] - zcc[c2]) * zfn[i]);

            vn = u[c2] * xfn[i] + v[c2] * yfn[i] + w[c2] * zfn[i];
            // if (!steady)
            //     vn -= vgn[i];
            vn = fabs(vn);
            eigv = vn + sqrt(c2tmp);

            // if (vis_run)
            // {
            //     muoopr = vis_l[c2] / prl + vis_t[c2] / prt;
            //     gam_tmp = gam;
            //     // eigv += C*gam_tmp/rho[c2]*muoopr/dn;
            //     eigv += C * gam_tmp / rho[c2] * muoopr * area[i] / vol[c2];
            // }
            dt[c2] = std::min(dt[c2], dn / eigv);
        }
    }

#elif (defined MF_OPENMP) && (defined OMP_FaceColor)
    IntType nIFace = grid->GetNIFace();
    IntType pfacenum = nBFace - nIFace;

    IntType bfacegroup_num, ifacegroup_num;
    // IntType *grid_bfacegroup, *grid_ifacegroup;
    ifacegroup_num = (*grid).ifacegroup.size(); ///< 内部面的颜色数
    bfacegroup_num = (*grid).bfacegroup.size(); ///< 物理边界面的颜色数
    // grid_bfacegroup = NULL;
    // grid_ifacegroup = NULL;
    // snew_array_1D(grid_bfacegroup, bfacegroup_num);
    // snew_array_1D(grid_ifacegroup, ifacegroup_num);
    // for (int i = 0; i < bfacegroup_num; i++)
    // {
    //     grid_bfacegroup[i] = (*grid).bfacegroup[i];
    // }
    // for (int i = 0; i < ifacegroup_num; i++)
    // {
    //     grid_ifacegroup[i] = (*grid).ifacegroup[i];
    // }

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

#pragma omp parallel for
        for (IntType i = startFace; i < endFace; i++)
        {
            IntType c1, c2;
            RealFlow eigv, dn, vn, c2tmp, gam_tmp;
            c1 = f2c[i + i];

            c2tmp = gam * (p[c1] + p_bar) / rho[c1];
            dn = fabs((xfc[i] - xcc[c1]) * xfn[i] + (yfc[i] - ycc[c1]) * yfn[i] + (zfc[i] - zcc[c1]) * zfn[i]);

            vn = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
            // if (!steady)
            //     vn -= vgn[i];
            vn = fabs(vn);
            eigv = vn + sqrt(c2tmp);

            // if (vis_run)
            // {
            //     muoopr = vis_l[c1] / prl + vis_t[c1] / prt;
            //     gam_tmp = gam;

            //     // eigv += C*gam_tmp/rho[c1]*muoopr/(dn+TINY);
            //     eigv += C * gam_tmp / rho[c1] * muoopr * area[i] / vol[c1];
            // }
            dt[c1] = std::min(dt[c1], dn / eigv);
        }
    }

    //!< Interface:
#ifdef MF_MPICH
    for (IntType i = pfacenum; i < nBFace; i++)
    {
        c1 = f2c[i + i];

        c2tmp = gam * (p[c1] + p_bar) / rho[c1];
        dn = fabs((xfc[i] - xcc[c1]) * xfn[i] + (yfc[i] - ycc[c1]) * yfn[i] + (zfc[i] - zcc[c1]) * zfn[i]);

        vn = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
        // if (!steady)
        //     vn -= vgn[i];
        vn = fabs(vn);
        eigv = vn + sqrt(c2tmp);

        // if (vis_run)
        // {
        //     muoopr = vis_l[c1] / prl + vis_t[c1] / prt;
        //     gam_tmp = gam;

        //     // eigv += C*gam_tmp/rho[c1]*muoopr/(dn+TINY);
        //     eigv += C * gam_tmp / rho[c1] * muoopr * area[i] / vol[c1];
        // }
        dt[c1] = std::min(dt[c1], dn / eigv);
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

#pragma omp parallel for
        for (IntType i = startFace; i < endFace; i++)
        {
            IntType c1, c2;
            RealFlow eigv, dn, vn, c2tmp, gam_tmp;
            c1 = f2c[i + i];
            c2 = f2c[i + i + 1];

            c2tmp = gam * (p[c1] + p_bar) / rho[c1];
            dn = fabs((xfc[i] - xcc[c1]) * xfn[i] + (yfc[i] - ycc[c1]) * yfn[i] + (zfc[i] - zcc[c1]) * zfn[i]);

            vn = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
            // if (!steady)
            //     vn -= vgn[i];
            vn = fabs(vn);
            eigv = vn + sqrt(c2tmp);

            // if (vis_run)
            // {
            //     muoopr = vis_l[c1] / prl + vis_t[c1] / prt;

            //     gam_tmp = gam;

            //     // eigv += C*gam_tmp/rho[c1]*muoopr/dn;
            //     eigv += C * gam_tmp / rho[c1] * muoopr * area[i] / vol[c1];
            // }
            dt[c1] = std::min(dt[c1], dn / eigv);

            c2tmp = gam * (p[c2] + p_bar) / rho[c2];
            dn = fabs((xfc[i] - xcc[c2]) * xfn[i] + (yfc[i] - ycc[c2]) * yfn[i] + (zfc[i] - zcc[c2]) * zfn[i]);

            vn = u[c2] * xfn[i] + v[c2] * yfn[i] + w[c2] * zfn[i];
            // if (!steady)
            //     vn -= vgn[i];
            vn = fabs(vn);
            eigv = vn + sqrt(c2tmp);

            // if (vis_run)
            // {
            //     muoopr = vis_l[c2] / prl + vis_t[c2] / prt;
            //     gam_tmp = gam;

            //     // eigv += C*gam_tmp/rho[c2]*muoopr/dn;
            //     eigv += C * gam_tmp / rho[c2] * muoopr * area[i] / vol[c2];
            // }
            dt[c2] = std::min(dt[c2], dn / eigv);
        }
    }
    // sdel_array_1D(grid_bfacegroup);
    // sdel_array_1D(grid_ifacegroup);

#elif (defined MF_OPENMP) && (defined OMP_Reduction)
    RealFlow *tmp_dt = NULL;
    IntType *nFPC = CalnFPC(grid); ///< number of faces in each real cell
    IntType **C2F = CalC2F(grid);  ///< real cell to face connection
    snew_array_1D(tmp_dt, 2 * nTFace);
    IntType j, face, count;

    //!< nBFace:
#pragma omp parallel for private(i, count, c1, eigv, dn, vn, c2tmp, gam_tmp, muoopr)
    for (i = 0; i < nBFace; i++)
    {
        count = 2 * i;
        c1 = f2c[count];

        c2tmp = gam * (p[c1] + p_bar) / rho[c1];
        dn = fabs((xfc[i] - xcc[c1]) * xfn[i] + (yfc[i] - ycc[c1]) * yfn[i] + (zfc[i] - zcc[c1]) * zfn[i]);

        vn = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
        // if (!steady)
        //     vn -= vgn[i];
        vn = fabs(vn);
        eigv = vn + sqrt(c2tmp);

        // if (vis_run)
        // {
        //     muoopr = vis_l[c1] / prl + vis_t[c1] / prt;
        //     gam_tmp = gam;

        //     // eigv += C*gam_tmp/rho[c1]*muoopr/(dn+TINY);
        //     eigv += C * gam_tmp / rho[c1] * muoopr * area[i] / vol[c1];
        // }
        tmp_dt[count] = dn / eigv;
    }

    //!< For interior faces
#pragma omp parallel for private(i, count, c1, c2, eigv, dn, vn, c2tmp, gam_tmp, muoopr)
    for (i = nBFace; i < nTFace; i++)
    {
        count = 2 * i;
        c1 = f2c[count];
        c2 = f2c[count + 1];

        c2tmp = gam * (p[c1] + p_bar) / rho[c1];
        dn = fabs((xfc[i] - xcc[c1]) * xfn[i] + (yfc[i] - ycc[c1]) * yfn[i] + (zfc[i] - zcc[c1]) * zfn[i]);

        vn = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
        // if (!steady)
        //     vn -= vgn[i];
        vn = fabs(vn);
        eigv = vn + sqrt(c2tmp);

        // if (vis_run)
        // {
        //     muoopr = vis_l[c1] / prl + vis_t[c1] / prt;

        //     gam_tmp = gam;

        //     // eigv += C*gam_tmp/rho[c1]*muoopr/dn;
        //     eigv += C * gam_tmp / rho[c1] * muoopr * area[i] / vol[c1];
        // }
        tmp_dt[count] = dn / eigv;

        c2tmp = gam * (p[c2] + p_bar) / rho[c2];
        dn = fabs((xfc[i] - xcc[c2]) * xfn[i] + (yfc[i] - ycc[c2]) * yfn[i] + (zfc[i] - zcc[c2]) * zfn[i]);

        vn = u[c2] * xfn[i] + v[c2] * yfn[i] + w[c2] * zfn[i];
        // if (!steady)
        //     vn -= vgn[i];
        vn = fabs(vn);
        eigv = vn + sqrt(c2tmp);

        // if (vis_run)
        // {
        //     muoopr = vis_l[c2] / prl + vis_t[c2] / prt;
        //     gam_tmp = gam;

        //     // eigv += C*gam_tmp/rho[c2]*muoopr/dn;
        //     eigv += C * gam_tmp / rho[c2] * muoopr * area[i] / vol[c2];
        // }
        tmp_dt[count + 1] = dn / eigv;
    }

#pragma omp parallel for private(i, j, count, c1, c2, face)
    for (i = 0; i < nTCell; i++)
    {
        for (j = 0; j < nFPC[i]; j++)
        {
            face = C2F[i][j];
            count = 2 * face;
            c1 = f2c[count];
            c2 = f2c[count + 1];
            if (i == c1)
            {
                dt[c1] = std::min(dt[c1], tmp_dt[count]);
            }
            else if (i == c2)
            {
                dt[c2] = std::min(dt[c2], tmp_dt[count + 1]);
            }
            else
            {
                // mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
                cerr << "Error in func[TimeStepNormal_new()] of file[compute_timestep.cpp]" << endl;
                exit(1);
            }
        }
    }
    sdel_array_1D(tmp_dt);

#elif (defined MF_OPENMP) && (defined OMP_DIVREP) ///< Division & replication
    IntType threads = grid->threads;
    IntType startFace, endFace, t, k, count, face;
    RealFlow tmp_dt1, tmp_dt2;
    if (grid->DivRepSuccess)
    {
#pragma omp parallel for private(t, i, k, startFace, endFace, count, c1, c2, face, eigv, dn, vn, c2tmp, gam_tmp, muoopr, tmp_dt1, tmp_dt2)
        for (t = 0; t < threads; t++)
        {
            //!< Boundary faces
            startFace = grid->idx_pthreads_bface[t];
            endFace = grid->idx_pthreads_bface[t + 1];
            for (i = startFace; i < endFace; i++)
            {
                face = grid->id_division_bface[i];
                count = 2 * face;
                c1 = f2c[count];
                c2tmp = gam * (p[c1] + p_bar) / rho[c1];
                dn = fabs((xfc[face] - xcc[c1]) * xfn[face] + (yfc[face] - ycc[c1]) * yfn[face] + (zfc[face] - zcc[c1]) * zfn[face]);

                vn = u[c1] * xfn[face] + v[c1] * yfn[face] + w[c1] * zfn[face];
                // if (!steady)
                //     vn -= vgn[face];
                vn = fabs(vn);
                eigv = vn + sqrt(c2tmp);

                // if (vis_run)
                // {
                //     muoopr = vis_l[c1] / prl + vis_t[c1] / prt;
                //     gam_tmp = gam;
                //     // eigv += C*gam_tmp/rho[c1]*muoopr/(dn+TINY);
                //     eigv += C * gam_tmp / rho[c1] * muoopr * area[face] / vol[c1];
                // }
                dt[c1] = std::min(dt[c1], dn / eigv);
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
                c2tmp = gam * (p[c1] + p_bar) / rho[c1];
                dn = fabs((xfc[face] - xcc[c1]) * xfn[face] + (yfc[face] - ycc[c1]) * yfn[face] + (zfc[face] - zcc[c1]) * zfn[face]);

                vn = u[c1] * xfn[face] + v[c1] * yfn[face] + w[c1] * zfn[face];
                // if (!steady)
                //     vn -= vgn[face];
                vn = fabs(vn);
                eigv = vn + sqrt(c2tmp);

                // if (vis_run)
                // {
                //     muoopr = vis_l[c1] / prl + vis_t[c1] / prt;
                //     gam_tmp = gam;
                //     // eigv += C*gam_tmp/rho[c1]*muoopr/dn;
                //     eigv += C * gam_tmp / rho[c1] * muoopr * area[face] / vol[c1];
                // }
                tmp_dt1 = dn / eigv;

                c2tmp = gam * (p[c2] + p_bar) / rho[c2];
                dn = fabs((xfc[face] - xcc[c2]) * xfn[face] + (yfc[face] - ycc[c2]) * yfn[face] + (zfc[face] - zcc[c2]) * zfn[face]);

                vn = u[c2] * xfn[face] + v[c2] * yfn[face] + w[c2] * zfn[face];
                // if (!steady)
                //     vn -= vgn[face];
                vn = fabs(vn);
                eigv = vn + sqrt(c2tmp);

                // if (vis_run)
                // {
                //     muoopr = vis_l[c2] / prl + vis_t[c2] / prt;
                //     gam_tmp = gam;
                //     // eigv += C*gam_tmp/rho[c2]*muoopr/dn;
                //     eigv += C * gam_tmp / rho[c2] * muoopr * area[face] / vol[c2];
                // }
                tmp_dt2 = dn / eigv;

                if (abs(k) < nTFace)
                {
                    // write back to c1 & c2
                    dt[c1] = MIN(dt[c1], tmp_dt1);
                    dt[c2] = MIN(dt[c2], tmp_dt2);
                }
                else
                {
                    if (k > 0)
                    {
                        // just write back to c1
                        dt[c1] = MIN(dt[c1], tmp_dt1);
                    }
                    else
                    {
                        // just write back to c2
                        dt[c2] = MIN(dt[c2], tmp_dt2);
                    }
                }
            }
        }
    }

#else

    //!< For boundary faces
    for (i = 0; i < nBFace; i++)
    {
        c1 = f2c[i + i];

        c2tmp = gam * (p[c1] + p_bar) / rho[c1];
        dn = fabs((xfc[i] - xcc[c1]) * xfn[i] + (yfc[i] - ycc[c1]) * yfn[i] + (zfc[i] - zcc[c1]) * zfn[i]);
        vn = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
        // if (!steady)
        //     vn -= vgn[i];
        vn = fabs(vn);
        eigv = vn + sqrt(c2tmp);

        // if (vis_run)
        // {
        //     muoopr = vis_l[c1] / prl + vis_t[c1] / prt;
        //     gam_tmp = gam;
        //     // eigv += C*gam_tmp/rho[c1]*muoopr/(dn+TINY);
        //     eigv += C * gam_tmp / rho[c1] * muoopr * area[i] / vol[c1];
        // }
        dt[c1] = std::min(dt[c1], dn / eigv);
    }
    //!< For interior faces
    for (i = nBFace; i < nTFace; i++)
    {
        c1 = f2c[i + i];
        c2 = f2c[i + i + 1];

        c2tmp = gam * (p[c1] + p_bar) / rho[c1];
        dn = fabs((xfc[i] - xcc[c1]) * xfn[i] + (yfc[i] - ycc[c1]) * yfn[i] + (zfc[i] - zcc[c1]) * zfn[i]);
        vn = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
        // if (!steady)
        //     vn -= vgn[i];
        vn = fabs(vn);
        eigv = vn + sqrt(c2tmp);

        /**
        if (vis_run)
        {
            muoopr = vis_l[c1] / prl + vis_t[c1] / prt;
            gam_tmp = gam;
            // eigv += C*gam_tmp/rho[c1]*muoopr/dn;
            eigv += C * gam_tmp / rho[c1] * muoopr * area[i] / vol[c1];
        }
        */
        dt[c1] = std::min(dt[c1], dn / eigv);

        c2tmp = gam * (p[c2] + p_bar) / rho[c2];
        dn = fabs((xfc[i] - xcc[c2]) * xfn[i] + (yfc[i] - ycc[c2]) * yfn[i] + (zfc[i] - zcc[c2]) * zfn[i]);
        vn = u[c2] * xfn[i] + v[c2] * yfn[i] + w[c2] * zfn[i];
        // if (!steady)
        //     vn -= vgn[i];
        vn = fabs(vn);
        eigv = vn + sqrt(c2tmp);

        /**
        if (vis_run)
        {
            muoopr = vis_l[c2] / prl + vis_t[c2] / prt;
            gam_tmp = gam;
            // eigv += C*gam_tmp/rho[c2]*muoopr/dn;
            eigv += C * gam_tmp / rho[c2] * muoopr * area[i] / vol[c2];
        }
        */
        dt[c2] = std::min(dt[c2], dn / eigv);
    }
#endif ///< ~MF_OPENMP
}

/*!
 * @brief       limit time step for robust
 * @param       grid
 * @param       dt
 * @remarks     modify according to the fun [void LimitTimeStep(PolyGrid *grid, RealFlow *dt)]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-13
 */
void LimitTimeStep(PolyGrid *grid, RealFlow *dt)
{
    IntType nTCell = grid->GetNTCell();
    IntType nBFace = grid->GetNBFace();
    IntType n = nTCell + nBFace;
    // IntType level = grid->GetLevel();
    RealFlow *p = grid->GetP();

    IntType i, count;
    RealFlow cfl, cfl_tmp;

    // for coarse grid, reduce cfl number using cfl_coeff
    // if (level > 0)
    // {
    //     cfl_start *= cfl_coeff;
    //     cfl_end *= cfl_coeff;
    // }

    // compute current step's cfl number
    // if (iter_done < 0)
    // { // 粗网格迭代
    //     cfl = cfl_start;
    // }
    // else if (iter_done > cfl_nstep)
    // {
    //     cfl = cfl_end;
    // }
    // else
    // {
    //     //!< modified from CFL3D, the ramping now is nonlinear, occurring slowerly at first and then increasing in rate.
    //     cfl = cfl_start * pow(cfl_ratio, (RealFlow)iter_done / cfl_nstep);
    // }
    if (iter_done > cfl_nstep)
        cfl = cfl_end;
    else
        cfl = cfl_start * pow(cfl_ratio, (RealFlow)iter_done / cfl_nstep); ///< modified from CFL3D, the ramping now is nonlinear, occurring slowerly at first and then increasing in rate.

    //!< limit cfl using gradient of p, decrease cfl in big gradient of p
    IntType *det = NULL;
    snew_array_1D(det, nTCell);
    CellIsMG(grid, det);

    cfl_min = 0.5 * cfl; ///< 在此处将 cfl_min 设为当前步 cfl 数乘以 0.5

#if (defined MF_OPENMP) || (defined TDTREE)
#pragma omp parallel for private(i, cfl_tmp)
#endif
    for (i = 0; i < nTCell; i++)
    {
        //!< 根据压力场来确定当地cfl数
        if (p[i] > p_break)
        {
            cfl_tmp = cfl;
        }
        else if (p[i] < p_min)
        {
            cfl_tmp = cfl_min;
        }
        else
        {
            cfl_tmp = (p[i] - p_min) / (p_break - p_min) * (cfl - cfl_min) + cfl_min;
        }
        //!< 根据压力梯度的极值来限制当地 cfl 数
        if (!det[i])
        {
            cfl_tmp = cfl_min;
        }

        dt[i] *= cfl_tmp;
    }

    sdel_array_1D(det);

    //!< Print out the maximum and minimun dt
    // RealFlow dt_max = 0.0, dt_min = BIG;
    // #ifdef MF_OPENMP
    // #pragma omp parallel for private(i)
    // #endif
    for (i = 0; i < nTCell; i++)
    {
        dt_max = std::max(dt_max, dt[i]);
        dt_min = std::min(dt_min, dt[i]);
    }

#ifdef MF_MPICH
    RealFlow dt_max_glb, dt_min_glb;
    MPI_Allreduce(&dt_max, &dt_max_glb, 1, MPIReal, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&dt_min, &dt_min_glb, 1, MPIReal, MPI_MIN, MPI_COMM_WORLD);
    dt_max = dt_max_glb;
    dt_min = dt_min_glb;
#endif

    //!< Now limit the dt to ratio_dtmax * dt_min
    RealFlow ratio_max = dt_max / dt_min;

    if (ratio_max > ratio_dtmax)
    {
        RealFlow dt_max_lim = ratio_dtmax * dt_min;
#ifdef MF_OPENMP
#pragma omp parallel for private(i)
#endif
        for (i = 0; i < nTCell; i++)
        {
            // dt[i] = std::min(dt[i], dt_max_lim);
            if (dt[i] > dt_max_lim)
            {
                dt[i] = dt_max_lim;
            }
        }
    }

#ifdef MF_DEBUG
    count = 0;
    if (ratio_max > ratio_dtmax)
    {
        RealFlow dt_max_lim = ratio_dtmax * dt_min;
        for (i = 0; i < nTCell; i++)
        {
            // dt[i] = std::min(dt[i], dt_max_lim);
            if (dt[i] > dt_max_lim)
            {
                count++;
            }
        }
#ifdef MF_MPICH
        Parallel::parallel_sum(count, MPI_COMM_WORLD);
        if (myid == 0)
#endif
        {
            cout << endl
                 << "dt_max/dt_min= "
                 << ratio_max << endl;
            cout << endl
                 << count << " cells are limited for dt too big." << endl;
        }
    }
#endif

#ifdef MF_DEBUG
    BCRecord **bcr = grid->Getbcr();
    IntType *f2c = grid->Getf2c();
    IntType type, c1;

    count = 0;
    RealFlow dt_avg = 0.0;
    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        if (type != WALL)
            continue;
        c1 = f2c[i + i];

        count++;
        dt_avg += dt[c1];
    }
#ifdef MF_MPICH
    Parallel::parallel_sum(dt_avg, MPI_COMM_WORLD);
    Parallel::parallel_sum(count, MPI_COMM_WORLD);
    if (myid == 0)
#endif
        cout << "First layer's average time step is: " << dt_avg / count << endl;
#endif
}

/*!
 * @brief       Determine one cell to use multigrid or not (based on cell)
 * @param       grid
 * @param       det
 * @remarks     modify according to the fun [void CellIsMG(PolyGrid *grid, IntType *det)]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-13
 */
void CellIsMG(PolyGrid *grid, IntType *det)
{
    IntType nTCell = grid->GetNTCell();
    IntType n = nTCell + grid->GetNBFace();
    IntType nTFace = grid->GetNTFace();
    IntType nBFace = grid->GetNBFace();
    IntType *f2c = grid->Getf2c();
    IntType c1, c2;

    RealFlow dp;
    RealFlow *p = grid->GetP();

    //!< 由于修改了激波探测的规则，由压力差比来流总压修改为压力差比压力和，stind参数需要调整。
    //!< 根据喷流数值试验，将 stind 参数固化为 0.1，这个值越大，判断出的激波单元越少，越小，越容易误判激波单元
    RealFlow stind = 0.1;

#if (defined MF_OPENMP) || (defined TDTREE)
#pragma omp parallel for
#endif
    for (IntType i = 0; i < nTCell; ++i)
        det[i] = 1;

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
#pragma omp parallel for private(c1, c2, dp) schedule(static, groupSize)
            for (IntType i = startFace; i < endFace; i++)
            {
                c1 = f2c[i + i];
                c2 = f2c[i + i + 1];
                dp = fabs(p[c2] - p[c1]) / (p[c2] + p[c1] + p_bar + p_bar);
                if (dp > stind)
                {
                    //!< 压力变化过大，不多重计算
                    det[c1] = 0;
                }
            }
        }

        //!< Iterface
#ifdef MF_MPICH
        for (IntType i = pfacenum; i < nBFace; i++)
        {
            c1 = f2c[i + i];
            c2 = f2c[i + i + 1];
            dp = fabs(p[c2] - p[c1]) / (p[c2] + p[c1] + p_bar + p_bar);
            if (dp > stind)
            {
                //!< 压力变化过大，不多重计算
                det[c1] = 0;
            }
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

#pragma omp parallel for private(c1, c2, dp) schedule(static, groupSize)
            for (IntType i = startFace; i < endFace; i++)
            {
                c1 = f2c[i + i];
                c2 = f2c[i + i + 1];
                dp = fabs(p[c2] - p[c1]) / (p[c2] + p[c1] + p_bar + p_bar);
                if (dp > stind)
                {
                    //!< 压力变化过大，不多重计算
                    det[c1] = 0;
                    det[c2] = 0;
                }
            }
        }
    }
    else
    {
        for (IntType i = 0; i < nBFace; ++i)
        {
            c1 = f2c[i + i];
            c2 = f2c[i + i + 1];
            dp = fabs(p[c2] - p[c1]) / (p[c2] + p[c1] + p_bar + p_bar);
            if (dp > stind)
            {
                //!< 压力变化过大，不多重计算
                det[c1] = 0;
            }
        }
        for (IntType i = nBFace; i < nTFace; ++i)
        {
            c1 = f2c[i + i];
            c2 = f2c[i + i + 1];
            dp = fabs(p[c2] - p[c1]) / (p[c2] + p[c1] + p_bar + p_bar);
            if (dp > stind)
            {
                //!< 压力变化过大，不多重计算
                det[c1] = 0;
                det[c2] = 0;
            }
        }
    }
#elif (defined MF_OPENMP) && (defined OMP_FaceColor)
    IntType nIFace = grid->GetNIFace();
    IntType pfacenum = nBFace - nIFace;

    IntType bfacegroup_num, ifacegroup_num;
    // IntType *grid_bfacegroup, *grid_ifacegroup;
    ifacegroup_num = (*grid).ifacegroup.size(); ///< 内部面的颜色数
    bfacegroup_num = (*grid).bfacegroup.size(); ///< 物理边界面的颜色数
    // grid_bfacegroup = NULL;
    // grid_ifacegroup = NULL;
    // snew_array_1D(grid_bfacegroup, bfacegroup_num);
    // snew_array_1D(grid_ifacegroup, ifacegroup_num);
    // for (int i = 0; i < bfacegroup_num; i++)
    // {
    //     grid_bfacegroup[i] = (*grid).bfacegroup[i];
    // }
    // for (int i = 0; i < ifacegroup_num; i++)
    // {
    //     grid_ifacegroup[i] = (*grid).ifacegroup[i];
    // }

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

#pragma omp parallel for private(c1, c2, dp)
        for (IntType i = startFace; i < endFace; i++)
        {
            c1 = f2c[i + i];
            c2 = f2c[i + i + 1];
            dp = fabs(p[c2] - p[c1]) / (p[c2] + p[c1] + p_bar + p_bar);
            if (dp > stind)
            {
                //!< 压力变化过大，不多重计算
                det[c1] = 0;
            }
        }
    }

    //!< Interface:
#ifdef MF_MPICH
    for (IntType i = pfacenum; i < nBFace; i++)
    {
        c1 = f2c[i + i];
        c2 = f2c[i + i + 1];
        dp = fabs(p[c2] - p[c1]) / (p[c2] + p[c1] + p_bar + p_bar);
        if (dp > stind)
        {
            //!< 压力变化过大，不多重计算
            det[c1] = 0;
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

#pragma omp parallel for private(c1, c2, dp)
        for (IntType i = startFace; i < endFace; i++)
        {
            c1 = f2c[i + i];
            c2 = f2c[i + i + 1];
            dp = fabs(p[c2] - p[c1]) / (p[c2] + p[c1] + p_bar + p_bar);
            if (dp > stind)
            {
                //!< 压力变化过大，不多重计算
                det[c1] = 0;
                det[c2] = 0;
            }
        }
    }

#elif (defined MF_OPENMP) && (defined OMP_Reduction)
    IntType *nFPC = CalnFPC(grid);
    IntType **C2F = CalC2F(grid);
    IntType face;
#pragma omp parallel for private(face, c1, c2, dp)
    for (IntType i = 0; i < nTCell; i++)
    {
        for (IntType j = 0; j < nFPC[i]; j++)
        {
            face = C2F[i][j];
            c1 = f2c[face + face];
            c2 = f2c[face + face + 1];
            dp = abs(p[c2] - p[c1]) / (p[c2] + p[c1] + p_bar + p_bar);
            if (dp > stind)
            {
                det[i] = 0;
                break;
            }
        }
    }
    // sdel_array_1D(grid_bfacegroup);
    // sdel_array_1D(grid_ifacegroup);

#elif (defined MF_OPENMP) && (defined OMP_DIVREP)
    IntType threads = grid->threads;
    IntType startFace, endFace, t, k, face, i;
    if (grid->DivRepSuccess)
    {
#pragma omp parallel for private(i, k, startFace, endFace, face, c1, c2, dp)
        for (t = 0; t < threads; t++)
        {
            //!< Boundary faces
            startFace = grid->idx_pthreads_bface[t];
            endFace = grid->idx_pthreads_bface[t + 1];
            for (i = startFace; i < endFace; i++)
            {
                face = grid->id_division_bface[i];
                c1 = f2c[face + face];
                c2 = f2c[face + face + 1];
                dp = abs(p[c2] - p[c1]) / (p[c2] + p[c1] + p_bar + p_bar);
                if (dp > stind)
                {
                    //!< 压力变化过大，不多重计算
                    det[c1] = 0;
                }
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
                c1 = f2c[face + face];
                c2 = f2c[face + face + 1];
                dp = abs(p[c2] - p[c1]) / (p[c2] + p[c1] + p_bar + p_bar);
                if (dp > stind)
                {
                    //!< 压力变化过大，不多重计算
                    if (abs(k) < nTFace)
                    {
                        //!< write back to c1 & c2
                        det[c1] = 0;
                        det[c2] = 0;
                    }
                    else
                    {
                        if (k > 0)
                        {
                            //!< just write back to c1
                            det[c1] = 0;
                        }
                        else
                        {
                            //!< just write back to c2
                            det[c2] = 0;
                        }
                    }
                }
            }
        }
    }

#else
    for (IntType i = 0; i < nBFace; ++i)
    {
        c1 = f2c[i + i];
        c2 = f2c[i + i + 1];
        dp = fabs(p[c2] - p[c1]) / (p[c2] + p[c1] + p_bar + p_bar);
        if (dp > stind)
        {
            //!< 压力变化过大，不多重计算
            det[c1] = 0;
        }
    }
    for (IntType i = nBFace; i < nTFace; ++i)
    {
        c1 = f2c[i + i];
        c2 = f2c[i + i + 1];
        dp = fabs(p[c2] - p[c1]) / (p[c2] + p[c1] + p_bar + p_bar);
        if (dp > stind)
        {
            //!< 压力变化过大，不多重计算
            det[c1] = 0;
            det[c2] = 0;
        }
    }
#endif ///< end of MF_OPENMP and OMP_*
}

#ifdef TDTREE

/*!
 * @brief
 * @param       grid
 * @param       dt
 * @param       vis_run
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-12
 */
void TimeStepNormal_new_TDTree(PolyGrid *grid, RealFlow *dt, IntType vis_run)
{
    IntType nTCell = grid->GetNTCell();
    IntType nBFace = grid->GetNBFace();
    IntType nTFace = grid->GetNTFace();
    IntType n = nTCell + nBFace;
    IntType *f2c = grid->Getf2c();
    RealGeom *xfn = grid->GetXfn();
    RealGeom *yfn = grid->GetYfn();
    RealGeom *zfn = grid->GetZfn();
    RealGeom *xfc = grid->GetXfc();
    RealGeom *yfc = grid->GetYfc();
    RealGeom *zfc = grid->GetZfc();
    RealGeom *xcc = grid->GetXcc();
    RealGeom *ycc = grid->GetYcc();
    RealGeom *zcc = grid->GetZcc();
    RealGeom *vgn = grid->GetFaceNormalVelocity();
    RealGeom *area = grid->GetFaceArea();
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

    TDTree *TDTreeRoot = grid->GetTDTree();

    // IntType steady;
    // RealFlow gam, p_bar;
    // grid->GetData(&steady, INT, 1, "steady");
    // grid->GetData(&gam, REAL_FLOW, 1, "gam");
    // grid->GetData(&p_bar, REAL_FLOW, 1, "p_bar");

    IntType i, c1, c2;
    RealFlow eigv, dn, vn, c2tmp, gam_tmp;

    RealFlow *vis_l, *vis_t;
    RealFlow prl, prt;
    // if (vis_run)
    // {
    //     vis_l = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "vis_l");
    //     grid->GetData(&prl, REAL_FLOW, 1, "prl");
    //     vis_t = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "vis_t");
    //     grid->GetData(&prt, REAL_FLOW, 1, "prt");
    // }

    // Set dt to BIG
#pragma omp parallel for
    for (i = 0; i < nTCell; i++)
        dt[i] = BIG;

    char *userArgs[28] = {(char *)grid, (char *)dt, (char *)f2c, (char *)xfn, (char *)yfn, (char *)zfn, (char *)xfc, (char *)yfc,
                          (char *)zfc, (char *)xcc, (char *)ycc, (char *)zcc, (char *)vgn, (char *)area, (char *)vol,
                          (char *)rho, (char *)u, (char *)v, (char *)w, (char *)p, (char *)&steady, (char *)&gam,
                          (char *)&p_bar, (char *)vis_l, (char *)vis_t, (char *)&prl, (char *)&prt, (char *)&vis_run};

    TDTreeRoot->task_traversal(TimeStepNormal_new_Kernel, NULL, userArgs, Forward);
}

/*!
 * @brief
 * @param       userArgs
 * @param       treeArgs
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-12
 */
void TimeStepNormal_new_Kernel(char **userArgs, TDTreeArgs *treeArgs)
{
    // #pragma omp critical
    // {
    IntType i, c1, c2, nMid;
    RealFlow eigv, dn, vn, c2tmp, gam_tmp;

    RealFlow C = 4.0;
    RealFlow muoopr;

    IntType ns = treeArgs->firstFace;
    IntType ne = treeArgs->lastFace + 1;
    // if (ns >= ne)
    //     return;

    PolyGrid *grid = (PolyGrid *)userArgs[0];
    RealFlow *dt = (RealFlow *)userArgs[1];
    IntType *f2c = (IntType *)userArgs[2];
    RealGeom *xfn = (RealGeom *)userArgs[3];
    RealGeom *yfn = (RealGeom *)userArgs[4];
    RealGeom *zfn = (RealGeom *)userArgs[5];
    RealGeom *xfc = (RealGeom *)userArgs[6];
    RealGeom *yfc = (RealGeom *)userArgs[7];
    RealGeom *zfc = (RealGeom *)userArgs[8];
    RealGeom *xcc = (RealGeom *)userArgs[9];
    RealGeom *ycc = (RealGeom *)userArgs[10];
    RealGeom *zcc = (RealGeom *)userArgs[11];
    RealGeom *vgn = (RealGeom *)userArgs[12];
    RealGeom *area = (RealGeom *)userArgs[13];
    RealGeom *vol = (RealGeom *)userArgs[14];

    RealFlow *rho = (RealFlow *)userArgs[15];
    RealFlow *u = (RealFlow *)userArgs[16];
    RealFlow *v = (RealFlow *)userArgs[17];
    RealFlow *w = (RealFlow *)userArgs[18];
    RealFlow *p = (RealFlow *)userArgs[19];

    IntType steady = (IntType)(*(IntType *)userArgs[20]);
    RealFlow gam = (RealFlow)(*(RealFlow *)userArgs[21]);
    RealFlow p_bar = (RealFlow)(*(RealFlow *)userArgs[22]);

    // RealFlow *vis_l = (RealFlow *)userArgs[23];
    // RealFlow *vis_t = (RealFlow *)userArgs[24];
    // RealFlow prl = (RealFlow)(*(RealFlow *)userArgs[25]);
    // RealFlow prt = (RealFlow)(*(RealFlow *)userArgs[26]);

    // IntType vis_run = *((IntType *)userArgs[27]);
    IntType nBFace = grid->GetNBFace();

    nMid = ns;
    if (ne <= nBFace)
    {
        // If all boundary faces
        nMid = ne;
    }
    else if (ns < nBFace)
    {
        // Part of them are boundary faces
        nMid = nBFace;
    }
    // For boundary faces
    for (i = ns; i < nMid; i++)
    {
        c1 = f2c[i + i];

        c2tmp = gam * (p[c1] + p_bar) / rho[c1];
        dn = fabs((xfc[i] - xcc[c1]) * xfn[i] + (yfc[i] - ycc[c1]) * yfn[i] + (zfc[i] - zcc[c1]) * zfn[i]);

        vn = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
        // if (!steady)
        //     vn -= vgn[i];
        vn = fabs(vn);
        eigv = vn + sqrt(c2tmp);

        // if (vis_run)
        // {
        //     muoopr = vis_l[c1] / prl + vis_t[c1] / prt;
        //     gam_tmp = gam;
        //     // eigv += C*gam_tmp/rho[c1]*muoopr/(dn+TINY);
        //     eigv += C * gam_tmp / rho[c1] * muoopr * area[i] / vol[c1];
        // }
        dt[c1] = std::min(dt[c1], dn / eigv);
    }
    // For interior faces
    for (i = nMid; i < ne; i++)
    {
        c1 = f2c[i + i];
        c2 = f2c[i + i + 1];

        c2tmp = gam * (p[c1] + p_bar) / rho[c1];
        dn = fabs((xfc[i] - xcc[c1]) * xfn[i] + (yfc[i] - ycc[c1]) * yfn[i] + (zfc[i] - zcc[c1]) * zfn[i]);

        vn = u[c1] * xfn[i] + v[c1] * yfn[i] + w[c1] * zfn[i];
        // if (!steady)
        //     vn -= vgn[i];
        vn = fabs(vn);
        eigv = vn + sqrt(c2tmp);

        // if (vis_run)
        // {
        //     muoopr = vis_l[c1] / prl + vis_t[c1] / prt;
        //     gam_tmp = gam;
        //     // eigv += C*gam_tmp/rho[c1]*muoopr/dn;
        //     eigv += C * gam_tmp / rho[c1] * muoopr * area[i] / vol[c1];
        // }
        dt[c1] = std::min(dt[c1], dn / eigv);

        c2tmp = gam * (p[c2] + p_bar) / rho[c2];
        dn = fabs((xfc[i] - xcc[c2]) * xfn[i] + (yfc[i] - ycc[c2]) * yfn[i] + (zfc[i] - zcc[c2]) * zfn[i]);

        vn = u[c2] * xfn[i] + v[c2] * yfn[i] + w[c2] * zfn[i];
        // if (!steady)
        //     vn -= vgn[i];
        vn = fabs(vn);
        eigv = vn + sqrt(c2tmp);

        // if (vis_run)
        // {
        //     muoopr = vis_l[c2] / prl + vis_t[c2] / prt;
        //     gam_tmp = gam;
        //     // eigv += C*gam_tmp/rho[c2]*muoopr/dn;
        //     eigv += C * gam_tmp / rho[c2] * muoopr * area[i] / vol[c2];
        // }
        dt[c2] = std::min(dt[c2], dn / eigv);
    }
    // }
}

#endif