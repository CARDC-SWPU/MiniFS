/*!
 * @file        cal_gradient.cpp
 * @brief
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @version     1.0
 * @date        2023-07-07
 * @copyright   Copyright (c) 2023  PERSONAL
 *
 * @par 修改日志:
 * <table>
 * <tr><th> Date        <th> Version  <th> Author  <th> Description
 * <tr><td> 2023-07-07  <td> 1.0      <td> Wisces  <td> 〔内容〕
 * <tr><td> 2023-07-12  <td> 2.0      <td> Wisces  <td> 添加了 TDTree 线程级并行（使用条件编译）
 * </table>
 */

//!< direct head file
#include "cal_gradient.h"

//!< C/C++ head files
#include <iostream>
#include <cstdlib>    ///< exit()
#include <cmath>      ///< sqrt()
#include <sys/time.h> ///< gettimeofday()

//!< user defined head files
#include "grid_polyhedra.h"
#include "memory_util.h"
#include "para_field_global.h"

//!< head file relying on condition-compiling
#ifdef MF_OPENMP
#include <omp.h>
#endif

#ifdef TDTREE
#include "TDTree.h"
#endif

using namespace std;

/*!
 * @brief       为 rho，u，v，w，p 梯度分配内存，并计算第一次的初始梯度值
 * @param       grid
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-07
 */
void AllocAndCalQuantityGradient(PolyGrid *grid)
{
    //!< 分配内存
    const IntType nTCell = grid->GetNTCell();
    const IntType n = nTCell + grid->GetNBFace();

    RealFlow *dqdx = NULL, *dqdy = NULL, *dqdz = NULL;

    snew_array_1D(dqdx, 5 * n);
    snew_array_1D(dqdy, 5 * n);
    snew_array_1D(dqdz, 5 * n);

#ifdef MF_OPENMP
#pragma omp parallel for
#endif
    for (IntType i = 0; i < 5 * n; ++i)
    {
        dqdx[i] = 0;
        dqdy[i] = 0;
        dqdz[i] = 0;
    }

    grid->SetDqdx(dqdx);
    grid->SetDqdy(dqdy);
    grid->SetDqdz(dqdz);

    //!< 调用梯度重构函数计算梯度值
    CalculateGradient(grid);
}

/*!
 * @brief       Zero dqdx, dqdy and dqdz. Drive the actual gradient calculation functions in 3D.
 * @param       grid
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-08
 */
void CalculateGradient(PolyGrid *grid /*, RealFlow *dqdx[5], RealFlow *dqdy[5], RealFlow *dqdz[5]*/)
{
    const IntType nTCell = grid->GetNTCell();
    const IntType n = nTCell + grid->GetNBFace();
    const IntType nTNode = grid->GetNTNode();

    // RealFlow *_dqdx[5], *_dqdy[5], *_dqdz[5];
    // _dqdx[0] = static_cast<RealFlow *>(dqdx);
    // _dqdy[0] = static_cast<RealFlow *>(dqdy);
    // _dqdz[0] = static_cast<RealFlow *>(dqdz);
    RealFlow *dqdx[5], *dqdy[5], *dqdz[5];
    dqdx[0] = static_cast<RealFlow *>(grid->GetDqdx());
    dqdy[0] = static_cast<RealFlow *>(grid->GetDqdy());
    dqdz[0] = static_cast<RealFlow *>(grid->GetDqdz());
    for (IntType i = 1; i < 5; ++i)
    {
        dqdx[i] = &dqdx[i - 1][n];
        dqdy[i] = &dqdy[i - 1][n];
        dqdz[i] = &dqdz[i - 1][n];
    }

    //!< Initialize dq
#ifdef MF_OPENMP
#pragma omp parallel for
#endif
    for (IntType i = 0; i < 5 * n; i++)
    {
        dqdx[0][i] = 0.;
        dqdy[0][i] = 0.;
        dqdz[0][i] = 0.;
    }

    RealFlow *q[5];
    q[0] = static_cast<RealFlow *>(grid->GetRho());
    q[1] = static_cast<RealFlow *>(grid->GetU());
    q[2] = static_cast<RealFlow *>(grid->GetV());
    q[3] = static_cast<RealFlow *>(grid->GetW());
    q[4] = static_cast<RealFlow *>(grid->GetP());

#ifdef TDTREE

    // CompGradientQ(grid, q, dqdx, dqdy, dqdz, 0, kNVar);
    CompGradientQ_Gauss_Node_TDTree(grid, q, dqdx, dqdy, dqdz, 0, 5);

#else

    RealFlow *u_n = NULL, *v_n = NULL, *w_n = NULL;
    snew_array_1D(u_n, nTNode);
    snew_array_1D(v_n, nTNode);
    snew_array_1D(w_n, nTNode);

    for (IntType i = 0; i < 5; ++i)
    {
        // CompGradientQ(grid, q[i], _dqdx[i], _dqdy[i], _dqdz[i], i);
        CompGradientQ_Gauss_Node(grid, q[i], dqdx[i], dqdy[i], dqdz[i], i, u_n, v_n, w_n);
    }

//!< 并行传值
#ifdef MF_MPICH
    RealFlow **grad_mpi = NULL;
    snew_array_1D(grad_mpi, 3 * 5);
    IntType count = 0;
    for (IntType i = 0; i < 5; ++i)
    {
        grad_mpi[count++] = dqdx[i];
        grad_mpi[count++] = dqdy[i];
        grad_mpi[count++] = dqdz[i];
    }
    grid->RecvSendVarNeighbor_Togeth(3 * 5, grad_mpi);
    // for (IntType i = 0; i < 3 * 5; ++i)
    // {
    //     grid->CommInterfaceDataMPI(grad_mpi[i]);
    // }
    sdel_array_1D(grad_mpi);
#endif

    sdel_array_1D(u_n);
    sdel_array_1D(v_n);
    sdel_array_1D(w_n);
#endif
}

/*!
 * @brief       Calculate the gradients of flow variable q in 3D.
 * @param       grid
 * @param       q
 * @param       dqdx
 * @param       dqdy
 * @param       dqdz
 * @param       name
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-08
 */
// void CompGradientQ(PolyGrid *grid, RealFlow *q, RealFlow *dqdx, RealFlow *dqdy, RealFlow *dqdz, IntType name)
// {
//     CompGradientQ_Gauss_Node(grid, q, dqdx, dqdy, dqdz, name);
// }

/*!
 * @brief       Calculate the gradients of flow variable q(rho, u, v, w, p) in 3D use Node-Green-Gauss Approach.
 * @param       grid
 * @param       q
 * @param       dqdx
 * @param       dqdy
 * @param       dqdz
 * @param       name
 * @remarks     modify according to the fun [void CompGradientQ_Gauss_Node(PolyGrid *grid, RealFlow *q, RealFlow *dqdx, RealFlow *dqdy, RealFlow *dqdz, IntType name)]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-08
 */
void CompGradientQ_Gauss_Node(PolyGrid *grid, RealFlow *q, RealFlow *dqdx, RealFlow *dqdy, RealFlow *dqdz, IntType name, RealFlow *u_n, RealFlow *v_n, RealFlow *w_n)
{
    IntType nTNode = grid->GetNTNode();
    IntType nTCell = grid->GetNTCell();
    IntType nTFace = grid->GetNTFace();
    IntType nBFace = grid->GetNBFace();
    IntType *f2c = grid->Getf2c();
    BCRecord **bcr = grid->Getbcr();
    IntType n = nTCell + nBFace;
    RealGeom *vol = grid->GetCellVol();
    RealGeom *area = grid->GetFaceArea();
    RealGeom *xfn = grid->GetXfn();
    RealGeom *yfn = grid->GetYfn();
    RealGeom *zfn = grid->GetZfn();
    IntType *nNPF = grid->GetnNPF();
    IntType **F2N = CalF2N(grid);
    IntType *nFPC = CalnFPC(grid);
    IntType **C2F = CalC2F(grid);

    IntType i, j, c1, c2, count, type, face;
    RealGeom tmpx, tmpy, tmpz;
    RealFlow qsum;

    //!< Initialize dq
    // #ifdef MF_OPENMP
    // #pragma omp parallel for
    // #endif
    // for (i = 0; i < n; i++)
    // {
    //     dqdx[i] = 0.;
    //     dqdy[i] = 0.;
    //     dqdz[i] = 0.;
    // }

    RealFlow *q_n = NULL;
    snew_array_1D(q_n, nTNode);
    CompNodeVar3D_dist(grid, q_n, q, name, u_n, v_n, w_n);

#ifdef MF_TIMING
    double time_tmp = 0.0;
#ifdef MF_MPICH
    MPI_Barrier(MPI_COMM_WORLD);
    time_tmp = -MPI_Wtime();
#else
    struct timeval starttimeTemRoe, endtimeTemRoe;
    gettimeofday(&starttimeTemRoe, 0);
#endif
#endif

#if (defined MF_OPENMP) && (defined OMP_GroupColor)
    IntType pfacenum = nBFace - grid->GetNIFace();
    if (grid->GroupColorSuccess)
    {
        IntType groupSize = grid->groupSize;
        IntType bfacegroup_num, ifacegroup_num;
        IntType startFace, endFace;
        bfacegroup_num = grid->bfacegroup.size();
        ifacegroup_num = grid->ifacegroup.size();

        // Boundary faces:
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

#pragma omp parallel for private(i, j, count, c1, c2, type, qsum) schedule(static, groupSize)
            for (i = startFace; i < endFace; i++)
            {
                count = 2 * i;
                c1 = f2c[count];
                c2 = f2c[count + 1];
                type = bcr[i]->GetType();
                qsum = 0.0;
                if (type == INTERFACE || type == SYMM)
                {
                    for (j = 0; j < nNPF[i]; j++)
                        qsum += q_n[F2N[i][j]];
                    qsum /= RealFlow(nNPF[i]);
                }
                else
                {
                    qsum = 0.5 * (q[c1] + q[c2]);
                }

                qsum *= area[i];
                dqdx[c1] += qsum * xfn[i];
                dqdy[c1] += qsum * yfn[i];
                dqdz[c1] += qsum * zfn[i];
            }
        }

        //!< Interfece
        count = 2 * pfacenum;
        for (i = pfacenum; i < nBFace; i++)
        {
            c1 = f2c[count++];
            c2 = f2c[count++];
            type = bcr[i]->GetType();
            qsum = 0.0;

            if (type == INTERFACE || type == SYMM)
            {
                for (j = 0; j < nNPF[i]; j++)
                    qsum += q_n[F2N[i][j]];
                qsum /= RealFlow(nNPF[i]);
            }
            else
            {
                qsum = 0.5 * (q[c1] + q[c2]);
            }

            qsum *= area[i];
            // tmpx = qsum * xfn[i];
            // tmpy = qsum * yfn[i];
            // tmpz = qsum * zfn[i];
            // dqdx[c1] += tmpx;
            // dqdy[c1] += tmpy;
            // dqdz[c1] += tmpz;
            dqdx[c1] += qsum * xfn[i];
            dqdy[c1] += qsum * yfn[i];
            dqdz[c1] += qsum * zfn[i];
        }

        // Interior faces
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

#pragma omp parallel for private(i, j, count, c1, c2, qsum, tmpx, tmpy, tmpz) schedule(static, groupSize)
            for (i = startFace; i < endFace; i++)
            {
                count = 2 * i;
                c1 = f2c[count];
                c2 = f2c[count + 1];
                qsum = 0.0;

                for (j = 0; j < nNPF[i]; j++)
                    qsum += q_n[F2N[i][j]];
                qsum /= RealFlow(nNPF[i]);

                qsum *= area[i];
                tmpx = qsum * xfn[i];
                tmpy = qsum * yfn[i];
                tmpz = qsum * zfn[i];

                // For cell c1
                dqdx[c1] += tmpx;
                dqdy[c1] += tmpy;
                dqdz[c1] += tmpz;

                // For cell c2
                dqdx[c2] -= tmpx;
                dqdy[c2] -= tmpy;
                dqdz[c2] -= tmpz;
            }
        }
    }
    else
    {
        for (i = 0; i < nBFace; i++)
        {
            count = 2 * i;
            c1 = f2c[count];
            c2 = f2c[count + 1];
            type = bcr[i]->GetType();
            qsum = 0.0;

            if (type == INTERFACE || type == SYMM)
            {
                for (j = 0; j < nNPF[i]; j++)
                    qsum += q_n[F2N[i][j]];
                qsum /= RealFlow(nNPF[i]);
            }
            else
            {
                qsum = 0.5 * (q[c1] + q[c2]);
            }

            qsum *= area[i];
            tmpx = qsum * xfn[i];
            tmpy = qsum * yfn[i];
            tmpz = qsum * zfn[i];
            dqdx[c1] += tmpx;
            dqdy[c1] += tmpy;
            dqdz[c1] += tmpz;
        }

        for (i = nBFace; i < nTFace; i++)
        {
            count = 2 * i;
            c1 = f2c[count];
            c2 = f2c[count + 1];
            qsum = 0.0;

            for (j = 0; j < nNPF[i]; j++)
                qsum += q_n[F2N[i][j]];
            qsum /= RealFlow(nNPF[i]);

            qsum *= area[i];
            tmpx = qsum * xfn[i];
            tmpy = qsum * yfn[i];
            tmpz = qsum * zfn[i];
            // For cell c1
            dqdx[c1] += tmpx;
            dqdy[c1] += tmpy;
            dqdz[c1] += tmpz;

            // For cell c2
            dqdx[c2] -= tmpx;
            dqdy[c2] -= tmpy;
            dqdz[c2] -= tmpz;
        }
    }
#elif (defined MF_OPENMP) && (defined OMP_FaceColor)
    IntType pfacenum = nBFace - grid->GetNIFace();
    IntType bfacegroup_num, ifacegroup_num;
    IntType startFace, endFace;
    ifacegroup_num = (*grid).ifacegroup.size();
    bfacegroup_num = (*grid).bfacegroup.size();
    // Boundary faces:
    for (IntType fcolor = 0; fcolor < bfacegroup_num; fcolor++)
    {
        if (fcolor == 0)
        {
            startFace = 0;
        }
        else
        {
            startFace = (*grid).bfacegroup[fcolor - 1];
        }
        endFace = (*grid).bfacegroup[fcolor];

#pragma omp parallel for private(i, j, c1, c2, count, type, qsum)
        for (i = startFace; i < endFace; i++)
        {
            // IntType j, c1, c2, count, type;
            // RealGeom tmpx, tmpy, tmpz;
            // RealFlow qsum;
            count = 2 * i;
            c1 = f2c[count];
            c2 = f2c[count + 1];
            type = bcr[i]->GetType();
            qsum = 0.0;

            if (type == INTERFACE || type == SYMM)
            {
                for (j = 0; j < nNPF[i]; j++)
                    qsum += q_n[F2N[i][j]];
                qsum /= RealFlow(nNPF[i]);
            }
            else
            {
                qsum = 0.5 * (q[c1] + q[c2]);
            }

            qsum *= area[i];
            // tmpx = qsum * xfn[i];
            // tmpy = qsum * yfn[i];
            // tmpz = qsum * zfn[i];
            dqdx[c1] += qsum * xfn[i];
            dqdy[c1] += qsum * yfn[i];
            dqdz[c1] += qsum * zfn[i];
        }
    }

    //!< Interface
    for (i = pfacenum; i < nBFace; i++)
    {
        count = 2 * i;
        c1 = f2c[count];
        c2 = f2c[count + 1];
        type = bcr[i]->GetType();
        qsum = 0.0;
        // RealGeom tmpx, tmpy, tmpz;

        if (type == INTERFACE || type == SYMM)
        {
            for (j = 0; j < nNPF[i]; j++)
                qsum += q_n[F2N[i][j]];
            qsum /= RealFlow(nNPF[i]);
        }
        else
        {
            qsum = 0.5 * (q[c1] + q[c2]);
        }

        qsum *= area[i];
        // tmpx = qsum * xfn[i];
        // tmpy = qsum * yfn[i];
        // tmpz = qsum * zfn[i];
        // dqdx[c1] += tmpx;
        // dqdy[c1] += tmpy;
        // dqdz[c1] += tmpz;
        dqdx[c1] += qsum * xfn[i];
        dqdy[c1] += qsum * yfn[i];
        dqdz[c1] += qsum * zfn[i];
    }

    // Interior faces:
    for (IntType fcolor = 0; fcolor < ifacegroup_num; fcolor++)
    {
        // IntType startFace, endFace;
        if (fcolor == 0)
        {
            startFace = nBFace;
        }
        else
        {
            startFace = (*grid).ifacegroup[fcolor - 1];
        }
        endFace = (*grid).ifacegroup[fcolor];

#pragma omp parallel for private(i, j, c1, c2, count, tmpx, tmpy, tmpz, qsum)
        for (i = startFace; i < endFace; i++)
        {
            // IntType j, c1, c2, count;
            // RealGeom tmpx, tmpy, tmpz;
            // RealFlow qsum;
            count = 2 * i;
            c1 = f2c[count];
            c2 = f2c[count + 1];
            qsum = 0.0;

            for (j = 0; j < nNPF[i]; j++)
                qsum += q_n[F2N[i][j]];
            qsum /= RealFlow(nNPF[i]);

            qsum *= area[i];
            tmpx = qsum * xfn[i];
            tmpy = qsum * yfn[i];
            tmpz = qsum * zfn[i];
            // For cell c1
            dqdx[c1] += tmpx;
            dqdy[c1] += tmpy;
            dqdz[c1] += tmpz;
            // For cell c2
            dqdx[c2] -= tmpx;
            dqdy[c2] -= tmpy;
            dqdz[c2] -= tmpz;
        }
    }

#elif (defined MF_OPENMP) && (defined OMP_Reduction)
    RealGeom *tmpxyz = NULL;
    snew_array_1D(tmpxyz, 3 * nTFace);
#pragma omp parallel for private(i, j, count, c1, c2, type, qsum)
    for (i = 0; i < nBFace; i++)
    {
        count = 2 * i;
        c1 = f2c[count];
        c2 = f2c[count + 1];
        type = bcr[i]->GetType();
        qsum = 0.0;

        if (type == INTERFACE || type == SYMM)
        {
            for (j = 0; j < nNPF[i]; j++)
                qsum += q_n[F2N[i][j]];
            qsum /= RealFlow(nNPF[i]);
        }
        else
        {
            qsum = 0.5 * (q[c1] + q[c2]);
        }
        j = 3 * i;
        qsum *= area[i];
        tmpxyz[j] = qsum * xfn[i];
        tmpxyz[j + 1] = qsum * yfn[i];
        tmpxyz[j + 2] = qsum * zfn[i];
    }
#pragma omp parallel for private(i, j, count, c1, c2, qsum)
    for (i = nBFace; i < nTFace; i++)
    {
        count = 2 * i;
        c1 = f2c[count];
        c2 = f2c[count + 1];
        qsum = 0.0;

        for (j = 0; j < nNPF[i]; j++)
            qsum += q_n[F2N[i][j]];
        qsum /= RealFlow(nNPF[i]);
        j = 3 * i;
        qsum *= area[i];
        tmpxyz[j] = qsum * xfn[i];
        tmpxyz[j + 1] = qsum * yfn[i];
        tmpxyz[j + 2] = qsum * zfn[i];
    }

#pragma omp parallel for private(i, j, count, c1, c2, face)
    for (i = 0; i < nTCell; i++)
    {
        for (j = 0; j < nFPC[i]; j++)
        {
            face = C2F[i][j];
            count = 3 * face;
            c1 = f2c[face + face];
            c2 = f2c[face + face + 1];
            if (i == c1)
            {
                dqdx[i] += tmpxyz[count];
                dqdy[i] += tmpxyz[count + 1];
                dqdz[i] += tmpxyz[count + 2];
            }
            else if (i == c2)
            {
                dqdx[i] -= tmpxyz[count];
                dqdy[i] -= tmpxyz[count + 1];
                dqdz[i] -= tmpxyz[count + 2];
            }
            else
            {
                // mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
                cerr << "Error in func[CompGradientQ_Gauss_Node()] of file[load_grid.cpp]!" << endl;
                exit(1);
            }
        }
    }
    sdel_array_1D(tmpxyz);

#elif (defined MF_OPENMP) && (defined OMP_DIVREP) // Division & replication
    IntType threads = grid->threads;
    IntType startFace, endFace, t, k;
    if (grid->DivRepSuccess)
    {
#pragma omp parallel for private(t, i, j, k, startFace, endFace, count, c1, c2, face, type, qsum, tmpx, tmpy, tmpz)
        for (t = 0; t < threads; t++)
        {
            // Boundary faces
            startFace = grid->idx_pthreads_bface[t];
            endFace = grid->idx_pthreads_bface[t + 1];
            for (i = startFace; i < endFace; i++)
            {
                face = grid->id_division_bface[i];
                count = 2 * face;
                c1 = f2c[count];
                c2 = f2c[count + 1];
                type = bcr[face]->GetType();
                qsum = 0.0;
                if (type == INTERFACE || type == SYMM)
                {
                    for (j = 0; j < nNPF[face]; j++)
                        qsum += q_n[F2N[face][j]];
                    qsum /= RealFlow(nNPF[face]);
                }
                else
                {
                    qsum = 0.5 * (q[c1] + q[c2]);
                }
                qsum *= area[face];
                tmpx = qsum * xfn[face];
                tmpy = qsum * yfn[face];
                tmpz = qsum * zfn[face];
                dqdx[c1] += tmpx;
                dqdy[c1] += tmpy;
                dqdz[c1] += tmpz;
            }
            // Interior faces
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
                qsum = 0.0;
                for (j = 0; j < nNPF[face]; j++)
                    qsum += q_n[F2N[face][j]];
                qsum /= RealFlow(nNPF[face]);
                qsum *= area[face];
                tmpx = qsum * xfn[face];
                tmpy = qsum * yfn[face];
                tmpz = qsum * zfn[face];
                if (abs(k) < nTFace)
                {
                    dqdx[c1] += tmpx;
                    dqdy[c1] += tmpy;
                    dqdz[c1] += tmpz;
                    dqdx[c2] -= tmpx;
                    dqdy[c2] -= tmpy;
                    dqdz[c2] -= tmpz;
                }
                else
                {
                    if (k > 0)
                    {
                        dqdx[c1] += tmpx;
                        dqdy[c1] += tmpy;
                        dqdz[c1] += tmpz;
                    }
                    else
                    {
                        dqdx[c2] -= tmpx;
                        dqdy[c2] -= tmpy;
                        dqdz[c2] -= tmpz;
                    }
                }
            }
        }
    }

    // #elif (defined MF_OPENMP) && (defined DIVCON) // D&C TREE
    //     RealGeom *tmpxyz = NULL;
    //     snew_array_1D(tmpxyz, 3 * (nTFace - nBFace) d);
    // #pragma omp parallel
    //     {
    // #pragma omp single nowait
    //         tree_traversal(grid->treeHead, dqdx, dqdy, dqdz, tmpxyz,
    //                        f2c, bcr, nNPF, F2N, q_n, q, area, xfn, yfn, zfn, nBFace);
    //     }
    //     sdel_array_1D(tmpxyz);

    // #elif (defined FS_SIMD) && (defined FaceColoring) && (!defined FS_SIMD_AVX) && (!defined Tile) // add by ruitian, 2021.11.30
    //     IntType bfacegroup_num, ifacegroup_num;
    //     IntType *grid_bfacegroup, *grid_ifacegroup;
    //     ifacegroup_num = (*grid).ifacegroup.size();
    //     bfacegroup_num = (*grid).bfacegroup.size();
    //     grid_bfacegroup = NULL;
    //     grid_ifacegroup = NULL;
    //     const IntType Vec = 4;
    //     snew_array_1D(grid_bfacegroup, bfacegroup_numd);
    //     snew_array_1D(grid_ifacegroup, ifacegroup_numd);
    //     for (int i = 0; i < bfacegroup_num; i++)
    //     {
    //         grid_bfacegroup[i] = (*grid).bfacegroup[i];
    //     }
    //     for (int i = 0; i < ifacegroup_num; i++)
    //     {
    //         grid_ifacegroup[i] = (*grid).ifacegroup[i];
    //     }
    //     // #ifdef _OPENMP
    //      // #pragma omp parallel for private(j,c1,c2,count,type,qsum,tmpx,tmpy,tmpz) \
//       reduction(+:dqdx[0:n],dqdy[0:n],dqdz[0:n])// schedule(static)
    //     for (IntType fcolor = 0; fcolor < bfacegroup_num; fcolor++)
    //     {
    //         IntType startFace, endFace;
    //         if (fcolor == 0)
    //         {
    //             startFace = 0;
    //         }
    //         else
    //         {
    //             startFace = grid_bfacegroup[fcolor - 1];
    //         }
    //         endFace = grid_bfacegroup[fcolor];
    //         // #pragma omp parallel for private(type,j,c1,c2,count,qsum,tmpx,tmpy,tmpz)
    //         for (IntType i = startFace; i < endFace; i++)
    //         {
    //             count = 2 * i;
    //             c1 = f2c[count];
    //             c2 = f2c[count + 1];
    //             type = bcr[i]->GetType();
    //             qsum = 0.0;
    //             if (type == INTERFACE || type == SYMM)
    //             {
    //                 for (j = 0; j < nNPF[i]; j++)
    //                     qsum += q_n[F2N[i][j]];
    //                 qsum /= RealFlow(nNPF[i]);
    //             }
    //             else
    //             {
    //                 qsum = 0.5 * (q[c1] + q[c2]);
    //             }
    //             qsum *= area[i];
    //             tmpx = qsum * xfn[i];
    //             tmpy = qsum * yfn[i];
    //             tmpz = qsum * zfn[i];
    //             dqdx[c1] += tmpx;
    //             dqdy[c1] += tmpy;
    //             dqdz[c1] += tmpz;
    //         }
    //     }
    //     for (IntType i = pfacenum; i < nBFace; i++)
    //     {
    //         IntType count = 2 * i;
    //         IntType c1 = f2c[count];
    //         IntType c2 = f2c[count + 1];
    //         IntType type = bcr[i]->GetType();
    //         RealFlow qsum = 0.0;
    //         RealGeom tmpx, tmpy, tmpz;
    //         if (type == INTERFACE || type == SYMM)
    //         {
    //             for (j = 0; j < nNPF[i]; j++)
    //                 qsum += q_n[F2N[i][j]];
    //             qsum /= RealFlow(nNPF[i]);
    //         }
    //         else
    //         {
    //             qsum = 0.5 * (q[c1] + q[c2]);
    //         }
    //         qsum *= area[i];
    //         tmpx = qsum * xfn[i];
    //         tmpy = qsum * yfn[i];
    //         tmpz = qsum * zfn[i];
    //         dqdx[c1] += tmpx;
    //         dqdy[c1] += tmpy;
    //         dqdz[c1] += tmpz;
    //     }
    //     // Interior faces
    //      // #pragma omp parallel for private(j,c1,c2,count,qsum,tmpx,tmpy,tmpz) \
//       reduction(+:dqdx[0:n],dqdy[0:n],dqdz[0:n])// schedule(static)
    //     for (IntType fcolor = 0; fcolor < ifacegroup_num; fcolor++)
    //     {
    //         IntType startFace, endFace;
    //         if (fcolor == 0)
    //         {
    //             startFace = nBFace;
    //         }
    //         else
    //         {
    //             startFace = grid_ifacegroup[fcolor - 1];
    //         }
    //         endFace = grid_ifacegroup[fcolor];
    //         // #pragma omp parallel for private(j,c1,c2,count,qsum,tmpx,tmpy,tmpz)
    //         // #pragma omp parallel for
    //         IntType k;
    //         for (k = startFace; k + Vec < endFace; k += Vec)
    //         {
    //             IntType jv[Vec], c1v[Vec], c2v[Vec], countv[Vec];
    //             RealGeom xfnv[Vec], yfnv[Vec], zfnv[Vec], areav[Vec];
    //             RealGeom qsumv[Vec];
    //             RealGeom dqdxc1v[Vec], dqdxc2v[Vec];
    //             RealGeom dqdyc1v[Vec], dqdyc2v[Vec];
    //             RealGeom dqdzc1v[Vec], dqdzc2v[Vec];
    //             RealGeom tmpxv[Vec], tmpyv[Vec], tmpzv[Vec];
    //             RealGeom nNPFv[Vec];
    //             // Load:
    //             // #pragma omp simd
    //             for (IntType iv = 0; iv < Vec; iv++)
    //             {
    //                 countv[iv] = 2 * (iv + k);
    //                 c1v[iv] = f2c[countv[iv]];
    //                 c2v[iv] = f2c[countv[iv] + 1];
    //                 qsumv[iv] = 0.0;
    //                 areav[iv] = area[iv + k];
    //                 xfnv[iv] = xfn[iv + k];
    //                 yfnv[iv] = yfn[iv + k];
    //                 zfnv[iv] = zfn[iv + k];
    //                 nNPFv[iv] = nNPF[iv + k];
    //                 // c1 cell:
    //                 dqdxc1v[iv] = dqdx[c1v[iv]];
    //                 dqdyc1v[iv] = dqdy[c1v[iv]];
    //                 dqdzc1v[iv] = dqdz[c1v[iv]];
    //                 // c2 cell:
    //                 dqdxc2v[iv] = dqdx[c2v[iv]];
    //                 dqdyc2v[iv] = dqdy[c2v[iv]];
    //                 dqdzc2v[iv] = dqdz[c2v[iv]];
    //             }
    //             // Computation:
    //             for (IntType iv = 0; iv < Vec; iv++)
    //             {
    //                 for (j = 0; j < nNPFv[iv]; j++)
    //                 {
    //                     qsumv[iv] += q_n[F2N[iv + k][j]];
    //                 }
    //             }
    // #pragma omp simd safelen(Vec)
    //             for (IntType iv = 0; iv < Vec; iv++)
    //             {
    //                 qsumv[iv] /= nNPFv[iv];
    //                 qsumv[iv] *= areav[iv];
    //                 tmpxv[iv] = qsumv[iv] * xfnv[iv];
    //                 tmpyv[iv] = qsumv[iv] * yfnv[iv];
    //                 tmpzv[iv] = qsumv[iv] * zfnv[iv];
    //                 // c1 cell:
    //                 dqdxc1v[iv] += tmpxv[iv];
    //                 dqdyc1v[iv] += tmpyv[iv];
    //                 dqdzc1v[iv] += tmpzv[iv];
    //                 // c2 cell:
    //                 dqdxc2v[iv] -= tmpxv[iv];
    //                 dqdyc2v[iv] -= tmpyv[iv];
    //                 dqdzc2v[iv] -= tmpzv[iv];
    //             }
    //             // Load Back:
    //             for (IntType iv = 0; iv < Vec; iv++)
    //             {
    //                 // c1 cell:
    //                 dqdx[c1v[iv]] = dqdxc1v[iv];
    //                 dqdy[c1v[iv]] = dqdyc1v[iv];
    //                 dqdz[c1v[iv]] = dqdzc1v[iv];
    //                 // c2 cell:
    //                 dqdx[c2v[iv]] = dqdxc2v[iv];
    //                 dqdy[c2v[iv]] = dqdyc2v[iv];
    //                 dqdz[c2v[iv]] = dqdzc2v[iv];
    //             }
    //         }
    //         for (IntType i = k; i < endFace; i++)
    //         {
    //             IntType j, c1, c2, count;
    //             RealGeom tmpx, tmpy, tmpz;
    //             RealFlow qsum;
    //             count = 2 * i;
    //             c1 = f2c[count];
    //             c2 = f2c[count + 1];
    //             qsum = 0.0;
    //             for (j = 0; j < nNPF[i]; j++)
    //                 qsum += q_n[F2N[i][j]];
    //             qsum /= RealFlow(nNPF[i]);
    //             qsum *= area[i];
    //             tmpx = qsum * xfn[i];
    //             tmpy = qsum * yfn[i];
    //             tmpz = qsum * zfn[i];
    //             // For cell c1
    //             dqdx[c1] += tmpx;
    //             dqdy[c1] += tmpy;
    //             dqdz[c1] += tmpz;
    //             // For cell c2
    //             dqdx[c2] -= tmpx;
    //             dqdy[c2] -= tmpy;
    //             dqdz[c2] -= tmpz;
    //         }
    //     }
    //     sdel_array_1D(grid_bfacegroup);
    //     sdel_array_1D(grid_ifacegroup);

    // #elif (defined FS_SIMD) && (defined Tile)                                  // add by ruitian, 2021.11.30, for AVX512 simd based on tile
    //     for (i = 0; i < nBFace; i++)
    //     {
    //         count = 2 * i;
    //         c1 = f2c[count];
    //         c2 = f2c[count + 1];
    //         type = bcr[i]->GetType();
    //         qsum = 0.0;
    //         if (type == INTERFACE || type == SYMM)
    //         {
    //             for (j = 0; j < nNPF[i]; j++)
    //                 qsum += q_n[F2N[i][j]];
    //             qsum /= RealFlow(nNPF[i]);
    //         }
    //         else
    //         {
    //             qsum = 0.5 * (q[c1] + q[c2]);
    //         }
    //         qsum *= area[i];
    //         tmpx = qsum * xfn[i];
    //         tmpy = qsum * yfn[i];
    //         tmpz = qsum * zfn[i];
    //         dqdx[c1] += tmpx;
    //         dqdy[c1] += tmpy;
    //         dqdz[c1] += tmpz;
    //     }
    //     // added by ruitianSIMD
    //     RealFlow *qsumiface = NULL;
    //     snew_array_1D(qsumiface, nTFaced);
    //     // transform the grid information into the tile:
    //     // ruitian, 2021.12.21
    //     // cout<<"start tile comp."<<endl;
    //     for (IntType i = 0; i < nTFace; i++)
    //     {
    //         qsumiface[i] = 0.;
    //         grid->qsumtile[i] = 0.;
    //     }
    //     for (IntType i = 0; i < nTCell; i++)
    //     {
    //         grid->dqdxtile[i] = dqdx[i];
    //         grid->dqdytile[i] = dqdy[i];
    //         grid->dqdztile[i] = dqdz[i];
    //     }
    //     RealFlow *qsumtile = NULL;
    //     snew_array_1D(qsumtile, nTFaced);
    //     for (IntType i = 0; i < nTFace; i++)
    //     {
    //         qsumtile[i] = 0.;
    //     }
    //     for (IntType i = nBFace; i < nTFace; i++)
    //     {
    //         for (IntType j = 0; j < nNPF[i]; j++)
    //             qsumtile[i] += q_n[F2N[i][j]];
    //         qsumtile[i] /= RealFlow(nNPF[i]);
    //     }
    //     for (IntType ii = 0; ii < grid->iSIMDnnz; ii++)
    //     {
    //         IntType i = grid->iSIMDval[ii];
    //         grid->qsumt[ii] = qsumtile[i];
    //     }
    //     for (IntType ii = 0; ii < (*grid->ioffsets).size(); ii++)
    //     {
    //         // #pragma omp parallel for num_threads(2)
    //         for (IntType jj = 0; jj < (*grid->ioffsets)[ii].size() - 1; jj++)
    //         {
    //             // this for cycle contains an independent tile
    //             IntType kk;
    //             __m256i vc1, vc2;
    //             __m512i vc1tem, vc2tem;
    //             __m512d vtmpx, vtmpy, vtmpz;
    //             __m512d vqsum, varea; // , vfacezero;
    //             __m512d vxfn, vyfn, vzfn;
    //             __m512d vdqdxc1, vdqdyc1, vdqdzc1;
    //             __m512d vdqdxc2, vdqdyc2, vdqdzc2;
    //             // for (kk = (*grid->ioffsets)[ii][jj]; kk< (*grid->ioffsets)[ii][jj + 1]; kk ++) {
    //             for (kk = (*grid->ioffsets)[ii][jj]; kk < (*grid->ioffsets)[ii][jj + 1]; kk += 8)
    //             {
    //                 ///*
    //                 // this for cycle contains 8 times faces
    //                 if (grid->ifacezero[kk + 7] == 1)
    //                 {
    //                     // load c1 and c2:
    //                     vc1tem = _mm512_load_epi64(&grid->iSIMDrow[kk]);
    //                     vc2tem = _mm512_load_epi64(&grid->iSIMDcol[kk]);
    //                     vc1 = _mm512_castsi512_si256(vc1tem);
    //                     vc2 = _mm512_castsi512_si256(vc2tem);
    //                     // load data on face:
    //                     vxfn = _mm512_load_pd(&grid->xfnt[kk]);
    //                     vyfn = _mm512_load_pd(&grid->yfnt[kk]);
    //                     vzfn = _mm512_load_pd(&grid->zfnt[kk]);
    //                     vqsum = _mm512_load_pd(&grid->qsumt[kk]);
    //                     varea = _mm512_load_pd(&grid->areat[kk]);
    //                     // comput.
    //                     vqsum = _mm512_mul_pd(vqsum, varea);
    //                     vtmpx = _mm512_mul_pd(vqsum, vxfn);
    //                     vtmpy = _mm512_mul_pd(vqsum, vyfn);
    //                     vtmpz = _mm512_mul_pd(vqsum, vzfn);
    //                     // gather c1:
    //                     vdqdxc1 = _mm512_i32gather_pd(vc1, grid->dqdxtile, 8);
    //                     vdqdyc1 = _mm512_i32gather_pd(vc1, grid->dqdytile, 8);
    //                     vdqdzc1 = _mm512_i32gather_pd(vc1, grid->dqdztile, 8);
    //                     // comput.
    //                     vdqdxc1 = _mm512_add_pd(vdqdxc1, vtmpx);
    //                     vdqdyc1 = _mm512_add_pd(vdqdyc1, vtmpy);
    //                     vdqdzc1 = _mm512_add_pd(vdqdzc1, vtmpz);
    //                     // scatter c1:
    //                     _mm512_i32scatter_pd(grid->dqdxtile, vc1, vdqdxc1, 8);
    //                     _mm512_i32scatter_pd(grid->dqdytile, vc1, vdqdyc1, 8);
    //                     _mm512_i32scatter_pd(grid->dqdztile, vc1, vdqdzc1, 8);
    //                     // gather c2:
    //                     vdqdxc2 = _mm512_i32gather_pd(vc2, grid->dqdxtile, 8);
    //                     vdqdyc2 = _mm512_i32gather_pd(vc2, grid->dqdytile, 8);
    //                     vdqdzc2 = _mm512_i32gather_pd(vc2, grid->dqdztile, 8);
    //                     vdqdxc2 = _mm512_sub_pd(vdqdxc2, vtmpx);
    //                     vdqdyc2 = _mm512_sub_pd(vdqdyc2, vtmpy);
    //                     vdqdzc2 = _mm512_sub_pd(vdqdzc2, vtmpz);
    //                     // scatter c2:
    //                     _mm512_i32scatter_pd(grid->dqdxtile, vc2, vdqdxc2, 8);
    //                     _mm512_i32scatter_pd(grid->dqdytile, vc2, vdqdyc2, 8);
    //                     _mm512_i32scatter_pd(grid->dqdztile, vc2, vdqdzc2, 8);
    //                 }
    //                 else
    //                 {
    //                     for (IntType ikk = kk; ikk < kk + 8; ikk++)
    //                     {
    //                         if (grid->ifacezero[kk] == 0)
    //                         {
    //                             break;
    //                         }
    //                         else
    //                         {
    //                             IntType c1 = grid->iSIMDrow[ikk];
    //                             IntType c2 = grid->iSIMDcol[ikk];
    //                             grid->qsumt[ikk] *= grid->areat[ikk];
    //                             RealGeom tmpx = grid->qsumt[ikk] * grid->xfnt[ikk] * grid->ifacezero[ikk];
    //                             RealGeom tmpy = grid->qsumt[ikk] * grid->yfnt[ikk] * grid->ifacezero[ikk];
    //                             RealGeom tmpz = grid->qsumt[ikk] * grid->zfnt[ikk] * grid->ifacezero[ikk];
    //                             // For cell c1
    //                             grid->dqdxtile[c1] += tmpx;
    //                             grid->dqdytile[c1] += tmpy;
    //                             grid->dqdztile[c1] += tmpz;
    //                             // For cell c2
    //                             grid->dqdxtile[c2] -= tmpx;
    //                             grid->dqdytile[c2] -= tmpy;
    //                             grid->dqdztile[c2] -= tmpz;
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    //     // transform back:
    //     for (IntType i = 0; i < nTCell; i++)
    //     {
    //         dqdx[i] = grid->dqdxtile[i];
    //         dqdy[i] = grid->dqdytile[i];
    //         dqdz[i] = grid->dqdztile[i];
    //     }
    //     sdel_array_1D(qsumtile);

    // #elif (defined FS_SIMD) && (defined FS_SIMD_AVX) && (defined FaceColoring) // for AVX512 simd based on face coloring
    //     IntType bfacegroup_num, ifacegroup_num;
    //     IntType *grid_bfacegroup, *grid_ifacegroup;
    //     ifacegroup_num = (*grid).ifacegroup.size();
    //     bfacegroup_num = (*grid).bfacegroup.size();
    //     grid_bfacegroup = NULL;
    //     grid_ifacegroup = NULL;
    //     snew_array_1D(grid_bfacegroup, bfacegroup_numd);
    //     snew_array_1D(grid_ifacegroup, ifacegroup_numd);
    //     for (int i = 0; i < bfacegroup_num; i++)
    //     {
    //         grid_bfacegroup[i] = (*grid).bfacegroup[i];
    //     }
    //     for (int i = 0; i < ifacegroup_num; i++)
    //     {
    //         grid_ifacegroup[i] = (*grid).ifacegroup[i];
    //     }
    //     // add by ruitian, for SIMD AVX facecoloring
    //     // #pragma omp parallel for
    //     for (IntType i = 0; i < nTCell; i++)
    //     {
    //         grid->dqdxtile[i] = dqdx[i];
    //         grid->dqdytile[i] = dqdy[i];
    //         grid->dqdztile[i] = dqdz[i];
    //     }
    //     RealFlow *bqsumtile = NULL;
    //     // snew_array_1D(qsumtile, nTFaced);
    //     bqsumtile = (RealGeom *)_mm_malloc(sizeof(RealGeom) * nBFace, 64);
    //     // #pragma omp parallel for
    //     for (IntType i = 0; i < nBFace; i++)
    //     {
    //         bqsumtile[i] = 0.;
    //     }
    //     // #pragma omp parallel for
    //     for (IntType i = 0; i < nBFace; i++)
    //     {
    //         IntType c1, c2, type;
    //         c1 = grid->f2c1[i];
    //         c2 = grid->f2c2[i];
    //         type = bcr[i]->GetType();
    //         if (type == INTERFACE || type == SYMM)
    //         {
    //             for (j = 0; j < nNPF[i]; j++)
    //                 bqsumtile[i] += q_n[F2N[i][j]];
    //             bqsumtile[i] /= RealFlow(nNPF[i]);
    //         }
    //         else
    //         {
    //             bqsumtile[i] = 0.5 * (q[c1] + q[c2]);
    //         }
    //         bqsumtile[i] *= area[i];
    //     }
    //     // Boundary faces:
    //     for (IntType fcolor = 0; fcolor < bfacegroup_num; fcolor++)
    //     {
    //         IntType startFace, endFace;
    //         if (fcolor == 0)
    //         {
    //             startFace = 0;
    //         }
    //         else
    //         {
    //             startFace = grid_bfacegroup[fcolor - 1];
    //         }
    //         endFace = grid_bfacegroup[fcolor];
    //         IntType i;
    // #ifdef MF_OPENMP
    // #pragma omp parallel for
    // #endif
    //         for (i = startFace; i < endFace - 8; i += 8)
    //         {
    //             __m256i vc1;
    //             __m512i vc1tem;
    //             __m512d vtmpx, vtmpy, vtmpz;
    //             __m512d vqsum;
    //             __m512d vxfn, vyfn, vzfn;
    //             __m512d vdqdxc1, vdqdyc1, vdqdzc1;
    //             // load c1 and c2:
    //             //_mm256_load_epi32
    //             // vc1 = _mm256_load_epi32(&f2c[i]);
    //             vc1tem = _mm512_load_epi64(&grid->f2c1[i]);
    //             vc1 = _mm512_castsi512_si256(vc1tem);
    //             // load data on face:
    //             vxfn = _mm512_load_pd(&grid->xfntile[i]);
    //             vyfn = _mm512_load_pd(&grid->yfntile[i]);
    //             vzfn = _mm512_load_pd(&grid->zfntile[i]);
    //             vqsum = _mm512_load_pd(&bqsumtile[i]);
    //             // gather:
    //             vdqdxc1 = _mm512_i32gather_pd(vc1, grid->dqdxtile, 8);
    //             vdqdyc1 = _mm512_i32gather_pd(vc1, grid->dqdytile, 8);
    //             vdqdzc1 = _mm512_i32gather_pd(vc1, grid->dqdztile, 8);
    //             // comput.
    //             vtmpx = _mm512_mul_pd(vqsum, vxfn);
    //             vtmpy = _mm512_mul_pd(vqsum, vyfn);
    //             vtmpz = _mm512_mul_pd(vqsum, vzfn);
    //             vdqdxc1 = _mm512_add_pd(vdqdxc1, vtmpx);
    //             vdqdyc1 = _mm512_add_pd(vdqdyc1, vtmpy);
    //             vdqdzc1 = _mm512_add_pd(vdqdzc1, vtmpz);
    //             // scatter:
    //             _mm512_i32scatter_pd(grid->dqdxtile, vc1, vdqdxc1, 8);
    //             _mm512_i32scatter_pd(grid->dqdytile, vc1, vdqdyc1, 8);
    //             _mm512_i32scatter_pd(grid->dqdztile, vc1, vdqdzc1, 8);
    //         }
    //         for (IntType vi = i; vi < endFace; vi++)
    //         {
    //             IntType c1, c2;
    //             RealGeom tmpx, tmpy, tmpz;
    //             c1 = grid->f2c1[vi];
    //             c2 = grid->f2c2[vi];
    //             tmpx = bqsumtile[vi] * xfn[vi];
    //             tmpy = bqsumtile[vi] * yfn[vi];
    //             tmpz = bqsumtile[vi] * zfn[vi];
    //             grid->dqdxtile[c1] += tmpx;
    //             grid->dqdytile[c1] += tmpy;
    //             grid->dqdztile[c1] += tmpz;
    //         }
    //     }
    //     for (IntType i = pfacenum; i < nBFace; i++)
    //     {
    //         IntType count = 2 * i;
    //         IntType c1 = f2c[count];
    //         IntType c2 = f2c[count + 1];
    //         IntType type = bcr[i]->GetType();
    //         RealFlow qsum = 0.0;
    //         RealGeom tmpx, tmpy, tmpz;
    //         if (type == INTERFACE || type == SYMM)
    //         {
    //             for (j = 0; j < nNPF[i]; j++)
    //                 qsum += q_n[F2N[i][j]];
    //             qsum /= RealFlow(nNPF[i]);
    //         }
    //         else
    //         {
    //             qsum = 0.5 * (q[c1] + q[c2]);
    //         }
    //         qsum *= area[i];
    //         tmpx = qsum * xfn[i];
    //         tmpy = qsum * yfn[i];
    //         tmpz = qsum * zfn[i];
    //         grid->dqdxtile[c1] += tmpx;
    //         grid->dqdytile[c1] += tmpy;
    //         grid->dqdztile[c1] += tmpz;
    //     }
    //     /*************************************************************************/
    //     // ruitianSIMD, 2021.12.28
    //     // added by ruitian, for SIMD
    //     // transform the grid information into the tile:
    //     // ruitian, 2021.12.21
    //     RealFlow *qsumtile = NULL;
    //     // snew_array_1D(qsumtile, nTFaced);
    //     qsumtile = (RealGeom *)_mm_malloc(sizeof(RealGeom) * nTFace, 64);
    // #ifdef MF_OPENMP
    // #pragma omp parallel for
    // #endif
    //     for (IntType i = 0; i < nTFace; i++)
    //     {
    //         qsumtile[i] = 0.;
    //     }
    // #ifdef MF_OPENMP
    // #pragma omp parallel for
    // #endif
    //     for (IntType i = nBFace; i < nTFace; i++)
    //     {
    //         for (j = 0; j < nNPF[i]; j++)
    //             qsumtile[i] += q_n[F2N[i][j]];
    //         qsumtile[i] /= RealFlow(nNPF[i]);
    //         qsumtile[i] *= area[i];
    //     }
    //     // Interior faces:
    //     for (IntType fcolor = 0; fcolor < ifacegroup_num; fcolor++)
    //     {
    //         IntType startFace, endFace;
    //         if (fcolor == 0)
    //         {
    //             startFace = nBFace;
    //         }
    //         else
    //         {
    //             startFace = grid_ifacegroup[fcolor - 1];
    //         }
    //         endFace = grid_ifacegroup[fcolor];
    //         IntType i;
    // #ifdef MF_OPENMP
    // #pragma omp parallel for
    // #endif
    //         for (i = startFace; i < endFace - 8; i += 8)
    //         {
    //             __m256i vc1, vc2;
    //             __m512i vc1tem, vc2tem;
    //             __m512d vtmpx, vtmpy, vtmpz;
    //             __m512d vqsum;
    //             __m512d vxfn, vyfn, vzfn;
    //             __m512d vdqdxc1, vdqdyc1, vdqdzc1;
    //             __m512d vdqdxc2, vdqdyc2, vdqdzc2;
    //             // load c1 and c2:
    //             //_mm256_load_epi32
    //             // vc1 = _mm256_load_epi32(&f2c[i]);
    //             vc1tem = _mm512_load_epi64(&grid->f2c1[i]);
    //             vc2tem = _mm512_load_epi64(&grid->f2c2[i]);
    //             vc1 = _mm512_castsi512_si256(vc1tem);
    //             vc2 = _mm512_castsi512_si256(vc2tem);
    //             // load data on face:
    //             vxfn = _mm512_load_pd(&grid->xfntile[i]);
    //             vyfn = _mm512_load_pd(&grid->yfntile[i]);
    //             vzfn = _mm512_load_pd(&grid->zfntile[i]);
    //             vqsum = _mm512_load_pd(&qsumtile[i]);
    //             // gather:
    //             vdqdxc1 = _mm512_i32gather_pd(vc1, grid->dqdxtile, 8);
    //             vdqdyc1 = _mm512_i32gather_pd(vc1, grid->dqdytile, 8);
    //             vdqdzc1 = _mm512_i32gather_pd(vc1, grid->dqdztile, 8);
    //             vdqdxc2 = _mm512_i32gather_pd(vc2, grid->dqdxtile, 8);
    //             vdqdyc2 = _mm512_i32gather_pd(vc2, grid->dqdytile, 8);
    //             vdqdzc2 = _mm512_i32gather_pd(vc2, grid->dqdztile, 8);
    //             // comput.
    //             vtmpx = _mm512_mul_pd(vqsum, vxfn);
    //             vtmpy = _mm512_mul_pd(vqsum, vyfn);
    //             vtmpz = _mm512_mul_pd(vqsum, vzfn);
    //             vdqdxc1 = _mm512_add_pd(vdqdxc1, vtmpx);
    //             vdqdyc1 = _mm512_add_pd(vdqdyc1, vtmpy);
    //             vdqdzc1 = _mm512_add_pd(vdqdzc1, vtmpz);
    //             vdqdxc2 = _mm512_sub_pd(vdqdxc2, vtmpx);
    //             vdqdyc2 = _mm512_sub_pd(vdqdyc2, vtmpy);
    //             vdqdzc2 = _mm512_sub_pd(vdqdzc2, vtmpz);
    //             // scatter:
    //             _mm512_i32scatter_pd(grid->dqdxtile, vc1, vdqdxc1, 8);
    //             _mm512_i32scatter_pd(grid->dqdytile, vc1, vdqdyc1, 8);
    //             _mm512_i32scatter_pd(grid->dqdztile, vc1, vdqdzc1, 8);
    //             _mm512_i32scatter_pd(grid->dqdxtile, vc2, vdqdxc2, 8);
    //             _mm512_i32scatter_pd(grid->dqdytile, vc2, vdqdyc2, 8);
    //             _mm512_i32scatter_pd(grid->dqdztile, vc2, vdqdzc2, 8);
    //         }
    //         for (IntType vi = i; vi < endFace; vi++)
    //         {
    //             IntType c1 = grid->f2c1[vi];
    //             IntType c2 = grid->f2c2[vi];
    //             RealGeom tmpx = qsumtile[vi] * grid->xfntile[vi];
    //             RealGeom tmpy = qsumtile[vi] * grid->yfntile[vi];
    //             RealGeom tmpz = qsumtile[vi] * grid->zfntile[vi];
    //             // For cell c1
    //             grid->dqdxtile[c1] += tmpx;
    //             grid->dqdytile[c1] += tmpy;
    //             grid->dqdztile[c1] += tmpz;
    //             // For cell c2
    //             grid->dqdxtile[c2] -= tmpx;
    //             grid->dqdytile[c2] -= tmpy;
    //             grid->dqdztile[c2] -= tmpz;
    //         }
    //     }
    //     sdel_array_1D(grid_bfacegroup);
    //     sdel_array_1D(grid_ifacegroup);
    //     // transform back:
    // #ifdef MF_OPENMP
    // #pragma omp parallel for
    // #endif
    //     for (IntType i = 0; i < nTCell; i++)
    //     {
    //         dqdx[i] = grid->dqdxtile[i];
    //         dqdy[i] = grid->dqdytile[i];
    //         dqdz[i] = grid->dqdztile[i];
    //     }

    // #elif (defined MF_OPENMP) && (defined FS_SIMD) && (defined DIVREP) && (defined BoundedColoring) // add by dingxin, 2021-11-30
    //     IntType threads = grid->threads;
    //     IntType startFace, endFace, t, k, iv;
    //     IntType endIndex_iFace_vec;
    //     const IntType Vec = VEC_SIZE;
    //     if (grid->DivRepSuccess)
    //     {
    // #pragma omp parallel for private(t, i, j, startFace, endFace, count, c1, c2, face, type, qsum, tmpx, tmpy, tmpz)
    //         for (t = 0; t < threads; t++)
    //         { // Boundary faces
    //             startFace = grid->idx_pthreads_bface[t];
    //             endFace = grid->idx_pthreads_bface[t + 1];
    //             for (i = startFace; i < endFace; i++)
    //             {
    //                 face = grid->id_division_bface[i];
    //                 count = 2 * face;
    //                 c1 = f2c[count];
    //                 c2 = f2c[count + 1];
    //                 type = bcr[face]->GetType();
    //                 qsum = 0.0;
    //                 if (type == INTERFACE || type == SYMM)
    //                 {
    //                     for (j = 0; j < nNPF[face]; j++)
    //                         qsum += q_n[F2N[face][j]];
    //                     qsum /= RealFlow(nNPF[face]);
    //                 }
    //                 else
    //                 {
    //                     qsum = 0.5 * (q[c1] + q[c2]);
    //                 }
    //                 qsum *= area[face];
    //                 tmpx = qsum * xfn[face];
    //                 tmpy = qsum * yfn[face];
    //                 tmpz = qsum * zfn[face];
    //                 dqdx[c1] += tmpx;
    //                 dqdy[c1] += tmpy;
    //                 dqdz[c1] += tmpz;
    //             }
    //         }
    // #pragma omp parallel for private(t, i, j, k, iv, startFace, endFace, endIndex_iFace_vec, count, c1, c2, face, type, qsum, tmpx, tmpy, tmpz)
    //         for (t = 0; t < threads; t++)
    //         {
    //             // Interior faces
    //             startFace = grid->idx_pthreads_iface[t];
    //             endFace = grid->idx_pthreads_iface[t + 1];
    //             endIndex_iFace_vec = grid->endIndex_iFace_vec[t];
    //             for (k = startFace; k + Vec <= endIndex_iFace_vec; k += Vec)
    //             {
    //                 IntType c1v[Vec], c2v[Vec], countv[Vec], tag[Vec], facev[Vec];
    //                 RealGeom xfnv[Vec], yfnv[Vec], zfnv[Vec], areav[Vec];
    //                 RealGeom qsumv[Vec];
    //                 RealGeom dqdxc1v[Vec], dqdxc2v[Vec];
    //                 RealGeom dqdyc1v[Vec], dqdyc2v[Vec];
    //                 RealGeom dqdzc1v[Vec], dqdzc2v[Vec];
    //                 RealGeom tmpxv[Vec], tmpyv[Vec], tmpzv[Vec];
    //                 RealGeom nNPFv[Vec];
    //                 // Load:
    //                 for (iv = 0; iv < Vec; iv++)
    //                 {
    //                     tag[iv] = grid->id_division_iface[iv + k];
    //                     if (abs(tag[iv]) < nTFace)
    //                         facev[iv] = tag[iv];
    //                     else
    //                         facev[iv] = abs(tag[iv]) - nTFace;
    //                 }
    // #pragma omp simd safelen(Vec)
    //                 for (iv = 0; iv < Vec; iv++)
    //                 {
    //                     countv[iv] = 2 * facev[iv];
    //                     c1v[iv] = f2c[countv[iv]];
    //                     c2v[iv] = f2c[countv[iv] + 1];
    //                     qsumv[iv] = 0.0;
    //                     areav[iv] = area[facev[iv]];
    //                     xfnv[iv] = xfn[facev[iv]];
    //                     yfnv[iv] = yfn[facev[iv]];
    //                     zfnv[iv] = zfn[facev[iv]];
    //                     nNPFv[iv] = nNPF[facev[iv]];
    //                     if (abs(tag[iv]) < nTFace)
    //                     {
    //                         dqdxc1v[iv] = dqdx[c1v[iv]];
    //                         dqdyc1v[iv] = dqdy[c1v[iv]];
    //                         dqdzc1v[iv] = dqdz[c1v[iv]];
    //                         dqdxc2v[iv] = dqdx[c2v[iv]];
    //                         dqdyc2v[iv] = dqdy[c2v[iv]];
    //                         dqdzc2v[iv] = dqdz[c2v[iv]];
    //                     }
    //                     else
    //                     {
    //                         if (tag[iv] > 0)
    //                         {
    //                             dqdxc1v[iv] = dqdx[c1v[iv]];
    //                             dqdyc1v[iv] = dqdy[c1v[iv]];
    //                             dqdzc1v[iv] = dqdz[c1v[iv]];
    //                         }
    //                         else
    //                         {
    //                             dqdxc2v[iv] = dqdx[c2v[iv]];
    //                             dqdyc2v[iv] = dqdy[c2v[iv]];
    //                             dqdzc2v[iv] = dqdz[c2v[iv]];
    //                         }
    //                     }
    //                 }
    //                 // Computation:
    //                 for (iv = 0; iv < Vec; iv++)
    //                 {
    //                     for (j = 0; j < nNPFv[iv]; j++)
    //                     {
    //                         qsumv[iv] += q_n[F2N[facev[iv]][j]];
    //                     }
    //                 }
    // #pragma omp simd safelen(Vec)
    //                 for (iv = 0; iv < Vec; iv++)
    //                 {
    //                     qsumv[iv] /= RealFlow(nNPFv[iv]);
    //                     qsumv[iv] *= areav[iv];
    //                     tmpxv[iv] = qsumv[iv] * xfnv[iv];
    //                     tmpyv[iv] = qsumv[iv] * yfnv[iv];
    //                     tmpzv[iv] = qsumv[iv] * zfnv[iv];
    //                     if (abs(tag[iv]) < nTFace)
    //                     {
    //                         dqdxc1v[iv] += tmpxv[iv];
    //                         dqdyc1v[iv] += tmpyv[iv];
    //                         dqdzc1v[iv] += tmpzv[iv];
    //                         dqdxc2v[iv] -= tmpxv[iv];
    //                         dqdyc2v[iv] -= tmpyv[iv];
    //                         dqdzc2v[iv] -= tmpzv[iv];
    //                     }
    //                     else
    //                     {
    //                         if (tag[iv] > 0)
    //                         {
    //                             dqdxc1v[iv] += tmpxv[iv];
    //                             dqdyc1v[iv] += tmpyv[iv];
    //                             dqdzc1v[iv] += tmpzv[iv];
    //                         }
    //                         else
    //                         {
    //                             dqdxc2v[iv] -= tmpxv[iv];
    //                             dqdyc2v[iv] -= tmpyv[iv];
    //                             dqdzc2v[iv] -= tmpzv[iv];
    //                         }
    //                     }
    //                 }
    //                 // Load Back:
    // #pragma omp simd safelen(Vec)
    //                 for (iv = 0; iv < Vec; iv++)
    //                 {
    //                     if (abs(tag[iv]) < nTFace)
    //                     {
    //                         dqdx[c1v[iv]] = dqdxc1v[iv];
    //                         dqdy[c1v[iv]] = dqdyc1v[iv];
    //                         dqdz[c1v[iv]] = dqdzc1v[iv];
    //                         dqdx[c2v[iv]] = dqdxc2v[iv];
    //                         dqdy[c2v[iv]] = dqdyc2v[iv];
    //                         dqdz[c2v[iv]] = dqdzc2v[iv];
    //                     }
    //                     else
    //                     {
    //                         if (tag[iv] > 0)
    //                         {
    //                             dqdx[c1v[iv]] = dqdxc1v[iv];
    //                             dqdy[c1v[iv]] = dqdyc1v[iv];
    //                             dqdz[c1v[iv]] = dqdzc1v[iv];
    //                         }
    //                         else
    //                         {
    //                             dqdx[c2v[iv]] = dqdxc2v[iv];
    //                             dqdy[c2v[iv]] = dqdyc2v[iv];
    //                             dqdz[c2v[iv]] = dqdzc2v[iv];
    //                         }
    //                     }
    //                 }
    //             }
    //             startFace = k;
    //             for (i = startFace; i < endFace; i++)
    //             {
    //                 k = grid->id_division_iface[i];
    //                 if (abs(k) < nTFace)
    //                     face = k;
    //                 else
    //                     face = abs(k) - nTFace;
    //                 count = 2 * face;
    //                 c1 = f2c[count];
    //                 c2 = f2c[count + 1];
    //                 qsum = 0.0;
    //                 for (j = 0; j < nNPF[face]; j++)
    //                     qsum += q_n[F2N[face][j]];
    //                 qsum /= RealFlow(nNPF[face]);
    //                 qsum *= area[face];
    //                 tmpx = qsum * xfn[face];
    //                 tmpy = qsum * yfn[face];
    //                 tmpz = qsum * zfn[face];
    //                 if (abs(k) < nTFace)
    //                 {
    //                     dqdx[c1] += tmpx;
    //                     dqdy[c1] += tmpy;
    //                     dqdz[c1] += tmpz;
    //                     dqdx[c2] -= tmpx;
    //                     dqdy[c2] -= tmpy;
    //                     dqdz[c2] -= tmpz;
    //                 }
    //                 else
    //                 {
    //                     if (k > 0)
    //                     {
    //                         dqdx[c1] += tmpx;
    //                         dqdy[c1] += tmpy;
    //                         dqdz[c1] += tmpz;
    //                     }
    //                     else
    //                     {
    //                         dqdx[c2] -= tmpx;
    //                         dqdy[c2] -= tmpy;
    //                         dqdz[c2] -= tmpz;
    //                     }
    //                 }
    //             }
    //         }
    //     }

    // #elif (defined MF_OPENMP) && (defined FS_SIMD) && (defined DIVREP) // add by dingxin, 2021-12-05
    //     IntType threads = grid->threads;
    //     IntType startFace, endFace, t, k, iv;
    //     const IntType Vec = VEC_SIZE;
    //     if (grid->DivRepSuccess)
    //     {
    // #pragma omp parallel for private(t, i, j, startFace, endFace, count, c1, c2, face, type, qsum, tmpx, tmpy, tmpz)
    //         for (t = 0; t < threads; t++)
    //         { // Boundary faces
    //             startFace = grid->idx_pthreads_bface[t];
    //             endFace = grid->idx_pthreads_bface[t + 1];
    //             for (i = startFace; i < endFace; i++)
    //             {
    //                 face = grid->id_division_bface[i];
    //                 count = 2 * face;
    //                 c1 = f2c[count];
    //                 c2 = f2c[count + 1];
    //                 type = bcr[face]->GetType();
    //                 qsum = 0.0;
    //                 if (type == INTERFACE || type == SYMM)
    //                 {
    //                     for (j = 0; j < nNPF[face]; j++)
    //                         qsum += q_n[F2N[face][j]];
    //                     qsum /= RealFlow(nNPF[face]);
    //                 }
    //                 else
    //                 {
    //                     qsum = 0.5 * (q[c1] + q[c2]);
    //                 }
    //                 qsum *= area[face];
    //                 tmpx = qsum * xfn[face];
    //                 tmpy = qsum * yfn[face];
    //                 tmpz = qsum * zfn[face];
    //                 dqdx[c1] += tmpx;
    //                 dqdy[c1] += tmpy;
    //                 dqdz[c1] += tmpz;
    //             }
    //         }
    // #pragma omp parallel for private(t, i, j, k, iv, startFace, endFace, count, c1, c2, face, type, qsum, tmpx, tmpy, tmpz)
    //         for (t = 0; t < threads; t++)
    //         {
    //             // Interior faces
    //             IntType c1v[Vec], c2v[Vec], countv[Vec], tag[Vec], facev[Vec];
    //             RealGeom xfnv[Vec], yfnv[Vec], zfnv[Vec], areav[Vec];
    //             RealGeom qsumv[Vec];
    //             RealGeom tmpxv[Vec], tmpyv[Vec], tmpzv[Vec];
    //             RealGeom nNPFv[Vec];
    //             startFace = grid->idx_pthreads_iface[t];
    //             endFace = grid->idx_pthreads_iface[t + 1];
    //             for (k = startFace; k + Vec < endFace; k += Vec)
    //             {
    //                 // Load:
    //                 for (iv = 0; iv < Vec; iv++)
    //                 {
    //                     tag[iv] = grid->id_division_iface[iv + k];
    //                     if (abs(tag[iv]) < nTFace)
    //                         facev[iv] = tag[iv];
    //                     else
    //                         facev[iv] = abs(tag[iv]) - nTFace;
    //                 }
    // #pragma omp simd safelen(Vec)
    //                 for (iv = 0; iv < Vec; iv++)
    //                 {
    //                     countv[iv] = 2 * facev[iv];
    //                     c1v[iv] = f2c[countv[iv]];
    //                     c2v[iv] = f2c[countv[iv] + 1];
    //                     qsumv[iv] = 0.0;
    //                     areav[iv] = area[facev[iv]];
    //                     xfnv[iv] = xfn[facev[iv]];
    //                     yfnv[iv] = yfn[facev[iv]];
    //                     zfnv[iv] = zfn[facev[iv]];
    //                     nNPFv[iv] = nNPF[facev[iv]];
    //                 }
    //                 // Computation:
    //                 for (iv = 0; iv < Vec; iv++)
    //                 {
    //                     for (j = 0; j < nNPFv[iv]; j++)
    //                     {
    //                         qsumv[iv] += q_n[F2N[facev[iv]][j]];
    //                     }
    //                 }
    // #pragma omp simd safelen(Vec)
    //                 for (iv = 0; iv < Vec; iv++)
    //                 {
    //                     qsumv[iv] /= nNPFv[iv];
    //                     qsumv[iv] *= areav[iv];
    //                     tmpxv[iv] = qsumv[iv] * xfnv[iv];
    //                     tmpyv[iv] = qsumv[iv] * yfnv[iv];
    //                     tmpzv[iv] = qsumv[iv] * zfnv[iv];
    //                 }
    //                 // Load Back: no color, exist conflicts
    //                 for (iv = 0; iv < Vec; iv++)
    //                 {
    //                     if (abs(tag[iv]) < nTFace)
    //                     {
    //                         dqdx[c1v[iv]] += tmpxv[iv];
    //                         dqdy[c1v[iv]] += tmpyv[iv];
    //                         dqdz[c1v[iv]] += tmpzv[iv];
    //                         dqdx[c2v[iv]] -= tmpxv[iv];
    //                         dqdy[c2v[iv]] -= tmpyv[iv];
    //                         dqdz[c2v[iv]] -= tmpzv[iv];
    //                     }
    //                     else
    //                     {
    //                         if (tag[iv] > 0)
    //                         {
    //                             dqdx[c1v[iv]] += tmpxv[iv];
    //                             dqdy[c1v[iv]] += tmpyv[iv];
    //                             dqdz[c1v[iv]] += tmpzv[iv];
    //                         }
    //                         else
    //                         {
    //                             dqdx[c2v[iv]] -= tmpxv[iv];
    //                             dqdy[c2v[iv]] -= tmpyv[iv];
    //                             dqdz[c2v[iv]] -= tmpzv[iv];
    //                         }
    //                     }
    //                 }
    //             }
    //             startFace = k;
    //             for (i = startFace; i < endFace; i++)
    //             {
    //                 k = grid->id_division_iface[i];
    //                 if (abs(k) < nTFace)
    //                     face = k;
    //                 else
    //                     face = abs(k) - nTFace;
    //                 count = 2 * face;
    //                 c1 = f2c[count];
    //                 c2 = f2c[count + 1];
    //                 qsum = 0.0;
    //                 for (j = 0; j < nNPF[face]; j++)
    //                     qsum += q_n[F2N[face][j]];
    //                 qsum /= RealFlow(nNPF[face]);
    //                 qsum *= area[face];
    //                 tmpx = qsum * xfn[face];
    //                 tmpy = qsum * yfn[face];
    //                 tmpz = qsum * zfn[face];
    //                 if (abs(k) < nTFace)
    //                 {
    //                     dqdx[c1] += tmpx;
    //                     dqdy[c1] += tmpy;
    //                     dqdz[c1] += tmpz;
    //                     dqdx[c2] -= tmpx;
    //                     dqdy[c2] -= tmpy;
    //                     dqdz[c2] -= tmpz;
    //                 }
    //                 else
    //                 {
    //                     if (k > 0)
    //                     {
    //                         dqdx[c1] += tmpx;
    //                         dqdy[c1] += tmpy;
    //                         dqdz[c1] += tmpz;
    //                     }
    //                     else
    //                     {
    //                         dqdx[c2] -= tmpx;
    //                         dqdy[c2] -= tmpy;
    //                         dqdz[c2] -= tmpz;
    //                     }
    //                 }
    //             }
    //         }
    //     }

#else
    for (i = 0; i < nBFace; i++)
    {
        count = 2 * i;
        c1 = f2c[count];
        c2 = f2c[count + 1];
        type = bcr[i]->GetType();
        qsum = 0.0;

        if (type == INTERFACE || type == SYMM)
        {
            for (j = 0; j < nNPF[i]; j++)
                qsum += q_n[F2N[i][j]];
            qsum /= RealFlow(nNPF[i]);
        }
        else
        {
            qsum = 0.5 * (q[c1] + q[c2]);
        }

        qsum *= area[i];
        // tmpx = qsum * xfn[i];
        // tmpy = qsum * yfn[i];
        // tmpz = qsum * zfn[i];
        dqdx[c1] += qsum * xfn[i];
        dqdy[c1] += qsum * yfn[i];
        dqdz[c1] += qsum * zfn[i];
    }

    // Interior faces
    for (i = nBFace; i < nTFace; i++)
    {
        count = 2 * i;
        c1 = f2c[count];
        c2 = f2c[count + 1];
        qsum = 0.0;

        for (j = 0; j < nNPF[i]; j++)
            qsum += q_n[F2N[i][j]];
        qsum /= RealFlow(nNPF[i]);

        qsum *= area[i];
        tmpx = qsum * xfn[i];
        tmpy = qsum * yfn[i];
        tmpz = qsum * zfn[i];
        // For cell c1
        dqdx[c1] += tmpx;
        dqdy[c1] += tmpy;
        dqdz[c1] += tmpz;

        // For cell c2
        dqdx[c2] -= tmpx;
        dqdy[c2] -= tmpy;
        dqdz[c2] -= tmpz;
    }
#endif

    // 如果单元含有一个以上的物面，该单元梯度采用Gauss求解
    //     IntType level = 0;
    //     if (vis_mode != INVISCID && level == 0)
    //     {
    //         IntType *cellwallnumber = grid->GetGridQualityCellWallNumber();
    // #ifdef MF_OPENMP
    // #pragma omp parallel for private(i, j, face, c1, c2, qsum, tmpx, tmpy, tmpz)
    // #endif
    //         for (i = 0; i < nTCell; i++)
    //         {
    //             if (cellwallnumber[i] < 2)
    //                 continue;
    //             dqdx[i] = 0.0;
    //             dqdy[i] = 0.0;
    //             dqdz[i] = 0.0;
    //             for (j = 0; j < nFPC[i]; j++)
    //             {
    //                 face = C2F[i][j];
    //                 c1 = f2c[face + face];
    //                 c2 = f2c[face + face + 1];
    //                 qsum = 0.5 * (q[c1] + q[c2]) * area[face];
    //                 tmpx = qsum * xfn[face];
    //                 tmpy = qsum * yfn[face];
    //                 tmpz = qsum * zfn[face];
    //                 if (i == c1)
    //                 {
    //                     dqdx[i] += tmpx;
    //                     dqdy[i] += tmpy;
    //                     dqdz[i] += tmpz;
    //                 }
    //                 else if (i == c2)
    //                 {
    //                     dqdx[i] -= tmpx;
    //                     dqdy[i] -= tmpy;
    //                     dqdz[i] -= tmpz;
    //                 }
    //                 else
    //                 {
    //                     cerr << endl
    //                          << "Error in function CompGradientQ_Gauss_Node_Dist!  i is not c1 or c2! " << i << "  " << c1 << "  " << c2 << endl;
    //                     // mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
    //                     exit(1);
    //                 }
    //             }
    //         }
    //     }

    // 边界层前n层采用Gauss方法
    // IntType GaussLayer = -1;
    // grid->GetData(&GaussLayer, INT, 1, "GaussLayer"); ///< = 5
    // IntType *CellLayerNo = (IntType *)grid->GetDataPtr(INT, n, "CellLayerNo");
    // if (level == 0 && GaussLayer > 0)
    // {
    //     for (i = 0; i < nTCell; i++)
    //     {
    //         if ((CellLayerNo[i] == -1) || (CellLayerNo[i] >= GaussLayer))
    //             continue;
    //         dqdx[i] = 0.0;
    //         dqdy[i] = 0.0;
    //         dqdz[i] = 0.0;
    //         for (j = 0; j < nFPC[i]; j++)
    //         {
    //             face = C2F[i][j];
    //             c1 = f2c[face + face];
    //             c2 = f2c[face + face + 1];
    //             qsum = 0.5 * (q[c1] + q[c2]) * area[face];
    //             tmpx = qsum * xfn[face];
    //             tmpy = qsum * yfn[face];
    //             tmpz = qsum * zfn[face];
    //             if (i == c1)
    //             {
    //                 dqdx[i] += tmpx;
    //                 dqdy[i] += tmpy;
    //                 dqdz[i] += tmpz;
    //             }
    //             else if (i == c2)
    //             {
    //                 dqdx[i] -= tmpx;
    //                 dqdy[i] -= tmpy;
    //                 dqdz[i] -= tmpz;
    //             }
    //             else
    //             {
    //                 cerr << endl
    //                      << "Error in function CompGradientQ_Gauss_Node_Dist!  i is not c1 or c2! " << i << "  " << c1 << "  " << c2 << endl;
    //                 // mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
    //                 exit(1);
    //             }
    //         }
    //     }
    // }

#ifdef MF_OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nTCell; i++)
    {
        dqdx[i] /= vol[i];
        dqdy[i] /= vol[i];
        dqdz[i] /= vol[i];
    }
    sdel_array_1D(q_n);

#ifdef MF_TIMING
#ifdef MF_MPICH
    total_t[2] += time_tmp + MPI_Wtime();
#else
    gettimeofday(&endtimeTemRoe, 0);
    time_tmp = (RealGeom)1000000 * (endtimeTemRoe.tv_sec - starttimeTemRoe.tv_sec) + endtimeTemRoe.tv_usec - starttimeTemRoe.tv_usec;
    total_t[2] += time_tmp;
#endif
#endif
}

/*!
 * @brief       Compute the node variable use the distance weight
 * @param       grid
 * @param       q_n
 * @param       q
 * @param       name
 * @remarks     modify according to the fun [void CompNodeVar3D_dist(PolyGrid *grid, RealFlow *q_n, RealFlow *q, IntType name)]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-08
 */
void CompNodeVar3D_dist(PolyGrid *grid, RealFlow *q_n, RealFlow *q, IntType name, RealFlow *u_n, RealFlow *v_n, RealFlow *w_n)
{
#ifdef MF_TIMING
    double time_tmp = 0.0;
#ifdef MF_MPICH
    MPI_Barrier(MPI_COMM_WORLD);
    time_tmp = -MPI_Wtime();
#else
    struct timeval starttimeTemRoe, endtimeTemRoe;
    gettimeofday(&starttimeTemRoe, 0);
#endif
#endif

    IntType i, j, c1, c2, type, p1;
    IntType nTNode = grid->GetNTNode();
    IntType nBFace = grid->GetNBFace();
    IntType nTCell = grid->GetNTCell();
    IntType n = nTCell + nBFace;

    IntType *f2c = grid->Getf2c();
    IntType *nNPF = grid->GetnNPF();
    IntType *nNPC = CalnNPC(grid);
    IntType **C2N = CalC2N(grid);
    IntType **F2N = CalF2N(grid);
    RealGeom *x = grid->GetX();
    RealGeom *y = grid->GetY();
    RealGeom *z = grid->GetZ();
    RealGeom *xfc = grid->GetXfc();
    RealGeom *yfc = grid->GetYfc();
    RealGeom *zfc = grid->GetZfc();
    RealGeom dx, dy, dz, wr;

    // node to cell connection:
    IntType *nCPN = CalnCPN(grid);
    IntType **N2C = CalN2C(grid);

    BCRecord **bcr = grid->Getbcr();
    // if (grid->GetWeightNodeDist() == NULL)
    // {
    //     ComputeWeight3D_Node(grid); // 距离分之一权
    // }
    RealGeom *WeightNode = grid->GetWeightNodeDist();
    IntType *Nmark = grid->GetNodeType();
    RealGeom **WeightNodeC2N = grid->GetWeightNodeC2N();
    RealGeom **WeightNodeN2C = grid->GetWeightNodeN2C();
    RealGeom **WeightNodeBFace2C = grid->GetWeightNodeBFace2C();

    // RealFlow *rho = grid->GetRho();
    RealFlow *u = grid->GetU();
    RealFlow *v = grid->GetV();
    RealFlow *w = grid->GetW();
    // RealFlow *p = grid->GetP();

    //!< 求顶点速度需要的量
    // RealFlow /* *u, *v, *w, *u_n, *v_n, *w_n,*/ *facu, *facv, *facw;
    IntType *nodesymm = grid->GetNodeSymm();
    if (!nodesymm)
    {
        FindNodeSYMM(grid);
        nodesymm = grid->GetNodeSymm();
    }
    RealGeom *xfn_n_symm, *yfn_n_symm, *zfn_n_symm;
    xfn_n_symm = grid->GetXfnNSymm();
    yfn_n_symm = grid->GetYfnNSymm();
    zfn_n_symm = grid->GetZfnNSymm();
    RealFlow *facu = NULL, *facv = NULL, *facw = NULL;

    //!< At u gradient calculation, work out u_n, v_n and w_n for the v and w gradient calculation
    // if (name > 0 && name < 4)
    if (name == 1)
    {
        // nodesymm = (IntType *)grid->GetDataPtr(INT, nTNode, "nodesymm");
        // if (!nodesymm)
        // {
        //     FindNodeSYMM(grid);
        //     nodesymm = (IntType *)grid->GetDataPtr(INT, nTNode, "nodesymm");
        // }
        // u_n = NULL;
        // v_n = NULL;
        // w_n = NULL;
        // snew_array_1D(u_n, nTNode);
        // snew_array_1D(v_n, nTNode);
        // snew_array_1D(w_n, nTNode);
#ifdef MF_OPENMP
#pragma omp parallel for private(i)
#endif
        for (i = 0; i < nTNode; i++)
        {
            u_n[i] = 0.0;
            v_n[i] = 0.0;
            w_n[i] = 0.0;
        }

        snew_array_1D(facu, nBFace);
        snew_array_1D(facv, nBFace);
        snew_array_1D(facw, nBFace);

        // xfn_n_symm = (RealFlow *)grid->GetDataPtr(REAL_GEOM, nTNode, "xfn_n_symm");
        // yfn_n_symm = (RealFlow *)grid->GetDataPtr(REAL_GEOM, nTNode, "yfn_n_symm");
        // zfn_n_symm = (RealFlow *)grid->GetDataPtr(REAL_GEOM, nTNode, "zfn_n_symm");
    }

#ifdef MF_OPENMP
#pragma omp parallel for private(i)
#endif
    for (i = 0; i < nTNode; i++)
        q_n[i] = 0.0;

    //!< 计算物理边界点的值，使用物理面的值进行加权计算，不包括并行边界和对称面边界
    //!< 计算物理边界面心的物理值
    RealFlow *facq = NULL;
    snew_array_1D(facq, nBFace);

#ifdef MF_OPENMP
#pragma omp parallel for private(i, type, c1, c2)
#endif
    for (i = 0; i < nBFace; i++)
    {
        facq[i] = 0.0;
        type = bcr[i]->GetType();
        if (type == INTERFACE || type == SYMM)
            continue;
        c1 = f2c[2 * i];
        c2 = f2c[2 * i + 1];
        facq[i] = q[c1] + q[c2];
        facq[i] *= 0.5;

        // if (name > 0 && name < 4)
        if (name == 1)
        {
            facu[i] = u[c1] + u[c2];
            facu[i] *= 0.5;
            facv[i] = v[c1] + v[c2];
            facv[i] *= 0.5;
            facw[i] = w[c1] + w[c2];
            facw[i] *= 0.5;
        }
    }

    //!< 利用物理边界面心的值，计算物理边界点

    // #ifdef MF_OPENMP
    // #pragma omp parallel
    //     {
    // #pragma omp for private(j, type, p1, dx, dy, dz, wr)
    //         for (i = 0; i < nBFace; i++)
    //         {
    //             type = bcr[i]->GetType();
    //             if (type != WALL)
    //                 continue;
    //             for (j = 0; j < nNPF[i]; j++)
    //             {
    //                 p1 = F2N[i][j];
    //                 dx = x[p1] - xfc[i];
    //                 dy = y[p1] - yfc[i];
    //                 dz = z[p1] - zfc[i];
    //                 wr = dx * dx + dy * dy + dz * dz;
    //                 wr = sqrt(wr);
    //                 wr = 1. / wr;
    // #pragma omp atomic
    //                 q_n[p1] += facq[i] * wr;
    //                 if (name > 0 && name < 4 && nodesymm[p1] == 1)
    //                 {
    // #pragma omp atomic
    //                     u_n[p1] += facu[i] * wr;
    // #pragma omp atomic
    //                     v_n[p1] += facv[i] * wr;
    // #pragma omp atomic
    //                     w_n[p1] += facw[i] * wr;
    //                 }
    //             }
    //         }
    // #pragma omp for private(j, type, p1, dx, dy, dz, wr)
    //         for (i = 0; i < nBFace; i++)
    //         {
    //             type = bcr[i]->GetType();
    //             if (type != FAR_FIELD && type != FAR_FIELD_MOD)
    //                 continue;
    //             for (j = 0; j < nNPF[i]; j++)
    //             {
    //                 p1 = F2N[i][j];
    //                 if (Nmark[p1] == WALL)
    //                     continue;
    //                 dx = x[p1] - xfc[i];
    //                 dy = y[p1] - yfc[i];
    //                 dz = z[p1] - zfc[i];
    //                 wr = dx * dx + dy * dy + dz * dz;
    //                 wr = sqrt(wr);
    //                 wr = 1. / wr;
    // #pragma omp atomic
    //                 q_n[p1] += facq[i] * wr;=
    //                 if (name > 0 && name < 4 && nodesymm[p1] == 1)
    //                 {
    // #pragma omp atomic
    //                     u_n[p1] += facu[i] * wr;
    // #pragma omp atomic
    //                     v_n[p1] += facv[i] * wr;
    // #pragma omp atomic
    //                     w_n[p1] += facw[i] * wr;
    //                 }
    //             }
    //         }
    // #pragma omp for private(j, type, p1, dx, dy, dz, wr)
    //         for (i = 0; i < nBFace; i++)
    //         {
    //             type = bcr[i]->GetType();
    //             if (type == WALL || type == SYMM || type == FAR_FIELD || type == FAR_FIELD_MOD || type == INTERFACE)
    //                 continue;
    //             for (j = 0; j < nNPF[i]; j++)
    //             {
    //                 p1 = F2N[i][j];
    //                 if (Nmark[p1] == WALL || Nmark[p1] == FAR_FIELD)
    //                     continue;
    //                 dx = x[p1] - xfc[i];
    //                 dy = y[p1] - yfc[i];
    //                 dz = z[p1] - zfc[i];
    //                 wr = dx * dx + dy * dy + dz * dz;
    //                 wr = sqrt(wr);
    //                 wr = 1. / wr;
    // #pragma omp atomic
    //                 q_n[p1] += facq[i] * wr;
    //                 if (name > 0 && name < 4 && nodesymm[p1] == 1)
    //                 {
    // #pragma omp atomic
    //                     u_n[p1] += facu[i] * wr;
    // #pragma omp atomic
    //                     v_n[p1] += facv[i] * wr;
    // #pragma omp atomic
    //                     w_n[p1] += facw[i] * wr;
    //                 }
    //             }
    //         }
    //     }
    // #else
    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        if (type != WALL)
            continue;
        for (j = 0; j < nNPF[i]; j++)
        {
            p1 = F2N[i][j];
            // dx = x[p1] - xfc[i];
            // dy = y[p1] - yfc[i];
            // dz = z[p1] - zfc[i];
            // wr = dx * dx + dy * dy + dz * dz;
            // wr = sqrt(wr);
            // wr = 1. / wr;
            wr = WeightNodeBFace2C[i][j];
            q_n[p1] += facq[i] * wr;

            // if (name > 0 && name < 4 && nodesymm[p1] == 1)
            if (name == 1 && nodesymm[p1] == 1)
            {
                u_n[p1] += facu[i] * wr;
                v_n[p1] += facv[i] * wr;
                w_n[p1] += facw[i] * wr;
            }
        }
    }
    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        if (type != FAR_FIELD)
            continue;
        for (j = 0; j < nNPF[i]; j++)
        {
            p1 = F2N[i][j];
            if (Nmark[p1] == WALL)
                continue;
            // dx = x[p1] - xfc[i];
            // dy = y[p1] - yfc[i];
            // dz = z[p1] - zfc[i];
            // wr = dx * dx + dy * dy + dz * dz;
            // wr = sqrt(wr);
            // wr = 1. / wr;
            wr = WeightNodeBFace2C[i][j];
            q_n[p1] += facq[i] * wr;

            // if (name > 0 && name < 4 && nodesymm[p1] == 1)
            if (name == 1 && nodesymm[p1] == 1)
            {
                u_n[p1] += facu[i] * wr;
                v_n[p1] += facv[i] * wr;
                w_n[p1] += facw[i] * wr;
            }
        }
    }
    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        if (type == WALL || type == SYMM || type == FAR_FIELD || type == INTERFACE)
            continue;
        for (j = 0; j < nNPF[i]; j++)
        {
            p1 = F2N[i][j];
            if (Nmark[p1] == WALL || Nmark[p1] == FAR_FIELD)
                continue;
            // dx = x[p1] - xfc[i];
            // dy = y[p1] - yfc[i];
            // dz = z[p1] - zfc[i];
            // wr = dx * dx + dy * dy + dz * dz;
            // wr = sqrt(wr);
            // wr = 1. / wr;
            wr = WeightNodeBFace2C[i][j];
            q_n[p1] += facq[i] * wr;

            // if (name > 0 && name < 4 && nodesymm[p1] == 1)
            if (name == 1 && nodesymm[p1] == 1)
            {
                u_n[p1] += facu[i] * wr;
                v_n[p1] += facv[i] * wr;
                w_n[p1] += facw[i] * wr;
            }
        }
    }
    // #endif

    //!< 计算其他点的物理值，使用与其相相邻的控制体体心值
#ifdef MF_OPENMP
#pragma omp parallel for private(i, j)
    for (i = 0; i < nTNode; i++)
    {
        if (Nmark[i] != 0)
            continue;
        for (j = 0; j < nCPN[i]; j++)
        {
            IntType cellx = N2C[i][j];
            q_n[i] += q[cellx] * WeightNodeN2C[i][j];
            // if (name > 0 && name < 4 && nodesymm[i] == 1)
            if (name == 1 && nodesymm[i] == 1)
            {
                u_n[i] += u[cellx] * WeightNodeN2C[i][j];
                v_n[i] += v[cellx] * WeightNodeN2C[i][j];
                w_n[i] += w[cellx] * WeightNodeN2C[i][j];
            }
        }
    }
#else
    for (i = 0; i < nTCell; i++)
    {
        for (j = 0; j < nNPC[i]; j++)
        {
            p1 = C2N[i][j];
            if (Nmark[p1] != 0)
                continue;

            q_n[p1] += q[i] * WeightNodeC2N[i][j];

            // if (name > 0 && name < 4 && nodesymm[p1] == 1)
            if (name == 1 && nodesymm[p1] == 1)
            {
                u_n[p1] += u[i] * WeightNodeC2N[i][j];
                v_n[p1] += v[i] * WeightNodeC2N[i][j];
                w_n[p1] += w[i] * WeightNodeC2N[i][j];
            }
        }
    }
#endif

#ifdef MF_TIMING
#ifdef MF_MPICH
    total_t[2] += time_tmp + MPI_Wtime();
#else
    gettimeofday(&endtimeTemRoe, 0);
    time_tmp = (RealGeom)1000000 * (endtimeTemRoe.tv_sec - starttimeTemRoe.tv_sec) + endtimeTemRoe.tv_usec - starttimeTemRoe.tv_usec;
    total_t[2] += time_tmp;
#endif
#endif

    //!< 传递并行边界点的加权值
#ifdef MF_MPICH
    grid->CommInternodeDataMPI(q_n);
    if (name == 1)
    {
        grid->CommInternodeDataMPI(u_n);
        grid->CommInternodeDataMPI(v_n);
        grid->CommInternodeDataMPI(w_n);
    }
#endif

#ifdef MF_TIMING
#ifdef MF_MPICH
    MPI_Barrier(MPI_COMM_WORLD);
    time_tmp = -MPI_Wtime();
#else
    gettimeofday(&starttimeTemRoe, 0);
#endif
#endif

#ifdef MF_OPENMP
#pragma omp parallel for private(i)
#endif
    for (i = 0; i < nTNode; i++)
    {
        q_n[i] /= (WeightNode[i] + TINY);
        // if (name > 0 && name < 4 && nodesymm[i] == 1)
        if (name == 1 && nodesymm[i] == 1)
        {
            u_n[i] /= (WeightNode[i] + TINY);
            v_n[i] /= (WeightNode[i] + TINY);
            w_n[i] /= (WeightNode[i] + TINY);
        }
    }

    //!< 修正对称面顶点的速度
    // xfn_n_symm = (RealFlow *)grid->GetDataPtr(REAL_GEOM, nTNode, "xfn_n_symm");
    // yfn_n_symm = (RealFlow *)grid->GetDataPtr(REAL_GEOM, nTNode, "yfn_n_symm");
    // zfn_n_symm = (RealFlow *)grid->GetDataPtr(REAL_GEOM, nTNode, "zfn_n_symm");

    if (name > 0 && name < 4)
    {
#ifdef MF_OPENMP
#pragma omp parallel for private(i)
#endif
        for (i = 0; i < nTNode; i++)
        {
            if (nodesymm[i] != 1)
                continue;

            RealGeom vn = u_n[i] * xfn_n_symm[i] + v_n[i] * yfn_n_symm[i] + w_n[i] * zfn_n_symm[i];
            if (name == 1)
                q_n[i] = u_n[i] - vn * xfn_n_symm[i];
            else if (name == 2)
                q_n[i] = v_n[i] - vn * yfn_n_symm[i];
            else if (name == 3)
                q_n[i] = w_n[i] - vn * zfn_n_symm[i];
        }

        // sdel_array_1D(u_n);
        // sdel_array_1D(v_n);
        // sdel_array_1D(w_n);
        if (name == 1)
        {
            sdel_array_1D(facu);
            sdel_array_1D(facv);
            sdel_array_1D(facw);
        }
    }
    sdel_array_1D(facq);

#ifdef MF_TIMING
#ifdef MF_MPICH
    total_t[2] += time_tmp + MPI_Wtime();
#else
    gettimeofday(&endtimeTemRoe, 0);
    time_tmp = (RealGeom)1000000 * (endtimeTemRoe.tv_sec - starttimeTemRoe.tv_sec) + endtimeTemRoe.tv_usec - starttimeTemRoe.tv_usec;
    total_t[2] += time_tmp;
#endif
#endif
}

#ifdef TDTREE

void CompGradientQ_Gauss_Node_TDTree(PolyGrid *grid, RealFlow **q, RealFlow **dqdx, RealFlow **dqdy, RealFlow **dqdz, IntType ns, IntType ne)
{
    IntType nTNode = grid->GetNTNode();
    IntType nTCell = grid->GetNTCell();
    IntType nTFace = grid->GetNTFace();
    IntType nBFace = grid->GetNBFace();
    IntType *f2c = grid->Getf2c();
    BCRecord **bcr = grid->Getbcr();
    IntType n = nTCell + nBFace;
    RealGeom *vol = grid->GetCellVol();
    RealGeom *area = grid->GetFaceArea();
    RealGeom *xfn = grid->GetXfn();
    RealGeom *yfn = grid->GetYfn();
    RealGeom *zfn = grid->GetZfn();
    IntType *nNPF = grid->GetnNPF();
    IntType **F2N = CalF2N(grid);
    IntType *nFPC = CalnFPC(grid);
    IntType **C2F = CalC2F(grid);

    IntType i, j, c1, c2, count, type, face;
    RealGeom tmpx, tmpy, tmpz;
    RealFlow qsum;

    // Initialize dq
    // #pragma omp parallel for collapse(2)
    //     for (int j = ns; j < ne; j++)
    //         for (i = 0; i < n; i++)
    //         {
    //             dqdx[j][i] = 0.;
    //             dqdy[j][i] = 0.;
    //             dqdz[j][i] = 0.;
    //         }

    RealFlow **q_n = NULL, **u_n = NULL, **v_n = NULL, **w_n = NULL;
    snew_array_2D(q_n, ne - ns, nTNode);
    snew_array_2D(u_n, ne - ns, nTNode);
    snew_array_2D(v_n, ne - ns, nTNode);
    snew_array_2D(w_n, ne - ns, nTNode);
    // CompNodeVar3D_dist(grid, q_n, q, name);
    CompNodeVar3D_dist_TDTree(grid, q_n, q, u_n, v_n, w_n, ns, ne);

    IntType *nodesymm = grid->GetNodeSymm();
    if (!nodesymm)
    {
        FindNodeSYMM(grid);
        nodesymm = grid->GetNodeSymm();
    }
    RealGeom *xfn_n_symm, *yfn_n_symm, *zfn_n_symm;
    xfn_n_symm = grid->GetXfnNSymm();
    yfn_n_symm = grid->GetYfnNSymm();
    zfn_n_symm = grid->GetZfnNSymm();

    // IntType vis_mode;
    // grid->GetData(&vis_mode, INT, 1, "vis_mode");

    // IntType GaussLayer = -1;
    // grid->GetData(&GaussLayer, INT, 1, "GaussLayer");
    // IntType *CellLayerNo = (IntType *)grid->GetDataPtr(INT, n, "CellLayerNo");

    // IntType level;
    // level = grid->GetLevel();

    // if (level == 0)

    // char *userArgs[21] = {(char *)grid, (char *)q, (char *)dqdx, (char *)dqdy, (char *)dqdz, (char *)&ns, (char *)&ne, (char *)q_n, (char *)F2N,
    //                       (char *)nFPC, (char *)C2F, (char *)&vis_mode, (char *)&GaussLayer, (char *)CellLayerNo,
    //                       (char *)u_n, (char *)v_n, (char *)w_n, (char *)xfn_n_symm, (char *)yfn_n_symm, (char *)zfn_n_symm, (char *)nodesymm};
    char *userArgs[21] = {(char *)grid, (char *)q, (char *)dqdx, (char *)dqdy, (char *)dqdz, (char *)&ns, (char *)&ne, (char *)q_n, (char *)F2N,
                          (char *)nFPC, (char *)C2F, (char *)&vis_mode, NULL, NULL,
                          (char *)u_n, (char *)v_n, (char *)w_n, (char *)xfn_n_symm, (char *)yfn_n_symm, (char *)zfn_n_symm, (char *)nodesymm};

    TDTree *TDTreeRoot = grid->GetTDTree();

#ifdef MF_TIMING
    double time_tmp = 0.0;
#ifdef MF_MPICH
    MPI_Barrier(MPI_COMM_WORLD);
    time_tmp = -MPI_Wtime();
#else
    struct timeval starttimeTemRoe, endtimeTemRoe;
    gettimeofday(&starttimeTemRoe, 0);
#endif
#endif

    TDTreeRoot->task_traversal(CompGradientQ_Gauss_Node_Kernel3, NULL, userArgs, Backward);

    TDTreeRoot->task_traversal(CompGradientQ_Gauss_Node_Kernel4, NULL, userArgs, Forward);

#ifdef MF_TIMING
#ifdef MF_MPICH
    total_t[2] += time_tmp + MPI_Wtime();
#else
    gettimeofday(&endtimeTemRoe, 0);
    time_tmp = (RealGeom)1000000 * (endtimeTemRoe.tv_sec - starttimeTemRoe.tv_sec) + endtimeTemRoe.tv_usec - starttimeTemRoe.tv_usec;
    total_t[2] += time_tmp;
#endif
#endif

#ifdef MF_MPICH
    RealFlow **grad_mpi = NULL;
    snew_array_1D(grad_mpi, 3 * (ne - ns));
    count = 0;
    for (IntType i = 0; i < (ne - ns); ++i)
    {
        grad_mpi[count++] = dqdx[i];
        grad_mpi[count++] = dqdy[i];
        grad_mpi[count++] = dqdz[i];
    }
    grid->RecvSendVarNeighbor_Togeth(3 * (ne - ns), grad_mpi);
    // for (IntType i = 0; i < 3 * (ne - ns); ++i)
    // {
    //     grid->CommInterfaceDataMPI(grad_mpi[i]);
    // }
    sdel_array_1D(grad_mpi);
#endif

    sdel_array_2D(q_n);
    sdel_array_2D(u_n);
    sdel_array_2D(v_n);
    sdel_array_2D(w_n);
}

/*!
 * @brief
 * @param       grid
 * @param       q_n
 * @param       q
 * @param       u_n
 * @param       v_n
 * @param       w_n
 * @param       ns
 * @param       ne
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-12
 */
void CompNodeVar3D_dist_TDTree(PolyGrid *grid, RealFlow **q_n, RealFlow **q, RealFlow **u_n, RealFlow **v_n, RealFlow **w_n, IntType ns, IntType ne)
{

    IntType i, j, c1, c2, type, p1, name;
    IntType nTNode = grid->GetNTNode();
    IntType nBFace = grid->GetNBFace();
    IntType nTCell = grid->GetNTCell();
    IntType n = nTCell + nBFace;

    IntType *f2c = grid->Getf2c();
    IntType *nNPF = grid->GetnNPF();
    IntType *nNPC = CalnNPC(grid);
    IntType **C2N = CalC2N(grid);
    IntType **F2N = CalF2N(grid);
    RealGeom *x = grid->GetX();
    RealGeom *y = grid->GetY();
    RealGeom *z = grid->GetZ();
    RealGeom *xfc = grid->GetXfc();
    RealGeom *yfc = grid->GetYfc();
    RealGeom *zfc = grid->GetZfc();
    RealGeom dx, dy, dz, wr;

    // node to cell connection:
    IntType *nCPN = CalnCPN(grid);
    IntType **N2C = CalN2C(grid);

    BCRecord **bcr = grid->Getbcr();
    // if (grid->GetWeightNodeDist() == NULL)
    // {
    //     ComputeWeight3D_Node(grid); // 距离分之一权
    // }
    RealGeom *WeightNode = grid->GetWeightNodeDist();
    IntType *Nmark = grid->GetNodeType();
    RealGeom **WeightNodeC2N = grid->GetWeightNodeC2N();
    RealGeom **WeightNodeN2C = grid->GetWeightNodeN2C();
    RealGeom **WeightNodeBFace2C = grid->GetWeightNodeBFace2C();

    // RealFlow *rho = grid->GetRho();
    RealFlow *u = grid->GetU();
    RealFlow *v = grid->GetV();
    RealFlow *w = grid->GetW();
    // RealFlow *p = grid->GetP();

    // 求顶点速度需要的量
    IntType *nodesymm = grid->GetNodeSymm();
    if (!nodesymm)
    {
        FindNodeSYMM(grid);
        nodesymm = grid->GetNodeSymm();
    }
    RealGeom *xfn_n_symm, *yfn_n_symm, *zfn_n_symm;
    xfn_n_symm = grid->GetXfnNSymm();
    yfn_n_symm = grid->GetYfnNSymm();
    zfn_n_symm = grid->GetZfnNSymm();

    // RealFlow *u, *v, *w, *facu, *facv, *facw;
    RealFlow *facu = NULL, *facv = NULL, *facw = NULL;
    RealFlow vn;
    // u = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "u");
    // v = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "v");
    // w = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "w");

    snew_array_1D(facu, nBFace);
    snew_array_1D(facv, nBFace);
    snew_array_1D(facw, nBFace);

    // 计算物理边界点的值，使用物理面的值进行加权计算，不包括并行边界和对称面边界
    // 计算物理边界面心的物理值
    RealFlow *facq = NULL;
    snew_array_1D(facq, nBFace);

    TDTree *TDTreeRoot = grid->GetTDTree();

    char *userArgs[21] = {(char *)grid, (char *)q, (char *)u, (char *)v, (char *)w, (char *)u_n, (char *)v_n, (char *)w_n, (char *)&ns, (char *)&ne,
                          (char *)q_n, (char *)nNPC, (char *)C2N, (char *)F2N, (char *)nCPN, (char *)N2C, (char *)facu, (char *)facv, (char *)facw,
                          (char *)facq, (char *)nodesymm};

#ifdef MF_TIMING
    double time_tmp = 0.0;
#ifdef MF_MPICH
    MPI_Barrier(MPI_COMM_WORLD);
    time_tmp = -MPI_Wtime();
#else
    struct timeval starttimeTemRoe, endtimeTemRoe;
    gettimeofday(&starttimeTemRoe, 0);
#endif
#endif

    TDTreeRoot->task_traversal(CompNodeVar3D_dist_Kernel2, NULL, userArgs, Backward);

#ifdef MF_TIMING
#ifdef MF_MPICH
    total_t[2] += time_tmp + MPI_Wtime();
#else
    gettimeofday(&endtimeTemRoe, 0);
    time_tmp = (RealGeom)1000000 * (endtimeTemRoe.tv_sec - starttimeTemRoe.tv_sec) + endtimeTemRoe.tv_usec - starttimeTemRoe.tv_usec;
    total_t[2] += time_tmp;
#endif
#endif

    // 传递并行边界点的加权值
#ifdef MF_MPICH
    for (name = ns; name < ne; name++)
    {
        grid->CommInternodeDataMPI(q_n[name]);
        if (name > 0 && name < 4)
        {
            grid->CommInternodeDataMPI(u_n[name]);
            grid->CommInternodeDataMPI(v_n[name]);
            grid->CommInternodeDataMPI(w_n[name]);
        }
    }
#endif

    sdel_array_1D(facu);
    sdel_array_1D(facv);
    sdel_array_1D(facw);
    sdel_array_1D(facq);
}

/*!
 * @brief
 * @param       userArgs
 * @param       treeArgs
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-12
 */
void CompNodeVar3D_dist_Kernel2(char **userArgs, TDTreeArgs *treeArgs)
{
    // #pragma omp critical
    //     {
    PolyGrid *grid = (PolyGrid *)userArgs[0];
    RealFlow **q = (RealFlow **)userArgs[1];
    RealFlow *u = (RealFlow *)userArgs[2];
    RealFlow *v = (RealFlow *)userArgs[3];
    RealFlow *w = (RealFlow *)userArgs[4];
    RealFlow **u_n = (RealFlow **)userArgs[5];
    RealFlow **v_n = (RealFlow **)userArgs[6];
    RealFlow **w_n = (RealFlow **)userArgs[7];
    IntType ns = *(IntType *)userArgs[8];
    IntType ne = *(IntType *)userArgs[9];
    RealFlow **q_n = (RealFlow **)userArgs[10];
    IntType *nNPC = (IntType *)userArgs[11];
    IntType **C2N = (IntType **)userArgs[12];
    IntType **F2N = (IntType **)userArgs[13];
    IntType *nCPN = (IntType *)userArgs[14];
    IntType **N2C = (IntType **)userArgs[15];
    RealFlow *facu = (RealFlow *)userArgs[16];
    RealFlow *facv = (RealFlow *)userArgs[17];
    RealFlow *facw = (RealFlow *)userArgs[18];
    RealFlow *facq = (RealFlow *)userArgs[19];
    IntType *nodesymm = (IntType *)userArgs[20];

    IntType i, j, c1, c2, type, p1, name;
    IntType nTNode = grid->GetNTNode();
    IntType nBFace = grid->GetNBFace();
    IntType nTCell = grid->GetNTCell();
    IntType n = nTCell + nBFace;

    IntType *f2c = grid->Getf2c();
    IntType *nNPF = grid->GetnNPF();
    RealGeom *x = grid->GetX();
    RealGeom *y = grid->GetY();
    RealGeom *z = grid->GetZ();
    RealGeom *xfc = grid->GetXfc();
    RealGeom *yfc = grid->GetYfc();
    RealGeom *zfc = grid->GetZfc();
    RealGeom dx, dy, dz, wr;

    IntType nsCell = treeArgs->firstCell;
    IntType neCell = treeArgs->lastCell + 1;
    IntType nsFace = treeArgs->firstFace;
    IntType neFace = treeArgs->lastFace + 1;
    IntType nsNode = treeArgs->firstNode;
    IntType neNode = treeArgs->lastNode + 1;

    BCRecord **bcr = grid->Getbcr();
    RealGeom *WeightNode = grid->GetWeightNodeDist();
    IntType *Nmark = grid->GetNodeType();
    RealGeom **WeightNodeC2N = grid->GetWeightNodeC2N();
    RealGeom **WeightNodeN2C = grid->GetWeightNodeN2C();
    RealGeom **WeightNodeBFace2C = grid->GetWeightNodeBFace2C();

    for (name = ns; name < ne; name++)
    {

        for (i = nsNode; i < neNode; i++)
        {
            u_n[name][i] = 0.0;
            v_n[name][i] = 0.0;
            w_n[name][i] = 0.0;
            q_n[name][i] = 0.0;
        }
        if (neFace <= nBFace)
        {

            for (i = nsFace; i < neFace; i++)
            {
                facq[i] = 0.0;
                type = bcr[i]->GetType();
                if (type == INTERFACE || type == SYMM)
                    continue;
                c1 = f2c[2 * i];
                c2 = f2c[2 * i + 1];
                facq[i] = q[name][c1] + q[name][c2];
                facq[i] *= 0.5;

                if (name > 0 && name < 4)
                {
                    facu[i] = u[c1] + u[c2];
                    facu[i] *= 0.5;
                    facv[i] = v[c1] + v[c2];
                    facv[i] *= 0.5;
                    facw[i] = w[c1] + w[c2];
                    facw[i] *= 0.5;
                }
            }

            // 利用物理边界面心的值，计算物理边界点
            for (i = nsFace; i < neFace; i++)
            {
                type = bcr[i]->GetType();
                if (type != WALL)
                    continue;
                for (j = 0; j < nNPF[i]; j++)
                {
                    p1 = F2N[i][j];
                    // dx = x[p1] - xfc[i];
                    // dy = y[p1] - yfc[i];
                    // dz = z[p1] - zfc[i];
                    // wr = dx * dx + dy * dy + dz * dz;
                    // wr = sqrt(wr);
                    // wr = 1. / wr;
                    wr = WeightNodeBFace2C[i][j];
                    q_n[name][p1] += facq[i] * wr;

                    if (name > 0 && name < 4 && nodesymm[p1] == 1)
                    {
                        u_n[name][p1] += facu[i] * wr;
                        v_n[name][p1] += facv[i] * wr;
                        w_n[name][p1] += facw[i] * wr;
                    }
                }
            }
            for (i = nsFace; i < neFace; i++)
            {
                type = bcr[i]->GetType();
                if (type != FAR_FIELD && type != FAR_FIELD_MOD)
                    continue;
                for (j = 0; j < nNPF[i]; j++)
                {
                    p1 = F2N[i][j];
                    if (Nmark[p1] == WALL)
                        continue;
                    // dx = x[p1] - xfc[i];
                    // dy = y[p1] - yfc[i];
                    // dz = z[p1] - zfc[i];
                    // wr = dx * dx + dy * dy + dz * dz;
                    // wr = sqrt(wr);
                    // wr = 1. / wr;
                    wr = WeightNodeBFace2C[i][j];
                    q_n[name][p1] += facq[i] * wr;

                    if (name > 0 && name < 4 && nodesymm[p1] == 1)
                    {
                        u_n[name][p1] += facu[i] * wr;
                        v_n[name][p1] += facv[i] * wr;
                        w_n[name][p1] += facw[i] * wr;
                    }
                }
            }
            for (i = nsFace; i < neFace; i++)
            {
                type = bcr[i]->GetType();
                if (type == WALL || type == SYMM || type == FAR_FIELD || type == FAR_FIELD_MOD || type == INTERFACE)
                    continue;
                for (j = 0; j < nNPF[i]; j++)
                {
                    p1 = F2N[i][j];
                    if (Nmark[p1] == WALL || Nmark[p1] == FAR_FIELD)
                        continue;
                    // dx = x[p1] - xfc[i];
                    // dy = y[p1] - yfc[i];
                    // dz = z[p1] - zfc[i];
                    // wr = dx * dx + dy * dy + dz * dz;
                    // wr = sqrt(wr);
                    // wr = 1. / wr;
                    wr = WeightNodeBFace2C[i][j];
                    q_n[name][p1] += facq[i] * wr;

                    if (name > 0 && name < 4 && nodesymm[p1] == 1)
                    {
                        u_n[name][p1] += facu[i] * wr;
                        v_n[name][p1] += facv[i] * wr;
                        w_n[name][p1] += facw[i] * wr;
                    }
                }
            }
        }

        // 计算其他点的物理值，使用与其相相邻的控制体体心值

        for (i = nsNode; i < neNode; i++)
        {
            if (Nmark[i] != 0)
                continue;
            for (j = 0; j < nCPN[i]; j++)
            {
                IntType cellx = N2C[i][j];
                q_n[name][i] += q[name][cellx] * WeightNodeN2C[i][j];
                if (name > 0 && name < 4 && nodesymm[i] == 1)
                {
                    u_n[name][i] += u[cellx] * WeightNodeN2C[i][j];
                    v_n[name][i] += v[cellx] * WeightNodeN2C[i][j];
                    w_n[name][i] += w[cellx] * WeightNodeN2C[i][j];
                }
            }
        }
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
void CompGradientQ_Gauss_Node_Kernel3(char **userArgs, TDTreeArgs *treeArgs)
{
    // #pragma omp critical
    //     {
    PolyGrid *grid = (PolyGrid *)userArgs[0];
    RealFlow **q = (RealFlow **)userArgs[1];
    RealFlow **dqdx = (RealFlow **)userArgs[2];
    RealFlow **dqdy = (RealFlow **)userArgs[3];
    RealFlow **dqdz = (RealFlow **)userArgs[4];
    IntType ns = *(IntType *)userArgs[5];
    IntType ne = *(IntType *)userArgs[6];
    RealFlow **q_n = (RealFlow **)userArgs[7];
    IntType **F2N = (IntType **)userArgs[8];
    IntType *nFPC = (IntType *)userArgs[9];
    IntType **C2F = (IntType **)userArgs[10];
    RealFlow **u_n = (RealFlow **)userArgs[14];
    RealFlow **v_n = (RealFlow **)userArgs[15];
    RealFlow **w_n = (RealFlow **)userArgs[16];
    RealFlow *xfn_n_symm = (RealFlow *)userArgs[17];
    RealFlow *yfn_n_symm = (RealFlow *)userArgs[18];
    RealFlow *zfn_n_symm = (RealFlow *)userArgs[19];
    IntType *nodesymm = (IntType *)userArgs[20];

    IntType nTNode = grid->GetNTNode();
    IntType nTCell = grid->GetNTCell();
    IntType nTFace = grid->GetNTFace();
    IntType nBFace = grid->GetNBFace();
    IntType *f2c = grid->Getf2c();
    BCRecord **bcr = grid->Getbcr();
    IntType n = nTCell + nBFace;
    RealGeom *vol = grid->GetCellVol();
    RealGeom *area = grid->GetFaceArea();
    RealGeom *xfn = grid->GetXfn();
    RealGeom *yfn = grid->GetYfn();
    RealGeom *zfn = grid->GetZfn();
    IntType *nNPF = grid->GetnNPF();
    RealGeom *WeightNode = grid->GetWeightNodeDist();

    IntType nsCell = treeArgs->firstCell;
    IntType neCell = treeArgs->lastCell + 1;
    IntType nsFace = treeArgs->firstFace;
    IntType neFace = treeArgs->lastFace + 1;
    IntType nsNode = treeArgs->firstNode;
    IntType neNode = treeArgs->lastNode + 1;

    IntType i, j, c1, c2, type, face, name;
    RealGeom tmpx, tmpy, tmpz;
    RealFlow qsum;
    RealFlow vn;

    for (name = ns; name < ne; name++)
    {
        for (i = nsNode; i < neNode; i++)
        {
            q_n[name][i] /= (WeightNode[i] + TINY);
            if (name > 0 && name < 4 && nodesymm[i] == 1)
            {
                u_n[name][i] /= (WeightNode[i] + TINY);
                v_n[name][i] /= (WeightNode[i] + TINY);
                w_n[name][i] /= (WeightNode[i] + TINY);
                vn = u_n[name][i] * xfn_n_symm[i] + v_n[name][i] * yfn_n_symm[i] + w_n[name][i] * zfn_n_symm[i];
                if (name == 1)
                    q_n[name][i] = u_n[name][i] - vn * xfn_n_symm[i];
                else if (name == 2)
                    q_n[name][i] = v_n[name][i] - vn * yfn_n_symm[i];
                else if (name == 3)
                    q_n[name][i] = w_n[name][i] - vn * zfn_n_symm[i];
            }
        }
        // count = 2*nsFace;
        if (neFace <= nBFace)
        {
            for (i = nsFace; i < neFace; i++)
            {
                c1 = f2c[i + i];
                c2 = f2c[i + i + 1];
                type = bcr[i]->GetType();
                qsum = 0.0;

                if (type == INTERFACE || type == SYMM)
                {
                    for (j = 0; j < nNPF[i]; j++)
                        qsum += q_n[name][F2N[i][j]];
                    qsum /= RealFlow(nNPF[i]);
                }
                else
                {
                    qsum = 0.5 * (q[name][c1] + q[name][c2]);
                }

                qsum *= area[i];
                tmpx = qsum * xfn[i];
                tmpy = qsum * yfn[i];
                tmpz = qsum * zfn[i];
                dqdx[name][c1] += tmpx;
                dqdy[name][c1] += tmpy;
                dqdz[name][c1] += tmpz;
            }
        }
        else
        {
            // Interior faces
            for (i = nsFace; i < neFace; i++)
            {
                c1 = f2c[i + i];
                c2 = f2c[i + i + 1];
                qsum = 0.0;

                for (j = 0; j < nNPF[i]; j++)
                    qsum += q_n[name][F2N[i][j]];
                qsum /= RealFlow(nNPF[i]);

                qsum *= area[i];
                tmpx = qsum * xfn[i];
                tmpy = qsum * yfn[i];
                tmpz = qsum * zfn[i];
                // For cell c1
                dqdx[name][c1] += tmpx;
                dqdy[name][c1] += tmpy;
                dqdz[name][c1] += tmpz;

                // For cell c2
                dqdx[name][c2] -= tmpx;
                dqdy[name][c2] -= tmpy;
                dqdz[name][c2] -= tmpz;
            }
        }
    }
    // }
}

void CompGradientQ_Gauss_Node_Kernel4(char **userArgs, TDTreeArgs *treeArgs)
{
    // #pragma omp critical
    // {
    PolyGrid *grid = (PolyGrid *)userArgs[0];
    // RealFlow **q = (RealFlow **)userArgs[1];
    RealFlow **dqdx = (RealFlow **)userArgs[2];
    RealFlow **dqdy = (RealFlow **)userArgs[3];
    RealFlow **dqdz = (RealFlow **)userArgs[4];
    IntType ns = *(IntType *)userArgs[5];
    IntType ne = *(IntType *)userArgs[6];
    // RealFlow **q_n = (RealFlow **)userArgs[7];
    // IntType **F2N = (IntType **)userArgs[8];
    // IntType *nFPC = (IntType *)userArgs[9];
    // IntType **C2F = (IntType **)userArgs[10];
    // IntType vis_mode = *(IntType *)userArgs[11];
    // IntType GaussLayer = *(IntType *)userArgs[12];
    // IntType *CellLayerNo = (IntType *)userArgs[13];
    // RealFlow *u_n  = (RealFlow *)userArgs[10];
    // RealFlow *v_n  = (RealFlow *)userArgs[11];
    // RealFlow *w_n  = (RealFlow *)userArgs[12];
    // RealFlow *xfn_n_symm  = (RealFlow *)userArgs[13];
    // RealFlow *yfn_n_symm  = (RealFlow *)userArgs[14];
    // RealFlow *zfn_n_symm  = (RealFlow *)userArgs[15];
    // IntType *nodesymm   = (IntType *)userArgs[16];

    // IntType nTNode = grid->GetNTNode();
    // IntType nTCell = grid->GetNTCell();
    // IntType nTFace = grid->GetNTFace();
    // IntType nBFace = grid->GetNBFace();
    // IntType *f2c = grid->Getf2c();
    // BCRecord **bcr = grid->Getbcr();
    // IntType n = nTCell + nBFace;
    RealGeom *vol = grid->GetCellVol();
    // RealGeom *area = grid->GetFaceArea();
    // RealGeom *xfn = grid->GetXfn();
    // RealGeom *yfn = grid->GetYfn();
    // RealGeom *zfn = grid->GetZfn();
    // IntType *nNPF = grid->GetnNPF();

    IntType nsCell = treeArgs->firstCell;
    IntType neCell = treeArgs->lastCell + 1;
    // IntType nsFace = treeArgs->firstFace;
    // IntType neFace = treeArgs->lastFace + 1;
    // IntType nsNode = treeArgs->firstNode;
    // IntType neNode = treeArgs->lastNode + 1;

    IntType i, j, c1, c2, type, face, name;
    RealGeom tmpx, tmpy, tmpz;
    RealFlow qsum;
    RealFlow vn;

    // IntType *cellwallnumber = grid->GetGridQualityCellWallNumber();

    for (name = ns; name < ne; name++)
    {
        // if (vis_mode != INVISCID)
        // {
        //     for (i = nsCell; i < neCell; i++)
        //     {
        //         if (cellwallnumber[i] < 2)
        //             continue;
        //         dqdx[name][i] = 0.0;
        //         dqdy[name][i] = 0.0;
        //         dqdz[name][i] = 0.0;
        //         for (j = 0; j < nFPC[i]; j++)
        //         {
        //             face = C2F[i][j];
        //             c1 = f2c[face + face];
        //             c2 = f2c[face + face + 1];
        //             qsum = 0.5 * (q[name][c1] + q[name][c2]) * area[face];
        //             tmpx = qsum * xfn[face];
        //             tmpy = qsum * yfn[face];
        //             tmpz = qsum * zfn[face];
        //             if (i == c1)
        //             {
        //                 dqdx[name][i] += tmpx;
        //                 dqdy[name][i] += tmpy;
        //                 dqdz[name][i] += tmpz;
        //             }
        //             else if (i == c2)
        //             {
        //                 dqdx[name][i] -= tmpx;
        //                 dqdy[name][i] -= tmpy;
        //                 dqdz[name][i] -= tmpz;
        //             }
        //             else
        //             {
        //                 std::cerr << endl
        //                           << "Error in function CompGradientQ_Gauss_Node_Dist!  i is not c1 or c2! " << i << "  " << c1 << "  " << c2 << endl;
        //                 mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
        //             }
        //         }
        //     }
        // }

        // // 边界层前n层采用Gauss方法
        // if (GaussLayer > 0)
        // {
        //     for (i = nsCell; i < neCell; i++)
        //     {
        //         if ((CellLayerNo[i] == -1) || (CellLayerNo[i] >= GaussLayer))
        //             continue;
        //         dqdx[name][i] = 0.0;
        //         dqdy[name][i] = 0.0;
        //         dqdz[name][i] = 0.0;
        //         for (j = 0; j < nFPC[i]; j++)
        //         {
        //             face = C2F[i][j];
        //             c1 = f2c[face + face];
        //             c2 = f2c[face + face + 1];
        //             qsum = 0.5 * (q[name][c1] + q[name][c2]) * area[face];
        //             tmpx = qsum * xfn[face];
        //             tmpy = qsum * yfn[face];
        //             tmpz = qsum * zfn[face];
        //             if (i == c1)
        //             {
        //                 dqdx[name][i] += tmpx;
        //                 dqdy[name][i] += tmpy;
        //                 dqdz[name][i] += tmpz;
        //             }
        //             else if (i == c2)
        //             {
        //                 dqdx[name][i] -= tmpx;
        //                 dqdy[name][i] -= tmpy;
        //                 dqdz[name][i] -= tmpz;
        //             }
        //             else
        //             {
        //                 std::cerr << endl
        //                           << "Error in function CompGradientQ_Gauss_Node_Dist!  i is not c1 or c2! " << i << "  " << c1 << "  " << c2 << endl;
        //                 mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
        //             }
        //         }
        //     }
        // }
        for (i = nsCell; i < neCell; i++)
        {
            dqdx[name][i] /= vol[i];
            dqdy[name][i] /= vol[i];
            dqdz[name][i] /= vol[i];
        }
    }
    // }
}
#endif