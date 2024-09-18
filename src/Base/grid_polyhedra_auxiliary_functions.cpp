/*!
 * @file        grid_polyhedra_auxiliary_functions.cpp
 * @brief       Auxiliary functions for PolyGrid
 *
 * @author      Wisces ��wwwangqs17@163.com��
 * @version     1.0
 * @date        2023-05-21
 * @copyright   Copyright (c) 2023  PERSONAL
 *
 * @par �޸���־:
 * <table>
 * <tr><th> Date        <th> Version  <th> Author  <th> Description
 * <tr><td> 2023-05-21  <td> 1.0      <td> Wisces  <td> �����ݡ�
 * </table>
 */

// direct head files
#include "grid_polyhedra.h"

// C/C++build-in head files
// #include <cstring>
// #include <fstream>
#include <iostream>
// #include <sstream>
// #include <string>
#include <cstdlib>
#include <cmath>
// #include <cassert>
// #include <deque>
// #include <list>
// #include <set>
// #include <map>
// #include <iomanip>
using namespace std;

// user defined head files
#include "number_type.h"
#include "constant.h"
#include "memory_util.h"
#include "parallel_base.h"
#include "para_field_global.h"
// #include "zone.h"
// #include "solver_ns.h"
// #include "utility_functions.h"
// #include "algm.h"
// #include "io_log.h"
// #include "io_base_format.h"
// #include "parallel_base_functions.h"
// #include "system_base_functions.h"
// #include "grid_patch_type.h"�����ɵ� constant.h��

//!< head file relying on condition-compiling
#ifdef MF_MPICH
#include <mpi.h>
#endif

#ifdef MF_MPICH
extern int myid;
extern int numprocs;
extern int myZone;        // = myid+1
extern MPI_Comm GridComm; ///< for each grid
#endif

/**********************************************************************************
    Compute the normal vector, the area and face center on each cell face in 3D
                            ~~ Modified by ZHYB ~~
**********************************************************************************/
void FaceAreaNormalCentroid_cycle(PolyGrid *grid, RealGeom *area,
                                  RealGeom *xfn, RealGeom *yfn, RealGeom *zfn,
                                  RealGeom *xfc, RealGeom *yfc, RealGeom *zfc)
{
    IntType i, j, k, count, num_cycle, p1, p2, p3, c1;
    RealGeom x21, y21, z21, x31, y31, z31;
    RealGeom xn, yn, zn, art, xt, yt, zt;
    IntType *f2n = grid->Getf2n();
    IntType *nNPF = grid->GetnNPF();
    IntType nTFace = grid->GetNTFace();
    IntType nBFace = grid->GetNBFace();
    IntType *f2c = grid->Getf2c();
    RealGeom *xcc = grid->GetXcc();
    RealGeom *ycc = grid->GetYcc();
    RealGeom *zcc = grid->GetZcc();
    RealGeom *x = grid->GetX();
    RealGeom *y = grid->GetY();
    RealGeom *z = grid->GetZ();

    RealGeom *xfct = NULL, *yfct = NULL, *zfct = NULL;
    snew_array_1D(xfct, nTFace);
    snew_array_1D(yfct, nTFace);
    snew_array_1D(zfct, nTFace);

    // ͨ����������ȷ�������
    // num_cycle = 100;
    num_cycle = 4;
    for (k = 0; k < num_cycle; k++)
    {
        // res = 0.0;
        count = 0;
        for (i = 0; i < nTFace; i++)
        {
            area[i] = 0.0;
            xfct[i] = 0.0;
            yfct[i] = 0.0;
            zfct[i] = 0.0;

            p1 = f2n[count];
            for (j = 0; j < nNPF[i]; j++)
            {
                p2 = f2n[count++];
                if (j == nNPF[i] - 1)
                    p3 = p1;
                else
                    p3 = f2n[count];

                x21 = x[p2] - xfc[i];
                y21 = y[p2] - yfc[i];
                z21 = z[p2] - zfc[i];
                x31 = x[p3] - xfc[i];
                y31 = y[p3] - yfc[i];
                z31 = z[p3] - zfc[i];

                xn = y21 * z31 - y31 * z21;
                yn = z21 * x31 - z31 * x21;
                zn = x21 * y31 - x31 * y21;

                art = sqrt(xn * xn + yn * yn + zn * zn);

                xt = x[p2] + x[p3] + xfc[i];
                yt = y[p2] + y[p3] + yfc[i];
                zt = z[p2] + z[p3] + zfc[i];

                area[i] += art;
                xfct[i] += art * xt;
                yfct[i] += art * yt;
                zfct[i] += art * zt;
            }

            if (area[i] > TINY)
            {
                art = 1.0 / area[i] / 3.0;
                xfct[i] *= art;
                yfct[i] *= art;
                zfct[i] *= art;

                xfc[i] = xfct[i];
                yfc[i] = yfct[i];
                zfc[i] = zfct[i];
            }
        }
    } // end k

    // ���һ�ε�����ͬʱ����������ĺ��淨��
    count = 0;
    for (i = 0; i < nTFace; i++)
    {
        xfn[i] = 0.0;
        yfn[i] = 0.0;
        zfn[i] = 0.0;
        area[i] = 0.0;
        xfct[i] = 0.0;
        yfct[i] = 0.0;
        zfct[i] = 0.0;

        p1 = f2n[count];
        for (j = 0; j < nNPF[i]; j++)
        {
            p2 = f2n[count++];
            if (j == nNPF[i] - 1)
                p3 = p1;
            else
                p3 = f2n[count];

            x21 = x[p2] - xfc[i];
            y21 = y[p2] - yfc[i];
            z21 = z[p2] - zfc[i];
            x31 = x[p3] - xfc[i];
            y31 = y[p3] - yfc[i];
            z31 = z[p3] - zfc[i];

            xn = y21 * z31 - y31 * z21;
            yn = z21 * x31 - z31 * x21;
            zn = x21 * y31 - x31 * y21;

            art = sqrt(xn * xn + yn * yn + zn * zn);

            xt = x[p2] + x[p3] + xfc[i];
            yt = y[p2] + y[p3] + yfc[i];
            zt = z[p2] + z[p3] + zfc[i];

            area[i] += art;
            xfct[i] += art * xt;
            yfct[i] += art * yt;
            zfct[i] += art * zt;

            xfn[i] += xn;
            yfn[i] += yn;
            zfn[i] += zn;
        }

        if (area[i] > TINY)
        {
            art = 1.0 / area[i] / 3.0;
            xfc[i] = xfct[i] * art;
            yfc[i] = yfct[i] * art;
            zfc[i] = zfct[i] * art;
        }

        art = sqrt(xfn[i] * xfn[i] + yfn[i] * yfn[i] + zfn[i] * zfn[i]);
        area[i] = art * 0.5;
        if (art > TINY)
        {
            xfn[i] /= art;
            yfn[i] /= art;
            zfn[i] /= art;
        }
        else
        {
            c1 = f2c[i + i];

            xn = xfc[i] - xcc[c1];
            yn = yfc[i] - ycc[c1];
            zn = zfc[i] - zcc[c1];
            art = sqrt(xn * xn + yn * yn + zn * zn);

            xfn[i] = xn / art;
            yfn[i] = yn / art;
            zfn[i] = zn / art;
        }
    }

    //!< Wisces: ���� area_scale û��ʵ����;����������������������������
    // // find area scale search wall boundary faces first.
    // BCRecord **bcr = grid->Getbcr();
    // IntType type;
    // RealGeom area_scale = 0.0;
    // if (bcr != NULL)
    // {
    //     for (i = 0; i < nBFace; i++)
    //     {
    //         type = bcr[i]->GetType();
    //         if (type != WALL)
    //             continue;
    //         area_scale = MAX(area_scale, area[i]);
    //     }
    //     area_scale *= 1.0e-14;
    // }
    // else
    // {
    //     // if there are not, search all boundary faces
    //     for (i = 0; i < nBFace; i++)
    //     {
    //         area_scale = MAX(area_scale, area[i]);
    //     }
    //     area_scale *= 1.0e-16;
    // }
    // // #ifdef PARALLEL_MPI
    // //     RealGeom area_scale_glb;
    // //     MPI_Allreduce(&area_scale, &area_scale_glb, 1, MPIReal, MPI_MAX, MPI_COMM_WORLD);
    // //     area_scale = area_scale_glb;
    // // #endif
    // grid->UpdateData(&area_scale, REAL_GEOM, 1, "area_scale");
    sdel_array_1D(xfct);
    sdel_array_1D(yfct);
    sdel_array_1D(zfct);
}

/************************************************************************
  Compute face and cell centroid by simple average all node value
                       ~~ Modified by ZHYB ~~
************************************************************************/
void FaceCellCenterbyAverage(PolyGrid *grid, RealGeom *xfc, RealGeom *yfc, RealGeom *zfc, RealGeom *xcc, RealGeom *ycc, RealGeom *zcc)
{
    IntType nBFace = grid->GetNBFace();
    IntType nTFace = grid->GetNTFace();
    IntType nTCell = grid->GetNTCell();
    IntType *f2c = grid->Getf2c();
    IntType *f2n = grid->Getf2n();
    IntType *nNPF = grid->GetnNPF();
    RealGeom *x = grid->GetX();
    RealGeom *y = grid->GetY();
    RealGeom *z = grid->GetZ();

    IntType i, j, p1, node, cell, c1, c2;

    // face center
    node = 0;
    for (i = 0; i < nTFace; i++)
    {
        xfc[i] = 0.0;
        yfc[i] = 0.0;
        zfc[i] = 0.0;

        for (j = 0; j < nNPF[i]; j++)
        {
            p1 = f2n[node++];
            xfc[i] += x[p1];
            yfc[i] += y[p1];
            zfc[i] += z[p1];
        }
        xfc[i] /= nNPF[i];
        yfc[i] /= nNPF[i];
        zfc[i] /= nNPF[i];
    }

    // cell center, the average of face center
    for (i = 0; i < nTCell; i++)
    {
        xcc[i] = 0.0;
        ycc[i] = 0.0;
        zcc[i] = 0.0;
    }

    IntType *nface = NULL;
    snew_array_1D(nface, nTCell);
    for (i = 0; i < nTCell; i++)
        nface[i] = 0;
    cell = 0;
    for (i = 0; i < nBFace; i++)
    {
        c1 = f2c[cell++];
        cell++;

        xcc[c1] += xfc[i];
        ycc[c1] += yfc[i];
        zcc[c1] += zfc[i];
        nface[c1]++;
    }
    for (i = nBFace; i < nTFace; i++)
    {
        c1 = f2c[cell++];
        c2 = f2c[cell++];

        xcc[c1] += xfc[i];
        ycc[c1] += yfc[i];
        zcc[c1] += zfc[i];
        xcc[c2] += xfc[i];
        ycc[c2] += yfc[i];
        zcc[c2] += zfc[i];
        nface[c1]++;
        nface[c2]++;
    }

    for (i = 0; i < nTCell; i++)
    {
        xcc[i] /= nface[i];
        ycc[i] /= nface[i];
        zcc[i] /= nface[i];
    }

    sdel_array_1D(nface);
}

/************************************************************************
              Compute Cell Center and Cell Volume in 3D
                       ~~ Modified by ZHYB ~~
************************************************************************/
void CellVolCentroid(PolyGrid *grid, RealGeom *vol, RealGeom *xcc, RealGeom *ycc, RealGeom *zcc)
{
    IntType i, j, cell, node, c1, c2, p1, p2, p3;
    IntType nBFace = grid->GetNBFace();
    IntType nTFace = grid->GetNTFace();
    IntType nTCell = grid->GetNTCell();
    IntType *f2c = grid->Getf2c();
    IntType *f2n = grid->Getf2n();
    IntType *nNPF = grid->GetnNPF();
    RealGeom *xfc = grid->GetXfc();
    RealGeom *yfc = grid->GetYfc();
    RealGeom *zfc = grid->GetZfc();
    RealGeom *xfn = grid->GetXfn();
    RealGeom *yfn = grid->GetYfn();
    RealGeom *zfn = grid->GetZfn();
    RealGeom *area = grid->GetFaceArea();
    RealGeom *x = grid->GetX();
    RealGeom *y = grid->GetY();
    RealGeom *z = grid->GetZ();
    RealGeom x21, y21, z21, x31, y31, z31, x41, y41, z41, nx, ny, nz;
    RealGeom volt, dot, xt, yt, zt;
    RealGeom minvol, maxvol, tmp;

    RealGeom *xcct = NULL, *ycct = NULL, *zcct = NULL;
    snew_array_1D(xcct, nTCell);
    snew_array_1D(ycct, nTCell);
    snew_array_1D(zcct, nTCell);
    for (i = 0; i < nTCell; i++)
    {
        xcct[i] = 0.0;
        ycct[i] = 0.0;
        zcct[i] = 0.0;
        vol[i] = 0.0;
    }

    cell = 0;
    node = 0;
    for (i = 0; i < nBFace; i++)
    {
        c1 = f2c[cell++];
        cell++;

        p1 = f2n[node];
        for (j = 0; j < nNPF[i]; j++)
        {
            p2 = f2n[node++];
            if (j == nNPF[i] - 1)
                p3 = p1;
            else
                p3 = f2n[node];

            x21 = x[p2] - xfc[i];
            y21 = y[p2] - yfc[i];
            z21 = z[p2] - zfc[i];
            x31 = x[p3] - xfc[i];
            y31 = y[p3] - yfc[i];
            z31 = z[p3] - zfc[i];
            nx = y21 * z31 - y31 * z21;
            ny = z21 * x31 - z31 * x21;
            nz = x21 * y31 - x31 * y21;

            x41 = xcc[c1] - xfc[i];
            y41 = ycc[c1] - yfc[i];
            z41 = zcc[c1] - zfc[i];

            xt = x[p2] + x[p3] + xfc[i] + xcc[c1];
            yt = y[p2] + y[p3] + yfc[i] + ycc[c1];
            zt = z[p2] + z[p3] + zfc[i] + zcc[c1];

            volt = nx * x41 + ny * y41 + nz * z41;
            volt = -volt;
            xcct[c1] += volt * xt;
            ycct[c1] += volt * yt;
            zcct[c1] += volt * zt;
            vol[c1] += volt;
        }
    }
    for (i = nBFace; i < nTFace; i++)
    {
        c1 = f2c[cell++];
        c2 = f2c[cell++];

        p1 = f2n[node];
        for (j = 0; j < nNPF[i]; j++)
        {
            p2 = f2n[node++];
            if (j == nNPF[i] - 1)
                p3 = p1;
            else
                p3 = f2n[node];

            x21 = x[p2] - xfc[i];
            y21 = y[p2] - yfc[i];
            z21 = z[p2] - zfc[i];
            x31 = x[p3] - xfc[i];
            y31 = y[p3] - yfc[i];
            z31 = z[p3] - zfc[i];
            nx = y21 * z31 - y31 * z21;
            ny = z21 * x31 - z31 * x21;
            nz = x21 * y31 - x31 * y21;

            x41 = xcc[c1] - xfc[i];
            y41 = ycc[c1] - yfc[i];
            z41 = zcc[c1] - zfc[i];
            xt = x[p2] + x[p3] + xfc[i] + xcc[c1];
            yt = y[p2] + y[p3] + yfc[i] + ycc[c1];
            zt = z[p2] + z[p3] + zfc[i] + zcc[c1];
            volt = nx * x41 + ny * y41 + nz * z41;
            volt = -volt;
            xcct[c1] += volt * xt;
            ycct[c1] += volt * yt;
            zcct[c1] += volt * zt;
            vol[c1] += volt;

            x41 = xcc[c2] - xfc[i];
            y41 = ycc[c2] - yfc[i];
            z41 = zcc[c2] - zfc[i];
            xt = x[p2] + x[p3] + xfc[i] + xcc[c2];
            yt = y[p2] + y[p3] + yfc[i] + ycc[c2];
            zt = z[p2] + z[p3] + zfc[i] + zcc[c2];
            volt = nx * x41 + ny * y41 + nz * z41;
            xcct[c2] += volt * xt;
            ycct[c2] += volt * yt;
            zcct[c2] += volt * zt;
            vol[c2] += volt;
        }
    }

    // mflog::log.set_all_processors_out();

    cell = 0;
    minvol = BIG;
    maxvol = -BIG;
    for (i = 0; i < nTCell; i++)
    {
        tmp = 1.0 / (4.0 * vol[i] + TINY);
        xcct[i] *= tmp;
        ycct[i] *= tmp;
        zcct[i] *= tmp;
        vol[i] /= 6.0;
        minvol = MIN(minvol, vol[i]);
        maxvol = MAX(maxvol, vol[i]);
        if (vol[i] < TINY)
        {
            cell++;
            cout << endl
                 << "cell volume is:"
                 << " " << vol[i] << endl;
            cout << "cell center is:"
                 << " " << xcct[i] << " " << ycct[i] << " " << zcct[i] << endl;
        }
    }

#ifdef MF_MPICH
    Parallel::parallel_sum(cell, MPI_COMM_WORLD);
    Parallel::parallel_min_max(minvol, maxvol, MPI_COMM_WORLD);
#endif

// mflog::log.set_one_processor_out();
#ifdef MF_DEBUG
#ifdef MF_MPICH
    if (myid == 0)
    {
#endif
        cout << "Min and max volumes are "
             << "\t" << minvol << "\t" << maxvol << endl;
#ifdef MF_MPICH
    }
#endif
#endif

    if (cell > 0)
    {
        // mflog::log << std::endl
        //            << "There are " << cell << " cell's vol too small! This must correct!" << std::endl;
        // mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
        cerr << endl
             << "There are " << cell << " cell's vol too small! This must correct!" << endl;
        exit(1);
    }

    // ֻ�������������ȷһ������ĲŸ��£��������Ĳ��ü�ƽ��
    IntType *mark = NULL;
    snew_array_1D(mark, nTCell);
    for (i = 0; i < nTCell; i++)
        mark[i] = 1;
    for (i = 0; i < nBFace; i++)
    {
        c1 = f2c[i + i];

        x21 = xcct[c1] - xfc[i];
        y21 = ycct[c1] - yfc[i];
        z21 = zcct[c1] - zfc[i];
        dot = x21 * xfn[i] + y21 * yfn[i] + z21 * zfn[i];
        if (dot > 0.0)
            mark[c1] = 0;
    }
    for (i = nBFace; i < nTFace; i++)
    {
        c1 = f2c[i + i];
        c2 = f2c[i + i + 1];

        x21 = xcct[c1] - xfc[i];
        y21 = ycct[c1] - yfc[i];
        z21 = zcct[c1] - zfc[i];
        dot = x21 * xfn[i] + y21 * yfn[i] + z21 * zfn[i];
        if (dot > 0.0)
            mark[c1] = 0;

        x21 = xcct[c2] - xfc[i];
        y21 = ycct[c2] - yfc[i];
        z21 = zcct[c2] - zfc[i];
        dot = x21 * xfn[i] + y21 * yfn[i] + z21 * zfn[i];
        if (dot < 0.0)
            mark[c2] = 0;
    }

    for (i = 0; i < nTCell; i++)
    {
        if (mark[i] == 1)
        {
            xcc[i] = xcct[i];
            ycc[i] = ycct[i];
            zcc[i] = zcct[i];
        }
    }

    // For ghost cells
    cell = 0;
    for (i = 0; i < nBFace; i++)
    {
        // ���б߽��ں���ĳ����е�����ֵ
        c1 = f2c[cell++];
        cell++;
        c2 = i + nTCell; // because no num. for ghost cell in some condition
        if (area[i] > TINY)
        {
            tmp = 2. * ((xcc[c1] - xfc[i]) * xfn[i] + (ycc[c1] - yfc[i]) * yfn[i] + (zcc[c1] - zfc[i]) * zfn[i]);
            xcc[c2] = xcc[c1] - xfn[i] * tmp;
            ycc[c2] = ycc[c1] - yfn[i] * tmp;
            zcc[c2] = zcc[c1] - zfn[i] * tmp;
        }
        else
        {
            // Degenerated faces
            xcc[c2] = -xcc[c1] + 2. * xfc[i];
            ycc[c2] = -ycc[c1] + 2. * yfc[i];
            zcc[c2] = -zcc[c1] + 2. * zfc[i];
        }

        vol[c2] = vol[c1];
    }
    sdel_array_1D(xcct);
    sdel_array_1D(ycct);
    sdel_array_1D(zcct);
    sdel_array_1D(mark);
}

/************************************************************************
              Correct face normal vector for face of TINY AREA
                       ~~ Builded by ZHYB ~~
************************************************************************/
void CorrectFaceNormal(PolyGrid *grid, RealGeom *xfn, RealGeom *yfn, RealGeom *zfn)
{
    IntType i, c1, p1, node;
    RealGeom x21, y21, z21, x31, y31, z31, dot, nn, ratio;
    IntType nTFace = grid->GetNTFace();
    IntType *f2n = grid->Getf2n();
    IntType *nNPF = grid->GetnNPF();
    IntType *f2c = grid->Getf2c();
    RealGeom *x = grid->GetX();
    RealGeom *y = grid->GetY();
    RealGeom *z = grid->GetZ();
    RealGeom *xfc = grid->GetXfc();
    RealGeom *yfc = grid->GetYfc();
    RealGeom *zfc = grid->GetZfc();
    RealGeom *xcc = grid->GetXcc();
    RealGeom *ycc = grid->GetYcc();
    RealGeom *zcc = grid->GetZcc();
    RealGeom *area = grid->GetFaceArea();

    node = 0;
    for (i = 0; i < nTFace; i++)
    {
        p1 = f2n[node];
        node += nNPF[i];
        if (area[i] < TINY)
        {
            c1 = f2c[i + i];

            x21 = xcc[c1] - xfc[i];
            y21 = ycc[c1] - yfc[i];
            z21 = zcc[c1] - zfc[i];
            x31 = x[p1] - xfc[i];
            y31 = y[p1] - yfc[i];
            z31 = z[p1] - zfc[i];

            dot = x21 * x31 + y21 * y31 + z21 * z31;
            nn = x31 * x31 + y31 * y31 + z31 * z31;

            ratio = 0.0;
            if (nn > TINY)
                ratio = dot / nn;

            x21 -= ratio * x31;
            y21 -= ratio * y31;
            z21 -= ratio * z31;

            nn = sqrt(x21 * x21 + y21 * y21 + z21 * z21);
            xfn[i] = -x21 / nn;
            yfn[i] = -y21 / nn;
            zfn[i] = -z21 / nn;
        }
    }
}

/************************************************************************
 Correct cell centroid if it is too close to wall through moving cell
 centroid in the direct of normal vector
 by ZHYB
 2019-11-01
************************************************************************/
void CorrectCellCentroid(PolyGrid *grid, RealGeom *xcc, RealGeom *ycc, RealGeom *zcc,
                         RealGeom *xfc, RealGeom *yfc, RealGeom *zfc,
                         RealGeom *xfn, RealGeom *yfn, RealGeom *zfn)
{
    // �޸����ģ������ĵ���ľ���С��minDn�����������ģ�ʹ�����ٴﵽminDn
    // ���ĵ���ľ����С�������ʱ�䲽����С����������������������

    IntType nBFace = grid->GetNBFace();
    IntType nTFace = grid->GetNTFace();
    IntType *f2c = grid->Getf2c();

    IntType iexp = 16;
    // grid->GetData(&iexp, INT, 1, "iexp");
    const RealGeom scaleXYZcc = pow(10.0, -iexp + 1); // ��ֹ�ƶ�����̫С���ü��㾫��������
    const RealGeom minDn = 1.0e-10;                   // �߶ȣ�С�ڸó߶ȵ���Ҫ����������

    IntType i, c1, c2, sign;
    RealGeom x21, y21, z21, dot, dx, dy, dz;

    IntType count = 0;
    for (i = 0; i < nBFace; ++i)
    {
        c1 = f2c[i + i];

        x21 = xcc[c1] - xfc[i];
        y21 = ycc[c1] - yfc[i];
        z21 = zcc[c1] - zfc[i];
        dot = x21 * xfn[i] + y21 * yfn[i] + z21 * zfn[i];

        if (fabs(dot) < minDn)
        {
#ifdef MF_DEBUG
            count++;
            cout << endl
                 << "dot=" << dot
                 << endl
                 << "xcc=" << xcc[c1] << "  " << ycc[c1] << "  " << zcc[c1]
                 << endl
                 << "xfc=" << xfc[i] << "  " << yfc[i] << "  " << zfc[i]
                 << endl
                 << "xfn=" << xfn[i] << "  " << yfn[i] << "  " << zfn[i]
                 << endl;
#endif

            dx = minDn * xfn[i];
            dy = minDn * yfn[i];
            dz = minDn * zfn[i];

            dx = AbsMaxSignFirst(dx, scaleXYZcc * xcc[c1]);
            dy = AbsMaxSignFirst(dy, scaleXYZcc * ycc[c1]);
            dz = AbsMaxSignFirst(dz, scaleXYZcc * zcc[c1]);

            sign = (dot > 0.0) ? 1 : -1;
            xcc[c1] += sign * dx;
            ycc[c1] += sign * dy;
            zcc[c1] += sign * dz;
        }
    }

    for (i = nBFace; i < nTFace; i++)
    {
        c1 = f2c[i + i];
        c2 = f2c[i + i + 1];

        x21 = xcc[c1] - xfc[i];
        y21 = ycc[c1] - yfc[i];
        z21 = zcc[c1] - zfc[i];
        dot = x21 * xfn[i] + y21 * yfn[i] + z21 * zfn[i];
        if (fabs(dot) < minDn)
        {
#ifdef MF_DEBUG
            count++;
            cout << endl
                 << "dot=" << dot
                 << endl
                 << "xcc=" << xcc[c1] << "  " << ycc[c1] << "  " << zcc[c1]
                 << endl
                 << "xfc=" << xfc[i] << "  " << yfc[i] << "  " << zfc[i]
                 << endl
                 << "xfn=" << xfn[i] << "  " << yfn[i] << "  " << zfn[i]
                 << endl;
#endif
            dx = minDn * xfn[i];
            dy = minDn * yfn[i];
            dz = minDn * zfn[i];

            dx = AbsMaxSignFirst(dx, scaleXYZcc * xcc[c1]);
            dy = AbsMaxSignFirst(dy, scaleXYZcc * ycc[c1]);
            dz = AbsMaxSignFirst(dz, scaleXYZcc * zcc[c1]);

            sign = (dot > 0.0) ? 1 : -1;
            xcc[c1] += sign * dx;
            ycc[c1] += sign * dy;
            zcc[c1] += sign * dz;
        }

        x21 = xcc[c2] - xfc[i];
        y21 = ycc[c2] - yfc[i];
        z21 = zcc[c2] - zfc[i];
        dot = x21 * xfn[i] + y21 * yfn[i] + z21 * zfn[i];
        if (fabs(dot) < minDn)
        {
#ifdef MF_DEBUG
            count++;
            cout << endl
                 << "dot=" << dot
                 << endl
                 << "xcc=" << xcc[c2] << "  " << ycc[c2] << "  " << zcc[c2]
                 << endl
                 << "xfc=" << xfc[i] << "  " << yfc[i] << "  " << zfc[i]
                 << endl
                 << "xfn=" << xfn[i] << "  " << yfn[i] << "  " << zfn[i]
                 << endl;
#endif
            dx = minDn * xfn[i];
            dy = minDn * yfn[i];
            dz = minDn * zfn[i];

            dx = AbsMaxSignFirst(dx, scaleXYZcc * xcc[c2]);
            dy = AbsMaxSignFirst(dy, scaleXYZcc * ycc[c2]);
            dz = AbsMaxSignFirst(dz, scaleXYZcc * zcc[c2]);

            sign = (dot > 0.0) ? 1 : -1;
            xcc[c2] += sign * dx;
            ycc[c2] += sign * dy;
            zcc[c2] += sign * dz;
        }
    }
}

/************************************************************************
  ����ֵȡ������������ֵ�еĴ�ֵ������ȡ��һ�������ķ���
************************************************************************/
RealGeom AbsMaxSignFirst(RealGeom a, RealGeom b)
{
    IntType sign;
    sign = (a > 0) ? 1 : -1;

    RealGeom c;
    c = MAX(fabs(a), fabs(b));
    c *= sign;

    return c;
}

/*!
 * @brief       ʹ�þ����Ȩ�ķ�������ڵ��ϵ�Ȩ���ӣ�������ڵ������ֵ
 * @param       grid
 *
 * @author      Wisces ��wwwangqs17@163.com��
 * @date        2023-07-09
 */
void ComputeWeight3D_Node(PolyGrid *grid)
{
    IntType i, j, type, p1;
    IntType nTNode = grid->GetNTNode();
    IntType nBFace = grid->GetNBFace();
    IntType nTCell = grid->GetNTCell();

    IntType *nNPF = grid->GetnNPF();
    // IntType *nNPC = grid->GetnNPC();
    IntType *nNPC = CalnNPC(grid);
    IntType *nCPN = CalnCPN(grid);
    IntType **N2C = CalN2C(grid);
    // IntType **C2N = grid->GetC2N();
    IntType **C2N = CalC2N(grid);
    IntType **F2N = CalF2N(grid);
    RealGeom *x = grid->GetX();
    RealGeom *y = grid->GetY();
    RealGeom *z = grid->GetZ();
    RealGeom *xcc = grid->GetXcc();
    RealGeom *ycc = grid->GetYcc();
    RealGeom *zcc = grid->GetZcc();
    RealGeom *xfc = grid->GetXfc();
    RealGeom *yfc = grid->GetYfc();
    RealGeom *zfc = grid->GetZfc();
    BCRecord **bcr = grid->Getbcr();
    RealGeom dx, dy, dz, wr;

    // if (grid->GetWeightNodeDist() == NULL)
    // {
    //     RealGeom *WeightNodeDist = NULL;
    //     snew_array_1D(WeightNodeDist, nTNode);
    //     grid->SetWeightNodeDist(WeightNodeDist);
    // }
    // if (grid->GetNodeType() == NULL)
    // {
    //     IntType *Nmark = NULL;
    //     snew_array_1D(Nmark, nTNode);
    //     grid->SetNodeType(Nmark);
    // }
    // if (grid->GetWeightNodeC2N() == NULL)
    // {
    //     RealGeom **WeightNodeC2N = NULL;
    //     snew_array_2D(WeightNodeC2N, nTCell, nNPC, true);
    //     grid->SetWeightNodeC2N(WeightNodeC2N);
    // }
    // if (grid->GetWeightNodeN2C() == NULL)
    // {
    //     RealGeom **WeightNodeN2C = NULL;
    //     snew_array_2D(WeightNodeN2C, nTNode, nCPN, true);
    //     grid->SetWeightNodeN2C(WeightNodeN2C);
    // }
    // RealGeom *WeightNode = grid->GetWeightNodeDist();
    // RealGeom **WeightNodeC2N = grid->GetWeightNodeC2N();
    // IntType *Nmark = grid->GetNodeType();
    // RealGeom **WeightNodeN2C = grid->GetWeightNodeN2C();

    //!< �����ڴ�
    RealGeom *WeightNode = NULL, **WeightNodeC2N = NULL, **WeightNodeN2C = NULL, **WeightNodeBFace2C = NULL;
    IntType *Nmark = NULL;
    snew_array_1D(WeightNode, nTNode);
    snew_array_1D(Nmark, nTNode);
    snew_array_2D(WeightNodeC2N, nTCell, nNPC, true);
    snew_array_2D(WeightNodeN2C, nTNode, nCPN, true);
    snew_array_2D(WeightNodeBFace2C, nBFace, nNPF, true);

    grid->SetWeightNodeDist(WeightNode);
    grid->SetNodeType(Nmark);
    grid->SetWeightNodeC2N(WeightNodeC2N);
    grid->SetWeightNodeN2C(WeightNodeN2C);
    grid->SetWeightNodeBFace2C(WeightNodeBFace2C);

    //!< ��ʼ��
    for (i = 0; i < nTNode; i++)
    {
        Nmark[i] = 0;
        WeightNode[i] = 0.;
    }
    for (i = 0; i < nTCell; i++)
    {
        for (j = 0; j < nNPC[i]; j++)
        {
            WeightNodeC2N[i][j] = 0.0;
        }
    }
    for (i = 0; i < nTNode; i++)
    {
        for (j = 0; j < nCPN[i]; j++)
        {
            WeightNodeN2C[i][j] = 0.0;
        }
    }
    for (i = 0; i < nBFace; i++)
    {
        for (j = 0; j < nNPF[i]; j++)
        {
            WeightNodeBFace2C[i][j] = 0.0;
        }
    }

    //!< ��������߽���������ֵ����������߽���Ȩϵ��
    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        if (type != WALL)
            continue;
        for (j = 0; j < nNPF[i]; j++)
        {
            p1 = F2N[i][j];
            dx = x[p1] - xfc[i];
            dy = y[p1] - yfc[i];
            dz = z[p1] - zfc[i];
            wr = dx * dx + dy * dy + dz * dz;
            wr = sqrt(wr);
            wr = 1. / wr;
            WeightNodeBFace2C[i][j] = wr;
            WeightNode[p1] += wr;
            Nmark[p1] = WALL;
        }
    }
#ifdef MF_MPICH
    grid->CommInternodeDataMPI2(Nmark);
#endif

    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        if (type != FAR_FIELD && type != FAR_FIELD_MOD)
            continue;
        for (j = 0; j < nNPF[i]; j++)
        {
            p1 = F2N[i][j];
            if (Nmark[p1] == WALL)
                continue;
            dx = x[p1] - xfc[i];
            dy = y[p1] - yfc[i];
            dz = z[p1] - zfc[i];
            wr = dx * dx + dy * dy + dz * dz;
            wr = sqrt(wr);
            wr = 1. / wr;
            WeightNodeBFace2C[i][j] = wr;
            WeightNode[p1] += wr;
            Nmark[p1] = FAR_FIELD;
        }
    }
#ifdef MF_MPICH
    grid->CommInternodeDataMPI2(Nmark);
#endif

    //!< for other boundary condition(symmetry not include), e.g. penliu (zhyb)
    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        if (type == WALL || type == SYMM || type == FAR_FIELD || type == FAR_FIELD_MOD || type == INTERFACE)
            continue;
        for (j = 0; j < nNPF[i]; j++)
        {
            p1 = F2N[i][j];
            if (Nmark[p1] == WALL || Nmark[p1] == FAR_FIELD)
                continue;
            dx = x[p1] - xfc[i];
            dy = y[p1] - yfc[i];
            dz = z[p1] - zfc[i];
            wr = dx * dx + dy * dy + dz * dz;
            wr = sqrt(wr);
            wr = 1. / wr;
            WeightNodeBFace2C[i][j] = wr;
            WeightNode[p1] += wr;
            Nmark[p1] = 101;
        }
    }

#ifdef MF_MPICH
    grid->CommInternodeDataMPI2(Nmark);
#endif

    //!< ���������������ֵ��ʹ�����������ڵĿ���������ֵ
    for (i = 0; i < nTCell; i++)
    {
        for (j = 0; j < nNPC[i]; j++)
        {
            p1 = C2N[i][j];
            if (Nmark[p1] != 0)
                continue;
            dx = x[p1] - xcc[i];
            dy = y[p1] - ycc[i];
            dz = z[p1] - zcc[i];
            wr = dx * dx + dy * dy + dz * dz;
            wr = sqrt(wr);
            wr = 1. / wr;
            WeightNodeC2N[i][j] = wr;
            WeightNode[p1] += wr;
        }
    }

    //!< Calculate WeightNodeN2C from WeightNodeC2N:
    for (i = 0; i < nTNode; i++)
    {
        if (Nmark[i] != 0)
            continue;
        for (IntType j = 0; j < nCPN[i]; j++)
        {
            IntType cellx = N2C[i][j];
            IntType k;
            for (k = 0; k < nNPC[cellx]; k++)
            {
                IntType nodex = C2N[cellx][k];
                if (nodex == i)
                {
                    break;
                }
            }
            WeightNodeN2C[i][j] = WeightNodeC2N[cellx][k];
        }
    }
//!< ���ݲ��б߽��ļ�Ȩֵ
#ifdef MF_MPICH
    grid->CommInternodeDataMPI(WeightNode);
#endif
}

/**
//!< Grid closure checking
void ClosureCheck(PolyGrid *grid, RealGeom *xfn, RealGeom *area)
{
    IntType nTCell = grid->GetNTCell(), *f2c = grid->Getf2c(),
            nTFace = grid->GetNTFace();
    IntType i, count, c1, c2;

    RealGeom *sum = NULL;
    RealGeom *asum = NULL;
    snew_array_1D(sum, nTCell);
    snew_array_1D(asum, nTCell);
    for (i = 0; i < nTCell; i++)
        sum[i] = 0.;
    for (i = 0; i < nTCell; i++)
        asum[i] = 0.;

    count = 0;
    for (i = 0; i < nTFace; i++)
    {
        c1 = f2c[count++];
        c2 = f2c[count++];

        sum[c1] += xfn[i] * area[i];
        asum[c1] += fabs(xfn[i] * area[i]);
        if (c2 < 0 || c2 >= nTCell)
            continue;
        sum[c2] -= xfn[i] * area[i];
        asum[c2] += fabs(xfn[i] * area[i]);
    }

    // total = 0.;
    // mflog::log.set_all_processors_out();
    for (i = 0; i < nTCell; i++)
    {
        // total += fabs(sum[i]);
        if (fabs(sum[i]) / asum[i] > 1.e-2)
            // mflog::log << "Area sum for cell " << i + 1 << " is " << sum[i] << " " << fabs(sum[i]) / asum[i] << std::endl;
            cout << "Area sum for cell " << i + 1 << " is " << sum[i] << " " << fabs(sum[i]) / asum[i] << std::endl;
    }

    sdel_array_1D(sum);
    sdel_array_1D(asum);
}

*/

/*!
 * @brief       Calculate nCPC(number of neighboring cells in each cell without including itself)
 * @param       grid
 * @return      IntType*
 * @attention   nCPC does not include cell itself and ghost cells on physical boundary. nCPC do include ghost cells on parallel boundary.
 *
 * @author      Wisces ��wwwangqs17@163.com��
 * @date        2023-07-06
 */
IntType *CalnCPC(PolyGrid *grid)
{
    IntType *nCPC = grid->GetnCPC();
    //!< IF nCPC has already existed
    if (nCPC)
        return nCPC;

    IntType i, c1, c2, count, type;
    IntType nBFace = grid->GetNBFace();
    IntType nTFace = grid->GetNTFace();
    IntType nTCell = grid->GetNTCell();
    IntType *f2c = grid->Getf2c();
    BCRecord **bcr = grid->Getbcr();

    //!< Allocate memories for number of cells per cell
    snew_array_1D(nCPC, nTCell);
    // assert(nCPC != 0);

    //!< set nCPC to be 0 without including itself
    for (i = 0; i < nTCell; i++)
        nCPC[i] = 0;

    //!< If boundary is an INTERFACE, need to count ghost cell
    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        if (type == INTERFACE)
        {
            count = 2 * i;
            c1 = f2c[count];
            nCPC[c1]++;
        }
    }

    //!< Interior faces
    count = 2 * nBFace;
    for (i = nBFace; i < nTFace; i++)
    {
        c1 = f2c[count++];
        c2 = f2c[count++];
        nCPC[c1]++;
        nCPC[c2]++;
    }

    //!< Attach nCPC to the grid and return it
    grid->SetnCPC(nCPC);
    return nCPC;
}

/*!
 * @brief       calculate c2c (Cell to cell connection)
 * @param       grid
 * @return      IntType**
 * @attention   c2c does not include cell itself and ghost cells on physical boundary. c2c do include ghost cells on parallel boundary.
 *
 * @author      Wisces ��wwwangqs17@163.com��
 * @date        2023-07-06
 */
IntType **CalC2C(PolyGrid *grid)
{
    IntType **c2c = grid->Getc2c();
    // IF c2c has already existed
    if (c2c)
        return c2c;

    IntType i, c1, c2, count, type;
    IntType nBFace = grid->GetNBFace();
    IntType nTFace = grid->GetNTFace();
    IntType nTCell = grid->GetNTCell();
    IntType *f2c = grid->Getf2c();
    BCRecord **bcr = grid->Getbcr();
    IntType *nCPC = grid->GetnCPC();

    // Check if nCPC is available
    if (!nCPC)
        nCPC = CalnCPC(grid);

    // Allocate memories for cell to cell connection
    snew_array_2D(c2c, nTCell, nCPC, true);
    // Need to reset nCPC to 0 and recover it later
    for (i = 0; i < nTCell; i++)
        nCPC[i] = 0;

    // If boundary is an INTERFACE, need to count ghost cell
    // Note: Symmetry???
    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        if (type == INTERFACE)
        {
            count = 2 * i;
            c1 = f2c[count++];
            c2 = f2c[count];
            c2c[c1][nCPC[c1]++] = c2;
        }
    }

    // Interior faces
    count = 2 * nBFace;
    for (i = nBFace; i < nTFace; i++)
    {
        c1 = f2c[count++];
        c2 = f2c[count++];
        c2c[c1][nCPC[c1]++] = c2;
        c2c[c2][nCPC[c2]++] = c1;
    }

    // Attach c2c to the grid and return it
    grid->Setc2c(c2c);
    return c2c;
}

/**
    //  Calculate fcptr (There are two cells for each cell face. Which cell
    //  number is a cell in its neighboring cell).

void CalCNNCF(PolyGrid *grid)
{
    IntType nTFace = grid->GetNTFace();
    // IF fcptr has already existed
    IntType *fcptr;
    // fcptr = (IntType *)grid->GetDataPtr(INT, 2 * nTFace, "fcptr");
    if (fcptr)
        return;

    IntType i, c1, c2, count, type;
    IntType nBFace = grid->GetNBFace();
    IntType nTCell = grid->GetNTCell();
    IntType *f2c = grid->Getf2c();
    BCRecord **bcr = grid->Getbcr();
    IntType *nCPC = grid->GetnCPC();

    // Allocate memories for number of cells per cell
    snew_array_1D(fcptr, 2 * nTFace);
    assert(fcptr != 0);
    for (i = 0; i < 2 * nTFace; i++)
        fcptr[i] = -1;

    if (!nCPC)
        nCPC = CalnCPC(grid);

    // Need to reset nCPC to zero and recover it later
    for (i = 0; i < nTCell; i++)
        nCPC[i] = 0;

    // If boundary is an INTERFACE, need to count ghost cell
    // Note: Symmetry???
    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        if (type == INTERFACE)
        {
            count = 2 * i;
            c1 = f2c[count];
            nCPC[c1]++;
            fcptr[count] = nCPC[c1];
        }
    }

    // Interior faces
    count = 2 * nBFace;
    for (i = nBFace; i < nTFace; i++)
    {
        c1 = f2c[count];
        nCPC[c1]++;
        fcptr[count++] = nCPC[c1];

        c2 = f2c[count];
        nCPC[c2]++;
        fcptr[count++] = nCPC[c2];
    }

    // Attach fcptr to the grid
    // grid->UpdateDataPtr(fcptr, INT, 2 * nTFace, "fcptr");
}
*/

/*!
 * @brief       Calculate nFPC (number of faces in each real cell). The real cell is "physical" cell, ghost cells not included
 * @param       grid
 * @return      IntType*
 * @note        This function is almost the same as CalnCPC. But we can't use it, because nCPC doesn't count boundary cells except for interfaces.
 *
 * @author      Wisces ��wwwangqs17@163.com��
 * @date
 */
IntType *CalnFPC(PolyGrid *grid)
{
    IntType *nFPC = (IntType *)grid->GetnFPC();
    //!< if nFPC has already existed
    if (nFPC)
        return nFPC;

    IntType i, c1, c2, count;
    IntType nTCell = grid->GetNTCell();
    IntType nBFace = grid->GetNBFace();
    IntType nTFace = grid->GetNTFace();
    IntType *f2c = grid->Getf2c();

    //!< Allocate memories for number of faces per cell
    // nFPC = NULL;
    snew_array_1D(nFPC, nTCell);
    // assert(nFPC != 0);

    //!< set nFPC to zero
    for (i = 0; i < nTCell; i++)
        nFPC[i] = 0;

    //!< Boundary faces
    for (i = 0; i < nBFace; i++)
    {
        c1 = f2c[2 * i];
        nFPC[c1]++;
    }
    //!< Interior faces
    count = 2 * nBFace;
    for (i = nBFace; i < nTFace; i++)
    {
        c1 = f2c[count++];
        c2 = f2c[count++];
        nFPC[c1]++;
        nFPC[c2]++;
    }

    //!< Attach nFPC to the grid and return it
    grid->SetnFPC(nFPC);
    return nFPC;
}

/*!
 * @brief       Calculate C2F (real cell to face connection). The real cell is "physical" cell, ghost cells not included
 * @param       grid
 * @return      IntType**
 * @note        This function is similar to CalC2C.
 *
 * @author      Wisces ��wwwangqs17@163.com��
 * @date
 */
IntType **CalC2F(PolyGrid *grid)
{
    IntType nTCell = grid->GetNTCell();
    IntType **C2F = (IntType **)grid->GetC2F();
    //!< if C2F has already existed
    if (C2F)
        return C2F;

    IntType i, c1, c2, count;
    IntType nBFace = grid->GetNBFace();
    IntType nTFace = grid->GetNTFace();
    IntType *f2c = grid->Getf2c();
    IntType *nFPC = (IntType *)grid->GetnFPC();

    //!< Check if nFPC is available
    if (!nFPC)
        nFPC = CalnFPC(grid);

    //!< Allocate memories for cell to face connection
    // C2F = NULL;
    snew_array_2D(C2F, nTCell, nFPC, true);

    //!< Need to reset nFPC to 0 and recover it later
    for (i = 0; i < nTCell; i++)
        nFPC[i] = 0;

    //!< Boundary faces
    for (i = 0; i < nBFace; i++)
    {
        c1 = f2c[2 * i];
        C2F[c1][nFPC[c1]++] = i;
    }
    //!< Interior faces
    count = 2 * nBFace;
    for (i = nBFace; i < nTFace; i++)
    {
        c1 = f2c[count++];
        c2 = f2c[count++];
        C2F[c1][nFPC[c1]++] = i;
        C2F[c2][nFPC[c2]++] = i;
    }

    //!< Attach C2F to the grid and return it
    grid->SetC2F(C2F);
    return C2F;
}

/*!
 * @brief       Calculate F2N (Face to Node connection).
 * @param       grid
 * @return      IntType**
 * @note        F2N is very special which is different to C2N or C2C, et al. Because f2n is the essential data that cannot be absent, so we can construct F2N as a reference to f2n. As a result, we can reduce memory use.
 *
 * @author      Wisces ��wwwangqs17@163.com��
 * @date        2023-07-06
 */
IntType **CalF2N(PolyGrid *grid)
{
    IntType **F2N = grid->GetF2N();

    // IF F2N has already existed
    if (F2N != NULL)
        return F2N;

    IntType nTFace = grid->GetNTFace();
    IntType *nNPF = grid->GetnNPF();
    IntType *f2n = grid->Getf2n();

    // Check if nNPF is available
    // assert(nNPF != 0);

    // Allocate memories for cell to cell connection
    snew_array_1D(F2N, nTFace);
    F2N[0] = f2n;
    for (IntType i = 1; i < nTFace; ++i)
    {
        F2N[i] = &(F2N[i - 1][nNPF[i - 1]]);
    }

    // Attach F2N to the grid and return it
    grid->SetF2N(F2N);
    return F2N;
}

/*!
 * @brief       Calculate C2N (Cell to node connection).
 * @param       grid
 * @return      IntType**
 * @note        The node order for regular cell definition complies the CGNS standard
 *
 * @author      Wisces ��wwwangqs17@163.com��
 * @date        2023-07-06
 */
IntType **CalC2N(PolyGrid *grid)
{
    IntType **C2N = (IntType **)grid->GetC2N();
    // IF C2N has already existed
    if (C2N)
        return C2N;

    IntType i, j, k, n, fac, pi, count;
    IntType nTCell = grid->GetNTCell();
    IntType *nNPF = grid->GetnNPF();
    IntType *f2c = grid->Getf2c();
    IntType *nFPC = CalnFPC(grid);
    IntType **F2N = CalF2N(grid);
    IntType **C2F = CalC2F(grid);
    IntType *nNPC = CalnNPC(grid);

    // Check if nFPC is available
    if (!nFPC)
        nFPC = CalnFPC(grid);

    // Allocate memories for cell to face connection
    // C2N = NULL;
    snew_array_2D(C2N, nTCell, nNPC, true);
    for (i = 0; i < nTCell; i++)
    {

        if (nNPC[i] == 8 && nFPC[i] == 6)
        { // Ϊ����������cgns��ʽ����
            IntType Nmark[8];
            IntType fac6, npnt;
            fac = C2F[i][0];
            if (f2c[fac << 1] == i)
            {
                C2N[i][0] = F2N[fac][0];
                C2N[i][1] = F2N[fac][3];
                C2N[i][2] = F2N[fac][2];
                C2N[i][3] = F2N[fac][1];
            }
            else
            {
                C2N[i][0] = F2N[fac][0];
                C2N[i][1] = F2N[fac][1];
                C2N[i][2] = F2N[fac][2];
                C2N[i][3] = F2N[fac][3];
            }
            Nmark[0] = 1;
            Nmark[1] = 2;
            Nmark[2] = 4;
            Nmark[3] = 8;
            for (j = 4; j < 8; j++)
                Nmark[j] = 0;
            npnt = 4;

            for (j = 1; j < nFPC[i]; j++)
            {
                n = 0;
                fac = C2F[i][j];
                for (k = 0; k < nNPF[fac]; k++)
                {
                    pi = F2N[fac][k];
                    for (count = 0; count < 4; count++)
                    {
                        if (C2N[i][count] == pi)
                        {
                            n += Nmark[count];
                            break;
                        }
                    }
                }
                // if(n==3||n==6||n== 9||n==12) {
                if (n != 0)
                {
                    for (k = 0; k < nNPF[fac]; k++)
                    {
                        pi = F2N[fac][k];
                        if (C2N[i][0] != pi && C2N[i][1] != pi && C2N[i][2] != pi && C2N[i][3] != pi)
                        {
                            IntType mm;
                            for (mm = 4; mm < npnt; mm++)
                            {
                                if (C2N[i][mm] == pi)
                                    break;
                            }
                            Nmark[mm] -= n;
                            if (mm == npnt)
                            {
                                C2N[i][npnt] = pi;
                                npnt++;
                            }
                        }
                    }
                }
                else
                {
                    fac6 = fac;
                }
            }
            IntType p4, p5, p6, p7;
            for (k = 0; k < nNPF[fac6]; k++)
            {
                pi = F2N[fac6][k];
                for (j = 4; j < 8; j++)
                {
                    if (C2N[i][j] == pi)
                        break;
                }
                // assert(j < 8);
                if (Nmark[j] == -12)
                    p4 = C2N[i][j];
                else if (Nmark[j] == -9)
                    p5 = C2N[i][j];
                else if (Nmark[j] == -18)
                    p6 = C2N[i][j];
                else if (Nmark[j] == -21)
                    p7 = C2N[i][j];
                else
                {
                    // mflog::log.set_all_processors_out();
                    // mflog::log << "fatal error hex-cell !!" << endl;
                    // mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
                    cerr << "fatal error hex-cell !!" << endl;
                    exit(1);
                }
            }
            C2N[i][4] = p4;
            C2N[i][5] = p5;
            C2N[i][6] = p6;
            C2N[i][7] = p7;
        }
        else if (nNPC[i] == 6 && nFPC[i] == 5)
        { // Prism
            // ����������������tf1��tf2����6����
            IntType tf1, tf2;
            tf1 = -1;
            tf2 = -1;
            for (j = 0; j < nFPC[i]; ++j)
            {
                fac = C2F[i][j];
                if (nNPF[fac] == 4)
                {
                    continue;
                }
                if (tf1 == -1)
                {
                    // tf1û���ҵ�����ֵtf1
                    tf1 = fac;
                }
                else
                {
                    // tf1�Ѿ��ҵ�����ֵtf2
                    tf2 = fac;
                    break;
                }
            }

            // ��tf1��f2c����ķ����f2n[tf1]����C2Nǰ��������
            if (f2c[tf1 << 1] == i)
            {
                C2N[i][0] = F2N[tf1][0];
                C2N[i][2] = F2N[tf1][1];
                C2N[i][1] = F2N[tf1][2];
            }
            else
            {
                C2N[i][0] = F2N[tf1][0];
                C2N[i][1] = F2N[tf1][1];
                C2N[i][2] = F2N[tf1][2];
            }
            // ��tf2��f2c����ķ����f2n[tf1]����C2Nt[3](����������)
            IntType C2Nt[3];
            if (f2c[tf2 << 1] == i)
            {
                C2Nt[0] = F2N[tf2][0];
                C2Nt[1] = F2N[tf2][1];
                C2Nt[2] = F2N[tf2][2];
            }
            else
            {
                C2Nt[0] = F2N[tf2][0];
                C2Nt[2] = F2N[tf2][1];
                C2Nt[1] = F2N[tf2][2];
            }

            // �ҵ�����0��Ӧ�Ķ���3
            // ���������ı��ε��棬C2Nt�еĶ���Ͷ���0����������Ϊ����3
            IntType fc[3]; // ��ǹ���Ĵ���
            for (j = 0; j < 3; ++j)
            {
                fc[j] = 0;
            }

            for (j = 0; j < nFPC[i]; ++j)
            {
                fac = C2F[i][j];
                if (nNPF[fac] == 3)
                {
                    continue;
                }

                count = 0;
                for (k = 0; k < 4; ++k)
                {
                    if (F2N[fac][k] == C2N[i][0])
                    {
                        count = 1;
                        break;
                    }
                }
                if (count == 0)
                {
                    continue;
                } // ���治��������0

                for (k = 0; k < 4; ++k)
                {
                    if (F2N[fac][k] == C2Nt[0])
                    {
                        ++fc[0];
                    }
                    else if (F2N[fac][k] == C2Nt[1])
                    {
                        ++fc[1];
                    }
                    else if (F2N[fac][k] == C2Nt[2])
                    {
                        ++fc[2];
                    }
                }
            }

            IntType nid = -1;
            for (j = 0; j < 3; ++j)
            {
                if (fc[j] == 2)
                {
                    nid = j;
                    break;
                }
            }

            // ����������3��4��5������
            for (j = 0; j < 3; ++j)
            {
                C2N[i][j + 3] = C2Nt[(j + nid) % 3];
            }
        }
        else if (nNPC[i] == 5 && nFPC[i] == 5)
        { // Pyramid
            // �����ı�����rfa����4����
            IntType rfa;
            rfa = -1;
            for (j = 0; j < nFPC[i]; ++j)
            {
                fac = C2F[i][j];
                if (nNPF[fac] == 4)
                {
                    rfa = fac;
                    break;
                }
            }

            // ��tf1��f2c����ķ����f2n[tf1]����C2Nǰ�ĸ�����
            if (f2c[rfa << 1] == i)
            {
                C2N[i][0] = F2N[rfa][0];
                C2N[i][3] = F2N[rfa][1];
                C2N[i][2] = F2N[rfa][2];
                C2N[i][1] = F2N[rfa][3];
            }
            else
            {
                C2N[i][0] = F2N[rfa][0];
                C2N[i][1] = F2N[rfa][1];
                C2N[i][2] = F2N[rfa][2];
                C2N[i][3] = F2N[rfa][3];
            }

            // �ҵ�5������
            // �������������ε��棬�ҵ���ǰ�ĸ��㶼��ͬ�ĵ㣬��Ϊ���������
            IntType nid;
            nid = -1;
            for (j = 0; j < nFPC[i]; ++j)
            {
                fac = C2F[i][j];
                if (nNPF[fac] == 4)
                {
                    continue;
                }

                for (k = 0; k < 3; ++k)
                {
                    for (n = 0; n < 4; ++n)
                    {
                        if (F2N[fac][k] == C2N[i][n])
                        {
                            break;
                        }
                    }
                    if (n == 4)
                    { // �ҵ���ǰ�ĸ��㶼���ȵĵ�
                        nid = k;
                        break;
                    }
                }

                if (nid != -1)
                {
                    C2N[i][4] = F2N[fac][nid];
                    break;
                }
            }
        }
        else if (nNPC[i] == 4 && nFPC[i] == 4)
        { // Tetra
            // ֱ��ʹ�õ�һ����ȷ��������ǰ��������
            fac = C2F[i][0];

            // ��fac��f2c����ķ����f2n[fac]����C2Nǰ��������
            if (f2c[fac << 1] == i)
            {
                C2N[i][0] = F2N[fac][0];
                C2N[i][2] = F2N[fac][1];
                C2N[i][1] = F2N[fac][2];
            }
            else
            {
                C2N[i][0] = F2N[fac][0];
                C2N[i][1] = F2N[fac][1];
                C2N[i][2] = F2N[fac][2];
            }

            // �ҵ�4������
            // ʹ�õڶ����棬�ҵ���ǰ�����㶼��ͬ�ĵ㣬��Ϊ���ĸ�����
            IntType nid = -1;
            fac = C2F[i][1];

            for (k = 0; k < 3; ++k)
            {
                for (n = 0; n < 3; ++n)
                {
                    if (F2N[fac][k] == C2N[i][n])
                    {
                        break;
                    }
                }
                if (n == 3)
                { // �ҵ���ǰ3���㶼���ȵĵ�
                    nid = k;
                    break;
                }
            }
            // assert(nid != -1);
            C2N[i][3] = F2N[fac][nid];
        }

        // other cell type may be builded by multigrid's coarsen
        // so the irregular cell can not be ordered.
        else
        {
            count = 0;
            for (j = 0; j < nFPC[i]; j++)
            {
                fac = C2F[i][j];
                for (k = 0; k < nNPF[fac]; k++)
                {
                    pi = F2N[fac][k];
                    for (n = 0; n < count; n++)
                    {
                        if (C2N[i][n] == pi)
                            break;
                    }
                    if (n == count)
                    {
                        C2N[i][count++] = pi;
                    }
                }
            }
        } // end cycle cell type
    }
    // Attach C2N to the grid and return it
    grid->SetC2N(C2N);
    return C2N;
}

/*!
 * @brief       Calculate nNPC (number of node in each real cell). Real cell is "physical" cell, ghost cell not included
 * @param       grid
 * @return      IntType*
 * @note        This function also calculates C2N, but the node index of one cell is arbitrary, the topology of cell to node is deserted.
 * @note        This function is almost the same as CalnCPC. But we can't use it, because nNPC doesn't count boundary cells except for interfaces.
 *
 * @author      Wisces ��wwwangqs17@163.com��
 * @date        2023-07-06
 */
IntType *CalnNPC(PolyGrid *grid)
{
    IntType *nNPC = (IntType *)grid->GetnNPC();
    // IF nNPC has already existed
    if (nNPC)
        return nNPC;

    IntType i, j, c1, c2, count;
    IntType nTCell = grid->GetNTCell();
    IntType nBFace = grid->GetNBFace();
    IntType nTFace = grid->GetNTFace();
    IntType *f2c = grid->Getf2c();
    IntType *nNPF = grid->GetnNPF();
    IntType **F2N = CalF2N(grid);

    // Allocate memories for number of faces per cell
    IntType *nNPC_tmp = NULL;
    snew_array_1D(nNPC_tmp, nTCell);
    // set nFPC to zero
    for (i = 0; i < nTCell; i++)
        nNPC_tmp[i] = 0;

    for (i = 0; i < nBFace; i++)
    {
        c1 = f2c[i + i];
        nNPC_tmp[c1] += nNPF[i];
    }
    count = nBFace * 2;
    for (i = nBFace; i < nTFace; i++)
    {
        c1 = f2c[count++];
        c2 = f2c[count++];
        nNPC_tmp[c1] += nNPF[i];
        nNPC_tmp[c2] += nNPF[i];
    }

    IntType **C2N_tmp = NULL;
    snew_array_2D(C2N_tmp, nTCell, nNPC_tmp, true);
    for (i = 0; i < nTCell; i++)
        nNPC_tmp[i] = 0;
    for (i = 0; i < nBFace; i++)
    {
        c1 = f2c[i + i];
        for (j = 0; j < nNPF[i]; j++)
        {
            C2N_tmp[c1][nNPC_tmp[c1]++] = F2N[i][j];
        }
    }
    count = nBFace * 2;
    for (i = nBFace; i < nTFace; i++)
    {
        c1 = f2c[count++];
        c2 = f2c[count++];
        for (j = 0; j < nNPF[i]; j++)
        {
            C2N_tmp[c1][nNPC_tmp[c1]++] = F2N[i][j];
            C2N_tmp[c2][nNPC_tmp[c2]++] = F2N[i][j];
        }
    }

    // nNPC = NULL;
    snew_array_1D(nNPC, nTCell);

    // set nNPC to 1
    for (i = 0; i < nTCell; i++)
        nNPC[i] = 1;
    // Boundary faces
    for (i = 0; i < nTCell; i++)
    {
        grid->quick_sort(C2N_tmp[i], 0, nNPC_tmp[i] - 1);
        for (j = 1; j < nNPC_tmp[i]; j++)
        {
            if (C2N_tmp[i][j] != C2N_tmp[i][j - 1])
                nNPC[i]++;
        }
    }

    // Attach nNPC to the grid and return it
    grid->SetnNPC(nNPC);
    sdel_array_1D(nNPC_tmp);
    sdel_array_2D(C2N_tmp);
    return nNPC;
}

/**
// ///////////////////////////////////////////////////
// ��һ�������Ž�����������
// ///////////////////////////////////////////////////
void ReorderPnts(IntType *pnt, IntType npt)
{
    IntType idx = 0;
    IntType is = 0, ie = npt;
    for (IntType i = is + 1; i < ie; i++)
    {
        if (pnt[i] < pnt[idx])
            idx = i;
    }
    IntType *itmp = NULL;
    snew_array_1D(itmp, npt);
    for (IntType i = 0; i < npt; i++)
        itmp[i] = pnt[i];
    for (IntType i = 0; i < npt; i++)
        pnt[i] = itmp[(i + idx) % npt]; // zhyb: �����������С���ҳ����ڵ�һλ,��˳��û�б�

    sdel_array_1D(itmp);
}

void BreakFaceLoop(PolyGrid *grid)
{
    if (!grid->GetZ())
        return;

    register IntType i, j, k, n = 1, p1;

    while (n > 0)
    {
        IntType *nNPF = grid->GetnNPF(), *f2n = grid->Getf2n();
        IntType nTFace = grid->GetNTFace();
        IntType *f2nind = NULL;
        snew_array_1D(f2nind, nTFace + 1);
        f2nind[0] = 0;
        for (i = 0; i < nTFace; i++)
            f2nind[i + 1] = f2nind[i] + nNPF[i];

        n = 0;
        for (i = 0; i < nTFace; i++)
        {
            p1 = f2n[f2nind[i]];
            for (j = f2nind[i] + 1; j < f2nind[i + 1]; j++)
            {
                if (f2n[j] == p1)
                {
                    n++;
                    break;
                }
            }
        }

        if (n == 0)
        {
            sdel_array_1D(f2nind);
            break;
        }

        // now break the multi-looped faces into two faces
        IntType n_face = nTFace + n;
        IntType *nnpf = NULL;
        snew_array_1D(nnpf, n_face);
        IntType *f2c = grid->Getf2c();
        IntType *f2cn = NULL;
        snew_array_1D(f2cn, n_face * 2);
        IntType ii, nf2;

        n = 0;
        for (i = 0; i < nTFace; i++)
        {
            nnpf[i] = nNPF[i];
            ii = i + i;
            f2cn[ii] = f2c[ii];
            f2cn[ii + 1] = f2c[ii + 1];
            p1 = f2n[f2nind[i]];
            for (j = f2nind[i] + 1; j < f2nind[i + 1]; j++)
            {
                if (f2n[j] == p1)
                {
                    nnpf[i] = j - f2nind[i];
                    nf2 = (nTFace + n) * 2;
                    f2cn[nf2] = f2c[ii];
                    f2cn[nf2 + 1] = f2c[ii + 1];
                    nnpf[nTFace + n++] = nNPF[i] - j + f2nind[i] - 2;
                    break;
                }
            }
        }

        // get new f2n
        IntType *f2nnind = NULL;
        snew_array_1D(f2nnind, n_face + 1);
        f2nnind[0] = 0;

        for (i = 0; i < n_face; i++)
            f2nnind[i + 1] = f2nnind[i] + nnpf[i];
        IntType *f2nn = NULL;
        snew_array_1D(f2nn, f2nnind[n_face]);

        n = 0;
        for (i = 0; i < nTFace; i++)
        {
            p1 = f2n[f2nind[i]];
            f2nn[f2nnind[i]] = p1;

            for (j = 1; j < nNPF[i]; j++)
            {
                if (f2n[f2nind[i] + j] == p1)
                {
                    // add the face to the end
                    assert(f2n[f2nind[i] + j + 1] == f2n[f2nind[i + 1] - 1]);

                    for (k = j + 1; k < nNPF[i] - 1; k++)
                        f2nn[f2nnind[nTFace + n] + k - j - 1] = f2n[f2nind[i] + k];
                    n++;
                    break;
                }
                else
                {
                    f2nn[f2nnind[i] + j] = f2n[f2nind[i] + j];
                }
            }
        }
        sdel_array_1D(f2nnind);
        sdel_array_1D(f2nind);
        grid->SetNTFace(n_face);
        grid->SetnNPF(nnpf);
        grid->Setf2n(f2nn);
        grid->Setf2c(f2cn);

        printf("No of faces having more than one loop is %ld\n", (long)n);
    }
}
*/

/*!
 * @brief       find pure symmetry boundary node for node Gauss gradient compute
 * @param       grid
 *
 * @author      Wisces ��wwwangqs17@163.com��
 * @date        2023-07-10
 */
void FindNodeSYMM(PolyGrid *grid)
{
    IntType i, j, p1, type;

    IntType nBFace = grid->GetNBFace();
    IntType nTNode = grid->GetNTNode();
    IntType *nNPF = grid->GetnNPF();
    IntType **F2N = CalF2N(grid);
    RealGeom *xfn = grid->GetXfn();
    RealGeom *yfn = grid->GetYfn();
    RealGeom *zfn = grid->GetZfn();
    BCRecord **bcr = grid->Getbcr();

    IntType *nodesymm = NULL;
    RealGeom *xfn_n_symm = NULL, *yfn_n_symm = NULL, *zfn_n_symm = NULL;

    snew_array_1D(nodesymm, nTNode);
    // grid->UpdateDataPtr(nodesymm, INT, nTNode, "nodesymm");
    // nodesymm=0 : not symmetry boundary node
    // nodesymm=1 : symmetry boundary node
    snew_array_1D(xfn_n_symm, nTNode);
    snew_array_1D(yfn_n_symm, nTNode);
    snew_array_1D(zfn_n_symm, nTNode);
    // grid->UpdateDataPtr(xfn_n_symm, REAL_GEOM, nTNode, "xfn_n_symm");
    // grid->UpdateDataPtr(yfn_n_symm, REAL_GEOM, nTNode, "yfn_n_symm");
    // grid->UpdateDataPtr(zfn_n_symm, REAL_GEOM, nTNode, "zfn_n_symm");

    RealGeom *nface = NULL;
    snew_array_1D(nface, nTNode);
    for (i = 0; i < nTNode; i++)
    {
        nodesymm[i] = 2;
        xfn_n_symm[i] = 0.0;
        yfn_n_symm[i] = 0.0;
        zfn_n_symm[i] = 0.0;
        nface[i] = 0.0;
    }

    for (i = 0; i < nBFace; i++)
    {
        type = bcr[i]->GetType();
        if (type != SYMM)
            continue;

        for (j = 0; j < nNPF[i]; j++)
        {
            p1 = F2N[i][j];
            nodesymm[p1] = 1;
            xfn_n_symm[p1] += xfn[i];
            yfn_n_symm[p1] += yfn[i];
            zfn_n_symm[p1] += zfn[i];
            nface[p1] += 1.0;
        }
    }
#ifdef MF_MPICH
    grid->CommInternodeDataMPI2(nodesymm);  // ȡСֵ
    grid->CommInternodeDataMPI(xfn_n_symm); // ���
    grid->CommInternodeDataMPI(yfn_n_symm);
    grid->CommInternodeDataMPI(zfn_n_symm);
    grid->CommInternodeDataMPI(nface);
#endif

    for (i = 0; i < nTNode; i++)
    {
        if (nodesymm[i] == 2)
            nodesymm[i] = 0;
        if (nodesymm[i] == 1)
        {
            xfn_n_symm[i] /= nface[i];
            yfn_n_symm[i] /= nface[i];
            zfn_n_symm[i] /= nface[i];
        }
    }
    sdel_array_1D(nface);

    grid->SetNodeSymm(nodesymm);
    grid->SetXfnNSymm(xfn_n_symm);
    grid->SetYfnNSymm(yfn_n_symm);
    grid->SetZfnNSymm(zfn_n_symm);
}
/**
/// \brief  ����ϸ����ڵ����ڵĴ�������ĸ���
/// \par    Update records:
/// <pre>
/// Date        Author      Description
///
/// </pre>
IntType *CalnCCPN(PolyGrid *grid)
{
    IntType *nCCPN = grid->GetnCCPN();
    if (nCCPN)
        return nCCPN;

    IntType **N2CC = CalN2CC(grid);
    nCCPN = grid->GetnCCPN();

    return nCCPN;
}

/// \brief  ����ϸ����ڵ����ڵĴ�����������
/// \par    Update records:
/// <pre>
/// Date        Author      Description
/// 20220217    ���½�        ɾ���� NinCN �ļ���
/// </pre>
IntType **CalN2CC(PolyGrid *grid)
{
    IntType **N2CC = grid->GetN2CC();
    if (N2CC)
        return N2CC;

    IntType *c2cc = grid->Getc2cc();
    assert(c2cc != 0);

    IntType i, j, k, c1, cc;
    IntType nTNode = grid->GetNTNode();

    IntType *nVPN = CalnVPN(grid);
    IntType **N2V = CalN2V(grid);

    // Allocate memories for number of CCell per node
    IntType *nCCPN = (IntType *)grid->GetnCCPN();
    if (nCCPN == 0)
        snew_array_1D(nCCPN, nTNode);

    // set nCCPN to zero
    for (i = 0; i < nTNode; i++)
        nCCPN[i] = 0;

    IntType **N2CC_tmp = NULL;
    snew_array_2D(N2CC_tmp, nTNode, nVPN, false);
    for (i = 0; i < nTNode; i++)
    {
        for (j = 0; j < nVPN[i]; j++)
        {
            N2CC_tmp[i][j] = -1;
        }
    }

    for (i = 0; i < nTNode; i++)
    {
        for (j = 0; j < nVPN[i]; j++)
        {
            c1 = N2V[i][j];
            cc = c2cc[c1];
            for (k = 0; k < nCCPN[i]; k++)
            {
                if (cc == N2CC_tmp[i][k])
                    break;
            }
            if (k == nCCPN[i])
            {
                N2CC_tmp[i][k] = cc;
                nCCPN[i]++;
            }
        }
    }
    // allocate memory for **N2CC
    N2CC = NULL;
    snew_array_2D(N2CC, nTNode, nCCPN, true);
    for (i = 0; i < nTNode; i++)
    {
        for (j = 0; j < nCCPN[i]; j++)
        {
            N2CC[i][j] = N2CC_tmp[i][j];
        }
    }
    sdel_array_2D(N2CC_tmp, nTNode, false); ////�Ƿ�ɾ��N2CC_tmp;lihuan-2018-11-15

    // Attach C2N to the grid and return it
    grid->SetnCCPN(nCCPN);
    grid->SetN2CC(N2CC);

    return N2CC;
}

//  Calculate nVPN (number of faces in each )
//    This function is almost the same as CalnPC.
//    And nVPN doesn't count ghost cells.
IntType *CalnVPN(PolyGrid *grid)
{
    IntType *nVPN = grid->GetnVPN();
    if (nVPN)
        return nVPN;

    IntType nTNode = grid->GetNTNode();
    IntType nTCell = grid->GetNTCell();

    nVPN = NULL;
    snew_array_1D(nVPN, nTNode);
    assert(nVPN != 0);
    grid->SetnVPN(nVPN);

    IntType *nNPC = grid->GetnNPC();
    assert(nNPC != 0);
    IntType **C2N = grid->GetC2N();
    assert(C2N != 0);

    IntType i, j;
    for (i = 0; i < nTNode; i++)
        nVPN[i] = 0;

    for (i = 0; i < nTCell; i++)
    {
        for (j = 0; j < nNPC[i]; j++)
        {
            IntType jN = C2N[i][j];
            nVPN[jN]++;
        }
    }

    return nVPN;
}

//  Calculate N2V (number of faces in each )
//    This function is almost the same as CalC2C.
//    And N2V doesn't count ghost cells.
IntType **CalN2V(PolyGrid *grid)
{
    IntType **N2V = grid->GetN2V();
    if (N2V)
        return N2V;

    IntType nTNode = grid->GetNTNode();
    IntType nTCell = grid->GetNTCell();

    IntType *nVPN = grid->GetnVPN();
    if (!nVPN)
        nVPN = CalnVPN(grid);

    IntType i, j;

    N2V = NULL;
    snew_array_2D(N2V, nTNode, nVPN, true);
    grid->SetN2V(N2V);
    IntType *nNPC = grid->GetnNPC();
    assert(nNPC != 0);
    IntType **C2N = grid->GetC2N();
    assert(C2N != 0);

    for (i = 0; i < nTNode; i++)
        nVPN[i] = 0;
    for (i = 0; i < nTCell; i++)
    {
        for (j = 0; j < nNPC[i]; j++)
        {
            IntType jN = C2N[i][j];
            N2V[jN][nVPN[jN]] = i;
            nVPN[jN]++;
        }
    }

    return N2V;
}
*/

/*!
 * @brief       calculate N2C(Node to cell connection)
 * @param       grid
 * @return      IntType**
 *
 * @author      Wisces ��wwwangqs17@163.com��
 * @date        2023-07-09
 */
IntType **CalN2C(PolyGrid *grid)
{
    IntType **N2C = grid->GetN2C();
    if (N2C)
        return N2C;

    IntType *nCPN = CalnCPN(grid);
    N2C = grid->GetN2C();
    return N2C;
}

/*!
 * @brief       Calculate nCPN(number of cells in each node). Not include all ghost cell!
 * @param       grid
 * @return      IntType*
 *
 * @author      Wisces ��wwwangqs17@163.com��
 * @date        2023-07-09
 */
IntType *CalnCPN(PolyGrid *grid)
{
    IntType *nCPN = (IntType *)grid->GetnCPN();
    // IF nCPN has already existed
    if (nCPN)
        return nCPN;

    IntType i, j, k, c1, c2, p1, mark;

    IntType nBFace = grid->GetNBFace();
    IntType nTFace = grid->GetNTFace();
    IntType nTNode = grid->GetNTNode();
    IntType *f2c = grid->Getf2c();
    IntType **F2N = CalF2N(grid);
    IntType *nNPF = grid->GetnNPF();

    // nCPN = NULL;
    snew_array_1D(nCPN, nTNode);
    IntType **n2c = NULL;
    snew_array_1D(n2c, nTNode);
    for (i = 0; i < nTNode; i++)
    {
        n2c[i] = NULL;
        snew_array_1D(n2c[i], 500);
        nCPN[i] = 0;
        for (j = 0; j < 500; j++)
        {
            n2c[i][j] = -1;
        }
    }

    for (i = 0; i < nBFace; i++)
    {
        c1 = f2c[i + i];

        for (j = 0; j < nNPF[i]; j++)
        {
            p1 = F2N[i][j];
            mark = 0;
            for (k = 0; k < nCPN[p1]; k++)
            {
                if (n2c[p1][k] == c1)
                {
                    mark = 1;
                    break;
                }
            }
            if (mark == 0)
            {
                n2c[p1][nCPN[p1]] = c1;
                nCPN[p1]++;
            }
        }
    }
    for (i = nBFace; i < nTFace; i++)
    {
        c1 = f2c[i + i];
        c2 = f2c[i + i + 1];

        for (j = 0; j < nNPF[i]; j++)
        {
            p1 = F2N[i][j];
            mark = 0;
            for (k = 0; k < nCPN[p1]; k++)
            {
                if (n2c[p1][k] == c1)
                {
                    mark = 1;
                    break;
                }
            }
            if (mark == 0)
            {
                n2c[p1][nCPN[p1]] = c1;
                nCPN[p1]++;
            }
            mark = 0;
            for (k = 0; k < nCPN[p1]; k++)
            {
                if (n2c[p1][k] == c2)
                {
                    mark = 1;
                    break;
                }
            }
            if (mark == 0)
            {
                n2c[p1][nCPN[p1]] = c2;
                nCPN[p1]++;
            }
        }
    }

    IntType **N2C = NULL;
    snew_array_2D(N2C, nTNode, nCPN, false);
    for (i = 0; i < nTNode; i++)
    {
        for (j = 0; j < nCPN[i]; j++)
        {
            N2C[i][j] = n2c[i][j];
        }
    }

    grid->SetnCPN(nCPN);
    grid->SetN2C(N2C);

    // IntType maxncpn = 0, minncpn = 10000;
    // for (i = 0; i < nTNode; i++)
    // {
    //     maxncpn = MAX(maxncpn, nCPN[i]);
    //     minncpn = MIN(minncpn, nCPN[i]);
    // }
    // #ifdef PARALLEL_MPI
    //     Parallel::parallel_min_max(minncpn, maxncpn, MPI_COMM_WORLD);
    // #endif

    // mflog::log.set_one_processor_out();
    // mflog::log << endl
    //            << "Max nCPN = " << maxncpn << "  Min nCPN = " << minncpn << endl;
    // cout << endl
    //      << "Max nCPN = " << maxncpn << "  Min nCPN = " << minncpn << endl;

    for (i = 0; i < nTNode; i++)
        sdel_array_1D(n2c[i]);
    sdel_array_1D(n2c);
    return nCPN;
}
