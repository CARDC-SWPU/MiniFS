/*!
 * @file        cal_flux.cpp
 * @brief       无粘通量计算的主要函数
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
 * <tr><td> 2023-06-11  <td> 1.1      <td> Wisces  <td> 添加了部分函数，同时对部分函数进行了修改
 * <tr><td> 2023-06-15  <td> 2.0      <td> Wisces  <td> 添加头文件 para_field_global.h，采用全局变量形式存储数据
 * <tr><td> 2023-07-03  <td> 3.0      <td> Wisces  <td> 添加了 OpenMP 并行（使用条件编译）
 * <tr><td> 2023-07-12  <td> 4.0      <td> Wisces  <td> 添加了 TDTree 线程级并行（使用条件编译）
 * </table>
 */

//!< C/C++ head files
#include <iostream>
// #include <cassert>
#include <cmath>
#include <cstdlib> ///< exit()

//!< direct head file
#include "cal_flux.h"

//!< user defined head files
#include "grid_polyhedra.h"
#include "memory_util.h"
#include "para_field_global.h"

//!< head file relying on condition-compiling
#if (defined MF_OPENMP) || (defined TDTREE)
#include <omp.h>
#endif

#ifdef TDTREE
#include "TDTree.h"
#endif

using namespace std;

/*!
 * @brief       Zero the residuals. Also allocate memory for the residuals if they had not been allocated
 * @param       grid
 * @remarks     modify according to the fun [void ZeroResiduals(PolyGrid *grid)]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-11
 */
void ZeroResiduals(PolyGrid *grid)
{
    IntType nTCell = grid->GetNTCell(), nT5 = 5 * nTCell;
    RealFlow *res = grid->GetRes();
    if (!res)
    {
        snew_array_1D(res, nT5);
        grid->SetRes(res);
        // assert(res != 0);
    }

#if (defined MF_OPENMP) || (defined TDTREE)
#pragma omp parallel for
#endif
    for (IntType i = 0; i < nT5; i++)
    {
        res[i] = 0.;
    }
}

/*!
 * @brief       Drive individual functions to calculate inviscid fluxes in 3D
 * @param       grid
 * @param       limit   (= NULL)
 * @param       level   (= 0)
 * @remarks     modify according to the fun [void InviscidFlux(PolyGrid *grid, RealFlow **limit, IntType level)]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-05-31
 */
void InviscidFlux(PolyGrid *grid, RealFlow **limit, IntType level)
{
    IntType i, ns, ne, len;
    IntType nTCell = grid->GetNTCell();
    IntType nBFace = grid->GetNBFace();
    IntType nTFace = grid->GetNTFace();
    IntType n = nTCell + nBFace;

    // Get metrics
    RealGeom *xfn = grid->GetXfn();
    RealGeom *yfn = grid->GetYfn();
    RealGeom *zfn = grid->GetZfn();
    RealGeom *area = grid->GetFaceArea();
    RealGeom *vgn = grid->GetFaceNormalVelocity();
    // for overlap
    // IntType *face_act = NULL;

    // Allocate temporary memories for ql, qr and flux
    RealFlow *ql[5], *qr[5], *flux[5], *q[5];

    // Get flow variables
    // q[0] = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "rho");
    // q[1] = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "u");
    // q[2] = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "v");
    // q[3] = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "w");
    // q[4] = (RealFlow *)grid->GetDataPtr(REAL_FLOW, n, "p");
    // RealFlow *rho = grid->GetRho();
    // RealFlow *u = grid->GetU();
    // RealFlow *v = grid->GetV();
    // RealFlow *w = grid->GetW();
    // RealFlow *p = grid->GetP();

    q[0] = (RealFlow *)grid->GetRho();
    q[1] = (RealFlow *)grid->GetU();
    q[2] = (RealFlow *)grid->GetV();
    q[3] = (RealFlow *)grid->GetW();
    q[4] = (RealFlow *)grid->GetP();

// const IntType kNVar = 5;
// RealFlow **dqdx = NULL, **dqdy = NULL, **dqdz = NULL;
// snew_array_1D(dqdx, kNVar);
// snew_array_1D(dqdy, kNVar);
// snew_array_1D(dqdz, kNVar);
// dqdx[0] = static_cast<RealFlow *>(
//     grid->GetDataPtr(REAL_FLOW, kNVar * n, "dqdx"));
// dqdy[0] = static_cast<RealFlow *>(
//     grid->GetDataPtr(REAL_FLOW, kNVar * n, "dqdy"));
// dqdz[0] = static_cast<RealFlow *>(
//     grid->GetDataPtr(REAL_FLOW, kNVar * n, "dqdz"));
// for (IntType i = 1; i < kNVar; ++i)
// {
//     dqdx[i] = &dqdx[i - 1][n];
//     dqdy[i] = &dqdy[i - 1][n];
//     dqdz[i] = &dqdz[i - 1][n];
// }
#ifdef MF_OPENMP
    len = nTFace;
#else
    len = SEG_LEN;
#endif
    ql[0] = NULL;
    qr[0] = NULL;
    flux[0] = NULL;
    snew_array_1D(ql[0], 5 * len);
    snew_array_1D(qr[0], 5 * len);
    snew_array_1D(flux[0], 5 * len);
    // assert(ql[0] != 0);
    // assert(qr[0] != 0);
    // assert(flux[0] != 0);

    for (i = 1; i < 5; i++)
    {
        ql[i] = &ql[i - 1][len];
        qr[i] = &qr[i - 1][len];
        flux[i] = &flux[i - 1][len];
    }

    ns = 0;
    do
    {
#ifdef MF_OPENMP
        ne = nTFace;
#else
        ne = ns + SEG_LEN;
        if (ne > nTFace)
            ne = nTFace;
#endif
        // Get left variables and right variables
        SetQlQrWithQ(grid, q, ql, qr, ns, ne);
        // if (limit != NULL)
        // {
        //     CalcuQlQr(grid, ql, qr, limit, dqdx, dqdy, dqdz, ns, ne);
        // }
        // ModQlQrBou(grid, ql, qr, ns, ne);

        // CompInvFlux(grid, ql, qr, flux, &xfn[ns], &yfn[ns], &zfn[ns], &area[ns], &vgn[ns],
        //             &face_act[ns], gam, p_bar, alf_l, alf_n, 0, ns, ne);
        CompInvFlux(grid, ql, qr, flux, &xfn[ns], &yfn[ns], &zfn[ns], &area[ns], &vgn[ns],
                    NULL, gam, p_bar, alf_l, alf_n, 0, ns, ne);

        // Load the fluxes to residuals
        LoadFlux(grid, flux, ns, ne);
        ns = ne;
    } while (ns < nTFace);

    sdel_array_1D(ql[0]);
    sdel_array_1D(qr[0]);
    sdel_array_1D(flux[0]);
}

/*!
 * @brief       Set the Ql Qr With Q object
 * @param       grid
 * @param       q
 * @param       ql
 * @param       qr
 * @param       ns
 * @param       ne
 * @remarks     modify according to the fun [void SetQlQrWithQ(PolyGrid *grid, RealFlow *q[], RealFlow *ql[], RealFlow *qr[], IntType ns, IntType ne)]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-05-31
 */
void SetQlQrWithQ(PolyGrid *grid, RealFlow *q[], RealFlow *ql[], RealFlow *qr[], IntType ns, IntType ne)
{
    IntType *f2c = grid->Getf2c();
    IntType nvar, i, c1, c2, count, n, face;

    nvar = 5;
#ifdef MF_OPENMP
#pragma omp parallel for private(face, count, i, c1, c2, n)
    for (face = ns; face < ne; face++)
    {
        count = 2 * face;
        i = face - ns;
        c1 = f2c[count++];
        c2 = f2c[count];
        for (n = 0; n < nvar; n++)
        {
            ql[n][i] = q[n][c1];
            qr[n][i] = q[n][c2];
        }
    }
#else
    for (n = 0; n < nvar; n++)
    {
        count = 2 * ns;
        i = 0;
        for (face = ns; face < ne; face++)
        {
            c1 = f2c[count++];
            c2 = f2c[count++];

            ql[n][i] = q[n][c1];
            qr[n][i] = q[n][c2];
            i++;
        }
    }
#endif
}

/*!
 * @brief       Compute fluxes in 3D according to its type for preconditioning
 * @param       grid
 * @param       ql
 * @param       qr
 * @param       flux
 * @param       xfn
 * @param       yfn
 * @param       zfn
 * @param       area
 * @param       vgn         （未使用）
 * @param       face_act
 * @param       gam
 * @param       p_bar
 * @param       alf_l
 * @param       alf_n
 * @param       type_flux
 * @param       ns
 * @param       ne
 * @remarks     modify according to the fun [void CompInvFlux(PolyGrid *grid, RealFlow *ql[5], RealFlow *qr[5], RealFlow *flux[5], RealGeom *xfn, RealGeom *yfn, RealGeom *zfn, RealGeom *area, RealGeom *vgn, IntType *face_act, RealFlow gam, RealFlow p_bar, RealFlow alf_l, RealFlow alf_n, IntType type_flux, IntType ns, IntType ne)]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-05-31
 */
void CompInvFlux(PolyGrid *grid, RealFlow *ql[5], RealFlow *qr[5], RealFlow *flux[5],
                 RealGeom *xfn, RealGeom *yfn, RealGeom *zfn, RealGeom *area,
                 RealGeom *vgn, IntType *face_act, RealFlow gam, RealFlow p_bar,
                 RealFlow alf_l, RealFlow alf_n, IntType type_flux,
                 IntType ns, IntType ne)
{
    RoeFlux_noprec(grid, ql, qr, flux, xfn, yfn, zfn,
                   area, face_act, gam, p_bar, alf_l, alf_n, ns, ne);
}

/*!
 * @brief       Compute inviscid fluxes Using Roe scheme
 * @param       grid
 * @param       ql
 * @param       qr
 * @param       flux
 * @param       xfn
 * @param       yfn
 * @param       zfn
 * @param       area
 * @param       face_act
 * @param       gam
 * @param       p_bar
 * @param       alf_l
 * @param       alf_n
 * @param       ns
 * @param       ne
 * @remarks     modify according to the fun [void RoeFlux_noprec(PolyGrid *grid, RealFlow *ql[5], RealFlow *qr[5], RealFlow *flux[5], RealGeom *xfn, RealGeom *yfn, RealGeom *zfn, RealGeom *area, IntType *face_act, RealFlow gam, RealFlow p_bar, RealFlow alf_l, RealFlow alf_n, IntType ns, IntType ne)]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-05-31
 */
void RoeFlux_noprec(PolyGrid *grid, RealFlow *ql[5], RealFlow *qr[5], RealFlow *flux[5],
                    RealGeom *xfn, RealGeom *yfn, RealGeom *zfn, RealGeom *area, IntType *face_act,
                    RealFlow gam, RealFlow p_bar, RealFlow alf_l, RealFlow alf_n,
                    IntType ns, IntType ne)
{
    // register IntType i, c1, c2, count, ni;
    IntType len;
    // RealFlow rho_a, u_a, v_a, w_a, h_a, c_a, c2_a, vn_a, q2;
    // RealFlow vn_l, et_l, ht_l, vn_r, et_r, ht_r /*, gamm1*/;
    // RealFlow tmp0, tmp1, tmp2, alpha1, alpha2, alpha3, eigv1, eigv2, eigv3;
    // RealFlow drho, du, dv, dw, dp, dvn, dq2;
    // RealGeom areax, areay, areaz;
    // RealFlow spectral, epsaa, epsbb, epscc, epsa_r;
    // RealFlow u_vgn, v_vgn, w_vgn;

    IntType nTCell = grid->GetNTCell();
    IntType nBFace = grid->GetNBFace();
    IntType nTFace = grid->GetNTFace();
    IntType n = nTCell + nBFace;
    IntType *f2c = grid->Getf2c();

    RealGeom *vgn = grid->GetFaceNormalVelocity();

    // IntType *IsNormalFace = 0;
    // IntType *IsShockFace = 0;
    //!< Wisces: 暂时简化为 EntropyCorType == 3
    /**
    if (EntropyCorType == 4)
    {
        IsNormalFace = (IntType *)grid->GetDataPtr(INT, nTFace, "IsNormalFace");
        if (!IsNormalFace)
        {
            grid->FindNormalFace();
            IsNormalFace = (IntType *)grid->GetDataPtr(INT, nTFace, "IsNormalFace");
        }

        // shock face or not
        IsShockFace = NULL;
        snew_array_1D(IsShockFace, ne - ns);
        // CalIsShockFace(grid, IsShockFace, ns, ne);
    }
    */

    len = ne - ns;

#ifdef FS_SIMD
    //!< @todo 暂时什么都不做
/**
    // containing SIMD
    // const IntType    Vec = 8;
    IntType i, ni;
    IntType k;
#ifdef MF_OPENMP // OpenMP && SIMD
#pragma omp parallel for private(ni)
#endif
    for (i = 0; i < len - Vec; i += Vec)
    {
        ni = ns + i;
        // IntType ni_v[Vec];//areax_v[Vec], areay_v[Vec], areaz_v[Vec];
        RealFlow et_l_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow et_r_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow ht_l_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow ht_r_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow vn_l_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow vn_r_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow tmp0_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow tmp1_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow tmp2_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow rho_a_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow u_a_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow v_a_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow w_a_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow vn_a_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow h_a_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow q2_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow c2_a_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow c_a_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow eigv1_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow eigv2_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow eigv3_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow epsa_r_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow spectral_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow u_vgn_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow v_vgn_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow w_vgn_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow epsaa_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow epsbb_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow epscc_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow drho_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow du_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow dv_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow dw_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow dp_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow dvn_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow dq2_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow alpha1_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow alpha2_v[Vec] __attribute__((aligned(ALIGN)));
        RealFlow alpha3_v[Vec] __attribute__((aligned(ALIGN)));
#pragma omp simd safelen(Vec)
        for (IntType iv = 0; iv < Vec; iv++)
        {
            // Total energy
            et_l_v[iv] = (ql[4][i + iv] + p_bar) / gamm1 + 0.5 * ql[0][i + iv] *
                                                               (ql[1][i + iv] * ql[1][i + iv] + ql[2][i + iv] * ql[2][i + iv] + ql[3][i + iv] * ql[3][i + iv]);
            et_r_v[iv] = (qr[4][i + iv] + p_bar) / gamm1 + 0.5 * qr[0][i + iv] *
                                                               (qr[1][i + iv] * qr[1][i + iv] + qr[2][i + iv] * qr[2][i + iv] + qr[3][i + iv] * qr[3][i + iv]);
            ht_l_v[iv] = et_l_v[iv] + ql[4][i + iv] + p_bar;
            ht_r_v[iv] = et_r_v[iv] + qr[4][i + iv] + p_bar;
            // Full flux
            vn_l_v[iv] = xfn[i + iv] * ql[1][i + iv] + yfn[i + iv] * ql[2][i + iv] + zfn[i + iv] * ql[3][i + iv];
            vn_r_v[iv] = xfn[i + iv] * qr[1][i + iv] + yfn[i + iv] * qr[2][i + iv] + zfn[i + iv] * qr[3][i + iv];
            if (!steady)
            { // unsteady
                vn_l_v[iv] -= vgn[ni + iv];
                vn_r_v[iv] -= vgn[ni + iv];
            }
            tmp0_v[iv] = vn_l_v[iv] * ql[0][i + iv];
            tmp1_v[iv] = vn_r_v[iv] * qr[0][i + iv];

            flux[0][i + iv] = tmp0_v[iv] + tmp1_v[iv];
            flux[1][i + iv] = tmp0_v[iv] * ql[1][i + iv] + xfn[i + iv] * ql[4][i + iv] + tmp1_v[iv] * qr[1][i + iv] + xfn[i + iv] * qr[4][i + iv];
            flux[2][i + iv] = tmp0_v[iv] * ql[2][i + iv] + yfn[i + iv] * ql[4][i + iv] + tmp1_v[iv] * qr[2][i + iv] + yfn[i + iv] * qr[4][i + iv];
            flux[3][i + iv] = tmp0_v[iv] * ql[3][i + iv] + zfn[i + iv] * ql[4][i + iv] + tmp1_v[iv] * qr[3][i + iv] + zfn[i + iv] * qr[4][i + iv];
            flux[4][i + iv] = ht_l_v[iv] * vn_l_v[iv] + ht_r_v[iv] * vn_r_v[iv];
            if (!steady)
                flux[4][i + iv] += (ql[4][i + iv] + qr[4][i + iv] + 2.0 * p_bar) * vgn[ni + iv]; // unsteady, 0.5?ú×?oó3????yμ?μ?·?

            // 2éó?roe???ù????μ￥?a??é?μ???àíá?
            tmp0_v[iv] = sqrt(qr[0][i + iv] / ql[0][i + iv]);
            tmp1_v[iv] = 1.0 / (1.0 + tmp0_v[iv]);
            rho_a_v[iv] = sqrt(qr[0][i + iv] * ql[0][i + iv]);
            u_a_v[iv] = (ql[1][i + iv] + qr[1][i + iv] * tmp0_v[iv]) * tmp1_v[iv];
            v_a_v[iv] = (ql[2][i + iv] + qr[2][i + iv] * tmp0_v[iv]) * tmp1_v[iv];
            w_a_v[iv] = (ql[3][i + iv] + qr[3][i + iv] * tmp0_v[iv]) * tmp1_v[iv];
            vn_a_v[iv] = u_a_v[iv] * xfn[i + iv] + v_a_v[iv] * yfn[i + iv] + w_a_v[iv] * zfn[i + iv];
            h_a_v[iv] = (ht_l_v[iv] / ql[0][i + iv] + ht_r_v[iv] / qr[0][i + iv] * tmp0_v[iv]) * tmp1_v[iv];
            q2_v[iv] = 0.5 * (u_a_v[iv] * u_a_v[iv] + v_a_v[iv] * v_a_v[iv] + w_a_v[iv] * w_a_v[iv]);
            c2_a_v[iv] = gamm1 * (h_a_v[iv] - q2_v[iv]);
            c2_a_v[iv] = fabs(c2_a_v[iv]);
            c_a_v[iv] = sqrt(c2_a_v[iv]);

            if (steady)
            {
                eigv1_v[iv] = fabs(vn_a_v[iv]);
                eigv2_v[iv] = fabs(vn_a_v[iv] + c_a_v[iv]);
                eigv3_v[iv] = fabs(vn_a_v[iv] - c_a_v[iv]);
            }
            else
            { // unsteady
                eigv1_v[iv] = fabs(vn_a_v[iv] - vgn[ni + iv]);
                eigv2_v[iv] = fabs(vn_a_v[iv] - vgn[ni + iv] + c_a_v[iv]);
                eigv3_v[iv] = fabs(vn_a_v[iv] - vgn[ni + iv] - c_a_v[iv]);
            }
            // Entropy fix
            if (EntropyCorType == 3)
            {
                epsa_r_v[iv] = alf_l;
            }
            else if (EntropyCorType == 4)
            {
                if (IsNormalFace[ni + iv] && IsShockFace[i + iv] == 0)
                {
                    epsa_r_v[iv] = 0.01 * alf_l;
                    // epsa_r = 0.0002;
                }
                else
                {
                    epsa_r_v[iv] = alf_l;
                }
            }
            else
            {
                mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
            }
            // cfl3d form
            if (steady)
            {
                spectral_v[iv] = fabs(u_a_v[iv]) + fabs(v_a_v[iv]) + fabs(w_a_v[iv]) + c_a_v[iv];
            }
            else
            {
                u_vgn_v[iv] = vgn[ni + iv] * xfn[i + iv];
                v_vgn_v[iv] = vgn[ni + iv] * yfn[i + iv];
                w_vgn_v[iv] = vgn[ni + iv] * zfn[i + iv];
                spectral_v[iv] = fabs(u_a_v[iv] - u_vgn_v[iv]) + fabs(v_a_v[iv] - v_vgn_v[iv]) + fabs(w_a_v[iv] - w_vgn_v[iv]) + c_a_v[iv];
            }
            epsaa_v[iv] = epsa_r_v[iv] * spectral_v[iv];
            epsbb_v[iv] = 0.25 / std::max(epsaa_v[iv], TINY);
            epscc_v[iv] = 2.0 * epsaa_v[iv];
            if (eigv1_v[iv] < epscc_v[iv])
                eigv1_v[iv] = eigv1_v[iv] * eigv1_v[iv] * epsbb_v[iv] + epsaa_v[iv];
            if (eigv2_v[iv] < epscc_v[iv])
                eigv2_v[iv] = eigv2_v[iv] * eigv2_v[iv] * epsbb_v[iv] + epsaa_v[iv];
            if (eigv3_v[iv] < epscc_v[iv])
                eigv3_v[iv] = eigv3_v[iv] * eigv3_v[iv] * epsbb_v[iv] + epsaa_v[iv];

            drho_v[iv] = qr[0][i + iv] - ql[0][i + iv];
            du_v[iv] = qr[1][i + iv] - ql[1][i + iv];
            dv_v[iv] = qr[2][i + iv] - ql[2][i + iv];
            dw_v[iv] = qr[3][i + iv] - ql[3][i + iv];
            dp_v[iv] = qr[4][i + iv] - ql[4][i + iv];
            dvn_v[iv] = vn_r_v[iv] - vn_l_v[iv];

            dq2_v[iv] = u_a_v[iv] * du_v[iv] + v_a_v[iv] * dv_v[iv] + w_a_v[iv] * dw_v[iv];

            tmp0_v[iv] = dp_v[iv] / c2_a_v[iv];
            tmp1_v[iv] = rho_a_v[iv] * dvn_v[iv] / c_a_v[iv];
            alpha1_v[iv] = (drho_v[iv] - tmp0_v[iv]) * eigv1_v[iv];
            alpha2_v[iv] = 0.5 * (tmp0_v[iv] + tmp1_v[iv]) * eigv2_v[iv];
            alpha3_v[iv] = 0.5 * (tmp0_v[iv] - tmp1_v[iv]) * eigv3_v[iv];

            tmp0_v[iv] = alpha1_v[iv] + alpha2_v[iv] + alpha3_v[iv];
            tmp1_v[iv] = eigv1_v[iv] * rho_a_v[iv];
            tmp2_v[iv] = -tmp1_v[iv] * dvn_v[iv] + (alpha2_v[iv] - alpha3_v[iv]) * c_a_v[iv];
            flux[0][i + iv] -= tmp0_v[iv];
            flux[1][i + iv] -= tmp0_v[iv] * u_a_v[iv] + tmp1_v[iv] * du_v[iv] + tmp2_v[iv] * xfn[i + iv];
            flux[2][i + iv] -= tmp0_v[iv] * v_a_v[iv] + tmp1_v[iv] * dv_v[iv] + tmp2_v[iv] * yfn[i + iv];
            flux[3][i + iv] -= tmp0_v[iv] * w_a_v[iv] + tmp1_v[iv] * dw_v[iv] + tmp2_v[iv] * zfn[i + iv];
            flux[4][i + iv] -= alpha1_v[iv] * q2_v[iv] + (alpha2_v[iv] + alpha3_v[iv]) * h_a_v[iv] + tmp1_v[iv] * dq2_v[iv] + tmp2_v[iv] * vn_a_v[iv];

            tmp0_v[iv] = 0.5 * area[i + iv];
            flux[0][i + iv] *= tmp0_v[iv];
            flux[1][i + iv] *= tmp0_v[iv];
            flux[2][i + iv] *= tmp0_v[iv];
            flux[3][i + iv] *= tmp0_v[iv];
            flux[4][i + iv] *= tmp0_v[iv];
        }
    }
    k = i;
    for (i = k; i < len; i++)
    {
        IntType ni;
        RealFlow rho_a, u_a, v_a, w_a, h_a, c_a, c2_a, vn_a, q2;
        RealFlow vn_l, et_l, ht_l, vn_r, et_r, ht_r;
        RealFlow tmp0, tmp1, tmp2, alpha1, alpha2, alpha3, eigv1, eigv2, eigv3;
        RealFlow drho, du, dv, dw, dp, dvn, dq2;
        RealGeom areax, areay, areaz;
        RealFlow spectral, epsaa, epsbb, epscc, epsa_r;
        RealFlow u_vgn, v_vgn, w_vgn;
        ni = ns + i;
        areax = xfn[i];
        areay = yfn[i];
        areaz = zfn[i];

        // Total energy
        et_l = (ql[4][i] + p_bar) / gamm1 + 0.5 * ql[0][i] *
                                                (ql[1][i] * ql[1][i] + ql[2][i] * ql[2][i] + ql[3][i] * ql[3][i]);
        et_r = (qr[4][i] + p_bar) / gamm1 + 0.5 * qr[0][i] *
                                                (qr[1][i] * qr[1][i] + qr[2][i] * qr[2][i] + qr[3][i] * qr[3][i]);
        ht_l = et_l + ql[4][i] + p_bar;
        ht_r = et_r + qr[4][i] + p_bar;

        // Full flux
        vn_l = areax * ql[1][i] + areay * ql[2][i] + areaz * ql[3][i];
        vn_r = areax * qr[1][i] + areay * qr[2][i] + areaz * qr[3][i];
        if (!steady)
        { // unsteady
            vn_l -= vgn[ni];
            vn_r -= vgn[ni];
        }

        tmp0 = vn_l * ql[0][i];
        tmp1 = vn_r * qr[0][i];
        flux[0][i] = tmp0 + tmp1;
        flux[1][i] = tmp0 * ql[1][i] + areax * ql[4][i] + tmp1 * qr[1][i] + areax * qr[4][i];
        flux[2][i] = tmp0 * ql[2][i] + areay * ql[4][i] + tmp1 * qr[2][i] + areay * qr[4][i];
        flux[3][i] = tmp0 * ql[3][i] + areaz * ql[4][i] + tmp1 * qr[3][i] + areaz * qr[4][i];
        flux[4][i] = ht_l * vn_l + ht_r * vn_r;
        if (!steady)
            flux[4][i] += (ql[4][i] + qr[4][i] + 2.0 * p_bar) * vgn[ni]; // unsteady, 0.5在最后乘面积的地方

        // 采用roe平均计算单元面上的物理量
        tmp0 = sqrt(qr[0][i] / ql[0][i]);
        tmp1 = 1.0 / (1.0 + tmp0);
        rho_a = sqrt(qr[0][i] * ql[0][i]);
        u_a = (ql[1][i] + qr[1][i] * tmp0) * tmp1;
        v_a = (ql[2][i] + qr[2][i] * tmp0) * tmp1;
        w_a = (ql[3][i] + qr[3][i] * tmp0) * tmp1;
        vn_a = u_a * areax + v_a * areay + w_a * areaz;
        h_a = (ht_l / ql[0][i] + ht_r / qr[0][i] * tmp0) * tmp1;

        q2 = 0.5 * (u_a * u_a + v_a * v_a + w_a * w_a);
        c2_a = gamm1 * (h_a - q2);
        c2_a = fabs(c2_a);
        c_a = sqrt(c2_a);

        if (steady)
        {
            eigv1 = fabs(vn_a);
            eigv2 = fabs(vn_a + c_a);
            eigv3 = fabs(vn_a - c_a);
        }
        else
        { // unsteady
            eigv1 = fabs(vn_a - vgn[ns + i]);
            eigv2 = fabs(vn_a - vgn[ns + i] + c_a);
            eigv3 = fabs(vn_a - vgn[ns + i] - c_a);
        }

        // Entropy fix
        if (EntropyCorType == 3)
        {
            epsa_r = alf_l;
        }
        else if (EntropyCorType == 4)
        {
            if (IsNormalFace[ni] && IsShockFace[i] == 0)
            {
                epsa_r = 0.01 * alf_l;
                // epsa_r = 0.0002;
            }
            else
            {
                epsa_r = alf_l;
            }
        }
        else
        {
            mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
        }

        // cfl3d form
        if (steady)
        {
            spectral = fabs(u_a) + fabs(v_a) + fabs(w_a) + c_a;
        }
        else
        {
            u_vgn = vgn[ni] * xfn[i];
            v_vgn = vgn[ni] * yfn[i];
            w_vgn = vgn[ni] * zfn[i];
            spectral = fabs(u_a - u_vgn) + fabs(v_a - v_vgn) + fabs(w_a - w_vgn) + c_a;
        }
        epsaa = epsa_r * spectral;
        epsbb = 0.25 / std::max(epsaa, TINY);
        epscc = 2.0 * epsaa;
        if (eigv1 < epscc)
            eigv1 = eigv1 * eigv1 * epsbb + epsaa;
        if (eigv2 < epscc)
            eigv2 = eigv2 * eigv2 * epsbb + epsaa;
        if (eigv3 < epscc)
            eigv3 = eigv3 * eigv3 * epsbb + epsaa;

        drho = qr[0][i] - ql[0][i];
        du = qr[1][i] - ql[1][i];
        dv = qr[2][i] - ql[2][i];
        dw = qr[3][i] - ql[3][i];
        dp = qr[4][i] - ql[4][i];
        dvn = vn_r - vn_l;

        dq2 = u_a * du + v_a * dv + w_a * dw;

        tmp0 = dp / c2_a;
        tmp1 = rho_a * dvn / c_a;
        alpha1 = (drho - tmp0) * eigv1;
        alpha2 = 0.5 * (tmp0 + tmp1) * eigv2;
        alpha3 = 0.5 * (tmp0 - tmp1) * eigv3;

        tmp0 = alpha1 + alpha2 + alpha3;
        tmp1 = eigv1 * rho_a;
        tmp2 = -tmp1 * dvn + (alpha2 - alpha3) * c_a;
        flux[0][i] -= tmp0;
        flux[1][i] -= tmp0 * u_a + tmp1 * du + tmp2 * areax;
        flux[2][i] -= tmp0 * v_a + tmp1 * dv + tmp2 * areay;
        flux[3][i] -= tmp0 * w_a + tmp1 * dw + tmp2 * areaz;
        flux[4][i] -= alpha1 * q2 + (alpha2 + alpha3) * h_a + tmp1 * dq2 + tmp2 * vn_a;

        tmp0 = 0.5 * area[i];
        flux[0][i] *= tmp0;
        flux[1][i] *= tmp0;
        flux[2][i] *= tmp0;
        flux[3][i] *= tmp0;
        flux[4][i] *= tmp0;
    }
*/
#else
    //!< not containing SIMD
#ifdef MF_OPENMP // only OpenMP
#pragma omp parallel for
#endif
    //!< else: serial code
    for (IntType i = 0; i < len; i++)
    {
        IntType ni;
        RealFlow rho_a, u_a, v_a, w_a, h_a, c_a, c2_a, vn_a, q2;
        RealFlow vn_l, et_l, ht_l, vn_r, et_r, ht_r;
        RealFlow tmp0, tmp1, tmp2, alpha1, alpha2, alpha3, eigv1, eigv2, eigv3;
        RealFlow drho, du, dv, dw, dp, dvn, dq2;
        RealGeom areax, areay, areaz;
        RealFlow spectral, epsaa, epsbb, epscc, epsa_r;
        // RealFlow u_vgn, v_vgn, w_vgn;
        ni = ns + i;
        areax = xfn[i];
        areay = yfn[i];
        areaz = zfn[i];

        //!< Total energy
        et_l = (ql[4][i] + p_bar) / gamm1 + 0.5 * ql[0][i] *
                                                (ql[1][i] * ql[1][i] + ql[2][i] * ql[2][i] + ql[3][i] * ql[3][i]);
        et_r = (qr[4][i] + p_bar) / gamm1 + 0.5 * qr[0][i] *
                                                (qr[1][i] * qr[1][i] + qr[2][i] * qr[2][i] + qr[3][i] * qr[3][i]);
        ht_l = et_l + ql[4][i] + p_bar;
        ht_r = et_r + qr[4][i] + p_bar;

        //!< Full flux
        vn_l = areax * ql[1][i] + areay * ql[2][i] + areaz * ql[3][i];
        vn_r = areax * qr[1][i] + areay * qr[2][i] + areaz * qr[3][i];
        // if (!steady) {   //!< unsteady
        //     vn_l -= vgn[ni];
        //     vn_r -= vgn[ni];
        // }

        tmp0 = vn_l * ql[0][i];
        tmp1 = vn_r * qr[0][i];
        flux[0][i] = tmp0 + tmp1;
        flux[1][i] = tmp0 * ql[1][i] + areax * ql[4][i] + tmp1 * qr[1][i] + areax * qr[4][i];
        flux[2][i] = tmp0 * ql[2][i] + areay * ql[4][i] + tmp1 * qr[2][i] + areay * qr[4][i];
        flux[3][i] = tmp0 * ql[3][i] + areaz * ql[4][i] + tmp1 * qr[3][i] + areaz * qr[4][i];
        flux[4][i] = ht_l * vn_l + ht_r * vn_r;
        // if (!steady) flux[4][i] += (ql[4][i] + qr[4][i] + 2.0 * p_bar) * vgn[ni];   ///< unsteady, 0.5 在最后乘面积的地方

        //!< 采用 roe 平均计算单元面上的物理量
        tmp0 = sqrt(qr[0][i] / ql[0][i]);
        tmp1 = 1.0 / (1.0 + tmp0);
        rho_a = sqrt(qr[0][i] * ql[0][i]);
        u_a = (ql[1][i] + qr[1][i] * tmp0) * tmp1;
        v_a = (ql[2][i] + qr[2][i] * tmp0) * tmp1;
        w_a = (ql[3][i] + qr[3][i] * tmp0) * tmp1;
        vn_a = u_a * areax + v_a * areay + w_a * areaz;
        h_a = (ht_l / ql[0][i] + ht_r / qr[0][i] * tmp0) * tmp1;

        q2 = 0.5 * (u_a * u_a + v_a * v_a + w_a * w_a);
        c2_a = gamm1 * (h_a - q2);
        c2_a = fabs(c2_a);
        c_a = sqrt(c2_a);

        // if (steady)
        // {
        eigv1 = fabs(vn_a);
        eigv2 = fabs(vn_a + c_a);
        eigv3 = fabs(vn_a - c_a);
        // }
        // else
        // {   //!< unsteady
        //     eigv1 = fabs(vn_a - vgn[ns + i]);
        //     eigv2 = fabs(vn_a - vgn[ns + i] + c_a);
        //     eigv3 = fabs(vn_a - vgn[ns + i] - c_a);
        // }

        //!< Entropy fix
        //!< Wisces: 暂时简化为 EntropyCorType == 3
        epsa_r = alf_l;
        /**
        if (EntropyCorType == 3)
        {
            epsa_r = alf_l;
        }
        else if (EntropyCorType == 4)
        {
            if (IsNormalFace[ni] && IsShockFace[i] == 0)
            {
                epsa_r = 0.01 * alf_l;
                // epsa_r = 0.0002;
            }
            else
            {
                epsa_r = alf_l;
            }
        }
        else
        {
            // mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
            exit(1);
        }
        */

        //!< cfl3d form
        //!< Wisces: steady == 1
        spectral = fabs(u_a) + fabs(v_a) + fabs(w_a) + c_a;
        /**
        if (steady)
        {
            spectral = fabs(u_a) + fabs(v_a) + fabs(w_a) + c_a;
        }
        else
        {
            u_vgn = vgn[ni] * xfn[i];
            v_vgn = vgn[ni] * yfn[i];
            w_vgn = vgn[ni] * zfn[i];
            spectral = fabs(u_a - u_vgn) + fabs(v_a - v_vgn) + fabs(w_a - w_vgn) + c_a;
        }
        */
        epsaa = epsa_r * spectral;
        epsbb = 0.25 / std::max(epsaa, TINY);
        epscc = 2.0 * epsaa;
        if (eigv1 < epscc)
            eigv1 = eigv1 * eigv1 * epsbb + epsaa;
        if (eigv2 < epscc)
            eigv2 = eigv2 * eigv2 * epsbb + epsaa;
        if (eigv3 < epscc)
            eigv3 = eigv3 * eigv3 * epsbb + epsaa;

        drho = qr[0][i] - ql[0][i];
        du = qr[1][i] - ql[1][i];
        dv = qr[2][i] - ql[2][i];
        dw = qr[3][i] - ql[3][i];
        dp = qr[4][i] - ql[4][i];
        dvn = vn_r - vn_l;

        dq2 = u_a * du + v_a * dv + w_a * dw;

        tmp0 = dp / c2_a;
        tmp1 = rho_a * dvn / c_a;
        alpha1 = (drho - tmp0) * eigv1;
        alpha2 = 0.5 * (tmp0 + tmp1) * eigv2;
        alpha3 = 0.5 * (tmp0 - tmp1) * eigv3;

        tmp0 = alpha1 + alpha2 + alpha3;
        tmp1 = eigv1 * rho_a;
        tmp2 = -tmp1 * dvn + (alpha2 - alpha3) * c_a;
        flux[0][i] -= tmp0;
        flux[1][i] -= tmp0 * u_a + tmp1 * du + tmp2 * areax;
        flux[2][i] -= tmp0 * v_a + tmp1 * dv + tmp2 * areay;
        flux[3][i] -= tmp0 * w_a + tmp1 * dw + tmp2 * areaz;
        flux[4][i] -= alpha1 * q2 + (alpha2 + alpha3) * h_a + tmp1 * dq2 + tmp2 * vn_a;

        tmp0 = 0.5 * area[i];
        flux[0][i] *= tmp0;
        flux[1][i] *= tmp0;
        flux[2][i] *= tmp0;
        flux[3][i] *= tmp0;
        flux[4][i] *= tmp0;
    }
#endif ///< ~ MF_SIMD / MF_OPENMP

    /**
    if (EntropyCorType == 4)
    {
        sdel_array_1D(IsShockFace);
    }
    */
}

/*!
 * @brief       Update residuals in cell with the fluxes at cell faces in 3D
 * @param       grid
 * @param       flux
 * @param       ns
 * @param       ne
 * @remarks     modify according to the fun [void LoadFlux(PolyGrid *grid, RealFlow *flux[], IntType ns, IntType ne)]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-05-31
 */
void LoadFlux(PolyGrid *grid, RealFlow *flux[], IntType ns, IntType ne)
{
    IntType face, i, c1, c2, count, nMid;
    IntType nTCell = grid->GetNTCell();
    IntType nBFace = grid->GetNBFace();
    IntType *f2c = grid->Getf2c();

    //!< Determine if there are boundary faces.
    nMid = ns;
    if (ne <= nBFace)
    {
        //!< If all boundary faces
        nMid = ne;
    }
    else if (ns < nBFace)
    {
        //!< Part of them are boundary faces
        nMid = nBFace;
    }

    //!< Get Residual
    RealFlow *res[5];
    res[0] = grid->GetRes();
    // RealFlow *_res[5];
    // _res[0] = (RealFlow *)res; // (RealFlow *)grid->GetDataPtr(REAL_FLOW, 5 * nTCell, "res");
    res[1] = &res[0][nTCell];
    res[2] = &res[1][nTCell];
    res[3] = &res[2][nTCell];
    res[4] = &res[3][nTCell];

//!< Group color openmp
#if (defined MF_OPENMP) && (defined OMP_GroupColor)
    if (grid->GroupColorSuccess)
    {
        IntType nIFace = grid->GetNIFace();
        IntType groupSize = grid->groupSize;
        IntType n_bcolor, n_icolor;
        n_bcolor = grid->bfacegroup.size();
        n_icolor = grid->ifacegroup.size();

        //!< physical face
        for (i = 0; i < n_bcolor; i++)
        {
            if (!i)
                ns = 0;
            else
                ns = grid->bfacegroup[i - 1];
            nMid = grid->bfacegroup[i];
#pragma omp parallel for private(face, c1) schedule(static, groupSize)
            for (face = ns; face < nMid; face++)
            {
                c1 = f2c[2 * face];

                res[0][c1] -= flux[0][face];
                res[1][c1] -= flux[1][face];
                res[2][c1] -= flux[2][face];
                res[3][c1] -= flux[3][face];
                res[4][c1] -= flux[4][face];
            }
        }

//!< Interface
#ifdef MF_MPICH
        count = 2 * (nBFace - nIFace);
        for (face = nBFace - nIFace; face < nBFace; face++)
        {
            c1 = f2c[count++];
            count++;

            res[0][c1] -= flux[0][face];
            res[1][c1] -= flux[1][face];
            res[2][c1] -= flux[2][face];
            res[3][c1] -= flux[3][face];
            res[4][c1] -= flux[4][face];
        }
#endif

        //!< Interior faces
        for (i = 0; i < n_icolor; i++)
        {
            if (!i)
                nMid = nBFace;
            else
                nMid = grid->ifacegroup[i - 1];
            ne = grid->ifacegroup[i];
#pragma omp parallel for private(face, c1, c2) schedule(static, groupSize)
            for (face = nMid; face < ne; face++)
            {
                c1 = f2c[2 * face];
                c2 = f2c[2 * face + 1];

                res[0][c1] -= flux[0][face];
                res[1][c1] -= flux[1][face];
                res[2][c1] -= flux[2][face];
                res[3][c1] -= flux[3][face];
                res[4][c1] -= flux[4][face];

                res[0][c2] += flux[0][face];
                res[1][c2] += flux[1][face];
                res[2][c2] += flux[2][face];
                res[3][c2] += flux[3][face];
                res[4][c2] += flux[4][face];
            }
        }
    }
    else
    {
        // For boundary faces, remember c2 is ghost cell
        // Determine if there are boundary faces.
        count = 2 * ns;
        i = 0;
        for (face = ns; face < nMid; face++)
        {
            c1 = f2c[count++];
            count++;

            res[0][c1] -= flux[0][i];
            res[1][c1] -= flux[1][i];
            res[2][c1] -= flux[2][i];
            res[3][c1] -= flux[3][i];
            res[4][c1] -= flux[4][i];
            i++;
        }

        // Interior faces
        for (face = nMid; face < ne; face++)
        {
            c1 = f2c[count++];
            c2 = f2c[count++];

            res[0][c1] -= flux[0][i];
            res[1][c1] -= flux[1][i];
            res[2][c1] -= flux[2][i];
            res[3][c1] -= flux[3][i];
            res[4][c1] -= flux[4][i];

            res[0][c2] += flux[0][i];
            res[1][c2] += flux[1][i];
            res[2][c2] += flux[2][i];
            res[3][c2] += flux[3][i];
            res[4][c2] += flux[4][i];
            i++;
        }
    }

#elif (defined MF_OPENMP) && (defined OMP_FaceColor)
    IntType nIFace = grid->GetNIFace();
    // IntType nTFace = grid->GetNTFace();
    // IntType ifacenum = nTFace - nBFace;
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
            startFace = 0; ///< for ns > 0 && ns < (*grid).bfacegroup[0]
        }
        else
        {
            startFace = (*grid).bfacegroup[fcolor - 1];
        }
        endFace = (*grid).bfacegroup[fcolor];

#pragma omp parallel for
        for (IntType face = startFace; face < endFace; face++)
        {
            IntType c1, c2, count;
            IntType i;
            count = 2 * face;
            c1 = f2c[count];
            i = face - ns; ///< ns = 0

            res[0][c1] -= flux[0][i];
            res[1][c1] -= flux[1][i];
            res[2][c1] -= flux[2][i];
            res[3][c1] -= flux[3][i];
            res[4][c1] -= flux[4][i];
        }
    }

    //!< Interface:
#ifdef MF_MPICH
    for (IntType face = pfacenum; face < nBFace; face++)
    {
        IntType count = 2 * face;
        IntType c1 = f2c[count];
        IntType i = face - ns;
        res[0][c1] -= flux[0][i];
        res[1][c1] -= flux[1][i];
        res[2][c1] -= flux[2][i];
        res[3][c1] -= flux[3][i];
        res[4][c1] -= flux[4][i];
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
        for (IntType face = startFace; face < endFace; face++)
        {
            IntType c1, c2, count;
            IntType i;
            count = 2 * face;
            c1 = f2c[count];
            c2 = f2c[count + 1];
            i = face - ns;

            res[0][c1] -= flux[0][i];
            res[1][c1] -= flux[1][i];
            res[2][c1] -= flux[2][i];
            res[3][c1] -= flux[3][i];
            res[4][c1] -= flux[4][i];

            res[0][c2] += flux[0][i];
            res[1][c2] += flux[1][i];
            res[2][c2] += flux[2][i];
            res[3][c2] += flux[3][i];
            res[4][c2] += flux[4][i];
        }
    }
    // sdel_array_1D(grid_bfacegroup);
    // sdel_array_1D(grid_ifacegroup);

#elif (defined MF_OPENMP) && (defined OMP_Reduction)
    IntType *nFPC = CalnFPC(grid); ///< number of faces in each real cell
    IntType **C2F = CalC2F(grid);  ///< real cell to face connection
    IntType j;
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
                res[0][i] -= flux[0][face];
                res[1][i] -= flux[1][face];
                res[2][i] -= flux[2][face];
                res[3][i] -= flux[3][face];
                res[4][i] -= flux[4][face];
            }
            else if (i == c2)
            {
                res[0][i] += flux[0][face];
                res[1][i] += flux[1][face];
                res[2][i] += flux[2][face];
                res[3][i] += flux[3][face];
                res[4][i] += flux[4][face];
            }
            else
            {
                // mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
                cerr << "Error in func[LoadFlux()] of file[cal_flux.cpp]" << endl;
                exit(1);
            }
        }
    }

#elif (defined MF_OPENMP) && (defined OMP_DIVREP) // Division & replication
    IntType threads = grid->threads;
    IntType nTFace = grid->GetNTFace();
    IntType startFace, endFace, t, k;
    if (grid->DivRepSuccess)
    {
#pragma omp parallel for private(t, i, k, startFace, endFace, c1, c2, face)
        for (t = 0; t < threads; t++)
        {
            // Boundary faces
            startFace = grid->idx_pthreads_bface[t];
            endFace = grid->idx_pthreads_bface[t + 1];
            for (i = startFace; i < endFace; i++)
            {
                face = grid->id_division_bface[i];
                c1 = f2c[2 * face];

                res[0][c1] -= flux[0][face];
                res[1][c1] -= flux[1][face];
                res[2][c1] -= flux[2][face];
                res[3][c1] -= flux[3][face];
                res[4][c1] -= flux[4][face];
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
                c1 = f2c[2 * face];
                c2 = f2c[2 * face + 1];
                if (abs(k) < nTFace)
                {
                    res[0][c1] -= flux[0][face];
                    res[1][c1] -= flux[1][face];
                    res[2][c1] -= flux[2][face];
                    res[3][c1] -= flux[3][face];
                    res[4][c1] -= flux[4][face];

                    res[0][c2] += flux[0][face];
                    res[1][c2] += flux[1][face];
                    res[2][c2] += flux[2][face];
                    res[3][c2] += flux[3][face];
                    res[4][c2] += flux[4][face];
                }
                else
                {
                    if (k > 0)
                    {
                        res[0][c1] -= flux[0][face];
                        res[1][c1] -= flux[1][face];
                        res[2][c1] -= flux[2][face];
                        res[3][c1] -= flux[3][face];
                        res[4][c1] -= flux[4][face];
                    }
                    else
                    {
                        res[0][c2] += flux[0][face];
                        res[1][c2] += flux[1][face];
                        res[2][c2] += flux[2][face];
                        res[3][c2] += flux[3][face];
                        res[4][c2] += flux[4][face];
                    }
                }
            }
        }
    }

#else
    //!< For boundary faces, remember c2 is ghost cell
    count = 2 * ns;
    i = 0;
    for (face = ns; face < nMid; face++)
    {
        c1 = f2c[count++];
        count++;

        res[0][c1] -= flux[0][i];
        res[1][c1] -= flux[1][i];
        res[2][c1] -= flux[2][i];
        res[3][c1] -= flux[3][i];
        res[4][c1] -= flux[4][i];
        i++;
    }

    //!< Interior faces
    for (face = nMid; face < ne; face++)
    {
        c1 = f2c[count++];
        c2 = f2c[count++];

        res[0][c1] -= flux[0][i];
        res[1][c1] -= flux[1][i];
        res[2][c1] -= flux[2][i];
        res[3][c1] -= flux[3][i];
        res[4][c1] -= flux[4][i];

        res[0][c2] += flux[0][i];
        res[1][c2] += flux[1][i];
        res[2][c2] += flux[2][i];
        res[3][c2] += flux[3][i];
        res[4][c2] += flux[4][i];
        i++;
    }
#endif ///< ~ MF_OPENMP
}

#ifdef TDTREE

/*!
 * @brief
 * @param       grid
 * @param       limit
 * @param       level
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-12
 */
void InviscidFlux_TDTree(PolyGrid *grid, RealFlow **limit, IntType level)
{

    IntType i, len;
    IntType nTCell = grid->GetNTCell();
    IntType nBFace = grid->GetNBFace();
    IntType nTFace = grid->GetNTFace();
    IntType n = nTCell + nBFace;
    IntType *f2c = grid->Getf2c();

    // // Get parameters
    // IntType steady;
    // RealFlow gam, p_bar, alf_l, alf_n, disFact = 1.;
    // grid->GetData(&steady, INT, 1, "steady");
    // grid->GetData(&gam, REAL_FLOW, 1, "gam");
    // grid->GetData(&p_bar, REAL_FLOW, 1, "p_bar");
    // grid->GetData(&alf_l, REAL_FLOW, 1, "alf_l");
    // grid->GetData(&alf_n, REAL_FLOW, 1, "alf_n");
    // grid->GetData(&disFact, REAL_FLOW, 1, "disFact", 0);

    // // for overlap
    // IntType *face_act = NULL;

    // Allocate temporary memories for q
    RealFlow *q[5];

    // Get flow variables
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

    // const IntType kNVar = 5;
    // RealFlow **dqdx = NULL, **dqdy = NULL, **dqdz = NULL;
    // snew_array_1D(dqdx, kNVar, dmrfl);
    // snew_array_1D(dqdy, kNVar, dmrfl);
    // snew_array_1D(dqdz, kNVar, dmrfl);
    // dqdx[0] = static_cast<RealFlow *>(
    //     grid->GetDataPtr(REAL_FLOW, kNVar * n, "dqdx"));
    // dqdy[0] = static_cast<RealFlow *>(
    //     grid->GetDataPtr(REAL_FLOW, kNVar * n, "dqdy"));
    // dqdz[0] = static_cast<RealFlow *>(
    //     grid->GetDataPtr(REAL_FLOW, kNVar * n, "dqdz"));
    // for (IntType i = 1; i < kNVar; ++i)
    // {
    //     dqdx[i] = &dqdx[i - 1][n];
    //     dqdy[i] = &dqdy[i - 1][n];
    //     dqdz[i] = &dqdz[i - 1][n];
    // }

    // RealFlow gascon;
    // grid->GetData(&gascon, REAL_FLOW, 1, "gascon");
    // IntType EntropyCorType = 4;
    // grid->GetData(&EntropyCorType, INT, 1, "EntropyCorType");

    // IntType *IsNormalFace = 0;
    // if (EntropyCorType == 4)
    // {
    //     IsNormalFace = (IntType *)grid->GetDataPtr(INT, nTFace, "IsNormalFace");
    //     if (!IsNormalFace)
    //     {
    //         grid->FindNormalFace();
    //         IsNormalFace = (IntType *)grid->GetDataPtr(INT, nTFace, "IsNormalFace");
    //     }
    // }

    // Get Residual
    RealFlow *res[5];
    res[0] = grid->GetRes();
    // RealFlow *res = grid->GetRes();
    // RealFlow *_res[5];
    // _res[0] = (RealFlow *)res; // (RealFlow *)grid->GetDataPtr(REAL_FLOW, 5 * nTCell, "res");
    res[1] = &res[0][nTCell];
    res[2] = &res[1][nTCell];
    res[3] = &res[2][nTCell];
    res[4] = &res[3][nTCell];

    RealFlow mach00 = mach;
    // grid->GetData(&mach00, REAL_FLOW, 1, "mach");

    // qleonardo: build usrArgs
    // char *userArgs[18] = {(char *)grid, (char *)limit, (char *)dqdx, (char *)dqdy, (char *)dqdz, (char *)q, (char *)res,
    //                       (char *)face_act, (char *)IsNormalFace, (char *)&steady, (char *)&gam, (char *)&p_bar, (char *)&alf_l,
    //                       (char *)&alf_n, (char *)&disFact, (char *)&gascon, (char *)&EntropyCorType, (char *)&mach00};
    char *userArgs[18] = {(char *)grid, (char *)limit, NULL, NULL, NULL, (char *)q, (char *)res,
                          NULL, NULL, (char *)&steady, (char *)&gam, (char *)&p_bar, (char *)&alf_l,
                          (char *)&alf_n, NULL, (char *)&gascon, (char *)&EntropyCorType, (char *)&mach00};

    TDTree *TDTreeRoot = grid->GetTDTree();
    TDTreeRoot->task_traversal(InviscidFlux_Kernel, NULL, userArgs, Forward);

    // sdel_array_1D(dqdx);
    // sdel_array_1D(dqdy);
    // sdel_array_1D(dqdz);
}

/*!
 * @brief
 * @param       userArgs
 * @param       treeArgs
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-12
 */
void InviscidFlux_Kernel(char **userArgs, TDTreeArgs *treeArgs)
{
    // #pragma omp critical
    //     {
    IntType ns = treeArgs->firstFace;
    IntType ne = treeArgs->lastFace + 1;

    // if (ns >= ne)
    //     return;

    PolyGrid *grid = (PolyGrid *)userArgs[0];
    RealFlow **limit = (RealFlow **)userArgs[1];
    // RealFlow **dqdx = (RealFlow **)userArgs[2];
    // RealFlow **dqdy = (RealFlow **)userArgs[3];
    // RealFlow **dqdz = (RealFlow **)userArgs[4];
    RealFlow **q = (RealFlow **)userArgs[5];
    RealFlow **res = (RealFlow **)userArgs[6];
    // IntType *face_act = ((IntType *)userArgs[7]);
    // IntType *IsNormalFace = ((IntType *)userArgs[8]);
    IntType steady = *((IntType *)userArgs[9]);
    RealFlow gam = *((RealFlow *)userArgs[10]);
    RealFlow p_bar = *((RealFlow *)userArgs[11]);
    RealFlow alf_l = *((RealFlow *)userArgs[12]);
    RealFlow alf_n = *((RealFlow *)userArgs[13]);
    // RealFlow disFact = *((RealFlow *)userArgs[14]);
    RealFlow gascon = *((RealFlow *)userArgs[15]);
    IntType EntropyCorType = *((IntType *)userArgs[16]);

    IntType nTCell = grid->GetNTCell();
    IntType nBFace = grid->GetNBFace();
    IntType n = nTCell + nBFace;

    // Get metrics
    RealGeom *xfn = grid->GetXfn();
    RealGeom *yfn = grid->GetYfn();
    RealGeom *zfn = grid->GetZfn();
    RealGeom *area = grid->GetFaceArea();
    RealGeom *vgn = grid->GetFaceNormalVelocity();

    RealFlow *flux[5], *ql[5], *qr[5];

    IntType len = ne - ns;
    ql[0] = new RealFlow[5 * len];
    qr[0] = new RealFlow[5 * len];
    flux[0] = new RealFlow[5 * len];

    for (IntType i = 1; i < 5; i++)
    {
        ql[i] = &ql[i - 1][len];
        qr[i] = &qr[i - 1][len];
        flux[i] = &flux[i - 1][len];
    }

    // Get left variables and right variables
    SetQlQrWithQ(grid, q, ql, qr, ns, ne);
    // if (limit != NULL)
    // {
    //     CalcuQlQr(grid, ql, qr, limit, dqdx, dqdy, dqdz, ns, ne, userArgs);
    // }
    // ModQlQrBou(grid, ql, qr, ns, ne, userArgs);

    // RoeFlux_noprec(grid, ql, qr, flux, IsNormalFace, &xfn[ns], &yfn[ns], &zfn[ns], &area[ns],
    //                &face_act[ns], gam, p_bar, alf_l, alf_n, gascon, EntropyCorType, steady, ns, ne, userArgs);
    RoeFlux_noprec_TDTree(grid, ql, qr, flux, NULL, &xfn[ns], &yfn[ns], &zfn[ns], &area[ns],
                          NULL, gam, p_bar, alf_l, alf_n, gascon, EntropyCorType, steady, ns, ne, userArgs);

    // Load the fluxes to residuals
    LoadFlux_TDTree(grid, flux, res, ns, ne);

    delete[] ql[0];
    delete[] qr[0];
    delete[] flux[0];
    // }
}

/*!
 * @brief
 * @param       grid
 * @param       ql
 * @param       qr
 * @param       flux
 * @param       IsNormalFace
 * @param       xfn
 * @param       yfn
 * @param       zfn
 * @param       area
 * @param       face_act
 * @param       gam
 * @param       p_bar
 * @param       alf_l
 * @param       alf_n
 * @param       gascon
 * @param       EntropyCorType
 * @param       steady
 * @param       ns
 * @param       ne
 * @param       userArgs
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-12
 */
void RoeFlux_noprec_TDTree(PolyGrid *grid, RealFlow *ql[5], RealFlow *qr[5], RealFlow *flux[5], IntType *IsNormalFace,
                           RealGeom *xfn, RealGeom *yfn, RealGeom *zfn, RealGeom *area, IntType *face_act, RealFlow gam,
                           RealFlow p_bar, RealFlow alf_l, RealFlow alf_n, RealFlow gascon, IntType EntropyCorType,
                           IntType steady, IntType ns, IntType ne, char **userArgs)
{
    // register IntType i, c1, c2, count, ni;
    IntType len;
    // RealFlow rho_a, u_a, v_a, w_a, h_a, c_a, c2_a, vn_a, q2;
    // RealFlow vn_l, et_l, ht_l, vn_r, et_r, ht_r /*, gamm1*/;
    // RealFlow tmp0, tmp1, tmp2, alpha1, alpha2, alpha3, eigv1, eigv2, eigv3;
    // RealFlow drho, du, dv, dw, dp, dvn, dq2;
    // RealGeom areax, areay, areaz;
    // RealFlow spectral, epsaa, epsbb, epscc, epsa_r;
    // RealFlow u_vgn, v_vgn, w_vgn;

    IntType nTCell = grid->GetNTCell();
    IntType nBFace = grid->GetNBFace();
    IntType nTFace = grid->GetNTFace();
    IntType n = nTCell + nBFace;
    IntType *f2c = grid->Getf2c();

    RealGeom *vgn = grid->GetFaceNormalVelocity();

    // IntType *IsShockFace = 0;
    // if (EntropyCorType == 4)
    // {
    //     // shock face or not
    //     IsShockFace = new int[ne - ns]();
    //     CalIsShockFace(grid, IsShockFace, ns, ne, userArgs);
    // }

    // gamm1 = gam - 1.0;
    len = ne - ns;
    // count = 2 * ns;

    for (IntType i = 0; i < len; i++)
    {
        IntType ni;
        RealFlow rho_a, u_a, v_a, w_a, h_a, c_a, c2_a, vn_a, q2;
        RealFlow vn_l, et_l, ht_l, vn_r, et_r, ht_r;
        RealFlow tmp0, tmp1, tmp2, alpha1, alpha2, alpha3, eigv1, eigv2, eigv3;
        RealFlow drho, du, dv, dw, dp, dvn, dq2;
        RealGeom areax, areay, areaz;
        RealFlow spectral, epsaa, epsbb, epscc, epsa_r;
        // RealFlow u_vgn, v_vgn, w_vgn;

        ni = ns + i;
        // c1 = f2c[count++];
        // c2 = f2c[count++];
        areax = xfn[i];
        areay = yfn[i];
        areaz = zfn[i];

        // Total energy
        et_l = (ql[4][i] + p_bar) / gamm1 + 0.5 * ql[0][i] *
                                                (ql[1][i] * ql[1][i] + ql[2][i] * ql[2][i] + ql[3][i] * ql[3][i]);
        et_r = (qr[4][i] + p_bar) / gamm1 + 0.5 * qr[0][i] *
                                                (qr[1][i] * qr[1][i] + qr[2][i] * qr[2][i] + qr[3][i] * qr[3][i]);
        ht_l = et_l + ql[4][i] + p_bar;
        ht_r = et_r + qr[4][i] + p_bar;

        // Full flux
        vn_l = areax * ql[1][i] + areay * ql[2][i] + areaz * ql[3][i];
        vn_r = areax * qr[1][i] + areay * qr[2][i] + areaz * qr[3][i];
        // if (!steady)
        // { // unsteady
        //     vn_l -= vgn[ni];
        //     vn_r -= vgn[ni];
        // }

        tmp0 = vn_l * ql[0][i];
        tmp1 = vn_r * qr[0][i];
        flux[0][i] = tmp0 + tmp1;
        flux[1][i] = tmp0 * ql[1][i] + areax * ql[4][i] + tmp1 * qr[1][i] + areax * qr[4][i];
        flux[2][i] = tmp0 * ql[2][i] + areay * ql[4][i] + tmp1 * qr[2][i] + areay * qr[4][i];
        flux[3][i] = tmp0 * ql[3][i] + areaz * ql[4][i] + tmp1 * qr[3][i] + areaz * qr[4][i];
        flux[4][i] = ht_l * vn_l + ht_r * vn_r;
        // if (!steady)
        //     flux[4][i] += (ql[4][i] + qr[4][i] + 2.0 * p_bar) * vgn[ni]; // unsteady, 0.5ĺ¨ćĺäšé˘ç§Żçĺ°ć

        // 采用 roe 平均计算单元面上的物理量
        tmp0 = sqrt(qr[0][i] / ql[0][i]);
        tmp1 = 1.0 / (1.0 + tmp0);
        rho_a = sqrt(qr[0][i] * ql[0][i]);
        u_a = (ql[1][i] + qr[1][i] * tmp0) * tmp1;
        v_a = (ql[2][i] + qr[2][i] * tmp0) * tmp1;
        w_a = (ql[3][i] + qr[3][i] * tmp0) * tmp1;
        vn_a = u_a * areax + v_a * areay + w_a * areaz;
        h_a = (ht_l / ql[0][i] + ht_r / qr[0][i] * tmp0) * tmp1;

        q2 = 0.5 * (u_a * u_a + v_a * v_a + w_a * w_a);
        c2_a = gamm1 * (h_a - q2);
        c2_a = fabs(c2_a);
        c_a = sqrt(c2_a);

        // if (steady)
        // {
        eigv1 = fabs(vn_a);
        eigv2 = fabs(vn_a + c_a);
        eigv3 = fabs(vn_a - c_a);
        // }
        // else
        // { // unsteady
        //     eigv1 = fabs(vn_a - vgn[ns + i]);
        //     eigv2 = fabs(vn_a - vgn[ns + i] + c_a);
        //     eigv3 = fabs(vn_a - vgn[ns + i] - c_a);
        // }

        // Entropy fix
        epsa_r = alf_l;
        // if (EntropyCorType == 3)
        // {
        //     epsa_r = alf_l;
        // }
        // else if (EntropyCorType == 4)
        // {
        //     if (IsNormalFace[ni] && IsShockFace[i] == 0)
        //     {
        //         epsa_r = 0.01 * alf_l;
        //         // epsa_r = 0.0002;
        //     }
        //     else
        //     {
        //         epsa_r = alf_l;
        //     }
        // }
        // else
        // {
        //     mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
        // }

        // cfl3d form
        spectral = fabs(u_a) + fabs(v_a) + fabs(w_a) + c_a;
        // if (steady)
        // {
        //     spectral = fabs(u_a) + fabs(v_a) + fabs(w_a) + c_a;
        // }
        // else
        // {
        //     u_vgn = vgn[ni] * xfn[i];
        //     v_vgn = vgn[ni] * yfn[i];
        //     w_vgn = vgn[ni] * zfn[i];
        //     spectral = fabs(u_a - u_vgn) + fabs(v_a - v_vgn) + fabs(w_a - w_vgn) + c_a;
        // }
        epsaa = epsa_r * spectral;
        epsbb = 0.25 / std::max(epsaa, TINY);
        epscc = 2.0 * epsaa;
        if (eigv1 < epscc)
            eigv1 = eigv1 * eigv1 * epsbb + epsaa;
        if (eigv2 < epscc)
            eigv2 = eigv2 * eigv2 * epsbb + epsaa;
        if (eigv3 < epscc)
            eigv3 = eigv3 * eigv3 * epsbb + epsaa;

        drho = qr[0][i] - ql[0][i];
        du = qr[1][i] - ql[1][i];
        dv = qr[2][i] - ql[2][i];
        dw = qr[3][i] - ql[3][i];
        dp = qr[4][i] - ql[4][i];
        dvn = vn_r - vn_l;

        dq2 = u_a * du + v_a * dv + w_a * dw;

        tmp0 = dp / c2_a;
        tmp1 = rho_a * dvn / c_a;
        alpha1 = (drho - tmp0) * eigv1;
        alpha2 = 0.5 * (tmp0 + tmp1) * eigv2;
        alpha3 = 0.5 * (tmp0 - tmp1) * eigv3;

        tmp0 = alpha1 + alpha2 + alpha3;
        tmp1 = eigv1 * rho_a;
        tmp2 = -tmp1 * dvn + (alpha2 - alpha3) * c_a;
        flux[0][i] -= tmp0;
        flux[1][i] -= tmp0 * u_a + tmp1 * du + tmp2 * areax;
        flux[2][i] -= tmp0 * v_a + tmp1 * dv + tmp2 * areay;
        flux[3][i] -= tmp0 * w_a + tmp1 * dw + tmp2 * areaz;
        flux[4][i] -= alpha1 * q2 + (alpha2 + alpha3) * h_a + tmp1 * dq2 + tmp2 * vn_a;

        tmp0 = 0.5 * area[i];
        flux[0][i] *= tmp0;
        flux[1][i] *= tmp0;
        flux[2][i] *= tmp0;
        flux[3][i] *= tmp0;
        flux[4][i] *= tmp0;
    }
    // if (EntropyCorType == 4)
    // {
    //     delete[] IsShockFace;
    // }
}

/*!
 * @brief
 * @param       grid
 * @param       flux
 * @param       res
 * @param       ns
 * @param       ne
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-12
 */
void LoadFlux_TDTree(PolyGrid *grid, RealFlow *flux[], RealFlow **res, IntType ns, IntType ne)
{
    IntType face, i, c1, c2, count, nMid;
    IntType nTCell = grid->GetNTCell();
    IntType nBFace = grid->GetNBFace();
    IntType *f2c = grid->Getf2c();

    // Determine if there are boundary faces.
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

    // For boundary faces, remember c2 is ghost cell
    count = 2 * ns;
    i = 0;
    for (face = ns; face < nMid; face++)
    {
        c1 = f2c[count++];
        count++;

        res[0][c1] -= flux[0][i];
        res[1][c1] -= flux[1][i];
        res[2][c1] -= flux[2][i];
        res[3][c1] -= flux[3][i];
        res[4][c1] -= flux[4][i];
        i++;
    }

    // Interior faces
    for (face = nMid; face < ne; face++)
    {
        c1 = f2c[count++];
        c2 = f2c[count++];

        res[0][c1] -= flux[0][i];
        res[1][c1] -= flux[1][i];
        res[2][c1] -= flux[2][i];
        res[3][c1] -= flux[3][i];
        res[4][c1] -= flux[4][i];

        res[0][c2] += flux[0][i];
        res[1][c2] += flux[1][i];
        res[2][c2] += flux[2][i];
        res[3][c2] += flux[3][i];
        res[4][c2] += flux[4][i];
        i++;
    }
}

#endif