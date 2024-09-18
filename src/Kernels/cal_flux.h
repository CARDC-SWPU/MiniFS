/*!
 * @file        cal_flux.h
 * @brief       无粘通量的主要函数接口
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @version     1.0
 * @date        2023-06-10
 * @copyright   Copyright (c) 2023  PERSONAL
 *
 * @par 修改日志:
 * <table>
 * <tr><th> Date        <th> Version  <th> Author  <th> Description
 * <tr><td> 2023-06-10  <td> 1.0      <td> Wisces  <td> 〔内容〕
 * <tr><td> 2023-06-15  <td> 2.0      <td> Wisces  <td> 采用全局变量形式存储数据
 * </table>
 */

#ifndef MF_CAL_FLUX_H
#define MF_CAL_FLUX_H

#include "grid_polyhedra.h"

//!< 残差置零（或在第一次迭代时初始化）
void ZeroResiduals(PolyGrid *grid);

//!< 计算无粘通量
void InviscidFlux(PolyGrid *grid, RealFlow **limit, IntType level);
void SetQlQrWithQ(PolyGrid *grid, RealFlow *q[], RealFlow *ql[], RealFlow *qr[], IntType ns, IntType ne);
void CompInvFlux(PolyGrid *grid, RealFlow *ql[5], RealFlow *qr[5], RealFlow *flux[5],
                 RealGeom *xfn, RealGeom *yfn, RealGeom *zfn, RealGeom *area,
                 RealGeom *vgn, IntType *face_act, RealFlow gam, RealFlow p_bar,
                 RealFlow alf_l, RealFlow alf_n, IntType type_flux,
                 IntType ns, IntType ne);
void RoeFlux_noprec(PolyGrid *grid, RealFlow *ql[5], RealFlow *qr[5], RealFlow *flux[5],
                    RealGeom *xfn, RealGeom *yfn, RealGeom *zfn, RealGeom *area, IntType *face_act,
                    RealFlow gam, RealFlow p_bar, RealFlow alf_l, RealFlow alf_n,
                    IntType ns, IntType ne);
void LoadFlux(PolyGrid *grid, RealFlow *flux[], IntType ns, IntType ne);

#ifdef TDTREE
void InviscidFlux_TDTree(PolyGrid *grid, RealFlow **limit, IntType level);
void InviscidFlux_Kernel(char **userArgs, TDTreeArgs *treeArgs);
void RoeFlux_noprec_TDTree(PolyGrid *grid, RealFlow *ql[5], RealFlow *qr[5], RealFlow *flux[5], IntType *IsNormalFace,
                           RealGeom *xfn, RealGeom *yfn, RealGeom *zfn, RealGeom *area, IntType *face_act, RealFlow gam,
                           RealFlow p_bar, RealFlow alf_l, RealFlow alf_n, RealFlow gascon, IntType EntropyCorType,
                           IntType steady, IntType ns, IntType ne, char **userArgs);
void LoadFlux_TDTree(PolyGrid *grid, RealFlow *flux[], RealFlow **res, IntType ns, IntType ne);
#endif

#endif ///< ~MF_CAL_FLUX_H