/*!
 * @file        forward_lusgs.h
 * @brief       LUSGS 求解的主要函数接口
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
 * </table>
 */

#ifndef MF_FORWARD_LUSGS_H
#define MF_FORWARD_LUSGS_H

#include "grid_polyhedra.h"

//!< Forward the solution one time step using LU-SGS method.
void ForwardLUSGS(PolyGrid *grid, IntType level);
void CalDiagLUSGS(PolyGrid *grid, RealFlow *Diag, IntType level);
void SolveLUSGS3D(PolyGrid *grid, RealFlow *Diag, RealFlow *DQ[5], IntType *nFPC, IntType **C2F, IntType level);
void FluxLUSGS3D(RealFlow flux[5], RealFlow q[5], RealFlow DQ[5], RealGeom fa_n[3], RealFlow gam, RealFlow p_bar, RealFlow lhs_omga);
void UpdateFlowField3D_CFL3d(PolyGrid *grid, RealFlow *DQ[5]);

#ifdef TDTREE
void CalDiagLUSGS_TDTree(PolyGrid *grid, RealFlow *Diag, IntType level);
void CalDiagLUSGS_Kernel(char **userArgs, TDTreeArgs *treeArgs);
void SolveLUSGS3D_TDTree(PolyGrid *grid, RealFlow *Diag, RealFlow *DQ[5], IntType *nFPC, IntType **C2F, IntType level);
void SolveLUSGS3D_forward_Kernel(char **userArgs, TDTreeArgs *treeArgs);
void SolveLUSGS3D_backward_Kernel(char **userArgs, TDTreeArgs *treeArgs);
#endif

#endif ///< ~MF_FORWARD_LUSGS_H