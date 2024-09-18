/*!
 * @file        cal_gradient.h
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
 * <tr><td> 2023-07-07 <td> 1.0      <td> Wisces  <td> 〔内容〕
 * </table>
 */

#ifndef MF_CAL_GRADIENT_H
#define MF_CAL_GRADIENT_H

#include "grid_polyhedra.h"

void AllocAndCalQuantityGradient(PolyGrid *grid);
void CalculateGradient(PolyGrid *grid /*, RealFlow *dqdx[5], RealFlow *dqdy[5], RealFlow *dqdz[5]*/);
// void CompGradientQ(PolyGrid *grid, RealFlow *q, RealFlow *dqdx, RealFlow *dqdy, RealFlow *dqdz, IntType name);
void CompGradientQ_Gauss_Node(PolyGrid *grid, RealFlow *q, RealFlow *dqdx, RealFlow *dqdy, RealFlow *dqdz, IntType name, RealFlow *u_n, RealFlow *v_n, RealFlow *w_n);
void CompNodeVar3D_dist(PolyGrid *grid, RealFlow *q_n, RealFlow *q, IntType name, RealFlow *u_n, RealFlow *v_n, RealFlow *w_n);

#ifdef TDTREE
void CompGradientQ_Gauss_Node_TDTree(PolyGrid *grid, RealFlow **q, RealFlow **dqdx, RealFlow **dqdy, RealFlow **dqdz, IntType ns, IntType ne);
void CompNodeVar3D_dist_TDTree(PolyGrid *grid, RealFlow **q_n, RealFlow **q, RealFlow **u_n, RealFlow **v_n, RealFlow **w_n, IntType ns, IntType ne);
void CompNodeVar3D_dist_Kernel2(char **userArgs, TDTreeArgs *treeArgs);
void CompGradientQ_Gauss_Node_Kernel3(char **userArgs, TDTreeArgs *treeArgs);
void CompGradientQ_Gauss_Node_Kernel4(char **userArgs, TDTreeArgs *treeArgs);
#endif

#endif ///< ~MF_CAL_GRADIENT_H