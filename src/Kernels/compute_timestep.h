/*!
 * @file        compute_timestep.h
 * @brief       计算时间步的函数接口
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
 * <tr><td> 2023-06-15  <td> 2.0      <td> Wisces  <td> 采用全局变量形式存储数据
 * </table>
 */

#ifndef MF_COMPUTE_TIMESTEP_H
#define MF_COMPUTE_TIMESTEP_H

#include "grid_polyhedra.h"

//!< Compute time step
void ComputeTimeStep(PolyGrid *grid);
void TimeStepNormal_new(PolyGrid *grid, RealFlow *dt, IntType vis_run);
void LimitTimeStep(PolyGrid *grid, RealFlow *dt);
void CellIsMG(PolyGrid *grid, IntType *det);

#ifdef TDTREE
void TimeStepNormal_new_TDTree(PolyGrid *grid, RealFlow *dt, IntType vis_run);
void TimeStepNormal_new_Kernel(char **userArgs, TDTreeArgs *treeArgs);
#endif

#endif ///< ~MF_COMPUTE_TIMESTEP_H