/*!
 * @file        load_grid.h
 * @brief       加载网格的主要函数接口
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @version     1.0
 * @date        2023-05-27
 * @copyright   Copyright (c) 2023  PERSONAL
 *
 * @par 修改日志:
 * <table>
 * <tr><th> Date        <th> Version  <th> Author  <th> Description
 * <tr><td> 2023-05-27  <td> 1.0      <td> Wisces  <td> 〔内容〕
 * <tr><td> 2023-06-13  <td> 1.1      <td> Wisces  <td> 添加了 InitParameter() 函数
 * <tr><td> 2023-06-15  <td> 2.0      <td> Wisces  <td> 删除了 InitParameter()、ComputeMachineZero() 函数
 * <tr><td> 2023-06-19  <td> 3.0      <td> Wisces  <td> 添加了 SpecifyNeighbors()、SpecifyBC()、PassInterfaceData() 和 SetGhostVariables() 函数
 * <tr><td> 2023-07-03  <td> 4.0      <td> Wisces  <td> 添加了 Decoupling()、ReorderCellforLUSGS() 函数
 * </table>
 */

#ifndef MF_LOAD_GRID_H
#define MF_LOAD_GRID_H

#include "grid_polyhedra.h"

//!< Load original grid
void LoadGrid(PolyGrid *grid, BCond *bc);
// void ReadMFlowGrid_binary(PolyGrid **grids, const IntType mg, const string &filename, ExtraGridData &extra_grid_data);
void ReadMFlowGrid_binary(PolyGrid *grid, const string &filename /*, ExtraGridData &extra_grid_data*/);
void ReadMFlowBoundaryCondition(BCond *bc, const char *filename);

//!< Do some preprocessing on grid data for subsequent solution calculations
void PreprocessGrid(PolyGrid *grid, BCond *bc);
void SpecifyNeighbors(PolyGrid *grid);
void CommGraph(PolyGrid *grid, IntType &nNeighbor, IntType *&nb);
void CommGraph_node(PolyGrid *grid, IntType &nNeighborN, IntType *&nbN);
void SpecifyBC(PolyGrid *grid, BCond *bc);
void AllocAndInitFlowfield(PolyGrid *grid);
void PassInterfaceData(PolyGrid *grid);
void SetGhostVariables(PolyGrid *grid);

//!< 设计实现了基于循环级与任务级两种共享内存并行编程模型的总共四种并行算法
void Decoupling(PolyGrid *grid);
//!< 针对 LUSGS，进行网格排序（内含 CellColoring）
void ReorderCellforLUSGS(PolyGrid *grid);

//!< 工具函数：int 转化为 string 类型函数
string int2str(const IntType source);

#endif // ~MF_ALGOMMFLOW_H