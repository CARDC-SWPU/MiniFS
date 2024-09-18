/*!
 * @file        para_field_global.cpp
 * @brief       全局参数和流场变量的定义
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @version     1.0
 * @date        2023-06-15
 * @copyright   Copyright (c) 2023  PERSONAL
 *
 * @par 修改日志:
 * <table>
 * <tr><th> Date        <th> Version  <th> Author  <th> Description
 * <tr><td> 2023-06-15  <td> 1.0      <td> Wisces  <td> 〔内容〕
 * <tr><td> 2023-06-20  <td> 2.0      <td> Wisces  <td> 添加了边界相关参数
 * <tr><td> 2023-07-06  <td> 3.0      <td> Wisces  <td> 添加了网格重排序相关参数
 * </table>
 */

//!< direct head file
#include "para_field_global.h"

//!< C++ build-in head files
#include <iostream>
#include <cmath>

//!< user defined head files
#include "constant.h"
#include "boundary_condition.h" ///< 边界相关参数
#include "grid_polyhedra.h"     ///< 边界类型声明

/*!
 * @brief       Parameters for simulation
 *
 */

//!< 气动力计算类型
IntType steady = 1; ///< steady or unsteady (0 -- unsteady; 1 -- steady)

//!< Parameters for iterations
IntType Unst_steps = 1; ///< total physical steps of iteration
IntType n_steps = 20;

/*!
 * @brief       Some input conditions for the computation
 *
 */

//!< CFL number, the global cfl number is compute using cfl_start, cfl_end and cfl_nstep
RealFlow cfl_start = 100.0; ///< for m6wing_12w / m6wing_370w
RealFlow cfl_end = 100.0;
// RealFlow cfl_start = 20.0; ///< for CHNT1_600w
// RealFlow cfl_end = 20.0;
IntType cfl_nstep = 2000;

//!< func-setted: Some parameters which are setted in func [void Zone::FixParameter()] and [void Zone::UpdateParameter()]
//!< func-setted: CFL number and relatives
RealFlow cfl_coeff = 0.5;           ///< cfl_coeff：粗网格 cfl 系数，与细网格的 cfl 数相乘得到粗网格的 cfl 数，建议取 0.5
RealFlow cfl_min = 0.5 * cfl_start; ///< cfl_min used for minish cfl number where p is mini. cfl_min=0.5*cfl_start is suggested.
RealFlow cfl_ratio = cfl_end / cfl_start;

//!< Incoming flow condition type
// IntType IncomingType = 1; ///< 1--Re and T; 2--p and T; 3--Altitude(unit:km,0~91)
RealFlow re = 1.8140449177e7; ///< for m6wing_12w / m6wing_370w
RealFlow T = 288.15;
// RealFlow re = 17040000.0; ///< for CHNT1_600w
// RealFlow T = 300.0;
// RealFlow p_bar = 93994.7; ///< 后面计算得最终值

//!< Mach number, Angle of attack and Angle of slide
RealFlow mach = 0.8395;
RealFlow alpha = 3.06;
RealFlow beita = 0.0;
RealFlow theta = 0.0;

/*!
 * @brief       Some parameters for the equation
 *
 */
//!< viscous model: 0--inviscid; 1--laminar; 2--SA model; 3--SST model
IntType vis_mode = INVISCID;

/*!
 * @brief       几何多重相关参数
 *
 */

/*!
 * @brief       Some parameters for the Invisflux computations
 *
 */
///< Entropy fix for Roe scheme
IntType EntropyCorType = 3; ///< = 3 (entropy correction is original harten's); = 4 (new modified method(4 is suggested))
RealFlow epsa_r = 0.1;      ///< entropy correction constants(0.0~0.3), 0.0 is no entropy correction(0.1 is suggested)

//!< func-setted - entropy correction constants for Roe flux
RealFlow alf_l = epsa_r; ///< 0.3 is suggested for alf_l, as 0.3 is suggested for alf_n, 0.0 is no entropy correction
RealFlow alf_n = epsa_r;

//!< order of spacial scheme
IntType order = 1; ///< 1--1st order; 2--2nd order no limiter; 3--2nd order(Barths); 4--2nd order(Venkatakrishnan)

/*!
 * @brief       Some parameters for the time iteration computations
 *
 */
IntType sweeps = 1;

//!< func-setted: limit for DQ in lu-sgs
IntType DQ_limit = 2; ///< 2 is suggested for non-precondition, as 1 is used for precondition

//!< func-setted: ratio for lhs of the lusgs iteration
RealFlow lhs_omga = 1.0; ///< 1.0 is suggested

//!< func-setted: limit max_dt for max_dt/min_dt<=ratio_dtmax
RealFlow ratio_dtmax = 1.0e80; ///< approach to no limit

/*!
 * @brief       Some parameters for the GMRES method
 *
 */
//!< mattype: NONE, MATRIX_FREE, MATRIX_FREE_PETSC, BLOCK_MATRIX_PETSC, BLOCK_MATRIX_LI
IntType mattype = 0;

/*!
 * @brief       Some parameters for the unsteady computation
 *
 */
//!< real time step for each time step of dual-time method
RealFlow real_dt = 0.010;

//!< func-setted: time_accuracy
RealFlow time_accuracy = 0.5; ///< 0 is the first order , 0.5 is the second order

/*!
 * @brief       Some parameters for the air
 *
 */
//!< func-setted: some constants for air(Do Not modify)
RealFlow gam = 1.4;
RealFlow gascon = 287.053;
// RealFlow prl = 0.72;
// RealFlow prt = 0.9;
// RealFlow cp = 1003.0;

/*!
 * @brief       Some parameters which are setted through simply calculation (in func [void Zone::UpdateParameter()])
 *
 */

RealFlow tref = 288.15;       ///< 中间值
RealFlow sref = 110.4;        ///< 中间值
RealFlow amuref = 1.78938e-5; ///< 中间值
RealFlow trat = (T / tref);   ///< 中间值
RealFlow amu = amuref * trat * sqrt(trat) * (tref + sref) / (T + sref);
RealFlow ainf = sqrt(gam * gascon * T);
RealFlow uqq = ainf * mach; ///< 中间值

RealFlow alpha_t = (alpha * PI / 180.);                 ///< 中间值
RealFlow beita_t = (beita * PI / 180.);                 ///< 中间值
RealFlow u_s = uqq * cos(alpha_t) * cos(-1. * beita_t); ///< 标量 u，区别于数组 *u
RealFlow v_s = uqq * sin(-1. * beita_t);                ///< 标量 v，区别于数组 *v
RealFlow w_s = uqq * sin(alpha_t) * cos(-1. * beita_t); ///< 标量 w，区别于数组 *w
RealFlow p_s = 0.00;                                    ///< 标量 p，区别于数组 *p

//!< for m6wing_12w
// T = 288.15
// re = 1.81404e7
// u_s = 285.27;
// v_s = 0;
// w_s = 15.2499
// p_s = 0
// rho_s = 1.12625
// amu = 1.78938e-5
// ainf = 340.294
// p_bar = 93984.5
// p_stag = 55120.5
// e_stag = 281327
// iexp = 16
// p_min = -93975.1
// p_max = 1.39707e6
// rho_min = 0.000113625
// rho_max = 15.7995
// p_break = -93044.7
// e_stag_max = 2.81327e10

//!< if(IncomingType == 1)
RealFlow rho_s = re * amu / uqq; ///< rho_s 表示标量 rho，区别于数组 *rho
RealFlow p_bar = rho_s * gascon * T;

RealFlow gamm1 = gam - 1.0;                         ///< 中间值
RealFlow temp = 1.0 + 0.5 * gamm1 * mach * mach;    ///< 中间值
RealFlow p_stag_t = pow(temp, gam / gamm1);         ///< 中间值
RealFlow rho_stag = pow(temp, 1.0 / gamm1) * rho_s; ///< 中间值
RealFlow p_stag = (p_stag_t - 1.0) * p_bar;

RealFlow e_stag = p_bar / gamm1 + 0.5 * rho_s * (u_s * u_s + v_s * v_s + w_s * w_s);

IntType iexp = 16; ///< [void ComputeMachineZero(PolyGrid *grid)]

RealFlow ratio_rhop_max = 10.0;   ///< 中间值
RealFlow ratio_rhop_min = 1.0e-4; ///< 中间值
RealFlow ratio_p_break = 0.01;    ///< 中间值

RealFlow p_min = ratio_rhop_min * p_bar - p_bar;
RealFlow p_max = ratio_rhop_max * (p_stag + p_bar) - p_bar;
RealFlow rho_min = ratio_rhop_min * rho_s;
RealFlow rho_max = ratio_rhop_max * rho_stag;
RealFlow p_break = ratio_p_break * p_bar - p_bar;
RealFlow e_stag_max = 1.0e5 * e_stag;

/*!
 * @brief       其他标量参数（对应于 gField）
 */

IntType iter_done = 0; ///< 用于统计计算次数（内外层循环总次数），在计算时间步的时候需要

/*!
 * @brief       其他标量参数（对应于 gPara）
 */
RealFlow dt_max = 0.0;
RealFlow dt_min = BIG;

/*!
 * @brief       Some parameters for the Gradient computations（改为 grid 成员变量）
 */
// IntType GaussLayer = 5; ///< if GradQ or GradQTurb = 6, then in boundary layer, the grid layer number of GaussLayer will use Gauss method
//                         ///< if GaussLayer <= 0, then all use Gauss node.
// IntType *CellLayerNo = NULL;

// IntType *nodesymm = NULL;
// RealGeom *xfn_n_symm = NULL;
// RealGeom *yfn_n_symm = NULL;
// RealGeom *zfn_n_symm = NULL;

/*!
 * @brief       其他数组参数（对应于 gField）（改为 grid 成员变量）
 */

// RealGeom *norm_dist_c2c = NULL; ///< store normal disstance of two cells（used for viscous spectral radii in LUSGS）

//!< 网格重排序
// IntType *LUSGSLayer = NULL;
// IntType *LUSGSCellOrder = NULL;
// IntType *LUSGScellsPerlayer = NULL; ///< Used for CellColor

//!< flow variables
// RealFlow *rho = NULL;
// RealFlow *u = NULL;
// RealFlow *v = NULL;
// RealFlow *w = NULL;
// RealFlow *p = NULL;

//!< 梯度计算
// RealFlow *dqdx = NULL;
// RealFlow *dqdy = NULL;
// RealFlow *dqdz = NULL;

//!< time step
// RealFlow *dt = NULL; ///< dt_timestep

//!< residuals
// RealFlow *res = NULL;

//!< DQ
// RealFlow *DQ = NULL;

/*!
 * @brief       其他数组参数（对应于 gPara）
 */

/*!
 * @brief       边界相关参数（input.para）
 */
//!< for m6wing_12w
// BC
// {
//     int n_patch_groups = 3;
//     {
//         int ids[] = {1};
//         int type = 6; // FAR_FIELD
//     }
//     {
//         int ids[] = {2};
//         int type = 4; // SYMM
//     }
//     {
//         int ids[] = {3};
//         int type = 3; // WALL;
//     }
// }
// IntType n_patch = 3;
// BCRecord *bcr1 = new BCRecord(1, FAR_FIELD, "far_field");
// BCRecord *bcr2 = new BCRecord(2, SYMM, "symm");
// BCRecord *bcr3 = new BCRecord(3, WALL, "wall");
// BCond *bc = new BCond();

//!< for m6wing_25w
// BC
// {
//     int n_patch_groups = 3;
//     {
//         int ids[] = {1};
//         int type = 6; // FAR_FIELD
//     }
//     {
//         int ids[] = {2};
//         int type = 4; // SYMM
//     }
//     {
//         int ids[] = {3:5};
//         int type = 3; // WALL;
//     }
// }
// IntType n_patch = 5;
// BCRecord *bcr1 = new BCRecord(1, FAR_FIELD, "far_field");
// BCRecord *bcr2 = new BCRecord(2, SYMM, "symm");
// BCRecord *bcr3 = new BCRecord(3, WALL, "wall");
// BCRecord *bcr4 = new BCRecord(4, WALL, "wall");
// BCRecord *bcr5 = new BCRecord(5, WALL, "wall");
// BCond *bc = new BCond();

// !< for CHNT1_600w
// BC
// {
//     int n_patch_groups = 3;
//     {
//         int ids[] = {6};
//         int type = 4; // SYMM
//     }
//     {
//         int ids[] = {3};
//         int type = 6; // FAR_FIELD
//     }
//     {
//         int ids[] = {1, 2, 4, 5, 7, 8};
//         int type = 3; // WALL;
//     }
// }
// IntType n_patch = 8;
// BCRecord *bcr1 = new BCRecord(1, WALL, "wall");
// BCRecord *bcr2 = new BCRecord(2, WALL, "wall");
// BCRecord *bcr3 = new BCRecord(3, FAR_FIELD, "far_field");
// BCRecord *bcr4 = new BCRecord(4, WALL, "wall");
// BCRecord *bcr5 = new BCRecord(5, WALL, "wall");
// BCRecord *bcr6 = new BCRecord(6, SYMM, "symm");
// BCRecord *bcr7 = new BCRecord(7, WALL, "wall");
// BCRecord *bcr8 = new BCRecord(8, WALL, "wall");
// BCond *bc = new BCond();

/*!
 * @brief       计时相关参数
 *
 */
double total_t[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
