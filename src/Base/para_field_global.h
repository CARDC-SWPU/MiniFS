/*!
 * @file        para_field_global.h
 * @brief       全局参数和流场变量的声明
 * @note        全局变量需要使用 extern 关键字进行声明，并在 .cpp 文件中进行定义
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

#ifndef MF_PARA_FIELD_GLOBAL_H
#define MF_PARA_FIELD_GLOBAL_H

#include "number_type.h"
#include "boundary_condition.h"

/*!
 * @brief       Parameters for simulation
 *
 */

//!< 气动力计算类型
extern IntType steady; ///< steady or unsteady (0 -- unsteady; 1 -- steady)

//!< Parameters for iterations
extern IntType Unst_steps; ///< total physical steps of iteration
extern IntType n_steps;

/*!
 * @brief       Some input conditions for the computation
 *
 */

//!< CFL number, the global cfl number is compute using cfl_start, cfl_end and cfl_nstep
extern RealFlow cfl_start;
extern RealFlow cfl_end;
extern IntType cfl_nstep;

//!< func-setted: Some parameters which are setted in func [void Zone::FixParameter()] and [void Zone::UpdateParameter()]
//!< func-setted: CFL number and relatives
extern RealFlow cfl_coeff; ///< cfl_coeff：粗网格 cfl 系数，与细网格的 cfl 数相乘得到粗网格的 cfl 数，建议取 0.5
extern RealFlow cfl_min;   ///< cfl_min used for minish cfl number where p is mini. cfl_min=0.5*cfl_start is suggested.
extern RealFlow cfl_ratio;

//!< Incoming flow condition type
// extern IntType IncomingType = 1; ///< 1--Re and T; 2--p and T; 3--Altitude(unit:km,0~91)
extern RealFlow re;
extern RealFlow T;
// extern RealFlow p_bar = 93994.7; ///< 后面计算得最终值

//!< Mach number, Angle of attack and Angle of slide
extern RealFlow mach;
extern RealFlow alpha;
extern RealFlow beita;
extern RealFlow theta;

/*!
 * @brief       Some parameters for the equation
 *
 */
//!< viscous model: 0--inviscid; 1--laminar; 2--SA model; 3--SST model
extern IntType vis_mode;

/*!
 * @brief       几何多重相关参数
 *
 */

/*!
 * @brief       Some parameters for the Invisflux computations
 *
 */
///< Entropy fix for Roe scheme
extern IntType EntropyCorType; ///< = 3 (entropy correction is original harten's); = 4 (new modified method(4 is suggested))
extern RealFlow epsa_r;        ///< entropy correction constants(0.0~0.3), 0.0 is no entropy correction(0.1 is suggested)

//!< func-setted - entropy correction constants for Roe flux
extern RealFlow alf_l; ///< 0.3 is suggested for alf_l, as 0.3 is suggested for alf_n, 0.0 is no entropy correction
extern RealFlow alf_n;

//!< order of spacial scheme
extern IntType order; ///< 1--1st order; 2--2nd order no limiter; 3--2nd order(Barths); 4--2nd order(Venkatakrishnan)

/*!
 * @brief       Some parameters for the time iteration computations
 *
 */
extern IntType sweeps;

//!< func-setted: limit for DQ in lu-sgs
extern IntType DQ_limit; ///< 2 is suggested for non-precondition, as 1 is used for precondition

//!< func-setted: ratio for lhs of the lusgs iteration
extern RealFlow lhs_omga; ///< 1.0 is suggested

//!< func-setted: limit max_dt for max_dt/min_dt<=ratio_dtmax
extern RealFlow ratio_dtmax; ///< approach to no limit

/*!
 * @brief       Some parameters for the GMRES method
 *
 */
//!< mattype: NONE, MATRIX_FREE, MATRIX_FREE_PETSC, BLOCK_MATRIX_PETSC, BLOCK_MATRIX_LI
extern IntType mattype;

/*!
 * @brief       Some parameters for the unsteady computation
 *
 */
//!< real time step for each time step of dual-time method
extern RealFlow real_dt;

//!< func-setted: time_accuracy
extern RealFlow time_accuracy; ///< 0 is the first order , 0.5 is the second order

/*!
 * @brief       Some parameters for the air
 *
 */
//!< func-setted: some constants for air(Do Not modify)
extern RealFlow gam;
extern RealFlow gascon;
// extern RealFlow prl = 0.72;
// extern RealFlow prt = 0.9;
// extern RealFlow cp = 1003.0;

/*!
 * @brief       Some parameters which are setted through simply calculation (in func [void Zone::UpdateParameter()])
 *
 */

// extern RealFlow tref;   ///< 中间值
// extern RealFlow sref;   ///< 中间值
// extern RealFlow amuref; ///< 中间值
// extern RealFlow trat;   ///< 中间值
extern RealFlow amu;
extern RealFlow ainf;
extern RealFlow uqq;

// extern RealFlow alpha_t; ///< 中间值
// extern RealFlow beita_t; ///< 中间值
extern RealFlow u_s; ///< 标量 u，区别于数组 *u
extern RealFlow v_s; ///< 标量 v，区别于数组 *v
extern RealFlow w_s; ///< 标量 w，区别于数组 *w
extern RealFlow p_s; ///< 标量 p，区别于数组 *p

//!< if(IncomingType == 1)
extern RealFlow rho_s; ///< rho_s 表示标量 rho，区别于数组 *rho
extern RealFlow p_bar;

extern RealFlow gamm1;
// extern RealFlow temp;     ///< 中间值
// extern RealFlow p_stag_t; ///< 中间值
// extern RealFlow rho_stag; ///< 中间值
extern RealFlow p_stag;

extern RealFlow e_stag;

extern IntType iexp; ///< [void ComputeMachineZero(PolyGrid *grid)]

// extern RealFlow ratio_rhop_max; ///< 中间值
// extern RealFlow ratio_rhop_min; ///< 中间值
// extern RealFlow ratio_p_break;  ///< 中间值

extern RealFlow p_min;
extern RealFlow p_max;
extern RealFlow rho_min;
extern RealFlow rho_max;
extern RealFlow p_break;
extern RealFlow e_stag_max;

/*!
 * @brief       其他标量参数（对应于 gField）
 */

extern IntType iter_done; ///< 用于统计计算次数（内外层循环总次数），在计算时间步的时候需要

/*!
 * @brief       其他标量参数（对应于 gPara）
 */
extern RealFlow dt_max;
extern RealFlow dt_min;

/*!
 * @brief       Some parameters for the Gradient computations（改为 grid 成员变量）
 */
// extern IntType GaussLayer;
// extern IntType *CellLayerNo;
// extern IntType *nodesymm;
// extern RealGeom *xfn_n_symm;
// extern RealGeom *yfn_n_symm;
// extern RealGeom *zfn_n_symm;

/*!
 * @brief       其他数组参数（对应于 gField）（改为 grid 成员变量）
 */

// extern RealGeom *norm_dist_c2c; ///< store normal disstance of two cells（used for viscous spectral radii in LUSGS）

//!< 网格重排序
// extern IntType *LUSGSLayer;
// extern IntType *LUSGSCellOrder;
// extern IntType *LUSGScellsPerlayer;

//!< flow variables
// extern RealFlow *rho;
// extern RealFlow *u;
// extern RealFlow *v;
// extern RealFlow *w;
// extern RealFlow *p;

//!< 梯度计算
// extern RealFlow *dqdx;
// extern RealFlow *dqdy;
// extern RealFlow *dqdz;

//!< time step
// extern RealFlow *dt; ///< dt_timestep

//!< residuals
// extern RealFlow *res;

//!< DQ
// extern RealFlow *DQ;

/*!
 * @brief       其他数组参数（对应于 gPara）
 */

/*!
 * @brief       边界相关参数（input.para）
 */
//!< for m6wing_12w
// extern IntType n_patch;
// extern BCRecord *bcr1;
// extern BCRecord *bcr2;
// extern BCRecord *bcr3;
// extern BCond *bc;

//!< for m6wing_25w
// extern IntType n_patch;
// extern BCRecord *bcr1;
// extern BCRecord *bcr2;
// extern BCRecord *bcr3;
// extern BCRecord *bcr4;
// extern BCRecord *bcr5;
// extern BCond *bc;

//!< for CHNT1_600w
// extern BCRecord *bcr1;
// extern BCRecord *bcr2;
// extern BCRecord *bcr3;
// extern BCRecord *bcr4;
// extern BCRecord *bcr5;
// extern BCRecord *bcr6;
// extern BCRecord *bcr7;
// extern BCRecord *bcr8;
// extern BCond *bc;

/*!
 * @brief       计时相关参数
 *
 */
extern double total_t[5]; ///< inviscid、lusgs、gradient、timestep、comm

#endif //~MF_PARA_FIELD_GLOBAL_H