/*!
 * @file        constant.h
 * @brief       define some constants used in MiniFS
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @version     1.0
 * @date        2023-05-21
 * @copyright   Copyright (c) 2023  PERSONAL
 *
 * @par 修改日志:
 * <table>
 * <tr><th> Date        <th> Version  <th> Author  <th> Description
 * <tr><td> 2023-05-21  <td> 1.0      <td> Wisces  <td> 〔内容〕
 * <tr><td> 2023-06-20  <td> 2.0      <td> Wisces  <td> 添加了边界条件类型
 * </table>
 */

#ifndef MF_CONSTANT_H
#define MF_CONSTANT_H
#include "number_type.h"

// Constant type
const RealGeom PI = 3.14159265358979323846;
const IntType MAXLINE = 1024;

// Face variable reconstruction order
const IntType FIRST_ORDER = 1;
const IntType SECOND_ORDER = 2;
const IntType LIMITED_VENCAT = 4;

// Flow type
const IntType INVISCID = 0;
const IntType LAMINAR = 1;
const IntType S_A_MODEL = 2;

// Attribute of cell/node for overlap grid
const IntType INACTIVE = 0;
const IntType CHANGACTIVE = 3;

const IntType SEG_LEN = 10240;

/// linear system methods
const IntType LUSGS = 0;
const IntType GMRES = 1;

// Boundary condition type
const IntType WALL = 3;
const IntType SYMM = 4;
const IntType FAR_FIELD = 6;
const IntType FAR_FIELD_MOD = 61;
const IntType INTERFACE = 10;

// Simple function
#ifndef MIN
#define MIN(X, Y) std::min((X), (Y))
#endif
#ifndef MAX
#define MAX(X, Y) std::max((X), (Y))
#endif

// #ifndef SIN
// #define SIN(X) (sin(X * PI / 180.0))
// #endif
// #ifndef COS
// #define COS(X) (cos(X * PI / 180.0))
// #endif

//!< TEST
#ifndef PR
#define PR(x) cout << #x << " = " << x << endl;
#endif

#endif //~MF_CONSTANT_H
