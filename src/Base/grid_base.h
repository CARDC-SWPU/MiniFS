/*!
 * @file        grid_base.h
 * @brief       Abstract base grid object
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
 * </table>
 */

#ifndef MF_GRID_BASE_H
#define MF_GRID_BASE_H

//!< C++ build-in head files
#include <iostream>

//!< user defined head files
#include "number_type.h"

/*!
 * @brief       Base Grid class
 *
 */
class Grid
{

private:
    IntType nTNode;      ///< no of total nodes
    RealGeom *x, *y, *z; ///< coordinates
    // DataStore *gField;   ///< store the field variables, anything
    // DataSafe *gPara;     ///< store the control parameters
    // IntType zn;          ///< the zone number

public:
    virtual void ComputeMetrics() = 0;

    // IntType GetZone() const;
    // void SetZone(const IntType znin);
    IntType GetNTNode() const;
    void SetNTNode(const IntType ntn);
    void SetX(RealGeom *xin);
    void SetY(RealGeom *yin);
    void SetZ(RealGeom *zin);
    RealGeom *GetX() const;
    RealGeom *GetY() const;
    RealGeom *GetZ() const;

    Grid();
    // explicit Grid(IntType in);
    virtual ~Grid();
};

//!< inline functions of class Grid
inline IntType Grid::GetNTNode() const
{
    return nTNode;
}

inline void Grid::SetNTNode(const IntType ntn)
{
    nTNode = ntn;
}

inline void Grid::SetX(RealGeom *xin)
{
    if (x != NULL && x != xin)
        delete[] x;
    x = xin;
    // sSet(x,xin);
}

inline void Grid::SetY(RealGeom *yin)
{
    if (y != NULL && y != yin)
        delete[] y;
    y = yin;
    // sSet(y, yin);
}

inline void Grid::SetZ(RealGeom *zin)
{
    if (z != NULL && z != zin)
        delete[] z;
    z = zin;
    // sSet(z, zin);
}

inline RealGeom *Grid::GetX() const
{
    return x;
}

inline RealGeom *Grid::GetY() const
{
    return y;
}

inline RealGeom *Grid::GetZ() const
{
    return z;
}

#endif //~MF_GRID_BASE_H
