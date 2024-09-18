/*!
 * @file        grid_base.cpp
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

//!< direct head file
#include "grid_base.h"

//!< C++ build-in head files
#include <iostream>

//!< user defined head files
#include "memory_util.h"

using namespace std;

Grid::Grid() : nTNode(0)
{
#ifdef MF_DEBUG
    cout << endl
         << "Constructor [Grid()] is called!" << endl;
#endif
    x = NULL;
    y = NULL;
    z = NULL;
    // gPara = NULL;
    // gField = NULL;
}

// Grid::Grid(IntType in) : nTNode(0) , zn(in)
// {
//     x = NULL;
//     y = NULL;
//     z = NULL;
//     // gPara = NULL;
//     // gField = NULL;
// }

Grid::~Grid()
{
    // if (x != NULL)
    // {
    //     delete[] x;
    //     x = NULL;
    // }
    // if (y != NULL)
    // {
    //     delete[] y;
    //     y = NULL;
    // }
    // if (z != NULL)
    // {
    //     delete[] z;
    //     z = NULL;
    // }

    sdel_array_1D(x);
    sdel_array_1D(y);
    sdel_array_1D(z);
}
