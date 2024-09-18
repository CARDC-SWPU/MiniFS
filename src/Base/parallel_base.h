/*!
 * @file        parallel_base.h
 * @brief       和 MPI 并行相关的函数接口
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @version     1.0
 * @date        2023-06-17
 * @copyright   Copyright (c) 2023  PERSONAL
 *
 * @par 修改日志:
 * <table>
 * <tr><th> Date        <th> Version  <th> Author  <th> Description
 * <tr><td> 2023-06-17  <td> 1.0      <td> Wisces  <td> 〔内容〕
 * </table>
 */

#ifndef MF_PARALLEL_BASE_H
#define MF_PARALLEL_BASE_H

//!< user defined head files
#include "number_type.h"

//!< head file relying on condition-compiling
#ifdef MF_MPICH
#include <mpi.h>
#endif

namespace Parallel
{
#ifdef MF_MPICH
    //!< Initialize MPI parallel enviroment
    void InitMPI(int argc, char *argv[]);

    //!< Finalize MPI parallel enviroment
    void FinalMPI();

    //!< Sum the data among all parallel processors in global_comm_world
    void parallel_sum(IntType &data, const MPI_Comm global_comm_world);

    //!< Sum each the item of data among all parallel processors in global_comm_world
    void parallel_sum(IntType *data, IntType n_data, const MPI_Comm global_comm_world);

    //!< Sum the data among all parallel processors in global_comm_world
    void parallel_sum(RealFlow &data, const MPI_Comm global_comm_world);

    //!< Find the minimum and maximum data among all parallel processors in global_comm_world
    void parallel_min_max(RealFlow &min_data, RealFlow &max_data, const MPI_Comm global_comm_world);

    //!< Find the maximum data among all parallel processors in global_comm_world
    void parallel_max(RealFlow *data, IntType n_data, const MPI_Comm global_comm_world);

#endif ///< ~MF_MPICH
} ///< ~namespace Parallel

#endif ///< ~MF_PARALLEL_BASE_H