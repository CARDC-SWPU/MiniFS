/*!
 * @file        parallel_base.cpp
 * @brief       和 MPI 并行相关的函数
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @version     1.0
 * @date        2023-06-18
 * @copyright   Copyright (c) 2023  PERSONAL
 *
 * @par 修改日志:
 * <table>
 * <tr><th> Date        <th> Version  <th> Author  <th> Description
 * <tr><td> 2023-06-18  <td> 1.0      <td> Wisces  <td> 〔内容〕
 * </table>
 */

//!< direct head file
#include "parallel_base.h"

//!< user difined head file
#include "number_type.h"
#include "grid_base.h"
#include "memory_util.h"

//!< head file relying on condition-compiling
#ifdef MF_MPICH
#include <mpi.h>
#endif

#ifdef MF_MPICH
int myid = 0;
int numprocs = 1;
int myZone = 1;    ///< start from 1 (= myid + 1)
MPI_Comm GridComm; ///< for each grid
#endif

namespace Parallel
{
#ifdef MF_MPICH

    /*!
     * @brief       Initialize MPI parallel enviroment
     * @param       argc
     * @param       argv
     *
     * @author      Wisces 〔wwwangqs17@163.com〕
     * @date        2023-06-17
     */
    void InitMPI(int argc, char *argv[])
    {
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);

        myZone = myid + 1;

        // // Initialize to global communication world, tangj add
        // // For overlap case, GridComm will be reset to correct
        // // value, so the ranks used for one same grid will be
        // // assigned to one same sub-communication.
        GridComm = MPI_COMM_WORLD;
    }

    /*!
     * @brief       Finalize MPI parallel enviroment
     *
     * @author      Wisces 〔wwwangqs17@163.com〕
     * @date        2023-06-17
     */
    void FinalMPI()
    {
        MPI_Finalize();
    }

    /*!
     * @brief       Sum the data among all parallel processors in global_comm_world
     * @param       data
     * @param       global_comm_world
     *
     * @author      Wisces 〔wwwangqs17@163.com〕
     * @date        2023-06-24
     */
    void parallel_sum(IntType &data, const MPI_Comm global_comm_world)
    {
        IntType data_init = data;
        MPI_Allreduce(&data_init, &data, 1, MPIIntType, MPI_SUM, global_comm_world);
    }

    /*!
     * @brief       Sum each the item of data among all parallel processors in global_comm_world
     * @param       data
     * @param       n_data
     * @param       global_comm_world
     *
     * @author      Wisces 〔wwwangqs17@163.com〕
     * @date        2023-06-24
     */
    void parallel_sum(IntType *data, IntType n_data, const MPI_Comm global_comm_world)
    {
        IntType *data_init = NULL;
        snew_array_1D(data_init, n_data);
        for (IntType n = 0; n < n_data; ++n)
            data_init[n] = data[n];
        MPI_Allreduce(data_init, data, n_data, MPIIntType, MPI_SUM, global_comm_world);
        sdel_array_1D(data_init);
    }

    /*!
     * @brief       Sum the data among all parallel processors in global_comm_world
     * @param       data
     * @param       global_comm_world
     *
     * @author      Wisces 〔wwwangqs17@163.com〕
     * @date        2023-06-29
     */
    void parallel_sum(RealFlow &data, const MPI_Comm global_comm_world)
    {
        RealFlow data_init = data;
        MPI_Allreduce(&data_init, &data, 1, MPIReal, MPI_SUM, global_comm_world);
    }

    /*!
     * @brief       Find the minimum and maximum data among all parallel processors in global_comm_world
     * @param       min_data
     * @param       max_data
     * @param       global_comm_world
     *
     * @author      Wisces 〔wwwangqs17@163.com〕
     * @date        2023-06-24
     */
    void parallel_min_max(RealFlow &min_data, RealFlow &max_data, const MPI_Comm global_comm_world)
    {
        RealFlow min_max[2] = {min_data, -max_data};
        RealFlow min_max_glb[2] = {min_max[0], min_max[1]};
        MPI_Allreduce(min_max, min_max_glb, 2, MPIReal, MPI_MIN, global_comm_world);
        min_data = min_max_glb[0];
        max_data = -min_max_glb[1];
    }

    /*!
     * @brief       ind the maximum data among all parallel processors in global_comm_world
     * @param       data
     * @param       n_data
     * @param       global_comm_world
     *
     * @author      Wisces 〔wwwangqs17@163.com〕
     * @date        2023-06-29
     */
    void parallel_max(RealFlow *data, IntType n_data, const MPI_Comm global_comm_world)
    {
        RealFlow *data_init = NULL;
        snew_array_1D(data_init, n_data);
        for (IntType n = 0; n < n_data; ++n)
            data_init[n] = data[n];
        MPI_Allreduce(data_init, data, n_data, MPIReal, MPI_MAX, global_comm_world);
        sdel_array_1D(data_init);
    }

#endif

} ///< ~namespace Parallel