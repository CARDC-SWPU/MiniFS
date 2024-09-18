#!/bin/bash
g++ -std=c++11 -O3 -o output/main_serial src/main.cpp src/Base/*.cpp src/Kernels/*.cpp -I src/Base/ -I src/Kernels/ \
-DMF_TIMING \
-DMF_DEBUG \

# -DMF_MPICH \

# -DMF_OPENMP \
    # -DOMP_CellColor \
    # -DOMP_FaceColor \
    # -DOMP_GroupColor \
    # -DOMP_Reduction \
    # -DOMP_DIVREP \
# -DTDTREE

# -DMF_DEBUG \
# -DNDEBUG \

## 串行
    g++ -std=c++11 -O3 -o output/main_serial src/main.cpp src/Base/*.cpp src/Kernels/*.cpp  -I src/Base/ -I src/Kernels/   -DMF_TIMING

## OpenMP 并行
    g++ -std=c++11 -O3 -fopenmp -o output/main_facecoloring src/main.cpp src/Base/*.cpp src/Kernels/*.cpp  -I src/Base/ -I src/Kernels/   -DMF_TIMING -DMF_OPENMP -DOMP_CellColor -DOMP_FaceColor

    g++ -std=c++11 -O3 -fopenmp -o output/main_groupcoloring src/main.cpp src/Base/*.cpp src/Kernels/*.cpp  -I src/Base/ -I src/Kernels/   -DMF_TIMING -DMF_OPENMP -DOMP_CellColor -DOMP_GroupColor

    g++ -std=c++11 -O3 -fopenmp -o output/main_reduction src/main.cpp src/Base/*.cpp src/Kernels/*.cpp  -I src/Base/ -I src/Kernels/   -DMF_TIMING -DMF_OPENMP -DOMP_CellColor -DOMP_Reduction

    g++ -std=c++11 -O3 -fopenmp -o output/main_divrep src/main.cpp src/Base/*.cpp src/Kernels/*.cpp -L src/Lib/ -lmetis  -I src/Base/ -I src/Kernels/   -I src/Lib/ -DMF_TIMING -DMF_OPENMP -DOMP_CellColor -DOMP_DIVREP

    g++ -std=c++11 -O3 -fopenmp -o output/main_task src/main.cpp src/Base/*.cpp src/Kernels/*.cpp src/Lib/libTDTree.so -L src/Lib/ -lmetis  -I src/Base/ -I src/Kernels/   -I src/Lib/ -DMF_TIMING -DTDTREE

## MPI 并行

    mpicxx -std=c++11 -O3 -o output/main_mpi src/main.cpp src/Base/*.cpp src/Kernels/*.cpp  -I src/Base/ -I src/Kernels/   -DMF_TIMING -DMF_MPICH

    mpicxx -std=c++11 -O3 -o output/main_gaspi src/main.cpp src/Base/*.cpp src/Kernels/*.cpp /public5/home/sch8427/GPI-2-next/install/lib64/libGPI2.so -I/public5/home/sch8427/GPI-2-next/install/include/  -I src/Base/ -I src/Kernels/   -DMF_TIMING -DMF_MPICH -DGASPI

## MPI+OpenMP 混合并行

     mpicxx -std=c++11 -O3 -fopenmp -o output/main_mpi_facecoloring src/main.cpp src/Base/*.cpp src/Kernels/*.cpp  -I src/Base/ -I src/Kernels/   -DMF_TIMING -DMF_MPICH -DMF_OPENMP -DOMP_CellColor -DOMP_FaceColor

    mpicxx -std=c++11 -O3 -fopenmp -o output/main_mpi_groupcoloring src/main.cpp src/Base/*.cpp src/Kernels/*.cpp  -I src/Base/ -I src/Kernels/   -DMF_TIMING -DMF_MPICH -DMF_OPENMP -DOMP_CellColor -DOMP_GroupColor

    mpicxx -std=c++11 -O3 -fopenmp -o output/main_mpi_reduction src/main.cpp src/Base/*.cpp src/Kernels/*.cpp  -I src/Base/ -I src/Kernels/   -DMF_TIMING -DMF_MPICH -DMF_OPENMP -DOMP_CellColor -DOMP_Reduction

    mpicxx -std=c++11 -O3 -fopenmp -o output/main_mpi_divrep src/main.cpp src/Base/*.cpp src/Kernels/*.cpp -L src/Lib/ -lmetis  -I src/Base/ -I src/Kernels/   -I src/Lib/ -DMF_TIMING -DMF_MPICH -DMF_OPENMP -DOMP_CellColor -DOMP_DIVREP

    mpicxx -std=c++11 -O3 -fopenmp -o output/main_mpi_task src/main.cpp src/Base/*.cpp src/Kernels/*.cpp src/Lib/libTDTree.so -L src/Lib/ -lmetis  -I src/Base/ -I src/Kernels/   -I src/Lib/ -DMF_TIMING -DMF_MPICH -DTDTREE

    mpicxx -std=c++11 -O3 -fopenmp -o output/main_gaspi_task src/main.cpp src/Base/*.cpp src/Kernels/*.cpp  src/Lib/libTDTree.so -L src/Lib/ -lmetis /public5/home/sch8427/GPI-2-next/install/lib64/libGPI2.so -I/public5/home/sch8427/GPI-2-next/install/include/  -I src/Base/ -I src/Kernels/  -I src/Lib/   -DMF_TIMING -DMF_MPICH -DGASPI -DTDTREE

    mpicxx -std=c++11 -O3 -fopenmp -o output/main_mpi_task src/main.cpp src/Base/*.cpp src/Kernels/*.cpp src/src/*.cpp -L src/Lib/ -lmetis  -I src/Base/ -I src/Kernels/ -I src/headers/   -I src/Lib/ -DMF_TIMING -DMF_MPICH -DTDTREE -DOMP

    mpicxx -std=c++11 -O3 -fopenmp -o output/main_gaspi_task src/main.cpp src/Base/*.cpp src/Kernels/*.cpp  src/src/*.cpp -L src/Lib/ -lmetis /public5/home/sch8427/GPI-2-next/install/lib64/libGPI2.so -I/public5/home/sch8427/GPI-2-next/install/include/  -I src/Base/ -I src/Kernels/  -I src/Lib/ -I src/headers/  -DMF_TIMING -DMF_MPICH -DGASPI -DTDTREE -DOMP
## Others

    export OMP_NUM_THREADS=32
    numactl --cpunodebind=0 --membind=0 ./output/main_facecoloring > result/370w20step/facecoloring1.out
    numactl --cpunodebind=0 --membind=0 ./output/main_reduction > result/370w20step/reduction1.out
    numactl --cpunodebind=0 --membind=0 ./output/main_divrep > result/370w20step/divrep1.out
    numactl --cpunodebind=0 --membind=0 ./output/main_task > result/370w20step/task1.out
    numactl --cpunodebind=0 --membind=0 ./output/main_forkjoin > result/370w20step/forkjoin1.out

    chmod a+x run.sh


g++ -std=c++11 -O3 -o output/main_serial src/main.cpp src/Base/*.cpp src/Kernels/*.cpp  -I src/Base/ -I src/Kernels/   -DMF_TIMING

g++ -std=c++11 -O3 -fopenmp -o output/main_facecoloring src/main.cpp src/Base/*.cpp src/Kernels/*.cpp  -I src/Base/ -I src/Kernels/   -DMF_TIMING -DMF_OPENMP -DOMP_CellColor -DOMP_FaceColor

g++ -std=c++11 -O3 -fopenmp -o output/main_task src/main.cpp src/Base/*.cpp src/Kernels/*.cpp src/Lib/libTDTree.so -L src/Lib/ -lmetis  -I src/Base/ -I src/Kernels/   -I src/Lib/ -DMF_TIMING -DTDTREE

g++ -std=c++11 -O3 -fopenmp -o output/main_reduction src/main.cpp src/Base/*.cpp src/Kernels/*.cpp  -I src/Base/ -I src/Kernels/   -DMF_TIMING -DMF_OPENMP -DOMP_CellColor -DOMP_Reduction

g++ -std=c++11 -O3 -fopenmp -o output/main_divrep src/main.cpp src/Base/*.cpp src/Kernels/*.cpp -L src/Lib/ -lmetis  -I src/Base/ -I src/Kernels/   -I src/Lib/ -DMF_TIMING -DMF_OPENMP -DOMP_CellColor -DOMP_DIVREP

--------------------------------------------------

g++ -std=c++11 -O3 -o output/main_ori_serial src/main.cpp src/Base/*.cpp src/Kernels/*.cpp  -I src/Base/ -I src/Kernels/   -DMF_TIMING

g++ -std=c++11 -O3 -fopenmp -o output/main_ori_facecoloring src/main.cpp src/Base/*.cpp src/Kernels/*.cpp  -I src/Base/ -I src/Kernels/   -DMF_TIMING -DMF_OPENMP -DOMP_CellColor -DOMP_FaceColor

g++ -std=c++11 -O3 -fopenmp -o output/main_ori_task src/main.cpp src/Base/*.cpp src/Kernels/*.cpp src/Lib/libTDTree.so -L src/Lib/ -lmetis  -I src/Base/ -I src/Kernels/   -I src/Lib/ -DMF_TIMING -DTDTREE

g++ -std=c++11 -O3 -fopenmp -o output/main_ori_reduction src/main.cpp src/Base/*.cpp src/Kernels/*.cpp  -I src/Base/ -I src/Kernels/   -DMF_TIMING -DMF_OPENMP -DOMP_CellColor -DOMP_Reduction

g++ -std=c++11 -O3 -fopenmp -o output/main_ori_divrep src/main.cpp src/Base/*.cpp src/Kernels/*.cpp -L src/Lib/ -lmetis  -I src/Base/ -I src/Kernels/   -I src/Lib/ -DMF_TIMING -DMF_OPENMP -DOMP_CellColor -DOMP_DIVREP


papi:
---------------------------------------------------------

g++ -std=c++11 -O3 -o output/main_serial src/main.cpp src/Base/*.cpp src/Kernels/*.cpp  src/Monitoring/papi_func.cpp -L /public5/home/sch8427/para/papi-install/lib/ -lpapi -I /public5/home/sch8427/para/papi-install/include/ -I src/Base/ -I src/Kernels/   -DMF_TIMING -DMF_PAPI

g++ -std=c++11 -O3 -fopenmp -o output/main_facecoloring src/main.cpp src/Base/*.cpp src/Kernels/*.cpp  src/Monitoring/papi_func.cpp -L /public5/home/sch8427/para/papi-install/lib/ -lpapi -I /public5/home/sch8427/para/papi-install/include/ -I src/Base/ -I src/Kernels/   -DMF_TIMING -DMF_OPENMP -DOMP_CellColor -DOMP_FaceColor -DMF_PAPI

g++ -std=c++11 -O3 -fopenmp -o output/main_task src/main.cpp src/Base/*.cpp src/Kernels/*.cpp src/Lib/libTDTree.so -L src/Lib/ -lmetis  src/Monitoring/papi_func.cpp -L /public5/home/sch8427/para/papi-install/lib/ -lpapi -I /public5/home/sch8427/para/papi-install/include/ -I src/Base/ -I src/Kernels/   -I src/Lib/ -DMF_TIMING -DTDTREE -DMF_PAPI

g++ -std=c++11 -O3 -fopenmp -o output/main_reduction src/main.cpp src/Base/*.cpp src/Kernels/*.cpp  src/Monitoring/papi_func.cpp -L /public5/home/sch8427/para/papi-install/lib/ -lpapi -I /public5/home/sch8427/para/papi-install/include/ -I src/Base/ -I src/Kernels/   -DMF_TIMING -DMF_OPENMP -DOMP_CellColor -DOMP_Reduction -DMF_PAPI

g++ -std=c++11 -O3 -fopenmp -o output/main_divrep src/main.cpp src/Base/*.cpp src/Kernels/*.cpp -L src/Lib/ -lmetis  src/Monitoring/papi_func.cpp -L /public5/home/sch8427/para/papi-install/lib/ -lpapi -I /public5/home/sch8427/para/papi-install/include/ -I src/Base/ -I src/Kernels/   -I src/Lib/ -DMF_TIMING -DMF_OPENMP -DOMP_CellColor -DOMP_DIVREP -DMF_PAPI

g++ -std=c++11 -O3 -o output/main_ori_serial src/main.cpp src/Base/*.cpp src/Kernels/*.cpp  src/Monitoring/papi_func.cpp -L /public5/home/sch8427/para/papi-install/lib/ -lpapi -I /public5/home/sch8427/para/papi-install/include/ -I src/Base/ -I src/Kernels/   -DMF_TIMING -DMF_PAPI

g++ -std=c++11 -O3 -fopenmp -o output/main_ori_facecoloring src/main.cpp src/Base/*.cpp src/Kernels/*.cpp  src/Monitoring/papi_func.cpp -L /public5/home/sch8427/para/papi-install/lib/ -lpapi -I /public5/home/sch8427/para/papi-install/include/ -I src/Base/ -I src/Kernels/   -DMF_TIMING -DMF_OPENMP -DOMP_CellColor -DOMP_FaceColor -DMF_PAPI

g++ -std=c++11 -O3 -fopenmp -o output/main_ori_task src/main.cpp src/Base/*.cpp src/Kernels/*.cpp src/Lib/libTDTree.so -L src/Lib/ -lmetis  src/Monitoring/papi_func.cpp -L /public5/home/sch8427/para/papi-install/lib/ -lpapi -I /public5/home/sch8427/para/papi-install/include/ -I src/Base/ -I src/Kernels/   -I src/Lib/ -DMF_TIMING -DTDTREE -DMF_PAPI

g++ -std=c++11 -O3 -fopenmp -o output/main_ori_reduction src/main.cpp src/Base/*.cpp src/Kernels/*.cpp  src/Monitoring/papi_func.cpp -L /public5/home/sch8427/para/papi-install/lib/ -lpapi -I /public5/home/sch8427/para/papi-install/include/ -I src/Base/ -I src/Kernels/   -DMF_TIMING -DMF_OPENMP -DOMP_CellColor -DOMP_Reduction -DMF_PAPI

g++ -std=c++11 -O3 -fopenmp -o output/main_ori_divrep src/main.cpp src/Base/*.cpp src/Kernels/*.cpp -L src/Lib/ -lmetis  src/Monitoring/papi_func.cpp -L /public5/home/sch8427/para/papi-install/lib/ -lpapi -I /public5/home/sch8427/para/papi-install/include/ -I src/Base/ -I src/Kernels/   -I src/Lib/ -DMF_TIMING -DMF_OPENMP -DOMP_CellColor -DOMP_DIVREP -DMF_PAPI