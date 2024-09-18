/*!
 * @file        parallel_mpi.cpp
 * @brief
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @version     1.0
 * @date        2023-06-23
 * @copyright   Copyright (c) 2023  PERSONAL
 *
 * @par 修改日志:
 * <table>
 * <tr><th> Date        <th> Version  <th> Author  <th> Description
 * <tr><td> 2023-06-23  <td> 1.0      <td> Wisces  <td> 添加了 MPI 并行（使用条件编译）
 * -# 添加了 SetUpComm()、SetUpComm_Node() 函数，用于建立边界面 / 边界点的传值机制
 * -# 添加了 CommInterfaceDataMPI()、RecvSendVarNeighbor() 函数，用于并行传值
 * </table>
 */
// direct head file
#include "grid_polyhedra.h"

//!< C/C++ head files
#include <cstdio>  ///printf
#include <cstdlib>
#include <cmath>   ///< ceil()
// #include <cassert>

//!< user difined head file
#include "memory_util.h"
#include "para_field_global.h"

//!< head file relying on condition-compiling
#ifdef MF_MPICH
#include <mpi.h>
#endif

#ifdef MF_MPICH
extern int myid;
extern int numprocs;
extern int myZone;        // = myid+1
extern MPI_Comm GridComm; ///< for each grid
#endif

#ifdef MF_MPICH
IntType *nbn;//对应的邻块中本地块的编号
/*!
 * @brief       建立边界面的传值机制
 * @note        modify according to the fun [void PolyGrid::SetUpComm()]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-23
 */
void PolyGrid::SetUpComm()
{
    IntType n, i, ni, g, np;
    BCRecord **bcr = Getbcr();
    IntType nbZone;
    MPI_Status status;
    MPI_Request *req_send, *req_recv;
    MPI_Status *status_array;

    IntType nNeighbor = this->GetNumberOfFaceNeighbors();

    if (nNeighbor == 0)
        return;
    status_array = NULL;
    req_send = NULL;
    req_recv = NULL;
    nZIFace = NULL; // No. of interfaces for each neighbor
    snew_array_1D(status_array, nNeighbor);
    snew_array_1D(req_send, nNeighbor);
    snew_array_1D(req_recv, nNeighbor);
    snew_array_1D(nZIFace, nNeighbor);

#ifdef GASPI
    // IntType maxnNeighbor;
    // MPI_Allreduce(&nNeighbor, &maxnNeighbor, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    cout << "init begin" << endl << flush;
    gaspi_config_t config;
    SUCCESS_OR_DIE (gaspi_config_get (&config));
    config.build_infrastructure = GASPI_TOPOLOGY_DYNAMIC;
    SUCCESS_OR_DIE (gaspi_config_set (config));
    SUCCESS_OR_DIE (gaspi_proc_init (GASPI_BLOCK));

    SUCCESS_OR_DIE (gaspi_proc_rank (&GASPI_iProc));
    SUCCESS_OR_DIE (gaspi_proc_num (&GASPI_nProc));

    unsigned int ss=0;
    gaspi_segment_max(&ss);
    cout << "segment_max ="<< ss << endl << flush;

    nC2B = new IntType[nTCell]();
    C2B = new pair<int,int>*[nTCell]();
    vector<pair<int,int>> *c2b = new vector<pair<int,int>>[nTCell]();
    //f2b = new vector<pair<int,int>>[nTFace]();
    // GASPI_sendDataSegmentId = new gaspi_segment_id_t[GASPI_nProc];
    // GASPI_recvDataSegmentId = new gaspi_segment_id_t[GASPI_nProc];
    // GASPI_sendOffsetSegmentId = new gaspi_segment_id_t[GASPI_nProc];
    // GASPI_recvOffsetSegmentId = new gaspi_segment_id_t[GASPI_nProc];
    // GASPI_sendDataPtr = new gaspi_pointer_t[GASPI_nProc];
    // GASPI_recvDataPtr = new gaspi_pointer_t[GASPI_nProc];
    // GASPI_sendOffsetPtr = new gaspi_pointer_t[GASPI_nProc];
    // GASPI_recvOffsetPtr = new gaspi_pointer_t[GASPI_nProc];
    // GASPI_sendData = new double*[GASPI_nProc];
    // GASPI_recvData = new double*[GASPI_nProc];
    // GASPI_sendOffset = new int*[GASPI_nProc];
    // GASPI_recvOffset = new int*[GASPI_nProc];
    // GASPI_lock = new omp_lock_t[GASPI_nProc];
    // GASPI_offset = new int[GASPI_nProc]();
    // GASPI_commed = new int[GASPI_nProc]();
    GASPI_mapping = new int[GASPI_nProc]();
    GASPI_mappingN = new int[GASPI_nProc]();
    GASPI_recvNotifications = new int[GASPI_nProc]();
    GASPI_need_Comm = new int[GASPI_nProc]();

    for(int i=0;i<GASPI_nProc;i++)
    {
        GASPI_need_Comm[i]=0;
    }

#ifdef TDTREE
    GASPI_blockSize = 4096;
    GASPI_lock = new omp_lock_t[GASPI_nProc];
    GASPI_offset = new int[GASPI_nProc]();
    GASPI_commed = new int[GASPI_nProc]();
    for(g=0; g<GASPI_nProc; g++)
    {
        omp_init_lock(&GASPI_lock[g]);
        GASPI_offset[g] = GASPI_commed[g] = 0;
    }    
#endif
    // for(g=0; g<GASPI_nProc; g++)
    // {
    //     // omp_init_lock(&GASPI_lock[g]);
    //     GASPI_sendDataSegmentId[g] = GASPI_recvDataSegmentId[g] = g;
    //     GASPI_sendOffsetSegmentId[g] = GASPI_recvOffsetSegmentId[g] = g + GASPI_nProc;
    // }
    GASPI_sendDataSegmentId = 0;
    GASPI_recvDataSegmentId = 1;
    GASPI_sendOffsetSegmentId = 2;
    GASPI_recvOffsetSegmentId = 3;
    // gaspi_number_t num;
    // gaspi_queue_num(&num);
    // cout << num << endl << flush;

    maxLength = 0;

    snew_array_1D(Ghost, nNeighbor);
#endif
    // bCNo = NULL; // cell number in the current zone(process) for each neighbor, each interface, ordered according to the current zone(process). Used to send data on cell center to other processor.
    snew_array_1D(bCNo, nNeighbor);
    IntType **bqs = NULL, **bqr = NULL;
    snew_array_1D(bqr, nNeighbor);
    snew_array_1D(bqs, nNeighbor);
    // ofstream fp ;
    // fp.open( "data" , ofstream::app|ofstream::out);
    for (g = 0; g < nNeighbor; g++)
    {
        nbZone = nb[g];

        // 计数与其相邻的各块分别有多少个并行边界面
        n = 0;
        ni = 0;
        for (i = 0; i < nBFace; i++)
        {
            if (bcr[i]->GetType() == INTERFACE)
            {
                if (nbZ[ni] == nbZone)
                {
                    n++;
                }
                ni++;
            }
        }
        nZIFace[g] = n;
        // printf("%d %d %d \n", myZone-1,nbZone,n);
        // cout << myZone-1 << "  "<< nbZone << "  " << n <<endl;
        // fp << myZone-1 << "  "<< nbZone << "  " << n <<endl;
#ifdef GASPI
        GASPI_mapping[nbZone] = g;
        if (n > maxLength)
            maxLength = n;
        Ghost[g] = NULL; 
        snew_array_1D(Ghost[g], n);
#endif         
        // 将当前块中的需要传递的体积号储存在bCNo[g][n]中:g代表将要传给第几块,n表示序列号
        bCNo[g] = NULL;
        bqs[g] = NULL;
        snew_array_1D(bCNo[g], n);
        snew_array_1D(bqs[g], n);
        n = 0;
        ni = 0;
        for (i = 0; i < nBFace; i++)
        {
            if (bcr[i]->GetType() == INTERFACE)
            {
                if (nbZ[ni] == nbZone)
                {
                    bCNo[g][n] = f2c[i * 2];
#ifdef GASPI
                    c2b[f2c[i*2]].push_back({g,n});
                    // c2b[f2c[i*2]].push_back({g,nbBF[ni]});
                    //f2b[i].push_back({g,n});
                    Ghost[g][n] = nbBF[ni];
#endif
                    bqs[g][n++] = nbBF[ni]; // parallel information for interfaces: face No. in that zone(process). Size is nIFace.
                }
                ni++;
            }
        }
    }
    // fp.close();
#ifdef GASPI
#ifdef TDTREE
    omp_set_max_active_levels(2);
    int blocks = ceil(1.0 * maxLength / GASPI_blockSize);
    MPI_Allreduce (&blocks, &GASPI_nbBlock, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    maxLength = GASPI_nbBlock * GASPI_blockSize * 4;  
#else  
    int blocks;
    MPI_Allreduce(&maxLength, &blocks, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    maxLength = blocks  * 4;
#endif
    cout << "maxLength = " << maxLength << endl << flush;  

    SUCCESS_OR_DIE(gaspi_segment_create
    (GASPI_sendDataSegmentId,
    20*maxLength*GASPI_nProc*sizeof(double),
    GASPI_GROUP_ALL,
    GASPI_BLOCK,
    GASPI_ALLOC_DEFAULT
    ));

    SUCCESS_OR_DIE(gaspi_segment_create
    (GASPI_recvDataSegmentId,
    20*maxLength*GASPI_nProc*sizeof(double),
    GASPI_GROUP_ALL,
    GASPI_BLOCK,
    GASPI_ALLOC_DEFAULT
    ));

    SUCCESS_OR_DIE(gaspi_segment_ptr (GASPI_sendDataSegmentId,
                &GASPI_sendDataPtr));

    SUCCESS_OR_DIE(gaspi_segment_ptr (GASPI_recvDataSegmentId,
                &GASPI_recvDataPtr));

    GASPI_sendData = (double *)GASPI_sendDataPtr;
    GASPI_recvData = (double *)GASPI_recvDataPtr;
    
    for (int i = 0 ; i < 20*maxLength*GASPI_nProc; i++)
    {
        GASPI_sendData[i] = 0.;
        GASPI_recvData[i] = 0.;
    }

    for(int i=0;i<nTCell;i++)
    {
        nC2B[i] = c2b[i].size();
        snew_array_1D(C2B[i], (int)c2b[i].size());
        for(int j=0;j<nC2B[i];j++)
            C2B[i][j] = c2b[i][j];
    }
    delete[] c2b;
#ifdef TDTREE
    SUCCESS_OR_DIE(gaspi_segment_create
        (GASPI_sendOffsetSegmentId,
        maxLength*GASPI_nProc*sizeof(int),
        GASPI_GROUP_ALL,
        GASPI_BLOCK,
        GASPI_ALLOC_DEFAULT
        ));
    SUCCESS_OR_DIE(gaspi_segment_create
        (GASPI_recvOffsetSegmentId,
        maxLength*GASPI_nProc*sizeof(int),
        GASPI_GROUP_ALL,
        GASPI_BLOCK,
        GASPI_ALLOC_DEFAULT
        ));

    SUCCESS_OR_DIE(gaspi_segment_ptr (GASPI_sendOffsetSegmentId,
                &GASPI_sendOffsetPtr));

    SUCCESS_OR_DIE(gaspi_segment_ptr (GASPI_recvOffsetSegmentId,
                &GASPI_recvOffsetPtr));

    GASPI_sendOffset = (int *)GASPI_sendOffsetPtr;
    GASPI_recvOffset = (int *)GASPI_recvOffsetPtr;

    for (int i = 0 ; i < maxLength*GASPI_nProc; i++)
    {
        GASPI_sendOffset[i] = 0;
        GASPI_recvOffset[i] = 0;
    }

    GASPI_notifications = 0;
    for(g=0; g<nNeighbor; g++)
        GASPI_notifications += nZIFace[g] / GASPI_blockSize + (nZIFace[g] % GASPI_blockSize > 0); 

#endif
    // MPI_Request *req_sendgaspi, *req_recvgaspi;
    // MPI_Status *status_arraygaspi;

    // status_arraygaspi = NULL;
    // req_sendgaspi = NULL;
    // req_recvgaspi = NULL;
    // snew_array_1D(status_arraygaspi, nNeighbor);
    // snew_array_1D(req_sendgaspi, nNeighbor);
    // snew_array_1D(req_recvgaspi, nNeighbor);
    // snew_array_1D(nbn, nNeighbor);
    // for (g = 0 ; g < nNeighbor ; g++ )
    // {
    //     nbZone = nb[g];
    //     MPI_Isend(&g , 1 ,MPIIntType, nbZone, level, MPI_COMM_WORLD, &req_sendgaspi[g]);
    //     MPI_Irecv(&nbn[g], 1, MPIIntType, nbZone, level, MPI_COMM_WORLD, &req_recvgaspi[g]);
    // }
    // MPI_Waitall(nNeighbor, req_sendgaspi, status_arraygaspi);
    // MPI_Waitall(nNeighbor, req_recvgaspi, status_arraygaspi);
    // sdel_array_1D(req_sendgaspi);
    // sdel_array_1D(req_recvgaspi);
    // sdel_array_1D(status_arraygaspi);
#endif
    // now receiving
    // for (g = 0; g < nNeighbor; g++)
    // {
    //     nbZone = nb[g];
    //     MPI_Send(&nZIFace[g], 1, MPIIntType, nbZone, level, MPI_COMM_WORLD);
    // }
    // for (g = 0; g < nNeighbor; g++)
    // {
    //     nbZone = nb[g];
    //     MPI_Recv(&n, 1, MPIIntType, nbZone, level, MPI_COMM_WORLD, &status);
    //     assert(n == nZIFace[g]);
    // }

    for (g = 0; g < nNeighbor; g++)
    {
        bqr[g] = NULL;
        snew_array_1D(bqr[g], nZIFace[g]);
    }

    for (np = 1; np <= numprocs; np++)
    // for (np = 0; np < numprocs; np++)
    {
        if (myZone == np)
        {
            //       send   to other
            for (g = 0; g < nNeighbor; g++)
            {
                nbZone = nb[g];
                MPI_Isend(bqs[g], nZIFace[g], MPIIntType, nbZone, level, MPI_COMM_WORLD, &req_send[g]);
            }
        }
        else
        {
            //       receive   from np
            for (g = 0; g < nNeighbor; g++)
            {
                nbZone = nb[g];
                if (nbZone == (np - 1))
                {
                    MPI_Irecv(bqr[g], nZIFace[g], MPIIntType, nbZone, level, MPI_COMM_WORLD, &req_recv[g]);
                }
            }
        }
    } // np

    MPI_Waitall(nNeighbor, req_recv, status_array);
    // bFNo = NULL; // boundary face number in the current zone(process) for each neighbor, each interface, ordered according to the neighbor zone(process).  Used to receive data on cell center from other processor.
    snew_array_2D(bFNo, nNeighbor, nZIFace, false);
    for (g = 0; g < nNeighbor; g++)
    {
        for (i = 0; i < nZIFace[g]; i++)
        {
            bFNo[g][i] = bqr[g][i];
        }
    }
    MPI_Waitall(nNeighbor, req_send, status_array);

    for (g = 0; g < nNeighbor; g++)
    {
        sdel_array_1D(bqr[g]);
        sdel_array_1D(bqs[g]);
    }
    sdel_array_1D(bqr);
    sdel_array_1D(bqs);
    sdel_array_1D(req_send);
    sdel_array_1D(req_recv);
    sdel_array_1D(status_array);

    // 检测并行分区是否存在非一一对应的现象
    IntType *count = NULL;
    snew_array_1D(count, nBFace);
    IntType mark = 1;
    for (i = 0; i < nBFace; i++)
        count[i] = 0;
    for (g = 0; g < nNeighbor; g++)
    {
        for (i = 0; i < nZIFace[g]; i++)
        {
            count[bFNo[g][i]]++;
        }
    }

    for (i = 0; i < nBFace; i++)
    {
        if (count[i] > 1)
        {
            cout << "myZone = " << myZone << ", and the face number is" << i << endl;
            cout << "The count is " << count[i] << endl;
            mark = 0;
        }
    }
    sdel_array_1D(count);
    if (!mark)
    {
        cerr << "Error in func[DivRep()] of file[parallel_omp.cpp]!" << endl;
        exit(1);
    }
}

#ifdef GASPI

#ifdef TDTREE

void PolyGrid::GASPIPrintf()
{
    String senddata_doc , recvdata_doc;
    // File *fp1 , *fp2;
    IntType nbZone;
    IntType nVar = 5;
    RealFlow data0;
    int i , j , k , g;

    sprintf(senddata_doc,"SendData%i.txt",GASPI_iProc);
    sprintf(recvdata_doc,"RecvData%i.txt",GASPI_iProc); 

    ofstream fp1 ,fp2 ;
    fp1.open( senddata_doc , ofstream::app|ofstream::out);

    // fp1=fopen(senddata_doc,"a");

    for (g = 0 ; g < nNeighbor ; g++)
    {
        nbZone = nb[g];
        for (i = 0 ; i < nZIFace[g] ; i++)
        {
            fp1 << nbZone << " " << i << " ";
            for (j = 0; j < nVar; j++)
            {
                fp1 << GASPI_sendData[i*nVar+j+maxLength*20*nbZone] << " ";
            }
            fp1 << endl ;
        }
    }
    // for (i = 0 ;i < GASPI_nProc ; i++)
    // {
    //     for (j = 0 ;j < 4 * maxLength ; j++)
    //     {
    //         if (GASPI_sendData[j*nVar+maxLength*20*i] != 100)
    //         {
    //             fp1 << i << " " << j <<" ";
    //             for (k = 0; k < nVar ; k++)
    //             {
    //                 fp1 << GASPI_sendData[j*nVar+k+maxLength*20*i] << " ";
    //             }
    //             fp1 << endl ;
    //         }
    //     }
    // }
    fp1 << endl ;
    fp1.close();

    fp2.open( recvdata_doc , ofstream::app|ofstream::out);

    for (g = 0 ; g < nNeighbor ; g++)
    {
        nbZone = nb[g];
        for (i = 0 ; i < nZIFace[g] ; i++)
        {
            fp2 << nbZone << " "  << i << " ";
            for (j = 0; j < nVar; j++)
            {
                fp2 << GASPI_recvData[i*nVar+j+maxLength*20*nbZone] << " ";
            }
            fp2 << endl ;
        }
    }
    // for (i = 0 ;i < GASPI_nProc ; i++)
    // {
    //     for (j = 0 ;j < 4 * maxLength ; j++)
    //     {
    //         if (GASPI_recvData[j*nVar+maxLength*20*i] != 100)
    //         {
    //             fp2 << i << " " << j <<" ";
    //             for (k = 0; k < nVar ; k++)
    //             {
    //                 fp2 << GASPI_recvData[j*nVar+k+maxLength*20*i] << " ";
    //             }
    //             fp2 << endl ;
    //         }
    //     }
    // }
    fp2 << endl ;
    fp2.close();    


}

void initBoundryEntity(char **userArgs, TDTreeArgs *treeArgs)
{	
	IntType nsCell = treeArgs->firstCell;
	IntType neCell = treeArgs->lastCell + 1;
	IntType nsNode = treeArgs->firstNode;
	IntType neNode = treeArgs->lastNode + 1;

    PolyGrid *grid = (PolyGrid *)userArgs;    
    int *nC2B = grid->GetnC2B();
    int *nN2B = grid->GetnN2B();

    (*treeArgs->boundryNbCell) = 0;
    (*treeArgs->boundryNbNode) = 0;

    IntType i, n=0;
    for(i=nsCell;i<neCell;i++)
        if (nC2B[i]>0)
            (*treeArgs->boundryNbCell)++;
    (*treeArgs->boundryCell) = new int[(*treeArgs->boundryNbCell)];
    for(i=nsCell;i<neCell;i++)
        if (nC2B[i]>0)
            (*treeArgs->boundryCell)[n++] = i;

    //////////////////add by zz/////////////////////////
    IntType **GHOST = grid->GetGhost();
    pair<int,int> **C2B = grid->GetC2B();
    IntType nNeighbor = grid->GetNumberOfFaceNeighbors();
    IntType nNo, No;
    (*treeArgs->boundryN) = new int[nNeighbor];
    for (i = 0 ; i < nNeighbor; i++){
        (*treeArgs->boundryN)[i] = 0;
    }
    IntType ilu, cell, ghost;
    for(ilu=0;ilu<(*treeArgs->boundryNbCell);ilu++){
        cell = (*treeArgs->boundryCell)[ilu];
        for(i=0;i<nC2B[cell];i++)
        {
            int neighbor = C2B[cell][i].first;
            No = C2B[cell][i].second;
            (*treeArgs->boundryN)[neighbor]++;
        }
    }
    int num=0;
    for (i = 0 ; i < nNeighbor; i++){
        num+=(*treeArgs->boundryN)[i];
    }
    (*treeArgs->boundryCellNo) = new int[num];
    IntType ni[nNeighbor];//计数
    for (i = 0 ; i < nNeighbor; i++){
        ni[i] = 0;
    } 
    for(ilu=0;ilu<(*treeArgs->boundryNbCell);ilu++){
        cell = (*treeArgs->boundryCell)[ilu];
        for(i=0;i<nC2B[cell];i++)
        {
            int neighbor = C2B[cell][i].first;
            No = C2B[cell][i].second;
            num = 0;
            for (int j = 0; j < neighbor; j++){
                num += (*treeArgs->boundryN)[j];
            }
            (*treeArgs->boundryCellNo)[num + ni[neighbor]] = No;
            // cout<< "r "<< GASPI_iProc <<" " <<num + ni[neighbor] << " "<< No <<endl;
            ni[neighbor]++;
        }
    }   
    ////////////////////////////////////////////////// 

    n=0;
    for(i=nsNode;i<neNode;i++)
        if (nN2B[i]>0)
            (*treeArgs->boundryNbNode)++;
    (*treeArgs->boundryNode) = new int[(*treeArgs->boundryNbNode)];
    for(i=nsNode;i<neNode;i++)
        if (nN2B[i]>0)
            (*treeArgs->boundryNode)[n++] = i;
}

void PolyGrid::GASPI_wait_for_queue_half_full (gaspi_queue_id_t queueID)
{
    gaspi_number_t queueSizeMax;
    gaspi_number_t queueSize;

    SUCCESS_OR_DIE (gaspi_queue_size_max (&queueSizeMax));
    SUCCESS_OR_DIE (gaspi_queue_size (queueID, &queueSize));

    // cout << queueSize << " " << queueSizeMax << endl << flush;
    // if (queueSize >= queueSizeMax)
    //     cout << "ERROR!!!!!!!!!!!!!!!" << endl << flush;

    if (queueSize >= queueSizeMax / 2) {
        SUCCESS_OR_DIE (gaspi_wait (queueID, GASPI_BLOCK));
    }
}

void PolyGrid::GASPI_wait_for_flush_queues ()
{
    gaspi_number_t queue_num;
    gaspi_queue_num(&queue_num);

    gaspi_queue_id_t queue = 0;
    while(queue < queue_num)
    {
        SUCCESS_OR_DIE( gaspi_wait(queue, GASPI_BLOCK) );
        queue++;
    }
}
#endif

/*
 多种通讯实现
*/

#ifdef TDTREE
// void PolyGrid::GASPIWaitAll(IntType nVar, RealFlow **q)
// {
//     IntType nNeighbor = this->GetNumberOfFaceNeighbors();

//     gaspi_barrier (GASPI_GROUP_ALL, GASPI_BLOCK);
//     // GASPI_wait_for_flush_queues ();

//     for (int i = 0; i < nNeighbor; i++)
//     {
//         int nbZone = nb[i];
// // #pragma omp parallel for
//         for(int j=GASPI_recvNotifications[nbZone]*GASPI_blockSize, offset =  j*nVar; j<nZIFace[i]; j++)
//         {
//             int No = GASPI_recvOffset[nbZone*maxLength+j];
//             int ghost    = nTCell + bFNo[i][No];
//             for(int k=0; k<nVar; k++)
//             {
//                 q[k][ghost] = GASPI_recvData[nbZone*maxLength*20+offset];
//                 offset++;
//             }
//         }
//         GASPI_offset[nbZone] = GASPI_commed[nbZone] = GASPI_recvNotifications[nbZone] = 0;
//     }
// }

void PolyGrid::GASPIWaitAll(IntType nVar, RealFlow **q)
{
    IntType nNeighbor = this->GetNumberOfFaceNeighbors();
    SUCCESS_OR_DIE(gaspi_wait(0,GASPI_BLOCK));
    SUCCESS_OR_DIE (gaspi_barrier (GASPI_GROUP_ALL, GASPI_BLOCK)); 
    // GASPI_wait_for_flush_queues ();
// #pragma omp parallel for num_threads(2)
#pragma omp parallel for
    for (int i = 0; i < nNeighbor; i++)
    {
        int nbZone = nb[i];
// #pragma omp parallel for
        for(int j=0, offset =  j*nVar; j<nZIFace[i]; j++)
        // for(int j=0; j<nZIFace[i]; j++)
        {
            // int offset =  j*nVar;
            int No = GASPI_recvOffset[nbZone*maxLength+j];
            int ghost    = nTCell + No;
            // int ghost    = nTCell + bFNo[i][No];
            for(int k=0; k<nVar; k++)
            {
                q[k][ghost] = GASPI_recvData[nbZone*maxLength*20+offset];
                offset++;
            }
        }
        GASPI_offset[nbZone] = GASPI_commed[nbZone] = 0;
    }
}
// void PolyGrid::GASPIUnzip(IntType nVar, RealFlow **q, IntType neighbor, IntType ns, IntType ne)
// {
//     ne += ns;
//     IntType j, k, ghost;
//     int nbZone = nb[neighbor];
//     for(j=ns; j<ne; j++) {
//         {
//             int No = GASPI_recvOffset[nbZone*maxLength+j];
//             ghost    = nTCell + bFNo[neighbor][No];
//             for(k=0; k<nVar; k++) {
//                 q[k][ghost] = GASPI_recvData[nbZone*maxLength*20+j*nVar+k];
//             }
//         }
//     }
// }

// void PolyGrid::GASPIWaitAll_Node(IntType nVar, RealFlow **q)
// {
// // #pragma omp parallel
// // #pragma omp single
// {    
// // #pragma omp parallel for
//     for(int g = 0; g < nNeighborN; g++)
//     {
//         gaspi_notification_id_t notifyID = 0;
//         gaspi_notification_t notifyValue = 0;

//         // #pragma omp task default(shared)
//         // #pragma omp task
//         {
//             int nbZone = nbN[g];
//             for(int i = 0; i < GASPI_notifications_node[g]; i++)
//             {
//                 gaspi_notification_id_t notifyID = 0;
//                 gaspi_notification_t notifyValue = 0;

//                 SUCCESS_OR_DIE (gaspi_notify_waitsome (GASPI_recvOffsetSegmentId,
//                     0,
//                     GASPI_nbBlock * GASPI_nProc,
//                     &notifyID,
//                     GASPI_BLOCK));

//                 SUCCESS_OR_DIE (gaspi_notify_reset (GASPI_recvOffsetSegmentId,
//                     notifyID,
//                     &notifyValue));

//                 notifyID /= GASPI_nProc;
//                 // cout << "recv " << nbZone << " " << GASPI_iProc << " " << notifyID * GASPI_blockSize << endl << flush;
                        
//                 // #pragma omp task default(shared)
//                 // #pragma omp task
//                 GASPIUnzip_Node(nVar, q, g, notifyID*GASPI_blockSize, notifyValue);
//             }
//             GASPI_offset[nbZone] = GASPI_commed[nbZone] = 0;
//         }

//     }
    
//     // #pragma omp taskwait
//     // SUCCESS_OR_DIE( gaspi_wait(0, GASPI_BLOCK) );
// }
//     // gaspi_barrier (GASPI_GROUP_ALL, GASPI_BLOCK);
//     GASPI_wait_for_flush_queues ();
// }

// void PolyGrid::GASPIUnzip_Node(IntType nVar, RealFlow **q, IntType neighbor, IntType ns, IntType ne)
// {
//     ne += ns;
//     IntType j, k, ghost;
//     int nbZone = nbN[neighbor];
//     // #pragma omp parallel for
//     for(j=ns; j<ne; j++) {
//         // #pragma omp task default(shared)
//         // #pragma omp task
//         {        
//             int No = GASPI_recvOffset[nbZone*maxLength+j];
//             int node = bNRNo[neighbor][No];
//             for(k=0; k<nVar; k++) {
//                 q[k][node] += GASPI_recvData[nbZone*maxLength+j*nVar+k];
//             }
//         }
//     }
// }

// void PolyGrid::GASPIComm(RealFlow **q, int *boundryCell, int boundryNbCell)
// {
//     pair<int,int> **C2B = this->GetC2B();
//     IntType *nC2B = this->GetnC2B();
//     IntType **GHOST = this->GetGhost();
//     IntType nNo, No, i;
//     IntType nVar = 5;
 
//     IntType ilu, cell, ghost;
//     for(ilu=0;ilu<boundryNbCell;ilu++){
//         cell = boundryCell[ilu];
//         for(i=0;i<nC2B[cell];i++)
//         {
//             int neighbor = C2B[cell][i].first;
//             No = C2B[cell][i].second;
//             ghost = GHOST[neighbor][No];

//             bool needComm = false;
//             int nbZone = nb[neighbor];

//             omp_set_lock(&GASPI_lock[nbZone]);

//             int offset = GASPI_offset[nbZone];
//             int commed = GASPI_commed[nbZone];
//             GASPI_offset[nbZone] ++;
//             if (GASPI_offset[nbZone] - GASPI_commed[nbZone] >= GASPI_blockSize)
//             {
//                 needComm = true;
//                 GASPI_commed[nbZone] += GASPI_blockSize;        
//             }
            
//             omp_unset_lock(&GASPI_lock[nbZone]);

//             GASPI_sendOffset[nbZone*maxLength+offset] = ghost;
//             for(int j = 0; j < nVar; j++)
//                 GASPI_sendData[nbZone*maxLength*20+offset*nVar+j] = q[j][bCNo[neighbor][No]];
                
//             if (needComm)
//             {
//                 GASPI_wait_for_queue_half_full (0);
//                 SUCCESS_OR_DIE(gaspi_write( GASPI_sendDataSegmentId, (nbZone*maxLength*20+commed * nVar )* sizeof(double), nbZone, 
//                                             GASPI_recvDataSegmentId, (GASPI_iProc*maxLength * 20 + commed * nVar) * sizeof(double), 
//                                             GASPI_blockSize * nVar * sizeof(double), 0, GASPI_BLOCK ));
                
//                 // GASPI_wait_for_queue_half_full (0);
//                 SUCCESS_OR_DIE(gaspi_write_notify( GASPI_sendOffsetSegmentId, (nbZone*maxLength+commed) * sizeof(int), nbZone, 
//                                                     GASPI_recvOffsetSegmentId, (GASPI_iProc*maxLength + commed) * sizeof(int), 
//                                                     GASPI_blockSize * sizeof (int), commed / GASPI_blockSize * GASPI_nProc + nbZone,
//                                                     GASPI_blockSize, 0, GASPI_BLOCK ));
                
//                 commed+=GASPI_blockSize;
//             }
            
//             if (offset + 1 == nZIFace[neighbor] && nZIFace[neighbor] % GASPI_blockSize > 0)
//             {
//                 // GASPI_wait_for_queue_half_full (0);
//                 SUCCESS_OR_DIE(gaspi_write( GASPI_sendDataSegmentId, (nbZone*maxLength*20+commed * nVar )* sizeof(double), nbZone, 
//                                             GASPI_recvDataSegmentId, (GASPI_iProc*maxLength * 20 + commed * nVar) * sizeof(double),  
//                                             nZIFace[neighbor] % GASPI_blockSize * nVar * sizeof (double), 0, GASPI_BLOCK ));
                
//                 // GASPI_wait_for_queue_half_full (0);
//                 SUCCESS_OR_DIE(gaspi_write_notify( GASPI_sendOffsetSegmentId, (nbZone*maxLength+commed) * sizeof(int), nbZone, 
//                                                     GASPI_recvOffsetSegmentId, (GASPI_iProc*maxLength + commed) * sizeof(int),
//                                                     nZIFace[neighbor] % GASPI_blockSize * sizeof (int), 
//                                                     commed / GASPI_blockSize * GASPI_nProc + nbZone,
//                                                     nZIFace[neighbor] % GASPI_blockSize, 0, GASPI_BLOCK ));
                
//                 // GASPI_need_Comm[nbZone]=0;
//             }            
//         }
//     }
// }

// void PolyGrid::GASPIComm(RealFlow **q, int *boundryCell, int boundryNbCell)
// {
//     pair<int,int> **C2B = this->GetC2B();
//     IntType *nC2B = this->GetnC2B();
//     IntType **GHOST = this->GetGhost();
//     IntType nNo, No, i;
//     IntType nVar = 5;

//     // IntType Dcell[nNeighbor][boundryNbCell];//中间数组
//     IntType n[nNeighbor];//计数
//     // cout <<GASPI_iProc <<endl;
//     for (i = 0 ; i < nNeighbor; i++){
//         n[i] = 0;
//     }
 
//     IntType ilu, cell, ghost;
//     for(ilu=0;ilu<boundryNbCell;ilu++){
//         cell = boundryCell[ilu];
//         for(i=0;i<nC2B[cell];i++)
//         {
//             int neighbor = C2B[cell][i].first;
//             No = C2B[cell][i].second;
//             // Dcell[neighbor][n[neighbor]] = No;
//             n[neighbor]++;
//             // ghost = GHOST[neighbor][No];

//             // bool needComm = false;
//             // int nbZone = nb[neighbor];

//             // omp_set_lock(&GASPI_lock[nbZone]);

//             // int offset = GASPI_offset[nbZone];
//             // int commed = GASPI_commed[nbZone];
//             // GASPI_offset[nbZone] ++;
//             // if (GASPI_offset[nbZone] - GASPI_commed[nbZone] >= GASPI_blockSize)
//             // {
//             //     needComm = true;
//             //     GASPI_commed[nbZone] += GASPI_blockSize;        
//             // }
            
//             // omp_unset_lock(&GASPI_lock[nbZone]);

//             // GASPI_sendOffset[nbZone*maxLength+offset] = ghost;
//             // for(int j = 0; j < nVar; j++)
//             //     GASPI_sendData[nbZone*maxLength*20+offset*nVar+j] = q[j][bCNo[neighbor][No]];
                
//             // if (needComm)
//             // {
//             //     GASPI_wait_for_queue_half_full (0);
//             //     SUCCESS_OR_DIE(gaspi_write( GASPI_sendDataSegmentId, (nbZone*maxLength*20+commed * nVar )* sizeof(double), nbZone, 
//             //                                 GASPI_recvDataSegmentId, (GASPI_iProc*maxLength * 20 + commed * nVar) * sizeof(double), 
//             //                                 GASPI_blockSize * nVar * sizeof(double), 0, GASPI_BLOCK ));
                
//             //     // GASPI_wait_for_queue_half_full (0);
//             //     SUCCESS_OR_DIE(gaspi_write_notify( GASPI_sendOffsetSegmentId, (nbZone*maxLength+commed) * sizeof(int), nbZone, 
//             //                                         GASPI_recvOffsetSegmentId, (GASPI_iProc*maxLength + commed) * sizeof(int), 
//             //                                         GASPI_blockSize * sizeof (int), commed / GASPI_blockSize * GASPI_nProc + nbZone,
//             //                                         GASPI_blockSize, 0, GASPI_BLOCK ));
                
//             //     commed+=GASPI_blockSize;
//             // }
            
//             // if (offset + 1 == nZIFace[neighbor] && nZIFace[neighbor] % GASPI_blockSize > 0)
//             // {
//             //     // GASPI_wait_for_queue_half_full (0);
//             //     SUCCESS_OR_DIE(gaspi_write( GASPI_sendDataSegmentId, (nbZone*maxLength*20+commed * nVar )* sizeof(double), nbZone, 
//             //                                 GASPI_recvDataSegmentId, (GASPI_iProc*maxLength * 20 + commed * nVar) * sizeof(double),  
//             //                                 nZIFace[neighbor] % GASPI_blockSize * nVar * sizeof (double), 0, GASPI_BLOCK ));
                
//             //     // GASPI_wait_for_queue_half_full (0);
//             //     SUCCESS_OR_DIE(gaspi_write_notify( GASPI_sendOffsetSegmentId, (nbZone*maxLength+commed) * sizeof(int), nbZone, 
//             //                                         GASPI_recvOffsetSegmentId, (GASPI_iProc*maxLength + commed) * sizeof(int),
//             //                                         nZIFace[neighbor] % GASPI_blockSize * sizeof (int), 
//             //                                         commed / GASPI_blockSize * GASPI_nProc + nbZone,
//             //                                         nZIFace[neighbor] % GASPI_blockSize, 0, GASPI_BLOCK ));
                
//             //     // GASPI_need_Comm[nbZone]=0;
//             // }            
//         }
//     }
//     int num=0;
//     for (i = 0 ; i < nNeighbor; i++){
//         num+=n[i];
//     }    
//     IntType Dcell[num];

//     IntType ni[nNeighbor];//计数
//     // cout <<GASPI_iProc <<endl;
//     for (i = 0 ; i < nNeighbor; i++){
//         ni[i] = 0;
//         // cout << "boudry num =" <<i << " "<<n[i]<<endl;
//     }
//     // cout<< "write data" << endl;
//     for(ilu=0;ilu<boundryNbCell;ilu++){
//         cell = boundryCell[ilu];
//         for(i=0;i<nC2B[cell];i++)
//         {
//             int neighbor = C2B[cell][i].first;
//             No = C2B[cell][i].second;
//             num = 0;
//             for (int j = 0; j < neighbor; j++){
//                 num += n[j];
//             }
//             Dcell[num + ni[neighbor]] = No;
//             // cout<< "r "<< GASPI_iProc <<" " <<num + ni[neighbor] << " "<< No <<endl;
//             ni[neighbor]++;
//         }
//     }    
//     // cout<< "get data"  << endl;
//     num=0;
//     for (int neighbor = 0 ; neighbor < nNeighbor; neighbor++){
//         // for (i = 0; i < n[neighbor]; i++)
//         if(n[neighbor])
//         {
//             // No = Dcell[num];
//             // // cout<< "g "<< GASPI_iProc <<" " <<num <<" " << No <<endl;
//             // num++;
//             // ghost = GHOST[neighbor][No];

//             bool needComm = false;
//             int nbZone = nb[neighbor];

//             omp_set_lock(&GASPI_lock[nbZone]);

//             int offset = GASPI_offset[nbZone];
//             int commed = GASPI_commed[nbZone];
//             GASPI_offset[nbZone] += n[neighbor];
//             if (GASPI_offset[nbZone] - GASPI_commed[nbZone] >= GASPI_blockSize)
//             {
//                 needComm = true;
//                 GASPI_commed[nbZone] += GASPI_blockSize;        
//             }
            
//             omp_unset_lock(&GASPI_lock[nbZone]);

//             for (i = 0; i < n[neighbor]; i++)
//             {
//                 No = Dcell[num];
//                 // cout<< "g "<< GASPI_iProc <<" " <<num <<" " << No <<endl;
//                 num++;
//                 GASPI_sendOffset[nbZone*maxLength+offset+i] = GHOST[neighbor][No];
//                 for(int j = 0; j < nVar; j++)
//                     GASPI_sendData[nbZone*maxLength*20+(offset+i)*nVar+j] = q[j][bCNo[neighbor][No]];
//             }  
//             if (needComm)
//             {
//                 GASPI_wait_for_queue_half_full (0);
//                 SUCCESS_OR_DIE(gaspi_write( GASPI_sendDataSegmentId, (nbZone*maxLength*20+commed * nVar )* sizeof(double), nbZone, 
//                                             GASPI_recvDataSegmentId, (GASPI_iProc*maxLength * 20 + commed * nVar) * sizeof(double), 
//                                             GASPI_blockSize * nVar * sizeof(double), 0, GASPI_BLOCK ));
                
//                 // GASPI_wait_for_queue_half_full (0);
//                 SUCCESS_OR_DIE(gaspi_write_notify( GASPI_sendOffsetSegmentId, (nbZone*maxLength+commed) * sizeof(int), nbZone, 
//                                                     GASPI_recvOffsetSegmentId, (GASPI_iProc*maxLength + commed) * sizeof(int), 
//                                                     GASPI_blockSize * sizeof (int), commed / GASPI_blockSize * GASPI_nProc + nbZone,
//                                                     GASPI_blockSize, 0, GASPI_BLOCK ));
                
//                 commed+=GASPI_blockSize;
//             }
            
//             if (offset + n[neighbor] == nZIFace[neighbor] && nZIFace[neighbor] % GASPI_blockSize > 0)
//             {
//                 // GASPI_wait_for_queue_half_full (0);
//                 SUCCESS_OR_DIE(gaspi_write( GASPI_sendDataSegmentId, (nbZone*maxLength*20+commed * nVar )* sizeof(double), nbZone, 
//                                             GASPI_recvDataSegmentId, (GASPI_iProc*maxLength * 20 + commed * nVar) * sizeof(double),  
//                                             nZIFace[neighbor] % GASPI_blockSize * nVar * sizeof (double), 0, GASPI_BLOCK ));
                
//                 // GASPI_wait_for_queue_half_full (0);
//                 SUCCESS_OR_DIE(gaspi_write_notify( GASPI_sendOffsetSegmentId, (nbZone*maxLength+commed) * sizeof(int), nbZone, 
//                                                     GASPI_recvOffsetSegmentId, (GASPI_iProc*maxLength + commed) * sizeof(int),
//                                                     nZIFace[neighbor] % GASPI_blockSize * sizeof (int), 
//                                                     commed / GASPI_blockSize * GASPI_nProc + nbZone,
//                                                     nZIFace[neighbor] % GASPI_blockSize, 0, GASPI_BLOCK ));
                
//                 // GASPI_need_Comm[nbZone]=0;
//             }   
//             /* code */
//         }
//     }
//     // cout<< "end" << endl;
// }

void PolyGrid::GASPIComm(RealFlow **q, int *n, int *Dcell)
{
    pair<int,int> **C2B = this->GetC2B();
    IntType *nC2B = this->GetnC2B();
    IntType **GHOST = this->GetGhost();
    IntType nNo, No, i;
    IntType nVar = 5;

    int num=0;

    for (int neighbor = 0 ; neighbor < nNeighbor; neighbor++){
        // for (i = 0; i < n[neighbor]; i++)
        if(n[neighbor])
        {
            // No = Dcell[num];
            // // cout<< "g "<< GASPI_iProc <<" " <<num <<" " << No <<endl;
            // num++;
            // ghost = GHOST[neighbor][No];

            bool needComm = false;
            int nbZone = nb[neighbor];

            omp_set_lock(&GASPI_lock[nbZone]);

            int offset = GASPI_offset[nbZone];
            int commed = GASPI_commed[nbZone];
            GASPI_offset[nbZone] += n[neighbor];
            if (GASPI_offset[nbZone] - GASPI_commed[nbZone] >= GASPI_blockSize)
            {
                needComm = true;
                GASPI_commed[nbZone] += GASPI_blockSize;        
            }
            
            omp_unset_lock(&GASPI_lock[nbZone]);

            for (i = 0; i < n[neighbor]; i++)
            {
                No = Dcell[num];
                // cout<< "g "<< GASPI_iProc <<" " <<num <<" " << No <<endl;
                num++;
                GASPI_sendOffset[nbZone*maxLength+offset+i] = GHOST[neighbor][No];
                for(int j = 0; j < nVar; j++)
                    GASPI_sendData[nbZone*maxLength*20+(offset+i)*nVar+j] = q[j][bCNo[neighbor][No]];
            }  
            if (needComm)
            {
                GASPI_wait_for_queue_half_full (0);
                SUCCESS_OR_DIE(gaspi_write( GASPI_sendDataSegmentId, (nbZone*maxLength*20+commed * nVar )* sizeof(double), nbZone, 
                                            GASPI_recvDataSegmentId, (GASPI_iProc*maxLength * 20 + commed * nVar) * sizeof(double), 
                                            GASPI_blockSize * nVar * sizeof(double), 0, GASPI_BLOCK ));
                
                // GASPI_wait_for_queue_half_full (0);
                SUCCESS_OR_DIE(gaspi_write_notify( GASPI_sendOffsetSegmentId, (nbZone*maxLength+commed) * sizeof(int), nbZone, 
                                                    GASPI_recvOffsetSegmentId, (GASPI_iProc*maxLength + commed) * sizeof(int), 
                                                    GASPI_blockSize * sizeof (int), commed / GASPI_blockSize * GASPI_nProc + nbZone,
                                                    GASPI_blockSize, 0, GASPI_BLOCK ));
                
                commed+=GASPI_blockSize;
            }
            
            if (offset + n[neighbor] == nZIFace[neighbor] && nZIFace[neighbor] % GASPI_blockSize > 0)
            {
                // GASPI_wait_for_queue_half_full (0);
                SUCCESS_OR_DIE(gaspi_write( GASPI_sendDataSegmentId, (nbZone*maxLength*20+commed * nVar )* sizeof(double), nbZone, 
                                            GASPI_recvDataSegmentId, (GASPI_iProc*maxLength * 20 + commed * nVar) * sizeof(double),  
                                            nZIFace[neighbor] % GASPI_blockSize * nVar * sizeof (double), 0, GASPI_BLOCK ));
                
                // GASPI_wait_for_queue_half_full (0);
                SUCCESS_OR_DIE(gaspi_write_notify( GASPI_sendOffsetSegmentId, (nbZone*maxLength+commed) * sizeof(int), nbZone, 
                                                    GASPI_recvOffsetSegmentId, (GASPI_iProc*maxLength + commed) * sizeof(int),
                                                    nZIFace[neighbor] % GASPI_blockSize * sizeof (int), 
                                                    commed / GASPI_blockSize * GASPI_nProc + nbZone,
                                                    nZIFace[neighbor] % GASPI_blockSize, 0, GASPI_BLOCK ));
                
                // GASPI_need_Comm[nbZone]=0;
            }   
            /* code */
        }
    }
    // cout<< "end" << endl;
}

// void PolyGrid::GASPIComm(RealFlow **q, int *boundryCell, int boundryNbCell)
// {
//     pair<int,int> **C2B = this->GetC2B();
//     IntType *nC2B = this->GetnC2B();
//     IntType **ghost = this->GetGhost();
//     IntType nNo, no, i;
//     IntType nVar = 5;
//     IntType Dcell[nNeighbor][boundryNbCell];
//     IntType n[nNeighbor];
//     // cout <<GASPI_iProc <<endl;
//     for (i = 0 ; i < nNeighbor; i++){
//         n[i] = 0;
//     }

//     IntType ilu, cell;
//     for(ilu=0;ilu<boundryNbCell;ilu++){
//         cell = boundryCell[ilu];
//         for(i=0;i<nC2B[cell];i++)
//         {
//             int neighbor = C2B[cell][i].first;
//             no = C2B[cell][i].second;
//             Dcell[neighbor][n[neighbor]]=no;
//             // GHOST[neighbor][n[neighbor]] = ghost[neighbor][no];
//             // for(int j = 0; j < 5; j++)
//             //     Dq[neighbor][n[neighbor]*5+j]=q[j][bCNo[neighbor][No]]
//             n[neighbor]++;
//         }
//     }

//     for (int neighbor = 0; neighbor < nNeighbor;neighbor++){
//         if (n[neighbor]){

//             bool needComm = false;
//             int nbZone = nb[neighbor];
//             omp_set_lock(&GASPI_lock[nbZone]);

//             int offset = GASPI_offset[nbZone];
//             int commed = GASPI_commed[nbZone];
//             GASPI_offset[nbZone] += n[neighbor];

//             if (GASPI_offset[nbZone] - GASPI_commed[nbZone] >= GASPI_blockSize)
//             {
//                 needComm = true;
//                 GASPI_commed[nbZone] += GASPI_blockSize;        
//             }

//             omp_unset_lock(&GASPI_lock[nbZone]);

//             for (i = 0 ; i < n[neighbor] ; i++){
//                 int No=Dcell[neighbor][i];
//                 GASPI_sendOffset[nbZone*maxLength+offset+i] = ghost[neighbor][No];
//                 for(int j = 0; j < nVar; j++)
//                     GASPI_sendData[nbZone*maxLength*20+offset*nVar+j+i] = q[j][bCNo[neighbor][No]];
//             }
//             if (needComm)
//             {
//                 GASPI_wait_for_queue_half_full (0);
//                 SUCCESS_OR_DIE(gaspi_write( GASPI_sendDataSegmentId, (nbZone*maxLength*20+commed * nVar )* sizeof(double), nbZone, 
//                                             GASPI_recvDataSegmentId, (GASPI_iProc*maxLength * 20 + commed * nVar) * sizeof(double), 
//                                             GASPI_blockSize * nVar * sizeof(double), 0, GASPI_BLOCK ));
                
//                 GASPI_wait_for_queue_half_full (0);
//                 SUCCESS_OR_DIE(gaspi_write_notify( GASPI_sendOffsetSegmentId, (nbZone*maxLength+commed) * sizeof(int), nbZone, 
//                                                     GASPI_recvOffsetSegmentId, (GASPI_iProc*maxLength + commed) * sizeof(int), 
//                                                     GASPI_blockSize * sizeof (int), commed / GASPI_blockSize * GASPI_nProc + nbZone,
//                                                     GASPI_blockSize, 0, GASPI_BLOCK ));
                
//                 commed+=GASPI_blockSize;
//             }
            
//             if (offset + n[neighbor] == nZIFace[neighbor] && nZIFace[neighbor] % GASPI_blockSize > 0)
//             {
//                 GASPI_wait_for_queue_half_full (0);
//                 SUCCESS_OR_DIE(gaspi_write( GASPI_sendDataSegmentId, (nbZone*maxLength*20+commed * nVar )* sizeof(double), nbZone, 
//                                             GASPI_recvDataSegmentId, (GASPI_iProc*maxLength * 20 + commed * nVar) * sizeof(double),  
//                                             nZIFace[neighbor] % GASPI_blockSize * nVar * sizeof (double), 0, GASPI_BLOCK ));
                
//                 GASPI_wait_for_queue_half_full (0);
//                 SUCCESS_OR_DIE(gaspi_write_notify( GASPI_sendOffsetSegmentId, (nbZone*maxLength+commed) * sizeof(int), nbZone, 
//                                                     GASPI_recvOffsetSegmentId, (GASPI_iProc*maxLength + commed) * sizeof(int),
//                                                     nZIFace[neighbor] % GASPI_blockSize * sizeof (int), 
//                                                     commed / GASPI_blockSize * GASPI_nProc + nbZone,
//                                                     nZIFace[neighbor] % GASPI_blockSize, 0, GASPI_BLOCK ));
                
//                 // GASPI_need_Comm[nbZone]=0;
//             }
//         }        
//     }
// }

// void PolyGrid::GASPIComm(IntType nVar, RealFlow **q, IntType neighbor, int No, int ghost)
// {
//     bool needComm = false;
//     int nbZone = nb[neighbor];
//     // cout << "check messge :" << omp_test_lock(&GASPI_lock[nbZone]) <<endl;
//     // bool needlock = false; 
//     // if (GASPI_need_Comm[nbZone]==0)
//     // {
//     //     needlock = true;
//     //     GASPI_need_Comm[nbZone]++;
//     // }

//     // if(needlock){omp_set_lock(&GASPI_lock[nbZone]);}

//     int offset = GASPI_offset[nbZone];
//     int commed = GASPI_commed[nbZone];
//     GASPI_offset[nbZone] ++;
//     if (GASPI_offset[nbZone] - GASPI_commed[nbZone] >= GASPI_blockSize)
//     {
//         needComm = true;
//         GASPI_commed[nbZone] += GASPI_blockSize;        
//     }
    
//     // if (needlock){omp_unset_lock(&GASPI_lock[nbZone]);}

//     GASPI_sendOffset[nbZone*maxLength+offset] = ghost;
//     for(int j = 0; j < nVar; j++)
//         GASPI_sendData[nbZone*maxLength*20+offset*nVar+j] = q[j][bCNo[neighbor][No]];
        
//     if (needComm)
//     {
//         GASPI_wait_for_queue_half_full (0);
//         SUCCESS_OR_DIE(gaspi_write( GASPI_sendDataSegmentId, (nbZone*maxLength*20+commed * nVar )* sizeof(double), nbZone, 
//                                     GASPI_recvDataSegmentId, (GASPI_iProc*maxLength * 20 + commed * nVar) * sizeof(double), 
//                                     GASPI_blockSize * nVar * sizeof(double), 0, GASPI_BLOCK ));
        
//         // GASPI_wait_for_queue_half_full (0);
//         SUCCESS_OR_DIE(gaspi_write_notify( GASPI_sendOffsetSegmentId, (nbZone*maxLength+commed) * sizeof(int), nbZone, 
//                                             GASPI_recvOffsetSegmentId, (GASPI_iProc*maxLength + commed) * sizeof(int), 
//                                             GASPI_blockSize * sizeof (int), commed / GASPI_blockSize * GASPI_nProc + nbZone,
//                                             GASPI_blockSize, 0, GASPI_BLOCK ));
        
//         commed+=GASPI_blockSize;
//     }
    
//     if (offset + 1 == nZIFace[neighbor] && nZIFace[neighbor] % GASPI_blockSize > 0)
//     {
//         // GASPI_wait_for_queue_half_full (0);
//         SUCCESS_OR_DIE(gaspi_write( GASPI_sendDataSegmentId, (nbZone*maxLength*20+commed * nVar )* sizeof(double), nbZone, 
//                                     GASPI_recvDataSegmentId, (GASPI_iProc*maxLength * 20 + commed * nVar) * sizeof(double),  
//                                     nZIFace[neighbor] % GASPI_blockSize * nVar * sizeof (double), 0, GASPI_BLOCK ));
        
//         // GASPI_wait_for_queue_half_full (0);
//         SUCCESS_OR_DIE(gaspi_write_notify( GASPI_sendOffsetSegmentId, (nbZone*maxLength+commed) * sizeof(int), nbZone, 
//                                             GASPI_recvOffsetSegmentId, (GASPI_iProc*maxLength + commed) * sizeof(int),
//                                             nZIFace[neighbor] % GASPI_blockSize * sizeof (int), 
//                                             commed / GASPI_blockSize * GASPI_nProc + nbZone,
//                                             nZIFace[neighbor] % GASPI_blockSize, 0, GASPI_BLOCK ));
        
//         // GASPI_need_Comm[nbZone]=0;
//     }
// }

// void PolyGrid::GASPIComm(IntType nVar, RealFlow **q, IntType neighbor, int No)
// {
//     bool needComm = false;
//     int nbZone = nb[neighbor];

//     omp_set_lock(&GASPI_lock[nbZone]);

//     int offset = GASPI_offset[nbZone];
//     int commed = GASPI_commed[nbZone];
//     GASPI_offset[nbZone] ++;
//     if (GASPI_offset[nbZone] - GASPI_commed[nbZone] >= GASPI_blockSize)
//     {
//         needComm = true;
//         GASPI_commed[nbZone] += GASPI_blockSize;        
//     }

//     omp_unset_lock(&GASPI_lock[nbZone]);

//     GASPI_sendOffset[nbZone*maxLength+offset] = No;
//     for(int j = 0; j < nVar; j++)
//         GASPI_sendData[nbZone*maxLength*20+offset*nVar+j] = q[j][bCNo[neighbor][No]];
        
//     if (needComm)
//     {
//         // cout << "Send " << GASPI_iProc << " " << nbZone << " " << commed << " " <<  GASPI_sendData[nbZone*maxLength*20+commed*nVar] << endl << flush;
//         // printf("Send 2 %d %d %d %f\n", GASPI_iProc, nbZone, commed, GASPI_sendData[nbZone*maxLength*20+commed*nVar]);
//         GASPI_wait_for_queue_half_full (0);
//         SUCCESS_OR_DIE(gaspi_write( GASPI_sendDataSegmentId, (nbZone*maxLength*20+commed * nVar )* sizeof(double), nbZone, 
//                                     GASPI_recvDataSegmentId, (GASPI_iProc*maxLength * 20 + commed * nVar) * sizeof(double), 
//                                     GASPI_blockSize * nVar * sizeof(double), 0, GASPI_BLOCK ));
        
//         // GASPI_wait_for_queue_half_full (0);
//         SUCCESS_OR_DIE(gaspi_write_notify( GASPI_sendOffsetSegmentId, (nbZone*maxLength+commed) * sizeof(int), nbZone, 
//                                             GASPI_recvOffsetSegmentId, (GASPI_iProc*maxLength + commed) * sizeof(int), 
//                                             GASPI_blockSize * sizeof (int), commed / GASPI_blockSize * GASPI_nProc + nbZone,
//                                             GASPI_blockSize, 0, GASPI_BLOCK ));
        
//         commed+=GASPI_blockSize;
//     }
    
//     if (offset + 1 == nZIFace[neighbor] && nZIFace[neighbor] % GASPI_blockSize > 0)
//     {
//         // cout << "SendLast " << GASPI_iProc << " " << nbZone << " " << commed << " " << GASPI_sendData[nbZone*maxLength*20+commed*nVar] << endl << flush;
//         // printf("SendLast 2 %d %d %d %f\n", GASPI_iProc, nbZone, commed, GASPI_sendData[nbZone*maxLength*20+commed*nVar]);
//         // GASPI_wait_for_queue_half_full (0);
//         SUCCESS_OR_DIE(gaspi_write( GASPI_sendDataSegmentId, (nbZone*maxLength*20+commed * nVar )* sizeof(double), nbZone, 
//                                     GASPI_recvDataSegmentId, (GASPI_iProc*maxLength * 20 + commed * nVar) * sizeof(double),  
//                                     nZIFace[neighbor] % GASPI_blockSize * nVar * sizeof (double), 0, GASPI_BLOCK ));
        
//         // GASPI_wait_for_queue_half_full (0);
//         SUCCESS_OR_DIE(gaspi_write_notify( GASPI_sendOffsetSegmentId, (nbZone*maxLength+commed) * sizeof(int), nbZone, 
//                                             GASPI_recvOffsetSegmentId, (GASPI_iProc*maxLength + commed) * sizeof(int),
//                                             nZIFace[neighbor] % GASPI_blockSize * sizeof (int), 
//                                             commed / GASPI_blockSize * GASPI_nProc + nbZone,
//                                             nZIFace[neighbor] % GASPI_blockSize, 0, GASPI_BLOCK ));
//     }
// }

// void PolyGrid::checkMessage(IntType nVar, RealFlow **q)
// {
//     gaspi_notification_id_t notifyID;
//     gaspi_notification_t notifyValue;
//     IntType nNeighbor = this->GetNumberOfFaceNeighbors();
//     do
//     {
//         notifyID = 0;
//         notifyValue = 0;

//         for (int g = 0; g < nNeighbor; g++)
//         {
//             int nbZone = nb[g];
//             if (gaspi_notify_waitsome (GASPI_recvOffsetSegmentId,
//                 0,
//                 GASPI_nbBlock * GASPI_nProc,
//                 &notifyID,
//                 GASPI_TEST) != GASPI_SUCCESS)
//                     continue;

//             SUCCESS_OR_DIE (gaspi_notify_reset (GASPI_recvOffsetSegmentId,
//                 notifyID,
//                 &notifyValue));
                
//             if (notifyValue) 
//             {
// // #pragma omp critical
//                 GASPI_recvNotifications[nbZone]++;
//                 notifyID /= GASPI_nProc;
//                 GASPIUnzip(nVar, q, g, notifyID*GASPI_blockSize, notifyValue);
//                 // cout << "recv " << nbZone << " " << GASPI_iProc << " " << notifyID * GASPI_blockSize << endl << flush;
//                 break;
//             }
//         }
//     } while(notifyValue);
// }

// void PolyGrid::GASPIComm_Node(IntType nVar, RealFlow **q, IntType neighbor, int No)
// {
//     // cout << "begin" << flush;
//     bool needComm = false;
//     int nbZone = nbN[neighbor];

//     omp_set_lock(&GASPI_lock[nbZone]);

//     int offset = GASPI_offset[nbZone];
//     int commed = GASPI_commed[nbZone];
//     GASPI_offset[nbZone] ++;
//     if (GASPI_offset[nbZone] - GASPI_commed[nbZone] >= GASPI_blockSize)
//     {
//         needComm = true;
//         GASPI_commed[nbZone] += GASPI_blockSize;        
//     }

//     omp_unset_lock(&GASPI_lock[nbZone]);

//     GASPI_sendOffset[nbZone*maxLength+offset] = No;
//     for(int j = 0; j < nVar; j++)
//         GASPI_sendData[nbZone*maxLength*20+offset*nVar+j] = q[j][bNSNo[neighbor][No]];
       
//     if (needComm)
//     {
//         GASPI_wait_for_queue_half_full (0);
//         // cout << "Send " << GASPI_iProc << " " << nbZone << " " << commed << endl << flush;
//         SUCCESS_OR_DIE(gaspi_write( GASPI_sendDataSegmentId, (nbZone * maxLength * 20 + commed * nVar )* sizeof(double), nbZone, 
//                                     GASPI_recvDataSegmentId, (nbZone * maxLength * 20 + commed * nVar) * sizeof(double), 
//                                     GASPI_blockSize * nVar * sizeof (double), 0, GASPI_BLOCK ));
        
//         GASPI_wait_for_queue_half_full (0);
//         SUCCESS_OR_DIE(gaspi_write_notify( GASPI_sendOffsetSegmentId, (nbZone*maxLength+commed) * sizeof(int), nbZone, 
//                                             GASPI_recvOffsetSegmentId, (GASPI_iProc*maxLength + commed)  * sizeof(int), 
//                                             GASPI_blockSize * sizeof (int), commed / GASPI_blockSize * GASPI_nProc + nbZone,
//                                             GASPI_blockSize, 0, GASPI_BLOCK ));
        
//         commed+=GASPI_blockSize;
//     }
    
//     if (offset + 1 == nZINode[neighbor] && nZINode[neighbor] % GASPI_blockSize > 0)
//     {
//         GASPI_wait_for_queue_half_full (0);
//         // cout << "Send " << GASPI_iProc << " " << nbZone << " " << commed << endl << flush;
//         SUCCESS_OR_DIE(gaspi_write( GASPI_sendDataSegmentId, (nbZone*maxLength*20+commed * nVar )* sizeof(double), nbZone, 
//                                     GASPI_recvDataSegmentId, (GASPI_iProc*maxLength * 20 + commed * nVar) * sizeof(double), 
//                                     nZINode[neighbor] % GASPI_blockSize * nVar * sizeof (double), 0, GASPI_BLOCK ));
        
//         GASPI_wait_for_queue_half_full (0);
//         SUCCESS_OR_DIE(gaspi_write_notify( GASPI_sendOffsetSegmentId, (nbZone*maxLength+commed) * sizeof(int), nbZone, 
//                                             GASPI_recvOffsetSegmentId, (GASPI_iProc*maxLength + commed) * sizeof(int), 
//                                             nZINode[neighbor] % GASPI_blockSize * sizeof (int), 
//                                             commed / GASPI_blockSize * GASPI_nProc + nbZone,
//                                             nZINode[neighbor] % GASPI_blockSize, 0, GASPI_BLOCK ));
//     }
// }

// #elif TDTREE1

// void PolyGrid::GASPIWaitAll(IntType nVar, RealFlow **q)
// {
//     SUCCESS_OR_DIE(gaspi_wait(0,GASPI_BLOCK));
//     for(int i = 0; i < GASPI_notifications; i++)
//     {
//         gaspi_notification_id_t notifyID = 0;
//         gaspi_notification_t notifyValue = 0;

//         int j = i % GASPI_nProc;
//         while(1)
//         {
//             j = (j+1) % GASPI_nProc;
//             if (j == GASPI_iProc || gaspi_notify_waitsome (GASPI_recvOffsetSegmentId,
//                 0,
//                 GASPI_nbBlock * GASPI_nProc,
//                 &notifyID,
//                 GASPI_TEST) != GASPI_SUCCESS)
//                     continue;

//             SUCCESS_OR_DIE (gaspi_notify_reset (GASPI_recvOffsetSegmentId,
//                 notifyID,
//                 &notifyValue));
                
//             if (notifyValue) break;
//         }
//         // cout << i << " " << GASPI_notifications << endl;
//         int sendID = j;
//         notifyID /= GASPI_nProc;
//         GASPIUnzip(nVar, q, GASPI_mapping[sendID], notifyID*GASPI_blockSize, notifyValue);
//     }
    
//     SUCCESS_OR_DIE (gaspi_barrier (GASPI_GROUP_ALL, GASPI_BLOCK));
//     for(int i=0;i<GASPI_nProc;i++)
//         GASPI_offset[i] = GASPI_commed[i] = 0 ;

//     GASPI_wait_for_queue_half_full (0);
//     // MPI_Barrier(MPI_COMM_WORLD);
// }

// void PolyGrid::GASPIUnzip(IntType nVar, RealFlow **q, IntType neighbor, IntType ns, IntType ne)
// {
//     ne += ns;
//     IntType j, k, ghost;
//     int nbZone = nb[neighbor];
//     for(j=ns; j<ne; j++) {
        
//         int No = GASPI_recvOffset[nbZone*maxLength+j];
//         ghost    = nTCell + bFNo[neighbor][No];
//         for(k=0; k<nVar; k++) {
//             q[k][ghost] = GASPI_recvData[nbZone*maxLength*20+j*nVar+k];
//         }
//     }
// }
// void PolyGrid::GASPIComm(IntType nVar, RealFlow **q, IntType neighbor, int No)
// {
//     bool needComm = false;
//     int nbZone = nb[neighbor];

//     omp_set_lock(&GASPI_lock[nbZone]);

//     int offset = GASPI_offset[nbZone];
//     int commed = GASPI_commed[nbZone];
//     GASPI_offset[nbZone] ++;
//     if (GASPI_offset[nbZone] - GASPI_commed[nbZone] >= GASPI_blockSize)
//     {
//         needComm = true;
//         GASPI_commed[nbZone] += GASPI_blockSize;        
//     }

//     omp_unset_lock(&GASPI_lock[nbZone]);

//     GASPI_sendOffset[nbZone*maxLength+offset] = No;
//     for(int j = 0; j < nVar; j++)
//         GASPI_sendData[nbZone*maxLength*20+offset*nVar+j] = q[j][bCNo[neighbor][No]];

//     if (needComm)
//     {
//         GASPI_wait_for_queue_half_full (0);
//         // cout << commed*nVar << " " << maxLength * 20 + commed * nVar << endl << flush;
//         SUCCESS_OR_DIE(gaspi_write( GASPI_sendDataSegmentId, (nbZone * maxLength * 20+commed * nVar) * sizeof(double), nbZone, 
//                                     GASPI_recvDataSegmentId, (GASPI_iProc * maxLength * 20 + commed * nVar) * sizeof(double), 
//                                     GASPI_blockSize*nVar*sizeof(double), 0, GASPI_BLOCK ));
//         GASPI_wait_for_queue_half_full (0);
        
//         SUCCESS_OR_DIE(gaspi_write_notify( GASPI_sendOffsetSegmentId, (nbZone * maxLength + commed) * sizeof(int), nbZone, 
//                                             GASPI_recvOffsetSegmentId, (GASPI_iProc * maxLength + commed) * sizeof(int), 
//                                             GASPI_blockSize * sizeof (int), commed / GASPI_blockSize * GASPI_nProc + nbZone,
//                                             GASPI_blockSize, 0, GASPI_BLOCK ));
        
//         commed+=GASPI_blockSize;
//     }
    
//     if (offset + 1 == nZIFace[neighbor] && nZIFace[neighbor] % GASPI_blockSize > 0)
//     {
//         GASPI_wait_for_queue_half_full (0);
//         SUCCESS_OR_DIE(gaspi_write( GASPI_sendDataSegmentId, (nbZone * maxLength * 20+commed * nVar) * sizeof(double), nbZone, 
//                                     GASPI_recvDataSegmentId, (GASPI_iProc * maxLength * 20 + commed * nVar) * sizeof(double), 
//                                     nZIFace[neighbor] % GASPI_blockSize * nVar * sizeof (double), 0, GASPI_BLOCK ));
//         GASPI_wait_for_queue_half_full (0);
//         SUCCESS_OR_DIE(gaspi_write_notify( GASPI_sendOffsetSegmentId, (nbZone*maxLength+commed) * sizeof(int), nbZone, 
//                                             GASPI_recvOffsetSegmentId, (GASPI_iProc*maxLength + commed) * sizeof(int), 
//                                             nZIFace[neighbor] % GASPI_blockSize * sizeof (int), 
//                                             commed / GASPI_blockSize * GASPI_nProc + nbZone,
//                                             nZIFace[neighbor] % GASPI_blockSize, 0, GASPI_BLOCK ));
//     }
    
// }
// #elif TDTREE1

// void PolyGrid::GASPIWaitAll(IntType nVar, RealFlow **q)
// {
//     SUCCESS_OR_DIE(gaspi_wait(0,GASPI_BLOCK));
//     SUCCESS_OR_DIE (gaspi_barrier (GASPI_GROUP_ALL, GASPI_BLOCK)); 
//     IntType ghost;
//     for (int i = 0; i < nNeighbor; i++)
//     {
//         int nbZone = nb[i];
//         for (int j = 0; j < nZIFace[i]; j++)
//         {
//             ghost = nTCell + bFNo[i][j];
//             for (int k = 0 ; k < nVar ; k++)
//             {
//                 q[k][ghost] = GASPI_recvData[j*nVar+k+maxLength*20*nbZone];
//             }
//         }
//     }
//     for(int i=0;i<GASPI_nProc;i++)
//         GASPI_offset[i] = 0 ;

// }

// void PolyGrid::GASPIComm(IntType nVar, RealFlow **q, IntType neighbor, int No)
// {
//     bool needComm = false;
//     int nbZone = nb[neighbor];

//     omp_set_lock(&GASPI_lock[nbZone]);

//     int offset = GASPI_offset[nbZone];
//     GASPI_offset[nbZone] ++;
//     if (GASPI_offset[nbZone] == nZIFace[neighbor])
//     {
//         needComm = true;    
//     }
//     omp_unset_lock(&GASPI_lock[nbZone]);  

//     for(int j = 0; j < nVar; j++)
//     {
//         GASPI_sendData[nbZone*maxLength*20+offset*nVar+j] = q[j][bCNo[neighbor][No]];
//     }

//     if (needComm)
//     {
//         SUCCESS_OR_DIE(gaspi_write( GASPI_sendDataSegmentId, (nbZone * maxLength * 20) * sizeof(double), nbZone, 
//                                     GASPI_recvDataSegmentId, (GASPI_iProc * maxLength * 20 ) * sizeof(double), 
//                                     nZIFace[neighbor]*nVar*sizeof(double), 0, GASPI_BLOCK ));
//     }
    

    
// }
#endif
/*!
 * @brief       gaspi非阻塞发送，非阻塞接收
 * @param       nvar
 * @param       q
 *
 * @author      zz
 * @date        2023-12-18
 */
void PolyGrid::GASPIRecvSendVarNeighbor_Togeth(IntType nvar, RealFlow **q)
{
#ifdef MF_TIMING
    double time_tmp = -MPI_Wtime();
#endif

    if (nNeighbor == 0)
        return;

    IntType i, j, k;
// #ifdef TDTREE
// #pragma omp parallel for
// #endif
    for (i = 0; i < nNeighbor; i++)
    {
        int nbZone = nb[i];
        for (j = 0; j < nZIFace[i]; j++)
        {
            for (k = 0; k < nvar; k++)
                GASPI_sendData[maxLength*20*nbZone+j*nvar+k] = q[k][bCNo[i][j]];
        }
    }

    GASPIRecvSendVarNeighbor_Over(nvar);
    gaspi_barrier (GASPI_GROUP_ALL, GASPI_BLOCK);

    IntType ghost;
// #ifdef TDTREE
// #pragma omp parallel for
// #endif
    for (i = 0; i < nNeighbor; i++)
    {
        int nbZone = nb[i];
        for (j = 0; j < nZIFace[i]; j++)
        {
            ghost = nTCell + bFNo[i][j];
            for (k = 0 ; k < nvar ; k++)
            {
                q[k][ghost] = GASPI_recvData[j*nvar+k+maxLength*20*nbZone];
            }
        }
    }

    // for (int i = 0 ; i < 20*maxLength*GASPI_nProc; i++)
    // {
    //     GASPI_sendData[i] = 100.;
    //     GASPI_recvData[i] = 100.;
    // }
#ifdef MF_TIMING
    total_t[4] += time_tmp + MPI_Wtime();
#endif
}
/*!
 * @brief
 * @param       nvar
 *
 * @author      zz
 * @date        2023-12-18
 */
void PolyGrid::GASPIRecvSendVarNeighbor_Over(IntType nvar)
{
    IntType i, nbZone;
            // // send to other
// #ifdef TDTREE
// #pragma omp parallel for
// #endif
    for (i = 0; i < nNeighbor; i++)
    {
        nbZone = nb[i];
        SUCCESS_OR_DIE(gaspi_write( GASPI_sendDataSegmentId, maxLength * 20 * nbZone * sizeof(double), nbZone, 
                            GASPI_recvDataSegmentId, (maxLength * 20 * GASPI_iProc) * sizeof(double), 
                            nZIFace[i] * nvar * sizeof (double), 0, GASPI_BLOCK ));
    }
    SUCCESS_OR_DIE(gaspi_wait(0,GASPI_BLOCK));
}
#endif

/*!
 * @brief       建立边界点的传值机制
 * @note        modify according to the fun [void PolyGrid::SetUpComm_Node()]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-23
 */

void PolyGrid::SetUpComm_Node()
{
    IntType n, i, g, np;
    IntType nbZone;
    MPI_Status status;
    MPI_Request *req_send, *req_recv;
    MPI_Status *status_array;

    IntType nNeighborN = this->GetNumberOfNodeNeighbors();

    if (nNeighborN == 0)
        return;
    status_array = NULL;
    req_send = NULL;
    req_recv = NULL;
    // nZINode = NULL;
    // bNSNo = NULL;
    snew_array_1D(status_array, nNeighborN);
    snew_array_1D(req_send, nNeighborN);
    snew_array_1D(req_recv, nNeighborN);
    snew_array_1D(nZINode, nNeighborN);
    snew_array_1D(bNSNo, nNeighborN);

    IntType **bqr = NULL;
    IntType **bqs = NULL;
    snew_array_1D(bqr, nNeighborN);
    snew_array_1D(bqs, nNeighborN);
#if (defined GASPI) && (defined TDTREE)
    // cout << "SetUpComm_Node bigin" << endl << flush;
    IntType nTNode = this->GetNTNode();
    vector<pair<int,int>> *n2b = new vector<pair<int,int>>[nTNode]();
    nN2B = new IntType[nTNode]();
    N2B = new pair<int,int>*[nTNode]();
#endif 
    for (g = 0; g < nNeighborN; g++)
    {
        nbZone = nbN[g];

        // 计数与其相邻的各块分别有多少个并行边点
        n = 0;
        for (i = 0; i < nINode; i++)
        {
            if (nbZN[i] == nbZone)
            {
                n++;
            }
        }
        nZINode[g] = n;

        // 将当前块中的需要传递的点号储存在bNSNo[g][n]中:g代表将要传给第几块,n表示序列号
        bNSNo[g] = NULL;
        bqs[g] = NULL;
        snew_array_1D(bNSNo[g], n);
        snew_array_1D(bqs[g], n);
        n = 0;
        for (i = 0; i < nINode; i++)
        {
            if (nbZN[i] == nbZone)
            {
                bNSNo[g][n] = nbSN[i];
#if (defined GASPI) && (defined TDTREE)
                n2b[nbSN[i]].push_back(make_pair(g, n));
#endif    
                bqs[g][n++] = nbRN[i];
            }
        }
    }

    // now receiving
    for (g = 0; g < nNeighborN; g++)
    {
        nbZone = nbN[g];
        MPI_Send(&nZINode[g], 1, MPIIntType, nbZone, level, MPI_COMM_WORLD);
    }
    for (g = 0; g < nNeighborN; g++)
    {
        nbZone = nbN[g];
        MPI_Recv(&n, 1, MPIIntType, nbZone, level, MPI_COMM_WORLD, &status);
        // assert(n == nZINode[g]);
    }

    for (g = 0; g < nNeighborN; g++)
    {
        bqr[g] = NULL;
        snew_array_1D(bqr[g], nZINode[g]);
    }
    for (np = 1; np <= numprocs; np++)
    {
        if (myZone == np)
        {
            //       send   to other
            for (g = 0; g < nNeighborN; g++)
            {
                nbZone = nbN[g];
                MPI_Isend(bqs[g], nZINode[g], MPIIntType, nbZone, level, MPI_COMM_WORLD, &req_send[g]);
            }
        }
        else
        {
            //       receive   from np
            for (g = 0; g < nNeighborN; g++)
            {
                nbZone = nbN[g];
                if (nbZone == (np - 1))
                {
                    MPI_Irecv(bqr[g], nZINode[g], MPIIntType, nbZone, level, MPI_COMM_WORLD, &req_recv[g]);
                }
            }
        }
    } // np

    MPI_Waitall(nNeighborN, req_recv, status_array);
    // bNRNo = NULL;
    snew_array_1D(bNRNo, nNeighborN);
    for (g = 0; g < nNeighborN; g++)
    {
        bNRNo[g] = NULL;
        snew_array_1D(bNRNo[g], nZINode[g]);
        for (i = 0; i < nZINode[g]; i++)
        {
            bNRNo[g][i] = bqr[g][i];
        }
    }
    MPI_Waitall(nNeighborN, req_send, status_array);

    for (g = 0; g < nNeighborN; g++)
    {
        sdel_array_1D(bqr[g]);
        sdel_array_1D(bqs[g]);
    }
#if (defined GASPI) && (defined TDTREE)
    // cout << "SetUpComm_Node bigin2" << endl << flush;
    GASPI_notifications_node = new IntType[nNeighborN];
    int aa=0;
    for(int g=0; g<nNeighborN; g++){
        nbZone = nbN[g];
        GASPI_mappingN[nbZone] = g;
        GASPI_notifications_node[g]= nZINode[g] / GASPI_blockSize + (nZINode[g] % GASPI_blockSize > 0); 
        aa+=GASPI_notifications_node[g];
    }
    cout << "GASPI_notifications_node =" << aa<< endl << flush;
    for(int i=0;i<nTNode;i++)
    {
        nN2B[i] = n2b[i].size();
        snew_array_1D(N2B[i], (int)n2b[i].size());
        for(int j=0;j<nN2B[i];j++)
            N2B[i][j] = n2b[i][j];
    }

    delete[] n2b;
    // cout<< "eeeeee"<< endl << flush;
    TDTree *TDTreeRoot = this->GetTDTree();
    TDTreeRoot->task_traversal(initBoundryEntity, NULL, (char **)this, Forward);
    cout << "init Finish!!!!!!!!" << endl << flush;
#endif   
    sdel_array_1D(bqr);
    sdel_array_1D(bqs);
    sdel_array_1D(req_send);
    sdel_array_1D(req_recv);
    sdel_array_1D(status_array);
}

/*!
 * @brief       并行传值
 * @param       q
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-27
 */
void PolyGrid::CommInterfaceDataMPI(RealFlow *q)
{
    if (nNeighbor == 0)
        return;
#ifdef MF_TIMING
    MPI_Barrier(MPI_COMM_WORLD);
    double time_tmp = -MPI_Wtime();
#endif
    if (q)
        RecvSendVarNeighbor(q);
#ifdef MF_TIMING
    total_t[4] += time_tmp + MPI_Wtime();
#endif
}

/*!
 * @brief       非阻塞发送，非阻塞接收
 * @param       q
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-06-27
 */
void PolyGrid::RecvSendVarNeighbor(RealFlow *q)
{
    IntType i, j, nbZone, ghost;
    IntType np;
    RealFlow **bqs = NULL;
    RealFlow **bqr = NULL;

    MPI_Request *req_send, *req_recv;
    MPI_Status *status_array;

    status_array = NULL;
    req_send = NULL;
    req_recv = NULL;
    snew_array_1D(status_array, nNeighbor);
    snew_array_1D(req_send, nNeighbor);
    snew_array_1D(req_recv, nNeighbor);
    snew_array_2D(bqr, nNeighbor, nZIFace, false);
    snew_array_2D(bqs, nNeighbor, nZIFace, false);

    // t1 = MPI_Wtime();

    for (np = 1; np <= numprocs; np++)
    {
        if (myZone == np)
        {
            // send to other
            for (i = 0; i < nNeighbor; i++)
            {
                nbZone = nb[i];
                for (j = 0; j < nZIFace[i]; j++)
                {
                    bqs[i][j] = q[bCNo[i][j]];
                }
                MPI_Isend(bqs[i], nZIFace[i], MPIReal, nbZone, level, MPI_COMM_WORLD, &req_send[i]);
            }
        }
        else
        {
            // receive from np
            for (i = 0; i < nNeighbor; i++)
            {
                nbZone = nb[i];
                if (nbZone == (np - 1))
                {
                    MPI_Irecv(bqr[i], nZIFace[i], MPIReal, nbZone, level, MPI_COMM_WORLD, &req_recv[i]);
                }
            }
        }
    } // np

    MPI_Waitall(nNeighbor, req_recv, status_array);
    for (i = 0; i < nNeighbor; i++)
    {
        for (j = 0; j < nZIFace[i]; j++)
        {
            ghost = nTCell + bFNo[i][j];
            q[ghost] = bqr[i][j];
        }
    }
    MPI_Waitall(nNeighbor, req_send, status_array);

    // t2 = MPI_Wtime();
    // mpi_time += t2 - t1;
    sdel_array_2D(bqr, nNeighbor, false);
    sdel_array_2D(bqs, nNeighbor, false);
    sdel_array_1D(req_recv);
    sdel_array_1D(req_send);
    sdel_array_1D(status_array);
}

/*!
 * @brief       并行传值
 * @param       q
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-10
 */
void PolyGrid::CommInternodeDataMPI(RealFlow *q)
{
    if (nNeighborN == 0)
        return;
#ifdef MF_TIMING
    MPI_Barrier(MPI_COMM_WORLD);
    double time_tmp = -MPI_Wtime();
#endif
    if (q)
        RecvSendVarNeighbor_Node(q);
#ifdef MF_TIMING
    total_t[4] += time_tmp + MPI_Wtime();
#endif
}

/*!
 * @brief       非阻塞发送，非阻塞接收
 * @param       q
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-10
 */
void PolyGrid::RecvSendVarNeighbor_Node(RealFlow *q)
{
    IntType i, j, nbZone;
    IntType np;
    RealFlow **bqs = NULL;
    RealFlow **bqr = NULL;

    MPI_Request *req_send, *req_recv;
    MPI_Status *status_array;
    status_array = NULL;
    req_send = NULL;
    req_recv = NULL;
    snew_array_1D(status_array, nNeighborN);
    snew_array_1D(req_send, nNeighborN);
    snew_array_1D(req_recv, nNeighborN);
    snew_array_2D(bqr, nNeighborN, nZINode, false);
    snew_array_2D(bqs, nNeighborN, nZINode, false);
    // t1 = MPI_Wtime();

    for (np = 1; np <= numprocs; np++)
    {
        if (myZone == np)
        {
            // send to other
            for (i = 0; i < nNeighborN; i++)
            {
                nbZone = nbN[i];
                for (j = 0; j < nZINode[i]; j++)
                {
                    bqs[i][j] = q[bNSNo[i][j]];
                }
                MPI_Isend(bqs[i], nZINode[i], MPIReal, nbZone, level, MPI_COMM_WORLD, &req_send[i]);
            }
        }
        else
        {
            //  receive from np
            for (i = 0; i < nNeighborN; i++)
            {
                nbZone = nbN[i];
                if (nbZone == (np - 1))
                {
                    MPI_Irecv(bqr[i], nZINode[i], MPIReal, nbZone, level, MPI_COMM_WORLD, &req_recv[i]);
                }
            }
        }
    } // np

    MPI_Waitall(nNeighborN, req_recv, status_array);
    for (i = 0; i < nNeighborN; i++)
    {
        for (j = 0; j < nZINode[i]; j++)
        {
            q[bNRNo[i][j]] += bqr[i][j];
        }
    }
    MPI_Waitall(nNeighborN, req_send, status_array);

    // t2 = MPI_Wtime();
    // mpi_time += t2 - t1;

    sdel_array_2D(bqr, nNeighborN, false);
    sdel_array_2D(bqs, nNeighborN, false);
    sdel_array_1D(req_recv);
    sdel_array_1D(req_send);
    sdel_array_1D(status_array);
}

/*!
 * @brief       并行传值
 * @param       q
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-10
 */
void PolyGrid::CommInternodeDataMPI2(IntType *q)
{
    if (nNeighborN == 0)
        return;
#ifdef MF_TIMING
    MPI_Barrier(MPI_COMM_WORLD);
    double time_tmp = -MPI_Wtime();
#endif
    if (q)
        RecvSendVarNeighbor_Node2(q);
#ifdef MF_TIMING
    total_t[4] += time_tmp + MPI_Wtime();
#endif
}

/*!
 * @brief       node variable communication and synchronize to the minimum value of multi-nodes among processors
 * @param       q
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-10
 */
void PolyGrid::RecvSendVarNeighbor_Node2(IntType *q)
{
    IntType i, j, nbZone;
    IntType np;
    IntType **bqs = NULL;
    IntType **bqr = NULL;

    MPI_Request *req_send, *req_recv;
    MPI_Status *status_array;

    status_array = NULL;
    req_send = NULL;
    req_recv = NULL;
    snew_array_1D(status_array, nNeighborN);
    snew_array_1D(req_send, nNeighborN);
    snew_array_1D(req_recv, nNeighborN);
    snew_array_2D(bqr, nNeighborN, nZINode, false);
    snew_array_2D(bqs, nNeighborN, nZINode, false);

    // t1 = MPI_Wtime();

    for (np = 1; np <= numprocs; np++)
    {
        if (myZone == np)
        {
            // send to other
            for (i = 0; i < nNeighborN; i++)
            {
                nbZone = nbN[i];
                for (j = 0; j < nZINode[i]; j++)
                {
                    bqs[i][j] = q[bNSNo[i][j]];
                }
                MPI_Isend(bqs[i], nZINode[i], MPIIntType, nbZone, level, MPI_COMM_WORLD, &req_send[i]);
            }
        }
        else
        {
            // receive from np
            for (i = 0; i < nNeighborN; i++)
            {
                nbZone = nbN[i];
                if (nbZone == (np - 1))
                {
                    MPI_Irecv(bqr[i], nZINode[i], MPIIntType, nbZone, level, MPI_COMM_WORLD, &req_recv[i]);
                }
            }
        }
    } // np

    MPI_Waitall(nNeighborN, req_recv, status_array);
    for (i = 0; i < nNeighborN; i++)
    {
        for (j = 0; j < nZINode[i]; j++)
        {
            if (bqr[i][j] == 0)
                continue;
            if (q[bNRNo[i][j]] == 0)
            {
                q[bNRNo[i][j]] = bqr[i][j];
            }
            else if (q[bNRNo[i][j]] > bqr[i][j])
            {
                q[bNRNo[i][j]] = bqr[i][j];
            }
        }
    }
    MPI_Waitall(nNeighborN, req_send, status_array);

    // t2 = MPI_Wtime();
    // mpi_time += t2 - t1;

    sdel_array_2D(bqr, nNeighborN, false);
    sdel_array_2D(bqs, nNeighborN, false);
    sdel_array_1D(req_recv);
    sdel_array_1D(req_send);
    sdel_array_1D(status_array);
}

/*!
 * @brief       非阻塞发送，非阻塞接收
 * @param       nvar
 * @param       q
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-16
 */
void PolyGrid::RecvSendVarNeighbor_Togeth(IntType nvar, RealFlow **q)
{
#ifdef MF_TIMING
    MPI_Barrier(MPI_COMM_WORLD);
    double time_tmp = -MPI_Wtime();
#endif

    if (nNeighbor == 0)
        return;

    RealFlow ***bqs = NULL, ***bqr = NULL;
    IntType i;

    MPI_Request *req_send = 0, *req_recv = 0;
    MPI_Status *status_array = 0;

    status_array = NULL;
    req_send = NULL;
    req_recv = NULL;
    snew_array_1D(status_array, nNeighbor);
    snew_array_1D(req_send, nNeighbor);
    snew_array_1D(req_recv, nNeighbor);
    snew_array_1D(bqr, nNeighbor);
    snew_array_1D(bqs, nNeighbor);

    Set_RecvSend(bqs, bqr, nvar);

    for (i = 0; i < nvar; i++)
        Add_RecvSend(bqs, q[i], i);

    RecvSendVarNeighbor_Over(bqs, bqr, req_send, req_recv, status_array, nvar);

    MPI_Waitall(nNeighbor, req_recv, status_array);
    MPI_Waitall(nNeighbor, req_send, status_array);

    sdel_array_1D(req_send);
    sdel_array_1D(req_recv);
    sdel_array_1D(status_array);
    for (i = 0; i < nvar; i++)
        Read_RecvSend(bqr, q[i], i);

    sdel_array_1D(bqr[0][0]);
    sdel_array_1D(bqr[0]);
    sdel_array_1D(bqr);
    sdel_array_1D(bqs[0][0]);
    sdel_array_1D(bqs[0]);
    sdel_array_1D(bqs);

#ifdef MF_TIMING
    total_t[4] += time_tmp + MPI_Wtime();
#endif
}

/*!
 * @brief
 * @param       bqs
 * @param       bqr
 * @param       req_send
 * @param       req_recv
 * @param       status_array
 * @param       nvar
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-16
 */
void PolyGrid::RecvSendVarNeighbor_Over(RealFlow ***bqs, RealFlow ***bqr, MPI_Request *req_send,
                                        MPI_Request *req_recv, MPI_Status *status_array, IntType nvar)
{
    IntType i, nbZone;
    IntType np;

    for (np = 1; np <= numprocs; np++)
    {
        if (myZone == np)
        {
            // send to other
            for (i = 0; i < nNeighbor; i++)
            {
                nbZone = nb[i];
                MPI_Isend(bqs[i][0], nZIFace[i] * nvar, MPIReal, nbZone, level, MPI_COMM_WORLD, &req_send[i]);
            }
        }
        else
        {
            // receive from np
            for (i = 0; i < nNeighbor; i++)
            {
                nbZone = nb[i];
                if (nbZone == (np - 1))
                {
                    MPI_Irecv(bqr[i][0], nZIFace[i] * nvar, MPIReal, nbZone, level, MPI_COMM_WORLD, &req_recv[i]);
                }
            }
        }
    }
}

/*!
 * @brief
 * @param       bqs
 * @param       bqr
 * @param       nvar
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-16
 */
void PolyGrid::Set_RecvSend(RealFlow ***bqs, RealFlow ***bqr, IntType nvar)
{
    IntType i, j, k, temp_nZIFace = 0;

    for (i = 0; i < nNeighbor; i++)
        temp_nZIFace += nZIFace[i];
    bqs[0] = NULL;
    bqr[0] = NULL;
    snew_array_1D(bqs[0], nvar * nNeighbor);
    snew_array_1D(bqr[0], nvar * nNeighbor);

    for (i = 1; i < nNeighbor; i++)
    {
        bqs[i] = &bqs[i - 1][nvar];
        bqr[i] = &bqr[i - 1][nvar];
    }
    bqs[0][0] = NULL;
    bqr[0][0] = NULL;
    snew_array_1D(bqs[0][0], nNeighbor * nvar * temp_nZIFace);
    snew_array_1D(bqr[0][0], nNeighbor * nvar * temp_nZIFace);
    for (i = 1; i < nNeighbor; i++)
    {
        bqs[i][0] = &bqs[i - 1][0][nvar * nZIFace[i - 1]];
        bqr[i][0] = &bqr[i - 1][0][nvar * nZIFace[i - 1]];
    }
    for (i = 0; i < nNeighbor; i++)
    {
        for (j = 1; j < nvar; j++)
        {
            bqs[i][j] = &bqs[i][j - 1][nZIFace[i]];
            bqr[i][j] = &bqr[i][j - 1][nZIFace[i]];
        }
    }

    for (i = 0; i < nNeighbor; i++)
    {
        for (j = 0; j < nvar; j++)
        {
            for (k = 0; k < nZIFace[i]; k++)
            {
                bqs[i][j][k] = 0.0;
                bqr[i][j][k] = 0.0;
            }
        }
    }
}

/*!
 * @brief
 * @param       bqs
 * @param       q
 * @param       num_var
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-16
 */
void PolyGrid::Add_RecvSend(RealFlow ***bqs, RealFlow *q, IntType num_var)
{
    IntType i, j;
    for (i = 0; i < nNeighbor; i++)
    {
        for (j = 0; j < nZIFace[i]; j++)
        {
            bqs[i][num_var][j] = q[bCNo[i][j]];
        }
    }
}

/*!
 * @brief
 * @param       bqr
 * @param       q
 * @param       num_var
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-16
 */
void PolyGrid::Read_RecvSend(RealFlow ***bqr, RealFlow *q, IntType num_var)
{
    IntType i, j, ghost;
    for (i = 0; i < nNeighbor; i++)
    {
        for (j = 0; j < nZIFace[i]; j++)
        {
            ghost = nTCell + bFNo[i][j];
            q[ghost] = bqr[i][num_var][j];
        }
    }
}

#endif