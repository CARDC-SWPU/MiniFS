/*!
 * @file        parallel_omp.cpp
 * @brief
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @version     1.0
 * @date        2023-07-03
 * @copyright   Copyright (c) 2023  PERSONAL
 *
 * @par 修改日志:
 * <table>
 * <tr><th> Date        <th> Version  <th> Author  <th> Description
 * <tr><td> 2023-07-03  <td> 1.0      <td> Wisces  <td> 实现了支持共享内存并行编程模型的四种并行算法
 * </table>
 */

//!< direct head file
#include "grid_polyhedra.h"

//!< C/C++ head files
#include <iostream>
#include <set>
#include <cstdlib> ///< exit()
#include <cmath>   ///< ceil()
#include<fstream>

//!< user defined head files
#include "memory_util.h"
#include "constant.h"
#include "number_type.h"
// #include "para_field_global.h"

//!< head file relying on condition-compiling
#ifdef MF_MPICH
#include <mpi.h>
#endif

#ifdef MF_OPENMP
#include <omp.h>
#endif

#ifdef OMP_DIVREP
#include "metis.h"
#endif

#ifdef MF_MPICH
extern int myid;
extern int numprocs;
extern int myZone;        // = myid+1
extern MPI_Comm GridComm; ///< for each grid
#endif

using namespace std;

#ifdef OMP_FaceColor
/*!
 * @brief       Colouring faces for fine-grained parallelization (i.e., openmp and simd)
 * @param       grid
 * @param       bgroupsize
 * @param       igroupsize
 * @note        参考文献：张健，李瑞田，邓亮，代喆，刘杰，徐传福．面向多核 CPU/众核 GPU 架构的非结构 CFD 共享内存并行计算技术研究[J/OL]．航空学报.
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-03
 */
void FaceColoring(PolyGrid *grid, IntType bgroupsize, IntType igroupsize)
{
    set<IntType> cellset;
    IntType nTCell, nBFace, nTFace, nIFace;
    IntType *f2c, *f2n, *nNPF;

    nTCell = grid->GetNTCell();
    nBFace = grid->GetNBFace();
    nIFace = grid->GetNIFace();
    nTFace = grid->GetNTFace();
    f2c = grid->Getf2c();
    f2n = grid->Getf2n();
    nNPF = grid->GetnNPF();

    IntType *index_bface = NULL; ///< 物理边界面
    IntType *index_bface_1 = NULL;
    // IntType *nNPF_bface = NULL;
    // IntType *f2n_bface = NULL;
    // IntType **F2N_bface = NULL;

    IntType *index_iface = NULL; ///< 内部面
    IntType *index_iface_1 = NULL;

    IntType ifacenum = nTFace - nBFace; ///< 非边界面（内部面）数量

    IntType n = nTCell + nBFace;
    IntType pfacenum = nBFace - nIFace; ///< 物理边界面数量

    //!< allocate memory
    snew_array_1D(index_bface, pfacenum);
    snew_array_1D(index_bface_1, pfacenum);

    snew_array_1D(index_iface, ifacenum);
    snew_array_1D(index_iface_1, ifacenum);

    //!< initialization
    for (IntType i = 0; i < pfacenum; i++)
        index_bface[i] = -1;

    for (IntType i = nBFace; i < nTFace; i++)
        index_iface[i - nBFace] = -1;

    IntType newindex = 0, bgroupnum = 0, igroupnum = 0;
    cellset.clear();
    //!< 2.2 colouring boundary faces not including interface faces（并行边界面单独处理）
    IntType bgrouplength = 0;
    do
    {
        bgrouplength = 0;
        for (IntType i = 0; i < pfacenum; i++)
        {
            if (index_bface[i] == -1) ///< unmarked face
            {
                if (cellset.count(f2c[2 * i]) == 0) // unmarked c1 cell, cellset容器里还没有 f2c[2*i], 该重循环以 f2c[2*i] 为 reference cell
                {
                    index_bface[i] = newindex;   // 将 newindex 存入 index_bface, 表明该面已着色, 注意 i 仍有信息
                    index_bface_1[newindex] = i; // index_bface[i] 转置存储格式，用于重排序，将相同的颜色面放在一起
                    newindex = newindex + 1;     // 此面已着色，着色数量加一
                    cellset.insert(f2c[2 * i]);  // mark c1 cell
                    bgrouplength = bgrouplength + 1;
                    if (bgrouplength >= bgroupsize)
                        break;
                }
            }
        }
        // cout << "newindex: " << newindex << endl;

        //!< edpas(ngroup) = newindex - 1;
        (*grid).bfacegroup.push_back(newindex); ///< newindex is also ok, 表明该颜色 group 内 c1 cell 个数
        cellset.clear();
        bgroupnum = bgroupnum + 1; ///< 最终得到 bgroupnum 个颜色，对于边界面

    } while (newindex < pfacenum); ///< until all boundary faces marked, 每重 do 循环新建一个 color group

    //!< 2.3 update f2c, nNPF, f2n for boundary faces
    UpdateFaceData(grid, index_bface_1, 0);
    // IntType *f2c_backup = NULL;
    // snew_array_1D(f2c_backup, 2 * pfacenum);
    // for (IntType i = 0; i < pfacenum; i++)
    // {
    //     f2c_backup[2 * i] = f2c[2 * i];
    //     f2c_backup[2 * i + 1] = f2c[2 * i + 1];
    // }
    // //!< 2.3.1 update f2c for boundary faces
    // for (IntType i = 0; i < pfacenum; i++)
    // {
    //     f2c[2 * i] = f2c_backup[index_bface_1[i] * 2];
    //     f2c[2 * i + 1] = f2c_backup[index_bface_1[i] * 2 + 1];
    // }

    // // IntType bfirst = (*grid).bfacegroup[0];

    // //!< update nNPF, f2n for boundary faces
    // //!< prepare information and data structures
    // IntType *nNPF_bface = NULL;
    // IntType *f2n_bface = NULL;
    // snew_array_1D(nNPF_bface, pfacenum);
    // IntType m = 0;
    // for (IntType i = 0; i < pfacenum; i++)
    // {
    //     m += nNPF[i];
    // }
    // snew_array_1D(f2n_bface, m);

    // m = 0;
    // for (IntType i = 0; i < pfacenum; i++)
    // {
    //     nNPF_bface[i] = nNPF[i];
    //     for (IntType j = 0; j < nNPF[i]; j++)
    //     {
    //         f2n_bface[m] = f2n[m];
    //         m++;
    //     }
    // }

    // IntType **F2N_bface = NULL;
    // snew_array_1D(F2N_bface, pfacenum);
    // F2N_bface[0] = f2n_bface;
    // for (IntType i = 1; i < pfacenum; i++)
    // {
    //     F2N_bface[i] = &(F2N_bface[i - 1][nNPF_bface[i - 1]]);
    // }

    // //!< 2.3.2. update nNPF for physical faces
    // for (IntType i = 0; i < pfacenum; i++)
    // {
    //     nNPF[i] = nNPF_bface[index_bface_1[i]];
    // }

    // //!< 2.3.3. update f2n for physical faces
    // m = 0;
    // for (IntType i = 0; i < pfacenum; i++)
    // {
    //     for (IntType j = 0; j < nNPF[i]; j++)
    //     {
    //         f2n[m] = F2N_bface[index_bface_1[i]][j];
    //         m++;
    //     }
    // }

    //!< 2.4 colouring interior faces
    newindex = nBFace;
    cellset.clear();
    IntType igrouplength = 0;
    do
    {
        igrouplength = 0;
        for (IntType i = nBFace; i < nTFace; i++)
        {
            if (index_iface[i - nBFace] == -1) ///< unmarked face
            {
                if (cellset.count(f2c[2 * i]) == 0) ///< unmarked c1 cell
                {
                    if (cellset.count(f2c[2 * i + 1]) == 0) ///< unmarked c2 cel
                    {
                        index_iface[i - nBFace] = newindex;
                        index_iface_1[newindex - nBFace] = i;
                        newindex = newindex + 1;
                        cellset.insert(f2c[2 * i]); ///< mark c1 cell
                        cellset.insert(f2c[2 * i + 1]);
                        igrouplength = igrouplength + 1;
                        if (igrouplength >= igroupsize)
                            break;
                    }
                }
            }
        }
        // edpas(ngroup)=newindex-1;
        (*grid).ifacegroup.push_back(newindex);
        cellset.clear();
        igroupnum = igroupnum + 1;

    } while (newindex < nTFace); ///< until all interior faces marked

    //!< 2.5 update f2c, nNPF, f2n for interior faces
    UpdateFaceData(grid, index_iface_1, nBFace);

    IntType bfacegroup_num, ifacegroup_num;
    ifacegroup_num = (*grid).ifacegroup.size(); ///< 内部面的颜色数
    bfacegroup_num = (*grid).bfacegroup.size(); ///< 物理边界面的颜色数

    ofstream fp1 ,fp2 ;
    fp1.open( "face.txt" , ofstream::app|ofstream::out);
    for (IntType fcolor = 0; fcolor < bfacegroup_num; fcolor++)
    {
        if (fcolor == 0)
        {
            fp1 << (*grid).bfacegroup[fcolor] << endl;
        } else{
            fp1 << (*grid).bfacegroup[fcolor] -(*grid).bfacegroup[fcolor-1]<< endl;
        }   
        
    }
    fp1.close();

    fp2.open( "face.txt" , ofstream::app|ofstream::out);
    for (IntType fcolor = 0; fcolor < ifacegroup_num; fcolor++)
    {
        if (fcolor == 0)
        {
            fp2 << (*grid).ifacegroup[fcolor] -(*grid).bfacegroup[bfacegroup_num-1]<< endl;
        } else{
            fp2 << (*grid).ifacegroup[fcolor] -(*grid).ifacegroup[fcolor-1]<< endl;
        }   
    }
    fp2.close();

    //!< deallocation
    sdel_array_1D(index_bface);
    // sdel_array_1D(nNPF_bface);
    // sdel_array_1D(f2n_bface);
    // sdel_array_1D(F2N_bface);
    sdel_array_1D(index_bface_1);

    sdel_array_1D(index_iface);
    // sdel_array_1D(f2c_backup);
    sdel_array_1D(index_iface_1);
}
#endif

/*!
 * @brief       Update boundary-face/interior-face conectivities after reordering.
 * @param       grid
 * @param       index_face
 * @param       flag
 * @remarks     modify according to the fun [void UpdateIface(PolyGrid* grid, IntType* index_iface)] of file[grid_reorder.cpp]
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-03
 */
void UpdateFaceData(PolyGrid *grid, IntType *index_face, IntType flag)
{
    //!< flag == 0 物理边界面
    //!< flag == nBFace 内部面
    IntType *f2c, *f2n, *nNPF;
    f2c = grid->Getf2c();
    f2n = grid->Getf2n();
    nNPF = grid->GetnNPF();

    IntType nBFace, nTFace, nIFace;
    nBFace = grid->GetNBFace();
    nTFace = grid->GetNTFace();
    nIFace = grid->GetNIFace();

    IntType *f2c_backup = NULL;
    IntType *nNPF_backup = NULL;
    IntType *f2n_backup = NULL;
    IntType **F2N_bakcup = NULL;
    IntType facenum, m, nodecount_bface;
    IntType start_I, end_I;

    facenum = nBFace - nIFace; ///< 物理边界面
    start_I = 0;
    end_I = facenum;
    if (flag != 0) ///< 内部面
    {
        facenum = nTFace - flag; ///< nTFace - nBFace;
        start_I = nBFace;
        end_I = nTFace;
    }

    //!< update f2c for boundary/interior faces
    snew_array_1D(f2c_backup, 2 * facenum);
    for (IntType i = start_I; i < end_I; i++)
    {
        f2c_backup[2 * (i - flag)] = f2c[2 * i];
        f2c_backup[2 * (i - flag) + 1] = f2c[2 * i + 1];
    }
    for (IntType i = start_I; i < end_I; i++)
    {
        m = index_face[i - flag] - flag;
        f2c[2 * i] = f2c_backup[m * 2];
        f2c[2 * i + 1] = f2c_backup[m * 2 + 1];
    }

    //!< update nNPF, f2n for boundary/interior faces
    snew_array_1D(nNPF_backup, facenum);
    m = 0;
    for (IntType i = start_I; i < end_I; i++)
    {
        m += nNPF[i];
    }
    snew_array_1D(f2n_backup, m);
    m = 0;
    nodecount_bface = 0;
    if (flag != 0)
    {
        for (IntType j = 0; j < nBFace; ++j) ///< 前面是所有边界面，后面是所有内部面
            nodecount_bface += nNPF[j];
    }
    for (IntType i = start_I; i < end_I; i++)
    {
        nNPF_backup[i - flag] = nNPF[i];
        for (IntType j = 0; j < nNPF[i]; j++)
        {
            f2n_backup[m] = f2n[nodecount_bface + m];
            m++;
        }
    }

    snew_array_1D(F2N_bakcup, facenum);
    F2N_bakcup[0] = f2n_backup;
    for (IntType j = 1; j < facenum; j++)
    {
        F2N_bakcup[j] = &(F2N_bakcup[j - 1][nNPF_backup[j - 1]]);
    }

    //!< update nNPF for boundary/interior faces
    for (IntType i = start_I; i < end_I; i++)
    {
        nNPF[i] = nNPF_backup[index_face[i - flag] - flag];
    }

    //!< update f2n for boundary/interior faces
    m = 0;
    for (IntType i = start_I; i < end_I; i++)
    {
        for (IntType j = 0; j < nNPF[i]; j++)
        {
            f2n[nodecount_bface + m] = F2N_bakcup[index_face[i - flag] - flag][j];
            ++m;
        }
    }
    sdel_array_1D(nNPF_backup);
    sdel_array_1D(f2n_backup);
    sdel_array_1D(F2N_bakcup);
    sdel_array_1D(f2c_backup);
}

#ifdef OMP_GroupColor
/*!
 * @brief       Color contiguous groups of faces
 * @param       grid
 * @param       balanceColors
 * @return      true
 * @return      false
 * @note        参考文献：张健，李瑞田，邓亮，代喆，刘杰，徐传福．面向多核 CPU/众核 GPU 架构的非结构 CFD 共享内存并行计算技术研究[J/OL]．航空学报.
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-03
 */
bool GroupColoring(PolyGrid *grid, bool balanceColors)
{
    IntType MaxColors = grid->MaxColors;
    IntType threads;
    IntType groupEnd, tmp, color, nColor;
    IntType groupSize = grid->groupSize; ///< = 128
    IntType *idxColor = NULL;

    IntType nTCell, nBFace, nTFace, nIFace;
    IntType *f2c, *f2n, *nNPF;
    IntType *index_bface = NULL;
    // IntType *f2c_backup = NULL;
    // IntType *nNPF_bface = NULL;
    // IntType *f2n_bface = NULL;
    // IntType **F2N_bface = NULL;

    IntType *index_iface = NULL;
    /*IntType *nNPF_iface = NULL;
    IntType *f2n_iface = NULL;
    IntType **F2N_iface = NULL;*/

    /*
    IntType* f2c_bface = NULL;
    IntType* index_bface_1 = NULL;
    IntType* f2c_iface = NULL;
    IntType* index_iface_1 = NULL;*/

    nTCell = grid->GetNTCell();
    nBFace = grid->GetNBFace();
    nIFace = grid->GetNIFace();
    nTFace = grid->GetNTFace();
    f2c = grid->Getf2c();
    f2n = grid->Getf2n();
    nNPF = grid->GetnNPF();

    IntType ifacenum = nTFace - nBFace; ///< 内部面
    IntType n = nTCell + nBFace;
    IntType pfacenum = nBFace - nIFace; ///< 物理边界面

#ifdef MF_OPENMP
    threads = omp_get_max_threads();
#else
    threads = 1;
#endif ///< ~MF_OPENMP
    if (pfacenum < groupSize * threads || ifacenum < groupSize * threads)
    {
#ifdef MF_MPICH
        cout << "Rank - " << myid << " skip group coloring." << endl;
#endif
        cout << "groupSize too big." << endl;
        return false;
    }

    //!< allocate memory
    snew_array_1D(idxColor, nTFace);
    snew_array_1D(index_bface, pfacenum);
    // snew_array_1D(f2c_backup, 2 * pfacenum);
    // snew_array_1D(nNPF_bface, pfacenum);
    // snew_array_1D(F2N_bface, pfacenum);
    snew_array_1D(index_iface, ifacenum);

    vector<IntType> colorSize(1, 0);
    vector<IntType> colorSize_bak(MaxColors, 0);
    vector<set<IntType>> innerInColor(MaxColors);
    vector<IntType> searchOrder(MaxColors);

    /* pface color */
    nColor = 1;
    for (IntType i = 0; i < pfacenum; i += groupSize)
    {
        groupEnd = MIN(i + groupSize, pfacenum);
        searchOrder.resize(nColor);
        for (IntType j = 0; j < nColor; j++)
            searchOrder[j] = j;

        //!< Balance sizes by looking for space in smaller colors first. 首先用较小的颜色寻找空间来平衡大小
        if (balanceColors)
        {
            colorSize_bak.resize(nColor);
            colorSize_bak.assign(colorSize.begin(), colorSize.end());
            quicksort_vecint(1, nColor, searchOrder, colorSize_bak);
        }

        vector<IntType>::iterator iter_color = searchOrder.begin();
        for (; iter_color != searchOrder.end(); ++iter_color)
        {
            bool free = true;
            tmp = 2 * i;

            //!< Traverse entire group as a large outer index. 遍历整个组作为一个大的外部索引
            for (IntType j = i; j < groupEnd && free; ++j)
            {
                IntType c = f2c[tmp];
                tmp += 2;
                // if (innerInColor[*iter_color].find(c) != innerInColor[*iter_color].end())
                if (innerInColor[*iter_color].count(c) != 0)
                    free = false;
            }
            //!< If none of the inner indices in the group appears in this color yet, it is assigned to the group.
            //!< 如果组中的内部索引还没有以这种颜色显示，则将其分配给该组
            if (free)
                break;
        }

        IntType color;
        if (iter_color != searchOrder.end())
        {
            //!< Found a color conflict-free. 找到了一个没有冲突的颜色
            color = *iter_color;
        }
        else
        {
            //!< No color was free, make space for a new one. 没有一种颜色是有空位，为新的颜色腾出空间
            color = nColor++;
            if (nColor == MaxColors)
            {
#ifdef MF_MPICH
                cout << "Rank - " << myid << " skip group coloring." << endl;
#endif
                cout << "pface colors over limit." << endl;

                return false;
            }
            colorSize.push_back(0);
        }
        //!< test
        // for (IntType j = i; j < groupEnd; ++j)
        // {
        //     assert(innerInColor[color].count(f2c[2 * j]) == 0);
        // }

        tmp = 2 * i;
        for (IntType j = i; j < groupEnd; ++j)
        {
            idxColor[j] = color;
            innerInColor[color].insert(f2c[tmp]);
            tmp += 2;
        }
        colorSize[color] += groupEnd - i;
    } ///< ~pface color

    tmp = 0;
    for (IntType i = 0; i < nColor; i++)
    {
        tmp += (IntType)colorSize[i];
        grid->bfacegroup.push_back(tmp);
    }
    tmp = 0;
    for (IntType i = 0; i < nColor; i++)
    {
        for (IntType j = 0; j < pfacenum; j++)
        {
            if (idxColor[j] == i)
            {
                index_bface[tmp++] = j;
            }
        }
    }
    // assert(tmp == pfacenum);

    //!< update f2c, nNPF, f2n for physical faces
    UpdateFaceData(grid, index_bface, 0);
    // // update f2c for physical faces
    // for (IntType i = 0; i < pfacenum; i++)
    // {
    //     f2c_backup[2 * i] = f2c[2 * i];
    //     f2c_backup[2 * i + 1] = f2c[2 * i + 1];
    // }
    // for (IntType i = 0; i < pfacenum; i++)
    // {
    //     f2c[2 * i] = f2c_backup[index_bface[i] * 2];
    //     f2c[2 * i + 1] = f2c_backup[index_bface[i] * 2 + 1];
    // }
    // // update nNPF, f2n for physical faces
    // tmp = 0;
    // for (IntType i = 0; i < pfacenum; i++)
    // {
    //     tmp += nNPF[i];
    // }
    // snew_array_1D(f2n_bface, tmp);
    // tmp = 0;
    // for (IntType i = 0; i < pfacenum; i++)
    // {
    //     nNPF_bface[i] = nNPF[i];
    //     for (IntType j = 0; j < nNPF[i]; j++)
    //     {
    //         f2n_bface[tmp] = f2n[tmp];
    //         tmp++;
    //     }
    // }
    // F2N_bface[0] = f2n_bface;
    // for (IntType i = 1; i < pfacenum; i++)
    // {
    //     F2N_bface[i] = &(F2N_bface[i - 1][nNPF_bface[i - 1]]);
    // }
    // for (IntType i = 0; i < pfacenum; i++)
    // {
    //     nNPF[i] = nNPF_bface[index_bface[i]];
    // }
    // tmp = 0;
    // for (IntType i = 0; i < pfacenum; i++)
    // {
    //     for (IntType j = 0; j < nNPF[i]; j++)
    //     {
    //         f2n[tmp] = F2N_bface[index_bface[i]][j];
    //         tmp++;
    //     }
    // }

    // // deallocation
    // sdel_array_1D(index_bface);
    // sdel_array_1D(nNPF_bface);
    // sdel_array_1D(f2n_bface);
    // sdel_array_1D(F2N_bface);

    /* iface color */
    nColor = 1;
    colorSize.resize(1);
    colorSize[0] = 0;
    for (vector<set<IntType>>::iterator iter_color = innerInColor.begin(); iter_color != innerInColor.end(); iter_color++)
        (*iter_color).clear();
    for (IntType i = 0; i < ifacenum; i += groupSize)
    {
        groupEnd = MIN(i + groupSize, ifacenum);
        searchOrder.resize(nColor);

        for (IntType j = 0; j < nColor; j++)
            searchOrder[j] = j;

        //!< Balance sizes by looking for space in smaller colors first.
        if (balanceColors)
        {
            colorSize_bak.resize(nColor);
            colorSize_bak.assign(colorSize.begin(), colorSize.end());
            quicksort_vecint(1, nColor, searchOrder, colorSize_bak);
        }

        vector<IntType>::iterator iter_color = searchOrder.begin();
        for (; iter_color != searchOrder.end(); ++iter_color)
        {
            bool free = true;
            tmp = 2 * (i + nBFace);

            //!< Traverse entire group as a large outer index.
            for (IntType j = i; j < groupEnd && free; ++j)
            {
                IntType c1 = f2c[tmp++];
                IntType c2 = f2c[tmp++];
                // if (innerInColor[*iter_color].find(c) != innerInColor[*iter_color].end())
                if (innerInColor[*iter_color].count(c1) != 0 || innerInColor[*iter_color].count(c2) != 0)
                    free = false;
            }
            //!< If none of the inner indices in the group appears in this color yet, it is assigned to the group.
            if (free)
                break;
        }

        IntType color;
        if (iter_color != searchOrder.end())
        {
            //!< Found a color conflict-free.
            color = *iter_color;
        }
        else
        {
            //!< No color was free, make space for a new one.
            color = nColor++;
            if (nColor == MaxColors)
            {
#ifdef MF_MPICH
                cout << "Rank - " << myid << " skip group coloring." << endl;
#endif
                cout << "iface colors over limit." << endl;
                return false;
            }
            colorSize.push_back(0);
        }
        // test
        // for (IntType j = i; j < groupEnd; ++j)
        // {
        //     assert(innerInColor[color].count(f2c[2 * (j + nBFace)]) == 0);
        //     assert(innerInColor[color].count(f2c[2 * (j + nBFace) + 1]) == 0);
        // }

        tmp = 2 * (i + nBFace);
        for (IntType j = i; j < groupEnd; ++j)
        {
            idxColor[j + nBFace] = color;
            innerInColor[color].insert(f2c[tmp++]);
            innerInColor[color].insert(f2c[tmp++]);
        }
        colorSize[color] += groupEnd - i;
    } ///< ~iface color

    tmp = nBFace;
    for (IntType i = 0; i < nColor; i++)
    {
        tmp += (IntType)colorSize[i];
        grid->ifacegroup.push_back(tmp);
    }
    tmp = 0;
    for (IntType i = 0; i < nColor; i++)
    {
        for (IntType j = nBFace; j < nTFace; j++)
        {
            if (idxColor[j] == i)
                index_iface[tmp++] = j;
        }
    }
    // assert(tmp == ifacenum);
    for (vector<set<IntType>>::iterator iter_color = innerInColor.begin(); iter_color != innerInColor.end(); iter_color++)
        (*iter_color).clear();
    sdel_array_1D(idxColor);

    //!< update f2c nNPF, f2n for interior faces
    // UpdateIface(grid, index_iface);
    UpdateFaceData(grid, index_iface, nBFace);

    //!< deallocation
    sdel_array_1D(index_iface);
    // sdel_array_1D(f2c_backup);
    return true;
}

/*!
 * @brief       To sort a list of points for IntType stored in vector<IntType>
 * @param       s
 * @param       e
 * @param       order
 * @param       vecint
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-03
 */
void quicksort_vecint(IntType s, IntType e, vector<IntType> &order, vector<IntType> &vecint)
{
    IntType i, last;

    if (e - s <= 0) // nothing to do
        return;

    last = s - 1;
    for (i = s; i < e; i++)
    {
        if (vecint[i] < vecint[s - 1]) // ascending
        {
            // xswap_vecint(order, vecint, ++last, i);
            ++last;
            std::swap(vecint[last], vecint[i]);
            std::swap(order[last], order[i]);
        }
    }

    // xswap_vecint(order, vecint, s - 1, last);
    std::swap(vecint[s - 1], vecint[last]);
    std::swap(order[s - 1], order[last]);
    quicksort_vecint(s, last, order, vecint);
    quicksort_vecint(last + 2, e, order, vecint);
}
#endif ///< ~OMP_GroupColor

#ifdef OMP_DIVREP
/*!
 * @brief       Divide subzone by Metis, and replicate face
 * @param       grid
 * @return      true
 * @return      false
 * @note        参考文献：张健，李瑞田，邓亮，代喆，刘杰，徐传福．面向多核 CPU/众核 GPU 架构的非结构 CFD 共享内存并行计算技术研究[J/OL]．航空学报.
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @date        2023-07-03
 */
bool DivRep(PolyGrid *grid)
{
    IntType i, j;
    IntType nTCell = grid->GetNTCell();
    IntType nBFace = grid->GetNBFace();
    IntType nTFace = grid->GetNTFace();
    IntType nSubZone = grid->threads;
    IntType *nCPC_tmp = NULL;
    IntType **c2c_tmp = NULL;
    IntType *f2c;
    IntType count, c1, c2, subzone;

    idx_t *xadj = NULL;   ///< csr, start and end index
    idx_t *adjncy = NULL; ///< csr, cell id
    idx_t nvtxs, ncon, nparts, objval;
    idx_t *part = NULL;
    IntType status;

    // IntType *idx;
    IntType *nface_pthread = NULL;
    if (nSubZone < 2)
    {
        // mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
        cerr << "Error in func[DivRep()] of file[parallel_omp.cpp]!" << endl;
        exit(1);
    }
    /* Division */
    {
        //!< get cell to cell, exclude ghost
        snew_array_1D(nCPC_tmp, nTCell);
        for (i = 0; i < nTCell; i++)
            nCPC_tmp[i] = 0;
        f2c = grid->Getf2c();
        count = 2 * nBFace;
        for (i = nBFace; i < nTFace; i++)
        {
            c1 = f2c[count++];
            c2 = f2c[count++];
            nCPC_tmp[c1]++;
            nCPC_tmp[c2]++;
        }
        snew_array_2D(c2c_tmp, nTCell, nCPC_tmp, true);
        for (i = 0; i < nTCell; i++)
            nCPC_tmp[i] = 0;
        count = 2 * nBFace;
        for (i = nBFace; i < nTFace; i++)
        {
            c1 = f2c[count++];
            c2 = f2c[count++];
            c2c_tmp[c1][nCPC_tmp[c1]++] = c2;
            c2c_tmp[c2][nCPC_tmp[c2]++] = c1;
        }

        //!< csr
        snew_array_1D(xadj, nTCell + 1);
        snew_array_1D(adjncy, 2 * (nTFace - nBFace));
        xadj[0] = 0;
        count = 0;
        for (i = 0; i < nTCell; i++)
        {
            xadj[i + 1] = xadj[i] + nCPC_tmp[i];
            for (j = 0; j < nCPC_tmp[i]; j++)
            {
                adjncy[count++] = c2c_tmp[i][j];
            }
        }
        sdel_array_1D(nCPC_tmp);
        sdel_array_2D(c2c_tmp);

        //!< divide by Metis
        nvtxs = (idx_t)nTCell;
        ncon = 1;
        nparts = (idx_t)nSubZone;
        snew_array_1D(part, static_cast<size_t>(nvtxs));
        for (i = 0; i < nvtxs; i++)
        {
            part[i] = 0;
        }
        status = METIS_PartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, NULL, NULL, NULL,
                                          &nparts, NULL, NULL, NULL, &objval, part);
        sdel_array_1D(xadj);
        sdel_array_1D(adjncy);
        xadj = NULL;
        adjncy = NULL;
        if (status != METIS_OK)
        {
            sdel_array_1D(part);
            // mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
            cerr << "Error in func[DivRep()] of file[parallel_omp.cpp]!" << endl;
            exit(1);
        }
    }

    /* Replication */
    {
        //!< Bface
        // grid->idx_pthreads_bface = NULL; ///< 前 n 个线程 Bface 的总数（第一个是 0，第二个是 tid0 的 Bface 个数，第三个是 tid0 和 tid1 的 Bface 个数之和，以此类推）
        // grid->id_division_bface = NULL;  ///< 每个线程 Bface 的编号（按照线程顺序排列）
        // nface_pthread = NULL; ///< 每个线程对应 face(Bface/Iface) 的个数
        snew_array_1D(grid->idx_pthreads_bface, nSubZone + 1);
        snew_array_1D(nface_pthread, nSubZone);
        for (i = 0; i < nSubZone + 1; i++)
            grid->idx_pthreads_bface[i] = 0;
        count = 0;
        for (i = 0; i < nBFace; i++)
        {
            c1 = f2c[count];
            count += 2;
            subzone = (IntType)part[c1];
            grid->idx_pthreads_bface[subzone + 1]++;
        }
        for (i = 1; i < nSubZone + 1; i++)
        {
            grid->idx_pthreads_bface[i] += grid->idx_pthreads_bface[i - 1];
            nface_pthread[i - 1] = 0;
        }
        snew_array_1D(grid->id_division_bface, grid->idx_pthreads_bface[nSubZone]);
        count = 0;
        for (i = 0; i < nBFace; i++)
        {
            c1 = f2c[count];
            count += 2;
            subzone = (IntType)part[c1];
            grid->id_division_bface[grid->idx_pthreads_bface[subzone] + nface_pthread[subzone]] = i;
            nface_pthread[subzone]++;
        }

        //!< Iface
        // grid->idx_pthreads_iface = NULL;
        // grid->id_division_iface = NULL;
        snew_array_1D(grid->idx_pthreads_iface, nSubZone + 1);
        for (i = 0; i < nSubZone + 1; i++)
            grid->idx_pthreads_iface[i] = 0;
        count = 2 * nBFace;
        for (i = nBFace; i < nTFace; i++)
        {
            c1 = f2c[count++];
            c2 = f2c[count++];
            subzone = (IntType)part[c1];
            grid->idx_pthreads_iface[subzone + 1]++;
            if (part[c1] != part[c2])
            {
                //!< subzone boundary face, replication
                subzone = (IntType)part[c2];
                grid->idx_pthreads_iface[subzone + 1]++;
            }
        }
        for (i = 1; i < nSubZone + 1; i++)
        {
            grid->idx_pthreads_iface[i] += grid->idx_pthreads_iface[i - 1];
            nface_pthread[i - 1] = 0;
        }
        // assert(grid->idx_pthreads_iface[nSubZone] - nTFace + nBFace == objval);
        snew_array_1D(grid->id_division_iface, grid->idx_pthreads_iface[nSubZone]);
        count = 2 * nBFace;
        for (i = nBFace; i < nTFace; i++)
        {
            c1 = f2c[count++];
            c2 = f2c[count++];
            if (part[c1] != part[c2])
            {
                //!< subzone boundary face, replication
                subzone = (IntType)part[c1];
                j = i + nTFace;
                grid->id_division_iface[grid->idx_pthreads_iface[subzone] + nface_pthread[subzone]] = j;
                nface_pthread[subzone]++;
                subzone = (IntType)part[c2];
                j = -1 * (i + nTFace);
                grid->id_division_iface[grid->idx_pthreads_iface[subzone] + nface_pthread[subzone]] = j;
                nface_pthread[subzone]++;
            }
            else
            {
                subzone = (IntType)part[c1];
                grid->id_division_iface[grid->idx_pthreads_iface[subzone] + nface_pthread[subzone]] = i;
                nface_pthread[subzone]++;
            }
        }
        sdel_array_1D(nface_pthread);
    }

    sdel_array_1D(part);
    return true;
}

// #ifdef BoundedColoring
// /************************************************************************
//                Color faces of each subzone after Divide & Replicate
//                         Add by dingxin 2021-11-30
// ************************************************************************/
// void colorAfterDivRep(PolyGrid *grid)
// {
//     IntType i, t, k, startFace, endFace;
//     IntType n_bface, n_iface;
//     IntType nTFace = grid->GetNTFace();
//     IntType *f2c = grid->Getf2c();
//     IntType bfacenum_vec, ifacenum_vec;

// #pragma omp parallel for private(t, i, k, startFace, endFace, n_bface, n_iface, bfacenum_vec, ifacenum_vec)
//     for (t = 0; t < grid->threads; t++)
//     {
//         IntType *bface = NULL;
//         IntType *iface = NULL;
//         IntType *order_b = NULL;
//         IntType *order_i = NULL;
//         // Boundary faces
//         startFace = grid->idx_pthreads_bface[t];
//         grid->endIndex_bFace_vec[t] = startFace;
//         endFace = grid->idx_pthreads_bface[t + 1];
//         n_bface = endFace - startFace;
//         snew_array_1D(bface, n_bface);
//         snew_array_1D(order_b, n_bface);
//         for (i = startFace; i < endFace; i++)
//             bface[i - startFace] = grid->id_division_bface[i];

//         // Interior faces
//         startFace = grid->idx_pthreads_iface[t];
//         grid->endIndex_iFace_vec[t] = startFace;
//         endFace = grid->idx_pthreads_iface[t + 1];
//         n_iface = endFace - startFace;
//         snew_array_1D(iface, n_iface);
//         snew_array_1D(order_i, n_iface);
//         for (i = startFace; i < endFace; i++)
//         {
//             k = grid->id_division_iface[i];
//             if (abs(k) < nTFace)
//                 iface[i - startFace] = k;
//             else
//                 iface[i - startFace] = abs(k) - nTFace;
//         }
//         BoundedColoring(f2c, bface, iface, order_b, order_i, n_bface, n_iface, bfacenum_vec, ifacenum_vec);
//         grid->endIndex_bFace_vec[t] += bfacenum_vec;
//         grid->endIndex_iFace_vec[t] += ifacenum_vec;
//         updateOrder(grid->id_division_bface, order_b, grid->idx_pthreads_bface[t], n_bface);
//         updateOrder(grid->id_division_iface, order_i, startFace, n_iface);
// #ifdef MF_DEBUG
//         {
//             updateOrder(bface, order_b, 0, n_bface);
//             if (!validation_decoupling(bface, 0, bfacenum_vec, f2c, BOUNDARYFACE))
//             {
//                 mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
//             }
//             updateOrder(iface, order_i, 0, n_iface);
//             if (!validation_decoupling(iface, 0, ifacenum_vec, f2c, INTERIORFACE))
//             {
//                 mflow_exit(mflow_error_flag(CPP_FILD_ID, CPP_LINE));
//             }
//         }
// #endif // MF_DEBUG
//         sdel_array_1D(bface);
//         sdel_array_1D(iface);
//         sdel_array_1D(order_b);
//         sdel_array_1D(order_i);
//     }
// }
// #endif ///< ~BoundedColoring
#endif ///< ~OMP_DIVREP
