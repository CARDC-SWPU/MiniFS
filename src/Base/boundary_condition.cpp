/*!
 * @file        boundary_condition.cpp
 * @brief       BCond Object
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @version     1.0
 * @date
 * @copyright   Copyright (c) 2023  PERSONAL
 *
 * @par 修改日志:
 * <table>
 * <tr><th> Date        <th> Version  <th> Author  <th> Description
 * <tr><td> 2023-06-15  <td> 1.0      <td> Wisces  <td> 〔内容〕
 * <tr><td> 2023-06-20  <td> 2.0      <td> Wisces  <td> 添加了构造函数 BCRecord(const IntType patch_id_, const IntType type, const ShortString type_symbols);
 * </table>
 */

// direct head file
#include "boundary_condition.h"

// build-in head files
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
// #include <cassert>
#include <iostream>
using namespace std;

// other user defined head file
#include "number_type.h"
#include "grid_base.h"
#include "memory_util.h"
// #include "memory_util.h"

BCond::~BCond()
{
    IntType i;
    for (i = 0; i < nBCRecord; i++)
    {
        // if (bcRecord[i] != NULL)
        // {
        //     delete bcRecord[i];
        //     bcRecord[i] = NULL;
        // }
        sdel_object(bcRecord[i]);
    }
    if (bcRecord != NULL)
    {
        delete[] bcRecord;
        bcRecord = NULL;
    }
    // sdel_array_1D(bcRecord);
}

void BCond::AddBCRecord(BCRecord *bcRecordin)
{
    if (nBCRecord == 0)
    {
        snew_array_1D(bcRecord, 1);
        bcRecord[nBCRecord++] = bcRecordin;
    }
    else
    {
        IntType i;
        BCRecord **bcRecordt = bcRecord;
        bcRecord = NULL;
        snew_array_1D(bcRecord, nBCRecord + 1);
        for (i = 0; i < nBCRecord; i++)
            bcRecord[i] = bcRecordt[i];
        bcRecord[nBCRecord++] = bcRecordin;
        sdel_array_1D(bcRecordt);
    }
}

// Copy boundary conditions from the other
// void BCond::CopyFrom(const BCond *other)
// {
//     IntType nbc = other->GetNoBCR();
//     for (IntType ibc = 0; ibc < nbc; ++ibc)
//     {
//         BCRecord *bcr_other = other->GetBCRecord(ibc);
//         BCRecord *bcr_this = this->FindBCRecord(bcr_other->GetPatchID());

//         if (bcr_this != NULL)
//         {
//             // bcr_this->CopyFrom(bcr_other);
//         }
//         else
//         {
//             BCRecord *bcr = new BCRecord();
//             // snew_object(bcr);
//             // bcr->CopyFrom(other->GetBCRecord(ibc));
//             this->AddBCRecord(bcr);
//         }
//     }
// }

/// Return the BCRecord which patch id is patch_id, return NULL if not exist.
// BCRecord *BCond::FindBCRecord(const IntType patch_id)
// {
//     for (IntType ibc = 0; ibc < this->GetNoBCR(); ++ibc)
//     {
//         BCRecord *current_bcr = this->GetBCRecord(ibc);
//         if (current_bcr->GetPatchID() == patch_id)
//         {
//             return current_bcr;
//         }
//     }

//     return NULL;
// }

// void BCRecord::UpdateData(void *data, IntType typein, IntType size, const ShortString name)
// {
//     bcData_.UpdateDataSafe(data, typein, size, name);
// }
// void BCRecord::GetBCVar(void *data, IntType type, const ShortString name, IntType messageOn)
// {
//     bcData_.GetDataByName(data, type, 1, name, messageOn);
// }

// void BCRecord::GetBCVar(void *data, IntType type, const ShortString name)
// {
//     bcData_.GetDataByName(data, type, 1, name);
// }

// // delete extra data
// void BCRecord::EraseExtraData(void)
// {
//     bcData_.DeleteAllData();
// }

// /// \brief Deeply copy data from the other
// void BCRecord::CopyFrom(const BCRecord *other)
// {
//     this->SetType(other->GetType());
//     this->SetPatchID(other->GetPatchID());
//     this->SetTypeSymbol(*(other->GetTypeSymbol()));
//     this->SetPatchName(*(other->GetPatchName()));
//     this->EraseExtraData();
//     this->bcData_.CopyDataFrom(other->GetDataPointer());
// }
