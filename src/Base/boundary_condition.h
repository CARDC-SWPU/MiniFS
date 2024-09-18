/*!
 * @file        boundary_condition.h
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

#ifndef MF_BOUNDARY_CONDITION_H
#define MF_BOUNDARY_CONDITION_H

#include <string.h>
#include "number_type.h"
#include "constant.h"

class BCRecord
{

private:
    IntType type_, patch_id_; // the bc key(or bc type), patch id (base 1)
    // ShortString type_symbols_; // the type symbol, size less than MAX_SHORT_STRING
    // String patch_name_;        // the type string, patch name

public:
    /// \brief Set the patch type string defined in MFlow,
    /// such as "wall", "symm", "far_field", et al.
    // void SetTypeSymbol(const ShortString type_symbols);

    // /// \brief Set the name of this boundary patch
    // void SetPatchName(const String patch_name);

    /// \brief Set the patch type of an integer
    void SetType(IntType type);

    /// \brief Set the id of this boundary patch
    void SetPatchID(IntType id);

    /// \brief Get the bc type string defined in MFlow,
    /// such as "wall", "symm", "far_field", et al.
    // const ShortString *GetTypeSymbol() const;

    // /// \brief Get the patch name
    // const String *GetPatchName() const;

    /// \brief Get the patch type of an integer
    IntType GetType() const;

    /// \brief Get the id of this boundary patch
    IntType GetPatchID() const;

    // /// \brief delete all data which are given with $ for this BCRecord.
    // void EraseExtraData(void);

    // /// \brief Deeply copy data from the other
    // void CopyFrom(const BCRecord *other);

    // constructor
    BCRecord(const IntType patch_id_, const IntType type);
    // BCRecord(const IntType type, const ShortString type_symbols, const String patch_name);
    // BCRecord(const ShortString type_symbols, const String patch_name);
    // explicit BCRecord(const ShortString type_symbols);
    // explicit BCRecord(const IntType type);
    BCRecord();

    // destructor
    ~BCRecord();
};

// inline functions for class BCRecord

// // Set the bc type with an string
// inline void BCRecord::SetTypeSymbol(const ShortString type_symbols)
// {
//     strcpy(type_symbols_, type_symbols);
// }

// // Set the name of this boundary patch
// inline void BCRecord::SetPatchName(const String patch_name)
// {
//     strcpy(patch_name_, patch_name);
// }

// Set the patch type with an integer
inline void BCRecord::SetType(IntType type)
{
    type_ = type;
};

// Set the id of this boundary patch
inline void BCRecord::SetPatchID(IntType id)
{
    patch_id_ = id;
}

// // Get the bc type with an string defined in MFlow,
// // such as "wall", "symm", "far_field", et al.
// inline const ShortString *BCRecord::GetTypeSymbol() const
// {
//     return &type_symbols_;
// }

// // Get the patch name
// inline const String *BCRecord::GetPatchName() const
// {
//     return &patch_name_;
// }

// Get the patch type with an integer
inline IntType BCRecord::GetType() const
{
    return type_;
}

// Get the id of this boundary patch
inline IntType BCRecord::GetPatchID() const
{
    return patch_id_;
}

inline BCRecord::BCRecord(const IntType patch_id, const IntType type) : patch_id_(patch_id), type_(type)
{
}

// inline BCRecord::BCRecord(const IntType type, const ShortString type_symbols, const String patch_name) : type_(type)
// {
//     strcpy(type_symbols_, type_symbols);
//     strcpy(patch_name_, patch_name);
// }

// inline BCRecord::BCRecord(const ShortString type_symbols, const String patch_name) : type_(0)
// {
//     strcpy(type_symbols_, type_symbols);
//     strcpy(patch_name_, patch_name);
// }

// inline BCRecord::BCRecord(const ShortString type_symbols) : type_(0)
// {
//     strcpy(type_symbols_, type_symbols);
//     patch_name_[0] = '\0';
// }

// inline BCRecord::BCRecord(const IntType type) : type_(type)
// {
//     type_symbols_[0] = '\0';
//     patch_name_[0] = '\0';
// }

inline BCRecord::BCRecord() : type_(0)
{
    // type_symbols_[0] = '\0';
    // patch_name_[0] = '\0';
}

inline BCRecord::~BCRecord()
{
}

//
// Container to save several BCRecord
class BCond
{
    IntType nBCRecord;
    BCRecord **bcRecord; // the bc rec. each boun. face is assoc.

public:
    /// \brief constructor
    BCond();

    /// \brief Get the number of BCRecords
    IntType GetNoBCR() const;

    /// \brief Add a BCRecord
    void AddBCRecord(BCRecord *bcRecord);

    /// \brief Get the i-th BCRecord item
    BCRecord *GetBCRecord(IntType i) const;

    // /// \brief Add or update boundary conditions from the other
    // void CopyFrom(const BCond *other);

    // /// \brief Return the BCRecord which patch id is patch_id, return NULL if not exist.
    // BCRecord *FindBCRecord(const IntType patch_id);

    /// \brief Destructor
    ~BCond();
};

// inline functions for class BCond

inline BCond::BCond() : nBCRecord(0)
{
    bcRecord = NULL;
}

inline IntType BCond::GetNoBCR() const
{
    return nBCRecord;
}

inline BCRecord *BCond::GetBCRecord(IntType i) const
{
    return bcRecord[i];
}

#endif //~MF_BOUNDARY_CONDITION_H
