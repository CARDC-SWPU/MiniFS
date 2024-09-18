/*!
 * @file        memory_util.h
 * @brief       A library for dynamic memory management in C++
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @version     1.0
 * @date        2023-06-28
 * @copyright   Copyright (c) 2023  PERSONAL
 *
 * @par 修改日志:
 * <table>
 * <tr><th> Date        <th> Version  <th> Author  <th> Description
 * <tr><td> 2023-06-28  <td> 1.0      <td> Wisces  <td> 将动态内存管理单独提取出来
 * </table>
 */

#ifndef MF_MEMORY_UTIL_H
#define MF_MEMORY_UTIL_H

//!< C++ build-in head files
#include <iostream>

//!< user defined head files
#include "number_type.h"

template <typename T>
inline void sdel_object(T *&ptr)
{
    // std::cout << "The func [sdel_object] is called!" << std::endl;

    if (ptr != NULL)
    {
        delete ptr;
        ptr = NULL;
    }
}

template <typename T>
inline void sdel_array_1D(T *&ptr)
{
    // std::cout << "The func [sdel_array_1D] is called!" << std::endl;

    if (ptr != NULL)
    {
        delete[] ptr;
        ptr = NULL;
    }
}

template <typename T>
inline void sdel_array_2D(T **&ptr, const std::size_t n = 1, const bool contiguous = true)
{
    // std::cout << "The func [sdel_array_2D] is called!" << std::endl;

    if (ptr != NULL)
    {
        // first delete the dynamic memory of the second dimension
        if (contiguous)
        {
            // if(ptr[0] != NULL)
            sdel_array_1D(ptr[0]);
        }
        else
        {
            for (std::size_t i = 0; i < n; ++i)
            {
                // if(ptr[i] != NULL)
                sdel_array_1D(ptr[i]);
            }
        }

        // delete the 'ptr' and erase it from the map
        sdel_array_1D(ptr);
    }
}

template <typename T>
inline void snew_array_1D(T *&ptr, const std::size_t n)
{
    // std::cout << "The func [snew_array_1D] is called!" << std::endl;

    sdel_array_1D<T>(ptr);
    ptr = new T[n];
}

template <typename T>
inline void snew_array_2D(T **&ptr, const std::size_t n, const std::size_t m, const bool contiguous = true)
{
    // std::cout << "The func [snew_array_2D] is called!" << std::endl;

    sdel_array_2D(ptr, n, contiguous);

    snew_array_1D(ptr, n);
    if (contiguous)
    {
        ptr[0] = NULL;
        snew_array_1D(ptr[0], n * m);
        for (std::size_t i = 1; i < n; ++i)
        {
            ptr[i] = &(ptr[i - 1][m]);
        }
    }
    else
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            ptr[i] = NULL;
            snew_array_1D(ptr[i], m);
        }
    }
}

template <typename T>
inline void snew_array_2D(T **&ptr, const std::size_t n, const IntType *size, const bool contiguous = true)
{
    // std::cout << "The func [snew_array_2D] is called!" << std::endl;

    sdel_array_2D(ptr, n, contiguous);

    snew_array_1D(ptr, n);

    if (contiguous)
    {
        std::size_t sum = 0;
        for (std::size_t i = 0; i < n; ++i)
            sum += size[i];

        ptr[0] = NULL;
        snew_array_1D(ptr[0], sum);

        for (std::size_t i = 1; i < n; ++i)
        {
            ptr[i] = &(ptr[i - 1][size[i - 1]]);
        }
    }
    else
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            ptr[i] = NULL;
            snew_array_1D(ptr[i], size[i]);
        }
    }
}

template <typename T>
inline void sSet(T *&ptr_old, T *ptr_new)
{
    if ((ptr_old != NULL) && (ptr_old != ptr_new))
    {
        sdel_array_1D(ptr_old);
    }
    ptr_old = ptr_new;
}

template <typename T>
inline void sSet(T **&ptr_old, T **ptr_new)
{
    if ((ptr_old != NULL) && (ptr_old != ptr_new))
    {
        sdel_array_2D(ptr_old); // contiguous mode
    }
    ptr_old = ptr_new;
}

template <typename T>
inline void sSet(T **&ptr_old, T **ptr_new, std::size_t n)
{
    if ((ptr_old != NULL) && (ptr_old != ptr_new))
    {
        sdel_array_2D(ptr_old, n, false); // non-contiguous mode
    }
    ptr_old = ptr_new;
}

#endif ///< ~MF_MEMORY_UTIL_H