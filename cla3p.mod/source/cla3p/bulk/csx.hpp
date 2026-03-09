/*
 * Copyright (c) 2025-2026 Simulisoft
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef CLA3P_BULK_CSX_HPP_
#define CLA3P_BULK_CSX_HPP_

/**
 * @file
 */

#include <iostream>

#include "cla3p/types/integer.hpp"
#include "cla3p/types/enums.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace blk {
namespace csx {
/*-------------------------------------------------*/

enum class csx_t {
    CSR   = 0,
    CSC       
};

void roll(int_t np, int_t *xxxptr);

void unroll(int_t np, int_t *xxxptr);

int_t maxXlen(int_t np, const int_t *xxxptr);

void check(prop_t ptype, uplo_t uplo, csx_t type, int_t np, int_t maxIdx, const int_t *xxxptr, const int_t *xxxidx);

void sort(int_t np, const int_t *xxxptr, int_t *xxxidx);

template <typename T_Scalar>
void sort(int_t np, const int_t *xxxptr, int_t *xxxidx, T_Scalar *values);

template <typename T_Scalar>
void print_to_stream(std::ostream& os, csx_t type, int_t np, const int_t *xxxptr, const int_t *xxxidx, const T_Scalar *values, std::streamsize prec = 0);

template <typename T_Scalar>
void transpose(int_t inp, const int_t *ixptr, const int_t *ixidx, const T_Scalar *ivals,
               int_t onp, int_t *oxptr, int_t *oxidx, T_Scalar *ovals, T_Scalar coeff = 1);

template <typename T_Scalar>
void conjugate_transpose(int_t inp, const int_t *ixptr, const int_t *ixidx, const T_Scalar *ivals,
                         int_t onp, int_t *oxptr, int_t *oxidx, T_Scalar *ovals, T_Scalar coeff = 1);

void uplo2ge_xxptr(uplo_t uplo, csx_t type, int_t np, const int_t *ixptr, const int_t *ixidx, int_t *oxptr);

template <typename T_Scalar>
void sy2ge(uplo_t uplo, csx_t type, int_t np, 
           const int_t *ixptr, const int_t *ixidx, const T_Scalar *ivals,
           int_t *oxptr, int_t *oxidx, T_Scalar *ovals);

template <typename T_Scalar>
void he2ge(uplo_t uplo, csx_t type, int_t np, 
           const int_t *ixptr, const int_t *ixidx, const T_Scalar *ivals,
           int_t *oxptr, int_t *oxidx, T_Scalar *ovals);

template <typename T_Scalar>
void remove_duplicates(int_t np, int_t *xxxptr, int_t *xxxidx, T_Scalar *values, dup_t op);

template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type 
idx_max_asum(prop_t ptype, int_t np, const int_t *xxxptr, const int_t *xxxidx, const T_Scalar *values);

template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type 
ptr_max_asum(prop_t ptype, int_t np, int_t maxIdx, const int_t *xxxptr, const int_t *xxxidx, const T_Scalar *values);

template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type 
norm_max(int_t np, const int_t *xxxptr, const T_Scalar *values);

template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type 
norm_fro(prop_t ptype, int_t np, const int_t *xxxptr, const int_t *xxxidx, const T_Scalar *values);

/*-------------------------------------------------*/
} // namespace csx
} // namespace blk
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_BULK_CSX_HPP_