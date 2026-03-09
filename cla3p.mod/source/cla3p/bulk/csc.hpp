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

#ifndef CLA3P_BULK_CSC_HPP_
#define CLA3P_BULK_CSC_HPP_

/**
 * @file
 */

#include <iostream>

#include "cla3p/bulk/csx.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace blk {
namespace csc {
/*-------------------------------------------------*/

inline void check(prop_t ptype, uplo_t uplo, int_t m, int_t n, const int_t *colptr, const int_t *rowidx)
{
    csx::check(ptype, uplo, csx::csx_t::CSC, n, m, colptr, rowidx);
}

template <typename T_Scalar>
void print_to_stream(std::ostream& os, int_t n, const int_t *colptr, const int_t *rowidx, const T_Scalar *values, std::streamsize prec = 0)
{
    csx::print_to_stream(os, csx::csx_t::CSC, n, colptr, rowidx, values, prec);
}

template <typename T_Scalar>
void print(int_t n, const int_t *colptr, const int_t *rowidx, const T_Scalar *values, std::streamsize prec = 0)
{
    print_to_stream(std::cout, n, colptr, rowidx, values, prec);
}

template <typename T_Scalar>
void transpose(int_t m, int_t n, const int_t *icolptr, const int_t *irowidx, const T_Scalar *ivalues,
               int_t *ocolptr, int_t *orowidx, T_Scalar *ovalues, T_Scalar coeff = 1)
{
    csx::transpose(n, icolptr, irowidx, ivalues,
                   m, ocolptr, orowidx, ovalues, coeff);
}

template <typename T_Scalar>
void conjugate_transpose(int_t m, int_t n, const int_t *icolptr, const int_t *irowidx, const T_Scalar *ivalues,
                         int_t *ocolptr, int_t *orowidx, T_Scalar *ovalues, T_Scalar coeff = 1)
{
    csx::conjugate_transpose(n, icolptr, irowidx, ivalues,
                             m, ocolptr, orowidx, ovalues, coeff);
}

inline void uplo2ge_colptr(uplo_t uplo, int_t n, const int_t *icolptr, const int_t *irowidx, int_t *ocolptr)
{
    csx::uplo2ge_xxptr(uplo, csx::csx_t::CSC, n, icolptr, irowidx, ocolptr);
}

template <typename T_Scalar>
void sy2ge(uplo_t uplo, int_t n, 
           const int_t *icolptr, const int_t *irowidx, const T_Scalar *ivalues,
           int_t *ocolptr, int_t *orowidx, T_Scalar *ovalues)
{
    csx::sy2ge(uplo, csx::csx_t::CSC, n, 
               icolptr, irowidx, ivalues,
               ocolptr, orowidx, ovalues);
}

template <typename T_Scalar>
void he2ge(uplo_t uplo, int_t n, 
           const int_t *icolptr, const int_t *irowidx, const T_Scalar *ivalues,
           int_t *ocolptr, int_t *orowidx, T_Scalar *ovalues)
{
    csx::he2ge(uplo, csx::csx_t::CSC, n, 
               icolptr, irowidx, ivalues,
               ocolptr, orowidx, ovalues);
}

template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type 
norm_one(prop_t ptype, int_t n, const int_t *colptr, const int_t *rowidx, const T_Scalar *values)
{
    return csx::idx_max_asum(ptype, n, colptr, rowidx, values);
}

template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type 
norm_inf(prop_t ptype, int_t m, int_t n, const int_t *colptr, const int_t *rowidx, const T_Scalar *values)
{
    return csx::ptr_max_asum(ptype, n, m, colptr, rowidx, values);
}

template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type 
norm_max(int_t n, const int_t *colptr, const T_Scalar *values)
{
    return csx::norm_max(n, colptr, values);
}

template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type 
norm_fro(prop_t ptype, int_t n, const int_t *colptr, const int_t *rowidx, const T_Scalar *values)
{
    return csx::norm_fro(ptype, n, colptr, rowidx, values);
}

template <typename T_Scalar>
void permute(prop_t ptype, uplo_t uplo, int_t m, int_t n, 
             const int_t *icolptr, const int_t *irowidx, const T_Scalar *ivalues,
             int_t *ocolptr, int_t *orowidx, T_Scalar *ovalues, const int_t *P, const int_t *Q);

/*-------------------------------------------------*/
} // namespace csc
} // namespace blk
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_BULK_CSC_HPP_
