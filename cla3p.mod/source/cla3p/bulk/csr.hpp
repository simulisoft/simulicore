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

#ifndef CLA3P_BULK_CSR_HPP_
#define CLA3P_BULK_CSR_HPP_

/**
 * @file
 */

#include <iostream>

#include "cla3p/bulk/csx.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace blk {
namespace csr {
/*-------------------------------------------------*/

inline void check(prop_t ptype, uplo_t uplo, int_t m, int_t n, const int_t *rowptr, const int_t *colidx)
{
    csx::check(ptype, uplo, csx::csx_t::CSR, m, n, rowptr, colidx);
}

template <typename T_Scalar>
void print_to_stream(std::ostream& os, int_t m, const int_t *rowptr, const int_t *colidx, const T_Scalar *values, std::streamsize prec = 0)
{
    csx::print_to_stream(os, csx::csx_t::CSR, m, rowptr, colidx, values, prec);
}

template <typename T_Scalar>
void print(int_t m, const int_t *rowptr, const int_t *colidx, const T_Scalar *values, std::streamsize prec = 0)
{
	print_to_stream(std::cout, m, rowptr, colidx, values, prec);
}

template <typename T_Scalar>
void transpose(int_t m, int_t n, 
               const int_t *irowptr, const int_t *icolidx, const T_Scalar *ivalues,
		       int_t *orowptr, int_t *ocolidx, T_Scalar *ovalues, T_Scalar coeff = 1)
{
    csx::transpose(m, irowptr, icolidx, ivalues,
                   n, orowptr, ocolidx, ovalues, coeff);
}

template <typename T_Scalar>
void conjugate_transpose(int_t m, int_t n, 
                         const int_t *irowptr, const int_t *icolidx, const T_Scalar *ivalues,
		                 int_t *orowptr, int_t *ocolidx, T_Scalar *ovalues, T_Scalar coeff = 1)
{
    csx::conjugate_transpose(m, irowptr, icolidx, ivalues,
                             n, orowptr, ocolidx, ovalues, coeff);
}

inline void uplo2ge_rowptr(uplo_t uplo, int_t m, const int_t *irowptr, const int_t *colxidx, int_t *orowptr)
{
    csx::uplo2ge_xxptr(uplo, csx::csx_t::CSR, m, irowptr, colxidx, orowptr);
}

template <typename T_Scalar>
void sy2ge(uplo_t uplo, int_t m, 
           const int_t *irowptr, const int_t *icolidx, const T_Scalar *ivalues,
		   int_t *orowptr, int_t *ocolidx, T_Scalar *ovalues)
{
    csx::sy2ge(uplo, csx::csx_t::CSR, m, 
               irowptr, icolidx, ivalues,
               orowptr, ocolidx, ovalues);
}

template <typename T_Scalar>
void he2ge(uplo_t uplo, int_t m, 
           const int_t *irowptr, const int_t *icolidx, const T_Scalar *ivalues,
		   int_t *orowptr, int_t *ocolidx, T_Scalar *ovalues)
{
    csx::he2ge(uplo, csx::csx_t::CSR, m, 
               irowptr, icolidx, ivalues,
               orowptr, ocolidx, ovalues);
}

template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type 
norm_one(prop_t ptype, int_t m, int_t n, const int_t *rowptr, const int_t *colidx, const T_Scalar *values)
{
    return csx::ptr_max_asum(ptype, m, n, rowptr, colidx, values);
}

template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type 
norm_inf(prop_t ptype, int_t m, const int_t *rowptr, const int_t *colidx, const T_Scalar *values)
{
    return csx::idx_max_asum(ptype, m, rowptr, colidx, values);
}

template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type 
norm_max(int_t m, const int_t *rowptr, const T_Scalar *values)
{
    return csx::norm_max(m, rowptr, values);
}

template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type 
norm_fro(prop_t ptype, int_t m, const int_t *rowptr, const int_t *colidx, const T_Scalar *values)
{
    return csx::norm_fro(ptype, m, rowptr, colidx, values);
}

template <typename T_Scalar>
void permute(prop_t ptype, uplo_t uplo, int_t m, int_t n, 
		     const int_t *irowptr, const int_t *icolidx, const T_Scalar *ivalues,
		     int_t *orowptr, int_t *ocolidx, T_Scalar *ovalues, const int_t *P, const int_t *Q);

/*-------------------------------------------------*/
} // namespace csr
} // namespace blk
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_BULK_CSR_HPP_
