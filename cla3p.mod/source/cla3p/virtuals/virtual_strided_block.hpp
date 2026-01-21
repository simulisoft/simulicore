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

#ifndef CLA3P_VIRTUAL_STRIDED_BLOCK_HPP_
#define CLA3P_VIRTUAL_STRIDED_BLOCK_HPP_

/**
 * @file
 */

#include "cla3p/types/integer.hpp"
#include "cla3p/generic/meta1d.hpp"
#include "cla3p/generic/matrix_meta.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

namespace dns { template <typename T_Scalar> class XxVector; }
namespace dns { template <typename T_Scalar> class XxMatrix; }
namespace csc { template <typename T_Int, typename T_Scalar> class XxMatrix; }

/*-------------------------------------------------*/

template <typename T_Scalar>
class VirtualStridedDnsVector : public Meta1D<int_t> {

    public:
        VirtualStridedDnsVector(int_t n, const T_Scalar *vals, int_t incv); 
        ~VirtualStridedDnsVector();

        void evaluateOnNew(dns::XxVector<T_Scalar>& dest) const;
		void evaluateOnExisting(dns::XxVector<T_Scalar>& dest) const;

    private:
        const T_Scalar *m_values;
        int_t m_incv;
};

/*-------------------------------------------------*/

template <typename T_Scalar>
class VirtualStridedDnsMatrix : public MatrixMeta<int_t> {

    public:
        VirtualStridedDnsMatrix(int_t m, int_t n, const T_Scalar *vals, int_t ldv, 
                           const Property& pr, int_t incv); 
        ~VirtualStridedDnsMatrix();

        void evaluateOnNew(dns::XxMatrix<T_Scalar>& dest) const;
		void evaluateOnExisting(dns::XxMatrix<T_Scalar>& dest) const;

    private:
        const T_Scalar *m_values;
        int_t m_incv;
        int_t m_ldv;
};

/*-------------------------------------------------*/

template <typename T_Int, typename T_Scalar>
class VirtualStridedCscMatrix : MatrixMeta<T_Int> {

    public:
        VirtualStridedCscMatrix(T_Int m, T_Int n, 
                           const T_Int *cptr, const T_Int *ridx, const T_Scalar *vals, 
                           const Property& pr, T_Int incv); 
        ~VirtualStridedCscMatrix();

        void evaluateOnNew(csc::XxMatrix<T_Int,T_Scalar>& dest) const;
		void evaluateOnExisting(csc::XxMatrix<T_Int,T_Scalar>& dest) const;

    private:
        const T_Int *m_colptr;
        const T_Int *m_rowidx;
        const T_Scalar *m_values;
        T_Int m_incv;

        T_Int nnz() const { return m_colptr[this->ncols()] ; }
};

/*-------------------------------------------------*/

namespace alias {

template <typename T_Scalar>
using VirtualStrided_vec = VirtualStridedDnsVector<T_Scalar>;

template <typename T_Scalar>
using VirtualStrided_mat = VirtualStridedDnsMatrix<T_Scalar>;

template <typename T_Int, typename T_Scalar>
using VirtualStrided_csc = VirtualStridedCscMatrix<T_Int,T_Scalar>;

} // namespace alias

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_VIRTUAL_STRIDED_BLOCK_HPP_
