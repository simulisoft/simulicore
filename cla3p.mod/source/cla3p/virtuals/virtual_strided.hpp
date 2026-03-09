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

#ifndef CLA3P_VIRTUAL_STRIDED_HPP_
#define CLA3P_VIRTUAL_STRIDED_HPP_

/**
 * @file
 */

#include "cla3p/types/integer.hpp"
#include "cla3p/generic/meta1d.hpp"
#include "cla3p/generic/matrix_meta.hpp"
#include "cla3p/generic/guard.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

namespace dns { template <typename T_Scalar> class XxVector; }
namespace dns { template <typename T_Scalar> class XxMatrix; }
namespace csr { template <typename T_Int, typename T_Scalar> class XxMatrix; }
namespace csc { template <typename T_Int, typename T_Scalar> class XxMatrix; }

/*-------------------------------------------------*/

template <typename T_Scalar>
class VirtualStridedDnsVector : public Meta1D<int_t> {

    public:
        VirtualStridedDnsVector();
        VirtualStridedDnsVector(int_t n, T_Scalar *vals, int_t incv); 
        ~VirtualStridedDnsVector();

        VirtualStridedDnsVector(const VirtualStridedDnsVector& other) = delete;
        VirtualStridedDnsVector(VirtualStridedDnsVector&& other) = default;
        VirtualStridedDnsVector& operator=(VirtualStridedDnsVector&& other) = default;

        VirtualStridedDnsVector& operator=(const VirtualStridedDnsVector& other);
        void operator=(const dns::XxVector<T_Scalar>& vec);
        void operator=(T_Scalar val);

        void clear();
        VirtualStridedDnsVector<T_Scalar> rcopy();

        static Guard<VirtualStridedDnsVector<T_Scalar>> 
        view(int_t n, const T_Scalar *vals, int_t incv);

        void evaluateOnNew(dns::XxVector<T_Scalar>& dest) const;
        void evaluateOnExisting(dns::XxVector<T_Scalar>& dest) const;

    private:
        T_Scalar *m_values;
        int_t m_incv;

        void defaults();
};

/*-------------------------------------------------*/

template <typename T_Scalar>
class VirtualStridedDnsMatrix : public MatrixMeta<int_t> {

    public:
        VirtualStridedDnsMatrix();
        VirtualStridedDnsMatrix(int_t m, int_t n, 
                                T_Scalar *vals, int_t ldv, 
                                const Property& pr, int_t incv); 
        ~VirtualStridedDnsMatrix();

        VirtualStridedDnsMatrix(const VirtualStridedDnsMatrix& other) = delete;
        VirtualStridedDnsMatrix(VirtualStridedDnsMatrix&& other) = default;
        VirtualStridedDnsMatrix& operator=(VirtualStridedDnsMatrix&& other) = default;

        VirtualStridedDnsMatrix& operator=(const VirtualStridedDnsMatrix& other);
        void operator=(const dns::XxMatrix<T_Scalar>& mat);
        void operator=(T_Scalar val);

        void clear();
        VirtualStridedDnsMatrix<T_Scalar> rcopy();

        static Guard<VirtualStridedDnsMatrix<T_Scalar>> 
        view(int_t m, int_t n, 
             const T_Scalar *vals, int_t ldv, 
             const Property& pr, int_t incv);

        void evaluateOnNew(dns::XxMatrix<T_Scalar>& dest) const;
        void evaluateOnExisting(dns::XxMatrix<T_Scalar>& dest) const;

    private:
        T_Scalar *m_values;
        int_t m_ldv;
        int_t m_incv;

        void defaults();
};

/*-------------------------------------------------*/

template <typename T_Int, typename T_Scalar>
class VirtualStridedCsrMatrix : MatrixMeta<T_Int> {

    public:
        VirtualStridedCsrMatrix();
        VirtualStridedCsrMatrix(T_Int m, T_Int n, 
                                T_Int *rptr, T_Int *cidx, T_Scalar *vals, 
                                const Property& pr, T_Int incv); 
        ~VirtualStridedCsrMatrix();

        VirtualStridedCsrMatrix(const VirtualStridedCsrMatrix& other) = delete;
        VirtualStridedCsrMatrix(VirtualStridedCsrMatrix&& other) = default;
        VirtualStridedCsrMatrix& operator=(VirtualStridedCsrMatrix&& other) = default;

        VirtualStridedCsrMatrix<T_Int,T_Scalar>& operator=(const VirtualStridedCsrMatrix<T_Int,T_Scalar>& other);
        void operator=(const csr::XxMatrix<T_Int,T_Scalar>& mat);
        void operator=(T_Scalar val);

        void clear();
        VirtualStridedCsrMatrix<T_Int,T_Scalar> rcopy();

        static Guard<VirtualStridedCsrMatrix<T_Int,T_Scalar>> 
        view(T_Int m, T_Int n, 
             const T_Int *rptr, const T_Int *cidx, const T_Scalar *vals, 
             const Property& pr, T_Int incv);

        void evaluateOnNew(csr::XxMatrix<T_Int,T_Scalar>& dest) const;
        void evaluateOnExisting(csr::XxMatrix<T_Int,T_Scalar>& dest) const;

    private:
        T_Int *m_rowptr;
        T_Int *m_colidx;
        T_Scalar *m_values;
        T_Int m_incv;

        void defaults();

        T_Int nnz() const { return m_rowptr[this->nrows()] ; }
};

/*-------------------------------------------------*/

template <typename T_Int, typename T_Scalar>
class VirtualStridedCscMatrix : MatrixMeta<T_Int> {

    public:
        VirtualStridedCscMatrix();
        VirtualStridedCscMatrix(T_Int m, T_Int n, 
                                T_Int *cptr, T_Int *ridx, T_Scalar *vals, 
                                const Property& pr, T_Int incv); 
        ~VirtualStridedCscMatrix();

        VirtualStridedCscMatrix(const VirtualStridedCscMatrix& other) = delete;
        VirtualStridedCscMatrix(VirtualStridedCscMatrix&& other) = default;
        VirtualStridedCscMatrix& operator=(VirtualStridedCscMatrix&& other) = default;

        VirtualStridedCscMatrix<T_Int,T_Scalar>& operator=(const VirtualStridedCscMatrix<T_Int,T_Scalar>& other);
        void operator=(const csc::XxMatrix<T_Int,T_Scalar>& mat);
        void operator=(T_Scalar val);

        void clear();
        VirtualStridedCscMatrix<T_Int,T_Scalar> rcopy();

        static Guard<VirtualStridedCscMatrix<T_Int,T_Scalar>> 
        view(T_Int m, T_Int n, 
             const T_Int *cptr, const T_Int *ridx, const T_Scalar *vals, 
             const Property& pr, T_Int incv);

        void evaluateOnNew(csc::XxMatrix<T_Int,T_Scalar>& dest) const;
        void evaluateOnExisting(csc::XxMatrix<T_Int,T_Scalar>& dest) const;

    private:
        T_Int *m_colptr;
        T_Int *m_rowidx;
        T_Scalar *m_values;
        T_Int m_incv;

        void defaults();

        T_Int nnz() const { return m_colptr[this->ncols()] ; }
};

/*-------------------------------------------------*/

namespace alias {

template <typename T_Scalar>
using VirtualStrided_vec = VirtualStridedDnsVector<T_Scalar>;
template <typename T_Scalar>
using GuardedStrided_vec = Guard<VirtualStrided_vec<T_Scalar>>;

template <typename T_Scalar>
using VirtualStrided_dns = VirtualStridedDnsMatrix<T_Scalar>;
template <typename T_Scalar>
using GuardedStrided_dns = Guard<VirtualStrided_dns<T_Scalar>>;

template <typename T_Int, typename T_Scalar>
using VirtualStrided_csr = VirtualStridedCsrMatrix<T_Int,T_Scalar>;
template <typename T_Int, typename T_Scalar>
using GuardedStrided_csr = Guard<VirtualStrided_csr<T_Int,T_Scalar>>;

template <typename T_Int, typename T_Scalar>
using VirtualStrided_csc = VirtualStridedCscMatrix<T_Int,T_Scalar>;
template <typename T_Int, typename T_Scalar>
using GuardedStrided_csc = Guard<VirtualStrided_csc<T_Int,T_Scalar>>;

} // namespace alias

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_VIRTUAL_STRIDED_HPP_
