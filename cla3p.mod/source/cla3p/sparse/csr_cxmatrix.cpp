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

// this file inc
#include <cla3p/sparse/csr_cxmatrix.hpp>

// system

// 3rd

// cla3p

/*-------------------------------------------------*/
namespace cla3p {
namespace csr {
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
CxMatrix<T_Int,T_Scalar>::CxMatrix(const XxMatrix<T_Int,T_Scalar>& other)
    : XxMatrix<T_Int,T_Scalar>(other)
{
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
CxMatrix<T_Int,T_Scalar>& CxMatrix<T_Int,T_Scalar>::operator=(const XxMatrix<T_Int,T_Scalar>& other)
{
    XxMatrix<T_Int,T_Scalar>::operator=(other);
    return *this;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
CxMatrix<T_Int,T_Scalar>::CxMatrix(XxMatrix<T_Int,T_Scalar>&& other)
    : XxMatrix<T_Int,T_Scalar>(std::move(other))
{
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
CxMatrix<T_Int,T_Scalar>& CxMatrix<T_Int,T_Scalar>::operator=(XxMatrix<T_Int,T_Scalar>&& other)
{
    XxMatrix<T_Int,T_Scalar>::operator=(std::move(other));
    return *this;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
CxMatrix<T_Int,T_Scalar>::CxMatrix()
{
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
CxMatrix<T_Int,T_Scalar>::CxMatrix(T_Int nr, T_Int nc, T_Int nz, const Property& pr)
    : CxMatrix<T_Int,T_Scalar>::XxMatrix(nr, nc, nz, pr)
{
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
CxMatrix<T_Int,T_Scalar>::CxMatrix(T_Int nr, T_Int nc, 
                                   T_Int *rptr, T_Int *cidx, T_Scalar *vals, 
                                   bool bind, const Property& pr)
    : CxMatrix<T_Int,T_Scalar>::XxMatrix(nr, nc, rptr, cidx, vals, bind, pr)
{
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
CxMatrix<T_Int,T_Scalar>::~CxMatrix()
{
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
alias::VirtualStrided_csr<T_Int, typename TypeTraits<T_Scalar>::real_type> 
CxMatrix<T_Int,T_Scalar>::real()
{
    T_RScalar *realValues = reinterpret_cast<T_RScalar*>(this->values());
    return alias::VirtualStrided_csr<T_Int,T_RScalar>(this->nrows(), 
                                                      this->ncols(),
                                                      this->rowptr(), 
                                                      this->colidx(),
                                                      realValues, 
                                                      this->prop(), 2);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
alias::VirtualStrided_csr<T_Int, typename TypeTraits<T_Scalar>::real_type> 
CxMatrix<T_Int,T_Scalar>::imag()
{
    T_RScalar *imagValues = reinterpret_cast<T_RScalar*>(this->values()) + 1;
    return alias::VirtualStrided_csr<T_Int,T_RScalar>(this->nrows(), 
                                                      this->ncols(),
                                                      this->rowptr(), 
                                                      this->colidx(),
                                                      imagValues, 
                                                      this->prop(), 2);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
alias::GuardedStrided_csr<T_Int, typename TypeTraits<T_Scalar>::real_type> 
CxMatrix<T_Int,T_Scalar>::real() const
{
    const T_RScalar *realValues = reinterpret_cast<const T_RScalar*>(this->values());
    return alias::VirtualStrided_csr<T_Int,T_RScalar>::view(this->nrows(), 
                                                            this->ncols(), 
                                                            this->rowptr(), 
                                                            this->colidx(),
                                                            realValues, 
                                                            this->prop(), 2);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
alias::GuardedStrided_csr<T_Int, typename TypeTraits<T_Scalar>::real_type> 
CxMatrix<T_Int,T_Scalar>::imag() const
{
    const T_RScalar *imagValues = reinterpret_cast<const T_RScalar*>(this->values()) + 1;
    return alias::VirtualStrided_csr<T_Int,T_RScalar>::view(this->nrows(), 
                                                            this->ncols(), 
                                                            this->rowptr(), 
                                                            this->colidx(),
                                                            imagValues, 
                                                            this->prop(), 2);
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template class CxMatrix<int_t,complex_t>;
template class CxMatrix<int_t,complex8_t>;
/*-------------------------------------------------*/
} // namespace csr
} // namespace cla3p
/*-------------------------------------------------*/
