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
#include "culite/algebra/functional_inner.hpp"

// system

// 3rd

// culite
#include "culite/dense/dns_xxvector.hpp"

// forwards
#include "culite/checks/cla3p_forwards.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace ops {
/*-------------------------------------------------*/
template <typename T_Scalar>
T_Scalar dot(const dns::XxVector<T_Scalar>& x, const dns::XxVector<T_Scalar>& y, CuBlasHandler& cublasHandler)
{
    ::cla3p::similarity_dim_check(x.size(), y.size());

    T_Scalar ret = makeScalar<T_Scalar>(0);
    cublasHandler.dot(x.size(), x.values(), 1, y.values(), 1, &ret);
    return ret;
}
/*-------------------------------------------------*/
#define instantiate_dot(T_Scl) \
template T_Scl dot(const dns::XxVector<T_Scl>&, const dns::XxVector<T_Scl>&, CuBlasHandler&)
instantiate_dot(real_t);
instantiate_dot(real4_t);
instantiate_dot(complex_t);
instantiate_dot(complex8_t);
#undef instantiate_dot
/*-------------------------------------------------*/
template <typename T_Scalar>
T_Scalar dotc(const dns::XxVector<T_Scalar>& x, const dns::XxVector<T_Scalar>& y, CuBlasHandler& cublasHandler)
{
    ::cla3p::similarity_dim_check(x.size(), y.size());

    T_Scalar ret = makeScalar<T_Scalar>(0);
    cublasHandler.dotc(x.size(), x.values(), 1, y.values(), 1, &ret);
    return ret;
}
/*-------------------------------------------------*/
#define instantiate_dotc(T_Scl) \
template T_Scl dotc(const dns::XxVector<T_Scl>&, const dns::XxVector<T_Scl>&, CuBlasHandler&)
instantiate_dotc(real_t);
instantiate_dotc(real4_t);
instantiate_dotc(complex_t);
instantiate_dotc(complex8_t);
#undef instantiate_dotc
/*-------------------------------------------------*/
} // namespace ops
} // namespace culite
/*-------------------------------------------------*/
