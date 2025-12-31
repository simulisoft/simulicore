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
#include "cla3p/algebra/functional_inner.hpp"

// system

// 3rd

// cla3p
#include "cla3p/checks/basic_checks.hpp"
#include "cla3p/proxies/blas_proxy.hpp"
#include "cla3p/dense/dns_xxvector.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace ops {
/*-------------------------------------------------*/
template <typename T_Scalar>
T_Scalar dot(const dns::XxVector<T_Scalar>& X, const dns::XxVector<T_Scalar>& Y)
{
	similarity_dim_check(X.size(), Y.size());

	return blas::dot(X.size(), X.values(), 1, Y.values(), 1);
}
/*-------------------------------------------------*/
#define instantiate_dot(T_Scl) \
template T_Scl dot(const dns::XxVector<T_Scl>&, const dns::XxVector<T_Scl>&)
instantiate_dot(real_t);
instantiate_dot(real4_t);
instantiate_dot(complex_t);
instantiate_dot(complex8_t);
#undef instantiate_dot
/*-------------------------------------------------*/
template <typename T_Scalar>
T_Scalar dotc(const dns::XxVector<T_Scalar>& X, const dns::XxVector<T_Scalar>& Y)
{
	similarity_dim_check(X.size(), Y.size());

	return blas::dotc(X.size(), X.values(), 1, Y.values(), 1);
}
/*-------------------------------------------------*/
#define instantiate_dotc(T_Scl) \
template T_Scl dotc(const dns::XxVector<T_Scl>&, const dns::XxVector<T_Scl>&)
instantiate_dotc(real_t);
instantiate_dotc(real4_t);
instantiate_dotc(complex_t);
instantiate_dotc(complex8_t);
#undef instantiate_dotc
/*-------------------------------------------------*/
} // namespace ops
} // namespace cla3p
/*-------------------------------------------------*/
