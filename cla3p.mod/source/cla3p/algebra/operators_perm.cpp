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
#include "cla3p/algebra/operators_perm.hpp"

// system

// 3rd

// cla3p

/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
cla3p::dns::XxVector<T_Scalar> operator*(
    const cla3p::prm::PxMatrix<T_Int>& P,
    const cla3p::dns::XxVector<T_Scalar>& X)
{
  return X.permuteLeft(P);
}
/*-------------------------------------------------*/
#define instantiate_perm(T_Int, T_Scl) \
template cla3p::dns::XxVector<T_Scl> operator*( \
		const cla3p::prm::PxMatrix<T_Int>&, \
		const cla3p::dns::XxVector<T_Scl>&)
instantiate_perm(cla3p::int_t, cla3p::real_t);
instantiate_perm(cla3p::int_t, cla3p::real4_t);
instantiate_perm(cla3p::int_t, cla3p::complex_t);
instantiate_perm(cla3p::int_t, cla3p::complex8_t);
#undef instantiate_perm
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scl>
cla3p::dns::XxMatrix<T_Scl> operator*(
    const cla3p::prm::PxMatrix<T_Int>& P,
    const cla3p::dns::XxMatrix<T_Scl>& A)
{
  return A.permuteLeft(P);
}
/*-------------------------------------------------*/
#define instantiate_perm(T_Int, T_Scl) \
template cla3p::dns::XxMatrix<T_Scl> operator*( \
		const cla3p::prm::PxMatrix<T_Int>&, \
		const cla3p::dns::XxMatrix<T_Scl>&)
instantiate_perm(cla3p::int_t, cla3p::real_t);
instantiate_perm(cla3p::int_t, cla3p::real4_t);
instantiate_perm(cla3p::int_t, cla3p::complex_t);
instantiate_perm(cla3p::int_t, cla3p::complex8_t);
#undef instantiate_perm
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scl>
cla3p::dns::XxMatrix<T_Scl> operator*(
		const cla3p::dns::XxMatrix<T_Scl>& A,
		const cla3p::prm::PxMatrix<T_Int>& P)
{
	return A.permuteRight(P);
}
/*-------------------------------------------------*/
#define instantiate_perm(T_Int, T_Scl) \
template cla3p::dns::XxMatrix<T_Scl> operator*( \
		const cla3p::dns::XxMatrix<T_Scl>&, \
		const cla3p::prm::PxMatrix<T_Int>&)
instantiate_perm(cla3p::int_t, cla3p::real_t);
instantiate_perm(cla3p::int_t, cla3p::real4_t);
instantiate_perm(cla3p::int_t, cla3p::complex_t);
instantiate_perm(cla3p::int_t, cla3p::complex8_t);
#undef instantiate_perm
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scl>
cla3p::csc::XxMatrix<T_Int, T_Scl> operator*(
		const cla3p::prm::PxMatrix<T_Int>& P,
		const cla3p::csc::XxMatrix<T_Int, T_Scl>& A)
{
	return A.permuteLeft(P);
}
/*-------------------------------------------------*/
#define instantiate_perm(T_Int, T_Scl) \
template cla3p::csc::XxMatrix<T_Int, T_Scl> operator*( \
		const cla3p::prm::PxMatrix<T_Int>&, \
		const cla3p::csc::XxMatrix<T_Int, T_Scl>&)
instantiate_perm(cla3p::int_t, cla3p::real_t);
instantiate_perm(cla3p::int_t, cla3p::real4_t);
instantiate_perm(cla3p::int_t, cla3p::complex_t);
instantiate_perm(cla3p::int_t, cla3p::complex8_t);
#undef instantiate_perm
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scl>
cla3p::csc::XxMatrix<T_Int, T_Scl> operator*(
		const cla3p::csc::XxMatrix<T_Int, T_Scl>& A,
		const cla3p::prm::PxMatrix<T_Int>& P)
{
	return A.permuteRight(P);
}
/*-------------------------------------------------*/
#define instantiate_perm(T_Int, T_Scl) \
template cla3p::csc::XxMatrix<T_Int, T_Scl> operator*( \
		const cla3p::csc::XxMatrix<T_Int, T_Scl>&, \
		const cla3p::prm::PxMatrix<T_Int>&)
instantiate_perm(cla3p::int_t, cla3p::real_t);
instantiate_perm(cla3p::int_t, cla3p::real4_t);
instantiate_perm(cla3p::int_t, cla3p::complex_t);
instantiate_perm(cla3p::int_t, cla3p::complex8_t);
#undef instantiate_perm
/*-------------------------------------------------*/
template <typename T_Int>
cla3p::prm::PxMatrix<T_Int> operator*(
    const cla3p::prm::PxMatrix<T_Int>& P,
    const cla3p::prm::PxMatrix<T_Int>& Q)
{
  return Q.permuteLeft(P);
}
/*-------------------------------------------------*/
#define instantiate_perm(T_Int) \
template cla3p::prm::PxMatrix<T_Int> operator*( \
		const cla3p::prm::PxMatrix<T_Int>&, \
		const cla3p::prm::PxMatrix<T_Int>&)
instantiate_perm(cla3p::int_t);
#undef instantiate_perm
/*-------------------------------------------------*/
