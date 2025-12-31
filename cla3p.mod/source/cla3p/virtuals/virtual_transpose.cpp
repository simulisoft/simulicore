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
#include "cla3p/virtuals/virtual_transpose.hpp"

// system

// 3rd

// cla3p
#include "cla3p/bulk/dns.hpp"
#include "cla3p/bulk/csc.hpp"
#include "cla3p/dense/dns_xxvector.hpp"
#include "cla3p/sparse/csc_xxmatrix.hpp"
#include "cla3p/algebra/functional_update.hpp"

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualTransposeEvaluateOnNewSpec(
		const dns::XxMatrix<T_Scalar>& src,
		bool conj,
		dns::XxMatrix<T_Scalar>& dest)
{
	dest.clear();
	dest = dns::XxMatrix<T_Scalar>(src.ncols(), src.nrows(), src.prop().transpose());
	VirtualTransposeEvaluateOnExistingSpec(src, conj, dest);
}
/*-------------------------------------------------*/
template void VirtualTransposeEvaluateOnNewSpec(const dns::XxMatrix<real_t    >& , bool, dns::XxMatrix<real_t    >&);
template void VirtualTransposeEvaluateOnNewSpec(const dns::XxMatrix<real4_t   >& , bool, dns::XxMatrix<real4_t   >&);
template void VirtualTransposeEvaluateOnNewSpec(const dns::XxMatrix<complex_t >& , bool, dns::XxMatrix<complex_t >&);
template void VirtualTransposeEvaluateOnNewSpec(const dns::XxMatrix<complex8_t>& , bool, dns::XxMatrix<complex8_t>&);
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualTransposeEvaluateOnExistingSpec(
		const dns::XxMatrix<T_Scalar>& src,
		bool conj,
		dns::XxMatrix<T_Scalar>& dest)
{
	if(conj) {
		blk::dns::conjugate_transpose(
				src.nrows(),
				src.ncols(),
				src.values(),
				src.ld(),
				dest.values(),
				dest.ld());
	} else {
		blk::dns::transpose(
				src.nrows(),
				src.ncols(),
				src.values(),
				src.ld(),
				dest.values(), 
				dest.ld());
	} // conj
}
/*-------------------------------------------------*/
template void VirtualTransposeEvaluateOnExistingSpec(const dns::XxMatrix<real_t    >&, bool, dns::XxMatrix<real_t    >&);
template void VirtualTransposeEvaluateOnExistingSpec(const dns::XxMatrix<real4_t   >&, bool, dns::XxMatrix<real4_t   >&);
template void VirtualTransposeEvaluateOnExistingSpec(const dns::XxMatrix<complex_t >&, bool, dns::XxMatrix<complex_t >&);
template void VirtualTransposeEvaluateOnExistingSpec(const dns::XxMatrix<complex8_t>&, bool, dns::XxMatrix<complex8_t>&);
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualTransposeAccumulateOnExistingSpec(
		const dns::XxMatrix<T_Scalar>& src,
		bool conj,
		dns::XxMatrix<T_Scalar>& dest,
		T_Scalar coeff)
{
	// TODO: add without evaluating to temp
	dns::XxMatrix<T_Scalar> tmp;
	VirtualTransposeEvaluateOnNewSpec(src, conj, tmp);
	ops::update(coeff, tmp, dest);
}
/*-------------------------------------------------*/
template void VirtualTransposeAccumulateOnExistingSpec(const dns::XxMatrix<real_t    >&, bool, dns::XxMatrix<real_t    >&, real_t    );
template void VirtualTransposeAccumulateOnExistingSpec(const dns::XxMatrix<real4_t   >&, bool, dns::XxMatrix<real4_t   >&, real4_t   );
template void VirtualTransposeAccumulateOnExistingSpec(const dns::XxMatrix<complex_t >&, bool, dns::XxMatrix<complex_t >&, complex_t );
template void VirtualTransposeAccumulateOnExistingSpec(const dns::XxMatrix<complex8_t>&, bool, dns::XxMatrix<complex8_t>&, complex8_t);
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void VirtualTransposeEvaluateOnNewSpec(
		const csc::XxMatrix<T_Int,T_Scalar>& src,
		bool conj,
		csc::XxMatrix<T_Int, T_Scalar>& dest)
{
	dest.clear();
	dest = csc::XxMatrix<T_Int,T_Scalar>(src.ncols(), src.nrows(), src.nnz(), src.prop().transpose());
	VirtualTransposeEvaluateOnExistingSpec(src, conj, dest);
}
/*-------------------------------------------------*/
template void VirtualTransposeEvaluateOnNewSpec(const csc::XxMatrix<int_t,real_t    >&, bool, csc::XxMatrix<int_t,real_t    >&);
template void VirtualTransposeEvaluateOnNewSpec(const csc::XxMatrix<int_t,real4_t   >&, bool, csc::XxMatrix<int_t,real4_t   >&);
template void VirtualTransposeEvaluateOnNewSpec(const csc::XxMatrix<int_t,complex_t >&, bool, csc::XxMatrix<int_t,complex_t >&);
template void VirtualTransposeEvaluateOnNewSpec(const csc::XxMatrix<int_t,complex8_t>&, bool, csc::XxMatrix<int_t,complex8_t>&);
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void VirtualTransposeEvaluateOnExistingSpec(
		const csc::XxMatrix<T_Int,T_Scalar>& src,
		bool conj,
		csc::XxMatrix<T_Int, T_Scalar>& dest)
{
	if(conj) {
		blk::csc::conjugate_transpose(
				src.nrows(),
				src.ncols(),
				src.colptr(),
				src.rowidx(),
				src.values(),
				dest.colptr(),
				dest.rowidx(),
				dest.values());
	} else {
		blk::csc::transpose(
				src.nrows(),
				src.ncols(),
				src.colptr(),
				src.rowidx(),
				src.values(),
				dest.colptr(),
				dest.rowidx(),
				dest.values());
	} // conj
}
/*-------------------------------------------------*/
template void VirtualTransposeEvaluateOnExistingSpec(const csc::XxMatrix<int_t,real_t    >&, bool, csc::XxMatrix<int_t,real_t    >&);
template void VirtualTransposeEvaluateOnExistingSpec(const csc::XxMatrix<int_t,real4_t   >&, bool, csc::XxMatrix<int_t,real4_t   >&);
template void VirtualTransposeEvaluateOnExistingSpec(const csc::XxMatrix<int_t,complex_t >&, bool, csc::XxMatrix<int_t,complex_t >&);
template void VirtualTransposeEvaluateOnExistingSpec(const csc::XxMatrix<int_t,complex8_t>&, bool, csc::XxMatrix<int_t,complex8_t>&);
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void VirtualTransposeAccumulateOnExistingSpec(
		const csc::XxMatrix<T_Int,T_Scalar>& src,
		bool conj,
		csc::XxMatrix<T_Int, T_Scalar>& dest,
		T_Scalar coeff)
{
	csc::XxMatrix<T_Int,T_Scalar> tmp;
	VirtualTransposeEvaluateOnNewSpec(src, conj, tmp);
	ops::update(coeff, tmp, dest);
}
/*-------------------------------------------------*/
template void VirtualTransposeAccumulateOnExistingSpec(const csc::XxMatrix<int_t,real_t    >&, bool, csc::XxMatrix<int_t,real_t    >&, real_t    );
template void VirtualTransposeAccumulateOnExistingSpec(const csc::XxMatrix<int_t,real4_t   >&, bool, csc::XxMatrix<int_t,real4_t   >&, real4_t   );
template void VirtualTransposeAccumulateOnExistingSpec(const csc::XxMatrix<int_t,complex_t >&, bool, csc::XxMatrix<int_t,complex_t >&, complex_t );
template void VirtualTransposeAccumulateOnExistingSpec(const csc::XxMatrix<int_t,complex8_t>&, bool, csc::XxMatrix<int_t,complex8_t>&, complex8_t);
/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/
