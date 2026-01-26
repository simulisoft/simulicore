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
#include "cla3p/bulk/csc.hpp"

// system

// 3rd

// cla3p
#include "cla3p/error/exceptions.hpp"
#include "cla3p/checks/basic_checks.hpp"
#include "cla3p/support/utils.hpp"


/*-------------------------------------------------*/
namespace cla3p {
namespace blk {
namespace csc {
/*-------------------------------------------------*/
template <typename T_Scalar>
static void permute_ge_both(int_t /*m*/, int_t n,
		const int_t *icolptr, const int_t *irowidx, const T_Scalar *ivalues,
		int_t *ocolptr, int_t *orowidx, T_Scalar *ovalues, const int_t *P, const int_t *Q)
{
	ocolptr[0] = 0;
	for(int_t j = 0; j < n; j++) {
		ocolptr[j + 1] = icolptr[Q[j]+1] - icolptr[Q[j]];
	} // j

	csx::roll(n, ocolptr);

	for(int_t j = 0; j < n; j++) {
		for(int_t irow = icolptr[Q[j]]; irow < icolptr[Q[j]+1]; irow++) {
			orowidx[ocolptr[j]] = P[irowidx[irow]];
			ovalues[ocolptr[j]] = ivalues[irow];
			ocolptr[j]++;
		} // irow
	} // j
	
	csx::unroll(n, ocolptr);

	csx::sort(n, ocolptr, orowidx, ovalues);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
static void permute_ge_left(int_t /*m*/, int_t n,
		const int_t *icolptr, const int_t *irowidx, const T_Scalar *ivalues,
		int_t *ocolptr, int_t *orowidx, T_Scalar *ovalues, const int_t *P)
{
	std::copy(icolptr, icolptr + n + 1, ocolptr);

	for(int_t j = 0; j < n; j++) {
		for(int_t irow = icolptr[j]; irow < icolptr[j+1]; irow++) {
			orowidx[ocolptr[j]] = P[irowidx[irow]];
			ovalues[ocolptr[j]] = ivalues[irow];
			ocolptr[j]++;
		} // irow
	} // j
	
	csx::unroll(n, ocolptr);

	csx::sort(n, ocolptr, orowidx, ovalues);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
static void permute_ge_right(int_t /*m*/, int_t n,
		const int_t *icolptr, const int_t *irowidx, const T_Scalar *ivalues,
		int_t *ocolptr, int_t *orowidx, T_Scalar *ovalues, const int_t *Q)
{
	ocolptr[0] = 0;
	for(int_t j = 0; j < n; j++) {
		ocolptr[j + 1] = icolptr[Q[j] + 1] - icolptr[Q[j]];
	} // j

	csx::roll(n, ocolptr);

	for(int_t j = 0; j < n; j++) {
		std::copy(irowidx + icolptr[Q[j]], irowidx + icolptr[Q[j] + 1], orowidx + ocolptr[j]);
		std::copy(ivalues + icolptr[Q[j]], ivalues + icolptr[Q[j] + 1], ovalues + ocolptr[j]);
	} // j
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void permute_xx_mirror(prop_t ptype, uplo_t uplo, int_t n,
		const int_t *icolptr, const int_t *irowidx, const T_Scalar *ivalues,
		int_t *ocolptr, int_t *orowidx, T_Scalar *ovalues, const int_t *P)
{
	int_t Pi;
	int_t Pj;

	for(int_t j = 0; j < n + 1; j++) 
		ocolptr[j] = 0;

	for(int_t j = 0; j < n; j++) {
		for(int_t irow = icolptr[j]; irow < icolptr[j+1]; irow++) {
			Pi = P[irowidx[irow]];
			Pj = P[j];
			if(uplo == uplo_t::Upper && Pj < Pi) 
				ocolptr[Pi+1]++;
			else if(uplo == uplo_t::Lower && Pj > Pi) 
				ocolptr[Pi+1]++;
			else
				ocolptr[Pj+1]++;
		} // irow
	} // j

	csx::roll(n, ocolptr);

	for(int_t j = 0; j < n; j++) {
		for(int_t irow = icolptr[j]; irow < icolptr[j+1]; irow++) {
			Pi = P[irowidx[irow]];
			Pj = P[j];
			if(uplo == uplo_t::Upper && Pj < Pi) {
				orowidx[ocolptr[Pi]] = Pj;
				ovalues[ocolptr[Pi]] = opposite_element(ivalues[irow],ptype);
				ocolptr[Pi]++;
			} else if(uplo == uplo_t::Lower && Pj > Pi) {
				orowidx[ocolptr[Pi]] = Pj;
				ovalues[ocolptr[Pi]] = opposite_element(ivalues[irow],ptype);
				ocolptr[Pi]++;
			} else {
				orowidx[ocolptr[Pj]] = Pi;
				ovalues[ocolptr[Pj]] = ivalues[irow];
				ocolptr[Pj]++;
			}
		} // irow
	} // j
	
	csx::unroll(n, ocolptr);

	csx::sort(n, ocolptr, orowidx, ovalues);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void permute(prop_t ptype, uplo_t uplo, int_t m, int_t n,
		const int_t *icolptr, const int_t *irowidx, const T_Scalar *ivalues,
		int_t *ocolptr, int_t *orowidx, T_Scalar *ovalues, const int_t *P, const int_t *Q)
{
	if(!m || !n) return;

	Property prop(ptype, uplo);

	if(prop.isSquare()) {
		square_check(m, n);
	}

	if(prop.isGeneral()) {

		if(P && Q) {

			permute_ge_both(m, n, icolptr, irowidx, ivalues, ocolptr, orowidx, ovalues, P, Q);

		} else if(P && !Q) {

			permute_ge_left(m, n, icolptr, irowidx, ivalues, ocolptr, orowidx, ovalues, P);

		} else if(!P && Q) { 

			permute_ge_right(m, n, icolptr, irowidx, ivalues, ocolptr, orowidx, ovalues, Q);

		} else {

			std::copy(icolptr, icolptr + n + 1    , ocolptr);
			std::copy(irowidx, irowidx + icolptr[n], orowidx);
			std::copy(ivalues, ivalues + icolptr[n], ovalues);

		} // P/Q

	} else if(prop.isSymmetric() || prop.isHermitian() || prop.isSkew()) {

		if(P) {

			permute_xx_mirror(ptype, uplo, n, icolptr, irowidx, ivalues, ocolptr, orowidx, ovalues, P);

		} else {

			std::copy(icolptr, icolptr + n + 1    , ocolptr);
			std::copy(irowidx, irowidx + icolptr[n], orowidx);
			std::copy(ivalues, ivalues + icolptr[n], ovalues);

		} // P

	} else {

		throw err::Exception("Invalid property: " + prop.name());

	} // prop
}
/*-------------------------------------------------*/
#define instantiate_permute(int_t, T_Scl) \
template void permute(prop_t, uplo_t, int_t, int_t, \
                      const int_t*, const int_t*, const T_Scl*, \
                      int_t*, int_t*, T_Scl*, const int_t*, const int_t*)
instantiate_permute(int_t, real_t    );
instantiate_permute(int_t, real4_t   );
instantiate_permute(int_t, complex_t );
instantiate_permute(int_t, complex8_t);
#undef instantiate_permute
/*-------------------------------------------------*/
} // namespace csc
} // namespace blk
} // namespace cla3p
/*-------------------------------------------------*/
