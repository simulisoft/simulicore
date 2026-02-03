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
#include "cla3p/bulk/csx.hpp"

// system
#include <algorithm>
#include <cmath>
#include <vector>

// 3rd

// cla3p
#include "cla3p/error/exceptions.hpp"
#include "cla3p/proxies/blas_proxy.hpp"
#include "cla3p/support/utils.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace blk {
namespace csx {
/*-------------------------------------------------*/
void roll(int_t np, int_t *xxxptr)
{
	for(int_t k = 0; k < np; k++) {
		xxxptr[k+1] += xxxptr[k];
	} // k
}
/*-------------------------------------------------*/
void unroll(int_t np, int_t *xxxptr)
{
	if(np <= 0) return;

	for(int_t k = np-1; k > 0; k--) {
		xxxptr[k] = xxxptr[k-1];
	} // k

	xxxptr[0] = 0;
}
/*-------------------------------------------------*/
int_t maxXlen(int_t np, const int_t *xxxptr)
{
	int_t ret = 0;

	for(int_t k = 0; k < np; k++) {
		ret = std::max(ret, xxxptr[k+1] - xxxptr[k]);
	} // k

	if(ret < 0) {
		throw err::NoConsistency("Negative length found in pointer array");
	}

	return ret;
}
/*-------------------------------------------------*/
static std::string indexDimTypeStr(csx_t type, bool capitalize = false)
{
    switch(type) {
        case csx_t::CSR: return capitalize ? "Column" : "column";
        case csx_t::CSC: return capitalize ? "Row" : "row";
        default:         return "Unknown";
    }
}
/*-------------------------------------------------*/
static std::string pointerDimTypeStr(csx_t type, bool capitalize = false)
{
    switch(type) {
        case csx_t::CSR: return capitalize ? "Row" : "row";
        case csx_t::CSC: return capitalize ? "Column" : "column";
        default:         return "Unknown";
    }
}
/*-------------------------------------------------*/
void check(prop_t ptype, uplo_t uplo, csx_t type, int_t np, int_t maxIdx, const int_t *xxxptr, const int_t *xxxidx)
{
	Property prop(ptype, uplo);

	if(xxxptr[0]) {
		throw err::NoConsistency("Pointer array must contain a zero at position 0");
	}

	std::vector<bool> mark(maxIdx, false);

	for(int_t k = 0; k < np; k++) {

		int_t xbgn = xxxptr[k];
		int_t xend = xxxptr[k+1];

		if(xend < xbgn) {
			throw err::NoConsistency(pointerDimTypeStr(type, true) + " pointer array must be in ascending order");
		}

		for(int_t xidx = xbgn; xidx < xend; xidx++) {

            int_t l = xxxidx[xidx];

			if(l < 0) {
				throw err::OutOfBounds("Negative " + indexDimTypeStr(type) + 
                                       " index detected: " + std::to_string(l));
			}

            if(l >= maxIdx) {
				throw err::OutOfBounds(indexDimTypeStr(type, true) + " index " + std::to_string(l) + 
                                       " is greater than " + indexDimTypeStr(type) + "s number " + std::to_string(maxIdx));
			}

            int_t row = type == csx_t::CSC ? l : k;
            int_t col = type == csx_t::CSC ? k : l;

			if(prop.isLower() && row < col) {
				throw err::NoConsistency("Found coordinate " + coordToString(row,col) + " in upper part");
			}

			if(prop.isUpper() && row > col) {
				throw err::NoConsistency("Found coordinate " + coordToString(row,col) + " in lower part");
			}

			// TODO: perhaps check skew diagonals for values

			if(xidx > xbgn) {
				if(l < xxxidx[xidx-1]) {
					throw err::NoConsistency(pointerDimTypeStr(type, true) + " " + std::to_string(k) + " is not sorted");
				}
			}

			if(mark[l]) {
				throw err::NoConsistency("Duplicate entry detected at " + coordToString(row,col));
			}

			mark[l] = true;

		} // xidx

        //
        // reset marks
        //
		for(int_t xidx = xbgn; xidx < xend; xidx++) {
			mark[xxxidx[xidx]] = false;
		} // xidx

	} // k
}
/*-------------------------------------------------*/
void sort(int_t np, const int_t *xxxptr, int_t *xxxidx)
{
	for(int_t k = 0; k < np; k++) {

		int_t xbgn = xxxptr[k];
		int_t xend = xxxptr[k+1];
		int_t xlen = xend - xbgn;

		if(xlen > 0) {
			std::sort(xxxidx + xbgn, xxxidx + xend);
		} // xlen

	} // k
}
/*-------------------------------------------------*/
template <typename T_Scalar>
class IdxValuePair {
	public:
		IdxValuePair() : i(0), v(0) {}
		~IdxValuePair() = default;
		IdxValuePair(const IdxValuePair<T_Scalar>& p) : i(p.i), v(p.v) {} // Removes ABI warning
		int_t    i;
		T_Scalar v;

		bool operator<(const IdxValuePair<T_Scalar>& other) const
		{
			return (i < other.i);
		}
};
/*-------------------------------------------------*/
// 
// TODO: implement without extra storage
//
template <typename T_Scalar>
void sort(int_t np, const int_t *xxxptr, int_t *xxxidx, T_Scalar *values)
{
	if(np <= 0) return;

	int_t mxlen = maxXlen(np, xxxptr);

	if(!mxlen) return;

	std::vector<IdxValuePair<T_Scalar>> idx_pairs;
	idx_pairs.reserve(mxlen);

	for(int_t k = 0; k < np; k++) {

		int_t xbgn = xxxptr[k];
		int_t xend = xxxptr[k+1];
		int_t xlen = xend - xbgn;

		if(xlen < 2) continue;

		idx_pairs.resize(xlen);

		for(int_t xidx = xbgn; xidx < xend; xidx++) {
			idx_pairs[xidx - xbgn].i = xxxidx[xidx];
			idx_pairs[xidx - xbgn].v = values[xidx];
		} // xidx

		std::sort(idx_pairs.begin(), idx_pairs.end());

		for(int_t xidx = xbgn; xidx < xend; xidx++) {
			xxxidx[xidx] = idx_pairs[xidx - xbgn].i;
			values[xidx] = idx_pairs[xidx - xbgn].v;
		} // xidx

	} // k
}
/*-------------------------------------------------*/
template void sort(int_t, const int_t*, int_t*, real_t    *);
template void sort(int_t, const int_t*, int_t*, real4_t   *);
template void sort(int_t, const int_t*, int_t*, complex_t *);
template void sort(int_t, const int_t*, int_t*, complex8_t*);
/*-------------------------------------------------*/
template <typename T_Scalar>
void print_to_stream(std::ostream& os, csx_t type, int_t np, const int_t *xxxptr, const int_t *xxxidx, const T_Scalar *values, std::streamsize prec)
{
	if(!np) return;

	int_t maxIdx = 0;
	for(int_t k = 0; k < xxxptr[np]; k++) 
		maxIdx = std::max(maxIdx, xxxidx[k]);

	ListPrinter listPrinter(os, maxIdx, np, xxxptr[np], prec);
	listPrinter.streamHeader();

	for(int_t k = 0, cnt = 0; k < np; k++) {

		int_t xbgn = xxxptr[k];
		int_t xend = xxxptr[k+1];

		for(int_t xidx = xbgn; xidx < xend; xidx++) {
            int_t row = type == csx_t::CSC ? xxxidx[xidx] : k;
            int_t col = type == csx_t::CSC ? k : xxxidx[xidx];
            listPrinter.streamTuple(cnt++, row, col, values[xidx]);
		} // xidx

	} // k
}
/*-------------------------------------------------*/
template void print_to_stream(std::ostream&, csx_t, int_t, const int_t*, const int_t*, const real_t    *, std::streamsize);
template void print_to_stream(std::ostream&, csx_t, int_t, const int_t*, const int_t*, const real4_t   *, std::streamsize);
template void print_to_stream(std::ostream&, csx_t, int_t, const int_t*, const int_t*, const complex_t *, std::streamsize);
template void print_to_stream(std::ostream&, csx_t, int_t, const int_t*, const int_t*, const complex8_t*, std::streamsize);
/*-------------------------------------------------*/
static void transpose_xxptr(int_t inp, const int_t *ixptr, const int_t *ixidx, 
                            int_t onp, int_t *oxptr)
{
    std::fill(oxptr, oxptr + (onp+1), 0);

	for(int_t k = 0; k < inp; k++) {

		int_t xbgn = ixptr[k];
		int_t xend = ixptr[k+1];

		for(int_t xidx = xbgn; xidx < xend; xidx++) {
			oxptr[ixidx[xidx] + 1]++;
		} // xidx

	} // k

	roll(onp, oxptr);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
static void hybrid_transpose_tmpl(int_t inp, const int_t *ixptr, const int_t *ixidx, const T_Scalar *ivals, 
		                          int_t onp, int_t *oxptr, int_t *oxidx, T_Scalar *ovals, 
                                  T_Scalar coeff, bool conjop) 
{
	transpose_xxptr(inp, ixptr, ixidx, onp, oxptr);

	for(int_t k = 0; k < inp; k++) {

		int_t xbgn = ixptr[k];
		int_t xend = ixptr[k+1];

		for(int_t xidx = xbgn; xidx < xend; xidx++) {
			oxidx[oxptr[ixidx[xidx]]] = k;
			ovals[oxptr[ixidx[xidx]]] = coeff * (conjop ? arith::conj(ivals[xidx]) : ivals[xidx]);
			oxptr[ixidx[xidx]]++;
		} // xidx

	} // k

	unroll(onp, oxptr);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void transpose(int_t inp, const int_t *ixptr, const int_t *ixidx, const T_Scalar *ivals,
		       int_t onp, int_t *oxptr, int_t *oxidx, T_Scalar *ovals, T_Scalar coeff)
{
	hybrid_transpose_tmpl(inp, ixptr, ixidx, ivals, 
                          onp, oxptr, oxidx, ovals, 
                          coeff, false);
}
/*-------------------------------------------------*/
template void transpose(int_t, const int_t*, const int_t*, const real_t    *, int_t, int_t*, int_t*, real_t    *, real_t    );
template void transpose(int_t, const int_t*, const int_t*, const real4_t   *, int_t, int_t*, int_t*, real4_t   *, real4_t   );
template void transpose(int_t, const int_t*, const int_t*, const complex_t *, int_t, int_t*, int_t*, complex_t *, complex_t );
template void transpose(int_t, const int_t*, const int_t*, const complex8_t*, int_t, int_t*, int_t*, complex8_t*, complex8_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
void conjugate_transpose(int_t inp, const int_t *ixptr, const int_t *ixidx, const T_Scalar *ivals,
		                 int_t onp, int_t *oxptr, int_t *oxidx, T_Scalar *ovals, T_Scalar coeff)
{
	hybrid_transpose_tmpl(inp, ixptr, ixidx, ivals, 
                          onp, oxptr, oxidx, ovals, 
                          coeff, true);
}
/*-------------------------------------------------*/
template void conjugate_transpose(int_t, const int_t*, const int_t*, const real_t    *, int_t, int_t*, int_t*, real_t    *, real_t    );
template void conjugate_transpose(int_t, const int_t*, const int_t*, const real4_t   *, int_t, int_t*, int_t*, real4_t   *, real4_t   );
template void conjugate_transpose(int_t, const int_t*, const int_t*, const complex_t *, int_t, int_t*, int_t*, complex_t *, complex_t );
template void conjugate_transpose(int_t, const int_t*, const int_t*, const complex8_t*, int_t, int_t*, int_t*, complex8_t*, complex8_t);
/*-------------------------------------------------*/
void uplo2ge_xxptr(uplo_t uplo, csx_t type, int_t np, const int_t *ixptr, const int_t *ixidx, int_t *oxptr)
{
	if(uplo == uplo_t::Full) {
		std::copy(ixptr, ixptr + (np+1), oxptr);
		return;
	}

	std::fill(oxptr, oxptr + (np+1), 0);

	for(int_t k = 0; k < np; k++) {

		int_t xbgn = ixptr[k];
		int_t xend = ixptr[k+1];

		for(int_t xidx = xbgn; xidx < xend; xidx++) {

            int_t row = (type == csx_t::CSC ? ixidx[xidx] : k);
            int_t col = (type == csx_t::CSC ? k : ixidx[xidx]);

			if(row == col) {
				oxptr[col+1]++;
			} else if(uplo == uplo_t::Lower && row > col) {
				oxptr[col+1]++;
				oxptr[row+1]++;
			} else if(uplo == uplo_t::Upper && row < col) {
				oxptr[col+1]++;
				oxptr[row+1]++;
			} // lower/upper

		} // xidx

	} // k

	roll(np, oxptr);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
static inline void append_to_selected_part(int_t l, int_t k, T_Scalar v, int_t *oxptr, int_t *oxidx, T_Scalar *ovals, bool conjop)
{
	// original part
	oxidx[oxptr[k]] = l;
	ovals[oxptr[k]] = v;
	oxptr[k]++;

	// opposite part
	if(l != k) {
		oxidx[oxptr[l]] = k;
		ovals[oxptr[l]] = (conjop ? arith::conj(v) : v);
		oxptr[l]++;
	} // strict part
}
/*-------------------------------------------------*/
template <typename T_Scalar>
static void xx2ge(uplo_t uplo, csx_t type, int_t np, 
                  const int_t *ixptr, const int_t *ixidx, const T_Scalar *ivals, 
		          int_t *oxptr, int_t *oxidx, T_Scalar *ovals, bool conjop) 
{
	uplo2ge_xxptr(uplo, type, np, ixptr, ixidx, oxptr);
    
	for(int_t k = 0; k < np; k++) {

		int_t xbgn = ixptr[k];
		int_t xend = ixptr[k+1];

		for(int_t xidx = xbgn; xidx < xend; xidx++) {

            int_t l = ixidx[xidx];
            int_t row = (type == csx_t::CSC ? l : k);
            int_t col = (type == csx_t::CSC ? k : l);
            T_Scalar val = ivals[xidx];

			if(uplo == uplo_t::Lower && row >= col) {
				//append_to_selected_part(type == csx_t::CSC ? row : col, 
                //                        type == csx_t::CSC ? col : row, 
                //                        val, oxptr, oxidx, ovals, conjop);
                append_to_selected_part(l, k, val, oxptr, oxidx, ovals, conjop);
			} else if(uplo == uplo_t::Upper && row <= col) {
				//append_to_selected_part(type == csx_t::CSC ? row : col, 
                //                        type == csx_t::CSC ? col : row, 
                //                        val, oxptr, oxidx, ovals, conjop);
                append_to_selected_part(l, k, val, oxptr, oxidx, ovals, conjop);
			} // lower/upper

		} // xidx

	} // k

	unroll(np, oxptr);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void sy2ge(uplo_t uplo, csx_t type, int_t np, 
           const int_t *ixptr, const int_t *ixidx, const T_Scalar *ivals,
           int_t *oxptr, int_t *oxidx, T_Scalar *ovals)
{
	xx2ge(uplo, type, np, ixptr, ixidx, ivals, oxptr, oxidx, ovals, false);
}
/*-------------------------------------------------*/
template void sy2ge(uplo_t, csx_t, int_t, const int_t*, const int_t*, const real_t    *, int_t*, int_t*, real_t    *);
template void sy2ge(uplo_t, csx_t, int_t, const int_t*, const int_t*, const real4_t   *, int_t*, int_t*, real4_t   *);
template void sy2ge(uplo_t, csx_t, int_t, const int_t*, const int_t*, const complex_t *, int_t*, int_t*, complex_t *);
template void sy2ge(uplo_t, csx_t, int_t, const int_t*, const int_t*, const complex8_t*, int_t*, int_t*, complex8_t*);
/*-------------------------------------------------*/
template <typename T_Scalar>
void he2ge(uplo_t uplo, csx_t type, int_t np, 
           const int_t *ixptr, const int_t *ixidx, const T_Scalar *ivals,
           int_t *oxptr, int_t *oxidx, T_Scalar *ovals)
{
	xx2ge(uplo, type, np, ixptr, ixidx, ivals, oxptr, oxidx, ovals, true);
}
/*-------------------------------------------------*/
template void he2ge(uplo_t, csx_t, int_t, const int_t*, const int_t*, const real_t    *, int_t*, int_t*, real_t    *);
template void he2ge(uplo_t, csx_t, int_t, const int_t*, const int_t*, const real4_t   *, int_t*, int_t*, real4_t   *);
template void he2ge(uplo_t, csx_t, int_t, const int_t*, const int_t*, const complex_t *, int_t*, int_t*, complex_t *);
template void he2ge(uplo_t, csx_t, int_t, const int_t*, const int_t*, const complex8_t*, int_t*, int_t*, complex8_t*);
/*-------------------------------------------------*/
template <typename T_Scalar>
static inline void apply_op(T_Scalar& v, const T_Scalar& u, dup_t op)
{
	switch(op) {
		case dup_t::Sum:  v += u; break;
		case dup_t::Prod: v *= u; break;
		case dup_t::Amax: v = std::max(std::abs(v), std::abs(u)); break;
		case dup_t::Amin: v = std::min(std::abs(v), std::abs(u)); break;
        default:
            throw err::Exception("Unknown duplicate handling operation");
	}
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void remove_duplicates(int_t np, int_t *xxxptr, int_t *xxxidx, T_Scalar *values, dup_t op)
{
	for(int_t k = 0; k < np; k++) {

		int_t xbgn = xxxptr[k];
		int_t xend = xxxptr[k+1];
		int_t xlen = xend - xbgn;

		xxxptr[k] = (k > 0 ? xxxptr[k-1] : 0);

		if(!xlen) continue;

		int_t xref = xxxidx[xbgn];

		xxxidx[xxxptr[k]] = xxxidx[xbgn];
		values[xxxptr[k]] = values[xbgn];
		xxxptr[k]++;

		for(int_t xidx = xbgn + 1; xidx < xend; xidx++) {

			int_t    l = xxxidx[xidx];
			T_Scalar v = values[xidx];

			if(l == xref) {
				apply_op<T_Scalar>(values[xxxptr[k]-1], v, op);
			} else {
				xxxidx[xxxptr[k]] = l;
				values[xxxptr[k]] = v;
				xxxptr[k]++;
				xref = l;
			} // dup check

		} // xidx

	} // k

	xxxptr[np] = xxxptr[np-1];

	unroll(np, xxxptr);
}
/*-------------------------------------------------*/
template void remove_duplicates(int_t, int_t*, int_t*, real_t    *, dup_t);
template void remove_duplicates(int_t, int_t*, int_t*, real4_t   *, dup_t);
template void remove_duplicates(int_t, int_t*, int_t*, complex_t *, dup_t);
template void remove_duplicates(int_t, int_t*, int_t*, complex8_t*, dup_t);
/*-------------------------------------------------*/
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type
idx_max_asum(prop_t ptype, int_t np, const int_t *xxxptr, const int_t *xxxidx, const T_Scalar *values)
{
	using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

  if(np <= 0) return 0;

	std::vector<T_RScalar> idx_norms(np, 0);

	for(int_t k = 0; k < np; k++) {

		int_t xbgn = xxxptr[k];
		int_t xend = xxxptr[k+1];

		for(int_t xidx = xbgn; xidx < xend; xidx++) {

			int_t l = xxxidx[xidx];
			T_RScalar av = std::abs(values[xidx]);

			idx_norms[k] += av;

			if(ptype == prop_t::Symmetric || ptype == prop_t::Hermitian || ptype == prop_t::Skew) {
				if(l != k) {
					idx_norms[l] += av;
				} // off diag
			} // ptype

		} // xidx

	} // k

	return idx_norms[blas::iamax(np, idx_norms.data(), 1)];
}
/*-------------------------------------------------*/
template real_t  idx_max_asum(prop_t, int_t, const int_t*, const int_t*, const real_t    *);
template real4_t idx_max_asum(prop_t, int_t, const int_t*, const int_t*, const real4_t   *);
template real_t  idx_max_asum(prop_t, int_t, const int_t*, const int_t*, const complex_t *);
template real4_t idx_max_asum(prop_t, int_t, const int_t*, const int_t*, const complex8_t*);
/*-------------------------------------------------*/
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type
ptr_max_asum(prop_t ptype, int_t np, int_t maxIdx, const int_t *xxxptr, const int_t *xxxidx, const T_Scalar *values)
{
	using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

  if(np <= 0 || maxIdx <= 0) return 0;

	std::vector<T_RScalar> ptr_norms(maxIdx, 0);

	for(int_t k = 0; k < np; k++) {

		int_t xbgn = xxxptr[k];
		int_t xend = xxxptr[k+1];

		for(int_t xidx = xbgn; xidx < xend; xidx++) {

			int_t l = xxxidx[xidx];
			T_RScalar av = std::abs(values[xidx]);

			ptr_norms[l] += av;

			if(ptype == prop_t::Symmetric || ptype == prop_t::Hermitian || ptype == prop_t::Skew) {
				if(l != k) {
					ptr_norms[k] += av;
				} // off diag
			} // ptype

		} // xidx

	} // k

	return ptr_norms[blas::iamax(maxIdx, ptr_norms.data(), 1)];
}
/*-------------------------------------------------*/
template real_t  ptr_max_asum(prop_t, int_t, int_t, const int_t*, const int_t*, const real_t    *);
template real4_t ptr_max_asum(prop_t, int_t, int_t, const int_t*, const int_t*, const real4_t   *);
template real_t  ptr_max_asum(prop_t, int_t, int_t, const int_t*, const int_t*, const complex_t *);
template real4_t ptr_max_asum(prop_t, int_t, int_t, const int_t*, const int_t*, const complex8_t*);
/*-------------------------------------------------*/
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type
norm_max(int_t np, const int_t *xxxptr, const T_Scalar *values)
{
    if(np <= 0) return 0;

	T_Scalar vamax = values[blas::iamax(xxxptr[np], values, 1)];
	return std::abs(vamax);
}
/*-------------------------------------------------*/
template real_t  norm_max(int_t, const int_t*, const real_t    *);
template real4_t norm_max(int_t, const int_t*, const real4_t   *);
template real_t  norm_max(int_t, const int_t*, const complex_t *);
template real4_t norm_max(int_t, const int_t*, const complex8_t*);
/*-------------------------------------------------*/
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type
norm_fro(prop_t ptype, int_t np, const int_t *xxxptr, const int_t *xxxidx, const T_Scalar *values)
{
	using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

  if(np <= 0) return 0;

	T_RScalar ret = 0;

	for(int_t k = 0; k < np; k++) {

		int_t xbgn = xxxptr[k];
		int_t xend = xxxptr[k+1];

		for(int_t xidx = xbgn; xidx < xend; xidx++) {

			T_RScalar av = std::abs(values[xidx]);
			T_RScalar av2 = av * av;

			if(ptype == prop_t::General || ptype == prop_t::Triangular) {

				ret += av2;

			} else if(ptype == prop_t::Symmetric || ptype == prop_t::Hermitian || ptype == prop_t::Skew) {

				ret += (xxxidx[xidx] == k ? av2 : T_RScalar(2) * av2);

			} // ptype

		} // xidx

	} // k

	return std::sqrt(ret);
}
/*-------------------------------------------------*/
template real_t  norm_fro(prop_t, int_t, const int_t*, const int_t*, const real_t    *);
template real4_t norm_fro(prop_t, int_t, const int_t*, const int_t*, const real4_t   *);
template real_t  norm_fro(prop_t, int_t, const int_t*, const int_t*, const complex_t *);
template real4_t norm_fro(prop_t, int_t, const int_t*, const int_t*, const complex8_t*);
/*-------------------------------------------------*/
} // namespace csx
} // namespace blk
} // namespace cla3p
/*-------------------------------------------------*/
