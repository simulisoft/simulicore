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
#include "cla3p/bulk/csr.hpp"

// system

// 3rd

// cla3p
#include "cla3p/error/exceptions.hpp"
#include "cla3p/checks/basic_checks.hpp"
#include "cla3p/support/utils.hpp"
#include "cla3p/support/heap_buffer.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace blk {
namespace csr {
/*-------------------------------------------------*/
template <typename T_Scalar>
static void permute_ge_both_invq(int_t m, int_t /*n*/,
        const int_t *irowptr, const int_t *icolidx, const T_Scalar *ivalues,
        int_t *orowptr, int_t *ocolidx, T_Scalar *ovalues, const int_t *P, const int_t *Q)
{
    orowptr[0] = 0;
    for(int_t i = 0; i < m; i++) {
        orowptr[P[i] + 1] = irowptr[i+1] - irowptr[i];
    } // i

    csx::roll(m, orowptr);

    for(int_t i = 0; i < m; i++) {
        int_t Pi = P[i];
        for(int_t jcol = irowptr[i]; jcol < irowptr[i+1]; jcol++) {
            ocolidx[orowptr[Pi]] = Q[icolidx[jcol]];
            ovalues[orowptr[Pi]] = ivalues[jcol];
            orowptr[Pi]++;
        } // jcol
    } // i

    csx::unroll(m, orowptr);

    csx::sort(m, orowptr, ocolidx, ovalues);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
static void permute_ge_both(int_t m, int_t n,
        const int_t *irowptr, const int_t *icolidx, const T_Scalar *ivalues,
        int_t *orowptr, int_t *ocolidx, T_Scalar *ovalues, const int_t *P, const int_t *Q)
{
    // TODO : do this without extra buffer if possible
    HeapBuffer<int_t> invQ(n);

    for(int_t j = 0; j < n; j++) {
        invQ.data()[Q[j]] = j;
    } // j

    permute_ge_both_invq(m, n,
        irowptr, icolidx, ivalues,
        orowptr, ocolidx, ovalues, P, invQ.data());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
static void permute_ge_left(int_t m, int_t /*n*/,
        const int_t *irowptr, const int_t *icolidx, const T_Scalar *ivalues,
        int_t *orowptr, int_t *ocolidx, T_Scalar *ovalues, const int_t *P)
{
    orowptr[0] = 0;
    for(int_t i = 0; i < m; i++) {
        orowptr[P[i] + 1] = irowptr[i + 1] - irowptr[i];
    } // i

    csx::roll(m, orowptr);

    for(int_t i = 0; i < m; i++) {
        std::copy(icolidx + irowptr[i], icolidx + irowptr[i + 1], ocolidx + orowptr[P[i]]);
        std::copy(ivalues + irowptr[i], ivalues + irowptr[i + 1], ovalues + orowptr[P[i]]);
    } // i
}
/*-------------------------------------------------*/
template <typename T_Scalar>
static void permute_ge_right_invq(int_t m, int_t /*n*/,
        const int_t *irowptr, const int_t *icolidx, const T_Scalar *ivalues,
        int_t *orowptr, int_t *ocolidx, T_Scalar *ovalues, const int_t *Q)
{
    std::copy(irowptr, irowptr + m + 1, orowptr);

    for(int_t i = 0; i < m; i++) {
        for(int_t jcol = irowptr[i]; jcol < irowptr[i+1]; jcol++) {
            ocolidx[orowptr[i]] = Q[icolidx[jcol]];
            ovalues[orowptr[i]] = ivalues[jcol];
            orowptr[i]++;
        } // jcol
    } // i
    
    csx::unroll(m, orowptr);

    csx::sort(m, orowptr, ocolidx, ovalues);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
static void permute_ge_right(int_t m, int_t n,
        const int_t *irowptr, const int_t *icolidx, const T_Scalar *ivalues,
        int_t *orowptr, int_t *ocolidx, T_Scalar *ovalues, const int_t *Q)
{
    // TODO : do this without extra buffer if possible
    HeapBuffer<int_t> invQ(n);

    for(int_t j = 0; j < n; j++) {
        invQ.data()[Q[j]] = j;
    } // j

    permute_ge_right_invq(m, n,
                          irowptr, icolidx, ivalues,
                          orowptr, ocolidx, ovalues, invQ.data());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void permute_xx_mirror(prop_t ptype, uplo_t uplo, int_t n,
        const int_t *irowptr, const int_t *icolidx, const T_Scalar *ivalues,
        int_t *orowptr, int_t *ocolidx, T_Scalar *ovalues, const int_t *P)
{
    int_t Pi;
    int_t Pj;

    for(int_t i = 0; i < n + 1; i++) 
        orowptr[i] = 0;

    for(int_t i = 0; i < n; i++) {
        for(int_t jcol = irowptr[i]; jcol < irowptr[i+1]; jcol++) {
            Pj = P[icolidx[jcol]];
            Pi = P[i];
            if(uplo == uplo_t::Upper && Pj < Pi) 
                orowptr[Pj+1]++;
            else if(uplo == uplo_t::Lower && Pj > Pi) 
                orowptr[Pj+1]++;
            else
                orowptr[Pi+1]++;
        } // jcol
    } // i

    csx::roll(n, orowptr);

    for(int_t i = 0; i < n; i++) {
        for(int_t jcol = irowptr[i]; jcol < irowptr[i+1]; jcol++) {
            Pj = P[icolidx[jcol]];
            Pi = P[i];
            if(uplo == uplo_t::Upper && Pj < Pi) {
                ocolidx[orowptr[Pj]] = Pi;
                ovalues[orowptr[Pj]] = opposite_element(ivalues[jcol],ptype);
                orowptr[Pj]++;
            } else if(uplo == uplo_t::Lower && Pj > Pi) {
                ocolidx[orowptr[Pj]] = Pi;
                ovalues[orowptr[Pj]] = opposite_element(ivalues[jcol],ptype);
                orowptr[Pj]++;
            } else {
                ocolidx[orowptr[Pi]] = Pj;
                ovalues[orowptr[Pi]] = ivalues[jcol];
                orowptr[Pi]++;
            }
        } // jcol
    } // i
    
    csx::unroll(n, orowptr);

    csx::sort(n, orowptr, ocolidx, ovalues);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void permute(prop_t ptype, uplo_t uplo, int_t m, int_t n,
        const int_t *irowptr, const int_t *icolidx, const T_Scalar *ivalues,
        int_t *orowptr, int_t *ocolidx, T_Scalar *ovalues, const int_t *P, const int_t *Q)
{
    if(!m || !n) return;

    Property prop(ptype, uplo);

    if(prop.isSquare()) {
        square_check(m, n);
    }

    if(prop.isGeneral()) {

        if(P && Q) {

            permute_ge_both(m, n, irowptr, icolidx, ivalues, orowptr, ocolidx, ovalues, P, Q);

        } else if(P && !Q) {

            permute_ge_left(m, n, irowptr, icolidx, ivalues, orowptr, ocolidx, ovalues, P);

        } else if(!P && Q) { 

            permute_ge_right(m, n, irowptr, icolidx, ivalues, orowptr, ocolidx, ovalues, Q);

        } else {

            std::copy(irowptr, irowptr + m + 1    , orowptr);
            std::copy(icolidx, icolidx + irowptr[m], ocolidx);
            std::copy(ivalues, ivalues + irowptr[m], ovalues);

        } // P/Q

    } else if(prop.isSymmetric() || prop.isHermitian() || prop.isSkew()) {

        if(P) {

            permute_xx_mirror(ptype, uplo, n, irowptr, icolidx, ivalues, orowptr, ocolidx, ovalues, P);

        } else {

            std::copy(irowptr, irowptr + m + 1     , orowptr);
            std::copy(icolidx, icolidx + irowptr[m], ocolidx);
            std::copy(ivalues, ivalues + irowptr[m], ovalues);

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
} // namespace csr
} // namespace blk
} // namespace cla3p
/*-------------------------------------------------*/
