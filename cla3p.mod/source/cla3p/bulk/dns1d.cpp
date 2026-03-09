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
#include "cla3p/bulk/dns1d.hpp"

// system
#include <functional>

// 3rd

// cla3p
#include "cla3p/types/scalar.hpp"
#include "cla3p/proxies/blas_proxy.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace blk {
namespace dns {
/*-------------------------------------------------*/
template <typename T_Scalar>
void fill(int_t n, T_Scalar *x, T_Scalar val, int_t incx)
{
	if(!n) return;

    if(!incx) {

        *x = val; 

    } else if(incx == 1) {

        std::fill_n(x, n, val);

    } else {

        for(int_t i = 0; i < n; i++, x += incx) {
            *x = val;
        } // i

    } // incx
}
/*-------------------------------------------------*/
template void fill(int_t n, int_t      *x, int_t      val, int_t incx);
template void fill(int_t n, real_t     *x, real_t     val, int_t incx);
template void fill(int_t n, real4_t    *x, real4_t    val, int_t incx);
template void fill(int_t n, complex_t  *x, complex_t  val, int_t incx);
template void fill(int_t n, complex8_t *x, complex8_t val, int_t incx);
/*-------------------------------------------------*/
template <>
void copy<int_t>(int_t n, const int_t *x, int_t *y, int_t incx, int_t incy)
{
    if(!n) return;

    for(int_t i = 0; i < n; i++, x += incx, y += incy) {
        *y = *x;
    } // i
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void copy(int_t n, const T_Scalar *x, T_Scalar* y, int_t incx, int_t incy)
{
	if(!n) return;

    blas::copy(n, x, incx, y, incy);
}
/*-------------------------------------------------*/
template void copy(int_t n, const real_t     *x, real_t     *y, int_t incx, int_t incy);
template void copy(int_t n, const real4_t    *x, real4_t    *y, int_t incx, int_t incy);
template void copy(int_t n, const complex_t  *x, complex_t  *y, int_t incx, int_t incy);
template void copy(int_t n, const complex8_t *x, complex8_t *y, int_t incx, int_t incy);
/*-------------------------------------------------*/
} // namespace dns
} // namespace blk
} // namespace cla3p
/*-------------------------------------------------*/
