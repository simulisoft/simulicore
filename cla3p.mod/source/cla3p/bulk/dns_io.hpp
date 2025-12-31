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

#ifndef CLA3P_BULK_DNS_IO_HPP_
#define CLA3P_BULK_DNS_IO_HPP_

/**
 * @file
 */

#include <string>
#include <iostream>
#include <sstream>

#include "cla3p/types.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace blk {
namespace dns {
/*-------------------------------------------------*/

void print_to_stream(std::ostream& os, uplo_t uplo, int_t m, int_t n, const int_t      *a, int_t lda, std::streamsize prec = 0, int_t cols_per_page = 10);
void print_to_stream(std::ostream& os, uplo_t uplo, int_t m, int_t n, const real_t     *a, int_t lda, std::streamsize prec = 0, int_t cols_per_page = 10);
void print_to_stream(std::ostream& os, uplo_t uplo, int_t m, int_t n, const real4_t    *a, int_t lda, std::streamsize prec = 0, int_t cols_per_page = 10);
void print_to_stream(std::ostream& os, uplo_t uplo, int_t m, int_t n, const complex_t  *a, int_t lda, std::streamsize prec = 0, int_t cols_per_page = 5 );
void print_to_stream(std::ostream& os, uplo_t uplo, int_t m, int_t n, const complex8_t *a, int_t lda, std::streamsize prec = 0, int_t cols_per_page = 5 );

/*-------------------------------------------------*/
} // namespace dns
} // namespace blk
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_BULK_DNS_IO_HPP_
