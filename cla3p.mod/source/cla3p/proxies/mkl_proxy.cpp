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
#include "cla3p/proxies/mkl_proxy.hpp"

// system

// 3rd
#include <mkl.h>

// cla3p

/*-------------------------------------------------*/
namespace cla3p {
namespace mkl {
/*-------------------------------------------------*/
std::string version()
{
	const nint_t len = 256;
	char buffer[len];
	mkl_get_version_string(buffer, len);
	return buffer;
}
/*-------------------------------------------------*/
#define MKL_ALLOC_ALIGNMENT 64
/*-------------------------------------------------*/
void* malloc(std::size_t size)
{
	return mkl_malloc(size, MKL_ALLOC_ALIGNMENT);
}
/*-------------------------------------------------*/
void* calloc(std::size_t size, std::size_t nmemb)
{
	return mkl_calloc(size, nmemb, MKL_ALLOC_ALIGNMENT);
}
/*-------------------------------------------------*/
#undef MKL_ALLOC_ALIGNMENT
/*-------------------------------------------------*/
void* realloc(void *ptr, std::size_t nmemb)
{
	return mkl_realloc(ptr, nmemb);
}
/*-------------------------------------------------*/
void free(void *ptr)
{
	mkl_free(ptr);
}
/*-------------------------------------------------*/
nint_t get_dynamic()
{
	return mkl_get_dynamic();
}
/*-------------------------------------------------*/
nint_t get_max_threads()
{
	return mkl_get_max_threads();
}
/*-------------------------------------------------*/
void set_num_threads(nint_t numThreads)
{
	mkl_set_num_threads(numThreads);
}
/*-------------------------------------------------*/
nint_t set_num_threads_local(nint_t numThreads)
{
	return mkl_set_num_threads_local(numThreads);
}
/*-------------------------------------------------*/
#define omatcopy_macro(typeout, typein, prefix) \
typeout omatcopy(char ordering, char trans, std::size_t rows, std::size_t cols, typein alpha, \
		const typein *a, std::size_t lda, \
		      typein *b, std::size_t ldb) \
{ \
	mkl_##prefix##omatcopy(ordering, trans, rows, cols, alpha, a, lda, b, ldb); \
}
omatcopy_macro(void, real_t    , d)
omatcopy_macro(void, real4_t   , s)
omatcopy_macro(void, complex_t , z)
omatcopy_macro(void, complex8_t, c)
#undef omatcopy_macro
/*-------------------------------------------------*/
#define omatcopy_macro(typeout, typein, prefix) \
typeout omatcopy(char ordering, char trans, std::size_t rows, std::size_t cols, typein alpha, \
		const typein *a, std::size_t lda, std::size_t stridea, \
		      typein *b, std::size_t ldb, std::size_t strideb) \
{ \
	mkl_##prefix##omatcopy2(ordering, trans, rows, cols, alpha, a, lda, stridea, b, ldb, strideb); \
}
omatcopy_macro(void, real_t    , d)
omatcopy_macro(void, real4_t   , s)
omatcopy_macro(void, complex_t , z)
omatcopy_macro(void, complex8_t, c)
#undef omatcopy_macro
/*-------------------------------------------------*/
#define imatcopy_macro(typeout, typein, prefix) \
typeout imatcopy(char ordering, char trans, std::size_t rows, std::size_t cols, typein alpha, typein *a, std::size_t lda, std::size_t ldb) \
{ \
	mkl_##prefix##imatcopy(ordering, trans, rows, cols, alpha, a, lda, ldb); \
}
imatcopy_macro(void, real_t    , d)
imatcopy_macro(void, real4_t   , s)
imatcopy_macro(void, complex_t , z)
imatcopy_macro(void, complex8_t, c)
#undef imatcopy_macro
/*-------------------------------------------------*/
#define omatadd_macro(typeout, typein, prefix) \
typeout omatadd(char ordering, char transa, char transb, std::size_t rows, std::size_t cols, \
    typein alpha, const typein *A, std::size_t lda, \
    typein beta , const typein *B, std::size_t ldb, typein *C, std::size_t ldc) \
{ \
	mkl_##prefix##omatadd(ordering, transa, transb, rows, cols, alpha, A, lda, beta , B, ldb, C, ldc); \
}
omatadd_macro(void, real_t    , d)
omatadd_macro(void, real4_t   , s)
omatadd_macro(void, complex_t , z)
omatadd_macro(void, complex8_t, c)
#undef omatadd_macro
/*-------------------------------------------------*/
} // namespace mkl
} // namespace cla3p
/*-------------------------------------------------*/
