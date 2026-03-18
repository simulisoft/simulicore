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
#include <cla3p/eigsol/lapack_syev.hpp>

// system

// 3rd

// cla3p
#include <cla3p/dense.hpp>

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/
template <typename T_Matrix>
LapackSyev<T_Matrix>::LapackSyev()
{
    defaults();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
LapackSyev<T_Matrix>::~LapackSyev()
{
    clear();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackSyev<T_Matrix>::defaults()
{
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void LapackSyev<T_Matrix>::clear()
{
    defaults();
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template class LapackSyev<dns::RdMatrix>;
template class LapackSyev<dns::RfMatrix>;
template class LapackSyev<dns::CdMatrix>;
template class LapackSyev<dns::CfMatrix>;
/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/
