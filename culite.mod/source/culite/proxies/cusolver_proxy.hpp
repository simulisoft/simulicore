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

#ifndef CULITE_CUSOLVER_PROXY_HPP_
#define CULITE_CUSOLVER_PROXY_HPP_

/**
 * @file
 */

#include <cusolverDn.h>

/*-------------------------------------------------*/
namespace culite {
namespace cusolver {
/*-------------------------------------------------*/

cusolverEigMode_t bool2cusolverEigMode(bool job);

/*-------------------------------------------------*/

class DnParams {

    public:
        DnParams();
        ~DnParams();

        cusolverDnParams_t get() const { return m_params; }

    private:
        cusolverDnParams_t m_params{nullptr};
};

/*-------------------------------------------------*/
} // namespace cusolver
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_CUSOLVER_PROXY_HPP_
