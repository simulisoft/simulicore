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

#ifndef CULITE_COMPLEX_HPP_
#define CULITE_COMPLEX_HPP_

/** 
 * @file
 */

 #include <iostream>
 #include <cuComplex.h>
 #include <culite/types/cuda_api.hpp>

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/

class DoubleComplex {

    private:
        using T_RScalar = double;

    public:
        CULITE_HOST CULITE_DEVICE DoubleComplex(T_RScalar re = 0.0, T_RScalar im = 0.0) 
        : m_data(cuDoubleComplex{re, im}) {}
        ~DoubleComplex() = default;

        // host only
        CULITE_HOST T_RScalar& rreal() { return m_data.x; }
        CULITE_HOST T_RScalar& rimag() { return m_data.y; }
        CULITE_HOST const T_RScalar& rreal() const { return m_data.x; }
        CULITE_HOST const T_RScalar& rimag() const { return m_data.y; }

        // generic
        CULITE_DEVICE CULITE_DEVICE T_RScalar real() const { return m_data.x; }
        CULITE_DEVICE CULITE_DEVICE T_RScalar imag() const { return m_data.y; }

        CULITE_HOST CULITE_DEVICE void setReal(const T_RScalar r) { m_data.x = r; }
        CULITE_HOST CULITE_DEVICE void setImag(const T_RScalar i) { m_data.y = i; }

        CULITE_HOST CULITE_DEVICE T_RScalar abs() const { return cuCabs(m_data); }
        CULITE_HOST CULITE_DEVICE T_RScalar abs2() const 
        { 
            return cuCreal(m_data)*cuCreal(m_data) + cuCimag(m_data)*cuCimag(m_data); 
        }

        CULITE_HOST CULITE_DEVICE DoubleComplex operator+(const DoubleComplex c) const {
            cuDoubleComplex res = cuCadd(m_data, c.m_data);
            return DoubleComplex(cuCreal(res), cuCimag(res));
        }

        CULITE_HOST CULITE_DEVICE DoubleComplex operator-(const DoubleComplex c) const {
            cuDoubleComplex res = cuCsub(m_data, c.m_data);
            return DoubleComplex(cuCreal(res), cuCimag(res));
        }

        CULITE_HOST CULITE_DEVICE DoubleComplex operator*(const DoubleComplex c) const {
            cuDoubleComplex res = cuCmul(m_data, c.m_data);
            return DoubleComplex(cuCreal(res), cuCimag(res));
        }

        CULITE_HOST CULITE_DEVICE DoubleComplex operator/(const DoubleComplex c) const {
            cuDoubleComplex res = cuCdiv(m_data, c.m_data);
            return DoubleComplex(cuCreal(res), cuCimag(res));
        }

    private:
        cuDoubleComplex m_data;

};

/*-------------------------------------------------*/

class SingleComplex {

    private:
        using T_RScalar = float;

    public:
        CULITE_HOST CULITE_DEVICE SingleComplex(T_RScalar re = 0.0f, T_RScalar im = 0.0f) 
        : m_data(cuFloatComplex{re, im}) {}
        ~SingleComplex() = default;

        // host only
        CULITE_HOST T_RScalar& rreal() { return m_data.x; }
        CULITE_HOST T_RScalar& rimag() { return m_data.y; }
        CULITE_HOST const T_RScalar& rreal() const { return m_data.x; }
        CULITE_HOST const T_RScalar& rimag() const { return m_data.y; }

        // generic
        CULITE_HOST CULITE_DEVICE T_RScalar real() const { return m_data.x; }
        CULITE_HOST CULITE_DEVICE T_RScalar imag() const { return m_data.y; }

        CULITE_HOST CULITE_DEVICE void setReal(const T_RScalar r) { m_data.x = r; }
        CULITE_HOST CULITE_DEVICE void setImag(const T_RScalar i) { m_data.y = i; }

        CULITE_HOST CULITE_DEVICE T_RScalar abs() const { return cuCabsf(m_data); }
        CULITE_HOST CULITE_DEVICE T_RScalar abs2() const 
        { 
            return cuCrealf(m_data)*cuCrealf(m_data) + cuCimagf(m_data)*cuCimagf(m_data); 
        }

        CULITE_HOST CULITE_DEVICE SingleComplex operator+(const SingleComplex c) const {
            cuFloatComplex res = cuCaddf(m_data, c.m_data);
            return SingleComplex(cuCrealf(res), cuCimagf(res));
        }

        CULITE_HOST CULITE_DEVICE SingleComplex operator-(const SingleComplex c) const {
            cuFloatComplex res = cuCsubf(m_data, c.m_data);
            return SingleComplex(cuCrealf(res), cuCimagf(res));
        }

        CULITE_HOST CULITE_DEVICE SingleComplex operator*(const SingleComplex c) const {
            cuFloatComplex res = cuCmulf(m_data, c.m_data);
            return SingleComplex(cuCrealf(res), cuCimagf(res));
        }

        CULITE_HOST CULITE_DEVICE SingleComplex operator/(const SingleComplex c) const {
            cuFloatComplex res = cuCdivf(m_data, c.m_data);
            return SingleComplex(cuCrealf(res), cuCimagf(res));
        }

    private:
        cuFloatComplex m_data;

};

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

CULITE_HOST inline std::ostream& operator<<(std::ostream& os, const culite::DoubleComplex& c) {
    os << "(" << c.rreal() << "," << c.rimag() << ")";
    return os;
}

CULITE_HOST inline std::ostream& operator<<(std::ostream& os, const culite::SingleComplex& c) {
    os << "(" << c.rreal() << "," << c.rimag() << ")";
    return os;
}

/*-------------------------------------------------*/

#endif // CULITE_COMPLEX_HPP_
