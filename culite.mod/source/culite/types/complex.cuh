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

#ifndef CULITE_COMPLEX_CUH_
#define CULITE_COMPLEX_CUH_

/** 
 * @file
 */

 #include <iostream>
 #include <cuComplex.h>

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/

class DoubleComplex {

	private:
		using T_RScalar = double;

	public:
		__host__ __device__ DoubleComplex(T_RScalar re = 0.0, T_RScalar im = 0.0) 
		: m_data(cuDoubleComplex{re, im}) {}
		~DoubleComplex() = default;

    __host__ __device__ T_RScalar& real() { return m_data.x; }
    __host__ __device__ T_RScalar& imag() { return m_data.y; }

    __host__ __device__ const T_RScalar& real() const { return m_data.x; }
    __host__ __device__ const T_RScalar& imag() const { return m_data.y; }

	__host__ __device__ DoubleComplex operator+(const DoubleComplex& c) const {
        cuDoubleComplex res = cuCadd(m_data, c.m_data);
        return DoubleComplex(cuCreal(res), cuCimag(res));
    }

	__host__ __device__ DoubleComplex operator-(const DoubleComplex& c) const {
        cuDoubleComplex res = cuCsub(m_data, c.m_data);
        return DoubleComplex(cuCreal(res), cuCimag(res));
    }

	__host__ __device__ DoubleComplex operator*(const DoubleComplex& c) const {
        cuDoubleComplex res = cuCmul(m_data, c.m_data);
        return DoubleComplex(cuCreal(res), cuCimag(res));
    }

	__host__ __device__ DoubleComplex operator/(const DoubleComplex& c) const {
        cuDoubleComplex res = cuCdiv(m_data, c.m_data);
        return DoubleComplex(cuCreal(res), cuCimag(res));
    }

	cuDoubleComplex* cudaPtr() { return &m_data; }
	const cuDoubleComplex* cudaPtr() const { return &m_data; }

	private:
		cuDoubleComplex m_data;

};

/*-------------------------------------------------*/

class SingleComplex {

	private:
		using T_RScalar = float;

	public:
		__host__ __device__ SingleComplex(T_RScalar re = 0.0f, T_RScalar im = 0.0f) 
		: m_data(cuFloatComplex{re, im}) {}
		~SingleComplex() = default;

	__host__ __device__ T_RScalar& real() { return m_data.x; }
	__host__ __device__ T_RScalar& imag() { return m_data.y; }

	__host__ __device__ const T_RScalar& real() const { return m_data.x; }
	__host__ __device__ const T_RScalar& imag() const { return m_data.y; }
	
	__host__ __device__ SingleComplex operator+(const SingleComplex& c) const {
		cuFloatComplex res = cuCaddf(m_data, c.m_data);
		return SingleComplex(cuCrealf(res), cuCimagf(res));
	}

	__host__ __device__ SingleComplex operator-(const SingleComplex& c) const {
		cuFloatComplex res = cuCsubf(m_data, c.m_data);
		return SingleComplex(cuCrealf(res), cuCimagf(res));
	}

	__host__ __device__ SingleComplex operator*(const SingleComplex& c) const {
		cuFloatComplex res = cuCmulf(m_data, c.m_data);
		return SingleComplex(cuCrealf(res), cuCimagf(res));
	}

	__host__ __device__ SingleComplex operator/(const SingleComplex& c) const {
		cuFloatComplex res = cuCdivf(m_data, c.m_data);
		return SingleComplex(cuCrealf(res), cuCimagf(res));
	}

	cuFloatComplex* cudaPtr() { return &m_data; }
	const cuFloatComplex* cudaPtr() const { return &m_data; }

	private:
		cuFloatComplex m_data;

};

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

__host__ std::ostream& operator<<(std::ostream& os, const culite::DoubleComplex& c) {
    os << "(" << c.real() << "," << c.imag() << ")";
    return os;
}

__host__ std::ostream& operator<<(std::ostream& os, const culite::SingleComplex& c) {
    os << "(" << c.real() << "," << c.imag() << ")";
    return os;
}

/*-------------------------------------------------*/

#endif // CULITE_COMPLEX_CUH_
