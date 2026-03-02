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

#ifndef CULITE_GRID_HPP_
#define CULITE_GRID_HPP_

/** 
 * @file
 */

#include <cstddef>
#include <cuda_runtime.h>

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief Base class for CUDA kernel grid configuration.
 * @details GridXD provides the foundation for managing CUDA grid and block dimensions.
 * It includes validation and overflow checking to ensure safe grid configurations.
 * This is a base class for specialized 1D and 2D grid configurations.
 */
class GridXD {

    protected:
        using cudaInt = unsigned int;  ///< Type used by CUDA dim3 dimensions

    public:
        /**
         * @brief Default constructor.
         * @details Creates a GridXD with default-initialized dimensions.
         */
        GridXD();
        
        /**
         * @brief Constructs a GridXD with explicit grid and block dimensions.
         * @param numBlocks Number of blocks in the grid.
         * @param threadsPerBlock Number of threads per block.
         */
        GridXD(dim3 numBlocks, dim3 threadsPerBlock);
        
        /**
         * @brief Default destructor.
         */
        ~GridXD();

        /**
         * @brief Gets the number of blocks in the grid.
         * @return Reference to the dim3 structure containing block counts.
         */
        const dim3& numBlocks() const;
        
        /**
         * @brief Gets the number of threads per block.
         * @return Reference to the dim3 structure containing threads per block.
         */
        const dim3& threadsPerBlock() const;

    protected:
        dim3 m_numBlocks;
        dim3 m_threadsPerBlock;

        void input_check(std::size_t n, std::size_t blockSize) const;
        void input_check(std::size_t m, std::size_t n, std::size_t blockSizeX, std::size_t blockSizeY) const;
        void overflow_check(std::size_t val) const;
};

/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief 1D CUDA kernel grid configuration.
 * @details Grid1D simplifies configuration for 1D kernels by automatically calculating
 * the number of blocks needed to cover n elements with the specified block size.
 */
class Grid1D : public GridXD {

    public:
        /**
         * @brief Default constructor.
         * @details Creates a Grid1D with default-initialized dimensions.
         */
        Grid1D();
        
        /**
         * @brief Constructs a 1D grid configuration.
         * @details Automatically calculates the number of blocks needed to cover n elements
         * with the given block size. Performs overflow and validation checks.
         * @param n Total number of elements to process.
         * @param blockSize Number of threads per block (default: 256).
         * @throws err::CudaException if blockSize is zero, if arithmetic overflow occurs,
         * or if calculated values exceed unsigned int maximum.
         */
        Grid1D(std::size_t n, std::size_t blockSize = 256);
        
        /**
         * @brief Default destructor.
         */
        ~Grid1D();
};

/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief 2D CUDA kernel grid configuration.
 * @details Grid2D simplifies configuration for 2D kernels (e.g., matrix operations) by
 * automatically calculating the number of blocks needed to cover m×n elements with the
 * specified block dimensions. Optimized for column-major matrices.
 */
class Grid2D : public GridXD {

    public:
        /**
         * @brief Default constructor.
         * @details Creates a Grid2D with default-initialized dimensions.
         */
        Grid2D();
        
        /**
         * @brief Constructs a 2D grid configuration with explicit block dimensions.
         * @details Automatically calculates the number of blocks needed in X and Y dimensions
         * to cover m×n elements with the given block sizes. Performs overflow and validation checks.
         * @param m Number of elements in the X dimension (rows).
         * @param n Number of elements in the Y dimension (columns).
         * @param blockSizeX Number of threads per block in X dimension (default: 16).
         * @param blockSizeY Number of threads per block in Y dimension (default: 16).
         * @throws err::CudaException if block sizes are zero, if arithmetic overflow occurs,
         * or if calculated values exceed unsigned int maximum.
         */
        Grid2D(std::size_t m, std::size_t n, std::size_t blockSizeX = 16, std::size_t blockSizeY = 16);
        
        /**
         * @brief Default destructor.
         */
        ~Grid2D();
};

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_GRID_HPP_
