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
#include <culite/support/grid.hpp>

// system
#include <limits>

// 3rd

// culite
#include <culite/error/exceptions.hpp>

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/
GridXD::GridXD() {}
GridXD::~GridXD() {}
/*-------------------------------------------------*/
GridXD::GridXD(dim3 numBlocks, dim3 threadsPerBlock)
: m_numBlocks(numBlocks), 
  m_threadsPerBlock(threadsPerBlock)
{
}
/*-------------------------------------------------*/
const dim3& GridXD::numBlocks() const { return m_numBlocks; }
const dim3& GridXD::threadsPerBlock() const { return m_threadsPerBlock; }
/*-------------------------------------------------*/
void GridXD::input_check(std::size_t n, std::size_t blockSize) const
{
    if (blockSize == 0) {
        throw err::CudaException("GridXD: blockSize cannot be zero");
    }
    if (n > std::numeric_limits<std::size_t>::max() - blockSize + 1) {
        throw err::CudaException("GridXD: arithmetic overflow in grid dimension calculation");
    }
}
/*-------------------------------------------------*/
void GridXD::input_check(std::size_t m, std::size_t n, std::size_t blockSizeX, std::size_t blockSizeY) const
{
    if (blockSizeX == 0 || blockSizeY == 0) {
        throw err::CudaException("GridXD: block sizes cannot be zero");
    }
    if (m > std::numeric_limits<std::size_t>::max() - blockSizeX + 1) {
        throw err::CudaException("GridXD: arithmetic overflow in X grid dimension calculation");
    }
    if (n > std::numeric_limits<std::size_t>::max() - blockSizeY + 1) {
        throw err::CudaException("GridXD: arithmetic overflow in Y grid dimension calculation");
    }
}
/*-------------------------------------------------*/
void GridXD::overflow_check(std::size_t val) const
{
    if (val > std::numeric_limits<cudaInt>::max()) {
        throw err::CudaException("GridXD: value in grid dimension calculation exceeds unsigned integer maximum");
    }
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
Grid1D::Grid1D() {}
Grid1D::~Grid1D() {}
/*-------------------------------------------------*/
Grid1D::Grid1D(std::size_t n, std::size_t blockSize)
{
    input_check(n, blockSize);
    
    std::size_t numBlocks = (n + blockSize - 1) / blockSize;
    
    overflow_check(numBlocks);
    overflow_check(blockSize);
    
    m_numBlocks = dim3(numBlocks);
    m_threadsPerBlock = dim3(blockSize);
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
Grid2D::Grid2D() {}
Grid2D::~Grid2D() {}
/*-------------------------------------------------*/
Grid2D::Grid2D(std::size_t m, std::size_t n, std::size_t blockSizeX, std::size_t blockSizeY)
 {
    input_check(m, n, blockSizeX, blockSizeY);
    
    std::size_t numBlocksX = (m + blockSizeX - 1) / blockSizeX;
    std::size_t numBlocksY = (n + blockSizeY - 1) / blockSizeY;
    
    overflow_check(numBlocksX);
    overflow_check(numBlocksY);
    overflow_check(blockSizeX);
    overflow_check(blockSizeY);
    
    m_numBlocks = dim3(numBlocksX, numBlocksY);
    m_threadsPerBlock = dim3(blockSizeX, blockSizeY);
 }
/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/
