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
#include <culite/generic/cusparse_handler.hpp>

// system

// 3rd

// culite
#include <culite/error/cuda.hpp>

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/
CuSparseSpm::CuSparseSpm(CuSparseHandler& cuSparseHandler, cusparsePointerMode_t mode)
    : m_handler(cuSparseHandler)
{
    m_oldMode = m_handler.setPointerMode(mode);
}
/*-------------------------------------------------*/
CuSparseSpm::~CuSparseSpm()
{
    m_handler.setPointerMode(m_oldMode);
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
CuSparseHandler::CuSparseHandler()
{
    defaults();
    cusparseStatus_t cusparseStatus = cusparseCreate(&m_handle);
    err::check_cusparse(cusparseStatus);
}
/*-------------------------------------------------*/
CuSparseHandler::~CuSparseHandler()
{
    clear();
    cusparseStatus_t cusparseStatus = cusparseDestroy(m_handle);
    err::check_cusparse(cusparseStatus);
}
/*-------------------------------------------------*/
void CuSparseHandler::defaults()
{
    m_workspaceInBytes = 0;
}
/*-------------------------------------------------*/
void CuSparseHandler::clear()
{
    deviceWork().clear();

    defaults();
}
/*-------------------------------------------------*/
cusparsePointerMode_t CuSparseHandler::setPointerMode(cusparsePointerMode_t mode)
{
    cusparsePointerMode_t ret = pointerMode();
    if (ret != mode) {
        cusparseStatus_t status = cusparseSetPointerMode(m_handle, mode);
        err::check_cusparse(status);
    }
    return ret;
}
/*-------------------------------------------------*/
cusparsePointerMode_t CuSparseHandler::pointerMode()
{
    cusparsePointerMode_t ret;
    cusparseStatus_t status = cusparseGetPointerMode(m_handle, &ret);
    err::check_cusparse(status);
    return ret;
}
/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/
