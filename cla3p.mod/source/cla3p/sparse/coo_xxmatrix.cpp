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
#include <cla3p/sparse/coo_xxmatrix.hpp>

// system
#include <algorithm>

// 3rd

// cla3p
#include <cla3p/sparse/csc_xxmatrix.hpp>

#include <cla3p/support/imalloc.hpp>
#include <cla3p/support/utils.hpp>
#include <cla3p/bulk/csc.hpp>

#include <cla3p/checks/basic_checks.hpp>
#include <cla3p/checks/coo_checks.hpp>

/*-------------------------------------------------*/
namespace cla3p {
namespace coo {
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar>::XxMatrix()
{
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar>::XxMatrix(T_Int nr, T_Int nc, const Property& pr)
    : MatrixMeta<T_Int>(nr, nc, sanitizeProperty<T_Scalar>(pr))
{
    if(nr > 0 && nc > 0) {
        checker();
    } else {
        clear();
    } // nr/nc
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
XxMatrix<T_Int,T_Scalar>::~XxMatrix()
{
    clear();
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void XxMatrix<T_Int,T_Scalar>::clear()
{
    MatrixMeta<T_Int>::clear();
    tupleVec().clear();
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
T_Int XxMatrix<T_Int,T_Scalar>::nnz() const
{
    return static_cast<T_Int>(tupleVec().size());
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
typename XxMatrix<T_Int,T_Scalar>::TupleVec& XxMatrix<T_Int,T_Scalar>::tupleVec()
{
    return m_tuples;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
const typename XxMatrix<T_Int,T_Scalar>::TupleVec& XxMatrix<T_Int,T_Scalar>::tupleVec() const
{
    return m_tuples;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void XxMatrix<T_Int,T_Scalar>::reserve(T_Int nz)
{
    tupleVec().reserve(nz);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void XxMatrix<T_Int,T_Scalar>::insert(const Tuple<T_Int,T_Scalar>& tuple)
{
    coo_check_triplet(this->nrows(), this->ncols(), this->prop(), tuple.row(), tuple.col(), tuple.val());

    tupleVec().push_back(tuple);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void XxMatrix<T_Int,T_Scalar>::insert(T_Int i, T_Int j, T_Scalar v)
{
    Tuple<T_Int,T_Scalar> tuple(i, j, v);
    insert(tuple);
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
std::string XxMatrix<T_Int,T_Scalar>::info(const std::string& header) const
{ 
    std::string top;
    std::string bottom;
    fill_info_margins(header, top, bottom);

    std::ostringstream ss;

    ss << top << "\n";

    ss << "  Datatype............. " << TypeTraits<T_Scalar>::type_name() << "\n";
    ss << "  Precision............ " << TypeTraits<T_Scalar>::prec_name() << "\n";
    ss << "  Index Precision...... " << TypeTraits<T_Int>::prec_name() << "\n";
    ss << "  Number of rows....... " << this->nrows() << "\n";
    ss << "  Number of columns.... " << this->ncols() << "\n";
    ss << "  Number of non zeros.. " << nnz() << "\n";
    ss << "  Property............. " << this->prop() << "\n";

    ss << bottom << "\n";

    return ss.str();
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void XxMatrix<T_Int,T_Scalar>::toStream(std::ostream& os, std::streamsize prec) const
{
    if(this->empty() || !nnz()) return;

    ListPrinter listPrinter(os, this->nrows(), this->ncols(), nnz(), prec);

    listPrinter.streamHeader();

    for(T_Int cnt = 0; cnt < nnz(); cnt++) {
            T_Int    i = tupleVec()[cnt].row();
            T_Int    j = tupleVec()[cnt].col();
            T_Scalar v = tupleVec()[cnt].val();
            listPrinter.streamTuple(cnt, i, j, v);
    } // cnt
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
csc::XxMatrix<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::toCsc(dup_t duplicatePolicy) const
{
    if(!this->nrows() || !this->ncols())
        return csc::XxMatrix<T_Int,T_Scalar>();

    T_Int *colptr = i_calloc_t<T_Int>(this->ncols() + 1);

    std::for_each(tupleVec().begin(), tupleVec().end(), 
            [&](const Tuple<T_Int,T_Scalar> &tuple) 
            { 
            colptr[tuple.col() + 1]++;
            });

    blk::csx::roll(this->ncols(), colptr);

    T_Int nnz = colptr[this->ncols()];

    T_Int    *rowidx = nullptr;
    T_Scalar *values = nullptr;

    if(nnz) {

        rowidx = i_malloc_t<T_Int>(nnz);
        values = i_malloc_t<T_Scalar>(nnz);

        std::for_each(tupleVec().begin(), tupleVec().end(), 
                [&](const Tuple<T_Int,T_Scalar> &tuple) 
                { 
                rowidx[colptr[tuple.col()]] = tuple.row();
                values[colptr[tuple.col()]] = tuple.val();
                colptr[tuple.col()]++;
                });

        blk::csx::unroll(this->ncols(), colptr);
        blk::csx::sort(this->ncols(), colptr, rowidx, values);
        blk::csx::remove_duplicates(this->ncols(), colptr, rowidx, values, duplicatePolicy);

        rowidx = i_realloc_t<T_Int>(rowidx, colptr[this->ncols()]);
        values = i_realloc_t<T_Scalar>(values, colptr[this->ncols()]);

    } // nnz

    csc::XxMatrix<T_Int,T_Scalar> ret(this->nrows(), this->ncols(), colptr, rowidx, values, true, this->prop());

    return ret;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
csr::XxMatrix<T_Int,T_Scalar> XxMatrix<T_Int,T_Scalar>::toCsr(dup_t duplicatePolicy) const
{
    if(!this->nrows() || !this->ncols())
        return csr::XxMatrix<T_Int,T_Scalar>();

    T_Int *rowptr = i_calloc_t<T_Int>(this->nrows() + 1);

    // TODO: group csr/csc functionality into a utility function

    std::for_each(tupleVec().begin(), tupleVec().end(), 
            [&](const Tuple<T_Int,T_Scalar> &tuple) 
            { 
            rowptr[tuple.row() + 1]++;
            });

    blk::csx::roll(this->nrows(), rowptr);

    T_Int nnz = rowptr[this->nrows()];

    T_Int    *colidx = nullptr;
    T_Scalar *values = nullptr;

    if(nnz) {

        colidx = i_malloc_t<T_Int>(nnz);
        values = i_malloc_t<T_Scalar>(nnz);

        std::for_each(tupleVec().begin(), tupleVec().end(), 
                [&](const Tuple<T_Int,T_Scalar> &tuple) 
                { 
                colidx[rowptr[tuple.row()]] = tuple.col();
                values[rowptr[tuple.row()]] = tuple.val();
                rowptr[tuple.row()]++;
                });

        blk::csx::unroll(this->nrows(), rowptr);
        blk::csx::sort(this->nrows(), rowptr, colidx, values);
        blk::csx::remove_duplicates(this->nrows(), rowptr, colidx, values, duplicatePolicy);

        colidx = i_realloc_t<T_Int>(colidx, rowptr[this->nrows()]);
        values = i_realloc_t<T_Scalar>(values, rowptr[this->nrows()]);

    } // nnz

    csr::XxMatrix<T_Int,T_Scalar> ret(this->nrows(), this->ncols(), rowptr, colidx, values, true, this->prop());
    return ret;
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void XxMatrix<T_Int,T_Scalar>::checker() const
{
    coo_consistency_check(this->prop(), this->nrows(), this->ncols());
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template class XxMatrix<int_t,real_t>;
template class XxMatrix<int_t,real4_t>;
template class XxMatrix<int_t,complex_t>;
template class XxMatrix<int_t,complex8_t>;
/*-------------------------------------------------*/
} // namespace coo
} // namespace cla3p
/*-------------------------------------------------*/
