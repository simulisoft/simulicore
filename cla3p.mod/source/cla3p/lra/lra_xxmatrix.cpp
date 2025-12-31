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
#include "cla3p/lra/lra_xxmatrix.hpp"

// system

// 3rd

// cla3p
#include "cla3p/support/utils.hpp"
#include "cla3p/dense.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace lra {
/*-------------------------------------------------*/
template <typename T_Matrix>
XxMatrix<T_Matrix>::XxMatrix()
{
}
/*-------------------------------------------------*/
template <typename T_Matrix>
XxMatrix<T_Matrix>::XxMatrix(int_t nr, int_t nc, int_t k)
{
	if(nr > 0 && nc > 0 && k > 0) {

		A() = T_Matrix(nr, k);
		B() = T_Matrix(nc, k);

	} else {

		// TODO: enable checks ?
		// A.checker();
		// B.checker();

	} // nr/nc/k
}
/*-------------------------------------------------*/
template <typename T_Matrix>
XxMatrix<T_Matrix>::~XxMatrix()
{
	clear();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
XxMatrix<T_Matrix>::XxMatrix(const XxMatrix<T_Matrix>& other)
	: XxMatrix(other.nrows(), other.ncols(), other.nrank())
{
	copyFromExisting(other);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
XxMatrix<T_Matrix>& XxMatrix<T_Matrix>::operator=(const XxMatrix<T_Matrix>& other)
{
  if(!(*this)) {
    *this = XxMatrix(other.nrows(), other.ncols(), other.nrank());
  }
  copyFromExisting(other);
  return *this;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
XxMatrix<T_Matrix>::XxMatrix(XxMatrix<T_Matrix>&& other)
{
	moveFrom(other);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
XxMatrix<T_Matrix>& XxMatrix<T_Matrix>::operator=(XxMatrix<T_Matrix>&& other)
{
	moveFrom(other);
	return *this;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void XxMatrix<T_Matrix>::operator=(T_Scalar val)
{
	fill(val);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void XxMatrix<T_Matrix>::clear()
{
	Meta2D::clear();

	A().clear();
	B().clear();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
int_t XxMatrix<T_Matrix>::nrows() const
{
	return A().nrows();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
int_t XxMatrix<T_Matrix>::ncols() const
{
	return B().nrows();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
int_t XxMatrix<T_Matrix>::nrank() const
{
	return A().ncols();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
const T_Matrix& XxMatrix<T_Matrix>::A() const
{
	return m_A;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
T_Matrix& XxMatrix<T_Matrix>::A()
{
	return m_A;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
const T_Matrix& XxMatrix<T_Matrix>::B() const
{
	return m_B;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
T_Matrix& XxMatrix<T_Matrix>::B()
{
	return m_B;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void XxMatrix<T_Matrix>::fill(T_Scalar val)
{
	A().fill(val);
	B().fill(val);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void XxMatrix<T_Matrix>::iscale(T_Scalar val)
{
	if(nrows() > ncols())
		B().iscale(val);
	else
		A().iscale(val);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
std::string XxMatrix<T_Matrix>::info(const std::string& header) const
{
	std::string top;
	std::string bottom;
	fill_info_margins(header, top, bottom);

	std::ostringstream ss;

	ss << top << "\n";

	ss << "  Datatype............. " << TypeTraits<typename T_Matrix::value_type>::type_name() << "\n";
	ss << "  Precision............ " << TypeTraits<typename T_Matrix::value_type>::prec_name() << "\n";
	ss << "  Number of rows....... " << nrows() << "\n";
	ss << "  Number of columns.... " << ncols() << "\n";
	ss << "  Rank................. " << nrank() << "\n";
	ss << "  A values / owner..... " << A().values() << " / " << boolToYesNo(A().owner()) << "\n";
	ss << "  B values / owner..... " << B().values() << " / " << boolToYesNo(B().owner()) << "\n";

	ss << bottom << "\n";

	return ss.str();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void XxMatrix<T_Matrix>::moveFrom(XxMatrix<T_Matrix>& other)
{
	if(this != &other) {

		if(*this) {
			*this = other;
		} else {
			Meta2D<int_t>::operator=(std::move(other));
			A() = std::move(other.A());
			B() = std::move(other.B());
		} // similar

		other.clear();

	} // do not apply on self
}
/*-------------------------------------------------*/
template <typename T_Matrix>
void XxMatrix<T_Matrix>::copyFromExisting(const XxMatrix<T_Matrix>& other)
{
	if(this != &other) {
		A() = other.A();
		B() = other.B();
	} // do not apply on self
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template class XxMatrix<dns::RdMatrix>;
template class XxMatrix<dns::RfMatrix>;
template class XxMatrix<dns::CdMatrix>;
template class XxMatrix<dns::CfMatrix>;
/*-------------------------------------------------*/
} // namespace lra
} // namespace cla3p
/*-------------------------------------------------*/

