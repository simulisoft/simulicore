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
#include "cla3p/sparse/csc_xxcontainer.hpp"

// system

// 3rd

// cla3p
#include "cla3p/support/imalloc.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace csc {
/*-------------------------------------------------*/
#define XxContainerTmpl XxContainer<T_Int,T_Scalar>
#define XxContainerTlst template <typename T_Int, typename T_Scalar>
/*-------------------------------------------------*/
XxContainerTlst
XxContainerTmpl::XxContainer()
{
	defaults();
}
/*-------------------------------------------------*/
XxContainerTlst
XxContainerTmpl::XxContainer(std::size_t nc, std::size_t nz)
{
	if(nc) {
		T_Int    *cptr = i_malloc<T_Int>(nc+1);
		T_Int    *ridx = i_malloc<T_Int>(nz);
		T_Scalar *vals = i_malloc<T_Scalar>(nz);

		Ownership::operator=(Ownership(true));

		cptr[nc] = nz;
		setColptr(cptr);
		setRowidx(ridx);
		setValues(vals);
	}
}
/*-------------------------------------------------*/
XxContainerTlst
XxContainerTmpl::XxContainer(T_Int *cptr, T_Int *ridx, T_Scalar *vals, bool bind)
{
	if(cptr) {
		Ownership::operator=(Ownership(bind));
		setColptr(cptr);
		setRowidx(ridx);
		setValues(vals);
	}
}
/*-------------------------------------------------*/
XxContainerTlst
XxContainerTmpl::~XxContainer()
{
  clear();
}
/*-------------------------------------------------*/
XxContainerTlst
void XxContainerTmpl::defaults()
{
	setColptr(nullptr);
	setRowidx(nullptr);
	setValues(nullptr);
}
/*-------------------------------------------------*/
XxContainerTlst
void XxContainerTmpl::setColptr(T_Int *cptr)
{
	m_colptr = cptr;
}
/*-------------------------------------------------*/
XxContainerTlst
void XxContainerTmpl::setRowidx(T_Int *ridx)
{
	m_rowidx = ridx;
}
/*-------------------------------------------------*/
XxContainerTlst
void XxContainerTmpl::setValues(T_Scalar *vals)
{
	m_values = vals;
}
/*-------------------------------------------------*/
XxContainerTlst
const T_Int* XxContainerTmpl::colptr() const
{
  return m_colptr;
}
/*-------------------------------------------------*/
XxContainerTlst
T_Int* XxContainerTmpl::colptr()
{
  return m_colptr;
}
/*-------------------------------------------------*/
XxContainerTlst
const T_Int* XxContainerTmpl::rowidx() const
{
  return m_rowidx;
}
/*-------------------------------------------------*/
XxContainerTlst
T_Int* XxContainerTmpl::rowidx()
{
  return m_rowidx;
}
/*-------------------------------------------------*/
XxContainerTlst
T_Scalar* XxContainerTmpl::values()
{
	return m_values;
}
/*-------------------------------------------------*/
XxContainerTlst
const T_Scalar* XxContainerTmpl::values() const
{
	return m_values;
}
/*-------------------------------------------------*/
XxContainerTlst
void XxContainerTmpl::clear()
{
	if(owner()) {
		i_free(colptr());
		i_free(rowidx());
		i_free(values());
	} // owner

	Ownership::clear();

	defaults();
}
/*-------------------------------------------------*/
#undef XxContainerTmpl
#undef XxContainerTlst
/*-------------------------------------------------*/
#define instantiate_xxcontainer(T_Scl) \
template class XxContainer<int_t,T_Scl>
instantiate_xxcontainer(real_t);
instantiate_xxcontainer(real4_t);
instantiate_xxcontainer(complex_t);
instantiate_xxcontainer(complex8_t);
#undef instantiate_xxcontainer
/*-------------------------------------------------*/
} // namespace csc
} // namespace cla3p
/*-------------------------------------------------*/
