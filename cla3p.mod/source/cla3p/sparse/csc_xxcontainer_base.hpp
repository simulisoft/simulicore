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

#ifndef CLA3P_CSC_XXCONTAINER_BASE_HPP_
#define CLA3P_CSC_XXCONTAINER_BASE_HPP_

/**
 * @file
 */

/*-------------------------------------------------*/
namespace cla3p { 
namespace csc {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The sparse container base class (compressed sparse column).
 */
template <typename T_Int, typename T_Scalar>
class XxContainerBase {

	public:
		using index_type = T_Int;
		using value_type = T_Scalar;

	public:
		XxContainerBase()
		{
			defaults();
		}

		explicit XxContainerBase(T_Int *cptr, T_Int *ridx, T_Scalar *vals)
		{
			setColptr(cptr);
			setRowidx(ridx);
			setValues(vals);
		}

		XxContainerBase(XxContainerBase<T_Int,T_Scalar>&) = delete;
		XxContainerBase<T_Int,T_Scalar>& operator=(XxContainerBase<T_Int,T_Scalar>&) = delete;

		XxContainerBase(XxContainerBase<T_Int,T_Scalar>&& other)
		{
			moveFrom(other);
		}

		XxContainerBase<T_Int,T_Scalar>& operator=(XxContainerBase<T_Int,T_Scalar>&& other)
		{
			moveFrom(other);
			return *this;
		}

		~XxContainerBase()
		{
			clear();
		}

		/**
		 * @copydoc standard_docs::colptr()
		 */
		T_Int* colptr() { return m_colptr; }

		/**
		 * @copydoc standard_docs::colptr()
		 */
		const T_Int* colptr() const { return m_colptr; }

		/**
		 * @copydoc standard_docs::rowidx()
		 */
		T_Int* rowidx() { return m_rowidx; }

		/**
		 * @copydoc standard_docs::rowidx()
		 */
		const T_Int* rowidx() const { return m_rowidx; }

		/**
		 * @copydoc standard_docs::values()
		 */
		T_Scalar* values() {  return m_values; }

		/**
		 * @copydoc standard_docs::values()
		 */
		const T_Scalar* values() const { return m_values; }

	protected:
		void clear()
		{
			defaults();
		}

	private:
		T_Int*    m_colptr;
		T_Int*    m_rowidx;
		T_Scalar* m_values;

		void setColptr(T_Int* cptr) { m_colptr = cptr; }
		void setRowidx(T_Int* ridx) { m_rowidx = ridx; }
		void setValues(T_Scalar* vals) { m_values = vals; }

		void defaults()
		{
			setColptr(nullptr);
			setRowidx(nullptr);
			setValues(nullptr);
		}

		void moveFrom(XxContainerBase<T_Int,T_Scalar>& other)
		{
			if(this != &other) {
				clear();
				setColptr(other.colptr());
				setRowidx(other.rowidx());
				setValues(other.values());
				other.clear();
			} // do not apply on self
		}
};

/*-------------------------------------------------*/
} // namespace csc
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_CSC_XXCONTAINER_BASE_HPP_
