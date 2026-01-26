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

#ifndef CLA3P_CSX_XXCONTAINER_BASE_HPP_
#define CLA3P_CSX_XXCONTAINER_BASE_HPP_

/**
 * @file
 */

/*-------------------------------------------------*/
namespace cla3p { 
namespace csx {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The sparse container base class (compressed sparse generic).
 */
template <typename T_Int, typename T_Scalar>
class XxContainerBase {

	public:
		using index_type = T_Int;
		using value_type = T_Scalar;

	public:
		XxContainerBase() { defaults(); }

		explicit XxContainerBase(T_Int *xptr, T_Int *xidx, T_Scalar *vals)
		{
			setXxxptr(xptr);
			setXxxidx(xidx);
			setValues(vals);
		}

		XxContainerBase(XxContainerBase<T_Int,T_Scalar>&) = delete;
		XxContainerBase<T_Int,T_Scalar>& operator=(XxContainerBase<T_Int,T_Scalar>&) = delete;

		XxContainerBase(XxContainerBase<T_Int,T_Scalar>&& other) { moveFrom(other); }
		XxContainerBase<T_Int,T_Scalar>& operator=(XxContainerBase<T_Int,T_Scalar>&& other) { return moveFrom(other); }

		~XxContainerBase() { clear(); }

        /**
		 * @brief Access the data buffer.
		 * @details Returns a pointer to the underlying data array.
		 * @return Pointer to the data buffer.
		 */
        T_Scalar* values() { return m_values; }

        /**
		 * @brief Access the data buffer.
		 * @details Returns a pointer to the underlying data array.
		 * @return Pointer to the data buffer.
		 */
        const T_Scalar* values() const { return m_values; }

   	protected:
		T_Int* xxxptr() { return m_xxxptr; }
		T_Int* xxxidx() { return m_xxxidx; }

   		const T_Int* xxxptr() const { return m_xxxptr; }
		const T_Int* xxxidx() const { return m_xxxidx; }

		void clear()
		{
			defaults();
		}

	private:
		T_Int*    m_xxxptr;
		T_Int*    m_xxxidx;
		T_Scalar* m_values;

		void setXxxptr(T_Int*    xptr) { m_xxxptr = xptr; }
		void setXxxidx(T_Int*    xidx) { m_xxxidx = xidx; }
		void setValues(T_Scalar* vals) { m_values = vals; }

		void defaults()
		{
			setXxxptr(nullptr);
			setXxxidx(nullptr);
			setValues(nullptr);
		}

		XxContainerBase<T_Int,T_Scalar>& moveFrom(XxContainerBase<T_Int,T_Scalar>& other)
		{
			if(this != &other) {
				clear();
				setXxxptr(other.xxxptr());
				setXxxidx(other.xxxidx());
				setValues(other.values());
				other.clear();
			} // do not apply on self
            return *this;
		}
};

/*-------------------------------------------------*/
} // namespace csx
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_CSX_XXCONTAINER_BASE_HPP_
