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

#ifndef CLA3P_GUARD_HPP_
#define CLA3P_GUARD_HPP_

/**
 * @file
 */

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/
/**
 * @brief Immutable object wrapper.
 * @tparam T_Object The object type to guard.
 * @details Provides a read-only view of an object, preventing modifications
 *          to its contents. Used when functions need to return objects with
 *          immutable data that must be protected from external changes.
 * @include ex04a_guard_const_vector.cpp
 */

template <typename T_Object>
class Guard {

    public:

        /**
         * @brief Default constructor.
         * @details Constructs an empty guard with no attached object.
         */
        Guard() {}

        /**
         * @brief Destructor.
         * @details Detaches the guarded object and destroys the guard.
         */
        ~Guard() { clear(); }

        /**
         * @brief Copy constructor.
         * @details Constructs a guard with a referenced copy of @p other.
         */
        Guard(const Guard<T_Object>& other) { setLocalObj(other.get()); }

        /**
         * @brief Copy assignment operator.
         * @details Replaces the guarded object with a referenced copy of @p other.
         */
        Guard<T_Object>& operator=(const Guard<T_Object>& other)
        {
            setLocalObj(other.get());
            return *this;
        }

        /**
         * @brief Move converting constructor.
         * @tparam U_Object The source object type.
         * @details Constructs a guard by converting and moving from @p other.
         */
        template <typename U_Object>
        Guard(Guard<U_Object>&& other) 
            : Guard(static_cast<const T_Object&>(other.get())) {} 

        /**
         * @brief Object constructor.
         * @details Constructs a guard with a referenced copy of @p obj.
         * @param[in] obj The object to guard.
         */
        explicit Guard(const T_Object& obj) { setLocalObj(obj); }

        /**
         * @brief Clears the guard.
         * @details Detaches the guarded object and resets the guard to an empty state.
         */
        void clear() { m_obj.clear(); }

        /**
         * @brief Retrieves the guarded object.
         * @details Returns a constant reference to the object being guarded.
         * @return A constant reference to the guarded object.
         */
        const T_Object& get() const { return m_obj; }

    private:
        T_Object m_obj;

        void setLocalObj(const T_Object& obj)
        {
            m_obj.clear();
            m_obj = const_cast<T_Object&>(obj).rcopy(); 
        }
};
/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_GUARD_HPP_
