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

#ifndef CLA3P_VIRTUAL_EXPRESSION_HPP_
#define CLA3P_VIRTUAL_EXPRESSION_HPP_

/**
 * @file
 */

#include <ostream>

/*-------------------------------------------------*/

namespace cla3p {
namespace dns { template <typename T_Scalar> class XxVector; }
namespace dns { template <typename T_Scalar> class XxMatrix; }
namespace csr { template <typename T_Int, typename T_Scalar> class XxMatrix; }
namespace csc { template <typename T_Int, typename T_Scalar> class XxMatrix; }
} // namespace cla3p

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief The virtual expression base class.
 */
template <typename T_Result, typename T_Virtual>
class VirtualExpression {

    private:
        using T_Scalar = typename T_Result::value_type;
        using virtual_type = T_Virtual;

    public:
        using result_type = T_Result;

    public:
        VirtualExpression() {}
        ~VirtualExpression() {}

        const T_Virtual& self() const { return static_cast<const T_Virtual&>(*this); }

        /**
         * @brief Evaluates the expression.
         * @details Evaluates the expression and stores the result in a new object.
         */
        virtual T_Result evaluate() const 
        {
            T_Result ret;
            evaluateOnNew(ret);
            return ret; 
        }

        /**
         * @brief Evaluates the expression on a return object.
         * @details Clears existing data in @p dest, reallocates and evaluates the expression.
         */
        virtual void evaluateOnNew(T_Result& dest) const = 0;

        /**
         * @brief Evaluates the expression on a return object.
         * @details Evaluates the expression on the pre-allocated & compatible @p dest.
         */
        virtual void evaluateOnExisting(T_Result& dest) const = 0;

        /**
         * @brief Adds the expression on a return object.
         * @details Adds the scaled expression result to the pre-allocated & compatible @p dest.
         */
        virtual void accumulateOnExisting(T_Result& dest, T_Scalar coeff) const = 0;
};

/*-------------------------------------------------*/

namespace alias { 

template <typename T_Scalar, typename T_Virtual>
using VirtualExpr_vec = VirtualExpression<dns::XxVector<T_Scalar>, T_Virtual>;

template <typename T_Scalar, typename T_Virtual>
using VirtualExpr_dns = VirtualExpression<dns::XxMatrix<T_Scalar>, T_Virtual>;

template <typename T_Int, typename T_Scalar, typename T_Virtual>
using VirtualExpr_csr = VirtualExpression<csr::XxMatrix<T_Int,T_Scalar>, T_Virtual>;

template <typename T_Int, typename T_Scalar, typename T_Virtual>
using VirtualExpr_csc = VirtualExpression<csc::XxMatrix<T_Int,T_Scalar>, T_Virtual>;

} // namespace alias

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_stream_operators
 * @brief Writes to os the contents of the evaluated v.
 */
template <typename T_Result, typename T_Virtual>
std::ostream& operator<<(std::ostream& os, const cla3p::VirtualExpression<T_Result, T_Virtual>& v)
{
    os << v.evaluate();
    return os;
}

/*-------------------------------------------------*/

#endif // CLA3P_VIRTUAL_EXPRESSION_HPP_
