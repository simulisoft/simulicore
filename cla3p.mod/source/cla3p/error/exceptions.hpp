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

#ifndef CLA3P_EXCEPTIONS_HPP_
#define CLA3P_EXCEPTIONS_HPP_

/** 
 * @file
 */

#include <string>
#include <exception>

#include <cla3p/types.hpp>

/*-------------------------------------------------*/
namespace cla3p {
namespace err {
/*-------------------------------------------------*/

/**
 * @brief Base exception class for cla3p library.
 * @details Extends @c std::exception to provide library-specific error reporting.
 */
class Exception : public std::exception {

    public:
        /**
         * @brief Message constructor.
         * @details Constructs an exception with the specified message.
         * @param[in] msg The error message.
         */
        Exception(const std::string& msg = "Untitled");

        /**
         * @brief Destructor.
         * @details Destroys the exception object.
         */
        ~Exception() throw();

        /**
         * @brief Returns the exception message.
         * @details Provides the explanatory error string.
         * @return Pointer to the null-terminated error message.
         */
        const char *what() const throw();

    private:
        std::string m_msg;
};

/*-------------------------------------------------*/

/**
 * @brief Consistency violation exception.
 * @details Thrown when an operation encounters inconsistent or incompatible data.
 */
class NoConsistency : public Exception {
    public:
        /** 
         * @copydoc Exception::Exception()
         */
        NoConsistency(const std::string& msg = "Untitled");

        /** 
         * @copydoc Exception::~Exception()
         */
        ~NoConsistency() throw();
};

/*-------------------------------------------------*/

/**
 * @brief Invalid operation exception.
 * @details Thrown when an operation is not supported or not applicable in the
 *          current context.
 */
class InvalidOp : public Exception {
    public:
        /** 
         * @copydoc Exception::Exception()
         */
        InvalidOp(const std::string& msg = "Untitled");

        /** 
         * @copydoc Exception::~Exception()
         */
        ~InvalidOp() throw();
};

/*-------------------------------------------------*/

/**
 * @brief Memory allocation failure exception.
 * @details Thrown when a memory allocation request cannot be satisfied.
 */
class OutOfMemory : public Exception {
    public:
        /** 
         * @copydoc Exception::Exception()
         */
        OutOfMemory(const std::string& msg = "Untitled");

        /** 
         * @copydoc Exception::~Exception()
         */
        ~OutOfMemory() throw();
};

/*-------------------------------------------------*/

/**
 * @brief Index out-of-bounds exception.
 * @details Thrown when an array or matrix access uses an invalid index.
 */
class OutOfBounds : public Exception {
    public:
        /** 
         * @copydoc Exception::Exception()
         */
        OutOfBounds(const std::string& msg = "Untitled");

        /** 
         * @copydoc Exception::~Exception()
         */
        ~OutOfBounds() throw();
};

/*-------------------------------------------------*/
} // namespace err
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_EXCEPTIONS_HPP_
