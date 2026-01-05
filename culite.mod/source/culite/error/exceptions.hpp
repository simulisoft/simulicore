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

#ifndef CULITE_EXCEPTIONS_HPP_
#define CULITE_EXCEPTIONS_HPP_

/** 
 * @file
 */

#include <string>
#include <exception>

/*-------------------------------------------------*/
namespace culite {
namespace err {
/*-------------------------------------------------*/

/**
 * @brief Basic culite exception class.
 */
class Exception : public std::exception {

	public:
		/** 
		 * @brief The default constructor.
		 * @details Constructs an exception with the given message.
		 */
		Exception(const std::string& msg = "Untitled");

		/** 
		 * @brief Destroys the exception.
		 */
		~Exception() throw();

		/** 
		 * @brief Returns the explanatory string.
		 */
		const char *what() const throw();

	private:
		std::string m_msg;
};

/*-------------------------------------------------*/

/**
 * @brief CUDA exception class.
 */
class CudaException : public Exception {
	public:
		/** 
		 * @copydoc Exception::Exception()
		 */
		CudaException(const std::string& msg = "Untitled");

		/** 
		 * @copydoc Exception::~Exception()
		 */
		~CudaException() throw();
};

/*-------------------------------------------------*/
} // namespace err
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_EXCEPTIONS_HPP_
