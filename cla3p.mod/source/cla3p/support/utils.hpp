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

#ifndef CLA3P_UTILS_HPP_
#define CLA3P_UTILS_HPP_

/** 
 * @file
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "cla3p/types.hpp"
#include "cla3p/error.hpp"

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/

/*
 * Returns the number of digits of n
 * Sign not included
 */
int_t inumlen(int_t n);

/*
 * Bytes to human readable string
 * KBytes to human readable string
 * MBytes to human readable string
 */
std::string bytesToString(std::size_t nbytes);
std::string kbytesToString(std::size_t nkbytes);
std::string mbytesToString(std::size_t nmbytes);

/*
 * Fills object info top & bottom margin strings
 */
void fill_info_margins(const std::string& header, std::string& top, std::string& bottom);

/*
 * Boolean to yes/no
 * Boolean to on/off
 */
std::string boolToYesNo(bool flg);
std::string boolToOnOff(bool flg);

/*
 * Given Aij returns Aji depending on property
 */
template <typename T_Scalar>
T_Scalar opposite_element(const T_Scalar&, const prop_t&);

/*
 * Coordinate to string (i,j)
 */
std::string coordToString(int_t i, int_t j);

template <typename T_Scalar>
void streamValue(std::ostream& os, const T_Scalar& v);

class OstreamSettings {

	public:
		OstreamSettings(std::ostream& os) : m_os(os) {}
		~OstreamSettings() {}

		void backup() 
		{
			m_flags = m_os.flags();
			m_precision = m_os.precision();
			m_fill = m_os.fill();
		}

		void restore() const {
			m_os.flags(m_flags);
			m_os.precision(m_precision);
			m_os.fill(m_fill);
		}

	private:
		std::ostream& m_os;
		std::ios_base::fmtflags m_flags;
		std::streamsize m_precision;
		char m_fill;
};

/*
 * Unified code for Matrix element print in list form
 */
class ListPrinter {
	public:
		ListPrinter(std::ostream& os, int_t maxRows, int_t maxCols, int_t maxNnz, std::streamsize prec = 0);
		~ListPrinter();

		void streamHeader() const;

		template <typename T_Scalar>
		void streamTuple(int_t k, int_t i, int_t j, const T_Scalar& v) const
		{
			m_os << std::setw(m_ndCount) << k;
			m_os << " | ";
			m_os << std::setw(m_ndRows) << i;
			m_os << ' ';
			m_os << std::setw(m_ndCols) << j;
			m_os << ' ';
			m_os << ' ';
			streamValue(m_os, v);
			m_os << '\n';
		}

	private:
		std::ostream& m_os;
		OstreamSettings m_settings;
		const int_t m_ndCount;
		const int_t m_ndRows;
		const int_t m_ndCols;
};

/*
 * Checks for auto & returns appropriate type
 */
decomp_t determineDecompType(decomp_t, const Property&);

/*-------------------------------------------------*/
typedef struct RowRange {
	int_t ibgn;
	int_t iend;
	int_t ilen;
} RowRange;
/*-------------------------------------------------*/
inline RowRange irange(uplo_t uplo, int_t m, int_t j)
{
	RowRange ret = {0, m, m};
	//ret.ibgn = 0;
	//ret.iend = m;
	//ret.ilen = m;

	if(!m) return ret;

	if(uplo == uplo_t::Upper) {
		ret.ibgn = 0;
		ret.iend = std::min(j+1,m);
	} else if(uplo == uplo_t::Lower) {
		ret.ibgn = std::min(j,m);
		ret.iend = m;
	} // uplo

	ret.ilen = ret.iend - ret.ibgn;

	return ret;
}
/*-------------------------------------------------*/
inline RowRange irange_strict(uplo_t uplo, int_t m, int_t j)
{
	RowRange ret = irange(uplo, m, j);

	if(!m) return ret;

	if(uplo == uplo_t::Upper) {
		ret.iend--;
	} else if(uplo == uplo_t::Lower) {
		ret.ibgn++;
	} // uplo

	ret.ilen = ret.iend - ret.ibgn;

	return ret;
}
/*-------------------------------------------------*/
inline RowRange irange_complement(uplo_t uplo, int_t m, int_t j)
{
	RowRange ret = {0, 0, 0};
	//ret.ibgn = 0;
	//ret.iend = 0;
	//ret.ilen = 0;

	if(!m) return ret;

	if(uplo == uplo_t::Upper) {
		ret.ibgn = std::min(j+1,m);
		ret.iend = m;
	} else if(uplo == uplo_t::Lower) {
		ret.ibgn = 0;
		ret.iend = std::min(j,m);
	} // uplo

	ret.ilen = ret.iend - ret.ibgn;

	return ret;
}
/*-------------------------------------------------*/
inline bool coord_in_range(uplo_t uplo, int_t i, int_t j)
{
	/**/ if(uplo == uplo_t::Upper && i > j) return false;
	else if(uplo == uplo_t::Lower && i < j) return false;

	return true;
}
/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_UTILS_HPP_
