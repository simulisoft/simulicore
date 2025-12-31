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
#include "cla3p/support/utils.hpp"

// system
#include <iterator>

// 3rd

// cla3p
#include "cla3p/error/exceptions.hpp"
#include "cla3p/support/rand.hpp"

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/
int_t inumlen(int_t n)
{
	int_t len = 1;

	while(n /= 10)
		len++;

	return len;
}
/*-------------------------------------------------*/
static constexpr std::size_t SIZEKB = 1024LLU;
static constexpr std::size_t SIZEMB = 1048576LLU;
static constexpr std::size_t SIZEGB = 1073741824LLU;
static constexpr std::size_t SIZETB = 1099511627776LLU;
/*-------------------------------------------------*/
static std::string bytesToString(std::size_t nbytes, std::size_t base, const char *suffix)
{
	real_t scaledBytes = static_cast<real_t>(nbytes) / static_cast<real_t>(base);
	std::stringstream ss;
	ss << std::fixed << std::setprecision(3) << scaledBytes << ' ' << suffix;
	return ss.str();
}
/*-------------------------------------------------*/
std::string bytesToString(std::size_t nbytes)
{
	if(nbytes >= SIZETB)
		return bytesToString(nbytes, SIZETB, "Tb");

	else if(nbytes >= SIZEGB)
		return bytesToString(nbytes, SIZEGB, "Gb");

	else if(nbytes >= SIZEMB)
		return bytesToString(nbytes, SIZEMB, "Mb");

	else if(nbytes >= SIZEKB)
		return bytesToString(nbytes, SIZEKB, "Kb");

	else
		return bytesToString(nbytes, 1, "b");
}
/*-------------------------------------------------*/
std::string kbytesToString(std::size_t nkbytes)
{
	return bytesToString(nkbytes * SIZEKB);
}
/*-------------------------------------------------*/
std::string mbytesToString(std::size_t nmbytes)
{
	return bytesToString(nmbytes * SIZEMB);
}
/*-------------------------------------------------*/
void fill_info_margins(const std::string& header, std::string& top, std::string& bottom)
{
	const std::size_t deflen = 44;
	const std::size_t msglen = header.size();

	top.clear();
	bottom.clear();

	if(!msglen) {
		top.resize(deflen, '=');
	} else if(msglen > deflen - 2) {
		top = header;
	} else {
		std::size_t lenmargin = (deflen - msglen - 2) / 2;
		top.append(lenmargin, '=');
		top.append(" ");
		top.append(header);
		top.append(" ");
		top.append(lenmargin, '=');
	} // msglen

	bottom.resize(top.size(), '=');
}
/*-------------------------------------------------*/
std::string boolToYesNo(bool flg)
{
	return (flg ? "Yes" : "No");
}
/*-------------------------------------------------*/
std::string boolToOnOff(bool flg)
{
	return (flg ? "On" : "Off");
}
/*-------------------------------------------------*/
template <typename T_Scalar>
T_Scalar opposite_element(const T_Scalar& x, const prop_t& ptype)
{
  if(ptype == prop_t::Symmetric)
		return x;

  if(ptype == prop_t::Hermitian)
		return arith::conj(x);

  if(ptype == prop_t::Skew)
		return (-x);

  throw err::Exception("Invalid Property");
  return x;
}
/*-------------------------------------------------*/
template int_t      opposite_element(const int_t     &, const prop_t&);
template real_t     opposite_element(const real_t    &, const prop_t&);
template real4_t    opposite_element(const real4_t   &, const prop_t&);
template complex_t  opposite_element(const complex_t &, const prop_t&);
template complex8_t opposite_element(const complex8_t&, const prop_t&);
/*-------------------------------------------------*/
std::string coordToString(int_t i, int_t j)
{
	std::stringstream ss;
	ss <<  "(" << i << ", " << j << ")";
  return ss.str();
}
/*-------------------------------------------------*/
static void valToStringInt(std::ostream& os, const int_t& v)
{
	os << v;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
static void valToStringReal(std::ostream& os, const T_Scalar& v)
{
	os << std::setw(os.precision() + 7) << v;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
static void valToStringComplex(std::ostream& os, const T_Scalar& v)
{
	os << "(";
	valToStringReal(os, v.real());
	os << ",";
	valToStringReal(os, v.imag());
	os << ")";
}
/*-------------------------------------------------*/
template <> void streamValue(std::ostream& os, const int_t&      v) { valToStringInt    (os, v); }
template <> void streamValue(std::ostream& os, const real_t&     v) { valToStringReal   (os, v); }
template <> void streamValue(std::ostream& os, const real4_t&    v) { valToStringReal   (os, v); }
template <> void streamValue(std::ostream& os, const complex_t&  v) { valToStringComplex(os, v); }
template <> void streamValue(std::ostream& os, const complex8_t& v) { valToStringComplex(os, v); }
/*-------------------------------------------------*/
ListPrinter::ListPrinter(std::ostream& os, int_t maxRows, int_t maxCols, int_t maxNnz, std::streamsize prec)
	: 
		m_os(os),
		m_settings(os),
		m_ndCount(std::max(int_t(5), inumlen(maxNnz))),
		m_ndRows(std::max(int_t(5), inumlen(maxRows))),
		m_ndCols(std::max(int_t(8), inumlen(maxCols)))
{
	m_settings.backup();
	m_os << std::scientific;
	m_os << std::setfill(' ');
	if(prec > 0) {
		m_os << std::setprecision(prec);
	}
}
/*-------------------------------------------------*/
ListPrinter::~ListPrinter()
{
	m_settings.restore();
}
/*-------------------------------------------------*/
void ListPrinter::streamHeader() const
{
	m_os << std::setw(m_ndCount) << "#nz";
	m_os << " | ";
	m_os << std::setw(m_ndRows) << "row";
	m_os << ' ';
	m_os << std::setw(m_ndCols) << "column";
	m_os << ' ';
	m_os << "  value";
	m_os << '\n';

	std::size_t headerLength = m_ndCount + m_ndRows + m_ndCols + 12;

	std::fill_n(std::ostream_iterator<char>(m_os), headerLength, '-');

	m_os << '\n';
}
/*-------------------------------------------------*/
decomp_t determineDecompType(decomp_t dtype, const Property& pr)
{
	decomp_t ret = dtype;

	if(ret == decomp_t::Auto) {

		if(pr.isSymmetric() || pr.isHermitian())
			ret = decomp_t::LDLT;

		else if(pr.isGeneral())
			ret = decomp_t::LU;

	} // auto

	return ret;
}
/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/
