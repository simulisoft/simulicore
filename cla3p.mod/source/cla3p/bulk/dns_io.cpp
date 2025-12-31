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
#include "cla3p/bulk/dns_io.hpp"

// system
#include <iomanip>
#include <algorithm>

// 3rd

// cla3p
#include "cla3p/bulk/dns.hpp"
#include "cla3p/support/utils.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace blk {
namespace dns {
/*-------------------------------------------------*/
class Printer {

	public:
		Printer(std::ostream& os, std::streamsize prec, int_t cols_per_page);
		~Printer();

		template <typename T_Scalar>
		void printToStream(uplo_t uplo, int_t m, int_t n, const T_Scalar *a, int_t lda);

	private:
		template <typename T_Scalar>
		int_t pure_element_length(uplo_t uplo, int_t m, int_t n, const T_Scalar *a, int_t lda) const;

		template <typename T_Scalar>
		void initialize(uplo_t uplo, int_t m, int_t n, const T_Scalar *a, int_t lda);

		void restore_settings();

		template <typename T_Scalar>
		void streamData(uplo_t uplo, int_t m, int_t n, const T_Scalar *a, int_t lda);

		template <typename T_Scalar>
		void streamPage(uplo_t uplo, int_t m, int_t jbgn, int_t jend, const T_Scalar *a, int_t lda);

		void streamPageHeader(int_t jbgn, int_t jend);

		template <typename T_Scalar>
		void streamIthRowOfPage(uplo_t uplo, int_t i, int_t jbgn, int_t jend, const T_Scalar *a, int_t lda);

		void streamEmpty();

		template <typename T_Scalar>
		void streamElement(const T_Scalar& a);

		std::ostream& m_os;
		std::streamsize m_prec;
		const int_t m_cols_per_page;
		const int_t m_separator_length;
		OstreamSettings m_settings;

		int_t m_row_numdigits;
		int_t m_col_numdigits;
		int_t m_record_length;
};
/*-------------------------------------------------*/
static int_t sanitized_cols_per_page(int_t cols_per_page)
{
	const int_t cols_per_page_min =  1;
	//const int_t cols_per_page_max = 50;

	cols_per_page = std::max(cols_per_page, cols_per_page_min);
	//cols_per_page = std::min(cols_per_page, cols_per_page_max);

	return cols_per_page;
}
/*-------------------------------------------------*/
Printer::Printer(std::ostream& os, std::streamsize prec, int_t cols_per_page)
	: 
		m_os(os),
		m_prec(prec > 0 ? prec : os.precision()), 
		m_cols_per_page(sanitized_cols_per_page(cols_per_page)),
		m_separator_length(1),
		m_settings(os)
{
	m_row_numdigits = 0;
	m_col_numdigits = 0;
	m_record_length = 0;
}
/*-------------------------------------------------*/
Printer::~Printer()
{
}
/*-------------------------------------------------*/
static int_t calc_max_ilen(uplo_t uplo, int_t m, int_t n, const int_t *a, int_t lda)
{
	int_t maxlen = 0;

	for(int_t j = 0; j < n; j++) {
		RowRange ir = irange(uplo, m, j);
		for(int_t i = ir.ibgn; i < ir.iend; i++) {
			maxlen = std::max(maxlen, inumlen(entry(lda,a,i,j)));
		} // i
	} // j

	return maxlen;
}
/*-------------------------------------------------*/
template <> int_t Printer::pure_element_length(uplo_t uplo, int_t m, int_t n, const int_t *a, int_t lda) const 
{ 
	return (calc_max_ilen(uplo,m,n,a,lda) + 1); 
}
/*-------------------------------------------------*/
static int_t pure_real_length(std::streamsize prec) { return (prec + 7); }
static int_t pure_complex_length(std::streamsize prec) { return (2 * pure_real_length(prec) + 3); }
/*-------------------------------------------------*/
template <> int_t Printer::pure_element_length(uplo_t, int_t, int_t, const real_t    *, int_t) const { return pure_real_length(m_os.precision()); }
template <> int_t Printer::pure_element_length(uplo_t, int_t, int_t, const real4_t   *, int_t) const { return pure_real_length(m_os.precision()); }
template <> int_t Printer::pure_element_length(uplo_t, int_t, int_t, const complex_t *, int_t) const { return pure_complex_length(m_os.precision()); }
template <> int_t Printer::pure_element_length(uplo_t, int_t, int_t, const complex8_t*, int_t) const { return pure_complex_length(m_os.precision()); }
/*-------------------------------------------------*/
template <typename T_Scalar>
void Printer::initialize(uplo_t uplo, int_t m, int_t n, const T_Scalar *a, int_t lda)
{
	//
	// Keep this first
	//
	m_settings.backup();
	m_os << std::scientific << std::setprecision(m_prec) << std::setfill(' ');

	m_row_numdigits = inumlen(m);
	m_col_numdigits = inumlen(n);

	// 
	// calculate pure element length (no spaces, with sign)
	//
	int_t element_length = pure_element_length(uplo, m, n, a, lda);
	m_record_length = std::max(element_length, m_col_numdigits) + m_separator_length;
}
/*-------------------------------------------------*/
void Printer::restore_settings()
{
	m_settings.restore();
}
/*-------------------------------------------------*/
void Printer::streamPageHeader(int_t jbgn, int_t jend)
{
	m_os << std::setw(m_row_numdigits + 2) << "";

	for(int_t j = jbgn; j < jend; j++) {
		m_os << std::setw(m_record_length) << j;
	} // j

	m_os << "\n";
}
/*-------------------------------------------------*/
template <> void Printer::streamElement(const int_t& a)
{
	m_os << std::setw(m_record_length) << a;
}
/*-------------------------------------------------*/
template <> void Printer::streamElement(const real_t& a)
{
	m_os << std::setw(m_record_length) << a;
}
/*-------------------------------------------------*/
template <> void Printer::streamElement(const real4_t& a)
{
	m_os << std::setw(m_record_length) << a;
}
/*-------------------------------------------------*/
template <typename T>
static void streamComplexElement(std::ostream& os, const T& a, int_t separator_length)
{
	int_t real_width = pure_real_length(os.precision());

	os << std::setw(separator_length) << "";
	os << "(";
	os << std::setw(real_width) << a.real();
	os << ",";
	os << std::setw(real_width) << a.imag();
	os << ")";
}
/*-------------------------------------------------*/
template <> void Printer::streamElement(const complex_t& a)
{
	streamComplexElement(m_os, a, m_separator_length);
}
/*-------------------------------------------------*/
template <> void Printer::streamElement(const complex8_t& a)
{
	streamComplexElement(m_os, a, m_separator_length);
}
/*-------------------------------------------------*/
void Printer::streamEmpty()
{
	m_os << std::setw(m_record_length) << "";
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void Printer::streamIthRowOfPage(uplo_t uplo, int_t i, int_t jbgn, int_t jend, const T_Scalar *a, int_t lda)
{
	m_os << std::setw(m_row_numdigits) << i << " |";

	for(int_t j = jbgn; j < jend; j++) {

		if(uplo == uplo_t::Upper && j < i) {
			streamEmpty();
		} else if(uplo == uplo_t::Lower && j > i) {
			break;
		} else {
			streamElement(entry(lda,a,i,j));
		}

	} // j

	m_os << "\n";
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void Printer::streamPage(uplo_t uplo, int_t m, int_t jbgn, int_t jend, const T_Scalar *a, int_t lda)
{
	streamPageHeader(jbgn, jend);
	
	for(int_t i = 0; i < m; i++) {
		streamIthRowOfPage(uplo, i, jbgn, jend, a, lda);
	} // i

	m_os << "\n";
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void Printer::streamData(uplo_t uplo, int_t m, int_t n, const T_Scalar *a, int_t lda)
{
	int_t columns_per_page = m_cols_per_page;

	int_t num_pages = n / columns_per_page;
	int_t rem_cols = n % columns_per_page;

	for(int_t ipage = 0; ipage < num_pages; ipage++) {
		int_t jbgn = ipage * columns_per_page;
		int_t jend = jbgn + columns_per_page;
		streamPage(uplo, m, jbgn, jend, a, lda);
	} // ipage

	if(rem_cols) {
		int_t jbgn = num_pages * columns_per_page;
		int_t jend = n;
		streamPage(uplo, m, jbgn, jend, a, lda);
	} // rem_cols
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void Printer::printToStream(uplo_t uplo, int_t m, int_t n, const T_Scalar *a, int_t lda)
{
	initialize(uplo, m, n, a, lda);
	streamData(uplo, m, n, a, lda);
	restore_settings();
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template <typename T_Scalar>
void print_to_stream_tmpl(std::ostream& os, uplo_t uplo, int_t m, int_t n, T_Scalar *a, int_t lda, std::streamsize prec, int_t cols_per_page)
{
	if(!(m > 0) || !(n > 0)) return;

	Printer printer(os, prec, cols_per_page);
	printer.printToStream(uplo, m, n, a, lda);
}
/*-------------------------------------------------*/
void print_to_stream(std::ostream& os, uplo_t uplo, int_t m, int_t n, const int_t *a, int_t lda, std::streamsize prec, int_t cols_per_page)
{
	return print_to_stream_tmpl(os, uplo, m, n, a, lda, prec, cols_per_page);
}
/*-------------------------------------------------*/
void print_to_stream(std::ostream& os, uplo_t uplo, int_t m, int_t n, const real_t *a, int_t lda, std::streamsize prec, int_t cols_per_page)
{
	return print_to_stream_tmpl(os, uplo, m, n, a, lda, prec, cols_per_page);
}
/*-------------------------------------------------*/
void print_to_stream(std::ostream& os, uplo_t uplo, int_t m, int_t n, const real4_t *a, int_t lda, std::streamsize prec, int_t cols_per_page)
{
	return print_to_stream_tmpl(os, uplo, m, n, a, lda, prec, cols_per_page);
}
/*-------------------------------------------------*/
void print_to_stream(std::ostream& os, uplo_t uplo, int_t m, int_t n, const complex_t *a, int_t lda, std::streamsize prec, int_t cols_per_page)
{
	return print_to_stream_tmpl(os, uplo, m, n, a, lda, prec, cols_per_page);
}
/*-------------------------------------------------*/
void print_to_stream(std::ostream& os, uplo_t uplo, int_t m, int_t n, const complex8_t *a, int_t lda, std::streamsize prec, int_t cols_per_page)
{
	return print_to_stream_tmpl(os, uplo, m, n, a, lda, prec, cols_per_page);
}
/*-------------------------------------------------*/
} // namespace dns
} // namespace blk
} // namespace cla3p
/*-------------------------------------------------*/
