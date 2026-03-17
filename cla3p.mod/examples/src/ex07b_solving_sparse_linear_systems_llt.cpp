/**
 * @example ex07b_solving_sparse_linear_systems_llt.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>
#include <cla3p/sparse.hpp>
#include <cla3p/linsol.hpp>
#include <cla3p/algebra.hpp>

/*-----------------------------------------------------*/
static cla3p::csr::RdMatrix DefaultSymmetricSparseMatrix()
{
    static cla3p::int_t  rowptr[] = {0, 3, 5, 7,  8, 9};
    static cla3p::int_t  colidx[] = { 0,  1,  3,  1, 4,  2, 3,  3,  4};
    static cla3p::real_t values[] = {10, -2, -4, 50, 8, 40, 2, 70, 50};

    return cla3p::csr::RdMatrix(5, 5, rowptr, colidx, values, false, cla3p::Property::SymmetricUpper());
}
/*-----------------------------------------------------*/
int main()
{
    const cla3p::csr::RdMatrix A = DefaultSymmetricSparseMatrix();
    const cla3p::dns::RdVector b = cla3p::dns::RdVector::random(5);
    const cla3p::dns::RdMatrix B = cla3p::dns::RdMatrix::random(5, 3);
    cla3p::dns::RdVector x; // x will be created in solve
    cla3p::dns::RdMatrix X(5, 3); // Preallocate space for X

    cla3p::PardisoLLt<cla3p::csr::RdMatrix> lltSolver;

    /*
     * Perform analysis & symbolic decomposition on A.
     */
    lltSolver.analysis(A);

    /*
     * Decompose A into LL' product.
     */
    lltSolver.decompose(A);

    {
        /*
         * Single column (vector) rhs.
         * Calculate x = (A^{-1} * b).
         */

        lltSolver.solve(b, x);
        std::cout << "Dense Vector rhs::Absolute Error: "
                  << (b - A * x).evaluate().normOne() << std::endl;
    }

    {
        /*
         * Multiple column (matrix) rhs.
         * Calculate X = (A^{-1} * B).
         */
        lltSolver.solve(B, X);
        std::cout << "Dense Matrix rhs::Absolute Error: "
                  << (B - A * X).evaluate().normOne() << std::endl;
    }

    return 0;
}
