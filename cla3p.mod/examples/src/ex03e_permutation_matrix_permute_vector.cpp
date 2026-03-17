/**
 * @example ex03e_permutation_matrix_permute_vector.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>
#include <cla3p/perms.hpp>
#include <cla3p/algebra.hpp>

int main()
{
    cla3p::dns::RdVector x(3);

    for(cla3p::uint_t i = 0; i < 3; i++)
        x(i) = 3 - i;

    cla3p::prm::PiMatrix P(3);
    P(0) = 1;
    P(1) = 2;
    P(2) = 0;

    std::cout << "x:\n" << x << "\n";
    std::cout << "P:\n" << P << "\n";

    /*
     * Perform the operation (P * x).
     * using operators and the permute member function respectively.
     */
    cla3p::dns::RdVector PX1 = P * x;
    std::cout << "PX1:\n" << PX1 << "\n";

    cla3p::dns::RdVector PX2 = x.permuteLeft(P);
    std::cout << "PX2:\n" << PX2 << "\n";

    return 0;
}
