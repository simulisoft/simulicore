/**
 * @example ex00_property_create.cpp
 */

#include <iostream>
#include <cla3p/types.hpp>

int main()
{
    /*
     * Instantiate various properties
     */

    cla3p::Property prGe = cla3p::Property::General(); // the default property
    cla3p::Property prSy = cla3p::Property::SymmetricLower(); // symmetric lower property
    cla3p::Property prHe = cla3p::Property(cla3p::prop_t::Hermitian, cla3p::uplo_t::Upper); // hermitian upper property

    std::cout << "Property G: " << prGe << std::endl;
    std::cout << "Property S: " << prSy << std::endl;
    std::cout << "Property H: " << prHe << std::endl;

    /*
     * Instantiate an invalid property
     */

    try {
        cla3p::Property invalidProperty(cla3p::prop_t::General, cla3p::uplo_t::Upper);
    } catch(...) {
        std::cout << "Exception thrown for invalid property" << std::endl;
    }

    return 0;
}

