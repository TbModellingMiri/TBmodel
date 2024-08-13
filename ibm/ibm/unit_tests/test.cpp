#include <iostream>
#include <cassert>

#include "test_io.hpp"
#include "test_common.hpp"

int main() {
    std::cout << "== io : \n";
    assert(test_io());
    std::cout << "-- done\n\n";

    std::cout << "== common : \n";
    test_rand();
    std::cout << "-- done\n\n";
}