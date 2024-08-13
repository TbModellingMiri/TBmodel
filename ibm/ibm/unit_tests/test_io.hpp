#include <iostream>
#include <vector>
#include <string>

#include "../inc/io.h"

bool test_io() {
    std::vector<std::vector<int>> vov {
        {1, 2, 3, 4},
        {0, 0},
        {10, 20, 100}
    };
    print_vov(vov);
    std::vector<std::string> vov_names {
        "v1234", "v00", "v1020100"
    };
    print_vov(vov, vov_names);
    
    return true;
}