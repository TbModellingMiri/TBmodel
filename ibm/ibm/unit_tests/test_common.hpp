#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>

#include "../inc/common.h"

void uniform_bin_print(int kNum_intervals, int kNum_draws) {
    std::vector<int> res(kNum_intervals, 0);
    for (int ii = 0; ii < kNum_draws; ii++) {
        res[std::floor(get_rand_uniform() * kNum_intervals)]++;
    }
    int ii = 0;
    for (int ee : res) {
        std::cout << "intv " << ii++ << ": " << ee << "\n";
    }
}

bool test_rand() {
    const int kNum_draws = kRandomness_buffer_size * 10;
    const int kNum_intervals = 10;

    std::cout << "testing get_rand_uniform():\n";

    // --

    set_rand_uniform_buffer(false);
    std::cout << "No buffer:\n";
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    begin = std::chrono::steady_clock::now();
    uniform_bin_print(kNum_intervals, kNum_draws);
    end = std::chrono::steady_clock::now();

    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[micro-s]" << std::endl;

    // --

    set_rand_uniform_buffer(true);
    std::cout << "With buffer:\n";

    begin = std::chrono::steady_clock::now();
    uniform_bin_print(kNum_intervals, kNum_draws);
    end = std::chrono::steady_clock::now();

    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[micro-s]" << std::endl;

    return true;
}