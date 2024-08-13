#ifndef COMMON_H
#define COMMON_H
#include <ostream>
#include <ctime>
// #include <chrono>

const double kNo_of_days_in_year = 365.0;

const int kRandomness_buffer_size = 10000000;
const float kEPS = 1.2e-7;





enum class HumanSeirStatus {
    SUCP,
    EXPD,
    INFT,
    LATN,
    TRTD,
    RCVD,
    DISD
};
std::ostream& operator<<(std::ostream& os, HumanSeirStatus seir_status);


void set_rand_uniform_buffer(bool use_buffer);
double get_rand_uniform();
void set_rand_uniform_buffered(double* result_array, int array_size);
void set_rand_uniform_switch(double* result_array, int array_size);
void set_rand_uniform_boost(double* result_array, int array_size);
void set_rand_uniform_std(double* result_array, int array_size);

void do_shuffle(std::vector<HumanSeirStatus>& human_seir_status);

int get_rand_discrete(std::vector<double> weights);
double get_rand_gamma(double scale, double rate);
int sample_integer(int N);
double get_rand_normal(double mean, double sd);
double get_rand_lnorm(double log_mean, double log_sd);
std::vector<int> get_hh_sizes(int m, double log_mean, double log_sd);
std::vector<int> sample_households(int n, int m, const std::vector<int>& household_sizes);

// class Rand {
//     // static std::random_device rd;
//     static std::mt19937 gen(std::random_device());

// public:
//     static double get_rand_uniform() {
//         static std::uniform_real_distribution<> uni_dis(0.0, 1.0);
//         return uni_dis(gen);
//     }
// };

/*
std::vector<std::vector<double>> get_euclidean_distance(std::vector<double> x, std::vector<double> y);
*/

std::vector<int> get_household_members(const std::vector<int>& household, int householdNum);

std::vector<std::vector<double>> get_gravity_distance(
    std::vector<int> pop,
    std::vector<double> x,
    std::vector<double> y,
    int regions);

    std::time_t get_unix_from_date(int year, int month, int day);

/*
std::string get_date_from_unix(std::time_t unix_time);
*/
#endif