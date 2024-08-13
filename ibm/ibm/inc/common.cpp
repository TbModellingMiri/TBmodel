#include <random>
#include <iostream>
#include <algorithm>
#include <cassert>

#include <cstring> // std::memcpy

#include <boost/random/random_device.hpp> //
#include <boost/random/mersenne_twister.hpp> // boost::mt19937
#include <boost/random/uniform_real_distribution.hpp> // boost::random::uniform_real_distribution
#include <boost/random/variate_generator.hpp>  // boost::variate_generator

#include "common.h"

bool rand_uniform_buffer = true;

std::ostream& operator<< (std::ostream& os, HumanSeirStatus seir_status) {
    switch(seir_status) {
        case HumanSeirStatus::SUCP : os << "sucp"; break;
        case HumanSeirStatus::EXPD : os << "expd"; break;
        case HumanSeirStatus::INFT : os << "inft"; break;
        case HumanSeirStatus::LATN : os << "latn"; break;
        case HumanSeirStatus::TRTD : os << "trtd"; break;
        case HumanSeirStatus::RCVD : os << "rcvd"; break;
        case HumanSeirStatus::DISD : os << "disd"; break;
        default : os << "Unknown HumanSeirStatus";
        // default    : os.setstate(std::ios_base::failbit);
    }
    return os;
}

void set_rand_uniform_buffer(bool use_buffer) {
    rand_uniform_buffer = use_buffer;
}

double get_rand_uniform() {
    // // const int kSeed = 12;
    // // std::mt19937 gen(kSeed);
    // static std::mt19937 gen(std::random_device{}());
    // static std::uniform_real_distribution<> uni_dis(0.0, 1.0);
    double rnd_number = 0.0;
    if (rand_uniform_buffer) {
        set_rand_uniform_buffered(&rnd_number, 1);
    } else {
        // return uni_dis(gen);
        set_rand_uniform_switch(&rnd_number, 1);
    }
    return rnd_number;
}

void set_rand_uniform_buffered(double* result_array, int array_size) {
    static int head(-1);
    static double random_number_buffer[kRandomness_buffer_size];

    if (head == -1) { // first time this function is called
        set_rand_uniform_switch(random_number_buffer, kRandomness_buffer_size);
        head = 0;
    }
    assert(array_size < kRandomness_buffer_size);
    assert(head <= kRandomness_buffer_size);

    if (array_size > (kRandomness_buffer_size - head)) { // once the buffer has been consumed
        set_rand_uniform_switch(random_number_buffer, kRandomness_buffer_size);
        head = 0;
    }

    std::memcpy(result_array, random_number_buffer + head, array_size * sizeof(double));
    head = head + array_size;

}

void set_rand_uniform_switch(double* result_array, int array_size) {
    // You may add a switch statement here, but we hard code the choice for the boost library here for simplicity
    set_rand_uniform_boost(result_array, array_size); // boost library
    // set_rand_uniform_std(result_array, array_size); // std library
}

void set_rand_uniform_boost(double* result_array, int array_size) {
    static boost::mt19937 b_gen;

    static bool seeded(false);
    if (!seeded) {
        b_gen.seed(time(NULL));
        seeded = true;
    }
    
    static boost::random::uniform_real_distribution<float> b_udist(0,1);
    static boost::variate_generator<boost::mt19937&, boost::random::uniform_real_distribution<float>> b_uniform_rvg(b_gen,b_udist);  
    
    float value = 0.0;
    for (int ii= 0; ii < array_size; ii++) {
        value = b_uniform_rvg();
        result_array[ii] = (value >= 1.0 ? 1.0-kEPS : value);
    }
}

void set_rand_uniform_std(double* result_array, int array_size) {
    // const int kSeed = 12;
    // std::mt19937 gen(kSeed);
    static std::mt19937 gen(std::random_device{}());
    static std::uniform_real_distribution<> uni_dis(0.0, 1.0);

    float value = 0.0;

    for( int ii= 0; ii < array_size; ii++ ){
        value = uni_dis(gen);
        result_array[ii] = (value >= 1.0 ? 1.0-kEPS : value);        
    }
}

int sample_integer(int N) {
    static std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<> dist(0, N);
    return dist(gen);
}

void do_shuffle(std::vector<HumanSeirStatus>& human_seir_status){
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(human_seir_status.begin(), human_seir_status.end(), g);
}

double get_rand_gamma(double scale, double rate) {
    static std::mt19937 gen(std::random_device{}());
    std::gamma_distribution<> gma_dis(scale, rate);
    return gma_dis(gen);
};

int get_rand_discrete(std::vector<double> weights) {
    static std::mt19937 gen(std::random_device{}());
    std::discrete_distribution<> dsc_dis(weights.begin(), weights.end());
    return dsc_dis(gen);
}

double get_rand_normal(double mean, double sd) {
    static std::mt19937 gen(std::random_device{}());
    std::normal_distribution<double> nrm_dis(mean, sd);
    return nrm_dis(gen);
}

double get_rand_lnorm(double log_mean, double log_sd) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::lognormal_distribution<> log_nrm(log_mean, log_sd);
    return log_nrm(gen);
}

std::vector<int> get_hh_sizes(int m, double log_mean, double log_sd) {
    std::vector<int> household_sizes(m);

    for (int i = 0; i < m; ++i) {
        household_sizes[i] = std::round(get_rand_lnorm(log_mean, log_sd));
    }
    
    std::replace(household_sizes.begin(), household_sizes.end(), 0, 1);

    return household_sizes;
}

// Function to sample household indices with given weights
std::vector<int> sample_households(int n, int m, const std::vector<int>& household_sizes) {
    std::vector<int> household(n);
    std::vector<double> weights(m);

    double sum = std::accumulate(household_sizes.begin(), household_sizes.end(), 0.0);
    for (int i = 0; i < m; ++i) {
        weights[i] = static_cast<double>(household_sizes[i]) / sum;
    }

    // Sample households
    for (int i = 0; i < n; ++i) {
        household[i] = get_rand_discrete(weights);
    }

    return household;
}

/*
std::vector<std::vector<double>> get_euclidean_distance(std::vector<double> x, std::vector<double> y)
{
    
    std::vector<std::vector<double>> res(x.size(), std::vector<double>(x.size(), 0));

    for (size_t i = 0; i < res.size(); i++)
    {
        for (size_t j = 0; j < res[i].size(); j++)
        {
            res[i][j] = std::sqrt(std::pow(x[i] - x[j], 2) + std::pow(y[i] - y[j], 2));
        }
    }
    return res;
}
*/

std::vector<int> get_household_members(const std::vector<int>& household, int householdNum) {
    std::vector<int> indices;
    for (int i = 0; i < household.size(); ++i) {
        if (household[i] == householdNum) {
            indices.push_back(i);
        }
    }
    return indices;
}

std::vector<std::vector<double>> get_gravity_distance(
    std::vector<int> pop,
    std::vector<double> x,
    std::vector<double> y,
    int regions)
{
    
    std::vector<std::vector<double>> res(regions, std::vector<double>(regions, 0));

    for (size_t i = 0; i < regions; i++)
    {
        for (size_t j = 0; j < regions; j++)
        {
            if(i==j){
                res[i][j] == 0.0;
            } else{
                res[i][j] = (pop[i]*pop[j])/sqrt(pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2));
            }
        }
    }
    return res;
}
std::time_t get_unix_from_date(int year, int month, int day) {
    std::tm timeinfo = {};
    timeinfo.tm_year = year - 1900;  // years since 1900
    timeinfo.tm_mon = month - 1;     // January = 0 correction
    timeinfo.tm_mday = day;

    std::time_t unix_time = std::mktime(&timeinfo);

    return unix_time;
}
/*
std::string get_date_from_unix(std::time_t unix_time) {
    std::tm* timeinfo = std::localtime(&unix_time);

    int year = timeinfo->tm_year + 1900; 
    int month = timeinfo->tm_mon + 1;    
    int day = timeinfo->tm_mday;

    std::string date = std::to_string(year) + "-" + std::to_string(month) + "-" + std::to_string(day);

    return date;
}
*/