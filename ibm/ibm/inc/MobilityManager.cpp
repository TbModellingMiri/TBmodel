#include <iostream>
#include <vector>
#include <cassert>
#include <random>
#include <cmath>
#include <chrono>

#include "io.h"
#include "common.h"
#include "RegionManager.h"
#include "MobilityManager.h"

int MobilityManager::pick_region_uniform(int home_region) {
    int r = uni_int_dis(gen);
    if (r >= home_region) {
        r = r + 1;
    }
    return(r);
}

int MobilityManager::pick_region_gravity(int home_region) {
    std::vector<double> weights = kGravity[home_region];
    int r = get_rand_discrete(weights);
    return(r);
}
// Mobility type 1
// short-term return-home travel

int MobilityManager::type1_gen_travel_length() {
    // return std::max(1, static_cast<int>(nrm_dis_type1(gen)));
    return std::max(1, static_cast<int>(std::round(nrm_dis_type1(gen))));
}

// Mobility type season
// same as type 1, but different values 
int MobilityManager::season_gen_travel_length() {
    // return std::max(1, static_cast<int>(nrm_dis_type1(gen)));
    return std::max(1, static_cast<int>(std::round(nrm_dis_season(gen))));
}

void MobilityManager::mobility_season(
    RegionManager& rgn_mgr
    ) {
    int hh = 0;
    for (int& dd: human_days_until_return_home) {
        assert(dd >= 0);
        if ( dd == 1) { // time to return home
            rgn_mgr.send_human_home(hh);
            dd = 0;
        } else if (dd == 0) { // time to travel?
            if (get_rand_uniform() < kProb_travel_season) {
                // how many days?
                dd = season_gen_travel_length();
                // where to?
                int to_region;
                if(kMobility_model==2) {
                    to_region = pick_region_gravity(rgn_mgr.get_human_current_region(hh));
                } else {
                    to_region = pick_region_uniform(rgn_mgr.get_human_current_region(hh));
                }
                rgn_mgr.send_human_to_region(hh, to_region);
            }
        } else { // spend time away from home
            dd--;
        }
        hh++;
    }
}

void MobilityManager::mobility_type1(
    RegionManager& rgn_mgr
    ) {
    int hh = 0;
    for (int& dd: human_days_until_return_home) {
        assert(dd >= 0);
        if ( dd == 1) { // time to return home
            rgn_mgr.send_human_home(hh);
            dd = 0;
        } else if (dd == 0) { // time to travel?
            if (get_rand_uniform() < kProb_travel) {
                // how many days?
                dd = type1_gen_travel_length();
                // where to?
                int to_region;
                if(kMobility_model==2) {
                    to_region = pick_region_gravity(rgn_mgr.get_human_current_region(hh));
                } else {
                    to_region = pick_region_uniform(rgn_mgr.get_human_current_region(hh));
                }
                rgn_mgr.send_human_to_region(hh, to_region);
            }
        } else { // spend time away from home
            dd--;
        }
        hh++;
    }
}


MobilityManager::MobilityManager(
        int no_of_humans,
        int no_of_regions,
        // type 1
        double prob_travel,
        double days_away_mean,
        double days_away_sd,

        double prob_travel_season,
        double days_away_mean_season,
        double days_away_sd_season,

        const std::vector<int>& population_col,
        const std::vector<double>& lon_col,
        const std::vector<double>& lat_col,
        const int mobility_model,
        const int year,
        const int month,
        const int day,
        const int holiday_start_y,
        const int holiday_start_m,
        const int holiday_start_d,
        const int holiday_length,
        const int recurrent

    ) : human_days_until_return_home(no_of_humans, 0),
        uni_int_dis(0, no_of_regions - 2),
        kProb_travel(prob_travel),
        kProb_travel_season(prob_travel_season),
        // kDays_away_mean(days_away_mean),
        // kDays_away_sd(days_away_sd),
        nrm_dis_type1(days_away_mean, days_away_sd),
        nrm_dis_season(days_away_mean_season, days_away_sd_season),

        kGravity(get_gravity_distance(population_col, lon_col, lat_col, no_of_regions)),
        kMobility_model(mobility_model),
        current_time(get_unix_from_date(year, month, day)),
        holiday_start(get_unix_from_date(holiday_start_y, holiday_start_m, holiday_start_d)),
        holiday_end(holiday_start+holiday_length*86400),
        kRecurrent(recurrent)
    {
}
MobilityManager::~MobilityManager() {
    std::cout << "~MobilityManager() called\n";
}

void MobilityManager::step(
    RegionManager& rgn_mgr
    ) {
    /*
    if(kRecurrent != -1 and current_time > holiday_start and current_time < holiday_end) {
        mobility_season(rgn_mgr);
    } else {
        mobility_type1(rgn_mgr);
    }

   if (current_time == holiday_end and kRecurrent == 1) {
    // if end of holiday and holiday is recurrent
    // add +1 year to holiday start and end date
    std::tm start_placeholder = *std::localtime(&holiday_start);
    start_placeholder.tm_year = start_placeholder.tm_year + 1;
    std::time_t updated_start_holiday = std::mktime(&start_placeholder);
    
    std::tm end_placeholder = *std::localtime(&holiday_end);
    end_placeholder.tm_year = end_placeholder.tm_year + 1;
    std::time_t updated_holiday_end = std::mktime(&end_placeholder);
    
    holiday_start = updated_start_holiday;
    holiday_end = updated_holiday_end;
}
    */
    // add 86400 seconds (one day) at the end of each day
    current_time += 86400;
}


