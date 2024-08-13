#include <iostream>
#include <string>
#include <vector>
#include <algorithm> //std::fill, std::shuffle
#include <random>
#include <cassert>
#include <chrono>

#include "common.h"
#include "io.h"
#include "RegionManager.h"
#include "MobilityManager.h"


void RegionManager::log_human_movement(int hid, int from_rid, int to_rid) {
    logbook_human_movement.push_back(
        "h" + std::to_string(hid)
        + kCsv_sep_str + "r" + std::to_string(human_home_region[hid])
        + kCsv_sep_str + "r" + std::to_string(from_rid)
        + kCsv_sep_str + "r" + std::to_string(to_rid)
    );
}

RegionManager::RegionManager (
        int no_of_regions,
        const std::vector<int>& population_col,
        int seed_region,
        int init_e
    ) : kNo_of_regions(no_of_regions),
        region_populations_home(population_col.begin(), population_col.begin() + kNo_of_regions),
        region_populations_current(region_populations_home.begin(), region_populations_home.end()),
        region_s_t(kNo_of_regions, -1),
        region_e_t(kNo_of_regions, -1),
        region_i_t(kNo_of_regions, -1),
        region_l_t(kNo_of_regions, -1),
        region_a_t(kNo_of_regions, -1),
        region_r_t(kNo_of_regions, -1),
        region_prev_t(kNo_of_regions, -1),
        kSeed_region(seed_region)
    {
    
    assert(kNo_of_regions <= static_cast<int>(population_col.size()));

    assert(kNo_of_regions == static_cast<int>(region_populations_home.size()));
    assert(kNo_of_regions == static_cast<int>(region_populations_current.size()));

    int pop_sum = 0;
    for (int pop : region_populations_home) {
        pop_sum += pop;
    }
    human_home_region.resize(pop_sum, -1);
    human_current_region.resize(pop_sum, -1);
    int hh = 0;
    int rr = 0;
    for (int pop : region_populations_home) {
        for (int ii = 0; ii < pop; ii++) {
            human_home_region[hh++] = rr;
        }
        rr++;
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    
    if(kSeed_region == -1) {
        // if setting is -1, do randomise
        std::shuffle(human_home_region.begin(), human_home_region.end(), gen);
    
    } else {
    assert(seed_region >= 0); // cannot seed in bad region index
    assert(seed_region < kNo_of_regions);
    
    // Create vector of exposed in seeding region
    std::vector<int> seeding(init_e, seed_region); 

    // find start of seeding region
    int rotate_start = 0;
    if(seed_region != 0) {
        // if seed_region == 0, rotate_start = 0.
        for(int i = 0 ; i < seed_region ; i++) {
            rotate_start += region_populations_home[i];
            }
        }
    int rotate_stop = rotate_start + init_e;
    
    // remove init_e number of people from seeding region to avoid double population
    human_home_region.erase(human_home_region.begin()+rotate_start, human_home_region.begin()+rotate_stop);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    // do randomise of remaining people
    std::shuffle(human_home_region.begin(), human_home_region.end(), gen);

    // concatenate exposed in seeding region and randomised remaining individuals 
    seeding.insert(seeding.end(), human_home_region.begin(),human_home_region.end());
    human_home_region = seeding;
    }

    hh = 0;
    for (int rr : human_home_region) {
        assert(rr >= 0);
        assert(rr < kNo_of_regions);
        human_current_region[hh++] = rr;
    }


    for (int rr = 0; rr < kNo_of_regions; rr++) {
        int pop_home = std::count(human_home_region.begin(), human_home_region.end(), rr);
       assert(pop_home == region_populations_home[rr]);
        int pop_current = std::count(human_current_region.begin(), human_current_region.end(), rr);
       assert(pop_current == region_populations_current[rr]);
    }

    logbook_human_movement.clear();
    logbook_human_movement.push_back(
        "hid" + kCsv_sep_str + "home" + kCsv_sep_str + "from" + kCsv_sep_str + "to"
    );
}

RegionManager::~RegionManager() {
    std::cout << "~RegionManager() called\n";
}

void RegionManager::update_region_seir_t(const std::vector<HumanSeirStatus>& human_seir_status) {
    assert(human_seir_status.size() == human_current_region.size());
    std::fill(region_s_t.begin(), region_s_t.end(), 0);
    std::fill(region_e_t.begin(), region_e_t.end(), 0);
    std::fill(region_i_t.begin(), region_i_t.end(), 0);
    std::fill(region_l_t.begin(), region_l_t.end(), 0);
    std::fill(region_a_t.begin(), region_a_t.end(), 0);
    std::fill(region_r_t.begin(), region_r_t.end(), 0);
    int hh = 0;
    for ( const HumanSeirStatus& ss: human_seir_status) {
        switch(ss) {
            case HumanSeirStatus::SUCP : region_s_t[human_current_region[hh++]]++; break;
            case HumanSeirStatus::EXPD : region_e_t[human_current_region[hh++]]++; break;
            case HumanSeirStatus::INFT : region_i_t[human_current_region[hh++]]++; break;
            case HumanSeirStatus::LATN : region_l_t[human_current_region[hh++]]++; break;
            case HumanSeirStatus::TRTD : region_a_t[human_current_region[hh++]]++; break;
            case HumanSeirStatus::RCVD : region_r_t[human_current_region[hh++]]++; break;
        }
    }
    for (int rr = 0; rr < kNo_of_regions; rr++) {
        region_prev_t[rr] = static_cast<double>(region_i_t[rr])
                        / static_cast<double>(region_populations_current[rr]);
    }

}

