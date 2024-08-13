#ifndef REGION_MAN_H
#define REGION_MAN_H

class RegionManager {
    const int kNo_of_regions;

    // size kNo_of_regions;
    std::vector<int> region_populations_home;
    std::vector<int> region_populations_current;

    std::vector<int> region_s_t;
    std::vector<int> region_e_t;
    std::vector<int> region_i_t;
    std::vector<int> region_l_t;
    std::vector<int> region_a_t;
    std::vector<int> region_r_t;

    std::vector<double> region_prev_t; //prevalence

    const int kSeed_region;


    // size kNo_of_humans
    std::vector<int> human_home_region;
    std::vector<int> human_current_region;

    // friend class MobilityManager;
    // friend void MobilityManager::step(RegionManager& rgn_mgr);
    // friend void MobilityManager::send_human_home(RegionManager& rgn_mgr);

    std::vector<std::string> logbook_human_movement;

// FUNCTIONS

    void log_human_movement(int hid, int from_rid, int to_rid);

public:
    RegionManager(
        int no_of_regions,
        const std::vector<int>& population_col,
        int seed_region,
        int init_e
    );
    ~RegionManager();
    
    void update_region_seir_t(
        const std::vector<HumanSeirStatus>& human_seir_status
    );
    void step();

    // called by MobilityManager
    inline void send_human_to_region(int hid, int rid) {
        log_human_movement(hid, human_current_region[hid], rid);
        region_populations_current[human_current_region[hid]]--;
        region_populations_current[rid]++;
        human_current_region[hid] = rid;
    }
    inline void send_human_home(int hid) {
        // region_populations_current[human_current_region[hid]]--;
        // region_populations_current[human_home_region[hid]]++;
        // human_current_region[hid] = human_home_region[hid];
        send_human_to_region(hid, human_home_region[hid]);
    }

    // modifying getters
    inline void clear_logbook_human_movement(){
        logbook_human_movement.clear();
    }

    // const getters
    inline const std::vector<std::string>& get_logbook_human_movement() const {
        return logbook_human_movement;
    }
    inline const std::vector<int>& get_region_populations_current() const {
        return region_populations_current;
    }
    inline const std::vector<int>& get_region_i_t() const {
        return region_i_t;
    }
    inline const std::vector<double>& get_region_prev_t() const {
        return region_prev_t;
    }
    inline const std::vector<int>& get_human_current_region() const {
        return human_current_region;
    }

    inline int get_no_of_regions() const {
        return kNo_of_regions;
    }
    inline int get_total_population() const {
        return static_cast<int>(human_home_region.size());
    }
    inline int get_human_current_region(int human_id) const {
        return human_current_region[human_id];
    }

    inline int get_human_home_region(int human_id) const {
        return human_home_region[human_id];
    }

};
#endif