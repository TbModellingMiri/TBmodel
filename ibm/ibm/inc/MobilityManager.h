#ifndef MOBILITY_MAN_H
#define MOBILITY_MAN_H

class MobilityManager {
    
    std::vector<int> human_days_until_return_home;
    
    std::mt19937 gen{std::random_device{}()};
    std::uniform_int_distribution<> uni_int_dis;
    std::vector<std::vector<double>> kGravity;
    int kMobility_model;
    std::time_t current_time;
    std::time_t holiday_start;
    std::time_t holiday_end;
    int kRecurrent;
        
    // sharable functions
    int pick_region_uniform(int home_region);
    int pick_region_gravity(int home_region);

    // type 1 : short-term return-home travel
    const double kProb_travel;
    const double kProb_travel_season;
    // const double kDays_away_mean;
    // const double kDays_away_sd;
    std::normal_distribution<double> nrm_dis_type1;
    std::normal_distribution<double> nrm_dis_season;

    // type 1 functions
    int type1_gen_travel_length();
    int season_gen_travel_length();
    
    void mobility_season(
    RegionManager& rgn_mgr
    );
    
    void mobility_type1(
        RegionManager& rgn_mgr
    );

public:
    MobilityManager(
        int no_of_humans,
        int no_of_regions,
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
        const int holiday_start_m,
        const int holiday_start_d,
        const int holiday_start_y,
        const int holiday_length,
        const int recurrent

        );
    ~MobilityManager();

    void step( // friend of RegionManager
        // std::vector<int>& human_home_region,
        // std::vector<int>& human_current_region,
        // std::vector<int>& region_populations_home,
        // std::vector<int>& region_populations_current
        RegionManager& rgn_mgr
    );

};

#endif