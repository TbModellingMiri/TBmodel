#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <string>
#include <chrono>
#include <ctime>

#include <filesystem> // std::filesystem::path, std::filesystem::exists

#include "inc/io.h"
#include "inc/common.h"
#include "inc/RegionManager.h"
#include "inc/HumanManager.h"
#include "inc/MobilityManager.h"

#define LOG(x) {std::cout << x << std::endl;}

void print_help() {
    std::cout << "Usage: command_line.exe -c config.json -o output_dir [-p prefix]"
                << "\n -c\t JSON formatted file containing simulation parameters"
                << "\n -o\t name of output directory"
                << "\n -p\t (optional) prefix of output, default = timestamp"
                << "\n";
}

int main(int argc, char** argv) {
    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "start: ok" << std::endl;
    ////
    // Inputs

    ///
    // Command line inputs
    InputParser input_parser(argc, argv);

    if (input_parser.cmdOptionExists("--help")) {
        print_help();
    }
    bool input_check_passed = true;
    const std::string& config_file_name = input_parser.getCmdOption("-c");
    if (config_file_name.empty()) {
        std::cout << "Error: missing config file. See `-c` option.\n";
        input_check_passed = false;
    }
    const std::string& output_dir_name = input_parser.getCmdOption("-o");
    if (output_dir_name.empty()) {
        std::cout << "Error: missing output directory name. See `-o` option.\n";
        input_check_passed = false;
    }
    if (!input_check_passed) {
        print_help();
        return EXIT_FAILURE;
    }
    if (!std::filesystem::exists(output_dir_name)) {
        std::filesystem::create_directory(output_dir_name);
        std::cout << "Created " << output_dir_name << " for outputs\n";
    }

    std::string output_prefix = input_parser.getCmdOption("-p");
    if (output_prefix.empty()) {
        output_prefix = get_time_stamp();
    }

    std::cout << "ibm.exe started with:"
              << "\n\tconfig :\t\t" << config_file_name
              << "\n\toutput directory : \t" << output_dir_name
              << "\n\toutput prefix : \t" << output_prefix
              << "\n";

    const std::filesystem::path kOut_path = std::filesystem::path(".")
                                        / std::filesystem::path(output_dir_name)
                                        / std::filesystem::path(output_prefix);

    ///
    // Parse JSON -> config
    rapidjson::Document config = get_json_from_file(config_file_name);
    const int kNo_of_timesteps = config["simulation"]["length_days"].GetInt();

    ///
    // Parse CSV
    const std::filesystem::path region_file = std::filesystem::path(".")
        / std::filesystem::path(config["data"]["path"].GetString())
        / std::filesystem::path(config["data"]["region"]["file"].GetString());

    rapidcsv::Document region_csv_doc(region_file.string());
    std::vector<int> region_id_col = region_csv_doc.GetColumn<int>("reg_id");
    std::vector<int> region_pops_col = region_csv_doc.GetColumn<int>("population");
    std::vector<double> region_lon_col = region_csv_doc.GetColumn<double>("lon");
    std::vector<double> region_lat_col = region_csv_doc.GetColumn<double>("lat");

    const std::filesystem::path age_dist_file = std::filesystem::path(".")
        / std::filesystem::path(config["data"]["path"].GetString())
        / std::filesystem::path(config["data"]["age_dist"]["file"].GetString());

    rapidcsv::Document age_dist_csv_doc(age_dist_file.string());
    std::vector<double> age_weights = age_dist_csv_doc.GetColumn<double>("population");

    const std::filesystem::path household_file = std::filesystem::path(".")
        / std::filesystem::path(config["data"]["path"].GetString())
        / std::filesystem::path(config["data"]["household_dist"]["file"].GetString());

    rapidcsv::Document household_dist_doc(household_file.string());
    std::vector<int> household = household_dist_doc.GetColumn<int>("household");
    
    if (config["data"]["region"]["no_of_regions"].GetInt() > static_cast<int>(region_pops_col.size())) {
        std::cout << "Error: no_of_regions is greater than the number of rows in region file.\n";
        return EXIT_FAILURE;
    }
    
    // inputs parsed.
    ////

    ////
    // Construct data managers
    RegionManager* rgn_mgr (
        new RegionManager(
            config["data"]["region"]["no_of_regions"].GetInt(),
            region_pops_col,
            config["initial_condition"]["region"]["region"].GetInt(),
            config["initial_condition"]["exposed"]["size"].GetInt()
        )
    );
    const int kTotal_population = rgn_mgr->get_total_population();
    const int kNo_of_regions = config["data"]["region"]["no_of_regions"].GetInt();
    
    /*
    AttributeManager* atr_mgr(
        new AttributeManager(
            kTotal_population,
            avg_household_size,
            age_weights
        )
    );
    */

    HumanManager* hmn_mgr (
        new HumanManager(
            kTotal_population,
            config["parameters"]["beta"].GetDouble(),
            config["parameters"]["gamma"].GetDouble(),
            config["parameters"]["hh_attack_rate"].GetDouble(),
            config["parameters"]["alpha"].GetDouble(),
            config["parameters"]["phi"].GetDouble(),
            config["parameters"]["delta1"].GetDouble(),
            config["parameters"]["delta2"].GetDouble(),

            config["parameters"]["tau_I_mu"].GetDouble(),
            config["parameters"]["tau_I_sd"].GetDouble(),
            config["parameters"]["tau_I_min"].GetInt(),
            /*
            config["parameters"]["tau_L_mu"].GetDouble(),
            config["parameters"]["tau_L_sd"].GetDouble(),
            config["parameters"]["tau_L_min"].GetInt(),
        
            config["parameters"]["omega_mu"].GetDouble(),
            config["parameters"]["omega_sd"].GetDouble(),
            config["parameters"]["omega_min"].GetInt(),
            
            config["parameters"]["eta_mu"].GetDouble(),
            config["parameters"]["eta_sd"].GetDouble(),
            config["parameters"]["eta_min"].GetInt(),
            */
            config["parameters"]["omega"].GetDouble(),
            config["parameters"]["eta"].GetDouble(),
            config["parameters"]["nu"].GetDouble(),
            config["parameters"]["mu"].GetDouble(),

            config["initial_condition"]["exposed"]["size"].GetInt(),
            config["initial_condition"]["infectious"]["size"].GetInt(),
            config["initial_condition"]["latent"]["size"].GetInt(),
            config["initial_condition"]["treated"]["size"].GetInt(),
            config["initial_condition"]["recovered"]["size"].GetInt(),

            config["parameters"]["no_of_households"].GetInt(),
            //age_weights,
            household,
            config["intervention"]["intervention_start"].GetInt(),
            config["intervention"]["intervention_stop"].GetInt(),
            config["intervention"]["intervention_phase"].GetDouble(),
            config["intervention"]["delta1_intervention"].GetDouble(),
            config["intervention"]["delta2_intervention"].GetDouble(),
            config["intervention"]["acf_screened"].GetInt(),
            config["intervention"]["HHS"].GetInt(),
            config["intervention"]["sensitivity"].GetDouble()
        )
    );

    rgn_mgr->update_region_seir_t(hmn_mgr->get_human_seir_status());

    std::cout << "Human 0 assigned to region " << rgn_mgr->get_human_current_region(0) << "\n";

    MobilityManager mbl_mgr (
        kTotal_population,
        kNo_of_regions,
        config["mobility"]["type_1"]["prob_travel_times_per_year"].GetDouble()
            / kNo_of_days_in_year,
        config["mobility"]["type_1"]["days_away_mean"].GetDouble(),
        config["mobility"]["type_1"]["days_away_sd"].GetDouble(),

        config["mobility"]["season"]["prob_travel_times_per_year"].GetDouble()
            / kNo_of_days_in_year,
        config["mobility"]["season"]["days_away_mean"].GetDouble(),
        config["mobility"]["season"]["days_away_sd"].GetDouble(),

        region_pops_col,
        region_lon_col, 
        region_lat_col,


        config["mobility"]["mobility_model"].GetInt(),
        config["simulation"]["year_start"].GetInt(),
        config["simulation"]["month_start"].GetInt(),
        config["simulation"]["day_start"].GetInt(),
        config["mobility"]["season"]["holiday_start_year"].GetInt(),
        config["mobility"]["season"]["holiday_start_month"].GetInt(),
        config["mobility"]["season"]["holiday_start_day"].GetInt(),
        config["mobility"]["season"]["length"].GetInt(),
        config["mobility"]["season"]["recurrent"].GetInt()
    );

    ////
    // Create logging / output data containers
    //std::vector<std::vector<int>> sim_table_infectious_households(kNo_of_timesteps, std::vector<int>(config["parameters"]["no_of_households"].GetInt()));
    std::vector<float> sim_table_day(kNo_of_timesteps, 0.0);
    std::vector<float> sim_table_s(kNo_of_timesteps, 0.0);
    std::vector<float> sim_table_e(kNo_of_timesteps, 0.0);
    std::vector<float> sim_table_i(kNo_of_timesteps, 0.0);
    std::vector<float> sim_table_l(kNo_of_timesteps, 0.0);
    std::vector<float> sim_table_a(kNo_of_timesteps, 0.0);
    std::vector<float> sim_table_r(kNo_of_timesteps, 0.0);
    std::vector<float> sim_table_d(kNo_of_timesteps, 0.0);
    std::vector<float> sim_table_inc(kNo_of_timesteps, 0.0);
    std::vector<float> sim_table_hh_inc(kNo_of_timesteps, 0.0);
    std::vector<double> rpt_prev(kNo_of_timesteps, -1.0);
    std::vector<std::vector<double>> rpt_prev_regions(
        kNo_of_regions,
        std::vector<double>(kNo_of_timesteps, -1.0)
    );

    std::vector<std::string> rpt_mobility;


    ////
    // start of simulation
    // Get settings for model expansions
    int mobility_model_setting = config["mobility"]["mobility_model"].GetInt();
    std::string mobility_model;
    if(mobility_model_setting == 2) {
        mobility_model = "Gravity model";
    } else {
        mobility_model = "Random uniform";
        }

     std::string start_date =
     std::to_string(config["simulation"]["year_start"].GetInt()) + "-" +
     std::to_string(config["simulation"]["month_start"].GetInt()) + "-" +
     std::to_string(config["simulation"]["day_start"].GetInt());
    
    std::cout << "set-up completed:\n"
                << "\tpopulation: \t"       << kTotal_population << "\n"
                << "\tregions: \t"          << kNo_of_regions    << "\n"
                << "\ttime: \t"             << kNo_of_timesteps  << "\n"
                << "\tmobility model: \t"   << mobility_model    << "\n"
                << "\tstart date: \t"       << start_date        << "\n"
                << "\n";

    for (int tt = 0; tt < kNo_of_timesteps; tt++) {
        ////
        // Logging for outputs
        //sim_table_infectious_households[tt] = hmn_mgr->get_infectious_households();
        sim_table_day[tt] = tt;
    
        sim_table_s[tt] = hmn_mgr->get_s_t();
        sim_table_e[tt] = hmn_mgr->get_e_t();
        sim_table_i[tt] = hmn_mgr->get_i_t();
        sim_table_l[tt] = hmn_mgr->get_l_t();
        sim_table_a[tt] = hmn_mgr->get_a_t();
        sim_table_r[tt] = hmn_mgr->get_r_t();
        sim_table_d[tt] = hmn_mgr->get_d_t();
        sim_table_inc[tt] = hmn_mgr->get_inc_t();
        sim_table_hh_inc[tt] = hmn_mgr->get_hh_inc_t();
        rpt_prev[tt] = static_cast<double>(hmn_mgr->get_i_t()) / static_cast<double>(kTotal_population);
        int rr = 0;
        for (const double& pp : rgn_mgr->get_region_prev_t()) {
            rpt_prev_regions[rr++][tt] = pp;
        }
        ////
        // start of time step

        // mobility
        
        // mbl_mgr.step(*rgn_mgr);
        // rgn_mgr->update_region_seir_t(hmn_mgr->get_human_seir_status());
        // s->e->i->r->s
        hmn_mgr->step(
            // rgn_mgr->get_region_populations_current(),
            // rgn_mgr->get_region_i_t(),
            // rgn_mgr->get_region_prev_t(),
            // rgn_mgr->get_human_current_region()
        );
        // rgn_mgr->update_region_seir_t(hmn_mgr->get_human_seir_status());
       
        ////
        // end of day report
        /*
        for ( const std::string& ss: rgn_mgr->get_logbook_human_movement()) {
            rpt_mobility.push_back(
                std::to_string(tt)
                + kCsv_sep_str
                + ss
            );
        }
        rgn_mgr->clear_logbook_human_movement();
        */
        
    }

    // outputs

    std::vector<std::vector<float>> sim_table {
        sim_table_day,
        sim_table_s,
        sim_table_e,
        sim_table_i,
        sim_table_l,
        sim_table_a,
        sim_table_r,
        sim_table_d,
        sim_table_inc,
        sim_table_hh_inc
    };
    std::vector<std::string> sim_table_headers {
        "day", "susceptible", "exposed", "infected", "latent", "treated", "recovered", "diseased", "inc", "hh_inc"
    };
    // print_vov(sim_table, sim_table_headers);
    // print_vov(sim_table);

    /*
    std::vector<int> age_group = atr_mgr->get_age_group();
    std::vector<int> households = atr_mgr->get_household();
    std::vector<std::string> headers = {"age","household","prevalence"};

    write_to_csv(kOut_path.string() + "_AGE_GROUPS.csv", age_group, headers[0]);
    write_to_csv(kOut_path.string() + "_HOUSEHOLDS.csv", households, headers[1]);
    
    */

    write_to_csv(kOut_path.string() + "_sim_table.csv", sim_table, sim_table_headers);
    //write_to_csv(kOut_path.string() + "_prev.csv", rpt_prev, "prevalence");
    //write_to_csv(kOut_path.string() + "_prev_regions.csv", rpt_prev_regions);
    //write_to_csv(kOut_path.string() + "_human_mobility.csv", rpt_mobility);
    //write_to_csv(kOut_path.string() + "_infectious_households.csv", sim_table_infectious_households);

    // print_v(rpt_mobility);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "exit: ok" << std::endl;
    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;
    return EXIT_SUCCESS;
}