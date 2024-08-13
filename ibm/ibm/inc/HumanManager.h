#ifndef HUMAN_MAN_H
#define HUMAN_MAN_H

class HumanManager {
    int tt = 0;
    const int kNo_of_humans;
    const int kIntervention_start;
    const int kIntervention_stop;
    const double kIntervention_phase;
    double coverage = 0.0;
    const int kACF_screened;
    int HHS;
    const double kSensitivity;
    std::vector<int> phi_destination;
    std::vector<int> tau_I;
    //std::vector<int> tau_L;
    //std::vector<int> eta;
    std::vector<int> omega_vec;
    int sampled_hh = 0;

    const double kBeta;
    const double kGamma;
    const double kHh_attack_rate;
    const double kAlpha;
    const double kPhi;
    double kDelta1;
    double kDelta2;
    
    const double kDelta1_intervention;
    const double kDelta2_intervention;

    const double kTau_I_mu;
    const double kTau_I_sd;
    const int kTau_I_min;

    /*
    const double kTau_L_mu;
    const double kTau_L_sd;
    const int kTau_L_min;
    */
   const double kOmega;
   const double kEta;
   const double kNu;
   const double kMu;
    
    /*
    const double kEta_mu;
    const double kEta_sd;
    const int kEta_min;
    */
    // double lambda_t = 0.0;
    int s_t = 0;
    int e_t = 0;
    int i_t = 0;
    int l_t = 0;
    int a_t = 0;
    int r_t = 0;
    int d_t = 0;
    int inc_t = 0;
    int hh_inc_t = 0;

    //std::vector<bool> sex;
    //std::vector<int> age_group;
    //std::vector<double> kAge_weights;
    std::vector<int> household;
    std::vector<int> was_latent_status;
    std::vector<int> time_to_hh_attack;

    int kAvg_household_size;
    int kNo_of_households;

    std::vector<int> infectious_households;
    
    std::vector<HumanSeirStatus> human_seir_status; //pop
public:
    HumanManager(
        int no_of_humans,
        double beta,
        double gamma,
        double hh_attack_rate,
        double alpha,
        double phi,
        double delta1,
        double delta2,

        double tau_I_mu,
        double tau_I_sd,
        int tau_I_min,
        /*
        double tau_L_mu,
        double tau_L_sd,
        int tau_L_min,
        */
        double omega,
        double eta,
        double nu,
        double mu,
        /*
        double eta_mu,
        double eta_sd,
        int eta_min,
        */
        
        int init_e,
        int init_i,
        int init_l,
        int init_a,
        int init_r,

        int no_of_households,
        //const std::vector<double>& age_weights,
        std::vector<int> input_household,
        int intervention_start,
        int intervention_stop,
        double intervention_phase,
        double delta1_intervention,
        double delta2_intervention,
        int acf_screened,
        int HHS,
        double sensitivity
        );
    ~HumanManager();
    void step(
        // const std::vector<int>& region_populations_current,
        // const std::vector<int>& region_i_t,
        //const std::vector<double>& region_prev_t,
        //const std::vector<int>& human_current_region
    );

    inline const std::vector<HumanSeirStatus>& get_human_seir_status() const {
        return human_seir_status;
    }

    inline int get_s_t() const {return s_t;}
    inline int get_e_t() const {return e_t;}
    inline int get_i_t() const {return i_t;}
    inline int get_l_t() const {return l_t;}
    inline int get_a_t() const {return a_t;}
    inline int get_r_t() const {return r_t;}
    inline int get_d_t() const {return d_t;}
    inline int get_inc_t() const {return inc_t;}
    inline int get_hh_inc_t() const {return hh_inc_t;}

    inline int get_no_of_humans() const {return kNo_of_humans;}
    inline int get_no_of_households() const {return kNo_of_households;}
    inline std::vector<int> get_infectious_households() const {return infectious_households;}

    };
#endif