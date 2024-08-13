#include <iostream>
#include <vector>
#include <algorithm> // std::count
#include <cassert>
#include <cmath>
#include "common.h"

#include "HumanManager.h"

HumanManager::HumanManager(
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
  
    int init_e,
    int init_i,
    int init_l,
    int init_a,
    int init_r,

    int no_of_households,

    std::vector<int> input_household,
    int intervention_start,
    int intervention_stop,
    double intervention_phase,
    double delta1_intervention,
    double delta2_intervention,
    int acf_screened,
    int HHS,
    double sensitivity

    ) : kNo_of_humans(no_of_humans),
        human_seir_status(kNo_of_humans, HumanSeirStatus::SUCP),

        phi_destination(kNo_of_humans),
        tau_I(kNo_of_humans),

        omega_vec(kNo_of_humans, kOmega),
        kBeta(beta),
        kGamma(gamma),
        kHh_attack_rate(hh_attack_rate),
        kAlpha(alpha),
        kPhi(phi),
        kDelta1(delta1),
        kDelta2(delta2),
        kTau_I_mu(tau_I_mu),
        kTau_I_sd(tau_I_sd),
        kTau_I_min(tau_I_min),
       
        kOmega(omega),
        kEta(eta),
        kNu(nu),
        kMu(mu),
        
        infectious_households(kNo_of_households, 0), // and ensure size > kHuman number or greater than csv
        was_latent_status(kNo_of_humans, 1),
        time_to_hh_attack(kNo_of_humans, -1),
        household(input_household),
        kNo_of_households(household.size()),
        kIntervention_start(intervention_start),
        kIntervention_stop(intervention_stop),
        kIntervention_phase(intervention_phase),
        kDelta1_intervention(delta1_intervention),
        kDelta2_intervention(delta2_intervention),
        coverage(0),
        kACF_screened(acf_screened),
        HHS(HHS),
        kSensitivity(sensitivity)
{

    for (int ee = 0; ee < init_e; ee++)
    {
        human_seir_status[ee] = HumanSeirStatus::EXPD;
    }
    for (int ii = 0; ii < init_i; ii++)
    {
        human_seir_status[init_e + ii] = HumanSeirStatus::INFT;
    }
    for (int ll = 0; ll < init_l; ll++)
    {
        human_seir_status[init_e + init_i + ll] = HumanSeirStatus::LATN;
    }
    for (int aa = 0; aa < init_a; aa++)
    {
        human_seir_status[init_e + init_i + init_l + aa] = HumanSeirStatus::TRTD;
    }
    for (int rr = 0; rr < init_r; rr++)
    {
        human_seir_status[init_e + init_i + init_l + init_a + rr] = HumanSeirStatus::RCVD;
    }

    s_t = std::count(human_seir_status.begin(), human_seir_status.end(), HumanSeirStatus::SUCP);
    e_t = std::count(human_seir_status.begin(), human_seir_status.end(), HumanSeirStatus::EXPD);
    i_t = std::count(human_seir_status.begin(), human_seir_status.end(), HumanSeirStatus::INFT);
    l_t = std::count(human_seir_status.begin(), human_seir_status.end(), HumanSeirStatus::LATN);
    a_t = std::count(human_seir_status.begin(), human_seir_status.end(), HumanSeirStatus::TRTD);
    r_t = std::count(human_seir_status.begin(), human_seir_status.end(), HumanSeirStatus::RCVD);
    d_t = std::count(human_seir_status.begin(), human_seir_status.end(), HumanSeirStatus::DISD);

    // do_shuffle(human_seir_status);

    assert(e_t == init_e);
    assert(i_t == init_i);
    assert(l_t == init_l);
    assert(a_t == init_a);
    assert(r_t == init_r);
    assert(s_t + e_t + i_t + l_t + a_t + r_t == kNo_of_humans);
  

    // assign attributes
    for (int hh = 0; hh < kNo_of_humans; hh++)
    {
        phi_destination[hh] = get_rand_uniform() > kPhi ? 1 : 0;
        tau_I[hh] = std::max(kTau_I_min, static_cast<int>(std::round(get_rand_normal(kTau_I_mu, kTau_I_sd))));
       
        // add infections to households
        if (human_seir_status[hh] == HumanSeirStatus::INFT)
        {
            infectious_households[household[hh]]++;
        }
    }
    assert(time_to_hh_attack.size() == kNo_of_humans);
    assert(was_latent_status.size() == kNo_of_humans);
    assert(household.size() == kNo_of_humans);
    assert(infectious_households.size() == kNo_of_households);
}

HumanManager::~HumanManager()
{
    std::cout << "~HumanManager() called\n";
}

void HumanManager::step(
    // const std::vector<int>& region_populations_current,
    // const std::vector<int>& region_i_t,
    // const std::vector<double> &region_prev_t,
    // const std::vector<int> &human_current_region
)
{
    double pop = kNo_of_humans - get_d_t();
    double my_prev = get_i_t() / pop;
    inc_t = 0;
    hh_inc_t = 0;
    if (tt >= kIntervention_start)
    {
        kDelta1 = kDelta1_intervention;
        kDelta2 = kDelta2_intervention;
    }

    if (tt >= kIntervention_start && tt < kIntervention_stop && sampled_hh <= kNo_of_humans )
    {
        int samples = std::min(kACF_screened, static_cast<int>(std::round((tt-kIntervention_start)*(kACF_screened/kIntervention_phase))));
        for (int j = 0; j < samples; j++)
        {
            //int sampled_hh = sample_integer(kNo_of_humans);
            if(sampled_hh >= kNo_of_humans)
            {
                break;
            }
            if(human_seir_status[sampled_hh] == HumanSeirStatus::DISD)
            {
                j--;
            }
            else if ((human_seir_status[sampled_hh] == HumanSeirStatus::EXPD ||
            human_seir_status[sampled_hh] == HumanSeirStatus::INFT ||
            human_seir_status[sampled_hh] == HumanSeirStatus::LATN) &&
            get_rand_uniform() <kSensitivity)
            {
                if (human_seir_status[sampled_hh] == HumanSeirStatus::INFT)
                {
                    was_latent_status[sampled_hh] = 0;
                    omega_vec[sampled_hh] = kOmega;
                    time_to_hh_attack[sampled_hh] = -1;
                    infectious_households[household[sampled_hh]]--;
                }
                else
                {
                    human_seir_status[sampled_hh] = HumanSeirStatus::TRTD;
                    was_latent_status[sampled_hh] = 1;
                    time_to_hh_attack[sampled_hh] = -1; 
                    omega_vec[sampled_hh] = kOmega;
                }
                // !! intervention
                if (HHS == 1)
                {
                    std::vector<int> members = get_household_members(household, household[sampled_hh]);

                    for (int i = 0; i < members.size(); i++)
                    {
                        if ((human_seir_status[members[i]] == HumanSeirStatus::EXPD ||
                             human_seir_status[members[i]] == HumanSeirStatus::INFT ||
                             human_seir_status[members[i]] == HumanSeirStatus::LATN) &&
                             get_rand_uniform() < kSensitivity)
                        {
                            // time_to_hh_attack[members[i]] = static_cast<int>(std::round(get_rand_lnorm(log(30), log(1.8))));
                            if (human_seir_status[members[i]] == HumanSeirStatus::INFT)
                            {
                                was_latent_status[members[i]] = 0;
                                infectious_households[household[members[i]]]--;
                            }
                            else
                            {
                                was_latent_status[members[i]] = 1;
                            }
                            human_seir_status[members[i]] = HumanSeirStatus::TRTD;
                            omega_vec[members[i]] = kOmega;
                            time_to_hh_attack[members[i]] = -1; // empty hh attack timer
                        }
                    }
                }
            }
            sampled_hh++;
        }
    }

    for (size_t hh = 0; hh < kNo_of_humans; hh++)
    {
        if (human_seir_status[hh] == HumanSeirStatus::SUCP)
        {
            if (get_rand_uniform() < 1 / (74.8 * 365))
            {
                // human_seir_status[hh] = HumanSeirStatus::DISD;
                human_seir_status[hh] = HumanSeirStatus::SUCP;
                time_to_hh_attack[hh] = -1;
                household[hh] = sample_integer(household.size());
            }
            else
            {
                // int hh_region = human_current_region[hh];
                // double lambda_hh_t = kBeta
                //             * static_cast<double>(region_i_t[hh_region])
                //             / static_cast<double>(region_populations_current[hh_region]);
                // double lambda_hh_t = kBeta * region_prev_t[hh_region];
                double lambda_hh_t = kBeta * my_prev;
                // insert
                // end insert
                if (get_rand_uniform() < lambda_hh_t)
                {
                    human_seir_status[hh] = HumanSeirStatus::EXPD;
                    phi_destination[hh] = get_rand_uniform() > kPhi ? 1 : 0;
                    time_to_hh_attack[hh] = -1; // empty hh attack timer
                }
                // }
                else if (time_to_hh_attack[hh] >= 0) // check if hh exposed
                {
                    if (time_to_hh_attack[hh] == 0) // check if time to transition
                    {
                        human_seir_status[hh] = HumanSeirStatus::EXPD;
                        phi_destination[hh] = get_rand_uniform() > kPhi ? 1 : 0;
                        time_to_hh_attack[hh] = -1; // empty hh attack timer
                        if (phi_destination[hh] == 0)
                        {
                            hh_inc_t++;
                        }
                    }
                    else
                    {
                        time_to_hh_attack[hh]--;
                    }
                }
            }
        }
        else if (human_seir_status[hh] == HumanSeirStatus::EXPD)
        {
            if (get_rand_uniform() < 1 / (74.8 * 365))
            {
                // human_seir_status[hh] = HumanSeirStatus::DISD;
                human_seir_status[hh] = HumanSeirStatus::SUCP;
                time_to_hh_attack[hh] = -1;
                household[hh] = sample_integer(household.size());
            }
            else
            {
                if (get_rand_uniform() < kGamma)
                {
                    if (phi_destination[hh] == 1)
                    {
                        human_seir_status[hh] = HumanSeirStatus::LATN;
                    }
                    else
                    {
                        human_seir_status[hh] = HumanSeirStatus::INFT;
                        tau_I[hh] = std::max(kTau_I_min, static_cast<int>(std::round(get_rand_normal(kTau_I_mu, kTau_I_sd))));
                        infectious_households[household[hh]]++;
                        inc_t++;

                        time_to_hh_attack[hh] = -1; // empty hh attack timer
                        if (infectious_households[household[hh]] == 1)
                        {
                            std::vector<int> members = get_household_members(household, household[hh]);
                            for (int i = 0; i < members.size(); i++)
                            {
                                if (human_seir_status[members[i]] != HumanSeirStatus::TRTD && human_seir_status[members[i]] != HumanSeirStatus::INFT && get_rand_uniform() < kHh_attack_rate)
                                {
                                    // time_to_hh_attack[members[i]] = static_cast<int>(std::round(get_rand_lnorm(log(30), log(1.8))));
                                    time_to_hh_attack[members[i]] = static_cast<int>(std::round(get_rand_gamma(4, 0.2)));
                                }
                            }
                        }
                    }
                }
                // insert
                if (time_to_hh_attack[hh] >= 0)
                {
                    if (time_to_hh_attack[hh] == 0)
                    {
                        human_seir_status[hh] = HumanSeirStatus::EXPD;
                        phi_destination[hh] = get_rand_uniform() > kPhi ? 1 : 0;
                        time_to_hh_attack[hh]--;
                        if (phi_destination[hh] == 0)
                        {
                            hh_inc_t++;
                        }
                    }
                    else
                    {
                        time_to_hh_attack[hh]--;
                    }
                }
                // end insert
            }
        }
        else if (human_seir_status[hh] == HumanSeirStatus::INFT)
        {
            if (get_rand_uniform() < 1 / (74.8 * 365))
            {
                // human_seir_status[hh] = HumanSeirStatus::DISD;
                infectious_households[household[hh]]--;
                human_seir_status[hh] = HumanSeirStatus::SUCP;
                time_to_hh_attack[hh] = -1;
                household[hh] = sample_integer(household.size());
            }
            else
            {
                if (get_rand_uniform() < kMu)
                {
                    human_seir_status[hh] = HumanSeirStatus::DISD;
                    infectious_households[household[hh]]--;
                    time_to_hh_attack[hh] = -1; // empty hh attack timer
                }
                else if (tau_I[hh] == 0)
                {
                    human_seir_status[hh] = HumanSeirStatus::RCVD;
                    infectious_households[household[hh]]--;
                }
                else if (get_rand_uniform() < kDelta1 + kDelta2) 
                {
                    human_seir_status[hh] = HumanSeirStatus::TRTD;
                    was_latent_status[hh] = 0;
                    omega_vec[hh] = kOmega;
                    time_to_hh_attack[hh] = -1; // empty hh attack timer
                    // DECREMENT INFECTIOUS HOUSEHOLD

                    infectious_households[household[hh]]--;
                    // !! intervention
                    if (tt >= kIntervention_start && HHS == 1)
                    {
                        if (tt <= kIntervention_start + kIntervention_phase)
                        {
                            coverage = (tt - kIntervention_start) / kIntervention_phase;
                            coverage = coverage > get_rand_uniform() ? 1.0 : 0.0;
                        }
                        else
                        {
                            coverage = 1.0;
                        }

                        if (coverage == 1.0)
                        {
                            std::vector<int> members = get_household_members(household, household[hh]);

                            for (int i = 0; i < members.size(); i++)
                            {
                                if ((human_seir_status[members[i]] == HumanSeirStatus::EXPD ||
                                     human_seir_status[members[i]] == HumanSeirStatus::INFT ||
                                     human_seir_status[members[i]] == HumanSeirStatus::LATN) &&
                                     get_rand_uniform() < kSensitivity)
                                {
                                    // time_to_hh_attack[members[i]] = static_cast<int>(std::round(get_rand_lnorm(log(30), log(1.8))));
                                    if (human_seir_status[members[i]] == HumanSeirStatus::INFT)
                                    {
                                        was_latent_status[members[i]] = 0;
                                        infectious_households[household[members[i]]]--;
                                    }
                                    else
                                    {
                                        was_latent_status[members[i]] = 1;
                                    }
                                    human_seir_status[members[i]] = HumanSeirStatus::TRTD;
                                    omega_vec[members[i]] = kOmega;
                                    time_to_hh_attack[members[i]] = -1; // empty hh attack timer
                                }
                            }
                        }
                    }
                    // !! intervention end
                }
                tau_I[hh] = tau_I[hh] - 1;
            }
        }
        else if (human_seir_status[hh] == HumanSeirStatus::LATN)
        {
            if (get_rand_uniform() < 1 / (74.8 * 365))
            {
                // human_seir_status[hh] = HumanSeirStatus::DISD;
                human_seir_status[hh] = HumanSeirStatus::SUCP;
                time_to_hh_attack[hh] = -1;
                household[hh] = sample_integer(household.size());
            }

            else
            {
                // int hh_region = human_current_region[hh];
                // double lambda_hh_t = kBeta * region_prev_t[hh_region];
                double lambda_hh_t = kBeta * my_prev;

                if (get_rand_uniform() < kEta)
                {
                    // else if(eta[hh] == 0) {
                    human_seir_status[hh] = HumanSeirStatus::INFT;
                    tau_I[hh] = std::max(kTau_I_min, static_cast<int>(std::round(get_rand_normal(kTau_I_mu, kTau_I_sd))));
                    infectious_households[household[hh]]++;
                    inc_t++;

                    if (infectious_households[household[hh]] == 1)
                    {
                        std::vector<int> members = get_household_members(household, household[hh]);
                        for (int i = 0; i < members.size(); i++)
                        {
                            if (human_seir_status[members[i]] != HumanSeirStatus::TRTD && human_seir_status[members[i]] != HumanSeirStatus::INFT && get_rand_uniform() < kHh_attack_rate)
                            {
                                // time_to_hh_attack[members[i]] = static_cast<int>(std::round(get_rand_lnorm(log(30), log(1.8))));
                                time_to_hh_attack[members[i]] = static_cast<int>(std::round(get_rand_gamma(4, 0.2)));
                            }
                        }
                    }
                }
                else if (get_rand_uniform() < kDelta1)
                {
                    human_seir_status[hh] = HumanSeirStatus::TRTD;
                    was_latent_status[hh] = 1;
                    time_to_hh_attack[hh] = -1; // empty hh attack
                    omega_vec[hh] = kOmega;
                    // !! intervention
                    if (tt >= kIntervention_start && HHS == 1)
                    {
                        if (tt <= kIntervention_start + kIntervention_phase)
                        {
                            coverage = (tt - kIntervention_start) / kIntervention_phase;
                            coverage = coverage > get_rand_uniform() ? 1.0 : 0.0;
                        }
                        else
                        {
                            coverage = 1.0;
                        }

                        if (coverage == 1.0)
                        {
                            std::vector<int> members = get_household_members(household, household[hh]);
                            for (int i = 0; i < members.size(); i++)
                            {
                                if ((human_seir_status[members[i]] == HumanSeirStatus::EXPD ||
                                     human_seir_status[members[i]] == HumanSeirStatus::INFT ||
                                     human_seir_status[members[i]] == HumanSeirStatus::LATN) &&
                                     get_rand_uniform() < kSensitivity)
                                {
                                    // time_to_hh_attack[members[i]] = static_cast<int>(std::round(get_rand_lnorm(log(30), log(1.8))));
                                    if (human_seir_status[members[i]] == HumanSeirStatus::INFT)
                                    {
                                        was_latent_status[members[i]] = 0;
                                        infectious_households[household[members[i]]]--;
                                    }
                                    else
                                    {
                                        was_latent_status[members[i]] = 1;
                                    }
                                    human_seir_status[members[i]] = HumanSeirStatus::TRTD;
                                    omega_vec[members[i]] = kOmega;
                                    time_to_hh_attack[members[i]] = -1;
                                }
                            }
                        }
                    }
                    // !! intervention end
                }
                // else if(infectious_households[household[hh]]>0) {
                //         int household_lambda = 0.355598/573.05;
                //         if(get_rand_uniform() < household_lambda) {
                //             human_seir_status[hh] = HumanSeirStatus::EXPD;
                //             phi_destination[hh] = false;
                //         }
                // }
                else if (get_rand_uniform() < (kNu * lambda_hh_t))
                {
                    human_seir_status[hh] = HumanSeirStatus::EXPD;
                    phi_destination[hh] = get_rand_uniform() > kPhi ? 1 : 0;
                }

                // tau_L[hh] = tau_L[hh]-1;
                // eta[hh] = eta[hh]-1;
                // insert
                else if (time_to_hh_attack[hh] >= 0)
                {
                    if (time_to_hh_attack[hh] == 0)
                    {
                        human_seir_status[hh] = HumanSeirStatus::EXPD;
                        phi_destination[hh] = get_rand_uniform() > kPhi ? 1 : 0;
                        time_to_hh_attack[hh] = -1; // empty hh attack timer
                        if (phi_destination[hh] == 0)
                        {
                            hh_inc_t++;
                        }
                    }
                    else
                    {
                        time_to_hh_attack[hh]--;
                    }
                }
                // end insert
            }
        }
        else if (human_seir_status[hh] == HumanSeirStatus::TRTD)
        {
            if (get_rand_uniform() < 1 / (74.8 * 365))
            {
                // human_seir_status[hh] = HumanSeirStatus::DISD;
                human_seir_status[hh] = HumanSeirStatus::SUCP;
                time_to_hh_attack[hh] = -1;
                household[hh] = sample_integer(household.size());
            }
            else
            {
                if (omega_vec[hh] == 0)
                {
                    human_seir_status[hh] = HumanSeirStatus::RCVD;
                }
                omega_vec[hh] = omega_vec[hh] - 1;
            }
        }
        else if (human_seir_status[hh] == HumanSeirStatus::RCVD)
        {
            if (get_rand_uniform() < 1 / (74.8 * 365))
            {
                // human_seir_status[hh] = HumanSeirStatus::DISD;
                human_seir_status[hh] = HumanSeirStatus::SUCP;
                time_to_hh_attack[hh] = -1;
                household[hh] = sample_integer(household.size());
            }
            else
            {
                // int hh_region = human_current_region[hh];
                // double lambda_hh_t = kBeta * region_prev_t[hh_region];
                double lambda_hh_t = kBeta * my_prev;

                if (was_latent_status[hh] == 0 && get_rand_uniform() < kAlpha)
                {
                    human_seir_status[hh] = HumanSeirStatus::INFT;
                    tau_I[hh] = std::max(kTau_I_min, static_cast<int>(std::round(get_rand_normal(kTau_I_mu, kTau_I_sd))));
                    infectious_households[household[hh]]++;
                    inc_t++;
                }

                // insert
                else if (get_rand_uniform() < (kNu * lambda_hh_t))
                {
                    human_seir_status[hh] = HumanSeirStatus::EXPD;
                    phi_destination[hh] = get_rand_uniform() > kPhi ? 1 : 0;
                }
                else if (time_to_hh_attack[hh] >= 0)
                {
                    if (time_to_hh_attack[hh] == 0)
                    {
                        human_seir_status[hh] = HumanSeirStatus::EXPD;
                        phi_destination[hh] = get_rand_uniform() > kPhi ? 1 : 0;
                        time_to_hh_attack[hh] = -1; // empty hh attack timer
                        if (phi_destination[hh] == 0)
                        {
                            hh_inc_t++;
                        }
                    }
                    else
                    {
                        time_to_hh_attack[hh]--;
                    }
                }
                // end insert
            }
        }
    }

    // int hh_region = human_current_region[1];
    // double lambda_hh_t = kBeta
    //             * static_cast<double>(region_i_t[hh_region])
    //             / static_cast<double>(region_populations_current[hh_region]);
    // double lambda_hh_t = kBeta * region_prev_t[hh_region];
    double lambda_hh_t = kBeta * my_prev;
    // std::cout << "Infected = " << get_i_t() << " population: " << pop << std::endl;
    // std::cout << "prev: " << my_prev << " lambda: " << lambda_hh_t << std::endl;

    tt++;
    s_t = std::count(human_seir_status.begin(), human_seir_status.end(), HumanSeirStatus::SUCP);
    e_t = std::count(human_seir_status.begin(), human_seir_status.end(), HumanSeirStatus::EXPD);
    i_t = std::count(human_seir_status.begin(), human_seir_status.end(), HumanSeirStatus::INFT);
    l_t = std::count(human_seir_status.begin(), human_seir_status.end(), HumanSeirStatus::LATN);
    a_t = std::count(human_seir_status.begin(), human_seir_status.end(), HumanSeirStatus::TRTD);
    r_t = std::count(human_seir_status.begin(), human_seir_status.end(), HumanSeirStatus::RCVD);
    d_t = std::count(human_seir_status.begin(), human_seir_status.end(), HumanSeirStatus::DISD);
}