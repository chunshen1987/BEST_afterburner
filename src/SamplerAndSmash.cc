//#include "SamplerAndSmash.h"

#include "smash/decaymodes.h"
#include "smash/inputfunctions.h"
#include "smash/particles.h"
#include "smash/setup_particles_decaymodes.h"

#include <string>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

/*
void SamplerAndSmash::InitSmash() {
  boost::filesystem::path input_config_path(smash_config_filename_);
  if (!boost::filesystem::exists(input_config_path)) {
    std::cout << "SMASH config file " << smash_config << " not found.";
    std::exit(-1);
  } else {
    std::cout << "Obtaining SMASH configuration from " << smash_config_filename << std::endl;
  }
  smash::Configuration config(input_config_path.parent_path(),
                              input_config_path.filename());
  const auto &log = logger<LogArea::Main>();
  // SMASH logging
  smash::set_default_loglevel(
      config.take({"Logging", "default"}, einhard::ALL));
  smash::create_all_loggers(configuration["Logging"]);

  auto particles_and_decays =
      load_particles_and_decaymodes(smash_particlelist_filename_.c_str(), smash_decaymodes_filename_.c_str());
  // Take care of the random seed. This will make SMASH results reproducible.
  int64_t seed = config.read({"General", "Randomseed"});
  if (seed < 0) {
    configuration["General"]["Randomseed"] = random::generate_63bit_seed();
  }
  // Read in the rest of configuration
  log.info("Seting up SMASH Experiment object");
  boost::filesystem::path output_path("./smash_output");
  smash_experiment_ =
      make_unique<smash::Experiment<AfterburnerModus>>(config, output_path);
  log.info("Finish initializing SMASH");
}

void SamplerAndSmash::Execute() {
  modus->jetscape_hadrons_ = soft_particlization_sampler_->Hadron_list_;
  const int n_events = modus->jetscape_hadrons_.size();
  JSINFO << "SMASH: obtained " << n_events << " events from particlization";
  smash::Particles *smash_particles = smash_experiment_->particles();
  for (unsigned int i = 0; i < n_events; i++) {
    JSINFO << "Event " << i << " SMASH starts with "
           << modus->jetscape_hadrons_[i].size() << " particles.";
    smash_experiment_->initialize_new_event();
    smash_experiment_->run_time_evolution();
    smash_experiment_->do_final_decays();
    smash_experiment_->final_output(i);
    smash_particles_to_JS_hadrons(*smash_particles,
                                  modus->jetscape_hadrons_[i]);
    JSINFO << modus->jetscape_hadrons_[i].size() << " hadrons from SMASH.";
  }
}*/

int main() {
  std::cout << "SMASH includes are working" << std::endl;
}
