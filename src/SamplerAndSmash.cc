#include "SamplerAndSmash.h"

#include "smash/decaymodes.h"
#include "smash/inputfunctions.h"
#include "smash/particles.h"
#include "smash/setup_particles_decaymodes.h"

#include "microcanonical_sampler/main.h"

#include <string>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

namespace bf = boost::filesystem;

namespace {
bf::path default_output_path() {
  const bf::path p = bf::absolute("smash_output");
  if (!bf::exists(p)) {
    return p / "0";
  }
  bf::path p2;
  for (int id = 0; id < std::numeric_limits<int>::max(); ++id) {
    p2 = p / std::to_string(id);
    if (!bf::exists(p2)) {
      break;
    }
  }
  if (p == p2) {
    throw std::runtime_error("no unique data subdir ID left");
  }
  return p2;
}

void ensure_path_is_valid(const bf::path &path) {
  if (bf::exists(path)) {
    if (!bf::is_directory(path)) {
      throw std::runtime_error("The given path (" + path.native() +
                               ") exists, but it is not a directory.");
    }
  } else {
    if (!bf::create_directories(path)) {
      throw std::runtime_error(
          "Race condition detected: The directory " + path.native() +
          " did not exist a few cycles ago, but was created in the meantime by "
          "a different process.");
    }
  }
}

} // namespace

void SamplerAndSmash::InitSmash() {
  bf::path input_config_path(smash_config_filename_);
  if (!bf::exists(input_config_path)) {
    std::cout << "SMASH config file " << smash_config_filename_
	      << " not found.";
    std::exit(-1);
  } else {
    std::cout << "Obtaining SMASH configuration from "
	      << smash_config_filename_ << std::endl;
  }
  smash::Configuration config(input_config_path.parent_path(),
                              input_config_path.filename());
  const auto &log = smash::logger<smash::LogArea::Main>();
  // SMASH logging
  smash::set_default_loglevel(
      config.take({"Logging", "default"}, einhard::ALL));
  smash::create_all_loggers(config["Logging"]);

  auto particles_and_decays =
      smash::load_particles_and_decaymodes(smash_particlelist_filename_,
		                           smash_decaymodes_filename_);
  smash::ParticleType::create_type_list(particles_and_decays.first);
  smash::DecayModes::load_decaymodes(particles_and_decays.second);
  smash::ParticleType::check_consistency();

  // Take care of the random seed. This will make SMASH results reproducible.
  int64_t seed = config.read({"General", "Randomseed"});
  if (seed < 0) {
    config["General"]["Randomseed"] = smash::random::generate_63bit_seed();
  }
  // Read in the rest of configuration
  log.info("Seting up SMASH Experiment object");
  bf::path output_path = default_output_path();
  ensure_path_is_valid(output_path);

  smash_experiment_ =
      smash::make_unique<smash::Experiment<AfterburnerModus>>(config, output_path);
  log.info("Finish initializing SMASH");
}
/*
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
  SamplerAndSmash sampler_and_smash;
  sampler_and_smash.InitSmash();
}
