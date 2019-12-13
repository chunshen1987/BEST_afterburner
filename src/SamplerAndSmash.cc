#include "SamplerAndSmash.h"

#include "smash/decaymodes.h"
#include "smash/inputfunctions.h"
#include "smash/particles.h"
#include "smash/setup_particles_decaymodes.h"

#include "microcanonical_sampler/hydro_cells.h"
#include "microcanonical_sampler/microcanonical_sampler.h"
#include "microcanonical_sampler/sampler_particletype_list.h"

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

SamplerAndSmash::SamplerAndSmash() {
  // Set up configuration
  bf::path input_config_path(smash_config_filename_);
  if (!bf::exists(input_config_path)) {
    std::cout << "SMASH config file " << smash_config_filename_
	      << " not found.";
    std::exit(-1);
  } else {
    std::cout << "Obtaining SMASH configuration from "
	      << smash_config_filename_ << std::endl;
  }
  smash::Configuration config = smash::Configuration(
      input_config_path.parent_path(), input_config_path.filename());
  // Set up logging
  smash::set_default_loglevel(
      config.take({"Logging", "default"}, einhard::ALL));
  smash::create_all_loggers(config["Logging"]);

  // Initialize SMASH particle types, which may be also used by the sampler
  std::string smash_particlelist_filename = config.take({"Particles"});
  std::string smash_decaymodes_filename = config.take({"DecayModes"});
  auto particles_and_decays =
      smash::load_particles_and_decaymodes(smash_particlelist_filename.c_str(),
		                           smash_decaymodes_filename.c_str());
  smash::ParticleType::create_type_list(particles_and_decays.first);
  smash::DecayModes::load_decaymodes(particles_and_decays.second);
  smash::ParticleType::check_consistency();

  // Take care of the random seed. This will make SMASH results reproducible.
  int64_t seed = config.read({"General", "Randomseed"});
  if (seed < 0) {
    config["General"]["Randomseed"] = smash::random::generate_63bit_seed();
  }
 
  const auto &log = smash::logger<smash::LogArea::Main>();

  log.info("Seting up SMASH Experiment object");
  bf::path output_path = default_output_path();
  ensure_path_is_valid(output_path);

  smash_experiment_ =
      smash::make_unique<smash::Experiment<AfterburnerModus>>(config, output_path);
  log.info("Finish initializing SMASH");

  /**
   *   Initialize sampler
   */
  log.info("Initializing microcanonical sampler");
  sampler::ParticleListFormat plist_format = sampler::ParticleListFormat::SMASH;
  sampler::read_particle_list(smash_particlelist_filename, plist_format);

  // Maximal hadron mass to be sampled
  smash::Configuration subconf = config["MicrocanonicalSampler"];
  const double max_mass = subconf.take({"MaxMass"}, 2.5);
  const double E_patch = subconf.take({"PatchEnergy"}, 10.0);
  const int N_warmup = subconf.take({"NumberOfWarmupSteps"}, 1E6);

  // Quantum statistics is not implemented properly in the microcanonical sampler
  constexpr bool quantum_statistics = false;

  /**
   * A function, which defines, which species will be sampled. For
   * a given species, if it returns true, the species will be sampled.
   */
  auto is_sampled_type = [&](ParticleTypePtr t) {
    return t->is_hadron() && t->mass() < max_mass;
  };

  std::string hypersurface_input_file = subconf.take({"HyperSurface"}); 
  HyperSurfacePatch::InputFormat hypersurface_format =
    HyperSurfacePatch::InputFormat::MUSIC_ASCII_3plus1D;

  // So far this is dummy, because we do not use 2+1D hydro hypersurface
  std::array<double, 3> eta_for_2Dhydro = {-2.0, 2.0, 0.4};

  HyperSurfacePatch hyper(hypersurface_input_file, hypersurface_format,
                        eta_for_2Dhydro, is_sampled_type, quantum_statistics);
  log.info("Full hypersurface: ", hyper);
  MicrocanonicalSampler sampler(is_sampled_type, 0, quantum_statistics);

  auto patches = hyper.split(E_patch);
  size_t number_of_patches = patches.size();

  std::vector<MicrocanonicalSampler::SamplerParticleList> particles;
  particles.resize(number_of_patches);

  for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
    std::cout << "Initializing patch " << i_patch << std::endl;
    sampler.initialize(patches[i_patch], particles[i_patch]);
  }

  std::cout << "Warming up." << std::endl;
  step_until_sufficient_decorrelation(sampler, patches, particles, N_warmup);
  std::cout << "Finished warming up." << std::endl;
  size_t total_particles = 0;
  for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
    total_particles += particles[i_patch].size();
  }
  std::cout << total_particles << " particles" << std::endl;
  sampler.print_rejection_stats();
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

  for (size_t j = 0; j < N_printout; j++) {
    step_until_sufficient_decorrelation(sampler, patches, particles,
                                        N_decorrelate);
    // print out
    if (j % 10000 == 0 && j != 0) {
      std::cout << "sample " << j << std::endl;
    }
  }

  for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
    MicrocanonicalSampler::QuantumNumbers cons(particles[i_patch]);
    assert((cons.momentum - patches[i_patch].pmu()).abs() < 1.e-4);
    assert(cons.B == patches[i_patch].B());
    assert(cons.S == patches[i_patch].S());
    assert(cons.Q == patches[i_patch].Q());
  }
  sampler.print_rejection_stats();

}*/

int main() {
  std::cout << "SMASH includes are working" << std::endl;
  SamplerAndSmash sampler_and_smash;
}
