#include "SamplerAndSmash.h"

#include "smash/decaymodes.h"
#include "smash/inputfunctions.h"
#include "smash/particles.h"
#include "smash/setup_particles_decaymodes.h"

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

  /**
   *   Set up configuration
   */
  std::string smash_config_filename = "../config.yaml";
  bf::path input_config_path(smash_config_filename);
  if (!bf::exists(input_config_path)) {
    std::cout << "SMASH config file " << smash_config_filename
	      << " not found.";
    std::exit(-1);
  } else {
    std::cout << "Obtaining SMASH configuration from "
	      << smash_config_filename << std::endl;
  }
  smash::Configuration config = smash::Configuration(
      input_config_path.parent_path(), input_config_path.filename());
  // Set up logging
  smash::set_default_loglevel(
      config.take({"Logging", "default"}, einhard::ALL));
  smash::create_all_loggers(config["Logging"]);
  const auto &log = smash::logger<smash::LogArea::Main>();

  // Take care of the random seed. This will make SMASH results reproducible.
  int64_t seed = config.read({"General", "Randomseed"});
  if (seed < 0) {
    config["General"]["Randomseed"] = smash::random::generate_63bit_seed();
  }
  std::string smash_particlelist_filename = config.take({"Particles"});
  std::string smash_decaymodes_filename = config.take({"DecayModes"});
  
  /**
   *   Initialize sampler
   */
  log.info("Initializing microcanonical sampler");
  sampler::ParticleListFormat plist_format = sampler::ParticleListFormat::SMASH;
  sampler::read_particle_list(smash_particlelist_filename, plist_format);

  smash::Configuration subconf = config["MicrocanonicalSampler"];
  // Maximal hadron mass to be sampled
  const double max_mass = subconf.take({"MaxMass"}, 2.5);
  const double E_patch = subconf.take({"PatchEnergy"}, 10.0);
  const size_t N_warmup = subconf.take({"WarmupSteps"}, 1E6);
  N_decorrelate_ = subconf.take({"DecorrelationSteps"}, 2E2);
  N_samples_per_hydro_ = subconf.take({"SamplesPerHydro"}, 100);
 
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
  microcanonical_sampler_ = smash::make_unique<MicrocanonicalSampler>(is_sampled_type, 0, quantum_statistics);

  microcanonical_sampler_patches_ = smash::make_unique<std::vector<HyperSurfacePatch>>(
      hyper.split(E_patch));
  size_t number_of_patches = microcanonical_sampler_patches_->size();

  microcanonical_sampler_particles_ = smash::make_unique<std::vector<MicrocanonicalSampler::SamplerParticleList>>();
  microcanonical_sampler_particles_->resize(number_of_patches);

  for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
    std::cout << "Initializing patch " << i_patch << std::endl;
    microcanonical_sampler_->initialize((*microcanonical_sampler_patches_)[i_patch], (*microcanonical_sampler_particles_)[i_patch]);
  }

  std::cout << "Warming up." << std::endl;
  step_until_sufficient_decorrelation(*microcanonical_sampler_, *microcanonical_sampler_patches_, *microcanonical_sampler_particles_, N_warmup);
  std::cout << "Finished warming up." << std::endl;
  size_t total_particles = 0;
  for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
    total_particles += (*microcanonical_sampler_particles_)[i_patch].size();
  }
  std::cout << total_particles << " particles" << std::endl;
  microcanonical_sampler_->print_rejection_stats();

  /**
   *   Initialize SMASH Experiment
   */

  // Initialize SMASH particle types, consistency with the sampler
  // requires the particles file is the same for sampler and SMASH
  auto particles_and_decays =
      smash::load_particles_and_decaymodes(smash_particlelist_filename.c_str(),
		                           smash_decaymodes_filename.c_str());
  smash::ParticleType::create_type_list(particles_and_decays.first);
  smash::DecayModes::load_decaymodes(particles_and_decays.second);
  smash::ParticleType::check_consistency();

  log.info("Seting up SMASH Experiment object");
  bf::path output_path = default_output_path();
  ensure_path_is_valid(output_path);

  smash_experiment_ =
      smash::make_unique<smash::Experiment<AfterburnerModus>>(config, output_path);
  smash_experiment_->modus()->set_sampler_type(SamplerType::Microcanonical);
  log.info("Finish initializing SMASH");

}

void AfterburnerModus::microcanonical_sampler_hadrons_to_smash_particles(
    const std::vector<MicrocanonicalSampler::SamplerParticleList> &sampler_hadrons,
    smash::Particles &smash_particles) {
  smash_particles.reset();
  const size_t n_patches = sampler_hadrons.size();
  for (size_t i_patch = 0; i_patch < n_patches; i_patch++) {
    for (const MicrocanonicalSampler::SamplerParticle &h : sampler_hadrons[i_patch]) {
      const smash::FourVector p = h.momentum;
      const double mass = p.abs();
      const smash::FourVector r =
          (*microcanonical_sampler_patches_)[i_patch].cells()[h.cell_index].r;
      this->try_create_particle(smash_particles, h.type->pdgcode(),
                              r.x0(), r.x1(), r.x2(), r.x3(),
                              mass, p.x0(), p.x1(), p.x2(), p.x3());
    }
  }
}

void SamplerAndSmash::Execute() {
  const auto &log = smash::logger<smash::LogArea::Experiment>();
  for (size_t j = 0; j < N_samples_per_hydro_; j++) {
    step_until_sufficient_decorrelation(*microcanonical_sampler_, *microcanonical_sampler_patches_, *microcanonical_sampler_particles_,
                                        N_decorrelate_);
    AfterburnerModus *modus = smash_experiment_->modus();
    modus->microcanonical_sampler_hadrons_ = microcanonical_sampler_particles_.get();
    modus->microcanonical_sampler_patches_ = microcanonical_sampler_patches_.get();
    log.info("Event ", j);
    smash_experiment_->initialize_new_event();
    smash_experiment_->run_time_evolution();
    smash_experiment_->do_final_decays();
    smash_experiment_->final_output(j);
  }

  // smash::Particles *smash_particles = smash_experiment_->particles();
}

int main() {
  SamplerAndSmash sampler_and_smash;
  sampler_and_smash.Execute();
}
