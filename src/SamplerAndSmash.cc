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
  sampler_type_ = SamplerType::Pratt; //Microcanonical;  // Todo: get it from config
  if (sampler_type_ == SamplerType::Microcanonical) {
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
    microcanonical_sampler_ = smash::make_unique<MicrocanonicalSampler>(
      is_sampled_type, 0, quantum_statistics);

    microcanonical_sampler_patches_ = smash::make_unique<std::vector<HyperSurfacePatch>>(
        hyper.split(E_patch));
    size_t number_of_patches = microcanonical_sampler_patches_->size();

    microcanonical_sampler_particles_ =
      smash::make_unique<std::vector<MicrocanonicalSampler::SamplerParticleList>>();
    microcanonical_sampler_particles_->resize(number_of_patches);

    for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
      std::cout << "Initializing patch " << i_patch << std::endl;
      microcanonical_sampler_->initialize((*microcanonical_sampler_patches_)[i_patch],
                                          (*microcanonical_sampler_particles_)[i_patch]);
    }

    std::cout << "Warming up." << std::endl;
    step_until_sufficient_decorrelation(*microcanonical_sampler_,
                                        *microcanonical_sampler_patches_,
                                        *microcanonical_sampler_particles_, N_warmup);
    std::cout << "Finished warming up." << std::endl;
    size_t total_particles = 0;
    for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
      total_particles += (*microcanonical_sampler_particles_)[i_patch].size();
    }
    std::cout << total_particles << " particles" << std::endl;
    microcanonical_sampler_->print_rejection_stats();
  }

  if (sampler_type_ == SamplerType::Pratt) {
    log.info("Initializing Pratt sampler");
    pratt_sampler_parameters_.ReadParsFromFile("parameters.dat");
    log.info("Creating particle list");
    pratt_sampler_particlelist_ = smash::make_unique<CpartList>(&pratt_sampler_parameters_);
    log.info("Creating mean field");
    pratt_sampler_meanfield_ = smash::make_unique<CmeanField_Simple>(&pratt_sampler_parameters_);
    log.info("Creating sampler");

    pratt_sampler_ = smash::make_unique<CmasterSampler>(&pratt_sampler_parameters_);
    pratt_sampler_->meanfield = pratt_sampler_meanfield_.get();
    pratt_sampler_->partlist = pratt_sampler_particlelist_.get();
    pratt_sampler_->randy->reset(time(NULL));
    pratt_sampler_->ReadHyper2D();
  }

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
  smash_experiment_->modus()->set_sampler_type(sampler_type_);
  log.info("Finish initializing SMASH");

}

void AfterburnerModus::sampler_hadrons_to_smash_particles(
    smash::Particles &smash_particles) {
  smash_particles.reset();
  const auto &log = smash::logger<smash::LogArea::Experiment>();

  if (sampler_type_ == SamplerType::Microcanonical) {
    const size_t n_patches = microcanonical_sampler_hadrons_->size();
    for (size_t i_patch = 0; i_patch < n_patches; i_patch++) {
      for (const MicrocanonicalSampler::SamplerParticle &h :
               (*microcanonical_sampler_hadrons_)[i_patch]) {
        const smash::FourVector p = h.momentum;
        const double mass = p.abs();
        const smash::FourVector r =
            (*microcanonical_sampler_patches_)[i_patch].cells()[h.cell_index].r;
        this->try_create_particle(smash_particles, h.type->pdgcode(),
                                r.x0(), r.x1(), r.x2(), r.x3(),
                                mass, p.x0(), p.x1(), p.x2(), p.x3());
      }
    }
  } else if (sampler_type_ == SamplerType::Pratt) {
    log.info("Transfering ", pratt_sampler_hadrons_->size(),
             " particles from sampler to smash.");
    for (const Cpart &h : *pratt_sampler_hadrons_) {
      const smash::FourVector p(h.p[0], h.p[1], h.p[2], h.p[3]),
                              r(h.r[0], h.r[1], h.r[2], h.r[3]);
      const double mass = p.abs();
      log.info("pdg = ", h.pid, ", p = ", p, ", r = ", r);
      this->try_create_particle(smash_particles,
                                smash::PdgCode::from_decimal(h.pid),
                                h.r[0], h.r[1], h.r[2], h.r[3],
                                mass, h.p[0], h.p[1], h.p[2], h.p[3]);
    }
  }
}

void SamplerAndSmash::Execute() {
  const auto &log = smash::logger<smash::LogArea::Experiment>();
  for (size_t j = 0; j < N_samples_per_hydro_; j++) {
    AfterburnerModus *modus = smash_experiment_->modus();

    if (sampler_type_ == SamplerType::Microcanonical) {
      step_until_sufficient_decorrelation(
        *microcanonical_sampler_, *microcanonical_sampler_patches_,
        *microcanonical_sampler_particles_, N_decorrelate_);
      modus->microcanonical_sampler_hadrons_ = microcanonical_sampler_particles_.get();
      modus->microcanonical_sampler_patches_ = microcanonical_sampler_patches_.get();
    }

    if (sampler_type_ == SamplerType::Pratt) {
      pratt_sampler_->MakeEvent();
      modus->pratt_sampler_hadrons_ = &(pratt_sampler_->partlist->partvec);
      // Fix a bug/feature of the wrong vector size in Pratt sampler
      modus->pratt_sampler_hadrons_->resize(pratt_sampler_->partlist->nparts);
    }

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
