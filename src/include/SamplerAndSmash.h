#ifndef SAMPLER_AND_SMASH_H
#define SAMPLER_AND_SMASH_H

#include "smash/configuration.h"
#include "smash/experiment.h"
#include "smash/listmodus.h"

/**
 * An Afterburner modus, which acts similarly to SMASH ListModus,
 * external particles of smash::Particles class, not from file.
 * This class is needed to use SMASH as a 3rd party afterburner.
 */
class AfterburnerModus : public smash::ListModus {
 public:
  // Unlike for ListModus there is no need to get any data from the config
  AfterburnerModus(smash::Configuration, const smash::ExperimentParameters &) {
    std::cout << "Constructing AfterburnerModus" << std::endl;
  }

  // The converter is not static, because modus holds int variables
  // for the number of warnings, which are used in try_create_particle,
  // called by this function. Maybe I (oliiny) will change this design in SMASH
  // later, but now I have to put this converter inside the AfterburnerModus.
  void sampler_hadrons_to_smash_particles(
      const std::vector<SamplerHadron> &sampler_hadrons,
      smash::Particles &smash_particles);

  // This function overrides the function from ListModus.
  double initial_conditions(smash::Particles *particles,
                            const smash::ExperimentParameters &) {
    sampler_hadrons_to_smash_particles(sampler_hadrons_, *particles);
    backpropagate_to_same_time(*particles);
    return start_time_;
  }
  std::vector<SamplerHadron> input_hadrons_;
};

class SamplerAndSmash {
private:

  std::string smash_config_filename_ = "./config.yaml",
              // empty string means default for SMASH
              smash_particlelist_filename_ = "",
              smash_decaymodes_filename_ = "";
  std::string sampler_config_filename_ = "",
              sampler_hypersurface_filename_ = "";
  std::unique_ptr<smash::Experiment<AfterburnerModus>> smash_experiment_;
public:
};

#endif // SAMPLER_AND_SMASH_H
