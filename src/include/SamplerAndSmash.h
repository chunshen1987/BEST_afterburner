#ifndef SAMPLER_AND_SMASH_H
#define SAMPLER_AND_SMASH_H

#include "smash/configuration.h"
#include "smash/experiment.h"
#include "smash/listmodus.h"

#include "microcanonical_sampler/hydro_cells.h"
#include "microcanonical_sampler/microcanonical_sampler.h"

#include "pratt_sampler/master.h"

enum class SamplerType {
 Microcanonical,
 Pratt,
 Sangwook,
 iSS,
};

/**
 * An Afterburner modus, which acts similarly to SMASH ListModus,
 * external particles of smash::Particles class, not from file.
 * This class is needed to use SMASH as a 3rd party afterburner.
 */
class AfterburnerModus : public smash::ListModus {
 public:
  // Unlike for ListModus there is no need to get any data from the config
  AfterburnerModus(smash::Configuration, const smash::ExperimentParameters &) {
    const auto &log = smash::logger<smash::LogArea::Main>();
    log.info("Constructing AfterburnerModus");
  }
  void set_sampler_type(SamplerType sampler_type) { sampler_type_ = sampler_type; }

  void microcanonical_sampler_hadrons_to_smash_particles(
      const std::vector<MicrocanonicalSampler::SamplerParticleList> &sampler_hadrons,
      smash::Particles &smash_particles);

  // This function overrides the function from ListModus.
  double initial_conditions(smash::Particles *particles,
                            const smash::ExperimentParameters &) {
    if (sampler_type_ == SamplerType::Microcanonical) {
      microcanonical_sampler_hadrons_to_smash_particles(
          *microcanonical_sampler_hadrons_, *particles);
    }
    backpropagate_to_same_time(*particles);
    return start_time_;
  }
  std::vector<MicrocanonicalSampler::SamplerParticleList> *microcanonical_sampler_hadrons_;
  std::vector<HyperSurfacePatch> *microcanonical_sampler_patches_;
 private:
  SamplerType sampler_type_;
};

class SamplerAndSmash {
 public:
  SamplerAndSmash();
  void Execute();
 private:
  SamplerType sampler_type_;

  // Microcanonical sampler
  std::unique_ptr<std::vector<HyperSurfacePatch>> microcanonical_sampler_patches_;
  std::unique_ptr<std::vector<MicrocanonicalSampler::SamplerParticleList>>
      microcanonical_sampler_particles_;
  std::unique_ptr<MicrocanonicalSampler> microcanonical_sampler_;
  size_t N_decorrelate_;
  size_t N_samples_per_hydro_;

  // Pratt sampler
  CparameterMap pratt_sampler_parameters_;
  std::unique_ptr<CmeanField_Simple> pratt_sampler_meanfield_;
  std::unique_ptr<CpartList> pratt_sampler_particlelist_;
  std::unique_ptr<CmasterSampler> pratt_sampler_;

  std::unique_ptr<smash::Experiment<AfterburnerModus>> smash_experiment_;
};

#endif // SAMPLER_AND_SMASH_H
