#ifndef SAMPLER_AND_SMASH_H
#define SAMPLER_AND_SMASH_H

#include "smash/configuration.h"
#include "smash/experiment.h"
#include "smash/listmodus.h"

#include "microcanonical_sampler/hydro_cells.h"
#include "microcanonical_sampler/microcanonical_sampler.h"

#include "msu_sampler/master.h"
#include "msu_sampler/part.h"

#ifdef iSSFlag
    #include "iSS.h"
#endif

enum class SamplerType {
    Microcanonical,
    MSU,
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
    AfterburnerModus(smash::Configuration,
                     const smash::ExperimentParameters &) {
        const auto &log = smash::logger<smash::LogArea::Main>();
        log.debug("Constructing AfterburnerModus");
    }

    void set_sampler_type(SamplerType sampler_type) {
        sampler_type_ = sampler_type;
    }

    void sampler_hadrons_to_smash_particles(smash::Particles &smash_particles);

    // This function overrides the function from ListModus.
    double initial_conditions(smash::Particles *particles,
                              const smash::ExperimentParameters &) {
        sampler_hadrons_to_smash_particles(*particles);
        backpropagate_to_same_time(*particles);
        return start_time_;
    }

    // Microcanonical sampler hook-up
    std::vector<MicrocanonicalSampler::SamplerParticleList>
        *microcanonical_sampler_hadrons_;
    std::vector<HyperSurfacePatch> *microcanonical_sampler_patches_;

    // MSU sampler hook-up
    std::vector<msu_sampler::Cpart> *msu_sampler_hadrons_;

    // iSS sampler hook-up
#ifdef iSSFlag
    std::vector<iSS_Hadron> *iSS_hadrons_;
#endif

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
    std::unique_ptr<std::vector<HyperSurfacePatch>>
        microcanonical_sampler_patches_;
    std::unique_ptr<std::vector<MicrocanonicalSampler::SamplerParticleList>>
        microcanonical_sampler_particles_;
    std::unique_ptr<MicrocanonicalSampler> microcanonical_sampler_;
    size_t N_decorrelate_;
    size_t N_samples_per_hydro_;

    // MSU sampler
    msu_sampler::CparameterMap msu_sampler_parameters_;
    std::unique_ptr<msu_sampler::CmeanField_Simple> msu_sampler_meanfield_;
    std::unique_ptr<msu_sampler::CpartList> msu_sampler_particlelist_;
    std::unique_ptr<msu_sampler::CmasterSampler> msu_sampler_;

    // iSS sampler
#ifdef iSSFlag
    std::unique_ptr<iSS> iSpectraSampler_ptr_;
#endif

    // SMASH
    std::unique_ptr<smash::Experiment<AfterburnerModus>> smash_experiment_;
};

#endif // SAMPLER_AND_SMASH_H
