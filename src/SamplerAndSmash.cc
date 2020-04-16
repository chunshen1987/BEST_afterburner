#include <getopt.h>

#include "SamplerAndSmash.h"

#include "smash/decaymodes.h"
#include "smash/inputfunctions.h"
#include "smash/particles.h"
#include "smash/setup_particles_decaymodes.h"
#include "smash/stringfunctions.h"

#include "microcanonical_sampler/sampler_particletype_list.h"

#include <string>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

namespace bf = boost::filesystem;

namespace {
bf::path default_output_path(std::string output_folder) {
    const bf::path p = bf::absolute(output_folder);
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
            throw std::runtime_error("Race condition detected: The directory " +
                                     path.native() +
                                     " did not exist a few cycles ago, but was "
                                     "created in the meantime by "
                                     "a different process.");
        }
    }
}

bool avoid_in_sampler_table(const smash::ParticleTypePtr ptype) {
    const smash::PdgCode pdg = ptype->pdgcode();
    std::array<int, 3> q = pdg.quark_content();
    bool has_heavy_quark =
        std::abs(q[0]) > 3 || std::abs(q[1]) > 3 || std::abs(q[2]) > 3;
    return pdg.is_lepton() ||
           (pdg.is_nucleus() && !ptype->is_stable()) ||  // Unstable nuclei, like t'
           (pdg.is_hadron() && has_heavy_quark);         // heavy flavor hadrons
}

bool avoid_sampling(const smash::ParticleTypePtr ptype) {
    // Photons should be in the table for decaymodes consistency,
    // but they should not be sampled
    return !ptype->is_hadron() || avoid_in_sampler_table(ptype);
}

void smash_particles_to_iSS_format(
    const std::string &sampler_particles_filename,
    const std::string &sampler_particles_selector_filename) {
    FILE *sampler_particles_file;
    FILE *sampler_particles_selector_file;
    sampler_particles_file = fopen(sampler_particles_filename.c_str(), "w");
    sampler_particles_selector_file =
        fopen(sampler_particles_selector_filename.c_str(), "w");

    smash::ParticleTypePtrList list;
    list.clear();
    for (const auto &ptype : smash::ParticleType::list_all()) {
        if (!ptype.pdgcode().is_lepton() && ptype.pdgcode().charmness() == 0 &&
            ptype.pdgcode().bottomness() == 0) {
            list.push_back(&ptype);
        }
    }
    std::sort(list.begin(), list.end(),
              [](smash::ParticleTypePtr a, smash::ParticleTypePtr b) {
                  return a->mass() < b->mass();
              });

    for (const smash::ParticleTypePtr ptype : list) {
        if (!avoid_sampling(ptype)) {
            fprintf(sampler_particles_selector_file, "%13i\n",
                    ptype->pdgcode().get_decimal());
        }
        if (ptype->baryon_number() < 0 || avoid_in_sampler_table(ptype)) {
            continue;
        }
        const auto &decay_modes = ptype->decay_modes();
        const auto &modelist = decay_modes.decay_mode_list();
        int ndecays = ptype->is_stable() ? 1 : modelist.size();
        fprintf(sampler_particles_file,
                "%13i %s %10.5f %10.5f %5i %5i %5i %5i %5i %5i %5i %5i\n",
                ptype->pdgcode().get_decimal(),
                smash::utf8::fill_left(ptype->name(), 12, ' ').c_str(),
                ptype->mass(), ptype->width_at_pole(),
                ptype->pdgcode().spin_degeneracy(), ptype->baryon_number(),
                ptype->strangeness(), ptype->pdgcode().charmness(),
                ptype->pdgcode().bottomness(), ptype->isospin() + 1, ptype->charge(),
                ndecays);
        if (!ptype->is_stable()) {
            for (const auto &decay : modelist) {
                auto ptypes = decay->particle_types();
                fprintf(sampler_particles_file,
                        "%13i %13i %20.5f %13i %13i %13i %13i %13i %13i\n",
                        ptype->pdgcode().get_decimal(), 2, decay->weight(),
                        ptypes[0]->pdgcode().get_decimal(),
                        ptypes[1]->pdgcode().get_decimal(), 0, 0, 0,
                        decay->angular_momentum());
            }
        } else {
            fprintf(sampler_particles_file,
                    "%13i %13i %20.5f %13i %13i %13i %13i %13i %13i\n",
                    ptype->pdgcode().get_decimal(), 1, 1.0,
                    ptype->pdgcode().get_decimal(), 0, 0, 0, 0, 0);
        }
    }
    fclose(sampler_particles_file);
    fclose(sampler_particles_selector_file);
}

}  // namespace

SamplerAndSmash::SamplerAndSmash(std::string config_filename,
                                 const std::vector<std::string> &extra_config) {
    /**
     *   Set up configuration
     */
    bf::path input_config_path(config_filename);
    if (!bf::exists(input_config_path)) {
        std::cout << "Config file " << config_filename << " not found." << std::endl;
        std::exit(-1);
    } else {
        std::cout << "Obtaining configuration from " << config_filename << std::endl;
    }
    smash::Configuration config = smash::Configuration(
        input_config_path.parent_path(), input_config_path.filename());
    for (const auto &c : extra_config) {
        config.merge_yaml(c);
    }

    // Set up logging
    smash::set_default_loglevel(config.take({"Logging", "default"}, einhard::ALL));
    smash::create_all_loggers(config["Logging"]);
    const auto &log = smash::logger<smash::LogArea::Main>();

    // Take care of the random seed. This will make SMASH results reproducible.
    int64_t seed = config.read({"General", "Randomseed"});
    if (seed < 0) {
        config["General"]["Randomseed"] = smash::random::generate_63bit_seed();
    }

    // Prepare output directory
    std::string output_dir =
        config.take({"Output_Directory"}, std::string("smash_output"));
    bf::path output_path = default_output_path(output_dir);
    ensure_path_is_valid(output_path);
    // Keep a copy of the used configuration in the output directory
    std::ofstream config_copy;
    config_copy.open((output_path / "config.yaml").c_str());
    config_copy << config.to_string();
    config_copy.close();

    // Initialize sampler type and its general parameters
    std::string sampler_type_str = config.take({"Sampler", "Type"});
    if (sampler_type_str == "Microcanonical") {
        sampler_type_ = SamplerType::Microcanonical;
    } else if (sampler_type_str == "MSU") {
        sampler_type_ = SamplerType::MSU;
    } else if (sampler_type_str == "iSS") {
        sampler_type_ = SamplerType::iSS;
    } else {
        log.error("Unknown sampler type: ", sampler_type_str);
        throw std::runtime_error("Unknown sampler type.");
    }

    N_samples_per_hydro_ = config.read({"General", "Nevents"});
    const bool shear_deltaf = config.take({"Sampler", "ViscousCorections", "Shear"});
    const bool bulk_deltaf = config.take({"Sampler", "ViscousCorections", "Bulk"});
    const bool diff_deltaf =
        config.take({"Sampler", "ViscousCorections", "Diffusion"});
    const bool quantum_stat = config.take({"Sampler", "QuantumCorrections"});
    const bool spectral_functions = config.take({"Sampler", "SpectralFunctions"});

    /** Initialize SMASH particle types.
     * A subset of these particles will be used for sampling.
     * This is necessary for consistency: SMASH should be able to handle
     * the particles sampler provides. It is also good for consistent comparison
     * between different samplers: they are given the same particle table.
     */
    std::string smash_particlelist_filename = config.take({"Particles"});
    std::string smash_decaymodes_filename = config.take({"DecayModes"});
    log.info("Initializing particles table from ", smash_particlelist_filename,
             " and ", smash_decaymodes_filename);
    auto particles_and_decays = smash::load_particles_and_decaymodes(
        smash_particlelist_filename.c_str(), smash_decaymodes_filename.c_str());
    smash::ParticleType::create_type_list(particles_and_decays.first);
    smash::DecayModes::load_decaymodes(particles_and_decays.second);
    smash::ParticleType::check_consistency();

    // Temporary file for particle table in iSS format. It is filled automatically
    // below and gaurantees consistency between different samplers, because
    // ALL samplers are going to be initialized from it.
    const std::string sampler_particles_filename(
        (output_path / "pdg-SMASH.dat").c_str());
    // Additional temporary file, specific for iSS sampler.
    // For our purposes this file should always contain all the particles
    // from sampler_particles.
    const std::string sampler_particles_selector_filename(
        (output_path / "chosen_particles_SMASH.dat").c_str());

    std::string external_codes_dir(
        (bf::path(__FILE__).parent_path() / "/../external_codes").string());
    log.info("External codes directory: ", external_codes_dir);

    // Convert SMASH particle table to iSS format and cut off particles
    // that definitely should not be sampled (such as photons, leptons,
    // open charm and beauty hadrons, etc).
    smash_particles_to_iSS_format(sampler_particles_filename,
                                  sampler_particles_selector_filename);

    const std::string hypersurface_input_file = config.take({"HyperSurface"});

    if (sampler_type_ == SamplerType::Microcanonical) {
        log.info("Initializing microcanonical sampler");
        sampler::ParticleListFormat plist_format = sampler::ParticleListFormat::iSS;
        sampler::read_particle_list(sampler_particles_filename, plist_format);

        // Maximal hadron mass to be sampled
        const double E_patch =
            config.take({"Sampler", "Microcanonical", "PatchEnergy"}, 20.0);
        const size_t N_warmup = 1E6;
        N_decorrelate_ = 2E2;

        if (shear_deltaf || bulk_deltaf || diff_deltaf) {
            log.error("Viscous corrections are not implemented",
                      " in microcanonical sampler");
        }
        if (spectral_functions) {
            log.error("Spectral functions are not implemented",
                      "in microcanonical sampler. Pole masses will be used.");
        }

        // Quantum statistics is not implemented properly in the microcanonical
        // sampler
        if (quantum_stat) {
            log.error("Quantum statistics is not implemented in the microcanonical",
                      " sampler.");
        }
        /**
         * A function, which defines, which species will be sampled. For
         * a given species, if it returns true, the species will be sampled.
         */
        auto is_sampled_type = [](sampler::ParticleTypePtr t) {
            const smash::PdgCode pdg = t->pdgcode();
            smash::ParticleTypePtr ptype = &smash::ParticleType::find(pdg);
            return !avoid_sampling(ptype);
        };

        HyperSurfacePatch::InputFormat hypersurface_format =
            HyperSurfacePatch::InputFormat::MUSIC_ASCII_3plus1D;

        // So far this is dummy, because we do not use 2+1D hydro hypersurface
        std::array<double, 3> eta_for_2Dhydro = {-2.0, 2.0, 0.4};

        HyperSurfacePatch hyper(hypersurface_input_file, hypersurface_format,
                                eta_for_2Dhydro, is_sampled_type, quantum_stat);
        log.info("Full hypersurface: ", hyper);
        microcanonical_sampler_ =
            smash::make_unique<MicrocanonicalSampler>(is_sampled_type, 0, false);

        microcanonical_sampler_patches_ =
            smash::make_unique<std::vector<HyperSurfacePatch>>(hyper.split(E_patch));
        size_t number_of_patches = microcanonical_sampler_patches_->size();

        microcanonical_sampler_particles_ = smash::make_unique<
            std::vector<MicrocanonicalSampler::SamplerParticleList>>();
        microcanonical_sampler_particles_->resize(number_of_patches);

        for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
            std::cout << "Initializing patch " << i_patch << std::endl;
            microcanonical_sampler_->initialize(
                (*microcanonical_sampler_patches_)[i_patch],
                (*microcanonical_sampler_particles_)[i_patch]);
        }

        std::cout << "Warming up." << std::endl;
        step_until_sufficient_decorrelation(
            *microcanonical_sampler_, *microcanonical_sampler_patches_,
            *microcanonical_sampler_particles_, N_warmup);
        std::cout << "Finished warming up." << std::endl;
        size_t total_particles = 0;
        for (size_t i_patch = 0; i_patch < number_of_patches; i_patch++) {
            total_particles += (*microcanonical_sampler_particles_)[i_patch].size();
        }
        std::cout << total_particles << " particles" << std::endl;
        microcanonical_sampler_->print_rejection_stats();
    }

    if (sampler_type_ == SamplerType::MSU) {
        log.info("Initializing MSU sampler");
        std::map<std::string, std::string> msu_parameters = {
            {"RESONANCES_INFO_FILE", sampler_particles_filename},
            {"RESONANCES_DECAYS_FILE", sampler_particles_filename},
            {"HYPER_INFO_FILE", hypersurface_input_file},
            {"SAMPLER_SFDIRNAME",
             external_codes_dir +
                 "/best_sampler/software/resinfo/spectralfunctions"},
            // Todo(MSU): make sure MSU sampler complains if hypersurface T is out of
            // this range
            {"SAMPLER_TFMIN", "0.110"},
            {"SAMPLER_TFMAX", "0.170"},
            {"SAMPLER_NTF", "60"},
            {"SAMPLER_SIGMAFMIN", "0.093"},
            {"SAMPLER_SIGMAFMAX", "0.093"},
            {"SAMPLER_NSIGMAF", "1"},
            {"SAMPLER_SETMU0", "false"},
            {"SAMPLER_BOSE_CORR", std::to_string(quantum_stat)},
            {"SAMPLER_N_BOSE_CORR", "5"},
            {"SAMPLER_USE_POLE_MASS", std::to_string(!spectral_functions)},
            {"VISCOUSCORRECTION", std::to_string(shear_deltaf)}};
        for (const auto &par : msu_parameters) {
            msu_sampler_parameters_.set(par.first, par.second);
        }

        if (bulk_deltaf || diff_deltaf) {
            log.error("MSU sampler does not include bulk or diffusion corrections.");
        }

        log.info("MSU sampler: creating particle list");
        msu_sampler_particlelist_ =
            (smash::make_unique<msu_sampler::CpartList>(&msu_sampler_parameters_));
        log.info("MSU sampler: creating mean field");
        msu_sampler_meanfield_ = (smash::make_unique<msu_sampler::CmeanField_Simple>(
            &msu_sampler_parameters_));
        log.info("MSU sampler: creating sampler object");

        msu_sampler_ = smash::make_unique<msu_sampler::CmasterSampler>(
            &msu_sampler_parameters_);
        msu_sampler_->meanfield = msu_sampler_meanfield_.get();
        msu_sampler_->partlist = msu_sampler_particlelist_.get();
        msu_sampler_->randy->reset(time(NULL));
        msu_sampler_->ReadHyper();
    }

    if (sampler_type_ == SamplerType::iSS) {
#ifdef iSSFlag
        log.info("Initializing the iSS sampler");
        bf::path iSS_dir(external_codes_dir + "/iSS");

        // Assume that music_input file with hydro parameter, that iSS reads in
        // is in the same directory with hypersurface
        iSpectraSampler_ptr_ = smash::make_unique<iSS>(
            bf::path(hypersurface_input_file).parent_path().string(),
            (iSS_dir / "iSS_tables").string(),
            output_path
                .string(),  // Particle tables were automatically created there
            (iSS_dir / "iSS_parameters.dat").string());

        iSpectraSampler_ptr_->paraRdr_ptr->setVal("number_of_repeated_sampling",
                                                  N_samples_per_hydro_);
        int64_t random_seed = config.read({"General", "Randomseed"});
        iSpectraSampler_ptr_->set_random_seed(random_seed);
        iSpectraSampler_ptr_->paraRdr_ptr->setVal("hydro_mode", 2);
        iSpectraSampler_ptr_->paraRdr_ptr->setVal("quantum_statistics",
                                                  quantum_stat);

        // set the hadronic afterburner to be SMASH
        iSpectraSampler_ptr_->paraRdr_ptr->setVal("afterburner_type", 2);

        iSpectraSampler_ptr_->paraRdr_ptr->setVal("output_samples_into_files", 0);
        iSpectraSampler_ptr_->paraRdr_ptr->setVal("use_OSCAR_format", 0);
        iSpectraSampler_ptr_->paraRdr_ptr->setVal("use_gzip_format", 0);
        iSpectraSampler_ptr_->paraRdr_ptr->setVal("store_samples_in_memory", 1);
        iSpectraSampler_ptr_->paraRdr_ptr->setVal("perform_decays", 0);

        iSpectraSampler_ptr_->paraRdr_ptr->setVal("include_deltaf_shear",
                                                  shear_deltaf);
        iSpectraSampler_ptr_->paraRdr_ptr->setVal("include_deltaf_bulk",
                                                  bulk_deltaf);
        iSpectraSampler_ptr_->paraRdr_ptr->setVal("bulk_deltaf_kind", 1);
        iSpectraSampler_ptr_->paraRdr_ptr->setVal("include_deltaf_diffusion",
                                                  diff_deltaf);
        iSpectraSampler_ptr_->paraRdr_ptr->setVal("restrict_deltaf", 0);
        iSpectraSampler_ptr_->paraRdr_ptr->setVal("deltaf_max_ratio", 1.0);
        iSpectraSampler_ptr_->paraRdr_ptr->setVal("f0_is_not_small", 1);

        iSpectraSampler_ptr_->paraRdr_ptr->setVal("calculate_vn", 0);
        iSpectraSampler_ptr_->paraRdr_ptr->setVal("MC_sampling", 2);
        iSpectraSampler_ptr_->paraRdr_ptr->setVal(
            "sample_upto_desired_particle_number", 0);
        if (spectral_functions) {
            log.error("iSS does not support spectral functions. ",
                      "Pole masses will be used.");
        }
        iSpectraSampler_ptr_->paraRdr_ptr->echo();
#else
        log.info("Please compile the BEST afterburner with the iSS sampler");
        log.info("cmake .. -DiSS=ON");
        exit(1);
#endif
    }

    // Get rid of configuration for unused samplers
    for (const std::string s : {"Microcanonical", "MSU", "iSS"}) {
        if (s != sampler_type_str) {
            if (config.has_value({"Sampler", s.c_str()})) {
                config.take({"Sampler", s.c_str()});
            }
        }
    }

    /**
     *   Initialize SMASH Experiment
     */

    log.info("Seting up SMASH Experiment object");
    smash_experiment_ = (smash::make_unique<smash::Experiment<AfterburnerModus>>(
        config, output_path));
    smash_experiment_->modus()->set_sampler_type(sampler_type_);
    log.info("Finish initializing SMASH");

    // Complain, if there are unused configuration values
    const std::string report = config.unused_values_report();
    if (report != "{}") {
        log.warn("The following configuration values were not used:\n", report);
    }
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
                    ((*microcanonical_sampler_patches_)[i_patch]
                         .cells()[h.cell_index]
                         .r);
                this->try_create_particle(smash_particles, h.type->pdgcode(), r.x0(),
                                          r.x1(), r.x2(), r.x3(), mass, p.x0(),
                                          p.x1(), p.x2(), p.x3());
            }
        }
    } else if (sampler_type_ == SamplerType::MSU) {
        log.info("Transfering ", msu_sampler_hadrons_->size(),
                 " particles from sampler to smash.");
        for (const msu_sampler::Cpart &h : *msu_sampler_hadrons_) {
            const smash::FourVector p(h.p[0], h.p[1], h.p[2], h.p[3]),
                r(h.r[0], h.r[1], h.r[2], h.r[3]);
            const double mass = p.abs();
            log.debug("pdg = ", h.pid, ", p = ", p, ", r = ", r);
            this->try_create_particle(
                smash_particles, smash::PdgCode::from_decimal(h.pid), h.r[0], h.r[1],
                h.r[2], h.r[3], mass, h.p[0], h.p[1], h.p[2], h.p[3]);
        }
    } else if (sampler_type_ == SamplerType::iSS) {
#ifdef iSSFlag
        log.info("Transfering ", iSS_hadrons_->size(),
                 " particles from sampler to smash.");
        for (const auto &ihadron : *iSS_hadrons_) {
            log.debug("pdg = ", ihadron.pid, ", mass = ", ihadron.mass,
                      " GeV, E = ", ihadron.E, " GeV");
            this->try_create_particle(
                smash_particles, smash::PdgCode::from_decimal(ihadron.pid),
                ihadron.t, ihadron.x, ihadron.y, ihadron.z, ihadron.mass, ihadron.E,
                ihadron.px, ihadron.py, ihadron.pz);
        }
#endif
    }
}

void SamplerAndSmash::Execute() {
    const auto &log = smash::logger<smash::LogArea::Experiment>();

    if (sampler_type_ == SamplerType::iSS) {
        // The iSS sample all events at once
#ifdef iSSFlag
        iSpectraSampler_ptr_->shell();
#endif
    }

    for (size_t j = 0; j < N_samples_per_hydro_; j++) {
        // event loop: pass events one by one from the sampler to SMASH
        AfterburnerModus *modus = smash_experiment_->modus();

        if (sampler_type_ == SamplerType::Microcanonical) {
            step_until_sufficient_decorrelation(
                *microcanonical_sampler_, *microcanonical_sampler_patches_,
                *microcanonical_sampler_particles_, N_decorrelate_);
            modus->microcanonical_sampler_hadrons_ =
                microcanonical_sampler_particles_.get();
            modus->microcanonical_sampler_patches_ =
                microcanonical_sampler_patches_.get();
        }

        if (sampler_type_ == SamplerType::MSU) {
            msu_sampler_->MakeEvent();
            modus->msu_sampler_hadrons_ = &(msu_sampler_->partlist->partvec);
            // Fix a bug/feature of the wrong vector size in MSU sampler
            modus->msu_sampler_hadrons_->resize(msu_sampler_->partlist->nparts);
        }

        if (sampler_type_ == SamplerType::iSS) {
#ifdef iSSFlag
            modus->iSS_hadrons_ = iSpectraSampler_ptr_->get_hadron_list_iev(j);
#endif
        }

        log.info("Event ", j);
        smash_experiment_->initialize_new_event();
        smash_experiment_->run_time_evolution();
        smash_experiment_->do_final_decays();
        smash_experiment_->final_output(j);
    }

    if (sampler_type_ == SamplerType::iSS) {
#ifdef iSSFlag
        iSpectraSampler_ptr_->clear();
#endif
    }
}

namespace {
void usage(const int rc, const std::string &progname) {
    std::printf("\nUsage: %s [option]\n\n", progname.c_str());
    std::printf(
        "  -h, --help              usage information\n"
        "  -c, --configfile <file> path to configuration file\n"
        "                          (default: ../config.yaml)\n"
        "  -o, --option <YAML>     override option from config file\n"
        "                          Can be used multiple times\n"
        "                          For example:\n"
        "                          ./%s\n"
        "                             -o \"Output: {Output_Interval: 10.0}\"\n"
        "                             -o \"Sampler: {Type: iSS}\"\n"
        "                             -o \"Output_Directory: "
        "iSS_plus_SMASH_results\"\n\n",
        progname.c_str());
    std::exit(rc);
}
}  // anonymous namespace

int main(int argc, char *argv[]) {
    constexpr option longopts[] = {{"configfile", required_argument, 0, 'c'},
                                   {"option", required_argument, 0, 'o'},
                                   {"help", no_argument, 0, 'h'},
                                   {nullptr, 0, 0, 0}};

    int opt;
    const std::string progname = bf::path(argv[0]).filename().native();
    std::string config_file_path = "../config.yaml";
    std::vector<std::string> extra_config;
    while ((opt = getopt_long(argc, argv, "c:o:h", longopts, nullptr)) != -1) {
        switch (opt) {
            case 'o':
                extra_config.emplace_back(optarg);
                break;
            case 'c':
                config_file_path = optarg;
                break;
            case 'h':
                usage(EXIT_SUCCESS, progname);
                break;
            default:
                usage(EXIT_FAILURE, progname);
        }
    }

    SamplerAndSmash sampler_and_smash(config_file_path, extra_config);
    sampler_and_smash.Execute();
}
