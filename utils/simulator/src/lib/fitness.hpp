#pragma once

#include "random.hpp"
#include "tclap/CmdLine.h"

class FitnessState {
  public:
    double drift{0};
    double fitness_optimum{0};
    explicit FitnessState() = default;
    ~FitnessState() = default;
};

class FitnessModel {
  protected:
    FitnessState state{};

  public:
    explicit FitnessModel() = default;
    virtual ~FitnessModel() = default;

    virtual double fitness(double v) const = 0;

    void set_state(FitnessState const &s) { state = s; }
    FitnessState get_state() const { return state; }
    virtual void update() = 0;
};

class NeutralModel : public FitnessModel {
  public:
    explicit NeutralModel() = default;
    ~NeutralModel() override = default;

    double fitness(double v) const override { return 1.0; }
    void update() override {};
};


class OrnsteinUhlenbeckBiasModel : public FitnessModel {
  protected:
    double peakness{};
    double epistasis{};
    double bias{};
    double bias_rate{};
    double sigma{};
    double theta{};
    std::exponential_distribution<double> exponential_distrib;

  public:
    static int nbr_changes;

    explicit OrnsteinUhlenbeckBiasModel() = default;
    explicit OrnsteinUhlenbeckBiasModel(double const &peakness, double const &epistasis,
        double const &bias_optimum, double const &bias_optimum_rate, double const &sigma_optimum,
        double const &theta_optimum)
        : peakness{peakness},
          epistasis{epistasis},
          bias{bias_optimum},
          bias_rate{bias_optimum_rate},
          sigma{sigma_optimum},
          theta{theta_optimum} {
        // Exponential distribution with rate bias
        if (bias_rate < 0.0 or bias_rate > 1.0) {
            throw std::invalid_argument("'bias_optimum_rate' must be in [0, 1]");
        }
        if (bias < 0.0) { throw std::invalid_argument("'bias_optimum' must be >= 0"); }
        if (sigma < 0.0) { throw std::invalid_argument("'sigma_optimum' must be >= 0"); }
        if (theta < 0.0 or theta >= 1.0) { throw std::invalid_argument("'theta_optimum' must be in [0, 1)"); }
        if (peakness < 0.0) { throw std::invalid_argument("'peakness' must be >= 0"); }
        if (epistasis < 0.0) { throw std::invalid_argument("'epistasis' must be >= 0"); }

        exponential_distrib = std::exponential_distribution<double>(1.0 / bias);
    }

    void update() override {
        if (bias_rate != 0.0) {
            if (uniform_distrib(generator) < bias_rate) {
                nbr_changes++;
                if (bernouilli_distrib(generator)) {
                    state.drift += exponential_distrib(generator);
                } else {
                    state.drift -= exponential_distrib(generator);
                }
            }
        } else {
            state.drift += bias;
        }
        state.fitness_optimum +=
            sigma * normal_distrib(generator) + theta * (state.drift - state.fitness_optimum);
    }

    double fitness(double v) const override {
        return exp(-peakness * pow((v - state.fitness_optimum), epistasis));
    }
    ~OrnsteinUhlenbeckBiasModel() override = default;
};

class OrnsteinUhlenbeckBiasArgParse {
  protected:
    TCLAP::CmdLine &cmd;

  public:
    explicit OrnsteinUhlenbeckBiasArgParse(TCLAP::CmdLine &cmd) : cmd{cmd} {}

    TCLAP::ValueArg<double> bias_optimum{"", "bias_optimum",
        "The Ornstein–Uhlenbeck bias (>=0) applied to the fitness optimum at each generation",
        false, 0.0, "double", cmd};
    TCLAP::ValueArg<double> bias_optimum_rate{"", "bias_optimum_rate",
        "The rate (probability, 0<=r<1) at which we change the bias_optimum", false, 0.0, "double",
        cmd};
    TCLAP::ValueArg<double> sigma_optimum{"", "sigma_optimum",
        "The Ornstein–Uhlenbeck sigma (>=0) applied to the fitness optimum at each generation",
        false, 0.0, "double", cmd};
    TCLAP::ValueArg<double> theta_optimum{"", "theta_optimum",
        "The Ornstein–Uhlenbeck theta (0<=theta<1) applied to the optimum at each generation",
        false, 0.0, "double", cmd};
    TCLAP::ValueArg<double> peakness{"", "peakness",
        "'alpha' parameter (peakness) of the fitness function "
        "(exp(-alpha*(phenotype^beta))",
        false, 1.0, "double", cmd};
    TCLAP::ValueArg<double> epistasis{"", "epistasis",
        "'beta' parameter (epistasis) of fitness "
        "function (exp(-alpha*(phenotype^beta))",
        false, 2.0, "double", cmd};


    OrnsteinUhlenbeckBiasModel get_model() {
        return OrnsteinUhlenbeckBiasModel(peakness.getValue(), epistasis.getValue(),
            bias_optimum.getValue(), bias_optimum_rate.getValue(), sigma_optimum.getValue(),
            theta_optimum.getValue());
    }
};
