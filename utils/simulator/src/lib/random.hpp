#pragma once

#include <random>

// Random generator engine with seed 0.
u_long seed{0};
std::default_random_engine generator(seed);
std::default_random_engine generator_pop_size(seed);
std::default_random_engine generator_mut_rate(seed);

// Random distributions.
// Normal distribution with mean 0 and standard deviation 1.
std::normal_distribution<double> normal_distrib(0.0, 1.0);
// Uniform distribution between 0 and 1.
std::uniform_real_distribution<double> uniform_distrib(0.0, 1.0);
// Bernoulli distribution with probability 0.5.
std::bernoulli_distribution bernouilli_distrib(0.5);