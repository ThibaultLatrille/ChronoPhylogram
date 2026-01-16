EXT=png
rm -rf ./manuscript/figures/simulation_*
# Constant generation time
cp ./results/simulated_mammals_no/inference_RevBayes.boxplot.is_nuc.${EXT} ./manuscript/figures/simulation_cstGT_phylogram_support.${EXT}
cp ./results/simulated_mammals_no/inference_RevBayes.boxplot.is_OU.${EXT} ./manuscript/figures/simulation_cstGT_OU.${EXT}
# Variable generation time
cp ./results/simulated_mammals_slow/inference_RevBayes.boxplot.is_nuc.${EXT} ./manuscript/figures/simulation_varGT_phylogram_support.${EXT}
cp ./results/simulated_mammals_slow/inference_RevBayes.boxplot.is_OU.${EXT} ./manuscript/figures/simulation_varGT_OU.${EXT}
cp ./results/simulated_mammals_slow/inference_RevBayes.boxplot.var_multiplier.${EXT} ./manuscript/figures/simulation_varGT_MultiBM.${EXT}
cp ./results/simulated_mammals_slow/inference_RevBayes.boxplot.num_rate_changes.${EXT} ./manuscript/figures/simulation_varGT_MultiBM_NSwitch.${EXT}
cp ./results/simulated_mammals_slow/inference_RevBayes.boxplot.num_theta_changes.${EXT} ./manuscript/figures/simulation_varGT_MultiOU_NSwitch.${EXT}
cp ./results/simulated_mammals_slow/plot_RevBayesDistance_Phylo_neutral.pdf ./manuscript/figures/simulation_varGT_dist_NeutralPhylo.pdf
cp ./results/simulated_mammals_slow/plot_RevBayesDistance_Phylo_moving_optimum.pdf ./manuscript/figures/simulation_varGT_dist_SelectionPhylo.pdf
cp ./results/simulated_mammals_slow/plot_RevBayesDistance_Chrono_neutral.pdf ./manuscript/figures/simulation_varGT_dist_NeutralChrono.pdf
cp ./results/simulated_mammals_slow/plot_RevBayesDistance_Chrono_moving_optimum.pdf ./manuscript/figures/simulation_varGT_dist_SelectionChrono.pdf
cp ./results/simulated_mammals_slow/plot_MaximumLikelihood.pdf ./manuscript/figures/simulation_varGT_ML_vs_Bayes.pdf
cp ./results/simulated_mammals_slow/plot_MaximumLikelihood_logit.pdf ./manuscript/figures/simulation_varGT_ML_vs_Bayes_logit.pdf
# Saturation
cp ./results/simulated_mammals_slow/plot_distance_Phylo_neutral_saturation.pdf ./manuscript/figures/simulation_varGT_dist_NeutralPhylo_saturation.pdf
cp ./results/simulated_mammals_slow_low_u/plot_distance_Phylo_neutral_saturation.pdf ./manuscript/figures/simulation_varGT_dist_NeutralPhylo_slowU_saturation.pdf
cp ./results/simulated_mammals_slow_high_u/plot_distance_Phylo_neutral_saturation.pdf ./manuscript/figures/simulation_varGT_dist_NeutralPhylo_fastU_saturation.pdf
cp ./results/simulated_mammals_slow_very_high_u/plot_distance_Phylo_neutral_saturation.pdf ./manuscript/figures/simulation_varGT_dist_NeutralPhylo_veryFastU_saturation.pdf
cp ./results/simulated_mammals_slow_low_u/inference_RevBayes.boxplot.is_OU.${EXT} ./manuscript/figures/simulation_varGT_slowU_OU.${EXT}
cp ./results/simulated_mammals_slow_high_u/inference_RevBayes.boxplot.is_OU.${EXT} ./manuscript/figures/simulation_varGT_fastU_OU.${EXT}
cp ./results/simulated_mammals_slow_very_high_u/inference_RevBayes.boxplot.is_OU.${EXT} ./manuscript/figures/simulation_varGT_veryFastU_OU.${EXT}