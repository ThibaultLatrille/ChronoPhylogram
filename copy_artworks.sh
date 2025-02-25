EXT=pdf
# Constant generation time
rm -rf ./manuscript/figures/simulation_*
cp ./results/mammals_no/inference_RevBayes.boxplot.is_nuc.${EXT} ./manuscript/figures/simulation_cstGT_phylogram_support.${EXT}
cp ./results/mammals_no/inference_RevBayes.boxplot.is_OU.${EXT} ./manuscript/figures/simulation_cstGT_OU.${EXT}
# Variable generation time
cp ./results/mammals_slow/inference_RevBayes.boxplot.is_nuc.${EXT} ./manuscript/figures/simulation_varGT_phylogram_support.${EXT}
cp ./results/mammals_slow/inference_RevBayes.boxplot.is_OU.${EXT} ./manuscript/figures/simulation_varGT_OU.${EXT}
cp ./results/mammals_slow/inference_RevBayes.boxplot.var_multiplier.${EXT} ./manuscript/figures/simulation_varGT_MultiBM.${EXT}
cp ./results/mammals_slow/inference_RevBayes.boxplot.num_rate_changes.${EXT} ./manuscript/figures/simulation_varGT_MultiBM_NSwitch.${EXT}
cp ./results/mammals_slow/inference_RevBayes.boxplot.num_theta_changes.${EXT} ./manuscript/figures/simulation_varGT_MultiOU_NSwitch.${EXT}
cp ./results/mammals_slow/plot_RevBayesDistance_Phylo_neutral.${EXT} ./manuscript/figures/simulation_varGT_dist_NeutralPhylo.${EXT}
cp ./results/mammals_slow/plot_RevBayesDistance_Phylo_moving_optimum.${EXT} ./manuscript/figures/simulation_varGT_dist_SelectionPhylo.${EXT}
cp ./results/mammals_slow/plot_RevBayesDistance_Chrono_neutral.${EXT} ./manuscript/figures/simulation_varGT_dist_NeutralChrono.${EXT}
cp ./results/mammals_slow/plot_RevBayesDistance_Chrono_moving_optimum.${EXT} ./manuscript/figures/simulation_varGT_dist_SelectionChrono.${EXT}
