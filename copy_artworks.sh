cp ./results/empiricalRevBayes_Mammal.pdf ./manuscript/figures/empirical_mammal.pdf
# Constant generation time
cp ./results/mammals_no/inference_RevBayes.boxplot.is_nuc.pdf ./manuscript/figures/simulation_cstGT_phylogram_support.pdf
cp ./results/mammals_no/inference_RevBayes.boxplot.is_OU.pdf ./manuscript/figures/simulation_cstGT_OU.pdf
# Variable generation time
cp ./results/mammals_slow/inference_RevBayes.boxplot.is_nuc.pdf ./manuscript/figures/simulation_varGT_phylogram_support.pdf
cp ./results/mammals_slow/inference_RevBayes.boxplot.is_OU.pdf ./manuscript/figures/simulation_varGT_OU.pdf
cp ./results/mammals_slow/inference_RevBayes.boxplot.var_multiplier.pdf ./manuscript/figures/simulation_varGT_MultiBM.pdf
cp ./results/mammals_slow/inference_RevBayes.boxplot.num_rate_changes.pdf ./manuscript/figures/simulation_varGT_MultiBM_NSwitch.pdf
cp ./results/mammals_slow/inference_RevBayes.boxplot.num_theta_changes.pdf ./manuscript/figures/simulation_varGT_MultiOU_NSwitch.pdf
cp ./results/mammals_slow/plot_RevBayesDistance_Phylo_neutral.pdf ./manuscript/figures/simulation_varGT_dist_NeutralPhylo.pdf
cp ./results/mammals_slow/plot_RevBayesDistance_Phylo_moving_optimum.pdf ./manuscript/figures/simulation_varGT_dist_SelectionPhylo.pdf
cp ./results/mammals_slow/plot_RevBayesDistance_Chrono_neutral.pdf ./manuscript/figures/simulation_varGT_dist_NeutralChrono.pdf
cp ./results/mammals_slow/plot_RevBayesDistance_Chrono_moving_optimum.pdf ./manuscript/figures/simulation_varGT_dist_SelectionChrono.pdf
