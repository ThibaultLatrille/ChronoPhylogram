EXT=png
#/Users/tlatrille/Documents/venv/py312stats/bin/python3 /Users/tlatrille/Documents/ChronoPhylogram/scripts/plot_simulations_RevBayes.py --folder /Users/tlatrille/Documents/ChronoPhylogram/data_RevBayes/mammals_no --output /Users/tlatrille/Documents/ChronoPhylogram/results/mammals_no/inference_RevBayes.tsv
#/Users/tlatrille/Documents/venv/py312stats/bin/python3 /Users/tlatrille/Documents/ChronoPhylogram/scripts/plot_simulations_RevBayes.py --folder /Users/tlatrille/Documents/ChronoPhylogram/data_RevBayes/mammals_slow --output /Users/tlatrille/Documents/ChronoPhylogram/results/mammals_slow/inference_RevBayes.tsv
rm -rf ./manuscript/figures/simulation_*
# Constant generation time
cp ./results/mammals_no/inference_RevBayes.boxplot.is_nuc.${EXT} ./manuscript/figures/simulation_cstGT_phylogram_support.${EXT}
cp ./results/mammals_no/inference_RevBayes.boxplot.is_OU.${EXT} ./manuscript/figures/simulation_cstGT_OU.${EXT}
# Variable generation time
cp ./results/mammals_slow/inference_RevBayes.boxplot.is_nuc.${EXT} ./manuscript/figures/simulation_varGT_phylogram_support.${EXT}
cp ./results/mammals_slow/inference_RevBayes.boxplot.is_OU.${EXT} ./manuscript/figures/simulation_varGT_OU.${EXT}
cp ./results/mammals_slow/inference_RevBayes.boxplot.var_multiplier.${EXT} ./manuscript/figures/simulation_varGT_MultiBM.${EXT}
cp ./results/mammals_slow/inference_RevBayes.boxplot.num_rate_changes.${EXT} ./manuscript/figures/simulation_varGT_MultiBM_NSwitch.${EXT}
cp ./results/mammals_slow/inference_RevBayes.boxplot.num_theta_changes.${EXT} ./manuscript/figures/simulation_varGT_MultiOU_NSwitch.${EXT}
cp ./results/mammals_slow/plot_RevBayesDistance_Phylo_neutral.pdf ./manuscript/figures/simulation_varGT_dist_NeutralPhylo.pdf
cp ./results/mammals_slow/plot_RevBayesDistance_Phylo_moving_optimum.pdf ./manuscript/figures/simulation_varGT_dist_SelectionPhylo.pdf
cp ./results/mammals_slow/plot_RevBayesDistance_Chrono_neutral.pdf ./manuscript/figures/simulation_varGT_dist_NeutralChrono.pdf
cp ./results/mammals_slow/plot_RevBayesDistance_Chrono_moving_optimum.pdf ./manuscript/figures/simulation_varGT_dist_SelectionChrono.pdf
