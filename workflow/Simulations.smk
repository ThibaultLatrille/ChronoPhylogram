import os
import numpy as np

FOLDER = os.path.abspath('.')

executable_list = ["nodetraits", "readnodetraits", "rb"]
exec_dico = {}

for executable in executable_list:
    exec_path = ""
    for b in ["BayesCode", "bayescode", "Bayescode", "bayesCode", "BAYESCODe", "BayesCODE"]:
        exec_path = os.path.join(FOLDER,f'utils/{b}/bin/{executable}')
        if os.path.exists(exec_path):
            break
    if executable == "rb":
        exec_path = "/opt/homebrew/Caskroom/miniforge/base/envs/osx-64/bin/rb"
    if not os.path.exists(exec_path):
        # Find executable in the path using whereis. If not found, raise an error.
        split = os.popen(f'whereis {executable}').read().split()
        if len(split) > 1:
            exec_path = split[1].strip()
        else:
            pg = "RevBayes" if executable == 'rb' else "BayesCode"
            raise FileNotFoundError(f'{executable} not found. Please install {pg} and add it to your path.')
    exec_dico[executable] = exec_path
    print(f"Found {executable} at {exec_path}")

configfile: 'config/config.yaml'

EXP = config["experiment_name"]
SEEDS = [int(i) for i in np.linspace(config["min_seed"],config["max_seed"],config["nb_seeds"])]
SIMULATOR_PARAMS = {s: " ".join(['--{0} {1}'.format(k,v) for k, v in d.items() if k != "model"]) for s, d in
                    config["simulators"].items()}
for tag in ["population_size", "mutation_rate"]:
    tag_sigma = f"{tag}_brownian_sigma"
    if f"scaled_{tag_sigma}" not in config:
        continue
    print(f"Order of magnitude in change of {tag}: {config[f'scaled_{tag_sigma}']}")
    config["core"][f"{tag_sigma}"] = config[f"scaled_{tag_sigma}"] / np.sqrt(config["core"]['number_of_generations'])
    print(f"{tag_sigma}: {config['core'][f'{tag_sigma}']}")
SIMULATOR_PARAMS["core"] = " ".join(['--{0} {1}'.format(k,v) for k, v in config["core"].items()])
SEED_POP_SIZE = 42
SEED_MUT_RATE = 24
GRAMS = ["Phylo", "Chrono"]

wildcard_constraints:
    simulator=r"[a-zA-Z_]+",gram=r"[a-zA-Z]+"

rule all:
    input:
        expand(f"{FOLDER}/results/{EXP}/simu_prob_{{simulator}}.pdf",simulator=config["simulators"]),
        expand(f"{FOLDER}/results/{EXP}/simu_wAIC_{{simulator}}.pdf",simulator=config["simulators"]),
        expand(f"{FOLDER}/results/{EXP}/branch_wAIC_{{simulator}}.pdf",simulator=config["simulators"]),
        expand(f"{FOLDER}/results/{EXP}/plot_ancestral_{{gram}}_{{simulator}}.pdf",
            gram=GRAMS,simulator=config["simulators"]),
        expand(f"{FOLDER}/results/{EXP}/plot_distance_{{gram}}_{{simulator}}.pdf",
            gram=GRAMS,simulator=config["simulators"]),
        expand(f"{FOLDER}/results/{EXP}/plot_distance_reconstructed_{{gram}}_{{simulator}}.pdf",
            gram=GRAMS,simulator=config["simulators"]),
        f"{FOLDER}/results/{EXP}/plot_div_comparison.pdf",
        f"{FOLDER}/results/{EXP}/inference_RevBayes.tsv"


def variance_env(nbr_loci, a, mut_rate, pop_size, h2):
    vG = 4 * float(a) ** 2 * float(nbr_loci) * float(mut_rate) * int(pop_size)
    h2 = float(h2)
    assert vG > 0
    assert 0 <= h2 <= 1
    return vG * (1 - h2) / h2


rule prepare_chronogram:
    input:
        script=f"{FOLDER}/scripts/prepare_chronogram.py",
        tree=f"{FOLDER}/{config['tree']}"
    output:
        tree=f"{FOLDER}/data_simulated/{EXP}/chronogram.tree"
    shell:
        'python3 {input.script} --input {input.tree} --output {output.tree}'

rule run_simulations:
    input:
        exec=lambda wildcards: f"{FOLDER}/utils/simulator/build/{config['simulators'][wildcards.simulator]['model']}",
        tree=rules.prepare_chronogram.output.tree
    output:
        nhx=f"{FOLDER}/data_simulated/{EXP}/{{simulator}}/replicate_seed{{seed}}.nhx.gz"
    params:
        file=lambda wildcards: f"{FOLDER}/data_simulated/{EXP}/{wildcards.simulator}/replicate_seed{wildcards.seed}",
        folder=lambda wildcards: f"{FOLDER}/data_simulated/{EXP}/{wildcards.simulator}",
        simulator=lambda wildcards: "{0} {1}".format(SIMULATOR_PARAMS["core"],SIMULATOR_PARAMS[wildcards.simulator]),
        ve=lambda wildcards: variance_env(
            nbr_loci=config["simulators"][wildcards.simulator]["number_loci"],
            a=config["core"]["mutation_mean_effect_size"],
            mut_rate=config["core"]["mutation_rate_per_loci_per_generation"],
            pop_size=config["core"]["population_size"],
            h2=config["heritability"])
    shell:
        'mkdir -p {params.folder} && {input.exec} --tree {input.tree} {params.simulator} --variance_environment {params.ve} --seed {wildcards.seed} --seed_pop_size {SEED_POP_SIZE} --seed_mut_rate {SEED_MUT_RATE} --output {params.file} && gzip {params.file}.nhx'

rule run_neutral_simulation:
    input:
        exec=lambda wildcards: f"{FOLDER}/utils/simulator/build/neutral",
        tree=rules.prepare_chronogram.output.tree
    output:
        nhx=f"{FOLDER}/data_simulated/{EXP}/neutral_tree.nhx.gz"
    params:
        file=f"{FOLDER}/data_simulated/{EXP}/neutral_tree",
        folder=lambda wildcards: f"{FOLDER}/data_simulated",
        simulator=SIMULATOR_PARAMS["core"]
    shell:
        'mkdir -p {params.folder} && {input.exec} --tree {input.tree} {params.simulator} --variance_environment 0.0 --number_loci 50000 --seed 42 --seed_pop_size {SEED_POP_SIZE} --seed_mut_rate {SEED_MUT_RATE} --output {params.file} && gzip {params.file}.nhx'

rule simulation_traits:
    input:
        script=f"{FOLDER}/scripts/pre_processed_traits_simulations.py",
        nhx=rules.run_simulations.output.nhx
    output:
        traits=f"{FOLDER}/data_simulated/{EXP}/{{simulator}}/replicate_seed{{seed}}.traits.tsv"
    shell:
        'python3 {input.script} --input {input.nhx} --traitsfile {output.traits}'

rule distance_tree:
    input:
        script=f"{FOLDER}/scripts/neutral_tree.py",
        nhx=rules.run_simulations.output.nhx
    output:
        tree=f"{FOLDER}/data_simulated/{EXP}/{{simulator}}/replicate_seed{{seed}}.d.tree"
    shell:
        'python3 {input.script} --neutral_tree {input.nhx} --tree {output.tree} '

rule plot_div_comparison:
    input:
        scripts=f"{FOLDER}/scripts/plot_div_comparison.py",
        distance_tree=expand(f"{FOLDER}/data_simulated/{EXP}/{{simulator}}/replicate_seed{{seed}}.d.tree",
            seed=SEEDS,simulator=config["simulators"]),
    output:
        pdf=f"{FOLDER}/results/{EXP}/plot_div_comparison.pdf"
    shell:
        'python3 {input.scripts} --distance_tree {input.distance_tree} --output {output.pdf}'

rule neutral_tree:
    input:
        script=f"{FOLDER}/scripts/neutral_tree.py",
        nhx=rules.run_neutral_simulation.output.nhx
    output:
        tree=f"{FOLDER}/data_simulated/{EXP}/neutral_tree.tree"
    shell:
        'python3 {input.script} --neutral_tree {input.nhx} --tree {output.tree} '

rule reconstructed_chronogram:
    input:
        script=f"{FOLDER}/scripts/calib_pl.R",
        tree=rules.neutral_tree.output.tree
    output:
        tree=f"{FOLDER}/data_simulated/{EXP}/reconstructed_chronogram.tree"
    shell:
        'Rscript {input.script} {input.tree} {output.tree}'

rule scale_tree:
    input:
        script=f"{FOLDER}/scripts/scale_tree.py",
        tree_1=rules.neutral_tree.output.tree,
        tree_2=rules.reconstructed_chronogram.output.tree
    output:
        tree_1=f"{FOLDER}/data_simulated/{EXP}/Phylo_scaled.tree",
        tree_2=f"{FOLDER}/data_simulated/{EXP}/Chrono_scaled.tree"
    shell:
        'python3 {input.script} --tree_1 {input.tree_1} --tree_2 {input.tree_2} --tree_output_1 {output.tree_1} --tree_output_2 {output.tree_2}'

rule plot_trait_distance:
    input:
        scripts=f"{FOLDER}/scripts/plot_trait_distance.py",
        distance_tree=f"{FOLDER}/data_simulated/{EXP}/{{gram}}_scaled.tree",
        annotated_tree=expand(f"{FOLDER}/data_simulated/{EXP}/{{{{simulator}}}}/replicate_seed{{seed}}.nhx.gz",seed=SEEDS)
    output:
        run=f"{FOLDER}/results/{EXP}/plot_distance_{{gram}}_{{simulator}}.pdf"
    shell:
        'python3 {input.scripts} --distance_tree {input.distance_tree} --annotated_tree {input.annotated_tree} --output {output.run}'

rule bayescode_inference:
    input:
        exec=exec_dico['nodetraits'],
        tree=f"{FOLDER}/data_simulated/{EXP}/{{gram}}_scaled.tree",
        traits=rules.simulation_traits.output.traits
    output:
        run=f"{FOLDER}/data_BayesCode/{EXP}/{{simulator}}/inference_{{gram}}_seed{{seed}}.run"
    params:
        chain=lambda wildcards: f"{FOLDER}/data_BayesCode/{EXP}/{wildcards.simulator}/inference_{wildcards.gram}_seed{wildcards.seed}",
        until=f"-u {config['bayes_until']}"
    shell:
        '{input.exec} {params.until} --uniq_kappa --df 1 --tree {input.tree} --traitsfile {input.traits} {params.chain}'

rule bayescode_read:
    input:
        exec=exec_dico['readnodetraits'],
        run=rules.bayescode_inference.output.run
    output:
        tree=f"{FOLDER}/data_BayesCode/{EXP}/{{simulator}}/inference_{{gram}}_seed{{seed}}.Phenotype_mean.nhx",
        wAIC=f"{FOLDER}/data_BayesCode/{EXP}/{{simulator}}/inference_{{gram}}_seed{{seed}}.wAIC.tsv",
    params:
        chain=rules.bayescode_inference.params.chain,
        until=rules.bayescode_inference.params.until,
        burnin=f"-b {config['bayes_burn_in']}"
    shell:
        'if [ -f {params.chain}.trace.gz ]; then gunzip -f {params.chain}.trace.gz; fi; if [ -f {params.chain}.chain.gz ]; then gunzip -f {params.chain}.chain.gz; fi;'
        '{input.exec} {params.until} {params.burnin} --wAIC {params.chain}; {input.exec} {params.until} {params.burnin} --newick {params.chain}; gzip -f {params.chain}.trace; gzip -f {params.chain}.chain'

rule plot_ancestral_traits:
    input:
        script=f"{FOLDER}/scripts/plot_ancestral_traits.py",
        inference_tree=expand(f"{FOLDER}/data_BayesCode/{EXP}/{{{{simulator}}}}/inference_{{{{gram}}}}_seed{{seed}}.Phenotype_mean.nhx",seed=SEEDS),
        simu_tree=expand(f"{FOLDER}/data_simulated/{EXP}/{{{{simulator}}}}/replicate_seed{{seed}}.nhx.gz",seed=SEEDS),
    output:
        pdf=f"{FOLDER}/results/{EXP}/plot_ancestral_{{gram}}_{{simulator}}.pdf"
    shell:
        'python3 {input.script} --inference_tree {input.inference_tree} --simu_tree {input.simu_tree} --output {output.pdf}'

rule plot_ancestral_trait_distance:
    input:
        scripts=f"{FOLDER}/scripts/plot_trait_distance.py",
        distance_tree=f"{FOLDER}/data_simulated/{EXP}/{{gram}}_scaled.tree",
        annotated_tree=expand(f"{FOLDER}/data_BayesCode/{EXP}/{{{{simulator}}}}/inference_{{{{gram}}}}_seed{{seed}}.Phenotype_mean.nhx",seed=SEEDS)
    output:
        run=f"{FOLDER}/results/{EXP}/plot_distance_reconstructed_{{gram}}_{{simulator}}.pdf"
    shell:
        'python3 {input.scripts} --distance_tree {input.distance_tree} --annotated_tree {input.annotated_tree} --output {output.run}'


rule merge_simulated_trace:
    input:
        script=f"{FOLDER}/scripts/merge_results_simulations.py",
        traces=expand(f"{FOLDER}/data_BayesCode/{EXP}/{{{{simulator}}}}/inference_{{gram}}_seed{{seed}}.Phenotype_mean.nhx",gram=GRAMS,seed=SEEDS),
        bayescode=expand(f"{FOLDER}/data_BayesCode/{EXP}/{{{{simulator}}}}/inference_{{gram}}_seed{{seed}}.run",gram=GRAMS,seed=SEEDS)
    output:
        tsv=f"{FOLDER}/results/{EXP}/simu_prob_{{simulator}}.pdf"
    shell:
        'python3 {input.script} --bayescode {input.bayescode} --output {output.tsv}'

rule merge_simulated_wAIC:
    input:
        script=f"{FOLDER}/scripts/merge_wAIC_simulations.py",
        bayescode=expand(f"{FOLDER}/data_BayesCode/{EXP}/{{{{simulator}}}}/inference_{{gram}}_seed{{seed}}.wAIC.tsv",gram=GRAMS,seed=SEEDS)
    output:
        tsv=f"{FOLDER}/results/{EXP}/simu_wAIC_{{simulator}}.pdf"
    shell:
        'python3 {input.script} --bayescode {input.bayescode} --output {output.tsv}'

rule plot_branch_wAIC:
    input:
        script=f"{FOLDER}/scripts/plot_branch_wAIC.py",
        tree_x=rules.scale_tree.output.tree_2,
        tree_y=rules.scale_tree.output.tree_1,
        bayescode=expand(f"{FOLDER}/data_BayesCode/{EXP}/{{{{simulator}}}}/inference_{{gram}}_seed{{seed}}.wAIC.tsv",gram=GRAMS,seed=SEEDS)
    output:
        tsv=f"{FOLDER}/results/{EXP}/branch_wAIC_{{simulator}}.pdf"
    shell:
        'python3 {input.script} --tree_x {input.tree_x} --tree_y {input.tree_y} --bayescode {input.bayescode} --output {output.tsv}'


rule convert_to_RevBayes:
    input:
        script=f"{FOLDER}/scripts/convert_to_nexus.py",
        tree=f"{FOLDER}/data_simulated/{EXP}/{{gram}}_scaled.tree",
        traits=rules.simulation_traits.output.traits
    output:
        nexus_tree=f"{FOLDER}/data_RevBayes/{EXP}/{{simulator}}/inference_{{gram}}_seed{{seed}}/tree.nex",
        nexus_traits=f"{FOLDER}/data_RevBayes/{EXP}/{{simulator}}/inference_{{gram}}_seed{{seed}}/traits.nex"
    shell:
        'python3 {input.script} --input_tree {input.tree} --input_traits {input.traits} --output_tree {output.nexus_tree} --output_traits {output.nexus_traits}'

rule run_RevBayes:
    input:
        exec=exec_dico['rb'],
        rev_file=f"{FOLDER}/scripts/mcmc_{{rb}}.Rev",
        tree=rules.convert_to_RevBayes.output.nexus_tree,
        traits=rules.convert_to_RevBayes.output.nexus_traits
    output:
        log=f"{FOLDER}/data_RevBayes/{EXP}/{{simulator}}/inference_{{gram}}_seed{{seed}}/{{rb}}.log"
    params:
        folder=lambda wildcards: f"{FOLDER}/data_RevBayes/{EXP}/{wildcards.simulator}/inference_{wildcards.gram}_seed{wildcards.seed}"
    shell:
        'cd {params.folder} && {input.exec} {input.rev_file}'

rule run_Both_RevBayes:
    input:
        exec=exec_dico['rb'],
        rev_file=f"{FOLDER}/scripts/mcmc_{{rb}}.Rev",
        timetree=f"{FOLDER}/data_RevBayes/{EXP}/{{simulator}}/inference_Chrono_seed{{seed}}/tree.nex",
        nuctree=f"{FOLDER}/data_RevBayes/{EXP}/{{simulator}}/inference_Phylo_seed{{seed}}/tree.nex",
        traits=f"{FOLDER}/data_RevBayes/{EXP}/{{simulator}}/inference_Chrono_seed{{seed}}/traits.nex"
    output:
        log=f"{FOLDER}/data_RevBayes/{EXP}/{{simulator}}/inference_Both_seed{{seed}}/{{rb}}.log"
    params:
        folder=lambda wildcards: f"{FOLDER}/data_RevBayes/{EXP}/{wildcards.simulator}/inference_Both_seed{wildcards.seed}"
    shell:
        'mkdir -p {params.folder};'
        'cp {input.timetree} {params.folder}/tree_time.nex;'
        'cp {input.nuctree} {params.folder}/tree_nuc.nex;'
        'cp {input.traits} {params.folder}/traits.nex;'
        'cd {params.folder} && {input.exec} {input.rev_file}'

rule compress_RevBayes_BM_log:
    input:
        log=lambda wildcards: rules.run_Both_RevBayes.output.log if wildcards.gram == "Both" else rules.run_RevBayes.output.log
    output:
        log=f"{FOLDER}/data_RevBayes/{EXP}/{{simulator}}/inference_{{gram}}_seed{{seed}}/{{rb}}.log.gz"
    params:
        prefix=f"{FOLDER}/data_RevBayes/{EXP}/{{simulator}}/inference_{{gram}}_seed{{seed}}/{{rb}}"
    shell:
        'for f in {params.prefix}*.log; do if [ ! -f $f.gz ]; then gzip $f; fi; done'

rule gather_RevBayes_log:
    input:
        script=f"{FOLDER}/scripts/plot_simulations_RevBayes.py",
        rb_log=expand(rules.compress_RevBayes_BM_log.output.log,gram=GRAMS,seed=SEEDS,
            simulator=config["simulators"],rb=["simple_BM_model", "simple_OU_RJ", "relaxed_BM_RJ"]),
        rb_bm_log=expand(rules.compress_RevBayes_BM_log.output.log,gram=["Both"],seed=SEEDS,
            simulator=config["simulators"],rb=["simple_BM_Switch"])
    output:
        plot=f"{FOLDER}/results/{EXP}/inference_RevBayes.tsv"
    params:
        folder=f"{FOLDER}/data_RevBayes/{EXP}"
    shell:
        'python3 {input.script} --folder {params.folder} --output {output.plot}'
