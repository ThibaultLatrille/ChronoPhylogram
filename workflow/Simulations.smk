import os
import numpy as np

FOLDER = os.path.abspath('.')

executable_list = ["nodetraits", "readnodetraits"]
exec_dico = {}

for executable in executable_list:
    exec_path = ""
    for b in ["BayesCode", "bayescode", "Bayescode", "bayesCode", "BAYESCODe", "BayesCODE"]:
        exec_path = os.path.join(FOLDER,f'utils/{b}/bin/{executable}')
        if os.path.exists(exec_path):
            break
    if not os.path.exists(exec_path):
        # Find executable in the path using whereis. If not found, raise an error.
        split = os.popen(f'whereis {executable}').read().split()
        if len(split) > 1:
            exec_path = split[1].strip()
        else:
            raise FileNotFoundError(f'{executable} not found. Please install BayesCode and add it to your path.')
    exec_dico[executable] = exec_path
    print(f"Found {executable} at {exec_path}")

configfile: 'config/config.yaml'

ChronoGram = f"{FOLDER}/utils/simulator/{config['tree']}"
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

rule all:
    input:
        expand(f"{FOLDER}/results/simu_{{simulator}}.pdf",simulator=config["simulators"])


def variance_env(nbr_loci, a, mut_rate, pop_size, h2):
    vG = 4 * float(a) ** 2 * float(nbr_loci) * float(mut_rate) * int(pop_size)
    h2 = float(h2)
    assert vG > 0
    assert 0 <= h2 <= 1
    return vG * (1 - h2) / h2


rule run_simulations:
    input:
        exec=lambda wildcards: f"{FOLDER}/utils/simulator/build/{config['simulators'][wildcards.simulator]['model']}",
        tree=ChronoGram
    output:
        nhx=f"{FOLDER}/data_simulated/{{simulator}}/replicate_seed{{seed}}.nhx.gz"
    params:
        file=lambda wildcards: f"{FOLDER}/data_simulated/{wildcards.simulator}/replicate_seed{wildcards.seed}",
        folder=lambda wildcards: f"{FOLDER}/data_simulated/{wildcards.simulator}",
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
        tree=f"{FOLDER}/utils/simulator/{config['tree']}"
    output:
        nhx=f"{FOLDER}/data_simulated/neutral_tree.nhx.gz"
    params:
        file=f"{FOLDER}/data_simulated/neutral_tree",
        folder=lambda wildcards: f"{FOLDER}/data_simulated",
        simulator=SIMULATOR_PARAMS["core"]
    shell:
        'mkdir -p {params.folder} && {input.exec} --tree {input.tree} {params.simulator} --variance_environment 0.0 --number_loci 50000 --seed 42 --seed_pop_size {SEED_POP_SIZE} --seed_mut_rate {SEED_MUT_RATE} --output {params.file} && gzip {params.file}.nhx'

rule simulation_traits:
    input:
        script=f"{FOLDER}/scripts/pre_processed_traits_simulations.py",
        nhx=rules.run_simulations.output.nhx
    output:
        traits=f"{FOLDER}/data_simulated/{{simulator}}/replicate_{{gram}}_seed{{seed}}.traits.tsv"
    shell:
        'python3 {input.script} --input {input.nhx} --traitsfile {output.traits}'

rule distance_tree:
    input:
        script=f"{FOLDER}/scripts/neutral_tree.py",
        nhx=rules.run_simulations.output.nhx
    output:
        tree=f"{FOLDER}/data_simulated/{{simulator}}/replicate_seed{{seed}}.d.tree"
    shell:
        'python3 {input.script} --neutral_tree {input.nhx} --tree {output.tree} '


rule neutral_tree:
    input:
        script=f"{FOLDER}/scripts/neutral_tree.py",
        nhx=rules.run_neutral_simulation.output.nhx
    output:
        tree=f"{FOLDER}/data_simulated/neutral_tree.tree"
    shell:
        'python3 {input.script} --neutral_tree {input.nhx} --tree {output.tree} '


rule scale_tree:
    input:
        script=f"{FOLDER}/scripts/scale_tree.py",
        tree_1=rules.neutral_tree.output.tree,
        tree_2=ChronoGram
    output:
        tree_1=f"{FOLDER}/data_simulated/Phylo_scaled.tree",
        tree_2=f"{FOLDER}/data_simulated/Chrono_scaled.tree"
    shell:
        'python3 {input.script} --tree_1 {input.tree_1} --tree_2 {input.tree_2} --tree_output_1 {output.tree_1} --tree_output_2 {output.tree_2}'

rule bayescode_inference:
    input:
        exec=exec_dico['nodetraits'],
        tree=f"{FOLDER}/data_simulated/{{gram}}_scaled.tree",
        traits=rules.simulation_traits.output.traits
    output:
        run=f"{FOLDER}/data_simulated/{{simulator}}/inference_{{gram}}_seed{{seed}}.run"
    params:
        chain=lambda wildcards: f"{FOLDER}/data_simulated/{wildcards.simulator}/inference_{wildcards.gram}_seed{wildcards.seed}",
        until=f"-u {config['bayes_until']}"
    shell:
        '{input.exec} {params.until} --uniq_kappa --df 1 --tree {input.tree} --traitsfile {input.traits} {params.chain}'

rule plot_trait_distance:
    input:
        scripts=f"{FOLDER}/scripts/plot_trait_distance.py",
        distance_tree=f"{FOLDER}/data_simulated/{{gram}}_scaled.tree",
        simu_tree=rules.run_simulations.output.nhx
    output:
        run=f"{FOLDER}/data_simulated/{{simulator}}/plot_{{gram}}_seed{{seed}}.pdf"
    shell:
        'python3 {input.scripts} --distance_tree {input.distance_tree} --simu_tree {input.simu_tree} --output {output.run}'

rule merge_simulated_trace:
    input:
        script=f"{FOLDER}/scripts/merge_results_simulations.py",
        bayescode=expand(f"{FOLDER}/data_simulated/{{{{simulator}}}}/inference_{{gram}}_seed{{seed}}.run",
            gram=["Phylo", "Chrono"],seed=SEEDS),
        pdf=expand(f"{FOLDER}/data_simulated/{{{{simulator}}}}/plot_{{gram}}_seed{{seed}}.pdf",
            gram=["Phylo", "Chrono"],seed=SEEDS),
        tree=expand(f"{FOLDER}/data_simulated/{{{{simulator}}}}/replicate_seed{{seed}}.d.tree",seed=SEEDS)
    output:
        tsv=f"{FOLDER}/results/simu_{{simulator}}.pdf"
    shell:
        'python3 {input.script} --bayescode {input.bayescode} --output {output.tsv}'
