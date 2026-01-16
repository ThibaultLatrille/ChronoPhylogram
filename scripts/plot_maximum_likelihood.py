import os
import argparse
from collections import defaultdict
from libraries_plot import *
from maximum_likelihood import *


def annotate_leaves_with_traits(tree: Tree, trait_leaves: dict[str, float], trait_name: str):
    for leaf in tree.get_leaves():
        if leaf.name in trait_leaves:
            setattr(leaf, trait_name, trait_leaves[leaf.name])
        else:
            raise ValueError(f"Leaf {leaf.name} not found in trait_leaves")


def main(chronogram_path: str, phylogram_path: str, rb_results_path: str, simu_trait_paths: list[str],
         output_path: str):
    chronogram = Tree(chronogram_path, format=1)
    phylogram = Tree(phylogram_path, format=1)

    output_dico = defaultdict(list)
    for simu_trait_path in simu_trait_paths:
        simu_name = os.path.basename(os.path.dirname(simu_trait_path))
        seed = "seed" + os.path.basename(simu_trait_path).split("_seed")[1].split(".")[0]
        print(f"Processing simulation: {simu_name}, seed: {seed}")
        traits = pd.read_csv(simu_trait_path, sep="\t", comment="#")
        trait_leaves = {row['TaxonName']: row['Phenotype_mean'] for _, row in traits.iterrows()}
        annotate_leaves_with_traits(chronogram, trait_leaves, "trait")
        annotate_leaves_with_traits(phylogram, trait_leaves, "trait")
        anc_z_chrono, var_chrono, ll_chrono = brownian_fitting(chronogram, "trait")
        anc_z_phylo, var_phylo, ll_phylo = brownian_fitting(phylogram, "trait")
        AIC_chrono = 2 * 2 - 2 * ll_chrono  # 2 parameters: anc_z and var
        AIC_phylo = 2 * 2 - 2 * ll_phylo
        Delta_AIC = AIC_chrono - AIC_phylo
        AIC_weight_phylo = 1 - np.exp(-0.5 * Delta_AIC) / (1 + np.exp(-0.5 * Delta_AIC))
        # Save to a dataframe
        output_dico["simu"].append(simu_name)
        output_dico["seed"].append(seed)
        output_dico["anc_z_chronogram"].append(anc_z_chrono)
        output_dico["var_chronogram"].append(var_chrono)
        output_dico["ll_chronogram"].append(ll_chrono)
        output_dico["anc_z_phylogram"].append(anc_z_phylo)
        output_dico["var_phylogram"].append(var_phylo)
        output_dico["ll_phylogram"].append(ll_phylo)
        output_dico["Delta_AIC"].append(Delta_AIC)
        output_dico["AIC_weight_phylogram"].append(AIC_weight_phylo)
    df_output = pd.DataFrame(output_dico)
    print(df_output.head())

    # Join with the RevBayes results
    rb_results = pd.read_csv(rb_results_path, sep="\t")
    # Filter only to "simple_BM_Switchnodes" model
    if "is_nuc" not in rb_results.columns:
        # Plot histogram of AIC weights to see available models
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        df_output.hist(column="AIC_weight_phylogram", by="simu", ax=ax)
        plt.tight_layout()
        plt.savefig(output_path)
        plt.close("all")
        return

    rb_results = rb_results[rb_results["model"] == "simple_BM_Switchnodes"][["simu", "seed", "is_nuc"]]
    print(f"Joining {len(df_output)} ML results with {len(rb_results)} RevBayes results")
    df_merged = pd.merge(df_output, rb_results, on=["simu", "seed"], how="left")
    df_merged.to_csv(output_path, sep="\t", index=False)
    assert len(df_merged) == len(df_output), "Merged dataframe has different length than ML results"
    # Clip Delta_AIC to -40, 40 for better visualization
    df_merged["Delta_AIC"] = df_merged["Delta_AIC"].clip(-40, 40)

    # Plot the results
    for t, l, f, m in [("", "", lambda x: x, 0.5), ("_logit", " (logit)", lambda x: np.log(x / (1 - x)), 0.0)]:
        # Scatter plot of Delta_AIC vs support for phylogram from RevBayes
        fig, ax = plt.subplots(1, 1, figsize=(7, 5))
        # Show in different colors the different simulations in the scatter plot, add legend
        # group by simulation
        for i, (simu, df_simu) in enumerate(df_merged.groupby("simu")):
            color, label = cs_simu_models[gr_simu_models(simu)]
            ax.scatter(f(df_simu["AIC_weight_phylogram"]), f(df_simu["is_nuc"]), label=label, color=color, alpha=0.7)
        # regression line
        x, y = f(df_merged["AIC_weight_phylogram"]), f(df_merged["is_nuc"])
        mask = np.isfinite(x) & np.isfinite(y)
        x, y = x[mask], y[mask]
        a, b = np.polyfit(x, y, 1)
        x_lim = np.array([min(0, min(x)), max(x)])
        r2 = np.corrcoef(x, y)[0, 1] ** 2
        ax.plot(x_lim, a * x_lim + b, color="black", linestyle="--", label=f"Regression line (R² = {r2:.4f})")
        ax.set_xlabel(f"AIC weight support for Phylogram{l}", fontsize=14)
        ax.set_ylabel(f"Bayesian posterior support for Phylogram{l}", fontsize=14)
        ax.axhline(m, color="grey", linestyle="--", label="Equal support lines")
        ax.axvline(m, color="grey", linestyle="--")
        ax.legend()
        plt.tight_layout()
        plt.savefig(output_path.replace(".pdf", f"{t}.pdf"))
        plt.close("all")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--chronogram", help="Input chronogram tree file", required=True)
    parser.add_argument("--phylogram", help="Input phylogram tree file", required=True)
    parser.add_argument("--rb_results", help="Input rb_results file", required=True)
    parser.add_argument("--simu_traits", help="Input simulated tree file", required=True, nargs="+")
    parser.add_argument("--output", help="Output file", required=True)
    args = parser.parse_args()
    main(args.chronogram, args.phylogram, args.rb_results, args.simu_traits, args.output)
