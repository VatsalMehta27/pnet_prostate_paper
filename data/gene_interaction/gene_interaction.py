import pandas as pd

from config_path import GENE_INTERACTION_PATH

def read_gene_interaction_data():
    # TODO: edit this for dataset
    gene_interactions = pd.read_csv(GENE_INTERACTION_PATH, sep="\t")

    return gene_interactions


def get_gene_interaction_layer(genes):
    mapping = {}

    gene_interactions = read_gene_interaction_data()

    missing_gene_interactions = []

    for gene in genes:
        interacts_with_gene = gene_interactions[gene_interactions["interactor_A"] == gene]["interactor_B"].unique()

        if len(interacts_with_gene) == 0:
            missing_gene_interactions.append(gene)

        mapping[gene] = interacts_with_gene

    print len(missing_gene_interactions), " genes have no interactors"

    return mapping