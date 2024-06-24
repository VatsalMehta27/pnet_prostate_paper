import re
import networkx as nx
from data.go_ontology.read_obo import read_obo
from config_path import GO_ONTOLOGY_PATH
from os.path import join
import pandas as pd


def add_edges(G, node, n_levels):
    edges = []
    source = node
    for l in range(n_levels):
        target = node + '_copy' + str(l + 1)
        edge = (source, target)
        source = target
        edges.append(edge)

    G.add_edges_from(edges)
    return G

def complete_network(G, n_leveles=4):
    sub_graph = nx.ego_graph(G, 'root', radius=n_leveles)
    terminal_nodes = [n for n, d in sub_graph.out_degree() if d == 0]
    distances = [len(nx.shortest_path(G, source='root', target=node)) for node in terminal_nodes]
    for node in terminal_nodes:
        distance = len(nx.shortest_path(sub_graph, source='root', target=node))
        if distance <= n_leveles:
            diff = n_leveles - distance + 1
            sub_graph = add_edges(sub_graph, node, diff)

    return sub_graph

def get_nodes_at_level(net, distance):
    # get all nodes within distance around the query node
    nodes = set(nx.ego_graph(net, 'root', radius=distance))

    # remove nodes that are not **at** the specified distance but closer
    if distance >= 1.:
        nodes -= set(nx.ego_graph(net, 'root', radius=distance - 1))

    return list(nodes)

def get_layers_from_net(net, n_levels):
    layers = []
    for i in range(n_levels):
        nodes = get_nodes_at_level(net, i)
        dict = {}
        for n in nodes:
            n_name = re.sub('_copy.*', '', n)
            next = net.successors(n)
            dict[n_name] = [re.sub('_copy.*', '', nex) for nex in next]
        layers.append(dict)
    return layers

class GO():

    def __init__(self):
        self.ontology = self.load_ontology()
        self.genes = self.load_genes()

    def load_genes(self):
        gene_to_go_term = pd.read_csv(join(GO_ONTOLOGY_PATH, "gene2go"), sep="\t")
        gene_to_go_term = gene_to_go_term[['#tax_id', 'GeneID', 'GO_ID']]
        
        gene_id_to_symbol = pd.read_csv(join(GO_ONTOLOGY_PATH, "gene_info"), sep="\t")
        gene_id_to_symbol = gene_id_to_symbol[['#tax_id', 'GeneID', 'Symbol']]
        
        combined = pd.merge(gene_to_go_term, gene_id_to_symbol, on=['#tax_id', 'GeneID'], how="inner")
        
        return combined

    def load_ontology(self):
        with open(join(GO_ONTOLOGY_PATH, "go-basic.obo")) as obofile:
            ontology = read_obo(obofile)
        return ontology

class GONetwork:
    def __init__(self):
        self.go = GO()
        self.netx = self.get_network()
        
    def get_network(self, edge_type="is_a"):
        if hasattr(self, 'netx'):
            return self.netx
        
        complete_graph = self.go.ontology
        filtered_graph = nx.DiGraph()
        
        for u, v, key in complete_graph.edges(keys=True):
            if key == edge_type:
                filtered_graph.add_edge(u, v, key=edge_type)
        
        filtered_graph.name = "go ontology"
        
        # add root node
        roots = [n for n, d in filtered_graph.in_degree() if d == 0]
        root_node = 'root'
        edges = [(root_node, n) for n in roots]
        filtered_graph.add_edges_from(edges)

        return filtered_graph
    
    def get_completed_network(self, n_levels):
        G = complete_network(self.netx, n_leveles=n_levels)
        return G

    def get_completed_tree(self, n_levels):
        G = self.get_tree()
        G = complete_network(G, n_leveles=n_levels)
        return G

    def get_layers(self, n_levels, direction='root_to_leaf'):
        if direction == 'root_to_leaf':
            net = self.get_completed_network(n_levels)
            layers = get_layers_from_net(net, n_levels)
        else:
            net = self.get_completed_network(5)
            layers = get_layers_from_net(net, 5)
            layers = layers[5 - n_levels:5]

        # get the last layer (genes level)
        terminal_nodes = [n for n, d in net.out_degree() if d == 0]  # set of terminal go terms
        # we need to find genes belonging to these go terms
        genes_df = self.go.genes

        go_term_gene_mapping = {}
        missing_ontology_terms = []
        print len(terminal_nodes)
        
        for term in terminal_nodes:
            term = re.sub('_copy.*', '', term)
            
            genes = genes_df[genes_df['GO_ID'] == term]['Symbol'].unique()
            
            if len(genes) == 0:
                missing_ontology_terms.append(term)
                
            go_term_gene_mapping[term] = genes
            
        print "num missing terms: ", str(len(missing_ontology_terms))

        layers.append(go_term_gene_mapping)
        
        return layers