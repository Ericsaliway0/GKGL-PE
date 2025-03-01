import pandas as pd
import networkx as nx
class Network:
    def __init__(self, p_value_df, kge):
        self.p_value_df = p_value_df
        self.kge = kge
        self.graph_nx = nx.DiGraph()  # Directed graph for protein network

    def load_protein_network(self, protein_file):
        """
        Loads a protein network from a CSV file.
        
        The CSV file should contain columns:
        - 'protein1': identifier for the first protein in the interaction
        - 'protein2': identifier for the second protein in the interaction
        - 'shared_partners_count': number of shared interaction partners
        - 'p-value': statistical p-value of the interaction
        - 'significance': significance level of the interaction
        """
        protein_df = pd.read_csv(protein_file)
        for _, row in protein_df.iterrows():
            p1, p2 = row['protein1'], row['protein2']
            shared_count = row['shared_partners_count']
            p_value = row['p-value']
            significance = row['significance']
            
            # Add protein nodes with attributes, where `stId` is equivalent to `protein2`
            if p1 not in self.graph_nx:
                self.graph_nx.add_node(p1, type='protein', stId=p1, weight=p_value, significance=significance)
            if p2 not in self.graph_nx:
                self.graph_nx.add_node(p2, type='protein', stId=p2, weight=p_value, significance=significance)
            
            # Add edge with attributes
            self.graph_nx.add_edge(p1, p2, shared_count=shared_count, weight=p_value)

    def set_significance_by_p_value(self, threshold=0.05):
        """
        Updates node significance based on a p-value threshold.
        
        Parameters:
        - threshold (float): significance threshold for p-values.
        """
        for node, attrs in self.graph_nx.nodes(data=True):
            if 'weight' in attrs:
                attrs['significance'] = 'significant' if attrs['weight'] < threshold else 'not-significant'

    def display_network_info(self):
        """
        Displays basic information about the network.
        """
        print("Number of nodes:", self.graph_nx.number_of_nodes())
        print("Number of edges:", self.graph_nx.number_of_edges())
        print("Node data:", dict(self.graph_nx.nodes(data=True)))
        print("Edge data:", dict(self.graph_nx.edges(data=True)))



'''protein_network = ProteinNetwork()
protein_network.load_protein_network('path/to/protein_network.csv')
protein_network.display_network_info()'''
