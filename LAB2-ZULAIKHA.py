import requests
import pandas as pd
import networkx as nx
import streamlit as st
import matplotlib.pyplot as plt

# Step 1: Functions to Retrieve Data

# a. BioGRID retrieval function
def retrieve_ppi_biogrid(EGFR):
    biogrid_url = "https://webservice.thebiogrid.org/interactions"
    params = {
        "accessKey": "5ad3f0aea24a429da6fa81594be47320",  # Replace with your BioGRID API key
        "format": "json",
        "searchNames": True,
        "geneList": EGFR,
        "organism": 9606,  # Human species ID
        "includeInteractors": True
    }
    response = requests.get(biogrid_url, params=params)
    network_data = response.json()
    df = pd.DataFrame.from_dict(network_data, orient='index')
    
    # Check if 'OFFICIAL_SYMBOL_A' and 'OFFICIAL_SYMBOL_B' exist
    if 'OFFICIAL_SYMBOL_A' in df.columns and 'OFFICIAL_SYMBOL_B' in df.columns:
        df['OFFICIAL_SYMBOL_A'] = df['OFFICIAL_SYMBOL_A'].str.upper()
        df['OFFICIAL_SYMBOL_B'] = df['OFFICIAL_SYMBOL_B'].str.upper()
    else:
        st.warning("Columns 'OFFICIAL_SYMBOL_A' or 'OFFICIAL_SYMBOL_B' not found in the data.")
        return pd.DataFrame()  # Return an empty DataFrame or handle accordingly

    return df[['OFFICIAL_SYMBOL_A', 'OFFICIAL_SYMBOL_B']]

# b. STRING retrieval function
def retrieve_ppi_string(target_protein):
    string_url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": "EGFR",
        "species": 9606
    }
    response = requests.get(string_url, params=params)
    network_data = response.json()
    df = pd.json_normalize(network_data)
    return df[['preferredName_A', 'preferredName_B']]

# Step 2: Network Generation Function
def generate_network(dataframe):
    network_graph = nx.from_pandas_edgelist(dataframe, source=dataframe.columns[0], target=dataframe.columns[1])
    return network_graph

# Step 3: Centrality Measures Function
def get_centralities(network_graph):
    degree_centrality = nx.degree_centrality(network_graph)
    closeness_centrality = nx.closeness_centrality(network_graph)
    betweenness_centrality = nx.betweenness_centrality(network_graph)
    eigenvector_centrality = nx.eigenvector_centrality(network_graph, max_iter=1000)
    pagerank_centrality = nx.pagerank(network_graph)
    
    centralities = {
        'Degree Centrality': degree_centrality,
        'Closeness Centrality': closeness_centrality,
        'Betweenness Centrality': betweenness_centrality,
        'Eigenvector Centrality': eigenvector_centrality,
        'PageRank Centrality': pagerank_centrality
    }
    return centralities

# Step 4: Streamlit App for User Interaction
st.title("Human Protein-Protein Interaction (PPI) Network")
target_protein = st.text_input("Enter Protein ID (e.g., TP53): EGFR")
database_choice = st.selectbox("Select Database", ("BioGRID", "STRING"))

if st.button("Retrieve PPI Data"):
    # Retrieve data based on user choice
    if database_choice == "BioGRID":
        ppi_data = retrieve_ppi_biogrid(target_protein)
    else:
        ppi_data = retrieve_ppi_string(target_protein)
    
    # Check if data is empty
    if ppi_data.empty:
        st.warning("No PPI data found. Please check the protein ID and try again.")
    else:
        # Generate network graph
        network_graph = generate_network(ppi_data)
        
        # Display PPI data information in the first column
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("PPI Data Information")
            st.dataframe(ppi_data)
            st.write("Number of nodes:", network_graph.number_of_nodes())
            st.write("Number of edges:", network_graph.number_of_edges())
            
            # Visualize the network
            if network_graph.number_of_nodes() == 0:
                st.warning("No interactions found for the given protein.")
            else:
                fig, ax = plt.subplots(figsize=(10, 10))  # Adjust size for larger graphs
                layout = nx.spring_layout(network_graph, seed=42, k=0.15)  # Adjust 'k' for better spacing
                nx.draw(network_graph, layout, ax=ax, with_labels=True, node_size=30, node_color="skyblue", font_size=8)
                st.pyplot(fig)

        # Display centrality measures in the second column
        with col2:
            st.subheader("Centrality Measures")
            centralities = get_centralities(network_graph)
            for centrality_name, centrality_values in centralities.items():
                sorted_nodes = sorted(centrality_values.items(), key=lambda x: -x[1])[:5]  # Top 5 nodes for each centrality
                st.write(f"{centrality_name}:")
                for node, centrality in sorted_nodes:
                    st.write(f"{node}: {centrality:.4f}")
