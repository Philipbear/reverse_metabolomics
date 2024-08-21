from matchms import calculate_scores, Spectrum
from matchms.similarity import ModifiedCosine
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pyteomics import mzml


def create_ms2_network(feature_df, df, mzml_file, out_dir):
    """
    Create a network of MS/MS spectra with nodes colored as in the scatter plot
    and edges based on modified cosine similarity > 0.7 and at least 6 matched peaks.
    """
    # Read mzML file
    with mzml.MzML(mzml_file) as reader:
        spectra_dict = {spec['index']: spec for spec in reader if spec['ms level'] == 2}

    # Combine feature_df and df, ensuring 'Selected' status is prioritized
    df['status'] = df['discard_reason'].apply(get_status)
    feature_df['status'] = 'm/z unmatched to target list'
    combined_df = pd.concat([df, feature_df])

    combined_df = combined_df[~pd.isnull(combined_df['best_MS2_scan_idx'])].reset_index(drop=True)

    # Sort by status, first 'Selected', then others
    combined_df['status_sort'] = combined_df['status'].map({'Selected': 0}).fillna(1)
    combined_df = combined_df.sort_values(['status_sort', 'status']).drop('status_sort', axis=1)
    combined_df['best_MS2_scan_idx'] = combined_df['best_MS2_scan_idx'].astype(int)
    # dereplicate
    combined_df = combined_df.drop_duplicates(subset='best_MS2_scan_idx', keep='first')

    # Create Spectrum objects with preprocessed spectra
    spectra = []
    for idx, row in combined_df.iterrows():
        scan = spectra_dict.get(row['best_MS2_scan_idx'] - 1)
        if scan:
            mz = scan['m/z array']
            intensities = scan['intensity array']
            precursor_mz = scan['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0][
                'selected ion m/z']

            # Preprocess the spectrum
            mz, intensities = preprocess_spectrum(mz, intensities, precursor_mz)

            # Create Spectrum object
            spectrum = Spectrum(mz=mz,
                                intensities=intensities,
                                metadata={"precursor_mz": precursor_mz,
                                          "RT": scan['scanList']['scan'][0]['scan start time'],
                                          "spectrum_id": row['best_MS2_scan_idx'],
                                          "status": row['status']})
            spectra.append(spectrum)

    # Calculate scores
    print("Calculating mod cos scores...")
    modified_cosine = ModifiedCosine(tolerance=0.02)
    scores = calculate_scores(spectra, spectra, modified_cosine)

    # Create network
    G = nx.Graph()

    # Add all nodes first
    for spectrum in spectra:
        G.add_node(spectrum.get('spectrum_id'), status=spectrum.get('status'))

    # Add edges
    for (i, score) in enumerate(scores):
        if score[0].get('spectrum_id') <= score[1].get('spectrum_id'):
            continue
        if score[2][0] >= 0.7 and score[2][1] >= 6:
            node1, node2 = score[0].get('spectrum_id'), score[1].get('spectrum_id')
            G.add_edge(node1, node2, weight=score[2][0])

    print("Number of nodes in the graph:", len(G.nodes()))
    print("Number of edges in the graph:", len(G.edges()))

    # Define node colors
    def get_node_color(node_id):
        status = G.nodes[node_id]['status']
        return distinct_colors.get(status, distinct_colors['Other'])

    # Set node colors
    node_colors = [get_node_color(node) for node in G.nodes()]

    # Create the plot
    plt.figure(figsize=(20, 15))  # Increased figure size for better visibility

    # Use a combination of layouts for connected and isolated nodes
    connected_nodes = [node for node in G.nodes() if G.degree(node) > 0]
    isolated_nodes = list(nx.isolates(G))

    # Position connected nodes using Fruchterman-Reingold layout
    pos_connected = nx.fruchterman_reingold_layout(G.subgraph(connected_nodes), k=2, iterations=50)

    # Position isolated nodes in a grid layout around the connected component
    pos_isolated = nx.spring_layout(G.subgraph(isolated_nodes), k=5, iterations=50)

    # Combine layouts
    pos = {**pos_connected, **pos_isolated}

    # Adjust positions to spread isolated nodes
    max_x = max(p[0] for p in pos.values())
    max_y = max(p[1] for p in pos.values())
    for node in isolated_nodes:
        pos[node] = (pos[node][0] * 1.5 + max_x, pos[node][1] * 1.5 + max_y)

    # Draw edges
    nx.draw_networkx_edges(G, pos, alpha=0.2, edge_color='gray')

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=30)

    # Add labels for legend
    unique_statuses = set(nx.get_node_attributes(G, 'status').values())
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', label=status,
                                  markerfacecolor=distinct_colors.get(status, distinct_colors['Other']), markersize=5)
                       for status in unique_statuses]
    plt.legend(handles=legend_elements, title="Status", loc="center left", bbox_to_anchor=(1, 0.5))

    plt.title(f'MS/MS Network - {mzml_file.split("/")[-1]}')
    plt.axis('off')

    # Save the plot
    output_path = f"{out_dir}/{mzml_file.split('/')[-1]}_ms2_network_with_isolated.svg"
    plt.savefig(output_path, format='svg', bbox_inches='tight')
    plt.close()

    print(f"MS/MS network plot with isolated nodes saved to {output_path}")
    print(f"Number of nodes: {G.number_of_nodes()}")
    print(f"Number of edges: {G.number_of_edges()}")


def preprocess_spectrum(mz, intensities, precursor_mz):
    """
    Preprocess the spectrum by removing ions with m/z larger than precursor - 0.5
    and intensities less than 1% of the base peak.
    """
    # Remove ions with m/z larger than precursor - 0.5
    mask = mz < (precursor_mz - 0.5)
    mz = mz[mask]
    intensities = intensities[mask]

    # Remove ions with intensities less than 1% of the base peak
    max_intensity = np.max(intensities)
    intensity_threshold = 0.01 * max_intensity
    mask = intensities > intensity_threshold
    mz = mz[mask]
    intensities = intensities[mask]

    return mz, intensities


# Helper function to determine status (same as in the scatter plot function)
def get_status(reasons):
    if pd.isnull(reasons) or reasons == '':
        return 'Selected'
    reasons_list = reasons.split(';')
    return 'Selected' if 'Selected' in reasons_list else reasons_list[0]


# Define distinct colors (same as in the scatter plot function)
distinct_colors = {
    'Selected': '#1f77b4',
    'No core adduct': '#d62728',
    'MS2 explanation below cutoff': '#2ca02c',
    'Matched to doubly charged ion': '#ff7f0e',
    'Other': '#9467bd',
    'm/z unmatched to target list': '#D3D3D3'
}
