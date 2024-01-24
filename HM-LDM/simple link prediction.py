import pandas as pd
import numpy as np
import random
import networkx as nx
from tqdm import tqdm
import re
import matplotlib.pyplot as plt
from node2vec import Node2Vec
import matplotlib.lines as mlines
from scipy.linalg import sqrtm

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import jaccard_score
from sklearn.cluster import KMeans
from sklearn.neighbors import kneighbors_graph
from sklearn.metrics import roc_curve, auc, precision_recall_curve

file_path = 'MN_edge.txt'

fb_df = pd.read_csv(file_path, sep=' ', header=None, names=['node_1', 'node_2'])

# create graph
G = nx.from_pandas_edgelist(fb_df, "node_1", "node_2", create_using=nx.Graph())

fb_df = pd.read_csv(file_path, sep=' ', header=None, names=['node_1', 'node_2'])
max_node = max(fb_df[['node_1', 'node_2']].values.flatten())
G = nx.Graph()
G.add_nodes_from(range(1, max_node + 1))
G.add_edges_from(fb_df[['node_1', 'node_2']].values)
nodes = sorted(G.nodes())
adj_G = nx.to_numpy_matrix(G, nodelist=nodes)

# combine all nodes in a list
node_list = pd.unique(fb_df['node_1']).tolist() + pd.unique(fb_df['node_2']).tolist()

# remove duplicate items from the list
node_list = list(dict.fromkeys(node_list))

# build adjacency matrix
adj_G = nx.to_numpy_array(G, nodelist = node_list)


# Get the Minimum Spanning Tree (MST)
mst = nx.minimum_spanning_tree(G)

# Create a graph from the MST edges
G_mst = nx.Graph()
G_mst.add_edges_from(mst.edges())

# Remove the edges that are not in the MST
G_removed_edges = G.copy()
G_removed_edges.remove_edges_from([edge for edge in G.edges() if edge not in G_mst.edges()])

# combine all nodes in a list
node_list = pd.unique(fb_df['node_1']).tolist() + pd.unique(fb_df['node_2']).tolist()

# remove duplicate items from the list
node_list = list(dict.fromkeys(node_list))

# build adjacency matrix
adj_G = nx.to_numpy_array(G, nodelist = node_list)

# get unconnected node-pairs
all_unconnected_pairs = []

# Iterate through the adjacency matrix
for i in range(adj_G.shape[0]):
    for j in range(i + 1, adj_G.shape[1]):  # Use i + 1 to avoid duplicate pairs and the diagonal
        if adj_G[i, j] == 0:
            all_unconnected_pairs.append([node_list[i], node_list[j]])


node_1_unlinked = [i[0] for i in all_unconnected_pairs]
node_2_unlinked = [i[1] for i in all_unconnected_pairs]

data = pd.DataFrame({'node_1':node_1_unlinked, 
                     'node_2':node_2_unlinked})

# add target variable 'link'
data['link'] = 0

initial_node_count = len(G.nodes)

fb_df_temp = fb_df.copy()

# list to store removable links
omissible_links_index = [i for i, edge in enumerate(fb_df[['node_1', 'node_2']].values) if tuple(edge) not in G_mst.edges()]

# Drop removable edges based on MST
half_omissible_links_index = omissible_links_index[:len(omissible_links_index)//2]
fb_df_ghost = fb_df.loc[half_omissible_links_index]

#fb_df_ghost = fb_df.loc[omissible_links_index]

# add the target variable 'link'
fb_df_ghost['link'] = 1

data = pd.concat([data, fb_df_ghost[['node_1', 'node_2', 'link']]], ignore_index=True)

# drop removable edges
fb_df_partial = fb_df.drop(index=fb_df_ghost.index.values)

# build graph
G_data = nx.from_pandas_edgelist(fb_df_partial, "node_1", "node_2", create_using=nx.Graph())

# Plot the original graph, MST, and the modified graph with removed edges
plt.figure(figsize=(15, 5))

# Plot Original Graph
plt.subplot(131)
pos = nx.spring_layout(G, seed=42)
nx.draw(G, with_labels=False, pos=pos, node_size=40, alpha=0.6, width=0.7)
plt.title('Original Graph')
plt.text(0.5, -0.1, f"Edges: {len(G.edges())}\nNon-Edges: {len(G.nodes()) * (len(G.nodes()) - 1) // 2 - len(G.edges())}", horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)

# Plot Minimum Spanning Tree
plt.subplot(132)
nx.draw(G_mst, with_labels=False, node_size=40, alpha=0.6, width=0.7)
plt.title('Minimum Spanning Tree')
plt.text(0.5, -0.1, f"Edges: {len(G_mst.edges())}\nNon-Edges: {len(G_mst.nodes()) * (len(G_mst.nodes()) - 1) // 2 - len(G_mst.edges())}", horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)

# Plot Modified Graph with Removed Edges
plt.subplot(133)
nx.draw(G_data, with_labels=False, node_size=40, alpha=0.6, width=0.7)
plt.title('Graph with Removed Edges')
plt.text(0.5, -0.1, f"Edges: {len(G_data.edges())}\nNon-Edges: {len(G_data.nodes()) * (len(G_removed_edges.nodes()) - 1) // 2 - len(G_removed_edges.edges())}", horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)

plt.show()



# Function to compute common neighbors
def common_neighbors(G, node1, node2):
    return len(list(nx.common_neighbors(G, node1, node2)))

# Function to compute the Jaccard score
def jaccard_similarity(G, node1, node2):
    neighbors_node1 = set(nx.neighbors(G, node1))
    neighbors_node2 = set(nx.neighbors(G, node2))
    intersection_size = len(neighbors_node1.intersection(neighbors_node2))
    union_size = len(neighbors_node1.union(neighbors_node2))
    jaccard_similarity = intersection_size / union_size if union_size != 0 else 0
    return jaccard_similarity

# Function to compute the preferential attachment
def preferential_attachment(G, node1, node2):
    degree_node1 = G.degree[node1]
    degree_node2 = G.degree[node2]
    return degree_node1 * degree_node2

# Function to compute Adamic-Adar
def adamic_adar(G, node1, node2):
    aa = sum(1 / np.log(len(list(nx.neighbors(G, common_neighbor)))) for common_neighbor in nx.common_neighbors(G, node1, node2))
    return aa


# Function to compute resource allocation
def resource_allocation(G, node1, node2):
    ra = sum(1 / (len(list(nx.neighbors(G, common_neighbor)))) for common_neighbor in nx.common_neighbors(G, node1, node2))
    return ra

def katz_index(G, node1, node2, alpha=0.01):
    try:
        # Calculate Katz centrality
        katz = nx.katz_centrality_numpy(G, alpha=alpha)
        # Return the sum of Katz centrality for the two nodes
        return katz[node1] + katz[node2]
    except Exception as e:
        print(f"Error calculating Katz index: {e}")
        return 0


# Function to compute the Personalized PageRank
def ppa(G, alpha):
    A = nx.to_numpy_array(G)
    D = np.diag(np.sum(A, axis=1))
    D_inv = np.linalg.inv(D)
    P = D_inv @ A
    num_rows, num_columns = A.shape
    ppr = np.linalg.inv(np.eye(num_rows) - alpha * P)
    ppr = ppr / np.sum(ppr, axis=1, keepdims=True)
    return ppr

# Function to compute features for link prediction
def compute_features(G, edges):
    features = []
    for edge in edges:
        node1, node2 = edge
        common_neighbors_val = common_neighbors(G, node1, node2)
        adamic_adar_val = adamic_adar(G, node1, node2)
        katz_index_val = katz_index(G, node1, node2)
        jaccard_similarity_val = jaccard_similarity(G, node1, node2)

        # Ensure the values are real (take the real part if complex)
        common_neighbors_val = common_neighbors_val.real if np.iscomplexobj(common_neighbors_val) else common_neighbors_val
        adamic_adar_val = adamic_adar_val.real if np.iscomplexobj(adamic_adar_val) else adamic_adar_val
        katz_index_val = katz_index_val.real if np.iscomplexobj(katz_index_val) else katz_index_val
        jaccard_similarity_val = jaccard_similarity_val.real if np.iscomplexobj(jaccard_similarity_val) else jaccard_similarity_val

        features.append([common_neighbors_val, adamic_adar_val, katz_index_val, jaccard_similarity_val])
    return np.array(features)

# Function to compute ROC curve and AUC for the binary classification problem
def compute_roc_curve(labels, scores, model_name):
    fpr, tpr, _ = roc_curve(labels, scores)
    roc_auc = auc(fpr, tpr)
    return fpr, tpr, roc_auc


# Create sparse_i_rem and sparse_j_rem for all links
sparse_i_rem = fb_df_ghost['node_1']
sparse_j_rem = fb_df_ghost['node_2']

# Create sparse_i and sparse_j for all links
sparse_i = fb_df_partial['node_1']
sparse_j = fb_df_partial['node_2']

# Create non_sparse_i and non_sparse_j for all links
num_nodes = adj_G.shape[0]
non_sparse_i, non_sparse_j = [], []
for i in range(num_nodes):
    for j in range(i + 1, num_nodes):
        if adj_G[i, j] == 0:
            non_sparse_i.append(i)
            non_sparse_j.append(j)
non_sparse_i = pd.DataFrame(np.array(non_sparse_i))
non_sparse_j = pd.DataFrame(np.array(non_sparse_j))

# Create non_sparse_i and non_sparse_j with the same number of elements as sparse_i_rem and sparse_j_rem
num_elements = len(sparse_i_rem)

# Choose random elements from the current content of non_sparse_i and non_sparse_j
random_indices = np.random.choice(len(non_sparse_i), num_elements, replace=False)

# Create non_sparse_i and non_sparse_j with the chosen random elements
non_sparse_i = pd.DataFrame(np.array(non_sparse_i)[random_indices])
non_sparse_j = pd.DataFrame(np.array(non_sparse_j)[random_indices])

# Positive samples (edges with links)
positive_samples_edges = pd.DataFrame({'node_1': sparse_i, 'node_2': sparse_j})
positive_samples_labels = pd.Series(np.ones(len(sparse_i)))  # Assuming 1 for links

# Missing links (edges with missing links)
missing_links_edges = pd.DataFrame({'node_1': sparse_i_rem, 'node_2': sparse_j_rem})
missing_links_labels = pd.Series(np.ones(len(sparse_i_rem)))  # Assuming 1 for missing links

# Negative samples (non-links)
negative_samples_edges = pd.DataFrame({'node_1': non_sparse_i.values.flatten(),'node_2': non_sparse_j.values.flatten()})
negative_samples_labels = pd.Series(np.zeros(len(non_sparse_i)), index=non_sparse_i.index)  # Assuming 0 for non-links

# Combine positive, missing links, and negative samples into a single dataset
all_samples = pd.concat([positive_samples_edges, missing_links_edges, negative_samples_edges], ignore_index=True)
all_labels = pd.concat([positive_samples_labels, missing_links_labels, negative_samples_labels], ignore_index=True)

# Split the data into training and testing sets
train_data, test_data, train_labels, test_labels = train_test_split(all_samples, all_labels, test_size=0.3, random_state=35, stratify=all_labels)

train_data = train_data[(train_data['node_1'] != 0) & (train_data['node_2'] != 0)]
train_labels = train_labels[train_data.index]
test_data = test_data[(test_data['node_1'] != 0) & (test_data['node_2'] != 0)]
test_labels = test_labels[test_data.index]

# Compute features for training and testing sets
features_train = compute_features(G, train_data[['node_1', 'node_2']].values)
features_test = compute_features(G, test_data[['node_1', 'node_2']].values)



# Train logistic regression model
lr = LogisticRegression(class_weight="balanced")
lr.fit(features_train, train_labels)

# Make predictions on the test set
predictions = lr.predict_proba(features_test)

# Evaluate the model using ROC-AUC score
roc_auc = roc_auc_score(test_labels, predictions[:, 1])
print("ROC-AUC Score:", roc_auc)


# Train logistic regression models for common neighbors, Adamic-Adar, Katz index, and Jaccard score
lr_common_neighbors = LogisticRegression(class_weight="balanced")
lr_common_neighbors.fit(features_train[:, 0].reshape(-1, 1), train_labels)
predictions_common_neighbors = lr_common_neighbors.predict_proba(features_test[:, 0].reshape(-1, 1))

lr_adamic_adar = LogisticRegression(class_weight="balanced")
lr_adamic_adar.fit(features_train[:, 1].reshape(-1, 1), train_labels)
predictions_adamic_adar = lr_adamic_adar.predict_proba(features_test[:, 1].reshape(-1, 1))

lr_katz_index = LogisticRegression(class_weight="balanced")
lr_katz_index.fit(features_train[:, 2].reshape(-1, 1), train_labels)
predictions_katz_index = lr_katz_index.predict_proba(features_test[:, 2].reshape(-1, 1))

lr_jaccard = LogisticRegression(class_weight="balanced")
lr_jaccard.fit(features_train[:, 3].reshape(-1, 1), train_labels)
predictions_jaccard = lr_jaccard.predict_proba(features_test[:, 3].reshape(-1, 1))

# Compute ROC curves for all models
fpr_common_neighbors, tpr_common_neighbors, roc_auc_common_neighbors = compute_roc_curve(test_labels, predictions_common_neighbors[:, 1], 'Common Neighbors')
fpr_adamic_adar, tpr_adamic_adar, roc_auc_adamic_adar = compute_roc_curve(test_labels, predictions_adamic_adar[:, 1], 'Adamic-Adar')
fpr_katz_index, tpr_katz_index, roc_auc_katz_index = compute_roc_curve(test_labels, predictions_katz_index[:, 1], 'Katz Index')
fpr_jaccard, tpr_jaccard, roc_auc_jaccard = compute_roc_curve(test_labels, predictions_jaccard[:, 1], 'Jaccard Similarity')


# Initialize the plot
plt.figure(figsize=(10, 6))

# Plot ROC curves
plt.plot(fpr_common_neighbors, tpr_common_neighbors, label=f'Common Neighbors (AUC = {roc_auc_common_neighbors:.2f})')
plt.plot(fpr_adamic_adar, tpr_adamic_adar, label=f'Adamic-Adar (AUC = {roc_auc_adamic_adar:.2f})')
plt.plot(fpr_katz_index, tpr_katz_index, label=f'Katz Index (AUC = {roc_auc_katz_index:.2f})')
plt.plot(fpr_jaccard, tpr_jaccard, label=f'Jaccard Similarity (AUC = {roc_auc_jaccard:.2f})')


# Plot random guessing line
plt.plot([0, 1], [0, 1], linestyle='--', color='gray', label='Random')

# Set labels and title
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend(loc='lower right')
plt.grid(True)
plt.show()











file_path = 'SN_edge.txt'
fb_df = pd.read_csv(file_path, sep=' ', header=None, names=['node_1', 'node_2'])
max_node = max(fb_df[['node_1', 'node_2']].values.flatten())
G = nx.Graph()
G.add_nodes_from(range(1, max_node + 1))
G.add_edges_from(fb_df[['node_1', 'node_2']].values)
nodes = sorted(G.nodes())
adj_G = nx.to_numpy_matrix(G, nodelist=nodes)

# combine all nodes in a list
node_list = pd.unique(fb_df['node_1']).tolist() + pd.unique(fb_df['node_2']).tolist()

# remove duplicate items from the list
node_list = list(dict.fromkeys(node_list))

# build adjacency matrix
adj_G = nx.to_numpy_array(G, nodelist = node_list)


num_edges = G.number_of_edges()
num_nodes = G.number_of_nodes()
isolated_nodes = [node for node, degree in G.degree() if degree == 0]
num_isolated_nodes = len(isolated_nodes)
num_components = nx.number_connected_components(G)
print("Number of Edges:", num_edges)
print("Number of Nodes:", num_nodes)
print("Number of Isolated Nodes:", num_isolated_nodes)
print("Number of Connected Components:", num_components)


# Save data for original data
# Create sparse_i and sparse_j for all links
sparse_i = fb_df['node_1']
sparse_j = fb_df['node_2']

# Create non_sparse_i and non_sparse_j for all links
num_nodes = adj_G.shape[0]
non_sparse_i, non_sparse_j = [], []
for i in range(num_nodes):
    for j in range(i + 1, num_nodes):
        if adj_G[i, j] == 0:
            non_sparse_i.append(i)
            non_sparse_j.append(j)
non_sparse_i = pd.DataFrame(np.array(non_sparse_i))
non_sparse_j = pd.DataFrame(np.array(non_sparse_j))

sparse_i.to_csv('datasets/Original data/sparse_i.txt', header=False, index=False)
sparse_j.to_csv('datasets/Original data/sparse_j.txt', header=False, index=False)
non_sparse_i.to_csv('datasets/Original data/non_sparse_i.txt', header=False, index=False)
non_sparse_j.to_csv('datasets/Original data/non_sparse_j.txt', header=False, index=False)


# Save data for negative samples
# fb_df_partial contains remaining graph after negative sample
# Create sparse_i_rem and sparse_j_rem for all links
sparse_i_rem = fb_df_ghost['node_1']
sparse_j_rem = fb_df_ghost['node_2']

# Create sparse_i and sparse_j for all links
sparse_i = fb_df_partial['node_1']
sparse_j = fb_df_partial['node_2']

# Create non_sparse_i and non_sparse_j for all links
num_nodes = adj_G.shape[0]
non_sparse_i, non_sparse_j = [], []
for i in range(num_nodes):
    for j in range(i + 1, num_nodes):
        if adj_G[i, j] == 0:
            non_sparse_i.append(i)
            non_sparse_j.append(j)
non_sparse_i = pd.DataFrame(np.array(non_sparse_i))
non_sparse_j = pd.DataFrame(np.array(non_sparse_j))

# Create non_sparse_i and non_sparse_j with the same number of elements as sparse_i_rem and sparse_j_rem
num_elements = len(sparse_i_rem)

# Choose random elements from the current content of non_sparse_i and non_sparse_j
random_indices = np.random.choice(len(non_sparse_i), num_elements, replace=False)

# Create non_sparse_i and non_sparse_j with the chosen random elements
non_sparse_i = pd.DataFrame(np.array(non_sparse_i)[random_indices])
non_sparse_j = pd.DataFrame(np.array(non_sparse_j)[random_indices])

# Save to .txt files 
sparse_i.to_csv('datasets/Negative samples/sparse_i.txt', header=False, index=False)
sparse_j.to_csv('datasets/Negative samples/sparse_j.txt', header=False, index=False)
non_sparse_i.to_csv('datasets/Negative samples/non_sparse_i.txt', header=False, index=False)
non_sparse_j.to_csv('datasets/Negative samples/non_sparse_j.txt', header=False, index=False)
sparse_i_rem.to_csv('datasets/Negative samples/sparse_i_rem.txt', header=False, index=False)
sparse_j_rem.to_csv('datasets/Negative samples/sparse_j_rem.txt', header=False, index=False)


# Save data for link prediction
num_nodes = adj_G.shape[0]
sparse_i, sparse_j, sparse_i_rem, sparse_j_rem = [], []
for i in range(num_nodes):
    for j in range(i + 1, num_nodes):
        if adj_G[i, j] == 1:
            sparse_i.append(i)
            sparse_j.append(j)
        if adj_G[i, j] == 0:
            sparse_i_rem.append(i)
            sparse_j_rem.append(j)
sparse_i = pd.DataFrame(np.array(sparse_i))
sparse_j = pd.DataFrame(np.array(sparse_j))
sparse_i_rem = pd.DataFrame(np.array(sparse_i_rem))
sparse_j_rem = pd.DataFrame(np.array(sparse_j_rem))
sparse_i.to_csv('sparse_i.txt', header=False, index=False)
sparse_j.to_csv('sparse_j.txt', header=False, index=False)
sparse_i_rem.to_csv('sparse_i_rem.txt', header=False, index=False)
sparse_j_rem.to_csv('sparse_j_rem.txt', header=False, index=False)













# Determine number of connected components
W = nx.adjacency_matrix(G).toarray()
D = np.diag(np.sum(W, axis=1))
L = D-W

# Calculate the different Laplacian matrices
L_unnormalized = L # unnormalized
L_normalized = np.identity(len(W)) - np.linalg.inv(D) @ W # normalized
L_symmetric = np.identity(len(W)) - np.diag(1 / np.sqrt(np.sum(W, axis=1))) @ W @ np.diag(1 / np.sqrt(np.sum(W, axis=1))) # normalized symmetric

eigenvalues_unnormalized = np.linalg.eigvals(L_unnormalized)
eigenvalues_normalized = np.linalg.eigvals(L_normalized)
eigenvalues_symmetric = np.linalg.eigvals(L_symmetric)

# Sort eigenvalues in ascending order
eigenvalues_unnormalized.sort()
eigenvalues_normalized.sort()
eigenvalues_symmetric.sort()

# Plot the 10 smallest eigenvalues
num_smallest_eigenvalues = 10

plt.figure(figsize=(15, 5))

# Plot for unnormalized Laplacian
plt.subplot(1, 3, 1)
plt.plot(eigenvalues_unnormalized[:num_smallest_eigenvalues], marker='o')
plt.title('Unnormalized Laplacian')
plt.xlabel('Eigenvalue Index')
plt.ylabel('Eigenvalue')

# Plot for normalized Laplacian
plt.subplot(1, 3, 2)
plt.plot(eigenvalues_normalized[:num_smallest_eigenvalues], marker='o')
plt.title('Normalized Laplacian')
plt.xlabel('Eigenvalue Index')
plt.ylabel('Eigenvalue')

# Plot for normalized symmetric Laplacian
plt.subplot(1, 3, 3)
plt.plot(eigenvalues_symmetric[:num_smallest_eigenvalues], marker='o')
plt.title('Normalized Symmetric Laplacian')
plt.xlabel('Eigenvalue Index')
plt.ylabel('Eigenvalue')

plt.tight_layout()
plt.show()

print(eigenvalues_normalized[:num_smallest_eigenvalues]-eigenvalues_symmetric[:num_smallest_eigenvalues])


# Create the graph and get the largest connected component
connected_components = list(nx.connected_components(G))
largest_connected_component = max(connected_components, key=len)
largest_connected_subgraph = G.subgraph(largest_connected_component)

# Compute the adjacency matrix of the largest connected component
adj_large = nx.adjacency_matrix(largest_connected_subgraph).toarray()

# Calculate the Laplacian matrices
D = np.diag(np.sum(adj_large, axis=1))
L_unnormalized = D - adj_large
L_normalized = np.identity(len(adj_large)) - np.linalg.inv(D) @ adj_large
L_symmetric = np.identity(len(adj_large)) - sqrtm(np.linalg.inv(D)) @ adj_large @ sqrtm(np.linalg.inv(D))

# Calculate the eigenvalues
eigenvalues_unnormalized = np.linalg.eigvals(L_unnormalized)
eigenvalues_normalized = np.linalg.eigvals(L_normalized)
eigenvalues_symmetric = np.linalg.eigvals(L_symmetric)

# Plot the 10 largest eigenvalues
num_eigenvalues = 10
fig, axs = plt.subplots(1, 3, figsize=(15, 5))

axs[0].plot(np.sort(eigenvalues_unnormalized)[-num_eigenvalues:], marker='o')
axs[0].set_title('Unnormalized Laplacian')

axs[1].plot(np.sort(eigenvalues_normalized)[-num_eigenvalues:], marker='o')
axs[1].set_title('Normalized Laplacian')

axs[2].plot(np.sort(eigenvalues_symmetric)[-num_eigenvalues:], marker='o')
axs[2].set_title('Normalized Symmetric Laplacian')

plt.show()

















import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import label_binarize
import numpy as np
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score, accuracy_score

def misclass_roc(memb, Y, a=1, b=1):
    z = pd.get_dummies(memb)
    H = z.shape[1]
    V = Y.shape[0]
    
    Abs_Freq = np.dot(z.T, np.dot(Y, z))
    np.fill_diagonal(Abs_Freq, np.diag(Abs_Freq) / 2)
    
    Tot = np.dot(z.T, np.dot(np.ones((V, V)), z))
    np.fill_diagonal(Tot, (np.diag(Tot) - np.bincount(memb - 1, minlength=H)) / 2)
    
    Rel_Freq = (a + Abs_Freq) / (a + b + Tot)
    pred = np.dot(z, np.dot(Rel_Freq, z.T))
    
    true_labels = Y[np.tril_indices(V, k=-1)]
    predicted_scores = pred[np.tril_indices(V, k=-1)]
    
    # ROC Curve
    fpr, tpr, _ = roc_curve(true_labels, predicted_scores)
    auc_roc = auc(fpr, tpr)
    
    # Precision-Recall Curve
    precision, recall, _ = precision_recall_curve(true_labels, predicted_scores)
    auc_pr = average_precision_score(true_labels, predicted_scores)
    
    # Accuracy
    predicted_labels = (predicted_scores >= 0.5).astype(int)
    acc = accuracy_score(true_labels, predicted_labels)
    
    return {'perf': {'fpr': fpr, 'tpr': tpr, 'precision': precision, 'recall': recall},
            'auc_roc': auc_roc, 'auc_pr': auc_pr, 'acc': acc}

# Load CSV files in Python
membership_1 = pd.read_csv("../ESBM ROC/memb_Z_DP_mis.csv")["x"]
membership_2 = pd.read_csv("../ESBM ROC/memb_Z_DM_mis.csv")["x"]
membership_3 = pd.read_csv("../ESBM ROC/memb_Z_PY_mis.csv")["x"]
membership_4 = pd.read_csv("../ESBM ROC/memb_Z_GN_mis.csv")["x"]
Y = pd.read_csv("../ESBM ROC/Y.csv").to_numpy()

# Plot ROC curves for each membership vector
plt.figure(figsize=(10, 6))

membership_list = [membership_1, membership_2, membership_3, membership_4]
auc_scores = []

for i, memb in enumerate(membership_list):
    result = misclass_roc(memb, Y)
    auc_score = result['auc_roc']
    auc_scores.append(auc_score)
    plt.plot(result['perf']['fpr'], result['perf']['tpr'], label=f'{["DP", "DM", "PY", "GN"][i]} (AUC = {auc_score:.2f})')

plt.plot([0, 1], [0, 1], 'k--', label='Random')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend(loc='lower right')
plt.grid(True)
plt.show()
