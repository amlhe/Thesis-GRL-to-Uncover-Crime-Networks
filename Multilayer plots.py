import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import networkx as nx
import matplotlib.patheffects as path_effects
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import scipy.stats as stats
from scipy.stats import gaussian_kde


# Read in adjacency matrix with network data
adj_MN = pd.read_csv('datasets/adj_MN.csv')
adj_PC = pd.read_csv('datasets/adj_PC.csv')
adj_SN = pd.read_csv('datasets/adj_SN.csv')
adj_WR = pd.read_csv('datasets/adj_WR.csv')
adj_AW = pd.read_csv('datasets/adj_AW.csv')
adj_JU = pd.read_csv('datasets/adj_JU.csv')

adj_MN.index = range(1, len(adj_MN.index) + 1)
adj_MN.columns = range(1, len(adj_MN.index) + 1)

adj_PC.index = range(1, len(adj_PC.index) + 1)
adj_PC.columns = range(1, len(adj_PC.index) + 1)

adj_SN.index = range(1, len(adj_SN.index) + 1)
adj_SN.columns = range(1, len(adj_SN.index) + 1)

adj_WR.index = range(1, len(adj_WR.index) + 1)
adj_WR.columns = range(1, len(adj_WR.index) + 1)

adj_AW.index = range(1, len(adj_AW.index) + 1)
adj_AW.columns = range(1, len(adj_AW.index) + 1)

adj_JU.index = range(1, len(adj_JU.index) + 1)
adj_JU.columns = range(1, len(adj_JU.index) + 1)

cols = ['steelblue', 'darksalmon', 'mediumseagreen']

# Set random seed for reproducability
np.random.seed(1)

# Defining networkx graph objects for the multilayer networks 
G_MN = nx.from_pandas_adjacency(adj_MN)
G_PC = nx.from_pandas_adjacency(adj_PC)

G_SN = nx.from_pandas_adjacency(adj_SN)

G_WR = nx.from_pandas_adjacency(adj_WR)
G_AW = nx.from_pandas_adjacency(adj_AW)
G_JU = nx.from_pandas_adjacency(adj_JU)

# Define network node positions
pos_MN = nx.spring_layout(G_MN) 
pos_PC = nx.spring_layout(G_PC) 
pos_WR = nx.spring_layout(G_WR) 
pos_AW = nx.spring_layout(G_AW) 
pos_JU = nx.spring_layout(G_JU) 

# Change between considered networks
graphs = [G_WR, G_AW, G_JU]
#graphs = [G_SN]
#graphs = [G_WR, G_AW, G_JU]
#graph_names = ["Meetings", "Phonecalls"]
#graph_names = ["Infinito"]
graph_names = ["Wiretap Records", "Arrest Warrants", "Judgments"]


# Plotting features
w = 10
h = 8

fig, ax = plt.subplots(1, 1, figsize=(w,h), dpi=200, subplot_kw={'projection':'3d'})

for gi, G in enumerate(graphs):
    if gi == 0:
        pos = pos_WR
    elif gi == 1:
        pos = pos_AW
    else:
        pos = pos_JU
    # node positions
    xs = list(list(zip(*list(pos.values())))[0])
    ys = list(list(zip(*list(pos.values())))[1])
    zs = [gi]*len(xs) # set a common z-position of the nodes 

    # node colors
    cs = [cols[gi]]*len(xs)
    
    # add within-layer edges 
    lines3d = [(list(pos[i])+[gi],list(pos[j])+[gi]) for i,j in G.edges()]
    line_collection = Line3DCollection(lines3d, zorder=gi-1, color=cols[gi], alpha=0.8)
    ax.add_collection3d(line_collection)
    
    # now add nodes
    ax.scatter(xs, ys, zs, c=cs, s=20, edgecolors='.2', marker='o', alpha=1, zorder=gi+1)

    
    # add a plane to designate the layer
    xdiff = max(xs)-min(xs)
    ydiff = max(ys)-min(ys)
    ymin = min(ys)-ydiff*0.1
    ymax = max(ys)+ydiff*0.1
    xmin = min(xs)-xdiff*0.1 * (w/h)
    xmax = max(xs)+xdiff*0.1 * (w/h)
    xx, yy = np.meshgrid([xmin, xmax],[ymin, ymax])
    zz = np.zeros(xx.shape)+gi
    ax.plot_surface(xx, yy, zz, color=cols[gi], alpha=0.1, zorder=gi)

    # add label
    layertext = ax.text(0.0, 1.1, gi*0.95+0.5, "Layer "+(graph_names[gi]),
                    color='.95', fontsize='large', zorder=1e5, ha='left', va='center',
                    path_effects=[path_effects.Stroke(linewidth=2, foreground=cols[gi]),
                                  path_effects.Normal()])

# set them all at the same x,y,zlims
ax.set_ylim(min(ys)-ydiff*0.1,max(ys)+ydiff*0.1)
ax.set_xlim(min(xs)-xdiff*0.1,max(xs)+xdiff*0.1)
ax.set_zlim(-0.1, len(graphs) - 1 + 0.1)

# select viewing angle
angle = 65
height_angle = 20
ax.view_init(height_angle, angle)

# how much do you want to zoom into the fig
ax.dist = 9.0

ax.set_axis_off()

plt.show()




# Create the aggregate networks of Montagna and Oversize
G_Montagna = nx.Graph()
G_Montagna.add_nodes_from(G_MN)
G_Montagna.add_nodes_from(G_PC)
G_Montagna.add_edges_from(G_MN.edges)
G_Montagna.add_edges_from(G_PC.edges)

G_Oversize = nx.Graph()
G_Oversize.add_nodes_from(G_WR)
G_Oversize.add_nodes_from(G_AW)
G_Oversize.add_nodes_from(G_JU)
G_Oversize.add_edges_from(G_WR.edges)
G_Oversize.add_edges_from(G_AW.edges)
G_Oversize.add_edges_from(G_JU.edges)

# Plot the aggregate networks
pos = nx.spring_layout(G_Montagna)
nx.draw(G_Montagna, pos, with_labels=True, node_size=300, node_color="skyblue", font_size=10, font_color="black")
plt.title("Aggregate Network of Montagna")
plt.show()

pos = nx.spring_layout(G_Oversize)
nx.draw(G_Oversize, pos, with_labels=True, node_size=300, node_color="skyblue", font_size=10, font_color="black")
plt.title("Aggregate Network of Oversize")
plt.show()


# Calculate the degree of each node
degrees_MN = dict(G_MN.degree)

# Determine the sizes of the nodes based on their degree
node_sizes_MN = [v * 50 for v in degrees_MN.values()]

# Plot the graph with node sizes based on degree
pos = nx.spring_layout(G_MN)  # positions for all nodes
nx.draw(G_MN, pos, with_labels=True, node_size=node_sizes_MN, node_color="skyblue", font_size=10, font_color="black")
plt.title("Graph with Node Sizes Based on Degree")
plt.show()




# Function to create histogram
def plot_degree_histogram(G):
    # Calculate the degrees
    degrees = [G.degree(n) for n in G.nodes()]

    # Calculate the average degree
    avg_degree = sum(degrees) / len(degrees)

    # Create a histogram of the degree distribution
    plt.hist(degrees, bins=max(degrees) - min(degrees) + 1, color='lightblue', edgecolor='black', alpha=0.7)

    # Add a red vertical line for the mean degree
    plt.axvline(x=avg_degree, color='red', linestyle='dashed', linewidth=2, label=f'Mean Degree: {avg_degree:.2f}')

    # Add labels and title
    plt.title('Degree Distribution Histogram')
    plt.xlabel('Degree')
    plt.ylabel('Density')
    plt.legend()
    plt.show()

plot_degree_histogram(G_SN)
