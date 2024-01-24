import argparse
import torch
import numpy as np
import torch.optim as optim
import sys
from tqdm import tqdm
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import pandas as pd

sys.path.append('./src/')

from HM_LDM import LSM


parser = argparse.ArgumentParser(description='Hybrid Membership-Latent Distance Model')

parser.add_argument('--epochs', type=int, default=10000, metavar='N',
                    help='number of epochs for training (default: 10K)')

parser.add_argument('--scaling_epochs', type=int, default=2000, metavar='N',
                    help='number of epochs for learning initial scale for the random effects (default: 2K)')

parser.add_argument('--delta_sq', type=float, default=1, metavar='N',
                    help='delta^2 hyperparameter controlling the volume of the simplex')

parser.add_argument('--cuda', type=eval, 
                      choices=[True, False],  default=True,
                    help='CUDA training')

parser.add_argument('--p', type=int, 
                      choices=[1, 2],  default=1,
                    help='L2 norm power (default: 1)')

parser.add_argument('--LP',type=eval, 
                      choices=[True, False], default=False,
                    help='performs link prediction')

parser.add_argument('--D', type=int, default=10, metavar='N',
                    help='dimensionality of the embeddings (default: 8)')

parser.add_argument('--lr', type=float, default=0.1, metavar='N',
                    help='learning rate for the ADAM optimizer, for large values of delta 0.01 is more stable (default: 0.1)')

parser.add_argument('--sample_percentage', type=float, default=1, metavar='N',
                    help='Sample size network percentage, it should be equal or less than 1 (default: 1)')

parser.add_argument('--dataset', type=str, default='SN',
                    help='dataset to apply HM-LDM on')


sys.argv = ['']
del sys
args, unknown = parser.parse_known_args()

device = torch.device("cpu")
#torch.set_default_tensor_type('torch.FloatTensor')
torch.set_default_dtype(torch.float32)
torch.set_default_device('cpu')


if __name__ == "__main__":
    latent_dims=[args.D]
    datasets=[args.dataset]
    for dataset in datasets:
        for latent_dim in latent_dims:
            # input data, link rows i positions with i<j
            sparse_i=torch.from_numpy(np.loadtxt("./datasets/Original data/"+dataset+'/sparse_i.txt')).long().to(device)
            # input data, link column positions with i<j
            sparse_j=torch.from_numpy(np.loadtxt("./datasets/Original data/"+dataset+'/sparse_j.txt')).long().to(device)
            
            if args.LP:
                # file denoting rows i of missing links, with i<j 
                sparse_i_rem=torch.from_numpy(np.loadtxt("./datasets/Negative samples/"+dataset+'/sparse_i_rem.txt')).long().to(device)
                # file denoting columns j of missing links, with i<j
                sparse_j_rem=torch.from_numpy(np.loadtxt("./datasets/Negative samples/"+dataset+'/sparse_j_rem.txt')).long().to(device)
                # file denoting negative sample rows i, with i<j
                non_sparse_i=torch.from_numpy(np.loadtxt("./datasets/Negative samples/"+dataset+'/non_sparse_i.txt')).long().to(device)
                # file denoting negative sample columns, with i<j
                non_sparse_j=torch.from_numpy(np.loadtxt("./datasets/Negative samples/"+dataset+'/non_sparse_j.txt')).long().to(device)
               
            else:
                non_sparse_i=None
                non_sparse_j=None
                sparse_i_rem=None
                sparse_j_rem=None
                
            N=int(sparse_j.max()+1)
            #Missing data here denoted if Marginalization is applied or not
            # In case missing data is set to True then the input should be the complete graph
            sample_size=int(args.sample_percentage*N)
            model = LSM(sparse_i,sparse_j,N,latent_dim=latent_dim,sample_size=sample_size,non_sparse_i=non_sparse_i,non_sparse_j=non_sparse_j,sparse_i_rem=sparse_i_rem,sparse_j_rem=sparse_j_rem,CVflag=True,graph_type='undirected',missing_data=False,device=device,p=args.p).to(device)
            optimizer = optim.Adam(model.parameters(), args.lr)  
            if args.p==1:
                model.reg_l=args.delta_sq**0.5
            else:
                model.reg_l=args.delta_sq
            elements=(N*(N-1))*0.5
            for epoch in tqdm(range(args.epochs),desc="HM-LDM is Running…",ascii=False, ncols=75):
                if epoch==args.scaling_epochs:
                    model.scaling=0

                loss=-model.LSM_likelihood_bias(epoch=epoch)/sample_size

                optimizer.zero_grad() # clear the gradients.   
                loss.backward() # backpropagate
                optimizer.step() # update the weights
                if epoch%1000==0:
                      print('Iteration Number:', epoch)
                      print('Negative Log-Likelihood:',(loss.item()*N)/elements)
                      if args.LP:
                          roc,pr=model.link_prediction() 
                          print('AUC-ROC:',roc)
                          print('AUC-PR:',pr)
                          
    plt.rcParams["figure.figsize"] = (10,10)
    
    z_idx=model.latent_z.argmax(1)
    w_idx=model.latent_z.argmax(1)
    
    f_z=z_idx.argsort()
    f_w=w_idx.argsort()
    
    new_i=torch.cat((sparse_i,sparse_j))
    new_j=torch.cat((sparse_j,sparse_i))
    
    D=csr_matrix((np.ones(new_i.shape[0]),(new_i.cpu().numpy(),new_j.cpu().numpy())),shape=(N,N))#.todense()
 
    
    D = D[:, f_w.cpu().numpy()][f_z.cpu().numpy()]
    
    
    plt.spy(D,markersize=4)
    plt.xticks([])
    plt.yticks([])
    #plt.savefig(f'adjacency_{args.dataset}.pdf')
    plt.show()


# Check for identifiable solutions and node champions
latent_z = model.latent_z
identifiable_count = 0
epsilon = 1e-4
observed_corners = set()
possible_corners = {(0,) * i + (1,) + (0,) * (args.D - i - 1) for i in range(args.D)}
expected_corners = args.D

# Function to check if a tensor is binary
def is_binary(tensor):
    return torch.all(torch.isclose(tensor, torch.tensor(0.0), atol=epsilon) | torch.isclose(tensor, torch.tensor(1.0), atol=epsilon))

# Check each node's latent representation
for node_embedding in latent_z:
    if is_binary(node_embedding):
        identifiable_count += 1
        rounded_embedding = (torch.round(node_embedding / epsilon) * epsilon).int()
        #observed_corners.add(tuple(rounded_embedding.tolist()))
        if not torch.all(torch.eq(rounded_embedding, torch.zeros_like(rounded_embedding))):
            observed_corners.add(tuple(rounded_embedding.tolist()))

# Calculate the fraction of node champions and print
fraction_node_champions = identifiable_count / len(latent_z)
print(f"Fraction of node champions: {fraction_node_champions:.2%}")

# Check if there is at least one champion node for each possible corner
if observed_corners.issuperset(possible_corners):
    print(f"At least one champion node for each corner of the {args.D}-dimensional simplex is observed.")
else:
    print(f"Identifiability condition is not met. Observed corners: {observed_corners}, Expected corners: {possible_corners}")

# Check if there is at least one champion node for each corner
if len(observed_corners) >= expected_corners:
    print(f"All {expected_corners} corners of the simplex are occupied by at least one champion node.")
else:
    print(f"Identifiability condition is not met. Observed corners: {len(observed_corners)}, Expected corners: {expected_corners}")


z_idx = (z_idx).cpu().numpy()
f_z = (f_z).cpu().numpy()
df = pd.DataFrame({'Node': np.arange(1, N+1), 'z_idx': z_idx, 'f_z': f_z})
df.to_csv(f'../ESBM-master/Community from HM-LDM/community_{args.dataset}_D{args.D}_p{args.p}.csv', index=False)


















from HM_LDM import LSM
import argparse
import torch
import numpy as np
import torch.optim as optim
import sys
from tqdm import tqdm
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import networkx as nx

sys.path.append('./src/')


parser = argparse.ArgumentParser(
    description='Hybrid Membership-Latent Distance Model')

parser.add_argument('--epochs', type=int, default=10000, metavar='N',
                    help='number of epochs for training (default: 10K)')

parser.add_argument('--scaling_epochs', type=int, default=2000, metavar='N',
                    help='number of epochs for learning initial scale for the random effects (default: 2K)')

parser.add_argument('--delta_sq', type=float, default=1000, metavar='N',
                    help='delta^2 hyperparameter controlling the volume of the simplex')

parser.add_argument('--cuda', type=eval,
                    choices=[True, False],  default=True,
                    help='CUDA training')

parser.add_argument('--p', type=int,
                    choices=[1, 2],  default=2,
                    help='L2 norm power (default: 1)')


parser.add_argument('--LP', type=eval,
                    choices=[True, False], default=False,
                    help='performs link prediction')

parser.add_argument('--D', type=int, default=8, metavar='N',
                    help='dimensionality of the embeddings (default: 8)')

parser.add_argument('--lr', type=float, default=0.1, metavar='N',
                    help='learning rate for the ADAM optimizer, for large values of delta 0.01 is more stable (default: 0.1)')

parser.add_argument('--sample_percentage', type=float, default=1, metavar='N',
                    help='Sample size network percentage, it should be equal or less than 1 (default: 1)')


parser.add_argument('--dataset', type=str, default='SN',
                    help='dataset to apply HM-LDM on')

def is_binary(tensor, epsilon=1e-4):
    return torch.all(torch.isclose(tensor, torch.tensor(0.0), atol=epsilon) | torch.isclose(tensor, torch.tensor(1.0), atol=epsilon))


sys.argv = ['']
del sys
args, unknown = parser.parse_known_args()

D = [5, 10, 25]

device = torch.device("cpu")
#torch.set_default_tensor_type('torch.FloatTensor')
torch.set_default_dtype(torch.float32)
torch.set_default_device('cpu')

if __name__ == "__main__":
    latent_dims = D
    datasets = [args.dataset]
    for dataset in datasets:
        D_col = []
        pr_large = []
        auc_scores = []
        fraction_node_champions_mat = []

        for latent_dim in latent_dims:
            # input data, link rows i positions with i<j
            sparse_i = torch.from_numpy(np.loadtxt(
                "./datasets/Negative samples/"+dataset+'/sparse_i.txt')).long().to(device)
            # input data, link column positions with i<j
            sparse_j = torch.from_numpy(np.loadtxt(
                "./datasets/Negative samples/"+dataset+'/sparse_j.txt')).long().to(device)

            if args.LP:
                # input data, link rows i positions with i<j
                sparse_i = torch.from_numpy(np.loadtxt(
                "./datasets/Original data/"+dataset+'/sparse_i.txt')).long().to(device)
                # input data, link column positions with i<j
                sparse_j = torch.from_numpy(np.loadtxt(
                "./datasets/Original data/"+dataset+'/sparse_j.txt')).long().to(device)
                # file denoting rows i of missing links, with i<j
                sparse_i_rem = torch.from_numpy(np.loadtxt(
                    "./datasets/Negative samples/"+dataset+'/sparse_i_rem.txt')).long().to(device)
                # file denoting columns j of missing links, with i<j
                sparse_j_rem = torch.from_numpy(np.loadtxt(
                    "./datasets/Negative samples/"+dataset+'/sparse_j_rem.txt')).long().to(device)
                # file denoting negative sample rows i, with i<j
                non_sparse_i = torch.from_numpy(np.loadtxt(
                    "./datasets/Negative samples/"+dataset+'/non_sparse_i.txt')).long().to(device)
                # file denoting negative sample columns, with i<j
                non_sparse_j = torch.from_numpy(np.loadtxt(
                    "./datasets/Negative samples/"+dataset+'/non_sparse_j.txt')).long().to(device)

            else:
                non_sparse_i = None
                non_sparse_j = None
                sparse_i_rem = None
                sparse_j_rem = None

            delta_sq_mult = [0.01, 0.1, 1, 10, 50, 100]
            #delta_sq_mult = [args.delta_sq]
            roc_col = []
            pr_col = []
            fraction_node_champions_list = []

            for delta_sq in delta_sq_mult:
                N = int(sparse_j.max()+1)
                # Missing data here denoted if Marginalization is applied or not
                # In case missing data is set to True then the input should be the complete graph
                sample_size = int(args.sample_percentage*N)
                model = LSM(sparse_i, sparse_j, N, latent_dim=latent_dim, sample_size=sample_size, non_sparse_i=non_sparse_i, non_sparse_j=non_sparse_j,
                            sparse_i_rem=sparse_i_rem, sparse_j_rem=sparse_j_rem, CVflag=True, graph_type='undirected', missing_data=False, device=device, p=args.p).to(device)
                optimizer = optim.Adam(model.parameters(), args.lr)
                if args.p == 1:
                    model.reg_l = delta_sq**0.5
                else:
                    model.reg_l = delta_sq
                elements = (N*(N-1))*0.5

                # Lists to store results
                roc_loc = []
                like = []
                roc_curve_data = []  
                pr_curve_data = []

                for epoch in tqdm(range(args.epochs), desc="HM-LDM is Running…", ascii=False, ncols=75):
                    if epoch == args.scaling_epochs:
                        model.scaling = 0

                    loss = -model.LSM_likelihood_bias(epoch=epoch)/sample_size

                    optimizer.zero_grad()  # clear the gradients
                    loss.backward()  # backpropagate
                    optimizer.step()  # update the weights

                    if epoch % 1000 == 0:
                        print('Iteration Number:', epoch)
                        print('Negative Log-Likelihood:', (loss.item()*N)/elements)
                        if args.LP:
                            roc, pr, roc_curve = model.link_prediction()
                            print('AUC-ROC:', roc)
                            print('AUC-PR:', pr)
                            roc_curve_data.append(roc_curve)

                # Plot ROC curve for the current D
                #fpr, tpr, _ = roc_curve_data[-1]
                #plt.plot(fpr, tpr, label=f'AUC of D={latent_dim}: {roc:.2f}')

                # Store AUC-ROC score for later
                auc_scores.append((latent_dim, roc))

        
            # Add legend and labels to the plot
            #plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
            #plt.xlabel('False Positive Rate')
            #plt.ylabel('True Positive Rate')
            #plt.legend(loc='lower right')
            #plt.show()
            #plt.savefig(f'roc_{args.dataset}_p={args.p}_delta={args.delta_sq}.jpg')

                # Check fraction of node champions at each epoch
                identifiable_count = 0
                latent_z = model.latent_z
                for node_embedding in latent_z:
                    if is_binary(node_embedding):
                        identifiable_count += 1

                fraction_node_champions_list.append(identifiable_count / len(latent_z))
                if args.LP:
                    roc_col.append(roc)
                    pr_col.append(pr)

            fraction_node_champions_mat.append(fraction_node_champions_list)
            
            if args.LP:
                D_col.append(roc_col)
                pr_large.append(pr_col)


    # Plot Fraction of Node Champions vs. Delta_sq for different dimensions D
    plt.figure(figsize=(10, 6))
    default_x_ticks = range(len(delta_sq_mult))
    for i, y_values in enumerate(fraction_node_champions_mat):
        plt.plot(default_x_ticks, y_values, label=f'D={latent_dims[i]}', marker='o')

    plt.xticks(default_x_ticks, delta_sq_mult)
    plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.xlabel('$\delta^2$')
    plt.ylabel('Fraction of Node Champions')
    plt.legend()
    plt.grid(True)
    plt.show()
    


    plt.figure(figsize=(10,6))
    default_x_ticks = range(len(delta_sq_mult))
    for i, y_values in enumerate(D_col):
        plt.plot(default_x_ticks, y_values, label=f'D={latent_dims[i]}', marker='o')
    plt.xticks(default_x_ticks, delta_sq_mult)
    plt.xlabel('$\delta^2$')
    plt.ylabel('AUC-ROC')
    plt.legend()
    plt.grid(True)
    #plt.savefig(f'roc_{args.dataset}_{args.p}_new.jpg')
    plt.show()


    plt.rcParams["figure.figsize"] = (10, 10)

    labels = ['D=2', 'D=8', 'D=16', 'D=25']

    plt.figure(figsize=(10,6))
    default_x_ticks = range(len(delta_sq_mult))
    max_auc_scores = []
    for i, y_values in enumerate(D_col):
        max_auc = max(y_values)
        max_auc_scores.append((f'{labels[i]}: {max_auc:.2f}', max_auc))
        plt.plot(default_x_ticks, y_values, label=f'{labels[i]} (Max AUC: {max_auc:.2f})', marker='o')
    plt.xticks(default_x_ticks, delta_sq_mult)
    plt.xlabel('$\delta^2$')
    plt.ylabel('AUC-ROC')
    plt.legend()
    plt.grid(True)
    #plt.savefig(f'roc_{args.dataset}_{args.p}_new.jpg')
    plt.show()

    plt.rcParams["figure.figsize"] = (10, 10)

    z_idx = model.latent_z.argmax(1)
    w_idx = model.latent_z.argmax(1)

    f_z = z_idx.argsort()
    f_w = w_idx.argsort()

    new_i = torch.cat((sparse_i, sparse_j))
    new_j = torch.cat((sparse_j, sparse_i))

    D = csr_matrix((np.ones(new_i.shape[0]), (new_i.cpu().numpy(), new_j.cpu().numpy())), shape=(N, N))  # .todense()

    D = D[:, f_w.cpu().numpy()][f_z.cpu().numpy()]

    plt.spy(D, markersize=1)
    plt.xticks([])
    plt.yticks([])
    #plt.savefig(f'adjacency_{args.dataset}.jpg')
    plt.show()






 










from HM_LDM import LSM
import argparse
import torch
import numpy as np
import torch.optim as optim
import sys
from tqdm import tqdm
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

sys.path.append('./src/')

parser = argparse.ArgumentParser(description='Hybrid Membership-Latent Distance Model')

parser.add_argument('--epochs', type=int, default=10000, metavar='N',
                    help='number of epochs for training (default: 10K)')

parser.add_argument('--scaling_epochs', type=int, default=2000, metavar='N',
                    help='number of epochs for learning initial scale for the random effects (default: 2K)')

parser.add_argument('--delta_sq', type=float, default=0.1, metavar='N',
                    help='delta^2 hyperparameter controlling the volume of the simplex')

parser.add_argument('--cuda', type=eval, 
                      choices=[True, False],  default=True,
                    help='CUDA training')

parser.add_argument('--p', type=int, 
                      choices=[1, 2],  default=2,
                    help='L2 norm power (default: 1)')

parser.add_argument('--LP',type=eval, 
                      choices=[True, False], default=False,
                    help='performs link prediction')

parser.add_argument('--D', type=int, default=8, metavar='N',
                    help='dimensionality of the embeddings (default: 8)')

parser.add_argument('--lr', type=float, default=0.1, metavar='N',
                    help='learning rate for the ADAM optimizer, for large values of delta 0.01 is more stable (default: 0.1)')

parser.add_argument('--sample_percentage', type=float, default=1, metavar='N',
                    help='Sample size network percentage, it should be equal or less than 1 (default: 1)')



parser.add_argument('--dataset', type=str, default='facebook',
                    help='dataset to apply HM-LDM on')


sys.argv = ['']
del sys
args, unknown = parser.parse_known_args()

device = torch.device("cpu")
#torch.set_default_tensor_type('torch.FloatTensor')
torch.set_default_dtype(torch.float32)
torch.set_default_device('cpu')

if __name__ == "__main__":
    latent_dims=[args.D]
    datasets=[args.dataset]
    for dataset in datasets:
        for latent_dim in latent_dims:
            # input data, link rows i positions with i<j
            sparse_i=torch.from_numpy(np.loadtxt("./datasets/Original data/"+dataset+'/sparse_i.txt')).long().to(device)
            # input data, link column positions with i<j
            sparse_j=torch.from_numpy(np.loadtxt("./datasets/Original data/"+dataset+'/sparse_j.txt')).long().to(device)
            
            if args.LP:
                # file denoting rows i of missing links, with i<j 
                sparse_i_rem=torch.from_numpy(np.loadtxt("./datasets/"+dataset+'/sparse_i_rem.txt')).long().to(device)
                # file denoting columns j of missing links, with i<j
                sparse_j_rem=torch.from_numpy(np.loadtxt("./datasets/"+dataset+'/sparse_j_rem.txt')).long().to(device)
                # file denoting negative sample rows i, with i<j
                non_sparse_i=torch.from_numpy(np.loadtxt("./datasets/"+dataset+'/non_sparse_i.txt')).long().to(device)
                # file denoting negative sample columns, with i<j
                non_sparse_j=torch.from_numpy(np.loadtxt("./datasets/"+dataset+'/non_sparse_j.txt')).long().to(device)
               
            else:
                non_sparse_i=None
                non_sparse_j=None
                sparse_i_rem=None
                sparse_j_rem=None   
            
            N=int(sparse_j.max()+1)
            #Missing data here denoted if Marginalization is applied or not
            # In case missing data is set to True then the input should be the complete graph
            sample_size=int(args.sample_percentage*N)
            model = LSM(sparse_i,sparse_j,N,latent_dim=latent_dim,sample_size=sample_size,non_sparse_i=non_sparse_i,non_sparse_j=non_sparse_j,sparse_i_rem=sparse_i_rem,sparse_j_rem=sparse_j_rem,CVflag=True,graph_type='undirected',missing_data=False,device=device,p=args.p).to(device)
            optimizer = optim.Adam(model.parameters(), args.lr)  
            if args.p==1:
                model.reg_l=args.delta_sq**0.5
            else:
                model.reg_l=args.delta_sq
            elements=(N*(N-1))*0.5
            for epoch in tqdm(range(args.epochs),desc="HM-LDM is Running…",ascii=False, ncols=75):
                if epoch==args.scaling_epochs:
                    model.scaling=0
                

                loss=-model.LSM_likelihood_bias(epoch=epoch)/sample_size        
             
                optimizer.zero_grad() # clear the gradients.   
                loss.backward() # backpropagate
                optimizer.step() # update the weights
                if epoch%1000==0:
                      print('Iteration Number:', epoch)
                      print('Negative Log-Likelihood:',(loss.item()*N)/elements)
                      if args.LP:
                          roc,pr=model.link_prediction() 
                          print('AUC-ROC:',roc)
                          print('AUC-PR:',pr)

                          
    plt.rcParams["figure.figsize"] = (10,10)
    
    z_idx=model.latent_z.argmax(1)
    w_idx=model.latent_z.argmax(1)

    # Export membership vector
    community_memberships = (z_idx+1).cpu().numpy()
    df = pd.DataFrame({'Node': np.arange(1, N + 1), 'CommunityMembership': community_memberships})
    #df.to_csv(f'../ESBM-master/Community from HM-LDM/community_{args.dataset}_D{args.D}_p{args.p}_delta{args.delta_sq}.csv', index=False)
    
    f_z=z_idx.argsort()
    f_w=w_idx.argsort()
    
    new_i=torch.cat((sparse_i,sparse_j))
    new_j=torch.cat((sparse_j,sparse_i))
    
    D=csr_matrix((np.ones(new_i.shape[0]),(new_i.cpu().numpy(),new_j.cpu().numpy())),shape=(N,N))#.todense()
 
    D = D[:, f_w.cpu().numpy()][f_z.cpu().numpy()]
    
    plt.spy(D,markersize=2, color='black')
    plt.xticks([])
    plt.yticks([])
    plt.axis('off')

    #plt.savefig(f'adjacency_{args.dataset}_D{args.D}_p{args.p}_delta{args.delta_sq}.jpg')
    plt.show()




    # Create a graph using NetworkX
    G = nx.Graph()
    G.add_edges_from(zip(new_i.cpu().numpy(), new_j.cpu().numpy()))

    # Create a dictionary to map nodes to their corresponding community
    community_mapping = {node: z_idx[i] for i, node in enumerate(G.nodes)}

    # Get a list of colors based on community assignments
    node_colors = [community_mapping[node] for node in G.nodes]

    # Get unique community assignments
    unique_communities = np.unique(z_idx.cpu().numpy())
    num_communities = len(unique_communities)

    # Plot the graph with node colors
    plt.figure(figsize=(10, 10))
    pos = nx.spring_layout(G)  # You can use other layouts as well
    nx.draw(G, pos, node_color=node_colors, cmap=plt.cm.rainbow, with_labels=False, node_size=50)
    # Additional settings for the plot
    plt.xticks([])
    plt.yticks([])
    plt.axis('off')

    # Add a legend with the number of communities
    plt.suptitle(f'Number of Communities: {num_communities}', y=0.95, fontsize=16)

    # Save or show the plot
    plt.savefig(f'community_graph_{args.dataset}_D{args.D}_p{args.p}_delta{args.delta_sq}.jpg')
    plt.show()
