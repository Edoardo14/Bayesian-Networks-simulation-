# Load the libraries
library(bnlearn)
library(gRain)
library(igraph)
library(Rgraphviz)





# Load the data of the Asia Network
data("asia")

# Change the names of the columns to match my paper
names(asia) <- c("V", "S", "T", "L", "B", "TL", "X", "D")

# Create the Asia network
asia_dag <- model2network("[V][S][T|V][L|S][B|S][TL|T:L][X|TL][D|TL:B]")

# Fit the model to the data
asia_model <- bn.fit(asia_dag, data = asia)

# Plot the fitted model
graphviz.chart(asia_model, layout="dot", type="barprob")





# CLIQUE TREE ALGORITHM
cliquetree <- function(dag, model) {
  
  # 1. Input:
  #           dag must be a "bn" object
  #           model must be a "bn.fit" object
  
  # 2. Moralisation
  moral_graph <- moral(dag)
  
  # 3. Triangulation
  triang_graph <- triangulation(moral_graph)
  
  # 4. Compute the junction tree
  # 5. Calculate the potentials
  junction_tree <- junctiontree_potentials(model)
  print(junction_tree)
  
  return(junction_tree)
}





# gRain internally handles triangulation
# But this function explains the process of triangulation

# Function for triangulation
triangulation <- function(moral_graph) {
  
  # Triangulation
  moral_igraph <- as.igraph(moral_graph)
  triang_igraph <- triangulate(moral_igraph)
  
  # Convert back to bnlearn format
  adj_matrix <- as_adjacency_matrix(triang_igraph, sparse = FALSE)
  node_names <- nodes(moral_graph)
  triang_graph <- empty.graph(node_names)
  
  #  Add edges and form the triangulated graph
  for (i in 1:(length(node_names) - 1)) {
    
    for (j in (i + 1):length(node_names)) {
      
      if (adj_matrix[i, j] == 1) {
        
        triang_graph$arcs <- rbind(triang_graph$arcs, c(node_names[i], node_names[j]))
        triang_graph$arcs <- rbind(triang_graph$arcs, c(node_names[j], node_names[i]))
      }
    }
  }
  
  return(triang_graph)
}





# Function to compute the junction tree and potentials
junctiontree_potentials <- function(model) {
  
  # Convert bn.fit object to grain object
  grain_obj <- as.grain(model)
  
  # Compile the junction tree to prepare it for inference
  compiled_grain <- compile(grain_obj)
  
  return(compiled_grain)
}





# Function to perform probabilistic inference
inference <- function(junction_tree, targets, evidence) {
  
  # targets must be a vector of character objects
  # evidence must be a "list" object
  # NULL evidence results in computing the marginal probability
  
  # Set evidence
  if (!is.null(evidence) && length(evidence) > 0) {
    
    junction_tree <- setEvidence(junction_tree, nodes = names(evidence), states = unlist(evidence))
  }
  
  # 6. Belief propagation
  propagated_junction_tree <- propagate(junction_tree)
  
  # Perform inference
  query_result <- querygrain(propagated_junction_tree, nodes = targets, type = "joint")
  return(query_result)
}





# Running the Clique Tree Algorithm
junction_tree <- cliquetree(asia_dag, asia_model)

# Perform an inference query on the joint probability of nodes "T" and "L"
result_joint_T_L <- inference(junction_tree, c("T", "L"), NULL)
print(result_joint_T_L)

# Perform an inference query on the joint probability of nodes "T", "L", and "S"
result_joint_T_L_S <- inference(junction_tree, c("T", "L", "S"), NULL)
print(result_joint_T_L_S)

# Perform an inference query on the joint probability of nodes "TL" given evidence T = "yes" and L = "yes"
result_joint_TL_given_T_and_L <- inference(junction_tree, c("TL"), list(T = "yes", L = "yes"))
print(result_joint_TL_given_T_and_L)

# Perform an inference query on node "L" given evidence S = "yes"
result_L_given_S <- inference(junction_tree, "L", list(S = "yes"))
print(result_L_given_S)

# Perform an inference query on node "D" given evidence T = "yes" and B = "no"
result_D_given_T_and_B <- inference(junction_tree, "D", list(T = "yes", B = "no"))
print(result_D_given_T_and_B)

# Perform an inference query on node "X"
result_X_no_evidence <- inference(junction_tree, "X", NULL)
print(result_X_no_evidence)