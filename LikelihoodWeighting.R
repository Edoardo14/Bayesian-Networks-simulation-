# Load the libraries
library(bnlearn)





# Load the data of the Asia Network
data("asia")

# Change the names of the columns to match the paper
names(asia) <- c("V", "S", "T", "L", "B", "TL", "X", "D")

# Create the Asia network
asia_dag <- model2network("[V][S][T|V][L|S][B|S][TL|T:L][X|TL][D|TL:B]")

# Fit the model to the data
asia_model <- bn.fit(asia_dag, data = asia)





# LIKELIHOOD WEIGHTING ALGORITHM
likelihood_weighting <- function(dag, model, evidence, n_samples = 1000) {
  
  # 1. Input:
  #           dag must be a "bn" object
  #           model must be a "bn.fit" object
  #           evidence must be a "list" object
  #           n_samples must be a "numeric" object
  
  # 2. Initialisation: I use a list and a data frame
  # The weights are already all set to 1
  weights <- rep(1, n_samples)
  samples <- data.frame(matrix(ncol = length(names(model)), nrow = n_samples))
  colnames(samples) <- names(model)
  
  # 3. Generating the samples
  for (i in 1:n_samples) {
    
    # Iterate the nodes in topological order
    for (node in node.ordering(dag)) {
      
      if (node %in% names(evidence)) {
        
        # If there is evidence, add it to the sample
        samples[i, node] <- evidence[[node]]
        
        # Check if the node has parents
        if (length(model[[node]]$parents) == 0) {
          
          # The node has no parents
          # Update the weight of the sample
          prob <- model[[node]]$prob[evidence[[node]]]
          weights[i] <- weights[i] * prob
        } else {
          
          # The node has parents
          # Sample from its conditional probabilities
          conditional_prob <- model[[node]]$prob
          
          # Get the values of the parents
          parent_values <- as.character(samples[i, model[[node]]$parents])
          
          # Transform to data frame object
          conditional_prob_df <- as.data.frame(conditional_prob)
          
          # Remove unnecessary parts of the data frame
          for (j in 1:length(parent_values)) {
            
            conditional_prob_df <- subset(conditional_prob_df, conditional_prob_df[j+1]==parent_values[j])
          }
          
          # Update the weight of the sample
          prob <- subset(conditional_prob_df, conditional_prob_df[node]==evidence[[node]])$Freq
          weights[i] <- weights[i] * prob
        }
      } else {
        
        # Sample the node based on the values of its parents
        if (length(model[[node]]$parents) == 0) {
          
          # The node has no parents
          # Sample from its marginal probabilities
          marginal_prob <- model[[node]]$prob
          
          # Transform to data frame object
          marginal_prob_df <- as.data.frame(marginal_prob)
          colnames(marginal_prob_df)[1] <- node
          
          # Get the probabilities to sample from
          probs <- marginal_prob_df[, "Freq", drop = FALSE]
          rownames(probs) <- c("no", "yes")
        } else {
          
          # The node has parents
          # Sample from the conditional probabilities
          conditional_prob <- model[[node]]$prob
          
          # Get the values of the parents
          parent_values <- as.character(samples[i, model[[node]]$parents])
          
          # Transform to data frame object
          conditional_prob_df <- as.data.frame(conditional_prob)
          
          # Remove unnecessary parts of the data frame
          for (j in 1:length(parent_values)) {
            
            conditional_prob_df <- subset(conditional_prob_df, conditional_prob_df[j+1]==parent_values[j])
          }
          
          # Get the probabilities to sample from
          probs <- conditional_prob_df[, "Freq", drop = FALSE]
          rownames(probs) <- c("no", "yes")
        }
        
        # Sample from the node's probability distribution
        sampled_value <- sample(rownames(probs), size = 1, prob = probs$Freq)
        samples[i, node] <- sampled_value
      }
    }
  }
  
  # 4. Returning the samples and the corresponding weights
  return(list(samples = samples, weights = weights))
}






# Function for probabilistic inference
inference <- function(weightedsamples, targets) {
  
  # targets must be a vector of character objects
  
  # Extracting samples and weights
  samples <- weightedsamples$samples
  weights <- weightedsamples$weights
  
  # Create an empty list to store results
  target_probs <- list()
  
  # Estimate the probability for each target
  for (target in targets) {
    
    # Get the unique values for the target variable
    target_values <- unique(samples[[target]])
    
    # Initialize a vector to store the probabilities
    probs <- numeric(length(target_values))
    
    for (j in seq_along(target_values)) {
      
      # Calculate the weighted probability of each target value
      probs[j] <- sum(weights[samples[[target]] == target_values[j]]) / sum(weights)
    }
    
    # Store the probabilities in the list
    target_probs[[target]] <- data.frame(value = target_values, prob = probs)
  }
  
  return(target_probs)
}


# Example 1
weightedsamples <- likelihood_weighting(asia_dag, asia_model, 
                                        evidence = list(T = "yes", B = "no"))
result <- inference(weightedsamples, targets = c("L", "D"))
print(result)

# Example 2
weightedsamples <- likelihood_weighting(asia_dag, asia_model, 
                                        evidence = list(TL = "no"))
result <- inference(weightedsamples, targets = c("T", "S"))
print(result)

# Example 3
weightedsamples <- likelihood_weighting(asia_dag, asia_model, 
                                        evidence = list(T = "yes", S = "yes"))
result <- inference(weightedsamples, targets = c("B", "D"))
print(result)
