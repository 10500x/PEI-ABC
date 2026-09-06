
### Matriz de juguete de las peliculas ###

R <- matrix(c(4, NA, NA, 2, NA, 5,
              NA, 1, 5, NA, 2, NA,
              2, 5, 3, 4, NA, NA,
              NA, NA, NA, 5, 4, 3,
              NA, 2, NA, NA, 3, 5,
              1, NA, 2, 4, 5, NA), nrow=6, byrow=TRUE)

### Obtengo códigos países
ISO3_country <- open_dataset("df_paises.parquet") %>%
  select(iso_wits) %>% 
  collect() %>% 
  pull(iso_wits)


open_dataset("df_paises.parquet") %>%
  collect() %>% 
  View()

### Obtengo años para esta revisión
Year <- open_dataset("df_impo.parquet") %>%
  filter(revision == "HS22") %>% 
  select(Year) %>% 
  distinct(Year) %>% 
  collect() %>%
  pull(Year)

open_dataset("df_product.parquet") %>%
  collect() %>% 
  View()

open_dataset("df_product.parquet") %>%
  collect() %>% 
  View()

Cod_prod <- open_dataset("df_product.parquet") %>%
  filter(!is.na(desc_prod_22)) %>% 
  select(Cod_prod) %>% 
  collect() %>% 
  pull(Cod_prod)

### Genero todas las combinaciones 
df_merge <- expand_grid(Year, ISO3_country, Cod_prod)

df_merge <- df_merge %>% 
  left_join(open_dataset("df_impo.parquet") %>%
              filter(PartnerISO3 == "ARG" & revision == "HS22") %>% 
              select(ReporterISO3,
                     ProductCode,
                     Year,
                     TradeValue.in.1000.USD,
                     NetWeight.in.KGM) %>%
              rename(ISO3_country = ReporterISO3,
                     Cod_prod = ProductCode,
                     MA_tij = TradeValue.in.1000.USD,
                     QA_tij = NetWeight.in.KGM) %>% 
              collect(), 
            by = c("Year", "ISO3_country", "Cod_prod"),
            relationship = "many-to-one")

df_merge <- df_merge %>% 
  left_join(open_dataset("df_impo.parquet") %>%
              filter(PartnerISO3 == "WLD" & revision == "HS22") %>% 
              select(ReporterISO3,
                     ProductCode,
                     Year,
                     TradeValue.in.1000.USD,
                     NetWeight.in.KGM) %>%
              rename(ISO3_country = ReporterISO3,
                     Cod_prod = ProductCode,
                     M_tij = TradeValue.in.1000.USD,
                     Q_tij = NetWeight.in.KGM) %>% 
              collect(), 
            by = c("Year", "ISO3_country", "Cod_prod"),
            relationship = "many-to-one")

df_merge = df_merge %>% 
  mutate(w_tij=MA_tij/M_tij)

df_merge <- df_merge %>%
  rename(year = Year,
         ISO3 = ISO3_country,
         posHS17 = Cod_prod)

df_merge <- df_merge %>%
  filter(ISO3 != "ARG" & ISO3 != "WLD")

list_countries <- open_dataset("df_distance.parquet") %>%
  distinct(iso3_o) %>%
  rename(iso_of = iso3_o) %>% 
  filter(iso_of != "ARG") %>% 
  collect() %>% 
  left_join(open_dataset("df_paises.parquet") %>%
              select(iso_of, iso_wits) %>% 
              collect(),
            by = "iso_of") %>% 
  pull(iso_wits)

df_merge = df_merge %>%
  filter(ISO3 %in% list_countries)

### A partir de acá pueden armar la matriz R, 
### usen df_HS22_PEI.csv

R = df_HS22_PEI %>% 
  filter(year==2022) %>% 
  select(ISO3,posHS17,w_tij) %>% 
  pivot_wider(names_from=posHS17,
              values_from=w_tij)

R <- R %>%
  column_to_rownames("ISO3") %>%
  as.matrix()

#write.csv2(df_merge, "df_HS22_PEI.csv", row.names = FALSE)

#==============================================================#
# 1) Matrix-based Batch Gradient Descent (Aggarwal 2016 style) #
#==============================================================#

matrix_batch_gd <- function(R, k=2, gamma=0.01, max_iter=1000, tol=1e-4) {
  # Initialize
  set.seed(42)
  P <- matrix(rnorm(nrow(R)*k, 0, 0.1), nrow=nrow(R))
  Q <- matrix(rnorm(ncol(R)*k, 0, 0.1), nrow=ncol(R))
  
  # Observed indices
  observed <- which(!is.na(R), arr.ind=TRUE)
  n_obs <- nrow(observed)
  
  # Track RMSE for convergence
  rmse_history <- numeric()
  prev_rmse <- Inf
  
  for(iter in 1:max_iter) {
    # Compute prediction and error matrix
    R_hat <- P %*% t(Q)
    E <- R - R_hat
    
    # Set missing entries to zero (they don't contribute to gradient)
    E_zero <- E
    E_zero[is.na(E_zero)] <- 0
    
    # Update rules from Aggarwal (2016):
    # p_il <- p_il + γ * Σ_{j:(i,j)∈Ω} e_ij * q_jl
    # In matrix form: P ← P + γ * (E_zero) %*% Q
    
    P_new <- P + gamma * (E_zero %*% Q)
    Q_new <- Q + gamma * (t(E_zero) %*% P)
    
    # Update matrices
    P <- P_new
    Q <- Q_new
    
    # Compute RMSE on observed entries only
    predictions <- sapply(1:n_obs, function(idx) {
      i <- observed[idx, 1]
      j <- observed[idx, 2]
      sum(P[i, ] * Q[j, ])
    })
    current_rmse <- sqrt(mean((R[observed] - predictions)^2))
    rmse_history <- c(rmse_history, current_rmse)
    
    # Check convergence condition
    if(iter > 1) {
      rmse_change <- abs(prev_rmse - current_rmse)
      if(rmse_change < tol) {
        cat(paste("Converged at iteration", iter, "- RMSE change:", round(rmse_change, 6), "\n"))
        break
      }
    }
    prev_rmse <- current_rmse
    
    # Print progress
    if(iter %% 100 == 0) {
      cat(paste("Iteration", iter, "- RMSE:", round(current_rmse, 4), "\n"))
    }
  }
  
  return(list(P=P, Q=Q, rmse_history=rmse_history, iterations=iter))
}


result = matrix_batch_gd(R)

round(result$P%*%t(result$Q),2)

rownames(R)=c("Usuario 1","Usuario 2","Usuario 3",
              "Usuario 4","Usuario 5","Usuario 6")
colnames(R)=c("Tonto y re tonto", "El padrino", "La pistola desnuda",  
"Buenos  Muchachos", "Casino", "Zoolander")

R_hat = result$P%*%t(result$Q)

rownames(R_hat)=rownames(R)
colnames(R_hat)=colnames(R)

# Simple function to list ALL unseen movies
recommend_all_unseen <- function(user_name, R_hat, R_actual) {
  pred <- R_hat[user_name, ]
  seen <- !is.na(R_actual[user_name, ])
  unseen <- pred[!seen]
  unseen_sorted <- sort(unseen, decreasing = TRUE)
  return(data.frame(Movie = names(unseen_sorted), Rating = round(unseen_sorted, 2)))
}

# Usage
recommend_all_unseen("Usuario 1", R_hat, R)

#===========================================================#
#     2)  SGD Algorithm based on Aggarwal (2016)            #
#===========================================================#

sgd_aggarwal <- function(R, k=2, gamma=0.01, max_iter=1000, tol=1e-4) {
  
  # Initialize matrices
  set.seed(42)
  P <- matrix(rnorm(nrow(R)*k, 0, 0.1), nrow=nrow(R))
  Q <- matrix(rnorm(ncol(R)*k, 0, 0.1), nrow=ncol(R))
  
  # Get observed entries
  observed <- which(!is.na(R), arr.ind=TRUE)
  n_obs <- nrow(observed)
  
  # Track RMSE for convergence
  prev_rmse <- Inf
  
  # SGD loop
  for(iter in 1:max_iter) {
    # Shuffle observed entries
    shuffle <- sample(1:n_obs)
    
    # Loop through each observed entry in shuffled order
    for(idx in shuffle) {
      i <- observed[idx, 1]
      j <- observed[idx, 2]
      
      # Compute error
      r_ij <- R[i, j]
      r_hat <- sum(P[i, ] * Q[j, ])
      e <- r_ij - r_hat
      
      # Vectorized updates (from Aggarwal)
      P[i, ] <- P[i, ] + gamma * e * Q[j, ]
      Q[j, ] <- Q[j, ] + gamma * e * P[i, ]
    }
    
    # Check convergence after each full pass (epoch)
    # Compute RMSE on observed entries
    predictions <- sapply(1:n_obs, function(idx) {
      i <- observed[idx, 1]
      j <- observed[idx, 2]
      sum(P[i, ] * Q[j, ])
    })
    current_rmse <- sqrt(mean((R[observed] - predictions)^2))
    
    # Check convergence condition
    if(abs(prev_rmse - current_rmse) < tol) {
      cat(paste("Converged at iteration", iter, "- RMSE change:", round(abs(prev_rmse - current_rmse), 6), "\n"))
      break
    }
    
    prev_rmse <- current_rmse
    
    # Optional: Print progress
    if(iter %% 10 == 0) {
      cat(paste("Iteration", iter, "- RMSE:", round(current_rmse, 4), "\n"))
    }
  }
  
  return(list(P=P, Q=Q, R_hat=P %*% t(Q)))
}

# Run the algorithm
result_sgd <- sgd_aggarwal(R)

# Get predictions
R_hat_sgd <- result_sgd$R_hat
rownames(R_hat_sgd) <- rownames(R)
colnames(R_hat_sgd) <- colnames(R)

# Show results
round(R_hat_sgd, 2)

recommend_all_unseen("Usuario 2", R_hat_sgd, R)

#==============================================================#
# 3) Matrix-based Batch Gradient Descent with regularization   #
#==============================================================#

matrix_batch_gd_reg <- function(R, k=2, gamma=0.01, lambda=0.1, max_iter=1000, tol=1e-4) {
  # Initialize
  set.seed(42)
  P <- matrix(rnorm(nrow(R)*k, 0, 0.1), nrow=nrow(R))
  Q <- matrix(rnorm(ncol(R)*k, 0, 0.1), nrow=ncol(R))
  
  # Observed indices
  observed <- which(!is.na(R), arr.ind=TRUE)
  n_obs <- nrow(observed)
  
  # Track RMSE for convergence
  rmse_history <- numeric()
  prev_rmse <- Inf
  
  for(iter in 1:max_iter) {
    # Compute prediction and error matrix
    R_hat <- P %*% t(Q)
    E <- R - R_hat
    
    # Set missing entries to zero (they don't contribute to gradient)
    E_zero <- E
    E_zero[is.na(E_zero)] <- 0
    
    # Update rules with regularization:
    # P^+ = P(1 - γλ) + γ * E_zero * Q
    # Q^+ = Q(1 - γλ) + γ * t(E_zero) * P
    
    P_new <- P * (1 - gamma * lambda) + gamma * (E_zero %*% Q)
    Q_new <- Q * (1 - gamma * lambda) + gamma * (t(E_zero) %*% P)
    
    # Update matrices
    P <- P_new
    Q <- Q_new
    
    # Compute RMSE on observed entries only
    predictions <- sapply(1:n_obs, function(idx) {
      i <- observed[idx, 1]
      j <- observed[idx, 2]
      sum(P[i, ] * Q[j, ])
    })
    current_rmse <- sqrt(mean((R[observed] - predictions)^2))
    rmse_history <- c(rmse_history, current_rmse)
    
    # Check convergence condition
    if(iter > 1) {
      rmse_change <- abs(prev_rmse - current_rmse)
      if(rmse_change < tol) {
        cat(paste("Converged at iteration", iter, "- RMSE change:", round(rmse_change, 6), "\n"))
        break
      }
    }
    prev_rmse <- current_rmse
    
    # Print progress
    if(iter %% 100 == 0) {
      cat(paste("Iteration", iter, "- RMSE:", round(current_rmse, 4), "\n"))
    }
  }
  
  return(list(P=P, Q=Q, rmse_history=rmse_history, iterations=iter))
}

result <- matrix_batch_gd_reg(R, k=2, gamma=0.01, lambda=0.1, max_iter=1000, tol=1e-4)

R_hat_gd_reg = result$P%*%t(result$Q)

recommend_all_unseen("Usuario 6", R_hat_gd_reg, R)


#===========================================================#
#     4)  SGD Algorithm with regularization                 #
#===========================================================#

sgd_aggarwal_reg <- function(R, k=2, gamma=0.01 , lambda=0.1, max_iter=1000, tol=1e-4) {
  
  # Initialize matrices
  set.seed(42)
  P <- matrix(rnorm(nrow(R)*k, 0, 0.1), nrow=nrow(R))
  Q <- matrix(rnorm(ncol(R)*k, 0, 0.1), nrow=ncol(R))
  
  # Get observed entries
  observed <- which(!is.na(R), arr.ind=TRUE)
  n_obs <- nrow(observed)
  
  # Track RMSE for convergence
  prev_rmse <- Inf
  
  # SGD loop
  for(iter in 1:max_iter) {
    # Shuffle observed entries
    shuffle <- sample(1:n_obs)
    
    # Loop through each observed entry in shuffled order
    for(idx in shuffle) {
      i <- observed[idx, 1]
      j <- observed[idx, 2]
      
      # Compute error
      r_ij <- R[i, j]
      r_hat <- sum(P[i, ] * Q[j, ])
      e <- r_ij - r_hat
      
      # Vectorized updates (from Aggarwal)
      P[i, ] <- P[i, ] + gamma*(e*Q[j, ]-lambda*P[i, ])
      Q[j, ] <- Q[j, ] + gamma*(e*P[i, ]-lambda*Q[j, ])
    }
    
    # Check convergence after each full pass (epoch)
    # Compute RMSE on observed entries
    predictions <- sapply(1:n_obs, function(idx) {
      i <- observed[idx, 1]
      j <- observed[idx, 2]
      sum(P[i, ] * Q[j, ])
    })
    current_rmse <- sqrt(mean((R[observed] - predictions)^2))
    
    # Check convergence condition
    if(abs(prev_rmse - current_rmse) < tol) {
      cat(paste("Converged at iteration", iter, "- RMSE change:", round(abs(prev_rmse - current_rmse), 6), "\n"))
      break
    }
    
    prev_rmse <- current_rmse
    
    # Optional: Print progress
    if(iter %% 10 == 0) {
      cat(paste("Iteration", iter, "- RMSE:", round(current_rmse, 4), "\n"))
    }
  }
  
  return(list(P=P, Q=Q, R_hat=P %*% t(Q)))
}

result <- sgd_aggarwal_reg(R)

R_hat_sgd_reg = result$P%*%t(result$Q)
rownames(R_hat_sgd_reg) <- rownames(R)
colnames(R_hat_sgd_reg) <- colnames(R)
recommend_all_unseen("Usuario 6", R_hat_sgd_reg, R)

result = sgd_aggarwal_reg(R, k=16, gamma=0.01, lambda=0.1, 
                 max_iter=1000, tol=1e-4)

R_hat = result$P%*%t(result$Q)

rownames(R_hat)=rownames(R)
colnames(R_hat)=colnames(R)

recommend_all_unseen <- function(user_name, R_hat, R_actual) {
  pred <- R_hat[user_name, ]
  seen <- !is.na(R_actual[user_name, ])
  unseen <- pred[!seen]
  unseen_sorted <- sort(unseen, decreasing = TRUE)
  return(data.frame(Product = names(unseen_sorted), Share = unseen_sorted))
}

recommend_all_unseen("USA", R_hat, R) %>% 
  View()

