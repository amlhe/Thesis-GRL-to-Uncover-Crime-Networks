rm(list=ls())

# Read in data
A <- read.csv(file="Application/NDRANGHETAMAFIA_2M.csv",header=TRUE, stringsAsFactors = TRUE)

# Manually insert missing data from reported attendance
A[23,20] <- 1
A[48,23] <- 1
A[78,35] <- 1
A[115,22] <- 1

# Nodes who never attended a summit
sel_empty <- which(apply(as.matrix(A[,-1]),1,sum)==0)

# Nodes not recognized during the investigation process
A[c(38,105,106,125,135),1]

# remove these suspects from the dataset
A <- A[-sel_empty,]

# create a vector with names
actors <- A[,1]
actors <- droplevels(actors)

A <- as.matrix(A[,-1])
A <- A%*%t(A)
A <- (A>0)*1
diag(A) <- 0
rownames(A) <- colnames(A) <- c(1:dim(A)[1])

# ----------------------------------------------------
# Locale membership ("OUT": Suspects not belonging to La Lombardia. "MISS": Information not available)

Locale <- c("C","OUT","A","MISS","O","A","MISS","D","D","D","D","D","C","P","L","L","Q","MISS","B","OUT","B","B","I","MISS","OUT","D","A","O","N","N","H","OUT","D","E","G","G","L","A","OUT","Q","C","OUT","Q","L","C","MISS","C","C","F","C","OUT","D","A","B","B","E","M","MISS","C","C","C","B","H","C","C","E","E","E","E","C","MISS","L","A","A","E","E","C","E","E","E","C","MISS","OUT","C","C","E","G","A","A","B","I","I","A","B","B","OUT","I","A","G","N","E","D","F","OUT","OUT","C","D","C","MISS","MISS","C","MISS","E","E","C","MISS","OUT","B","L","A","D","D","O","MISS","B","D","O","D","D","A","A","I","C","MISS","MISS","MISS","A","A","F","E","C","Q","H","B","B","B")   

# ----------------------------------------------------
# Leadership role ("miss": Information not available)

Role <- c("aff","aff","aff","miss","aff","boss","miss","boss","boss","aff","aff","aff","aff","aff","aff","aff","aff","miss","aff","boss","boss","boss","boss","miss","boss","boss","aff","aff","aff","boss","aff","boss","aff","aff","aff","aff","aff","aff","boss","aff","aff","aff","boss","aff","aff","miss","aff","aff","aff","aff","boss","aff","aff","aff","aff","aff","boss","miss","aff","aff","aff","aff","boss","boss","boss","aff","boss","aff","aff","boss","miss","aff","aff","boss","boss","aff","aff","aff","aff","aff","aff","miss","boss","aff","aff","aff","aff","aff","aff","boss","aff","boss","aff","aff","aff","aff","boss","boss","boss","boss","aff","aff","aff","aff","aff","aff","aff","boss","miss","miss","aff","miss","aff","aff","aff","miss","aff","aff","boss","aff","aff","aff","aff","miss","aff","aff","boss","aff","boss","aff","aff","aff","aff","miss","miss","miss","aff","aff","boss","aff","aff","aff","boss","aff","aff","boss")

# Suspects known not to be part of the organization La Lombardia 
sel_miss <- which(Locale=="MISS" | Locale=="OUT")

# Clean dataset
Locale_temp <- Locale[-sel_miss]
Role_temp <- Role[-sel_miss]
actors_temp <- actors[-sel_miss]
A <- A[-sel_miss,-sel_miss]

SN <- A
write.csv(SN, file="Adjacency matrices/adj_SN.csv", row.names = FALSE, col.names = FALSE)
save(SN, file="Adjacency matrices/adj_SN.RData")
