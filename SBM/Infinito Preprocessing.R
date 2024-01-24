rm(list=ls())

A <- read.csv(file="Application/NDRANGHETAMAFIA_2M.csv",header=TRUE, stringsAsFactors = TRUE)

# Manually insert missing data from reported attendance
A[23,20] <- 1
A[48,23] <- 1
A[78,35] <- 1
A[115,22] <- 1

# suspects who never attended a summit
sel_empty <- which(apply(as.matrix(A[,-1]),1,sum)==0)

# suspects not recognized during the investigation process
A[c(38,105,106,125,135),1]

# indicators of the two groups of suspects to be excluded
sel_empty <- c(c(sel_empty),c(38,105,106,125,135))

# remove these suspects from the dataset
A <- A[-sel_empty,]

# create a vector with the actors' names
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

# indicators of those suspects who are not known to be part of the organization La Lombardia 
sel_miss <- which(Locale=="MISS" | Locale=="OUT")

# remove such suspects from the dataset
Locale_temp <- Locale[-sel_miss]
Role_temp <- Role[-sel_miss]
actors_temp <- actors[-sel_miss]
A <- A[-sel_miss,-sel_miss]

SN <- A
write.csv(SN, file="Adjacency matrices/adj_SN.csv", row.names = FALSE, col.names = FALSE)
save(SN, file="Adjacency matrices/adj_SN.RData")

# Create the dataset used for modeling and inference
sel <- which(Locale_temp=="A" | Locale_temp=="B" | Locale_temp=="C" | Locale_temp=="D" | Locale_temp=="E")
Y <- A[sel,sel]
RoleLocale <- RoleLocale_temp[sel]
Role <- Role_temp[sel]
Locale <- Locale_temp[sel]
actors <- actors_temp[sel]
rownames(Y) <- colnames(Y) <- c(1:dim(Y)[1])

# Create the test dataset used for assessing predictive performance
Y_test <- A[-sel,sel]
rownames(Y_test) <- c(1:dim(Y_test)[1])
colnames(Y_test) <- c(1:dim(Y_test)[2])
RoleLocale_test <- RoleLocale_temp[-sel]
Role_test <- Role_temp[-sel]
Locale_test <- Locale_temp[-sel]
actors_test <- actors_temp[-sel]

save(Y,Y_test,Locale,Locale_test,Role,Role_test,RoleLocale,RoleLocale_test,file="Application/crime_net.RData")
save(Y, file="Adjacency matrices/adj_SN.csv")

write.csv(Y, file = "Adjacency matrices/adj_SN.csv", row.names = FALSE)

rm(list=ls())
