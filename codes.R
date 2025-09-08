
####################################################################################
############################### Preparing the dataset ##############################
####################################################################################

data <- read.csv("~/data.csv")

data[data == -99] <- NA
data[data == -88] <- NA
data[data == -77] <- NA

data <- data %>% drop_na(CAPE0)
data <- data %>% drop_na(eloc_nos_ctxlhbankssts_1)

data$CAPE0_group <- ifelse(data$CAPE0 > median(data$CAPE0), 1, 0)
data <- data %>% select(-CAPE0)


####################################################################################
############################## Performing GEE analysis #############################
####################################################################################

store <- NULL
for(i in 3:(dim(data)[2]-1)){ 
  
  tryCatch({
    
    model1 <- gee(CAPE0_group ~ data[,i], 
                  data = data, 
                  id = famid, 
                  family = binomial,
                  corstr = "exchangeable", na.action = na.omit)
    summary(model1)
    
    robustPval <- 2 * pnorm(-abs(summary(model1)$coefficients[, 5]))
    result_table <- cbind(round(exp(coef(model1)),3), 
                          round(exp(summary(model1)$coefficients[,1] - qnorm(0.975)*summary(model1)$coefficients[,4]),3), 	
                          round(exp(summary(model1)$coefficients[,1] + qnorm(0.975)*summary(model1)$coefficients[,4]),3),
                          round(robustPval,3))
    rownames(result_table)
    print(names(data)[i])
    colnames(result_table) <- c("Odds_Ratio", "95%_LL", "95%_UL", "p-value")
    print(result_table)
    
    capture.output(cat("Variable:"), paste(as.data.table(names(data)[i])), cat("\n"), cat("Summary_statistics: \n"), print(summary(data[,i])), cat("\n"),
                   cat("Marginal_model: \n"), result_table, cat("\n\n"), append=TRUE,
                   file = "results_CAPE0.txt") },
    error = function(e) {
      capture.output(cat("Variable:", names(data)[i], "\n\n"), cat("Summary_statistics: \n"),print(summary(data[,i])), cat("\n"),cat("Analysis could not be performed for this variable.\n\n"),
                     append=TRUE,
                     file = "results_CAPE0.txt")
      return(NULL)
    })
  if(dim(result_table)[1]==2){
    if(as.data.frame(result_table)$'p-value'[2]<=0.05) store <- c(store, names(data)[i])
  } else if(dim(result_table)[1]==3){
    if(as.data.frame(result_table)$'p-value'[2]<=0.05 || as.data.frame(result_table)$'p-value'[3]<=0.05) store <- c(store, names(data)[i])
  }
  
  print(i)
}

length(store)
store

store_last <- c("CAPE0_group", "gender", "age", "WASI_PRI_Composite", store) 
length(store_last)

df <- subset(data, select = c(store_last))
summary(df)
df <- na.omit(df)
dim(df)

summary(as.factor(data$CAPE0_group))
summary(as.factor(df$CAPE0_group))


####################################################################################
################################ Correlation Filter ################################
####################################################################################

corr_mat <- correlation(df, method = "spearman", p_adjust = "none")
View(corr_mat[corr_mat$rho>0.7,])

data_last <- df %>%
  select(-corr_mat[corr_mat$rho>0.7,]$Parameter2)    

dim(data_last)


####################################################################################
############################# Machine learning analysis ############################
####################################################################################

########## Dummy variable coding

data_last <- dummy_cols(data_last, select_columns = "gender")
data_last <- subset(data_last, select=-c(gender, gender_female))

data_last <- dummy_cols(data_last, select_columns = "SES")
data_last <- subset(data_last, select=-c(SES, SES_high)) 

data_last <- dummy_cols(data_last, select_columns = "FIGS")
data_last <- subset(data_last, select=-c(FIGS, FIGS_no_psychiatric_disorder_history_in_the_first_family_degree))


########## Train-test sets split

set.seed(123)
index=sample.split(data_last$CAPE0_group, SplitRatio=0.70)
train=data_last[index,]
test=data_last[!index,]
summary(as.factor(train$CAPE0_group))
summary(as.factor(test$CAPE0_group))

train$CAPE0_group <- as.factor(train$CAPE0_group)
test$CAPE0_group <- as.factor(test$CAPE0_group)

train$CAPE0_group <- plyr::revalue(train$CAPE0_group, 
                                   c("1" = "High", 
                                     "0" = "Low"))

test$CAPE0_group <- plyr::revalue(test$CAPE0_group, 
                                  c("1" = "High", 
                                    "0" = "Low"))


########## Standardization

preprocess_values_train <- preProcess(train, method = c("center", "scale"))
train.st <- predict(preprocess_values_train, train)
test.st <- predict(preprocess_values_train, test)


########## LASSO feature selection

train.st_lasso <- subset(train.st, select=-c(age, WASI_PRI_Composite, gender_male))

dim(train.st) 
dim(train.st_lasso)

tuneGrid_lasso <- expand.grid(
  .alpha=1,
  .lambda=10^seq(-3, 3, length = 100))

set.seed(123)
model_lasso <- train(CAPE0_group~., data = train.st_lasso, method = "glmnet",
                     tuneGrid  = tuneGrid_lasso, family="binomial")

coef_lasso <- as.data.frame.matrix(coef(model_lasso$finalModel, model_lasso$bestTune$lambda))

var.imp_lasso <- varImp(model_lasso)
plot(var.imp_lasso)

rownames(coef_lasso)[order(coef_lasso$s1, decreasing=TRUE)]
df_coef_lasso <- as.data.frame(cbind(rownames(coef_lasso), coef_lasso$s1))
non_zero_coefs <- df_coef_lasso[df_coef_lasso$V2 != 0, ]

non_zero_coefs <- non_zero_coefs[-1,]
dim(non_zero_coefs)

train_en <- subset(train.st, select = c("CAPE0_group", non_zero_coefs$V1, "age", "WASI_PRI_Composite", "gender_male"))


########## Elastic net

set.seed(123)

model_en <- train(CAPE0_group ~., data = train_en, method = "glmnet", 
                  tuneLength=20, family="binomial")

coef_en <- as.data.frame.matrix(coef(model_en$finalModel, model_en$bestTune$lambda))

var.imp <- varImp(model_en)
plot(var.imp) 

pred_en <- predict(model_en, test.st[,-1])
perf_en <- confMat(pred_en, test.st$CAPE0_group, verbose = TRUE, positive = "High")
roc(response = as.numeric(test.st$CAPE0_group), predictor = as.numeric(pred_en), ci=TRUE)

