# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Evaluating sufficiency and necessity in model building with Coincidence     #
# Analysis                                                                    #
# ---------------------                                                       #
# Replication material                                                        #
# [anonymized]                                                                #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Install development version of cna R package, i.e. version 3.5.3.4
# Other versions of the package will not work because they do not have
# the new evaluation measures required for the replication.
install.packages("cna_3.5.3.4.tar.gz", repos = NULL, type="source")
# Restart R
rstudioapi::restartSession()
# Check cna package version
packageVersion("cna") # ==> must be version 3.5.3.4
# Set working directory again
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load required packages
library(parallel)
library(MASS)
library(cna)
library(frscore)
library(ggplot2)
library(reshape2)
library(tikzDevice)
library(dplyr)
library(gridExtra)
library(cowplot)
library(ggpubr)

# source auxiliary functions
source("Auxfuncs.R")
# Detect cores for parallel processing
numCores <- detectCores() 

###################################################################
# Replicate Tables 1b and 1c and their CNA analyses (see pages 5-9)
###################################################################
# Set the ground truth.
gt <- "A*C + a*c + D <-> B"
tab1b <- selectCases(gt,full.ct("A*C*D*B")) |> ct2df()
rownames(tab1b) <-  paste0("c",1:nrow(tab1b))
# Table 1b
tab1b
# Analyze Table 1b with cna, at con=1 and cov=1.
cna(tab1b)
# Table 1c
tab1c <- rbind(tab1b[-3,],data.frame(A=c(0,0),C=c(0,1),D=c(1,0),B=c(0,1), 
            row.names = c("c9","c10")))
# Analyze Table 1c with cna, at con=1 and cov=1.
cna(tab1c) # no model
# Analyze Table 1c with cna, at con=.7 and cov=.8.
cna(tab1c, con=.7, cov=.8)


###########################################
# Simulation experiments (see pages 20-27)
###########################################

# Set replication seed.
set.seed(6771) 
# Set the number of datasets to be generated.
n <-1000 # For less resource intensive experiments, set e.g. n <- 100.
# Sample the seeds for the experiment.
seeds <- sample(.Machine$integer.max, n)
# Detect the number of CPUs available.
numCores <- detectCores()

### Data generation, as described on pages 20-21. ###
#####################################################
# Generate GTs and corresponding datasets from ground truths and seeds
datasets <- vector("list", n)
for(i in 1:n){
  # cat(i,"\n")
  set.seed(seeds[i])
  n_fac <-  sample(5:7,1)
  n_size <- sample(((2^{n_fac-1})/2):200,1)
  GT <- generate_GTs(num_factors = n_fac, compl=list(sample(1:3,1),sample(3:5,1)))
  datasets[[i]] <- suppressMessages(generate_data(GT, type="cs", sample_size=n_size,
                                      num_factors=n_fac))
  attr(datasets[[i]][[1]],"GT") <- GT
}

# Adjust prevalence.
# Set variation sequence.
var_seq <- seq(0.1,0.9,0.1)
# loop over the variation sequence to adjust prevalence
prev_data <- vector("list", length(var_seq))
for (i in 1:length(var_seq)){
  cat(i, "\n")
  for(j in 1:length(datasets)){
  set.seed(seeds[j])
  prev_data[[i]][[j]] <-  introInflation.with.constNoise(datasets[[j]], 
                                      attr(datasets[[j]][[1]],"GT"),
                                      t = list("A"), inflate.ratio = var_seq[i], 
                                      constant = "fragmentation")
  attr(prev_data[[i]][[j]][[1]],"GT") <- attr(datasets[[j]][[1]],"GT")
   }
}
# Save prevalence manipulated data.
saveRDS(prev_data, file = "data_cs_1000.RData")


### Data analysis and robustness scoring, as described on pages 22-24. ###
##########################################################################
# Source definitions of tested evaluation measures.
source("score_defs.R")
# Define analysis function.
analysis <- function(x,  ccDef){
 out <-  quiet(rean_cna(x, attempt = seq(0.95, 0.65, -0.1), maxstep = c(3, 5, 16), 
               outcome= "A", ccDef = ccDef))
 attr(out,"GT") <- attr(x,"GT")
 return(out)
 }
# Analysis loop.
results <- vector("list", length(score_defs))
for(k in 1:length(score_defs)){
  cat(k, "score \n")
  rean_series <- vector("list",length(prev_data))
  for(i in 1:length(prev_data)){
    cat(i, "\n")
    rean_series[[i]] <- mclapply(prev_data[[i]], function(x) analysis(x[[1]],   
    ccDef=score_defs[[k]]), mc.cores=numCores)   
  }
results[[k]] <- rean_series  
}
# Save results.
saveRDS(results, file = "results.RData")

# Robustness scoring
# Define coring function.
frscoring <- function(x){
  xx <- do.call(rbind,x)
  if(nrow(xx)>0){
           fr.scores <- quiet(frscore(xx$condition, maxsols=150))
          }else{
            fr.scores <- xx
          }
  attr(fr.scores,"GT") <- attr(x, "GT")
  return(fr.scores)
}

# Robustness scoring loop.
frscores.per.definition <- vector("list", length(results))
for(i in 1:length(results)){
  cat(i,"score \n")
  aux <- vector("list", length(results[[i]]))
  for(j in 1:length(results[[i]])){
    cat(j,"prev \n")
      aux[[j]]  <- mclapply(results[[i]][[j]],function(x)frscoring(x),
                            mc.cores = numCores)
  }
  frscores.per.definition[[i]] <- aux
}
# Save frscores.per.definition
saveRDS(frscores.per.definition, "frscores.per.definition.RData")

# Select the top 95the percentile.
results.select <- vector("list", length(frscores.per.definition))
for(i in 1:length(frscores.per.definition)){
  aux1 <- vector("list", length(frscores.per.definition[[i]]))
  for(j in 1:length(frscores.per.definition[[i]])){
    aux2 <- vector("list", length(frscores.per.definition[[i]][[j]]))
     for(k in 1:length(frscores.per.definition[[i]][[j]])){
       x <- frscores.per.definition[[i]][[j]][[k]]$models
       if(!is.null(x)){
       # select the 95 percentile on the norm.score column in x
        x.95 <- x[which(x$norm.score>=quantile(x$norm.score, 0.95)),]
         aux2[[k]] <-  list(quant95=x.95)
         attr(aux2[[k]], "GT") <- attr(frscores.per.definition[[i]][[j]][[k]], "GT")
       }else{
         x.95 <- frscores.per.definition[[i]][[j]][[k]]
         x.95$model <- character(0)
         aux2[[k]] <-  list(quant95=x.95)
         attr(aux2[[k]], "GT") <- attr(frscores.per.definition[[i]][[j]][[k]], "GT")
       }
     }
     aux1[[j]] <- aux2
  }
  results.select[[i]] <- aux1
  
}

### Benchmark scoring ###
#########################
results.scored95 <- vector("list", length(results.select))
for(i in 1:length(results.select)){
  aux95 <- vector("list", length(results.select[[i]]))
  for(j in 1:length(results.select[[i]])){
    out.scored.95 <- list()
       out.scored.95$correct <- lapply(results.select[[i]][[j]], function(x){
                                         corr_calc_frscore(x$quant95, attr(x, "GT"))})
      out.scored.95$complete <- lapply(results.select[[i]][[j]], function(x){ # nolint
                                         complete_calc_frscore(x$quant95, attr(x, "GT"))})
      out.scored.95$empty <- lapply(results.select[[i]][[j]], function(x)
                                                     ifelse(length(x$quant95) == 0, 1,
                                                        as.numeric(nrow(x$quant95) == 0)))
      out.scored.95$corr_prop <- lapply(results.select[[i]][[j]], function(x)
                                           corr_prop_calc_frscore(x$quant95, attr(x, "GT")))
      yy <- lapply(results[[i]][[j]], function(x)do.call(rbind,x))
      yy <- lapply(yy, function(x) x[which(!duplicated(x$condition)),])
      out.scored.95$ambiguity_select <- lapply(results.select[[i]][[j]], 
                                    function(x) amb_calc_frscore(x$quant95) )
    aux95[[j]] <- out.scored.95
  }
  results.scored95[[i]] <- aux95
}

# Fill out the scoreboards.
scoreboards.per.defscore <- vector("list", length(results.scored95))
for(i in 1:length(results.scored95)){
   scoreboard <- expand.grid(95,seq(0.1, 0.9, 0.1))
   colnames(scoreboard) <- c("quantile","prevalence")
   scoreboard <- cbind(scoreboard, correct=NA,complete=NA,empty=NA,corr_prop=NA,ambiguity_select=NA)
     for(j in 1:length(results.scored95[[i]])){
     scoreboard[j,3:7] <-   unlist(lapply(results.scored95[[i]][[j]], 
                            function(x) if(all(is.na(unlist(x)))){NA}else{mean(unlist(x), na.rm = T)}))
   }
scoreboards.per.defscore[[i]] <- scoreboard
}

# Overall aggregation.
overall <- vector("list", length(scoreboards.per.defscore))
for(i in 1:length(scoreboards.per.defscore)){
  aux <- vector("list", nrow(scoreboards.per.defscore[[i]]))
  for(j in 1:length(aux)){
    x <- scoreboards.per.defscore[[i]][j,]
    aux[[j]] <- (x$correct * x$complete * (1 - x$empty)^(1/2))
  }
 scoreboards.per.defscore[[i]]$overall <-  unlist(aux)  
 
split.by.setting <- split(scoreboards.per.defscore[[i]], list(scoreboards.per.defscore[[i]]$quantile))
 x <- unlist(lapply(split.by.setting, function(x)mean(replace(x$overall, is.na(x$overall),0))))      
 xx <- split.by.setting[which(x==max(x))]
 xx$mean.overall <- max(x)
 overall[[i]] <- xx
 }
overall

# Save scoreboards.
saveRDS(scoreboards.per.defscore, "scoreboards.per.defscore.RData")


### Plots ###
#############
n <-  1:5  
title <- c(paste("consistency / coverage"), 
           paste("PA-consistency / AA-coverage"),
           paste("AAC-consistency / AA-coverage"),
           paste("PA-consistency / PAC-coverage"),
           paste("AAC-consistency / PAC-coverage")
           )
plotting <- data.frame(n,title)

plot <- vector("list", length(plotting))
for(i in 1:nrow(plotting)){
  n <- plotting$n[i]
  title <- plotting$title[i]
scoresheet <- scoreboards.per.defscore[[n]]

selection1 <- as.data.frame(scoresheet)
selection1 <- selection1 %>% group_by(quantile, prevalence, correct, empty, complete, overall, 
                                      corr_prop, ambiguity_select) 
k1 <- melt(selection1, id.vars = c( "quantile","prevalence")) 
k11 <- k1[46:54,]
k11 <- rbind(k11,k11,k11)
colnames(k11) <- c("quantile", "prevalence", "variable1","value1")
k12 <- k1[1:27,]
colnames(k12) <- c("quantile", "prevalence", "variable2","value2")
k13 <- k1[37:45,]
k13 <- rbind(k13, k13,k13)
colnames(k13) <- c("quantile", "prevalence", "variable3","value3")
k1 <- cbind(k11, k12[,3:4],k13[,3:4])
k1$variable1 <- factor(k1$variable1,levels = c("overall"))
k1$variable2 <- factor(k1$variable2,levels = c("empty", "correct", "complete"))
k1$variable3 <- factor(k1$variable3,levels = c("ambiguity_select"))#, "ambiguity_total"))
k1$prevalence <- factor(k1$prevalence,levels = seq(0.1,0.9,0.1))
k1$quantile<- factor(k1$quantile,levels = selection1$quantile[1])

plot[[i]] <- ggplot(k1,aes(x = prevalence),na.rm=T) +
  geom_bar(aes(y = value2,fill = variable2), stat="identity",position = position_dodge(width=0.8), width=0.7) +
   scale_fill_manual(values=c('#aec9aa', '#221f1f','#df9199')) + facet_wrap(vars(quantile)) +
  theme_bw() + ggtitle(title) +
  geom_line(aes(y = value1,group=variable1,color=variable1,linetype=variable1),lwd=1, stat="identity") +
   geom_line(aes(y = value3/5,group=variable3,color=variable3,linetype=variable3), lwd=0.8, stat="identity")+
   scale_color_manual(values=c('#5233ff',"#fcdd3f")) + #,'#C70039','#0066ff'))+  # '#C70039',
      scale_linetype_manual(values = c("solid","solid","solid","dashed"))+
 theme(plot.title = element_text(size = 10), axis.text.x = element_text( size = 7 ),
        axis.text.y = element_text( size = 7 ),axis.title = element_text( size = 7),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 5, l = 5)),
        legend.spacing.x = unit(0.35, 'cm'),
        legend.title = element_blank(),strip.text = element_text(size = 7),
       legend.text=element_text(size=7,margin = margin(l = -9, r=0,unit = "pt")))+
   scale_x_discrete(name ="prevalence")+ theme(legend.position="right") +
  theme(legend.justification = "top")+
   scale_y_continuous(name = "score / ratio",limits = c(0,1),
    sec.axis = sec_axis(~.*5, name = "ambiguity", 
      labels = function(b) { round(b , 1)})
      )

}
ns <- t(data.frame(A =round(overall[[1]]$mean.overall,2),B=round(overall[[2]]$mean.overall,2), C=round(overall[[3]]$mean.overall,2),
                   D=round(overall[[4]]$mean.overall,2), E=round(overall[[5]]$mean.overall,2)))
ns <- cbind(rownames(ns),ns)
colnames(ns) <- c("Experiment","overall")
tt1 <- ttheme_default(base_size =8)
ngrob <- draw_grob(tableGrob(ns, rows = NULL,theme=tt1), x=0.72, y=0.25, width=0.001, height=0.001,scale = 0.2)

figure <- ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], plot[[5]],
                    labels = c("A", "B", "C", "D", "E"),
                    ncol = 2, nrow = 3, common.legend = TRUE, legend = "bottom") + ngrob 
figure
