library(testthat)
library(ggplot2)
library(overlapping)
library(factoextra)
  
#goes into each folder, reads in data
experiments <- list.dirs('.', recursive = FALSE)
for (e in experiments){
    setwd(e)
    table <- list.files(pattern = "*%.txt")
    temp <- read.table(table, header = TRUE)
    variable <- paste(gsub("./","",e), "aminoacids", sep = "_")
    assign(variable, temp)
    setwd('..')
  }
viewsamples <- as.character(sapply(experiments, function(x) gsub("./", "",x)))
viewsamples
#overview of all data
#data will be called according to the numbering in "viewsamples".
  
  
#run this function for each substrate selection (plus, minus, input replicates) to generate a master table containing all necessary data
#filtering and visualisation of the data will come at a later stage
# synthax is m1: vector with positive replicates, m2: vector with negative replicates, m3: vector with input replicates, 
#name: numeric identifier of substrate, w: placeholder value (recommended is a value close to 1 e.g. 0.9999)
master_analysis <- function(m1,m2,m3,name,w) {
    #read all POSITIVE TABLES
    for(i in 1:length(m1)) {
      temp <- viewsamples[m1[i]]
      assign(paste("plustable", i, sep=""), eval(parse(text=paste("`",temp, "_aminoacids","`", sep=""))))
    }

    #merge positive tables (USING AND operator)
    plus <- plustable1
    for(i in 1:(length(m1)-1)) {
      tabletomerge <- eval(parse(text=paste("plustable", i+1, sep="")))
      plus <- merge(tabletomerge, plus, by = "Member")
    }
    
    
    #read all NEGATIVE TABLES
    for(i in 1:length(m2)) {
      temp <- viewsamples[m2[i]]
      assign(paste("minustable", i, sep=""), eval(parse(text=paste("`",temp, "_aminoacids","`", sep=""))))
    }
       
    #merge all NEGATIVE TABLES WITH OR operator according to the previosuly merged plus table
    for(e in 1:(length(m2))) {
      table_x <- eval(parse(text=paste("minustable", e, sep="")))
      ###1 and 2
      present_in_both <- na.exclude(match(table_x$Member, plus$Member))
      #define exclusive table No 1
      exclusive_in <- plus[-present_in_both,c(1,2,3)]
      exclusive_in[,2] <- w
      exclusive_in[,3] <- (w*(table_x[1,3]/table_x[1,2]))
      colnames(exclusive_in) <- colnames(table_x)
      merged_full <-rbind(table_x, exclusive_in)
      assign(paste("minustable", e, sep = ""), merged_full)
    }
    minus <- minustable1
    for(i in 1:(length(m2)-1)) {
      tabletomerge <- eval(parse(text=paste("minustable", i+1, sep="")))
      minus <- merge(tabletomerge, minus, by = "Member")
    }

    
    
    #read all INPUT TABLES
    for(i in 1:length(m3)) {
      temp <- viewsamples[m3[i]]
      assign(paste("naivetable", i, sep=""), eval(parse(text=paste("`",temp, "_aminoacids","`", sep=""))))
    }
    
    
    #define inall (blackdots that are both in input and output in all replicates)
    in_plus_and_input <- naivetable1
    for(i in 1:(length(m3)-1)) {
      tabletomerge <- eval(parse(text=paste("naivetable", i+1, sep="")))
      in_plus_and_input <- merge(tabletomerge,  in_plus_and_input, by = "Member")
    }
    in_plus_and_input <- merge(plus, in_plus_and_input, by = "Member")
    in_plus_and_input <- in_plus_and_input[,c(1,2)]
    
    #merge all NAIVE TABLES WITH OR operator to previosuly created plus table
    for(e in 1:(length(m3))) {
      table_x <- eval(parse(text=paste("naivetable", e, sep="")))
      ###1 and 2
      present_in_both <- na.exclude(match(table_x$Member, plus$Member))
      #define exclusive table No 1
      exclusive_in <- plus[-present_in_both,c(1,2,3)]
      exclusive_in[,2] <- w
      exclusive_in[,3] <- (w*(table_x[1,3]/table_x[1,2]))
      colnames(exclusive_in) <- colnames(table_x)
      merged_full <-rbind(table_x, exclusive_in)
      assign(paste("naivetable", e, sep = ""), merged_full)
    }
    naive <- naivetable1
    for(i in 1:(length(m3)-1)) {
      tabletomerge <- eval(parse(text=paste("naivetable", i+1, sep="")))
      naive <- merge(tabletomerge, naive, by = "Member")
    }

    #merge everything
    merged <- merge(plus, minus, by = "Member")
    merged <- merge(merged, naive, by = "Member")
    
    #define coordinates to navigate table
    cord_minus <- length(m1)*2+1
    cord_naive <- cord_minus + length(m2)*2
    
    #calculate mean_fractions (for plus, negative and input conditions)
    merged$mean_plus <- 1
    merged$sd_plus <- 1
    merged$mean_minus <- 1
    merged$sd_minus <- 1
    merged$mean_naive <- 1
    merged$sd_naive <- 1
    for(i in 1:nrow(merged)) {
      plus_vector <- c()
      for(e in 1:length(m1)) {
        plus_vector <- append(plus_vector, merged[i,1+e*2])
      }
      minus_vector <- c()
      for(e in 1:length(m2)) {
        minus_vector <- append(minus_vector, merged[i,cord_minus+e*2])
      }
      naive_vector <- c()
      for(e in 1:length(m3)) {
        naive_vector <- append(naive_vector, merged[i,cord_naive+e*2])
      }
      #calculate means and SDs (relatives)
      merged[i,"mean_plus"] <- mean(plus_vector)
      merged[i,"sd_plus"] <- sd(plus_vector)/mean(plus_vector)
      merged[i,"mean_minus"] <- mean(minus_vector)
      merged[i,"sd_minus"] <- sd(minus_vector)/mean(minus_vector)
      merged[i,"mean_naive"] <- mean(naive_vector)
      merged[i,"sd_naive"] <- sd(naive_vector)/mean(naive_vector)
    }
    merged$enrichment_plus <- merged$mean_plus/merged$mean_naive
    merged$enrichment_minus <- merged$mean_minus/merged$mean_naive
    merged$selectivity <- merged$enrichment_plus/merged$enrichment_minus
    
    #order according to selectivity
    merged <- merged[order(merged$"selectivity", decreasing = TRUE),]

    
    #sums counts and sets minimums in different conditions. This is used for differentiation of data (classification into different bins) exploiting the fact that placeholder value can be differentiated from 1 which is minimal count
    merged$sum_plus <- 1
    merged$sum_minus <- 1
    merged$sum_naive <- 1
    merged$min_plus <- 1
    merged$min_minus <- 1
    merged$min_naive <- 1
    for(i in 1:nrow(merged)) {
      plus_vector <- c()
      for(e in 1:length(m1)) {
        plus_vector <- append(plus_vector, merged[i,e*2])
      }
      minus_vector <- c()
      for(e in 1:length(m2)) {
        minus_vector <- append(minus_vector, merged[i,cord_minus-1+e*2])
      }
      naive_vector <- c()
      for(e in 1:length(m3)) {
        naive_vector <- append(naive_vector, merged[i,cord_naive-1+e*2])
      }
      merged[i,"sum_plus"] <- sum(plus_vector)
      merged[i,"sum_minus"] <- sum(minus_vector)
      merged[i,"sum_naive"] <- sum(naive_vector)
      merged[i,"min_plus"] <- min(plus_vector)
      merged[i,"min_minus"] <- min(minus_vector)
      merged[i,"min_naive"] <- min(naive_vector)
    }
    
    #defines different groups of data that have different levels of coverage (depending on whether they do not contain any placeholder values (in_all), they contain one or more placeholders (on verge), or they only contain place holders.
    #This can be done both for the negative and the input samples, therefore creating a multidimensional Venn diagramm.)
    #The groups are differentiated using exploiting that placeholder values are different from 1, the minimal count:
    
    #group that has no placeholder values in the minus samples
    in_plus_and_minus <- merged[(merged[,"min_minus"]) > w,]
    
    #group that has no placeholder in the input samples
    in_plus_and_naive <- merged[(merged[,"min_naive"]) > w,]
    
    #group that has no placeholder values at all
    in_all <- in_plus_and_minus[(in_plus_and_minus[,"min_naive"]) > w,]
    
    #group that has one or more (but not only) placeholder values in the minus samples
    on_verge_minus <- merged[(merged[,"sum_minus"]) > length(m2)*w,]
    on_verge_minus <- on_verge_minus[(on_verge_minus[,"min_minus"]) == w,]
    
    #group that has one or more (but not only) placeholder values in the input samples
    on_verge_naive <- merged[(merged[,"sum_naive"]) > length(m3)*w,]
    on_verge_naive <- on_verge_naive[(on_verge_naive[,"min_naive"]) == w,]
    
    
    on_verge <- on_verge_minus[(on_verge_minus[,"sum_naive"]) > length(m3)*w,]
    on_verge <- on_verge[(on_verge[,"min_naive"]) == w,]
    
    #group that only has placeholder values in the minus samples
    none_minus <- merged[(merged[,"sum_minus"]) == length(m2)*w,]
    
    #group that only has placeholder values in the input samples
    none_naive <- merged[(merged[,"sum_naive"]) == length(m3)*w,]
    
    #group that only has placeholder values in both the minus and the input samples
    none <- none_minus[(none_minus[,"sum_naive"]) == length(m3)*w,]
    
    
    #define name and write tables for all group. Table that will be processed subsequently is the master table
    name_output <- paste(name, "master_table.csv", sep="_")
    write.csv(merged, name_output, quote = FALSE, row.names = FALSE)
    name_output <- paste(name, "table_in_plus_and_minus.csv", sep="_")
    write.csv(in_plus_and_minus, name_output, quote = FALSE, row.names = FALSE)
    name_output <- paste(name, "table_in_plus_and_naive.csv", sep="_")
    write.csv(in_plus_and_naive, name_output, quote = FALSE, row.names = FALSE)
    name_output <- paste(name, "table_in_all.csv", sep="_")
    write.csv(in_all, name_output, quote = FALSE, row.names = FALSE)
    name_output <- paste(name, "table_on_verge_minus.csv", sep="_")
    write.csv(on_verge_minus, name_output, quote = FALSE, row.names = FALSE)
    name_output <- paste(name, "table_on_verge_naive.csv", sep="_")
    write.csv(on_verge_naive, name_output, quote = FALSE, row.names = FALSE)
    name_output <- paste(name, "table_on_verge.csv", sep="_")
    write.csv(on_verge, name_output, quote = FALSE, row.names = FALSE)
    name_output <- paste(name, "table_none_minus.csv", sep="_")
    write.csv(none_minus, name_output, quote = FALSE, row.names = FALSE)
    name_output <- paste(name, "table_none_naive.csv", sep="_")
    write.csv(none_naive, name_output, quote = FALSE, row.names = FALSE)
    name_output <- paste(name, "table_none.csv", sep="_")
    write.csv(none, name_output, quote = FALSE, row.names = FALSE)
    }
    

#this analysis takes a master table (as well as the other subgroup tables) as input, and filters and plots it for hit discovery
#synthax is: name: nummerical identifier used as a prefix before, r_plus:enrichment cut off for the detail analysis, s: selectivity cut off for the detail analysis, e2: error cut off, q and p: coordinates for the spindle plot
detail_analysis <- function(name,r_plus,s,e2,q,p) {
  
  #read in the tables
  table_name <- paste(name, "master_table.csv", sep = "_")
  merged <- read.csv(table_name, header = TRUE)
  table_name <- paste(name, "table_in_plus_and_minus.csv", sep = "_")
  in_plus_and_minus <- read.csv(table_name, header = TRUE)
  table_name <- paste(name, "table_on_verge_minus.csv", sep = "_")
  on_verge_minus <- read.csv(table_name, header = TRUE)
  table_name <- paste(name, "table_none_minus.csv", sep = "_")
  none_minus <- read.csv(table_name, header = TRUE)
  
  #filter according to the defined values
  filtered <- merged[(merged$sd_plus) < e2,]
  filtered_in_plus_and_minus <- in_plus_and_minus[(in_plus_and_minus$sd_plus) < e2,]
  filtered_on_verge_minus <- on_verge_minus[(on_verge_minus$sd_plus) < e2,]
  filtered_none_minus <- none_minus[(none_minus$sd_plus) < e2,]
  
  #creates spindle plot
  
  #creates spindle plot
  plot(log(filtered$enrichment_plus), log(filtered$selectivity), cex=1, col="#17BFB3", pch=20, ylim=c(-q,q), xlim=c(-3,p))
  
  #verges (points that can be found at the margin of sequencing, that contain one ore more placeholder values but are still somehow covered)
  points(log(filtered_on_verge_minus$enrichment_plus), log(filtered_on_verge_minus$selectivity), cex=1, col="#17BFB3", pch=20, ylim=c(-q,q), xlim=c(-3,p))
  #IN ALLs (points that can be found in all replicates)
  points(log(filtered_in_plus_and_minus$enrichment_plus), log(filtered_in_plus_and_minus$selectivity), cex=1, col="#111111", pch=20, ylim=c(-q,q), xlim=c(-3,p))
  
  
  #NONES (points that can't found in none rplicate and only contain placeholder values)
  points(log(filtered_none_minus$enrichment_plus), log(filtered_none_minus$selectivity), cex=1, col="#C30000", pch=20, ylim=c(-q,q), xlim=c(-3,p))
  
  
  #cut offs define for the output
  abline(a = NULL, b = NULL, h = log(s))
  abline(a = NULL, b = NULL, v = log(r_plus))
  
  #define name and write table
  creme <- filtered[(filtered$enrichment_plus) > r_plus,]
  creme <- creme[(creme$selectivity) > s,]
  name_output <- paste(name, "creme.csv", sep="_")
  write.csv(creme, name_output, quote = FALSE, row.names = FALSE)
}
  