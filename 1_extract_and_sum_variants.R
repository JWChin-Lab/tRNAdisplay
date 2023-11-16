#1. Read aligned sequences file

aligned <- read.table("2_aligned.txt", header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
aligned <- aligned[-(1:3),]


#2. Alignment score is given in column 5, and a CIGAR string indicating deletions is given 
#in column 6. For now - seems reasonable to take all 453M (453 matches, i.e. no deletions) and
#potentially trim the lowest alignment scores. For now, low alignments scores are left in
#and processed downstream. 

full_length <- aligned[grep("453M", aligned[,6]),]
quality_cut <- max(as.numeric(full_length[,5]))*2/3  
full_length <- full_length[(as.numeric(full_length[,5]) >= quality_cut),]
truncated <- aligned[-grep("453M", aligned[,6]),]
truncated <- truncated[(as.numeric(truncated[,5]) >= quality_cut),]

#2a Processing the truncated reads to include deletions after 300 bp 
CIGAR <- substring(truncated$V6, 1, 3)
late_deletions <- which(as.numeric(CIGAR) >= 300 )
truncated <- truncated[late_deletions,]
truncated<- truncated[which(nchar(truncated[,6]) <= 10),]

#2b Recombine full_length + truncated 
full_data <- rbind(full_length,truncated)


#3 Analyse sequences
sequences <- as.data.frame(full_data[,10]) 
colnames(sequences) <- "Full sequence"
sequences$"length" <- nchar(as.character(sequences[,1]))
sequences$"trim" <- sub(".*CGTGTTGACAAAAA", "", sequences$`Full sequence` )
sequences$"trim_length" <- nchar(sequences$trim)
sequences$"library1" <- substring(sequences$trim, 14, 49) #library positions
sequences$"library2" <- substring(sequences$trim, 146, 166) #library positions
sequences$"library3" <- substring(sequences$trim, 266, 272) #library positions
sequences$"library4" <- substring(sequences$trim, 311, 324) #library positions
sequences$"library5" <- substring(sequences$trim, 365, 385) #library positions

print(table(sequences$trim_length))
sequences <- sequences[which(sequences$trim_length < 440),] #removing sequences with mutated positions

sequences$"full library" <- paste(substring(sequences$library1, 4,6),
                                  substring(sequences$library1, 10,12),
                                  substring(sequences$library1, 19,21),
                                  substring(sequences$library1, 22,24),
                                  substring(sequences$library1, 31,33),

                                  
                                  
                                  substring(sequences$library2, 4,6),
                                  substring(sequences$library2, 10,12),
                                  substring(sequences$library2, 16,18), 
                                  
                                  substring(sequences$library3, 4,6),
                                  
                                  substring(sequences$library4, 4,6),
                                  substring(sequences$library4, 10,12),
                                  
                                  substring(sequences$library5, 4,6),
                                  substring(sequences$library5, 16,18),
                                  sep="")

#At this point - discard mutations in main sequence; assuming that they do not lead
#to useful diversity.


#4. Count library members (removes potential seq mistakes)
libraries <- sequences$`full library`
libraries <- as.data.frame(table(libraries))
libraries <- libraries[order(libraries[,2], decreasing = TRUE),]
colnames(libraries) <- c("Member", "Count")
libraries$"Proportion(%)" <- libraries$Count/sum(libraries$Count)


write.table(libraries, "4_sequence_counts.txt", row.names = FALSE, quote = FALSE)


