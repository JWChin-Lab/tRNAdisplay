#1. Read back sequence file

sequences <- read.table("4_sequence_counts.txt", header = TRUE)

#2. Define genetic decryption function - compare a list of codons with corresponding amino acids. Note - TAG (Amber) is translated as Glutamine (N)
#in E. Coli

geneticdecryption <- function(input=c("ATG", "GCA")){
  codons <- c("TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TAA","TAG","TGT","TGC",
              "TGA","TGG","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG","CAT","CAC","CAA","CAG",
              "CGT","CGC","CGA","CGG","ATT","ATC","ATA","ATG","ACT","ACC","ACA","ACG","AAT","AAC",
              "AAA","AAG","AGT","AGC","AGA","AGG","GTT","GTC","GTA","GTG","GCT","GCC","GCA","GCG",
              "GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG")
  
  aminoacids <- c("F","F","L","L","S","S","S","S","Y","Y","*","*","C","C","*","W","L","L","L","L","P",
                  "P","P","P","H","H","Q","Q","R","R","R","R","I","I","I","M","T","T","T","T","N","N",
                  "K","K","S","S","R","R","V","V","V","V","A","A","A","A","D","D","E","E","G","G","G","G","X")
  return(aminoacids[match(input, codons, nomatch = "65")]) # "X" is returned as warning (item 65 in aminoacids)
}


#Define sequence translation function
translate <- function(insert){
  insert <- substring(insert, seq(1,nchar(insert)-1,3), seq(3, nchar(insert),3)) #Cuts insert into codons
  insert <- insert[insert!=""] #Removes overhanging bases
  #insert <- insert[insert!="---"] #Removes codon deletions, if desired
  return(paste(geneticdecryption(insert), collapse="")) #applies geneticdecryption, and removes spaces
}

#3 Apply translation and scale by Count
sequences$"Translated" <- sapply(sequences$Member, FUN=translate, 
                                 USE.NAMES=FALSE, simplify = TRUE)

translated <- as.data.frame(sequences[,c(4,2)])
translated$Cumul.Sum <- cumsum(sequences$Count) 
translated_seqlogo <- rep(NA, sum(translated$Count))
translated_seqlogo[1:translated$Cumul.Sum[1]] <- translated[1,1]
for(i in 2:nrow(translated)){
  translated_seqlogo[(translated$Cumul.Sum[i-1]+1):(translated$Cumul.Sum[i])] <- translated$Translated[i]
}

#4. Tabulate amino acids 

aminoacids <- as.data.frame(table(translated_seqlogo))
aminoacids <- aminoacids[order(aminoacids[,2], decreasing = TRUE),]  
colnames(aminoacids) <- c("Member", "Count")
aminoacids$"Proportion(%)" <- aminoacids$Count/sum(aminoacids$Count)
diversity <- round(nrow(aminoacids)/sum(aminoacids$Count)*100,1)
outputname <- paste("5_aminoacids_",diversity, "%.txt",sep="")

write.table(aminoacids, outputname, row.names = FALSE, quote = FALSE)

  
  



