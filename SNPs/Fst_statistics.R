
##################### prepare input table for hierfstat, with samples names and genotype of the individuals -one locus per column

library(tidyr)
library(dplyr)

SNPs_locus <- read.delim("Home/snps_/statistics/SNPs_locus_table.txt")

locus_transpose = t(SNPs_locus)

locus_transpose2=as.data.frame(locus_transpose)

#make first row as headers
names(locus_transpose2) <-locus_transpose2[1,]
locus_transpose2 <- locus_transpose2[-1,]

#add column with samples names
vec <- c("A1","A4","A7","B1","B4","B7",
         "C3","C4","C7","VI2_1","VI_2_4","VI_2_7",
         "V_2_4","V_2_7","V2_1","X_2_1","X_2_4", "X_2_8") 

locus_transpose3 <- cbind(locus_transpose2, sample = vec)  

#remove row names as headers
rownames(locus_transpose3)<-NULL

#move last column with samples names to first
locus_transpose4 <- locus_transpose3 %>%
  select(sample, everything())

#add column with depths and one with colonies

#shallow is '5', mesophotic is '45'

vec2 <- c("5","5","5","5","5","5",
         "5","5","5","45","45","45",
         "45","45","45","45","45",
         "45") 

#vec3 <- c("A","A","A","B","B","B","C","C","C",
          #"VI","VI","VI","V","V","V","X","X","X") #substitued letters with numbers

vec3 <- c("1","1","1","2","2","2","3","3","3",
          "4","5","4","5","4","5","6","6","6")

locus_transpose5 <- cbind(locus_transpose4, depth = vec2) 
locus_transpose6 <- cbind(locus_transpose5, colony = vec3)

#move last columns with depth and colony to 2nd and 3rd places

moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", 
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}

#move colony column after sample column using the moveme function above
locus_transpose7=locus_transpose6[moveme(names(locus_transpose6), "colony after sample")]

#move depth column after colony column using the moveme function above
locus_transpose7=locus_transpose7[moveme(names(locus_transpose7), "depth after colony")]

locus_transpose7$sample <- NULL #don't need this column

##################### Fst with package hierfstat, Method: Weir & Cockerham, 1984

library(hierfstat)

locus_transpose8 <- lapply(locus_transpose7,as.numeric) #pairwise.WCfst only works with numerics
locus_transpose9=as.data.frame(locus_transpose8)

Fst_between_colonies = pairwise.WCfst(locus_transpose9[,-2],diploid=TRUE) #Fst between colonies

Fst_between_depths = pairwise.WCfst(locus_transpose9[,-1],diploid=TRUE) #Fst between depths
ncol(locus_transpose9)

#### Tests the significance of the effect of depth on genetic differentiation

#test.between(data, test.lev, rand.unit, nperm, ...)
#test.lev: A vector containing the units from which to construct the contingency tables
#rand.unit: A vector containing the assignment of each observation to the units to be permutted

#### example from package
#data(gtrunchier)
#attach(gtrunchier)
#test whether the locality level has a significant effect on genetic structuring
#test.between(gtrunchier[,-c(1,2)],test.lev=Locality,rand.unit=Patch)

effect_of_depth = test.between(locus_transpose9, test.lev=vec2, rand.unit=vec3, nperm=999)
effect_of_depth

effect_of_colony = test.within(locus_transpose9[,-c(1,2)],within=vec2,test.lev=vec3, nperm=999)
effect_of_colony

