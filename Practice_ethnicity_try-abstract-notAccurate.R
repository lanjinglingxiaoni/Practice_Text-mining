# init
libs <- c("tm", "plyr", "RISmed", "class", "SnowballC")
lapply(libs, require, character.only = TRUE)

# set parameter
eth <- c("African", "America", "Asian", "Caucasian", "MiddleEast", "notClear", "oceania")
file <- "C:/Users/tingjing/Desktop/Trainset_ethnicy.txt"

# read list and retrieve title+abstract from pubmed
pmid <- read.table(file, stringsAsFactors = TRUE, header = TRUE)

retrieve <-function(eth){
  pmid.eth <- pmid[pmid$Ethnicity == eth, "PMID"]
  ref <- EUtilsGet(pmid.eth)
  ref.df <- data.frame("name" = eth, "text" = paste0(ArticleTitle(ref),AbstractText(ref)))
}

results <- lapply(eth, retrieve)


# clean corpus
cleanCorpus <- function(corpus){
  corpus.tmp <- tm_map(corpus, removePunctuation)
  corpus.tmp <- tm_map(corpus.tmp, stripWhitespace)
  corpus.tmp <- tm_map(corpus.tmp, removeNumbers)
  corpus.tmp <- tm_map(corpus.tmp, tolower)
  corpus.tmp <- tm_map(corpus.tmp, removeWords, stopwords("english"))
  corpus.tmp <- tm_map(corpus.tmp, removeWords, c("one", "two", "three", "four", "five","six", 
                                                  "seven", "eight", "nine", "ten"))
  corpus.tmp <- tm_map(corpus.tmp, stemDocument) # remove common word endings
  return(corpus.tmp)
}


# build TDM
generateTDM <- function(results){
  res.cor <- Corpus(VectorSource(results[,"text"]))
  res.cl <- cleanCorpus(res.cor)
  res.tdm <- TermDocumentMatrix(res.cl)
  
  res.tdm <- removeSparseTerms(res.tdm, 0.7)
  output <- list(name = unique(results["name"]), tdm = res.tdm)
}

tdm <- lapply(results, generateTDM)

# attach ethnicity
bindEthnicity <- function(tdm){
  tdm.t <- t(data.matrix(tdm[["tdm"]]))
  tdm.df <- as.data.frame(tdm.t, stringsAsFactors = TRUE)
  tdm.df <- cbind(tdm.df, tdm[["name"]][[1]])
  colnames(tdm.df)[ncol(tdm.df)] <- "ethnicity"
  return(tdm.df)
}
refTDM <- lapply(tdm, bindEthnicity)

# stack
refTDM.stack <- do.call(rbind.fill, refTDM)
refTDM.stack[is.na(refTDM.stack)] <- 0

# hold out training + testing
train.idx <- sample(nrow(refTDM.stack), ceiling(nrow(refTDM.stack) * 0.7))
test.idx <- (1:nrow(refTDM.stack))[-train.idx]

# model
refTDM.stack.eth <- refTDM.stack[, "ethnicity"]
refTDM.stack.term <- refTDM.stack[, !colnames(refTDM.stack) %in% "ethnicity"]
knn.pred <- knn(refTDM.stack.term[train.idx,], refTDM.stack.term[test.idx,], refTDM.stack.eth[train.idx])

# accuracy
conf.mat <- table("prediction" = knn.pred, "actual" = refTDM.stack.eth[test.idx])
accuracy <- sum(diag(conf.mat))/length(test.idx)
