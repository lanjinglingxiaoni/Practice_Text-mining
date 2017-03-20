# init
libs <-c("tm", "class", "plyr", "RISmed", "SnowballC")
lapply(libs, require, character.only = TRUE)

# set options
options(stringsAsFactors = FALSE)

# set category
category <- c("drugResponse", "riskPrediction")
  
# list of PMID
drugResponse <- c("22270724","21841502","9590747","18309959","25514114","15118073","15696205","15710947","16011858","16115929","24025416","20479403","26051236","16187797","26508876","26314834","18626007","17317677","26206867","18372921")
riskPrediction <- c("7493024","22798144","16949048","17100994","22970155","16683254","22923021","22762150","18566915","24799981","15366000","12362047","15926618","12658575","16807412","16736289","18759827","24078348","21344236","18572020")


# obtain the abstract from RISmed
ref <-EUtilsGet(c(drugResponse, riskPrediction))
ref.df <- data.frame("pmid"=PMID(ref), "title"=ArticleTitle(ref), 
                   "abstract"=AbstractText(ref), "text" = paste0(ArticleTitle(ref),AbstractText(ref)))

# clean corpus
cleanCorpus <- function(corpus){
  corpus.tmp <- tm_map(corpus, removePunctuation)
  corpus.tmp <- tm_map(corpus.tmp, stripWhitespace)
  corpus.tmp <- tm_map(corpus.tmp, removeNumbers)
  corpus.tmp <- tm_map(corpus.tmp, tolower)
  corpus.tmp <- tm_map(corpus.tmp, removeWords, stopwords("english"))
  corpus.tmp <- tm_map(corpus.tmp, stemDocument) # remove common word endings
  return(corpus.tmp)
}
  
# generate term matrix (TDM)
generateTDM <- function(cat){
  cat.value <- eval(parse(text = cat))
  ref.cor <- Corpus(VectorSource(ref.df[ref.df$pmid %in% cat.value, "text"]))
  ref.cl <- cleanCorpus(ref.cor)
  ref.tdm <- TermDocumentMatrix(ref.cl)
  
  #ref.tdm <- removeSparseTerms(ref.tdm, 0.7)
  result <- list(name = cat, tdm= ref.tdm)
}

tdm <- lapply(category, generateTDM)


# attach name
bindCategoryToTDM <- function(tdm){
  ref.mat <- t(data.matrix(tdm[["tdm"]]))
  ref.mat.df <- as.data.frame(ref.mat, stringsAsFactors = FALSE)
  
  ref.mat.df <- cbind(ref.mat.df, rep(tdm[["name"]], nrow(ref.mat.df)))
  colnames(ref.mat.df)[ncol(ref.mat.df)] <- "refCategory"
  return(ref.mat.df) #????????????????????????????????????????????????????????????????????? R??????????????????????????????????????????????????????????????????????????????????????????,??????????????????????????????????????????return()??????????????? ??????????????????????????????,????????????return(),??????????????????, ???return()????????????????????????????????????????????????
}
refTDM <- lapply(tdm, bindCategoryToTDM)

# stack
refTDM.stack <-do.call(rbind.fill, refTDM)
refTDM.stack[is.na(refTDM.stack)] <- 0

# hold-out (randomly generate training and testing dataset)
train.idx <- sample(nrow(refTDM.stack), ceiling(nrow(refTDM.stack) * 0.7))
test.idx <- (1:nrow(refTDM.stack))[-train.idx]

# black list of terms
black <- c("gefitinib", "erlotinib", "egfr", "cetuximab", "ten", "eighteen")

# model
tdm.cat <- refTDM.stack[,"refCategory"]
tdm.term <- refTDM.stack[, ! colnames(refTDM.stack) %in% c("refCategory", black)]

knn.pred <- knn(tdm.term[train.idx,], tdm.term[test.idx,],tdm.cat[train.idx])

# accuracy
conf.mat <- table("prediction" = knn.pred, Actual = tdm.cat[test.idx])
accuracy <- sum(diag(conf.mat))/length(test.idx)

# summarize frequency of term
freq <- colSums(tdm.term)

# generate histogram for top 20 terms
ord <- order(freq)
top20 <- data.frame(name = names(freq[tail(ord,20)]), freq = freq[tail(ord,20)])
library(ggplot2)
p <- ggplot(top20, aes(name,freq)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle=45, hjust =1))
p

# generate word clound
library(wordcloud)
set.seed(120)
wordcloud(names(freq), freq, min.freq = 15, scale = c(5,.1), colors = brewer.pal(6, "Dark2"))

# generate hierarchal cluster
library(cluster)
d <- dist(t(tdm.term), method="euclidian") # e.g d <- dist(t(tdm.term[1:20]), method="euclidian")
fit <- hclust(d = d, method="ward")
plot(fit, hang= -1)

# generate k-means clustering
library(fpc)
library(cluster)
d <- dist(t(tdm.term), method = "euclidian")
kfit <- kmeans(d, 2)
clusplot(as.matrix(d), kfit$cluster, color = TRUE, shade = TRUE, labels= 2, lines=0)
