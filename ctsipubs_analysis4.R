#load the following libraries for this analysis
library(bibliometrix)
library(tidyr) #tidy messy data
library(tidytext) 
library(dplyr)
library(ggplot2)
library(topicmodels)
library(textmineR)
library(stringr)
library(tm)
library(topicdoc)

setwd("~/R")
file <- "wos-ctsi-pubs.txt" #search results exported from Web of Science
#convert the exported data file from WoS to a dataframe
#convert2df creates a bibliographic data frame with cases corresponding
#to manuscripts and variables to Field Tag in the original export file.
ctsipubs = convert2df(file, dbsource="wos", format="plaintext")
names(ctsipubs) # displays columns of the data frame 

#split documents into train and test st using 80/20 rule
smp_siz = floor(0.8*nrow(ctsipubs))
set.seed(123)
train_ind = sample(seq_len(nrow(ctsipubs)),size = smp_siz)
train =ctsipubs[train_ind,] 
test=ctsipubs[-train_ind,]  

#descriptive analysis of "ctsipubs" bibliographic dataframe
results = biblioAnalysis(ctsipubs)
options(width=100)
report = summary(object = results, k = 10, pause = FALSE)

#Data cleaning and pre-processing
#preprocessing the data to remove delimiters from keywords
#tokenization

#Step 1. Create tidy data frame for entire abstract dataset, train and test sets
ctsi_abst <- tibble(id = ctsipubs$UT, 
                   abstract = ctsipubs$AB)
ctsi_abst
ctsi_abst %>% 
  na.omit()
train_abst <- tibble(id = train$UT, 
                    abstract = train$AB)
train_abst
train_abst %>% 
  na.omit()
test_abst <- tibble(id = test$UT, 
                    abstract = test$AB)
test_abst
test_abst %>% 
  na.omit()

#Step 2. Create tidy data frame for author keywords, train and test sets 
ctsi_keyword <- tibble(id = ctsipubs$UT, 
                       keyword = ctsipubs$DE)
ctsi_keyword
ctsi_keyword %>% 
  na.omit()

train_keyword <- tibble(id = train$UT, 
                       keyword = train$DE)
train_keyword
train_keyword %>% 
  na.omit()

test_keyword <- tibble(id = test$UT, 
                       keyword = test$DE)
test_keyword
test_keyword %>% 
  na.omit()

#Step 3. Tokenization to create one row per term or keyword for full dataset, train and test sets
#remove stopwords
ctsi_abst <- ctsi_abst %>%
  unnest_tokens(word, abstract) %>%
  anti_join(stop_words)
ctsi_abst
train_abst <- train_abst %>%
  unnest_tokens(word, abstract) %>%
  anti_join(stop_words)
train_abst
test_abst <- test_abst %>%
  unnest_tokens(word, abstract) %>%
  anti_join(stop_words)
test_abst

ctsi_keyword <- ctsi_keyword %>%
  unnest_tokens("keyword", "keyword", token = "regex", pattern = ";")
ctsi_keyword
train_keyword <- train_keyword %>%
  unnest_tokens("keyword", "keyword", token = "regex", pattern = ";")
train_keyword
test_keyword <- test_keyword %>%
  unnest_tokens("keyword", "keyword", token = "regex", pattern = ";")
test_keyword

#check most common words in abstract and keyword
ctsi_abst %>% 
  count(word, sort = TRUE) 
ctsi_abst_withcount <- ctsi_abst %>% 
  count(word, sort = TRUE)  
#write.csv(ctsi_abst_withcount, "ctsi_abst_withcount.csv")
ctsi_abst_withcount
ctsi_keyword %>% 
  count(keyword, sort = TRUE) %>% 
  na.omit()

#remove custom stop words from abstract 
my_stopwords <- read.csv("words.csv", header = FALSE)
my_stopwords <- as.character(my_stopwords$V1)

custom_stopwords <- tibble(word = c(as.character(1:10), 
                                    my_stopwords))
ctsi_abst <- ctsi_abst %>% 
  anti_join(custom_stopwords) 
ctsi_abst %>% 
  count(word, sort = TRUE) 
train_abst <- train_abst %>% 
  anti_join(custom_stopwords) 
train_abst %>% 
  count(word, sort = TRUE) 
test_abst <- test_abst %>% 
  anti_join(custom_stopwords) 
test_abst %>% 
  count(word, sort = TRUE) 

#convert keywords to lowercase 
ctsi_keyword <- ctsi_keyword %>% 
  mutate(keyword = tolower(keyword))
train_keyword <- train_keyword %>% 
  mutate(keyword = tolower(keyword))
test_keyword <- test_keyword %>% 
  mutate(keyword = tolower(keyword))

############################
##TOPIC MODELING USING LDA##
############################

#create a Document Term Matrix or DTM for full dataset, train and test sets
ctsi_word_counts <- ctsi_abst %>%
  count(id, word, sort = TRUE) %>%
  ungroup()
ctsi_word_counts

dtm_abst <- ctsi_word_counts %>% 
  cast_dtm(id, word, n)
dtm_abst
inspect(dtm_abst)

train_word_counts <- train_abst %>%
  count(id, word, sort = TRUE) %>%
  ungroup()
train_word_counts

dtm_train_abst <- train_word_counts %>% 
  cast_dtm(id, word, n)
dtm_train_abst

test_word_counts <- test_abst %>%
  count(id, word, sort = TRUE) %>%
  ungroup()
test_word_counts

dtm_test_abst <- test_word_counts %>% 
  cast_dtm(id, word, n)
dtm_test_abst

#LDA model with 10-100 topics increasing by 10
mod_perplexity = c()
for(i in seq(from = 5, to = 50, by = 5)) {
mod_lda <- LDA(dtm_train_abst, k = i, control = list(seed = 1234))
mod_perplexity[i] = perplexity(mod_lda, dtm_test_abst)
}
plot(mod_perplexity, xlab = "number of topics k", ylab = "perplexity")

#re-run model with entire data set with k = 20
abstlda20 <- LDA(dtm_abst, k = 20, method = "VEM", 
                 control = list(seed = 1234))
abstlda20


#model interpretation
#beta is the probability of that term belonging to that topic
tidy_lda20<- tidy(abstlda20)
tidy_lda20

#top 10 terms for each topic 
top_ten_terms <- tidy_lda20 %>% 
  group_by(topic) %>% 
  slice_max(beta, n = 10, with_ties = FALSE) %>% 
  ungroup() %>% 
  arrange (topic, -beta)
top_ten_terms
plot(top_ten_terms)

#plot top ten terms per topic
top_ten_terms %>%
  mutate(term = reorder_within(term, beta, topic)) %>%
  group_by(topic, term) %>%    
  arrange(desc(beta)) %>%  
  ungroup() %>%
  ggplot(aes(beta, term, fill = as.factor(topic))) +
  geom_col(show.legend = FALSE) +
  scale_y_reordered() +
  labs(title = "Top 10 terms in each LDA topic",
       x = expression(beta), y = NULL) +
  theme(axis.text = element_text(size = 6))+  
  facet_wrap(~ topic, ncol = 4, scales = "free")+
  theme(strip.text.x = element_text(size = 6))

ggplot(top_ten_terms[top_ten_terms$topic %in% "6", ], 
       aes(beta, reorder(term, beta), fill = as.factor(topic))) + 
  geom_col(show.legend = FALSE, fill = "brown") +
  scale_y_reordered() +
  labs(title = "Top 10 terms in topic 6",
       x = expression(beta), y = NULL) +
  theme(axis.text = element_text(size = 10))

#gamma is the probability of a document belonging to a topic
lda_doc <- tidy(abstlda20, matrix = "gamma")
lda_doc
ggplot(lda_doc, aes(gamma)) +
  geom_histogram(alpha = 0.8, bins = 15) +
  scale_y_log10() +
  labs(title = "Distribution of probabilities for all topics",
       y = "# of documents (log scale)", x = expression(gamma))

ggplot(lda_doc[lda_doc$topic %in% "19", ], aes(gamma)) +
  geom_histogram(alpha = 0.8, bins = 15) +
  scale_y_log10() +
  labs(title = "Distribution of probabilities for Topic19",
       y = "# of documents (log scale)", x = expression(gamma))

#connecting topics to article ids - topic most associated with each article
article_classification <- lda_doc %>%
  group_by(document) %>%
  slice_max(gamma) %>%
  ungroup()
article_classification

#connecting topic models with author keywords
#identify author keywords associated with each topic
lda_doc <- full_join(article_classification, ctsi_keyword, by = c("document" = "id"))
lda_doc1<- na.omit(lda_doc)
lda_doc1

#filter for document-topic prob >0.9
top_keywords <- lda_doc1 %>% 
  filter(gamma > 0.9) %>% 
  count(topic, keyword, sort = TRUE)
top_keywords

#plot top keywords associated with each topic
top_keywords %>%
  group_by(topic) %>%
  slice_max(n, n = 5, with_ties = FALSE) %>%
  ungroup %>%
  mutate(keyword = reorder_within(keyword, n, topic)) %>%
  ggplot(aes(n, keyword, fill = as.factor(topic))) +
  geom_col(show.legend = FALSE) +
  labs(title = "Top keywords for each LDA topic",
       x = "Number of documents", y = NULL) +
  theme(axis.text = element_text(size = 4))+
  scale_y_reordered() +
  facet_wrap(~ topic, ncol = 4, scales = "free")

top_keywords %>%
  group_by(topic) %>%
  slice_max(n, n = 5, with_ties = FALSE) %>%
  ungroup %>%
  mutate(keyword = reorder_within(keyword, n, topic))


###########################
#LDA using textminer#
###########################

# create dtm
set.seed(2100)
dtm <- CreateDtm(doc_vec = ctsi_abst$word,
                 doc_names = ctsi_abst$id,
                 stopword_vec = stopwords("en"),
                 verbose = F)

dtm <- dtm[,colSums(dtm)>2]

set.seed(1900)
mod_lda <- FitLdaModel(dtm = dtm,
                       k = 20, # number of topic
                       iterations = 500,
                       burnin = 180,
                       alpha = 0.1,beta = 0.05,
                       optimize_alpha = T,
                       calc_likelihood = T,
                       calc_coherence = T,
                       calc_r2 = F)

plot(mod_lda$log_likelihood,type = "l")
mod_lda$top_terms <- GetTopTerms(phi = mod_lda$phi,M = 15)
data.frame(mod_lda$top_terms)
mod_lda$coherence
mod_lda$prevalence <- colSums(mod_lda$theta)/sum(mod_lda$theta)*100
mod_lda$prevalence
mod_lda$summary <- data.frame(topic = rownames(mod_lda$phi),
                              coherence = round(mod_lda$coherence,3),
                              prevalence = round(mod_lda$prevalence,3),
                              top_terms = apply(mod_lda$top_terms,2,function(x){paste(x,collapse = ", ")}))

modsum <- mod_lda$summary %>%
  `rownames<-`(NULL)
modsum
modsum %>% pivot_longer(cols = c(coherence,prevalence)) %>%
  ggplot(aes(x = factor(topic,levels = unique(topic)), y = value, group = 1)) +
  geom_point() + geom_line() +
  facet_wrap(~name,scales = "free_y",nrow = 2) +
  theme_minimal() +
  labs(title = "Best topics by coherence and prevalence score",
       subtitle = "Publication abstracts",
       x = "Topics", y = "Value")

mod_lda$linguistic <- CalcHellingerDist(mod_lda$phi)
mod_lda$hclust <- hclust(as.dist(mod_lda$linguistic),"ward.D")
mod_lda$hclust$labels <- paste(mod_lda$hclust$labels, mod_lda$labels[,1])
plot(mod_lda$hclust)
rect.hclust(mod_lda$hclust, k = 5, border = 2:6) %>% 
  abline(h = 0.665, col = 'red')
plot(mod_lda$hclust)


