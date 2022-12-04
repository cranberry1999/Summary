## ----------------------------------------------
#
# HW3. Association and Text analysis of society data
#
# COPYRIGHT (c) 2020 AJOU University
# Author	: Hwang, Hyereen
# History : 2020/06/01
#
## ----------------------------------------------


# 0. 데이터 불러오기 
data <- read.csv("Society Data2.csv"); data <- data[,-1] #분석에 사용하지 않을 컬럼 삭제
str(data) #데이터를 확인
summary(data)
# 데이터를 불러오는 과정에서 변형된 변수속성 다시 factor로 변환
for (i in 9:24) {
  data[,i] <- as.factor(data[,i])
} 


### 1. 연관분석 
# 연관분석 할 데이터 정제 
adata <- data; adata <-adata[,-c(1,3:8,12:13)] #분석에 사용하지 않을 컬럼 제외 
# 범주형으로 데이터 정제
adata$일자리만족도 <- factor(adata$일자리만족도, levels=c(1:5),
                       labels=c('매우만족', '약간만족', '보통', '약간불만족', '매우불만족'))
# 사용할 패키지 불러오기 
library(arules)
library(arulesViz)
rule <- apriori(adata, parameter=list(minlen=2)) #처음 코드를 실행할 때 support=0.025, confidence=0.2 로 실행해보았으나 연관규칙이 너무 많이 나와서 디폴트값으로 코드를 실행하였습니다
summary(rule)
# 시각화 
plot(rule)
# 응답패턴 분석 
rule.df <- as(rule, 'data.frame')

## 1-1. 일자리 만족도에 영향을 주는 근로요건 분석 
rule1 <- apriori(adata, parameter = list(minlen=2, support=0.04, confidence=0.5))# 디폴트값으로 코드를 실행하였더니 연관규칙을 찾을 수 없어서 support, confidence값을 낮추었습니다 
work.df <- inspect(sort(subset(rule1, subset=rhs %in% c('일자리만족도=매우만족'))), by=lift)
# 일자리 만족도와 강한 연관관계를 가지는 응답규칙을 발견하기 위하여 일자리만족도를 매우만족에 두고 연관분석을 하였습니다.


### 2. 텍스트 분석 
#패키지 불러오기 
library(KoNLP)
library(stringr)
library(wordcloud2)
library(qgraph)
library(tm)

options(mc.cores=1)

dics <- c('sejong','woorimalsam') #사전을 선택 

# 한글 형태소 분석 함수
ko.words <- function(doc) {
  d <- as.character(doc)
  ## 형태소 중 명사, 용언 추출하기    
  pos <- paste(SimplePos09(d))
  ex <- str_match(pos, '([가-힣]+)/[NP]')
  keyword <- ex[,2]
  ## 결측값(NA) 정제하기    
  keyword <- keyword[!is.na(keyword)]
  paste(keyword, collapse = ' ')
}

# 텍스트 파일 불러오기 
txt <- readLines("work.txt")
#각 줄마다 형태소 분석 
words <- lapply(txt, ko.words) 

# 단어 정제
words <- gsub("위하", NA, words)
words <- gsub("따르", NA, words)
words <- gsub("불리", NA, words)
words <- gsub("아니", NA, words)

# tdm을 생성 
cps <- Corpus(VectorSource(words))
tdm <- TermDocumentMatrix(cps, control = list(removePunctuation=T, removeNumbers=T))
tdm.matrix <- as.matrix(tdm)


# 단어의 빈도수를 분석하여 키워드를 찾습니다 
word.count <- rowSums(tdm.matrix)
word.order <- order(word.count, decreasing = T)
freq.words <- tdm.matrix[word.order[1:50],] #빈도수 상위 50개의 단어를 추출하였습니다. 
freq.words <- subset(freq.words, str_length(row.names(freq.words))>=2) #단어의 길이가 2이상인 단어를 출력 
keyword.df <- data.frame(rownames(freq.words), rowSums(freq.words))


## 2-1. 워드클라우드 시각화
library(RColorBrewer) #RColorBrewer패키지를 통해 워드클라우드에 색을 입혀주었습니다 
wordcloud2(keyword.df, minRotation=0, maxRotation=0, color=rep(brewer.pal(8, "Dark2"))) # 단어가 회전하지 않게 파라미터에 minRotation=0, maxRotation=0 을 추가 

## 2-2. 네트워크 시각화
co.matrix <- freq.words %*% t(freq.words)
qgraph(co.matrix, labels=rownames(co.matrix), diag=F, layout='spring', edge.color='blue')       
   
