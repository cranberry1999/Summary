## ----------------------------------------------
#
# Society research
#
# COPYRIGHT (c) 2020 AJOU University
# Author	: Hwang, Hyereen
# History : 2020/04/20
#
## ----------------------------------------------


# 1. 데이터 불러오기 
data <- read.csv("사회조사_복지.사회참여.문화와여가.소득과소비.노동(제공)_2019_20200429_95711.csv", header=F)

# 2. 분석에 사용할 데이터 
## 사용할 데이터만 선택하기 
data <- data[,c(4,10,11,58,194,195,196,206:207,210:223)]

## 머리글행 이름 부여
colnames(data) <- c("성별", "주관적만족감", "성취만족도", "인간관계만족도", "소득만족도", "소비생활만족도", "직업선택요인", "일과가정_우선도", "고용의안정성", "종사상의지위", "임금근로자유형", "근로여건만족도_하는일", "근로여건만족도_인사관리" , "근로여건만족도_임금", "근로여건만족도_복리후생", "근로여건만족도_장래성", "근로여건만족도_근무환경", "근로여건만족도_인간관계","근로여건만족도_근로시간", "근로여건만족도_일-가정양립직장문화", "근로여건만족도_폭력방지노력", "근로여건만족도_교육훈련기회","일자리만족도")


# 3. 범주형 데이터 전처리
data$성별<-factor(data$성별, levels=c(1,2), labels=c("남자","여자"))
data$직업선택요인 <- factor(data$직업선택요인, levels=c(1:9), labels=c("명예", "안정성", "수입", "적성", "자아실현", "발전성", "근무시간", "기타", "잘모르겠다"))
data$일과가정_우선도 <- factor(data$일과가정_우선도, levels=c(1:5), labels=c("일", "대부분일", "둘다", "대부분가정", "가정"))
data$고용의안정성 <- factor(data$고용의안정성, levels=c(1:4), labels=c("매우많이느낌", "느끼는편", "느끼지않는편", "전혀느끼지않음"))
data$종사상의지위 <- factor(data$종사상의지위, levels=c(1:4), labels=c("임금근로자","고용원있는자영업자" , "고용원없는자영업자", "무급가족봉사자"))
data$임금근로자유형 <- factor(data$임금근로자유형, levels=c(1:3), labels=c("상용근로자", "임시근로자", "일용근로자"))



## NA값 제거 
data <- subset(data, !is.na(임금근로자유형)) #직장을 가지지 않은 사람들에 의해 NA로 표시된 값들은 분석과 연관이 없으므로 제거
data <- subset(data, !is.na(성취만족도)) #필수항목에서 NA로 나타난 값은 불성실한 응답에 의한 것이므로 제거


# 4. 요약 통계 확인
summary(data)

#일과 가정의 우선도를 묻는 문항- 응답이 "일"에 많이 치우치긴 했지만, 둘 다 중요하다는 응답이 6763으로 가장 많음
#대부분의 문항들과는 다르게 “소득만족도”와 “근로여건만족도_임금”에 관한 항목에서는 응답자들이 만족하지 않는다고 답함
#일자리에 대한 만족도- 만족과 보통에 많은 응답이 몰려 있음
#따라서 2019년 사회조사에서는 보다 많은 사람들이 본인의 일자리에 만족하고 있음을 알 수 있음
#앞으로의 과제에서 일자리만족과 가장 상관관계가 높은 변수가 무엇일지 찾아보아야 함 



## 데이터 저장
#write.csv(data, file = "society data.csv")

---

# 1. 분석환경 만들기 
# 1-1. 데이터 불러오기
data <- read.csv("society data.csv", header = T)

# 1-2. 패키지 불러오기 
library(ggplot2)
library(psych)
library(gmodels)




# 2. 데이터 시각화
## 2-1. 일과 가정_우선도 시각화 
ggplot(data, aes(일과가정_우선도)) + geom_bar(fill="#a3c4dc") + 
  labs(title="일과가정_우선도 시각화", subtitle="일과 가정생활의 균형을 중요시하는 직장인이 다수")
### 주제 선정 배경에서  "일`가정 모두 중요" 응답자가 "일이 우선" 응답자를 처음으로 넘으면서 "일을 우선시하던 사회에서 일과 가정생활의 균형을 중요시하는 사회로 변화하고 있다"고 서술한 내용을 직접 그래프로 나타내어 보았습니다. 일과 가정 둘 다 중요하다고 선택한 직장인의 수가 다른 항목보다 두배가량 높은 것을 알 수 있습니다. 


## 2-2. 개인생활 만족도와 일자리만족도의 관련성 시각화 
data$개인생활만족도 <- (data$주관적만족감+ data$성취만족도+ data$인간관계만족도+ data$소득만족도+ data$소비생활만족도)/5 # 분석에 사용하기 위해 개인의 생활에 관련된 컬럼들의 평균을 내어 "개인생활만족도"라는 새로운 컬럼을 만듭니다.
data <- data[,c(1:7,25,8:24)] # 순서를 다시 정렬합니다.

ggplot(data, aes(x=개인생활만족도, y=X, size = 일자리만족도)) +
  geom_point(alpha=0.2, aes(color=일자리만족도)) # 버블의 색깔과 크기를 통해 개인생활만족도가 낮아질수록 일자리만족도가 낮아지는 경향을 확인할 수 있으나, 데이터 수가 너무 많아 분석하기 힘듭니다.

# 총 14713개의 데이터 중에서 1000개를 임의추출하여 그래프를 다시 그려봅니다. 
set.seed(1004)
A <- sample(x=c(1:14713), size = 1000, replace = F)
data_point <- data[data$X %in% A,]
ggplot(data_point, aes(x=개인생활만족도, y=X, size = 일자리만족도)) +
  geom_point(alpha=0.2, aes(color=일자리만족도)) # 개인생활만족도에 해당하는 숫자가 커질수록 일자리만족도 버블의 색깔이 옅어지고 크기가 커지는 경향을 확인 할 수 있습니다. 따라서 개인생활 만족도와 일자리만족도는 연관이 있습니다. 

 
## 2-3. 근로 여건의 분야에 따른 근로 여건 만족도 분포
## 2-3-0. 그래프를 분석할 때 한눈에 파악하기 위해 데이터를 추가로 정제
### 수치형 속성을 범주형 속성으로 바꿔주었습니다.
data2 <- data
data2$근로여건만족도_하는일 <- factor(data2$근로여건만족도_하는일, levels=c(1:6), labels=c("매우만족", "약간만족", "보통", "약간불만족", "매우불만족", "모르겠다"))
data2$근로여건만족도_인사관리 <- factor(data2$근로여건만족도_인사관리, levels=c(1:7), labels=c("매우만족", "약간만족", "보통", "약간불만족", "매우불만족", "모르겠다", "해당없다"))
data2$근로여건만족도_임금 <- factor(data2$근로여건만족도_임금, levels=c(1:6), labels=c("매우만족", "약간만족", "보통", "약간불만족", "매우불만족", "모르겠다"))
data2$근로여건만족도_복리후생 <- factor(data2$근로여건만족도_복리후생, levels=c(1:7), labels=c("매우만족", "약간만족", "보통", "약간불만족", "매우불만족", "모르겠다", "해당없다"))
data2$근로여건만족도_장래성 <- factor(data2$근로여건만족도_장래성, levels=c(1:7), labels=c("매우만족", "약간만족", "보통", "약간불만족", "매우불만족", "모르겠다", "해당없다"))
data2$근로여건만족도_근무환경 <- factor(data2$근로여건만족도_근무환경, levels=c(1:6), labels=c("매우만족", "약간만족", "보통", "약간불만족", "매우불만족", "모르겠다"))
data2$근로여건만족도_인간관계 <- factor(data2$근로여건만족도_인간관계, levels=c(1:7), labels=c("매우만족", "약간만족", "보통", "약간불만족", "매우불만족", "모르겠다", "해당없다"))
data2$근로여건만족도_근로시간 <- factor(data2$근로여건만족도_근로시간, levels=c(1:6), labels=c("매우만족", "약간만족", "보통", "약간불만족", "매우불만족", "모르겠다"))
data2$근로여건만족도_일.가정양립직장문화 <- factor(data2$근로여건만족도_일.가정양립직장문화, levels=c(1:7), labels=c("매우만족", "약간만족", "보통", "약간불만족", "매우불만족", "모르겠다", "해당없다"))
data2$근로여건만족도_폭력방지노력 <- factor(data2$근로여건만족도_폭력방지노력, levels=c(1:7), labels=c("매우만족", "약간만족", "보통", "약간불만족", "매우불만족", "모르겠다", "해당없다"))
data2$근로여건만족도_교육훈련기회 <- factor(data2$근로여건만족도_교육훈련기회, levels=c(1:7), labels=c("매우만족", "약간만족", "보통", "약간불만족", "매우불만족", "모르겠다", "해당없다"))

## 2-3-1. 정제한 데이터를 바탕으로 데이터 시각화 
# 그래프에 대한 설명은 밑에 있습니다. 
ggplot(data2, aes(x=근로여건만족도_하는일, y=일자리만족도, color=근로여건만족도_하는일)) +
  geom_boxplot() + scale_x_discrete(limits=c("매우만족", "약간만족", "보통", "약간불만족", "매우불만족", "모르겠다"))

ggplot(data2, aes(x=근로여건만족도_인사관리, y=일자리만족도, color=근로여건만족도_인사관리)) + 
  geom_boxplot() + scale_x_discrete(limits=c("매우만족", "약간만족", "보통", "약간불만족", "매우불만족", "모르겠다", "해당없다"))

ggplot(data2, aes(x=근로여건만족도_임금, y=일자리만족도, color=근로여건만족도_임금)) + 
  geom_boxplot() + scale_x_discrete(limits=c("매우만족", "약간만족", "보통", "약간불만족", "매우불만족", "모르겠다"))

ggplot(data2, aes(x=근로여건만족도_복리후생, y=일자리만족도, color=근로여건만족도_복리후생)) + 
  geom_boxplot() + scale_x_discrete(limits=c("매우만족", "약간만족", "보통", "약간불만족", "매우불만족", "모르겠다", "해당없다"))

ggplot(data2, aes(x=근로여건만족도_장래성, y=일자리만족도, color=근로여건만족도_장래성)) + 
  geom_boxplot() + scale_x_discrete(limits=c("매우만족", "약간만족", "보통", "약간불만족", "매우불만족", "모르겠다", "해당없다"))

ggplot(data2, aes(x=근로여건만족도_근무환경, y=일자리만족도, color=근로여건만족도_근무환경)) + 
  geom_boxplot() + scale_x_discrete(limits=c("매우만족", "약간만족", "보통", "약간불만족", "매우불만족", "모르겠다"))

ggplot(data2, aes(x=근로여건만족도_인간관계, y=일자리만족도, color=근로여건만족도_인간관계)) + 
  geom_boxplot() + scale_x_discrete(limits=c("매우만족", "약간만족", "보통", "약간불만족", "매우불만족", "모르겠다", "해당없다"))

ggplot(data2, aes(x=근로여건만족도_근로시간, y=일자리만족도, color=근로여건만족도_근로시간)) + 
  geom_boxplot() + scale_x_discrete(limits=c("매우만족", "약간만족", "보통", "약간불만족", "매우불만족", "모르겠다"))

ggplot(data2, aes(x=근로여건만족도_일.가정양립직장문화, y=일자리만족도, color=근로여건만족도_일.가정양립직장문화)) + 
  geom_boxplot() + scale_x_discrete(limits=c("매우만족", "약간만족", "보통", "약간불만족", "매우불만족", "모르겠다", "해당없다"))

ggplot(data2, aes(x=근로여건만족도_폭력방지노력, y=일자리만족도, color=근로여건만족도_폭력방지노력)) + 
  geom_boxplot() + scale_x_discrete(limits=c("매우만족", "약간만족", "보통", "약간불만족", "매우불만족", "모르겠다", "해당없다"))

ggplot(data2, aes(x=근로여건만족도_교육훈련기회, y=일자리만족도, color=근로여건만족도_교육훈련기회)) + 
  geom_boxplot() + scale_x_discrete(limits=c("매우만족", "약간만족", "보통", "약간불만족", "매우불만족", "모르겠다", "해당없다"))
### 그래프로 나타내어 보니 그래프가 크게 두가지 분포로 나뉜다는 것을 알게 되었습니다. 근로여건 중에서 인사관리, 임금, 복리후생, 장래성, 근무환경, 일`가정 양립 직장문화, 교육훈련기회에 해당하는 그래프들은  근로여건에 대한 만족도가 상승하면서 일자리만족도 또한 상승하지만, 근로여건만족도가 매우 불만족인 지점에서는 일자리만족도가 보통에서 매우 불만족 사이에 분포하는 것을 알 수 있습니다. 
### 근로여건 중에서 근로시간, 하는 일, 인간관계, 폭력방지노력에 해당하는 그래프들은 앞서 설명한 그래프들보다 분포가 제각각이지만 근로여건 만족도가 상승함에 따라 일자리만족도도 상승하는 모습을 보입니다. 


# 3. 가설 검정
## 3-1. 상관분석 
corr.test(data$근로여건만족도_일.가정양립직장문화, data$주관적만족감)
corr.test(data$근로여건만족도_일.가정양립직장문화, data$일자리만족도)
### 일과 가정 모두 중요시하는 직장문화에 대한 만족도는 주관적만족감과 일자리만족도와 양의 상관관계가 있습니다. 앞서 그래프 1번에서 일과 가정 모두 중요하다고 선택한 직장인의 수가 많았던 만큼 그것을 존중해주는 직장문화가 조성되면 개인의 주관적만족감과 일자리만족도가 높아질 것입니다. 

corr.test(data$개인생활만족도, data$일자리만족도)
### 데이터 탐색과정에서 개인생활만족도와 일자리만족도의 연관성을 발견하였습니다. 따라서 두가지 요인을 상관분석 해본 결과 높은 양의 상관관계를 가지는 것으로 분석되었습니다. 이는 일자리만족도가 높아질 수록 개인생활만족도 또한 높아진다는 것을 의미합니다.


data_verify <- data[,c(8, 14:25)] # 한 번에 상관분석을 하기 위해서 필요한 컬럼을 선택합니다. 
corr.test(data_verify[,c(2:9,11:13)])
### 일자리만족도와 각 항목별 근로여건만족도는 모두 양의 상관관계를 가지는 것으로 분석되었습니다. 일자리만족도와 보다 높은 상관관계를 가지는 근로여건만족도는  순서대로 하는일, 근무환경, 폭력방지노력=임금, 근로시간, 인간관계, 복리후생=장래성, 교육훈련기회, 인사관리 입니다. 기업에서는 분석한 자료를 참고하여 근로여건을 개선한다면 직장인들의 일자리만족도를 향상시키는데 도움이 될 것입니다. 앞서 일자리만족도는 개인생활만족도와 양의 상관관계를 가진다는 것을 분석하였으므로 기업이 일자리만족도를 향상시키면 개인생활만족도도 향상되어 전보다 높은 워라밸을 누릴 수 있습니다.

#pairs.panels(data_analysis, scale=T, lm=T, stars=T)  코드를 읽는데 시간이 너무 오래걸려 주석처리 하였습니다. 내용은 보고서에 포함하였습니다. 

## 3-2. 편상관분석 
partial.r(cor(data_verify), c(2:13), 1)
### 앞의 분석에서 개인생활만족도가 일자리만족도와 상관관계가 있음을 알게 되었으므로 편상관분석을 시행하였습니다. 상관계수가 작아졌지만 편상관분석을 거쳐도 개인생활만족도와 일자리만족도는 양의 상관관계에 있음을 알 수 있습니다. 

# 다음 분석에 사용하기 위해 data저장
#write.csv(data, "Society Data2.csv")

---

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
   
