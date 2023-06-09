library(Biobase)
library(corrplot)

# Type of data
df <- antiProfilesData::apColonData
class(df)

pdata = pData(df)
edata = as.data.frame(exprs(df))
fdata = fData(df)

#1.a
str(edata)
str(pdata)
str(fdata)

#1.b
colnames(edata)
row.names(edata)

colnames(pdata)
row.names(pdata)

#1.c
summary(edata)
summary(pdata$Status)

#1.d
for (col in colnames(pdata[,1:6])) {
  print(col)
  tbl <- table(pdata[ ,col], useNA='always')
  print(tbl)
}

#1.e
cov(edata[,1:10], method= 'spearman')
correlation <- cor(edata[,1:10], method="spearman")
correlation

corrplot(correlation, type = "upper", 
         tl.col = "black", tl.srt = 85)

col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = correlation, col = col, symm = TRUE)

#1.f
plot(edata$GSM95478~edata$GSM95473, col=5, main = "Relation betweeen GSM95478 and GSM95473", pch=19, xlab="GSM95463", ylab='GSM95478')
fit1 <- lm(edata$GSM95478~edata$GSM95473)
abline(fit1, type='l')

#2
edata = log2(abs(edata))*sign(edata)
pc1 = prcomp(edata)

edata_centered = t(t(edata) - colMeans(edata))
svd_s = svd(edata_centered)

svd_s$v[1:5,1:5]
pc1$rotation[1:5,1:5]

plot(pc1$rotation[,1],svd_s$v[,1],col=2)

#3
visual_artists <- as.factor(c(rep('Aries', 29), rep('Taurus', 24), rep('Gemini', 22), rep("Cancer", 19), rep("Leo", 21),
                              rep("Virgo", 18), rep("Libra", 19), rep("Scorpio", 20), rep("Sagittarius", 23), rep("Capricorn", 18),
                              rep("Aquarius", 20), rep("Pisces", 23)))
p = c(1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1/12)
p = p / sum(p)
chisq.test(table(visual_artists), p=p)
print("P value = 0.9265 greater than 0.05 so we will accept null hypothesis and reject alternative hypothesis")

#4
dist1 = dist(t(edata[,1:10]))
hclust1 = hclust(dist1)
plot(hclust1, hang = -1)

kmeans1 = kmeans(edata, centers = 3)
table(kmeans1$cluster)

kmeans2 = kmeans(edata, centers = dim(edata)[2])
table(kmeans2$cluster)


