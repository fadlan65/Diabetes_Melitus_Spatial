library(sf)
library(tidyverse)
library(tidycensus)
library(tmap)
library(spdep)
library(tigris)
library(rmapshaper)
library(broom)
library(car)
library(spatialreg)
library(knitr)
library(stargazer)
library(spdep)
library(haven)
library(rgdal)
library(dplyr)
library(mapview)
library(stringr)
library(lubridate)
library(readr)

train_fktp <- data.frame(read_csv("D:/spasial/train_fktp.txt", 
                                  col_names = FALSE))
test_fktp <- data.frame(read_csv("D:/spasial/test_fktp.txt", 
                                 col_names = FALSE))
train_fktp<-train_fktp[,-26]
fktp <- rbind(train_fktp,test_fktp)

#Filter x7=provinsi
fktp_kaltengtimsel=filter(fktp, X7 %in% c(62,63,64))

#Mengambil data DMT2
fktp_kaltengtimsel_ICD10_E11 <- fktp_kaltengtimsel %>%
  filter(str_detect(X16, "E11"))
fktp_kaltengtimsel_ICD10_E11 <- fktp_kaltengtimsel_ICD10_E11[!duplicated(fktp_kaltengtimsel_ICD10_E11$X1),] %>%
  count(X8)
View(fktp_kaltengtimsel_ICD10_E11)
library(writexl)
#write_xlsx(fktp_kaltengtimsel_ICD10_E11,"D:/OneDrive/Documents/KULON/SMT 6/Analisis Data Spasial/dataprake11.xlsx")


#untuk uji moran variabel y
Indonesia <- st_read("D:/spasial/idn_adm_bps_20200401_shp/idn_admbnda_adm2_bps_20200401.shp")

library(readxl)
Data <- read_excel("D:/OneDrive/Documents/KULON/SMT 6/Analisis Data Spasial/dataprake11.xlsx")

kaltengtimsel <- subset(Indonesia, ADM1_EN %in%c('Kalimantan Tengah',
                                                 'Kalimantan Selatan','Kalimantan Timur'
))

kaltengtimsel <- kaltengtimsel %>%
  dplyr::select(ADM2_EN, geometry)

Prevalance_ss<- left_join(kaltengtimsel,Data, by="ADM2_EN")

#Peta
Indonesia.shp <- readOGR(dsn ="D:/spasial/idn_adm_bps_20200401_shp", 
                         layer="idn_admbnda_adm2_bps_20200401")

kaltengtimsel.shp <- Indonesia.shp[which(Indonesia.shp@data$ADM1_EN%in%c('Kalimantan Tengah','Kalimantan Selatan','Kalimantan Timur'
)),]
D <- as.matrix(dist(coordinates(kaltengtimsel.shp), method = "euclidean"))
head(D)

#Invers Distance Weights
#Invers weight matrix 
WM <- 1/D
#Inverse weight matrix - row-normalized 
diag(WM) <- 0
rtot <- rowSums(WM, na.rm = T)
WM_std <- WM/rtot
rowSums(WM_std, na.rm = T)
invers.WM <- mat2listw(WM_std, style = "W")
summary(invers.WM)

#Eksponential Distance Weights
#Eksponential weight matrix 
alpha <- 2
WM2 <- exp(-alpha*D)
#Eksponential weight matrix - row normalized 
diag(WM2) <- 0
rtot2 <- rowSums(WM2, na.rm = T)
rtot2
WM2_std <- WM2/rtot
rowSums(WM2_std, na.rm = T)
eksp.WM <- mat2listw(WM2_std, style = "W")
summary(eksp.WM)

# Calculate Moran's I
#Uji Autokorelasi 
moran.test(Prevalance_ss$n, invers.WM)
moran.test(Prevalance_ss$n, eksp.WM)
moran.test(Prevalance_ss$n, knn.MW)

#Plot Data
tmap_mode("view")
tmap_mode("plot")

View_I10<-tm_shape(Prevalance_ss,unit="mi",simplify = 1) +
  tm_basemap(leaflet::providers$Esri.WorldImagery)+
  tm_polygons(col = "n", style = "jenks",palette = "Reds", n=5,
              border.alpha = 0.3, title = "", size=0.5)+
  tm_text("ADM2_EN",size = 0.3, legend.size.show = T,col="black")+
  tm_view(set.view = c(118.641492,-2.943041,5),text.size.variable=T)+
  tm_minimap(leaflet::providers$CartoDB.VoyagerNoLabels,toggle = F)

lf<-tmap_leaflet(View_I10)
View_I10

tmap_save(View_I10,"h.png",width = 900, height = 900,asp=0)

#Data
library(readxl)
datalagmodel <- read_excel("KULON/SMT 6/Analisis Data Spasial/spasial lag model/datakaltengkaltimkalsel.xlsx", 
                           sheet = "Sheet4")
summary(datalagmodel)

#model regresi berganda biasa
reg.ols=lm(Y~X1+X2+X3+X4+X5+X6+X7+X8+X9, data=datalagmodel)
summary(reg.ols)
AIC(reg.ols)

#uji asumsi klasik regresi berganda
resid=reg.ols$residuals

library(nortest)
lillie.test(resid)

library(car)
vif(reg.ols)

library(lmtest)
bptest(reg.ols)

library(lmtest)
dwtest(reg.ols)

#Uji lagrange multiplier
reg.lm=lm.LMtests(reg.ols,eksp.WM,test=c("LMlag","LMerr"))
summary(reg.lm)

#model SLM
lagmodel=lagsarlm(Y~X1+X2+X3+X4+X5+X6+X7+X8+X9,data=datalagmodel,eksp.WM)
summary(lagmodel)

#uji asumsi model SLM
shapiro.test(lagmodel$residuals)
bptest.Sarlm(lagmodel)

