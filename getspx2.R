library(lubridate)
library(xts)

##Download the data from 
##http://realized.oxford-man.ox.ac.uk/media/1366/oxfordmanrealizedvolatilityindices.zip
##It contains the file OxfordManRealizedVolatilityIndices.csv.

rvi <- read.csv("OxfordManRealizedVolatilityIndices.csv",check.names=FALSE,skip=2)

rvi_date <- function(x) {
    year <- substring(x,1,4)
    month <- substring(x,5,6)
    day <- substring(x,7,8)
    as.Date(paste(year,month,day,sep="-"))
}

rvix <- xts(rvi[,-1],rvi_date(rvi$DateID))

SPX2<-rvix[,"SPX2.rv"]
save(SPX2,file=data="data/spx2.RData")
