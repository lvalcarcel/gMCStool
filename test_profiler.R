gc()
cat("\f")

# Clean all object from workspace
rm(list = ls())
# browser()


library(rstudioapi) # Include it in helpers.R
setwd(dirname(getActiveDocumentContext()$path)) # Change into file directory

source("app.R")

runApp(app)



# DepMap.info.all.2 <- DepMap.info.all
# 
# DepMap.info.all.2$DepMapEssentialityByGene <- DepMap.info.all.2$DepMapEssentialityByGene %>% 
#   mutate(DepMap_ID = as.factor(DepMap_ID)) %>% 
#   mutate(ENSEMBL = as.factor(ENSEMBL)) %>% 
#   mutate(essentiality.database = as.factor(essentiality.database))
# 
# DepMap.info.all.2$DepMapCorrelationByGene <- DepMap.info.all.2$DepMapCorrelationByGene %>% 
#   mutate(gMCS.database = as.factor(gMCS.database)) %>% 
#   mutate(ENSEMBL = as.factor(ENSEMBL)) %>% 
#   mutate(essentiality.database = as.factor(essentiality.database)) %>% 
#   mutate(UNIT = as.factor(UNIT))
# 
# DepMap.info.all.2$DepMapGeneExpression <- DepMap.info.all.2$DepMapGeneExpression %>% 
#   mutate(ENSEMBL = as.factor(ENSEMBL)) %>% 
#   mutate(SYMBOL = as.factor(SYMBOL)) %>% 
#   mutate(DepMap_ID = as.factor(DepMap_ID)) %>% 
#   mutate(UNIT = as.factor(UNIT))
# 
# DepMap.info.all.2$DepMapExpressionByGene <- DepMap.info.all.2$DepMapExpressionByGene %>% 
#   mutate(gmcs.idx = as.factor(gmcs.idx)) %>% 
#   mutate(ENSEMBL = as.factor(ENSEMBL)) %>% 
#   mutate(DepMap_ID = as.factor(DepMap_ID)) %>% 
#   mutate(UNIT = as.factor(UNIT)) %>% 
#   mutate(gMCS.database = as.factor(gMCS.database))

sizes <- sapply(ls(), function(n) object.size(get(n)), simplify = FALSE)
print(sapply(sizes[order(as.integer(sizes))], function(s) format(s, unit = 'auto')))


sizes <- sapply(gMCS.info.raw$EssentialTasks_CultureMedium, function(n) object.size(n), simplify = FALSE)
print(sapply(sizes[order(as.integer(sizes))], function(s) format(s, unit = 'auto')))

sizes <- sapply(zz, function(n) object.size(n), simplify = FALSE); print(sapply(sizes[order(as.integer(sizes))], function(s) format(s, unit = 'auto')))

# Just for testing purposes --> Remove it when done
# install.packages("profvis")

library(profvis)
profvis({
  runApp(app)
})

# z <- lapply(dir()[grepl('.R', dir(), fixed = T) & !grepl("copy", dir(), ignore.case = T)], NCmisc::list.functions.in.file)
# zz <- unique(unlist(lapply(z,names))); zz <- lapply(zz, function(x)return(c())); names(zz) <- unique(unlist(lapply(z,names)));
# for(n in names(zz)){for(i in 1:length(z)){try({zz[[n]] <- union(zz[[n]], z[[i]][[n]])})}}
# zz
# 
# names(zz)


