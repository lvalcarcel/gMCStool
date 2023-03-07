SaveResultsInExcel <- function(Results, gMCS.info, filename){
  
  wb <- createWorkbook()
  ## Add worksheets
  addWorksheet(wb, "single met tasks")
  addWorksheet(wb, "ratio met tasks")
  addWorksheet(wb, "gMCS_single")
  ## Add the data
  writeDataTable(wb, "single met tasks", Results$num.essential.gene, tableStyle = "TableStyleLight1")
  writeDataTable(wb, "ratio met tasks", Results$ratio.essential.gene, tableStyle = "TableStyleLight1")
  aux <- Results$list.gMCS.essential
  aux$gMCS_ENSEMBL <- gMCS.info$gMCSs.ENSEMBL.txt[as.numeric(as.character(aux$gMCS))]
  aux$gMCS_SYMBOL <- gMCS.info$gMCSs.SYMBOL.txt[as.numeric(as.character(aux$gMCS))]
  writeDataTable(wb, "gMCS_single", aux, tableStyle = "TableStyleLight1", colNames = T, rowNames = F)
  
  saveWorkbook(wb, filename, overwrite = TRUE)
  
}
