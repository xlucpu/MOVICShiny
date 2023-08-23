read_maf <- function (maf, clinicalData = NULL, removeDuplicatedVariants = TRUE, 
          useAll = TRUE, gisticAllLesionsFile = NULL, gisticAmpGenesFile = NULL, 
          gisticDelGenesFile = NULL, gisticScoresFile = NULL, cnLevel = "all", 
          cnTable = NULL, isTCGA = FALSE, vc_nonSyn = NULL, verbose = TRUE) 
{
  start_time = proc.time()
  if (is.data.frame(x = maf)) {
    maf = data.table::as.data.table(maf)
  }else {
    if (verbose) {
      cat("-Reading\n")
    }
    maf <- data.table::fread(file = maf, sep = "\t", stringsAsFactors = FALSE, 
                             verbose = FALSE, data.table = TRUE, showProgress = TRUE, 
                             header = TRUE, fill = TRUE, skip = "Hugo_Symbol", 
                             quote = "")
  }
  if (verbose) {
    cat("-Validating\n")
  }
  maf = maftools:::validateMaf(maf = maf, isTCGA = isTCGA, rdup = removeDuplicatedVariants, 
                    chatty = verbose)
  if (!useAll) {
    cat("--Using only `Somatic` variants from Mutation_Status. Set useAll = TRUE to include everything.")
    if (length(colnames(maf)[colnames(x = maf) %in% "Mutation_Status"]) > 
        0) {
      maf = maf[Mutation_Status %in% "Somatic"]
      if (nrow(maf) == 0) {
        stop("No more Somatic mutations left after filtering for Mutation_Status! Maybe set useAll to TRUE ?")
      }
    }else {
      cat("Mutation_Status not found. Assuming all variants are Somatic and validated\n")
    }
  }
  if (is.null(vc_nonSyn)) {
    vc.nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins", 
                     "Splice_Site", "Translation_Start_Site", "Nonsense_Mutation", 
                     "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins", 
                     "Missense_Mutation")
  }else {
    vc.nonSilent = vc_nonSyn
  }
  maf.silent = maf[!Variant_Classification %in% vc.nonSilent]
  if (nrow(maf.silent) > 0) {
    maf.silent.vc = maf.silent[, .N, .(Tumor_Sample_Barcode, 
                                       Variant_Classification)]
    maf.silent.vc.cast = data.table::dcast(data = maf.silent.vc, 
                                           formula = Tumor_Sample_Barcode ~ Variant_Classification, 
                                           fill = 0, value.var = "N")
    ##下面在得到summary.silent的过程中colSums函数会报错,提示maf.silent.vc.cast从第2列开始的列名太长了,需要在计算前进行修改,待计算后再把对应的列名改回来(variable names are limited to 10000 bytes)
    colnames_tmp <- colnames(maf.silent.vc.cast)
    colnames(maf.silent.vc.cast) <- c("Tumor_Sample_Barcode",paste0("V",2:ncol(maf.silent.vc.cast)))
    summary.silent = data.table::data.table(ID = c("Samples", 
                                                   colnames(maf.silent.vc.cast)[2:ncol(maf.silent.vc.cast)]), 
                                            N = c(nrow(maf.silent.vc.cast), colSums(maf.silent.vc.cast[, 
                                                                                                       2:ncol(maf.silent.vc.cast), with = FALSE])))
    summary.silent$ID <- c("Samples",colnames_tmp[2:length(colnames_tmp)])
    colnames(maf.silent.vc.cast) <- colnames_tmp
    maf = maf[Variant_Classification %in% vc.nonSilent]
    if (verbose) {
      cat(paste0("-Silent variants: ", nrow(maf.silent)), 
          "\n")
    }
  }
  if (nrow(maf) == 0) {
    stop("No non-synonymous mutations found\nCheck `vc_nonSyn`` argumet in `read.maf` for details")
  }
  if (!is.null(gisticAllLesionsFile)) {
    gisticIp = readGistic(gisticAllLesionsFile = gisticAllLesionsFile, 
                          gisticAmpGenesFile = gisticAmpGenesFile, gisticDelGenesFile = gisticDelGenesFile, 
                          isTCGA = isTCGA, gisticScoresFile = gisticScoresFile, 
                          cnLevel = cnLevel, verbose = verbose)
    gisticIp = gisticIp@data
    suppressWarnings(gisticIp[, `:=`(id, paste(Hugo_Symbol, 
                                               Tumor_Sample_Barcode, sep = ":"))])
    gisticIp = gisticIp[!duplicated(id)]
    gisticIp[, `:=`(id, NULL)]
    maf = data.table::rbindlist(list(maf, gisticIp), fill = TRUE, 
                                use.names = TRUE)
    maf$Tumor_Sample_barcode = factor(x = maf$Tumor_Sample_barcode, 
                                      levels = unique(c(levels(maf$Tumor_Sample_barcode), 
                                                        unique(as.character(gisticIp$Tumor_Sample_barcode)))))
  }else if (!is.null(cnTable)) {
    if (verbose) {
      cat("-Processing copy number data\n")
    }
    if (is.data.frame(cnTable)) {
      cnDat = data.table::copy(cnTable)
      data.table::setDT(x = cnDat)
    }else {
      cnDat = data.table::fread(input = cnTable, sep = "\t", 
                                stringsAsFactors = FALSE, header = TRUE, colClasses = "character")
    }
    colnames(cnDat) = c("Hugo_Symbol", "Tumor_Sample_Barcode", 
                        "Variant_Classification")
    if (isTCGA) {
      cnDat[, `:=`(Tumor_Sample_Barcode, substr(x = cnDat$Tumor_Sample_Barcode, 
                                                start = 1, stop = 12))]
    }
    cnDat$Variant_Type = "CNV"
    suppressWarnings(cnDat[, `:=`(id, paste(Hugo_Symbol, 
                                            Tumor_Sample_Barcode, sep = ":"))])
    cnDat = cnDat[!duplicated(id)]
    cnDat[, `:=`(id, NULL)]
    maf = data.table::rbindlist(l = list(maf, cnDat), fill = TRUE, 
                                use.names = TRUE)
    maf$Tumor_Sample_barcode = factor(x = maf$Tumor_Sample_barcode, 
                                      levels = unique(c(levels(maf$Tumor_Sample_barcode), 
                                                        unique(as.character(cnDat$Tumor_Sample_barcode)))))
  }
  maf$Variant_Classification = as.factor(as.character(maf$Variant_Classification))
  maf$Variant_Type = as.factor(as.character(maf$Variant_Type))
  if (verbose) {
    cat("-Summarizing\n")
  }
  mafSummary = maftools:::summarizeMaf(maf = maf, anno = clinicalData, 
                            chatty = verbose)
  #MAF在新版本的maftools(比如说2.12.0)中已经变成函数了,会报错；因此须使用版本为2.6.05的maftools包,这时的MAF定义为一个对象
  m = maftools:::MAF(data = maf, variants.per.sample = mafSummary$variants.per.sample, 
          variant.type.summary = mafSummary$variant.type.summary, 
          variant.classification.summary = mafSummary$variant.classification.summary, 
          gene.summary = mafSummary$gene.summary, summary = mafSummary$summary, 
          maf.silent = maf.silent, clinical.data = mafSummary$sample.anno)
  if (verbose) {
    cat("-Finished in", data.table::timetaken(start_time), 
        "\n")
  }
  return(m)
}
