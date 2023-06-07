dat<- read_csv("~/Documents/MWAS_home/cor_filter/1-Hydroxynapthalene.csv")

res<- read_csv("~/Documents/MWAS_home/pj_results/pah_mom/c18neg/c18neg_mom.csv")|>
  filter(Variable=="1-Hydroxynapthalene") |>
  filter(beta_dir %in% c("positive-significant", "negative-significant"))


filtered_res <- res %>% 
  inner_join(dat, by = c("mz", "time"))

#============================================================================


runMummichog <- function(input_folder, output_folder) {
  # create object for storing data
  mSet3 <- InitDataObjects("mass_all", "mummichog", FALSE)
  
  # set peak format
  mSet3 <- SetPeakFormat(mSet3, "rmp")
  mSet3 <- UpdateInstrumentParameters(mSet3, 5.0, "mixed", "no")
  
  # get a list of text files in the input folder
  input_files <- list.files(input_folder, pattern = "\\.txt$", full.names = TRUE)
  
  # process each input file
  for (input_file in input_files) {
    # read the peak list data
    mSet3 <- Read.PeakListData(mSet3, input_file)
    mSet3 <- SanityCheckMummichogData(mSet3)
    
    # map selected adducts to current data
    mSet3 <- Setup.AdductData(mSet3, add.vec)
    mSet3 <- PerformAdductMapping(mSet3, add.mode="mixed")
    
    # perform mummichog algorithm using selected adducts, using version 2 of the mummichog algorithm
    mSet3 <- SetPeakEnrichMethod(mSet3, algOpt="mum", version="v2")
    mSet3 <- SetMummichogPval(mSet3, cutoff=0.05) # pval
    
    # the next step takes three or four minutes to run
    mSet3 <- PerformPSEA(mSet3, "hsa_mfn", "current", 3, 10000)
    
    # store the results as a data frame
    mummi_results <- as.data.frame(mSet3$mummi.resmat)
    
    # generate output file path and name
    output_file <- file.path(output_folder, paste0(tools::file_path_sans_ext(basename(input_file)), ".csv"))
    
    # save results as a CSV file
    write.csv(mummi_results, file = output_file, row.names = T)
    
    # print the output file path for each processed file
    cat("Processed file:", input_file, "Output file:", output_file, "\n")
  }
}