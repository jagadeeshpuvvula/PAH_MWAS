dat<- read_csv("~/Documents/MWAS_home/cor_filter/1-Hydroxynapthalene.csv")

res<- read_csv("~/Documents/MWAS_home/pj_results/pah_mom/c18neg/c18neg_mom.csv")|>
  filter(Variable=="1-Hydroxynapthalene") |>
  filter(beta_dir %in% c("positive-significant", "negative-significant"))


filtered_res <- res %>% 
  inner_join(dat, by = c("mz", "time"))
