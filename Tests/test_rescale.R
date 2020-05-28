dat <- do.call(rbind, 
               lapply(X = paste0("~sarahlotspeich/Downloads/rescale50/hn1/", list.files("~sarahlotspeich/Downloads/rescale50/hn1/")), FUN = read.csv)) %>% 
  dplyr::mutate(rescale_if = 0.5, hn = 1)

dat %<>% rbind(
  do.call(rbind, 
          lapply(X = paste0("~sarahlotspeich/Downloads/rescale50/hn2/", list.files("~sarahlotspeich/Downloads/rescale50/hn2/")), FUN = read.csv)) %>% dplyr::mutate(rescale_if = 0.5, hn = 2)
) 

dat %<>% rbind(
  do.call(rbind, 
          lapply(X = paste0("~sarahlotspeich/Downloads/rescale25/", list.files("~sarahlotspeich/Downloads/rescale25/")), FUN = read.csv)) %>% dplyr::mutate(rescale_if = 0.25, hn = 1)
) 

dat %<>% rbind(
  do.call(rbind, 
          lapply(X = paste0("~sarahlotspeich/Downloads/rescale90/hn1/", list.files("~sarahlotspeich/Downloads/rescale90/hn1/")), FUN = read.csv)) %>% dplyr::mutate(rescale_if = 0.90, hn = 1)
) 

dat %<>% rbind(
  do.call(rbind, 
          lapply(X = paste0("~sarahlotspeich/Downloads/rescale90/hn2/", list.files("~sarahlotspeich/Downloads/rescale90/hn2/")), FUN = read.csv)) %>% dplyr::mutate(rescale_if = 0.90, hn = 2)
) 

dat %>% dplyr::group_by(method, rescale_if, hn) %>% 
  dplyr::summarize(avg_beta0 = mean(beta0_est, na.rm = TRUE),
                   avg_beta1 = mean(beta1_est, na.rm = TRUE),
                   avg_beta2 = mean(beta2_est, na.rm = TRUE))

dat %>% dplyr::group_by(method, rescale_if, hn) %>%
  dplyr::filter(!(!conv & method == "SMLE"), !(is.na(conv) & method == "SMLE")) %>% 
  dplyr::summarize(avg_beta0 = mean(beta0_est, na.rm = TRUE),
                   avg_beta1 = mean(beta1_est, na.rm = TRUE),
                   avg_beta2 = mean(beta2_est, na.rm = TRUE),
                   se_beta0 = sd(beta0_est, na.rm = TRUE),
                   see_beta0 = mean(beta0_se, na.rm = TRUE),
                   se_beta1 = sd(beta1_est, na.rm = TRUE),
                   see_beta1 = mean(beta1_se, na.rm = TRUE),
                   se_beta2 = sd(beta2_est, na.rm = TRUE),
                   see_beta2 = mean(beta2_se, na.rm = TRUE),
                   sims = n(),
                   se_conv = mean(!is.na(beta1_se))) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(method = "") %>% 
  kableExtra::kable(digits = 3, row.names = FALSE, format = "latex") %>% 
  kableExtra::group_rows(group_label = "All validated", start_row = 1, end_row = 5) %>% 
  kableExtra::group_rows(group_label = "Complete-data", start_row = 6, end_row = 10) %>% 
  kableExtra::group_rows(group_label = "Naive", start_row = 11, end_row = 15) %>%
  kableExtra::group_rows(group_label = "SMLE", start_row = 16, end_row = 20)
  

