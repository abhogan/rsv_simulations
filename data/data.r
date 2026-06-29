################
##### data #####
################

# https://www.nature.com/articles/s41598-021-88524-w
# https://pubmed.ncbi.nlm.nih.gov/26732801/

############################
##### attack rate data #####
############################

# https://pmc.ncbi.nlm.nih.gov/articles/PMC11431910/
# https://publications.ersnet.org/content/erj/57/4/2002688.abstract?implicit-login=true%26428

ages_ar <- c(1, 1, 2, 3, 1, 2, 1, 2, 70, 70)
ages_chr_ar <- c("0-1", "0-1", "0-2", "0-3", "0-2", "0-3", "0-2", "0-3", "median 75", "median 75")
attack_rate <- c(0.534, 0.37, 0.68, 0.86, 0.441, 0.846, 0.38, 0.68, 22/527, 37/513)
ref_ar <- c("Cacho et al", "Kutsaya et al", "Kutsaya et al", "Kutsaya et al", "Andeweg et al", "Andeweg et al", "Kazakova et al", "Kazakova et al", "Korsten et al", "Korsten et al")

#################################################
##### time of peak incidence and amplitude ######
#################################################

# https://www.health.gov.au/resources/publications/australian-respiratory-surveillance-report-8-september-to-21-september-2025?language=en
# figure 10 - extracted using webplotdigitizer for year 2024
# all ages

case_times_pia <- 1:53

cases_pia <- c(961, 1030, 1099, 1133, 1442, 1717, 2094, 2575,
               3159, 3571, 3983, 4395, 4498, 5116, 5494, 6009,
               5391, 5562, 5803, 6661, 7382, 7210, 6970, 6180,
               6318, 6112, 5700, 5185, 4807, 4326, 3983, 3708,
               3330, 3159, 3090, 2712, 2300, 2129, 1785, 1579,
               1236, 1202, 1133, 1030, 1099, 1064, 1202, 1236,
               1270, 1442, 1511, 1099, 1202)

amplitude <- (max(cases_pia) - min(cases_pia))/max(cases_pia)

peak_time <- case_times_pia[cases_pia == max(cases_pia)]

###########################
##### total incidence #####
###########################

# year 2024
# https://nindss.health.gov.au/pbi-dashboard/

ages_ti <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85)

cases_ti <- c(86286, 14799, 7578, 4031, 2901, 3313, 4266, 4296, 3699, 3607, 4536, 4444, 5176, 5036, 5209, 5345, 4567, 6822)

##### attack rate in older adults ######

