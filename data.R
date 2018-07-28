library(Hmisc)
library(knitr)
library(bizdays)
library(tidyverse)

#getwd()

setwd("/Volumes/yilunsun/Datasets/Feb2018csv")

L4 = read_csv("L4.csv")

L5 = read_csv("L5.csv")

LRTD = read_csv("LRTD.csv")

L6 = read_csv("L6.csv")

PT_KEEP = read_csv("PT_KEEP.csv")

L6_KPS = read_csv("L6_ARCHIVE.csv")

DVH = read_csv("DVH.csv")

DICOM = read_csv("DICOM_TRANS.csv")

DVH_LONG_LU = read_csv("DVH_LONG_LU.csv")

L8 = read_csv("L8.csv")

L9 = read_csv("L9.csv")

###########################################################################################
### data manipulation
###########################################################################################
L6_KPS = L6_KPS %>%
  select(nid, Karnofsky_KPS_L6) %>%
  rename(kps = Karnofsky_KPS_L6)


L4 = L4 %>% 
  mutate(comorbid_count = (Does_HYP_L4 == "Yes") + (Does_DM_L4 == "Yes") +
           (Does_SCL_L4 == "Yes") + (Does_RA_L4 == "Yes") + (Does_LUP_L4 == "Yes") + 
           (Does_CD_L4 == "Yes") + (Does_CPD_L4 == "Yes") + (Does_CHF_L4 == "Yes") +
           (Does_CTD_L4 == "Yes") + (Does_DEM_L4 == "Yes") + (Does_HEM_L4 == "Yes") +
           (Does_LEUK_L4 == "Yes") + (Does_ML_L4 == "Yes") + (Does_MI_L4 == "Yes") +
           (Does_PVD_L4 == "Yes") + (Does_UD_L4 == "Yes") + (Does_LD_L4 == "Yes") +
           (Does_RD_L4 == "Yes") + (Does_MST_L4 == "Yes")) %>%
  mutate(comorbid_count = ifelse(is.na(comorbid_count),0,comorbid_count)) %>%
  mutate(hist_collapsed = ifelse(Histology_L4 == "SCLC", "SCLC", ifelse(is.na(Histology_L4), NA, "NSCLC")))


L6 = L6 %>%
  select(nid, WeightLoss_L6, AdverseEvent_11_L6) %>%
  rename(ecog_L6 = AdverseEvent_11_L6) %>%
  full_join(L6_KPS, by = "nid") %>%
  mutate(ecog_L6 = ifelse(!is.na(ecog_L6), ecog_L6, 
                          5 - cut(kps, breaks = c(seq(5,85,by=20), Inf),labels = FALSE, include.lowest = T, right= F))) 


PT = PT_KEEP %>% 
  select(nid, age, gender, race, insur_1, insur_2, marriage_1, marriage_2)


PT2 = PT_KEEP %>% 
  filter(cancer =="lung") %>%
  mutate(date_first_data_entry = as.Date(date_first_data_entry, format= "%m/%d/%y")) %>%
  filter(date_first_data_entry > as.Date("01/01/2015", format = "%m/%d/%Y"), date_first_data_entry < as.Date("12/31/2015", format = "%m/%d/%Y"))


QUARTILE = PT2 %>%
  inner_join(L4 %>% 
               select(nid, hist_collapsed), by="nid") %>%
  group_by(hosp) %>%
  summarise(vol_quartile = n()/12) %>%
  mutate(volume_quartile = ntile(vol_quartile,4) -1)


DVH_L = DVH %>%
  select(nid, cancer, strctname, voltotal, meangy) %>%
  filter(cancer == "lung", strctname %in% c("PrimaryTumor_PTV", "normallung", "Esophagus", "heart"), !is.na(meangy)) %>%
  select(-cancer) %>%
  distinct(nid, strctname, .keep_all = TRUE)


DVH_L = DVH_L %>% 
  filter(strctname == "PrimaryTumor_PTV") %>% 
  rename(ptv_voltotal = voltotal, meangy_ptv = meangy) %>%
  select(-strctname) %>%
  full_join(DVH_L %>% 
              filter(strctname == "normallung") %>% 
              select(-voltotal) %>%
              rename(meangy_lung = meangy) %>%
              select(-strctname), by = "nid") %>%
  full_join(DVH_L %>% 
              filter(strctname == "Esophagus") %>% 
              rename(esoph_voltotal = voltotal, meangy_esoph = meangy) %>%
              select(-strctname), by = "nid") %>%
  full_join(DVH_L %>% 
              select(-voltotal) %>%
              filter(strctname == "heart") %>% 
              rename(meangy_heart = meangy) %>%
              select(-strctname), by = "nid") 

### check DVH_L
#psych::describe(DVH_L[,c("meangy_ptv", "meangy_lung", "meangy_esoph", "meangy_heart")])


DICOM = DICOM %>% 
  rename(TPN3Dnum = tpn3d) %>%
  select(nid, ppn3d_cat, TPN3Dnum, tnfxp)


DVH_L2 = DVH %>% 
  filter(cancer == "lung", voltotal != "") %>% ### check voltotal value
  select(nid, cancer, struclocal, voltotal) %>%
  select(-cancer) %>%
  distinct(nid, struclocal, .keep_all = TRUE)


D_and_V = DVH_LONG_LU %>%
  select(nid, struclocal, dose, volume) %>%
  filter(struclocal %in% c("primarytumor_ptv", "Heart", "lungs_r+l-gtv", "lungs_r+l-itv", "lungs_rl-ctv", "Esophagus"), !is.na(volume), !is.na(dose)) %>%
  full_join(DVH_L2, by = c("nid", "struclocal")) %>%
  mutate(dose = ifelse(dose >200, dose / 100, dose),
         vol95_diff = abs(volume - 95), 
         dose5_diff = abs(dose - 5),
         dose20_diff = abs(dose - 20),
         vol2cc_diff = abs(volume * voltotal / 100 - 2))


D95 = D_and_V %>%
  group_by(nid, struclocal) %>%
  slice(which.min(vol95_diff)) %>%
  rename(D95 = dose) %>%
  select(nid, struclocal, D95)


V5 = D_and_V %>%
  group_by(nid, struclocal) %>%
  slice(which.min(dose5_diff)) %>%
  rename(V5 = volume) %>%
  select(nid, struclocal, V5)


V20 = D_and_V %>%
  group_by(nid, struclocal) %>%
  slice(which.min(dose20_diff)) %>%
  rename(V20 = volume) %>%
  select(nid, struclocal, V20)


D2CC = D_and_V %>%
  group_by(nid, struclocal) %>%
  slice(which.min(vol2cc_diff)) %>%
  rename(D2CC = dose) %>%
  select(nid, struclocal, D2CC)


DOSIMETRICS = D95 %>%
  filter(struclocal == "primarytumor_ptv") %>%
  ungroup() %>%
  rename(D95_ptv = D95) %>%
  select(-struclocal) %>%
  full_join(D2CC %>% 
              filter(struclocal == "Esophagus") %>% 
              ungroup() %>%
              rename(D2CC_esoph = D2CC) %>% 
              select(-struclocal), by = "nid") %>%
  full_join(V5 %>%
              filter(struclocal == "Heart") %>% 
              ungroup() %>%
              rename(V5_heart = V5) %>% 
              select(-struclocal), by = "nid") %>%
  full_join(V5 %>%
              filter(struclocal == "lungs_r+l-gtv") %>% 
              ungroup() %>%
              rename(V5_lung_gtv = V5) %>% 
              select(-struclocal), by = "nid") %>%
  full_join(V5 %>%
              filter(struclocal == "lungs_r+l-itv") %>% 
              ungroup() %>%
              rename(V5_lung_itv = V5) %>% 
              select(-struclocal), by = "nid") %>%
  full_join(V5 %>%
              filter(struclocal == "lungs_rl-ctv") %>% 
              ungroup() %>%
              rename(V5_lung_ctv = V5) %>% 
              select(-struclocal), by = "nid") %>%
  full_join(V20 %>%
              filter(struclocal == "lungs_r+l-gtv") %>% 
              ungroup() %>%
              rename(V20_lung_gtv = V20) %>% 
              select(-struclocal), by = "nid") %>%
  full_join(V20 %>%
              filter(struclocal == "lungs_r+l-itv") %>% 
              ungroup() %>%
              rename(V20_lung_itv = V20) %>% 
              select(-struclocal), by = "nid") %>%
  full_join(V20 %>%
              filter(struclocal == "lungs_rl-ctv") %>% 
              ungroup() %>%
              rename(V20_lung_ctv = V20) %>% 
              select(-struclocal), by = "nid") %>%
  mutate(V5_lung = coalesce(V5_lung_gtv, V5_lung_itv, V5_lung_ctv),
         V20_lung = coalesce(V20_lung_gtv, V20_lung_itv, V20_lung_ctv))

### still need to check why number of observations is different from SAS ver.
create.calendar(name='WeekendsOnly', weekdays=c('sunday', 'saturday'))

A = L4 %>%
  full_join(L5, by="nid") %>%
  full_join(LRTD %>%
              rename(submitdate_lrtd = submitdate), by="nid") %>%
  left_join(PT_KEEP %>%
              select(nid, rt_st_Date, rt_end_Date) %>%
              mutate(rt_st_Date = as.Date(as.character(rt_st_Date),format="%d%b%Y"),
                     rt_end_Date = as.Date(as.character(rt_end_Date),format="%d%b%Y")) %>%
              full_join(DVH_L %>% 
                          select(nid, meangy_ptv, meangy_lung), by="nid") %>%
              full_join(DICOM %>%
                          select(nid, TPN3Dnum, tnfxp), by="nid"), by="nid") %>%
  filter(!is.na(rt_end_Date)) %>%
  distinct(nid, .keep_all = T) %>%
  mutate(hist = substr(Histology_L4,1,5),
         rt_days = bizdays(rt_st_Date,rt_end_Date,"WeekendsOnly"),
         BID_plans = gsub("NA","",paste0(BID1_LRTD, BID2_LRTD, BID3_LRTD, BID4_LRTD, BID5_LRTD))) %>%
  rename(BID = BID1_LRTD) %>%
  mutate(BID = ifelse(!is.na(BID), BID, ifelse(tnfxp <= rt_days, "N", ifelse(tnfxp/rt_days > 1.75, "Y", "N"))),
         BID = ifelse(!is.na(BID), BID, "N"),
         BID_plans = ifelse(BID_plans != '', BID_plans, BID),
         excl_rs = ifelse(is.na(hist_collapsed), "1. Histology missing", NA),
         excl_rs = ifelse(is.na(excl_rs) & is.na(SurgicalResection_L4), "2. Surgery", excl_rs),
         excl_rs = ifelse(is.na(excl_rs) & SurgicalResection_L4 != "N" & !is.na(SurgicalResection_L4), "2. Surgery", excl_rs),
         excl_rs = ifelse(is.na(excl_rs) & HeterogeneityCorrect_LRTD == "N", "3. No hetero correction", excl_rs),
         excl_rs = ifelse(is.na(excl_rs) & is.na(group_stage), "4b. Stage unknown", excl_rs),
         excl_rs = ifelse(is.na(excl_rs) & !group_stage %in% c('IIA', 'IIB', 'IIIA', 'IIIB'), "4a. Stage I", excl_rs),
         excl_rs = ifelse(is.na(excl_rs) & is.na(meangy_ptv), "5. Dose data missing or anomalous", excl_rs),
         excl_rs = ifelse(is.na(excl_rs) & !is.na(meangy_ptv) & (meangy_ptv > 200 | meangy_ptv < 30 | meangy_lung > 50) , "5. Dose data missing or anomalous", excl_rs),
         excl_rs = ifelse(is.na(excl_rs) & is.na(TPN3Dnum), "6. Unknown treatment plan", excl_rs),
         excl_rs = ifelse(is.na(excl_rs) & BID == 'Y', "7. BID Yes", excl_rs)
  )


### this data has excluded ineligible patients*/
SUBGROUP = A %>%
  filter(is.na(excl_rs)) %>%
  select(nid, hosp.x, hist_collapsed, SurgicalResection_L4, HeterogeneityCorrect_LRTD,
         Chemotherapy_L5, IfYesChemo_Con_L5, IfYesChemo_Adj_L5, IfYesChemo_Neo_L5,
         comorbid_count, Does_HYP_L4, Does_DM_L4, Does_SCL_L4, Does_RA_L4, Does_LUP_L4,
         Does_CD_L4, Does_CPD_L4, Does_CHF_L4, Does_CTD_L4, Does_DEM_L4, Does_HEM_L4,
         Does_LEUK_L4, Does_ML_L4, Does_MI_L4, Does_PVD_L4, Does_UD_L4, Does_LD_L4, 
         Does_RD_L4, Does_MST_L4, BMI_L4, CurrentSmoker_L4, FormerSmoker_L4, Prior_02_L4,
         Liters_L4, PFTDone_L4, PFTValues_1_L4, PFTValues_2_L4, DLCDone_L4, DLCValues_1_L4,
         DLCValues_2_L4, TumorStage_T_L4, TumorStage_N_L4, group_stage, MotionManagement_LRTD,
         AvoidanceStructures__5_LRTD,  ### heart
         AvoidanceStructures__7_LRTD,  ### esophagus
         PTVStructures_1_LRTD, ### Spinal cord within 2cm of PTV 
         PTVStructures_2_LRTD, ###  heart within 2cm of PTV 
         PTVStructures_3_LRTD, ### esoph within 2cm of PTV 
         PTVStructures_4_LRTD, ###  Brachial plexus within 2cm of PTV 
         PTVStructures_5_LRTD, ###  other within 2cm of PTV
         OtherSpecify_3_LRTD, ###  what other structure within 2cm of PTV
         BID, BID_plans, tnfxp, rt_days) %>%
  rename(hosp = hosp.x) %>%
  mutate(academic = ifelse(hosp %in% c(151,152,153,154,155), 1, 0),
         sub = 1) %>%
  rename(ptv_2cm_spinalcord = PTVStructures_1_LRTD,
         ptv_2cm_heart = PTVStructures_2_LRTD,
         ptv_2cm_esoph = PTVStructures_3_LRTD,
         ptv_2cm_brachialplexus = PTVStructures_4_LRTD,
         ptv_2cm_other = PTVStructures_5_LRTD,
         ptv_2cm_otherstruct = OtherSpecify_3_LRTD) %>%
  left_join(PT %>%
              full_join(DVH_L, by='nid') %>%
              full_join(DOSIMETRICS, by='nid') %>%
              full_join(L6, by='nid') %>%
              full_join(DICOM, by = 'nid'), by = "nid") %>%
  distinct(nid, .keep_all=T)


### consort diagram
cat('Reasons patients are excluded')
table(A$excl_rs)


###########################################################################################
### Esoph Tox
###########################################################################################

TOX = read_csv("LU_TOXICITY_ALL.csv")

ESOPH_LONG = TOX %>%
  filter(!is.na(Esophagitis)) %>%
  mutate(dateofeval = as.Date(dateofeval, format='%d%b%Y'),
         rt_st_Date = as.Date(rt_st_Date, format='%d%b%Y'),
         rt_end_Date = as.Date(rt_end_Date, format='%d%b%Y'))


G2onset = ESOPH_LONG %>%
  filter(Esophagitis>1) %>%
  group_by(nid) %>%
  summarise(freq = n(),
            G2onset = min(dateofeval)) %>%
  arrange(nid)


ESOPH_END = ESOPH_LONG %>%
  filter(timepoint == '3-EOT' & !is.na(Esophagitis)) %>%
  mutate(esoph_end = Esophagitis) %>%
  select(nid, timepoint, rt_st_Date, dateofeval, esoph_end) 


TEMP2 = ESOPH_LONG %>%
  group_by(nid) %>%
  summarise(freq_temp2 = n(),
            rtend = max(rt_end_Date),
            lastfu = max(dateofeval),
            maxesoph = max(Esophagitis))


ESOPH2 = TEMP2 %>%
  full_join(ESOPH_END, by = 'nid') %>%
  full_join(G2onset, by = 'nid') %>%
  mutate(time_postrt = as.numeric(lastfu-rtend)/(365/12),
         keep1mo = ifelse((!is.na(esoph_end)) | (!is.na(time_postrt) & time_postrt >= .5 & time_postrt <= 4), 1, 0),
         g2onset_week = as.numeric(G2onset-rt_st_Date)/7)


###########################################################################################
### Pneumonitis Tox
###########################################################################################
### caution! check time_postrt definition
TEMP1 = TOX %>%
  filter(!is.na(Pneumonitis)) %>%
  mutate(rt_st_Date = as.Date(rt_st_Date, format="%d%b%Y"),
         dateofeval = as.Date(dateofeval, format="%d%b%Y"),
         pneum2plus = ifelse(Pneumonitis >=2,1,0),
         time_postrt = as.numeric(dateofeval-rt_st_Date)/(365/12),
         pneum_time = ifelse(time_postrt >= (1.5+.5) & time_postrt < (1.5+2), 1, 
                             ifelse(time_postrt >= (1.5+2) & time_postrt < (1.5+4.5), 3,
                                    ifelse(time_postrt >= (1.5+4.5), 6, NA)))) %>%
  select(nid, rt_end_Date, rt_st_Date, dateofeval, Pneumonitis, pneum2plus, time_postrt, pneum_time)


M1 = TEMP1 %>%
  filter(pneum_time==1 & !is.na(Pneumonitis)) %>%
  mutate(rp_1mo = Pneumonitis,
         have_rp_1mo = 1,
         time2 = 1) %>%
  select(nid, rp_1mo, have_rp_1mo, time2)


M3 = TEMP1 %>%
  filter(pneum_time==3 & !is.na(Pneumonitis)) %>%
  mutate(rp_3mo = Pneumonitis,
         have_rp_3mo = 1,
         time2 = 3) %>%
  select(nid, rp_3mo, have_rp_3mo, time2)


M6 = TEMP1 %>%
  filter(pneum_time==6 & !is.na(Pneumonitis)) %>%
  mutate(rp_6mo = Pneumonitis,
         have_rp_6mo = 1,
         time2 = 6) %>%
  select(nid, rp_6mo, have_rp_6mo, time2)


MAXPNEU = TEMP1 %>%
  filter(!is.na(Pneumonitis)) %>%
  group_by(nid) %>%
  summarise(max_fu_pneum = max(time_postrt),
            maxpneum = max(Pneumonitis)) %>%
  arrange(nid) %>%
  mutate(pneum3plus = ifelse(maxpneum >= 3, 1, 0))


EVENT = TEMP1 %>%
  filter(pneum2plus == 1)


ETIME = EVENT %>% 
  group_by(nid) %>%
  summarise(rp2_time = min(time_postrt)) %>%
  mutate(rp2 = 1)
  

NEVENT = TEMP1 %>%
  filter(pneum2plus == 0)


NETIME = NEVENT %>% 
  group_by(nid) %>%
  summarise(fu_time = max(time_postrt))


RP2 = ETIME %>%
  full_join(NETIME, by = 'nid') %>%
  mutate(rp2_time = ifelse(is.na(rp2_time), fu_time, rp2_time),
         rp2 = ifelse(is.na(rp2), 0, rp2)) %>%
  full_join(MAXPNEU, by = 'nid')

ALL = SUBGROUP %>%
  full_join(RP2, by = 'nid') %>%
  full_join(ESOPH2, by = 'nid') %>%
  full_join(QUARTILE, by = 'hosp') %>%
  filter(sub == 1) %>%
  mutate(esoph2plus = ifelse(maxesoph >= 2, 1, 0),
         esoph3plus = ifelse(maxesoph >= 3, 1, 0)) %>%
  rename(ecog = ecog_L6) %>%
  distinct(nid, .keep_all = T)


ALL2 = ALL %>%
  filter(meangy_ptv <= 200 | meangy_ptv >= 30 | meangy_lung <= 50) %>%
  mutate(RT = ifelse(TPN3Dnum > 0, 'IMRT', '3D'),
         ptv_2cm_spinalcord = ifelse(is.na(ptv_2cm_spinalcord), 'N', ptv_2cm_spinalcord),
         ptv_2cm_heart = ifelse(is.na(ptv_2cm_heart), 'N', ptv_2cm_heart),
         ptv_2cm_esoph = ifelse(is.na(ptv_2cm_esoph), 'N', ptv_2cm_esoph),
         ptv_2cm_brachialplexus = ifelse(is.na(ptv_2cm_brachialplexus), 'N', ptv_2cm_brachialplexus),
         ptv_2cm_other = ifelse(is.na(ptv_2cm_other), 'N', ptv_2cm_other),
         stage = ifelse(group_stage %in% c('IA', 'IB'), 1,
                        ifelse(group_stage %in% c('IIA', 'IIB'), 2, 
                               ifelse(group_stage %in% c('IIIA', 'IIIB'), 3, NA))),
         chemo_con = ifelse(IfYesChemo_Con_L5 == 'Y', 'Y', 'N'),
         chemo_con = ifelse(is.na(chemo_con), 'N', chemo_con),
         chemo_adj = ifelse(IfYesChemo_Adj_L5 == 'Y', 'Y', 'N'),
         chemo_adj = ifelse(is.na(chemo_adj), 'N', chemo_adj),
         chemo_neo = ifelse(IfYesChemo_Neo_L5 == 'Y', 'Y', 'N'),
         chemo_neo = ifelse(is.na(chemo_neo), 'N', chemo_neo),
         married_current = ifelse(marriage_1 == 'married', 'Y', 'N'),
         insur_type = ifelse(substr(insur_1,2,7)=='edicar', 'Medicare',
                             ifelse(insur_1 == 'none', 'None', 'Non-Medicare')),
         insur_secondary = ifelse(!is.na(insur_2), 'Y', 'N'),
         race2 = ifelse(race == 'white', 'White', 'AfAm/Other'),
         academic = ifelse(is.na(academic), 0, academic),
         CurrentSmoker_L4 = ifelse(CurrentSmoker_L4 == 'Unknown', NA, CurrentSmoker_L4),
         motion_mgmt = ifelse(MotionManagement_LRTD %in% c('Abdominal compression', 'Breath hold with device', 'Gating of radiotherapy (RPM, AlignRT, etc..)'), 'Breathhold','ITV/Other'),
         age = ifelse(age <= 0, NA, age)) %>%
  rename(ptv_vol = ptv_voltotal) %>%
  filter(stage %in% c(2,3)) %>%
  full_join(M1, by = 'nid') %>%
  full_join(M3, by = 'nid') %>%
  full_join(M6, by = 'nid') %>%
  filter(sub == 1) %>%
  mutate(have_rp_1mo = ifelse(is.na(have_rp_1mo), 0, have_rp_1mo),
         have_rp_3mo = ifelse(is.na(have_rp_3mo), 0, have_rp_3mo),
         have_rp_6mo = ifelse(is.na(have_rp_6mo), 0, have_rp_6mo)) %>%
  distinct(nid, .keep_all=T)
 
  
TEMP2 = L9 %>%
  select(nid, DiseaseStatus_L9, DateOfEval_l9) %>%
  mutate(DateOfEval_l9 = as.Date(DateOfEval_l9, format = '%d%b%Y')) %>%
  full_join(ALL2 %>%
              select(nid, rt_st_Date, sub), by='nid') %>%
  filter(sub == 1) %>%
  mutate(time = as.numeric(DateOfEval_l9 - rt_st_Date)/(365/12),
         timepoint2 = ifelse(time < 2, 'EOT', 
                             ifelse(time < 3.5, '1m',
                                    ifelse(time < 6, '3m', 
                                           ifelse(time < 12, '6m', NA)))))
  

### prog2 as in SAS file
PROG = TEMP2 %>%
  filter(DiseaseStatus_L9 %in% c('Local-Regional Progression', 'Distant Progression', 'Both Local/Regional Progression and Distant Progression')) %>%
  mutate(progression = 1,
         time_prog = time,
         lp = ifelse(DiseaseStatus_L9 %in% c('Local-Regional Progression',  'Both Local/Regional Progression and Distant Progression'), 1, 0),
         dp = ifelse(DiseaseStatus_L9 %in% c('Distant Progression',  'Both Local/Regional Progression and Distant Progression'), 1, 0)) %>%
  select(nid, timepoint2, time, progression, DiseaseStatus_L9, time_prog, lp, dp) %>%
  group_by(nid,timepoint2) %>%
  filter(row_number()==1) 
  

NOPROG = TEMP2 %>%
  filter(DiseaseStatus_L9 %in% c('No evidence of disease', 'No evidence of progression of disease')) %>%
  mutate(progression = 0,
         time_prog_fu = time) %>%
  select(nid, time_prog_fu) %>%
  distinct(nid, .keep_all = T)
  
  
PROG_TIME = PROG %>%
  full_join(NOPROG, by = 'nid') %>%
  mutate(time_prog = ifelse(is.na(time_prog), time_prog_fu, time_prog),
         progression = ifelse(!is.na(time_prog) & is.na(progression), 0, progression),
         lp = ifelse(!is.na(time_prog) & is.na(progression), 0, lp),
         dp = ifelse(!is.na(time_prog) & is.na(progression), 0, dp)) %>%
  distinct(nid, .keep_all=T) 


ALL2 = ALL2 %>%
  full_join(PROG_TIME, by = 'nid') %>%
  distinct(nid, .keep_all=T) %>%
  mutate(ptv_2cm_all = 0,
         ptv_2cm_all = ptv_2cm_all + ifelse(ptv_2cm_spinalcord == "Y" & !is.na(ptv_2cm_spinalcord), 1, 0),
         ptv_2cm_all = ptv_2cm_all + ifelse(ptv_2cm_heart == "Y" & !is.na(ptv_2cm_heart), 1, 0),
         ptv_2cm_all = ptv_2cm_all + ifelse(ptv_2cm_esoph == "Y" & !is.na(ptv_2cm_esoph), 1, 0),
         ptv_2cm_all = ptv_2cm_all + ifelse(ptv_2cm_brachialplexus == "Y" & !is.na(ptv_2cm_brachialplexus), 1, 0),
         ptv_2cm_all = ptv_2cm_all + ifelse(ptv_2cm_other == "Y" & !is.na(ptv_2cm_other), 1, 0),
         decile = factor(ntile(ptv_vol, 10)),
         ptv_vol_sq = ptv_vol ^ 2,
         row_num = row_number()
  ) 


