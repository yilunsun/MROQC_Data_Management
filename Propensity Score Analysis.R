rm(list = ls())

setwd('~/Dropbox/GSRA/MROQC/');

source('data.R');

###########################################################################################
### Propensity Score
###########################################################################################
### fit model stage 1 

m1_ps = glm(RT=='IMRT' ~ age + race2 + married_current + CurrentSmoker_L4 + comorbid_count +
            group_stage + ptv_vol + ptv_vol_sq + academic + chemo_con + chemo_neo + chemo_adj +
            TumorStage_N_L4 + TumorStage_T_L4 + hist_collapsed + ptv_2cm_spinalcord +
            ptv_2cm_heart + ptv_2cm_esoph + ptv_2cm_brachialplexus + ptv_2cm_other + ptv_2cm_all + factor(decile),
            family = binomial(), data = ALL2)

summary(m1_ps)

prediction = data.frame(row_num = as.numeric(names(m1_ps$fitted.values)),
                        phat = m1_ps$fitted.values)

W = ALL2 %>%
  left_join(prediction, by= 'row_num') %>%
  mutate(wt = ifelse(RT == 'IMRT', 1/phat, 1/(1-phat))) %>%
  group_by(RT) %>%
  mutate(wt = wt * length(wt)/sum(wt, na.rm=T)) %>%
  ungroup()
  
### check common support
ggplot(data = W) +
  geom_histogram(aes(phat, fill = RT), bins = 30, alpha=0.5)

### Consortium Diagram
EXCLUDE = ALL2 %>%
  select(nid, meangy_ptv, RT, BID) %>%
  full_join(W %>%
              select(nid, wt), by = 'nid') %>%
  mutate(excl_rs = ifelse(is.na(meangy_ptv), "1. Dose missing",
                          ifelse(is.na(RT), "2. RT missing", 
                                 ifelse(is.na(wt), "3. Prop variable", NA))))

table(EXCLUDE$excl_rs)

W = W %>%
  filter(!is.na(wt))

  
### check covariate balance
library(flextable)
library(officer)
library(questionr)

### unweighted
num_of_digits = 0

cat_var = c('gender', 'race2', 'ecog', 'academic', 'chemo_neo', 'chemo_adj', 'chemo_con',
              'comorbid_count', 'TumorStage_T_L4', 'TumorStage_N_L4', 'group_stage', 
              'ptv_2cm_spinalcord', 'ptv_2cm_heart', 'ptv_2cm_esoph', 'ptv_2cm_brachialplexus',
              'ptv_2cm_other', 'ptv_2cm_all', 'hist_collapsed')

con_var = c('age', 'BMI_L4', 'ptv_vol', 'meangy_ptv', 'D95_ptv', 'meangy_lung', 'V5_lung',
            'V20_lung', 'meangy_esoph', 'meangy_heart', 'V5_heart', 'WeightLoss_L6',
            'PFTValues_2_L4', 'DLCValues_2_L4', 'D2CC_esoph')

############################################3######################3######################3######################3

cat('Unweighted:\n')

for (i in 1:length(cat_var))
{
  cat(paste0(cat_var[i],':\n\n'))
  eval(parse(text=paste0('print(round(prop.table(with(W, questionr::wtd.table(',cat_var[i], ', RT)), margin = 2) * 100, num_of_digits))')))
  cat('\n###########################\n\n')
}

for (i in 1:length(con_var))
{
  cat(paste0(con_var[i],':\n\n'))
 
  cat('3D:\t')
  eval(parse(text=paste0('temp_mean1=round(with(W[W$RT==\'3D\',], questionr::wtd.mean(', con_var[i], ')), 2)')))
  eval(parse(text=paste0('temp_sd1=round(sqrt(with(W[W$RT==\'3D\',], questionr::wtd.var(', con_var[i], '))), 2)')))
  cat(paste0(temp_mean1, ' (', temp_sd1, ')'))
  
  cat('\nIMRT:\t')
  eval(parse(text=paste0('temp_mean2=round(with(W[W$RT==\'IMRT\',], questionr::wtd.mean(', con_var[i], ')), 2)')))
  eval(parse(text=paste0('temp_sd2=round(sqrt(with(W[W$RT==\'IMRT\',], questionr::wtd.var(', con_var[i], '))), 2)')))
  cat(paste0(temp_mean2, ' (', temp_sd2, ')'))
  
  cat('\n###########################\n\n')
}

cat('Unweighted done.')

############################################3######################3######################3######################3
#weighted
############################################3######################3######################3######################3

cat('Weighted:\n')

for (i in 1:length(cat_var))
{
  cat(paste0(cat_var[i],':\n\n'))
  eval(parse(text=paste0('print(round(prop.table(with(W, questionr::wtd.table(',cat_var[i], ', RT, weights = wt)), margin = 2) * 100, num_of_digits))')))
  cat('\n###########################\n\n')
}

for (i in 1:length(con_var))
{
  cat(paste0(con_var[i],':\n\n'))
  
  cat('3D:\t')
  eval(parse(text=paste0('temp_mean1=round(with(W[W$RT==\'3D\',], questionr::wtd.mean(', con_var[i], ', weights = wt)), 2)')))
  eval(parse(text=paste0('temp_sd1=round(sqrt(with(W[W$RT==\'3D\',], questionr::wtd.var(', con_var[i], ', weights = wt))), 2)')))
  cat(paste0(temp_mean1, ' (', temp_sd1, ')'))
  
  cat('\nIMRT:\t')
  eval(parse(text=paste0('temp_mean2=round(with(W[W$RT==\'IMRT\',], questionr::wtd.mean(', con_var[i], ', weights = wt)), 2)')))
  eval(parse(text=paste0('temp_sd2=round(sqrt(with(W[W$RT==\'IMRT\',], questionr::wtd.var(', con_var[i], ', weights = wt))), 2)')))
  cat(paste0(temp_mean2, ' (', temp_sd2, ')'))
  
  cat('\n###########################\n\n')
}

cat('Weighted done.')

#########################################################################################################
### histogram and stageN
#########################################################################################################
ggplot(data = W) + 
  geom_histogram(aes(x = ptv_vol, color = RT), bins = 50, alpha=0.5)


W = W %>%
  mutate(stageN = ifelse(TumorStage_N_L4 != 'X', TumorStage_N_L4, NA))

#########################################################################################################
### ESOPHAGITIS ANALYSIS
#########################################################################################################

### ESoph plot
plotdata_unwt = W %>%
  select("RT", "esoph2plus") %>%
  filter(complete.cases(.)) %>%
  group_by(RT) %>%
  summarise(Prob = mean(esoph2plus) * 100,
            upper = Hmisc::binconf(sum(esoph2plus), n(), method="wilson")[3] * 100, 
            lower = Hmisc::binconf(sum(esoph2plus), n(), method="wilson")[2] * 100) %>%
  mutate(cat = 'Unadjusted')

plotdf = W %>%
  select("RT", "esoph2plus" , 'wt') %>%
  filter(complete.cases(.)) %>%
  mutate(weight = wt*n()/sum(wt)) %>%
  group_by(RT) %>%
  summarise(Prob = sum(esoph2plus*weight)/sum(weight) * 100, 
            upper = binconf(sum(esoph2plus*weight),sum(weight),method="wilson")[3] * 100, 
            lower = binconf(sum(esoph2plus*weight),sum(weight),method="wilson")[2] * 100) %>%
  mutate(cat = 'Adjusted') %>%
  union(plotdata_unwt)

### start plot unweighted
ggplot(data = plotdf, aes(x=factor(cat, levels = c('Unadjusted', 'Adjusted')), y=Prob, color=RT, group = RT)) +
  geom_point(position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin=lower, ymax=upper),width=.2,
                position=position_dodge(0.5)) +
  ylab("Probability of Gr 2 or Higher Esophagitis (%)") +
  xlab('') +
  ylim(c(0, 70)) +
  theme(axis.title.x = element_text(face="bold", colour="black",size=18),axis.text.x=element_text(size=20,face="bold",colour="black"))+
  theme(axis.title.y = element_text(face="bold", colour="black",size=18),axis.text.y = element_text(face="bold", colour="black",size=20))+
  theme(legend.title=element_blank())+
  theme(legend.text = element_text(size = 16, face = "bold"))+
  theme(legend.key.size=unit(1.2,"cm"))


### GLMM for Esoph2plus

proc glimmix data=w empirical method=laplace;
class rt(ref='3D') group_stage hosp  chemo_con chemo_adj chemo_neo ptv_2cm_esoph ptv_2cm_heart ptv_2cm_spinalcord ptv_2cm_brachialplexus;
model esoph2plus(event='1') = rt meangy_ptv group_stage chemo_con chemo_adj chemo_neo ptv_vol stageN ptv_2cm_esoph ptv_2cm_heart ptv_2cm_spinalcord ptv_2cm_brachialplexus
/ solution dist=binomial link=logit;
weight w;
random intercept/subject=hosp;
estimate 'OR(IMRT vs 3D)' rt 1 -1/ilink cl exp;
covtest glm;
where keep1mo=1;
ods output estimates=esoph_or;
ods output ParameterEstimates=parm_esoph;
run; title;


title 'Univairtae GLMM for Esoph2plus';
proc glimmix data=w empirical method=laplace;
class rt(ref='3D');
model esoph2plus(event='1') = rt / solution dist=binomial link=logit;
weight w;
random intercept/subject=hosp;
estimate 'OR(IMRT vs 3D)' rt 1 -1/ilink cl exp;
covtest glm;
where keep1mo=1;
ods output estimates=esoph_or;
ods output ParameterEstimates=parm_esoph;
run; title;

data parm_esoph;
set parm_esoph;
value = coalescec(of RT--chemo_neo);
or = exp(estimate);
or_lower = exp(estimate - 1.96*stderr);
or_upper = exp(estimate + 1.96*stderr);
run;

proc print data=parm_esoph noobs label;
var effect value or or_lower or_upper probt;
label or="OR" or_lower="95% Lower" or_upper="95% Upper" probt="P-value";
run;

%macro glm_uni_esoph(VAR, CLASS=);
title "GLMM for Esoph2plus - &VAR";
proc glimmix data=w empirical method=laplace;
class rt(ref='3D') &CLASS;
model esoph2plus(event='1') = &VAR / solution dist=binomial link=logit;
weight w;
random intercept/subject=hosp;
where keep1mo=1;
ods output ParameterEstimates=parm_esoph_&VAR;
run; title;

data parm_esoph_&VAR;
set parm_esoph_&VAR;
where effect ne "Intercept" and probt ne .;
length value $20;
value = &VAR;
or = exp(estimate);
or_lower = exp(estimate - 1.96*stderr);
or_upper = exp(estimate + 1.96*stderr);
keep effect value or: probt;
run;

proc print data=parm_esoph_&VAR noobs label;
var effect value or or_lower or_upper probt;
label or="OR" or_lower="95% Lower" or_upper="95% Upper" probt="P-value";
run;
%mend glm_uni_esoph;

%glm_uni_esoph(rt, CLASS=rt)
%glm_uni_esoph(meangy_ptv)
%glm_uni_esoph(ptv_2cm_esoph, CLASS=ptv_2cm_esoph)
%glm_uni_esoph(group_stage, CLASS=group_stage)
%glm_uni_esoph(chemo_con, CLASS=chemo_con)
%glm_uni_esoph(chemo_adj, CLASS=chemo_adj)
%glm_uni_esoph(chemo_neo, CLASS=chemo_neo)
%glm_uni_esoph(ptv_vol)
%glm_uni_esoph(stageN)

data parm_esoph_uni;
set parm_esoph_rt
parm_esoph_meangy_ptv
parm_esoph_ptv_2cm_esoph
parm_esoph_group_stage
parm_esoph_chemo_con
parm_esoph_chemo_adj
parm_esoph_chemo_neo
parm_esoph_ptv_vol
parm_esoph_stageN
;
run;

proc print data=parm_esoph_uni noobs;
var effect value or or_lower or_upper probt;
run;

#########################################################################################################
### PNEUMONITIS ANALYSIS
###########################################################################################################

* PNEUMONITIS ANALYSIS;
* create dataset that all Pneumonitis analyses will be based on;
* add weight according to number FU points (1,3,6);

proc sort data=temp1;
by pneum_time;
run;

proc freq data=temp1;
tables pneum2plus;
by pneum_time;
ods output onewayfreqs=pneum_time_prob;
where pneum_time ne .;
run;

proc sql noprint;
select percent into :PP1-:PP3
from pneum_time_prob
where pneum2plus = 1;
quit;

%put &PP1, &PP2, &PP3;

/* number of Pneumonitis is # of obs in w_rp */
data w_rp;
set w;
where max_fu_pneum ge 2 or maxpneum ge 2;  * check;
if have_rp_1mo=. then have_rp_1mo=0;
if have_rp_3mo=. then have_rp_3mo=0;
if have_rp_6mo=. then have_rp_6mo=0;
w2 = (have_rp_1mo*&PP1 + have_rp_3mo*&PP2 + have_rp_6mo*&PP3) / (&PP1 + &PP2 + &PP3);
if (rp_1mo=1 or rp_3mo=1 or rp_6mo=1) then w2=1;
if w2 in (. 0) then delete;
drop w;
run;

proc freq data=w;
where w ne .;
tables have_rp_1mo have_rp_3mo have_rp_6mo;
run;

%propensity(IN=w_rp, OUT=w_rp2, PNEUM_ONLY=Y);

data w_rp2;
set w_rp2;
w_rp2 = w * w2;
run;

ods rtf file='baseline char - pneum only.rtf';
%baseline_stats(w_rp2, Unweighted);
%baseline_stats(w_rp2, Weighted);
ods rtf close;

proc means data=w_rp2 n mean std min q1 median q3 max ndec=3;
var w_rp2;
run;

* FIRST, SHOW RAW DATA;

proc freq data=w_rp2;
tables have_rp_1mo*rt
have_rp_3mo*rt
have_rp_6mo*rt /norow nocum nopercent;
run;

* by timepoint;
data plot;
set m1(rename=rp_1mo=rp_bytime) m3(rename=rp_3mo=rp_bytime) m6(rename=rp_6mo=rp_bytime);
run;

proc sort data=w; by nid; run;

proc sort data=plot; by nid; run;

data plot2;
merge plot w;
by nid;
run;

data plot2;
set plot2;
where sub=1 and w ne .;
if rp_bytime ge 2 then rp2_bytime=1;
else if rp_bytime in (0,1) then rp2_bytime=0;
run;

ods graphics on;
goptions reset=all htext=1.75;
axis1 label=(angle=90 'Grade 2+ Pneumonitis') width=2 order=0 to 1 by .1 color=black;* repeat=1;
axis2 label= ('') width=2 offset=(10 pct 10pct);

proc gchart data=plot2;
vbar rt /sumvar=rp2_bytime group=time2 type=mean outside=mean space=2  errorbar=both raxis=axis1 ;
run; quit;

proc freq data=w;
tables rt*rp_1mo /nocol nocum nopercent;
tables rt*rp_3mo /nocol nocum nopercent;
tables rt*rp_6mo /nocol nocum nopercent;
run;

* make plot;
proc sort data=w_rp2; by rt; run;
*ods trace on;

proc freq data=w_rp2;
tables rp2 /nocol nocum nopercent binomial(score level='1');
by rt;
weight w2; * note, w2 is not the IPTW but the weight to accomodate varying FU;
ods output BinomialCLs=pneum_prob;
run;
*ods trace off;

proc freq data=w_rp2;
tables rp2 /nocol nocum nopercent binomial(score level='1');
by rt;
weight w_RP2; * note, w_RP2 is PRODUCT OF IPTW AND FU WEIGHT;
ods output BinomialCLs=pneum_prob_w;
run;

*include both results in plot;
data plot;
set pneum_prob(in=in1) pneum_prob_w;
if in1 then x = _N_;
else x = _N_ + 1;
y = 100 * proportion; output;
y = 100 * lowerCL; output;
y = 100 * upperCL; output;
keep rt x y;
run;

data anno;
length xsys ysys position $ 1
color function $ 8 text $ 50;
retain xsys '1' ysys '1' position '6'
color 'black' function 'label' text '+';
x=23; y=5; function='label'; color='black'; size=2.5; font='swiss'; text='Unadjusted'; output;
x=63; y=5; function='label'; color='black'; size=2.5; font='swiss'; text='Adjusted'; output;
run;

goptions reset=all htext=2.5;
axis1 label=(angle=90 'Grade 2+ Pneumonitis (%)') width=2 order=0 to 30 by 5;
axis2 label= ('') width=2 order=0 to 6 by 1 value=('' '3D  ' 'IMRT' '' '3D  ' 'IMRT' '') offset=(10 pct 10pct) minor=none;
symbol1 interpol=hiloc width=3 color=black; 

proc gplot data=plot;
plot y*x /vaxis=axis1 haxis=axis2 annotate=anno vref=5 to 25 by 5 cvref=ligr; 
run; quit;

proc freq data=w_rp2; tables rp2*hosp; run;

proc univariate data=w_rp2;
var rtend;
run;
/*min time is 19068*/
  
  data w_rp2;
set w_rp2;
time = rtend - 19068;
run;

title 'GLMM for RP2';
proc glimmix data=w_rp2 empirical;
class rt(ref='3D') group_stage hosp chemo_con chemo_adj chemo_neo ptv_2cm_esoph ptv_2cm_heart ptv_2cm_spinalcord ptv_2cm_brachialplexus;
model rp2(event='1') = rt meangy_ptv group_stage chemo_con chemo_adj chemo_neo ptv_vol ptv_2cm_esoph ptv_2cm_heart ptv_2cm_spinalcord ptv_2cm_brachialplexus
/*rt_st_yr rt_st_yr*rt*/
  / solution dist=binomial link=logit ;
weight w_rp2;
*weight w2;
random intercept / subject=hosp;
covtest glm;
where hosp not in (976 5175 13840);
estimate 'OR(IMRT vs 3D)' rt 1 -1/ilink cl exp;
*ods output estimates=pneum_or;
ods output ParameterEstimates=parm_pneum;
run; title;


title 'Univariate GLMM for RP2';
proc glimmix data=w_rp2 empirical;
class rt(ref='3D');
model rp2(event='1') = rt
/*rt_st_yr rt_st_yr*rt*/
  / solution dist=binomial link=logit ;
weight w_rp2;
*weight w2;
random intercept / subject=hosp;
covtest glm;
where hosp not in (976 5175 13840);
estimate 'OR(IMRT vs 3D)' rt 1 -1/ilink cl exp;
*ods output estimates=pneum_or;
ods output ParameterEstimates=parm_pneum;
run; title;


data parm_pneum;
set parm_pneum;
value = coalescec(of RT--chemo_neo);
or = exp(estimate);
or_lower = exp(estimate - 1.96*stderr);
or_upper = exp(estimate + 1.96*stderr);
run;

proc print data=parm_pneum noobs label;
var effect value or or_lower or_upper probt;
label or="OR" or_lower="95% Lower" or_upper="95% Upper" probt="P-value";
run;

%macro glm_uni_pneum(VAR, CLASS=);
title "GLMM for RP2 - &VAR";
proc glimmix data=w_rp2 empirical;
class rt(ref='3D') &CLASS;
model rp2(event='1') = &VAR
/ solution dist=binomial link=logit ;
weight w2;
random intercept / subject=hosp;
covtest glm;
where hosp not in (976 5175 13840);
ods output ParameterEstimates=parm_pneum_&VAR;
run; title;

data parm_pneum_&VAR;
set parm_pneum_&VAR;
where effect ne "Intercept" and probt ne .;
length value $20;
value = &VAR;
or = exp(estimate);
or_lower = exp(estimate - 1.96*stderr);
or_upper = exp(estimate + 1.96*stderr);
keep effect value or: probt;
run;

proc print data=parm_pneum_&VAR noobs label;
var effect value or or_lower or_upper probt;
label or="OR" or_lower="95% Lower" or_upper="95% Upper" probt="P-value";
run;
%mend glm_uni_pneum;

%glm_uni_pneum(rt, CLASS=rt)
%glm_uni_pneum(meangy_ptv)
%glm_uni_pneum(group_stage, CLASS=group_stage)
%glm_uni_pneum(chemo_con, CLASS=chemo_con)
%glm_uni_pneum(chemo_adj, CLASS=chemo_adj)
%glm_uni_pneum(chemo_neo, CLASS=chemo_neo)
%glm_uni_pneum(ptv_vol)

data parm_pneum_uni;
set parm_pneum_rt
parm_pneum_meangy_ptv
parm_pneum_group_stage
parm_pneum_chemo_con
parm_pneum_chemo_adj
parm_pneum_chemo_neo
parm_pneum_ptv_vol
;
run;

proc print data=parm_pneum_uni noobs;
var effect value or or_lower or_upper probt;
run;


%macro print_or(DATA, EXTRAVAR=);
proc print data=&DATA label noobs;
var &EXTRAVAR label ExpEstimate ExpLower ExpUpper Probt;
label ExpEstimate = "Odds Ratio"
ExpLower = "95% CI Lower"
ExpUpper = "95% CI Upper"
Probt = "P-Value";
run; title;
%mend print_or;

ods rtf file="Odds Ratios from GLMs.rtf" bodytitle;
title 'OR for Grade 2+ Esophagitis';
%print_or(esoph_or);
title 'OR for Grade 2+ Pneumonitis';
%print_or(pneum_or);
ods rtf close;

