%Script to run mixed effects fitting for the diffusion-limited model on the
%Novartis PDX data (open source)

% Author: Adam Nasim,
%
% part of https://github.com/DrAdamNasim/Diffusion_Limited_Cancer_Growth_Model
% If using this or related code please cite 
% Nasim, A.; Yates, J.; Derks, G.; Dunlop, C. 
%     Mechanistic mathematical model of tumour growth and inhibition (diffusion-limited model)
%     (Manuscript submitted for publication).

%Split data by cancer type
tdata_PDAC = xlsread('Gao_Untreated_Vol2.xlsm','PDAC','b2:b454');
ydata_PDAC = xlsread('Gao_Untreated_Vol2.xlsm','PDAC','c2:c454');
subject_PDAC = xlsread('Gao_Untreated_Vol2.xlsm','PDAC','d2:d454');
no_of_subjects_PDAC = subject_PDAC(end);
V_PDAC = linspace(1,no_of_subjects_PDAC,no_of_subjects_PDAC)';

tdata_GC = xlsread('Gao_Untreated_Vol2.xlsm','GC','b2:b481');
ydata_GC = xlsread('Gao_Untreated_Vol2.xlsm','GC','c2:c481');
subject_GC = xlsread('Gao_Untreated_Vol2.xlsm','GC','d2:d481');
no_of_subjects_GC = subject_GC(end);
V_GC = linspace(1,no_of_subjects_GC,no_of_subjects_GC)';

tdata_CRC = xlsread('Gao_Untreated_Vol2.xlsm','CRC','b2:b533');
ydata_CRC = xlsread('Gao_Untreated_Vol2.xlsm','CRC','c2:c533');
subject_CRC = xlsread('Gao_Untreated_Vol2.xlsm','CRC','d2:d533');
no_of_subjects_CRC = subject_CRC(end);
V_CRC = linspace(1,no_of_subjects_CRC,no_of_subjects_CRC)';

tdata_CM = xlsread('Gao_Untreated_Vol2.xlsm','CM','b2:b269');
ydata_CM = xlsread('Gao_Untreated_Vol2.xlsm','CM','c2:c269');
subject_CM = xlsread('Gao_Untreated_Vol2.xlsm','CM','d2:d269');
no_of_subjects_CM = subject_CM(end);
V_CM = linspace(1,no_of_subjects_CM,no_of_subjects_CM)';

tdata_BRCA = xlsread('Gao_Untreated_Vol2.xlsm','BRCA','b2:b490');
ydata_BRCA = xlsread('Gao_Untreated_Vol2.xlsm','BRCA','c2:c490');
subject_BRCA = xlsread('Gao_Untreated_Vol2.xlsm','BRCA','d2:d490');
no_of_subjects_BRCA = subject_BRCA(end);
V_BRCA = linspace(1,no_of_subjects_BRCA,no_of_subjects_BRCA)';

tdata_NSCLC = xlsread('Gao_Untreated_Vol2.xlsm','NSCLC','b2:b279');
ydata_NSCLC = xlsread('Gao_Untreated_Vol2.xlsm','NSCLC','c2:c279');
subject_NSCLC = xlsread('Gao_Untreated_Vol2.xlsm','NSCLC','d2:d279');
no_of_subjects_NSCLC = subject_NSCLC(end);
V_NSCLC = linspace(1,no_of_subjects_NSCLC,no_of_subjects_NSCLC)';

%Call to the diffusion-limited model without treatment
nlme_model = @Diffusion_Limited_Growth_Function;
%Set initial parameters
phi0 = [-1 -1 -1 -1];

[PHI_PDAC,PSI_PDAC,stats_PDAC,res_PDAC] = nlmefitsa(tdata_PDAC,ydata_PDAC,subject_PDAC,V_PDAC,nlme_model,phi0,'ParamTransform',[2 2 2 2],'LogLikMethod','lin','ComputeStdErrors',true);
[PHI_GC,PSI_GC,stats_GC,res_GC] = nlmefitsa(tdata_GC,ydata_GC,subject_GC,V_GC,nlme_model,phi0,'ParamTransform',[2 2 2 2],'LogLikMethod','lin','ComputeStdErrors',true);
[PHI_CRC,PSI_CRC,stats_CRC,res_CRC] = nlmefitsa(tdata_CRC,ydata_CRC,subject_CRC,V_CRC,nlme_model,phi0,'ParamTransform',[2 2 2 2],'LogLikMethod','lin','ComputeStdErrors',true);
[PHI_CM,PSI_CM,stats_CM,res_CM] = nlmefitsa(tdata_CM,ydata_CM,subject_CM,V_CM,nlme_model,phi0,'ParamTransform',[2 2 2 2],'LogLikMethod','lin','ComputeStdErrors',true);
[PHI_BRCA,PSI_BRCA,stats_BRCA,res_BRCA] = nlmefitsa(tdata_BRCA,ydata_BRCA,subject_BRCA,V_BRCA,nlme_model,phi0,'ParamTransform',[2 2 2 2],'LogLikMethod','lin','ComputeStdErrors',true);
[PHI_NSCLC,PSI_NSCLC,stats_NSCLC,res_NSCLC] = nlmefitsa(tdata_NSCLC,ydata_NSCLC,subject_NSCLC,V_NSCLC,nlme_model,phi0,'ParamTransform',[2 2 2 2],'LogLikMethod','lin','ComputeStdErrors',true);


