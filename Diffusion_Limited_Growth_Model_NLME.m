%Script to run mixed effects fitting for the diffusion-limited model on
%mock data

% Author: Adam Nasim,
%
% part of https://github.com/DrAdamNasim/Diffusion_Limited_Cancer_Growth_Model
% If using this or related code please cite 
% Nasim, A.; Yates, J.; Derks, G.; Dunlop, C. 
%     Mechanistic mathematical model of tumour growth and inhibition (diffusion-limited model)
%     (Manuscript submitted for publication).

%Call to mock growth data generator, as an example we simulate 10 subjects
%on time data for a typical preclinical xenograft experiment
Time_Data = [0 3 7 10 14 17 21];
No_Subjects = 10;
[Tdata, Volume, Subject,Group] = Mock_Growth_Data_Generator(Time_Data,No_Subjects);

%Call to the diffusion-limited model without treatment
nlme_model = @Diffusion_Limited_Growth_Function;
%Set initial parameters
phi0 = [-1 -1 -1 -1];

[PHI,PSI,stats,res] = nlmefitsa(Tdata,Volume,Subject,Group,nlme_model,phi0,'ParamTransform',[2 2 2 2],'LogLikMethod','lin','ComputeStdErrors',true);


