%Function to generate mock growth data

% Author: Adam Nasim,
%
% part of https://github.com/DrAdamNasim/Diffusion_Limited_Cancer_Growth_Model
% If using this or related code please cite 
% Nasim, A.; Yates, J.; Derks, G.; Dunlop, C. 
%     Mechanistic mathematical model of tumour growth and inhibition (diffusion-limited model)
%     (Manuscript submitted for publication).
function [Tdata, Volume, Subject, Group] = Mock_Growth_Data_Generator(Time_Data,No_Subjects)

%Example inputs
%Time_Data = [0 3 7 10 14 17 21];
%No_Subjects = 2;

linspace(1,No_Subjects,No_Subjects)
Volume = [];
Tdata = [];
Subject = [];


Random_No = rand(length(Time_Data),No_Subjects);%Randomly generated from a uniform distribution between 0-1
Volume_Sort = sort(Random_No);%Sort in ascending order to mimick growth data
%Extract data for each mock subject and align in format for NLME fitting
for i = 1:No_Subjects
    Volume = [Volume; Volume_Sort(:,i)];
    Subject = [Subject; repmat(i,length(Volume_Sort(:,i)),1)];
end
Tdata = repmat(Time_Data',No_Subjects,1);
Group = linspace(1,No_Subjects,No_Subjects)';
end