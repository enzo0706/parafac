% Post-test analysis for the ALS-MPPM algorithm of different configurations.
%
%  ----------    --      --    --------          --        --      --
% |          |  |  |    |  |  |        \        /  \      |  \    |  |
% |   -------   |  |    |  |  |   ---\  \      /    \     |   \   |  |
% |  |          |  |    |  |  |  |    \  \    /      \    |    \  |  |
% |   -------   |  |    |  |  |  |    |  |   /   /\   \   |     \ |  |
% |          |  |  |    |  |  |  |    |  |  |   /  \   |  |      \|  |
% |   -------   |  |    |  |  |  |    |  |  |   ----   |  |  |\      |
% |  |          |  |    |  |  |  |    |  |  |          |  |  | \     |
% |  |          |  |    |  |  |  |    /  /  |   ----   |  |  |  \    |
% |  |          \  \----/  /  |   ---/  /   |  |    |  |  |  |   \   |
% |  |           \        /   |        /    |  |    |  |  |  |    \  |
%  --             --------     --------      --      --    --      --
%  
% |    | |\    |  ---  |    |  -----  ----   ----   ---  ----- \    /
% |    | | \   |   |   |    | |      |    | |        |     |    \  /
% |    | |  \  |   |   \    / |----- |----   ----    |     |     \/
% |    | |   \ |   |    \  /  |      |  \        |   |     |     /
%  ----  |    \|  ---    \/    ----- |   \   ----   ---    |    /
% 
% ----------------------------------------------------------------------
% 
%  ----------    --      --    ----------   |           /\     ------
% |          |  |  \    /  |  |          |  |          /  \   |      \
% |   -------   |   \  /   |   ---    ---   |         /    \  |      |
% |  |          |    \/    |      |  |      |        /      \ |      /
% |   -------   |          |      |  |      |        |      | |------
% |          |  |  |\  /|  |      |  |      |        |------| |      \
% |   -------   |  | \/ |  |      |  |      |        |      | |      |
% |  |          |  |    |  |      |  |      |        |      | |      /
% |   -------   |  |    |  |   ---    ---    ------- |      |  ------
% |          |  |  |    |  |  |          |  
%  ----------    --      --    ----------   
%
% History:
%   2017.12.24 initial committed by Ming Shi
%
% Copyright (C) 2017 Ming Shi
% Copyright (C) 2017 EMI Lab, Dept. of E.E., Fudan University
%

clear all
close all
clc

[folder,~,~] = fileparts(which(mfilename));
load([folder,'/results/ALS-MPPM_COMPARISON.mat'])

SNR = 0:5:50;
N = 100;
Alg = {'L=0, M=1','L=1, M=2','L=2, M=3','L=1,M=4','L=3,M=4'};
Init = {'random','dtld'};
Trunc = 0.05;   % 5% truncated mean
N_trunc = N*Trunc;

for s = 1:length(SNR)
    for i = 1:length(Init)
        [~,idx] = sort(Err1(s,:,i),'descend');
        idx = idx(1:N_trunc);
        Res1(s,idx,i) = NaN;
        Iter1(s,idx,i) = NaN;
        Time1(s,idx,i) = NaN;
        Fail1(s,idx,i) = NaN;
        
        [~,idx] = sort(Err2(s,:,i),'descend');
        idx = idx(1:N_trunc);
        Res2(s,idx,i) = NaN;
        Iter2(s,idx,i) = NaN;
        Time2(s,idx,i) = NaN;
        Fail2(s,idx,i) = NaN;
        
        [~,idx] = sort(Err3(s,:,i),'descend');
        idx = idx(1:N_trunc);
        Res3(s,idx,i) = NaN;
        Iter3(s,idx,i) = NaN;
        Time3(s,idx,i) = NaN;
        Fail3(s,idx,i) = NaN;
        
        [~,idx] = sort(Err4(s,:,i),'descend');
        idx = idx(1:N_trunc);
        Res4(s,idx,i) = NaN;
        Iter4(s,idx,i) = NaN;
        Time4(s,idx,i) = NaN;
        Fail4(s,idx,i) = NaN;
        
        [~,idx] = sort(Err5(s,:,i),'descend');
        idx = idx(1:N_trunc);
        Res5(s,idx,i) = NaN;
        Iter5(s,idx,i) = NaN;
        Time5(s,idx,i) = NaN;
        Fail5(s,idx,i) = NaN;
    end
end

Iter = zeros(length(SNR),length(Alg),length(Init));
Err = zeros(length(SNR),length(Alg),length(Init));
Time = zeros(length(SNR),length(Alg),length(Init));
Fail = zeros(length(SNR),length(Alg),length(Init));
Acce_ratio = zeros(length(SNR),length(Alg),length(Init));

Iter(:,1,:) = mean(Iter1,2,'omitnan');
Err(:,1,:) = mean(Err1,2,'omitnan');
Time(:,1,:) = mean(Time1,2,'omitnan');
Fail(:,1,:) = mean(Fail1./Iter1,2,'omitnan');

Iter(:,2,:) = mean(Iter2,2,'omitnan');
Err(:,2,:) = mean(Err2,2,'omitnan');
Time(:,2,:) = mean(Time2,2,'omitnan');
Fail(:,2,:) = mean(Fail2./Iter2,2,'omitnan');

Iter(:,3,:) = mean(Iter3,2,'omitnan');
Err(:,3,:) = mean(Err3,2,'omitnan');
Time(:,3,:) = mean(Time3,2,'omitnan');
Fail(:,3,:) = mean(Fail3./Iter3,2,'omitnan');

Iter(:,4,:) = mean(Iter4,2,'omitnan');
Err(:,4,:) = mean(Err4,2,'omitnan');
Time(:,4,:) = mean(Time4,2,'omitnan');
Fail(:,4,:) = mean(Fail4./Iter4,2,'omitnan');

Iter(:,5,:) = mean(Iter5,2,'omitnan');
Err(:,5,:) = mean(Err5,2,'omitnan');
Time(:,5,:) = mean(Time5,2,'omitnan');
Fail(:,5,:) = mean(Fail5./Iter5,2,'omitnan');

Acce_ratio(:,1,:) = Iter(:,1,:)./Iter(:,1,:);
Acce_ratio(:,2,:) = Iter(:,2,:)./Iter(:,1,:);
Acce_ratio(:,3,:) = Iter(:,3,:)./Iter(:,1,:);
Acce_ratio(:,4,:) = Iter(:,4,:)./Iter(:,1,:);
Acce_ratio(:,5,:) = Iter(:,5,:)./Iter(:,1,:);
Acce_ratio = 1./Acce_ratio;

clear -regexp Iter[0-9] Err[0-9] Time[0-9] Fail[0-9]
save([folder,'/results/ALS-MPPM_COMPARISON_RESULT.mat'])