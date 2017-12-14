% Test code for the ALS-MPPM algorithm of different configurations.
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

clear all
close all
clc

[folder,~,~] = fileparts(which(mfilename));

% Basic test conditions
D = [4,4,4];            % Dimensions of the tensor
Rt = 3;                 % Rank of the tensor
IterMax = 5000;         % Maximum number of iterations
Eps = 1e-8;             % Threshold of the normalized squared error
SNR = 0:5:50;           % Signal-to-noise ratio in dB
N = 2000;               % Number of random samples

% Additional conditions
Type = 'real';          % Data type 'real' or 'complex'
Init_vector = {'rand','dtld'};
R = 3;                  % Initialized rank
% UI control
h = waitbar(0,'Initializing...');
S = length(SNR);
I = length(Init_vector);
Step = 1/(N*S*I);

% Results pre-allocation
Err1 = zeros(S,N,I);    % ALS-MPPM (M=0, L=1)
Iter1 = zeros(S,N,I);
Time1 = zeros(S,N,I);
Fail1 = zeros(S,N,I);
Err2 = zeros(S,N,I);    % ALS-MPPM (M=1, L=2)
Iter2 = zeros(S,N,I);
Time2 = zeros(S,N,I);
Fail2 = zeros(S,N,I);
Err3 = zeros(S,N,I);    % ALS-MPPM (M=2, L=3)
Iter3 = zeros(S,N,I);
Time3 = zeros(S,N,I);
Fail3 = zeros(S,N,I);
Err4 = zeros(S,N,I);    % ALS-MPPM (M=1, L=4)
Iter4 = zeros(S,N,I);
Time4 = zeros(S,N,I);
Fail4 = zeros(S,N,I);
Err5 = zeros(S,N,I);    % ALS-MPPM (M=3, L=4)
Iter5 = zeros(S,N,I);
Time5 = zeros(S,N,I);
Fail5 = zeros(S,N,I);

% Test
for s = 1:S
    for n = 1:N
        %-------------------- Generating Data --------------------
        if strcmpi(Type,'real')
            At = rand(D(1),Rt);
            Bt = rand(D(2),Rt);
            Ct = rand(D(3),Rt);
        elseif strcmpi(Type,'complex')
            At = rand(D(1),Rt)+1i*rand(D(1),Rt);
            Bt = rand(D(2),Rt)+1i*rand(D(2),Rt);
            Ct = rand(D(3),Rt)+1i*rand(D(3),Rt);
        else
            error('Wrong data type! Must be ''real'' or ''complex''');
        end
        X1 = At*kat_rao(Ct,Bt).';
        X = reshape(X1,D(1),D(2),D(3));
        X = reshape(awgn(X(:),SNR(s),'measured'),D(1),D(2),D(3));
        X1 = reshape(X,D(1),D(2)*D(3));
        X2 = reshape(permute(X,[2,1,3]),D(2),D(1)*D(3));
        X3 = reshape(permute(X,[3,1,2]),D(3),D(1)*D(2));
        for i = 1:I
            waitbar(((s-1)*N*I+(n-1)*I+i)*Step,h,['SNR: ',num2str(SNR(s)),'dB, Tensor: ', num2str(n)])
            Init = Init_vector{i};
            %-------------------- Initialization --------------------
            if strcmpi(Init,'rand')     % Random
                A = randn(D(1),R);
                B = randn(D(2),R);
                C = randn(D(3),R);
            elseif strcmpi(Init,'eig')  % Eigenvectors
                if R<=D(1)
                    [A,~] = eigs(X1*X1',R,'lm');
                else
                    A = randn(D(1),R);
                    [A(:,1:D(1)),~] = eigs(Y1*Y1',D(1),'lm');
                end
                if R<=D(2)
                    [B,~] = eigs(X2*X2',R,'lm');
                else
                    B = randn(D(2),R);
                    [B(:,1:D(2)),~] = eigs(Y2*Y2',D(2),'lm');
                end
                if R<=D(3)
                    [C,~] = eigs(X3*X3',R,'lm');
                else
                    C = randn(D(3),R);
                    [C(:,1:D(3)),~] = eigs(Y3*Y3',D(3),'lm');
                end
            elseif strcmpi(Init,'svd')  % Singular vectors
                [U,~,~] = svd(X1);
                if R<=D(1)
                    A = U(:,1:R);
                else
                    A = randn(D(1),R);
                    A(:,1:D(1)) = U;
                end
                [U,~,~] = svd(X2);
                if R<=D(2)
                    B = U(:,1:R);
                else
                    B = randn(D(2),R);
                    B(:,1:D(2)) = U;
                end
                [U,~,~] = svd(X3);
                if R<=D(3)
                    C = U(:,1:R);
                else
                    C = randn(D(3),R);
                    C(:,1:D(3)) = U;
                end
            elseif strcmpi(Init,'dtld') % DTLD
                [A,B,C] = cp3_dtld(X,R);
            else
                error('Wrong initialization method!');
            end
            
            %-------------------- Running Test --------------------
            % ALS-MPPM (M=0, L=1)
            tic
            [~,~,~,iter,err,fail] = als_mppm_cpd3(X,R,Eps,IterMax,A,B,C,0,1);
            Time1(s,n,i) = toc;
            Iter1(s,n,i) = iter;
            Err1(s,n,i) = err(end);
            Fail1(s,n,i) = numel(find(fail == 1));
            % ALS-MPPM (M=1, L=2)
            tic
            [~,~,~,iter,err,fail] = als_mppm_cpd3(X,R,Eps,IterMax,A,B,C,1,2);
            Time2(s,n,i) = toc;
            Iter2(s,n,i) = iter;
            Err2(s,n,i) = err(end);
            Fail2(s,n,i) = numel(find(fail == 1));
            % ALS-MPPM (M=2, L=3)
            tic
            [~,~,~,iter,err,fail] = als_mppm_cpd3(X,R,Eps,IterMax,A,B,C,2,3);
            Time3(s,n,i) = toc;
            Iter3(s,n,i) = iter;
            Err3(s,n,i) = err(end);
            Fail3(s,n,i) = numel(find(fail == 1));
            % ALS-MPPM (M=1, L=4)
            tic
            [~,~,~,iter,err,fail] = als_mppm_cpd3(X,R,Eps,IterMax,A,B,C,1,4);
            Time4(s,n,i) = toc;
            Iter4(s,n,i) = iter;
            Err4(s,n,i) = err(end);
            Fail4(s,n,i) = numel(find(fail == 1));
            % ALS-MPPM (M=3, L=4)
            tic
            [~,~,~,iter,err,fail] = als_mppm_cpd3(X,R,Eps,IterMax,A,B,C,3,4);
            Time5(s,n,i) = toc;
            Iter5(s,n,i) = iter;
            Err5(s,n,i) = err(end);
            Fail5(s,n,i) = numel(find(fail == 1));
        end
    end
end
close(h)

%Save data
save([folder,'/results/ALS-MPPM_COMPARISON.mat'],'Iter*','Err*','Time*','Fail*','Type')