function [A,B,C,Iter,Err,Fail] = als_mppm_cpd3(X,R,Eps,IterMax,A,B,C,M,L)
% ALS_MPPM_CPD3 PARAFAC decomposition using the ALS_MPPM algorithm.
%
% als_mppm_cpd3(X,R,Eps,IterMax,A,B,C,M,L)
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
% This function performs the CANDECOMP/PARAFAC Decomposition of a third-order 
% tensor (CPD3) using the Alternating Least Squares (ALS) algorithm with the 
% Matrix Polynomial Predictive Model (MPPM) acceleration.
%
% The CPD3 has the matricized form of
% 
%   X1 = A * kat_rao(C,B).'
%   X2 = B * kat_rao(C,A).'
%   X3 = C * kat_rao(B,A).'
%
% where X1, X2, X3 are the mode-1,2,3 matricizations of tensor X, and A, B, C 
% are the unknown loading/factor matrices to be estimated. 
%
% In MPPM, a matrix A is predicted by its previous values as
%
%   A(k) = sum(h(l) * A(k-l))
%
% where h(l) are the coefficients of the matrix-valued FIR predictive filter. 
% 
% The ALS-MPPM algorithm stops when any of the following convergence criteria is
% met:
%
%   - The normalized squared error drops below the preset threshold;
%   - The relative squared error drops below the preset threshold;
%   - The algorithm exceeds the maximum number of iterations.
%
% All inputs are explained as follows:
%   X       - I*J*K tensor
%   R       - scalar, rank of X
%   Eps     - scalar, threshold of the normalized squared error
%             (optional, default 1e-8)
%   IterMax - scalar, maximum number of the iterations
%             (optional, default 5000)
%   A       - I*R matrix, initial estimate of loading matrix A
%             (optional, default randomly initialized)
%   B       - J*R matrix, initial estimate of loading matrix B
%             (optional, default randomly initialized)
%   C       - K*R matrix, initial estimate of loading matrix C
%             (optional, default randomly initialized)
%   M       - scalar, degree of the assumed polynomial
%             (optional, default 0)
%   L       - scalar, length of the FIR predictive filter
%             (optional, default 1)
%
% All outputs are explained as follows:
%   A       - I*R final estimate of loading matrix A
%   B       - J*R final estimate of loading matrix B
%   C       - K*R final estiamte of loading matrix C 
%   Iter    - scalar, actual number of iterations taken before convergence
%   Err     - Iter*1 vector, normalized squared error of each iteration
%   Fail    - Iter*1 vector, indicator of prediction failure
%
% For detailed information of the ALS-MPPM algorithm, please refer to 
%
%   M. Shi, J. Zhang, B. Hu, B. Wang, and Q. Lu, “Convergence acceleration of 
% alternating least squares with a matrix polynomial predictive model for 
% PARAFAC decomposition of a tensor,” in Signal Processing Conference (EUSIPCO),
% 2017 25th European, pp. 1115–1119, IEEE, 2017.
%
% History:
%
%   2017.12.14 initial committed by Ming Shi
%
% Copyright (C) 2017 Ming Shi
% Copyright (C) 2017 EMI Lab, Dept. of E.E., Fudan University
%

    %% Check input parameters
    if nargin<2;    error('At least 2 parameters are needed.'); end
    if nargin<3;    Eps     = [];                               end
    if nargin<4;    IterMax = [];                               end
    if nargin<5;    A       = [];                               end
    if nargin<6;    B       = [];                               end
    if nargin<7;    C       = [];                               end
    if nargin<8;    M       = [];                               end
    if nargin<9;    L       = [];                               end
    if ndims(X)~=3; error('X must be a third-order tensor.');   end
    if L<=M;        error('L must be greater than M.');         end
    %% Default values
    if isempty(Eps);        Eps     = 1e-8;                     end
    if isempty(IterMax);    IterMax = 5000;                     end
    if isempty(M);          M       = 0;                        end
    if isempty(L);          L       = 1;                        end
    
    %% Loading matrices initialization
    [I,J,K] = size(X);
    if isempty(A)
        if isreal(X)
            A = randn(I,R);
        else
            A = randn(I,R)+1i*randn(I,R);
        end
    end
    if isempty(B)
        if isreal(X)
            B = randn(J,R);
        else
            B = randn(J,R)+1i*randn(J,R);
        end
    end
    if isempty(C)
        if isreal(X)
            C = randn(K,R);
        else
            C = randn(K,R)+1i*randn(K,R);
        end
    end
    
    %% Tracking initialization
    isPredictionOn = (L~=1 && M~=0);
    if isPredictionOn
        Ac = {A};
        Bc = {B};
        Cc = {C};
        for l = 2:L
            Ac = [{A},Ac];
            Bc = [{B},Bc];
            Cc = [{C},Cc];
        end
    end
    h0 = @(l) ones(1,l)./l;
    h1 = @(l) (4*l-6*(1:l)+2)/(l*(l-1));
    h2 = @(l) (9*l^2+(9-36*(1:l))*l+30*(1:l).^2-18*(1:l)+6)/(l*(l-1)*(l-2));
    h3 = @(l) (16*l^3+(24-120*(1:l))*l^2+(240*(1:l).^2-120*(1:l)+56)*l...
              -140*(1:l).^3+120*(1:l).^2-100*(1:l)+24)/(l*(l-1)*(l-2)*(l-3));
    switch(M)
        case 0
            h = h0(L);
        case 1
            h = h1(L);
        case 2
            h = h2(L);
        case 3
            h = h3(L);
        otherwise
            error('Not supported yet.')
    end
    Err = inf;
    Iter = 0;
    Fail = [];
    X1 = reshape(X,I,J*K);
    X2 = reshape(permute(X,[2,1,3]),J,I*K);
    X3 = reshape(permute(X,[3,1,2]),K,I*J);
    
    %% Tracking
    while Iter<IterMax && Err(end)>Eps
        % Prediction using MPPM
        if isPredictionOn
            Ap = zeros(I,R);
            Bp = zeros(J,R);
            Cp = zeros(K,R);
            for l = 1:L
                Ap = Ap+Ac{l}*h(l);
                Bp = Bp+Bc{l}*h(l);
                Cp = Cp+Cc{l}*h(l);
            end
            Fail = [Fail,0];
            % Handling the prediction failure
            if norm(X1-Ap*kat_rao(Cp,Bp).','fro')^2/norm(X1,'fro')^2>Err(end)
                Ap = A;
                Bp = B;
                Cp = C;
                Fail(end) = 1;
            end
        else
            Ap = A;
            Bp = B;
            Cp = C;
        end
        % ALS update
        A = X1*kat_rao(Cp,Bp)/((Cp'*Cp).*(Bp'*Bp));
        B = X2*kat_rao(Cp,A)/((Cp'*Cp).*(A'*A));
        C = X3*kat_rao(B,A)/((B'*B).*(A'*A));
        % Postprocess - Rescaling (Normalizing)
        DA = sqrt(sum(A.*conj(A)));
        DB = sqrt(sum(B.*conj(B)));
        DC = sqrt(sum(C.*conj(C)));
        D = diag((DA.*DB.*DC).^(1/3));
        A = A*diag(1./DA)*D;
        B = B*diag(1./DB)*D;
        C = C*diag(1./DC)*D;
        % Postprocess - Parameters update
        Iter = Iter+1;
        Err = [Err; norm(X1-A*kat_rao(C,B).','fro')^2/norm(X1,'fro')^2];
        if isPredictionOn
            Ac = [{A},Ac{1:end-1}];
            Bc = [{B},Bc{1:end-1}];
            Cc = [{C},Cc{1:end-1}];
        end
        if abs((Err(end)-Err(end-1))/Err(end)) < 1e-6
            break
        end
    end
    Err = Err(2:end);
end