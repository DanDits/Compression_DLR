function [Frames, ranks] = DLR(StartA, frameCount, GetDelta, MakeFrame, approx, approxAsRank)
%DLR Executes the dynamic low rank algorithm using the given approximation
%level. Based on the paper "A projector-splitting integrator for dynamical
%low-rank approximation" by Christian Lubich and Ivan V. Oseledets in 2013
%published by Springer.
%INPUT:
%StartA: Start frame, a matrix defining the size of the whole problem.
%frameCount: Will perform frameCount - 1 DLR steps.
%GetDelta: Function that takes the current Frame and the index of
%   current time step starting from 1. Can calculate a difference to the 
%   next frame or evaluate some function on the current frame and returns a 
%   matrix with the same size as StartA.
%MakeFrame: The function MakeFrame can be the
%   default @Make_Frame or used to alter resulting Frame for better
%   visualization. Takes matrices U, S, V as arguments, default: Y=U*S*V'
%approx: An approximation level between 0.0 and 2.0 indicating the amount
%   of data saved. This is only an estimation. For compression the value
%   should be lower than 1.0. Regulates the rank used for approximation.
%   Higher values result in a higher rank.
%approxAsRank: If equal to 'asrank' the given approximation argument is interpreted
%   as the approximation rank to use. In this case approx should be an int
%   >=1. If ending with a '!' in any case, only one fix rank will be used,
%   else the rank can vary adaptively.
%OUTPUT:
%Frames: A tensor of size Frame dimension times framesCount holding the
%   approximations calculated by DLR and frames produced by MakeFrame.
%ranks: The approximation ranks used.

%% Setting up parameters for the Dynamic low rank approximation
frameDim = size(StartA);
normStartA = norm(StartA, 'fro');

if nargin >= 6 && strcmp(approxAsRank(1:min(end, 6)), 'asrank') && approx >= 1
    maxR = approx;
else
    % and r = approx/2 * sqrt(n*m), the factor 1/2 since we need to 'store' two
    % matrices U and V, we ignore the matrix of size r^2/2 for the approxfactor
    % we could perform a svd on this small matrix, push orthogonal side
    % matrices to U and V and only need to save r singular values
    % additionally
    maxR = approx/2*sqrt(prod(frameDim));
end
if nargin >= 6 && ~isempty(approxAsRank) && strcmp(approxAsRank(end), '!')
    fixedRank = true;
else
    fixedRank = false;
end

%fixed maximum approximation rank r >= 1, r <= n,m (dimensions of single frames)
maxR = floor(min(frameDim(2), min(frameDim(1), max(1, maxR))));
disp(strcat('Starting to calculate DLR with max rank=', num2str(maxR), '; true rank of starting matrix A=', num2str(rank(StartA))));

%% Starting value setup
ranks = zeros(frameCount, 1);
currR = maxR; 

%Given a mxn matrix StartA
[U,S,V] = GetRankApprox(StartA, currR);

%S is now diagonal, though it would only need to be invertible
Frames = zeros(frameDim(1), frameDim(2), frameCount);

%Now do the iteration and perform a step with each delta 
for i = 1:frameCount
    %Save to Frames
    ranks(i) = currR;
    Frames(:,:,i) = MakeFrame(U, S, V);
    
    %except for the last frame, perform the DLR step to the next frame
    if i < frameCount
        %GetDelta is a given function, can be a delta to the next frame or
        %a function that depends on the current frame and timestep
        DeltaA = GetDelta(Frames(:,:,i), i);
        [U,S,V] = DLR_Step(U, S, V, DeltaA);
        if ~fixedRank
            D = abs(diag(S)); % S is upper triangle matrix from default DLR_Step
            % If S does not have such a simple form we would need to check its SVD
            rankS  = nnz(D > eps * normStartA);
            if rankS < currR
                
                % we overapproximated! we can go smaller
                currR = rankS;
                [US, S, VS] = GetRankApprox(S, currR);
                U = U * US;
                V = V * VS;
            elseif rankS == currR && currR < maxR
                % we got full rank but can go bigger
                oldR = currR;
                currR = round(currR + (maxR - currR) / 2);
                currR = min(currR, maxR);
                UNew = zeros(frameDim(1), currR);UNew(:,1:oldR) = U;
                VNew = zeros(frameDim(2), currR);VNew(:,1:oldR) = V;
                SNew = zeros(currR);SNew(1:oldR,1:oldR) = S;
                U = UNew; V = VNew; S = SNew;
            end
        end
    end
end

end

function [U,S,V] = GetRankApprox(A, r)
% Calculate singular value decomposition and use spectral cutoff to get
% rank r

[U,S,V] = svd(A, 'econ');

%Sort to cut off greatest eigenvalues (absolute value) 
% and make matrix rank <= r, eigenvalues of svd are already sorted
%[S,SortInd]=sort(abs(diag(S)), 'descend');
%SortInd = SortInd(1:r);

S = diag(S);
SortInd = 1:r;

%Spectral cut off and set new smaller U,S,V
U = U(:,SortInd);
S = diag(S(SortInd));
V = V(:,SortInd);

end
