function [ranks] = DLR(StartA, frameCount, Get_Delta, Make_Frame, varargin)
%DLR Executes the dynamic low rank algorithm using the given approximation
%level. Based on the paper "A projector-splitting integrator for dynamical
%low-rank approximation" by Christian Lubich and Ivan V. Oseledets in 2013
%published by Springer.
%INPUT:
%StartA: Start frame, a matrix defining the size of the whole problem.
%frameCount: Will perform frameCount - 1 DLR steps and invoke Make_Frame
%exactly frameCount times.
%Get_Delta: Function that takes the current Frame and the index of
%   current time step starting from 1. Can calculate a difference to the 
%   next frame or evaluate some function on the current frame and returns a 
%   matrix with the same size as StartA.
%Make_Frame: The function MakeFrame creates a frame out of the approximation and
%   alters it depending on the client requirements (like normalization). 
%   Takes matrices U, S, V as arguments, default: Y=U*S*V' and an index for
%   the frame to make starting from 1. Y can be interpreted and saved as
%   the current frame of interest. For indexed images you will want to
%   ensure every entry is >=1.
%approx: Optional numeric parameter. Defaults to 0.33. 
%   An approximation level between 0.0 and 2.0 indicating the amount
%   of data saved. This is only an estimation. For compression the value
%   should be lower than 1.0. Regulates the rank used for approximation.
%   Higher values result in a higher rank.
%rank: Optional numeric parameter. Only used if positive, overwrites set approx.
%fixed: Optional logic parameter. If true (default) then the approximation rank 
%   will stay fixed during the process, else it will be adapated to
%   previous calculations trying to make it as small as possible.
%OUTPUT:
%ranks: The approximation ranks used. If fixed is true this is a vector
%   with framesCount entries of the same value of the used rank.

%% Setting up parameters for the Dynamic low rank approximation
frameDim = size(StartA);
normStartA = norm(StartA, 'fro');


inp = inputParser;
addRequired(inp, 'StartA');
addRequired(inp, 'frameCount', @isnumeric);
addRequired(inp, 'Get_Delta', @ (f) isa(f, 'function_handle'));
addRequired(inp, 'Make_Frame', @ (f) isa(f, 'function_handle'));
addParameter(inp,'approx',0.33,@isnumeric);
addParameter(inp,'rank',0,@isnumeric);
addParameter(inp,'fixed',true,@islogical);
addParameter(inp,'showVTV',false,@islogical);

parse(inp,StartA,frameCount,Get_Delta,Make_Frame,varargin{:});
showVTV = inp.Results.showVTV;

% and r = approx/2 * sqrt(n*m), the factor 1/2 since we need to 'store' two
% matrices U and V, we ignore the matrix of size r^2/2 for the approxfactor
% we could perform a svd on this small matrix, push orthogonal side
% matrices to U and V and only need to save r singular values
% additionally
maxR = inp.Results.approx/2*sqrt(prod(frameDim));

if inp.Results.rank > 0
    %explicitly wants to use given rank for approximation
    maxR = inp.Results.rank;
end
fixedRank = inp.Results.fixed;

%fixed maximum approximation rank r >= 1, r <= n,m (dimensions of single frames)
maxR = floor(min(frameDim(2), min(frameDim(1), max(1, maxR))));
disp(strcat('Starting to calculate DLR with max rank=', num2str(maxR), '; true rank of starting matrix A=', num2str(rank(StartA))));

%% Starting value setup
ranks = zeros(frameCount, 1);
currR = maxR; 

%Given a mxn matrix StartA
[U,S,V] = Get_Rank_Approx(StartA, currR);
%S is now diagonal, though it would only need to be invertible

%Now do the iteration and perform a step with each delta
global FRAMESIN %DEBUGGING VTV
for i = 1:frameCount
    %Save to Frames
    ranks(i) = currR;
    CurrentFrame = Make_Frame(U, S, V, i);
    
    %FOR TESTING SOME PROPERTIES OF VVS FOR SOME FRAMES ONLY:
    if showVTV
        if i >= 1 && i <= 17

            aRank = rank(FRAMESIN(:,:,i+1));
            [~,~,AV] = Get_Rank_Approx(FRAMESIN(:,:,i+1), currR);
            VVS = AV'*V;

            disp(strcat('L=',num2str(i+1), '_: rankA=', num2str(aRank), ' norm vv_s=', num2str(norm(VVS)), ...
                ' cond vv_s=', num2str(cond(VVS)), ' rank vv_s=', num2str(rank(VVS)), ...
                ' dim vv_s=', num2str(size(VVS, 1))));

        end
    end

    %except for the last frame, perform the DLR step to the next frame
    if i < frameCount
        %GetDelta is a given function, can be a delta to the next frame or
        %a function that depends on the current frame and timestep
        DeltaA = Get_Delta(CurrentFrame, i);
        
        [U,S,V] = DLR_Step(U, S, V, DeltaA);
        
        %% If rank is not fixed we adaptively make it smaller or bigger
        if ~fixedRank
            D = abs(diag(S)); % S is upper triangle matrix from default DLR_Step
            % If S does not have such a simple form we would need to check its SVD
            rankS  = nnz(D > eps * normStartA);
            if rankS < currR
                
                % we overapproximated! we can go smaller
                currR = rankS;
                [US, S, VS] = Get_Rank_Approx(S, currR);
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
