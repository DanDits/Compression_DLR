function [Frames, FramesBestApprox] = DLR_Frames(FramesIn, varargin)
%DLR_FRAMES Performs the DLR algorithm on a list of frames using the given
%approximation level. 

%% Setup and parse parameters to Demo function

%DEBUG:
global FRAMESIN;
FRAMESIN=FramesIn;

inp = inputParser;
addRequired(inp, 'FramesIn');
addParameter(inp,'approx',0.33,@isnumeric);
addParameter(inp,'rank',0,@isnumeric);
addParameter(inp,'fixed',true,@islogical);
addParameter(inp,'ploterror',false,@islogical);
addParameter(inp,'plotranks',false,@islogical);
addParameter(inp,'showbestapprox',false,@islogical);
addParameter(inp,'showVTV',false,@islogical);

parse(inp,FramesIn, varargin{:});
showVTV = inp.Results.showVTV;
approx = inp.Results.approx;
approxAsRank = approx > 2;
if inp.Results.rank > 0
    %explicitly wants to use given rank for approximation
    approx = inp.Results.rank;
    approxAsRank = true;
end
approxRankFixed = inp.Results.fixed;
plotError = inp.Results.ploterror;
plotRanks = inp.Results.plotranks;
showBestApprox = inp.Results.showbestapprox;

%% Extract amount of frames and the dimension, giving it to the DLR routine
frameCount = size(FramesIn, 3);%assume >= 1

function DeltaA = Get_Delta(Current, i)
    DeltaA = FramesIn(:,:,i+1) - Current;
end
StartA = FramesIn(:,:,1); %Starting value

approxParam = 'approx';
if approxAsRank
    approxParam = 'rank';
end
Frames = zeros(size(StartA, 1), size(StartA, 2), frameCount);
function Frame = Make_Frame(U, S, V, frameIndex)
    %MAKE_FRAME Out of the decomposition Y = U*S*V' generate a valid frame
    %   This implements the default behavior, that is returning Y and ensuring
    %   that no entry is lower than 1 by setting those to exactly 1. This is
    %   expected by indexed bitmaps. Values are not guaranteed to be integers.

    Frame = U * S * V';
    %normalize if value is out of bounds
    %for indexed image all values need to be >=1
    Frame(Frame < 1) = 1;
    if frameIndex >= 1
        Frames(:,:,frameIndex) = Frame;
    end
end
approxRanks = DLR(StartA, frameCount, @Get_Delta, @Make_Frame, approxParam, approx, 'fixed', approxRankFixed, 'showVTV', showVTV);

%% If requested setup best approximation in rank r manifold
FramesBestApprox = [];
if showBestApprox
   FramesBestApprox = zeros(size(Frames));
   for i = 1:frameCount
       [U,S,V] = Get_Rank_Approx(FramesIn(:,:,i), approxRanks(i));
       FramesBestApprox(:,:,i) = Make_Frame(U, S, V, -1);
   end
end

%% Only responsible for optionally plotting the error
if plotError
    errors = zeros(frameCount, 1);
    for i = 1:frameCount
       errors(i) = norm(Frames(:,:,i) - FramesIn(:,:,i), 'fro'); 
    end
    figure();
    plot(1:frameCount, errors, 'o');
    legendText = 'Error of approximation frame to real frame';
    if showBestApprox 
        errorsBestApprox = zeros(frameCount, 1);
        for i = 1:frameCount
           errorsBestApprox(i) = norm(FramesBestApprox(:,:,i) - FramesIn(:,:,i), 'fro'); 
        end
        hold on
        plot(1:frameCount, errorsBestApprox, 'x');
        hold off
        legend(legendText, 'Error of rank r best approximation frame to real frame');
    else
        legend(legendText);
    end
    title('Error in Frobenius norm');
    xlabel('Frame number')
    ylabel('Error')
end

%% Only responsible for optionally creating a figure of involved ranks
if plotRanks
    %if requested calculate and plot ranks of all frames
    ranks = zeros(frameCount, 1);
    for i = 1:frameCount
       ranks(i) = rank(FramesIn(:,:,i)); 
    end
    figure();
    plot(1:frameCount, ranks, 'o', 1:frameCount, approxRanks);
    title('Involved ranks');
    legend('Rank of given frame', 'Used approximation rank r');
    xlabel('Frame number')
    ylabel('Rank')
end
end

