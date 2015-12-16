function [Frames, FramesBestApprox] = DLR_Frames(FramesIn, varargin)
%DLR_FRAMES Performs the DLR algorithm on a list of frames using the given
%approximation level. 

%% Setup and parse parameters to Demo function

inp = inputParser;
addRequired(inp, 'FramesIn');
addParameter(inp,'approx',0.33,@isnumeric);
addParameter(inp,'rank',0,@isnumeric);
addParameter(inp,'fixed',true,@islogical);
addParameter(inp,'ploterror',false,@islogical);
addParameter(inp,'plotranks',false,@islogical);
addParameter(inp,'showbestapprox',false,@islogical);

parse(inp,FramesIn, varargin{:});
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

approxParam = '';
if approxAsRank
    approxParam = 'asrank';
end
if approxRankFixed
    approxParam = strcat(approxParam, '!');
end
Frame_Maker = @Make_Frame;
[Frames, approxRanks] = DLR(StartA, frameCount, @Get_Delta, Frame_Maker, approx, approxParam);

%% If requested setup best approximation in rank r manifold
FramesBestApprox = [];
if showBestApprox
   FramesBestApprox = zeros(size(Frames));
   for i = 1:frameCount
       [U,S,V] = Get_Rank_Approx(FramesIn(:,:,i), approxRanks(i));
       FramesBestApprox(:,:,i) = Frame_Maker(U, S, V);
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

