function Frames = DLR_Frames(FramesIn, approx, approxAsRank, plotError, plotRanks)
%DLR_FRAMES Performs the DLR algorithm on a list of frames using the given
%approximation level. Optionally plots the errors ('errors') and/or ranks
%('ranks' or 'ranks!' for big matrices) in separate figures.

%% Extract amount of frames and the dimension, giving it to the DLR routine
frameCount = size(FramesIn, 3);%assume >= 1
frameDim = [size(FramesIn, 1), size(FramesIn, 2)];

function DeltaA = Get_Delta(Current, i)
    DeltaA = FramesIn(:,:,i+1) - Current;
end
StartA = FramesIn(:,:,1); %Starting value
if nargin <= 2
    approxAsRank = '';
end
[Frames, approxRanks] = DLR(StartA, frameCount, @Get_Delta, @Make_Frame, approx, approxAsRank);

%% Only responsible for optionally plotting the error
if nargin >= 4 && strcmp(plotError, 'errors')
    errors = zeros(frameCount, 1);
    for i = 1:frameCount
       errors(i) = norm(Frames(:,:,i) - FramesIn(:,:,i), 'fro'); 
    end
    figure();
    plot(1:frameCount, errors, 'o');
    title('Error in Frobenius norm');
    legend('Error of approximation frame to real frame');
    xlabel('Frame number')
    ylabel('Error')
end

%% Only responsible for optionally creating a figure of involved ranks
if nargin >= 5 && strcmp(plotRanks(1:5), 'ranks') ...
        && (max(frameDim)<600 || strcmp(plotRanks, 'ranks!'))
    %if not too big and requested calculate and plot ranks of all frames
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

