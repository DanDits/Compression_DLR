function [Movie, Frames] = GenerateSimpleMovie()

framesCount = 40; %amount of frames to calculate
frameDim = [256; 128]; %dimension of each frame

Frames = zeros(frameDim(1),frameDim(2),framesCount);

%Generate simple figure that moves
figureDim = [floor(frameDim(1) / 8), floor(frameDim(2) / 8)];
Figure = zeros(figureDim);

for i = 1:figureDim(2)
    for j = 1:figureDim(1)
        Figure(j, i)=floor(1 + max(abs(i - figureDim(2) / 2), abs(j - figureDim(1) / 2)));
    end
end

%Generate frames
%Start and dest mark top left points of figure
figureStart = zeros(2,1);
figureStart(1) = frameDim(1) - figureDim(1) + 1;
figureStart(2) = 1;
figureDest = zeros(2,1);
figureDest(1) = 1; % end at first row
figureDest(2) = frameDim(2) - figureDim(2) + 1; % end so that it touches last column
for i = 1:framesCount
    %Draw figure in current Frame, let it move from bottom left to top
    %right, interpolate position linearly
    figureTopLeftRowColumn = (figureDest - figureStart) / ...
        (framesCount - 1) * (i - 1) + figureStart;
    row = floor(figureTopLeftRowColumn(1));
    column = floor(figureTopLeftRowColumn(2));
    Frames(row:row+figureDim(1)-1, column:column+figureDim(2)-1,i) = Figure; 
end
Frames(Frames < 1) = 1; % indexed image, set missing entries
Movie = FramesToMovie(Frames);
end