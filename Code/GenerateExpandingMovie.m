function [Movie, Frames] = GenerateExpandingMovie()

framesCount = 50; %amount of frames to calculate
frameDim = [256; 256]; %dimension of each frame

Frames = zeros(frameDim(1),frameDim(2),framesCount);

dx = 1 / frameDim(2);
dy = 1/ frameDim(1);
cycles = 1;

for i = 1:framesCount
    radius = (1 + sin(cycles * (i - 1) * 2 * pi / (framesCount - 1))) / 4;
    for m = 1: frameDim(1)
        y = dy * (m - 1);
        for n = 1: frameDim(2)
            x = dx * (n - 1);
            %dist = sqrt((x - 0.5)^2 + (y - 0.5)^2);%for expanding circle
            dist = max(abs(x - 0.5), abs(y - 0.5));
            if dist <= radius
                Frames(m, n, i) = 256 * (dist / radius);
            end
        end
    end
end

Frames(Frames < 1) = 1; % indexed image, set missing entries
Movie = FramesToMovie(Frames);
end

