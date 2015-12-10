function [Movie, Frames, map] = GenerateMovieOfGif(file)
    [Frames, map] = imread(file);
    Frames = im2double(Frames,'indexed');
    Movie = FramesToMovie(Frames, map);
    Frames = reshape(Frames, size(Frames,1),size(Frames,2),size(Frames,4)); % third single dimension not needed
end