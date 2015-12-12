function Movie = FramesToMovie(Frames, map)
if (length(size(Frames))==3)
    % the one dimension in the third variable is required by the immovie
    % function
    Frames = reshape(Frames, size(Frames,1),size(Frames,2),1,size(Frames,3));
end
if nargin <= 1 || isempty(map)
    %map = colormap('jet'); % this would open an empty figure
    map = jet(size(get(0,'defaultfigurecolormap'),1));
end
Movie = immovie(Frames,map);

end

