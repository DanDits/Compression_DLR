function [] = Animated_Surf(Frames, X, Y, tau, slowDownFactor)
%% Animate frames in a surface plot
s = surf(X, Y, Frames(:,:,1),'EdgeColor','none','LineStyle','none');
colormap('hot')
colorbar
totalMin = min(min(min(Frames)));
totalMax =  max(totalMin + 1, max(max(max(Frames))));
caxis([totalMin totalMax] ) ; 
mode manual;
zlim([totalMin totalMax]);
view(45,45);
T = tau * (size(Frames,3) - 1);
for j = 2:size(Frames,3)
   s.ZData = Frames(:,:,j); 
   title(strcat('Time passed:  ', num2str((j - 1) * tau, '%.1f'), '/', num2str(T, '%.1f')));
   pause(tau * slowDownFactor);
end

end

