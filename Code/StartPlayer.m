function [] = StartPlayer(PlayerHandle, autoReverse)
    if nargin <=1
        autoReverse = 3; % by default reverse and repeat
    end
    PlayerHandle.DataSource.Controls.AutoReverse=autoReverse;
    PlayerHandle.DataSource.Controls.Repeat=1;
    play(PlayerHandle.DataSource.Controls);
end