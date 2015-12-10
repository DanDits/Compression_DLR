function Frame = Make_Frame(U, S, V)
%MAKE_FRAME Out of the decomposition Y = U*S*V' generate a valid frame
%   This implements the default behavior, that is returning Y and ensuring
%   that no entry is lower than 1 by setting those to exactly 1. This is
%   expected by indexed bitmaps. Values are not guaranteed to be integers.
   
Frame = U * S * V';
%normalize if value is out of bounds
%for indexed image all values need to be >=1
Frame(Frame < 1) = 1;

end

