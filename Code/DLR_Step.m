function [U, S, V] = DLR_Step(U, S, V, DeltaA)
%DLR_STEP Performs a DLR step starting from an approximation Y=U*S*V' using
% a DeltaA matrix that describes the current delta to the next
% approximation. Returns the next approximation in the series U*S*V'

%% Perform a Dynamic low rank step, using first order splitting method
%see 3.2 in paper, given Y = USV' and DeltaA 
% calculate same decomposition for next A, exact if rank A <= r

%Step1:
K = U * S + DeltaA * V;
[U, S] = qr(K, 0); %U*S = K, U orthonormal, S upper tridiagonal

%Step2:
S = S - U' * DeltaA * V;

%Step3:
L = V * S' + DeltaA' * U;
[V,S] = qr(L, 0); %V*S = L, V orthonormal, S upper tridiagonal
S = S';

end