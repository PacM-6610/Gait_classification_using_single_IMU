function [Entropy] = spectral_entropy(signal)
% This function clalculate the spectral entropy of a signal in the
% frequency domain.
P=sum(abs(fft(signal)).^2);

%Normalization
d=P(:);
d=d/sum(d+ 1e-12);

%Entropy Calculation
logd = log2(d + 1e-12);
Entropy = -sum(d.*logd)/log2(length(d));

end

