function [nb_zc] = zero_crossing_count(signal)
%This function calculate the zero crossing rate (ZCR)which indicate the number 
% of times the signal crosses zero and changes its signs. The ZCR is initialised 
%as zero and is incremented when the signal pass zero.
temp=find(signal(:).*circshift(signal(:), [-1 0]) <= 0);
nb_zc=size(temp,1);
end

