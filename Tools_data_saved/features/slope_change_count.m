function [nb_slope_change] = slope_change_count(signal)
% This function counts the number of times the slope of the signal changes its sign
nb_slope_change=sum(ischange(signal));
end

