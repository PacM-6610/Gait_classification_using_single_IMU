function [cross_corr_coef] = cross_correlation(a,b,i,matrice, window, overlap)
% This function calculate the cross-correlation coefficient which measures the linear 
% dependence between two axis. 
cross_corr_table=corrcoef(matrice{1,a}{1}(1+window*overlap*(i-1):window*overlap*(i+1)), matrice{1,b}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
cross_corr_coef=cross_corr_table(2,1);
end

