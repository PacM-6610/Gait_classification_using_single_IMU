function [features_reduced] = pca(features_normalized, variance)
% The principle composant analysis (PCA) is a well-known dimensionality reduction technique
% which projects the data onto a new orthogonal basis of smaller dimensions
feat_norm=features_normalized{:,1:end-1};

[m, ~] = size(feat_norm); % m - no of egs; n - no of features

Sigma = (1/m)*(feat_norm')*(feat_norm); % Covariance matrix

[U, S, ~] = svd(Sigma); % Perform Singular Value Decomposition

traceS = trace(S); % Calculate sum of diagonal elements of S

for i=1:size(S,2)
    tempS = sum(diag(S(1:i,1:i))); % sum of K diagonal elements
    if ((tempS/traceS) >= (variance/100))
        break;
    end
end
K = i;

UReduced = U(:,1:K); % Consider only first K features
features_reduced = feat_norm*UReduced; % Get reduced data 

features_name=string(zeros(1,size(features_reduced,2)));
for i=1:size(features_reduced,2)
    features_name(i)= append(['feature_', num2str(i)]);
end
features_reduced=array2table(features_reduced, 'VariableNames', features_name);
temp=array2table(features_normalized.features_type, 'VariableNames', "features_type");
features_reduced=[features_reduced, temp];

end



