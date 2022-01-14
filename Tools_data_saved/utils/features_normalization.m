function [features_normalized] = features_normalization(features, type)
% This function normalised the feature set. Two different normalization
% were implemented; MinMax and standardScaler. Choose the feature 
% normalization type with 'type'.
    
    if strcmp(type, 'minMax')
        nb_features=size(features,2)-1;
        features_normalized=features;
        for i = 1:nb_features
            features_temp=features{:,i};
            features_normalized{:,i} = 2.*(features_temp-min(features_temp))./(max(features_temp)-min(features_temp))-1 ;  % normalizing the data p
        end
    elseif strcmp(type, 'standardScaler')
        nb_features=size(features,2)-1;
        features_normalized=features;
        for i = 1:nb_features
            features_temp=features{:,i};
            mu = mean(features_temp);
            stddev = std(features_temp);
            features_normalized{:,i} = bsxfun(@minus,features_temp,mu); % subtract mean of each feature from original value
            features_normalized{:,i} = bsxfun(@rdivide,features_normalized{:,i},stddev); % divide by standard deviation
        end
    else
        disp("  'type' not valid, please choose minMax or standardScaler");
        features_normalized=features;
    end
    
end

