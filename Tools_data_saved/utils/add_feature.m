
function [features_name, feature_temp] = add_feature(features_name, name_new_features, nb_chunk ,bool_freq_domain)
% This function initialised the name of the features and initialised a zero matrix to 
% store the value calculated later.

if nargin <= 3 
        bool_freq_domain=false;
    end
    
    features_name_temp={};
    if bool_freq_domain 
        if strcmp(name_new_features, 'cross_correlation')
            features_name_temp{end+1}= append(name_new_features,'_acc_xy');
            features_name_temp{end+1}= append(name_new_features,'_acc_xz');
            features_name_temp{end+1}= append(name_new_features,'_acc_yz');
            features_name_temp{end+1}= append(name_new_features,'_rate_xy');
            features_name_temp{end+1}= append(name_new_features,'_rate_xz');
            features_name_temp{end+1}= append(name_new_features,'_rate_yz');
        else
            features_name_temp{end+1}= append(name_new_features,'_acc_x');
            features_name_temp{end+1}= append(name_new_features,'_acc_y');
            features_name_temp{end+1}= append(name_new_features,'_acc_z');
            features_name_temp{end+1}= append(name_new_features,'_rate_x');
            features_name_temp{end+1}= append(name_new_features,'_rate_y');
            features_name_temp{end+1}= append(name_new_features,'_rate_z');
        end
        feature_temp=zeros(nb_chunk,6);
    
    else
        features_name_temp{end+1}= append(name_new_features,'_acc_x');
        features_name_temp{end+1}= append(name_new_features,'_acc_y');
        features_name_temp{end+1}= append(name_new_features,'_acc_z');
        features_name_temp{end+1}= append(name_new_features,'_rate_x');
        features_name_temp{end+1}= append(name_new_features,'_rate_y');
        features_name_temp{end+1}= append(name_new_features,'_rate_z');
        features_name_temp{end+1}= append(name_new_features,'_norm_acc');
        features_name_temp{end+1}= append(name_new_features,'_norm_rate');

        feature_temp=zeros(nb_chunk,8);
    end
    features_name=[features_name, features_name_temp];
end



