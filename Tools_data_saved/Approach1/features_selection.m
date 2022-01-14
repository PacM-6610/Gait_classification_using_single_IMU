function [features_selected_table] = features_selection(features, nb_features_selected, VERBOSE, nb_feature_plot)
% This function targets the selection of the best features to maximise the classification performance for a 
% given model. It first ranked every feature with a score (From three Tests; T-Test (TT), Chi Square Test (CST), 
% Minimum Redundancy Maximum Relevance Test (MRMRT)) and selected the right number to keep to optimise performance and computational cost. 
% This approach is categorised as a " filter method" since it selected the feature independently of the model.
%% feature selection T-Test
rng(8000,'twister');
obs=table2array(features(:,1:end-1));
grp=table2array(features(:, end));

% Divide data into a training set and a testing set
testing_percent=0.3;
holdoutCVP = cvpartition(grp,'holdout', round(testing_percent*length(obs)));

dataTrain = obs(holdoutCVP.training,:);
grpTrain = grp(holdoutCVP.training);

% Filter approach
dataTrainG1 = dataTrain(grp2idx(grpTrain)==1,:); % Flat
dataTrainG2 = dataTrain(grp2idx(grpTrain)==2,:); % Downsteps
dataTrainG3 = dataTrain(grp2idx(grpTrain)==3,:); % Upsteps

[~,p1,~,~] = ttest2(dataTrainG1,dataTrainG3,'Vartype','unequal'); % Upsteps-Flat
[~,p2,~,~] = ttest2(dataTrainG2,dataTrainG3,'Vartype','unequal'); % Upsteps-Downsteps
[~,p3,~,~] = ttest2(dataTrainG1,dataTrainG2,'Vartype','unequal'); % Flat-Downsteps

% Plot empirical cumulative distribution of the p-values
if VERBOSE==1
    figure;
    subplot(1,3,1);
    ecdf(p1);
    xlabel('P value'); 
    ylabel('CDF value')
    title('Upsteps and Flat');

    subplot(1,3,2);
    ecdf(p2);
    xlabel('P value'); 
    ylabel('CDF value')
    title('Upsteps and Downsteps');
    
    subplot(1,3,3);
    ecdf(p3);
    xlabel('P value'); 
    ylabel('CDF value')
    title('Flat and Downsteps');
    
end

 % Misclassification error MCE and resubstitution MCE for Upsteps-Flat
[~,idx_tt_uf] = sort(p1,2); % sort the features
 features_selected=idx_tt_uf(1:15);
 strrep(features.Properties.VariableNames(features_selected),'_','\_');
 
% Misclassification error MCE and resubstitution MCE for Upsteps-Downsteps
[~,idx_tt_ud] = sort(p2,2); % sort the features
 features_selected=idx_tt_ud(1:15);
 strrep(features.Properties.VariableNames(features_selected),'_','\_');
 
% Misclassification error MCE and resubstitution MCE for Flat-Downsteps
[~,idx_tt_fd] = sort(p3,2); % sort the features
 features_selected=idx_tt_fd(1:15);
 strrep(features.Properties.VariableNames(features_selected),'_','\_');
 
 
%% feature selection MRMR and FSC-Test
[idx_mrmr,scores_mrmrt] = fscmrmr(features,'features_type');
[idx_sct,score_sct] = fscchi2(features,'features_type');


% plot 
if VERBOSE==1
    figure;
    subplot(2,1,1);
    idx_mrmr_plot=idx_mrmr(1:nb_feature_plot);
    bar(scores_mrmrt(idx_mrmr_plot))
    xlabel('Predictor rank')
    ylabel('Predictor importance score')
    title('Minimum Redundancy Maximum Relevance')  
    set(gca, 'XTick', linspace(1,length(idx_mrmr_plot), length(idx_mrmr_plot)), 'XTickLabels', strrep(features.Properties.VariableNames(idx_mrmr_plot),'_','\_'));
    xtickangle(45)

    subplot(2,1,2);
    score_sct(isinf(score_sct))= nan;
    score_sct(isnan(score_sct))=max(score_sct)*1.3; % Replace Inf by an high value
    idx_sct_plot=idx_sct(1:nb_feature_plot);
    bar(score_sct(idx_sct_plot));
    xlabel('Predictor rank')
    ylabel('Predictor importance score')
    title('Chi-square test statistics')  
    set(gca, 'XTick', linspace(1,length(idx_sct_plot), length(idx_sct_plot)), 'XTickLabels', strrep(features.Properties.VariableNames(idx_sct_plot),'_','\_'));
    xtickangle(45)
end

%% Feature selection - Ranking and extraction of the chosen features
% Select best features
feature_score = zeros(1, length(idx_mrmr)); % initialisation
feature_score_temp = zeros(1, length(idx_mrmr)); % initialisation
threshold = find(scores_mrmrt(idx_mrmr) < 0.01 , 1); % Find the first score under 0.01
feature_score(idx_mrmr(1:threshold)) = threshold; % replace all score over 0.01 by a fixed weight
feature_score_temp(idx_sct) = linspace(length(idx_mrmr), 1, length(idx_mrmr))*threshold/(size(features, 2)-1); % Adjust second score
feature_score =  feature_score + feature_score_temp; % Weighted sum of our two scores
feature_score_temp(idx_tt_uf) = linspace(length(idx_mrmr), 1, length(idx_mrmr))*threshold/(size(features, 2)-1)/3; % Adjust third score
feature_score =  feature_score + feature_score_temp; % Weighted sum of all scores
feature_score_temp(idx_tt_ud) = linspace(length(idx_mrmr), 1, length(idx_mrmr))*threshold/(size(features, 2)-1)/3; % Adjust third score
feature_score =  feature_score + feature_score_temp; % Weighted sum of all scores
feature_score_temp(idx_tt_fd) = linspace(length(idx_mrmr), 1, length(idx_mrmr))*threshold/(size(features, 2)-1)/3; % Adjust third score
feature_score =  feature_score + feature_score_temp; % Weighted sum of all scores

[~,indexes] = sort(feature_score, 'descend');
features_selected_name = strrep(features.Properties.VariableNames(indexes(1:nb_features_selected)),'_','\_');
if VERBOSE==1
    disp(features_selected_name);
end
feature_score(indexes(1:nb_features_selected));
features_selected_table = features(:, indexes(1:nb_features_selected));
features_selected_table = [features_selected_table, features(:, end)];

if VERBOSE==1
    indexes_plot=indexes(1:nb_feature_plot);
    figure;
    bar(feature_score(indexes_plot));
    xlabel('Predictor rank')
    ylabel('Predictor importance score')
    title('Best features')  
    set(gca, 'XTick', linspace(1,length(indexes_plot), length(indexes_plot)), 'XTickLabels', strrep(features.Properties.VariableNames(indexes_plot),'_','\_'));
    xtickangle(45)


     % Misclassification error MCE and resubstitution MCE 

    testMCE = zeros(1,nb_feature_plot);
    resubMCE = zeros(1,nb_feature_plot);
    nfs = 1:1:nb_feature_plot;
    classf = @(xtrain,ytrain,xtest,ytest) ...
                 sum(~strcmp(ytest,classify(xtest,xtrain,ytrain,'quadratic'))); % if ERROR : The covariance matrix of each group in TRAINING must be positive definite.
                                                                                % Change 'quadratic' with 'diagQuadratic'


    for i = 1:nb_feature_plot
       fs = indexes(1:nfs(i));
       testMCE(i) = crossval(classf,obs(:,fs),grp,'partition',holdoutCVP)...
           /holdoutCVP.TestSize;
    end

    if VERBOSE==1
        figure;
        plot(nfs, testMCE,'-o');
        xlabel('Number of Features');
        ylabel('MCE');
        legend({'MCE on the test set'},'location','NW');
        title('Simple Filter Feature Selection Method');
    end

     min_MCE=min(testMCE);

     features_selected=indexes(1:nb_features_selected);
     strrep(features.Properties.VariableNames(features_selected),'_','\_');

    for i=1:nb_features_selected
        temp=strrep(features.Properties.VariableNames(features_selected(i)),'_',' ');
        if VERBOSE==1
            disp(temp{1});
        end
    end
end
end

