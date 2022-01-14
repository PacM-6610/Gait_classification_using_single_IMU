%% Feature selection algorithm (Approach I)
% This is the 'main' which contains all the step applied for the first
% method. For more information go see link: 
% https://github.com/PacM-6610/Gait_classification_using_single_IMU

clear all;
close all;
clc;

%% Create file 'processed data'
export_data('C:\Users\paco-\OneDrive\Bureau\Projet_de_semestre\DataBase\input_data_SD', 'on');

%% Add to path features function and all the models
% Determine the current path (It should be in the 'Tools_data_saved' folder)
folder = pwd;  
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

% Remove 'Approach2' from path
rmv_folder=strcat(folder, '\Approach2');
rmpath(rmv_folder)
clc;clear all;

display ('All subfolder were added to the current global path' );

%% Import and plot resulted processed data - (Quick inspection of a specific signal)
load ('../processed data/8');  %choose candidate
accelerations=trunk(8,(10:12)); %choose imu and Try number

% Plot
figure();
hold on
for i=1:3
    plot(accelerations{1,i}{1}) 
end
xlabel('counter')
ylabel('Acceleration in m/s^2')
title('Pre-processed signal');

 
%% Dataset of features creation 
% Init
window=640;
overlap=0.5;
nb_participant=30;
imuPosition='trunk'; % Trunk, wrist, shankL, shankR, thighL, thighR
outliers=[18]; % list of outlier in the dataset

[features, participant_chunk_index] = create_dataset(window, overlap, nb_participant, imuPosition, false, outliers);
%[features,participant_chunk_index] = create_dataset_old(window, overlap, nb_participant, imuPosition, outliers); % old version with less features (For midterm)

display('Features extracted');

%% Feature verification
% is there nan or inf
nb_nan=0;
nb_inf=0;
for f=1:size(features, 2)-1
    nb_nan=sum(isnan(features{:,f}), 'all')+nb_nan;
    nb_inf=sum(isinf(features{:,f}), 'all')+nb_inf;
    features{:,f}(isnan(features{:,f}))=0; % replace nan by 0
    features{:,f}(isinf(features{:,f}))=0; % replace inf by 0
end
disp(nb_nan) % display total number of nan found
disp(nb_inf) % display total number of inf found
%% Features normalization (Optionnal)
features = features_normalization(features, 'standardScaler');
display('Features normalized');

%% feature selection 
close all;
nb_features_selected = 32; % choose number of feature to keep
plot = false;
nb_feature_plot = 16; % Only the score and the result of the 'nb_feature_plot' first features with the best rank are ploted
[features_selected_table] = features_selection(features, nb_features_selected, plot, nb_feature_plot);
display('Features selection done - the results is stocked in features_selected_table');

 %% Models - Classification
% Every model is derived and exported from MATLAB Classification Learner 
 % MidTerm
%cubic_SVM_v1; % first try
%trainClassifier_Midterm; % extended window for first try 
%trainClassifier; % Second try (more features)

 % Report
%trainClassifier_approach1_w640;                % Baseline: 6,4s window, 13 features, 3 classes, trunk, gaussian SVM
%trainClassifier__approach1_w640_comparison;    % Baseline: 6,4s window, 32 features, 3 classes, trunk, cubic SVM

%% Leave one out - Cross-Validation
model=@trainClassifier__approach1_w640_comparison; % chose model to test from above
outliers=[18];
[MCE, cmV_participants] = leave_one_out(model, features_selected_table, nb_participant, participant_chunk_index, outliers);

% Final MCE
MCE_final=0;
for participant = 1:nb_participant
    if ismember(participant, outliers)
        continue;
    end
    MCE_final= MCE_final+MCE(participant);
end
MCE_final=MCE_final/(nb_participant-size(outliers, 2))*100;
fprintf('LOO :  %f \n',100-MCE_final);

%% Confusion matrix metrics
% For each participant
precision_participants=zeros(3, nb_participant);
recall_participants=zeros(3, nb_participant);
f1score_participants=zeros(3, nb_participant);

for participant = 1:nb_participant
    if ismember(participant, outliers)
        continue;
    end
    [precision_participants(:, participant), recall_participants(:, participant), ~, f1score_participants(:, participant)] = metrics_confusion_matrix(cmV_participants{participant});
end

% Calculate total mean per participant for each metrics
precision_participants_mean=mean(precision_participants)'; % mean of precision for all classes
recall_participants_mean=mean(recall_participants)';       % mean of recall for all classes
f1score_participants_mean=mean(f1score_participants)';     % mean of f1 for all classes

% All together
cmV_tot=combine_confusion_matrix(cmV_participants, outliers); % Aggregate confusion matrix of all participants
[precision, recall, accuracy, f1score]=metrics_confusion_matrix(cmV_tot);

display('Confusion matrix done');

%% Plot Results
% Confusion matrix and table(precision - recall - f1)
Class = {'Flat'; 'DownSteps'; 'UpSteps'}; figure;
confusionchart(cmV_tot, Class, 'DiagonalColor', [0.9290, 0.6940, 0.1250, 0.95]); % Plot confusion matrice

% Write table 'metric' (Recall - Precision - F1) in excel file
Precision = precision;
Recall = recall;
F1 = f1score;
metrics = table(Class,Precision,Recall,F1);

filename = 'Metrics.xlsx';
sheet_number=6; % Excel sheet number
writetable(metrics,filename,'Sheet', sheet_number,'Range','D2'); % Save Results in excel sheet

%% Plot Histogram
y=[precision_participants_mean, recall_participants_mean, f1score_participants_mean];
nb_metrics=size(y,2);
index=1:nb_participant;
for outlier=1:length(outliers)
    index=index(find(index~=outliers(outlier)));
end

y_without_outlier=y( index,:);
mean_y_wo=mean(y_without_outlier);
std_y_wo=std(y_without_outlier);

err=ones(nb_participant+1, size(y, 2)).*std_y_wo;
y=[y;mean_y_wo];

for outlier=1:length(outliers)
    err(outliers(outlier),:)=zeros(1,nb_metrics);
end

figure
bar(y)
x_tick=1:nb_participant;
x_tick=[string(x_tick), 'Mean'];
set(gca, 'XTickLabel', x_tick);

hold on

ngroups = size(y, 1);
nbars = size(y, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:,i), err(:,i), '.k');
end
legend('Precision', 'Recall', 'F1');
xticklabels(x_tick)
xticks(x)
xtickangle(90)
xlabel('Participants');
ylabel('Mean Error');
ylim([0 1.1]);
xlim([0, 32.1])

% Save graphe
name='Final_results_boxPlot';
saveas(gcf,[name, '.png'])

