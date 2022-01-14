%% Feature reduction algorithm (Approach II)
% This is the 'main' which contains all the step applied for the second
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

% Remove 'Approach1' from path
rmv_folder=strcat(folder, '\Approach1');
rmpath(rmv_folder)
clc; clear all;
disp('All subfolder were added to the current global path' );

%% Import and plot resulted processed data - (Quick inspection of a specific signal)
load ('../processed data/1');  %choose candidate
accelerations=trunk(4,(10:12)); %choose imu and Try number

% Plot
figure();
hold on
for i=1:3
    plot(accelerations{1,i}{1}) 
end
xlabel('counter')
ylabel('Acceleration in m/s^2')
title('Pre-processed signal');

%% Features extraction
% Init
window=100; %640 or 100;
overlap=0.5;
nb_participant=30;
moreClass=false; % true or false
imuPosition='trunk'; % trunk, wrist, shankL, shankR, thighL, thighR
outliers=[18]; % list of outlier in the dataset

[features, participant_chunk_index] = create_dataset(window,overlap, nb_participant, imuPosition, moreClass, outliers);

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

%% Features normalization
features_normalized = features_normalization(features, 'standardScaler');
display('Features normalized');

%% Features reduction
explained_variance = 92; % Enter a number between 1 and 99 both inclusive based on
                         % how much variance to preserve. (Explained variance)
features_reduced = pca(features_normalized, explained_variance);
display('Features reduction done - the results is stocked in features_reduced');

%% Models - Classification
% Every model is derived and exported from MATLAB Classification Learner 
 % Endterm
%trainClassifier_approach2_w640                 % Baseline: 6,4s window, 32 features,  3 classes, trunk,  polynomial SVM, explained variance 92%
%trainClassifier_approach2_w100                 % Window:     1s window, 41 features,  3 classes, trunk,  quadratic SVM,  explained variance 92%
%trainClassifier_approach2_w640_MoreClasses     % Class:    6,4s window, 34 features,  5 classes, trunk,  quadratic SVM,  explained variance 92%
%trainClassifier_approach2_w640_wrist           % ImuPos:   6,4s window, 30 features,  3 classes, wrist,  quadratic SVM,  explained variance 92%
%trainClassifier_approach2_w640_shankR          % ImuPos:   6,4s window, 29 features,  3 classes, shankR, quadratic SVM,  explained variance 92%

%% Leave one out - Cross-Validation
model=@trainClassifier_approach2_w100; % chose model to test from above
[MCE, cmV_participants] = leave_one_out(model, features_reduced, nb_participant, participant_chunk_index, outliers);

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
if moreClass
    nb_class=5;
    precision_participants=zeros(nb_class, nb_participant);
    recall_participants=zeros(nb_class, nb_participant);
    f1score_participants=zeros(nb_class, nb_participant);
else
    precision_participants=zeros(3, nb_participant);
    recall_participants=zeros(3, nb_participant);
    f1score_participants=zeros(3, nb_participant);
end


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
if moreClass
    Class = {'UpSteps'; 'DownSteps'; 'DownSlope'; 'Flat'; 'UpSlope'};
else
    Class = {'Flat'; 'DownSteps'; 'UpSteps'};
end
figure;
confusionchart(cmV_tot, Class, 'DiagonalColor', [0.9290, 0.6940, 0.1250, 0.95]); % Plot confusion matrice

% Write table 'metric' (Recall - Precision - F1) in excel file
Precision = precision;
Recall = recall;
F1 = f1score;
metrics = table(Class,Precision,Recall,F1);

filename = 'Metrics.xlsx';
sheet_number=6; % Excel sheet number (For moreclass sheet number 3)
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

