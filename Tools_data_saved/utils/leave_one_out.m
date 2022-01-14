function [MCE, cmV_tot] = leave_one_out(model, features_reduced, nb_participant, participant_chunk_index, outliers)
% Apply the leave-one-out crossvalidation. 

%% Dataset of features decomposition
% Init
participant_features_table=cell(nb_participant,1);
lower_index=1;
 for k=1:nb_participant
    if ismember(k, outliers)
        continue;
    end
    upper_index=lower_index+participant_chunk_index(k)-1;
    participant_features_table{k}=features_reduced(lower_index:upper_index,:);
    lower_index=lower_index+participant_chunk_index(k);
 end

%% Leave one out Cross-Validation

validationAccuracy=zeros(nb_participant,1);
MCE=zeros(nb_participant,1);
cmV_tot=cell(nb_participant,1);
for participant=1:nb_participant
    
    if ismember(participant, outliers)
        continue;
    end
    
    % Create training dataset leaving out only 'participant'
    trainingData=[];
    for k=1:nb_participant 
        if k~=participant
            trainingData = table_fusion(participant_features_table{k},trainingData);
        else
            testData = participant_features_table{k};
        end
        
    end
    [trainedClassifier, validationAccuracy(k)] = model(trainingData);
    
    prediction=trainedClassifier.predictFcn(testData);
    
    for pred = 1:size(prediction,1)
        if ~strcmp(prediction(pred), testData{pred, end})
            MCE(participant)=MCE(participant)+1;
        end
    end
    
    MCE(participant)=MCE(participant)/size(prediction,1);
    [cmV_tot{participant}, order] = confusionmat(prediction, testData{:, end});
    
end


end

