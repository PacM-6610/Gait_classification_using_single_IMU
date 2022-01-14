function [features, participant_chunk_index] = create_dataset(window,overlap, nb_participant, imuPosition, moreClass, outliers)
% This function create the dataset use for training and testing. It extract
% the feature of the signal for a particular 'imuposition'. By default it
% extract the feature of the three class 'Upsetp', 'Flat' and 'Downstep'
% but 'Slope up and 'slope down' can be added with 'Moreclass' equal to
% true. The participant in the matrix 'ouliers' are excluded.
% 'participant_chunk_index' contains the number of sample for each
% participant to easily sepearte the dataset between participants.
%% Dataset of features creation
participant_chunk_index=zeros(nb_participant,1);
 for k=1:nb_participant
    if ismember(k, outliers)
        continue;
    end
     
    load(['../processed data/', num2str(k)]); 
    if strcmp(imuPosition, 'trunk')
        dataset=trunk;  
    elseif strcmp(imuPosition, 'wrist')
        dataset=wrist;
        
    elseif strcmp(imuPosition, 'shankL')
        dataset=shankL;
        
    elseif strcmp(imuPosition, 'shankR')
        dataset=shankR;
        
    elseif strcmp(imuPosition, 'thighL')
        dataset=thighL;
        
    elseif strcmp(imuPosition, 'thighR')
        dataset=thighR;
        
    else
        error('Invalid IMU position')
    end
    
    
    % Extract
    % Upstep and downstep
    try_number=16;
    if k==1
        [features, nb_chunk] = features_extraction(dataset, try_number, window, overlap);
        participant_chunk_index(k)=participant_chunk_index(k)+nb_chunk;
    else
        [features_temp, nb_chunk] = features_extraction(dataset, try_number, window, overlap);
        [features] = table_fusion(features, features_temp);
        participant_chunk_index(k)=participant_chunk_index(k)+nb_chunk;
    end
    for try_number =17:27
        [features_temp, nb_chunk] = features_extraction(dataset, try_number, window, overlap);
        [features] = table_fusion(features, features_temp);
        participant_chunk_index(k)=participant_chunk_index(k)+nb_chunk;
    end

    % Flat
    for try_number = 4:9
        [features_temp, nb_chunk] = features_extraction(dataset, try_number, window, overlap);
        [features] = table_fusion(features, features_temp);
        participant_chunk_index(k)=participant_chunk_index(k)+nb_chunk;
    end
    
    % UpSlop and DownSlope
    if moreClass    
        for try_number = 28:39
            [features_temp, nb_chunk] = features_extraction(dataset, try_number, window, overlap);
            [features] = table_fusion(features, features_temp);
            participant_chunk_index(k)=participant_chunk_index(k)+nb_chunk;
        end
    end
 end
 
end



