function [features,participant_chunk_index] = create_dataset_old(window, overlap, nb_participant, imuPosition, outliers)
%create_dataset_old is the older version of create_dataset which create the
%dataset with features_extraction_old which contains less feature than it
%is newer version.

participant_chunk_index=zeros(nb_participant,1);
 for k=1:nb_participant
    if ismember(k, outliers)
        continue;
    end
    load (['../processed data/', num2str(k)]); 
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
        [features, nb_chunk] = features_extraction_old(dataset, try_number, window, overlap);
        participant_chunk_index(k)=participant_chunk_index(k)+nb_chunk;
    else
        [features_temp, nb_chunk] = features_extraction_old(dataset, try_number, window, overlap);
        [features] = table_fusion(features, features_temp);
        participant_chunk_index(k)=participant_chunk_index(k)+nb_chunk;
    end
    for try_number =17:27
        [features_temp, nb_chunk] = features_extraction_old(dataset, try_number, window, overlap);
        [features] = table_fusion(features, features_temp);
        participant_chunk_index(k)=participant_chunk_index(k)+nb_chunk;
    end

    % Flat
    for try_number = 4:9
        [features_temp, nb_chunk] = features_extraction_old(dataset, try_number, window, overlap);
        [features] = table_fusion(features, features_temp);
        participant_chunk_index(k)=participant_chunk_index(k)+nb_chunk;
    end
 end
end

