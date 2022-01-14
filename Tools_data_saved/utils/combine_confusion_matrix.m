function [cmV_tot] = combine_confusion_matrix(cmV_participants, outliers)
% This function aggregate the confusion matrix of all participant in one 
% unique confusion matrix.
nb_participant=length(cmV_participants);
nb_class=length(cmV_participants{1});

cmV_tot=zeros(nb_class);
for participant=1:nb_participant
    if ismember(participant, outliers)
        continue;
    end
    cmV_tot=cmV_tot+cmV_participants{participant};
end

end

