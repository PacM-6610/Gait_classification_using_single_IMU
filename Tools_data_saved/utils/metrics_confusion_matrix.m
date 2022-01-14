function [precision, recall, accuracy, f1score] = metrics_confusion_matrix(confusion_matrix)
% Calculate the Precision, Recall and F1-measure from the
% 'confusion_matrix'.

nb_class=length(confusion_matrix);

% Precision - Recall - F1
precision=zeros(nb_class,1);
recall=zeros(nb_class,1);
f1score=zeros(nb_class,1);
for class=1:nb_class
    if sum(confusion_matrix(class,:))==0
        precision(class)=1;
    else
        precision(class)=confusion_matrix(class,class)/sum(confusion_matrix(class,:));
    end
    
    if sum(confusion_matrix(:, class))==0
        recall(class)=1;
    else
        recall(class)=confusion_matrix(class,class)/sum(confusion_matrix(:, class));
    end
    
    f1score(class)= 2*precision(class)*recall(class)/(precision(class)+recall(class));
end
accuracy=trace(confusion_matrix)/sum(confusion_matrix, 'all');
end

