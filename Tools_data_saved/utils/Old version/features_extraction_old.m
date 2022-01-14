function [features, nb_chunk] = features_extraction(trunk, try_number, window, overlap, offset, VERBOSE)

% Input
if nargin <= 2
    window=100;
    overlap=0.5;
end

if nargin <= 4
    offset=200;
end

if nargin <= 5
    VERBOSE = false;
end

if nargin < 2
    error('Not enough input argument')
end

if overlap>0.75
    error('Overlap should not be greater than 0.75 (actual value: %.2f)', overlap);
end

if try_number >=16  && try_number <=27
    if mod(try_number,2) == 0 
        features_type = 'upStep';
    else
        features_type = 'downStep';
    end
elseif try_number >=4  && try_number <=9
    features_type = 'flat';
else
    error('Invalid Try number')
end


% Init features extraction
if VERBOSE
    disp('Begin feature extraction ...');
end
accelerations=trunk(try_number,(10:12));
angular_rate=trunk(try_number,(16:18));
for k=1:3
    accelerations{1,k}{1}=accelerations{1,k}{1}(offset+1:end);
    angular_rate{1,k}{1}=angular_rate{1,k}{1}(offset+1:end);
end
N=length(accelerations{1,1}{1});

% Preparation
modulus=mod(N, window*overlap);
nb_chunk=(N-modulus)/(window*overlap)-1;

% Features
features_name={};
[features_name, cross_temp] = add_feature(features_name, 'cross_correlation', nb_chunk, true);
[features_name, max_temp] = add_feature(features_name, 'max', nb_chunk);
[features_name, median_temp] = add_feature(features_name, 'median', nb_chunk);
[features_name, std_temp] = add_feature(features_name, 'std', nb_chunk);
[features_name, mean_temp] = add_feature(features_name, 'mean', nb_chunk);
[features_name, rms_temp] = add_feature(features_name, 'rms', nb_chunk);
[features_name, var_temp] = add_feature(features_name, 'var', nb_chunk);
[features_name, fft_specCentroid] = add_feature(features_name, 'specCentroid', nb_chunk, true);
[features_name, fft_specSpread] = add_feature(features_name, 'specSpread', nb_chunk, true);
features_type_temp=strings(nb_chunk,1);

% Extraction
for i=1:nb_chunk % nb_of_sample
    % Cross-correlation
    cross_temp(i,1)=cross_correlation(1,2,i,accelerations, window, overlap);
    cross_temp(i,2)=cross_correlation(1,3,i,accelerations, window, overlap);
    cross_temp(i,3)=cross_correlation(2,3,i,accelerations, window, overlap);
    cross_temp(i,4)=cross_correlation(1,2,i,angular_rate, window, overlap);
    cross_temp(i,5)=cross_correlation(1,3,i,angular_rate, window, overlap);
    cross_temp(i,6)=cross_correlation(2,3,i,angular_rate, window, overlap);

    
    for j=1:3 % nb_of_axis
        % Time domain
        max_temp(i,j)=max(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        max_temp(i,j+3)=max(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        
        
        median_temp(i,j)=median(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        median_temp(i,j+3)=median(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        
        std_temp(i,j)=std(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        std_temp(i,j+3)=std(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        
        mean_temp(i,j)=mean(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        mean_temp(i,j+3)=mean(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
    
        rms_temp(i,j)=rms(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        rms_temp(i,j+3)=rms(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        
        var_temp(i,j)=var(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        var_temp(i,j+3)=var(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        
        % Frequency domain
        %%Time specifications:
        fs = 100;                      % samples per second
        N = window;

        %%Fourier Transform:
        accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1))=accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1))-mean_temp(i,j);
        angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1))=angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1))-mean_temp(i,j+3);
        
        accelerations_fft = abs(fft(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1))));
        accelerations_fft = accelerations_fft(1:N/2+1);
        
        angular_rate_fft = abs(fft(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1))));
        angular_rate_fft = angular_rate_fft(1:N/2+1);
        
        freq_bins = [0:N/2]*fs/N;
        
        %spectral centroid
        fft_specCentroid(i,j) = (freq_bins*accelerations_fft)/sum(accelerations_fft);
        fft_specCentroid(i,j+3) = (freq_bins*angular_rate_fft)/sum(angular_rate_fft);

        %spectral spread
        fft_specSpread(i,j) = (((freq_bins - fft_specCentroid(i,j)).^2 * accelerations_fft)/sum(accelerations_fft))^0.5;
        fft_specSpread(i,j+3) = (((freq_bins - fft_specCentroid(i,j+3)).^2 * angular_rate_fft)/sum(angular_rate_fft))^0.5;
    end
    
    % Norm
    acceleration_norm = sqrt(accelerations{1,1}{1}(1+window*overlap*(i-1):window*overlap*(i+1)).^2 + ...
        accelerations{1,2}{1}(1+window*overlap*(i-1):window*overlap*(i+1)).^2+ ...
        accelerations{1,3}{1}(1+window*overlap*(i-1):window*overlap*(i+1)).^2);
    angular_rate_norm = sqrt(angular_rate{1,1}{1}(1+window*overlap*(i-1):window*overlap*(i+1)).^2 + ...
        angular_rate{1,2}{1}(1+window*overlap*(i-1):window*overlap*(i+1)).^2+ ...
        angular_rate{1,3}{1}(1+window*overlap*(i-1):window*overlap*(i+1)).^2);
    
    median_temp(i,7)=median(acceleration_norm);
    max_temp(i,7)=max(acceleration_norm);
    std_temp(i,7)=std(acceleration_norm);
    mean_temp(i,7)=mean(acceleration_norm);
    rms_temp(i,7)=rms(acceleration_norm);
    var_temp(i,7)=var(acceleration_norm);    
    
    median_temp(i,8)=median(angular_rate_norm);
    max_temp(i,8)=max(angular_rate_norm);
    std_temp(i,8)=std(angular_rate_norm);
    mean_temp(i,8)=mean(angular_rate_norm);
    rms_temp(i,8)=rms(angular_rate_norm);
    var_temp(i,8)=var(angular_rate_norm);  
    
    features_type_temp(i) = features_type;
end

% Output
features = horzcat(cross_temp, max_temp, median_temp, std_temp, mean_temp, rms_temp, var_temp, fft_specCentroid, fft_specSpread);
features = array2table(features, 'VariableNames', features_name);
features = addvars(features, features_type_temp, 'NewVariableNames', 'features_type');
if VERBOSE
    disp('End of feature extraction');
end

end
