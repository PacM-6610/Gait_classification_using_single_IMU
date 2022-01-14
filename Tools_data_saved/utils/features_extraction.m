function [features, nb_chunk] = features_extraction(signal, try_number, window, overlap, offset, VERBOSE)
% Extract all the feature in the continuous and frequency domain from the signal.

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
elseif try_number >=28  && try_number <=39
    if mod(try_number,2) == 0 
        features_type = 'upSlope';
    else
        features_type = 'downSlope';
    end
else
    error('Invalid Try number')
end


% Init features extraction
if VERBOSE
    disp('Begin feature extraction ...');
end
accelerations=signal(try_number,(10:12));
angular_rate=signal(try_number,(16:18));
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

% Continuous time
[features_name, cross_temp] = add_feature(features_name, 'cross_correlation', nb_chunk, true);
[features_name, range_temp] = add_feature(features_name, 'range', nb_chunk);
[features_name, max_temp] = add_feature(features_name, 'max', nb_chunk);
[features_name, min_temp] = add_feature(features_name, 'min', nb_chunk);
[features_name, median_temp] = add_feature(features_name, 'median', nb_chunk);
[features_name, std_temp] = add_feature(features_name, 'std', nb_chunk);
[features_name, iqr_temp] = add_feature(features_name, 'iqr', nb_chunk);
[features_name, mean_temp] = add_feature(features_name, 'mean', nb_chunk);
[features_name, mean_harmonic_temp] = add_feature(features_name, 'mean_harmonic', nb_chunk);
[features_name, mean_absolute_value_temp] = add_feature(features_name, 'mean_absolute_value', nb_chunk);
[features_name, mad_temp] = add_feature(features_name, 'mad', nb_chunk);
[features_name, rms_temp] = add_feature(features_name, 'rms', nb_chunk);
[features_name, var_temp] = add_feature(features_name, 'var', nb_chunk);
[features_name, skewness_temp] = add_feature(features_name, 'skewness', nb_chunk);
[features_name, kurtosis_temp] = add_feature(features_name, 'kurtosis', nb_chunk);
[features_name, zero_crossing_temp] = add_feature(features_name, 'zero_crossing', nb_chunk, true);
[features_name, slope_change_temp] = add_feature(features_name, 'slope_change', nb_chunk);

% Frequency domain
[features_name, fft_specCentroid] = add_feature(features_name, 'specCentroid', nb_chunk, true);
[features_name, fft_specSpread] = add_feature(features_name, 'specSpread', nb_chunk, true);
%[features_name, fft_specEntropy] = add_feature(features_name, 'specEntropy', nb_chunk);

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
        range_temp(i,j)=range(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        range_temp(i,j+3)=range(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
      
        max_temp(i,j)=max(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        max_temp(i,j+3)=max(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        
        min_temp(i,j)=min(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        min_temp(i,j+3)=min(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        
        median_temp(i,j)=median(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        median_temp(i,j+3)=median(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        
        std_temp(i,j)=std(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        std_temp(i,j+3)=std(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        
        iqr_temp(i,j)=iqr(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        iqr_temp(i,j+3)=iqr(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        
        mean_temp(i,j)=mean(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        mean_temp(i,j+3)=mean(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        
        mean_harmonic_temp(i,j)=harmmean(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        mean_harmonic_temp(i,j+3)=harmmean(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        
        mean_absolute_value_temp(i,j)=mean(abs(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1))));
        mean_absolute_value_temp(i,j+3)=mean(abs(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1))));
    
        mad_temp(i,j)=mad(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        mad_temp(i,j+3)=mad(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        
        rms_temp(i,j)=rms(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        rms_temp(i,j+3)=rms(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        
        var_temp(i,j)=var(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        var_temp(i,j+3)=var(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        
        skewness_temp(i,j)=skewness(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        skewness_temp(i,j+3)=skewness(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        
        kurtosis_temp(i,j)=kurtosis(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        kurtosis_temp(i,j+3)=kurtosis(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        
        zero_crossing_temp(i,j)=zero_crossing_count(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        zero_crossing_temp(i,j+3)=zero_crossing_count(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        
        slope_change_temp(i,j)=slope_change_count(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        slope_change_temp(i,j+3)=slope_change_count(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)));
        
        
        % Frequency domain
        %%Time specifications:
        fs = 100;                      % samples per second
        N = window;

        %%Fourier Transform:
        acceleration_zm=accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1))-mean_temp(i,j);
        angular_rate_zm=angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1))-mean_temp(i,j+3);
        
        acceleration_zm_windowed = hann(window).*acceleration_zm;
        angular_zm_windowed = hann(window).*angular_rate_zm;
        
        accelerations_fft = abs(fft(acceleration_zm_windowed));
        accelerations_fft = accelerations_fft(1:N/2+1);
        
        angular_rate_fft = abs(fft(angular_zm_windowed));
        angular_rate_fft = angular_rate_fft(1:N/2+1);
        
        freq_bins = [0:N/2]*fs/N;
        
        %spectral centroid
        fft_specCentroid(i,j) = (freq_bins*accelerations_fft)/sum(accelerations_fft);
        fft_specCentroid(i,j+3) = (freq_bins*angular_rate_fft)/sum(angular_rate_fft);

        %spectral spread
        fft_specSpread(i,j) = (((freq_bins - fft_specCentroid(i,j)).^2 * accelerations_fft)/sum(accelerations_fft))^0.5;
        fft_specSpread(i,j+3) = (((freq_bins - fft_specCentroid(i,j+3)).^2 * angular_rate_fft)/sum(angular_rate_fft))^0.5;
   
        %spectral entropy
        %fft_specEntropy(i,j) = spectral_entropy(accelerations{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)),1/fs);
        %fft_specEntropy(i,j+3) = spectral_entropy(angular_rate{1,j}{1}(1+window*overlap*(i-1):window*overlap*(i+1)),1/fs);

    end
    
    % Norm
    acceleration_norm = sqrt(accelerations{1,1}{1}(1+window*overlap*(i-1):window*overlap*(i+1)).^2 + ...
        accelerations{1,2}{1}(1+window*overlap*(i-1):window*overlap*(i+1)).^2+ ...
        accelerations{1,3}{1}(1+window*overlap*(i-1):window*overlap*(i+1)).^2);
    angular_rate_norm = sqrt(angular_rate{1,1}{1}(1+window*overlap*(i-1):window*overlap*(i+1)).^2 + ...
        angular_rate{1,2}{1}(1+window*overlap*(i-1):window*overlap*(i+1)).^2+ ...
        angular_rate{1,3}{1}(1+window*overlap*(i-1):window*overlap*(i+1)).^2);
    
    range_temp(i,7)=range(acceleration_norm);
    max_temp(i,7)=max(acceleration_norm);
    min_temp(i,7)=min(acceleration_norm);
    median_temp(i,7)=median(acceleration_norm);
    std_temp(i,7)=std(acceleration_norm);
    iqr_temp(i,7)=iqr(acceleration_norm);
    mean_temp(i,7)=mean(acceleration_norm);
    mean_harmonic_temp(i,7)=harmmean(acceleration_norm);
    mean_absolute_value_temp(i,7)=mean(abs(acceleration_norm));
    mad_temp(i,7)=mad(acceleration_norm);
    rms_temp(i,7)=rms(acceleration_norm);
    var_temp(i,7)=var(acceleration_norm);  
    skewness_temp(i,7)=skewness(acceleration_norm);
    kurtosis_temp(i,7)=kurtosis(acceleration_norm);
    slope_change_temp(i,7)=slope_change_count(acceleration_norm);
    %fft_specEntropy(i,7)=spectral_entropy(acceleration_norm);
    
    
    range_temp(i,8)=range(angular_rate_norm);
    max_temp(i,8)=max(angular_rate_norm);
    min_temp(i,8)=min(angular_rate_norm);
    median_temp(i,8)=median(angular_rate_norm);
    std_temp(i,8)=std(angular_rate_norm);
    iqr_temp(i,8)=iqr(angular_rate_norm);
    mean_temp(i,8)=mean(angular_rate_norm);
    mean_harmonic_temp(i,8)=harmmean(angular_rate_norm);
    mean_absolute_value_temp(i,8)=mean(abs(angular_rate_norm));
    mad_temp(i,8)=mad(angular_rate_norm);
    rms_temp(i,8)=rms(angular_rate_norm);
    var_temp(i,8)=var(angular_rate_norm);  
    skewness_temp(i,8)=skewness(angular_rate_norm);
    kurtosis_temp(i,8)=kurtosis(angular_rate_norm);
    slope_change_temp(i,8)=slope_change_count(angular_rate_norm);
    %fft_specEntropy(i,8)=spectral_entropy(angular_rate_norm);
    
    features_type_temp(i) = features_type;
end

% Output
features = horzcat(cross_temp, range_temp, max_temp, min_temp, median_temp,...
    std_temp, iqr_temp, mean_temp, mean_harmonic_temp, mean_absolute_value_temp, ...
    mad_temp, rms_temp, var_temp, skewness_temp, kurtosis_temp, zero_crossing_temp, ...
    slope_change_temp, fft_specCentroid, fft_specSpread);
features = array2table(features, 'VariableNames', features_name);
features = addvars(features, features_type_temp, 'NewVariableNames', 'features_type');
if VERBOSE
    disp('End of feature extraction');
end

end
