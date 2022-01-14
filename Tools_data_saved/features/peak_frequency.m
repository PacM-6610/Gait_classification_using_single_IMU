function [pF] = peak_frequency(signal, N, fs)
% Find the peak frequency of a signal
    [freq_response,freq_index] = freqz(signal,1,N,fs);  %N is the number of samples
    pM = max(abs(freq_response)); %magnitude
    pF = freq_index(abs(freq_response)==pM); %frequency
    plot(freq_index,10*log10(abs(freq_response)),'b',pF,10*log10(pM),'r*')
    title(sprintf('Peak Frequency = %f Hz (%.1f dB)',pF,10*log10(pM)));
    xlabel('Frequency [Hz]')
    ylabel('Power [dB]')
end

