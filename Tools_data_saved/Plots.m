%% Plot results of a participant for stair up
% This regroups a lot of the different plot of the signal to explore more
% in details the result or the raw/preprocessed dataset.

 load ('../processed data/1');  
 
figure();

 for i=16:2:26
    accelerations=trunk(i,(10:12));
    p=i/2-7;
    subplot(3,2,p)
    title(['Try ',num2str(i)])
    hold on;
    xlabel('counter')
    ylabel('Acceleration in m/s^2')
    xlim([0 1400])
    ylim([-5 20])
    for j=1:3
        plot(accelerations{1,j}{1})
    end
 end
 sgtitle('Participant 1') 
 
 %% Plot 3 axis Observation for flat,upstairs,downstaire
close all
clc
load ('../processed data/1'); 
w_subplot=3;
h_subplot=3;
try_number_array=[4, 16, 17];
try_number_name = {'Flat', 'Cobble stone', 'Upsatairs', 'Downstairs', 'Slope Up', 'Slope Down', 'Bank left', 'Bank right', 'Grass'};
fig1 = figure;

for scenario_i=1:length(try_number_array)
    try_number=try_number_array(scenario_i);
    accelerations=trunk(try_number,(10:12));
    angular_rate = trunk(try_number,(16:18));
    for axes_i=1:3
        subplot(h_subplot, w_subplot, (scenario_i-1)*3+axes_i);
        plot(accelerations{1,axes_i}{1})
        if axes_i==1
            ylim([0 20])
            if scenario_i==1
                title('Vertical')
                ylabel('Flat','fontweight','bold')
            elseif scenario_i==2
                ylabel('Upstairs', 'fontweight','bold')
            elseif scenario_i==3
                ylabel('Downstairs','fontweight','bold')
            end
        
        elseif axes_i==2
            ylim([-3 3]) 
            if scenario_i==1
                title('Medio-lateral')
            end
        elseif axes_i==3
            ylim([-5 4])
            if scenario_i==1
                title('Anterior-posterior')
            end
        end
        yyaxis right
        plot(angular_rate{1,axes_i}{1}, ':');
        if axes_i==1
            ylim([-2 3])
        elseif axes_i==2
            ylim([-1 1])
        elseif axes_i==3
            ylim([-1 0.6])
        end
        xlim([100 600])
    end
end
han=axes(fig1,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Acceleration Amplitude (m/s^2)');
yyaxis right
han.YLabel.Visible='on';
ylabel(han,'Angular Velocity Amplitude (rad/s)');
xlabel(han,'Counter');

%% Plot results of a participant for different scenario
close all
clc
load ('../processed data/1'); 
w_subplot=3;
h_subplot=3;
try_number_array=[4, 10, 16, 17, 28, 29, 40, 41, 52];
try_number_name = {'Flat', 'Cobble stone', 'Upsatairs', 'Downstairs', 'Slope Up', 'Slope Down', 'Bank left', 'Bank right', 'Grass'};
fig = figure;

for i=1:length(try_number_array)
    
    subplot(h_subplot,w_subplot,i)
    try_number=try_number_array(i);
    accelerations=trunk(try_number,(10:12));
    angular_rate = trunk(try_number,(16:18));
    
    acceleration_amplitude=sqrt(accelerations{1,1}{1}.^2 + ...
    accelerations{1,2}{1}.^2 + ...
    accelerations{1,3}{1}.^2);

    angular_rate_amplitude=sqrt(angular_rate{1,1}{1}.^2 + ...
    angular_rate{1,2}{1}.^2 + ...
    angular_rate{1,3}{1}.^2);

    plot(acceleration_amplitude)
    xlim([100 600])
    ylim([0 20])
    yyaxis right
    plot(angular_rate_amplitude, ':');
    xlim([100 600])
    ylim([0 1.5])
    title(try_number_name{i})
end

han=axes(fig,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Resultant Acceleration Amplitude (m/s^2)');
yyaxis right
han.YLabel.Visible='on';
ylabel(han,'Resultant Angular Velocity Amplitude (rad/s)');
xlabel(han,'Counter');
 
%% Plot results of a participant for different sensor location
close all
clc
load ('../processed data/1'); 
w_subplot=2;
h_subplot=3;
try_number=4;
try_number_name = {'Flat'};
fig = figure;
data={trunk, wrist, shankL, shankR, thighL, thighR};
data_name = {'Trunk', 'Wrist', 'Shank left', 'Shank right', 'Thigh left', 'Thigh right'};
for i=1:length(data)
    
    subplot(h_subplot,w_subplot,i)
    accelerations = data{i}(try_number,(10:12));
    angular_rate = data{i}(try_number,(16:18));
    
    acceleration_amplitude=sqrt(accelerations{1,1}{1}.^2 + ...
    accelerations{1,2}{1}.^2 + ...
    accelerations{1,3}{1}.^2);

    angular_rate_amplitude=sqrt(angular_rate{1,1}{1}.^2 + ...
    angular_rate{1,2}{1}.^2 + ...
    angular_rate{1,3}{1}.^2);

    plot(acceleration_amplitude)
    xlim([100 600])
    ylim([0 23])
    yyaxis right
    plot(angular_rate_amplitude, ':');
    xlim([100 600])
    ylim([0 5])
    title(data_name{i})
end

han=axes(fig,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Resultant Acceleration Amplitude (m/s^2)');
yyaxis right
han.YLabel.Visible='on';
ylabel(han,'Resultant Angular Velocity Amplitude (rad/s)');
xlabel(han,'Counter');
 
 %% Plot resuls of different participant for specific try and sensor location
w_subplot=2;
h_subplot=3;
try_number=4;
fig = figure;

for i=1:6
    load (['../processed data/', num2str(i)]); 
    subplot(h_subplot,w_subplot,i)
    accelerations = trunk(try_number,(10:12));
    angular_rate = trunk(try_number,(16:18));
    
    acceleration_amplitude=sqrt(accelerations{1,1}{1}.^2 + ...
    accelerations{1,2}{1}.^2 + ...
    accelerations{1,3}{1}.^2);

    angular_rate_amplitude=sqrt(angular_rate{1,1}{1}.^2 + ...
    angular_rate{1,2}{1}.^2 + ...
    angular_rate{1,3}{1}.^2);

    plot(acceleration_amplitude)
    xlim([100 600])
    ylim([0 18])
    yyaxis right
    plot(angular_rate_amplitude, ':');
    xlim([100 600])
    ylim([0 2])
    title(['Participant ', num2str(i)]) 
end

han=axes(fig,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Resultant Acceleration Amplitude (m/s^2)');
yyaxis right
han.YLabel.Visible='on';
ylabel(han,'Resultant Angular Velocity Amplitude (rad/s)');
xlabel(han,'Counter');
 
 %% Plot results of different trials for a specific participant, sensor locationt and scenario
close all
clc
load ('../processed data/1'); 
w_subplot=2;
h_subplot=3;
fig = figure;

for try_number=16:2:26
    p=try_number/2-7;
    subplot(h_subplot,w_subplot,p)
    accelerations=trunk(try_number,(10:12));
    angular_rate = trunk(try_number,(16:18));
    
    acceleration_amplitude=sqrt(accelerations{1,1}{1}.^2 + ...
    accelerations{1,2}{1}.^2 + ...
    accelerations{1,3}{1}.^2);

    angular_rate_amplitude=sqrt(angular_rate{1,1}{1}.^2 + ...
    angular_rate{1,2}{1}.^2 + ...
    angular_rate{1,3}{1}.^2);

    plot(acceleration_amplitude)
    xlim([100 600])
    ylim([0 20])
    yyaxis right
    plot(angular_rate_amplitude, ':');
    xlim([100 600])
    ylim([0 1.5])
    title(['Try ',num2str(try_number)])
end

han=axes(fig,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Resultant Acceleration Amplitude (m/s^2)');
yyaxis right
han.YLabel.Visible='on';
ylabel(han,'Resultant Angular Velocity Amplitude (rad/s)');
xlabel(han,'Counter'); 
 
 %% Plot signal
load ('../processed data/1');  
try_number=4;
accelerations=trunk(try_number,(10:12));
angular_rate = trunk(try_number,(16:18));
acceleration_amplitude=sqrt(accelerations{1,1}{1}.^2 + ...
    accelerations{1,2}{1}.^2 + ...
    accelerations{1,3}{1}.^2);
angular_rate_amplitude=sqrt(angular_rate{1,1}{1}.^2 + ...
    angular_rate{1,2}{1}.^2 + ...
    angular_rate{1,3}{1}.^2);

fig = figure;
    plot(acceleration_amplitude)
    ylim([0 20])
    yyaxis right
    plot(angular_rate_amplitude, ':');
    ylim([0 1.5])
plot(acceleration_amplitude)
han=axes(fig,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Resultant Acceleration Amplitude (m/s^2)');
yyaxis right
han.YLabel.Visible='on';
ylabel(han,'Resultant Angular Velocity Amplitude (rad/s)');
xlabel(han,'Counter');
 
 %% Plot resuls of all the participant for stair up
 for k=1:30
    load (['../processed data/', num2str(k)]);  
    figure();

     for i=16:2:26
        accelerations=trunk(i,(10:12));
        p=i/2-7;
        subplot(3,2,p)
        title(['Try ',num2str(i)])
        hold on;
        xlabel('counter')
        ylabel('Acceleration in m/s^2')
        xlim([0 1400])
        ylim([-5 20])
        for j=1:3
            plot(accelerations{1,j}{1})
        end
     end
     sgtitle(['Participant ', num2str(k)]) 
 end
 
 
%% Plot the frequency domain of a participant
figure;   

for i=1:3
    %%Time specifications:
    Fs = 100;                      % samples per second
    dt = 1/Fs;                     % seconds per sample
    N = size(accelerations{1,i}{1},1);

    %%Fourier Transform:
    X = fftshift(fft(accelerations{1,i}{1}));

    %%Frequency specifications:
    dF = Fs/N;                      % hertz
    f = -Fs/2:dF:Fs/2-dF;           % hertz

    %%Plot the spectrum:
    subplot(3,1,i)
    nb_removed=1;
    X(end/2-nb_removed:end/2+nb_removed)=0;
    plot(f,abs(X)/N);
    ylabel('Amplitude')
    xlabel('Frequency (in hertz)');
    ylim([0 2])
    xlim([-50 50])
end

%% Plot the frequency domain of a participant
figure;   

for i=1:3
    %%Time specifications:
    Fs = 100;                      % samples per second
    dt = 1/Fs;                     % seconds per sample
    N = size(accelerations{1,i}{1},1);

    %%Fourier Transform:
    X = fftshift(fft(accelerations{1,i}{1}));

    %%Frequency specifications:
    dF = Fs/N;                      % hertz
    f = -Fs/2:dF:Fs/2-dF;           % hertz

    %%Plot the spectrum:
    subplot(3,1,i)
    nb_removed=1;
    X(end/2-nb_removed:end/2+nb_removed)=0;
    plot(f,abs(X)/N);
    ylabel('Amplitude')
    xlabel('Frequency (in hertz)');
    ylim([0 2])
    xlim([-50 50])
end

%% Plot the frequency domain 
 for k=1:1
    load (['../processed data/', num2str(k)]); 
    accelerations=trunk(16,(10:12)); 
    figure;   
    for i=1:3
        %%Time specifications:
        Fs = 100;                      % samples per second
        dt = 1/Fs;                     % seconds per sample
        N = size(accelerations{1,i}{1},1);
        
        %%Frequency specifications:
        dF = Fs/N;                      % hertz
        f = 0:dF:Fs/2-dF;           % hertz

        %%Fourier Transform:
        X = abs(fft(accelerations{1,i}{1}));
        X= X(1:size(f,2));


        %%Plot the spectrum
        subplot(3,1,i)
        plot(f,abs(X)/N);
        ylabel('Amplitude')
        xlabel('Frequency (in hertz)');
        ylim([0 2])
        xlim([-1 25])
    end
    sgtitle(['Participant ', num2str(k)]) 
 end
 