%% Wasserstein Modulation Index (wMI)
% This code ("wMI_CallerRoutine.m") demonstrates the basic usage of the wMI (Ohki, 2022, J Neurosci Methods, DOI: 10.1016/j.jneumeth.2022.109578). 

% You will need the following files in the current directory or Matlab path:
%  eegfilt.m
%  Wasserstein_MI.m
%  Questions -> takefumi2ohki@gmail.com


clear all, close all

% Constructing an Example
% parameters for simulation data
srate=1000;
data_length = 10000; %( i.e., 10 sec)
t=1:1:data_length;
tVec =1/srate:1:data_length;

% parameter for the coupling strength
nonmodulatedamplitude = .3; % increase this parameters to get less modulation; you'll see that this is reflected in the wMI value

% lower frex (i.e., nesting phase frex) and higher frex (i.e., nested amplitude frex) 
Phase_Modulating_Freq = 5;
Amp_Modulated_Freq     = 40;

% creation simulation data based on Tort et al. (2010)
lfp=(0.2*(sin(2*pi*t*Phase_Modulating_Freq/srate)+1)+nonmodulatedamplitude*0.1).*sin(2*pi*t*Amp_Modulated_Freq/srate)+sin(2*pi*t*Phase_Modulating_Freq/srate);

% add random noise
lfp=lfp+0.3*randn(1,length(lfp));


% Plot for check
figure(1), clf
subplot(2, 2, 1)
plot(lfp, 'k', 'linew', 1.5), xlim([0 1000])
set(gca,'fontsize',12), xlabel('Time (ms)'), ylabel('Amplitude (a.u.)')
title('Simulation data')


%% Define the range of amplitude and phase frequencies that you want to analyze 
PhaseFreqVector =   2 :1:15;
AmpFreqVector    = 10 :5:100;

PhaseFreq_BandWidth = 2;
AmpFreq_BandWidth    = 5;


%% Create phase bins for storing amplitude values in each phase position
nbin = 10; % we are binning 0-360 degrees based on he number of bin, (e.g., each bin has 20 degrees if we set 18)
position=zeros(1, nbin); % Note: this variable will get the beginning (not the center) of each phase bin (in rads)
winsize = 2*pi/nbin;
for j=1:nbin 
     position(j) = -pi+(j-1)*winsize; 
end

%% Filtering and Hilbert transform
% Data box
Comodulogram  = zeros(length(PhaseFreqVector), length(AmpFreqVector));
[AmpFreqTransformed, PhaseFreqTransformed]   =  deal(zeros(length(AmpFreqVector), data_length), zeros(length(PhaseFreqVector), data_length));

% Loop for extracting phase info
for amp_idx = 1:length(AmpFreqVector)
    
     Af1 = AmpFreqVector(amp_idx);  Af2 = Af1+AmpFreq_BandWidth;
     AmpFreq = eegfilt(lfp, srate, Af1, Af2); % filtering via eegfilt.m
     AmpFreqTransformed(amp_idx, :) = abs(hilbert(AmpFreq)); % getting the amplitude envelope

end

% Loop for extracting amplitude info
for phase_idx =1:length(PhaseFreqVector)
    
    Pf1 = PhaseFreqVector(phase_idx);  Pf2 = Pf1 + PhaseFreq_BandWidth;
    PhaseFreq=eegfilt(lfp, srate, Pf1, Pf2); % filtering via eegfilt.m
    PhaseFreqTransformed(phase_idx , :) = angle(hilbert(PhaseFreq)); % getting the phase time series

end



%% Calculate the wMI and produce comodulogram
tic
counter1=0;
for phase_idx = 1 : length(PhaseFreqVector)

     counter1 = counter1+1;

     % define the target phase frequency bands
     Pf1 = PhaseFreqVector(phase_idx);  Pf2 = Pf1+ PhaseFreq_BandWidth;
    
     counter2=0;
   
    for amp_idx = 1 : length(AmpFreqVector)
         
        counter2=counter2+1;
    
        % define the target amplitude frequency bands 
        Af1 = AmpFreqVector(amp_idx); Af2 = Af1+AmpFreq_BandWidth;
        
        %Calucaltion the wMI via L1wasserstein_MI_test 
        [wMI, optpath, MeanAmp] = Wasserstein_MI(PhaseFreqTransformed(phase_idx, :), AmpFreqTransformed(amp_idx, :), position);
         Comodulogram(counter1, counter2) = wMI;
         
         % just keepling one MeanAmp and coupling phase matrix (i.e., optimal transport route) showing the maximum coupling strength
         if phase_idx == 4 && amp_idx == 7 % 5 Hz for phase and 40 Hz for amplitude
            MeanAmp4plot = MeanAmp;
            optpath4plot = optpath;
            
         end
    
    end
end
toc


%% Plot
figure(1)
subplot(222)
bar(MeanAmp4plot, 'k') %  Phase: 4 Hz  Amplitude: 40 Hz
set(gca,'fontsize',14), ylabel('Amplitude'), xlabel('Phase (binned)')
xticklabels({[]})
title('Amplitude distirbution for phase')

subplot(223)
contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,Comodulogram', 30,'lines','none')
set(gca,'fontsize',14), ylabel('Amplitude (Hz)'), xlabel('Phase(Hz)'), 
colormap turbo,  caxis([.2 .7]), title('wMI'), colorbar

subplot(224)
imagesc(diag(diag(optpath4plot)))
set(gca,'fontsize',14)
title('Coupling phase matrix'), colorbar


