clear;clc;  

Tau=1000:1000:20000;
Tau1=1000:1000:10000;
LTP=3;
Pro=0;

 directory=[cd,'\image\'];
 EEG=cell(length(LTP),length(Pro));
%  addpath([pwd, '\peaks-detection-master\']);
  for k=1:length(Tau)
     for l=1:length(Tau1)
         parset_tau_rec_e=Tau(k);
         parset_tau_rec_i=Tau1(l);
         parset_ltp_factor=LTP;
         parset_propofol=Pro;
         % parset_titletext='Example 5: irregular discharges';
         parset_ylimits=[-60 40];
         EEG{k,l}= mmc6(parset_tau_rec_e,parset_tau_rec_i,parset_ltp_factor,parset_propofol);
         EEG_new = downsample(EEG{k,l},70);
         fs=8960/70;						% sample frequency after downsampling
         EEG_new_1 = EEG_new((20*fs)+1:end); 	% remove transient effects
         filtOrder=3;
         filtFreq= [0.5 30];
         [b,a]= butter(filtOrder,filtFreq/(0.5*fs),'bandpass');
         EEG_filtered = filtfilt(b,a,EEG_new_1); % apply bandpass filter
    
%         
    h = figure;
    ty = 1/fs:1/fs:length(EEG_filtered)/fs;
    plot(ty, EEG_filtered,'color','k');
    set(gca, 'FontSize',14)

     saveas(h,fullfile(directory,['Im_', num2str(k),'_', num2str(l), '.fig']));


     end

  end
  
