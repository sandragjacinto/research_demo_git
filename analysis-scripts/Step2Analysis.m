

%% loop all participants

subjects = {'EEG_01'}; %Define Subjects

cond ={'_Pred_Ac','_NPred_Ac','_Pred_In','_NPred_In'};
ext={'.mat'}; % extension

% variable number of wavelet cycles & wavelet parameters
%% 
num_frex = 150;
min_freq = 15;
max_freq = 90;
% set a few different wavelet widths (number of wavelet cycles)
range_cycles = [ 3 12 ];

% other wavelet parameters
nCycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex);

%baseline window
baseline_window = [ -500 -100 ];

[~,baseidx(1)] = min(abs(EEG.times-baseline_window(1)));
[~,baseidx(2)] = min(abs(EEG.times-baseline_window(2)));

%Folder='/run/media/bouth/MACBOU_D/data/Noise'; %get path to the main directory


% initialize output time-frequency data

% FFT of data (doesn't change on frequency iteration), here we concatenate
% trials
% data=reshape(EEG.data(strcmpi(channel2use,{EEG.chanlocs.labels}),:,:),1,[]);
% dataX=fft(data,nConv);


%% Scipt (1) all trials

for jj = 1 :1% numel(subjects)  %
   % subjectPath= eval(['''' Folder '/' subjects{jj} '''']);
    %cd(subjectPath); % go to subject directory
    for ii = 1 :1% numel(cond)
      %  name= strcat(char(subjects{jj}),(cond{ii}), ext);
     %   eval([ 'load ' subjects{jj} cond{ii} '.mat']); %Victor %Define name of current mat file
        % variable number of wavelet cycles & wavelet parameters
        frex = linspace(min_freq,max_freq,num_frex);
        time = -2:1/EEG.srate:2;
        half_wave = (length(time)-1)/2;

        % FFT parameters
        nKern = length(time);
        nData = EEG.pnts*EEG.trials;
        nConv = nKern+nData-1;
        tf_chan = zeros(length(EEG.chanlocs),length(frex),EEG.pnts);

            for chan=1:64
                dataX = fft(reshape(EEG.data(chan,:,:),1,[]),nConv);
                tf = zeros(length(frex),EEG.pnts,EEG.trials);
                    for fi=1:length(frex)
                         s = nCycles(fi)/(2*pi*frex(fi));
                         cmw = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));    
                         cmwX = fft(cmw,nConv);   
                      % run convolution
                         as = ifft(cmwX.*dataX,nConv);
                         as = as(half_wave+1:end-half_wave);
                         as = reshape(as,EEG.pnts,EEG.trials);
                         tf(fi,:,:)=as;
                    end
        % normalize and baseline correction
                 tmpbase=1:size(tf,2);
                 mstd = std(tf(:,tmpbase,:),[],2);
                 mbase = mean(tf(:,tmpbase,:),2);
                 P = bsxfun(@rdivide, bsxfun(@minus, tf, mbase), mstd);
                 tfCohen = mean(abs(P).^2,3);
                 tfCohen = 10*log10( bsxfun(@rdivide, tfCohen, mean(tfCohen(:,baseidx(1):baseidx(2)),2)) );
                 tf_chan(chan,:,:)=tfCohen; 
            end                   
          % FILENAME= strcat(char(subjects{jj}), (cond{ii}), '_TFRpow'); % here we try to create a *.mat file name holding the subject's number, the condition and type of analysis(TF here)
%            eval([FILENAME '= tf_chan;']);  % it is like renaming the template TF by using the new *.mat filename (FILENAME1)
%            eval(['save ' FILENAME ' ' FILENAME]); % here we save the *.mat filename
%            clear FILENAME a0* a1* c0* c1* as0* as*   
    end
end
%% Plot

figure;

%% Scipt (1)  non-phase-locked activity

subjects = {'a01','a02','a03','a04','a05','a06','a07',...
    'a08','a09','a10','a11','a12','a13','a14','c02','c03','c04',...
    'c05','c06','c07','c08','c09','c10','c11','c14','as01',...
    'as02','as03','as04','as05','as06','as07','as08','as09',...
    'as10','as11','as12','as13','as14'}; %Define Subjects
subjects = {'as07','as08','as09',...
    'as10','as11','as12','as13','as14'}; %Define Subjects

cond ={'_SignalCued','_NoSignalCued','_SignalNoCue','_NoSignalNoCue'};
ext={'.mat'}; % extension

% variable number of wavelet cycles & wavelet parameters
% 
num_frex = 70;
min_freq =  2;
max_freq = 80;
% set a few different wavelet widths (number of wavelet cycles)
range_cycles = [ 3 12 ];

% other wavelet parameters
nCycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex);

%baseline window
baseline_window = [ -300 -100 ];
[~,baseidx(1)] = min(abs(EEG.times-baseline_window(1)));
[~,baseidx(2)] = min(abs(EEG.times-baseline_window(2)));


Folder='/run/media/bouth/MACBOU_D/data/Noise'; %get path to the main directory

%% loop Scipt (1)  non-phase-locked activity

for jj = 1 : numel(subjects)  %
    subjectPath= eval(['''' Folder '/' subjects{jj} '''']);
    cd(subjectPath); % go to subject directory
    for ii = 1 : numel(cond)
        name= strcat(char(subjects{jj}),(cond{ii}), ext);
        eval([ 'load ' subjects{jj} cond{ii} '.mat']); %Victor %Define name of current mat file
        % variable number of wavelet cycles & wavelet parameters
        frex = linspace(min_freq,max_freq,num_frex);
        time = -2:1/EEG.srate:2;
        half_wave = (length(time)-1)/2;

        % FFT parameters
        nKern = length(time);
        nData = EEG.pnts*EEG.trials;
        nConv = nKern+nData-1;
        % compute ERP
        erp = squeeze(mean(EEG.data(:,:,:),3));

        % compute non-phase-locked power by subtracting ERP from each trial
        nonphase_EEG = squeeze( bsxfun(@minus,EEG.data(:,:,:),erp) );
       
        tf_chan = zeros(length(EEG.chanlocs),length(frex),EEG.pnts);

            for chan=1:58
                dataX = fft(reshape(nonphase_EEG(chan,:,:),1,[]),nConv);
                tf = zeros(length(frex),EEG.pnts,EEG.trials);
                    for fi=1:length(frex)
                         s = nCycles(fi)/(2*pi*frex(fi));
                         cmw = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));    
                         cmwX = fft(cmw,nConv);   
                      % run convolution
                         as = ifft(cmwX.*dataX,nConv);
                         as = as(half_wave+1:end-half_wave);
                         as = reshape(as,EEG.pnts,EEG.trials);
                         tf(fi,:,:)=as;
                    end
        % normalize and baseline correction
                 tmpbase=1:size(tf,2);
                 mstd = std(tf(:,tmpbase,:),[],2);
                 mbase = mean(tf(:,tmpbase,:),2);
                 P = bsxfun(@rdivide, bsxfun(@minus, tf, mbase), mstd);
                 tfCohen = mean(abs(P).^2,3);
                 tfCohen = 10*log10( bsxfun(@rdivide, tfCohen, mean(tfCohen(:,baseidx(1):baseidx(2)),2)) );
                 tf_chan(chan,:,:)=tfCohen; 
            end                   
           FILENAME= strcat(char(subjects{jj}), (cond{ii}), '_TFRpowInd'); % here we try to create a *.mat file name holding the subject's number, the condition and type of analysis(TF here)
           eval([FILENAME '= tf_chan;']);  % it is like renaming the template TF by using the new *.mat filename (FILENAME1)
           eval(['save ' FILENAME ' ' FILENAME]); % here we save the *.mat filename
           clear FILENAME a0* a1* c0* c1* as0* as*   
    end
end



%% Scipt (2)  non-phase-locked activity

for jj = 1 : numel(subjects)  % 
    for ii = 1 : numel(cond)
        name= strcat(char(subjects{jj}),(cond{ii}), ext);
        eval([ 'load ' subjects{jj} cond{ii} '.mat']); %Victor %Define name of current mat file
        % variable number of wavelet cycles & wavelet parameters
        frex = linspace(min_freq,max_freq,num_frex);
        time = -2:1/EEG.srate:2;
        half_wave = (length(time)-1)/2;

        % FFT parameters
        nKern = length(time);
        nData = EEG.pnts*EEG.trials;
        nConv = nKern+nData-1;

        % convert baseline time into indices
        [~,baseidx(1)] = min(abs(EEG.times-baseline_window(1)));
        [~,baseidx(2)] = min(abs(EEG.times-baseline_window(2)));
        % compute ERP
        erp = squeeze(mean(EEG.data(:,:,:),3));

        % compute non-phase-locked power by subtracting ERP from each trial
        nonphase_EEG = squeeze( bsxfun(@minus,EEG.data(:,:,:),erp) );


        % FFT of data (doesn't change on frequency iteration), here we concatenate
        % trials
        % data=reshape(EEG.data(strcmpi(channel2use,{EEG.chanlocs.labels}),:,:),1,[]);
        % dataX=fft(data,nConv);
        
            for chan=1:58
                dataX = fft(reshape(nonphase_EEG(chan,:,:),1,[]),nConv);
                tf = zeros(length(frex),EEG.pnts);
                    for fi=1:length(frex)
                         s = nCycles(fi)/(2*pi*frex(fi));
                         cmw = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));    
                         cmwX = fft(cmw,nConv);   
                      % run convolution
                         as = ifft(cmwX.*dataX,nConv);
                         as = as(half_wave+1:end-half_wave);
                         as = reshape(as,EEG.pnts,EEG.trials);
                       % put power data into big matrix
                         tf(fi,:) = mean(abs(as).^2,2);
                    end
                 tf = 10*log10( bsxfun(@rdivide, tf, mean(tf(:,baseidx(1):baseidx(2)),2)) ); %db conversion
                 tf_chan(chan,:,:)= tf;
            end   
          FILENAME= strcat(char(subjects{jj}), (cond{ii}), '_TFRComplexInd'); % here we try to create a *.mat file name holding the subject's number, the condition and type of analysis(TF here)
          eval([FILENAME '= tf_chan;']);  % it is like renaming the template TF by using the new *.mat filename (FILENAME1)
          eval(['save ' FILENAME ' ' FILENAME]); % here we save the *.mat filename
    end
end

%% 




%% Scipt (1)

for jj = 1 : numel(subjects)  % 
    for ii = 1 : numel(cond)
        name= strcat(char(subjects{jj}),(cond{ii}), ext);
        eval([ 'load ' subjects{jj} cond{ii} '.mat']); %Victor %Define name of current mat file
        % variable number of wavelet cycles & wavelet parameters
        frex = linspace(min_freq,max_freq,num_frex);
        time = -2:1/EEG.srate:2;
        half_wave = (length(time)-1)/2;

        % FFT parameters
        nKern = length(time);
        nData = EEG.pnts*EEG.trials;
        nConv = nKern+nData-1;

        % convert baseline time into indices
        [~,baseidx(1)] = min(abs(EEG.times-baseline_window(1)));
        [~,baseidx(2)] = min(abs(EEG.times-baseline_window(2)));
        tf_chan = zeros(58,length(frex),EEG.pnts);
        
            for chan=1:58
                dataX = fft(reshape(EEG.data(chan,:,:),1,[]),nConv);
                tf = zeros(length(frex),EEG.pnts);
                    for fi=1:length(frex)
                         s = nCycles(fi)/(2*pi*frex(fi));
                         cmw = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));    
                         cmwX = fft(cmw,nConv);   
                      % run convolution
                         as = ifft(cmwX.*dataX,nConv);
                         as = as(half_wave+1:end-half_wave);
                         as = reshape(as,EEG.pnts,EEG.trials);
                       % put power data into big matrix
                         tf(fi,:) = mean(abs(as).^2,2);
                    end
                 tf = 10*log10( bsxfun(@rdivide, tf, mean(tf(:,baseidx(1):baseidx(2)),2)) ); %db conversion
                 tf_chan(chan,:,:)= tf;
            end   
          FILENAME= strcat(char(subjects{jj}), (cond{ii}), '_TFRpow'); % here we try to create a *.mat file name holding the subject's number, the condition and type of analysis(TF here)
          eval([FILENAME '= tf_chan;']);  % it is like renaming the template TF by using the new *.mat filename (FILENAME1)
          eval(['save ' FILENAME ' ' FILENAME]); % here we save the *.mat filename
    end
end

%% Scipt (2)  non-phase-locked activity

for jj = 1 : numel(subjects)  % 
    for ii = 1 : numel(cond)
        name= strcat(char(subjects{jj}),(cond{ii}), ext);
        eval([ 'load ' subjects{jj} cond{ii} '.mat']); %Victor %Define name of current mat file
        % variable number of wavelet cycles & wavelet parameters
        frex = linspace(min_freq,max_freq,num_frex);
        time = -2:1/EEG.srate:2;
        half_wave = (length(time)-1)/2;

        % FFT parameters
        nKern = length(time);
        nData = EEG.pnts*EEG.trials;
        nConv = nKern+nData-1;

        % convert baseline time into indices
        [~,baseidx(1)] = min(abs(EEG.times-baseline_window(1)));
        [~,baseidx(2)] = min(abs(EEG.times-baseline_window(2)));
        % compute ERP
        erp = squeeze(mean(EEG.data(:,:,:),3));

        % compute non-phase-locked power by subtracting ERP from each trial
        nonphase_EEG = squeeze( bsxfun(@minus,EEG.data(:,:,:),erp) );


        % FFT of data (doesn't change on frequency iteration), here we concatenate
        % trials
        % data=reshape(EEG.data(strcmpi(channel2use,{EEG.chanlocs.labels}),:,:),1,[]);
        % dataX=fft(data,nConv);
        
            for chan=1:58
                dataX = fft(reshape(nonphase_EEG(chan,:,:),1,[]),nConv);
                tf = zeros(length(frex),EEG.pnts);
                    for fi=1:length(frex)
                         s = nCycles(fi)/(2*pi*frex(fi));
                         cmw = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));    
                         cmwX = fft(cmw,nConv);   
                      % run convolution
                         as = ifft(cmwX.*dataX,nConv);
                         as = as(half_wave+1:end-half_wave);
                         as = reshape(as,EEG.pnts,EEG.trials);
                       % put power data into big matrix
                         tf(fi,:) = mean(abs(as).^2,2);
                    end
                 tf = 10*log10( bsxfun(@rdivide, tf, mean(tf(:,baseidx(1):baseidx(2)),2)) ); %db conversion
                 tf_chan(chan,:,:)= tf;
            end   
          FILENAME= strcat(char(subjects{jj}), (cond{ii}), '_TFRpowind'); % here we try to create a *.mat file name holding the subject's number, the condition and type of analysis(TF here)
          eval([FILENAME '= tf_chan;']);  % it is like renaming the template TF by using the new *.mat filename (FILENAME1)
          eval(['save ' FILENAME ' ' FILENAME]); % here we save the *.mat filename
    end
end

%% 