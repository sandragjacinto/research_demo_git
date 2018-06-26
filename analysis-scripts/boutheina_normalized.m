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
                chan=30;
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
                 tf_chan(chan,:,:) = tfCohen; 
            end                   
           FILENAME= strcat(char(subjects{jj}), (cond{ii}), '_TFRpowInd'); % here we try to create a *.mat file name holding the subject's number, the condition and type of analysis(TF here)
           eval([FILENAME '= tf_chan;']);  % it is like renaming the template TF by using the new *.mat filename (FILENAME1)
           eval(['save ' FILENAME ' ' FILENAME]); % here we save the *.mat filename
           clear FILENAME a0* a1* c0* c1* as0* as*   
    end
end
