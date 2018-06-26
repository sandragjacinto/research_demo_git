%% TF using trial concatenation
%sujets = {'s01','s02','s03','s04','s05','s06','s07','s08','s09'}; %inserer le noms des fichier des sujets
sujets = {'s11', 's12'};
eventsOI = {'_Scramb', '_Pred', '_UnPred'};


baseline_windows = [ -0.500 0];


min_freq =  2;
max_freq = 80;
num_frex = 40;
frex = linspace(min_freq,max_freq,num_frex);

% other wavelet parameters
range_cycles = [ 10 16 ];

s = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex) ./ (2*pi*frex);


tic; % restart matlab timer
for ss = 1:length(sujets)
    for tt = 1:length(eventsOI)
        load (strcat('EEG_',sujets{ss},(eventsOI{tt}),'.mat'));
        avance = strcat('EEG_',sujets{ss},(eventsOI{tt}),'.mat')
        % convert baseline time into indices
        baseidx = reshape( dsearchn(EEG.times',baseline_windows(:)), [],2);
        wavtime = -2:1/EEG.srate:2;
        half_wave = (length(wavtime)-1)/2;
        
        
        for chanInd = 1:64
            % which channel to plot

            % FFT parameters
            nWave = length(wavtime);
            nData = EEG.pnts * EEG.trials; % This line is different from above!!
            nConv = nWave + nData - 1;
            
            % now compute the FFT of all trials concatenated
            alldata = reshape( EEG.data(chanInd,:,:) ,1,[]);
            dataX   = fft( alldata ,nConv );
            %tf = zeros(length(frex),EEG.pnts);
            tf= zeros(length(frex),EEG.pnts,EEG.trials);
            
            % loop over frequencies
            for fi=1:length(frex)
                
                % create wavelet and get its FFT
                % the wavelet doesn't change on each trial...
                wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
                waveletX = fft(wavelet,nConv);
                waveletX = waveletX ./ max(waveletX);
                
                % now run convolution in one step
                as = ifft(waveletX .* dataX);
                as = as(half_wave+1:end-half_wave);
                
                % and reshape back to time X trials
                as = reshape( as, EEG.pnts, EEG.trials );
                tf(fi,:,:) =as;
                % compute power and average over trials
                % tf(fi,:) = mean( abs(as).^2 ,2);
                
            end
            tmpbase=1:size(tf,2);
            mstd = std(tf(:,tmpbase,:),[],2);
            mbase = mean(tf(:,tmpbase,:),2);
            P = bsxfun(@rdivide, bsxfun(@minus, tf, mbase), mstd);
            tfCohen = mean(abs(P).^2,3);
            tfCohen = 10*log10( bsxfun(@rdivide, tfCohen, mean(tfCohen(:,baseidx(1):baseidx(2)),2)) );
            TF_data(chanInd,:,:)=tfCohen;
            logichan = chanInd
        end
       
        name = strcat('EEG_TF_',sujets{ss},(eventsOI{tt}));
        eval([name '= TF_data;'])
        save ((strcat(name,'.mat')), (name)) 
        clear TF_data tfCohen as
        
    end
    
end
 computationTime(2) = toc
% % plot results
% figure(2), clf
% contourf(EEG.times,frex,tf,40,'linecolor','none')
% set(gca,'ydir','normal','xlim',[-0.300 1.3000])
%%

tmpbase=1:size(tf,2);
mstd = std(tf(:,tmpbase,:),[],2);
mbase = mean(tf(:,tmpbase,:),2);
P = bsxfun(@rdivide, bsxfun(@minus, tf, mbase), mstd);
tfCohen = mean(abs(P).^2,3);
tfCohen = 10*log10( bsxfun(@rdivide, tfCohen, mean(tfCohen(:,baseidx(1):baseidx(2)),2)) );
%tf_chan(:,:)=tfCohen;

figure, clf
contourf(EEG.times,frex,tfCohen,40,'linecolor','none')
%% db conversion and plot results
%%Baseline

% define color limits
climdb  = [-3 3];

% create new matrix for percent change
tfpct = zeros(size(tf));

for basei=1:size(tf,1)
    
    activity = tf(4,:,:);
    baseline = mean( tf(4,:,baseidx(basei,1):baseidx(basei,2)) ,3);
    
    % decibel
    tf(basei,:,:) = 10*log10( bsxfun(@rdivide, activity, baseline) );
    
end


% plot results
figure(2), clf
contourf(EEG.times,frex,tf,40,'linecolor','none')
%set(gca,'clim',[0 5],'ydir','normal','xlim',[-300 1000])

%% show computation time differences

figure(3), clf
bar(computationTime)
set(gca,'xtick',1:2,'xticklabel',{'Each trial';'Concatenated'},'xlim',[0 3])
ylabel('Computation time (s)')

%% show differences

figure(4), clf
plot(EEG.times,tfTrialAve(5,:)), hold on
plot(EEG.times,tf(5,:),'r')
set(gca,'xlim',[-500 1200])

title([ 'Power at ' num2str(round(10*frex(5))/10) ' Hz' ])
xlabel('Time (ms)'), ylabel('Power (\muV^2)')
legend({'each trial';'concatenated'})

%% end
