%%% 1 - Prepocessing 
% script pour preparer les données aux analyses
% ce script reçoit des donnees cnt, déjà nettoyé (blink, corrections)
% des variables clés se trouvent au debut du script de façon a ce que les modifications soient simplifiées
% il est important de garder le memes noms de variables
% les donnees sont traités en boucles, simplifiant la tâches
% les donnees sont garder dans des structures
% Les analyses utilisent fieldtrip toolbox

sujets = {'s02', 's03', 's04', 's05', 's06', 's07', 's08', 's09', 's11','s12'};%,'s02','s03','s04','s06','s07','s08','s09', 's11', 's12'}; %inserer le noms des fichier des sujets
%sujets = {'s11', 's12'};
%input_path = ('/Users/sandrajacinto/Documents/MATLAB/ex/Projet01/cnt/08_mai'); %le path jusqu'aux fichiers de donnees raw
%output_path = strcat(input_path,'/fieldtrip_output'); %le path pour les resultats
%eventsOI = [410,411,510,511,630,631,420,421,520,521]; % numero de trigger a analyser 410 - Sc Neutre Predictive Bonne Réponse 630 Sc Scrambled Bonne Réponse
eventsOI = [200,201,210,211,230,231];%,420,421,520,521];
prestim_Time = 0.5; % temps avant stim en sec (baseline)
poststim_Time = 1.3; %temps après stim
load elec

for ss=1:length(sujets);                             
    %if ~exist([output_path '\' sujets{ss}], 'dir'), mkdir([output_path '\' sujets{ss}]); end;  %% créer le dossier 'non_du_sujet' dans output s'il n'existe pas
    for tt= 1:length (eventsOI) % boucle pour chaque trigger
       trig=eventsOI(tt);    %% trig = 

        %% extraire\downsampler les données
        % definir les trials (ft_definetrial) 
        cfg = [];
        cfg.dataset             = (strcat('sandra_',sujets{ss},'.cnt'))   % nom complet du fichier . ds ex ('users\Backup\work/.../raw\s1\s1.cnt')
        cfg.trialdef.eventtype  = 'trigger';       % nom de la voie de trigger
        cfg.trialdef.eventvalue = trig;            % trigger a extraire
        cfg.trialdef.prestim	= prestim_Time;    % information definie plus haut
        cfg.trialdef.poststim	= poststim_Time;   % information definie plus haut
        cfg.channel	   	        = {'EEG' 'EOG'};   	% canaux a conserver 
        cfg = ft_definetrial(cfg);                 

   		% decouper le fichier .cnt, importer dans matlab, a partir des elements definits par ft_definetrial       
        cfg.demean     = 'yes';                    % ligne de base a 0
        cfg.hpfilter      = 'yes' ;
        cfg.hpfreq        = 1;
%        cfg.dftfilter  = 'yes';                    % supprime le bruit a 50 100 150hz
        % Baseline-correction options
        cfg.demean          = 'yes';
        cfg.baselinewindow  = [-0.5 0];
        data       = ft_preprocessing(cfg);        % data = structure avec essais découpés
        
        %resample in a way that all data is 256
        cfg.resamplefs = 1000;
        cfg.detrend = 'no';
        data = ft_resampledata(cfg, data);
        
        switch ss
            case {2}
        
            elec = ft_read_sens('standard_1020.elc');
            cfg = [];
            cfg.method        = 'triangulation';
            cfg.layout        = 'biosemi64.lay';
            cfg.channel       = 'eeg';   %channels for which neighbours should be found
            cfg.feedback      = 'no';
            neighbours = ft_prepare_neighbours(cfg, data); 

            cfg = [];
            cfg.method         = 'average';%, 'average', 'spline' or 'slap' (default = 'weighted')
            cfg.badchannel     = {'F07', 'F08'}; %cell-array, see FT_CHANNELSELECTION for details
           % cfg.missingchannel = cell-array, see FT_CHANNELSELECTION for details
            cfg.neighbours     = neighbours
            cfg.trials         = 'all';% or a selection given as a 1xN vector (default = 'all')
            cfg.lambda         =  1e-5; %default 1e-5 not for method 'distance')
            cfg.order          = 4; %order of the polynomial interpolation (default = 4, not for m
            cfg.elec      = elec;
            ss
            data = ft_channelrepair(cfg, data)
            
            case{3}
                
            cfg = [];
            cfg.method         = 'average';%, 'average', 'spline' or 'slap' (default = 'weighted')
            cfg.badchannel     = {'T7','T8'}; %cell-array, see FT_CHANNELSELECTION for details
           % cfg.missingchannel = cell-array, see FT_CHANNELSELECTION for details
            cfg.neighbours     = neighbours
            cfg.trials         = 'all';% or a selection given as a 1xN vector (default = 'all')
            cfg.lambda         =  1e-5; %default 1e-5 not for method 'distance')
            cfg.order          = 4; %order of the polynomial interpolation (default = 4, not for m
            cfg.elec      = elec;
            ss
            data = ft_channelrepair(cfg, data)
            
            case{4}
            cfg = [];
            cfg.method         = 'average';%, 'average', 'spline' or 'slap' (default = 'weighted')
            cfg.badchannel     = {'M1','T7'}; %cell-array, see FT_CHANNELSELECTION for details
           % cfg.missingchannel = cell-array, see FT_CHANNELSELECTION for details
            cfg.neighbours     = neighbours
            cfg.trials         = 'all';% or a selection given as a 1xN vector (default = 'all')
            cfg.lambda         =  1e-5; %default 1e-5 not for method 'distance')
            cfg.order          = 4; %order of the polynomial interpolation (default = 4, not for m
            cfg.elec      = elec;
            ss
            data = ft_channelrepair(cfg, data)
            
            case{6}
            cfg = [];
            cfg.method         = 'average';%, 'average', 'spline' or 'slap' (default = 'weighted')
            cfg.badchannel     = {'FT7','T7', 'T8'}; %cell-array, see FT_CHANNELSELECTION for details
           % cfg.missingchannel = cell-array, see FT_CHANNELSELECTION for details
            cfg.neighbours     = neighbours
            cfg.trials         = 'all';% or a selection given as a 1xN vector (default = 'all')
            cfg.lambda         =  1e-5; %default 1e-5 not for method 'distance')
            cfg.order          = 4; %order of the polynomial interpolation (default = 4, not for m
            cfg.elec      = elec;
            ss
            data = ft_channelrepair(cfg, data)
            
             case{9}
            cfg = [];
            cfg.method         = 'average';%, 'average', 'spline' or 'slap' (default = 'weighted')
            cfg.badchannel     = {'T7'}; %cell-array, see FT_CHANNELSELECTION for details
           % cfg.missingchannel = cell-array, see FT_CHANNELSELECTION for details
            cfg.neighbours     = neighbours
            cfg.trials         = 'all';% or a selection given as a 1xN vector (default = 'all')
            cfg.lambda         =  1e-5; %default 1e-5 not for method 'distance')
            cfg.order          = 4; %order of the polynomial interpolation (default = 4, not for m
            cfg.elec      = elec;
            ss
            data = ft_channelrepair(cfg, data)
         end
        
        %cd (output_path)
        save(strcat(sujets{ss},'_',num2str(trig)), 'data')    % sauvegarde data avec le nom du trigger ex: save('G:\MEG\dimitri\analysis\benba\data_clean_trig2', 'data')
        %cd ..
        clear data
    end
end
clear

%this script appends triggs so we have pred unpred and scramb conditions

%sujets = {'s01','s02', 's03', 's04', 's05', 's06', 's07', 's08', 's09'};%,'s02','s03','s04','s06','s07','s08','s09'}; %inserer le noms des fichier des sujets
sujets = {'s11', 's12'};
%input_path = ('/Users/sandrajacinto/Documents/MATLAB/ex/Projet01/cnt/08_mai'); %le path jusqu'aux fichiers de donnees raw
%output_path = strcat(input_path,'/fieldtrip_output'); %le path pour les resultats
eventsOI = [630, 631];%[420,421,520,521] %[410,411,510,511]; % [630, 631]; numero de trigger a analyser 410 - Sc Neutre Predictive Bonne Réponse 630 Sc Scrambled Bonne Réponse

for ss=1:length(sujets);                             
    %if ~exist([output_path '\' sujets{ss}], 'dir'), mkdir([output_path '\' sujets{ss}]); end;  %% créer le dossier 'non_du_sujet' dans output s'il n'existe pas
    for tt= 1:length (eventsOI) % boucle pour chaque trigger
        name = (strcat(sujets{ss},'_', int2str(eventsOI(tt)),'.mat'));
        load (name);
        data_s{tt} = data;
    end
    cfg= [];
    data_Scramb =  ft_appenddata(cfg, data_s{1}, data_s{2})%, data_s{3}, data_s{4});
    save ((strcat(sujets{ss},'_Scramb.mat')), 'data_Scramb');
    clear data_s data_Scramb data
end

eventsOI = [420,421,520,521]; %[410,411,510,511]; % [630, 631]; numero de trigger a analyser 410 - Sc Neutre Predictive Bonne Réponse 630 Sc Scrambled Bonne Réponse

for ss=1:length(sujets);                             
    %if ~exist([output_path '\' sujets{ss}], 'dir'), mkdir([output_path '\' sujets{ss}]); end;  %% créer le dossier 'non_du_sujet' dans output s'il n'existe pas
    for tt= 1:length (eventsOI) % boucle pour chaque trigger
        name = (strcat(sujets{ss},'_', int2str(eventsOI(tt)),'.mat'));
        load (name);
        data_s{tt} = data;
    end
    cfg= [];
    data_UnPred =  ft_appenddata(cfg, data_s{1}, data_s{2}, data_s{3},data_s{4})%, data_s{3}, data_s{4});
    save ((strcat(sujets{ss},'_UnPred.mat')), 'data_UnPred');
    clear data_s data_UnPred data
end

eventsOI = [410,411,510,511]; % [630, 631]; numero de trigger a analyser 410 - Sc Neutre Predictive Bonne Réponse 630 Sc Scrambled Bonne Réponse

for ss=1:length(sujets);                             
    %if ~exist([output_path '\' sujets{ss}], 'dir'), mkdir([output_path '\' sujets{ss}]); end;  %% créer le dossier 'non_du_sujet' dans output s'il n'existe pas
    for tt= 1:length (eventsOI) % boucle pour chaque trigger
        name = (strcat(sujets{ss},'_', int2str(eventsOI(tt)),'.mat'));
        load (name);
        data_s{tt} = data;
    end
    cfg= [];
    data_Pred =  ft_appenddata(cfg, data_s{1}, data_s{2}, data_s{3},data_s{4})%, data_s{3}, data_s{4});
    save ((strcat(sujets{ss},'_Pred.mat')), 'data_Pred');
    clear data_s data_Pred data
end
%% Timelockanalysis

sujets = {'s01', 's02', 's03','s04', 's06', 's07', 's08', 's09', 's11', 's12'};%,'s02','s03','s04','s06','s07','s08','s09'}; %inserer le noms des fichier des sujets
%input_path = ('/Users/sandrajacinto/Documents/MATLAB/ex/Projet01/cnt/08_mai'); %le path jusqu'aux fichiers de donnees raw
%output_path = strcat(input_path,'/fieldtrip_output'); %le path pour les resultats
eventsOI = {'_UnPred','_Pred','_Scramb'};%,510,511]; % numero de trigger a analyser 410 - Sc Neutre Predictive Bonne Réponse 630 Sc Scrambled Bonne Réponse

for ss=1:length(sujets);                             
    %if ~exist([output_path '\' sujets{ss}], 'dir'), mkdir([output_path '\' sujets{ss}]); end;  %% créer le dossier 'non_du_sujet' dans output s'il n'existe pas
    for tt= 1:length (eventsOI) % boucle pour chaque trigger
        data = load (strcat(sujets{ss},(eventsOI{tt}),'.mat')); data = data.(strcat('data',eventsOI{tt}));
        
    cfg = [];
    cfg.keeptrials = 'no';
    dataAv = ft_timelockanalysis(cfg, data)
    
    save ((strcat(sujets{ss},(eventsOI{tt}),'_Av')), 'dataAv')
    cfg=[];
    cfg.layout='biosemi64.lay';
    figure; ft_multiplotER (cfg, dataAv)
    end
end



%%
load('s01_Pred_Av.mat')
data1 = dataAv;
load('s02_Pred_Av.mat')
data2 = dataAv;
load('s03_Pred_Av.mat')
data3 = dataAv;
load('s04_Pred_Av.mat')
data4 = dataAv;
load('s06_Pred_Av.mat')
data6 = dataAv;
load('s07_Pred_Av.mat')
data7 = dataAv;
load('s08_Pred_Av.mat')
data8 = dataAv;
load('s09_Pred_Av.mat')
data9 = dataAv;

cfg = [];
datAv = ft_timelockgrandaverage (cfg, data1, data2, data3,data4,data6,data7,data8,data9)

cfg = [];
cfg.layout='biosemi64.lay';
figure; ft_multiplotER (cfg, data)



