%sujets = { 's01','s02','s03','s04','s05', 's06', 's07', 's08', 's09'};%,'s02','s03','s04','s06','s07','s08','s09'}; %inserer le noms des fichier des sujets
sujets = {'s11', 's12'};

eventsOI = {'_UnPred','_Pred','_Scramb'};



min_freq =  2;
max_freq = 80;
num_frex = 40;
freq = linspace(min_freq,max_freq,num_frex);
load label
load dimord

for ss = 1:length(sujets)
    for ee = 1:length (eventsOI)
        load(strcat ('EEG_',sujets{ss},eventsOI{ee}));
        name = strcat ('EEG_TF_',sujets{ss},eventsOI{ee});
        data = load(name); data=data.(name);
        TF = [];
        TF.label = label;
        TF.dimord = dimord;
        TF.freq = freq;
        TF.time = EEG.times;
        TF.powspctrm = data;
        TF.cfg = 'imported from cohen TF';
        namesave = strcat ('TFp_',sujets{ss},eventsOI{ee});
        save ((strcat(namesave,'.mat')), 'TF')
    end
    
end
