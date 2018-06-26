%Grand Average
sujets = { 's01','s02','s03','s04', 's05','s06', 's07', 's08', 's09'};%,'s02','s03','s04','s06','s07','s08','s09'}; %inserer le noms des fichier des sujets
eventsOI = {'_UnPred','_Pred','_Scramb'};
prename = 'TFP_';

for ee = 1:length(eventsOI)
for ss = 2:length(sujets)
    name = strcat (prename, sujets{ss}, eventsOI{ee});
    temp = load(name); temp = temp.TF;
    tempname = strcat('freq', num2str(ss));
    eval ([tempname '=temp']);
end
cfg = [];
cfg.keepindividual = 'no';
cfg.foilim         = 'all';
cfg.toilim         = 'all';
cfg.channel        = 'all';

GA = ft_freqgrandaverage(cfg, freq2, freq3, freq4, freq5, freq6, freq7, freq8, freq9);
GAname = strcat ('GA', eventsOI{ee});

eval ([GAname '=GA']);
end

%%
cfg = [];
cfg.layout='biosemi64.lay';
figure; ft_multiplotTFR (cfg, GA_Pred)
