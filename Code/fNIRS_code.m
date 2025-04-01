clear
load('rhythm_R.mat') % R: right channels / L: left channels
load('arrhythm_R.mat') % R: right channels / L: left channels

load('distancesRight.mat')
load('distancesLeft.mat')
dist = distR; % distR for right side or distL for left side
nch = 45; % 45 for right or 44 for left

fs = 10.1215;
HbOfinaltot_rhythm(:,:,10) = [];
HbOfinaltot_arrhythm(:,:,10) = [];
for k=1:11  
    for i = 1:nch
        HbOfinaltot_rhythm(:,i,k) = HbOfinaltot_rhythm(:,i,k) ./ dist(1,i);
        HbOfinaltot_arrhythm(:,i,k) = HbOfinaltot_arrhythm(:,i,k) ./ dist(1,i);
    end
end
temp = HbOfinaltot_rhythm;
temp = permute(temp,[2,1,3]);
temp = pop_importdata('data',temp,'srate',fs);
temp = pop_resample(temp,1);
HbOfinaltot_rhythm = permute(temp.data,[2,1,3]);
time = temp.times/1000 - 5;
temp = HbOfinaltot_arrhythm;
temp = permute(temp,[2,1,3]);
temp = pop_importdata('data',temp,'srate',fs);
temp = pop_resample(temp,1);
HbOfinaltot_arrhythm = permute(temp.data,[2,1,3]);
%baseline
nsub = 11;
for k =1:nsub
for i =1:nch

    y1 = HbOfinaltot_rhythm(:,i,k);
    yy = nanmean(y1(1:5));
    x1b(:,i,k) = yy*ones(1,31);

    y2 = HbOfinaltot_arrhythm(:,i,k);
    yy = nanmean(y2(1:5));
    x2b(:,i,k) = yy*ones(1,31);

end
end
x1 = HbOfinaltot_rhythm;
x1 = permute(x1,[2 1 3]);
x1b = permute(x1b,[2 1 3]);
x2 = HbOfinaltot_arrhythm;
x2 = permute(x2,[2 1 3]);
x2b = permute(x2b,[2 1 3]);
fs = 1;
%%
nch = 45;
dist = distR; % distR for right side or distL for left side
x1([9 13 14 15 17 18 19 20 21 22 23 24 25 31 38 39 40 43],:,:) = nan;
x1b([9 13 14 15 17 18 19 20 21 22 23 24 25 31 38 39 40 43],:,:) = nan;
x2([9 13 14 15 17 18 19 20 21 22 23 24 25 31 38 39 40 43],:,:) = nan;
x2b([9 13 14 15 17 18 19 20 21 22 23 24 25 31 38 39 40 43],:,:) = nan;
load('Neighbours_Right.mat');

% nch = 44;
% dist = distL; % distR for right side or distL for left side
% x1([6 9 13 16 17 18 21 23 30 33 40],:,:) = nan;
% x1b([6 9 13 16 17 18 21 23 30 33 40],:,:) = nan;
% x2([6 9 13 16 17 18 21 23 30 33 40],:,:) = nan;
% x2b([6 9 13 16 17 18 21 23 30 33 40],:,:) = nan;
% load('Neighbours_Left.mat');
%% cluster-based statistical learning
sub = 1:11;
load('OldFieldLayout.mat');
load('chan128.mat');

sig{1} = x1; %x1,x2
sig{2} = x2; %x2,x1b,x2b

Rhyt = cell(length(sub),1);
for i=1:length(Rhyt)
    EEG = pop_importdata('data',sig{1}(:,:,i),'srate',fs);
    EEG.chanlocs = chan(1:nch);
    
    Rhyt{i} = eeglab2fieldtrip(EEG,'preprocessing');
    Rhyt{i}.trial = EEG.data;    
    Rhyt{i}.avg   = EEG.data;
    Rhyt{i}.dimord = 'chan_time';
    Rhyt{i}.time = time;  % time or 1  
    tmp = Rhyt{i}.elec.pnt(:,1);
    Rhyt{i}.elec.pnt(:,1) = -Rhyt{i}.elec.pnt(:,2);
    Rhyt{i}.elec.pnt(:,2) = tmp;             
end

cRhyt = cell(length(sub),1);
for i=1:length(cRhyt)
    EEG = pop_importdata('data',sig{2}(:,:,i),'srate',fs);
    EEG.chanlocs = chan(1:nch);
    
    cRhyt{i} = eeglab2fieldtrip(EEG,'preprocessing');
    cRhyt{i}.trial = EEG.data;    
    cRhyt{i}.avg   = EEG.data;
    cRhyt{i}.dimord = 'chan_time';
    cRhyt{i}.time = time; % time or 1       
    tmp = cRhyt{i}.elec.pnt(:,1);
    cRhyt{i}.elec.pnt(:,1) = -cRhyt{i}.elec.pnt(:,2);
    cRhyt{i}.elec.pnt(:,2) = tmp;             
end

% for k=1:nch  % channel by channel : uncomment
    cfg = [];
    cfg.channel     = [{chan(1:nch).labels}]; % chan(k).labels;  % [{chan(1:nch).labels}];
    cfg.neighbours  = neighbours; % defined as above   % channel by channel : comment 
    cfg.latency     = [5 15]; % [5 15] or [1]
    cfg.avgovertime = 'no';
    cfg.parameter   = 'avg';
    cfg.method      = 'montecarlo';
    cfg.statistic   = 'ft_statfun_depsamplesT';
    cfg.alpha       = 0.05;
    cfg.clusteralpha= 0.05;
    cfg.correctm    = 'cluster';
    cfg.correcttail = 'prob';
    cfg.tail = 0;
    cfg.numrandomization = 5000;
    cfg.minnbchan        = 1; % minimal neighbouring channels           % channel by channel : comment 

    Nsub = length(sub);
    cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
    cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
    cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
    cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

    stat = ft_timelockstatistics(cfg,Rhyt{:},cRhyt{:}); % channel by channel : stat{k} = ft_timelockstatistics(cfg,Rhyt{:},cRhyt{:});  % stat = ft_timelockstatistics(cfg,Rhyt{:},cRhyt{:});
% end     % channel by channel : uncomment

% % specify significant channels
% for k = 1:nch   % channel by channel : uncomment
%     cfg = [];
%     cfg.layout = layout;
%     cfg.channel = 'all';
%     cfg.colorbar = 'yes';
%     cfg.renderer = 'painters';
%     cfg.parameter = 'stat';
%     cfg.maskparameter = 'mask';
%     cfg.maskalpha = 1;
%     cfg.maskstyle = 'box';
%     cfg.colorbartext = 'Power (dB)';
%     cfg.showlabels = 'yes';
%     cfg.showoutline = 'yes';
%     cfg.colormap = 'jet';
%     cfg.masknans = 'yes';
%     figure
%     ft_multiplotER(cfg,stat{k}); % stat{k} % stat
% end    % channel by channel : uncomment

%% plot for clusters

% R
temp1 = HbOfinaltot_rhythm;
temp1 = temp1(:,[27 30 41 42 44 45],:); % 27 30 41 42 44 45  % 3 4 16 35
HbOfinaltot_rhythm = nansum(temp1,2) ./ sum(1 ./ dist(1,[27 30 41 42 44 45]));
temp2 = HbOfinaltot_arrhythm;
temp2 = temp2(:,[27 30 41 42 44 45],:); % 27 30 41 42 44 45  % 3 4 16 35
HbOfinaltot_arrhythm = nansum(temp2,2) ./ sum(1 ./ dist(1,[27 30 41 42 44 45]));

% % L
% temp1 = HbOfinaltot_rhythm;
% temp1 = temp1(:,[12 35 38],:); % 12 35 38  % 1 2 14 26 32
% HbOfinaltot_rhythm = nansum(temp1,2) ./ sum(1 ./ dist(1,[12 35 38]));
% 
% temp2 = HbOfinaltot_arrhythm;
% temp2 = temp2(:,[12 35 38],:); % 12 35 38  % 1 2 14 26 32
% HbOfinaltot_arrhythm = nansum(temp2,2) ./ sum(1 ./ dist(1,[12 35 38]));

HbO_rhythm = nanmean(HbOfinaltot_rhythm,3); 
HbO_arrhythm = nanmean(HbOfinaltot_arrhythm,3);

figure
i = 1;
ch = 1;
% for ch = 1:nch
    HbOstd_rhythm(:,ch) = std(squeeze(HbOfinaltot_rhythm(:,ch,:)),[],2,'omitnan')/sqrt(sum(~isnan(HbOfinaltot_rhythm(1,ch,:))));
    numSUB_avg(ch) = sum(~isnan(HbOfinaltot_rhythm(1,ch,:))); % avg over this subjects for channel ch
    HbOstd_arrhythm(:,ch) = std(squeeze(HbOfinaltot_arrhythm(:,ch,:)),[],2,'omitnan')/sqrt(sum(~isnan(HbOfinaltot_arrhythm(1,ch,:))));
    numSUB_avg2(ch) = sum(~isnan(HbOfinaltot_arrhythm(1,ch,:)));
    shadedErrorBar(time',HbO_rhythm(:,ch),HbOstd_rhythm(:,ch),[0.39 0.83 0.07]); 
    hold on
    shadedErrorBar(time',HbO_arrhythm(:,ch),HbOstd_arrhythm(:,ch),[0.91 0.47 0.04]);

    ylim([-0.15 0.3]) 
    box off
    set(gca,'FontSize',45)
    set(gca,'LineWidth',3)
    set(gca,'FontWeight','Bold')
    
% end
%% plot RC and AC for all channels 

HbO_rhythm = nanmean(HbOfinaltot_rhythm,3); 
HbR_rhythm = nanmean(HbRfinaltot_rhythm,3); 
HbT_rhythm = nanmean(HbTfinaltot_rhythm,3); 

HbO_arrhythm = nanmean(HbOfinaltot_arrhythm,3);
HbR_arrhythm = nanmean(HbRfinaltot_arrhythm,3);
HbT_arrhythm = nanmean(HbTfinaltot_arrhythm,3);

for i = 1:nch
    statNEW(i,:) = [nan(1,10) stat{1, i}.mask*0.29 nan(1,10)];   
    statNEW1(i,:) = [nan(1,10) statb1{1, i}.mask*(-0.1) nan(1,10)];
    statNEW2(i,:) = [nan(1,10) statb2{1, i}.mask*(-0.13) nan(1,10)];
end

statNEW(statNEW==0) = NaN;
statNEW1(statNEW1==0) = NaN;
statNEW2(statNEW2==0) = NaN;

for ch = 1:nch
    figure
    HbOstd_rhythm(:,ch) = std(squeeze(HbOfinaltot_rhythm(:,ch,:)),[],2,'omitnan')/sqrt(sum(~isnan(HbOfinaltot_rhythm(1,ch,:))));
    numSUB_avg(ch) = sum(~isnan(HbOfinaltot_rhythm(1,ch,:))); % avg over this subjects for channel ch
    HbOstd_arrhythm(:,ch) = std(squeeze(HbOfinaltot_arrhythm(:,ch,:)),[],2,'omitnan')/sqrt(sum(~isnan(HbOfinaltot_arrhythm(1,ch,:))));
    numSUB_avg2(ch) = sum(~isnan(HbOfinaltot_arrhythm(1,ch,:)));

    if (numSUB_avg(1,ch)<4)
        HbO_rhythm(:,ch) = nan;
        HbO_arrhythm(:,ch) = nan;
        statNEW(ch,:) = nan;
    end

    if (dist(1,ch)>40)
        HbO_rhythm(:,ch) = nan;
        HbO_arrhythm(:,ch) = nan;
        statNEW(ch,:) = nan;
    end

    shadedErrorBar(time',HbO_rhythm(:,ch),HbOstd_rhythm(:,ch),[0.39 0.83 0.07]); 
    hold on
    shadedErrorBar(time',HbO_arrhythm(:,ch),HbOstd_arrhythm(:,ch),[0.91 0.47 0.04]);
    hold on
    plot(time',statNEW(ch,:),'linewidth',12,'color','k')
    hold on
    plot(time',statNEW1(ch,:),'linewidth',12,'color',[0.39 0.83 0.07])
    hold on
    plot(time',statNEW2(ch,:),'linewidth',12,'color',[0.91 0.47 0.04])
    fulltitle =  ['channel ',num2str(ch)];
    title(fulltitle)
    ylim([-0.15 0.3]) 

    box off
    set(gca,'FontSize',45)
    set(gca,'LineWidth',5)
    set(gca,'FontWeight','Bold')
    xticklabels({})
    yticklabels({})

end