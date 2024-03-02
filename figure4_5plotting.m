%% CpGtoLPS screps 
clc
clear
colormap parula
load("./single_cell_features/fulltraces_CpGLPSscrep.mat")
rng(3)

lig1 = 1; %100 nM CpG
lig2 = 2; %1 ng LPS
feature = 3; %Plot AUC
ngroups = 5; %# of groups for kmeans

conditions = fulltraces_CpGLPSscrep( :, 14);
conditions = cell2mat(conditions);
conditions = conditions(cell2mat( fulltraces_CpGLPSscrep(conditions, 10) ) == lig1 & cell2mat( fulltraces_CpGLPSscrep(conditions, 11)) == lig2);

temp1 = [];
maxnorm1 = [];
clust1 = [];
clust2 = [];

for aa = 1:length(conditions)
    clust1 = [clust1; cell2mat( fulltraces_CpGLPSscrep(cell2mat( fulltraces_CpGLPSscrep(:, 14) ) == conditions(aa), [2, 3,  4])) ];
    clust2 = [clust2; cell2mat( fulltraces_CpGLPSscrep(cell2mat( fulltraces_CpGLPSscrep(:, 14) ) == conditions(aa), [6, 7, 8])) ];
    temp1 = [temp1; fulltraces_CpGLPSscrep{cell2mat( fulltraces_CpGLPSscrep(:, 14) ) == conditions(aa), 1} ];
    maxnorm1 = [maxnorm1; fulltraces_CpGLPSscrep{cell2mat( fulltraces_CpGLPSscrep(:, 14) ) == conditions(aa), feature} ];

end
maxnorm1 = mean(maxnorm1);

opts = statset('Display','final');
groups1 = kmeans(clust1(:, [1 3]),ngroups, 'Distance','cityblock', 'Replicates',5,'Options',opts);

%%%Scatterplot of peak amplitude (y axis) and time (x axis)%%%
figure(1)
clf
colors = parula(ngroups+1);

colors1 = ...
[128/244, 0, 0;
 192/244, 57/244, 43/244;
 244/244, 67/244, 54/244;
 229/244, 115/244, 115/244;
 244/244, 163/244, 187/244];
plotter = [4, 1, 2, 5, 3];
groupcol = [];
for aa = 1:length(groups1)
    groupcol(aa, :) = colors1(plotter == groups1(aa), :);
end

hold on
scatter(clust1(:, 3), clust1(:,1), 3, groupcol, 'filled', 'jitter', 'on', 'jitterAmount', .5)
set(gca, 'Xlim', [0 40])
    xticks(0:5:40);
    xticklabels(0:30:240);
hold off

%%%Heatmaps based on group%%%
figure(2) %plotting heatmap
clf
plotter = [4, 1, 2, 5, 3]; %FOR PK AND TIME 
for aa = 1:ngroups
    counter = plotter(aa);
    temptr = temp1(groups1==counter, :);
    [nr, ~] = size(temptr);
    subplot(1, ngroups, aa)
    nrused = min(nr, 100);
    rows = datasample(1:nr, nrused, 'Replace',false);
    maxes = max(temptr(rows, 1:40), [], 2);
    [~, idx] = sort(maxes, 'descend');
    imagesc(temptr(rows(idx), :), [0 5])
    %axis off;
    nr

end

%%%Violin plots based on group for first interval%%%
figure(3) %plotting group features
clf
plotter = [4, 1, 2, 5, 3]; %FOR PK AND TIME
for aa = 1:ngroups
    counter = plotter(aa);
    hold on
    cond = cell2mat( fulltraces_CpGLPSscrep(cell2mat( fulltraces_CpGLPSscrep(:, 14) ) == conditions(1), 11:13) );
        temp1 = [];
        temp2 = [];
        for bb = 1:length(conditions)
            temp1 = [temp1; fulltraces_CpGLPSscrep{cell2mat( fulltraces_CpGLPSscrep(:, 14) ) == conditions(bb), feature} ];
        end

        temp1 =  max( temp1./maxnorm1, 0);
        Violin( temp1(groups1 == counter), aa, 'Bandwidth', 0.4,  'MedianColor', [1,1,1], 'ViolinColor', colors1(aa, :), 'ViolinAlpha', .5, 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32, ...
            'BoxWidth', 0.1, 'ShowMean', true);

    set(gca, 'XLim', [0.5, 5.5], 'YLim', [-.1, 2.5]);
    xticks([1 2 3 4 5]);

end

%%%Violin plots based on group for second interval%%%
figure(4) %plotting group features
clf
maxnorms2 = mean(cell2mat( fulltraces_CpGLPSscrep(1:4, feature+4)));
plotter = [4, 1, 2, 5, 3]; %FOR PK AND TIME  
% plotter = 1:5;
for aa = 1:ngroups
    counter = plotter(aa);
    hold on
    cond = cell2mat( fulltraces_CpGLPSscrep(cell2mat( fulltraces_CpGLPSscrep(:, 14) ) == conditions(1), 11:13) );
        
        if isequal( cond,  [2, 4, 1] )
            maxnorm2 = maxnorms2(1);
        else
            maxnorm2 = maxnorms2(2);
        end
        temp1 = [];
        temp2 = [];
        for bb = 1:length(conditions)
            temp2 = [temp2; fulltraces_CpGLPSscrep{cell2mat( fulltraces_CpGLPSscrep(:, 14) ) == conditions(bb), feature+4} ];
        end
        temp2 =  max( temp2./maxnorm2, 0);
        Violin( temp2(groups1 == counter), aa, 'Bandwidth', 0.4,  'MedianColor', [1,1,1], 'ViolinColor', colors1(aa, :), 'ViolinAlpha', .5, 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32, ...
            'BoxWidth', 0.1, 'ShowMean', true);
    set(gca, 'XLim', [0.5, 5.5], 'YLim', [-.1, 2.5]);
    xticks([1 2 3 4 5]);
end

%%%Histogram of peak time%%%
figure(5)
clf
temp = clust1(:,3);
histogram(temp,20, 'Normalization', 'probability')
set(gca, 'XLim', [0, 40])

%%%Histogram of peak amplitude%%%
figure(6)
temp = clust1(:,1);
histogram(temp,40, 'Normalization', 'probability')


%% PICtoTNF screps
clc
clear
colormap parula
load("./single_cell_features/fulltraces_PICTNFscrep.mat")
rng(8)

lig1 = 3; %1000 ng/mL PIC
lig2 = 6; %10 ng/mL TNF
feature = 3; %AUC

maxnorms2 = mean(cell2mat( fulltraces_PICTNFscrep(1:4, feature+4)));

ngroups = 5;

conditions = fulltraces_PICTNFscrep( :, 14);
conditions = cell2mat(conditions);
conditions = conditions(cell2mat( fulltraces_PICTNFscrep(conditions, 10) ) == lig1 & cell2mat( fulltraces_PICTNFscrep(conditions, 11)) == lig2);

temp1 = [];
maxnorm1 = [];
clust1 = [];
clust2 = [];


for aa = 1:length(conditions)
    clust1 = [clust1; cell2mat( fulltraces_PICTNFscrep(cell2mat( fulltraces_PICTNFscrep(:, 14) ) == conditions(aa), [2, 3, 4])) ];
    clust2 = [clust2; cell2mat( fulltraces_PICTNFscrep(cell2mat( fulltraces_PICTNFscrep(:, 14) ) == conditions(aa), [6,7, 8])) ];
    temp1 = [temp1; cell2mat( fulltraces_PICTNFscrep(cell2mat( fulltraces_PICTNFscrep(:, 14) ) == conditions(aa), 1)) ];
    maxnorm1 = [maxnorm1; fulltraces_PICTNFscrep{cell2mat( fulltraces_PICTNFscrep(:, 14) ) == conditions(aa), feature} ];


end

maxnorm1 = mean(maxnorm1);

opts = statset('Display','final');
groups1 = kmeans( clust1(:, [1 3]) ,ngroups, 'Distance','cityblock', 'Replicates',5,'Options',opts);
groups1old = groups1;

%%%Scatterplot of peak amplitude (y axis) and time (x axis)%%%
figure(1)
clf
colors = parula(ngroups+1);
colors1 = ...
[13/244, 71/244, 161/244;
 2/244, 119/244, 189/244;
 30/244, 136/244, 229/244;
 66/244, 165/244, 244/244;
 144/244, 202/244, 244/244];
groupcol = [];
plotter = [2 4 5 1 3];
for aa = 1:length(groups1)
    groupcol(aa, :) = colors1(plotter == groups1(aa), :);
end
rng(1)
hold on
scatter(clust1(:, 3), clust1(:,1), 3, groupcol, 'filled', 'jitter', 'on', 'jitterAmount', .5)
    xticks(0:5:40);
    xticklabels(0:30:240);
hold off


%%%Heatmaps based on group%%%
figure(2) %plotting heatmap
clf
rng(22)
plotter = [2 4 5 1 3]; %FOR PK PLUS TIME %1:ngroups; 
for aa=1:ngroups
    counter = plotter(aa);
    temptr = temp1(groups1==counter, :);
    [nr, ~] = size(temptr);
    subplot(1, ngroups, aa)
    nrused = min(nr, 100);
    rows = datasample(1:nr, nrused, 'Replace',false);
    maxes = max(temptr(rows, 1:40), [], 2);
    [~, idx] = sort(maxes, 'descend');
    imagesc(temptr(rows(idx), :), [0 5])
    hold off

end

%%%Violin plots based on group for first interval%%%
figure(3) %plotting group features
clf
plotter = [2 4 5 1 3]; %FOR PK PLUS TIME 
for aa = 1:ngroups
    counter = plotter(aa);
    hold on
    maxnorm2 = maxnorms2(1);
    temp1 = [];
    temp2 = [];

        for bb = 1:length(conditions)
            temp1 = [temp1; fulltraces_PICTNFscrep{cell2mat( fulltraces_PICTNFscrep(:, 14) ) == conditions(bb), feature} ];
        end
        temp1 =  max( temp1./maxnorm1, 0);
        Violin( temp1(groups1 == counter), aa, 'Bandwidth', 0.4,  'MedianColor', [1,1,1], 'ViolinColor', colors(counter, :), 'ViolinAlpha', .5, 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32, ...
            'BoxWidth', 0.1, 'ShowMean', true);
    set(gca, 'XLim', [0.5, ngroups+.5], 'YLim', [-.1, 5]);
    xticks([1 2 3 4 5]);
    xticklabels({'1','2', '3', '4', '5'});

end

%%%Violin plots based on group for second interval%%%
figure(4) %plotting group features
clf
plotter = [2 4 5 1 3]; %FOR PK PLUS TIME 
for aa = 1:ngroups
    counter = plotter(aa);
    hold on
    maxnorm2 = maxnorms2(1);
    temp2 = [];
    for bb = 1:length(conditions)
        temp2 = [temp2; fulltraces_PICTNFscrep{cell2mat( fulltraces_PICTNFscrep(:, 14) ) == conditions(bb), feature+4} ];
    end
    temp2 =  max( temp2./maxnorm2, 0);
    Violin( temp2(groups1 == counter), aa, 'Bandwidth', 0.4,  'MedianColor', [1,1,1], 'ViolinColor', colors(counter, :), 'ViolinAlpha', .5, 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32, ...
        'BoxWidth', 0.1, 'ShowMean', true);
    set(gca, 'XLim', [0.5, ngroups+.5], 'YLim', [-.1, 5]);
    xticks([1 2 3 4 5]);
    xticklabels({'1','2', '3', '4', '5'});

end

%%%Histogram of peak time%%%
figure(5)
clf
temp = clust1(:,3);
histogram(temp,40, 'Normalization', 'probability')

%%%Histogram of peak amplitude%%%
figure(6)
temp = clust1(:,1);
histogram(temp,40, 'Normalization', 'probability')


%% TNFblock PIC features
clc
clear
load("./single_cell_features/peakfeatures_PICfirst_TNFblock.mat")

%Feature 2-5 = Stim A max, width, AUC, and late AUC
%Feature 7-10 = Stim B max, AUC, early AUC, and late AUC
feature = 8;

%6 and 4 for TNF, 7 and 5 for LPS
conditions = peakfeatures_PICfirst_TNFblock( ([peakfeatures_PICfirst_TNFblock{:, 12}] == 6 )', 11:15);
conditions = [conditions; peakfeatures_PICfirst_TNFblock( ([peakfeatures_PICfirst_TNFblock{:, 12}] == 4 )', 11:15)];
conditions = cell2mat(conditions);
%For PIC
inputs = [17, 19, 17, 19, 17, 19; 0, 1, 0, 1, 0, 1; 2, 3, 2, 3, 0, 1; 4, 5, 4, 5, 0, 1 ];


maxnorms = [mean(cell2mat( peakfeatures_PICfirst_TNFblock(ismember( cell2mat( peakfeatures_PICfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 17 & conditions(:, 4) == 1, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_PICfirst_TNFblock(ismember( cell2mat( peakfeatures_PICfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 19 & conditions(:, 4) == 1, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_PICfirst_TNFblock(ismember( cell2mat( peakfeatures_PICfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 17 & conditions(:, 4) == 2, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_PICfirst_TNFblock(ismember( cell2mat( peakfeatures_PICfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 19 & conditions(:, 4) == 2, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_PICfirst_TNFblock(ismember( cell2mat( peakfeatures_PICfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 17 & conditions(:, 4) == 3, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_PICfirst_TNFblock(ismember( cell2mat( peakfeatures_PICfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 19 & conditions(:, 4) == 3, 5) ), feature) ))];

figure(3) %plotting TNFblock effects
alltraces = {};
clf
%For PIC
cond = conditions(ismember( conditions(:, 1), inputs(3, :)), :);
for bb = 1:2
    temp = [];
    for cc = 1:3
        temp1 = cell2mat(  peakfeatures_PICfirst_TNFblock( ismember( cell2mat( peakfeatures_PICfirst_TNFblock(:, 15) ), cond(cond(:, 1) == inputs(3, bb+(cc-1)*2) & cond(:, 4) == cc, 5) ), feature) ) ;
        %Norm to no block PIC
        temp1 =  temp1./mean(cell2mat( peakfeatures_PICfirst_TNFblock(ismember( cell2mat( peakfeatures_PICfirst_TNFblock(:, 15) ), cond(cond(:, 1) == inputs(3, 1+(cc-1)*2) & cond(:, 4) == cc, 5) ), feature) ));
        temp = [temp;temp1];
    end
    length(temp)
    Violin( max( temp, 0), bb, 'Bandwidth', 0.3,  'MedianColor', [1,1,1], 'ViolinColor', [0 0 0], 'ViolinAlpha', .3, 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32, ...
        'BoxWidth', 0.1, 'ShowMean', true);
    alltraces{2, bb} = max( temp, 0);
end
set(gca, 'XLim', [0.25, 2.75], 'YLim', [-0.5, 4])
title( inputs(2, 1) );
xticks([1 2])
xticklabels({'Con','Blk'})
ylabel('S2 LPS AUC')

%% TNFblock PIC tracemap
rng(12)
load("./single_cell_features/peakfeatures_PICfirst_TNFblock.mat")
feature = 1;
conditions = peakfeatures_PICfirst_TNFblock( ([peakfeatures_PICfirst_TNFblock{:, 12}] == 6 )', 11:15);
conditions = [conditions; peakfeatures_PICfirst_TNFblock( ([peakfeatures_PICfirst_TNFblock{:, 12}] == 4 )', 11:15)];
conditions = cell2mat(conditions);
inputs = [17, 19, 17, 19, 17, 19; 0, 1, 0, 1, 0, 1; 2, 3, 2, 3, 0, 1; 4, 5, 4, 5, 0, 1, ];

figure(3) %plotting TNFblock effects
ylim1 = [0 4];

clf
counter = 1;

cond = conditions(ismember( conditions(:, 1), inputs(3, :)), :);
for bb = 1:2
    temp = [];
    tempt = [];
    for cc = 1:3
        temp1 = cell2mat(  peakfeatures_PICfirst_TNFblock( ismember( cell2mat( peakfeatures_PICfirst_TNFblock(:, 15) ), cond(cond(:, 1) == inputs(3, bb+(cc-1)*2) & cond(:, 4) == cc, 5) ), feature) ) ;
        tempt1 = cell2mat(  peakfeatures_PICfirst_TNFblock( ismember( cell2mat( peakfeatures_PICfirst_TNFblock(:, 15) ), cond(cond(:, 1) == inputs(3, bb+(cc-1)*2) & cond(:, 4) == cc, 5) ), 3) );
        temp = [temp;temp1];
        tempt = [tempt;tempt1];
    end
    [nr, ~] = size(temp);
    rows = datasample(1:nr, 100, 'replace', false);
    traces = temp(rows, 1:40);
    tracest = tempt(rows);

    [tracest,sortIdx]  = sort(tracest,'ascend');
    subplot(1,2, counter)

    imagesc(traces(sortIdx, :), ylim1)
    length(temp)
    counter = counter + 1;
    axis off
    if bb == 1
        title("Control")
    else
        title("TNF block")
    end

end


%% TNFblock CpG features
clc
clear
load("./single_cell_features/peakfeatures_CPGfirst_TNFblock.mat")

%Feature 2-5 = Stim A max, width, AUC, and late AUC
%Feature 7-10 = Stim B max, AUC, early AUC, and late AUC
feature = 8;

%6 and 4 for TNF, 7 and 5 for LPS
conditions = peakfeatures_CPGfirst_TNFblock( ([peakfeatures_CPGfirst_TNFblock{:, 12}] == 6 )', 11:15);
conditions = [conditions; peakfeatures_CPGfirst_TNFblock( ([peakfeatures_CPGfirst_TNFblock{:, 12}] == 4 )', 11:15)];
conditions = cell2mat(conditions);

%For CpG
inputs = [17, 19, 17, 19; 2, 3, 3, 4];
maxnorms = [mean(cell2mat( peakfeatures_CPGfirst_TNFblock(ismember( cell2mat( peakfeatures_CPGfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 17 & conditions(:, 4) == 1, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_CPGfirst_TNFblock(ismember( cell2mat( peakfeatures_CPGfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 19 & conditions(:, 4) == 1, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_CPGfirst_TNFblock(ismember( cell2mat( peakfeatures_CPGfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 17 & conditions(:, 4) == 2, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_CPGfirst_TNFblock(ismember( cell2mat( peakfeatures_CPGfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 19 & conditions(:, 4) == 2, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_CPGfirst_TNFblock(ismember( cell2mat( peakfeatures_CPGfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 17 & conditions(:, 4) == 3, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_CPGfirst_TNFblock(ismember( cell2mat( peakfeatures_CPGfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 19 & conditions(:, 4) == 3, 5) ), feature) ))];

figure(3) %plotting TNFblock effects
clf
hold on
cond = conditions(ismember( conditions(:, 1), inputs(2, :)), :);
for bb = 1:2
    temp = [];
    for cc = 1:2
        temp1 = cell2mat(  peakfeatures_CPGfirst_TNFblock( ismember( cell2mat( peakfeatures_CPGfirst_TNFblock(:, 15) ), cond(cond(:, 1) == inputs(2, bb+(cc-1)*2) & cond(:, 4) == cc, 5) ), feature) ) ;
        %Norm to no block PIC
        temp1 =  temp1./mean(cell2mat( peakfeatures_CPGfirst_TNFblock(ismember( cell2mat( peakfeatures_CPGfirst_TNFblock(:, 15) ), cond(cond(:, 1) == inputs(2, 1+(cc-1)*2) & cond(:, 4) == cc, 5) ), feature) ));
        temp = [temp;temp1];
    end
    length(temp)
    if bb == 1
        tempcon = temp;
    else
        tempblk = temp;
    end
    Violin( max( temp, 0), bb, 'Bandwidth', 0.3,  'MedianColor', [1,1,1], 'ViolinColor', [0 0 0], 'ViolinAlpha', .3, 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32, ...
        'BoxWidth', 0.1, 'ShowMean', true);
    alltraces{2, bb} = max( temp, 0);
end
set(gca, 'XLim', [0.25, 2.75], 'YLim', [-0.5, 4])
title( inputs(2, 1) );
xticks([1 2])
xticklabels({'Con','Blk'})
ylabel('S2 LPS AUC')

%% TNFblock CpG heatmap
clc
clear
rng(12)
load("./single_cell_features/peakfeatures_CPGfirst_TNFblock.mat")

feature = 1;
conditions = peakfeatures_CPGfirst_TNFblock( ([peakfeatures_CPGfirst_TNFblock{:, 12}] == 6 )', 11:15);
conditions = [conditions; peakfeatures_CPGfirst_TNFblock( ([peakfeatures_CPGfirst_TNFblock{:, 12}] == 4 )', 11:15)];
conditions = cell2mat(conditions);
%For CpG
inputs = [17, 19, 17, 19; 2, 3, 3, 4];

figure(3) %plotting TNFblock effects
ylim1 = [0 4];
clf
counter = 1;
cond = conditions(ismember( conditions(:, 1), inputs(2, :)), :);
for bb = 1:2
    temp = [];
    tempt = [];
    for cc = 1:2
        temp1 = cell2mat(  peakfeatures_CPGfirst_TNFblock( ismember( cell2mat( peakfeatures_CPGfirst_TNFblock(:, 15) ), cond(cond(:, 1) == inputs(2, bb+(cc-1)*2) & cond(:, 4) == cc, 5) ), feature) ) ;
        tempt1 = cell2mat(  peakfeatures_CPGfirst_TNFblock( ismember( cell2mat( peakfeatures_CPGfirst_TNFblock(:, 15) ), cond(cond(:, 1) == inputs(2, bb+(cc-1)*2) & cond(:, 4) == cc, 5) ), 3) );
        temp = [temp;temp1];
        tempt = [tempt;tempt1];
    end
    [nr, ~] = size(temp);
    rows = datasample(1:nr, 100, 'replace', false);
    traces = temp(rows, 1:40);
    tracest = tempt(rows);

    [tracest,sortIdx]  = sort(tracest,'ascend');
    subplot(1,2, counter)
    imagesc(traces(sortIdx, :), ylim1)
    length(temp)
    counter = counter + 1;
    axis off
    if bb == 1
        title("Control")
    else
        title("TNF block")
    end

end

%% IL10block CpG features
clc
clear

load("./single_cell_features/peakfeatures_CPGfirst_IL10block.mat")

%Feature 2-5 = Stim A max, width, AUC, and late AUC
%Feature 7-10 = Stim B max, AUC, early AUC, and late AUC
feature = 8;

%6 and 4 for TNF, 7 and 5 for LPS
conditions = peakfeatures_CPGfirst_IL10block( ([peakfeatures_CPGfirst_IL10block{:, 12}] == 6 )', 11:15);
conditions = [conditions; peakfeatures_CPGfirst_IL10block( ([peakfeatures_CPGfirst_IL10block{:, 12}] == 4 )', 11:15)];
conditions = cell2mat(conditions);
%For CpG
inputs = [17, 19, 17, 19; 3, 21, 2, 7];

maxnorms = [mean(cell2mat( peakfeatures_CPGfirst_IL10block(ismember( cell2mat( peakfeatures_CPGfirst_IL10block(:, 15) ), conditions(conditions(:, 1) == 17 & conditions(:, 4) == 1, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_CPGfirst_IL10block(ismember( cell2mat( peakfeatures_CPGfirst_IL10block(:, 15) ), conditions(conditions(:, 1) == 19 & conditions(:, 4) == 1, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_CPGfirst_IL10block(ismember( cell2mat( peakfeatures_CPGfirst_IL10block(:, 15) ), conditions(conditions(:, 1) == 17 & conditions(:, 4) == 2, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_CPGfirst_IL10block(ismember( cell2mat( peakfeatures_CPGfirst_IL10block(:, 15) ), conditions(conditions(:, 1) == 19 & conditions(:, 4) == 2, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_CPGfirst_IL10block(ismember( cell2mat( peakfeatures_CPGfirst_IL10block(:, 15) ), conditions(conditions(:, 1) == 17 & conditions(:, 4) == 3, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_CPGfirst_IL10block(ismember( cell2mat( peakfeatures_CPGfirst_IL10block(:, 15) ), conditions(conditions(:, 1) == 19 & conditions(:, 4) == 3, 5) ), feature) ))];

figure(3) %plotting TNFblock effects
clf
%For PIC
    hold on
    cond = conditions(ismember( conditions(:, 1), inputs(2, :)), :);
        for bb = 1:2
            temp = [];
            for cc = 1:2
                temp1 = cell2mat(  peakfeatures_CPGfirst_IL10block( ismember( cell2mat( peakfeatures_CPGfirst_IL10block(:, 15) ), cond(cond(:, 1) == inputs(2, bb+(cc-1)*2) & cond(:, 4) == cc, 5) ), feature) ) ;

                %Norm to no block PIC
                temp1 =  temp1./mean(cell2mat( peakfeatures_CPGfirst_IL10block(ismember( cell2mat( peakfeatures_CPGfirst_IL10block(:, 15) ), cond(cond(:, 1) == inputs(2, 1+(cc-1)*2) & cond(:, 4) == cc, 5) ), feature) ));
                temp = [temp;temp1];
            end
            length(temp)
            if bb == 1
                tempcon = temp;
            else
                tempblk = temp;
            end
            Violin( max( temp, 0), bb, 'Bandwidth', 0.3,  'MedianColor', [1,1,1], 'ViolinColor', [0 0 0], 'ViolinAlpha', .3, 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32, ...
                'BoxWidth', 0.1, 'ShowMean', true);
            alltraces{2, bb} = max( temp, 0);
        end
        set(gca, 'XLim', [0.25, 2.75], 'YLim', [-0.5, 4])
        xticks([1 2])
        xticklabels({'Con','Blk'})
        ylabel('S2 TNF AUC')

