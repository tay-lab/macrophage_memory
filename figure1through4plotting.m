%% PIC single cell heatmap for F1
%TITLE IN N1 N2 N3 format. N1 = ligand 1, N2 = ligand 2, N3 = time
clear
load("./single_cell_features/peakfeatures_firstPICCPG.mat")
rng(1)
conditions = peakfeatures_PICCPGfirst( ([peakfeatures_PICCPGfirst{:, 11}] < 3 )', 15);
conditions = cell2mat(conditions);
conditions = conditions(~ismember(conditions, [29, 31, 61, 63] ));
conditions = sort([conditions; (1:4)'; ((1:4)+32)']);
position = [9, 1, 25, 17, 11, 3, 27, 19, 13, 5, 29, 21, 15, 7, 31, 23];
position = [position, position+32];
xlim = [0, 234];

figure(1)
conditions = conditions(conditions < 33);
ylim1 = [0 1.5];
clf
normfactors = [];
for aa=1:length(conditions)  
    temp1 = [peakfeatures_PICCPGfirst{conditions(aa), 1}; peakfeatures_PICCPGfirst{conditions(aa)+32, 1}] ;
    maxPICs = cell2mat( peakfeatures_PICCPGfirst( cell2mat( peakfeatures_PICCPGfirst(:, 11) ) == 2, 2) );
    normfactors(1) = mean(maxPICs);
    if peakfeatures_PICCPGfirst{conditions(aa), 12} == 6 && peakfeatures_PICCPGfirst{conditions(aa), 13} == 4
        normfactors(2) = mean([ peakfeatures_PICCPGfirst{conditions(2), 7}; peakfeatures_PICCPGfirst{conditions(2)+32, 7} ]);
    elseif peakfeatures_PICCPGfirst{conditions(aa), 12} == 6 && peakfeatures_PICCPGfirst{conditions(aa), 13} == 8
        normfactors(2) = mean([ peakfeatures_PICCPGfirst{conditions(1), 7}; peakfeatures_PICCPGfirst{conditions(1)+32, 7} ]);
    elseif peakfeatures_PICCPGfirst{conditions(aa), 12} == 7 && peakfeatures_PICCPGfirst{conditions(aa), 13} == 4
        normfactors(2) = mean([ peakfeatures_PICCPGfirst{conditions(4), 7}; peakfeatures_PICCPGfirst{conditions(4)+32, 7} ]);
    else
        normfactors(2) = mean([ peakfeatures_PICCPGfirst{conditions(3), 7}; peakfeatures_PICCPGfirst{conditions(3)+32, 7} ]);
    end
    [nr, ~] = size(temp1);
    subplot(4, 8, position(aa) )
    rows = datasample(1:nr, 100, "replace", false);
    traces = temp1(rows, 1:40)./normfactors(1);
    maxs = max(traces, [], 2);
    [maxs,sortIdx]  = sort(maxs,'descend');
    imagesc(traces(sortIdx, :), ylim1)
    foo = string( cell2mat( peakfeatures_PICCPGfirst(conditions(aa), 11:13) ) );
    title(strjoin(foo ) )  
    axis off
    temp2 = [peakfeatures_PICCPGfirst{conditions(aa), 6}; peakfeatures_PICCPGfirst{conditions(aa)+32, 6}] ;
    [nr, ~] = size(temp2);
    rows = datasample(1:nr, 100, "replace", false);
    subplot(4, 8, position(aa)+1)
    traces = temp2(rows, 1:40)./normfactors(2);
    maxs = max(traces, [], 2);
    [maxs,sortIdx]  = sort(maxs,'descend');
    imagesc(traces(sortIdx, :), ylim1)
    foo = string( cell2mat( peakfeatures_PICCPGfirst(conditions(aa), 11:13) ) );
    title(strjoin(foo ) )  
    axis off
end
hold off

%% PIC features for F2
%Subplot1 4 hour PIC to TNF, Subplot2 8 hours PIC to TNF, Subplot3 4 hours PIC to LPS, subplot4 8 hours PIC to LPS 
clc
clear
load("./single_cell_features/peakfeatures_firstPICCPG.mat")

%Feature 2-5 = Stim A max, width, AUC, and late AUC
%Feature 7-10 = Stim B max, AUC, early AUC, and late AUC
feature = 10; 

conditions = peakfeatures_PICCPGfirst( ([peakfeatures_PICCPGfirst{:, 11}] < 3 )', 15);
conditions = cell2mat(conditions);
conditions = conditions(~ismember(conditions, [29, 31, 61, 63] ));
conditions = sort([conditions', (1:4), ((1:4)+32)]);
positions = [5, 1, 13, 9, 6, 2, 14, 10, 7, 3, 15, 11, 8, 4, 16, 12];
positions = [positions, positions+16];
conditions(positions) = conditions;
maxnorms = [mean(peakfeatures_PICCPGfirst{1, feature}), mean(peakfeatures_PICCPGfirst{2, feature}), ...
    mean(peakfeatures_PICCPGfirst{3, feature}), mean(peakfeatures_PICCPGfirst{4, feature}),...
    mean(peakfeatures_PICCPGfirst{33, feature}), mean(peakfeatures_PICCPGfirst{34, feature}), ...
    mean(peakfeatures_PICCPGfirst{35, feature}), mean(peakfeatures_PICCPGfirst{36, feature})];

figure(2) %plotting PIC effects
clf
counter = 1;
for aa = 1:4
    subplot(4,1, aa)
    y = [];
    hold on
    for bb = 1:4
        cond1 = cell2mat( peakfeatures_PICCPGfirst(conditions(counter), 12:14) );
        cond2 = cell2mat( peakfeatures_PICCPGfirst(conditions(counter)+32, 12:14) );
        if isequal( cond1,  [6, 8, 1] )
            maxnorm1 = maxnorms(1);
            maxnorm2 = maxnorms(5);
        elseif isequal( cond1,  [6, 4, 1] )
            maxnorm1 = maxnorms(2);
            maxnorm2 = maxnorms(6);
        elseif isequal( cond1,  [7, 8, 1] )
            maxnorm1 = maxnorms(3);
            maxnorm2 = maxnorms(7);
        else
            maxnorm1 = maxnorms(4);
            maxnorm2 = maxnorms(8);
        end
        temp1 = peakfeatures_PICCPGfirst{conditions(counter), feature}./maxnorm1 ;
        temp2 = peakfeatures_PICCPGfirst{conditions(counter)+32, feature}./maxnorm2 ;
        Violin( [temp1; temp2], bb/2, 'Bandwidth', 0.3,  'MedianColor', [1,1,1], 'ViolinColor', [0 0 0], 'ViolinAlpha', .3, 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32/2, ...
            'BoxWidth', 0.1, 'ShowMean', true);
        y = [y, nanmean([temp1;temp2])];
        counter = counter + 1;
        controldist = [peakfeatures_PICCPGfirst{conditions((aa-1)*4+1), feature}./maxnorm1; peakfeatures_PICCPGfirst{conditions((aa-1)*4+1)+32, feature}./maxnorm2];
        p = ranksum(controldist, [temp1; temp2])*3
        alltraces{aa, bb} = [temp1; temp2];
    end
    x=0:0.5:1;
    P = polyfit(x,y(2:4),1);
    y = polyval(P, -.5/2:.1:2.5/2);
    plot(1.5/2:.1:4.5/2,y,'r-', 'LineWidth',2);
    if ismember(aa, 1:2)
        text(1/2, 4.75, string(round(P(1), 2)+ "*log(x) + " + round(P(2), 2) ),'fontsize',12)
        set(gca, 'XLim', [0.25, 4.75/2], 'YLim', [-0.5, 5])
    else
        text(1/2, 3.25, string(round(P(1), 2)+ "*log(x) + " + round(P(2), 2) ),'fontsize',12)
        set(gca, 'XLim', [0.25, 4.75/2], 'YLim', [-0.5, 3.5])
    end
    title( strjoin( string( cell2mat( peakfeatures_PICCPGfirst(conditions(counter-1), 12:13) ) ) ) );
    xticks([1/2 2/2 3/2 4/2])
    xticklabels({'Con','0.1 ug','0.3 ug','1 ug'})
    set(gca,'XTickLabel',{'0 ng/mL','100 ng/mL','300 ng/mL','1000 ng/mL'},'fontsize',12)
    set(gca,'fontsize',12)
end

%% CPG single cell heatmap for F1
%TITLE IN N1 N2 N3 format. N1 = ligand 1, N2 = ligand 2, N3 = time
clear
load("./single_cell_features/peakfeatures_firstPICCPG.mat")
rng(5)
conditions = peakfeatures_PICCPGfirst( ([peakfeatures_PICCPGfirst{:, 11}] > 2 )', 15);
conditions = cell2mat(conditions);
conditions = conditions(~ismember(conditions, [30, 32, 62, 64] ));
position = [9, 1, 25, 17, 11, 3, 27, 19, 13, 5, 29, 21, 15, 7, 31, 23];
position = [position, position+32];
xlim = [0, 234];
figure(2) %plotting PIC effects
clf
conditions = conditions(conditions < 33);
ylim1 = [0 1.5];
normfactors = [];
for aa=1:length(conditions)  
    temp1 = [peakfeatures_PICCPGfirst{conditions(aa), 1}; peakfeatures_PICCPGfirst{conditions(aa)+32, 1}] ;
    maxCPGs = cell2mat( peakfeatures_PICCPGfirst( cell2mat( peakfeatures_PICCPGfirst(:, 11) ) == 5, 2) );
    normfactors(1) = mean(maxCPGs);
    if peakfeatures_PICCPGfirst{conditions(aa), 12} == 6 && peakfeatures_PICCPGfirst{conditions(aa), 13} == 4
        normfactors(2) = mean([ peakfeatures_PICCPGfirst{conditions(2), 7}; peakfeatures_PICCPGfirst{conditions(2)+32, 7} ]);
    elseif peakfeatures_PICCPGfirst{conditions(aa), 12} == 6 && peakfeatures_PICCPGfirst{conditions(aa), 13} == 8
        normfactors(2) = mean([ peakfeatures_PICCPGfirst{conditions(1), 7}; peakfeatures_PICCPGfirst{conditions(1)+32, 7} ]);
    elseif peakfeatures_PICCPGfirst{conditions(aa), 12} == 7 && peakfeatures_PICCPGfirst{conditions(aa), 13} == 4
        normfactors(2) = mean([ peakfeatures_PICCPGfirst{conditions(4), 7}; peakfeatures_PICCPGfirst{conditions(4)+32, 7} ]);
    else
        normfactors(2) = mean([ peakfeatures_PICCPGfirst{conditions(3), 7}; peakfeatures_PICCPGfirst{conditions(3)+32, 7} ]);
    end
    [nr, ~] = size(temp1);
    subplot(4, 8, position(aa) )
    rows = datasample(1:nr, 100, "replace", false);
    traces = temp1(rows, 1:40)./normfactors(1);
    maxs = max(traces, [], 2);
    [maxs,sortIdx]  = sort(maxs,'descend');
    imagesc(traces(sortIdx, :), ylim1)
    foo = string( cell2mat( peakfeatures_PICCPGfirst(conditions(aa), 11:14) ) );
    title(strjoin(foo ) )  
    axis off
    temp2 = [peakfeatures_PICCPGfirst{conditions(aa), 6}; peakfeatures_PICCPGfirst{conditions(aa)+32, 6}] ;
    [nr, ~] = size(temp2);
    rows = datasample(1:nr, 100, "replace", false);
    subplot(4, 8, position(aa)+1)
    traces = temp2(rows, 1:40)./normfactors(2);
    maxs = max(traces, [], 2);
    [maxs,sortIdx]  = sort(maxs,'descend');
    imagesc(traces(sortIdx, :), ylim1)
    foo = string( cell2mat( peakfeatures_PICCPGfirst(conditions(aa), 11:13) ) );
    title(strjoin(foo ) )  
    axis off
end
hold off

%% CPG features for F2
%Subplot1 4 hour CpG to TNF, Subplot2 8 hours CpG to TNF, Subplot3 4 hours CpG to LPS, subplot4 8 hours CpG to LPS 
clc
clear
load("./single_cell_features/peakfeatures_firstPICCPG.mat"),

%Feature 2-5 = Stim A max, width, AUC, and late AUC
%Feature 7-10 = Stim B max, AUC, early AUC, and late AUC
feature = 7;

conditions = peakfeatures_PICCPGfirst( ([peakfeatures_PICCPGfirst{:, 11}] > 2 )', 15);
conditions = cell2mat(conditions);
conditions = conditions(~ismember(conditions, [30, 32, 62, 64] ));
positions = [5, 1, 13, 9, 6, 2, 14, 10, 7, 3, 15, 11, 8, 4, 16, 12];
positions = [positions, positions+16];
conditions(positions) = conditions;
maxnorms = [mean(peakfeatures_PICCPGfirst{1, feature}), mean(peakfeatures_PICCPGfirst{2, feature}), ...
    mean(peakfeatures_PICCPGfirst{3, feature}), mean(peakfeatures_PICCPGfirst{4, feature}),...
    mean(peakfeatures_PICCPGfirst{33, feature}), mean(peakfeatures_PICCPGfirst{34, feature}), ...
    mean(peakfeatures_PICCPGfirst{35, feature}), mean(peakfeatures_PICCPGfirst{36, feature})];
figure(4) %plotting CPG effects
clf
counter = 1;
for aa = 1:4
    y = [];
    subplot(4,1, aa)
    hold on
    for bb = 1:4
        cond1 = cell2mat( peakfeatures_PICCPGfirst(conditions(counter), 12:14) );
        cond2 = cell2mat( peakfeatures_PICCPGfirst(conditions(counter)+32, 12:14) );
        if isequal( cond1,  [6, 8, 1] )
            maxnorm1 = maxnorms(1);
            maxnorm2 = maxnorms(5);
        elseif isequal( cond1,  [6, 4, 1] )
            maxnorm1 = maxnorms(2);
            maxnorm2 = maxnorms(6);
        elseif isequal( cond1,  [7, 8, 1] )
            maxnorm1 = maxnorms(3);
            maxnorm2 = maxnorms(7);
        else
            maxnorm1 = maxnorms(4);
            maxnorm2 = maxnorms(8);
        end
        temp1 = peakfeatures_PICCPGfirst{conditions(counter), feature}./maxnorm1 ;
        temp2 = peakfeatures_PICCPGfirst{conditions(counter)+32, feature}./maxnorm2 ;
        Violin( [temp1; temp2], bb/2, 'Bandwidth', 0.3,  'MedianColor', [1,1,1], 'ViolinColor', [0 0 0], 'ViolinAlpha', .3, 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32/2, ...
            'BoxWidth', 0.1, 'ShowMean', true);
        y = [y, nanmean([temp1;temp2])];
        counter = counter + 1;
        controldist = [peakfeatures_PICCPGfirst{conditions((aa-1)*4+1), feature}./maxnorm1; peakfeatures_PICCPGfirst{conditions((aa-1)*4+1)+32, feature}./maxnorm2];
        p = ranksum(controldist, [temp1; temp2])*3
        alltraces{aa, bb} = [temp1; temp2];
    end
    x=0:0.5:1;
    P = polyfit(x,y(2:4),1);
    y = polyval(P, -.5/2:.1:2.5/2);
    plot(1.5/2:.1:4.5/2,y,'r-', 'LineWidth',2);
    if ismember(aa, [1,2])
        set(gca, 'XLim', [0.25, 4.75/2], 'YLim', [-0.25, 4])
        text(1/2, 3.75, string(round(P(1), 2)+ "*log(x) + " + round(P(2), 2) ),'fontsize',12)
    else
        set(gca, 'XLim', [0.25, 4.75/2], 'YLim', [-0.25, 3])
        text(1/2, 2.75, string(round(P(1), 2)+ "*log(x) + " + round(P(2), 2) ),'fontsize',12)
    end
    title( strjoin( string( cell2mat( peakfeatures_PICCPGfirst(conditions(counter-1), 12:13) ) ) ) );
    xticks([1/2 2/2 3/2 4/2])
    xticklabels({'Con','10 nM','30 nM','100 nM'})
    set(gca,'fontsize',12)
end

%% TNF single cell heatmap for F1
%TITLE IN N1 N2 N3 format. N1 = ligand 1, N2 = ligand 2, N3 = time
clear
clf
load("./single_cell_features/peakfeatures_firstTNFLPS.mat")
rng(6)
conditions = peakfeatures_TNFLPSfirst( ([peakfeatures_TNFLPSfirst{:, 11}] == 17 )', 15);
conditions = [conditions; peakfeatures_TNFLPSfirst( ([peakfeatures_TNFLPSfirst{:, 11}] < 4 )', 15)];
conditions = cell2mat(conditions);
position = [1, 11, 21, 31, 13, 3, 33, 23, 15, 5, 35, 25, 17, 37, 27, 29, 39, 19, 7];
position = [position, position+40, 27, 66];
xlim = [0, 234];
figure(2) %plotting LPS effects
clf
conditions = [conditions(conditions < 33); 65];
ylim1 = [0 1.5];
clf
normfactors = [];
for aa=1:length(conditions) 
    if aa ~= 19
        temp1 = [peakfeatures_TNFLPSfirst{conditions(aa), 1}; peakfeatures_TNFLPSfirst{conditions(aa)+32, 1}] ;
        temp2 = [peakfeatures_TNFLPSfirst{conditions(aa), 6}; peakfeatures_TNFLPSfirst{conditions(aa)+32, 6}] ;
    else
        temp1 = [peakfeatures_TNFLPSfirst{conditions(aa), 1}; peakfeatures_TNFLPSfirst{conditions(aa)+2, 1}] ;
        temp2 = [peakfeatures_TNFLPSfirst{conditions(aa), 6}; peakfeatures_TNFLPSfirst{conditions(aa)+2, 6}] ;
    end
    maxTNFs = cell2mat( peakfeatures_TNFLPSfirst( cell2mat( peakfeatures_TNFLPSfirst(:, 11) ) == 2, 2) );
    normfactors(1) = mean(maxTNFs);
    if peakfeatures_TNFLPSfirst{conditions(aa), 12} == 2 && peakfeatures_TNFLPSfirst{conditions(aa), 13} == 4
        normfactors(2) = mean([ peakfeatures_TNFLPSfirst{conditions(2), 7}; peakfeatures_TNFLPSfirst{conditions(2)+32, 7} ]);
    elseif peakfeatures_TNFLPSfirst{conditions(aa), 12} == 2 && peakfeatures_TNFLPSfirst{conditions(aa), 13} == 8
        normfactors(2) = mean([ peakfeatures_TNFLPSfirst{conditions(1), 7}; peakfeatures_TNFLPSfirst{conditions(1)+32, 7} ]);
    elseif peakfeatures_TNFLPSfirst{conditions(aa), 12} == 6 && peakfeatures_TNFLPSfirst{conditions(aa), 13} == 4
        normfactors(2) = mean([ peakfeatures_TNFLPSfirst{conditions(4), 7}; peakfeatures_TNFLPSfirst{conditions(4)+32, 7} ]);
    else
        normfactors(2) = mean([ peakfeatures_TNFLPSfirst{conditions(3), 7}; peakfeatures_TNFLPSfirst{conditions(3)+32, 7} ]);
    end
    [nr, ~] = size(temp1);
    subplot(4, 10, position(aa) )
    rows = datasample(1:nr, 100, "replace", false);
    traces = temp1(rows, 1:40)./normfactors(1);
    maxs = max(traces, [], 2);
    [maxs,sortIdx]  = sort(maxs,'descend');
    imagesc(traces(sortIdx, :), ylim1)
    foo = string( cell2mat( peakfeatures_TNFLPSfirst(conditions(aa), 11:14) ) );
    title(strjoin(foo ) )  
    axis off
    [nr, ~] = size(temp2);
    rows = datasample(1:nr, 100, "replace", false);
    subplot(4, 10, position(aa)+1)
    traces = temp2(rows, 1:40)./normfactors(2);
    maxs = max(traces, [], 2);
    [maxs,sortIdx]  = sort(maxs,'descend');
    imagesc(traces(sortIdx, :), ylim1)
    foo = string( cell2mat( peakfeatures_TNFLPSfirst(conditions(aa), 11:13) ) );
    title(strjoin(foo ) )  
    axis off
end
hold off
%% TNF features for F2
%Subplot1 4 hour TNF to TNF, Subplot2 8 hours TNF to TNF, Subplot3 4 hours TNF to LPS, subplot4 8 hours TNF to LPS 
clc
clear
load("./single_cell_features/peakfeatures_firstTNFLPS.mat")

%Feature 2-5 = Stim A max, width, AUC, and late AUC
%Feature 7-10 = Stim B max, AUC, early AUC, and late AUC
feature = 8;

conditions = peakfeatures_TNFLPSfirst( ([peakfeatures_TNFLPSfirst{:, 11}] < 4 )', 15);
conditions = cell2mat(conditions);
conditions = sort([conditions; (1:4)'; ((1:4)+32)']);
positions = [1, 6, 11, 16, 7, 2, 17, 12, 8, 3, 18, 13, 9, 19, 14, 15, 20, 10];
positions = [positions, positions+20, 4, 24];
newconds = NaN([40, 1]);
newconds(positions) = conditions;
maxnorms = [mean(peakfeatures_TNFLPSfirst{1, feature}), mean(peakfeatures_TNFLPSfirst{2, feature}), ...
    mean(peakfeatures_TNFLPSfirst{3, feature}), mean(peakfeatures_TNFLPSfirst{4, feature}),...
    mean(peakfeatures_TNFLPSfirst{33, feature}), mean(peakfeatures_TNFLPSfirst{34, feature}), ...
    mean(peakfeatures_TNFLPSfirst{35, feature}), mean(peakfeatures_TNFLPSfirst{36, feature})];
figure(1) 
clf
counter = 1;
for aa = 1:4
    y = [];
    subplot(4,1, aa)
    hold on
    for bb = 1:5
        if ~isnan(newconds(counter))
            cond1 = cell2mat( peakfeatures_TNFLPSfirst(newconds(counter), 12:14) );
            if isequal( cond1,  [2, 4, 1] )
                maxnorm1 = maxnorms(1);
                maxnorm2 = maxnorms(5);
            elseif isequal( cond1,  [2, 8, 1] )
                maxnorm1 = maxnorms(2);
                maxnorm2 = maxnorms(6);
            elseif isequal( cond1,  [6, 4, 1] )
                maxnorm1 = maxnorms(3);
                maxnorm2 = maxnorms(7);
            else
                maxnorm1 = maxnorms(4);
                maxnorm2 = maxnorms(8);
            end
            if newconds(counter) == 65
                temp1 = peakfeatures_TNFLPSfirst{newconds(counter), feature}./maxnorm1 ;
                temp2 = peakfeatures_TNFLPSfirst{newconds(counter)+2, feature}./maxnorm2 ;
            else
                temp1 = peakfeatures_TNFLPSfirst{newconds(counter), feature}./maxnorm1 ;
                temp2 = peakfeatures_TNFLPSfirst{newconds(counter)+32, feature}./maxnorm2 ;
            end
            Violin( [temp1; temp2], bb/2, 'Bandwidth', 0.3,  'MedianColor', [1,1,1], 'ViolinColor', [0 0 0], 'ViolinAlpha', .3, 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32/2, ...
                'BoxWidth', 0.1, 'ShowMean', true);
        end
        counter = counter + 1;
        y = [y, mean([temp1;temp2])];
        controldist = [peakfeatures_TNFLPSfirst{newconds((aa-1)*5+1), feature}./maxnorm1; peakfeatures_TNFLPSfirst{newconds((aa-1)*5+1)+32, feature}./maxnorm2];
        p = ranksum(controldist, [temp1; temp2])*4
    end
    x=0:0.5:1.5;
    P = polyfit(x,y(2:5),1);
    y = polyval(P, -.5/2:.1:3.5/2);
    plot(1.5/2:.1:5.5/2,y,'r-', 'LineWidth',2);
    if ismember(aa, [1,2])
        set(gca, 'XLim', [0.25, 5.75/2], 'YLim', [-0.25, 5])
        text(1/2, 4.75, string(round(P(1), 2)+ "*log(x) + " + round(P(2), 2) ),'fontsize',12)
    else
        set(gca, 'XLim', [0.25, 5.75/2], 'YLim', [-0.25, 3.5])
        text(1/2, 3.25, string(round(P(1), 2)+ "*log(x) + " + round(P(2), 2) ),'fontsize',12)
    end
    title( strjoin( string( cell2mat( peakfeatures_TNFLPSfirst(newconds(counter-2), 12:13) ) ) ) );
    xticks([1/2 2/2 3/2 4/2 5/2])
    set(gca,'XTickLabel',{'0 ng/mL','1 ng/mL','3 ng/mL','10 ng/mL','30 ng/mL'},'fontsize',12)
    set(gca,'fontsize',12)
end

%% LPS whole cell heatmap for F1
%TITLE IN N1 N2 N3 format. N1 = ligand 1, N2 = ligand 2, N3 = time
clear
load("./single_cell_features/peakfeatures_firstTNFLPS.mat")
rng(6)
conditions = peakfeatures_TNFLPSfirst( ([peakfeatures_TNFLPSfirst{:, 11}] > 3 )', 15);
conditions = cell2mat(conditions);
position = [1, 11, 21, 31, 13, 3, 33, 23, 15, 5, 35, 25, 17, 7, 37, 19, 9, 39, 27];
position = [position, position+40, 27, 67];
xlim = [0, 234];
figure(2) %plotting LPS effects
clf
conditions = [conditions(conditions < 33); 66];
ylim1 = [0 1.5];
clf
normfactors = [];
for aa=1:length(conditions) 
    if aa ~= 19
        temp1 = [peakfeatures_TNFLPSfirst{conditions(aa), 1}; peakfeatures_TNFLPSfirst{conditions(aa)+32, 1}] ;
        temp2 = [peakfeatures_TNFLPSfirst{conditions(aa), 6}; peakfeatures_TNFLPSfirst{conditions(aa)+32, 6}] ;
    else
        temp1 = [peakfeatures_TNFLPSfirst{conditions(aa), 1}; peakfeatures_TNFLPSfirst{conditions(aa)+2, 1}] ;
        temp2 = [peakfeatures_TNFLPSfirst{conditions(aa), 6}; peakfeatures_TNFLPSfirst{conditions(aa)+2, 6}] ;
    end
    maxLPSs = cell2mat( peakfeatures_TNFLPSfirst( cell2mat( peakfeatures_TNFLPSfirst(:, 11) ) == 7, 2) );
    normfactors(1) = mean(maxLPSs);
    if peakfeatures_TNFLPSfirst{conditions(aa), 12} == 2 && peakfeatures_TNFLPSfirst{conditions(aa), 13} == 4
        normfactors(2) = mean([ peakfeatures_TNFLPSfirst{conditions(2), 7}; peakfeatures_TNFLPSfirst{conditions(2)+32, 7} ]);
    elseif peakfeatures_TNFLPSfirst{conditions(aa), 12} == 2 && peakfeatures_TNFLPSfirst{conditions(aa), 13} == 8
        normfactors(2) = mean([ peakfeatures_TNFLPSfirst{conditions(1), 7}; peakfeatures_TNFLPSfirst{conditions(1)+32, 7} ]);
    elseif peakfeatures_TNFLPSfirst{conditions(aa), 12} == 6 && peakfeatures_TNFLPSfirst{conditions(aa), 13} == 4
        normfactors(2) = mean([ peakfeatures_TNFLPSfirst{conditions(4), 7}; peakfeatures_TNFLPSfirst{conditions(4)+32, 7} ]);
    else
        normfactors(2) = mean([ peakfeatures_TNFLPSfirst{conditions(3), 7}; peakfeatures_TNFLPSfirst{conditions(3)+32, 7} ]);
    end
    [nr, ~] = size(temp1);
    subplot(4, 10, position(aa) )
    rows = datasample(1:nr, 100, "replace", false);
    traces = temp1(rows, 1:40)./normfactors(1);
    maxs = max(traces, [], 2);
    [maxs,sortIdx]  = sort(maxs,'descend');
    imagesc(traces(sortIdx, :), ylim1)
    foo = string( cell2mat( peakfeatures_TNFLPSfirst(conditions(aa), 11:14) ) );
    title(strjoin(foo ) )  
    axis off
    [nr, ~] = size(temp2);
    rows = datasample(1:nr, 100, "replace", false);
    subplot(4, 10, position(aa)+1)
    traces = temp2(rows, 1:40)./normfactors(2);
    maxs = max(traces, [], 2);
    [maxs,sortIdx]  = sort(maxs,'descend');
    imagesc(traces(sortIdx, :), ylim1)
    foo = string( cell2mat( peakfeatures_TNFLPSfirst(conditions(aa), 11:13) ) );
    title(strjoin(foo ) )  
    axis off
end
hold off

%% LPS features for F2
%Subplot1 4 hour LPS to TNF, Subplot2 8 hours LPS to TNF, Subplot3 4 hours LPS to LPS, subplot4 8 hours LPS to LPS 
clc
clear
load("./single_cell_features/peakfeatures_firstTNFLPS.mat")

%Feature 2-5 = Stim A max, width, AUC, and late AUC
%Feature 7-10 = Stim B max, AUC, early AUC, and late AUC
feature = 8;

conditions = peakfeatures_TNFLPSfirst( ([peakfeatures_TNFLPSfirst{:, 11}] > 3 )', 15);
conditions = cell2mat(conditions);
position = [1, 6, 11, 16, 7, 2, 17, 12, 8, 3, 18, 13, 9, 19, 14, 19, 9, 39];
position = [position, position+40, 27, 67];
positions = [1, 6, 11, 16, 7, 2, 17, 12, 8, 3, 18, 13, 9, 4, 19, 10, 5, 20];
positions = [positions, positions+20, 14, 34];
newconds = NaN([40, 1]);
newconds(positions) = conditions;
maxnorms = [mean(peakfeatures_TNFLPSfirst{1, feature}), mean(peakfeatures_TNFLPSfirst{2, feature}), ...
    mean(peakfeatures_TNFLPSfirst{3, feature}), mean(peakfeatures_TNFLPSfirst{4, feature}),...
    mean(peakfeatures_TNFLPSfirst{33, feature}), mean(peakfeatures_TNFLPSfirst{34, feature}), ...
    mean(peakfeatures_TNFLPSfirst{35, feature}), mean(peakfeatures_TNFLPSfirst{36, feature})];
figure(1) %plotting LPS effects
clf
counter = 1;
for aa = 1:4
    y = [];
    subplot(4,1, aa)
    hold on
    for bb = 1:5
        if ~isnan(newconds(counter))
            cond1 = cell2mat( peakfeatures_TNFLPSfirst(newconds(counter), 12:14) );
            if isequal( cond1,  [2, 4, 1] )
                maxnorm1 = maxnorms(1);
                maxnorm2 = maxnorms(5);
            elseif isequal( cond1,  [2, 8, 1] )
                maxnorm1 = maxnorms(2);
                maxnorm2 = maxnorms(6);
            elseif isequal( cond1,  [6, 4, 1] )
                maxnorm1 = maxnorms(3);
                maxnorm2 = maxnorms(7);
            else
                maxnorm1 = maxnorms(4);
                maxnorm2 = maxnorms(8);
            end
            if newconds(counter) == 66
                temp1 = peakfeatures_TNFLPSfirst{newconds(counter), feature}./maxnorm1 ;
                temp2 = peakfeatures_TNFLPSfirst{newconds(counter)+2, feature}./maxnorm2 ;
            else
                temp1 = peakfeatures_TNFLPSfirst{newconds(counter), feature}./maxnorm1 ;
                temp2 = peakfeatures_TNFLPSfirst{newconds(counter)+32, feature}./maxnorm2 ;
            end
            temp1 = max(temp1, 0);
            temp2 = max(temp2, 0);
            Violin( [temp1; temp2], bb/2, 'Bandwidth', 0.3,  'MedianColor', [1,1,1], 'ViolinColor', [0 0 0], 'ViolinAlpha', .3, 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32/2, ...
                'BoxWidth', 0.1, 'ShowMean', true);
        end
        counter = counter + 1;
        y = [y, mean([temp1;temp2])];
        controldist = [peakfeatures_TNFLPSfirst{newconds((aa-1)*5+1), feature}./maxnorm1; peakfeatures_TNFLPSfirst{newconds((aa-1)*5+1)+32, feature}./maxnorm2];
        p = ranksum(controldist, [temp1; temp2])*4
    end
    x=0:0.5:1.5;
    P = polyfit(x,y(2:5),1);
    y = polyval(P, -0.5/2:.1:3.5/2);
    plot(1.5/2:.1:5.5/2,y,'r-', 'LineWidth',2);    
    if ismember(aa, [1,2])
        set(gca, 'XLim', [0.25, 5.75/2], 'YLim', [-0.25, 5])
        text(1/2, 4.75, string(round(P(1), 2)+ "*log(x) + " + round(P(2), 2) ),'fontsize',12)
    else
        set(gca, 'XLim', [0.25, 5.75/2], 'YLim', [-0.25, 3.5])
        text(1/2, 3.25, string(round(P(1), 2)+ "*log(x) + " + round(P(2), 2) ),'fontsize',12)
    end
    title( strjoin( string( cell2mat( peakfeatures_TNFLPSfirst(newconds(counter-2), 12:13) ) ) ) );

    xticks([1/2 2/2 3/2 4/2 5/2])
    set(gca,'XTickLabel',{'0 ng/mL','0.1 ng/mL','0.3 ng/mL','1 ng/mL','3 ng/mL'},'fontsize',12)
    set(gca,'fontsize',12)
end

%% PAM single cell heatmap for F1
%TITLE IN N1 N2 N3 format. N1 = ligand 1, N2 = ligand 2, N3 = time
clear
load("./single_cell_features/peakfeatures_firstPAMIL1")
rng(1)
conditions = peakfeatures_PAMIL1first( ([peakfeatures_PAMIL1first{:, 11}] < 3 )', 15);
conditions = cell2mat(conditions);
conditions = conditions(~ismember(conditions, [29, 31, 61, 63] ));
conditions = sort([conditions; (1:4)'; ((1:4)+32)']);
position = [9, 1, 25, 17, 11, 3, 27, 19, 13, 5, 29, 21, 15, 7, 31, 23];
position = [position, position+32];
xlim = [0, 234];
figure(2) %plotting PAM effects
conditions = conditions(conditions < 33);
ylim1 = [0 1.5];
clf
normfactors = [];
for aa=1:length(conditions)  
    temp1 = [peakfeatures_PAMIL1first{conditions(aa), 1}; peakfeatures_PAMIL1first{conditions(aa)+32, 1}] ;
    maxPAMs = cell2mat( peakfeatures_PAMIL1first( cell2mat( peakfeatures_PAMIL1first(:, 11) ) == 2, 2) );
    normfactors(1) = mean(maxPAMs);
    if peakfeatures_PAMIL1first{conditions(aa), 12} == 6 && peakfeatures_PAMIL1first{conditions(aa), 13} == 4
        normfactors(2) = mean([ peakfeatures_PAMIL1first{conditions(2), 7}; peakfeatures_PAMIL1first{conditions(2)+32, 7} ]);
    elseif peakfeatures_PAMIL1first{conditions(aa), 12} == 6 && peakfeatures_PAMIL1first{conditions(aa), 13} == 8
        normfactors(2) = mean([ peakfeatures_PAMIL1first{conditions(1), 7}; peakfeatures_PAMIL1first{conditions(1)+32, 7} ]);
    elseif peakfeatures_PAMIL1first{conditions(aa), 12} == 7 && peakfeatures_PAMIL1first{conditions(aa), 13} == 4
        normfactors(2) = mean([ peakfeatures_PAMIL1first{conditions(4), 7}; peakfeatures_PAMIL1first{conditions(4)+32, 7} ]);
    else
        normfactors(2) = mean([ peakfeatures_PAMIL1first{conditions(3), 7}; peakfeatures_PAMIL1first{conditions(3)+32, 7} ]);
    end
    [nr, ~] = size(temp1);
    subplot(4, 8, position(aa) )
    rows = datasample(1:nr, 100, "replace", false);
    traces = temp1(rows, 1:40)./normfactors(1);
    maxs = max(traces, [], 2);
    [maxs,sortIdx]  = sort(maxs,'descend');
    imagesc(traces(sortIdx, :), ylim1)
    foo = string( cell2mat( peakfeatures_PAMIL1first(conditions(aa), 11:14) ) );
    title(strjoin(foo ) )  
    axis off
    temp2 = [peakfeatures_PAMIL1first{conditions(aa), 6}; peakfeatures_PAMIL1first{conditions(aa)+32, 6}] ;
    [nr, ~] = size(temp2);
    rows = datasample(1:nr, 100, "replace", false);
    subplot(4, 8, position(aa)+1)
    traces = temp2(rows, 1:40)./normfactors(2);
    maxs = max(traces, [], 2);
    [maxs,sortIdx]  = sort(maxs,'descend');
    imagesc(traces(sortIdx, :), ylim1)
    foo = string( cell2mat( peakfeatures_PAMIL1first(conditions(aa), 11:13) ) );
    title(strjoin(foo ) )  
    axis off
end
hold off

%% PAM features for F2
%Subplot1 4 hour PAM to TNF, Subplot2 8 hours PAM to TNF, Subplot3 4 hours PAM to LPS, subplot4 8 hours PAM to LPS 
clc
clear
load("./single_cell_features/peakfeatures_firstPAMIL1")

%Feature 2-5 = Stim A max, width, AUC, and late AUC
%Feature 7-10 = Stim B max, AUC, early AUC, and late AUC
feature = 7;

conditions = peakfeatures_PAMIL1first( ([peakfeatures_PAMIL1first{:, 11}] < 3 )', 15);
conditions = cell2mat(conditions);
conditions = conditions(~ismember(conditions, [28:32, 60:64] ));
conditions = sort([conditions', (1:4), ((1:4)+32)]);
positions = [5, 1, 13, 9, 6, 2, 14, 10, 7, 3, 15, 11, 8, 4, 16, 12];
positions = [positions, positions+16];
conditions(positions) = conditions;
maxnorms = [mean(peakfeatures_PAMIL1first{1, feature}), mean(peakfeatures_PAMIL1first{2, feature}), ...
    mean(peakfeatures_PAMIL1first{3, feature}), mean(peakfeatures_PAMIL1first{4, feature}) ,...
    mean(peakfeatures_PAMIL1first{33, feature}), mean(peakfeatures_PAMIL1first{34, feature}), ...
    mean(peakfeatures_PAMIL1first{35, feature}), mean(peakfeatures_PAMIL1first{36, feature})];
figure(1)
clf
counter = 1;
for aa = 1:4
    subplot(4,1, aa)
    y = [];
    hold on
    for bb = 1:4
        cond1 = cell2mat( peakfeatures_PAMIL1first(conditions(counter), 12:14) );
        cond2 = cell2mat( peakfeatures_PAMIL1first(conditions(counter)+32, 12:14) );
        if isequal( cond1,  [6, 8, 1] )
            maxnorm1 = maxnorms(1);
            maxnorm2 = maxnorms(5);
        elseif isequal( cond1,  [6, 4, 1] )
            maxnorm1 = maxnorms(2);
            maxnorm2 = maxnorms(6);
        elseif isequal( cond1,  [7, 8, 1] )
            maxnorm1 = maxnorms(3);
            maxnorm2 = maxnorms(7);
        else
            maxnorm1 = maxnorms(4);
            maxnorm2 = maxnorms(8);
        end
        temp1 = peakfeatures_PAMIL1first{conditions(counter), feature}./maxnorm1 ;
        temp1 =  temp1( temp1 > prctile(temp1, .5) & temp1 < prctile(temp1, 99.5) );
        temp2 = peakfeatures_PAMIL1first{conditions(counter)+32, feature}./maxnorm2 ;
        temp2 =  temp2( temp2 > prctile(temp2, .5) & temp2 < prctile(temp2, 99.5) );
        Violin( [temp1; temp2], bb/2, 'Bandwidth', 0.3,  'MedianColor', [1,1,1], 'ViolinColor', [0 0 0], 'ViolinAlpha', .3, 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32/2, ...
            'BoxWidth', 0.1, 'ShowMean', true);
        y = [y, mean([temp1;temp2])];
        counter = counter + 1;
        controldist = [peakfeatures_PAMIL1first{conditions((aa-1)*4+1), feature}./maxnorm1; peakfeatures_PAMIL1first{conditions((aa-1)*4+1)+32, feature}./maxnorm2];
        p = ranksum(controldist, [temp1; temp2])*3
    end
    x=0:0.5:1;
    P = polyfit(x,y(2:4),1);
    y = polyval(P, -.5/2:.1:2.5/2);
    plot(1.5/2:.1:4.5/2,y,'r-', 'LineWidth',2);

    if ismember(aa, [1,2])
        set(gca, 'XLim', [0.25, 4.75/2], 'YLim', [-0.25, 5])
        text(1/2, 4.75, string(round(P(1), 2)+ "*log(x) + " + round(P(2), 2) ),'fontsize',12)
    else
        set(gca, 'XLim', [0.25, 4.75/2], 'YLim', [-0.25, 3.5])
        text(1/2, 3.25, string(round(P(1), 2)+ "*log(x) + " + round(P(2), 2) ),'fontsize',12)
    end

    title( strjoin( string( cell2mat( peakfeatures_PAMIL1first(conditions(counter-1), 12:14) ) ) ) );
    xticks([1/2 2/2 3/2 4/2])
    set(gca,'XTickLabel',{'0 ng/mL','0.03 ng/mL','0.1 ng/mL','0.3 ng/mL'},'fontsize',12)
    set(gca,'fontsize',12)
end

%% IL1 whole cell heatmap for F1
%TITLE IN N1 N2 N3 format. N1 = ligand 1, N2 = ligand 2, N3 = time
clear
load("./single_cell_features/peakfeatures_firstPAMIL1")
rng(5)
conditions = peakfeatures_PAMIL1first( ([peakfeatures_PAMIL1first{:, 11}] > 2 )', 15);
conditions = cell2mat(conditions);
conditions = conditions(~ismember(conditions, [30, 32, 62, 64] ));
position = [9, 1, 25, 17, 11, 3, 27, 19, 13, 5, 29, 21, 15, 7, 31, 23];
position = [position, position+32];
xlim = [0, 234];
figure(2) %plotting PIC effects
clf
conditions = conditions(conditions < 33);
ylim1 = [0 1.5];
clf
normfactors = [];
for aa=1:length(conditions)  
    temp1 = [peakfeatures_PAMIL1first{conditions(aa), 1}; peakfeatures_PAMIL1first{conditions(aa)+32, 1}] ;
    maxIL1s = cell2mat( peakfeatures_PAMIL1first( cell2mat( peakfeatures_PAMIL1first(:, 11) ) == 5, 2) );
    normfactors(1) = mean(maxIL1s);
    if peakfeatures_PAMIL1first{conditions(aa), 12} == 6 && peakfeatures_PAMIL1first{conditions(aa), 13} == 4
        normfactors(2) = mean([ peakfeatures_PAMIL1first{conditions(2), 7}; peakfeatures_PAMIL1first{conditions(2)+32, 7} ]);
    elseif peakfeatures_PAMIL1first{conditions(aa), 12} == 6 && peakfeatures_PAMIL1first{conditions(aa), 13} == 8
        normfactors(2) = mean([ peakfeatures_PAMIL1first{conditions(1), 7}; peakfeatures_PAMIL1first{conditions(1)+32, 7} ]);
    elseif peakfeatures_PAMIL1first{conditions(aa), 12} == 7 && peakfeatures_PAMIL1first{conditions(aa), 13} == 4
        normfactors(2) = mean([ peakfeatures_PAMIL1first{conditions(4), 7}; peakfeatures_PAMIL1first{conditions(4)+32, 7} ]);
    else
        normfactors(2) = mean([ peakfeatures_PAMIL1first{conditions(3), 7}; peakfeatures_PAMIL1first{conditions(3)+32, 7} ]);
    end
    [nr, ~] = size(temp1);
    subplot(4, 8, position(aa) )
    ncells = min(nr, 100);
    rows = datasample(1:nr, ncells, "replace", false);
    traces = temp1(rows, 1:40)./normfactors(1);
    maxs = max(traces, [], 2);
    [maxs,sortIdx]  = sort(maxs,'descend');
    imagesc(traces(sortIdx, :), ylim1)
    foo = string( cell2mat( peakfeatures_PAMIL1first(conditions(aa), 11:14) ) );
    title(strjoin(foo ) )  
    axis off
    temp2 = [peakfeatures_PAMIL1first{conditions(aa), 6}; peakfeatures_PAMIL1first{conditions(aa)+32, 6}] ;
    [nr, ~] = size(temp2);
    ncells = min(nr, 100);
    rows = datasample(1:nr, ncells, "replace", false);
    subplot(4, 8, position(aa)+1)
    traces = temp2(rows, 1:40)./normfactors(2);
    maxs = max(traces, [], 2);
    [maxs,sortIdx]  = sort(maxs,'descend');
    imagesc(traces(sortIdx, :), ylim1)
    foo = string( cell2mat( peakfeatures_PAMIL1first(conditions(aa), 11:13) ) );
    title(strjoin(foo ) )  
    axis off
end
hold off
%% IL-1 features for F2
%Subplot1 4 hour IL-1 to TNF, Subplot2 8 hours IL-1 to TNF, Subplot3 4 hours IL-1 to LPS, subplot4 8 hours IL-1 to LPS 
clc
clear
load("./single_cell_features/peakfeatures_firstPAMIL1")

%Feature 2-5 = Stim A max, width, AUC, and late AUC
%Feature 7-10 = Stim B max, AUC, early AUC, and late AUC
feature = 8;

conditions = peakfeatures_PAMIL1first( ([peakfeatures_PAMIL1first{:, 11}] > 2 )', 15);
conditions = cell2mat(conditions);
conditions = conditions(~ismember(conditions, [29:32, 61:64] ));
positions = [5, 1, 13, 9, 6, 2, 14, 10, 7, 3, 15, 11, 8, 4, 16, 12];
positions = [positions, positions+16];
conditions(positions) = conditions;
maxnorms = [mean(peakfeatures_PAMIL1first{1, feature}), mean(peakfeatures_PAMIL1first{2, feature}), ...
    mean(peakfeatures_PAMIL1first{3, feature}), mean(peakfeatures_PAMIL1first{4, feature}) ,...
    mean(peakfeatures_PAMIL1first{33, feature}), mean(peakfeatures_PAMIL1first{34, feature}), ...
    mean(peakfeatures_PAMIL1first{35, feature}), mean(peakfeatures_PAMIL1first{36, feature})];
figure(3)
clf
counter = 1;
for aa = 1:4
    subplot(4,1, aa)
    y = [];
    hold on
    for bb = 1:4
        cond1 = cell2mat( peakfeatures_PAMIL1first(conditions(counter), 12:14) );
        cond2 = cell2mat( peakfeatures_PAMIL1first(conditions(counter)+32, 12:14) );
        if isequal( cond1,  [6, 8, 1] )
            maxnorm1 = maxnorms(1);
            maxnorm2 = maxnorms(5);
        elseif isequal( cond1,  [6, 4, 1] )
            maxnorm1 = maxnorms(2);
            maxnorm2 = maxnorms(6);
        elseif isequal( cond1,  [7, 8, 1] )
            maxnorm1 = maxnorms(3);
            maxnorm2 = maxnorms(7);
        else
            maxnorm1 = maxnorms(4);
            maxnorm2 = maxnorms(8);
        end
        temp1 = peakfeatures_PAMIL1first{conditions(counter), feature}./maxnorm1 ;
        temp2 = peakfeatures_PAMIL1first{conditions(counter)+32, feature}./maxnorm2 ;
        Violin( [temp1; temp2], bb/2, 'Bandwidth', 0.3,  'MedianColor', [1,1,1], 'ViolinColor', [0 0 0], 'ViolinAlpha', .3, 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32/2, ...
            'BoxWidth', 0.1, 'ShowMean', true);
        y = [y, mean([temp1;temp2])];
        counter = counter + 1;
        controldist = [peakfeatures_PAMIL1first{conditions((aa-1)*4+1), feature}./maxnorm1; peakfeatures_PAMIL1first{conditions((aa-1)*4+1)+32, feature}./maxnorm2];
        p = ranksum(controldist, [temp1; temp2])*3
    end
    x=0:0.5:1;
    P = polyfit(x,y(2:4),1);
    y = polyval(P, -.5/2:.1:2.5/2);
    plot(1.5/2:.1:4.5/2,y,'r-', 'LineWidth',2);
    set(gca, 'XLim', [0.25, 4.75/2], 'YLim', [-0.25, 3.5])
    text(1/2, 3.25, string(round(P(1), 2)+ "*log(x) + " + round(P(2), 2) ),'fontsize',12)
    set(gca, 'XLim', [0.25/2, 4.75/2], 'YLim', [-0.5, 3.5])
    title( strjoin( string( cell2mat( peakfeatures_PAMIL1first(conditions(counter-1), 12:14) ) ) ) );
    xticks([1/2 2/2 3/2 4/2])
    set(gca,'XTickLabel',{'0 ng/mL','3 ng/mL','10 ng/mL','30 ng/mL'},'fontsize',12)
    set(gca,'fontsize',12)
end


%% sepsis single cell traces for F4
%TITLE IN N1 N2 format. N1 = ligand, N2 = control (0) or sepsis (1), N3 = time
clc
clear
load("./single_cell_features/peakfeatures_sepsis.mat") 
rng(1)
conditions = cell2mat( peakfeatures_sepsis( :, 6:8));
xlim = [0, 234];
ylim1 = [0 3];
figure(1) %plotting PIC effects
clf
rep = unique(conditions(:, 3));
treat = unique(conditions(:, 2));
lig = unique(conditions(:, 1));
counter = 1;
position = [1, 3, 5, 2, 4, 6];
for bb = 1:length(treat)
    for cc = 1:length(lig)
        [~, temp1]  = ismember(conditions(:,1:2), [lig(cc), treat(bb)], 'rows');
        foo = 1:length(temp1);
        temp1 = foo(temp1==1);
        temp1 = cell2mat( peakfeatures_sepsis(temp1, 1) ) ;
        [nr, ~] = size(temp1);
        subplot(1, 6, position(counter) )
        rows = datasample(1:nr, min(nr, 110), "Replace", false);
        maxs = max(temp1(rows, 1:40), [], 2);
        [maxs,sortIdx]  = sort(maxs,'descend');
        imagesc(temp1(rows(sortIdx(5:end-5)), 1:40), ylim1)
        foo2 = string( [lig(cc), treat(bb)] ); %string( [lig(cc), treat(bb), rep(aa)] );
        title(strjoin(foo2 ) )
        counter = counter + 1;
        axis off
    end
end

%% sepsis features
clc

%Feature 2-5 = Stim max, early AUC, AUC, and time to activation
feature = 4;

conditions = cell2mat( peakfeatures_sepsis( :, 6:8));
treat = unique(conditions(:, 2));
lig = unique(conditions(:, 1));

figure(2)
clf
counter = 1;
pos = [1, 3, 5, 2, 4, 6];
tempall = [];
for aa = 1:length(treat)
    for bb = 1:length(lig)
        [~, temp1]  = ismember(conditions(:,1:2), [lig(bb), treat(aa)], 'rows');
        foo = 1:length(temp1);
        temp1 = foo(temp1==1);
        temp1 = cell2mat( peakfeatures_sepsis(temp1, feature) ) ;
        bot5 = prctile(temp1, 5);
        top5 = prctile(temp1, 95);
        temp1 = temp1( temp1 <= top5 & temp1 >= bot5);
        Violin( temp1, pos(counter), 'Bandwidth', 3,  'MedianColor', [1,1,1], 'ViolinColor', [0 0 0], 'ViolinAlpha', .3, 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32, ...
            'BoxWidth', 0.1, 'ShowMean', true);
        set(gca, 'XLim', [0.25, 6.75], 'YLim', [-0.5, 80])
        set(gca,'XTickLabel',{'TNF-sepsis','TNF-control','LPS-sepsis','LPS-control','Media-sepsis','Media-control'},'fontsize',12)
        tempall{counter} = temp1;
        counter = counter + 1;
    end
end

ranksum(tempall{1}, tempall{4})*3
mean(tempall{1} )/mean(tempall{4})

ranksum(tempall{2}, tempall{5})*3
mean(tempall{2} )/mean(tempall{5})

ranksum(tempall{3}, tempall{6})*3
mean(tempall{3} )/mean(tempall{6})

