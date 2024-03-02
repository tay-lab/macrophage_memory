%% load files CPG PIC
clc
clear
load("./processedtraces/PICCPGfirst_traces.mat")
scmat = scmat_PICCPG;
clear('scmat_PICCPG')
%% finding features CPG PIC
peakfeatures = cell(11);
uniques = unique(scmat(:, 131));

for aa = 1:length(uniques)
    temp = scmat(scmat(:, 131) == uniques(aa), : );

    %%%Finds all cells tracked for entire interval of first 4 hour stimulus%%%
    peakfeatures{aa, 1} = temp(:, 3:42);
    peakfeatures{aa, 1} = peakfeatures{aa, 1}(~any(isnan( peakfeatures{aa, 1} ), 2 ), :);
    peakfeatures{aa, 1} = peakfeatures{aa, 1}(:, : ) - min( peakfeatures{aa, 1}(:, 1:2)' )';

    [nr1, ~] = size(peakfeatures{aa, 1} );
    pkmax = []; %Max amplitude
    pkwidth = []; %Peak width
    auc = []; %AUC over 240 minutes
    auc2 = []; %AUC over last 120 minutes

    %%%Finds peak features for each cell tracked in the entire first interval%%%
    for cc = 1:nr1
        sm_tr = smoothdata(peakfeatures{aa, 1}(cc, :), 'loess', 3);
        [~, idx, w, p] = findpeaks(sm_tr', 'MinPeakDistance',5 ,'MinPeakProminence', 0.1 );
        if isempty(idx)
            pkmax = [pkmax; 0];
            pkwidth = [pkwidth; 0 ];
        else
            pks = sortrows( [idx peakfeatures{aa, 1}(cc, idx )' ], 2, 'desc');
            pkmax = [pkmax; max(pks(:, 2))];
            pkwidth = [pkwidth; max(w) ];
        end
        auc = [auc; trapz(peakfeatures{aa, 1}(cc, 1:20)) ];
        auc2 = [auc2; trapz(peakfeatures{aa, 1}(cc, 11:20)) ];
    end
    peakfeatures{aa, 2} = pkmax;
    peakfeatures{aa, 3} = pkwidth;
    peakfeatures{aa, 4} = auc;
    peakfeatures{aa, 5} = auc2;

    %%%Finds all cells tracked for entire interval of second 4 hour stimulus%%%
    time = temp(1, 129);
    if time == 4
        peakfeatures{aa, 6} = temp(:, 43:82);
    else
        peakfeatures{aa, 6} = temp(:, 83:122);
    end
    peakfeatures{aa, 6} = peakfeatures{aa, 6}(~any(isnan( peakfeatures{aa, 6} ), 2 ), :);
    peakfeatures{aa, 6} = peakfeatures{aa, 6}(:, : ) -  min( peakfeatures{aa, 6}(:, 1:2)')';

    [nr1, ~] = size(peakfeatures{aa, 6} );
    pkmax = []; %Max amplitude
    auc = []; %AUC over 240 minutes
    auc1 = []; %AUC over first 120 minutes
    auc2 = []; %AUC over last 120 minutes
    
    %%%Finds peak features for each cell tracked in the entire second interval%%%
    for cc = 1:nr1
        sm_tr = smoothdata(peakfeatures{aa, 6}(cc, :), 'loess', 3);
        [~, idx, w, p] = findpeaks(sm_tr', 'MinPeakDistance',5 ,'MinPeakProminence', 0.1 );
        if isempty(idx)
            pkmax = [pkmax; 0];
        else
            pks = sortrows( [idx peakfeatures{aa, 6}(cc, idx )' ], 2, 'desc');
            pkmax = [pkmax; max(pks(:, 2))];
        end
        auc = [auc; trapz(peakfeatures{aa, 6}(cc, 1:40)) ];
        auc1 = [auc1; trapz(peakfeatures{aa, 6}(cc, 1:20)) ];
        auc2 = [auc2; trapz(peakfeatures{aa, 6}(cc, 11:40)) ];
    end
    peakfeatures{aa, 7} = max(pkmax, 0);
    peakfeatures{aa, 8} = max(auc,0);
    peakfeatures{aa, 9} = max(auc1,0); 
    peakfeatures{aa, 10} = max(auc2,0);


    %%%Carries over metadata%%%
    peakfeatures{aa, 11} = temp(1, 127);
    peakfeatures{aa, 12} = temp(1, 128);
    peakfeatures{aa, 13} = temp(1, 129);
    peakfeatures{aa, 14} = temp(1, 130);
    peakfeatures{aa, 15} = temp(1, 131);
end

peakfeatures_PICCPGfirst = peakfeatures;
save('./single_cell_features/peakfeatures_firstPICCPG', 'peakfeatures_PICCPGfirst' );

%% load files TNF LPS
clc
clear
load("./processedtraces/TNFLPSfirst_traces.mat")
scmat = scmat_TNFLPS;
clear('scmat_TNFLPS')

%% finding features TNF LPS
peakfeatures = cell(11);
counter = 1;
uniques = unique(scmat(:, 131));

for aa = 1:length(uniques)
    temp = scmat(scmat(:, 131) == uniques(aa), : );
    
    peakfeatures{aa, 1} = temp(:, 4:43);
    peakfeatures{aa, 1} = peakfeatures{aa, 1}(~any(isnan( peakfeatures{aa, 1} ), 2 ), :);
    peakfeatures{aa, 1} = peakfeatures{aa, 1}(:, : ) -  min( peakfeatures{aa, 1}(:, 1:2)' )';

    [nr1, ~] = size(peakfeatures{aa, 1} );
    pkmax = []; %Max amplitude
    pkwidth = []; %Peak width
    auc = []; %AUC over 240 minutes
    auc2 = []; %AUC over last 120 minutes

    for cc = 1:nr1
        sm_tr = smoothdata(peakfeatures{aa, 1}(cc, :), 'loess', 3);
        [~, idx, w, p] = findpeaks(sm_tr', 'MinPeakDistance',5 ,'MinPeakProminence', 0.1 );
        if isempty(idx)
            pkmax = [pkmax; 0];
            pkwidth = [pkwidth; 0 ];
        else
            pks = sortrows( [idx peakfeatures{aa, 1}(cc, idx )' ], 2, 'desc');
            pkmax = [pkmax; max(pks(:, 2))];
            pkwidth = [pkwidth; max(w) ];
        end
        auc = [auc; trapz(peakfeatures{aa, 1}(cc, :)) ];
        auc2 = [auc2; trapz(peakfeatures{aa, 1}(cc, 11:20)) ];
    end
    peakfeatures{aa, 2} = pkmax;
    peakfeatures{aa, 3} = pkwidth;
    peakfeatures{aa, 4} = auc;
    peakfeatures{aa, 5} = auc2;

    
    time = temp(1, 129);
    if time == 4
        peakfeatures{aa, 6} = temp(:, 44:83);
    else
        peakfeatures{aa, 6} = temp(:, 84:123);
    end
    peakfeatures{aa, 6} = peakfeatures{aa, 6}(~any(isnan( peakfeatures{aa, 6} ), 2 ), :);
    peakfeatures{aa, 6} = peakfeatures{aa, 6}(:, : ) -  min( peakfeatures{aa, 6}(:, 1:2)' )';

    [nr1, ~] = size(peakfeatures{aa, 6} );
    pkmax = [];
    auc = [];
    auc1 = [];
    auc2 = [];

    for cc = 1:nr1
        sm_tr = smoothdata(peakfeatures{aa, 6}(cc, :), 'loess', 3);
        [~, idx, w, p] = findpeaks(sm_tr', 'MinPeakDistance',5 ,'MinPeakProminence', 0.1 );
        if isempty(idx)
            pkmax = [pkmax; 0];
        else
            pks = sortrows( [idx peakfeatures{aa, 6}(cc, idx )' ], 2, 'desc');
            pkmax = [pkmax; max(pks(:, 2))];
        end
        auc = [auc; trapz(peakfeatures{aa, 6}(cc, 1:40)) ];
        auc1 = [auc1; trapz(peakfeatures{aa, 6}(cc, 1:15)) ];
        auc2 = [auc2; trapz(peakfeatures{aa, 6}(cc, 16:40)) ];
    end
    peakfeatures{aa, 7} = pkmax;
    peakfeatures{aa, 8} = max(auc,0);
    peakfeatures{aa, 9} = max(auc1,0);
    peakfeatures{aa, 10} = max(auc2,0);

    peakfeatures{aa, 11} = temp(1, 127);
    peakfeatures{aa, 12} = temp(1, 128);
    peakfeatures{aa, 13} = temp(1, 129);
    peakfeatures{aa, 14} = temp(1, 130);
    peakfeatures{aa, 15} = temp(1, 131);

    %%%Dealing with repeated stimulus with same dose%%%
    if ismember(aa, [13, 29, 13+32, 29+32])
        newindex = 64+counter;
        counter = counter + 1;

        [peakfeatures{newindex, 1}, peakfeatures{newindex, 2}, peakfeatures{newindex, 3}, peakfeatures{newindex, 4}, peakfeatures{newindex, 5}] = peakfeatures{aa, 1:5};
        peakfeatures{newindex, 6} = temp(:, 44:83);
        peakfeatures{newindex, 6} = peakfeatures{newindex, 6}(~any(isnan( peakfeatures{newindex, 6} ), 2 ), :);
        peakfeatures{newindex, 6} = peakfeatures{newindex, 6}(:, : ) -  mean( peakfeatures{newindex, 6}(:, 1:2), 2 );

        [nr1, ~] = size(peakfeatures{newindex, 6} );
        pkp = [];
        auc = [];
        auc1 = [];
        auc2 = [];

        for cc = 1:nr1
            sm_tr = smoothdata(peakfeatures{newindex, 6}(cc, :), 'loess', 3);
            [~, idx, w, p] = findpeaks(sm_tr', 'MinPeakDistance',5 ,'MinPeakProminence', 0.1 );
            if isempty(idx)
                pkmax = [pkmax; 0];
            else
                pks = sortrows( [idx peakfeatures{newindex, 6}(cc, idx )' ], 2, 'desc');
                pkmax = [pkmax; max(pks(:, 2))];
            end
            auc = [auc; trapz(peakfeatures{newindex, 6}(cc, 1:40)) ];
            auc1 = [auc1; trapz(peakfeatures{newindex, 6}(cc, 1:15)) ];
            auc2 = [auc2; trapz(peakfeatures{newindex, 6}(cc, 16:40)) ];
        end
        peakfeatures{newindex, 7} = pkmax;
        peakfeatures{newindex, 8} = max(auc,0);
        peakfeatures{newindex, 9} = max(auc1,0);
        peakfeatures{newindex, 10} = max(auc2,0);


        peakfeatures{newindex, 11} = temp(1, 127);
        peakfeatures{newindex, 12} = temp(1, 128);
        peakfeatures{newindex, 13} = 4;
        peakfeatures{newindex, 14} = temp(1, 130);
        peakfeatures{newindex, 15} = newindex;
    end
end

peakfeatures_TNFLPSfirst = peakfeatures;
save('./single_cell_features/peakfeatures_firstTNFLPS', 'peakfeatures_TNFLPSfirst' );

%% load files PAM IL1
clc
clear
load("./processedtraces/PAMIL1first_traces.mat")
scmat = scmat_PAMIL1;
clear('scmat_PAMIL1')
%% finding features PAM IL1

peakfeatures = cell(11);

uniques = unique(scmat(:, 131));
for aa = 1:length(uniques)
    temp = scmat(scmat(:, 131) == uniques(aa), : );
    peakfeatures{aa, 1} = temp(:, 4:43);
    peakfeatures{aa, 1} = peakfeatures{aa, 1}(~any(isnan( peakfeatures{aa, 1} ), 2 ), :);
    peakfeatures{aa, 1} = peakfeatures{aa, 1}(:, : ) - min( peakfeatures{aa, 1}(:, 1:2)' )';

    [nr1, ~] = size(peakfeatures{aa, 1} );
    pkmax = []; %Max amplitude
    pkwidth = []; %Peak width
    auc = []; %AUC over 240 minutes
    auc2 = []; %AUC over last 120 minutes

    for cc = 1:nr1
        sm_tr = smoothdata(peakfeatures{aa, 1}(cc, :), 'loess', 3);
        [~, idx, w, p] = findpeaks(sm_tr', 'MinPeakDistance',5 ,'MinPeakProminence', 0.1 );
        if isempty(idx)
            pkmax = [pkmax; 0];
        else
            pks = sortrows( [idx peakfeatures{aa, 1}(cc, idx )' ], 2, 'desc');
            pkmax = [pkmax; max(pks(:, 2))];
        end
        auc = [auc; trapz(peakfeatures{aa, 1}(cc, :)) ];
        auc2 = [auc2; trapz(peakfeatures{aa, 1}(cc, 11:20)) ];
    end
    peakfeatures{aa, 2} = pkmax;
    peakfeatures{aa, 3} = pkwidth;
    peakfeatures{aa, 4} = auc;
    peakfeatures{aa, 5} = auc2;

    
    time = temp(1, 129);
    if time == 4
        peakfeatures{aa, 6} = temp(:, 44:83);
    else
        peakfeatures{aa, 6} = temp(:, 84:123);
    end
    peakfeatures{aa, 6} = peakfeatures{aa, 6}(~any(isnan( peakfeatures{aa, 6} ), 2 ), :);
    peakfeatures{aa, 6} = peakfeatures{aa, 6}(:, : ) -  min( peakfeatures{aa, 6}(:, 1:2)' )';

    [nr1, ~] = size(peakfeatures{aa, 6} );
    pkmax = [];
    auc = [];
    auc1 = [];
    auc2 = [];

    for cc = 1:nr1
        sm_tr = smoothdata(peakfeatures{aa, 6}(cc, :), 'loess', 3);
        [~, idx, w, p] = findpeaks(sm_tr', 'MinPeakDistance',5 ,'MinPeakProminence', 0.1 );
        if isempty(idx)
            pkmax = [pkmax; 0];
        else
            pks = sortrows( [idx peakfeatures{aa, 6}(cc, idx )' ], 2, 'desc');
            pkmax = [pkmax; max(pks(:, 2))];
        end
        auc = [auc; trapz(peakfeatures{aa, 6}(cc, 1:40)) ];
        auc1 = [auc1; trapz(peakfeatures{aa, 6}(cc, 1:15)) ];
        auc2 = [auc2; trapz(peakfeatures{aa, 6}(cc, 16:40)) ];
    end
    peakfeatures{aa, 7} = pkmax;
    peakfeatures{aa, 8} = max(auc,0);
    peakfeatures{aa, 9} = max(auc1,0);
    peakfeatures{aa, 10} = max(auc2,0);

    peakfeatures{aa, 11} = temp(1, 127);
    peakfeatures{aa, 12} = temp(1, 128);
    peakfeatures{aa, 13} = temp(1, 129);
    peakfeatures{aa, 14} = temp(1, 130);
    peakfeatures{aa, 15} = temp(1, 131);
end

peakfeatures_PAMIL1first = peakfeatures;
save('./single_cell_features/peakfeatures_firstPAMIL1', 'peakfeatures_PAMIL1first' );

%% load files CPG to LPS single cell replicates
clc
clear
load("./processedtraces/CPGLPSscrep_traces.mat")
scmat = scmat_screps3;
clear('scmat_screps3')

%% finding full cells CPG to LPS single cell replicates
wdw = 3; %smoothing window

peakfeatures = cell(11);

uniques = unique(scmat(:, 131));
counter = 1;
for aa = 1:length(uniques)
    temp = scmat(scmat(:, 131) == uniques(aa), : );
    time = temp(1, 129);
    if time == 4

        %%%Finds all cells tracked for entire interval of all 8 hours%%%
        peakfeatures{counter, 1} = temp(:, 3:82);
        peakfeatures{counter, 1} = peakfeatures{counter, 1}(~any(isnan( peakfeatures{counter, 1} ), 2 ), :);
        peakfeatures{counter, 1} = peakfeatures{counter, 1}(:, : )- min( peakfeatures{counter, 1}(:, 1:2)' )';  

        [nr1, ~] = size(peakfeatures{counter, 1} );
        pkmax = []; %max amplitude
        pkwidth = []; %time between first and second half max
        pktime = []; %time to half max
        auc = []; %AUC of first interval

        %%%Finds features for first 4 hour interval from completely tracked cells%%%
        for cc = 1:nr1
            sm_tr = smoothdata(peakfeatures{counter, 1}(cc, 1:40), 'loess', wdw);
            [smax, idx, w, p] = findpeaks(sm_tr', 'MinPeakDistance',5 ,'MinPeakProminence', 0.1 );
            if isempty(idx)
                pkmax = [pkmax; 0]; 
                pkwidth = [pkwidth; 0 ]; 
                pktime = [pktime; 0]; 
            else
                pks = sortrows( [idx peakfeatures{counter, 1}(cc, idx )' ], 2, 'desc');
                [pkmaxtemp, pktimetemp] = max(pks(:, 2));
                if pkmaxtemp < 0
                    pkmax = [pkmax; 0];
                    pkwidth = [pkwidth; 40 ];
                    pktime = [pktime; 40];
                else

                    hpktimetemp = find( peakfeatures{counter, 1}(cc, 1:39) < pkmaxtemp/2 & peakfeatures{counter, 1}(cc,2:40) > pkmaxtemp/2);
                    hpkwidthtemp = find( peakfeatures{counter, 1}(cc, 1:39) > pkmaxtemp/2 & peakfeatures{counter, 1}(cc,2:40) < pkmaxtemp/2);
                    pkmax = [pkmax; pkmaxtemp];
                    pktime = [pktime; hpktimetemp(1)];
                    len = max(hpkwidthtemp);
                    if len < pks(1)
                        pkwidth = [pkwidth; 40];
                    else
                        widthidx = find(hpkwidthtemp > pks(1));
                        if isempty( widthidx )
                            pkwidth = [pkwidth; 40-hpktimetemp(1)];
                        else
                            pkwidth = [pkwidth; hpkwidthtemp(min(widthidx))-hpktimetemp(1)];
                        end
                        
                    end

                end
            end
            auc = [auc; trapz(peakfeatures{counter, 1}(cc, 1:40)) ];
        end
        peakfeatures{counter, 2} = max( pkmax, 0);
        peakfeatures{counter, 3} = max( auc, 0);
        peakfeatures{counter, 4} = pktime;
        peakfeatures{counter, 5} = pkwidth;

        pkmax = [];
        pkwidth = [];
        pktime = [];
        pkp = [];
        auc = [];
        auc1 = [];
        auc2 = [];

        %%%Finds features for second 4 hour interval from completely tracked cells%%%
        temptr = peakfeatures{counter, 1}(:, 41:80) - min( peakfeatures{counter, 1}(:, 41:42)' )';

        for cc = 1:nr1
            sm_tr = smoothdata(temptr(cc, :), 'loess', wdw);
            [~, idx, w, p] = findpeaks(sm_tr', 'MinPeakDistance',5 ,'MinPeakProminence', 0.1 );
            if isempty(idx)
                pkmax = [pkmax; 0];
                pkwidth = [pkwidth; 0 ];
                pktime = [pktime; 0];
                pkp = [pkp; 0];
            else
                pks = sortrows( [idx temptr(cc, idx )' ], 2, 'desc');
                pkmax = [pkmax; max(pks(:, 2))];
                pkwidth = [pkwidth; max(w) ];
                pktime = [pktime; idx(1)];
                pkp = [pkp; p];
            end
            auc = [auc; trapz(temptr(cc, 1:40))];
            auc1 = [auc1; trapz(temptr(cc, 1:15)) ];
            auc2 = [auc2; trapz(temptr(cc, 16:40)) ];
        end
        peakfeatures{counter, 6} = max( pkmax, 0);
        peakfeatures{counter, 7} = max( auc, 0);
        peakfeatures{counter, 8} = max( auc1, 0);
        peakfeatures{counter, 9} = max( auc2, 0);

        peakfeatures{counter, 10} = temp(1, 127);
        peakfeatures{counter, 11} = temp(1, 128);
        peakfeatures{counter, 12} = temp(1, 129);
        peakfeatures{counter, 13} = temp(1, 130);
        peakfeatures{counter, 14} = temp(1, 131);
        counter = counter + 1;
    end


end

fulltraces_CpGLPSscrep = peakfeatures;
save('./single_cell_features/fulltraces_CpGLPSscrep', 'fulltraces_CpGLPSscrep' );

%% load files polyI:C to TNF single cell replicates
clc
clear
load("./processedtraces/PICTNFscrep_traces.mat")
scmat = scmat_screps5;
clear('scmat_screps5')

%% finding full cells PIC screps5
peakfeatures = cell(11);

uniques = unique(scmat(:, 131));
counter = 1;
for aa = 1:length(uniques)
    temp = scmat(scmat(:, 131) == uniques(aa), : );
    time = temp(1, 129);
    if time == 4

        %%%Finds all cells tracked for entire interval of all 8 hours%%%
        peakfeatures{counter, 1} = temp(:, 3:82);
        peakfeatures{counter, 1} = peakfeatures{counter, 1}(~any(isnan( peakfeatures{counter, 1} ), 2 ), :);
        peakfeatures{counter, 1} = peakfeatures{counter, 1}(:, : )- min( peakfeatures{counter, 1}(:, 1:2)' )';  % -peakfeatures{counter, 1}(:, 1);

        [nr1, ~] = size(peakfeatures{counter, 1} );
        pkmax = []; %max amplitude
        pkwidth = []; %time between first and second half max
        pktime = []; %time to half max
        auc = []; %AUC of first interval
        
        %%%Finds features for first 4 hour interval from completely tracked cells%%%
        for cc = 1:nr1
            sm_tr = smoothdata(peakfeatures{counter, 1}(cc, 1:40), 'loess', 3);
            [~, idx, w, p] = findpeaks(sm_tr', 'MinPeakDistance',5 ,'MinPeakProminence', 0.1 );
            if isempty(idx)
                pkmax = [pkmax; 0];
                pkwidth = [pkwidth; 0 ];
                pktime = [pktime; 40];
            else
                pks = sortrows( [idx peakfeatures{counter, 1}(cc, idx )' ], 2, 'desc');
                [pkmaxtemp, pktimetemp] = max(pks(:, 2));
                if pkmaxtemp < 0 
                    pkmax = [pkmax; 0];
                    pkwidth = [pkwidth; 0 ];
                    pktime = [pktime; 40];
                else

                    hpktimetemp = find( peakfeatures{counter, 1}(cc, 1:end-1) < pkmaxtemp/2 & peakfeatures{counter, 1}(cc,2:end) > pkmaxtemp/2);
                    hpkwidthtemp = find( peakfeatures{counter, 1}(cc, 1:39) > pkmaxtemp/2 & peakfeatures{counter, 1}(cc,2:40) < pkmaxtemp/2);
                    pkmax = [pkmax; pkmaxtemp];
                    pktime = [pktime; hpktimetemp(1)];
                    len = max(hpkwidthtemp);
                    if len < pks(1)
                        pkwidth = [pkwidth; 40];
                    else
                        widthidx = find(hpkwidthtemp > pks(1));
                        if isempty( widthidx )
                            pkwidth = [pkwidth; 40-hpktimetemp(1)];
                        else
                            pkwidth = [pkwidth; hpkwidthtemp(min(widthidx))-hpktimetemp(1)];
                        end
                    end
                end

            end
            auc = [auc; trapz(peakfeatures{counter, 1}(cc, 1:40)) ];

        end
        peakfeatures{counter, 2} = max( pkmax, 0);
        peakfeatures{counter, 3} = max( auc, 0);
        peakfeatures{counter, 4} = pktime; 
        peakfeatures{counter, 5} = pkwidth;

        pkmax = [];
        pktime = [];
        auc = [];
        auc2 = [];

        %%%Finds features for second 4 hour interval from completely tracked cells%%%
        temptr = peakfeatures{counter, 1}(:, 41:80) - min( peakfeatures{counter, 1}(:, 41:42)' )'; 
        for cc = 1:nr1
            sm_tr = smoothdata(temptr(cc, :), 'loess', 3);
            [~, idx, w, p] = findpeaks(sm_tr', 'MinPeakDistance',5 ,'MinPeakProminence', 0.1 );
            if isempty(idx)
                pkmax = [pkmax; 0];
                pkwidth = [pkwidth; 0 ];
                pktime = [pktime; 0];
            else
                pks = sortrows( [idx temptr(cc, idx )' ], 2, 'desc');
                [pkmaxtemp, pktimetemp] = max(pks(:, 2));
                pkmax = [pkmax; pkmaxtemp];
                pkwidth = [pkwidth; w(pktimetemp) ];
                pktime = [pktime; idx(pktimetemp)];

            end
            auc = [auc; trapz(temptr(cc, 1:40))];
            auc2 = [auc2; trapz(temptr(cc, 21:40)) ];
        end
        peakfeatures{counter, 6} = max( pkmax, 0);
        peakfeatures{counter, 7} = max( auc, 0);
        peakfeatures{counter, 8} = pktime; 
        peakfeatures{counter, 9} = max( auc2, 0);     



        peakfeatures{counter, 10} = temp(1, 127);
        peakfeatures{counter, 11} = temp(1, 128);
        peakfeatures{counter, 12} = temp(1, 129);
        peakfeatures{counter, 13} = temp(1, 130);
        peakfeatures{counter, 14} = temp(1, 131);
        counter = counter + 1;
    end


end

fulltraces_PICTNFscrep = peakfeatures;
save('./single_cell_features/fulltraces_PICTNFscrep', 'fulltraces_PICTNFscrep' );

%% load files TNFblock PIC
clc
clear
load("./processedtraces/PICfirst_TNFblock_traces.mat")
scmat = scmat_TNFblock;
clear('scmat_TNFblock')
%% finding features TNFblock

peakfeatures = cell(11);

uniques = unique(scmat(:, 131));
for aa = 1:length(uniques)
    temp = scmat(scmat(:, 131) == uniques(aa), : );

    %%%Finds trace and features for first 4 hour interval%%%
    peakfeatures{aa, 1} = temp(:, 3:42);
    peakfeatures{aa, 1} = peakfeatures{aa, 1}(~any(isnan( peakfeatures{aa, 1} ), 2 ), :);
    peakfeatures{aa, 1} = peakfeatures{aa, 1}(:, : ) - min( peakfeatures{aa, 1}(:, 1:2)' )';

    [nr1, ~] = size(peakfeatures{aa, 1} );
    pkmax = [];
    pktime = [];
    auc = [];
    auc2 = [];

    for cc = 1:nr1
        sm_tr = smoothdata(peakfeatures{aa, 1}(cc, :), 'loess', 3);
        [~, idx, w, p] = findpeaks(sm_tr', 'MinPeakDistance',5 ,'MinPeakProminence', 0.1 );
        if isempty(idx)
            pkmax = [pkmax; 0];
            pktime = [pktime; 0];
        else
            pks = sortrows( [idx peakfeatures{aa, 1}(cc, idx )' ], 2, 'desc');
            [pkmaxtemp, pktimetemp] = max(pks(:, 2));
            if pkmaxtemp < 0.75
                pktime = [pktime; 40];
            else
                hpktimetemp = find( peakfeatures{aa, 1}(cc, 1:end-1) < pkmaxtemp/2 & peakfeatures{aa, 1}(cc,2:end) > pkmaxtemp/2);
                pktime = [pktime; hpktimetemp(1)];
            end
            pkmax = [pkmax; pkmaxtemp];
        end
        auc = [auc; trapz(peakfeatures{aa, 1}(cc, 1:40)) ];
        auc2 = [auc2; trapz(peakfeatures{aa, 1}(cc, 21:40)) ];
    end
    peakfeatures{aa, 2} = pkmax;
    peakfeatures{aa, 3} = pktime;
    peakfeatures{aa, 4} = auc;
    peakfeatures{aa, 5} = auc2;

    %%%Finds trace and features for second 4 hour interval%%%
    time = temp(1, 129);
    if time == 4
        peakfeatures{aa, 6} = temp(:, 43:82);
    else
        peakfeatures{aa, 6} = temp(:, 83:122);
    end
    peakfeatures{aa, 6} = peakfeatures{aa, 6}(~any(isnan( peakfeatures{aa, 6} ), 2 ), :);
    peakfeatures{aa, 6} = peakfeatures{aa, 6}(:, : ) -  min( peakfeatures{aa, 6}(:, 1:2)')';

    [nr1, ~] = size(peakfeatures{aa, 6} );
    pkmax = [];
    auc = [];
    auc1 = [];

    for cc = 1:nr1
        sm_tr = smoothdata(peakfeatures{aa, 6}(cc, :), 'loess', 3);
        [~, idx, w, p] = findpeaks(sm_tr', 'MinPeakDistance',5 ,'MinPeakProminence', 0.1 );
        if isempty(idx)
            pkmax = [pkmax; 0];
        else
            pks = sortrows( [idx peakfeatures{aa, 6}(cc, idx )' ], 2, 'desc');
            pkmax = [pkmax; max(pks(:, 2))];
        end
        auc = [auc; trapz(peakfeatures{aa, 6}(cc, 1:40)) ];
        auc1 = [auc1; trapz(peakfeatures{aa, 6}(cc, 1:20)) ];
    end
    peakfeatures{aa, 7} = max(pkmax, 0);
    peakfeatures{aa, 8} = max(auc,0);
    peakfeatures{aa, 9} = max(auc1,0);
    peakfeatures{aa, 10} = max(auc,0)./(max(pkmax,0));


    peakfeatures{aa, 11} = temp(1, 127);
    peakfeatures{aa, 12} = temp(1, 128);
    peakfeatures{aa, 13} = temp(1, 129);
    peakfeatures{aa, 14} = temp(1, 130);
    peakfeatures{aa, 15} = temp(1, 131);
end

peakfeatures_PICfirst_TNFblock = peakfeatures;
%save('./single_cell_features/peakfeatures_PICfirst_TNFblock', 'peakfeatures_PICfirst_TNFblock' );

%% load files TNFblock LPS
clc
clear
load("./processedtraces/CPGfirst_TNFblock_traces.mat")
scmat = scmat_TNFblock2;
clear('scmat_TNFblock2')
%% finding features TNFblock
peakfeatures = cell(11);

uniques = unique(scmat(:, 131));
for aa = 1:length(uniques)
    temp = scmat(scmat(:, 131) == uniques(aa), : );
    peakfeatures{aa, 1} = temp(:, 3:42);
    peakfeatures{aa, 1} = peakfeatures{aa, 1}(~any(isnan( peakfeatures{aa, 1} ), 2 ), :);
    peakfeatures{aa, 1} = peakfeatures{aa, 1}(:, : ) - min( peakfeatures{aa, 1}(:, 1:2)' )';

    [nr1, ~] = size(peakfeatures{aa, 1} );
    pkmax = [];
    pktime = [];
    auc = [];
    auc2 = [];

    for cc = 1:nr1
        sm_tr = smoothdata(peakfeatures{aa, 1}(cc, :), 'loess', 3);
        [~, idx, w, p] = findpeaks(sm_tr', 'MinPeakDistance',5 ,'MinPeakProminence', 0.1 );
        if isempty(idx)
            pkmax = [pkmax; 0];
            pktime = [pktime; 40];
        else
            pks = sortrows( [idx peakfeatures{aa, 1}(cc, idx )' ], 2, 'desc');
            pkmaxtemp = max(pks(:, 2));
            if pkmaxtemp < 0.75
                pktime = [pktime; 40];
            else
                hpktimetemp = find( peakfeatures{aa, 1}(cc, 1:end-1) < pkmaxtemp/2 & peakfeatures{aa, 1}(cc,2:end) > pkmaxtemp/2);
                pktime = [pktime; hpktimetemp(1)];
            end
            pkmax = [pkmax; pkmaxtemp];
        end
        auc = [auc; trapz(peakfeatures{aa, 1}(cc, 1:40)) ];
        auc2 = [auc2; trapz(peakfeatures{aa, 1}(cc, 21:40)) ];
    end
    peakfeatures{aa, 2} = pkmax;
    peakfeatures{aa, 3} = pktime;
    peakfeatures{aa, 4} = auc;
    peakfeatures{aa, 5} = auc2;

    
    time = temp(1, 129);
    if time == 4
        peakfeatures{aa, 6} = temp(:, 43:82);
    else
        peakfeatures{aa, 6} = temp(:, 83:122);
    end
    peakfeatures{aa, 6} = peakfeatures{aa, 6}(~any(isnan( peakfeatures{aa, 6} ), 2 ), :);
    peakfeatures{aa, 6} = peakfeatures{aa, 6}(:, : ) -  min( peakfeatures{aa, 6}(:, 1:2)')';

    [nr1, ~] = size(peakfeatures{aa, 6} );
    pkmax = [];
    auc = [];
    auc1 = [];

    for cc = 1:nr1
        sm_tr = smoothdata(peakfeatures{aa, 6}(cc, :), 'loess', 3);
        [~, idx, w, p] = findpeaks(sm_tr', 'MinPeakDistance',5 ,'MinPeakProminence', 0.1 );
        if isempty(idx)
            pkmax = [pkmax; 0];
        else
            pks = sortrows( [idx peakfeatures{aa, 6}(cc, idx )' ], 2, 'desc');
            pkmax = [pkmax; max(pks(:, 2))];
        end
        auc = [auc; trapz(peakfeatures{aa, 6}(cc, 1:40)) ];
        auc1 = [auc1; trapz(peakfeatures{aa, 6}(cc, 1:20)) ];
        auc2 = [auc2; trapz(peakfeatures{aa, 6}(cc, 21:40)) ];
    end
    peakfeatures{aa, 7} = max(pkmax, 0);
    peakfeatures{aa, 8} = max(auc,0);
    peakfeatures{aa, 9} = max(auc1,0);
    peakfeatures{aa, 10} = max(auc,0)./(max(pkmax,0));



    peakfeatures{aa, 11} = temp(1, 127);
    peakfeatures{aa, 12} = temp(1, 128);
    peakfeatures{aa, 13} = temp(1, 129);
    peakfeatures{aa, 14} = temp(1, 130);
    peakfeatures{aa, 15} = temp(1, 131);
end

peakfeatures_CPGfirst_TNFblock = peakfeatures;
save('./single_cell_features/peakfeatures_CPGfirst_TNFblock', 'peakfeatures_CPGfirst_TNFblock' );

%% load files IL10block
clc
clear
load("./processedtraces/CpGfirst_IL10block_traces.mat")
scmat = scmat_IL10block;
clear('scmat_IL10block')
%% finding features IL10block
peakfeatures = cell(11);

uniques = unique(scmat(:, 131));
for aa = 1:length(uniques)
    temp = scmat(scmat(:, 131) == uniques(aa), : );
    peakfeatures{aa, 1} = temp(:, 3:42);
    peakfeatures{aa, 1} = peakfeatures{aa, 1}(~any(isnan( peakfeatures{aa, 1} ), 2 ), :);
    peakfeatures{aa, 1} = peakfeatures{aa, 1}(:, : ) - min( peakfeatures{aa, 1}(:, 1:2)' )';

    [nr1, ~] = size(peakfeatures{aa, 1} );
    pkmax = [];
    pktime = [];
    auc = [];
    auc2 = [];

    for cc = 1:nr1
        sm_tr = smoothdata(peakfeatures{aa, 1}(cc, :), 'loess', 3);
        [~, idx, w, p] = findpeaks(sm_tr', 'MinPeakDistance',5 ,'MinPeakProminence', 0.1 );
        if isempty(idx)
            pkmax = [pkmax; 0];
            pktime = [pktime; 0];
        else
            pks = sortrows( [idx peakfeatures{aa, 1}(cc, idx )' ], 2, 'desc');
            [pkmaxtemp, pktimetemp] = max(pks(:, 2));
            if pkmaxtemp < 0.75
                pktime = [pktime; 40];
            else
                hpktimetemp = find( peakfeatures{aa, 1}(cc, 1:end-1) < pkmaxtemp/2 & peakfeatures{aa, 1}(cc,2:end) > pkmaxtemp/2);
                pktime = [pktime; hpktimetemp(1)];
            end
            pkmax = [pkmax; pkmaxtemp];
        end
        auc = [auc; trapz(peakfeatures{aa, 1}(cc, 1:40)) ];
        auc2 = [auc2; trapz(peakfeatures{aa, 1}(cc, 21:40)) ];
    end
    peakfeatures{aa, 2} = pkmax;
    peakfeatures{aa, 3} = pktime;
    peakfeatures{aa, 4} = auc;
    peakfeatures{aa, 5} = auc2;

    
    time = temp(1, 129);
    if time == 4
        peakfeatures{aa, 6} = temp(:, 43:82);
    else
        peakfeatures{aa, 6} = temp(:, 83:122);
    end
    peakfeatures{aa, 6} = peakfeatures{aa, 6}(~any(isnan( peakfeatures{aa, 6} ), 2 ), :);
    peakfeatures{aa, 6} = peakfeatures{aa, 6}(:, : ) -  min( peakfeatures{aa, 6}(:, 1:2)')';

    [nr1, ~] = size(peakfeatures{aa, 6} );
    pkmax = [];
    auc = [];
    auc1 = [];

    for cc = 1:nr1
        sm_tr = smoothdata(peakfeatures{aa, 6}(cc, :), 'loess', 3);
        [~, idx, w, p] = findpeaks(sm_tr', 'MinPeakDistance',5 ,'MinPeakProminence', 0.1 );
        if isempty(idx)
            pkmax = [pkmax; 0];
        else
            pks = sortrows( [idx peakfeatures{aa, 6}(cc, idx )' ], 2, 'desc');
            pkmax = [pkmax; max(pks(:, 2))];
        end
        auc = [auc; trapz(peakfeatures{aa, 6}(cc, 1:40)) ];
        auc1 = [auc1; trapz(peakfeatures{aa, 6}(cc, 1:20)) ];
    end
    peakfeatures{aa, 7} = max(pkmax, 0);
    peakfeatures{aa, 8} = max(auc,0);
    peakfeatures{aa, 9} = max(auc1,0); 
    peakfeatures{aa, 10} = max(auc,0)./(max(pkmax,0));

    peakfeatures{aa, 11} = temp(1, 127);
    peakfeatures{aa, 12} = temp(1, 128);
    peakfeatures{aa, 13} = temp(1, 129);
    peakfeatures{aa, 14} = temp(1, 130);
    peakfeatures{aa, 15} = temp(1, 131);
end

peakfeatures_CPGfirst_IL10block = peakfeatures;
save('./single_cell_features/peakfeatures_CPGfirst_IL10block', 'peakfeatures_CPGfirst_IL10block' );

%% load files sepsis
clc
clear
load("./processedtraces/sepsis_traces.mat")
scmat = scmat_PMP;
clear('scmat_PMP')
%% finding features PMP

peakfeatures = cell(9);

uniques = unique(scmat(:, 109));
for aa = 1:length(uniques)
    temp = scmat(scmat(:, 109) == uniques(aa), : );

    %%%Finds  cells tracked for 4 hours%%%
    peakfeatures{aa, 1} = temp(:, 3:42);
    peakfeatures{aa, 1} = peakfeatures{aa, 1}(~any(isnan( peakfeatures{aa, 1} ), 2 ), :);
    peakfeatures{aa, 1} = peakfeatures{aa, 1}(:, : ) - min( peakfeatures{aa, 1}(:, 1:2)' )';
    peakfeatures{aa, 1} = smoothdata( peakfeatures{aa, 1}(:, : ) , 'loess', 5);

    [nr1, ~] = size(peakfeatures{aa, 1} );
    pkmax = [];
    pktime = [];
    auc = [];
    auc2 = []; %AUC from first interval not second

    for cc = 1:nr1
        sm_tr = peakfeatures{aa, 1}(cc, :);
        [~, idx, w, p] = findpeaks(sm_tr', 'MinPeakDistance',5 ,'MinPeakProminence', 0.1 );
        if isempty(idx)
            pkmax = [pkmax; 0];
            pktime = [pktime; 40];
        else
            pks = sortrows( [idx peakfeatures{aa, 1}(cc, idx )' ], 2, 'desc');
            hpktimetemp = find( peakfeatures{aa, 1}(cc, 1:39) < max(pks(:, 2))/2 & peakfeatures{aa, 1}(cc,2:40) > max(pks(:, 2))/2);
            if max(pks(:, 2)) < 0.5 || min(hpktimetemp) <2
                pkmax = [pkmax; 0];
                pktime = [pktime; NaN];
            else
                pkmax = [pkmax; max(pks(:, 2))];
                pktime = [pktime; hpktimetemp(1)];
            end

        end
        auc = [auc; trapz(peakfeatures{aa, 1}(cc, 1:40)) ];
        auc2 = [auc2; trapz(peakfeatures{aa, 1}(cc, 1:20)) ];
    end
    peakfeatures{aa, 2} = pkmax;
    peakfeatures{aa, 3} = auc2;
    peakfeatures{aa, 4} = max(auc, 0);
    peakfeatures{aa, 5} = pktime;

    if temp(1, 108) == 1 && temp(1, 106) < 16
        peakfeatures{aa, 6} = temp(1, 106)-1;
    else
        peakfeatures{aa, 6} = temp(1, 106);
    end
    peakfeatures{aa, 7} = temp(1, 107);
    peakfeatures{aa, 8} = temp(1, 108);
    peakfeatures{aa, 9} = temp(1, 109);
end

peakfeatures_sepsis = peakfeatures;
save('./single_cell_features/peakfeatures_sepsis', 'peakfeatures_sepsis' );


