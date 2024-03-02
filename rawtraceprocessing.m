%% trace compilation ALL FILES
clc 
clear

%%%MODIFY THIS SECTION%%%
loadfolder = './rawtraces/';
savefolder = './processedtraces';

% savename = 'PICCPGfirst_traces.mat';
% savename = 'TNFLPSfirst_traces.mat';
% savename = 'PAMIL1first_traces.mat';
% savename = 'CPGLPSscrep_traces.mat';
% savename = 'PICTNFscrep_traces.mat';
% savename = 'PICfirst_TNFblock_traces.mat';
savename = 'CPGfirst_TNFblock_traces.mat';
% savename = 'CPGfirst_IL10block_traces.mat';
% savename = 'sepsis_traces.mat';

% rawfiles = {'20210825_FirstPICCPG_rep1.mat', '20210908_FirstPICCPG_rep2.mat'};
% rawfiles = {'20210824_FirstLPSTNF_rep1.mat', '20210907_FirstLPSTNF_rep2.mat'};
% rawfiles = {'20211006_FirstPAMIL1_rep1.mat', '20211013_FirstPAMIL1_rep2.mat'};
% rawfiles = {'20220407_CPGLPSscrep_traces.mat'};
% rawfiles = {'20220427_PICTNFscrep_traces.mat'};
% rawfiles = {'20220914_Cytoblock_rep1.mat', '20220920_Cytoblock_rep3.mat', '20220921_Cytoblock_rep4.mat', '20220927_Cytoblock_rep5.mat', '20220928_Cytoblock_rep6.mat'};
rawfiles = {'20220914_Cytoblock_rep1.mat', '20220920_Cytoblock_rep3.mat'};
% rawfiles = {'20220915_Cytoblock_rep2.mat', '20220921_Cytoblock_rep4.mat'};
% rawfiles = {'20220720_sepsis_rep3.mat', '20220707_sepsis_rep2.mat', '20220705_sepsis_rep1.mat'};

scmat = [];

counter = 1;
for cc = 1: length(rawfiles)


    load([loadfolder rawfiles{cc}]);
    [nrow, ~] = size(R);

    for aa = 1:nrow
        if ~isempty(R{aa, 15})
            startframe = round( R{aa, 15}(1, 1) )-3;
            endframe = startframe + 126;
            [nrow2, ~] = size(R{aa, 2});
            for bb = 1:nrow2
                if length(R{aa, 2}{bb, 2}(:, 1)) > 1
                    temp = NaN(1, length(startframe:(endframe+4)));
                    celltrace = R{aa, 2}{bb, 2}(:, [1, 9]);
                    if length(celltrace(:,1)) < length(celltrace(1,1):celltrace(end,1))
                        celltrace = [ [celltrace(1,1):celltrace(end,1)]',  interp1(celltrace(:, 1), celltrace(:,2), celltrace(1,1):celltrace(end,1) )'];
                    end
                    if celltrace(1, 1) >= startframe
                        startpos = celltrace(1,1);
                    else
                        startpos = startframe;
                    end
                    if celltrace(end, 1) <= endframe
                        endpos = celltrace(end, 1);
                    else
                        endpos = endframe;
                    end
                    temp(1, (startpos-startframe+1):(endpos-startframe+1)) = celltrace(find(celltrace(:,1) == startpos):find(celltrace(:,1) == endpos), 2)';
                    if isempty(scmat)
                        scmat = temp;
                    else
                        scmat = [scmat; temp];
                    end
                end

            end
            if counter == 1 && cc == 1
                lastcell = 0;
            end

            ligs = R{aa, 15}(:,2);

            scmat((lastcell+1):end,(end-4)) = ligs(1);
            scmat((lastcell+1):end,(end-3)) = ligs(end);
            if length( unique(ligs) ) < length( ligs)
                scmat((lastcell+1):end,end-2) = 8;
            else
                scmat((lastcell+1):end,end-2) = 4;
            end
                scmat((lastcell+1):end,end-1) = cc;
            scmat((lastcell+1):end,end) = counter;
            aa
            lastcell=length(scmat(:, 1));
            counter = counter + 1;
        end
    end

end

% save([savefolder savename], 'scmat' );
