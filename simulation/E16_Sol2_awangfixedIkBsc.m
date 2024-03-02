function [y_resall, y_resall2] = E16_Sol2_awangfixedIkBsc(TNF, TNF4ini, ParIdx, NumPar)

[~,nc1]= size(TNF4ini);
tspan1 = 0:1:nc1-1;
nTime1 = length(tspan1);
[nr,nc]=size(TNF);
tspan = 0:1:nc-1; % mins
nTime     = length(tspan);
y_res = [];
y_resall = [];
y_resall2 = [];
for aa = 1:3
    if aa == 1
        StimPar = [0, 0, 0];
    elseif aa == 2
        StimPar = [0, 1, 0];
    else
        StimPar = [0, 0, 1];
    end
    %%% Initial Condition %%%

    
    % y_list    = zeros(nTime1, length(y0));
    % y_list(1, :) = y0;
    % for iTime = 2:nTime1
    %     [t_, y_] = ode15s(@(t,y) E10_Mod0(t,y,TNF4ini(1,iTime),ParP,NumPar), tspan1(iTime-1:iTime), y0);
    %     y0       = y_(end, :);
    %     y_list(iTime, :) = y0;
    % end
    % ini = y_list(end,:);


    ParP = ParIdx';
    y0 = [0, 0, 0, 100, 0, ParP(23)*.15, ParP(23), 0, 0, 0, 0, 0, 1, 0, 0, 0];
    ini = y0;
    y_list = zeros(nTime, length(y0));
    y_list(1, :) = y0;


    for iTime = 2:nTime
        [t_, y_] = ode15s(@(t,y) E16_Mod2_awang(t,y,TNF(iTime),ParP,NumPar), tspan(iTime-1:iTime), y0);
        y0       = y_(end, :);
        y_list(iTime, :) = y0;

        if iTime == 480
            y0(14) = StimPar(1);
            y0(15) = StimPar(2);
            y0(16) = StimPar(3);
            %ParP(22) = 0;
        elseif iTime == 720
            y0(14) = 1;
            y0(15) = 0;
            y0(16) = 0;
        end
    end
    y_res = y_list(:,1) - y_list(1,1);

    y_res2 = y_list(:,6)  + y_list(:,7);

    y_resall = [y_resall;y_res'];

    y_resall2 = [y_resall2; y_res2' ];


end


end