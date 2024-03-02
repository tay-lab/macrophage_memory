function [y_resall, y_resall2, y_resall3, y_resall4] = E16_Sol2_plotter_awang(TNF, TNF4ini, ParIdx,NumPar)

[~,nc1]= size(TNF4ini);
tspan1 = 0:1:nc1-1;
nTime1 = length(tspan1);
[nr,nc]=size(TNF);
tspan = 0:1:nc-1; % mins
nTime     = length(tspan);
y_res = [];
y_resall = [];
y_res2 = [];
y_resall2 = [];
y_res3 = [];
y_resall3 = [];
y_res4 = [];
y_resall4 = [];
%ParA
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

    parfor bb = 1:3
        ParP = ParIdx';
        foo = [.6, 1, 1.8]

        %initial conditions
        y0 = [0, 0, 0, 10, 0, ParP(23)*0.15, ParP(23), 0, 0, 0, 0, 0, 1, 0, 0, 0];

        ini = y0;
        y_list = zeros(nTime, length(y0));
        %y0(5) = y0(5)*10^((bb-2)*.5);
        y0(6) = y0(6)*foo(bb);
        y0(7) = y0(7)*foo(bb);
        y_list(1, :) = y0;


        for iTime = 2:nTime
            [t_, y_] = ode15s(@(t,y) E16_Mod2_awang(t,y,TNF(iTime),ParP,NumPar), tspan(iTime-1:iTime), y0);
            y0       = y_(end, :);
            y_list(iTime, :) = y0;

            if iTime == 721
                y0(14) = StimPar(1);
                y0(15) = StimPar(2);
                y0(16) = StimPar(3);
                %ParP(22) = 0;
            elseif iTime == 960
                y0(14) = 1;
                y0(15) = 0;
                y0(16) = 0;
            end
        end
         y_res(bb,:) = y_list(:,1) - y_list(1,1);
         y_res2(bb,:) = y_list(:,6)  + y_list(:,7);
         y_res3(bb,:) =  y_list(:,11); %ParP(19)*y_list(:,11)./10;
         y_res4(bb,:) = y_list(:,12) + 1000; %6000-y_list(:,9)-y_list(:,10); %
    end

    y_resall = [y_resall;y_res];
    y_resall2 = [y_resall2;y_res2];
    y_resall3 = [y_resall3;y_res3];
    y_resall4 = [y_resall4;y_res4];


end


end