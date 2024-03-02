function dy = E16_Mod2_awang(t,y,TNF,Par,FPar)

dy = zeros(16,1);

dy(1)= Par(4)*y(2) - Par(5)*y(1) - Par(10)*FPar(1)*y(5)*y(1);                                             % NFkBn
dy(2)= -Par(4)*y(2) + Par(5)*y(1) - Par(10)*y(2)*y(6) + Par(17)*y(9)*y(7); % +Par(17)*y(7);               % NFkBc

%dy(3)= Par(13)*(y(1)^FPar(2))/((Par(11)*FPar(3))^FPar(2)+y(1)^FPar(2)) - Par(14)*y(3);                    % mRNA induced
dy(3)= Par(13)*(y(1)^FPar(2))/((Par(11)*35000)^FPar(2)+y(1)^FPar(2)) - Par(14)*y(3);                    % mRNA induced

dy(4)= Par(12) - Par(14)*y(4);%+ Par(13)*(y(1)^FPar(2))/((Par(11)*FPar(3))^FPar(2)+y(1)^FPar(2))          % mRNA ikba

dy(5)= -Par(10)*FPar(1)*y(2)*y(1)+Par(6)*y(6)-Par(7)*y(5);                                                 % IkBn
dy(6)= -Par(10)*y(6)*y(2) + Par(16)*y(4) -Par(6)*y(6) +Par(7)*y(5) -Par(18)*y(6) + Par(16)*y(3);                         % IkBc

dy(7)= Par(10)*y(6)*y(2) -Par(17)*y(9)*y(7) +Par(9)*y(8) -Par(8)*y(7);                                    % NFkB:IkBc
dy(8)= Par(10)*FPar(1)*y(5)*y(1) -Par(9)*y(8) +Par(8)*y(7);                                                % NFkB:IkBn

%IKK terms NEW
dy(9)= Par(1)*y(15)*(FPar(5)-y(9)-y(10))*Par(22)^2/(y(11)^2+Par(22)^2) + ...                       % IKKa CpG component 
    Par(1)*y(16)*(FPar(5)-y(9)-y(10))*Par(22)^2/(y(11)^2+Par(22)^2) + ...                         % IKKa PIC component
    Par(1)*y(14)*(1 + Par(21)*y(12))*(FPar(5)-y(9)-y(10))*Par(22)^2/(y(11)^2+Par(22)^2) - ...       % IKKa TNF component
    Par(2)*y(9) ;                                                                                         % IKKa inactivation

% dy(10)=Par(2)*y(9)*(y(16)^2+Par(20)^2)/Par(20)^2 - Par(5)*y(10); % IKKi
dy(10)=Par(2)*y(9) - Par(3)*y(10);                                                                        % IKKi

dy(11)= Par(15)*y(3)-Par(20)*y(11);                                                                       % A20 

dy(12)=Par(15)*max( sign(y(16)), 0)*y(3)-Par(20)*y(12);                                                   % Additional TNFR

dy(13)= Par(15)*max( sign(y(16)), 0)*y(3)-Par(20)*y(13);                                                  % cIAP1/2

dy(14)=-Par(19)*y(14);                                                                                    % TNF

dy(15)=-Par(19)*y(15);                                                                                    % CpG

dy(16)=-Par(19)*y(16);                                                                                    % PIC



end