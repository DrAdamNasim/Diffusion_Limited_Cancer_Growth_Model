function Volume = Two_Phase_Spatial_Drug_Hyperbolic_Function(Params,tdata)

% Author: Adam Nasim,
%
% part of https://github.com/DrAdamNasim/Diffusion_Limited_Cancer_Growth_Model
% If using this or related code please cite 
% Nasim, A.; Yates, J.; Derks, G.; Dunlop, C. 
%     Mechanistic mathematical model of tumour growth and inhibition (diffusion-limited model)
%     (Manuscript submitted for publication).


%Plots drug concentration heatmaps as well as can
%track drug at depth throughout experiment (can be set by user).
%Currently heat maps only supported when in diffusion-limited phase at the
%time of plotting
%%
Dose_Strength = 50;
No_Runs = size(Dose_Strength);%Make vector to run for different dose strengths
No_Runs = No_Runs(2);
DPI = 60; %How clear to make heatmaps
F_Vec = [0.001 2.5 5 7.5 10]; %Spatial drug parameter. Heat maps will be plotted for these values
for j = 1:length(F_Vec)
    F = F_Vec(j);
    for i = 1:No_Runs
        Dose = Dose_Strength(i);
        Dose_i = Dose_Strength(i);
        Outer_Radius = [];
        Plot_Time = [];
        Td = [1 8 15 38];
        tdata = linspace(Td(1),Td(end));
        %Set parameters
        alpha = 1.39;
        E = sqrt(alpha/((1.8*10^(-5))*60*60*24));%sqrt(wash out rate of drug in tissue/Diffusion coefficient of drug (per day)) - Assume wash out the same as in plasma (ASSUMPTION)
        Params = [0.131 0.02 8e-3 0.333 0.413];
        L = Params(1); %Growth rate
        mu = Params(2); %Loss rate from necrotic region
        Kill = Params(3); %Drug potency
        rstar = Params(4); %Critical radius
        r0 = Params(5); %Initial radius
        w_ext = [];
        noDose=size(Td);
        noDose=noDose(2);
        n=1;
        t=0;
        YI = [];
        %If we start in a diffusion limited state r_0>rstar we find the initial
        %conditions using the algebraic constraint equaiton
        if r0>=rstar
            PC = [2, -3, 0, 1-(rstar/r0).^2];%Phase 2 starting value r__T quiescent value if r1<r0<r3
            Q0 = roots(PC);
            r01_2 = Q0(2)*r0;
            r04 = [r0; r01_2^2];
            Fl=1;
            YI = r04;
            GF = (r04(1)^3-r04(2)^3)/(r04(1)^3);%Initial GF to print only
        else
            Fl=0;
        end
        m=1;
        M1 = [1 0 ; 0 0];
        options1 = odeset('Events', @events1);
        options2 = odeset('Mass',M1,'Events', @events2);
        options3 = odeset('Events', @events3);
        r = [];
        rd = [];
        while m<noDose
            if  t(end)<Td(1) && Fl==0    
                Kill=0;%No dosing in this phase
                w_ext = 0;
                Sol1 = ode45(@(t,r) L.*r./3, [t(end),Td(1)], r0,options1);
                %For plotting outer radius/Volume
                Outer_Radius = [Outer_Radius Sol1.y];
                Plot_Time = [Plot_Time Sol1.x];
                %
                ie = Sol1.ie;
                time_end = Sol1.x(end);
                timed1 = tdata(tdata<time_end);
                S1 = size(timed1);
                if S1(2)==0
                else
                Deval_time1 = tdata(tdata<time_end);
                ad = deval(Sol1,Deval_time1);
                elementsd = size(Deval_time1);
                rQd = zeros(1,elementsd(2));
                rd = [ad;rQd];
                end
                a = real(Sol1.y);
                at = Sol1.x;
                t = at;
                %Make rN artificially 0 so concatenation is possible
                elements = size(t); %how many elements of zero we need
                rQ = zeros(1,elements(2));
                r = [a;rQ];
                YI=a(:,end);
                r04 = [YI;0];

            elseif t(end)<Td(1) && Fl==1
                Kill=0;%No dosing in this phase
                w_ext = 0;
                Sol2 = ode15s(@Phase2_DAE,[t(end),Td(1)],r04,options2);
                t2 = Sol2.x(1);
                t2end = Sol2.x(end);
                %For plotting outer radius/Volume
                Outer_b = Sol2.y;
                Outer_b = Outer_b(1,:);
                Outer_Radius = [Outer_Radius Outer_b];
                Plot_Time = [Plot_Time Sol2.x];
                ie = Sol2.ie;
                %DEVAL
                timed = tdata(tdata>=t(end) & tdata<Td(1));
                S1 = size(timed);
                if S1(2) ==0
                else 
                    Deval_time1 = timed;
                    rd = deval(Sol2,Deval_time1);
                end
                bt = Sol2.x;
                b = real(Sol2.y);
                if t ==0
                    t = bt;
                else
                    t = [at bt];
                end
                YI = b(:,end);
                r = [r b];

            elseif t(end)>=Td(1) && Fl==0
                Kill = Params(3);
                if Dose_i == 0
                    Kill=0;
                end
                if Td(m)==0
                    YI = r0;
                end
                Sol3 = ode45(@DifEqn, [t(end)-Td(m),Td(m+1)-Td(m)], YI(1),options3);
                ie = Sol3.ie;
                %For plotting
                Outer_c = Sol3.y;
                Outer_c = Outer_c(1,:);
                Outer_Radius = [Outer_Radius Outer_c];
                Plot_Time = [Plot_Time Sol3.x+Td(m)];
                %DEVAL
                 time1 = Sol3.x(end)+Td(m);
                 if m == 4
                    timed = tdata(tdata>=t(end) & tdata<=time1);
                 else
                    timed = tdata(tdata>=t(end) & tdata<time1);
                 end
                if isempty(timed)  

                else
                    Deval_time2 = timed-Td(m);
                    cd = deval(Sol3,Deval_time2);
                    elementsd = size(Deval_time2);
                    rQd = zeros(1,elementsd(2));
                    cd = [cd;rQd];
                    rd = [rd cd];
                end   
                time_t = Sol3.x;
                elements2 = size(time_t);
                rQ = zeros(1,elements2(2));
                ct = Sol3.x + Td(m);
                c = real(Sol3.y);
                c = [c;rQ];
                YI= c(:,end);
                t = [t ct];
                r = [r c];
               YI=YI(:,end);
            elseif t(end)>=Td(1) && Fl==1
                Kill=Params(3);
                if Dose_i == 0
                    Kill=0;
                end            
                Sol4 = ode15s(@Phase2_DAE,[t(end)-Td(m),Td(m+1)-Td(m)],YI,options2);
                ie = Sol4.ie;
                %For plotting
                Outer_d = Sol4.y;
                Outer_d = Outer_d(1,:);
                Outer_Radius = [Outer_Radius Outer_d];
                Plot_Time = [Plot_Time Sol4.x+Td(m)];
                %Deval
                 time1 = Sol4.x(end)+Td(m);
                 if m == 4
                    timed = tdata(tdata>=t(end) & tdata<=time1);
                 else
                    timed = tdata(tdata>=t(end) & tdata<time1);
                 end
                if isempty(timed)  

                else
                    Deval_time2 = timed-Td(m);
                    cd = deval(Sol4,Deval_time2);              
                    rd = [rd cd];
                end 
                dt = Sol4.x + Td(m);
                d = real(Sol4.y);
                YI=d(:,end);
                t = [t dt];
                r = [r d];      
                YI=YI(:,end);
            end
            %% toggle flag
            s=size(ie);
            s=s(1);
            if ie==1 %i.e. an event happened change flag
                if Fl==0
                    Fl=1; %adjust for the cases you set
                else
                    Fl=0;
                end
            elseif s==0 && t(end)==Td(1) %i.e. no event happened just got to end of first phase
                %dont do anything
            elseif s==0 && t(end)>Td(1) %i.e. just got to end of one phase
                m=m+1;
                Dose = Dose_i + w_ext(end);
            end
             n=n+1;
             ie=[]; 
        end
        R = r(1,:);
        RN = sqrt(r(2,:));%The sqrt is used due to the variable substitution see the Phase2_DAE function below
        GF_Size = size(R);
        size(RN);
        GF_Vec = ones(1,GF_Size(2));
        GF = GF_Vec - (RN./R).^3;
        Volume = (4/3)*pi*Outer_Radius.^3;
        i = i+1;
    end
    %Define which time point to evaluate the drug concentration at
    Day_of_Dose = Td(3);
    largeEnough = t >= Day_of_Dose;
    smallEnough = t <= Day_of_Dose+1;
    result1 = t(largeEnough & smallEnough);
    m_I = find(t == result1(1));
    m_I = m_I(1);
    Extra_Space = R(m_I)/10;
    lower = -R(m_I)-Extra_Space;
    upper = R(m_I)+Extra_Space;
    X = linspace(lower,upper,DPI);
    Y = linspace(lower,upper,DPI);
    Tumour_Radius = R(m_I);
    WW=zeros(DPI,DPI);
    W_ext_after_1st_Dose = Dose_i*(exp(-alpha.*(t(m_I)-Day_of_Dose)));
    w_ext = W_ext_after_1st_Dose;
    % Here we use the hyperbolic form of the drug concentration w, in the functions below we
    % choose to use Bessel functions instead. Both are equivalent.
    M = (E.*cosh(E*RN(m_I)).*cosh(F*RN(m_I))-F.*sinh(F*RN(m_I)).*sinh(E.*RN(m_I)))./(E*cosh(E*RN(m_I)).*sinh(F*RN(m_I))-F*cosh(F*RN(m_I)).*sinh(E*RN(m_I)));
    A = -M*F*R(m_I)./(cosh(F*R(m_I))-M*sinh(F*R(m_I)));
    B = -A./M;
    A_Hat = (A.*sinh(F*RN(m_I))+B.*cosh(F.*RN(m_I)))./(sinh(E.*RN(m_I)));
    for HH=1:DPI
       for YY = 1:DPI
            Radius = sqrt(X(HH).^2+Y(YY).^2);
            if Radius<=RN(m_I)
                    WW(HH,YY) = A_Hat.*w_ext.*sinh(E.*Radius)./(F.*Radius);
            elseif Radius<=R(m_I) && Radius>=RN(m_I)
                    WW(HH,YY) = A.*w_ext.*sinh(F*Radius)./(F*Radius) + B.*w_ext.*cosh(F*Radius)./(F*Radius);
            else
                WEXT = max(WW);
                WVAL = max(WEXT);
                WW(HH,YY) = WVAL;
                WW(HH,YY) = W_ext_after_1st_Dose;
            end 
       end

    end
    Radius_Drug_Con = linspace(0,0.3,DPI);
    for DD=1:DPI
        Radius = Radius_Drug_Con(DD);
        if Radius<=rstar
            R_TT = Radius_Drug_Con(end);
                Drug_Con(DD) = (w_ext.*(R_TT)./(besseli(0.5,F.*R_TT))).*(besseli(0.5,F.*Radius)./sqrt(Radius));
        elseif Radius<=R(m_I) && Radius>=RN(m_I)
                Drug_Con(DD) = sqrt(pi./(2.*F.*Radius)).*((C_t.*besseli(0.5,F.*Radius) + D_t.*besselk(0.5,F.*Radius)));
        else
            DEXT = max(WW);
            WVAL = max(DEXT);
            Drug_Con(HH,YY) = WVAL;
            Drug_Con(HH,YY) = W_ext_after_1st_Dose;
        end 
    end
    figure(length(F_Vec)+1)
    plot(Radius_Drug_Con,Drug_Con./max(Drug_Con).*100,'LineWidth',3)
    legend('F=0','F = 2.5','F = 5','F = 7.5','F = 10')
    xlabel('Distance From Tumour Centre (cm)')
    ylabel('Drug Concentration (%)')
    xlim([0 .3])
    set(gca,'FontSize',60,'XDir','reverse')
    box off
    legend box off
    hold on

    %Heatmap plots
    Necrotic_Radius = RN(m_I);
    figure(j)
    WW_Max_Vec = max(WW);
    WW_Max = WW_Max_Vec(1);
    hold on
    surf(X, Y, WW./WW_Max.*100, 'edgecolor', 'none','FaceColor','interp');
    hold on
    colorbar()
    h = colorbar;
    %For necrotic overlay on heatmap
    circle(0,0,Necrotic_Radius,'--',[0.65 0.65 0.65])
    circle(0,0,Tumour_Radius,'-',[0 0 0])
    set(get(h,'label'),'string','Drug Concentration (%)');
    title(['F = ',num2str(F),''])
    set(gca,'FontSize',60)
    xlim([lower,upper])
    ylim([lower,upper])
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    caxis([10 100])
    axis square
    
    %Plot solution 
    figure(length(F_Vec)+2)
    yyaxis left
    plot(Plot_Time,Volume,'LineWidth',2)
    ylabel('Tumour Volume (cm^3)')
    yyaxis right
    plot(Plot_Time,GF*100,'LineWidth',3)
    xlabel('Time(days)')
    ylabel('Growth fraction %')
    hold on
end
set(gca,'FontSize',60)
box off
%%%%%%%%%%%%%%%____________________%%%%%%%%%%%%%%%
function f = Phase2_DAE(t,r)
%We use a change of variables R_N->sqrt(R_N) (abusing notation) to
%allow us to satisfy the Index 1 property of algebraic differential equations
R_T = r(1);
R_N = r(2);
w_ext = Dose*exp(-alpha.*t);
if r(2) == 0
    s = linspace(0,R_T,50);
    w = s.^(3/2).*Kill.*w_ext*sqrt(R_T).*besseli(1/2,F.*s)./(besseli(1/2,F.*R_T));
    f = L.*R_T./3 - (1./(R_T.^2)).*trapz(s,w);
else
    M = (F*besseli(3/2,F*sqrt(abs(R_N))).*besseli(1/2,E*sqrt(abs(R_N)))-E*besseli(1/2,F*sqrt(abs(R_N))).*besseli(3/2,E*sqrt(abs(R_N))))./((E*besselk(1/2,F*sqrt(abs(R_N))).*besseli(3/2,E*sqrt(abs(R_N)))+F*besselk(3/2,F*sqrt(abs(R_N))).*besseli(1/2,E*sqrt(abs(R_N)))));            
    s = linspace(sqrt(abs(R_N)),R_T,50);
    w = s.^2.*Kill.*w_ext.*sqrt(R_T).*(besseli(1/2,F.*s)+M.*besselk(1/2,F.*s))./(sqrt(s).*(besseli(1/2,F.*R_T)+M.*besselk(1/2,F.*R_T)));
    f =  [(L.*(R_T.^3-sqrt(abs(R_N)).^3)-mu.*(sqrt(abs(R_N))).^3)./(3.*R_T.^2) - (1./(R_T.^2)).*trapz(s,w);
    -(rstar./r(1)).^2 + (1 + 2.*sqrt(abs(R_N))./r(1)).*(1 - sqrt(abs(R_N))./r(1)).^2];
end
end

function f = DifEqn(t,r) 
 R_T = r;
 w_ext = Dose*exp(-alpha.*t);
 s = linspace(0,R_T,50);
 w = s.^(3/2).*Kill.*w_ext*sqrt(R_T).*besseli(1/2,F.*s)./(besseli(1/2,F.*R_T));
 f = L.*R_T./3 - (1./(R_T.^2)).*trapz(s,w);
end

    function [check, isterminal, direction] = events2(~,r)
       direction = -1;   % Negative direction only
      isterminal = 1;%terminate integration when event occurs
      check = double( r(1,:) > rstar);
    end

    function [check, isterminal, direction] = events1(~,r)
        direction = [];
        isterminal = 1;%terminate integration when event occurs
        check = double(r(1,:) < rstar); %terminate second phase when critical size is reached
    end
    function [check, isterminal, direction] = events3(~,r)
        direction = [];
        isterminal = 1;%terminate integration when event occurs
        check = double(r(1,:) > rstar); %terminate second phase when critical size is reached
    end
function circle(x,y,r,line_style,colour)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle, 0.01 is the anglular step
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
ZL = length(xp);
z = repmat(100,1,ZL);
plot3(x+xp,y+yp,z,'Linewidth',3,'LineStyle',line_style,'Color', colour);
end
end
