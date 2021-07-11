function V = Diffusion_Limited_Time_Dependent_Drug_Fitting_Function(Params,tdata,Subject)
% Author: Adam Nasim,
%
% part of https://github.com/DrAdamNasim/Diffusion_Limited_Cancer_Growth_Model
% If using this or related code please cite 
% Nasim, A.; Yates, J.; Derks, G.; Dunlop, C. 
%     Mechanistic mathematical model of tumour growth and inhibition (diffusion-limited model)
%     (Manuscript submitted for publication).


tdata = tdata'; %Change input time data to column data 
%If not being used to fit to data the variables Outer_Volume and Plot_Time
%will produce smooth curves of the predicted tumour volume over the
%timeframe of the experiment.

%Data specific, way to specify dose timings depending on subject, e.g., 
% % % % if Subject<10
% % % %     Td = [1 8 15 tdata(end)];
% % % % elseif 10<=Subject && Subject<30
% % % %     Td = [1 8 15 22 tdata(end)];
% % % % else
% % % %     Td = [4 11 18 25 tdata(end)];
% % % % end
Td = [1 8 15 tdata(end)]; %Sufficient if dosing same for all subjects
Dose_Index = find(Td<=tdata(end));
Td = Td(Dose_Index);

%Set parameters
lambda = Params(1); %kg-kp Net growth rate
kq = Params(2);     %necrotic loss rate
Kill = Params(3);   %drug potency parameter
Vstar = Params(4);  %Critical volume
V0 = Params(5);     %Initial volume

CL = 1.39;
GF0 = 1;        %If exponential growth, if not see below
C_0 = 50; %Initial dose concentration
C = 0;
Ini_Cond = [V0;GF0];

%initialise arrays
Outer_Volume = [];
Plot_Time = [];
Concentration_Drug = [];
GF = [];
Vd = [];   
GFd = [];
Fitting_Time = [];
at = [];
%Switch to radial values to calculate initial GF conditions
r0 = (3*V0/(4*pi))^(1/3);
rstar = (3*Vstar/(4*pi))^(1/3);
% Find initial value of the growth fraction if initial volume greater than
% critical volume
if r0>=rstar
    PC = [2, -3, 0, 1-(rstar/r0).^2];%Phase 2 starting value r_T quiescent value if r1<r0<r3
    if all(isfinite(rstar/r0))
        Q0 = roots(PC);
        r01_2 = Q0(2)*r0;
        r04 = [r0; r01_2];
        GF0 = (r04(1)^3-r04(2)^3)/(r04(1)^3);%Initial GF when GF ~= 1
        Fl=1;
        Ini_Cond = [V0;GF0];
    else

    end
else
    Fl=0;
end

options1 = odeset('Events', @events1);
options2 = odeset('Events', @events2);
noDose=size(Td);
noDose=noDose(2);
n=1;
t=0;
m = 1;
 warning('') % Clear last warning message
 %While loop allows to solve model depending on phase (fully proliferating
 %or diffusion-limited) as well as take into account the dose timings
while m<noDose
    if  t(end)<Td(1) && Fl==0    
        Kill=0; % No dosing in this phase therefore no drug effect
        Sol1 = ode45(@Exp_Growth, [t(end),Td(1)], Ini_Cond(1),options1);
%         Outer_Volume = [Outer_Volume Sol1.y];
%         Plot_Time = [Plot_Time Sol1.x];
        %
        ie = Sol1.ie; % flag that tells us whether or not we change growth phase
        time_end = Sol1.x(end);
        timed1 = tdata(t(end)<=tdata & tdata<time_end);
        S1 = size(timed1);
        if S1(2)==0
        else
            Deval_time1 = timed1;
            Fitting_Time = [Fitting_Time Deval_time1];
            Vd = [Vd deval(Sol1,Deval_time1)];
            elementsd = size(Deval_time1);
            GFd = [GFd ones(1,elementsd(2))];
        end
        a = real(Sol1.y);
        at = Sol1.x;
        t = at;
        % Make rN artificially 0 so concatenation is possible
        elements = size(t); 
        GFa = ones(1,elements(2));
        GF = [GF GFa];
        Ini_Cond = [a(end);1];
    elseif t(end)<Td(1) && Fl==1
        Kill=0; % No dosing in this phase therefore no drug effect
        Sol2 = ode15s(@DB_fun,[t(end),Td(1)],Ini_Cond,options2);
        t2 = Sol2.x(1);
        t2end = Sol2.x(end);
        Outer_b = Sol2.y;
        Outer_b = Outer_b(1,:);
        Outer_Volume = [Outer_Volume Outer_b];
        Plot_Time = [Plot_Time Sol2.x];
        ie = Sol2.ie;
        timed = tdata(tdata>=t(end) & tdata<Td(1));
        S1 = size(timed);
        if S1(2) ==0
        else 
            Fitting_Time = [Fitting_Time timed];
            Vd_Array = deval(Sol2,timed);
            Vd = [Vd Vd_Array(1,:)];
            GFd = [GFd Vd_Array(2,:)];
        end
        bt = Sol2.x;
        GFb = real(Sol2.y(2,:));
        GF = [GF GFb];
        if t ==0
            t = bt;
        else
            t = [at bt];
        end
        Ini_Cond = [Outer_Volume(end) GF(end)];
    elseif t(end)>=Td(1) && Fl==0
        Kill = Params(3);
        if t(end)-Td(m) == Td(m+1)-Td(m)
        else
        Sol3 = ode45(@Exp_Growth, [t(end)-Td(m),Td(m+1)-Td(m)], Ini_Cond(1),options1);
        ie = Sol3.ie;
        Outer_c = Sol3.y;
        Outer_c = Outer_c(1,:);
        Outer_Volume = [Outer_Volume Outer_c];
        Outer_Volume = real(Outer_Volume);
        Plot_Time = [Plot_Time Sol3.x+Td(m)];
        %DEVAL
         time1 = Sol3.x(end)+Td(m);
         if m == noDose-1
            timed = tdata(tdata>=t(end) & tdata<=time1);
         else
            timed = tdata(tdata>=t(end) & tdata<time1);
         end
        if isempty(timed)  

        else
            Deval_time2 = timed-Td(m);
            Fitting_Time = [Fitting_Time Deval_time2+Td(m)];
            Vd = [Vd deval(Sol3,Deval_time2)];
            elementsd = size(Deval_time2);
            GFd = [GFd ones(1,elementsd(2))];
        end   
        time_t = Sol3.x;
        elements2 = size(time_t);
        GFc = ones(1,elements2(2));         
        GF = [GF GFc];
        ct = Sol3.x + Td(m);
        t = [t ct];
        Ini_Cond = [Outer_Volume(end);1];
        end
    elseif t(end)>=Td(1) && Fl==1
        Kill=Params(3);
        if t(end)-Td(m) == Td(m+1)-Td(m)
        else
            Sol4 = ode15s(@DB_fun,[t(end)-Td(m),Td(m+1)-Td(m)],Ini_Cond,options2);
            ie = Sol4.ie;
%             %For plotting
            Outer_d = Sol4.y;
            GF = [GF Outer_d(2,:)];
            Outer_Volume = [Outer_Volume Outer_d(1,:)];
            Outer_Volume = real(Outer_Volume);
            Plot_Time = [Plot_Time Sol4.x+Td(m)];
            %Deval
             time1 = Sol4.x(end)+Td(m);
             if m == noDose-1
                timed = tdata(tdata>=t(end) & tdata<=time1);
             else
                timed = tdata(tdata>=t(end) & tdata<time1);
             end
            if isempty(timed)  

            else
                Deval_time2 = timed-Td(m);
                Fitting_Time = [Fitting_Time timed];
                Vd_Array1 = deval(Sol4,Deval_time2); 
                Vd = [Vd Vd_Array1(1,:)];
                GFd = [GFd Vd_Array1(2,:)];
            end 
            dt = Sol4.x + Td(m);
            t = [t dt];  
            Ini_Cond = [Outer_Volume(end) GF(end)];
        end
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
        m = m+1;
        C_0 = C_0 + Concentration_Drug(end);
    end
     n=n+1;
     ie=[]; 
end

Check_Size_tdata = size(tdata);
Check_Size_tdata = Check_Size_tdata(2);
Check_Size_V = size(Vd);
Check_Size_V = Check_Size_V(2);
if Check_Size_tdata == Check_Size_V
    V = Vd;
else
    V = zeros(1,Check_Size_tdata);
end
V = Vd';
[warnMsg, warnId] = lastwarn;
if ~isempty(warnMsg)
    V = zeros(length(tdata),1);
end
 %plot(Plot_Time,Outer_Volume)
function Vol = Exp_Growth(t,x)
    C = C_0*exp(-CL*t);
    Concentration_Drug = [Concentration_Drug C];
    Deff = Kill*C;
    Vol = (lambda-Deff)*x;
end    
function f = DB_fun(t,x)
    x(2) = x(2);
    C = C_0*exp(-CL*t);
    Concentration_Drug = [Concentration_Drug C];
    Deff = Kill*C;
    Q = nthroot(1-x(2), 3);
    K_V = (lambda-Deff)*x(2)-kq*(1-x(2));
    f = [K_V.*x(1);-(1/3)*K_V*Q*(1+2*Q)*(1-Q)];
end

%Event solvers to allow switching between exponential and non-exp phases
    function [check, isterminal, direction] = events2(~,x)
       direction = -1;   % Negative direction only
      isterminal = 1;%terminate integration when event occurs
      check = double(x(1,:) > Vstar);
    end

    function [check, isterminal, direction] = events1(~,x)
        direction = [];
        isterminal = 1;%terminate integration when event occurs
        check = double(x(1,:) < Vstar); %terminate second phase when critical size is reached
    end
end