function Vol_DL = Diffusion_Limited_Growth_Function(PHI,tdata,V)
%Solves the diffusion-limited growth model in the absence of drug input.

%Solver written in terms of the radial formulation, returns tumour volume
%given inputs parameters (PHI), time data (tdata) and the subject (leave blank if not fitting to multiple subjects. Functionality for nonlinear mixed effects fitting)

%Parameters
lambda = PHI(1);
mu = PHI(2);
Vstar = PHI(3);
V0 = PHI(4);
r0 = (3*V0./(4*pi)).^(1/3);
rstar = (3*Vstar./(4*pi)).^(1/3);
tswitch = 3*(log(rstar/r0)/lambda);

%Find time at which to solve exponential and diffusion-limited phases (if
%applicable)
index1 = find(tdata<tswitch);
tt=tdata(index1);
tt = tt';
index2 = find(tdata>=tswitch);
asizer=size(index2);
ttt=tdata(index2);
ttt = ttt';
if asizer(1)==0
    P=0;
    R2 = [];
    P1 = [];
else
    
    if tswitch >= 0
        P0 = 0;
        tstart = tswitch;
    else
        tstart = 0;
        PC = [2, -3, 0, 1-(rstar/r0).^2];
        Q0 = roots(PC);
        P0 = (Q0(2)).^2;
    end
    sol = ode45(DifEqn, [tstart tdata(end)], P0);
    P1 = sol.y;
    P=deval(sol,ttt);
    R2 = rstar.*sqrt(1./((1+2.*P.^(1/2)).*(1-P.^(1/2)).^2)); % P = Q^2 where Q = rN/rT 
end

R1=r0*exp(lambda.*tt./3);
TT = [tt ttt];
% convert to volume
R1 = (4/3)*pi*R1.^3;
R2 = (4/3)*pi*R2.^3;
Vol_DL = [R1 R2]';%Volume
function f = DifEqn(T,P)
        f = @(t,P) [2*lambda.*(1-((1+mu/lambda).*(sqrt(P)).^3)).*(1-sqrt(P)).*(1+2.*sqrt(P))./9];
end

end