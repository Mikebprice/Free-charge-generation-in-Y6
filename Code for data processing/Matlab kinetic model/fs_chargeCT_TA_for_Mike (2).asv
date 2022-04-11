clear all
time_raw    =importdata('Time_data_kinetics for Paul.xlsx');
A           =importdata('Exciton kinetics for Paul.xlsx');
B           =importdata('Charge kinetics for Paul.xlsx');
time        =time_raw.*10^-12;
extra_time  =[(-10^-12):(10^-13):(-2*10^-13)]';
extra_time  =extra_time.*ones(length(extra_time),4);
time        =[extra_time; time];
A           =[zeros(size(extra_time));A];
B           =[zeros(size(extra_time));B];
fluences    =[0.5 2.5 10 60]'*10^17;
n_time_points=length(time);
linewidth=1;

%% Start and Finish Times for Fitting
start_time  =200*10^-15;
end_time    =1600*10^-12;
ni          =zeros(1,4);
nf          =zeros(1,4);
NA          =zeros(1,4);
NB          =zeros(1,4);

for fluence=1:4
    for i=1:length(time)
        if  abs(time(i,fluence)-start_time)==min(abs(time(:,fluence)-start_time))
            ni(1,fluence)=i;
        end
        if  abs(time(i,fluence)-end_time)==min(abs(time(:,fluence)-end_time))
            nf(1,fluence)=i;
        end
    end
end
% nf=[length(A) length(B) 102 128];


fluences    =[0.5 2.5 10 60]'*10^17;
n_time_points=length(time);
linewidth=1;

%% Exciton, Total Charge Normalisation
A           =A./[0.9 1   1   1];
B           =B./(max(B).*[0.84 0.9 0.8 1]);

%% Check Kinetics
for fluence     =1:3
    fig
    timeindices =ni(fluence):nf(fluence);
    timerange   =time(timeindices,fluence);
    semilogx(timerange,A(timeindices,fluence),'LineWidth',1)
    hold on
    semilogx(timerange,B(timeindices,fluence),'LineWidth',1)
end

%% Integrate Rate Equations
na=1; nb=1; nd=1; ne=1;
fluences=[0.5 2.5 10 60]'*10^17;
n_fluences=3;
pb = CmdLineProgressBar('Iterations:');
best_global_err=inf;
n_time_points=length(A);
S_temp      =zeros(n_time_points,1);    S      =zeros(n_time_points,1);    
CT_temp     =zeros(n_time_points,1);    CT     =zeros(n_time_points,1);
C_temp      =zeros(n_time_points,1);    C      =zeros(n_time_points,1);
Ctotal_temp =zeros(n_time_points,1);    Ctotal =zeros(n_time_points,1);

% fig
% semilogx(time(:,fluence),B(:,fluence),'LineWidth',linewidth)
% hold on
for a=1:na
    for b=1:nb
        for d=1:nd
            for e=1:ne
                c0   = [0;0;0];               k=zeros(1,8);
                krad = 4*10^8;
                knr  = 1/(260*10^-12)-krad;
                kradnr=krad + knr;            k(1)=kradnr;
                kct  = 7.2*10^12;               k(2)=kct;
                kcr  = 4.2*10^12;               k(3)=kcr;
                knrct= 0*10^-10;              k(4)=knrct;
                kb   = 1.36*10^-7;            k(5)=kb;
                kcs  = 4.1*10^9;                k(6)=kcs;
                kenc = 1.95*10^-8;             k(7)=kenc;
                knrc = 1*10^8;                k(8)=knrc;
                err=inf(n_fluences,1);
                opts = odeset('RelTol',1e-5);
                for fluence=1:n_fluences
                    timeindices =ni(fluence):nf(fluence);
                    timerange   =time(timeindices,fluence);
                    i0=fluences(fluence); tspan=time(:,fluence);
                    [t,c]=ode15s(@(t,c) diffun_chargesCT_200fs(t,c,k,i0),tspan(1:n_time_points),c0,opts);
                    
                    S_temp(:,fluence) =c(:,1);
                    CT_temp(:,fluence)=c(:,2);
                    C_temp(:,fluence) =c(:,3);
                    Ctotal_temp(:,fluence)=C_temp(:,fluence)+CT_temp(:,fluence);
              
%%%                    
%                     Fit to CT + Charge
                    Ns  = max(S_temp(:,fluence)); 
                    Nc  = max(Ctotal_temp(:,fluence));
                    
                    combined=[S_temp(timeindices,fluence)./Ns,...
                              Ctotal_temp(timeindices,fluence)./Nc];
                    dM=combined-[A(timeindices,fluence),...
                                 B(timeindices,fluence)];
                    err(fluence) = sqrt(sum(sum(dM.^2)));
                end
                    global_err = norm(err);
                if  global_err < best_global_err
                    best_global_err=global_err;
                    bestk=k;
                    bestc=c;
                    S=S_temp; CT=CT_temp; C=C_temp; Ctotal=Ctotal_temp;
                end
            end 
        end
    %semilogx(time(1:n_time_points,1),Ctotal_temp(:,fluence)./Nc,'LineWidth',linewidth)
    end
    pb.print(a,na)
end
% legend

%% Plot Exciton Fits
fig
for fluence=1:3
    fig
    semilogx(time(1:n_time_points,fluence),S(:,fluence)./(max(S(:,fluence)*((fluence~=1)+(fluence==1)*0.9))),'LineWidth',linewidth)
    hold on
    semilogx(time(:,fluence),A(:,fluence),'LineWidth',linewidth)
    xlim([time(ni(fluence),fluence) time(nf(fluence),fluence)])
end

%% Plot Total Charge Fits
for fluence=1:n_fluences
    fig
    semilogx(time(1:n_time_points,fluence),Ctotal(:,fluence)./(max(Ctotal(:,fluence))),'LineWidth',linewidth)
    hold on
    semilogx(time(:,fluence),B(:,fluence),'LineWidth',linewidth)
    xlim([time(ni(fluence),fluence) time(nf(fluence),fluence)])
end

%% Plot Exciton/CT/Charge Fits
for fluence=1:3
    fig
    semilogx(time(1:n_time_points,fluence),S(:,fluence),'LineWidth',linewidth)
    hold on
    semilogx(time(1:n_time_points,fluence),CT(:,fluence),'LineWidth',linewidth)
    semilogx(time(1:n_time_points,fluence),C(:,fluence),'LineWidth',linewidth)
%     semilogx(time(1:n_time_points,fluence),CT(:,fluence)./S(:,fluence),'LineWidth',linewidth)
end
