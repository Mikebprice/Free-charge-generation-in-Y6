clear all
time_raw    =importdata('Time_data_kinetics for Paul.xlsx');
A           =importdata('Exciton kinetics for Paul.xlsx');
B           =importdata('Charge kinetics for Paul.xlsx');
PL          =importdata('Y6_combined.xlsx');
PL_fluences_all = PL(:,1);
PL_fluences = [PL_fluences_all(1) PL_fluences_all(7)];
PL_int      = PL(:,2);
time        =time_raw.*10^-12;
extra_time  =[(-10^-12):(10^-13):(-2*10^-13)]';
extra_time  =extra_time.*ones(length(extra_time),4);
time        =[extra_time; time];
A           =[zeros(size(extra_time));A];
B           =[zeros(size(extra_time));B];
datarange   =1:2;
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
% for fluence     =1:3
%     fig
%     timeindices =ni(fluence):nf(fluence);
%     timerange   =time(timeindices,fluence);
%     semilogx(timerange,A(timeindices,fluence),'LineWidth',1)
%     hold on
%     semilogx(timerange,B(timeindices,fluence),'LineWidth',1)
% end

%% Integrate Rate Equations
na=3; nb=3; nc=5; nd=3; ne=3;
kct_range = 10*10^12; kcr_range = 10*10^12; kcs_range = 10*10^9; knrc_range = 10*10^8; kenc_range = 10*10^(-8);
fluences=[0.5 2.5 10 60]'*10^17;
n_fluences=3;
pb = CmdLineProgressBar('Iterations:');
best_global_err=inf;
n_time_points=length(A);
S_temp      =zeros(n_time_points,1);    S      =zeros(n_time_points,1);    
CT_temp     =zeros(n_time_points,1);    CT     =zeros(n_time_points,1);
C_temp      =zeros(n_time_points,1);    C      =zeros(n_time_points,1);
Ctotal_temp =zeros(n_time_points,1);    Ctotal =zeros(n_time_points,1);
sun = [1.9401*10^22,1.9401*10^25];
global_err = zeros(na,nb,nd,ne);
error_to_plot = [0];
C_frac_to_plot = [0];
C_frac_to_plot2 = [0];
% fig
% semilogx(time(:,fluence),B(:,fluence),'LineWidth',linewidth)
% hold on
for a=1:na
    for b=1:nb
        for cc=1:nc
            for d=1:nd
                for e=1:ne
                    c0   = [0;0;0];               k=zeros(1,8);
                    krad = 2.3*10^8;
                    knr  = 1/(260*10^-12)-krad;
                    kradnr=krad + knr;            
                    k(1) =kradnr;
                    kct  = 7.2*10^12-.5*kct_range+a*kct_range/na;               k(2)=kct;
                    kcr  = 4.2*10^12-.5*kcr_range+b*kcr_range/nb;               k(3)=kcr;
                    knrct= 0*10^9;                k(4)=knrct;
                    kb   = 1.36*10^-7;            k(5)=kb;
                    kcs  = (3.9)*10^9-.5*kcs_range+d*kcs_range/nd;                k(6)=kcs;
                    kenc = (3.95)*10^(-8)-.5*kenc_range+cc*kenc_range/nc;             k(7)=kenc;
                    knrc = (2.5)*10^8-.5*knrc_range+e*knrc_range/ne;                k(8)=knrc;
                    err=inf(n_fluences,1);
                    opts = odeset('RelTol',1e-5);
                     for fluence=datarange
                        i0 =PL_fluences(fluence);
                        [t,c]=ode15s(@(t,c) diffun_chargesCT_600ps(t,c,k,i0),[0:(10^-12):1*10^-8],c0,opts);
                        PLQEraw(fluence)=trapz(t,krad*c(:,1))/i0;
    %                     S_temp(:,fluence) =c(:,1);
    %                     CT_temp(:,fluence)=c(:,2);
    %                     C_temp(:,fluence) =c(:,3);
    %                     
    %                     %C_fraction_temp(fluence) = C_temp(1200,fluence)./(C_temp(1200,fluence)+CT_temp(1200,fluence)+S_temp(1200,fluence));
    %                     Ctotal_temp(:,fluence)=C_temp(:,fluence)+CT_temp(:,fluence);
    %                     Ctotal_temp(:,fluence)=Ctotal_temp(:,fluence);
                     end

                     if PLQEraw(2)/PLQEraw(1)>1.1
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
                        global_err(a,b,d,e) = norm(err);

                        for n_sun=1:2
                            i0=sun(n_sun); tspan=[0:(10^-12):1*10^-8];
                            [t,c]=ode15s(@(t,c) diffun_chargesCT_steady(t,c,k,i0),tspan,c0,opts);
                            S_temp_stdy(n_sun) =c(end,1);
                            CT_temp_stdy(n_sun)=c(end,2);
                            C_temp_stdy(n_sun) =c(end,3);

                            C_fraction(n_sun) = C_temp_stdy(n_sun)/(S_temp_stdy(n_sun)+CT_temp_stdy(n_sun)+C_temp_stdy(n_sun));
                            C_fraction2(n_sun) = 2*C_temp_stdy(n_sun)/(S_temp_stdy(n_sun)+CT_temp_stdy(n_sun)+2*C_temp_stdy(n_sun));
                        end
                        error_to_plot = [error_to_plot norm(err)];
                        C_frac_to_plot = [C_frac_to_plot C_fraction(n_sun)];
                        C_frac_to_plot2 = [C_frac_to_plot2 C_fraction2(n_sun)];
    %                     figure(1)
    %                     plot(norm(err),C_fraction(1))
    %                     hold on
    %                     figure(2)
    %                     plot(norm(err),C_fraction(2))
    %                     hold on
    %                     figure(3)
    %                     plot(norm(err),C_fraction2(1))
    %                     hold on
    %                     figure(4)
    %                     plot(norm(err),C_fraction2(2))
    %                     hold on
    %                     if  global_err < best_global_err
    %                         best_global_err=global_err;
    %                         bestk=k;
    %                         bestc=c;
    %                         S=S_temp; CT=CT_temp; C=C_temp; Ctotal=Ctotal_temp;
    %                     end
                     end
                 end
            end 
        end
    %semilogx(time(1:n_time_points,1),Ctotal_temp(:,fluence)./Nc,'LineWidth',linewidth)
    end
    pb.print(a,na)
end
% legend
fig
plot(C_frac_to_plot(2:end),error_to_plot(2:end),'.')

%% Plot Exciton Fits
fig
for fluence=1:3
    fig
    semilogx(time(1:n_time_points,fluence),S_temp(:,fluence)./(max(S_temp(:,fluence)*((fluence~=1)+(fluence==1)*0.9))),'LineWidth',linewidth)
    hold on
    semilogx(time(:,fluence),A(:,fluence),'LineWidth',linewidth)
    xlim([time(ni(fluence),fluence) time(nf(fluence),fluence)])
end

%% Plot Total Charge Fits
for fluence=1:n_fluences
    fig
    semilogx(time(1:n_time_points,fluence),Ctotal_temp(:,fluence)./(max(Ctotal_temp(:,fluence))),'LineWidth',linewidth)
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
