clear all
time_raw=importdata('Time_data_kinetics for Paul.xlsx');
A       =importdata('Y6_combined.xlsx');
k       =importdata('bestk.mat');
fluences=A(:,1);  I=A(:,2)/max(A(:,2));
datarange=1:16;
datarange_sim = 1:30;
fluences_sim = 1e13*logspace(1,10,50);
time=[0:(10^-12):1*10^-8]; n_time_points=length(time);
linewidth=1;

%% Integrate Rate Equations
na=1; nb=1; nd=1; ne=1;
pb = CmdLineProgressBar('Iterations:');
besterr=inf;
S_temp=zeros(n_time_points,1);    S =zeros(n_time_points,1);    
CT_temp=zeros(n_time_points,1);   CT=zeros(n_time_points,1);
C_temp=zeros(n_time_points,1);    C =zeros(n_time_points,1);
fig
semilogx(fluences,I*2,'o','MarkerEdgeColor',[0 0 0],'MarkerSize',6)
hold on
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
                k=[kradnr kct kcr knrct kb kcs kenc knrc];
                opts   = odeset('RelTol',1e-3);
                for fluence=datarange_sim
                    i0 =fluences_sim(fluence);
                    [t,c]=ode15s(@(t,c) diffun_chargesCT_600ps(t,c,k,i0),time,c0,opts);
                    PLQEraw(fluence)=trapz(t,krad*c(:,1))/i0;
                    S_temp(:,fluence) =c(:,1);
                    CT_temp(:,fluence)=c(:,2);
                    C_temp(:,fluence) =c(:,3);
                    Ctotal_temp(:,fluence)=C_temp(:,fluence)+CT_temp(:,fluence);
                    Ctotal_temp(:,fluence)=Ctotal_temp(:,fluence);
                end
                %err=norm(PLQEraw./max(PLQEraw)-I(datarange,:)');
%                 if  err<besterr
%                     besterr=err;
%                     bestk=k; bestc=c;
%                     S=S_temp; CT=CT_temp; C=C_temp; Ctotal=Ctotal_temp;
%                     PLQE=PLQEraw;
                %end
            end
%             pb.print(d,nd)  
        end
        semilogx(fluences_sim(datarange_sim)',PLQEraw*100,'LineWidth',1)
    end
    pb.print(a,na)
end
legend
% plot(fluences(datarange)',PLQE./max(PLQE),'LineWidth',1)

%% Plot Exciton/CT/Charge Fits
for fluence=1:13
    fig
    semilogx(time,S(:,fluence),'LineWidth',linewidth)
    hold on
    semilogx(time,CT(:,fluence),'LineWidth',linewidth)
    semilogx(time,C(:,fluence),'LineWidth',linewidth)
end
