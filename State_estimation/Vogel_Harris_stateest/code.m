%Zachary Vogel
%Aubrey Harris
%For the State Estimation Class

%Quadrotor
clc, clear, close all

%Question 2
A=[0 1 0 0 0 0;
   0 -.0104 0 0 0 0;
   0 0 0 1 0 0;
   0 0 0 -.0104 0 0;
   0 0 0 0 0 1;
   0 0 0 0 0 -.0208];
B=[0 0 0 0;
   -.04167 0 .04167 0;
   0 0 0 0;
   0 -.04167 0 .04167;
   0 0 0 0;
   .4 .4 .4 .4];
C=[1 0 0 0 0 0;
   0 0 1 0 0 0;
   0 0 0 0 1 0];
D=[0 0 0 0;
   0 0 0 0;
   0 0 0 0];

T=[-.979884756356448, -4.82522265458613, -.230866541974081, -.612220975306226, .690563551134410, 1.22339914364605;
    6.89796334757239e-15, 1.48692340850896e-14, -.377537284162289, -.938512134420510, .714974800240852, 1.28716254599862;
    .979884756356460, 4.82522265458616, -.230866541974082, -.612220975306231, .690563551134411, 1.22339914364605;
    -2.50083851737725e-15, -1.30321138661878e-14, 1.23584087990799, 2.65069061583658, .446451060069997, .585765120120332];

%%
%Question 5


sysa=ss(A,B,C,D);
nat=abs(min(eig(A)));
freq=20*nat/(2*pi);

dT = 0.1; %sec
tvec = 0:dT:40;
DTsysa=c2d(sysa,dT,'zoh');
F = DTsysa.A;
G = DTsysa.B;
H = DTsysa.C;

Gamma   = [0 0 0 
          1 0 0;
          0 0 0;
          0 1 0;
          0 0 0;
          0 0 1];
       
Qtilde = [0.02 0.005 0.002 ;   %%My thoughts are that vertical winds should be weakest
         0.005 0.02  0.002;
         0.002 0.002 0.005]; 
 
M=dT*[-A, Gamma*Qtilde*Gamma';
     zeros(size(A)), A'];
E=expm(M);
Phi=E(7:12,7:12)';
Q=Phi*E(1:6,7:12);
%Q=diag([.001,.001,1,.001,.001,10]);
R=diag([.09325,.08853,.245]); %referenced scholarly article

%up forward left down back right
len=floor(length(tvec)/6);
u1= 5*[ones(1,len);ones(1,len);ones(1,len);ones(1,len)]; %straight up
u2= 5*[ones(1,len);0*ones(1,len);-1*ones(1,len);0*ones(1,len)];
u3= 5*[0*ones(1,len);1*ones(1,len);0*ones(1,len);-1*ones(1,len)];
u4=-u1;
u5=-u2;
u6=-u3;
len=length(tvec)-6*len;
u7=[zeros(1,len);zeros(1,len);zeros(1,len);zeros(1,len)];
u=[u1,u4,u2,u5,u3,u6,u7];

%%
%Question 6 Observability
%observability test is below

mu0 = [1 .5 1 0 1 0];

Hp1=[0 1 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 1];
Hp2=[0 0 1 0 0 0;1 0 0 0 0 0];
%actual observability test
obs=zeros(3*6,6);
obs1=zeros(3*6,6);
obs2=zeros(2*6,6);
for k=0:5
    obs((1+(k*3)):(3+(k*3)),:)=H*F^k;
    obs1((1+(k*3)):(3+(k*3)),:)=Hp1*F^k;
    obs2((1+(k*2)):(2+(k*2)),:)=Hp2*F^k;
end
rank(obs)
rank(obs1)
rank(obs2)

%%
%Question 7 Monte Carlo
mu0 = [1 0 1 0 1 0]; %%initial mean of x at k=0
P0 = eye(6); %%initial covar of x at k=0

muHistpred = zeros(6,length(tvec));
    muHistpred(:,1) = mu0;
PHistpred = zeros(6,6,length(tvec));
    PHistpred(:,:,1) = P0;

%%Draw samples to predict through state dynamics
nsamps = 100;
x_samps = mvnrnd(mu0',P0,nsamps)';    
x_sampHist = zeros(6,nsamps,length(tvec));
    x_sampHist(:,:,1) = x_samps;
    
%%Intialize plot for Quad    
figure(),
[Xell,Yell] = calc_gsigma_ellipse_plotpoints...
        ([mu0(1),mu0(3)]',[P0(1,1),P0(1,3);P0(3,1),P0(3,3);],2,100);
plot(Xell,Yell,'b.-'),hold on
plot(x_samps(1,:),x_samps(3,:),'r.')
xlabel('x position (m)')
ylabel('y position (m)')
title('Monte Carlo, Potential 2D Initial Position Distribution')
legend('Iterative Variance','Position')

   
for kk=1:length(tvec)-1
    %%Pure prediction update for A
    muHistpred(:,kk+1) = F*muHistpred(:,kk)+G*-T*muHistpred(:,kk); %ignore constant Guk;
    PHistpred(:,:,kk+1) = F*PHistpred(:,:,kk)*F'+Q;
    
    %%Sample update
    %%simulate process noise vector samples
    wkk_samp = mvnrnd([0,0,0,0,0,0]',Q,nsamps);
    x_sampHist(:,:,kk+1) = F*x_sampHist(:,:,kk)+wkk_samp';
    
    %%update plot for A
    if mod(kk,(length(tvec)-1)/8)==0
    [Xell,Yell] = calc_gsigma_ellipse_plotpoints...
    ([muHistpred(1,kk+1),muHistpred(3,kk+1)]',...
     [PHistpred(1,1,kk+1),PHistpred(1,3,kk+1);PHistpred(3,1,kk+1),PHistpred(3,3,kk+1)],...
      2,100);
   plot(Xell,Yell,'b.-')
    %%plot new ellipse   
    %%plot new Monte Carlo sample points
    plot(x_sampHist(1,:,kk+1),x_sampHist(3,:,kk+1),'r.')
    end
    
end

for k=0:2
    figure()
    subplot(211)
    plot(tvec(1:end),muHistpred(1+2*k,:),'b'), hold on
    plot(tvec(1:end),muHistpred(1+2*k,:)+ 2*sqrt(squeeze(PHistpred(1+2*k,1+2*k,:)))','b--')
    plot(tvec(1:end),muHistpred(1+2*k,:)- 2*sqrt(squeeze(PHistpred(1+2*k,1+2*k,:)))','b--')
    %xlabel('time (sec)')
     set(gca,'XTick',[]);
    if k==0
        title('Quadrotor Pure State Prediction')
        ylabel('x (m)')
    elseif k==1
        title('Quadrotor Pure State Prediction')
        ylabel('y (m)')
    else
        title('Quadrotor Pure State Prediction')
        ylabel('z (m)')
    end
    subplot(212)
    plot(tvec(1:end),muHistpred(2+2*k,:),'r'), hold on
    plot(tvec(1:end),muHistpred(2+2*k,:)+ 2*sqrt(squeeze(PHistpred(2+2*k,2+2*k,:)))','r--')
    plot(tvec(1:end),muHistpred(2+2*k,:)- 2*sqrt(squeeze(PHistpred(2+2*k,2+2*k,:)))','r--')
       xlabel('time (sec)')
    if k==0
        ylabel('dot x (m/s)')
    elseif k==1
        ylabel('dot y (m/s)')
    else
        ylabel('dot z (m/s)')
    end
end


%%
%Question 8
%%Initialize vanilla Kalman filter updates and...

%%
%Question 9
%NEES NIS Chisquare test


mk = mu0;
Pk = P0;
mk_filt_hist = zeros(6,length(tvec));
Pk_filt_hist = zeros(6,6,length(tvec));
Qkf = Q;


Nsimruns = 50;
NEESsamps = zeros(Nsimruns,length(tvec));
NISsamps = zeros(Nsimruns,length(tvec));


for ss=1:Nsimruns
    %%%1. %%Generate true trajectory and measurements from system
    
    xk_truehist = zeros(6,length(tvec));
    ykhist = zeros(3,length(tvec));
    xk = mvnrnd(mu0,P0); %sample initial robot state
    for jj=1:length(tvec)
        
        wk = mvnrnd(zeros(1,6),Q);
        xkp1 = F*xk' + G*(u(:,jj)-T*xk') + wk';
        vkp1 = mvnrnd(zeros(1,3),R)';
        ykp1 = H*xkp1 + vkp1;
        
        xk_truehist(:,jj) = xkp1;
        ykhist(:,jj) = ykp1;
        xk = xkp1';
    end
    xk_truehist5=xk_truehist;
    ykhist5=ykhist;
    %%% 2. Kalman Filter equations with simple NCV model
    
    %%Run the Kalman filter updates
    mk = mu0;
    Pk = P0;
    mk_filt_hist1 = zeros(6,length(tvec));
    Pk_filt_hist1 = zeros(6,6,length(tvec));
    innovk_hist = zeros(3,length(tvec));
    Pyyk_hist = zeros(3,3,length(tvec)); %%store measurement innovation covar
    Pyykp1_hist=zeros(3,3,length(tvec));
    NEESsshist = zeros(1,length(tvec));
    NISsshist = zeros(1,length(tvec));
    mkp1_plus=mu0';
    for jj=1:length(tvec)
        %I made a function to do the kalman filter, it's easier
        uk=u(:,jj)+-T*mkp1_plus;
        [K,mkp1_plus,Pkp1_plus,innov_kp1,Pyykp1]=KAL_FILT(F,G,H,uk,Pk,mk',Qkf,R,ykhist(:,jj));

        Pyykp1_hist(:,:,jj)=Pyykp1;
        mk = mkp1_plus';
        mk_filt_hist1(:,jj) = mkp1_plus;
        Pk = Pkp1_plus;
        Pk_filt_hist1(:,:,jj)= Pkp1_plus;
        innovk_hist(:,jj) = innov_kp1;
        
        %%Compute and store NEES and NIS statistics:
        invPkp1 = inv(Pkp1_plus);
        invPyykp1 = inv(Pyykp1);
        NEESsshist(jj) = ...
            (xk_truehist(:,jj) - mkp1_plus)'*invPkp1*(xk_truehist(:,jj) - mkp1_plus);
        NISsshist(jj) = innovk_hist(:,jj)'*invPyykp1*innovk_hist(:,jj);                
    end
    
    NEESsamps(ss,:) = NEESsshist;
    NISsamps(ss,:) = NISsshist;
    
end
 
figure(),

for k=1:6
  eval(['subplot(61',num2str(k),')'])
  plot(tvec,xk_truehist(k,:)-mk_filt_hist1(k,:),'r'), hold on
  plot(tvec,2*sqrt(squeeze(Pk_filt_hist1(k,k,:))'),'b--')
  plot(tvec,-2*sqrt(squeeze(Pk_filt_hist1(k,k,:))'),'b--')
  %xlabel('Time (sec)')
  if k<6
     set(gca,'XTick',[]);
  end
  if k==1
    title('State Estimation Errors vs Time')
    ylabel('X')
  elseif k==2
      ylabel('dot X')
  elseif k==3
      ylabel('Y')
  elseif k==4
      ylabel('dot Y')
  elseif k==5
      ylabel('Z')
  else
      ylabel('dot Z')
      xlabel('Time (sec)')
  end
end
  


error_Hist = xk_truehist-mk_filt_hist1; %Ex,k


figure()
plot(tvec,error_Hist)
title('E[exk]')
xlabel('time (s)')
ylabel('Truth Data Error')
legend('X','dot X','Y','dot Y', 'Z','dot Z')


%%
%Question 11 Batch least squares
nmeas=40; %number of measurements to use
sizing=size(ykhist(:,1));
nouts=sizing(1);%how many measurement devices
yki=zeros(nouts*nmeas,1);
sizing=size(F);
nstates=sizing(1);%number of states
Hki=zeros(nmeas*nouts,nstates);
sizing=size(u(:,1));
ninputs=sizing(1);
uki=zeros(nmeas*ninputs,1);
Gki=zeros(nmeas*nouts,nmeas*ninputs);
Cki=zeros(nmeas*nouts,nmeas*nstates);
Rki=zeros(nmeas*nouts,nmeas*nouts);
Qki=zeros(nmeas*nstates,nmeas*nstates);
for k=1:nmeas
   yki(1+nouts*(k-1):nouts+nouts*(k-1),1)=ykhist(:,k);
   Hki(1+nouts*(k-1):nouts+nouts*(k-1),:)=H*F^k;
   uki(1+ninputs*(k-1):ninputs+ninputs*(k-1))=u(:,1);
   Rki(1+nouts*(k-1):nouts+nouts*(k-1),1+nouts*(k-1):nouts+nouts*(k-1))=R;
   Qki(1+nstates*(k-1):nstates+nstates*(k-1),1+nstates*(k-1):nstates+nstates*(k-1))=Q;
   
   for kk=k:nmeas-1
      Gki(1+(kk-1)*nouts:nouts+(kk-1)*nouts,1+(k-1)*ninputs:ninputs+(k-1)*ninputs)=H*F^(kk-k)*G;
      Cki(1+(kk-1)*nouts:nouts+(kk-1)*nouts,1+(k-1)*nstates:nstates+(k-1)*nstates)=H*F^(kk-k);
   end
end
zetaki=yki-Gki*uki;
Lki=[Hki,Cki];


BLS1=(Hki'*Rki^(-1)*Hki)^(-1)*Hki'*Rki^(-1);
BLS2=(Cki'*Rki^(-1)*Cki+Qki^(-1))^(-1)*Cki'*Rki^(-1);
Abls=BLS1*Cki;
Bbls=BLS1*(yki-Gki*uki);
Cbls=BLS2*Hki;
Dbls=BLS2*(yki-Gki*uki);

x0bls=(eye(size(Abls*Cbls))-Abls*Cbls)^(-1)*(Abls*Dbls+Bbls);
wbls=(eye(size(Cbls*Abls))-Cbls*Abls)^(-1)*(Cbls*Bbls+Dbls);

%%
%Sim 1
%%DO NEES TEST:


epsNEESbar = mean(NEESsamps,1);
alphaNEES = 0.05; %%significance level
Nnx = Nsimruns*length(F);
%%compute intervals:
r1x = chi2inv(alphaNEES/2, Nnx )./ Nsimruns;
r2x = chi2inv(1-alphaNEES/2, Nnx )./ Nsimruns;
figure()
plot(epsNEESbar,'ro','MarkerSize',6,'LineWidth',2),hold on
plot(r1x*ones(size(epsNEESbar)),'r--','LineWidth',2)
plot(r2x*ones(size(epsNEESbar)),'r--','LineWidth',2)
ylabel('NEES statistic, bar epsilon x','FontSize',14)
xlabel('time step, k','FontSize',14)
title('50 run, .05 NEES Estimation Results','FontSize',14)
legend('NEES @ time k', 'r_1 bound', 'r_2 bound')


%%DO NIS TEST:
epsNISbar = mean(NISsamps,1);
alphaNIS = 0.05;
Nny = Nsimruns*size(H,1);
%%compute intervals:
r1y = chi2inv(alphaNIS/2, Nny )./ Nsimruns;
r2y = chi2inv(1-alphaNIS/2, Nny )./ Nsimruns;
figure()
plot(epsNISbar,'bo','MarkerSize',6,'LineWidth',2),hold on
plot(r1y*ones(size(epsNISbar)),'b--','LineWidth',2)
plot(r2y*ones(size(epsNISbar)),'b--','LineWidth',2)
ylabel('NIS statistic, bar epsilon y','FontSize',14)
xlabel('time step, k','FontSize',14)
title('50 run, .05 NIS Estimation Results','FontSize',14)
legend('NIS @ time k', 'r_1 bound', 'r_2 bound')
%%
%Sim #2


Nsimruns = 100;
NEESsamps = zeros(Nsimruns,length(tvec));
NISsamps = zeros(Nsimruns,length(tvec));


for ss=1:Nsimruns
    %%%1. %%Generate true trajectory and measurements from system
    
    xk_truehist = zeros(6,length(tvec));
    ykhist = zeros(3,length(tvec));
    xk = mvnrnd(mu0,P0); %sample initial robot state
    for jj=1:length(tvec)
        
        wk = mvnrnd(zeros(1,6),Q);
        xkp1 = F*xk' + G*(u(:,jj)-T*xk') + wk';
        vkp1 = mvnrnd(zeros(1,3),R)';
        ykp1 = H*xkp1 + vkp1;
        
        xk_truehist(:,jj) = xkp1;
        ykhist(:,jj) = ykp1;
        xk = xkp1';
    end
    
    %%% 2. Kalman Filter equations with simple NCV model
    
    %%Run the Kalman filter updates
    mk = mu0;
    Pk = P0;
    mk_filt_hist = zeros(6,length(tvec));
    Pk_filt_hist = zeros(6,6,length(tvec));
    innovk_hist = zeros(3,length(tvec));
    Pyyk_hist = zeros(3,3,length(tvec)); %%store measurement innovation covar
    Pyykp1_hist=zeros(3,3,length(tvec));
    NEESsshist = zeros(1,length(tvec));
    NISsshist = zeros(1,length(tvec));
    mkp1_plus=mu0';
    
    for jj=1:length(tvec)
        
        uk=u(:,jj)+-T*mkp1_plus;
        
        [K,mkp1_plus,Pkp1_plus,innov_kp1,Pyykp1]=KAL_FILT(F,G,H,uk,Pk,mk',Qkf,R,ykhist(:,jj));

        Pyykp1_hist(:,:,jj)=Pyykp1;
        mk = mkp1_plus';
        mk_filt_hist(:,jj) = mkp1_plus;
        Pk = Pkp1_plus;
        Pk_filt_hist(:,:,jj)= Pkp1_plus;
        innovk_hist(:,jj) = innov_kp1;
        
        %%Compute and store NEES and NIS statistics:
        invPkp1 = inv(Pkp1_plus);
        invPyykp1 = inv(Pyykp1);
        NEESsshist(jj) = ...
            (xk_truehist(:,jj) - mkp1_plus)'*invPkp1*(xk_truehist(:,jj) - mkp1_plus);
        NISsshist(jj) = innovk_hist(:,jj)'*invPyykp1*innovk_hist(:,jj);                
    end
    
    NEESsamps(ss,:) = NEESsshist;
    NISsamps(ss,:) = NISsshist;
    
end


%%DO NEES TEST:


epsNEESbar = mean(NEESsamps,1);
alphaNEES = 0.01; %%significance level
Nnx = Nsimruns*length(F);
%%compute intervals:
r1x = chi2inv(alphaNEES/2, Nnx )./ Nsimruns;
r2x = chi2inv(1-alphaNEES/2, Nnx )./ Nsimruns;
figure()
plot(epsNEESbar,'ro','MarkerSize',6,'LineWidth',2),hold on
plot(r1x*ones(size(epsNEESbar)),'r--','LineWidth',2)
plot(r2x*ones(size(epsNEESbar)),'r--','LineWidth',2)
ylabel('NEES statistic, bar epsilon x','FontSize',14)
xlabel('time step, k','FontSize',14)
title('100 run, .01 NEES Estimation Results','FontSize',14)
legend('NEES @ time k', 'r_1 bound', 'r_2 bound')


%%DO NIS TEST:
epsNISbar = mean(NISsamps,1);
alphaNIS = 0.01;
Nny = Nsimruns*size(H,1);
%%compute intervals:
r1y = chi2inv(alphaNIS/2, Nny )./ Nsimruns;
r2y = chi2inv(1-alphaNIS/2, Nny )./ Nsimruns;
figure()
plot(epsNISbar,'bo','MarkerSize',6,'LineWidth',2),hold on
plot(r1y*ones(size(epsNISbar)),'b--','LineWidth',2)
plot(r2y*ones(size(epsNISbar)),'b--','LineWidth',2)
ylabel('NIS statistic, bar epsilon y','FontSize',14)
xlabel('time step, k','FontSize',14)
title('100 run, .01 NIS Estimation Results','FontSize',14)
legend('NIS @ time k', 'r_1 bound', 'r_2 bound')


%Sim #3


Nsimruns = 200;
NEESsamps = zeros(Nsimruns,length(tvec));
NISsamps = zeros(Nsimruns,length(tvec));


for ss=1:Nsimruns
    %%%1. %%Generate true trajectory and measurements from system
    
    xk_truehist = zeros(6,length(tvec));
    ykhist = zeros(3,length(tvec));
    xk = mvnrnd(mu0,P0); %sample initial robot state
    for jj=1:length(tvec)
        
        wk = mvnrnd(zeros(1,6),Q);
        xkp1 = F*xk' + G*(u(:,jj)-T*xk') + wk';
        vkp1 = mvnrnd(zeros(1,3),R)';
        ykp1 = H*xkp1 + vkp1;
        
        xk_truehist(:,jj) = xkp1;
        ykhist(:,jj) = ykp1;
        xk = xkp1';
    end
    
    %%% 2. Kalman Filter equations with simple NCV model
    
    %%Run the Kalman filter updates
    mk = mu0;
    Pk = P0;
    mk_filt_hist = zeros(6,length(tvec));
    Pk_filt_hist = zeros(6,6,length(tvec));
    innovk_hist = zeros(3,length(tvec));
    Pyyk_hist = zeros(3,3,length(tvec)); %%store measurement innovation covar
    Pyykp1_hist=zeros(3,3,length(tvec));
    NEESsshist = zeros(1,length(tvec));
    NISsshist = zeros(1,length(tvec));
    mkp1_plus=mu0';
    
    for jj=1:length(tvec)
        
        uk=-T*mkp1_plus+u(:,jj);
        [K,mkp1_plus,Pkp1_plus,innov_kp1,Pyykp1]=KAL_FILT(F,G,H,uk,Pk,mk',Qkf,R,ykhist(:,jj));

        Pyykp1_hist(:,:,jj)=Pyykp1;
        mk = mkp1_plus';
        mk_filt_hist(:,jj) = mkp1_plus;
        Pk = Pkp1_plus;
        Pk_filt_hist(:,:,jj)= Pkp1_plus;
        innovk_hist(:,jj) = innov_kp1;
        
        %%Compute and store NEES and NIS statistics:
        invPkp1 = inv(Pkp1_plus);
        invPyykp1 = inv(Pyykp1);
        NEESsshist(jj) = ...
            (xk_truehist(:,jj) - mkp1_plus)'*invPkp1*(xk_truehist(:,jj) - mkp1_plus);
        NISsshist(jj) = innovk_hist(:,jj)'*invPyykp1*innovk_hist(:,jj);                
    end
    
    NEESsamps(ss,:) = NEESsshist;
    NISsamps(ss,:) = NISsshist;
    
end

% figure
% plot3(mk_filt_hist(1,:),mk_filt_hist(3,:),mk_filt_hist(5,:))
% title('estimated trajectory')
% xlabel('x')
% ylabel('y')
% zlabel('z')
%%DO NEES TEST:



epsNEESbar = mean(NEESsamps,1);
alphaNEES = 0.003; %%significance level
Nnx = Nsimruns*length(F);
%%compute intervals:
r1x = chi2inv(alphaNEES/2, Nnx )./ Nsimruns;
r2x = chi2inv(1-alphaNEES/2, Nnx )./ Nsimruns;
figure()
plot(epsNEESbar,'ro','MarkerSize',6,'LineWidth',2),hold on
plot(r1x*ones(size(epsNEESbar)),'r--','LineWidth',2)
plot(r2x*ones(size(epsNEESbar)),'r--','LineWidth',2)
ylabel('NEES statistic, bar epsilon x','FontSize',14)
xlabel('time step, k','FontSize',14)
title('200 run, .003 NEES Estimation Results','FontSize',14)
legend('NEES @ time k', 'r_1 bound', 'r_2 bound')


%%DO NIS TEST:
epsNISbar = mean(NISsamps,1);
alphaNIS = 0.003;
Nny = Nsimruns*size(H,1);
%%compute intervals:
r1y = chi2inv(alphaNIS/2, Nny )./ Nsimruns;
r2y = chi2inv(1-alphaNIS/2, Nny )./ Nsimruns;
figure()
plot(epsNISbar,'bo','MarkerSize',6,'LineWidth',2),hold on
plot(r1y*ones(size(epsNISbar)),'b--','LineWidth',2)
plot(r2y*ones(size(epsNISbar)),'b--','LineWidth',2)
ylabel('NIS statistic, bar epsilon y','FontSize',14)
xlabel('time step, k','FontSize',14)
title('200 run, .003 NIS Estimation Results','FontSize',14)
legend('NIS @ time k', 'r_1 bound', 'r_2 bound')

%%
%Question 10 Qkf not equal to Q


Nsimruns = 50;
NEESsamps = zeros(Nsimruns,length(tvec));
NISsamps = zeros(Nsimruns,length(tvec));
Qkf=5000*Q;
for ss=1:Nsimruns
    %%%1. %%Generate true trajectory and measurements from system
    
    xk_truehist = zeros(6,length(tvec));
    ykhist = zeros(3,length(tvec));
    xk = mvnrnd(mu0,P0); %sample initial robot state
    for jj=1:length(tvec)
        
        wk = mvnrnd(zeros(1,6),Qkf);
        xkp1 = F*xk' + G*(u(:,jj)-T*xk') + wk';
        vkp1 = mvnrnd(zeros(1,3),R)';
        ykp1 = H*xkp1 + vkp1;
        
        xk_truehist(:,jj) = xkp1;
        ykhist(:,jj) = ykp1;
        xk = xkp1';
    end
    
    %%% 2. Kalman Filter equations with simple NCV model
    
    %%Run the Kalman filter updates
    mk = mu0;
    Pk = P0;
    mk_filt_hist = zeros(6,length(tvec));
    Pk_filt_hist = zeros(6,6,length(tvec));
    innovk_hist = zeros(3,length(tvec));
    Pyyk_hist = zeros(3,3,length(tvec)); %%store measurement innovation covar
    Pyykp1_hist=zeros(3,3,length(tvec));
    NEESsshist = zeros(1,length(tvec));
    NISsshist = zeros(1,length(tvec));
    mkp1_plus=mu0';  
    
    for jj=1:length(tvec)
        uk=u(:,jj)-T*mkp1_plus;
        [K,mkp1_plus,Pkp1_plus,innov_kp1,Pyykp1]=KAL_FILT(F,G,H,uk,Pk,mk',Qkf,R,ykhist(:,jj));

        Pyykp1_hist(:,:,jj)=Pyykp1;
        mk = mkp1_plus';
        mk_filt_hist(:,jj) = mkp1_plus;
        Pk = Pkp1_plus;
        Pk_filt_hist(:,:,jj)= Pkp1_plus;
        innovk_hist(:,jj) = innov_kp1;
        
        %%Compute and store NEES and NIS statistics:
        invPkp1 = inv(Pkp1_plus);
        invPyykp1 = inv(Pyykp1);
        NEESsshist(jj) = ...
            (xk_truehist(:,jj) - mkp1_plus)'*invPkp1*(xk_truehist(:,jj) - mkp1_plus);
        NISsshist(jj) = innovk_hist(:,jj)'*invPyykp1*innovk_hist(:,jj);                
    end
    
    NEESsamps(ss,:) = NEESsshist;
    NISsamps(ss,:) = NISsshist;
    
end


%%DO NEES TEST:


epsNEESbar = mean(NEESsamps,1);
alphaNEES = 0.003; %%significance level
Nnx = Nsimruns*length(F);
%%compute intervals:
r1x = chi2inv(alphaNEES/2, Nnx )./ Nsimruns;
r2x = chi2inv(1-alphaNEES/2, Nnx )./ Nsimruns;
figure()
plot(epsNEESbar,'ro','MarkerSize',6,'LineWidth',2),hold on
plot(r1x*ones(size(epsNEESbar)),'r--','LineWidth',2)
plot(r2x*ones(size(epsNEESbar)),'r--','LineWidth',2)
ylabel('NEES statistic, bar epsilon x','FontSize',14)
xlabel('time step, k','FontSize',14)
title('5000*Q .003 NEES Estimation Results','FontSize',14)
legend('NEES @ time k', 'r_1 bound', 'r_2 bound')


%%DO NIS TEST:
epsNISbar = mean(NISsamps,1);
alphaNIS = 0.003;
Nny = Nsimruns*size(H,1);
%%compute intervals:
r1y = chi2inv(alphaNIS/2, Nny )./ Nsimruns;
r2y = chi2inv(1-alphaNIS/2, Nny )./ Nsimruns;
figure()
plot(epsNISbar,'bo','MarkerSize',6,'LineWidth',2),hold on
plot(r1y*ones(size(epsNISbar)),'b--','LineWidth',2)
plot(r2y*ones(size(epsNISbar)),'b--','LineWidth',2)
ylabel('NIS statistic, bar epsilon y','FontSize',14)
xlabel('time step, k','FontSize',14)
title('5000*Q .003 NIS Estimation Results','FontSize',14)
legend('NIS @ time k', 'r_1 bound', 'r_2 bound')




%%
%Question 13 Information Filter


mu_info_hist=zeros(6,length(tvec));
p_info_hist=zeros(6,6,length(tvec));
Ikp=ones([6,6]);
ikp=[1,1,1,1,1,1]';

mkp1_plus=mu0';

for c=1:length(tvec)

uc=u(:,c)-T*mkp1_plus;

Zk=inv(F)'*Ikp*inv(F);
%Prediction
ikp1_minus=(eye(6)-Zk*inv(Zk*inv(Q)))*(inv(F)'*ikp+Zk*G*uc); %mean 6x1
Ikp1_minus=Zk-Zk*inv(Zk*inv(Q))*Zk; %covar 6x6
%Measurement
ikp1p=ikp1_minus+(H'*inv(R)*ykhist(:,c));
Ikp1p=Ikp1_minus+(H'*inv(R)*H);
%Recovery
mkp1_plus=inv(Ikp1p)*ikp1p;
%mu_info_hist(:,c)=mkp1_plus;
pkp1_plus=inv(Ikp1p);
%p_info_hist(:,:,c)=pkp1_plus;


ikp=ikp1p;
Ikp=Ikp1p;
end
%%

Ykp=inv(P0);
ykp=inv(P0)*mu0';
xk = mvnrnd(mu0,P0)';

mkp1_plus=mu0';

for f=1:length(tvec)

uc=u(:,f)-T*mkp1_plus;
    
wk1 = mvnrnd(zeros(1,6),Q)';
xkp = F*xk + G*uc + wk1;
vkp = mvnrnd(zeros(1,3),R)';
zk = H*xkp + vkp;

%predict
Mk=inv(F)'*Ykp*inv(F);
Ck=Mk*inv(Mk+inv(Q));
Lk=eye(6)-Ck;
Ykminus=Lk*Mk*Lk'+Ck*inv(Q)*Ck';
ykminus=Lk*inv(F)'*ykp;
%measurement
Ik=H'*inv(R)*H;
ik=H'*inv(R)*zk;
%update
Ykp1=Ykminus+Ik;
ykp1=ykminus+ik;
%Recovery
mkp1_plus=inv(Ykp1)*ykp1;
mu_info_hist(:,f)=mkp1_plus;
pkp1_plus=inv(Ykp1);
p_info_hist(:,:,f)=pkp1_plus;

xk=xkp1;
Yk=Ykp1;
yk=ykp1;

end

figure(),
for k=1:6
  eval(['subplot(61',num2str(k),')'])
  plot(tvec,xk_truehist(k,:)-mu_info_hist(k,:),'r'), hold on
  plot(tvec,2*sqrt(squeeze(p_info_hist(k,k,:))'),'b--')
  plot(tvec,-2*sqrt(squeeze(p_info_hist(k,k,:))'),'b--')
  %xlabel('Time (sec)')
  if k<6
     set(gca,'XTick',[]);
  end
  if k==1
    title('Information Filter Errors vs Time')
    ylabel('X')
  elseif k==2
      ylabel('dot X')
  elseif k==3
      ylabel('Y')
  elseif k==4
      ylabel('dot Y')
  elseif k==5
      ylabel('Z')
  else
      ylabel('dot Z')
      xlabel('Time (sec)')
  end
end


%%
%I'll probably do 16 cause it's interesting.
a1=F';
a2=F'*H';
scale=[0.9 0.5 0.2];
a=[-1/(sqrt(2))-1/(sqrt(2))*i (-1/(sqrt(2))+i/(sqrt(2))) -1/(2*sqrt(2))+i/(2*sqrt(2)) (-1/(2*sqrt(2))-i/(2*sqrt(2))) 1 0.5];
a31=scale(1)*a;
a32=scale(2)*a;
a33=scale(3)*a;
Kplace1=place(a1,a2,a31)';
Kplace2=place(a1,a2,a32)';
Kplace3=place(a1,a2,a33)';

xk_truehist=xk_truehist5;
ykhist=ykhist5;

xk1 = mu0';
Pk1 = P0;
xk2=xk1;
xk3=xk1;
Pk2=Pk1;
Pk3=Pk1;
xhat_hist1 = zeros(6,length(tvec));
Plo_hist1 = zeros(6,6,length(tvec));
xhat_hist2=xhat_hist1;
Plo_hist2=Plo_hist1;
xhat_hist3=xhat_hist1;
Plo_hist3=Plo_hist1;
xkplus1=zeros(6,1);
xkplus2=zeros(6,1);
xkplus3=zeros(6,1);

for k=1:length(tvec)
    uk=u(:,k)-T*xkplus1;
    [xkplus1,Pkplus1,yerror]=Burger_Observer(F,G,H,uk,Pk1,xk1,Kplace1,Q,R,ykhist(:,k));
    xk1 = xkplus1;
    xhat_hist1(:,k) = xkplus1;
    Pk1 = Pkplus1;
    Plo_hist1(:,:,k)= Pkplus1;
    uk=u(:,k)-T*xkplus2;
    [xkplus2,Pkplus2,yerror]=Burger_Observer(F,G,H,uk,Pk2,xk2,Kplace2,Q,R,ykhist(:,k));
    xk2 = xkplus2;
    xhat_hist2(:,k) = xkplus2;
    Pk2 = Pkplus2;
    Plo_hist2(:,:,k)= Pkplus2;
    uk=u(:,k)-T*xkplus3;
    [xkplus3,Pkplus3,yerror]=Burger_Observer(F,G,H,uk,Pk3,xk3,Kplace3,Q,R,ykhist(:,k));
    xk3 = xkplus3;
    xhat_hist3(:,k) = xkplus3;
    Pk3 = Pkplus3;
    Plo_hist3(:,:,k)= Pkplus3;
end

for ll=1:3
    figure(),
    if ll==1
        xkhat_hist=xhat_hist1;
        Plo_hist=Plo_hist1;
    elseif ll==2
        xkhat_hist=xhat_hist2;
        Plo_hist=Plo_hist2;
    else
        xkhat_hist=xhat_hist3;
        Plo_hist=Plo_hist3;
    end
    for k=1:6
        eval(['subplot(61',num2str(k),')'])
        plot(tvec,xk_truehist(k,:)-xkhat_hist(k,:),'r'), hold on
        plot(tvec,2*sqrt(squeeze(Plo_hist(k,k,:))'),'b--')
        plot(tvec,-2*sqrt(squeeze(Plo_hist(k,k,:))'),'b--')
        %xlabel('Time (sec)')
        if k<6
            set(gca,'XTick',[]);
        end
        if k==1
            if ll==1
                title('State Estimation Errors vs Time Pole Placement 1')
            elseif ll==2
                title('State Estimation Errors vs Time Pole Placement 2')
            else
                title('State Estimation Errors vs Time Pole Placement 3')
            end 
            ylabel('X')
        elseif k==2
            ylabel('dot X')
        elseif k==3
            ylabel('Y')
        elseif k==4
            ylabel('dot Y')
        elseif k==5
            ylabel('Z')
        else
            ylabel('dot Z')
            xlabel('Time (sec)')
        end
    end

end

%%
%  %Question 14 Square Root Filter

xsr_truehist = zeros(6,length(tvec));
     ysrhist = zeros(3,length(tvec));
xsr = mvnrnd(mu0,P0);
Qsr=sqrt(Q);
Rsr=chol(R,'lower');

for jj=1:length(tvec)
    wsr = mvnrnd(zeros(1,6),Qsr);
         xkp1 = F*xsr' + G*u(:,jj) + wsr';
         vkp1 = mvnrnd(zeros(1,3),Rsr)';
         ykp1 = H*xkp1 + vkp1;
         
         xsr_truehist(:,jj) = xkp1;
         ysrhist(:,jj) = ykp1;
         xk = xkp1';
end

%Run the Square Root Kalman filter updates
     msr = sqrt(mu0);
     Psr = P0;
     msr_filt_hist = zeros(6,length(tvec));
     Psr_filt_hist = zeros(6,6,length(tvec));
     innovk_hist = zeros(3,length(tvec));
     S=chol(Psr,'lower');
     Psr=S*S';
     for jj=1:length(tvec)
         
         msr=sqrt(msr);
         %%Perform prediction step
         mkp1_minus = F*msr' + G*u(:,jj);
         Pkp1_minus = F*Psr*F' + Qsr;
         
         %%Compute Kalman gain
         Kkp1 = Pkp1_minus*H'/(H*Pkp1_minus*H' + Rsr);
         %%Perform measurement update step
         ykp1_report = ysrhist(:,jj); %simulate the reporting of data from sensor
         ykp1_pred = H*mkp1_minus; %predicted measurement
         innov_kp1 = ykp1_report - ykp1_pred; %compute meas innovation
         mkp1_plus = mkp1_minus + Kkp1*innov_kp1; %compute update to state mean
         Pkp1_plus = (eye(6) - Kkp1*H)*Pkp1_minus; %compute update to covar
         
         msr = mkp1_plus';
         msr_filt_hist(:,jj) = mkp1_plus;
         Psr = Pkp1_plus;
         Psr_filt_hist(:,:,jj)= Pkp1_plus;
         innovk_hist(:,jj) = innov_kp1;
         
     end
     
     figure(),
   subplot(611)
   plot(tvec,xsr_truehist(1,:)-msr_filt_hist(1,:),'r'), hold on
   plot(tvec,2*sqrt(squeeze(Psr_filt_hist(1,1,:))'),'b--')
   plot(tvec,-2*sqrt(squeeze(Psr_filt_hist(1,1,:))'),'b--')
     set(gca,'XTick',[]);
   %xlabel('Time (sec)')
   ylabel('X')
   title('Square Root Filter Errors vs Time'),
   subplot(612)
   plot(tvec,xsr_truehist(2,:)-msr_filt_hist(2,:),'r'), hold on
   plot(tvec,2*sqrt(squeeze(Psr_filt_hist(2,2,:))'),'b--')
   plot(tvec,-2*sqrt(squeeze(Psr_filt_hist(2,2,:))'),'b--')
   %xlabel('Time (sec)')
     set(gca,'XTick',[]);
   ylabel('dot X')
   subplot(613)
   plot(tvec,xsr_truehist(3,:)-msr_filt_hist(3,:),'r'), hold on
   plot(tvec,2*sqrt(squeeze(Psr_filt_hist(3,3,:))'),'b--')
   plot(tvec,-2*sqrt(squeeze(Psr_filt_hist(3,3,:))'),'b--')
   %xlabel('Time (sec)')
     set(gca,'XTick',[]);
   ylabel('Y')
   subplot(614)
   plot(tvec,xsr_truehist(4,:)-msr_filt_hist(4,:),'r'), hold on
   plot(tvec,2*sqrt(squeeze(Psr_filt_hist(4,4,:))'),'b--')
   plot(tvec,-2*sqrt(squeeze(Psr_filt_hist(4,4,:))'),'b--')
   %xlabel('Time (sec)')
     set(gca,'XTick',[]);
   ylabel('dot Y')
   subplot(615)
   plot(tvec,xsr_truehist(5,:)-msr_filt_hist(5,:),'r'), hold on
   plot(tvec,2*sqrt(squeeze(Psr_filt_hist(5,5,:))'),'b--')
   plot(tvec,-2*sqrt(squeeze(Psr_filt_hist(5,5,:))'),'b--')
   %xlabel('Time (sec)')
   ylabel('Z')
     set(gca,'XTick',[]);
   subplot(616)
   plot(tvec,xsr_truehist(6,:)-msr_filt_hist(6,:),'r'), hold on
   plot(tvec,2*sqrt(squeeze(Psr_filt_hist(6,6,:))'),'b--')
   plot(tvec,-2*sqrt(squeeze(Psr_filt_hist(6,6,:))'),'b--')
   xlabel('Time (sec)')
   ylabel('dot Z')
   
   
   %%
   %save all the figures
   str={'MC','MCx','MCy','MCz','KF','errors','50NEES','50NIS','100NEES','100NIS','200NEES','200NIS','5000NEES','5000NIS','IF','pole1','pole2','pole3','SRKF3'};
   
   h = get(0,'children');
   hold on
    for i=1:length(h)
        saveas(h(i), str{length(str)+1-i},'png');
    end