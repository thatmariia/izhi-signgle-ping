clear all; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change the input strength and see how the spike raster chanhe and
% the frequency of the network oscillation (see TFR main power band)
% you can try input_net1 =  4/6/10/15/20
input_net1 = 5; % excitatory input to network 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of excitatory and inhibitory neurons
Ne=200;     Ni=50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=[0.02*ones(Ne,1);        0.1*ones(Ni,1)]; % a =tiemscale of recover varibale u
b=[0.2*ones(Ne,1);         0.2*ones(Ni,1)]; % b= sensitivity of  u to subthreshold oscillations
c=[-65*ones(Ne,1);        -65*ones(Ni,1)]; % c= membrane voltage after spike (reset)
d=[8*ones(Ne,1);           2*ones(Ni,1)]; % d= spike reset of recover varibale u
v=-65*ones(Ne+Ni,1);    % Initial values of v = voltage
u=b.*v;                 % Initial values of u= membrane recovery variable
firings=[];             % spike timings
simulation_time=3000 ;
dt=1;
Ntot=(Ne+Ni);
%%%%%%%%%%%%%% input strength  gaussian input%%%%%%%%%%%%
% to excitatory neurons  %%%%
mean_E= [ input_net1*ones(Ne,1) ];
var_E= 2;
% to inhibitory neurons
mean_I= [ 4*ones(Ni,1) ];
var_I= 1.5;

%%%%%%%%%%%%%%%%%%%%
gampa=zeros(Ne,1);
gaba=zeros(Ni,1);
decay_ampa =1;decay_gaba =7;
rise_ampa =0.15;rise_gaba =0.2;

%% Constructing conenctivity matrix %%%%%%%
EE = 0.05 ;% ecitatory to excitatory
EI = 0.4;% ecitatory to inhibitory
IE = 0.3 ;% inhbitory to excitatory
II = 0.2; % inhibitory to inhibitory


S=zeros(Ntot);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% WITHIN NETWORK CONNECTIVITY %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E1_ind=1:Ne;
I1_ind=Ne+1:Ntot;
%%%  E - E  %%%%%%%%%
S(E1_ind,E1_ind) =   EE*rand(Ne);
%%%  E - I   %%%%%%%%%
S(I1_ind,E1_ind) =  EI*rand(Ni,Ne);
%%%  I - E  %%%%%%%%%
S(  E1_ind ,I1_ind) =     -IE*rand(Ne,Ni);
%%%  I - I %%%%%%%%%
S( I1_ind,I1_ind ) =  -II*rand(Ni);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:dt:simulation_time          % simulation
    I=[var_E*randn(Ne,1)+mean_E;var_I*randn(Ni,1)+mean_I]; % thalamic input
    fired=find(v>=30);    % indices of spikes
    firings=[firings; t+0*fired,fired];
    v(fired)=c(fired);
    u(fired)=u(fired)+d(fired);
    gampa=gampa + dt* (0.3*(((1+tanh((v(1:Ne)/10)  +2 ))/2).*(1-gampa)/rise_ampa - gampa/decay_ampa));
    gaba=gaba + dt* (0.3*(((1+tanh((v(Ne+1:end)/10)  +2 ))/2).*(1-gaba)/rise_gaba - gaba/decay_gaba));
    
    I=I+S*[gampa ;gaba]; % integrate input from other neurons
    
    
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
    u=u+a.*(b.*v-u);                 % stability
    
end
toc

%%% spike probability for given time bin
excpop1=find(firings(:,2) >0 & firings(:,2) <= Ne);
inhpop1=find(firings(:,2) >Ne & firings(:,2) <= Ntot);
% population signal
[t1,t2] = hist(firings(excpop1,1),0:1:t);
signal1=t1(300:end);

figure('Color','w','Position' ,[ 200 200  950 450],'name',[' network input: '  num2str(input_net1)]),
subplot(2,1,1,'Fontsize',17) %
firing1=firings(find(firings(:,2) <=Ne),:);
firing2=firings(find(firings(:,2) > Ne &  firings(:,2) < Ntot  ),:);
plot(firing1(:,1),firing1(:,2),'.','Color', [ 0.8 0.2 0.2]);
hold on, plot(firing2(:,1),firing2(:,2),'.','Color', [ 0.2 0.2 0.8]);
axis tight
title('spike raster');xlim([ 700 1500]), ylabel('neuron N')

subplot(2,1,2,'Fontsize',17) %
t=(1:simulation_time)./0.001; t=t(299:end);
plot(t,signal1);
xlabel('Time (ms) ');axis tight
%%%%%%%%%% make a TFR using what you have learned %%%%%%%%%%%%
%%






































































%%
% If you get stuck!
% 
% Wt=-1:0.001:1;
% 
% half_wave=floor((length(Wt)-1)/2); %% half the size of the wavelet
% nConv = length(Wt)+length(signal1)-1; %% the size of the data + zero padding
% FFT_data=fft(signal1,nConv); 
% 
% freq_of_interest=[20:1:60];
% TFR=NaN(length(freq_of_interest),length(t));
% 
% SD=[.05]; %% set the width of the Gaussian
% for freq=freq_of_interest
%     G=exp(-(Wt).^2./(2*SD^2));
%     complexSine = exp(1i*2*pi*freq*Wt );
%     complex_wavelet=complexSine.*G;
%     FFT_wavelet=fft(complex_wavelet,nConv); FFT_wavelet=FFT_wavelet./max(FFT_wavelet);
%     
%     tmp = ifft(FFT_wavelet.*FFT_data,nConv);
%     TFR(find(freq_of_interest==freq),:)  = tmp((half_wave)+1:end-half_wave); %%% trim the edges, these are the bits we included by zero padding
% end
% 
% imagesc(t,freq_of_interest,abs(TFR)); axis tight; axis xy
% [mx,mx_i]=max(nanmean(abs(TFR),2));
% text(max(t),freq_of_interest(mx_i),[' \leftarrow ' num2str(freq_of_interest(mx_i)) 'Hz'],'fontsize',12,'fontweight','bold')
% title('TFR');xlabel('Time (ms) ');ylabel('Frequency (Hz')
