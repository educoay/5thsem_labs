% Generation of Line Codes, PSD and Probability of Error

clc
close all

% input paramater

N=10;
a = round(rand(1,N))
A=5;
Tb=1;
fs=100;

% unipolar NRZ
u=[];
for k=1:N
    u = [u A*a(k)*ones(1,fs)];
end
t = linspace(0, N*Tb, length(u));

%unipolar RZ
u_rz=[];
for k=1:N
    c=ones(1,fs/2);
    b=zeros(1,fs/2);
    p=[c b];
    u_rz=[u_rz A*a(k)*p];
end

%polar NRZ
P=[];
for k=1:N
    P=[P ((-1)^(a(k)+1))*A*ones(1,fs)];
end

%polar rz
p_rz=[];
for k=1:N
    c=ones(1,fs/2);
    b=zeros(1,fs/2);
    p=[c b];
    p_rz=[p_rz ((-1)^(a(k)+1))*A*p];
end

%bipolar NRZ
B=[];
count=-1;
for k=1:N;
   if a(k)==1
        if count==-1
            B=[B A*a(k)*ones(1,fs)];
            count=-1;
        end
   else
       B=[B A*a(k)*ones(1,fs)];
   end
end

%bipolar RZ/AMI RZ
b_rz=[];
count=-1;
for k=1:N
   if a(k)==1
        if count==-1
            b_rz=[b_rz A*a(k)*ones(1,fs/2) zeros(1,fs/2)];
            count=-1;
   else
       b_rz=[b_rz A*a(k)*ones(1,fs) zeros(1,fs/2)];
       count=-1
   end
else
    b_rz=[b_rz A*a(k)*ones(1,fs)];
end
end

%split- phase or manscher code
 
m=[];
for k=1:N
    c=ones(1,fs/2);
    b=zeros(1,fs/2);
    p=[c b];
    m=[m ((-1)^(a(k)+1))*A*p];
end


t=linspace(0,N*Tb,length(u));
 
figure(1)
subplot (7,1,1);
plot(t,u)
axis([0 N*Tb -6 6])
title('unipolar NRZ')
grid on

% figure(2)
subplot(7,1,2);
plot(t, u_rz)
axis([0 N*Tb -6 6])
title('unipolar RZ')
grid on

%figure(3)
subplot (7,1,3);
plot(t,P)
axis([0 N*Tb -6 6])
title('polar NRZ')
grid on
 
 
%figure(4)
subplot (7,1,4);
plot(t,p_rz)
axis([0 N*Tb -6 6])
title('polar rz')
grid on
 
%figure(5)
subplot (7,1,5);
plot(t,B)
axis([0 N*Tb -6 6])
title('bipolar NRZ')
grid on
 
%figure(6)
subplot (7,1,6);
plot(t,b_rz)
axis([0 N*Tb -6 6])
title('bipolar RZ/AMI RZ')
grid on
 
%figure(7)
subplot (7,1,7);
plot(t,m)
axis([0 N*Tb -6 6])
title('split- phase or manscher code')
grid on

% PSD:
v=1; 
R=1; 
T=1/R; 
f=0:0.001*R:2*R;  
f= f+1e-10;

%Unipolar NRZ
s=((v^2*T/4).*(sin(pi.*f*T)./(pi.*f*T)).^2);
s(1)=s(1)+(v^2/4);
ff=0;
stem(ff,s(1),'*r','LineWidth',4)
hold on;
figure(2)
plot(f,s,'-r','LineWidth',2);
hold on;


%Manchester code
s=(v.^2.*T).*((sin(pi.*f*T/2)./(pi.*f*T/2)).^2).*(sin(pi.*f*T/2).^2);
plot(f,s,'--g','LineWidth',2);
hold on;

%Polar NRZ
s=((v^2*T).*(sin(pi.*f*T)./(pi.*f*T)).^2);
plot(f,s,'--b','LineWidth',2);
hold on;

%Bipolar RZ 
s=(v.^2.*T/4).*((sin(pi.*f*T/2)./(pi.*f*T/2)).^2).*(sin(pi.*f*T).^2);
% figure(3)
plot(f,s,'--k','LineWidth',2);
legend('Unipolar NRZ: impulse at at f=0','Unipolar NRZ','Manchestercode','PolarNRZ','Bipolar RZ/ RZ-AMI');
xlabel('Normalised frequency)');
ylabel('Power spectral density');


% PE ERROR:
 
%Unipolar NRZ
E=[0:1:25]; 
 
%Unipolar NRZ
P1=(1/2)*erfc(sqrt(E/2));
 
%polar NRZ and Manchester code has same Pe for equiprobable 1's and 0's
P2=(1/2)*erfc(sqrt(E));
 
%Bipolar RZ/ RZ-AMI
P3=(3/4)*erfc(sqrt(E/2));
figure(4)
E=10*log10(E); 
semilogy(E,P1,'-k',E,P2,'-r',E,P3,'-b','LineWidth',2)
legend('Unipolar NRZ','Polar NRZ and Manchester','Bipolar RZ/ RZ-AMI','Location','best');
xlabel('SNR per bit, Eb/No(dB)');
ylabel('Bit error probality Pe');

