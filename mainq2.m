%some constants describing the model
gL = 30*10^-9 ;% mhos
EL = -70*10^-3 ;% Volts
C  = 300*10^-12 ;% Farad
VT = 20*10^-3 ;% Volts

% inputs
M = 500/0.1 ; % no. of columns (no. of delta_T time intervals that amount to T)  
N = 2 ; % no. of rows (no. of neurons) 
T = 0.5 ; % the time upto which the expt is performed
h = 0.0001 ; % delta_T

% input current matrix
Ic = gL*(VT-EL) ; %considering Ic = 1
alpha = 0.1;
I = zeros(N,M);
for i = 1:N
    I(i,1:M) = (1 + alpha*i)*Ic;
end

A = kz/C ; 
B = gL*EL ;

y = zeros(N,M);
% init condition
y(:,1) = EL; 
F = @(t,V,U,I) A*(V-Et)*(V-Er)+(-U+I)/C ; 
G = @(t,V,U) a*(b*(V-Er)-U) ; 

%runge kutta implementation
for i = 1:M-1 
    k1 = F(i,       y(:,i),          U(:,i),          (I(:,i)+I(:,i+1))/2 );
    L1 = G(i,       y(:,i),          U(:,i));
    k2 = F(i+0.5*h, y(:,i)+0.5*h*k1, U(:,i)+0.5*h*L1, (I(:,i)+I(:,i+1))/2 );
    L2 = G(i+0.5*h, y(:,i)+0.5*h*k1, U(:,i)+0.5*h*L1);
    k3 = F(i+0.5*h, y(:,i)+0.5*h*k2, U(:,i)+0.5*h*L2, (I(:,i)+I(:,i+1))/2 );
    L3 = G(i+0.5*h, y(:,i)+0.5*h*k2, U(:,i)+0.5*h*L2);
    k4 = F(i+h,     y(:,i)+h*k3,     U(:,i)+h*L3,     (I(:,i)+I(:,i+1))/2 );
    L4 = G(i+h,     y(:,i)+h*k3,     U(:,i)+h*L3);
    
    y(:,i+1) = y(:,i) + (1/6)*(k1+2*(k2+k3)+k4)*h ;
    U(:,i+1) = U(:,i) + (1/6)*(L1+2*(L2+L3)+L4)*h ;
    
    y((y(:,i+1)>=Vpeak),i+1) = c;
    U((y(:,i+1)>=Vpeak),i+1) = U((y(:,i+1)>=Vpeak),i+1) + d;
end

t = 0.1:0.1:500;

figure(1)
for i=1:N
    plot(t,y(i,:))
    hold on
end
title('Spiking Neuron')
xlabel('Time (in ms)')
ylabel('Membrane Potential (in Volts)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1 = find(y(2,4000:5000)==-0.07);
t2 = find(y(4,4000:5000)==-0.07);
t3 = find(y(6,4000:5000)==-0.07);
t4 = find(y(8,4000:5000)==-0.07);

avg_t1 = t1(2)-t1(1);
avg_t2 = t2(2)-t2(1);
avg_t3 = t3(2)-t3(1);
avg_t4 = t4(2)-t4(1);

avg_t = [avg_t1 avg_t2 avg_t3 avg_t4];
Iappk = [I(2) I(4) I(6) I(8)];

figure(2)
plot(Iappk,avg_t)

title('Avg time of spiking vs time')
ylabel('Time (in ms)')
xlabel('Iappk (in nA)')
