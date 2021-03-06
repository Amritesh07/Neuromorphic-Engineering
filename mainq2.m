% Init Neuron Type values
%C(pF) kz(?S/V) Er(mV) Et(mV) a (KHz) b(nS) c(mV) d(pA) vpeak(mV)
type = 'RS';
if type == 'RS'
   C =  100E-12;   Kz=  0.7E-6;   Er=  -60E-3;   Et=  -40E-3;   a =  0.03E3;
   b =  -2E-9;   c =  -50E-3;   d =  100E-12;   Vpeak=35E-3;
elseif type == 'IB'
   C =  150E-12;   Kz=  1.2E-6;   Er=  -75E-3;   Et=  -45E-3;   a =  0.01E3;
   b =  5E-9;   c =  -56E-3;   d =  130E-12;   Vpeak=50E-3;
elseif type == 'CH'
   C =  50E-12;   Kz=  1.5E-6;   Er=  -60E-3;   Et=  -40E-3;   a =  0.03E3;
   b =  1E-9;   c =  -40E-3;   d =  150E-12;   Vpeak=25E-3;
end

N = 10; % No of neurons
h = 0.1E-3; % time step
Tmax = 0.5; %max time
M = Tmax/dt; % max time index i.e no. of columns (no. of delta_T time intervals that amount to T)

% input current matrix
    %Ic = gL*(VT-EL) ; %considering Ic = 1
    %alpha = 0.1;
I = zeros(N,M);
for i = 1:N
    I(i,1:M) = 400*10^-12 ; %(1 + alpha*i)*Ic;
    % I(i,1:M) = 500*10^-12 ;
    % I(i,1:M) = 600*10^-12 ;
end

y = zeros(N,M);
U = zeros(N,M);
% init condition
y(:,1) = Er;
U(:,1) = 0;


F = @(t,V,U,I) (Kz*(V-Et).*(V-Er)-U+I)/C ; 
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