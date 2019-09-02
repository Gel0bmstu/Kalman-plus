clear all;

% Row count
rowKol = 17875;
colNum = 0;

% Open source file
f = dlmread('6.csv', ';', [0 colNum rowKol colNum])';
dt = 1;

% Setup initialize data
eOptX(1) = f(1);
eOptX(2) = f(2);

wp = f(1);
wpp = f(1);


wf(1) = f(1);
wf(2) = f(2);

we(1) = f(1);
we(2) = f(2);

wn(1) = f(1);
wn(2) = f(2);

% Kalman's coefficient for default filter
K2 = 0.9995;
K = 0.997;

wn_2(1) = f(1);
wn_2(2) = f(2) * (1 - K2) + wn_2(1) * K2;

% Gyro error [sigma]
sigmaAx = 0.27*4.848e-6;

sp = 0;
sr = f(1) + f(2);

% Step to build line for complex
step = 50;
a = 0;

eOptX(1) =  0.27*4.848e-13;
eOptX(2) =  0.27*4.848e-13;

Karr(1) = 0.27*4.848e-6;
Karr(2) = 0.27*4.848e-6;

Karr1(1) = 0.27*4.848e-6;
Karr1(2) = 0.27*4.848e-6;

zap = 5;
buf = 0;
buf1 = 0;

for t = 3:dt:rowKol
    % At the beginnig algorothm working like a normal Kamlman filter
    if (t < step)
      wn(t) = f(t) * (1 - K) + wn(t - 1) * K;
      we(t) = wn(t);
      wn_2(t) = f(t) * (1 - K2) + wn_2(t-1) * K2;
      wn_3(t) = wn(t);
      eOptX(t)=sqrt((sigmaAx^4)*(eOptX(t-1)^2+2e-7^4)/(2e-7^4+eOptX(t-1)^2+sigmaAx^4));
      eOptX1(t)=sqrt((sigmaAx^4)*(eOptX(t-1)^2+2e5^4)/(2e5^4+eOptX(t-1)^2+sigmaAx^4));
      continue;
    end;
      
    % Then, if the current time is a multiple of "Step" (50 for default)
    if ( mod(t, step) == 0)
      wp = wn(t - 1);
      wpp = wn(t - zap);
      
      % We are building a line and consider it a reference
      y = @(x) (-((t - zap * dt) * wp - (t - dt) * wpp) - (wpp - wp) * x) / ((dt * (zap - 1)) );
    end;
      
    sigmaU = abs(f(t) - y(t)) / 100;
    sigmaU2 = abs(f(t) - y(t)) / 1000;
    
    eOptX(t)=(sigmaU^4)*(eOptX(t-1) + sigmaAx^4)/(sigmaU^4+eOptX(t-1)+sigmaAx^4);
    K=eOptX(t)/sigmaU^4;
    
    Karr(t) = abs(f(t) - y(t));
    Karr1(t) = abs(f(t) - y(t)) / 700;
    
    wn(t) = f(t) * (1 - K) + wn(t - 1) * K;
    wn_2(t) = f(t) * (1 - K2) + wn_2(t-1) * K2;
    we(t) = y(t);
end;

% crutches :)
f(rowKol) = f(rowKol);
we(rowKol + 1) = f(rowKol);
wn(rowKol + 1) = f(rowKol);
wn_2(rowKol + 1) = f(rowKol);
wn_3(rowKol + 1) = f(rowKol);
eOptX(rowKol + 1) = eOptX(rowKol);
k = 1:rowKol + 1;

Karr(rowKol + 1) = Karr(rowKol);
Karr1(rowKol + 1) = Karr1(rowKol);

% Plot the result
plot(k, f, k, wn, k, wn_2)
xlabel('Время t, [сек]');
ylabel('Приращение угла ψ, [рад / сек]');
legend('Сигнал с гироскопа', 'Калман модифицированный', 'Калман обычный');
grid on;
