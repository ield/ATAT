clear;
f = 2.45e9;
lambda = 3e8 / f;
k = 2*pi/lambda;
%% Initial restriction
% We start fixing b/lambda at 0.65 so that it is posible to use the grpahs
% in fig 15.17
b = 0.65*lambda;
% k*b % to check 3.83 < kb < 5.33

% However, phi_ph will be the last to calculate, because it only deppends
% on l. It will be adjusted with this restriction phi_ex and with the
% corrugated horn data phi_fl

%% Calculating phi_ex
% It is known that 2b/lambda = 1.3. Therefore, it is seen that for the
% majority of the b/a relations, phi_ph = 320º = -40º.
phi_ex = -40;

% It is chosen the ratio b/a = 1.27, so that
a = b/1.27;
% k*a % to check 1.84 < ka < 3.83

% This gives a c = -7.5. If it were that case that the TE was dominant, it
% would be necessary to increase the b/a relation. If the TM was dominant,
% it woulld be necessary to decrease the b/a relation.
%% Calculating phi_fl
% Parameters from diagram of corrugated horn (diapo 16)
% It is wanted an intensity of -12 dB at 35º. 
relFielInt = 10^(-12/20);
% It is obtained that for a s = 0.8, 2pia0/lambda*sin(35) = 8
s = 0.8;
theta_35 = 35*pi/180;
a_0 = 8/sin(theta_35)*lambda/(2*pi);
% For this value of a_0, it is calculated L_long, the L in the diapo 16,
% useful to calculate theta
L_long = a_0^2/(2*lambda*s);
theta = atan(a_0 / L_long);
% Therefore, a_0*lambda = 2.2, so looking at 15.17.a, for b/lambda = 0.065,
% it is obtained y_axis = 0.184
phi_fl = 0.184*360/tan(theta);

%% Calculating phi_ph and l
% Therefore, it is obtained that 
phi_ph = 360 + phi_ex - phi_fl;

% Therefore, from the graph 15.17.b
l_p = phi_ph/360*lambda/0.548;

%% Print all the values
a = a/lambda
b = b/lambda
a_0 = a_0/lambda
l_p = l_p/lambda
lambda
theta = theta*180/pi

%%
% How to calculate the lambda of the guide for a given r for te11 and tm11
r = 0.6*lambda;             % Radius of the guide

p_te = 1.841;               % p11 (taken from Pozar)
p_tm = 3.832;               % p11 (taken from Pozar)

lg_te = 2*pi*r/p_te;        % lambda g for te11       
lg_tm = 2*pi*r/p_tm;        % lambda g for tm11



%%
zero_1 = 3.8317;
b_a = 1.27;
lambda = 122.4;
kb = 4;
a = 20*log10(abs(besselj(1, zero_1/b_a))^2 / sqrt(abs(1-(zero_1/kb)^2)))
