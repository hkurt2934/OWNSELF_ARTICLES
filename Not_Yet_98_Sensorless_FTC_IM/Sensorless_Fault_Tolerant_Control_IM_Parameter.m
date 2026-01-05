f_base   = 200;             % Hz electrical
p        = 4;               % poles
wb       = 2*pi*f_base;     % base electrical rad/s
Ts       = 1e-5;

% Stator/rotor (phase, per unit-length)
Rs   = 0.435;               % ohm
Rr   = 0.816;               % ohm (referred to stator)
Ls   = 0.002;               % H   (stator self)
Lr   = 0.002;               % H   (rotor self, referred)
Lm   = 0.0693;              % H   magnetizing
J    = 0.015;               % kg·m^2
B    = 0.0001;              % N·m·s/rad viscous
Vdc  = 300;                 % V
fsw  = 10e3;                % Hz PWM
Ts   = 100e-6;              % s control & observer sample time

sigma = 1 - (Lm^2)/(Ls*Lr);
Tr    = Lr / Rr;                    % rotor time constant
Leak  = Ls - (Lm^2)/Lr;             % leakage

Param_Limit = [1,1]; %???
Ki = 100; %???
Kp = 5; %???
Zp = 1; %??? Number of Poles
