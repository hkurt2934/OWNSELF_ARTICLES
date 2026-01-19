function out = run_IM_case(enable_FTC)

%% PARAMETERS
Ts = 1e-4;  Tsim = 3;  N = Tsim/Ts;
t = (0:N-1)*Ts;

Rs=1.405; Rr=1.395; Ls=0.0058; Lr=0.0058; Lm=0.0055;
P=2; J=0.02; B=0.001;

w_ref = 150;     % rad/s
id_ref = 2;
t_fault = 1;

Kp_w=0.4; Ki_w=15;
Kp_i=5;  Ki_i=300;

%% STATES
id=0; iq=0; psi_r=0.8; w=0;
int_w=0; int_id=0; int_iq=0;

%% LOGS
w_log=zeros(1,N); iq_log=w_log; id_log=w_log;
Te_log=w_log; psi_log=w_log;

%% LOOP
for k=1:N

    % Rotor resistance fault
    if t(k) > t_fault
        Rr_eff = 3*Rr;
    else
        Rr_eff = Rr;
    end

    % FTC selection
    if enable_FTC
        Rr_hat = Rr_eff;
    else
        Rr_hat = Rr;
    end

    % Speed loop
    w_err = w_ref - w;
    int_w = int_w + Ki_w*Ts*w_err;
    iq_ref = Kp_w*w_err + int_w;

    % Current loops
    vd = Kp_i*(id_ref-id) + int_id;
    vq = Kp_i*(iq_ref-iq) + int_iq;

    id = id + Ts*(vd - Rs*id)/Ls;
    iq = iq + Ts*(vq - Rs*iq)/Ls;

    % Flux observer (affected by fault)
    psi_r = psi_r + Ts*(Lm/Lr*id - (Rr_hat/Lr)*psi_r);

    % Torque & mechanics
    Te = (3/2)*P*(Lm/Lr)*psi_r*iq;
    w = w + Ts*(Te - B*w)/J;

    % Logs
    w_log(k)=w;
    iq_log(k)=iq;
    id_log(k)=id;
    Te_log(k)=Te;
    psi_log(k)=psi_r;
end

out.t=t;
out.w=w_log;
out.iq=iq_log;
out.id=id_log;
out.Te=Te_log;
out.psi=psi_log;
end
