function [K,MU,P,yerr,Psup]=KAL_FILT(F,G,H,U,Pold,MUold,Q,R,ymeas)    
%vanilla Kalman filter

    %%Prediction step
    mu_minus = F*MUold+G*U;
    P_minus = F*Pold*F'+Q;

    %%Kalman gain
    Psup=H*P_minus*H'+R;
    Psup=0.5*(Psup+Psup');
    K = P_minus*H'*(Psup)^(-1);

    %%Measurement update
    yerr = ymeas-H*mu_minus;
    MU = mu_minus+K*yerr;
    P = (eye(size(K*H))-K*H)*P_minus;
end

