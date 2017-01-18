function [xhat,Pest,yerr]=Burger_Observer(F,G,H,U,Pold,xhatold,K,Q,R,ymeas)
    %prediction
    xhatminus=F*xhatold+G*U;
    P_minus=F*Pold*F'+Q;
    
    %measurement update
    yerr=ymeas-H*xhatminus;
    xhat=xhatminus+K*yerr;
    Pest=P_minus-K*H*P_minus-P_minus*H'*K'+K*(H*P_minus*H'+R)*K';
end