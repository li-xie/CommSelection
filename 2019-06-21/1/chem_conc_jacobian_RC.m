function J=chem_conc_jacobian_RC(t, y, para)

B = y(1);
R = y(2);

gM_max = para(1);
gH_max=para(2);
c_RM = para(3);
c_RH = para(4);
K_MR = para(5);
K_HR =para(6);
K_MB = para(7);
M = para(8);
H = para(9);
c_BM_M = para(10);
r_B_H = para(11);

BN = B/K_MB;
RN_M = R/K_MR;
RN_H = R/K_HR;

% M_coef = BN/(RN_M+BN).*RN_M./(RN_M+1) + RN_M./(RN_M+BN).*BN./(BN+1);
% H_coef=RN_H./(RN_H+1);

dM_coefdB = RN_M^2/(RN_M+BN)^2/(RN_M+1) ...
    + RN_M*(1-BN^2)/(RN_M+BN)^2/(BN+1)^2;
dM_coefdR = BN^2/(RN_M+BN)^2/(BN+1) ...
    + BN*(1-RN_M^2)/(RN_M+BN)^2/(RN_M+1)^2;
dH_coefdR = 1/(RN_H+1)^2/K_HR;
dM_coefdB = dM_coefdB/K_MB;
dM_coefdR = dM_coefdR/K_MR;

% dB = -c_BM_M*gM_max*M_coef + r_B_H*gH_max*H_coef;
% dR = -c_RM*gM_max*M_coef*M - c_RH*gH_max*H_coef*H;
J11 = c_BM_M*gM_max*dM_coefdB;
J12 = -c_BM_M*gM_max*dM_coefdR + r_B_H*gH_max*dH_coefdR;
J21 = -c_RM*gM_max*dM_coefdB*M;
J22 = -c_RM*gM_max*dM_coefdR*M - c_RH*gH_max*dH_coefdR*H;
J=[J11 J12; J21 J22];