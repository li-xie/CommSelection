function dy=chem_conc_RC(t,y,para)

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

% the term gM(R, B)/gM_max
M_coef = BN/(RN_M+BN).*RN_M./(RN_M+1) + RN_M./(RN_M+BN).*BN./(BN+1);
% the term gH(R)/gH_max
H_coef=RN_H./(RN_H+1);
dB = -c_BM_M*gM_max*M_coef + r_B_H*gH_max*H_coef;
dR = -c_RM*gM_max*M_coef*M - c_RH*gH_max*H_coef*H;

dy=[dB; dR];