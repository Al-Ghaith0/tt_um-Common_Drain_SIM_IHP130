eps0 = 8.854187e-12;
eps_ox = 3.9;
%NMOS
Cjorgatn = 1.3846e-11;
Cjorgatdn = 1e-9;
Cjorstin = 2.4557e-11;
Cjorstidn = 1e-9;
Cjorbotn = 0.85856;
Cjorbotdn = 0.001;
toxon = 7.4256e-9;
toxovon = 7.4256e-9;
meun = 0.46687 * 0.85058;
Lovn = 6.3866e-8;
Coxn = eps0 * eps_ox / toxon;
Coxovn = eps0 * eps_ox / toxovon;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pmos
Cjorgatdp = 1e-9;
Cjorgatp = 3.1e-11;
Cjorstip = 4.2583e-11;
Cjorstidp = 1e-9;
Cjorbotp = 0.0008085;
Cjorbotdp = 0.001;
toxop = 6.945e-9;
toxovop = 6.945e-9;
meup = 0.018428 * 1.8162;
Lovp = 3.471e-8;
Coxp = eps0 * eps_ox / toxop;
Coxovp = eps0 * eps_ox / toxovop;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NMOS XM17
W17 = 1e-5;
L17 = 7.2e-7;
As17 = 3.4e-12;
Ps17 = 2.068e-5;
Ad17 = 3.4e-12;
Pd17 = 2.068e-5;
Cgsn17 = 2/3 * W17 * L17 * Coxn + W17 * Lovn * Coxovn
Cgdn17 = W17 * Lovn * Coxovn
Cjsn17 = Cjorbotn / (sqrt(1 + 0))
Csbn17 = (As17 + W17 * L17) * Cjsn17 + Ps17 * Cjorstin + W17 * Cjorgatn
Cjdn17 = Cjorbotdn / (sqrt(1 + 2.056))
Cdbn17 = (Ad17 * Cjdn17 + Pd17 * Cjorstidn + W17 * Cjorgatdn)
I17 = 2.872e-6;
Veffe17 = sqrt (I17 * 2 / (meun * Coxn * W17 / L17))
gm = 2 * I17 / Veffe17
rds17 = 1.75/I17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NMOS XM18
W18 = 1e-5;
L18 = 7.2e-7;
As18 = 3.4e-12;
Ps18 = 2.068e-5;
Ad18 = 3.4e-12;
Pd18 = 2.068e-5;
Cgsn18 = 2/3 * W18 * L18 * Coxn + W18 * Lovn * Coxovn
Cgdn18 = W18 * Lovn * Coxovn
Cjsn18 = Cjorbotn / (sqrt(1 + 0))
Csbn18 = (As18 + W18 * L17) * Cjsn17 + Ps17 * Cjorstin + W18 * Cjorgatn
Cjdn18 = Cjorbotdn / (sqrt(1 + 2.056))
Cdbn18 = (Ad18 * Cjdn18 + Pd18 * Cjorstidn + W18 * Cjorgatdn)
I18 = 2.872e-6;
Veffe18 = sqrt (I18 * 2 / (meun * Coxn * W18 / L18))
rds18 = 1.75/I18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TF derivation
C1 = Cgdn17
C2 = Cgsn17
C3 = Cgdn18
C4 = Cgsn18
C5 = Csbn17 + Cdbn18
R20 = (20/gm+rds18)/(1+gm*rds18);
R30 = 10/gm+rds18/(1+rds18*gm)+10/gm*rds18/(1+rds18*gm)*gm;
R32 = 10/gm+(rds18/(1+rds18*gm)*20/gm)/(rds18/(1+rds18*gm)+20/gm)
b1 = 20/gm*C1 + R20*C2 + R30*C3 + 10/gm*C4 + C5*rds18/(1+rds18*gm);
b2 = 20/gm*rds18/(1+rds18*gm)*C1*C2 + 20/gm*R30*C1*C3 + 20/gm*10/gm*C1*C4 + (20*rds18)/(rds18+1)*C1*C5 + R20*R32*C2*C3 + R20*20/gm*C2*C4 + R20*((rds18/(1+gm*rds18))*20/gm)/((rds18/(1+gm*rds18)+20/gm))*C2*C5 + R30*((rds18/(1+gm*rds18))*10/gm)/((rds18/(1+gm*rds18)+10/gm))*C3*C4 + R30*((rds18/(1+gm*rds18))*10/gm)/((rds18/(1+gm*rds18)+10/gm))*C3*C5 + 10/gm*rds18/(1+rds18*gm)*C4*C5;
a0 = rds18/(1/gm+rds18)
a1 = R20*C2*rds18/(1/gm+rds18) + R30*C3*rds18/(1/gm+(rds18*10/gm)/(rds18+10/gm)) + rds18/(1/gm+rds18)*C4
nm = [a1,a0];
de = [b2,b1,1];
TF = tf(nm,de)
bode(TF);
margin(TF);
