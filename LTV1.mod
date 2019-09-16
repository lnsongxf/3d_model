

var 
b_e        // entrepreneurial debt
g_H        // housing investment adjustment cost
g_H_1      // first derivative of g_H with respect to I_H
IH         // housing investment 
PH         // Profit of housing capital producing firm
A          // productivity
b_m        // mortgage debt
C          // aggregate consumption
C_m        // consumption of borrowers (impatient households)
C_s        // consumption of savers (patient households)
D          // aggregate deposits
def_rate_e // corporate default rate
def_rate_m // mortgage default rate//def_rate_F // default rate of corporate banks
def_rate_H // default rate of mortgage banks//E_F        // equity invested in corporate banks
def_rate_F
epsilonH
epsilonF
G_e        // share of entrepreneurial capital belonging to firms that default (BGG parameter)
G_e_1      // First derivative of G_e with respect to omega_e (BGG parameter)//G_F        // share of corporate loans belonging to corporate banks that default (BGG parameter)
G_F        // share of mortgage loans belonging to mortgage banks that default (BGG parameter)
F_pF        //cdf for bank default
G_H        // share of mortgage loans belonging to mortgage banks that default (BGG parameter)
F_pH
g_I        // Capital investment adjustment cost 
g_I_1      // First derivative of g_I with respect to I_K
G_m        // share of housing belonging to households that default (BGG parameter)
G_m_1      // First derivative of G_m with respect to omega_m (BGG parameter)
Gamma_e    // Share of gross corporate revenues going to the bank (BGG parameter)
Gamma_e_1  // First derivative of Gamma_e with respect to omega_e(BGG parameter)//Gamma_F    // Share of gross corporate bank revenues going to depositors(BGG parameter)//Gamma_F_1  // First derivative of Gamma_F with respect to omega_F(BGG parameter)
Gamma_F
Gamma_F_1
Gamma_H    // Share of gross mortgage bank revenues going to depositors(BGG parameter)
Gamma_H_1  // First derivative of Gamma_H with respect to omega_H(BGG parameter)
Gamma_m    // Share of gross returns from housing going to the bank(BGG parameter)
Gamma_m_1  // First derivative of Gamma_m with respect to omega_m(BGG parameter)
F_pi       // cdf of probability of default
F_pe        //probability of defaultng households
//F_pi_1
//F_pe_1
H          // Aggregate housing supply
H_m        // Housing used by borrowers
H_s        // Housing used by savers
I          // Aggregate corporate investment
K          // Aggregate capital stock 
L          // Aggregate labour supply
L_m        // Borrowers' labour supply
L_s        // Savers' labour supply       
Lambda_m   // LM on the budget constraint in the household borrower's problem
Lambda_s   // LM on the budget constraint in the household saver's problem
n_b        // Bankers' net worth
n_e        // Entrepreneurs' net worth
omega_bar_e	// Idiosyncratic productivity shock below which the corporate borrower defaults
omega_bar_F // Idiosyncratic corporate bank loan return shock below which the corporate bank defaults	
omega_bar_H	// Idiosyncratic mortgage bank loan return shock below which the mortgage bank defaults
omega_bar_m // Idiosyncratic housing return shock below which the mortgage borrower defaults	
PI          // Profit of business capital producing firm 
q_H         // Housing price
q_K         // Capital price
R_D         // Deposit interest rate
R_F         // Corporate loan interest rate NOMINAL
R_H         // Aggregate financial return on housing (i.e. excluding imputed rents)
r_K         // Capital rental rate
R_K         // Capital rate of return
R_m         // Mortgage interest rate NOMINAL

R_tilde_F 	// Aggregate return on a diversified corporate loan portfolio (i.e. portfolio return after accounting for loan losses)
R_tilde_H 	// Aggregate return on a diversified housing loan portfolio (i.e. portfolio return after accounting for loan losses)
//rho_F       // Corporate bank return on equity
//rho_H       // Mortgage bank return on equity
Tr_H
Tr_F
Tr          // Total deposit insurance payments to banks//Tr_F        // deposit insurance payments to corporate banksTr_H        // deposit insurance payments to mortgage banks
UC_m        // Borrowers' utility from consumption
UC_m_1      // Borrowers' marginal utility from consumption
UC_s        // Savers' utility from consumption
UC_s_1      // Savers' marginal utility from consumption
UH_m        // Borrowers' utility from housing services
UH_m_1      // Borrowers' marginal utility from housing services
UH_s        // Savers' utility from housing services
UH_s_1      // Savers' marginal utility from housing services
UL_m        // Borrowers' disutility from work
UL_m_1      // Borrowers' marginal disutility from work
UL_s        // Savers' disutility from work
UL_s_1      // Savers' marginal disutility from work
Util_m      // Total borrower utility
Util_s      // Total saver utility
w           // Wage rate
W_b         // Bankers' wealth (pre-dividend)
Pr_H        //Profit of the bank
//rho         //return on bak equity
W_e         // Entrepreneurs' wealth (pre-dividend)
x_e         // Corporate leverage
x_m         // Household leverage
xi_e        // LM on bank's participation constraint in the entrepreneurs' problem
xi_m        // LM on bank's participation constraint in the household borrower's problem
//xi_b          // LM on bank's BS constraint 
Y           // Output
Y_obs       // All _obs variables are transformations to aid plotting. Transformations may differ across variables depending on what makes sense.
R_D_obs     // Please look at particular _obs variable definition in the code to see the exact transformation.
R_m_obs
R_H_obs
R_F_obs
H_m_obs
H_s_obs
b_m_obs
C_obs
C_m_obs
C_s_obs
D_obs
//E_F_obs
I_obs
K_obs
L_obs
L_m_obs
L_s_obs
n_b_obs
n_e_obs
q_H_obs
q_K_obs
r_K_obs
R_K_obs
R_tilde_F_obs
R_tilde_H_obs
Tr_obs
//Tr_F_obs
//Tr_H_obs
w_obs
x_e_obs
x_m_obs
Vs           // Value function of household savers
Vm           // Value function of household borrowers
Ve           // Value function of entrepreneurs
Vb           // Value function of bankers
EJ           // Shock (housing preference)
EK           // Shock (Capital Investment)
EH           // Shock (Housing Investment)
ESe          // Shock (entrepreneur risk)
ESm          // Shock (housing risk)
ESH          // Shock (mortgage bank risk)
ESF          // Shock (corporate bank risk)
EWe          // Shock (Entrepreneur net-worth)
EWb          // Shock (Banker net-worth)
EdH          // Shock (housing depreciation)
EdK          // Shock (capital depreciation)
//ERW        //RW shock
phi_F        // Minimum capital ratio for corporate banks
phi_H        // Minimum capital ratio for mortgage banks
phi
b_e_obs
bsp_F        // Corporate loan return spread
bsp_H        // Household loan return spread
Y_net        // Output net of default costs
Y_net_obs
Y_net_2
Y_net_2_obs
res_chk
deltaH       // Housing depreciation rate
deltaK       // Capital depreciation rate
GDP          // GDP
GDP_obs
H_obs
IH_obs
res_H
m_e          // Verification cost of entrepreneurs
m_m          // Verification cost of borrowing households
av_def       // Average default for all banks
RDsp         // Deposit rate spread over the risk free rate
R_DD         // Deposit rate return net of bank default costs
W_m          // Household borrower net worth
W_m_obs
b_tot        // Total debt in the economy
b_tot_obs
C_b          // Dividend payments by bankers
C_e          // Dividend payments by entrepreneurs
A1
B1
C1
D1
R_mi
R_Fi
phib
rho_F
rho_H
rho_F_obs
rho_H_obs
DY_NET
DY_NET_OBS

credit_gap_obs
bsp_H_obs
bsp_F_obs

//variables for measurement equations
dy_data
dq_H_data
dbe_data
dbm_data
bank_rate_data
b_to_Y_data
d_b_to_Y_data
int_rate_HH_data
bsp_H_data
bsp_F_data
dw_data
dinve_data
dc_data
CR_data
;

parameters sigma_epsiHd sigma_epsiHk sigma_epsiA sigma_epsiJ sigma_epsiK sigma_epsiH 
sigma_epsiSe sigma_epsiSm sigma_epsiSF sigma_epsiSH sigma_epsiWb sigma_epsiRW sigma_epsiWe pp  
hab Cyphi_H Cyphi_F phi_Fs phi_Hs alphaa delta_K delta_H betta_m betta_s mu_m mu_e mu_F mu_H 
sigma_e1 sigma_m1 sigma_F sigma_H varphi_s varphi_m v_s v_m chi_b chi_e eta  a_e a_s a_b   psi_i psi_h 
rhoA rhoJ rhoK rhoH rhoSe rhoSm rhoSF rhoSH rhoWb rhoWe rhoHd rhoHk rhoRW zeta1 zetae epsilonH1s epsilonF1s tau taue rp rpe phiinf1 kappa nu psib phis LTVHrule LTVFrule
gamma_y gamma_w gamma_inve gamma_c gamma_db gamma_bspH


;
    
varexo epsiA  epsiJ epsiK epsiSe epsiSm epsiSF epsiSH epsiWb  epsiWe epsiH epsiHd epsiHk epsilon_phi;
//===============================
//SET PARAMETER VALUES
load LTV1_parameter_values.mat;
//M=mortgage bank, E=entrepreneur, H=household, F=corporate bank
sigma_epsiHd=set_sigma_epsiHd;  //  shock to housing depreciation cost in EdH and EdK
sigma_epsiHk=set_sigma_epsiHk; // REDUNDANT, same as Hd in H and K sectors
 sigma_epsiA= set_sigma_epsiA; // productivity shock in A
 sigma_epsiJ= set_sigma_epsiJ; // housing preference shock in EJ
 sigma_epsiK = set_sigma_epsiK; // capital investment shock in EK
sigma_epsiH =set_sigma_epsiH ; // housing investment shock in EH
sigma_epsiSe=set_sigma_epsiSe; // entrepreneur risk shock in ESe
sigma_epsiSm=set_sigma_epsiSm; // housing risk shock in ESm
  sigma_epsiSF=  set_sigma_epsiSF; // corporate bank risk shock in ESF===> why is this so large?R 
 sigma_epsiSH= set_sigma_epsiSH; // mortgage bank risk shock is ESH
  sigma_epsiWb=set_sigma_epsiWb;// banker net worth shock in EWb
  sigma_epsiRW=set_sigma_epsiRW;// REDUNDANT, RW shock in ERW? what is this? 
  sigma_epsiWe=  set_sigma_epsiWe; //  entrepreneur net worth shock in EWe
sigma_epsilon_phi=set_sigma_epsilon_phi;
 pp  = set_pp ; // transaction cost when recovering funds from failed banks
betta_s=set_betta_s; // patient HH discount factor
betta_m=set_betta_m; // impatient HH discount factor
hab =set_hab; // habit formation
Cyphi_H =set_Cyphi_H ; // REDUNDANT? parameter on capital requirement on mortgage banks
phi_Hs=set_phi_Hs; // REDUNDANT? same equation as above, equation shut off? 
Cyphi_F=set_Cyphi_F; // REDUNDANT? parameter on capital requirement on corporate banks
phi_Fs =set_phi_Fs  ; // REDUNDANT? appears in the same equation as above, the equation is shut off. 
alphaa=set_alphaa; // share of capital in output 
rp=set_rp; // loan repayment rate of impatient households
rpe=set_rpe; // loan repayment rate of entrepreneurs
delta_K= set_delta_K; // depreciation rate of housing & capital shocks
delta_H=set_delta_H; // same as above
mu_m=set_mu_m;
mu_e=set_mu_e;
mu_F=set_mu_F;
mu_H=set_mu_H;
sigma_e1=set_sigma_e1;
sigma_m1=set_sigma_m1;
sigma_F=set_sigma_F;
sigma_H=set_sigma_H;// all 8 are hyperparameters appearing in default cdfs?
varphi_s=set_varphi_s; //  patient HH preference parameter in utility
 varphi_m= set_varphi_m; // impatient HH preference parameter in utility
v_s=set_v_s; // patient HH preference parameter in utility 
v_m=set_v_m; // impatient HH preference parameter in utility
chi_b=set_chi_b; // banker preference parameter in utility
chi_e=set_chi_e; // entrepreneur preference parameter in utility 
  eta = set_eta; // patient HH  inverse frisch elasticity of labor supply 
 a_e =set_a_e; // ??? 
a_s=set_a_s; // ??? appears in the budget constraint of savers
 a_b = set_a_b; //  ??? appears in the budget constraint of borrowers
psi_i=set_psi_i; // parameter is capital adjustment cost function 
psi_h=set_psi_h; // same as above. Why are there two of these? 
rhoA =set_rhoA;
rhoJ =set_rhoJ ;
rhoK =set_rhoK;
 rhoSe = set_rhoSe;
rhoSm =set_rhoSm; 
rhoSF =set_rhoSF; 
rhoSH =set_rhoSH; 
 rhoWb = set_rhoWb; 
rhoWe =set_rhoWe; 
rhoHd =set_rhoHd; 
rhoH=set_rhoH;
rhoHk =set_rhoHk ; 
rhoRW =set_rhoRW ; //shock persistence parameters, same notation as standard deviation
zeta1=set_zeta1; // Interest rate stickiness, same as below, appears in  R_m R_F
zetae=set_zetae; // Interest rate stickiness,  appears in 51 FOC for bank business lending, 51 is not in use 
epsilonH1s=set_epsilonH1s;//this is the part of LTV rule that is not procyclical
epsilonF1s=set_epsilonF1s;// ??? =epsilonH  & =epsilonF respectively. What are these two? 
tau=set_tau; // appears in A1, B1, C1, D1, R_mi R_Fi
taue=set_taue;
phiinf1=set_phiinf1; //reaction to inflation--why is this so low? 
 kappa = set_kappa  ; // interest rate smoothing
nu=set_nu; // ??? appears in A1, C1
 psib = set_psib; //??? appears in A1, C1
phis=set_phis; /// ???? this is equated to phi, phi is not used anywhere??  is it the same as phi_F and phi_H maybe? 
LTVHrule=set_LTVHrule;
LTVFrule=set_LTVFrule;
gamma_y=set_gamma_y;
gamma_w=set_gamma_w;
gamma_inve=set_gamma_inve;
gamma_c=set_gamma_c;
gamma_db=set_gamma_db;
gamma_bspH=set_gamma_bspH;
//************************************************************
model;
//************************************************************
// Households 
//************************************************************
//*******************************
// Utility functions
//*******************************
//1
UC_s = log(C_s-hab*C_s(-1)); 
//2
UC_m = log(C_m-hab*C_m(-1)); 
//3
UL_s = varphi_s*L_s^(1+eta)/(1+eta); 
//4
UL_m = varphi_m*L_m^(1+eta)/(1+eta); 
//5
UH_s = EJ*v_s*log(H_s(-1)); 
//6
UH_m = EJ*v_m*log(H_m(-1)); 
//7
Util_s = UC_s - UL_s + UH_s;
//8
Util_m = UC_m - UL_m + UH_m;
//9
UC_s_1 = 1/(C_s-hab*C_s(-1));
//10
UC_m_1 = 1/(C_m-hab*C_m(-1));
//11
UL_s_1 = varphi_s*L_s^(eta);
//12
UL_m_1 = varphi_m*L_m^(eta);
//13
UH_s_1 = EJ*v_s/(H_s(-1));
//14
UH_m_1 = EJ*v_m/(H_m(-1));
//*******************************
// WELFARE
//*******************************
// 15
Vs=Util_s+betta_s*Vs(+1);
// 16
Vm=Util_m+betta_m*Vm(+1);
// 17
Ve=n_e;
// 18
Vb=n_b;
//*******************************
// Savers
//*******************************
// 19 foc C_s
Lambda_s = UC_s_1;
// 20 foc L_s
UL_s_1 = w*Lambda_s;
// 21 foc d
Lambda_s = betta_s*Lambda_s(1)*R_DD(1);
//22 Definition of R_DD (deposit return net of bank default costs)
R_DD=R_D(-1)*(1-pp*(av_def));
// 23 foc H_s
Lambda_s*(q_H) = betta_s*UH_s_1(1) + betta_s*Lambda_s(1)*(1-deltaH(1))*q_H(1);
// 24 Budget Constraint
  C_s - C_e - C_b + q_H*(H_s-(1-deltaH)*H_s(-1))+ D = w*L_s +R_DD*D(-1) - Tr*a_s + PI + PH +(1-EWe)*(1-Gamma_e)*R_K*q_K(-1)*K(-1) + (1-EWb)*(W_b); 
//*******************************
// Borrowers
//*******************************
//25 Default cut off borrowers depends on leverage
omega_bar_m = x_m(-1)/R_H;
//26 Household leverage definition
x_m = (R_m)*b_m/(H_m*q_H);
//27 Housing rate of return
R_H = (1-deltaH)*q_H/q_H(-1);
//28 foc C_m    
Lambda_m = UC_m_1;
//29 foc L_m
UL_m_1 = w*Lambda_m;
// 30--ltv limit borrowers
(b_m-b_m(-1)*(1-rp)*(1-F_pi))*R_m=epsilonH*(H_m-H_m(-1)*(1-delta_H))*q_H(+1);
// 31--ltv limit entrepreneurs
(b_e-b_e(-1)*(1-rpe)*(1-F_pe))*R_F=epsilonF*(K-K(-1)*(1-delta_K))*q_K(+1);
// 32--ltv rule borrowers
epsilonH=epsilonH1s;
//exp(epsilonH-epsilonH1s)=(b_m/steady_state(b_m)^(-LTVHrule));
// 33--ltv rule entrepreneurs
epsilonF=epsilonF1s;
//exp(epsilonF-epsilonF1s)=(b_e/steady_state(b_e)^(-LTVFrule));
// 34 foc H_m
betta_m*UH_m_1(1) - Lambda_m*(q_H) + betta_m*Lambda_m(1)*((1-G_m(1))*R_H(1)*q_H)+(xi_m)*(epsilonH*(q_H(+1)))-xi_m(+1)*(epsilonH*(q_H(+2)))*(1-deltaH) = 0;  //needs to be changed                                       

// 35 foc b_m 
Lambda_m-betta_m*Lambda_m(1)*(R_m*(1-F_pi(+1)))-xi_m*R_m+xi_m(+1)*(1-rp)*(1-F_pi(+1))*R_m(+1)=0;///needs to be changed


// 36 Budget Constraint Borrowers
C_m + q_H*H_m - (1-Gamma_m)*R_H*q_H(-1)*H_m(-1) = w*L_m + b_m-Tr*a_b;
//************************************************************
// Entrepreneurs
//************************************************************
// 37 Rate of return to capital
R_K = (r_K+(1-deltaK)*q_K)/q_K(-1);
// 38 Default cut off entrepreneurs depends on leverage
omega_bar_e = x_e(-1)/R_K;
// 39 Firm leverage definition
x_e = (R_F)*(q_K*K-n_e)/(q_K*K);
//  FOC x_e
//- Gamma_e_1(1) + xi_e*((1-Gamma_F(1))*(Gamma_e_1(1)- m_e*G_e_1(1)));             
// 40 Foc k
(((1-G_e(1))*R_K(1)*q_K))*betta_s*(Lambda_s(+1)/Lambda_s)-(q_K*R_F*(1-F_pe(+1)))+xi_e*(epsilonF*q_K(+1)-R_F)-xi_e(+1)*(epsilonF*q_K(+2)*(1-deltaK)-(1-rpe)*(1-F_pe(+1))*R_F(+1))=0 ;//need to be changed for pc

// 41 Wealth dynamics entrepreneurs  
W_e= EWe*(1-Gamma_e)*R_K*q_K(-1)*K(-1);
// 42 Evolution individual net worth 
n_e = (1-chi_e)*W_e;
// 43 Corporate dividends
C_e=chi_e*W_e;


//***********************************************************
// Bankers
//************************************************************
// 42 Equalization of expected rates of return on equity invested in corporate and mortgage banks
W_b=Pr_H;

// 45 Banker net worth allocated to lending

n_b = (1-chi_b)*W_b;

// 46 Bank dividends
C_b=chi_b*(W_b);
//************************************************************
// Banks
//************************************************************
// Definitions
// 47 Default threshold corporate bank
omega_bar_F = (1-phi_F(-1))*(R_D(-1))/R_tilde_F;
// 48 Default threshold mortgage bank
omega_bar_H = (1-phi_H(-1))*R_D(-1)/R_tilde_H;




// 49
Pr_H(+1) =((1-G_H(+1))*((1-F_pi(+1))*b_m*R_m+G_m(+1)*(1 - mu_m)*(b_m*(R_m)/omega_bar_m(+1))))+(1-G_F(+1))*(((1-F_pe(+1))*b_e*R_F+G_e(+1)*(1 - mu_e)*(b_e*(R_F)/omega_bar_e)))-(1-F_pH(+1))*(b_m*(1-phib))*R_D-(1-F_pF(+1))*(b_e*(1-phib))*R_D;

//50 FOC for bank mortgage lending


phib=n_b/((((R_mi/R_m)^(-tau))*b_m+((R_Fi/R_F)^(-tau))*b_e));

A1 = (betta_s*zeta1)*(Lambda_s(+1)*((1-F_pH(+1))*(R_D)+nu*(phib/phi)^(1-psib)/(((R_mi/R_m)^(-tau))*b_m+((R_Fi/R_F)^(-tau))*b_e))*((R_m)^(tau))*b_m)+(betta_s*zeta1)*A1(+1);
B1=(betta_s*zeta1)*Lambda_s(+1)*((1-G_H(+1))*(((1-F_pi(+1)))+(G_m(+1)*(1 - mu_m)/omega_bar_m(+1))))*((((R_m)^(tau))*b_m))+(betta_s*zeta1)*B1(+1);

C1=(betta_s*zeta1)*(Lambda_s(+1)*((1-F_pF(+1))*(R_D)+(nu*(phib/phi)^(1-psib)/(((R_mi/R_m)^(-tau))*b_m+((R_Fi/R_F)^(-tau))*b_e)))*((R_F)^(tau))*b_e)+(betta_s*zeta1)*C1(+1);

D1=(betta_s*zeta1)*Lambda_s(+1)*((1-G_F(+1))*(((1-F_pe(+1)))+(G_e(+1)*(1 - mu_e)/omega_bar_e(+1))))*((((R_F)^(tau))*b_e))+(betta_s*zeta1)*D1(+1);

R_mi=(tau/(tau-1))*(A1)/(B1);
R_Fi=(tau/(tau-1))*(C1)/(D1);
R_m=((1-zeta1)*R_mi^(1-tau)+zeta1*R_m(-1)^(1-tau))^(1/(1-tau));
R_F=((1-zeta1)*R_Fi^(1-tau)+zeta1*R_F(-1)^(1-tau))^(1/(1-tau));



// 53 Rate of return on household mortgage loans
R_tilde_H = (Gamma_m - mu_m*G_m)*R_H*q_H(-1)*H_m(-1)/b_m(-1);
// 54 Rate of return on corporate loans
R_tilde_F = (Gamma_e-mu_e*G_e)*R_K*q_K(-1)*K(-1)/(q_K(-1)*K(-1)-n_e(-1));
// 55 Balance Sheet Bank
n_b + D = b_m + (q_K*K - n_e); 
//************************************************************
// Consumption good production
//************************************************************
// 56 Output
Y = A*K(-1)^(alphaa)*L^(1-alphaa);
// 57 Capital rental rate
r_K = alphaa*Y/K(-1);
// 58 Wage rate
w = (1-alphaa)*Y/L;
//************************************************************
// Capital good production with CEE adjustment costs 
// (depend on I/I(-1)) - Christiano, Eichenbaum and Evans (1995) - JPE
//************************************************************
// 59 Capital adj cost function
g_I = EK*psi_i/2*(I/I(-1)-1)^2;
// 60 foc of (57) wrt I
g_I_1 = EK*psi_i*(I/I(-1)-1);
// 61 Capital stock evolution
K = (1-deltaK)*K(-1) + I*(1-g_I );
// 62 Capital price
q_K = 1 + g_I + I/I(-1)*g_I_1 - betta_s*Lambda_s(1)/Lambda_s*(I(1)/I)^2*g_I_1(1);
// 63 Flow profit of capital producers
PI = q_K*I - (1 + g_I)*I; 
//*******************************
// Inside equity market
//*******************************
//*******************************
// Goods market
//*******************************
// 65 Aggregate consumption definition
C = C_s + C_m;
//*******************************
// Labour market 
//*******************************
// 66 Aggregate labour supply definition
L = L_s + L_m;
//*******************************
// Housing supply with CEE Adj costs
//*******************************
// 67 Capital adj cost function
g_H = EH*psi_h/2*(IH/IH(-1)-1)^2;
// 68 foc of (63) wrt IH
g_H_1 = EH*psi_h*(IH/IH(-1)-1);
// 69 Housing stock evolution
H = (1-deltaH)*H(-1) + IH*(1-g_H);
// 70 Housing price with CEE Adj costs
q_H = 1 + g_H + IH/IH(-1)*g_H_1 - betta_s*Lambda_s(1)/Lambda_s*(IH(1)/IH)^2*g_H_1(1);
// 71 Flow profits by housing producers 
PH= q_H*IH - (1 + g_H)*IH; 
// 72 Aggregate housing market clearing (Supply (LHS) = Demand (RHS))
H = H_m + H_s; 
//*******************************
// Deposit insurance agency
//*******************************
//  Deposit insurance transfers to corporate bank depositors whose banks have gone bankrupt
//Tr_F = (omega_bar_F - Gamma_F + mu_F*G_F)*R_tilde_F*(q_K(-1)*K(-1)-(1-chi_e)*W_e(-1));
// 70 Deposit insurance transfers to household bank depositors whose banks have gone bankrupt
Tr_H = (omega_bar_H - Gamma_H + mu_H*G_H)*((((1-F_pi)*b_m(-1)*R_m(-1)+G_m*(1 - mu_m)*(b_m(-1)*(R_m(-1))/omega_bar_m))));
Tr_F = (omega_bar_F - Gamma_F + mu_F*G_F)*(((1-F_pe)*b_e(-1)*R_F(-1)+G_e*(1 - mu_e)*(b_e(-1)*(R_F(-1))/omega_bar_e)));

// 73
// 71 Total deposit insurance transfers
Tr = Tr_F + Tr_H;



//**************************************************************
//Taylor rule
//**************************************************************



//************************************************************
// Distribution functions
//************************************************************
// 74
Gamma_m   = normcdf((log(omega_bar_m)- (ESm*sigma_m1)^2/2)/(ESm*sigma_m1)) + omega_bar_m*(1-normcdf((log(omega_bar_m)+ESm*sigma_m1^2/2)/(ESm*sigma_m1)));
// 75
Gamma_e   = normcdf((log(omega_bar_e) - (ESe*sigma_e1)^2/2)/(ESe*sigma_e1)) + omega_bar_e*(1-normcdf((log(omega_bar_e)+ESe*sigma_e1 ^2/2)/(ESe*sigma_e1)));
// 76
Gamma_F   = normcdf((log(omega_bar_F) - (ESF*sigma_F)^2/2)/(ESF*sigma_F)) + omega_bar_F*(1-normcdf((log(omega_bar_F)+ESF*sigma_F^2/2)/(ESF*sigma_F)));
// 77
Gamma_H   = normcdf((log(omega_bar_H)-(ESH*sigma_H)^2/2)/(ESH*sigma_H)) + omega_bar_H*(1-normcdf((log(omega_bar_H)+ESH*sigma_H^2/2)/(ESH*sigma_H)));
// 78
Gamma_m_1 = (normpdf((log(omega_bar_m)-(ESm*sigma_m1)^2/2)/(ESm*sigma_m1))-omega_bar_m*normpdf((log(omega_bar_m)+ESm*sigma_m1^2/2)/(ESm*sigma_m1)))/(ESm*sigma_m1*omega_bar_m) + (1-normcdf((log(omega_bar_m)+(ESm*sigma_m1)^2/2)/(ESm*sigma_m1)));
// 79
Gamma_e_1 = (normpdf((log(omega_bar_e)-(ESe*sigma_e1) ^2/2)/(ESe*sigma_e1))-omega_bar_e*normpdf((log(omega_bar_e)+ESe*sigma_e1 ^2/2)/(ESe*sigma_e1)))/(ESe*sigma_e1 *omega_bar_e) + (1-normcdf((log(omega_bar_e)+(ESe*sigma_e1)^2/2)/(ESe*sigma_e1)));
// 80
Gamma_F_1 = (normpdf((log(omega_bar_F)-(ESF*sigma_F)^2/2)/(ESF*sigma_F))-omega_bar_F*normpdf((log(omega_bar_F)+ESF*sigma_F^2/2)/(ESF*sigma_F)))/(ESF*sigma_F*omega_bar_F) + (1-normcdf((log(omega_bar_F)+(ESF*sigma_F)^2/2)/(ESF*sigma_F)));
// 81
Gamma_H_1 = (normpdf((log(omega_bar_H)-(ESH*sigma_H)^2/2)/(ESH*sigma_H))-omega_bar_H*normpdf((log(omega_bar_H)+ESH*sigma_H^2/2)/(ESH*sigma_H)))/(ESH*sigma_H*omega_bar_H) + (1-normcdf((log(omega_bar_H)+(ESH*sigma_H)^2/2)/(ESH*sigma_H)));
// 82
G_m   = normcdf((log(omega_bar_m)-ESm*sigma_m1^2/2)/(ESm*sigma_m1));
// 83
G_e   = normcdf((log(omega_bar_e)-ESe*sigma_e1 ^2/2)/(ESe*sigma_e1));
// 84
G_F   = normcdf((log(omega_bar_F)-ESF*sigma_F^2/2)/(ESF*sigma_F));
// 85
G_H   = normcdf((log(omega_bar_H)-ESH*sigma_H^2/2)/(ESH*sigma_H));
// 86
G_m_1 = normpdf((log(omega_bar_m)-ESm*sigma_m1^2/2)/(ESm*sigma_m1))/(ESm*sigma_m1*omega_bar_m);
// 87
F_pi = normcdf((log(omega_bar_m)+ESm*sigma_m1^2/2)/(ESm*sigma_m1));//probability of defaultng households
//F_pi_1 = normpdf((log(omega_bar_m)+ESm*sigma_m1^2/2)/(ESm*sigma_m1))/(ESm*sigma_m1*omega_bar_m);
// 88
F_pe = normcdf((log(omega_bar_e)+ESm*sigma_e1^2/2)/(ESm*sigma_e1));//probability of defaultng households
//F_pe_1 = normpdf((log(omega_bar_e)+ESm*sigma_e1^2/2)/(ESm*sigma_e1))/(ESe*sigma_e1 *omega_bar_e);//probability of defaultng households
// 89
F_pH   = normcdf((log(omega_bar_H)+sigma_H^2/2)/sigma_H);
// 90
F_pF   = normcdf((log(omega_bar_F)+sigma_F^2/2)/sigma_F);
// 91
G_e_1 = normpdf((log(omega_bar_e)-ESe*sigma_e1 ^2/2)/(ESe*sigma_e1))/(ESe*sigma_e1 *omega_bar_e);

// 92 Capital requirement on mortgage banks
exp(phi_H-phi_Hs) = (b_m/steady_state(b_m))^(Cyphi_H)+epsilon_phi(-1);


// 93 Capital requirement on corporate banks
exp(phi_F-phi_Fs) = (b_e/steady_state(b_e))^(Cyphi_F)+epsilon_phi(-1);

phi = phis;

// Capital depreciation
// 94
deltaK = delta_K+EdK;
// 95 Housing depreciation
deltaH = delta_H+EdH;

//*********************
//****Other definitions
//*********************
// NFC loans and total loans
//96
b_e=(q_K*K-n_e);
//97
b_tot=b_e+b_m;
// Return on loan portfolio spreads
//98
bsp_H = 400*(R_tilde_H(1)-R_D);
//99
bsp_F = 400*(R_tilde_F(1)-R_D);
//GDP DEFINITIONS (Y_net is the main definition of GDP used in the model)
//100
Y_net = C + (chi_b)*W_b + (chi_e)*W_e + I + IH;
//101
Y_net_2 = Y - (R_D(-1)*pp*av_def*D(-1)+m_e*G_e*R_K*q_K(-1)*K(-1) + m_m*G_m*R_H*q_H(-1)*H_m(-1)+ mu_F*G_H*(Gamma_e - m_e*G_e)*R_K*q_K(-1)*K(-1) + mu_H*G_H*(Gamma_m - m_m*G_m)*R_H*q_H(-1)*H_m(-1));
GDP = Y + UL_s_1*H_s/Lambda_s + UL_m_1*H_m/Lambda_m;
// Monitoring cost
m_e = mu_e;
m_m = mu_m;
//Average default of banks
av_def = ((1-phi_F)*b_e*def_rate_F + (1-phi_H)*b_m*def_rate_H)/D;

// Deposit spread
RDsp = 400*(R_D-R_DD(+1));
// Borrowers wealth
W_m = (1-Gamma_m)*R_H*q_H(-1)*H_m(-1) + w*L_m;

//************************************************************
// Reporting adjunct variables
//************************************************************
def_rate_m = normcdf((log(omega_bar_m)+ESm*sigma_m1^2/2)/(ESm*sigma_m1))*400;
def_rate_e = normcdf((log(omega_bar_e)+ESe*sigma_e1^2/2)/(ESe*sigma_e1))*400;
def_rate_H = normcdf((log(omega_bar_H)+ESH*sigma_H^2/2)/(ESH*sigma_H))*400;
def_rate_F = normcdf((log(omega_bar_F)+ESF*sigma_F^2/2)/(ESF*sigma_F))*400;
Y_obs = log(Y/steady_state(Y))*100; 
R_D_obs = (R_D-steady_state(R_D))*400; 
R_m_obs = (R_m-steady_state(R_m))*400;
R_H_obs = log(R_H/steady_state(R_H))*400;
R_F_obs = (R_F-steady_state(R_F))*400;
H_m_obs = log(H_m/steady_state(H_m))*100;
H_s_obs = log(H_s/steady_state(H_s))*100;
b_m_obs = log(b_m/steady_state(b_m))*100;
b_e_obs = log(b_e/steady_state(b_e))*100;
C_obs = log(C/steady_state(C))*100;
C_m_obs = log(C_m/steady_state(C_m))*100;
C_s_obs = log(C_s/steady_state(C_s))*100;
D_obs = log(D/steady_state(D))*100;
//E_F_obs = log(E_F/steady_state(E_F))*100;
I_obs = log(I/steady_state(I))*100;
K_obs = log(K/steady_state(K))*100;
L_obs = log(L/steady_state(L))*100;
L_m_obs = log(L_m/steady_state(L_m))*100;
L_s_obs = log(L_s/steady_state(L_s))*100;
n_b_obs = log(n_b/steady_state(n_b))*100;
n_e_obs = log(n_e/steady_state(n_e))*100;
q_H_obs = log(q_H/steady_state(q_H))*100;
q_K_obs = log(q_K/steady_state(q_K))*100;
r_K_obs = log(r_K/steady_state(r_K))*400;
R_K_obs = (R_K-steady_state(R_K))*400;
R_tilde_F_obs = (R_tilde_F-steady_state(R_tilde_F))*400;
R_tilde_H_obs = (R_tilde_H-steady_state(R_tilde_H))*400;

Tr_obs = log(Tr/steady_state(Tr))*100;
//Tr_F_obs = log(Tr_F/steady_state(Tr_F))*100;
//Tr_H_obs = log(Tr_H/steady_state(Tr_H))*100;
w_obs = log(w/steady_state(w))*100;
x_e_obs = log(x_e/steady_state(x_e))*100;
x_m_obs = log(x_m/steady_state(x_m))*100;
b_tot_obs=log(b_tot/steady_state(b_tot))*100;
W_m_obs = log(W_m/steady_state(W_m))*100;
GDP_obs = log(GDP/steady_state(GDP))*100;
H_obs = log(H/steady_state(H))*100;
IH_obs = log(IH/steady_state(IH))*100;
Y_net_2_obs = log(Y_net_2/steady_state(Y_net_2))*100;
Y_net_obs = log(Y_net/steady_state(Y_net))*100;


// Residuals
res_chk = (R_D(-1)*pp*av_def*D(-1)+C + chi_b*W_b + chi_e*W_e + I + IH + (m_e*G_e*R_K*q_K(-1)*K(-1) + m_m*G_m*R_H*q_H(-1)*H_m(-1)+ mu_F*G_H*(Gamma_e - m_e*G_e)*R_K*q_K(-1)*K(-1) + mu_H*G_H*(Gamma_m - m_m*G_m)*R_H*q_H(-1)*H_m(-1)))/Y;
res_H = 100*(H - H_s - H_m)/H;


//************************************************************
// SHOCKS 
// Note that each of these shocks potentially has a component 
// that hits immediately and a 'news shock' component
//************************************************************
log(A)   = rhoA*log(A(-1)) - epsiA;
log(EJ)  = rhoJ*log(EJ(-1)) + epsiJ;
log(EK)  = rhoK*log(EK(-1)) + epsiK;
log(EH)  = rhoH*log(EH(-1)) + epsiH;
log(ESe) = rhoSe*log(ESe(-1)) + epsiSe;
log(ESm) = rhoSm*log(ESm(-1)) + epsiSm;
//Bank risk shocks are perfectly correlated across banks--> this also has been changed it seems?
log(ESF) = rhoSF*log(ESF(-1)) + epsiSF;
log(ESH) = rhoSH*log(ESH(-1)) + epsiSH;
log(EWe) = rhoWe*log(EWe(-1)) - epsiWe;
log(EWb) = rhoWb*log(EWb(-1)) - epsiWb;
//log(ERW) = rhoRW*log(ERW(-1))+epsiRW;
//Depreciation shocks are perfectly correlated across H and K sectors--> this has been changed?
EdH      = rhoHd*EdH(-1) + epsiHd;
EdK      = rhoHk*EdK(-1) + epsiHk;

// 49 Rate of return on corporate bank equity
rho_F = (1-Gamma_F)*R_tilde_F/phi_F(-1);
// 50 Rate of return on mortgage bank equity
rho_H =(1-Gamma_H)*R_tilde_H/phi_H(-1);
rho_F_obs = (rho_F-steady_state(rho_F))*400;
rho_H_obs = (rho_H-steady_state(rho_H))*400;


//DY_NET=0;
DY_NET_OBS=0;
DY_NET=Y_net-Y_net(-1);
//DY_NET_OBS=log(DY_NET/steady_state(DY_NET))*100;


credit_gap_obs=100*log( (b_tot/steady_state(b_tot)) / (Y_net/steady_state(Y_net)));
//bsp_H_obs=100*log(bsp_H/steady_state(bsp_H));
bsp_H_obs=400*log((R_tilde_H(1)-R_D)/(steady_state(R_tilde_H)-steady_state(R_D)));
bsp_F_obs=400*log((R_tilde_F(1)-R_D)/(steady_state(R_tilde_F)-steady_state(R_D)));

//MEASUREMENT EQUATIONS
dy_data=Y_net_obs-Y_net_obs(-1)+gamma_y;
dq_H_data=q_H_obs-q_H_obs(-1);
dbe_data=b_e_obs-b_e_obs(-1)+gamma_db;
dbm_data=b_m_obs-b_m_obs(-1)+gamma_db;
bank_rate_data=R_D_obs+4*steady_state(R_D);
b_to_Y_data-b_to_Y_data(-1)=credit_gap_obs-credit_gap_obs(-1);
d_b_to_Y_data=b_to_Y_data-b_to_Y_data(-1);
int_rate_HH_data=R_m_obs+4*steady_state(R_m);
//d_bsp_H_data=(bsp_H_obs-bsp_H_obs(-1))/4;
bsp_H_data=bsp_H+gamma_bspH;
bsp_F_data-bsp_F_data(-1)= bsp_F_obs-bsp_F_obs(-1);
dw_data=w_obs-w_obs(-1)+gamma_w;
//dinve_data=(I_obs-I_obs(-1))+(IH_obs-IH_obs(-1))+gamma_y;
dinve_data=I_obs-I_obs(-1)+gamma_inve;
dc_data=C_obs-C_obs(-1)+gamma_c;
CR_data=100*phi_H;
//rho_F=roe;
end; 
//************************************************************
// Here specify the innovations
// Note that it is possible to shut down a shock by assigning it 
// a zero variance.
//************************************************************
shocks;
var epsiA   = sigma_epsiA;
var epsiJ   = sigma_epsiJ;
var epsiK   = sigma_epsiK;
var epsiH   = sigma_epsiSH;
var epsiSe  = sigma_epsiSe;
var epsiSm  = sigma_epsiSm;
var epsiSF  = sigma_epsiSF;
var epsiSH  = sigma_epsiSH;
var epsiWb  = sigma_epsiWb;
var epsiWe  = sigma_epsiWe; 

var epsiHd  = sigma_epsiHd;
var epsiHk  = sigma_epsiHk;
//var epsiRW  = sigma_epsiRW;
//var epsilon_phi = sigma_epsilon_phi;
end;

//resid;
//steady;
//check;


estimated_params;
betta_m,   0.9842  ,0.95,0.99,BETA_PDF,0.97,0.005;
rp, 0.035,0,1,BETA_PDF,0.05,0.025;
rpe,   0.04,0,1,BETA_PDF,0.05,0.025;
delta_H, 0.01 ,0,1,BETA_PDF,0.05,0.025;
delta_K,0.03 ,0,1,BETA_PDF,0.05,0.025;
v_s,0.1  ,0.01,1,BETA_PDF,0.5,0.25;
v_m,  0.273 ,0,1,BETA_PDF,0.5,0.25;
//chi_e, 0.0760,0,1,UNIFORM_PDF,0.5, 0.288675134594813;
//chi_b, 0.005 ,0.01,1,UNIFORM_PDF,0.5, 0.288675134594813;
gamma_y,  0.21 ,1e-6,5,UNIFORM_PDF,0,5;
gamma_w,  0.8 ,1e-6,5,UNIFORM_PDF,0,5;
gamma_inve,  0.17 ,1e-6,5,UNIFORM_PDF,0,5;
gamma_c,  0.24 ,1e-6,5,UNIFORM_PDF,0,5;

hab,0.4167,0.2,0.99,BETA_PDF,0.5,0.2;
psi_i,  1.4340,1,5,NORMAL_PDF,2,2;
psi_h, 4.9446 ,1,5,NORMAL_PDF,2,2;
varphi_s,  1.2989 ,0,5,NORMAL_PDF,1,0.5;
varphi_m, 0.7290 ,0,5,NORMAL_PDF,1,0.5;
sigma_e1,0.318,1e-5,1,         UNIFORM_PDF,5, 2.88675134595;
sigma_m1,0.1,1e-5,1,         UNIFORM_PDF,5, 2.88675134595;
sigma_F,0.0363,1e-5,1,         UNIFORM_PDF,5, 2.88675134595;
sigma_H,0.0363,1e-5,1,         UNIFORM_PDF,5, 2.88675134595;
alphaa,0.3 ,0,1,BETA_PDF,0.3,0.05;
mu_m, 0.3 ,0,1,BETA_PDF,0.3,0.05;
mu_e, 0.3,0,1,BETA_PDF,0.3,0.05;
mu_F,0.3,0,1,BETA_PDF,0.3,0.05;
mu_H, 0.3,0,1,BETA_PDF,0.3,0.05;




stderr epsiHd,0.1,0,99,          UNIFORM_PDF,5, 2.88675134595;
stderr epsiHk,0.01,0,99,    UNIFORM_PDF,5, 2.88675134595;
stderr epsiA,0.1,0,99,     UNIFORM_PDF,5, 2.88675134595;
stderr epsiJ,0.1,0,99,   UNIFORM_PDF,5, 2.88675134595;
stderr epsiK,0.1,0,99,    UNIFORM_PDF,5, 2.88675134595;
//stderr epsiH,0.1,0,99,    UNIFORM_PDF,5, 2.88675134595;
stderr epsiSe,0.1,0,99,       UNIFORM_PDF,5, 2.88675134595;
stderr epsiSm,0.1,0,99,  UNIFORM_PDF,5, 2.88675134595;
stderr epsiSF,0.1,0,99,     UNIFORM_PDF,5, 2.88675134595;
stderr epsiSH,0.1,0,99,   UNIFORM_PDF,5, 2.88675134595;
stderr epsiWe,0.1,0,99, UNIFORM_PDF,5, 2.88675134595;
stderr epsiWb,0.1,0,99,  UNIFORM_PDF,5, 2.88675134595;
//stderr epsiWb,0.1,0,99,  INV_GAMMA2_PDF,3,0.0005;
//stderr epsilon_phi,0.1,0,99, UNIFORM_PDF,5, 2.88675134595;

rhoHd,0.01,0,0.9999,BETA_PDF,0.5,0.2;
rhoHk,0.01,0,0.9999,BETA_PDF,0.5,0.2;
rhoA,0.01,0,0.9999,BETA_PDF,0.5,0.2;
rhoJ,0.01,0,0.9999,BETA_PDF,0.5,0.2;
rhoK,0.01,0,0.9999,BETA_PDF,0.5,0.2;
//rhoH,0.01,0,0.99,BETA_PDF,0.5,0.2;
rhoSe,0.01,0,0.9999,BETA_PDF,0.5,0.2;
rhoSm,0.01,0,0.99,BETA_PDF,0.5,0.2;                                                                                                                                                          
rhoSF,0.01,0,0.99,BETA_PDF,0.5,0.2; 
rhoSH,0.01,0,0.99,BETA_PDF,0.5,0.2;
rhoWe,0.01,0,0.99,BETA_PDF,0.5,0.2;
rhoWb,0.01,0,0.99,BETA_PDF,0.5,0.2;


//zeta1,0.2723,0,1,BETA_PDF,0.5,0.2;
end;
varobs dy_data,dq_H_data,int_rate_HH_data,bank_rate_data,bsp_H_data,dw_data,dinve_data,dbm_data,dc_data;//-->dc_data should be taken out potentially    ,dinve_data,dw_data,;,bsp_H_data,,dinve_data,;,dc_data, dbm_data;,, ;int_rate_HH_data , bank_rate_data ,dbe_data bank_rate_data  , int_rate_HH_data;bank_rate_data  ,dbe_data; dbe_data,dbm_data bsp_H_data;; , //////////////dq_H_data, b_to_Y_data, int_rate_HH_data, bank_rate_data,bsp_H_data;,b_to_Y_data ,b_to_Y_data,bsp_H_data


estimation(datafile='estimation_dataset_quarterly.mat', 
mode_compute=0,
mode_file='LTV1_mode.mat',
//optim=('TolFun',1e-4),
optim=('TolFun',1e-5,'TolX',1e-5,'Display','iter','MaxIter',1000,'Hessian','bfgs'),
//optimset=('TolFun',1e-5,'TolX',1e-5,'Display','iter','MaxIter',5000),
//optim=('TolFun',10e-4,'TolX',10e-4,'Display','iter','UseParallel',1,'MaxIter',9999999,'MaxFunEvals',9999),
kalman_algo=4,
//lik_init=1,
nograph,nodiagnostic,
first_obs=2,
nobs=75,
presample=2,
mh_nblocks=1,
mh_replic=0,
mh_jscale=0.28,
mh_drop=0.2);


//calib_smoother(datafile='estimation_dataset_quarterly.mat')dy_data, Y_net;
//
//shock_decomposition(parameter_set=calibration,datafile='estimation_dataset_quarterly.mat')dy_data, dq_H_data,b_to_Y_data;




//varobs dy_data, dq_H_data, b_to_Y_data,bank_rate_data,int_rate_HH_data;  
//shock_decomposition(datafile='estimation_dataset_quarterly.mat',parameter_set=posterior_mode)dy_data, dq_H_data,d_b_to_Y_data,int_rate_HH_data,bank_rate_data,bsp_H_data,dbe_data,dbm_data,dw_data,dc_data,dinve_data;// , dq_H_data,d_b_to_Y_data,int_rate_HH_data,bank_rate_data,bsp_H_data ,d_b_to_Y_data; 
//shock_decomposition(datafile='estimation_dataset_quarterly.mat',parameter_set=posterior_mode)dy_data,dw_data,dc_data,dinve_data,bsp_H_data;// , dq_H_data,d_b_to_Y_data,int_rate_HH_data,bank_rate_data,bsp_H_data ,d_b_to_Y_data; 



stoch_simul(order=1,irf=40,nograph,periods=10000,drop=5000);//Y_net_obs,dy_data,dw_data,dinve_data,dq_H_data,bank_rate_data,bsp_H_data,b_m_obs;

//R_D,R_D_obs,dy_data, dq_H_data,dbe_data,dbm_data,bank_rate_data,def_rate_e,def_rate_m,def_rate_H,def_rate_F,b_tot,b_tot_obs,Y_net,Y_net_obs,credit_gap_obs,w_obs,L_obs,b_m_obs,Tr_obs,w,L,b_m,Tr;//,Y_net,Y_net_obs b_tot b_tot_obs,b_m,b_m_obs,b_e,b_e_obs,dbe_data,dbm_data,bsp_H_data,bsp_F_data;

 


//varobs dy_data, dq_H_data, b_to_Y_data; // dbm_data;,dbe_data , ; , bank_rate_data; 
//shock_decomposition(datafile='estimation_dataset_quarterly.mat', parameter_set=calibration)dy_data;//, dq_H_data,b_to_Y_data;
