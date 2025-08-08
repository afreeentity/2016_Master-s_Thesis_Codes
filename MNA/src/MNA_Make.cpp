/*

	 Class Make as a member of MNA 

	written by Ehsan Nedaaee Oskoee*/

#include "mnamatrix.h"

void MNA_Matrix::Make(double *Vx){

 	double idt;
 	if(Delt > 0.)  idt = 1./Delt;
	for(int i=1; i<=size; i++) {
		for(int j=1; j<=mna[i].size(); j++) mna[i](mna[i].row(j))=0.;
		rval[i] = 0.;
		fvec[i]=0.;
	}

 /* Elements of  V1 */ 
	auto Vt01 = [] (double t, double vs){
		double val;
		val = 100.0;
		return val*vs;
	};
	rval[41] = rval[41] + Vt01(tm,Vs[1]);
	fvec[41] -= Vt01(tm,Vs[1]);
	mna[8](41) = mna[8][41] + 1.;
	fvec[8] +=  Vx[41];
	mna[41](8) = mna[41][8] + 1.;
	fvec[41] +=  Vx[8];

 /* Elements of  NM1*/ 

	double Un_M, Uo_M, Qn_M, Qo_M;
	double Geq_M, Ieq_M, Sr_M;
	double xp_M, gq_M,Kh_M,Kp_M;
	static int tt_M = 0;

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[2];
	Un_M += Vx[1];
	Uo_M += -uxtold[2];
	Uo_M += uxtold[1];
	Qo_M = uxtold[42];
	Qn_M = Vx[42];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[42](1) = mna[42][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](2) = mna[1][2] - idt*Geq_M;
	mna[2](1) = mna[2][1] - idt*Geq_M;
	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[42](2) = mna[42][2] +Geq_M;
	rval[2] += (Ieq_M - Qo_M)*idt;
	fvec[2] += -(Sr_M);
	mna[42](42) = mna[42][42] + 1.;
	rval[42] += Ieq_M;
	fvec[42] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM2*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[3];
	Un_M += Vx[1];
	Uo_M += -uxtold[3];
	Uo_M += uxtold[1];
	Qo_M = uxtold[43];
	Qn_M = Vx[43];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[43](1) = mna[43][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](3) = mna[1][3] - idt*Geq_M;
	mna[3](1) = mna[3][1] - idt*Geq_M;
	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[43](3) = mna[43][3] +Geq_M;
	rval[3] += (Ieq_M - Qo_M)*idt;
	fvec[3] += -(Sr_M);
	mna[43](43) = mna[43][43] + 1.;
	rval[43] += Ieq_M;
	fvec[43] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM3*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[4];
	Un_M += Vx[1];
	Uo_M += -uxtold[4];
	Uo_M += uxtold[1];
	Qo_M = uxtold[44];
	Qn_M = Vx[44];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[44](1) = mna[44][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](4) = mna[1][4] - idt*Geq_M;
	mna[4](1) = mna[4][1] - idt*Geq_M;
	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[44](4) = mna[44][4] +Geq_M;
	rval[4] += (Ieq_M - Qo_M)*idt;
	fvec[4] += -(Sr_M);
	mna[44](44) = mna[44][44] + 1.;
	rval[44] += Ieq_M;
	fvec[44] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM4*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[5];
	Un_M += Vx[1];
	Uo_M += -uxtold[5];
	Uo_M += uxtold[1];
	Qo_M = uxtold[45];
	Qn_M = Vx[45];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[45](1) = mna[45][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](5) = mna[1][5] - idt*Geq_M;
	mna[5](1) = mna[5][1] - idt*Geq_M;
	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[45](5) = mna[45][5] +Geq_M;
	rval[5] += (Ieq_M - Qo_M)*idt;
	fvec[5] += -(Sr_M);
	mna[45](45) = mna[45][45] + 1.;
	rval[45] += Ieq_M;
	fvec[45] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM5*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[6];
	Un_M += Vx[1];
	Uo_M += -uxtold[6];
	Uo_M += uxtold[1];
	Qo_M = uxtold[46];
	Qn_M = Vx[46];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[46](1) = mna[46][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](6) = mna[1][6] - idt*Geq_M;
	mna[6](1) = mna[6][1] - idt*Geq_M;
	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[46](6) = mna[46][6] +Geq_M;
	rval[6] += (Ieq_M - Qo_M)*idt;
	fvec[6] += -(Sr_M);
	mna[46](46) = mna[46][46] + 1.;
	rval[46] += Ieq_M;
	fvec[46] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM6*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[7];
	Un_M += Vx[1];
	Uo_M += -uxtold[7];
	Uo_M += uxtold[1];
	Qo_M = uxtold[47];
	Qn_M = Vx[47];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[47](1) = mna[47][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](7) = mna[1][7] - idt*Geq_M;
	mna[7](1) = mna[7][1] - idt*Geq_M;
	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[47](7) = mna[47][7] +Geq_M;
	rval[7] += (Ieq_M - Qo_M)*idt;
	fvec[7] += -(Sr_M);
	mna[47](47) = mna[47][47] + 1.;
	rval[47] += Ieq_M;
	fvec[47] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM7*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[8];
	Un_M += Vx[1];
	Uo_M += -uxtold[8];
	Uo_M += uxtold[1];
	Qo_M = uxtold[48];
	Qn_M = Vx[48];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[48](1) = mna[48][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](8) = mna[1][8] - idt*Geq_M;
	mna[8](1) = mna[8][1] - idt*Geq_M;
	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[48](8) = mna[48][8] +Geq_M;
	rval[8] += (Ieq_M - Qo_M)*idt;
	fvec[8] += -(Sr_M);
	mna[48](48) = mna[48][48] + 1.;
	rval[48] += Ieq_M;
	fvec[48] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM8*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[9];
	Un_M += Vx[1];
	Uo_M += -uxtold[9];
	Uo_M += uxtold[1];
	Qo_M = uxtold[49];
	Qn_M = Vx[49];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[49](1) = mna[49][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](9) = mna[1][9] - idt*Geq_M;
	mna[9](1) = mna[9][1] - idt*Geq_M;
	mna[9](9) = mna[9][9] + idt*Geq_M;
	mna[49](9) = mna[49][9] +Geq_M;
	rval[9] += (Ieq_M - Qo_M)*idt;
	fvec[9] += -(Sr_M);
	mna[49](49) = mna[49][49] + 1.;
	rval[49] += Ieq_M;
	fvec[49] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM9*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[11];
	Un_M += Vx[1];
	Uo_M += -uxtold[11];
	Uo_M += uxtold[1];
	Qo_M = uxtold[50];
	Qn_M = Vx[50];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[50](1) = mna[50][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](11) = mna[1][11] - idt*Geq_M;
	mna[11](1) = mna[11][1] - idt*Geq_M;
	mna[11](11) = mna[11][11] + idt*Geq_M;
	mna[50](11) = mna[50][11] +Geq_M;
	rval[11] += (Ieq_M - Qo_M)*idt;
	fvec[11] += -(Sr_M);
	mna[50](50) = mna[50][50] + 1.;
	rval[50] += Ieq_M;
	fvec[50] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM10*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[12];
	Un_M += Vx[1];
	Uo_M += -uxtold[12];
	Uo_M += uxtold[1];
	Qo_M = uxtold[51];
	Qn_M = Vx[51];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[51](1) = mna[51][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](12) = mna[1][12] - idt*Geq_M;
	mna[12](1) = mna[12][1] - idt*Geq_M;
	mna[12](12) = mna[12][12] + idt*Geq_M;
	mna[51](12) = mna[51][12] +Geq_M;
	rval[12] += (Ieq_M - Qo_M)*idt;
	fvec[12] += -(Sr_M);
	mna[51](51) = mna[51][51] + 1.;
	rval[51] += Ieq_M;
	fvec[51] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM11*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[13];
	Un_M += Vx[1];
	Uo_M += -uxtold[13];
	Uo_M += uxtold[1];
	Qo_M = uxtold[52];
	Qn_M = Vx[52];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[52](1) = mna[52][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](13) = mna[1][13] - idt*Geq_M;
	mna[13](1) = mna[13][1] - idt*Geq_M;
	mna[13](13) = mna[13][13] + idt*Geq_M;
	mna[52](13) = mna[52][13] +Geq_M;
	rval[13] += (Ieq_M - Qo_M)*idt;
	fvec[13] += -(Sr_M);
	mna[52](52) = mna[52][52] + 1.;
	rval[52] += Ieq_M;
	fvec[52] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM12*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[14];
	Un_M += Vx[1];
	Uo_M += -uxtold[14];
	Uo_M += uxtold[1];
	Qo_M = uxtold[53];
	Qn_M = Vx[53];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[53](1) = mna[53][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](14) = mna[1][14] - idt*Geq_M;
	mna[14](1) = mna[14][1] - idt*Geq_M;
	mna[14](14) = mna[14][14] + idt*Geq_M;
	mna[53](14) = mna[53][14] +Geq_M;
	rval[14] += (Ieq_M - Qo_M)*idt;
	fvec[14] += -(Sr_M);
	mna[53](53) = mna[53][53] + 1.;
	rval[53] += Ieq_M;
	fvec[53] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM13*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[15];
	Un_M += Vx[1];
	Uo_M += -uxtold[15];
	Uo_M += uxtold[1];
	Qo_M = uxtold[54];
	Qn_M = Vx[54];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[54](1) = mna[54][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](15) = mna[1][15] - idt*Geq_M;
	mna[15](1) = mna[15][1] - idt*Geq_M;
	mna[15](15) = mna[15][15] + idt*Geq_M;
	mna[54](15) = mna[54][15] +Geq_M;
	rval[15] += (Ieq_M - Qo_M)*idt;
	fvec[15] += -(Sr_M);
	mna[54](54) = mna[54][54] + 1.;
	rval[54] += Ieq_M;
	fvec[54] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM14*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[16];
	Un_M += Vx[1];
	Uo_M += -uxtold[16];
	Uo_M += uxtold[1];
	Qo_M = uxtold[55];
	Qn_M = Vx[55];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[55](1) = mna[55][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](16) = mna[1][16] - idt*Geq_M;
	mna[16](1) = mna[16][1] - idt*Geq_M;
	mna[16](16) = mna[16][16] + idt*Geq_M;
	mna[55](16) = mna[55][16] +Geq_M;
	rval[16] += (Ieq_M - Qo_M)*idt;
	fvec[16] += -(Sr_M);
	mna[55](55) = mna[55][55] + 1.;
	rval[55] += Ieq_M;
	fvec[55] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM15*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[17];
	Un_M += Vx[1];
	Uo_M += -uxtold[17];
	Uo_M += uxtold[1];
	Qo_M = uxtold[56];
	Qn_M = Vx[56];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[56](1) = mna[56][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](17) = mna[1][17] - idt*Geq_M;
	mna[17](1) = mna[17][1] - idt*Geq_M;
	mna[17](17) = mna[17][17] + idt*Geq_M;
	mna[56](17) = mna[56][17] +Geq_M;
	rval[17] += (Ieq_M - Qo_M)*idt;
	fvec[17] += -(Sr_M);
	mna[56](56) = mna[56][56] + 1.;
	rval[56] += Ieq_M;
	fvec[56] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM16*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[18];
	Un_M += Vx[1];
	Uo_M += -uxtold[18];
	Uo_M += uxtold[1];
	Qo_M = uxtold[57];
	Qn_M = Vx[57];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[57](1) = mna[57][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](18) = mna[1][18] - idt*Geq_M;
	mna[18](1) = mna[18][1] - idt*Geq_M;
	mna[18](18) = mna[18][18] + idt*Geq_M;
	mna[57](18) = mna[57][18] +Geq_M;
	rval[18] += (Ieq_M - Qo_M)*idt;
	fvec[18] += -(Sr_M);
	mna[57](57) = mna[57][57] + 1.;
	rval[57] += Ieq_M;
	fvec[57] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM17*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[19];
	Un_M += Vx[1];
	Uo_M += -uxtold[19];
	Uo_M += uxtold[1];
	Qo_M = uxtold[58];
	Qn_M = Vx[58];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[58](1) = mna[58][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](19) = mna[1][19] - idt*Geq_M;
	mna[19](1) = mna[19][1] - idt*Geq_M;
	mna[19](19) = mna[19][19] + idt*Geq_M;
	mna[58](19) = mna[58][19] +Geq_M;
	rval[19] += (Ieq_M - Qo_M)*idt;
	fvec[19] += -(Sr_M);
	mna[58](58) = mna[58][58] + 1.;
	rval[58] += Ieq_M;
	fvec[58] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM18*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[20];
	Un_M += Vx[1];
	Uo_M += -uxtold[20];
	Uo_M += uxtold[1];
	Qo_M = uxtold[59];
	Qn_M = Vx[59];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[59](1) = mna[59][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](20) = mna[1][20] - idt*Geq_M;
	mna[20](1) = mna[20][1] - idt*Geq_M;
	mna[20](20) = mna[20][20] + idt*Geq_M;
	mna[59](20) = mna[59][20] +Geq_M;
	rval[20] += (Ieq_M - Qo_M)*idt;
	fvec[20] += -(Sr_M);
	mna[59](59) = mna[59][59] + 1.;
	rval[59] += Ieq_M;
	fvec[59] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM19*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[21];
	Un_M += Vx[1];
	Uo_M += -uxtold[21];
	Uo_M += uxtold[1];
	Qo_M = uxtold[60];
	Qn_M = Vx[60];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[60](1) = mna[60][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](21) = mna[1][21] - idt*Geq_M;
	mna[21](1) = mna[21][1] - idt*Geq_M;
	mna[21](21) = mna[21][21] + idt*Geq_M;
	mna[60](21) = mna[60][21] +Geq_M;
	rval[21] += (Ieq_M - Qo_M)*idt;
	fvec[21] += -(Sr_M);
	mna[60](60) = mna[60][60] + 1.;
	rval[60] += Ieq_M;
	fvec[60] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM20*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[22];
	Un_M += Vx[1];
	Uo_M += -uxtold[22];
	Uo_M += uxtold[1];
	Qo_M = uxtold[61];
	Qn_M = Vx[61];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[61](1) = mna[61][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](22) = mna[1][22] - idt*Geq_M;
	mna[22](1) = mna[22][1] - idt*Geq_M;
	mna[22](22) = mna[22][22] + idt*Geq_M;
	mna[61](22) = mna[61][22] +Geq_M;
	rval[22] += (Ieq_M - Qo_M)*idt;
	fvec[22] += -(Sr_M);
	mna[61](61) = mna[61][61] + 1.;
	rval[61] += Ieq_M;
	fvec[61] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM21*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[23];
	Un_M += Vx[1];
	Uo_M += -uxtold[23];
	Uo_M += uxtold[1];
	Qo_M = uxtold[62];
	Qn_M = Vx[62];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[62](1) = mna[62][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](23) = mna[1][23] - idt*Geq_M;
	mna[23](1) = mna[23][1] - idt*Geq_M;
	mna[23](23) = mna[23][23] + idt*Geq_M;
	mna[62](23) = mna[62][23] +Geq_M;
	rval[23] += (Ieq_M - Qo_M)*idt;
	fvec[23] += -(Sr_M);
	mna[62](62) = mna[62][62] + 1.;
	rval[62] += Ieq_M;
	fvec[62] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM22*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[24];
	Un_M += Vx[1];
	Uo_M += -uxtold[24];
	Uo_M += uxtold[1];
	Qo_M = uxtold[63];
	Qn_M = Vx[63];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[63](1) = mna[63][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](24) = mna[1][24] - idt*Geq_M;
	mna[24](1) = mna[24][1] - idt*Geq_M;
	mna[24](24) = mna[24][24] + idt*Geq_M;
	mna[63](24) = mna[63][24] +Geq_M;
	rval[24] += (Ieq_M - Qo_M)*idt;
	fvec[24] += -(Sr_M);
	mna[63](63) = mna[63][63] + 1.;
	rval[63] += Ieq_M;
	fvec[63] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM23*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[25];
	Un_M += Vx[1];
	Uo_M += -uxtold[25];
	Uo_M += uxtold[1];
	Qo_M = uxtold[64];
	Qn_M = Vx[64];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[64](1) = mna[64][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](25) = mna[1][25] - idt*Geq_M;
	mna[25](1) = mna[25][1] - idt*Geq_M;
	mna[25](25) = mna[25][25] + idt*Geq_M;
	mna[64](25) = mna[64][25] +Geq_M;
	rval[25] += (Ieq_M - Qo_M)*idt;
	fvec[25] += -(Sr_M);
	mna[64](64) = mna[64][64] + 1.;
	rval[64] += Ieq_M;
	fvec[64] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM24*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[26];
	Un_M += Vx[1];
	Uo_M += -uxtold[26];
	Uo_M += uxtold[1];
	Qo_M = uxtold[65];
	Qn_M = Vx[65];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[65](1) = mna[65][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](26) = mna[1][26] - idt*Geq_M;
	mna[26](1) = mna[26][1] - idt*Geq_M;
	mna[26](26) = mna[26][26] + idt*Geq_M;
	mna[65](26) = mna[65][26] +Geq_M;
	rval[26] += (Ieq_M - Qo_M)*idt;
	fvec[26] += -(Sr_M);
	mna[65](65) = mna[65][65] + 1.;
	rval[65] += Ieq_M;
	fvec[65] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM25*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[28];
	Un_M += Vx[1];
	Uo_M += -uxtold[28];
	Uo_M += uxtold[1];
	Qo_M = uxtold[66];
	Qn_M = Vx[66];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[66](1) = mna[66][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](28) = mna[1][28] - idt*Geq_M;
	mna[28](1) = mna[28][1] - idt*Geq_M;
	mna[28](28) = mna[28][28] + idt*Geq_M;
	mna[66](28) = mna[66][28] +Geq_M;
	rval[28] += (Ieq_M - Qo_M)*idt;
	fvec[28] += -(Sr_M);
	mna[66](66) = mna[66][66] + 1.;
	rval[66] += Ieq_M;
	fvec[66] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM26*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[32];
	Un_M += Vx[1];
	Uo_M += -uxtold[32];
	Uo_M += uxtold[1];
	Qo_M = uxtold[67];
	Qn_M = Vx[67];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[67](1) = mna[67][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](32) = mna[1][32] - idt*Geq_M;
	mna[32](1) = mna[32][1] - idt*Geq_M;
	mna[32](32) = mna[32][32] + idt*Geq_M;
	mna[67](32) = mna[67][32] +Geq_M;
	rval[32] += (Ieq_M - Qo_M)*idt;
	fvec[32] += -(Sr_M);
	mna[67](67) = mna[67][67] + 1.;
	rval[67] += Ieq_M;
	fvec[67] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM27*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[33];
	Un_M += Vx[1];
	Uo_M += -uxtold[33];
	Uo_M += uxtold[1];
	Qo_M = uxtold[68];
	Qn_M = Vx[68];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[68](1) = mna[68][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](33) = mna[1][33] - idt*Geq_M;
	mna[33](1) = mna[33][1] - idt*Geq_M;
	mna[33](33) = mna[33][33] + idt*Geq_M;
	mna[68](33) = mna[68][33] +Geq_M;
	rval[33] += (Ieq_M - Qo_M)*idt;
	fvec[33] += -(Sr_M);
	mna[68](68) = mna[68][68] + 1.;
	rval[68] += Ieq_M;
	fvec[68] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM28*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[35];
	Un_M += Vx[1];
	Uo_M += -uxtold[35];
	Uo_M += uxtold[1];
	Qo_M = uxtold[69];
	Qn_M = Vx[69];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[69](1) = mna[69][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](35) = mna[1][35] - idt*Geq_M;
	mna[35](1) = mna[35][1] - idt*Geq_M;
	mna[35](35) = mna[35][35] + idt*Geq_M;
	mna[69](35) = mna[69][35] +Geq_M;
	rval[35] += (Ieq_M - Qo_M)*idt;
	fvec[35] += -(Sr_M);
	mna[69](69) = mna[69][69] + 1.;
	rval[69] += Ieq_M;
	fvec[69] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM29*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[37];
	Un_M += Vx[1];
	Uo_M += -uxtold[37];
	Uo_M += uxtold[1];
	Qo_M = uxtold[70];
	Qn_M = Vx[70];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[70](1) = mna[70][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](37) = mna[1][37] - idt*Geq_M;
	mna[37](1) = mna[37][1] - idt*Geq_M;
	mna[37](37) = mna[37][37] + idt*Geq_M;
	mna[70](37) = mna[70][37] +Geq_M;
	rval[37] += (Ieq_M - Qo_M)*idt;
	fvec[37] += -(Sr_M);
	mna[70](70) = mna[70][70] + 1.;
	rval[70] += Ieq_M;
	fvec[70] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM30*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[38];
	Un_M += Vx[1];
	Uo_M += -uxtold[38];
	Uo_M += uxtold[1];
	Qo_M = uxtold[71];
	Qn_M = Vx[71];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[71](1) = mna[71][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](38) = mna[1][38] - idt*Geq_M;
	mna[38](1) = mna[38][1] - idt*Geq_M;
	mna[38](38) = mna[38][38] + idt*Geq_M;
	mna[71](38) = mna[71][38] +Geq_M;
	rval[38] += (Ieq_M - Qo_M)*idt;
	fvec[38] += -(Sr_M);
	mna[71](71) = mna[71][71] + 1.;
	rval[71] += Ieq_M;
	fvec[71] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM31*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[40];
	Un_M += Vx[1];
	Uo_M += -uxtold[40];
	Uo_M += uxtold[1];
	Qo_M = uxtold[72];
	Qn_M = Vx[72];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[1](1) = mna[1][1] + idt*Geq_M;
	mna[72](1) = mna[72][1] -Geq_M;
	rval[1] += -(Ieq_M - Qo_M)*idt;
	fvec[1] += Sr_M;
	mna[1](40) = mna[1][40] - idt*Geq_M;
	mna[40](1) = mna[40][1] - idt*Geq_M;
	mna[40](40) = mna[40][40] + idt*Geq_M;
	mna[72](40) = mna[72][40] +Geq_M;
	rval[40] += (Ieq_M - Qo_M)*idt;
	fvec[40] += -(Sr_M);
	mna[72](72) = mna[72][72] + 1.;
	rval[72] += Ieq_M;
	fvec[72] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM32*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[3];
	Un_M += Vx[2];
	Uo_M += -uxtold[3];
	Uo_M += uxtold[2];
	Qo_M = uxtold[73];
	Qn_M = Vx[73];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[73](2) = mna[73][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](3) = mna[2][3] - idt*Geq_M;
	mna[3](2) = mna[3][2] - idt*Geq_M;
	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[73](3) = mna[73][3] +Geq_M;
	rval[3] += (Ieq_M - Qo_M)*idt;
	fvec[3] += -(Sr_M);
	mna[73](73) = mna[73][73] + 1.;
	rval[73] += Ieq_M;
	fvec[73] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM33*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[4];
	Un_M += Vx[2];
	Uo_M += -uxtold[4];
	Uo_M += uxtold[2];
	Qo_M = uxtold[74];
	Qn_M = Vx[74];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[74](2) = mna[74][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](4) = mna[2][4] - idt*Geq_M;
	mna[4](2) = mna[4][2] - idt*Geq_M;
	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[74](4) = mna[74][4] +Geq_M;
	rval[4] += (Ieq_M - Qo_M)*idt;
	fvec[4] += -(Sr_M);
	mna[74](74) = mna[74][74] + 1.;
	rval[74] += Ieq_M;
	fvec[74] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM34*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[5];
	Un_M += Vx[2];
	Uo_M += -uxtold[5];
	Uo_M += uxtold[2];
	Qo_M = uxtold[75];
	Qn_M = Vx[75];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[75](2) = mna[75][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](5) = mna[2][5] - idt*Geq_M;
	mna[5](2) = mna[5][2] - idt*Geq_M;
	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[75](5) = mna[75][5] +Geq_M;
	rval[5] += (Ieq_M - Qo_M)*idt;
	fvec[5] += -(Sr_M);
	mna[75](75) = mna[75][75] + 1.;
	rval[75] += Ieq_M;
	fvec[75] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM35*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[6];
	Un_M += Vx[2];
	Uo_M += -uxtold[6];
	Uo_M += uxtold[2];
	Qo_M = uxtold[76];
	Qn_M = Vx[76];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[76](2) = mna[76][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](6) = mna[2][6] - idt*Geq_M;
	mna[6](2) = mna[6][2] - idt*Geq_M;
	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[76](6) = mna[76][6] +Geq_M;
	rval[6] += (Ieq_M - Qo_M)*idt;
	fvec[6] += -(Sr_M);
	mna[76](76) = mna[76][76] + 1.;
	rval[76] += Ieq_M;
	fvec[76] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM36*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[7];
	Un_M += Vx[2];
	Uo_M += -uxtold[7];
	Uo_M += uxtold[2];
	Qo_M = uxtold[77];
	Qn_M = Vx[77];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[77](2) = mna[77][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](7) = mna[2][7] - idt*Geq_M;
	mna[7](2) = mna[7][2] - idt*Geq_M;
	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[77](7) = mna[77][7] +Geq_M;
	rval[7] += (Ieq_M - Qo_M)*idt;
	fvec[7] += -(Sr_M);
	mna[77](77) = mna[77][77] + 1.;
	rval[77] += Ieq_M;
	fvec[77] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM37*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[8];
	Un_M += Vx[2];
	Uo_M += -uxtold[8];
	Uo_M += uxtold[2];
	Qo_M = uxtold[78];
	Qn_M = Vx[78];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[78](2) = mna[78][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](8) = mna[2][8] - idt*Geq_M;
	mna[8](2) = mna[8][2] - idt*Geq_M;
	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[78](8) = mna[78][8] +Geq_M;
	rval[8] += (Ieq_M - Qo_M)*idt;
	fvec[8] += -(Sr_M);
	mna[78](78) = mna[78][78] + 1.;
	rval[78] += Ieq_M;
	fvec[78] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM38*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[9];
	Un_M += Vx[2];
	Uo_M += -uxtold[9];
	Uo_M += uxtold[2];
	Qo_M = uxtold[79];
	Qn_M = Vx[79];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[79](2) = mna[79][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](9) = mna[2][9] - idt*Geq_M;
	mna[9](2) = mna[9][2] - idt*Geq_M;
	mna[9](9) = mna[9][9] + idt*Geq_M;
	mna[79](9) = mna[79][9] +Geq_M;
	rval[9] += (Ieq_M - Qo_M)*idt;
	fvec[9] += -(Sr_M);
	mna[79](79) = mna[79][79] + 1.;
	rval[79] += Ieq_M;
	fvec[79] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM39*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[10];
	Un_M += Vx[2];
	Uo_M += -uxtold[10];
	Uo_M += uxtold[2];
	Qo_M = uxtold[80];
	Qn_M = Vx[80];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[80](2) = mna[80][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](10) = mna[2][10] - idt*Geq_M;
	mna[10](2) = mna[10][2] - idt*Geq_M;
	mna[10](10) = mna[10][10] + idt*Geq_M;
	mna[80](10) = mna[80][10] +Geq_M;
	rval[10] += (Ieq_M - Qo_M)*idt;
	fvec[10] += -(Sr_M);
	mna[80](80) = mna[80][80] + 1.;
	rval[80] += Ieq_M;
	fvec[80] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM40*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[11];
	Un_M += Vx[2];
	Uo_M += -uxtold[11];
	Uo_M += uxtold[2];
	Qo_M = uxtold[81];
	Qn_M = Vx[81];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[81](2) = mna[81][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](11) = mna[2][11] - idt*Geq_M;
	mna[11](2) = mna[11][2] - idt*Geq_M;
	mna[11](11) = mna[11][11] + idt*Geq_M;
	mna[81](11) = mna[81][11] +Geq_M;
	rval[11] += (Ieq_M - Qo_M)*idt;
	fvec[11] += -(Sr_M);
	mna[81](81) = mna[81][81] + 1.;
	rval[81] += Ieq_M;
	fvec[81] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM41*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[12];
	Un_M += Vx[2];
	Uo_M += -uxtold[12];
	Uo_M += uxtold[2];
	Qo_M = uxtold[82];
	Qn_M = Vx[82];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[82](2) = mna[82][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](12) = mna[2][12] - idt*Geq_M;
	mna[12](2) = mna[12][2] - idt*Geq_M;
	mna[12](12) = mna[12][12] + idt*Geq_M;
	mna[82](12) = mna[82][12] +Geq_M;
	rval[12] += (Ieq_M - Qo_M)*idt;
	fvec[12] += -(Sr_M);
	mna[82](82) = mna[82][82] + 1.;
	rval[82] += Ieq_M;
	fvec[82] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM42*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[13];
	Un_M += Vx[2];
	Uo_M += -uxtold[13];
	Uo_M += uxtold[2];
	Qo_M = uxtold[83];
	Qn_M = Vx[83];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[83](2) = mna[83][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](13) = mna[2][13] - idt*Geq_M;
	mna[13](2) = mna[13][2] - idt*Geq_M;
	mna[13](13) = mna[13][13] + idt*Geq_M;
	mna[83](13) = mna[83][13] +Geq_M;
	rval[13] += (Ieq_M - Qo_M)*idt;
	fvec[13] += -(Sr_M);
	mna[83](83) = mna[83][83] + 1.;
	rval[83] += Ieq_M;
	fvec[83] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM43*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[14];
	Un_M += Vx[2];
	Uo_M += -uxtold[14];
	Uo_M += uxtold[2];
	Qo_M = uxtold[84];
	Qn_M = Vx[84];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[84](2) = mna[84][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](14) = mna[2][14] - idt*Geq_M;
	mna[14](2) = mna[14][2] - idt*Geq_M;
	mna[14](14) = mna[14][14] + idt*Geq_M;
	mna[84](14) = mna[84][14] +Geq_M;
	rval[14] += (Ieq_M - Qo_M)*idt;
	fvec[14] += -(Sr_M);
	mna[84](84) = mna[84][84] + 1.;
	rval[84] += Ieq_M;
	fvec[84] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM44*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[16];
	Un_M += Vx[2];
	Uo_M += -uxtold[16];
	Uo_M += uxtold[2];
	Qo_M = uxtold[85];
	Qn_M = Vx[85];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[85](2) = mna[85][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](16) = mna[2][16] - idt*Geq_M;
	mna[16](2) = mna[16][2] - idt*Geq_M;
	mna[16](16) = mna[16][16] + idt*Geq_M;
	mna[85](16) = mna[85][16] +Geq_M;
	rval[16] += (Ieq_M - Qo_M)*idt;
	fvec[16] += -(Sr_M);
	mna[85](85) = mna[85][85] + 1.;
	rval[85] += Ieq_M;
	fvec[85] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM45*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[17];
	Un_M += Vx[2];
	Uo_M += -uxtold[17];
	Uo_M += uxtold[2];
	Qo_M = uxtold[86];
	Qn_M = Vx[86];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[86](2) = mna[86][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](17) = mna[2][17] - idt*Geq_M;
	mna[17](2) = mna[17][2] - idt*Geq_M;
	mna[17](17) = mna[17][17] + idt*Geq_M;
	mna[86](17) = mna[86][17] +Geq_M;
	rval[17] += (Ieq_M - Qo_M)*idt;
	fvec[17] += -(Sr_M);
	mna[86](86) = mna[86][86] + 1.;
	rval[86] += Ieq_M;
	fvec[86] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM46*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[18];
	Un_M += Vx[2];
	Uo_M += -uxtold[18];
	Uo_M += uxtold[2];
	Qo_M = uxtold[87];
	Qn_M = Vx[87];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[87](2) = mna[87][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](18) = mna[2][18] - idt*Geq_M;
	mna[18](2) = mna[18][2] - idt*Geq_M;
	mna[18](18) = mna[18][18] + idt*Geq_M;
	mna[87](18) = mna[87][18] +Geq_M;
	rval[18] += (Ieq_M - Qo_M)*idt;
	fvec[18] += -(Sr_M);
	mna[87](87) = mna[87][87] + 1.;
	rval[87] += Ieq_M;
	fvec[87] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM47*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[19];
	Un_M += Vx[2];
	Uo_M += -uxtold[19];
	Uo_M += uxtold[2];
	Qo_M = uxtold[88];
	Qn_M = Vx[88];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[88](2) = mna[88][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](19) = mna[2][19] - idt*Geq_M;
	mna[19](2) = mna[19][2] - idt*Geq_M;
	mna[19](19) = mna[19][19] + idt*Geq_M;
	mna[88](19) = mna[88][19] +Geq_M;
	rval[19] += (Ieq_M - Qo_M)*idt;
	fvec[19] += -(Sr_M);
	mna[88](88) = mna[88][88] + 1.;
	rval[88] += Ieq_M;
	fvec[88] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM48*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[20];
	Un_M += Vx[2];
	Uo_M += -uxtold[20];
	Uo_M += uxtold[2];
	Qo_M = uxtold[89];
	Qn_M = Vx[89];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[89](2) = mna[89][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](20) = mna[2][20] - idt*Geq_M;
	mna[20](2) = mna[20][2] - idt*Geq_M;
	mna[20](20) = mna[20][20] + idt*Geq_M;
	mna[89](20) = mna[89][20] +Geq_M;
	rval[20] += (Ieq_M - Qo_M)*idt;
	fvec[20] += -(Sr_M);
	mna[89](89) = mna[89][89] + 1.;
	rval[89] += Ieq_M;
	fvec[89] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM49*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[21];
	Un_M += Vx[2];
	Uo_M += -uxtold[21];
	Uo_M += uxtold[2];
	Qo_M = uxtold[90];
	Qn_M = Vx[90];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[90](2) = mna[90][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](21) = mna[2][21] - idt*Geq_M;
	mna[21](2) = mna[21][2] - idt*Geq_M;
	mna[21](21) = mna[21][21] + idt*Geq_M;
	mna[90](21) = mna[90][21] +Geq_M;
	rval[21] += (Ieq_M - Qo_M)*idt;
	fvec[21] += -(Sr_M);
	mna[90](90) = mna[90][90] + 1.;
	rval[90] += Ieq_M;
	fvec[90] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM50*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[23];
	Un_M += Vx[2];
	Uo_M += -uxtold[23];
	Uo_M += uxtold[2];
	Qo_M = uxtold[91];
	Qn_M = Vx[91];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[91](2) = mna[91][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](23) = mna[2][23] - idt*Geq_M;
	mna[23](2) = mna[23][2] - idt*Geq_M;
	mna[23](23) = mna[23][23] + idt*Geq_M;
	mna[91](23) = mna[91][23] +Geq_M;
	rval[23] += (Ieq_M - Qo_M)*idt;
	fvec[23] += -(Sr_M);
	mna[91](91) = mna[91][91] + 1.;
	rval[91] += Ieq_M;
	fvec[91] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM51*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[24];
	Un_M += Vx[2];
	Uo_M += -uxtold[24];
	Uo_M += uxtold[2];
	Qo_M = uxtold[92];
	Qn_M = Vx[92];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[92](2) = mna[92][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](24) = mna[2][24] - idt*Geq_M;
	mna[24](2) = mna[24][2] - idt*Geq_M;
	mna[24](24) = mna[24][24] + idt*Geq_M;
	mna[92](24) = mna[92][24] +Geq_M;
	rval[24] += (Ieq_M - Qo_M)*idt;
	fvec[24] += -(Sr_M);
	mna[92](92) = mna[92][92] + 1.;
	rval[92] += Ieq_M;
	fvec[92] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM52*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[26];
	Un_M += Vx[2];
	Uo_M += -uxtold[26];
	Uo_M += uxtold[2];
	Qo_M = uxtold[93];
	Qn_M = Vx[93];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[93](2) = mna[93][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](26) = mna[2][26] - idt*Geq_M;
	mna[26](2) = mna[26][2] - idt*Geq_M;
	mna[26](26) = mna[26][26] + idt*Geq_M;
	mna[93](26) = mna[93][26] +Geq_M;
	rval[26] += (Ieq_M - Qo_M)*idt;
	fvec[26] += -(Sr_M);
	mna[93](93) = mna[93][93] + 1.;
	rval[93] += Ieq_M;
	fvec[93] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM53*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[30];
	Un_M += Vx[2];
	Uo_M += -uxtold[30];
	Uo_M += uxtold[2];
	Qo_M = uxtold[94];
	Qn_M = Vx[94];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[94](2) = mna[94][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](30) = mna[2][30] - idt*Geq_M;
	mna[30](2) = mna[30][2] - idt*Geq_M;
	mna[30](30) = mna[30][30] + idt*Geq_M;
	mna[94](30) = mna[94][30] +Geq_M;
	rval[30] += (Ieq_M - Qo_M)*idt;
	fvec[30] += -(Sr_M);
	mna[94](94) = mna[94][94] + 1.;
	rval[94] += Ieq_M;
	fvec[94] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM54*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[31];
	Un_M += Vx[2];
	Uo_M += -uxtold[31];
	Uo_M += uxtold[2];
	Qo_M = uxtold[95];
	Qn_M = Vx[95];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[95](2) = mna[95][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](31) = mna[2][31] - idt*Geq_M;
	mna[31](2) = mna[31][2] - idt*Geq_M;
	mna[31](31) = mna[31][31] + idt*Geq_M;
	mna[95](31) = mna[95][31] +Geq_M;
	rval[31] += (Ieq_M - Qo_M)*idt;
	fvec[31] += -(Sr_M);
	mna[95](95) = mna[95][95] + 1.;
	rval[95] += Ieq_M;
	fvec[95] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM55*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[32];
	Un_M += Vx[2];
	Uo_M += -uxtold[32];
	Uo_M += uxtold[2];
	Qo_M = uxtold[96];
	Qn_M = Vx[96];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[96](2) = mna[96][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](32) = mna[2][32] - idt*Geq_M;
	mna[32](2) = mna[32][2] - idt*Geq_M;
	mna[32](32) = mna[32][32] + idt*Geq_M;
	mna[96](32) = mna[96][32] +Geq_M;
	rval[32] += (Ieq_M - Qo_M)*idt;
	fvec[32] += -(Sr_M);
	mna[96](96) = mna[96][96] + 1.;
	rval[96] += Ieq_M;
	fvec[96] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM56*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[34];
	Un_M += Vx[2];
	Uo_M += -uxtold[34];
	Uo_M += uxtold[2];
	Qo_M = uxtold[97];
	Qn_M = Vx[97];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[97](2) = mna[97][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](34) = mna[2][34] - idt*Geq_M;
	mna[34](2) = mna[34][2] - idt*Geq_M;
	mna[34](34) = mna[34][34] + idt*Geq_M;
	mna[97](34) = mna[97][34] +Geq_M;
	rval[34] += (Ieq_M - Qo_M)*idt;
	fvec[34] += -(Sr_M);
	mna[97](97) = mna[97][97] + 1.;
	rval[97] += Ieq_M;
	fvec[97] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM57*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[36];
	Un_M += Vx[2];
	Uo_M += -uxtold[36];
	Uo_M += uxtold[2];
	Qo_M = uxtold[98];
	Qn_M = Vx[98];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[98](2) = mna[98][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](36) = mna[2][36] - idt*Geq_M;
	mna[36](2) = mna[36][2] - idt*Geq_M;
	mna[36](36) = mna[36][36] + idt*Geq_M;
	mna[98](36) = mna[98][36] +Geq_M;
	rval[36] += (Ieq_M - Qo_M)*idt;
	fvec[36] += -(Sr_M);
	mna[98](98) = mna[98][98] + 1.;
	rval[98] += Ieq_M;
	fvec[98] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM58*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[37];
	Un_M += Vx[2];
	Uo_M += -uxtold[37];
	Uo_M += uxtold[2];
	Qo_M = uxtold[99];
	Qn_M = Vx[99];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[99](2) = mna[99][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](37) = mna[2][37] - idt*Geq_M;
	mna[37](2) = mna[37][2] - idt*Geq_M;
	mna[37](37) = mna[37][37] + idt*Geq_M;
	mna[99](37) = mna[99][37] +Geq_M;
	rval[37] += (Ieq_M - Qo_M)*idt;
	fvec[37] += -(Sr_M);
	mna[99](99) = mna[99][99] + 1.;
	rval[99] += Ieq_M;
	fvec[99] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM59*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[38];
	Un_M += Vx[2];
	Uo_M += -uxtold[38];
	Uo_M += uxtold[2];
	Qo_M = uxtold[100];
	Qn_M = Vx[100];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[100](2) = mna[100][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](38) = mna[2][38] - idt*Geq_M;
	mna[38](2) = mna[38][2] - idt*Geq_M;
	mna[38](38) = mna[38][38] + idt*Geq_M;
	mna[100](38) = mna[100][38] +Geq_M;
	rval[38] += (Ieq_M - Qo_M)*idt;
	fvec[38] += -(Sr_M);
	mna[100](100) = mna[100][100] + 1.;
	rval[100] += Ieq_M;
	fvec[100] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM60*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[40];
	Un_M += Vx[2];
	Uo_M += -uxtold[40];
	Uo_M += uxtold[2];
	Qo_M = uxtold[101];
	Qn_M = Vx[101];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[2](2) = mna[2][2] + idt*Geq_M;
	mna[101](2) = mna[101][2] -Geq_M;
	rval[2] += -(Ieq_M - Qo_M)*idt;
	fvec[2] += Sr_M;
	mna[2](40) = mna[2][40] - idt*Geq_M;
	mna[40](2) = mna[40][2] - idt*Geq_M;
	mna[40](40) = mna[40][40] + idt*Geq_M;
	mna[101](40) = mna[101][40] +Geq_M;
	rval[40] += (Ieq_M - Qo_M)*idt;
	fvec[40] += -(Sr_M);
	mna[101](101) = mna[101][101] + 1.;
	rval[101] += Ieq_M;
	fvec[101] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM61*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[4];
	Un_M += Vx[3];
	Uo_M += -uxtold[4];
	Uo_M += uxtold[3];
	Qo_M = uxtold[102];
	Qn_M = Vx[102];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[102](3) = mna[102][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](4) = mna[3][4] - idt*Geq_M;
	mna[4](3) = mna[4][3] - idt*Geq_M;
	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[102](4) = mna[102][4] +Geq_M;
	rval[4] += (Ieq_M - Qo_M)*idt;
	fvec[4] += -(Sr_M);
	mna[102](102) = mna[102][102] + 1.;
	rval[102] += Ieq_M;
	fvec[102] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM62*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[5];
	Un_M += Vx[3];
	Uo_M += -uxtold[5];
	Uo_M += uxtold[3];
	Qo_M = uxtold[103];
	Qn_M = Vx[103];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[103](3) = mna[103][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](5) = mna[3][5] - idt*Geq_M;
	mna[5](3) = mna[5][3] - idt*Geq_M;
	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[103](5) = mna[103][5] +Geq_M;
	rval[5] += (Ieq_M - Qo_M)*idt;
	fvec[5] += -(Sr_M);
	mna[103](103) = mna[103][103] + 1.;
	rval[103] += Ieq_M;
	fvec[103] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM63*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[6];
	Un_M += Vx[3];
	Uo_M += -uxtold[6];
	Uo_M += uxtold[3];
	Qo_M = uxtold[104];
	Qn_M = Vx[104];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[104](3) = mna[104][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](6) = mna[3][6] - idt*Geq_M;
	mna[6](3) = mna[6][3] - idt*Geq_M;
	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[104](6) = mna[104][6] +Geq_M;
	rval[6] += (Ieq_M - Qo_M)*idt;
	fvec[6] += -(Sr_M);
	mna[104](104) = mna[104][104] + 1.;
	rval[104] += Ieq_M;
	fvec[104] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM64*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[7];
	Un_M += Vx[3];
	Uo_M += -uxtold[7];
	Uo_M += uxtold[3];
	Qo_M = uxtold[105];
	Qn_M = Vx[105];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[105](3) = mna[105][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](7) = mna[3][7] - idt*Geq_M;
	mna[7](3) = mna[7][3] - idt*Geq_M;
	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[105](7) = mna[105][7] +Geq_M;
	rval[7] += (Ieq_M - Qo_M)*idt;
	fvec[7] += -(Sr_M);
	mna[105](105) = mna[105][105] + 1.;
	rval[105] += Ieq_M;
	fvec[105] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM65*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[8];
	Un_M += Vx[3];
	Uo_M += -uxtold[8];
	Uo_M += uxtold[3];
	Qo_M = uxtold[106];
	Qn_M = Vx[106];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[106](3) = mna[106][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](8) = mna[3][8] - idt*Geq_M;
	mna[8](3) = mna[8][3] - idt*Geq_M;
	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[106](8) = mna[106][8] +Geq_M;
	rval[8] += (Ieq_M - Qo_M)*idt;
	fvec[8] += -(Sr_M);
	mna[106](106) = mna[106][106] + 1.;
	rval[106] += Ieq_M;
	fvec[106] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM66*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[9];
	Un_M += Vx[3];
	Uo_M += -uxtold[9];
	Uo_M += uxtold[3];
	Qo_M = uxtold[107];
	Qn_M = Vx[107];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[107](3) = mna[107][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](9) = mna[3][9] - idt*Geq_M;
	mna[9](3) = mna[9][3] - idt*Geq_M;
	mna[9](9) = mna[9][9] + idt*Geq_M;
	mna[107](9) = mna[107][9] +Geq_M;
	rval[9] += (Ieq_M - Qo_M)*idt;
	fvec[9] += -(Sr_M);
	mna[107](107) = mna[107][107] + 1.;
	rval[107] += Ieq_M;
	fvec[107] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM67*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[10];
	Un_M += Vx[3];
	Uo_M += -uxtold[10];
	Uo_M += uxtold[3];
	Qo_M = uxtold[108];
	Qn_M = Vx[108];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[108](3) = mna[108][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](10) = mna[3][10] - idt*Geq_M;
	mna[10](3) = mna[10][3] - idt*Geq_M;
	mna[10](10) = mna[10][10] + idt*Geq_M;
	mna[108](10) = mna[108][10] +Geq_M;
	rval[10] += (Ieq_M - Qo_M)*idt;
	fvec[10] += -(Sr_M);
	mna[108](108) = mna[108][108] + 1.;
	rval[108] += Ieq_M;
	fvec[108] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM68*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[11];
	Un_M += Vx[3];
	Uo_M += -uxtold[11];
	Uo_M += uxtold[3];
	Qo_M = uxtold[109];
	Qn_M = Vx[109];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[109](3) = mna[109][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](11) = mna[3][11] - idt*Geq_M;
	mna[11](3) = mna[11][3] - idt*Geq_M;
	mna[11](11) = mna[11][11] + idt*Geq_M;
	mna[109](11) = mna[109][11] +Geq_M;
	rval[11] += (Ieq_M - Qo_M)*idt;
	fvec[11] += -(Sr_M);
	mna[109](109) = mna[109][109] + 1.;
	rval[109] += Ieq_M;
	fvec[109] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM69*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[12];
	Un_M += Vx[3];
	Uo_M += -uxtold[12];
	Uo_M += uxtold[3];
	Qo_M = uxtold[110];
	Qn_M = Vx[110];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[110](3) = mna[110][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](12) = mna[3][12] - idt*Geq_M;
	mna[12](3) = mna[12][3] - idt*Geq_M;
	mna[12](12) = mna[12][12] + idt*Geq_M;
	mna[110](12) = mna[110][12] +Geq_M;
	rval[12] += (Ieq_M - Qo_M)*idt;
	fvec[12] += -(Sr_M);
	mna[110](110) = mna[110][110] + 1.;
	rval[110] += Ieq_M;
	fvec[110] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM70*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[14];
	Un_M += Vx[3];
	Uo_M += -uxtold[14];
	Uo_M += uxtold[3];
	Qo_M = uxtold[111];
	Qn_M = Vx[111];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[111](3) = mna[111][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](14) = mna[3][14] - idt*Geq_M;
	mna[14](3) = mna[14][3] - idt*Geq_M;
	mna[14](14) = mna[14][14] + idt*Geq_M;
	mna[111](14) = mna[111][14] +Geq_M;
	rval[14] += (Ieq_M - Qo_M)*idt;
	fvec[14] += -(Sr_M);
	mna[111](111) = mna[111][111] + 1.;
	rval[111] += Ieq_M;
	fvec[111] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM71*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[17];
	Un_M += Vx[3];
	Uo_M += -uxtold[17];
	Uo_M += uxtold[3];
	Qo_M = uxtold[112];
	Qn_M = Vx[112];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[112](3) = mna[112][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](17) = mna[3][17] - idt*Geq_M;
	mna[17](3) = mna[17][3] - idt*Geq_M;
	mna[17](17) = mna[17][17] + idt*Geq_M;
	mna[112](17) = mna[112][17] +Geq_M;
	rval[17] += (Ieq_M - Qo_M)*idt;
	fvec[17] += -(Sr_M);
	mna[112](112) = mna[112][112] + 1.;
	rval[112] += Ieq_M;
	fvec[112] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM72*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[18];
	Un_M += Vx[3];
	Uo_M += -uxtold[18];
	Uo_M += uxtold[3];
	Qo_M = uxtold[113];
	Qn_M = Vx[113];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[113](3) = mna[113][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](18) = mna[3][18] - idt*Geq_M;
	mna[18](3) = mna[18][3] - idt*Geq_M;
	mna[18](18) = mna[18][18] + idt*Geq_M;
	mna[113](18) = mna[113][18] +Geq_M;
	rval[18] += (Ieq_M - Qo_M)*idt;
	fvec[18] += -(Sr_M);
	mna[113](113) = mna[113][113] + 1.;
	rval[113] += Ieq_M;
	fvec[113] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM73*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[19];
	Un_M += Vx[3];
	Uo_M += -uxtold[19];
	Uo_M += uxtold[3];
	Qo_M = uxtold[114];
	Qn_M = Vx[114];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[114](3) = mna[114][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](19) = mna[3][19] - idt*Geq_M;
	mna[19](3) = mna[19][3] - idt*Geq_M;
	mna[19](19) = mna[19][19] + idt*Geq_M;
	mna[114](19) = mna[114][19] +Geq_M;
	rval[19] += (Ieq_M - Qo_M)*idt;
	fvec[19] += -(Sr_M);
	mna[114](114) = mna[114][114] + 1.;
	rval[114] += Ieq_M;
	fvec[114] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM74*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[22];
	Un_M += Vx[3];
	Uo_M += -uxtold[22];
	Uo_M += uxtold[3];
	Qo_M = uxtold[115];
	Qn_M = Vx[115];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[115](3) = mna[115][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](22) = mna[3][22] - idt*Geq_M;
	mna[22](3) = mna[22][3] - idt*Geq_M;
	mna[22](22) = mna[22][22] + idt*Geq_M;
	mna[115](22) = mna[115][22] +Geq_M;
	rval[22] += (Ieq_M - Qo_M)*idt;
	fvec[22] += -(Sr_M);
	mna[115](115) = mna[115][115] + 1.;
	rval[115] += Ieq_M;
	fvec[115] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM75*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[26];
	Un_M += Vx[3];
	Uo_M += -uxtold[26];
	Uo_M += uxtold[3];
	Qo_M = uxtold[116];
	Qn_M = Vx[116];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[116](3) = mna[116][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](26) = mna[3][26] - idt*Geq_M;
	mna[26](3) = mna[26][3] - idt*Geq_M;
	mna[26](26) = mna[26][26] + idt*Geq_M;
	mna[116](26) = mna[116][26] +Geq_M;
	rval[26] += (Ieq_M - Qo_M)*idt;
	fvec[26] += -(Sr_M);
	mna[116](116) = mna[116][116] + 1.;
	rval[116] += Ieq_M;
	fvec[116] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM76*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[27];
	Un_M += Vx[3];
	Uo_M += -uxtold[27];
	Uo_M += uxtold[3];
	Qo_M = uxtold[117];
	Qn_M = Vx[117];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[117](3) = mna[117][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](27) = mna[3][27] - idt*Geq_M;
	mna[27](3) = mna[27][3] - idt*Geq_M;
	mna[27](27) = mna[27][27] + idt*Geq_M;
	mna[117](27) = mna[117][27] +Geq_M;
	rval[27] += (Ieq_M - Qo_M)*idt;
	fvec[27] += -(Sr_M);
	mna[117](117) = mna[117][117] + 1.;
	rval[117] += Ieq_M;
	fvec[117] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM77*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[28];
	Un_M += Vx[3];
	Uo_M += -uxtold[28];
	Uo_M += uxtold[3];
	Qo_M = uxtold[118];
	Qn_M = Vx[118];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[118](3) = mna[118][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](28) = mna[3][28] - idt*Geq_M;
	mna[28](3) = mna[28][3] - idt*Geq_M;
	mna[28](28) = mna[28][28] + idt*Geq_M;
	mna[118](28) = mna[118][28] +Geq_M;
	rval[28] += (Ieq_M - Qo_M)*idt;
	fvec[28] += -(Sr_M);
	mna[118](118) = mna[118][118] + 1.;
	rval[118] += Ieq_M;
	fvec[118] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM78*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[29];
	Un_M += Vx[3];
	Uo_M += -uxtold[29];
	Uo_M += uxtold[3];
	Qo_M = uxtold[119];
	Qn_M = Vx[119];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[119](3) = mna[119][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](29) = mna[3][29] - idt*Geq_M;
	mna[29](3) = mna[29][3] - idt*Geq_M;
	mna[29](29) = mna[29][29] + idt*Geq_M;
	mna[119](29) = mna[119][29] +Geq_M;
	rval[29] += (Ieq_M - Qo_M)*idt;
	fvec[29] += -(Sr_M);
	mna[119](119) = mna[119][119] + 1.;
	rval[119] += Ieq_M;
	fvec[119] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM79*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[30];
	Un_M += Vx[3];
	Uo_M += -uxtold[30];
	Uo_M += uxtold[3];
	Qo_M = uxtold[120];
	Qn_M = Vx[120];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[120](3) = mna[120][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](30) = mna[3][30] - idt*Geq_M;
	mna[30](3) = mna[30][3] - idt*Geq_M;
	mna[30](30) = mna[30][30] + idt*Geq_M;
	mna[120](30) = mna[120][30] +Geq_M;
	rval[30] += (Ieq_M - Qo_M)*idt;
	fvec[30] += -(Sr_M);
	mna[120](120) = mna[120][120] + 1.;
	rval[120] += Ieq_M;
	fvec[120] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM80*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[31];
	Un_M += Vx[3];
	Uo_M += -uxtold[31];
	Uo_M += uxtold[3];
	Qo_M = uxtold[121];
	Qn_M = Vx[121];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[121](3) = mna[121][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](31) = mna[3][31] - idt*Geq_M;
	mna[31](3) = mna[31][3] - idt*Geq_M;
	mna[31](31) = mna[31][31] + idt*Geq_M;
	mna[121](31) = mna[121][31] +Geq_M;
	rval[31] += (Ieq_M - Qo_M)*idt;
	fvec[31] += -(Sr_M);
	mna[121](121) = mna[121][121] + 1.;
	rval[121] += Ieq_M;
	fvec[121] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM81*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[32];
	Un_M += Vx[3];
	Uo_M += -uxtold[32];
	Uo_M += uxtold[3];
	Qo_M = uxtold[122];
	Qn_M = Vx[122];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[122](3) = mna[122][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](32) = mna[3][32] - idt*Geq_M;
	mna[32](3) = mna[32][3] - idt*Geq_M;
	mna[32](32) = mna[32][32] + idt*Geq_M;
	mna[122](32) = mna[122][32] +Geq_M;
	rval[32] += (Ieq_M - Qo_M)*idt;
	fvec[32] += -(Sr_M);
	mna[122](122) = mna[122][122] + 1.;
	rval[122] += Ieq_M;
	fvec[122] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM82*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[35];
	Un_M += Vx[3];
	Uo_M += -uxtold[35];
	Uo_M += uxtold[3];
	Qo_M = uxtold[123];
	Qn_M = Vx[123];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[123](3) = mna[123][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](35) = mna[3][35] - idt*Geq_M;
	mna[35](3) = mna[35][3] - idt*Geq_M;
	mna[35](35) = mna[35][35] + idt*Geq_M;
	mna[123](35) = mna[123][35] +Geq_M;
	rval[35] += (Ieq_M - Qo_M)*idt;
	fvec[35] += -(Sr_M);
	mna[123](123) = mna[123][123] + 1.;
	rval[123] += Ieq_M;
	fvec[123] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM83*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[37];
	Un_M += Vx[3];
	Uo_M += -uxtold[37];
	Uo_M += uxtold[3];
	Qo_M = uxtold[124];
	Qn_M = Vx[124];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[124](3) = mna[124][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](37) = mna[3][37] - idt*Geq_M;
	mna[37](3) = mna[37][3] - idt*Geq_M;
	mna[37](37) = mna[37][37] + idt*Geq_M;
	mna[124](37) = mna[124][37] +Geq_M;
	rval[37] += (Ieq_M - Qo_M)*idt;
	fvec[37] += -(Sr_M);
	mna[124](124) = mna[124][124] + 1.;
	rval[124] += Ieq_M;
	fvec[124] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM84*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[39];
	Un_M += Vx[3];
	Uo_M += -uxtold[39];
	Uo_M += uxtold[3];
	Qo_M = uxtold[125];
	Qn_M = Vx[125];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[3](3) = mna[3][3] + idt*Geq_M;
	mna[125](3) = mna[125][3] -Geq_M;
	rval[3] += -(Ieq_M - Qo_M)*idt;
	fvec[3] += Sr_M;
	mna[3](39) = mna[3][39] - idt*Geq_M;
	mna[39](3) = mna[39][3] - idt*Geq_M;
	mna[39](39) = mna[39][39] + idt*Geq_M;
	mna[125](39) = mna[125][39] +Geq_M;
	rval[39] += (Ieq_M - Qo_M)*idt;
	fvec[39] += -(Sr_M);
	mna[125](125) = mna[125][125] + 1.;
	rval[125] += Ieq_M;
	fvec[125] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM85*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[5];
	Un_M += Vx[4];
	Uo_M += -uxtold[5];
	Uo_M += uxtold[4];
	Qo_M = uxtold[126];
	Qn_M = Vx[126];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[126](4) = mna[126][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](5) = mna[4][5] - idt*Geq_M;
	mna[5](4) = mna[5][4] - idt*Geq_M;
	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[126](5) = mna[126][5] +Geq_M;
	rval[5] += (Ieq_M - Qo_M)*idt;
	fvec[5] += -(Sr_M);
	mna[126](126) = mna[126][126] + 1.;
	rval[126] += Ieq_M;
	fvec[126] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM86*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[6];
	Un_M += Vx[4];
	Uo_M += -uxtold[6];
	Uo_M += uxtold[4];
	Qo_M = uxtold[127];
	Qn_M = Vx[127];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[127](4) = mna[127][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](6) = mna[4][6] - idt*Geq_M;
	mna[6](4) = mna[6][4] - idt*Geq_M;
	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[127](6) = mna[127][6] +Geq_M;
	rval[6] += (Ieq_M - Qo_M)*idt;
	fvec[6] += -(Sr_M);
	mna[127](127) = mna[127][127] + 1.;
	rval[127] += Ieq_M;
	fvec[127] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM87*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[7];
	Un_M += Vx[4];
	Uo_M += -uxtold[7];
	Uo_M += uxtold[4];
	Qo_M = uxtold[128];
	Qn_M = Vx[128];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[128](4) = mna[128][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](7) = mna[4][7] - idt*Geq_M;
	mna[7](4) = mna[7][4] - idt*Geq_M;
	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[128](7) = mna[128][7] +Geq_M;
	rval[7] += (Ieq_M - Qo_M)*idt;
	fvec[7] += -(Sr_M);
	mna[128](128) = mna[128][128] + 1.;
	rval[128] += Ieq_M;
	fvec[128] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM88*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[8];
	Un_M += Vx[4];
	Uo_M += -uxtold[8];
	Uo_M += uxtold[4];
	Qo_M = uxtold[129];
	Qn_M = Vx[129];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[129](4) = mna[129][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](8) = mna[4][8] - idt*Geq_M;
	mna[8](4) = mna[8][4] - idt*Geq_M;
	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[129](8) = mna[129][8] +Geq_M;
	rval[8] += (Ieq_M - Qo_M)*idt;
	fvec[8] += -(Sr_M);
	mna[129](129) = mna[129][129] + 1.;
	rval[129] += Ieq_M;
	fvec[129] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM89*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[9];
	Un_M += Vx[4];
	Uo_M += -uxtold[9];
	Uo_M += uxtold[4];
	Qo_M = uxtold[130];
	Qn_M = Vx[130];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[130](4) = mna[130][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](9) = mna[4][9] - idt*Geq_M;
	mna[9](4) = mna[9][4] - idt*Geq_M;
	mna[9](9) = mna[9][9] + idt*Geq_M;
	mna[130](9) = mna[130][9] +Geq_M;
	rval[9] += (Ieq_M - Qo_M)*idt;
	fvec[9] += -(Sr_M);
	mna[130](130) = mna[130][130] + 1.;
	rval[130] += Ieq_M;
	fvec[130] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM90*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[12];
	Un_M += Vx[4];
	Uo_M += -uxtold[12];
	Uo_M += uxtold[4];
	Qo_M = uxtold[131];
	Qn_M = Vx[131];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[131](4) = mna[131][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](12) = mna[4][12] - idt*Geq_M;
	mna[12](4) = mna[12][4] - idt*Geq_M;
	mna[12](12) = mna[12][12] + idt*Geq_M;
	mna[131](12) = mna[131][12] +Geq_M;
	rval[12] += (Ieq_M - Qo_M)*idt;
	fvec[12] += -(Sr_M);
	mna[131](131) = mna[131][131] + 1.;
	rval[131] += Ieq_M;
	fvec[131] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM91*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[13];
	Un_M += Vx[4];
	Uo_M += -uxtold[13];
	Uo_M += uxtold[4];
	Qo_M = uxtold[132];
	Qn_M = Vx[132];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[132](4) = mna[132][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](13) = mna[4][13] - idt*Geq_M;
	mna[13](4) = mna[13][4] - idt*Geq_M;
	mna[13](13) = mna[13][13] + idt*Geq_M;
	mna[132](13) = mna[132][13] +Geq_M;
	rval[13] += (Ieq_M - Qo_M)*idt;
	fvec[13] += -(Sr_M);
	mna[132](132) = mna[132][132] + 1.;
	rval[132] += Ieq_M;
	fvec[132] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM92*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[14];
	Un_M += Vx[4];
	Uo_M += -uxtold[14];
	Uo_M += uxtold[4];
	Qo_M = uxtold[133];
	Qn_M = Vx[133];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[133](4) = mna[133][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](14) = mna[4][14] - idt*Geq_M;
	mna[14](4) = mna[14][4] - idt*Geq_M;
	mna[14](14) = mna[14][14] + idt*Geq_M;
	mna[133](14) = mna[133][14] +Geq_M;
	rval[14] += (Ieq_M - Qo_M)*idt;
	fvec[14] += -(Sr_M);
	mna[133](133) = mna[133][133] + 1.;
	rval[133] += Ieq_M;
	fvec[133] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM93*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[15];
	Un_M += Vx[4];
	Uo_M += -uxtold[15];
	Uo_M += uxtold[4];
	Qo_M = uxtold[134];
	Qn_M = Vx[134];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[134](4) = mna[134][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](15) = mna[4][15] - idt*Geq_M;
	mna[15](4) = mna[15][4] - idt*Geq_M;
	mna[15](15) = mna[15][15] + idt*Geq_M;
	mna[134](15) = mna[134][15] +Geq_M;
	rval[15] += (Ieq_M - Qo_M)*idt;
	fvec[15] += -(Sr_M);
	mna[134](134) = mna[134][134] + 1.;
	rval[134] += Ieq_M;
	fvec[134] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM94*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[16];
	Un_M += Vx[4];
	Uo_M += -uxtold[16];
	Uo_M += uxtold[4];
	Qo_M = uxtold[135];
	Qn_M = Vx[135];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[135](4) = mna[135][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](16) = mna[4][16] - idt*Geq_M;
	mna[16](4) = mna[16][4] - idt*Geq_M;
	mna[16](16) = mna[16][16] + idt*Geq_M;
	mna[135](16) = mna[135][16] +Geq_M;
	rval[16] += (Ieq_M - Qo_M)*idt;
	fvec[16] += -(Sr_M);
	mna[135](135) = mna[135][135] + 1.;
	rval[135] += Ieq_M;
	fvec[135] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM95*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[17];
	Un_M += Vx[4];
	Uo_M += -uxtold[17];
	Uo_M += uxtold[4];
	Qo_M = uxtold[136];
	Qn_M = Vx[136];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[136](4) = mna[136][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](17) = mna[4][17] - idt*Geq_M;
	mna[17](4) = mna[17][4] - idt*Geq_M;
	mna[17](17) = mna[17][17] + idt*Geq_M;
	mna[136](17) = mna[136][17] +Geq_M;
	rval[17] += (Ieq_M - Qo_M)*idt;
	fvec[17] += -(Sr_M);
	mna[136](136) = mna[136][136] + 1.;
	rval[136] += Ieq_M;
	fvec[136] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM96*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[18];
	Un_M += Vx[4];
	Uo_M += -uxtold[18];
	Uo_M += uxtold[4];
	Qo_M = uxtold[137];
	Qn_M = Vx[137];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[137](4) = mna[137][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](18) = mna[4][18] - idt*Geq_M;
	mna[18](4) = mna[18][4] - idt*Geq_M;
	mna[18](18) = mna[18][18] + idt*Geq_M;
	mna[137](18) = mna[137][18] +Geq_M;
	rval[18] += (Ieq_M - Qo_M)*idt;
	fvec[18] += -(Sr_M);
	mna[137](137) = mna[137][137] + 1.;
	rval[137] += Ieq_M;
	fvec[137] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM97*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[19];
	Un_M += Vx[4];
	Uo_M += -uxtold[19];
	Uo_M += uxtold[4];
	Qo_M = uxtold[138];
	Qn_M = Vx[138];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[138](4) = mna[138][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](19) = mna[4][19] - idt*Geq_M;
	mna[19](4) = mna[19][4] - idt*Geq_M;
	mna[19](19) = mna[19][19] + idt*Geq_M;
	mna[138](19) = mna[138][19] +Geq_M;
	rval[19] += (Ieq_M - Qo_M)*idt;
	fvec[19] += -(Sr_M);
	mna[138](138) = mna[138][138] + 1.;
	rval[138] += Ieq_M;
	fvec[138] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM98*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[20];
	Un_M += Vx[4];
	Uo_M += -uxtold[20];
	Uo_M += uxtold[4];
	Qo_M = uxtold[139];
	Qn_M = Vx[139];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[139](4) = mna[139][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](20) = mna[4][20] - idt*Geq_M;
	mna[20](4) = mna[20][4] - idt*Geq_M;
	mna[20](20) = mna[20][20] + idt*Geq_M;
	mna[139](20) = mna[139][20] +Geq_M;
	rval[20] += (Ieq_M - Qo_M)*idt;
	fvec[20] += -(Sr_M);
	mna[139](139) = mna[139][139] + 1.;
	rval[139] += Ieq_M;
	fvec[139] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM99*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[22];
	Un_M += Vx[4];
	Uo_M += -uxtold[22];
	Uo_M += uxtold[4];
	Qo_M = uxtold[140];
	Qn_M = Vx[140];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[140](4) = mna[140][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](22) = mna[4][22] - idt*Geq_M;
	mna[22](4) = mna[22][4] - idt*Geq_M;
	mna[22](22) = mna[22][22] + idt*Geq_M;
	mna[140](22) = mna[140][22] +Geq_M;
	rval[22] += (Ieq_M - Qo_M)*idt;
	fvec[22] += -(Sr_M);
	mna[140](140) = mna[140][140] + 1.;
	rval[140] += Ieq_M;
	fvec[140] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM100*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[23];
	Un_M += Vx[4];
	Uo_M += -uxtold[23];
	Uo_M += uxtold[4];
	Qo_M = uxtold[141];
	Qn_M = Vx[141];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[141](4) = mna[141][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](23) = mna[4][23] - idt*Geq_M;
	mna[23](4) = mna[23][4] - idt*Geq_M;
	mna[23](23) = mna[23][23] + idt*Geq_M;
	mna[141](23) = mna[141][23] +Geq_M;
	rval[23] += (Ieq_M - Qo_M)*idt;
	fvec[23] += -(Sr_M);
	mna[141](141) = mna[141][141] + 1.;
	rval[141] += Ieq_M;
	fvec[141] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM101*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[24];
	Un_M += Vx[4];
	Uo_M += -uxtold[24];
	Uo_M += uxtold[4];
	Qo_M = uxtold[142];
	Qn_M = Vx[142];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[142](4) = mna[142][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](24) = mna[4][24] - idt*Geq_M;
	mna[24](4) = mna[24][4] - idt*Geq_M;
	mna[24](24) = mna[24][24] + idt*Geq_M;
	mna[142](24) = mna[142][24] +Geq_M;
	rval[24] += (Ieq_M - Qo_M)*idt;
	fvec[24] += -(Sr_M);
	mna[142](142) = mna[142][142] + 1.;
	rval[142] += Ieq_M;
	fvec[142] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM102*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[25];
	Un_M += Vx[4];
	Uo_M += -uxtold[25];
	Uo_M += uxtold[4];
	Qo_M = uxtold[143];
	Qn_M = Vx[143];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[143](4) = mna[143][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](25) = mna[4][25] - idt*Geq_M;
	mna[25](4) = mna[25][4] - idt*Geq_M;
	mna[25](25) = mna[25][25] + idt*Geq_M;
	mna[143](25) = mna[143][25] +Geq_M;
	rval[25] += (Ieq_M - Qo_M)*idt;
	fvec[25] += -(Sr_M);
	mna[143](143) = mna[143][143] + 1.;
	rval[143] += Ieq_M;
	fvec[143] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM103*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[27];
	Un_M += Vx[4];
	Uo_M += -uxtold[27];
	Uo_M += uxtold[4];
	Qo_M = uxtold[144];
	Qn_M = Vx[144];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[144](4) = mna[144][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](27) = mna[4][27] - idt*Geq_M;
	mna[27](4) = mna[27][4] - idt*Geq_M;
	mna[27](27) = mna[27][27] + idt*Geq_M;
	mna[144](27) = mna[144][27] +Geq_M;
	rval[27] += (Ieq_M - Qo_M)*idt;
	fvec[27] += -(Sr_M);
	mna[144](144) = mna[144][144] + 1.;
	rval[144] += Ieq_M;
	fvec[144] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM104*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[28];
	Un_M += Vx[4];
	Uo_M += -uxtold[28];
	Uo_M += uxtold[4];
	Qo_M = uxtold[145];
	Qn_M = Vx[145];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[145](4) = mna[145][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](28) = mna[4][28] - idt*Geq_M;
	mna[28](4) = mna[28][4] - idt*Geq_M;
	mna[28](28) = mna[28][28] + idt*Geq_M;
	mna[145](28) = mna[145][28] +Geq_M;
	rval[28] += (Ieq_M - Qo_M)*idt;
	fvec[28] += -(Sr_M);
	mna[145](145) = mna[145][145] + 1.;
	rval[145] += Ieq_M;
	fvec[145] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM105*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[31];
	Un_M += Vx[4];
	Uo_M += -uxtold[31];
	Uo_M += uxtold[4];
	Qo_M = uxtold[146];
	Qn_M = Vx[146];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[146](4) = mna[146][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](31) = mna[4][31] - idt*Geq_M;
	mna[31](4) = mna[31][4] - idt*Geq_M;
	mna[31](31) = mna[31][31] + idt*Geq_M;
	mna[146](31) = mna[146][31] +Geq_M;
	rval[31] += (Ieq_M - Qo_M)*idt;
	fvec[31] += -(Sr_M);
	mna[146](146) = mna[146][146] + 1.;
	rval[146] += Ieq_M;
	fvec[146] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM106*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[32];
	Un_M += Vx[4];
	Uo_M += -uxtold[32];
	Uo_M += uxtold[4];
	Qo_M = uxtold[147];
	Qn_M = Vx[147];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[147](4) = mna[147][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](32) = mna[4][32] - idt*Geq_M;
	mna[32](4) = mna[32][4] - idt*Geq_M;
	mna[32](32) = mna[32][32] + idt*Geq_M;
	mna[147](32) = mna[147][32] +Geq_M;
	rval[32] += (Ieq_M - Qo_M)*idt;
	fvec[32] += -(Sr_M);
	mna[147](147) = mna[147][147] + 1.;
	rval[147] += Ieq_M;
	fvec[147] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM107*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[34];
	Un_M += Vx[4];
	Uo_M += -uxtold[34];
	Uo_M += uxtold[4];
	Qo_M = uxtold[148];
	Qn_M = Vx[148];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[148](4) = mna[148][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](34) = mna[4][34] - idt*Geq_M;
	mna[34](4) = mna[34][4] - idt*Geq_M;
	mna[34](34) = mna[34][34] + idt*Geq_M;
	mna[148](34) = mna[148][34] +Geq_M;
	rval[34] += (Ieq_M - Qo_M)*idt;
	fvec[34] += -(Sr_M);
	mna[148](148) = mna[148][148] + 1.;
	rval[148] += Ieq_M;
	fvec[148] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM108*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[35];
	Un_M += Vx[4];
	Uo_M += -uxtold[35];
	Uo_M += uxtold[4];
	Qo_M = uxtold[149];
	Qn_M = Vx[149];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[149](4) = mna[149][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](35) = mna[4][35] - idt*Geq_M;
	mna[35](4) = mna[35][4] - idt*Geq_M;
	mna[35](35) = mna[35][35] + idt*Geq_M;
	mna[149](35) = mna[149][35] +Geq_M;
	rval[35] += (Ieq_M - Qo_M)*idt;
	fvec[35] += -(Sr_M);
	mna[149](149) = mna[149][149] + 1.;
	rval[149] += Ieq_M;
	fvec[149] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM109*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[36];
	Un_M += Vx[4];
	Uo_M += -uxtold[36];
	Uo_M += uxtold[4];
	Qo_M = uxtold[150];
	Qn_M = Vx[150];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[150](4) = mna[150][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](36) = mna[4][36] - idt*Geq_M;
	mna[36](4) = mna[36][4] - idt*Geq_M;
	mna[36](36) = mna[36][36] + idt*Geq_M;
	mna[150](36) = mna[150][36] +Geq_M;
	rval[36] += (Ieq_M - Qo_M)*idt;
	fvec[36] += -(Sr_M);
	mna[150](150) = mna[150][150] + 1.;
	rval[150] += Ieq_M;
	fvec[150] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM110*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[37];
	Un_M += Vx[4];
	Uo_M += -uxtold[37];
	Uo_M += uxtold[4];
	Qo_M = uxtold[151];
	Qn_M = Vx[151];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[151](4) = mna[151][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](37) = mna[4][37] - idt*Geq_M;
	mna[37](4) = mna[37][4] - idt*Geq_M;
	mna[37](37) = mna[37][37] + idt*Geq_M;
	mna[151](37) = mna[151][37] +Geq_M;
	rval[37] += (Ieq_M - Qo_M)*idt;
	fvec[37] += -(Sr_M);
	mna[151](151) = mna[151][151] + 1.;
	rval[151] += Ieq_M;
	fvec[151] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM111*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[38];
	Un_M += Vx[4];
	Uo_M += -uxtold[38];
	Uo_M += uxtold[4];
	Qo_M = uxtold[152];
	Qn_M = Vx[152];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[152](4) = mna[152][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](38) = mna[4][38] - idt*Geq_M;
	mna[38](4) = mna[38][4] - idt*Geq_M;
	mna[38](38) = mna[38][38] + idt*Geq_M;
	mna[152](38) = mna[152][38] +Geq_M;
	rval[38] += (Ieq_M - Qo_M)*idt;
	fvec[38] += -(Sr_M);
	mna[152](152) = mna[152][152] + 1.;
	rval[152] += Ieq_M;
	fvec[152] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM112*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[39];
	Un_M += Vx[4];
	Uo_M += -uxtold[39];
	Uo_M += uxtold[4];
	Qo_M = uxtold[153];
	Qn_M = Vx[153];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[153](4) = mna[153][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](39) = mna[4][39] - idt*Geq_M;
	mna[39](4) = mna[39][4] - idt*Geq_M;
	mna[39](39) = mna[39][39] + idt*Geq_M;
	mna[153](39) = mna[153][39] +Geq_M;
	rval[39] += (Ieq_M - Qo_M)*idt;
	fvec[39] += -(Sr_M);
	mna[153](153) = mna[153][153] + 1.;
	rval[153] += Ieq_M;
	fvec[153] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM113*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[40];
	Un_M += Vx[4];
	Uo_M += -uxtold[40];
	Uo_M += uxtold[4];
	Qo_M = uxtold[154];
	Qn_M = Vx[154];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[4](4) = mna[4][4] + idt*Geq_M;
	mna[154](4) = mna[154][4] -Geq_M;
	rval[4] += -(Ieq_M - Qo_M)*idt;
	fvec[4] += Sr_M;
	mna[4](40) = mna[4][40] - idt*Geq_M;
	mna[40](4) = mna[40][4] - idt*Geq_M;
	mna[40](40) = mna[40][40] + idt*Geq_M;
	mna[154](40) = mna[154][40] +Geq_M;
	rval[40] += (Ieq_M - Qo_M)*idt;
	fvec[40] += -(Sr_M);
	mna[154](154) = mna[154][154] + 1.;
	rval[154] += Ieq_M;
	fvec[154] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM114*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[6];
	Un_M += Vx[5];
	Uo_M += -uxtold[6];
	Uo_M += uxtold[5];
	Qo_M = uxtold[155];
	Qn_M = Vx[155];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[155](5) = mna[155][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](6) = mna[5][6] - idt*Geq_M;
	mna[6](5) = mna[6][5] - idt*Geq_M;
	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[155](6) = mna[155][6] +Geq_M;
	rval[6] += (Ieq_M - Qo_M)*idt;
	fvec[6] += -(Sr_M);
	mna[155](155) = mna[155][155] + 1.;
	rval[155] += Ieq_M;
	fvec[155] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM115*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[7];
	Un_M += Vx[5];
	Uo_M += -uxtold[7];
	Uo_M += uxtold[5];
	Qo_M = uxtold[156];
	Qn_M = Vx[156];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[156](5) = mna[156][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](7) = mna[5][7] - idt*Geq_M;
	mna[7](5) = mna[7][5] - idt*Geq_M;
	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[156](7) = mna[156][7] +Geq_M;
	rval[7] += (Ieq_M - Qo_M)*idt;
	fvec[7] += -(Sr_M);
	mna[156](156) = mna[156][156] + 1.;
	rval[156] += Ieq_M;
	fvec[156] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM116*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[8];
	Un_M += Vx[5];
	Uo_M += -uxtold[8];
	Uo_M += uxtold[5];
	Qo_M = uxtold[157];
	Qn_M = Vx[157];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[157](5) = mna[157][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](8) = mna[5][8] - idt*Geq_M;
	mna[8](5) = mna[8][5] - idt*Geq_M;
	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[157](8) = mna[157][8] +Geq_M;
	rval[8] += (Ieq_M - Qo_M)*idt;
	fvec[8] += -(Sr_M);
	mna[157](157) = mna[157][157] + 1.;
	rval[157] += Ieq_M;
	fvec[157] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM117*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[11];
	Un_M += Vx[5];
	Uo_M += -uxtold[11];
	Uo_M += uxtold[5];
	Qo_M = uxtold[158];
	Qn_M = Vx[158];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[158](5) = mna[158][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](11) = mna[5][11] - idt*Geq_M;
	mna[11](5) = mna[11][5] - idt*Geq_M;
	mna[11](11) = mna[11][11] + idt*Geq_M;
	mna[158](11) = mna[158][11] +Geq_M;
	rval[11] += (Ieq_M - Qo_M)*idt;
	fvec[11] += -(Sr_M);
	mna[158](158) = mna[158][158] + 1.;
	rval[158] += Ieq_M;
	fvec[158] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM118*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[12];
	Un_M += Vx[5];
	Uo_M += -uxtold[12];
	Uo_M += uxtold[5];
	Qo_M = uxtold[159];
	Qn_M = Vx[159];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[159](5) = mna[159][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](12) = mna[5][12] - idt*Geq_M;
	mna[12](5) = mna[12][5] - idt*Geq_M;
	mna[12](12) = mna[12][12] + idt*Geq_M;
	mna[159](12) = mna[159][12] +Geq_M;
	rval[12] += (Ieq_M - Qo_M)*idt;
	fvec[12] += -(Sr_M);
	mna[159](159) = mna[159][159] + 1.;
	rval[159] += Ieq_M;
	fvec[159] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM119*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[13];
	Un_M += Vx[5];
	Uo_M += -uxtold[13];
	Uo_M += uxtold[5];
	Qo_M = uxtold[160];
	Qn_M = Vx[160];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[160](5) = mna[160][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](13) = mna[5][13] - idt*Geq_M;
	mna[13](5) = mna[13][5] - idt*Geq_M;
	mna[13](13) = mna[13][13] + idt*Geq_M;
	mna[160](13) = mna[160][13] +Geq_M;
	rval[13] += (Ieq_M - Qo_M)*idt;
	fvec[13] += -(Sr_M);
	mna[160](160) = mna[160][160] + 1.;
	rval[160] += Ieq_M;
	fvec[160] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM120*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[15];
	Un_M += Vx[5];
	Uo_M += -uxtold[15];
	Uo_M += uxtold[5];
	Qo_M = uxtold[161];
	Qn_M = Vx[161];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[161](5) = mna[161][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](15) = mna[5][15] - idt*Geq_M;
	mna[15](5) = mna[15][5] - idt*Geq_M;
	mna[15](15) = mna[15][15] + idt*Geq_M;
	mna[161](15) = mna[161][15] +Geq_M;
	rval[15] += (Ieq_M - Qo_M)*idt;
	fvec[15] += -(Sr_M);
	mna[161](161) = mna[161][161] + 1.;
	rval[161] += Ieq_M;
	fvec[161] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM121*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[16];
	Un_M += Vx[5];
	Uo_M += -uxtold[16];
	Uo_M += uxtold[5];
	Qo_M = uxtold[162];
	Qn_M = Vx[162];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[162](5) = mna[162][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](16) = mna[5][16] - idt*Geq_M;
	mna[16](5) = mna[16][5] - idt*Geq_M;
	mna[16](16) = mna[16][16] + idt*Geq_M;
	mna[162](16) = mna[162][16] +Geq_M;
	rval[16] += (Ieq_M - Qo_M)*idt;
	fvec[16] += -(Sr_M);
	mna[162](162) = mna[162][162] + 1.;
	rval[162] += Ieq_M;
	fvec[162] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM122*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[19];
	Un_M += Vx[5];
	Uo_M += -uxtold[19];
	Uo_M += uxtold[5];
	Qo_M = uxtold[163];
	Qn_M = Vx[163];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[163](5) = mna[163][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](19) = mna[5][19] - idt*Geq_M;
	mna[19](5) = mna[19][5] - idt*Geq_M;
	mna[19](19) = mna[19][19] + idt*Geq_M;
	mna[163](19) = mna[163][19] +Geq_M;
	rval[19] += (Ieq_M - Qo_M)*idt;
	fvec[19] += -(Sr_M);
	mna[163](163) = mna[163][163] + 1.;
	rval[163] += Ieq_M;
	fvec[163] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM123*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[21];
	Un_M += Vx[5];
	Uo_M += -uxtold[21];
	Uo_M += uxtold[5];
	Qo_M = uxtold[164];
	Qn_M = Vx[164];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[164](5) = mna[164][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](21) = mna[5][21] - idt*Geq_M;
	mna[21](5) = mna[21][5] - idt*Geq_M;
	mna[21](21) = mna[21][21] + idt*Geq_M;
	mna[164](21) = mna[164][21] +Geq_M;
	rval[21] += (Ieq_M - Qo_M)*idt;
	fvec[21] += -(Sr_M);
	mna[164](164) = mna[164][164] + 1.;
	rval[164] += Ieq_M;
	fvec[164] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM124*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[23];
	Un_M += Vx[5];
	Uo_M += -uxtold[23];
	Uo_M += uxtold[5];
	Qo_M = uxtold[165];
	Qn_M = Vx[165];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[165](5) = mna[165][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](23) = mna[5][23] - idt*Geq_M;
	mna[23](5) = mna[23][5] - idt*Geq_M;
	mna[23](23) = mna[23][23] + idt*Geq_M;
	mna[165](23) = mna[165][23] +Geq_M;
	rval[23] += (Ieq_M - Qo_M)*idt;
	fvec[23] += -(Sr_M);
	mna[165](165) = mna[165][165] + 1.;
	rval[165] += Ieq_M;
	fvec[165] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM125*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[24];
	Un_M += Vx[5];
	Uo_M += -uxtold[24];
	Uo_M += uxtold[5];
	Qo_M = uxtold[166];
	Qn_M = Vx[166];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[166](5) = mna[166][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](24) = mna[5][24] - idt*Geq_M;
	mna[24](5) = mna[24][5] - idt*Geq_M;
	mna[24](24) = mna[24][24] + idt*Geq_M;
	mna[166](24) = mna[166][24] +Geq_M;
	rval[24] += (Ieq_M - Qo_M)*idt;
	fvec[24] += -(Sr_M);
	mna[166](166) = mna[166][166] + 1.;
	rval[166] += Ieq_M;
	fvec[166] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM126*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[25];
	Un_M += Vx[5];
	Uo_M += -uxtold[25];
	Uo_M += uxtold[5];
	Qo_M = uxtold[167];
	Qn_M = Vx[167];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[167](5) = mna[167][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](25) = mna[5][25] - idt*Geq_M;
	mna[25](5) = mna[25][5] - idt*Geq_M;
	mna[25](25) = mna[25][25] + idt*Geq_M;
	mna[167](25) = mna[167][25] +Geq_M;
	rval[25] += (Ieq_M - Qo_M)*idt;
	fvec[25] += -(Sr_M);
	mna[167](167) = mna[167][167] + 1.;
	rval[167] += Ieq_M;
	fvec[167] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM127*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[26];
	Un_M += Vx[5];
	Uo_M += -uxtold[26];
	Uo_M += uxtold[5];
	Qo_M = uxtold[168];
	Qn_M = Vx[168];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[168](5) = mna[168][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](26) = mna[5][26] - idt*Geq_M;
	mna[26](5) = mna[26][5] - idt*Geq_M;
	mna[26](26) = mna[26][26] + idt*Geq_M;
	mna[168](26) = mna[168][26] +Geq_M;
	rval[26] += (Ieq_M - Qo_M)*idt;
	fvec[26] += -(Sr_M);
	mna[168](168) = mna[168][168] + 1.;
	rval[168] += Ieq_M;
	fvec[168] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM128*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[27];
	Un_M += Vx[5];
	Uo_M += -uxtold[27];
	Uo_M += uxtold[5];
	Qo_M = uxtold[169];
	Qn_M = Vx[169];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[169](5) = mna[169][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](27) = mna[5][27] - idt*Geq_M;
	mna[27](5) = mna[27][5] - idt*Geq_M;
	mna[27](27) = mna[27][27] + idt*Geq_M;
	mna[169](27) = mna[169][27] +Geq_M;
	rval[27] += (Ieq_M - Qo_M)*idt;
	fvec[27] += -(Sr_M);
	mna[169](169) = mna[169][169] + 1.;
	rval[169] += Ieq_M;
	fvec[169] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM129*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[29];
	Un_M += Vx[5];
	Uo_M += -uxtold[29];
	Uo_M += uxtold[5];
	Qo_M = uxtold[170];
	Qn_M = Vx[170];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[170](5) = mna[170][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](29) = mna[5][29] - idt*Geq_M;
	mna[29](5) = mna[29][5] - idt*Geq_M;
	mna[29](29) = mna[29][29] + idt*Geq_M;
	mna[170](29) = mna[170][29] +Geq_M;
	rval[29] += (Ieq_M - Qo_M)*idt;
	fvec[29] += -(Sr_M);
	mna[170](170) = mna[170][170] + 1.;
	rval[170] += Ieq_M;
	fvec[170] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM130*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[31];
	Un_M += Vx[5];
	Uo_M += -uxtold[31];
	Uo_M += uxtold[5];
	Qo_M = uxtold[171];
	Qn_M = Vx[171];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[171](5) = mna[171][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](31) = mna[5][31] - idt*Geq_M;
	mna[31](5) = mna[31][5] - idt*Geq_M;
	mna[31](31) = mna[31][31] + idt*Geq_M;
	mna[171](31) = mna[171][31] +Geq_M;
	rval[31] += (Ieq_M - Qo_M)*idt;
	fvec[31] += -(Sr_M);
	mna[171](171) = mna[171][171] + 1.;
	rval[171] += Ieq_M;
	fvec[171] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM131*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[32];
	Un_M += Vx[5];
	Uo_M += -uxtold[32];
	Uo_M += uxtold[5];
	Qo_M = uxtold[172];
	Qn_M = Vx[172];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[172](5) = mna[172][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](32) = mna[5][32] - idt*Geq_M;
	mna[32](5) = mna[32][5] - idt*Geq_M;
	mna[32](32) = mna[32][32] + idt*Geq_M;
	mna[172](32) = mna[172][32] +Geq_M;
	rval[32] += (Ieq_M - Qo_M)*idt;
	fvec[32] += -(Sr_M);
	mna[172](172) = mna[172][172] + 1.;
	rval[172] += Ieq_M;
	fvec[172] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM132*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[33];
	Un_M += Vx[5];
	Uo_M += -uxtold[33];
	Uo_M += uxtold[5];
	Qo_M = uxtold[173];
	Qn_M = Vx[173];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[173](5) = mna[173][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](33) = mna[5][33] - idt*Geq_M;
	mna[33](5) = mna[33][5] - idt*Geq_M;
	mna[33](33) = mna[33][33] + idt*Geq_M;
	mna[173](33) = mna[173][33] +Geq_M;
	rval[33] += (Ieq_M - Qo_M)*idt;
	fvec[33] += -(Sr_M);
	mna[173](173) = mna[173][173] + 1.;
	rval[173] += Ieq_M;
	fvec[173] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM133*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[34];
	Un_M += Vx[5];
	Uo_M += -uxtold[34];
	Uo_M += uxtold[5];
	Qo_M = uxtold[174];
	Qn_M = Vx[174];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[174](5) = mna[174][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](34) = mna[5][34] - idt*Geq_M;
	mna[34](5) = mna[34][5] - idt*Geq_M;
	mna[34](34) = mna[34][34] + idt*Geq_M;
	mna[174](34) = mna[174][34] +Geq_M;
	rval[34] += (Ieq_M - Qo_M)*idt;
	fvec[34] += -(Sr_M);
	mna[174](174) = mna[174][174] + 1.;
	rval[174] += Ieq_M;
	fvec[174] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM134*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[36];
	Un_M += Vx[5];
	Uo_M += -uxtold[36];
	Uo_M += uxtold[5];
	Qo_M = uxtold[175];
	Qn_M = Vx[175];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[175](5) = mna[175][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](36) = mna[5][36] - idt*Geq_M;
	mna[36](5) = mna[36][5] - idt*Geq_M;
	mna[36](36) = mna[36][36] + idt*Geq_M;
	mna[175](36) = mna[175][36] +Geq_M;
	rval[36] += (Ieq_M - Qo_M)*idt;
	fvec[36] += -(Sr_M);
	mna[175](175) = mna[175][175] + 1.;
	rval[175] += Ieq_M;
	fvec[175] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM135*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[37];
	Un_M += Vx[5];
	Uo_M += -uxtold[37];
	Uo_M += uxtold[5];
	Qo_M = uxtold[176];
	Qn_M = Vx[176];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[176](5) = mna[176][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](37) = mna[5][37] - idt*Geq_M;
	mna[37](5) = mna[37][5] - idt*Geq_M;
	mna[37](37) = mna[37][37] + idt*Geq_M;
	mna[176](37) = mna[176][37] +Geq_M;
	rval[37] += (Ieq_M - Qo_M)*idt;
	fvec[37] += -(Sr_M);
	mna[176](176) = mna[176][176] + 1.;
	rval[176] += Ieq_M;
	fvec[176] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM136*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[38];
	Un_M += Vx[5];
	Uo_M += -uxtold[38];
	Uo_M += uxtold[5];
	Qo_M = uxtold[177];
	Qn_M = Vx[177];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[177](5) = mna[177][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](38) = mna[5][38] - idt*Geq_M;
	mna[38](5) = mna[38][5] - idt*Geq_M;
	mna[38](38) = mna[38][38] + idt*Geq_M;
	mna[177](38) = mna[177][38] +Geq_M;
	rval[38] += (Ieq_M - Qo_M)*idt;
	fvec[38] += -(Sr_M);
	mna[177](177) = mna[177][177] + 1.;
	rval[177] += Ieq_M;
	fvec[177] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM137*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[40];
	Un_M += Vx[5];
	Uo_M += -uxtold[40];
	Uo_M += uxtold[5];
	Qo_M = uxtold[178];
	Qn_M = Vx[178];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[5](5) = mna[5][5] + idt*Geq_M;
	mna[178](5) = mna[178][5] -Geq_M;
	rval[5] += -(Ieq_M - Qo_M)*idt;
	fvec[5] += Sr_M;
	mna[5](40) = mna[5][40] - idt*Geq_M;
	mna[40](5) = mna[40][5] - idt*Geq_M;
	mna[40](40) = mna[40][40] + idt*Geq_M;
	mna[178](40) = mna[178][40] +Geq_M;
	rval[40] += (Ieq_M - Qo_M)*idt;
	fvec[40] += -(Sr_M);
	mna[178](178) = mna[178][178] + 1.;
	rval[178] += Ieq_M;
	fvec[178] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM138*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[7];
	Un_M += Vx[6];
	Uo_M += -uxtold[7];
	Uo_M += uxtold[6];
	Qo_M = uxtold[179];
	Qn_M = Vx[179];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[179](6) = mna[179][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](7) = mna[6][7] - idt*Geq_M;
	mna[7](6) = mna[7][6] - idt*Geq_M;
	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[179](7) = mna[179][7] +Geq_M;
	rval[7] += (Ieq_M - Qo_M)*idt;
	fvec[7] += -(Sr_M);
	mna[179](179) = mna[179][179] + 1.;
	rval[179] += Ieq_M;
	fvec[179] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM139*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[8];
	Un_M += Vx[6];
	Uo_M += -uxtold[8];
	Uo_M += uxtold[6];
	Qo_M = uxtold[180];
	Qn_M = Vx[180];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[180](6) = mna[180][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](8) = mna[6][8] - idt*Geq_M;
	mna[8](6) = mna[8][6] - idt*Geq_M;
	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[180](8) = mna[180][8] +Geq_M;
	rval[8] += (Ieq_M - Qo_M)*idt;
	fvec[8] += -(Sr_M);
	mna[180](180) = mna[180][180] + 1.;
	rval[180] += Ieq_M;
	fvec[180] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM140*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[10];
	Un_M += Vx[6];
	Uo_M += -uxtold[10];
	Uo_M += uxtold[6];
	Qo_M = uxtold[181];
	Qn_M = Vx[181];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[181](6) = mna[181][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](10) = mna[6][10] - idt*Geq_M;
	mna[10](6) = mna[10][6] - idt*Geq_M;
	mna[10](10) = mna[10][10] + idt*Geq_M;
	mna[181](10) = mna[181][10] +Geq_M;
	rval[10] += (Ieq_M - Qo_M)*idt;
	fvec[10] += -(Sr_M);
	mna[181](181) = mna[181][181] + 1.;
	rval[181] += Ieq_M;
	fvec[181] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM141*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[11];
	Un_M += Vx[6];
	Uo_M += -uxtold[11];
	Uo_M += uxtold[6];
	Qo_M = uxtold[182];
	Qn_M = Vx[182];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[182](6) = mna[182][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](11) = mna[6][11] - idt*Geq_M;
	mna[11](6) = mna[11][6] - idt*Geq_M;
	mna[11](11) = mna[11][11] + idt*Geq_M;
	mna[182](11) = mna[182][11] +Geq_M;
	rval[11] += (Ieq_M - Qo_M)*idt;
	fvec[11] += -(Sr_M);
	mna[182](182) = mna[182][182] + 1.;
	rval[182] += Ieq_M;
	fvec[182] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM142*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[12];
	Un_M += Vx[6];
	Uo_M += -uxtold[12];
	Uo_M += uxtold[6];
	Qo_M = uxtold[183];
	Qn_M = Vx[183];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[183](6) = mna[183][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](12) = mna[6][12] - idt*Geq_M;
	mna[12](6) = mna[12][6] - idt*Geq_M;
	mna[12](12) = mna[12][12] + idt*Geq_M;
	mna[183](12) = mna[183][12] +Geq_M;
	rval[12] += (Ieq_M - Qo_M)*idt;
	fvec[12] += -(Sr_M);
	mna[183](183) = mna[183][183] + 1.;
	rval[183] += Ieq_M;
	fvec[183] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM143*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[13];
	Un_M += Vx[6];
	Uo_M += -uxtold[13];
	Uo_M += uxtold[6];
	Qo_M = uxtold[184];
	Qn_M = Vx[184];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[184](6) = mna[184][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](13) = mna[6][13] - idt*Geq_M;
	mna[13](6) = mna[13][6] - idt*Geq_M;
	mna[13](13) = mna[13][13] + idt*Geq_M;
	mna[184](13) = mna[184][13] +Geq_M;
	rval[13] += (Ieq_M - Qo_M)*idt;
	fvec[13] += -(Sr_M);
	mna[184](184) = mna[184][184] + 1.;
	rval[184] += Ieq_M;
	fvec[184] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM144*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[14];
	Un_M += Vx[6];
	Uo_M += -uxtold[14];
	Uo_M += uxtold[6];
	Qo_M = uxtold[185];
	Qn_M = Vx[185];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[185](6) = mna[185][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](14) = mna[6][14] - idt*Geq_M;
	mna[14](6) = mna[14][6] - idt*Geq_M;
	mna[14](14) = mna[14][14] + idt*Geq_M;
	mna[185](14) = mna[185][14] +Geq_M;
	rval[14] += (Ieq_M - Qo_M)*idt;
	fvec[14] += -(Sr_M);
	mna[185](185) = mna[185][185] + 1.;
	rval[185] += Ieq_M;
	fvec[185] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM145*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[15];
	Un_M += Vx[6];
	Uo_M += -uxtold[15];
	Uo_M += uxtold[6];
	Qo_M = uxtold[186];
	Qn_M = Vx[186];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[186](6) = mna[186][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](15) = mna[6][15] - idt*Geq_M;
	mna[15](6) = mna[15][6] - idt*Geq_M;
	mna[15](15) = mna[15][15] + idt*Geq_M;
	mna[186](15) = mna[186][15] +Geq_M;
	rval[15] += (Ieq_M - Qo_M)*idt;
	fvec[15] += -(Sr_M);
	mna[186](186) = mna[186][186] + 1.;
	rval[186] += Ieq_M;
	fvec[186] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM146*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[18];
	Un_M += Vx[6];
	Uo_M += -uxtold[18];
	Uo_M += uxtold[6];
	Qo_M = uxtold[187];
	Qn_M = Vx[187];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[187](6) = mna[187][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](18) = mna[6][18] - idt*Geq_M;
	mna[18](6) = mna[18][6] - idt*Geq_M;
	mna[18](18) = mna[18][18] + idt*Geq_M;
	mna[187](18) = mna[187][18] +Geq_M;
	rval[18] += (Ieq_M - Qo_M)*idt;
	fvec[18] += -(Sr_M);
	mna[187](187) = mna[187][187] + 1.;
	rval[187] += Ieq_M;
	fvec[187] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM147*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[19];
	Un_M += Vx[6];
	Uo_M += -uxtold[19];
	Uo_M += uxtold[6];
	Qo_M = uxtold[188];
	Qn_M = Vx[188];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[188](6) = mna[188][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](19) = mna[6][19] - idt*Geq_M;
	mna[19](6) = mna[19][6] - idt*Geq_M;
	mna[19](19) = mna[19][19] + idt*Geq_M;
	mna[188](19) = mna[188][19] +Geq_M;
	rval[19] += (Ieq_M - Qo_M)*idt;
	fvec[19] += -(Sr_M);
	mna[188](188) = mna[188][188] + 1.;
	rval[188] += Ieq_M;
	fvec[188] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM148*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[20];
	Un_M += Vx[6];
	Uo_M += -uxtold[20];
	Uo_M += uxtold[6];
	Qo_M = uxtold[189];
	Qn_M = Vx[189];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[189](6) = mna[189][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](20) = mna[6][20] - idt*Geq_M;
	mna[20](6) = mna[20][6] - idt*Geq_M;
	mna[20](20) = mna[20][20] + idt*Geq_M;
	mna[189](20) = mna[189][20] +Geq_M;
	rval[20] += (Ieq_M - Qo_M)*idt;
	fvec[20] += -(Sr_M);
	mna[189](189) = mna[189][189] + 1.;
	rval[189] += Ieq_M;
	fvec[189] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM149*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[22];
	Un_M += Vx[6];
	Uo_M += -uxtold[22];
	Uo_M += uxtold[6];
	Qo_M = uxtold[190];
	Qn_M = Vx[190];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[190](6) = mna[190][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](22) = mna[6][22] - idt*Geq_M;
	mna[22](6) = mna[22][6] - idt*Geq_M;
	mna[22](22) = mna[22][22] + idt*Geq_M;
	mna[190](22) = mna[190][22] +Geq_M;
	rval[22] += (Ieq_M - Qo_M)*idt;
	fvec[22] += -(Sr_M);
	mna[190](190) = mna[190][190] + 1.;
	rval[190] += Ieq_M;
	fvec[190] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM150*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[23];
	Un_M += Vx[6];
	Uo_M += -uxtold[23];
	Uo_M += uxtold[6];
	Qo_M = uxtold[191];
	Qn_M = Vx[191];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[191](6) = mna[191][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](23) = mna[6][23] - idt*Geq_M;
	mna[23](6) = mna[23][6] - idt*Geq_M;
	mna[23](23) = mna[23][23] + idt*Geq_M;
	mna[191](23) = mna[191][23] +Geq_M;
	rval[23] += (Ieq_M - Qo_M)*idt;
	fvec[23] += -(Sr_M);
	mna[191](191) = mna[191][191] + 1.;
	rval[191] += Ieq_M;
	fvec[191] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM151*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[24];
	Un_M += Vx[6];
	Uo_M += -uxtold[24];
	Uo_M += uxtold[6];
	Qo_M = uxtold[192];
	Qn_M = Vx[192];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[192](6) = mna[192][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](24) = mna[6][24] - idt*Geq_M;
	mna[24](6) = mna[24][6] - idt*Geq_M;
	mna[24](24) = mna[24][24] + idt*Geq_M;
	mna[192](24) = mna[192][24] +Geq_M;
	rval[24] += (Ieq_M - Qo_M)*idt;
	fvec[24] += -(Sr_M);
	mna[192](192) = mna[192][192] + 1.;
	rval[192] += Ieq_M;
	fvec[192] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM152*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[25];
	Un_M += Vx[6];
	Uo_M += -uxtold[25];
	Uo_M += uxtold[6];
	Qo_M = uxtold[193];
	Qn_M = Vx[193];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[193](6) = mna[193][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](25) = mna[6][25] - idt*Geq_M;
	mna[25](6) = mna[25][6] - idt*Geq_M;
	mna[25](25) = mna[25][25] + idt*Geq_M;
	mna[193](25) = mna[193][25] +Geq_M;
	rval[25] += (Ieq_M - Qo_M)*idt;
	fvec[25] += -(Sr_M);
	mna[193](193) = mna[193][193] + 1.;
	rval[193] += Ieq_M;
	fvec[193] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM153*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[26];
	Un_M += Vx[6];
	Uo_M += -uxtold[26];
	Uo_M += uxtold[6];
	Qo_M = uxtold[194];
	Qn_M = Vx[194];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[194](6) = mna[194][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](26) = mna[6][26] - idt*Geq_M;
	mna[26](6) = mna[26][6] - idt*Geq_M;
	mna[26](26) = mna[26][26] + idt*Geq_M;
	mna[194](26) = mna[194][26] +Geq_M;
	rval[26] += (Ieq_M - Qo_M)*idt;
	fvec[26] += -(Sr_M);
	mna[194](194) = mna[194][194] + 1.;
	rval[194] += Ieq_M;
	fvec[194] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM154*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[27];
	Un_M += Vx[6];
	Uo_M += -uxtold[27];
	Uo_M += uxtold[6];
	Qo_M = uxtold[195];
	Qn_M = Vx[195];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[195](6) = mna[195][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](27) = mna[6][27] - idt*Geq_M;
	mna[27](6) = mna[27][6] - idt*Geq_M;
	mna[27](27) = mna[27][27] + idt*Geq_M;
	mna[195](27) = mna[195][27] +Geq_M;
	rval[27] += (Ieq_M - Qo_M)*idt;
	fvec[27] += -(Sr_M);
	mna[195](195) = mna[195][195] + 1.;
	rval[195] += Ieq_M;
	fvec[195] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM155*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[31];
	Un_M += Vx[6];
	Uo_M += -uxtold[31];
	Uo_M += uxtold[6];
	Qo_M = uxtold[196];
	Qn_M = Vx[196];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[196](6) = mna[196][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](31) = mna[6][31] - idt*Geq_M;
	mna[31](6) = mna[31][6] - idt*Geq_M;
	mna[31](31) = mna[31][31] + idt*Geq_M;
	mna[196](31) = mna[196][31] +Geq_M;
	rval[31] += (Ieq_M - Qo_M)*idt;
	fvec[31] += -(Sr_M);
	mna[196](196) = mna[196][196] + 1.;
	rval[196] += Ieq_M;
	fvec[196] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM156*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[32];
	Un_M += Vx[6];
	Uo_M += -uxtold[32];
	Uo_M += uxtold[6];
	Qo_M = uxtold[197];
	Qn_M = Vx[197];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[197](6) = mna[197][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](32) = mna[6][32] - idt*Geq_M;
	mna[32](6) = mna[32][6] - idt*Geq_M;
	mna[32](32) = mna[32][32] + idt*Geq_M;
	mna[197](32) = mna[197][32] +Geq_M;
	rval[32] += (Ieq_M - Qo_M)*idt;
	fvec[32] += -(Sr_M);
	mna[197](197) = mna[197][197] + 1.;
	rval[197] += Ieq_M;
	fvec[197] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM157*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[34];
	Un_M += Vx[6];
	Uo_M += -uxtold[34];
	Uo_M += uxtold[6];
	Qo_M = uxtold[198];
	Qn_M = Vx[198];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[198](6) = mna[198][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](34) = mna[6][34] - idt*Geq_M;
	mna[34](6) = mna[34][6] - idt*Geq_M;
	mna[34](34) = mna[34][34] + idt*Geq_M;
	mna[198](34) = mna[198][34] +Geq_M;
	rval[34] += (Ieq_M - Qo_M)*idt;
	fvec[34] += -(Sr_M);
	mna[198](198) = mna[198][198] + 1.;
	rval[198] += Ieq_M;
	fvec[198] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM158*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[35];
	Un_M += Vx[6];
	Uo_M += -uxtold[35];
	Uo_M += uxtold[6];
	Qo_M = uxtold[199];
	Qn_M = Vx[199];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[199](6) = mna[199][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](35) = mna[6][35] - idt*Geq_M;
	mna[35](6) = mna[35][6] - idt*Geq_M;
	mna[35](35) = mna[35][35] + idt*Geq_M;
	mna[199](35) = mna[199][35] +Geq_M;
	rval[35] += (Ieq_M - Qo_M)*idt;
	fvec[35] += -(Sr_M);
	mna[199](199) = mna[199][199] + 1.;
	rval[199] += Ieq_M;
	fvec[199] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM159*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[36];
	Un_M += Vx[6];
	Uo_M += -uxtold[36];
	Uo_M += uxtold[6];
	Qo_M = uxtold[200];
	Qn_M = Vx[200];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[200](6) = mna[200][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](36) = mna[6][36] - idt*Geq_M;
	mna[36](6) = mna[36][6] - idt*Geq_M;
	mna[36](36) = mna[36][36] + idt*Geq_M;
	mna[200](36) = mna[200][36] +Geq_M;
	rval[36] += (Ieq_M - Qo_M)*idt;
	fvec[36] += -(Sr_M);
	mna[200](200) = mna[200][200] + 1.;
	rval[200] += Ieq_M;
	fvec[200] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM160*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[38];
	Un_M += Vx[6];
	Uo_M += -uxtold[38];
	Uo_M += uxtold[6];
	Qo_M = uxtold[201];
	Qn_M = Vx[201];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[201](6) = mna[201][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](38) = mna[6][38] - idt*Geq_M;
	mna[38](6) = mna[38][6] - idt*Geq_M;
	mna[38](38) = mna[38][38] + idt*Geq_M;
	mna[201](38) = mna[201][38] +Geq_M;
	rval[38] += (Ieq_M - Qo_M)*idt;
	fvec[38] += -(Sr_M);
	mna[201](201) = mna[201][201] + 1.;
	rval[201] += Ieq_M;
	fvec[201] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM161*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[39];
	Un_M += Vx[6];
	Uo_M += -uxtold[39];
	Uo_M += uxtold[6];
	Qo_M = uxtold[202];
	Qn_M = Vx[202];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[202](6) = mna[202][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](39) = mna[6][39] - idt*Geq_M;
	mna[39](6) = mna[39][6] - idt*Geq_M;
	mna[39](39) = mna[39][39] + idt*Geq_M;
	mna[202](39) = mna[202][39] +Geq_M;
	rval[39] += (Ieq_M - Qo_M)*idt;
	fvec[39] += -(Sr_M);
	mna[202](202) = mna[202][202] + 1.;
	rval[202] += Ieq_M;
	fvec[202] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM162*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[40];
	Un_M += Vx[6];
	Uo_M += -uxtold[40];
	Uo_M += uxtold[6];
	Qo_M = uxtold[203];
	Qn_M = Vx[203];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[6](6) = mna[6][6] + idt*Geq_M;
	mna[203](6) = mna[203][6] -Geq_M;
	rval[6] += -(Ieq_M - Qo_M)*idt;
	fvec[6] += Sr_M;
	mna[6](40) = mna[6][40] - idt*Geq_M;
	mna[40](6) = mna[40][6] - idt*Geq_M;
	mna[40](40) = mna[40][40] + idt*Geq_M;
	mna[203](40) = mna[203][40] +Geq_M;
	rval[40] += (Ieq_M - Qo_M)*idt;
	fvec[40] += -(Sr_M);
	mna[203](203) = mna[203][203] + 1.;
	rval[203] += Ieq_M;
	fvec[203] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM163*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[8];
	Un_M += Vx[7];
	Uo_M += -uxtold[8];
	Uo_M += uxtold[7];
	Qo_M = uxtold[204];
	Qn_M = Vx[204];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[204](7) = mna[204][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](8) = mna[7][8] - idt*Geq_M;
	mna[8](7) = mna[8][7] - idt*Geq_M;
	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[204](8) = mna[204][8] +Geq_M;
	rval[8] += (Ieq_M - Qo_M)*idt;
	fvec[8] += -(Sr_M);
	mna[204](204) = mna[204][204] + 1.;
	rval[204] += Ieq_M;
	fvec[204] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM164*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[9];
	Un_M += Vx[7];
	Uo_M += -uxtold[9];
	Uo_M += uxtold[7];
	Qo_M = uxtold[205];
	Qn_M = Vx[205];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[205](7) = mna[205][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](9) = mna[7][9] - idt*Geq_M;
	mna[9](7) = mna[9][7] - idt*Geq_M;
	mna[9](9) = mna[9][9] + idt*Geq_M;
	mna[205](9) = mna[205][9] +Geq_M;
	rval[9] += (Ieq_M - Qo_M)*idt;
	fvec[9] += -(Sr_M);
	mna[205](205) = mna[205][205] + 1.;
	rval[205] += Ieq_M;
	fvec[205] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM165*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[11];
	Un_M += Vx[7];
	Uo_M += -uxtold[11];
	Uo_M += uxtold[7];
	Qo_M = uxtold[206];
	Qn_M = Vx[206];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[206](7) = mna[206][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](11) = mna[7][11] - idt*Geq_M;
	mna[11](7) = mna[11][7] - idt*Geq_M;
	mna[11](11) = mna[11][11] + idt*Geq_M;
	mna[206](11) = mna[206][11] +Geq_M;
	rval[11] += (Ieq_M - Qo_M)*idt;
	fvec[11] += -(Sr_M);
	mna[206](206) = mna[206][206] + 1.;
	rval[206] += Ieq_M;
	fvec[206] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM166*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[12];
	Un_M += Vx[7];
	Uo_M += -uxtold[12];
	Uo_M += uxtold[7];
	Qo_M = uxtold[207];
	Qn_M = Vx[207];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[207](7) = mna[207][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](12) = mna[7][12] - idt*Geq_M;
	mna[12](7) = mna[12][7] - idt*Geq_M;
	mna[12](12) = mna[12][12] + idt*Geq_M;
	mna[207](12) = mna[207][12] +Geq_M;
	rval[12] += (Ieq_M - Qo_M)*idt;
	fvec[12] += -(Sr_M);
	mna[207](207) = mna[207][207] + 1.;
	rval[207] += Ieq_M;
	fvec[207] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM167*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[15];
	Un_M += Vx[7];
	Uo_M += -uxtold[15];
	Uo_M += uxtold[7];
	Qo_M = uxtold[208];
	Qn_M = Vx[208];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[208](7) = mna[208][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](15) = mna[7][15] - idt*Geq_M;
	mna[15](7) = mna[15][7] - idt*Geq_M;
	mna[15](15) = mna[15][15] + idt*Geq_M;
	mna[208](15) = mna[208][15] +Geq_M;
	rval[15] += (Ieq_M - Qo_M)*idt;
	fvec[15] += -(Sr_M);
	mna[208](208) = mna[208][208] + 1.;
	rval[208] += Ieq_M;
	fvec[208] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM168*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[17];
	Un_M += Vx[7];
	Uo_M += -uxtold[17];
	Uo_M += uxtold[7];
	Qo_M = uxtold[209];
	Qn_M = Vx[209];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[209](7) = mna[209][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](17) = mna[7][17] - idt*Geq_M;
	mna[17](7) = mna[17][7] - idt*Geq_M;
	mna[17](17) = mna[17][17] + idt*Geq_M;
	mna[209](17) = mna[209][17] +Geq_M;
	rval[17] += (Ieq_M - Qo_M)*idt;
	fvec[17] += -(Sr_M);
	mna[209](209) = mna[209][209] + 1.;
	rval[209] += Ieq_M;
	fvec[209] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM169*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[20];
	Un_M += Vx[7];
	Uo_M += -uxtold[20];
	Uo_M += uxtold[7];
	Qo_M = uxtold[210];
	Qn_M = Vx[210];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[210](7) = mna[210][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](20) = mna[7][20] - idt*Geq_M;
	mna[20](7) = mna[20][7] - idt*Geq_M;
	mna[20](20) = mna[20][20] + idt*Geq_M;
	mna[210](20) = mna[210][20] +Geq_M;
	rval[20] += (Ieq_M - Qo_M)*idt;
	fvec[20] += -(Sr_M);
	mna[210](210) = mna[210][210] + 1.;
	rval[210] += Ieq_M;
	fvec[210] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM170*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[21];
	Un_M += Vx[7];
	Uo_M += -uxtold[21];
	Uo_M += uxtold[7];
	Qo_M = uxtold[211];
	Qn_M = Vx[211];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[211](7) = mna[211][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](21) = mna[7][21] - idt*Geq_M;
	mna[21](7) = mna[21][7] - idt*Geq_M;
	mna[21](21) = mna[21][21] + idt*Geq_M;
	mna[211](21) = mna[211][21] +Geq_M;
	rval[21] += (Ieq_M - Qo_M)*idt;
	fvec[21] += -(Sr_M);
	mna[211](211) = mna[211][211] + 1.;
	rval[211] += Ieq_M;
	fvec[211] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM171*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[22];
	Un_M += Vx[7];
	Uo_M += -uxtold[22];
	Uo_M += uxtold[7];
	Qo_M = uxtold[212];
	Qn_M = Vx[212];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[212](7) = mna[212][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](22) = mna[7][22] - idt*Geq_M;
	mna[22](7) = mna[22][7] - idt*Geq_M;
	mna[22](22) = mna[22][22] + idt*Geq_M;
	mna[212](22) = mna[212][22] +Geq_M;
	rval[22] += (Ieq_M - Qo_M)*idt;
	fvec[22] += -(Sr_M);
	mna[212](212) = mna[212][212] + 1.;
	rval[212] += Ieq_M;
	fvec[212] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM172*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[23];
	Un_M += Vx[7];
	Uo_M += -uxtold[23];
	Uo_M += uxtold[7];
	Qo_M = uxtold[213];
	Qn_M = Vx[213];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[213](7) = mna[213][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](23) = mna[7][23] - idt*Geq_M;
	mna[23](7) = mna[23][7] - idt*Geq_M;
	mna[23](23) = mna[23][23] + idt*Geq_M;
	mna[213](23) = mna[213][23] +Geq_M;
	rval[23] += (Ieq_M - Qo_M)*idt;
	fvec[23] += -(Sr_M);
	mna[213](213) = mna[213][213] + 1.;
	rval[213] += Ieq_M;
	fvec[213] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM173*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[24];
	Un_M += Vx[7];
	Uo_M += -uxtold[24];
	Uo_M += uxtold[7];
	Qo_M = uxtold[214];
	Qn_M = Vx[214];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[214](7) = mna[214][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](24) = mna[7][24] - idt*Geq_M;
	mna[24](7) = mna[24][7] - idt*Geq_M;
	mna[24](24) = mna[24][24] + idt*Geq_M;
	mna[214](24) = mna[214][24] +Geq_M;
	rval[24] += (Ieq_M - Qo_M)*idt;
	fvec[24] += -(Sr_M);
	mna[214](214) = mna[214][214] + 1.;
	rval[214] += Ieq_M;
	fvec[214] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM174*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[25];
	Un_M += Vx[7];
	Uo_M += -uxtold[25];
	Uo_M += uxtold[7];
	Qo_M = uxtold[215];
	Qn_M = Vx[215];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[215](7) = mna[215][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](25) = mna[7][25] - idt*Geq_M;
	mna[25](7) = mna[25][7] - idt*Geq_M;
	mna[25](25) = mna[25][25] + idt*Geq_M;
	mna[215](25) = mna[215][25] +Geq_M;
	rval[25] += (Ieq_M - Qo_M)*idt;
	fvec[25] += -(Sr_M);
	mna[215](215) = mna[215][215] + 1.;
	rval[215] += Ieq_M;
	fvec[215] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM175*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[28];
	Un_M += Vx[7];
	Uo_M += -uxtold[28];
	Uo_M += uxtold[7];
	Qo_M = uxtold[216];
	Qn_M = Vx[216];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[216](7) = mna[216][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](28) = mna[7][28] - idt*Geq_M;
	mna[28](7) = mna[28][7] - idt*Geq_M;
	mna[28](28) = mna[28][28] + idt*Geq_M;
	mna[216](28) = mna[216][28] +Geq_M;
	rval[28] += (Ieq_M - Qo_M)*idt;
	fvec[28] += -(Sr_M);
	mna[216](216) = mna[216][216] + 1.;
	rval[216] += Ieq_M;
	fvec[216] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM176*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[29];
	Un_M += Vx[7];
	Uo_M += -uxtold[29];
	Uo_M += uxtold[7];
	Qo_M = uxtold[217];
	Qn_M = Vx[217];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[217](7) = mna[217][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](29) = mna[7][29] - idt*Geq_M;
	mna[29](7) = mna[29][7] - idt*Geq_M;
	mna[29](29) = mna[29][29] + idt*Geq_M;
	mna[217](29) = mna[217][29] +Geq_M;
	rval[29] += (Ieq_M - Qo_M)*idt;
	fvec[29] += -(Sr_M);
	mna[217](217) = mna[217][217] + 1.;
	rval[217] += Ieq_M;
	fvec[217] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM177*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[30];
	Un_M += Vx[7];
	Uo_M += -uxtold[30];
	Uo_M += uxtold[7];
	Qo_M = uxtold[218];
	Qn_M = Vx[218];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[218](7) = mna[218][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](30) = mna[7][30] - idt*Geq_M;
	mna[30](7) = mna[30][7] - idt*Geq_M;
	mna[30](30) = mna[30][30] + idt*Geq_M;
	mna[218](30) = mna[218][30] +Geq_M;
	rval[30] += (Ieq_M - Qo_M)*idt;
	fvec[30] += -(Sr_M);
	mna[218](218) = mna[218][218] + 1.;
	rval[218] += Ieq_M;
	fvec[218] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM178*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[31];
	Un_M += Vx[7];
	Uo_M += -uxtold[31];
	Uo_M += uxtold[7];
	Qo_M = uxtold[219];
	Qn_M = Vx[219];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[219](7) = mna[219][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](31) = mna[7][31] - idt*Geq_M;
	mna[31](7) = mna[31][7] - idt*Geq_M;
	mna[31](31) = mna[31][31] + idt*Geq_M;
	mna[219](31) = mna[219][31] +Geq_M;
	rval[31] += (Ieq_M - Qo_M)*idt;
	fvec[31] += -(Sr_M);
	mna[219](219) = mna[219][219] + 1.;
	rval[219] += Ieq_M;
	fvec[219] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM179*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[33];
	Un_M += Vx[7];
	Uo_M += -uxtold[33];
	Uo_M += uxtold[7];
	Qo_M = uxtold[220];
	Qn_M = Vx[220];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[220](7) = mna[220][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](33) = mna[7][33] - idt*Geq_M;
	mna[33](7) = mna[33][7] - idt*Geq_M;
	mna[33](33) = mna[33][33] + idt*Geq_M;
	mna[220](33) = mna[220][33] +Geq_M;
	rval[33] += (Ieq_M - Qo_M)*idt;
	fvec[33] += -(Sr_M);
	mna[220](220) = mna[220][220] + 1.;
	rval[220] += Ieq_M;
	fvec[220] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM180*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[34];
	Un_M += Vx[7];
	Uo_M += -uxtold[34];
	Uo_M += uxtold[7];
	Qo_M = uxtold[221];
	Qn_M = Vx[221];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[221](7) = mna[221][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](34) = mna[7][34] - idt*Geq_M;
	mna[34](7) = mna[34][7] - idt*Geq_M;
	mna[34](34) = mna[34][34] + idt*Geq_M;
	mna[221](34) = mna[221][34] +Geq_M;
	rval[34] += (Ieq_M - Qo_M)*idt;
	fvec[34] += -(Sr_M);
	mna[221](221) = mna[221][221] + 1.;
	rval[221] += Ieq_M;
	fvec[221] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM181*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[36];
	Un_M += Vx[7];
	Uo_M += -uxtold[36];
	Uo_M += uxtold[7];
	Qo_M = uxtold[222];
	Qn_M = Vx[222];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[222](7) = mna[222][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](36) = mna[7][36] - idt*Geq_M;
	mna[36](7) = mna[36][7] - idt*Geq_M;
	mna[36](36) = mna[36][36] + idt*Geq_M;
	mna[222](36) = mna[222][36] +Geq_M;
	rval[36] += (Ieq_M - Qo_M)*idt;
	fvec[36] += -(Sr_M);
	mna[222](222) = mna[222][222] + 1.;
	rval[222] += Ieq_M;
	fvec[222] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM182*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[37];
	Un_M += Vx[7];
	Uo_M += -uxtold[37];
	Uo_M += uxtold[7];
	Qo_M = uxtold[223];
	Qn_M = Vx[223];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[223](7) = mna[223][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](37) = mna[7][37] - idt*Geq_M;
	mna[37](7) = mna[37][7] - idt*Geq_M;
	mna[37](37) = mna[37][37] + idt*Geq_M;
	mna[223](37) = mna[223][37] +Geq_M;
	rval[37] += (Ieq_M - Qo_M)*idt;
	fvec[37] += -(Sr_M);
	mna[223](223) = mna[223][223] + 1.;
	rval[223] += Ieq_M;
	fvec[223] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM183*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[38];
	Un_M += Vx[7];
	Uo_M += -uxtold[38];
	Uo_M += uxtold[7];
	Qo_M = uxtold[224];
	Qn_M = Vx[224];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[224](7) = mna[224][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](38) = mna[7][38] - idt*Geq_M;
	mna[38](7) = mna[38][7] - idt*Geq_M;
	mna[38](38) = mna[38][38] + idt*Geq_M;
	mna[224](38) = mna[224][38] +Geq_M;
	rval[38] += (Ieq_M - Qo_M)*idt;
	fvec[38] += -(Sr_M);
	mna[224](224) = mna[224][224] + 1.;
	rval[224] += Ieq_M;
	fvec[224] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM184*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[39];
	Un_M += Vx[7];
	Uo_M += -uxtold[39];
	Uo_M += uxtold[7];
	Qo_M = uxtold[225];
	Qn_M = Vx[225];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[225](7) = mna[225][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](39) = mna[7][39] - idt*Geq_M;
	mna[39](7) = mna[39][7] - idt*Geq_M;
	mna[39](39) = mna[39][39] + idt*Geq_M;
	mna[225](39) = mna[225][39] +Geq_M;
	rval[39] += (Ieq_M - Qo_M)*idt;
	fvec[39] += -(Sr_M);
	mna[225](225) = mna[225][225] + 1.;
	rval[225] += Ieq_M;
	fvec[225] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM185*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[40];
	Un_M += Vx[7];
	Uo_M += -uxtold[40];
	Uo_M += uxtold[7];
	Qo_M = uxtold[226];
	Qn_M = Vx[226];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[7](7) = mna[7][7] + idt*Geq_M;
	mna[226](7) = mna[226][7] -Geq_M;
	rval[7] += -(Ieq_M - Qo_M)*idt;
	fvec[7] += Sr_M;
	mna[7](40) = mna[7][40] - idt*Geq_M;
	mna[40](7) = mna[40][7] - idt*Geq_M;
	mna[40](40) = mna[40][40] + idt*Geq_M;
	mna[226](40) = mna[226][40] +Geq_M;
	rval[40] += (Ieq_M - Qo_M)*idt;
	fvec[40] += -(Sr_M);
	mna[226](226) = mna[226][226] + 1.;
	rval[226] += Ieq_M;
	fvec[226] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM186*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[9];
	Un_M += Vx[8];
	Uo_M += -uxtold[9];
	Uo_M += uxtold[8];
	Qo_M = uxtold[227];
	Qn_M = Vx[227];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[227](8) = mna[227][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](9) = mna[8][9] - idt*Geq_M;
	mna[9](8) = mna[9][8] - idt*Geq_M;
	mna[9](9) = mna[9][9] + idt*Geq_M;
	mna[227](9) = mna[227][9] +Geq_M;
	rval[9] += (Ieq_M - Qo_M)*idt;
	fvec[9] += -(Sr_M);
	mna[227](227) = mna[227][227] + 1.;
	rval[227] += Ieq_M;
	fvec[227] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM187*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[10];
	Un_M += Vx[8];
	Uo_M += -uxtold[10];
	Uo_M += uxtold[8];
	Qo_M = uxtold[228];
	Qn_M = Vx[228];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[228](8) = mna[228][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](10) = mna[8][10] - idt*Geq_M;
	mna[10](8) = mna[10][8] - idt*Geq_M;
	mna[10](10) = mna[10][10] + idt*Geq_M;
	mna[228](10) = mna[228][10] +Geq_M;
	rval[10] += (Ieq_M - Qo_M)*idt;
	fvec[10] += -(Sr_M);
	mna[228](228) = mna[228][228] + 1.;
	rval[228] += Ieq_M;
	fvec[228] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM188*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[11];
	Un_M += Vx[8];
	Uo_M += -uxtold[11];
	Uo_M += uxtold[8];
	Qo_M = uxtold[229];
	Qn_M = Vx[229];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[229](8) = mna[229][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](11) = mna[8][11] - idt*Geq_M;
	mna[11](8) = mna[11][8] - idt*Geq_M;
	mna[11](11) = mna[11][11] + idt*Geq_M;
	mna[229](11) = mna[229][11] +Geq_M;
	rval[11] += (Ieq_M - Qo_M)*idt;
	fvec[11] += -(Sr_M);
	mna[229](229) = mna[229][229] + 1.;
	rval[229] += Ieq_M;
	fvec[229] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM189*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[12];
	Un_M += Vx[8];
	Uo_M += -uxtold[12];
	Uo_M += uxtold[8];
	Qo_M = uxtold[230];
	Qn_M = Vx[230];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[230](8) = mna[230][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](12) = mna[8][12] - idt*Geq_M;
	mna[12](8) = mna[12][8] - idt*Geq_M;
	mna[12](12) = mna[12][12] + idt*Geq_M;
	mna[230](12) = mna[230][12] +Geq_M;
	rval[12] += (Ieq_M - Qo_M)*idt;
	fvec[12] += -(Sr_M);
	mna[230](230) = mna[230][230] + 1.;
	rval[230] += Ieq_M;
	fvec[230] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM190*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[14];
	Un_M += Vx[8];
	Uo_M += -uxtold[14];
	Uo_M += uxtold[8];
	Qo_M = uxtold[231];
	Qn_M = Vx[231];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[231](8) = mna[231][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](14) = mna[8][14] - idt*Geq_M;
	mna[14](8) = mna[14][8] - idt*Geq_M;
	mna[14](14) = mna[14][14] + idt*Geq_M;
	mna[231](14) = mna[231][14] +Geq_M;
	rval[14] += (Ieq_M - Qo_M)*idt;
	fvec[14] += -(Sr_M);
	mna[231](231) = mna[231][231] + 1.;
	rval[231] += Ieq_M;
	fvec[231] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM191*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[15];
	Un_M += Vx[8];
	Uo_M += -uxtold[15];
	Uo_M += uxtold[8];
	Qo_M = uxtold[232];
	Qn_M = Vx[232];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[232](8) = mna[232][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](15) = mna[8][15] - idt*Geq_M;
	mna[15](8) = mna[15][8] - idt*Geq_M;
	mna[15](15) = mna[15][15] + idt*Geq_M;
	mna[232](15) = mna[232][15] +Geq_M;
	rval[15] += (Ieq_M - Qo_M)*idt;
	fvec[15] += -(Sr_M);
	mna[232](232) = mna[232][232] + 1.;
	rval[232] += Ieq_M;
	fvec[232] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM192*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[16];
	Un_M += Vx[8];
	Uo_M += -uxtold[16];
	Uo_M += uxtold[8];
	Qo_M = uxtold[233];
	Qn_M = Vx[233];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[233](8) = mna[233][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](16) = mna[8][16] - idt*Geq_M;
	mna[16](8) = mna[16][8] - idt*Geq_M;
	mna[16](16) = mna[16][16] + idt*Geq_M;
	mna[233](16) = mna[233][16] +Geq_M;
	rval[16] += (Ieq_M - Qo_M)*idt;
	fvec[16] += -(Sr_M);
	mna[233](233) = mna[233][233] + 1.;
	rval[233] += Ieq_M;
	fvec[233] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM193*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[18];
	Un_M += Vx[8];
	Uo_M += -uxtold[18];
	Uo_M += uxtold[8];
	Qo_M = uxtold[234];
	Qn_M = Vx[234];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[234](8) = mna[234][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](18) = mna[8][18] - idt*Geq_M;
	mna[18](8) = mna[18][8] - idt*Geq_M;
	mna[18](18) = mna[18][18] + idt*Geq_M;
	mna[234](18) = mna[234][18] +Geq_M;
	rval[18] += (Ieq_M - Qo_M)*idt;
	fvec[18] += -(Sr_M);
	mna[234](234) = mna[234][234] + 1.;
	rval[234] += Ieq_M;
	fvec[234] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM194*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[19];
	Un_M += Vx[8];
	Uo_M += -uxtold[19];
	Uo_M += uxtold[8];
	Qo_M = uxtold[235];
	Qn_M = Vx[235];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[235](8) = mna[235][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](19) = mna[8][19] - idt*Geq_M;
	mna[19](8) = mna[19][8] - idt*Geq_M;
	mna[19](19) = mna[19][19] + idt*Geq_M;
	mna[235](19) = mna[235][19] +Geq_M;
	rval[19] += (Ieq_M - Qo_M)*idt;
	fvec[19] += -(Sr_M);
	mna[235](235) = mna[235][235] + 1.;
	rval[235] += Ieq_M;
	fvec[235] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM195*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[20];
	Un_M += Vx[8];
	Uo_M += -uxtold[20];
	Uo_M += uxtold[8];
	Qo_M = uxtold[236];
	Qn_M = Vx[236];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[236](8) = mna[236][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](20) = mna[8][20] - idt*Geq_M;
	mna[20](8) = mna[20][8] - idt*Geq_M;
	mna[20](20) = mna[20][20] + idt*Geq_M;
	mna[236](20) = mna[236][20] +Geq_M;
	rval[20] += (Ieq_M - Qo_M)*idt;
	fvec[20] += -(Sr_M);
	mna[236](236) = mna[236][236] + 1.;
	rval[236] += Ieq_M;
	fvec[236] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM196*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[22];
	Un_M += Vx[8];
	Uo_M += -uxtold[22];
	Uo_M += uxtold[8];
	Qo_M = uxtold[237];
	Qn_M = Vx[237];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[237](8) = mna[237][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](22) = mna[8][22] - idt*Geq_M;
	mna[22](8) = mna[22][8] - idt*Geq_M;
	mna[22](22) = mna[22][22] + idt*Geq_M;
	mna[237](22) = mna[237][22] +Geq_M;
	rval[22] += (Ieq_M - Qo_M)*idt;
	fvec[22] += -(Sr_M);
	mna[237](237) = mna[237][237] + 1.;
	rval[237] += Ieq_M;
	fvec[237] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM197*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[23];
	Un_M += Vx[8];
	Uo_M += -uxtold[23];
	Uo_M += uxtold[8];
	Qo_M = uxtold[238];
	Qn_M = Vx[238];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[238](8) = mna[238][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](23) = mna[8][23] - idt*Geq_M;
	mna[23](8) = mna[23][8] - idt*Geq_M;
	mna[23](23) = mna[23][23] + idt*Geq_M;
	mna[238](23) = mna[238][23] +Geq_M;
	rval[23] += (Ieq_M - Qo_M)*idt;
	fvec[23] += -(Sr_M);
	mna[238](238) = mna[238][238] + 1.;
	rval[238] += Ieq_M;
	fvec[238] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM198*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[24];
	Un_M += Vx[8];
	Uo_M += -uxtold[24];
	Uo_M += uxtold[8];
	Qo_M = uxtold[239];
	Qn_M = Vx[239];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[239](8) = mna[239][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](24) = mna[8][24] - idt*Geq_M;
	mna[24](8) = mna[24][8] - idt*Geq_M;
	mna[24](24) = mna[24][24] + idt*Geq_M;
	mna[239](24) = mna[239][24] +Geq_M;
	rval[24] += (Ieq_M - Qo_M)*idt;
	fvec[24] += -(Sr_M);
	mna[239](239) = mna[239][239] + 1.;
	rval[239] += Ieq_M;
	fvec[239] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM199*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[26];
	Un_M += Vx[8];
	Uo_M += -uxtold[26];
	Uo_M += uxtold[8];
	Qo_M = uxtold[240];
	Qn_M = Vx[240];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[240](8) = mna[240][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](26) = mna[8][26] - idt*Geq_M;
	mna[26](8) = mna[26][8] - idt*Geq_M;
	mna[26](26) = mna[26][26] + idt*Geq_M;
	mna[240](26) = mna[240][26] +Geq_M;
	rval[26] += (Ieq_M - Qo_M)*idt;
	fvec[26] += -(Sr_M);
	mna[240](240) = mna[240][240] + 1.;
	rval[240] += Ieq_M;
	fvec[240] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM200*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[27];
	Un_M += Vx[8];
	Uo_M += -uxtold[27];
	Uo_M += uxtold[8];
	Qo_M = uxtold[241];
	Qn_M = Vx[241];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[241](8) = mna[241][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](27) = mna[8][27] - idt*Geq_M;
	mna[27](8) = mna[27][8] - idt*Geq_M;
	mna[27](27) = mna[27][27] + idt*Geq_M;
	mna[241](27) = mna[241][27] +Geq_M;
	rval[27] += (Ieq_M - Qo_M)*idt;
	fvec[27] += -(Sr_M);
	mna[241](241) = mna[241][241] + 1.;
	rval[241] += Ieq_M;
	fvec[241] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM201*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[28];
	Un_M += Vx[8];
	Uo_M += -uxtold[28];
	Uo_M += uxtold[8];
	Qo_M = uxtold[242];
	Qn_M = Vx[242];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[242](8) = mna[242][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](28) = mna[8][28] - idt*Geq_M;
	mna[28](8) = mna[28][8] - idt*Geq_M;
	mna[28](28) = mna[28][28] + idt*Geq_M;
	mna[242](28) = mna[242][28] +Geq_M;
	rval[28] += (Ieq_M - Qo_M)*idt;
	fvec[28] += -(Sr_M);
	mna[242](242) = mna[242][242] + 1.;
	rval[242] += Ieq_M;
	fvec[242] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM202*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[30];
	Un_M += Vx[8];
	Uo_M += -uxtold[30];
	Uo_M += uxtold[8];
	Qo_M = uxtold[243];
	Qn_M = Vx[243];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[243](8) = mna[243][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](30) = mna[8][30] - idt*Geq_M;
	mna[30](8) = mna[30][8] - idt*Geq_M;
	mna[30](30) = mna[30][30] + idt*Geq_M;
	mna[243](30) = mna[243][30] +Geq_M;
	rval[30] += (Ieq_M - Qo_M)*idt;
	fvec[30] += -(Sr_M);
	mna[243](243) = mna[243][243] + 1.;
	rval[243] += Ieq_M;
	fvec[243] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM203*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[31];
	Un_M += Vx[8];
	Uo_M += -uxtold[31];
	Uo_M += uxtold[8];
	Qo_M = uxtold[244];
	Qn_M = Vx[244];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[244](8) = mna[244][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](31) = mna[8][31] - idt*Geq_M;
	mna[31](8) = mna[31][8] - idt*Geq_M;
	mna[31](31) = mna[31][31] + idt*Geq_M;
	mna[244](31) = mna[244][31] +Geq_M;
	rval[31] += (Ieq_M - Qo_M)*idt;
	fvec[31] += -(Sr_M);
	mna[244](244) = mna[244][244] + 1.;
	rval[244] += Ieq_M;
	fvec[244] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM204*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[32];
	Un_M += Vx[8];
	Uo_M += -uxtold[32];
	Uo_M += uxtold[8];
	Qo_M = uxtold[245];
	Qn_M = Vx[245];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[245](8) = mna[245][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](32) = mna[8][32] - idt*Geq_M;
	mna[32](8) = mna[32][8] - idt*Geq_M;
	mna[32](32) = mna[32][32] + idt*Geq_M;
	mna[245](32) = mna[245][32] +Geq_M;
	rval[32] += (Ieq_M - Qo_M)*idt;
	fvec[32] += -(Sr_M);
	mna[245](245) = mna[245][245] + 1.;
	rval[245] += Ieq_M;
	fvec[245] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM205*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[34];
	Un_M += Vx[8];
	Uo_M += -uxtold[34];
	Uo_M += uxtold[8];
	Qo_M = uxtold[246];
	Qn_M = Vx[246];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[246](8) = mna[246][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](34) = mna[8][34] - idt*Geq_M;
	mna[34](8) = mna[34][8] - idt*Geq_M;
	mna[34](34) = mna[34][34] + idt*Geq_M;
	mna[246](34) = mna[246][34] +Geq_M;
	rval[34] += (Ieq_M - Qo_M)*idt;
	fvec[34] += -(Sr_M);
	mna[246](246) = mna[246][246] + 1.;
	rval[246] += Ieq_M;
	fvec[246] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM206*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[35];
	Un_M += Vx[8];
	Uo_M += -uxtold[35];
	Uo_M += uxtold[8];
	Qo_M = uxtold[247];
	Qn_M = Vx[247];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[247](8) = mna[247][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](35) = mna[8][35] - idt*Geq_M;
	mna[35](8) = mna[35][8] - idt*Geq_M;
	mna[35](35) = mna[35][35] + idt*Geq_M;
	mna[247](35) = mna[247][35] +Geq_M;
	rval[35] += (Ieq_M - Qo_M)*idt;
	fvec[35] += -(Sr_M);
	mna[247](247) = mna[247][247] + 1.;
	rval[247] += Ieq_M;
	fvec[247] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM207*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[36];
	Un_M += Vx[8];
	Uo_M += -uxtold[36];
	Uo_M += uxtold[8];
	Qo_M = uxtold[248];
	Qn_M = Vx[248];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[248](8) = mna[248][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](36) = mna[8][36] - idt*Geq_M;
	mna[36](8) = mna[36][8] - idt*Geq_M;
	mna[36](36) = mna[36][36] + idt*Geq_M;
	mna[248](36) = mna[248][36] +Geq_M;
	rval[36] += (Ieq_M - Qo_M)*idt;
	fvec[36] += -(Sr_M);
	mna[248](248) = mna[248][248] + 1.;
	rval[248] += Ieq_M;
	fvec[248] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM208*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[37];
	Un_M += Vx[8];
	Uo_M += -uxtold[37];
	Uo_M += uxtold[8];
	Qo_M = uxtold[249];
	Qn_M = Vx[249];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[249](8) = mna[249][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](37) = mna[8][37] - idt*Geq_M;
	mna[37](8) = mna[37][8] - idt*Geq_M;
	mna[37](37) = mna[37][37] + idt*Geq_M;
	mna[249](37) = mna[249][37] +Geq_M;
	rval[37] += (Ieq_M - Qo_M)*idt;
	fvec[37] += -(Sr_M);
	mna[249](249) = mna[249][249] + 1.;
	rval[249] += Ieq_M;
	fvec[249] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM209*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[38];
	Un_M += Vx[8];
	Uo_M += -uxtold[38];
	Uo_M += uxtold[8];
	Qo_M = uxtold[250];
	Qn_M = Vx[250];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[250](8) = mna[250][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](38) = mna[8][38] - idt*Geq_M;
	mna[38](8) = mna[38][8] - idt*Geq_M;
	mna[38](38) = mna[38][38] + idt*Geq_M;
	mna[250](38) = mna[250][38] +Geq_M;
	rval[38] += (Ieq_M - Qo_M)*idt;
	fvec[38] += -(Sr_M);
	mna[250](250) = mna[250][250] + 1.;
	rval[250] += Ieq_M;
	fvec[250] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM210*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[39];
	Un_M += Vx[8];
	Uo_M += -uxtold[39];
	Uo_M += uxtold[8];
	Qo_M = uxtold[251];
	Qn_M = Vx[251];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[251](8) = mna[251][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](39) = mna[8][39] - idt*Geq_M;
	mna[39](8) = mna[39][8] - idt*Geq_M;
	mna[39](39) = mna[39][39] + idt*Geq_M;
	mna[251](39) = mna[251][39] +Geq_M;
	rval[39] += (Ieq_M - Qo_M)*idt;
	fvec[39] += -(Sr_M);
	mna[251](251) = mna[251][251] + 1.;
	rval[251] += Ieq_M;
	fvec[251] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM211*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[40];
	Un_M += Vx[8];
	Uo_M += -uxtold[40];
	Uo_M += uxtold[8];
	Qo_M = uxtold[252];
	Qn_M = Vx[252];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[8](8) = mna[8][8] + idt*Geq_M;
	mna[252](8) = mna[252][8] -Geq_M;
	rval[8] += -(Ieq_M - Qo_M)*idt;
	fvec[8] += Sr_M;
	mna[8](40) = mna[8][40] - idt*Geq_M;
	mna[40](8) = mna[40][8] - idt*Geq_M;
	mna[40](40) = mna[40][40] + idt*Geq_M;
	mna[252](40) = mna[252][40] +Geq_M;
	rval[40] += (Ieq_M - Qo_M)*idt;
	fvec[40] += -(Sr_M);
	mna[252](252) = mna[252][252] + 1.;
	rval[252] += Ieq_M;
	fvec[252] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM212*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[11];
	Un_M += Vx[9];
	Uo_M += -uxtold[11];
	Uo_M += uxtold[9];
	Qo_M = uxtold[253];
	Qn_M = Vx[253];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[9](9) = mna[9][9] + idt*Geq_M;
	mna[253](9) = mna[253][9] -Geq_M;
	rval[9] += -(Ieq_M - Qo_M)*idt;
	fvec[9] += Sr_M;
	mna[9](11) = mna[9][11] - idt*Geq_M;
	mna[11](9) = mna[11][9] - idt*Geq_M;
	mna[11](11) = mna[11][11] + idt*Geq_M;
	mna[253](11) = mna[253][11] +Geq_M;
	rval[11] += (Ieq_M - Qo_M)*idt;
	fvec[11] += -(Sr_M);
	mna[253](253) = mna[253][253] + 1.;
	rval[253] += Ieq_M;
	fvec[253] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM213*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[14];
	Un_M += Vx[9];
	Uo_M += -uxtold[14];
	Uo_M += uxtold[9];
	Qo_M = uxtold[254];
	Qn_M = Vx[254];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[9](9) = mna[9][9] + idt*Geq_M;
	mna[254](9) = mna[254][9] -Geq_M;
	rval[9] += -(Ieq_M - Qo_M)*idt;
	fvec[9] += Sr_M;
	mna[9](14) = mna[9][14] - idt*Geq_M;
	mna[14](9) = mna[14][9] - idt*Geq_M;
	mna[14](14) = mna[14][14] + idt*Geq_M;
	mna[254](14) = mna[254][14] +Geq_M;
	rval[14] += (Ieq_M - Qo_M)*idt;
	fvec[14] += -(Sr_M);
	mna[254](254) = mna[254][254] + 1.;
	rval[254] += Ieq_M;
	fvec[254] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM214*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[16];
	Un_M += Vx[9];
	Uo_M += -uxtold[16];
	Uo_M += uxtold[9];
	Qo_M = uxtold[255];
	Qn_M = Vx[255];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[9](9) = mna[9][9] + idt*Geq_M;
	mna[255](9) = mna[255][9] -Geq_M;
	rval[9] += -(Ieq_M - Qo_M)*idt;
	fvec[9] += Sr_M;
	mna[9](16) = mna[9][16] - idt*Geq_M;
	mna[16](9) = mna[16][9] - idt*Geq_M;
	mna[16](16) = mna[16][16] + idt*Geq_M;
	mna[255](16) = mna[255][16] +Geq_M;
	rval[16] += (Ieq_M - Qo_M)*idt;
	fvec[16] += -(Sr_M);
	mna[255](255) = mna[255][255] + 1.;
	rval[255] += Ieq_M;
	fvec[255] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM215*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[17];
	Un_M += Vx[9];
	Uo_M += -uxtold[17];
	Uo_M += uxtold[9];
	Qo_M = uxtold[256];
	Qn_M = Vx[256];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[9](9) = mna[9][9] + idt*Geq_M;
	mna[256](9) = mna[256][9] -Geq_M;
	rval[9] += -(Ieq_M - Qo_M)*idt;
	fvec[9] += Sr_M;
	mna[9](17) = mna[9][17] - idt*Geq_M;
	mna[17](9) = mna[17][9] - idt*Geq_M;
	mna[17](17) = mna[17][17] + idt*Geq_M;
	mna[256](17) = mna[256][17] +Geq_M;
	rval[17] += (Ieq_M - Qo_M)*idt;
	fvec[17] += -(Sr_M);
	mna[256](256) = mna[256][256] + 1.;
	rval[256] += Ieq_M;
	fvec[256] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM216*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[20];
	Un_M += Vx[9];
	Uo_M += -uxtold[20];
	Uo_M += uxtold[9];
	Qo_M = uxtold[257];
	Qn_M = Vx[257];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[9](9) = mna[9][9] + idt*Geq_M;
	mna[257](9) = mna[257][9] -Geq_M;
	rval[9] += -(Ieq_M - Qo_M)*idt;
	fvec[9] += Sr_M;
	mna[9](20) = mna[9][20] - idt*Geq_M;
	mna[20](9) = mna[20][9] - idt*Geq_M;
	mna[20](20) = mna[20][20] + idt*Geq_M;
	mna[257](20) = mna[257][20] +Geq_M;
	rval[20] += (Ieq_M - Qo_M)*idt;
	fvec[20] += -(Sr_M);
	mna[257](257) = mna[257][257] + 1.;
	rval[257] += Ieq_M;
	fvec[257] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM217*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[30];
	Un_M += Vx[9];
	Uo_M += -uxtold[30];
	Uo_M += uxtold[9];
	Qo_M = uxtold[258];
	Qn_M = Vx[258];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[9](9) = mna[9][9] + idt*Geq_M;
	mna[258](9) = mna[258][9] -Geq_M;
	rval[9] += -(Ieq_M - Qo_M)*idt;
	fvec[9] += Sr_M;
	mna[9](30) = mna[9][30] - idt*Geq_M;
	mna[30](9) = mna[30][9] - idt*Geq_M;
	mna[30](30) = mna[30][30] + idt*Geq_M;
	mna[258](30) = mna[258][30] +Geq_M;
	rval[30] += (Ieq_M - Qo_M)*idt;
	fvec[30] += -(Sr_M);
	mna[258](258) = mna[258][258] + 1.;
	rval[258] += Ieq_M;
	fvec[258] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM218*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[31];
	Un_M += Vx[9];
	Uo_M += -uxtold[31];
	Uo_M += uxtold[9];
	Qo_M = uxtold[259];
	Qn_M = Vx[259];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[9](9) = mna[9][9] + idt*Geq_M;
	mna[259](9) = mna[259][9] -Geq_M;
	rval[9] += -(Ieq_M - Qo_M)*idt;
	fvec[9] += Sr_M;
	mna[9](31) = mna[9][31] - idt*Geq_M;
	mna[31](9) = mna[31][9] - idt*Geq_M;
	mna[31](31) = mna[31][31] + idt*Geq_M;
	mna[259](31) = mna[259][31] +Geq_M;
	rval[31] += (Ieq_M - Qo_M)*idt;
	fvec[31] += -(Sr_M);
	mna[259](259) = mna[259][259] + 1.;
	rval[259] += Ieq_M;
	fvec[259] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM219*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[32];
	Un_M += Vx[9];
	Uo_M += -uxtold[32];
	Uo_M += uxtold[9];
	Qo_M = uxtold[260];
	Qn_M = Vx[260];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[9](9) = mna[9][9] + idt*Geq_M;
	mna[260](9) = mna[260][9] -Geq_M;
	rval[9] += -(Ieq_M - Qo_M)*idt;
	fvec[9] += Sr_M;
	mna[9](32) = mna[9][32] - idt*Geq_M;
	mna[32](9) = mna[32][9] - idt*Geq_M;
	mna[32](32) = mna[32][32] + idt*Geq_M;
	mna[260](32) = mna[260][32] +Geq_M;
	rval[32] += (Ieq_M - Qo_M)*idt;
	fvec[32] += -(Sr_M);
	mna[260](260) = mna[260][260] + 1.;
	rval[260] += Ieq_M;
	fvec[260] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM220*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[34];
	Un_M += Vx[9];
	Uo_M += -uxtold[34];
	Uo_M += uxtold[9];
	Qo_M = uxtold[261];
	Qn_M = Vx[261];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[9](9) = mna[9][9] + idt*Geq_M;
	mna[261](9) = mna[261][9] -Geq_M;
	rval[9] += -(Ieq_M - Qo_M)*idt;
	fvec[9] += Sr_M;
	mna[9](34) = mna[9][34] - idt*Geq_M;
	mna[34](9) = mna[34][9] - idt*Geq_M;
	mna[34](34) = mna[34][34] + idt*Geq_M;
	mna[261](34) = mna[261][34] +Geq_M;
	rval[34] += (Ieq_M - Qo_M)*idt;
	fvec[34] += -(Sr_M);
	mna[261](261) = mna[261][261] + 1.;
	rval[261] += Ieq_M;
	fvec[261] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM221*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[11];
	Un_M += Vx[10];
	Uo_M += -uxtold[11];
	Uo_M += uxtold[10];
	Qo_M = uxtold[262];
	Qn_M = Vx[262];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[10](10) = mna[10][10] + idt*Geq_M;
	mna[262](10) = mna[262][10] -Geq_M;
	rval[10] += -(Ieq_M - Qo_M)*idt;
	fvec[10] += Sr_M;
	mna[10](11) = mna[10][11] - idt*Geq_M;
	mna[11](10) = mna[11][10] - idt*Geq_M;
	mna[11](11) = mna[11][11] + idt*Geq_M;
	mna[262](11) = mna[262][11] +Geq_M;
	rval[11] += (Ieq_M - Qo_M)*idt;
	fvec[11] += -(Sr_M);
	mna[262](262) = mna[262][262] + 1.;
	rval[262] += Ieq_M;
	fvec[262] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM222*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[12];
	Un_M += Vx[10];
	Uo_M += -uxtold[12];
	Uo_M += uxtold[10];
	Qo_M = uxtold[263];
	Qn_M = Vx[263];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[10](10) = mna[10][10] + idt*Geq_M;
	mna[263](10) = mna[263][10] -Geq_M;
	rval[10] += -(Ieq_M - Qo_M)*idt;
	fvec[10] += Sr_M;
	mna[10](12) = mna[10][12] - idt*Geq_M;
	mna[12](10) = mna[12][10] - idt*Geq_M;
	mna[12](12) = mna[12][12] + idt*Geq_M;
	mna[263](12) = mna[263][12] +Geq_M;
	rval[12] += (Ieq_M - Qo_M)*idt;
	fvec[12] += -(Sr_M);
	mna[263](263) = mna[263][263] + 1.;
	rval[263] += Ieq_M;
	fvec[263] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM223*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[15];
	Un_M += Vx[10];
	Uo_M += -uxtold[15];
	Uo_M += uxtold[10];
	Qo_M = uxtold[264];
	Qn_M = Vx[264];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[10](10) = mna[10][10] + idt*Geq_M;
	mna[264](10) = mna[264][10] -Geq_M;
	rval[10] += -(Ieq_M - Qo_M)*idt;
	fvec[10] += Sr_M;
	mna[10](15) = mna[10][15] - idt*Geq_M;
	mna[15](10) = mna[15][10] - idt*Geq_M;
	mna[15](15) = mna[15][15] + idt*Geq_M;
	mna[264](15) = mna[264][15] +Geq_M;
	rval[15] += (Ieq_M - Qo_M)*idt;
	fvec[15] += -(Sr_M);
	mna[264](264) = mna[264][264] + 1.;
	rval[264] += Ieq_M;
	fvec[264] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM224*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[23];
	Un_M += Vx[10];
	Uo_M += -uxtold[23];
	Uo_M += uxtold[10];
	Qo_M = uxtold[265];
	Qn_M = Vx[265];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[10](10) = mna[10][10] + idt*Geq_M;
	mna[265](10) = mna[265][10] -Geq_M;
	rval[10] += -(Ieq_M - Qo_M)*idt;
	fvec[10] += Sr_M;
	mna[10](23) = mna[10][23] - idt*Geq_M;
	mna[23](10) = mna[23][10] - idt*Geq_M;
	mna[23](23) = mna[23][23] + idt*Geq_M;
	mna[265](23) = mna[265][23] +Geq_M;
	rval[23] += (Ieq_M - Qo_M)*idt;
	fvec[23] += -(Sr_M);
	mna[265](265) = mna[265][265] + 1.;
	rval[265] += Ieq_M;
	fvec[265] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM225*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[37];
	Un_M += Vx[10];
	Uo_M += -uxtold[37];
	Uo_M += uxtold[10];
	Qo_M = uxtold[266];
	Qn_M = Vx[266];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[10](10) = mna[10][10] + idt*Geq_M;
	mna[266](10) = mna[266][10] -Geq_M;
	rval[10] += -(Ieq_M - Qo_M)*idt;
	fvec[10] += Sr_M;
	mna[10](37) = mna[10][37] - idt*Geq_M;
	mna[37](10) = mna[37][10] - idt*Geq_M;
	mna[37](37) = mna[37][37] + idt*Geq_M;
	mna[266](37) = mna[266][37] +Geq_M;
	rval[37] += (Ieq_M - Qo_M)*idt;
	fvec[37] += -(Sr_M);
	mna[266](266) = mna[266][266] + 1.;
	rval[266] += Ieq_M;
	fvec[266] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM226*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[12];
	Un_M += Vx[11];
	Uo_M += -uxtold[12];
	Uo_M += uxtold[11];
	Qo_M = uxtold[267];
	Qn_M = Vx[267];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[11](11) = mna[11][11] + idt*Geq_M;
	mna[267](11) = mna[267][11] -Geq_M;
	rval[11] += -(Ieq_M - Qo_M)*idt;
	fvec[11] += Sr_M;
	mna[11](12) = mna[11][12] - idt*Geq_M;
	mna[12](11) = mna[12][11] - idt*Geq_M;
	mna[12](12) = mna[12][12] + idt*Geq_M;
	mna[267](12) = mna[267][12] +Geq_M;
	rval[12] += (Ieq_M - Qo_M)*idt;
	fvec[12] += -(Sr_M);
	mna[267](267) = mna[267][267] + 1.;
	rval[267] += Ieq_M;
	fvec[267] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM227*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[13];
	Un_M += Vx[11];
	Uo_M += -uxtold[13];
	Uo_M += uxtold[11];
	Qo_M = uxtold[268];
	Qn_M = Vx[268];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[11](11) = mna[11][11] + idt*Geq_M;
	mna[268](11) = mna[268][11] -Geq_M;
	rval[11] += -(Ieq_M - Qo_M)*idt;
	fvec[11] += Sr_M;
	mna[11](13) = mna[11][13] - idt*Geq_M;
	mna[13](11) = mna[13][11] - idt*Geq_M;
	mna[13](13) = mna[13][13] + idt*Geq_M;
	mna[268](13) = mna[268][13] +Geq_M;
	rval[13] += (Ieq_M - Qo_M)*idt;
	fvec[13] += -(Sr_M);
	mna[268](268) = mna[268][268] + 1.;
	rval[268] += Ieq_M;
	fvec[268] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM228*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[14];
	Un_M += Vx[11];
	Uo_M += -uxtold[14];
	Uo_M += uxtold[11];
	Qo_M = uxtold[269];
	Qn_M = Vx[269];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[11](11) = mna[11][11] + idt*Geq_M;
	mna[269](11) = mna[269][11] -Geq_M;
	rval[11] += -(Ieq_M - Qo_M)*idt;
	fvec[11] += Sr_M;
	mna[11](14) = mna[11][14] - idt*Geq_M;
	mna[14](11) = mna[14][11] - idt*Geq_M;
	mna[14](14) = mna[14][14] + idt*Geq_M;
	mna[269](14) = mna[269][14] +Geq_M;
	rval[14] += (Ieq_M - Qo_M)*idt;
	fvec[14] += -(Sr_M);
	mna[269](269) = mna[269][269] + 1.;
	rval[269] += Ieq_M;
	fvec[269] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM229*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[21];
	Un_M += Vx[11];
	Uo_M += -uxtold[21];
	Uo_M += uxtold[11];
	Qo_M = uxtold[270];
	Qn_M = Vx[270];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[11](11) = mna[11][11] + idt*Geq_M;
	mna[270](11) = mna[270][11] -Geq_M;
	rval[11] += -(Ieq_M - Qo_M)*idt;
	fvec[11] += Sr_M;
	mna[11](21) = mna[11][21] - idt*Geq_M;
	mna[21](11) = mna[21][11] - idt*Geq_M;
	mna[21](21) = mna[21][21] + idt*Geq_M;
	mna[270](21) = mna[270][21] +Geq_M;
	rval[21] += (Ieq_M - Qo_M)*idt;
	fvec[21] += -(Sr_M);
	mna[270](270) = mna[270][270] + 1.;
	rval[270] += Ieq_M;
	fvec[270] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM230*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[25];
	Un_M += Vx[11];
	Uo_M += -uxtold[25];
	Uo_M += uxtold[11];
	Qo_M = uxtold[271];
	Qn_M = Vx[271];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[11](11) = mna[11][11] + idt*Geq_M;
	mna[271](11) = mna[271][11] -Geq_M;
	rval[11] += -(Ieq_M - Qo_M)*idt;
	fvec[11] += Sr_M;
	mna[11](25) = mna[11][25] - idt*Geq_M;
	mna[25](11) = mna[25][11] - idt*Geq_M;
	mna[25](25) = mna[25][25] + idt*Geq_M;
	mna[271](25) = mna[271][25] +Geq_M;
	rval[25] += (Ieq_M - Qo_M)*idt;
	fvec[25] += -(Sr_M);
	mna[271](271) = mna[271][271] + 1.;
	rval[271] += Ieq_M;
	fvec[271] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM231*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[29];
	Un_M += Vx[11];
	Uo_M += -uxtold[29];
	Uo_M += uxtold[11];
	Qo_M = uxtold[272];
	Qn_M = Vx[272];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[11](11) = mna[11][11] + idt*Geq_M;
	mna[272](11) = mna[272][11] -Geq_M;
	rval[11] += -(Ieq_M - Qo_M)*idt;
	fvec[11] += Sr_M;
	mna[11](29) = mna[11][29] - idt*Geq_M;
	mna[29](11) = mna[29][11] - idt*Geq_M;
	mna[29](29) = mna[29][29] + idt*Geq_M;
	mna[272](29) = mna[272][29] +Geq_M;
	rval[29] += (Ieq_M - Qo_M)*idt;
	fvec[29] += -(Sr_M);
	mna[272](272) = mna[272][272] + 1.;
	rval[272] += Ieq_M;
	fvec[272] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM232*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[30];
	Un_M += Vx[11];
	Uo_M += -uxtold[30];
	Uo_M += uxtold[11];
	Qo_M = uxtold[273];
	Qn_M = Vx[273];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[11](11) = mna[11][11] + idt*Geq_M;
	mna[273](11) = mna[273][11] -Geq_M;
	rval[11] += -(Ieq_M - Qo_M)*idt;
	fvec[11] += Sr_M;
	mna[11](30) = mna[11][30] - idt*Geq_M;
	mna[30](11) = mna[30][11] - idt*Geq_M;
	mna[30](30) = mna[30][30] + idt*Geq_M;
	mna[273](30) = mna[273][30] +Geq_M;
	rval[30] += (Ieq_M - Qo_M)*idt;
	fvec[30] += -(Sr_M);
	mna[273](273) = mna[273][273] + 1.;
	rval[273] += Ieq_M;
	fvec[273] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM233*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[35];
	Un_M += Vx[11];
	Uo_M += -uxtold[35];
	Uo_M += uxtold[11];
	Qo_M = uxtold[274];
	Qn_M = Vx[274];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[11](11) = mna[11][11] + idt*Geq_M;
	mna[274](11) = mna[274][11] -Geq_M;
	rval[11] += -(Ieq_M - Qo_M)*idt;
	fvec[11] += Sr_M;
	mna[11](35) = mna[11][35] - idt*Geq_M;
	mna[35](11) = mna[35][11] - idt*Geq_M;
	mna[35](35) = mna[35][35] + idt*Geq_M;
	mna[274](35) = mna[274][35] +Geq_M;
	rval[35] += (Ieq_M - Qo_M)*idt;
	fvec[35] += -(Sr_M);
	mna[274](274) = mna[274][274] + 1.;
	rval[274] += Ieq_M;
	fvec[274] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM234*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[36];
	Un_M += Vx[11];
	Uo_M += -uxtold[36];
	Uo_M += uxtold[11];
	Qo_M = uxtold[275];
	Qn_M = Vx[275];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[11](11) = mna[11][11] + idt*Geq_M;
	mna[275](11) = mna[275][11] -Geq_M;
	rval[11] += -(Ieq_M - Qo_M)*idt;
	fvec[11] += Sr_M;
	mna[11](36) = mna[11][36] - idt*Geq_M;
	mna[36](11) = mna[36][11] - idt*Geq_M;
	mna[36](36) = mna[36][36] + idt*Geq_M;
	mna[275](36) = mna[275][36] +Geq_M;
	rval[36] += (Ieq_M - Qo_M)*idt;
	fvec[36] += -(Sr_M);
	mna[275](275) = mna[275][275] + 1.;
	rval[275] += Ieq_M;
	fvec[275] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM235*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[15];
	Un_M += Vx[12];
	Uo_M += -uxtold[15];
	Uo_M += uxtold[12];
	Qo_M = uxtold[276];
	Qn_M = Vx[276];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[12](12) = mna[12][12] + idt*Geq_M;
	mna[276](12) = mna[276][12] -Geq_M;
	rval[12] += -(Ieq_M - Qo_M)*idt;
	fvec[12] += Sr_M;
	mna[12](15) = mna[12][15] - idt*Geq_M;
	mna[15](12) = mna[15][12] - idt*Geq_M;
	mna[15](15) = mna[15][15] + idt*Geq_M;
	mna[276](15) = mna[276][15] +Geq_M;
	rval[15] += (Ieq_M - Qo_M)*idt;
	fvec[15] += -(Sr_M);
	mna[276](276) = mna[276][276] + 1.;
	rval[276] += Ieq_M;
	fvec[276] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM236*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[19];
	Un_M += Vx[12];
	Uo_M += -uxtold[19];
	Uo_M += uxtold[12];
	Qo_M = uxtold[277];
	Qn_M = Vx[277];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[12](12) = mna[12][12] + idt*Geq_M;
	mna[277](12) = mna[277][12] -Geq_M;
	rval[12] += -(Ieq_M - Qo_M)*idt;
	fvec[12] += Sr_M;
	mna[12](19) = mna[12][19] - idt*Geq_M;
	mna[19](12) = mna[19][12] - idt*Geq_M;
	mna[19](19) = mna[19][19] + idt*Geq_M;
	mna[277](19) = mna[277][19] +Geq_M;
	rval[19] += (Ieq_M - Qo_M)*idt;
	fvec[19] += -(Sr_M);
	mna[277](277) = mna[277][277] + 1.;
	rval[277] += Ieq_M;
	fvec[277] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM237*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[26];
	Un_M += Vx[12];
	Uo_M += -uxtold[26];
	Uo_M += uxtold[12];
	Qo_M = uxtold[278];
	Qn_M = Vx[278];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[12](12) = mna[12][12] + idt*Geq_M;
	mna[278](12) = mna[278][12] -Geq_M;
	rval[12] += -(Ieq_M - Qo_M)*idt;
	fvec[12] += Sr_M;
	mna[12](26) = mna[12][26] - idt*Geq_M;
	mna[26](12) = mna[26][12] - idt*Geq_M;
	mna[26](26) = mna[26][26] + idt*Geq_M;
	mna[278](26) = mna[278][26] +Geq_M;
	rval[26] += (Ieq_M - Qo_M)*idt;
	fvec[26] += -(Sr_M);
	mna[278](278) = mna[278][278] + 1.;
	rval[278] += Ieq_M;
	fvec[278] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM238*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[29];
	Un_M += Vx[12];
	Uo_M += -uxtold[29];
	Uo_M += uxtold[12];
	Qo_M = uxtold[279];
	Qn_M = Vx[279];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[12](12) = mna[12][12] + idt*Geq_M;
	mna[279](12) = mna[279][12] -Geq_M;
	rval[12] += -(Ieq_M - Qo_M)*idt;
	fvec[12] += Sr_M;
	mna[12](29) = mna[12][29] - idt*Geq_M;
	mna[29](12) = mna[29][12] - idt*Geq_M;
	mna[29](29) = mna[29][29] + idt*Geq_M;
	mna[279](29) = mna[279][29] +Geq_M;
	rval[29] += (Ieq_M - Qo_M)*idt;
	fvec[29] += -(Sr_M);
	mna[279](279) = mna[279][279] + 1.;
	rval[279] += Ieq_M;
	fvec[279] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM239*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[30];
	Un_M += Vx[12];
	Uo_M += -uxtold[30];
	Uo_M += uxtold[12];
	Qo_M = uxtold[280];
	Qn_M = Vx[280];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[12](12) = mna[12][12] + idt*Geq_M;
	mna[280](12) = mna[280][12] -Geq_M;
	rval[12] += -(Ieq_M - Qo_M)*idt;
	fvec[12] += Sr_M;
	mna[12](30) = mna[12][30] - idt*Geq_M;
	mna[30](12) = mna[30][12] - idt*Geq_M;
	mna[30](30) = mna[30][30] + idt*Geq_M;
	mna[280](30) = mna[280][30] +Geq_M;
	rval[30] += (Ieq_M - Qo_M)*idt;
	fvec[30] += -(Sr_M);
	mna[280](280) = mna[280][280] + 1.;
	rval[280] += Ieq_M;
	fvec[280] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM240*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[33];
	Un_M += Vx[12];
	Uo_M += -uxtold[33];
	Uo_M += uxtold[12];
	Qo_M = uxtold[281];
	Qn_M = Vx[281];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[12](12) = mna[12][12] + idt*Geq_M;
	mna[281](12) = mna[281][12] -Geq_M;
	rval[12] += -(Ieq_M - Qo_M)*idt;
	fvec[12] += Sr_M;
	mna[12](33) = mna[12][33] - idt*Geq_M;
	mna[33](12) = mna[33][12] - idt*Geq_M;
	mna[33](33) = mna[33][33] + idt*Geq_M;
	mna[281](33) = mna[281][33] +Geq_M;
	rval[33] += (Ieq_M - Qo_M)*idt;
	fvec[33] += -(Sr_M);
	mna[281](281) = mna[281][281] + 1.;
	rval[281] += Ieq_M;
	fvec[281] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM241*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[36];
	Un_M += Vx[12];
	Uo_M += -uxtold[36];
	Uo_M += uxtold[12];
	Qo_M = uxtold[282];
	Qn_M = Vx[282];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[12](12) = mna[12][12] + idt*Geq_M;
	mna[282](12) = mna[282][12] -Geq_M;
	rval[12] += -(Ieq_M - Qo_M)*idt;
	fvec[12] += Sr_M;
	mna[12](36) = mna[12][36] - idt*Geq_M;
	mna[36](12) = mna[36][12] - idt*Geq_M;
	mna[36](36) = mna[36][36] + idt*Geq_M;
	mna[282](36) = mna[282][36] +Geq_M;
	rval[36] += (Ieq_M - Qo_M)*idt;
	fvec[36] += -(Sr_M);
	mna[282](282) = mna[282][282] + 1.;
	rval[282] += Ieq_M;
	fvec[282] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM242*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[14];
	Un_M += Vx[13];
	Uo_M += -uxtold[14];
	Uo_M += uxtold[13];
	Qo_M = uxtold[283];
	Qn_M = Vx[283];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[13](13) = mna[13][13] + idt*Geq_M;
	mna[283](13) = mna[283][13] -Geq_M;
	rval[13] += -(Ieq_M - Qo_M)*idt;
	fvec[13] += Sr_M;
	mna[13](14) = mna[13][14] - idt*Geq_M;
	mna[14](13) = mna[14][13] - idt*Geq_M;
	mna[14](14) = mna[14][14] + idt*Geq_M;
	mna[283](14) = mna[283][14] +Geq_M;
	rval[14] += (Ieq_M - Qo_M)*idt;
	fvec[14] += -(Sr_M);
	mna[283](283) = mna[283][283] + 1.;
	rval[283] += Ieq_M;
	fvec[283] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM243*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[15];
	Un_M += Vx[13];
	Uo_M += -uxtold[15];
	Uo_M += uxtold[13];
	Qo_M = uxtold[284];
	Qn_M = Vx[284];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[13](13) = mna[13][13] + idt*Geq_M;
	mna[284](13) = mna[284][13] -Geq_M;
	rval[13] += -(Ieq_M - Qo_M)*idt;
	fvec[13] += Sr_M;
	mna[13](15) = mna[13][15] - idt*Geq_M;
	mna[15](13) = mna[15][13] - idt*Geq_M;
	mna[15](15) = mna[15][15] + idt*Geq_M;
	mna[284](15) = mna[284][15] +Geq_M;
	rval[15] += (Ieq_M - Qo_M)*idt;
	fvec[15] += -(Sr_M);
	mna[284](284) = mna[284][284] + 1.;
	rval[284] += Ieq_M;
	fvec[284] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM244*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[16];
	Un_M += Vx[13];
	Uo_M += -uxtold[16];
	Uo_M += uxtold[13];
	Qo_M = uxtold[285];
	Qn_M = Vx[285];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[13](13) = mna[13][13] + idt*Geq_M;
	mna[285](13) = mna[285][13] -Geq_M;
	rval[13] += -(Ieq_M - Qo_M)*idt;
	fvec[13] += Sr_M;
	mna[13](16) = mna[13][16] - idt*Geq_M;
	mna[16](13) = mna[16][13] - idt*Geq_M;
	mna[16](16) = mna[16][16] + idt*Geq_M;
	mna[285](16) = mna[285][16] +Geq_M;
	rval[16] += (Ieq_M - Qo_M)*idt;
	fvec[16] += -(Sr_M);
	mna[285](285) = mna[285][285] + 1.;
	rval[285] += Ieq_M;
	fvec[285] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM245*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[20];
	Un_M += Vx[13];
	Uo_M += -uxtold[20];
	Uo_M += uxtold[13];
	Qo_M = uxtold[286];
	Qn_M = Vx[286];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[13](13) = mna[13][13] + idt*Geq_M;
	mna[286](13) = mna[286][13] -Geq_M;
	rval[13] += -(Ieq_M - Qo_M)*idt;
	fvec[13] += Sr_M;
	mna[13](20) = mna[13][20] - idt*Geq_M;
	mna[20](13) = mna[20][13] - idt*Geq_M;
	mna[20](20) = mna[20][20] + idt*Geq_M;
	mna[286](20) = mna[286][20] +Geq_M;
	rval[20] += (Ieq_M - Qo_M)*idt;
	fvec[20] += -(Sr_M);
	mna[286](286) = mna[286][286] + 1.;
	rval[286] += Ieq_M;
	fvec[286] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM246*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[24];
	Un_M += Vx[13];
	Uo_M += -uxtold[24];
	Uo_M += uxtold[13];
	Qo_M = uxtold[287];
	Qn_M = Vx[287];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[13](13) = mna[13][13] + idt*Geq_M;
	mna[287](13) = mna[287][13] -Geq_M;
	rval[13] += -(Ieq_M - Qo_M)*idt;
	fvec[13] += Sr_M;
	mna[13](24) = mna[13][24] - idt*Geq_M;
	mna[24](13) = mna[24][13] - idt*Geq_M;
	mna[24](24) = mna[24][24] + idt*Geq_M;
	mna[287](24) = mna[287][24] +Geq_M;
	rval[24] += (Ieq_M - Qo_M)*idt;
	fvec[24] += -(Sr_M);
	mna[287](287) = mna[287][287] + 1.;
	rval[287] += Ieq_M;
	fvec[287] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM247*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[28];
	Un_M += Vx[13];
	Uo_M += -uxtold[28];
	Uo_M += uxtold[13];
	Qo_M = uxtold[288];
	Qn_M = Vx[288];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[13](13) = mna[13][13] + idt*Geq_M;
	mna[288](13) = mna[288][13] -Geq_M;
	rval[13] += -(Ieq_M - Qo_M)*idt;
	fvec[13] += Sr_M;
	mna[13](28) = mna[13][28] - idt*Geq_M;
	mna[28](13) = mna[28][13] - idt*Geq_M;
	mna[28](28) = mna[28][28] + idt*Geq_M;
	mna[288](28) = mna[288][28] +Geq_M;
	rval[28] += (Ieq_M - Qo_M)*idt;
	fvec[28] += -(Sr_M);
	mna[288](288) = mna[288][288] + 1.;
	rval[288] += Ieq_M;
	fvec[288] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM248*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[34];
	Un_M += Vx[13];
	Uo_M += -uxtold[34];
	Uo_M += uxtold[13];
	Qo_M = uxtold[289];
	Qn_M = Vx[289];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[13](13) = mna[13][13] + idt*Geq_M;
	mna[289](13) = mna[289][13] -Geq_M;
	rval[13] += -(Ieq_M - Qo_M)*idt;
	fvec[13] += Sr_M;
	mna[13](34) = mna[13][34] - idt*Geq_M;
	mna[34](13) = mna[34][13] - idt*Geq_M;
	mna[34](34) = mna[34][34] + idt*Geq_M;
	mna[289](34) = mna[289][34] +Geq_M;
	rval[34] += (Ieq_M - Qo_M)*idt;
	fvec[34] += -(Sr_M);
	mna[289](289) = mna[289][289] + 1.;
	rval[289] += Ieq_M;
	fvec[289] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM249*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[36];
	Un_M += Vx[13];
	Uo_M += -uxtold[36];
	Uo_M += uxtold[13];
	Qo_M = uxtold[290];
	Qn_M = Vx[290];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[13](13) = mna[13][13] + idt*Geq_M;
	mna[290](13) = mna[290][13] -Geq_M;
	rval[13] += -(Ieq_M - Qo_M)*idt;
	fvec[13] += Sr_M;
	mna[13](36) = mna[13][36] - idt*Geq_M;
	mna[36](13) = mna[36][13] - idt*Geq_M;
	mna[36](36) = mna[36][36] + idt*Geq_M;
	mna[290](36) = mna[290][36] +Geq_M;
	rval[36] += (Ieq_M - Qo_M)*idt;
	fvec[36] += -(Sr_M);
	mna[290](290) = mna[290][290] + 1.;
	rval[290] += Ieq_M;
	fvec[290] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM250*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[26];
	Un_M += Vx[14];
	Uo_M += -uxtold[26];
	Uo_M += uxtold[14];
	Qo_M = uxtold[291];
	Qn_M = Vx[291];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[14](14) = mna[14][14] + idt*Geq_M;
	mna[291](14) = mna[291][14] -Geq_M;
	rval[14] += -(Ieq_M - Qo_M)*idt;
	fvec[14] += Sr_M;
	mna[14](26) = mna[14][26] - idt*Geq_M;
	mna[26](14) = mna[26][14] - idt*Geq_M;
	mna[26](26) = mna[26][26] + idt*Geq_M;
	mna[291](26) = mna[291][26] +Geq_M;
	rval[26] += (Ieq_M - Qo_M)*idt;
	fvec[26] += -(Sr_M);
	mna[291](291) = mna[291][291] + 1.;
	rval[291] += Ieq_M;
	fvec[291] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM251*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[27];
	Un_M += Vx[14];
	Uo_M += -uxtold[27];
	Uo_M += uxtold[14];
	Qo_M = uxtold[292];
	Qn_M = Vx[292];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[14](14) = mna[14][14] + idt*Geq_M;
	mna[292](14) = mna[292][14] -Geq_M;
	rval[14] += -(Ieq_M - Qo_M)*idt;
	fvec[14] += Sr_M;
	mna[14](27) = mna[14][27] - idt*Geq_M;
	mna[27](14) = mna[27][14] - idt*Geq_M;
	mna[27](27) = mna[27][27] + idt*Geq_M;
	mna[292](27) = mna[292][27] +Geq_M;
	rval[27] += (Ieq_M - Qo_M)*idt;
	fvec[27] += -(Sr_M);
	mna[292](292) = mna[292][292] + 1.;
	rval[292] += Ieq_M;
	fvec[292] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM252*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[29];
	Un_M += Vx[14];
	Uo_M += -uxtold[29];
	Uo_M += uxtold[14];
	Qo_M = uxtold[293];
	Qn_M = Vx[293];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[14](14) = mna[14][14] + idt*Geq_M;
	mna[293](14) = mna[293][14] -Geq_M;
	rval[14] += -(Ieq_M - Qo_M)*idt;
	fvec[14] += Sr_M;
	mna[14](29) = mna[14][29] - idt*Geq_M;
	mna[29](14) = mna[29][14] - idt*Geq_M;
	mna[29](29) = mna[29][29] + idt*Geq_M;
	mna[293](29) = mna[293][29] +Geq_M;
	rval[29] += (Ieq_M - Qo_M)*idt;
	fvec[29] += -(Sr_M);
	mna[293](293) = mna[293][293] + 1.;
	rval[293] += Ieq_M;
	fvec[293] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM253*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[32];
	Un_M += Vx[14];
	Uo_M += -uxtold[32];
	Uo_M += uxtold[14];
	Qo_M = uxtold[294];
	Qn_M = Vx[294];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[14](14) = mna[14][14] + idt*Geq_M;
	mna[294](14) = mna[294][14] -Geq_M;
	rval[14] += -(Ieq_M - Qo_M)*idt;
	fvec[14] += Sr_M;
	mna[14](32) = mna[14][32] - idt*Geq_M;
	mna[32](14) = mna[32][14] - idt*Geq_M;
	mna[32](32) = mna[32][32] + idt*Geq_M;
	mna[294](32) = mna[294][32] +Geq_M;
	rval[32] += (Ieq_M - Qo_M)*idt;
	fvec[32] += -(Sr_M);
	mna[294](294) = mna[294][294] + 1.;
	rval[294] += Ieq_M;
	fvec[294] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM254*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[33];
	Un_M += Vx[14];
	Uo_M += -uxtold[33];
	Uo_M += uxtold[14];
	Qo_M = uxtold[295];
	Qn_M = Vx[295];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[14](14) = mna[14][14] + idt*Geq_M;
	mna[295](14) = mna[295][14] -Geq_M;
	rval[14] += -(Ieq_M - Qo_M)*idt;
	fvec[14] += Sr_M;
	mna[14](33) = mna[14][33] - idt*Geq_M;
	mna[33](14) = mna[33][14] - idt*Geq_M;
	mna[33](33) = mna[33][33] + idt*Geq_M;
	mna[295](33) = mna[295][33] +Geq_M;
	rval[33] += (Ieq_M - Qo_M)*idt;
	fvec[33] += -(Sr_M);
	mna[295](295) = mna[295][295] + 1.;
	rval[295] += Ieq_M;
	fvec[295] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM255*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[39];
	Un_M += Vx[14];
	Uo_M += -uxtold[39];
	Uo_M += uxtold[14];
	Qo_M = uxtold[296];
	Qn_M = Vx[296];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[14](14) = mna[14][14] + idt*Geq_M;
	mna[296](14) = mna[296][14] -Geq_M;
	rval[14] += -(Ieq_M - Qo_M)*idt;
	fvec[14] += Sr_M;
	mna[14](39) = mna[14][39] - idt*Geq_M;
	mna[39](14) = mna[39][14] - idt*Geq_M;
	mna[39](39) = mna[39][39] + idt*Geq_M;
	mna[296](39) = mna[296][39] +Geq_M;
	rval[39] += (Ieq_M - Qo_M)*idt;
	fvec[39] += -(Sr_M);
	mna[296](296) = mna[296][296] + 1.;
	rval[296] += Ieq_M;
	fvec[296] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM256*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[17];
	Un_M += Vx[15];
	Uo_M += -uxtold[17];
	Uo_M += uxtold[15];
	Qo_M = uxtold[297];
	Qn_M = Vx[297];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[15](15) = mna[15][15] + idt*Geq_M;
	mna[297](15) = mna[297][15] -Geq_M;
	rval[15] += -(Ieq_M - Qo_M)*idt;
	fvec[15] += Sr_M;
	mna[15](17) = mna[15][17] - idt*Geq_M;
	mna[17](15) = mna[17][15] - idt*Geq_M;
	mna[17](17) = mna[17][17] + idt*Geq_M;
	mna[297](17) = mna[297][17] +Geq_M;
	rval[17] += (Ieq_M - Qo_M)*idt;
	fvec[17] += -(Sr_M);
	mna[297](297) = mna[297][297] + 1.;
	rval[297] += Ieq_M;
	fvec[297] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM257*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[20];
	Un_M += Vx[15];
	Uo_M += -uxtold[20];
	Uo_M += uxtold[15];
	Qo_M = uxtold[298];
	Qn_M = Vx[298];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[15](15) = mna[15][15] + idt*Geq_M;
	mna[298](15) = mna[298][15] -Geq_M;
	rval[15] += -(Ieq_M - Qo_M)*idt;
	fvec[15] += Sr_M;
	mna[15](20) = mna[15][20] - idt*Geq_M;
	mna[20](15) = mna[20][15] - idt*Geq_M;
	mna[20](20) = mna[20][20] + idt*Geq_M;
	mna[298](20) = mna[298][20] +Geq_M;
	rval[20] += (Ieq_M - Qo_M)*idt;
	fvec[20] += -(Sr_M);
	mna[298](298) = mna[298][298] + 1.;
	rval[298] += Ieq_M;
	fvec[298] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM258*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[26];
	Un_M += Vx[15];
	Uo_M += -uxtold[26];
	Uo_M += uxtold[15];
	Qo_M = uxtold[299];
	Qn_M = Vx[299];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[15](15) = mna[15][15] + idt*Geq_M;
	mna[299](15) = mna[299][15] -Geq_M;
	rval[15] += -(Ieq_M - Qo_M)*idt;
	fvec[15] += Sr_M;
	mna[15](26) = mna[15][26] - idt*Geq_M;
	mna[26](15) = mna[26][15] - idt*Geq_M;
	mna[26](26) = mna[26][26] + idt*Geq_M;
	mna[299](26) = mna[299][26] +Geq_M;
	rval[26] += (Ieq_M - Qo_M)*idt;
	fvec[26] += -(Sr_M);
	mna[299](299) = mna[299][299] + 1.;
	rval[299] += Ieq_M;
	fvec[299] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM259*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[27];
	Un_M += Vx[15];
	Uo_M += -uxtold[27];
	Uo_M += uxtold[15];
	Qo_M = uxtold[300];
	Qn_M = Vx[300];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[15](15) = mna[15][15] + idt*Geq_M;
	mna[300](15) = mna[300][15] -Geq_M;
	rval[15] += -(Ieq_M - Qo_M)*idt;
	fvec[15] += Sr_M;
	mna[15](27) = mna[15][27] - idt*Geq_M;
	mna[27](15) = mna[27][15] - idt*Geq_M;
	mna[27](27) = mna[27][27] + idt*Geq_M;
	mna[300](27) = mna[300][27] +Geq_M;
	rval[27] += (Ieq_M - Qo_M)*idt;
	fvec[27] += -(Sr_M);
	mna[300](300) = mna[300][300] + 1.;
	rval[300] += Ieq_M;
	fvec[300] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM260*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[29];
	Un_M += Vx[15];
	Uo_M += -uxtold[29];
	Uo_M += uxtold[15];
	Qo_M = uxtold[301];
	Qn_M = Vx[301];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[15](15) = mna[15][15] + idt*Geq_M;
	mna[301](15) = mna[301][15] -Geq_M;
	rval[15] += -(Ieq_M - Qo_M)*idt;
	fvec[15] += Sr_M;
	mna[15](29) = mna[15][29] - idt*Geq_M;
	mna[29](15) = mna[29][15] - idt*Geq_M;
	mna[29](29) = mna[29][29] + idt*Geq_M;
	mna[301](29) = mna[301][29] +Geq_M;
	rval[29] += (Ieq_M - Qo_M)*idt;
	fvec[29] += -(Sr_M);
	mna[301](301) = mna[301][301] + 1.;
	rval[301] += Ieq_M;
	fvec[301] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM261*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[31];
	Un_M += Vx[15];
	Uo_M += -uxtold[31];
	Uo_M += uxtold[15];
	Qo_M = uxtold[302];
	Qn_M = Vx[302];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[15](15) = mna[15][15] + idt*Geq_M;
	mna[302](15) = mna[302][15] -Geq_M;
	rval[15] += -(Ieq_M - Qo_M)*idt;
	fvec[15] += Sr_M;
	mna[15](31) = mna[15][31] - idt*Geq_M;
	mna[31](15) = mna[31][15] - idt*Geq_M;
	mna[31](31) = mna[31][31] + idt*Geq_M;
	mna[302](31) = mna[302][31] +Geq_M;
	rval[31] += (Ieq_M - Qo_M)*idt;
	fvec[31] += -(Sr_M);
	mna[302](302) = mna[302][302] + 1.;
	rval[302] += Ieq_M;
	fvec[302] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM262*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[17];
	Un_M += Vx[16];
	Uo_M += -uxtold[17];
	Uo_M += uxtold[16];
	Qo_M = uxtold[303];
	Qn_M = Vx[303];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[16](16) = mna[16][16] + idt*Geq_M;
	mna[303](16) = mna[303][16] -Geq_M;
	rval[16] += -(Ieq_M - Qo_M)*idt;
	fvec[16] += Sr_M;
	mna[16](17) = mna[16][17] - idt*Geq_M;
	mna[17](16) = mna[17][16] - idt*Geq_M;
	mna[17](17) = mna[17][17] + idt*Geq_M;
	mna[303](17) = mna[303][17] +Geq_M;
	rval[17] += (Ieq_M - Qo_M)*idt;
	fvec[17] += -(Sr_M);
	mna[303](303) = mna[303][303] + 1.;
	rval[303] += Ieq_M;
	fvec[303] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM263*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[19];
	Un_M += Vx[16];
	Uo_M += -uxtold[19];
	Uo_M += uxtold[16];
	Qo_M = uxtold[304];
	Qn_M = Vx[304];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[16](16) = mna[16][16] + idt*Geq_M;
	mna[304](16) = mna[304][16] -Geq_M;
	rval[16] += -(Ieq_M - Qo_M)*idt;
	fvec[16] += Sr_M;
	mna[16](19) = mna[16][19] - idt*Geq_M;
	mna[19](16) = mna[19][16] - idt*Geq_M;
	mna[19](19) = mna[19][19] + idt*Geq_M;
	mna[304](19) = mna[304][19] +Geq_M;
	rval[19] += (Ieq_M - Qo_M)*idt;
	fvec[19] += -(Sr_M);
	mna[304](304) = mna[304][304] + 1.;
	rval[304] += Ieq_M;
	fvec[304] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM264*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[21];
	Un_M += Vx[16];
	Uo_M += -uxtold[21];
	Uo_M += uxtold[16];
	Qo_M = uxtold[305];
	Qn_M = Vx[305];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[16](16) = mna[16][16] + idt*Geq_M;
	mna[305](16) = mna[305][16] -Geq_M;
	rval[16] += -(Ieq_M - Qo_M)*idt;
	fvec[16] += Sr_M;
	mna[16](21) = mna[16][21] - idt*Geq_M;
	mna[21](16) = mna[21][16] - idt*Geq_M;
	mna[21](21) = mna[21][21] + idt*Geq_M;
	mna[305](21) = mna[305][21] +Geq_M;
	rval[21] += (Ieq_M - Qo_M)*idt;
	fvec[21] += -(Sr_M);
	mna[305](305) = mna[305][305] + 1.;
	rval[305] += Ieq_M;
	fvec[305] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM265*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[27];
	Un_M += Vx[16];
	Uo_M += -uxtold[27];
	Uo_M += uxtold[16];
	Qo_M = uxtold[306];
	Qn_M = Vx[306];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[16](16) = mna[16][16] + idt*Geq_M;
	mna[306](16) = mna[306][16] -Geq_M;
	rval[16] += -(Ieq_M - Qo_M)*idt;
	fvec[16] += Sr_M;
	mna[16](27) = mna[16][27] - idt*Geq_M;
	mna[27](16) = mna[27][16] - idt*Geq_M;
	mna[27](27) = mna[27][27] + idt*Geq_M;
	mna[306](27) = mna[306][27] +Geq_M;
	rval[27] += (Ieq_M - Qo_M)*idt;
	fvec[27] += -(Sr_M);
	mna[306](306) = mna[306][306] + 1.;
	rval[306] += Ieq_M;
	fvec[306] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM266*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[28];
	Un_M += Vx[16];
	Uo_M += -uxtold[28];
	Uo_M += uxtold[16];
	Qo_M = uxtold[307];
	Qn_M = Vx[307];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[16](16) = mna[16][16] + idt*Geq_M;
	mna[307](16) = mna[307][16] -Geq_M;
	rval[16] += -(Ieq_M - Qo_M)*idt;
	fvec[16] += Sr_M;
	mna[16](28) = mna[16][28] - idt*Geq_M;
	mna[28](16) = mna[28][16] - idt*Geq_M;
	mna[28](28) = mna[28][28] + idt*Geq_M;
	mna[307](28) = mna[307][28] +Geq_M;
	rval[28] += (Ieq_M - Qo_M)*idt;
	fvec[28] += -(Sr_M);
	mna[307](307) = mna[307][307] + 1.;
	rval[307] += Ieq_M;
	fvec[307] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM267*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[30];
	Un_M += Vx[16];
	Uo_M += -uxtold[30];
	Uo_M += uxtold[16];
	Qo_M = uxtold[308];
	Qn_M = Vx[308];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[16](16) = mna[16][16] + idt*Geq_M;
	mna[308](16) = mna[308][16] -Geq_M;
	rval[16] += -(Ieq_M - Qo_M)*idt;
	fvec[16] += Sr_M;
	mna[16](30) = mna[16][30] - idt*Geq_M;
	mna[30](16) = mna[30][16] - idt*Geq_M;
	mna[30](30) = mna[30][30] + idt*Geq_M;
	mna[308](30) = mna[308][30] +Geq_M;
	rval[30] += (Ieq_M - Qo_M)*idt;
	fvec[30] += -(Sr_M);
	mna[308](308) = mna[308][308] + 1.;
	rval[308] += Ieq_M;
	fvec[308] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM268*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[18];
	Un_M += Vx[17];
	Uo_M += -uxtold[18];
	Uo_M += uxtold[17];
	Qo_M = uxtold[309];
	Qn_M = Vx[309];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[17](17) = mna[17][17] + idt*Geq_M;
	mna[309](17) = mna[309][17] -Geq_M;
	rval[17] += -(Ieq_M - Qo_M)*idt;
	fvec[17] += Sr_M;
	mna[17](18) = mna[17][18] - idt*Geq_M;
	mna[18](17) = mna[18][17] - idt*Geq_M;
	mna[18](18) = mna[18][18] + idt*Geq_M;
	mna[309](18) = mna[309][18] +Geq_M;
	rval[18] += (Ieq_M - Qo_M)*idt;
	fvec[18] += -(Sr_M);
	mna[309](309) = mna[309][309] + 1.;
	rval[309] += Ieq_M;
	fvec[309] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM269*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[19];
	Un_M += Vx[17];
	Uo_M += -uxtold[19];
	Uo_M += uxtold[17];
	Qo_M = uxtold[310];
	Qn_M = Vx[310];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[17](17) = mna[17][17] + idt*Geq_M;
	mna[310](17) = mna[310][17] -Geq_M;
	rval[17] += -(Ieq_M - Qo_M)*idt;
	fvec[17] += Sr_M;
	mna[17](19) = mna[17][19] - idt*Geq_M;
	mna[19](17) = mna[19][17] - idt*Geq_M;
	mna[19](19) = mna[19][19] + idt*Geq_M;
	mna[310](19) = mna[310][19] +Geq_M;
	rval[19] += (Ieq_M - Qo_M)*idt;
	fvec[19] += -(Sr_M);
	mna[310](310) = mna[310][310] + 1.;
	rval[310] += Ieq_M;
	fvec[310] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM270*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[22];
	Un_M += Vx[17];
	Uo_M += -uxtold[22];
	Uo_M += uxtold[17];
	Qo_M = uxtold[311];
	Qn_M = Vx[311];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[17](17) = mna[17][17] + idt*Geq_M;
	mna[311](17) = mna[311][17] -Geq_M;
	rval[17] += -(Ieq_M - Qo_M)*idt;
	fvec[17] += Sr_M;
	mna[17](22) = mna[17][22] - idt*Geq_M;
	mna[22](17) = mna[22][17] - idt*Geq_M;
	mna[22](22) = mna[22][22] + idt*Geq_M;
	mna[311](22) = mna[311][22] +Geq_M;
	rval[22] += (Ieq_M - Qo_M)*idt;
	fvec[22] += -(Sr_M);
	mna[311](311) = mna[311][311] + 1.;
	rval[311] += Ieq_M;
	fvec[311] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM271*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[24];
	Un_M += Vx[17];
	Uo_M += -uxtold[24];
	Uo_M += uxtold[17];
	Qo_M = uxtold[312];
	Qn_M = Vx[312];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[17](17) = mna[17][17] + idt*Geq_M;
	mna[312](17) = mna[312][17] -Geq_M;
	rval[17] += -(Ieq_M - Qo_M)*idt;
	fvec[17] += Sr_M;
	mna[17](24) = mna[17][24] - idt*Geq_M;
	mna[24](17) = mna[24][17] - idt*Geq_M;
	mna[24](24) = mna[24][24] + idt*Geq_M;
	mna[312](24) = mna[312][24] +Geq_M;
	rval[24] += (Ieq_M - Qo_M)*idt;
	fvec[24] += -(Sr_M);
	mna[312](312) = mna[312][312] + 1.;
	rval[312] += Ieq_M;
	fvec[312] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM272*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[28];
	Un_M += Vx[17];
	Uo_M += -uxtold[28];
	Uo_M += uxtold[17];
	Qo_M = uxtold[313];
	Qn_M = Vx[313];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[17](17) = mna[17][17] + idt*Geq_M;
	mna[313](17) = mna[313][17] -Geq_M;
	rval[17] += -(Ieq_M - Qo_M)*idt;
	fvec[17] += Sr_M;
	mna[17](28) = mna[17][28] - idt*Geq_M;
	mna[28](17) = mna[28][17] - idt*Geq_M;
	mna[28](28) = mna[28][28] + idt*Geq_M;
	mna[313](28) = mna[313][28] +Geq_M;
	rval[28] += (Ieq_M - Qo_M)*idt;
	fvec[28] += -(Sr_M);
	mna[313](313) = mna[313][313] + 1.;
	rval[313] += Ieq_M;
	fvec[313] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM273*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[30];
	Un_M += Vx[17];
	Uo_M += -uxtold[30];
	Uo_M += uxtold[17];
	Qo_M = uxtold[314];
	Qn_M = Vx[314];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[17](17) = mna[17][17] + idt*Geq_M;
	mna[314](17) = mna[314][17] -Geq_M;
	rval[17] += -(Ieq_M - Qo_M)*idt;
	fvec[17] += Sr_M;
	mna[17](30) = mna[17][30] - idt*Geq_M;
	mna[30](17) = mna[30][17] - idt*Geq_M;
	mna[30](30) = mna[30][30] + idt*Geq_M;
	mna[314](30) = mna[314][30] +Geq_M;
	rval[30] += (Ieq_M - Qo_M)*idt;
	fvec[30] += -(Sr_M);
	mna[314](314) = mna[314][314] + 1.;
	rval[314] += Ieq_M;
	fvec[314] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM274*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[37];
	Un_M += Vx[17];
	Uo_M += -uxtold[37];
	Uo_M += uxtold[17];
	Qo_M = uxtold[315];
	Qn_M = Vx[315];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[17](17) = mna[17][17] + idt*Geq_M;
	mna[315](17) = mna[315][17] -Geq_M;
	rval[17] += -(Ieq_M - Qo_M)*idt;
	fvec[17] += Sr_M;
	mna[17](37) = mna[17][37] - idt*Geq_M;
	mna[37](17) = mna[37][17] - idt*Geq_M;
	mna[37](37) = mna[37][37] + idt*Geq_M;
	mna[315](37) = mna[315][37] +Geq_M;
	rval[37] += (Ieq_M - Qo_M)*idt;
	fvec[37] += -(Sr_M);
	mna[315](315) = mna[315][315] + 1.;
	rval[315] += Ieq_M;
	fvec[315] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM275*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[40];
	Un_M += Vx[17];
	Uo_M += -uxtold[40];
	Uo_M += uxtold[17];
	Qo_M = uxtold[316];
	Qn_M = Vx[316];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[17](17) = mna[17][17] + idt*Geq_M;
	mna[316](17) = mna[316][17] -Geq_M;
	rval[17] += -(Ieq_M - Qo_M)*idt;
	fvec[17] += Sr_M;
	mna[17](40) = mna[17][40] - idt*Geq_M;
	mna[40](17) = mna[40][17] - idt*Geq_M;
	mna[40](40) = mna[40][40] + idt*Geq_M;
	mna[316](40) = mna[316][40] +Geq_M;
	rval[40] += (Ieq_M - Qo_M)*idt;
	fvec[40] += -(Sr_M);
	mna[316](316) = mna[316][316] + 1.;
	rval[316] += Ieq_M;
	fvec[316] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM276*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[24];
	Un_M += Vx[18];
	Uo_M += -uxtold[24];
	Uo_M += uxtold[18];
	Qo_M = uxtold[317];
	Qn_M = Vx[317];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[18](18) = mna[18][18] + idt*Geq_M;
	mna[317](18) = mna[317][18] -Geq_M;
	rval[18] += -(Ieq_M - Qo_M)*idt;
	fvec[18] += Sr_M;
	mna[18](24) = mna[18][24] - idt*Geq_M;
	mna[24](18) = mna[24][18] - idt*Geq_M;
	mna[24](24) = mna[24][24] + idt*Geq_M;
	mna[317](24) = mna[317][24] +Geq_M;
	rval[24] += (Ieq_M - Qo_M)*idt;
	fvec[24] += -(Sr_M);
	mna[317](317) = mna[317][317] + 1.;
	rval[317] += Ieq_M;
	fvec[317] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM277*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[29];
	Un_M += Vx[18];
	Uo_M += -uxtold[29];
	Uo_M += uxtold[18];
	Qo_M = uxtold[318];
	Qn_M = Vx[318];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[18](18) = mna[18][18] + idt*Geq_M;
	mna[318](18) = mna[318][18] -Geq_M;
	rval[18] += -(Ieq_M - Qo_M)*idt;
	fvec[18] += Sr_M;
	mna[18](29) = mna[18][29] - idt*Geq_M;
	mna[29](18) = mna[29][18] - idt*Geq_M;
	mna[29](29) = mna[29][29] + idt*Geq_M;
	mna[318](29) = mna[318][29] +Geq_M;
	rval[29] += (Ieq_M - Qo_M)*idt;
	fvec[29] += -(Sr_M);
	mna[318](318) = mna[318][318] + 1.;
	rval[318] += Ieq_M;
	fvec[318] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM278*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[31];
	Un_M += Vx[18];
	Uo_M += -uxtold[31];
	Uo_M += uxtold[18];
	Qo_M = uxtold[319];
	Qn_M = Vx[319];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[18](18) = mna[18][18] + idt*Geq_M;
	mna[319](18) = mna[319][18] -Geq_M;
	rval[18] += -(Ieq_M - Qo_M)*idt;
	fvec[18] += Sr_M;
	mna[18](31) = mna[18][31] - idt*Geq_M;
	mna[31](18) = mna[31][18] - idt*Geq_M;
	mna[31](31) = mna[31][31] + idt*Geq_M;
	mna[319](31) = mna[319][31] +Geq_M;
	rval[31] += (Ieq_M - Qo_M)*idt;
	fvec[31] += -(Sr_M);
	mna[319](319) = mna[319][319] + 1.;
	rval[319] += Ieq_M;
	fvec[319] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM279*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[33];
	Un_M += Vx[18];
	Uo_M += -uxtold[33];
	Uo_M += uxtold[18];
	Qo_M = uxtold[320];
	Qn_M = Vx[320];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[18](18) = mna[18][18] + idt*Geq_M;
	mna[320](18) = mna[320][18] -Geq_M;
	rval[18] += -(Ieq_M - Qo_M)*idt;
	fvec[18] += Sr_M;
	mna[18](33) = mna[18][33] - idt*Geq_M;
	mna[33](18) = mna[33][18] - idt*Geq_M;
	mna[33](33) = mna[33][33] + idt*Geq_M;
	mna[320](33) = mna[320][33] +Geq_M;
	rval[33] += (Ieq_M - Qo_M)*idt;
	fvec[33] += -(Sr_M);
	mna[320](320) = mna[320][320] + 1.;
	rval[320] += Ieq_M;
	fvec[320] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM280*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[22];
	Un_M += Vx[19];
	Uo_M += -uxtold[22];
	Uo_M += uxtold[19];
	Qo_M = uxtold[321];
	Qn_M = Vx[321];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[19](19) = mna[19][19] + idt*Geq_M;
	mna[321](19) = mna[321][19] -Geq_M;
	rval[19] += -(Ieq_M - Qo_M)*idt;
	fvec[19] += Sr_M;
	mna[19](22) = mna[19][22] - idt*Geq_M;
	mna[22](19) = mna[22][19] - idt*Geq_M;
	mna[22](22) = mna[22][22] + idt*Geq_M;
	mna[321](22) = mna[321][22] +Geq_M;
	rval[22] += (Ieq_M - Qo_M)*idt;
	fvec[22] += -(Sr_M);
	mna[321](321) = mna[321][321] + 1.;
	rval[321] += Ieq_M;
	fvec[321] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM281*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[24];
	Un_M += Vx[19];
	Uo_M += -uxtold[24];
	Uo_M += uxtold[19];
	Qo_M = uxtold[322];
	Qn_M = Vx[322];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[19](19) = mna[19][19] + idt*Geq_M;
	mna[322](19) = mna[322][19] -Geq_M;
	rval[19] += -(Ieq_M - Qo_M)*idt;
	fvec[19] += Sr_M;
	mna[19](24) = mna[19][24] - idt*Geq_M;
	mna[24](19) = mna[24][19] - idt*Geq_M;
	mna[24](24) = mna[24][24] + idt*Geq_M;
	mna[322](24) = mna[322][24] +Geq_M;
	rval[24] += (Ieq_M - Qo_M)*idt;
	fvec[24] += -(Sr_M);
	mna[322](322) = mna[322][322] + 1.;
	rval[322] += Ieq_M;
	fvec[322] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM282*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[28];
	Un_M += Vx[19];
	Uo_M += -uxtold[28];
	Uo_M += uxtold[19];
	Qo_M = uxtold[323];
	Qn_M = Vx[323];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[19](19) = mna[19][19] + idt*Geq_M;
	mna[323](19) = mna[323][19] -Geq_M;
	rval[19] += -(Ieq_M - Qo_M)*idt;
	fvec[19] += Sr_M;
	mna[19](28) = mna[19][28] - idt*Geq_M;
	mna[28](19) = mna[28][19] - idt*Geq_M;
	mna[28](28) = mna[28][28] + idt*Geq_M;
	mna[323](28) = mna[323][28] +Geq_M;
	rval[28] += (Ieq_M - Qo_M)*idt;
	fvec[28] += -(Sr_M);
	mna[323](323) = mna[323][323] + 1.;
	rval[323] += Ieq_M;
	fvec[323] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM283*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[38];
	Un_M += Vx[19];
	Uo_M += -uxtold[38];
	Uo_M += uxtold[19];
	Qo_M = uxtold[324];
	Qn_M = Vx[324];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[19](19) = mna[19][19] + idt*Geq_M;
	mna[324](19) = mna[324][19] -Geq_M;
	rval[19] += -(Ieq_M - Qo_M)*idt;
	fvec[19] += Sr_M;
	mna[19](38) = mna[19][38] - idt*Geq_M;
	mna[38](19) = mna[38][19] - idt*Geq_M;
	mna[38](38) = mna[38][38] + idt*Geq_M;
	mna[324](38) = mna[324][38] +Geq_M;
	rval[38] += (Ieq_M - Qo_M)*idt;
	fvec[38] += -(Sr_M);
	mna[324](324) = mna[324][324] + 1.;
	rval[324] += Ieq_M;
	fvec[324] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM284*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[24];
	Un_M += Vx[20];
	Uo_M += -uxtold[24];
	Uo_M += uxtold[20];
	Qo_M = uxtold[325];
	Qn_M = Vx[325];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[20](20) = mna[20][20] + idt*Geq_M;
	mna[325](20) = mna[325][20] -Geq_M;
	rval[20] += -(Ieq_M - Qo_M)*idt;
	fvec[20] += Sr_M;
	mna[20](24) = mna[20][24] - idt*Geq_M;
	mna[24](20) = mna[24][20] - idt*Geq_M;
	mna[24](24) = mna[24][24] + idt*Geq_M;
	mna[325](24) = mna[325][24] +Geq_M;
	rval[24] += (Ieq_M - Qo_M)*idt;
	fvec[24] += -(Sr_M);
	mna[325](325) = mna[325][325] + 1.;
	rval[325] += Ieq_M;
	fvec[325] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM285*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[26];
	Un_M += Vx[20];
	Uo_M += -uxtold[26];
	Uo_M += uxtold[20];
	Qo_M = uxtold[326];
	Qn_M = Vx[326];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[20](20) = mna[20][20] + idt*Geq_M;
	mna[326](20) = mna[326][20] -Geq_M;
	rval[20] += -(Ieq_M - Qo_M)*idt;
	fvec[20] += Sr_M;
	mna[20](26) = mna[20][26] - idt*Geq_M;
	mna[26](20) = mna[26][20] - idt*Geq_M;
	mna[26](26) = mna[26][26] + idt*Geq_M;
	mna[326](26) = mna[326][26] +Geq_M;
	rval[26] += (Ieq_M - Qo_M)*idt;
	fvec[26] += -(Sr_M);
	mna[326](326) = mna[326][326] + 1.;
	rval[326] += Ieq_M;
	fvec[326] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM286*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[27];
	Un_M += Vx[20];
	Uo_M += -uxtold[27];
	Uo_M += uxtold[20];
	Qo_M = uxtold[327];
	Qn_M = Vx[327];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[20](20) = mna[20][20] + idt*Geq_M;
	mna[327](20) = mna[327][20] -Geq_M;
	rval[20] += -(Ieq_M - Qo_M)*idt;
	fvec[20] += Sr_M;
	mna[20](27) = mna[20][27] - idt*Geq_M;
	mna[27](20) = mna[27][20] - idt*Geq_M;
	mna[27](27) = mna[27][27] + idt*Geq_M;
	mna[327](27) = mna[327][27] +Geq_M;
	rval[27] += (Ieq_M - Qo_M)*idt;
	fvec[27] += -(Sr_M);
	mna[327](327) = mna[327][327] + 1.;
	rval[327] += Ieq_M;
	fvec[327] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM287*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[33];
	Un_M += Vx[20];
	Uo_M += -uxtold[33];
	Uo_M += uxtold[20];
	Qo_M = uxtold[328];
	Qn_M = Vx[328];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[20](20) = mna[20][20] + idt*Geq_M;
	mna[328](20) = mna[328][20] -Geq_M;
	rval[20] += -(Ieq_M - Qo_M)*idt;
	fvec[20] += Sr_M;
	mna[20](33) = mna[20][33] - idt*Geq_M;
	mna[33](20) = mna[33][20] - idt*Geq_M;
	mna[33](33) = mna[33][33] + idt*Geq_M;
	mna[328](33) = mna[328][33] +Geq_M;
	rval[33] += (Ieq_M - Qo_M)*idt;
	fvec[33] += -(Sr_M);
	mna[328](328) = mna[328][328] + 1.;
	rval[328] += Ieq_M;
	fvec[328] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM288*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[37];
	Un_M += Vx[20];
	Uo_M += -uxtold[37];
	Uo_M += uxtold[20];
	Qo_M = uxtold[329];
	Qn_M = Vx[329];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[20](20) = mna[20][20] + idt*Geq_M;
	mna[329](20) = mna[329][20] -Geq_M;
	rval[20] += -(Ieq_M - Qo_M)*idt;
	fvec[20] += Sr_M;
	mna[20](37) = mna[20][37] - idt*Geq_M;
	mna[37](20) = mna[37][20] - idt*Geq_M;
	mna[37](37) = mna[37][37] + idt*Geq_M;
	mna[329](37) = mna[329][37] +Geq_M;
	rval[37] += (Ieq_M - Qo_M)*idt;
	fvec[37] += -(Sr_M);
	mna[329](329) = mna[329][329] + 1.;
	rval[329] += Ieq_M;
	fvec[329] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM289*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[39];
	Un_M += Vx[20];
	Uo_M += -uxtold[39];
	Uo_M += uxtold[20];
	Qo_M = uxtold[330];
	Qn_M = Vx[330];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[20](20) = mna[20][20] + idt*Geq_M;
	mna[330](20) = mna[330][20] -Geq_M;
	rval[20] += -(Ieq_M - Qo_M)*idt;
	fvec[20] += Sr_M;
	mna[20](39) = mna[20][39] - idt*Geq_M;
	mna[39](20) = mna[39][20] - idt*Geq_M;
	mna[39](39) = mna[39][39] + idt*Geq_M;
	mna[330](39) = mna[330][39] +Geq_M;
	rval[39] += (Ieq_M - Qo_M)*idt;
	fvec[39] += -(Sr_M);
	mna[330](330) = mna[330][330] + 1.;
	rval[330] += Ieq_M;
	fvec[330] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM290*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[22];
	Un_M += Vx[21];
	Uo_M += -uxtold[22];
	Uo_M += uxtold[21];
	Qo_M = uxtold[331];
	Qn_M = Vx[331];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[21](21) = mna[21][21] + idt*Geq_M;
	mna[331](21) = mna[331][21] -Geq_M;
	rval[21] += -(Ieq_M - Qo_M)*idt;
	fvec[21] += Sr_M;
	mna[21](22) = mna[21][22] - idt*Geq_M;
	mna[22](21) = mna[22][21] - idt*Geq_M;
	mna[22](22) = mna[22][22] + idt*Geq_M;
	mna[331](22) = mna[331][22] +Geq_M;
	rval[22] += (Ieq_M - Qo_M)*idt;
	fvec[22] += -(Sr_M);
	mna[331](331) = mna[331][331] + 1.;
	rval[331] += Ieq_M;
	fvec[331] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM291*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[23];
	Un_M += Vx[21];
	Uo_M += -uxtold[23];
	Uo_M += uxtold[21];
	Qo_M = uxtold[332];
	Qn_M = Vx[332];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[21](21) = mna[21][21] + idt*Geq_M;
	mna[332](21) = mna[332][21] -Geq_M;
	rval[21] += -(Ieq_M - Qo_M)*idt;
	fvec[21] += Sr_M;
	mna[21](23) = mna[21][23] - idt*Geq_M;
	mna[23](21) = mna[23][21] - idt*Geq_M;
	mna[23](23) = mna[23][23] + idt*Geq_M;
	mna[332](23) = mna[332][23] +Geq_M;
	rval[23] += (Ieq_M - Qo_M)*idt;
	fvec[23] += -(Sr_M);
	mna[332](332) = mna[332][332] + 1.;
	rval[332] += Ieq_M;
	fvec[332] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM292*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[26];
	Un_M += Vx[21];
	Uo_M += -uxtold[26];
	Uo_M += uxtold[21];
	Qo_M = uxtold[333];
	Qn_M = Vx[333];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[21](21) = mna[21][21] + idt*Geq_M;
	mna[333](21) = mna[333][21] -Geq_M;
	rval[21] += -(Ieq_M - Qo_M)*idt;
	fvec[21] += Sr_M;
	mna[21](26) = mna[21][26] - idt*Geq_M;
	mna[26](21) = mna[26][21] - idt*Geq_M;
	mna[26](26) = mna[26][26] + idt*Geq_M;
	mna[333](26) = mna[333][26] +Geq_M;
	rval[26] += (Ieq_M - Qo_M)*idt;
	fvec[26] += -(Sr_M);
	mna[333](333) = mna[333][333] + 1.;
	rval[333] += Ieq_M;
	fvec[333] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM293*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[27];
	Un_M += Vx[21];
	Uo_M += -uxtold[27];
	Uo_M += uxtold[21];
	Qo_M = uxtold[334];
	Qn_M = Vx[334];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[21](21) = mna[21][21] + idt*Geq_M;
	mna[334](21) = mna[334][21] -Geq_M;
	rval[21] += -(Ieq_M - Qo_M)*idt;
	fvec[21] += Sr_M;
	mna[21](27) = mna[21][27] - idt*Geq_M;
	mna[27](21) = mna[27][21] - idt*Geq_M;
	mna[27](27) = mna[27][27] + idt*Geq_M;
	mna[334](27) = mna[334][27] +Geq_M;
	rval[27] += (Ieq_M - Qo_M)*idt;
	fvec[27] += -(Sr_M);
	mna[334](334) = mna[334][334] + 1.;
	rval[334] += Ieq_M;
	fvec[334] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM294*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[37];
	Un_M += Vx[21];
	Uo_M += -uxtold[37];
	Uo_M += uxtold[21];
	Qo_M = uxtold[335];
	Qn_M = Vx[335];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[21](21) = mna[21][21] + idt*Geq_M;
	mna[335](21) = mna[335][21] -Geq_M;
	rval[21] += -(Ieq_M - Qo_M)*idt;
	fvec[21] += Sr_M;
	mna[21](37) = mna[21][37] - idt*Geq_M;
	mna[37](21) = mna[37][21] - idt*Geq_M;
	mna[37](37) = mna[37][37] + idt*Geq_M;
	mna[335](37) = mna[335][37] +Geq_M;
	rval[37] += (Ieq_M - Qo_M)*idt;
	fvec[37] += -(Sr_M);
	mna[335](335) = mna[335][335] + 1.;
	rval[335] += Ieq_M;
	fvec[335] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM295*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[38];
	Un_M += Vx[21];
	Uo_M += -uxtold[38];
	Uo_M += uxtold[21];
	Qo_M = uxtold[336];
	Qn_M = Vx[336];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[21](21) = mna[21][21] + idt*Geq_M;
	mna[336](21) = mna[336][21] -Geq_M;
	rval[21] += -(Ieq_M - Qo_M)*idt;
	fvec[21] += Sr_M;
	mna[21](38) = mna[21][38] - idt*Geq_M;
	mna[38](21) = mna[38][21] - idt*Geq_M;
	mna[38](38) = mna[38][38] + idt*Geq_M;
	mna[336](38) = mna[336][38] +Geq_M;
	rval[38] += (Ieq_M - Qo_M)*idt;
	fvec[38] += -(Sr_M);
	mna[336](336) = mna[336][336] + 1.;
	rval[336] += Ieq_M;
	fvec[336] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM296*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[39];
	Un_M += Vx[21];
	Uo_M += -uxtold[39];
	Uo_M += uxtold[21];
	Qo_M = uxtold[337];
	Qn_M = Vx[337];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[21](21) = mna[21][21] + idt*Geq_M;
	mna[337](21) = mna[337][21] -Geq_M;
	rval[21] += -(Ieq_M - Qo_M)*idt;
	fvec[21] += Sr_M;
	mna[21](39) = mna[21][39] - idt*Geq_M;
	mna[39](21) = mna[39][21] - idt*Geq_M;
	mna[39](39) = mna[39][39] + idt*Geq_M;
	mna[337](39) = mna[337][39] +Geq_M;
	rval[39] += (Ieq_M - Qo_M)*idt;
	fvec[39] += -(Sr_M);
	mna[337](337) = mna[337][337] + 1.;
	rval[337] += Ieq_M;
	fvec[337] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM297*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[40];
	Un_M += Vx[21];
	Uo_M += -uxtold[40];
	Uo_M += uxtold[21];
	Qo_M = uxtold[338];
	Qn_M = Vx[338];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[21](21) = mna[21][21] + idt*Geq_M;
	mna[338](21) = mna[338][21] -Geq_M;
	rval[21] += -(Ieq_M - Qo_M)*idt;
	fvec[21] += Sr_M;
	mna[21](40) = mna[21][40] - idt*Geq_M;
	mna[40](21) = mna[40][21] - idt*Geq_M;
	mna[40](40) = mna[40][40] + idt*Geq_M;
	mna[338](40) = mna[338][40] +Geq_M;
	rval[40] += (Ieq_M - Qo_M)*idt;
	fvec[40] += -(Sr_M);
	mna[338](338) = mna[338][338] + 1.;
	rval[338] += Ieq_M;
	fvec[338] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM298*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[29];
	Un_M += Vx[22];
	Uo_M += -uxtold[29];
	Uo_M += uxtold[22];
	Qo_M = uxtold[339];
	Qn_M = Vx[339];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[22](22) = mna[22][22] + idt*Geq_M;
	mna[339](22) = mna[339][22] -Geq_M;
	rval[22] += -(Ieq_M - Qo_M)*idt;
	fvec[22] += Sr_M;
	mna[22](29) = mna[22][29] - idt*Geq_M;
	mna[29](22) = mna[29][22] - idt*Geq_M;
	mna[29](29) = mna[29][29] + idt*Geq_M;
	mna[339](29) = mna[339][29] +Geq_M;
	rval[29] += (Ieq_M - Qo_M)*idt;
	fvec[29] += -(Sr_M);
	mna[339](339) = mna[339][339] + 1.;
	rval[339] += Ieq_M;
	fvec[339] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM299*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[30];
	Un_M += Vx[22];
	Uo_M += -uxtold[30];
	Uo_M += uxtold[22];
	Qo_M = uxtold[340];
	Qn_M = Vx[340];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[22](22) = mna[22][22] + idt*Geq_M;
	mna[340](22) = mna[340][22] -Geq_M;
	rval[22] += -(Ieq_M - Qo_M)*idt;
	fvec[22] += Sr_M;
	mna[22](30) = mna[22][30] - idt*Geq_M;
	mna[30](22) = mna[30][22] - idt*Geq_M;
	mna[30](30) = mna[30][30] + idt*Geq_M;
	mna[340](30) = mna[340][30] +Geq_M;
	rval[30] += (Ieq_M - Qo_M)*idt;
	fvec[30] += -(Sr_M);
	mna[340](340) = mna[340][340] + 1.;
	rval[340] += Ieq_M;
	fvec[340] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM300*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[33];
	Un_M += Vx[22];
	Uo_M += -uxtold[33];
	Uo_M += uxtold[22];
	Qo_M = uxtold[341];
	Qn_M = Vx[341];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[22](22) = mna[22][22] + idt*Geq_M;
	mna[341](22) = mna[341][22] -Geq_M;
	rval[22] += -(Ieq_M - Qo_M)*idt;
	fvec[22] += Sr_M;
	mna[22](33) = mna[22][33] - idt*Geq_M;
	mna[33](22) = mna[33][22] - idt*Geq_M;
	mna[33](33) = mna[33][33] + idt*Geq_M;
	mna[341](33) = mna[341][33] +Geq_M;
	rval[33] += (Ieq_M - Qo_M)*idt;
	fvec[33] += -(Sr_M);
	mna[341](341) = mna[341][341] + 1.;
	rval[341] += Ieq_M;
	fvec[341] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM301*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[39];
	Un_M += Vx[22];
	Uo_M += -uxtold[39];
	Uo_M += uxtold[22];
	Qo_M = uxtold[342];
	Qn_M = Vx[342];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[22](22) = mna[22][22] + idt*Geq_M;
	mna[342](22) = mna[342][22] -Geq_M;
	rval[22] += -(Ieq_M - Qo_M)*idt;
	fvec[22] += Sr_M;
	mna[22](39) = mna[22][39] - idt*Geq_M;
	mna[39](22) = mna[39][22] - idt*Geq_M;
	mna[39](39) = mna[39][39] + idt*Geq_M;
	mna[342](39) = mna[342][39] +Geq_M;
	rval[39] += (Ieq_M - Qo_M)*idt;
	fvec[39] += -(Sr_M);
	mna[342](342) = mna[342][342] + 1.;
	rval[342] += Ieq_M;
	fvec[342] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM302*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[27];
	Un_M += Vx[23];
	Uo_M += -uxtold[27];
	Uo_M += uxtold[23];
	Qo_M = uxtold[343];
	Qn_M = Vx[343];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[23](23) = mna[23][23] + idt*Geq_M;
	mna[343](23) = mna[343][23] -Geq_M;
	rval[23] += -(Ieq_M - Qo_M)*idt;
	fvec[23] += Sr_M;
	mna[23](27) = mna[23][27] - idt*Geq_M;
	mna[27](23) = mna[27][23] - idt*Geq_M;
	mna[27](27) = mna[27][27] + idt*Geq_M;
	mna[343](27) = mna[343][27] +Geq_M;
	rval[27] += (Ieq_M - Qo_M)*idt;
	fvec[27] += -(Sr_M);
	mna[343](343) = mna[343][343] + 1.;
	rval[343] += Ieq_M;
	fvec[343] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM303*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[28];
	Un_M += Vx[23];
	Uo_M += -uxtold[28];
	Uo_M += uxtold[23];
	Qo_M = uxtold[344];
	Qn_M = Vx[344];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[23](23) = mna[23][23] + idt*Geq_M;
	mna[344](23) = mna[344][23] -Geq_M;
	rval[23] += -(Ieq_M - Qo_M)*idt;
	fvec[23] += Sr_M;
	mna[23](28) = mna[23][28] - idt*Geq_M;
	mna[28](23) = mna[28][23] - idt*Geq_M;
	mna[28](28) = mna[28][28] + idt*Geq_M;
	mna[344](28) = mna[344][28] +Geq_M;
	rval[28] += (Ieq_M - Qo_M)*idt;
	fvec[28] += -(Sr_M);
	mna[344](344) = mna[344][344] + 1.;
	rval[344] += Ieq_M;
	fvec[344] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM304*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[30];
	Un_M += Vx[23];
	Uo_M += -uxtold[30];
	Uo_M += uxtold[23];
	Qo_M = uxtold[345];
	Qn_M = Vx[345];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[23](23) = mna[23][23] + idt*Geq_M;
	mna[345](23) = mna[345][23] -Geq_M;
	rval[23] += -(Ieq_M - Qo_M)*idt;
	fvec[23] += Sr_M;
	mna[23](30) = mna[23][30] - idt*Geq_M;
	mna[30](23) = mna[30][23] - idt*Geq_M;
	mna[30](30) = mna[30][30] + idt*Geq_M;
	mna[345](30) = mna[345][30] +Geq_M;
	rval[30] += (Ieq_M - Qo_M)*idt;
	fvec[30] += -(Sr_M);
	mna[345](345) = mna[345][345] + 1.;
	rval[345] += Ieq_M;
	fvec[345] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM305*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[31];
	Un_M += Vx[23];
	Uo_M += -uxtold[31];
	Uo_M += uxtold[23];
	Qo_M = uxtold[346];
	Qn_M = Vx[346];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[23](23) = mna[23][23] + idt*Geq_M;
	mna[346](23) = mna[346][23] -Geq_M;
	rval[23] += -(Ieq_M - Qo_M)*idt;
	fvec[23] += Sr_M;
	mna[23](31) = mna[23][31] - idt*Geq_M;
	mna[31](23) = mna[31][23] - idt*Geq_M;
	mna[31](31) = mna[31][31] + idt*Geq_M;
	mna[346](31) = mna[346][31] +Geq_M;
	rval[31] += (Ieq_M - Qo_M)*idt;
	fvec[31] += -(Sr_M);
	mna[346](346) = mna[346][346] + 1.;
	rval[346] += Ieq_M;
	fvec[346] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM306*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[37];
	Un_M += Vx[23];
	Uo_M += -uxtold[37];
	Uo_M += uxtold[23];
	Qo_M = uxtold[347];
	Qn_M = Vx[347];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[23](23) = mna[23][23] + idt*Geq_M;
	mna[347](23) = mna[347][23] -Geq_M;
	rval[23] += -(Ieq_M - Qo_M)*idt;
	fvec[23] += Sr_M;
	mna[23](37) = mna[23][37] - idt*Geq_M;
	mna[37](23) = mna[37][23] - idt*Geq_M;
	mna[37](37) = mna[37][37] + idt*Geq_M;
	mna[347](37) = mna[347][37] +Geq_M;
	rval[37] += (Ieq_M - Qo_M)*idt;
	fvec[37] += -(Sr_M);
	mna[347](347) = mna[347][347] + 1.;
	rval[347] += Ieq_M;
	fvec[347] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM307*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[38];
	Un_M += Vx[23];
	Uo_M += -uxtold[38];
	Uo_M += uxtold[23];
	Qo_M = uxtold[348];
	Qn_M = Vx[348];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[23](23) = mna[23][23] + idt*Geq_M;
	mna[348](23) = mna[348][23] -Geq_M;
	rval[23] += -(Ieq_M - Qo_M)*idt;
	fvec[23] += Sr_M;
	mna[23](38) = mna[23][38] - idt*Geq_M;
	mna[38](23) = mna[38][23] - idt*Geq_M;
	mna[38](38) = mna[38][38] + idt*Geq_M;
	mna[348](38) = mna[348][38] +Geq_M;
	rval[38] += (Ieq_M - Qo_M)*idt;
	fvec[38] += -(Sr_M);
	mna[348](348) = mna[348][348] + 1.;
	rval[348] += Ieq_M;
	fvec[348] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM308*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[39];
	Un_M += Vx[23];
	Uo_M += -uxtold[39];
	Uo_M += uxtold[23];
	Qo_M = uxtold[349];
	Qn_M = Vx[349];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[23](23) = mna[23][23] + idt*Geq_M;
	mna[349](23) = mna[349][23] -Geq_M;
	rval[23] += -(Ieq_M - Qo_M)*idt;
	fvec[23] += Sr_M;
	mna[23](39) = mna[23][39] - idt*Geq_M;
	mna[39](23) = mna[39][23] - idt*Geq_M;
	mna[39](39) = mna[39][39] + idt*Geq_M;
	mna[349](39) = mna[349][39] +Geq_M;
	rval[39] += (Ieq_M - Qo_M)*idt;
	fvec[39] += -(Sr_M);
	mna[349](349) = mna[349][349] + 1.;
	rval[349] += Ieq_M;
	fvec[349] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM309*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[31];
	Un_M += Vx[24];
	Uo_M += -uxtold[31];
	Uo_M += uxtold[24];
	Qo_M = uxtold[350];
	Qn_M = Vx[350];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[24](24) = mna[24][24] + idt*Geq_M;
	mna[350](24) = mna[350][24] -Geq_M;
	rval[24] += -(Ieq_M - Qo_M)*idt;
	fvec[24] += Sr_M;
	mna[24](31) = mna[24][31] - idt*Geq_M;
	mna[31](24) = mna[31][24] - idt*Geq_M;
	mna[31](31) = mna[31][31] + idt*Geq_M;
	mna[350](31) = mna[350][31] +Geq_M;
	rval[31] += (Ieq_M - Qo_M)*idt;
	fvec[31] += -(Sr_M);
	mna[350](350) = mna[350][350] + 1.;
	rval[350] += Ieq_M;
	fvec[350] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM310*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[37];
	Un_M += Vx[24];
	Uo_M += -uxtold[37];
	Uo_M += uxtold[24];
	Qo_M = uxtold[351];
	Qn_M = Vx[351];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[24](24) = mna[24][24] + idt*Geq_M;
	mna[351](24) = mna[351][24] -Geq_M;
	rval[24] += -(Ieq_M - Qo_M)*idt;
	fvec[24] += Sr_M;
	mna[24](37) = mna[24][37] - idt*Geq_M;
	mna[37](24) = mna[37][24] - idt*Geq_M;
	mna[37](37) = mna[37][37] + idt*Geq_M;
	mna[351](37) = mna[351][37] +Geq_M;
	rval[37] += (Ieq_M - Qo_M)*idt;
	fvec[37] += -(Sr_M);
	mna[351](351) = mna[351][351] + 1.;
	rval[351] += Ieq_M;
	fvec[351] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM311*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[28];
	Un_M += Vx[25];
	Uo_M += -uxtold[28];
	Uo_M += uxtold[25];
	Qo_M = uxtold[352];
	Qn_M = Vx[352];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[25](25) = mna[25][25] + idt*Geq_M;
	mna[352](25) = mna[352][25] -Geq_M;
	rval[25] += -(Ieq_M - Qo_M)*idt;
	fvec[25] += Sr_M;
	mna[25](28) = mna[25][28] - idt*Geq_M;
	mna[28](25) = mna[28][25] - idt*Geq_M;
	mna[28](28) = mna[28][28] + idt*Geq_M;
	mna[352](28) = mna[352][28] +Geq_M;
	rval[28] += (Ieq_M - Qo_M)*idt;
	fvec[28] += -(Sr_M);
	mna[352](352) = mna[352][352] + 1.;
	rval[352] += Ieq_M;
	fvec[352] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM312*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[30];
	Un_M += Vx[25];
	Uo_M += -uxtold[30];
	Uo_M += uxtold[25];
	Qo_M = uxtold[353];
	Qn_M = Vx[353];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[25](25) = mna[25][25] + idt*Geq_M;
	mna[353](25) = mna[353][25] -Geq_M;
	rval[25] += -(Ieq_M - Qo_M)*idt;
	fvec[25] += Sr_M;
	mna[25](30) = mna[25][30] - idt*Geq_M;
	mna[30](25) = mna[30][25] - idt*Geq_M;
	mna[30](30) = mna[30][30] + idt*Geq_M;
	mna[353](30) = mna[353][30] +Geq_M;
	rval[30] += (Ieq_M - Qo_M)*idt;
	fvec[30] += -(Sr_M);
	mna[353](353) = mna[353][353] + 1.;
	rval[353] += Ieq_M;
	fvec[353] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM313*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[28];
	Un_M += Vx[26];
	Uo_M += -uxtold[28];
	Uo_M += uxtold[26];
	Qo_M = uxtold[354];
	Qn_M = Vx[354];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[26](26) = mna[26][26] + idt*Geq_M;
	mna[354](26) = mna[354][26] -Geq_M;
	rval[26] += -(Ieq_M - Qo_M)*idt;
	fvec[26] += Sr_M;
	mna[26](28) = mna[26][28] - idt*Geq_M;
	mna[28](26) = mna[28][26] - idt*Geq_M;
	mna[28](28) = mna[28][28] + idt*Geq_M;
	mna[354](28) = mna[354][28] +Geq_M;
	rval[28] += (Ieq_M - Qo_M)*idt;
	fvec[28] += -(Sr_M);
	mna[354](354) = mna[354][354] + 1.;
	rval[354] += Ieq_M;
	fvec[354] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM314*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[29];
	Un_M += Vx[26];
	Uo_M += -uxtold[29];
	Uo_M += uxtold[26];
	Qo_M = uxtold[355];
	Qn_M = Vx[355];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[26](26) = mna[26][26] + idt*Geq_M;
	mna[355](26) = mna[355][26] -Geq_M;
	rval[26] += -(Ieq_M - Qo_M)*idt;
	fvec[26] += Sr_M;
	mna[26](29) = mna[26][29] - idt*Geq_M;
	mna[29](26) = mna[29][26] - idt*Geq_M;
	mna[29](29) = mna[29][29] + idt*Geq_M;
	mna[355](29) = mna[355][29] +Geq_M;
	rval[29] += (Ieq_M - Qo_M)*idt;
	fvec[29] += -(Sr_M);
	mna[355](355) = mna[355][355] + 1.;
	rval[355] += Ieq_M;
	fvec[355] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM315*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[32];
	Un_M += Vx[26];
	Uo_M += -uxtold[32];
	Uo_M += uxtold[26];
	Qo_M = uxtold[356];
	Qn_M = Vx[356];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[26](26) = mna[26][26] + idt*Geq_M;
	mna[356](26) = mna[356][26] -Geq_M;
	rval[26] += -(Ieq_M - Qo_M)*idt;
	fvec[26] += Sr_M;
	mna[26](32) = mna[26][32] - idt*Geq_M;
	mna[32](26) = mna[32][26] - idt*Geq_M;
	mna[32](32) = mna[32][32] + idt*Geq_M;
	mna[356](32) = mna[356][32] +Geq_M;
	rval[32] += (Ieq_M - Qo_M)*idt;
	fvec[32] += -(Sr_M);
	mna[356](356) = mna[356][356] + 1.;
	rval[356] += Ieq_M;
	fvec[356] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM316*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[33];
	Un_M += Vx[26];
	Uo_M += -uxtold[33];
	Uo_M += uxtold[26];
	Qo_M = uxtold[357];
	Qn_M = Vx[357];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[26](26) = mna[26][26] + idt*Geq_M;
	mna[357](26) = mna[357][26] -Geq_M;
	rval[26] += -(Ieq_M - Qo_M)*idt;
	fvec[26] += Sr_M;
	mna[26](33) = mna[26][33] - idt*Geq_M;
	mna[33](26) = mna[33][26] - idt*Geq_M;
	mna[33](33) = mna[33][33] + idt*Geq_M;
	mna[357](33) = mna[357][33] +Geq_M;
	rval[33] += (Ieq_M - Qo_M)*idt;
	fvec[33] += -(Sr_M);
	mna[357](357) = mna[357][357] + 1.;
	rval[357] += Ieq_M;
	fvec[357] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM317*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[29];
	Un_M += Vx[27];
	Uo_M += -uxtold[29];
	Uo_M += uxtold[27];
	Qo_M = uxtold[358];
	Qn_M = Vx[358];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[27](27) = mna[27][27] + idt*Geq_M;
	mna[358](27) = mna[358][27] -Geq_M;
	rval[27] += -(Ieq_M - Qo_M)*idt;
	fvec[27] += Sr_M;
	mna[27](29) = mna[27][29] - idt*Geq_M;
	mna[29](27) = mna[29][27] - idt*Geq_M;
	mna[29](29) = mna[29][29] + idt*Geq_M;
	mna[358](29) = mna[358][29] +Geq_M;
	rval[29] += (Ieq_M - Qo_M)*idt;
	fvec[29] += -(Sr_M);
	mna[358](358) = mna[358][358] + 1.;
	rval[358] += Ieq_M;
	fvec[358] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM318*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[33];
	Un_M += Vx[27];
	Uo_M += -uxtold[33];
	Uo_M += uxtold[27];
	Qo_M = uxtold[359];
	Qn_M = Vx[359];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[27](27) = mna[27][27] + idt*Geq_M;
	mna[359](27) = mna[359][27] -Geq_M;
	rval[27] += -(Ieq_M - Qo_M)*idt;
	fvec[27] += Sr_M;
	mna[27](33) = mna[27][33] - idt*Geq_M;
	mna[33](27) = mna[33][27] - idt*Geq_M;
	mna[33](33) = mna[33][33] + idt*Geq_M;
	mna[359](33) = mna[359][33] +Geq_M;
	rval[33] += (Ieq_M - Qo_M)*idt;
	fvec[33] += -(Sr_M);
	mna[359](359) = mna[359][359] + 1.;
	rval[359] += Ieq_M;
	fvec[359] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM319*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[37];
	Un_M += Vx[27];
	Uo_M += -uxtold[37];
	Uo_M += uxtold[27];
	Qo_M = uxtold[360];
	Qn_M = Vx[360];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[27](27) = mna[27][27] + idt*Geq_M;
	mna[360](27) = mna[360][27] -Geq_M;
	rval[27] += -(Ieq_M - Qo_M)*idt;
	fvec[27] += Sr_M;
	mna[27](37) = mna[27][37] - idt*Geq_M;
	mna[37](27) = mna[37][27] - idt*Geq_M;
	mna[37](37) = mna[37][37] + idt*Geq_M;
	mna[360](37) = mna[360][37] +Geq_M;
	rval[37] += (Ieq_M - Qo_M)*idt;
	fvec[37] += -(Sr_M);
	mna[360](360) = mna[360][360] + 1.;
	rval[360] += Ieq_M;
	fvec[360] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM320*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[31];
	Un_M += Vx[28];
	Uo_M += -uxtold[31];
	Uo_M += uxtold[28];
	Qo_M = uxtold[361];
	Qn_M = Vx[361];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[28](28) = mna[28][28] + idt*Geq_M;
	mna[361](28) = mna[361][28] -Geq_M;
	rval[28] += -(Ieq_M - Qo_M)*idt;
	fvec[28] += Sr_M;
	mna[28](31) = mna[28][31] - idt*Geq_M;
	mna[31](28) = mna[31][28] - idt*Geq_M;
	mna[31](31) = mna[31][31] + idt*Geq_M;
	mna[361](31) = mna[361][31] +Geq_M;
	rval[31] += (Ieq_M - Qo_M)*idt;
	fvec[31] += -(Sr_M);
	mna[361](361) = mna[361][361] + 1.;
	rval[361] += Ieq_M;
	fvec[361] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM321*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[32];
	Un_M += Vx[28];
	Uo_M += -uxtold[32];
	Uo_M += uxtold[28];
	Qo_M = uxtold[362];
	Qn_M = Vx[362];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[28](28) = mna[28][28] + idt*Geq_M;
	mna[362](28) = mna[362][28] -Geq_M;
	rval[28] += -(Ieq_M - Qo_M)*idt;
	fvec[28] += Sr_M;
	mna[28](32) = mna[28][32] - idt*Geq_M;
	mna[32](28) = mna[32][28] - idt*Geq_M;
	mna[32](32) = mna[32][32] + idt*Geq_M;
	mna[362](32) = mna[362][32] +Geq_M;
	rval[32] += (Ieq_M - Qo_M)*idt;
	fvec[32] += -(Sr_M);
	mna[362](362) = mna[362][362] + 1.;
	rval[362] += Ieq_M;
	fvec[362] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM322*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[33];
	Un_M += Vx[28];
	Uo_M += -uxtold[33];
	Uo_M += uxtold[28];
	Qo_M = uxtold[363];
	Qn_M = Vx[363];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[28](28) = mna[28][28] + idt*Geq_M;
	mna[363](28) = mna[363][28] -Geq_M;
	rval[28] += -(Ieq_M - Qo_M)*idt;
	fvec[28] += Sr_M;
	mna[28](33) = mna[28][33] - idt*Geq_M;
	mna[33](28) = mna[33][28] - idt*Geq_M;
	mna[33](33) = mna[33][33] + idt*Geq_M;
	mna[363](33) = mna[363][33] +Geq_M;
	rval[33] += (Ieq_M - Qo_M)*idt;
	fvec[33] += -(Sr_M);
	mna[363](363) = mna[363][363] + 1.;
	rval[363] += Ieq_M;
	fvec[363] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM323*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[36];
	Un_M += Vx[28];
	Uo_M += -uxtold[36];
	Uo_M += uxtold[28];
	Qo_M = uxtold[364];
	Qn_M = Vx[364];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[28](28) = mna[28][28] + idt*Geq_M;
	mna[364](28) = mna[364][28] -Geq_M;
	rval[28] += -(Ieq_M - Qo_M)*idt;
	fvec[28] += Sr_M;
	mna[28](36) = mna[28][36] - idt*Geq_M;
	mna[36](28) = mna[36][28] - idt*Geq_M;
	mna[36](36) = mna[36][36] + idt*Geq_M;
	mna[364](36) = mna[364][36] +Geq_M;
	rval[36] += (Ieq_M - Qo_M)*idt;
	fvec[36] += -(Sr_M);
	mna[364](364) = mna[364][364] + 1.;
	rval[364] += Ieq_M;
	fvec[364] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM324*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[37];
	Un_M += Vx[28];
	Uo_M += -uxtold[37];
	Uo_M += uxtold[28];
	Qo_M = uxtold[365];
	Qn_M = Vx[365];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[28](28) = mna[28][28] + idt*Geq_M;
	mna[365](28) = mna[365][28] -Geq_M;
	rval[28] += -(Ieq_M - Qo_M)*idt;
	fvec[28] += Sr_M;
	mna[28](37) = mna[28][37] - idt*Geq_M;
	mna[37](28) = mna[37][28] - idt*Geq_M;
	mna[37](37) = mna[37][37] + idt*Geq_M;
	mna[365](37) = mna[365][37] +Geq_M;
	rval[37] += (Ieq_M - Qo_M)*idt;
	fvec[37] += -(Sr_M);
	mna[365](365) = mna[365][365] + 1.;
	rval[365] += Ieq_M;
	fvec[365] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM325*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[32];
	Un_M += Vx[29];
	Uo_M += -uxtold[32];
	Uo_M += uxtold[29];
	Qo_M = uxtold[366];
	Qn_M = Vx[366];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[29](29) = mna[29][29] + idt*Geq_M;
	mna[366](29) = mna[366][29] -Geq_M;
	rval[29] += -(Ieq_M - Qo_M)*idt;
	fvec[29] += Sr_M;
	mna[29](32) = mna[29][32] - idt*Geq_M;
	mna[32](29) = mna[32][29] - idt*Geq_M;
	mna[32](32) = mna[32][32] + idt*Geq_M;
	mna[366](32) = mna[366][32] +Geq_M;
	rval[32] += (Ieq_M - Qo_M)*idt;
	fvec[32] += -(Sr_M);
	mna[366](366) = mna[366][366] + 1.;
	rval[366] += Ieq_M;
	fvec[366] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM326*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[35];
	Un_M += Vx[29];
	Uo_M += -uxtold[35];
	Uo_M += uxtold[29];
	Qo_M = uxtold[367];
	Qn_M = Vx[367];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[29](29) = mna[29][29] + idt*Geq_M;
	mna[367](29) = mna[367][29] -Geq_M;
	rval[29] += -(Ieq_M - Qo_M)*idt;
	fvec[29] += Sr_M;
	mna[29](35) = mna[29][35] - idt*Geq_M;
	mna[35](29) = mna[35][29] - idt*Geq_M;
	mna[35](35) = mna[35][35] + idt*Geq_M;
	mna[367](35) = mna[367][35] +Geq_M;
	rval[35] += (Ieq_M - Qo_M)*idt;
	fvec[35] += -(Sr_M);
	mna[367](367) = mna[367][367] + 1.;
	rval[367] += Ieq_M;
	fvec[367] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM327*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[37];
	Un_M += Vx[29];
	Uo_M += -uxtold[37];
	Uo_M += uxtold[29];
	Qo_M = uxtold[368];
	Qn_M = Vx[368];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[29](29) = mna[29][29] + idt*Geq_M;
	mna[368](29) = mna[368][29] -Geq_M;
	rval[29] += -(Ieq_M - Qo_M)*idt;
	fvec[29] += Sr_M;
	mna[29](37) = mna[29][37] - idt*Geq_M;
	mna[37](29) = mna[37][29] - idt*Geq_M;
	mna[37](37) = mna[37][37] + idt*Geq_M;
	mna[368](37) = mna[368][37] +Geq_M;
	rval[37] += (Ieq_M - Qo_M)*idt;
	fvec[37] += -(Sr_M);
	mna[368](368) = mna[368][368] + 1.;
	rval[368] += Ieq_M;
	fvec[368] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM328*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[40];
	Un_M += Vx[30];
	Uo_M += -uxtold[40];
	Uo_M += uxtold[30];
	Qo_M = uxtold[369];
	Qn_M = Vx[369];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[30](30) = mna[30][30] + idt*Geq_M;
	mna[369](30) = mna[369][30] -Geq_M;
	rval[30] += -(Ieq_M - Qo_M)*idt;
	fvec[30] += Sr_M;
	mna[30](40) = mna[30][40] - idt*Geq_M;
	mna[40](30) = mna[40][30] - idt*Geq_M;
	mna[40](40) = mna[40][40] + idt*Geq_M;
	mna[369](40) = mna[369][40] +Geq_M;
	rval[40] += (Ieq_M - Qo_M)*idt;
	fvec[40] += -(Sr_M);
	mna[369](369) = mna[369][369] + 1.;
	rval[369] += Ieq_M;
	fvec[369] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM329*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[33];
	Un_M += Vx[31];
	Uo_M += -uxtold[33];
	Uo_M += uxtold[31];
	Qo_M = uxtold[370];
	Qn_M = Vx[370];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[31](31) = mna[31][31] + idt*Geq_M;
	mna[370](31) = mna[370][31] -Geq_M;
	rval[31] += -(Ieq_M - Qo_M)*idt;
	fvec[31] += Sr_M;
	mna[31](33) = mna[31][33] - idt*Geq_M;
	mna[33](31) = mna[33][31] - idt*Geq_M;
	mna[33](33) = mna[33][33] + idt*Geq_M;
	mna[370](33) = mna[370][33] +Geq_M;
	rval[33] += (Ieq_M - Qo_M)*idt;
	fvec[33] += -(Sr_M);
	mna[370](370) = mna[370][370] + 1.;
	rval[370] += Ieq_M;
	fvec[370] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM330*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[34];
	Un_M += Vx[31];
	Uo_M += -uxtold[34];
	Uo_M += uxtold[31];
	Qo_M = uxtold[371];
	Qn_M = Vx[371];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[31](31) = mna[31][31] + idt*Geq_M;
	mna[371](31) = mna[371][31] -Geq_M;
	rval[31] += -(Ieq_M - Qo_M)*idt;
	fvec[31] += Sr_M;
	mna[31](34) = mna[31][34] - idt*Geq_M;
	mna[34](31) = mna[34][31] - idt*Geq_M;
	mna[34](34) = mna[34][34] + idt*Geq_M;
	mna[371](34) = mna[371][34] +Geq_M;
	rval[34] += (Ieq_M - Qo_M)*idt;
	fvec[34] += -(Sr_M);
	mna[371](371) = mna[371][371] + 1.;
	rval[371] += Ieq_M;
	fvec[371] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM331*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[35];
	Un_M += Vx[31];
	Uo_M += -uxtold[35];
	Uo_M += uxtold[31];
	Qo_M = uxtold[372];
	Qn_M = Vx[372];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[31](31) = mna[31][31] + idt*Geq_M;
	mna[372](31) = mna[372][31] -Geq_M;
	rval[31] += -(Ieq_M - Qo_M)*idt;
	fvec[31] += Sr_M;
	mna[31](35) = mna[31][35] - idt*Geq_M;
	mna[35](31) = mna[35][31] - idt*Geq_M;
	mna[35](35) = mna[35][35] + idt*Geq_M;
	mna[372](35) = mna[372][35] +Geq_M;
	rval[35] += (Ieq_M - Qo_M)*idt;
	fvec[35] += -(Sr_M);
	mna[372](372) = mna[372][372] + 1.;
	rval[372] += Ieq_M;
	fvec[372] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM332*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[35];
	Un_M += Vx[32];
	Uo_M += -uxtold[35];
	Uo_M += uxtold[32];
	Qo_M = uxtold[373];
	Qn_M = Vx[373];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[32](32) = mna[32][32] + idt*Geq_M;
	mna[373](32) = mna[373][32] -Geq_M;
	rval[32] += -(Ieq_M - Qo_M)*idt;
	fvec[32] += Sr_M;
	mna[32](35) = mna[32][35] - idt*Geq_M;
	mna[35](32) = mna[35][32] - idt*Geq_M;
	mna[35](35) = mna[35][35] + idt*Geq_M;
	mna[373](35) = mna[373][35] +Geq_M;
	rval[35] += (Ieq_M - Qo_M)*idt;
	fvec[35] += -(Sr_M);
	mna[373](373) = mna[373][373] + 1.;
	rval[373] += Ieq_M;
	fvec[373] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM333*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[38];
	Un_M += Vx[32];
	Uo_M += -uxtold[38];
	Uo_M += uxtold[32];
	Qo_M = uxtold[374];
	Qn_M = Vx[374];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[32](32) = mna[32][32] + idt*Geq_M;
	mna[374](32) = mna[374][32] -Geq_M;
	rval[32] += -(Ieq_M - Qo_M)*idt;
	fvec[32] += Sr_M;
	mna[32](38) = mna[32][38] - idt*Geq_M;
	mna[38](32) = mna[38][32] - idt*Geq_M;
	mna[38](38) = mna[38][38] + idt*Geq_M;
	mna[374](38) = mna[374][38] +Geq_M;
	rval[38] += (Ieq_M - Qo_M)*idt;
	fvec[38] += -(Sr_M);
	mna[374](374) = mna[374][374] + 1.;
	rval[374] += Ieq_M;
	fvec[374] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM334*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[40];
	Un_M += Vx[32];
	Uo_M += -uxtold[40];
	Uo_M += uxtold[32];
	Qo_M = uxtold[375];
	Qn_M = Vx[375];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[32](32) = mna[32][32] + idt*Geq_M;
	mna[375](32) = mna[375][32] -Geq_M;
	rval[32] += -(Ieq_M - Qo_M)*idt;
	fvec[32] += Sr_M;
	mna[32](40) = mna[32][40] - idt*Geq_M;
	mna[40](32) = mna[40][32] - idt*Geq_M;
	mna[40](40) = mna[40][40] + idt*Geq_M;
	mna[375](40) = mna[375][40] +Geq_M;
	rval[40] += (Ieq_M - Qo_M)*idt;
	fvec[40] += -(Sr_M);
	mna[375](375) = mna[375][375] + 1.;
	rval[375] += Ieq_M;
	fvec[375] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM335*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[37];
	Un_M += Vx[34];
	Uo_M += -uxtold[37];
	Uo_M += uxtold[34];
	Qo_M = uxtold[376];
	Qn_M = Vx[376];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[34](34) = mna[34][34] + idt*Geq_M;
	mna[376](34) = mna[376][34] -Geq_M;
	rval[34] += -(Ieq_M - Qo_M)*idt;
	fvec[34] += Sr_M;
	mna[34](37) = mna[34][37] - idt*Geq_M;
	mna[37](34) = mna[37][34] - idt*Geq_M;
	mna[37](37) = mna[37][37] + idt*Geq_M;
	mna[376](37) = mna[376][37] +Geq_M;
	rval[37] += (Ieq_M - Qo_M)*idt;
	fvec[37] += -(Sr_M);
	mna[376](376) = mna[376][376] + 1.;
	rval[376] += Ieq_M;
	fvec[376] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM336*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[38];
	Un_M += Vx[35];
	Uo_M += -uxtold[38];
	Uo_M += uxtold[35];
	Qo_M = uxtold[377];
	Qn_M = Vx[377];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[35](35) = mna[35][35] + idt*Geq_M;
	mna[377](35) = mna[377][35] -Geq_M;
	rval[35] += -(Ieq_M - Qo_M)*idt;
	fvec[35] += Sr_M;
	mna[35](38) = mna[35][38] - idt*Geq_M;
	mna[38](35) = mna[38][35] - idt*Geq_M;
	mna[38](38) = mna[38][38] + idt*Geq_M;
	mna[377](38) = mna[377][38] +Geq_M;
	rval[38] += (Ieq_M - Qo_M)*idt;
	fvec[38] += -(Sr_M);
	mna[377](377) = mna[377][377] + 1.;
	rval[377] += Ieq_M;
	fvec[377] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM337*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[37];
	Un_M += Vx[36];
	Uo_M += -uxtold[37];
	Uo_M += uxtold[36];
	Qo_M = uxtold[378];
	Qn_M = Vx[378];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[36](36) = mna[36][36] + idt*Geq_M;
	mna[378](36) = mna[378][36] -Geq_M;
	rval[36] += -(Ieq_M - Qo_M)*idt;
	fvec[36] += Sr_M;
	mna[36](37) = mna[36][37] - idt*Geq_M;
	mna[37](36) = mna[37][36] - idt*Geq_M;
	mna[37](37) = mna[37][37] + idt*Geq_M;
	mna[378](37) = mna[378][37] +Geq_M;
	rval[37] += (Ieq_M - Qo_M)*idt;
	fvec[37] += -(Sr_M);
	mna[378](378) = mna[378][378] + 1.;
	rval[378] += Ieq_M;
	fvec[378] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM338*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[38];
	Un_M += Vx[37];
	Uo_M += -uxtold[38];
	Uo_M += uxtold[37];
	Qo_M = uxtold[379];
	Qn_M = Vx[379];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[37](37) = mna[37][37] + idt*Geq_M;
	mna[379](37) = mna[379][37] -Geq_M;
	rval[37] += -(Ieq_M - Qo_M)*idt;
	fvec[37] += Sr_M;
	mna[37](38) = mna[37][38] - idt*Geq_M;
	mna[38](37) = mna[38][37] - idt*Geq_M;
	mna[38](38) = mna[38][38] + idt*Geq_M;
	mna[379](38) = mna[379][38] +Geq_M;
	rval[38] += (Ieq_M - Qo_M)*idt;
	fvec[38] += -(Sr_M);
	mna[379](379) = mna[379][379] + 1.;
	rval[379] += Ieq_M;
	fvec[379] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM339*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[39];
	Un_M += Vx[37];
	Uo_M += -uxtold[39];
	Uo_M += uxtold[37];
	Qo_M = uxtold[380];
	Qn_M = Vx[380];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[37](37) = mna[37][37] + idt*Geq_M;
	mna[380](37) = mna[380][37] -Geq_M;
	rval[37] += -(Ieq_M - Qo_M)*idt;
	fvec[37] += Sr_M;
	mna[37](39) = mna[37][39] - idt*Geq_M;
	mna[39](37) = mna[39][37] - idt*Geq_M;
	mna[39](39) = mna[39][39] + idt*Geq_M;
	mna[380](39) = mna[380][39] +Geq_M;
	rval[39] += (Ieq_M - Qo_M)*idt;
	fvec[39] += -(Sr_M);
	mna[380](380) = mna[380][380] + 1.;
	rval[380] += Ieq_M;
	fvec[380] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM340*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[40];
	Un_M += Vx[37];
	Uo_M += -uxtold[40];
	Uo_M += uxtold[37];
	Qo_M = uxtold[381];
	Qn_M = Vx[381];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[37](37) = mna[37][37] + idt*Geq_M;
	mna[381](37) = mna[381][37] -Geq_M;
	rval[37] += -(Ieq_M - Qo_M)*idt;
	fvec[37] += Sr_M;
	mna[37](40) = mna[37][40] - idt*Geq_M;
	mna[40](37) = mna[40][37] - idt*Geq_M;
	mna[40](40) = mna[40][40] + idt*Geq_M;
	mna[381](40) = mna[381][40] +Geq_M;
	rval[40] += (Ieq_M - Qo_M)*idt;
	fvec[40] += -(Sr_M);
	mna[381](381) = mna[381][381] + 1.;
	rval[381] += Ieq_M;
	fvec[381] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  NM341*/ 

	Un_M = 0.;
	Uo_M = 0.;
	Un_M += -Vx[40];
	Un_M += Vx[39];
	Uo_M += -uxtold[40];
	Uo_M += uxtold[39];
	Qo_M = uxtold[382];
	Qn_M = Vx[382];

	xp_M = exp(-1.1*(Qn_M-10.0));
        gq_M = (1.e1-5.e2)/(1.+100.0*xp_M)+5.e2;
        Kh_M = idt*(Qn_M-Qo_M)*gq_M;
        Kp_M = idt*gq_M + idt*(Qn_M-Qo_M)*(1.1*100.0*(1.e1-5.e2)*xp_M/((1.+100.0*xp_M)*(1.+100.0*xp_M)));

	Geq_M = 1./Kp_M;
	Ieq_M = -Kh_M/Kp_M + Qn_M;
	Sr_M = Un_M + (Ieq_M-Qn_M)/Geq_M;

	mna[39](39) = mna[39][39] + idt*Geq_M;
	mna[382](39) = mna[382][39] -Geq_M;
	rval[39] += -(Ieq_M - Qo_M)*idt;
	fvec[39] += Sr_M;
	mna[39](40) = mna[39][40] - idt*Geq_M;
	mna[40](39) = mna[40][39] - idt*Geq_M;
	mna[40](40) = mna[40][40] + idt*Geq_M;
	mna[382](40) = mna[382][40] +Geq_M;
	rval[40] += (Ieq_M - Qo_M)*idt;
	fvec[40] += -(Sr_M);
	mna[382](382) = mna[382][382] + 1.;
	rval[382] += Ieq_M;
	fvec[382] += -Geq_M*Un_M+Qn_M-Ieq_M;

 /* Elements of  R1*/ 
	mna[13](383) = mna[13][383] + 1.;
	fvec[13] +=  Vx[383];
	mna[383](13) = mna[383][13] + 1.;
	fvec[383] += Vx[13];
	mna[383](383) = mna[383][383] - R[1];
	fvec[383] -= R[1]*Vx[383];

 /* Elements of  R2*/ 
	mna[17](384) = mna[17][384] + 1.;
	fvec[17] +=  Vx[384];
	mna[384](17) = mna[384][17] + 1.;
	fvec[384] += Vx[17];
	mna[384](384) = mna[384][384] - R[2];
	fvec[384] -= R[2]*Vx[384];

 /* Elements of  R3*/ 
	mna[21](385) = mna[21][385] + 1.;
	fvec[21] +=  Vx[385];
	mna[385](21) = mna[385][21] + 1.;
	fvec[385] += Vx[21];
	mna[385](385) = mna[385][385] - R[3];
	fvec[385] -= R[3]*Vx[385];

 /* Elements of  R4*/ 
	mna[25](386) = mna[25][386] + 1.;
	fvec[25] +=  Vx[386];
	mna[386](25) = mna[386][25] + 1.;
	fvec[386] += Vx[25];
	mna[386](386) = mna[386][386] - R[4];
	fvec[386] -= R[4]*Vx[386];

 /* Elements of  R5*/ 
	mna[33](387) = mna[33][387] + 1.;
	fvec[33] +=  Vx[387];
	mna[387](33) = mna[387][33] + 1.;
	fvec[387] += Vx[33];
	mna[387](387) = mna[387][387] - R[5];
	fvec[387] -= R[5]*Vx[387];

 /* Elements of  R6*/ 
	mna[29](388) = mna[29][388] + 1.;
	fvec[29] +=  Vx[388];
	mna[388](29) = mna[388][29] + 1.;
	fvec[388] += Vx[29];
	mna[388](388) = mna[388][388] - R[6];
	fvec[388] -= R[6]*Vx[388];
}