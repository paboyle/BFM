#ifndef _BFM_QDP_G5DPARAMS_H_
#define _BFM_QDP_G5DPARAMS_H_

#include <zolotarev.h>

enum BfmPrecisionType{ 
  Bfm32bit,
  Bfm64bit,
  Bfm32bitDefectCorrect
};

enum BfmSolver { 
                 DWF,           // CPS style preconditioning
 		 DWFrb4d,       // Chroma style preconditioning
		 WilsonFermion, // Wilson
		 // Clover,     // Unimplemented
		 WilsonTM,      // Twisted mass wilson
		 // NextNearest neigh Wilson (Eguchi-Kawamoto-Hamber-Wu)
		 WilsonNN,       
		 CloverFermion, //5
		 ////////////////////////////////////
		 // Group of generalised 5d actions
		 // Kerne/Representation/approximation Combos
		 ////////////////////////////////////
		 HwPartFracZolo, // KEK's approach
		 HwContFracZolo, // Edwards, Kennedy et al prefer this
		 HwPartFracTanh, // 
		 HwContFracTanh, // 
		 HwCayleyZolo, // Chiu Optimal
		 HtCayleyZolo, // 
		 HwCayleyTanh, // Scaled shamir
		 HmCayleyTanh, // Scaled shamir 13
		 HtCayleyTanh, // Plain old DWF.
		 DWFTransfer,  // DWF transfer matrix 
		 DWFTransferInv,// DWF transfer inverse matrix 
		 HtContFracTanh,
		 HtContFracZolo
};

// Ugly descriptor for the action
class bfmActionParams {
 public:

  double mass;
  double M5;
  double twistedmass;/*twisted fermion mass ( +/-ig5 mu)*/
  double Csw;        /*Clover coefficient - unimplemented as yet*/
  BfmSolver solver;
  static double mobius_scale;
  double zolo_hi;
  double zolo_lo;
  double zolo_delta; // minimax error. Not a parameter, set by G5D
  int Ls;
  int precon_5d;  
  int solveMobiusDminus;

  bfmActionParams(){
    precon_5d=0;
    mass=Csw=twistedmass=M5=0.0;
    Ls=1;
    zolo_hi=zolo_lo;
    solveMobiusDminus=0;
  }
  void pWilson(double _mass){
    Ls=1;
    precon_5d=0;
    mass=_mass;
    solver =WilsonFermion;
  }
  void pClover(void){
    Ls=1;
    precon_5d=0;
    mass=0.0; // Mass and so on included in A matrices
    solver =CloverFermion;
  }
  void pWilsonTM(double _mass,double _mu){
    Ls=1;
    precon_5d=0;
    mass=_mass;
    twistedmass=_mu;
    solver =WilsonTM;
  }
  void pDWF(double _mass, double _M5,int _Ls){
    mass=_mass; M5=_M5; Ls=_Ls;
    zolo_hi=0;
    zolo_lo=1;
    precon_5d=1;
    mobius_scale = 0.0;
    solver = DWF;
  }
  void pDWFrb4d(double _mass, double _M5,int _Ls){
    mass=_mass; M5=_M5; Ls=_Ls;
    zolo_hi=0;
    zolo_lo=1;
    mobius_scale = 0.0;
    solver = DWFrb4d;
  }

 void ScaledShamirCayleyTanh(double _mass, double _M5, int _Ls,double _mobius_scale){
    precon_5d=0;
    mass=_mass; M5=_M5; Ls=_Ls;
    zolo_hi=0.0;
    zolo_lo=0.0;
    mobius_scale = _mobius_scale;
    solver = HmCayleyTanh;
  }
  void ShamirCayleyTanh(double _mass, double _M5, int _Ls){
    precon_5d=0;
    mass=_mass; M5=_M5; Ls=_Ls;
    zolo_hi=0.0;
    zolo_lo=0.0;
    mobius_scale = 1.0;
    solver = HtCayleyTanh;
  }
  void ShamirCayleyZolo(double _mass, double _M5, int _Ls,double _zolo_lo,double _zolo_hi){
    precon_5d=0;
    mass=_mass; M5=_M5; Ls=_Ls;
    zolo_hi=_zolo_hi;
    zolo_lo=_zolo_lo;
    mobius_scale = 0.0;
    solver = HtCayleyZolo;
  }
  void WilsonCayleyZolo(double _mass, double _M5, int _Ls,double _zolo_lo,double _zolo_hi){
    precon_5d=0;
    mass=_mass; M5=_M5; Ls=_Ls;
    zolo_hi=_zolo_hi;
    zolo_lo=_zolo_lo;
    mobius_scale = 0.0;
    solver = HwCayleyZolo;
  }
  void WilsonCayleyTanh(double _mass, double _M5, int _Ls,double scale){
    precon_5d=0;
    mass=_mass; M5=_M5; Ls=_Ls;
    zolo_hi=0.0;
    zolo_lo=0.0;
    mobius_scale = scale;
    solver = HwCayleyTanh;
  }
  void WilsonContFracZolo(double _mass, double _M5, int _Ls,double _zolo_lo, double _zolo_hi){
    mass=_mass; M5=_M5; Ls=_Ls;
    precon_5d=0;
    zolo_hi=_zolo_hi;
    zolo_lo=_zolo_lo;
    mobius_scale = 0.0;
    solver = HwContFracZolo;
  }
  void WilsonPartFracZolo(double _mass, double _M5, int _Ls,double _zolo_lo, double _zolo_hi){
    precon_5d=0;
    mass=_mass; M5=_M5; Ls=_Ls;
    zolo_hi=_zolo_hi;
    zolo_lo=_zolo_lo;
    mobius_scale = 0.0;
    solver = HwPartFracZolo;
  }
  void VranasMobiusScale(void) {
    mobius_scale = 1+Ls*(1.0/8.0); // Formula from Brower & Vranas  
  }

};
typedef bfmActionParams g5dParams;

#endif
