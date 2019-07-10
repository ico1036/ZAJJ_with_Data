#ifndef NpKNU_HH
#define NpKNU_HH

#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <string>
#include "TObject.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TEfficiency.h"

using std::cout; using std::endl;

namespace npknu {
int CompareStringWild(const char *string, const char *wild);
};

namespace npknu {

class Evt : public TObject {
public:
   unsigned int      run         ; 
   unsigned int      lumi        ; 
   unsigned long int event       ; 
	bool isRealData               ;
	std::vector<double> doubleVec ;
	std::vector<int> intVec       ;
   Evt() {};
   virtual ~Evt() {};
   Evt(unsigned int _run, unsigned int _lumi, unsigned long int _event, bool _isRealData) : run(_run), lumi(_lumi), event(_event), isRealData(_isRealData) {};
  	virtual void Clear(Option_t* option ="") {
		std::vector<double>().swap(doubleVec ) ;
		std::vector<int>().swap(intVec ) ;
	};
	void Print(const char* str = "") {
		std::cout << "NpKNUEventInfo " << run << " " << lumi << " " << event << " " << str << " " ;
	};
	friend std::ostream& operator<<(std::ostream& os, const Evt& evt); 
	friend std::ostream& operator<<(std::ostream& os, const Evt* evt); 
	ClassDef(Evt,1);
};

class P4 : public TObject  {
public:
   double pt     ;
   double eta    ;
   double phi    ;
   double energy ;
	std::vector<double> doubleVec ;
	std::vector<int> intVec       ;

   P4() {};
   virtual ~P4() {};
	virtual void Clear(Option_t* option ="") {
		std::vector<double>().swap(doubleVec ) ;
		std::vector<int>().swap(intVec ) ;
	}; 
	P4(const P4& _p1, const P4& _p2) {
		TLorentzVector P4Vec1 = _p1.GetP4();
		TLorentzVector P4Vec2 = _p2.GetP4();
		pt  = (P4Vec1 + P4Vec2).Pt();
		eta = (P4Vec1 + P4Vec2).Eta();
		phi = (P4Vec1 + P4Vec2).Phi();
		energy = (P4Vec1 + P4Vec2).Energy();
	}
   void SetPtEtaPhiE(double _pt, double _eta, double _phi, double _energy) { pt = _pt; eta = _eta; phi = _phi; energy = _energy;  };
   void SetPxPyPzE(double _px, double _py, double _pz, double _energy) { 
		TLorentzVector p4Vec(_px, _py, _pz, _energy);
		pt = p4Vec.Pt();
		eta = p4Vec.Eta();
		phi = p4Vec.Phi();
		energy = p4Vec.Energy();
	};
	TLorentzVector GetP4() const { TLorentzVector _p4;  _p4.SetPtEtaPhiE(pt, eta, phi, energy); return _p4; }
   double Pt() const { return pt; }
   double Eta() const { return eta; }
   double Phi() const {return phi; }
   double Energy() const { return energy; }
   double E() const { return energy; }
	double AEta() const { return fabs(eta); }
   double DeltaR(const TLorentzVector& _p4) {
      TLorentzVector thisP4 = GetP4();
      return thisP4.DeltaR(_p4);
   };
   double DeltaR(const P4& _p4) {
      TLorentzVector thisP4 = GetP4();
      TLorentzVector otherP4 = _p4.GetP4();
      return thisP4.DeltaR(otherP4);
   };
	double DeltaEta(const P4& _p4) { return std::fabs(eta - _p4.eta) ;  }
	double DeltaPhi(const P4& _p4) {
		TLorentzVector thisP4 = GetP4();
		TLorentzVector otherP4 = _p4.GetP4();
		return fabs(thisP4.DeltaPhi(otherP4));
	}
	double Et() const { return (GetP4().CosTheta() * energy) ; }
   double DeltaR(const P4* _p4Ptr) { return DeltaR(*_p4Ptr); };
   double DeltaEta(const P4* _p4Ptr) { return DeltaEta(*_p4Ptr); };
   double DeltaPhi(const P4* _p4Ptr) { return DeltaPhi(*_p4Ptr); };
	double M(const P4& _p4) { return (GetP4() + _p4.GetP4()).M() ; }
	double M(const P4* _p4Ptr) { return M(*_p4Ptr); }

	double GetMinDRinTCA(TClonesArray* inTCA) {
   	double returnVal = 10000.0;
	   for(int i=0; i<inTCA->GetEntries(); i++) {
   	   npknu::P4* ptr = (npknu::P4*)inTCA->At(i);
      	if(DeltaR(ptr) < returnVal) returnVal = DeltaR(ptr);
   	}
   	return returnVal;
	};
	bool IsPassDRinTCA(const double _MinDR, TClonesArray* inTCA) {
		if(GetMinDRinTCA(inTCA) < _MinDR) return true;
		return false;
	};

   virtual void Print(const char* space = "") const { cout << space << " P4Obj pt eta phi energy : " << Pt() << " " << Eta() << " " << Phi() << " " << Energy() << endl; };
   ClassDef(P4,1);
};

class GenParticle : public P4 {
public:
	unsigned int idx   ;
	int    status      ;
	int    pdgId       ;
	double mass        ;
   unsigned int nDau  ; 
   int mIdx           ;
   GenParticle() {};
   virtual ~GenParticle() {}; 
	virtual void Clear(Option_t* option ="") { P4::Clear("C"); }; 
	void Print() {
		cout << "GenParticle " << idx << " " << pdgId << " " << mass << " " << mIdx << " " << nDau << " " << Pt() << " " << Eta() << " " << Phi() << " " << Energy()  << endl;
	};
   ClassDef(GenParticle,1);
};

class Photon : public P4 {
public:
   double superCluster_eta                                  ;
   double superCluster_phi                                  ;
   double hadTowOverEm                                      ;
   bool   passEleVeto                                       ;
   bool   hasPixelSeed                                      ;

   double full5x5_sigmaIetaIeta                             ;
   double phoChargedIsolationMap                            ;
   double phoESEffSigmaRRMap                                ;
   double phoFull5x5E1x3Map                                 ;
   double phoFull5x5E2x2Map                                 ;
   double phoFull5x5E2x5MaxMap                              ;
   double phoFull5x5E5x5Map                                 ;
   double phoFull5x5SigmaIEtaIEtaMap                        ;
   double phoFull5x5SigmaIEtaIPhiMap                        ;

   double phoNeutralHadronIsolationMap                      ;
   double phoPhotonIsolationMap                             ;
   double phoWorstChargedIsolationMap                       ;

   double chaIso                                            ; 
   double neuIso                                            ; 
   double phoIso                                            ; 
	double sumPUPt                                           ;

   double isoPhotonsWithEA                                  ;
   double isoNeutralHadronsWithEA                           ;
   double isoChargedHadronsWithEA                           ;

   int    matchToTruth                                      ;
   double rho                                               ;
   double effAreaCha                                        ;
   double effAreaNeu                                        ;
   double effAreaPho                                        ;

   bool   vidIsPassTight                                    ;  
   bool   vidIsPassMedium                                   ;  
   bool   vidIsPassLoose                                    ;  
	std::vector<std::pair<int, double> > cutsResultVecMedium ;

	bool   mva80Id                                           ;
	bool   mva90Id                                           ;
	double mvaValue                                          ;
	int    mvaCateg                                          ;

	double EnCorr                                            ;
	int    cutCheckMedium  ;
	
	unsigned int  superCluster_seed_seed_rawId        ;

   double correctedEt      ;
   unsigned int gainSeedSC ;
   double full5x5_r9       ;
   double scale        ;
   double smear        ;
   float scale_stat_up ;
   float scale_stat_dn ;
   float scale_syst_up ;
   float scale_syst_dn ;
   float scale_gain_up ;
   float scale_gain_dn ;
   float resol_rho_up  ;
   float resol_rho_dn  ;
   float resol_phi_up  ;
   float resol_phi_dn  ;

	int    genPdgId;
	double genPt   ;
	double genEta  ;
	double genPhi  ;
	double genEn   ;

   double scEta() { return superCluster_eta; };
   double scAEta() { return std::fabs(superCluster_eta); };
   bool   isGap() { return ((scAEta() > 1.4442) and (scAEta() < 1.566)) ; };
   bool   isEB() { return (scAEta() <= 1.479) ; };
   bool   isEE() { return ((scAEta() > 1.479) and (scAEta() < 2.5)); };

	double PUCha() { return (rho * effAreaCha); }
	double PUNeu() { return (rho * effAreaNeu); }
	double PUPho() { return (rho * effAreaPho); }

	double PUdBeta() { return 0.5 * sumPUPt; }
	double ChaIsoPUdBeta() { return std::max(0.0, phoChargedIsolationMap       - PUCha()) ; }
	double NeuIsoPUdBeta() { return std::max(0.0, phoNeutralHadronIsolationMap - PUNeu()) ; }
	double PhoIsoPUdBeta() { return std::max(0.0, phoPhotonIsolationMap        - PUPho()) ; }

	double CutIsoCha() { return isEB() ? (isoChargedHadronsWithEA / 0.441) : (isoChargedHadronsWithEA / 0.442) ; }
	double CutIsoNeu() { return isEB() ? (isoNeutralHadronsWithEA / (2.725 + 0.0148 * pt + 0.000017 * pt * pt)) : (isoNeutralHadronsWithEA / (1.715 + 0.0163 * pt + 0.000014 * pt * pt)) ; }
	double CutIsoPho() { return isEB() ? (isoPhotonsWithEA / (2.571 + 0.0047 * pt)) : (isoPhotonsWithEA / (3.863 + 0.0034 * pt)) ; }

	void MyPrintCutResults() {
		if(isEB()) {
      	cout << "MyCutPhoEB pt                      0 " << (pt                      > 0.0                                      ) << " " << pt                      << endl;
      	cout << "MyCutPhoEB eta                     1 " << (scAEta()                < 1.479                                    ) << " " << scAEta()                << endl;
      	cout << "MyCutPhoEB HadOverEM               1 " << (hadTowOverEm            < 0.0396                                   ) << " " << hadTowOverEm            << endl;
      	cout << "MyCutPhoEB full5x5_sigmaIetaIeta   3 " << (full5x5_sigmaIetaIeta   < 0.01022                                  ) << " " << full5x5_sigmaIetaIeta   << endl;
      	cout << "MyCutPhoEB isoChargedHadronsWithEA 4 " << (isoChargedHadronsWithEA < 0.441                                    ) << " " << isoChargedHadronsWithEA << endl;
      	cout << "MyCutPhoEB isoNeutralHadronsWithEA 5 " << (isoNeutralHadronsWithEA < 2.725 + 0.0148 * pt + 0.000017 * pt * pt ) << " " << isoNeutralHadronsWithEA << endl;
      	cout << "MyCutPhoEB isoPhotonsWithEA        6 " << (isoPhotonsWithEA        < 2.571 + 0.0047 * pt                      ) << " " << isoPhotonsWithEA        << endl;
		} else {
      	cout << "MyCutPhoEE pt                      0 " << (pt                      > 0.0                                      ) << " " << pt                      << endl;
      	cout << "MyCutPhoEE eta                     1 " << (scAEta()                > 1.479                                    ) << " " << scAEta()                << endl;
      	cout << "MyCutPhoEE HadOverEM               1 " << (hadTowOverEm            < 0.0219                                   ) << " " << hadTowOverEm            << endl;
      	cout << "MyCutPhoEE full5x5_sigmaIetaIeta   3 " << (full5x5_sigmaIetaIeta   < 0.03001                                  ) << " " << full5x5_sigmaIetaIeta   << endl;
      	cout << "MyCutPhoEE isoChargedHadronsWithEA 4 " << (isoChargedHadronsWithEA < 0.442                                    ) << " " << isoChargedHadronsWithEA << endl;
      	cout << "MyCutPhoEE isoNeutralHadronsWithEA 5 " << (isoNeutralHadronsWithEA < 1.715 + 0.0163 * pt + 0.000014 * pt * pt ) << " " << isoNeutralHadronsWithEA << endl;
      	cout << "MyCutPhoEE isoPhotonsWithEA        6 " << (isoPhotonsWithEA        < 3.863 + 0.0034 * pt                      ) << " " << isoPhotonsWithEA        << endl;
		}
	}

	bool MyMediumPUdBetaCut() {
		unsigned int myBit = 0;
		if(isEB()) {
			if(pt                      > 0.0      ) myBit += (1<<0);
			if(scAEta()                < 1.479    ) myBit += (1<<1);
			if(hadTowOverEm            < 0.0396   ) myBit += (1<<2); 
			if(full5x5_sigmaIetaIeta   < 0.01022  ) myBit += (1<<3); 
			if(isoChargedHadronsWithEA < 0.441    ) myBit += (1<<4); 
			if(isoNeutralHadronsWithEA < (2.725 + 0.0148 * pt + 0.000017 * pt * pt) ) myBit += (1<<5); 
			if(isoPhotonsWithEA        < (2.571 + 0.0047 * pt)  ) myBit += (1<<6); 
		} else {
			if(pt                      > 0.0      ) myBit += (1<<0);
			if(scAEta()                > 1.479    ) myBit += (1<<1);
			if(hadTowOverEm            < 0.0219   ) myBit += (1<<2); 
			if(full5x5_sigmaIetaIeta   < 0.03001  ) myBit += (1<<3); 
			if(isoChargedHadronsWithEA < 0.442    ) myBit += (1<<4); 
			if(isoNeutralHadronsWithEA < (1.715 + 0.0163 * pt + 0.000014 * pt * pt) ) myBit += (1<<5); 
			if(isoPhotonsWithEA        < (3.863 + 0.0034 * pt)  ) myBit += (1<<6); 
		}
		if(myBit == (1<<7)-1) return true;
		return false;	
	}
	unsigned int bitMyMediumCut() {
		unsigned int myBit = 0;
		if(isEB()) {
			if(pt                      > 0.0                                      ) myBit |= (1<<0);
			if(scAEta()                < 1.479                                    ) myBit |= (1<<1);
			if(hadTowOverEm            < 0.0396                                   ) myBit |= (1<<2); 
			if(full5x5_sigmaIetaIeta   < 0.01022                                  ) myBit |= (1<<3); 
			if(isoChargedHadronsWithEA < 0.441                                    ) myBit |= (1<<4); 
			if(isoNeutralHadronsWithEA < 2.725 + 0.0148 * pt + 0.000017 * pt * pt ) myBit |= (1<<5); 
			if(isoPhotonsWithEA        < 2.571 + 0.0047 * pt                      ) myBit |= (1<<6); 
		} else {
			if(pt                      > 0.0                                      ) myBit |= (1<<0);
			if(scAEta()                > 1.479                                    ) myBit |= (1<<1);
			if(hadTowOverEm            < 0.0219                                   ) myBit |= (1<<2); 
			if(full5x5_sigmaIetaIeta   < 0.03001                                  ) myBit |= (1<<3); 
			if(isoChargedHadronsWithEA < 0.442                                    ) myBit |= (1<<4); 
			if(isoNeutralHadronsWithEA < 1.715 + 0.0163 * pt + 0.000014 * pt * pt ) myBit |= (1<<5); 
			if(isoPhotonsWithEA        < 3.863 + 0.0034 * pt                      ) myBit |= (1<<6); 
		}
		return myBit;
	}

	unsigned int bitMyTightCut() {
		unsigned int myBit = 0;
		if(isEB()) {
			if(pt                      > 0.0                                      ) myBit |= (1<<0);
			if(scAEta()                < 1.479                                    ) myBit |= (1<<1);
			if(hadTowOverEm            < 0.0269                                   ) myBit |= (1<<2); 
			if(full5x5_sigmaIetaIeta   < 0.00994                                  ) myBit |= (1<<3); 
			if(isoChargedHadronsWithEA < 0.202                                    ) myBit |= (1<<4); 
			if(isoNeutralHadronsWithEA < 0.264 + 0.0148 * pt + 0.000017 * pt * pt ) myBit |= (1<<5); 
			if(isoPhotonsWithEA        < 2.362 + 0.0047 * pt                      ) myBit |= (1<<6); 
		} else {
			if(pt                      > 0.0                                      ) myBit |= (1<<0);
			if(scAEta()                > 1.479                                    ) myBit |= (1<<1);
			if(hadTowOverEm            < 0.0213                                   ) myBit |= (1<<2); 
			if(full5x5_sigmaIetaIeta   < 0.03000                                  ) myBit |= (1<<3); 
			if(isoChargedHadronsWithEA < 0.034                                    ) myBit |= (1<<4); 
			if(isoNeutralHadronsWithEA < 0.586 + 0.0163 * pt + 0.000014 * pt * pt ) myBit |= (1<<5); 
			if(isoPhotonsWithEA        < 2.617 + 0.0034 * pt                      ) myBit |= (1<<6); 
		}
		return myBit;
	}

	bool IsMediumNoSigmaNoIso() {
		unsigned int checkBit = (1<<7) - 1 - (1<<3) - (1<<4) - (1<<5) - (1<<6) ; // remove 
		if((bitMyMediumCut() & checkBit) == checkBit) return true;
		return false;
	}


	bool IsMediumNoSigmaNoCha() {
		unsigned int checkBit = (1<<7) - 1 - (1<<3) - (1<<4) ; // remove 
		if((bitMyMediumCut() & checkBit) == checkBit) return true;
		return false;
	}

	bool IsMediumNoSigmaNoNeu() {
		unsigned int checkBit = (1<<7) - 1 - (1<<3) - (1<<5) ; // remove 
		if((bitMyMediumCut() & checkBit) == checkBit) return true;
		return false;
	}

	bool IsMediumNoSigmaNoPho() {
		unsigned int checkBit = (1<<7) - 1 - (1<<3) - (1<<6) ; // remove 
		if((bitMyMediumCut() & checkBit) == checkBit) return true;
		return false;
	}


	bool IsIdMedium() {
		unsigned int checkBit = (1<<7) - 1; // remove 
		if((bitMyMediumCut() & checkBit) == checkBit) return true;
		return false;
	}

	bool IsTrueMedium() {
		unsigned int checkBit = (1<<7) - 1 - (1<<3); // remove sigmaIEtaIEtaCutBit
		if((bitMyMediumCut() & checkBit) == checkBit) return true;
		return false;
	}
	bool IsFakeMedium() {
		unsigned int checkBit = (1<<7) - 1 - (1<<3) - (1<<4); // remove sigmaIetaIeta and Charged Hadron Isolation;
		if((bitMyMediumCut() & checkBit) == checkBit) {
			return true;
		}
		return false;
	}


	bool IsIdTight() {
		unsigned int checkBit = (1<<7) - 1; 
		if((bitMyTightCut() & checkBit) == checkBit) return true;
		return false;
	}

	bool IsTrueTight() {
		unsigned int checkBit = (1<<7) - 1 - (1<<3); // remove sigmaIEtaIEtaCutBit
		if((bitMyTightCut() & checkBit) == checkBit) return true;
		return false;
	}

	bool IsFakeTight() {
		unsigned int checkBit = (1<<7) - 1 - (1<<3) - (1<<4); // remove sigmaIEtaIEtaCutBit
		if((bitMyTightCut() & checkBit) == checkBit) return true;
		return false;
	}

	bool MyMediumCut() {
		if(bitMyMediumCut() == (1<<7)-1) return true;
		return false;	
	}

	bool IsNtPhotonMedium() {
		unsigned int checkBit = (1<<7) - 1 - (1<<3) - (1<<4); // remove sigmaIetaIeta and ChaHadIso
		if((bitMyMediumCut() & checkBit) == checkBit) { 
			return true;
		}
		return false;
		
	}

	int CompareCutBit() {
		unsigned int myBit = bitMyMediumCut() ;
		unsigned int vidBit = 0;
		unsigned int bitIndex = 0 ;
		for(auto& it : cutsResultVecMedium) {
			if(it.first) { vidBit |= (1<<bitIndex); }
			bitIndex++;
		}
		if(myBit != vidBit) {
			cout << "### DifferentPhotonIdBit MyBit " << myBit << " vidBit " << vidBit << " DiffBit " << ((int)myBit - (int)vidBit) << endl;
			MyPrintCutResults();
			for(auto& it : cutsResultVecMedium) {
				cout << "vidPhoton " << it.first << " " << it.second << endl; 
			}
		}
		return (int)(myBit - vidBit);
	}
   double scDeltaR(double eta, double phi) {
      TLorentzVector p1; p1.SetPtEtaPhiM(1.0, superCluster_eta, superCluster_phi, 0.0);
      TLorentzVector p2; p2.SetPtEtaPhiM(1.0, eta, phi, 0.0);
      return p1.DeltaR(p2);
   };

	Photon() {};
	virtual ~Photon() { ClearVector(); };
	virtual void Clear(Option_t* option ="") { ClearVector(); }; 
	void ClearVector() {
		P4::Clear();
		std::vector<std::pair<int, double> >().swap(cutsResultVecMedium) ;
	}
	ClassDef(Photon,1);
};


class PhotonSlim : public P4 {
public:
   double superCluster_eta        ;
	double hadTowOverEm            ; 
	double full5x5_sigmaIetaIeta   ; 
	double isoChargedHadronsWithEA ; 
	double isoNeutralHadronsWithEA ; 
	double isoPhotonsWithEA        ; 
	bool   hasPixelSeed                 ;
	double phoWorstChargedIsolationMap  ;
   int    matchToTruth     ;
   double full5x5_r9       ;
   bool   vidIsPassTight                                    ;  
   bool   vidIsPassMedium                                   ;  
   bool   vidIsPassLoose                                    ;  
	
   double scEta() { return superCluster_eta; };
   double scAEta() { return std::fabs(superCluster_eta); };
   bool   isGap() { return ((scAEta() > 1.4442) and (scAEta() < 1.566)) ; };
   bool   isEB() { return (scAEta() <= 1.479) ; };
   bool   isEE() { return ((scAEta() > 1.479) and (scAEta() < 2.5)); };

	PhotonSlim() {};
	PhotonSlim(Photon* phoPtr) {
		pt     = phoPtr->pt     ;
		eta    = phoPtr->eta    ;
		phi    = phoPtr->phi    ;
		energy = phoPtr->energy ;

   	superCluster_eta        = phoPtr->superCluster_eta        ;
		hadTowOverEm            = phoPtr->hadTowOverEm            ; 
		full5x5_sigmaIetaIeta   = phoPtr->full5x5_sigmaIetaIeta   ; 
		isoChargedHadronsWithEA = phoPtr->isoChargedHadronsWithEA ; 
		isoNeutralHadronsWithEA = phoPtr->isoNeutralHadronsWithEA ; 
		isoPhotonsWithEA        = phoPtr->isoPhotonsWithEA        ; 
   	matchToTruth            = phoPtr->matchToTruth            ;
   	full5x5_r9              = phoPtr->full5x5_r9              ;
		hasPixelSeed            = phoPtr->hasPixelSeed            ;
		phoWorstChargedIsolationMap  = phoPtr->phoWorstChargedIsolationMap ;

		vidIsPassTight          = phoPtr->vidIsPassTight          ;  
		vidIsPassMedium         = phoPtr->vidIsPassMedium         ;  
		vidIsPassLoose          = phoPtr->vidIsPassLoose          ;  
	};
	unsigned int bitMyMediumCut() {
		unsigned int myBit = 0;
		if(isEB()) {
			if(pt                      > 0.0                                      ) myBit |= (1<<0);
			if(scAEta()                < 1.479                                    ) myBit |= (1<<1);
			if(hadTowOverEm            < 0.0396                                   ) myBit |= (1<<2); 
			if(full5x5_sigmaIetaIeta   < 0.01022                                  ) myBit |= (1<<3); 
			if(isoChargedHadronsWithEA < 0.441                                    ) myBit |= (1<<4); 
			if(isoNeutralHadronsWithEA < 2.725 + 0.0148 * pt + 0.000017 * pt * pt ) myBit |= (1<<5); 
			if(isoPhotonsWithEA        < 2.571 + 0.0047 * pt                      ) myBit |= (1<<6); 
		} else {
			if(pt                      > 0.0                                      ) myBit |= (1<<0);
			if(scAEta()                > 1.479                                    ) myBit |= (1<<1);
			if(hadTowOverEm            < 0.0219                                   ) myBit |= (1<<2); 
			if(full5x5_sigmaIetaIeta   < 0.03001                                  ) myBit |= (1<<3); 
			if(isoChargedHadronsWithEA < 0.442                                    ) myBit |= (1<<4); 
			if(isoNeutralHadronsWithEA < 1.715 + 0.0163 * pt + 0.000014 * pt * pt ) myBit |= (1<<5); 
			if(isoPhotonsWithEA        < 3.863 + 0.0034 * pt                      ) myBit |= (1<<6); 
		}
		return myBit;
	}
	virtual ~PhotonSlim() { };
	ClassDef(PhotonSlim,1);
};


class Muon : public P4 {
public:
   int    type                                                 ;
   double trackIso                                             ;
   double caloIso                                              ;
   double ecalIso                                              ;
   double hcalIso                                              ;

   double isolationR03_emEt                                    ;
   double isolationR03_emVetoEt                                ;
   double isolationR03_hadEt                                   ;
   double isolationR03_hoEt                                    ;
   double isolationR03_hoVetoEt                                ;
   int    isolationR03_nJets                                   ;
   int    isolationR03_nTracks                                 ;
   double isolationR03_sumPt                                   ;
   double isolationR03_trackerVetoPt                           ;

   double isolationR05_emEt                                    ;
   double isolationR05_emVetoEt                                ;
   double isolationR05_hadEt                                   ;
   double isolationR05_hoEt                                    ;
   double isolationR05_hoVetoEt                                ;
   int    isolationR05_nJets                                   ;
   int    isolationR05_nTracks                                 ;
   double isolationR05_sumPt                                   ;
   double isolationR05_trackerVetoPt                           ;

   double pfIsolationR03_sumChargedHadronPt                    ;
   double pfIsolationR03_sumNeutralHadronEt                    ;
   double pfIsolationR03_sumPhotonEt                           ;
   double pfIsolationR03_sumPUPt                               ;
   double pfIsolationR03_sumChargedParticlePt                  ;
   double pfIsolationR03_sumNeutralHadronEtHighThreshold       ;
   double pfIsolationR03_sumPhotonEtHighThreshold              ;

   double pfIsolationR04_sumChargedHadronPt                    ;
   double pfIsolationR04_sumNeutralHadronEt                    ;
   double pfIsolationR04_sumPhotonEt                           ;
   double pfIsolationR04_sumPUPt                               ;
   double pfIsolationR04_sumChargedParticlePt                  ;
   double pfIsolationR04_sumNeutralHadronEtHighThreshold       ;
   double pfIsolationR04_sumPhotonEtHighThreshold              ;

   bool  globalTrack_isNull ;
   bool  innerTrack_isNull  ;

   double globalTrack_dxyPV                                    ;
   double globalTrack_dzPV                                     ;
   double globalTrack_dxyGV                                    ;
   double globalTrack_dzGV                                     ;
   double globalTrack_vx                                       ;
   double globalTrack_vy                                       ;
   double globalTrack_vz                                       ;
   double globalTrack_normalizedChi2                           ;
   double globalTrack_hitPattern_numberOfValidTrackerHits      ;
   double globalTrack_hitPattern_numberOfValidMuonHits         ;
   double globalTrack_validFraction                            ;

   double innerTrack_dxyPV                                     ;
   double innerTrack_dzPV                                      ;
   double innerTrack_dxyGV                                     ;
   double innerTrack_dzGV                                      ;
   double innerTrack_pt                                        ;
   double innerTrack_ptError                                   ;
   double innerTrack_normalizedChi2                            ;
   double innerTrack_hitPattern_trackerLayersWithMeasurement   ;
   double innerTrack_hitPattern_pixelLayersWithMeasurement     ;
   double innerTrack_hitPattern_numberOfValidPixelHits         ;
   double innerTrack_validFraction                             ;

   int    numberOfMatchedStations                              ;
   int    numberOfMatches                                      ;

   bool   isPFMuon                                             ;
   bool   isGlobalMuon                                         ;
   bool   isTrackerMuon                                        ;
   bool   isLooseMuon                                          ;
   bool   isMediumMuon                                         ;
   bool   isTightMuonPV                                        ;
   bool   isTightMuonGV                                        ;
   bool   isSoftMuonPV                                         ;
   bool   isSoftMuonGV                                         ;
   bool   isHighPtMuonPV                                       ;
   bool   isHighPtMuonGV                                       ;

   double combinedQuality_chi2LocalPosition                    ;
   double combindQuality_trkKink                               ;
   double segmentCompatibility                                 ;
   double dB                                                   ;
   double muonBestTrack_dxyPV                                  ;
   double muonBestTrack_dzPV                                   ;
   double muonBestTrack_dxyGV                                  ;
   double muonBestTrack_dzGV                                   ;
   double muonBestTrack_pt                                     ;
   double muonBestTrack_ptError                                ;
   double pfBaseCRIso                                          ;
   double trkBaseIso                                           ;

	bool isBad1 ;
	bool isBad2 ;

	int charge;

	int    genPdgId;
	double genPt   ;
	double genEta  ;
	double genPhi  ;
	double genEn   ;


	Muon() {};
	virtual ~Muon() { ClearVector(); };
	virtual void Clear(Option_t* option ="") { ClearVector(); };

   double MuIsoPF() {
      double b = pfIsolationR04_sumNeutralHadronEt + pfIsolationR04_sumPhotonEt - 0.5 * pfIsolationR04_sumPUPt ;
      double v = pfIsolationR04_sumChargedHadronPt + std::max(0.0, b) ;
      return v / pt;
   }
   double MuIsoTrk() {
      return isolationR03_sumPt / pt;
   }

	bool MediumID2016(bool isBCDEF = false) {
		bool goodGlob =  	isGlobalMuon 
								&& (globalTrack_normalizedChi2 < 3) 
								&& (combinedQuality_chi2LocalPosition < 12) 
								&& (combindQuality_trkKink < 20) ;
		bool isPass =     isLooseMuon 
							&& ( (innerTrack_validFraction > 0.8) || (isBCDEF && innerTrack_validFraction > 0.49) ) 
							&& (segmentCompatibility > (goodGlob ? 0.303 : 0.451) )  ;
		return isPass ;
	}

   double DDeltaR(double eta2, double phi2) {
      TLorentzVector p1; p1.SetPtEtaPhiM(1.0, eta, phi, 0.0);
      TLorentzVector p2; p2.SetPtEtaPhiM(1.0, eta2, phi2, 0.0);
      return p1.DeltaR(p2);
   }
 
	void ClearVector() {
		P4::Clear();
	};
	ClassDef(Muon,1);
};


class Electron : public P4 {
public:
   double gsfTrack_Px                          ;
   double gsfTrack_Py                          ;
   double gsfTrack_Pz                          ;
   double gsfTrack_Pt                          ;
   double superCluster_x                       ;
   double superCluster_y                       ;
   double superCluster_z                       ;
   double superCluster_eta                     ;
   double superCluster_phi                     ;
   double superCluster_seed_eta                ;
   double superCluster_seed_phi                ;

   unsigned int  superCluster_seed_seed_rawId  ;

   double superCluster_energy                  ;
   double ecalEnergy                           ;
   double caloEnergy                           ;

   int    charge                               ;
   double scE1x5                               ;
   double scE2x5Max                            ;
   double scE5x5                               ;
   double scPixCharge                          ;
   double scSigmaEtaEta                        ;
   double scSigmaIEtaIEta                      ;

   bool   ecalDrivenSeed                       ;
   double deltaEtaSuperClusterTrackAtVtx       ;
   double deltaPhiSuperClusterTrackAtVtx       ;
   double ecalPFClusterIso                     ;
   double hcalPFClusterIso                     ;

   double hcalOverEcal                         ;
   double hadronicOverEm                       ;
   double e2x5Max                              ;
   double e5x5                                 ;
   double e1x5                                 ;
   double dr03EcalRecHitSumEt                  ;
   double dr03HcalDepth1TowerSumEt             ;
   double dr03HcalTowerSumEt                   ;
   double dr03TkSumPt                          ;
   double gsfTrack_dxy                         ;
   double gsfTrack_dz                          ;
   double gsfTrack_dxyPVtx                     ;
   double gsfTrack_dzPVtx                      ;
   double gsfTrack_dxyGVtx                     ;
   double gsfTrack_dzGVtx                      ;
   int    nMissingHits                         ;

   double caloIso                              ;
   double ecalIso                              ;
   double hcalIso                              ;
   double trackIso                             ;
   double full5x5_sigmaIetaIeta                ;
   double sigmaIetaIphi                        ;
   double sigmaIphiIphi                        ;
   bool   isPF                                 ;
   bool   passConversionVeto                   ;
   double r9                                   ;
   double fbrem                                ;

   double eSuperClusterOverP                   ;
   double ooEmooP                              ;

   double vtxFitConversion                     ;

   double pfIso_sumChargedHadronPt             ;
   double pfIso_sumNeutralHadronEt             ;
   double pfIso_sumPhotonEt                    ;
   double pfIso_sumPUPt                        ;
   double absIsoWithDBeta                      ;
   double effAreaPFIso                         ;

   double effArea                              ;
   double rho                                  ;

   bool vidIsPassHEEP      ;
	std::vector<std::pair<int, double> > cutsResultVecHEEP      ;


   Electron() {};
   virtual ~Electron() { ClearVector(); };
	virtual void Clear(Option_t* option ="") { ClearVector(); }; 
	void ClearVector() {
		P4::Clear();
		std::vector<std::pair<int, double> >().swap( cutsResultVecHEEP       ) ;
	};

   ClassDef(Electron,1);
};


class Jet : public P4 {
public:
   float neutralHadronEnergyFraction                              ;
   float neutralEmEnergyFraction                                  ;
   float chargedMultiplicity                                      ;
   float neutralMultiplicity                                      ;
   float muonEnergyFraction                                       ;
   float chargedHadronEnergyFraction                              ;
   float chargedEmEnergyFraction                                  ;

   bool passPFLooseId                                             ;
   bool passPFTightId                                             ;
   bool passPFTightLVetoId                                        ;

   double   bDiscri_pfTrackCountingHighEffBJetTags                ;
   double   bDiscri_pfTrackCountingHighPurBJetTags                ;
   double   bDiscri_pfJetProbabilityBJetTags                      ;
   double   bDiscri_pfJetBProbabilityBJetTags                     ;
   double   bDiscri_pfSimpleSecondaryVertexHighEffBJetTags        ;
   double   bDiscri_pfSimpleSecondaryVertexHighPurBJetTags        ;
   double   bDiscri_pfCombinedSecondaryVertexV2BJetTags           ;
   double   bDiscri_pfCombinedInclusiveSecondaryVertexV2BJetTags  ;
   double   bDiscri_pfCombinedMVABJetTags                         ;

   double ptRaw                                                   ;
   double uncRaw                                                  ;
   double uncCor                                                  ;
   double ptUncCorDo                                              ;
   double ptUncCorUp                                              ;

	int    genPdgId;
	double genPt   ;
	double genEta  ;
	double genPhi  ;
	double genEn   ;

   double genParton_pt                                            ;
   int    genHadronFlavour                                        ;
   int    genPartonFlavour                                        ;
   std::vector<double> jecFactorVec                               ;
   unsigned int bitCutCheck                                       ;
   std::vector<std::pair<int, int> > genStatusPidVec        ;
   std::vector<P4>                   genP4Vec               ;
	Jet() {};
   virtual ~Jet() { ClearVector(); };
   virtual void Clear(Option_t* option ="") { ClearVector(); };
   void ClearVector() {
		P4::Clear();
      std::vector<double>              ().swap( jecFactorVec   )  ;
      std::vector<std::pair<int, int> >().swap( genStatusPidVec)  ;
      std::vector<P4>                  ().swap( genP4Vec       )  ;
   };
   bool cutPFLooseID() {
      double NHF                 =  neutralHadronEnergyFraction    ;
      double NEMF                =  neutralEmEnergyFraction        ;
      double CHF                 =  chargedHadronEnergyFraction    ;
      double CEMF                =  chargedEmEnergyFraction        ;
      double NumConst            =  chargedMultiplicity + neutralMultiplicity ;
      double NumNeutralParticle  =  neutralMultiplicity            ;
      double CHM                 =  chargedMultiplicity            ;
      bool looseJetID1 = ((NHF<0.99 && NEMF<0.99 && NumConst>1) && ((std::abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || std::abs(eta)>2.4) && std::abs(eta)<=2.7)    ;
      bool looseJetID2 =  (NHF<0.98 && NEMF>0.01 && NumNeutralParticle>2 && std::abs(eta)>2.7 && std::abs(eta)<=3.0 ) ;
      bool looseJetID3 =  (NEMF<0.90 && NumNeutralParticle>10 && std::abs(eta)>3.0 ) ;
      return (looseJetID1 || looseJetID2 || looseJetID3) ;
   };
   bool cutPFTightID() {
      double NHF                 =  neutralHadronEnergyFraction    ;
      double NEMF                =  neutralEmEnergyFraction        ;
      double CHF                 =  chargedHadronEnergyFraction    ;
      double CEMF                =  chargedEmEnergyFraction        ;
      double NumConst            =  chargedMultiplicity + neutralMultiplicity ;
      double NumNeutralParticle  =  neutralMultiplicity            ;
      double CHM                 =  chargedMultiplicity            ;
      bool tightJetID1 = ((NHF<0.90 && NEMF<0.90 && NumConst>1) && ((std::abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || std::abs(eta)>2.4) && std::abs(eta)<=2.7) ;
      bool tightJetID2 =  (NHF<0.98 && NEMF>0.01 && NumNeutralParticle>2 && std::abs(eta)>2.7 && std::abs(eta)<=3.0 ) ;
      bool tightJetID3 =  (NEMF<0.90 && NumNeutralParticle>10 && std::abs(eta)>3.0 ) ;
      return (tightJetID1 || tightJetID2 || tightJetID3);
   };
   bool cutPFTightLVetoID() {
      double NHF                 =  neutralHadronEnergyFraction    ;
      double NEMF                =  neutralEmEnergyFraction        ;
      double CHF                 =  chargedHadronEnergyFraction    ;
      double MUF                 =  muonEnergyFraction             ;
      double CEMF                =  chargedEmEnergyFraction        ;
      double NumConst            =  chargedMultiplicity + neutralMultiplicity ;
      double CHM                 =  chargedMultiplicity            ;
      bool tightLepVetoJetID = ((NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((std::abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || std::abs(eta)>2.4) && std::abs(eta)<=2.7) ;
      return tightLepVetoJetID ;
   };
	ClassDef(Jet,1)	
};


class MET : public P4 {
public:
   bool   isCaloMET             ;
   bool   isPFMET               ;
   bool   isRecoMET             ;
	double pt                   ;
	double phi                  ;
   double sumEt                ; 
   double genMET_pt            ; 
   double genMET_phi           ; 
   double genMET_sumEt         ; 
   double shiftedPt_JetEnUp    ;         
   double shiftedPt_JetEnDown  ;     
	double metSignificance      ;
   MET() {};
   virtual ~MET() {};
	virtual void Clear(Option_t* option ="") { 
		P4::Clear();
	}; 
   ClassDef(MET,1)
};


class Trigger : public TObject  {
public:
   int nameIdx             ;
   std::string name        ;
   bool accept             ;
   double prescaleForIndex ;
   Trigger() {};
   virtual ~Trigger() {};
   virtual void Clear(Option_t* option ="") { };

//   bool Pass(int index) { return ((index == nameIndex) && (accept)) ; };
//   bool Pass(std::vector<int> indexVec) {
//      for(std::vector<int>::iterator it=indexVec.begin(); it!=indexVec.end(); it++) { if(Pass(*it)) return true; }
//      return false;
//   }
	bool IsPassHLT(const int _idx) {
		if(_idx == nameIdx and accept) return true;
		return false;
	};
	bool IsPassHLT(std::vector<int> _idxVec) {
		if(std::find(_idxVec.begin(), _idxVec.end(), nameIdx) != _idxVec.end() and accept) return true;
		return false;
	};
   ClassDef(Trigger,1)
};

class TriggerObject : public P4 {
public:
   std::vector<std::string>  filterLabelVec    ;
   std::vector<int>  pathNameIdxVec             ;
   std::vector<std::string>  pathNameVec        ;
   std::vector<unsigned int> pathBitVec         ;
   TriggerObject() {};
   virtual ~TriggerObject() { ClearVector(); };
   virtual void Clear(Option_t* option ="") { ClearVector(); };
   void ClearVector() {
		P4::Clear();
      std::vector<std::string>        ().swap(filterLabelVec)    ;
      std::vector<std::string>        ().swap(pathNameVec   )    ;
      std::vector<unsigned int>       ().swap(pathBitVec    ) ;
   };
	void Print(bool doAll = false) {
		P4::Print();
		std::cout << "TriggerObject pathNameIdxPassHLT " << " : " ;
		for(unsigned int i=0; i<pathNameIdxVec.size(); i++) {
			if(pathBitVec[i] == 15) {
				cout << pathNameIdxVec[i] << " " ;
			}
		}
		std::cout << std::endl;
		if(doAll) {
			std::cout << "TriggerObject pathNameIdxAll " << pathNameIdxVec.size() << " : " ;
			std::copy(pathNameIdxVec.begin(), pathNameIdxVec.end(), std::ostream_iterator<int>(std::cout, " "));
			std::cout << std::endl;
			std::cout << "TriggerObject pathBitAll " << pathBitVec.size() << " : " ;
			std::copy(pathBitVec.begin(), pathBitVec.end(), std::ostream_iterator<unsigned int>(std::cout, " "));
			std::cout << std::endl;
		}
	};
	bool IsPassHLT(const int _idx) {
		for(unsigned int i=0; i<pathNameIdxVec.size(); i++) {
			if((pathNameIdxVec[i] == _idx) and (pathBitVec[i] == 15)) return true;
		}
		return false;
	};
	bool IsPassHLT(std::vector<int> _idxVec) {
		for(unsigned int i=0; i<pathNameIdxVec.size(); i++) {
			if((std::find(_idxVec.begin(), _idxVec.end(), pathNameIdxVec[i]) != _idxVec.end()) and (pathBitVec[i] == 15)) return true;
		}
		return false;
	};
   ClassDef(TriggerObject,1)
};



class GenInfo : public TObject  {
public:
	double weight   ;
   double scalePDF ;
   int pdg1        ;
   int pdg2        ;
   double x1       ;
   double x2       ;
   double xpdf1    ;
   double xpdf2    ;
  	double qScale   ;
  	double alphaQCD ;
  	double alphaQED ;

	GenInfo() {};
	virtual ~GenInfo() {};
	virtual void Clear(Option_t* option ="") { }; 
	ClassDef(GenInfo,1);
};

class Pileup : public TObject  {
public:
	int    BunchCrossing       ;
	double TrueNumInteractions ;
	int    PU_NumInteractions  ;
   Pileup() {};
   virtual ~Pileup() {}; 
	virtual void Clear(Option_t* option ="") { }; 
   ClassDef(Pileup,1);
};

class Vertex : public TObject {
public:
	double x              ;  
	double y              ;  
	double z              ;  
	double xError         ;  
	double yError         ;  
	double zError         ;  
	int    tracksSize     ;  
	int    nTracks        ;  
	bool   isFake         ;  
	double ndof           ;  
	double position_rho   ;  
	double chi2           ;  
	double normalizedChi2 ;  

   Vertex() {};
   virtual ~Vertex() {}; 
	virtual void Clear(Option_t* option ="") { }; 
	bool isGoodVertex() { //https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_7.4.12/ElectronNtupler/plugins/ElectronNtuplerVIDDemo.cc
		return ((!isFake) and (ndof >= 4.0) and (position_rho <= 2.0) and (std::fabs(z) <= 24.0));
   };
   ClassDef(Vertex,1);
};

class NumIdObjects : public TObject {
public:
   int numGenPho    ; 
   int numGenPho10  ; 
   int numGenPho20  ; 
   int numGenPho25  ; 
   int numSelEleNC1 ; 
   int numSelEleNC2 ; 
   int numSelEleNC3 ; 
   int numSelEleNC4 ; 
   int numSelEleNC5 ; 
   int numSelEleNC6 ; 
   int numSelEleEC1 ; 
   int numSelEleEC2 ; 
   int numSelEleEC3 ; 
   int numSelEleEC4 ; 
   int numSelEleEC5 ; 
   int numSelEleEC6 ; 
   int numSelPhoNC1 ; 
   int numSelPhoNC2 ; 
   int numSelPhoNC3 ; 
   int numSelPhoNC4 ; 
   int numSelPhoNC5 ; 
   int numSelPhoEC1 ; 
   int numSelPhoEC2 ; 
   int numSelPhoEC3 ; 
   int numSelPhoEC4 ; 
   int numSelPhoEC5 ; 

   NumIdObjects() {};
   virtual ~NumIdObjects() {}; 
	virtual void Clear(Option_t* option ="") { }; 
   ClassDef(NumIdObjects,1);
};


class PUReweight {
private:
   long unsigned int events    ;
   double sumWeight   ;
public:
   TH1D* h1_data  ;
   TH1F* h1_mc    ;
   TH1D* PUweightHist ;
   PUReweight(TString fileName, TString fileNameMC, TString mcHistName) {
      std::cout << "### dataPUfile " << fileName << std::endl;
      std::cout << "### mcPUfile " << fileNameMC << std::endl;
      std::cout << "### mcPUhist " << mcHistName << std::endl;
      events = 0;
      sumWeight = 0;
      TFile* pileupFile = new TFile(fileName,"READ");
      PUweightHist = (TH1D*)pileupFile->Get("pileup");
      PUweightHist->SetDirectory(0);
      pileupFile->Close();
      double PUweightInt = PUweightHist->Integral();

      TH1F* mcPU=NULL;
      TFile* mcFile = new TFile(fileNameMC,"READ");
      if( mcPU==NULL) mcPU = (TH1F*)mcFile->Get(mcHistName);
      mcPU->SetDirectory(0);
      mcFile->Close();

      PUweightHist->Scale(1.0/PUweightInt);
      mcPU->Scale(1.0/mcPU->Integral());
      PUweightHist->Divide(mcPU);

      delete mcPU;
   };
   virtual ~PUReweight() {};
   double getWeight(double nPU) {
      events++;
      double PUweight = 0.0;
      PUweight = PUweightHist->GetBinContent(PUweightHist->GetXaxis()->FindBin(nPU));
      sumWeight+= PUweight;
      return PUweight;
   };
   double getAvgWeight() { return (sumWeight / events); };
};


typedef std::vector<npknu::GenInfo>       GenInfoCollection        ;
typedef std::vector<npknu::GenParticle>   GenParticleCollection    ;
typedef std::vector<npknu::Pileup>        PileupCollection         ;
typedef std::vector<npknu::Vertex>        VertexCollection         ;
typedef std::vector<npknu::Muon>          MuonCollection           ;
typedef std::vector<npknu::Electron>      ElectronCollection       ;
typedef std::vector<npknu::Photon>        PhotonCollection         ;
typedef std::vector<npknu::MET>           METCollection            ;
typedef std::vector<npknu::Jet>           JetCollection            ;
typedef std::vector<npknu::Trigger>       TriggerCollection        ;
typedef std::vector<npknu::TriggerObject> TriggerObjectCollection  ;

}; // End namespace npknu 

namespace npknu {
void GetHEEPElectronTCA(TClonesArray* inTCA, TClonesArray* outTCA);
void GetTriggerTCA(TClonesArray* inTCA, TClonesArray* outTCA, const int _thisNameIdx);
void GetTriggerTCA(TClonesArray* inTCA, TClonesArray* outTCA, std::vector<int> _thisNameIdxVec);
void GetTriggerObjectTCA(TClonesArray* inTCA, TClonesArray* outTCA, const int _thisNameIdx);
void GetTriggerObjectTCA(TClonesArray* inTCA, TClonesArray* outTCA, std::vector<int> _thisNameIdxVec);
//double GetMinDRinTCA(npknu::P4* p4Ptr, TClonesArray* inTCA);

int ElectronIdTCA(TClonesArray* inTCA, TClonesArray* outTCA, int type);
int PhotonIdTCA(TClonesArray* inTCA, TClonesArray* outTCA, int type);

int ElectronSelection(TClonesArray* inTCA, TClonesArray* outTCA, double pt1Cut, double pt2Cut, int type = 3);
int NumVetoEleSelection(TClonesArray* inTCA, TClonesArray* outTCA, double ptCut, npknu::Electron* ele1Ptr, npknu::Electron* ele2Ptr);
int NumVetoMuonSelection(TClonesArray* inTCA, TClonesArray* outTCA, double ptCut);
//void NtPhotonSelection(TClonesArray* inTCA, TClonesArray* outTCA, double pt1Cut, int type = 2);
int PhotonSelection(TClonesArray* inTCA, TClonesArray* outTCA, double pt1Cut, int type = 2, npknu::P4* p1Ptr =NULL, npknu::P4* p2Ptr = NULL);
void NtJetSelection(TClonesArray* inTCA, TClonesArray* outTCA, double ptCut);
double UnCorrPt(npknu::P4* p4Ptr, TClonesArray* inTCA);


void SortElectronPt(TClonesArray* inTCA);
void SortPhotonPt(TClonesArray* inTCA);
void SortJetPt(TClonesArray* inTCA);
double MCWeight(TClonesArray* genInfoTCA);
bool PassingTrigger(int id, TClonesArray* inTCA);
int NumGoodVertex(TClonesArray* inTCA);
double GetPileupTrue(TClonesArray* inTCA);
double MinDeltaR(TClonesArray* inTCA);
void RebinTH1F(TH1F* h1_in, TH1F* h1_out);

// class NpKNURun {
// 
// 	NpKNURun(TChain* inChain) {
// 	}
// 	void Run() 
// };

};


void WelcomeNpKNU();


#endif
