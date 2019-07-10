#include "NpKNU.hh"
#include <algorithm>

std::vector<std::pair<std::string, int> > errorMSGvec;

ClassImp(npknu::Evt)
ClassImp(npknu::GenInfo)
ClassImp(npknu::P4)
ClassImp(npknu::GenParticle)
ClassImp(npknu::Photon)
ClassImp(npknu::Pileup)
ClassImp(npknu::Vertex)
ClassImp(npknu::Electron)
ClassImp(npknu::Trigger)
ClassImp(npknu::TriggerObject)
ClassImp(npknu::Muon)
ClassImp(npknu::MET)
ClassImp(npknu::Jet)
ClassImp(npknu::PhotonSlim)
ClassImp(npknu::NumIdObjects)

namespace npknu {

void GetHEEPElectronTCA(TClonesArray* inTCA, TClonesArray* outTCA) {
	outTCA->Clear();
	for(int i=0; i < inTCA->GetEntries(); i++) {
		npknu::Electron* elePtr = (npknu::Electron*)inTCA->At(i);
		if(elePtr->vidIsPassHEEP) new ((*outTCA)[(int)outTCA->GetEntries()]) npknu::Electron(*elePtr);
	}
}

void GetTriggerObjectTCA(TClonesArray* inTCA, TClonesArray* outTCA, const int _thisNameIdx) {
	outTCA->Clear();
	for(int i=0; i < inTCA->GetEntries(); i++) {
		npknu::TriggerObject* triObjPtr = (npknu::TriggerObject*)inTCA->At(i);
		if(triObjPtr->IsPassHLT(_thisNameIdx)) new ((*outTCA)[(int)outTCA->GetEntries()]) npknu::TriggerObject(*triObjPtr);
	}
}
void GetTriggerObjectTCA(TClonesArray* inTCA, TClonesArray* outTCA, std::vector<int> _thisNameIdxVec) {
	outTCA->Clear();
	for(int i=0; i < inTCA->GetEntries(); i++) {
		npknu::TriggerObject* triObjPtr = (npknu::TriggerObject*)inTCA->At(i);
		if(triObjPtr->IsPassHLT(_thisNameIdxVec)) new ((*outTCA)[(int)outTCA->GetEntries()]) npknu::TriggerObject(*triObjPtr);
	}
}

void GetTriggerTCA(TClonesArray* inTCA, TClonesArray* outTCA, const int _thisNameIdx) {
	outTCA->Clear();
	for(int i=0; i < inTCA->GetEntries(); i++) {
		npknu::Trigger* triPtr = (npknu::Trigger*)inTCA->At(i);
		if(triPtr->IsPassHLT(_thisNameIdx)) new ((*outTCA)[(int)outTCA->GetEntries()]) npknu::Trigger(*triPtr);
	}
}
void GetTriggerTCA(TClonesArray* inTCA, TClonesArray* outTCA, std::vector<int> _thisNameIdxVec) {
	outTCA->Clear();
	for(int i=0; i < inTCA->GetEntries(); i++) {
		npknu::Trigger* triPtr = (npknu::Trigger*)inTCA->At(i);
		if(triPtr->IsPassHLT(_thisNameIdxVec)) new ((*outTCA)[(int)outTCA->GetEntries()]) npknu::Trigger(*triPtr);
	}
}

//double GetMinDRinTCA(npknu::P4* p4Ptr, TClonesArray* inTCA) {
//	double returnVal = 10000.0;
//	for(int i=0; i<inTCA->GetEntries(); i++) {
//		npknu::P4* ptr = (npknu::P4*)inTCA->At(i);
//		if(p4Ptr->DeltaR(ptr) < returnVal) returnVal = p4Ptr->DeltaR(ptr);
//	}
//	return returnVal;
//}


} // end namespace npknu

namespace npknu {
std::ostream& operator<<(std::ostream& os, const Evt& evt) {
	os << "npknu::Evt " << evt.run << " " << evt.lumi << " " << evt.event << " " << evt.isRealData ;
	return os;
}

std::ostream& operator<<(std::ostream& os, const Evt* evt) {
	os << "npknu::Evt " << evt->run << " " << evt->lumi << " " << evt->event << " " << evt->isRealData ;
	return os;
}

}


void WelcomeNpKNU() { std::cout << "Welcome to npknu " << std::endl; }

