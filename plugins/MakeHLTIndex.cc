#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"


class MakeHLTIndex : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MakeHLTIndex(const edm::ParameterSet&);
      ~MakeHLTIndex();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken;
      std::vector<std::string> triggerLabelsName ;
      int TriggerNameToInt(std::string name, unsigned int runNumber = 0);
        unsigned int oldRun;
};


// --- 1. Constructer: Scan and index my own HLTnames.txt
//		If the HLTnames.txt is null, output is not created 
MakeHLTIndex::MakeHLTIndex(const edm::ParameterSet& iConfig) {

    oldRun = 0;
   triggerResultsToken = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults")); //Get Trigger Result
    triggerLabelsName = iConfig.getParameter<std::vector<std::string> >("triggerLabelsName"); // Get offline TriggerLabelsName( data/HLTnames.txt )
    int ThisGlobalTriggerIndex = 0; // index number 
   for(std::vector<std::string>::iterator it = triggerLabelsName.begin(); it!= triggerLabelsName.end(); it++) { //HLTnames.txt loop start
        if((it->compare(0,1,"#")) == 0) continue;  // If the # is in the first space, continue
      ThisGlobalTriggerIndex++; 
      size_t lastSize = it->find("_v"); // fine the location of "_v"
      if(lastSize != std::string::npos) it->resize(lastSize+2);  // resize the string:  the words after "_v" are deleted
      std::cout << "### IndexTriggerLabelName " << ThisGlobalTriggerIndex << " " << *it << std::endl; 
   }
}

MakeHLTIndex::~MakeHLTIndex() {}
void MakeHLTIndex::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;

// START Event
   unsigned int runNum = (unsigned int) iEvent.id().run();
   unsigned int lumiNum = (unsigned int) iEvent.id().luminosityBlock();
   unsigned long long evtNum = (unsigned long long)iEvent.id().event();
   bool isRealData = iEvent.isRealData();
// END Event

    if(runNum == oldRun) return;
    oldRun = runNum;


// ---2. START Trigger
//		 triggerResults are converted TriggerNames(online) 
//		 output the online TriggerNames
//		 2_1 comparing my HLTnames.txt(offline) and output missing trigger names
   edm::Handle<edm::TriggerResults> triggerResults ; iEvent.getByToken(triggerResultsToken, triggerResults);
   const edm::TriggerNames &names = iEvent.triggerNames(*triggerResults);
    std::cout << "\n === TRIGGER PATHS === " << " " << (unsigned int) iEvent.id().run() << std::endl;
    for (unsigned int i = 0, n = triggerResults->size(); i < n; ++i) {
        int thisTriIdx = TriggerNameToInt(names.triggerName(i), (unsigned int) iEvent.id().run());
        std::cout << "#TriggerNameAll " << " isReal " << isRealData << " RunInfo " << runNum << " " << lumiNum << " " << evtNum << " TriggerName " << names.triggerName(i) << " " << thisTriIdx << std::endl;
    }
// END Trigger
}


// ---2_1 comparing my HLTnames.txt(offline) and output missing trigger names
int MakeHLTIndex::TriggerNameToInt(std::string name, unsigned int runNumber) {
   unsigned int index=0;
   TString nameStr(name);
   for(std::vector<std::string>::iterator it= triggerLabelsName.begin(); it!= triggerLabelsName.end(); it++) {
      index++;
      if(nameStr.Contains(*it)) return index;
   }
   std::cout << "### NotFound TriggerLabel " << name << " Run " << runNumber << std::endl;
   return -1;
}


void MakeHLTIndex::beginJob() {}
void MakeHLTIndex::endJob() {}
void MakeHLTIndex::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(MakeHLTIndex);

