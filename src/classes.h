#include "DataFormats/Common/interface/Wrapper.h"
#include "MiniAnalyzer/MiniAnalyzer/src/NpKNU.hh"
#include <vector>

namespace {
   struct dictionary {
		std::vector<npknu::Evt> yc_v_evt;
		edm::Wrapper<std::vector<npknu::Evt> > yac_v_evt_edm;
		edm::Wrapper<npknu::Evt> yac_evt_edm;

		std::vector<npknu::P4> yc_v_p4;
		edm::Wrapper<std::vector<npknu::P4> > yac_v_p4_edm;
		edm::Wrapper<npknu::P4> yac_p4_edm;

		std::vector<npknu::GenInfo> yc_v_geninfo;
		edm::Wrapper<std::vector<npknu::GenInfo> > yac_v_geninfo_edm;
		edm::Wrapper<npknu::GenInfo> yac_geninfo_edm;

		std::vector<npknu::GenParticle> yc_v_genparticle;
		edm::Wrapper<std::vector<npknu::GenParticle> > yac_v_genparticle_edm;
		edm::Wrapper<npknu::GenParticle> yac_genparticle_edm;

		std::vector<npknu::Vertex> yc_v_vertex;
		edm::Wrapper<std::vector<npknu::Vertex> > yac_v_vertex_edm;
		edm::Wrapper<npknu::Vertex> yac_vertex_edm;

		std::vector<npknu::Pileup> yc_v_pileup;
		edm::Wrapper<std::vector<npknu::Pileup> > yac_v_pileup_edm;
		edm::Wrapper<npknu::Pileup> yac_pileup_edm;

		std::vector<npknu::Electron> yc_v_electron;
		edm::Wrapper<std::vector<npknu::Electron> > yac_v_electron_edm;
		edm::Wrapper<npknu::Electron> yac_electron_edm;

		std::vector<npknu::Photon> yc_v_photon;
		edm::Wrapper<std::vector<npknu::Photon> > yac_v_photon_edm;
		edm::Wrapper<npknu::Photon> yac_photon_edm;

		std::vector<npknu::MET> yc_v_met;
		edm::Wrapper<std::vector<npknu::MET> > yac_v_met_edm;
		edm::Wrapper<npknu::MET> yac_met_edm;

		std::vector<npknu::Jet> yc_v_jet;
		edm::Wrapper<std::vector<npknu::Jet> > yac_v_jet_edm;
		edm::Wrapper<npknu::Jet> yac_jet_edm;

		std::vector<npknu::Trigger> yc_v_trigger;
		edm::Wrapper<std::vector<npknu::Trigger> > yac_v_trigger_edm;
		edm::Wrapper<npknu::Trigger> yac_trigger_edm;

		std::vector<npknu::TriggerObject> yc_v_triggerobj;
		edm::Wrapper<std::vector<npknu::TriggerObject> > yac_v_triggerobj_edm;
		edm::Wrapper<npknu::TriggerObject> yac_triggerobj_edm;

		std::vector<npknu::Muon> yc_v_muon;
		edm::Wrapper<std::vector<npknu::Muon> > yac_v_muon_edm;
		edm::Wrapper<npknu::Muon> yac_muon_edm;

	};
}
