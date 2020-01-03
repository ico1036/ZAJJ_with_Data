void draw(){

	double XMIN = -1;
	double XMAX = -1;
	long YMAX   = 100;
	
	double YratioMin = 0.1;
	double YratioMax = 2.1;

	//const double Lumi = 0.062138580*1000     ;
	//const double xsec_DYjet = 22880			 ;
		 


	int rebin=1; 
	
	// ---Data
	TFile *fData  = TFile::Open("/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_9_cand2/src/MiniAnalyzer/Analysis/result_data/Data_ele_pho_sel.root") ;

	// ---MC
	//TFile *fDYjet = TFile::Open("/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_9_cand2/src/MiniAnalyzer/Analysis/result_mc/DYjet_ele_sel.root") ;
	//TFile *fDYjet = TFile::Open("/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_9_cand2/src/MiniAnalyzer/Analysis/result_mc/DYjet_ele_sel_Z60_120.root") ;
	//TFile *fDYjet = TFile::Open("/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_9_cand2/src/MiniAnalyzer/Analysis/result_mc/DYjet_ele_sel_Z70_110.root") ;
	
	TFile *fQCD = TFile::Open("/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_9_cand2/src/MiniAnalyzer/Analysis/result_mc/QCDLLAJJ_ele_pho_sel.root") ;

	//TFile *fLLAJJ = TFile::Open("/hcp/data/data02/jwkim2/WORK/CMSSW_9_4_9_cand2/src/MiniAnalyzer/Analysis/result_mc/LLAJJ_ele_sel.root") ;
	
	//TString histname = "h1_Mee"; XMAX=111; XMIN=69; rebin=20; YMAX=100; TString title_name = "Mass_{ee}";
	//TString histname = "h1_e1PT"; XMAX=300; XMIN=0; rebin=10; YMAX=100; TString title_name = "Electron1 p_{T}"; YratioMin = 0.1; YratioMax = 5.1;
	//TString histname = "h1_e2PT"; XMAX=150; XMIN=0; rebin=5; YMAX=50; TString title_name = "Electron2 p_{T}";  YratioMin = 0.1; YratioMax = 5.1;
	TString histname = "h1_phoPT"; XMAX=250; XMIN=0; rebin=20; YMAX=250; TString title_name = "Photon p_{T}"; YratioMin = 0.1; YratioMax = 3.1;


	TH1F *hData		= (TH1F*)fData	  ->Get(histname); 
	TH1F *hQCD		= (TH1F*)fQCD	  ->Get(histname); 
	//TH1F *hDYjet	= (TH1F*)fDYjet	  ->Get(histname); 
	//TH1F *hLLAJJ	= (TH1F*)fLLAJJ	  ->Get(histname); 
	
	cout << "### Before Normalize ###" << endl;
	cout << hData->Integral() << endl;
	cout << hQCD->Integral() << endl;
	//cout << hDYjet->Integral() << endl;
	//cout << hLLAJJ->Integral() << endl;


	hQCD->SetLineWidth(1); hQCD->SetLineColor(kOrange-3); hQCD->SetFillColor(kOrange-3);
	//hDYjet->SetLineWidth(1); hDYjet->SetLineColor(kOrange-3); hDYjet->SetFillColor(kOrange-3);
	//hLLAJJ->SetLineWidth(3); hLLAJJ->SetLineColor(kRed-4);  
	
	//Find bin range
	//double x_start_bin = hDYjet->GetXaxis()->FindBin(60);
	//double x_end_bin = hDYjet->GetXaxis()->FindBin(120);
	//cout << "start bin: " <<x_start_bin << endl;
	//cout << "end bin: "   <<x_end_bin << endl;

	//Normalize
	//hDYjet->Scale(hData->Integral(x_start_bin,x_end_bin) / hDYjet->Integral(x_start_bin,x_end_bin));
	hQCD->Scale(hData->Integral() / hQCD->Integral());
	//hLLAJJ->Scale(hData->Integral() / hLLAJJ->Integral());
	
	//Rebin
	hData->Rebin(rebin);
	hQCD->Rebin(rebin);
	//hDYjet->Rebin(rebin);
	//hLLAJJ->Rebin(rebin);
	
	cout << "### After Normalize ###" << endl;
	cout << hData->Integral() << endl;
	cout << hQCD->Integral() << endl;
	//cout << hDYjet->Integral() << endl;
	//cout << hLLAJJ->Integral() << endl;

	double binwidth = hQCD->GetBinWidth(1);

	TH1F * hRatio = new TH1F(*hData);
	hRatio->Divide(hQCD);


	gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(0);
	gStyle->SetFrameBorderMode(0);

	TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
	TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.0001, 1.0, 1.0);
	pad1->SetBottomMargin(0.01);
	pad1->SetGrid();
	pad1->SetLogy();
	pad1->Draw();
    
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.009, 1, 0.299);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.35);
    pad2->SetGrid();
    pad2->Draw();

		pad1->cd();
		TH2F *null1 = new TH2F("null1","", 2, XMIN, XMAX, 2, 0.09,YMAX);
		null1->GetYaxis()->SetTitle(Form("Number of events / %3.1f GeV",binwidth));
		null1->GetXaxis()->SetTitle(title_name);
		null1->GetYaxis()->SetTitleOffset(1.2);
		null1->GetXaxis()->SetTitleOffset(1.2);
		null1->GetYaxis()->SetTitleSize(0.03);
		null1->GetYaxis()->SetLabelSize(0.03);
		null1->Draw();
		 hQCD->Draw("same hist");
		 hData ->Draw("E1 same");

	

	// --Legend and Latex	
	TLegend *l0 = new TLegend(0.65,0.89,0.90,0.65);
		l0->SetFillStyle(0);
		l0->SetBorderSize(0);
		l0->SetTextSize(0.03);

		TLegendEntry* l01 = l0->AddEntry(hQCD,"QCD LLAJJ"   ,"f"  );	l01->SetTextColor(hQCD->GetLineColor());  
		TLegendEntry* l03 = l0->AddEntry(hData, "Data"    ,"lep"  );			
		l0->Draw();


	pad2->cd();
	TH2F *null2 = new TH2F("null2", "null2", 2, XMIN, XMAX, 2, YratioMin, YratioMax);
	TLine *line = new TLine(XMIN,1,XMAX,1);
	
	null2->SetTitle("");
	null2->GetYaxis()->SetTitle("Data/MC(BKG)");
	null2->GetYaxis()->CenterTitle(true);
	null2->GetYaxis()->SetTitleSize(0.1) ;
	null2->GetYaxis()->SetTitleOffset(0.4);
	null2->GetYaxis()->SetNdivisions(504);
	null2->GetYaxis()->SetLabelSize(0.09) ;
	null2->GetXaxis()->SetLabelSize(0.09) ;
	null2->GetXaxis()->SetTitle(title_name) ;
	null2->GetXaxis()->SetTitleSize(0.1) ;
	null2->GetXaxis()->SetTitleOffset(1.1) ;
	null2->Draw();
	hRatio->SetMarkerStyle(21);
	hRatio->Draw("ep same");
	line->SetLineColor(kRed);
	line->SetLineStyle(9);
	line->Draw();
	 pad1->cd();
	 	TLatex latex;
	 	latex.SetNDC();
	 	latex.SetTextSize(0.04);
	 	latex.SetTextAlign(11);
	 	//latex.DrawLatex(0.6,0.91,Form("%.3f fb^{-1} (13 TeV)", Lumi/1000.0));
	 	
	 
	TString pngname=histname + ".png";
	c1->Print(pngname);
}
