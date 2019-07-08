void draw(){

	double XMIN = -1;
	double XMAX = -1;
	long YMAX   = 100;
	
	//const double Lumi = 0.062138580*1000     ;
	//const double xsec_DYjet = 22880			 ;
		 



	int rebin=1; 
	TFile *fData  = TFile::Open("hist/DoubleEG_GT_Run2016B/hist_Data_Z60_120.root") ;
	TFile *fDYjet = TFile::Open("hist/MC/hist_DYjet_Z60_120.root") ;

	
	TString histname = "h1_Mee"; XMAX=120; XMIN=60; rebin=10; YMAX=10000; TString title_name = "Mass_{ee}";


	TH1F *hData		= (TH1F*)fData	  ->Get(histname); 
	TH1F *hDYjet	= (TH1F*)fDYjet	  ->Get(histname); 
	
	cout << "### Before Normalize ###" << endl;
	cout << hData->Integral() << endl;
	cout << hDYjet->Integral() << endl;


	//hDYjet->SetLineWidth(3); hDYjet->SetLineColor(46);  hDYjet->Scale(Lumi * xsec_DYjet / 25000);
	hDYjet->SetLineWidth(3); hDYjet->SetLineColor(46);  
	
	//Find bin range
	double x_start_bin = hDYjet->GetXaxis()->FindBin(60);
	double x_end_bin = hDYjet->GetXaxis()->FindBin(120);
	
	cout << "start bin: " <<x_start_bin << endl;
	cout << "end bin: "   <<x_end_bin << endl;

	//Normalize
	hDYjet->Scale(hData->Integral(x_start_bin,x_end_bin) / hDYjet->Integral(x_start_bin,x_end_bin));
	//hDYjet->Scale(hData->Integral() / hDYjet->Integral());
	hDYjet->Rebin(rebin);
	
	hData->Rebin(rebin);
	
	cout << "### After Normalize ###" << endl;
	cout << hData->Integral() << endl;
	cout << hDYjet->Integral() << endl;

	double binwidth = hDYjet->GetBinWidth(1);

	gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(0);
	gStyle->SetFrameBorderMode(0);

	TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
	TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.0001, 1.0, 1.0);
		//   pad1->SetBottomMargin(0.01);
		pad1->SetGrid();
		pad1->SetLogy();
		pad1->Draw();
		pad1->cd();
		TH2F *null1 = new TH2F("null1","", 2, XMIN, XMAX, 2, 0.09,YMAX);
		null1->GetYaxis()->SetTitle(Form("Number of events / %3.1f GeV",binwidth));
		null1->GetXaxis()->SetTitle(title_name);
		null1->GetYaxis()->SetTitleOffset(1.2);
		null1->GetXaxis()->SetTitleOffset(1.2);
		null1->GetYaxis()->SetTitleSize(0.03);
		null1->GetYaxis()->SetLabelSize(0.03);
		null1->Draw();
		 hDYjet->Draw("same");
		 hData ->Draw("E1 same");

// --Legend and Latex	
	TLegend *l0 = new TLegend(0.65,0.89,0.90,0.65);
		l0->SetFillStyle(0);
		l0->SetBorderSize(0);
		l0->SetTextSize(0.03);

		TLegendEntry* l01 = l0->AddEntry(hDYjet,"DYjet"   ,"l"  );	l01->SetTextColor(hDYjet->GetLineColor());  
		TLegendEntry* l02 = l0->AddEntry(hData, "Data"    ,"lep"  );			
		l0->Draw();

	pad1->cd();
		TLatex latex;
		latex.SetNDC();
		latex.SetTextSize(0.04);
		latex.SetTextAlign(11);
		//latex.DrawLatex(0.6,0.91,Form("%.3f fb^{-1} (13 TeV)", Lumi/1000.0));
		
		TString pngname=histname + ".png";
		c1->Print("MCvsDAta6_"+pngname);

}























