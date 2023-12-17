float r_min = 40.;
float r_max = 90.;
float r_sigma = 2.;
int r_bins = (r_max - r_min) / r_sigma;

float xy_min = -30.;
float xy_max = 30.;
float xy_sigma = 2.;
int xy_bins = (xy_max - xy_min) / xy_sigma;

void
hough(std::string recodata_infilename, std::string ringdata_outfilename, int sev = 0, int nev = kMaxInt, bool display = false)
{

  /** create QA graphs and histograms **/
  
  TCanvas *c = nullptr;
  auto gXY = new TGraph;
  auto gXY_sel = new TGraph;

  auto hXYR = new TH3F("hXYR", ";x (mm);y (mm); R (mm)",
		       xy_bins + 1, xy_min - 0.5 * xy_sigma, xy_max + 0.5 * xy_sigma,
		       xy_bins + 1, xy_min - 0.5 * xy_sigma, xy_max + 0.5 * xy_sigma,
		       r_bins + 1, r_min - 0.5 * r_sigma, r_max + 0.5 * r_sigma);
  
  auto hXY = new TH2F("hMap", ";x (mm);y (mm)",
		       xy_bins + 1, xy_min - 0.5 * xy_sigma, xy_max + 0.5 * xy_sigma,
		      xy_bins + 1, xy_min - 0.5 * xy_sigma, xy_max + 0.5 * xy_sigma);
  
  auto hR = new TH1F("hR", ";x (mm);y (mm)",
		     r_bins + 1, r_min - 0.5 * r_sigma, r_max + 0.5 * r_sigma);
		     
  auto hDR = new TH1F("hDR", ";x (mm);y (mm)", 50, -25., 25.);


  /** link to input reconstructed data tree **/
  unsigned short n;
  float x[65534];
  float y[65534];
  float t[65534];
  auto fin = TFile::Open(recodata_infilename.c_str());
  auto tin = (TTree *)fin->Get("recodata");
  nev = nev < tin->GetEntries() ? nev : tin->GetEntries();
  tin->SetBranchAddress("n", &n);
  tin->SetBranchAddress("x", &x);
  tin->SetBranchAddress("y", &y);
  tin->SetBranchAddress("t", &t);

  /** create output ring data tree **/
  auto fout = TFile::Open(ringdata_outfilename.c_str(), "RECREATE");
  auto tout = new TTree("ringdata", "ringdata");
  unsigned short N;
  float X0[256];
  float Y0[256];
  float R[256];
  tout->Branch("N", &N, "N/s");
  tout->Branch("X0", &X0, "X0[N]/F");
  tout->Branch("Y0", &Y0, "Y0[N]/F");
  tout->Branch("R", &R, "R[N]/F");

  /** loop over events **/
  for (int iev = 0; iev < nev; ++iev) {
    if (iev % 100 == 0)
      std::cout << " --- done " << iev << " / " << nev << " events " << std::endl;
    tin->GetEntry(iev + sev);

    /** reset ring data **/
    N = 0;
    
    /** full 3D transform **/
    hXYR->Reset();
    for (int i = 0 ; i < n; ++i) {
      gXY->SetPoint(i, x[i], y[i]);
      for (int ix = 0; ix < hXYR->GetNbinsX(); ++ix) {
	auto cx = hXYR->GetXaxis()->GetBinCenter(ix + 1);
	for (int iy = 0; iy < hXYR->GetNbinsY(); ++iy) {
	  auto cy = hXYR->GetYaxis()->GetBinCenter(iy + 1);
	  for (int ir = 0; ir < hXYR->GetNbinsZ(); ++ir) {
	    auto r = hXYR->GetZaxis()->GetBinCenter(ir + 1);
	    auto dx = cx - x[i];
	    auto dy = cy - y[i];
	    auto eta = TMath::Sqrt(dx * dx + dy * dy) - r;
	    auto w = TMath::Gaus(eta, 0., 3.5, false);
	    hXYR->Fill(cx, cy, r, w);
	    
	  }}}}
    
    int locmax, locmay, locmaz;
    hXYR->GetMaximumBin(locmax, locmay, locmaz);
    X0[N] = hXYR->GetXaxis()->GetBinCenter(locmax);
    Y0[N] = hXYR->GetYaxis()->GetBinCenter(locmay);
    R[N] = hXYR->GetZaxis()->GetBinCenter(locmaz);
    //    std::cout << " --- after 3D iteration: " << CX << " " << CY << " " << R << std::endl;

    for (int iter = 0; iter < 0; ++iter) {
      
      /** 1D transormation to find radius at fixed center **/
      hR->Reset();
      for (int i = 0 ; i < n; ++i) {
	for (int ir = 0; ir < hR->GetNbinsX(); ++ir) {
	  auto r = hR->GetXaxis()->GetBinCenter(ir + 1);
	  auto dx = X0[N] - x[i];
	  auto dy = Y0[N] - y[i];
	  auto delta = TMath::Sqrt(dx * dx + dy * dy) - r;
	  auto w = TMath::Gaus(delta, 0., r_sigma, false);
	  hR->Fill(r, w);	
	}
      }
      int locmar = hR->GetMaximumBin();
      R[N] = hR->GetBinCenter(locmar);
    
      /** 2D transformation to find center at fixed radius **/
      hXY->Reset();
      for (int i = 0 ; i < n; ++i) {
	for (int ix = 0; ix < hXY->GetNbinsX(); ++ix) {
	  for (int iy = 0; iy < hXY->GetNbinsY(); ++iy) {
	    auto x0 = hXY->GetXaxis()->GetBinCenter(ix + 1);
	    auto y0 = hXY->GetYaxis()->GetBinCenter(iy + 1);
	    auto dx = x0 - x[i];
	    auto dy = y0 - y[i];
	    auto delta = TMath::Sqrt(dx * dx + dy * dy) - R[N];
	    auto w = TMath::Gaus(delta, 0., xy_sigma, false);
	    hXY->Fill(x0, y0, w);
	  }
	}
      }
      int locmax, locmay, locmaz;
      hXY->GetMaximumBin(locmax, locmay, locmaz);
      X0[N] = hXY->GetXaxis()->GetBinCenter(locmax);
      Y0[N] = hXY->GetYaxis()->GetBinCenter(locmay);

      //      std::cout << " --- after 1+2D iteration #" << iter << ": " << X0 << " " << Y0 << " " << R << std::endl;

    }

    /** fit circle **/

    gXY_sel->Set(0);
    for (int i = 0 ; i < n; ++i) {
      auto dx = X0[N] - x[i];
      auto dy = Y0[N] - y[i];
      auto delta = TMath::Sqrt(dx * dx + dy * dy) - R[N];
      if (fabs(delta) < 5.) gXY_sel->AddPoint(x[i], y[i]);
    }
  
    auto chi2 = [&](const double *par) {
      auto np = gXY_sel->GetN();
      double f = 0.;
      auto x = gXY_sel->GetX();
      auto y = gXY_sel->GetY();
      for (int i = 0; i < np; ++i) {
	double dx = x[i] - par[0];
	double dy = y[i] - par[1];
	double delta = TMath::Sqrt(dx * dx + dy * dy) - par[2];
	f += delta * delta;
      }
      return f;
    };
    ROOT::Math::Functor fcn(chi2, 3);
    ROOT::Fit::Fitter fitter;
    double pStart[3] = { X0[N], Y0[N], R[N] };
    fitter.SetFCN(fcn, pStart);
    fitter.Config().ParSettings(0).SetName("X0");
    fitter.Config().ParSettings(1).SetName("Y0");
    fitter.Config().ParSettings(2).SetName("R");
    bool ok = fitter.FitFCN();
    auto result = fitter.Result();
    X0[N] = result.Parameter(0);
    Y0[N] = result.Parameter(1);
    R[N] = result.Parameter(2);
    
    /** event display and QA plots **/    
    if (display) {
      if (!c) {
	c = new TCanvas("c", "c", 800, 800);
	c->Divide(2, 2);
      }
      
      /** prepare QA plots **/
      hXY->Reset();
      hR->Reset();
      hDR->Reset();
      gXY->Set(0);
      gXY_sel->Set(0);
      for (int i = 0 ; i < n; ++i) {
	auto dx = X0[N] - x[i];
	auto dy = Y0[N] - y[i];
	auto delta = TMath::Sqrt(dx * dx + dy * dy) - R[N];
	gXY->SetPoint(i, x[i], y[i]);
	if (fabs(delta) < 5.) gXY_sel->AddPoint(x[i], y[i]);
	hDR->Fill(delta);
	
	/** 1D transormation to find radius at fixed center **/
	for (int ir = 0; ir < hR->GetNbinsX(); ++ir) {
	  auto r = hR->GetXaxis()->GetBinCenter(ir + 1);
	  auto dx = X0[N] - x[i];
	  auto dy = Y0[N] - y[i];
	  auto delta = TMath::Sqrt(dx * dx + dy * dy) - r;
	  auto w = TMath::Gaus(delta, 0., r_sigma, false);
	  hR->Fill(r, w);	
	}
	
	/** 2D transformation of XY **/
	for (int ix = 0; ix < hXY->GetNbinsX(); ++ix) {
	  for (int iy = 0; iy < hXY->GetNbinsY(); ++iy) {
	    auto x0 = hXY->GetXaxis()->GetBinCenter(ix + 1);
	    auto y0 = hXY->GetYaxis()->GetBinCenter(iy + 1);
	    auto dx = x0 - x[i];
	    auto dy = y0 - y[i];
	    auto delta = TMath::Sqrt(dx * dx + dy * dy) - R[N];
	    auto w = TMath::Gaus(delta, 0., xy_sigma, false);
	    hXY->Fill(x0, y0, w);
	  }
	}
      }
      
      /** display **/
      c->cd(1)->DrawFrame(-100., -100., 100., 100., ";x (mm); y (mm)");
      gXY->Draw("samep*");
      gXY_sel->SetMarkerStyle(24);
      gXY_sel->Draw("samep");
      TMarker *marker = new TMarker(X0[N], Y0[N], 2);
      marker->SetMarkerColor(2);
      marker->Draw("same");
      TEllipse *circle = new TEllipse(X0[N], Y0[N], R[N]);
      circle->SetFillStyle(0);
      circle->SetLineColor(kRed);
      circle->Draw("same");
      circle = new TEllipse(X0[N], Y0[N], R[N] - 5.);
      circle->SetFillStyle(0);
      circle->SetLineColor(kRed);
      circle->SetLineStyle(2);    
      circle->Draw("same");
      circle = new TEllipse(X0[N], Y0[N], R[N] + 5.);
      circle->SetFillStyle(0);
      circle->SetLineColor(kRed);
      circle->SetLineStyle(2);    
      circle->Draw("same");
      
      c->cd(2);
      hDR->Draw("histo");
      
      c->cd(3);
      hR->Draw("histo");
      
      c->cd(4);
      hXY->SetStats(false);
      hXY->Draw("colz");
      marker->Draw("same");
      
      c->Update();
      gSystem->Sleep(500);

    } /** and of event display and QA plots **/    

    /** fill tree with ring data **/
    ++N;
    tout->Fill();    
    
  }

  /** write output and close **/
  fout->cd();
  tout->Write();
  fout->Close();
  fin->Close();
}
