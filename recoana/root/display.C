void
display(std::string recodata_infilename, std::string ringdata_infilename)
{

  /** link to input reconstructed data tree **/
  unsigned short n;
  float x[65534];
  float y[65534];
  float t[65534];
  auto fin = TFile::Open(recodata_infilename.c_str());
  auto tin = (TTree *)fin->Get("recodata");
  auto nev = tin->GetEntries();
  tin->SetBranchAddress("n", &n);
  tin->SetBranchAddress("x", &x);
  tin->SetBranchAddress("y", &y);
  tin->SetBranchAddress("t", &t);

  /** link to input ring data tree **/
  auto frin = TFile::Open(ringdata_infilename.c_str());
  auto trin = (TTree *)frin->Get("ringdata");
  unsigned short N;
  float X0[256];
  float Y0[256];
  float R[256];
  trin->SetBranchAddress("N", &N);
  trin->SetBranchAddress("X0", &X0);
  trin->SetBranchAddress("Y0", &Y0);
  trin->SetBranchAddress("R", &R);

  auto c = new TCanvas("c", "c", 1000, 500);
  c->Divide(2, 1);

  auto gXY = new TGraph;
  gXY->SetMarkerStyle(5);
  auto gXY_sel = new TGraph;
  gXY_sel->SetMarkerStyle(24);
  gXY_sel->SetMarkerColor(kRed+1);
  auto circle = new TEllipse(X0[N], Y0[N], R[N]);
  circle->SetFillStyle(0);
  circle->SetLineColor(kRed);
  auto hDeltaR = new TH1F("hDeltaR", "", 50, -50., 50.);
  auto fGaus = (TF1 *)gROOT->GetFunction("gaus");
  fGaus->SetRange(-50., 50.);
  fGaus->SetNpx(1000);

  /** loop over events **/
  for (int iev = 0; iev < nev; ++iev) {
    tin->GetEntry(iev);
    trin->GetEntry(iev);

    /** fill all hits **/
    gXY->Set(0);
    hDeltaR->Reset();
    for (int i = 0 ; i < n; ++i) {
      auto dx = X0[0] - x[i];
      auto dy = Y0[0] - y[i];
      auto deltaR = TMath::Sqrt(dx * dx + dy * dy) - R[0];
      hDeltaR->Fill(deltaR);
      gXY->AddPoint(x[i], y[i]);
    }

    /** fit deltaR distribution */
    hDeltaR->Fit(fGaus, "IWW0+");
    auto deltaR_mean = fGaus->GetParameter(1);
    auto deltaR_sigma = fGaus->GetParameter(2);

    /** fill selected hits **/
    gXY_sel->Set(0);
    for (int i = 0 ; i < n; ++i) {
      auto dx = X0[0] - x[i];
      auto dy = Y0[0] - y[i];
      auto deltaR = TMath::Sqrt(dx * dx + dy * dy) - R[0];
      if (fabs(deltaR - deltaR_mean) < 3. * deltaR_sigma)
	gXY_sel->AddPoint(x[i], y[i]);
    }

    /** draw **/
    c->cd(1)->DrawFrame(-100., -100., 100., 100., ";x (mm); y (mm)");
    gXY->Draw("samep");
    gXY_sel->Draw("samep");
    circle->DrawEllipse(X0[0], Y0[0], R[0], 0, 0, 360, 0, "same");
    c->cd(2);
    hDeltaR->Draw();
    fGaus->Draw("same");
    
    c->Update();
    gSystem->Sleep(100);
    getchar();
    
  }
  
}
