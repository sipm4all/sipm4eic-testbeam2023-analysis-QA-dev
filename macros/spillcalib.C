void
spillcalib(std::string refinedata_infilename, std::string spillcalib_outfilename)
{
  auto fGaus = (TF1 *)gROOT->GetFunction("gaus");
  std::map<int, TH1 *> hSPILL;

  auto fin = TFile::Open(refinedata_infilename.c_str());
  auto hin = (THnSparse *)fin->Get("hRefine");
  
  /** loop over devices **/
  for (int idevice = 192; idevice < 208; ++idevice) {
    auto di = idevice - 192;
    hin->GetAxis(0)->SetRange(di + 1, di + 1);
    auto hspill = hin->Projection(3, 4);
    if (hspill->Integral() < 100) {
      delete hspill;
      continue;
    }
    
    std::cout << " --- processing device: " << idevice << std::endl;
    hSPILL[idevice] = hspill->ProjectionX(Form("hSPILL_%d", idevice));
    hSPILL[idevice]->Reset();
    for (int ispill = 0; ispill < hspill->GetNbinsX(); ++ispill) {
      auto hpy = hspill->ProjectionY("hpy", ispill + 1, ispill + 1);
      if (hpy->Integral() < 100) {
        continue;
        delete hpy;
      }
      hpy->Rebin(4);
      auto maxbin = hpy->GetMaximumBin();
      auto maxval = hpy->GetBinContent(maxbin);
      auto maxcen = hpy->GetBinCenter(maxbin);
      fGaus->SetParameter(0, maxval);
      fGaus->SetParameter(1, maxcen);
      fGaus->SetParameter(2, 0.5);
      fGaus->SetParLimits(2., 0.1, 1.);
      hpy->Fit(fGaus, "q0B", "", maxcen - 1., maxcen + 1.);

      hSPILL[idevice]->SetBinContent(ispill + 1, std::round(fGaus->GetParameter(1)));
      hSPILL[idevice]->SetBinError(ispill + 1, fGaus->GetParError(1));

      delete hpy;
    }

#if 0
    hspill->GetXaxis()->SetRange(1, 26);
    hspill->Draw("colz");
    hSPILL[idevice]->Draw("same");
    gPad->Update();
    getchar();
#endif
    
    delete hspill;
  }
  
  auto fout = TFile::Open(spillcalib_outfilename.c_str(), "RECREATE");
  for (auto &[idevice, h] : hSPILL)
    h->Write();
  fout->Close();

}
