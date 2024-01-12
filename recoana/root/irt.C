#include "invraytrac.C"

void
irt(std::string recodata_infilename)
{

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

  auto hXY = new TH2F("hMap", ";x (mm);y (mm)", 396, -99, 99, 396, -99, 99);
  auto hT = new TH1F("hT", ";t (ns);", 50, -78.125, 78.125);
  auto hAngle = new TH1F("hAngle", ";#vartheta_{Ch} (rad);", 100, 0., 0.5);

  TVector3 E(0., 0., -61.);  // [mm] emission point, nominal (0., 0., -61.)
  TVector3 P(0., 0., 1.);    // [mm] direction vector, nominal (0., 0., 1.)
  TVector3 C(0., 0., -378.); // [mm] centre of curvature, nominal = (0., 0., -378.)
  TVector3 D(0., 0., 0.);    // [mm] detection point
  double R = 700.;           // [mm] radius of aerogel mirror, nominal = 700.
  
  for (int iev = 0; iev < nev; ++iev) {
    tin->GetEntry(iev);
    for (int i = 0 ; i < n; ++i) {

      hXY->Fill(gRandom->Uniform(x[i] - 1.5, x[i] + 1.5), gRandom->Uniform(y[i] - 1.5, y[i] + 1.5));
      hT->Fill(t[i]);

      /** inverse ray-tracking **/
      D.SetXYZ(x[i], y[i], -16.);
      auto S = invraytrac(E, D, C, R);
      S = S - E;
      auto angle = P.Angle(S);
      hAngle->Fill(angle);
      
    } 
  }

  auto c = new TCanvas("c", "c", 1500, 500);
  c->Divide(3, 1);
  c->cd(1)->SetLogz();
  hXY->Draw("col");
  c->cd(2)->SetLogy();
  hT->Draw();
  c->cd(3);
  hAngle->Draw();
  
}
