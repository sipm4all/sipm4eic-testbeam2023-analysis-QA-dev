#include <boost/program_options.hpp>
#include <iostream>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"

extern void hough_init(float *cpu_xmap, float *cpu_ymap, float *cpu_rmap, int Nx, int Ny, int Nr);
extern void hough_transform(float *cpu_x, float *cpu_y, float *cpu_h, int cpu_n, int Nx, int Ny, int Nr);
extern void hough_free();

struct program_options_t {
  std::string recodata, ringdata;
};

void
process_program_options(int argc, char *argv[], program_options_t &opt)
{
  /** process arguments **/
  namespace po = boost::program_options;
  po::options_description desc("Options");
  try {
    desc.add_options()
      ("help"             , "Print help messages")
      ("recodata"         , po::value<std::string>(&opt.recodata)->required(), "Reconstructed data input filename")
      ("ringdata"         , po::value<std::string>(&opt.ringdata)->required(), "Ring data output filename")
      ;
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
      std::cout << desc << std::endl;
      exit(1);
    }
  }
  catch(std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    std::cout << desc << std::endl;
    exit(1);
  }
}

int
main(int argc, char *argv[])
{

  program_options_t opt;
  process_program_options(argc, argv, opt);

  /** link to input reconstructed data tree **/
  unsigned short n;
  float x[65534];
  float y[65534];
  float t[65534];
  auto fin = TFile::Open(opt.recodata.c_str());
  auto tin = (TTree *)fin->Get("recodata");
  auto nev = tin->GetEntries();
  tin->SetBranchAddress("n", &n);
  tin->SetBranchAddress("x", &x);
  tin->SetBranchAddress("y", &y);
  tin->SetBranchAddress("t", &t);

  /** create output ring data tree **/
  auto fout = TFile::Open(opt.ringdata.c_str(), "RECREATE");
  auto tout = new TTree("ringdata", "ringdata");
  unsigned short N;
  float X0[256];
  float Y0[256];
  float R[256];
  tout->Branch("N", &N, "N/s");
  tout->Branch("X0", &X0, "X0[N]/F");
  tout->Branch("Y0", &Y0, "Y0[N]/F");
  tout->Branch("R", &R, "R[N]/F");

  /** initialise device **/
  const int Nx = 4;
  const int Ny = 4;
  const int Nr = 16;
  const int Nh = 256 * Nx * Ny * Nr;
  auto xmap = new float[Nh];
  auto ymap = new float[Nh];
  auto rmap = new float [Nh];
  hough_init(xmap, ymap, rmap, Nx, Ny, Nr);

  /** loop over events **/
  auto hough = new float[Nh];
  for (int iev = 0; iev < nev; ++iev) {
    if (iev % 1000 == 0) std::cout << iev << " / " << nev << std::endl;
    tin->GetEntry(iev);

    /** reset ring data **/
    N = 0;
    
    /** hough transform **/
    hough_transform(x, y, hough, n, Nx, Ny, Nr);

    /** get maximum **/
    int imax = std::distance(hough, std::max_element(hough, hough + Nh));
    X0[N] = xmap[imax];
    Y0[N] = ymap[imax];
    R[N] = rmap[imax];

    /** fill tree with ring data **/
    ++N;
    tout->Fill();
  }

  /** free **/
  delete [] xmap;
  delete [] ymap;
  delete [] rmap;
  delete [] hough;
  
  /** free device memory **/
  hough_free();
  
  /** write output and close **/
  fout->cd();
  tout->Write();
  fout->Close();
  fin->Close();

  return 0;
}
