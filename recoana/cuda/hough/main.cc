#include <boost/program_options.hpp>
#include <iostream>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "common.h"

extern void hough_init(data_t data);
extern void hough_transform(data_t data);
extern void hough_free();

struct program_options_t {
  std::string recodata, ringdata;
  float xmin, xmax, ymin, ymax, rmin, rmax, tmin, tmax;
  int xbins, ybins, rbins, tbins;
  float tsigma;
};

/**
const float x_min = -7.5;
const float x_stp = 1.;
const float y_min = -7.5;
const float y_stp = 1.;
const float r_min = 32.;
const float r_stp = 1.;
const float t_min = -787.5;
const float t_stp = 25.;
**/

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
      ("xmin"             , po::value<float>(&opt.xmin)->default_value(-10.), "Hough space xmin")
      ("xmax"             , po::value<float>(&opt.xmax)->default_value(10.), "Hough space xmax")
      ("xbins"            , po::value<int>(&opt.xbins)->default_value(20), "Hough space xbins")
      ("ymin"             , po::value<float>(&opt.ymin)->default_value(-10.), "Hough space ymin")
      ("ymax"             , po::value<float>(&opt.ymax)->default_value(10.), "Hough space ymax")
      ("ybins"            , po::value<int>(&opt.ybins)->default_value(20), "Hough space ybins")
      ("rmin"             , po::value<float>(&opt.rmin)->default_value(30.), "Hough space rmin")
      ("rmax"             , po::value<float>(&opt.rmax)->default_value(100.), "Hough space rmax")
      ("rbins"            , po::value<int>(&opt.rbins)->default_value(70), "Hough space rbins")
      ("tmin"             , po::value<float>(&opt.tmin)->default_value(-800), "Hough space tmin")
      ("tmax"             , po::value<float>(&opt.tmax)->default_value(800.), "Hough space tmax")
      ("tbins"            , po::value<int>(&opt.tbins)->default_value(32), "Hough space tbins")
      ("tsigma"            , po::value<float>(&opt.tsigma)->default_value(15.), "Hough space tsigma")
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
  float T[256];
  tout->Branch("N", &N, "N/s");
  tout->Branch("X0", &X0, "X0[N]/F");
  tout->Branch("Y0", &Y0, "Y0[N]/F");
  tout->Branch("R", &R, "R[N]/F");
  tout->Branch("T", &T, "T[N]/F");

  /** initialise **/
  data_t data;
  data.min.x = opt.xmin;
  data.max.x = opt.xmax;
  data.bins.x = opt.xbins;
  data.min.y = opt.ymin;
  data.max.y = opt.ymax;
  data.bins.y = opt.ybins;
  data.min.r = opt.rmin;
  data.max.r = opt.rmax;
  data.bins.r = opt.rbins;
  data.min.t = opt.tmin;
  data.max.t = opt.tmax;
  data.bins.t = opt.tbins;
  data.sigma.t = opt.tsigma;
  
  const int size = data.bins.x * data.bins.y * data.bins.r * data.bins.t;
  const int grid_size = 1 + (size - 1) / 256;
  data.map.x = new float[size];
  data.map.y = new float[size];
  data.map.r = new float[size];
  data.map.t = new float[size];
  data.hough.h = new float[size];
  data.hough.rh = new float[grid_size];
  data.hough.rhi = new int[grid_size];
  hough_init(data);

  /** loop over events **/
  for (int iev = 0; iev < nev; ++iev) {
    tin->GetEntry(iev);
    if (iev % 1000 == 0)
      std::cout << "processing event: " << iev << " / " << nev << std::endl;
    
    /** set data points **/
    data.points.n = n;
    data.points.x = x;
    data.points.y = y;
    data.points.t = t;
    
    /** reset ring data **/
    N = 0;
    
    /** hough transform **/
    hough_transform(data);

    /** get maximum **/
    int rimax = std::distance(data.hough.rh, std::max_element(data.hough.rh, data.hough.rh + grid_size));
    int imax = data.hough.rhi[rimax];
    X0[N] = data.map.x[imax];
    Y0[N] = data.map.y[imax];
    R[N] = data.map.r[imax];
    T[N] = data.map.t[imax];

    /** fill tree with ring data **/
    ++N;
    tout->Fill();
  }

  /** free **/
  delete [] data.map.x;
  delete [] data.map.y;
  delete [] data.map.r;
  delete [] data.map.t;
  delete [] data.hough.h;
  delete [] data.hough.rh;
  delete [] data.hough.rhi;
  
  /** free device memory **/
  hough_free();
  
  /** write output and close **/
  fout->cd();
  tout->Write();
  fout->Close();
  fin->Close();

  return 0;
}
