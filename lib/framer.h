#pragma once

#include "data.h"

namespace sipm4eic {

/*******************************************************************************/
  
class framer {

public:
  
  framer(std::vector<std::string> filenames, int frame_size = 1024) :
    _filenames(filenames), _frame_size(frame_size), _verbose(false) {
    for (const auto &filename : filenames)
      _next_spill[filename] = 0;
  };
  
  bool next_spill();
  void verbose(bool flag = true) { _verbose = flag; };
  
  typedef std::vector<data> channel_t;
  typedef std::map<int, channel_t> chip_t;
  typedef struct {
    std::map<int, chip_t> hits;
    std::vector<data> triggers;
  } device_t;
  typedef std::map<int, device_t> frame_t;
  
  std::map<int, frame_t> &frames() { return _frames; };
  std::map<int, unsigned int> &part_mask() { return _part_mask; };
  std::map<int, unsigned int> &dead_mask() { return _dead_mask; };
  
  void set_trigger_coarse_offset(int device, int offset) { _trigger_coarse_offset[device] = offset; };
  
private:
  
  bool _verbose;
  int _frame_size;
  std::vector<std::string> _filenames;
  std::map<std::string, int> _next_spill;
  std::map<int, frame_t> _frames;
  std::map<int, unsigned int> _part_mask;
  std::map<int, unsigned int> _dead_mask;

  std::map<int, int> _trigger_coarse_offset;
  
};
  
/*******************************************************************************/
  
bool framer::next_spill()
{
  bool has_data = false;
  _frames.clear();
  _part_mask.clear();
  _dead_mask.clear();
  
  /** loop over input file list **/
  for (const auto filename : _filenames) {
    
    /** open file **/
    if (_verbose) std::cout << " --- opening decoded file: " << filename << std::endl;
    if (gSystem->AccessPathName(filename.c_str())) {
      if (_verbose) std::cout << "     file does not exist: " << filename << std::endl;    continue;
    }
    auto fin = TFile::Open(filename.c_str());
    if (!fin || !fin->IsOpen()) continue;
    
    /** retrieve tree and link it **/
    auto tin = (TTree *)fin->Get("alcor");
    auto nev = tin->GetEntries();
    if (_verbose) std::cout << " --- found " << nev << " entries in tree " << std::endl;
    sipm4eic::data data;
    data.link_to_tree(tin);
    
    /** loop over events in tree **/
    for (int iev = _next_spill[filename]; iev < nev; ++iev) {
      tin->GetEntry(iev);
      
      /** start of spill **/
      if (data.is_start_spill()) {
        has_data = true;
        if (_verbose) std::cout << " --- start of spill found: event " << iev << std::endl;
	auto device = data.device;
	auto fifo = data.fifo;
	if (!_part_mask.count(device)) _part_mask[device] = 0x0;
	_part_mask[device] |= (1 << fifo);
	if (data.coarse_time_clock() == 0xdeadbeef) {
	  if (!_dead_mask.count(device)) _dead_mask[device] = 0x0;
	  _dead_mask[device] |= (1 << fifo);
	}
      }            
      
      /** ALCOR hit **/
      if (data.is_alcor_hit()) {
        auto device = data.device;
        auto chip = data.chip();
        auto channel = data.eo_channel();
        auto frame = data.coarse_time_clock() / _frame_size;
        //        if (_verbose) std::cout << " --- ALCOR hit: device=" << device << " chip=" << chip << " channel=" << channel << " frame=" << frame << std::endl;
        _frames[frame][device].hits[chip][channel].push_back(data);
      }

      /** trigger tag **/
      if (data.is_trigger_tag()) {
        auto device = data.device;
	if (_trigger_coarse_offset.count(device))
	  data.coarse -= _trigger_coarse_offset[device];
        auto frame = data.coarse_time_clock() / _frame_size;
        if (_verbose) std::cout << " --- trigger hit: device=" << device << " frame=" << frame << std::endl;
        _frames[frame][device].triggers.push_back(data);
      }
      
      /** end of spill **/
      if (data.is_end_spill()) {
        if (_verbose) std::cout << " --- end of spill found: event " << iev << std::endl;
        _next_spill[filename] = iev + 1;
        break;
      }
      
    } /** end of loop over events in tree **/
    
    fin->Close();
    
  } /** end of loop over input file list **/
  
  return has_data;
  
}
  
} /** namespace sipm4eic **/

