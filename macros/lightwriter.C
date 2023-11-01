#include "../lib/framer.h"
#include "../lib/lightio.h"

const int frame_size = 256;

std::vector<std::string> devices = {
  "kc705-192",
  "kc705-193",
  "kc705-194",
  "kc705-195",
  "kc705-196",
  "kc705-197",
  "kc705-198",
  "kc705-207"
};

void
lightwriter(std::string dirname, std::string outfilename = "lightdata.root", unsigned int max_spill = kMaxUInt)
{

  /**
   ** CREATE OUTPUT TREE
   **/

  sipm4eic::lightio io;
  io.write_to_tree(outfilename);
  
  /** 
   ** BUILD INPUT FILE LIST
   **/

  std::vector<std::string> filenames;
  for (auto device : devices) {
    for (int ififo = 0; ififo < 25; ++ififo) {
      std::string filename = dirname + "/" + device + "/decoded/alcdaq.fifo_" + std::to_string(ififo) + ".root";
      filenames.push_back(filename);
    }
  }

  /** 
   ** INITIALIZE FRAMER AND PROCESS
   **/

  std::cout << " --- initialize framer: frame size = " << frame_size << std::endl;
  sipm4eic::framer framer(filenames, frame_size);
  framer.set_trigger_coarse_offset(192, 112);
  
  /** loop over spills **/
  int n_spills = 0, n_frames = 0;
  for (int ispill = 0; ispill < max_spill && framer.next_spill(); ++ispill) {

    io.new_spill(ispill);

    for (auto &[idevice, amask] : framer.part_mask())
      io.add_part(idevice, amask);
    for (auto &[idevice, amask] : framer.dead_mask())
      io.add_dead(idevice, amask);

    /** loop over frames **/
    for (auto &[iframe, aframe] : framer.frames()) {

      io.new_frame(iframe);
      
      /** selection on Luca's trigger, device 192 **/
      if (aframe[192].triggers.size() != 1) continue;
      
      /** selection on timing scintillators, device 207 **/
      auto nsipm4 = aframe[207].hits[4].size();
      auto nsipm5 = aframe[207].hits[5].size();
      if (nsipm4 == 0 && nsipm5 == 0) continue;

      /** fill trigger0 hits **/
      auto trigger0 = aframe[192].triggers;
      for (auto &trigger : trigger0)
	io.add_trigger0(trigger.coarse_time_clock() - iframe * frame_size);

      /** fill timing hits **/
      for (auto &[ichip, achip] : aframe[207].hits) {
	for (auto &[ichannel, hits] : achip) {
	  for (auto &hit : hits) {
	    auto coarse = hit.coarse_time_clock() - iframe * frame_size;
	    io.add_timing(207, hit.device_index(), coarse, hit.fine, hit.tdc);
	  }}}
		
      /** fill cherenkov hits **/
      for (auto &[idevice, adevice] : aframe) {
	if (idevice == 207) continue;
	for (auto &[ichip, achip] : adevice.hits) {
	  for (auto &[ichannel, hits] : achip) {
	    for (auto &hit : hits) {
	      auto coarse = hit.coarse_time_clock() - iframe * frame_size;
	      io.add_cherenkov(idevice, hit.device_index(), coarse, hit.fine, hit.tdc);
	    }}}
	
      } /** end of loop over devices and hits **/

      io.add_frame();
      
    } /** end of loop over frames **/

    io.fill();

  } /** end of loop over spills **/

  /** 
   ** WRITE OUTPUT TO FILE
   **/

  io.write_and_close();

}
