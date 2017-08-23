#include "filemanager.hpp"

namespace fs = boost::filesystem;

// Given a tile, advance it to the next tile number.
// (I.e. 1216 -> 1301 or 1316 -> 2101)
int getNextTile(int tile) {
  if (tile == 0) {
    return 0;
  }

  // Get the components of the tile number
  int i1 = tile / 1000;               // Flowcell side
  int i2 = (tile - 1000 * i1) / 100;    // Swath of the cell on this side
  int i3 = tile - 1000 * i1 - 100 * i2;   // One of 16 positions on this swath

  int overflow = 0;

  i3 += 1;
  if (i3 > 16) {
    i3 = 1;
    overflow = 1;
  }

  i2 += overflow;
  overflow = 0;
  if (i2 > 3) {
    i2 = 1;
    overflow = 1;
  }

  i1 += overflow;
  if (i1 > 2) {
    // end the iteration here
    return 0;
  }

  return i1 * 1000 + i2 * 100 + i3;
}

// Comparator function for sorting the cycles (C1.1,..,C201.1,..)
bool cyclePathCompare(const fs::path &p1, const fs::path &p2) {
  std::string s1 = p1.filename().string();
  std::string s2 = p2.filename().string();

  int i1 = std::stoi(s1.substr(1, s1.find(".") - 1));
  int i2 = std::stoi(s2.substr(1, s2.find(".") - 1));

  return i1 < i2;
}

namespace kraken {

  BCLFileManager::BCLFileManager(std::string basecalls_folder, int length, int start_cycle, int step,
                                 std::vector<int> _target_tiles, std::vector<int> _target_lanes)
      : length(length), target_tiles(_target_tiles), target_lanes(_target_lanes), target_cycle(start_cycle),
        step(step) {
    basecalls_path = fs::path(basecalls_folder);

    if (!fs::exists(basecalls_path)) {
      err(EX_NOINPUT, "%s not found.", basecalls_folder.c_str());
    }

    if (fs::is_regular_file(basecalls_path)) {
      valid = false;
      err(EX_NOINPUT, "File given. 'BaseCalls' directory expected.");
    }

    // Scan target directory and get all paths for each of the lane folders.
    this->getFilePaths();

  }

  bool BCLFileManager::is_valid() {
    return valid;
  }

  // Return the next tile we want to scan and return the cycle tiles
  // we have not read before. We return the bases from the last read tile
  // to the most recent available tile.
  TileInfo BCLFileManager::getTile() {
    if (!valid)
      return TileInfo();

    // when we have read all tiles
    if (active_tile == target_tiles.end()) {
      // advance to next lane and reset tile
      active_lane++;
      active_tile = target_tiles.begin();

      // went through all lanes, check if we are done
      // if not, reset to first lane.
      if (active_lane == target_lanes.end()) {
        if (end_reached) {
          valid = false;
          return TileInfo();
        }

        // reset lane and set next target cycle
        active_lane = target_lanes.begin();
        cycle_position = target_cycle;
        target_cycle = std::min(target_cycle + step, length - 1);
        end_reached = (target_cycle == length - 1);
      }
    }

    // wait for the target cycle for this tile and lane to appear
    while (!fs::exists(cyclePaths[*active_lane - 1][target_cycle])) {
      std::this_thread::sleep_for(std::chrono::seconds(5));
    }

    // return the [start,end] cycle indices as TileInfo, with tile and lane
    TileInfo ret = {*active_lane, *active_tile, cycle_position, target_cycle + end_reached};

    active_tile++;

    return ret;
  }

  void BCLFileManager::getFilePaths() {
    // copy those directories to the lanePaths vector that start with 'L'
    copy_if(directory_iterator(basecalls_path), directory_iterator(), back_inserter(lanePaths),
            [&](const fs::path &p) {
                return (p.filename().string()[0] == 'L');
            }
    );

    sort(lanePaths.begin(), lanePaths.end());

    // Get paths to all cycles in all lanes.
    for (fs::path &lane_path : lanePaths) {
      cyclePaths.push_back(std::vector<fs::path>());

      for (int i = 1; i <= length; ++i) {
        fs::path tmp(lane_path);
        tmp += ("/C" + std::to_string(i) + ".1");
        cyclePaths.back().push_back(tmp);
      }
    }

    // collect the lane numbers into a set
    // assumption: LXXX where XXX is the lane number
    std::unordered_map<int, int> lane_numbers;
    for (int i=0; i < lanePaths.size(); ++i) {
      int lane_num = std::stoi(lanePaths[i].filename().string().substr(1));
      lane_numbers[lane_num] = i+1;
    }

    std::vector<int> lane_selection;

    // check if target lanes exist
    for (int lane_num : target_lanes) {
      if (lane_numbers.count(lane_num) == 0) {
        copy(lanePaths.begin(), lanePaths.end(), std::ostream_iterator<fs::path>(std::cout, "\n"));
        errx(EX_NOINPUT, "Lane %i was not found. Please only select from existing lanes shown above", lane_num);
      }
      else{
        lane_selection.push_back(lane_numbers[lane_num]);
      }
    }

    // if no target lanes were specified, use all existing ones
    if (lane_selection.empty()) {
      for (int i = 1; i <= lanePaths.size(); ++i)
        lane_selection.push_back(i);
    }

    target_lanes = lane_selection;

    // if no target tiles given, iterate through all tiles
    if (target_tiles.empty()) {
      int tile = 1101;
      do {
        target_tiles.push_back(tile);
      } while (tile = getNextTile(tile));
    }


    // Output the Lane numbers we are going to analyze.
    int lanes = target_lanes.size();
    std::cout << "Using " << lanes << (lanes == 1 ? " lane directory." : " lane directories.") << "\n";
    copy(target_lanes.begin(), target_lanes.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;

    // initialize file manager state by pointing at the first lanes/tiles to analyze
    active_lane = target_lanes.begin();
    active_tile = target_tiles.begin();

/*    std::cout << "Setting up paths for temporary directories.\n";
    for (size_t i = 0; i < lanePaths.size(); ++i) {
      // store paths to the tiles in LaneX -> Tiles map
      tmpPaths.insert(std::make_pair<int, TTilePathMap>(i + 1, TTilePathMap()));

      for (int tile : target_tiles) {
        std::string s(lanePaths[i].string() + "/" + std::to_string(tile) + ".filter");
        tmpPaths[i + 1][tile] = path(s);

        // if temporary filter files exist already, delete them
        if (fs::exists(tmpPaths[i + 1][tile]))
          fs::remove(tmpPaths[i + 1][tile]);

        std::cout << tmpPaths[i+1][tile].string() << "\n";
      }
    }*/
  }
}
