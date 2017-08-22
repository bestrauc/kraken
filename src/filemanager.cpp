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

  BCLFileManager::BCLFileManager(std::string basecalls_folder, int length, std::vector<int> target_tiles) {
    basecalls_path = fs::path(basecalls_folder);
    this->length = length;
    this->target_tiles = target_tiles;

    if (!fs::exists(basecalls_path)) {
      err(EX_NOINPUT, "%s not found.", basecalls_folder.c_str());
    }

    if (fs::is_regular_file(basecalls_path)) {
      valid = false;
      err(EX_NOINPUT, "File given. 'BaseCalls' directory expected.");
    }

    // Scan target directory and get all paths for each of the lane folders.
    this->getFilePaths();
    for (int i = 1; i <= lanePaths.size(); ++i)
      active_lanes.push(i);
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
      std::cout << "END: " << end_reached << "\n";

      // if we haven't reached the end in this lane, queue it again
      if (!end_reached)
        active_lanes.push(active_lane);

      // reset tile to start and go to next lane
      active_tile = target_tiles.begin();
      active_lane = active_lanes.front();
      active_lanes.pop();
      std::cout << "Switched to lane " << active_lane << "\n";
      end_reached = true;

      // if no lanes are left to process, end
      if (active_lanes.empty()) {
        valid = false;
        return TileInfo();
      }

    }

    // check if the next 'step' number of cycles were generated for the active_tile
    // otherwise we wait (blocks successive tiles too, but we read all together)
    int next_step = step;
    int last_cycle = lastTileCycle[*active_tile];
    if (last_cycle == 0)
      next_step += start_step;

    // don't overshoot the length of the sequences
    int target_cycle = std::min(last_cycle + next_step, length - 1);

    // wait for the tiles to be written
    while (!fs::exists(cyclePaths[active_lane - 1][target_cycle])) {
      // std::cout << cyclePaths[active_lane-1][target_cycle] << "\n";
      std::this_thread::sleep_for(std::chrono::seconds(5));
    }

    // end_reached will be set to true if all tiles have been read to the end
    bool end_cycle = target_cycle == length - 1;
    end_reached = end_reached && end_cycle;

    // workaround because seqreader checks < target_cycle. TODO: FIX
    if (end_cycle)
      target_cycle++;

    // return the [start,end] cycle indices as TileInfo, with tile and lane
    TileInfo ret = {active_lane, *active_tile, last_cycle, target_cycle};

    // remember tile progress for next iteration
    lastTileCycle[*active_tile] = target_cycle;

    std::cout << "Lane: " << active_lane <<" " << *active_tile <<" " << last_cycle <<" " << target_cycle << "\n";

    active_tile++;
    return ret;
  }

  void BCLFileManager::getFilePaths() {
    copy_if(directory_iterator(basecalls_path), directory_iterator(), back_inserter(lanePaths),
            [&](const fs::path &p) {
                return (p.filename().string()[0] == 'L');
            }
    );

    std::cout << "Following target tiles: \n";
//    copy(this->target_tiles.begin(), this->target_tiles.end(), std::ostream_iterator<int>(std::cout, "\n"));

    sort(lanePaths.begin(), lanePaths.end());

    for (fs::path &path : lanePaths) {
      cyclePaths.push_back(std::vector<fs::path>());

      for (int i = 1; i <= length; ++i) {
        fs::path tmp(path);
        tmp += ("/C" + std::to_string(i) + ".1");
        std::cout << tmp << "\n";
        cyclePaths.back().push_back(tmp);
      }
    }

    int lanes = lanePaths.size();
    std::cout << "Found " << lanes << (lanes == 1 ? " lane directory." : " lane directories.") << "\n";
    copy(lanePaths.begin(), lanePaths.end(), std::ostream_iterator<fs::path>(std::cout, "\n"));

    // if no target tiles given, iterate through all tiles
    if (target_tiles.empty()) {
      int tile = 1101;
      do {
        target_tiles.push_back(tile);
      } while (tile = getNextTile(tile));
    }

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
