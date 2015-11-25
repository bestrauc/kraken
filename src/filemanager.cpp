#include "filemanager.hpp"

namespace fs = boost::filesystem;

// Given a tile, advance it to the next tile number.
// (I.e. 1216 -> 1301 or 1316 -> 2101)
int getNextTile(int tile){
	// Get the components of the tile number
	int i1 = tile / 1000;               // Flowcell side
	int i2 = (tile - 1000*i1) / 100;    // Swath of the cell on this side
	int i3 = tile - 1000*i1 - 100*i2;   // One of 16 positions on this swath

	int overflow = 0;

	i3 +=1;
	if (i3 > 16){
		i3 = 1;
		overflow = 1;
	}

	i2 += overflow;
	overflow = 0;
	if (i2 > 3){
		i2 = 1;
		overflow = 1;
	}

	i1 += overflow;
	if (i1 > 2){
		return 0;
	}

	return i1*1000+i2*100+i3;
}

// Comparator function for sorting the cycles (C1.1,..,C201.1,..)
bool cyclePathCompare(const fs::path &p1, const fs::path &p2){
	std::string s1 = p1.filename().string();
	std::string s2 = p2.filename().string();

	int i1 = std::stoi(s1.substr(1, s1.find(".")-1));
	int i2 = std::stoi(s2.substr(1, s2.find(".")-1));

	return i1 < i2;
}

namespace kraken {

	BCLFileManager::BCLFileManager(std::string basecalls_folder, int length, int max_tile){
		basecalls_path = fs::path(basecalls_folder);
		this->length = length;
		this->max_tile = max_tile;

		if (!fs::exists(basecalls_path)){
			err(EX_NOINPUT, "%s not found.", basecalls_folder.c_str());
		}

		if (fs::is_regular_file(basecalls_path)){
			valid = false;
			err(EX_NOINPUT, "File given. 'BaseCalls' directory expected.");
		}

		// Scan target directory and get all paths for each of the lane folders.
		this->getFilePaths();
		for (int i=1; i <=lanePaths.size(); ++i)
			active_lanes.push(i);
	}

	bool BCLFileManager::is_valid(){
		return valid;
	}

	// Return the next tile we want to scan and return the cycle tiles
	// we have not read before. We return the bases from the last read tile
	// to the most recent available tile.
	TileInfo BCLFileManager::getTile(){
		if (!valid)
			return TileInfo();

		// when we have read all tiles
		if (active_tile == 0){
			std::cout << "END: " << end_reached << "\n";

			// if we haven't reached the end in this lane, queue it again
			if (!end_reached)
				active_lanes.push(active_lane);

			// reset tile to start and go to next lane
			active_tile = 1101;
			active_lane = active_lanes.front();
			active_lanes.pop();
			end_reached = true;

			// if no lanes are left to process, end
			if (active_lanes.empty()){
				valid = false;
				return TileInfo();
			}

		}

		// check if the next 'step' number of cycles were generated for the active_tile
		// otherwise we wait (blocks successive tiles too, but we read all together)
		int last_cycle = lastTileCycle[active_tile];
		int target_cycle = std::min(last_cycle + step, length-1);

		while (!fs::exists(cyclePaths[active_lane-1][target_cycle])){
			std::this_thread::sleep_for(std::chrono::seconds(5));
		}

		std::cout << active_lane-1 << " " << target_cycle << " " << cyclePaths[active_lane-1].size() << "\n";

		std::cout << "Target found: " << target_cycle << "\n";

		// add cycles that exist beyond our target cycle
		//while (target_cycle < length && fs::exists(cyclePaths[active_lane-1][target_cycle]))
		//	target_cycle += 1;

		//std::cout << "Target extended: " << target_cycle << "\n";

		// end_reached will be set to true if all tiles have been read to the end
		bool end_cycle = target_cycle == length-1;
		end_reached = end_reached && end_cycle;

		// workaround because seqreader checks < target_cycle. TODO: FIX
		if (end_cycle)
			target_cycle++;

		// return the [start,end] cycle indices as TileInfo, with tile and lane
		TileInfo ret = {active_lane, active_tile, last_cycle, target_cycle};

		// remember tile progress for next iteration
		lastTileCycle[active_tile] = target_cycle;

		// end earlier if we reached the max tile number in this cycle
		if (active_tile == max_tile)
			active_tile = 0;
		else
			active_tile = getNextTile(active_tile);

		return ret;
	}

	void BCLFileManager::getFilePaths(){
		copy_if(directory_iterator(basecalls_path), directory_iterator(), back_inserter(lanePaths),
				[&](const fs::path& p){
			return (p.filename().string()[0] == 'L');
		}
		);

		sort(lanePaths.begin(), lanePaths.end());

		for (fs::path &path : lanePaths){
			cyclePaths.push_back(std::vector<fs::path>());
			//copy_if(directory_iterator(path), directory_iterator(), back_inserter(cyclePaths.back()),
			//		[&](const fs::path& p){
			//			return (fs::is_directory(p) && p.filename().string()[0] == 'C');
			//		}
			//);

			for (int i=1; i <=length; ++i){
				fs::path tmp(path);
				tmp += ("/C" + std::to_string(i) + ".1");
				cyclePaths.back().push_back(tmp);
			}


			// sort cycle paths by keys C1.1, C2.1,...,Cl.1..
			//sort(cyclePaths.back().begin(), cyclePaths.back().end(), cyclePathCompare);
			// write to standard output
			//copy(cyclePaths.back().begin(), cyclePaths.back().end(), std::ostream_iterator<fs::path>(std::cout, "\n"));
		}

		std::cout << "Found " << lanePaths.size() << " lane directories.\n";
		copy(lanePaths.begin(), lanePaths.end(), std::ostream_iterator<fs::path>(std::cout, "\n"));

		std::cout << "Setting up paths for temporary directories.\n";
		for (size_t i=0; i < lanePaths.size(); ++i){
			tmpPaths.insert(std::make_pair<int, TTilePathMap>(i+1, TTilePathMap()));

			size_t j=1101;
			do{
				std::string s(lanePaths[i].string() + "/" + std::to_string(j) + ".filter");
				tmpPaths[i+1][j] = path(s);

				// if temporary filter files exist already, delete them
				if (fs::exists(tmpPaths[i+1][j]))
					fs::remove(tmpPaths[i+1][j]);

				//std::cout << tmpPaths[i+1][j].string() << "\n";
			} while( (j = getNextTile(j)) );
		}

		std::cout << "SCAN FINNISCHED \n";
	}
}
