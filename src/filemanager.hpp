#ifndef FILEMANAGER_HPP
#define FILEMANAGER_HPP

#include "kraken_headers.hpp"
#include <boost/filesystem.hpp>

using namespace boost::filesystem;

namespace kraken {

struct TileInfo{
	//std::vector<path> tile_paths;
	int lane_num;
	int tile_num;
	int first_cycle; // first cycle to process
	int last_cycle;  // last valid cycle available
};

class BCLFileManager {
public:
	BCLFileManager(std::string basecalls_folder, int length, int max_tile);
	TileInfo getTile();

	bool is_valid();

	typedef std::unordered_map<int, path> TTilePathMap;

	std::vector<std::vector<path> > cyclePaths;
	std::unordered_map<int, TTilePathMap> tmpPaths;

	// remember the last cycle we processed the tile in
	std::unordered_map<int, int> lastTileCycle;
	std::queue<int> active_lanes;
private:

	path basecalls_path;
	std::vector<path> lanePaths;
	int active_tile = 1101;
	int active_lane = 1;
	int length;
	int max_tile;
	int start_step = 21;
	int step = 10;
	//int step = 100;

	// set to true initially, will be made false when we don't reach the end
	bool end_reached = true;

	bool valid = true;

	void getFilePaths();

};

} // namespace

#endif
