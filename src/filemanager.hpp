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
	int first_tile; // first tile to process
	int last_tile;  // last valid tile available
};

class BCLFileManager {
public:
	BCLFileManager(std::string basecalls_folder, int length);
	TileInfo getTile();

	bool is_valid();

	typedef std::unordered_map<int, path> TTilePathMap;

	std::vector<std::vector<path> > cyclePaths;
	std::unordered_map<int, TTilePathMap> tmpPaths;

	// remember the last cycle we processed the tile in
	std::unordered_map<int, int> lastTileCycle;
private:

	path basecalls_path;
	std::vector<path> lanePaths;
	int active_tile = 1101;
	int active_lane = 1;
	int length;
	int step = 10;

	// set to true initially, will be made false when we don't reach the end
	bool end_reached = true;

	bool valid = true;


	void getFilePaths();

};

} // namespace

#endif
