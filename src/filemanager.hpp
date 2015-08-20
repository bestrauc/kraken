#ifndef FILEMANAGER_HPP
#define FILEMANAGER_HPP

#include "kraken_headers.hpp"
#include <boost/filesystem.hpp>

using namespace boost::filesystem;

namespace kraken {

struct TileInfo{
	std::vector<path> tile_paths;
	int first_tile; // first tile to process
	int last_tile;  // last valid tile available
};

class BCLFileManager {
public:
	BCLFileManager(std::string basecalls_folder);
	TileInfo getNextTile();

private:
	// Tile info by lane

};

} // namespace

#endif
