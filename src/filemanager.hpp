#ifndef FILEMANAGER_HPP
#define FILEMANAGER_HPP

#include "kraken_headers.hpp"
#include <boost/filesystem.hpp>

using namespace boost::filesystem;

namespace kraken {

struct TileInfo{
	std::vector<path> tile_paths;
	//int first_tile; // first tile to process
	//int last_tile;  // last valid tile available
};

class BCLFileManager {
public:
	BCLFileManager(std::string basecalls_folder);
	TileInfo getTile();

	bool is_valid();

private:
	typedef std::unordered_map<int, path> TTilePathMap;

	path basecalls_path;

	std::vector<path> lanePaths;
	std::vector<std::vector<path> > cyclePaths;
	std::unordered_map<int, TTilePathMap> tmpPaths;

	bool valid;


	void getFilePaths();

};

} // namespace

#endif
