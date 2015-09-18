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

	BCLFileManager::BCLFileManager(std::string basecalls_folder){
		basecalls_path = fs::path(basecalls_folder);
		if (!fs::exists(basecalls_path)){
			err(EX_NOINPUT, "%s not found.", basecalls_folder.c_str());
		}

		if (fs::is_regular_file(basecalls_path)){
			valid = false;
			err(EX_NOINPUT, "File given. 'BaseCalls' directory expected.");
		}

		// Scan target directory and get all paths for each of the lane folders.
		this->getFilePaths();
	}

	bool BCLFileManager::is_valid(){
		return valid;
	}

	TileInfo BCLFileManager::getTile(){
		return TileInfo();
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
			copy_if(directory_iterator(path), directory_iterator(), back_inserter(cyclePaths.back()),
					[&](const fs::path& p){
						return (fs::is_directory(p) && p.filename().string()[0] == 'C');
					}
			);

			sort(cyclePaths.back().begin(), cyclePaths.back().end(), cyclePathCompare);
			copy(cyclePaths.back().begin(), cyclePaths.back().end(), std::ostream_iterator<fs::path>(std::cout, "\n"));
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
				std::cout << tmpPaths[i+1][j].string() << "\n";
			} while( (j = getNextTile(j)) );
		}
	}
}
