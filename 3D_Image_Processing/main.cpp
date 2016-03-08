#include <iostream>
#include <string>
#include <vector>

using namespace std;

//Return paths of refurred image
vector<string> SearchImageMode(vector<string> all_str,string searchType)
{
	vector<string> refurredImagePaths;
	int size = all_str.size();
	for(int i = 0; i < size; ++i) {
		size_t pos = all_str[i].find(searchType);
		if (pos != string::npos) {
			refurredImagePaths.push_back(all_str[i]);
		}
	}
	return refurredImagePaths;
}

void main(int argc, char** argv)
{
	if (argc < 4){
		cerr << "Required both Doppler, Bmode and BmodeOrigin volume." << endl;
		system("PAUSE");
		return;
	}
	vector<string> inputPaths;
	for (int i = 0; i < argc; ++i) {
		inputPaths.push_back((string)argv[i]);
	}

	vector<string> bmodePaths = SearchImageMode(inputPaths, "bmode");
	vector<string> bmodeOriginalPaths = SearchImageMode(inputPaths, "original");
	vector<string> dopplerPaths = SearchImageMode(inputPaths, "doppler");
	
}
