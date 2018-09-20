#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <vector>


template <typename Int> Int ufind(Int x, Int* uf) {
    if (uf[x] != uf[uf[x]]) uf[x] = ufind(uf[x], uf);
    return uf[x];
} // ufind

template <typename Int> void uunion(Int x, Int y, Int* uf) {
    Int xx = ufind(x, uf);
    Int yy = ufind(y, uf);
    if (xx != yy ) if (xx < yy) uf[yy] = xx; else uf[xx] = yy;
} // uunion


int main(int argc, char* argv[]) {
    if (argc != 3) {
	std::cout << "input output\n";
	return 0;
    }

    std::vector<std::string> input;

    std::ifstream f(argv[1]);
    if (!f) return 0;
    for (std::string line; std::getline(f, line);) input.push_back(line);
    f.close();

    if (input.size() < 2) return 0;

    std::string head = input[0];
    input.erase(input.begin());

    unsigned int n = input.size();
    std::vector<std::string> alias(n);
    std::vector<unsigned int> consider(n, 0);

    for (unsigned int i = 0; i < n; ++i) {
	std::string s;
	std::stringstream ss(input[i]);
	std::getline(ss, s, '\t');
	if (s == "Description") {
	    alias[i] = "Description";
	    consider[i] = n;
	    continue;
	}
	std::getline(ss, s, '\t');
	alias[i].swap(s);
	if (alias[i] == "no_match") consider[i] = n;
	else consider[i] = std::count(alias[i].begin(), alias[i].end(), ';');
    } // for i

    std::vector<int> uf(n);
    for (unsigned int i = 0; i < n; ++i) uf[i] = i;

    const unsigned int C = 4;

    for (unsigned int i = 0; i < n - 1; ++i) {
	if (consider[i] < C) {
	    for (unsigned int j = i + 1; j < n; ++j) {
		if (consider[j] < C) {
		    bool res = false;
		    if (alias[i].size() < alias[j].size()) res = (alias[j].find(alias[i]) != std::string::npos);
		    else res = (alias[i].find(alias[j]) != std::string::npos);
		    if (res == true) uunion<int>(i, j, &uf[0]);
		}
	    } // for j
	}
    } // for i

    std::map<int, std::vector<int> > aset;
    std::map<int, std::vector<int> >::const_iterator iter;

    for (unsigned int i = 0; i < n; ++i) {
	if (i != uf[i]) aset[uf[i]].push_back(i);
    }

    std::ofstream os(argv[2]);
    os << head << "\n";
    for (unsigned int i = 0; i < n; ++i) {
	if (alias[i] == "Description") os << input[i] << "\n";
	else break;
    }

    for (iter = aset.begin(); iter != aset.end(); ++iter) {
	int pos = -1;
	std::string cur_str = alias[iter->first];
	unsigned int cur_len = cur_str.size();

	unsigned int cc = consider[iter->first];
	consider[iter->first] = n;

	for (unsigned int i = 0; i < iter->second.size(); ++i) {
	    cc += consider[iter->second[i]];
	    consider[iter->second[i]] = n;
	    unsigned int len = alias[iter->second[i]].size();
	    if (len < cur_len) {
		len = cur_len;
		cur_str = alias[iter->second[i]];
		pos = i;
	    }  else if (len == cur_len) {
		if (alias[iter->second[i]].compare(cur_str) < 0) {
		    cur_str = alias[iter->second[i]];
		    pos = i;
		}
	    }
	}

	if (cc != 0) {
	    if (pos == -1) os << input[iter->first]; else os << input[iter->second[pos]];
	    os << "\n";

	    std::cout << "Processing: "
		      << "{" << alias[iter->first] << ", " << iter->first << "} ->";
	    for (unsigned int i = 0; i < iter->second.size(); ++i) {
		std::cout << " {" << alias[iter->second[i]] << ", " << iter->second[i] << "}";
	    }
	    std::cout << "\nSelected: ";
	    if (pos == -1) std::cout << alias[iter->first];
	    else std::cout << alias[iter->second[pos]];
	    std::cout << "\n" << std::endl;
	} else {
	    os << input[iter->first] << "\n";
	    for (unsigned int i = 0; i < iter->second.size(); ++i) {
		os << input[iter->second[i]] << "\n";
	    }
	}
    } // for iter

    for (unsigned int i = 0; i < n; ++i) {
	if (consider[i] < C) os << input[i] << "\n";
    }

    os.close();

    return 0;
} // main
