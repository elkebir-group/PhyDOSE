/*
 * enumerateDF.cpp
 *
 *  Created on: 20-mar-2019
 *      Author: N. Aguse
 */

#include "Rcpp.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <iostream>
#include <set>
#include <algorithm>
#include <sstream>
#include <regex>
#include <string.h>

typedef std::map<std::string,std::string> StringMap;
typedef std::map<std::string,int> StringIntMap;
typedef std::map<int,std::string> IntStringMap;
typedef std::vector<StringMap> StringMapVector;
typedef std::set<std::string> StringSet;
typedef std::vector<StringSet> StringSetVector;
typedef std::set<StringSet> StringSetSet;
typedef std::bitset<10000> Subset;
typedef std::vector<Subset> Family;
typedef std::vector<int> Feature;
typedef std::vector<Feature> FeatureFamily;
typedef std::vector<int> IntVector;
typedef std::set<int> IntSet;
typedef std::set<IntSet> IntSetSet;
typedef std::vector<IntSet> IntSetVector;

bool is_cover(Feature pi, Family& Fam, Subset Universe){
    Subset cover;
    for (int mut : pi){
        cover |= Fam[mut-1];
    }
    return cover == Universe;
}

void getDFF(Family& fam, Subset U, int m, FeatureFamily& Phi, FeatureFamily& allFeatures){
    //  FeatureFamily Phi;
    for (int i = 0; i < allFeatures.size(); ++i){
        Feature afeature = allFeatures[i];
        bool c = is_cover(afeature, fam, U);
        if (c){
            bool to_add = true;
            for (int j = 0; j < Phi.size(); ++j){
                Feature feature = Phi[j];
                if (std::includes(afeature.begin(),afeature.end(),feature.begin(),feature.end())){
                    to_add = false;
                    break;
                }
            }
            if (to_add){
                Phi.push_back(afeature);
            }
        }
        
    }
}
/// Get line from stream in a platform independent manner
std::istream& getline(std::istream& is, std::string& t);
int g_lineNumber = 0;
std::istream& getline(std::istream& is, std::string& t)
{
    ++g_lineNumber;
    
    // source: http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
    t.clear();
    
    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.
    
    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();
    
    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if(sb->sgetc() == '\n')
                sb->sbumpc();
            return is;
        case EOF:
            // Also handle the case when the last line has no line ending
            if(t.empty())
                is.setstate(std::ios::eofbit);
            return is;
        default:
            t += (char)c;
        }
    }
}

// [[Rcpp::export]]
Rcpp::List enumerateDF(std::string filename){
    // vector of maps. each entry in the vector correspond to a tree
    // map maps each node to its parent
    // for each tree, we obtain the mutation set by going from each node to the root
    // convert vector of maps to vector of intsets (mutation set)
    // convert vector of intsets to vector of vector of bitsets
    //std::cout << "enumerating..." << std::endl;
    // Step 1. Parse file to vector of int maps (i.e. vector of edges)
    int nrTrees = -1;
    
    std::string line;
    std::ifstream inFile (filename.c_str());

    if (!inFile.good()){
        Rcpp::stop("Error: could not open file");
    }
    
    getline(inFile, line);
    while (line.empty() || line[0] == '#')
    {
        getline(inFile, line);
    }
    
    std::stringstream ss(line);
    ss >> nrTrees;
    
    if (nrTrees < 0)
    {
        Rcpp::stop("Error: number of trees should be nonnegative");
    }
    StringIntMap mutToIdxMap;
    IntStringMap idxToMutMap;
    //  StringMapVector imv;
    std::vector<IntVector> edgeListsVector;
    int ct = 0;
    int m = 0; // number of mutations = number of edges + 1
    for (int i = 0; i < nrTrees; ++i)
    {
        
        int nrEdges = -1;
        getline(inFile, line);
        std::size_t found = line.find("#edges");
        while (line.empty() || found == std::string::npos)
        {
            getline(inFile, line);
            found = line.find("#edges");
        }
        ss.clear();
        ss.str(line);
        ss >> nrEdges;
        if (nrEdges < 0)
        {
            Rcpp::stop("Error: number of edges should be nonnegative");
        }
        m = nrEdges + 1;
        IntVector edgeList(m, -1);
        for (int j = 0; j < nrEdges; ++j)
        {
            ss.clear();
            ss.str("");
            getline(inFile, line);
            ss << line << std::endl;
            std::string from;
            std::string to;
            std::getline(ss, from, ' ');
            std::getline(ss, to, ' ');
            std::regex newlines_re("\n");
            auto to2 = std::regex_replace(to, newlines_re, "");
            if (mutToIdxMap.find(from) == mutToIdxMap.end() && i == 0){
                idxToMutMap[ct] = from;
                mutToIdxMap[from] = ct++;
            }
            if (mutToIdxMap.find(to2) == mutToIdxMap.end() && i == 0){
                idxToMutMap[ct] = to2;
                mutToIdxMap[to2] = ct++;
            }
            edgeList[mutToIdxMap[to2]] = mutToIdxMap[from];
        }
        edgeListsVector.push_back(edgeList);
    }
    inFile.close();
    // Step 1.5 enumerate all combinations
    FeatureFamily allFeatures;
    for (int i = 1; i <= m; ++i){
        Feature new_feature;
        for (int j = 1; j <= i; ++j){
            new_feature.push_back(j);
        }
        Feature add_feature2(new_feature);
        allFeatures.push_back(add_feature2);
        while(true){
            bool nothing_new = true;
            for (int idx = i-1; idx >= 0; --idx){
                if (new_feature[idx] < m && (idx == i-1 || new_feature[idx] < new_feature[idx+1]-1)){
                    new_feature[idx]+= 1;
                    int inc = new_feature[idx]+1;
                    // Reset all entries to the right of the incremented number
                    if(idx < i-1){
                        for(int idx2 = idx+1; idx2 < i; ++idx2){
                            new_feature[idx2] = inc++;
                        }
                    }
                    nothing_new = false;
                    break;
                }
            }
            if (nothing_new){
                break;
            }
            Feature add_feature2(new_feature);
            allFeatures.push_back(add_feature2);
        }
        
    }
    
    // Step 2. Convert to vector of mutation sets
    
    std::vector<IntSetSet> treesFeaturesInt;
    for (int i = 0; i < nrTrees; ++i){
        IntSetSet featuresSet;
        // iterate all edges in the tree to create features of the tree
        for (int j = 0; j < m; ++j){
            IntSet feature;
            int currMut = j;
            while(true){
                feature.insert(currMut);
                if (edgeListsVector[i][currMut] != -1){
                    currMut = edgeListsVector[i][currMut];
                }
                else{
                    break;
                }
            }
            featuresSet.insert(feature);
        }
        treesFeaturesInt.push_back(featuresSet);
    }
    
    // Step 3. Convert to bitsets
    std::vector<IntSetVector> idxToFeatureMap2;
    std::vector<Family> FamilyVector2;
    for (int i = 0; i < nrTrees; ++i){
        Family fam;
        IntSetVector itfMap2(m);
        int idx = 0;
        for (auto it = treesFeaturesInt[i].begin(); it != treesFeaturesInt[i].end(); ++it){
            // For every mutation set (feature) of tree i, check if other trees have that feature
            Subset feat;
            
            itfMap2[idx++] = *it;
            for (int j = 0; j < nrTrees; ++j){
                if (i == j){
                    feat[j] = 1;
                }
                else if (treesFeaturesInt[j].count(*it) == 0){
                    feat[j] = 1;
                }
            }
            fam.push_back(feat);
        }
        FamilyVector2.push_back(fam);
        idxToFeatureMap2.push_back(itfMap2);
    }
    
    // Step 4. Define the universe
    Subset Universe;
    for (int i = 0; i < nrTrees; ++i){
        Universe[i] = 1;
    }
    
    std::vector<FeatureFamily> PhiVector;
    for (auto fam : FamilyVector2){
        FeatureFamily Phi;
        getDFF(fam, Universe, m, Phi, allFeatures);
        PhiVector.push_back(Phi);
    }
    
    // Step 5. Print the DFF for each tree
    Rcpp::List ret = Rcpp::List::create();
    for (int i = 0; i < nrTrees; ++i){
        Rcpp::List treeDFF = Rcpp::List::create();
        for (auto Feature : PhiVector[i]){
            std::string featureStr = "";
            for (int mutIdx : Feature){
                for (auto mut : idxToFeatureMap2[i][mutIdx-1]){
                    featureStr += idxToMutMap[mut];
                    featureStr += " ";
                }
                featureStr += ",";
            }
            featureStr.pop_back();
            featureStr.pop_back();
            treeDFF.push_back(featureStr);
        }
        ret.push_back(treeDFF);
    }
    return ret;
}