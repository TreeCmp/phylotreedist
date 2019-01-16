/* 
 * File:   main.cpp
 * Author: ana
 *
 * Created on 19 maj 2011, 12:12
 */

/*
This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <iostream>
#include <streambuf>
#include <PhylotreeDist.h>
#include <Phyl/Newick.h>
#include <Phyl/Tree.h>
#include <unistd.h>

using namespace std;
using namespace dist;
using namespace bpp;


string countDistance_int(int (*metricFun)(const TreeTemplate<Node>&, const TreeTemplate<Node>&, bool, bool), TreeTemplate<Node> *t1, TreeTemplate<Node> *t2, bool constr)
{
    stringstream ss;
    try {
        ss << metricFun(*t1, *t2, false, constr);
    } catch (bpp::Exception e) { 
        ss << e.what();
    } catch (exception e) { 
        ss << e.what();
    }
    return ss.str(); 
}

string countDistance_double(double (*metricFun)(const TreeTemplate<Node>&, const TreeTemplate<Node>&, bool, bool), TreeTemplate<Node> *t1, TreeTemplate<Node> *t2, bool constr)
{
    stringstream ss;    
    try {
        ss << metricFun(*t1, *t2, false, constr);
    } catch (bpp::Exception e) { 
        ss << e.what();
    } catch (exception e) { 
        ss << e.what();
    }
    return ss.str(); 
}

void print(int i, string result, ofstream& ofs)
{
   // cout << endl << i << "\t"<< result;
    ofs << endl << i << "\t"<< result;    
}

int main(int argc, char** argv) 
{            
    /*** Getting the commandline arguments ***/ 
    stringstream message;
    int opt;
    
    string inFile = "";
    string outFile = "";
    ofstream ofs;
    
    int compareMode = 0;
    string modeName = "pairs";
    bool checkConstraints = false;
        
    int (*metricFun_int)(const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2, bool setLeavesId, bool checkNames) 
        = PhylotreeDist::robinsonFoulds;        
    double (*metricFun_double)(const TreeTemplate<Node>& trIn1, const TreeTemplate<Node>& trIn2, bool setLeavesId, bool checkNames);
    string metricName = "Robinson-Foulds";    
    int doubleRes = 0;
    
    string info = "********************************\n"
                  "         PhylotreeDist\n"
                  "********************************\n"
            "The application for comparing phylogenetic trees.\n"
            "The metrics with their trees constraints \n"
            "are listed below, under option -d description.\n"
            "Usage: \n";
    info += argv[0];
    info += " -i inputFileInNewick -o outputFile [OPTIONS]\n"
            "OPTIONS:\n"
            "-m [p|m]  comparison mode\n"
            "\tp - pair\n "
            "\tm - matrix\n"
            "-d [ms|mc|mp|rf|t|n]  distance choice\n"
            "\tms - matching: split (unrooted trees)\n"
            "\tmc - matching: clusters (rooted trees)\n"
            "\tmp - matching: pairs (rooted trees)\n"
            "\trf - Robinson-Foulds  (default if no metric is chosen)\n"
            "\trfw - Robinson-Foulds with branch weights. (unrooted, branch-weighted trees)\n"
            "\tt  - triplets (rooted binary trees)\n"
            "\tn  - nodal (manhattan metric) (unrooted trees)\n"
            "\tnw  - nodal with branch weights (manhattan metric) (unrooted trees)\n"
            "\t\tnm - manhattan metric (unrooted trees)\n"
            "\t\tnp - pythagorean metric (unrooted trees)\n"
            "\t\tnmw - manhattan metric with branch weights (unrooted, branch-weighted trees)\n"
            "\t\tnpw - pythagorean metric with branch weights (unrooted, branch-weighted trees)\n"
            "-c  check trees constraints (un/rooted, bi/multifurcating,\n"
            "    the same leaves sets) and throw exception if\n"
            "    anything is incorrect."
            "\n";
    
    while ((opt = getopt(argc, argv, "i:o:m:d:c")) != -1) {
        switch (opt) {
            case 'i':
                inFile = optarg; break;
            case 'o':
                outFile = optarg; 
                ofs.open(optarg);
                break;
            case 'm':
                if (strncmp(optarg, "p", 2) == 0) { compareMode = 0; modeName = "pairs"; }
                else if (strncmp(optarg, "m", 2) == 0) { compareMode = 1; modeName = "matrix"; }
                else {
                    cout << "Wrong compare mode (-m). The program will terminate.\n" << info; 
                    return 0;
                }
                break;
            case 'd':
                if (strncmp(optarg, "ms", 3) == 0) { metricFun_int = PhylotreeDist::perfectMatching_splits; metricName = "Matching-Splits";}
                else if (strncmp(optarg, "mc", 3) == 0) { metricFun_int = PhylotreeDist::perfectMatching_clusters; metricName = "Matching-Clusters"; }
                else if (strncmp(optarg, "mp", 3) == 0) { metricFun_int = PhylotreeDist::perfectMatching_pairs; metricName = "Matching-Pairs"; }
                else if (strncmp(optarg, "rf", 3) == 0) { metricFun_int = PhylotreeDist::robinsonFoulds; metricName = "Robinson-Foulds"; }
                else if (strncmp(optarg, "rfw", 4) == 0) { doubleRes = 1; metricFun_double = PhylotreeDist::robinsonFouldsW; metricName = "Robinson-Foulds branch weighted"; }
                else if (strncmp(optarg, "q", 3) == 0) { metricFun_int = PhylotreeDist::quartetDistance; metricName = "Quartets"; }
                else if (strncmp(optarg, "t", 3) == 0) { metricFun_int = PhylotreeDist::tripletsDistance; metricName = "Triplets"; }
                else if (strncmp(optarg, "npw", 4) == 0) { doubleRes = 1; metricFun_double = PhylotreeDist::nodalDistanceW_pythagorean; metricName = "Nodal-Pythagorean branch weighted"; }
                else if (strncmp(optarg, "np", 3) == 0) { doubleRes = 1; metricFun_double = PhylotreeDist::nodalDistance_pythagorean; metricName = "Nodal-Pythagorean"; }
                else if (strncmp(optarg, "nw", 3) == 0 || strncmp(optarg, "nmw", 4) == 0) {  doubleRes = 1; metricFun_double =  PhylotreeDist::nodalDistanceW;  metricName = "Nodal-Manhattan branch weighted"; } // default for nodal: manhattan
                else if (strncmp(optarg, "nm", 1) == 0) { metricFun_int = PhylotreeDist::nodalDistance; metricName = "Nodal-Manhattan"; } // default for nodal: manhattan
                else {
                    cout << "Wrong metric choice (-d). The program will terminate.\n" << info; 
                    return 0;
                }
                break;
            case 'c':
                checkConstraints = true;
                break;
            default:
                cout << info;                        
        }
    }
    if (inFile == "" || outFile == "") {
        cout << "ERROR\n INFO: Incorrect parameters.\n You need to give input and output filenames.\nThe program will terminate.\n\n"
                << info;
        return 0;
    }
        
    
    message << "Metric:\t\t " << metricName << endl
        << "Comparison mode: " << modeName << endl
        << "Input file:\t " << inFile << endl 
        << "Output file:\t " << outFile << endl
        << "------------" << endl;
    cout << message.str();
    ofs << message.str();
    
    /*** Reading the trees ***/ 
    cout << "Scanning input file... " << flush;
    vector<TreeTemplate<Node> *> treesIn;
    vector<TreeTemplate<Node> *> trees;
    Newick newickReader(false);    
    try {
        newickReader.read(inFile, (vector<Tree *>&) (treesIn));
    } catch (exception& e) {
        cout << "Error when reading trees. Application terminated.\n";
        return 0;
    }
    cout << treesIn.size() << " trees found." << endl;    
    
    trees.resize(treesIn.size());    
    for (int i = 0; i < treesIn.size(); i++) {
            trees[i] = tools::TreesManip::createOrderedTrees(*treesIn[i]);
            delete treesIn[i];
    }    
    
    /*** Calling distance method (depending on the mode) ***/ 
    cout << "Counting the distances: PROCESSING: "; 
    int totalTime = 0;
    totalTime = clock();
    // The Nodal metric has different signature
    if (doubleRes) {
        if (compareMode == 0) {
            cout << trees.size() -1 << "calculations";
            for (int i = 1; i < trees.size(); i++) {
               print(i, countDistance_double(metricFun_double, trees[i - 1], trees[i], checkConstraints), ofs);
            }
        } else if (compareMode == 1) {
            cout << ((trees.size() * (trees.size() -1)) / 2)  << " calculations";
            int k = 0;
            for (int i = 0; i < trees.size(); i++) {
                for (int j = i + 1; j < trees.size(); j++) {
                    print(k++, countDistance_double(metricFun_double, trees[i], trees[j], checkConstraints), ofs);
                }
            }
        }    
    // All the rest metrics have the same signature - they are caled through a delegate    
    } else {    
        if (compareMode == 0) {
            cout << trees.size() -1 << " calculations";
            for (int i = 1; i < trees.size(); i++) {
                print(i, countDistance_int(metricFun_int, trees[i - 1], trees[i], checkConstraints), ofs);
            }
        } else if (compareMode == 1) {
            cout << ((trees.size() * (trees.size() -1)) / 2)  << " calculations";
            int k = 0;
            for (int i = 0; i < trees.size(); i++) {
                for (int j = i + 1; j < trees.size(); j++) {
                    print(k++, countDistance_int(metricFun_int , trees[i], trees[j], checkConstraints), ofs);
                }
            }
        }
    }
    cout << endl
        << "Counting the distances: FINISHED." << endl << endl
        << "Total calculation time: " << (clock() - totalTime) / 1000000.0F << endl
        << "For the results see the file: " << outFile << endl;
    ofs.close();
    cout << endl;
    /*** Tiding up ***/ 
    
    for (int i = 0; i < trees.size(); i++) {
        delete trees.at(i);  
    }
    return 0;
}