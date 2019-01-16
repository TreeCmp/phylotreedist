/* 
 * File:   rfTest.cpp
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
#include <vector>
#include <PhylotreeDist.h>
#include <Phyl/Newick.h>
#include <Phyl/TreeTemplate.h>
#include <unistd.h>
#include "NaiveAlgorithms.h"

using namespace std;
using namespace dist;
using namespace bpp;

void readFile(string f, vector<TreeTemplate<Node> *>& trees)
{
    cout << "\n\t\t**********" << f << endl;
    vector<TreeTemplate<Node> *> treesIn;
    Newick newickReader(false);    
    try {
        newickReader.read(f, (vector<Tree *>&) (treesIn));
    } catch (exception& e) {
        cout << "Error when reading trees. Application terminated.\n";
        return;
    }
    
    trees.resize(treesIn.size());    
    for (int i = 0; i < treesIn.size(); i++) {
            trees[i] = tools::TreesManip::createOrderedTrees(*treesIn[i]);
            delete treesIn[i];
    }
}

void compare_rf(vector<TreeTemplate<Node> *>& trees)
{
    int a, b;
    cout << trees.size() << "Trees " << endl;
    for (int i = 0; i < trees.size() - 1; i++) {
            try {
                a = NaiveAlgorithms::robinsonFoulds(*trees.at(i), *trees.at(i+1));
                b = PhylotreeDist::robinsonFoulds(*trees.at(i), *trees.at(i+1), false);
                if (a!=b) cout << endl << "trees " << i <<"and" << i+1 
                        << "  DIFFERENT res:" << a << "  " << b << endl;
            } catch (Exception e) {                
                cout << endl << "trees " << i <<"and" << i+1 
                        << e.what();
            }
    }
}

void compare_t(vector<TreeTemplate<Node> *>& trees)
{
    int a, b;
    cout << trees.size() << "Trees " << endl;
    for (int i = 0; i < trees.size() - 1; i++) {
            try {
                a = NaiveAlgorithms::tripletsDistance(*trees.at(i), *trees.at(i+1));
                b = PhylotreeDist::tripletsDistance(*trees.at(i), *trees.at(i+1), false);
                if (a!=b) cout << endl << "trees " << i <<"and" << i+1 
                        << "  DIFFERENT res:" << a << "  " << b << endl;
            } catch (Exception e) {                
                cout << endl << "trees " << i <<"and" << i+1 
                        << e.what();
            }
    }
}

void compare_mp(vector<TreeTemplate<Node> *>& trees, int* res)
{
    int b;
    cout << trees.size() << "Trees " << endl;
    for (int i = 0; i < trees.size() - 1; i++) {
            try {
                b = PhylotreeDist::perfectMatching_pairs (*trees.at(i), *trees.at(i+1), false);
                if (res[i] != b) cout << endl << "trees " << i <<"and" << i+1 
                        << "  DIFFERENT res:" << res[i] << "  " << b << endl;
            } catch (Exception e) {                
                cout << endl << "trees " << i <<"and" << i+1 
                        << e.what();
            }
    }
}



void rf()
{
    // Rooted trees return diff res (by 2)
    vector<TreeTemplate<Node> *> trees;
    readFile("../data/u5.newick", trees);
    compare_rf(trees);
        
    readFile("../data/ub17.newick", trees);
    compare_rf(trees);
    
    readFile("../data/u17-PMSplits_BogdamowiczJava.newick", trees);
    compare_rf(trees);
    
    readFile("../data/yule2_1250u_200.trees", trees);
    compare_rf(trees);    
}

void t()
{
    vector<TreeTemplate<Node> *> trees;
            
    readFile("../data/u17-PMSplits_BogdamowiczJava.newick", trees);
    compare_t(trees);
    
    readFile("../data/rb7-tripl.newick", trees);
    compare_t(trees);   
    
    readFile("../data/rb8.newick", trees);
    compare_t(trees);  
    
    readFile("../data/rb16.newick", trees);
    compare_t(trees);  
    
    readFile("../data/r55-treeBase.newick", trees);
    compare_t(trees);  
}

void n()
{
    vector<TreeTemplate<Node> *> trees;
            
    readFile("../data/u17-PMSplits_BogdamowiczJava.newick", trees);
    compare_t(trees);
    
    readFile("../data/rb7-tripl.newick", trees);
    compare_t(trees);   
    
    readFile("../data/rb8.newick", trees);
    compare_t(trees);  
    
    readFile("../data/rb16.newick", trees);
    compare_t(trees);  
    
    readFile("../data/r55-treeBase.newick", trees);
    compare_t(trees);  
}

void mp()
{
    vector<TreeTemplate<Node> *> trees;
            
    readFile("../data/rb4.newick", trees);
    int res[] = {0,8,6,4};
    compare_mp(trees, res);
}
int main()
{
    //rf();
    //t();
    //n();
    mp();
    cout << endl;
    return 0;
}