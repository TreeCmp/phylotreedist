/* 
 * File:   TestedTreesInstances.h
 * Author: ana
 *
 * Created on 20 maj 2011, 16:48
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

#ifndef TESTEDTREESINSTANCES_H
#define	TESTEDTREESINSTANCES_H
#include <Phyl/TreeTemplate.h>
#include <Phyl/Newick.h>
#include <string>
#include <vector>
#include <sstream>
using namespace bpp;
using namespace std;

class Reader
{    
public:        
        
	static void getTreesFromFile(string pathIn, vector<TreeTemplate<Node> *>& trees)
	{
            Trees treesTmp;
		Newick newickReader(false);
		try {
			newickReader.read(pathIn, trees);
                        trees = Trees::createOrdered(trees);
		 } catch (Exception e) {
			throw Exception("Error when reading trees from " + pathIn );
		 }
	}
	static void getTrees(string inTrees, vector<TreeTemplate<Node> *>& trees)
	{
		Newick n;
		istringstream istr(inTrees);
		n.read(istr, trees);
	}
};
class Trees
{
protected:
        vector<TreeTemplate<Node> *> trees;
public:
        Trees(){        }
        TreeTemplate<Node>* operator[](int i)
        {
                return trees.at(i);
        }
        TreeTemplate<Node>* at(int i)
        {
                return trees.at(i);
        }
        int size()
        {
                return trees.size();
        }
        
        static TreeTemplate<Node>& createOrdered(vector<TreeTemplate<Node> *>& treesIn)
        {
            vector<TreeTemplate<Node> *> trees;
            trees.resize(treesIn.size()); 
            for (int i = 0; i < treesIn.size(); i++) {
                trees[i] = tools::TreesManip::createOrderedTrees(*treesIn[i]);
            }
            return trees;
        }
    
        
};

class UnrootedTrees : public Trees
{
public:
        UnrootedTrees()
        {
                Reader::getTreesFromFile("data/u17-PMSplits_BogdamowiczJava.newick", trees);
        }
};

class RootedTrees : public Trees
{
public:
        RootedTrees()
        {
                Reader::getTreesFromFile("data/rb7-tripl.newick", trees);
        }
};


#endif	/* TESTEDTREESINSTANCES_H */

