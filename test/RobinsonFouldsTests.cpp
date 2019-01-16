/* 
 * File:   newsimpletest.cpp
 * Author: ana
 *
 * Created on 2011-05-20, 13:34:11
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

#define BOOST_TEST_MODULE RobinsonFoulds
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <Phyl/TreeTemplate.h>
#include <Phyl/BipartitionList.h>
#include <NumCalc/VectorTools.h>
#include <sstream>

#include "TestedTreesInstances.h"
#include "PhylotreeDist.h"
#include "NaiveAlgorithms.h"


using namespace bpp;
using namespace dist;

UnrootedTrees unrootedTrees;
RootedTrees rootedTrees;

BOOST_AUTO_TEST_SUITE( Correctness )

BOOST_AUTO_TEST_CASE( UnrootedTrees )
{    
	for (int i = 0; i < unrootedTrees.size() - 1; i++) {
                cout << endl << "unrootedTrees " << i <<"and" << i+1 << " ";
                try {
                        BOOST_CHECK_EQUAL(NaiveAlgorithms::robinsonFoulds(*unrootedTrees.at(i), *unrootedTrees.at(i+1)),
                                PhylotreeDist::robinsonFoulds(*unrootedTrees.at(i), *unrootedTrees.at(i+1),true));
                } catch (Exception e) {
                        if (strcmp(e.what(), "Trees have different sets of leaves.\n") != 0) {
                                ostringstream message;
                                message << "exception thrown by PhylotreeDist::robinsonFoulds(*unrootedTrees.at(" << i << "), *unrootedTrees.at(" << i+1 << ")): " << e.what();
                                BOOST_ERROR(message.str());                
                        }
                }
        }
	
}
BOOST_AUTO_TEST_CASE( RootedTrees )
{
	for (int i = 0; i < rootedTrees.size()-1; i++) {
                cout << endl << "rootedTrees " << i <<"and" << i+1 << " ";
                try {
                        BOOST_CHECK_EQUAL(NaiveAlgorithms::robinsonFoulds(*rootedTrees.at(i), *rootedTrees.at(i+1)),
                                PhylotreeDist::robinsonFoulds(*rootedTrees.at(i), *rootedTrees.at(i+1),true));
                } catch (Exception e) {
                        if (strcmp(e.what(), "Trees have different sets of leaves.\n") != 0) {
                                ostringstream message;
                                message << "exception thrown by PhylotreeDist::robinsonFoulds(*rootedTrees.at(" << i << "), *rootedTrees.at(" << i+1 << ")): " << e.what();
                                BOOST_ERROR(message.str());               
                        }              
                }
        }	
}

BOOST_AUTO_TEST_SUITE_END() //Correctness

BOOST_AUTO_TEST_SUITE( ConstraintsChecking )

BOOST_AUTO_TEST_CASE( ComparingTreesWithDifferentLeavesNumberThrowsException )
{
        vector<Tree*> trees;
        Reader::getTrees ("(a,(d,e));\n"
                        "(a,(c,(d,e)));\n", trees);
	BOOST_CHECK_THROW(PhylotreeDist::robinsonFoulds(*trees.at(0), *trees.at(1),true), bpp::Exception);
}
BOOST_AUTO_TEST_CASE( ComparingTreesWithDifferentLeavesSetsThrowsException )
{
        vector<Tree*> trees;
        Reader::getTrees ("(a,b,c);\n"
                        "(f,b,c);\n", trees);
	BOOST_CHECK_THROW(PhylotreeDist::robinsonFoulds(*trees.at(0), *trees.at(1),true), bpp::Exception);
}
BOOST_AUTO_TEST_CASE( ComparingRootedTreeWithUnroodedThrowsException )
{
        vector<Tree*> trees;
        Reader::getTrees ("((a,b),(c,d));\n"
                        "(a,b,(c,d));\n", trees);
	BOOST_CHECK_THROW(PhylotreeDist::robinsonFoulds(*trees.at(0), *trees.at(1),true), bpp::Exception);
	BOOST_CHECK_THROW(PhylotreeDist::robinsonFoulds(*trees.at(0), *trees.at(1), true), bpp::Exception);
}

BOOST_AUTO_TEST_CASE( ComparingTreebase55 )
{
        vector<Tree*> trees;
        Reader::getTreesFromFile("../data/rootedBifurcating.newick", trees);
	BOOST_CHECK_EQUAL(NaiveAlgorithms::robinsonFoulds(*rootedTrees.at(0), *rootedTrees.at(1)) + 1,
                                PhylotreeDist::robinsonFoulds(*rootedTrees.at(0), *rootedTrees.at(1),true));
}

BOOST_AUTO_TEST_SUITE_END() //ConstraintsChecking
