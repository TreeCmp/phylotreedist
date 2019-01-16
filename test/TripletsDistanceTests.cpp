/* 
 * File:   TripletsDistanceTests.h
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

#define BOOST_TEST_MODULE tripletsDistance
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <Phyl/Tree.h>
#include <Phyl/BipartitionList.h>
#include <NumCalc/VectorTools.h>

#include <PhylotreeDist.h>
#include "TestedTreesInstances.h"
#include "NaiveAlgorithms.h"


#include <iostream>
using namespace bpp;
using namespace dist;

RootedTrees trees;

BOOST_AUTO_TEST_SUITE( Correctness )

BOOST_AUTO_TEST_CASE( Correctness )
{
	for (int i = 0; i < trees.size() - 1; i++) {
                cout << endl << "Trees " << i <<"and" << i+1 << " ";
                try {
                        BOOST_CHECK_EQUAL(NaiveAlgorithms::tripletsDistance(*trees.at(i), *trees.at(i+1)),
                                PhylotreeDist::tripletsDistance(*trees.at(i), *trees.at(i+1), true));
                } catch (Exception e) {
                        if (strcmp(e.what(), "Trees have different sets of leaves.\n") != 0) {
                                ostringstream message;
                                message << "exception thrown by PhylotreeDist::tripletsDistance(*trees.at(" << i << "), *trees.at(" << i+1 << ")): " << e.what();
                                BOOST_ERROR(message.str());                
                        }
                }
        }	
}

BOOST_AUTO_TEST_SUITE_END() //Correctness


BOOST_AUTO_TEST_SUITE( ConstraintsChecking )

BOOST_AUTO_TEST_CASE( TreesWithDifferentLeavesNumberThrowException )
{
        vector<Tree*> trees;
        Reader::getTrees ("(a,(d,e));\n"
                        "(a,(c,(d,e)));\n", trees);
	BOOST_CHECK_THROW(PhylotreeDist::tripletsDistance(*trees.at(0), *trees.at(1), true), Exception);
}
BOOST_AUTO_TEST_CASE( TreesWithDifferentLeavesSetsThrowException )
{
        vector<Tree*> trees;
        Reader::getTrees ("(a,b,c);\n"
                        "(f,b,c);\n", trees);
	BOOST_CHECK_THROW(PhylotreeDist::tripletsDistance(*trees.at(0), *trees.at(1), true), Exception);
}
BOOST_AUTO_TEST_CASE( UnrootedTreesThrowException )
{
        vector<Tree*> trees;
        Reader::getTrees ("((a,b),(c,d));\n"
                        "(a,b,(c,d));\n", trees);
	BOOST_CHECK_THROW(PhylotreeDist::tripletsDistance(*trees.at(0), *trees.at(1), true), Exception);
}
BOOST_AUTO_TEST_CASE( MultifurcatingTreesThrowException )
{
        vector<Tree*> trees;
        Reader::getTrees ("((a,b),(c,d));\n"
                        "(a,(b,c,d));\n", trees);
	BOOST_CHECK_THROW(PhylotreeDist::tripletsDistance(*trees.at(0), *trees.at(1), true), Exception);
}

BOOST_AUTO_TEST_SUITE_END() //ConstraintsChecking

/*
BOOST_AUTO_TEST_SUITE( CorrectTreeStructureRecognition )

BOOST_AUTO_TEST_CASE( TheSameTreeGives0 )
{
	Tree* treeA = treesInstances.Rab_c;
	BOOST_CHECK_EQUAL(NaiveAlgorithms::tripletsDistance(*treeA, *treeA, Phylometrics::manhattan), 0);
}
BOOST_AUTO_TEST_CASE( TheSameTreesInRevertedInputSequenceOfLeavesGives0  )
{
	Tree* treeA = treesInstances.Rab_c;
	Tree* treeB = treesInstances.Rc_ab;
	BOOST_CHECK_EQUAL(NaiveAlgorithms::tripletsDistance(*treeA, *treeB, Phylometrics::manhattan), 0);
}

BOOST_AUTO_TEST_SUITE_END() // CorrectTreeStructureRecognition


BOOST_AUTO_TEST_SUITE( CorrectMetricValues )

BOOST_AUTO_TEST_CASE( ab_c__vs__ac_b  )
{
	Tree* treeA = treesInstances.Rab_c;
	Tree* treeB = treesInstances.Rac_b;
	BOOST_CHECK_EQUAL(NaiveAlgorithms::tripletsDistance(*treeA, *treeB), 1);
}
// Sample from "The triplets Distance for Rooted Bifurcating Phylogenetic Trees" Critchlow, Pearl, Qian
BOOST_AUTO_TEST_CASE( leaves7 )
{
	Tree* treeA = treesInstances.Rbin_7A;
	Tree* treeB = treesInstances.Rbin_7B;
	BOOST_CHECK_EQUAL(NaiveAlgorithms::tripletsDistance(*treeA, *treeB), 15);
}

BOOST_AUTO_TEST_SUITE_END() // CorrectMetricValues


*/
