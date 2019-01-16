/* 
 * File:   NodalDistanceTests.h
 * Author: ana
 *
 * Created on 20 maj 2011
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

#define BOOST_TEST_MODULE NodalDistance
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
/*
BOOST_AUTO_TEST_SUITE( ConstraintsChecking )

BOOST_AUTO_TEST_CASE( ComparingTreesWithDifferentLeavesNumberThrowsException )
{
	Tree* treeA = treesInstances.Uabc;
	Tree* treeB = treesInstances.Uabcd;
	BOOST_CHECK_THROW(Phylometrics::nodalDistance(*treeA, *treeB, Phylometrics::manhattan, true), DifferentLeavesSetsException);
}
BOOST_AUTO_TEST_CASE( ComparingTreesWithDifferentLeavesSetsThrowsException )
{
	Tree* treeA = treesInstances.Uabc;
	Tree* treeB = treesInstances.Uabd;
	BOOST_CHECK_THROW(Phylometrics::nodalDistance(*treeA, *treeB, Phylometrics::manhattan, true), DifferentLeavesSetsException);
}
BOOST_AUTO_TEST_CASE( ComparingTreesWithSameLeavesSetDoesntThrowException )
{
	Tree* treeA = treesInstances.Uabc;
	Tree* treeB = treesInstances.Rab_c;
	BOOST_CHECK_NO_THROW(Phylometrics::nodalDistance(*treeA, *treeB, Phylometrics::manhattan, true));
}

BOOST_AUTO_TEST_SUITE_END() //ConstraintsChecking


BOOST_AUTO_TEST_SUITE( CorrectTreeStructureRecognition )

BOOST_AUTO_TEST_SUITE( Manhattan_metric)

BOOST_AUTO_TEST_CASE( TheSameTreeGives0 )
{
	Tree* treeA = treesInstances.Uabcd;
	//std::cout << "Result = " << Phylometrics::nodalDistance(*treeA, *treeA, Phylometrics::manhattan, false) << endl << flush;
	BOOST_CHECK_EQUAL(Phylometrics::nodalDistance(*treeA, *treeA, Phylometrics::manhattan), 0);
}
BOOST_AUTO_TEST_CASE( TheSameTreesInRevertedInputSequenceOfLeavesGives0  )
{
	Tree* treeA = treesInstances.Rab_c;
	Tree* treeB = treesInstances.Rc_ab;
	BOOST_CHECK_EQUAL(Phylometrics::nodalDistance(*treeA, *treeB, Phylometrics::manhattan), 0);
}

BOOST_AUTO_TEST_SUITE_END() // Manhattan_metric

BOOST_AUTO_TEST_SUITE( Pythagorean_metric)

BOOST_AUTO_TEST_CASE( TheSameTreeGives0 )
{
	Tree* treeA = treesInstances.Uabcd;
	//std::cout << "Result = " << Phylometrics::nodalDistance(*treeA, *treeA, Phylometrics::pythagorean, false) << endl << flush;
	BOOST_CHECK_EQUAL(Phylometrics::nodalDistance(*treeA, *treeA, Phylometrics::manhattan), 0);
}
BOOST_AUTO_TEST_CASE( TheSameTreesInRevertedInputSequenceOfLeavesGives0  )
{
	Tree* treeA = treesInstances.Rab_c;
	Tree* treeB = treesInstances.Rc_ab;
	BOOST_CHECK_EQUAL(Phylometrics::nodalDistance(*treeA, *treeB, Phylometrics::pythagorean), 0);
}

BOOST_AUTO_TEST_SUITE_END() //Pythagorean_metric

BOOST_AUTO_TEST_SUITE_END() //CorrectTreeStructureReognition

*/
BOOST_AUTO_TEST_SUITE( Correctness )

BOOST_AUTO_TEST_SUITE( Manhattan_metric)

BOOST_AUTO_TEST_CASE( Correctness )
{
	for (int i = 0; i < trees.size() - 1; i++) {
                cout << endl << "Trees " << i <<"and" << i+1 << " ";
                try {
                        BOOST_CHECK_EQUAL(NaiveAlgorithms::nodalDistance(*trees.at(i), *trees.at(i+1), NaiveAlgorithms::manhattan, true),
                                PhylotreeDist::nodalDistance(*trees.at(i), *trees.at(i+1), PhylotreeDist::manhattan, true));
                } catch (Exception e) {
                        if (strcmp(e.what(), "Trees have different sets of leaves.\n") != 0) {
                                ostringstream message;
                                message << "exception thrown by PhylotreeDist::tripletsDistance(*trees.at(" << i << "), *trees.at(" << i+1 << ")): " << e.what();
                                BOOST_ERROR(message.str());                
                        }
                }
        }
}
        
BOOST_AUTO_TEST_SUITE_END() //Manhattan_metric

BOOST_AUTO_TEST_SUITE( Pythagorean_metric)

BOOST_AUTO_TEST_CASE( Correctness )
{
	for (int i = 0; i < trees.size() - 1; i++) {
                cout << endl << "Trees " << i <<"and" << i+1 << " ";
                try {
                        BOOST_CHECK_EQUAL(NaiveAlgorithms::nodalDistance(*trees.at(i), *trees.at(i+1), NaiveAlgorithms::pythagorean, true),
                                PhylotreeDist::nodalDistance(*trees.at(i), *trees.at(i+1), PhylotreeDist::pythagorean, true));
                } catch (Exception e) {
                        if (strcmp(e.what(), "Trees have different sets of leaves.\n") != 0) {
                                ostringstream message;
                                message << "exception thrown by PhylotreeDist::tripletsDistance(*trees.at(" << i << "), *trees.at(" << i+1 << ")): " << e.what();
                                BOOST_ERROR(message.str());                
                        }
                }
        }	
}
BOOST_AUTO_TEST_SUITE_END() // Pythagorean_metric

BOOST_AUTO_TEST_SUITE_END() //Correctness 


