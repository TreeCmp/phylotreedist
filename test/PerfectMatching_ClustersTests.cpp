/* 
 * File:   PerfectMatching_SplitsTests.cpp
 * Author: ana
 *
 * Created on 2011-05-23, 10:58:44
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
#define BOOST_TEST_MODULE PerfectMatchingDistance
#define BOOST_TEST_DYN_LINK
#include <PhylotreeDist.h>
#include <boost/test/unit_test.hpp>
#include <Phyl/Tree.h>
#include <vector>
#include "TestedTreesInstances.h"

using namespace dist;
using namespace std;

RootedTrees trees;

BOOST_AUTO_TEST_SUITE( Correctness )
BOOST_AUTO_TEST_CASE( SameTrees )
{
	BOOST_CHECK_EQUAL(
		PhylotreeDist::perfectMatching_clusters(*trees[0], *trees[0]),
		0);
}
BOOST_AUTO_TEST_SUITE_END() //Correctness

        
BOOST_AUTO_TEST_SUITE( ConstraintsChecking )

BOOST_AUTO_TEST_CASE( ComparingTreesWithDifferentLeavesNumberThrowsException )
{
        vector<Tree*> trees;
        Reader::getTrees ("(a,(d,e));\n"
                        "((a,c),(d,e));\n", trees);
	BOOST_CHECK_THROW(PhylotreeDist::perfectMatching_splits(*trees.at(0), *trees.at(1), true), bpp::Exception);
}
BOOST_AUTO_TEST_CASE( ComparingTreesWithDifferentLeavesSetsThrowsException )
{
        vector<Tree*> trees;
        Reader::getTrees ("(a,(b,c));\n"
                        "(f,(b,c));\n", trees);
	BOOST_CHECK_THROW(PhylotreeDist::perfectMatching_splits(*trees.at(0), *trees.at(1), true), bpp::Exception);
}

BOOST_AUTO_TEST_SUITE_END() //ConstraintsChecking
