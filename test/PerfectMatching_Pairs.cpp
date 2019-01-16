/* 
 * File:   PerfectMatching_Pairs.h
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
#define BOOST_TEST_MODULE PerfectMatching_Pairs
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <PhylotreeDist.h>
#include <Phyl/Tree.h>
#include "TestedTreesInstances.h"

using namespace dist;

RootedTrees trees;

BOOST_AUTO_TEST_SUITE( Correctness )
BOOST_AUTO_TEST_CASE( SameTree )
{
    BOOST_CHECK_EQUAL(
            PhylotreeDist::perfectMatching_pairs(*trees[35], *trees[35]), 
            0);
}
        
        
BOOST_AUTO_TEST_SUITE_END()