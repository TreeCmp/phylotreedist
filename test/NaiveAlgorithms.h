/* 
 * File:   NaiveAlgorithms.h
 * Author: ana
 *
 * Created on 21 maj 2011, 20:10
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

#ifndef NAIVEALGORITHMS_H
#define	NAIVEALGORITHMS_H

#include <Phyl/Tree.h>
#include <vector>
#include <string>
using namespace bpp;
using namespace std;

class NaiveAlgorithms
{
public:    
        static int robinsonFoulds(const Tree& tr1, const Tree& tr2)
                        throw (Exception)
        {
                int dist = 0;								// Returned value
                try {
                        dist = TreeTools::robinsonFouldsDistance(tr1, tr2, true);
                } catch (Exception) {
                        //throw Exception("Different leaves sets");
                }
                return dist;
        }
        
        
        enum nodalMetric {
                /**
                 * @brief d(T1, T2) = sqrt( add(( dist_a-b_Tr1 - dist_a-b_Tr2 ) ^ 2) )
                 */
                pythagorean,
                /**
                 * @brief d(T1, T2) = add( abs( dist_a-b_Tr1 - dist_a-b_Tr2 ) )
                 */
                manhattan
        };
        static int nodalDistance(const Tree & tr1, const Tree & tr2, nodalMetric metric, bool checkNames)
                throw (Exception)
        {
                if(!VectorTools::haveSameElements(tr1.getLeavesNames(), tr2.getLeavesNames()))
                        throw Exception("DiffetentLeavesSets [NaiveAlg]");

                double result = 0;								// Returned value
                vector<string> names = tr1.getLeavesNames();	// Names of taxons (both trees have the same leaves set which stand for taxons)
                map<string, int> leavesMap1, leavesMap2;		// Container for (name, id) pair of every leaf for tr1 and tr2.
                                                                                                                // (tr1 and tr2 have the same leaves set, but leaf's id may differ).
                int v1Id_t1, v2Id_t1, v1Id_t2, v2Id_t2;			// Ids of vertices of the same two taxons for trees (t1 or t2)
                createTaxonsMap(tr1, leavesMap1);
                createTaxonsMap(tr2, leavesMap2);
                for (unsigned int v1 = 0; v1 < names.size(); v1++) {
                        for (unsigned int v2 = v1+1; v2 < names.size(); v2++) {
                                v1Id_t1 = leavesMap1.find(names[v1])->second;
                                v2Id_t1 = leavesMap1.find(names[v2])->second;
                                v1Id_t2 = leavesMap2.find(names[v1])->second;
                                v2Id_t2 = leavesMap2.find(names[v2])->second;
                                int distT1 = (int)TreeTools::getPathBetweenAnyTwoNodes(tr1, v1Id_t1, v2Id_t1).size();
                                int distT2 = (int)TreeTools::getPathBetweenAnyTwoNodes(tr2, v1Id_t2, v2Id_t2).size();
                                
                                if (metric == NaiveAlgorithms::pythagorean) {
                                        result += pythagoreanDistance(distT1, distT2);
                                } else {
                                        result += manhattanDistance(distT1, distT2);
                                }
                        }
                }
                if (metric == NaiveAlgorithms::pythagorean) {
                        result = std::sqrt(result);
                }

                return (int)result;       
        }
        
        static int tripletsDistance(const Tree& tr1, const Tree& tr2)
                        throw (Exception)
        {
                if(!VectorTools::haveSameElements(tr1.getLeavesNames(), tr2.getLeavesNames()))
                        throw Exception("DiffetentLeavesSets [NaiveAlg]");
               
                int result = 0;									// Returned value
                vector<string> names = tr1.getLeavesNames();	// Names of taxons (both trees have the same leaves set which stand for taxons)
                int taxonsNumber = names.size();				// Number of taxons in tree
                map<string, int> taxonsTr1, taxonsTr2;			// A set of pairs <name, id_in_tree_structure> of every taxon in tree (tr1, tr2)
                createTaxonsMap(tr1, taxonsTr1);
                createTaxonsMap(tr2, taxonsTr2);

                for (int a = 0; a < taxonsNumber-2; a++) {
                          for (int b = a+1; b < taxonsNumber-1; b++)
                                  for (int c = b+1; c < taxonsNumber; c++) {
                                          //TODO zabezp przed second
                                          if (theFarestNode(tr1,
                                                          taxonsTr1.find(names[a]),
                                                          taxonsTr1.find(names[b]),
                                                          taxonsTr1.find(names[c]))
                                                        != theFarestNode(tr2,
                                                          taxonsTr2.find(names[a]),
                                                          taxonsTr2.find(names[b]),
                                                          taxonsTr2.find(names[c])))
                                                  result++;
                                  }
                        }
                 return result;
        }
        
private:

        static double pythagoreanDistance(int a, int b)
        {
                double x = (double)a;
                double y = (double)b;
                return (int)std::sqrt( std::pow(x - y, 2) );
        }

        static double manhattanDistance(int a, int b)
        {
                return std::abs(a - b);
        }


        static void createTaxonsMap(const Tree & tr, map<string, int> & leavesMap)
        {
                vector<int> ids = tr.getLeavesId();
                string name;
                for(vector<int>::iterator it = ids.begin(); it != ids.end(); it++) {
                        name = tr.getNodeName(*it);
                        leavesMap.insert(pair<string, int>(name, *it));
                }
        }
        
        static string theFarestNode(const Tree & tree, map<string, int>::const_iterator pa,
		map<string, int>::const_iterator pb, map<string, int>::const_iterator pc)
        {
                 int a = pa->second;
                 int b = pb->second;
                 int c = pc->second;

                 int ab, bc, ac;
                 vector<int> v;

                 ab = getCommonAncestorDepth(tree, a, b);
                 bc = getCommonAncestorDepth(tree, b, c);
                 ac = getCommonAncestorDepth(tree, a, c);

                 string ret;
                 if (ab >= bc && ac >= bc) ret = pa->first;
                 if (ac >= ab && bc >= ab) ret = pc->first;
                 if (bc >= ac && ab >= ac) ret = pb->first;
                 return ret;
        }
        static int getCommonAncestorDepth(const Tree & tree, int a, int b)
        {
                 vector<int> v;
                 v.push_back(a);
                 v.push_back(b);
                 int ancestor = TreeTools::getLastCommonAncestor(tree, v);
                 return TreeTools::getDepth(tree, ancestor);
        }
};
#endif	/* NAIVEALGORITHMS_H */
