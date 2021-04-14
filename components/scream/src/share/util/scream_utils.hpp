#ifndef SCREAM_UTILS_HPP
#define SCREAM_UTILS_HPP

#include <ekat/ekat_assert.hpp>

#include <list>
#include <algorithm>
#include <map>

namespace scream {

// This routine tries to find an arrangment of elements that allow each
// of the input groups to be a contiguous subarray of the global arrangement.
// E.g., given the groups of elments
//    G1=(A,B,C), G2=(A,B,C,D,E), G3=(C,D), G4=(C,D,E,F), G5=((D,E,F,G).
// the ordering of elements (A,B,C,D,E,F,G) clearly has all G's as contiguous
// subset, though it might not be obvious if one scrambled the entris within
// each input group.
// If such an ordering cannot be found, we return an empty list.

// The algorithm is the following.
//  - We proceed iteratively, trying to add one group at a time to the result.
//  - Internally, we store the result as a list of lists (LoL), like
//        [ [A,B], [C], [D,E,F], [G,H], [I], ... ]
//    If the algorithm succeeds, the output of the function is a single list
//    obtainied by splicing together all the inner lists.
//  - We are only allowed to modify the outer list in the following ways:
//      a) rearrange the elements in an inner list: [D,E,F]->[F,D,E]
//      b) split an inner list into multiple lists: [[A,B,C],[D,E]]->[[A,B],[C],[D,E}]
//      c) append a new list (that does not overlap with existing sublists):
//        [[A,B],[C,D]] -> [[A,B],[C,D],[E,F,G]]
//  - The 1st group is added as a single sublist, that is, given G1,..,G5 above,
//    after the first iteration, the list of lists would be [[A,B,C]].
//  - At the generic k-th step of the algorithm, we try to add group G to the LoL,
//    doing the following:
//     - compute intersection of G with each of the sublists of LoL
//     - check that all non-empty intersections are 'contiguous'. If not, the
//       overall algorithm failed, and we can return.
//     - if there's a part of G that does not intersect with any sublist, then
//       there are two scenarios:
//        i) all intersections are at the beginning or end of LoL: good.
//        ii) intersections are neither at the beginning nor at the end of LoL.
//            In this case, the overall algorithm failed, and we can return.
//     - if G has a non-empty intersection with sublist S, there are two scenarios:
//        i) the intersection is the whole S: nothing to do.
//        ii) the intersection is smaller than S: two sub-scenarios:
//          *) S is the first or last sublist with non-empty intersection:
//             split S into [intersection][rest] or [rest][intersection]
//             (the former if it's the first intersection, the former if it's the last).
//             This would correspond to operation (a+b) on the LoL.
//          *) S is neither the first nor the last. In this case the overall
//             algorithm has failed, and we can return.
//     - if G has a portion that does not intersect with any sublist S, add it
//       as a new sublist at the end. This corresponds to operation (c) from above.

// Regardless of the order in which we process the individual lists, if there is
// an arrangement A of T's that allows to have G1,...,Gn as a contiguous sublist
// of A, the above algorithm is guaranteed to find it.

template<typename T>
std::list<T> contiguous_superset (const std::list<std::list<T>>& groups)
{
  // Avoid accessing empty lists.
  if (groups.size()==0) {
    return std::list<T>();
  }

  // Intersect two lists. Simply wraps std::set_intersection, allowing a lighter syntax.
  auto intersect = [] (const std::list<T>& lhs,
                       const std::list<T>& rhs) -> std::list<T> {
    std::list<T> out;
    std::set_intersection(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(),
                          std::back_inserter(out));
    return out;
  };

  auto difference = [] (const std::list<T>& lhs,
                       const std::list<T>& rhs) -> std::list<T> {
    std::list<T> out;
    std::set_difference(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(),
                        std::back_inserter(out));
    return out;
  };

  // Would be nice if list::sort and list::unique returned the list,
  // so one can chain them into operations. Since they don't, use a lambda.
  auto sort = [] (const std::list<T>& l) {
    auto copy = l;
    copy.sort();
    return copy;
  };
  auto unique = [] (const std::list<T>& l) {
    auto copy = l;
    copy.unique();
    return copy;
  };

  std::list<std::list<T>> lol;
  // Process the remaining groups
  for (const auto& g : groups) {
    EKAT_REQUIRE_MSG(sort(g)==g,
        "Error! Individual input lists must already be sorted.\n");
    EKAT_REQUIRE_MSG(unique(g).size()==g.size(),
        "Error! Individual input lists must not contain repeated elements.\n");

    // Keep track of where non-empty intersections happen.
    // Assuming they happen in a contiguous subset of groups (othewise the algo fails),
    // they happen in the [first_pos,last_pos) range.
    // Note: unfortunately, std::list's bidirectional iterator cannot be used as key of std::map.
    std::map<std::list<T>,std::list<T>> intersections;
    auto first = lol.end();
    auto end   = lol.end();
    bool first_found = false;
    bool end_found = false;
    auto remainder = g;
    for (auto it=lol.begin(); it!=lol.end(); ++it) {
      intersections[*it] = intersect(*it,g);
      if (intersections[*it].size()>0) {
        if (not first_found) {
          first_found = true;
          first = it;
        } else {
          // If we already foudn a "last" non-empty intersection,
          // it means we found empty intersections after the 1st one.
          // As described above, if this happens, the overall algorithm has failed.
          if (end_found) {
            return std::list<T>();
          }
        }
        remainder = difference(remainder,intersections[*it]);
      } else {
        if (first_found && not end_found) {
          end_found = true;
          end = it;
        }
      }
    }

    // If first>lol.begin and last<lol.end(), then if remainder.size()>0 the overall algo failed.
    if (first!=lol.begin() && end!=lol.end() && remainder.size()>0) {
      return std::list<T>();
    }

    auto last = std::prev(end);
    // Now, go over each non-empty intersection, and, if its size is smaller than
    // the size of the sublist of lol, we split such sublist in intersection+rest
    for (auto it = first; it!=end; ++it) {
      auto diff = difference(*it,intersections[*it]);
      if (diff.size()>0) {
        // Only the first and last non-empty intersections can be less than lol's entry.
        if (not (it==first || it==last) ) {
          return std::list<T>();
        }
            
        // The intersection is not the whole sublist of lol, so replace the sublist
        // with the intersection, and add the remainder as a sublist immediately before.
        *it = intersections[*it];
        if (it==first) {
          lol.insert(it,diff);
        } else {
          lol.insert(std::next(it),diff);
          // We added a new entry after it, but we should not iterate over that entry,
          // so simply advance it.
          std::advance(it,1);
        }
      }
    }

    // If there is a remainder, we need to add it at the front or end (depending where the intersections are)
    if (remainder.size()>0) {
      if (end_found) {
        // We found an empty intersection after non-empty ones, so intersections must be at the front.
        lol.push_front(remainder);
      } else {
        // We never found an empty intersection after non-empty ones, so intersections must be at the back.
        lol.push_back(remainder);
      }
    }
  }

  // We made it! Not splice all sublists together.
  std::list<T> out;
  for (auto& l : lol) {
    out.splice(out.end(),l);
  }
  return out;
}

} // namespace scream

#endif // SCREAM_UTILS_HPP
