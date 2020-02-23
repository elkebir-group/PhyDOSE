/*
 * clonetree.h
 *
 *  Created on: 19-oct-2016
 *      Author: M. El-Kebir
 */

#ifndef CLONETREE_H
#define CLONETREE_H

#include "utils.h"
#include "basetree.h"
#include <unordered_map>
#include <unordered_set>

/// This class models a clone tree
class CloneTree : public BaseTree
{
public:
  /// Default constructor
  CloneTree();
  
  /// Copy constructor
  CloneTree(const CloneTree& other);
  
  /// Assignment operator
  CloneTree& operator =(const CloneTree& other);
  
  /// Constructor
  ///
  /// @param T Directed graph
  /// @param root Root node
  /// @param id Node identifier
  /// @param l Leaf label
  CloneTree(const Digraph& T,
            Node root,
            const StringNodeMap& id,
            const StringNodeMap& l);
  
  CloneTree(const BaseTree& tree)
    : BaseTree(tree)
    , _l(_tree)
  {
  }
    
  /// Overload constructor
  ///
  /// @param clusterTree Clone tree to expand
  /// @param expansionMap mapping of cluster labels to expanded node order
  CloneTree(const CloneTree& clusterTree,
            const std::unordered_map<std::string, StringVector>& expansionMap);
  
  /// Constructor to randomly generate expanded tree
  ///
  /// @param clusteredCloneTree Clone tree with clusters
  /// @param sep Cluster separator
  CloneTree(const CloneTree& clusteredCloneTree,
            const std::string& sep);
  
  /// Return the leaf label of the given node
  ///
  /// @param u Node
  const std::string& l(Node u) const
  {
    assert(_isLeaf[u]);
    return _l[u];
  }
  
  /// Return the set of leaf labels of subtree rooted at the given node
  ///
  /// @param u Node
  StringSet ll(Node u) const
  {
    StringSet res;
    for (Node v : _leafSubset[u])
    {
      res.insert(l(v));
    }
    return res;
  }
  
  /// Read leaf labeling
  ///
  /// @param in Input stream
  virtual bool readLeafLabeling(std::istream& in);
  
  /// Write leaf labeling
  ///
  /// @param out Output stream
  virtual void writeLeafLabeling(std::ostream& out) const;
  
  /// Write vertex labeling
  ///
  /// @param out Output stream
  /// @param lPlus Vertex labeling
  void writeVertexLabeling(std::ostream& out,
                           const StringNodeMap& lPlus) const;
  
  /// Print tree in DOT format
  ///
  /// @param out Output stream
  void writeDOT(std::ostream& out) const;
  
  /// Print tree in DOT format using the given color map
  ///
  /// @param out Output stream
  /// @param colorMap Color map
  virtual void writeDOT(std::ostream& out,
                        const StringToIntMap& colorMap) const;
  
  /// Print tree in DOT format using the given vertex labeling
  ///
  /// @param out Output stream
  /// @param lPlus Labeling of each node by a location
  void writeDOT(std::ostream& out,
                const StringNodeMap& lPlus) const;
  
  /// Print tree in DOT format using the given vertex labeling and color map
  ///
  /// @param out Output stream
  /// @param lPlus Labeling of each node by a location
  /// @param colorMap Color map
  virtual void writeDOT(std::ostream& out,
                        const StringNodeMap& lPlus,
                        const StringToIntMap& colorMap) const;
  
  /// Print tree in DOT format using the given vertex labeling and color map
  ///
  /// @param out Output stream
  /// @param lPlus Labeling of each node by a location
  /// @param colorMap Color map
  /// @param F Frequencies
  /// @param U Mixing proportions of the leaves
  void writeDOT(std::ostream& out,
                const StringNodeMap& lPlus,
                const StringToIntMap& colorMap,
                const DoubleVectorNodeMap& F,
                const DoubleNodeMap& U) const;
  
  /// Print tree in DOT format using the given color map
  ///
  /// @param out Output stream
  /// @param colorMap Color map
  /// @param U Mixing proportions of the leaves
  void writeDOT(std::ostream& out,
                const StringToIntMap& colorMap,
                const DoubleNodeMap& U) const;
  
  /// Print tree in DOT format using the given color map
  ///
  /// @param out Output stream
  /// @param lPlus Vertex labeling
  /// @param colorMap Color map
  /// @param U Mixing proportions of the leaves
  void writeDOT(std::ostream& out,
                const StringNodeMap& lPlus,
                const StringToIntMap& colorMap,
                const DoubleVectorNodeMap& U) const;
  
  /// Print tree in DOT format using the given color map
  ///
  /// @param out Output stream
  /// @param colorMap Color map
  /// @param U Mixing proportions of the leaves
  void writeDOT(std::ostream& out,
                const StringToIntMap& colorMap,
                const DoubleVectorNodeMap& U) const;

  /// Print tree in DOT format using the given vertex labeling and color map
  ///
  /// @param out Output stream
  /// @param lPlus Labeling of each node by a location
  /// @param colorMap Color map
  /// @param U Mixing proportions of the leaves
  /// @param characterLabel Character index of the character
  ///        introduced on the incoming edge of a node
  void writeDOT(std::ostream& out,
                const StringNodeMap& lPlus,
                const StringToIntMap& colorMap,
                const DoubleNodeMap& U,
                const IntNodeMap& characterLabel) const;
  
  /// Print tree in DOT format using the given vertex labeling and color map
  ///
  /// @param out Output stream
  /// @param lPlus Labeling of each node by a location
  /// @param colorMap Color map
  /// @param U Mixing proportions of the leaves
  /// @param characterLabel Character index of the character
  ///        introduced on the incoming edge of a node
  void writeDOT(std::ostream& out,
                const StringNodeMap& lPlus,
                const StringToIntMap& colorMap,
                const DoubleVectorNodeMap& U,
                const StringNodeMap& characterLabel) const;
  
  /// Print tree in DOT format using the given vertex labeling and color map
  ///
  /// @param out Output stream
  /// @param lPlus Labeling of each node by a location
  /// @param colorMap Color map
  /// @param U Mixing proportions of the leaves
  /// @param characterLabel Character index of the character
  ///        introduced on the incoming edge of a node
  void writeDOT(std::ostream& out,
                const StringNodeMap& lPlus,
                const StringToIntMap& colorMap,
                const DoubleNodeMap& U,
                const StringNodeMap& characterLabel) const;
  
  /// Return number of locations
  int getNrLocations() const
  {
    return getLocations().size();
  }
  
  /// Return location labels
  StringSet getLocations() const
  {
    StringSet res;
    
    for (NodeIt v(_tree); v != lemon::INVALID; ++v)
    {
      if (isLeaf(v))
      {
        res.insert(_l[v]);
      }
    }
    
    return res;
  }
  
  /// Return all migration edges of the provided vertex labeling
  ///
  /// @param lPlus Labeling of each node by a location
  ArcSet getMigrationEdges(const StringNodeMap& lPlus) const
  {
    ArcSet res;
    
    for (ArcIt a(_tree); a != lemon::INVALID; ++a)
    {
      Node u = _tree.source(a);
      Node v = _tree.target(a);
      if (lPlus[u] != lPlus[v])
      {
        res.insert(a);
      }
    }
    
    return res;
  }
  
  /// Merge sibling leaves labeled by the same locations
  void mergeSameSiblingLeaves();
  
  /// Return leaf labeling
  const StringNodeMap& getLeafLabeling() const
  {
    return _l;
  }
  
  /// Return split set
  SplitSet getSplits() const;
    
  /// Return vector containing all expansions of given node
  ///
  /// @param label String of cluster node label
  /// @param sep Separator for the cluster node label
  /// @param expansion Vector passed in to fill with expansions of nodes
  static void expandNode(const std::string& label,
                         const std::string& sep,
                         std::vector<StringVector>& expansion);
  
  /// Return random node expansion of given node
  ///
  /// @param label String of cluster node label
  /// @param sep Separator for the cluster node label
  /// @param expansion Vector passed in to fill with expansions of nodes
  static void randomNodeExpansion(const std::string& label,
                                  const std::string& sep,
                                  StringVector& expansion);
  
  /// Return random node expansion of given node
  ///
  /// @param label String of cluster node label
  /// @param sep Separator for the cluster node label
  /// @param expansion Vector passed in to fill with expansions of nodes
  static void levelNodeExpansion(const std::string& label,
                                 const std::string& sep,
                                 StringVector& expansion,
                                 const std::map<std::string, int> mutationToLevelMap);
  
  /// Return vector containing all expansions of given tree
  ///
  /// @param C Clone tree containing clusters to expand
  /// @param sep Separator for the cluster node label
  /// @param expansion Vector passed in to fill with expansions of tree
  static void expandTree(const CloneTree& C,
                         const std::string& sep,
                         std::vector<CloneTree>& expansion);
  
  /// Get level-order expansion of tree with respect to consensus
  static CloneTree levelExpandTree(const std::string sep, const std::map<std::string, int>& mutationToLevelMap, const CloneTree& clusterTree);
    
  /// Return Cartesian product
  /// https://www.rosettacode.org/wiki/Cartesian_product_of_two_or_more_lists#C.2B.2B
  static std::vector<std::vector<int>> product(const std::vector<std::vector<int>>& lists);
  
  /// Expand clone tree vector
  ///
  /// @param ctv Clone tree vector
  /// @param sep Separator
  /// @param expCtv Expanded clone tree vector (output)
  static void expandCloneTreeVector(const std::vector<CloneTree>& ctv,
                                    const std::string& sep,
                                    std::vector<CloneTree>& expCtv);
  
protected:
  /// Leaf labeling L(T) -> Sigma
  StringNodeMap _l;
};

typedef std::vector<CloneTree> CloneTreeVector;

std::ostream& operator<<(std::ostream& out, const CloneTree& T);

std::istream& operator>>(std::istream& in, CloneTree& T);

std::ostream& operator<<(std::ostream& out, const CloneTreeVector& trees);

std::istream& operator>>(std::istream& in, CloneTreeVector& trees);

#endif // CLONETREE_H
