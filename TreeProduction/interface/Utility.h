#ifndef __VHTauTau_TreeMaker_Utility_h
#define __VHTauTau_TreeMaker_Utility_h

#include <string>

class TTree;

namespace vhtm {
  class Utility {
  public:
    static TTree* getTree(const std::string& tree_name="vhtree");
  };
}
#endif
