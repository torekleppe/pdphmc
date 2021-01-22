
#ifndef __errorLog_HPP__
#define __errorLog_HPP__
#include <vector>
#include <string>

namespace PDPHMC_errorLog{
class errorLog{
  std::string name_;
  std::vector<std::string> theLog_;
  
public:
  errorLog(const std::string &name) : name_(name) {}
  void push(const std::string &message){theLog_.push_back(message);}
  bool isEmpty() const {return theLog_.size()==0;}
  void dump(const int maxEntries=50) const {
    if(! isEmpty()){
      std::cout << "contents of log: " << name_ << ":" << std::endl;  
      int c=1;
      for(std::vector<std::string>::const_iterator it = theLog_.begin(); it != theLog_.end(); ++it){
        std::cout << "#" << c << ": "<< *it << std::endl;
        c++;
        if(c>maxEntries) break;
      }
      if(theLog_.size()>maxEntries){
        std::cout << "further entries in log skipped" << std::endl;
      }
    }
  }
};
}
#endif




