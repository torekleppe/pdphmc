
//{[[BUILD_FLAGS]]}

#include "PDPHMC.hpp"
#include <algorithm>


//{[[INCLUDE_BLOC]]}

//{[[DATA_BLOC]]}

json_wrap __jw(
    //{[[JSON_FILE_NAME]]}
);

class __model{
  
public:
  __model() {}
  
  void setup(){
    if(__jw.subtreeHasMembers("data")){ //account for data=list() in R ( resulting in "data":[])
      //{[[READDATA]]} 
    }
    //{[[SETUP_BLOCK]]}
    
  }
  inline int dim() const {return(
      //{[[PARAMETER_DIM]]}
  );}
  inline int dimGenerated() const {return(
      //{[[GENERATED_DIM]]}
  );}
  
  inline void atEvent() {
    //{[[ATEVENT_BLOCK]]}
  }
  
  template <class var>
  var eval(const Eigen::Matrix<var,Eigen::Dynamic,1> &__par, 
           Eigen::Ref< Eigen::VectorXd > &__gen) const {
    
    var target = 0.0;
    
    //{[[VARIABLES_BLOCK]]}
    
    //{[[CP_FROM_PAR]]}
    
    //{[[MODEL_BLOCK]]}
    
    //{[[CP_TO_GEN]]}
    
    return(target);
  }
};


class __VARIABLES_Dimensions{
  std::vector<std::string> par_names_;
  std::vector<std::string> gen_names_;
  
  Eigen::MatrixXi par_info_;
  Eigen::MatrixXi gen_info_;
  
  std::vector<std::string> expanded_par_names_;
  std::vector<std::string> expanded_gen_names_;
  
  int dim_;
  int targetCopies_;
  
  void expandedNames(){
    /*
     *  make the complete set names for individual dimensions of the problem
     */
    for(int __i=0; __i<par_names_.size();__i++){
      if(par_info_(__i,2)==0){
        // scalar
        expanded_par_names_.push_back(par_names_[__i]);
      } else if(par_info_(__i,2)==1){
        // vector
        for(int __j=0; __j < par_info_(__i,0); __j++){
          expanded_par_names_.push_back(
            par_names_[__i] + "(" + std::to_string(__j) + ")"
          );
        }
      } else if(par_info_(__i,2)==2){
        // matrix (notice column major)
        std::string colInd;
        for(int __j=0; __j < par_info_(__i,1); __j++){ // columns
          colInd = "," + std::to_string(__j) + ")";
          for(int __k=0; __k < par_info_(__i,0); __k++){
            expanded_par_names_.push_back(
              par_names_[__i] + "(" + std::to_string(__k) + colInd
            );
          }
        }
      } else {
        std::cout << "Unknown storage type !!!" << std::endl;
      }
    }
    
    for(int __i=0; __i<gen_names_.size();__i++){
      if(gen_info_(__i,2)==0){
        // scalar
        expanded_gen_names_.push_back(gen_names_[__i]);
      } else if(gen_info_(__i,2)==1){
        // vector
        for(int __j=0; __j < gen_info_(__i,0); __j++){
          expanded_gen_names_.push_back(
            gen_names_[__i] + "(" + std::to_string(__j) + ")"
          );
        }
      } else if(gen_info_(__i,2)==2){
        // matrix (notice column major)
        std::string colInd;
        for(int __j=0; __j < gen_info_(__i,1); __j++){ // columns
          colInd = "," + std::to_string(__j) + ")";
          for(int __k=0; __k < gen_info_(__i,0); __k++){
            expanded_gen_names_.push_back(
              gen_names_[__i] + "(" + std::to_string(__k) + colInd
            );
          }
        }
      } else {
        std::cout << "Unknown storage type !!!" << std::endl;
      }
    }
  }
  
  
  void fill_dim_vec(std::string subtree,
                    Eigen::Ref<Eigen::VectorXd > q0,
                    bool checkOther) const {
    if(__jw.subtreeHasMembers(subtree)){
      Eigen::MatrixXd tmpInitM;
      Eigen::VectorXd tmpInitV;
      double tmpInitS;
      int first;
      int last;
      int len;
      int cc;
      for(int i = 0; i<par_names_.size();i++){
        first = par_info_(i,3);
        last = par_info_(i,4);
        len = last-first+1;
        if(__jw.getNumeric(subtree,par_names_[i],tmpInitS)) {
          for(int i=first;i<=last;i++) q0.coeffRef(i) = tmpInitS;
        } else if(__jw.getNumeric(subtree,par_names_[i],tmpInitV)) {
          if(tmpInitV.size() == len){
            q0.segment(first,len) = tmpInitV;
          } else {
            std::cout << "WARNING : wrong dimension of " << subtree << "::" << par_names_[i] << " (vector) , ignored" << std::endl;
          }
        } else if(__jw.getNumeric(subtree,par_names_[i],tmpInitM)) {
          if(tmpInitM.cols()*tmpInitM.rows() == len){
            cc = first;
            for(int j = 0; j<tmpInitM.cols(); j++){
              for(int k = 0; k<tmpInitM.rows(); k++){
                q0.coeffRef(cc) = tmpInitM.coeff(k,j);
                cc++;
              }
            }
            
          } else {
            std::cout << "WARNING : wrong dimension of " << subtree << "::" 
                      << par_names_[i] << " (matrix) , ignored" << std::endl;
          }
        } else {
          
          //std::cout << "no " << subtree << " found for parameter " << par_names_[i] << std::endl;
        }
      }
      if(checkOther){
        std::vector<std::string> namesInJSON;
        __jw.getAllNames(subtree,namesInJSON);
        for(int ii=0;ii<namesInJSON.size();ii++){
          //std::cout << " checking if " << namesInJSON[ii] << " is a parameter" << std::endl;
          auto its = std::find(par_names_.begin(),par_names_.end(),namesInJSON[ii]);
          if(its == par_names_.end()){
            std::cout << "WARNING : " << subtree << "::" << namesInJSON[ii] << " not a PARAMETER, ignored" << std::endl;
          }
        }
      }
    } // else{
      //std::cout << "subtree " << subtree << " empty" << std::endl;
    //}
  }
  
  
public:
  __VARIABLES_Dimensions(const int targetCopies) : targetCopies_(targetCopies) {
    
    //{[[PARAMETER_NAMES_PUSH]]}
    //{[[PARAMETER_INFO]]}
    
    //{[[GENERATED_NAMES_PUSH]]}
    //{[[GENERATED_INFO]]}
    
    /*
     * The last columns in {par,gen}_info_ give the first and last element of
     * the location of the different fields in the par and gen vectors
     */
    int __c = 0;
    for(int __i=0; __i<par_info_.rows(); __i++){
      par_info_(__i,3) = __c;
      __c += par_info_(__i,0)*par_info_(__i,1);
      par_info_(__i,4) = __c-1;
    }
    __c = 0;
    for(int __i=0; __i<gen_info_.rows(); __i++){
      gen_info_(__i,3) = __c;
      __c += gen_info_(__i,0)*gen_info_(__i,1);
      gen_info_(__i,4) = __c-1;
    }
    
    dim_ = par_info_(par_info_.rows()-1,4)+1;
    
    expandedNames();
    
  } // end constructor
  
  /*
   * Find the dimensions of par-vector corresponding to the
   * values passed using the pars-argument of the R-run method
   */
  Eigen::VectorXi getStorePars() const {
    
    
    
    int dim = par_info_(par_info_.rows()-1,4)+1;
    
    if(! __jw.subtreeHasMembers("pars")){
      return(Eigen::VectorXi::LinSpaced(dim*targetCopies_,0,dim*targetCopies_-1));
    }
    else{
      std::vector<std::string> pars;
      __jw.getStringVec("pars",pars);
      std::vector<int> tmpVec;
      int ind;
      for(int p = 0; p < pars.size(); p++){
        // try first full PARAMETER types
        //std::vector<std::string>::iterator 
        auto it  = std::find(par_names_.begin(),par_names_.end(),pars[p]);
        if(it != par_names_.end()){
          //std::cout << "found " << pars[p] << std::endl;
          ind = std::distance(par_names_.begin(), it);
          for(int parElem = par_info_(ind,3); parElem <= par_info_(ind,4); parElem++){
            tmpVec.push_back(parElem);
          }
        } else {
          // else, see if element matches a single dimension
          //std::vector<std::string>::iterator 
          auto its = std::find(expanded_par_names_.begin(),expanded_par_names_.end(),pars[p]);
          if(its != expanded_par_names_.end()){
            //std::cout << "found " << pars[p] << std::endl;
            tmpVec.push_back(std::distance(expanded_par_names_.begin(), its));
          } else {
            std::cout << "WARNING : " << pars[p] << " not a PARAMETER, ignored" << std::endl;
          }
        }
      }
      // sort and remove duplicates
      std::sort(tmpVec.begin(),tmpVec.end());
      tmpVec.erase(std::unique(tmpVec.begin(),tmpVec.end()),tmpVec.end());
      //for(int i=0;i<tmpVec.size();i++) std::cout << tmpVec[i] << "  " << expanded_par_names_[tmpVec[i]] << std::endl;
      int nstore = tmpVec.size();
      Eigen::VectorXi ret(nstore*targetCopies_);
      for(int i=0;i<nstore;i++){
        for(int cp=0;cp<targetCopies_;cp++){
          ret.coeffRef(i+cp*nstore) = tmpVec[i]+cp*dim;
        }
      }
      return(ret);
    }
  }
  
  Eigen::VectorXd getQ0(const int seed) const {
    
    rng r(seed);
    Eigen::VectorXd q0(dim_*targetCopies_);
    r.rnorm(q0.head(dim_));
    fill_dim_vec("init",q0.head(dim_),true);
    for(int cp = 1; cp < targetCopies_; cp++){
      q0.segment(cp*dim_,dim_) = q0.head(dim_);
    }
    return(q0);
  }
  
  Eigen::VectorXd getFixedMass() const {
    
    
    Eigen::VectorXd fm(dim_*targetCopies_);
    fm.setConstant(dim_*targetCopies_,-1.0);
    fill_dim_vec("fixedMass",fm.head(dim_),true);
    for(int cp = 1; cp < targetCopies_; cp++){
      fm.segment(cp*dim_,dim_) = fm.head(dim_);
    }
    
    
    return(fm);
  }
  
  std::string expanded_par_names(const int i) const { return(expanded_par_names_[i]);}
  
  std::vector<std::string> expanded_gen_names() const { return(expanded_gen_names_);}
  std::vector<std::string> expanded_par_names() const { return(expanded_par_names_);}
  
  std::vector<std::string> storeColHeaderPoint(const Eigen::VectorXi& storePars) const {
    std::vector<std::string> ret;
    if(targetCopies_==1){
    for(int i=0;i<storePars.size();i++){
      ret.push_back(expanded_par_names_[storePars(i)]);
    }
    } else {
      int origStoreN = storePars.size()/targetCopies_;
      for(int cp=0; cp < targetCopies_; cp++){
        for(int i=0; i<origStoreN;i++){
          ret.push_back(expanded_par_names_[storePars(i)]+"_cp"+std::to_string(cp));
        }
      }
    }
    if(targetCopies_==1){
    for(int i=0;i<expanded_gen_names_.size();i++){
      ret.push_back(expanded_gen_names_[i]);
    }
    } else {
      for(int cp=0; cp<targetCopies_;cp++){
        for(int i=0;i<expanded_gen_names_.size();i++){
          ret.push_back(expanded_gen_names_[i]+"_cp"+std::to_string(cp));
        }
      }
    }
    return(ret);
  } 
  
  
  std::vector<std::string> storeColHeaderInt(const Eigen::VectorXi& storePars) const {
    
    if(targetCopies_==1){
      return(expanded_gen_names_);
    } else {
      std::vector<std::string> ret;
      for(int cp=0; cp<targetCopies_;cp++){
        for(int i=0;i<expanded_gen_names_.size();i++){
          ret.push_back(expanded_gen_names_[i]+"_cp"+std::to_string(cp));
        }
      }
      return(ret);
    }
  } 
  
  void dump(){
    std::cout << "PARAMETERS : " << std::endl;
    for(int i=0;i<par_names_.size();i++ ){
      std::cout << par_names_[i] << " : " << par_info_.row(i) << std::endl;
    }
    for(int i=0;i<expanded_par_names_.size();i++){
      std::cout << i << "  " << expanded_par_names_[i] << std::endl;
    }
    if(gen_names_.size()>0){
      std::cout << "GENERATED : " << std::endl;
      for(int i=0;i<gen_names_.size();i++ ){
        std::cout << gen_names_[i] << " : " << gen_info_.row(i) << std::endl;
      }
      for(int i=0;i<expanded_gen_names_.size();i++){
        std::cout << i << "  " << expanded_gen_names_[i] << std::endl;
      }
    } else {
      std::cout << "No GENERATED variables found! " << std::endl;
    }
  }
};




int main(int argc, char** argv){
  
  int chain = 0;
  if(argc>1){
    chain = atoi(argv[1]);
    std::cout << "**** chain # " << chain << " ****" << std::endl;
  }
  
  
  
  PDPsampler< __model, {[[AD_WRAPPER_CLASS]]}, extRKN64,{[[LAMBDA_CLASS]]},{[[MASS_CLASS]]}> s_;
  
  s_.setup();
  
  
  
  // get the dimension of all variables
  __VARIABLES_Dimensions vd(s_.getTargetCopies());
  
  // get basic sampler parameters
  double tmp_control;
  
  // number of samples
  __jw.getNumeric("basic.params","samples",tmp_control);
  int samples = static_cast<int>(tmp_control);
  //std::cout << "samples : " << samples << std::endl;
  
  // Trajectory length
  double Tmax;
  __jw.getNumeric("basic.params","Tmax",Tmax);
  //std::cout << "Tmax : " << Tmax << std::endl;
  
  // fraction of trajectory time used for warmup
  double warmupFrac;
  __jw.getNumeric("basic.params","warmupFrac",warmupFrac);
  //std::cout << "warmupFrac : " << warmupFrac << std::endl;
  
  // seed used both by sampler and for the inital configuration
  int seed = 1;
  if(__jw.getNumeric("basic.params","seed",tmp_control)) seed = static_cast<int>(tmp_control);
  s_.setProperty("seed",seed+13*chain);
  
  // determine which elements of the parameters to store
  Eigen::VectorXi storePars = vd.getStorePars();
#ifdef __DEBUG__
  std::cout << " storePars : " << std::endl;
  int origDim = s_.getDim()/s_.getTargetCopies();
  for(int i=0; i< storePars.size();i++){
    std::cout << storePars[i] << "  :  " << vd.expanded_par_names(storePars[i] % origDim ) <<  std::endl;
  }
#endif
  
  // determine initial state of sampler
  Eigen::VectorXd q0 = vd.getQ0(3*seed+13*chain-1);
#ifdef __DEBUG__
  std::cout << " init : " << std::endl;
  for(int i=0; i< q0.size();i++){
    std::cout << vd.expanded_par_names(i % origDim) << " : " << q0(i) << std::endl;
  }
#endif
  
  // determine if any part of the mass matrix diagonal should be kept fixed
  Eigen::VectorXd fixedMass = vd.getFixedMass();
  Eigen::VectorXd allowsFixedMass = s_.getProperty("massallowsFixedSubvector");
  if(allowsFixedMass(0)>0.5){
    s_.setProperty("massFixedSubvector",fixedMass);
  }
  
#ifdef __DEBUG__
  std::cout << " fixedMass : " << std::endl;
  for(int i=0; i< fixedMass.size();i++){
    std::cout << vd.expanded_par_names(i % origDim) << " : " << fixedMass(i) << std::endl;
  }
#endif 
  
  // get tuning parameters for the sampler 
  std::vector<std::string> controlNames; 
  double propTmpd;
  Eigen::VectorXd propTmpv;
  __jw.getAllNames("control",controlNames);
  for(int i=0;i<controlNames.size();i++){
    if(__jw.getNumeric("control",controlNames[i],propTmpd)){
      // control field is scalar
      s_.setProperty(controlNames[i],propTmpd);
    } else if(__jw.getNumeric("control",controlNames[i],propTmpv)){
      // control field is vector
      s_.setProperty(controlNames[i],propTmpv);
    } else {
      std::cout << "WARNING : bad format on control::" << controlNames[i] << std::endl;
    }
  }
  
  s_.setPrintPrefix("chain # " + std::to_string(chain));
  
  
  // run the actual simulations
  s_.run(samples,Tmax,warmupFrac,storePars,q0);
  
  
  
  // dump samples
  
  int csvPrec = 8;
  double dtmp;
  if(__jw.getNumeric("out.files","csv.prec",dtmp)) csvPrec = static_cast<int>(dtmp);
  
  std::string fileNameBase;
  if(! __jw.getString("file.name.base",fileNameBase)){
    std::cout << "file.name.base missing in JSON file" << std::endl; 
    exit (EXIT_FAILURE);
  }
  
  
  std::vector<std::string> colHeaderPoint = vd.storeColHeaderPoint(storePars);
  std::vector<std::string> colHeaderInt = vd.storeColHeaderInt(storePars);
  
  
  // point samples
  s_.samplesToFile(csvPrec,true,colHeaderPoint,fileNameBase+"_"+std::to_string(chain)+"_point.csv");
  
  // integrated samples
  s_.samplesToFile(csvPrec,false,colHeaderInt,fileNameBase+"_"+std::to_string(chain)+"_int.csv");
  
  // diagnostics info
  s_.diagnosticsToFile(fileNameBase+"_"+std::to_string(chain)+"_diag.csv");
  
  // other info in JSON file
  jsonOut outf;
  outf.push("chain",chain);
  outf.push("CPUtime",s_.getProperty("CPUtime"));
  outf.push("lastMiDiag",s_.getProperty("lastMiDiag"));
  outf.push("fixedMi",fixedMass);
  outf.push("parNames",vd.expanded_par_names());
#ifdef __STORE_RNGS__ 
  outf.push("storeRngs_n",__rng_store_n);
  outf.push("storeRngs_u",__rng_store_u);
#endif
  outf.toFile(fileNameBase+"_"+std::to_string(chain)+"_misc.json");
  
#ifdef __STORE_EVENT_Q__  
  s_.eventSamplesToFile(fileNameBase+"_"+std::to_string(chain)+"_events.csv");
#endif
}


