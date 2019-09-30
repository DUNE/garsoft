// To glob.  Or not to glob?  What a stupid question.



#ifndef globbing_h
#define globbing_h



#include <glob.h>
#include <vector>
#include <string>



inline std::vector<std::string> globbing(const std::string& pat){

    glob_t glob_result;
    std::vector<std::string> ret;
	
    glob(pat.c_str(),GLOB_TILDE,NULL,&glob_result);

    for (unsigned int i=0; i<glob_result.gl_pathc; ++i){
        ret.push_back(std::string(glob_result.gl_pathv[i]));
    }

    globfree(&glob_result);
	
    return ret;
}



#endif
