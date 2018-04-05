//
//  java_interface.cpp
//  indigox
//
//  Created by Welsh, Ivan on 22/11/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//
#include <iostream>
#include <string>
#include <vector>
#ifdef BUILD_JAVA
#include <jni.h>
#endif

#include "api.hpp"

#include "utils/java_interface.hpp"
#include "utils/options.hpp"

namespace indigox {
  namespace utils {
    
    String GetEliminationOrdering(String& dgf_graph) {
#ifdef BUILD_JAVA
      static JavaVM* jvm;
      static JNIEnv* env;
      static jclass cls2 = nullptr;
      static jmethodID mid = nullptr;
      static bool ran_init = false;
      String r_str;
      jstring rv;
      String jar = Options::DATA_DIRECTORY;
      if (jar.back() != '/') jar.append("/");
      jar.append(Options::AssignElectrons::FPT::LIBTW_JAR_FILE);
      
      if (!ran_init) {
        JavaVMOption* jvm_options = new JavaVMOption[1];
        String cp = String("-Djava.class.path=") + jar;
        std::vector<char> cp_output(cp.c_str(), cp.c_str() + cp.size() + 1);
        jvm_options[0].optionString = cp_output.data();
        JavaVMInitArgs jvm_args;                    // JVM init arguments
        jvm_args.version = JNI_VERSION_1_6;
        jvm_args.nOptions = 1;
        jvm_args.options = &(*jvm_options);
        jvm_args.ignoreUnrecognized = true;
        
        jint rc = JNI_CreateJavaVM(&jvm, (void**)&env, &jvm_args);
        delete [] jvm_options;
        
        if (rc != JNI_OK) r_str = "ERROR: JVM failed to create";
        else {
          cls2 = env->FindClass("nl/uu/cs/treewidth/TDPrint");
          if (cls2 == nullptr) {
            if (env->ExceptionOccurred()) env->ExceptionDescribe();
            r_str = "ERROR: class not found!";
          } else {
            mid = env->GetStaticMethodID(cls2, "nativeMain", "(Ljava/lang/String;)Ljava/lang/String;");
          }
        }
        ran_init = true;
      }
      if (mid == nullptr) {
        if (r_str.empty()) r_str = "ERROR: method not found!";
      } else {
        rv = (jstring)env->CallStaticObjectMethod(cls2, mid, env->NewStringUTF(dgf_graph.c_str()));
        r_str = env->GetStringUTFChars(rv, 0);
      }
      
      //jvm->DestroyJavaVM();
      //std::cout << r_str << std::endl;
      return r_str;
#else
      return dgf_graph;
#endif
    }
  }
}
