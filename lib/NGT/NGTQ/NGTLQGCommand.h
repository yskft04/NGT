//
// Copyright (C) 2021 Yahoo Japan Corporation
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#pragma once

#include	"NGT/NGTQ/LargeQuantizedGraph.h"
#include	"NGT/Command.h"
#include	"NGT/NGTQ/NGTQGCommand.h"

namespace NGTLQG {
  
  class Command : public NGTQ::Command {
  public:
    class CreateParameters : public NGTQ::Command::CreateParameters {
    public:
    };

    class SearchParameters : public NGT::Command::SearchParameters {
    public:
  };

#if defined(NGT_SHARED_MEMORY_ALLOCATOR) || defined(NGTQ_SHARED_INVERTED_INDEX)
    void create(NGT::Args &args) {};
    void build(NGT::Args &args) {};
    void expand(NGT::Args &args) {};
    void distortion(NGT::Args &args) {};
    void hierarchicalKmeans(NGT::Args &args) {};
    void quantize(NGT::Args &args) {};
    void search(NGT::Args &args) {};
    void assign(NGT::Args &args) {};
    void info(NGT::Args &args) {};
#else
    void create(NGT::Args &args);
    void build(NGT::Args &args);
    void expand(NGT::Args &args);
    void distortion(NGT::Args &args);
    void hierarchicalKmeans(NGT::Args &args);
    void quantize(NGT::Args &args);
    void search(NGT::Args &args);
    void assign(NGT::Args &args);
    void info(NGT::Args &args);
#endif
    
    void setDebugLevel(int level) { debugLevel = level; }
    int getDebugLevel() { return debugLevel; }

    void help() {
      cerr << "Usage : ngtlqg command database [data]" << endl;
      cerr << "           command : create build quantize search" << endl;
    }

    void execute(NGT::Args args) {
      string command;
      try {
	command = args.get("#0");
      } catch(...) {
	help();
	return;
      }

      debugLevel = args.getl("X", 0);

      try {
	if (debugLevel >= 1) {
	  cerr << "ngt::command=" << command << endl;
	}
	if (command == "search") {
	  search(args);
	} else if (command == "create") {
	  create(args);
	} else if (command == "build") {
	  build(args);
	} else if (command == "expand") {
	  expand(args);
	} else if (command == "distortion") {
	  distortion(args);
	} else if (command == "kmeans") {
	  hierarchicalKmeans(args);
	} else if (command == "assign") {
	  assign(args);
	} else {
	  cerr << "Illegal command. " << command << endl;
	}
      } catch(NGT::Exception &err) {
	cerr << "ngtlqg: Fatal error: " << err.what() << endl;
      }
    }

  };

}; // NGTLQG
