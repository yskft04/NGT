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

#include	"NGT/NGTQ/LargeQuantizedGraph.h"

using namespace std;

class Loader {
public:
  Loader(const std::string &path) {
    if (path.find(".u8bin") != std::string::npos) {
      std::cerr << "type=u8bin" << std::endl;
      type = "u8";
    } else if (path.find(".i8bin") != std::string::npos) {
      std::cerr << "type=i8bin" << std::endl;
      type = "i8";
    } else if (path.find(".fbin") != std::string::npos) {
      std::cerr << "type=fbin" << std::endl;
      type = "f";
    } else {
      std::cerr << "Fatal error!!!" << std::endl;
      exit(1);
    }
    stream.open(path, std::ios::in | std::ios::binary);
    if (!stream) {
      std::cerr << "Error!" << std::endl;
      return;
    }
    stream.read(reinterpret_cast<char*>(&noOfObjects), sizeof(noOfObjects));
    stream.read(reinterpret_cast<char*>(&noOfDimensions), sizeof(noOfDimensions));
    std::cerr << "# of objects=" << noOfObjects << std::endl;
    std::cerr << "# of dimensions=" << noOfDimensions << std::endl;
    counter = 0;
  }
  std::vector<float> getObject() {
    vector<float> object;
    if (isEmpty()) {
      return object;
    }
    if (type == "u8") {
      for (uint32_t dimidx = 0; dimidx < noOfDimensions; dimidx++) {
	uint8_t v;
	stream.read(reinterpret_cast<char*>(&v), sizeof(v));
	if (stream.eof()) {
	  std::cerr << "something wrong" << std::endl;
	  return object;
	}
	object.push_back(v);
      }
    } else if (type == "i8") {
      for (uint32_t dimidx = 0; dimidx < noOfDimensions; dimidx++) {
	int8_t v;
	stream.read(reinterpret_cast<char*>(&v), sizeof(v));
	if (stream.eof()) {
	  std::cerr << "something wrong" << std::endl;
	  return object;
	}
	object.push_back(v);
      }
    } else if (type == "f") {
      for (uint32_t dimidx = 0; dimidx < noOfDimensions; dimidx++) {
	float v;
	stream.read(reinterpret_cast<char*>(&v), sizeof(v));
	if (stream.eof()) {
	  std::cerr << "something wrong" << std::endl;
	  return object;
	}
	object.push_back(v);
      }
    } else {
      std::cerr << "Fatal error!!!" << std::endl;
      exit(1);
    }
    counter++;
    return object;
  }
  bool isEmpty() {
    return noOfObjects <= counter;
  }
  std::ifstream stream;
  uint32_t	noOfObjects;
  uint32_t	noOfDimensions;
  uint32_t	counter;
  std::string	type;
};

void insert(NGT::Args &args)
{ 
  
  string	indexPath;
  string	objectPath;
  try {
    indexPath = args.get("#1");
  } catch (...) {
    std::cerr << "An index is not specified." << std::endl;
    return;
  }
  try {
    objectPath = args.get("#2");
  } catch (...) {
    std::cerr << "An object file is not specified." << std::endl;
    return;
  }

  size_t n = args.getl("n", 100);
  
  NGT::Property	property;
  property.dimension		= 128;
  property.objectType		= NGT::ObjectSpace::ObjectType::Float;
  property.distanceType	= NGT::Index::Property::DistanceType::DistanceTypeCosine;

  try {
    NGTQ::Index index(indexPath);
    Loader loader(objectPath);

    size_t cnt = 0;
    while (!loader.isEmpty()) {
      auto object = loader.getObject();
      //std::cerr << object.size() << std::endl;
      cnt++;
      if (cnt % 100000 == 0) {
	std::cerr << "loaded " << static_cast<float>(cnt) / 1000000.0 << "M objects." << std::endl;
      }
      index.insertIntoObjectRepository(object);
      if (cnt >= n) {
	break;
      }
    }
    std::cerr << "# of the objects=" << cnt << ":" << loader.noOfObjects << std::endl;

  } catch (NGT::Exception &err) {
    cerr << "Error " << err.what() << endl;
    return;
  } catch (...) {
    cerr << "Error" << endl;
    return;
  }
  return;
}


void extract(NGT::Args &args)
{ 
  
  string	objectPath;

  try {
    objectPath = args.get("#1");
  } catch (...) {
    std::cerr << "An object file is not specified." << std::endl;
    return;
  }

  size_t n = args.getl("n", 100);
  size_t dim = args.getl("d", 0);
  
  try {
    Loader loader(objectPath);

    size_t cnt = 0;
    while (!loader.isEmpty()) {
      auto object = loader.getObject();
      cnt++;
      if (cnt % 100000 == 0) {
	std::cerr << "loaded " << static_cast<float>(cnt) / 1000000.0 << "M objects." << std::endl;
      }
      if (dim != 0) {
	object.resize(dim, 0.0);
      }
      for (auto v = object.begin(); v != object.end(); ++v) {
	if (v + 1 != object.end()) {
	  std::cout << *v << "\t";
	} else {
	  std::cout << *v << std::endl;;
	}
      }
      if (n > 0 && cnt >= n) {
	break;
      }
    }
  } catch (...) {
    cerr << "Error" << endl;
    return;
  }
  return;
}

void gt(NGT::Args &args)
{ 
  string	path;

  try {
    path = args.get("#1");
  } catch (...) {
    std::cerr << "An object file is not specified." << std::endl;
    return;
  }

  std::ifstream stream;
  stream.open(path, std::ios::in | std::ios::binary);
  if (!stream) {
    std::cerr << "Error!" << std::endl;
    return;
  }
  uint32_t numQueries;
  uint32_t k;
  
  stream.read(reinterpret_cast<char*>(&numQueries), sizeof(numQueries));
  stream.read(reinterpret_cast<char*>(&k), sizeof(k));
  std::cerr << "# of queries=" << numQueries << std::endl;
  std::cerr << "k=" << k << std::endl;

  {
    std::ofstream idf;
    idf.open(path + "_gt.tsv");
    for (uint32_t qidx = 0; qidx < numQueries; qidx++) {
      for (uint32_t rank = 0; rank < k; rank++) {
	uint32_t id;
	stream.read(reinterpret_cast<char*>(&id), sizeof(id));
	idf << id;
	if (rank + 1 == k) {
	  idf << std::endl;
	} else {
	  idf << "\t";
	}
      }
    }
  }
  
  {
    std::ofstream df;
    df.open(path + "_gtdist.tsv");
    for (uint32_t qidx = 0; qidx < numQueries; qidx++) {
      for (uint32_t rank = 0; rank < k; rank++) {
	float distance;
	stream.read(reinterpret_cast<char*>(&distance), sizeof(distance));
	df << distance;
	if (rank + 1 == k) {
	  df << std::endl;
	} else {
	  df << "\t";
	}
      }
    }
  }

}


void gtRange(NGT::Args &args)
{ 
  string	path;

  try {
    path = args.get("#1");
  } catch (...) {
    std::cerr << "An object file is not specified." << std::endl;
    return;
  }

  std::ifstream stream;
  stream.open(path, std::ios::in | std::ios::binary);
  if (!stream) {
    std::cerr << "Error!" << std::endl;
    return;
  }
  uint32_t numQueries;
  uint32_t totalRes;
  uint32_t k;
  
  stream.read(reinterpret_cast<char*>(&numQueries), sizeof(numQueries));
  stream.read(reinterpret_cast<char*>(&totalRes), sizeof(totalRes));
  std::cerr << "# of queries=" << numQueries << std::endl;
  std::cerr << "totalRes=" << totalRes << std::endl;

  std::vector<int32_t> numResultsPerQuery(numQueries);
  for (size_t qidx = 0; qidx < numQueries; qidx++) {
    uint32_t v;
    stream.read(reinterpret_cast<char*>(&v), sizeof(v));
    //std::cerr << qidx << ":" << v << std::endl;
    numResultsPerQuery[qidx] = v;
  }
  {
    std::ofstream idf;
    idf.open(path + "_gt.tsv");
    std::ofstream df;
    df.open(path + "_gtdist.tsv");
    size_t count = 0;
    for (size_t qidx = 0; qidx < numQueries; qidx++) {
      if (numResultsPerQuery[qidx] == 0) {
	idf << std::endl;
	df << std::endl;
	continue;
      }
      for (size_t rank = 0; rank < numResultsPerQuery[qidx]; rank++) {
	uint32_t v;
	stream.read(reinterpret_cast<char*>(&v), sizeof(v));
	idf << v;
	df << 0.0;
	count++;
	if (rank + 1 == numResultsPerQuery[qidx]) {
	  idf << std::endl;
	  df << std::endl;
	} else {
	  idf << "\t";
	  df << "\t";
	}
      }
    }
    if (count != totalRes) {
      std::cerr << "Fatal error. " << count << ":" << totalRes << std::endl;
    }
  }
}


int
main(int argc, char **argv)
{
  NGT::Args args(argc, argv);

  string command;
  try {
    command = args.get("#0");
  } catch(...) {
    return 0;
  }
  
  try {
    if (command == "insert") {    
      insert(args);
    } else if (command == "extract") {
      extract(args);
    } else if (command == "gt") {
      gt(args);
    } else if (command == "gt-range") {
      gtRange(args);
    } else {
      std::cerr << "command error" << std::endl;
    }
  } catch(NGT::Exception &err) {
    cerr << "ngt: Error: " << err.what() << endl;
  }

  return 0;
}


