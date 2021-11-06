//
// Copyright (C) 2020 Yahoo Japan Corporation
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


#include	"NGT/NGTQ/NGTLQGCommand.h"

#if !defined(NGT_SHARED_MEMORY_ALLOCATOR) && !defined(NGTQ_SHARED_INVERTED_INDEX)

class PseudoRandomizedObjectList {
public:
  PseudoRandomizedObjectList(NGTQ::Quantizer::ObjectList &olist, size_t s = 10000):objectList(olist), nOfObjects(olist.size()), step(s) {}
  size_t getPseudoRandomizedID(size_t id) {
    return id;
  }
  bool multipleOpen(size_t nthreads) {
    return objectList.multipleOpen(nthreads);
  }
  bool get(size_t id, NGT::Object &data, NGT::ObjectSpace *objectSpace = 0) {
    return objectList.get(getPseudoRandomizedID(id), data, objectSpace);
  }
  bool get(const size_t streamID, size_t id, NGT::Object &data, NGT::ObjectSpace *objectSpace = 0) {
    return objectList.get(streamID, getPseudoRandomizedID(id), data, objectSpace);
  }

  NGTQ::Quantizer::ObjectList &getObjectList() { return objectList; };
  
  size_t size() { return objectList.size(); }

  NGTQ::Quantizer::ObjectList &objectList;
  const size_t nOfObjects;
  const size_t step;
};

void 
NGTLQG::Command::create(NGT::Args &args)
{
  const string usage = "Usage: ngtq create "
    " -d dimension [-o object-type (f:float|c:unsigned char)] [-D distance-function] [-n data-size] "
    "[-p #-of-thread] [-R global-codebook-range] [-r local-codebook-range] "
    "[-C global-codebook-size-limit] [-c local-codebook-size-limit] [-N local-division-no] "
    "[-T single-local-centroid (t|f)] [-e epsilon] [-i index-type (t:Tree|g:Graph)] "
    "[-M global-centroid-creation-mode (d|s)] [-L local-centroid-creation-mode (d|k|s)] "
    "[-s local-sample-coefficient] "
    "index(OUT) data.tsv(IN) rotation(IN)";

  try {
    NGTQ::Command::CreateParameters createParameters(args);

    cerr << "ngtq: Create" << endl;
#ifdef NGTQ_NEURIPS21
    std::vector<float> r;
    auto *rotation = &r;
    {
      try {
	std::string rotationPath = args.get("#3");
	cerr << "rotation is " << rotationPath << "." << endl;
	std::ifstream stream(rotationPath);
	if (!stream) {
	  std::cerr << "Cannot open the rotation. " << rotationPath << std::endl;
	  cerr << usage << endl;
	  return;
	}
	std::string line;
	while (getline(stream, line)) {
	  std::vector<std::string> tokens;
	  NGT::Common::tokenize(line, tokens, " \t");
	  for (auto &token : tokens) {
	    r.push_back(NGT::Common::strtof(token));
	  }
	}
      } catch (...) {
	rotation = 0;
      }
      std::cerr << "rotation matrix size=" << r.size() << std::endl;
    }
    
    NGTLQG::Index::create(createParameters.index, createParameters.property, createParameters.globalProperty,
			  createParameters.localProperty, rotation, createParameters.objectPath);
#else
    NGTLQG::Index::create(createParameters.index, createParameters.property, createParameters.globalProperty,
			  createParameters.localProperty);
#endif
    
#ifndef NGTQ_NEURIPS21
    if (!createParameters.objectPath.empty()) {
      cerr << "ngtlqg: Append" << endl;
      NGTQ::Index::append(createParameters.index, createParameters.objectPath, createParameters.numOfObjects);
    }
#endif
  } catch(NGT::Exception &err) {
    std::cerr << err.what() << std::endl;
    cerr << usage << endl;
  }
}

void
NGTLQG::Command::search(NGT::Args &args)
{
  
  const string usage = "Usage: ngtq search [-i g|t|s] [-n result-size] [-e epsilon] [-m mode(r|l|c|a)] "
    "[-E edge-size] [-o output-mode] [-b result expansion(begin:end:[x]step)] "
    "index(input) query.tsv(input)";
  string indexPath;
  try {
    indexPath = args.get("#1");
  } catch (...) {
    cerr << "DB is not specified" << endl;
    cerr << usage << endl;
    return;
  }

  string query;
  try {
    query = args.get("#2");
  } catch (...) {
    cerr << "Query is not specified" << endl;
    cerr << usage << endl;
    return;
  }

  size_t size		= args.getl("n", 20);
  auto exactResultExpansion	= args.getf("R", 0.0);
  char outputMode	= args.getChar("o", '-');
  float epsilon	= 0.1;

  char searchMode	= args.getChar("M", 'g');
  if (args.getString("e", "none") == "-") {
    // linear search
    epsilon = FLT_MAX;
  } else {
    epsilon = args.getf("e", 0.1);
  }
  float blobEpsilon = args.getf("B", 0.0);
  size_t edgeSize = args.getl("E", 0);
  float cutback = args.getf("C", 0.0);
  size_t explorationSize = args.getf("N", 256);

  size_t beginOfResultExpansion, endOfResultExpansion, stepOfResultExpansion;
  bool mulStep = false;
  {
    beginOfResultExpansion = stepOfResultExpansion = 1;
    endOfResultExpansion = 0;
    string str = args.getString("b", "16");
    vector<string> tokens;
    NGT::Common::tokenize(str, tokens, ":");
    if (tokens.size() >= 1) { beginOfResultExpansion = NGT::Common::strtod(tokens[0]); }
    if (tokens.size() >= 2) { endOfResultExpansion = NGT::Common::strtod(tokens[1]); }
    if (tokens.size() >= 3) { 
      if (tokens[2][0] == 'x') {
	mulStep = true;
	stepOfResultExpansion = NGT::Common::strtod(tokens[2].substr(1)); 
      } else {
	stepOfResultExpansion = NGT::Common::strtod(tokens[2]); 
      }
    }
  }
  if (debugLevel >= 1) {
    cerr << "size=" << size << endl;
    cerr << "result expansion=" << beginOfResultExpansion << "->" << endOfResultExpansion << "," << stepOfResultExpansion << endl;
  }

  NGTLQG::Index index(indexPath, true);
  std::cerr << "NGTLQGCommand::The index is open." << std::endl;
  std::cerr << "  vmsize==" << NGT::Common::getProcessVmSizeStr() << std::endl;
  std::cerr << "  peak vmsize==" << NGT::Common::getProcessVmPeakStr() << std::endl;
  auto dimension = index.getQuantizer().globalCodebookIndex.getObjectSpace().getDimension();
  try {
    ifstream		is(query);
    if (!is) {
      cerr << "Cannot open the specified file. " << query << endl;
      return;
    }
    if (outputMode == 's') { cout << "# Beginning of Evaluation" << endl; }
    string line;
    double totalTime = 0;
    int queryCount = 0;
    while(getline(is, line)) {
      NGT::Object *query;
      {
	vector<float>	queryVector;
	stringstream	linestream(line);
	while (!linestream.eof()) {
	  float value;
	  linestream >> value;
	  queryVector.push_back(value);
	}
	queryVector.resize(dimension);
	query = index.allocateObject(queryVector);
      }
      queryCount++;
      size_t resultExpansion = 0;
      for (size_t base = beginOfResultExpansion; 
	   resultExpansion <= endOfResultExpansion; 
	   base = mulStep ? base * stepOfResultExpansion : base + stepOfResultExpansion) {
	resultExpansion = base;
	NGT::ObjectDistances objects;
	NGT::Timer timer;
	timer.start();
	NGTLQG::SearchContainer searchContainer(*query);
	searchContainer.setResults(&objects);
	if (exactResultExpansion >= 1.0) {
	  searchContainer.setSize(static_cast<float>(size) * exactResultExpansion);
	  searchContainer.setExactResultSize(size);
	} else {
	  searchContainer.setSize(size);
	  searchContainer.setExactResultSize(0);
	}
	searchContainer.setEpsilon(epsilon);
	searchContainer.setBlobEpsilon(blobEpsilon);
	searchContainer.setEdgeSize(edgeSize);
	searchContainer.setCutback(cutback);
	searchContainer.setGraphExplorationSize(explorationSize);
	switch (searchMode) {
	case 'b':
	  index.searchBlobGraphNaively(searchContainer);
	  break;
	case 'g':
	  index.searchBlobGraph(searchContainer);
	  break;
	case 's':
	default:
	  index.searchBlobNaively(searchContainer);
	  break;
	}
	if (objects.size() > size) {
	  objects.resize(size);
	}
	timer.stop();
	totalTime += timer.time;
	if (outputMode == 'e') {
	  cout << "# Query No.=" << queryCount << endl;
	  cout << "# Query=" << line.substr(0, 20) + " ..." << endl;
	  cout << "# Index Type=" << "----" << endl;
	  cout << "# Size=" << size << endl;
	  cout << "# Epsilon=" << epsilon << endl;
	  cout << "# Result expansion=" << resultExpansion << endl;
	  cout << "# Distance Computation=" << index.getQuantizer().distanceComputationCount << endl;
	  cout << "# Query Time (msec)=" << timer.time * 1000.0 << endl;
	} else {
	  cout << "Query No." << queryCount << endl;
	  cout << "Rank\tIN-ID\tID\tDistance" << endl;
	}

	for (size_t i = 0; i < objects.size(); i++) {
	  cout << i + 1 << "\t" << objects[i].id << "\t";
	  cout << objects[i].distance << endl;
	}

	if (outputMode == 'e') {
	  cout << "# End of Search" << endl;
	} else {
	  cout << "Query Time= " << timer.time << " (sec), " << timer.time * 1000.0 << " (msec)" << endl;
	}
      } 
      if (outputMode == 'e') {
	cout << "# End of Query" << endl;
      }
      index.deleteObject(query);
    } 
    if (outputMode == 'e') {
      cout << "# Average Query Time (msec)=" << totalTime * 1000.0 / (double)queryCount << endl;
      cout << "# Number of queries=" << queryCount << endl;
      cout << "# End of Evaluation" << endl;
    } else {
      cout << "Average Query Time= " << totalTime / (double)queryCount  << " (sec), " 
	   << totalTime * 1000.0 / (double)queryCount << " (msec), (" 
	   << totalTime << "/" << queryCount << ")" << endl;
    }
  } catch (NGT::Exception &err) {
    cerr << "Error " << err.what() << endl;
    cerr << usage << endl;
  } catch (...) {
    cerr << "Error" << endl;
    cerr << usage << endl;
  }
  std::cerr << "NGTLQGCommand::The end of search" << std::endl;
  std::cerr << "  vmsize==" << NGT::Common::getProcessVmSizeStr() << std::endl;
  std::cerr << "  peak vmsize==" << NGT::Common::getProcessVmPeakStr() << std::endl;
  index.close();
}

void 
NGTLQG::Command::build(NGT::Args &args)
{
  const std::string usage = "Usage: ngtlqg build  [-Q dimension-of-subvector] [-E max-number-of-edges] index";
  string indexPath;
  try {
    indexPath = args.get("#1");
  } catch (...) {
    cerr << "An index is not specified." << endl;
    cerr << usage << endl;
    return;
  }
  char mode = args.getChar("m", '-');

  if (mode == 'g') {
    NGTLQG::Index::buildNGTLQG(indexPath);
    return;
  }

  std::vector<std::vector<float>> quantizerCodebook;

  {
    try {
      const std::string codebookPath = args.get("#2");
      cerr << "codebook is " << codebookPath << "." << endl;
      std::ifstream stream(codebookPath);
      if (!stream) {
	std::cerr << "Cannot open the codebook. " << codebookPath << std::endl;
	cerr << usage << endl;
	return;
      }
      std::string line;
      while (getline(stream, line)) {
	std::vector<std::string> tokens;
	NGT::Common::tokenize(line, tokens, " \t");
	std::vector<float> object;
	for (auto &token : tokens) {
	  object.push_back(NGT::Common::strtof(token));
	}
	if (!quantizerCodebook.empty() && quantizerCodebook[0].size() != object.size()) {
	  cerr << "The specified quantizer codebook is invalid. " << quantizerCodebook[0].size()
	       << ":" << object.size() << ":" << quantizerCodebook.size() << ":" << line << endl;
	  cerr << usage << endl;
	  return;
	}
	if (!object.empty()) {
	  quantizerCodebook.push_back(object);
	}
      }
    } catch (...) {}
  }

  std::vector<uint32_t> codebookIndex;
  {
    try {
      std::string codebookIndexPath = args.get("#3");
      cerr << "codebook index is " << codebookIndexPath << "." << endl;
      std::ifstream stream(codebookIndexPath);
      if (!stream) {
	std::cerr << "Cannot open the codebook index. " << codebookIndexPath << std::endl;
	cerr << usage << endl;
	return;
      }
      std::string line;
      while (getline(stream, line)) {
	std::vector<std::string> tokens;
	NGT::Common::tokenize(line, tokens, " \t");
	std::vector<float> object;
	if (tokens.size() != 1) {
	  cerr << "The specified codebook index is invalid. " << line << std::endl;
	  cerr << usage << endl;
	  return;
	}
	codebookIndex.push_back(NGT::Common::strtol(tokens[0]));
      }

    } catch (...) {}
  }

  std::vector<uint32_t> objectIndex;
  {
    try {
      std::string objectIndexPath = args.get("#4");
      cerr << "object index is " << objectIndexPath << "." << endl;
      std::ifstream stream(objectIndexPath);
      if (!stream) {
	std::cerr << "Cannot open the codebook index. " << objectIndexPath << std::endl;
	cerr << usage << endl;
	return;
      }
      std::string line;
      while (getline(stream, line)) {
	std::vector<std::string> tokens;
	NGT::Common::tokenize(line, tokens, " \t");
	std::vector<float> object;
	if (tokens.size() != 1) {
	  cerr << "The specified codebook index is invalid. " << line << std::endl;
	  cerr << usage << endl;
	  return;
	}
	objectIndex.push_back(NGT::Common::strtol(tokens[0]));
      }

    } catch (...) {}
  }

  size_t beginID = args.getl("s", 1);
  size_t size = args.getl("n", 0);
  size_t endID = beginID + size - 1;

  std::cerr << "quantizer codebook size=" << quantizerCodebook.size() << std::endl;
  std::cerr << "codebook index size=" << codebookIndex.size() << std::endl;
  std::cerr << "object index size=" << objectIndex.size() << std::endl;

  if (mode == 'q') {
    NGTLQG::Index::buildNGTQ(indexPath, quantizerCodebook, codebookIndex, objectIndex, beginID, endID);
    return;
  }
  NGTLQG::Index::build(indexPath, quantizerCodebook, codebookIndex, objectIndex, beginID, endID);
}

void
NGTLQG::Command::expand(NGT::Args &args)
{
  const std::string usage = "Usage: ngtlqg prebuild";
  string indexPath;
  try {
    indexPath = args.get("#1");
  } catch (...) {
    cerr << "An index is not specified." << endl;
    cerr << usage << endl;
    return;
  }

  NGTLQG::Index index(indexPath, true);
  NGT::Index globalIndex(indexPath + "/global");
  
  auto &quantizer = static_cast<NGTQ::QuantizerInstance<uint8_t>&>(index.getQuantizer());
  auto &invertedIndex = quantizer.invertedIndex;
  auto &globalGraph = static_cast<NGT::GraphAndTreeIndex&>(globalIndex.getIndex());

  std::cerr << "invertedIndex size=" << invertedIndex.size() - 1 << std::endl;

  size_t numOfObjects = 0;
  for (size_t idx = 1; idx < invertedIndex.size(); idx++) {
    numOfObjects += invertedIndex.at(idx)->size();
  }
  std::cerr << "# of objects=" << numOfObjects << std::endl;

  vector<uint32_t> iviid(numOfObjects + 1);
  for (size_t idx = 1; idx < invertedIndex.size(); idx++) {
    auto &entry = *invertedIndex.at(idx);
    for (size_t i = 0; i < entry.size(); i++) {
      iviid[entry[i].id] = idx;
      if (entry[i].id == 0) {
	std::cerr << "Fatal error" << std::endl;
	abort();
      }
    }
  }

  size_t addedCount = 0;
  size_t notAddedCount = 0;
  size_t max = 10;
  size_t k = 10;
  float epsilon = 0.0;
  NGT::Object query(&quantizer.globalCodebookIndex.getObjectSpace());
  for (size_t no = 0; no < max; no++) {
    for (size_t idx = 1; idx < invertedIndex.size(); idx++) {
      auto &node = *globalGraph.NGT::GraphIndex::NeighborhoodGraph::getNode(idx);
      auto &entry = *invertedIndex.at(idx);
      if (no < entry.size()) {
	auto id = entry[no].id;
	quantizer.objectList.get(id, query, &quantizer.globalCodebookIndex.getObjectSpace());
	NGT::ObjectDistances objects;
	NGTLQG::SearchContainer searchContainer(query);
	searchContainer.setSize(k);
	searchContainer.setEpsilon(epsilon);
	searchContainer.setResults(&objects);
	index.searchBlobGraphNaively(searchContainer);
	for (size_t i = 0; i < objects.size(); i++) {
	  auto dstidx = iviid[objects[i].id];
	  if (dstidx != idx) {
	    auto s = node.size();
	    globalGraph.addEdge(node, dstidx, objects[i].distance, false);
	    if (s != node.size()) {
	      addedCount++;
	    } else {
	      notAddedCount++;
	    }
	    auto &dstnode = *globalGraph.NGT::GraphIndex::NeighborhoodGraph::getNode(dstidx);
	    s = dstnode.size();
	    globalGraph.addEdge(dstnode, idx, objects[i].distance, false);
	    if (s != dstnode.size()) {
	      //std::cerr << "added!!" << std::endl;
	      addedCount++;
	    } else {
	      //std::cerr << "not added!!" << std::endl;
	      notAddedCount++;
	    }
	  }
	}
      }
    }
    std::cerr << "added=" << addedCount << " not added=" << notAddedCount << std::endl;
    globalIndex.save();
  }
  
}


void
NGTLQG::Command::distortion(NGT::Args &args)
{
  const std::string usage = "Usage: ngtlqg distortion";
  string indexPath;
  try {
    indexPath = args.get("#1");
  } catch (...) {
    cerr << "An index is not specified." << endl;
    cerr << usage << endl;
    return;
  }

  NGTLQG::Index index(indexPath);
  auto &quantizer = static_cast<NGTQ::QuantizerInstance<uint8_t>&>(index.getQuantizer());
  auto &invertedIndex = quantizer.invertedIndex;
  std::cerr << "invertedIndex size=" << invertedIndex.size() << std::endl;

  std::vector<std::vector<size_t>> subspaces;
  for (size_t idx = 1; idx < invertedIndex.size(); idx++) {
    auto &entry = *invertedIndex.at(idx);
    if (entry.subspaceID >= subspaces.size()) {
      subspaces.resize(entry.subspaceID + 1);
    }
    subspaces[entry.subspaceID].push_back(idx);
  }

  vector<pair<double, size_t>> distortions(subspaces.size(), make_pair(0.0, 0));
  std::vector<size_t> subspaceIDs;
  std::cerr << distortions.size() << ":" << distortions[0].first << ":" << distortions[0].second << std::endl;
  for (size_t ssidx = 0; ssidx < subspaces.size(); ssidx++) {
    for (auto iviid : subspaces[ssidx]) {
      auto &entry = *invertedIndex.at(iviid);
      std::cerr << ssidx << " entry size=" << entry.size() << " subspace ID=" << entry.subspaceID << std::endl;
      for (size_t i = 0; i < entry.size(); i++) {
	if (subspaceIDs.size() <= entry[i].id) {
	  subspaceIDs.resize(entry[i].id + 1);
	}
	subspaceIDs[entry[i].id] = ssidx;
	NGT::Object object(&quantizer.globalCodebookIndex.getObjectSpace());
	quantizer.objectList.get(entry[i].id, object, &quantizer.globalCodebookIndex.getObjectSpace());
	std::vector<float> v = quantizer.globalCodebookIndex.getObjectSpace().getObject(object);
	float subspaceObject[quantizer.globalCodebookIndex.getObjectSpace().getPaddedDimension()];
	(*quantizer.generateResidualObject)(object, // object
					    ssidx, // subspace ID
					    subspaceObject); // subspace objects
	double distortion = 0.0;
	for (size_t j = 0; j < v.size(); j++) {
	  auto d = static_cast<double>(subspaceObject[j]);
	  distortion += d * d;
	}
	distortion = sqrt(distortion);
	distortions[ssidx].first += distortion;
	distortions[ssidx].second++;
      }
    }
  }

  
  double distortion = 0.0;
  size_t count = 0;
  for (size_t ssidx = 0; ssidx < distortions.size(); ssidx++) {
    distortion += distortions[ssidx].first;
    count += distortions[ssidx].second;
  }
  distortion /= count;
  std::cerr << "# of subspaces=" << subspaces.size() << std::endl;
  std::cerr << "distortion=" << distortion << std::endl;
}

class HKNode {
public:
  bool leaf;
};

class HKLeafNode : public HKNode {
public:
  HKLeafNode():id(0){ leaf = true; }
  std::vector<uint32_t> members;
  uint32_t id;
};

class HKInternalNode : public HKNode {
public:
  HKInternalNode() { leaf = false; }
  std::vector<std::pair<uint32_t, std::vector<float>>> children;
};


int32_t hierarchcalKmeansSearch(std::vector<HKNode*> &nodes, int32_t rootID, float *object) {
  auto nodeID = rootID;
  while (true) {
    auto *node = nodes[nodeID];
    if (node->leaf) {
      return nodeID;
    } else {
      HKInternalNode &internalNode = static_cast<HKInternalNode&>(*node);      
      float min = std::numeric_limits<float>::max();
      int32_t minid = 0;
      for (auto &c : internalNode.children) {
	auto d = NGT::PrimitiveComparator::compareL2(reinterpret_cast<float*>(&object[0]),
						     c.second.data(), c.second.size());
	if (d < min) {
	  min = d;
	  minid = c.first;
	}
      }
      nodeID = minid;
    }
  }
  return -1;
}

void aggregateObjects(HKLeafNode &leafNode, std::vector<std::vector<float>> &vectors,
		      NGT::ObjectSpace &objectSpace, PseudoRandomizedObjectList &objectList)
{
  vectors.reserve(leafNode.members.size() + 1);
  NGT::Object obj(&objectSpace);
  for (auto &m : leafNode.members) {
    objectList.get(m, obj, &objectSpace);
    vectors.push_back(objectSpace.getObject(obj));
  }
}

void aggregateObjects(HKLeafNode &leafNode, std::vector<std::vector<float>> &vectors,
		      NGT::ObjectSpace &objectSpace, PseudoRandomizedObjectList &objectList,
		      std::vector<float> &object)
{
  aggregateObjects(leafNode, vectors, objectSpace, objectList);
  vectors.push_back(std::move(object));
}

void hierarchicalKmeansSplit(uint32_t id, std::vector<std::vector<float>> &vectors,
			     std::vector<HKNode*> &nodes, int32_t leafNodeID, NGT::Clustering &clustering)
{
  HKLeafNode &leafNode = static_cast<HKLeafNode&>(*nodes[leafNodeID]);
  std::vector<NGT::Clustering::Cluster> clusters;
  clustering.kmeans(vectors, clusters);
  auto *newNode = new HKInternalNode;
  for (auto &cluster : clusters) {
    auto centroid = std::move(cluster.centroid);
    centroid.resize(vectors[0].size());
    newNode->children.push_back(std::make_pair(nodes.size(), std::move(centroid)));
    auto *cnode = new HKLeafNode;
    nodes.push_back(cnode);
    for (auto &member : cluster.members) {
      if (member.vectorID > leafNode.members.size()) {
	std::cerr << "Fatal error. member:" << member.vectorID << ":" << leafNode.members.size() << std::endl;
	abort();
      }
      if (member.vectorID == leafNode.members.size()) {
	cnode->members.push_back(id);
      } else {
	cnode->members.push_back(leafNode.members[member.vectorID]);
      }
    }
  }
  delete nodes[leafNodeID]; 
  nodes[leafNodeID] = newNode; 
}

double computeError(std::vector<HKNode*> &nodes, NGT::ObjectSpace &objectSpace, NGTQ::Quantizer::ObjectList &objectList) {
  std::cerr << "node size=" << nodes.size() << std::endl;
  double distance = 0.0;
  size_t dcount = 0;
  for (auto *node : nodes) {
    if (node->leaf) {
    } else {
      HKInternalNode &internalNode = static_cast<HKInternalNode&>(*node);
      NGT::Object obj(&objectSpace);
      for (auto &child : internalNode.children) {
	if (nodes[child.first]->leaf) {
	  if (dcount % 100000 == 0) {
	    std::cerr << "Processed leaves=" << dcount << std::endl;
	  }
	  auto centroid = child.second;
	  HKLeafNode &leafNode = static_cast<HKLeafNode&>(*nodes[child.first]);
	  for (auto &m : leafNode.members) {
	    objectList.get(m, obj, &objectSpace);
	    std::vector<float> o = objectSpace.getObject(obj);
	    distance += NGT::Clustering::distanceL2(centroid, o);
	    dcount++;
	  }
	}
      }
    } 
  }
  distance /= dcount;
  std::cout << "# of vectors=" << dcount << std::endl;
  std::cout << "Quantization error=" << distance << std::endl;
  return distance;
}


size_t extractCentroids(std::ostream &oStream, std::vector<HKNode*> &nodes) {
  std::cerr << "node size=" << nodes.size() << std::endl;
  size_t clusterCount = 0;
  size_t objectCount = 0;
  size_t leafID = 0;
  for (auto *node : nodes) {
    if (node->leaf) {
      HKLeafNode &leafNode = static_cast<HKLeafNode&>(*node);
      objectCount += leafNode.members.size();
    } else {
      HKInternalNode &internalNode = static_cast<HKInternalNode&>(*node);
      for (auto &child : internalNode.children) {
	if (nodes[child.first]->leaf) {
	  if (static_cast<HKLeafNode*>(nodes[child.first])->id == 0) {
	    static_cast<HKLeafNode*>(nodes[child.first])->id = leafID;
	  } else if (static_cast<HKLeafNode*>(nodes[child.first])->id != leafID) {
	    std::cerr << "leaf ID is invalid?" << std::endl;
	  }
	  leafID++;
	  size_t count = 0;
	  clusterCount++;
	  for (auto &v : child.second) {
	    oStream << v;
	    if (++count == child.second.size()) {
	      oStream << std::endl;
	    } else {
	      oStream << "\t";
	    }
	  }
	}
      }
    } 
  }
  std::cerr << "# of clusters=" << clusterCount << std::endl;
  return objectCount;
}

size_t extractIndex(std::ostream &oStream, std::vector<HKNode*> &nodes, size_t numOfObjects) {
  std::vector<int32_t> clusterID(numOfObjects, -1);
  std::cerr << "numOfObjects=" << numOfObjects << std::endl;
  std::cerr << "node size=" << nodes.size() << std::endl;
  for (auto *node : nodes) {
    if (node->leaf) {
      HKLeafNode &leafNode = static_cast<HKLeafNode&>(*node);
      for (auto &member : leafNode.members) {
	if (member > numOfObjects) {
	  std::cerr << "output index: Internal fatal error. " << member << ":" << numOfObjects - 1 << std::endl;
	  abort();
	}
	if (member == 0) {
	  std::cerr << "output index: Internal fatal error. Invalid ID" << std::endl;
	  abort();
	}
	clusterID[member - 1] = leafNode.id;
      }
    }
  }
  std::cerr << "clusterID.size=" << clusterID.size() << std::endl;
  size_t count = 0;
  for (auto cid : clusterID) {
    count++;
    oStream << cid << std::endl;
  }
  std::cerr << "# of id=" << count << std::endl;
  return count;
}

void extractBtoQAndQCentroid(std::ostream &btoqStream, std::ostream &qStream,
			       std::vector<HKNode*> &nodes, size_t numOfThirdClusters) {
  std::cerr << "extractBtoQ" << std::endl;
  std::vector<int32_t> btoq(numOfThirdClusters);
  std::cerr << "numOfThirdClusters=" << numOfThirdClusters << std::endl;
  std::cerr << "node size=" << nodes.size() << std::endl;
  size_t rootID = 0;
  HKInternalNode &root = static_cast<HKInternalNode&>(*nodes[rootID]);
  std::cerr << "first=" << root.children.size() << std::endl;
  size_t secondCount = 0;	
  size_t thirdCount = 0;	
  size_t objectCount = 0;
  size_t leafID = 0;
  size_t qID = 0;
  for (auto &c1 : root.children) {
    HKInternalNode &node1 = static_cast<HKInternalNode&>(*nodes[c1.first]);
    std::cerr << "second=" << node1.children.size() << std::endl;
    secondCount += node1.children.size();
    for (auto &c2 : node1.children) {
      HKInternalNode &node2 = static_cast<HKInternalNode&>(*nodes[c2.first]);
      std::cerr << "third=" << node2.children.size() << std::endl;
      thirdCount += node2.children.size();
      size_t count = 0;
      for (auto &v : c2.second) {
	qStream << v;
	if (++count == c2.second.size()) {
	  qStream << std::endl;
	} else {
	  qStream << "\t";
	}
      }
      for (auto &c3 : node2.children) {
	btoqStream << qID << std::endl;
	HKLeafNode &leaf = static_cast<HKLeafNode&>(*nodes[c3.first]);
	objectCount += leaf.members.size();
	if (leaf.id != leafID++) {
	  std::cerr << "leaf is invalid" << leaf.id << ":" << leafID << std::endl;
	  abort();
	}
      }
      qID++;
    }
  }
  std::cerr << "second=" << secondCount << std::endl;
  std::cerr << "third=" << thirdCount << std::endl;
  std::cerr << "object=" << objectCount << std::endl;
}

void extractRandomObjectsFromEachBlob(std::ostream &oStream, std::vector<HKNode*> &nodes, size_t numOfObjects,
					size_t numOfRandomObjects, NGTQ::QuantizerInstance<uint8_t>& quantizer, char rmode) {
  std::cerr << "node size=" << nodes.size() << std::endl;
  std::vector<std::vector<std::vector<float>>> randomObjects(numOfObjects);
  std::vector<std::vector<float>> centroids(numOfObjects);
  for (auto *node : nodes) {
    if (node->leaf) {
      HKLeafNode &leafNode = static_cast<HKLeafNode&>(*node);
      std::vector<uint32_t> randomObjectIDXs;
      if (numOfRandomObjects >= leafNode.members.size()) {
	randomObjectIDXs = leafNode.members;
	while (randomObjectIDXs.size() < numOfRandomObjects) {
	  double random = ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
	  uint32_t idx = floor(leafNode.members.size() * random);
	  if (idx >= leafNode.members.size()) {
	    std::cerr << "Internal error. " << idx << ":" << leafNode.members.size() << std::endl;
	    abort();
	  }
	  randomObjectIDXs.push_back(leafNode.members[idx]);	  
	}
      } else {
	srand(leafNode.id);
	while (randomObjectIDXs.size() < numOfRandomObjects) {
	  uint32_t idx = 0;
	  do {
	    double random = ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
	    idx = floor(leafNode.members.size() * random);
	    if (idx >= leafNode.members.size()) {
	      std::cerr << "Internal error. " << idx << ":" << leafNode.members.size() << std::endl;
	      abort();
	    }
	  } while (std::find(randomObjectIDXs.begin(), randomObjectIDXs.end(), leafNode.members[idx]) != randomObjectIDXs.end());
	  std::cerr << "IDX=" << idx << "/" << leafNode.members.size() << std::endl;
	  randomObjectIDXs.push_back(leafNode.members[idx]);
	}
      }
      std::cerr << "randomObjectIDXs=" << randomObjectIDXs.size() << std::endl;
      for (auto member : randomObjectIDXs) {
	if (member == 0) {
	  std::cerr << "output index: Internal fatal error. Invalid ID. " <<  member << std::endl;
	  abort();
	}
	NGT::Object object(&quantizer.globalCodebookIndex.getObjectSpace());
	quantizer.objectList.get(member, object, &quantizer.globalCodebookIndex.getObjectSpace());
	std::vector<float> ro = quantizer.globalCodebookIndex.getObjectSpace().getObject(object);
	if (leafNode.id >= numOfObjects) {
	  std::cerr << "Internal error! Wrong leaf ID. " << leafNode.id << ":" << numOfObjects << std::endl;
	  abort();
	}
	randomObjects[leafNode.id].push_back(ro);
      }
    } else { 
      if (rmode == 'c') {
	HKInternalNode &internalNode = static_cast<HKInternalNode&>(*node);
	for (auto &child : internalNode.children) {
	  if (nodes[child.first]->leaf) {
	    HKLeafNode &leafNode = static_cast<HKLeafNode&>(*nodes[child.first]);
	    centroids[leafNode.id] = child.second;
	  }
	}
      }
    } 
  }
  for (size_t idx = 0; idx < centroids.size(); idx++) {
    auto &c = centroids[idx];
    if (rmode == 'c' && c.empty()) {
      std::cerr << "NGTLQGCommand: Fatal error! The centroid is empty." << std::endl;
      abort();
    }
    for (size_t i = 0; i < c.size(); i++) {
      oStream << c[i];
      if (i + 1 != c.size()) {
	oStream << "\t";
      } else {
	oStream << std::endl;;
      }
    }
    auto &ros = randomObjects[idx];
    for (auto &ro : ros) {
      if (ro.empty()) {
	std::cerr << "NGTLQGCommand: Fatal error! The random object vector is empty." << std::endl;
	abort();
      }
      for (size_t i = 0; i < ro.size(); i++) {
	oStream << ro[i];
	if (i + 1 != ro.size()) {
	  oStream << "\t";
	} else {
	  oStream << std::endl;;
	}
      }
    }
  }
}

void extractBtoQIndex(std::ofstream &of, std::vector<HKNode*> &nodes, std::vector<uint32_t> &qNodeIDs) {
  size_t leafID = 0;
  for (size_t qnidx = 0; qnidx < qNodeIDs.size(); qnidx++) {
    if (nodes[qNodeIDs[qnidx]]->leaf) {
      std::cerr << "Fatal error. this should be an internal node." << std::endl;
      abort();
    }
    HKInternalNode &inode = static_cast<HKInternalNode&>(*nodes[qNodeIDs[qnidx]]);
    for (auto &c : inode.children) {
      if (!nodes[c.first]->leaf) {
	std::cerr << "Fatal error. this should be a leaf." << std::endl;
	abort();
      }
      HKLeafNode &leaf = static_cast<HKLeafNode&>(*nodes[c.first]);
      if (leaf.id == 0) {
	leaf.id = leafID;
      }
      of << qnidx << std::endl;
      leafID++;
    }
  }
}

void hierarchicalKmeans(uint32_t id, int32_t rootID, NGT::Object &object,
			PseudoRandomizedObjectList &objectList, NGT::ObjectSpace &objectSpace,
			std::vector<HKNode*> &nodes,   NGT::Clustering &clustering, size_t maxSize) {
  NGT::Timer timer;
  objectList.get(id, object, &objectSpace);
  int32_t nodeID = hierarchcalKmeansSearch(nodes, rootID, reinterpret_cast<float*>(&object[0]));
  if (nodeID < 0) {
    std::cerr << "Fatal inner error! node ID=" << nodeID << std::endl;
    exit(1);
  }
  auto *node = nodes[nodeID];
  HKLeafNode &leafNode = static_cast<HKLeafNode&>(*node);
  if (leafNode.members.size() >= maxSize) {
    NGT::Timer subtimer;
    subtimer.start();
    auto objectVector = objectSpace.getObject(object);
    std::vector<std::vector<float>> vectors;
    aggregateObjects(leafNode, vectors, objectSpace, objectList, objectVector);
    subtimer.stop();
    std::cerr << "aggregate time=" << subtimer << std::endl;
    subtimer.start();
    hierarchicalKmeansSplit(id, vectors, nodes, nodeID, clustering);
    subtimer.stop();
    std::cerr << "split time=" << subtimer << std::endl;
  } else {
    leafNode.members.push_back(id);
  }
}


void hierarchicalKmeansBatch(std::vector<uint32_t> &batch, std::vector<pair<uint32_t, uint32_t>> &exceededLeaves,
			     int32_t rootID, NGT::Object &object,
			     PseudoRandomizedObjectList &objectList, NGT::ObjectSpace &objectSpace,
			     std::vector<HKNode*> &nodes, NGT::Clustering &clustering, size_t maxSize, size_t &nleaves,
			     size_t maxExceededLeaves) {

  if (batch.size() == 0) {
    return;
  }

  int32_t nodeIDs[batch.size()];

#pragma omp parallel for
  for (size_t idx = 0; idx < batch.size(); idx++) {
    auto id = batch[idx];
    #pragma omp critical
    objectList.get(id, object, &objectSpace);
    int32_t nodeID = hierarchcalKmeansSearch(nodes, rootID, reinterpret_cast<float*>(&object[0]));
    if (nodeID < 0) {
      std::cerr << "Fatal inner error! node ID=" << nodeID << std::endl;
      exit(1);
    }
    nodeIDs[idx] = nodeID;
  }

  
  for (size_t idx = 0; idx < batch.size(); idx++) {
    auto id = batch[idx];
    HKLeafNode &leafNode = static_cast<HKLeafNode&>(*nodes[nodeIDs[idx]]);
    leafNode.members.push_back(id);
    if (leafNode.members.size() > maxSize) {
      auto i = exceededLeaves.begin();
      for (; i != exceededLeaves.end(); i++) {
	if (static_cast<int32_t>((*i).second) == nodeIDs[idx]) break;
      }
      if (i == exceededLeaves.end()) {
	exceededLeaves.push_back(std::make_pair(batch[idx], nodeIDs[idx]));
      }
    }
  }

  batch.clear();

  if (exceededLeaves.size() < maxExceededLeaves) {
    return;
  }
	
  std::vector<std::vector<NGT::Clustering::Cluster>> clusters(exceededLeaves.size());
#pragma omp parallel for
  for (size_t idx = 0; idx < exceededLeaves.size(); idx++) {
    HKLeafNode &leafNode = static_cast<HKLeafNode&>(*nodes[exceededLeaves[idx].second]);
    std::vector<std::vector<float>> vectors;
    #pragma omp critical
    aggregateObjects(leafNode, vectors, objectSpace, objectList);
    clustering.kmeans(vectors, clusters[idx]);
  }

  std::cerr << "exceeded leaves=" << exceededLeaves.size() << std::endl;
  for (size_t idx = 0; idx < exceededLeaves.size(); idx++) {
    auto leafNodeID = exceededLeaves[idx].second;
    HKLeafNode &leafNode = static_cast<HKLeafNode&>(*nodes[leafNodeID]);
    auto *newNode = new HKInternalNode;
    for (auto &cluster : clusters[idx]) {
      newNode->children.push_back(std::make_pair(nodes.size(), std::move(cluster.centroid)));
      auto *cnode = new HKLeafNode;
      nodes.push_back(cnode);
      for (auto &member : cluster.members) {
	cnode->members.push_back(leafNode.members[member.vectorID]);
      }
    }
    nleaves += clusters[idx].size() - 1;
    delete nodes[leafNodeID]; 
    nodes[leafNodeID] = newNode; 
  }
  exceededLeaves.clear();

}

void hierarchicalKmeansWithNumberOfClusters(size_t numOfTotalClusters, size_t numOfObjects, size_t numOfLeaves, 
					    PseudoRandomizedObjectList &objectList, NGT::ObjectSpace &objectSpace, 
					    std::vector<HKNode*> &nodes, NGT::Clustering::InitializationMode initMode){
  std::cerr << "numOfTotalClusters=" << numOfTotalClusters << std::endl;
  std::cerr << "numOfLeaves=" << numOfLeaves << std::endl;
  if (numOfLeaves > numOfTotalClusters) {
    std::cerr << "# of clusters is invalid. " << numOfLeaves << ":" << numOfTotalClusters << std::endl;
    abort();
  }
  auto numOfRemainingClusters = numOfTotalClusters;
  auto numOfRemainingVectors = numOfObjects;
  size_t leafCount = 0;
  size_t nodeSize = nodes.size(); 
  for (size_t nidx = 0; nidx < nodeSize; nidx++) {
    if (nodes[nidx]->leaf) {
      leafCount++;
      if (numOfLeaves >= 100 && leafCount % (numOfLeaves / 100) == 0) {
	std::cerr << "Processed leaves: " << leafCount << " " << leafCount * 100 / numOfLeaves << "%" << std::endl;
      }
      HKLeafNode &leafNode = static_cast<HKLeafNode&>(*nodes[nidx]);
      std::vector<std::vector<float>> vectors;
      aggregateObjects(leafNode, vectors, objectSpace, objectList);
      size_t nClusters = round(static_cast<float>(leafNode.members.size()) / numOfRemainingVectors * numOfRemainingClusters);
      nClusters = nClusters == 0 ? 1 : nClusters;
      numOfRemainingVectors -= leafNode.members.size();
      numOfRemainingClusters -= nClusters;
      NGT::Clustering clustering(initMode, NGT::Clustering::ClusteringTypeKmeansWithoutNGT, 1000, nClusters);
      NGT::Timer timer;
      timer.start();
      hierarchicalKmeansSplit(0, vectors, nodes, nidx, clustering);
      timer.stop();
      if (nodes[nidx]->leaf) {
	std::cerr << "At this moment, the second node should be an internal" << std::endl;
	abort();
      }
    }
  } 
}

void hierarchicalKmeansWithNumberOfClustersInParallel(size_t numOfTotalClusters, size_t numOfObjects, size_t numOfLeaves, 
						      PseudoRandomizedObjectList &objectList, NGT::ObjectSpace &objectSpace, 
						      std::vector<HKNode*> &nodes, NGT::Clustering::InitializationMode initMode){
  NGT::Timer timer;
  timer.start();
  auto numOfRemainingClusters = numOfTotalClusters;
  auto numOfRemainingVectors = numOfObjects;
  size_t leafCount = 0;

  std::vector<pair<uint32_t, size_t>> leafNodes;
  leafNodes.reserve(numOfLeaves);
  for (size_t nidx = 0; nidx < nodes.size(); nidx++) {
    if (nodes[nidx]->leaf) {
      leafCount++;
      {
	size_t step = 10;
	if (numOfLeaves >= step && leafCount % (numOfLeaves / step) == 0) {
	  std::cerr << "Processed leaves: " << leafCount << " " << leafCount * step / numOfLeaves << "%" << std::endl;
	}
      }
      HKLeafNode &leafNode = static_cast<HKLeafNode&>(*nodes[nidx]);
      size_t nClusters = round(static_cast<float>(leafNode.members.size()) / numOfRemainingVectors * numOfRemainingClusters);
      nClusters = nClusters == 0 ? 1 : nClusters;
      numOfRemainingVectors -= leafNode.members.size();
      numOfRemainingClusters -= nClusters;
      leafNodes.push_back(std::make_pair(nidx, nClusters));
    }
  }
  timer.stop();
  std::cerr << "hierarchicalKmeansWithNumberOfClustersInParallel: extract leaves. Time=" << timer << std::endl;
  timer.start();

  std::cerr << "start kmeans..." << std::endl;
  std::vector<std::vector<NGT::Clustering::Cluster>> clusters(leafNodes.size());
#pragma omp parallel for
  for (size_t nidx = 0; nidx < leafNodes.size(); nidx++) {
    HKLeafNode &leafNode = static_cast<HKLeafNode&>(*nodes[leafNodes[nidx].first]);
    std::vector<std::vector<float>> vectors;
    #pragma omp critical
    aggregateObjects(leafNode, vectors, objectSpace, objectList);
    NGT::Clustering clustering(initMode, NGT::Clustering::ClusteringTypeKmeansWithoutNGT, 1000, leafNodes[nidx].second);
    clustering.kmeans(vectors, clusters[nidx]);
  }

  timer.stop();
  std::cerr << "hierarchicalKmeansWithNumberOfClustersInParallel: kmeans. Time=" << timer << std::endl;
  timer.start();

  std::cerr << "add nodes..." << std::endl;
  for (size_t idx = 0; idx < leafNodes.size(); idx++) {
    auto leafNodeID = leafNodes[idx].first;
    HKLeafNode &leafNode = static_cast<HKLeafNode&>(*nodes[leafNodeID]);
    auto *newNode = new HKInternalNode;
    for (auto &cluster : clusters[idx]) {
      newNode->children.push_back(std::make_pair(nodes.size(), std::move(cluster.centroid)));
      auto *cnode = new HKLeafNode;
      nodes.push_back(cnode);
      for (auto &member : cluster.members) {
	cnode->members.push_back(leafNode.members[member.vectorID]);
      }
    }
    delete nodes[leafNodeID]; 
    nodes[leafNodeID] = newNode; 
  }
  timer.stop();
  std::cerr << "hierarchicalKmeansWithNumberOfClustersInParallel: add nodes. Time=" << timer << std::endl;

}

void flattenClusters(std::vector<NGT::Clustering::Cluster> &upperClusters, 
		     std::vector<std::vector<NGT::Clustering::Cluster>> &lowerClusters,
		     size_t numOfLowerClusters,
		     std::vector<NGT::Clustering::Cluster> &flatClusters) {


  flatClusters.clear();
  flatClusters.reserve(numOfLowerClusters);

  for (size_t idx1 = 0; idx1 < lowerClusters.size(); idx1++) {
    for (size_t idx2 = 0; idx2 < lowerClusters[idx1].size(); idx2++) {
      for (auto &m : lowerClusters[idx1][idx2].members) {
	m.vectorID = upperClusters[idx1].members[m.vectorID].vectorID;
      }
      flatClusters.push_back(lowerClusters[idx1][idx2]);
    }
  }

}

#ifndef MULTIPLE_OBJECT_LISTS
void subclustering(std::vector<NGT::Clustering::Cluster> &upperClusters, size_t numOfLowerClusters, size_t numOfObjects,
		   NGT::ObjectSpace &objectSpace, PseudoRandomizedObjectList &objectList, 
		   NGT::Clustering::InitializationMode initMode, std::vector<std::vector<NGT::Clustering::Cluster>> &lowerClusters) {
  std::vector<uint32_t> nPartialClusters(upperClusters.size());
  auto numOfRemainingClusters = numOfLowerClusters;
  auto numOfRemainingVectors = numOfObjects;
  size_t ts = 0;
  for (size_t idx = 0; idx < upperClusters.size(); idx++) {
    size_t ncs = round(static_cast<float>(upperClusters[idx].members.size()) / numOfRemainingVectors * 
		       numOfRemainingClusters);
    ncs = ncs == 0 ? 1 : ncs;
    numOfRemainingVectors -= upperClusters[idx].members.size();
    if (numOfRemainingClusters >= ncs) {
      numOfRemainingClusters -= ncs;
    }
    nPartialClusters[idx] = ncs;
    ts += ncs;
  }
  std::cerr << "numOfRemainingClusters=" << numOfRemainingClusters << std::endl;
  std::cerr << "numOfRemainingVectors=" << numOfRemainingVectors << std::endl;
  std::cerr << "upperClusters=" << upperClusters.size() << std::endl;
  std::cerr << "total=" << ts << ":" << numOfLowerClusters << std::endl;
  if (ts < numOfLowerClusters || numOfRemainingClusters != 0) {
    std::cerr << "subclustering: Internal error! " << std::endl;
    exit(1);
  }

  lowerClusters.resize(upperClusters.size());
#pragma omp parallel for schedule(dynamic)
  for (size_t idx = 0; idx < upperClusters.size(); idx++) {
    std::cerr << "cluster idx=" << idx << " # of clusters=" << nPartialClusters[idx] << " # of members=" << upperClusters[idx].members.size() << std::endl;
    std::vector<std::vector<float>> partialVectors;
    partialVectors.reserve(upperClusters[idx].members.size());
    NGT::Object obj(&objectSpace);
#pragma omp critical
    {
      for (auto &m : upperClusters[idx].members) {
	objectList.get(m.vectorID + 1, obj, &objectSpace);
	partialVectors.push_back(objectSpace.getObject(obj));
      }
    }
    if (upperClusters[idx].members.size() != partialVectors.size()) {
      std::cerr << "the sizes of members are not consistent" << std::endl;
      abort();
    }
    NGT::Clustering lowerClustering(initMode, NGT::Clustering::ClusteringTypeKmeansWithoutNGT, 1000);
    lowerClustering.kmeans(partialVectors, nPartialClusters[idx], lowerClusters[idx]);
    if (nPartialClusters[idx] != lowerClusters[idx].size()) {
      std::cerr << "the sizes of cluster members are not consistent" << std::endl;	    
      abort();
    }
    std::cerr << "# of clusters=" << lowerClusters[idx].size() << std::endl;	  
  }
  size_t nc = 0;
  size_t mc = 0;
  for (auto &cs : lowerClusters) {
    nc += cs.size();
    for (auto &c : cs) {
      mc += c.members.size();
    }
  }
  std::cerr << "# of clusters=" << nc << " # of members=" << mc << std::endl;
}
void subclustering(std::vector<NGT::Clustering::Cluster> &upperClusters, size_t numOfLowerClusters, size_t numOfObjects,
		   NGT::ObjectSpace &objectSpace, PseudoRandomizedObjectList &objectList, 
		   NGT::Clustering::InitializationMode initMode, std::vector<NGT::Clustering::Cluster> &flatLowerClusters) {

  std::vector<std::vector<NGT::Clustering::Cluster>> lowerClusters;
  subclustering(upperClusters, numOfLowerClusters, numOfObjects, objectSpace, objectList, initMode, lowerClusters);

  flattenClusters(upperClusters, lowerClusters, numOfLowerClusters, flatLowerClusters);

}

#else 
void subclustering(std::vector<NGT::Clustering::Cluster> &upperClusters, size_t numOfLowerClusters, size_t numOfObjects,
		   NGT::ObjectSpace &objectSpace, PseudoRandomizedObjectList &objectList, 
		   NGT::Clustering::InitializationMode initMode, std::vector<std::vector<NGT::Clustering::Cluster>> &lowerClusters) {
  std::vector<uint32_t> nPartialClusters(upperClusters.size());
  auto numOfRemainingClusters = numOfLowerClusters;
  auto numOfRemainingVectors = numOfObjects;
  size_t ts = 0;
  for (size_t idx = 0; idx < upperClusters.size(); idx++) {
    size_t ncs = round(static_cast<float>(upperClusters[idx].members.size()) / numOfRemainingVectors * 
		       numOfRemainingClusters);
    ncs = ncs == 0 ? 1 : ncs;
    numOfRemainingVectors -= upperClusters[idx].members.size();
    if (numOfRemainingClusters >= ncs) {
      numOfRemainingClusters -= ncs;
    }
    nPartialClusters[idx] = ncs;
    ts += ncs;
  }

  std::cerr << "numOfRemainingClusters=" << numOfRemainingClusters << std::endl;
  std::cerr << "numOfRemainingVectors=" << numOfRemainingVectors << std::endl;
  std::cerr << "upperClusters=" << upperClusters.size() << std::endl;
  std::cerr << "total=" << ts << ":" << numOfLowerClusters << std::endl;
  if (ts < numOfLowerClusters || numOfRemainingClusters != 0) {
    std::cerr << "subclustering: Internal error! " << std::endl;
    exit(1);
  }

  auto nthreads = omp_get_max_threads();
  if (!objectList.multipleOpen(nthreads)) {
    std::cerr << "Cannot open multiple streams." << std::endl;
    abort();
  }

  lowerClusters.resize(upperClusters.size());
#pragma omp parallel for schedule(dynamic)
  for (size_t idx = 0; idx < upperClusters.size(); idx++) {
    std::vector<std::vector<float>> partialVectors;
    partialVectors.reserve(upperClusters[idx].members.size());
    NGT::Object obj(&objectSpace);
    auto threadid = omp_get_thread_num();
//#pragma omp critical
    {
      for (auto &m : upperClusters[idx].members) {
	if (threadid >= nthreads) {
	  std::cerr << "inner fatal error. # of threads=" << nthreads << ":" << threadid << std::endl;
	  exit(1);
	}
	if (!objectList.get(threadid, m.vectorID + 1, obj, &objectSpace)) {
	  std::cerr << "subclustering: Fatal error! cannot get!!!! " << m.vectorID + 1 << std::endl;
	  abort();
	}
	partialVectors.push_back(objectSpace.getObject(obj));
      }
    }
    if (upperClusters[idx].members.size() != partialVectors.size()) {
      std::cerr << "the sizes of members are not consistent" << std::endl;
      abort();
    }
    NGT::Clustering lowerClustering(initMode, NGT::Clustering::ClusteringTypeKmeansWithoutNGT, 1000);
    lowerClustering.kmeans(partialVectors, nPartialClusters[idx], lowerClusters[idx]);
    if (nPartialClusters[idx] != lowerClusters[idx].size()) {
      std::cerr << "the sizes of cluster members are not consistent" << std::endl;	    
      abort();
    }
  }
  size_t nc = 0;
  size_t mc = 0;
  for (auto &cs : lowerClusters) {
    nc += cs.size();
    for (auto &c : cs) {
      mc += c.members.size();
    }
  }
  std::cerr << "# of clusters=" << nc << " # of members=" << mc << std::endl;
}

void subclustering(std::vector<NGT::Clustering::Cluster> &upperClusters, size_t numOfLowerClusters, size_t numOfObjects,
		   NGT::ObjectSpace &objectSpace, PseudoRandomizedObjectList &objectList,
		   NGT::Clustering::InitializationMode initMode, std::vector<NGT::Clustering::Cluster> &flatLowerClusters) {

  std::vector<std::vector<NGT::Clustering::Cluster>> lowerClusters;
  subclustering(upperClusters, numOfLowerClusters, numOfObjects, objectSpace, objectList, initMode, lowerClusters);

  flattenClusters(upperClusters, lowerClusters, numOfLowerClusters, flatLowerClusters);

}

#endif 


void assign(std::vector<NGT::Clustering::Cluster> &clusters, size_t beginID, size_t endID,
	    NGT::ObjectSpace &objectSpace, PseudoRandomizedObjectList &objectList) {

  if (!objectList.multipleOpen(omp_get_max_threads())) {
    std::cerr << "Cannot open multiple streams." << std::endl;
    abort();
  }

  size_t count = 0;
#pragma omp parallel for
  for (size_t id = beginID; id <= endID; id++) {
    NGT::Object obj(&objectSpace);
    //#pragma omp critical
    objectList.get(omp_get_thread_num(), id, obj, &objectSpace);
    vector<float> v = objectSpace.getObject(obj);
    float min = std::numeric_limits<float>::max();
    int minidx = -1;
    for (size_t cidx = 0; cidx != clusters.size(); cidx++) {
      auto d = NGT::PrimitiveComparator::compareL2(reinterpret_cast<float*>(v.data()),
						   clusters[cidx].centroid.data(), v.size());
      if (d < min) {
	min = d;
	minidx = cidx;
      }
    }
    if (minidx < 0) {
      std::cerr << "assign: Fatal error!" << std::endl;
      abort();
    }
#pragma omp critical
    {
      clusters[minidx].members.push_back(NGT::Clustering::Entry(id - 1, minidx, min));
      count++;
      if (count % 1000000 == 0) {
	std::cerr << "# of assigned objects=" << count << std::endl;
      }
    }
  }

}

void assignWithNGT(std::vector<NGT::Clustering::Cluster> &clusters, size_t beginID, size_t endID,
		   NGT::ObjectSpace &objectSpace, PseudoRandomizedObjectList &objectList) {

  NGT::Property prop;
  prop.dimension = objectSpace.getDimension();
  prop.objectType = NGT::Index::Property::ObjectType::Float;
  prop.distanceType = NGT::Property::DistanceType::DistanceTypeL2;
  prop.edgeSizeForCreation = 10;
  prop.edgeSizeForSearch = 40;

  std::cerr << "dim=" << prop.dimension << std::endl;
  NGT::Index index(prop);
  std::cerr << "appending objects into NGT.." << std::endl;
  for (size_t cidx = 0; cidx < clusters.size(); cidx++) {
    if (cidx % 100000 == 0) {
      std::cerr << "# of appended cluster objects=" << cidx << std::endl;
    }
    index.append(clusters[cidx].centroid);
  }
  std::cerr << "createIndex..." << std::endl;
  index.createIndex(500);

  std::cerr << "assign with NGT..." << std::endl;
  endID++;
  size_t batchSize = 10000;
  std::vector<pair<uint32_t, float>> clusterIDs(endID - beginID);

  size_t nOfProcessed = 0;
#pragma omp parallel for
  for (size_t batchidx = beginID; batchidx < endID; batchidx += batchSize) {
    std::vector<std::vector<float>> objects(batchSize);
    size_t batchEnd = endID < batchidx + batchSize ? endID : batchidx + batchSize;
#pragma omp critical
    {
      for (size_t id = batchidx; id < batchEnd; id++) {
	NGT::Object obj(&objectSpace);
	objectList.get(id, obj, &objectSpace);
        auto v = objectSpace.getObject(obj);
	objects[id - batchidx] = v;
      }
      nOfProcessed += batchSize;
      if (nOfProcessed % 1000000 == 0) {
	std::cerr << "# of assigned objects (NGT)=" << nOfProcessed << "/" << endID - beginID << std::endl;
      }
    }
    for (size_t id = batchidx; id < batchEnd; id++) {
      NGT::SearchQuery sc(objects[id - batchidx]);
      NGT::ObjectDistances	objects;
      sc.setResults(&objects);
      sc.setSize(10);
      sc.setEpsilon(0.12);
      index.search(sc);
      clusterIDs[id - beginID] = make_pair(objects[0].id - 1, objects[0].distance);
    }
  }

  std::cerr << "pushing..." << std::endl;
  for (size_t id = beginID; id < endID; id++) {
    auto cid = clusterIDs[id - beginID].first;
    auto cdistance = clusterIDs[id - beginID].second;
    clusters[cid].members.push_back(NGT::Clustering::Entry(id - 1, cid, cdistance));
  }
}

void
NGTLQG::Command::hierarchicalKmeans(NGT::Args &args)
{
  const std::string usage = "ngtlqg ";
  std::string indexPath;

  try {
    indexPath = args.get("#1");
  } catch (...) {
    cerr << "DB is not specified" << endl;
    cerr << usage << endl;
    return;
  }

  std::string prefix;
  try {
    prefix = args.get("#2");
    std::cerr << "prefix=" << prefix << std::endl;
  } catch (...) {}

  std::string objectIDsFile;
  try {
    objectIDsFile = args.get("#3");
    std::cerr << "object IDs=" << objectIDsFile << std::endl;
  } catch (...) {}
  
  size_t maxSize = args.getl("m", 1000);
  size_t numOfObjects = args.getl("n", 0);
  size_t numOfClusters = args.getl("N", 2);
  size_t numOfTotalClusters = args.getl("C", 0);
  size_t numOfTotalBlobs = args.getl("b", 0);
  int32_t clusterID = args.getl("c", -1);
  
  char iMode = args.getChar("i", '-');
  NGT::Clustering::InitializationMode initMode = NGT::Clustering::InitializationModeHead;
  switch (iMode) {
  case '-':
  case 'l':
  case 'h': initMode = NGT::Clustering::InitializationModeHead; break;
  case 'r': initMode = NGT::Clustering::InitializationModeRandom; break;
  case 'p': initMode = NGT::Clustering::InitializationModeKmeansPlusPlus; break;
  }

  size_t numOfRandomObjects = args.getl("r", 0);
  char rmode = args.getChar("R", '-');

  size_t numOfFirstObjects = 0;
  size_t numOfFirstClusters = 0;
  size_t numOfSecondObjects = 0;
  size_t numOfSecondClusters = 0;
  size_t numOfThirdClusters = 0;
  std::string blob = args.getString("B", "");
  if (!blob.empty()) {
    std::vector<std::string> tokens;
    NGT::Common::tokenize(blob, tokens, ",");
    if (tokens.size() > 0) {
      std::vector<std::string> ftokens;
      NGT::Common::tokenize(tokens[0], ftokens, ":");
      if (ftokens.size() >= 2) {
	numOfFirstClusters = NGT::Common::strtof(ftokens[0]);
	numOfFirstObjects = NGT::Common::strtof(ftokens[1]);
      }
    }
    if (tokens.size() > 1) {
      std::vector<std::string> ftokens;
      NGT::Common::tokenize(tokens[1], ftokens, ":");
      if (ftokens.size() >= 1) {
	numOfSecondClusters = NGT::Common::strtof(ftokens[0]);
      }
      if (ftokens.size() >= 2) {
	numOfSecondObjects = NGT::Common::strtof(ftokens[1]);
      }
    }
    if (tokens.size() > 2) {
      std::vector<std::string> ftokens;
      NGT::Common::tokenize(tokens[2], ftokens, ":");
      if (ftokens.size() >= 1) {
	numOfThirdClusters = NGT::Common::strtof(ftokens[0]);
      }
    }
    std::cerr << "blob param=:" << std::endl;
    std::cerr << "numOfFirstObjects=" << numOfFirstObjects << std::endl;
    std::cerr << "numOfFirstClusters=" << numOfFirstClusters << std::endl;
    std::cerr << "numOfSecondClusters=" << numOfSecondClusters << std::endl;
    std::cerr << "numOfSecondObjects=" << numOfSecondObjects << std::endl;
    std::cerr << "numOfThirdClusters=" << numOfThirdClusters << std::endl;
    if (numOfFirstObjects < numOfFirstClusters) {
      std::cerr << "# of objects for the first should be larger than # of the first clusters. " << numOfFirstObjects << ":" << numOfFirstClusters << std::endl;
      abort();
    }
    if (numOfFirstClusters > numOfSecondClusters) {
      std::cerr << "# of the first clusters should be larger than # of the second clusters. " << numOfFirstClusters << ":" << numOfSecondClusters << std::endl;
      abort();
    }
    if (numOfSecondClusters > numOfThirdClusters) {
      std::cerr << "# of the first clusters should be larger than # of the second clusters. " << numOfFirstClusters << ":" << numOfSecondClusters << std::endl;
      abort();
    }
  }

  NGT::Clustering::ClusteringType clusteringType = NGT::Clustering::ClusteringTypeKmeansWithoutNGT;
  
  std::vector<HKNode*> nodes;
  
  NGTLQG::Index index(indexPath, true);
  auto &quantizer = static_cast<NGTQ::QuantizerInstance<uint8_t>&>(index.getQuantizer());
  auto &objectSpace = quantizer.globalCodebookIndex.getObjectSpace();
  PseudoRandomizedObjectList objectList(quantizer.objectList);
  size_t paddedDimension = objectSpace.getPaddedDimension();
  size_t dimension = objectSpace.getDimension();
  if (paddedDimension != dimension) {
    std::cerr << "NGTLQGCommand: Fatal error! Dimensions are inconsistent. Dimension=" << paddedDimension << "" << dimension << std::endl;
  }
  uint32_t rootID = 0;
  nodes.push_back(new HKLeafNode);
  NGT::Object object(&objectSpace);
  size_t iteration = 1000;
  NGT::Clustering clustering(initMode, clusteringType, iteration, numOfClusters);
  numOfObjects = numOfObjects == 0 ? objectList.size() : numOfObjects;
  std::cerr << "object list size=" << objectList.size() << std::endl;
  std::cerr << "maxSize=" << maxSize << std::endl;
  std::cerr << "# of objects=" << numOfObjects << std::endl;
  if (objectIDsFile.empty()) {
    if (numOfFirstObjects == 0) {  
      NGT::Timer timer;
      timer.start();
      std::vector<uint32_t> batch;
      std::vector<pair<uint32_t, uint32_t>> exceededLeaves;
      size_t nleaves = 1;
      size_t nOfThreads = 32;
      for (size_t id = 1; id <= numOfObjects; id++) {
	if (id % (numOfObjects / 100) == 0) {
	  timer.stop();
	  std::cerr << "# of processed objects=" << id << " " << id * 100 / numOfObjects << "% " << timer << " # of leaves=" << nleaves << std::endl;
	  timer.start();
	}
	batch.push_back(id);
	if (batch.size() > 100000) {
	  size_t kmeansBatchSize = nleaves < nOfThreads ? nleaves : nOfThreads;
	  hierarchicalKmeansBatch(batch, exceededLeaves, rootID, object, objectList, objectSpace, nodes, 
				  clustering, maxSize, nleaves, kmeansBatchSize);

	}
      }
      hierarchicalKmeansBatch(batch, exceededLeaves, rootID, object, objectList, objectSpace, nodes,
			      clustering, maxSize, nleaves, 0);

      if (numOfTotalClusters != 0) {
	NGT::Timer timer;
	timer.start();
	size_t numOfLeaves = 0;
	for (auto node : nodes) {
	  if (node->leaf) {
	    numOfLeaves++;
	  }
	}
	std::cerr << "# of nodes=" << nodes.size() << std::endl;
	std::cerr << "# of leaves=" << numOfLeaves << std::endl;
	std::cerr << "clustering for quantization." << std::endl;
	hierarchicalKmeansWithNumberOfClustersInParallel(numOfTotalClusters, numOfObjects, numOfLeaves, 
							 objectList, objectSpace, nodes, initMode);
	if (numOfTotalBlobs != 0) {
	  NGT::Timer timer;
	  timer.start();
	  size_t numOfLeaves = 0;
	  for (auto node : nodes) {
	    if (node->leaf) {
	      numOfLeaves++;
	    }
	  }
	  std::cerr << "# of leaves=" << numOfLeaves << ":" << numOfTotalClusters << std::endl;
	  if (numOfLeaves != numOfTotalClusters) {
	    std::cerr << "# of leaves is invalid " << numOfLeaves << ":" << numOfTotalClusters << std::endl;
	    abort();
	  }
	  {
	    std::ofstream of(prefix + "_qcentroid.tsv");
	    extractCentroids(of, nodes);
	  }
	  std::vector<uint32_t> qNodeIDs;
	  for (uint32_t nid = 0; nid < nodes.size(); nid++) {
	    if (nodes[nid]->leaf) {
	      qNodeIDs.push_back(nid);
	    }
	  }
	  std::cerr << "clustering to make blobs." << std::endl;
	  hierarchicalKmeansWithNumberOfClustersInParallel(numOfTotalBlobs, numOfObjects, numOfTotalClusters, 
							   objectList, objectSpace, nodes, initMode);
	  {
	    std::ofstream of(prefix + "_btoq_index.tsv");
	    extractBtoQIndex(of, nodes, qNodeIDs);
	  }
	}
      } 
      
    } else { 
      {
	std::cerr << "3 layers clustering..." << std::endl;
	std::cerr << "The first layer. " << numOfFirstClusters << ":" << numOfFirstObjects << std::endl;
	if (numOfFirstClusters >= numOfSecondClusters) {
	  std::cerr << "# of the second clusters should be larger than # of the first clusters." << std::endl;
	  return;
	}
	NGT::Clustering firstClustering(initMode, NGT::Clustering::ClusteringTypeKmeansWithoutNGT, 300);
	float clusterSizeConstraint = 5.0;
	firstClustering.setClusterSizeConstraintCoefficient(clusterSizeConstraint);
	std::cerr << "size constraint=" << clusterSizeConstraint << std::endl;
	std::vector<std::vector<float>> vectors;
	vectors.reserve(numOfFirstObjects);
	NGT::Object obj(&objectSpace);
	for (size_t id = 1; id <= numOfFirstObjects; id++) {
	  if (id % 1000000 == 0) {
	    std::cerr << "# of prcessed objects is " << id << std::endl;
	  }
	  objectList.get(id, obj, &objectSpace);
	  vectors.push_back(objectSpace.getObject(obj));
	}
	std::cerr << "Kmeans... " << vectors.size() << std::endl;
	std::vector<NGT::Clustering::Cluster> firstClusters;
	NGT::Timer timer;

	timer.start();
	firstClustering.kmeans(vectors, numOfFirstClusters, firstClusters);
	timer.stop();
	std::cerr << "# of clusters=" << firstClusters.size() << " time=" << timer << std::endl;
	std::cerr << "vmsize=" << NGT::Common::getProcessVmSize() << std::endl;
	for (size_t i = 0; i < firstClusters.size(); i++) {
	  std::cerr << "1st clusters ID=" << i << " # of members=" << firstClusters[i].members.size() << std::endl;
	}

	std::vector<std::vector<float>> otherVectors;
	timer.start();
	std::cerr << "Assign for the second. (" << numOfFirstObjects << "-" << numOfSecondObjects << ")..." << std::endl;
	::assign(firstClusters, numOfFirstObjects + 1, numOfSecondObjects, objectSpace, objectList);
	timer.stop();
	std::cerr << "Assign(1) time=" << timer << std::endl;
	////////////////////////////////////////
	std::cerr << "subclustering for the second." << std::endl;
	std::vector<NGT::Clustering::Cluster> secondClusters;
	timer.start();
#ifdef MULTIPLE_OBJECT_LISTS
	subclustering(firstClusters, numOfSecondClusters, numOfSecondObjects, objectSpace, objectList, initMode, secondClusters);
#else
	subclustering(firstClusters, numOfSecondClusters, numOfSecondObjects, objectSpace, objectList, initMode, secondClusters);
#endif
	timer.stop();
	std::cerr << "subclustering(1) time=" << timer << std::endl;
	std::cerr << "save quantization centroid" << std::endl;
	NGT::Clustering::saveClusters(prefix + "_qcentroid.tsv", secondClusters);
	timer.start();
	std::cerr << "Assign for the third. (" << numOfSecondObjects << "-" << numOfObjects << ")..." << std::endl;
	::assignWithNGT(secondClusters, numOfSecondObjects + 1, numOfObjects, objectSpace, objectList);
	timer.stop();
	std::cerr << "Assign(2) time=" << timer << std::endl;
	std::cerr << "subclustering for the third." << std::endl;
	std::vector<std::vector<NGT::Clustering::Cluster>> thirdClusters;
	timer.start();
#ifdef MULTIPLE_OBJECT_LISTS
	subclustering(secondClusters, numOfThirdClusters, numOfObjects, objectSpace, objectList, initMode, thirdClusters);
#else
	subclustering(secondClusters, numOfThirdClusters, numOfObjects, objectSpace, objectList, initMode, thirdClusters);
#endif
	timer.stop();
	std::cerr << "subclustering(2) time=" << timer << std::endl;
	{
	  std::vector<size_t> bqindex;
	  for (size_t idx1 = 0; idx1 < thirdClusters.size(); idx1++) {
	    for (size_t idx2 = 0; idx2 < thirdClusters[idx1].size(); idx2++) {
	      bqindex.push_back(idx1);
	    }
	  }
	  std::cerr << "save bqindex..." << std::endl;
	  NGT::Clustering::saveVector(prefix + "_bqindex.tsv", bqindex);
	}
	
	std::vector<NGT::Clustering::Cluster> thirdFlatClusters;
	flattenClusters(secondClusters, thirdClusters, numOfThirdClusters, thirdFlatClusters);

	std::cerr << "save centroid..." << std::endl;
	NGT::Clustering::saveClusters(prefix + "_centroid.tsv", thirdFlatClusters);

	{
	  std::vector<size_t> cindex(numOfObjects);
	  for (size_t cidx = 0; cidx < thirdFlatClusters.size(); cidx++) {
	    for (auto mit = thirdFlatClusters[cidx].members.begin(); mit != thirdFlatClusters[cidx].members.end(); ++mit) {
	      size_t vid = (*mit).vectorID;
	      cindex[vid] = cidx;
	    }
	  }
	  std::cerr << "save index..." << std::endl;
	  NGT::Clustering::saveVector(prefix + "_index.tsv", cindex);

	}
	std::cerr << "end of clustering" << std::endl;
	return;
      }
    }
  } else {
    std::cerr << "Cluster ID=" << clusterID << std::endl;
    if (clusterID < 0) {
      std::cerr << "A target cluster ID is not specified." << std::endl;
      std::cerr << usage << std::endl;
      return;
    }
    std::ifstream objectIDs(objectIDsFile);
    if (!objectIDs) {
      std::cerr << "Cannot open the object id file. " << objectIDsFile << std::endl;
    }
    uint32_t id = 1;
    int32_t cid;
    size_t ccount = 0;
    while (objectIDs >> cid) {
      std::cerr << cid << std::endl;
      if (id % 100000 == 0) {
	std::cerr << "# of processed objects=" << id << std::endl;
      }
      if (cid == -1) {
	continue;
      }
      if (cid == clusterID) {
	ccount++;
	::hierarchicalKmeans(id, rootID, object, objectList, objectSpace, nodes, clustering, maxSize);
      }
      id++;
    }
  }
  size_t objectCount = 0;
  if (prefix.empty()) {
    objectCount = extractCentroids(std::cout, nodes);
  } else {
    {
      std::ofstream of(prefix + "_centroid.tsv");
      objectCount = extractCentroids(of, nodes);
    }
    {
      std::ofstream of(prefix + "_index.tsv");
      extractIndex(of, nodes, numOfObjects);
    }
    if (numOfFirstObjects > 0) {
      std::ofstream btoqof(prefix + "_btoq.tsv");
      std::ofstream qcof(prefix + "_qcentroid.tsv");
      extractBtoQAndQCentroid(btoqof, qcof, nodes, numOfThirdClusters);
    }
    if (numOfRandomObjects > 0) {
      std::ofstream of(prefix + "_random_object.tsv");
      if (rmode == 'c') {
	extractRandomObjectsFromEachBlob(of, nodes, numOfObjects, numOfRandomObjects - 1, quantizer, rmode);
      } else {
	extractRandomObjectsFromEachBlob(of, nodes, numOfObjects, numOfRandomObjects, quantizer, rmode);
      }
    }
  }
  if (objectCount != numOfObjects) {
    std::cerr << "# of objects is invalid. " << objectCount << ":" << numOfObjects << std::endl;
  }
}

void
NGTLQG::Command::assign(NGT::Args &args)
{
  const std::string usage = "ngtlqg assign";
  std::string indexPath;

  try {
    indexPath = args.get("#1");
  } catch (...) {
    cerr << "Any index is not specified" << endl;
    cerr << usage << endl;
    return;
  }

  std::string queryPath;

  try {
    queryPath = args.get("#2");
  } catch (...) {
    cerr << "Any query is not specified" << endl;
    cerr << usage << endl;
    return;
  }

  auto epsilon = args.getf("e", 0.1);
  auto numOfObjects = args.getl("n", 10);
  auto mode = args.getChar("m", '-');

  try {
    NGT::Index		index(indexPath);
    NGT::Property	property;
    index.getProperty(property);
    ifstream		is(queryPath);
    if (!is) {
      std::cerr << "Cannot open the query file. " << queryPath << std::endl;
      return;
    }
    string		line;
    while (getline(is, line)) {
      vector<float>	query;
      {
	stringstream	linestream(line);
	while (!linestream.eof()) {
	  float value;
	  linestream >> value;
	  query.push_back(value);
	}
	if (static_cast<size_t>(property.dimension) != query.size()) {
	  std::cerr << "Dimension is invallid. " << property.dimension << ":" << query.size() << std::endl;
	  return;
	}
      }
      NGT::SearchQuery		sc(query);
      NGT::ObjectDistances	objects;
      sc.setResults(&objects);
      sc.setSize(numOfObjects);
      sc.setEpsilon(epsilon);

      index.search(sc);
      if (objects.size() == 0) {
	std::cerr << "The result is empty. Something wrong." << std::endl;
	return;
      }
      if (mode == '-') {
	std::cout << objects[0].id - 1 << std::endl;
      } else {
	std::cout << objects[0].id << std::endl;
      }
    }
  } catch (NGT::Exception &err) {
    cerr << "Error " << err.what() << endl;
    return;
  } catch (...) {
    cerr << "Error" << endl;
    return;
  }


}


#endif 
