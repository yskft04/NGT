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

#include	"NGT/Index.h"
#include	"NGT/NGTQ/Quantizer.h"
#include	"NGT/NGTQ/QuantizedGraph.h"


namespace NGTLQG {
  class SearchContainer : public NGT::SearchContainer {
  public:
    SearchContainer(NGT::Object &f): NGT::SearchContainer(f),
      cutback(0.0), graphExplorationSize(50), exactResultSize(0),
      blobExplorationCoefficient(0.0) {}
    SearchContainer &operator=(SearchContainer &sc) {
      NGT::SearchContainer::operator=(sc);
      cutback = sc.cutback;
      graphExplorationSize = sc.graphExplorationSize;
      exactResultSize = sc.exactResultSize;
      blobExplorationCoefficient = sc.blobExplorationCoefficient;
      return *this;
    }
    void setCutback(float c) { cutback = c; }
    void setGraphExplorationSize(size_t size) { graphExplorationSize = size; }
    void setExactResultSize(size_t esize) { exactResultSize = esize; }
    void setBlobEpsilon(float c) { blobExplorationCoefficient = c + 1.0; }
    float       cutback;
    size_t      graphExplorationSize;
    size_t      exactResultSize;
    float       blobExplorationCoefficient;
  };

  class LargeQuantizedBlobGraphRepository : public NGTQG::QuantizedGraphRepository {
  public:
    LargeQuantizedBlobGraphRepository(NGTQ::Index &quantizedIndex): NGTQG::QuantizedGraphRepository(quantizedIndex){
    }
    
    void construct(NGTQ::Index &quantizedIndex) {
      
      std::cerr << "# of size=" << quantizedIndex.getGlobalCodebookSize() << " # of subspaces=" << numOfSubspaces << std::endl;
      (*this).resize(quantizedIndex.getGlobalCodebookSize());
#pragma omp parallel for
      for (size_t gid = 1; gid < quantizedIndex.getGlobalCodebookSize(); gid++) {
	if (gid % 10000 == 0) {
	  std::cerr << "The number of processed blobs=" << gid << std::endl;
	}
	NGTQ::InvertedIndexEntry<uint16_t> invertedIndexObjects(numOfSubspaces);
	quantizedIndex.getQuantizer().extractInvertedIndexObject(invertedIndexObjects, gid);
	quantizedIndex.getQuantizer().eraseInvertedIndexObject(gid);
	NGTQ::QuantizedObjectProcessingStream quantizedStream(quantizedIndex.getQuantizer(), invertedIndexObjects.size());
	(*this)[gid].ids.reserve(invertedIndexObjects.size());
	for (size_t oidx = 0; oidx < invertedIndexObjects.size(); oidx++) {
	  (*this)[gid].ids.push_back(invertedIndexObjects[oidx].id);
	  for (size_t idx = 0; idx < numOfSubspaces; idx++) {
#ifdef NGTQ_UINT8_LUT
#ifdef NGTQ_SIMD_BLOCK_SIZE
            size_t dataNo = oidx;  
#if defined(NGT_SHARED_MEMORY_ALLOCATOR)
	    abort();
#else
	    quantizedStream.arrangeQuantizedObject(dataNo, idx, invertedIndexObjects[oidx].localID[idx] - 1);
#endif
#else 
	    objectData[idx * noobjs + dataNo] = invertedIndexObjects[oidx].localID[idx] - 1;
#endif 
#else  
	    objectData[idx * noobjs + dataNo] = invertedIndexObjects[oidx].localID[idx];
#endif 
	  }
	}

	(*this)[gid].subspaceID = invertedIndexObjects.subspaceID;
	(*this)[gid].objects = quantizedStream.compressIntoUint4();
      } 
#ifdef LQG_ERROR_EVALUATION
      std::cerr << "Distance(error)=" << sqrt(distance / distanceCount) << " # of objects=" << distanceCount << std::endl;
#endif
    }

  };
  
  class Index : public NGTQ::Index {
  public:
  Index(const std::string &indexPath, bool readOnly = false) :
    NGTQ::Index(indexPath, readOnly), path(indexPath), quantizedBlobGraph(*this) {
      try {
	load();
      } catch (NGT::Exception &err) {
	if (readOnly) {
	  std::cerr << "NGTLQG::Index: No quantized blob graph. " << err.what() << std::endl;
	} else {
	  quantizedBlobGraph.construct(*this);
	}
      }
    }

    ~Index() {}
    
    static void create(const std::string &index, NGTQ::Property &property, 
		       NGT::Property &globalProperty,
#ifdef NGTQ_NEURIPS21		       
		       NGT::Property &localProperty,
		       std::vector<float> *rotation,
		       const std::string &objectFile) {
#else
		       NGT::Property &localProperty) {
#endif
      property.quantizerType = NGTQ::QuantizerTypeLQG;
#ifdef NGTQ_NEURIPS21
      NGTQ::Index::create(index, property, globalProperty, localProperty, rotation, objectFile);
#else
      NGTQ::Index::create(index, property, globalProperty, localProperty);
#endif
    }
    
    void getSeeds(NGT::Index &index, NGT::Object *object, NGT::ObjectDistances &seeds, size_t noOfSeeds) {
      auto &graph = static_cast<NGT::GraphAndTreeIndex&>(index.getIndex());
      NGT::SearchContainer sc(*object);
      sc.setResults(&seeds);
      sc.setSize(noOfSeeds);
      sc.setEpsilon(0.0);
      sc.setEdgeSize(-2);
      graph.search(sc);
    }

    NGT::Distance getDistance(void *objects, std::vector<float> &distances, size_t noOfObjects, NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8 &lut
			      ) {
      auto &quantizedObjectDistance = getQuantizer().getQuantizedObjectDistance();
#ifdef NGTLQG_MIN
      auto min = quantizedObjectDistance(objects, distances.data(), noOfObjects, lut);
#else
      quantizedObjectDistance(objects, distances.data(), noOfObjects, lut);
#endif
#ifdef NGTLQG_MIN
      return min;
#endif
    }

    std::tuple<NGT::Distance, NGT::Distance>
      judge(NGTQG::QuantizedNode &ivi, size_t k, NGT::Distance radius,
	    NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8 &lut,
	    NGT::NeighborhoodGraph::ResultSet &result, size_t &foundCount
	    ) {
      auto noOfObjects = ivi.ids.size();
      float distances[NGTQ::QuantizedObjectProcessingStream::getNumOfAlignedObjects(noOfObjects)];
      auto &quantizedObjectDistance = getQuantizer().getQuantizedObjectDistance();
#ifdef NGTLQG_MIN
      float distance = quantizedObjectDistance(ivi.objects, &distances[0], noOfObjects, lut);
#else
      quantizedObjectDistance(ivi.objects, &distances[0], noOfObjects, lut);
#endif

#ifdef NGTLQG_MIN
      if (distance >= radius) {
	return std::make_pair(distance, radius);
      }
#endif
      bool found = false;
      for (size_t i = 0; i < noOfObjects; i++) {
	if (distances[i] <= radius) {
	  result.push(NGT::ObjectDistance(ivi.ids[i], distances[i]));
	  found = true;
	  if (result.size() > k) {
	    result.pop();
	  }
	  if (result.size() >= k) {
	    radius = result.top().distance;
	  }
	}
      }
      if (found) foundCount++;
#ifdef NGTLQG_MIN
      return std::make_pair(distance, radius);
#else
      return std::make_pair(0.0, radius);
#endif
    }

    ///////////////////
    void searchBlobGraphNaively(NGTLQG::SearchContainer &searchContainer) {
      NGT::Object *query = &searchContainer.object;
      
      auto &quantizer = getQuantizer();
      auto &globalIndex = quantizer.globalCodebookIndex;
      auto &globalGraph = static_cast<NGT::GraphAndTreeIndex&>(globalIndex.getIndex());
      NGT::ObjectDistances seeds;
      getSeeds(globalIndex, query, seeds, 5);

      if (seeds.empty()) {
	std::cerr << "something wrong." << std::endl;
	return;
      }
      size_t seedBlobID = seeds[0].id;
      auto &quantizedObjectDistance = quantizer.getQuantizedObjectDistance();
#ifdef NGTQG_ZERO_GLOBAL
      NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8 lut;
      quantizedObjectDistance.initialize(lut);
      quantizedObjectDistance.createDistanceLookup(*query, 1, lut); 
#else
#if defined(NGTQG_ROTATION)
      quantizedObjectDistance.rotation->mul(static_cast<float*>(query->getPointer()));
#endif
      uint32_t subspaceID = quantizedBlobGraph[seedBlobID].subspaceID;
      
      std::unordered_map<size_t, NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8> luts;
      luts.insert(std::make_pair(subspaceID, NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8()));
      auto lutfi = luts.find(subspaceID);
      quantizedObjectDistance.initialize((*lutfi).second);
      quantizedObjectDistance.createDistanceLookup(*query, subspaceID, (*lutfi).second);

#endif
      size_t visitCount = 1;
      size_t foundCount = 0;

      size_t k = searchContainer.size;
      NGT::Distance radius = FLT_MAX;
      NGT::Distance distance;
      NGT::NeighborhoodGraph::ResultSet result;
#ifdef NGTQG_ZERO_GLOBAL
      std::tie(distance, radius) = judge(quantizedBlobGraph[seedBlobID], k, radius, lut, result);
#else
      std::tie(distance, radius) = judge(quantizedBlobGraph[seedBlobID], k, radius, (*lutfi).second, result, foundCount);
#endif
      NGT::NeighborhoodGraph::UncheckedSet uncheckedBlobs;
      NGT::NeighborhoodGraph::DistanceCheckedSet distanceCheckedBlobs(globalGraph.repository.size());
      distanceCheckedBlobs.insert(seedBlobID);
      if (globalGraph.searchRepository.size() == 0) {
	std::cerr << "graph is empty! Is it read only?" << std::endl;
	abort();
      }
      auto *nodes = globalGraph.searchRepository.data();
      uncheckedBlobs.push(NGT::ObjectDistance(seedBlobID, distance));
      float explorationRadius = radius * searchContainer.explorationCoefficient;
      while (!uncheckedBlobs.empty()) {
	auto targetBlob = uncheckedBlobs.top();
	if (targetBlob.distance > explorationRadius) {
	  break;
	}
	uncheckedBlobs.pop();
	auto *neighbors = nodes[targetBlob.id].data();
	auto noOfEdges = (searchContainer.edgeSize == 0 || searchContainer.edgeSize > static_cast<int64_t>(nodes[targetBlob.id].size())) ?
	                 nodes[targetBlob.id].size() : searchContainer.edgeSize;
	auto neighborend = neighbors + noOfEdges;
;
	for (auto neighbor = neighbors; neighbor < neighborend; neighbor++) {
	  NGT::ObjectID neighborID = neighbor->first;
	  if (distanceCheckedBlobs[neighborID]) {
	    continue;
	  }
	  visitCount++;
	  distanceCheckedBlobs.insert(neighborID);
#ifdef NGTQG_ZERO_GLOBAL
	  std::tie(distance, radius) = judge(quantizedBlobGraph[neighborID], k, radius, lut, result);
#else
          auto luti = luts.find(subspaceID);
          if (luti == luts.end()) {
	    luts.insert(std::make_pair(subspaceID, NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8()));
	    luti = luts.find(subspaceID);
	    quantizedObjectDistance.initialize((*luti).second);
	    quantizedObjectDistance.createDistanceLookup(*query, subspaceID, (*luti).second);
	  }
	  std::tie(distance, radius) = judge(quantizedBlobGraph[neighborID], k, radius, (*luti).second, result, foundCount);
#endif 
	  if (static_cast<float>(foundCount) / visitCount < searchContainer.cutback) {
	    uncheckedBlobs = NGT::NeighborhoodGraph::UncheckedSet();
	    break;
	  }
	  if (distance <= explorationRadius) {
	    uncheckedBlobs.push(NGT::ObjectDistance(neighborID, distance));
	  }
	}
      }

      if (searchContainer.resultIsAvailable()) { 
	searchContainer.getResult().clear();
	searchContainer.getResult().moveFrom(result);
      } else {
	searchContainer.workingResult = result;
      }


    }


    void searchBlobNaively(NGTLQG::SearchContainer &searchContainer) {
      NGT::Object *query = &searchContainer.object;
      
      auto &quantizer = getQuantizer();
      auto &globalIndex = quantizer.globalCodebookIndex;
      NGT::ObjectDistances seeds;
      getSeeds(globalIndex, query, seeds, 1000);
      seeds.resize(10);

      if (seeds.empty()) {
	std::cerr << "something wrong." << std::endl;
	return;
      }

      auto &quantizedObjectDistance = quantizer.getQuantizedObjectDistance();
#if defined(NGTQG_ROTATION)
      quantizedObjectDistance.rotation->mul(static_cast<float*>(query->getPointer()));
#endif
      std::unordered_map<size_t, NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8> luts;
      size_t foundCount = 0;

      size_t k = searchContainer.size;
      NGT::Distance radius = FLT_MAX;
      NGT::Distance distance;
      NGT::NeighborhoodGraph::ResultSet result;

      for (size_t idx = 0; idx < seeds.size(); idx++) {
	auto blobID = seeds[idx].id;
	auto subspaceID = quantizedBlobGraph[blobID].subspaceID;
	auto luti = luts.find(subspaceID);
	//if (!availableLUTs[subspaceID]) {
	if (luti == luts.end()) {
	  luts.insert(std::make_pair(subspaceID, NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8()));
	  luti = luts.find(subspaceID);
	  quantizedObjectDistance.initialize((*luti).second);
	  quantizedObjectDistance.createDistanceLookup(*query, subspaceID, (*luti).second);
	}
	std::tie(distance, radius) = judge(quantizedBlobGraph[blobID], k, radius, (*luti).second, result, foundCount);

      }

      if (searchContainer.resultIsAvailable()) { 
	searchContainer.getResult().clear();
	searchContainer.getResult().moveFrom(result);
      } else {
	searchContainer.workingResult = result;
      }


    }

   void searchBlobGraph(NGTLQG::SearchContainer &searchContainer) {
     auto &globalIndex = getQuantizer().globalCodebookIndex;
     auto &globalGraph = static_cast<NGT::GraphAndTreeIndex&>(globalIndex.getIndex());
     NGT::ObjectDistances	seeds;
     globalGraph.getSeedsFromTree(searchContainer, seeds);
     if (seeds.empty()) {
       globalGraph.getRandomSeeds(globalGraph.repository, seeds, 20);
     }
     searchBlobGraph(searchContainer, seeds);
   }

   void searchBlobGraph(NGTLQG::SearchContainer &searchContainer, NGT::ObjectDistances &seeds) {

    auto &quantizer = getQuantizer();
    auto &globalIndex = quantizer.globalCodebookIndex;
    auto &globalGraph = static_cast<NGT::GraphAndTreeIndex&>(globalIndex.getIndex());
    auto &objectSpace = globalIndex.getObjectSpace();

    if (searchContainer.explorationCoefficient == 0.0) {
      searchContainer.explorationCoefficient = NGT_EXPLORATION_COEFFICIENT;
    }

    const auto requestedSize = searchContainer.size;
    searchContainer.size = std::numeric_limits<uint32_t>::max();

    // setup edgeSize
    size_t edgeSize = globalGraph.getEdgeSize(searchContainer);

    NGT::NeighborhoodGraph::UncheckedSet untracedNodes;

    NGT::NeighborhoodGraph::DistanceCheckedSet distanceChecked(globalGraph.searchRepository.size());
    NGT::NeighborhoodGraph::ResultSet results;

    globalGraph.setupDistances(searchContainer, seeds, NGT::PrimitiveComparator::L2Float::compare);
    std::sort(seeds.begin(), seeds.end());
    NGT::ObjectDistance currentNearestBlob = seeds.front();
    NGT::Distance explorationRadius = searchContainer.blobExplorationCoefficient * currentNearestBlob.distance;
    std::priority_queue<NGT::ObjectDistance, std::vector<NGT::ObjectDistance>, std::greater<NGT::ObjectDistance>> discardedObjects;
    untracedNodes.push(seeds.front());
    distanceChecked.insert(seeds.front().id);
    for (size_t i = 1; i < seeds.size(); i++) {
      untracedNodes.push(seeds[i]);
      distanceChecked.insert(seeds[i].id);
      discardedObjects.push(seeds[i]);
    }
    size_t explorationSize = 1;
    auto &quantizedObjectDistance = quantizer.getQuantizedObjectDistance();
    std::unordered_map<size_t, NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8> luts;
    NGT::Object rotatedQuery(&objectSpace);
    objectSpace.copy(rotatedQuery, searchContainer.object);
    quantizedObjectDistance.rotation->mul(static_cast<float*>(rotatedQuery.getPointer()));
    NGT::Distance radius = std::numeric_limits<NGT::Distance>::max();

    const size_t dimension = objectSpace.getPaddedDimension();
    NGT::ReadOnlyGraphNode *nodes = globalGraph.searchRepository.data();
    NGT::ReadOnlyGraphNode *neighbors = 0;
    NGT::ObjectDistance target;
    const size_t prefetchSize = objectSpace.getPrefetchSize();
    const size_t prefetchOffset = objectSpace.getPrefetchOffset();
    pair<uint64_t, NGT::PersistentObject*> *neighborptr;
    pair<uint64_t, NGT::PersistentObject*> *neighborendptr;
    for (;;) {
      if (untracedNodes.empty() || untracedNodes.top().distance > explorationRadius) {
	explorationSize++;
	auto blobID = currentNearestBlob.id;
	auto subspaceID = quantizedBlobGraph[blobID].subspaceID;
	//if (!availableLUTs[subspaceID]) {
	auto luti = luts.find(subspaceID);
	if (luti == luts.end()) {
	  luts.insert(std::make_pair(subspaceID, NGTQ::QuantizedObjectDistance::DistanceLookupTableUint8()));
	  luti = luts.find(subspaceID);
	  quantizedObjectDistance.initialize((*luti).second);
	  quantizedObjectDistance.createDistanceLookup(rotatedQuery, subspaceID, (*luti).second);
	}
	NGT::Distance blobDistance;
	size_t foundCount;
	std::tie(blobDistance, radius) = judge(quantizedBlobGraph[blobID], requestedSize, 
					       radius, (*luti).second, results, foundCount);

#ifdef NGTLQG_MIN
	if (blobDistance > radius * searchContainer.explorationCoefficient) {
	  break;
        }
#endif // NGTLQG_MIN
	if (explorationSize > searchContainer.graphExplorationSize) {
	  break;
	}
	if (discardedObjects.empty()) {
	  break;
	}
	currentNearestBlob = discardedObjects.top();
	discardedObjects.pop();
	explorationRadius = searchContainer.blobExplorationCoefficient * currentNearestBlob.distance;
	continue;
      }

      target = untracedNodes.top();
      untracedNodes.pop();


      neighbors = &nodes[target.id];
      neighborptr = &(*neighbors)[0];
      size_t neighborSize = neighbors->size() < edgeSize ? neighbors->size() : edgeSize;
      neighborendptr = neighborptr + neighborSize;

      pair<uint64_t, NGT::PersistentObject*>* nsPtrs[neighborSize];
      size_t nsPtrsSize = 0;
#ifndef PREFETCH_DISABLE
      for (; neighborptr < neighborendptr; ++neighborptr) {
#ifdef NGT_VISIT_COUNT
	searchContainer.visitCount++;
#endif
	if (!distanceChecked[(*(neighborptr)).first]) {
	  distanceChecked.insert((*(neighborptr)).first);
          nsPtrs[nsPtrsSize] = neighborptr;
          if (nsPtrsSize < prefetchOffset) {
            unsigned char *ptr = reinterpret_cast<unsigned char*>((*(neighborptr)).second);
	    NGT::MemoryCache::prefetch(ptr, prefetchSize);
          }
          nsPtrsSize++;
        }
      }
#endif
#ifdef PREFETCH_DISABLE
      for (; neighborptr < neighborendptr; ++neighborptr) {
#else
      for (size_t idx = 0; idx < nsPtrsSize; idx++) {
#endif
#ifdef PREFETCH_DISABLE
	if (distanceChecked[(*(neighborptr)).first]) {
	  continue;
	}
	distanceChecked.insert((*(neighborptr)).first);
#else
	neighborptr = nsPtrs[idx]; 
	if (idx + prefetchOffset < nsPtrsSize) {
	  unsigned char *ptr = reinterpret_cast<unsigned char*>((*(nsPtrs[idx + prefetchOffset])).second);
	  NGT::MemoryCache::prefetch(ptr, prefetchSize);
	}
#endif
#ifdef NGT_DISTANCE_COMPUTATION_COUNT
	searchContainer.distanceComputationCount++;
#endif

	NGT::Distance distance = NGT::PrimitiveComparator::L2Float::compare(searchContainer.object.getPointer(), 
									    neighborptr->second->getPointer(), dimension);
	NGT::ObjectDistance r;
	r.set(neighborptr->first, distance);
	untracedNodes.push(r);
	if (distance < currentNearestBlob.distance) {
	  discardedObjects.push(currentNearestBlob);
	  currentNearestBlob = r;
	  explorationRadius = searchContainer.blobExplorationCoefficient * currentNearestBlob.distance;
	} else {
	  discardedObjects.push(r);
	}
      } 
    } 


    if (searchContainer.resultIsAvailable()) { 
      if (searchContainer.exactResultSize > 0) {
	NGT::ObjectDistances &qresults = searchContainer.getResult();
	auto threadid = omp_get_thread_num();
	auto paddedDimension = getQuantizer().globalCodebookIndex.getObjectSpace().getPaddedDimension();
	NGT::ResultPriorityQueue	rs;
	NGT::Object object(&quantizer.globalCodebookIndex.getObjectSpace());
	qresults.resize(results.size());
	size_t idx = results.size();
	while (!results.empty()) {
	  auto r = results.top();
	  results.pop();
	  quantizer.objectList.get(threadid, r.id, object, &quantizer.globalCodebookIndex.getObjectSpace());
	  r.distance = NGT::PrimitiveComparator::compareL2(static_cast<float*>(searchContainer.object.getPointer()),
							   static_cast<float*>(object.getPointer()), paddedDimension);
	  qresults[--idx] = r;
	}
	std::sort(qresults.begin(), qresults.end());
	qresults.resize(searchContainer.exactResultSize);
      } else {
	NGT::ObjectDistances &qresults = searchContainer.getResult();
	qresults.moveFrom(results);
      }
    } else {
      if (searchContainer.exactResultSize > 0) {
	auto threadid = omp_get_thread_num();
	auto paddedDimension = getQuantizer().globalCodebookIndex.getObjectSpace().getPaddedDimension();
	NGT::ResultPriorityQueue	rs;
	NGT::Object object(&quantizer.globalCodebookIndex.getObjectSpace());
	while (!results.empty()) {
	  auto r = results.top();
	  results.pop();
	  quantizer.objectList.get(threadid, r.id, object, &quantizer.globalCodebookIndex.getObjectSpace());
	  r.distance = NGT::PrimitiveComparator::compareL2(static_cast<float*>(searchContainer.object.getPointer()),
							   static_cast<float*>(object.getPointer()), paddedDimension);
	  rs.push(r);
	}
	results = std::move(rs);
      } else {
	searchContainer.workingResult = std::move(results);
      }
    }

   }

    void save() {
      quantizedBlobGraph.save(path);
    }

    void load() {
      if (quantizedBlobGraph.stat(path)) {
	quantizedBlobGraph.load(path);
      } else {
	NGTThrowException("Not found the rearranged inverted index. [" + path + "]");
      }
    }
    
    static void build(const std::string &indexPath,
		      std::vector<std::vector<float>> &quantizerCodebook,
		      std::vector<uint32_t> &codebookIndex,
		      std::vector<uint32_t> &objectIndex,
		      size_t beginID, size_t endID) {
      buildNGTQ(indexPath, quantizerCodebook, codebookIndex, objectIndex, beginID, endID);
      buildNGTLQG(indexPath);
      std::cerr << "NGTQ index is completed. Vmsize=" << std::endl;
      std::cerr << "  vmsize=" << NGT::Common::getProcessVmSizeStr() << std::endl;
      std::cerr << "  peak vmsize=" << NGT::Common::getProcessVmPeakStr() << std::endl;
    }

    static void buildNGTQ(const std::string &indexPath,
		      std::vector<std::vector<float>> &quantizerCodebook,
		      std::vector<uint32_t> &codebookIndex,
		      std::vector<uint32_t> &objectIndex,
		      size_t beginID, size_t endID) {
      NGTQ::Index index(indexPath);
      if (quantizerCodebook.size() == 0 || codebookIndex.size() == 0) {
	if (quantizerCodebook.size() != codebookIndex.size()) {
	  stringstream msg;
	  msg << "The specified codebooks or indexes are invalild " << quantizerCodebook.size() << ":" << codebookIndex.size();
	  NGTThrowException(msg);
	}
	index.createIndex(beginID, endID);
      } else {
	index.createIndex(quantizerCodebook, codebookIndex, objectIndex, beginID, endID);
      }
      std::cerr << "NGTQ index is completed. Vmsize=" << std::endl;
      std::cerr << "  vmsize=" << NGT::Common::getProcessVmSizeStr() << std::endl;
      std::cerr << "  peak vmsize=" << NGT::Common::getProcessVmPeakStr() << std::endl;
      std::cerr << "saving..." << std::endl;
      index.save();
    }

    static void buildNGTLQG(const std::string &indexPath) {
      std::cerr << "build NGTLQG" << std::endl;
      NGTLQG::Index index(indexPath);
      std::cerr << "NGTLQG index is completed." << std::endl;
      std::cerr << "  vmsize=" << NGT::Common::getProcessVmSizeStr() << std::endl;
      std::cerr << "  peak vmsize=" << NGT::Common::getProcessVmPeakStr() << std::endl;
      std::cerr << "saving..." << std::endl;
      index.save();
    }

    const std::string path;
    
    LargeQuantizedBlobGraphRepository quantizedBlobGraph;


  };

} 
