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
using namespace std;

#define NGT_CLUSTERING

#ifdef NGT_CLUSTERING
#include "NGT/Clustering.h"
#else
#include "Cluster.h"
#endif

#include "Matrix.h"


class Args : public map<string, string>
{
public:
  Args(int argc, char **argv):
    option("a:b:c:d:e:f:g:hi:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:"
	   "A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P:Q:R:S:T:U:V:W:X:Y:Z:")
  {
    int opt;
    while ((opt = getopt(argc, argv, option)) != -1) {
      if ((char)opt == 'h') {
	string str;
	str.append(1, (char)opt);
	insert(pair<string, string>(str, ""));
	continue;
      }
      string str;
      str.append(1, (char)opt);
      insert(pair<string, string>(str, string(optarg)));
    }
    for (int i = 0; optind < argc; optind++, i++) {
      stringstream ss;
      ss << "#" << i;
      insert(pair<string, string>(ss.str(), string(argv[optind])));
    }
  }
  string &find(const char *s) { return get(s); }
  char getChar(const char *s, char v) {
    try {
      return get(s)[0];
    } catch (...) {
      return v;
    }
  }
  string getString(const char *s, const char *v) {
    try {
      return get(s);
    } catch (...) {
      return v;
    }
  }
  string &get(const char *s) {
    Args::iterator ai;
    ai = map<string, string>::find(string(s));
    if (ai == this->end()) {
      stringstream msg;
      msg << s << ": Not specified" << endl;
      throw invalid_argument(msg.str());
    }
    return ai->second;
  }
  long getl(const char *s, long v) {
    char *e;
    long val;
    try {
      val = strtol(get(s).c_str(), &e, 10);
    } catch (...) {
      return v;
    }
    if (*e != 0) {
      stringstream msg;
      msg << "ARGS::getl: Illegal string. Option=-" << s << " Specified value=" << get(s) 
	  << " Illegal string=" << e << endl;
      cerr << msg.str() << endl;
      throw invalid_argument(msg.str());
    }
    return val;
  }
  float getf(const char *s, float v) {
    char *e;
    float val;
    try {
      val = strtof(get(s).c_str(), &e);
    } catch (...) {
      return v;
    }
    if (*e != 0) {
      stringstream msg;
      msg << "ARGS::getf: Illegal string. Option=-" << s << " Specified value=" << get(s) 
	  << " Illegal string=" << e << endl;
      cerr << msg.str() << endl;
      throw invalid_argument(msg.str());
    }
    return val;
  }
  const char *option;
};


void
extractSubvector(vector<vector<float>> &vectors, vector<vector<float>> &subvectors, size_t start , size_t size)
{
  size_t vsize = vectors.size();
  subvectors.clear();
  subvectors.resize(vsize);
  for (size_t vidx = 0; vidx < vsize; vidx++) {
    subvectors[vidx].reserve(size);
    for (size_t i = 0; i < size; i++) {
      subvectors[vidx].push_back(vectors[vidx][start + i]);
    }
  }
}

void
catSubvector(vector<vector<float>> &vectors, vector<vector<float>> &subvectors)
{
  if (vectors.size() == 0) {
    vectors.resize(subvectors.size());
  }
  assert(vectors.size() == subvectors.size());
  size_t vsize = vectors.size();
  size_t subvsize = subvectors[0].size();
  size_t size = vectors[0].size() + subvsize;
  for (size_t vidx = 0; vidx < vsize; vidx++) {
    subvectors[vidx].reserve(size);
    for (size_t i = 0; i < subvsize; i++) {
      vectors[vidx].push_back(subvectors[vidx][i]);
    }
  }
}

void
#ifdef NGT_CLUSTERING
extractQuantizedVector(vector<vector<float>> &qvectors, vector<NGT::Clustering::Cluster> &clusters) 
#else
extractQuantizedVector(vector<vector<float>> &qvectors, vector<Cluster> &clusters) 
#endif
{
  for (size_t cidx = 0; cidx < clusters.size(); ++cidx) {
    for (auto mit = clusters[cidx].members.begin(); mit != clusters[cidx].members.end(); ++mit) {
      size_t vid = (*mit).vectorID;
      if (vid >= qvectors.size()) {
	qvectors.resize(vid + 1);
      }
      qvectors[vid] = clusters[cidx].centroid;
    }
  }
}

double
squareDistance(vector<vector<float>> &va, vector<vector<float>> &vb)
{
  assert(va.size() == vb.size());
  size_t vsize = va.size();
  assert(va[0].size() == vb[0].size());
  size_t dim = va[0].size();
  double distance = 0;
  for (size_t vidx = 0; vidx < vsize; vidx++) {
    for (size_t i = 0; i < dim; i++) {
      double d = va[vidx][i] - vb[vidx][i];
      distance += d * d;
    }
  }
  distance /= dim * vsize;
  return distance;
}

void evaluate(string global, vector<vector<float>> &vectors, char clusteringType, string &ofile, size_t &numberOfSubspaces, size_t &subvectorSize)
{
  vector<vector<float>> residualVectors;
  {
    // compute residual vectors by global centroids.
#ifdef NGT_CLUSTERING
    vector<NGT::Clustering::Cluster> globalCentroid;
#else
    vector<Cluster> globalCentroid;
#endif
    if (!global.empty()) {
      cerr << "generate residual vectors." << endl;
      try {
#ifdef NGT_CLUSTERING
	NGT::Clustering::loadClusters(global, globalCentroid);
#else
	loadClusters(global, globalCentroid);
#endif
      } catch (...) {
	cerr << "Cannot load vectors. " << global << endl;
	return;
      }
      if (clusteringType == 'k') {
#ifdef NGT_CLUSTERING
	NGT::Clustering::assign(vectors, globalCentroid);
#else
	assign(vectors, globalCentroid);
#endif
      } else {
	cerr << "Using NGT" << endl;
#ifdef NGT_CLUSTERING
	std::cerr << "Not implemented" << std::endl;
	abort();
#else
	assignWithNGT(vectors, globalCentroid);
#endif
      }
      residualVectors.resize(vectors.size());
      cerr << "global centroid size=" << globalCentroid.size() << endl;
      for (size_t cidx = 0; cidx < globalCentroid.size(); ++cidx) {
	for (auto mit = globalCentroid[cidx].members.begin(); mit != globalCentroid[cidx].members.end(); ++mit) {
	  size_t vid = (*mit).vectorID;
	  residualVectors[vid] = vectors[vid];
#ifdef NGT_CLUSTERING
	  NGT::Clustering::subtract(residualVectors[vid], globalCentroid[cidx].centroid);
#else
	  subtract(residualVectors[vid], globalCentroid[cidx].centroid);
#endif
	}
      }
    }
  }

  Matrix<float> R;
  Matrix<float>::load(ofile + "_R.tsv", R);    
  vector<vector<float>> qv(vectors.size());	// quantized vector
  vector<vector<float>> xp;	// residual vector
  if (residualVectors.empty()) {
    xp = vectors;
  } else {
    xp = residualVectors;
  }
  Matrix<float>::mulSquare(xp, R);
  for (size_t m = 0; m < numberOfSubspaces; m++) {
#ifdef NGT_CLUSTERING
    vector<NGT::Clustering::Cluster> subClusters;
#else
    vector<Cluster> subClusters;
#endif
    stringstream str;
    str << ofile << "-" << m << ".tsv";
#ifdef NGT_CLUSTERING
    NGT::Clustering::loadClusters(str.str(), subClusters);
#else
    loadClusters(str.str(), subClusters);
#endif
    vector<vector<float>> subVectors;
    extractSubvector(xp, subVectors, m * subvectorSize, subvectorSize);
    if (clusteringType == 'k') {
#ifdef NGT_CLUSTERING
      NGT::Clustering::assign(subVectors, subClusters);
#else
      assign(subVectors, subClusters);
#endif
    } else {
      cerr << "Using NGT for subvector" << endl;
#ifdef NGT_CLUSTERING
      std::cerr << "not implemented" << std::endl;
      abort();
#else
      assignWithNGT(subVectors, subClusters);
#endif
    }
#ifdef NGT_CLUSTERING
    double distortion = NGT::Clustering::calculateML2(subVectors, subClusters);
#else
    double distortion = calculateML2(subVectors, subClusters);
#endif
    cout << "distortion[" << m << "]=" << distortion << endl;
    vector<vector<float>> subCentroids(vectors.size());
    for (size_t cidx = 0; cidx < subClusters.size(); ++cidx) {
#ifdef NGT_CLUSTERING
      vector<NGT::Clustering::Entry> &members = subClusters[cidx].members;
#else
      vector<Entry> &members = subClusters[cidx].members;
#endif
      for (size_t eidx = 0; eidx < members.size(); ++eidx) {
#ifdef NGT_CLUSTERING
	NGT::Clustering::Entry &entry = members[eidx];
#else
	Entry &entry = members[eidx];
#endif
	assert(cidx == entry.centroidID);
	subCentroids[entry.vectorID] = subClusters[cidx].centroid;
      }
    }
    catSubvector(qv, subCentroids);
  }
#ifdef NGT_CLUSTERING
  double distortion = NGT::Clustering::distanceL2(qv, xp);
#else
  double distortion = distanceL2(qv, xp);
#endif
  cout << "distortion=" << distortion << endl;
  return;
}

void evaluate(vector<vector<float>> &vectors, string &ofile, size_t &numberOfSubspaces, size_t &subvectorSize) {
  cerr << "Evaluate" << endl;
  Matrix<float> R;
  Matrix<float>::load(ofile + "_R.tsv", R);    
  vector<vector<float>> xp = vectors;
  Matrix<float>::mulSquare(xp, R);
  for (size_t m = 0; m < numberOfSubspaces; m++) {
#ifdef NGT_CLUSTERING
    vector<NGT::Clustering::Cluster> subClusters;
#else
    vector<Cluster> subClusters;
#endif
    stringstream str;
    str << ofile << "-" << m << ".tsv";
#ifdef NGT_CLUSTERING
    NGT::Clustering::loadClusters(str.str(), subClusters);
#else
    loadClusters(str.str(), subClusters);
#endif
    vector<vector<float>> subVectors;
    extractSubvector(xp, subVectors, m * subvectorSize, subvectorSize);
#ifdef NGT_CLUSTERING
    NGT::Clustering::assign(subVectors, subClusters);
    double distortion = NGT::Clustering::calculateML2(subVectors, subClusters);
#else
    assign(subVectors, subClusters);
    double distortion = calculateML2(subVectors, subClusters);
#endif
    cout << "distortion[" << m << "]=" << distortion << endl;
    for (size_t cidx = 0; cidx < subClusters.size(); ++cidx) {
      cout << "  members[" << cidx << "]=" << subClusters[cidx].members.size() << endl;
    }
  }
  return;
}

void
generateResidualObjects(string global, vector<vector<float>> &vectors, 
			NGT::Clustering::ClusteringType clusteringType,
			size_t &numberOfSubspaces, size_t &subvectorSize)
{
  if (global.empty()) {
    std::cerr << "A global codebook is not specified!" << std::endl;
    exit(1);
    return;
  }
  vector<vector<float>> residualVectors;
  vector<vector<float>> globalCentroid;
  try {
    NGT::Clustering::loadVectors(global, globalCentroid);
  } catch (...) {
    cerr << "Cannot load vectors. " << global << endl;
    abort();
  }

  NGT::Property property;  property.objectType = NGT::Index::Property::ObjectType::Float;
  property.distanceType = NGT::Index::Property::DistanceType::DistanceTypeL2;
  property.dimension = globalCentroid[0].size();
  NGT::Index index(property);
  for (auto &c : globalCentroid) {
    index.append(c);
  }
  index.createIndex(100);

  for (size_t idx = 0; idx < vectors.size(); idx++) {
    auto &v = vectors[idx];
    if (idx % 10000 == 0) {
      std::cerr << "opq: Processing object ID=" << idx << std::endl;
    }
    NGT::ObjectDistances gc;
    NGT::SearchQuery query(v);
    query.setResults(&gc);
    query.setSize(10);
    query.setEpsilon(0.1);
    index.search(query);
    if (gc.empty()) {
      std::cerr << "inner fatal error. no results! something wrong." << std::endl;
      abort();
    }
    if (gc[0].id == 0 || gc[0].id > globalCentroid.size()) {
      std::cerr << "wrong id " << gc[0].id << ":" <<  globalCentroid.size() << std::endl;
      abort();
    }
    auto gcidx = gc[0].id - 1;
    try {
      NGT::Clustering::subtract(v, globalCentroid[gcidx]);
    } catch (NGT::Exception &err) {
      std::cerr << err.what() << ":" << v.size() << "x" << globalCentroid[gcidx].size() << std::endl;
      abort();
    }
  }
}

void optimizeRotation(
		      size_t iteration,
		      vector<vector<float>> &vectors,
		      Matrix<float> &xt,
		      Matrix<float> &R,
		      Matrix<float> &minR,
		      vector<vector<NGT::Clustering::Cluster>> &minLocalClusters,
		      NGT::Clustering::ClusteringType clusteringType,
		      NGT::Clustering::InitializationMode initMode,
		      size_t numberOfClusters,
		      size_t numberOfSubspaces,
		      size_t subvectorSize,
		      size_t clusterIteration,
		      bool clusterSizeConstraint,
		      size_t convergenceLimitTimes,
		      double &minDistortion,
		      NGT::Timer &timelimitTimer, float timelimit
		      ) {

  minDistortion = DBL_MAX;

  int minIt = 0;
  for (size_t it = 0; it < iteration; it++) {
    vector<vector<float>> xp = vectors;
    Matrix<float>::mulSquare(xp, R);
    float distance = 0.0;
#ifdef NGT_CLUSTERING
    vector<vector<NGT::Clustering::Cluster>> localClusters(numberOfSubspaces);
#else
    vector<vector<Cluster>> localClusters(numberOfSubspaces);
#endif
    vector<vector<float>> subQuantizedVectors[numberOfSubspaces];
#define ERROR_CALCULATION  
#ifdef ERROR_CALCULATION
    vector<float> subvectorDistances(numberOfSubspaces);
#endif
#pragma omp parallel for
    for (size_t m = 0; m < numberOfSubspaces; m++) {
      vector<vector<float>> subVectors;
      extractSubvector(xp, subVectors, m * subvectorSize, subvectorSize);
#ifdef NGT_CLUSTERING
      vector<NGT::Clustering::Cluster> &clusters = localClusters[m];
#else
      vector<Cluster> &clusters = localClusters[m];
#endif
#ifdef NGT_CLUSTERING
      NGT::Clustering clustering(initMode, clusteringType, clusterIteration);
      clustering.clusterSizeConstraint = clusterSizeConstraint;
      clustering.kmeans(subVectors, numberOfClusters, clusters);
#else
      size_t reassign;
      if (clusteringType == 'k') {
	reassign = kmeansClustering(initMode, subVectors, numberOfClusters, clusters, clusterSizeConstraint, 0);
      } else {
	reassign = kmeansClusteringWithNGT(initMode, subVectors, numberOfClusters, clusters, clusterSizeConstraint, 0);    
      }
#endif
      extractQuantizedVector(subQuantizedVectors[m], clusters);
      assert(subQuantizedVectors[m].size() == vectors.size());
      assert(subQuantizedVectors[m][0].size() == subvectorSize);
      // 入力部分ベクトルと量子化部分ベクトルを比較して量子化誤差の計算
#ifdef ERROR_CALCULATION
      double d = squareDistance(subQuantizedVectors[m], subVectors);
      subvectorDistances[m] = d;
#endif
    }
    distance = 0.0;
#ifdef ERROR_CALCULATION
    for (size_t m = 0; m < numberOfSubspaces; m++) {
      distance += subvectorDistances[m];
    }
#endif
    vector<vector<float>> quantizedVectors;
    for (size_t m = 0; m < numberOfSubspaces; m++) {
      catSubvector(quantizedVectors, subQuantizedVectors[m]);
    }
    distance = sqrt(distance / numberOfSubspaces);
    if (minDistortion > distance) {
      minDistortion = distance;
      minR = R;
      minIt = it;
      minLocalClusters = localClusters;
    }
    if (it + 1 > iteration || it - minIt > convergenceLimitTimes) {
      break;
    }
    timelimitTimer.stop();
    if (timelimitTimer.time > timelimit) {
      std::cerr << "opq: Warning. The elapsed time exceeds the limit-time. " << timelimit << ":" << timelimitTimer.time << std::endl;
      timelimitTimer.restart();
      break;
    }
    timelimitTimer.restart();
    Matrix<float> a(xt);
    a.mul(quantizedVectors);
    Matrix<float> u, s, v;
    Matrix<float>::svd(a, u, s, v);
    v.transpose();
    u.mul(v);
    R = u;
  }

}

void optimize(
	      size_t iteration,
	      vector<vector<float>> &vectors,
	      string global,
	      string ofile,
	      NGT::Clustering::ClusteringType clusteringType,
	      NGT::Clustering::InitializationMode initMode,
	      size_t dim,
	      size_t numberOfClusters,
	      size_t numberOfSubspaces,
	      size_t subvectorSize,
	      size_t clusterIteration,
	      bool clusterSizeConstraint,
	      size_t convergenceLimitTimes,
	      vector<Matrix<float>> &rs,
	      vector<vector<vector<NGT::Clustering::Cluster>>> &localClusters,
	      vector<double> &errors,
	      NGT::Timer &timelimitTimer, float timelimit
	      ) {

  generateResidualObjects(global, vectors, clusteringType, numberOfSubspaces, subvectorSize);

  Matrix<float> xt(vectors);
#ifndef BLAS_MATRIX
  xt.transpose();
#endif
  localClusters.resize(rs.size());
  errors.resize(rs.size());
  for (size_t ri = 0; ri < rs.size(); ri++) {
    Matrix<float> optr;
    optimizeRotation(
		     iteration,			
		     vectors,			
		     xt,			
		     rs[ri],			
		     optr,			
		     localClusters[ri],
		     clusteringType,
		     initMode,
		     numberOfClusters,
		     numberOfSubspaces,		
		     subvectorSize,		
		     clusterIteration,		
		     clusterSizeConstraint,
		     convergenceLimitTimes,
		     errors[ri],
		     timelimitTimer, timelimit
		     );
    std::cerr << "opt: R it=" << ri << " error=" << errors[ri] << std::endl;
    rs[ri] = optr;
  }
}

int
main(int argc, char **argv)
{

  string usage = "USAGE: opq -n number-of-clusters -m number-of subspaces [-O t|f] [-s t|f] [-I cluster-iteration] [-t R-max-iteration] [-c convergence-limit-times] vector-file [output-file-prefix]\n"
    "       ope -e E -n number-of-clusters -m number-of subspaces vector-file local-centroid-file [global-centroid-file]";

  Args args(argc, argv);

  string invector;
  try {
    invector = args.get("#0");
  } catch(...) {
    cerr << "opq: input vector is not specified." << endl;
    cerr << usage << endl;
    return 1;
  }

  char execType = args.getChar("e", '-');

  string ofile;
  size_t numberOfClusters = args.getl("n", 0);
  size_t numberOfSubspaces = args.getl("m", 0);
  try {
    ofile = args.get("#1");
  } catch(...) {
    cerr << "output file is not specified." << endl;
    cerr << usage << endl;
    return 1;
  }
  if (numberOfClusters == 0) {
    cerr << "-n adequate number-of-clusters should be specified." << endl;
    cerr << usage << endl;
    return 1;
  };
  if (numberOfSubspaces == 0) {
    cerr << "-m adequate number-of-subspaces should be specified." << endl;
    cerr << usage << endl;
    return 1;
  };

#ifdef NGT_CLUSTERING
  string cType;
  try {
    cType = args.getString("C", "k");
  } catch(...) {}

  NGT::Clustering::ClusteringType clusteringType = NGT::Clustering::ClusteringTypeKmeansWithNGT;;
  if (cType == "k") {
    clusteringType = NGT::Clustering::ClusteringTypeKmeansWithoutNGT;
  } else if (cType == "KS") {
    clusteringType = NGT::Clustering::ClusteringTypeKmeansWithNGT;
  } else if (cType == "i") {
    clusteringType = NGT::Clustering::ClusteringTypeKmeansWithIteration;
  } else {
    cerr << "invalid clustering type. " << cType << endl;
    cerr << usage << endl;
    return 1;
  }
#else
  char clusteringType;
  try {
    clusteringType = args.getChar("C", 'k');
  } catch(...) {}
#endif

  
#ifdef NGT_CLUSTERING
  char iMode = args.getChar("i", '-');
  NGT::Clustering::InitializationMode initMode = NGT::Clustering::InitializationModeHead;
  switch (iMode) {
  case '-':
  case 'l':
  case 'h': initMode = NGT::Clustering::InitializationModeHead; break;
  case 'r': initMode = NGT::Clustering::InitializationModeRandom; break;
  case 'p': initMode = NGT::Clustering::InitializationModeKmeansPlusPlus; break;
  }
#else
  char initMode = args.getChar("i", '-');
#endif
  
  size_t convergenceLimitTimes = args.getl("c", 5);
  size_t iteration = args.getl("t", 100);
  size_t clusterIteration = args.getl("I", 100);

  bool clusterSizeConstraint = false;
  if (args.getChar("s", 'f') == 't') {
    clusterSizeConstraint = true;
  }

  if (args.getChar("O", 'f') != 't') {
    iteration = 1;
  }

  string global;
  try {
    global = args.get("#2");
  } catch(...) {}

  size_t nOfMatrices = args.getl("M", 10);
  float sizeSeed = args.getf("S", 0.1);
  size_t sizeStep = args.getl("X", 2);
  float reject = args.getf("R", 0.9);
  float timelimit = args.getf("L", 24 * 1); 
  timelimit *= 60.0 * 60.0; 

  vector<vector<float>> vectors;

  cerr << "opq: Loading..." << endl;
  try {
#ifdef NGT_CLUSTERING
    NGT::Clustering::loadVectors(invector, vectors);
#else
    loadVectors(invector, vectors);
#endif
  } catch (...) {
    cerr << "opq: Cannot load vectors. " << invector << endl;
    return 1;
  }
  cerr << "opq: End of loading." << endl;
  cerr << "opq: invector=" << invector << endl;
  cerr << "opq: Size=" << vectors.size() << endl;
  if (vectors.size() == 0) {
    std::cerr << "opq: error! the specified vetor file is empty. " << invector << ". the size=" << vectors.size() << std::endl;
    exit(1);
  }
  size_t dim = vectors[0].size();
  size_t subvectorSize = dim / numberOfSubspaces;
  if (dim % numberOfSubspaces != 0) {
    cerr << "# of subspaces (m) is illegal. " << numberOfSubspaces << endl;
    return 0;
  }

  switch (execType) {
  case 'E':
    evaluate(global, vectors, clusteringType, ofile, numberOfSubspaces, subvectorSize);
    return 1;
  case 'e':
    evaluate(vectors, ofile, numberOfSubspaces, subvectorSize);
    return 1;
  }

  //////////////////

  NGT::Timer timelimitTimer;
  timelimitTimer.start();
  vector<vector<vector<NGT::Clustering::Cluster>>> localClusters;
  vector<double> errors;
  bool useEye = false;
  if (nOfMatrices == 0) {
    useEye = true;
    nOfMatrices = 1;
  }
  vector<Matrix<float>> rs(nOfMatrices);
  std::cerr << "opq: # of matrices=" << rs.size() << std::endl;
  if (useEye) {
    rs[0].eye(dim);
  } else {
    for (auto &r: rs) {
      r.randomRotation(dim);
    }
  }
  for (size_t vsize = static_cast<float>(vectors.size()) * sizeSeed; ; vsize *= sizeStep) {

    std::cerr << "opq: vsize=" << vsize << " nOfMatrices=" << nOfMatrices << std::endl;
    auto partialVectors = vectors;
    if (vsize < vectors.size()) {
      partialVectors.resize(vsize);
    }

    optimize(
	     iteration,
	     partialVectors,
	     global,
	     ofile,
	     clusteringType,
	     initMode,
	     dim,
	     numberOfClusters,
	     numberOfSubspaces,
	     subvectorSize,
	     clusterIteration,
	     clusterSizeConstraint,
	     convergenceLimitTimes,
	     rs,
	     localClusters,
	     errors,
	     timelimitTimer, timelimit);
    if (rs.size() != 1) {
      nOfMatrices = static_cast<float>(nOfMatrices) * (1.0 - reject);
      nOfMatrices = nOfMatrices == 0 ? 1 : nOfMatrices;
      vector<pair<double, Matrix<float>*>> sortedErrors;
      for (size_t idx = 0; idx < errors.size(); idx++) {
	sortedErrors.emplace_back(make_pair(errors[idx], &rs[idx]));
      }
      sort(sortedErrors.begin(), sortedErrors.end());
      vector<Matrix<float>> tmp;
      for (size_t idx = 0; idx < nOfMatrices; idx++) {
	tmp.emplace_back(*sortedErrors[idx].second);
      }
      if (tmp.size() != nOfMatrices) {
	std::cerr << "something strange. " << tmp.size() << ":" << nOfMatrices << std::endl;
      }
      rs = std::move(tmp);
    }
    if (vsize >= vectors.size()) {
      break;
    }
  }

  if (rs.size() != 1) {
    std::cerr << "opq: Warning. rs.size=" << rs.size() << std::endl;
  }
  auto minR = std::move(rs[0]);
  auto minLocalClusters = std::move(localClusters[0]);
  /////////////////
  size_t pos = std::distance(std::find(ofile.rbegin(), ofile.rend(), '.'), ofile.rend()) - 1;
  string of = ofile.substr(0, pos);
  Matrix<float>::save(of + "_R.tsv", minR);
  for (size_t m = 0; m < numberOfSubspaces; m++) {
    stringstream str;
    str << of << "-" << m << ".tsv";
#ifdef NGT_CLUSTERING
    NGT::Clustering::saveClusters(str.str(), minLocalClusters[m]);
#else
    saveClusters(str.str(), minLocalClusters[m]);
#endif
  }
  std::cerr << "opt: time=" << timelimitTimer << std::endl;
  return 0;

}

