#!/bin/bash

elapse_time() {
    END_TIME=`date +%s`
    ES=$(($END_TIME - $START_TIME))
    ((S=ES%60, M=(ES/60)%60, H=(ES/3600)%60))
    echo [ELAPSED TIME] H:M:S=$H:$M:$S
    ELAPSE_TIME="$H:$M:$S"
}

execute() {
    if [ -z $2 ]; then
	echo "execute: $1 &>> $LOG"
	`$1 &>> $LOG`
    else
	echo "execute: $1 &> ${LOG_ROOT}/$2"
	`$1 &> ${LOG_ROOT}/$2`
    fi
}

execute_search() {
    if [ -z $2 ]; then
	echo "execute_search: $1 >> $LOG"
	`$1 &>> $LOG`
    else
	echo "execute_search: $1 > $2"
	`$1 &> $2`
    fi
}


extract_objects() {
    echo [start] extract_objects: `date`
    START_TIME=`date +%s`
    OBJECTS=${WORKSPACE_ROOT}/base.$1.u8.tsv
    if [ -n "$BINOBJECTS" ]; then
	echo "extract objects $BINOBJECTS -> $OBJECTS for benchmarks"
	$NEURIPS21 extract -n $1 -d $DIMENSION $BINOBJECTS > $OBJECTS
    fi
    echo "objects $OBJECTS"
    if [ ! -f $OBJECTS ]; then
	if [ ! -f ${DATASET_BASE} ]; then
	    echo "Not found the object file. ${DATASET_BASE}"
	    exit
	fi
	echo "extract objects ${DATASET_BASE} -> $OBJECTS"
	echo "$NEURIPS21 extract -n $1 -d $DIMENSION ${DATASET_BASE} > $OBJECTS"
	$NEURIPS21 extract -n $1 -d $DIMENSION ${DATASET_BASE} > $OBJECTS
    fi
    RETURN=$OBJECTS
    elapse_time
    echo [end] extract_objects: `date`
    EO_ELAPSE_TIME=$ELAPSE_TIME
}

generate_blobs_with_hierarchical_clustering() {
    echo [start] generate_blobs_with_hierarchical_clustering: `date`
    START_TIME=`date +%s`
    CLUSTER=${CLUSTER_ROOT}/hierarchy-cluster-$NSTR-$BLOB
    CINDEX=${CLUSTER}_index.tsv
    CLUSTER_CENTROID=${CLUSTER}_centroid.tsv
    if [ ! -f ${CLUSTER_CENTROID} ]; then
	echo "clustering CLUSTER=$CLUSTER"
	echo "clustering CINDEX=$CINDEX"
	TMP_INDEX="index.$$"
	rm -rf $TMP_INDEX
	echo "$NGTLQG create -d $DIMENSION -o f -D 2 -C $BLOB -c 16 -N $NUM_OF_SUBVECTORS -M s -L k -B 1 $TMP_INDEX ${DATASET_BASE}"
	$NGTLQG create -d $DIMENSION -o f -D 2 -C $BLOB -c 16 -N $NUM_OF_SUBVECTORS -M s -L k -B 1 $TMP_INDEX ${DATASET_BASE}
	MAX=`expr $N / $BLOB \* 2`
	echo "The max number of objects in each clusters=$MAX"
	COM="$NGTLQG kmeans -n $N -m $MAX $TMP_INDEX ${CLUSTER}"
	COM_LOG="hcluster.log"
	execute "$COM" $COM_LOG
	rm -rf $TMP_INDEX
    fi
    elapse_time
    echo [end] generate_blobs_with_hierarchical_clustering: `date`
    GBHC_ELAPSE_TIME=$ELAPSE_TIME
}


generate_blobs_and_quantization_clusters_with_kmeans() {
    echo [start] generate_blobs_and_quantization_clusters_with_kmeans: `date`
    START_TIME=`date +%s`
    CLUSTER=${CLUSTER_ROOT}/kmeans-cluster-$NSTR-$BLOB-$QUANTIZATION
    CINDEX=${CLUSTER}_index.tsv
    CLUSTER_CENTROID=${CLUSTER}_centroid.tsv
    if [ ! -f ${CLUSTER_CENTROID} ]; then
	echo "${CLUSTER_CENTROID} is not generated yet."
	echo "clustering CLUSTER=$CLUSTER"
	echo "clustering CINDEX=$CINDEX"
	TMP_INDEX="index.$$"
	rm -rf $TMP_INDEX
	echo "$NGTLQG create -d $DIMENSION -o f -D 2 -C $BLOB -c 16 -N $NUM_OF_SUBVECTORS -M s -L k -B 1 $TMP_INDEX ${DATASET_BASE}"
	$NGTLQG create -d $DIMENSION -o f -D 2 -C $BLOB -c 16 -N $NUM_OF_SUBVECTORS -M s -L k -B 1 $TMP_INDEX ${DATASET_BASE}
	MAX=`expr $N \* $NUM_OF_CLUSTERS \* $NUM_OF_CLUSTERS_FOR_ADJUST / $QUANTIZATION`
	echo "The max number of objects in each cluster=$MAX"
	COM_LOG="clustering.log"
	FC=`echo "sqrt($QUANTIZATION)"|bc`
	FS=`expr $FC \* 2000`
	QS=`expr $QUANTIZATION \* $QUANTIZATION_SAMPLE_EXPANSION`
	COM="$NGTLQG kmeans -n $N -m $MAX -N $NUM_OF_CLUSTERS -B $FC:$FS,$QUANTIZATION:$QS,$BLOB $TMP_INDEX ${CLUSTER}"
	execute "$COM" $COM_LOG
	rm -rf $TMP_INDEX
    fi
    elapse_time
    echo [end] generate_blobs_quantization_clusters_with_kmeans: `date`
    GBHC_ELAPSE_TIME=$ELAPSE_TIME
}


generate_R_and_local_codebooks() {
    echo [start] generate_R_and_local_codebooks: `date`
    START_TIME=`date +%s`
    ITER=$R_ITERATION
    CITER=100
    SEED="h"
    CMODE="k"
    LC=${CLUSTER_RESIDUAL%.*}_opt
    OPTR="${LC}_R.tsv"
    echo "** Generate R and local codebooks. ${CLUSTER_RESIDUAL} -> $LC"
    if [ ! -f ${LC}-0_centroid.tsv ]; then
	RESIDUAL_TMP_FILE="residual_file.$$"
	echo "Optimization! ${R_SAMPLE} ${RESIDUAL_TMP_FILE}"
	if [ -n "$R_SAMPLE" -a $R_SAMPLE -ne 0 ]; then
	    head -${R_SAMPLE} ${CLUSTER_RESIDUAL} > ${RESIDUAL_TMP_FILE}
	    echo "Cut for the optimization. ${R_SAMPLE}, `wc -l ${RESIDUAL_TMP_FILE}`"
	else
	    ln -s ${CLUSTER_RESIDUAL} ${RESIDUAL_TMP_FILE}
	fi
	COM="$OPQ -n 16 -m ${NUM_OF_SUBVECTORS} -I $CITER -i $SEED -O t -s f -C $CMODE -t $ITER ${RESIDUAL_TMP_FILE} ${CLUSTER_RESIDUAL%.*}_opt.tsv"
	execute "$COM"
	rm -f ${RESIDUAL_TMP_FILE}
	CUTLAST=`expr $DIMENSION - ${SUBVECTOR_DIMENSION} + 1`
	for S in `seq 1 ${SUBVECTOR_DIMENSION} $CUTLAST`
	do
	    IDX=`expr \( $S - 1 \) / 2`
	    mv ${CLUSTER_RESIDUAL%.*}_opt-$IDX.tsv ${CLUSTER_RESIDUAL%.*}_opt-${IDX}_centroid.tsv
	done
    fi
    LOCAL_CREATION_MODE=s
    elapse_time
    echo [end] generate_R_and_local_codebooks: `date`
    GRLC_ELAPSE_TIME=$ELAPSE_TIME
}


generate_R_and_local_codebooks_from_kmeans() {
    echo [start] generate_R_and_local_codebooks_from_kmeans: `date`
    START_TIME=`date +%s`
    ITER=$R_ITERATION
    CITER=500
    SEED="h"
    CMODE="k"
    echo "** Generate R and local codebooks from kmeans. ${CLUSTER} -> $OPTR"
    OPT_DIR=${CLUSTER}_opt
    if [ ! -d $OPT_DIR ]; then
	mkdir $OPT_DIR
    fi
    OPTR="${OPT_DIR}/opt_R.tsv"
    if [ ! -f ${OPT_DIR}/opt-0_centroid.tsv ]; then
	COM="$OPQ -L 24 -M $R_N -S 0.02 -X $R_STEP -R 0.9 -e r -H ${R_SAMPLE} -n 16 -m ${NUM_OF_SUBVECTORS} -I $CITER -i $SEED -O t -s f -C $CMODE -t $ITER ${RANDOM_OBJECT} ${OPT_DIR}/opt.tsv $GC"
	execute "$COM" opq.log
	CUTLAST=`expr ${NUM_OF_SUBVECTORS} - 1`
	for IDX in `seq 0 $CUTLAST`
	do
	    mv ${OPT_DIR}/opt-$IDX.tsv ${OPT_DIR}/opt-${IDX}_centroid.tsv
	done
    fi
    LOCAL_CREATION_MODE=s
    elapse_time
    echo [end] generate_R_and_local_codebooks_from_kmeans: `date`
    GRLC_ELAPSE_TIME=$ELAPSE_TIME
}

write_settings() {
    SETTINGS=$INDEX/SETTINGS
    echo "TITLE=$TITLE" > $SETTINGS
    echo "TYPE=$TYPE" >> $SETTINGS
    echo "HOST=$HOSTNAME" >> $SETTINGS
}

build() {
    echo [start] build: `date`
    START_TIME=`date +%s`
    if [ ! -d $INDEX_DIR ]; then
	mkdir $INDEX_DIR
    fi
    if [ ! -e $INDEX/ivt ]; then
	echo "create index ($INDEX) ..."
	echo "$NGTLQG create -d $DIMENSION -o f -D 2 -E 10 -C $QUANTIZATION -c 16 -N $NUM_OF_SUBVECTORS -M $CREATION_MODE -L ${LOCAL_CREATION_MODE} -B 1 $INDEX ${DATASET_BASE}"
	if [ -n "$OPTR" ]; then
	    $NGTLQG create -d $DIMENSION -o f -D 2 -E 10 -C $QUANTIZATION -c 16 -N $NUM_OF_SUBVECTORS -M $CREATION_MODE -L ${LOCAL_CREATION_MODE} -B 1 -s 800 $INDEX ${DATASET_BASE} ${OPTR}
	else
	    $NGTLQG create -d $DIMENSION -o f -D 2 -E 10 -C $QUANTIZATION -c 16 -N $NUM_OF_SUBVECTORS -M $CREATION_MODE -L ${LOCAL_CREATION_MODE} -B 1 -s 800 $INDEX ${DATASET_BASE}
	fi
	echo "$NGT append ${INDEX}/global/ $GC"
	$NGT append ${INDEX}/global/ $GC
	if [ -n "$LC" ]; then
	    echo "set local centroids $LC"
	    CUTLAST=`expr ${NUM_OF_SUBVECTORS} - 1`
	    for IDX in `seq 0 $CUTLAST`
	    do
		echo "append $S : $IDX"
		echo "$NGT append ${INDEX}/local-$IDX ${LC}-${IDX}_centroid.tsv"
		$NGT append ${INDEX}/local-$IDX ${LC}-${IDX}_centroid.tsv
	    done
	fi
	echo "Build: Quantization centroids=QCC:$QCC"
	echo "Build: Blob to quantization index=QCI::$QCI"
	echo "Build: Object to blob index=QBI::$OBI"
	COM="$NGTLQG build -m q -n $N $INDEX $QCC $QCI $OBI"
	COM_LOG="build_q.log"
	execute "$COM" $COM_LOG
    fi
    if [ ! -e $INDEX/grp ]; then
	COM="$NGTLQG build -m g -n $N $INDEX"
	COM_LOG="build_g.log"
	execute "$COM" $COM_LOG
    else
	echo "Found $INDEX"
    fi

    if [ ! -f $INDEX/onng ]; then
	if [ $BLOB -ge 500000 ]; then
	    COM="$NGT reconstruct-graph -m S -s p -o 10 -i 120 $INDEX/global $INDEX/onng"
	    COM_LOG="reconstruct.log"
	    execute "$COM" $COM_LOG
	    mv $INDEX/global $INDEX/anng
	    ln -s onng $INDEX/global
	fi
    fi

    elapse_time
    write_settings
    echo [end] build: `date`
    B_ELAPSE_TIME=$ELAPSE_TIME
}

generate_data() {
    echo "Create ..."
    case $MODE in
	HIERARCHY_STATIC_OPT)
	    echo "$MODE ${GC%.*}"
	    TYPE=$MODE-$NSTR-x-x-$BLOB-$QUANTIZATION-$SUBVECTOR_DIMENSION
	    TITLE="$MODE	$NSTR	$SAMPLE	$QUANTIZATION	$BLOB_SAMPLE	$BLOB	$SUBVECTOR_DIMENSION	$R_SAMPLE	$R_ITERATION"
	    if [ -n "$SEARCH" ]; then
		return
	    fi
	    
	    generate_blobs_with_hierarchical_clustering
	    GC=${CLUSTER}_centroid.tsv
	    echo "**GC=$GC"

	    generate_quantization_clusters_from_blobs

	    generate_R_and_local_codebooks

	    echo "** CLUSTER=${CLUSTER} CLUSTER_CENTROID=${CLUSTER_CENTROID} CLUSTER_DIR=${CLUSTER_DIR}"
	    QCC=${CLUSTER_CENTROID}
	    QCI=${CLUSTER}_index.tsv
	    CREATION_MODE=l
	    ;;
	KMEANS_QUANTIZATION)
	    echo "$MODE ${GC%.*}"
	    TYPE=$DATASET-$MODE-$NSTR-$BLOB-$QUANTIZATION-$SUBVECTOR_DIMENSION-$R_SAMPLE-$R_ITERATION-$QUANTIZATION_SAMPLE_EXPANSION
	    TITLE="$MODE	$NSTR	$QUANTIZATION	$BLOB	$SUBVECTOR_DIMENSION	$R_SAMPLE	$R_ITERATION	$QUANTIZATION_SAMPLE_EXPANSION"
	    if [ -n "$SEARCH" ]; then
		return
	    fi
	    
	    NUM_OF_CLUSTERS_FOR_ADJUST=2

	    generate_blobs_and_quantization_clusters_with_kmeans

	    GC=${CLUSTER}_qcentroid.tsv
	    echo "GC=${CLUSTER}_qcentroid.tsv"
	    echo "$GC"
	    

	    extract_objects ${R_SAMPLE}
	    RANDOM_OBJECT=$RETURN
	    generate_R_and_local_codebooks_from_kmeans

	    GC=${CLUSTER}_centroid.tsv
	    echo "**GC=$GC"
	    OBI=${CLUSTER}_index.tsv
	    echo "**OBI=$OBI"
	    QCC=${CLUSTER}_qcentroid.tsv
	    echo "**QCC=$QCC"
	    QCI=${CLUSTER}_bqindex.tsv
	    echo "**QCI=$QCI"
	    LC=${OPT_DIR}/opt
	    CREATION_MODE=l
	    ;;
	*)
	    echo "mode error! $MODE"
	    exit 1
    esac
}

echo Start: `date`
N=1000000
BLOB=10000
QUANTIZATION=1024
SUBVECTOR_DIMENSION=1
MODE=HEAD
BLOB_SAMPLE=0
SAMPLE=0
DATASET=dataset
R_ITERATION=20
R_N=5
BINOBJECTS=
NUM_OF_CLUSTERS=10
NUM_OF_RANDOM_OBJECTS=2
RANDOM_OBJECTS_MODE=c
QUANTIZATION_SAMPLE_EXPANSION=50

NGT_ROOT=.
NGTLQG=./ngtlqg
NGTC=./ngtc
NGT=./ngt
OPQ=./opq
NEURIPS21=./neurips21

FORCE=false
while [ $# -gt 0 ]; do
    case "$1" in
	--root=*)
	    NGT_ROOT="${1#*=}"
	    ;;
	--dataset=*)
	    DATASET="${1#*=}"
	    ;;
	--#-of-objects=* | -n=*)
	    N="${1#*=}"
	    ;;
	--#-of-blobs=* | -b=*)
	    BLOB="${1#*=}"
	    ;;
	--#-of-clusters=* | -C=*)
	    NUM_OF_CLUSTERS="${1#*=}"
	    ;;
	--#-of-random-objects=* | -M=*)
	    NUM_OF_RANDOM_OBJECTS="${1#*=}"
	    ;;
	--#-of-centroids=* | -q=*)
	    QUANTIZATION="${1#*=}"
	    ;;
	--quantization-sample-expansion=* | -E=*)
	    QUANTIZATION_SAMPLE_EXPANSION="${1#*=}"
	    ;;
	--#-of-dimension=* | -d=*)
	    DIMENSION="${1#*=}"
	    ;;
	--#-of-suvbector-dimension=* | -D=*)
	    SUBVECTOR_DIMENSION="${1#*=}"
	    ;;
	--#-of-samples=* | -s=*)
	    SAMPLE="${1#*=}"
	    ;;
	--#-of-blob-samples=* | -B=*)
	    BLOB_SAMPLE="${1#*=}"
	    ;;
	--#-of-r-matrices=* | -N=*)
	    R_N="${1#*=}"
	    ;;
	--#-of-r-samples=* | -r=*)
	    R_SAMPLE="${1#*=}"
	    ;;
	--r-iteration=* | -R=*)
	    R_ITERATION="${1#*=}"
	    ;;
	--object=*)
	    BINOBJECTS="${1#*=}"
	    ;;
	--base=*)
	    BASE="${1#*=}"
	    ;;
	--query=*)
	    QUERY="${1#*=}"
	    ;;
	--benchmark=*)
	    BENCHMARK="${1#*=}"
	    ;;
	-X=*)
	    R_STEP="${1#*=}"
	    ;;
	--search)
	    SEARCH=true
	    ;;
	-m=*)
	    MODE="${1#*=}"
	    ;;
	-f)
	    FORCE=true
	    ;;
	*)
	    echo "Error! invalid argument. $1"
	    echo "  Usage: eval.sh --n #-of-objects --blob #-of-blobs --q #-of-centroids"
	    exit 1
    esac
    shift
done

if [ -n "$BENCHMARK" ]; then
    echo "Benchmark settings. [$BENCHMARK]"
    if [ -z $BINOBJECTS ]; then
	echo "BINOBJECTS should be specified"
	exit 1
    fi
    DATASET_BASE=$BINOBJECTS
    NGTLQG=ngtlqg
    NGTC=ngtc
    NGT=ngt
    OPQ=opq
    NEURIPS21=neurips21
    if [ -z $NGT_ROOT ]; then
	NGT_ROOT=ngt
    fi
    if [ ! -d $NGT_ROOT ]; then
	echo "No directory! $NGT_ROOT"
	exit 1
    fi
else
    DATASET_BASE=${DATASET}/${BASE}
fi
LOG_ROOT=${NGT_ROOT}/logs
RESULT_ROOT=${NGT_ROOT}/results
WORKSPACE_ROOT=${NGT_ROOT}/workspace
CLUSTER_ROOT=$WORKSPACE_ROOT

echo "BINOBJECTS=$BINOBJECTS"
if [ -z "$BINOBJECTS" -a ! -d $DATASET ]; then
    echo "No dataset. Make a symbolic link to the dataset. $DATASET"
    exit 1
fi

if [ `expr $N / 1000000` -ge 1000 ]; then
  NSTR=`expr $N / 1000000000`B
  SSTR=`expr $SAMPLE / 1000000000`B
else
  NSTR=`expr $N / 1000000`M
  SSTR=`expr $SAMPLE / 1000000`M
fi
CREATION_MODE=s
LOCAL_CREATION_MODE=k

NUM_OF_SUBVECTORS=`expr $DIMENSION / $SUBVECTOR_DIMENSION`
echo "N=$N ($NSTR)"
echo "BLOB=$BLOB"
echo "QUANTIZATION=$QUANTIZATION"
echo "SUBVECTOR_DIMENSION=$SUBVECTOR_DIMENSION"
echo "NUM_OF_SUBVECTORS=$NUM_OF_SUBVECTORS"

if [ -n "$BENCHMARK" ]; then
    if [ -d $BENCHMARK ]; then
	echo "index already exists. $INDEX"
	exit 0
    fi
fi

if [ ! -d $NGT_ROOT ]; then
    echo "Create the NGT root directory. $NGT_ROOT"
    mkdir $NGT_ROOT
fi
if [ ! -d $WORKSPACE_ROOT ]; then
    echo "Create the workspace directory. $WORKSPACE_ROOT"
    mkdir $WORKSPACE_ROOT
fi
if [ ! -d $LOG_ROOT ]; then
    echo "Create the log directory. $LOG_ROOT"
    mkdir $LOG_ROOT
fi
if [ ! -d $RESULT_ROOT ]; then
    echo "Create the result directory. $RESULT_ROOT"
    mkdir $RESULT_ROOT
fi

LOG=${LOG_ROOT}/log.$$
rm -f $LOG

generate_data

INDEX_DIR=${NGT_ROOT}/indexes
INDEX_TITLE=index-$TYPE
if [ -n "$BENCHMARK" ]; then
    INDEX=$BENCHMARK
else
    INDEX=$INDEX_DIR/$INDEX_TITLE
fi

echo "INDEX=$INDEX"

if [ -z "$SEARCH" ]; then
    build
fi

if [ -n "$BENCHMARK" ]; then
    echo "bench mark index=$BENCHMARK"
    echo "NGT index=$INDEX_DIR/$INDEX_TITLE"
    exit 0
fi

rm -f NGTLQG_ERROR

SEARCH_RESULT=${RESULT_ROOT}/search.txt
echo "search index=$INDEX result=$SEARCH_RESULT"
ETIME="EO=$EO_ELAPSE_TIME	GQC=$GQC_ELAPSE_TIME	GBQC=$GBQC_ELAPSE_TIME	GBHC=$GBHC_ELAPSE_TIME	GQCLB=$GQCLB_ELAPSE_TIME	GRLC=$GRLC_ELAPSE_TIME	B=$B_ELAPSE_TIME"
echo "@ELAPSE@ $ETIME"
if [ -n "$BUILD_LOG" ]; then
    if [ -f $BUILD_LOG ]; then
	echo "@QUANTIZATION ERROR@ `grep Distance\(error\)= $BUILD_LOG`" 
    fi
fi
HOST=`hostname --short`
TITLE="@EVAL@	$HOST	$TITLE"
CUTBACK=0.0
NUM_OF_EDGES=100
#MAX_NUM_OF_BLOBS=100
for MAX_NUM_OF_BLOBS in 10000
do
    for EPSILON in 0.0 0.06 0.1 0.12 0.14 0.16 0.18 0.20 0.25 0.30 0.35 0.40 0.50 0.60
    do
	COM="$NGTLQG search -M g -o e -n 10 -N $MAX_NUM_OF_BLOBS -B 0.0 -e $EPSILON -E ${NUM_OF_EDGES} -C $CUTBACK $INDEX ${DATASET}/${QUERY}.tsv"
	execute_search "$COM" $SEARCH_RESULT
	$NGT eval -m a -t "$TITLE	${NUM_OF_EDGES}	$MAX_NUM_OF_BLOBS	$CUTBACK" ${DATASET}/${QUERY}-$NSTR.gt $SEARCH_RESULT
    done
done

if [ -f NGTLQG_ERROR ]; then
    cat NGTLQG_ERROR 
fi

echo End: `date`
