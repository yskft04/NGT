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

#include <fstream>
#include <string>
#include <cstddef>
#include <stdint.h>
#include <iostream>
#include <stdexcept>
#include <cerrno>
#include <cstring>

namespace NGT {
  class ObjectSpace;
};

template <class TYPE>
class ObjectFile {
 private:
  enum Type {
    TypeFloat	= 0,
    TypeUint8	= 1,
    TypeInt8	= 2
  };
  struct FileHeadStruct {
    uint32_t noOfObjects;
    uint32_t noOfDimensions;
  };

  bool			_isOpen;
  std::ifstream		_stream;
  FileHeadStruct	_fileHead;

  uint32_t		_recordSize;
  uint32_t		_sizeOfElement;
  std::string		_typeName;
  Type			_type;
  std::string		_objectPath;
  std::string		_fileName;

  bool			_readFileHead();

  size_t		_pseudoDimension;
  std::vector<ObjectFile<TYPE>*>	_objectFiles;

 public:
  ObjectFile();
  ~ObjectFile();
  bool create(const std::string &file, const std::string &objectPath);
  bool open(const std::string &file, const size_t pseudoDimension = 0);
  void close();
  size_t insert(TYPE &data, NGT::ObjectSpace *objectSpace = 0) {std::cerr << "insert: not implemented."; abort();}
  void put(const size_t id, TYPE &data, NGT::ObjectSpace *objectSpace = 0) {std::cerr << "put: not implemented."; abort();}
  bool get(size_t id, TYPE &data, NGT::ObjectSpace *objectSpace = 0);
  bool get(const size_t streamID, size_t id, TYPE &data, NGT::ObjectSpace *objectSpace = 0);
  void remove(const size_t id) {std::cerr << "remove: not implemented."; abort();}
  bool isOpen() const;
  size_t size();
  size_t getRecordSize() { return _recordSize; }
  bool multipleOpen(const size_t nOfStreams);
  void multipleClose();
};


// constructor 
template <class TYPE>
ObjectFile<TYPE>::ObjectFile()
  : _isOpen(false) {
}

// destructor
template <class TYPE>
ObjectFile<TYPE>::~ObjectFile() {
  close();
}

template <class TYPE>
bool ObjectFile<TYPE>::create(const std::string &file, const std::string &objectPath) {
  std::cerr << "ObjectFile create " << file << ":" << objectPath << std::endl;
  {
    _stream.open(objectPath, std::ios::in);
    if (!_stream) {
      std::cerr << "Cannot open " << objectPath << std::endl;
      return false;
    }
    bool ret = _readFileHead();
    if (ret == false) {
      return false;
    }
    _stream.close();
  }
  {
    std::ofstream tmpstream;
    tmpstream.open(file);
    std::vector<std::string> tokens;
    NGT::Common::tokenize(objectPath, tokens, ".");
    tmpstream << _fileHead.noOfObjects << std::endl;
    tmpstream << _fileHead.noOfDimensions << std::endl;
    if (tokens.size() <= 1) {
      std::cerr << "The specifiled object file name has no proper extensions. " << objectPath << " use fbin." << std::endl;
      tmpstream << "fbin" << std::endl;
    } else {
      tmpstream << tokens.back() << std::endl;
    }
    tmpstream << objectPath << std::endl;
  }
  return true;
}

template <class TYPE>
bool ObjectFile<TYPE>::open(const std::string &file, size_t pseudoDimension) {
  _pseudoDimension = pseudoDimension;
  uint32_t noOfObjects;
  uint32_t noOfDimensions;
  _fileName = file;
  {
    std::ifstream tmpstream;
    tmpstream.open(file);
    if (!tmpstream) {
      std::cerr << "Cannot open " << file << std::endl;
      abort();
    }
    tmpstream >> noOfObjects;
    tmpstream >> noOfDimensions;
    tmpstream >> _typeName;
    tmpstream >> _objectPath;
  }

  if (_typeName == "fbin") {
    _sizeOfElement = 4;
    _type = TypeFloat;
  } else if (_typeName == "u8bin") {
    _sizeOfElement = 1;
    _type = TypeUint8;
  } else if (_typeName == "i8bin") {
    _sizeOfElement = 1;
    _type = TypeInt8;
  } else {
    _sizeOfElement = 4;
    _type = TypeFloat;
  }

  _stream.open(_objectPath, std::ios::in);
  if(!_stream){
    _isOpen = false;    
    return false;
  }
  _isOpen = true;

  bool ret = _readFileHead();
  if (_fileHead.noOfObjects != noOfObjects) {
    std::cerr << "Invalid # of objects=" << _fileHead.noOfObjects << ":" << noOfObjects << std::endl;
    abort();
  }
  if (_fileHead.noOfDimensions != noOfDimensions) {
    std::cerr << "Invalid # of dimensions=" << _fileHead.noOfDimensions << ":" << noOfDimensions << std::endl;
    abort();
  }
  _recordSize = _sizeOfElement * _fileHead.noOfDimensions;
  return ret;
}

template <class TYPE>
bool ObjectFile<TYPE>::multipleOpen(const size_t nOfStreams) {
  if (!isOpen()) {
    return false;
  }
  if (!_objectFiles.empty()) {
    std::cerr << "ObjectFile : already opened multiple streams. close and reopen. # of streams=" << nOfStreams << std::endl;
    multipleClose();
  }
  for (size_t i = 0; i < nOfStreams; i++) {
    auto *of = new ObjectFile<TYPE>;
    if (!of->open(_fileName, _pseudoDimension)) {
      return false;
    }
    _objectFiles.push_back(of);
  }
  return true;
}

template <class TYPE>
bool ObjectFile<TYPE>::get(const size_t streamID, size_t id, TYPE &data, NGT::ObjectSpace *objectSpace) {
  if (streamID >= _objectFiles.size()) {
    std::cerr << "streamID is invalid. " << streamID << ":" << _objectFiles.size() << std::endl;
    return false;
  }
  if (!_objectFiles[streamID]->get(id, data, objectSpace)) {
    return false;
  }
  return true;
}

template <class TYPE>
void ObjectFile<TYPE>::close(){
  _stream.close();
  _isOpen = false;
  multipleClose();
}

template <class TYPE>
void ObjectFile<TYPE>::multipleClose() {
  for (auto *i : _objectFiles) {
    i->close();
    delete i;
    i = 0;
  }
  _objectFiles.clear();
}

template <class TYPE>
bool ObjectFile<TYPE>::get(size_t id, TYPE &data, NGT::ObjectSpace *objectSpace) {
  id--;
  if( size() <= id ){
    return false;
  }
  
  //uint64_t offset_pos = (id * (sizeof(RecordStruct) + _fileHead.recordSize)) + sizeof(FileHeadStruct);
  uint64_t offset_pos = id * _recordSize + sizeof(FileHeadStruct);
  //offset_pos += sizeof(RecordStruct);  
  _stream.seekg(offset_pos, std::ios::beg);
  if (!_stream.fail()) {
    switch (_type) {
    case TypeFloat:
      if (_pseudoDimension == 0) {
	data.deserialize(_stream, objectSpace);
      } else {
	auto dim = _pseudoDimension;
	if (dim == 0) {
	  dim = _recordSize / sizeof(float);
	}
	if (_recordSize > dim * sizeof(float)) {
	  abort();
	}
	std::vector<float> record(dim);
	_stream.read(reinterpret_cast<char*>(record.data()), _recordSize);
	std::stringstream object;
	object.write(reinterpret_cast<char*>(record.data()), dim * sizeof(float));
	data.deserialize(object, objectSpace);
      }
      break;
    case TypeUint8:
    case TypeInt8:
      {
	auto dim = _pseudoDimension;
	if (dim == 0) {
	  dim = _recordSize;
	}
	uint8_t src[_recordSize];
	_stream.read(reinterpret_cast<char*>(src), _recordSize);
	std::vector<float> dst(dim);
	if (_type == TypeUint8) {
	  for (size_t i = 0; i < _recordSize; i++) {
	    dst[i] = static_cast<float>(src[i]);
	  }
	} else {
	  int8_t *intsrc = reinterpret_cast<int8_t*>(&src[0]);
	  for (size_t i = 0; i < _recordSize; i++) {
	    dst[i] = static_cast<float>(intsrc[i]);
	  }
	}
	std::stringstream object;
	object.write(reinterpret_cast<char*>(dst.data()), dim * sizeof(float));
	//delete[] dst;
	data.deserialize(object, objectSpace);
      }
      break;
    }
  } else {
    std::cerr << "ObjectFile::get something wrong! id=" << id << " type=" << _type << std::endl;
    abort();
#if 0

    const int trialCount = 10;
    for (int tc = 0; tc < trialCount; tc++) {
      _stream.clear();
      _stream.seekg(offset_pos, std::ios::beg);
      if (_stream.fail()) {
	continue;
      }
      data.deserialize(_stream, objectSpace);
      if (_stream.fail()) {
	continue;
      } else {
	break;
      }
    }
    if (_stream.fail()) {
      throw std::runtime_error("ObjectFile::get: Error!");
    }
#endif
  }

  return true;
}

template <class TYPE>
bool ObjectFile<TYPE>::isOpen() const
{
  return _isOpen;
}

template <class TYPE>
size_t ObjectFile<TYPE>::size()
{
  _stream.seekg(0, std::ios::end);
  int64_t offset_pos = _stream.tellg();
  offset_pos -= sizeof(FileHeadStruct);
  size_t num = offset_pos / _recordSize;
  num++; 
  return num;
}

template <class TYPE>
bool ObjectFile<TYPE>::_readFileHead() {
  _stream.seekg(0, std::ios::beg);
  _stream.read((char *)(&_fileHead), sizeof(FileHeadStruct));
  if(_stream.bad()){
    return false;
  }
  return true;
}

