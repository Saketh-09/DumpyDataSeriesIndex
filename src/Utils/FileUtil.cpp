#include "../../include/Utils/FileUtil.h"
#include "../../include/DataStructures/TimeSeries.h"
#include <fstream>
#include <cstdio>
#include <filesystem>
#include <dirent.h>
#include <cstdlib>


bool FileUtil::checkFileExists(const std::string& name) {
    std::ifstream f(name);
    return f.good();
}

bool FileUtil::fileRemove(const std::string& fname) {
    return std::remove(fname.c_str());
}

void FileUtil::deleteFiles(const std::vector<std::string>& files) {
    for(const auto& file : files) {
        fileRemove(file);
    }
}

long FileUtil::getFileSize(const std::string& fname) {
    std::filesystem::path filePath(fname);
    if(std::filesystem::exists(filePath)) {
        return std::filesystem::file_size(filePath);
    }
    return -1;
}

std::string FileUtil::getFilePath(const std::string& path, const std::string& filename) {
    std::filesystem::path filePath = std::filesystem::path(path) / filename;
    return filePath.string();
}

bool FileUtil::createDir(const std::string& path) {
    return std::filesystem::create_directory(path);
}

bool FileUtil::checkDirClean(const std::string& path) {
    if (!checkFileExists(path)) {
        createDir(path);
    }
    for (const auto& entry : std::filesystem::directory_iterator(path)) {
        if (entry.is_regular_file()) {
            std::filesystem::remove(entry.path());
        } else if (entry.is_directory()) {
            checkDirClean(entry.path().string());
            std::filesystem::remove(entry.path());
        }
    }
    return false;
}

void FileUtil::mergeFiles(const string sources[], const string& dest, int num) {
    FILE * out = fopen(dest.c_str(), "wb");
    int i =0;
    while( i < num) {
        long size = getFileSize(sources[i].c_str());
        char * tmp = new char[size];
        FILE *in = fopen(sources[i].c_str(), "rb");
        fread(tmp, 1, size, in);
        fclose(in);
        fwrite(tmp, 1, size, out);
        free(tmp);
        i++;
    }
    fclose(out);
}

float* FileUtil::readQueries() {
    auto queries = new float[Const::tsLength * Const::query_num];
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    int i = 0;
    while (i < Const::query_num) {
        fread(queries + i * Const::tsLength, sizeof(float ), Const::tsLength, f);
        i++;
    }
    fclose(f);
    return queries;
}

float* FileUtil::readSeries(std::FILE* f) {
    auto *ts = new float[Const::tsLength];
    std::fread(ts, sizeof(float), Const::tsLength, f);
    return ts;;
}

void FileUtil::getFiles(const std::string& path, std::vector<std::string>& files) {
    files.clear();
    for (const auto& entry : std::filesystem::directory_iterator(path)) {
        if (entry.is_regular_file()) {
            files.push_back(entry.path().string());
        }
    }
}