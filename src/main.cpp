#include <iostream>
#include <cstdio>
#include <vector>
#include <thread>
#include <chrono>
#include "../include/DataStructures/Node.h"
#include "../include/DataStructures/GraphConstruction.h"
#include "../include/Searchers/Searcher.h"
#include "../include/Utils/FileUtil.h"
#include "../include/Utils/MathUtil.h"
#include "../include/Utils/TimeSeriesUtil.h"
using namespace std;

void constructGraph(){
    GraphConstruction::buildAndSave2Disk();
}

vector<vector<int>>* loadGraphSkeleton(){
    int vd = 0;
    int i=1;
    while(i<=Const::bitsReserve) {
        vd += MathUtil::nChooseK(Const::segmentNum, i);
        i++;
    }
    auto nnList = new vector<vector<int>>(Const::vertexNum, vector<int>(vd, -1));

    if(!FileUtil::checkFileExists(Const::graphfn.c_str())){
        cout << "File not exists!" << Const::graphfn << endl;
        exit(-1);
    }
    FILE *f = fopen(Const::graphfn.c_str(), "rb");
    i=1;
    while(i<Const::vertexNum) {
        fread(&((*nnList)[i][0]), sizeof(int), vd, f);
        i++;
    }
    return nnList;

}

void buildDumpy(){
    auto g = loadGraphSkeleton();
    Node* root = Node::BuildingIndex(Const::datafn, Const::saxfn);
    root->save2Disk(Const::idxfn + "root.idx");
}

void approxSearchOneNode() {
    Node *root = Node::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *approxKnn = Searcher::approxSearch(root, queries + i * Const::tsLength, Const::k,
                                                                        g, Const::idxfn);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j) {
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
        }
    }
}

void approxSearchMoreNode() {
    Node *root = Node::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *approxKnn = Searcher::approxIncSearch(root, queries + i * Const::tsLength,
                                                                           Const::k, Const::idxfn,
                                                                           Const::visited_node_num);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
    }
}

void approxSearchOneNodeDTW() {
    Node *root = Node::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *approxKnn = Searcher::approxSearchDTW(root, queries + i * Const::tsLength, Const::k,
                                                                        g, Const::idxfn);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j) {
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
        }
    }
}

void approxSearchMoreNodeDTW() {
    Node *root = Node::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *approxKnn = Searcher::approxIncSearchDTW(root, queries + i * Const::tsLength,
                                                                           Const::k, Const::idxfn,
                                                                           Const::visited_node_num);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
    }
}

void buildFuzzy(){
    auto g = loadGraphSkeleton();
    int bound = Const::fuzzy_f * 100;
    Const::fuzzyidxfn += "/" + to_string(bound) + "-" + to_string(Const::delta) + "/";
    Node* root = Node::BuildIndexFuzzy(Const::datafn, Const::saxfn, Const::paafn, g);
    root->save2Disk(Const::fuzzyidxfn + "root.idx");
}

void approxSearchOneNodeFuzzy() {
    int bound = Const::fuzzy_f * 100;
    Const::fuzzyidxfn += "/" + to_string(bound) + "-" + to_string(Const::delta) + "/";
    Node *root = Node::loadFromDisk(Const::saxfn, Const::fuzzyidxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *approxKnn = Searcher::approxSearch(root, queries + i * Const::tsLength, Const::k,
                                                                        g, Const::fuzzyidxfn);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
    }
}

void approxSearchMoreNodeFuzzy() {
    int bound = Const::fuzzy_f * 100;
    Const::fuzzyidxfn += "/" + to_string(bound) + "-" + to_string(Const::delta) + "/";
    Node *root = Node::loadFromDisk(Const::saxfn, Const::fuzzyidxfn + "root.idx", false);
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        auto start = chrono::system_clock::now();
        vector<PqItemSeries *> *approxKnn = Searcher::approxIncSearchFuzzy(root, queries + i * Const::tsLength,
                                                                                Const::k, Const::fuzzyidxfn,
                                                                                Const::visited_node_num);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
    }
}

void exactSearchDumpy() {
    Node *root = Node::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *exactKnn = Searcher::exactSearch(root, queries + i * Const::tsLength, Const::k, g);
        Const::logPrint("Results:");
        for (int j = 0; j < exactKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*exactKnn)[j]->ts) << endl;
    }
}

void exactSearchDumpyDTW() {
    Node *root = Node::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *exactKnn = Searcher::exactSearchDTW(root, queries + i * Const::tsLength, Const::k, g);
        Const::logPrint("Results:");
        for (int j = 0; j < exactKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*exactKnn)[j]->ts) << endl;
    }
}

void ngSearchDumpy() {
    Node *root = Node::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    root->assignLeafNum();
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *approxKnn = Searcher::ngSearch(root, queries + i * Const::tsLength,
                                                                    Const::k, Const::nprobes);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
    }
}

void ngSearchFuzzy() {
    Node *root = Node::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    root->assignLeafNum();
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *approxKnn = Searcher::ngSearchFuzzy(root, queries + i * Const::tsLength,
                                                                    Const::k, Const::nprobes);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
    }
}

void statIndexDumpy(){
    Node* root = Node::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    root->getIndexStats();
}

void statIndexFuzzy(){
    int bound = Const::fuzzy_f * 100;
    Const::fuzzyidxfn += "/" + to_string(bound) + "-" + to_string(Const::delta) + "/";
    Node* root = Node::loadFromDisk(Const::saxfn, Const::fuzzyidxfn + "root.idx", false);
    root->getIndexStats();
}

int main() {
    Const::readConfig();
    if(Const::index == 0) {
        constructGraph();
    } else if(Const::index == 1) {
        int ops = Const::ops;
        if(ops == 0) {
            buildDumpy();
        } else if (ops == 1) {
            approxSearchOneNode();
        } else if (ops == 2) {
            exactSearchDumpy();
        } else if (ops == 3) {
            statIndexDumpy();
        } else if (ops == 4) {
            approxSearchMoreNode();
        } else if (ops == 5) {
            approxSearchOneNodeDTW();
        } else if (ops == 6) {
            approxSearchMoreNodeDTW();
        } else if (ops == 7) {
            ngSearchDumpy();
        } else if (ops == 8) {
            exactSearchDumpyDTW();
        }
    } else if(Const::index == 2) {
        int ops = Const::ops;
        if(ops == 0) {
            buildFuzzy();
        } else if (ops == 1) {
            approxSearchOneNodeFuzzy();
        } else if (ops == 3) {
            statIndexFuzzy();
        } else if (ops == 4) {
            approxSearchMoreNodeFuzzy();
        } else if(ops == 7) {
            ngSearchFuzzy();
        }
    }
}
