

#ifndef DUMPY_Searcher_H
#define DUMPY_Searcher_H
#include <vector>
#include "../DataStructures/Node.h"

class Searcher {
public:
    static void approxSearchInterNode(Node *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                      vector<PqItemSeries *> *heap, const string &index_dir);

    static vector<PqItemSeries *> *
    approxSearch(Node *root, float *query, int k, vector<vector<int>> *g, const string &index_dir);

    static vector<PqItemSeries *> *exactSearch(Node *root, float *query, int k, vector<vector<int>> *g);

    static vector<PqItemSeries*>* exactSearchDTW(Node* root, float *query, int k, vector<vector<int>> *g);

    static vector<PqItemSeries *> *ngSearch(Node *root, float *query, int k, int nprobes);

    static vector<PqItemSeries *> *ngSearchFuzzy(Node *root, float *query, int k, int nprobes);

    static vector<PqItemSeries *> *
    approxIncSearch(Node *root, float *query, int k, const string &index_dir, int node_num);

    static void
    approxIncSearchInterNode(Node *root, TimeSeries *queryTs, unsigned short *sax, int k,
                             vector<PqItemSeries *> *heap, const string &index_dir,
                             int &node_num);

    static void approxIncSearchInterNodeFuzzy(Node *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                              vector<PqItemSeries *> *heap, const string &index_dir, int &node_num,
                                              unordered_set<float *, createhash, isEqual> *hash_set);

    static void approxIncSearchInterNode(Node *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                                 vector<PqItemSeries *> *heap, const string &index_dir,
                                                 int &node_num, unordered_set<Node*>&visit,
                                                 unordered_set<float*, createhash, isEqual>*hash_set);

    static void approxIncSearchInterNode(Node *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                                        vector<PqItemSeries *> *heap, const string &index_dir,
                                                        int &node_num, unordered_set<Node*>&visit);

    static vector<PqItemSeries *> *
    approxIncSearchFuzzy(Node *root, float *query, int k, const string &index_dir, int node_num);

    static vector<PqItemSeries *> *
    approxSearchDTW(Node *root, float *query, int k, vector<vector<int>> *g, const string &index_dir);

    static void
    approxSearchInterNodeDTW(Node *root, TimeSeries *queryTs, unsigned short *sax, int k,
                             vector<PqItemSeries *> *heap,
                             const string &index_dir);

    static void approxIncSearchInterNodeDTW(Node *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                                           vector<PqItemSeries *> *heap, const string &index_dir,int &node_num);

    static vector<PqItemSeries *> * approxIncSearchDTW(Node *root, float *query, int k, const string &index_dir,
                                                                      int node_num);
};


#endif //DUMPY_Searcher_H
