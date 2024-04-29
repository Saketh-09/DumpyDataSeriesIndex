

#include <set>
#include <unordered_set>
#include <chrono>
#include <cmath>
#include <queue>
#include "../../include/Searchers/Searcher.h"
#include "../../include/DataStructures/PqItemSeries.h"
#include "../../include/DataStructures/TimeSeries.h"
#include "../../include/Utils/ConversionUtil.h"
#include "../../include/Utils/MathUtil.h"
#include "../../include/Const.h"

static Node* targetNode;
vector<PqItemSeries *> * Searcher::approxSearch(Node *root, float *query, int k, vector<vector<int>> *g,
                                                     const string &index_dir) {
    auto* queryTs = new TimeSeries(query);
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    int head = ConversionUtil::invSaxHeadFromSax(sax, Const::bitsCardinality, Const::segmentNum);
    Node *current = (root->ch)[head];
    if(current == nullptr){
        Node *node = nullptr;
        for(int i=0;i<Node::a + Node::b + Node::c;++i){
            if(root->ch[(*g)[head][i]] != nullptr){
                node = root->ch[(*g)[head][i]];
                break;
            }
        }
        // we only deal if the nearest node is a leaf or an internal node
        if(node->isInternalNode()){
            approxSearchInterNode(node, queryTs, sax, k, heap, index_dir);
        } else { node->search(k, queryTs, *heap, index_dir); targetNode = node;}
    } else if(!current->isInternalNode()){
        { current->search(k, queryTs, *heap, index_dir); targetNode = current;}
    } else approxSearchInterNode(current, queryTs, sax, k, heap, index_dir);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

void Searcher::approxSearchInterNode(Node *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                          vector<PqItemSeries *> *heap, const string &index_dir) {
    Node *current = root->route(sax);
    if(!current->isInternalNode()){
        current->search(k, queryTs, *heap, index_dir);
        targetNode = current;
        return;
    }

    // below is only for a empty pointer target leaf node, then we find the closest sibling
    double minimum_distance = numeric_limits<double>::max(), maximum_size = 0;
    Node *node;
    for(int index=0;index<current->ch.size();++index){
        if(current->ch[index] == nullptr)  continue;
        double distance;
        if(!current->ch[index]->isInternalNode())
            distance = ConversionUtil::LowerBound_Paa_iSax(queryTs->paa, current->sax, current->bits_cardinality, current->ch[index]->chosenSegments, index);
//        if(cur->children[i]->isLeafNode())  dist = ConversionUtil::LowerBound_Paa_iSax(queryTs->paa, cur->children[i]->sax, cur->children[i]->bits_cardinality);
        else distance = ConversionUtil::LowerBound_Paa_iSax(queryTs->paa, current->ch[index]->sax, current->ch[index]->bits_cardinality);
        if(distance < minimum_distance){
            minimum_distance = distance;
            maximum_size = current->ch[index]->size;
            node = current->ch[index];
        }else if(distance == minimum_distance && current->ch[index]->size > maximum_size){
            maximum_size = current->ch[index]->size;
            node = current->ch[index];
        }
    }

    // we only deal if the nearest node is a leaf or an internal node
    if(node->isInternalNode()){
        approxSearchInterNode(node, queryTs, sax, k, heap, index_dir);
        return;
    }else { node->search(k, queryTs, *heap, index_dir); targetNode = node;}
}

vector<PqItemSeries *> * Searcher::approxSearchDTW(Node *root, float *query, int k, vector<vector<int>> *g,
                                                        const string &index_dir) {
    auto* queryTs = new TimeSeries(query);
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    int head = ConversionUtil::invSaxHeadFromSax(sax, Const::bitsCardinality, Const::segmentNum);
    Node *cur = (root->ch)[head];
    if(cur == nullptr){
        Node *node = nullptr;
        for(int i=0;i<Node::a + Node::b + Node::c;++i){
            if(root->ch[(*g)[head][i]] != nullptr){
                node = root->ch[(*g)[head][i]];
                break;
            }
        }
        // we only concern whether the nearest node is a leaf or an internal node
        if(node->isInternalNode()){
            approxSearchInterNodeDTW(node, queryTs, sax, k, heap, index_dir);
        }else { node->searchDTW(k, queryTs, *heap, index_dir); targetNode = node;}
    }else if(!cur->isInternalNode()){
        { cur->searchDTW(k, queryTs, *heap, index_dir); targetNode = cur;}
    }else approxSearchInterNodeDTW(cur, queryTs, sax, k, heap, index_dir);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

void Searcher::approxSearchInterNodeDTW(Node *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                             vector<PqItemSeries *> *heap, const string &index_dir) {
    Node *cur = root->route(sax);
    if(!cur->isInternalNode()){
        cur->searchDTW(k, queryTs, *heap, index_dir);
        targetNode = cur;
        return;
    }

    auto* lowerLemire = new float [Const::tsLength];
    auto* upperLemire = new float [Const::tsLength];
    ConversionUtil::lower_upper_lemire(queryTs->ts, Const::tsLength, Const::dtw_window_size, lowerLemire, upperLemire);
    double *lowerPaa = ConversionUtil::paaFromTs(lowerLemire, Const::tsLengthPerSegment, Const::segmentNum);
    double *upperPaa = ConversionUtil::paaFromTs(upperLemire, Const::tsLengthPerSegment, Const::segmentNum);


    // below is only for a empty pointer target leaf node, then we find the closest sibling
    double min_dist = numeric_limits<double>::max(), max_size = 0;
    Node *node;
    for(auto & i : cur->ch){
        if(i == nullptr)  continue;
        double dist = ConversionUtil::minidist_paa_to_isax_DTW(upperPaa, lowerPaa, i->sax, i->bits_cardinality);
        if(dist < min_dist){
            min_dist = dist;
            max_size = i->size;
            node = i;
        }else if(dist == min_dist && i->size > max_size){
            max_size = i->size;
            node = i;
        }
    }

    // we only deal if the closest node is a leaf or an internal node
    if(node->isInternalNode()){
        approxSearchInterNodeDTW(node, queryTs, sax, k, heap, index_dir);
        return;
    }else { node->searchDTW(k, queryTs, *heap, index_dir); targetNode = node;}
}
static double *low_paa, *up_paa;
static float *t_paa;
bool comp_Dumpy_dtw(const Node* x, const Node* y){
    if(x == nullptr)    return false;
    if(y == nullptr)    return true;
    return ConversionUtil::minidist_paa_to_isax_DTW(up_paa, low_paa, x->sax, x->bits_cardinality) < ConversionUtil::minidist_paa_to_isax_DTW(up_paa, low_paa, y->sax, y->bits_cardinality);
}
vector<PqItemSeries *> * Searcher::approxIncSearchDTW(Node *root, float *query, int k, const string &index_dir,
                                                           int node_num) {
    auto* queryTs = new TimeSeries(query);
    auto* lowerLemire = new float [Const::tsLength];
    auto* upperLemire = new float [Const::tsLength];
    ConversionUtil::lower_upper_lemire(query, Const::tsLength, Const::dtw_window_size, lowerLemire, upperLemire);
    low_paa= ConversionUtil::paaFromTs(lowerLemire, Const::tsLengthPerSegment, Const::segmentNum);
    up_paa = ConversionUtil::paaFromTs(upperLemire, Const::tsLengthPerSegment, Const::segmentNum);

    t_paa = queryTs->paa;
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    approxIncSearchInterNodeDTW(root, queryTs, sax, k, heap, index_dir, node_num);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;

}

void Searcher::approxIncSearchInterNodeDTW(Node *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                                vector<PqItemSeries *> *heap, const string &index_dir,int &node_num) {
    if(!root->isInternalNode() || node_num <= 0)  return;
    Node *current = root->route1step(sax), *parent = root;
    while (current!= nullptr && current->isInternalNode() && current->getLeafNodeNum() > node_num) {
        parent = current;
        current = current->route(sax);
    }

    if(current!= nullptr){
        if(!current->isInternalNode()){
            current->searchDTW(k, queryTs, *heap, index_dir);
            --node_num;
        }else{
            approxIncSearchInterNodeDTW(current, queryTs, sax, k, heap, index_dir, node_num);
        }
    }

    if(node_num <=0)    return;

    vector<Node*>candidates;
    unordered_set<Node*>cands;
    for(Node *node: parent->ch)
        if(node != nullptr && node!=current && cands.find(node) == cands.end()) {
            candidates.push_back(node);
            cands.insert(node);
        }
    cands.clear();
    sort(candidates.begin(), candidates.end(), comp_Dumpy_dtw);



    for(int i=0;i<candidates.size() && node_num > 0;++i){
        if(!candidates[i]->isInternalNode()) {
            candidates[i]->searchDTW(k, queryTs, *heap, index_dir);
            --node_num;
        }
        else {
            approxIncSearchInterNodeDTW(candidates[i], queryTs, sax, k, heap, index_dir, node_num);
        }
    }

}

struct PqItem{
    Node* parent;
    int id;
    double distance{};

    PqItem(Node* _parent, int _id, double d){ parent = _parent;id= _id; distance = d;}
    PqItem(){ id=-1; distance = 0;parent = nullptr;}

    bool operator <(const PqItem & pt) const{
        if(distance != pt.distance)
            return distance < pt.distance;
        if(parent->layer != pt.parent->layer)
            return parent->layer > pt.parent->layer;
        return parent < pt.parent;
    }

    bool operator >(const PqItem& pt) const{
        if(distance != pt.distance)
            return distance > pt.distance;
        if(parent->layer != pt.parent->layer)
            return parent->layer < pt.parent->layer;
        return parent > pt.parent;
    }
};

struct cmp_PqItem{
    bool operator()(const PqItem& a, const PqItem& pt) const{
        if(a.distance != pt.distance)
            return a.distance > pt.distance;
        if(a.parent->layer != pt.parent->layer)
            return a.parent->layer < pt.parent->layer;
        return a.parent > pt.parent;
    }
};

vector<PqItemSeries*>*Searcher::exactSearch(Node* root, float *query, int k, vector<vector<int>> *g){

    vector<PqItemSeries*>* heap = approxSearch(root, query, k, g, Const::idxfn);
    unordered_set<Node*>visited;
    visited.insert(targetNode);
    make_heap(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    double bsf = (*heap)[0]->dist;
    auto *queryTs = new TimeSeries(query);

    set<PqItem>pq;
    for(int i =0;i<Const::vertexNum;++i){
        if(root->ch[i] == nullptr || visited.count(root->ch[i]) > 0)    continue;
        double dist  = ConversionUtil::getMinDist1stLayer(queryTs->paa , i);
        pq.insert(PqItem(root, i, dist));
    }

    double top_dist;
    Node* node;
    int length;
    while(!pq.empty()){
        top_dist = pq.begin()->distance;
        if(top_dist > bsf)  break;
        node = pq.begin()->parent->ch[pq.begin()->id];
        pq.erase(pq.begin());
        if(visited.count(node) > 0) continue;
        visited.insert(node);
        if(node->isInternalNode()){
            length = (1 << (node->chosenSegments.size()));
            for(int i =0;i<length;++i){
                if(node->ch[i] == nullptr)    continue;
                double dist = ConversionUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax,
                                                           node->bits_cardinality, node->chosenSegments, i);

                if(dist < bsf){
                    pq.insert(PqItem(node, i, dist));
                }
            }
        }else{
            node->search(k,queryTs,*heap,Const::idxfn);
            bsf = (*heap)[0]->dist;
        }
    }
    pq.clear();
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

vector<PqItemSeries*>*Searcher::exactSearchDTW(Node* root, float *query, int k, vector<vector<int>> *g){
    vector<PqItemSeries*>* heap = approxSearchDTW(root, query, k, g, Const::idxfn);
    unordered_set<Node*>visited;
    char _ = 0;
    visited.insert(targetNode);
    make_heap(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    double bsf = (*heap)[0]->dist;
    auto *queryTs = new TimeSeries(query);
    auto* lowerLemire = new float [Const::tsLength];
    auto* upperLemire = new float [Const::tsLength];
    ConversionUtil::lower_upper_lemire(query, Const::tsLength, Const::dtw_window_size, lowerLemire, upperLemire);
    double *lowerPaa = ConversionUtil::paaFromTs(lowerLemire, Const::tsLengthPerSegment, Const::segmentNum);
    double *upperPaa = ConversionUtil::paaFromTs(upperLemire, Const::tsLengthPerSegment, Const::segmentNum);


    priority_queue<PqItem, vector<PqItem>, cmp_PqItem> pq;
    for(int i =0;i<Const::vertexNum;++i){
        if(root->ch[i] == nullptr || visited.count(root->ch[i]) > 0)    continue;
        double dist = ConversionUtil::getMinDist1stLayerDTW(upperPaa, lowerPaa, i);
        pq.emplace(root, i, dist);
    }

    int len;
    PqItem cur; double top_dist, local_bsf; Node* node;

    while(!pq.empty()){
        cur = pq.top();
        pq.pop();
        top_dist = cur.distance;
        if(top_dist >= bsf)  break;
        node = cur.parent->ch[cur.id];
        if(visited.count(node) > 0) continue;
        visited.insert(node);
        if(node->isInternalNode()){
            len = (1 << (node->chosenSegments.size()));
            for(int i =0;i<len;++i){
                if(node->ch[i] == nullptr)    continue;
                double dist = ConversionUtil::minidist_paa_to_isax_DTW( upperPaa, lowerPaa, node->sax,
                                             node->bits_cardinality, node->chosenSegments, i);
                if(dist < bsf)
                    pq.emplace(node, i, dist);
            }
        }else{
            node->searchDTWSIMD(k, queryTs, *heap, Const::idxfn, upperLemire, lowerLemire);
            bsf = (*heap)[0]->dist;
        }
    }

    delete queryTs;
    sort_heap(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}


void searchSubTree(Node *root, TimeSeries *queryTs, int k, vector<PqItemSeries *> *heap, const string &index_dir,
                   int &node_num){
    unordered_set<Node*>visited;
    if(!root->isInternalNode()){
        if(root != targetNode) {
            root->search(k, queryTs, *heap, index_dir);
            --node_num;
        }
        return;
    }
    for(auto child:root->ch){
        if(child == nullptr || child == targetNode || visited.count(child) > 0)    continue;
        visited.insert(child);
        if(!child->isInternalNode()){
            child->search(k,queryTs,*heap, index_dir);
            --node_num;
        }else{
            searchSubTree(child, queryTs, k, heap, index_dir, node_num);
        }
    }
}

void searchSubTree(Node *root, TimeSeries *queryTs, int k, vector<PqItemSeries *> *heap, const string &index_dir,
                   int &node_num, unordered_set<float*, createhash, isEqual>*hash_set){
    unordered_set<Node*>visited;
    if(!root->isInternalNode()){
        if(root != targetNode) {
            root->search(k, queryTs, *heap, index_dir, hash_set);
            --node_num;
        }
        return;
    }
    for(auto child:root->ch){
        if(child == nullptr || child == targetNode || visited.count(child) > 0)    continue;
        visited.insert(child);
        if(!child->isInternalNode()){
            child->search(k,queryTs,*heap, index_dir, hash_set);
            --node_num;
        }else{
            searchSubTree(child, queryTs, k, heap, index_dir, node_num, hash_set);
        }
    }
}


vector<PqItemSeries *> *Searcher::ngSearch(Node *root, float *query, int k, int nprobes){
    auto *queryTs = new TimeSeries(query);
    auto heap = new vector<PqItemSeries*>();
    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    Node* root_subtree = root->route1step(sax);
    unordered_set<Node*>visited;
    if(root_subtree){
        if(!root_subtree->isInternalNode()){
            visited.insert(root_subtree);
            root_subtree->search(k, queryTs, *heap, Const::idxfn);
            --nprobes;
        }else if(root_subtree->leaf_num <= Const::pre_read){
            // a small subtree
            if(nprobes >= root_subtree->leaf_num){
                int _ =root_subtree->leaf_num;
                visited.insert(root_subtree);
                searchSubTree(root_subtree, queryTs, k, heap, Const::idxfn, _);
                nprobes -= root_subtree->leaf_num;
            }else{
                int rest = nprobes;
                approxIncSearchInterNode(root_subtree, queryTs, sax, k, heap, Const::idxfn, rest,  visited);
                nprobes = rest;
            }
        }else{
            // a large subtree
            int to_search = min(nprobes, Const::pre_read);
            int _ = to_search;
            approxIncSearchInterNode(root_subtree, queryTs, sax, k, heap, Const::idxfn, to_search,  visited);
            nprobes = nprobes - _ + to_search;
        }
    }

    if(nprobes <= 0){
        delete queryTs;
        sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
        return heap;
    }

    double bsf = heap->size() < k ? numeric_limits<double>::max(): (*heap)[0]->dist;

    set<PqItem>pq;
    for(int i =0;i<Const::vertexNum;++i){
        if(root->ch[i] == nullptr || visited.count(root->ch[i]) > 0)    continue;
        double dist  = ConversionUtil::getMinDist1stLayer(queryTs->paa , i);
        pq.insert(PqItem(root,i, dist));
    }
    int cur_probe = 0;
    while(!pq.empty() && cur_probe < nprobes){
        double top_dist;
        Node* node;
        top_dist = pq.begin()->distance;
        if(top_dist > bsf)  break;
        node = pq.begin()->parent->ch[pq.begin()->id];
        pq.erase(pq.begin());
        if(visited.count(node) > 0) continue;
        visited.insert(node);

        if(node->isInternalNode()){
            int len = (1 << (node->chosenSegments.size()));
            for(int i =0;i<len;++i){
                if(node->ch[i] == nullptr)    continue;
                double dist = ConversionUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax, node->bits_cardinality, node->chosenSegments, i);
                if(dist < bsf){
                    pq.insert(PqItem(node, i, dist));
                }
            }
        }else{
            node->search(k, queryTs, *heap, Const::idxfn);
            ++cur_probe;
            bsf = (*heap)[0]->dist;
        }
    }
    pq.clear();
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;

}

vector<PqItemSeries *> *Searcher::ngSearchFuzzy(Node *root, float *query, int k, int nprobes){
    auto *queryTs = new TimeSeries(query);
    auto heap = new vector<PqItemSeries*>();
    auto*hash_set = new unordered_set<float*, createhash, isEqual>();
    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    Node* root_subtree = root->route1step(sax);
    unordered_set<Node*>visited;
    if(root_subtree){
        if(!root_subtree->isInternalNode()){
            visited.insert(root_subtree);
            root_subtree->search(k, queryTs, *heap, Const::fuzzyidxfn, hash_set);
            --nprobes;
        }else if(root_subtree->leaf_num <= Const::pre_read){
            // a tiny subtree
            if(nprobes >= root_subtree->leaf_num) {
                int _ = root_subtree->leaf_num;
                visited.insert(root_subtree);
                searchSubTree(root_subtree, queryTs, k, heap, Const::fuzzyidxfn, _, hash_set);
                nprobes -= root_subtree->leaf_num;
            }else{
                int rest = nprobes;
                approxIncSearchInterNode(root_subtree, queryTs, sax, k, heap, Const::fuzzyidxfn, rest,  visited,
                                         hash_set);
                nprobes = rest;
            }
        }else{
            // a large subtree
            int to_search = min(nprobes, Const::pre_read);
            int _ = to_search;
            approxIncSearchInterNode(root_subtree, queryTs, sax, k, heap, Const::fuzzyidxfn, to_search,  visited,
                                     hash_set);
            nprobes = nprobes - _ + to_search;
        }
    }

    if(nprobes <= 0){
        delete queryTs;
        sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
        return heap;
    }

    double bsf = heap->size() < k ? numeric_limits<double>::max(): (*heap)[0]->dist;

    set<PqItem>pq;
    for(int i =0;i<Const::vertexNum;++i){
        if(root->ch[i] == nullptr || visited.count(root->ch[i]) > 0)    continue;
        double dist  = ConversionUtil::getMinDist1stLayer(queryTs->paa , i);
        pq.insert(PqItem(root,i, dist));
    }
    int cur_probe = 0;
    while(!pq.empty() && cur_probe < nprobes){
        double top_dist;
        Node* node;
        top_dist = pq.begin()->distance;
        if(top_dist > bsf)  break;
        node = pq.begin()->parent->ch[pq.begin()->id];
        pq.erase(pq.begin());
        if(visited.count(node) > 0) continue;
        visited.insert(node);

        if(node->isInternalNode()){
            int len = (1 << (node->chosenSegments.size()));
            for(int i =0;i<len;++i){
                if(node->ch[i] == nullptr)    continue;
                double dist = ConversionUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax, node->bits_cardinality, node->chosenSegments, i);
                if(dist < bsf){
                    pq.insert(PqItem(node, i, dist));
                }
            }
        }else{
            node->search(k, queryTs, *heap, Const::fuzzyidxfn, hash_set);
            ++cur_probe;
            bsf = (*heap)[0]->dist;
        }
    }
    pq.clear();
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;

}

bool comp_Dumpy(const Node* x, const Node* y){
    if(x == nullptr)    return false;
    if(y == nullptr)    return true;
    return ConversionUtil::LowerBound_Paa_iSax(t_paa, x->sax, x->layer) < ConversionUtil::LowerBound_Paa_iSax(t_paa, y->sax, y->layer);
}

vector<PqItemSeries *> * Searcher::approxIncSearch(Node *root, float *query, int k, const string &index_dir,
                                                        int node_num) {
    auto* queryTs = new TimeSeries(query);
    t_paa = queryTs->paa;
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    approxIncSearchInterNode(root, queryTs, sax, k, heap, index_dir, node_num);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;

}

void Searcher::approxIncSearchInterNode(Node *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                             vector<PqItemSeries *> *heap, const string &index_dir,int &node_num) {
    if(!root->isInternalNode() || node_num <= 0)  return;
    Node *cur = root->route1step(sax), *parent = root;
    while (cur!= nullptr && cur->isInternalNode() && cur->getLeafNodeNum() > node_num) {
        parent = cur;
        cur = cur->route(sax);
    }

    if(cur!= nullptr){
        if(!cur->isInternalNode()){
            cur->search(k, queryTs, *heap, index_dir);
            --node_num;
        }else{
            approxIncSearchInterNode(cur, queryTs, sax, k, heap, index_dir, node_num);
        }
    }

    if(node_num <=0)    return;

    vector<Node*>candidates;
    unordered_set<Node*>cands;
    for(Node *node: parent->ch)
        if(node != nullptr && node!=cur && cands.find(node) == cands.end()) {
            candidates.push_back(node);
            cands.insert(node);
        }
    cands.clear();
    sort(candidates.begin(), candidates.end(), comp_Dumpy);


    for(int i=0;i<candidates.size() && node_num > 0;++i){
        if(!candidates[i]->isInternalNode()) {
            candidates[i]->search(k, queryTs, *heap, index_dir);
            --node_num;
        }
        else {
            approxIncSearchInterNode(candidates[i], queryTs, sax, k, heap, index_dir, node_num);
        }
    }

}


vector<PqItemSeries *> * Searcher::approxIncSearchFuzzy(Node *root, float *query, int k, const string &index_dir,
                                                             int node_num) {
    auto* queryTs = new TimeSeries(query);
    t_paa = queryTs->paa;
    auto*heap = new vector<PqItemSeries*>();
    auto*hash_set = new unordered_set<float*, createhash, isEqual>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    approxIncSearchInterNodeFuzzy(root, queryTs, sax, k, heap, index_dir, node_num, hash_set);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;

}


void Searcher::approxIncSearchInterNodeFuzzy(Node *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                                  vector<PqItemSeries *> *heap, const string &index_dir, int &node_num,
                                                  unordered_set<float*, createhash, isEqual>*hash_set) {
    if(root->isLeafNode() || node_num <= 0)  return;
    Node *current = root->route1step(sax), *parent = root;
    while (current!= nullptr && !current->isLeafNode() && current->getLeafNodeNum() > node_num) {
        parent = current;
        current = current->route(sax);
    }

    if(current!= nullptr){
        if(current->isLeafNode()){
            current->search(k, queryTs, *heap, index_dir);
            --node_num;
        }else{
            approxIncSearchInterNode(current, queryTs, sax, k, heap, index_dir, node_num);
        }
    }

    if(node_num <=0)    return;

    vector<Node*>candidates;
    unordered_set<Node*>cands;
    for(Node *node: parent->ch)
        if(node != nullptr && cands.find(node) == cands.end()) {
            candidates.push_back(node);
            cands.insert(node);
        }
    cands.clear();

    sort(candidates.begin(), candidates.end(), comp_Dumpy);

    for(int i=0;i<candidates.size() && node_num > 0;++i){
        if(candidates[i]->isLeafNode()) {
            candidates[i]->search(k, queryTs, *heap, index_dir, hash_set);
            --node_num;
        }
        else {
            approxIncSearchInterNodeFuzzy(candidates[i], queryTs, sax, k, heap, index_dir, node_num, hash_set);
        }
    }

}

// for ng-search fuzzy
void Searcher::approxIncSearchInterNode(Node *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                             vector<PqItemSeries *> *heap, const string &index_dir,
                                             int &node_num, unordered_set<Node*>&visit,
                                             unordered_set<float*, createhash, isEqual>*hash_set) {
    if(node_num <= 0)  return;
    Node *cur = root->route1step(sax), *parent = root;
    while (cur!= nullptr && cur->isInternalNode() && cur->leaf_num > node_num) {
        parent = cur;
        cur = cur->route1step(sax);
    }

    if(cur!= nullptr){
        visit.insert(cur);
        if(!cur->isInternalNode()){
            cur->search(k, queryTs, *heap, index_dir, hash_set);
            --node_num;
        }else{
            searchSubTree(cur, queryTs,k, heap, index_dir, node_num, hash_set);
        }
    }

    if(node_num <=0)    return;

    double bsf = (*heap)[0]->dist;
    vector<PqItem>candidates;
    int len = (1 << (parent->chosenSegments.size()));
    for(int i =0;i<len;++i){
        if(parent->ch[i] == nullptr || parent->ch[i] == cur)    continue;
        double dist = ConversionUtil::LowerBound_Paa_iSax(queryTs->paa, parent->sax, parent->bits_cardinality, parent->chosenSegments, i);
        if(dist < bsf){
            candidates.emplace_back(parent, i, dist);
        }
    }
    sort(candidates.begin(),  candidates.end());

    for(int i=0;i<candidates.size() && node_num > 0;++i){
        Node* node = parent->ch[candidates[i].id];
        if(candidates[i].distance > bsf)    break;
        if(visit.count(node) > 0) continue;
        visit.insert(node);
        if(!node->isInternalNode()) {
            node->search(k, queryTs, *heap, index_dir, hash_set);
            --node_num;
        }
        else {
            searchSubTree(node, queryTs,k, heap, index_dir, node_num, hash_set);
        }
    }

}


// for ng-search
void Searcher::approxIncSearchInterNode(Node *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                             vector<PqItemSeries *> *heap, const string &index_dir,
                                             int &node_num, unordered_set<Node*>&visit) {
    if(node_num <= 0)  return;
    Node *current = root->route1step(sax), *parent = root;
    while (current!= nullptr && current->isInternalNode() && current->leaf_num > node_num) {
        parent = current;
        current = current->route1step(sax);
    }

    if(current!= nullptr){
        visit.insert(current);
        if(!current->isInternalNode()){
            current->search(k, queryTs, *heap, index_dir);
            --node_num;
        }else{
            searchSubTree(current, queryTs,k, heap, index_dir, node_num);
        }
    }

    if(node_num <=0)    return;

    double bsf = (*heap)[0]->dist;
    vector<PqItem>candidates;
    int len = (1 << (parent->chosenSegments.size()));
    for(int i =0;i<len;++i){
        if(parent->ch[i] == nullptr || parent->ch[i] == current)    continue;
        double dist = ConversionUtil::LowerBound_Paa_iSax(queryTs->paa, parent->sax, parent->bits_cardinality, parent->chosenSegments, i);
        if(dist < bsf){
            candidates.emplace_back(parent, i, dist);
        }
    }
    sort(candidates.begin(),  candidates.end());

    for(int i=0;i<candidates.size() && node_num > 0;++i){
        Node* node = parent->ch[candidates[i].id];
        if(candidates[i].distance > bsf)    break;
        if(visit.count(node) > 0) continue;
        visit.insert(node);
        if(!node->isInternalNode()) {
            node->search(k, queryTs, *heap, index_dir);
            --node_num;
        }
        else {
            searchSubTree(node, queryTs,k, heap, index_dir, node_num);
        }
    }

}
