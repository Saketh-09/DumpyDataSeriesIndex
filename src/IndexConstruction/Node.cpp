

#include <thread>
#include <cassert>
#include <cmath>
#include "../../include/DataStructures/Node.h"
#include "../../include/Utils/FileUtil.h"
#include "../../include/Utils/MathUtil.h"
#include "../../include/Utils/TimeSeriesUtil.h"
#include "../../include/Utils/ConversionUtil.h"

unsigned short *Node::saxes = nullptr;
float *Node::paas = nullptr;
int Node::a = MathUtil::nChooseK(Const::segmentNum, 1), Node::b = MathUtil::nChooseK(Const::segmentNum, 2), Node::c = MathUtil::nChooseK(Const::segmentNum, 3);
const int Node::power_2[]{1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536};
int * Node::combined_number = nullptr;
int*** Node::combines = nullptr;

struct pack{
    vector<bool> current_bits;
    vector<bool> current_mask;
    int tot_size;
    int process_id{};
    int masked_bits_num;
    bool disabled;

    pack(){
        tot_size = 0;
        masked_bits_num = 0;
        disabled = false;
    }
    pack(partUnit* node, int chosen_segment_number, int _pid){
        tot_size = node->s;
        masked_bits_num = 0;
        current_bits.resize(chosen_segment_number, false);
        current_mask.resize(chosen_segment_number, false);
        int _id =  node->id;
        for(int i=0;i<chosen_segment_number;++i){
            current_bits[chosen_segment_number - 1- i] = _id % 2;
            _id >>= 1;
        }
        process_id = _pid;
        node->process_id = process_id;
        disabled = false;
    }

    int calc_cost(int _id, int chosen_seg_num){
        int cost = 0;
        for(int index=0;index<chosen_seg_num;++index){
            if(!current_mask[chosen_seg_num-1-index] && current_bits[chosen_seg_num - 1 -index]!= (_id %2)){
                ++cost;
            }
            _id >>=1;
        }
        return cost;
    }

    int calc_pack_merge_cost(const pack & p, int chosen_seg_num, int *cur_cost, int *tar_cost){
        *cur_cost = 0; *tar_cost = 0;
        int cost = 0;
        for(int i=0;i<chosen_seg_num;++i){
            if(current_mask[i] && p.current_mask[i])    continue;
            if(current_mask[i] && !p.current_mask[i]){
                (*tar_cost)++;
                ++cost;
            }
            else if (!current_mask[i] && p.current_mask[i]) { (*cur_cost)++; ++cost;}
            else if(current_bits[i] != p.current_bits[i])   {
                (*cur_cost)++;
                (*tar_cost)++;
                ++cost;
            }
        }
        return cost;
    }

    void merge_pack(pack* p, int chosen_seg_num){
        pack * dis_one, * res_one;
        if(process_id < p->process_id){
            dis_one = p;
            res_one = this;
        }else{
            dis_one = this;
            res_one = p;
        }
        dis_one->disabled = true;
        res_one->tot_size = tot_size + p->tot_size;
        for(int i=0;i<chosen_seg_num;++i){
            if(dis_one->current_mask[i] && res_one->current_mask[i])    continue;
            if(dis_one->current_mask[i] && !res_one->current_mask[i]) { res_one->current_mask[i] = true; res_one->masked_bits_num++;}
            else if(!dis_one->current_mask[i] && res_one->current_mask[i]) continue;
            else if(dis_one->current_bits[i] != res_one->current_bits[i]){
                res_one->masked_bits_num++;
                res_one->current_mask[i] = true;
            }
        }

    }

    void insert(partUnit* node, int chosen_seg_num){
        node->process_id = process_id;
        int _id = node->id;
        tot_size += node->s;
        for(int i=0;i<chosen_seg_num;++i){
            if(!current_mask[chosen_seg_num-1-i] && current_bits[chosen_seg_num - 1 -i]!= (_id %2)){
                current_mask[chosen_seg_num-1-i] = true;
                masked_bits_num++;
            }
            _id>>=1;
        }
    }

    static bool compare_size(const pack&x, const pack&y){
        return x.tot_size < y.tot_size;
    }
};

void Node::loadCombines(){
    string base = "../combines/" + to_string(Const::segmentNum) + "-";
    auto ret = new int**[Const::segmentNum + 1];
    combined_number = new int[Const::segmentNum];
    ifstream ff("../combines/cnum-"+ to_string(Const::segmentNum) + ".txt", ios::in);
    for(int i=0;i<Const::segmentNum;++i){
        ff >> combined_number[i];
    }
    ff.close();

    for(int i=1;i<=Const::segmentNum - 1;++i){
        ret[i] = new int*[combined_number[i]];
        ifstream f(base + to_string(i) + ".txt", ios::in);
        for(int j=0;j<combined_number[i];++j){
            ret[i][j] = new int[i];
            for(int k=0;k<i;++k) {
                f >> ret[i][j][k];
            }
        }
        f.close();
    }
    ret[Const::segmentNum] = new int*[1];
    ret[Const::segmentNum][0] = new int[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i){
        ret[Const::segmentNum][0][i] = i;
    }
    combines = ret;

}

void materializeAllLeavesWithSax(string datafn, Node* root, int *navids, string index_dir, unsigned short*sax_tbl){
    auto start_sax = chrono::system_clock::now();
    Const::logPrint("Start move sax to disk file in 1st layer.");

    unordered_map<Node*, vector<unsigned short *>>sax_buffer;
    for(int i=0;i<root->size;++i){
        auto * sax = sax_tbl + i * Const::segmentNum;
        Node* node = root->route(sax);
        sax_buffer[node].push_back(sax);
    }
    for(auto &[node, buffer]:sax_buffer){
        string outfile = Const::idxfn + node->getFileName() + "_sax";
        if(node->partition_id == -1)    outfile += "_L";
        FILE *outf = fopen(outfile.c_str(), "a");
        for(auto sax:buffer)
            fwrite(sax, sizeof(unsigned short ), Const::segmentNum, outf);
        fclose(outf);
    }
    sax_buffer.clear();

    auto start_t = chrono::system_clock::now();
    Const::logPrint("Start move data to disk file in 1st layer.");
    FILE *f = fopen(datafn.c_str(), "r");
    long rest = root->size, total = root->size, current = 0;
    unordered_map<Node*, LBL_UNIT>lbl;

    // There is one more implementation method that fbl in each node adds a pointer vector where each pointer points to a series.
    // This may result many write calls.
    while(rest > 0){
        lbl.clear();
        long n;
        if(rest > Const::fbl_series_num)    n = Const::fbl_series_num;
        else n = rest;
        auto *tss = new float[n * Const::tsLength];

        auto end = chrono::system_clock::now();
        fread(tss, sizeof(float),n * Const::tsLength,  f);
        auto start = chrono::system_clock::now();

        //copy series to node to make sure write is only invoked one time upon a node fbl
        for(long index = current; index<current+n;++index){
            Node* node = root->route(sax_tbl + index * Const::segmentNum);
            lbl[node].buffer.push_back(tss + (index-current) * Const::tsLength);
        }

        // write series in order to node file from node fbl
        for(auto & [node,lbl_unit]:lbl){
            string out_file = Const::idxfn + node->getFileName();
            if(node->partition_id == -1)  out_file += "_L";
            FILE *outf = fopen(out_file.c_str(), "a");

            for(float *dat:lbl_unit.buffer)
                fwrite(dat, sizeof(float), Const::tsLength, outf);
            fclose(outf);
        }
        delete[] tss;

        rest-=n;
        current += n;
        Const::logPrint("Now Materialize all leaves. Progress: " + to_string((double)current / (double)total * 100) + "%");

    }

    fclose(f);
    delete[] navids;
}


Node* Node::route(const unsigned short *_sax){
    if(isLeafNode())
        return this;
    int nav_id = ConversionUtil::extendSax(_sax, bits_cardinality, chosenSegments);
    if(ch[nav_id] == nullptr) return this;
    return ch[nav_id]->route(_sax);
}

Node* Node::route1step(const unsigned short *_sax){
    assert(!isLeafNode());
    int nav_id;
    if(layer >= 1)
        nav_id = ConversionUtil::extendSax(_sax, bits_cardinality, chosenSegments);
    else
        nav_id = ConversionUtil::invSaxHeadFromSax(_sax, Const::bitsCardinality, Const::segmentNum);
    return ch[nav_id];
}

void Node::search(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir) const{
    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
    string fn = index_dir+ getFileName();


    FILE *f = fopen(fn.c_str(), "rb");
    struct timeval io{};
    Const::timerStart(&io);
    auto *ts = new float[size * Const::tsLength];
    fread(ts, sizeof(float), size * Const::tsLength, f);

    for(int i=0;i<size;++i){
        double dist = TimeSeriesUtil::euclideanDist(queryTs->ts, ts + i * Const::tsLength, Const::tsLength, bsf);

        if(heap.size() < k){
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }else if(dist < bsf){
            pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            delete heap.back();
            heap.pop_back();
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }

        if(heap.size() >= k)    bsf = heap[0]->dist;
    }

    for(PqItemSeries*s: heap){
        if(s->needDeepCopy) s->copyData();
    }
    delete[]ts;
    fclose(f);
}

void Node::search(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir, unordered_set<float*, createhash, isEqual>*hash_set) const{
    assert(isLeafNode());
    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
    string fn = index_dir+ getFileName();

    FILE *f = fopen(fn.c_str(), "rb");
    auto *ts = new float[size * Const::tsLength];
    fread(ts, sizeof(float), size * Const::tsLength, f);

    for(int i=0;i<size;++i){
        if(hash_set->find(ts + i * Const::tsLength) != hash_set->end()) continue;
        double dist = TimeSeriesUtil::euclideanDist(queryTs->ts, ts + i * Const::tsLength, Const::tsLength, bsf);

        if(heap.size() < k){
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            hash_set->insert(ts + i * Const::tsLength);
        }else if(dist < bsf){
            pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            hash_set->erase(heap.back()->ts);
            delete heap.back();
            heap.pop_back();
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            hash_set->insert(ts + i * Const::tsLength);
        }

        if(heap.size() >= k)    bsf = heap[0]->dist;
    }

    for(PqItemSeries*s: heap){
        if(s->needDeepCopy) s->copyData();
    }
    delete[]ts;
    fclose(f);
}

void Node::searchDTW(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir) const{
    assert(!isInternalNode());
    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
    string fn = index_dir+ getFileName();

    FILE *f = fopen(fn.c_str(), "rb");
    struct timeval io{};
    Const::timerStart(&io);
    auto *ts = new float[size * Const::tsLength];
    for(int i=0;i<size;++i)
        fread(ts + i * Const::tsLength, sizeof(float), Const::tsLength, f);

    for(int i=0;i<size;++i){
        double dist = TimeSeriesUtil::dtw(queryTs->ts, ts + i * Const::tsLength, Const::tsLength,Const::dtw_window_size, bsf);

        if(heap.size() < k){
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }else if(dist < bsf){
            pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            delete heap.back();
            heap.pop_back();
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }

        if(heap.size() >= k)    bsf = heap[0]->dist;
    }

    for(PqItemSeries*s: heap){
        if(s->needDeepCopy) s->copyData();
    }
    delete[]ts;
    fclose(f);
}

void Node::searchDTWSIMD(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir,
                              float* upperLemire, float* lowerLemire) const{
    assert(!isInternalNode());
    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
    string fn = index_dir+ getFileName();

    FILE *f = fopen(fn.c_str(), "rb");
    struct timeval io{};
    Const::timerStart(&io);
    auto *ts = new float[size * Const::tsLength];
    for(int i=0;i<size;++i)
        fread(ts + i * Const::tsLength, sizeof(float), Const::tsLength, f);

    float cb[Const::tsLength];
    float cb1[Const::tsLength];

    int length = 2*Const::dtw_window_size+1;
    float tSum[length];
    // pre_cost
    float pCost[length];
    // raw distance
    float rDist[length];
    double dist;


    for(int i=0;i<size;++i){
        float dist2=TimeSeriesUtil::lb_keogh_data_bound(ts + i * Const::tsLength, upperLemire,lowerLemire,
                                                        cb1, Const::tsLength, bsf);
        if(dist2 < bsf) {
            cb[Const::tsLength - 1] = cb1[Const::tsLength - 1];
            for (int ii = Const::tsLength - 2; ii >= 0; ii--)
                cb[ii] = cb[ii + 1] + cb1[ii];

            dist = TimeSeriesUtil::dtwsimd(queryTs->ts, ts + i * Const::tsLength, cb, Const::tsLength,
                                           Const::dtw_window_size, bsf,
                                           tSum, pCost, rDist);

            if(heap.size() < k){
                heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
                push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            }else if(dist < bsf){
                pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
                delete heap.back();
                heap.pop_back();
                heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
                push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            }

            if(heap.size() >= k)    bsf = heap[0]->dist;
        }
    }

    for(PqItemSeries*s: heap){
        if(s->needDeepCopy) s->copyData();
    }
    delete[]ts;
    fclose(f);
}


Node *Node::BuildingIndex(string &datafn, string &saxfn) {
    Const::logPrint("Start building index.");
    FileUtil::checkDirClean(Const::idxfn.c_str());
    loadCombines();
    long series_number = generateSaxTbl();

    auto* root = new Node();
    root->size = series_number;
    for(int &i:root->bits_cardinality)  i=0;
    partUnit nodeIn1stLayer[Const::vertexNum];
    int *navids = new int[series_number];
    for(int i=0;i<Const::vertexNum;++i)
        nodeIn1stLayer[i].id = i, nodeIn1stLayer[i].s=0, nodeIn1stLayer[i].process_id = -1;

    // obtain initial layer node size
    for(long index=0;index<series_number;++index){
        unsigned short *asax = saxes + index * Const::segmentNum;
        int nav_id = ConversionUtil::invSaxHeadFromSax(asax, Const::bitsCardinality, Const::segmentNum);
        navids[index] = nav_id;
        nodeIn1stLayer[nav_id].s++;
    }

    Const::logPrint("Finish statistic size of nodes in the 1st layer.");

    // partition 1st layer
    int partNum = partition(nodeIn1stLayer, Const::segmentNum);
    Const::logPrint("Finish partition");
    Node* childrenList[partNum];
    for(int index=0;index<partNum;++index)  childrenList[index] = new Node(1, index);
    root->ch.resize(Const::vertexNum);
    for(int i=0;i<Const::vertexNum;++i){
        if(nodeIn1stLayer[i].s <= 0) continue;
        if(nodeIn1stLayer[i].s > Const::th) {
            root->ch[i] = new Node(1, nodeIn1stLayer[i].s, i);
            root->ch[i]->generateSaxAndCardIn1stLayer(i);
        }else if(nodeIn1stLayer[i].process_id == -1){
            root->ch[i] = new Node(1, nodeIn1stLayer[i].s, i);
            root->ch[i]->generateSaxAndCardIn1stLayer(i);
        }
        else{
            int pid = nodeIn1stLayer[i].process_id;
            root->ch[i] = childrenList[pid];
            childrenList[pid]->size += nodeIn1stLayer[i].s;
            childrenList[pid]->generateSaxAndCardIn1stLayer4LeafNode(i);
        }
    }
    Const::logPrint("Finish build index structure 1st layer.");

    // add data offsets to internal nodes in first layer
    for(int i=0;i<Const::vertexNum;++i)
        if(nodeIn1stLayer[i].s > Const::th)
            root->ch[i]->offsets.reserve(nodeIn1stLayer[i].s);
    for(int i=0;i<series_number;++i){
        int nav_id = navids[i];
        root->ch[nav_id]->offsets.push_back(i);
    }
    Const::logPrint("data offsets have been put into nodes in the 1st layer.");

    int j = 0;
    int milestone = 0.1 * Const::vertexNum;
    Const::logPrint("start grow the index structure");
    for(int i=0;i<Const::vertexNum;++i){
        if(nodeIn1stLayer[i].s > Const::th) {
            root->ch[i]->growIndex();
        }
        if(++j%milestone == 0)
            Const::logPrint(to_string(j) + " nodes in the 1st layer has been processed.");
    }

    Const::logPrint("build index skeleton finished.");

    Const::logPrint("Start materialize leaves");
    materializeAllLeavesWithSax(datafn, root, navids, Const::idxfn, saxes);
    Const::logPrint("build index successfully!");
    delete[] saxes;

    return root;
}

void Node::growIndex() {
    if(size <= Const::th)   return;    
    determineSegments();
    int chosen_num = chosenSegments.size();
    
    // statistic children information in order to partition
    partUnit nodes[1<<chosen_num];
    for(int i=0;i<(1<<chosen_num);++i)
        nodes[i].id = i, nodes[i].s=0, nodes[i].process_id = -1;
    vector<vector<int>>node_offsets(1<<chosen_num, vector<int>());

    for(int i=0;i<size;++i){
        int new_id = ConversionUtil::extendSax(Node::saxes + (long)offsets[i] * (Const::segmentNum), bits_cardinality, chosenSegments);
        nodes[new_id].s++;
        node_offsets[new_id].push_back(offsets[i]);
    }

    if(this->layer > 1) vector<int>().swap(offsets);

    int partNum = partition(nodes, chosen_num);     
    
    Node* childrenList[partNum];
    for(int i=0;i<partNum;++i)  childrenList[i] = new Node(this, i);
    ch.resize(1 << chosen_num);
    for(int i=0;i<(1 << chosen_num);++i){
        if(nodes[i].s <= 0)  continue;
        else if(nodes[i].s > Const::th) {
            ch[i] = new Node(this, nodes[i].s, i);
            generateSaxAndCardinality(ch[i], i);
            ch[i]->offsets.resize(nodes[i].s);
            copy(node_offsets[i].begin(),  node_offsets[i].end(), ch[i]->offsets.begin());
            vector<int>().swap(node_offsets[i]);
        }else if(partition_id == -1){
            ch[i] = new Node(this, nodes[i].s, i);
            generateSaxAndCardinality(ch[i], i);
            vector<int>().swap(node_offsets[i]);
        }
        else{
            int _pid = nodes[i].process_id;
            ch[i] = childrenList[_pid];
            childrenList[_pid]->size += nodes[i].s;
            generateSaxAndCardinality4LeafNode(ch[i], i);
            vector<int>().swap(node_offsets[i]);
        }
    }
    
    vector<vector<int>>().swap(node_offsets);
    
    for(auto &child: ch){
        if(child!= nullptr && child->size > Const::th){
            child->growIndex();
        }
    }

}

void Node::determineFanout(int *lambda_minimum, int * lambda_max) const{
    if(size < 2 * Const::th)    {
        *lambda_minimum = 1;
        *lambda_max = 1;
        return;
    }
    *lambda_minimum = -1;
    *lambda_max = -1;
    double _min = size / (Const::th * Const::f_high);
    double _max = size / (Const::th * Const::f_low);
    for(int i = 1; i <= Const::segmentNum; ++i){
        if(*lambda_minimum == -1){
            if((1<< i) >= _min){
                *lambda_minimum = i;
            }
        }else{
            if((1<<i) == _max){
                *lambda_max = i;
                break;
            }else if((1<<i) > _max){
                *lambda_max = max(i-1,*lambda_minimum);
                break;
            }
        }
    }
}

void Node::determineSegments() {
    int lambda_min, lambda_max;
    determineFanout(&lambda_min, &lambda_max);

    vector<int>unit_size(Const::vertexNum, 0);

    vector<unordered_map<unsigned short, int>>data_seg_symbols(Const::segmentNum);
    for(int offset:offsets){
        unsigned short* cur_sax = saxes + offset * Const::segmentNum;
        for(int i=0;i<Const::segmentNum;++i){
            data_seg_symbols[i][cur_sax[i]]++;
        }
        int head = ConversionUtil::extendSax(cur_sax, bits_cardinality);
        unit_size[head]++;
    }

    // calculate stdev of each part
    vector<double>data_seg_mean(Const::segmentNum, 0);
    vector<double>data_seg_stdev(Const::segmentNum, 0);
    for(int i=0;i<Const::segmentNum;++i){
        auto& map = data_seg_symbols[i];
        for(auto &iter:map){
            unsigned short symbol = iter.first;
            data_seg_mean[i] += (ConversionUtil::getMidLineFromSaxSymbolbc8(symbol) * iter.second);
        }
        data_seg_mean[i] /= size;
        for(auto &iter:map){
            unsigned short symbol = iter.first;
            double mid_value = ConversionUtil::getMidLineFromSaxSymbolbc8(symbol);
            data_seg_stdev[i] += (iter.second * ((mid_value - data_seg_mean[i]) * (mid_value - data_seg_mean[i])));
        }
        data_seg_stdev[i] /= size;
    }

    vector<double>().swap(data_seg_mean);
    vector<unordered_map<unsigned short, int>>().swap(data_seg_symbols);

    // start to find the size of each node in each plan
    int plan_num;
    if(lambda_max < Const::segmentNum)
        plan_num = combined_number[lambda_max];
    else
        plan_num = 1;
    unordered_set<int>visited;
    double max_score = 0;
    vector<int> best_plan;
    for(int i=0;i<plan_num;++i){
        int *plan = combines[lambda_max][i];
        // first evaluate the whole plan
        vector<int>plan_node_sizes(1<<lambda_max, 0);
        int mask_code = MathUtil::generateMaskSettingKbits(plan, lambda_max, Const::segmentNum);
        map<int,int>max_node_size;
        for(int j=0;j<Const::vertexNum;++j){
            max_node_size[mask_code & j] += unit_size[j];
        }

        int _ = 0;
        for(auto & iter: max_node_size){
            plan_node_sizes[_++] = iter.second;
        }
        map<int,int>().swap(max_node_size);

        double score = compute_score(plan_node_sizes, plan, lambda_max, data_seg_stdev);
        if(score > max_score){
            max_score = score;
            best_plan.clear();
            for(int j = 0; j<lambda_max;++j)
                best_plan.push_back(plan[j]);
        }

        if(lambda_min <= lambda_max - 1)
            visitPlanFromBaseTable(visited, lambda_max - 1, plan, plan_node_sizes,
                                   &max_score, best_plan, lambda_min, mask_code, data_seg_stdev, score);
        vector<int>().swap(plan_node_sizes);
    }

    unordered_set<int>().swap(visited);
    vector<int>().swap(unit_size);
    chosenSegments = best_plan;
}

double Node::compute_score(vector<int>&node_sizes, int *plan, int lambda, vector<double>&data_seg_stdev) const{
    if(size < 2*Const::th){
        if(node_sizes[0] > Const::th || node_sizes[1] > Const::th)
            return (double)min(node_sizes[0],node_sizes[1]) / Const::th;
        return data_seg_stdev[plan[0]] * 100;
    }
    int over_th_nodes_no = 0;
    for(int _:node_sizes){
        if(_ > Const::th)
            over_th_nodes_no++;
    }
    double w = ((double)over_th_nodes_no) / ((int)node_sizes.size());
    double sum_seg = 0;
    for(int i=0;i<lambda;++i){
        sum_seg += data_seg_stdev[plan[i]];
    }
    sum_seg  = sqrt(sum_seg / lambda);
    sum_seg = exp(1+sum_seg);

    auto *tmp = new double[node_sizes.size()];
    for(int i=0;i<node_sizes.size();++i){
        tmp[i] = ((double)node_sizes[i]) / Const::th;
    }
    double stdev_fill_factor = MathUtil::deviation(tmp, node_sizes.size());

    double balance = exp(-(1+w) * stdev_fill_factor);
    double ret = sum_seg + Const::alpha* balance;
    delete[] tmp;
    return ret;
}

void
Node::visitPlanFromBaseTable(unordered_set<int> &visited, int cur_lambda, const int *plan, vector<int> &base_tbl,
                                  double *max_score, vector<int> &best_plan, int lambda_min, int mask_code,
                                  vector<double> &data_seg_stdev, double base_score) {
    // base mask is used to find base tbl
    int base_mask = 1;
    for(int i=0;i<cur_lambda;++i)
        base_mask = (base_mask << 1) +1;

    for(int i=0;i<cur_lambda + 1 ;++i){
        int reset_pos = plan[i];
        // find the whole plan code
        int cur_whole_mask = mask_code - (1 << (Const::segmentNum - 1  - reset_pos));
        if(visited.contains(cur_whole_mask))
            continue;
        visited.insert(cur_whole_mask);
        // find the new plan
        int *new_plan = new int[cur_lambda];
        int _=0;
        for(int j=0;j<cur_lambda+1;++j){
            if(i != j)  new_plan[_++] = plan[j];
        }

        // find the current base mask
        int cur_base_mask = base_mask - (1 << (cur_lambda - i));
        map<int,int>node_size_map;
        for(int j=0;j<base_tbl.size();++j)
            node_size_map[cur_base_mask &j] += base_tbl[j];

        vector<int>new_tbl(1<<cur_lambda, 0);
        _ = 0;
        for(auto &iter:node_size_map)
            new_tbl[_++] = iter.second;
        map<int,int>().swap(node_size_map);
        double score  = compute_score(new_tbl, new_plan, cur_lambda, data_seg_stdev);
        if(score > *max_score){
            *max_score = score;
            best_plan.clear();
            for(_ = 0; _<cur_lambda;++_)
                best_plan.push_back(new_plan[_]);
        }

        if(cur_lambda > lambda_min)
            visitPlanFromBaseTable(visited, cur_lambda - 1, new_plan, new_tbl,
                                   max_score, best_plan, lambda_min, cur_whole_mask, data_seg_stdev, score);
        vector<int>().swap(new_tbl);
        delete[] new_plan;
    }
}


PAA_INFO* Node::statPaa(){
    auto* r = new PAA_INFO();
    double split_line[Const::segmentNum],paa_max[Const::segmentNum],paa_min[Const::segmentNum],paa_mu[Const::segmentNum];

    double lb;  // ub is the new split line
    for(int i=0; i < Const::segmentNum; ++i)
        ConversionUtil::getValueRange(sax[i] << 1, bits_cardinality[i] + 1, &lb, &split_line[i]);
    for(auto &i:paa_max) i = - numeric_limits<double>::max();
    for(auto &i:paa_min) i = numeric_limits<double>::max();
    for(auto &i:r->paa_up_size) i = 0;
    for(auto &i:r->paa_below_size) i = 0;
    for(auto &i:r->paa_variance) i = 0;
    for(auto &i:paa_mu) i=0;
    for(long offset:offsets){
        float * start = paas + offset * (Const::segmentNum);
        for(int i=0;i<Const::segmentNum;++i){
            double value = *(start + i);
            paa_mu[i] += value;
            paa_min[i] = min(paa_min[i], value);
            paa_max[i] = max(paa_max[i], value);
            if(value > split_line[i]) {
                r->paa_up_size[i]++;
            }
            else {
                r->paa_below_size[i]++;
            }
        }
    }
    for(double & index : paa_mu) {
        index /= size;
    }

    for(long offset:offsets){
        float * start = paas + offset * (Const::segmentNum);
        for(int i=0;i<Const::segmentNum;++i){
            double value = *(start + i);
            r->paa_variance[i] += (value - paa_mu[i]) * (value - paa_mu[i]);
        }
    }
    return r;
}

struct temporary{
    int i{};
    double score{};
    temporary(int _i, double _score){i=_i;score = _score;}
    temporary(){;}

    static bool order(temporary a,temporary b){
        return a.score < b.score;
    }

    static bool orderdesc(temporary a,temporary b){
        return a.score > b.score;
    }
};

int Node::chooseOneSegment(PAA_INFO* node){
    int minimum = numeric_limits<int>::max(), min_index = -1;
    for(int i = 0; i < Const::segmentNum;++i){
        int large = max(node->paa_up_size[i], node->paa_below_size[i]);
        if(large < minimum){
            minimum = large;
            min_index = i;
        }
    }
    return min_index;
}

void Node::chooseSegment(PAA_INFO *paa, int chosen_number) {
    chosenSegments.resize(chosen_number);
    if(chosen_number == 1) { chosenSegments[0]=chooseOneSegment(paa) ; return;}

    temporary scores[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        if(bits_cardinality[i] >= Const::bitsCardinality)
            scores[i] = temporary(i, -1);
        else
            scores[i] = temporary(i, paa->paa_variance[i]);
    sort(scores, scores+Const::segmentNum, temporary::orderdesc);
    for(int i=0;i<chosen_number;++i)
        chosenSegments[i] = scores[i].i;
    sort(chosenSegments.begin(), chosenSegments.end());
}

void Node::generateSaxAndCardIn1stLayer(int new_id){
    for(int i = Const::segmentNum - 1; i >=0 ;--i){
        unsigned short t = new_id % 2 ;
        new_id >>= 1;
        sax[i] =  t;
    }
}

void Node::generateSaxAndCardinality(Node* node, int new_id){
    copy(sax, sax + Const::segmentNum, node->sax);
    copy(bits_cardinality, bits_cardinality + Const::segmentNum, node->bits_cardinality);
    for(int i = chosenSegments.size() - 1; i >=0 ;--i){
        int seg = (chosenSegments)[i];
        node->bits_cardinality[seg]++;
        int t = new_id % 2 ;
        new_id >>= 1;
        node->sax[seg] = (node->sax[seg] << 1) + t;
    }
}

void Node::generateSaxAndCardIn1stLayer4LeafNode(int new_id){
    if(bits_cardinality[0] == -1){
        for(int &i:bits_cardinality)    i=1;
        generateSaxAndCardIn1stLayer(new_id);
        return;
    }
    for(int i = Const::segmentNum - 1; i >=0 ;--i){
        int t = new_id % 2 ;
        new_id >>= 1;
        if(bits_cardinality[i] == 1 && sax[i] != t){
            bits_cardinality[i] = 0;
        }
    }
}

void Node::generateSaxAndCardinality4LeafNode(Node* node, int new_id){
    if(node->bits_cardinality[0] == -1){
        generateSaxAndCardinality(node, new_id);
        return;
    }
    for(int i = chosenSegments.size() - 1; i >=0 ;--i){
        int seg = chosenSegments[i];
        int t = new_id % 2 ;
        new_id >>= 1;
        if(node->bits_cardinality[seg] == bits_cardinality[seg] + 1 && node->sax[seg] % 2 != t){
            node->bits_cardinality[seg]--;
            node->sax[seg] >>= 1;
        }
    }
}

struct failUnit{
    partUnit *node{};
    int neighbor_size{};

    failUnit(partUnit* a, int b){ node = a; neighbor_size = b;}
};

static bool comp_fail_node(failUnit *x, failUnit *y){
    return x->neighbor_size > y->neighbor_size;
}

int Node::partition(partUnit* nodes_map, int chosen_segment_number){
    vector<partUnit*>nodes;
    int input_node_number = 1 << chosen_segment_number;
    int total_size = 0, node_number;
    for(int i=0;i<input_node_number;++i){
        if(nodes_map[i].s < Const::th * Const::small_perc && nodes_map[i].s > 0) {
            nodes.push_back(&nodes_map[i]);
            total_size += nodes_map[i].s;
        }
    }

    node_number = nodes.size();
    if(node_number <= 3) return 0;
    int pid = 0;
    int first_round_pack_num = (int)floor(total_size / (Const::th));
    if(node_number <= first_round_pack_num) return 0;
    int max_mask_num = floor(chosen_segment_number * Const::max_mask_bit_percentage);
    // 0.5 is a small number, generate more packs in the first round
    vector<pack>packs(first_round_pack_num);
    sort(nodes.begin(), nodes.end(), partUnit::comp_size);
    for(int i=0;i<first_round_pack_num;++i){
        packs[i] = pack(nodes[i], chosen_segment_number, i);
    }

    for(partUnit* cur_node:nodes) {
        if (cur_node->process_id != -1) continue;   // the node has been packed
        int cur_id = cur_node->id;

        int minimum_cost = chosen_segment_number;
        int pack_id = -1;
        for (pack &p: packs) {
            if (p.tot_size + cur_node->s > Const::th || p.masked_bits_num >= max_mask_num) continue;
            int cost = p.calc_cost(cur_id, chosen_segment_number);
            if(cost + p.masked_bits_num >= max_mask_num)    continue;
            if (cost < minimum_cost) {
                minimum_cost = cost;
                pack_id = p.process_id;
            }
        }
        if (pack_id == -1) {
            packs.emplace_back(cur_node, chosen_segment_number, packs.size());
        }else {
            packs[pack_id].insert(cur_node, chosen_segment_number);
        }
    }

    // merge packs
    unordered_map<int,int>process_id_map;
    for(int i = 0;i<packs.size();++i)
        process_id_map[i] = i;
    sort(packs.begin(),packs.end(), pack::compare_size);
    for(int i=0;i<packs.size();++i){
        pack& cur_pack = packs[i];
        if(cur_pack.process_id != process_id_map[cur_pack.process_id])    continue;

        int minimum_cost = chosen_segment_number;
        int minimum_size = numeric_limits<int>::max();
        int min_pack_id = -1;
        int cur_cost, tar_cost;
        for(int j=0;j<packs.size();++j){
            pack&target_pack = packs[j];
            if(target_pack.disabled || cur_pack.process_id == target_pack.process_id || cur_pack.tot_size + target_pack.tot_size > Const::th
               || cur_pack.masked_bits_num >= max_mask_num || target_pack.masked_bits_num >= max_mask_num) continue;
            cur_pack.calc_pack_merge_cost(target_pack, chosen_segment_number, &cur_cost, &tar_cost);
            if(cur_cost + cur_pack.masked_bits_num >= max_mask_num ||
               tar_cost + target_pack.masked_bits_num >= max_mask_num) continue;
            int cost = cur_cost + tar_cost;
            if(cost < minimum_cost || (cost == minimum_cost && cur_pack.tot_size < minimum_size)){
                minimum_cost = cost;
                minimum_size = target_pack.tot_size;
                min_pack_id = j;
            }
        }

        if(minimum_size < numeric_limits<int>::max()){
            cur_pack.merge_pack(&packs[min_pack_id], chosen_segment_number);
        }
    }

    // re-assign the process ids to the nodes
    int maximum_pid = 0;
    for(partUnit* node:nodes) {
        node->process_id = process_id_map[node->process_id];
        maximum_pid = max(maximum_pid, node->process_id);
    }

    return maximum_pid + 1;
}

void Node::getIndexStats(){
    int total_leaf_node_num = getLeafNodeNum();
    int total_size = getTotalSize();
    cout << "Total size = " << total_size << endl;
    cout <<"Total nodes number = " << getNodeNum() << endl;
    cout << "Leaf node number = " << total_leaf_node_num << endl;
    cout << "1st layer node number = " << get1stLayerNodesNo() <<endl;
    cout << "1st layer internal node number = " << get1stLayerInterNodesNo() << endl;
    cout << "1st layer internal series number = " << get1stLayerInterNodeSeriesNo() << endl;
    cout << "Max. height = " << getMaxHeight() - 1 <<endl;
    cout << "Avg. Height = " << getSumHeight() / (double) total_leaf_node_num << endl;
    cout <<"Avg. Filling Factor = "<< total_size / (double)total_leaf_node_num / Const::th << endl;
    cout << "Bias leaf node ratio = " << (double)getBiasLeafNodeNum() / total_leaf_node_num << endl;
}

int Node::get1stLayerInterNodesNo(){
    unordered_set<Node*>node;
    for(Node* child:ch){
        if(child== nullptr || child->size <= Const::th || node.find(child) != node.end())
            continue;
        node.insert(child);
    }
    return node.size();
}

int Node::get1stLayerNodesNo(){
    unordered_set<Node*>node;
    for(Node* child:ch){
        if(child== nullptr || node.find(child) != node.end())
            continue;
        node.insert(child);
    }
    return node.size();
}

int Node::get1stLayerInterNodeSeriesNo(){
    unordered_set<Node*>node;
    int ret = 0;
    for(Node* child:ch){
        if(child== nullptr || child->size <= Const::th || node.find(child) != node.end())
            continue;
        node.insert(child);
        ret += child->size;
    }
    return ret;
}

int Node::getMaxHeight(){
    if(isLeafNode())    return 1;
    int max_height = 0;
    unordered_set<Node*>hash_map;
    for(Node*child:ch){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        max_height = max(child->getMaxHeight(), max_height);
        hash_map.insert(child);
    }
    return max_height + 1;
}

int Node::getLeafNodeNum(){
    if(isLeafNode())    return 1;
    int sum = 0;
    unordered_set<Node*>hash_map;
    for(Node*child:ch){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        sum += child->getLeafNodeNum();
        hash_map.insert(child);
    }
    return sum;
}

// isax bias leaf nodes num
int Node::getBiasLeafNodeNum(){
    if(isLeafNode()){
        int max_b = 0, min_b=8;
        for(int &bc:bits_cardinality){
            max_b = max(max_b, bc);
            min_b = min(min_b, bc);
        }
        return (max_b - min_b >= 4);
    }
    int sum = 0;
    unordered_set<Node*>hash_map;
    for(Node*child:ch){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        sum += child->getBiasLeafNodeNum();
        hash_map.insert(child);
    }
    return sum;
}

int Node::getTotalSize(){
    if(isLeafNode())    return size;
    int sum = 0;
    unordered_set<Node*>hash_map;
    for(Node*child:ch){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        sum += child->getTotalSize();
        hash_map.insert(child);
    }
    return sum;
}

int Node::getNodeNum(){
    if(isLeafNode())    return 1;
    int sum = 0;
    unordered_set<Node*>hash_map;
    for(Node*child:ch){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        sum += child->getNodeNum();
        hash_map.insert(child);
    }
    return sum + 1;
}

int Node::getSumHeight(){
    if(isLeafNode())    return layer;
    int sum_height = 0;
    unordered_set<Node*>hash_map;
    for(Node*child:ch){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        sum_height += child->getSumHeight();
        hash_map.insert(child);
    }
    return sum_height;
}

int Node::loadSax(const string & saxfn){
    long f_size = FileUtil::getFileSize(saxfn.c_str()), series_num = f_size / (sizeof(unsigned short) * Const::segmentNum);
    saxes = new unsigned short [f_size / sizeof(unsigned short )];
    FILE *f = fopen(saxfn.c_str(), "rb");
    fread(saxes, sizeof(unsigned short ), f_size / sizeof(unsigned short), f);
    fclose(f);
    Const::logPrint("Finish loading sax");
    return series_num;
}

void Node::loadPaa(const string & paafn){
    long f_size = FileUtil::getFileSize(paafn.c_str());
    paas = new float [f_size / sizeof(float )];
    FILE *f = fopen(paafn.c_str(), "rb");
    fread(paas, sizeof(float ), f_size / sizeof(float ), f);
    fclose(f);
    Const::logPrint( "Finish loading paa");
}

long Node::generateSaxAndPaaTbl(){
    string fn = Const::datafn;
    long fs = FileUtil::getFileSize(fn.c_str());
    long series_num = fs / Const::tsLengthBytes;
    cout << "Total Series Number is "<<series_num <<endl;
    float ts[Const::tsLength];
    saxes = new unsigned short[series_num * Const::segmentNum];
    paas = new float[series_num * Const::segmentNum];
    long rest = series_num, cur = 0;
    FILE *f = fopen(fn.c_str(), "rb");


    while(rest > 0){
        unsigned num;
        if(rest > 4000000)    num = 4000000;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        fread(tss, sizeof(float),num * Const::tsLength,  f);

        for(int i=0;i<num;++i){
            if(isnan(tss[i * Const::tsLength])){
                for(int j = 0; j < Const::segmentNum; ++j)
                    saxes[i * Const::segmentNum + j] = 0;
                cout << "Dirty data: "<<i << "," <<endl;
            }
            else{
                ConversionUtil::paaAndSaxFromTs(tss + i * Const::tsLength,
                                         paas + (cur+ i) * Const::segmentNum, saxes + (cur+ i) * Const::segmentNum,
                                         Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
            }
        }
        delete[] tss;
        rest -=num;
        cur+=num;
    }

    fclose(f);
    return series_num;
}

long Node::generateSaxTbl(){
    string fn = Const::datafn;
    long series_num;
    if(Const::series_num == -1) {
        long fs = FileUtil::getFileSize(fn.c_str());
        series_num = fs / Const::tsLengthBytes;
    }else{
        series_num = Const::series_num;
    }
    cout << "Total Series Number is "<<series_num <<endl;
    float ts[Const::tsLength];
    saxes = new unsigned short[series_num * Const::segmentNum];
    long rest = series_num, cur = 0;
    FILE *f = fopen(fn.c_str(), "rb");


    while(rest > 0){
        unsigned num;
        if(rest > 2000000)    num = 2000000;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        fread(tss, sizeof(float),num * Const::tsLength,  f);

        for(int i=0;i<num;++i){
            if(isnan(tss[i * Const::tsLength])){
                for(int j = 0; j < Const::segmentNum; ++j)
                    saxes[i * Const::segmentNum + j] = 0;
                cout << "Dirty data: "<<i << "," <<endl;
            }
            else{
                ConversionUtil::saxFromTs(tss + i * Const::tsLength, saxes + (cur+ i) * Const::segmentNum,
                                   Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
            }
        }
        delete[] tss;
        rest -=num;
        cur+=num;
    }

    fclose(f);
    return series_num;
}

Node *Node::loadFromDisk(const string &saxfn, const string &idxfn, bool need_sax) {
    if(need_sax)
        loadSax(saxfn);
    ifstream ifs(idxfn, ios::binary);
    boost::archive::binary_iarchive ia(ifs);
    auto *g = new Node();
    ia >> (*g);
    ifs.close();
    return g;
}

int Node::assignLeafNum() {
    if(!isInternalNode()) {
        leaf_num = 1;
        return 1;
    }

    unordered_set<Node*>visited;
    for(Node* child: ch){
        if(child == nullptr || visited.count(child) > 0)    continue;
        visited.insert(child);
        leaf_num += child->assignLeafNum();
    }

    return leaf_num;
}