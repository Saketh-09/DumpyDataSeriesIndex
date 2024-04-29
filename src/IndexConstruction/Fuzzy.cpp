#include "../../include/DataStructures/Node.h"
#include "../../include/Utils/FileUtil.h"
#include "../../include/Utils/MathUtil.h"
#include <thread>

// We suggest you to  build Dumpy-fuzzy with the table of PAA as given.
// Only with SAX table is also viable.

static int fuzzy_number = 0;
int *Node::mask = nullptr;

// add actual series in the disk file of nodes in first layer
void materialize1stLayerFuzzy(string datafn, Node *root, int *navids, string index_dir,
                              unordered_map<Node *, NODE_RECORDER> *navigating_tbl) {
    Const::logPrint("Start move data to disk file in 1st layer.");
    FILE *f = fopen(datafn.c_str(), "r");
    long rest = root->size, total = root->size, current = 0;
    unordered_map<Node *, LBL_UNIT> fbl;

    // There is one more implementation method that fbl in every node stores a pointer vector where each pointer points to a series.
    // This may be incurring many write calls.
    while (rest > 0) {
        fbl.clear();
        long number;
        if (rest > Const::fbl_series_num) number = Const::fbl_series_num;
        else number = rest;
        auto *tss = new float[number * Const::tsLength];
        fread(tss, sizeof(float), number * Const::tsLength, f);

        long bound = current + number;
        // statistic the length of every node fbl size, and assign storage for them
        for (long i = current; i < bound; ++i) {
            fbl[root->ch[navids[i]]].buffer.push_back(tss + (i - current) * Const::tsLength);
        }

        for (auto &iter: *navigating_tbl) {
            if (iter.first->partition_id != -1) {
                NODE_RECORDER &recorder = iter.second;
                int &cur_pos = recorder.actual_size;
                while (cur_pos < recorder.series_index_list.size() && recorder.series_index_list[cur_pos] < bound) {
                    int index = recorder.series_index_list[cur_pos];
                    fbl[iter.first].buffer.push_back(tss + (index - current) * Const::tsLength);
                    ++cur_pos;
                }
            }
        }

        // add series in order to node file from node fbl
        for (auto &iter: fbl) {
            string outfile = index_dir;
            if (iter.first->partition_id == -1)
                outfile += "U_" + to_string(iter.first->id);
            else outfile += to_string(iter.first->layer) + "_" + to_string(iter.first->partition_id);
            FILE *outf = fopen(outfile.c_str(), "a");

            for (int i = 0; i < iter.second.buffer.size(); ++i)
                fwrite(iter.second.buffer[i], sizeof(float), Const::tsLength, outf);
            fclose(outf);
        }
        delete[] tss;

        rest -= number;
        current += number;
        Const::logPrint(
            "Now in 1st layer " + to_string((double) current / (double) total * 100) +
            "% series have been written to disk.(FIRST STAGE)");
    }

    fseek(f, 0, SEEK_SET);
    rest = root->size;
    current = 0;
    while (rest > 0) {
        fbl.clear();
        long number;
        if (rest > Const::fbl_series_num) number = Const::fbl_series_num;
        else number = rest;
        auto *tss = new float[number * Const::tsLength];
        fread(tss, sizeof(float), number * Const::tsLength, f);
        long bound = current + number;

        for (auto &iter: *navigating_tbl) {
            // internal node
            if (iter.first->partition_id == -1) {
                NODE_RECORDER &recorder = iter.second;
                int &current_position = recorder.actual_size;
                while (current_position < recorder.series_index_list.size() && recorder.series_index_list[
                           current_position] < bound) {
                    int i = recorder.series_index_list[current_position];
                    //                    fbl[iter.first].offsets.push_back(index);
                    fbl[iter.first].buffer.push_back(tss + (i - current) * Const::tsLength);
                    ++current_position;
                }
            }
        }

        for (auto &iter: fbl) {
            string out_file = index_dir;
            assert(iter.first->partition_id == -1);
            out_file += "U_" + to_string(iter.first->id);
            FILE *outf = fopen(out_file.c_str(), "a");


            for (int i = 0; i < iter.second.buffer.size(); ++i)
                fwrite(iter.second.buffer[i], sizeof(float), Const::tsLength, outf);
            fclose(outf);
        }

        delete[] tss;

        rest -= number;
        current += number;
        Const::logPrint(
            "Now in 1st layer " + to_string((double) current / (double) total * 100) +
            "% series have been written to disk.(SECOND STAGE)");
    }

    fclose(f);
    delete[] navids;
    navigating_tbl->clear();
}

void materializeInterNodeFuzzy(Node *node, unsigned short *saxes, int actual_size,
                               unordered_map<Node *, NODE_RECORDER> &navigating_tbl) {
    FILE *f = fopen((Const::fuzzyidxfn + "U_" + to_string(node->id)).c_str(), "r");
    long rest = node->size, cur = 0, number, bound;
    unordered_map<Node *, LBL_UNIT> lbl;

    while (rest > 0) {
        lbl.clear();
        if (rest > Const::fbl_series_num) number = Const::fbl_series_num;
        else number = rest;
        auto *tss = new float[number * Const::tsLength];
        fread(tss, sizeof(float), number * Const::tsLength, f);
        bound = cur + number;

        // actual series routing and inserting to lbl
        for (long i = cur; i < bound && i < actual_size; ++i) {
            Node *target = node->route(saxes + node->offsets[i] * Const::segmentNum);
            lbl[target].buffer.push_back(tss + (i - cur) * Const::tsLength);
        }

        // fuzzy series searching
        for (auto &iter: navigating_tbl) {
            NODE_RECORDER &recorder = iter.second;
            int &current_position = recorder.actual_size;
            while (current_position < recorder.series_index_list.size() && recorder.series_index_list[current_position]
                   < bound) {
                int i = recorder.series_index_list[current_position];
                lbl[iter.first].buffer.push_back(tss + (i - cur) * Const::tsLength);
                ++current_position;
            }
        }

        // add one series at a time
        for (auto &iter: lbl) {
            string out_file = Const::fuzzyidxfn + iter.first->getFileName();
            FILE *outf = fopen(out_file.c_str(), "a");


            for (int i = 0; i < iter.second.buffer.size(); ++i)
                fwrite(iter.second.buffer[i], sizeof(float), Const::tsLength, outf);

            fclose(outf);
        }

        delete[]tss;
        rest -= number;
        cur += number;
    }
    fclose(f);
    FileUtil::fileRemove((Const::fuzzyidxfn + "U_" + to_string(node->id)).c_str());
}

Node *Node::BuildIndexFuzzy(const string &datafn, const string &saxfn, const string &paafn,
                                      vector<vector<int> > *g) {
    FileUtil::checkDirClean(Const::fuzzyidxfn.c_str());
    Const::logPrint("start building index.");

    long series_number = generateSaxAndPaaTbl();

    mask = MathUtil::generateMask(Const::segmentNum);
    auto *root = new Node();
    root->size = series_number;
    for (int &i: root->bits_cardinality) i = 0;
    partUnit nodeIn1stLayer[Const::vertexNum];

    int *navids = new int[series_number];
    for (int i = 0; i < Const::vertexNum; ++i)
        nodeIn1stLayer[i].id = i, nodeIn1stLayer[i].s = 0, nodeIn1stLayer[i].process_id = -1;

    auto *paa_mu_part_units = new vector<vector<double> >(Const::vertexNum, vector<double>(Const::segmentNum, 0));
    // fetch initial layer node length
    for (long i = 0; i < series_number; ++i) {
        unsigned short *asax = saxes + i * Const::segmentNum;
        int navId = ConversionUtil::invSaxHeadFromSax(asax, Const::bitsCardinality, Const::segmentNum);
        float *paa = paas + i * Const::segmentNum;
        for (int j = 0; j < Const::segmentNum; ++j)
            (*paa_mu_part_units)[navId][j] += paa[j];
        navids[i] = navId;
        nodeIn1stLayer[navId].s++;
    }

    Const::logPrint("Finish statistic size of nodes in the 1st layer.");

    // partition first layer
    int part_num = partition(nodeIn1stLayer, Const::segmentNum);
    Const::logPrint("Finish partition");
    Node *children_list[part_num];
    for (int index = 0; index < part_num; ++index) children_list[index] = new Node(1, index);
    root->ch.resize(Const::vertexNum);
    int internalSize = 0;
    for (int index = 0; index < Const::vertexNum; ++index) {
        if (nodeIn1stLayer[index].s <= 0) continue;
        for (int j = 0; j < Const::segmentNum; ++j)
            (*paa_mu_part_units)[index][j] /= nodeIn1stLayer[index].s;
        if (nodeIn1stLayer[index].s > Const::th) {
            root->ch[index] = new Node(1, nodeIn1stLayer[index].s, index);
            root->ch[index]->generateSaxAndCardIn1stLayer(index);
            internalSize += root->ch[index]->size;
        } else if (nodeIn1stLayer[index].process_id == -1) {
            root->ch[index] = new Node(1, nodeIn1stLayer[index].s, index);
            root->ch[index]->generateSaxAndCardIn1stLayer(index);
        } else {
            int pid = nodeIn1stLayer[index].process_id;
            root->ch[index] = children_list[pid];
            children_list[pid]->size += nodeIn1stLayer[index].s;
            children_list[pid]->generateSaxAndCardIn1stLayer4LeafNode(index);
        }
    }
    Const::logPrint("Finish build index structure 1st layer.");

    // add data offsets to internal nodes in first layer
    for (int i = 0; i < Const::vertexNum; ++i) {
        if (root->ch[i] != nullptr) {
            root->ch[i]->offsets.reserve(nodeIn1stLayer[i].s);
        }
    }

    for (int i = 0; i < series_number; ++i) {
        int nav_id = navids[i];
        root->ch[nav_id]->offsets.push_back(i);
    }
    Const::logPrint("Start fuzzy 1st layer.");

    unordered_map<Node *, NODE_RECORDER> FLNT;
    root->fuzzyFirstLayer(nodeIn1stLayer, navids, FLNT, *paa_mu_part_units);
    delete paa_mu_part_units;
    Const::logPrint("1st layer fuzzy number is " + to_string(fuzzy_number));


    thread IO(materialize1stLayerFuzzy, datafn, root, navids, Const::fuzzyidxfn, &FLNT);

    int finished_series = 0, finished_percent = 0;
    Const::logPrint("start to grow the index structure");
    unordered_map<int, unordered_map<Node *, NODE_RECORDER> *> DPNT;
    for (int index = 0; index < Const::vertexNum; ++index) {
        if (nodeIn1stLayer[index].s > Const::th) {
            Node *ch = root->ch[index];
            // make sure the offset and series index list in DPLT is in sync with FLNT(file)
            sort(ch->offsets.begin() + nodeIn1stLayer[index].s, ch->offsets.end());
            auto *navigating_tbl = new unordered_map<Node *, NODE_RECORDER>();
            navigating_tbl->emplace(root->ch[index], NODE_RECORDER(nodeIn1stLayer[index].s, root->ch[index]));
            DPNT.emplace(index, navigating_tbl);
            root->ch[index]->growIndexFuzzy(*navigating_tbl, g);

            finished_series += nodeIn1stLayer[index].s;
            double percentage = (double) finished_series / (double) internalSize;
            if (percentage >= (finished_percent + 1) * 0.1) {
                Const::logPrint(to_string(percentage * 100) + "% internal series in the 1st layer has been processed.");
                finished_percent = percentage * 10;
            }
        }
    }
    Const::logPrint("build index skeleton finished.");

    IO.join();
    // navids have been removed in the process of materializing the 1st layer

    Const::logPrint("Start materialize internal nodes in the 1st layer");
    finished_series = 0;
    finished_percent = 0;
    for (int index = 0; index < Const::vertexNum; ++index) {
        if (nodeIn1stLayer[index].s > Const::th) {
            materializeInterNodeFuzzy(root->ch[index], saxes, nodeIn1stLayer[index].s, *DPNT[index]);

            finished_series += nodeIn1stLayer[index].s;
            double percent = (double) finished_series / (double) internalSize;
            if (percent >= (finished_percent + 1) * 0.1) {
                Const::logPrint(
                    to_string(percent * 100) + "% internal series in the 1st layer has been written to disk.");
                finished_percent = percent * 10;
            }
        }
    }

    Const::logPrint("Fuzzy series number is " + to_string(fuzzy_number));
    Const::logPrint("build index successfully!");

    return root;
}

void Node::growIndexFuzzy(unordered_map<Node *, NODE_RECORDER> &navigating_table, vector<vector<int> > *g) {
    if (size <= Const::th) return;

    int chosen_number = ConversionUtil::findFirstGE(power_2, 1, Const::segmentNum + 1, size / Const::th + 1);
    NODE_RECORDER &parentRecorder = navigating_table[this];

    determineSegments();

    // statistic children info in order to partition
    partUnit nodes[power_2[chosen_number]];
    for (int index = 0; index < power_2[chosen_number]; ++index)
        nodes[index].id = index, nodes[index].s = 0, nodes[index].process_id = -1;
    vector<vector<int> > node_offsets(power_2[chosen_number], vector<int>());
    vector<vector<int> > series_index_list(power_2[chosen_number], vector<int>());
    vector<int> node_actual_size(power_2[chosen_number], 0);

    for (int i = 0; i < parentRecorder.actual_size; ++i) {
        int new_id = ConversionUtil::extendSax(Node::saxes + (long) offsets[i] * Const::segmentNum, bits_cardinality,
                                        chosenSegments);
        nodes[new_id].s++;
        node_offsets[new_id].push_back(offsets[i]);
        series_index_list[new_id].push_back(parentRecorder.series_index_list[i]);
    }

    bool isLeaf[power_2[chosen_number]];
    for (int i = 0; i < power_2[chosen_number]; ++i) {
        isLeaf[i] = false;
        if (nodes[i].s <= Const::th) isLeaf[i] = true;
        node_actual_size[i] = nodes[i].s;
    }

    // divide fuzzy series
    for (int i = parentRecorder.actual_size; i < size; ++i) {
        int new_id = ConversionUtil::extendSax(Node::saxes + (long) offsets[i] * Const::segmentNum, bits_cardinality,
                                        chosenSegments, sax);
        // if the target node needn't have split, then the fuzzy series should not be inserted into it which will result in unnecessary split
        if (isLeaf[new_id] && nodes[new_id].s >= Const::th) continue;
        nodes[new_id].s++;
        node_offsets[new_id].push_back(offsets[i]);
        series_index_list[new_id].push_back(parentRecorder.series_index_list[i]);
    }
    if (layer > 1) vector<int>().swap(offsets);

    int partNum = partition(nodes, chosen_number);

    Node *leaf_children_list[partNum];
    for (int i = 0; i < partNum; ++i) {
        leaf_children_list[i] = new Node(this, i);
        navigating_table.emplace(leaf_children_list[i], NODE_RECORDER());
    }
    ch.resize(power_2[chosen_number]);
    for (int i = 0; i < power_2[chosen_number]; ++i) {
        if (nodes[i].s <= 0) continue;
        else if (nodes[i].s > Const::th || nodes[i].process_id == -1) {
            ch[i] = new Node(this, nodes[i].s, i);
            generateSaxAndCardinality(ch[i], i);
            ch[i]->offsets.resize(nodes[i].s);
            copy(node_offsets[i].begin(), node_offsets[i].end(), ch[i]->offsets.begin());
            navigating_table.emplace(ch[i], NODE_RECORDER(node_actual_size[i], series_index_list[i]));
        } else {
            int _pid = nodes[i].process_id;
            ch[i] = leaf_children_list[_pid];
            leaf_children_list[_pid]->size += nodes[i].s;
            generateSaxAndCardinality4LeafNode(ch[i], i);
            navigating_table[ch[i]].actual_size += node_actual_size[i];
            for (int j = node_actual_size[i]; j < nodes[i].s; ++j)
                navigating_table[ch[i]].series_index_list.push_back(series_index_list[i][j]);
        }
    }


    // note that for internal nodes, table stores all series index list inside while for leaf nodes, only fuzzy series index list are in the tbl
    fuzzy(nodes, node_actual_size, node_offsets, series_index_list, chosen_number, navigating_table);

    vector<vector<int> >().swap(node_offsets);
    vector<vector<int> >().swap(series_index_list);
    vector<int>().swap(node_actual_size);
    navigating_table.erase(this);

    for (auto &child: ch) {
        if (child != nullptr) {
            if (child->size > Const::th)
                child->growIndexFuzzy(navigating_table, g);
            else {
                // this may be executed many times for the same node, but no problem
                if (navigating_table[child].series_index_list.empty())
                    navigating_table.erase(child);
                else {
                    navigating_table[child].actual_size = 0;
                    sort(navigating_table[child].series_index_list.begin(),
                         navigating_table[child].series_index_list.end());
                }
            }
        }
    }
}

void Node::fuzzy(partUnit *part_units, vector<int> &actual_sizes,
                      vector<vector<int> > &node_offsets, vector<vector<int> > &series_index_list,
                      int chosen_num, unordered_map<Node *, NODE_RECORDER> &navigating_tbl) const {
    for (int i = 0; i < power_2[chosen_num]; ++i) {
        fuzzySeriesInPartUnit(part_units, actual_sizes[i], chosen_num, node_offsets[i], series_index_list[i],
                              navigating_tbl, i);
    }
}

struct CAND {
    int id;
    double score;

    CAND(int _i, double _score) {
        id = _i;
        score = _score;
    }

    static bool order(CAND &a, CAND &b) { return a.score < b.score; }
};


void Node::fuzzySeriesInPartUnit(partUnit *part_units, int actual_size, int chosen_num, vector<int> &node_offsets,
                                      vector<int> &series_index_list,
                                      unordered_map<Node *, NODE_RECORDER> &navigating_tbl, int _id) const {
    float *paa;
    double range, lb, ub;
    Node *temporary_node;
    vector<CAND> candidates; // max is chosen num
    int newId, seg, sax_symbol_last_bit, sax_symbol;

    for (int j = 0; j < actual_size; ++j) {
        paa = Node::paas + (long) node_offsets[j] * Const::segmentNum;
        for (int i = 0; i < chosenSegments.size(); ++i) {
            seg = chosenSegments[i];
            newId = _id ^ mask[Const::segmentNum - chosen_num + i];
            sax_symbol_last_bit = (_id & mask[Const::segmentNum - chosen_num + i]) > 0 ? 1 : 0;
            sax_symbol = (sax[seg] << 1) + sax_symbol_last_bit;
            if (part_units[newId].s > 0 && ch[_id] != ch[newId]) {
                temporary_node = ch[newId];
                if (temporary_node->partition_id != -1 && temporary_node->size >= Const::th) continue;

                // find the range of new partition unit in this segment
                // if any bound is infinity, then select the neighbor range as range
                if (sax_symbol == (power_2[bits_cardinality[seg] + 1] - 1))
                    ConversionUtil::getValueRange(sax[seg] << 1, bits_cardinality[seg] + 1, &lb, &ub);
                else if (sax[seg] == 0 && sax_symbol_last_bit == 0)
                    ConversionUtil::getValueRange(1, bits_cardinality[seg] + 1, &lb, &ub);
                else
                    ConversionUtil::getValueRange(sax_symbol, bits_cardinality[seg] + 1, &lb, &ub);
                range = (ub - lb) * Const::fuzzy_f;
                // find whether this series can be added into the candidates list
                // Single precision paa may lead to conflicts with sax, so we need 1e-4 and abs. series are routed by sax to target node.
                if (sax_symbol == (power_2[bits_cardinality[seg] + 1] - 1)) {
                    if (abs(paa[seg] - ub) <= range) candidates.emplace_back(newId, abs(paa[seg] - ub));
                } else if (sax_symbol == 0) {
                    if (abs(lb - paa[seg]) <= range) candidates.emplace_back(newId, abs(lb - paa[seg]));
                } else if (sax_symbol_last_bit == 0) {
                    //                     TODO: check
                    if (abs(ub - paa[seg]) <= range) candidates.emplace_back(newId, abs(ub - paa[seg]));
                } else {
                    if (abs(paa[seg] - lb) <= range) candidates.emplace_back(newId, abs(paa[seg] - lb));
                }
            }
        }
        sort(candidates.begin(), candidates.end(), CAND::order);
        int n = 0;
        // need not update node_offsets and series_index_lists
        for (int index = 0; index < candidates.size() && n < Const::delta; ++index) {
            temporary_node = ch[candidates[index].id];
            if (temporary_node->partition_id != -1) {
                // temp node is a leaf node
                if (temporary_node->size < Const::th) {
                    temporary_node->size++;
                    navigating_tbl[temporary_node].series_index_list.push_back(series_index_list[j]);
                    ++n;
                }
            } else {
                // temporary node is an internal node
                assert(navigating_tbl[temporary_node].actual_size > Const::th);
                temporary_node->size++;
                temporary_node->offsets.push_back(node_offsets[j]);
                navigating_tbl[temporary_node].series_index_list.push_back(series_index_list[j]);
                ++n;
            }
        }
        candidates.clear();
        fuzzy_number += n;
    }
}

void Node::fuzzySeriesInPartUnitInFirstLayer(partUnit *part_units, vector<int> &node_offsets, int _id,
                                                  unordered_map<Node *, NODE_RECORDER> &navigating_tbl,
                                                  vector<vector<double> > &paa_mu_part_units) const {
    float *paa, range;
    Node *temporary_node;
    int newId;
    vector<CAND> candidates;
    for (int j = 0; j < part_units[_id].s; ++j) {
        paa = Node::paas + (long) node_offsets[j] * Const::segmentNum;
        for (int i = 0; i < Const::segmentNum; ++i) {
            newId = _id ^ mask[i];
            if (part_units[newId].s > 0 && ch[_id] != ch[newId]) {
                if (part_units[newId].s < Const::th && ch[newId]->size >= Const::th) continue;
                if (paa[i] > 0) {
                    range = paa_mu_part_units[_id][i] * Const::fuzzy_f_1st;
                    if (paa[i] <= range) candidates.emplace_back(newId, paa[i]);
                } else {
                    range = (-paa_mu_part_units[_id][i]) * Const::fuzzy_f_1st;
                    if (-paa[i] <= range) candidates.emplace_back(newId, -paa[i]);
                }
            }
        }
        sort(candidates.begin(), candidates.end(), CAND::order);
        int n = 0;
        for (int i = 0; i < candidates.size() && n < Const::delta; ++i) {
            temporary_node = ch[candidates[i].id];
            if (temporary_node->partition_id != -1) {
                // temp node is a leaf node
                if (temporary_node->size < Const::th) {
                    temporary_node->size++;
                    navigating_tbl[temporary_node].series_index_list.push_back(node_offsets[j]);
                    ++n;
                }
            } else {
                // temp node is an internal node
                temporary_node->size++;
                temporary_node->offsets.push_back(node_offsets[j]);
                navigating_tbl[temporary_node].series_index_list.push_back(node_offsets[j]);
                ++n;
            }
        }
        candidates.clear();
        fuzzy_number += n;
    } //end loop
}


void Node::fuzzyFirstLayer(partUnit *part_units, const int *nav_ids,
                                unordered_map<Node *, NODE_RECORDER> &navigating_tbl,
                                vector<vector<double> > &paa_mu_part_units) const {
    // data offsets are added in internal node instead of leaf part units
    // node_offsets adds data offsets of leaf part units
    unordered_map<int, vector<int> > node_offsets;
    for (int i = 0; i < Const::vertexNum; ++i) {
        if (part_units[i].s <= Const::th)
            node_offsets[i].reserve(part_units[i].s);
    }
    for (int index = 0; index < this->size; ++index) {
        if (part_units[nav_ids[index]].s <= Const::th)
            node_offsets[nav_ids[index]].push_back(index);
    }
    for (int index = 0; index < Const::vertexNum; ++index) {
        if (ch[index] == nullptr) continue;
        if (this->ch[index]->partition_id == -1)
            fuzzySeriesInPartUnitInFirstLayer(part_units, this->ch[index]->offsets, index, navigating_tbl,
                                              paa_mu_part_units);
        else
            fuzzySeriesInPartUnitInFirstLayer(part_units, node_offsets[index], index, navigating_tbl,
                                              paa_mu_part_units);
    }
    node_offsets.clear();
    for (auto &iter: navigating_tbl) {
        iter.second.actual_size = 0;
        sort(iter.second.series_index_list.begin(), iter.second.series_index_list.end());
    }
}
