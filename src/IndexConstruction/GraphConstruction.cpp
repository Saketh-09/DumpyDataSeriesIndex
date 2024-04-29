

#include <cstdio>
#include <vector>
#include <iostream>
#include <thread>
#include "../../include/DataStructures/GraphConstruction.h"
#include "../../include/Utils/MathUtil.h"
#include "../../include/Utils/FileUtil.h"
int segment_number  = 16, bitsReserve =3;
string graphFileName =  "../RawGraph_" + to_string(segment_number) + "_" + to_string(bitsReserve) + ".bin";

void func(int s, int e, int neighbor_number, int neighborBits, int segment_num, const string & graphFileName){
    int graphSkeleton[neighbor_number];
    FILE *f = fopen((graphFileName + to_string(s)).c_str(), "wb");
    int a = MathUtil::nChooseK(segment_num, 1), b = MathUtil::nChooseK(segment_num, 2), c = MathUtil::nChooseK(segment_num, 3);
    for(int i=s; i< e; ++i)
    {
        if((i - s) % 10000 == 0) cout <<s << ": " << (i - s) <<endl;
        int s = 0;
        if(neighborBits >= 1) {
            MathUtil::get1BitChangedNums(i, segment_num, graphSkeleton, s);
            s += a;
        }
        if(neighborBits >= 2){
            MathUtil::get2BitsChangedNums(i, segment_num, graphSkeleton, s);
            s += b;
        }
        if(neighborBits >= 3){
            MathUtil::get3BitsChangedNums(i, segment_num, graphSkeleton, s);
            s += c;
        }
        if(neighborBits >= 4){
            MathUtil::get4BitsChangedNums(i, segment_num, graphSkeleton, s);
        }
        fwrite(graphSkeleton, sizeof(int), neighbor_number, f);
    }
    fclose(f);

}

void GraphConstruction::buildAndSave2Disk() {

    int arr_length = 1 << segment_number, neighborBits = bitsReserve;
    int neighbor_num = 0;
    for(int index=1;index<=neighborBits;++index)
        neighbor_num += MathUtil::nChooseK(segment_number, index);
    cout << "neighbor number = " << neighbor_num << endl;
    int thread_number = 20, chunk_size = arr_length / thread_number;
    thread threads[thread_number];
    for(int i=0;i<thread_number - 1;++i)
        threads[i] = thread(func, i * chunk_size, (i + 1) * chunk_size, neighbor_num, neighborBits, segment_number, graphFileName);
    threads[thread_number - 1] = thread(func, (thread_number - 1) * chunk_size, arr_length, neighbor_num, neighborBits, segment_number, graphFileName);

    for(int i=0;i<thread_number;++i)
        threads[i].join();

    string sources[thread_number];
    for(int i=0;i<thread_number;++i)
        sources[i] = (graphFileName + to_string(i * chunk_size));
    FileUtil::mergeFiles(sources, graphFileName, thread_number);


}

