#ifndef DUMPY_FILEUTIL_H
#define DUMPY_FILEUTIL_H

#include <string>
#include <vector>

class FileUtil {
public:
    static bool checkFileExists(const std::string& name);

    static std::string getFilePath(const std::string& path, const std::string& filename);

    static bool checkDirClean(const std::string& path);

    static long getFileSize(const std::string& fname);

    static void getFiles(const std::string& path, std::vector<std::string>& files);

    static float* readSeries(FILE* f);

    static bool fileRemove(const std::string& fname);

    static void mergeFiles(const std::string sources[], const std::string &dest, int num);
    static void deleteFiles(const std::vector<std::string>& fs);

    static bool createDir(const std::string& path);

    static float* readQueries();
};

#endif //DUMPY_FILEUTIL_H
