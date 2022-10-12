#include "all.h"

int main(){
    std::string matfile_name = "matrix/matrix0.txt";
    std::string outfile_name = "output/result1.txt";
    std::string testfile_name = "test/test0.csv";
    
    dns_eigen('T', 'V', 'y', 'y', 4, 4, matfile_name, outfile_name,testfile_name,10);
}
