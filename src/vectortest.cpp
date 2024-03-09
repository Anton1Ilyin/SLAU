#include <iostream>
#include<vector>


int main() {
    std::vector<int> vec1={1,2,3,4};
    for(int i=0; i<vec1.size(); i++)
    std::cout<<vec1[i]<<std::endl;
    return 0;
}