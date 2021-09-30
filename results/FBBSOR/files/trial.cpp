#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>

Eigen::VectorXf CSVtoVector(std::string file, int i, int n) {

    double data[i*n];

    Eigen::VectorXf vector = Eigen::VectorXf::Zero(i*n - 1);

    std::ifstream input(file);

    for (int j = 0; j < i*n - 1; j++) {
        input >> data[j];
        vector[j] = data[j];

    }
    
    std::cout << vector << std::endl;

    return vector;
}

Eigen::VectorXf write_w(std::string file, int i, int n) {

    double data[2*(n+1)*(i-1)];

    Eigen::VectorXf vector = Eigen::VectorXf::Zero(2*(n+1)*(i-1));

    std::ifstream input(file);

    for (int j = 0; j < 2*(n+1)*(i-1) - 1; ++j) {
        input >> data[j];
        vector[j] = data[j];

    }
    
    // std::cout << vector << std::endl;

    return vector;
}

// Saving the vector to the file
void WriteFile(Eigen::VectorXf vector, std::string filename) {
    
    std::ofstream myfile (filename.append(".txt"), std::ios_base::app);

    if (myfile.is_open()){

            myfile << vector << "\n";
    }
    myfile.close();
}


int main(){

    int n;
    int i;
    
    std::cout << "Enter the number of time step: " << std::endl;

    std::cin >> n;

    std::cout << "Enter the number of space grids: " << std::endl;

    std::cin >> i;

    Eigen::VectorXf Y = CSVtoVector("Y.txt", i-1, n+1);

    std::cout << Y;
    WriteFile(Y, "new_y");

    return 0;
}