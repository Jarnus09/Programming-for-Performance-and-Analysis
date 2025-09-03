#include <iostream>
#include <fstream>
#include <vector>
#include <string>

int main() {
    // File names
    std::vector<std::string> filenames = {
        "file1.txt",
        "file2.txt",
        "file3.txt",
        "file4.txt",
        "file5.txt"
    };

    // Create input.txt and write the file names to it
    std::ofstream inputFile("input.txt");
    if (!inputFile.is_open()) {
        std::cerr << "Error opening input.txt for writing!\n";
        return 1;
    }
    for (const auto& name : filenames) {
        inputFile << name << "\n";
    }
    inputFile.close();

    // Create output file streams
    std::vector<std::ofstream> files;
    for (const auto& name : filenames) {
        files.emplace_back(name);
        if (!files.back().is_open()) {
            std::cerr << "Error opening " << name << " for writing!\n";
            return 1;
        }
    }

    std::string text = "I am god heheh lmao";
    const long long repetitions = 1'000'000;

    // Write the string 1 million times to all files
    for (long long i = 0; i < repetitions; ++i) {
        for (auto& ofs : files) {
            ofs << text << "\n";  // using "\n" instead of std::endl for speed
        }
    }

    // Close all output files
    for (auto& ofs : files) {
        ofs.close();
    }

    std::cout << "All files generated and input.txt written successfully.\n";
    return 0;
}
