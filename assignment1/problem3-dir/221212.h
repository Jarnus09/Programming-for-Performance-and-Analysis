#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
struct Args
{
    std::string input_path;  // R
    int num_producers;       // T
    int min_lines;           // Lmin
    int max_lines;           // Lmax
    int buffer_size;         // M
    std::string output_path; // W
};

bool parseArgs(int argc, char *argv[], Args &args)
{
    if (argc != 7)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <input_file_path> <num_producers> <min_lines_per_thread> "
                  << "<max_lines_per_thread> <buffer_size_in_lines> <output_file_path>\n";
        return false;
    }

    try
    {
        args.input_path = argv[1];
        args.num_producers = std::stoi(argv[2]);
        args.min_lines = std::stoi(argv[3]);
        args.max_lines = std::stoi(argv[4]);
        args.buffer_size = std::stoi(argv[5]);
        args.output_path = argv[6];
    }
    catch (const std::exception e)
    {
        std::cerr << "Erro parsing arguments: " << e.what() << "\n";
        return false;
    }

    if (args.num_producers <= 0)
    {
        std::cerr << "Error: num_producers must be > 0\n";
        return false;
    }
    if (args.max_lines <= 0)
    {
        std::cerr << "Error: max lines must be > 0\n";
        return false;
    }
    if (args.min_lines > args.max_lines)
    {
        std::cerr << "Error: min_lines_per_thread cannot be greater than max_lines_per_thread\n";
        return false;
    }
    if (args.buffer_size <= 0)
    {
        std::cerr << "Error: buffer_size must be > 0\n";
        return false;
    }

    return true;
}

int randomBetween(int a, int b)
{
    srand(time(0));
    return a + rand() % (b - a + 1);
}