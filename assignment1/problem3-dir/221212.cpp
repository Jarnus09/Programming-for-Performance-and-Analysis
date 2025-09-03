#include <iostream>
#include <string>
#include <cstdlib>
#include <thread>
#include <mutex>
#include <queue>
#include <fstream>
#include <atomic>
#include <condition_variable>

#include "221212.h"

std::queue<std::string> buffer_queue;

std::mutex buffer_mutex, file_mutex;
std::condition_variable cv_writer, cv_reader;

std::atomic<bool> done_reading = false;
std::atomic<int> active_producers = 0;

Args args;
const int ANY_PRODUCER_THREAD = -1;

std::atomic<int> writing_thread = ANY_PRODUCER_THREAD;

void produce(std::ifstream &input_stream, int thread_id)
{

    while (true)
    {
        std::queue<std::string> local_buffer;
        {
            std::unique_lock lock(file_mutex);
            active_producers++;

            int number_of_lines_read = 0;
            int lines_to_be_read = randomBetween(args.min_lines, args.max_lines);
            std::string line;

            while (number_of_lines_read != lines_to_be_read && std::getline(input_stream, line))
            {
                number_of_lines_read++;
                local_buffer.push(line);
            }
        }
        {
            std::unique_lock lock(buffer_mutex);

            cv_writer.wait(lock, []
                           { return writing_thread == ANY_PRODUCER_THREAD; });

            writing_thread = thread_id;
            std::string line;

            while (!local_buffer.empty())
            {

                std::string line;
                line = local_buffer.front();
                local_buffer.pop();

                buffer_queue.push(line);

                if (buffer_queue.size() == args.buffer_size)
                {
                    cv_reader.notify_all();
                    cv_writer.wait(lock, [thread_id]
                                   { return writing_thread == thread_id && buffer_queue.empty(); });
                }
            }

            writing_thread = ANY_PRODUCER_THREAD;

            cv_reader.notify_all();

            active_producers--;

            if (done_reading || input_stream.eof())
            {
                done_reading = true;
                cv_reader.notify_all();
                break;
            }
        }
    }

    // Notify all the waiting writers that the entire file has been read
    cv_writer.notify_all();
}

void consume(std::ofstream &output_stream)
{

    while ((!done_reading) || (!buffer_queue.empty() || active_producers != 0))
    {
        
        std::unique_lock lock(buffer_mutex);

        cv_reader.wait(lock, []
                       { return (done_reading) || (!buffer_queue.empty()); });

        if (!buffer_queue.empty())
        {
            std::string line = buffer_queue.front();
            buffer_queue.pop();
            output_stream << line << std::endl;
        }

        if (buffer_queue.empty())
        {
            cv_writer.notify_all();
        }
    }
}

int main(int argc, char *argv[])
{

    if (!parseArgs(argc, argv, args))
    {
        return EXIT_FAILURE;
    }

    std::ifstream input_stream(args.input_path);
    if (!input_stream)
    {
        std::cerr << "Error: Could not open " << args.input_path << " for reading.\n";
        return 1;
    }

    std::ofstream output_stream(args.output_path);
    if (!output_stream)
    {
        std::cerr << "Error: Could not open " << args.output_path << " for writing.\n";
        return 1;
    }

    std::vector<std::thread> producer_threads;
    int num_consumers = std::max(1, args.num_producers / 2);
    std::vector<std::thread> consumer_threads;
    for (int i = 0; i < args.num_producers; i++)
    {
        producer_threads.emplace_back(produce, std::ref(input_stream), i);
    }

    for (int i = 0; i < num_consumers; i++)
    {
        consumer_threads.emplace_back(consume, std::ref(output_stream));
    }

    for (auto &t : producer_threads)
    {
        t.join();
    }
    for (auto &t : consumer_threads)
    {
        t.join();
    }
    return 0;
}