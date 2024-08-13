#ifndef IO_H
#define IO_H

// #include "thirdparty/rapidjson/fwd.h"
#include "thirdparty/rapidjson/document.h"
#include "thirdparty/rapidjson/error/en.h"

#include "thirdparty/rapidcsv.hpp"

rapidjson::Document get_json_from_file(std::string file_name);

const char kCsv_sep = ',';
const std::string kCsv_sep_str = std::string(1, kCsv_sep);

std::string get_time_stamp();

class InputParser {
    private:
        std::vector <std::string> tokens;

    public:
        InputParser(int argc, char **argv);
        
        const std::string& getCmdOption(const std::string& option) const;
        bool cmdOptionExists(const std::string& option) const;
};

template <typename T>
void write_to_csv(
    std::string file_name,
    std::vector<std::vector<T>> vov,
    std::vector<std::string> v_names = std::vector<std::string>{}
); 

template <typename T>
void write_to_csv(
    std::string file_name,
    std::vector<T> v,
    std::string v_name = ""
);

template <typename T>
void print_v(std::vector<T> v, std::string v_name = "");

template <typename T>
void print_vov(
    std::vector<std::vector<T>> vov,
    std::vector<std::string> v_names = std::vector<std::string>{}
);

#endif